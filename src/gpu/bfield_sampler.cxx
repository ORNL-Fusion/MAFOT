// CPU-side B-field pre-sampling for GPU field-line tracing.
// Compiled by g++.  No CUDA headers required.

#include <bfield_sampler.hxx>
#include <mafot.hxx>      // getBfield(), bndy[]
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#ifdef USE_NETCDF
#include <netcdf.h>
#endif

// ---------------------------------------------------------------------------
// Internal helper: allocate a zero-initialised FieldGrid3D
// ---------------------------------------------------------------------------
static FieldGrid3D* alloc_grid(int NR, int Nphi, int NZ)
{
    FieldGrid3D* g = new FieldGrid3D();
    memset(g, 0, sizeof(FieldGrid3D));
    long Ntotal = (long)NR * Nphi * NZ;
    g->BR    = new double[Ntotal]();
    g->Bphi  = new double[Ntotal]();
    g->BZ    = new double[Ntotal]();
    g->psi_norm = new double[(long)NR * NZ]();
    g->NR    = NR;
    g->Nphi  = Nphi;
    g->NZ    = NZ;
    // wall filled separately
    g->wall_R = nullptr;
    g->wall_Z = nullptr;
    g->Nwall  = 0;
    return g;
}

// ---------------------------------------------------------------------------
// sample_bfield
// ---------------------------------------------------------------------------
FieldGrid3D* sample_bfield(EFIT& EQD, IO& PAR,
                            int NR,   double Rmin,   double Rmax,
                            int Nphi, double phimin, double phimax,
                            int NZ,   double Zmin,   double Zmax)
{
    FieldGrid3D* g = alloc_grid(NR, Nphi, NZ);

    g->Rmin   = Rmin;   g->Rmax   = Rmax;
    g->phimin = phimin; g->phimax = phimax;
    g->Zmin   = Zmin;   g->Zmax   = Zmax;

    g->dR   = (NR   > 1) ? (Rmax   - Rmin)   / (NR   - 1) : 0.0;
    g->dphi = (Nphi > 1) ? (phimax - phimin)  / (Nphi - 1) : 0.0;
    g->dZ   = (NZ   > 1) ? (Zmax   - Zmin)    / (NZ   - 1) : 0.0;

    g->RmAxis = EQD.RmAxis;
    g->ZmAxis = EQD.ZmAxis;

    for(int i = 0; i < 4; i++) g->bndy[i] = bndy[i];
    g->simpleBndy = simpleBndy;

    printf("GPU: sampling B-field on %d x %d x %d (R x phi x Z) grid ...\n", NR, Nphi, NZ);
    fflush(stdout);

    // OpenMP parallel loop over grid points
    #pragma omp parallel for schedule(dynamic,4) collapse(2)
    for(int iR = 0; iR < NR; iR++)
    {
        for(int iphi = 0; iphi < Nphi; iphi++)
        {
            double R       = Rmin + iR   * g->dR;
            double phi_deg = (Nphi > 1) ? (phimin + iphi * g->dphi) : PAR.phistart;
            double phi_rad = phi_deg * 0.017453292519943;  // getBfield expects radians

            for(int iZ = 0; iZ < NZ; iZ++)
            {
                double Z = Zmin + iZ * g->dZ;
                double BR_val, BZ_val, Bphi_val;

                int chk = getBfield(R, Z, phi_rad, BR_val, BZ_val, Bphi_val, EQD, PAR);
                if(chk != 0) { BR_val = 0.0; BZ_val = 0.0; Bphi_val = 1.0; }  // fallback avoids /0 in kernel

                long idx = (long)iR * (Nphi * NZ) + (long)iphi * NZ + iZ;
                g->BR  [idx] = BR_val;
                g->Bphi[idx] = Bphi_val;
                g->BZ  [idx] = BZ_val;
            }
        }
    }

    printf("GPU: B-field sampling complete.\n");
    fflush(stdout);
    return g;
}

// ---------------------------------------------------------------------------
// read_bfield_netcdf
// ---------------------------------------------------------------------------
FieldGrid3D* read_bfield_netcdf(const std::string& filename)
{
#ifndef USE_NETCDF
    fprintf(stderr, "ERROR: MAFOT was not compiled with netCDF support (USE_NETCDF).\n");
    fprintf(stderr, "       Rebuild with VMEC=True or add -DUSE_NETCDF to CFLAGS.\n");
    return nullptr;
#else
    int ncid, retval;

    if((retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid)) != NC_NOERR) {
        fprintf(stderr, "ERROR: cannot open netCDF file '%s': %s\n",
                filename.c_str(), nc_strerror(retval));
        return nullptr;
    }

    // --- Read dimension sizes ---
    int nr_dimid, nphi_dimid, nz_dimid;
    size_t NR, Nphi, NZ;

    #define NC_CHK(call) do { if((retval=(call))!=NC_NOERR) { \
        fprintf(stderr,"netCDF error: %s\n",nc_strerror(retval)); nc_close(ncid); return nullptr; } } while(0)

    NC_CHK(nc_inq_dimid(ncid, "nr",   &nr_dimid));
    NC_CHK(nc_inq_dimid(ncid, "nphi", &nphi_dimid));
    NC_CHK(nc_inq_dimid(ncid, "nz",   &nz_dimid));
    NC_CHK(nc_inq_dimlen(ncid, nr_dimid,   &NR));
    NC_CHK(nc_inq_dimlen(ncid, nphi_dimid, &Nphi));
    NC_CHK(nc_inq_dimlen(ncid, nz_dimid,   &NZ));

    FieldGrid3D* g = alloc_grid((int)NR, (int)Nphi, (int)NZ);

    // --- Read coordinate arrays ---
    double *R_arr   = new double[NR];
    double *phi_arr = new double[Nphi];
    double *Z_arr   = new double[NZ];

    int varid;
    NC_CHK(nc_inq_varid(ncid, "R",   &varid)); NC_CHK(nc_get_var_double(ncid, varid, R_arr));
    NC_CHK(nc_inq_varid(ncid, "phi", &varid)); NC_CHK(nc_get_var_double(ncid, varid, phi_arr));
    NC_CHK(nc_inq_varid(ncid, "Z",   &varid)); NC_CHK(nc_get_var_double(ncid, varid, Z_arr));

    g->Rmin   = R_arr[0];     g->Rmax   = R_arr[NR-1];
    g->phimin = phi_arr[0];   g->phimax = phi_arr[Nphi-1];
    g->Zmin   = Z_arr[0];     g->Zmax   = Z_arr[NZ-1];
    g->dR   = (NR   > 1) ? (g->Rmax   - g->Rmin)   / (NR   - 1) : 0.0;
    g->dphi = (Nphi > 1) ? (g->phimax - g->phimin)  / (Nphi - 1) : 0.0;
    g->dZ   = (NZ   > 1) ? (g->Zmax   - g->Zmin)    / (NZ   - 1) : 0.0;

    delete[] R_arr;
    delete[] phi_arr;
    delete[] Z_arr;

    // --- Read B-field arrays ---
    NC_CHK(nc_inq_varid(ncid, "BR",   &varid)); NC_CHK(nc_get_var_double(ncid, varid, g->BR));
    NC_CHK(nc_inq_varid(ncid, "Bphi", &varid)); NC_CHK(nc_get_var_double(ncid, varid, g->Bphi));
    NC_CHK(nc_inq_varid(ncid, "BZ",   &varid)); NC_CHK(nc_get_var_double(ncid, varid, g->BZ));

    // --- Optional psi array ---
    if(nc_inq_varid(ncid, "psi", &varid) == NC_NOERR) {
        NC_CHK(nc_get_var_double(ncid, varid, g->psi_norm));
    } else {
        // Mark unavailable — laminar will skip psimin/psiav columns
        long Npsi = (long)NR * NZ;
        for(long k = 0; k < Npsi; k++) g->psi_norm[k] = -1.0;
    }

    nc_close(ncid);
    #undef NC_CHK

    printf("GPU: read B-field from '%s' (%zu x %zu x %zu grid).\n",
           filename.c_str(), NR, Nphi, NZ);
    return g;
#endif // USE_NETCDF
}

// ---------------------------------------------------------------------------
// extract_wall
// ---------------------------------------------------------------------------
void extract_wall(EFIT& EQD, FieldGrid3D* grid)
{
    delete[] grid->wall_R;
    delete[] grid->wall_Z;

    int N = EQD.Nwall;
    grid->Nwall  = N;
    grid->wall_R = new double[N];
    grid->wall_Z = new double[N];

    // EQD.wall is 1-indexed blitz, alternating R,Z: wall(1)=R1, wall(2)=Z1, ...
    for(int k = 0; k < N; k++) {
        grid->wall_R[k] = EQD.wall(2*k + 1);
        grid->wall_Z[k] = EQD.wall(2*k + 2);
    }
}

// ---------------------------------------------------------------------------
// fill_psi_grid
// ---------------------------------------------------------------------------
void fill_psi_grid(EFIT& EQD, FieldGrid3D* grid)
{
    int NR = grid->NR;
    int NZ = grid->NZ;

    #pragma omp parallel for collapse(2) schedule(static)
    for(int iR = 0; iR < NR; iR++) {
        for(int iZ = 0; iZ < NZ; iZ++) {
            double R = grid->Rmin + iR * grid->dR;
            double Z = grid->Zmin + iZ * grid->dZ;
            double psi_val, dummy1, dummy2;
            int chk = EQD.get_psi(R, Z, psi_val, dummy1, dummy2);
            grid->psi_norm[(long)iR * NZ + iZ] = (chk == 0) ? psi_val : 1.5;
        }
    }
}

// ---------------------------------------------------------------------------
// free_field_grid
// ---------------------------------------------------------------------------
void free_field_grid(FieldGrid3D* grid)
{
    if(!grid) return;
    delete[] grid->BR;
    delete[] grid->Bphi;
    delete[] grid->BZ;
    delete[] grid->psi_norm;
    delete[] grid->wall_R;
    delete[] grid->wall_Z;
    delete grid;
}
