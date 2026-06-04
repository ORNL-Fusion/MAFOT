// CPU-side B-field pre-sampling for GPU field-line tracing.
// Compiled by g++.  No CUDA headers required.

#include <mafot.hxx>      // getBfield(), bndy[]  (must precede bfield_sampler.hxx:
                          //   sets up <vector>/std + la_string before efit_class.hxx)
#include <bfield_sampler.hxx>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

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
