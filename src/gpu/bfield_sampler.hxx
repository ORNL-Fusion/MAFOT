// CPU-side helpers that build a FieldGrid3D ready to be handed to the GPU kernels.
//
// Header-only (like the rest of MAFOT): these functions are `inline` and compile
// directly into the including translation unit. There is intentionally NO
// bfield_sampler.cxx / bfield_sampler.o -- a separate object that #included
// mafot.hxx would emit its own copy of every (non-inline) MAFOT symbol and
// collide with structure.o / laminar_mpi.o at link time ("multiple definition").
//
// MUST be included AFTER <mafot.hxx> (it uses getBfield(), bndy[], simpleBndy).
// Compiled by g++, not nvcc -- no CUDA headers required here.

#ifndef BFIELD_SAMPLER_HXX
#define BFIELD_SAMPLER_HXX

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <gpu_fields.hxx>
#include <efit_class.hxx>
#include <io_class.hxx>

// ---------------------------------------------------------------------------
// Internal helper: allocate a zero-initialised FieldGrid3D
// ---------------------------------------------------------------------------
inline FieldGrid3D* alloc_grid(int NR, int Nphi, int NZ)
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
// Sample getBfield() on a uniform (R, phi, Z) grid using OpenMP.
//
// Uses whatever equilibrium backend is already loaded (EFIT, M3DC1, GPEC,
// VMEC/XFIELD, …) via the same getBfield() function the CPU path calls.  All
// backends therefore work transparently.
//
// If Nphi == 1 the single phi plane is PAR.phistart; use this for axisymmetric
// equilibria (EFIT-only, response_field < 0 without 3D perturbations).
// For 3D cases pass Nphi >= 16 (64 is a practical default).
//
// Returns a heap-allocated FieldGrid3D.  Caller must call free_field_grid().
// ---------------------------------------------------------------------------
inline FieldGrid3D* sample_bfield(EFIT& EQD, IO& PAR,
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
// Fill grid->wall_R / wall_Z / Nwall from EQD.wall (2-D wall polygon).
// grid must already be allocated; its wall pointers are replaced.
// ---------------------------------------------------------------------------
inline void extract_wall(EFIT& EQD, FieldGrid3D* grid)
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
// Fill grid->psi_norm (NR × NZ) by calling EQD.get_psi() on the same
// (R, Z) grid used for the B-field.  Call after sample_bfield().
// ---------------------------------------------------------------------------
inline void fill_psi_grid(EFIT& EQD, FieldGrid3D* grid)
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
// Free all heap memory owned by grid, then delete grid itself.
// ---------------------------------------------------------------------------
inline void free_field_grid(FieldGrid3D* grid)
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

#endif // BFIELD_SAMPLER_HXX
