// CPU-side helpers that build a FieldGrid3D ready to be handed to the GPU kernels.
// Compiled by g++, not nvcc — no CUDA headers required here.

#ifndef BFIELD_SAMPLER_HXX
#define BFIELD_SAMPLER_HXX

#include <string>
#include <gpu_fields.hxx>
#include <efit_class.hxx>
#include <io_class.hxx>

// ---------------------------------------------------------------------------
// Sample getBfield() on a uniform (R, phi, Z) grid using OpenMP.
//
// Uses whatever equilibrium backend is already loaded (EFIT, M3DC1, GPEC,
// VMEC, …) via the same getBfield() function the CPU path calls.  All
// backends therefore work transparently.
//
// If Nphi == 1 the single phi plane is PAR.phistart; use this for axisymmetric
// equilibria (EFIT-only, response_field < 0 without 3D perturbations).
// For 3D cases pass Nphi ≥ 16 (64 is a practical default).
//
// Returns a heap-allocated FieldGrid3D.  Caller must call free_field_grid().
// ---------------------------------------------------------------------------
FieldGrid3D* sample_bfield(EFIT& EQD, IO& PAR,
                            int NR,   double Rmin,   double Rmax,
                            int Nphi, double phimin, double phimax,
                            int NZ,   double Zmin,   double Zmax);

// ---------------------------------------------------------------------------
// Fill grid->wall_R / wall_Z / Nwall from EQD.wall (2-D wall polygon).
// grid must already be allocated; its wall pointers are replaced.
// ---------------------------------------------------------------------------
void extract_wall(EFIT& EQD, FieldGrid3D* grid);

// ---------------------------------------------------------------------------
// Fill grid->psi_norm (NR × NZ) by calling EQD.get_psi() on the same
// (R, Z) grid used for the B-field.  Call after sample_bfield().
// ---------------------------------------------------------------------------
void fill_psi_grid(EFIT& EQD, FieldGrid3D* grid);

// ---------------------------------------------------------------------------
// Free all heap memory owned by grid, then delete grid itself.
// ---------------------------------------------------------------------------
void free_field_grid(FieldGrid3D* grid);

#endif // BFIELD_SAMPLER_HXX
