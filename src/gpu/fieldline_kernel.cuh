// Host-callable declarations for the CUDA field-line tracing kernels.
// Include this from structure.cxx / laminar_mpi.cxx under #ifdef USE_GPU.
// Compiled by nvcc as part of fieldline_kernel.cu.

#ifndef FIELDLINE_KERNEL_CUH
#define FIELDLINE_KERNEL_CUH

#include <gpu_fields.hxx>

// ---------------------------------------------------------------------------
// Laminar tracing — one GPU thread per (R,Z) point.
//
// init[N]:     initial conditions (R,Z,phi) for every point
// results[N]:  output scalars (ntor, Lc, psimin, psiav) per point
// grid:        pre-sampled B-field and psi on a 3-D (R,phi,Z) grid
// params:      integration parameters (dpinit, itt, MapDirection, …)
//
// Returns 0 on success, -1 on CUDA error.
// ---------------------------------------------------------------------------
int gpu_trace_laminar(const FieldlineInit* init,
                      LaminarResult*       results,
                      int                  N,
                      const FieldGrid3D*   grid,
                      const GPUTraceParams& params);

// ---------------------------------------------------------------------------
// Structure tracing — one GPU thread per field line, per direction.
//
// init[N]:            initial conditions
// steps_bwd[N]:       number of valid steps in backward direction per field line
// steps_fwd[N]:       same for forward direction
// results_bwd[N*max]: StructureStep array, layout [iLine * max_steps + iStep]
// results_fwd[N*max]: same for forward direction
// max_steps:          caller-allocated capacity per field line per direction
// grid / params:      same as laminar
//
// Returns 0 on success, -1 on CUDA error.
// ---------------------------------------------------------------------------
int gpu_trace_structure(const FieldlineInit* init,
                        int                  N,
                        int                  max_steps,
                        StructureStep*       results_bwd,
                        StructureStep*       results_fwd,
                        int*                 steps_bwd,
                        int*                 steps_fwd,
                        const FieldGrid3D*   grid,
                        const GPUTraceParams& params);

// ---------------------------------------------------------------------------
// Number of CUDA-capable devices visible to the runtime, or -1 if the CUDA
// runtime reports an error (e.g. no driver / no GPU). Used to fail clearly
// when -g is requested on a host without a usable GPU.
// ---------------------------------------------------------------------------
int gpu_device_count(void);

#endif // FIELDLINE_KERNEL_CUH
