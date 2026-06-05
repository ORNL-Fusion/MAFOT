# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

MAFOT (MAnifold & FOotprint Tracer / MAgnetic Fieldlines Of Tokamaks) is a C++/Fortran MPI+OpenMP code for tracing magnetic field lines and computing connection lengths, Poincaré maps, and divertor heat footprints in tokamak plasmas. It supports multiple machines (DIII-D, ITER, NSTX, MAST, C-Mod, TCABR, and a generic "any" machine), with a Python/Tkinter GUI frontend. The current branch (`HEATv4.1`) adds a HEAT target for use as a field-line tracing backend for the HEAT heat-flux analysis code.

## Build system

All builds require a `make.inc` in the repo root. Copy and adapt from `install/make.inc.*` for your platform; the HEAT/Docker variant is `install/make.inc.HEAT` (also the current root `make.inc`).

```bash
# Build all tools for all machines + GUI + xpand
make all

# Build only HEAT targets (heatstructure + heatlaminar_mpi)
make heat

# Build a specific machine (e.g. DIII-D)
make d3d

# Build a single tool type across all machines
make laminar   # all *laminar_mpi binaries
make structure # all *structure binaries
make foot      # all *foot_mpi binaries

# Build the optional GUI symlink
make gui

# Clean compiled objects
make clean
```

Compiled binaries go to `$(BIN_DIR)` and libraries to `$(LIB_DIR)` as defined in `make.inc`. The HEAT make.inc uses `/root/source/MAFOT/build/{bin,lib}`.

Machine selection is controlled by preprocessor defines (`-DD3D`, `-DITER`, `-DNSTX`, `-DMAST`, `-DCMOD`, `-DTCABR`, `-DHEAT`, `-DANYM`). Optional 3rd-party backends are controlled by `make.inc` variables:
- `M3DC1 = True` — enables M3D-C1 resistive MHD perturbations (`-Dm3dc1`)
- `VMEC = True` — enables VMEC/XFIELD flux-surface coordinates (`-DUSE_XFIELD`)
- `SIESTA = True` — enables SIESTA (implies VMEC)

## Dependencies

- **blitz++** — array library, bundled in `blitz++/` (must be built first)
- **openmpi** — MPI parallelism
- **HDF5** — required for M3DC1 (openmpi-flavored HDF5)
- **netcdf** — required for VMEC
- **gfortran** — Fortran 77 coil-geometry subroutines in `src/libtrip3d/`
- **Python 3 + Tkinter/ttk** — optional GUI

## Running a case

Each tool takes a *control/parameter file* (filename conventionally starts with `_`) that specifies grid, field options, shot/time, and coil settings. Example: `tests/_plot.dat`, `examples/testcase/_plot.dat`.

```bash
# Serial tools (plot, fix, man, structure, lcfs)
dtstructure _lam.dat [tag]

# Parallel tools (laminar, foot, plot_mpi, trace)
mpirun -np <N> -host <hostname> dtlaminar_mpi _lam.dat [tag]
mpirun -np <N> -host <hostname> heatlaminar_mpi _lam.dat [tag]
```

The optional trailing `tag` string is appended to output filenames.

Common command-line flags (supported by `laminar_mpi`, `structure`, and most tools):
- `-h` — show help
- `-b` — use a simple bounding box instead of the g-file wall outline
- `-i dpinit` — set integration step size in degrees (default 1.0)
- `-B Rmin,Rmax,Zmin,Zmax` — override boundary box
- `-W wall_file` — use custom wall geometry file
- `-I island` — specify fake-island perturbation file
- `-P points` — specify initial-points file (structure tool)

The g-file (EFIT equilibrium) path and shot/time are set inside the control file. The supplemental coil-current file (e.g. `diiidsup.in`) must be in the working directory for machines that use it.

## Architecture

### Header-only class hierarchy (`include/`)

All core logic lives in `.hxx` headers; there are no separate `.cpp` implementation files for the classes.

| Header | Role |
|---|---|
| `mafot.hxx` | Top-level include that pulls in all classes and machine header; defines global `getBfield_general()` and boundary helpers |
| `efit_class.hxx` | Reads g-files; bicubic interpolation of ψ, B-field, flux surfaces |
| `io_class.hxx` | Reads the control parameter file; holds all run parameters |
| `particle_class.hxx` | Field-line / particle state; Runge-Kutta integrator (`dpinit` step in φ) |
| `mafot.hxx` (machine selection) | Conditionally includes one of `d3d.hxx`, `iter.hxx`, `nstx.hxx`, `mast.hxx`, `anymachine.hxx`, `tcabr.hxx`, `heat.hxx` |
| `heat.hxx` | Machine-agnostic HEAT backend; sets a wide `bndy[]`, stubs `point_along_target()`, omits machine-specific coil calls |
| `m3dc1_class.hxx` | M3D-C1 3D-field interface (when `M3DC1=True`) |
| `gpec_class.hxx` | GPEC 3D-field interface |
| `xfield_class.hxx`, `vmec_class.hxx` | VMEC/XFIELD flux-coordinate interface (when `VMEC=True`) |
| `splines.hxx` | Bicubic spline routines used by `efit_class` |
| `andi.hxx` | Utility macros, constants, string helpers |

### Executable sources (`src/`)

Each `.cxx` file compiles into multiple binaries — one per machine — via the `$(DEFINES) -D<MACHINE>` flag:

| Source | Parallelism | Output |
|---|---|---|
| `laminar_mpi.cxx` | MPI + OpenMP | Connection-length Lc and ψ maps on an R-Z grid |
| `foot_mpi.cxx` | MPI + OpenMP | Heat/particle footprints on target plates |
| `plot_mpi.cxx` | MPI + OpenMP | Poincaré / field-line plots |
| `structure.cxx` | Serial | Full 3D field-line paths (for visualization) |
| `trace.cxx` | MPI + OpenMP | General-purpose field-line tracing |
| `lcfs.cxx` | Serial | Last closed flux surface finder |
| `fix.cxx` | Serial | X-point / O-point finder |
| `man.cxx` | Serial | Manifold tracing |
| `xpand_mpi.cxx` | MPI + OpenMP | VMEC flux-surface expansion (requires `VMEC=True`) |

### Fortran coil libraries (`src/libtrip3d/`)

Biot-Savart and coil-geometry routines in Fortran 77, compiled into `libtrip3d.a`. Machine-specific geometry files: `d3dicoils.f`, `itericoilsgeom.f`, `nstxecgeom.f`, `masteccoilsgeom.f`, `tcabr*.f`, etc. The HEAT machine uses `nstxecgeom.f` stubs (geometry is effectively disabled — only EFIT 2D equilibrium is used).

### MPI+OpenMP parallelism pattern

The `_mpi` tools use a master/slave pattern: MPI rank 0 runs two OpenMP threads — a master thread dispatching work rows (in Z or t) to MPI slaves, and a slave thread doing local computation. All other MPI ranks compute and return results. This pattern is replicated identically across all `_mpi` tools.

### HEAT integration notes

The `-DHEAT` define routes through `heat.hxx` instead of any real machine header. Key differences from standard machines:
- `bndy[]` covers all tokamaks up to DEMO size (Rmax=13 m, |Z|≤10 m)
- No machine-specific coil perturbations; only EFIT 2D + optional M3DC1 3D fields
- `point_along_target()` is a no-op stub (HEAT handles geometry itself)
- Use `-b` flag on the command line if you want simple bounding-box limits instead of g-file wall

## GPU support (branch `mafot_gpu`)

The GPU path pre-samples whichever equilibrium backend is loaded (EFIT, M3DC1, VMEC, GPEC, SIESTA) onto a flat 3-D `(R, phi, Z)` host array, uploads it to the GPU, and runs a massively parallel RK4 integrator — one CUDA thread per field line. All existing CPU code paths are unchanged; the GPU is opt-in at both compile time and runtime.

### Build

In `make.inc` (or `install/make.inc.HEAT`):

```makefile
GPU = True
CUDA_PATH = /usr/local/cuda
NVCCFLAGS = -O3 -std=c++14 -arch=sm_89   # sm_89=Ada (RTX 2000/40xx); sm_80=A100; sm_70=V100; sm_86=RTX30xx; sm_90=H100
```

```bash
make heat    # GPU-enabled heatstructure + heatlaminar_mpi (one binary per tool)
```

There are **no separate `_gpu` binaries**: with `GPU=True`, GPU support is compiled
directly into `heatstructure` and `heatlaminar_mpi`. The same binary runs on the GPU
when given `-g` and on the CPU when `-g` is omitted. With `GPU=False`, the identical
target names build CPU-only binaries that need no CUDA toolchain (and print an error
if given `-g`).

GPU builds add `-DUSE_GPU` and link `libcudart` **statically** (`-lcudart_static`), so
the runtime/deployment image needs no CUDA libraries — only the driver (`libcuda.so`,
injected by `--gpus`) when `-g` is actually used. Source files involved:
- `include/gpu_fields.hxx` — POD structs shared between host and device
- `src/gpu/bfield_sampler.hxx/.cxx` — CPU-side grid sampling (samples `getBfield()` onto a 3-D grid)
- `src/gpu/fieldline_kernel.cuh/.cu` — CUDA kernels and host wrappers

### Runtime

Add `-g` to use the GPU; omit it to use the normal CPU/MPI path (same binary):

```bash
heatstructure -g _str.dat [tag]                    # GPU; samples the active field source
heatstructure _str.dat [tag]                       # CPU (no -g)
mpirun -n 1 heatlaminar_mpi -g _lam.dat [tag]      # GPU laminar (single rank)
mpirun -n N heatlaminar_mpi _lam.dat [tag]         # CPU laminar (no -g)
```

If `-g` is given but no CUDA-capable GPU is present, the tool prints
`ERROR: -g requires a CUDA-capable GPU, but none was found.` and exits.

The GPU samples whatever field source `response_field` selects in the control file (EFIT, M3D-C1, XFIELD/XPAND, SIESTA, GPEC) onto a uniform 3-D `(R, phi, Z)` grid, then runs the RK4 kernel on that grid. To trace a **user-supplied 3-D field**, set `response_field = -3` and provide an `xpand.dat` ASCII grid file — the existing `XFIELD` reader loads it and the GPU samples it via `getBfield()` like any other source (this needs `VMEC = True` in `make.inc`, which compiles in the XFIELD/XPAND reader). `response_field = -1` (axisymmetric EFIT) uses a single phi plane; all other sources use a full 3-D grid.

Key limitations:
- Structure tool (`heatstructure -g`) bypasses the GPU path when `sigma != 0` (heat-flux weighting), falling back to the CPU automatically.
- Laminar GPU mode (`heatlaminar_mpi -g`) runs entirely on MPI rank 0; passing `-n 1` is sufficient (other ranks exit immediately).
- GPU + M3DC1 works only on rank 0; if `prepare_common_perturbations` requires MPI coordination in your M3DC1 configuration, use CPU mode.

### Python GUI (`python/`)

`mafot_gui.py` is a Tkinter frontend that constructs control files and launches MPI jobs. `d3dplot.py` handles post-processing visualization. `use_xpand.py` is a helper for VMEC expansion. These are symlinked into `$(BIN_DIR)` by `make gui` / `make xpand` / `make d3dplot`.
