# README.md

This file gives a brief overview of the MAFOT code and its history. It also provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

MAFOT (MAnifold & FOotprint Tracer / MAgnetic Fieldlines Of Tokamaks) is a C++/Fortran MPI+OpenMP code for tracing magnetic field lines and computing connection lengths, Poincaré maps, and divertor heat footprints in tokamak plasmas. It supports multiple machines (DIII-D, ITER, NSTX, MAST, C-Mod, TCABR, and a generic "any" machine), with a Python/Tkinter GUI frontend. The current branch (`HEATv4.1`) adds a HEAT target for use as a field-line tracing backend for the HEAT heat-flux analysis code.

## Developer and References

Code Author and primary developer: Andreas Wingen   
Email: wingen@fusion.gat.com

To cite MAFOT please use:   
A. Wingen, T. E. Evans, and K. H. Spatschek, “High resolution numerical studies of separatrix splitting due to non-axisymmetric perturbation in DIII-D,” Nuclear Fusion 49, 55027, (2009), doi: 10.1088/0029-5515/49/5/055027.

Some other references that use the MAFOT code: 
- A. Wingen et al., Physics of Plasmas 16, 42504, (2009), doi: 10.1063/1.3099053.
- A. Wingen et al., Physical Review Letters 104, 175001, (2010), doi: 10.1103/physrevlett.104.175001.
- A. Wingen et al., Physics of Plasmas 21, 12509, (2014), doi: 10.1063/1.4862034.
- A. Wingen et al., Nuclear Fusion 54, 64007, (2014), doi: 10.1088/0029-5515/54/6/064007.
- A. Wingen et al., Nuclear Fusion 61, 16018, (2021), doi: 10.1088/1741-4326/abbfe9.
- T. Looby et al., Fusion Science and Technology 78, 10–27, (2022), doi: 10.1080/15361055.2021.1951532.

## Current Version

- MAFOT 6.1
- GUI 3.21

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
- `src/gpu/bfield_sampler.hxx` — CPU-side grid sampling (samples `getBfield()` onto a 3-D grid); header-only (inline), compiled into structure.o / laminar_mpi.o
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

## Python GUI (`python/`)

`mafot_gui.py` is a Tkinter frontend that constructs control files and launches MPI jobs. `d3dplot.py` handles post-processing visualization. `use_xpand.py` is a helper for VMEC expansion. These are symlinked into `$(BIN_DIR)` by `make gui` / `make xpand` / `make d3dplot`.

## Version Notes

MAFOT 6.1 -- July 2026
- GPU support added for HEAT specific binaries only for now.
- matplotlib get_cmap deprecation fixed in d3dplot.py

MAFOT 6.0 -- July 2026
- OpenMPI 5.0 compatibility patch: All C++ bindings (deprecated in 5.0) removed and replaced with C bindings.
- libmpi_cxx removed from all make files, as no longer needed.
- Backwards compatible with OpenMPI 4.x

GUI 3.21 -- Aug 2025
- Added Omega cluster Slurm Queue batch file support

MAFOT 5.7 -- May 2025
- Removed the simpleBndy exception for HEAT. Use the -b command line option instead.

MAFOT 5.6 -- Nov 2024
- Added GPEC support
- Added higher order derivatives to bicubic interpolation. Note that interpolation itself is still cubic, so accuracy is only slightly improved. 
- Enable variable integrator step size via command line flag in plot, laminar, foot and structure.

GUI 3.2 -- Nov 2024
- Added GPEC support

GUI 3.1.1 -- Oct 2024
- Bug fix: GUI would not start if gPath in control file does not exist. This is now fixed.

MAFOT 5.5.2 -- Feb 2024
- for compile setups without 3rd party code supports, like M3D-C1, MAFOT now defaults back to using the gFile mode even if other (unavailable) Field modes are selected in the control file. Before, this was not captured and a run would fail doing nothing. 

MAFOT 5.5.1 -- Dec 2023
- laminar now returns B-field components instead of previously obsolete variables (pitch, yaw amd psimax)

MAFOT 5.5 -- May 2023
- Code can now read arbitrarily named gFiles
- bug fix in structure, where toroidal angle was wrong for dpinit != 1

GUI 3.1 -- Jan 2023
- gFile is now a full pathname, selectable from drop down menu, and sets shot & time, if applicable

MAFOT 5.4.1 -- Feb 2023
- Enable any machine without specific 3D coils

GUI 3.0 -- Jan 2023
- converted to Python3
- PPPL cluster support added

MAFOT 5.4 -- Oct 2021
- Collision class (SULI project) merged into main

MAFOT 5.3 -- Oct 2021
- HEAT inluded in main file structure

MAFOT 5.2 -- Aug 2021
- TCABR tokamak support added, also to d3dplot.py

MAFOT 5.11 -- Nov 2019
- sheath model added

MAFOT 5.1 -- July 2019
- new tools trace and lcfs added
- bug fixes

MAFOT 5.0 -- Jan 2019
- Zeff added to all input files. This required some minor restructuring, which makes old ITER and NSTX input files incompatible.
  For D3D, only the position of useFilaments shifted by one. MAST and C-Mod are unchanged
- electric field for ExB particle drifts added

GUI 2.1 -- Jan 2019
- Zeff added
- electric field input added

MAFOT 4.2 -- Sep 2018
- restart option added to xpand, foot and laminar
- VMEC field line tracing now posible inside s = 1 without xpand
- d3dplot uses now pcolormesh instead of imshow

GUI 2.1 -- Jul 2018
- support for full wall outline in foot, using Swall
- Mapdirection now selectable in foot, because of reversed Bt cases

MAFOT 4.11 -- Jul 2018
- special 2D bisection method added to support s,u calculation from R,Z in xpand

MAFOT 4.1 -- May 2018
- laminar now set up to work with Sheft toolkit on NSTX. Several fixes for NSTX
- 3D wall with any number of toroidal slices
- check for wall crossing within first integration step

MAFOT 4.03 -- Mar 2018
- full wall now accessible for foot_mpi using the swall coordinate and target = 0

GUI 2.0 -- Feb 2018
- complete restructuring of ther code to unify common code segments
- several new features, like use of shell Flags for VMEC/SIESTA, or dtlaminar
- new autocomplete Entry widget for file names, that searches through the directory structure

MAFOT 4.02 -- Feb 2018
- B-coil shift&tilt error added for DIII-D

MAFOT 4.01 -- Jan 2018
- Buswork error field added for DIII-D
- For backward compatibility: c++11 standard no longer required, but recommended
- structure can now use the 3D wall
- compile targets for "each tool in all machines" added to makefile

MAFOT 4.0 -- Dec 2017
- unification of M3DC1 and regular machine file -> *_m3dc1.hxx machine files no longer supported
- simplification of machine specific files by moving larg pieces of common code to mafot.hxx
- bug fix for xpand
- CMOD now fully available (except target plate specification for footprints)
- improved command line options for all tools, including file name support for SIESTA, VMEC and XPAND
- new read-in routine for parameterfile
- 3D field line tracing with structure tool now fully supported
- c++ 2011 standard now required -> compiler flag -std=c++1 added in general

MAFOT 3.9 -- July 2017
- laminar can now compute Lc, psimin, etc. along an EFIT flux surface of constant psi
- d3dplot.py can plot the new laminar plots

MAFOT 3.8 -- May 2017
- C-Mod support added (only plot_mpi for now)
- 3D wall added
- Python "== None" warning fixed

MAFOT 3.75 -- Jul 2016
- divB correction in xpand does not work properly - DO NOT USE
- use of multiple M3DC1 files now allows arbitrary relative phase
- fake island perturbation added (D3D only) 

MAFOT 3.74 -- Apr 2016
- splines separated from EFIT_Class into own header
- function added to read & store header lines from files
- bug fix in VMEC class
- structure tool updated
- divB correction added to xpand
- minor fixes in python tools

GUI 1.43 -- Feb 2016
- bug fix in all plot calls for Drop cluster only where it calls the wrong job file
- replaced all ${NSLOTS} with str(nproc) in Drop cluster job files
- add „source /etc/profile.d/modules.sh“ to the mpirun job batch files for Drop as a workaround for the module load bug in bash

MAFOT 3.73 -- Oct 2015
- shell command line flags enabled in all tools, use -h for details

GUI 1.42 -- Oct 2015
- shell command flags can be entered through the "File Tag" input now

MAFOT 3.72 -- Aug 2015
- Bug fix in EFIT_class
- fix now searches on RZ grid as well; for period = 1: only first result within grid is returned

GUI 1.41 -- Aug 2015
- changed fix to use RZ grid in D3D and NSTX

MAFOT 3.71 -- Aug 2015
- in xpand: enable force axisymmetriy in VMECthrough n0only
- in xpand: use full VMEC B-field for virtual casing; interpolation error in vacuum field much larger 
            than any noise from including the vacuum field

MAFOT 3.7 -- Aug 2015
- NSTX-U divertor targets and NSTX outer inclined wall added

GUI 1.4 -- Aug 2015
- support for new NSTX & NSTX-U walls added

MAFOT 3.61 -- July 2015
- readfile now counts rows and allocates total array at once
	-> XFIELD read in works now with very large files
- M3DC1 response > 1 now possible -> reads proper time_xxx.h5

GUI 1.3 -- July 2015
- new entry field added for M3DC1 response_time > 1

MAFOT 3.6 -- July 2015
- M3DC1 interface is now a separate class

MAFOT 3.52 -- July 2015
- Bug fix for M3DC1: psi_n is now calculated correctly from M3DC1 equilibrium

GUI 1.2 -- May 2015
- All NSTX tools added

MAFOT 3.51 -- May 2015
- Bug fix for NSTX: M3DC1 settings now enabled

MAFOT 3.5 -- Feb 2015
- Include the support for external fields from xpand_mpi. This replaces Diagno.
  code looks for file "xpand.dat" ("diagno.dat" no longer recognized, 
  "xpand.dat" can come from the DIAGNO code)

MAFOT 3.41 -- Feb 2015
- BUG fixes in M3DC1 multimode capability for all machines:
    Equilibrium only used from first modefile, not from all
    machinesup.in file no longer required, if m3dc1sup.in is used
    field eval fails are now properly handled

GUI 1.15 -- Oct 2014
- SIESTA support added; plot und laminar only

MAFOT 3.4 -- Oct 2014
- SIESTA support added; code looks for file "siesta.dat"

MAFOT 3.3 -- Oct 2014
- VMEC/DIAGNO support added; code looks for files "diagno.dat" and "wout.nc"
- pragma OMP barrier removed to comply with latest compiler rules

GUI 1.14 -- Oct 2014
- VMEC/DIAGNO support added

GUI 1.131 -- Sep 2014
- g-file path can now be a relative path too

MAFOT 3.2 -- Jul 2014
- M3D-C1 interface updated to use Nate's libfusionio -> resistive wall and multiple M3D-C1 files now enabled
- Tools now unified across machines
- plot is now parallel -> plot_mpi (old, serial plot version still exists)

GUI 1.13 -- Jul 2014
- GUI now uses plot_mpi instead of plot

GUI 1.12 -- Feb 2014
- NSTX foot_mpi added

GUI 1.11 -- Nov 2013
- Bug fixes for ITER

GUI 1.1 -- Nov 2013
- ITER added

GUI 1.0 -- Oct 2013
- first release for D3D MAFOT control

MAFOT 3.1 -- Apr 2012
- M3D-C1 interface added; code looks for file "C1.h5"

MAFOT 3.0 -- Jan 2011
- EFIT, IO and Particle Classes introduced -> replace old includes
- includes unified across machines; the parts that are unique to each machine are in a machine include file

MAFOT 2.0 -- March 2009
- current filaments enabled
- support for multiple machines enabled
- ITER, MAST & NSTX added

MAFOT 1.0 -- May 2008
- first release of D3D MAFOT code
