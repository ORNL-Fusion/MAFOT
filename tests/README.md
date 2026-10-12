# MAFOT test data & procedures

This directory holds a small DIII-D test case (shot 126006 @ 3600 ms) plus control
files for regression-testing the CPU/MPI tools and the GPU path. There is no
automated test runner; the procedures below are run by hand (or inside the HEAT
builder docker image) and compared as described.

## Files

| File | Purpose |
|---|---|
| `g126006.03600` | EFIT g-file (129x129), DIII-D shot 126006 @ 3600 ms |
| `diiidsup.in` | DIII-D supplemental coil currents (needed by the `dt*` tools) |
| `_plot.dat` | Control file for `dtplot`: 40-line Poincaré plot, 300 turns, I-coil on |
| `plot_ref.dat` | Reference `dtplot` output (upstream, 2022) for `_plot.dat` |
| `points_test.dat` | Reference points output belonging to the `_plot.dat` case |
| `wout_akima.126006.03600_fixed.nc` | VMEC wout file for xpand tests (`xpand_mpi`) |
| `_struct_HEAT.dat` | Control file for `heatstructure` (axisymmetric, MapDirection=0, dpinit=1) |
| `points_HEAT.dat` | 5 start points (R [m], phi [deg], Z [m]) for `_struct_HEAT.dat` |
| `_lam_HEAT.dat` | Control file for `heatlaminar_mpi`: 30x30 R-Z grid, 150 turns, MapDirection=0 |

## 1. Poincaré regression vs. reference (`dtplot`)

Validates the equilibrium reader, splines, RK4 integrator, and the trip3d coil
library against a known-good 2022 output. Requires the DIII-D build (`make d3d`).

```bash
cd tests/
dtplot _plot.dat test        # writes plot_test.dat (~16 s)
```

Chaotic field lines diverge exponentially, so outputs are NOT bitwise identical
across compilers/hardware. Instead check that (a) the first ~200 rows agree to
high precision and (b) the column statistics agree closely:

```python
import numpy as np
d = np.genfromtxt('plot_test.dat', comments='#')
r = np.genfromtxt('plot_ref.dat', comments='#')
n = 200
assert np.abs((d[:n]-r[:n])/(np.abs(r[:n])+1e-12)).max() < 1e-6   # pre-divergence
for c in range(d.shape[1]):                                        # statistics
    assert abs(d[:,c].mean()-r[:,c].mean()) < 5e-2*max(1.0, abs(r[:,c].mean()))
```

Typical good result: first-200-row max relative difference ~1e-8; means agree
to ~0.3 %.

## 2. HEAT structure — CPU vs. GPU parity (`heatstructure`)

Requires the HEAT build (`make GPU=True heat`) and a CUDA-capable GPU for the
`-g` run. Axisymmetric field lines must stay on their flux surface, so the psi
column is the key check (see also the psi-conservation note below).

```bash
cd tests/
heatstructure    -P points_HEAT.dat _struct_HEAT.dat cpu   # struct_cpu.dat
heatstructure -g -P points_HEAT.dat _struct_HEAT.dat gpu   # struct_gpu.dat
```

Compare the 8 columns (X, Y, Z, R, phi, psi, Lc, dpsidLc) row-wise; the GPU uses
trilinear interpolation on a 300x300 sampled grid instead of bicubic splines, so
expect small (~1e-3 psi) but not zero differences. Structural checks that must
hold exactly:

- Same row ordering: backward trace reversed, then start point, then forward.
- With MapDirection=0 the Lc column is continuous and increasing through the
  start point (the start row carries the backward-total Lc, as on the CPU).
- psi along each axisymmetric trace stays ~constant on both CPU and GPU.

### 2b. Drift-term parity (newest GPU kernel code)

The relativistic guiding-center drift (`sigma != 0`) and the prescribed filament
drift velocities (control-file entries 26-28) are separate code paths in the GPU
kernel; test them the same way with modified control files:

```bash
# guiding-center drift: co-passing 10 keV ions
sed -e 's/^ParticleDirection.*/ParticleDirection(1=co-pass,-1=ctr-pass,0=field-lines)=\t1/' \
    _struct_HEAT.dat > _struct_drift.dat
heatstructure    -P points_HEAT.dat _struct_drift.dat dcpu
heatstructure -g -P points_HEAT.dat _struct_drift.dat dgpu

# filament drift velocities
sed -e 's/^v_parallel.*/v_parallel=\t100000.0/' -e 's/^v_radial.*/v_radial=\t500.0/' \
    _struct_HEAT.dat > _struct_vrad.dat
heatstructure    -P points_HEAT.dat _struct_vrad.dat vcpu
heatstructure -g -P points_HEAT.dat _struct_vrad.dat vgpu
```

Expected agreement (2026-07 baseline): identical row counts; position columns
within ~1e-3 m, psi within ~1e-4, Lc within ~1e-3 km — the same
trilinear-vs-bicubic interpolation floor as the pure-field-line test.

## 3. HEAT laminar — CPU vs. GPU parity (`heatlaminar_mpi`)

```bash
cd tests/
mpirun -n 4 heatlaminar_mpi    _lam_HEAT.dat cpu   # CPU/MPI path (~30 s on 4 ranks)
mpirun -n 1 heatlaminar_mpi -g _lam_HEAT.dat gpu   # GPU path (rank 0 only)
```

Both write an 11-column map over the same 30x30 grid (R, Z, N_toroidal,
Lc [km], psimin, psiav, B_R, B_Z, B_phi, theta, psi). Compare column-wise;
Lc/psimin of chaotic edge points differ more than core points. The B_* and
theta/psi start-point columns should agree to interpolation accuracy.

## 4. Structure no-pointfile mode sanity

Without `-P`, structure builds start points on a line from (Rmin,Zmin) to
(Rmax,Zmax) (control-file entries; set entry 0 "Nphi/free_Parameter" to 0 so the
parabolic deformation `verschieb` is off). Verify the rows with Lc == 0 sit
exactly on that line. (This mode silently produced garbage Z values before the
`initial` array indexing fix in structure.cxx.)

## Known coverage gaps

- **3-D field path untested**: every GPU test above is axisymmetric (`Field = -1`,
  single phi plane). The 3-D interpolation branch (Nphi=64 grid, phi periodicity
  wrap in `interp3D`) is not exercised by any test here. The cleanest way to close
  this: generate an `xpand.dat` from `wout_akima.126006.03600_fixed.nc` with
  `xpand_mpi`, then run `heatstructure`/`heatlaminar_mpi` with
  `response_field = -3` CPU vs GPU (requires `VMEC = True`).
- E-field / sheath drift terms are intentionally NOT on the GPU (guarded: structure
  falls back to CPU, laminar errors) — the guard behavior itself is easy to check
  but the terms have no GPU parity test by design.
- `foot_mpi`, `plot_mpi`, `trace`, `xpand_mpi` MPI ports are compile-verified only;
  the runtime MPI test coverage comes from `heatlaminar_mpi` (same master/slave
  pattern replicated across all `_mpi` tools).

## GPU psi-conservation smoke test

The fastest GPU sanity check (catches interpolation/indexing regressions in the
kernel): trace a single point one full turn with `heatstructure -g` on an
axisymmetric equilibrium (`Field = -1`) and confirm the psi column is constant
to ~1e-3 and the output has the expected number of rows. On the CPU the same
trace conserves psi to ~1e-6.

## Running inside the HEAT builder image

The HEAT docker image builds MAFOT with GPU support; to test local (unpushed)
changes, mount the repo into a `heatbuilder` container, copy `src/ include/
makefile install/` over `/root/source/MAFOT/`, and rebuild:

```bash
cd /path/to/HEAT && docker build -f github/docker/Dockerfile --target heatbuilder -t heatbuilder:local github
docker run -d --name mafot_test --gpus all -v /path/to/MAFOT:/mnt/MAFOT heatbuilder:local sleep infinity
docker exec mafot_test bash -c "cp -r /mnt/MAFOT/{src,include,makefile,install} /root/source/MAFOT/ && \
    cd /root/source/MAFOT && make GPU=True heat && make dtplot dtstructure dtlaminar_mpi"
```

MPI runs inside the container need `--allow-run-as-root` (and `--oversubscribe`
if ranks exceed cores). If `-g` reports "no CUDA-capable GPU" although
`nvidia-smi` works, the host CUDA/UVM state is stale (common after
suspend/resume): close GPU apps and reload with
`sudo modprobe -r nvidia_uvm && sudo modprobe nvidia_uvm`.
