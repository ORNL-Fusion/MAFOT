"""Minimal Python wrapper for the persistent MAFOT C++ map."""

import ctypes
from dataclasses import dataclass
import os
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class MAFOTMapConfig:
    run_directory: Path
    gfile: Path
    cpp_library: Path
    phistart: float = 0.0
    field: int = -1
    plasma_response: int = 0
    use_fcoil: bool = False
    use_ccoil: bool = False
    use_icoil: bool = False
    use_filament: bool = False
    use_bus_error: bool = False
    use_bcoil_error: bool = False
    particle_direction: int = 0
    particle_charge: int = 1
    energy_kev: float = 100.0
    pitch: float = 0.1
    mass: int = 2
    temperature_profile: Path = Path("prof_t.dat")
    density_profile: Path = Path("prof_n.dat")
    boundary1: Path | None = None
    boundary2: Path | None = None
    collisions: bool = False
    f: float = 0.0
    zbar: float = 2.0
    rc: float = 1.0
    mc: float = 1.0


class MAFOTMap:
    def __init__(self, config: MAFOTMapConfig) -> None:
        self.config = config
        self._write_iter_dat()
        self._library = ctypes.CDLL(os.fspath(config.cpp_library))

        self._create = self._library.mafot_engine_create
        self._create.argtypes = [
            ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p,
            ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
            ctypes.c_char_p, ctypes.c_char_p, ctypes.c_bool,
        ]
        self._create.restype = ctypes.c_void_p

        self._map_once = self._library.mafot_engine_map_once
        self._map_once.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_int,
            ctypes.POINTER(ctypes.c_double),
        ]
        self._map_once.restype = ctypes.c_int

        self._destroy = self._library.mafot_engine_destroy
        self._destroy.argtypes = [ctypes.c_void_p]
        self._destroy.restype = None

        self._last_error = self._library.mafot_engine_last_error
        self._last_error.argtypes = []
        self._last_error.restype = ctypes.c_char_p

        def path_arg(path: Path | None) -> bytes:
            return b"" if path is None else os.fsencode(config.run_directory / path)

        self._handle = self._create(
            os.fsencode(config.run_directory),
            path_arg(config.temperature_profile),
            path_arg(config.density_profile),
            config.f, config.zbar, config.rc, config.mc,
            path_arg(config.boundary1), path_arg(config.boundary2),
            config.collisions,
        )
        if not self._handle:
            error = self._last_error()
            raise RuntimeError(error.decode(errors="replace") if error else "Could not create MAFOT engine")

    def _write_iter_dat(self) -> None:
        c = self.config
        gfile = c.gfile if c.gfile.is_absolute() else c.run_directory / c.gfile
        (c.run_directory / "_iter.dat").write_text(
            f"""# Parameterfile for DIII-D Programs
# Path: {gfile}
free_parameter=0
itt=0
rmin=0
rmax=0
thmin=0
thmax=0
N=1
phistart(deg)={c.phistart}
MapDirection=1
PlasmaResponse(0=no,>1=yes)={c.plasma_response}
Field={c.field}
target=0
createPoints=0
useFcoil={int(c.use_fcoil)}
useCcoil={int(c.use_ccoil)}
useIcoil={int(c.use_icoil)}
ParticleDirection={c.particle_direction}
ParticleCharge={c.particle_charge}
Ekin[keV]={c.energy_kev}
lambda={c.pitch}
Mass={c.mass}
useFilament={int(c.use_filament)}
useBusError={int(c.use_bus_error)}
useBcoilError={int(c.use_bcoil_error)}
"""
        )

    def map_once(self, point: np.ndarray, direction: int) -> np.ndarray:
        if self._handle is None:
            raise RuntimeError("The MAFOT engine is closed")
        r, z = point
        output = np.empty(6, dtype=np.float64)
        status = self._map_once(
            self._handle, float(r), float(z), direction,
            output.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        )
        if status:
            error = self._last_error()
            raise RuntimeError(error.decode(errors="replace") if error else "MAFOT map failed")
        return output[:2].copy()

    def forward(self, point: np.ndarray) -> np.ndarray:
        return self.map_once(point, 1)

    def inverse(self, point: np.ndarray) -> np.ndarray:
        return self.map_once(point, -1)

    def poincare(
        self, psi_min: float, psi_max: float, surfaces: int = 16, turns: int = 300
    ) -> tuple[np.ndarray, np.ndarray]:
        rz, psi_theta = [], []
        for psi in np.linspace(psi_min, psi_max, surfaces):
            point = self.map_once(np.array((psi, 0.0)), 0)
            for _ in range(turns):
                output = np.empty(6)
                if self._map_once(
                    self._handle, *point, 1,
                    output.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                ):
                    break
                point = output[:2].copy()
                rz.append(point)
                psi_theta.append(output[2:4].copy())
        return np.array(rz), np.array(psi_theta)

    def close(self) -> None:
        if self._handle is not None:
            self._destroy(self._handle)
            self._handle = None
