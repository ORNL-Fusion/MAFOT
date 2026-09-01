"""Visual smoke test for the persistent MAFOT map."""

from pathlib import Path
import sys

import matplotlib.pyplot as plt

from mafot_map import MAFOTMap, MAFOTMapConfig


MAFOT_ROOT = Path(__file__).resolve().parents[1]
LIBRARY_NAME = "libmafot_cpp_bridge.dylib" if sys.platform == "darwin" else "libmafot_cpp_bridge.so"


def main() -> None:
    mafot = MAFOTMap(MAFOTMapConfig(
        run_directory=MAFOT_ROOT / "examples/testcase",
        gfile=Path("g126006.03600"),
        cpp_library=MAFOT_ROOT / "build/lib" / LIBRARY_NAME,
        use_icoil=True,
    ))

    rz, psi_theta = mafot.poincare(0.8, 1.0)
    mafot.close()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    ax1.plot(rz[:, 0], rz[:, 1], "k.", markersize=0.5)
    ax1.set_xlabel("R (m)")
    ax1.set_ylabel("Z (m)")
    ax2.plot(psi_theta[:, 1], psi_theta[:, 0], "k.", markersize=0.5)
    ax2.set_xlabel(r"$\theta$")
    ax2.set_ylabel(r"$\psi$")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
