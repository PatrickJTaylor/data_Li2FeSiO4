# SPDX-FileCopyrightText: Copyright (C) 2026 Benjamin J. Morgan
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Extract the total energies of the relaxed (quenched) AIMD frames from the raw VASP
output, and write them to extracted_data/molecular_dynamics/<T>/relaxed_frame_energies.npy.

This script documents how the archived energies were obtained. It requires the raw
dataset (University of Bath Research Data Archive, doi:10.15125/BATH-01647), which is
not included in this repository.

Each selected frame was quenched by geometry optimisation (relaxation/), then a static
calculation was run at the relaxed geometry (lobster/: same cutoff and k-points,
tetrahedron-method Brillouin-zone integration, any O-O dimer atoms initialised with a
small negative magnetic moment). The archived energies, magnetic moments, ICOBIs and
Wannier oxidation states for each frame all come from this static calculation, so the
energy taken here is its TOTEN (for the tetrahedron method this equals energy(sigma->0)).

Usage:
    python scripts/extract_relaxed_frame_energies.py /path/to/raw/molecular_dynamics
"""

import argparse
import re
from pathlib import Path

import numpy as np

from result_regeneration.paths import EXTRACTED_DATA

FRAME_SETS = {"500K": (0, 6881, 16602, 39366), "1000K": (0, 2581, 8315, 12276)}
ENERGY_RE = re.compile(r"free  energy   TOTEN\s*=\s*(-?\d+\.\d+)")


def final_energy(outcar: Path) -> float:
    """Return the final TOTEN in an OUTCAR."""
    matches = ENERGY_RE.findall(outcar.read_text())
    if not matches:
        raise ValueError(f"No energies found in {outcar}")

    return float(matches[-1])


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "raw_data",
        type=Path,
        help="path to the raw molecular_dynamics directory (containing 500K/ and 1000K/)",
    )
    args = parser.parse_args()

    for temperature, frames in FRAME_SETS.items():
        energies = np.array(
            [
                final_energy(
                    args.raw_data
                    / temperature
                    / "selected_frames"
                    / str(frame)
                    / "lobster"
                    / "OUTCAR"
                )
                for frame in frames
            ]
        )

        output = EXTRACTED_DATA / "molecular_dynamics" / temperature
        np.save(output / "relaxed_frame_energies.npy", energies)

        print(f"{temperature}: {energies} -> {output / 'relaxed_frame_energies.npy'}")


if __name__ == "__main__":
    main()
