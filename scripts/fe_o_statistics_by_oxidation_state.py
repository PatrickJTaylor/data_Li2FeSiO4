# SPDX-FileCopyrightText: Copyright (C) 2026 Benjamin J. Morgan
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Mean Fe-O ICOBI and mean Fe-O bond length for each Wannier-assigned Fe oxidation
state, pooled over the relaxed AIMD frames.

These are the summary statistics quoted in the Supporting Information as independent
validation of the Wannier oxidation-state assignments, and they correspond to the
distributions plotted in Fig. S3.

Convention (matching Fig. S3): for each Fe, the mean is taken over that Fe's own Fe-O
bonds; the reported value is then the mean over all Fe of a given oxidation state,
across all eight relaxed frames.  Bonds are those listed in ICOBILIST.lobster, so the
bond set is identical for the ICOBI and bond-length statistics.
"""

import numpy as np

from result_regeneration.errors import DependencyError
from result_regeneration.paths import EXTRACTED_DATA

TEMPERATURES = ("500K", "1000K")
FRAME_SETS = ((0, 6881, 16602, 39366), (0, 2581, 8315, 12276))
N_FE = 36
OXIDATION_STATES = (3, 4, 5)
NUMERALS = {3: "III", 4: "IV", 5: "V"}


def fe_o_bonds(path):
    """
    Parse Fe-O bonds from an ICOBILIST.lobster file.

    Returns
    -------
    bonds:
        {(fe_index, oxygen_label): (distance, icobi_summed_over_spins)}
    """
    lines = path.read_text().splitlines()
    headings = [i for i, line in enumerate(lines) if "COBI#" in line]

    bonds: dict[tuple[int, str], list[float]] = {}
    for start, end in zip(headings, headings[1:] + [len(lines)]):
        for line in lines[start + 1 : end]:
            fields = line.split()
            if len(fields) < 8:
                continue
            atom_i, atom_j = fields[1], fields[2]
            if not (atom_i.startswith("Fe") and atom_j.startswith("O")):
                continue
            key = (int(atom_i[2:]) - 1, atom_j)
            if key not in bonds:
                bonds[key] = [float(fields[3]), 0.0]
            bonds[key][1] += float(fields[7])

    return {k: tuple(v) for k, v in bonds.items()}


def main() -> None:
    prefix = EXTRACTED_DATA / "molecular_dynamics"

    icobis = {state: [] for state in OXIDATION_STATES}
    lengths = {state: [] for state in OXIDATION_STATES}

    for temperature, frames in zip(TEMPERATURES, FRAME_SETS):
        for frame_idx in frames:
            frame_dir = prefix / temperature / "selected_frames" / str(frame_idx)

            try:
                oxidation_states = np.load(frame_dir / "oxidation_states.npy")
            except OSError as base_error:
                raise DependencyError from base_error

            bonds = fe_o_bonds(frame_dir / "ICOBILIST.lobster")

            per_fe: dict[int, list[tuple[float, float]]] = {i: [] for i in range(N_FE)}
            for (fe_idx, _), (distance, icobi) in bonds.items():
                per_fe[fe_idx].append((distance, icobi))

            for fe_idx, fe_bonds in per_fe.items():
                if not fe_bonds:
                    continue
                state = int(oxidation_states[fe_idx])
                distances, values = zip(*fe_bonds)
                lengths[state].append(np.mean(distances))
                icobis[state].append(np.mean(values))

    print("Mean Fe-O ICOBI and Fe-O bond length by Wannier-assigned oxidation state")
    print("(per-Fe means, pooled over all eight relaxed AIMD frames)\n")
    print(f"  {'':<8} {'n(Fe)':>6} {'mean ICOBI':>12} {'mean Fe-O / A':>16}")
    for state in OXIDATION_STATES:
        print(
            f"  Fe({NUMERALS[state]:<3}) {len(icobis[state]):>6} "
            f"{np.mean(icobis[state]):>12.3f} {np.mean(lengths[state]):>16.3f}"
        )


if __name__ == "__main__":
    main()
