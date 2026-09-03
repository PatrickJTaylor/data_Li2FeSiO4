# SPDX-FileCopyrightText: Copyright (C) 2026 Benjamin J. Morgan
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Full cation coordination shell of the peroxide-forming O atoms (O52, O80) at 500 K.

For each of the 4 quenched 500K frames, report all Fe and Si neighbours of
O52 and O80 within bonding (2.5 A) and extended (3.0 A) cutoffs.

Atom indexing (0-based): Fe 0-35, O 36-179, Si 180-215 (216 atoms total).
O52 = structure index 52, O80 = structure index 80.
"""

import numpy as np
from pymatgen.io.vasp import Poscar
import os

from result_regeneration.errors import DependencyError
from result_regeneration.paths import EXTRACTED_DATA

# --- Configuration ---
BASE = EXTRACTED_DATA / "molecular_dynamics"

FRAMES = [0, 6881, 16602, 39366]
FRAME_LABELS = ["i", "ii", "iii", "iv"]

TARGET_O = [52, 80]
BONDING_CUTOFF = 2.5   # A
EXTENDED_CUTOFF = 3.0  # A

N_FE = 36
SI_START = 180


def load_frame(frame_idx):
    frame_path = os.path.join(BASE, "500K", "selected_frames", str(frame_idx))
    poscar = Poscar.from_file(os.path.join(frame_path, "POSCAR"))
    structure = poscar.structure
    try:
        ox_states = np.load(os.path.join(frame_path, "oxidation_states.npy"))
    except OSError as base_error:
        raise DependencyError from base_error
    return structure, ox_states


def species_label(idx):
    if idx < N_FE:
        return f"Fe{idx}"
    elif idx < SI_START:
        return f"O{idx}"
    else:
        return f"Si{idx}"


def get_cation_neighbours(structure, o_idx, cutoff):
    """Get all Fe and Si neighbours of an O atom within cutoff."""
    neighbours = structure.get_neighbors(structure[o_idx], r=cutoff)
    cations = []
    for n in neighbours:
        sp = n.species_string
        if sp in ("Fe", "Si"):
            cations.append((n.index, sp, n.nn_distance))
    cations.sort(key=lambda x: x[2])
    return cations


def main():
    print("=" * 75)
    print("R3.3b: Full cation coordination of dimer-forming O atoms")
    print("=" * 75)
    print(f"Target O atoms: O52 (index 52), O80 (index 80)")
    print(f"Bonding cutoff: {BONDING_CUTOFF} A")
    print(f"Extended cutoff: {EXTENDED_CUTOFF} A")
    print(f"Temperature: 500K")

    # Also compute O52-O80 distance in each frame
    print(f"\n{'=' * 75}")
    print("O52-O80 distance across frames")
    print(f"{'=' * 75}")
    for label, step in zip(FRAME_LABELS, FRAMES):
        structure, _ = load_frame(step)
        d = structure.get_distance(52, 80)
        print(f"  Frame {label} (step {step}): O52-O80 = {d:.4f} A")

    # Main analysis: full coordination shell per frame
    for label, step in zip(FRAME_LABELS, FRAMES):
        structure, ox_states = load_frame(step)

        print(f"\n{'=' * 75}")
        print(f"  Frame {label} (step {step})")
        print(f"{'=' * 75}")

        fe_ox = ox_states[:N_FE]

        # Fe oxidation state summary
        unique, counts = np.unique(fe_ox.astype(int), return_counts=True)
        ox_str = ", ".join(f"Fe^{u}: {c}" for u, c in zip(unique, counts))
        print(f"  Fe oxidation states: {ox_str}")

        for o_idx in TARGET_O:
            o_ox = int(ox_states[o_idx])
            print(f"\n  --- O{o_idx} (O oxidation state: {o_ox:+d}) ---")

            # Bonding shell
            bonding = get_cation_neighbours(structure, o_idx, BONDING_CUTOFF)
            # Extended shell (beyond bonding)
            extended_all = get_cation_neighbours(structure, o_idx, EXTENDED_CUTOFF)
            extended_only = [c for c in extended_all if c[2] > BONDING_CUTOFF]

            print(f"  Cation neighbours within {BONDING_CUTOFF} A:")
            if bonding:
                n_fe = sum(1 for _, sp, _ in bonding if sp == "Fe")
                n_si = sum(1 for _, sp, _ in bonding if sp == "Si")
                print(f"    Count: {n_fe} Fe + {n_si} Si = {len(bonding)} total")
                for idx, sp, d in bonding:
                    if sp == "Fe":
                        ox = int(fe_ox[idx])
                        print(f"    {species_label(idx):>6}  d = {d:.4f} A  "
                              f"(Fe^{ox})")
                    else:
                        print(f"    {species_label(idx):>6}  d = {d:.4f} A")
            else:
                print(f"    None")

            if extended_only:
                print(f"  Additional cations in {BONDING_CUTOFF}-{EXTENDED_CUTOFF} A:")
                for idx, sp, d in extended_only:
                    if sp == "Fe":
                        ox = int(fe_ox[idx])
                        print(f"    {species_label(idx):>6}  d = {d:.4f} A  "
                              f"(Fe^{ox})")
                    else:
                        print(f"    {species_label(idx):>6}  d = {d:.4f} A")
            else:
                print(f"  Additional cations in {BONDING_CUTOFF}-{EXTENDED_CUTOFF} A: "
                      f"none")

    # =====================================================================
    # Summary: track how coordination changes across frames
    # =====================================================================
    print(f"\n{'=' * 75}")
    print("Summary: cation coordination evolution")
    print(f"{'=' * 75}")

    for o_idx in TARGET_O:
        print(f"\n  O{o_idx} coordination (within {BONDING_CUTOFF} A):")
        header = f"    {'Neighbour':<10}"
        for label in FRAME_LABELS:
            header += f"{'Frame ' + label:>16}"
        print(header)
        print("    " + "-" * (10 + 16 * len(FRAME_LABELS)))

        # Collect all cation indices that appear in any frame
        all_cations = set()
        frame_data = {}
        for label, step in zip(FRAME_LABELS, FRAMES):
            structure, ox_states = load_frame(step)
            bonding = get_cation_neighbours(structure, o_idx, BONDING_CUTOFF)
            frame_data[label] = {idx: (sp, d, ox_states) for idx, sp, d in bonding}
            for idx, sp, d in bonding:
                all_cations.add((idx, sp))

        # Also check extended shell for completeness
        for label, step in zip(FRAME_LABELS, FRAMES):
            structure, ox_states = load_frame(step)
            extended = get_cation_neighbours(structure, o_idx, EXTENDED_CUTOFF)
            for idx, sp, d in extended:
                all_cations.add((idx, sp))

        # Sort: Fe first, then Si, by index
        cation_list = sorted(all_cations, key=lambda x: (0 if x[1] == "Fe" else 1, x[0]))

        for idx, sp in cation_list:
            row = f"    {species_label(idx):<10}"
            for label, step in zip(FRAME_LABELS, FRAMES):
                structure, ox_states = load_frame(step)
                # Get distance to this specific cation
                d = structure.get_distance(o_idx, idx)
                if sp == "Fe":
                    ox = int(ox_states[idx])
                    if d <= BONDING_CUTOFF:
                        row += f"  {d:.3f} Fe^{ox:d}"
                    elif d <= EXTENDED_CUTOFF:
                        row += f" ({d:.3f} Fe^{ox:d})"
                    else:
                        row += f"        --     "
                else:
                    if d <= BONDING_CUTOFF:
                        row += f"  {d:.3f}      "
                    elif d <= EXTENDED_CUTOFF:
                        row += f" ({d:.3f})     "
                    else:
                        row += f"        --     "
            print(row)

        print(f"    (values in parentheses are in the {BONDING_CUTOFF}-{EXTENDED_CUTOFF} A "
              f"extended shell)")


if __name__ == "__main__":
    main()
