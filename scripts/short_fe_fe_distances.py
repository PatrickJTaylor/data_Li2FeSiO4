# SPDX-FileCopyrightText: Copyright (C) 2026 Benjamin J. Morgan
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Cross-reference short Fe-Fe distances with coordination environments
and oxidation states.

For each quenched AIMD frame, compute all 630 pairwise Fe-Fe distances
(minimum image convention), identify pairs < 3.2 A, and report the oxidation
state and coordination number of each Fe in the pair.

Atom indexing: Fe 0-35, O 36-179, Si 180-215 (216 atoms total).
"""

import numpy as np
from pymatgen.io.vasp import Poscar
from collections import defaultdict
import os

from result_regeneration.errors import DependencyError
from result_regeneration.paths import EXTRACTED_DATA

# --- Configuration ---
BASE = EXTRACTED_DATA / "molecular_dynamics"

TEMPERATURES = ("500K", "1000K")
FRAME_SETS = {
    "500K": (0, 6881, 16602, 39366),
    "1000K": (0, 2581, 8315, 12276),
}
FRAME_LABELS = ("i", "ii", "iii", "iv")

N_FE = 36
FE_INDICES = range(N_FE)
FE_O_CUTOFF = 2.4       # A, for Fe-O coordination number
FE_FE_CUTOFF = 3.2      # A, generous cutoff for short Fe-Fe pairs
FE2O3_NN = 2.97         # A, nearest-neighbour Fe-Fe in corundum Fe2O3


def load_frame(temperature, frame_idx):
    """Load structure and oxidation states for a frame."""
    frame_path = os.path.join(BASE, temperature, "selected_frames", str(frame_idx))
    poscar = Poscar.from_file(os.path.join(frame_path, "POSCAR"))
    structure = poscar.structure
    try:
        ox_states = np.load(os.path.join(frame_path, "oxidation_states.npy"))
    except OSError as base_error:
        raise DependencyError from base_error
    return structure, ox_states


def get_fe_cn(structure, cutoff=FE_O_CUTOFF):
    """Compute Fe-O coordination number for each Fe atom."""
    cn_list = []
    for fe_idx in FE_INDICES:
        neighbors = structure.get_neighbors(structure[fe_idx], r=cutoff)
        n_o = sum(1 for n in neighbors if n.species_string == "O")
        cn_list.append(n_o)
    return cn_list


def classify_cn(cn):
    """Short label for coordination number."""
    labels = {3: "3", 4: "4(tet)", 5: "5", 6: "6(oct)", 7: "7", 8: "8"}
    return labels.get(cn, str(cn))


def get_short_fefe_pairs(structure, cutoff=FE_FE_CUTOFF):
    """
    Compute all pairwise Fe-Fe distances using minimum image convention.
    Return list of (i, j, distance) for pairs with distance < cutoff,
    sorted by distance.
    """
    pairs = []
    for i in range(N_FE):
        for j in range(i + 1, N_FE):
            d = structure.get_distance(i, j)
            if d < cutoff:
                pairs.append((i, j, d))
    pairs.sort(key=lambda x: x[2])
    return pairs


def ox_pair_label(ox_a, ox_b):
    """Canonical label for an oxidation-state pair (sorted so Fe^III-Fe^IV, not Fe^IV-Fe^III)."""
    lo, hi = sorted([int(ox_a), int(ox_b)])
    return f"Fe^{lo}-Fe^{hi}"


def main():
    print("=" * 80)
    print("R3.2b: Short Fe-Fe Distances vs Coordination and Oxidation State")
    print("=" * 80)
    print(f"Fe-Fe cutoff: {FE_FE_CUTOFF} A")
    print(f"Fe-O coordination cutoff: {FE_O_CUTOFF} A")
    print(f"Fe2O3 nearest-neighbour reference: {FE2O3_NN} A")

    for temperature in TEMPERATURES:
        frames = FRAME_SETS[temperature]
        print(f"\n{'#' * 80}")
        print(f"# {temperature}")
        print(f"{'#' * 80}")

        for frame_num, frame_idx in enumerate(frames):
            label = FRAME_LABELS[frame_num]
            structure, ox_states = load_frame(temperature, frame_idx)
            fe_ox = ox_states[:N_FE]
            fe_cn = get_fe_cn(structure)

            # Ox state summary
            unique, counts = np.unique(fe_ox.astype(int), return_counts=True)
            ox_summary = ", ".join(f"Fe^{u}: {c}" for u, c in zip(unique, counts))

            # CN summary
            cn_unique, cn_counts = np.unique(fe_cn, return_counts=True)
            cn_summary = ", ".join(
                f"CN={c}: {n}" for c, n in zip(cn_unique, cn_counts)
            )

            print(f"\n{'=' * 70}")
            print(f"  Frame {label} (step {frame_idx}), {temperature}")
            print(f"{'=' * 70}")
            print(f"  Fe oxidation states: {ox_summary}")
            print(f"  Fe coordination numbers: {cn_summary}")

            # Find short Fe-Fe pairs
            short_pairs = get_short_fefe_pairs(structure)

            if not short_pairs:
                print(f"  Short Fe-Fe pairs (< {FE_FE_CUTOFF} A): none")
                # Report the shortest Fe-Fe distance for reference
                min_d = float("inf")
                for i in range(N_FE):
                    for j in range(i + 1, N_FE):
                        d = structure.get_distance(i, j)
                        if d < min_d:
                            min_d = d
                print(f"  Shortest Fe-Fe distance: {min_d:.3f} A")
                continue

            print(f"  Short Fe-Fe pairs (< {FE_FE_CUTOFF} A): {len(short_pairs)}")
            print()

            # Detailed table
            print(f"    {'Pair':<12} {'d(A)':>7} {'Fe_a ox':>8} {'CN_a':>8} "
                  f"{'Fe_b ox':>8} {'CN_b':>8} {'Ox pair':>14}")
            print("    " + "-" * 70)

            ox_pair_counts = defaultdict(int)
            cn_pair_counts = defaultdict(int)
            both_non_tet = 0
            at_least_one_non_tet = 0

            for fe_a, fe_b, d in short_pairs:
                ox_a = int(fe_ox[fe_a])
                ox_b = int(fe_ox[fe_b])
                cn_a = fe_cn[fe_a]
                cn_b = fe_cn[fe_b]
                pair_label = ox_pair_label(ox_a, ox_b)

                cn_a_str = classify_cn(cn_a)
                cn_b_str = classify_cn(cn_b)

                print(f"    Fe{fe_a:<2}-Fe{fe_b:<2}  {d:>7.3f} "
                      f"{'Fe^' + str(ox_a):>8} {cn_a_str:>8} "
                      f"{'Fe^' + str(ox_b):>8} {cn_b_str:>8} "
                      f"{pair_label:>14}")

                ox_pair_counts[pair_label] += 1

                # CN pair (sorted)
                cn_lo, cn_hi = sorted([cn_a, cn_b])
                cn_pair_label = f"CN{cn_lo}-CN{cn_hi}"
                cn_pair_counts[cn_pair_label] += 1

                if cn_a != 4 and cn_b != 4:
                    both_non_tet += 1
                if cn_a != 4 or cn_b != 4:
                    at_least_one_non_tet += 1

            # Summary tables
            n_pairs = len(short_pairs)
            print(f"\n  Oxidation state pair counts:")
            for pair_label in sorted(ox_pair_counts.keys()):
                c = ox_pair_counts[pair_label]
                print(f"    {pair_label:<16} {c:>4}  ({100*c/n_pairs:.0f}%)")
            print(f"    {'Total':<16} {n_pairs:>4}")

            print(f"\n  Coordination number pair counts:")
            for cn_label in sorted(cn_pair_counts.keys()):
                c = cn_pair_counts[cn_label]
                print(f"    {cn_label:<16} {c:>4}  ({100*c/n_pairs:.0f}%)")
            print(f"    {'Total':<16} {n_pairs:>4}")

            print(f"\n  Both Fe non-tetrahedral: {both_non_tet}/{n_pairs} "
                  f"({100*both_non_tet/n_pairs:.0f}%)")
            print(f"  At least one non-tetrahedral: {at_least_one_non_tet}/{n_pairs} "
                  f"({100*at_least_one_non_tet/n_pairs:.0f}%)")

            # Check for shared O neighbours (face-sharing / edge-sharing)
            print(f"\n  Shared O-neighbour analysis (Fe-O cutoff {FE_O_CUTOFF} A):")
            for fe_a, fe_b, d in short_pairs:
                nbrs_a = set()
                for n in structure.get_neighbors(structure[fe_a], r=FE_O_CUTOFF):
                    if n.species_string == "O":
                        nbrs_a.add(n.index)
                nbrs_b = set()
                for n in structure.get_neighbors(structure[fe_b], r=FE_O_CUTOFF):
                    if n.species_string == "O":
                        nbrs_b.add(n.index)
                shared = nbrs_a & nbrs_b
                n_shared = len(shared)
                if n_shared >= 3:
                    sharing = "face-sharing"
                elif n_shared == 2:
                    sharing = "edge-sharing"
                elif n_shared == 1:
                    sharing = "corner-sharing"
                else:
                    sharing = "no shared O"
                shared_str = ", ".join(f"O{idx}" for idx in sorted(shared))
                print(f"    Fe{fe_a}-Fe{fe_b} (d={d:.3f} A): "
                      f"{n_shared} shared O [{sharing}]"
                      + (f" ({shared_str})" if shared else ""))


if __name__ == "__main__":
    main()
