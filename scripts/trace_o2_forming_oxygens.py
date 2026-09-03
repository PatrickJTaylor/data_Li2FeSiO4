# SPDX-FileCopyrightText: Copyright (C) 2026 Benjamin J. Morgan
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Trace the O atoms that form O-O dimers through the relaxed AIMD frames.

For each trajectory:
  - Scans all relaxed frames for short O-O contacts and classifies them (peroxide / superoxide / O2).
  - Traces the Fe coordination of dimer-forming O atoms back through the preceding frames.
  - Compares the oxidation states of dimer-parent Fe against the background distribution.
  - Checks the 500 K trajectory for precursor O-O shortening before dimer formation.

Supports the O-O dimer and Fe-O distance statements in Sec. III C and III D of the manuscript.

Atom indexing: Fe 0-35, O 36-179, Si 180-215 (216 atoms total).
LOBSTER uses 1-based indexing: Fe1 = index 0, O37 = index 36, etc.
"""

import numpy as np
from pymatgen.io.vasp import Poscar
from collections import defaultdict
import os

from result_regeneration.errors import DependencyError
from result_regeneration.paths import EXTRACTED_DATA

# --- Configuration ---
BASE = EXTRACTED_DATA / "molecular_dynamics"

N_FE = 36
N_O = 144
FE_INDICES = range(N_FE)
O_INDICES = range(36, 180)

FE_O_CUTOFF = 2.4  # A, for Fe-O coordination
OO_DIMER_CUTOFF = 1.6  # A, for identifying O-O dimers (peroxide ~1.4, O2 ~1.2)
OO_SHORT_CUTOFF = 2.0  # A, for identifying precursor short O-O contacts

FRAME_SETS = {
    "500K": [0, 6881, 16602, 39366],
    "1000K": [0, 2581, 8315, 12276],
}
FRAME_LABELS = ("i", "ii", "iii", "iv")


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


def find_oo_contacts(structure, max_dist):
    """
    Find all O-O pairs with distance < max_dist.
    Returns list of (o_idx_1, o_idx_2, distance) sorted by distance.
    """
    contacts = []
    o_list = list(O_INDICES)
    for i_pos, i_idx in enumerate(o_list):
        for j_pos in range(i_pos + 1, len(o_list)):
            j_idx = o_list[j_pos]
            d = structure.get_distance(i_idx, j_idx)
            if d < max_dist:
                contacts.append((i_idx, j_idx, d))
    contacts.sort(key=lambda x: x[2])
    return contacts


def find_fe_neighbors_of_o(structure, o_idx, cutoff=FE_O_CUTOFF):
    """Find Fe atoms coordinating a given O atom within cutoff."""
    neighbors = structure.get_neighbors(structure[o_idx], r=cutoff)
    fe_neighbors = [
        (n.index, n.nn_distance)
        for n in neighbors
        if n.species_string == "Fe"
    ]
    fe_neighbors.sort(key=lambda x: x[1])
    return fe_neighbors


def classify_oo(dist):
    """Classify O-O contact type based on distance."""
    if dist < 1.3:
        return "molecular O2 (~1.21 A)"
    elif dist < 1.55:
        return "peroxide (~1.4-1.5 A)"
    elif dist < 1.7:
        return "possible superoxide/long peroxide"
    else:
        return "short contact (pre-dimer?)"


def fe_ox_distribution(ox_states):
    """Return dict of Fe oxidation state counts."""
    fe_ox = ox_states[:N_FE].astype(int)
    unique, counts = np.unique(fe_ox, return_counts=True)
    return dict(zip(unique.tolist(), counts.tolist()))


def print_comparison(label, dimer_fe_ox_counts, background_dist, n_fe_total=N_FE):
    """
    Print side-by-side comparison of dimer-parent Fe ox states vs background.
    """
    all_ox = sorted(set(list(dimer_fe_ox_counts.keys()) + list(background_dist.keys())))
    n_dimer = sum(dimer_fe_ox_counts.values())

    print(f"\n  {label}")
    print(f"    {'Ox state':<10} {'Dimer-parent Fe':>20} {'Background (all Fe)':>25}")
    print("    " + "-" * 58)
    for ox in all_ox:
        dc = dimer_fe_ox_counts.get(ox, 0)
        bc = background_dist.get(ox, 0)
        d_pct = 100 * dc / n_dimer if n_dimer > 0 else 0
        b_pct = 100 * bc / n_fe_total if n_fe_total > 0 else 0
        print(f"    Fe^{ox:<5} {dc:>6}/{n_dimer:<4} ({d_pct:5.1f}%)   "
              f"{bc:>6}/{n_fe_total:<4} ({b_pct:5.1f}%)")


def analyze_oo_contacts_all_frames(temperature):
    """Scan all frames for O-O short contacts."""
    frames = FRAME_SETS[temperature]

    print(f"\n{'=' * 70}")
    print(f"  {temperature}: O-O contacts across all frames")
    print(f"{'=' * 70}")

    for frame_num, frame_idx in enumerate(frames):
        label = FRAME_LABELS[frame_num]
        structure, ox_states = load_frame(temperature, frame_idx)

        # Find O-O contacts at different thresholds
        dimers = find_oo_contacts(structure, OO_DIMER_CUTOFF)
        short_contacts = find_oo_contacts(structure, OO_SHORT_CUTOFF)

        print(f"\n  Frame {label} (step {frame_idx}):")

        # Ox state summary
        dist = fe_ox_distribution(ox_states)
        ox_str = ", ".join(f"Fe^{k}: {v}" for k, v in sorted(dist.items()))
        print(f"    Fe ox states: {ox_str}")

        o_ox = ox_states[36:180].astype(int)
        o_dist = dict(zip(*np.unique(o_ox, return_counts=True)))
        o_str = ", ".join(f"O^{int(k)}: {int(v)}" for k, v in sorted(o_dist.items()))
        print(f"    O ox states: {o_str}")

        if dimers:
            print(f"    O-O dimers (< {OO_DIMER_CUTOFF} A): {len(dimers)}")
            for o1, o2, d in dimers:
                species = classify_oo(d)
                ox_o1 = int(ox_states[o1])
                ox_o2 = int(ox_states[o2])
                print(f"      O{o1}--O{o2}: d={d:.4f} A [{species}], "
                      f"O ox states: ({ox_o1:+d}, {ox_o2:+d})")
        else:
            print(f"    O-O dimers (< {OO_DIMER_CUTOFF} A): none")

        # Short contacts beyond dimer cutoff
        precursors = [c for c in short_contacts if c[2] >= OO_DIMER_CUTOFF]
        if precursors:
            print(f"    Short O-O contacts ({OO_DIMER_CUTOFF}-{OO_SHORT_CUTOFF} A): "
                  f"{len(precursors)}")
            for o1, o2, d in precursors:
                ox_o1 = int(ox_states[o1])
                ox_o2 = int(ox_states[o2])
                print(f"      O{o1}--O{o2}: d={d:.4f} A, "
                      f"O ox states: ({ox_o1:+d}, {ox_o2:+d})")
        else:
            print(f"    Short O-O contacts ({OO_DIMER_CUTOFF}-{OO_SHORT_CUTOFF} A): none")


def trace_dimer_oxygens(temperature):
    """
    For O atoms that form dimers, trace their Fe coordination in frame i and frame iii,
    and compare dimer-parent Fe ox states against background.
    """
    frames = FRAME_SETS[temperature]

    print(f"\n{'=' * 70}")
    print(f"  {temperature}: Tracing dimer-forming O atoms")
    print(f"{'=' * 70}")

    # Load all frames
    structures = {}
    ox_arrays = {}
    for frame_num, frame_idx in enumerate(frames):
        label = FRAME_LABELS[frame_num]
        s, ox = load_frame(temperature, frame_idx)
        structures[label] = s
        ox_arrays[label] = ox

    # Identify the earliest frame containing dimers
    first_dimer_label = None
    for label in FRAME_LABELS:
        if find_oo_contacts(structures[label], OO_DIMER_CUTOFF):
            first_dimer_label = label
            break

    if first_dimer_label is None:
        print("  No O-O dimers found in any frame.")
        return

    first_dimer_num = FRAME_LABELS.index(first_dimer_label)
    first_dimer_idx = frames[first_dimer_num]
    dimers = find_oo_contacts(structures[first_dimer_label], OO_DIMER_CUTOFF)

    print(f"\n  First frame with dimers: {first_dimer_label} (step {first_dimer_idx})")
    for o1, o2, d in dimers:
        species = classify_oo(d)
        ox_o1 = int(ox_arrays[first_dimer_label][o1])
        ox_o2 = int(ox_arrays[first_dimer_label][o2])
        print(f"  Dimer: O{o1}--O{o2}, d={d:.4f} A [{species}], "
              f"O ox states: ({ox_o1:+d}, {ox_o2:+d})")

    # Collect all dimer O atoms
    dimer_o_atoms = set()
    for o1, o2, _ in dimers:
        dimer_o_atoms.add(o1)
        dimer_o_atoms.add(o2)

    # Trace Fe coordination of dimer O atoms through ALL frames
    for trace_num, trace_label in enumerate(FRAME_LABELS):
        trace_idx = frames[trace_num]
        struct_trace = structures[trace_label]
        ox_trace = ox_arrays[trace_label]

        is_dimer_frame = (trace_label == first_dimer_label)
        tag = " [dimer frame]" if is_dimer_frame else ""
        print(f"\n  --- Frame {trace_label} (step {trace_idx}){tag} ---")

        # Background Fe ox distribution
        bg_dist = fe_ox_distribution(ox_trace)
        bg_str = ", ".join(f"Fe^{k}: {v}" for k, v in sorted(bg_dist.items()))
        print(f"  Background Fe ox: {bg_str}")

        dimer_parent_fe_ox = defaultdict(int)

        for o_idx in sorted(dimer_o_atoms):
            fe_neighbors = find_fe_neighbors_of_o(struct_trace, o_idx)
            o_ox = int(ox_trace[o_idx])
            if fe_neighbors:
                fe_str = ", ".join(
                    f"Fe{fe_idx} (d={d:.3f} A, Fe^{int(ox_trace[fe_idx])})"
                    for fe_idx, d in fe_neighbors
                )
                for fe_idx, _ in fe_neighbors:
                    dimer_parent_fe_ox[int(ox_trace[fe_idx])] += 1
            else:
                fe_str = "no Fe within cutoff"
            print(f"    O{o_idx} (O^{o_ox:+d}) -> {fe_str}")

        # Background comparison
        print_comparison(
            f"Dimer-parent Fe vs background (frame {trace_label}):",
            dict(dimer_parent_fe_ox),
            bg_dist,
        )

        # Per-dimer shared-Fe analysis
        print(f"\n  Per-dimer shared-Fe (frame {trace_label}):")
        for o1, o2, d in dimers:
            species = classify_oo(d)
            fe_n1 = find_fe_neighbors_of_o(struct_trace, o1)
            fe_n2 = find_fe_neighbors_of_o(struct_trace, o2)
            fe_set1 = set(f for f, _ in fe_n1)
            fe_set2 = set(f for f, _ in fe_n2)
            shared = fe_set1 & fe_set2

            parts1 = ", ".join(f"Fe{f}(Fe^{int(ox_trace[f])})" for f, _ in fe_n1) or "none"
            parts2 = ", ".join(f"Fe{f}(Fe^{int(ox_trace[f])})" for f, _ in fe_n2) or "none"
            shared_str = (", ".join(f"Fe{f}" for f in shared) + " (shared)") if shared else "none"
            print(f"    O{o1}--O{o2}: O{o1}->[{parts1}], O{o2}->[{parts2}], shared=[{shared_str}]")


def check_500k_precursor():
    """
    Special check: in the 500K trajectory, are there any O-O distances in frame iii
    that are shorter than typical (even if above the dimer cutoff)?
    Look at the shortest O-O distances to check for precursor shortening.
    """
    print(f"\n{'=' * 70}")
    print("  500K: Precursor O-O shortening analysis")
    print(f"{'=' * 70}")

    frames = FRAME_SETS["500K"]

    # Check O-O distances in frame iii for the O atoms that eventually form dimers in iv
    struct_iv, ox_iv = load_frame("500K", frames[3])
    dimers_iv = find_oo_contacts(struct_iv, OO_DIMER_CUTOFF)

    if not dimers_iv:
        print("  No dimers in 500K frame iv to trace back.")
        return

    dimer_o_atoms = set()
    for o1, o2, _ in dimers_iv:
        dimer_o_atoms.add(o1)
        dimer_o_atoms.add(o2)

    print(f"  O atoms forming dimers in frame iv: {sorted(dimer_o_atoms)}")

    # For each earlier frame, check the O-O distance between these specific atoms
    for frame_num in range(4):
        label = FRAME_LABELS[frame_num]
        structure, ox_states = load_frame("500K", frames[frame_num])

        print(f"\n  Frame {label} (step {frames[frame_num]}):")
        for o1, o2, _ in dimers_iv:
            d = structure.get_distance(o1, o2)
            print(f"    O{o1}--O{o2}: d = {d:.4f} A")

    # Also check: in frame iii, what are the 10 shortest O-O distances overall?
    struct_iii, _ = load_frame("500K", frames[2])
    print(f"\n  10 shortest O-O distances in 500K frame iii (step {frames[2]}):")
    all_short = find_oo_contacts(struct_iii, 2.5)
    for o1, o2, d in all_short[:10]:
        marker = " <-- dimer pair" if (o1 in dimer_o_atoms and o2 in dimer_o_atoms) else ""
        print(f"    O{o1}--O{o2}: {d:.4f} A{marker}")

    # Same for 1000K
    print(f"\n{'=' * 70}")
    print("  1000K: Precursor O-O shortening analysis")
    print(f"{'=' * 70}")

    frames_1k = FRAME_SETS["1000K"]
    struct_iv_1k, ox_iv_1k = load_frame("1000K", frames_1k[3])
    dimers_iv_1k = find_oo_contacts(struct_iv_1k, OO_DIMER_CUTOFF)

    if not dimers_iv_1k:
        print("  No dimers in 1000K frame iv.")
        return

    dimer_o_1k = set()
    for o1, o2, _ in dimers_iv_1k:
        dimer_o_1k.add(o1)
        dimer_o_1k.add(o2)

    print(f"  O atoms forming dimers in frame iv: {sorted(dimer_o_1k)}")

    # Trace specific dimer pairs through all frames
    for frame_num in range(4):
        label = FRAME_LABELS[frame_num]
        structure, ox_states = load_frame("1000K", frames_1k[frame_num])

        print(f"\n  Frame {label} (step {frames_1k[frame_num]}):")
        for o1, o2, _ in dimers_iv_1k:
            d = structure.get_distance(o1, o2)
            print(f"    O{o1}--O{o2}: d = {d:.4f} A")

    # 10 shortest O-O in frame iii
    struct_iii_1k, _ = load_frame("1000K", frames_1k[2])
    print(f"\n  10 shortest O-O distances in 1000K frame iii (step {frames_1k[2]}):")
    all_short_1k = find_oo_contacts(struct_iii_1k, 2.5)
    for o1, o2, d in all_short_1k[:10]:
        marker = " <-- dimer pair" if (o1 in dimer_o_1k and o2 in dimer_o_1k) else ""
        print(f"    O{o1}--O{o2}: {d:.4f} A{marker}")


def main():
    print("=" * 70)
    print("R3.3: Tracing O2-Forming Oxygens Through the AIMD Trajectory")
    print("=" * 70)
    print(f"Fe-O coordination cutoff: {FE_O_CUTOFF} A")
    print(f"O-O dimer cutoff: {OO_DIMER_CUTOFF} A")
    print(f"O-O short contact cutoff: {OO_SHORT_CUTOFF} A")

    # Part 1: Scan all frames for O-O contacts
    print("\n\n" + "#" * 70)
    print("# PART 1: O-O contacts across all frames")
    print("#" * 70)
    for temp in ("500K", "1000K"):
        analyze_oo_contacts_all_frames(temp)

    # Part 2: Trace dimer-forming O atoms
    print("\n\n" + "#" * 70)
    print("# PART 2: Tracing dimer O atoms to their Fe coordination")
    print("#" * 70)
    for temp in ("500K", "1000K"):
        trace_dimer_oxygens(temp)

    # Part 3: Precursor shortening analysis
    print("\n\n" + "#" * 70)
    print("# PART 3: Precursor O-O shortening analysis")
    print("#" * 70)
    check_500k_precursor()


if __name__ == "__main__":
    main()
