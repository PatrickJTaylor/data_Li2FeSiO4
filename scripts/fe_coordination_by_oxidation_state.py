# SPDX-FileCopyrightText: Copyright (C) 2026 Benjamin J. Morgan
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Fe coordination environments by Wannier-assigned oxidation state in the relaxed AIMD frames.

For each relaxed frame, reports:
  - Fe-O coordination numbers, with cutoff sensitivity (2.2-2.5 A; 2.4 A is used in the manuscript).
  - Effective coordination numbers (ECoN, Hoppe's method) as a continuous metric.
  - Cross-tabulations of oxidation state against coordination number (counts and percentages).
  - Fe-O distance distributions by oxidation state.

Supports the Fe coordination statistics quoted in Sec. III C and III D of the manuscript.

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
N_O = 144
FE_INDICES = range(N_FE)
O_INDICES = range(36, 180)

CUTOFFS = [2.2, 2.3, 2.4, 2.5]
DEFAULT_CUTOFF = 2.4
ECON_SEARCH_RADIUS = 3.5  # generous radius for ECoN calculation


def classify_cn(cn):
    """Classify integer coordination number."""
    labels = {
        2: "CN=2",
        3: "CN=3",
        4: "CN=4 (tet)",
        5: "CN=5",
        6: "CN=6 (oct)",
        7: "CN=7",
        8: "CN=8",
    }
    return labels.get(cn, f"CN={cn}")


def get_fe_o_data(structure, search_radius=ECON_SEARCH_RADIUS):
    """
    For each Fe atom, get all O neighbours within search_radius.

    Returns list of dicts with 'distances' and 'o_indices' (sorted by distance).
    """
    results = []
    for fe_idx in FE_INDICES:
        neighbors = structure.get_neighbors(structure[fe_idx], r=search_radius)
        o_neighbors = [
            (n.index, n.nn_distance)
            for n in neighbors
            if n.species_string == "O"
        ]
        o_neighbors.sort(key=lambda x: x[1])
        results.append({
            "distances": np.array([d for _, d in o_neighbors]),
            "o_indices": [idx for idx, _ in o_neighbors],
        })
    return results


def coordination_number(fe_data, cutoff):
    """Hard-cutoff coordination number."""
    return int(np.sum(fe_data["distances"] <= cutoff))


def effective_coordination_number(fe_data, cutoff_for_pool=ECON_SEARCH_RADIUS):
    """
    Hoppe's effective coordination number (ECoN).

    ECoN = sum_i exp(1 - (d_i / d_av)^6)
    where d_av is the weighted average bond length, iterated to self-consistency.

    Uses all O within cutoff_for_pool.
    """
    dists = fe_data["distances"][fe_data["distances"] <= cutoff_for_pool]
    if len(dists) == 0:
        return 0.0

    d_min = dists[0]
    # Initial d_av = d_min
    d_av = d_min
    for _ in range(20):  # iterate to convergence
        weights = np.exp(1.0 - (dists / d_av) ** 6)
        d_av_new = np.sum(weights * dists) / np.sum(weights)
        if abs(d_av_new - d_av) < 1e-6:
            break
        d_av = d_av_new
    econ = np.sum(np.exp(1.0 - (dists / d_av) ** 6))
    return econ


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


def print_crosstab(title, row_label, col_label, table, row_keys, col_keys):
    """Print a cross-tabulation with counts and row/column totals."""
    col_w = max(10, max((len(ck) for ck in col_keys), default=10) + 2)
    row_w = max(14, max((len(rk) for rk in row_keys), default=14) + 2)

    print(f"\n  {title}")
    header = f"    {row_label:<{row_w}}"
    for ck in col_keys:
        header += f"{ck:>{col_w}}"
    header += f"{'Total':>{col_w}}"
    print(header)
    print("    " + "-" * (row_w + col_w * (len(col_keys) + 1)))

    col_totals = {ck: 0 for ck in col_keys}
    grand = 0
    for rk in row_keys:
        row_total = sum(table[rk].values())
        row = f"    {rk:<{row_w}}"
        for ck in col_keys:
            count = table[rk].get(ck, 0)
            row += f"{count:>{col_w}}"
            col_totals[ck] += count
        row += f"{row_total:>{col_w}}"
        grand += row_total
        print(row)

    row = f"    {'Total':<{row_w}}"
    for ck in col_keys:
        row += f"{col_totals[ck]:>{col_w}}"
    row += f"{grand:>{col_w}}"
    print(row)


def print_distance_stats(fe_data_list, fe_ox, cutoff):
    """Print Fe-O distance statistics grouped by oxidation state."""
    # Collect distances within cutoff for each ox state
    by_ox = defaultdict(list)
    for fe_idx in FE_INDICES:
        ox = int(fe_ox[fe_idx])
        dists = fe_data_list[fe_idx]["distances"]
        dists_in = dists[dists <= cutoff]
        by_ox[ox].extend(dists_in.tolist())

    print(f"\n  Fe-O distance statistics (within {cutoff} A cutoff):")
    print(f"    {'Ox state':<10} {'N_bonds':>8} {'mean':>7} {'std':>7} "
          f"{'min':>7} {'Q25':>7} {'median':>7} {'Q75':>7} {'max':>7}")
    print("    " + "-" * 78)
    for ox in sorted(by_ox.keys()):
        d = np.array(by_ox[ox])
        if len(d) == 0:
            continue
        print(f"    Fe^{ox:<5} {len(d):>8} {d.mean():>7.3f} {d.std():>7.3f} "
              f"{d.min():>7.3f} {np.percentile(d, 25):>7.3f} "
              f"{np.median(d):>7.3f} {np.percentile(d, 75):>7.3f} {d.max():>7.3f}")

    # Also report the full distance distribution (including beyond cutoff) up to 3.0 A
    print(f"\n  All Fe-O distances up to 3.0 A (for cutoff context):")
    by_ox_full = defaultdict(list)
    for fe_idx in FE_INDICES:
        ox = int(fe_ox[fe_idx])
        dists = fe_data_list[fe_idx]["distances"]
        dists_in = dists[dists <= 3.0]
        by_ox_full[ox].extend(dists_in.tolist())

    for ox in sorted(by_ox_full.keys()):
        d = np.array(by_ox_full[ox])
        # Histogram bins
        bins = [1.6, 1.8, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.8, 3.0]
        hist, _ = np.histogram(d, bins=bins)
        print(f"    Fe^{ox} (n_Fe={sum(1 for i in FE_INDICES if int(fe_ox[i]) == ox)}, "
              f"n_bonds={len(d)}):")
        for k in range(len(bins) - 1):
            bar = "#" * hist[k]
            print(f"      {bins[k]:.1f}-{bins[k+1]:.1f} A: {hist[k]:>4}  {bar}")


def main():
    print("=" * 90)
    print("R3.2: Fe Coordination Environments by Wannier Oxidation State")
    print("=" * 90)

    # =========================================================================
    # PART 1: Cutoff sensitivity analysis (compact summary)
    # =========================================================================
    print("\n" + "=" * 90)
    print("PART 1: Cutoff sensitivity analysis")
    print("=" * 90)
    print("For each frame, how does the ox-state vs CN cross-tab change with cutoff?")

    for temperature in TEMPERATURES:
        frames = FRAME_SETS[temperature]
        print(f"\n{'=' * 70}")
        print(f"  {temperature}")
        print(f"{'=' * 70}")

        for frame_num, frame_idx in enumerate(frames):
            label = FRAME_LABELS[frame_num]
            structure, ox_states = load_frame(temperature, frame_idx)
            fe_ox = ox_states[:N_FE]
            fe_data = get_fe_o_data(structure)

            print(f"\n  --- Frame {label} (step {frame_idx}) ---")
            # Ox state summary
            unique, counts = np.unique(fe_ox.astype(int), return_counts=True)
            ox_summary = ", ".join(f"Fe^{u}={c}" for u, c in zip(unique, counts))
            print(f"  Fe oxidation states: {ox_summary}")

            for cutoff in CUTOFFS:
                # Build ox -> CN table
                table = defaultdict(lambda: defaultdict(int))
                for fe_idx in FE_INDICES:
                    ox = f"Fe^{int(fe_ox[fe_idx])}"
                    cn = coordination_number(fe_data[fe_idx], cutoff)
                    cn_label = classify_cn(cn)
                    table[ox][cn_label] += 1

                ox_keys = sorted(table.keys())
                cn_keys = sorted(set(c for row in table.values() for c in row))

                print(f"\n  Cutoff = {cutoff} A:")
                header = f"    {'':>10}"
                for ck in cn_keys:
                    header += f"  {ck:>12}"
                print(header)
                for ok in ox_keys:
                    row = f"    {ok:>10}"
                    for ck in cn_keys:
                        row += f"  {table[ok].get(ck, 0):>12}"
                    print(row)

    # =========================================================================
    # PART 2: Detailed frame-by-frame analysis (default cutoff + ECoN)
    # =========================================================================
    print("\n\n" + "=" * 90)
    print(f"PART 2: Detailed frame-by-frame analysis (cutoff = {DEFAULT_CUTOFF} A + ECoN)")
    print("=" * 90)

    for temperature in TEMPERATURES:
        frames = FRAME_SETS[temperature]
        print(f"\n{'#' * 80}")
        print(f"# {temperature}")
        print(f"{'#' * 80}")

        for frame_num, frame_idx in enumerate(frames):
            label = FRAME_LABELS[frame_num]
            structure, ox_states = load_frame(temperature, frame_idx)
            fe_ox = ox_states[:N_FE]
            fe_data = get_fe_o_data(structure)

            print(f"\n{'=' * 70}")
            print(f"  Frame {label} (step {frame_idx}), {temperature}")
            print(f"{'=' * 70}")

            # Ox state summary
            unique, counts = np.unique(fe_ox.astype(int), return_counts=True)
            ox_summary = ", ".join(f"Fe^{u}: {c}" for u, c in zip(unique, counts))
            print(f"  Fe oxidation state counts: {ox_summary}")

            # Compute CN and ECoN for each Fe
            cn_list = []
            econ_list = []
            for fe_idx in FE_INDICES:
                cn = coordination_number(fe_data[fe_idx], DEFAULT_CUTOFF)
                econ = effective_coordination_number(fe_data[fe_idx])
                cn_list.append(cn)
                econ_list.append(econ)

            # (a) Oxidation state -> coordination environment
            ox_to_cn = defaultdict(lambda: defaultdict(int))
            for fe_idx in FE_INDICES:
                ox_key = f"Fe^{int(fe_ox[fe_idx])}"
                cn_key = classify_cn(cn_list[fe_idx])
                ox_to_cn[ox_key][cn_key] += 1

            ox_keys = sorted(ox_to_cn.keys())
            cn_keys = sorted(set(c for row in ox_to_cn.values() for c in row))
            print_crosstab(
                "(a) Given oxidation state, what coordination?",
                "Ox state", "CN", ox_to_cn, ox_keys, cn_keys
            )

            # (b) Coordination environment -> oxidation state
            cn_to_ox = defaultdict(lambda: defaultdict(int))
            for fe_idx in FE_INDICES:
                ox_key = f"Fe^{int(fe_ox[fe_idx])}"
                cn_key = classify_cn(cn_list[fe_idx])
                cn_to_ox[cn_key][ox_key] += 1

            cn_keys_b = sorted(cn_to_ox.keys())
            ox_keys_b = sorted(set(o for row in cn_to_ox.values() for o in row))
            print_crosstab(
                "(b) Given coordination, what oxidation state?",
                "CN", "Ox state", cn_to_ox, cn_keys_b, ox_keys_b
            )

            # ECoN statistics by oxidation state
            print(f"\n  Effective coordination number (ECoN) by oxidation state:")
            print(f"    {'Ox state':<10} {'N':>4} {'mean':>7} {'std':>7} "
                  f"{'min':>7} {'max':>7}  individual values")
            print("    " + "-" * 80)
            for ox_val in sorted(unique):
                ox_val = int(ox_val)
                mask = fe_ox.astype(int) == ox_val
                econs = np.array(econ_list)[mask]
                vals = ", ".join(f"{v:.2f}" for v in sorted(econs))
                if len(econs) > 1:
                    print(f"    Fe^{ox_val:<5} {len(econs):>4} {econs.mean():>7.2f} "
                          f"{econs.std():>7.2f} {econs.min():>7.2f} {econs.max():>7.2f}  "
                          f"[{vals}]")
                else:
                    print(f"    Fe^{ox_val:<5} {len(econs):>4} {econs.mean():>7.2f} "
                          f"{'--':>7} {econs.min():>7.2f} {econs.max():>7.2f}  "
                          f"[{vals}]")

            # Fe-O distance distributions by oxidation state
            print_distance_stats(fe_data, fe_ox, DEFAULT_CUTOFF)

            # Detailed listing of non-tetrahedral Fe atoms
            non_tet = [(i, cn_list[i], econ_list[i]) for i in FE_INDICES if cn_list[i] != 4]
            if non_tet:
                print(f"\n  Non-tetrahedral Fe atoms ({len(non_tet)} of {N_FE}):")
                for fe_idx, cn, econ in non_tet:
                    dists = fe_data[fe_idx]["distances"]
                    dists_in = dists[dists <= DEFAULT_CUTOFF]
                    d_str = ", ".join(f"{d:.3f}" for d in dists_in)
                    # Also show next distance beyond cutoff
                    dists_out = dists[dists > DEFAULT_CUTOFF]
                    next_d = f", next: {dists_out[0]:.3f}" if len(dists_out) > 0 else ""
                    print(f"    Fe{fe_idx}: Fe^{int(fe_ox[fe_idx])}, CN={cn}, "
                          f"ECoN={econ:.2f}, d=[{d_str}]{next_d}")
            else:
                print(f"\n  All {N_FE} Fe atoms are tetrahedral (CN=4).")


if __name__ == "__main__":
    main()
