# SPDX-FileCopyrightText: Copyright (C) 2026 Patrick J. Taylor
# SPDX-License-Identifier: GPL-3.0-or-later

import json
import matplotlib.pyplot as plt

from result_regeneration.colours import COLOURS
from result_regeneration.formatting import set_formatting
from result_regeneration.paths import EXTRACTED_DATA, FIGURES


def main() -> None:
    prefix = EXTRACTED_DATA / "hull"

    with open(prefix / "hull_data.json", "r") as stream:
        hull_data = json.load(stream)

    li2_stable = hull_data[2]["energy"]
    li1_stable = hull_data[1]["energy"]
    li0_stable = hull_data[0]["energy"]

    li2_per_atom = li2_stable / hull_data[2]["n_atoms"]
    li0_per_atom = li0_stable / hull_data[0]["n_atoms"]

    compositions, mixing_energies = [], []
    for config_data in hull_data:
        composition = (config_data["n_li"] / config_data["n_o"]) * 2
        energy = config_data["energy"] / config_data["n_atoms"]

        mixing_energy = energy - (
            composition * li2_per_atom + (1 - composition) * li0_per_atom
        )

        compositions.append(composition)
        mixing_energies.append(mixing_energy)

    li_compositions = tuple(composition * 2 for composition in compositions)
    li1_stable_mixing_energy = mixing_energies[1]

    hull_compositions = (0, 1, 2)
    hull_energies = (0, li1_stable_mixing_energy, 0)

    li2_per_fu = li2_stable / 2
    li1_per_fu = li1_stable / 2
    li0_per_fu = li0_stable / 2

    with open(prefix / "li_energy.json", "r") as stream:
        li_energy = json.load(stream)

    v1 = -1 * (li2_per_fu - li1_per_fu - li_energy)
    v2 = -1 * (li1_per_fu - li0_per_fu - li_energy)

    x = [2, 1, 1, 0]
    voltage = [v1, v1, v2, v2]

    fig_width = 3.15
    fig_height = fig_width * 1.75
    fig_size = (fig_width, fig_height)
    set_formatting(fig_size)

    fig, axes = plt.subplots(2, 1)

    axes[0].plot(
        li_compositions,
        mixing_energies,
        "o",
        color=COLOURS["coral_glow"],
        zorder=2,
        alpha=0.7,
        markeredgecolor=None,
        markeredgewidth=0,
    )
    axes[0].plot(
        hull_compositions, hull_energies, color=COLOURS["grape_soda"], zorder=1
    )

    axes[1].plot(x, voltage, color=COLOURS["grape_soda"])

    axes[0].set_ylim(-0.135, 0.01)
    axes[0].set_yticks([-0.10, -0.05, 0.00])
    axes[0].set_ylabel("mixing energy /\n eV per atom")

    axes[1].set_xlabel(r"$x$ in Li$_{x}$FeSiO$_{4}$")

    axes[1].set_ylim(2.2, 5.1)
    axes[1].set_yticks([3.0, 3.5, 4.0, 4.5, 5.0])
    axes[1].set_ylabel("voltage / V")

    for ax in axes:
        ax.invert_xaxis()

        for direction in ("top", "right"):
            ax.spines[direction].set_visible(False)

    axes[0].annotate("(a)", (-0.35, 1.15), xycoords="axes fraction")
    axes[1].annotate("(b)", (-0.35, 1.15), xycoords="axes fraction")

    fig.subplots_adjust(hspace=0.6)
    fig.savefig(FIGURES / "figure_S2.pdf")


if __name__ == "__main__":
    main()
