# SPDX-FileCopyrightText: Copyright (C) 2026 Patrick J. Taylor
# SPDX-License-Identifier: GPL-3.0-or-later

import matplotlib.pyplot as plt

from result_regeneration.colours import COLOURS
from result_regeneration.dos import parse_pdos
from result_regeneration.formatting import set_formatting
from result_regeneration.paths import EXTRACTED_DATA, FIGURES


def main() -> None:
    prefix = EXTRACTED_DATA / "sequential_delithiation"

    compositions = ("Li2FeSiO4", "LiFeSiO4", "FeSiO4")
    x_shifts = (-0.10, -0.05, 0.68)
    sigma = 4
    y_shift = 9

    fig_width = 3.15
    fig_height = fig_width * 2.1
    fig_size = (fig_width, fig_height)
    set_formatting(fig_size)

    fig, axes = plt.subplot_mosaic(
        [["Li2FeSiO4"], ["LiFeSiO4"], ["FeSiO4"]],
        constrained_layout=True,
        gridspec_kw={"hspace": 0.15},
    )

    for composition, x_shift in zip(compositions, x_shifts):
        energies, pdos = parse_pdos(prefix / composition, sigma, x_shift, y_shift)
        fe_up, fe_down, o_up, o_down = pdos

        ax = axes[composition]

        ax.fill_between(
            energies,
            fe_up,
            y2=y_shift,
            lw=0,
            facecolor=COLOURS["grape_soda"],
            alpha=0.6,
        )
        ax.fill_between(
            energies,
            fe_down,
            y2=y_shift,
            lw=0,
            facecolor=COLOURS["grape_soda"],
            alpha=0.6,
        )
        ax.plot(energies, fe_up, lw=0.5, color=COLOURS["grape_soda"])
        ax.plot(energies, fe_down, lw=0.5, color=COLOURS["grape_soda"])

        ax.fill_between(
            energies,
            o_up,
            y2=-y_shift,
            lw=0,
            facecolor=COLOURS["coral_glow"],
            alpha=0.6,
        )
        ax.fill_between(
            energies,
            o_down,
            y2=-y_shift,
            lw=0,
            facecolor=COLOURS["coral_glow"],
            alpha=0.6,
        )
        ax.plot(energies, o_up, lw=0.5, color=COLOURS["coral_glow"])
        ax.plot(energies, o_down, lw=0.5, color=COLOURS["coral_glow"])

        ax.axvline(x=0, lw=0.5, color=COLOURS["silver"])

    for ax in axes.values():
        ax.annotate("Fe", (0.93, 0.75), xycoords="axes fraction")
        ax.annotate("O", (0.945, 0.30), xycoords="axes fraction")

        ax.set_xlim(-10, 5)
        ax.set_xticks([-10, -5, 0, 5])

        ax.set_ylim(-20, 20)
        ax.set_yticks([])
        ax.set_ylabel("DOS")

        for direction in ("left", "top", "right"):
            ax.spines[direction].set_visible(False)

    axes["FeSiO4"].set_xlabel(r"$E - E_{\mathrm{F}}$ / eV")

    axes["Li2FeSiO4"].annotate(
        r"(a) Li$_{2}$FeSiO$_{4}$", (-0.1, 1.05), xycoords="axes fraction"
    )
    axes["LiFeSiO4"].annotate(
        r"(b) LiFeSiO$_{4}$", (-0.1, 1.05), xycoords="axes fraction"
    )
    axes["FeSiO4"].annotate(r"(c) FeSiO$_{4}$", (-0.1, 1.05), xycoords="axes fraction")

    fig.savefig(FIGURES / "figure_2.pdf")


if __name__ == "__main__":
    main()
