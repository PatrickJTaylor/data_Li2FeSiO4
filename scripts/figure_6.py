# SPDX-FileCopyrightText: Copyright (C) 2026 Patrick J. Taylor
# SPDX-License-Identifier: GPL-3.0-or-later

import matplotlib.pyplot as plt
import numpy as np
from pymatgen.core import Structure
from vasppy.rdf import RadialDistributionFunction

from result_regeneration.colours import COLOURS
from result_regeneration.formatting import set_formatting
from result_regeneration.paths import EXTRACTED_DATA, FIGURES


def main() -> None:
    prefix = EXTRACTED_DATA / "molecular_dynamics"

    rdf_species = (("Fe", "Fe"), ("Fe", "Si"), ("O", "O"))
    sigma = 0.05

    plasma = plt.get_cmap("plasma")
    colour_indices = (0, 33, 66, 100, -56)
    plasma_colours = [plasma.colors[idx] for idx in colour_indices]

    frame_indices_500 = (0, 6881, 16602, 39366)
    frames = [
        Structure.from_file(prefix / "500K" / "selected_frames" / f"{idx}" / "POSCAR")
        for idx in frame_indices_500
    ]
    final_frame_1000 = Structure.from_file(
        prefix / "1000K" / "selected_frames" / "12276" / "POSCAR"
    )
    frames.append(final_frame_1000)

    inset_bounds = (0.00, 0.48, 0.20, 0.52)
    inset_500 = {
        "xlim": (0.9, 1.9),
        "xticks": ([1.4], ["1.40"]),
        "ylim": (0, 0.05),
        "rectangle_bounds": (0.9, -0.1, 1.0, 0.3),
    }
    inset_1000 = {
        "xlim": (0.74, 1.74),
        "xticks": ([1.24], ["1.24"]),
        "ylim": (0, 0.2),
        "rectangle_bounds": (0.74, -0.1, 1.0, 0.5),
    }
    inset_params_set = (inset_500, inset_1000)

    fe2o3_fe_fe_distance = 2.97

    golden_ratio = (1 + np.sqrt(5)) / 2

    fig_width = 6.698
    fig_height = fig_width / golden_ratio
    fig_size = (fig_width, fig_height)
    set_formatting(fig_size)

    fig, axes = plt.subplots(5, 3)

    for frame_idx in range(5):
        for species_idx, species in enumerate(rdf_species):
            species_i, species_j = species
            rdf = RadialDistributionFunction.from_species_strings(
                structures=[frames[frame_idx]], species_i=species_i, species_j=species_j
            )

            axes[frame_idx][species_idx].plot(
                rdf.r,
                rdf.smeared_rdf(sigma),
                linewidth=0.75,
                color=plasma_colours[frame_idx],
            )

            # Indicate Fe-Fe distance in Fe2O3 for Fe-Fe RDFs
            if species_idx == 0:
                axes[frame_idx][species_idx].axvline(
                    x=fe2o3_fe_fe_distance, linewidth=0.75, color=COLOURS["silver"]
                )

            # Generate insets to highlight O-O dimers
            if frame_idx in (3, 4) and species_idx == 2:
                inset_params = inset_params_set[frame_idx - 3]

                inset = axes[frame_idx][species_idx].inset_axes(inset_bounds)
                inset.plot(
                    rdf.r,
                    rdf.smeared_rdf(),
                    linewidth=0.75,
                    color=plasma_colours[frame_idx],
                )
                inset.set_xlim(inset_params["xlim"])
                inset.set_xticks(*inset_params["xticks"])
                inset.set_ylim(inset_params["ylim"])
                inset.set_yticks([])
                inset.set_yticks([], minor=True)

                for direction in ("left", "top", "right"):
                    inset.spines[direction].set_visible(False)

                inset_indicator = axes[frame_idx][species_idx].indicate_inset(
                    inset_params["rectangle_bounds"], inset, edgecolor=COLOURS["silver"]
                )
                for line in inset_indicator.connectors:
                    line.set(linewidth=0.75)

    frame_labels = ("i", "ii", "iii", "iv", "iv")

    for i in range(5):
        for j in range(3):
            ax = axes[i][j]

            ax.set_xlim(0, 10)
            ax.set_xticks([0.0, 2.5, 5.0, 7.5, 10.0])

            ax.set_yticks([])
            ax.set_yticks([], minor=True)

            for direction in ("left", "top", "right"):
                ax.spines[direction].set_visible(False)

            if j == 0:
                ax.annotate(frame_labels[i], (-0.1, 0.4), xycoords="axes fraction")

            if i < 4:
                ax.tick_params(axis="x", which="major", labelsize=0)

            else:
                ax.set_xlabel(r"$r$ / Å")

    axes[0][0].set_title("Fe-Fe")
    axes[0][1].set_title("Fe-Si")
    axes[0][2].set_title("O-O")

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.05)
    fig.savefig(FIGURES / "figure_6.pdf")


if __name__ == "__main__":
    main()
