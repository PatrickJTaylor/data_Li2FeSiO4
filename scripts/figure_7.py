import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from figure_regeneration.colours import COLOURS
from figure_regeneration.formatting import set_formatting
from figure_regeneration.paths import EXTRACTED_DATA, FIGURES
from figure_regeneration.wannier import correct_dimer_oxidation_states


def main() -> None:
    prefix = EXTRACTED_DATA / "molecular_dynamics"

    temperatures = ("500K", "1000K")
    frame_sets = ((0, 6881, 16602, 39366), (0, 2581, 8315, 12276))

    fe_states = (3, 4, 5)
    fe_cmap = {
        3: COLOURS["stormy_teal"],
        4: COLOURS["grape_soda"],
        5: COLOURS["amber_flame"],
    }

    # Note that -1 here refers only to individual O^I- ions (i.e., localised oxygen
    # holes). Oxygen atoms belonging to peroxide dimers, despite formally also
    # having oxidation states of -1, are labelled by -0.75 so as to differentiate
    # the two species in the final plots.
    o_states = (-2, -1, -0.75, -0.5, 0)
    o_cmap = {
        -2: COLOURS["coral_glow"],
        -1: COLOURS["jungle_green"],
        -0.75: COLOURS["powder_blue"],
        -0.5: COLOURS["vintage_lavender"],
        0: COLOURS["punch_red"],
    }

    n_fe = 36
    n_o = 144

    fig_width = 6.698
    fig_height = fig_width * 1.5
    fig_size = (fig_width, fig_height)
    set_formatting(fig_size)

    fig, axes = plt.subplots(2, 2)

    for idx, temperature in enumerate(temperatures):
        for col_idx, frame_idx in enumerate(frame_sets[idx]):
            oxidation_states = correct_dimer_oxidation_states(
                prefix / f"{temperature}" / "selected_frames" / f"{frame_idx}"
            )

            fe_oxidation_states = oxidation_states[:36]
            o_oxidation_states = oxidation_states[36:180]

            cumulative_fraction = 0
            for state_i in fe_states:
                n = 0
                for state_j in fe_oxidation_states:
                    if state_i == state_j:
                        n += 1

                fraction = n / n_fe

                axes[0][idx].bar(
                    col_idx,
                    fraction,
                    bottom=cumulative_fraction,
                    width=0.65,
                    color=fe_cmap[state_i],
                )
                cumulative_fraction += fraction

            cumulative_fraction = 0
            for state_i in o_states:
                n = 0
                for state_j in o_oxidation_states:
                    if state_i == state_j:
                        n += 1

                fraction = n / n_o

                axes[1][idx].bar(
                    col_idx,
                    fraction,
                    bottom=cumulative_fraction,
                    width=0.65,
                    color=o_cmap[state_i],
                )
                cumulative_fraction += fraction

    subplot_labels = (
        r"(a) $500\,$K",
        r"(b) $1000\,$K",
        r"(c) $500\,$K",
        r"(d) $1000\,$K",
    )

    for label, ax in zip(subplot_labels, axes.flatten()):
        ax.set_xticks([0, 1, 2, 3], ["i", "ii", "iii", "iv"])

        ax.set_ylim(0, 1)
        ax.set_ylabel("fraction")

        for direction in ("top", "right"):
            ax.spines[direction].set_visible(False)

        ax.annotate(label, (-0.275, 1.100), xycoords="axes fraction")

    fe_labels = {
        3: r"Fe$^{\mathrm{III}}$",
        4: r"Fe$^{\mathrm{IV}}$",
        5: r"Fe$^{\mathrm{V}}$",
    }
    fe_patches = []
    for oxidation_state in fe_states:
        patch = Patch(color=fe_cmap[oxidation_state], label=fe_labels[oxidation_state])

        fe_patches.append(patch)

    o_labels = {
        -2: r"O$^{\mathrm{II}-}$",
        -1: r"O$^{\mathrm{I}-}$",
        -0.75: r"O$^{\mathrm{II}-}_{2}$",
        -0.5: r"O$^{\mathrm{I}-}_{2}$",
        0: r"O$_{2}$",
    }
    o_patches = []
    for oxidation_state in o_states:
        patch = Patch(color=o_cmap[oxidation_state], label=o_labels[oxidation_state])

        o_patches.append(patch)

    axes[0][1].legend(
        handles=fe_patches,
        loc="center right",
        ncols=1,
        handlelength=0.5,
        handleheight=0.5,
        bbox_to_anchor=(1.4, 0.5),
    )
    axes[1][1].legend(
        handles=o_patches,
        loc="center right",
        ncols=1,
        handlelength=0.5,
        handleheight=0.5,
        bbox_to_anchor=(1.4, 0.5),
    )

    fig.subplots_adjust(hspace=0.3, wspace=0.4)
    fig.savefig(FIGURES / "figure_7.pdf")


if __name__ == "__main__":
    main()
