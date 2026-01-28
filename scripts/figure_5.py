import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from figure_regeneration.colours import COLOURS
from figure_regeneration.formatting import set_formatting
from figure_regeneration.paths import EXTRACTED_DATA, FIGURES
from figure_regeneration.trajectory import moving_average


def main() -> None:
    prefix = EXTRACTED_DATA / "molecular_dynamics"

    temperatures = ("500K", "1000K")
    frame_sets = ((0, 6881, 16602, 39366), (0, 2581, 8315, 12276))
    window_size = 400

    relaxed_frame_energies_500 = np.load(prefix / "500K" / "relaxed_frame_energies.npy")
    pristine_energy = relaxed_frame_energies_500[0]

    fig_width = 3.15
    fig_height = fig_width * 1.75
    fig_size = (fig_width, fig_height)
    set_formatting(fig_size)

    fig, axes = plt.subplot_mosaic(
        [["500K"], ["1000K"], ["1000K"]], constrained_layout=True
    )

    for frame_indices, temperature in zip(frame_sets, temperatures):
        raw_energies = np.load(prefix / temperature / "energies.npy")
        relaxed_frame_energies = np.load(
            prefix / temperature / "relaxed_frame_energies.npy"
        )

        raw_energies -= pristine_energy
        relaxed_frame_energies -= pristine_energy

        average_energies = moving_average(raw_energies, window_size)

        time = tuple((idx * 2) / 1000 for idx in range(len(raw_energies)))
        average_time = time[: len(average_energies)]

        selected_frame_time = tuple(time[idx] for idx in frame_indices)

        pre_relaxed_frame_energies = []
        for t_i in selected_frame_time:
            for idx, t_j in enumerate(average_time):
                if t_i == t_j:
                    pre_relaxed_frame_energies.append(average_energies[idx])

        ax = axes[temperature]

        ax.plot(
            time,
            raw_energies,
            color=COLOURS["grape_soda"],
            alpha=0.2,
            zorder=1,
            linewidth=1,
        )
        ax.plot(average_time, average_energies, color=COLOURS["grape_soda"], zorder=2)
        ax.scatter(
            selected_frame_time,
            pre_relaxed_frame_energies,
            color=COLOURS["grape_soda"],
            zorder=3,
        )
        ax.scatter(selected_frame_time, relaxed_frame_energies, color="black", zorder=4)

        for idx, delta_t in enumerate(selected_frame_time):
            start = pre_relaxed_frame_energies[idx]
            end = relaxed_frame_energies[idx] + 2

            ax.annotate(
                "",
                xytext=(delta_t, start),
                xy=(delta_t, end),
                arrowprops=dict(arrowstyle="->"),
                zorder=2,
            )

        min_time, max_time = ax.get_xlim()
        for delta_t, energy in zip(selected_frame_time, relaxed_frame_energies):
            ax.plot(
                [min_time, delta_t],
                [energy, energy],
                linewidth=0.8,
                color=COLOURS["silver"],
            )

        ax.set_xlim(min_time, max_time)

    axes["500K"].set_xticks([0, 20, 40, 60, 80])
    axes["1000K"].set_xticks([0, 5, 10, 15, 20, 25])

    axes["500K"].set_ylim(-15, 20)
    axes["1000K"].set_ylim(-25, 30)

    for ax in axes.values():
        ax.set_xlabel(r"$\Delta t$ / ps")
        ax.set_ylabel(r"$\Delta E$ / eV")

        for direction in ("top", "right"):
            ax.spines[direction].set_visible(False)

    axes["500K"].annotate("i", (0.04, 0.85), xycoords="axes fraction")
    axes["500K"].annotate("ii", (0.16, 0.72), xycoords="axes fraction")
    axes["500K"].annotate("iii", (0.34, 0.71), xycoords="axes fraction")
    axes["500K"].annotate("iv", (0.76, 0.64), xycoords="axes fraction")

    axes["1000K"].annotate("i", (0.04, 0.99), xycoords="axes fraction")
    axes["1000K"].annotate("ii", (0.22, 0.94), xycoords="axes fraction")
    axes["1000K"].annotate("iii", (0.62, 0.83), xycoords="axes fraction")
    axes["1000K"].annotate("iv", (0.90, 0.79), xycoords="axes fraction")

    axes["500K"].annotate(r"(a) $500\,$K", (-0.2, 1.1), xycoords="axes fraction")
    axes["1000K"].annotate(r"(b) $1000\,$K", (-0.2, 1.1), xycoords="axes fraction")

    handles = (
        Line2D(
            [0],
            [0],
            marker="o",
            color="white",
            label="unrelaxed",
            markerfacecolor=COLOURS["grape_soda"],
            markersize=8,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="white",
            label="relaxed",
            markerfacecolor="black",
            markersize=8,
        ),
    )
    fig.legend(
        handles=handles,
        loc="upper right",
        ncols=1,
        bbox_to_anchor=(1.10, 1.05),
        handlelength=0.5,
    )

    fig.savefig(FIGURES / "figure_5.pdf")


if __name__ == "__main__":
    main()
