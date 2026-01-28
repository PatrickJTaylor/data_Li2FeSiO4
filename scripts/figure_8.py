import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from scipy.stats import gaussian_kde

from figure_regeneration.colours import COLOURS
from figure_regeneration.dos import parse_icobi
from figure_regeneration.formatting import set_formatting
from figure_regeneration.paths import EXTRACTED_DATA, FIGURES


def main() -> None:
    prefix = EXTRACTED_DATA / "molecular_dynamics"

    temperatures = ("500K", "1000K")
    frame_sets = ((0, 6881, 16602, 39366), (0, 2581, 8315, 12276))

    n_fe = 36

    mean_fe_icobis = {3: [], 4: [], 5: []}
    for idx, temperature in enumerate(temperatures):
        for frame_idx in frame_sets[idx]:
            oxidation_states = np.load(
                prefix
                / f"{temperature}"
                / "selected_frames"
                / f"{frame_idx}"
                / "oxidation_states.npy"
            )
            icobis = parse_icobi(
                prefix
                / f"{temperature}"
                / "selected_frames"
                / f"{frame_idx}"
                / "ICOBILIST.lobster"
            )

            fe_icobis = [[] for _ in range(n_fe)]
            for symbols, icobi in icobis.items():
                symbol_i, symbol_j = symbols

                if "Fe" in symbol_i:
                    fe_idx = int(symbol_i[2:]) - 1

                    fe_icobis[fe_idx].append(icobi)

            fe_oxidation_states = oxidation_states[:n_fe]
            for oxidation_state, icobis in zip(fe_oxidation_states, fe_icobis):
                mean_icobi = sum(icobis) / len(icobis)

                mean_fe_icobis[oxidation_state].append(mean_icobi)

    x_range = np.linspace(0, 1, 1000)
    baseline_offset = 0.35
    bandwidth = 0.04
    kde_threshold = 0.01
    kde_init_scale = 6e-04

    kde_scaling = {3: 3.20, 4: 0.85, 5: 9.50}

    cmap = {
        3: COLOURS["stormy_teal"],
        4: COLOURS["grape_soda"],
        5: COLOURS["amber_flame"],
    }

    fig_width = 3.15
    fig_height = fig_width * 1.25
    fig_size = (fig_width, fig_height)
    set_formatting(fig_size)

    fig, ax = plt.subplots()

    for baseline_idx, oxidation_state in enumerate(mean_fe_icobis):
        mean_icobis = mean_fe_icobis[oxidation_state]
        y_scatter = (
            np.full_like(mean_icobis, baseline_idx) * 2
            + np.random.random(size=len(mean_icobis)) * 0.15
        )

        kde = gaussian_kde(mean_icobis, bandwidth / np.std(mean_icobis))
        baseline = baseline_idx * 2 + baseline_offset
        scaled_kde = (
            kde(x_range)
            * len(mean_icobis)
            * kde_init_scale
            * kde_scaling[oxidation_state]
            + baseline
        )

        colour = cmap[oxidation_state]

        ax.plot(
            mean_icobis,
            y_scatter,
            "o",
            color=colour,
            alpha=0.7,
            markeredgecolor=None,
            markeredgewidth=0,
        )
        ax.axhline((baseline_idx * 2) + 0.075, color=COLOURS["silver"], alpha=0.2)
        ax.fill_between(
            x_range,
            scaled_kde,
            baseline,
            color=colour,
            where=kde(x_range) > kde_threshold,
            interpolate=False,
        )

        percentiles = np.percentile(mean_icobis, [0, 25, 50, 75, 100])
        base_y = baseline_idx * 2 + 1.35

        box = Rectangle(
            xy=(percentiles[1], base_y - 0.15),
            width=percentiles[3] - percentiles[1],
            height=0.3,
            linewidth=1.0,
            edgecolor=colour,
            fill=False,
        )
        ax.add_patch(box)

        ax.vlines(
            (percentiles[0], percentiles[-1]),
            base_y - 0.075,
            base_y + 0.075,
            linewidth=1.0,
            color=colour,
        )
        ax.vlines(
            percentiles[2], base_y - 0.15, base_y + 0.15, linewidth=1.0, color=colour
        )

        ax.hlines(base_y, percentiles[0], percentiles[1], linewidth=1.0, color=colour)
        ax.hlines(base_y, percentiles[3], percentiles[-1], linewidth=1.0, color=colour)

    ax.set_xlabel("mean Fe-O ICOBI per Fe")

    ax.set_ylim(-0.5, 5.6)
    ax.set_yticks([])

    for direction in ("left", "top", "right"):
        ax.spines[direction].set_visible(False)

    ax.annotate(r"Fe$^{\mathrm{III}}$", (0.41, 0.21), xycoords="axes fraction")
    ax.annotate(r"Fe$^{\mathrm{IV}}$", (0.61, 0.54), xycoords="axes fraction")
    ax.annotate(r"Fe$^{\mathrm{V}}$", (0.66, 0.87), xycoords="axes fraction")

    fig.tight_layout()
    fig.savefig(FIGURES / "figure_8.pdf")


if __name__ == "__main__":
    main()
