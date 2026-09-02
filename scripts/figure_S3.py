# SPDX-FileCopyrightText: Copyright (C) 2026 Patrick J. Taylor
# SPDX-License-Identifier: GPL-3.0-or-later

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from scipy.stats import gaussian_kde

from result_regeneration.colours import COLOURS
from result_regeneration.dos import parse_icobi
from result_regeneration.errors import DependencyError
from result_regeneration.formatting import set_formatting
from result_regeneration.paths import EXTRACTED_DATA, FIGURES

TEMPERATURES = ("500K", "1000K")
FRAME_SETS = ((0, 6881, 16602, 39366), (0, 2581, 8315, 12276))
N_FE = 36

CMAP = {
    3: COLOURS["stormy_teal"],
    4: COLOURS["grape_soda"],
    5: COLOURS["amber_flame"],
}


def load_per_fe_data() -> tuple[dict[int, list[float]], dict[int, list[float]]]:
    """
    Collect, for every Fe in every relaxed AIMD frame, the absolute magnetic moment and
    the mean Fe-O ICOBI, grouped by Wannier-assigned oxidation state.
    """
    prefix = EXTRACTED_DATA / "molecular_dynamics"

    abs_fe_magmoms = {3: [], 4: [], 5: []}
    mean_fe_icobis = {3: [], 4: [], 5: []}

    for temperature, frames in zip(TEMPERATURES, FRAME_SETS):
        for frame_idx in frames:
            frame_dir = prefix / temperature / "selected_frames" / str(frame_idx)

            try:
                oxidation_states = np.load(frame_dir / "oxidation_states.npy")

            except OSError as base_error:
                raise DependencyError from base_error

            magmoms = np.load(frame_dir / "magmoms.npy")
            icobis = parse_icobi(frame_dir / "ICOBILIST.lobster")

            fe_icobis = [[] for _ in range(N_FE)]
            for (symbol_i, _), icobi in icobis.items():
                if "Fe" in symbol_i:
                    fe_icobis[int(symbol_i[2:]) - 1].append(icobi)

            for fe_idx, oxidation_state in enumerate(oxidation_states[:N_FE]):
                abs_fe_magmoms[oxidation_state].append(abs(magmoms[fe_idx]))
                mean_fe_icobis[oxidation_state].append(
                    sum(fe_icobis[fe_idx]) / len(fe_icobis[fe_idx])
                )

    return abs_fe_magmoms, mean_fe_icobis


def raincloud(
    ax: plt.Axes,
    data: dict[int, list[float]],
    x_range: np.ndarray,
    kde_threshold: float,
    kde_init_scale: float,
    kde_scaling: dict[int, float],
    rng: np.random.Generator,
    bandwidth: float = 0.04,
    baseline_offset: float = 0.35,
) -> None:
    """Draw one raincloud (KDE + box plot + jittered points) per oxidation state."""
    for baseline_idx, oxidation_state in enumerate(data):
        values = np.asarray(data[oxidation_state])
        colour = CMAP[oxidation_state]

        y_scatter = baseline_idx * 2 + rng.random(len(values)) * 0.15

        kde = gaussian_kde(values, bandwidth / np.std(values))
        baseline = baseline_idx * 2 + baseline_offset
        scaled_kde = (
            kde(x_range) * len(values) * kde_init_scale * kde_scaling[oxidation_state]
            + baseline
        )

        ax.plot(values, y_scatter, "o", color=colour, alpha=0.7, markeredgewidth=0)
        ax.axhline((baseline_idx * 2) + 0.075, color=COLOURS["silver"], alpha=0.2)
        ax.fill_between(
            x_range,
            scaled_kde,
            baseline,
            color=colour,
            where=kde(x_range) > kde_threshold,
            interpolate=False,
        )

        percentiles = np.percentile(values, [0, 25, 50, 75, 100])
        base_y = baseline_idx * 2 + 1.35

        ax.add_patch(
            Rectangle(
                xy=(percentiles[1], base_y - 0.15),
                width=percentiles[3] - percentiles[1],
                height=0.3,
                linewidth=1.0,
                edgecolor=colour,
                fill=False,
            )
        )
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

    ax.set_ylim(-0.5, 5.6)
    ax.set_yticks([])

    for direction in ("left", "top", "right"):
        ax.spines[direction].set_visible(False)


def main() -> None:
    abs_fe_magmoms, mean_fe_icobis = load_per_fe_data()

    rng = np.random.default_rng(seed=0)

    fig_width = 7.0
    fig_height = 3.15 * 1.25
    set_formatting((fig_width, fig_height))

    fig, (ax_a, ax_b) = plt.subplots(1, 2)

    # (a) absolute Fe magnetic moments
    raincloud(
        ax_a,
        abs_fe_magmoms,
        x_range=np.linspace(0.0, 4.5, 1000),
        kde_threshold=0.002,
        kde_init_scale=8e-04,
        kde_scaling={3: 2.60, 4: 0.65, 5: 16.50},
        rng=rng,
    )
    ax_a.set_xlabel(r"absolute Fe magnetic moment / $\mu_{\mathrm{B}}$")
    _, x_max = ax_a.get_xlim()
    ax_a.set_xlim(0.7, x_max)
    ax_a.annotate(r"Fe$^{\mathrm{III}}$", (0.77, 0.21), xycoords="axes fraction")
    ax_a.annotate(r"Fe$^{\mathrm{IV}}$", (0.56, 0.53), xycoords="axes fraction")
    ax_a.annotate(r"Fe$^{\mathrm{V}}$", (0.31, 0.87), xycoords="axes fraction")

    # (b) mean Fe-O ICOBI per Fe
    raincloud(
        ax_b,
        mean_fe_icobis,
        x_range=np.linspace(0.0, 1.0, 1000),
        kde_threshold=0.01,
        kde_init_scale=6e-04,
        kde_scaling={3: 3.20, 4: 0.85, 5: 9.50},
        rng=rng,
    )
    ax_b.set_xlabel("mean Fe-O ICOBI per Fe")
    ax_b.annotate(r"Fe$^{\mathrm{III}}$", (0.41, 0.21), xycoords="axes fraction")
    ax_b.annotate(r"Fe$^{\mathrm{IV}}$", (0.61, 0.54), xycoords="axes fraction")
    ax_b.annotate(r"Fe$^{\mathrm{V}}$", (0.66, 0.87), xycoords="axes fraction")

    ax_a.annotate("(a)", (-0.04, 1.02), xycoords="axes fraction")
    ax_b.annotate("(b)", (-0.04, 1.02), xycoords="axes fraction")

    fig.tight_layout()
    fig.savefig(FIGURES / "figure_S3.pdf")


if __name__ == "__main__":
    main()
