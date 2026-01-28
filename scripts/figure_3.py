import matplotlib.pyplot as plt
import numpy as np

from result_regeneration.colours import COLOURS
from result_regeneration.dos import parse_icobi
from result_regeneration.formatting import set_formatting
from result_regeneration.paths import EXTRACTED_DATA, FIGURES


def main() -> None:
    prefix = EXTRACTED_DATA / "sequential_delithiation"

    compositions = ("Li2FeSiO4", "LiFeSiO4", "FeSiO4")
    xs = (2, 1, 0)

    fig_width = 3.15
    fig_height = fig_width * 1.6
    fig_size = (fig_width, fig_height)
    set_formatting(fig_size)

    fig, axes = plt.subplot_mosaic(
        [["magmoms"], ["icobis"]], constrained_layout=True, gridspec_kw={"hspace": 0.10}
    )

    average_magmoms, average_fe_o_icobis, average_si_o_icobis = [], [], []
    for x, composition in zip(xs, compositions):
        magmoms = np.load(prefix / composition / "magmoms.npy")
        icobis = parse_icobi(prefix / composition / "ICOBILIST.lobster")

        fe_o_icobis, si_o_icobis = [], []
        for symbols, icobi in icobis.items():
            symbol_i, symbol_j = symbols

            if "Fe" in symbol_i:
                fe_o_icobis.append(icobi)

            elif "Si" in symbol_j:
                si_o_icobis.append(icobi)

        average_magmoms.append(np.mean(magmoms))
        average_fe_o_icobis.append(np.mean(fe_o_icobis))
        average_si_o_icobis.append(np.mean(si_o_icobis))

        axes["magmoms"].scatter(
            np.repeat(x, len(magmoms)), magmoms, color=COLOURS["grape_soda"]
        )
        axes["icobis"].scatter(
            np.repeat(x, len(fe_o_icobis)), fe_o_icobis, color=COLOURS["grape_soda"]
        )
        axes["icobis"].scatter(
            np.repeat(x, len(si_o_icobis)), si_o_icobis, color=COLOURS["azure_blue"]
        )

    axes["magmoms"].plot(
        xs, average_magmoms, "o--", color=COLOURS["grape_soda"], markeredgecolor="black"
    )
    axes["icobis"].plot(
        xs,
        average_fe_o_icobis,
        "o--",
        color=COLOURS["grape_soda"],
        markeredgecolor="black",
    )
    axes["icobis"].plot(
        xs,
        average_si_o_icobis,
        "o--",
        color=COLOURS["azure_blue"],
        markeredgecolor="black",
    )

    for ax in axes.values():
        ax.invert_xaxis()
        ax.set_xticks([2, 1, 0])
        ax.set_xlabel(r"$x$ in Li$_{x}$FeSiO$_{4}$")

        for direction in ("top", "right"):
            ax.spines[direction].set_visible(False)

    axes["magmoms"].set_ylabel(r"$\mu_{\mathrm{B}}$")
    axes["icobis"].set_ylabel("ICOBI")

    axes["magmoms"].annotate(
        r"Fe$^{\mathrm{II}}$", (0.1, 0.35), xycoords="axes fraction"
    )
    axes["magmoms"].annotate(
        r"Fe$^{\mathrm{III}}$", (0.55, 0.9), xycoords="axes fraction"
    )
    axes["magmoms"].annotate(
        r"Fe$^{\mathrm{IV}}$", (1.0, 0.025), xycoords="axes fraction"
    )

    axes["icobis"].annotate("Fe-O", (0.15, 0.2), xycoords="axes fraction")
    axes["icobis"].annotate("Si-O", (0.15, 1.0), xycoords="axes fraction")

    axes["magmoms"].annotate(
        "(a) Fe magnetic moments", (-0.25, 1.10), xycoords="axes fraction"
    )
    axes["icobis"].annotate(
        "(b) (Fe/Si)-O ICOBI", (-0.25, 1.20), xycoords="axes fraction"
    )

    fig.savefig(FIGURES / "figure_3.pdf")


if __name__ == "__main__":
    main()
