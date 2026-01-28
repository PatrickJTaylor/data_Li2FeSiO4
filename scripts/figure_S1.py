import json
import matplotlib.pyplot as plt

from result_regeneration.colours import COLOURS
from result_regeneration.formatting import set_formatting
from result_regeneration.paths import EXTRACTED_DATA, FIGURES


def main() -> None:
    prefix = EXTRACTED_DATA / "polymorph_stability"

    with open(prefix / "li2_data.json", "r") as stream:
        li2_data = json.load(stream)

    with open(prefix / "li0_data.json", "r") as stream:
        li0_data = json.load(stream)

    li2_volumes, li2_energies = [], []
    for polymorph, data in li2_data.items():
        volume = data["volume"] / data["n_atoms"]
        energy = (data["energy"] * 1000) / data["n_atoms"]

        li2_volumes.append(volume)
        li2_energies.append(energy)

    li0_volumes, li0_energies = [], []
    for polymorph, data in li0_data.items():
        volume = data["volume"] / data["n_atoms"]
        energy = (data["energy"] * 1000) / data["n_atoms"]

        li0_volumes.append(volume)
        li0_energies.append(energy)

    min_li2_energy = min(li2_energies)
    min_li0_energy = min(li0_energies)

    li2_energies = [energy - min_li2_energy for energy in li2_energies]
    li0_energies = [energy - min_li0_energy for energy in li0_energies]

    fig_width = 3.15
    fig_height = fig_width * 2.5
    fig_size = (fig_width, fig_height)
    set_formatting(fig_size)

    fig, axes = plt.subplots(2, 1)

    axes[0].scatter(li2_volumes, li2_energies, color=COLOURS["grape_soda"])
    axes[1].scatter(li0_volumes, li0_energies, color=COLOURS["grape_soda"])

    for ax in axes:
        ax.set_xlabel("volume / Å$^{3}$ per atom")
        ax.set_ylabel("relative energy / meV per atom")

        for direction in ("top", "right", "bottom", "left"):
            ax.spines[direction].set_visible(False)

        ax.set_facecolor(COLOURS["bright_snow"])

        ax.set_box_aspect(1)

    axes[0].annotate(r"Pmn2$_{1}$", (0.09, 0.025), xycoords="axes fraction")
    axes[0].annotate(
        "inverse\n " + r"Pmn2$_{1}$", (0.33, 0.75), xycoords="axes fraction"
    )
    axes[0].annotate(r"Pbn2$_{1}$", (0.52, 0.935), xycoords="axes fraction")
    axes[0].annotate(r"P2$_{1}$/n", (0.835, 0.86), xycoords="axes fraction")
    axes[0].annotate("Pmnb", (0.65, 0.425), xycoords="axes fraction")
    axes[0].annotate(r"P2$_{1}$", (0.535, 0.135), xycoords="axes fraction")

    axes[1].annotate(r"Pmn2$_{1}$", (0.245, 0.505), xycoords="axes fraction")
    axes[1].annotate(
        "inverse\n " + r"Pmn2$_{1}$", (0.75, 0.09), xycoords="axes fraction"
    )
    axes[1].annotate(r"Pbn2$_{1}$", (0.715, 0.935), xycoords="axes fraction")
    axes[1].annotate(r"P2$_{1}$/n", (0.75, 0.625), xycoords="axes fraction")
    axes[1].annotate("Pmnb", (0.31, 0.045), xycoords="axes fraction")
    axes[1].annotate(r"P2$_{1}$", (0.09, 0.065), xycoords="axes fraction")

    axes[0].annotate(r"(a) Li$_{2}$FeSiO$_{4}$", (-0.2, 1.08), xycoords="axes fraction")
    axes[1].annotate(r"(b) FeSiO$_{4}$", (-0.2, 1.08), xycoords="axes fraction")

    fig.subplots_adjust(hspace=0.5)
    fig.savefig(FIGURES / "figure_S1.pdf")


if __name__ == "__main__":
    main()
