import numpy as np
from pymatgen.core import Structure

from figure_regeneration.paths import EXTRACTED_DATA
from figure_regeneration.wannier import assign_wannier_centres


def main() -> None:
    prefix = EXTRACTED_DATA / "sequential_delithiation"

    compositions = ("Li2FeSiO4", "LiFeSiO4", "FeSiO4")
    valence_counts = {"Li": 3, "Fe": 8, "Si": 4, "O": 6}

    for composition in compositions:
        structure = Structure.from_file(prefix / composition / "POSCAR")
        electrons = assign_wannier_centres(prefix / composition)

        print(composition)
        for idx, (site, electron_count) in enumerate(zip(structure, electrons)):
            oxidation_state = valence_counts[site.species_string] - electron_count

            print(f"{site.species_string}{idx}: {oxidation_state}")


if __name__ == "__main__":
    main()
