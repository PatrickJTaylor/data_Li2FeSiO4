# SPDX-FileCopyrightText: Copyright (C) 2026 Patrick J. Taylor
# SPDX-License-Identifier: GPL-3.0-or-later

import json

from result_regeneration.paths import EXTRACTED_DATA


def main() -> None:
    prefix = EXTRACTED_DATA / "thermodynamic_stability"

    with open(prefix / "stability_data.json", "r") as stream:
        stability_data = json.load(stream)

    energies_per_fu = {}
    for formula, data in stability_data.items():
        energy_per_fu = data["energy"] / data["n_fu"]

        energies_per_fu[formula] = energy_per_fu

    reactants = energies_per_fu["FeSiO4"]
    products = (
        0.5 * energies_per_fu["Fe2O3"]
        + energies_per_fu["SiO2"]
        + 0.25 * energies_per_fu["O2"]
    )
    energy_above_hull_per_fu = products - reactants

    n_atoms_per_fu_fesio4 = 6
    energy_above_hull_per_atom = energy_above_hull_per_fu / n_atoms_per_fu_fesio4

    print(f"FeSiO4 is {energy_above_hull_per_atom} eV / atom above the hull")


if __name__ == "__main__":
    main()
