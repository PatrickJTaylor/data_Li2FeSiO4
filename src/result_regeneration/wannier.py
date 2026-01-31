# SPDX-FileCopyrightText: Copyright (C) 2026 Patrick J. Taylor
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np
from ase import Atoms
from ase.build.tools import sort
from ase.io import read
from numpy.typing import NDArray
from pymatgen.core import Molecule

from result_regeneration.dos import parse_icobi


def assign_wannier_centres(path: str) -> NDArray[np.float64]:
    """
    Assign Wannier centres to atoms to compute an effective "electron count".

    Parameters
    ----------
    path:
        Filepath to a directory containing the xyz files output by Wannier90, as well as
        a POSCAR from which to read the cell vectors.

    Returns
    -------
    electrons:
        The computed "electron count" for each atom in the structure.
    """
    base_atoms = read(f"{path}/POSCAR")
    n_atoms = len(base_atoms)
    electrons = np.zeros(n_atoms)

    for spin in (1, 2):
        molecule = Molecule.from_file(f"{path}/wannier90.{spin}_centres.xyz")

        symbols = []
        for site in molecule:
            symbol = "X" if site.species_string == "X0+" else site.species_string

            symbols.append(symbol)

        atoms = Atoms(
            cell=base_atoms.cell,
            symbols=symbols,
            positions=molecule.cart_coords,
            pbc=True,
        )
        atoms = sort(atoms)

        distance_matrix = atoms.get_all_distances(mic=True)
        wannier_indices = range(n_atoms, len(atoms))
        for i in wannier_indices:
            distances = distance_matrix[i, :n_atoms]
            min_idx = distances.argmin()

            electrons[min_idx] += 1

    return electrons


def correct_dimer_oxidation_states(path: str) -> NDArray[np.float64]:
    """
    Correct/adjust a set of Wannier-assigned oxidation states to account for the
    existence of dimeric species.

    This function effectively averages the oxidation states assigned to atoms in O-O
    dimers, to account for the possibility of uneven splitting of near-equidistant
    Wannier centres.

    Parameters
    ----------
    path:
        Filepath to a directory containing oxidation_states.npy and an associated
        ICOBILIST.lobster.

    Returns
    -------
    oxidation_states:
        The corrected/adjusted oxidation states.

    Notes
    -----
    As described in the supplementary information, O-O dimers are formally identified
    here by verifying that their associated ICOBIs are >= 0.95.
    """
    oxidation_states = np.load(f"{path}/oxidation_states.npy")
    icobis = parse_icobi(f"{path}/ICOBILIST.lobster")

    dimer_indices = []
    for symbols, icobi in icobis.items():
        symbol_i, symbol_j = symbols

        if "O" in symbol_i and "O" in symbol_j and icobi >= 0.95:
            i = int(symbol_i[1:]) - 1
            j = int(symbol_j[1:]) - 1

            dimer_indices.append((i, j))

    for dimer in dimer_indices:
        dimer_oxidation_states = tuple(
            oxidation_states[dimer_idx] for dimer_idx in dimer
        )
        average_oxidation_state = sum(dimer_oxidation_states) / 2

        # -0.75 is assigned to the constituent O atoms of peroxide dimers, so as to
        # differentiate them from isolated O^I- ions.
        new_oxidation_state = (
            -0.75 if average_oxidation_state == -1 else average_oxidation_state
        )
        for dimer_idx in dimer:
            oxidation_states[dimer_idx] = new_oxidation_state

    return oxidation_states
