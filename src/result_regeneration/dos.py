# SPDX-FileCopyrightText: Copyright (C) 2026 Patrick J. Taylor
# SPDX-License-Identifier: GPL-3.0-or-later

import numpy as np
from numpy.typing import NDArray
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.lobster import Doscar
from scipy.ndimage import gaussian_filter1d


def parse_pdos(
    path: str, sigma: float | None = None, x_shift: float = 0.0, y_shift: float = 0.0
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Parse the projected density of states (PDOS) for configurations of Li_xFeSiO4.

    This function extracts the Fe-3d and O-2p PDOS for both spin channels alongside
    the energies at which they have been evaluated.

    Parameters
    ----------
    path:
        Filepath to a directory containing a DOSCAR.lobster and its associated
        POSCAR.
    sigma:
        The smearing width. If set to None, the PDOS is not broadened.
    x_shift:
        An offset to apply along the x-axis, effectively shifting the Fermi level.
    y_shift:
        An offset to apply along the y-axis, separating the Fe and O PDOS for sake of
        plotting them side by side.

    Returns
    -------
    energies:
        The energies at which the PDOS has been evaluated.
    pdos:
        The PDOS, with rows ordered as follows: (fe_up, fe_down, o_up, o_down).
    """
    doscar = Doscar(doscar=f"{path}/DOSCAR.lobster", structure_file=f"{path}/POSCAR")
    structure = doscar._final_structure

    energies = doscar.energies + x_shift
    fe_pdos, o_pdos = [], []
    for idx, pdos in enumerate(doscar.pdos):
        if structure[idx].species_string == "Fe":
            fe_pdos.append(pdos)

        elif structure[idx].species_string == "O":
            o_pdos.append(pdos)

    num_energies = len(energies)
    pdos = np.zeros((4, num_energies))

    # Sum Fe and O PDOS over orbitals and atoms.
    for species_idx, species_pdos in enumerate((fe_pdos, o_pdos)):
        for atom_pdos in species_pdos:
            l = "d" if species_idx == 0 else "p"

            up_idx = 2 * species_idx
            down_idx = up_idx + 1

            pdos[up_idx] += np.sum(
                [
                    atom_pdos[orbital][Spin.up]
                    for orbital in atom_pdos.keys()
                    if l in orbital
                ],
                axis=0,
            )
            pdos[down_idx] -= np.sum(
                [
                    atom_pdos[orbital][Spin.down]
                    for orbital in atom_pdos.keys()
                    if l in orbital
                ],
                axis=0,
            )

    pdos[:2] += y_shift
    pdos[2:] -= y_shift

    if sigma is not None:
        gaussian_filter1d(pdos, sigma, output=pdos)

    return energies, pdos


def parse_icobi(path: str) -> dict[tuple[str, str], float]:
    """
    Parse the integrated crystal orbital bond indices (ICOBIs) for configurations of
    Li_xFeSiO4.

    Parameters
    ----------
    path:
        Filepath to ICOBILIST.lobster.

    Returns
    -------
    icobis:
        The parsed ICOBIs, summed over spin channels.

    Notes
    -----
    The dictionary returned by this function mirrors the structure of the the original
    ICOBILIST.lobster file, except that the spin channels have been summed:

    {
        ("Fe1", "O3"): 0.41621,
        ("Fe1", "O4"): 0.41692,
        ...
    }
    """
    with open(path, "r") as stream:
        lines = stream.readlines()

    heading_indices = []
    for idx, line in enumerate(lines):
        if "COBI#" in line:
            heading_indices.append(idx)

    start, end = heading_indices
    icobis = {}
    for idx, line in enumerate(lines[start + 1 : end]):
        up = line.split()
        down = lines[idx + 1 + end].split()

        symbols = tuple(sorted([up[1], up[2]]))
        icobi = float(up[-1]) + float(down[-1])

        icobis[symbols] = icobi

    return icobis
