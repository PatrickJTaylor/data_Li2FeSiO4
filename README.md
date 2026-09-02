# Supporting Data for "Multi-Step Oxygen Redox Mechanism in the Polyanionic Lithium-Rich Cathode Li<sub>2</sub>FeSiO<sub>4</sub>"

[![License: GPL v3](https://img.shields.io/badge/license-GPLv3-blue)](https://www.gnu.org/licenses/gpl-3.0)
[![License: CC BY 4.0](https://img.shields.io/badge/license-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

Authors:
- Patrick J. Taylor, ORCID: [0009-0003-6511-6442](https://orcid.org/0009-0003-6511-6442)
- Kit McColl, ORCID: [0000-0002-7794-8276](https://orcid.org/0000-0002-7794-8276)
- M. Saiful Islam, ORCID: [0000-0002-8077-6241](https://orcid.org/0000-0002-8077-6241)
- Benjamin J. Morgan, ORCID: [0000-0002-3056-8233](https://orcid.org/0000-0002-3056-8233)

This is the Supporting Data for the manuscript "Multi-Step Oxygen Redox Mechanism in the Polyanionic Lithium-Rich Cathode Li<sub>2</sub>FeSiO<sub>4</sub>", accepted for publication in the *Journal of the American Chemical Society* (DOI to follow). Preprint: [10.26434/chemrxiv.10001983/v1](https://doi.org/10.26434/chemrxiv.10001983/v1).

More precisely, this repository contains the [minimal dataset](extracted_data) and accompanying [analysis code](scripts) required to regenerate all of the results and figures[^1] in the manuscript.

[^1]: Note that figures 1 and 4 were not generated programmatically, so they cannot be reproduced in an automated manner.
That being said, the raw data used to produce these figures are present in the dataset included in this archive.
The geometry of Li<sub>2</sub>FeSiO<sub>4</sub> required to plot figure 1, as well as the `xsf` files necessary to plot the Wannier isosurfaces of Figure 4 can all be found in [extracted_data/sequential_delithiation](data_Li2FeSiO4/tree/main/extracted_data/sequential_delithiation).

## Usage 📝

Practically speaking, to reproduce any one result or figure in the manuscript, all that is required is to run the relevant [Python script](scripts).
Every figure that can be programmatically generated has an associated script, which when run will deposit the regenerated pdf in [figures/regenerated](figures/regenerated).
There are also several scripts that do not regenerate figures, but reproduce other results that are discussed in the text of the manuscript:
- [assign_delithiation_oxidation_states.py](scripts/assign_delithiation_oxidation_states.py): assign Wannier oxidation states for sequentially delithiated Li<sub>$x$</sub>FeSiO<sub>4</sub> at $x = 2, 1, 0$.
- [thermodynamic_stability.py](scripts/thermodynamic_stability.py): determine the thermodynamic stability of FeSiO<sub>4</sub>.
- [fe_coordination_by_oxidation_state.py](scripts/fe_coordination_by_oxidation_state.py): Fe coordination environments by Wannier-assigned oxidation state, for the selected AIMD frames.
- [short_fe_fe_distances.py](scripts/short_fe_fe_distances.py): short Fe&ndash;Fe distances cross-referenced with Fe coordination and oxidation state.
- [trace_o2_forming_oxygens.py](scripts/trace_o2_forming_oxygens.py): O&ndash;O dimer formation and the Fe coordination of the dimer-forming O atoms.
- [dimer_oxygen_coordination.py](scripts/dimer_oxygen_coordination.py): full cation coordination shell of the peroxide-forming O atoms in the 500 K trajectory.

The relaxed AIMD frame energies in `extracted_data/molecular_dynamics/*/relaxed_frame_energies.npy` were extracted from the raw VASP output with [extract_relaxed_frame_energies.py](scripts/extract_relaxed_frame_energies.py); this script requires the raw dataset (University of Bath Research Data Archive, [10.15125/BATH-01647](https://doi.org/10.15125/BATH-01647)) and is included to document the extraction.

Whilst one could technically run the relevant scripts in any old Python environment that has access to the necessary dependencies, this repository provides several convenient methods for building bespoke and version-pinned environments for the sake of maximal reproducibility.

### Option 1 (recommended) - Docker 🐋

All figures and results can be regenerated via the [regen_results.sh](regen_results.sh) script, which uses [Docker](https://www.docker.com/resources/what-container/) to run all of the relevant Python scripts in a containerised manner.
For guidance on installing Docker, please see the official [documentation](https://docs.docker.com/engine/install/).

To regenerate everything in one go, simply run:

```shell
./regen_results.sh all
```

Executing the above will build a Docker image, spin up a container, run the full suite of analysis code, display the regenerated results and then copy the resulting figures to `./regenerated_figures`.

Rather than regenerating everything in one go, it is also possible to select a subset of figures/results to be reproduced - see the help menu for further details:

```
./regen_results.sh --help
```

### Option 2 - uv ⚡

This repository is technically structured as a [uv](https://docs.astral.sh/uv/) project, so in addition to the Docker workflow, it is also possible to run any of the analysis code in the following manner:

```shell
uv run scripts/figure_6.py
```

Running the analysis scripts through `uv` ensures that an appropriate Python environment is used (much like the Docker solution, minus the additional reproducibility afforded by containerisation).
To install uv, see the official [documentation](https://docs.astral.sh/uv/getting-started/installation/).

### Option 3 - Other Python environments 🐍

If you do not have Docker or uv installed and would rather stick with the Python tooling you already possess, you can (at the cost of some stability/reproducibility) navigate to the repository root and attempt to:

```shell
pip install -e .
```

This will effectively install the analysis code and its dependencies into the current Python environment, after which you should be able to run any of the relevant scripts directly:

```shell
python scripts/figure_S2.py
```

It should be noted that due to the strict version pinning in this repository's [pyproject.toml](pyproject.toml), it may be challenging for `pip` to resolve its dependency graph without any conflicts if the environment that you are using already has many other packages installed.
For this reason, it is advisable to start from a fresh virtual environment if you intend to follow this workflow.

## License 🖋

### Code ⌨

All of the source code included in this repository is licensed under the GNU General Public License, version 3.0 or (at your option) any later version (GPL-3.0-or-later): https://www.gnu.org/licenses/gpl-3.0

See `LICENSE` for the full license text.

### Data and figures 📈

All of the data and figures included in this repository (as well as figures that may be generated by execution of the source code) are licensed under the Creative Commons Attribution 4.0 International License (CC-BY-4.0): https://creativecommons.org/licenses/by/4.0/

See `LICENSE_CC-BY-4.0` for the full license text.
