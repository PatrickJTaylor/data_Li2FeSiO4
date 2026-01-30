# Supporting Data for "Multi-Step Oxygen Redox Mechanism in the Polyanionic Lithium-Rich Cathode Li<sub>2</sub>FeSiO<sub>4</sub>"

Authors:
- Patrick J. Taylor, ORCID: [0009-0003-6511-6442](https://orcid.org/0009-0003-6511-6442)
- Kit McColl, ORCID: [0000-0002-7794-8276](https://orcid.org/0000-0002-7794-8276)
- M. Saiful Islam, ORCID: [0000-0002-8077-6241](https://orcid.org/0000-0002-8077-6241)
- Benjamin J. Morgan, ORCID: [0000-0002-3056-8233](https://orcid.org/0000-0002-3056-8233)

This is the Supporting Data for the manuscript "Multi-Step Oxygen Redox Mechanism in the Polyanionic Cathode Li<sub>2</sub>FeSiO<sub>4</sub>".
More precisely, this repository contains the [minimal dataset](extracted_data) and accompanying [analysis code](scripts) required to regenerate all of the results and figures[^1] in the manuscript.

[^1]: Note that figures 1 and 4 were not generated programmatically, so they cannot be reproduced in an automated manner.
That being said, the raw data used to produce these figures are present in the dataset included in this archive.
The geometry of Li<sub>2</sub>FeSiO<sub>4</sub> required to plot figure 1, as well as the `xsf` files necessary to plot the Wannier isosurfaces of Figure 4 can all be found in [extracted_data/sequential_delithiation](extracted_data/sequential_delithiation).

Figure 1 was produced in [VESTA]() based on the [pristine geometry](extracted_data/sequential_delithiation/Li2FeSiO4/POSCAR) of Li<sub>2</sub>FeSiO<sub>4</sub>.
Figure 4 was produced 

## Usage 📝

Practically speaking, to reproduce any one result or figure in the manuscript, all that is required is to run the relevant [Python script](scripts).
Every figure that can be programmatically generated has an associated script, which when run will deposit the regenerated pdf in [figures/regenerated](figures/regenerated).
There are also several scripts that do not regenerate figures, but reproduce other results that are discussed in the text of the manuscript:
- [assign_delithiation_oxidation_states.py](scripts/assign_delithiation_oxidation_states.py): assign Wannier oxidation states for sequentially delithiated Li<sub>$x$</sub>FeSiO<sub>4</sub> at $x = 2, 1, 0$.
- [thermodynamic_stability.py](scripts/thermodynamic_stability.py): determine the thermodynamic stability of FeSiO<sub>4</sub>.

Whilst one could run the relevant scripts in any old Python environment that has access to the necessary dependencies, this repository provides two convenient methods for building bespoke and version-pinned environments for the sake of maximal reproducibility.

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

This repository is technically structured as a [uv](https://docs.astral.sh/uv/) project, so in addition to the Docker workflow, it is also possible to directly run any of the analysis code in the following manner:

```shell
uv run scripts/figure_6.py
```

Running the analysis scripts through `uv` ensures that an appropriate Python environment is used (much like the Docker solution, minus the additional reproducibility afforded by containerisation).
To install `uv`, see the official [documentation](https://docs.astral.sh/uv/getting-started/installation/).
