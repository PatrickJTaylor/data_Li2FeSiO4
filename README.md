# Supporting Data for "Multi-Step Oxygen Redox Mechanism in the Polyanionic Lithium-Rich Cathode Li<sub>2</sub>FeSiO<sub>4</sub>"

This is the Supporting Data for the manuscript "Multi-Step Oxygen Redox Mechanism in the Polyanionic Cathode Li<sub>2</sub>FeSiO<sub>4</sub>".

## Usage

This archive provides two methods by which the figures/results discussed in the manuscript can be regenerated from the underlying data.

### Option 1 (recommended) - Docker

For maximal reproducibility, all figures/results can be regenerated via the `regen_results.sh` script, which uses Docker to run all the relevant Python scripts in a containerised manner.

To regenerate everything in one go:

```shell
./regen_results.sh all
```

Running the above will build a Docker image, spin up a container, run all the necessary Python and then copy any resulting figures to `./regenerated_figures`.

For additional options, see the help menu built into the bash script:

```
./regen_results.sh --help
```

### Option 2 - uv

All of the Python scripts required for result regeneration are located in `./scripts`.
For example, every figure that can be programmatically regenerated has a corresponding `figure_x.py` script in `./scripts` - there are also other scripts that do not generate figures, but do reproduce other results discussed in the text of the manuscript (e.g., the energy of FeSiO<sub>4</sub> above the hull).
Technically these scripts can be run with any Python environment that contains the required packages, but for sake of reproducibility, it is better to use `uv`:

```
uv run scripts/figure_2.py
```

The lockfile that ships with this archive will ensure that the above is executed with a highly reproducible Python environment, which should mitigate the vast majority of possible problems with breaking changes in dependencies etc.
