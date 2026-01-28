import math

from matplotlib import rcParams

golden_ratio = (1 + math.sqrt(5)) / 2
fig_width = 426 / 72
fig_height = fig_width / golden_ratio

rc_parameters = {
    "backend": "pdf",
    "figure.figsize": (fig_width, fig_height),
    "font.family": "sans-serif",
    "font.size": 10,
    "legend.facecolor": "white",
    "legend.fancybox": False,
    "legend.frameon": False,
    "savefig.bbox": "tight",
}


def set_formatting(fig_size: tuple[float, float] | None = None) -> None:
    """
    Adjust matplotlib's rcParams to best fit the manuscript, with a given figure size.

    Parameters
    ----------
    fig_size:
        The desired figure size.
    """
    if fig_size is not None:
        rc_parameters["figure.figsize"] = fig_size

    for rc_parameter, value in rc_parameters.items():
        rcParams[rc_parameter] = value
