import numpy as np
from numpy.typing import NDArray


def moving_average(data: NDArray[np.float64], window_size: int) -> NDArray[np.float64]:
    """
    Compute a moving average over a numpy array.

    Parameters
    ----------
    data:
        The numpy array to average.
    window_size:
        The size of the moving window used to average the data.

    Returns
    -------
    moving_average:
        The computed moving average.
    """
    moving_average = []

    idx = 0
    while idx < len(data) - window_size + 1:
        window = data[idx : idx + window_size]

        average = window.mean()

        moving_average.append(average)
        idx += 1

    return np.array(moving_average)
