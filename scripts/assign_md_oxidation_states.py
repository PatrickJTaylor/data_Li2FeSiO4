import numpy as np

from figure_regeneration.paths import EXTRACTED_DATA
from figure_regeneration.wannier import assign_wannier_centres

prefix = EXTRACTED_DATA / "molecular_dynamics"

valence_counts = {"Fe": 8, "Si": 4, "O": 6}

temperatures = ("500K", "1000K")
frame_indices = ((0, 6881, 16602, 39366), (0, 2581, 8315, 12276))
n_atoms = 216

valence = np.zeros(n_atoms, dtype=np.float64)
valence[:36] += valence_counts["Fe"]
valence[36:180] += valence_counts["O"]
valence[180:] += valence_counts["Si"]

for idx, temperature in enumerate(temperatures):
    for frame_idx in frame_indices[idx]:
        electrons = assign_wannier_centres(
            f"{prefix}/{temperature}/selected_frames/{frame_idx}"
        )
        oxidation_states = valence - electrons

        np.save(
            f"{prefix}/{temperature}/selected_frames/{frame_idx}/oxidation_states.npy",
            oxidation_states,
        )
