# SPDX-FileCopyrightText: Copyright (C) 2026 Patrick J. Taylor
# SPDX-License-Identifier: GPL-3.0-or-later


class DependencyError(Exception):
    """
    A trivial custom exception to explicitly highlight the dependency of certain figure
    plotting scripts on the pre-calculation of Wannier-assigned MD oxidation states.
    """

    def __init__(self) -> None:
        message = "scripts/assign_md_oxidation_states.py must be run to regenerate Wannier-assigned oxidation states for selected AIMD frames."

        super().__init__(message)
