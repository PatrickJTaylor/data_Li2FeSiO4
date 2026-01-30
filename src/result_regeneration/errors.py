class DependencyError(Exception):
    def __init__(self) -> None:
        message = "scripts/assign_md_oxidation_states.py must be run to regenerate Wannier-assigned oxidation states for selected AIMD frames."

        super().__init__(message)
