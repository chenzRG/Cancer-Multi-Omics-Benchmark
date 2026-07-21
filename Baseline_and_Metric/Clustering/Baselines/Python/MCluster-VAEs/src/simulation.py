"""
simulation.py — stub for SimulatorMorgane.

The original simulation.py is not part of this release.
MCluster-VAEs on MLOmics data does not use simulation mode,
so this stub satisfies the import in src/__init__.py without
breaking run_mlomics.py.
"""


class SimulatorMorgane:
    """Stub — raises NotImplementedError if actually instantiated."""

    def __init__(self, seed=0):
        raise NotImplementedError(
            "SimulatorMorgane simulation data is not available in this installation. "
            "Use a real dataset (e.g. run run_mlomics.py) instead."
        )

    def get_data(self, data_index, n_mult):
        raise NotImplementedError("SimulatorMorgane is not available.")
