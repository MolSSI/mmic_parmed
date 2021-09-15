"""
mmic_parmed
A short description of the project.
"""

# Add imports here
from . import components, models
from .mmic_parmed import (
    molread_ext_maps,
    molwrite_ext_maps,
    ffread_ext_maps,
    ffwrite_ext_maps,
    __version__,
)


_classes_map = {
    "Molecule": models.ParmedMol,
    # "Trajectory": models.ParmedTraj,
    "ForceField": models.ParmedFF,
}
