"""
mmic_parmed
A short description of the project.
"""

# Add imports here
from . import components, models

# Handle versioneer
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions

# Need to update these lists
molread_ext_maps = {
    ".gro": "gro",
    ".psf": "psf",
    ".pdb": "pdb",
    ".sdf": "sdf",
    ".mol": "mol",
    ".mol2": "mol2",
}

molwrite_ext_maps = {".gro": "gro", ".pdb": "pdb"}

ffread_ext_maps = {".psf": "psf", ".top": "top", ".prm": "prm"}

_classes_map = {
    "Molecule": models.ParmedMol,
    # "Trajectory": models.ParmedTraj,
    "ForceField": models.ParmedFF,
}
