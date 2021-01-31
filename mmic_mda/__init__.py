"""
mmic_mda
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
    ".pdb": "pdb",
    ".top": "top",
    ".psf": "psf",
}

molwrite_ext_maps = {".gro": "gro", ".pdb": "pdb"}

units = {
    "length": "angstrom",
    "time": "ps",
    "energy": "kJ/mol",
    "charge": "e",
    "speed": "angstrom/ps",
    "force": "kJ/(mol*angstrom)",
    "mass": "amu",
    "angle": "degrees",
}

_classes_map = {"Mol": models.MdaMol, "Traj": models.MdaTraj}
