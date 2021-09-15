"""
mmic_parmed.py
A short description of the project.

Handles the primary functions
"""


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
    ".inpcrd": "inpcrd",
}

molwrite_ext_maps = {
    ".gro": "gro",
    ".pdb": "pdb",
    ".inpcrd": "inpcrd",
}

ffread_ext_maps = {
    ".psf": "psf",
    ".top": "top",
    ".prm": "prm",
    ".prmtop": "prmtop",
}
ffwrite_ext_maps = {
    ".psf": "psf",
    ".top": "top",
    ".prm": "prm",
    ".prmtop": "prmtop",
}
