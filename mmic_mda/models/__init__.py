from .mol import *
from .traj import *

molread_ext_maps = {
    ".gro": "gro",
    ".pdb": "pdb",
    ".top": "top",
    ".psf": "psf",
}

molwrite_ext_maps = {
    ".gro": "gro",
    ".pdb": "pdb",
    ".top": "top",
    ".psf": "psf",
}

classes_map = {"Mol": MdaMol, "Traj": MdaTraj}
