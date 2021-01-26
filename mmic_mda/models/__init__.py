from .mol import *
from .traj import *

# Need to update these lists
molread_ext_maps = {
    ".gro": "gro",
    ".pdb": "pdb",
    ".top": "top",
    ".psf": "psf",
}

molwrite_ext_maps = {".gro": "gro", ".pdb": "pdb"}

classes_map = {"Mol": MdaMol, "Traj": MdaTraj}
