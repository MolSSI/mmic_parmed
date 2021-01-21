from pydantic import Field
from typing import Dict, Any
from mmelemental.models.molecule.gen_mol import ToolkitMol

try:
    import MDAnalysis
except:
    raise ModuleNotFoundError('Make sure MDAnalysis is installed.')

class MdaMol(ToolkitMol):
    mol: MDAnalysis.Universe = Field(..., description = 'MDAnalysis molecule object.')

    @property
    def dtype(self):
        return 'mdanalysis'   

    @classmethod
    def build(cls, inputs: Dict[str, Any], dtype: str) -> "MdaMol":
        """
        Creates an instance of MdaMol object storing MDAnalysis.Universe. 
        This is done by parsing an input file (pdb, gro, ...).
        """
        if inputs.file:
            coords_fname = inputs.file.path
            if inputs.top_file:
                top_fname = inputs.top_file.path
            else:
                top_fname = None
            try:
                if top_fname:
                    mmol = MDAnalysis.Universe(top_fname, coords_fname)
                else:
                    mmol = MDAnalysis.Universe(coords_fname)
            except:
                raise ValueError(f"File type not supported: {inputs.file.ext}")

        elif inputs.code:
            raise NotImplementedError('No support for Chemical codes with MDAnalysis.')

        return cls(mol=mmol)
