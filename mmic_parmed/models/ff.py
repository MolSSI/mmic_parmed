from pydantic import Field, validator
from typing import Dict, Any, Optional
from mmelemental.models.base import ToolkitModel
from mmelemental.models.molecule.mm_ff import ForceField
from mmelemental.util.decorators import require


__all__ = ["ParmedFF"]


class ParmedFF(ToolkitModel):
    """ A model for ParmEd.Universe storing an MM molecule. """

    @property
    @require("parmed")
    def dtype(self):
        """ Returns the fundamental molecule object type. """
        from parmed.structure import Structure

        return Structure

    @validator("data")
    def valid_mol(cls, data):
        """ Makes sure the Structure object stores atoms. """
        if hasattr(data, "atoms"):
            if len(data.atoms):
                return data

        raise ValueError("ParmEd molecule object does not contain any atoms!")

    @classmethod
    @require("parmed")
    def from_file(
        cls, filename: str = None, top_filename: str = None, dtype: str = None, **kwargs
    ) -> "ParmedFF":
        """
        Constructs an ParmedMol object from file(s).

        Parameters
        ----------
        filename : str, optional
            The atomic positions filename to read
        top_filename: str, optional
            The topology filename to read
        dtype: str, optional
            The type of file to interpret. If unset, ParmEd attempts to discover dtype from the file extension.
        **kwargs
            Any additional keywords to pass to the constructor
        Returns
        -------
        ParmedMol
            A constructed ParmedMol class.
        """
        import parmed

        if dtype:
            raise ValueError(f"This argument is not supported: dtype = {dtype}.")

        if filename and top_filename:
            mol = parmed.load_file(filename=top_filename, xyz=filename, **kwargs)
        elif filename:
            mol = parmed.load_file(filename, **kwargs)
        elif top_filename:
            mol = parmed.load_file(top_filename, **kwargs)
        else:
            raise TypeError(
                "You must supply at least one of the following: filename or top_filename."
            )

        return cls(data=mol)

    @classmethod
    def from_schema(
        cls,
        data: Mol,
        version: Optional[str] = None,
        **kwargs: Dict[str, Any],
    ) -> "ParmedMol":
        """
        Constructs an ParmedMol object from an MMSchema Mol object.
        Parameters
        ----------
        data: Mol
            Data to construct Molecule from.
        version: str, optional
            Schema version e.g. 1.0.1
        **kwargs
            Additional kwargs to pass to the constructors.
        Returns
        -------
        ParmedMol
            A constructed ParmedMol class.
        """
        from mmic_parmed.components.mol_component import MolToParmedComponent

        return MolToParmedComponent.compute(data)

    def to_file(self, filename: str, dtype: str = None, **kwargs):
        """Writes the molecule to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
            File format
        **kwargs
            Additional kwargs to pass to the constructors.
        """
        if dtype:
            kwargs["format"] = dtype
        self.data.save(filename, **kwargs)

    def to_schema(self, version: Optional[str] = None, **kwargs) -> Mol:
        """Converts the molecule to MMSchema molecule.
        Parameters
        ----------
        version: str, optional
            Schema specification version to comply with e.g. 1.0.1.
        **kwargs
            Additional kwargs to pass to the constructor.
        """
        from mmic_parmed.components.mol_component import ParmedToMolComponent

        return ParmedToMolComponent.compute(self)
