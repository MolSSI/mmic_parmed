from pydantic import Field, validator
from typing import Dict, Any, Optional
from mmic_translator.models.base import ToolkitModel
from mmelemental.models.molecule import Molecule
import parmed

# ParmEd converter components
from mmic_parmed.components.mol_component import ParmedToMolComponent
from mmic_parmed.components.mol_component import MolToParmedComponent


__all__ = ["ParmedMol"]


class ParmedMol(ToolkitModel):
    """ A model for ParmEd.Universe storing an MM molecule. """

    @property
    def dtype(self):
        """ Returns the fundamental molecule object type. """
        return parmed.structure.Structure

    @classmethod
    def isvalid(cls, data):
        """ Makes sure the Structure object stores atoms. """
        if hasattr(data, "atoms"):
            if len(data.atoms):
                return data

        raise ValueError("ParmEd molecule object does not contain any atoms!")

    @classmethod
    def from_file(
        cls, filename: str = None, top_filename: str = None, **kwargs
    ) -> "ParmedMol":
        """
        Constructs an ParmedMol object from file(s).

        Parameters
        ----------
        filename : str, optional
            The atomic positions filename to read
        top_filename: str, optional
            The topology filename to read
        **kwargs
            Any additional keywords to pass to the constructor
        Returns
        -------
        ParmedMol
            A constructed ParmedMol class.
        """

        kwargs.pop(
            "dtype", None
        )  # load_file doesn't seem to support specifying file formats

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
        cls, data: Molecule, version: Optional[str] = None, **kwargs: Dict[str, Any]
    ) -> "ParmedMol":
        """
        Constructs an ParmedMol object from an MMSchema Molecule object.
        Parameters
        ----------
        data: Molecule
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
        inputs = {"schema_object": data, "schema_version": version}
        out = MolToParmedComponent.compute(inputs)
        return cls(data=out.tk_object, units=out.tk_units)

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

    def to_schema(self, version: Optional[str] = None, **kwargs) -> Molecule:
        """Converts the molecule to MMSchema molecule.
        Parameters
        ----------
        version: str, optional
            Schema specification version to comply with e.g. 1.0.1.
        **kwargs
            Additional kwargs to pass to the constructor.
        """
        inputs = {"tk_object": self.data, "schema_version": version, "kwargs": kwargs}
        out = ParmedToMolComponent.compute(inputs)
        if version:
            assert version == out.schema_version
        return out.schema_object
