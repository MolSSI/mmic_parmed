from typing import Dict, Any, Optional
from mmic_translator.models import ToolkitModel
from mmelemental.models import Molecule
import parmed

# ParmEd converter components
from mmic_parmed.components.mol_component import ParmedToMolComponent
from mmic_parmed.components.mol_component import MolToParmedComponent


__all__ = ["ParmedMol"]


class ParmedMol(ToolkitModel):
    """A model for ParmEd Structure storing a molecule."""

    @classmethod
    def engine(cls):
        return "parmed", parmed.__version__

    @classmethod
    def dtype(cls):
        """Returns the fundamental molecule object type."""
        return parmed.structure.Structure

    @classmethod
    def isvalid(cls, data):
        """Makes sure the Structure object stores atoms."""
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
        cls, data: Molecule, version: Optional[int] = None, **kwargs: Dict[str, Any]
    ) -> "ParmedMol":
        """
        Constructs an ParmedMol object from an MMSchema Molecule object.
        Parameters
        ----------
        data: Molecule
            Data to construct Molecule from.
        version: int, optional
            Schema version e.g. 1. Overwrites data.schema_version.
        **kwargs
            Additional kwargs to pass to the constructors.
        Returns
        -------
        ParmedMol
            A constructed ParmedMol class.
        """
        inputs = {
            "schema_object": data,
            "schema_version": version or data.schema_version,
            "schema_name": data.schema_name,
        }
        out = MolToParmedComponent.compute(inputs)
        return cls(data=out.data_object, data_units=out.data_units)

    def to_file(self, filename: str, dtype: str = None, mode: str = "w", **kwargs):
        """Writes the molecule to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
            File format
        **kwargs
            Additional kwargs to pass to the constructors. kwargs takes precedence over  data.
        """
        if dtype:
            kwargs["format"] = kwargs.get("format", dtype)
        if mode == "w":
            kwargs["overwrite"] = True
        elif mode == "a":
            kwargs["overwrite"] = False
        else:
            raise NotImplementedError(
                "File write mode can be either 'w' (write) or 'a' (append) for now."
            )

        self.data.save(filename, **kwargs)

    def to_schema(self, version: Optional[int] = 0, **kwargs) -> Molecule:
        """Converts the molecule to MMSchema molecule.
        Parameters
        ----------
        version: str, optional
            Schema specification version to comply with.
        **kwargs
            Additional kwargs to pass to the constructor.
        """
        inputs = {
            "data_object": self.data,
            "schema_name": kwargs.pop("schema_name", Molecule.default_schema_name),
            "schema_version": version,
            "keywords": kwargs,
        }
        out = ParmedToMolComponent.compute(inputs)
        if version:
            assert version == out.schema_version
        return out.schema_object
