from typing import Dict, Any, Optional
from mmic_translator.models import ToolkitModel
from mmelemental.models.forcefield import ForceField
import parmed

from mmic_parmed.components.ff_component import FFToParmedComponent
from mmic_parmed.components.ff_component import ParmedToFFComponent


__all__ = ["ParmedFF"]


class ParmedFF(ToolkitModel):
    """A model for ParmEd Structure storing FF object."""

    @classmethod
    def engine(cls):
        return "parmed", parmed.__version__

    @classmethod
    def dtype(cls):
        """Returns the fundamental FF object type."""
        return parmed.structure.Structure

    @classmethod
    def isvalid(cls, data):
        """Makes sure the Structure object stores atoms."""
        if hasattr(data, "atoms"):
            if len(data.atoms):
                return data
        raise ValueError("ParmEd forcefield object does not contain any atoms!")

    @classmethod
    def from_file(cls, filename: str, **kwargs) -> "ParmedFF":
        """
        Constructs an ParmedFF object from file(s).

        Parameters
        ----------
        filename : str
            The forcefield filename to read
        **kwargs
            Any additional keywords to pass to the constructor
        Returns
        -------
        ParmedFF
            A constructed ParmedFF object.
        """
        kwargs.pop(
            "dtype", None
        )  # load_file doesn't seem to support specifying file formats
        ff = parmed.load_file(filename=filename, **kwargs)

        return cls(data=ff)

    @classmethod
    def from_schema(
        cls, data: ForceField, version: Optional[int] = None, **kwargs: Dict[str, Any]
    ) -> "ParmedFF":
        """
        Constructs an ParmedFF object from an MMSchema ForceField object.
        Parameters
        ----------
        data: ForceField
            Data to construct the forcefield object from.
        version: int, optional
            Schema version e.g. 1. Overwrites data.schema_version.
        **kwargs
            Additional kwargs to pass to the constructors.
        Returns
        -------
        ParmedFF
            A constructed ParmedFF object.
        """
        inputs = {
            "schema_object": data,
            "schema_version": version or data.schema_version,
            "schema_name": data.schema_name,
        }
        out = FFToParmedComponent.compute(inputs)
        return cls(data=out.data_object, units=out.data_units)

    def to_file(self, filename: str, dtype: str = None, **kwargs):
        """Writes the forcefield to a file.
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

    def to_schema(self, version: Optional[int] = 0, **kwargs) -> ForceField:
        """Converts the forcefield to MMSchema ForceField object.
        Parameters
        ----------
        version: int, optional
            Schema specification version to comply with e.g. 1.
        **kwargs
            Additional kwargs to pass to the constructor.
        """
        inputs = {
            "data_object": self.data,
            "schema_version": version,
            "schema_name": kwargs.pop("schema_name", ForceField.default_schema_name),
            "keywords": kwargs,
        }
        out = ParmedToFFComponent.compute(inputs)
        if version:
            assert version == out.schema_version
        return out.schema_object
