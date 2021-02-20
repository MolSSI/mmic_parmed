from pydantic import Field, validator
from typing import Dict, Any, Optional
from mmelemental.models.base import ToolkitModel
from mmelemental.models.forcefield.mm_ff import ForceField
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

    @classmethod
    def isvalid(cls, data):
        """ Makes sure the Structure object stores atoms. """
        if hasattr(data, "atoms"):
            if len(data.atoms):
                return data

        raise ValueError("ParmEd forcefield object does not contain any atoms!")

    @classmethod
    @require("parmed")
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
        import parmed

        kwargs.pop(
            "dtype", None
        )  # load_file doesn't seem to support specifying file formats
        ff = parmed.load_file(filename=filename, **kwargs)

        return cls(data=ff)

    @classmethod
    def from_schema(
        cls, data: ForceField, version: Optional[str] = None, **kwargs: Dict[str, Any]
    ) -> "ParmedFF":
        """
        Constructs an ParmedFF object from an MMSchema ForceField object.
        Parameters
        ----------
        data: ForceField
            Data to construct the forcefield object from.
        version: str, optional
            Schema version e.g. 1.0.1
        **kwargs
            Additional kwargs to pass to the constructors.
        Returns
        -------
        ParmedFF
            A constructed ParmedFF object.
        """
        from mmic_parmed.components.ff_component import FFToParmedComponent

        return FFToParmedComponent.compute(data)

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

    def to_schema(self, version: Optional[str] = None, **kwargs) -> ForceField:
        """Converts the molecule to MMSchema molecule.
        Parameters
        ----------
        version: str, optional
            Schema specification version to comply with e.g. 1.0.1.
        **kwargs
            Additional kwargs to pass to the constructor.
        """
        from mmic_parmed.components.ff_component import ParmedToFFComponent

        return ParmedToFFComponent.compute(self)
