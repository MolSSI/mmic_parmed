from mmelemental.models import Molecule
from mmic.components import TacticComponent
from typing import List, Tuple, Optional, Set
from mmelemental.util.units import convert
from cmselemental.util.decorators import classproperty
from mmic_parmed.mmic_parmed import __version__
import parmed

from mmic_translator import (
    InputTrans,
    OutputTrans,
)

provenance_stamp = {
    "creator": "mmic_parmed",
    "version": __version__,
    "routine": __name__,
}

__all__ = ["MolToParmedComponent", "ParmedToMolComponent"]


class MolToParmedComponent(TacticComponent):
    """A component for converting Molecule to ParmEd molecule object."""

    @classproperty
    def input(cls):
        return InputTrans

    @classproperty
    def output(cls):
        return OutputTrans

    @classproperty
    def version(cls) -> str:
        """Finds program, extracts version, returns normalized version string.
        Returns
        -------
        str
            Return a valid, safe python version string.
        """
        ...

    @classproperty
    def strategy_comps(cls) -> Set[str]:
        """Returns the strategy component(s) this (tactic) component belongs to.
        Returns
        -------
        Set[str]
        """
        return {"mmic_translator"}

    def execute(
        self,
        inputs: InputTrans,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, OutputTrans]:
        """
        Works for writing PDB files e.g. pmol.save("file.pdb") but fails for gro files
        TODO: need to investigate this more. Routine is also very slow. Try to vectorize.
        """
        if isinstance(inputs, dict):
            inputs = self.input(**inputs)

        mmol = inputs.schema_object
        ndim = mmol.ndim

        if ndim != 3:
            raise NotImplementedError("mmic_parmed supports only 3D molecules.")

        pmol = parmed.structure.Structure()
        natoms = len(mmol.symbols)
        atom_empty = parmed.topologyobjects.Atom()

        if mmol.masses is not None:
            masses = convert(
                mmol.masses, mmol.masses_units, atom_empty.umass.unit.get_name()
            )
        else:
            masses = None

        if mmol.partial_charges is not None:
            charges = convert(
                mmol.partial_charges,
                mmol.partial_charges_units,
                atom_empty.ucharge.unit.get_symbol(),
            )
        else:
            charges = None

        if mmol.atom_labels is not None:
            atom_labels = mmol.atom_labels
        else:
            atom_labels = None

        if mmol.atomic_numbers is not None:
            atomic_numbers = mmol.atomic_numbers
        else:
            raise NotImplementedError(
                "mmic_parmed is supported only for atomic/molecular systems. Molecule.atomic_numbers must be defined."
            )

        for index, symb in enumerate(mmol.symbols):

            label = atom_labels[index] if atom_labels is not None else None
            # name = ToolkitMolecule.check_name(name)
            atomic_number = (
                atomic_numbers[index] if atomic_numbers is not None else None
            )
            mass = None if masses is None else masses[index]
            charge = None if charges is None else charges[index]

            # Will likely lose FF-related info ... but then Molecule is not supposed to store any params specific to FFs
            atom = parmed.topologyobjects.Atom(
                list=None,
                atomic_number=atomic_number,
                name=label or "",  # do not use symb?,
                type=label,
                mass=mass,
                charge=charge,
                nb_idx=0,
                solvent_radius=0.0,
                screen=0.0,
                tree="BLA",
                join=0.0,
                irotat=0.0,
                occupancy=1.0,
                bfactor=0.0,
                altloc="",
                number=-1,
                rmin=None,
                epsilon=None,
                rmin14=None,
                epsilon14=None,
            )

            if mmol.substructs is not None:
                resname, resnum = mmol.substructs[index]
            else:
                resname, resnum = "UNK", 0
            # classparmed.Residue(name, number=- 1, chain='', insertion_code='', segid='', list=None)[source]
            pmol.add_atom(atom, resname, resnum, chain="", inscode="", segid="")

        if mmol.geometry is not None:
            pmol.coordinates = mmol.geometry.reshape(natoms, ndim)
            pmol.coordinates = convert(
                pmol.coordinates, mmol.geometry_units, pmol.positions.unit.get_name()
            )

        if mmol.velocities is not None:
            units_speed = "angstrom/picosecond"  # hard-coded in parmed 3.4.0
            pmol.velocities = mmol.velocities.reshape(natoms, ndim)
            pmol.velocities = convert(
                pmol.velocities, mmol.velocities_units, units_speed
            )

        if mmol.connectivity is not None:
            for (
                i,
                j,
                order,
            ) in mmol.connectivity:
                # pmol.atoms[i].bond_to(pmol.atoms[j])
                pmol.bonds.append(
                    parmed.topologyobjects.Bond(pmol.atoms[i], pmol.atoms[j], order)
                )
                # both implementations seem to perform almost the same

        # if mmol.angles:
        #    for i, j, k in mmol.angles:
        #        pmol.angles.append(
        #            parmed.topologyobjects.Angle(
        #                pmol.atoms[i], pmol.atoms[j], pmol.atoms[k]
        #            )
        #        )

        # if mmol.dihedrals:
        #    for i, j, k, l in mmol.dihedrals:
        #        pmol.dihedrals.append(
        #            parmed.topologyobjects.Dihedral(
        #                pmol.atoms[i], pmol.atoms[j], pmol.atoms[k], pmol.atoms[l]
        #            )
        #        )

        return True, OutputTrans(
            schema_version=inputs.schema_version,
            schema_name=inputs.schema_name,
            proc_input=inputs,
            data_object=pmol,
            success=True,
            provenance=provenance_stamp,
        )


class ParmedToMolComponent(TacticComponent):
    """A component for converting ParmEd molecule to Molecule object."""

    @classproperty
    def input(cls):
        return InputTrans

    @classproperty
    def output(cls):
        return OutputTrans

    @classproperty
    def version(cls) -> str:
        """Finds program, extracts version, returns normalized version string.

        Returns
        -------
        str
            Return a valid, safe python version string.

        """
        ...

    @classproperty
    def strategy_comps(cls) -> Set[str]:
        """Returns the strategy component(s) this (tactic) component belongs to.

        Returns
        -------
        Set[str]

        """
        return {"mmic_translator"}

    def execute(
        self,
        inputs: InputTrans,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, OutputTrans]:

        if isinstance(inputs, dict):
            inputs = self.input(**inputs)

        # I think parmed.Structure does not store forces
        pmol = inputs.data_object
        geo_units = Molecule.default_units[
            "geometry_units"
        ]  # update later if specified
        vel_units = Molecule.default_units["velocities_units"]

        geo = getattr(pmol, "coordinates", None)
        if geo is not None:  # General enough? hackish?
            geo = geo.flatten()
            geo_units = pmol.positions.unit.get_name()

        vel = getattr(pmol, "velocities", None)
        if vel is not None:
            vel = vel.flatten()
            vel_units = "angstrom/picosecond"  # hard-coded in ParmEd

        atomic_nums = [atom.atomic_number for atom in pmol.atoms]
        names = [atom.name for atom in pmol.atoms]
        element_names = [atom.element_name for atom in pmol.atoms]

        masses = [atom.mass for atom in pmol.atoms]
        masses_units = pmol.atoms[0].umass.unit.get_name()

        charges = [atom.charge for atom in pmol.atoms]
        charges_units = pmol.atoms[0].ucharge.unit.get_symbol()

        # If bond order is none, set it to 1.
        if hasattr(pmol, "bonds"):
            connectivity = [
                (bond.atom1.idx, bond.atom2.idx, bond.order or 1) for bond in pmol.bonds
            ]
        else:
            connectivity = None

        if hasattr(pmol, "residues"):
            residues = [(atom.residue.name, atom.residue.idx) for atom in pmol.atoms]

        input_dict = {
            "atomic_numbers": atomic_nums,
            "symbols": element_names,
            "atom_labels": names,
            "geometry": geo,
            "geometry_units": geo_units,
            "velocities": vel,
            "velocities_units": vel_units,
            "substructs": residues,
            "connectivity": connectivity,
            "masses": masses,
            "masses_units": masses_units,
            "partial_charges": charges,
            "partial_charges_units": charges_units,
            "extras": inputs.keywords.get("extras"),
        }

        return True, OutputTrans(
            schema_version=inputs.schema_version,
            schema_name=inputs.schema_name,
            proc_input=inputs,
            schema_object=Molecule(**input_dict),
            success=True,
            provenance=provenance_stamp,
        )
