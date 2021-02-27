from mmelemental.models.util.output import FileOutput
from mmelemental.models.molecule import Molecule
from typing import Dict, Any, List, Tuple, Optional
from mmelemental.util.units import convert
import parmed

from mmic_translator import TransComponent
from mmic_translator.models.mol_io import (
    MolInToSchema,
    MolInFromSchema,
    MolOutToSchema,
    MolOutFromSchema,
)

__all__ = ["MolToParmedComponent", "ParmedToMolComponent"]


class MolToParmedComponent(TransComponent):
    """ A component for converting Molecule to ParmEd molecule object. """

    @classmethod
    def input(cls):
        return MolInFromSchema

    @classmethod
    def output(cls):
        return MolOutFromSchema

    def execute(
        self,
        inputs: Molecule,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, MolOutFromSchema]:
        """
        Works for writing PDB files e.g. pmol.save("file.pdb") but fails for gro files
        TODO: need to investigate this more. Routine is also very slow. Try to vectorize.
        """
        if isinstance(inputs, dict):
            inputs = self.input()(**inputs)

        mmol = inputs.schema_object
        pmol = parmed.structure.Structure()
        natoms = len(mmol.symbols)

        # Use dalton as the default unit for mass in ParmEd
        units = {"mass": "dalton"}

        if mmol.masses is not None:
            masses = convert(mmol.masses, mmol.masses_units, units["mass"])
        else:
            masses = [0] * natoms

        for index, symb in enumerate(mmol.symbols):

            name = mmol.atom_labels[index]
            # name = ToolkitMolecule.check_name(name)

            atomic_number = mmol.atomic_numbers[index]

            # Will likely lose FF-related info ... but then Molecule is not supposed to store any params specific to FFs
            atom = parmed.topologyobjects.Atom(
                list=None,
                atomic_number=atomic_number,
                name=name,
                type=symb,
                mass=masses[index],
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

            if mmol.residues:
                resname, resnum = mmol.residues[index]
            else:
                resname, resnum = "UNK", 0
            # classparmed.Residue(name, number=- 1, chain='', insertion_code='', segid='', list=None)[source]
            pmol.add_atom(atom, resname, resnum + 1, chain="", inscode="", segid="")

        if mmol.geometry is not None:
            pmol.coordinates = mmol.geometry.reshape(natoms, 3)
            pmol.coordinates = convert(
                pmol.coordinates, mmol.geometry_units, pmol.positions.unit.get_name()
            )
            units["length"] = pmol.positions.unit.get_name()

        if mmol.velocities is not None:
            units["speed"] = "angstrom/picosecond"
            pmol.velocities = mmol.velocities.reshape(natoms, 3)
            pmol.velocities = convert(
                pmol.velocities, mmol.velocities_units, units["speed"]
            )

        if mmol.connectivity:
            for (
                i,
                j,
                _,
            ) in mmol.connectivity:  # _: bond order ... can we use it in ParmEd?
                # pmol.atoms[i].bond_to(pmol.atoms[j])
                pmol.bonds.append(
                    parmed.topologyobjects.Bond(pmol.atoms[i], pmol.atoms[j])
                )
                # both implementations seem to perform almost the same

        if mmol.angles:
            for i, j, k in mmol.angles:
                pmol.angles.append(
                    parmed.topologyobjects.Angle(
                        pmol.atoms[i], pmol.atoms[j], pmol.atoms[k]
                    )
                )

        if mmol.dihedrals:
            for i, j, k, l in mmol.dihedrals:
                pmol.dihedrals.append(
                    parmed.topologyobjects.Dihedral(
                        pmol.atoms[i], pmol.atoms[j], pmol.atoms[k], pmol.atoms[l]
                    )
                )

        return True, MolOutFromSchema(tk_object=pmol, tk_units=units)


class ParmedToMolComponent(TransComponent):
    """ A component for converting ParmEd molecule to Molecule object. """

    @classmethod
    def input(cls):
        return MolInToSchema

    @classmethod
    def output(cls):
        return MolOutToSchema

    def execute(
        self,
        inputs: MolInToSchema,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, MolOutToSchema]:

        if isinstance(inputs, dict):
            inputs = self.input()(**inputs)

        # I think parmed.Structure does not store forces
        pmol = inputs.tk_object
        geo_units, vel_units = None, None

        geo = TransComponent.get(pmol, "coordinates")
        if geo is not None:  # General enough? hackish?
            geo = geo.flatten()
            geo_units = pmol.positions.unit.get_name()

        vel = TransComponent.get(pmol, "velocities")
        if vel is not None:
            vel = vel.flatten()
            vel_units = "angstrom/picosecond"  # hard-coded in ParmEd

        atomic_nums = [atom.atomic_number for atom in pmol.atoms]
        names = [atom.name for atom in pmol.atoms]

        masses = [atom.mass for atom in pmol.atoms]
        masses_units = pmol.atoms[0].umass.unit.get_name()

        # If bond order is none, set it to 1.
        if hasattr(pmol, "bonds"):
            connectivity = [
                (bond.atom1.idx, bond.atom2.idx, bond.order or 1) for bond in pmol.bonds
            ]
        else:
            connectivity = None

        if hasattr(pmol, "residues"):
            residues = [
                (atom.residue.name, atom.residue.idx + 1) for atom in pmol.atoms
            ]

        input_dict = {
            "atomic_numbers": atomic_nums,
            "geometry": geo,
            "geometry_units": geo_units,
            "velocities": vel,
            "velocities_units": vel_units,
            "residues": residues,
            "connectivity": connectivity,
            "masses": masses,
            "masses_units": masses_units,
            "atom_labels": names,
        }

        out = MolOutToSchema(schema_object=Molecule(**input_dict))
        return True, out
