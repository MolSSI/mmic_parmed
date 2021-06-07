from mmelemental.models.molecule import Molecule
from typing import List, Tuple, Optional
from mmelemental.util.units import convert
import parmed

from mmic_translator import (
    TransComponent,
    TransInput,
    TransOutput,
)

__all__ = ["MolToParmedComponent", "ParmedToMolComponent"]


class MolToParmedComponent(TransComponent):
    """A component for converting Molecule to ParmEd molecule object."""

    def execute(
        self,
        inputs: TransInput,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, TransOutput]:
        """
        Works for writing PDB files e.g. pmol.save("file.pdb") but fails for gro files
        TODO: need to investigate this more. Routine is also very slow. Try to vectorize.
        """
        if isinstance(inputs, dict):
            inputs = self.input()(**inputs)

        mmol = inputs.schema_object
        pmol = parmed.structure.Structure()
        natoms = len(mmol.symbols)
        atom_empty = parmed.topologyobjects.Atom()

        if mmol.masses is not None:
            masses = convert(
                mmol.masses, mmol.masses_units, atom_empty.umass.unit.get_symbol()
            )
        else:
            masses = None

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
            mass = masses[index] if masses is not None else None

            # Will likely lose FF-related info ... but then Molecule is not supposed to store any params specific to FFs
            atom = parmed.topologyobjects.Atom(
                list=None,
                atomic_number=atomic_number,
                name=label or "",  # do not use symb?,
                type=label,
                mass=mass,
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

            if mmol.substructs:
                resname, resnum = mmol.substructs[index]
            else:
                resname, resnum = "UNK", 0
            # classparmed.Residue(name, number=- 1, chain='', insertion_code='', segid='', list=None)[source]
            pmol.add_atom(atom, resname, resnum, chain="", inscode="", segid="")

        if mmol.geometry is not None:
            pmol.coordinates = mmol.geometry.reshape(natoms, 3)
            pmol.coordinates = convert(
                pmol.coordinates, mmol.geometry_units, pmol.positions.unit.get_name()
            )

        if mmol.velocities is not None:
            units_speed = "angstrom/picosecond"  # hard-coded in parmed 3.4.0
            pmol.velocities = mmol.velocities.reshape(natoms, 3)
            pmol.velocities = convert(
                pmol.velocities, mmol.velocities_units, units_speed
            )

        if mmol.connectivity:
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

        return True, TransOutput(proc_input=inputs, data_object=pmol)


class ParmedToMolComponent(TransComponent):
    """A component for converting ParmEd molecule to Molecule object."""

    def execute(
        self,
        inputs: TransInput,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, TransOutput]:

        if isinstance(inputs, dict):
            inputs = self.input()(**inputs)

        # I think parmed.Structure does not store forces
        pmol = inputs.data_object
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
        element_names = [atom.element_name for atom in pmol.atoms]

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
        }

        return True, TransOutput(
            proc_input=inputs, schema_object=Molecule(**input_dict)
        )
