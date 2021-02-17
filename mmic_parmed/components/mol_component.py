from mmelemental.components.trans.template_component import TransComponent
from mmelemental.models.util.output import FileOutput
from mmelemental.models.molecule import Molecule
from mmic_parmed.models import ParmedMol
from typing import Dict, Any, List, Tuple, Optional
from mmelemental.util.decorators import require
from mmelemental.util.units import convert

__all__ = ["MolToParmedComponent", "ParmedToMolComponent"]


class MolToParmedComponent(TransComponent):
    """ A component for converting Molecule to ParmEd molecule object. """

    @classmethod
    def input(cls):
        return Molecule

    @classmethod
    def output(cls):
        return ParmedMol

    @require("parmed")
    def execute(
        self,
        inputs: Molecule,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, ParmedMol]:
        """
        Works for writing PDB files e.g. pmol.save("file.pdb") but fails for gro files
        TODO: need to investigate this more. Routine is also very slow. Try to vectorize.
        """
        import parmed

        pmol = parmed.structure.Structure()
        natoms = len(inputs.symbols)

        # Use dalton as the default unit for mass in ParmEd
        units = {"mass": "dalton"}

        if inputs.masses is not None:
            masses = convert(inputs.masses, inputs.masses_units, units["mass"])
        else:
            masses = [0] * natoms

        for index, symb in enumerate(inputs.symbols):

            name = inputs.names[index]
            # name = ToolkitMolecule.check_name(name)

            atomic_number = inputs.atomic_numbers[index]

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

            if inputs.residues:
                resname, resnum = inputs.residues[index]
            else:
                resname, resnum = "UNK", 0
            # classparmed.Residue(name, number=- 1, chain='', insertion_code='', segid='', list=None)[source]
            pmol.add_atom(atom, resname, resnum + 1, chain="", inscode="", segid="")

        if inputs.geometry is not None:
            pmol.coordinates = inputs.geometry.reshape(natoms, 3)
            pmol.coordinates = convert(
                pmol.coordinates, inputs.geometry_units, pmol.positions.unit.get_name()
            )
            units["length"] = pmol.positions.unit.get_name()

        if inputs.velocities is not None:
            units["speed"] = "angstrom/picosecond"
            pmol.velocities = inputs.velocities.reshape(natoms, 3)
            pmol.velocities = convert(
                pmol.velocities, inputs.velocities_units, units["speed"]
            )

        if inputs.connectivity:
            for i, j, btype in inputs.connectivity:
                pmol.atoms[i].bond_to(pmol.atoms[j])

        return True, ParmedMol(data=pmol, units=units)


class ParmedToMolComponent(TransComponent):
    """ A component for converting ParmEd molecule to Molecule object. """

    @classmethod
    def input(cls):
        return ParmedMol

    @classmethod
    def output(cls):
        return Molecule

    def execute(
        self,
        inputs: ParmedMol,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Molecule]:

        # I think parmed.Structure does not store forces
        pmol = inputs.data

        geo = TransComponent.get(pmol, "coordinates")
        geo_units = pmol.positions.unit.get_name() if geo is not None else None

        vel = TransComponent.get(pmol, "velocities")
        vel_units = "angstrom/picosecond" if vel is not None else None

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
            "names": names,
        }

        return True, Molecule(**input_dict)
