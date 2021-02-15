from mmelemental.components.trans import TransComponent
from mmelemental.models.util.output import FileOutput
from mmelemental.models.forcefield.mm_ff import ForceField
from mmic_parmed.models import ParmedFF
from typing import Dict, Any, List, Tuple, Optional
from mmelemental.util.decorators import require
from mmelemental.util.units import convert

__all__ = ["FFToParmedComponent", "ParmedToFFComponent"]


class FFToParmedComponent(TransComponent):
    """ A component for converting Molecule to ParmEd molecule object. """

    @classmethod
    def input(cls):
        return ForceField

    @classmethod
    def output(cls):
        return ParmedFF

    @require("parmed")
    def execute(
        self,
        inputs: ForceField,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, ParmedFF]:
        """
        Works for writing PDB files e.g. pmol.save("file.pdb") but fails for gro files
        TODO: need to investigate this more. Routine is also very slow. Try to vectorize.
        """
        import parmed
        import mmic_parmed

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
            pmol.velocities = inputs.velocities.reshape(natoms, 3)
            pmol.velocities = convert(
                pmol.velocities,
                inputs.velocities_units,
                mmic_parmed.units["speed"],
            )
            units["speed"] = mmic_parmed.units["speed"]

        if inputs.connectivity:
            for i, j, btype in inputs.connectivity:
                pmol.atoms[i].bond_to(pmol.atoms[j])

        return True, ParmedFF(data=pmol, units=units)


class ParmedToFFComponent(TransComponent):
    """ A component for converting ParmEd molecule to Molecule object. """

    @classmethod
    def input(cls):
        return ParmedFF

    @classmethod
    def output(cls):
        return ForceField

    def execute(
        self,
        inputs: ParmedFF,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, ForceField]:

        ff = inputs.data
        atom = ff.atoms[0]
        bond = ff.bonds[0] if len(ff.bonds) else None
        angle = ff.angles[0] if len(ff.angles) else None
        dihedral = ff.dihedrals[0] if len(ff.dihedrals) else None

        data = [(atom.name, atom.charge) for atom in ff.atoms]
        names, charges = zip(*data)
        charges = convert(
            charges, atom.ucharge.unit.get_symbol(), mmic_parmed.units["charge"]
        )

        rmin_unit = atom.urmin.unit.get_name()
        rmin_14_factor = convert(
            1.0, atom.urmin_14.unit.get_name(), mmic_parmed.units["length"]
        )
        uepsilon_factor = convert(
            1.0, atom.uepsilon.unit.get_name(), mmic_parmed.units["energy"]
        )
        epsilon_14_factor = convert(
            1.0, atom.uepsilon_14.unit.get_name(), mmic_parmed.units["energy"]
        )
        params = [
            (
                atom.rmin * rmin_factor,
                atom.rmin_14 * rmin_14_factor,
                atom.epsilon * uepsilon_factor,
                atom.epsilon_14 * epsilon_14_factor,
            )
            for atom in ff.atoms
        ]

        if bond:
            bond_req_factor = convert(
                1.0, bond.type.ureq.unit.get_name(), mmic_parmed.units["length"]
            )
            bond_k_factor = convert(
                1.0, bond.type.uk.unit.get_name(), mmic_parmed.units["energy"]
            )
            bonds = [
                (bond.type.req * bond_req_factor, bond.type.k * bond_k_factor)
                for bond in ff.bonds
            ]
            bonds_type = [bond.funct for bond in ff.bonds]

        if angle:
            angle_theta_factor = convert(
                1.0, angle.type.utheteq.unit.get_name(), mmic_parmed.units["angle"]
            )
            angle_k_factor = convert(
                1.0, angle.type.uk.unit.get_name(), mmic_parmed.units["energy"]
            )
            angles = [(angle.type.theteq, angle.type.k) for angle in ff.angles]
            angles_type = [angle.funct for angle in ff.angles]

        if dihedral:
            # Why the hell does parmed interpret dihedral.type as a list sometimes? See dialanine.top, last 8 dihedral entries. WEIRD!
            # This is a hack around this
            dtype = dihedral.type[0] if isinstance(dihedral.type, list) else dihedral.type

            dihedral_phi_factor = convert(
                1.0, dtype.uphi_k.unit.get_name(), mmic_parmed.units["energy"]
            )
            dihedral_phase_factor = convert(
                1.0, dtype.uphase.unit.get_name(), mmic_parmed.units["degree"]
            )

            # Need to look into per, scee, and scnb ... their meaning, units, etc.
            dihedrals = [
                (
                    dihedral.type[0].phase * dihedral_phase_factor ,
                    dihedral.type[0].phi_k * dihedral_phi_factor,
                    dihedral.type[0].per,
                    dihedral.type[0].scee,
                    dihedral.type[0].scnb,
                )
                if isinstance(dihedral.type, list)
                else (
                    dihedral.type.phase * dihedral_phase_factor,
                    dihedral.type.phi_k * dihedral_phi_factor,
                    dihedral.type.per,
                    dihedral.type.scee,
                    dihedral.type.scnb,
                )
                for dihedral in ff.dihedrals
            ]
            dihedrals_type = [dihedral.funct for dihedral in ff.dihedrals]

        improper_dihedrals = None
        charge_groups = None
        exclusions = ff.nrexcl
        inclusions = None

        input_dict = {
            "bonds": bonds,
            "angles": angles,
            "dihedrals": dihedrals,
            "improper_dihedrals": improper_dihedrals,
            "charges": charges,
            "charge_groups": charge_groups,
            "nonbonded": params,
            "exclusions": exclusions,
            "inclusions": inclusions,
            "types": names,
        }

        return True, ForceField(**input_dict)
