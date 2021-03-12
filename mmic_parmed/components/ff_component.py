from mmelemental.models.util.output import FileOutput
from mmelemental.models import forcefield
from mmic_translator.models.io import (
    TransInput,
    TransOutput,
)
from mmic_translator.components import TransComponent
from typing import Dict, Any, List, Tuple, Optional
from mmelemental.util.units import convert
import parmed
from enum import Enum

__all__ = ["FFToParmedComponent", "ParmedToFFComponent"]

bondTypes = {
    1: forcefield.bonded.bonds.potentials.Harmonic,
    2: forcefield.bonded.bonds.potentials.Gromos96,
}


class FFToParmedComponent(TransComponent):
    """ A component for converting Molecule to ParmEd molecule object. """

    @classmethod
    def input(cls):
        return TransInput

    @classmethod
    def output(cls):
        return TransOutput

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

        ff = parmed.structure.Structure()
        params = parmed.parameters.ParameterSet()

        # Use dalton as the default unit for mass in ParmEd
        units = {"mass": "dalton"}
        mm_ff = inputs.schema_object
        natoms = len(mm_ff.symbols)
        masses = mm_ff.masses

        rmin_factor = convert(1.0, atom.urmin.unit.get_name(), mm_ff.bonds.lengths_units)
        rmin_14_factor = convert(
            1.0, atom.urmin_14.unit.get_name(), mm_ff.nonbonded.sigma_units
        )
        epsilon_factor = convert(
            1.0, atom.uepsilon.unit.get_name(), mm_ff.nonbonded.epsilon_units
        )
        epsilon_14_factor = convert(
            1.0, atom.uepsilon_14.unit.get_name(), mm_ff.nonbonded.epsilon_units
        )
        scaling_factor = 2 ** (1.0 / 6.0)  # rmin = 2^(1/6) sigma
        sigma = mm_ff.nonbonded.sigma * rmin_factor * scaling_factor
        sigma_14 = mm_ff.nonbonded.sigma14 * rmin_14_factor * scaling_factor
        epsilon = mm_ff.nonbonded.epsilon * epsilon_factor
        epsilon_14 = mm_ff.nonbonded.epsilon * epsilon_14_factor
        sigma, sigma_14, epsilon, epsilon_14 = zip(*params)

        for index, symbol in enumerate(mm_ff.symbols):

            atype = AtomType(
                name=name, number=idx, mass=mass, atomic_number=atomic_number
            )

            # atomic_number = mm_ff.atomic_numbers[index]
            type = mm_ff.types[index]

            # Will likely lose FF-related info ... but then Molecule is not supposed to store any params specific to FFs
            atom = parmed.topologyobjects.Atom(
                list=None,
                # atomic_number=atomic_number,
                name=type,
                type=symbol,
                # mass=masses[index],
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

        return True, TransOutput(tk_object=pmol, units=units)


class ParmedToFFComponent(TransComponent):
    """ A component for converting ParmEd molecule to Molecule object. """

    @classmethod
    def input(cls):
        return TransInput

    @classmethod
    def output(cls):
        return TransOutput

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

        ff = inputs.tk_object
        mm_units = forcefield.ForceField.get_units()

        # Need to map these to potential types
        bonds_units = forcefield.bonded.bonds.potentials.harmonic.Harmonic.get_units()
        angles_units = forcefield.bonded.angles.potentials.harmonic.Harmonic.get_units()
        dihedrals_units = forcefield.bonded.dihedrals.potentials.harmonic.Harmonic.get_units()
        lj_units = forcefield.nonbonded.potentials.lenjones.LennardJones.get_units()

        atom = ff.atoms[0]
        bond = ff.bonds[0] if len(ff.bonds) else None
        angle = ff.angles[0] if len(ff.angles) else None
        dihedral = ff.dihedrals[0] if len(ff.dihedrals) else None

        data = [(atom.name, atom.type, atom.charge, atom.mass) for atom in ff.atoms]
        names, symbols, charges, masses = zip(*data)
        charges = convert(
            charges, atom.ucharge.unit.get_symbol(), mm_units["charges_units"]
        )
        masses = convert(
            masses, atom.umass.unit.get_symbol(), mm_units["masses_units"]
        )

        rmin_factor = convert(1.0, atom.urmin.unit.get_name(), bonds_units["lengths_units"])
        rmin_14_factor = convert(
            1.0, atom.urmin_14.unit.get_name(), lj_units["sigma_units"]
        )
        uepsilon_factor = convert(
            1.0, atom.uepsilon.unit.get_name(), lj_units["epsilon_units"]
        )
        epsilon_14_factor = convert(
            1.0, atom.uepsilon_14.unit.get_name(), lj_units["epsilon_units"]
        )
        scaling_factor = 1.0 / 2 ** (1.0 / 6.0)  # rmin = 2^(1/6) sigma
        params = [
            (
                atom.rmin * rmin_factor * scaling_factor,
                atom.rmin_14 * rmin_14_factor * scaling_factor,
                atom.epsilon * uepsilon_factor,
                atom.epsilon_14 * epsilon_14_factor,
            )
            for atom in ff.atoms
        ]
        sigma, sigma_14, epsilon, epsilon_14 = zip(*params)

        lj = forcefield.nonbonded.potentials.LennardJones(sigma=sigma, epsilon=epsilon)
        # Need to include sigma_14 and epsilon_14
        nonbonded = forcefield.nonbonded.NonBonded(params=lj)

        if bond:
            bond_req_factor = convert(
                1.0, bond.type.ureq.unit.get_name(), bonds_units["lengths_units"]
            )
            bond_k_factor = convert(
                1.0, bond.type.uk.unit.get_name(), bonds_units["spring_units"]
            )
            bonds_lengths = [bond.type.req * bond_req_factor for bond in ff.bonds]
            bonds_k = [bond.type.k * bond_k_factor for bond in ff.bonds]
            bonds_type = [bond.funct for bond in ff.bonds]

            unique_bonds_type = set(bonds_type)
            if len(unique_bonds_type) > 1:
                params = [
                    bondTypes.get(btype)(spring=bonds_k[bonds_type == btype])
                    for btype in unique_bonds_type if bondTypes.get(btype)
                ]
            else:
                params = forcefield.bonded.bonds.potentials.Harmonic(spring=bonds_k, lengths=bonds_lengths)

            bonds = forcefield.bonded.Bonds(params=params)
        else:
            bonds = None

        if angle:
            angle_theta_factor = convert(
                1.0, angle.type.utheteq.unit.get_name(), angles_units["angles_units"]
            )
            angle_k_factor = convert(
                1.0, angle.type.uk.unit.get_name(), angles_units["spring_units"]
            )
            angles_ = [angle.type.theteq for angle in ff.angles]
            angles_k = [angle.type.k for angle in ff.angles]
            angles_type = [angle.funct for angle in ff.angles]

            unique_angles_type = set(angles_type)
            if len(unique_angles_type) > 1:
                params = [
                    angleTypes.get(btype)(spring=angles_k[angles_type == btype])
                    for btype in unique_angles_type if angleTypes.get(btype)
                ]
            else:
                params = forcefield.bonded.angles.potentials.Harmonic(spring=angles_k, angles=angles_)

            angles = forcefield.bonded.Angles(params=params)
        else:
            angles = None

        if False:  # dihedral:
            # Why the hell does parmed interpret dihedral.type as a list sometimes? See dialanine.top, last 8 dihedral entries. WEIRD!
            # This is a hack around this
            dtype = (
                dihedral.type[0] if isinstance(dihedral.type, list) else dihedral.type
            )

            dihedral_phi_factor = convert(
                1.0, dtype.uphi_k.unit.get_name(), dihedrals_units["angles_units"]
            )
            dihedral_phase_factor = convert(
                1.0, dtype.uphase.unit.get_name(), dihedrals_units["phases_units"]
            )

            # Need to look into per, scee, and scnb ... their meaning, units, etc.
            dihedrals = [
                (
                    dihedral.type[0].phase * dihedral_phase_factor,
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

        charge_groups = None
        exclusions = ff.nrexcl
        inclusions = None

        input_dict = {
            "masses": masses,
            "charges": charges,
            "bonds": bonds,
            "angles": angles,
            "dihedrals": None,
            "im_dihedrals": None,
            "nonbonded": nonbonded,
            "exclusions": exclusions,
            "inclusions": inclusions,
            "types": names,
            "symbols": symbols,
        }

        return True, TransOutput(schema_object=forcefield.ForceField(**input_dict))
