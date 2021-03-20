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

angleTypes = {
    1: forcefield.bonded.angles.potentials.Harmonic,
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

        empty_atom = parmed.topologyobjects.Atom()
        mmff = inputs.schema_object
        pff = parmed.structure.Structure()
        natoms = len(mmff.symbols)

        masses = TransComponent.get(mmff, "masses")
        masses = convert(
            mmff.masses, mmff.masses_units, empty_atom.umass.unit.get_symbol()
        )

        charges = TransComponent.get(mmff, "charges")
        charges = convert(
            charges, mmff.charges_units, empty_atom.ucharge.unit.get_symbol()
        )

        atomic_numbers = TransComponent.get(mmff, "atomic_numbers")
        atom_types = TransComponent.get(mmff, "types")

        # Non-bonded
        assert (
            mmff.nonbonded.form == "LennardJones",
            "Only LJ potential supported for now",
        )

        lj_units = forcefield.nonbonded.potentials.lenjones.LennardJones.get_units()
        scaling_factor = 2 ** (1.0 / 6.0)  # rmin = 2^(1/6) sigma

        rmin = mmff.nonbonded.params.sigma * scaling_factor
        rmin = convert(rmin, lj_units["sigma_units"], empty_atom.urmin.unit.get_name())
        # atom.rmin_14 * rmin_14_factor * scaling_factor,
        epsilon = convert(
            mmff.nonbonded.params.epsilon,
            lj_units["epsilon_units"],
            empty_atom.uepsilon.unit.get_name(),
        )
        # atom.epsilon_14 * epsilon_14_factor,

        # Bonds
        if False:  # mmff.bonds is not None:
            assert (
                mmff.bonds.form == "Harmonic",
                "Only Harmonic potential supported for now",
            )

            for (
                i,
                j,
                _,
            ) in mmff.bonds.indices:  # _: bond order ... can we use it in ParmEd?
                # pmol.atoms[i].bond_to(pmol.atoms[j])
                pmol.bonds.append(
                    parmed.topologyobjects.Bond(pmol.atoms[i], pmol.atoms[j])
                )
                # both implementations seem to perform almost the same

        # Angles
        if False:  # mmff.angles is not None:
            assert (
                mmff.angles.form == "Harmonic",
                "Only Harmonic potential supported for now",
            )

            for i, j, k in mmff.angles.indices:
                pmol.angles.append(
                    parmed.topologyobjects.Angle(
                        pmol.atoms[i], pmol.atoms[j], pmol.atoms[k]
                    )
                )

        # Dihedrals
        if False:  # mmff.dihedrals is not None:
            assert (
                mmff.dihedrals.form == "Harmonic",
                "Only Harmonic potential supported for now",
            )
            for i, j, k, l in mmff.dihedrals.indices:
                pmol.dihedrals.append(
                    parmed.topologyobjects.Dihedral(
                        pmol.atoms[i], pmol.atoms[j], pmol.atoms[k], pmol.atoms[l]
                    )
                )

        for index, symb in enumerate(mmff.symbols):

            # Will likely lose FF-related info ... but then Molecule is not supposed to store any params specific to FFs
            if atomic_numbers is not None:
                atomic_number = atomic_numbers[index]
            else:
                atomic_number = None

            if atom_types is not None:
                atom_type = atom_types[index]
            else:
                atom_type = None

            if masses is not None:
                mass = masses[index]
            else:
                mass = None

            if charges is not None:
                charge = charges[index]
            else:
                charge = None

            atom = parmed.topologyobjects.Atom(
                list=None,
                atomic_number=atomic_number,
                name=symb,
                type=atom_type,
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
                rmin=rmin[index],
                epsilon=epsilon[index],
                rmin14=None,
                epsilon14=None,
                # bonds=...,
                # angles=...,
                # dihedrals=...,
                # impropers
                # polarizable
            )

            pff.add_atom(atom, "", 0, chain="", inscode="", segid="")

        return True, TransOutput(tk_object=pff)


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
        masses = convert(masses, atom.umass.unit.get_symbol(), mm_units["masses_units"])

        rmin_factor = convert(1.0, atom.urmin.unit.get_name(), lj_units["sigma_units"])
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
        nonbonded = forcefield.nonbonded.NonBonded(params=lj, form="LennardJones")

        if bond:
            bonds_units = (
                forcefield.bonded.bonds.potentials.harmonic.Harmonic.get_units()
            )
            bonds_units.update(forcefield.bonded.Bonds.get_units())

            bond_req_factor = convert(
                1.0, bond.type.ureq.unit.get_name(), bonds_units["lengths_units"]
            )
            bond_k_factor = convert(
                1.0, bond.type.uk.unit.get_name(), bonds_units["spring_units"]
            )
            bonds_lengths = [bond.type.req * bond_req_factor for bond in ff.bonds]

            bonds_type = [bond.funct for bond in ff.bonds]
            bonds_k = [bond.type.k * bond_k_factor for bond in ff.bonds]

            unique_bonds_type = set(bonds_type)

            if len(unique_bonds_type) > 1:
                raise NotImplementedError("Multiple bond types not yet supported.")
                params = [
                    bondTypes.get(btype)(spring=bonds_k[bonds_type == btype])
                    for btype in unique_bonds_type
                    if bondTypes.get(btype)
                ]
            else:
                params = forcefield.bonded.bonds.potentials.Harmonic(spring=bonds_k)

            bonds = forcefield.bonded.Bonds(
                params=params, lengths=bonds_lengths, form="Harmonic"
            )
        else:
            bonds = None

        if angle:
            angles_units = (
                forcefield.bonded.angles.potentials.harmonic.Harmonic.get_units()
            )
            angles_units.update(forcefield.bonded.Angles.get_units())

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
                raise NotImplementedError("Multiple angle types not yet supported.")
                params = [
                    angleTypes.get(btype)(spring=angles_k[angles_type == btype])
                    for btype in unique_angles_type
                    if angleTypes.get(btype)
                ]
            else:
                params = forcefield.bonded.angles.potentials.Harmonic(spring=angles_k)

            angles = forcefield.bonded.Angles(
                params=params, angles=angles_, form="Harmonic"
            )
        else:
            angles = None

        if False:  # dihedral:
            dihedrals_units = (
                forcefield.bonded.dihedrals.potentials.harmonic.Harmonic.get_units()
            )
            dihedrals_units.update(forcefield.bonded.Dihedrals.get_units())

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
