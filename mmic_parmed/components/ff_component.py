from mmelemental.models import forcefield
from mmic_translator.models.io import (
    TransInput,
    TransOutput,
)
from mmic_translator.components import TransComponent
from typing import List, Tuple, Optional
from mmelemental.util.units import convert
import parmed

__all__ = ["FFToParmedComponent", "ParmedToFFComponent"]

bondTypes = {
    1: forcefield.bonded.bonds.potentials.Harmonic,
    2: forcefield.bonded.bonds.potentials.Gromos96,
}

angleTypes = {
    1: forcefield.bonded.angles.potentials.Harmonic,
}

# Need to fix: fudgeLJ, fudgeQQ
# Need to handle dihedrals properly


class FFToParmedComponent(TransComponent):
    """ A component for converting Molecule to ParmEd molecule object. """

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
        Too many for loops. Can easily reduce these esp in the ParmedToFFComponent.
        """
        if isinstance(inputs, dict):
            inputs = self.input()(**inputs)

        empty_atom = parmed.topologyobjects.Atom()
        mmff = inputs.schema_object
        pff = parmed.structure.Structure()

        masses = convert(
            mmff.masses, mmff.masses_units, empty_atom.umass.unit.get_symbol()
        )

        charges = TransComponent.get(mmff, "charges")
        charges = convert(
            charges, mmff.charges_units, empty_atom.ucharge.unit.get_symbol()
        )

        atomic_numbers = TransComponent.get(mmff, "atomic_numbers")
        atom_types = TransComponent.get(mmff, "defs")

        rmin, epsilon = self.get_nonbonded(mmff, empty_atom)

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
                # bonds=..., faster than connecting atoms one by one as done below?
                # angles=...,
                # dihedrals=...,
                # impropers=...,
                # polarizable=...,
            )

            residues = TransComponent.get(mmff, "substructs")

            if residues:
                resname, resnum = residues[index]
            else:
                raise NotImplementedError(
                    "Residues must be supplied for forcefields based on atom typing."
                )

            pff.add_atom(atom, resname, resnum, chain="", inscode="", segid="")

        # Bonds
        bonds = TransComponent.get(mmff, "bonds")
        if bonds is not None:
            assert (
                mmff.bonds.form == "Harmonic"
            ), "Only Harmonic potential supported for now"

            spring = convert(
                bonds.params.spring, bonds.params.spring_units, "kcal/mol/angstroms**2"
            )
            req = convert(bonds.lengths, bonds.lengths_units, "angstroms")

            for (
                bi,
                (
                    i,
                    j,
                    order,
                ),
            ) in enumerate(mmff.bonds.indices):
                btype = parmed.topologyobjects.BondType(
                    k=spring[bi], req=req[bi], list=pff.bond_types
                )
                pff.bonds.append(
                    parmed.topologyobjects.Bond(
                        pff.atoms[i], pff.atoms[j], order=order, type=btype
                    )
                )
                pff.bond_types.append(btype)
                # both implementations seem to perform almost the same:
                # pff.atoms[i].bond_to(pff.atoms[j])

        # Angles
        angles = TransComponent.get(mmff, "angles")
        if angles is not None:
            assert (
                mmff.angles.form == "Harmonic"
            ), "Only Harmonic potential supported for now"

            spring = convert(
                angles.params.spring, angles.params.spring_units, "kcal/mol/radians^2"
            )
            angles_eq = convert(angles.angles, angles.angles_units, "degrees")

            for ai, (i, j, k) in enumerate(mmff.angles.indices):
                atype = parmed.topologyobjects.AngleType(
                    k=spring[ai], theteq=angles_eq[ai], list=pff.angle_types
                )
                pff.angles.append(
                    parmed.topologyobjects.Angle(
                        pff.atoms[i], pff.atoms[j], pff.atoms[k], type=atype
                    )
                )
                pff.angle_types.append(atype)

        # Dihedrals
        dihedrals = TransComponent.get(mmff, "dihedrals")
        if dihedrals is not None:
            assert (
                mmff.dihedrals.form == "Charmm"
            ), "Only Charmm-style potential supported for now"

            energy = convert(
                dihedrals.params.energy, dihedrals.params.energy_units, "kcal/mol"
            )
            phase = convert(
                dihedrals.params.phase, dihedrals.params.phase_units, "degrees"
            )
            periodicity = dihedrals.params.periodicity

            for di, (i, j, k, l) in enumerate(mmff.dihedrals.indices):
                dtype = parmed.topologyobjects.DihedralType(
                    phi_k=energy[di],
                    per=periodicity[di],
                    phase=phase[di],
                    # scee
                    # scnb
                    list=pff.dihedral_types,
                )
                # assert:
                # dtype.funct = (
                #    9  # hackish: assume all dihedrals are proper and charmm-style
                # )
                pff.dihedrals.append(
                    parmed.topologyobjects.Dihedral(
                        pff.atoms[i],
                        pff.atoms[j],
                        pff.atoms[k],
                        pff.atoms[l],
                        improper=False,
                        type=dtype,
                    )
                )
                pff.dihedral_types.append(dtype)

        return True, TransOutput(proc_input=inputs, data_object=pff)

    def get_nonbonded(
        self,
        mmff: forcefield.ForceField,
        empty_atom: parmed.topologyobjects.Atom,
    ) -> Tuple["numpy.ndarray", "numpy.ndarray"]:

        assert (
            mmff.nonbonded.form == "LennardJones"
        ), "Only LJ potential supported for now"

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
        return rmin, epsilon


class ParmedToFFComponent(TransComponent):
    """ A component for converting ParmEd molecule to Molecule object. """

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

        ff = inputs.data_object
        mm_units = forcefield.ForceField.get_units()

        # Need to map these to potential types
        lj_units = forcefield.nonbonded.potentials.lenjones.LennardJones.get_units()

        atom = ff.atoms[0]
        bond = ff.bonds[0] if len(ff.bonds) else None
        angle = ff.angles[0] if len(ff.angles) else None
        dihedral = ff.dihedrals[0] if len(ff.dihedrals) else None

        data = [
            (
                atom.name,
                atom.type,
                atom.charge,
                atom.mass,
                atom.atomic_number,
                atom.element_name,
            )
            for atom in ff.atoms
        ]
        names, types, charges, masses, atomic_numbers, elements = zip(*data)
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
            connectivity = [
                (bond.atom1.idx, bond.atom2.idx, bond.order or 1) for bond in ff.bonds
            ]

            unique_bonds_type = set(bonds_type)

            if len(unique_bonds_type) > 1:
                raise NotImplementedError("Multiple bond types not yet supported.")
                # params = [
                #    bondTypes.get(btype)(spring=bonds_k[bonds_type == btype])
                #    for btype in unique_bonds_type
                #    if bondTypes.get(btype)
                # ]
            else:
                params = forcefield.bonded.bonds.potentials.Harmonic(spring=bonds_k)

            bonds = forcefield.bonded.Bonds(
                params=params,
                lengths=bonds_lengths,
                indices=connectivity,
                form="Harmonic",
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
            angles_ = [angle.type.theteq * angle_theta_factor for angle in ff.angles]
            angles_k = [angle.type.k * angle_k_factor for angle in ff.angles]
            angles_type = [angle.funct for angle in ff.angles]

            angles_indices = [
                (angle.atom1.idx, angle.atom2.idx, angle.atom3.idx)
                for angle in ff.angles
            ]

            unique_angles_type = set(angles_type)
            if len(unique_angles_type) > 1:
                raise NotImplementedError("Multiple angle types not yet supported.")
                # params = [
                #    angleTypes.get(btype)(spring=angles_k[angles_type == btype])
                #    for btype in unique_angles_type
                #    if angleTypes.get(btype)
                # ]
            else:
                params = forcefield.bonded.angles.potentials.Harmonic(spring=angles_k)

            angles = forcefield.bonded.Angles(
                params=params, angles=angles_, indices=angles_indices, form="Harmonic"
            )
        else:
            angles = None

        if dihedral:
            dihedrals_units = forcefield.bonded.dihedrals.potentials.Charmm.get_units()
            dihedrals_units.update(forcefield.bonded.Dihedrals.get_units())

            # Why the hell does parmed interpret dihedral.type as a list sometimes? See dialanine.top, last 8 dihedral entries. WEIRD!
            # This is a hack around this
            # dihedrals_type = (
            #    dihedral.type[0] if isinstance(dihedral.type, list) else dihedral.type
            # )

            dihedrals_type = [di.funct for di in ff.dihedrals if di.funct != 4]
            # funct == 4: periodic improper dihedralm see https://manual.gromacs.org/documentation/2019/reference-manual/topologies/topology-file-formats.html#tab-topfile2

            dihedral_phi_factor = convert(
                1.0,
                "kcal/mol",  # dihedrals_type.uphi_k.unit.get_name(),
                dihedrals_units["energy_units"],
            )
            dihedral_phase_factor = convert(
                1.0,
                "degrees",  # dihedrals_type.uphase.unit.get_name(),
                dihedrals_units["phase_units"],
            )

            dihedrals_indices = [
                (
                    dihedral.atom1.idx,
                    dihedral.atom2.idx,
                    dihedral.atom3.idx,
                    dihedral.atom4.idx,
                )
                for dihedral in ff.dihedrals
            ]

            # Need to look into per, scee, and scnb ... their meaning, units, etc.
            dihedrals_params = [
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

            phase, energy, per, _, _ = zip(*dihedrals_params)

            unique_dihedrals_type = set(dihedrals_type)
            if len(unique_dihedrals_type) > 1:
                raise NotImplementedError(
                    "Multiple dihedral types not yet supported.", unique_dihedrals_type
                )
            else:
                params = forcefield.bonded.dihedrals.potentials.Charmm(
                    phase=phase,
                    energy=energy,
                    periodicity=per,
                )

            dihedrals = forcefield.bonded.Dihedrals(
                params=params,
                angles=angles_,
                indices=dihedrals_indices,
                form="Charmm",
            )
        else:
            dihedrals = None

        if hasattr(ff, "residues"):
            residues = [(atom.residue.name, atom.residue.idx) for atom in ff.atoms]

        # charge_groups = None ... should support charge_groups?
        exclusions = ff.nrexcl
        inclusions = None

        input_dict = {
            "masses": masses,
            "charges": charges,
            "bonds": bonds,
            "angles": angles,
            "dihedrals": dihedrals,
            "im_dihedrals": None,
            "nonbonded": nonbonded,
            "exclusions": exclusions,
            "inclusions": inclusions,
            "defs": types,  # or names?
            "symbols": elements,
            "substructs": residues,
            "atomic_numbers": atomic_numbers,
        }

        ff = forcefield.ForceField(**input_dict)
        return True, TransOutput(proc_input=inputs, schema_object=ff)
