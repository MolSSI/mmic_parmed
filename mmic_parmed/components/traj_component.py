from mmelemental.components.trans.template_component import TransComponent
from mmelemental.models.util.output import FileOutput
from mmelemental.models.trajectory.mm_traj import Traj, Frame
from mmic_parmed.models import MdaMol, MdaTraj
from typing import Dict, Any, List, Tuple, Optional
from mmelemental.util.decorators import require

__all__ = ["MdaToTrajComponent"]


class MolToMdaComponent(TransComponent):
    """ A component for converting Molecule to MDAnalysis molecule object. """

    @classmethod
    def input(cls):
        return Mol

    @classmethod
    def output(cls):
        return MdaTraj

    @require("MDAnalysis")
    def execute(
        self,
        inputs: Traj,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, MdaTraj]:

        import MDAnalysis

        natoms = len(inputs.masses)

        if inputs.residues:
            nres = len(inputs.residues)
            resnames, resids = zip(*inputs.residues)

        # Must account for segments as well
        segindices = None

        mda_mol = Universe.empty(
            natoms,
            n_residues=nres,
            atom_resindex=resids,
            residue_segindex=segindices,
            trajectory=True,
        )

        mda_mol.add_TopologyAttr("type", inputs.symbols)
        mda_mol.add_TopologyAttr("mass", inputs.masses)

        if inputs.names is not None:
            mda_mol.add_TopologyAttr("name", inputs.names)

        if inputs.residues is not None:
            mda_mol.add_TopologyAttr("resname", resnames)
            mda_mol.add_TopologyAttr("resid", resids)

        # mda_mol.add_TopologyAttr('segid', ['SOL'])

        if inputs.geometry is not None:
            mda_mol.atoms.positions = inputs.geometry.reshape(natoms, 3)

        if inputs.velocities is not None:
            mda_mol.atoms.velocities = inputs.velocities.reshape(natoms, 3)

        if inputs.forces is not None:
            mda_mol.atoms.positions = inputs.forces.reshape(natoms, 3)

        if inputs.connectivity:
            bonds = [(bond[0], bond[1]) for bond in inputs.connectivity]
            mda_mol.add_TopologyAttr("bonds", bonds)
            # How to load bond order?

        return True, MdaMol(data=mda_mol)


class MdaToTrajComponent(TransComponent):
    """ A component for converting MDAnalysis molecule to Molecule object. """

    @classmethod
    def input(cls):
        return MdaTraj

    @classmethod
    def output(cls):
        return Traj

    def execute(
        self,
        inputs: MdaTraj,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, Traj]:

        mol = None
        orient, validate, kwargs = False, None, None
        uni = inputs.data

        if hasattr(uni.atoms, "names"):
            mda_mol = MdaMol(data=uni)
            mol = mda_mol.to_schema()

        frames = [
            Frame(
                geometry=frame.positions if frame.has_positions else None,
                geometry_units="A",
                velocities=frame.velocities if frame.has_velocities else None,
                velocities_units="A/ps",
                forces=frame.forces if frame.has_forces else None,
                forces_units="kJ/(mol*A)",
                timestep=frame.dt,
                timestep_units="ps",
            )
            for frame in uni.trajectory
        ]

        if kwargs:
            input_dict.update(kwargs)

        return True, Traj(top=mol, frames=frames)
