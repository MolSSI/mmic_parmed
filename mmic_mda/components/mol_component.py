from mmelemental.components.trans.template_component import TransComponent
from mmelemental.models.util.output import FileOutput
from mmelemental.models.molecule.mm_mol import Mol
from mmic_mda.models import MdaMol
from typing import Dict, Any, List, Tuple, Optional

try:
    from MDAnalysis import Universe
except:
    raise ModuleNotFoundError('Make sure MDAnalysis is installed.')

class MolToMdaComponent(TransComponent):
    """ A component for converting Molecule to MDAnalysis molecule object. """
    @classmethod
    def input(cls):
        return Mol

    @classmethod
    def output(cls):
        return MdaMol

    def execute(
        self,
        inputs: Mol,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, MdaMol]:

        natoms = len(inputs.masses)

        if inputs.residues:
            nres = len(inputs.residues)
            resnames, resids = zip(*inputs.residues)
        
        # Must account for segments as well
        segindices = None

        mda_mol = Universe.empty(natoms, n_residues=nres, atom_resindex=resids,
            residue_segindex=segindices, trajectory=True)

        mda_mol.add_TopologyAttr('type', inputs.symbols)
        mda_mol.add_TopologyAttr('mass', inputs.masses)

        if inputs.names:
            mda_mol.add_TopologyAttr('name', inputs.names)
        
        if inputs.residues:
            mda_mol.add_TopologyAttr('resname', resnames)
            mda_mol.add_TopologyAttr('resid', resids)

        #mda_mol.add_TopologyAttr('segid', ['SOL'])

        if inputs.geometry is not None:
            mda_mol.atoms.positions = inputs.geometry.reshape(natoms,3)

        if inputs.velocities is not None:
            mda_mol.atoms.velocities = inputs.velocities.reshape(natoms,3)

        if inputs.forces is not None:
            mda_mol.atoms.positions = inputs.forces.reshape(natoms,3)

        if inputs.connectivity:
            bonds = [(bond[0], bond[1]) for bond in inputs.connectivity]
            mda_mol.add_TopologyAttr('bonds', bonds)
            # How to load bond order?

        return True, MdaMol(mol=mda_mol)

class MdaToMolComponent(TransComponent):
    """ A component for converting MDAnalysis molecule to Molecule object. """

    @classmethod
    def input(cls):
        return MdaMol

    @classmethod
    def output(cls):
        return Mol

    def execute(
        self,
        inputs: MdaMol,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None) -> Tuple[bool, Mol]:
        
        assert inputs.dtype == 'mdanalysis'
        orient, validate, kwargs = False, None, None

        # get all properties + more from Universe?
        uni = inputs.mol
        geo = TransComponent.get(uni.atoms, 'positions')
        vel = TransComponent.get(uni.atoms, 'velocities')
        forces = TransComponent.get(uni.atoms, 'forces')
        symbs = TransComponent.get(uni.atoms, 'types')
        names = TransComponent.get(uni.atoms, 'names')
        masses = TransComponent.get(uni.atoms, 'masses')

        # If bond order is none, set it to 1.
        connectivity = [(bond.indices[0], bond.indices[1], bond.order or 1) for bond in uni.atoms.bonds]
        residues = [(atom.resname, atom.resnum) for atom in uni.atoms]

        input_dict = {'symbols': symbs, 
                      'geometry': geo, 
                      'velocities': vel,
                      'forces': forces,
                      'residues': residues, 
                      'connectivity': connectivity,
                      'masses': masses,
                      'names': names}

        if kwargs:
            input_dict.update(kwargs)

        return True, Mol(orient=orient, validate=validate, **input_dict)
