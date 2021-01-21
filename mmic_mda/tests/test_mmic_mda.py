"""
Unit and regression test for the mmic_mda package.
"""

# Import package, test suite, and other packages as needed
import mmic_mda
import pytest
import sys
import MDAnalysis as mda
import mmelemental as mm

def test_mmic_mda_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmic_mda" in sys.modules

def test_mda_to_mol(guess_bonds):
    uni = mda.Universe('mmic_mda/data/1dzl_fixed.pdb', guess_bonds=guess_bonds)
    mda_mol = mmic_mda.models.MdaMolecule(mol=uni)
    mm_mol = mmic_mda.components.MdaToMolComponent.compute(mda_mol)
    mm_mol.to_file('1dzl_mm.pdb')

    return mm_mol

def test_mol_to_mda():
    with mm.models.util.output.FileOutput(path='1dzl_mm.pdb', clean=True) as fo:
        mm_mol = mm.models.molecule.mm_molecule.Molecule.from_file(fo.abs_path)

    return mmic_mda.components.MolToMdaComponent.compute(mm_mol)

mm_mol = test_mda_to_mol(guess_bonds=False)
mda_mol = test_mol_to_mda()

print(mm_mol)
print(mda_mol)
print(mm_mol.geometry)

mm_mol = test_mda_to_mol(guess_bonds=True)
mda_mol = test_mol_to_mda()

print(mm_mol)
print(mda_mol)
