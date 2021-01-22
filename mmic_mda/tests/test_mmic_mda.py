"""
Unit and regression test for the mmic_mda package.
"""

# Import package, test suite, and other packages as needed
import mmic_mda
import pytest
import sys
import MDAnalysis as mda
import mmelemental as mm


def pytest_generate_tests(metafunc):
    if "guess_bonds" in metafunc.fixturenames:
        metafunc.parametrize("guess_bonds", ["False", "True"])


def test_mmic_mda_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmic_mda" in sys.modules


def test_mda_to_mol(guess_bonds):
    uni = mda.Universe("mmic_mda/data/1dzl_fixed.pdb", guess_bonds=guess_bonds)
    mda_mol = mmic_mda.models.MdaMol(mol=uni)
    mm_mol = mmic_mda.components.MdaToMolComponent.compute(mda_mol)

    return mm_mol


def test_mol_to_mda(guess_bonds):
    mm_mol = mm.models.molecule.mm_mol.Mol.from_file("mmic_mda/data/1dzl_fixed.pdb")

    return mmic_mda.components.MolToMdaComponent.compute(mm_mol)


def test_io_methods(guess_bonds):
    mda_mol = mmic_mda.models.MdaMol.from_file(
        "mmic_mda/data/1dzl_fixed.pdb", guess_bonds=guess_bonds
    )
    assert isinstance(mda_mol.mol, mda_mol.dtype)

    mm_mol = mda_mol.to_data()
    assert isinstance(mm_mol, mm.models.molecule.Mol)
