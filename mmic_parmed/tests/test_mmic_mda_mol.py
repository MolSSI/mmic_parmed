"""
Unit and regression test for the mmic_parmed package.
"""

# Import package, test suite, and other packages as needed
import mmic_parmed
import pytest
import sys
import os
import MDAnalysis as mda
import mmelemental as mm


data_dir = os.path.join("mmic_parmed", "data")
top_file = os.path.join(data_dir, "molecules", "1dzl_fixed.pdb")


def pytest_generate_tests(metafunc):
    if "guess_bonds" in metafunc.fixturenames:
        metafunc.parametrize("guess_bonds", ["False", "True"])


def test_mmic_parmed_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmic_parmed" in sys.modules


def test_mda_to_mol(guess_bonds):
    uni = mda.Universe(top_file, guess_bonds=guess_bonds)
    mda_mol = mmic_parmed.models.MdaMol(data=uni)
    mm_mol = mmic_parmed.components.MdaToMolComponent.compute(mda_mol)

    return mm_mol


def test_mol_to_mda(guess_bonds):
    mm_mol = mm.models.molecule.mm_mol.Mol.from_file(top_file)

    return mmic_parmed.components.MolToMdaComponent.compute(mm_mol)


def test_io_methods(guess_bonds):
    mda_mol = mmic_parmed.models.MdaMol.from_file(top_file, guess_bonds=guess_bonds)
    assert isinstance(mda_mol.data, mda_mol.dtype)

    mm_mol = mda_mol.to_schema()
    assert isinstance(mm_mol, mm.models.molecule.Mol)
