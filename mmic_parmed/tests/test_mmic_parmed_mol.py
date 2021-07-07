"""
Unit and regression test for the mmic_parmed package.
"""

# Import package, test suite, and other packages as needed
import mmic_parmed
import pytest
import sys
import os
import parmed
import mmelemental as mm
import mm_data

cfile = lambda ext: mm_data.mols[f"1dzl_fixed.{ext}"]
exts = ["pdb", "gro"]


def pytest_generate_tests(metafunc):
    if "file" in metafunc.fixturenames:
        metafunc.parametrize("file", [cfile(ext) for ext in exts])


def test_mmic_parmed_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmic_parmed" in sys.modules


def test_parmed_to_mol(file, **kwargs):
    struct = parmed.load_file(file)
    inputs = {"data_object": struct, "kwargs": kwargs}
    return mmic_parmed.components.ParmedToMolComponent.compute(inputs)


def test_mol_to_parmed(file):
    mmol = mm.models.molecule.mm_mol.Molecule.from_file(file)
    inputs = {"schema_object": mmol}
    return mmic_parmed.components.MolToParmedComponent.compute(inputs)


def test_io_methods(file):
    pmol = mmic_parmed.models.ParmedMol.from_file(file)
    assert isinstance(pmol.data, pmol.dtype)

    mmol = pmol.to_schema()
    assert isinstance(mmol, mm.models.molecule.Molecule)
