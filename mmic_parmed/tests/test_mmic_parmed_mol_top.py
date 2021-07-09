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

cfile = lambda ext: mm_data.mols[f"dialanine.{ext}"]
ffile = lambda ext: mm_data.ffs[f"dialanine.{ext}"]
cexts = ["pdb", "gro", "inpcrd"]
fexts = ["top", "top", "prmtop"]


def pytest_generate_tests(metafunc):
    if "cfile" in metafunc.fixturenames:
        metafunc.parametrize("cfile", [cfile(ext) for ext in cexts])
    if "ffile" in metafunc.fixturenames:
        metafunc.parametrize("ffile", [ffile(ext) for ext in fexts])


def test_mmic_parmed_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmic_parmed" in sys.modules


def test_parmed_to_mol(cfile, ffile, **kwargs):
    struct = parmed.load_file(filename=ffile, xyz=cfile)
    inputs = {"data_object": struct, "keywords": kwargs}
    return mmic_parmed.components.ParmedToMolComponent.compute(inputs)


def test_mol_to_parmed(cfile, ffile):
    mmol = mm.models.molecule.mm_mol.Molecule.from_file(cfile, ffile)
    inputs = {"schema_object": mmol}
    return mmic_parmed.components.MolToParmedComponent.compute(inputs)


def test_io_methods(cfile, ffile):
    pmol = mmic_parmed.models.ParmedMol.from_file(cfile, ffile)
    assert isinstance(pmol.data, pmol.dtype)

    mmol = pmol.to_schema()
    assert isinstance(mmol, mm.models.molecule.Molecule)
