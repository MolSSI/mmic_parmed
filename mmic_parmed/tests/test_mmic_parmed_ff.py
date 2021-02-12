"""
Forcefield tests for the mmic_parmed package.
"""

# Import package, test suite, and other packages as needed
import mmic_parmed
import pytest
import sys
import os
import parmed
import mmelemental as mm


data_dir = os.path.join("mmic_parmed", "data")
top_file = lambda ext: os.path.join(data_dir, "forcefields", f"1dzl_gro.{ext}")


def pytest_generate_tests(metafunc):
    if "file" in metafunc.fixturenames:
        metafunc.parametrize("file", [top_file("top")])


def test_mmic_parmed_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmic_parmed" in sys.modules


def test_parmed_to_ff(file):
    ff = parmed.load_file(file)
    pff = mmic_parmed.models.ParmedFF(data=ff)
    return mmic_parmed.components.ParmedToFFComponent.compute(pff)


def test_ff_to_parmed(file):
    mm_ff = mm.models.forcefield.mm_ff.ForceField.from_file(file)
    # return mmic_parmed.components.FFToParmedComponent.compute(mm_ff)


def test_io_methods(file):
    pff = mmic_parmed.models.ParmedFF.from_file(file)
    assert isinstance(pff.data, pff.dtype)

    pff.to_file("ff.psf")
    os.remove("ff.psf")

    mff = pff.to_schema()
    assert isinstance(mff, mm.models.forcefield.ForceField)
