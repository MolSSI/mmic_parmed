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
top_files = lambda ext: os.path.join(data_dir, "forcefields", f"1dzl_gro.{ext}")
json_files = os.path.join(data_dir, "forcefields", "forcefield-single.json")


def pytest_generate_tests(metafunc):
    if "top_file" in metafunc.fixturenames:
        metafunc.parametrize("top_file", [top_files("top")])

    if "json_file" in metafunc.fixturenames:
        metafunc.parametrize("json_file", [json_files])


def test_mmic_parmed_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmic_parmed" in sys.modules


def test_parmed_to_ff(top_file, **kwargs):
    ff = parmed.load_file(top_file)
    inputs = {"tk_object": ff, "kwargs": kwargs}
    return mmic_parmed.components.ParmedToFFComponent.compute(inputs)


def test_ff_to_parmed(json_file, **kwargs):
    mm_ff = mm.models.forcefield.mm_ff.ForceField.from_file(json_file)
    inputs = {"schema_object": mm_ff, "kwargs": kwargs}
    return mmic_parmed.components.FFToParmedComponent.compute(inputs)


def test_io_methods(top_file):
    pff = mmic_parmed.models.ParmedFF.from_file(top_file)
    assert isinstance(pff.data, pff.dtype)

    pff.to_file("ff.psf")
    os.remove("ff.psf")

    mff = pff.to_schema()
    assert isinstance(mff, mm.models.forcefield.ForceField)
