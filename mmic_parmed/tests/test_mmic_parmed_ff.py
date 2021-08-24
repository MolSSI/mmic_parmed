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
import mm_data
from pathlib import Path


top_files = lambda ext: [
    mm_data.ffs[f"1dzl_gro.{ext}"],
    mm_data.ffs[f"dialanine.{ext}"],
]
json_files = [mm_data.ffs["water-ff.json"]]  # "alanine.json"]]


def pytest_generate_tests(metafunc):
    if "top_file" in metafunc.fixturenames:
        metafunc.parametrize("top_file", top_files("top"))

    if "json_file" in metafunc.fixturenames:
        metafunc.parametrize("json_file", json_files)


def test_mmic_parmed_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmic_parmed" in sys.modules


def test_parmed_to_ff(top_file, **kwargs):
    ff = parmed.load_file(top_file)
    inputs = {
        "data_object": ff,
        "keywords": kwargs,
        "schema_name": "mmel_input",
        "schema_version": 1,
    }
    return mmic_parmed.components.ParmedToFFComponent.compute(inputs)


def test_ff_to_parmed(json_file, **kwargs):
    mm_ff = mm.models.forcefield.mm_ff.ForceField.from_file(json_file)
    inputs = {
        "schema_object": mm_ff,
        "keywords": kwargs,
        "schema_name": "mmel_input",
        "schema_version": 1,
    }
    return mmic_parmed.components.FFToParmedComponent.compute(inputs)


def test_io_methods(top_file):

    pff = mmic_parmed.models.ParmedFF.from_file(top_file)
    assert isinstance(pff.data, pff.dtype())

    top_filename = Path("tmp.top")

    pff.to_file(top_filename.name)
    assert top_filename.is_file()
    # assert filecmp.cmp(top_filename, top_file, shallow=False), f"ff.top != {top_file}"
    top_filename.unlink()

    psf_filename = Path("tmp.psf")
    pff.to_file(psf_filename.name)
    assert psf_filename.is_file()
    psf_filename.unlink()

    mff = pff.to_schema()
    assert isinstance(mff, mm.models.forcefield.ForceField)
