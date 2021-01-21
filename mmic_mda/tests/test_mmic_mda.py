"""
Unit and regression test for the mmic_mda package.
"""

# Import package, test suite, and other packages as needed
import mmic_mda
import pytest
import sys

def test_mmic_mda_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmic_mda" in sys.modules
