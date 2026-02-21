import iDiffIR
import pytest


@pytest.mark.unit
def test_version_present():
    assert hasattr(iDiffIR, "__version__")
