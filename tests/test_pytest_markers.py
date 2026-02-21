from __future__ import annotations

import pytest


def test_expected_pytest_markers_are_registered(pytestconfig: pytest.Config) -> None:
    markers = pytestconfig.getini("markers")
    marker_names = {marker.split(":", maxsplit=1)[0].strip() for marker in markers}
    assert {"performance", "smoke", "integration", "unit"} <= marker_names
