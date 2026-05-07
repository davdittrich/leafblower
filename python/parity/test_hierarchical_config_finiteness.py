"""Tests: HierarchicalConfig rejects non-finite outer_tol (FIX-11 / leafblower-6ycz.1.21)."""
import pytest
from leafblower import HierarchicalConfig


@pytest.mark.parametrize("bad", [float("inf"), float("-inf"), float("nan")])
def test_outer_tol_non_finite_rejected(bad):
    with pytest.raises(ValueError, match="finite"):
        HierarchicalConfig(coarse_margins=["a"], outer_tol=bad)
