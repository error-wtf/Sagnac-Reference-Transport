"""Gate tolerance dosage: the SAG-S9/S10 gates must PASS at the frozen
resolution and FAIL when the grid is deliberately starved.  A gate that
passes at every resolution measures nothing."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from sagnac_reference.analytic import SagnacConfig
from sagnac_reference.diagnostics import evaluate


def test_gate_passes_at_frozen_resolution():
    cfg = SagnacConfig(1.0, 1.0, 0.2)
    verdict = evaluate(cfg, pde_nx=2400)
    assert verdict["status"] == "SAGNAC_REFERENCE_CLOSURE_PASS"


def test_gate_fails_at_starved_resolution():
    """nx=200 (the minimum the domain guard allows) starves the PDE routes:
    SAG-S9 must FAIL — proves the tolerances are dosage-sensitive, not
    rubber stamps.  (nx<200 is correctly rejected by the guards.)"""
    cfg = SagnacConfig(1.0, 1.0, 0.2)
    verdict = evaluate(cfg, pde_nx=200)
    s9 = verdict["gates"]["SAG-S9-pde-independent"]
    assert s9 is False, "gate passed even with a starved grid — tolerance is decorative"
    assert verdict["status"] == "SAGNAC_REFERENCE_CLOSURE_FAIL"


def test_domain_guard_rejects_sub_minimum_grid():
    """The fail-closed grid guard itself: nx < 200 must raise."""
    import pytest
    cfg = SagnacConfig(1.0, 1.0, 0.2)
    with pytest.raises(ValueError, match="nx"):
        evaluate(cfg, pde_nx=60)
