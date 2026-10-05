"""Dimensional scaling and unit-invariance contracts (Phase 6)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from sagnac_reference.analytic import (
    SagnacConfig,
    delta_tau_detector,
    directed_times,
    normalized_asymmetry,
)
from sagnac_reference.true_chain import return_time_true_chain


@pytest.mark.parametrize("lam", [10.0, 100.0])
def test_time_scales_linearly_with_L(lam):
    base = SagnacConfig(1.0, 1.0, 0.2)
    scaled = SagnacConfig(lam, 1.0, 0.2)
    tp0, tm0 = directed_times(base)
    tp1, tm1 = directed_times(scaled)
    assert tp1 == pytest.approx(lam * tp0, rel=1e-15)
    assert tm1 == pytest.approx(lam * tm0, rel=1e-15)
    chain = return_time_true_chain(scaled, +1, steps=32000).return_time
    assert abs(chain - tp1) < 1e-10 * lam


def test_dimensionless_invariance_under_unit_rescale():
    """Rescaling all speeds by k leaves beta, gamma and asymmetry invariant."""
    k = 100.0
    a = SagnacConfig(1.0, 1.0, 0.2)
    b = SagnacConfig(1.0, 1.0 * k, 0.2 * k)
    assert normalized_asymmetry(a) == pytest.approx(normalized_asymmetry(b), abs=1e-15)
    assert a.beta == pytest.approx(b.beta)
    assert a.gamma == pytest.approx(b.gamma)


def test_proper_time_scales_linearly_too():
    a = SagnacConfig(1.0, 1.0, 0.2)
    b = SagnacConfig(50.0, 1.0, 0.2)
    assert delta_tau_detector(b) == pytest.approx(50.0 * delta_tau_detector(a), rel=1e-15)
