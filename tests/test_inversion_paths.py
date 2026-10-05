import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from sagnac_reference.inversion import velocity_from_delta_t, velocity_from_times


@pytest.mark.parametrize("beta", [0.0, 0.01, 0.1, 0.5, 0.9, -0.5, -0.9])
def test_inversion_round_trip_stable_root(beta):
    c, L = 1.0, 1.0
    v = beta * c
    dt = 2.0 * L * v / (c * c - v * v)
    v_rec = velocity_from_delta_t(dt, L, c)
    assert abs(v_rec - v) < 1e-12


def test_inversion_times_path():
    c = 1.0
    tp, tm = 1.25, 0.8333333333333334
    assert abs(velocity_from_times(tp, tm, c) - 0.2) < 1e-12
