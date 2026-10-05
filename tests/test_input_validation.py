"""Input-validation guards on all public entry points."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from sagnac_reference.analytic import SagnacConfig
from sagnac_reference.inversion import (
    velocity_from_delta_t,
    velocity_from_times,
)
from sagnac_reference.phase import detector_phase


def test_velocity_from_times_rejects_bad_input():
    with pytest.raises(ValueError):
        velocity_from_times(1.0, 1.0, -1.0)   # c < 0
    with pytest.raises(ValueError):
        velocity_from_times(0.0, 1.0, 1.0)    # t_plus <= 0
    with pytest.raises(ValueError):
        velocity_from_times(1.0, -1.0, 1.0)   # t_minus <= 0


def test_velocity_from_delta_t_rejects_bad_geometry():
    with pytest.raises(ValueError):
        velocity_from_delta_t(0.1, 0.0, 1.0)  # L <= 0
    with pytest.raises(ValueError):
        velocity_from_delta_t(0.1, 1.0, -1.0) # c <= 0


def test_detector_phase_rejects_nonfinite_omega():
    cfg = SagnacConfig(1.0, 1.0, 0.2)
    with pytest.raises(ValueError):
        detector_phase(cfg, float("inf"))
    with pytest.raises(ValueError):
        detector_phase(cfg, float("nan"))


def test_detector_phase_positive_for_positive_omega():
    cfg = SagnacConfig(1.0, 1.0, 0.2)
    assert detector_phase(cfg, 10.0) > 0.0
