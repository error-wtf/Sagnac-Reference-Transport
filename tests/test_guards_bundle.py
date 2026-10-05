"""Fail-closed guards across all modules (the untested ValueError paths)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from sagnac_reference.analytic import SagnacConfig, series_partial
from sagnac_reference.inversion import (
    velocity_from_delta_t,
    velocity_from_times,
)
from sagnac_reference.phase import detector_phase
from sagnac_reference.segment_chain import return_time_segment_chain
from sagnac_reference.transport_pde import return_time_pde
from sagnac_reference.transport_pde_upwind import return_time_pde_upwind
from sagnac_reference.true_chain import return_time_true_chain


def test_inversion_guards():
    with pytest.raises(ValueError):
        velocity_from_times(1.0, 1.0, -1.0)
    with pytest.raises(ValueError):
        velocity_from_times(0.0, 1.0, 1.0)
    with pytest.raises(ValueError):
        velocity_from_times(1.0, -1.0, 1.0)
    with pytest.raises(ValueError):
        velocity_from_delta_t(0.1, 0.0, 1.0)
    with pytest.raises(ValueError):
        velocity_from_delta_t(0.1, 1.0, -1.0)


def test_phase_guard():
    cfg = SagnacConfig(1.0, 1.0, 0.2)
    with pytest.raises(ValueError):
        detector_phase(cfg, float("inf"))


def test_series_direction_guard():
    cfg = SagnacConfig(1.0, 1.0, 0.2)
    with pytest.raises(ValueError):
        series_partial(cfg, 0, 5)
    with pytest.raises(ValueError):
        series_partial(cfg, +1, -1)


def test_segment_chain_direction_guard():
    cfg = SagnacConfig(1.0, 1.0, 0.2)
    with pytest.raises(ValueError):
        return_time_segment_chain(cfg, 7)


def test_pde_guards_fail_closed():
    cfg = SagnacConfig(1.0, 1.0, 0.2)
    for solver in (return_time_pde, return_time_pde_upwind):
        with pytest.raises(ValueError):  # direction
            solver(cfg, 0, nx=400)
        with pytest.raises(ValueError):  # nx too small
            solver(cfg, +1, nx=100)
        with pytest.raises(ValueError):  # sigma too narrow for grid
            solver(cfg, +1, nx=4000, pulse_sigma=1e-4)


def test_true_chain_direction_guard():
    cfg = SagnacConfig(1.0, 1.0, 0.2)
    with pytest.raises(ValueError):
        return_time_true_chain(cfg, 0, steps=1000)
