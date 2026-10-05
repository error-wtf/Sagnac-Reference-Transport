import pytest
from sagnac_reference.analytic import SagnacConfig, directed_times, delta_t_axle, delta_tau_detector
from sagnac_reference.phase import detector_phase
from sagnac_reference.inversion import velocity_from_times, velocity_from_delta_t


def test_frame_and_phase_chain():
    cfg = SagnacConfig(2.0, 3.0, 0.8)
    assert delta_tau_detector(cfg) == pytest.approx(delta_t_axle(cfg) / cfg.gamma)
    omega = 19.0
    assert detector_phase(cfg, omega) == pytest.approx(omega * delta_tau_detector(cfg))


def test_inverse_closure():
    cfg = SagnacConfig(2.0, 3.0, -0.8)
    tp, tm = directed_times(cfg)
    assert velocity_from_times(tp, tm, cfg.c) == pytest.approx(cfg.v, abs=1e-13)
    assert velocity_from_delta_t(delta_t_axle(cfg), cfg.L, cfg.c) == pytest.approx(cfg.v, abs=1e-13)
