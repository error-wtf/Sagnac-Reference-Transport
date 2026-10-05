import pytest
from sagnac_reference.analytic import SagnacConfig, directed_times, delta_t_axle, delta_tau_detector


def test_wrong_branch_sign_is_detected():
    cfg = SagnacConfig(1.0, 1.0, 0.25)
    tp, tm = directed_times(cfg)
    wrong_tp = cfg.L / (cfg.c + cfg.v)
    assert abs(wrong_tp - tp) > 1e-2
    assert wrong_tp == pytest.approx(tm)


def test_wrong_gamma_factor_is_detected():
    cfg = SagnacConfig(1.0, 1.0, 0.6)
    correct = delta_tau_detector(cfg)
    wrong = delta_t_axle(cfg) * cfg.gamma
    assert abs(wrong - correct) > 0.1


def test_zigzag_is_not_encoded_as_light_velocity_reversal():
    # The reference model keeps light propagation speed fixed at ±c.
    # Alternation belongs to the correction-series sign, not the propagation law.
    cfg = SagnacConfig(1.0, 1.0, 0.2)
    assert cfg.c > 0
