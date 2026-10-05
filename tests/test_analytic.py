import math
import pytest

from sagnac_reference.analytic import (
    SagnacConfig, directed_times, delta_t_axle, delta_t_axle_closed,
    normalized_asymmetry, series_partial,
)


def test_closed_forms_and_parity():
    cfg = SagnacConfig(3.0, 5.0, 0.7)
    tp, tm = directed_times(cfg)
    assert delta_t_axle(cfg) == pytest.approx(delta_t_axle_closed(cfg), rel=1e-14)
    assert normalized_asymmetry(cfg) == pytest.approx(cfg.beta, rel=1e-14)
    rev = SagnacConfig(cfg.L, cfg.c, -cfg.v)
    rtp, rtm = directed_times(rev)
    assert rtp == pytest.approx(tm)
    assert rtm == pytest.approx(tp)
    assert delta_t_axle(rev) == pytest.approx(-delta_t_axle(cfg))


def test_static_limit():
    cfg = SagnacConfig(2.0, 7.0, 0.0)
    tp, tm = directed_times(cfg)
    assert tp == tm == pytest.approx(cfg.L / cfg.c)
    assert delta_t_axle(cfg) == 0.0


def test_series_converges():
    cfg = SagnacConfig(1.0, 1.0, 0.2)
    tp, tm = directed_times(cfg)
    assert series_partial(cfg, +1, 40) == pytest.approx(tp, abs=1e-14)
    assert series_partial(cfg, -1, 40) == pytest.approx(tm, abs=1e-14)


def test_domain_guards():
    with pytest.raises(ValueError):
        SagnacConfig(1.0, 1.0, 1.0)
    with pytest.raises(ValueError):
        SagnacConfig(-1.0, 1.0, 0.1)
