"""Phase 2+3: high-precision oracle tests + explicit series analysis.

Ground truth is verified against `decimal` (50-digit precision) so the float
implementation is anchored to an independent arithmetic, not to itself.
"""
from __future__ import annotations

import sys
from decimal import Decimal, getcontext
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from sagnac_reference.analytic import (
    SagnacConfig,
    delta_t_axle,
    delta_t_axle_closed,
    directed_times,
    normalized_asymmetry,
    series_partial,
    series_remainder_bound,
)

getcontext().prec = 50


def _decimal_dt(L: float, c: float, v: float) -> Decimal:
    """50-digit reference for dt = 2Lv/(c^2-v^2)."""
    Ld, cd, vd = (Decimal(x) for x in (L, c, v))
    return (2 * Ld * vd) / (cd * cd - vd * vd)


@pytest.mark.parametrize("beta", [1e-3, 0.1, 0.2, 0.5, 0.8, -0.5, -0.8])
def test_oracle_against_decimal_50_digit(beta):
    cfg = SagnacConfig(L=1.0, c=1.0, v=beta)
    ref = float(_decimal_dt(1.0, 1.0, beta))
    got = float(delta_t_axle_closed(cfg))
    # float has ~2.2e-16 relative eps; a correct formula must agree to that.
    assert abs(got - ref) <= 1e-15 * max(abs(ref), 1e-300)


def test_small_v_limit_matches_linear_term():
    """lim v->0: dt = 2Lv/c^2 (the linear term dominates)."""
    L, c = 1.0, 1.0
    for v in (1e-3, 1e-6, 1e-9):
        cfg = SagnacConfig(L, c, v)
        linear = 2.0 * L * v / (c * c)
        full = float(delta_t_axle_closed(cfg))
        # deviation from the linear term is the beta^2 correction ~ v^2 relative;
        # require it to be at most 1% for v <= 1e-3 (beta^2 <= 1e-6 << 1%).
        assert abs(full - linear) <= 0.01 * abs(linear)


def test_direction_swap_identity():
    """t_+(v) == t_-(-v) exactly (float-level)."""
    for v in (0.1, 0.37, 0.8):
        cfg_p = SagnacConfig(1.0, 1.0, v)
        cfg_m = SagnacConfig(1.0, 1.0, -v)
        tp, _ = directed_times(cfg_p)
        _, tm = directed_times(cfg_m)
        assert tp == tm


def test_dt_antisymmetry():
    for v in (0.1, 0.37, 0.8):
        cfg_p = SagnacConfig(1.0, 1.0, v)
        cfg_m = SagnacConfig(1.0, 1.0, -v)
        assert delta_t_axle(cfg_p) == -delta_t_axle(cfg_m)


def test_near_c_horizon_behaviour():
    """|v| -> c: t_+ diverges, dt -> finite-free divergence, no NaN/inf leak
    for the largest beta the domain guard allows."""
    cfg = SagnacConfig(1.0, 1.0, 1.0 - 1e-12)
    tp, tm = directed_times(cfg)
    assert tp > 1e12
    assert 0.0 < tm < 1.0


# ---------------- Phase 3: series with exact remainders ----------------

@pytest.mark.parametrize("beta", [0.1, 0.3, 0.5, 0.8])
@pytest.mark.parametrize("N", [0, 1, 5, 10, 30])
def test_series_partial_matches_exact_up_to_remainder(beta, N):
    cfg = SagnacConfig(1.0, 1.0, beta)
    for direction in (+1, -1):
        partial = series_partial(cfg, direction, N)
        tp, tm = directed_times(cfg)
        exact = tp if direction == +1 else tm
        bound = series_remainder_bound(cfg, direction, N)
        # bound is mathematically sharp; allow float rounding headroom.
        assert abs(partial - exact) <= bound * (1.0 + 1e-12) + 1e-15


def test_series_remainder_bound_is_tight():
    """The remainder bound must be sharp: actual error < bound, same order."""
    cfg = SagnacConfig(1.0, 1.0, 0.5)
    tp, tm = directed_times(cfg)
    for N in (5, 10, 20):
        err_p = abs(series_partial(cfg, +1, N) - tp)
        bound_p = series_remainder_bound(cfg, +1, N)
        assert err_p <= bound_p
        # geometric decay: err at N+5 is ~beta^5 of err at N
        err_p5 = abs(series_partial(cfg, +1, N + 5) - tp)
        assert err_p5 < err_p * cfg.beta**5 * 1.01


def test_normalized_asymmetry_identity():
    """(t+ - t-)/(t+ + t-) == v/c exactly (the paper's key identity)."""
    for v in (0.05, 0.2, 0.5, 0.9):
        cfg = SagnacConfig(1.0, 1.0, v)
        assert normalized_asymmetry(cfg) == pytest.approx(cfg.beta, abs=1e-15)


def test_even_odd_decomposition():
    """Sum of series isolates even terms; difference isolates odd terms."""
    cfg = SagnacConfig(1.0, 1.0, 0.4)
    tp, tm = directed_times(cfg)
    s, d = tp + tm, tp - tm
    # even part: 2*T0*(1 + b^2 + b^4 + ...) = 2*T0/(1-b^2)
    even_expected = 2.0 * cfg.T0 / (1.0 - cfg.beta**2)
    assert abs(s - even_expected) < 1e-15
    # odd part: 2*T0*b/(1-b^2)
    odd_expected = 2.0 * cfg.T0 * cfg.beta / (1.0 - cfg.beta**2)
    assert abs(d - odd_expected) < 1e-15
