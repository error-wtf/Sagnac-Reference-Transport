import sys
from pathlib import Path
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from sagnac_reference.analytic import SagnacConfig, directed_times
from sagnac_reference.transport_pde import return_time_pde


def test_pde_convergence_order_about_two():
    """Frozen from the 2026-10-05 convergence campaign: with grid-independent
    pulse width the Lax-Wendroff route converges at ~2nd order. Gate: err at
    nx=2400 below 5e-3 (campaign max 4.225e-04, ~12x margin)."""
    cfg = SagnacConfig(1.0, 1.0, 0.2)
    tp, tm = directed_times(cfg)
    pp = return_time_pde(cfg, +1, nx=2400).return_time
    pm = return_time_pde(cfg, -1, nx=2400).return_time
    assert abs(pp - tp) < 5e-3
    assert abs(pm - tm) < 5e-3


def test_pde_pulse_width_guard_fail_closed():
    cfg = SagnacConfig(1.0, 1.0, 0.2)
    # explicit narrow pulse: sigma < 8dx must fail closed
    with pytest.raises(ValueError):
        return_time_pde(cfg, +1, nx=2400, pulse_sigma=1e-4)  # 1e-4 < 8/2400=3.3e-3
