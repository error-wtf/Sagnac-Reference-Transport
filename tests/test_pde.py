import pytest

from sagnac_reference.analytic import SagnacConfig, directed_times
from sagnac_reference.transport_pde import return_time_pde


def test_pde_matches_oracle_moderate_beta():
    cfg = SagnacConfig(1.0, 1.0, 0.2)
    tp, tm = directed_times(cfg)
    pp = return_time_pde(cfg, +1, nx=4800)
    pm = return_time_pde(cfg, -1, nx=4800)
    assert pp.return_time == pytest.approx(tp, abs=3.5e-3)
    assert pm.return_time == pytest.approx(tm, abs=3.5e-3)
