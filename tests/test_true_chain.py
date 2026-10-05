import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from sagnac_reference.analytic import SagnacConfig, directed_times
from sagnac_reference.true_chain import return_time_true_chain


@pytest.mark.parametrize("beta", [0.0, 0.1, 0.2, 0.5, -0.2, -0.5])
def test_true_chain_matches_oracle(beta):
    cfg = SagnacConfig(1.0, 1.0, beta)
    tp, tm = directed_times(cfg)
    rp = return_time_true_chain(cfg, +1, steps=32000).return_time
    rm = return_time_true_chain(cfg, -1, steps=32000).return_time
    assert abs(rp - tp) < 1e-10
    assert abs(rm - tm) < 1e-10


def test_true_chain_direction_reversal():
    cfg = SagnacConfig(1.0, 1.0, 0.3)
    rp = return_time_true_chain(cfg, +1, steps=32000).return_time
    rm = return_time_true_chain(cfg, -1, steps=32000).return_time
    assert rp > rm  # co-rotating return is slower
