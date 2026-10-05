import pytest

from sagnac_reference.analytic import SagnacConfig, directed_times
from sagnac_reference.segment_chain import return_time_segment_chain


def test_chain_matches_oracle():
    cfg = SagnacConfig(1.0, 1.0, 0.23)
    tp, tm = directed_times(cfg)
    cp = return_time_segment_chain(cfg, +1, 5000).return_time
    cm = return_time_segment_chain(cfg, -1, 5000).return_time
    assert cp == pytest.approx(tp, abs=2e-12)
    assert cm == pytest.approx(tm, abs=2e-12)
