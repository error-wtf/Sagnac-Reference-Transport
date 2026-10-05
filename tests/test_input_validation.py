"""Input-validation guards: untested ValueError paths on public entry points."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from sagnac_reference.analytic import (
    SagnacConfig,
    directed_times,
    series_partial,
)
from sagnac_reference.segment_chain import (
    return_time_segment_chain,
)


def test_config_rejects_bad_geometry():
    with pytest.raises(ValueError, match="L"):
        SagnacConfig(L=0.0, c=1.0, v=0.1)
    with pytest.raises(ValueError, match="c"):
        SagnacConfig(L=1.0, c=-1.0, v=0.1)
    with pytest.raises(ValueError, match=r"\|v\| < c"):
        SagnacConfig(L=1.0, c=1.0, v=1.0)


def test_series_partial_rejects_bad_direction_and_n():
    cfg = SagnacConfig(1.0, 1.0, 0.2)
    with pytest.raises(ValueError, match="direction"):
        series_partial(cfg, 0, 5)
    with pytest.raises(ValueError, match="N"):
        series_partial(cfg, +1, -1)


def test_segment_chain_rejects_bad_direction():
    cfg = SagnacConfig(1.0, 1.0, 0.2)
    with pytest.raises(ValueError):
        return_time_segment_chain(cfg, 7)


def test_directed_times_domains():
    cfg = SagnacConfig(1.0, 1.0, 0.2)
    tp, tm = directed_times(cfg)
    assert tp > tm > 0  # co-rotating takes longer
