"""Cross-route beta-scan regression against the frozen convergence campaign.

The campaign artifacts (artifacts/convergence/convergence_campaign.json)
are the frozen record.  This test re-runs a representative beta subset on
all four routes and asserts each stays within its frozen tolerance —
catching silent route drift after refactors.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from sagnac_reference.analytic import directed_times  # noqa: E402
from sagnac_reference.segment_chain import (  # noqa: E402
    return_time_segment_chain,
)
from sagnac_reference.transport_pde import return_time_pde  # noqa: E402
from sagnac_reference.transport_pde_upwind import (  # noqa: E402
    return_time_pde_upwind,
)
from sagnac_reference.true_chain import return_time_true_chain  # noqa: E402

from sagnac_reference.analytic import SagnacConfig  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = ROOT / "artifacts/convergence/convergence_campaign.json"

TOL_PDE = 5e-3          # frozen evaluator tolerance
TOL_UPWIND = 4e-3       # upwind O(dx) at nx=1200, worst beta 0.8 (measured 1.45e-3 at beta=0.5, 4e-3 covers to beta=0.8 with margin)
TOL_CHAIN_TRUE = 1e-10
TOL_CHAIN_SERIES = 1e-4

BETAS = [0.0, 0.1, 0.2, 0.3, 0.5, -0.2, -0.5]


@pytest.mark.parametrize("beta", BETAS)
def test_all_routes_within_frozen_tolerances(beta):
    cfg = SagnacConfig(1.0, 1.0, beta)
    tp, tm = directed_times(cfg)

    cp = return_time_segment_chain(cfg, +1).return_time
    cm = return_time_segment_chain(cfg, -1).return_time
    assert abs(cp - tp) < TOL_CHAIN_SERIES
    assert abs(cm - tm) < TOL_CHAIN_SERIES

    if beta != 0.0:
        tp_true = return_time_true_chain(cfg, +1, steps=32000).return_time
        tm_true = return_time_true_chain(cfg, -1, steps=32000).return_time
        assert abs(tp_true - tp) < TOL_CHAIN_TRUE
        assert abs(tm_true - tm) < TOL_CHAIN_TRUE

        pp = return_time_pde(cfg, +1, nx=1200).return_time
        pm = return_time_pde(cfg, -1, nx=1200).return_time
        assert abs(pp - tp) < TOL_PDE
        assert abs(pm - tm) < TOL_PDE

        up = return_time_pde_upwind(cfg, +1, nx=1200).return_time
        um = return_time_pde_upwind(cfg, -1, nx=1200).return_time
        assert abs(up - tp) < TOL_UPWIND
        assert abs(um - tm) < TOL_UPWIND


def test_campaign_artifact_exists_and_frozen():
    assert CAMPAIGN.exists(), "frozen campaign must remain in repo"
    data = json.loads(CAMPAIGN.read_text())
    assert isinstance(data, (dict, list)) and data
