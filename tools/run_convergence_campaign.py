#!/usr/bin/env python3
"""Phase 6 convergence campaign (regenerates artifacts/convergence/)."""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from sagnac_reference.analytic import SagnacConfig, directed_times
from sagnac_reference.transport_pde import return_time_pde
from sagnac_reference.true_chain import return_time_true_chain

BETAS = [0.0, 1e-3, -1e-3, 1e-2, -1e-2, 0.1, -0.1, 0.2, -0.2, 0.5, -0.5, 0.8, -0.8]


def main() -> int:
    out = Path("artifacts/convergence")
    out.mkdir(parents=True, exist_ok=True)
    L, c = 1.0, 1.0
    rows = []
    for b in BETAS:
        cfg = SagnacConfig(L, c, b * c)
        tp, tm = directed_times(cfg)
        row = {"beta": b, "t_plus_exact": tp, "t_minus_exact": tm}
        for n in (600, 1200, 2400):
            pp = return_time_pde(cfg, +1, nx=n).return_time
            pm = return_time_pde(cfg, -1, nx=n).return_time
            row[f"pde_plus_err_n{n}"] = abs(pp - tp)
            row[f"pde_minus_err_n{n}"] = abs(pm - tm)
        for s in (4000, 8000, 16000, 64000):
            rp = return_time_true_chain(cfg, +1, steps=s).return_time
            rm = return_time_true_chain(cfg, -1, steps=s).return_time
            row[f"chain_plus_err_s{s}"] = abs(rp - tp)
            row[f"chain_minus_err_s{s}"] = abs(rm - tm)
        rows.append(row)
        print(f"beta={b:+.3f} done")
    (out / "convergence_campaign.json").write_text(json.dumps(rows, indent=1))
    with open(out / "convergence_campaign.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("written:", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
