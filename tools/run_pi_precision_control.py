#!/usr/bin/env python3
"""SAGNAC PI-PRECISION CONTROL (numerical precision, NOT SSZ physics).

Lino's task 2: at an analytically known runtime (Sagnac reference),
truncate pi to d decimal places everywhere it enters the calculation and
measure the induced error against the full-precision reference.

Physics stays SSZ-independent. This is a NUMERICS control:
  - the Sagnac loop uses pi only via geometry (L = 2 pi R for a circular
    loop of radius R) — the standard construction in the repo's docs.
  - we evaluate the exact reference (mpmath, 50 digits) and then the same
    formula with pi truncated to d digits, for d = 2..30, and measure
    |delta_t(d) - delta_t(exact)| / delta_t(exact).
  - expected scaling: relative error ~ 10^-d / (2 ln-ish) i.e. one digit
    of pi per decade of precision. Verify the empirically measured
    digits-vs-error law against the analytic prediction |dU|/U = |dpi|/pi.
  - cross-check against SSZ_PI_DIGITS_RESOLUTION_LAW_V1: the digit count
    needed to push the runtime error below a given threshold epsilon is
    d > log10(U/epsilon_target)-ish — confirming that NUMERICAL precision
    is a free parameter, categorically different from physical resolution.
"""
import datetime
import json
import math
import subprocess
from pathlib import Path

from mpmath import mp, mpf
from mpmath import pi as mpi_pi

ROOT = Path("/home/error/physics/clones/Sagnac-Reference-Transport")

def head(p):
    return subprocess.run(["git", "-C", str(p), "rev-parse", "HEAD"],
                          capture_output=True, text=True).stdout.strip()

mp.dps = 50

R = mpf(1000.0)          # 1 km loop radius
c = mpf("299792458.0")
v = mpf("0.0")           # non-rotating reference: T_p = T_m = L/c (cleanest pi-test)
# L = 2 pi R (the ONLY place pi enters)
L_exact = 2 * mpi_pi * R
T_exact = L_exact / c

def pi_trunc(d):
    """pi truncated (not rounded) to d decimals via mpmath string ops."""
    mp.dps = d + 5
    s = mp.nstr(mpi_pi, d + 2)
    mp.dps = 50
    return mpf(s[: d + 2])  # "3." + d digits

rows = []
for d in range(1, 31):
    pi_d = pi_trunc(d)
    L_d = 2 * pi_d * R
    T_d = L_d / c
    rel = abs(T_d - T_exact) / T_exact
    # analytic prediction: relative runtime error = |dpi|/pi
    dpi = abs(pi_d - mpi_pi) / mpi_pi
    rows.append({"d": d,
                  "pi_used": float(pi_d),
                  "rel_err_runtime": float(rel),
                  "analytic_rel_err": float(dpi),
                  "match_order": bool(abs(math.log10(float(rel)) - math.log10(float(dpi))) < 0.5)})

# digits needed to push runtime error below thresholds
def digits_for(rel_target):
    # rel ~ |dpi|/pi ~ 10^-d / pi  -> d ~ log10(1/(pi*rel))
    return math.ceil(math.log10(1.0 / (math.pi * rel_target)))

thresholds = {"1e-6": digits_for(1e-6), "1e-12": digits_for(1e-12),
               "1e-16 (float64 eps)": digits_for(1e-16)}

out = {
    "artifact": "SAGNAC_PI_PRECISION_CONTROL_V1",
    "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
    "sagnac_head": head(ROOT),
    "scope": ("NUMERICS CONTROL ONLY. SSZ-independent: measures how pi-digit "
               "truncation propagates into the analytically known Sagnac/loop "
               "runtime. No segmentation claim, no theory content."),
    "setup": {"loop": "circular, L = 2 pi R", "R_m": float(R), "v_m_s": 0.0,
               "reference": "mpmath 50 digits"},
    "scaling_law": {
        "analytic": "rel_err(T) = |dpi|/pi  ~ 10^-d / pi",
        "empirical_matches_analytic_all_rows": all(r["match_order"] for r in rows),
    },
    "digits_for_target_relative_error": thresholds,
    "rows": rows,
    "interpretation": {
        "numerical_precision_is_free": ("digits of pi cost nothing computationally; "
                                         "float64 (16 digits) already puts runtime "
                                         "error at 1e-16 level for this geometry"),
        "distinct_from_physical_resolution": ("the physical-resolution law "
                                               "(SSZ_PI_DIGITS_RESOLUTION_LAW_V1: "
                                               "d > log10(R/l_min)) bounds what NATURE "
                                               "resolves; this control bounds what the "
                                               "CODE needs. Two different d's, two "
                                               "different questions — never mix them."),
        "crossover_example": ("a 1-km loop computed at float64 precision has runtime "
                               "error ~1e-16, while the physical-resolution floor for "
                               "R=1 km at l_min=lp would be d ~ 38 — the numerical "
                               "requirement is far below any physical claim"),
    },
    "no_observational_input": True,
    "observational_seal": "maintained",
}
out["closure_head"] = head("/home/error/physics/clones/SSZ_FULL_CLOSURE")
out_path = ROOT / "artifacts" / "SAGNAC_PI_PRECISION_CONTROL_V1.json"
out_path.write_text(json.dumps(out, indent=1) + "\n")
print("rows:", len(rows))
print("scaling law holds for all rows:", out["scaling_law"]["empirical_matches_analytic_all_rows"])
print("digits needed:", thresholds)
