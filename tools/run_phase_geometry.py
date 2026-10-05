#!/usr/bin/env python3
"""PHASE 10 - MATHEMATICAL_CONVERGENCE_GEOMETRY / NOT_PHYSICAL_LIGHT_OSCILLATION

Explorative, NICHT-normative Analyse der alternierenden Korrekturreihe als
gedämpfte diskrete Phasor-Dynamik:

    a_n = T0 * beta^n * exp(i*pi*n)      (Amplitude * Phase, Phase = n*pi)

Diese Geometrie lebt im ITERATIONSINDEX n, nicht in physikalischer Zeit.
Sie validiert Konvergenzstruktur — nichts physikalisch Schwingendes.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sagnac_reference.analytic import SagnacConfig  # noqa: E402

OUT = Path("artifacts/phase_geometry")


def analyse(beta: float, T0: float = 1.0, n_max: int = 24) -> dict:
    n = np.arange(n_max + 1)
    amp = T0 * beta**n
    phi = np.pi * n
    z = amp * np.exp(1j * phi)
    partial = np.cumsum(z)
    S_inf = T0 / (1.0 + beta)
    err = S_inf - partial
    return {
        "beta": beta,
        "S_inf_analytic": S_inf,
        "partial_sums_real": partial.real.tolist(),
        "partial_sums_imag": partial.imag.tolist(),
        "remainder_abs": np.abs(err).tolist(),
        "remainder_arg_pi": (np.angle(err) / np.pi).tolist(),
        "log_damping_check": {
            "ratio_last_two_terms": float(amp[-1] / amp[-2]),
            "equals_beta": bool(abs(amp[-1] / amp[-2] - beta) < 1e-15),
        },
        "phase_step_is_pi": bool(np.allclose(np.diff(phi) % (2 * np.pi), np.pi)),
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    out = {
        "label": "MATHEMATICAL_CONVERGENCE_GEOMETRY / NOT_PHYSICAL_LIGHT_OSCILLATION",
        "note": ("The alternating correction series is analysed as damped phasor "
                 "dynamics in the ITERATION index. This is convergence geometry, "
                 "not a claim about physical light motion."),
        "cases": {},
    }
    for beta in (0.2, 0.5, 0.8):
        out["cases"][f"beta_{beta}"] = analyse(beta)
        d = out["cases"][f"beta_{beta}"]
        ok = (d["log_damping_check"]["equals_beta"]
              and d["phase_step_is_pi"]
              and abs(float(d["remainder_abs"][-1]) - float(beta)**25 / (1 + float(beta))) < 1e-9)
        out["cases"][f"beta_{beta}"]["checks_pass"] = bool(ok)
    all_ok = all(c["checks_pass"] for c in out["cases"].values())
    out["all_checks_pass"] = all_ok
    (OUT / "phase_geometry.json").write_text(json.dumps(out, indent=1))
    print("label:", out["label"])
    print("all_checks_pass:", all_ok)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
