#!/usr/bin/env python3
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from sagnac_reference.diagnostics import write_verdict

out = write_verdict(ROOT / "artifacts" / "SAGNAC_REFERENCE_VERDICT.json")
print(out["status"])
for name, ok in out["gates"].items():
    print(f"{name}: {'PASS' if ok else 'FAIL'}")
raise SystemExit(0 if out["status"].endswith("PASS") else 1)
