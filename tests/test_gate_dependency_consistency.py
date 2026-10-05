"""Gate dependencies must be transitively satisfied: a dependent gate may
only PASS if all its declared dependencies PASS (metadata becomes a
contract, not decoration)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from sagnac_reference.diagnostics import write_verdict


def test_dependencies_are_transitively_satisfied(tmp_path):
    verdict = write_verdict(tmp_path / "verdict.json")
    assert verdict["status"] == "SAGNAC_REFERENCE_CLOSURE_PASS"
    gates = verdict["gates"]
    graph = verdict["dependency_graph"]
    detailed = verdict["gates_detailed"]

    def passes(short: str) -> bool:
        return any(g.startswith(short + "-") and v for g, v in gates.items())

    for short, deps in graph.items():
        for dep in deps:
            assert passes(dep), (
                f"{short} declares dependency {dep} but {dep} does not pass")
        # cross-check metadata consistency: the depends_on in gates_detailed
        # matches the graph
        key = next(k for k in detailed if k.startswith(short + "-"))
        assert detailed[key].get("depends_on") == deps, (
            f"{key}: gates_detailed.depends_on != dependency_graph entry")
