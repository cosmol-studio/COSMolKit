from __future__ import annotations

import importlib.util
import json
import tempfile
from pathlib import Path

TOOL_DIR = Path(__file__).resolve().parents[1]
ROOT = TOOL_DIR.parents[2]


def load_runner():
    specification = importlib.util.spec_from_file_location(
        "chembl_parity_tautomer_runner", TOOL_DIR / "run.py"
    )
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


runner = load_runner()


def load_auditor():
    specification = importlib.util.spec_from_file_location(
        "chembl_parity_tautomer_auditor", TOOL_DIR / "audit_tautomer.py"
    )
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


audit_tautomer = load_auditor()


def test_complete_profile_owns_full_tautomer_phase() -> None:
    profile = runner.load_json(TOOL_DIR / "profiles/complete.json")
    runner.validate_profile(profile)
    phase = {item["name"]: item for item in profile["phases"]}["tautomer"]
    assert phase["script"] == "audit_tautomer.py"
    assert phase["mode"] == "tautomer"
    assert phase["expected_processed"] == profile["corpus_records"]
    assert phase["expected_profiles"] == {
        "tautomer_parse": 1,
        "tautomer_default": 9,
        "tautomer_v1": 9,
        "tautomer_max_tautomers_1": 9,
        "tautomer_max_transforms_1": 9,
        "tautomer_retain_sp3_stereo": 9,
        "tautomer_retain_bond_stereo": 9,
        "tautomer_retain_isotopic_hydrogens": 9,
        "tautomer_no_reassign_stereo": 9,
    }


def test_runner_passes_pinned_profile_and_rust_oracle() -> None:
    profile = runner.load_json(TOOL_DIR / "profiles/complete.json")
    phase = {item["name"]: item for item in profile["phases"]}["tautomer"]
    command = runner.command_for(
        ROOT,
        phase,
        Path("corpus/shard-003.jsonl"),
        Path("run/tautomer/shard-003.json"),
        100,
        5,
    )
    assert command[command.index("--mode") + 1] == "tautomer"
    assert Path(command[command.index("--tautomer-profile") + 1]).name == "tautomer_profile.json"
    assert Path(command[command.index("--cosmolkit-oracle") + 1]).name.startswith(
        "cosmolkit-chembl-tautomer-oracle"
    )


def test_audit_compares_every_declared_field_and_records_no_hidden_skip() -> None:
    record = {"row": 0, "chembl_id": "HARNESS", "smiles": "CC(C)=O"}
    branch = {
        "ok": True,
        "error": None,
        "ordered_smiles": ["C=C(C)O", "CC(C)=O"],
        "molecule_states": [{"isomeric_smiles": "C=C(C)O"}],
        "status": "Completed",
        "modified_atoms": [0, 1, 3],
        "modified_bonds": [0, 1, 2],
        "scores": [{"total": 1}],
        "canonical_smiles": "CC(C)=O",
        "canonical_state": {"isomeric_smiles": "CC(C)=O"},
    }
    reference = {"parse": {"ok": True}, "branches": {"default": branch}}
    actual = {
        "parse": {"ok": True},
        "branches": {"default": {"result": branch}},
    }
    fields = [
        "ordered_smiles",
        "molecule_states",
        "status",
        "modified_atoms",
        "modified_bonds",
        "scores",
        "canonical_smiles",
        "canonical_state",
    ]
    with tempfile.TemporaryDirectory() as temporary:
        audit = audit_tautomer.Audit(Path(temporary) / "summary.json", 4)
        audit_tautomer.compare_record(
            audit, record, reference, actual, ["default"], fields
        )
        assert sum(audit.counts.values()) == 10
        assert set(audit.counts) == {
            "match.tautomer.parse",
            "match.tautomer.default.enumeration_outcome",
            *(f"match.tautomer.default.{field}" for field in fields),
        }
        audit.finish(1, {"tautomer_parse": 1, "tautomer_default": 9})
        summary = json.loads((Path(temporary) / "summary.json").read_text())
        assert summary["processed"] == 1
        assert summary["profiles"] == {
            "tautomer_parse": 1,
            "tautomer_default": 9,
        }


def test_auditor_streams_the_rust_oracle_without_a_temporary_result_file() -> None:
    source = (TOOL_DIR / "audit_tautomer.py").read_text(encoding="utf-8")
    assert '"--output",\n        "-"' in source
    assert "subprocess.Popen" in source
    assert "read_wire_message" in source
    assert "TemporaryDirectory" not in source
