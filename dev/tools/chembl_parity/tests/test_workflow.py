from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from dev.tools.chembl_parity import audit_surfaces


TOOL_DIR = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


runner = load_module("chembl_parity_run", TOOL_DIR / "run.py")
prepare = load_module("chembl_parity_prepare", TOOL_DIR / "prepare_corpus.py")


class WorkflowTests(unittest.TestCase):
    def test_profile_is_valid_and_uses_owned_scripts(self) -> None:
        profile = runner.load_json(TOOL_DIR / "profiles/complete.json")
        runner.validate_profile(profile)
        scripts = {phase["script"] for phase in profile["phases"]}
        self.assertEqual(scripts, runner.KNOWN_SCRIPTS)
        for script in scripts:
            self.assertTrue((TOOL_DIR / script).is_file())

    def test_complete_profile_contains_extended_chembl_surfaces(self) -> None:
        profile = runner.load_json(TOOL_DIR / "profiles/complete.json")
        phases = {phase["name"]: phase for phase in profile["phases"]}
        expected = {
            "topology-operations": "topology_operations",
            "fragments": "fragments",
            "coordinates-2d": "coordinates_2d",
            "binary-roundtrip": "binary_roundtrip",
        }
        for name, mode in expected.items():
            self.assertEqual(phases[name]["script"], "audit_surfaces.py")
            self.assertEqual(phases[name]["mode"], mode)
            self.assertIn(mode, runner.SCRIPT_MODES["audit_surfaces.py"])
        self.assertEqual(
            phases["topology-operations"]["expected_processed"], 2_854_376
        )

    def test_observation_errors_are_structured_values(self) -> None:
        def fail() -> None:
            raise ValueError("deliberate observation failure")

        self.assertEqual(
            audit_surfaces.capture_observation(fail),
            {
                "status": "error",
                "type": "ValueError",
                "message": "deliberate observation failure",
            },
        )

    def test_topology_mode_reaches_raw_parse_sanitize_rejection(self) -> None:
        record = {"chembl_id": "HARNESS", "row": 0, "smiles": "c"}
        rd_mol, ck_mol = audit_surfaces.parse_pair(record, sanitize=False)
        self.assertIsNotNone(rd_mol)
        self.assertIsNotNone(ck_mol)
        with tempfile.TemporaryDirectory() as temporary:
            audit = audit_surfaces.Audit(Path(temporary) / "summary.json", 5)
            audit_surfaces.audit_topology_operations(
                audit, record, rd_mol, ck_mol
            )
            self.assertEqual(audit.counts["match.topology.raw_parse.graph_state"], 1)
            self.assertEqual(audit.counts["match.topology.sanitize.rejected"], 1)
            self.assertFalse(
                any(key.startswith("mismatch.") for key in audit.counts)
            )
            audit.finish(1, "topology_operations", 0, 1)

    def test_prepare_corpus_is_deterministic_and_checksummed(self) -> None:
        rows = (
            "chembl_id\tcanonical_smiles\tstandard_inchi\tstandard_inchi_key\n"
            "CHEMBL1\tCCO\tI1\tK1\n"
            "CHEMBL2\t\tI2\tK2\n"
            "CHEMBL3\tN\tI3\tK3\n"
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "chembl.tsv.gz"
            with gzip.open(source, "wt", encoding="utf-8", newline="") as stream:
                stream.write(rows)
            source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
            output = root / "corpus"
            with mock.patch.object(
                prepare, "CHEMBL_37_SOURCE_SHA256", source_sha
            ), mock.patch(
                "sys.argv",
                [
                    "prepare_corpus.py",
                    "--input",
                    str(source),
                    "--output-dir",
                    str(output),
                    "--shards",
                    "2",
                ],
            ):
                prepare.main()
            manifest = json.loads((output / "manifest.json").read_text())
            self.assertEqual(manifest["records"], 2)
            self.assertEqual(manifest["missing_smiles"], 1)
            profile = {
                "corpus_source_sha256": source_sha,
                "corpus_records": 2,
            }
            validated, _ = runner.validate_corpus(output, profile, verify_files=True)
            self.assertEqual(validated["records"], 2)

    def test_aggregate_blocks_only_covered_mismatches(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "shard-000.json"
            output.write_text(
                json.dumps(
                    {
                        "processed": 3,
                        "counts": {
                            "match.graph": 3,
                            "mismatch.external.chembl_inchi_vs_rdkit": 1,
                        },
                    }
                )
            )
            result = {
                "0": {
                    "returncode": 0,
                    "output": str(output),
                }
            }
            aggregate = runner.aggregate_phase(
                root,
                result,
                1,
                {"mismatch.external.chembl_inchi_vs_rdkit"},
                3,
            )
            self.assertTrue(aggregate["passed"])
            summary = json.loads(output.read_text())
            summary["counts"]["mismatch.graph"] = 1
            output.write_text(json.dumps(summary))
            aggregate = runner.aggregate_phase(
                root,
                result,
                1,
                {"mismatch.external.chembl_inchi_vs_rdkit"},
                3,
            )
            self.assertFalse(aggregate["passed"])

    def test_corpus_validation_rejects_corrupt_shard(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            shard = root / "shard-000.jsonl"
            shard.write_text('{"chembl_id":"CHEMBL1","smiles":"CCO"}\n')
            checksum = hashlib.sha256(shard.read_bytes()).hexdigest()
            manifest = {
                "schema": runner.CORPUS_SCHEMA,
                "source_sha256": "source",
                "records": 1,
                "files": [
                    {
                        "shard": 0,
                        "file": shard.name,
                        "records": 1,
                        "bytes": shard.stat().st_size,
                        "sha256": checksum,
                    }
                ],
            }
            (root / "manifest.json").write_text(json.dumps(manifest))
            runner.validate_corpus(
                root,
                {"corpus_source_sha256": "source", "corpus_records": 1},
                verify_files=True,
            )
            shard.write_text("corrupt\n")
            with self.assertRaisesRegex(ValueError, "size mismatch|checksum mismatch"):
                runner.validate_corpus(
                    root,
                    {"corpus_source_sha256": "source", "corpus_records": 1},
                    verify_files=True,
                )

    def test_completed_task_requires_command_and_output_identity(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "summary.json"
            output.write_text(json.dumps({"processed": 1, "counts": {}}))
            command = ["python", "audit.py"]
            result = {
                "returncode": 0,
                "command": command,
                "output": str(output),
                "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
            }
            findings = output.with_suffix(".findings.jsonl")
            findings.write_text("")
            result["findings_sha256"] = hashlib.sha256(findings.read_bytes()).hexdigest()
            self.assertTrue(runner.completed_task_is_valid(result, command))
            self.assertFalse(
                runner.completed_task_is_valid(result, ["python", "other.py"])
            )
            output.write_text("{}")
            self.assertFalse(runner.completed_task_is_valid(result, command))

    def test_repository_outputs_must_be_ignored(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with mock.patch.object(runner.subprocess, "run") as run:
                run.return_value.returncode = 1
                with self.assertRaisesRegex(ValueError, "not ignored by Git"):
                    runner.require_external_or_ignored(root, root / "results", "run")
            with mock.patch.object(runner.subprocess, "run") as run:
                runner.require_external_or_ignored(
                    root, root.parent / "external-results", "run"
                )
                run.assert_not_called()


if __name__ == "__main__":
    unittest.main()
