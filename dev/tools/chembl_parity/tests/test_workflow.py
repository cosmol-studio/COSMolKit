from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from dev.tools.chembl_parity import audit_fingerprints, audit_stereo, audit_surfaces


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

    def test_molalign_phase_uses_the_complete_supported_coordinate_boundary(self) -> None:
        profile = runner.load_json(TOOL_DIR / "profiles/complete.json")
        phase = {item["name"]: item for item in profile["phases"]}["molalign"]
        self.assertEqual(phase["script"], "audit_surfaces.py")
        self.assertEqual(phase["mode"], "molalign")
        self.assertEqual(phase["max_atoms"], 80)
        self.assertEqual(phase["expected_processed"], 2_854_362)
        command = runner.command_for(
            TOOL_DIR.parents[2],
            phase,
            Path("corpus/shard-017.jsonl"),
            Path("run/molalign/shard-017.json"),
            123,
            4,
        )
        self.assertEqual(command[command.index("--mode") + 1], "molalign")
        self.assertEqual(command[command.index("--max-atoms") + 1], "80")
        self.assertEqual(command[command.index("--limit") + 1], "123")

    def test_molalign_auditor_exercises_every_source_call_family(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            audit = audit_surfaces.Audit(Path(temporary) / "summary.json", 4)
            for row in range(6):
                record = {
                    "row": row,
                    "chembl_id": f"HARNESS{row}",
                    "smiles": "CCCO",
                }
                rd_mol = audit_surfaces.Chem.MolFromSmiles(record["smiles"])
                ck_mol = audit_surfaces.cosmolkit.Molecule.from_smiles(
                    record["smiles"]
                )
                self.assertIsNotNone(rd_mol)
                audit_surfaces.audit_molalign(audit, record, rd_mol, ck_mol)
            self.assertFalse(
                [key for key in audit.counts if key.startswith("mismatch.")]
            )
            required = {
                "match.molalign.alignment_transform.rmsd",
                "match.molalign.best_alignment.rmsd",
                "match.molalign.coordinate_rmsd.rmsd",
                "match.molalign.align_to.coordinates.conformer_0",
                "match.molalign.all_conformer_best_rms.rmsds",
                "match.molalign.align_conformers.rmsds",
            }
            self.assertTrue(required.issubset(audit.counts))
            audit.finish(6, "molalign", 0, 1)

    def test_stereo_phase_is_complete_source_bounded_and_unfiltered(self) -> None:
        profile = runner.load_json(TOOL_DIR / "profiles/complete.json")
        phase = {item["name"]: item for item in profile["phases"]}[
            "stereo-enumeration"
        ]
        self.assertEqual(phase["script"], "audit_stereo.py")
        self.assertEqual(phase["mode"], "stereo")
        self.assertEqual(phase["expected_processed"], profile["corpus_records"])
        self.assertEqual(
            phase["expected_profiles"],
            {"potential_stereo": 4, "enumeration": 4},
        )
        self.assertNotIn("max_atoms", phase)
        self.assertNotIn("records_per_shard", phase)

        self.assertEqual(
            audit_stereo.POTENTIAL_PROFILES,
            (
                {"id": "preserve_possible", "clean": False, "flag_possible": True},
                {"id": "clean_possible", "clean": True, "flag_possible": True},
                {"id": "preserve_known", "clean": False, "flag_possible": False},
                {"id": "clean_known", "clean": True, "flag_possible": False},
            ),
        )
        self.assertEqual(
            [item["id"] for item in audit_stereo.ENUMERATION_PROFILES],
            [
                "default_bounded",
                "all_assigned_bounded",
                "non_unique_bounded",
                "seeded_three",
            ],
        )
        self.assertEqual(
            [item["max_isomers"] for item in audit_stereo.ENUMERATION_PROFILES],
            [8, 8, 8, 3],
        )

        command = runner.command_for(
            TOOL_DIR.parents[2],
            phase,
            Path("corpus/shard-017.jsonl"),
            Path("run/stereo-enumeration/shard-017.json"),
            123,
            4,
        )
        self.assertEqual(command[command.index("--mode") + 1], "stereo")
        self.assertNotIn("--limit", command)
        self.assertNotIn("--max-atoms", command)

    def test_stereo_auditor_compares_every_profile_and_exact_state_field(self) -> None:
        record = {
            "row": 0,
            "chembl_id": "HARNESS",
            "smiles": "CC(F)C(Cl)Br",
        }
        with tempfile.TemporaryDirectory() as temporary:
            audit = audit_stereo.Audit(Path(temporary) / "summary.json", 4)
            audit_stereo.audit_record(audit, record)
            audit.finish(1, "stereo", 0, 1)
            self.assertFalse(
                any(key.startswith("mismatch.") for key in audit.counts),
                audit.counts,
            )
            self.assertEqual(sum(audit.counts.values()), 31)
            expected = {
                *(f"match.potential.{profile['id']}.records" for profile in audit_stereo.POTENTIAL_PROFILES),
                *(f"match.potential.{profile['id']}.state" for profile in audit_stereo.POTENTIAL_PROFILES),
                *(f"match.enumeration.{profile['id']}.theoretical_count" for profile in audit_stereo.ENUMERATION_PROFILES),
                *(f"match.enumeration.{profile['id']}.bounded_out" for profile in audit_stereo.ENUMERATION_PROFILES),
                *(f"match.enumeration.{profile['id']}.outputs" for profile in audit_stereo.ENUMERATION_PROFILES),
            }
            self.assertTrue(expected <= set(audit.counts))

    def test_stereo_auditor_treats_source_none_and_structured_rejection_equally(self) -> None:
        record = {
            "row": 0,
            "chembl_id": "SOURCE_REJECTED",
            "smiles": "F[PH](F)(F)(F)(F)F",
        }
        with tempfile.TemporaryDirectory() as temporary:
            audit = audit_stereo.Audit(Path(temporary) / "summary.json", 4)
            audit_stereo.audit_record(audit, record)
            audit.finish(1, "stereo", 0, 1)

            self.assertEqual(audit.counts["match.parse.accepted"], 1)
            self.assertEqual(audit.counts["parse.both_rejected"], 1)
            self.assertFalse(
                any(key.startswith("mismatch.") for key in audit.counts),
                audit.counts,
            )

    def test_stereo_aggregate_requires_all_shards_profiles_and_rows(self) -> None:
        expected_profiles = {"potential_stereo": 4, "enumeration": 4}
        with tempfile.TemporaryDirectory() as temporary:
            phase_dir = Path(temporary)
            results = {}
            for shard in range(3):
                output = phase_dir / f"shard-{shard:03d}.json"
                output.write_text(
                    json.dumps(
                        {
                            "processed": shard + 1,
                            "profiles": expected_profiles,
                            "counts": {"match.enumeration.default_bounded.outputs": shard + 1},
                        }
                    )
                )
                output.with_suffix(".findings.jsonl").write_text("")
                results[str(shard)] = {"returncode": 0, "output": str(output)}

            aggregate = runner.aggregate_phase(
                phase_dir, results, 3, set(), 6, expected_profiles
            )
            self.assertTrue(aggregate["complete"])
            self.assertTrue(aggregate["passed"])

            partial = phase_dir / "shard-002.json"
            summary = json.loads(partial.read_text())
            summary["processed"] = 0
            partial.write_text(json.dumps(summary))
            aggregate = runner.aggregate_phase(
                phase_dir, results, 3, set(), 6, expected_profiles
            )
            self.assertFalse(aggregate["complete"])
            self.assertFalse(aggregate["passed"])

            summary["processed"] = 3
            summary["profiles"] = {"potential_stereo": 4, "enumeration": 3}
            partial.write_text(json.dumps(summary))
            aggregate = runner.aggregate_phase(
                phase_dir, results, 3, set(), 6, expected_profiles
            )
            self.assertFalse(aggregate["complete"])
            self.assertEqual(aggregate["invalid_profile_tasks"], 1)

            del results["2"]
            aggregate = runner.aggregate_phase(
                phase_dir, results, 3, set(), 6, expected_profiles
            )
            self.assertFalse(aggregate["complete"])
            self.assertEqual(aggregate["failed_or_missing_tasks"], 1)

    def test_stereo_script_and_profile_are_part_of_resume_identity(self) -> None:
        root = TOOL_DIR.parents[2]
        profile = runner.load_json(TOOL_DIR / "profiles/complete.json")
        with mock.patch.object(runner, "extension_artifacts", return_value=[]), mock.patch.object(
            runner, "git_head", return_value="HEAD"
        ), mock.patch.object(runner, "distribution_version", return_value="0"):
            identity = runner.build_identity(
                root,
                TOOL_DIR,
                TOOL_DIR / "profiles/complete.json",
                profile,
                {"source_sha256": "source"},
                "manifest-sha",
                b"",
                b"",
            )
        self.assertIn("audit_stereo.py", identity["scripts"])
        self.assertEqual(
            identity["selected_phases"][-1], "stereo-enumeration"
        )

        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "shard-000.json"
            findings = output.with_suffix(".findings.jsonl")
            output.write_text(json.dumps({"processed": 1, "counts": {}}))
            findings.write_text("")
            command = ["python", "-m", "audit_stereo"]
            result = {
                "returncode": 0,
                "command": command,
                "output": str(output),
                "output_sha256": runner.sha256_file(output),
                "findings_sha256": runner.sha256_file(findings),
            }
            self.assertTrue(runner.completed_task_is_valid(result, command))
            output.write_text(json.dumps({"processed": 0, "counts": {}}))
            self.assertFalse(runner.completed_task_is_valid(result, command))

    def test_ciplabeler_phase_is_complete_and_source_backed(self) -> None:
        profile = runner.load_json(TOOL_DIR / "profiles/complete.json")
        phase = {item["name"]: item for item in profile["phases"]}["ciplabeler"]
        self.assertEqual(phase["script"], "audit_core.py")
        self.assertEqual(phase["mode"], "ciplabeler")
        self.assertEqual(phase["expected_processed"], 2_897_819)
        command = runner.command_for(
            TOOL_DIR.parents[2],
            phase,
            Path("corpus/shard-017.jsonl"),
            Path("run/ciplabeler/shard-017.json"),
            123,
            4,
        )
        self.assertEqual(command[command.index("--mode") + 1], "ciplabeler")
        self.assertNotIn("--max-atoms", command)

    def test_atom_pair_phase_is_complete_and_source_profiled(self) -> None:
        profile = runner.load_json(TOOL_DIR / "profiles/complete.json")
        phases = {phase["name"]: phase for phase in profile["phases"]}
        phase = phases["atom-pair"]
        self.assertEqual(phase["script"], "audit_fingerprints.py")
        self.assertEqual(phase["mode"], "atom_pair")
        self.assertEqual(phase["expected_processed"], profile["corpus_records"])
        self.assertEqual(
            phase["expected_profiles"],
            {
                "atom_pair_sparse_count": 10,
                "atom_pair_count": 10,
                "atom_pair_sparse_bit": 10,
                "atom_pair_explicit_bit": 10,
                "atom_pair_additional_output": 1,
            },
        )
        command = runner.command_for(
            TOOL_DIR.parents[2],
            phase,
            Path("corpus/shard-017.jsonl"),
            Path("run/atom-pair/shard-017.json"),
            123,
            4,
        )
        self.assertIn("corpus/shard-017.jsonl", command)
        self.assertIn("run/atom-pair/shard-017.json", command)
        self.assertEqual(command[command.index("--mode") + 1], "atom_pair")
        atom_profile = Path(
            command[command.index("--atom-pair-profile") + 1]
        )
        self.assertEqual(
            atom_profile.name, "atom_pair_fingerprint_profile.json"
        )

    def test_atom_pair_auditor_compares_every_configured_exact_field(self) -> None:
        profile = runner.load_json(
            TOOL_DIR.parents[2]
            / "tools/testdata/rdkit/atom_pair_fingerprint_profile.json"
        )
        by_name = {branch["name"]: branch for branch in profile["branches"]}
        branches = [by_name["default"], by_name["additional_output"]]
        generators = audit_fingerprints.make_atom_pair_generators(branches)
        record = {"row": 0, "chembl_id": "HARNESS", "smiles": "CCCO"}
        rd_mol = audit_fingerprints.Chem.MolFromSmiles(record["smiles"])
        ck_mol = audit_fingerprints.cosmolkit.Molecule.from_smiles(record["smiles"])
        self.assertIsNotNone(rd_mol)
        with tempfile.TemporaryDirectory() as temporary:
            audit = audit_fingerprints.Audit(Path(temporary) / "summary.json", 4)
            audit_fingerprints.compare_atom_pair(
                audit, record, rd_mol, ck_mol, branches, generators
            )
            self.assertEqual(sum(audit.counts.values()), 9)
            self.assertTrue(
                all(key.startswith("match.atom_pair.") for key in audit.counts)
            )
            audit.finish(
                1,
                {
                    "atom_pair_sparse_count": 2,
                    "atom_pair_count": 2,
                    "atom_pair_sparse_bit": 2,
                    "atom_pair_explicit_bit": 2,
                    "atom_pair_additional_output": 1,
                },
            )

    def test_topological_torsion_phase_is_complete_and_source_profiled(self) -> None:
        profile = runner.load_json(TOOL_DIR / "profiles/complete.json")
        phases = {phase["name"]: phase for phase in profile["phases"]}
        phase = phases["topological-torsion"]
        self.assertEqual(phase["script"], "audit_fingerprints.py")
        self.assertEqual(phase["mode"], "topological_torsion")
        self.assertEqual(phase["expected_processed"], profile["corpus_records"])
        self.assertEqual(
            phase["expected_profiles"],
            {
                "topological_torsion_sparse_count": 9,
                "topological_torsion_count": 9,
                "topological_torsion_sparse_bit": 9,
                "topological_torsion_explicit_bit": 9,
                "topological_torsion_sparse_count_additional_output": 2,
                "topological_torsion_count_additional_output": 2,
                "topological_torsion_sparse_bit_additional_output": 2,
                "topological_torsion_explicit_bit_additional_output": 2,
            },
        )
        command = runner.command_for(
            TOOL_DIR.parents[2],
            phase,
            Path("corpus/shard-017.jsonl"),
            Path("run/topological-torsion/shard-017.json"),
            123,
            4,
        )
        self.assertEqual(
            command[command.index("--mode") + 1], "topological_torsion"
        )
        torsion_profile = Path(
            command[command.index("--topological-torsion-profile") + 1]
        )
        self.assertEqual(
            torsion_profile.name, "topological_torsion_fingerprint_profile.json"
        )

    def test_layered_phase_is_complete_and_source_profiled(self) -> None:
        profile = runner.load_json(TOOL_DIR / "profiles/complete.json")
        phase = {item["name"]: item for item in profile["phases"]}["layered"]
        self.assertEqual(phase["script"], "audit_fingerprints.py")
        self.assertEqual(phase["mode"], "layered")
        self.assertEqual(phase["expected_processed"], profile["corpus_records"])
        self.assertEqual(phase["expected_profiles"], {"layered": 18})
        command = runner.command_for(
            TOOL_DIR.parents[2],
            phase,
            Path("corpus/shard-017.jsonl"),
            Path("run/layered/shard-017.json"),
            123,
            4,
        )
        self.assertEqual(command[command.index("--mode") + 1], "layered")
        layered_profile = Path(command[command.index("--layered-profile") + 1])
        self.assertEqual(layered_profile.name, "layered_fingerprint_profile.json")

    def test_layered_auditor_compares_every_configured_exact_field(self) -> None:
        profile = runner.load_json(
            TOOL_DIR.parents[2]
            / "tools/testdata/rdkit/layered_fingerprint_profile.json"
        )
        branches = audit_fingerprints.layered_branches(profile)
        self.assertEqual(len(branches), 18)
        by_name = {branch["name"]: branch for branch in branches}
        selected = [
            by_name["default"],
            by_name["rooted_terminal_pair_branched"],
            by_name["mod_three_mask_seeded_counts"],
            by_name["high_source_bits_only"],
        ]
        record = {"row": 0, "chembl_id": "HARNESS", "smiles": "c1ccccc1O"}
        rd_mol = audit_fingerprints.Chem.MolFromSmiles(record["smiles"])
        ck_mol = audit_fingerprints.cosmolkit.Molecule.from_smiles(record["smiles"])
        self.assertIsNotNone(rd_mol)
        with tempfile.TemporaryDirectory() as temporary:
            audit = audit_fingerprints.Audit(Path(temporary) / "summary.json", 4)
            audit_fingerprints.compare_layered(
                audit, record, rd_mol, ck_mol, selected
            )
            self.assertEqual(sum(audit.counts.values()), len(selected))
            self.assertTrue(
                all(key.startswith("match.layered.") for key in audit.counts)
            )
            audit.finish(1, {"layered": len(selected)})

    def test_pattern_phase_is_complete_source_profiled_and_sharded(self) -> None:
        profile = runner.load_json(TOOL_DIR / "profiles/complete.json")
        phase = {item["name"]: item for item in profile["phases"]}["pattern"]
        self.assertEqual(phase["script"], "audit_fingerprints.py")
        self.assertEqual(phase["mode"], "pattern")
        self.assertEqual(phase["expected_processed"], profile["corpus_records"])
        self.assertEqual(phase["expected_profiles"], {"pattern": 10})
        command = runner.command_for(
            TOOL_DIR.parents[2],
            phase,
            Path("corpus/shard-017.jsonl"),
            Path("run/pattern/shard-017.json"),
            123,
            4,
        )
        self.assertEqual(command[command.index("--mode") + 1], "pattern")
        pattern_profile = Path(command[command.index("--pattern-profile") + 1])
        self.assertEqual(pattern_profile.name, "pattern_fingerprint_profile.json")

    def test_pattern_auditor_compares_every_complete_corpus_profile_exactly(self) -> None:
        profile = runner.load_json(
            TOOL_DIR.parents[2]
            / "tools/testdata/rdkit/pattern_fingerprint_profile.json"
        )
        branches = audit_fingerprints.validate_pattern_profile(profile)
        self.assertEqual(len(branches), 10)
        self.assertNotIn(
            "set_only_bits_wrong_width",
            {branch["name"] for branch in branches},
        )
        record = {"row": 0, "chembl_id": "HARNESS", "smiles": "CCCO"}
        rd_mol = audit_fingerprints.Chem.MolFromSmiles(record["smiles"])
        ck_mol = audit_fingerprints.cosmolkit.Molecule.from_smiles(record["smiles"])
        self.assertIsNotNone(rd_mol)
        with tempfile.TemporaryDirectory() as temporary:
            audit = audit_fingerprints.Audit(Path(temporary) / "summary.json", 4)
            audit_fingerprints.compare_pattern(
                audit, record, rd_mol, ck_mol, branches
            )
            self.assertEqual(sum(audit.counts.values()), 10)
            self.assertTrue(
                all(key.startswith("match.pattern.") for key in audit.counts)
            )
            audit.finish(1, {"pattern": 10})

    def test_pattern_aggregate_rejects_partial_profile_and_missing_shard(self) -> None:
        expected_profiles = {"pattern": 10}
        with tempfile.TemporaryDirectory() as temporary:
            phase_dir = Path(temporary)
            results = {}
            for shard in range(2):
                output = phase_dir / f"shard-{shard:03d}.json"
                output.write_text(
                    json.dumps(
                        {
                            "processed": 1,
                            "profiles": expected_profiles,
                            "counts": {"match.pattern.default": 1},
                        }
                    )
                )
                output.with_suffix(".findings.jsonl").write_text("")
                results[str(shard)] = {"returncode": 0, "output": str(output)}
            aggregate = runner.aggregate_phase(
                phase_dir, results, 2, set(), 2, expected_profiles
            )
            self.assertTrue(aggregate["passed"])

            summary = json.loads((phase_dir / "shard-001.json").read_text())
            summary["profiles"] = {"pattern": 9}
            (phase_dir / "shard-001.json").write_text(json.dumps(summary))
            aggregate = runner.aggregate_phase(
                phase_dir, results, 2, set(), 2, expected_profiles
            )
            self.assertFalse(aggregate["complete"])
            self.assertEqual(aggregate["invalid_profile_tasks"], 1)

            del results["1"]
            aggregate = runner.aggregate_phase(
                phase_dir, results, 2, set(), 2, expected_profiles
            )
            self.assertFalse(aggregate["complete"])

    def test_pattern_profile_is_part_of_resume_identity(self) -> None:
        root = TOOL_DIR.parents[2]
        profile = runner.load_json(TOOL_DIR / "profiles/complete.json")
        with mock.patch.object(runner, "extension_artifacts", return_value=[]), mock.patch.object(
            runner, "git_head", return_value="HEAD"
        ), mock.patch.object(runner, "distribution_version", return_value="0"):
            identity = runner.build_identity(
                root,
                TOOL_DIR,
                TOOL_DIR / "profiles/complete.json",
                profile,
                {"source_sha256": "source"},
                "manifest-sha",
                b"",
                b"",
            )
        self.assertIn(
            "tools/testdata/rdkit/pattern_fingerprint_profile.json",
            identity["reference_profiles"],
        )

    def test_topological_torsion_auditor_compares_all_vectors_and_provenance(self) -> None:
        profile = runner.load_json(
            TOOL_DIR.parents[2]
            / "tools/testdata/rdkit/topological_torsion_fingerprint_profile.json"
        )
        by_name = {branch["name"]: branch for branch in profile["corpus_branches"]}
        branches = [by_name["default"], by_name["all_provenance"]]
        rdkit_generators, cosmolkit_generators = (
            audit_fingerprints.make_topological_torsion_generators(branches)
        )
        record = {"row": 0, "chembl_id": "HARNESS", "smiles": "CCCO"}
        rd_mol = audit_fingerprints.Chem.MolFromSmiles(record["smiles"])
        ck_mol = audit_fingerprints.cosmolkit.Molecule.from_smiles(record["smiles"])
        self.assertIsNotNone(rd_mol)
        with tempfile.TemporaryDirectory() as temporary:
            audit = audit_fingerprints.Audit(Path(temporary) / "summary.json", 4)
            audit_fingerprints.compare_topological_torsion(
                audit,
                record,
                rd_mol,
                ck_mol,
                branches,
                rdkit_generators,
                cosmolkit_generators,
            )
            self.assertEqual(sum(audit.counts.values()), 12)
            self.assertTrue(
                all(
                    key.startswith("match.topological_torsion.")
                    for key in audit.counts
                )
            )
            audit.finish(
                1,
                {
                    "topological_torsion_sparse_count": 2,
                    "topological_torsion_count": 2,
                    "topological_torsion_sparse_bit": 2,
                    "topological_torsion_explicit_bit": 2,
                    "topological_torsion_sparse_count_additional_output": 1,
                    "topological_torsion_count_additional_output": 1,
                    "topological_torsion_sparse_bit_additional_output": 1,
                    "topological_torsion_explicit_bit_additional_output": 1,
                },
            )

    def test_atom_pair_aggregate_rejects_filtered_or_partial_shards(self) -> None:
        expected_profiles = {
            "atom_pair_sparse_count": 10,
            "atom_pair_count": 10,
            "atom_pair_sparse_bit": 10,
            "atom_pair_explicit_bit": 10,
            "atom_pair_additional_output": 1,
        }
        with tempfile.TemporaryDirectory() as temporary:
            phase_dir = Path(temporary)
            results = {}
            for shard in range(2):
                output = phase_dir / f"shard-{shard:03d}.json"
                output.write_text(
                    json.dumps(
                        {
                            "processed": 1,
                            "profiles": expected_profiles,
                            "counts": {"match.atom_pair.default.explicit_bit": 1},
                        }
                    )
                )
                output.with_suffix(".findings.jsonl").write_text("")
                results[str(shard)] = {"returncode": 0, "output": str(output)}
            aggregate = runner.aggregate_phase(
                phase_dir, results, 2, set(), 2, expected_profiles
            )
            self.assertTrue(aggregate["complete"])
            self.assertTrue(aggregate["passed"])

            filtered = phase_dir / "shard-001.json"
            summary = json.loads(filtered.read_text())
            summary["profiles"] = {"atom_pair_explicit_bit": 1}
            filtered.write_text(json.dumps(summary))
            aggregate = runner.aggregate_phase(
                phase_dir, results, 2, set(), 2, expected_profiles
            )
            self.assertFalse(aggregate["complete"])
            self.assertFalse(aggregate["passed"])
            self.assertEqual(aggregate["invalid_profile_tasks"], 1)

            del results["1"]
            aggregate = runner.aggregate_phase(
                phase_dir, results, 2, set(), 2, expected_profiles
            )
            self.assertFalse(aggregate["complete"])

    def test_atom_pair_run_resume_requires_identical_phase_identity(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run_dir = Path(temporary) / "run"
            identity = {
                "selected_phases": ["atom-pair"],
                "reference_profiles": {
                    "tools/testdata/rdkit/atom_pair_fingerprint_profile.json": "abc"
                },
            }
            profile = {"name": "atom-pair-harness"}
            manifest = runner.initialize_manifest(
                run_dir, identity, profile, workers=2, resume=False
            )
            runner.atomic_json(run_dir / "manifest.json", manifest)
            resumed = runner.initialize_manifest(
                run_dir, identity, profile, workers=4, resume=True
            )
            self.assertEqual(resumed["identity"], identity)
            self.assertEqual(resumed["workers"], 4)
            with self.assertRaisesRegex(ValueError, "run identity changed"):
                runner.initialize_manifest(
                    run_dir,
                    {**identity, "selected_phases": ["other"]},
                    profile,
                    workers=4,
                    resume=True,
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
