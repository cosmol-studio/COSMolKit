"""Original 152 conformer library rows through the public Python API.

COSMOLKIT_CONFORMER_ORACLE_152 selects pinned JSONL inputs;
the default is the project's smiles_small oracle.
"""
import hashlib
import json
import math
import os
import sys
import tempfile
from pathlib import Path

import cosmolkit as ck

WORKSPACE = Path(__file__).resolve().parents[2]


def _oracle(count, corpus):
    key = f"COSMOLKIT_CONFORMER_ORACLE_{count}"
    path = Path(os.environ.get(
        key,
        WORKSPACE / "testdata/conformer/expected/rdkit" / corpus / "conformer_generation_library.jsonl",
    ))
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing original{count} oracle: {path}; set {key} to the pinned "
            "original JSONL file. Do not regenerate the reference."
        )
    return path


COORD_TOLERANCE = 1.0e-6


def _execute(profile, path, digest, count):
    raw = path.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == digest
    rows = [json.loads(line) for line in raw.splitlines()]
    assert len(rows) == count
    # Preflight all source-defined settings before any chemistry call.
    for row in rows:
        assert (row["seed"], row["preset"], row["max_iterations"], row["timeout"]) == (61453, "ETKDGv3", 3, 0)
        assert row["rdkit_add_hs_ok"] or not row["rdkit_parse_ok"]
        assert row["rdkit_embed_ok"] == (row["status"] == 0)
    output = (Path(os.environ["COSMOLKIT_CONFORMER_NATIVE_OUTPUT"])
              if "COSMOLKIT_CONFORMER_NATIVE_OUTPUT" in os.environ
              else Path(tempfile.mkdtemp(prefix=f"cosmolkit-{profile}-")))
    output.mkdir(exist_ok=True)
    actual_path = output / (profile + "-actual.jsonl")
    assert not actual_path.exists()
    failures = []
    stages = {}
    with actual_path.open("x") as stream:
        for index, expected in enumerate(rows):
            actual = {"row": index + 1, "smiles": expected["smiles"], "parse_ok": False, "add_hs_ok": False, "embed_ok": False, "status": None, "coords": None, "error_stage": None, "error": None, "error_type": None}
            stage = "parse"
            try:
                molecule = ck.Molecule.from_smiles(expected["smiles"])
                actual["parse_ok"] = True
                stage = "add_hs"
                molecule = molecule.with_hydrogens()
                actual["add_hs_ok"] = True
                stage = "embed"
                params = ck.EmbedParams.etkdg_v3()
                params = configured(params, random_seed=61453)
                params = configured(params, num_threads=1)
                params = configured(params, max_iterations=3)
                params = configured(params, timeout=0)
                result = molecule.with_3d_conformer_result(params)
                actual["status"] = result.conf_id()
                actual["embed_ok"] = actual["status"] == 0
                actual["failure_counts"] = result.params().failures
                if actual["embed_ok"]:
                    conformers = result.molecule().conformers_3d()
                    actual["conformer_ids"] = [conformer.id() for conformer in conformers]
                    actual["coords"] = conformers[0].coordinates()
                else:
                    actual["error_stage"] = "embed"
                    actual["error"] = "EmbedMolecule returned " + str(actual["status"])
            except Exception as error:
                # Record and compare every exception; never skip the row or call it unsupported.
                actual["error_stage"] = stage
                actual["error"] = str(error)
                actual["error_type"] = type(error).__module__ + "." + type(error).__qualname__
            mismatch = []
            for field in ("parse_ok", "add_hs_ok", "embed_ok"):
                if actual[field] != expected["rdkit_" + field]:
                    mismatch.append(field)
            if actual["status"] != expected["status"]:
                mismatch.append("status")
            if actual["error_stage"] != expected["error_stage"]:
                mismatch.append("error_stage")
            # Preserve original parse-failure condition (error exists, not equal text).
            if not expected["rdkit_parse_ok"]:
                if not actual["error"] or not expected["error"]:
                    mismatch.append("missing_parse_error")
            elif actual["error"] != expected["error"]:
                mismatch.append("error")
            if expected["coords"] is not None:
                xyz = actual["coords"]
                if xyz is None or len(xyz) != len(expected["coords"]):
                    mismatch.append("coordinate_shape")
                else:
                    differences = [(i, axis, a[axis], e[axis]) for i, (a, e) in enumerate(zip(xyz, expected["coords"])) for axis in range(3) if not math.isfinite(a[axis]) or abs(a[axis] - e[axis]) > COORD_TOLERANCE]
                    actual["coordinate_mismatch_count"] = len(differences)
                    actual["first_coordinate_mismatch"] = differences[0] if differences else None
                    if differences:
                        mismatch.append("coordinates")
            elif actual["coords"] is not None:
                mismatch.append("unexpected_coordinates")
            actual["mismatches"] = mismatch
            actual["source_error"] = expected["error"]
            stream.write(json.dumps(actual, allow_nan=False) + "\n")
            stream.flush()
            stages[actual["error_stage"] or "success"] = stages.get(actual["error_stage"] or "success", 0) + 1
            if mismatch:
                failures.append({"row": index + 1, "smiles": expected["smiles"], "mismatches": mismatch})
            if (index + 1) % 100 == 0:
                print(profile, "actual rows", index + 1, "mismatches", len(failures), flush=True)
    native = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(ck.__file__).parent.glob("*.so")}
    summary = {"profile": profile, "rows_executed": count, "matched": count - len(failures), "failed": len(failures), "actual_stage_counts": stages, "failures": failures, "original_reference": str(path), "original_reference_sha256": digest, "actual_sha256": hashlib.sha256(actual_path.read_bytes()).hexdigest(), "native_sha256": native, "interpreter": sys.version, "source_pin": "351f8f378f8ad6bbd517980c38896e66bf907af8", "tolerance": COORD_TOLERANCE, "whole_block_ACCEPTED": False, "independent_review_required": True}
    (output / (profile + "-summary.json")).write_text(json.dumps(summary, indent=2) + "\n")
    assert not failures, f"{profile}: {len(failures)}/{count} native mismatches; all actual outcomes retained at {actual_path}"


def test_original_library_native_152():
    _execute("original152", _oracle(152, "smiles_small"), "536fe7d81cfcf913ad3901645bcf52ba1078f2e1faac7f897e861e3cdffcc7b3", 152)


def configured(params, **changes):
    values = {field: getattr(params, field) for field in ['max_iterations', 'num_threads', 'random_seed', 'clear_confs', 'use_random_coords', 'box_size_mult', 'rand_neg_eig', 'num_zero_fail', 'coord_map', 'optimizer_force_tol', 'ignore_smoothing_failures', 'enforce_chirality', 'use_exp_torsion_angle_prefs', 'use_basic_knowledge', 'verbose', 'basin_thresh', 'prune_rms_thresh', 'only_heavy_atoms_for_rms', 'et_version', 'embed_fragments_separately', 'use_small_ring_torsions', 'use_macrocycle_torsions', 'use_macrocycle14config', 'timeout', 'cpci', 'force_trans_amides', 'use_symmetry_for_pruning', 'bounds_mat_force_scaling', 'track_failures', 'enable_sequential_random_seeds', 'symmetrize_conjugated_terminal_groups_for_pruning']}
    values.update(changes)
    return ck.EmbedParams(**values)
