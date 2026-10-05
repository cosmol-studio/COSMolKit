"""All original Pattern focused18/small152/5000 rows and all11 branches."""
import hashlib
import json
import os
from pathlib import Path
import cosmolkit as ck
import pytest

ROOT = Path(__file__).resolve().parents[2]
PROFILES = {
    "pattern_focused": (18, "COSMOLKIT_PATTERN_ORIGINAL_FOCUSED_GOLDEN", "0e5a58bc9cf7ce07836f71fc32dca80b25ae8b3cc038f4e729cf976b8528cd74", "testdata/fingerprint/fixtures/rdkit/pattern_fingerprint_focused.smi"),
    "smiles_small": (152, "COSMOLKIT_PATTERN_ORIGINAL_SMALL_GOLDEN", "05908b2765e0d7fa8665f5ac421dd82c8c3ad1f8e8799acde9d5e916c4199eab", "testdata/smiles/corpus/smiles_small.smi"),
    "smiles_5000": (5000, "COSMOLKIT_PATTERN_ORIGINAL5000_GOLDEN", "b4df21695a5b0caf71b5b4505f43c25f67782c4c72a248c9a7abce7f27767a85", "testdata/smiles/corpus/smiles_5000.smi"),
}


def input_rows(path):
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("smarts\t"):
            rows.append(("smarts", line.split("\t", 1)[1]))
        elif line.startswith("smiles\t"):
            rows.append(("smiles", line.split("\t", 1)[1]))
        else:
            rows.append(("smiles", line))
    return rows


@pytest.fixture(scope="session")
def original_pattern_reference():
    profile_path = ROOT / "tools/testdata/rdkit/pattern_fingerprint_profile.json"
    raw = profile_path.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == "52b748e5397ea81ff77ead0b7c7c9cf42bdfe076410821e304a9c4f0a84a197a"
    branch_profiles = {b["name"]: b for b in json.loads(raw)["branches"]}
    assert len(branch_profiles) == 11
    references = {}
    for profile, (count, variable, sha, relative_input) in PROFILES.items():
        raw = Path(os.environ[variable]).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == sha
        records = [json.loads(line) for line in raw.splitlines()]
        rows = input_rows(ROOT / relative_input)
        assert len(records) == len(rows) == count
        for row, (record, (kind, text)) in enumerate(zip(records, rows, strict=True)):
            assert (record["row"], record["input_kind"], record["smiles"]) == (row, kind, text)
        references[profile] = records
    with Path(os.environ["COSMOLKIT_PATTERN_NATIVE_OBSERVATIONS"]).open("x") as observations:
        yield references, branch_profiles, observations


@pytest.mark.parametrize("profile,row", [(p, i) for p, (n, *_rest) in PROFILES.items() for i in range(n)])
def test_original_pattern_all_rows_all11_complete_results(profile, row, original_pattern_reference):
    references, profiles, observations = original_pattern_reference
    expected = references[profile][row]
    failures = []
    observed = {"profile": profile, "row": row, "input_kind": expected["input_kind"], "smiles": expected["smiles"], "branches": {}, "source_error": expected["error"]}
    try:
        value = (ck.Molecule.from_smiles(expected["smiles"]) if expected["input_kind"] == "smiles" else ck.parse_smarts(expected["smiles"]))
    except Exception as error:
        observed.update(ok=False, error={"type": type(error).__name__, "message": str(error), "domain": getattr(error, "domain", None), "kind": getattr(error, "kind", None)})
        if expected["rdkit_ok"]:
            failures.append(("parse", observed["error"]))
        else:
            assert not expected["branches"] and expected["error"] is not None
    else:
        observed.update(ok=True, error=None)
        assert expected["rdkit_ok"] and expected["error"] is None
        assert set(expected["branches"]) == set(profiles)
        before = value.to_smiles() if expected["input_kind"] == "smiles" else ck.write_smarts(value, ck.SmartsWriteParams())
        for name, branch in expected["branches"].items():
            p = branch["parameters"]
            assert p == profiles[name]
            if p["atomCounts"] == "none":
                assert branch["atom_counts_before"] is branch["atom_counts_after"] is None
            else:
                assert branch["atom_counts_before"] == branch["atom_counts_after"]
            if p["setOnlyBits"] == "wrong_width":
                assert not branch["ok"] and "bad setOnlyBits size" in branch["error"]
                observed["branches"][name] = {"parameters": p, "public_arguments_omitted_as_source_inert": True, "source_validation": branch}
                continue
            assert branch["ok"] and branch["error"] is None
            assert (branch["set_only_bits"] is None) == (p["setOnlyBits"] == "none")
            params = ck.PatternFingerprintParams(n_bits=p["fpSize"], tautomeric=p["tautomericFingerprint"])
            result = {"parameters": p, "atom_counts_before": branch["atom_counts_before"], "atom_counts_after": branch["atom_counts_after"], "set_only_bits": branch["set_only_bits"], "source_inert_arguments_omitted": True}
            try:
                fp = (value.pattern_fingerprint_with_params(params) if expected["input_kind"] == "smiles" else ck.pattern_query_fingerprint_with_params(value, params))
                result.update(ok=True, n_bits=fp.n_bits(), on_bits=fp.on_bits(), error=None)
                assert (params.n_bits, params.tautomeric) == (p["fpSize"], p["tautomericFingerprint"])
                if name == "default":
                    default = value.pattern_fingerprint() if expected["input_kind"] == "smiles" else ck.pattern_query_fingerprint(value)
                    assert (default.n_bits(), default.on_bits()) == (fp.n_bits(), fp.on_bits())
            except Exception as error:
                result.update(ok=False, n_bits=None, on_bits=None, error={"type": type(error).__name__, "message": str(error), "domain": getattr(error, "domain", None), "kind": getattr(error, "kind", None)})
            observed["branches"][name] = result
            for field in ("ok", "n_bits", "on_bits"):
                if result[field] != branch[field]:
                    failures.append((name, field, result[field], branch[field]))
        after = value.to_smiles() if expected["input_kind"] == "smiles" else ck.write_smarts(value, ck.SmartsWriteParams())
        if after != before:
            failures.append(("input changed", before, after))
    observations.write(json.dumps(observed, sort_keys=True, separators=(",", ":")) + "\n")
    observations.flush()
    assert not failures, (profile, row, expected["smiles"], failures)
