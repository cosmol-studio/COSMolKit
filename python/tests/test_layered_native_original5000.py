"""All original Layered18 configurations and all original5000 inputs, without filters."""
import hashlib
import json
import os
from pathlib import Path

import cosmolkit as ck
import pytest

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "testdata/smiles/corpus/smiles_5000.smi"
PROFILE = ROOT / "tools/testdata/rdkit/layered_fingerprint_profile.json"
GOLDEN_SHA = "68a2029b115c7b962971317c3ee11af21626d71ce5b98327601bc30190377482"
CORPUS_SHA = "a4d579cd72621af27772256bb23ba796452276bb924fd20aac83625ffa67d849"


@pytest.fixture(scope="session")
def original_layered_reference():
    golden = Path(os.environ["COSMOLKIT_LAYERED_ORIGINAL5000_GOLDEN"])
    raw = golden.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == GOLDEN_SHA
    assert hashlib.sha256(CORPUS.read_bytes()).hexdigest() == CORPUS_SHA
    rows = [line.strip() for line in CORPUS.read_text().splitlines()
            if line.strip() and not line.strip().startswith("#")]
    records = [json.loads(line) for line in raw.splitlines()]
    profiles = {p["name"]: p for p in json.loads(PROFILE.read_text())["branches"]}
    assert len(rows) == len(records) == 5000
    assert len(profiles) == 18
    with Path(os.environ["COSMOLKIT_LAYERED_ORIGINAL5000_OBSERVATIONS"]).open("x") as observations:
        yield records, rows, profiles, observations


@pytest.mark.parametrize("row", range(5000))
def test_original5000_all_layered18_complete_results(row, original_layered_reference):
    records, corpus, profiles, observations = original_layered_reference
    expected = records[row]
    assert expected["row"] == row
    assert expected["smiles"] == corpus[row]
    assert set(expected["branches"]) == set(profiles)
    failures = []
    observed = {"row": row, "smiles": corpus[row], "ok": False, "branches": {}}
    try:
        molecule = ck.Molecule.from_smiles(corpus[row])
    except Exception as error:
        observed["error"] = {"type": type(error).__name__, "message": str(error)}
        failures.append(("parse", observed["error"]))
    else:
        observed["ok"] = True
        observed["error"] = None
        before = molecule.to_smiles()
        assert expected["rdkit_ok"] and expected["error"] is None
        for name, branch in expected["branches"].items():
            assert branch["parameters"] == profiles[name]
            p = branch["parameters"]
            a = branch["resolved_arguments"]
            assert a["layerFlags"] == int(p["layerFlags"], 16)
            for key in ("fromAtoms", "atomCounts", "setOnlyBits"):
                assert (key in p) == (a[key] is not None)
            params = ck.LayeredFingerprintParams(
                layers=a["layerFlags"], min_path=p["minPath"], max_path=p["maxPath"],
                fp_size=p["fpSize"], branched_paths=p["branchedPaths"],
                from_atoms=a["fromAtoms"], atom_counts=a["atomCounts"],
                set_only_bits=(None if a["setOnlyBits"] is None else
                               ck.Fingerprint.from_on_bits(p["fpSize"], a["setOnlyBits"])),
            )
            result = {"parameters": p, "resolved_arguments": a}
            try:
                output = molecule.layered_fingerprint_with_output_with_params(params)
                fp = output.fingerprint()
                result.update(ok=True, num_bits=fp.n_bits(), on_bits=fp.on_bits(),
                              atom_counts=output.atom_counts(), error=None)
                # Each original scalar bits entry must agree with the complete output entry.
                bits = molecule.layered_fingerprint_with_params(params)
                assert (bits.n_bits(), bits.on_bits()) == (result["num_bits"], result["on_bits"])
                assert params.atom_counts == a["atomCounts"]
                assert params.from_atoms == a["fromAtoms"]
                mask = params.set_only_bits
                assert (None if mask is None else mask.on_bits()) == a["setOnlyBits"]
            except Exception as error:
                result.update(ok=False, num_bits=None, on_bits=None, atom_counts=None,
                              error={"type": type(error).__name__, "message": str(error)})
            observed["branches"][name] = result
            for field in ("ok", "num_bits", "on_bits", "atom_counts"):
                if result[field] != branch[field]:
                    failures.append((name, field, result[field], branch[field]))
            if branch["ok"]:
                assert branch["error"] is None
        if molecule.to_smiles() != before:
            failures.append(("molecule changed", before, molecule.to_smiles()))
    observations.write(json.dumps(observed, sort_keys=True, separators=(",", ":")) + "\n")
    observations.flush()
    assert not failures, (row, corpus[row], failures)
