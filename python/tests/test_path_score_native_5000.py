"""Full original5000, immutable source-prepared path values."""
import hashlib
import json
import os
from pathlib import Path
import cosmolkit as ck
import pytest

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "testdata/smiles/corpus/smiles_5000.smi"
GOLDEN = Path(os.environ["COSMOLKIT_PATH_SOURCE_GOLDEN"])
assert hashlib.sha256(CORPUS.read_bytes()).hexdigest() == "a4d579cd72621af27772256bb23ba796452276bb924fd20aac83625ffa67d849"
assert hashlib.sha256(GOLDEN.read_bytes()).hexdigest() == "f919a277cce6741dec5b792b1dacc6d17cb90deef70bf0ef88158449ab26da03"
SMILES = [s.strip() for s in CORPUS.read_text().splitlines() if s.strip() and not s.strip().startswith("#")]
ROWS = [json.loads(s) for s in GOLDEN.read_text().splitlines()]
assert len(SMILES) == len(ROWS) == 5000
ACTUAL = Path(os.environ["COSMOLKIT_PATH_ACTUAL_OUTPUT"]).open("x")

@pytest.mark.parametrize("index", range(5000))
def test_all_original_inputs_all_four_profiles_score_and_explanation(index):
    expected = ROWS[index]
    assert expected["index"] == index and expected["smiles"] == SMILES[index]
    observed = {"index": index, "smiles": SMILES[index], "profiles": {}}
    differences = []
    try:
        molecule = ck.Molecule.from_smiles(SMILES[index])
        observed["atoms"] = molecule.num_atoms()
        assert observed["atoms"] == expected["atoms"]
        for name, source in expected["profiles"].items():
            value = {k: source[k] for k in ("path", "size", "atom_codes")}
            try:
                score = molecule.topological_torsion_path_score(source["path"], source["size"], source["atom_codes"])
                value["score"] = score
                explanation = ck.explain_path_score(score, source["size"])
                assert type(explanation) is tuple and all(type(row) is tuple for row in explanation)
                value["explanation"] = [list(row) for row in explanation]
                value["error"] = None
            except Exception as error:
                value["error"] = {"type": type(error).__name__, "message": str(error), "domain": getattr(error, "domain", None), "kind": getattr(error, "kind", None)}
            observed["profiles"][name] = value
            if value != source:
                differences.append({"profile": name, "source": source, "actual": value})
    except Exception as error:
        observed["parse_or_transport_error"] = {"type": type(error).__name__, "message": str(error)}
        differences.append(observed["parse_or_transport_error"])
    finally:
        ACTUAL.write(json.dumps(observed, sort_keys=True) + "\n")
        ACTUAL.flush()
    assert not differences, differences
