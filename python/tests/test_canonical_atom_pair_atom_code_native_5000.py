"""Complete immutable native default and legacy chiral atom-code helper rows."""
import json
from pathlib import Path
import pytest
import cosmolkit as ck
from test_canonical_topological_torsion_reusable_native_5000 import CORPUS, GOLDEN, sha
ROOT = Path(__file__).resolve().parents[2]

@pytest.fixture(scope="module")
def original_atom_codes():
    assert sha(CORPUS) == "a4d579cd72621af27772256bb23ba796452276bb924fd20aac83625ffa67d849"
    assert sha(GOLDEN) == "28d91a8ca475ffedf43ac0e50621be2b0b6ae394225dffbe76b6fff061affe7d"
    smiles = [s.strip() for s in CORPUS.read_text().splitlines() if s.strip() and not s.strip().startswith("#")]
    assert len(smiles) == 5000
    positions = []
    with GOLDEN.open("rb") as stream:
        while True:
            offset = stream.tell()
            if not stream.readline(): break
            positions.append(offset)
    assert len(positions) == 5000
    actual = (ROOT / "target/agent-handoff/fp-graph-block/atom-pair-atom-code-original5000-public-observations-v1.jsonl").open("x")
    try: yield smiles, positions, actual
    finally: actual.close()

@pytest.mark.parametrize("index", range(5000))
def test_original_every_atom_default_and_chiral_source_code(index, original_atom_codes):
    smiles, positions, actual = original_atom_codes
    with GOLDEN.open("rb") as stream:
        stream.seek(positions[index]); expected = json.loads(stream.readline())
    assert expected["smiles"] == smiles[index] and expected["rdkit_ok"]
    observed = {"index": index, "smiles": smiles[index], "atoms": []}
    try:
        molecule = ck.Molecule.from_smiles(smiles[index])
        for atom_index in range(len(expected["helpers"]["atom_codes"])):
            result = {"atom": atom_index}
            try:
                result["code"] = molecule.with_atom_pair_atom_code(atom_index).code
                result["chiral_code"] = molecule.with_atom_pair_atom_code(atom_index, include_chirality=True).code
            except Exception as error:
                result["error"] = {"type": type(error).__name__, "message": str(error), "domain": getattr(error, "domain", None), "kind": getattr(error, "kind", None), "cause": str(error.__cause__)}
            observed["atoms"].append(result)
    except Exception as error:
        observed["parse_error"] = {"type": type(error).__name__, "message": str(error)}
    actual.write(json.dumps(observed, ensure_ascii=True) + "\n"); actual.flush()
    assert "parse_error" not in observed, observed
    assert len(observed["atoms"]) == len(expected["helpers"]["atom_codes"])
    for result, atom in zip(observed["atoms"], expected["helpers"]["atom_codes"], strict=True):
        assert "error" not in result, f"original successful row {index}: {result}"
        assert result["code"] == atom["code"]
        assert result["chiral_code"] == atom["chiral_code"]
