"""All original native helper atom codes/explanations, with no case filtering."""
import json
from pathlib import Path
import pytest
import cosmolkit as ck
from test_canonical_atom_code_explanation import fields
from test_canonical_topological_torsion_reusable_native_5000 import CORPUS, GOLDEN, sha
ROOT = Path(__file__).resolve().parents[2]

@pytest.fixture(scope="module")
def original_codes():
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
    actual = (ROOT / "target/agent-handoff/fp-graph-block/atom-code-explanation-original5000-public-observations-v1.jsonl").open("x")
    try: yield smiles, positions, actual
    finally: actual.close()

@pytest.mark.parametrize("index", range(5000))
def test_original_all_atom_explanations(index, original_codes):
    smiles, positions, actual = original_codes
    with GOLDEN.open("rb") as stream:
        stream.seek(positions[index]); expected = json.loads(stream.readline())
    assert expected["smiles"] == smiles[index] and expected["rdkit_ok"]
    observed = {"index": index, "smiles": smiles[index], "atoms": []}
    for atom_index, atom in enumerate(expected["helpers"]["atom_codes"]):
        result = {"atom": atom_index}
        try:
            result["explanation"] = fields(ck.AtomCodeExplanation.from_code(atom["code"]))
            result["chiral_explanation"] = fields(ck.AtomCodeExplanation.from_code(atom["chiral_code"], include_chirality=True))
            result["branch_subtract_ignored"] = fields(ck.AtomCodeExplanation.from_code(atom["chiral_code"], 2, True))
        except Exception as error:
            result["error"] = {"type": type(error).__name__, "message": str(error), "domain": getattr(error, "domain", None), "kind": getattr(error, "kind", None)}
        observed["atoms"].append(result)
    actual.write(json.dumps(observed, ensure_ascii=True) + "\n"); actual.flush()
    for result, atom in zip(observed["atoms"], expected["helpers"]["atom_codes"], strict=True):
        assert "error" not in result, f"original successful row {index}: {result}"
        assert result["explanation"] == atom["explanation"]
        assert result["chiral_explanation"] == atom["chiral_explanation"]
        assert result["branch_subtract_ignored"] == atom["chiral_explanation"]
