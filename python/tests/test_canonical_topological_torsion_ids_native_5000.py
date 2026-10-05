"""Source IDs helper on all immutable original native 5000 records."""
import json
from pathlib import Path
import pytest
import cosmolkit as ck
from test_canonical_topological_torsion_reusable_native_5000 import CORPUS, GOLDEN, sha
ROOT = Path(__file__).resolve().parents[2]

@pytest.fixture(scope="module")
def original_ids():
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
    actual = (ROOT / "target/agent-handoff/fp-graph-block/topological-torsion-native-original5000-ids-public-observations-v1.jsonl").open("x")
    try: yield smiles, positions, actual
    finally: actual.close()

@pytest.mark.parametrize("index", range(5000))
def test_original_native_sorted_ids_with_all_multiplicities(index, original_ids):
    smiles, positions, actual = original_ids
    with GOLDEN.open("rb") as stream:
        stream.seek(positions[index]); expected = json.loads(stream.readline())
    assert expected["smiles"] == smiles[index] and expected["rdkit_ok"]
    observed = {"index": index, "smiles": smiles[index]}
    try:
        molecule = ck.Molecule.from_smiles(smiles[index])
        ids = molecule.topological_torsion_ids()
        observed["ids"] = ids
        observed["explicit_default_ids"] = molecule.topological_torsion_ids_with_params(4)
    except Exception as error:
        observed["error"] = {"type": type(error).__name__, "message": str(error), "domain": getattr(error, "domain", None), "kind": getattr(error, "kind", None), "cause": str(error.__cause__)}
    actual.write(json.dumps(observed, ensure_ascii=True) + "\n"); actual.flush()
    assert "error" not in observed, f"original successful row {index}: {observed}"
    assert observed["ids"] == expected["helpers"]["ids"]
    assert observed["explicit_default_ids"] == expected["helpers"]["ids"]
