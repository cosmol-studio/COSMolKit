"""All ten original native legacy torsion records on all unchanged 5000 rows.

Author test proposal; p1 reviews conditions and ROOT adjudicates acceptance.
"""
import json
from pathlib import Path
import pytest
import cosmolkit as ck
from test_canonical_topological_torsion_reusable_native_5000 import CORPUS, GOLDEN, sha

ROOT = Path(__file__).resolve().parents[2]

@pytest.fixture(scope="module")
def original_legacy():
    assert sha(CORPUS) == "a4d579cd72621af27772256bb23ba796452276bb924fd20aac83625ffa67d849"
    assert sha(GOLDEN) == "28d91a8ca475ffedf43ac0e50621be2b0b6ae394225dffbe76b6fff061affe7d"
    smiles = [line.strip() for line in CORPUS.read_text().splitlines() if line.strip() and not line.strip().startswith("#")]
    assert len(smiles) == 5000
    positions = []
    with GOLDEN.open("rb") as stream:
        while True:
            offset = stream.tell()
            if not stream.readline(): break
            positions.append(offset)
    assert len(positions) == 5000
    actual = (ROOT / "target/agent-handoff/fp-graph-block/topological-torsion-native-original5000-legacy-public-observations-v1.jsonl").open("x")
    try: yield smiles, positions, actual
    finally: actual.close()

@pytest.mark.parametrize("index", range(5000))
def test_original_native_legacy_all_ten_calls(index, original_legacy):
    smiles, positions, actual = original_legacy
    with GOLDEN.open("rb") as stream:
        stream.seek(positions[index]); expected = json.loads(stream.readline())
    assert expected["smiles"] == smiles[index] and expected["rdkit_ok"]
    observed = {"index": index, "smiles": smiles[index], "legacy": {}, "calls": 0, "differences": []}
    try: molecule = ck.Molecule.from_smiles(smiles[index])
    except Exception as error:
        observed["parse_error"] = {"type": type(error).__name__, "message": str(error)}
        actual.write(json.dumps(observed) + "\n"); actual.flush()
        pytest.fail(f"Original successful row failed to parse: {observed}")
    calls = [
        ("unfolded", "legacy_topological_torsion_sparse_count_fingerprint", None),
        ("unfolded_chiral", "legacy_topological_torsion_sparse_count_fingerprint_with_params", ck.LegacyTopologicalTorsionParams(include_chirality=True)),
        ("unfolded_custom", "legacy_topological_torsion_sparse_count_fingerprint_with_params", ck.LegacyTopologicalTorsionParams(custom_atom_invariants=[i + 17 for i in range(molecule.num_atoms())])),
        ("hashed", "legacy_topological_torsion_count_fingerprint_with_params", ck.LegacyTopologicalTorsionParams(fp_size=1000)),
        ("hashed_rooted", "legacy_topological_torsion_count_fingerprint_with_params", ck.LegacyTopologicalTorsionParams(fp_size=1000, from_atoms=[0] if molecule.num_atoms() else [])),
        ("hashed_ignored", "legacy_topological_torsion_count_fingerprint_with_params", ck.LegacyTopologicalTorsionParams(fp_size=1000, ignore_atoms=[0] if molecule.num_atoms() else [])),
    ] + [("bits_" + str(entry), "legacy_topological_torsion_fingerprint_with_params", ck.LegacyTopologicalTorsionParams(fp_size=256, bits_per_entry=entry)) for entry in [1, 2, 4, 6]]
    for name, method, params in calls:
        observed["calls"] += 1
        target = expected["legacy"]["bit_vectors"][name[5:]] if name.startswith("bits_") else expected["legacy"][name]
        try:
            value = getattr(molecule, method)() if params is None else getattr(molecule, method)(params)
            result = {"size": value.n_bits(), "on_bits": value.on_bits()} if name.startswith("bits_") else {"size": value.length(), "nonzero_elements": {str(k): v for k, v in value.nonzero_elements().items()}}
            observed["legacy"][name] = result
            if result != target: observed["differences"].append(name)
        except Exception as error:
            observed["legacy"][name] = {"error": {"type": type(error).__name__, "message": str(error), "domain": getattr(error, "domain", None), "kind": getattr(error, "kind", None), "cause": str(error.__cause__)}}
            observed["differences"].append(name + "/exception")
    actual.write(json.dumps(observed, ensure_ascii=True) + "\n"); actual.flush()
    assert observed["calls"] == 10
    assert not observed["differences"], f"row {index} {smiles[index]}: {observed['differences']}; full actual bytes retained"
