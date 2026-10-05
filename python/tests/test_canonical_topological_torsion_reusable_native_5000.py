"""Unmodified pinned native rows, all nine profiles and reusable public calls.

Test/comparison additions are delivery proposals for independent p1 review.
Legacy and helper records stay in the original golden; this test covers modern
single calls, all provenance fields, bulk, metadata and restored JSON calls.
"""
import hashlib
import json
from pathlib import Path

import cosmolkit as ck
import pytest

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "testdata/smiles/corpus/smiles_5000.smi"
GOLDEN = ROOT / "target/agent-handoff/structure-oracle-native-env/generated-source-native-v1/topological_torsion_fingerprint.jsonl"
PROFILE = ROOT / "tools/testdata/rdkit/topological_torsion_fingerprint_profile.json"
SOURCE_GENERATOR = ROOT / "tools/testdata/rdkit/_generate_topological_torsion_fingerprint_golden.py"
METHODS = [
    ("sparse_count", "topological_torsion_sparse_count_fingerprint_with_generator", "sparse_counts"),
    ("sparse_bit", "topological_torsion_sparse_fingerprint_with_generator", "sparse_fingerprints"),
    ("count", "topological_torsion_count_fingerprint_with_generator", "counts"),
    ("bit", "topological_torsion_fingerprint_with_generator", "fingerprints"),
]


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            digest.update(block)
    return digest.hexdigest()


def record(fp, form):
    if form in ("count", "sparse_count"):
        return {"size": fp.length(), "nonzero_elements": {str(k): v for k, v in fp.nonzero_elements().items()}}
    # Only the source SparseBitVect signed-index projection is u32. AO is u64.
    bits = sorted(bit & 0xFFFFFFFF for bit in fp.on_bits()) if form == "sparse_bit" else fp.on_bits()
    return {"size": fp.n_bits(), "on_bits": bits}


def metadata(output):
    return {
        "atom_counts": output.atom_counts(),
        "atom_to_bits": output.atom_to_bits(),
        "bit_info_map": {str(k): [list(pair) for pair in v] for k, v in output.bit_info_map().items()},
        "bit_paths": {str(k): [list(path) for path in v] for k, v in output.bit_paths().items()},
        "atoms_per_bit": {str(k): [list(atoms) for atoms in v] for k, v in output.atoms_per_bit().items()},
    }


def error_record(error):
    return {"type": type(error).__name__, "message": str(error), "domain": getattr(error, "domain", None), "kind": getattr(error, "kind", None), "cause": str(error.__cause__)}


@pytest.fixture(scope="module")
def original_and_reusable():
    for path, expected in [
        (CORPUS, "a4d579cd72621af27772256bb23ba796452276bb924fd20aac83625ffa67d849"),
        (GOLDEN, "28d91a8ca475ffedf43ac0e50621be2b0b6ae394225dffbe76b6fff061affe7d"),
        (PROFILE, "f03a8b5a1e9b9f05347c0058761f7cac2efed4ea0baafce4bbcb967e44040daa"),
        (SOURCE_GENERATOR, "820ca13d717f691890bbc3e80cc803096df3a1515358bf83ecf47c89b6cb0ee4"),
    ]:
        assert sha(path) == expected
    smiles = [line.strip() for line in CORPUS.read_text().splitlines() if line.strip() and not line.strip().startswith("#")]
    assert len(smiles) == 5000
    molecules = []
    parse_errors = {}
    for index, text in enumerate(smiles):
        try:
            molecules.append(ck.Molecule.from_smiles(text))
        except Exception as error:
            molecules.append(None)
            parse_errors[index] = error_record(error)
    branches = json.loads(PROFILE.read_bytes())["corpus_branches"]
    assert len(branches) == 9
    generators, restored, bulk, infos, payloads = {}, {}, {}, {}, {}
    for branch in branches:
        name = branch["name"]
        generator = ck.TopologicalTorsionFingerprintGenerator(params=ck.TopologicalTorsionParams(
            torsion_atom_count=branch.get("torsionAtomCount", 4),
            include_chirality=branch.get("includeChirality", False),
            count_simulation=branch.get("countSimulation", True),
            count_bounds=branch.get("countBounds"), fp_size=branch.get("fpSize", 2048),
        ))
        settings = generator.settings()
        settings.only_shortest_paths = branch.get("onlyShortestPaths", False)
        settings.bits_per_feature = branch.get("numBitsPerFeature", 1)
        generators[name] = generator
        infos[name] = generator.info_string()
        payloads[name] = json.loads(generator.to_json())
        restored[name] = ck.TopologicalTorsionFingerprintGenerator.from_json(generator.to_json())
        bulk[name] = {}
        for form, _, method in METHODS:
            try:
                values = getattr(generator, method)(molecules, num_threads=2)
                assert len(values) == 5000
                bulk[name][form] = [None if value is None else record(value, form) for value in values]
            except Exception as error:
                bulk[name][form] = {"error": error_record(error)}
    positions = []
    with GOLDEN.open("rb") as stream:
        while True:
            offset = stream.tell()
            if not stream.readline():
                break
            positions.append(offset)
    assert len(positions) == 5000
    actual = (ROOT / "target/agent-handoff/fp-graph-block/topological-torsion-native-5000-reusable-state-bulk-json-observations-v1.jsonl").open("x")
    try:
        yield smiles, molecules, parse_errors, branches, generators, restored, bulk, infos, payloads, positions, actual
    finally:
        actual.close()


@pytest.mark.parametrize("index", range(5000))
def test_reusable_modern_full_native_row(index, original_and_reusable):
    smiles, molecules, parse_errors, branches, generators, restored, bulk, infos, payloads, positions, actual = original_and_reusable
    with GOLDEN.open("rb") as stream:
        stream.seek(positions[index])
        expected = json.loads(stream.readline())
    assert expected["smiles"] == smiles[index] and expected["rdkit_ok"]
    observed = {"index": index, "smiles": smiles[index], "profiles": {}, "single_calls": 0, "bulk_values": 0, "json_calls": 0}
    differences = []
    if index in parse_errors:
        observed["parse_error"] = parse_errors[index]
        actual.write(json.dumps(observed) + "\n")
        actual.flush()
        pytest.fail(f"original successful input failed to parse: {observed}")
    molecule = molecules[index]
    for branch in branches:
        name = branch["name"]
        target = expected["profiles"][name]
        generator = generators[name]
        call = ck.TopologicalTorsionCallParams(
            from_atoms=[0] if branch.get("fromAtoms") == "first" and molecule.num_atoms() else None,
            ignore_atoms=[0] if branch.get("ignoreAtoms") == "first" and molecule.num_atoms() else None,
            custom_atom_invariants=[i + 17 for i in range(molecule.num_atoms())] if branch.get("customAtomInvariants") == "index_plus_17" else None,
        )
        result = observed["profiles"][name] = {"info_string": infos[name], "json": payloads[name], "bulk": {}}
        if infos[name] != target["info_string"]:
            differences.append(f"{name}/info_string")
        if payloads[name] != target["json"]:
            differences.append(f"{name}/json")
        for form, method, _ in METHODS:
            output = None
            if branch.get("additionalOutput"):
                output = ck.FingerprintAdditionalOutput()
                output.allocate_atom_counts(); output.allocate_atom_to_bits(); output.allocate_bit_info_map(); output.allocate_bit_paths(); output.allocate_atoms_per_bit()
            observed["single_calls"] += 1
            try:
                value = record(getattr(molecule, method)(generator, params=call, output=output), form)
                ao = None if output is None else metadata(output)
                result[form] = {"value": value, "additional_output": ao}
                if value != target[form]:
                    differences.append(f"{name}/{form}/value")
                if ao is not None and ao != target["additional_output"][form]:
                    differences.append(f"{name}/{form}/AO")
            except Exception as error:
                result[form] = {"error": error_record(error)}
                differences.append(f"{name}/{form}/exception")
            values = bulk[name][form]
            result["bulk"][form] = values[index] if isinstance(values, list) else values
            observed["bulk_values"] += 1
            if result["bulk"][form] != target["bulk"][form]:
                differences.append(f"{name}/{form}/bulk")
        observed["json_calls"] += 1
        try:
            result["json_restored_count"] = record(molecule.topological_torsion_count_fingerprint_with_generator(restored[name], params=call), "count")
            if result["json_restored_count"] != target["json_restored_count"]:
                differences.append(f"{name}/json_restored_count")
        except Exception as error:
            result["json_restored_count"] = {"error": error_record(error)}
            differences.append(f"{name}/json_restored_count/exception")
    observed["differences"] = differences
    actual.write(json.dumps(observed, ensure_ascii=True) + "\n")
    actual.flush()
    assert observed["single_calls"] == 36 and observed["bulk_values"] == 36 and observed["json_calls"] == 9
    assert not differences, f"row {index} {smiles[index]} differs: {differences}; full actual bytes retained"
