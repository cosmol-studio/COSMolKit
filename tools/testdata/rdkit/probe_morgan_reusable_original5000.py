"""Pinned native reusable Morgan supplement; preserves original 17 profiles.

Imports the frozen preparation generator only, never ordinary test modules.
Bulk source API has no per-call roots/custom invariants; its independent values
are retained separately from original scalar golden rows.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator as fg

FORMS = [("sparse_count", "GetSparseCountFingerprints"),
         ("sparse_bit", "GetSparseFingerprints"),
         ("hashed_count", "GetCountFingerprints"),
         ("explicit_bit", "GetFingerprints")]

def record(value):
    if value is None:
        return None
    if hasattr(value, "GetNonzeroElements"):
        return {"size": int(value.GetLength()), "nonzero": sorted([int(k), int(v)] for k, v in value.GetNonzeroElements().items())}
    return {"size": int(value.GetNumBits()), "native_bits": [int(bit) for bit in value.GetOnBits()]}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--source-generator", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert hashlib.sha256(args.corpus.read_bytes()).hexdigest() == "a4d579cd72621af27772256bb23ba796452276bb924fd20aac83625ffa67d849"
    assert hashlib.sha256(args.source_generator.read_bytes()).hexdigest() == "925b88bec3db4a2f47f7d68acc4bab67d8796049fdd021e0061cf066f10eb94a"
    spec = importlib.util.spec_from_file_location("frozen_source_generator", args.source_generator)
    original = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(original)
    assert len(original.BRANCHES) == 17
    smiles = list(original.iter_smiles(args.corpus))
    assert len(smiles) == 5000
    Chem.SetUseLegacyStereoPerception(True)
    molecules = [Chem.MolFromSmiles(text) for text in smiles]
    assert all(molecule is not None for molecule in molecules)
    rows = [None, *molecules[:2500], None, *molecules[2500:], None]
    assert len(rows) == 5003
    args.output.mkdir()
    metadata = []
    total_bulk = total_json = 0
    for branch in original.BRANCHES:
        name = branch["name"]
        generator = original.make_generator(branch)
        payload = generator.ToJSON()
        restored = fg.FingerprintGeneratorFromJSON(payload)
        bulk = {form: getattr(generator, method)(rows, numThreads=2) for form, method in FORMS}
        for values in bulk.values():
            assert len(values) == 5003 and values[0] is None and values[2501] is None and values[-1] is None
        metadata.append({"name": name, "source_branch": branch, "info": generator.GetInfoString(), "json": json.loads(payload),
                         "restored_info": restored.GetInfoString(), "restored_json": json.loads(restored.ToJSON()),
                         "bulk_rows": 5003, "none_indices": [0, 2501, 5002]})
        with (args.output / (name + ".jsonl")).open("x") as stream:
            for index, (text, molecule) in enumerate(zip(smiles, molecules)):
                position = index + 1 + (index >= 2500)
                call = original.fingerprint_kwargs(molecule, branch)
                value = {"index": index, "smiles": text, "bulk": {form: record(values[position]) for form, values in bulk.items()},
                         "restored_sparse_count": record(restored.GetSparseCountFingerprint(molecule, **call))}
                total_bulk += 4
                total_json += 1
                stream.write(json.dumps(value, sort_keys=True) + "\n")
        print(json.dumps({"profile": name, "bulk_values": 20000, "json_calls": 5000}), flush=True)
    with (args.output / "metadata.json").open("x") as stream:
        json.dump(metadata, stream, indent=2, sort_keys=True)
    print(json.dumps({"profiles": len(metadata), "inputs": len(smiles), "bulk_values": total_bulk,
                      "source_none_slots": 204, "json_calls": total_json, "test_count": 0}), flush=True)

if __name__ == "__main__":
    main()
