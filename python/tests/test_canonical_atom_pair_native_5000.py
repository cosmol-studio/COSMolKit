"""Complete original 5000-row/10-branch/four-form native AtomPair delivery proposal.

Reference bytes come from the pinned source-built private RDKit service.
All original inputs/branches/outputs are preserved; independent review pending.
"""
import ast
import hashlib
import json
import os
from pathlib import Path
import pytest
import cosmolkit as ck

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "testdata/smiles/corpus/smiles_5000.smi"
REFERENCE_ROOT = Path(os.environ["COSMOLKIT_FP_SOURCE_REFERENCE_ROOT"]).resolve()
GOLDEN = REFERENCE_ROOT / "target/agent-handoff/structure-oracle-native-env/generated-source-native-v1/atom_pair_fingerprint.jsonl"
GENERATOR = ROOT / "tools/testdata/rdkit/_generate_atom_pair_fingerprint_golden.py"
assert hashlib.sha256(CORPUS.read_bytes()).hexdigest() == "a4d579cd72621af27772256bb23ba796452276bb924fd20aac83625ffa67d849"
assert hashlib.sha256(GOLDEN.read_bytes()).hexdigest() == "5ad83a60e5109fe21510390f2a7fb21a4eaa72de9aa911c35dcd7a1a3663765a"
assert hashlib.sha256(GENERATOR.read_bytes()).hexdigest() == "6bfc31870e263d6e16f4b5a1d36c0559263e92c98f62e6c7d84fe5014856d41c"
PROFILE = ROOT / "tools/testdata/rdkit/atom_pair_fingerprint_profile.json"
assert hashlib.sha256(PROFILE.read_bytes()).hexdigest() == "b104ba22a87409b3631f1dd10c0bb5aec775115d8d2c3203ad2d37e9d504d261"
BRANCHES = json.loads(PROFILE.read_bytes())["branches"]
assert len(BRANCHES) == 10
SMILES = [line.strip() for line in CORPUS.read_text().splitlines() if line.strip() and not line.strip().startswith("#")]
assert len(SMILES) == 5000
POSITIONS = []
with GOLDEN.open("rb") as handle:
    while True:
        position = handle.tell()
        if not handle.readline(): break
        POSITIONS.append(position)
assert len(POSITIONS) == 5000
ACTUAL = Path(os.environ["COSMOLKIT_FP_ACTUAL_OUTPUT"]).open("x")
METHODS = [("sparse_count", "atom_pair_sparse_count_fingerprint_with_params"), ("sparse_bit", "atom_pair_sparse_fingerprint_with_params"), ("count", "atom_pair_count_fingerprint_with_params"), ("explicit_bit", "atom_pair_fingerprint_with_params")]


def params_for(molecule, branch):
    atoms = molecule.num_atoms()
    return ck.AtomPairFingerprintParams(generator=ck.AtomPairParams(min_distance=branch["minDistance"], max_distance=branch["maxDistance"], include_chirality=branch["includeChirality"], use_2d=branch["use2D"], count_simulation=branch["countSimulation"], fp_size=branch["fpSize"], count_bounds=branch["countBounds"], bits_per_feature=branch["numBitsPerFeature"]), from_atoms=[0] if branch.get("fromAtoms")=="first" and atoms else None, ignore_atoms=[0] if branch.get("ignoreAtoms")=="first" and atoms else None, custom_atom_invariants=[i+11 for i in range(atoms)] if branch.get("customAtomInvariants")=="index_plus_11" else None)


def metadata_record(output):
    return {"atom_counts":output.atom_counts(), "atom_to_bits":output.atom_to_bits(), "bit_info_map":{str(key):[list(pair) for pair in pairs] for key,pairs in output.bit_info_map().items()}, "atoms_per_bit":{str(key):[list(atoms) for atoms in groups] for key,groups in output.atoms_per_bit().items()}}


@pytest.mark.parametrize("index", range(5000))
def test_current_atom_pair_matches_original_native_all_branches_and_outputs(index):
    with GOLDEN.open("rb") as handle:
        handle.seek(POSITIONS[index]); expected = json.loads(handle.readline())
    assert expected["smiles"] == SMILES[index]
    errors = []
    observed = {"index":index,"smiles":SMILES[index],"branches":{},"calls":0}
    try:
        molecule = ck.Molecule.from_smiles(SMILES[index])
    except Exception as error:
        observed["parse_error"]={"type":type(error).__name__,"message":str(error)}
        ACTUAL.write(json.dumps(observed)+"\n"); ACTUAL.flush()
        pytest.fail(str(observed["parse_error"]))
    for branch in BRANCHES:
        name = branch["name"]; observed["branches"][name] = {}
        params = params_for(molecule, branch)
        for output_name, method in METHODS:
            target = expected["branches"][name][output_name]
            output = None
            if branch.get("additionalOutput"):
                output = ck.FingerprintAdditionalOutput(); output.allocate_atom_counts(); output.allocate_atom_to_bits(); output.allocate_bit_info_map(); output.allocate_atoms_per_bit()
            observed["calls"] += 1
            try:
                fingerprint = getattr(molecule, method)(params, output)
                value = {"ok":True,"error":None,"length":fingerprint.length() if output_name in ("sparse_count","count") else fingerprint.n_bits()}
                if output_name in ("sparse_count","count"):
                    value["nonzero_elements"]={str(bit):count for bit,count in fingerprint.nonzero_elements().items()}
                else:
                    bits=sorted(x & 0xFFFFFFFF for x in fingerprint.on_bits()); value.update(on_bits=bits)
                if output is not None: value["additional_output"] = metadata_record(output)
            except Exception as error:
                value={"ok":False,"type":type(error).__name__,"error":str(error),"domain":getattr(error,"domain",None),"kind":getattr(error,"kind",None),"cause":str(error.__cause__),"cause_type":type(error.__cause__).__name__}
            observed["branches"][name][output_name] = value
            if target["ok"]:
                matched = value == target
            else:
                # New proposal mapping: preserve the original IndexError and
                # exact index, then check the canonical typed error/cause.
                # Original expected bytes remain unchanged; p1 review pending.
                assert target["error"].startswith("IndexError: ")
                bit = int(target["error"].split(": ",1)[1])
                expected_cause = f"fingerprint index {bit} is outside vector length {1 << (27 if branch['includeChirality'] else 23)}"
                matched = (not value["ok"] and value["type"] == "AtomPairReadError" and value["domain"] == "fingerprints" and value["kind"] == "Generator" and value["cause"] == expected_cause and value["cause_type"] == "ValueError" and value["error"] == "AtomPair generation failed: " + expected_cause)
            if not matched: errors.append(f"{name}/{output_name}")
    ACTUAL.write(json.dumps(observed,ensure_ascii=True)+"\n"); ACTUAL.flush()
    assert observed["calls"] == 40
    assert not errors, f"row {index} {SMILES[index]} differs: {errors} (full actual bytes retained in observations JSONL)"
