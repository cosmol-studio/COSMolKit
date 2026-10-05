"""Complete original 5000-row/17-branch/four-form native Morgan delivery proposal.

Reference bytes come from the pinned source-built private RDKit service.
All original inputs/branches/outputs are preserved; independent review pending.
"""
import ast
import hashlib
import json
from pathlib import Path
import pytest
import cosmolkit as ck

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "testdata/smiles/corpus/smiles_5000.smi"
GOLDEN = ROOT / "target/agent-handoff/structure-oracle-native-env/generated-source-native-v1/morgan_fingerprint.jsonl"
GENERATOR = ROOT / "tools/testdata/rdkit/_generate_morgan_fingerprint_golden.py"
assert hashlib.sha256(CORPUS.read_bytes()).hexdigest() == "a4d579cd72621af27772256bb23ba796452276bb924fd20aac83625ffa67d849"
assert hashlib.sha256(GOLDEN.read_bytes()).hexdigest() == "6f2c92c1827fde8ec893b9fddec2d7130d327d199db6cd860facebdfe3406196"
assert hashlib.sha256(GENERATOR.read_bytes()).hexdigest() == "925b88bec3db4a2f47f7d68acc4bab67d8796049fdd021e0061cf066f10eb94a"
BRANCHES = next(ast.literal_eval(node.value) for node in ast.parse(GENERATOR.read_text()).body if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id == "BRANCHES" for target in node.targets))
assert len(BRANCHES) == 17
SMILES = [line.strip() for line in CORPUS.read_text().splitlines() if line.strip() and not line.strip().startswith("#")]
assert len(SMILES) == 5000
POSITIONS = []
with GOLDEN.open("rb") as handle:
    while True:
        position = handle.tell()
        if not handle.readline(): break
        POSITIONS.append(position)
assert len(POSITIONS) == 5000
ACTUAL = (ROOT / "target/agent-handoff/fp-graph-block/morgan-native-5000-u64-full-observations-v3.jsonl").open("x")
METHODS = [("sparse_count", "morgan_sparse_count_fingerprint_with_params"), ("sparse_bit", "morgan_sparse_fingerprint_with_params"), ("hashed_count", "morgan_count_fingerprint_with_params"), ("explicit_bit", "morgan_fingerprint_with_params")]


def params_for(molecule, branch):
    atoms = molecule.num_atoms()
    bonds = molecule.num_bonds()
    custom = branch.get("customAtomInvariants")
    atom_values = None
    if custom == "index_plus_one": atom_values = [i + 1 for i in range(atoms)]
    elif custom == "even_atoms_zero_odd_index_plus_one": atom_values = [0 if i % 2 == 0 else i + 1 for i in range(atoms)]
    invariant = ck.MorganInvariants.features() if branch.get("atomInvariantsGenerator") == "morgan_feature_default" else ck.MorganInvariants.connectivity()
    ring_membership = False if branch.get("atomInvariantsGenerator") == "morgan_ring_false" else branch.get("includeRingMembership", True)
    use_bond_types = False if branch.get("bondInvariantsGenerator") == "morgan_no_bond_types" else branch["useBondTypes"]
    return ck.MorganFingerprintParams(generator=ck.MorganParams(radius=branch["radius"], fp_size=branch["nBits"], include_chirality=branch["includeChirality"], use_bond_types=use_bond_types, count_simulation=branch.get("countSimulation",False), count_bounds=branch.get("countBounds"), only_nonzero_invariants=branch.get("onlyNonzeroInvariants",False), include_ring_membership=ring_membership, include_redundant_environments=branch.get("includeRedundantEnvironments",False), bits_per_feature=branch.get("numBitsPerFeature",1)), from_atoms=[0] if branch.get("fromAtoms")=="first" and atoms else None, ignore_atoms=[0] if branch.get("ignoreAtoms")=="first" and atoms else None, custom_atom_invariants=atom_values, custom_bond_invariants=[i+7 for i in range(bonds)] if branch.get("customBondInvariants")=="index_plus_seven" else None, invariants=invariant)


def metadata_record(output):
    return {"atom_counts":output.atom_counts(), "atom_to_bits":output.atom_to_bits(), "bit_info_map":{str(key):[list(pair) for pair in pairs] for key,pairs in output.bit_info_map().items()}, "atoms_per_bit":{str(key):[list(atoms) for atoms in groups] for key,groups in output.atoms_per_bit().items()}}


@pytest.mark.parametrize("index", range(5000))
def test_current_morgan_matches_original_native_all_branches_and_outputs(index):
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
                output = ck.AdditionalOutput(); output.allocate_atom_counts(); output.allocate_atom_to_bits(); output.allocate_bit_info_map(); output.allocate_atoms_per_bit()
            observed["calls"] += 1
            try:
                fingerprint = getattr(molecule, method)(params, output)
                value = {"ok":True,"error":None}
                if output_name in ("sparse_count","hashed_count"):
                    value["nonzero_elements"]={str(bit):count for bit,count in fingerprint.nonzero_elements().items()}
                else:
                    bits=sorted(x & 0xFFFFFFFF for x in fingerprint.on_bits()); value.update(on_bits=bits,num_on_bits=len(bits))
                if output is not None: value["additional_output"] = metadata_record(output)
            except Exception as error:
                value={"ok":False,"type":type(error).__name__,"error":str(error),"domain":getattr(error,"domain",None),"kind":getattr(error,"kind",None)}
            observed["branches"][name][output_name] = value
            if value != target: errors.append(f"{name}/{output_name}")
    ACTUAL.write(json.dumps(observed,ensure_ascii=True)+"\n"); ACTUAL.flush()
    assert observed["calls"] == 68
    assert not errors, f"row {index} {SMILES[index]} differs: {errors} (full actual bytes retained in observations JSONL)"
