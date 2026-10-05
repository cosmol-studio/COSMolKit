"""Original full 5000 inputs, nine profiles, four single-call forms/AO proposal.
Bulk, generator metadata, JSON/legacy/helper calls are retained in the reference
and remain separately pending. Independent p1 review required.
"""
import hashlib
import json
from pathlib import Path
import pytest
import cosmolkit as ck
ROOT = Path(__file__).resolve().parents[2]
def sha(path):
    digest=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):digest.update(b)
    return digest.hexdigest()
CORPUS=ROOT/"testdata/smiles/corpus/smiles_5000.smi"
GOLDEN=ROOT/"target/agent-handoff/structure-oracle-native-env/generated-source-native-v1/topological_torsion_fingerprint.jsonl"
GENERATOR=ROOT/"tools/testdata/rdkit/_generate_topological_torsion_fingerprint_golden.py"
PROFILE=ROOT/"tools/testdata/rdkit/topological_torsion_fingerprint_profile.json"
assert sha(CORPUS)=="a4d579cd72621af27772256bb23ba796452276bb924fd20aac83625ffa67d849"
assert sha(GOLDEN)=="28d91a8ca475ffedf43ac0e50621be2b0b6ae394225dffbe76b6fff061affe7d"
assert sha(GENERATOR)=="820ca13d717f691890bbc3e80cc803096df3a1515358bf83ecf47c89b6cb0ee4"
assert sha(PROFILE)=="f03a8b5a1e9b9f05347c0058761f7cac2efed4ea0baafce4bbcb967e44040daa"
BRANCHES=json.loads(PROFILE.read_bytes())["corpus_branches"]
assert len(BRANCHES)==9
SMILES=[s.strip() for s in CORPUS.read_text().splitlines() if s.strip() and not s.strip().startswith("#")]
assert len(SMILES)==5000
POSITIONS=[]
with GOLDEN.open("rb") as f:
    while True:
        p=f.tell()
        if not f.readline():break
        POSITIONS.append(p)
assert len(POSITIONS)==5000
ACTUAL=(ROOT/"target/agent-handoff/fp-graph-block/topological-torsion-native-5000-single-current-observations-v1.jsonl").open("x")
METHODS=[("sparse_count","topological_torsion_sparse_count_fingerprint_with_params"),("sparse_bit","topological_torsion_sparse_fingerprint_with_params"),("count","topological_torsion_count_fingerprint_with_params"),("bit","topological_torsion_fingerprint_with_params")]
def params_for(m,b):
    atoms=m.num_atoms()
    return ck.TopologicalTorsionFingerprintParams(generator=ck.TopologicalTorsionParams(torsion_atom_count=b.get("torsionAtomCount",4),include_chirality=b.get("includeChirality",False),only_shortest_paths=b.get("onlyShortestPaths",False),count_simulation=b.get("countSimulation",True),count_bounds=b.get("countBounds"),fp_size=b.get("fpSize",2048),bits_per_feature=b.get("numBitsPerFeature",1)),from_atoms=[0] if b.get("fromAtoms")=="first" and atoms else None,ignore_atoms=[0] if b.get("ignoreAtoms")=="first" and atoms else None,custom_atom_invariants=[i+17 for i in range(atoms)] if b.get("customAtomInvariants")=="index_plus_17" else None)
def metadata_record(o):
    return {"atom_counts":o.atom_counts(),"atom_to_bits":o.atom_to_bits(),"bit_info_map":{str(k):[list(p) for p in v] for k,v in o.bit_info_map().items()},"bit_paths":{str(k):[list(p) for p in v] for k,v in o.bit_paths().items()},"atoms_per_bit":{str(k):[list(p) for p in v] for k,v in o.atoms_per_bit().items()}}
@pytest.mark.parametrize("index",range(5000))
def test_current_torsion_single_forms_match_original_native_all9profiles(index):
    with GOLDEN.open("rb") as f:
        f.seek(POSITIONS[index]); expected=json.loads(f.readline())
    assert expected["smiles"]==SMILES[index] and expected["rdkit_ok"]
    observed={"index":index,"smiles":SMILES[index],"profiles":{},"calls":0};errors=[]
    try:m=ck.Molecule.from_smiles(SMILES[index])
    except Exception as e:
        observed["parse_error"]={"type":type(e).__name__,"message":str(e)}
        ACTUAL.write(json.dumps(observed)+"\n");ACTUAL.flush();pytest.fail(str(observed))
    for branch in BRANCHES:
        name=branch["name"];target=expected["profiles"][name];observed["profiles"][name]={};params=params_for(m,branch)
        for form,method in METHODS:
            output=None
            if branch.get("additionalOutput"):
                output=ck.AdditionalOutput()
                output.allocate_atom_to_bits();output.allocate_atom_counts();output.allocate_bit_info_map();output.allocate_bit_paths();output.allocate_atoms_per_bit()
            observed["calls"]+=1
            try:
                fp=getattr(m,method)(params,output)
                value={"size":fp.length() if form in ("sparse_count","count") else fp.n_bits()}
                if form in ("sparse_count","count"):value["nonzero_elements"]={str(k):v for k,v in fp.nonzero_elements().items()}
                else:value["on_bits"]=sorted(k & 0xFFFFFFFF for k in fp.on_bits()) if form=="sparse_bit" else fp.on_bits()
                meta=metadata_record(output) if output is not None else None
                observed["profiles"][name][form]={"value":value,"additional_output":meta}
                if value!=target[form]:errors.append(f"{name}/{form}/value")
                if meta is not None and meta!=target["additional_output"][form]:errors.append(f"{name}/{form}/AO")
            except Exception as e:
                observed["profiles"][name][form]={"error":{"type":type(e).__name__,"message":str(e),"domain":getattr(e,"domain",None),"kind":getattr(e,"kind",None),"cause":str(e.__cause__)}};errors.append(f"{name}/{form}/exception")
    ACTUAL.write(json.dumps(observed,ensure_ascii=True)+"\n");ACTUAL.flush()
    assert observed["calls"]==36
    assert not errors,f"row {index} {SMILES[index]} differs {errors}; full actual bytes retained"
