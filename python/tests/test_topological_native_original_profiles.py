"""Complete original5152 rows/all14 RDK branches, including source requested provenance.
References are preflighted in full before CK calls; tests never prepare an oracle.
Supply the two raw references through COSMOLKIT_RDK_ORIGINAL_SMALL_GOLDEN and
COSMOLKIT_RDK_ORIGINAL5000_GOLDEN, as the Rust corpus test does. Supply source
requested atomBits/bitInfo through COSMOLKIT_RDK_ORIGINAL_SMALL_REQUESTED_OUTPUT_GOLDEN
and COSMOLKIT_RDK_ORIGINAL5000_REQUESTED_OUTPUT_GOLDEN. Missing/stale references fail.
References come from tools/testdata/rdkit/_generate_rdkit_topological_fingerprint_golden.py
and its original14 profile, with source provenance enabled for requested outputs;
the pinned source/parameters, complete input identities and expected hashes stay fixed.
"""
from pathlib import Path
import hashlib
import json
import os
import cosmolkit as ck
import pytest

ROOT = Path(__file__).resolve().parents[2]
REFERENCES = [{'name': 'small152', 'input_file': 'smiles_small.smi', 'count': 152, 'input_sha256': '47380e477dc2ab4b3c2b7cd62754e52b718bb9f3f0b977d47eb97db45074442e', 'raw': os.environ['COSMOLKIT_RDK_ORIGINAL_SMALL_GOLDEN'], 'raw_sha256': '1a1e1fd9bcb60e702814d984434c68a1bb5726c045f7bb1bf2b7af7ed3f11fb1', 'requested_output': os.environ['COSMOLKIT_RDK_ORIGINAL_SMALL_REQUESTED_OUTPUT_GOLDEN'], 'requested_output_sha256': 'b2180c7266bde9b6118f820446d3ffe9062e42d01a4b6fb0b90561f6192c8ed7', 'requested_output_bytes': 132819946}, {'name': 'original5000', 'input_file': 'smiles_5000.smi', 'count': 5000, 'input_sha256': 'a4d579cd72621af27772256bb23ba796452276bb924fd20aac83625ffa67d849', 'raw': os.environ['COSMOLKIT_RDK_ORIGINAL5000_GOLDEN'], 'raw_sha256': '0d443196f862ad5d2fa0c9be3c16b8fc3464800d0ea64f87f894b2eba1f4cc69', 'requested_output': os.environ['COSMOLKIT_RDK_ORIGINAL5000_REQUESTED_OUTPUT_GOLDEN'], 'requested_output_sha256': 'ce11fef76dd0f950d412e8c72b8e78f9f2dd7819e7acf61fb7438e1c24f1856d', 'requested_output_bytes': 6690736867}]
PROFILE_SHA = "227d51b15c6387885fd38502ed69b8a48691700d343f2fd3103ae9e0738f34ce"

def file_sha(path):
    digest=hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda:handle.read(1<<20),b""):digest.update(block)
    return digest.hexdigest()

@pytest.fixture(scope="session")
def preflight_all_selected_original_rows():
    profile_path=ROOT/"tools/testdata/rdkit/rdkit_topological_fingerprint_profile.json"
    assert file_sha(profile_path)==PROFILE_SHA
    branches=json.loads(profile_path.read_text())["branches"]
    assert len(branches)==14
    index=[]
    for reference in REFERENCES:
        input_path=ROOT/"testdata/smiles/corpus"/reference["input_file"]
        assert file_sha(input_path)==reference["input_sha256"]
        corpus=[line.strip() for line in input_path.read_text().splitlines() if line.strip() and not line.strip().startswith("#")]
        assert len(corpus)==reference["count"]
        assert file_sha(reference["raw"])==reference["raw_sha256"]
        assert file_sha(reference["requested_output"])==reference["requested_output_sha256"]
        with Path(reference["raw"]).open("rb") as raw,Path(reference["requested_output"]).open("rb") as output:
            for row,smiles in enumerate(corpus):
                raw_offset=raw.tell();output_offset=output.tell()
                ordinary=json.loads(raw.readline());requested=json.loads(output.readline())
                for record in [ordinary,requested]:
                    assert record["row"]==row
                    assert record["smiles"]==smiles
                    assert record["rdkit_ok"]==ordinary["rdkit_ok"]
                if ordinary["rdkit_ok"]:
                    for record in [ordinary,requested]:
                        assert record["error"] is None
                        assert len(record["branches"])==14
                    for branch in branches:
                        name=branch["name"];a=ordinary["branches"][name];b=requested["branches"][name]
                        assert a["parameters"]==branch
                        assert b["parameters"]==dict(branch,provenance=True)
                        for field in ["ok","on_bits","num_bits","num_on_bits","error"]:assert a[field]==b[field]
                        assert a["ok"] and a["error"] is None
                        assert a["num_on_bits"]==len(a["on_bits"])
                        assert isinstance(b["atom_bits"],list) and isinstance(b["bit_info"],dict)
                else:
                    assert ordinary["error"]==requested["error"]
                    assert not ordinary["branches"] and not requested["branches"]
                index.append((reference,row,smiles,raw_offset,output_offset))
            assert not raw.read(1) and not output.read(1)
    assert len(index)==5152
    destination=Path(os.environ["COSMOLKIT_RDK_NATIVE_OBSERVATIONS"])
    with destination.open("x") as observations:
        yield index,observations


def params_for_branch(branch,atom_count):
    mapping={"minPath":"min_path","maxPath":"max_path","fpSize":"fp_size","nBitsPerHash":"num_bits_per_feature","useHs":"use_hs","tgtDensity":"target_density","minSize":"min_size","branchedPaths":"branched_paths","useBondOrder":"use_bond_order"}
    kwargs={target:branch[source] for source,target in mapping.items()}
    if branch.get("fromAtoms")=="first" and atom_count>0:kwargs["from_atoms"]=[0]
    if branch.get("atomInvariants")=="index_plus_one":kwargs["atom_invariants"]=list(range(1,atom_count+1))
    return ck.TopologicalFingerprintParams(**kwargs)

@pytest.mark.parametrize("selected_row",range(5152))
def test_complete_original14_scalar_and_requested_metadata_every_row(selected_row,preflight_all_selected_original_rows):
    index,observations=preflight_all_selected_original_rows
    reference,row,smiles,raw_offset,output_offset=index[selected_row]
    with Path(reference["requested_output"]).open("rb") as handle:
        handle.seek(output_offset);expected=json.loads(handle.readline())
    actual={"profile":reference["name"],"row":row,"smiles":smiles,"branches":{}}
    if not expected["rdkit_ok"]:
        with pytest.raises(ck.SmilesError) as raised:ck.Molecule.from_smiles(smiles)
        actual.update(ok=False,error=str(raised.value),domain=raised.value.domain,kind=raised.value.kind)
        observations.write(json.dumps(actual,sort_keys=True)+"\n");return
    molecule=ck.Molecule.from_smiles(smiles);before=molecule.to_smiles()
    for name,branch in expected["branches"].items():
        context=f"{reference['name']} row {row} {smiles} / {name}"
        params=params_for_branch(branch["parameters"],molecule.num_atoms())
        scalar=molecule.topological_fingerprint_with_params(params)
        result=molecule.topological_fingerprint_with_output_with_params(params,ck.TopologicalFingerprintOutputRequest(atom_bits=True,bit_info=True))
        fp=result.fingerprint();atom_bits=result.atom_bits();bit_info=result.bit_info()
        actual["branches"][name]={"scalar_on_bits":scalar.on_bits(),"num_bits":fp.n_bits(),"on_bits":fp.on_bits(),"num_on_bits":len(fp.on_bits()),"atom_bits":atom_bits,"bit_info":{str(bit):paths for bit,paths in bit_info.items()}}
        assert scalar.n_bits()==branch["num_bits"],context
        assert scalar.on_bits()==branch["on_bits"],context
        assert fp.n_bits()==branch["num_bits"],context
        assert fp.on_bits()==branch["on_bits"],context
        assert len(fp.on_bits())==branch["num_on_bits"],context
        assert atom_bits==branch["atom_bits"],context
        assert bit_info=={int(bit):paths for bit,paths in branch["bit_info"].items()},context
        assert params.from_atoms==([0] if branch["parameters"].get("fromAtoms")=="first" and molecule.num_atoms()>0 else None),context
    assert molecule.to_smiles()==before
    actual["ok"]=True
    observations.write(json.dumps(actual,sort_keys=True)+"\n")
