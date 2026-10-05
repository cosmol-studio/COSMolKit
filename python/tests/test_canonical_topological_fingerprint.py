"""Complete original3 RDK scalar conditions, frozen source observations and canonical projection."""
import cosmolkit as ck
import pytest

def test_topological_fingerprint_matches_rdkit_exact_bits_and_is_deterministic():
    molecule = ck.Molecule.from_smiles("CCO")
    params = ck.TopologicalFingerprintParams(fp_size=64,num_bits_per_feature=1)
    fingerprint = molecule.topological_fingerprint_with_params(params)
    assert fingerprint.on_bits() == [0,28,59]
    assert molecule.topological_fingerprint_with_params(params).on_bits() == fingerprint.on_bits()

def test_topological_fingerprint_with_output_matches_rdkit_provenance():
    molecule = ck.Molecule.from_smiles("CCO")
    params = ck.TopologicalFingerprintParams(fp_size=64,num_bits_per_feature=1)
    result = molecule.topological_fingerprint_with_output_with_params(params,ck.TopologicalFingerprintOutputRequest(atom_bits=True,bit_info=True))
    assert result.fingerprint().on_bits() == [0,28,59]
    assert result.atom_bits() == [[28,0],[28,59,0],[59,0]]
    assert result.bit_info() == {0:[[0,1]],28:[[0]],59:[[1]]}

def test_topological_fingerprint_rejects_source_precondition_ranges():
    molecule = ck.Molecule.from_smiles("CCO")
    with pytest.raises(ValueError,match="minPath==0"):
        molecule.topological_fingerprint_with_params(ck.TopologicalFingerprintParams(min_path=0))

def test_original11_defaults_and_source_absent_output_getters():
    params = ck.TopologicalFingerprintParams()
    assert (params.min_path,params.max_path,params.fp_size,params.num_bits_per_feature,params.use_hs,params.target_density,params.min_size,params.branched_paths,params.use_bond_order,params.atom_invariants,params.from_atoms)==(1,7,2048,2,True,0.0,128,True,True,None,None)
    request = ck.TopologicalFingerprintOutputRequest()
    assert (request.atom_bits,request.bit_info)==(False,False)
    for value,field in [(params,"min_path"),(request,"atom_bits")]:
        with pytest.raises(AttributeError): setattr(value,field,2)
    molecule = ck.Molecule.from_smiles("CCO")
    default = molecule.topological_fingerprint()
    assert default.on_bits()==[562,1183,1308,1339,1728,1772]
    result = molecule.topological_fingerprint_with_output()
    assert result.fingerprint().on_bits()==default.on_bits()
    assert repr(result)=="TopologicalFingerprintResult(n_bits=2048, has_atom_bits=false, has_bit_info=false)"
    for field in ["atom_bits","bit_info"]:
        with pytest.raises(ck.TopologicalFingerprintError) as raised: getattr(result,field)()
        assert str(raised.value)==f"topological {field} output was not requested for this fingerprint result"
        assert raised.value.domain=="Fingerprint"
        assert raised.value.kind=="OutputNotRequested"
        assert raised.value.field==field

@pytest.mark.parametrize("kwargs,reason",[({"min_path":0},"minPath==0"),({"max_path":0},"maxPath<minPath"),({"fp_size":0},"fpSize==0"),({"num_bits_per_feature":0},"nBitsPerHash==0")])
def test_original4_preconditions_structured_projection(kwargs,reason):
    with pytest.raises(ck.TopologicalFingerprintError) as raised:
        ck.Molecule.from_smiles("CCO").topological_fingerprint_with_params(ck.TopologicalFingerprintParams(**kwargs))
    assert str(raised.value)==reason
    assert raised.value.reason==reason
    assert raised.value.domain=="Fingerprint"
    assert raised.value.kind=="InvalidArguments"

def test_original_source_pre_fold_metadata_and_detached_native_values():
    molecule = ck.Molecule.from_smiles("CCO");before=molecule.to_smiles()
    params = ck.TopologicalFingerprintParams(fp_size=64,num_bits_per_feature=1,target_density=0.2,min_size=16)
    result = molecule.topological_fingerprint_with_output_with_params(params,ck.TopologicalFingerprintOutputRequest(atom_bits=True,bit_info=True))
    assert result.fingerprint().n_bits()==16
    assert result.atom_bits()==[[28,0],[28,59,0],[59,0]]
    assert result.bit_info()=={0:[[0,1]],28:[[0]],59:[[1]]}
    atom_copy=result.atom_bits();atom_copy[0].clear()
    bit_copy=result.bit_info();bit_copy[0][0].clear()
    assert result.atom_bits()==[[28,0],[28,59,0],[59,0]]
    assert result.bit_info()=={0:[[0,1]],28:[[0]],59:[[1]]}
    roots=[];params=ck.TopologicalFingerprintParams(from_atoms=roots);roots.append(0)
    assert params.from_atoms==[]
    assert molecule.topological_fingerprint_with_params(params).on_bits()==[]
    assert molecule.to_smiles()==before


def test_original_complex_query_bond_skips_environments_through_canonical_projection():
    query = ck.QueryGraph.from_smarts("[#6]~[#6]")
    params = ck.TopologicalFingerprintParams()
    fingerprint = ck.topological_query_fingerprint_with_params(query,params)
    assert fingerprint.n_bits()==2048
    assert fingerprint.on_bits()==[]
    result = ck.topological_query_fingerprint_with_output_with_params(query,params,ck.TopologicalFingerprintOutputRequest(atom_bits=True,bit_info=True))
    assert result.fingerprint().on_bits()==[]
    assert result.atom_bits()==[[],[]]
    assert result.bit_info()=={}
