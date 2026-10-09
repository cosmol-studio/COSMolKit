"""Original Layered scalar and query observations through canonical projections."""
from pathlib import Path
import json
import cosmolkit as ck
import pytest

DEFAULT_BITS = [92, 360, 596, 610, 611, 674, 867, 1044, 1111, 1783, 1784]
FIXTURE = Path(__file__).resolve().parents[2] / "testdata/fingerprint/fixtures/rdkit/layered_fingerprint_query_cases.json"
SOURCE_CASES = json.loads(FIXTURE.read_text())
QUERY_CASES = SOURCE_CASES["complexity_masks"] + SOURCE_CASES["aromaticity_branches"]

def test_original_defaults_counts_masks_roots_and_immutable_results():
    mol = ck.Molecule.from_smiles("CCO")
    before = mol.to_smiles()
    assert mol.fingerprint_layered().on_bits() == DEFAULT_BITS
    assert mol.fingerprint_layered().n_bits() == 2048
    absent = mol.fingerprint_layered_with_output()
    assert absent.atom_counts() is None
    assert repr(absent) == "LayeredFingerprintResult(n_bits=2048, has_atom_counts=false)"
    seed = [10, 20, 30]
    p = ck.LayeredFingerprintParams(atom_counts=seed, set_only_bits=ck.Fingerprint.from_on_bits(2048, [674]))
    result = mol.fingerprint_layered_with_output_with_params(p)
    assert result.fingerprint().on_bits() == [674]
    assert result.atom_counts() == [11, 22, 31]
    assert repr(result) == "LayeredFingerprintResult(n_bits=2048, has_atom_counts=true)"
    copy = result.atom_counts()
    copy.clear()
    assert result.atom_counts() == [11, 22, 31]
    assert seed == p.atom_counts == [10, 20, 30]
    assert mol.fingerprint_layered_with_params(ck.LayeredFingerprintParams(from_atoms=[])).on_bits() == []
    assert mol.fingerprint_layered_with_params(ck.LayeredFingerprintParams(
        from_atoms=[0], branched_paths=False)).on_bits() == [360,596,610,611,674,867,1044,1111,1783,1784]
    assert mol.to_smiles() == before

def test_original_unsigned_layer_flags_and_configuration_defaults():
    p = ck.LayeredFingerprintParams()
    assert (p.layers,p.min_path,p.max_path,p.fp_size,p.branched_paths) == (0xffff_ffff,1,7,2048,True)
    assert p.atom_counts is p.set_only_bits is p.from_atoms is None
    assert ck.LayeredFingerprintLayers.from_bits_retain(0xffff_ffc0).bits() == 0xffff_ffc0
    mol = ck.Molecule.from_smiles("CCO")
    result = mol.fingerprint_layered_with_output_with_params(
        ck.LayeredFingerprintParams(layers=0xffff_ffc0, atom_counts=[5,6,7]))
    assert result.fingerprint().on_bits() == []
    assert result.atom_counts() == [5,6,7]
    p.min_path = 2
    assert p.min_path == 2

@pytest.mark.parametrize("kwargs,reason", [
    ({"min_path":0},"minPath==0"),
    ({"min_path":3,"max_path":2},"maxPath<minPath"),
    ({"fp_size":0},"fpSize==0"),
    ({"atom_counts":[0,0]},"bad atomCounts size"),
    ({"from_atoms":[3]},"fromAtoms contains atom index out of range"),
])
def test_original_preconditions_remain_structured_errors(kwargs,reason):
    mol = ck.Molecule.from_smiles("CCO")
    with pytest.raises(ck.LayeredFingerprintError) as raised:
        mol.fingerprint_layered_with_params(ck.LayeredFingerprintParams(**kwargs))
    assert str(raised.value) == reason
    assert raised.value.reason == reason
    assert raised.value.domain == "Fingerprint"
    assert raised.value.kind == "InvalidArguments"

def test_original_mask_width_failure_and_unsigned_transport():
    mol = ck.Molecule.from_smiles("CCO")
    with pytest.raises(ck.LayeredFingerprintError,match="bad setOnlyBits size"):
        mol.fingerprint_layered_with_params(ck.LayeredFingerprintParams(set_only_bits=ck.Fingerprint.from_on_bits(64, [])))
    for kwargs in [{"layers":-1},{"layers":1<<32},{"from_atoms":[-1]},{"atom_counts":[-1]}]:
        with pytest.raises(OverflowError):
            ck.LayeredFingerprintParams(**kwargs)

@pytest.mark.parametrize("case",QUERY_CASES,ids=[c["case_id"] for c in QUERY_CASES])
def test_all_original_query_masks_aromaticity_and_complete_on_bits(case):
    params = ck.LayeredFingerprintParams(layers=0x3f, fp_size=512)
    if case["notation"] == "smarts":
        query = ck.QueryGraph.from_smarts(case["input"])
        result = ck.fingerprint_layered_query_with_params(query,params)
    else:
        mol = ck.Molecule.from_smiles(case["input"])
        result = mol.fingerprint_layered_with_params(params)
    assert result.n_bits() == 512
    assert result.on_bits() == case["on_bits"]

def test_canonical_query_owned_output_protocol():
    query = ck.QueryGraph.from_smarts("[C,N]~[C,N]")
    params = ck.LayeredFingerprintParams(fp_size=512,atom_counts=[5,6])
    result = ck.fingerprint_layered_query_with_output_with_params(query,params)
    assert result.fingerprint().on_bits() == [162]
    assert result.atom_counts() == [6,7]
    assert params.atom_counts == [5,6]


def test_original_source_layered_scalar_defaults_types_counts_masks_roots_and_immutability():
    molecule = ck.Molecule.from_smiles('CCO')
    before = molecule.to_smiles()
    default = molecule.fingerprint_layered()
    assert default.n_bits() == 2048
    assert default.on_bits() == [92, 360, 596, 610, 611, 674, 867, 1044, 1111, 1783, 1784]
    assert molecule.to_smiles() == before
    topology_mask = molecule.fingerprint_layered_with_params(ck.LayeredFingerprintParams(layers=1))
    assert topology_mask.on_bits() == [674, 867]
    masked = molecule.fingerprint_layered_with_output_with_params(ck.LayeredFingerprintParams(atom_counts=[10, 20, 30], set_only_bits=topology_mask))
    assert masked.fingerprint().on_bits() == [674, 867]
    assert masked.atom_counts() == [12, 23, 32]
    no_counts = molecule.fingerprint_layered_with_output()
    assert no_counts.atom_counts() is None
    assert no_counts.fingerprint().on_bits() == default.on_bits()
    assert molecule.fingerprint_layered_with_params(ck.LayeredFingerprintParams(from_atoms=[])).on_bits() == []
    assert molecule.fingerprint_layered_with_params(ck.LayeredFingerprintParams(branched_paths=False, from_atoms=[0])).on_bits() == [360, 596, 610, 611, 674, 867, 1044, 1111, 1783, 1784]
    assert molecule.to_smiles() == before

def test_original_source_layered_scalar_preserves_source_errors_and_unsigned_layer_bits():
    molecule = ck.Molecule.from_smiles('CCO')
    assert molecule.fingerprint_layered_with_params(ck.LayeredFingerprintParams(layers=4294967232, atom_counts=[5, 6, 7])).on_bits() == []
    with pytest.raises(ValueError, match='minPath==0'):
        molecule.fingerprint_layered_with_params(ck.LayeredFingerprintParams(min_path=0))
    with pytest.raises(ValueError, match='maxPath<minPath'):
        molecule.fingerprint_layered_with_params(ck.LayeredFingerprintParams(min_path=3, max_path=2))
    with pytest.raises(ValueError, match='bad atomCounts size'):
        molecule.fingerprint_layered_with_params(ck.LayeredFingerprintParams(atom_counts=[0, 0]))
    with pytest.raises(ValueError, match='bad setOnlyBits size'):
        molecule.fingerprint_layered_with_params(ck.LayeredFingerprintParams(fp_size=64, set_only_bits=molecule.fingerprint_layered()))
    with pytest.raises(ValueError, match='fromAtoms contains atom index out of range'):
        molecule.fingerprint_layered_with_params(ck.LayeredFingerprintParams(from_atoms=[3]))
