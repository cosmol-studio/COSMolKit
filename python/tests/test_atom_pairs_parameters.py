"""Pinned RDKit351f8f FingerprintUtil.h and original d892ec3 vocabulary proposal."""
import pytest
import cosmolkit as ck

def fields(value):
    result=(value.symbol,value.branch_count,value.pi_electrons)
    return result if value.chirality is None else result+(value.chirality,)

@pytest.mark.parametrize("name,expected", [
    ("version","1.1.0"),("num_type_bits",4),("num_pi_bits",2),
    ("num_branch_bits",3),("num_chiral_bits",2),("code_size",9),
    ("num_path_bits",5),("max_path_length",31),("num_atom_pair_fingerprint_bits",23),
])
def test_exact_original_constant_values(name,expected):
    assert getattr(ck.AtomPairsParameters,name)==expected


def test_atom_types_complete_order_and_implicit_zero():
    assert ck.AtomPairsParameters.atom_types == [5,6,7,8,9,14,15,16,17,33,34,35,51,52,53,0]
    detached=list(ck.AtomPairsParameters.atom_types)
    detached[0]=99
    assert ck.AtomPairsParameters.atom_types[0]==5


@pytest.mark.parametrize("index,symbol", enumerate(["B","C","N","O","F","Si","P","S","Cl","As","Se","Br","Sb","Te","I","*"]))
def test_type_order_drives_existing_explanation(index,symbol):
    p=ck.AtomPairsParameters
    code=(index << (p.num_branch_bits+p.num_pi_bits)) | 7 | (3 << p.num_branch_bits)
    assert fields(ck.AtomCodeExplanation.from_code(code))==(symbol,7,3)
    for bits,label in [(0,""),(1,"R"),(2,"S")]:
        assert fields(ck.AtomCodeExplanation.from_code(code | (bits << p.code_size),include_chirality=True))==(symbol,7,3,label)


def test_namespace_has_no_constructor():
    with pytest.raises(TypeError):
        ck.AtomPairsParameters()


def test_existing_atom_code_producer_uses_source_vocabulary():
    molecule=ck.Molecule.from_smiles("CCO")
    p=ck.AtomPairsParameters
    for index,expected in enumerate([("C",1,0),("C",2,0),("O",1,0)]):
        code=molecule.with_atom_pair_atom_code(index).code
        assert code < (1 << p.code_size)
        assert fields(ck.AtomCodeExplanation.from_code(code))==expected
