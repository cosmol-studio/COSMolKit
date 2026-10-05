"""Original modern u64 binding and pinned RDKit Torsions.py fixed proposals."""
import pytest
import cosmolkit as ck

@pytest.mark.parametrize("smiles,path,expected", [
    ("C=CC",[0,1,2],(("C",1,0),("C",2,1),("C",1,1))),
    ("C=CO",[0,1,2],(("C",1,1),("C",2,1),("O",1,0))),
    ("OC=CO",[0,1,2,3],(("O",1,0),("C",2,1),("C",2,1),("O",1,0))),
    ("CC=CO",[0,1,2,3],(("C",1,0),("C",2,1),("C",2,1),("O",1,0))),
    ("C=CC(=O)O",[0,1,2,3],(("C",1,1),("C",2,1),("C",3,1),("O",1,1))),
    ("C=CC(=O)O",[0,1,2,4],(("C",1,1),("C",2,1),("C",3,1),("O",1,0))),
    ("OOOO",[0,1,2],(("O",1,0),("O",2,0),("O",2,0))),
    ("OOOO",[0,1,2,3],(("O",1,0),("O",2,0),("O",2,0),("O",1,0))),
])
def test_complete_source_doc_examples_direction_and_exact_nested_tuple(smiles,path,expected):
    molecule=ck.Molecule.from_smiles(smiles)
    score=molecule.topological_torsion_path_score(path,len(path))
    assert ck.explain_path_score(score,len(path))==expected
    assert type(ck.explain_path_score(score,len(path))) is tuple
    assert all(type(x) is tuple for x in ck.explain_path_score(score,len(path)))
    assert molecule.topological_torsion_path_score(path[::-1],len(path))==score
    assert molecule.topological_torsion_path_score(path+[999],len(path),[])==score


def test_fixed_packing_values_custom_codes_and_source_nonpath_inputs():
    m=ck.Molecule.from_smiles("C=CC")
    assert m.topological_torsion_path_score([0,1,2],3)==10506272
    assert ck.Molecule.from_smiles("C=CO").topological_torsion_path_score([0,1,2],3)==25186344
    m=ck.Molecule.from_smiles("C.C.O")
    expected=32 | (96<<9) | (32<<18) | (96<<27)
    assert m.topological_torsion_path_score([0,2,0,2],4)==expected
    assert m.topological_torsion_path_score([2,0,2,0],4)==expected
    m=ck.Molecule.from_smiles("CC")
    assert m.topological_torsion_path_score([0,1],2,[1,512])==261632
    assert ck.explain_path_score(261632,2)==(("B",1,0),("*",8,3))

@pytest.mark.parametrize("path,size,codes,kind,details", [
    ([],0,None,"ZeroSize",{}),([0],2,None,"ShortPath",{"actual":1,"required":2}),
    ([9],1,[4],"ShortAtomCodes",{"actual":1,"required":3}),
    ([3],1,[4,4,4],"AtomIndexOutOfRange",{"index":3,"atom_count":3}),
    ([0,1,2],3,[10,1,10],"AtomCodeUnderflow",{"index":1,"code":1,"subtract":2}),
    ([0]*9,9,[10]*3,"PackedCode",{}),
])
def test_original_error_categories_order_and_typed_context(path,size,codes,kind,details):
    m=ck.Molecule.from_smiles("CCC")
    category=IndexError if kind=="AtomIndexOutOfRange" else ck.TopologicalTorsionPathScoreError
    with pytest.raises(category) as raised:m.topological_torsion_path_score(path,size,codes)
    assert raised.value.domain=="Fingerprint"
    assert raised.value.kind==kind
    for field,value in details.items():assert getattr(raised.value,field)==value
    if kind=="PackedCode":assert raised.value.__cause__ is not None


def test_default_size_zero_and_extra_zero_chunks_unsigned_score_boundary():
    assert ck.explain_path_score(0)==(("B",1,0),("B",2,0),("B",2,0),("B",1,0))
    assert ck.explain_path_score((1<<64)-1,0)==()
    decoded=ck.explain_path_score((1<<64)-1,10)
    assert decoded[0]==("*",8,3) and decoded[1]==("*",9,3)
    assert decoded[8:]==(("B",2,0),("B",1,0))
    for invalid in [-1,1<<64]:
        with pytest.raises(OverflowError):ck.explain_path_score(invalid)
    with pytest.raises(OverflowError):ck.explain_path_score(0,-1)
    with pytest.raises(TypeError):ck.explain_path_score("0")
    m=ck.Molecule.from_smiles("CCC")
    for path,size,codes in [([-1],1,None),([0],-1,None),([0],1,[-1]*3),([0],1,[1<<32]*3)]:
        with pytest.raises(OverflowError):m.topological_torsion_path_score(path,size,codes)


def test_generated_path_error_class_matches_native_public_fields():
    import ast
    from pathlib import Path

    stub = ast.parse((Path(__file__).resolve().parents[2] / "python/cosmolkit.pyi").read_text())
    declarations = [node for node in stub.body if isinstance(node, ast.ClassDef) and node.name == "TopologicalTorsionPathScoreError"]
    assert len(declarations) == 1
    declaration = declarations[0]
    assert [ast.unparse(base) for base in declaration.bases] == ["builtins.ValueError"]
    assert {node.target.id: ast.unparse(node.annotation) for node in declaration.body if isinstance(node, ast.AnnAssign)} == {
        "domain": "builtins.str", "kind": "builtins.str",
        "actual": "builtins.int", "required": "builtins.int",
        "index": "builtins.int", "atom_count": "builtins.int",
        "code": "builtins.int", "subtract": "builtins.int",
    }
    exported = next(ast.literal_eval(node.value) for node in stub.body if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id == "__all__" for target in node.targets))
    assert exported.count("TopologicalTorsionPathScoreError") == 1
    assert issubclass(ck.TopologicalTorsionPathScoreError, ValueError)
