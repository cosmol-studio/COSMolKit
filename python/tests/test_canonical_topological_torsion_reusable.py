"""Source-fixed reusable torsion public-chain tests; independent review proposal."""
import ast
import json
from pathlib import Path
import pytest
import cosmolkit as ck
SINGLE=[("fingerprint_topological_torsion_sparse_count_with_generator","sparse_counts"),("fingerprint_topological_torsion_sparse_with_generator","sparse_fingerprints"),("fingerprint_topological_torsion_count_with_generator","counts"),("fingerprint_topological_torsion_with_generator","fingerprints")]

def record(value):
    if hasattr(value,"nonzero_elements"):return {"size":value.length(),"nonzero_elements":value.nonzero_elements()}
    return {"size":value.n_bits(),"on_bits":value.on_bits()}

@pytest.mark.parametrize("single,bulk",SINGLE)
def test_reusable_scalar_source_values_and_bulk_none_order(single,bulk):
    g=ck.TopologicalTorsionFingerprintGenerator()
    a,b,c=(ck.Molecule.from_smiles(s) for s in ["CCCCO","CCCC","C1CC1"])
    expected=[getattr(m,single)(g) for m in [a,b,c]]
    for threads in [1,2,7]:
        values=getattr(g,bulk)([a,None,b,c,None],num_threads=threads)
        assert len(values)==5 and values[1] is None and values[4] is None
        assert [record(values[i]) for i in [0,2,3]]==[record(v) for v in expected]
    if bulk=="sparse_counts":
        assert [v.nonzero_elements() for v in expected]==[{4437590048:1,12893306913:1},{4303372320:1},{4437590049:1}]
        assert all(v.length()==68719476736 for v in expected)
    if bulk=="counts":assert [v.nonzero_elements() for v in expected]==[{0:2},{15:1},{71:1}]
    if bulk=="fingerprints":assert [v.on_bits() for v in expected]==[[0,1],[60],[284]]
    assert getattr(g,bulk)([],num_threads=2)==[]
    assert getattr(g,bulk)([None,None],num_threads=7)==[None,None]


def test_shared_settings_snapshot_lifetime_and_independent_atom_generator():
    g=ck.TopologicalTorsionFingerprintGenerator();s=g.settings();alias=g.settings();m=ck.Molecule.from_smiles("CCCCO")
    snapshot=s.params();assert snapshot.torsion_atom_count==4
    s.torsion_atom_count=3;s.fp_size=1000;s.only_shortest_paths=True
    assert alias.torsion_atom_count==3 and alias.fp_size==1000
    assert m.fingerprint_topological_torsion_sparse_count_with_generator(g).length()==1<<27
    assert sum(m.fingerprint_topological_torsion_sparse_count_with_generator(g).nonzero_elements().values())==3
    assert m.fingerprint_topological_torsion_count_with_generator(g).length()==1000
    assert snapshot.torsion_atom_count==4 and snapshot.fp_size==2048
    s.include_chirality=True
    value=json.loads(g.to_json());assert value["fingerprintArguments"]["includeChirality"]=="true"
    assert value["atomInvariantsGenerator"]["includeChirality"]=="false"
    assert m.fingerprint_topological_torsion_sparse_count_with_generator(g).length()==1<<33
    bounds=[1,3,5];s.count_bounds=bounds;bounds.append(7);assert alias.count_bounds==[1,3,5]
    copied=s.count_bounds;copied.append(9);assert s.count_bounds==[1,3,5]
    s.count_bounds=(1, 2, 4, 8);assert alias.count_bounds==[1, 2, 4, 8]
    del g;alias.fp_size=777;assert s.fp_size==777
    with pytest.raises(TypeError):ck.TopologicalTorsionSettings()


def test_source_metadata_json_roundtrip_and_call_selections():
    g=ck.TopologicalTorsionFingerprintGenerator();expected="Common arguments : countSimulation=1 fpSize=2048 bitsPerFeature=1 includeChirality=0 --- TopologicalTorsionArguments torsionAtomCount=4 onlyShortestPaths=0 --- TopologicalTorsionEnvGenerator --- AtomPairInvariantGenerator topologicalTorsionCorrection=1 --- No bond invariants generator"
    assert g.info_string()==expected;assert repr(g)==f"TopologicalTorsionFingerprintGenerator({expected})"
    value=json.loads(g.to_json());assert value["fingerprintArguments"]["countBounds"]==["1","2","4","8"]
    restored=ck.TopologicalTorsionFingerprintGenerator.from_json(g.to_json());assert json.loads(restored.to_json())==value
    m=ck.Molecule.from_smiles("CCCCO");before=m.to_smiles();custom=[17,18,19,20,21];params=ck.TopologicalTorsionCallParams(custom_atom_invariants=custom)
    custom.append(99);assert params.custom_atom_invariants==[17,18,19,20,21]
    for method,_ in SINGLE:
        assert record(getattr(m,method)(restored,params=params))==record(getattr(m,method)(g,params=params))
        empty = record(getattr(m,method)(g,params=ck.TopologicalTorsionCallParams(from_atoms=[])))
        assert empty["nonzero_elements"] == {} if "nonzero_elements" in empty else empty["on_bits"] == []
    assert m.to_smiles()==before
    params.conformer_id=0
    assert params.conformer_id == 0
    with pytest.raises(ck.TopologicalTorsionReadError):ck.TopologicalTorsionFingerprintGenerator.from_json("not json")


def test_source_live_dense_bounds_error_preserves_usable_state_and_errors():
    g=ck.TopologicalTorsionFingerprintGenerator();m=ck.Molecule.from_smiles("CCCC");s=g.settings();s.count_bounds=[]
    with pytest.raises(ck.TopologicalTorsionReadError,match="Count bounds are empty") as caught:g.fingerprints([None,m],num_threads=2)
    assert caught.value.kind=="Generator" and caught.value.__cause__ is not None
    assert sum(g.sparse_counts([m],num_threads=2)[0].nonzero_elements().values())==1
    s.count_bounds=[1,2,4,8];assert g.fingerprints([m],num_threads=2)[0].on_bits()==[60]
    with pytest.raises(ck.TopologicalTorsionReadError,match="INT_MIN"):g.counts([m],num_threads=-2147483648)
    for name,value,error in [("fp_size",-1,OverflowError),("torsion_atom_count",4294967296,OverflowError),("bits_per_feature",1.5,TypeError),("count_bounds",[-1],OverflowError)]:
        old=getattr(s,name)
        with pytest.raises(error):setattr(s,name,value)
        assert getattr(s,name)==old


def test_source_reusable_output_is_unique_and_empty_roots_are_distinct():
    m=ck.Molecule.from_smiles("CCCCO");g=ck.TopologicalTorsionFingerprintGenerator();output=ck.FingerprintAdditionalOutput();output.allocate_bit_paths();output.allocate_atom_to_bits();output.allocate_atom_counts()
    m.fingerprint_topological_torsion_sparse_count_with_generator(g,output=output)
    assert set(output.bit_paths())=={4437590048,12893306913}
    copied=output.bit_paths();copied.clear();assert output.bit_paths()
    m.fingerprint_topological_torsion_sparse_count_with_generator(g,params=ck.TopologicalTorsionCallParams(from_atoms=[]),output=output)
    assert output.bit_paths()=={} and output.atom_counts()==[0]*5
    assert output.atom_to_bits()==[[]]*5
    assert ck.TopologicalTorsionCallParams().from_atoms is None
    assert ck.TopologicalTorsionCallParams(from_atoms=[]).from_atoms==[]


def test_generated_reusable_constructor_bound_settings_and_bulk_protocols():
    tree=ast.parse((Path(__file__).resolve().parents[1]/"cosmolkit.pyi").read_text());classes={n.name:n for n in tree.body if isinstance(n,ast.ClassDef)}
    for name in ["TopologicalTorsionFingerprintGenerator","TopologicalTorsionSettings","TopologicalTorsionCallParams"]:assert name in classes
    generator={n.name:n for n in classes["TopologicalTorsionFingerprintGenerator"].body if isinstance(n,ast.FunctionDef)}
    assert set(generator)=={"__new__","new","from_json","settings","info_string","to_json","__repr__","fingerprints","sparse_fingerprints","counts","sparse_counts"}
    settings={}
    for n in classes["TopologicalTorsionSettings"].body:
        if isinstance(n,ast.FunctionDef): settings.setdefault(n.name, []).append(n)
    fields = {"torsion_atom_count", "only_shortest_paths", "include_chirality", "count_simulation", "fp_size", "bits_per_feature", "count_bounds"}
    assert set(settings) == fields | {"set_" + name for name in fields} | {"params", "__repr__"}
    for name in ["torsion_atom_count","only_shortest_paths","include_chirality","count_simulation","fp_size","bits_per_feature","count_bounds"]:
        assert len(settings[name]) == 2
        getter,setter=settings[name]
        assert [ast.unparse(d) for d in getter.decorator_list] == ["property"]
        assert [ast.unparse(d) for d in setter.decorator_list] == [f"{name}.setter"]
        assert [a.arg for a in getter.args.args] == ["self"]
        assert [a.arg for a in setter.args.args] == ["self","value"]
        expected_type = "builtins.list[builtins.int]" if name=="count_bounds" else ("builtins.bool" if name in ("only_shortest_paths","include_chirality","count_simulation") else "builtins.int")
        assert ast.unparse(getter.returns) == expected_type
        assert ast.unparse(setter.args.args[1].annotation) == ("typing.Sequence[builtins.int]" if name == "count_bounds" else expected_type)
        assert ast.unparse(setter.returns) == "None"
