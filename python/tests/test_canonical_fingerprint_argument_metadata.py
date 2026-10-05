"""Proposed canonical argument metadata projections, pinned RDKit351f/modern d892."""
import ast
import gc
import json
from pathlib import Path

import cosmolkit as ck
import pytest

TYPES = [ck.MorganParams, ck.AtomPairParams, ck.TopologicalTorsionParams]


@pytest.mark.parametrize("kind", TYPES)
def test_immutable_update_owns_values_and_survives_receiver_lifetime(kind):
    original = kind()
    updated = original.with_json('{"fpSize":"1000","includeChirality":"true","countBounds":[1,"3",7]}')
    assert type(updated) is kind
    assert updated is not original
    assert original.fp_size == 2048
    assert original.include_chirality is False
    assert original.count_bounds == [1, 2, 4, 8]
    bounds = updated.count_bounds
    bounds.append(99)
    assert updated.count_bounds == [1, 3, 7]
    with pytest.raises(AttributeError):
        updated.fp_size = 12
    del original
    gc.collect()
    assert updated.fp_size == 1000
    assert updated.include_chirality is True
    assert kind().with_json(updated.to_json()).to_json() == updated.to_json()


@pytest.mark.parametrize("kind", TYPES)
def test_json_leaf_types_empty_and_missing_bounds(kind):
    original = kind()
    encoded = json.loads(original.to_json())
    assert encoded["fpSize"] == "2048"
    assert encoded["countBounds"] == ["1", "2", "4", "8"]
    for text in ('{}', '{"countBounds":[]}', '{"countBounds":""}'):
        result = original.with_json(text)
        assert result.count_bounds == []
        assert json.loads(result.to_json())["countBounds"] == ""
    for text in ('', ' \n\t'):
        assert original.with_json(text).to_json() == original.to_json()
    assert original.count_bounds == [1, 2, 4, 8]


@pytest.mark.parametrize("kind", TYPES)
@pytest.mark.parametrize("text,expected_kind", [
    ('{', 'Parse'),
    ('[1]', 'SourceSuccess'),
    ('{"fpSize":-1}', 'SourceSuccess'),
    ('{"includeChirality":"invalid"}', 'SourceSuccess'),
    ('{"countBounds":[1,-1]}', 'SourceSuccess'),
])
def test_typed_errors_retain_receiver_and_source_cause(kind, text, expected_kind):
    original = kind()
    before = original.to_json()
    # Every original input is retained. Actual pinned source succeeds for the
    # previous false Invalid cases; only parser failures retain Parse here.
    if expected_kind == 'SourceSuccess':
        updated = original.with_json(text)
        assert updated is not original
        assert type(updated) is kind
        assert updated.include_chirality is False
        assert updated.fp_size == (4294967295 if 'fpSize' in text else 2048)
        assert updated.count_bounds == ([1, 4294967295] if 'countBounds' in text else [])
    else:
        with pytest.raises(ck.FingerprintJsonError) as caught:
            original.with_json(text)
        assert isinstance(caught.value, ValueError)
        assert caught.value.domain == 'fingerprints'
        assert caught.value.kind == expected_kind
        assert str(caught.value)
        assert isinstance(caught.value.__cause__, ValueError)
        assert str(caught.value.__cause__) == str(caught.value)
    assert original.to_json() == before
    with pytest.raises(TypeError):
        original.with_json(None)


def test_derived_fields_and_source_omitted_morgan_policies():
    m = ck.MorganParams(use_bond_types=False, include_ring_membership=False,
                        include_redundant_environments=True)
    result = m.with_json('{"radius":"4294967295","onlyNonzeroInvariants":"1","useBondTypes":true,"includeRedundantEnvironments":false}')
    assert result.radius == 4294967295
    assert result.only_nonzero_invariants is True
    assert result.use_bond_types is False
    assert result.include_ring_membership is False
    assert result.include_redundant_environments is True
    assert result.info_string() == 'MorganArguments onlyNonzeroInvariants=1 radius=4294967295'
    assert not {'useBondTypes', 'includeRingMembership', 'includeRedundantEnvironments'} & json.loads(result.to_json()).keys()
    ap = ck.AtomPairParams().with_json('{"use2D":false,"minDistance":30,"maxDistance":1,"numBitsPerFeature":0}')
    assert ap.info_string() == 'AtomPairArguments use2D=0 minDistance=30 maxDistance=1'
    assert ap.bits_per_feature == 0
    tt = ck.TopologicalTorsionParams().with_json('{"torsionAtomCount":8,"onlyShortestPaths":"true","numBitsPerFeature":0}')
    assert tt.info_string() == 'TopologicalTorsionArguments torsionAtomCount=8 onlyShortestPaths=1'
    assert tt.bits_per_feature == 0


def test_original_generated_stub_projects_value_update_and_structured_error():
    stub = Path(__file__).resolve().parents[1] / 'cosmolkit.pyi'
    classes = {node.name: node for node in ast.parse(stub.read_text()).body
               if isinstance(node, ast.ClassDef)}
    for name in ('MorganParams', 'AtomPairParams', 'TopologicalTorsionParams'):
        methods = {node.name: node for node in classes[name].body
                   if isinstance(node, ast.FunctionDef)}
        assert {'info_string', 'to_json', 'with_json'} <= methods.keys()
        method = methods['with_json']
        assert [arg.arg for arg in method.args.args] == ['self', 'json']
        assert ast.unparse(method.returns) == name
        assert not any(isinstance(decorator, ast.Attribute) and decorator.attr == 'setter'
                       for node in methods.values() for decorator in node.decorator_list)
    error = classes['FingerprintJsonError']
    assert ast.unparse(error.bases[0]) == 'builtins.ValueError'
    assert {node.target.id for node in error.body if isinstance(node, ast.AnnAssign)} >= {'domain', 'kind'}


@pytest.mark.parametrize("kind", TYPES)
@pytest.mark.parametrize("text,size,bounds", [
    ('{"includeChirality":"invalid"}', 2048, []),
    ('{"fpSize":-1}', 4294967295, []),
    ('{"fpSize":"invalid"}', 2048, []),
    ('{"fpSize":4294967296}', 2048, []),
    ('{"countBounds":{"x":1,"y":3}}', 2048, [1,3]),
    ('{"countBounds":"bad"}', 2048, []),
    ('{"fpSize":1,"fpSize":2}', 1, []),
])
def test_independent_source_seven_cases_first_duplicate_and_all_child_iteration(kind, text, size, bounds):
    original = kind()
    updated = original.with_json(text)
    assert updated.fp_size == size
    assert updated.include_chirality is False
    assert updated.count_bounds == bounds
    assert original.fp_size == 2048
    assert original.count_bounds == [1,2,4,8]
    same = kind().with_json(updated.to_json())
    assert same.to_json() == updated.to_json()


@pytest.mark.parametrize("kind", TYPES)
@pytest.mark.parametrize("text", [
    '{"countBounds":[1,"invalid"]}',
    '{"countBounds":[{}]}',
    '{"countBounds":[4294967296]}',
])
def test_strict_bounds_child_errors_retain_receiver_and_structured_class(kind, text):
    original = kind()
    before = original.to_json()
    with pytest.raises(ck.FingerprintJsonError) as caught:
        original.with_json(text)
    error = caught.value
    assert isinstance(error, ValueError)
    assert error.domain == 'fingerprints'
    assert error.kind == 'Invalid'
    assert error.__cause__ is None
    assert str(error)
    assert original.to_json() == before


@pytest.mark.parametrize("kind", TYPES)
def test_source_lexical_stream_retry_and_current_default_in_binding(kind):
    original = kind()
    for text in ('{"fpSize":1e0,"includeChirality":1e0}',
                 '{"fpSize":1.0,"includeChirality":1.0}'):
        updated = original.with_json(text)
        assert updated.fp_size == 2048
        assert updated.include_chirality is False
    for text, expected in [('2true',True),('2 true',True),('-1true',True),
                           ('+true',True),('+false',False),('1true',False),
                           ('+1',True),('01',True),('-0',False)]:
        updated = original.with_json(json.dumps({'includeChirality':text}))
        assert updated.include_chirality is expected
    current = kind(fp_size=19, include_chirality=True)
    updated = current.with_json('{"fpSize":"bad","includeChirality":"bad","countBounds":{"z":3,"a":1,"z":7}}')
    assert updated.fp_size == 19
    assert updated.include_chirality is True
    assert updated.count_bounds == [3,1,7]
    assert current.count_bounds == [1,2,4,8]
