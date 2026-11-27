"""Constructor input ownership, live configuration and report snapshots.

The compiled registry checker owns complete field/signature agreement. Keep a
fixed set of scalar protocol samples here; do not grow a second field inventory.
"""
import inspect
import json
import cosmolkit as ck
import pytest

DEFAULTS = {
    'max_iterations':0, 'use_basic_knowledge':False, 'optimizer_force_tol':1e-3,
}
DISTINCT = {
    'max_iterations':19, 'use_basic_knowledge':True, 'optimizer_force_tol':0.125,
}


def test_constructor_preserves_explicit_source_defaults():
    signature = inspect.signature(ck.EmbedParams)
    params = ck.EmbedParams()
    assert params.failures == []
    for field, expected in DEFAULTS.items():
        assert getattr(params, field) == expected, field
        assert signature.parameters[field].default is Ellipsis, field
        assert signature.parameters[field].kind == inspect.Parameter.KEYWORD_ONLY
    with pytest.raises(TypeError):
        ck.EmbedParams(failures=[1])


@pytest.mark.parametrize('field', DEFAULTS)
def test_constructor_scalar_protocol_preserves_other_defaults(field):
    params = ck.EmbedParams(**{field:DISTINCT[field]})
    assert getattr(params, field) == DISTINCT[field]
    for other, expected in DEFAULTS.items():
        if other != field:
            assert getattr(params, other) == expected, other


def test_constructor_maps_are_detached_and_configuration_views_are_live():
    maps={0:[1.,2.,3.]};cpci={(0,1):0.125}
    params=ck.EmbedParams(coord_map=maps,cpci=cpci)
    maps[0][0]=9.;cpci[(0,1)]=9.
    assert params.coord_map == {0:(1.,2.,3.)}
    assert params.cpci == {(0,1):0.125}
    changed=params.with_json('{"randomSeed":42}')
    assert params.random_seed == -1 and changed.random_seed == 42
    assert changed.coord_map == params.coord_map and changed.cpci == params.cpci
    restored=changed.with_json(changed.to_json())
    assert json.loads(restored.to_json()) == json.loads(changed.to_json())
    restored.coord_map = {0: (99.,2.,3.)}
    assert restored.coord_map == {0:(99.,2.,3.)}
    assert changed.coord_map == params.coord_map == {0:(1.,2.,3.)}
    failures=restored.failures;failures.append(1)
    assert restored.failures == []
