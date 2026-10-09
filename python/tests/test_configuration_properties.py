"""Small native-binding regressions, not a second corpus parity suite."""
import json
import gc
import weakref

import pytest

import cosmolkit as ck


def test_configuration_assignment_reaches_the_actual_operation():
    params = ck.SubstructMatchParams()
    params.max_matches = 1
    params.uniquify = False
    molecule = ck.Molecule.from_smiles("CCC")
    query = ck.parse_smarts("C")
    assert len(molecule.substruct_matches_with_params(query, params)) == 1
    with pytest.raises((TypeError, OverflowError)):
        params.max_matches = -1
    assert params.max_matches == 1
    assert params.uniquify is False


def test_assignment_preserves_nested_configuration_and_owned_lists():
    generator = ck.MorganParams(radius=2, include_chirality=True)
    params = ck.MorganFingerprintParams(generator=generator, from_atoms=[0])
    params.from_atoms = [1, 2]
    assert params.from_atoms == [1, 2]
    assert params.generator.radius == 2
    assert params.generator.include_chirality is True
    params.generator = ck.MorganParams(radius=4)
    assert params.generator.radius == 4
    assert params.from_atoms == [1, 2]


def test_nested_edits_reach_native_calls_without_aliasing_constructor_inputs():
    original = ck.MorganParams(radius=3)
    roots = [0]
    params = ck.MorganFingerprintParams(generator=original, from_atoms=roots)
    generator = params.generator
    generator.radius = 0
    params.from_atoms.append(1)
    assert params.generator is generator
    assert params.generator.radius == 0
    assert params.from_atoms == [0, 1]
    assert original.radius == 3 and roots == [0]
    molecule = ck.Molecule.from_smiles("CCO")
    explicit = ck.MorganFingerprintParams(generator=ck.MorganParams(radius=0), from_atoms=[0, 1])
    assert molecule.fingerprint_morgan_with_params(params, None).on_bits() == molecule.fingerprint_morgan_with_params(explicit, None).on_bits()
    params.generator = ck.MorganParams(radius=2)
    assert generator.radius == 2
    assert ck.MorganFingerprintParams(generator=generator).generator.radius == 2
    params.generator.count_bounds.append(16)
    assert params.generator.count_bounds == [1, 2, 4, 8, 16]


def test_nested_and_list_conversion_failures_are_atomic():
    params = ck.MorganFingerprintParams(from_atoms=[0])
    child, roots = params.generator, params.from_atoms
    with pytest.raises((OverflowError, TypeError)):
        child.radius = -1
    assert child.radius == params.generator.radius == 3
    with pytest.raises((OverflowError, TypeError)):
        roots.extend([1, -1])
    assert roots == params.from_atoms == [0]
    with pytest.raises((OverflowError, TypeError)):
        params.generator.count_bounds.append(-1)
    assert child.count_bounds == [1, 2, 4, 8]
    params.from_atoms = None
    roots.append(2)
    assert params.from_atoms is None and roots == [0, 2]


@pytest.mark.parametrize("method,args,expected", [
    ("append", (3,), [0, 1, 2, 3]),
    ("extend", ([3, 4],), [0, 1, 2, 3, 4]),
    ("insert", (1, 4), [0, 4, 1, 2]),
    ("pop", (), [0, 1]),
    ("remove", (1,), [0, 2]),
    ("reverse", (), [2, 1, 0]),
    ("sort", (), [0, 1, 2]),
    ("clear", (), []),
    ("__setitem__", (slice(1, 3), [4]), [0, 4]),
    ("__delitem__", (0,), [1, 2]),
    ("__iadd__", ([3],), [0, 1, 2, 3]),
    ("__imul__", (2,), [0, 1, 2, 0, 1, 2]),
])
def test_list_mutations_write_through(method, args, expected):
    params = ck.MorganFingerprintParams(from_atoms=[0, 1, 2])
    view = params.from_atoms
    getattr(view, method)(*args)
    assert view == params.from_atoms == expected


def test_nested_alignment_list_edits_are_validated_and_written_back():
    params = ck.BestAlignmentParameters(atom_maps=[[ck.AlignmentAtomMap(0, 0)]])
    params.atom_maps[0].append(ck.AlignmentAtomMap(1, 1))
    assert len(params.atom_maps[0]) == 2
    with pytest.raises(TypeError):
        params.atom_maps[0].append(42)
    assert len(params.atom_maps[0]) == 2
    # Passing the view to a native constructor reads the committed contents.
    copied = ck.BestAlignmentParameters(atom_maps=params.atom_maps)
    assert len(copied.atom_maps[0]) == 2


def test_manually_defined_setters_also_support_live_configuration_views():
    params = ck.AlignmentParameters(atom_map=[ck.AlignmentAtomMap(0, 0)])
    params.atom_map.append(ck.AlignmentAtomMap(1, 1))
    assert len(ck.AlignmentParameters(atom_map=params.atom_map).atom_map) == 2
    callback = lambda *args: False
    tautomer = ck.TautomerParams.v1()
    tautomer.callback = callback
    transforms = tautomer.transform_count()
    tautomer.score_params.terms = []
    tautomer.score_params.terms.append(ck.TautomerScoreTerm("carbon", "[#6]", 3))
    assert [term.score() for term in tautomer.score_params.terms] == [3]
    assert tautomer.callback is callback and tautomer.transform_count() == transforms
    with pytest.raises(TypeError):
        tautomer.score_params.terms.append(42)
    assert [term.score() for term in tautomer.score_params.terms] == [3]


def test_mapping_and_nested_coordinate_edits_preserve_preset_state():
    params = ck.EmbedParams.etkdg_v3()
    before = json.loads(params.to_json())
    params.coord_map = {0: [1., 2., 3.]}
    params.coord_map[0][0] = 4.
    params.coord_map.update({1: [5., 6., 7.]})
    assert params.coord_map == {0: [4., 2., 3.], 1: [5., 6., 7.]}
    with pytest.raises((TypeError, ValueError)):
        params.coord_map[0].append(8.)
    assert params.coord_map[0] == [4., 2., 3.]
    assert params.coord_map.pop(1) == [5., 6., 7.]
    params.coord_map.clear()
    assert params.coord_map == {}
    after = json.loads(params.to_json())
    # Preset-specific state outside the changed field is not reconstructed.
    assert {key: value for key, value in before.items() if "coord" not in key.lower()} == {key: value for key, value in after.items() if "coord" not in key.lower()}


def test_views_do_not_keep_owners_alive_and_survive_owner_release():
    params = ck.MorganFingerprintParams(from_atoms=[0])
    reference = weakref.ref(params)
    child, roots = params.generator, params.from_atoms
    del params
    gc.collect()
    assert reference() is None
    child.radius = 2
    roots.append(1)
    assert child.radius == 2 and roots == [0, 1]
    child = ck.MorganFingerprintParams().generator
    child.radius = 0
    assert child.radius == 0


def test_setter_reuses_constructor_validation_and_is_atomic():
    callback = lambda *args: None
    params = ck.BatchQueryParams(n_jobs=3, progress_callback=callback)
    with pytest.raises(TypeError, match="callable"):
        ck.BatchQueryParams(progress_callback=42)
    with pytest.raises(TypeError, match="callable"):
        params.progress_callback = 42
    assert params.progress_callback is callback
    assert params.n_jobs == 3
    params.progress_bar = True
    assert params.progress_callback is callback
    assert params.progress_bar is True


def test_fixed_settings_and_fragment_required_inputs_remain_writable():
    params = ck.ForceFieldMinimizeParams()
    params.energy_tolerance = 1e-8
    assert params.energy_tolerance == 1e-8
    assert params.max_iterations == 200
    fragment = ck.FragmentSmilesWriteParams([0, 1])
    fragment.atoms = [0]
    assert fragment.atoms == [0]


def test_batch_scalar_slice_and_mask_have_distinct_runtime_results():
    batch = ck.MoleculeBatch.from_smiles_list(["C", "CC"])
    assert isinstance(batch[0], ck.Molecule)
    assert isinstance(batch[-1], ck.Molecule)
    assert isinstance(batch[:1], ck.MoleculeBatch)
    assert isinstance(batch[[0]], ck.MoleculeBatch)
    assert isinstance(batch[[True, False]], ck.MoleculeBatch)
    assert len(batch[:1]) == 1
    with pytest.raises(TypeError):
        batch[True]
    with pytest.raises(IndexError):
        batch[2]


def test_conformer_assignment_retains_all_other_preset_state():
    params = ck.EmbedParams.etkdg_v3()
    before = json.loads(params.to_json())
    params.random_seed = 37
    after = json.loads(params.to_json())
    changed = {key for key in before if before[key] != after[key]}
    assert changed == {"randomSeed"}
    assert params.random_seed == 37


def test_tautomer_assignment_retains_catalog_and_callback_identity():
    params = ck.TautomerParams.v1()
    count = params.transform_count()
    callback = lambda *args: False
    params.callback = callback
    params.max_tautomers = 9
    assert params.callback is callback
    assert params.transform_count() == count
    assert params.max_tautomers == 9
    with pytest.raises(TypeError, match="callable"):
        params.callback = 42
    assert params.callback is callback
    assert params.transform_count() == count


def test_bio_group_switch_remains_writable_and_preserves_output_formatting():
    params = ck.BioMmcifWriteParams(compact=True, auth_all=True, align_pairs=7)
    assert params.all_groups is True
    params.all_groups = False
    assert params.atoms is False
    assert params.entry is False
    assert params.software is False
    assert params.compact is True
    assert params.auth_all is True
    assert params.align_pairs == 7
    params.atoms = True
    assert params.all_groups is False
    params.all_groups = True
    assert params.entry is True
    assert params.software is True
    params.entry = None
    assert params.entry is True
    params.auth_all = None
    assert params.auth_all is False
    assert params.compact is True


@pytest.mark.parametrize("name", ["BestAlignmentParameters", "CoordinateRmsdParameters", "AllConformerRmsdParameters"])
def test_alignment_optional_maps_assignment_matches_constructor(name):
    cls = getattr(ck, name)
    params = cls(atom_maps=[[ck.AlignmentAtomMap(0, 0)]], max_matches=123)
    assert len(params.atom_maps) == 1
    params.atom_maps = None
    assert params.atom_maps == cls(atom_maps=None).atom_maps == []
    assert params.max_matches == 123
    with pytest.raises(TypeError):
        params.atom_maps = [42]
    assert params.atom_maps == []
    assert params.max_matches == 123
