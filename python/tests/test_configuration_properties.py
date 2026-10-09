"""Small native-binding regressions, not a second corpus parity suite."""
import json

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
