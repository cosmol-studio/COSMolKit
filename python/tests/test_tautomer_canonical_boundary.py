"""Canonical TAU language-boundary conditions from the pinned RDKit wrapper.

Chemistry corpus parity belongs to parity-tests. These fixed cases cover Python
ownership, conversions, cancellation, exception identity and exact ordering.
"""
import gc
from typing import Any, cast

import pytest
import cosmolkit as ck


def acetone():
    return ck.Molecule.from_smiles("CC(C)=O")


def test_default_and_v1_parameter_factories():
    current, v1 = ck.TautomerParams(), ck.TautomerParams.v1()
    assert current.transform_count() == 37
    assert v1.transform_count() == 36
    for params in current, v1:
        assert params.max_tautomers() == params.max_transforms() == 1000
        assert params.remove_sp3_stereo()
        assert params.remove_bond_stereo()
        assert params.remove_isotopic_hydrogens()
        assert params.reassign_stereo()
        assert params.callback() is params.scorer() is None


@pytest.mark.parametrize("field,value", [
    ("max_tautomers", 7), ("max_transforms", 0),
    ("remove_sp3_stereo", False), ("remove_bond_stereo", False),
    ("remove_isotopic_hydrogens", False), ("reassign_stereo", False),
])
def test_parameter_updates_preserve_value_semantics(field, value):
    params = ck.TautomerParams()
    before = getattr(params, field)()
    updated = getattr(params, "with_" + field)(value)
    assert getattr(updated, field)() == value
    assert getattr(params, field)() == before
    getattr(params, "set_" + field)(value)
    assert getattr(params, field)() == value


@pytest.mark.parametrize("field", ["max_tautomers", "max_transforms"])
@pytest.mark.parametrize("value", [-1, 2**32])
def test_unsigned_limit_conversion_raises(field, value):
    with pytest.raises(OverflowError):
        ck.TautomerParams(**{field: value})


def test_rich_result_order_index_and_independent_ownership():
    source = acetone()
    result = source.enumerate_tautomers()
    assert result.status() == ck.TautomerEnumerationStatus.Completed
    assert len(result) == result.len() == 2
    assert not result.is_empty()
    assert result.canonical_smiles() == ["C=C(C)O", "CC(C)=O"]
    assert result.modified_atoms() == [0, 2, 3]
    assert result.modified_bonds() == [0, 1, 2]
    assert [m.to_smiles() for m in result] == result.canonical_smiles()
    assert [(k, m.to_smiles()) for k, m in result.entries()] == [
        ("C=C(C)O", "C=C(C)O"), ("CC(C)=O", "CC(C)=O")
    ]
    assert result[-1].to_smiles() == "CC(C)=O"
    assert result[-2].to_smiles() == "C=C(C)O"
    for index in [-3, 2]:
        with pytest.raises(IndexError, match="index out of bounds"):
            result[index]
    assert result.get(2) is None
    selected = result[0]
    iterator = result.iter()
    del result, source
    gc.collect()
    assert selected.to_smiles() == "C=C(C)O"
    assert [m.to_smiles() for m in iterator] == ["C=C(C)O", "CC(C)=O"]


def test_callback_pre_application_snapshot_survives_return_and_cancellation():
    saved = []
    def cancel(source, progress):
        saved.append((source, progress))
        assert isinstance(source, ck.TautomerMoleculeView)
        assert isinstance(progress, ck.TautomerProgress)
        return False
    params = ck.TautomerParams(callback=cancel)
    assert params.callback() is cancel
    source = acetone()
    result = source.enumerate_tautomers_with_params(params)
    assert result.status() == ck.TautomerEnumerationStatus.Canceled
    assert len(result) == len(saved) == 1
    assert source.to_smiles() == "CC(C)=O"
    del source, result, params
    gc.collect()
    source, progress = saved[0]
    assert source.to_smiles() == "CC(C)=O"
    assert source.to_owned().tautomer_score().total() == 5
    assert progress.status() == ck.TautomerEnumerationStatus.Completed
    assert progress.num_transforms() == 1
    assert len(progress) == 1 and not progress.is_empty()
    assert progress.modified_atoms() == progress.modified_bonds() == []
    assert [(key, value.to_smiles()) for key, value in progress.to_owned().entries()] == [
        ("CC(C)=O", "CC(C)=O")
    ]
    atoms, bonds = source.atoms(), source.bonds()
    assert all(isinstance(atom, ck.Atom) for atom in atoms)
    assert all(isinstance(bond, ck.Bond) for bond in bonds)
    assert [atom.atomic_number() for atom in atoms] == [6, 6, 6, 8]
    assert [atom.degree() for atom in atoms] == [1, 3, 1, 1]
    assert [row.total_hydrogens() for row in source.atom_metadata()] == [3, 0, 3, 0]
    assert source.atom(4) is None and source.bond(3) is None
    assert source.atom_degree(4) is None
    assert source.atom(0).total_hydrogens() == 3
    assert source.bond(2).stereo_code() == 0
    with pytest.raises(AttributeError):
        source.num_atoms = 0


@pytest.mark.parametrize("entry", ["enumerate_tautomers_with_params", "canonical_tautomer_with_params"])
def test_callback_exception_identity_and_source_atomicity(entry):
    error = RuntimeError("original callback error")
    def fail(source, progress):
        raise error
    source = acetone()
    params = ck.TautomerParams(callback=fail)
    with pytest.raises(RuntimeError) as caught:
        getattr(source, entry)(params)
    assert caught.value is error
    assert source.to_smiles() == "CC(C)=O"
    params.set_callback(None)
    assert source.canonical_tautomer_with_params(params).to_smiles() == "CC(C)=O"


@pytest.mark.parametrize("entry", ["source", "result"])
def test_custom_signed_scorer_exact_inputs_and_retention(entry):
    saved = []
    def score(molecule):
        saved.append(molecule)
        return 100 if molecule.to_smiles() == "C=C(C)O" else -100
    params = ck.TautomerParams(scorer=score)
    assert params.scorer() is score
    source = acetone()
    receiver = source if entry == "source" else source.enumerate_tautomers()
    selected = receiver.canonical_tautomer_with_params(params)
    assert selected.to_smiles() == "C=C(C)O"
    assert [m.to_smiles() for m in saved] == ["C=C(C)O", "CC(C)=O"]
    del receiver, selected, params
    gc.collect()
    assert saved[0].tautomer_score().total() == 1
    assert source.to_smiles() == "CC(C)=O"


@pytest.mark.parametrize("score,exception", [(2**31, OverflowError), ("wrong", TypeError)])
def test_scorer_return_conversion_propagates_original_error(score, exception):
    source = acetone()
    with pytest.raises(exception):
        source.canonical_tautomer_with_params(ck.TautomerParams(scorer=lambda _: score))
    assert source.to_smiles() == "CC(C)=O"


def test_scorer_exception_identity_is_preserved():
    error = LookupError("original scorer error")
    def fail(molecule):
        raise error
    with pytest.raises(LookupError) as caught:
        acetone().enumerate_tautomers().canonical_tautomer_with_params(ck.TautomerParams(scorer=fail))
    assert caught.value is error


def test_minimum_score_error_retains_structured_domain_and_cause():
    with pytest.raises(ck.OperationError) as caught:
        acetone().canonical_tautomer_with_params(ck.TautomerParams(scorer=lambda _: -(2**31)))
    assert caught.value.domain == "operation"
    assert caught.value.kind == "Tautomer"
    cause = caught.value.__cause__
    assert isinstance(cause, ck.TautomerRunError)
    assert cause.domain == "tautomer" and cause.kind == "NoCanonicalTautomer"


def test_single_result_skips_scorer_and_result_selection_does_not_enumerate():
    def unexpected(*args):
        raise AssertionError("unexpected re-enumeration/scoring")
    result = ck.Molecule.from_smiles("C").enumerate_tautomers()
    selected = result.canonical_tautomer_with_params(ck.TautomerParams(callback=unexpected, scorer=unexpected))
    assert selected.to_smiles() == "C"


def test_tie_break_uses_retained_lexical_order():
    result = acetone().enumerate_tautomers()
    assert result.canonical_tautomer_with_params(ck.TautomerParams(scorer=lambda _: -7)).to_smiles() == "C=C(C)O"
    assert result.canonical_tautomer().to_smiles() == "CC(C)=O"


def test_explicit_score_terms_are_detached_and_signed():
    source = acetone()
    default = ck.default_tautomer_score_terms()
    assert len(default) == 12
    term = ck.TautomerScoreTerm("negative carbonyl", "[C]=[O]", -9)
    assert term.name() == "negative carbonyl" and term.smarts() == "[C]=[O]" and term.score() == -9
    assert term == ck.TautomerScoreTerm("negative carbonyl", "[C]=[O]", -9)
    params = ck.TautomerScoreParams([term])
    score = source.tautomer_score_with_params(params)
    assert (score.ring(), score.substructure(), score.hetero_hydrogen(), score.total()) == (0, -9, 0, -9)
    assert source.tautomer_score_with_params(ck.TautomerScoreParams([])).total() == 0
    assert source.tautomer_score().total() == 5
    assert ck.TautomerScoreParams().terms is None


def test_catalog_file_error_is_structured(tmp_path):
    with pytest.raises(ck.TautomerCatalogError) as caught:
        ck.TautomerParams.from_transform_file(str(tmp_path / "missing"))
    assert caught.value.domain == "tautomer_catalog"
    assert caught.value.kind == "BadInputFile"
    assert caught.value.__cause__ is not None


@pytest.mark.parametrize("field", ["callback", "scorer"])
def test_noncallable_parameter_is_rejected_without_changing_configuration(field):
    with pytest.raises(TypeError, match="callable"):
        # Deliberately invalid runtime input, outside the public static type.
        ck.TautomerParams(**{field: cast(Any, 3)})
    params = ck.TautomerParams()
    with pytest.raises(TypeError, match="callable"):
        getattr(params, "set_" + field)(3)
    assert getattr(params, field)() is None


def test_iterable_selection_preserves_input_order_duplicates_and_first_tie():
    values = [ck.Molecule.from_smiles("C first"), ck.Molecule.from_smiles("C second"), ck.Molecule.from_smiles("CC third")]
    observed = []
    def score(value):
        observed.append(value.properties().name())
        return -9
    selected = ck.canonical_tautomer_from_molecules_with_params(iter(values), ck.TautomerParams(scorer=score))
    assert observed == ["first", "second", "third"]
    assert selected.to_smiles() == "C"
    assert selected.properties().name() == "first"
    assert [m.properties().name() for m in values] == observed
    result = acetone().enumerate_tautomers()
    assert ck.canonical_tautomer_from_molecules(result).to_smiles() == result.canonical_tautomer().to_smiles()


def test_iterable_selection_empty_and_wrong_types_remain_errors():
    with pytest.raises(ck.OperationError) as caught:
        ck.canonical_tautomer_from_molecules([])
    cause = caught.value.__cause__
    assert isinstance(cause, ck.TautomerRunError)
    assert cause.kind == "NoCanonicalTautomer"
    with pytest.raises(TypeError):
        ck.canonical_tautomer_from_molecules(cast(Any, [acetone(), 42]))
    with pytest.raises(TypeError):
        ck.canonical_tautomer_from_molecules(cast(Any, 42))


def test_iterable_generator_exception_keeps_identity():
    error = RuntimeError("original iteration error")
    def values():
        yield acetone()
        raise error
    with pytest.raises(RuntimeError) as caught:
        ck.canonical_tautomer_from_molecules(values())
    assert caught.value is error


def test_callback_property_values_are_owned_read_only_and_keep_computed_state():
    captured = []
    def callback(source, progress):
        captured.append(source.properties())
        return False
    source = ck.Molecule.from_smiles("CC(C)=O original")
    source.enumerate_tautomers_with_params(ck.TautomerParams(callback=callback))
    props = captured[0]
    assert isinstance(props, ck.MoleculeProperties)
    assert props.name() == "original"
    assert props.prop("_StereochemDone") == "1"
    assert props.is_prop_computed("_StereochemDone")
    assert "_StereochemDone" in props.computed_prop_names()
    assert props.sdf_data_fields() == props.sdf_property_lists() == []
    returned = props.props()
    returned["_StereochemDone"] = "changed copy"
    del source, captured
    gc.collect()
    assert props.prop("_StereochemDone") == "1"
    with pytest.raises(AttributeError):
        setattr(props, "name", "changed")


def test_atom_hybridization_is_the_shared_read_only_source_vocabulary():
    source = acetone()
    assert [atom.hybridization() for atom in source.atoms()] == [
        ck.Hybridization.SP3, ck.Hybridization.SP2,
        ck.Hybridization.SP3, ck.Hybridization.SP2,
    ]
    assert [atom.hybridization().name for atom in source.atoms()] == ["SP3", "SP2", "SP3", "SP2"]
    assert source.to_smiles() == "CC(C)=O"
