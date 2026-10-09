"""Local batch policy regressions; no corpus or reference preparation."""

from pathlib import Path
from collections.abc import Callable
from typing import cast

import pytest
import cosmolkit as ck


@pytest.mark.parametrize("name", ["with_hydrogens", "without_hydrogens", "sanitize", "with_kekulized_bonds", "with_2d_coordinates"])
@pytest.mark.parametrize("form", ["default", "keywords", "params"])
def test_transforms_inherit_keep_in_every_call_form(name: str, form: str):
    batch = ck.MoleculeBatch.from_smiles_list(["CCO", "C1CC", "O"], errors="keep")
    before = batch.to_smiles_list()
    errors = [error.as_dict() for error in batch.errors()]
    transform = cast(Callable[..., ck.MoleculeBatch], getattr(batch, name))
    if form == "default":
        result = transform()
    elif form == "keywords":
        result = transform(n_jobs=1, errors=None)
    else:
        result = transform(params=ck.BatchParams(n_jobs=1))
    assert result.valid_mask() == [True, False, True]
    assert result.error_mode() == ck.BatchErrorMode.KEEP
    assert [error.as_dict() for error in result.errors()] == errors
    assert result.with_hydrogens().without_hydrogens().valid_mask() == result.valid_mask()
    with pytest.raises(ck.BatchValidationError) as caught:
        _ = transform(errors="raise")
    assert [error.as_dict() for error in caught.value.errors()] == errors
    assert batch.to_smiles_list() == before
    assert batch.error_mode() == ck.BatchErrorMode.KEEP


def test_new_errors_are_retained_and_explicit_override_changes_only_returned_chain():
    batch = ck.MoleculeBatch.from_smiles_list(["CCO", "C1CC", "O"], errors="keep")
    bad = ck.AddHsParams(only_on_atoms=[99])
    failed = batch.with_hydrogens(options=bad)
    assert failed.valid_mask() == [False, False, False]
    assert [error.index() for error in failed.errors()] == [0, 1, 2]
    assert failed.errors()[1].as_dict() == batch.errors()[0].as_dict()
    assert [error.operation() for error in failed.errors()] == [
        "batch.with_hydrogens", "batch.from_smiles_list", "batch.with_hydrogens"
    ]
    assert failed.without_hydrogens().invalid_count() == 3
    with pytest.raises(ck.BatchValidationError) as caught:
        _ = batch.with_hydrogens(options=bad, params=ck.BatchParams(errors="raise"))
    assert [error.as_dict() for error in caught.value.errors()] == [
        error.as_dict() for error in failed.errors()
    ]
    strict = batch.with_valid_records().with_hydrogens(errors="raise")
    assert strict.error_mode() == ck.BatchErrorMode.RAISE
    with pytest.raises(ck.BatchValidationError):
        _ = strict.with_hydrogens(options=bad)
    assert batch.error_mode() == ck.BatchErrorMode.KEEP
    assert batch.valid_mask() == [True, False, True]
    assert ck.MoleculeBatch.from_smiles_list(["CCO"]).error_mode() == ck.BatchErrorMode.RAISE


def test_selection_and_configuration_reset_preserve_policy():
    batch = ck.MoleculeBatch.from_smiles_list(["CCO", "C1CC"], errors="keep")
    for selected in [batch[:], batch[[1, 0]], batch[[True, True]], batch.with_valid_records()]:
        assert selected.error_mode() == ck.BatchErrorMode.KEEP
        assert selected.with_hydrogens().error_mode() == ck.BatchErrorMode.KEEP
    strict = ck.MoleculeBatch.from_smiles_list(["CCO"])
    assert strict[:].error_mode() == ck.BatchErrorMode.RAISE
    for params in [ck.BatchParams(), ck.BatchExportParams()]:
        assert cast(object, params.errors) is None
        params.errors = "raise"
        assert cast(object, params.errors) == ck.BatchErrorMode.RAISE
        params.errors = None
        assert cast(object, params.errors) is None
        with pytest.raises(ValueError):
            params.errors = "invalid"
        assert cast(object, params.errors) is None


@pytest.mark.parametrize("kind", ["sdf", "sdf_files", "images"])
def test_exports_inherit_keep_and_strict_override_keeps_all_error_details(kind: str, tmp_path: Path):
    batch = ck.MoleculeBatch.from_smiles_list(["CCO", "C1CC"], errors="keep")
    if kind == "sdf":
        report = batch.write_sdf(tmp_path / "kept.sdf", n_jobs=1)
        with pytest.raises(ck.BatchValidationError) as caught:
            _ = batch.write_sdf(tmp_path / "strict.sdf", errors="raise")
    elif kind == "sdf_files":
        report = batch.write_sdf_files(tmp_path / "kept", params=ck.BatchExportParams())
        with pytest.raises(ck.BatchValidationError) as caught:
            _ = batch.write_sdf_files(tmp_path / "strict", errors="raise")
    else:
        report = batch.write_images(str(tmp_path / "kept"), format="svg")
        with pytest.raises(ck.BatchValidationError) as caught:
            _ = batch.write_images(str(tmp_path / "strict"), format="svg", execution=ck.BatchParams(errors="raise"))
    assert (report.total(), report.success(), report.failed()) == (2, 1, 1)
    expected = [error.as_dict() for error in batch.errors()]
    assert [error.as_dict() for error in report.errors()] == expected
    assert [error.as_dict() for error in caught.value.errors()] == expected
