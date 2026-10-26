"""Typed configuration extraction and native batch execution validation."""

import pytest
import cosmolkit as ck


@pytest.mark.parametrize("params_type", [ck.BatchParams, ck.BatchQueryParams, ck.BatchExportParams])
def test_worker_defaults_and_valid_assignments_match_construction(
    params_type: type[ck.BatchParams] | type[ck.BatchQueryParams] | type[ck.BatchExportParams],
):
    params = params_type(progress_bar=True)
    assert params.n_jobs is None
    for value in (1, 2, None):
        params.n_jobs = value
        assert params.n_jobs == params_type(n_jobs=value).n_jobs
        assert params.progress_bar is True


@pytest.mark.parametrize("params_type", [ck.BatchParams, ck.BatchQueryParams, ck.BatchExportParams])
def test_zero_workers_are_rejected_at_native_execution(
    params_type: type[ck.BatchParams] | type[ck.BatchQueryParams] | type[ck.BatchExportParams],
    tmp_path,
):
    constructed = params_type(n_jobs=0)
    params = params_type(n_jobs=2, progress_bar=True)
    params.n_jobs = 0
    batch = ck.mols_from_smiles_list(["CCO"])
    before = batch.to_smiles()
    for configuration in (constructed, params):
        with pytest.raises(ck.BatchValidationError, match="n_jobs must be >= 1"):
            if params_type is ck.BatchParams:
                batch.with_hydrogens(params=configuration)
            elif params_type is ck.BatchQueryParams:
                batch.to_smiles(params=configuration)
            else:
                batch.write_sdf(tmp_path / "not-written.sdf", params=configuration)
        assert configuration.n_jobs == 0
    assert batch.to_smiles() == before
    assert not (tmp_path / "not-written.sdf").exists()
    assert params.progress_bar is True


@pytest.mark.parametrize("params_type", [ck.BatchParams, ck.BatchQueryParams, ck.BatchExportParams])
def test_wrong_worker_types_are_rejected_without_changing_configuration(
    params_type: type[ck.BatchParams] | type[ck.BatchQueryParams] | type[ck.BatchExportParams],
):
    params = params_type(n_jobs=2, progress_bar=True)
    for invalid in (-1, 1.5, "2"):
        # Deliberately invalid arguments must also fail at runtime.
        with pytest.raises((TypeError, OverflowError)):
            _ = params_type(n_jobs=invalid)  # pyright: ignore[reportArgumentType]
        with pytest.raises((TypeError, OverflowError)):
            params.n_jobs = invalid  # pyright: ignore[reportAttributeAccessIssue]
        assert params.n_jobs == 2
        assert params.progress_bar is True


def test_failed_query_worker_assignment_preserves_callback_identity():
    def callback() -> None:
        pass

    params = ck.BatchQueryParams(n_jobs=2, progress_bar=True, progress_callback=callback)
    with pytest.raises(OverflowError):
        params.n_jobs = -1
    assert params.n_jobs == 2
    assert params.progress_bar is True
    assert params.progress_callback is callback
