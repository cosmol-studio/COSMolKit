"""Fixed configuration-display checks; no chemistry or corpus preparation."""

from typing_extensions import override

import pytest
import cosmolkit as ck


@pytest.mark.parametrize("params,expected", [
    (ck.BatchParams(), "BatchParams(errors=None, n_jobs=None, progress_bar=None)"),
    (ck.BatchQueryParams(), "BatchQueryParams(n_jobs=None, progress_bar=None, progress_callback=None)"),
    (ck.BatchExportParams(), "BatchExportParams(format='v2000', errors=None, n_jobs=None, progress_bar=None)"),
    (ck.BatchImageParams(), "BatchImageParams(format='png', width=300, height=300, execution=BatchParams(errors=None, n_jobs=None, progress_bar=None), filenames=None, report_path=None)"),
])
def test_batch_configuration_repr_displays_all_defaults(
    params: ck.BatchParams | ck.BatchQueryParams | ck.BatchExportParams | ck.BatchImageParams,
    expected: str,
):
    assert repr(params) == expected
    assert str(params) == expected


def test_batch_configuration_repr_updates_after_assignment_without_mutation():
    params = ck.BatchParams(errors="keep", n_jobs=2, progress_bar=True)
    assert repr(params) == "BatchParams(errors=<BatchErrorMode.KEEP: 2>, n_jobs=2, progress_bar=True)"
    params.n_jobs = 3
    assert repr(params) == "BatchParams(errors=<BatchErrorMode.KEEP: 2>, n_jobs=3, progress_bar=True)"
    assert params == ck.BatchParams(errors="keep", n_jobs=3, progress_bar=True)
    assert params.n_jobs == 3
    assert params.progress_bar is True
    with pytest.raises(ValueError):
        params.n_jobs = 0
    assert repr(params) == "BatchParams(errors=<BatchErrorMode.KEEP: 2>, n_jobs=3, progress_bar=True)"


def test_export_repr_uses_python_values_and_reflects_every_field():
    params = ck.BatchExportParams(format="v3000", errors="keep", n_jobs=2, progress_bar=True)
    assert repr(params) == "BatchExportParams(format='v3000', errors=<BatchErrorMode.KEEP: 2>, n_jobs=2, progress_bar=True)"
    params.format = "v2000"
    params.errors = "raise"
    params.n_jobs = None
    params.progress_bar = False
    assert repr(params) == "BatchExportParams(format='v2000', errors=<BatchErrorMode.RAISE: 1>, n_jobs=None, progress_bar=False)"


def test_image_repr_includes_nested_configuration_and_escapes_text():
    filenames = ["water's.sdf", None, "line\nbreak"]
    path = 'reports/"counts".json'
    params = ck.BatchImageParams(
        format="svg", width=100, height=200,
        execution=ck.BatchParams(errors="keep", n_jobs=2),
        filenames=filenames, report_path=path,
    )
    expected = (
        "BatchImageParams(format='svg', width=100, height=200, "
        "execution=BatchParams(errors=<BatchErrorMode.KEEP: 2>, n_jobs=2, progress_bar=None), "
        f"filenames={filenames!r}, report_path={path!r})"
    )
    assert repr(params) == expected
    assert repr(params) == expected
    assert params.execution.n_jobs == 2
    assert params.execution == ck.BatchParams(errors="keep", n_jobs=2)
    assert params.filenames == filenames
    assert params.report_path == path
    params.width = 101
    assert repr(params) == expected.replace("width=100", "width=101")


class Callback:
    def __call__(self) -> None:
        raise AssertionError("repr must not execute the progress callback")

    @override
    def __repr__(self) -> str:
        return "Callback()"


def test_query_repr_displays_callback_without_invoking_or_replacing_it():
    callback = Callback()
    params = ck.BatchQueryParams(n_jobs=2, progress_bar=True, progress_callback=callback)
    assert repr(params) == "BatchQueryParams(n_jobs=2, progress_bar=True, progress_callback=Callback())"
    assert params.progress_callback is callback
    params.progress_bar = False
    assert repr(params) == "BatchQueryParams(n_jobs=2, progress_bar=False, progress_callback=Callback())"
    assert params.progress_callback is callback


def test_callback_repr_failure_is_not_silently_replaced():
    error = RuntimeError("callback repr failed")

    class BadCallback(Callback):
        @override
        def __repr__(self) -> str:
            raise error

    callback = BadCallback()
    params = ck.BatchQueryParams(n_jobs=2, progress_callback=callback)
    with pytest.raises(RuntimeError) as caught:
        _ = repr(params)
    assert caught.value is error
    assert params.n_jobs == 2
    assert params.progress_callback is callback
