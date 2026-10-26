"""Configuration values and side-effect-free display; no custom repr contract."""

from typing_extensions import override
from pathlib import Path

import pytest
import cosmolkit as ck


@pytest.mark.parametrize("params,expected", [
    (ck.BatchParams(), dict(errors=None, n_jobs=None, progress_bar=None)),
    (ck.BatchQueryParams(), dict(n_jobs=None, progress_bar=None, progress_callback=None)),
    (ck.BatchExportParams(), dict(format=ck.SdfFormat.V2000, errors=None, n_jobs=None, progress_bar=None)),
    (ck.BatchImageParams(), dict(format="png", width=300, height=300, filenames=None, report_path=None)),
])
def test_batch_configuration_repr_displays_all_defaults(
    params: ck.BatchParams | ck.BatchQueryParams | ck.BatchExportParams | ck.BatchImageParams,
    expected: dict,
):
    assert {name: getattr(params, name) for name in expected} == expected
    assert type(params).__name__ in repr(params)
    assert {name: getattr(params, name) for name in expected} == expected


def test_batch_configuration_repr_updates_after_assignment_without_mutation():
    params = ck.BatchParams(errors="keep", n_jobs=2, progress_bar=True)
    assert (params.errors, params.n_jobs, params.progress_bar) == (ck.BatchErrorMode.KEEP, 2, True)
    params.n_jobs = 3
    assert type(params).__name__ in repr(params)
    assert params.n_jobs == 3
    assert params.progress_bar is True
    with pytest.raises(OverflowError):
        params.n_jobs = -1
    assert (params.errors, params.n_jobs, params.progress_bar) == (ck.BatchErrorMode.KEEP, 3, True)


def test_export_repr_uses_python_values_and_reflects_every_field():
    params = ck.BatchExportParams(format="v3000", errors="keep", n_jobs=2, progress_bar=True)
    assert (params.format, params.errors, params.n_jobs, params.progress_bar) == (ck.SdfFormat.V3000, ck.BatchErrorMode.KEEP, 2, True)
    params.format = "v2000"
    params.errors = "raise"
    params.n_jobs = None
    params.progress_bar = False
    assert (params.format, params.errors, params.n_jobs, params.progress_bar) == (ck.SdfFormat.V2000, ck.BatchErrorMode.RAISE, None, False)


def test_image_repr_includes_nested_configuration_and_escapes_text():
    filenames = ["water's.sdf", None, "line\nbreak"]
    path = 'reports/"counts".json'
    params = ck.BatchImageParams(
        format="svg", width=100, height=200,
        execution=ck.BatchParams(errors="keep", n_jobs=2),
        filenames=filenames, report_path=path,
    )
    assert (params.format, params.width, params.height) == ("svg", 100, 200)
    assert type(params).__name__ in repr(params)
    assert params.execution.n_jobs == 2
    assert params.execution.errors == ck.BatchErrorMode.KEEP
    assert params.execution.progress_bar is None
    assert params.filenames == tuple(filenames)
    assert params.report_path == Path(path)
    params.width = 101
    assert params.width == 101
    assert params.height == 200


class Callback:
    def __call__(self) -> None:
        raise AssertionError("repr must not execute the progress callback")

    @override
    def __repr__(self) -> str:
        return "Callback()"


def test_query_repr_displays_callback_without_invoking_or_replacing_it():
    callback = Callback()
    params = ck.BatchQueryParams(n_jobs=2, progress_bar=True, progress_callback=callback)
    assert "BatchQueryParams" in repr(params)
    assert params.n_jobs == 2 and params.progress_bar is True
    assert params.progress_callback is callback
    params.progress_bar = False
    assert "BatchQueryParams" in repr(params)
    assert params.n_jobs == 2 and params.progress_bar is False
    assert params.progress_callback is callback


def test_configuration_display_does_not_execute_callback_repr():
    error = RuntimeError("callback repr failed")

    class BadCallback(Callback):
        @override
        def __repr__(self) -> str:
            raise error

    callback = BadCallback()
    params = ck.BatchQueryParams(n_jobs=2, progress_callback=callback)
    assert "BatchQueryParams" in repr(params)
    with pytest.raises(RuntimeError) as caught:
        repr(callback)
    assert caught.value is error
    assert params.n_jobs == 2
    assert params.progress_callback is callback
