"""Local filesystem-argument regressions; no prepared corpus or reference."""

import json
import os
from pathlib import Path
from typing import cast
from typing_extensions import override

import pytest
import cosmolkit as ck


class TextPath(os.PathLike[str]):
    def __init__(self, path: Path):
        self.path: Path = path

    @override
    def __fspath__(self) -> str:
        return str(self.path)


def path_input(path: Path, kind: str) -> str | os.PathLike[str]:
    if kind == "str":
        return str(path)
    if kind == "custom":
        return TextPath(path)
    return path


@pytest.mark.parametrize("kind", ["str", "pathlib", "custom"])
@pytest.mark.parametrize("form", ["default", "params", "explicit"])
def test_sdf_paths_support_all_call_forms(kind: str, form: str, tmp_path: Path):
    batch = ck.MoleculeBatch.from_smiles_list(["CCO", "O"])
    path = path_input(tmp_path / "分子.sdf", kind)
    counts = path_input(tmp_path / "计数.json", kind)
    directory = path_input(tmp_path / "分子", kind)
    options = ck.BatchExportParams()
    if form == "explicit":
        report = batch.write_sdf_with_params(path, options, counts)
        files = batch.write_sdf_files_with_params(directory, options, ["ethanol.sdf", "water.sdf"], None)
        result = ck.MoleculeBatch.read_sdf_with_params(path, ck.SdfReadParams())
    elif form == "params":
        report = batch.write_sdf(path, options, report_path=counts)
        files = batch.write_sdf_files(directory, options, filenames=["ethanol.sdf", "water.sdf"])
        result = ck.MoleculeBatch.read_sdf(path, ck.SdfReadParams())
    else:
        report = batch.write_sdf(path, report_path=counts)
        files = batch.write_sdf_files(directory, filenames=["ethanol.sdf", "water.sdf"])
        result = ck.MoleculeBatch.read_sdf(path)
    assert result.to_smiles_list() == ["CCO", "O"]
    assert (report.total(), report.success(), report.failed()) == (2, 2, 0)
    assert (files.total(), files.success(), files.failed()) == (2, 2, 0)
    assert json.loads((tmp_path / "计数.json").read_text()) == {"written": 2, "failed": 0}
    for filename, expected in [("ethanol.sdf", "CCO"), ("water.sdf", "O")]:
        assert ck.MoleculeBatch.read_sdf(path_input(tmp_path / "分子" / filename, kind)).to_smiles_list() == [expected]
    report.write_report(path_input(tmp_path / "report.csv", kind))
    assert (tmp_path / "report.csv").read_text().splitlines() == ["written,failed", "2,0"]
    assert batch.to_smiles_list() == ["CCO", "O"]


@pytest.mark.parametrize("kind", ["str", "pathlib", "custom"])
@pytest.mark.parametrize("configured", [False, True])
def test_sdf_dataset_and_reader_accept_path_protocol(kind: str, configured: bool, tmp_path: Path):
    path = tmp_path / "input.sdf"
    batch = ck.MoleculeBatch.from_smiles_list(["CCO", "O"])
    _ = batch.write_sdf(str(path))
    value = path_input(path, kind)
    if configured:
        dataset = ck.SdfDataset.open_with_params(value, ck.SdfReadParams())
        reader = ck.SdfReader.open_with_params(value, ck.SdfReadParams())
    else:
        dataset = ck.SdfDataset.open(value)
        reader = ck.SdfReader.open(value)
    assert len(dataset) == 2
    assert dataset[:].to_smiles_list() == ["CCO", "O"]
    assert reader.path() == path
    assert [item.to_smiles_list() for item in reader.batches(size=1)] == [["CCO"], ["O"]]


class BytesPath:
    def __fspath__(self) -> bytes:
        return b"not-a-text-path.sdf"


class InvalidPath:
    def __fspath__(self) -> int:
        return 42


@pytest.mark.parametrize("invalid", [b"file.sdf", BytesPath(), InvalidPath(), object()])
def test_non_text_paths_fail_before_export(invalid: object, tmp_path: Path):
    batch = ck.MoleculeBatch.from_smiles_list(["O"])
    output = tmp_path / "must-not-be-written.sdf"
    # Deliberately invalid argument types exercise runtime extraction.
    bad = cast(str, invalid)
    with pytest.raises(TypeError):
        _ = ck.MoleculeBatch.read_sdf(bad)
    with pytest.raises(TypeError):
        _ = batch.write_sdf(bad)
    with pytest.raises(TypeError):
        _ = batch.write_sdf_files(bad)
    with pytest.raises(TypeError):
        _ = batch.write_sdf(output, report_path=bad)
    assert not output.exists()
    assert batch.to_smiles_list() == ["O"]


def test_fspath_error_is_preserved_without_string_fallback(tmp_path: Path):
    error = RuntimeError("filesystem protocol failed")

    class BrokenPath(os.PathLike[str]):
        @override
        def __fspath__(self) -> str:
            raise error

        @override
        def __str__(self) -> str:
            return str(tmp_path / "must-not-be-written.sdf")

    batch = ck.MoleculeBatch.from_smiles_list(["O"])
    with pytest.raises(RuntimeError) as caught:
        _ = batch.write_sdf(BrokenPath())
    assert caught.value is error
    assert list(tmp_path.iterdir()) == []
