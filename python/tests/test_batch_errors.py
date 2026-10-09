"""Local regressions for Python batch error projections."""

import json
from pathlib import Path

import pytest
import cosmolkit as ck


@pytest.mark.parametrize("mode", ["keep", "raise"])
def test_batch_error_as_dict_preserves_field_types_and_is_detached(mode: str):
    inputs = ["CCO", "C1CC", "N1"]
    if mode == "keep":
        errors = ck.MoleculeBatch.from_smiles_list(inputs, errors="keep").errors()
    else:
        with pytest.raises(ck.BatchValidationError) as caught:
            _ = ck.MoleculeBatch.from_smiles_list(inputs, errors="raise")
        errors = caught.value.errors()

    assert [error.index() for error in errors] == [1, 2]
    for error in errors:
        expected = {
            "index": error.index(),
            "operation": error.operation(),
            "message": error.message(),
        }
        result = error.as_dict()
        assert type(result) is dict
        assert result == expected
        assert type(result["index"]) is int
        assert type(result["operation"]) is str
        assert type(result["message"]) is str
        index = result["index"]
        assert isinstance(index, int)
        assert inputs[index] in ("C1CC", "N1")
        assert json.loads(json.dumps(result)) == expected

        result["index"] = -1
        result["message"] = "changed by caller"
        assert error.as_dict() == expected


@pytest.mark.parametrize("kind", ["sdf", "sdf_files", "images"])
@pytest.mark.parametrize("inputs", [["CCO", "C1CC", "O"], ["C1CC", "N1"], []])
def test_export_report_counts_all_errors_without_a_skipped_state(
    kind: str, inputs: list[str], tmp_path: Path
):
    batch = ck.MoleculeBatch.from_smiles_list(inputs, errors="keep")
    original = [error.as_dict() for error in batch.errors()]
    if kind == "sdf":
        report = batch.write_sdf(str(tmp_path / "out.sdf"), errors="keep")
    elif kind == "sdf_files":
        report = batch.write_sdf_files(str(tmp_path / "files"), errors="keep")
    else:
        report = batch.write_images(
            str(tmp_path / "images"), format="svg", execution=ck.BatchParams(errors="keep")
        )
    assert report.total() == len(inputs)
    assert report.success() == batch.valid_count()
    assert report.failed() == len(original)
    assert report.total() == report.success() + report.failed()
    assert [error.as_dict() for error in report.errors()] == original
    assert report.failed() == len(report.errors())
    assert not hasattr(report, "skipped")
    assert "skipped" not in repr(report)
    path = tmp_path / "counts.json"
    report.write_report(str(path))
    assert json.loads(path.read_text()) == {"written": report.success(), "failed": report.failed()}
    assert [error.as_dict() for error in batch.errors()] == original


def test_export_report_keeps_input_and_write_errors_in_order(tmp_path: Path):
    batch = ck.MoleculeBatch.from_smiles_list(["CCO", "C1CC", "O"], errors="keep")
    original = batch.errors()[0].as_dict()
    (tmp_path / "blocked.sdf").mkdir()
    report = batch.write_sdf_files(
        str(tmp_path), errors="keep", filenames=["blocked.sdf", None, "water.sdf"]
    )
    assert (report.total(), report.success(), report.failed()) == (3, 1, 2)
    errors = report.errors()
    assert [error.index() for error in errors] == [0, 1]
    assert errors[0].operation() == "batch.write_sdf_files"
    assert errors[0].message()
    assert errors[1].as_dict() == original
    assert (tmp_path / "water.sdf").is_file()
