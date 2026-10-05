import json
from pathlib import Path

import pytest

import _generate_tautomer_special_regression as generator


def test_fixed_selection_reuses_existing_oracle_without_expanding_corpus(tmp_path, monkeypatch):
    fixture = Path(__file__).resolve().parents[3] / "testdata/tautomer/fixtures/rdkit/long_conjugated_cases.json"
    document = json.loads(fixture.read_text())
    calls = []
    monkeypatch.setattr(generator, "assert_rdkit_version", lambda: None)

    def build_record(case, profile, branches):
        calls.append((case, branches))
        return {"case_id": case["case_id"], "branches": branches}

    monkeypatch.setattr(generator, "build_record", build_record)
    output = tmp_path / "reference.jsonl"
    generator.generate(fixture, output)
    assert calls == [(document["cases"][0], ["default", "v1"])]
    assert json.loads(output.read_text()) == {
        "case_id": "smiles_5000:1399", "branches": ["default", "v1"]
    }
    assert output.read_bytes().endswith(b"\n")
    document["branches"][0]["max_transforms"] = 1
    invalid = tmp_path / "invalid.json"
    invalid.write_text(json.dumps(document))
    with pytest.raises(ValueError, match="pinned profile mismatch"):
        generator.generate(invalid, output)
    assert len(calls) == 1
