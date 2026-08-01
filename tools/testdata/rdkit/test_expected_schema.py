from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from _expected_schema import SCHEMAS, validate_jsonl_output
from generate_all import GENERATOR_SPECS


class ExpectedSchemaTests(unittest.TestCase):
    def test_every_generator_output_has_exactly_one_schema(self) -> None:
        self.assertEqual({item.output for item in GENERATOR_SPECS}, set(SCHEMAS))

    def test_validation_rejects_missing_unknown_and_wrong_type_fields(self) -> None:
        valid = {
            "error": None,
            "inchi": "InChI=1S/CH4/h1H4",
            "rdkit_ok": True,
            "row": 1,
            "smiles": "C",
        }
        invalid_records = [
            {key: value for key, value in valid.items() if key != "row"},
            {**valid, "unknown": 1},
            {**valid, "row": "1"},
        ]
        for record in invalid_records:
            with self.subTest(record=record):
                with tempfile.TemporaryDirectory() as directory:
                    path = Path(directory) / "inchi.jsonl"
                    path.write_text(json.dumps(record) + "\n", encoding="utf-8")
                    with self.assertRaises(ValueError):
                        validate_jsonl_output(path, "inchi.jsonl")

    def test_validation_returns_checksum_and_noncomment_record_count(self) -> None:
        record = {
            "error": None,
            "inchi": "InChI=1S/CH4/h1H4",
            "rdkit_ok": True,
            "row": 1,
            "smiles": "C",
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "inchi.jsonl"
            path.write_text(
                "# generated\n" + json.dumps(record) + "\n\n", encoding="utf-8"
            )
            checksum, records = validate_jsonl_output(path, "inchi.jsonl")
            self.assertEqual(records, 1)
            self.assertEqual(len(checksum), 64)


if __name__ == "__main__":
    unittest.main()
