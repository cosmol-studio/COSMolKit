from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from _generate_builtin_fixture_migration_golden import build_records


class BuiltinFixtureGenerationTests(unittest.TestCase):
    def test_metadata_is_excluded_without_dropping_unknown_fixture_formats(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "README.md").write_text("fixture documentation\n", encoding="utf-8")
            (root / "source_manifest.jsonl").write_text("{}\n", encoding="utf-8")
            (root / "upstream_case.csv").write_bytes(b"case,value\n1,x\n")

            records = build_records(root)

        self.assertEqual(
            records,
            [
                {
                    "kind": "inventory",
                    "fixture": "upstream_case.csv",
                    "byte_len": 15,
                    "nonempty": True,
                }
            ],
        )


if __name__ == "__main__":
    unittest.main()
