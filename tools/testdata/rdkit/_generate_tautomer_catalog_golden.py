#!/usr/bin/env python3
"""Generate exact current and V1 tautomer catalog records from pinned RDKit."""

from __future__ import annotations

import argparse
import ast
import json
import re
from pathlib import Path

from _tautomer_oracle import assert_rdkit_version


REPO_ROOT = Path(__file__).resolve().parents[3]
CATALOGS = {
    "current": REPO_ROOT / "third_party/rdkit/Code/GraphMol/MolStandardize/TautomerCatalog/tautomerTransforms.in",
    "v1": REPO_ROOT / "third_party/rdkit/Code/GraphMol/MolStandardize/TautomerCatalog/tautomerTransforms.v1.in",
}


def parse_catalog(name: str, path: Path) -> list[dict[str, object]]:
    source = path.read_text(encoding="utf-8")
    records: list[dict[str, object]] = []
    cursor = 0
    marker = "std::make_tuple("
    string_argument = re.compile(
        r'std::string\(((?:\s*"(?:\\.|[^"\\])*"\s*)*)\)', re.DOTALL
    )
    string_literal = re.compile(r'"(?:\\.|[^"\\])*"')
    while (start := source.find(marker, cursor)) >= 0:
        depth = 1
        position = start + len(marker)
        in_string = False
        escaped = False
        while depth and position < len(source):
            character = source[position]
            if in_string:
                if escaped:
                    escaped = False
                elif character == "\\":
                    escaped = True
                elif character == '"':
                    in_string = False
            elif character == '"':
                in_string = True
            elif character == "(":
                depth += 1
            elif character == ")":
                depth -= 1
            position += 1
        if depth:
            raise ValueError(f"{path}: unterminated std::make_tuple at byte {start}")
        expression = source[start:position]
        columns = []
        for match in string_argument.finditer(expression):
            columns.append(
                "".join(ast.literal_eval(token) for token in string_literal.findall(match.group(1)))
            )
        if len(columns) != 4:
            source_line = source.count("\n", 0, start) + 1
            raise ValueError(
                f"{path}:{source_line}: expected four std::string arguments, got {len(columns)}"
            )
        source_line = source.count("\n", 0, start) + 1
        records.append(
            {
                "schema_version": 1,
                "catalog": name,
                "index": len(records),
                "source_path": path.relative_to(REPO_ROOT / "third_party/rdkit").as_posix(),
                "source_line": source_line,
                "name": columns[0],
                "smarts": columns[1],
                "bonds": columns[2],
                "charges": columns[3],
            }
        )
        cursor = position
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert_rdkit_version()
    records = [record for name, path in CATALOGS.items() for record in parse_catalog(name, path)]
    counts = {name: sum(record["catalog"] == name for record in records) for name in CATALOGS}
    if counts != {"current": 37, "v1": 36}:
        raise ValueError(f"pinned tautomer catalog counts changed: {counts}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as output:
        for record in records:
            output.write(json.dumps(record, ensure_ascii=True, sort_keys=True, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
