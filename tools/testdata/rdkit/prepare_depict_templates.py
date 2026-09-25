"""Reproduce RDKit 2026.03.1's generated depictor template rows.

Input: third_party/rdkit/Code/GraphMol/Depictor/TemplateSmarts.h
Output: crates/cosmolkit-depict/assets/default_templates.cxsmarts

The source notice and license are recorded beside the output. This is an
explicit preparation command, never run by ordinary tests.
"""

from __future__ import annotations

import ast
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "third_party/rdkit/Code/GraphMol/Depictor/TemplateSmarts.h"
TARGET = ROOT / "crates/cosmolkit-depict/assets/default_templates.cxsmarts"


def main() -> None:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    start = lines.index("const std::vector<std::string> TEMPLATE_SMARTS = {") + 1
    end = lines.index("};", start)
    rows: list[str] = []
    for source_line in lines[start:end]:
        entry = source_line.strip()
        if not entry.startswith('"') or not entry.endswith('",'):
            raise ValueError(f"unexpected template source row: {source_line!r}")
        row = ast.literal_eval(entry[:-1])
        if not isinstance(row, str) or "\n" in row or "\r" in row:
            raise ValueError("template row is not a single text line")
        rows.append(row)
    if len(rows) != 578:
        raise ValueError(f"unexpected pinned corpus size: {len(rows)}")
    TARGET.parent.mkdir(parents=True, exist_ok=True)
    TARGET.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"prepared {len(rows)} depictor templates at {TARGET}")


if __name__ == "__main__":
    main()
