#!/usr/bin/env python3
"""Generate RDKit kekulized molblock golden data for bond-block parity tests."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

from rdkit import Chem
from rdkit.Chem import AllChem


def iter_smiles(path: Path) -> Iterable[str]:
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        yield line


def molblock_body(block: str) -> str:
    lines = block.splitlines()
    if len(lines) < 4:
        return ""
    return "\n".join(lines[3:])


def render_molblock_variant(
    mol: Chem.Mol, *, force_v3000: bool
) -> tuple[bool, str | None, str | None]:
    try:
        mb = Chem.MolToMolBlock(mol, forceV3000=force_v3000)
    except Exception as exc:
        return False, None, str(exc)
    body = molblock_body(mb)
    first = body.splitlines()[0] if body else ""
    expected_tag = "V3000" if force_v3000 else "V2000"
    if expected_tag not in first:
        return False, None, f"RDKit returned {first!r} instead of {expected_tag}"
    return True, body, None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("tests/smiles.smi"),
        help="input smiles corpus (default: tests/smiles.smi)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/golden/molblock_v2000_kekulized.jsonl"),
        help="output jsonl (default: tests/golden/molblock_v2000_kekulized.jsonl)",
    )
    args = parser.parse_args()

    rows = list(iter_smiles(args.input))
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w", encoding="utf-8") as f:
        for smi in rows:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                rec = {
                    "smiles": smi,
                    "parse_ok": False,
                    "parse_error": "MolFromSmiles returned None",
                    "kekulize_ok": False,
                    "kekulize_error": "parse failed",
                    "v2000_ok": False,
                    "v2000_body": None,
                    "v2000_error": "parse failed",
                    "v3000_ok": False,
                    "v3000_body": None,
                    "v3000_error": "parse failed",
                }
            else:
                parse_ok = True
                parse_error = None
                try:
                    Chem.Kekulize(mol, clearAromaticFlags=True)
                    kek_ok = True
                    kek_error = None
                except Exception as exc:
                    kek_ok = False
                    kek_error = str(exc)

                if kek_ok:
                    AllChem.Compute2DCoords(mol)
                    v2000_ok, v2000_body, v2000_error = render_molblock_variant(
                        mol, force_v3000=False
                    )
                    v3000_ok, v3000_body, v3000_error = render_molblock_variant(
                        mol, force_v3000=True
                    )
                else:
                    v2000_ok, v2000_body, v2000_error = (
                        False,
                        None,
                        "kekulize failed",
                    )
                    v3000_ok, v3000_body, v3000_error = (
                        False,
                        None,
                        "kekulize failed",
                    )

                rec = {
                    "smiles": smi,
                    "parse_ok": parse_ok,
                    "parse_error": parse_error,
                    "kekulize_ok": kek_ok,
                    "kekulize_error": kek_error,
                    "v2000_ok": v2000_ok,
                    "v2000_body": v2000_body,
                    "v2000_error": v2000_error,
                    "v3000_ok": v3000_ok,
                    "v3000_body": v3000_body,
                    "v3000_error": v3000_error,
                }

            f.write(json.dumps(rec, ensure_ascii=True))
            f.write("\n")

    print(f"Wrote {len(rows)} records to {args.output}")


if __name__ == "__main__":
    main()
