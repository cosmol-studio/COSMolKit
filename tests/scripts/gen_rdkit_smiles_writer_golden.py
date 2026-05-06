#!/usr/bin/env python3
"""Generate RDKit MolToSmiles golden data across writer parameter branches."""

from __future__ import annotations

import argparse
import itertools
import json
from importlib.metadata import version
from pathlib import Path
from typing import Iterable

from rdkit import Chem, RDLogger

EXPECTED_RDKIT_VERSION = "2026.3.1"

BOOL_PARAMS = [
    "do_isomeric_smiles",
    "do_kekule",
    "canonical",
    "clean_stereo",
    "all_bonds_explicit",
    "all_hs_explicit",
    "include_dative_bonds",
    "ignore_atom_map_numbers",
]


def iter_branches() -> Iterable[dict[str, object]]:
    for values in itertools.product([False, True], repeat=len(BOOL_PARAMS)):
        base = dict(zip(BOOL_PARAMS, values, strict=True))
        for rooted_at_atom in [None, "first", "last"]:
            params = {**base, "rooted_at_atom": rooted_at_atom}
            name_bits = [
                f"iso{int(params['do_isomeric_smiles'])}",
                f"kek{int(params['do_kekule'])}",
                f"can{int(params['canonical'])}",
                f"clean{int(params['clean_stereo'])}",
                f"bond{int(params['all_bonds_explicit'])}",
                f"hs{int(params['all_hs_explicit'])}",
                f"dat{int(params['include_dative_bonds'])}",
                f"map{int(params['ignore_atom_map_numbers'])}",
                f"root_{rooted_at_atom or 'none'}",
            ]
            yield {"name": "_".join(name_bits), "params": params}


def assert_rdkit_version() -> None:
    actual = version("rdkit")
    assert actual == EXPECTED_RDKIT_VERSION, (
        f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, got {actual}"
    )


def iter_smiles(path: Path) -> Iterable[str]:
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        yield line


def rdkit_params(mol: Chem.Mol, branch_params: dict[str, object]) -> Chem.SmilesWriteParams:
    params = Chem.SmilesWriteParams()
    params.doIsomericSmiles = bool(branch_params["do_isomeric_smiles"])
    params.doKekule = bool(branch_params["do_kekule"])
    params.canonical = bool(branch_params["canonical"])
    params.cleanStereo = bool(branch_params["clean_stereo"])
    params.allBondsExplicit = bool(branch_params["all_bonds_explicit"])
    params.allHsExplicit = bool(branch_params["all_hs_explicit"])
    params.includeDativeBonds = bool(branch_params["include_dative_bonds"])
    params.ignoreAtomMapNumbers = bool(branch_params["ignore_atom_map_numbers"])
    rooted_at_atom = branch_params["rooted_at_atom"]
    if rooted_at_atom is None:
        params.rootedAtAtom = -1
    elif rooted_at_atom == "first":
        params.rootedAtAtom = 0 if mol.GetNumAtoms() else -1
    elif rooted_at_atom == "last":
        params.rootedAtAtom = mol.GetNumAtoms() - 1
    else:
        raise AssertionError(f"unknown rooted_at_atom value: {rooted_at_atom!r}")
    params.doRandom = False
    return params


def branch_result(mol: Chem.Mol, branch: dict[str, object]) -> dict[str, object]:
    params = branch["params"]
    try:
        return {
            "params": params,
            "ok": True,
            "smiles": Chem.MolToSmiles(mol, rdkit_params(mol, params)),
            "error": None,
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "params": params,
            "ok": False,
            "smiles": None,
            "error": f"{type(exc).__name__}: {exc}",
        }


def build_record(smiles: str) -> dict[str, object]:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return {
            "smiles": smiles,
            "rdkit_ok": False,
            "branches": {},
            "error": "MolFromSmiles returned None",
        }
    branches = {branch["name"]: branch_result(mol, branch) for branch in iter_branches()}
    return {
        "smiles": smiles,
        "rdkit_ok": True,
        "branches": branches,
        "error": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("tests/smiles.smi"),
        help="input SMILES file (default: tests/smiles.smi)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/golden/smiles_writer.jsonl"),
        help="output JSONL path (default: tests/golden/smiles_writer.jsonl)",
    )
    args = parser.parse_args()

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    records = [build_record(smiles) for smiles in iter_smiles(args.input)]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as f:
        for record in records:
            f.write(json.dumps(record, ensure_ascii=True))
            f.write("\n")
    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
