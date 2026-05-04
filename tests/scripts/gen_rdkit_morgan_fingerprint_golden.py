#!/usr/bin/env python3
"""Generate RDKit Morgan fingerprint golden data across parameter branches."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
from typing import Iterable

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

EXPECTED_RDKIT_VERSION = "2026.3.1"

BRANCHES = [
    {
        "name": "r2_n2048_chiral_false_bonds_true",
        "radius": 2,
        "nBits": 2048,
        "includeChirality": False,
        "useBondTypes": True,
    },
    {
        "name": "r2_n2048_chiral_true_bonds_true",
        "radius": 2,
        "nBits": 2048,
        "includeChirality": True,
        "useBondTypes": True,
    },
    {
        "name": "r3_n2048_chiral_false_bonds_true",
        "radius": 3,
        "nBits": 2048,
        "includeChirality": False,
        "useBondTypes": True,
    },
    {
        "name": "r2_n1024_chiral_false_bonds_false",
        "radius": 2,
        "nBits": 1024,
        "includeChirality": False,
        "useBondTypes": False,
    },
    {
        "name": "r2_n2048_count_sim_bounds_1_2_4_8",
        "radius": 2,
        "nBits": 2048,
        "includeChirality": False,
        "useBondTypes": True,
        "countSimulation": True,
        "countBounds": [1, 2, 4, 8],
    },
    {
        "name": "r2_n2048_redundant_true",
        "radius": 2,
        "nBits": 2048,
        "includeChirality": False,
        "useBondTypes": True,
        "includeRedundantEnvironments": True,
    },
    {
        "name": "r2_n2048_ring_membership_false",
        "radius": 2,
        "nBits": 2048,
        "includeChirality": False,
        "useBondTypes": True,
        "includeRingMembership": False,
    },
    {
        "name": "r2_n2048_from_atom_0",
        "radius": 2,
        "nBits": 2048,
        "includeChirality": False,
        "useBondTypes": True,
        "fromAtoms": "first",
    },
    {
        "name": "r2_n2048_ignore_atom_0",
        "radius": 2,
        "nBits": 2048,
        "includeChirality": False,
        "useBondTypes": True,
        "ignoreAtoms": "first",
    },
    {
        "name": "r2_n2048_custom_atom_invariants",
        "radius": 2,
        "nBits": 2048,
        "includeChirality": False,
        "useBondTypes": True,
        "customAtomInvariants": "index_plus_one",
    },
    {
        "name": "r2_n2048_custom_bond_invariants",
        "radius": 2,
        "nBits": 2048,
        "includeChirality": False,
        "useBondTypes": True,
        "customBondInvariants": "index_plus_seven",
    },
    {
        "name": "r2_n2048_only_nonzero_custom_atom",
        "radius": 2,
        "nBits": 2048,
        "includeChirality": False,
        "useBondTypes": True,
        "onlyNonzeroInvariants": True,
        "customAtomInvariants": "even_atoms_zero_odd_index_plus_one",
    },
    {
        "name": "r2_n2048_atom_inv_gen_ring_false",
        "radius": 2,
        "nBits": 2048,
        "includeChirality": False,
        "useBondTypes": True,
        "atomInvariantsGenerator": "morgan_ring_false",
    },
    {
        "name": "r2_n2048_feature_atom_inv_gen",
        "radius": 2,
        "nBits": 2048,
        "includeChirality": False,
        "useBondTypes": True,
        "atomInvariantsGenerator": "morgan_feature_default",
    },
    {
        "name": "r2_n2048_bond_inv_gen_no_bond_types",
        "radius": 2,
        "nBits": 2048,
        "includeChirality": False,
        "useBondTypes": True,
        "bondInvariantsGenerator": "morgan_no_bond_types",
    },
    {
        "name": "r2_n2048_num_bits_per_feature_2",
        "radius": 2,
        "nBits": 2048,
        "includeChirality": False,
        "useBondTypes": True,
        "numBitsPerFeature": 2,
    },
    {
        "name": "r2_n2048_additional_output",
        "radius": 2,
        "nBits": 2048,
        "includeChirality": False,
        "useBondTypes": True,
        "additionalOutput": True,
    },
]


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


def make_generator(branch: dict[str, object]):
    atom_inv_gen = None
    if branch.get("atomInvariantsGenerator") == "morgan_ring_false":
        atom_inv_gen = rdFingerprintGenerator.GetMorganAtomInvGen(False)
    elif branch.get("atomInvariantsGenerator") == "morgan_feature_default":
        atom_inv_gen = rdFingerprintGenerator.GetMorganFeatureAtomInvGen()

    bond_inv_gen = None
    if branch.get("bondInvariantsGenerator") == "morgan_no_bond_types":
        bond_inv_gen = rdFingerprintGenerator.GetMorganBondInvGen(False, False)

    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=branch["radius"],
        fpSize=branch["nBits"],
        includeChirality=branch["includeChirality"],
        useBondTypes=branch["useBondTypes"],
        countSimulation=branch.get("countSimulation", False),
        countBounds=branch.get("countBounds"),
        onlyNonzeroInvariants=branch.get("onlyNonzeroInvariants", False),
        includeRingMembership=branch.get("includeRingMembership", True),
        includeRedundantEnvironments=branch.get("includeRedundantEnvironments", False),
        atomInvariantsGenerator=atom_inv_gen,
        bondInvariantsGenerator=bond_inv_gen,
    )
    if "numBitsPerFeature" in branch:
        generator.GetOptions().numBitsPerFeature = branch["numBitsPerFeature"]
    return generator


def fingerprint_kwargs(mol: Chem.Mol, branch: dict[str, object]) -> dict[str, object]:
    kwargs: dict[str, object] = {}
    if branch.get("fromAtoms") == "first" and mol.GetNumAtoms() > 0:
        kwargs["fromAtoms"] = [0]
    if branch.get("ignoreAtoms") == "first" and mol.GetNumAtoms() > 0:
        kwargs["ignoreAtoms"] = [0]
    custom_atom = branch.get("customAtomInvariants")
    if custom_atom == "index_plus_one":
        kwargs["customAtomInvariants"] = [idx + 1 for idx in range(mol.GetNumAtoms())]
    elif custom_atom == "even_atoms_zero_odd_index_plus_one":
        kwargs["customAtomInvariants"] = [
            0 if idx % 2 == 0 else idx + 1 for idx in range(mol.GetNumAtoms())
        ]
    if branch.get("customBondInvariants") == "index_plus_seven":
        kwargs["customBondInvariants"] = [idx + 7 for idx in range(mol.GetNumBonds())]
    return kwargs


def build_records(smiles_values: list[str]) -> list[dict[str, object]]:
    previous_fps: dict[str, object] = {}
    records: list[dict[str, object]] = []
    generators = {branch["name"]: make_generator(branch) for branch in BRANCHES}

    for smiles in smiles_values:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol is None:
            records.append(
                {
                    "smiles": smiles,
                    "rdkit_ok": False,
                    "branches": {},
                    "error": "MolFromSmiles returned None",
                }
            )
            previous_fps.clear()
            continue

        branches = {}
        for branch in BRANCHES:
            name = branch["name"]
            try:
                additional_output = None
                if branch.get("additionalOutput"):
                    additional_output = rdFingerprintGenerator.AdditionalOutput()
                    additional_output.AllocateAtomCounts()
                    additional_output.AllocateAtomToBits()
                    additional_output.AllocateBitInfoMap()
                    additional_output.AllocateAtomsPerBit()
                kwargs = fingerprint_kwargs(mol, branch)
                if additional_output is not None:
                    kwargs["additionalOutput"] = additional_output
                fp = generators[name].GetFingerprint(mol, **kwargs)
                previous = previous_fps.get(name)
                tanimoto = (
                    None
                    if previous is None
                    else DataStructs.TanimotoSimilarity(fp, previous)
                )
                branch_record = {
                    "ok": True,
                    "on_bits": list(fp.GetOnBits()),
                    "num_on_bits": fp.GetNumOnBits(),
                    "tanimoto_to_previous": tanimoto,
                    "error": None,
                }
                if additional_output is not None:
                    branch_record["additional_output"] = {
                        "atom_counts": list(additional_output.GetAtomCounts()),
                        "atom_to_bits": [
                            list(bits) for bits in additional_output.GetAtomToBits()
                        ],
                        "bit_info_map": {
                            str(bit): [list(pair) for pair in pairs]
                            for bit, pairs in sorted(additional_output.GetBitInfoMap().items())
                        },
                        "atoms_per_bit": {
                            str(bit): [list(atoms) for atoms in atoms_per_bit]
                            for bit, atoms_per_bit in sorted(additional_output.GetAtomsPerBit().items())
                        },
                    }
                branches[name] = branch_record
                previous_fps[name] = fp
            except Exception as exc:  # noqa: BLE001
                branches[name] = {
                    "ok": False,
                    "on_bits": None,
                    "num_on_bits": None,
                    "tanimoto_to_previous": None,
                    "error": f"{type(exc).__name__}: {exc}",
                }
                previous_fps.pop(name, None)

        records.append(
            {
                "smiles": smiles,
                "rdkit_ok": True,
                "branches": branches,
                "error": None,
            }
        )

    return records


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
        default=Path("tests/golden/morgan_fingerprint.jsonl"),
        help="output JSONL path (default: tests/golden/morgan_fingerprint.jsonl)",
    )
    args = parser.parse_args()

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    records = build_records(list(iter_smiles(args.input)))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as f:
        for record in records:
            f.write(json.dumps(record, ensure_ascii=True))
            f.write("\n")
    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
