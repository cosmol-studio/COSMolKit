#!/usr/bin/env python3
"""Generate exact RDKit AtomPair fingerprint expected data for a corpus."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
from typing import Any, Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator

EXPECTED_RDKIT_VERSION = "2026.3.1"
PROFILE_PATH = Path(__file__).with_name("atom_pair_fingerprint_profile.json")
UINT32_MASK = (1 << 32) - 1


def assert_rdkit_version() -> None:
    actual = version("rdkit")
    if actual != EXPECTED_RDKIT_VERSION:
        raise SystemExit(
            f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, got {actual}"
        )


def iter_smiles(path: Path) -> Iterable[str]:
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line and not line.startswith("#"):
            yield line


def load_profile() -> dict[str, Any]:
    profile = json.loads(PROFILE_PATH.read_text(encoding="utf-8"))
    if profile.get("schema_version") != 1:
        raise SystemExit("AtomPair fingerprint profile schema_version must be 1")
    if profile.get("rdkit_version") != EXPECTED_RDKIT_VERSION:
        raise SystemExit("AtomPair fingerprint profile has the wrong RDKit version")
    source_paths = profile.get("source_paths")
    if not isinstance(source_paths, list) or not source_paths:
        raise SystemExit("AtomPair fingerprint profile must identify its RDKit sources")
    comparison_fields = profile.get("comparison_fields")
    if not isinstance(comparison_fields, dict) or not comparison_fields:
        raise SystemExit("AtomPair fingerprint profile must define comparison fields")
    branches = profile.get("branches")
    if not isinstance(branches, list) or not branches:
        raise SystemExit("AtomPair fingerprint profile must contain branches")
    names = [branch.get("name") for branch in branches]
    if any(not isinstance(name, str) or not name for name in names):
        raise SystemExit("AtomPair fingerprint branch names must be non-empty strings")
    if len(set(names)) != len(names):
        raise SystemExit("AtomPair fingerprint branch names must be unique")
    return profile


def make_generator(branch: dict[str, Any]):
    generator = rdFingerprintGenerator.GetAtomPairGenerator(
        minDistance=branch["minDistance"],
        maxDistance=branch["maxDistance"],
        includeChirality=branch["includeChirality"],
        use2D=branch["use2D"],
        countSimulation=branch["countSimulation"],
        countBounds=branch["countBounds"],
        fpSize=branch["fpSize"],
    )
    generator.GetOptions().numBitsPerFeature = branch["numBitsPerFeature"]
    return generator


def resolved_arguments(mol: Chem.Mol, branch: dict[str, Any]) -> dict[str, list[int]]:
    arguments: dict[str, list[int]] = {}
    if branch.get("fromAtoms") == "first":
        arguments["fromAtoms"] = [0] if mol.GetNumAtoms() else []
    if branch.get("ignoreAtoms") == "first":
        arguments["ignoreAtoms"] = [0] if mol.GetNumAtoms() else []
    if branch.get("customAtomInvariants") == "index_plus_11":
        arguments["customAtomInvariants"] = [
            index + 11 for index in range(mol.GetNumAtoms())
        ]
    return arguments


def make_additional_output() -> rdFingerprintGenerator.AdditionalOutput:
    output = rdFingerprintGenerator.AdditionalOutput()
    output.AllocateAtomCounts()
    output.AllocateAtomToBits()
    output.AllocateBitInfoMap()
    output.AllocateAtomsPerBit()
    return output


def unsigned_bit_id(bit_id: int) -> int:
    return int(bit_id) & UINT32_MASK


def additional_output_record(
    output: rdFingerprintGenerator.AdditionalOutput | None,
) -> dict[str, Any] | None:
    if output is None:
        return None
    return {
        "atom_counts": list(output.GetAtomCounts()),
        "atom_to_bits": [
            [unsigned_bit_id(bit) for bit in bits] for bits in output.GetAtomToBits()
        ],
        "bit_info_map": {
            str(unsigned_bit_id(bit)): [list(pair) for pair in pairs]
            for bit, pairs in sorted(output.GetBitInfoMap().items())
        },
        "atoms_per_bit": {
            str(unsigned_bit_id(bit)): [list(atoms) for atoms in entries]
            for bit, entries in sorted(output.GetAtomsPerBit().items())
        },
    }


def call_output(
    generator: Any,
    method_name: str,
    mol: Chem.Mol,
    arguments: dict[str, list[int]],
    capture_additional_output: bool,
) -> tuple[Any, dict[str, Any] | None]:
    output = make_additional_output() if capture_additional_output else None
    kwargs: dict[str, Any] = dict(arguments)
    if output is not None:
        kwargs["additionalOutput"] = output
    fingerprint = getattr(generator, method_name)(mol, **kwargs)
    return fingerprint, additional_output_record(output)


def output_record(
    generator: Any,
    method_name: str,
    output_kind: str,
    mol: Chem.Mol,
    arguments: dict[str, list[int]],
    capture_additional_output: bool,
) -> dict[str, Any]:
    try:
        fingerprint, additional_output = call_output(
            generator,
            method_name,
            mol,
            arguments,
            capture_additional_output,
        )
        if output_kind == "count":
            result: dict[str, Any] = {
                "ok": True,
                "length": int(fingerprint.GetLength()),
                "nonzero_elements": {
                    str(int(bit)): int(count)
                    for bit, count in sorted(
                        fingerprint.GetNonzeroElements().items()
                    )
                },
                "error": None,
            }
        else:
            on_bits = sorted(unsigned_bit_id(bit) for bit in fingerprint.GetOnBits())
            result = {
                "ok": True,
                "length": int(fingerprint.GetNumBits()),
                "on_bits": on_bits,
                "error": None,
            }
        if additional_output is not None:
            result["additional_output"] = additional_output
        return result
    except Exception as exc:  # noqa: BLE001
        return {
            "ok": False,
            "length": None,
            "nonzero_elements" if output_kind == "count" else "on_bits": None,
            "error": f"{type(exc).__name__}: {exc}",
        }


def build_records(
    smiles_values: list[str], profile: dict[str, Any]
) -> list[dict[str, Any]]:
    branches = profile["branches"]
    generators = {branch["name"]: make_generator(branch) for branch in branches}
    records: list[dict[str, Any]] = []
    for row, smiles in enumerate(smiles_values):
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol is None:
            records.append(
                {
                    "row": row,
                    "smiles": smiles,
                    "rdkit_ok": False,
                    "branches": {},
                    "error": "MolFromSmiles returned None",
                }
            )
            continue

        branch_records: dict[str, Any] = {}
        for branch in branches:
            name = branch["name"]
            arguments = resolved_arguments(mol, branch)
            capture = bool(branch.get("additionalOutput", False))
            outputs = {
                "sparse_count": output_record(
                    generators[name],
                    "GetSparseCountFingerprint",
                    "count",
                    mol,
                    arguments,
                    capture,
                ),
                "count": output_record(
                    generators[name],
                    "GetCountFingerprint",
                    "count",
                    mol,
                    arguments,
                    capture,
                ),
                "sparse_bit": output_record(
                    generators[name],
                    "GetSparseFingerprint",
                    "bit",
                    mol,
                    arguments,
                    capture,
                ),
                "explicit_bit": output_record(
                    generators[name],
                    "GetFingerprint",
                    "bit",
                    mol,
                    arguments,
                    capture,
                ),
            }
            ok = all(output["ok"] for output in outputs.values())
            branch_records[name] = {
                "parameters": branch,
                "resolved_arguments": arguments,
                **outputs,
                "ok": ok,
                "error": None
                if ok
                else "one or more AtomPair fingerprint calls failed",
            }

        records.append(
            {
                "row": row,
                "smiles": smiles,
                "rdkit_ok": True,
                "branches": branch_records,
                "error": None,
            }
        )
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    profile = load_profile()
    records = build_records(list(iter_smiles(args.input)), profile)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=True, sort_keys=True))
            handle.write("\n")
    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
