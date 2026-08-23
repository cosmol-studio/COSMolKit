#!/usr/bin/env python3
"""Generate exact pinned-RDKit Topological Torsion fingerprint expected data."""

from __future__ import annotations

import argparse
import json
import warnings
from importlib.metadata import version
from pathlib import Path
from typing import Any, Callable, Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.AtomPairs import Torsions, Utils

EXPECTED_RDKIT_VERSION = "2026.3.1"
UINT32_MASK = (1 << 32) - 1

PROFILE_PATH = Path(__file__).with_name("topological_torsion_fingerprint_profile.json")


def load_profile_matrix(name: str) -> tuple[dict[str, Any], ...]:
    profile = json.loads(PROFILE_PATH.read_text(encoding="utf-8"))
    branches = profile[name]
    if not isinstance(branches, list) or not branches:
        raise SystemExit(f"invalid Topological Torsion profile matrix {name!r}")
    return tuple(branches)


# Expected-data generation and the ChEMBL audit share these exact source
# branches, so there is no independent test-only parameter matrix to drift.
FOCUSED_PROFILES = load_profile_matrix("focused_branches")
CORPUS_PROFILES = load_profile_matrix("corpus_branches")


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


def unsigned_bit_id(bit_id: int) -> int:
    return int(bit_id) & UINT32_MASK


def exception_record(call: Callable[[], object]) -> dict[str, str | None]:
    try:
        call()
    except Exception as error:  # noqa: BLE001 - exception identity is oracle data
        return {"type": type(error).__name__, "message": str(error)}
    return {"type": None, "message": None}


def sparse_count_record(fingerprint: object) -> dict[str, object]:
    return {
        "size": int(fingerprint.GetLength()),
        "nonzero_elements": {
            str(int(bit)): int(count)
            for bit, count in sorted(fingerprint.GetNonzeroElements().items())
        },
    }


def sparse_bit_record(fingerprint: object) -> dict[str, object]:
    on_bits = sorted(unsigned_bit_id(bit) for bit in fingerprint.GetOnBits())
    return {
        "size": int(fingerprint.GetNumBits()),
        "on_bits": on_bits,
    }


def explicit_bit_record(fingerprint: object) -> dict[str, object]:
    return {
        "size": int(fingerprint.GetNumBits()),
        "on_bits": list(fingerprint.GetOnBits()),
    }


def make_additional_output() -> rdFingerprintGenerator.AdditionalOutput:
    output = rdFingerprintGenerator.AdditionalOutput()
    output.AllocateAtomToBits()
    output.AllocateBitInfoMap()
    output.AllocateBitPaths()
    output.AllocateAtomCounts()
    output.AllocateAtomsPerBit()
    return output


def additional_output_record(
    output: rdFingerprintGenerator.AdditionalOutput,
) -> dict[str, object]:
    return {
        "atom_to_bits": [[int(bit) for bit in bits] for bits in output.GetAtomToBits()],
        "bit_info_map": {
            str(int(bit)): [list(pair) for pair in pairs]
            for bit, pairs in sorted(output.GetBitInfoMap().items())
        },
        "bit_paths": {
            str(int(bit)): [list(path) for path in paths]
            for bit, paths in sorted(output.GetBitPaths().items())
        },
        "atom_counts": list(output.GetAtomCounts()),
        "atoms_per_bit": {
            str(int(bit)): [list(atoms) for atoms in paths]
            for bit, paths in sorted(output.GetAtomsPerBit().items())
        },
    }


def make_generator(profile: dict[str, Any]):
    generator = rdFingerprintGenerator.GetTopologicalTorsionGenerator(
        includeChirality=profile.get("includeChirality", False),
        torsionAtomCount=profile.get("torsionAtomCount", 4),
        countSimulation=profile.get("countSimulation", True),
        countBounds=profile.get("countBounds"),
        fpSize=profile.get("fpSize", 2048),
    )
    options = generator.GetOptions()
    options.onlyShortestPaths = profile.get("onlyShortestPaths", False)
    options.numBitsPerFeature = profile.get("numBitsPerFeature", 1)
    return generator


def call_arguments(mol: Chem.Mol, profile: dict[str, Any]) -> dict[str, object]:
    result: dict[str, object] = {}
    if profile.get("fromAtoms") == "first" and mol.GetNumAtoms():
        result["fromAtoms"] = [0]
    if profile.get("ignoreAtoms") == "first" and mol.GetNumAtoms():
        result["ignoreAtoms"] = [0]
    if profile.get("customAtomInvariants") == "index_plus_17":
        result["customAtomInvariants"] = [
            index + 17 for index in range(mol.GetNumAtoms())
        ]
    return result


def call_with_optional_output(
    generator: object,
    method: str,
    mol: Chem.Mol,
    profile: dict[str, Any],
) -> tuple[object, dict[str, object] | None]:
    kwargs = call_arguments(mol, profile)
    output = make_additional_output() if profile.get("additionalOutput") else None
    if output is not None:
        kwargs["additionalOutput"] = output
    fingerprint = getattr(generator, method)(mol, **kwargs)
    return fingerprint, None if output is None else additional_output_record(output)


def modern_profile_record(mol: Chem.Mol, profile: dict[str, Any]) -> dict[str, object]:
    generator = make_generator(profile)
    sparse_count, sparse_count_output = call_with_optional_output(
        generator, "GetSparseCountFingerprint", mol, profile
    )
    sparse_bit, sparse_bit_output = call_with_optional_output(
        generator, "GetSparseFingerprint", mol, profile
    )
    count, count_output = call_with_optional_output(
        generator, "GetCountFingerprint", mol, profile
    )
    bit, bit_output = call_with_optional_output(
        generator, "GetFingerprint", mol, profile
    )
    bulk = {
        "sparse_count": sparse_count_record(
            generator.GetSparseCountFingerprints([mol], numThreads=2)[0]
        ),
        "sparse_bit": sparse_bit_record(
            generator.GetSparseFingerprints([mol], numThreads=2)[0]
        ),
        "count": sparse_count_record(
            generator.GetCountFingerprints([mol], numThreads=2)[0]
        ),
        "bit": explicit_bit_record(generator.GetFingerprints([mol], numThreads=2)[0]),
    }
    payload = generator.ToJSON()
    restored = rdFingerprintGenerator.FingerprintGeneratorFromJSON(payload)
    record: dict[str, object] = {
        "parameters": profile,
        "info_string": generator.GetInfoString(),
        "json": json.loads(payload),
        "json_restored_count": sparse_count_record(
            restored.GetCountFingerprint(mol, **call_arguments(mol, profile))
        ),
        "sparse_count": sparse_count_record(sparse_count),
        "sparse_bit": sparse_bit_record(sparse_bit),
        "count": sparse_count_record(count),
        "bit": explicit_bit_record(bit),
        "bulk": bulk,
    }
    if profile.get("additionalOutput"):
        record["additional_output"] = {
            "sparse_count": sparse_count_output,
            "sparse_bit": sparse_bit_output,
            "count": count_output,
            "bit": bit_output,
        }
    return record


def legacy_record(mol: Chem.Mol) -> dict[str, object]:
    atom_invariants = [index + 17 for index in range(mol.GetNumAtoms())]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        unfolded = rdMolDescriptors.GetTopologicalTorsionFingerprint(mol)
        unfolded_chiral = rdMolDescriptors.GetTopologicalTorsionFingerprint(
            mol, includeChirality=True
        )
        unfolded_custom = rdMolDescriptors.GetTopologicalTorsionFingerprint(
            mol, atomInvariants=atom_invariants
        )
        hashed = rdMolDescriptors.GetHashedTopologicalTorsionFingerprint(
            mol, nBits=1000
        )
        hashed_rooted = rdMolDescriptors.GetHashedTopologicalTorsionFingerprint(
            mol, nBits=1000, fromAtoms=[0] if mol.GetNumAtoms() else []
        )
        hashed_ignored = rdMolDescriptors.GetHashedTopologicalTorsionFingerprint(
            mol, nBits=1000, ignoreAtoms=[0] if mol.GetNumAtoms() else []
        )
        bit_vectors = {
            str(bits_per_entry): explicit_bit_record(
                rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(
                    mol, nBits=256, nBitsPerEntry=bits_per_entry
                )
            )
            for bits_per_entry in (1, 2, 4, 6)
        }
    return {
        "unfolded": sparse_count_record(unfolded),
        "unfolded_chiral": sparse_count_record(unfolded_chiral),
        "unfolded_custom": sparse_count_record(unfolded_custom),
        "hashed": sparse_count_record(hashed),
        "hashed_rooted": sparse_count_record(hashed_rooted),
        "hashed_ignored": sparse_count_record(hashed_ignored),
        "bit_vectors": bit_vectors,
    }


def helper_record(mol: Chem.Mol) -> dict[str, object]:
    atom_codes = []
    for atom in mol.GetAtoms():
        code = rdMolDescriptors.GetAtomPairAtomCode(atom)
        chiral_code = rdMolDescriptors.GetAtomPairAtomCode(
            atom, includeChirality=True
        )
        atom_codes.append(
            {
                "code": int(code),
                "explanation": list(Utils.ExplainAtomCode(code)),
                "chiral_code": int(chiral_code),
                "chiral_explanation": list(
                    Utils.ExplainAtomCode(chiral_code, includeChirality=True)
                ),
            }
        )
    paths = [
        list(path) for path in Chem.FindAllPathsOfLengthN(mol, 4, useBonds=False)
    ]
    path_scores = []
    for path in paths:
        score = int(Torsions.pyScorePath(mol, path, 4))
        path_scores.append(
            {
                "path": path,
                "score": score,
                "explanation": [list(entry) for entry in Torsions.ExplainPathScore(score)],
            }
        )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        ids = [int(value) for value in Torsions.GetTopologicalTorsionFingerprintAsIds(mol)]
    return {"atom_codes": atom_codes, "paths": path_scores, "ids": ids}


def error_record(mol: Chem.Mol) -> dict[str, object]:
    def empty_count_bounds() -> object:
        generator = rdFingerprintGenerator.GetTopologicalTorsionGenerator(countBounds=[])
        return generator.GetFingerprint(mol)

    def short_custom_invariants() -> object:
        generator = rdFingerprintGenerator.GetTopologicalTorsionGenerator()
        return generator.GetCountFingerprint(mol, customAtomInvariants=[1])

    def short_path() -> object:
        return Torsions.pyScorePath(mol, [0], 4)

    return {
        "oversized_torsion": exception_record(
            lambda: rdFingerprintGenerator.GetTopologicalTorsionGenerator(
                torsionAtomCount=8
            )
        ),
        "empty_count_bounds": exception_record(empty_count_bounds),
        "short_custom_invariants": exception_record(short_custom_invariants),
        "invalid_json": exception_record(
            lambda: rdFingerprintGenerator.FingerprintGeneratorFromJSON("not json")
        ),
        "short_path": exception_record(short_path),
    }


def build_records(
    smiles_values: list[str], profile_matrix: tuple[dict[str, Any], ...]
) -> list[dict[str, object]]:
    records = []
    for row, smiles in enumerate(smiles_values):
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol is None:
            records.append(
                {
                    "row": row,
                    "smiles": smiles,
                    "rdkit_ok": False,
                    "profiles": {},
                    "legacy": None,
                    "helpers": None,
                    "errors": None,
                    "error": "MolFromSmiles returned None",
                }
            )
            continue
        try:
            profile_records = {
                profile["name"]: modern_profile_record(mol, profile)
                for profile in profile_matrix
            }
            records.append(
                {
                    "row": row,
                    "smiles": smiles,
                    "rdkit_ok": True,
                    "profiles": profile_records,
                    "legacy": legacy_record(mol),
                    "helpers": helper_record(mol),
                    "errors": error_record(mol),
                    "error": None,
                }
            )
        except Exception as error:  # noqa: BLE001 - failed rows remain explicit
            records.append(
                {
                    "row": row,
                    "smiles": smiles,
                    "rdkit_ok": False,
                    "profiles": {},
                    "legacy": None,
                    "helpers": None,
                    "errors": None,
                    "error": f"{type(error).__name__}: {error}",
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
    profiles = (
        CORPUS_PROFILES
        if args.input.stem in {"smiles_5000", "5000"}
        else FOCUSED_PROFILES
    )
    records = build_records(list(iter_smiles(args.input)), profiles)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=True, sort_keys=True))
            handle.write("\n")
    print(
        f"Wrote {len(records)} deterministic Topological Torsion records "
        f"across {len(profiles)} source profiles "
        f"to {args.output}"
    )


if __name__ == "__main__":
    main()
