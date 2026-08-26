#!/usr/bin/env python3
# pyright: reportAttributeAccessIssue=false
"""Exact ChEMBL parity audit for repository-supported fingerprint families."""

from __future__ import annotations

import argparse
import json
import os
import time
from collections import Counter
from pathlib import Path
from typing import Any

import cosmolkit
from rdkit import Chem, DataStructs, RDLogger, rdBase
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import rdFingerprintGenerator


EXPECTED_RDKIT_VERSION = "2026.03.1"
EXPECTED_RDKIT_SOURCE_REVISION = "351f8f378f8ad6bbd517980c38896e66bf907af8"
EXPECTED_LAYERED_VERSION = "0.7.0"
UINT32_MASK = (1 << 32) - 1
RDLogger.DisableLog("rdApp.*")


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def atomic_json(path: Path, value: Any) -> None:
    temporary = path.with_suffix(path.suffix + f".tmp-{os.getpid()}")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def error_value(error: BaseException) -> dict[str, str]:
    return {"type": type(error).__name__, "message": str(error)}


class Audit:
    def __init__(self, output: Path, max_examples: int) -> None:
        self.output = output
        self.findings_path = output.with_suffix(".findings.jsonl")
        self.max_examples = max_examples
        self.counts: Counter[str] = Counter()
        self.example_counts: Counter[str] = Counter()
        self.started = time.monotonic()
        output.parent.mkdir(parents=True, exist_ok=True)
        self.findings = self.findings_path.open("w", encoding="utf-8")

    def match(self, feature: str) -> None:
        self.counts[f"match.{feature}"] += 1

    def count(self, feature: str, amount: int = 1) -> None:
        self.counts[feature] += amount

    def mismatch(
        self,
        feature: str,
        record: dict[str, Any],
        rdkit_value: Any,
        cosmolkit_value: Any,
        detail: str | None = None,
    ) -> None:
        key = f"mismatch.{feature}"
        self.counts[key] += 1
        if self.example_counts[key] >= self.max_examples:
            return
        self.example_counts[key] += 1
        self.findings.write(
            json.dumps(
                {
                    "feature": feature,
                    "row": record["row"],
                    "chembl_id": record["chembl_id"],
                    "smiles": record["smiles"],
                    "rdkit": rdkit_value,
                    "cosmolkit": cosmolkit_value,
                    "detail": detail,
                },
                ensure_ascii=True,
                default=str,
                separators=(",", ":"),
            )
            + "\n"
        )
        self.findings.flush()

    def compare(
        self,
        feature: str,
        record: dict[str, Any],
        rdkit_value: Any,
        cosmolkit_value: Any,
    ) -> None:
        if rdkit_value == cosmolkit_value:
            self.match(feature)
        else:
            self.mismatch(feature, record, rdkit_value, cosmolkit_value)

    def finish(self, processed: int, profiles: dict[str, int]) -> None:
        self.findings.close()
        atomic_json(
            self.output,
            {
                "processed": processed,
                "elapsed_seconds": time.monotonic() - self.started,
                "profiles": profiles,
                "counts": dict(sorted(self.counts.items())),
            },
        )


def topological_kwargs(
    branch: dict[str, Any], atom_count: int
) -> tuple[dict[str, Any], dict[str, Any]]:
    rdkit_kwargs = {
        "minPath": branch["minPath"],
        "maxPath": branch["maxPath"],
        "fpSize": branch["fpSize"],
        "nBitsPerHash": branch["nBitsPerHash"],
        "useHs": branch["useHs"],
        "tgtDensity": branch["tgtDensity"],
        "minSize": branch["minSize"],
        "branchedPaths": branch["branchedPaths"],
        "useBondOrder": branch["useBondOrder"],
    }
    cosmolkit_kwargs = {
        "min_path": branch["minPath"],
        "max_path": branch["maxPath"],
        "fp_size": branch["fpSize"],
        "num_bits_per_feature": branch["nBitsPerHash"],
        "use_hs": branch["useHs"],
        "target_density": branch["tgtDensity"],
        "min_size": branch["minSize"],
        "branched_paths": branch["branchedPaths"],
        "use_bond_order": branch["useBondOrder"],
    }
    if branch.get("fromAtoms") == "first" and atom_count:
        rdkit_kwargs["fromAtoms"] = [0]
        cosmolkit_kwargs["from_atoms"] = [0]
    if branch.get("atomInvariants") == "index_plus_one":
        invariants = list(range(1, atom_count + 1))
        rdkit_kwargs["atomInvariants"] = invariants
        cosmolkit_kwargs["atom_invariants"] = invariants
    return rdkit_kwargs, cosmolkit_kwargs


def compare_topological(
    audit: Audit,
    record: dict[str, Any],
    rd_mol: Chem.Mol,
    ck_mol: Any,
    branches: list[dict[str, Any]],
) -> None:
    for branch in branches:
        name = str(branch["name"])
        rdkit_kwargs, cosmolkit_kwargs = topological_kwargs(
            branch, rd_mol.GetNumAtoms()
        )
        try:
            rd_fp = Chem.RDKFingerprint(rd_mol, **rdkit_kwargs)
            rd_value: Any = {
                "n_bits": rd_fp.GetNumBits(),
                "on_bits": list(rd_fp.GetOnBits()),
            }
        except Exception as error:  # noqa: BLE001
            rd_value = {"error": error_value(error)}
        try:
            ck_fp = ck_mol.topological_fingerprint(**cosmolkit_kwargs)
            ck_value: Any = {
                "n_bits": ck_fp.n_bits(),
                "on_bits": ck_fp.on_bits(),
            }
        except Exception as error:  # noqa: BLE001
            ck_value = {"error": error_value(error)}
        audit.compare(f"topological.{name}", record, rd_value, ck_value)


def normalize_bit_info(value: dict[Any, Any]) -> dict[int, list[list[int]]]:
    return {
        int(bit): [list(path) for path in paths]
        for bit, paths in sorted(value.items(), key=lambda item: int(item[0]))
    }


def compare_topological_provenance(
    audit: Audit,
    record: dict[str, Any],
    rd_mol: Chem.Mol,
    ck_mol: Any,
    branches: list[dict[str, Any]],
) -> None:
    by_name = {str(branch["name"]): branch for branch in branches}
    for name in ("default", "density_folded"):
        branch = by_name[name]
        rdkit_kwargs, cosmolkit_kwargs = topological_kwargs(
            branch, rd_mol.GetNumAtoms()
        )
        atom_bits: list[list[int]] = []
        bit_info: dict[int, list[list[int]]] = {}
        try:
            rd_fp = Chem.RDKFingerprint(
                rd_mol,
                atomBits=atom_bits,
                bitInfo=bit_info,
                **rdkit_kwargs,
            )
            rd_value: Any = {
                "n_bits": rd_fp.GetNumBits(),
                "on_bits": list(rd_fp.GetOnBits()),
                "atom_bits": [list(bits) for bits in atom_bits],
                "bit_info": normalize_bit_info(bit_info),
            }
        except Exception as error:  # noqa: BLE001
            rd_value = {"error": error_value(error)}
        try:
            ck_result = ck_mol.topological_fingerprint_with_output(
                atom_bits=True,
                bit_info=True,
                **cosmolkit_kwargs,
            )
            ck_fp = ck_result.fingerprint()
            ck_value: Any = {
                "n_bits": ck_fp.n_bits(),
                "on_bits": ck_fp.on_bits(),
                "atom_bits": ck_result.atom_bits(),
                "bit_info": normalize_bit_info(ck_result.bit_info()),
            }
        except Exception as error:  # noqa: BLE001
            ck_value = {"error": error_value(error)}
        audit.compare(f"topological_provenance.{name}", record, rd_value, ck_value)


def compare_avalon(
    audit: Audit,
    record: dict[str, Any],
    rd_mol: Chem.Mol,
    ck_mol: Any,
    branches: list[dict[str, Any]],
) -> None:
    for branch in branches:
        name = str(branch["name"])
        n_bits = int(branch["n_bits"])
        is_query = bool(branch["is_query"])
        bit_flags = int(str(branch["bit_flags"]), 16)
        try:
            rd_fp = pyAvalonTools.GetAvalonFP(
                rd_mol, nBits=n_bits, isQuery=is_query, bitFlags=bit_flags
            )
            rd_value: Any = {
                "n_bits": rd_fp.GetNumBits(),
                "on_bits": list(rd_fp.GetOnBits()),
            }
        except Exception as error:  # noqa: BLE001
            rd_value = {"error": error_value(error)}
        try:
            ck_fp = ck_mol.avalon_fingerprint(
                n_bits=n_bits, is_query=is_query, bit_flags=bit_flags
            )
            ck_value: Any = {
                "n_bits": ck_fp.n_bits(),
                "on_bits": ck_fp.on_bits(),
            }
        except Exception as error:  # noqa: BLE001
            ck_value = {"error": error_value(error)}
        audit.compare(f"avalon.{name}", record, rd_value, ck_value)


def layered_branches(profile: dict[str, Any]) -> list[dict[str, Any]]:
    if profile.get("schema_version") != 1:
        raise ValueError("Layered profile schema_version must be 1")
    if profile.get("rdkit_version") != "2026.3.1":
        raise ValueError("Layered profile RDKit version does not match the pin")
    if profile.get("rdkit_source_revision") != EXPECTED_RDKIT_SOURCE_REVISION:
        raise ValueError("Layered profile source revision does not match the pin")
    if profile.get("algorithm_version") != EXPECTED_LAYERED_VERSION:
        raise ValueError("Layered profile algorithm version does not match the source")
    if profile.get("comparison_fields") != ["num_bits", "on_bits", "atom_counts"]:
        raise ValueError("Layered profile comparison_fields must cover every exact output")
    branches = profile.get("branches")
    if not isinstance(branches, list) or not branches:
        raise ValueError("Layered profile branches must be a nonempty list")
    names: set[str] = set()
    for branch in branches:
        if not isinstance(branch, dict):
            raise ValueError("every Layered branch must be an object")
        name = branch.get("name")
        if not isinstance(name, str) or not name or name in names:
            raise ValueError(f"Layered branch name is empty or duplicated: {name!r}")
        names.add(name)
        for key in ("layerFlags", "minPath", "maxPath", "fpSize", "branchedPaths"):
            if key not in branch:
                raise ValueError(f"Layered branch {name!r} is missing {key}")
        try:
            int(str(branch["layerFlags"]), 0)
        except ValueError as error:
            raise ValueError(f"Layered branch {name!r} has invalid layerFlags") from error
        if not isinstance(branch["minPath"], int) or not isinstance(
            branch["maxPath"], int
        ):
            raise ValueError(f"Layered branch {name!r} has invalid path bounds")
        if not isinstance(branch["fpSize"], int) or branch["fpSize"] < 1:
            raise ValueError(f"Layered branch {name!r} has invalid fpSize")
        if not isinstance(branch["branchedPaths"], bool):
            raise ValueError(f"Layered branch {name!r} has invalid branchedPaths")
        if branch.get("fromAtoms") not in (None, "first", "terminal_pair"):
            raise ValueError(f"Layered branch {name!r} has invalid fromAtoms selector")
        if branch.get("atomCounts") not in (None, "zeros", "index_plus_10"):
            raise ValueError(f"Layered branch {name!r} has invalid atomCounts selector")
        if branch.get("setOnlyBits") not in (None, "even", "mod_three"):
            raise ValueError(f"Layered branch {name!r} has invalid setOnlyBits selector")
    return branches


def layered_from_atoms(atom_count: int, selector: str | None) -> list[int] | None:
    if selector is None:
        return None
    if selector == "first":
        return [0] if atom_count else []
    if atom_count == 0:
        return []
    if atom_count == 1:
        return [0]
    return [0, atom_count - 1]


def layered_atom_counts(atom_count: int, selector: str | None) -> list[int] | None:
    if selector is None:
        return None
    if selector == "zeros":
        return [0] * atom_count
    return [index + 10 for index in range(atom_count)]


def layered_masks(
    fp_size: int, selector: str | None
) -> tuple[Any | None, Any | None]:
    if selector is None:
        return None, None
    step = 2 if selector == "even" else 3
    on_bits = list(range(0, fp_size, step))
    rdkit_mask = DataStructs.ExplicitBitVect(fp_size)
    for bit in on_bits:
        rdkit_mask.SetBit(bit)
    return rdkit_mask, cosmolkit.Fingerprint.from_on_bits(fp_size, on_bits)


def compare_layered(
    audit: Audit,
    record: dict[str, Any],
    rd_mol: Chem.Mol,
    ck_mol: Any,
    branches: list[dict[str, Any]],
) -> None:
    atom_count = rd_mol.GetNumAtoms()
    for branch in branches:
        name = str(branch["name"])
        fp_size = int(branch["fpSize"])
        from_atoms = layered_from_atoms(atom_count, branch.get("fromAtoms"))
        rdkit_counts = layered_atom_counts(atom_count, branch.get("atomCounts"))
        cosmolkit_counts = (
            None if rdkit_counts is None else list(rdkit_counts)
        )
        rdkit_mask, cosmolkit_mask = layered_masks(
            fp_size, branch.get("setOnlyBits")
        )
        rdkit_kwargs: dict[str, Any] = {
            "layerFlags": int(str(branch["layerFlags"]), 0),
            "minPath": int(branch["minPath"]),
            "maxPath": int(branch["maxPath"]),
            "fpSize": fp_size,
            "branchedPaths": bool(branch["branchedPaths"]),
        }
        cosmolkit_kwargs: dict[str, Any] = {
            "layers": int(str(branch["layerFlags"]), 0),
            "min_path": int(branch["minPath"]),
            "max_path": int(branch["maxPath"]),
            "fp_size": fp_size,
            "branched_paths": bool(branch["branchedPaths"]),
        }
        if from_atoms is not None:
            rdkit_kwargs["fromAtoms"] = from_atoms
            cosmolkit_kwargs["from_atoms"] = from_atoms
        if rdkit_counts is not None:
            rdkit_kwargs["atomCounts"] = rdkit_counts
            cosmolkit_kwargs["atom_counts"] = cosmolkit_counts
        if rdkit_mask is not None:
            rdkit_kwargs["setOnlyBits"] = rdkit_mask
            cosmolkit_kwargs["set_only_bits"] = cosmolkit_mask
        try:
            rdkit_fp = Chem.LayeredFingerprint(rd_mol, **rdkit_kwargs)
            rdkit_value: Any = {
                "n_bits": rdkit_fp.GetNumBits(),
                "on_bits": list(rdkit_fp.GetOnBits()),
                "atom_counts": rdkit_counts,
            }
        except Exception as error:  # noqa: BLE001
            rdkit_value = {"error": error_value(error)}
        try:
            cosmolkit_result = ck_mol.fingerprint_layered_with_output(
                **cosmolkit_kwargs
            )
            cosmolkit_fp = cosmolkit_result.fingerprint()
            cosmolkit_value: Any = {
                "n_bits": cosmolkit_fp.n_bits(),
                "on_bits": cosmolkit_fp.on_bits(),
                "atom_counts": cosmolkit_result.atom_counts(),
            }
        except Exception as error:  # noqa: BLE001
            cosmolkit_value = {"error": error_value(error)}
        audit.compare(f"layered.{name}", record, rdkit_value, cosmolkit_value)


def atom_pair_kwargs(
    branch: dict[str, Any], atom_count: int
) -> tuple[dict[str, Any], dict[str, Any]]:
    rdkit_kwargs: dict[str, Any] = {}
    cosmolkit_kwargs: dict[str, Any] = {
        "n_bits": branch["fpSize"],
        "min_distance": branch["minDistance"],
        "max_distance": branch["maxDistance"],
        "use_2d": branch["use2D"],
        "include_chirality": branch["includeChirality"],
        "count_simulation": branch["countSimulation"],
        "count_bounds": branch["countBounds"],
        "num_bits_per_feature": branch["numBitsPerFeature"],
    }
    if branch.get("fromAtoms") == "first":
        atoms = [0] if atom_count else []
        rdkit_kwargs["fromAtoms"] = atoms
        cosmolkit_kwargs["from_atoms"] = atoms
    if branch.get("ignoreAtoms") == "first":
        atoms = [0] if atom_count else []
        rdkit_kwargs["ignoreAtoms"] = atoms
        cosmolkit_kwargs["ignore_atoms"] = atoms
    if branch.get("customAtomInvariants") == "index_plus_11":
        invariants = [index + 11 for index in range(atom_count)]
        rdkit_kwargs["customAtomInvariants"] = invariants
        cosmolkit_kwargs["custom_atom_invariants"] = invariants
    return rdkit_kwargs, cosmolkit_kwargs


def make_atom_pair_generators(
    branches: list[dict[str, Any]],
) -> dict[str, Any]:
    generators: dict[str, Any] = {}
    for branch in branches:
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
        generators[str(branch["name"])] = generator
    return generators


def atom_pair_additional_output() -> rdFingerprintGenerator.AdditionalOutput:
    output = rdFingerprintGenerator.AdditionalOutput()
    output.AllocateAtomCounts()
    output.AllocateAtomToBits()
    output.AllocateBitInfoMap()
    output.AllocateAtomsPerBit()
    return output


def normalize_atom_pair_additional_output(output: Any) -> dict[str, Any]:
    return {
        "atom_counts": list(output.GetAtomCounts()),
        "atom_to_bits": [list(bits) for bits in output.GetAtomToBits()],
        "bit_info_map": {
            int(bit): [list(pair) for pair in pairs]
            for bit, pairs in sorted(output.GetBitInfoMap().items())
        },
        "atoms_per_bit": {
            int(bit): [list(atoms) for atoms in entries]
            for bit, entries in sorted(output.GetAtomsPerBit().items())
        },
    }


def normalize_cosmolkit_atom_pair_additional_output(output: Any) -> dict[str, Any]:
    return {
        "atom_counts": output.atom_counts(),
        "atom_to_bits": output.atom_to_bits(),
        "bit_info_map": {
            int(bit): [list(pair) for pair in pairs]
            for bit, pairs in sorted(output.bit_info_map().items())
        },
        "atoms_per_bit": {
            int(bit): [list(atoms) for atoms in entries]
            for bit, entries in sorted(output.atoms_per_bit().items())
        },
    }


def capture_count_fingerprint(call: Any) -> dict[str, Any]:
    try:
        fingerprint = call()
        return {
            "status": "ok",
            "length": fingerprint.GetLength()
            if hasattr(fingerprint, "GetLength")
            else fingerprint.size(),
            "nonzero_elements": {
                int(bit): int(count)
                for bit, count in sorted(
                    (
                        fingerprint.GetNonzeroElements()
                        if hasattr(fingerprint, "GetNonzeroElements")
                        else fingerprint.nonzero_elements()
                    ).items()
                )
            },
        }
    except Exception:  # noqa: BLE001
        return {"status": "error"}


def capture_bit_fingerprint(call: Any, *, unsigned_bits: bool = False) -> dict[str, Any]:
    try:
        fingerprint = call()
        if hasattr(fingerprint, "GetNumBits"):
            length = fingerprint.GetNumBits()
            on_bits = list(fingerprint.GetOnBits())
        else:
            length = fingerprint.n_bits() if hasattr(fingerprint, "n_bits") else fingerprint.size()
            on_bits = fingerprint.on_bits()
        if unsigned_bits:
            on_bits = sorted(int(bit) & UINT32_MASK for bit in on_bits)
        return {"status": "ok", "length": length, "on_bits": on_bits}
    except Exception:  # noqa: BLE001
        return {"status": "error"}


def compare_atom_pair(
    audit: Audit,
    record: dict[str, Any],
    rd_mol: Chem.Mol,
    ck_mol: Any,
    branches: list[dict[str, Any]],
    generators: dict[str, Any],
) -> None:
    for branch in branches:
        name = str(branch["name"])
        generator = generators[name]
        rdkit_kwargs, cosmolkit_kwargs = atom_pair_kwargs(
            branch, rd_mol.GetNumAtoms()
        )
        audit.compare(
            f"atom_pair.{name}.sparse_count",
            record,
            capture_count_fingerprint(
                lambda: generator.GetSparseCountFingerprint(rd_mol, **rdkit_kwargs)
            ),
            capture_count_fingerprint(
                lambda: ck_mol.fingerprint_atom_pair_sparse_count(**cosmolkit_kwargs)
            ),
        )
        audit.compare(
            f"atom_pair.{name}.count",
            record,
            capture_count_fingerprint(
                lambda: generator.GetCountFingerprint(rd_mol, **rdkit_kwargs)
            ),
            capture_count_fingerprint(
                lambda: ck_mol.fingerprint_atom_pair_count(**cosmolkit_kwargs)
            ),
        )
        audit.compare(
            f"atom_pair.{name}.sparse_bit",
            record,
            capture_bit_fingerprint(
                lambda: generator.GetSparseFingerprint(rd_mol, **rdkit_kwargs)
            ),
            capture_bit_fingerprint(
                lambda: ck_mol.fingerprint_atom_pair_sparse_bits(**cosmolkit_kwargs)
            ),
        )
        audit.compare(
            f"atom_pair.{name}.explicit_bit",
            record,
            capture_bit_fingerprint(
                lambda: generator.GetFingerprint(rd_mol, **rdkit_kwargs)
            ),
            capture_bit_fingerprint(
                lambda: ck_mol.fingerprint_atom_pair(**cosmolkit_kwargs)
            ),
        )
        if branch.get("additionalOutput"):
            try:
                rd_output = atom_pair_additional_output()
                rd_fingerprint = generator.GetFingerprint(
                    rd_mol, additionalOutput=rd_output, **rdkit_kwargs
                )
                rd_value: Any = {
                    "status": "ok",
                    "length": rd_fingerprint.GetNumBits(),
                    "on_bits": list(rd_fingerprint.GetOnBits()),
                    "additional_output": normalize_atom_pair_additional_output(
                        rd_output
                    ),
                }
            except Exception:  # noqa: BLE001
                rd_value = {"status": "error"}
            try:
                ck_result = ck_mol.fingerprint_atom_pair_with_output(
                    **cosmolkit_kwargs
                )
                ck_fingerprint = ck_result.fingerprint()
                ck_value: Any = {
                    "status": "ok",
                    "length": ck_fingerprint.n_bits(),
                    "on_bits": ck_fingerprint.on_bits(),
                    "additional_output": normalize_cosmolkit_atom_pair_additional_output(
                        ck_result.additional_output()
                    ),
                }
            except Exception:  # noqa: BLE001
                ck_value = {"status": "error"}
            audit.compare(
                f"atom_pair.{name}.additional_output",
                record,
                rd_value,
                ck_value,
            )


def topological_torsion_kwargs(
    branch: dict[str, Any], atom_count: int
) -> tuple[dict[str, Any], dict[str, Any]]:
    rdkit_kwargs: dict[str, Any] = {}
    cosmolkit_kwargs: dict[str, Any] = {}
    if branch.get("fromAtoms") == "first":
        atoms = [0] if atom_count else []
        rdkit_kwargs["fromAtoms"] = atoms
        cosmolkit_kwargs["from_atoms"] = atoms
    if branch.get("ignoreAtoms") == "first":
        atoms = [0] if atom_count else []
        rdkit_kwargs["ignoreAtoms"] = atoms
        cosmolkit_kwargs["ignore_atoms"] = atoms
    if branch.get("customAtomInvariants") == "index_plus_17":
        invariants = [index + 17 for index in range(atom_count)]
        rdkit_kwargs["customAtomInvariants"] = invariants
        cosmolkit_kwargs["custom_atom_invariants"] = invariants
    return rdkit_kwargs, cosmolkit_kwargs


def make_topological_torsion_generators(
    branches: list[dict[str, Any]],
) -> tuple[dict[str, Any], dict[str, Any]]:
    rdkit_generators: dict[str, Any] = {}
    cosmolkit_generators: dict[str, Any] = {}
    for branch in branches:
        name = str(branch["name"])
        rdkit_generator = rdFingerprintGenerator.GetTopologicalTorsionGenerator(
            includeChirality=branch.get("includeChirality", False),
            torsionAtomCount=branch.get("torsionAtomCount", 4),
            countSimulation=branch.get("countSimulation", True),
            countBounds=branch.get("countBounds"),
            fpSize=branch.get("fpSize", 2048),
        )
        rdkit_options = rdkit_generator.GetOptions()
        rdkit_options.onlyShortestPaths = branch.get("onlyShortestPaths", False)
        rdkit_options.numBitsPerFeature = branch.get("numBitsPerFeature", 1)

        cosmolkit_generator = cosmolkit.get_topological_torsion_generator(
            include_chirality=branch.get("includeChirality", False),
            torsion_atom_count=branch.get("torsionAtomCount", 4),
            count_simulation=branch.get("countSimulation", True),
            count_bounds=branch.get("countBounds"),
            fp_size=branch.get("fpSize", 2048),
        )
        cosmolkit_options = cosmolkit_generator.get_options()
        cosmolkit_options.only_shortest_paths = branch.get("onlyShortestPaths", False)
        cosmolkit_options.num_bits_per_feature = branch.get("numBitsPerFeature", 1)
        rdkit_generators[name] = rdkit_generator
        cosmolkit_generators[name] = cosmolkit_generator
    return rdkit_generators, cosmolkit_generators


def topological_torsion_additional_output() -> tuple[Any, Any]:
    rdkit_output = rdFingerprintGenerator.AdditionalOutput()
    cosmolkit_output = cosmolkit.AdditionalOutput()
    for rdkit_method, cosmolkit_method in (
        ("AllocateAtomCounts", "allocate_atom_counts"),
        ("AllocateAtomToBits", "allocate_atom_to_bits"),
        ("AllocateBitInfoMap", "allocate_bit_info_map"),
        ("AllocateBitPaths", "allocate_bit_paths"),
        ("AllocateAtomsPerBit", "allocate_atoms_per_bit"),
    ):
        getattr(rdkit_output, rdkit_method)()
        getattr(cosmolkit_output, cosmolkit_method)()
    return rdkit_output, cosmolkit_output


def normalize_topological_torsion_additional_output(output: Any) -> dict[str, Any]:
    if hasattr(output, "GetAtomCounts"):
        atom_counts = output.GetAtomCounts()
        atom_to_bits = output.GetAtomToBits()
        bit_info_map = output.GetBitInfoMap()
        bit_paths = output.GetBitPaths()
        atoms_per_bit = output.GetAtomsPerBit()
    else:
        atom_counts = output.atom_counts()
        atom_to_bits = output.atom_to_bits()
        bit_info_map = output.bit_info_map()
        bit_paths = output.bit_paths()
        atoms_per_bit = output.atoms_per_bit()
    return {
        "atom_counts": list(atom_counts),
        "atom_to_bits": [list(bits) for bits in atom_to_bits],
        "bit_info_map": {
            int(bit): [list(pair) for pair in pairs]
            for bit, pairs in sorted(bit_info_map.items())
        },
        "bit_paths": {
            int(bit): [list(path) for path in paths]
            for bit, paths in sorted(bit_paths.items())
        },
        "atoms_per_bit": {
            int(bit): [list(atoms) for atoms in entries]
            for bit, entries in sorted(atoms_per_bit.items())
        },
    }


def compare_topological_torsion(
    audit: Audit,
    record: dict[str, Any],
    rd_mol: Chem.Mol,
    ck_mol: Any,
    branches: list[dict[str, Any]],
    rdkit_generators: dict[str, Any],
    cosmolkit_generators: dict[str, Any],
) -> None:
    forms = (
        (
            "sparse_count",
            "GetSparseCountFingerprint",
            "get_sparse_count_fingerprint",
            capture_count_fingerprint,
            False,
        ),
        (
            "count",
            "GetCountFingerprint",
            "get_count_fingerprint",
            capture_count_fingerprint,
            False,
        ),
        (
            "sparse_bit",
            "GetSparseFingerprint",
            "get_sparse_fingerprint",
            capture_bit_fingerprint,
            True,
        ),
        (
            "explicit_bit",
            "GetFingerprint",
            "get_fingerprint",
            capture_bit_fingerprint,
            False,
        ),
    )
    for branch in branches:
        name = str(branch["name"])
        rdkit_generator = rdkit_generators[name]
        cosmolkit_generator = cosmolkit_generators[name]
        rdkit_kwargs, cosmolkit_kwargs = topological_torsion_kwargs(
            branch, rd_mol.GetNumAtoms()
        )
        for form, rdkit_method, cosmolkit_method, capture, unsigned_bits in forms:
            audit.compare(
                f"topological_torsion.{name}.{form}",
                record,
                capture(
                    lambda method=rdkit_method: getattr(rdkit_generator, method)(
                        rd_mol, **rdkit_kwargs
                    ),
                    **(
                        {"unsigned_bits": unsigned_bits}
                        if capture is capture_bit_fingerprint
                        else {}
                    ),
                ),
                capture(
                    lambda method=cosmolkit_method: getattr(cosmolkit_generator, method)(
                        ck_mol, **cosmolkit_kwargs
                    ),
                    **(
                        {"unsigned_bits": unsigned_bits}
                        if capture is capture_bit_fingerprint
                        else {}
                    ),
                ),
            )
            if branch.get("additionalOutput"):
                rdkit_output, cosmolkit_output = topological_torsion_additional_output()
                rdkit_value = capture(
                    lambda method=rdkit_method: getattr(rdkit_generator, method)(
                        rd_mol, additionalOutput=rdkit_output, **rdkit_kwargs
                    ),
                    **(
                        {"unsigned_bits": unsigned_bits}
                        if capture is capture_bit_fingerprint
                        else {}
                    ),
                )
                cosmolkit_value = capture(
                    lambda method=cosmolkit_method: getattr(cosmolkit_generator, method)(
                        ck_mol, additional_output=cosmolkit_output, **cosmolkit_kwargs
                    ),
                    **(
                        {"unsigned_bits": unsigned_bits}
                        if capture is capture_bit_fingerprint
                        else {}
                    ),
                )
                rdkit_value["additional_output"] = normalize_topological_torsion_additional_output(rdkit_output)
                cosmolkit_value["additional_output"] = normalize_topological_torsion_additional_output(cosmolkit_output)
                audit.compare(
                    f"topological_torsion.{name}.{form}.additional_output",
                    record,
                    rdkit_value,
                    cosmolkit_value,
                )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--topological-profile", type=Path, required=True)
    parser.add_argument("--avalon-profile", type=Path, required=True)
    parser.add_argument("--atom-pair-profile", type=Path, required=True)
    parser.add_argument("--topological-torsion-profile", type=Path, required=True)
    parser.add_argument("--layered-profile", type=Path, required=True)
    parser.add_argument(
        "--mode",
        choices=("topological_avalon", "atom_pair", "topological_torsion", "layered"),
        required=True,
    )
    parser.add_argument("--limit", type=int)
    parser.add_argument("--max-examples", type=int, default=12)
    args = parser.parse_args()

    if rdBase.rdkitVersion != EXPECTED_RDKIT_VERSION:
        raise SystemExit(
            f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, "
            f"got {rdBase.rdkitVersion}"
        )
    expected_cosmolkit = os.environ.get("COSMOLKIT_EXPECTED_VERSION")
    if expected_cosmolkit and cosmolkit.__version__ != expected_cosmolkit:
        raise SystemExit(
            "COSMolKit version mismatch: expected "
            f"{expected_cosmolkit}, got {cosmolkit.__version__}"
        )

    topological_profile = load_json(args.topological_profile)
    avalon_profile = load_json(args.avalon_profile)
    atom_pair_profile = load_json(args.atom_pair_profile)
    topological_torsion_profile = load_json(args.topological_torsion_profile)
    layered_profile = load_json(args.layered_profile)
    topological_branches = topological_profile["branches"]
    avalon_branches = avalon_profile["corpus_branches"]
    atom_pair_branches = atom_pair_profile["branches"]
    atom_pair_generators = make_atom_pair_generators(atom_pair_branches)
    topological_torsion_branches = topological_torsion_profile["corpus_branches"]
    complete_layered_branches = layered_branches(layered_profile)
    rdkit_torsion_generators, cosmolkit_torsion_generators = (
        make_topological_torsion_generators(topological_torsion_branches)
    )
    audit = Audit(args.output, args.max_examples)
    processed = 0
    try:
        with args.input.open(encoding="utf-8") as source:
            for line in source:
                if args.limit is not None and processed >= args.limit:
                    break
                record = json.loads(line)
                processed += 1
                try:
                    rd_mol = Chem.MolFromSmiles(record["smiles"], sanitize=True)
                except Exception as error:  # noqa: BLE001
                    rd_mol = None
                    rd_error: Any = error_value(error)
                else:
                    rd_error = None
                try:
                    ck_mol = cosmolkit.Molecule.from_smiles(record["smiles"])
                except Exception as error:  # noqa: BLE001
                    ck_mol = None
                    ck_error = error_value(error)
                else:
                    ck_error = None
                if rd_mol is None or ck_mol is None:
                    if rd_mol is None and ck_mol is None:
                        audit.count("parse.both_rejected")
                    else:
                        audit.mismatch(
                            "parse",
                            record,
                            "ok" if rd_mol is not None else rd_error or "rejected",
                            "ok" if ck_mol is not None else ck_error,
                        )
                    continue
                audit.count("parse.both_ok")
                if args.mode == "topological_avalon":
                    compare_topological(
                        audit, record, rd_mol, ck_mol, topological_branches
                    )
                    compare_topological_provenance(
                        audit, record, rd_mol, ck_mol, topological_branches
                    )
                    compare_avalon(audit, record, rd_mol, ck_mol, avalon_branches)
                elif args.mode == "atom_pair":
                    compare_atom_pair(
                        audit,
                        record,
                        rd_mol,
                        ck_mol,
                        atom_pair_branches,
                        atom_pair_generators,
                    )
                elif args.mode == "topological_torsion":
                    compare_topological_torsion(
                        audit,
                        record,
                        rd_mol,
                        ck_mol,
                        topological_torsion_branches,
                        rdkit_torsion_generators,
                        cosmolkit_torsion_generators,
                    )
                else:
                    compare_layered(
                        audit, record, rd_mol, ck_mol, complete_layered_branches
                    )
    finally:
        audit.finish(
            processed,
            (
                {
                    "topological": len(topological_branches),
                    "topological_provenance": 2,
                    "avalon": len(avalon_branches),
                }
                if args.mode == "topological_avalon"
                else {
                    "atom_pair_sparse_count": len(atom_pair_branches),
                    "atom_pair_count": len(atom_pair_branches),
                    "atom_pair_sparse_bit": len(atom_pair_branches),
                    "atom_pair_explicit_bit": len(atom_pair_branches),
                    "atom_pair_additional_output": sum(
                        bool(branch.get("additionalOutput"))
                        for branch in atom_pair_branches
                    ),
                }
                if args.mode == "atom_pair"
                else {
                    "topological_torsion_sparse_count": len(topological_torsion_branches),
                    "topological_torsion_count": len(topological_torsion_branches),
                    "topological_torsion_sparse_bit": len(topological_torsion_branches),
                    "topological_torsion_explicit_bit": len(topological_torsion_branches),
                    "topological_torsion_sparse_count_additional_output": sum(
                        bool(branch.get("additionalOutput"))
                        for branch in topological_torsion_branches
                    ),
                    "topological_torsion_count_additional_output": sum(
                        bool(branch.get("additionalOutput"))
                        for branch in topological_torsion_branches
                    ),
                    "topological_torsion_sparse_bit_additional_output": sum(
                        bool(branch.get("additionalOutput"))
                        for branch in topological_torsion_branches
                    ),
                    "topological_torsion_explicit_bit_additional_output": sum(
                        bool(branch.get("additionalOutput"))
                        for branch in topological_torsion_branches
                    ),
                }
                if args.mode == "topological_torsion"
                else {"layered": len(complete_layered_branches)}
            ),
        )


if __name__ == "__main__":
    main()
