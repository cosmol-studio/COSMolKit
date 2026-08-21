#!/usr/bin/env python3
# pyright: reportAttributeAccessIssue=false
"""Exact ChEMBL parity audit for topological and Avalon fingerprints."""

from __future__ import annotations

import argparse
import json
import os
import time
from collections import Counter
from pathlib import Path
from typing import Any

import cosmolkit
from rdkit import Chem, RDLogger, rdBase
from rdkit.Avalon import pyAvalonTools


EXPECTED_RDKIT_VERSION = "2026.03.1"
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--topological-profile", type=Path, required=True)
    parser.add_argument("--avalon-profile", type=Path, required=True)
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
    topological_branches = topological_profile["branches"]
    avalon_branches = avalon_profile["corpus_branches"]
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
                compare_topological(
                    audit, record, rd_mol, ck_mol, topological_branches
                )
                compare_topological_provenance(
                    audit, record, rd_mol, ck_mol, topological_branches
                )
                compare_avalon(audit, record, rd_mol, ck_mol, avalon_branches)
    finally:
        audit.finish(
            processed,
            {
                "topological": len(topological_branches),
                "topological_provenance": 2,
                "avalon": len(avalon_branches),
            },
        )


if __name__ == "__main__":
    main()
