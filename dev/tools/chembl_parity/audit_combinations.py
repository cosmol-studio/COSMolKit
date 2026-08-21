#!/usr/bin/env python3
# pyright: reportAttributeAccessIssue=false
"""Combination, order, cache, batch, and shared-object stress audit."""

from __future__ import annotations

import argparse
import hashlib
import json
import struct
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Callable

import cosmolkit
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, MACCSkeys, QED, rdDistGeom
from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors

from .audit_core import (
    Audit,
    ck_atom_state,
    ck_bond_state,
    compare_value,
    f64_bits,
    rd_atom_state,
    rd_bond_state,
    validate_runtime_versions,
)
from .audit_surfaces import normalize_svg, rd_svg


RDLogger.DisableLog("rdApp.*")


def iter_records(path: Path, shard: int, shards: int, limit: int | None):
    selected = 0
    with path.open(encoding="utf-8") as source:
        for index, line in enumerate(source):
            if index % shards != shard:
                continue
            if limit is not None and selected >= limit:
                break
            selected += 1
            yield json.loads(line)


def digest_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def digest_array(value: Any) -> str:
    array = np.ascontiguousarray(value, dtype=np.float64)
    digest = hashlib.sha256()
    digest.update(struct.pack(">I", array.ndim))
    for dimension in array.shape:
        digest.update(struct.pack(">Q", int(dimension)))
    digest.update(array.astype(">f8", copy=False).tobytes())
    return digest.hexdigest()


def rd_graph_state(mol: Chem.Mol) -> dict[str, Any]:
    return {
        "atoms": [rd_atom_state(atom) for atom in mol.GetAtoms()],
        "bonds": [rd_bond_state(bond) for bond in mol.GetBonds()],
    }


def ck_graph_state(mol: Any) -> dict[str, Any]:
    return {
        "atoms": [ck_atom_state(atom) for atom in mol.atoms()],
        "bonds": [ck_bond_state(bond) for bond in mol.bonds()],
    }


def rd_observations(mol: Chem.Mol) -> dict[str, Any]:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        inchi = Chem.MolToInchi(mol)
        key = Chem.MolToInchiKey(mol)
    return {
        "graph": rd_graph_state(mol),
        "smiles": Chem.MolToSmiles(mol),
        "inchi": inchi,
        "inchi_key": key,
        "mol_wt": f64_bits(Descriptors.MolWt(Chem.Mol(mol))),
        "qed": f64_bits(QED.qed(Chem.Mol(mol))),
        "tpsa": f64_bits(rdMolDescriptors.CalcTPSA(Chem.Mol(mol))),
        "morgan": list(generator.GetFingerprint(mol).GetOnBits()),
        "maccs": [
            bit - 1 for bit in MACCSkeys.GenMACCSKeys(mol).GetOnBits() if bit > 0
        ],
        "svg": digest_text(normalize_svg(rd_svg(mol, 300, 300))),
        "explicit_h": rd_graph_state(Chem.AddHs(mol)),
    }


def ck_observations(mol: Any) -> dict[str, Any]:
    return {
        "graph": ck_graph_state(mol),
        "smiles": mol.to_smiles(),
        "inchi": mol.to_inchi(),
        "inchi_key": mol.to_inchi_key(),
        "mol_wt": f64_bits(cosmolkit.calc_mol_wt(mol)),
        "qed": f64_bits(cosmolkit.calc_qed(mol)),
        "tpsa": f64_bits(cosmolkit.calc_tpsa(mol)),
        "morgan": mol.fingerprint_morgan().on_bits(),
        "maccs": mol.maccs_fingerprint().on_bits(),
        "svg": digest_text(normalize_svg(mol.to_svg(300, 300))),
        "explicit_h": ck_graph_state(mol.with_hydrogens()),
    }


def ck_operation_functions(mol: Any) -> dict[str, Callable[[], Any]]:
    return {
        "smiles": mol.to_smiles,
        "inchi": mol.to_inchi,
        "inchi_key": mol.to_inchi_key,
        "mol_wt": lambda: f64_bits(cosmolkit.calc_mol_wt(mol)),
        "qed": lambda: f64_bits(cosmolkit.calc_qed(mol)),
        "tpsa": lambda: f64_bits(cosmolkit.calc_tpsa(mol)),
        "morgan": lambda: mol.fingerprint_morgan().on_bits(),
        "maccs": lambda: mol.maccs_fingerprint().on_bits(),
        "bounds": lambda: digest_array(mol.dg_bounds_matrix()),
        "svg": lambda: digest_text(normalize_svg(mol.to_svg(300, 300))),
        "explicit_h": lambda: ck_graph_state(mol.with_hydrogens()),
    }


def run_ck_sequence(mol: Any, order: tuple[str, ...]) -> dict[str, Any]:
    functions = ck_operation_functions(mol)
    return {name: functions[name]() for name in order}


def compare_bounds(audit: Audit, record: dict[str, Any], rd_mol: Chem.Mol, ck_mol: Any):
    rd_bounds = np.asarray(rdDistGeom.GetMoleculeBoundsMatrix(rd_mol), dtype=np.float64)
    ck_bounds = np.asarray(ck_mol.dg_bounds_matrix(), dtype=np.float64)
    if rd_bounds.shape != ck_bounds.shape:
        audit.mismatch("combination.bounds.shape", record, rd_bounds.shape, ck_bounds.shape)
        return
    max_error = float(np.abs(rd_bounds - ck_bounds).max(initial=0.0))
    audit.count("combination.bounds.entries", rd_bounds.size)
    if max_error <= 1.0e-8:
        audit.count("match.combination.bounds")
    else:
        audit.mismatch(
            "combination.bounds", record, "within 1e-8", max_error, f"max_error={max_error}"
        )


def compare_cache_sequences(audit: Audit, record: dict[str, Any], rd_mol: Chem.Mol, ck_mol: Any):
    rd_crippen_mol = Chem.Mol(rd_mol)
    rd_crippen = [
        f64_bits(Crippen.MolLogP(rd_crippen_mol, includeHs=False, force=True)),
        f64_bits(Crippen.MolLogP(rd_crippen_mol, includeHs=True, force=False)),
    ]
    ck_crippen = [
        f64_bits(cosmolkit.calc_crippen_descriptors(ck_mol, include_hs=False, force=True)[0]),
        f64_bits(cosmolkit.calc_crippen_descriptors(ck_mol, include_hs=True, force=False)[0]),
    ]
    compare_value(
        audit, "combination.crippen_cache_sequence", record, rd_crippen, ck_crippen
    )

    rd_tpsa_mol = Chem.Mol(rd_mol)
    rd_tpsa = [
        f64_bits(rdMolDescriptors.CalcTPSA(rd_tpsa_mol, force=True, includeSandP=False)),
        f64_bits(rdMolDescriptors.CalcTPSA(rd_tpsa_mol, force=False, includeSandP=True)),
    ]
    ck_tpsa = [
        f64_bits(cosmolkit.calc_tpsa(ck_mol, force=True, include_sandp=False)),
        f64_bits(cosmolkit.calc_tpsa(ck_mol, force=False, include_sandp=True)),
    ]
    compare_value(audit, "combination.tpsa_cache_sequence", record, rd_tpsa, ck_tpsa)


def compare_morgan_option_sequence(
    audit: Audit, record: dict[str, Any], rd_mol: Chem.Mol, ck_mol: Any
):
    rd_default = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    rd_two = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    rd_two.GetOptions().numBitsPerFeature = 2
    rd_values = [
        list(rd_default.GetFingerprint(rd_mol).GetOnBits()),
        list(rd_two.GetFingerprint(rd_mol).GetOnBits()),
        list(rd_default.GetFingerprint(rd_mol).GetOnBits()),
    ]
    ck_values = [
        ck_mol.fingerprint_morgan().on_bits(),
        ck_mol.fingerprint_morgan(num_bits_per_feature=2).on_bits(),
        ck_mol.fingerprint_morgan().on_bits(),
    ]
    compare_value(audit, "combination.morgan_option_sequence", record, rd_values, ck_values)


def audit_scalar(audit: Audit, record: dict[str, Any]):
    rd_mol = Chem.MolFromSmiles(record["smiles"])
    if rd_mol is None:
        audit.count("skip.rdkit_parse")
        return
    try:
        ck_mol = cosmolkit.Molecule.from_smiles(record["smiles"])
    except Exception as error:
        audit.mismatch("combination.parse", record, "success", repr(error))
        return

    rd_values = rd_observations(rd_mol)
    ck_values = ck_observations(ck_mol)
    for name, rd_value in rd_values.items():
        compare_value(audit, f"combination.scalar.{name}", record, rd_value, ck_values[name])
    compare_bounds(audit, record, rd_mol, ck_mol)

    functions = ck_operation_functions(ck_mol)
    baseline = {name: function() for name, function in functions.items()}
    names = tuple(functions)
    orders = (
        names,
        tuple(reversed(names)),
        names[::2] + names[1::2],
        names[1::2] + names[::2],
    )
    graph_before = ck_graph_state(ck_mol)
    for order_index, order in enumerate(orders):
        observed = run_ck_sequence(ck_mol, order)
        for name in order:
            compare_value(
                audit,
                f"combination.order_{order_index}.{name}",
                record,
                baseline[name],
                observed[name],
            )
        compare_value(
            audit,
            f"combination.order_{order_index}.source_graph",
            record,
            graph_before,
            ck_graph_state(ck_mol),
        )

    rd_roundtrip = Chem.RemoveHs(Chem.AddHs(rd_mol))
    ck_roundtrip = ck_mol.with_hydrogens().without_hydrogens()
    compare_value(
        audit,
        "combination.add_remove_h.graph",
        record,
        rd_graph_state(rd_roundtrip),
        ck_graph_state(ck_roundtrip),
    )
    compare_cache_sequences(audit, record, rd_mol, ck_mol)
    compare_morgan_option_sequence(audit, record, rd_mol, ck_mol)


def audit_concurrent(audit: Audit, record: dict[str, Any]):
    try:
        mol = cosmolkit.Molecule.from_smiles(record["smiles"])
    except Exception as error:
        audit.mismatch("concurrent.parse", record, "success", repr(error))
        return
    baseline = ck_observations(mol)

    def worker(_: int) -> tuple[dict[str, Any], dict[str, Any]]:
        return ck_observations(mol), ck_graph_state(mol)

    try:
        with ThreadPoolExecutor(max_workers=8) as executor:
            results = list(executor.map(worker, range(16)))
    except Exception as error:
        audit.mismatch("concurrent.shared_object_error", record, "success", repr(error))
        return
    for result, graph in results:
        compare_value(audit, "concurrent.shared_object_values", record, baseline, result)
        compare_value(audit, "concurrent.shared_object_graph", record, baseline["graph"], graph)


def batch_compare_value(
    audit: Audit,
    feature: str,
    records: list[dict[str, Any]],
    expected: list[Any],
    actual: list[Any],
):
    if len(expected) != len(actual):
        audit.mismatch(feature + ".length", records[0], len(expected), len(actual))
        return
    for record, expected_value, actual_value in zip(records, expected, actual):
        compare_value(audit, feature, record, expected_value, actual_value)


def audit_batch_chunk(audit: Audit, records: list[dict[str, Any]]):
    smiles = [record["smiles"] for record in records]
    scalar: list[Any] = []
    valid_records: list[dict[str, Any]] = []
    valid_smiles: list[str] = []
    for record, text in zip(records, smiles):
        try:
            scalar.append(cosmolkit.Molecule.from_smiles(text))
            valid_records.append(record)
            valid_smiles.append(text)
        except Exception as error:
            audit.mismatch("batch.scalar_parse", record, "success", repr(error))
    if not scalar:
        return
    batch = cosmolkit.MoleculeBatch.from_smiles_list(
        valid_smiles, errors="raise", n_jobs=8
    )
    raw_scalar = [
        cosmolkit.Molecule.from_smiles(text, sanitize=False) for text in valid_smiles
    ]
    raw_batch = cosmolkit.MoleculeBatch.from_smiles_list(
        valid_smiles, sanitize=False, errors="raise", n_jobs=8
    )
    compare_value(
        audit,
        "batch.valid_count",
        valid_records[0],
        len(scalar),
        batch.valid_count(),
    )
    expected_smiles = [mol.to_smiles() for mol in scalar]
    expected_fps = [mol.fingerprint_morgan().on_bits() for mol in scalar]
    expected_outputs = [
        {
            "fingerprint": result.fingerprint().on_bits(),
            "atom_counts": result.additional_output().atom_counts(),
            "atom_to_bits": result.additional_output().atom_to_bits(),
            "bit_info_map": result.additional_output().bit_info_map(),
            "atoms_per_bit": result.additional_output().atoms_per_bit(),
        }
        for result in (mol.fingerprint_morgan_with_output() for mol in scalar)
    ]
    expected_bounds = [digest_array(mol.dg_bounds_matrix()) for mol in scalar]
    expected_svgs = [
        digest_text(normalize_svg(mol.to_svg(300, 300))) for mol in scalar
    ]
    expected_h_roundtrip = [
        mol.with_hydrogens().without_hydrogens().to_smiles() for mol in scalar
    ]
    expected_sanitized = [mol.sanitize().to_smiles() for mol in raw_scalar]
    expected_kekulized = {
        clear_aromatic_flags: [
            mol.with_kekulized_bonds(clear_aromatic_flags).to_smiles(kekule=True)
            for mol in scalar
        ]
        for clear_aromatic_flags in (False, True)
    }
    expected_2d = [
        digest_array(mol.with_2d_coordinates().coordinates_2d()) for mol in scalar
    ]
    for jobs in (1, 8):
        batch_compare_value(
            audit,
            f"batch.jobs_{jobs}.smiles",
            valid_records,
            expected_smiles,
            batch.to_smiles_list(n_jobs=jobs),
        )
        batch_compare_value(
            audit,
            f"batch.jobs_{jobs}.morgan",
            valid_records,
            expected_fps,
            [
                fp.on_bits() if fp is not None else None
                for fp in batch.fingerprint_morgan_list(n_jobs=jobs)
            ],
        )
        batch_outputs = []
        for result in batch.fingerprint_morgan_with_output_list(n_jobs=jobs):
            if result is None:
                batch_outputs.append(None)
                continue
            output = result.additional_output()
            batch_outputs.append(
                {
                    "fingerprint": result.fingerprint().on_bits(),
                    "atom_counts": output.atom_counts(),
                    "atom_to_bits": output.atom_to_bits(),
                    "bit_info_map": output.bit_info_map(),
                    "atoms_per_bit": output.atoms_per_bit(),
                }
            )
        batch_compare_value(
            audit,
            f"batch.jobs_{jobs}.morgan_output",
            valid_records,
            expected_outputs,
            batch_outputs,
        )
        batch_compare_value(
            audit,
            f"batch.jobs_{jobs}.bounds",
            valid_records,
            expected_bounds,
            [
                digest_array(value) if value is not None else None
                for value in batch.dg_bounds_matrix_list(n_jobs=jobs)
            ],
        )
        batch_compare_value(
            audit,
            f"batch.jobs_{jobs}.svg",
            valid_records,
            expected_svgs,
            [
                digest_text(normalize_svg(value)) if value is not None else None
                for value in batch.to_svg_list(n_jobs=jobs)
            ],
        )
        transformed = batch.with_hydrogens(n_jobs=jobs).without_hydrogens(n_jobs=jobs)
        batch_compare_value(
            audit,
            f"batch.jobs_{jobs}.add_remove_h",
            valid_records,
            expected_h_roundtrip,
            transformed.to_smiles_list(n_jobs=jobs),
        )
        batch_compare_value(
            audit,
            f"batch.jobs_{jobs}.sanitize",
            valid_records,
            expected_sanitized,
            raw_batch.sanitize(errors="raise", n_jobs=jobs).to_smiles_list(),
        )
        for clear_aromatic_flags in (False, True):
            branch = f"clear_aromatic_flags_{str(clear_aromatic_flags).lower()}"
            batch_compare_value(
                audit,
                f"batch.jobs_{jobs}.kekulize.{branch}",
                valid_records,
                expected_kekulized[clear_aromatic_flags],
                batch.with_kekulized_bonds(
                    clear_aromatic_flags, errors="raise", n_jobs=jobs
                ).to_smiles_list(kekule=True),
            )
        batch_2d = batch.with_2d_coordinates(errors="raise", n_jobs=jobs).to_list()
        batch_compare_value(
            audit,
            f"batch.jobs_{jobs}.coordinates_2d",
            valid_records,
            expected_2d,
            [
                digest_array(molecule.coordinates_2d())
                if molecule is not None
                else None
                for molecule in batch_2d
            ],
        )


def main() -> None:
    validate_runtime_versions()
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--mode", choices=("scalar", "concurrent", "batch"), required=True
    )
    parser.add_argument("--limit", type=int)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--max-atoms", type=int, default=80)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--max-examples", type=int, default=300)
    args = parser.parse_args()

    audit = Audit(args.output, args.max_examples)
    processed = 0
    batch_records: list[dict[str, Any]] = []
    try:
        for record in iter_records(args.input, args.shard, args.shards, args.limit):
            try:
                rd_mol = Chem.MolFromSmiles(record["smiles"])
            except Exception:
                rd_mol = None
            if rd_mol is None:
                try:
                    cosmolkit.Molecule.from_smiles(record["smiles"])
                except Exception:
                    audit.count("skip.both_parse_rejected")
                else:
                    audit.mismatch(
                        "combination.parse", record, "rejected", "COSMolKit accepted"
                    )
                continue
            if rd_mol.GetNumAtoms() > args.max_atoms:
                audit.count("skip.selection")
                continue
            processed += 1
            if args.mode == "scalar":
                audit_scalar(audit, record)
            elif args.mode == "concurrent":
                audit_concurrent(audit, record)
            else:
                batch_records.append(record)
                if len(batch_records) >= args.batch_size:
                    audit_batch_chunk(audit, batch_records)
                    batch_records.clear()
            if processed % 100 == 0:
                print(
                    json.dumps(
                        {
                            "mode": args.mode,
                            "processed": processed,
                            "mismatches": sum(
                                value
                                for key, value in audit.counts.items()
                                if key.startswith("mismatch.")
                            ),
                        }
                    ),
                    flush=True,
                )
        if batch_records:
            audit_batch_chunk(audit, batch_records)
    finally:
        audit.finish(processed, args.mode, args.shard, args.shards)


if __name__ == "__main__":
    main()
