#!/usr/bin/env python3
"""Generate RDKit UFF/MMFF molecule-parameter coverage golden data."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
from typing import Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

EXPECTED_RDKIT_VERSION = "2026.3.1"
FORCEFIELD_PARITY_SEED = 61453
FORCEFIELD_NONBONDED_THRESH = 100.0
FORCEFIELD_OPT_MAX_ITERS = 200
FORCEFIELD_MULTI_CONFS = 2


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


def forcefield_result(fn, mol: Chem.Mol) -> dict[str, object]:
    try:
        probe_mol = Chem.Mol(mol)
        return {"ok": True, "has_all": bool(fn(probe_mol)), "error": None}
    except Exception as err:
        return {"ok": False, "has_all": None, "error": str(err)}


def mmff_result(mol: Chem.Mol) -> dict[str, object]:
    try:
        probe_mol = Chem.Mol(mol)
        has_all = bool(AllChem.MMFFHasAllMoleculeParams(probe_mol))
        props = AllChem.MMFFGetMoleculeProperties(probe_mol)
        atom_types = (
            [int(props.GetMMFFAtomType(i)) for i in range(probe_mol.GetNumAtoms())]
            if props is not None
            else None
        )
        formal_charges = (
            [
                float(props.GetMMFFFormalCharge(i))
                for i in range(probe_mol.GetNumAtoms())
            ]
            if props is not None
            else None
        )
        partial_charges = (
            [
                float(props.GetMMFFPartialCharge(i))
                for i in range(probe_mol.GetNumAtoms())
            ]
            if props is not None
            else None
        )
        return {
            "ok": True,
            "has_all": has_all,
            "atom_types": atom_types,
            "formal_charges": formal_charges,
            "partial_charges": partial_charges,
            "error": None,
        }
    except Exception as err:
        return {
            "ok": False,
            "has_all": None,
            "atom_types": None,
            "formal_charges": None,
            "partial_charges": None,
            "error": str(err),
        }


def embedded_forcefield_result(smiles: str) -> dict[str, object]:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return {
            "ok": False,
            "cxsmiles": None,
            "coords": None,
            "uff": {"ok": False, "needs_more": None, "energy": None, "error": "MolFromSmiles returned None"},
            "mmff": {"ok": False, "needs_more": None, "energy": None, "error": "MolFromSmiles returned None"},
            "uff_optimized": {
                "ok": False,
                "needs_more": None,
                "energy": None,
                "coords": None,
                "error": "MolFromSmiles returned None",
            },
            "mmff_optimized": {
                "ok": False,
                "needs_more": None,
                "energy": None,
                "coords": None,
                "error": "MolFromSmiles returned None",
            },
            "uff_multi_optimized": {
                "ok": False,
                "conformer_results": None,
                "initial_coords": None,
                "error": "MolFromSmiles returned None",
            },
            "mmff_multi_optimized": {
                "ok": False,
                "conformer_results": None,
                "initial_coords": None,
                "error": "MolFromSmiles returned None",
            },
            "error": "MolFromSmiles returned None",
        }
    try:
        status = AllChem.EmbedMolecule(
            mol,
            randomSeed=FORCEFIELD_PARITY_SEED,
            useRandomCoords=True,
        )
        if status != 0:
            return {
                "ok": False,
                "cxsmiles": None,
                "coords": None,
                "uff": {
                    "ok": False,
                    "needs_more": None,
                    "energy": None,
                    "error": f"EmbedMolecule returned {status}",
                },
                "mmff": {
                    "ok": False,
                    "needs_more": None,
                    "energy": None,
                    "error": f"EmbedMolecule returned {status}",
                },
                "uff_optimized": {
                    "ok": False,
                    "needs_more": None,
                    "energy": None,
                    "coords": None,
                    "error": f"EmbedMolecule returned {status}",
                },
                "mmff_optimized": {
                    "ok": False,
                    "needs_more": None,
                    "energy": None,
                    "coords": None,
                    "error": f"EmbedMolecule returned {status}",
                },
                "uff_multi_optimized": {
                    "ok": False,
                    "conformer_results": None,
                    "initial_coords": None,
                    "error": f"EmbedMolecule returned {status}",
                },
                "mmff_multi_optimized": {
                    "ok": False,
                    "conformer_results": None,
                    "initial_coords": None,
                    "error": f"EmbedMolecule returned {status}",
                },
                "error": f"EmbedMolecule returned {status}",
            }
        cxsmiles = Chem.MolToCXSmiles(mol, isomericSmiles=True)
        parity_mol = Chem.MolFromSmiles(cxsmiles, sanitize=True)
        if parity_mol is None:
            raise ValueError("RDKit failed to parse its own forcefield CXSMILES")
        conf = parity_mol.GetConformer(0)
        coords = [
            [float(conf.GetAtomPosition(i)[axis]) for axis in range(3)]
            for i in range(parity_mol.GetNumAtoms())
        ]
        return {
            "ok": True,
            "cxsmiles": cxsmiles,
            "coords": coords,
            "uff": forcefield_initial_energy_result(
                lambda m: AllChem.UFFGetMoleculeForceField(
                    m,
                    vdwThresh=FORCEFIELD_NONBONDED_THRESH,
                    confId=0,
                    ignoreInterfragInteractions=True,
                ),
                parity_mol,
            ),
            "mmff": forcefield_initial_energy_result(
                lambda m: AllChem.MMFFGetMoleculeForceField(
                    m,
                    AllChem.MMFFGetMoleculeProperties(m, mmffVariant="MMFF94"),
                    nonBondedThresh=FORCEFIELD_NONBONDED_THRESH,
                    confId=0,
                    ignoreInterfragInteractions=True,
                ),
                parity_mol,
            ),
            "uff_optimized": forcefield_optimized_result(
                lambda m: AllChem.UFFGetMoleculeForceField(
                    m,
                    vdwThresh=FORCEFIELD_NONBONDED_THRESH,
                    confId=0,
                    ignoreInterfragInteractions=True,
                ),
                parity_mol,
                FORCEFIELD_OPT_MAX_ITERS,
            ),
            "mmff_optimized": forcefield_optimized_result(
                lambda m: AllChem.MMFFGetMoleculeForceField(
                    m,
                    AllChem.MMFFGetMoleculeProperties(m, mmffVariant="MMFF94"),
                    nonBondedThresh=FORCEFIELD_NONBONDED_THRESH,
                    confId=0,
                    ignoreInterfragInteractions=True,
                ),
                parity_mol,
                FORCEFIELD_OPT_MAX_ITERS,
            ),
            "uff_multi_optimized": forcefield_multi_optimized_result(
                lambda m: AllChem.UFFOptimizeMoleculeConfs(
                    m,
                    numThreads=1,
                    maxIters=FORCEFIELD_OPT_MAX_ITERS,
                    vdwThresh=FORCEFIELD_NONBONDED_THRESH,
                    ignoreInterfragInteractions=True,
                ),
                parity_mol,
                FORCEFIELD_MULTI_CONFS,
            ),
            "mmff_multi_optimized": forcefield_multi_optimized_result(
                lambda m: AllChem.MMFFOptimizeMoleculeConfs(
                    m,
                    numThreads=1,
                    maxIters=FORCEFIELD_OPT_MAX_ITERS,
                    mmffVariant="MMFF94",
                    nonBondedThresh=FORCEFIELD_NONBONDED_THRESH,
                    ignoreInterfragInteractions=True,
                ),
                parity_mol,
                FORCEFIELD_MULTI_CONFS,
            ),
            "error": None,
        }
    except Exception as err:
        return {
            "ok": False,
            "cxsmiles": None,
            "coords": None,
            "uff": {"ok": False, "needs_more": None, "energy": None, "error": str(err)},
            "mmff": {"ok": False, "needs_more": None, "energy": None, "error": str(err)},
            "uff_optimized": {
                "ok": False,
                "needs_more": None,
                "energy": None,
                "coords": None,
                "error": str(err),
            },
            "mmff_optimized": {
                "ok": False,
                "needs_more": None,
                "energy": None,
                "coords": None,
                "error": str(err),
            },
            "uff_multi_optimized": {
                "ok": False,
                "conformer_results": None,
                "initial_coords": None,
                "error": str(err),
            },
            "mmff_multi_optimized": {
                "ok": False,
                "conformer_results": None,
                "initial_coords": None,
                "error": str(err),
            },
            "error": str(err),
        }


def forcefield_initial_energy_result(fn, mol: Chem.Mol) -> dict[str, object]:
    try:
        probe_mol = Chem.Mol(mol)
        ff = fn(probe_mol)
        if ff is None:
            return {
                "ok": True,
                "needs_more": -1,
                "energy": -1.0,
                "error": None,
            }
        ff.Initialize()
        needs_more = int(ff.Minimize(maxIts=0))
        return {
            "ok": True,
            "needs_more": needs_more,
            "energy": float(ff.CalcEnergy()),
            "gradient": [float(value) for value in ff.CalcGrad()],
            "error": None,
        }
    except Exception as err:
        return {
            "ok": False,
            "needs_more": None,
            "energy": None,
            "gradient": None,
            "error": str(err),
        }


def forcefield_optimized_result(fn, mol: Chem.Mol, max_iters: int) -> dict[str, object]:
    try:
        opt_mol = Chem.Mol(mol)
        ff = fn(opt_mol)
        if ff is None:
            return {
                "ok": True,
                "needs_more": -1,
                "energy": -1.0,
                "coords": None,
                "error": None,
            }
        ff.Initialize()
        needs_more = int(ff.Minimize(maxIts=max_iters))
        conf = opt_mol.GetConformer(0)
        coords = [
            [float(conf.GetAtomPosition(i)[axis]) for axis in range(3)]
            for i in range(opt_mol.GetNumAtoms())
        ]
        return {
            "ok": True,
            "needs_more": needs_more,
            "energy": float(ff.CalcEnergy()),
            "coords": coords,
            "error": None,
        }
    except Exception as err:
        return {
            "ok": False,
            "needs_more": None,
            "energy": None,
            "coords": None,
            "error": str(err),
        }


def forcefield_multi_optimized_result(fn, mol: Chem.Mol, num_confs: int) -> dict[str, object]:
    try:
        opt_mol = Chem.Mol(mol)
        opt_mol.RemoveAllConformers()
        conf_ids = AllChem.EmbedMultipleConfs(
            opt_mol,
            numConfs=num_confs,
            randomSeed=FORCEFIELD_PARITY_SEED,
            useRandomCoords=True,
        )
        if len(conf_ids) != num_confs:
            return {
                "ok": False,
                "conformer_results": None,
                "initial_coords": None,
                "error": f"EmbedMultipleConfs returned {len(conf_ids)} conformers",
            }
        initial_coords = [
            [
                [float(conf.GetAtomPosition(i)[axis]) for axis in range(3)]
                for i in range(opt_mol.GetNumAtoms())
            ]
            for conf in opt_mol.GetConformers()
        ]
        results = fn(opt_mol)
        conformer_results = []
        for conf, result in zip(opt_mol.GetConformers(), results):
            conformer_results.append(
                {
                    "ok": True,
                    "needs_more": int(result[0]),
                    "energy": float(result[1]),
                    "coords": [
                        [float(conf.GetAtomPosition(i)[axis]) for axis in range(3)]
                        for i in range(opt_mol.GetNumAtoms())
                    ],
                    "error": None,
                }
            )
        return {
            "ok": True,
            "conformer_results": conformer_results,
            "initial_coords": initial_coords,
            "error": None,
        }
    except Exception as err:
        return {
            "ok": False,
            "conformer_results": None,
            "initial_coords": None,
            "error": str(err),
        }


def build_record(smiles: str) -> dict[str, object]:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return {
            "smiles": smiles,
            "rdkit_ok": False,
            "uff": {"ok": False, "has_all": None, "error": "MolFromSmiles returned None"},
            "mmff": {
                "ok": False,
                "has_all": None,
                "atom_types": None,
                "formal_charges": None,
                "partial_charges": None,
                "error": "MolFromSmiles returned None",
            },
            "uff_explicit_h": {
                "ok": False,
                "has_all": None,
                "error": "MolFromSmiles returned None",
            },
            "mmff_explicit_h": {
                "ok": False,
                "has_all": None,
                "atom_types": None,
                "formal_charges": None,
                "partial_charges": None,
                "error": "MolFromSmiles returned None",
            },
            "embedded": embedded_forcefield_result(smiles),
            "error": "MolFromSmiles returned None",
        }

    explicit_h_mol = Chem.AddHs(mol)
    return {
        "smiles": smiles,
        "rdkit_ok": True,
        "uff": forcefield_result(AllChem.UFFHasAllMoleculeParams, mol),
        "mmff": mmff_result(mol),
        "uff_explicit_h": forcefield_result(
            AllChem.UFFHasAllMoleculeParams, explicit_h_mol
        ),
        "mmff_explicit_h": mmff_result(explicit_h_mol),
        "embedded": embedded_forcefield_result(smiles),
        "error": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("testdata/smiles/corpus/smiles_small.smi"),
        help="input SMILES file (default: testdata/smiles/corpus/smiles_small.smi)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("testdata/forcefield/expected/rdkit/smiles_small/forcefield_params.jsonl"),
        help="output JSONL path (default: testdata/forcefield/expected/rdkit/smiles_small/forcefield_params.jsonl)",
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
