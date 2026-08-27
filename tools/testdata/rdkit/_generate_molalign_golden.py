#!/usr/bin/env python3
"""Generate complete focused RDKit MolAlign observable-state records."""

from __future__ import annotations

import argparse
import copy
import json
import math
from importlib.metadata import version
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import rdMolAlign


EXPECTED_RDKIT_VERSION = "2026.3.1"


def assert_rdkit_version() -> None:
    actual = version("rdkit")
    if actual != EXPECTED_RDKIT_VERSION:
        raise SystemExit(
            f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, got {actual}"
        )


def molecule_from_source(source: dict[str, Any]) -> Chem.Mol:
    molecule = Chem.MolFromSmiles(source["smiles"])
    if molecule is None:
        raise ValueError(f"failed to parse {source['smiles']!r}")
    molecule.RemoveAllConformers()
    for item in source["conformers"]:
        coordinates = item["coordinates"]
        if len(coordinates) != molecule.GetNumAtoms():
            raise ValueError(
                f"conformer {item['id']} has {len(coordinates)} coordinates for "
                f"{molecule.GetNumAtoms()} atoms"
            )
        conformer = Chem.Conformer(molecule.GetNumAtoms())
        conformer.SetId(int(item["id"]))
        conformer.Set3D(True)
        for atom_index, point in enumerate(coordinates):
            conformer.SetAtomPosition(atom_index, tuple(float(value) for value in point))
        molecule.AddConformer(conformer, assignId=False)
    return molecule


def snapshot(molecule: Chem.Mol) -> list[dict[str, Any]]:
    return [
        {
            "id": conformer.GetId(),
            "is_3d": conformer.Is3D(),
            "coordinates": [
                list(conformer.GetAtomPosition(atom_index))
                for atom_index in range(molecule.GetNumAtoms())
            ],
        }
        for conformer in molecule.GetConformers()
    ]


def atom_map(value: Any) -> list[tuple[int, int]]:
    return [(int(pair[0]), int(pair[1])) for pair in value]


def atom_maps(value: Any) -> list[list[tuple[int, int]]] | None:
    if value is None:
        return None
    return [atom_map(mapping) for mapping in value]


def common_best_arguments(call: dict[str, Any]) -> dict[str, Any]:
    return {
        "prbCid": int(call.get("probe_conformer_id", -1)),
        "refCid": int(call.get("reference_conformer_id", -1)),
        "map": atom_maps(call.get("atom_maps")),
        "maxMatches": int(call.get("max_matches", 1_000_000)),
        "symmetrizeConjugatedTerminalGroups": bool(
            call.get("symmetrize_conjugated_terminal_groups", True)
        ),
        "weights": [float(value) for value in call.get("weights", [])],
        "numThreads": int(call.get("num_threads", 1)),
    }


def classify_error(operation: str, error: Exception) -> str:
    message = str(error)
    if "Bad Conformer Id" in message:
        return "conformer_not_found"
    if (
        "Mismatch in number of weights" in message
        or "Incorrect number of weights specified" in message
    ):
        return "weight_count_mismatch"
    if "No sub-structure match found" in message:
        return "no_substructure_match"
    return f"unclassified:{operation}:{type(error).__name__}"


def execute_pair_call(
    operation: str,
    call: dict[str, Any],
    probe: Chem.Mol,
    reference: Chem.Mol,
) -> dict[str, Any]:
    if operation == "alignment_transform":
        mapping = atom_map(call.get("atom_map", []))
        rmsd, transform = rdMolAlign.GetAlignmentTransform(
            probe,
            reference,
            prbCid=int(call.get("probe_conformer_id", -1)),
            refCid=int(call.get("reference_conformer_id", -1)),
            atomMap=mapping,
            weights=[float(value) for value in call.get("weights", [])],
            reflect=bool(call.get("reflect", False)),
            maxIters=int(call.get("max_iterations", 50)),
        )
        return {
            "rmsd": rmsd,
            "transform": transform.tolist(),
            "atom_map": [list(pair) for pair in mapping],
        }
    if operation == "best_alignment":
        arguments = common_best_arguments(call)
        arguments.update(
            {
                "reflect": bool(call.get("reflect", False)),
                "maxIters": int(call.get("max_iterations", 50)),
            }
        )
        rmsd, transform, mapping = rdMolAlign.GetBestAlignmentTransform(
            probe, reference, **arguments
        )
        return {
            "rmsd": rmsd,
            "transform": transform.tolist(),
            "atom_map": [list(pair) for pair in mapping],
        }
    if operation == "coordinate_rmsd":
        arguments = common_best_arguments(call)
        arguments.pop("numThreads")
        arguments["prbId"] = arguments.pop("prbCid")
        arguments["refId"] = arguments.pop("refCid")
        return {"rmsd": rdMolAlign.CalcRMS(probe, reference, **arguments)}
    if operation == "align_to":
        mapping = atom_map(call.get("atom_map", []))
        arguments = {
            "prbCid": int(call.get("probe_conformer_id", -1)),
            "refCid": int(call.get("reference_conformer_id", -1)),
            "atomMap": mapping,
            "weights": [float(value) for value in call.get("weights", [])],
            "reflect": bool(call.get("reflect", False)),
            "maxIters": int(call.get("max_iterations", 50)),
        }
        expected_rmsd, transform = rdMolAlign.GetAlignmentTransform(
            probe, reference, **arguments
        )
        rmsd = rdMolAlign.AlignMol(probe, reference, **arguments)
        if rmsd != expected_rmsd:
            raise AssertionError(
                "AlignMol and GetAlignmentTransform returned different RMSD values"
            )
        return {
            "rmsd": rmsd,
            "transform": transform.tolist(),
            "atom_map": [list(pair) for pair in mapping],
        }
    raise ValueError(f"unsupported pair operation: {operation}")


def execute_single_call(
    operation: str, call: dict[str, Any], molecule: Chem.Mol
) -> dict[str, Any]:
    if operation == "all_conformer_best_rms":
        maps = atom_maps(call.get("atom_maps"))
        values = rdMolAlign.GetAllConformerBestRMS(
            molecule,
            numThreads=int(call.get("num_threads", 1)),
            map=maps,
            maxMatches=int(call.get("max_matches", 1_000_000)),
            symmetrizeConjugatedTerminalGroups=bool(
                call.get("symmetrize_conjugated_terminal_groups", True)
            ),
            weights=[float(value) for value in call.get("weights", [])],
        )
        ids = [conformer.GetId() for conformer in molecule.GetConformers()]
        pairs = [
            [ids[probe_index], ids[reference_index]]
            for probe_index in range(len(ids))
            for reference_index in range(probe_index)
        ]
        return {"rmsds": list(values), "conformer_pairs": pairs}
    if operation == "align_conformers":
        rmsds: list[float] = []
        rdMolAlign.AlignMolConformers(
            molecule,
            atomIds=[int(value) for value in call.get("atom_indices", [])],
            confIds=[int(value) for value in call.get("conformer_ids", [])],
            weights=[float(value) for value in call.get("weights", [])],
            reflect=bool(call.get("reflect", False)),
            maxIters=int(call.get("max_iterations", 50)),
            RMSlist=rmsds,
        )
        return {"rmsds": rmsds}
    raise ValueError(f"unsupported single-molecule operation: {operation}")


def record_for_call(case: dict[str, Any], call_index: int, call: dict[str, Any]) -> dict[str, Any]:
    operation = call["operation"]
    if operation == "input_parse":
        return {
            "schema_version": 1,
            "case_id": case["case_id"],
            "call_index": call_index,
            "operation": operation,
            "source": {"input_smiles": case["input_smiles"]},
            "parameters": copy.deepcopy(call),
            "status": "error",
            "result": None,
            "error_kind": "input_parse_error",
            "error_type": "MolFromSmiles",
            "error_message": "RDKit returned no molecule",
            "before": {"probe": None, "reference": None, "molecule": None},
            "after": {"probe": None, "reference": None, "molecule": None},
        }
    is_pair = "probe" in case
    source = (
        {"probe": case["probe"], "reference": case["reference"]}
        if is_pair
        else {"molecule": case["molecule"]}
    )
    probe = molecule_from_source(case["probe"]) if is_pair else None
    reference = molecule_from_source(case["reference"]) if is_pair else None
    molecule = molecule_from_source(case["molecule"]) if not is_pair else None
    before = {
        "probe": snapshot(probe) if probe is not None else None,
        "reference": snapshot(reference) if reference is not None else None,
        "molecule": snapshot(molecule) if molecule is not None else None,
    }
    try:
        result = (
            execute_pair_call(operation, call, probe, reference)
            if is_pair
            else execute_single_call(operation, call, molecule)
        )
        status = "ok"
        error_kind = None
        error_type = None
        error_message = None
    except Exception as error:  # noqa: BLE001 - exception state is oracle output
        result = None
        status = "error"
        error_kind = classify_error(operation, error)
        error_type = type(error).__name__
        error_message = str(error)
    after = {
        "probe": snapshot(probe) if probe is not None else None,
        "reference": snapshot(reference) if reference is not None else None,
        "molecule": snapshot(molecule) if molecule is not None else None,
    }
    return {
        "schema_version": 1,
        "case_id": case["case_id"],
        "call_index": call_index,
        "operation": operation,
        "source": source,
        "parameters": copy.deepcopy(call),
        "status": status,
        "result": result,
        "error_kind": error_kind,
        "error_type": error_type,
        "error_message": error_message,
        "before": before,
        "after": after,
    }


def corpus_coordinates(row_index: int, atom_count: int) -> list[list[float]]:
    phase = (row_index % 97) * 0.013
    return [
        [
            atom_index * 0.37 + 0.05 * math.sin(atom_index * 0.19 + phase),
            math.sin(atom_index * 0.71 + phase),
            math.cos(atom_index * 0.43 - phase) + 0.03 * (atom_index % 5),
        ]
        for atom_index in range(atom_count)
    ]


def transform_coordinates(
    coordinates: list[list[float]],
    row_index: int,
    *,
    variant: int,
) -> list[list[float]]:
    angle = 0.19 + 0.017 * ((row_index + variant * 11) % 23)
    cosine = math.cos(angle)
    sine = math.sin(angle)
    translation = [
        0.7 + 0.03 * (row_index % 17),
        -0.4 + 0.02 * (row_index % 13),
        0.2 + 0.01 * (variant % 7),
    ]
    transformed = []
    for atom_index, (x, y, z) in enumerate(coordinates):
        perturbation = 0.0025 * (((atom_index + row_index + variant) % 7) - 3)
        transformed.append(
            [
                cosine * x - sine * y + translation[0] + perturbation,
                sine * x + cosine * y + translation[1] - perturbation * 0.5,
                z + translation[2] + perturbation * 0.25,
            ]
        )
    return transformed


def corpus_case(row_index: int, smiles: str) -> dict[str, Any]:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return {
            "case_id": f"smiles-{row_index + 1:05d}",
            "input_smiles": smiles,
            "calls": [{"operation": "input_parse"}],
        }
    atom_count = molecule.GetNumAtoms()
    identity_map = [[atom_index, atom_index] for atom_index in range(atom_count)]
    coordinates = corpus_coordinates(row_index, atom_count)
    weights = [1.0 + 0.25 * ((atom_index + row_index) % 4) for atom_index in range(atom_count)]
    operation = row_index % 6
    case_id = f"smiles-{row_index + 1:05d}"
    if operation < 4:
        probe = transform_coordinates(coordinates, row_index, variant=1)
        call: dict[str, Any]
        if operation == 0:
            call = {
                "operation": "alignment_transform",
                "probe_conformer_id": 7,
                "reference_conformer_id": 17,
                "atom_map": identity_map,
                "weights": weights,
                "reflect": row_index % 12 == 0,
                "max_iterations": 50,
            }
        elif operation == 1:
            alternate_map = list(reversed(identity_map))
            call = {
                "operation": "best_alignment",
                "probe_conformer_id": 7,
                "reference_conformer_id": 17,
                "atom_maps": [alternate_map, identity_map],
                "weights": weights,
                "reflect": row_index % 12 == 1,
                "num_threads": 2,
            }
        elif operation == 2:
            call = {
                "operation": "coordinate_rmsd",
                "probe_conformer_id": 7,
                "reference_conformer_id": 17,
                "atom_maps": [identity_map],
                "weights": weights,
            }
        else:
            call = {
                "operation": "align_to",
                "probe_conformer_id": 7,
                "reference_conformer_id": 17,
                "atom_map": identity_map,
                "weights": weights,
                "reflect": row_index % 12 == 3,
                "max_iterations": 50,
            }
        return {
            "case_id": case_id,
            "probe": {
                "smiles": smiles,
                "conformers": [{"id": 7, "coordinates": probe}],
            },
            "reference": {
                "smiles": smiles,
                "conformers": [{"id": 17, "coordinates": coordinates}],
            },
            "calls": [call],
        }

    conformers = [
        {"id": 7, "coordinates": coordinates},
        {
            "id": 27,
            "coordinates": transform_coordinates(coordinates, row_index, variant=2),
        },
        {
            "id": 17,
            "coordinates": transform_coordinates(coordinates, row_index, variant=3),
        },
    ]
    if operation == 4:
        call = {
            "operation": "all_conformer_best_rms",
            "atom_maps": [identity_map],
            "weights": weights,
            "num_threads": 2,
        }
    else:
        call = {
            "operation": "align_conformers",
            "atom_indices": list(range(atom_count)),
            "conformer_ids": [27, 17, 7],
            "weights": weights,
            "reflect": row_index % 12 == 5,
            "max_iterations": 50,
        }
    return {
        "case_id": case_id,
        "molecule": {"smiles": smiles, "conformers": conformers},
        "calls": [call],
    }


def load_cases(input_path: Path) -> list[dict[str, Any]]:
    if input_path.suffix.lower() == ".json":
        fixture = json.loads(input_path.read_text(encoding="utf-8"))
        if fixture.get("schema_version") != 1:
            raise SystemExit("unsupported focused MolAlign fixture schema")
        return fixture["cases"]
    return [
        corpus_case(row_index, line.strip())
        for row_index, line in enumerate(input_path.read_text(encoding="utf-8").splitlines())
        if line.strip()
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert_rdkit_version()
    cases = load_cases(args.input)
    records = [
        record_for_call(case, call_index, call)
        for case in cases
        for call_index, call in enumerate(case["calls"])
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True, separators=(",", ":")))
            handle.write("\n")


if __name__ == "__main__":
    main()
