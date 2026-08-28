#!/usr/bin/env python3
"""Differential ChEMBL audit for potential stereo and stereoisomer enumeration."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Callable

import cosmolkit
from rdkit import Chem, RDLogger
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    GetStereoisomerCount,
    StereoEnumerationOptions,
)

from .audit_core import Audit, atomic_json, validate_runtime_versions


RDLogger.DisableLog("rdApp.*")

POTENTIAL_PROFILES = (
    {"id": "preserve_possible", "clean": False, "flag_possible": True},
    {"id": "clean_possible", "clean": True, "flag_possible": True},
    {"id": "preserve_known", "clean": False, "flag_possible": False},
    {"id": "clean_known", "clean": True, "flag_possible": False},
)

ENUMERATION_PROFILES = (
    {
        "id": "default_bounded",
        "only_unassigned": True,
        "only_stereo_groups": False,
        "max_isomers": 8,
        "rand": None,
        "unique": True,
    },
    {
        "id": "all_assigned_bounded",
        "only_unassigned": False,
        "only_stereo_groups": False,
        "max_isomers": 8,
        "rand": None,
        "unique": True,
    },
    {
        "id": "non_unique_bounded",
        "only_unassigned": True,
        "only_stereo_groups": False,
        "max_isomers": 8,
        "rand": None,
        "unique": False,
    },
    {
        "id": "seeded_three",
        "only_unassigned": True,
        "only_stereo_groups": False,
        "max_isomers": 3,
        "rand": 61453,
        "unique": True,
    },
)

EXPECTED_PROFILES = {"potential_stereo": 4, "enumeration": 4}
NO_ATOM = int(Chem.StereoInfo().NOATOM)


def normalize_chiral_tag(value: Any) -> str:
    return {
        "CHI_UNSPECIFIED": "Unspecified",
        "CHI_TETRAHEDRAL_CW": "TetrahedralCw",
        "CHI_TETRAHEDRAL_CCW": "TetrahedralCcw",
        "CHI_OTHER": "Other",
        "CHI_TETRAHEDRAL": "Tetrahedral",
        "CHI_ALLENE": "Allene",
        "CHI_SQUAREPLANAR": "SquarePlanar",
        "CHI_TRIGONALBIPYRAMIDAL": "TrigonalBipyramidal",
        "CHI_OCTAHEDRAL": "Octahedral",
    }[str(value)]


def normalize_bond_direction(value: Any) -> str:
    return {
        "NONE": "None",
        "BEGINWEDGE": "BeginWedge",
        "BEGINDASH": "BeginDash",
        "ENDDOWNRIGHT": "EndDownRight",
        "ENDUPRIGHT": "EndUpRight",
        "EITHERDOUBLE": "EitherDouble",
        "UNKNOWN": "Unknown",
    }[str(value)]


def normalize_bond_stereo(value: Any) -> str:
    return {
        "STEREONONE": "None",
        "NONE": "None",
        "STEREOANY": "Any",
        "ANY": "Any",
        "STEREOZ": "Z",
        "Z": "Z",
        "STEREOE": "E",
        "E": "E",
        "STEREOCIS": "Cis",
        "CIS": "Cis",
        "STEREOTRANS": "Trans",
        "TRANS": "Trans",
        "STEREOATROPCW": "AtropCw",
        "ATROPCW": "AtropCw",
        "STEREOATROPCCW": "AtropCcw",
        "ATROPCCW": "AtropCcw",
    }[str(value)]


def normalize_stereo_type(value: Any) -> str:
    return {
        "Unspecified": "unspecified",
        "Atom_Tetrahedral": "atom_tetrahedral",
        "Atom_SquarePlanar": "atom_square_planar",
        "Atom_TrigonalBipyramidal": "atom_trigonal_bipyramidal",
        "Atom_Octahedral": "atom_octahedral",
        "Bond_Double": "bond_double",
        "Bond_Cumulene_Even": "bond_even_cumulene",
        "Bond_Atropisomer": "bond_atropisomer",
    }[str(value)]


def rd_stereo_info(info: Any) -> dict[str, Any]:
    stereo_type = normalize_stereo_type(info.type)
    center = int(info.centeredOn)
    return {
        "stereo_type": stereo_type,
        "specified": str(info.specified).lower(),
        "center_kind": "atom" if stereo_type.startswith("atom_") else "bond",
        "center_index": None if center == NO_ATOM else center,
        "descriptor": {
            "NoValue": "none",
            "Tet_CW": "tetrahedral_clockwise",
            "Tet_CCW": "tetrahedral_counterclockwise",
            "Bond_Cis": "bond_cis",
            "Bond_Trans": "bond_trans",
            "Bond_AtropCW": "bond_atrop_clockwise",
            "Bond_AtropCCW": "bond_atrop_counterclockwise",
        }[str(info.descriptor)],
        "permutation": int(info.permutation),
        "controlling_atoms": [
            None if int(atom) == NO_ATOM else int(atom)
            for atom in info.controllingAtoms
        ],
    }


def ck_stereo_info(info: Any) -> dict[str, Any]:
    return {
        "stereo_type": info.stereo_type,
        "specified": info.specified,
        "center_kind": info.center_kind,
        "center_index": info.center_index,
        "descriptor": info.descriptor,
        "permutation": info.permutation,
        "controlling_atoms": info.controlling_atoms,
    }


def normalized_direction_gauge(bonds: list[dict[str, Any]]) -> dict[int, str]:
    directed = {
        bond["index"]: bond
        for bond in bonds
        if bond["direction"] in {"EndDownRight", "EndUpRight"}
    }
    parent = {index: index for index in directed}
    participating: set[int] = set()

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for double_bond in bonds:
        if double_bond["stereo"] not in {"Cis", "Trans", "Z", "E"}:
            continue
        component = [
            bond["index"]
            for bond in bonds
            if bond["index"] in directed
            and (
                bond["begin"] in {double_bond["begin"], double_bond["end"]}
                or bond["end"] in {double_bond["begin"], double_bond["end"]}
            )
        ]
        if not component:
            continue
        participating.update(component)
        for other in component[1:]:
            union(component[0], other)

    normalized = {
        index: directed[index]["direction"] for index in participating
    }
    components: dict[int, list[int]] = {}
    for index in participating:
        components.setdefault(find(index), []).append(index)
    for component in components.values():
        if normalized[min(component)] == "EndDownRight":
            for index in component:
                normalized[index] = {
                    "EndDownRight": "EndUpRight",
                    "EndUpRight": "EndDownRight",
                }[normalized[index]]
    return normalized


def rd_bonds(molecule: Chem.Mol) -> list[dict[str, Any]]:
    return [
        {
            "index": bond.GetIdx(),
            "begin": bond.GetBeginAtomIdx(),
            "end": bond.GetEndAtomIdx(),
            "direction": normalize_bond_direction(bond.GetBondDir()),
            "stereo": normalize_bond_stereo(bond.GetStereo()),
            "stereo_atoms": [int(atom) for atom in bond.GetStereoAtoms()],
        }
        for bond in molecule.GetBonds()
    ]


def ck_bonds(molecule: Any) -> list[dict[str, Any]]:
    return [
        {
            "index": bond.idx(),
            "begin": bond.begin_atom_idx(),
            "end": bond.end_atom_idx(),
            "direction": normalize_bond_direction(bond.bond_dir_name()),
            "stereo": normalize_bond_stereo(bond.stereo_name()),
            "stereo_atoms": bond.stereo_atoms(),
        }
        for bond in molecule.bonds()
    ]


def molecule_state(
    canonical_smiles: str,
    atom_tags: list[str],
    bonds: list[dict[str, Any]],
    conformer_count: int,
    *,
    normalize_directions: bool,
) -> dict[str, Any]:
    gauge = normalized_direction_gauge(bonds) if normalize_directions else {}
    return {
        "canonical_smiles": canonical_smiles,
        "atom_chiral_tags": atom_tags,
        "bonds": [
            {
                "direction": gauge.get(bond["index"], bond["direction"]),
                "stereo": bond["stereo"],
                "stereo_atoms": bond["stereo_atoms"],
            }
            for bond in bonds
        ],
        "conformer_count": conformer_count,
    }


def rd_molecule_state(
    molecule: Chem.Mol, *, normalize_directions: bool = False
) -> dict[str, Any]:
    return molecule_state(
        Chem.MolToSmiles(molecule, isomericSmiles=True),
        [normalize_chiral_tag(atom.GetChiralTag()) for atom in molecule.GetAtoms()],
        rd_bonds(molecule),
        molecule.GetNumConformers(),
        normalize_directions=normalize_directions,
    )


def ck_molecule_state(
    molecule: Any, *, normalize_directions: bool = False
) -> dict[str, Any]:
    return molecule_state(
        molecule.to_smiles(),
        [normalize_chiral_tag(atom.chiral_tag_name()) for atom in molecule.atoms()],
        ck_bonds(molecule),
        molecule.num_conformers(),
        normalize_directions=normalize_directions,
    )


def capture(function: Callable[[], Any]) -> dict[str, Any]:
    try:
        return {"status": "ok", "value": function()}
    except Exception as error:  # noqa: BLE001 - preserve source error boundary
        return {"status": "error", "type": type(error).__name__}


def potential_rd(source: Chem.Mol, profile: dict[str, Any]) -> dict[str, Any]:
    def run() -> dict[str, Any]:
        working = Chem.Mol(source)
        records = Chem.FindPotentialStereo(
            working,
            cleanIt=profile["clean"],
            flagPossible=profile["flag_possible"],
        )
        return {
            "records": [rd_stereo_info(info) for info in records],
            "state": rd_molecule_state(working),
        }

    return capture(run)


def potential_ck(source: Any, profile: dict[str, Any]) -> dict[str, Any]:
    def run() -> dict[str, Any]:
        analysis = source.analyze_potential_stereo(
            clean=profile["clean"], flag_possible=profile["flag_possible"]
        )
        return {
            "records": [ck_stereo_info(info) for info in analysis.stereo_info],
            "state": ck_molecule_state(analysis.molecule),
        }

    return capture(run)


def rd_options(profile: dict[str, Any]) -> StereoEnumerationOptions:
    return StereoEnumerationOptions(
        tryEmbedding=False,
        onlyUnassigned=profile["only_unassigned"],
        onlyStereoGroups=profile["only_stereo_groups"],
        maxIsomers=profile["max_isomers"],
        rand=profile["rand"],
        unique=profile["unique"],
    )


def ck_options(profile: dict[str, Any]) -> Any:
    return cosmolkit.StereoisomerOptions(
        try_embedding=False,
        only_unassigned=profile["only_unassigned"],
        only_stereo_groups=profile["only_stereo_groups"],
        max_isomers=profile["max_isomers"],
        rand=profile["rand"],
        unique=profile["unique"],
    )


def enumeration_rd(source: Chem.Mol, profile: dict[str, Any]) -> dict[str, Any]:
    def run() -> dict[str, Any]:
        options = rd_options(profile)
        count = int(GetStereoisomerCount(source, options))
        bounded_out = profile["id"] != "seeded_three" and count > 8
        outputs = (
            []
            if bounded_out
            else [
                rd_molecule_state(isomer, normalize_directions=True)
                for isomer in EnumerateStereoisomers(source, options)
            ]
        )
        return {
            "theoretical_count": str(count),
            "bounded_out": bounded_out,
            "outputs": outputs,
        }

    return capture(run)


def enumeration_ck(source: Any, profile: dict[str, Any]) -> dict[str, Any]:
    def run() -> dict[str, Any]:
        options = ck_options(profile)
        count = int(source.stereoisomer_count(options))
        bounded_out = profile["id"] != "seeded_three" and count > 8
        outputs = (
            []
            if bounded_out
            else [
                ck_molecule_state(isomer, normalize_directions=True)
                for isomer in source.stereoisomers(options)
            ]
        )
        return {
            "theoretical_count": str(count),
            "bounded_out": bounded_out,
            "outputs": outputs,
        }

    return capture(run)


def compare(
    audit: Audit,
    feature: str,
    record: dict[str, Any],
    rd_value: Any,
    ck_value: Any,
) -> None:
    if rd_value == ck_value:
        audit.count(f"match.{feature}")
    else:
        audit.mismatch(feature, record, rd_value, ck_value)


def compare_capture(
    audit: Audit,
    feature: str,
    record: dict[str, Any],
    rd_value: dict[str, Any],
    ck_value: dict[str, Any],
) -> None:
    compare(audit, f"{feature}.status", record, rd_value["status"], ck_value["status"])
    if rd_value["status"] != "ok" or ck_value["status"] != "ok":
        compare(audit, f"{feature}.error_type", record, rd_value.get("type"), ck_value.get("type"))
        return
    rd_inner = rd_value["value"]
    ck_inner = ck_value["value"]
    fields = (
        ("records", "state")
        if "records" in rd_inner
        else ("theoretical_count", "bounded_out", "outputs")
    )
    for field in fields:
        compare(audit, f"{feature}.{field}", record, rd_inner[field], ck_inner[field])


def audit_record(audit: Audit, record: dict[str, Any]) -> None:
    rd_parse = capture(lambda: Chem.MolFromSmiles(record["smiles"], sanitize=True))
    ck_parse = capture(lambda: cosmolkit.Molecule.from_smiles(record["smiles"]))
    rd_source = rd_parse.get("value") if rd_parse["status"] == "ok" else None
    ck_source = ck_parse.get("value") if ck_parse["status"] == "ok" else None
    rd_accepted = rd_source is not None
    ck_accepted = ck_source is not None
    compare(audit, "parse.accepted", record, rd_accepted, ck_accepted)
    if not rd_accepted or not ck_accepted:
        if not rd_accepted and not ck_accepted:
            audit.count("parse.both_rejected")
        return
    rd_before = rd_molecule_state(rd_source)
    ck_before = ck_molecule_state(ck_source)
    for profile in POTENTIAL_PROFILES:
        compare_capture(
            audit,
            f"potential.{profile['id']}",
            record,
            potential_rd(rd_source, profile),
            potential_ck(ck_source, profile),
        )
    for profile in ENUMERATION_PROFILES:
        compare_capture(
            audit,
            f"enumeration.{profile['id']}",
            record,
            enumeration_rd(rd_source, profile),
            enumeration_ck(ck_source, profile),
        )
    compare(audit, "source.rdkit_unchanged", record, rd_before, rd_molecule_state(rd_source))
    compare(audit, "source.cosmolkit_unchanged", record, ck_before, ck_molecule_state(ck_source))


def main() -> None:
    validate_runtime_versions()
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mode", choices=("stereo",), required=True)
    parser.add_argument("--max-examples", type=int, default=12)
    args = parser.parse_args()

    audit = Audit(args.output, args.max_examples)
    processed = 0
    try:
        with args.input.open(encoding="utf-8") as source:
            for processed, line in enumerate(source, start=1):
                audit_record(audit, json.loads(line))
                if processed % 1000 == 0:
                    print(json.dumps({"mode": "stereo", "processed": processed}), flush=True)
    finally:
        audit.finish(processed, "stereo", 0, 1)
        summary = json.loads(args.output.read_text(encoding="utf-8"))
        summary["profiles"] = EXPECTED_PROFILES
        atomic_json(args.output, summary)


if __name__ == "__main__":
    main()
