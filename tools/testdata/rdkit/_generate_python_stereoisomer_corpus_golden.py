#!/usr/bin/env python3
"""Generate compact pinned-RDKit stereoisomer corpus oracle records."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
from typing import Any, Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    GetStereoisomerCount,
    StereoEnumerationOptions,
)


EXPECTED_RDKIT_VERSION = "2026.3.1"
NO_ATOM = int(Chem.StereoInfo().NOATOM)

POTENTIAL_PROFILES = (
    {"id": "preserve_possible", "clean_it": False, "flag_possible": True},
    {"id": "clean_possible", "clean_it": True, "flag_possible": True},
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


def iter_smiles(path: Path) -> Iterable[tuple[int, str]]:
    row = 0
    for raw in path.read_text(encoding="utf-8").splitlines():
        smiles = raw.strip()
        if not smiles or smiles.startswith("#"):
            continue
        yield row, smiles
        row += 1


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
        "STEREOANY": "Any",
        "STEREOZ": "Z",
        "STEREOE": "E",
        "STEREOCIS": "Cis",
        "STEREOTRANS": "Trans",
        "STEREOATROPCW": "AtropCw",
        "STEREOATROPCCW": "AtropCcw",
    }[str(value)]


def _enumeration_direction_gauge(molecule: Chem.Mol) -> dict[int, str]:
    """Return stable directions modulo RDKit's component-wide slash gauge."""
    directed = {
        bond.GetIdx(): bond
        for bond in molecule.GetBonds()
        if str(bond.GetBondDir()) in {"ENDDOWNRIGHT", "ENDUPRIGHT"}
    }
    participating: set[int] = set()
    parent = {index: index for index in directed}

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

    for double_bond in molecule.GetBonds():
        if str(double_bond.GetStereo()) not in {
            "STEREOCIS",
            "STEREOTRANS",
            "STEREOZ",
            "STEREOE",
        }:
            continue
        component = []
        for atom in (double_bond.GetBeginAtom(), double_bond.GetEndAtom()):
            component.extend(
                bond.GetIdx()
                for bond in atom.GetBonds()
                if bond.GetIdx() != double_bond.GetIdx()
                and bond.GetIdx() in directed
            )
        participating.update(component)
        for other in component[1:]:
            union(component[0], other)

    components: dict[int, list[int]] = {}
    for index in participating:
        components.setdefault(find(index), []).append(index)
    normalized = {
        index: normalize_bond_direction(bond.GetBondDir())
        for index, bond in directed.items()
        if index in participating
    }
    for component in components.values():
        first = min(component)
        if normalized[first] == "EndDownRight":
            for index in component:
                normalized[index] = {
                    "EndDownRight": "EndUpRight",
                    "EndUpRight": "EndDownRight",
                }[normalized[index]]
    return normalized


def molecule_stereo_state(
    molecule: Chem.Mol, *, normalize_enumeration_directions: bool = False
) -> dict[str, Any]:
    writer_copy = Chem.Mol(molecule)
    direction_gauge = (
        _enumeration_direction_gauge(molecule)
        if normalize_enumeration_directions
        else {}
    )
    return {
        "canonical_smiles": Chem.MolToSmiles(writer_copy, isomericSmiles=True),
        "atom_chiral_tags": [
            normalize_chiral_tag(atom.GetChiralTag()) for atom in molecule.GetAtoms()
        ],
        "bonds": [
            {
                "direction": direction_gauge.get(
                    bond.GetIdx(), normalize_bond_direction(bond.GetBondDir())
                ),
                "stereo": normalize_bond_stereo(bond.GetStereo()),
                "stereo_atoms": [int(atom) for atom in bond.GetStereoAtoms()],
            }
            for bond in molecule.GetBonds()
        ],
        "conformer_count": molecule.GetNumConformers(),
    }


def stereo_info_state(info: Any) -> dict[str, Any]:
    stereo_type = {
        "Unspecified": "unspecified",
        "Atom_Tetrahedral": "atom_tetrahedral",
        "Atom_SquarePlanar": "atom_square_planar",
        "Atom_TrigonalBipyramidal": "atom_trigonal_bipyramidal",
        "Atom_Octahedral": "atom_octahedral",
        "Bond_Double": "bond_double",
        "Bond_Cumulene_Even": "bond_even_cumulene",
        "Bond_Atropisomer": "bond_atropisomer",
    }[str(info.type)]
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


def potential_run(source: Chem.Mol, profile: dict[str, Any]) -> dict[str, Any]:
    working = Chem.Mol(source)
    infos = Chem.FindPotentialStereo(
        working,
        cleanIt=profile["clean_it"],
        flagPossible=profile["flag_possible"],
    )
    return {
        "profile": profile,
        "records": [stereo_info_state(info) for info in infos],
        "state": molecule_stereo_state(working),
    }


def enumeration_options(profile: dict[str, Any]) -> StereoEnumerationOptions:
    return StereoEnumerationOptions(
        tryEmbedding=False,
        onlyUnassigned=profile["only_unassigned"],
        onlyStereoGroups=profile["only_stereo_groups"],
        maxIsomers=profile["max_isomers"],
        rand=profile["rand"],
        unique=profile["unique"],
    )


def enumeration_run(source: Chem.Mol, profile: dict[str, Any]) -> dict[str, Any]:
    options = enumeration_options(profile)
    theoretical_count = GetStereoisomerCount(source, options)
    bounded_out = profile["id"] != "seeded_three" and theoretical_count > 8
    outputs = (
        []
        if bounded_out
        else [
            molecule_stereo_state(isomer, normalize_enumeration_directions=True)
            for isomer in EnumerateStereoisomers(source, options)
        ]
    )
    return {
        "profile": profile,
        "theoretical_count": str(theoretical_count),
        "bounded_out": bounded_out,
        "outputs": outputs,
    }


def build_record(row: int, smiles: str) -> dict[str, Any]:
    try:
        source = Chem.MolFromSmiles(smiles, sanitize=True)
        if source is None:
            return {
                "schema_version": 2,
                "row": row,
                "smiles": smiles,
                "parse_status": "none",
                "error_type": None,
                "error_message": None,
                "source_state": None,
                "potential_stereo": [],
                "enumeration": [],
            }
        source_state = molecule_stereo_state(source)
        potential = [potential_run(source, profile) for profile in POTENTIAL_PROFILES]
        enumeration = [
            enumeration_run(source, profile) for profile in ENUMERATION_PROFILES
        ]
        if molecule_stereo_state(source) != source_state:
            raise RuntimeError("reference operations mutated the source molecule")
        return {
            "schema_version": 2,
            "row": row,
            "smiles": smiles,
            "parse_status": "ok",
            "error_type": None,
            "error_message": None,
            "source_state": source_state,
            "potential_stereo": potential,
            "enumeration": enumeration,
        }
    except Exception as error:  # noqa: BLE001 - reference failure is oracle data
        return {
            "schema_version": 2,
            "row": row,
            "smiles": smiles,
            "parse_status": "error",
            "error_type": type(error).__name__,
            "error_message": str(error),
            "source_state": None,
            "potential_stereo": [],
            "enumeration": [],
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    actual_version = version("rdkit")
    if actual_version != EXPECTED_RDKIT_VERSION:
        raise SystemExit(
            f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, got {actual_version}"
        )
    RDLogger.DisableLog("rdApp.*")
    rows = list(iter_smiles(args.input))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as output:
        for row, smiles in rows:
            output.write(
                json.dumps(
                    build_record(row, smiles),
                    sort_keys=True,
                    separators=(",", ":"),
                )
            )
            output.write("\n")
    print(f"Wrote {len(rows)} records to {args.output}")


if __name__ == "__main__":
    main()
