#!/usr/bin/env python3
"""Generate pinned-RDKit Python stereoisomer-enumeration oracle records."""

from __future__ import annotations

import argparse
import json
import random
from importlib.metadata import version
from pathlib import Path
from typing import Any

from rdkit import Chem, RDLogger
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    GetStereoisomerCount,
    StereoEnumerationOptions,
)


EXPECTED_RDKIT_VERSION = "2026.3.1"
REPO_ROOT = Path(__file__).resolve().parents[3]


def assert_rdkit_version() -> None:
    actual = version("rdkit")
    if actual != EXPECTED_RDKIT_VERSION:
        raise SystemExit(
            f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, got {actual}"
        )


def json_value(value: Any) -> Any:
    if isinstance(value, float):
        return value.hex()
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    if isinstance(value, (tuple, list)):
        return [json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): json_value(item) for key, item in sorted(value.items())}
    if hasattr(value, "__iter__"):
        return [json_value(item) for item in value]
    return str(value)


def property_snapshot(owner: Any) -> dict[str, Any]:
    persistent_names = sorted(
        owner.GetPropNames(includePrivate=True, includeComputed=False)
    )
    complete_names = sorted(
        owner.GetPropNames(includePrivate=True, includeComputed=True)
    )
    values = owner.GetPropsAsDict(includePrivate=True, includeComputed=True)
    return {
        "persistent_names": persistent_names,
        "computed_names": sorted(set(complete_names) - set(persistent_names)),
        "values": json_value(values),
    }


def molecule_snapshot(mol: Chem.Mol) -> dict[str, Any]:
    atoms = []
    for atom in mol.GetAtoms():
        atoms.append(
            {
                "index": atom.GetIdx(),
                "atomic_number": atom.GetAtomicNum(),
                "isotope": atom.GetIsotope(),
                "formal_charge": atom.GetFormalCharge(),
                "explicit_hydrogens": atom.GetNumExplicitHs(),
                "implicit_hydrogens": atom.GetNumImplicitHs(),
                "no_implicit": atom.GetNoImplicit(),
                "radical_electrons": atom.GetNumRadicalElectrons(),
                "aromatic": atom.GetIsAromatic(),
                "hybridization": str(atom.GetHybridization()),
                "chiral_tag": str(atom.GetChiralTag()),
                "properties": property_snapshot(atom),
            }
        )
    bonds = []
    for bond in mol.GetBonds():
        bonds.append(
            {
                "index": bond.GetIdx(),
                "begin": bond.GetBeginAtomIdx(),
                "end": bond.GetEndAtomIdx(),
                "bond_type": str(bond.GetBondType()),
                "aromatic": bond.GetIsAromatic(),
                "conjugated": bond.GetIsConjugated(),
                "direction": str(bond.GetBondDir()),
                "stereo": str(bond.GetStereo()),
                "stereo_atoms": [int(item) for item in bond.GetStereoAtoms()],
                "properties": property_snapshot(bond),
            }
        )
    stereo_groups = []
    for group in mol.GetStereoGroups():
        stereo_groups.append(
            {
                "type": str(group.GetGroupType()),
                "read_id": group.GetReadId(),
                "write_id": group.GetWriteId(),
                "atoms": [atom.GetIdx() for atom in group.GetAtoms()],
                "bonds": [bond.GetIdx() for bond in group.GetBonds()],
            }
        )
    conformers = []
    for conformer in mol.GetConformers():
        coordinates = []
        for atom_index in range(mol.GetNumAtoms()):
            point = conformer.GetAtomPosition(atom_index)
            coordinates.append([point.x.hex(), point.y.hex(), point.z.hex()])
        conformers.append(
            {
                "id": conformer.GetId(),
                "is_3d": conformer.Is3D(),
                "coordinates": coordinates,
                "properties": property_snapshot(conformer),
            }
        )

    # Writers can populate computed properties, so serialize observable text from a
    # clone and leave the state captured above untouched.
    writer_copy = Chem.Mol(mol)
    try:
        canonical_smiles = Chem.MolToSmiles(writer_copy, isomericSmiles=True)
        canonical_smiles_error = None
    except Exception as error:  # noqa: BLE001 - exception identity is oracle data
        canonical_smiles = None
        canonical_smiles_error = {
            "type": type(error).__name__,
            "message": str(error),
        }
    return {
        "atom_count": mol.GetNumAtoms(),
        "bond_count": mol.GetNumBonds(),
        "atoms": atoms,
        "bonds": bonds,
        "stereo_groups": stereo_groups,
        "conformers": conformers,
        "properties": property_snapshot(mol),
        "canonical_isomeric_smiles": canonical_smiles,
        "canonical_isomeric_smiles_error": canonical_smiles_error,
    }


def materialize_input(input_spec: dict[str, Any]) -> tuple[str, str]:
    kind = input_spec["kind"]
    if kind in {"smiles", "cxsmiles"}:
        return kind, input_spec["value"]
    if kind == "repeat_smiles":
        value = (
            input_spec["prefix"]
            + input_spec["unit"] * int(input_spec["repeat"])
            + input_spec["suffix"]
        )
        return "smiles", value
    if kind == "mol_file":
        return kind, input_spec["path"]
    raise ValueError(f"unsupported input kind: {kind!r}")


def parse_input(input_spec: dict[str, Any]) -> Chem.Mol | None:
    kind, source = materialize_input(input_spec)
    if kind in {"smiles", "cxsmiles"}:
        params = Chem.SmilesParserParams()
        params.sanitize = True
        params.removeHs = True
        params.allowCXSMILES = kind == "cxsmiles"
        params.parseName = False
        return Chem.MolFromSmiles(source, params)
    if kind == "mol_file":
        return Chem.MolFromMolFile(
            str(REPO_ROOT / source),
            sanitize=True,
            removeHs=True,
            strictParsing=True,
        )
    raise AssertionError(f"unhandled input kind: {kind!r}")


CHIRAL_TAGS = {
    "TetrahedralCw": Chem.ChiralType.CHI_TETRAHEDRAL_CW,
    "TetrahedralCcw": Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
    "Unspecified": Chem.ChiralType.CHI_UNSPECIFIED,
}
BOND_STEREO = {
    "Any": Chem.BondStereo.STEREOANY,
    "None": Chem.BondStereo.STEREONONE,
    "Cis": Chem.BondStereo.STEREOCIS,
    "Trans": Chem.BondStereo.STEREOTRANS,
    "E": Chem.BondStereo.STEREOE,
    "Z": Chem.BondStereo.STEREOZ,
}


def apply_mutations(mol: Chem.Mol, mutations: list[dict[str, Any]]) -> None:
    for mutation in mutations:
        kind = mutation["kind"]
        if kind == "set_atom_chiral_tag":
            mol.GetAtomWithIdx(int(mutation["atom"])).SetChiralTag(
                CHIRAL_TAGS[mutation["value"]]
            )
        elif kind == "set_bond_stereo":
            mol.GetBondWithIdx(int(mutation["bond"])).SetStereo(
                BOND_STEREO[mutation["value"]]
            )
        else:
            raise ValueError(f"unsupported mutation kind: {kind!r}")


def stereo_info_snapshot(info: Any) -> dict[str, Any]:
    return {
        "centered_on": int(info.centeredOn),
        "type": str(info.type),
        "specified": str(info.specified),
        "descriptor": str(info.descriptor),
        "permutation": int(info.permutation),
        "controlling_atoms": [int(item) for item in info.controllingAtoms],
    }


def potential_stereo_run(
    source: Chem.Mol, profile: dict[str, Any]
) -> dict[str, Any]:
    working = Chem.Mol(source)
    before = molecule_snapshot(working)
    try:
        infos = Chem.FindPotentialStereo(
            working,
            cleanIt=bool(profile["clean_it"]),
            flagPossible=bool(profile["flag_possible"]),
        )
        result = {
            "status": "ok",
            "error_type": None,
            "error_message": None,
            "stereo_info": [stereo_info_snapshot(info) for info in infos],
        }
    except Exception as error:  # noqa: BLE001 - exception identity is oracle data
        result = {
            "status": "error",
            "error_type": type(error).__name__,
            "error_message": str(error),
            "stereo_info": [],
        }
    return {
        "options": {
            "clean_it": bool(profile["clean_it"]),
            "flag_possible": bool(profile["flag_possible"]),
        },
        "before": before,
        "result": result,
        "after": molecule_snapshot(working),
    }


class CounterRandom(random.Random):
    """Declared fixture random source: successive integers from getrandbits()."""

    def __init__(self) -> None:
        super().__init__(0)
        self.next_value = 0

    def getrandbits(self, bit_count: int) -> int:
        value = self.next_value
        self.next_value += 1
        if value >= 1 << bit_count:
            raise AssertionError("counter random source exhausted the declared bit width")
        return value


def enumeration_options(profile: dict[str, Any]) -> StereoEnumerationOptions:
    random_spec = profile["rand"]
    kind = random_spec["kind"]
    if kind == "default":
        random_source: Any = None
    elif kind == "integer":
        random_source = int(random_spec["value"])
    elif kind == "counter_getrandbits":
        random_source = CounterRandom()
    else:
        raise ValueError(f"unsupported random source kind: {kind!r}")
    return StereoEnumerationOptions(
        tryEmbedding=bool(profile["try_embedding"]),
        onlyUnassigned=bool(profile["only_unassigned"]),
        maxIsomers=int(profile["max_isomers"]),
        rand=random_source,
        unique=bool(profile["unique"]),
        onlyStereoGroups=bool(profile["only_stereo_groups"]),
    )


def profile_snapshot(profile: dict[str, Any]) -> dict[str, Any]:
    return {
        "id": profile["id"],
        "try_embedding": bool(profile["try_embedding"]),
        "only_unassigned": bool(profile["only_unassigned"]),
        "only_stereo_groups": bool(profile["only_stereo_groups"]),
        "max_isomers": int(profile["max_isomers"]),
        "rand": profile["rand"],
        "unique": bool(profile["unique"]),
    }


def enumeration_run(source: Chem.Mol, profile: dict[str, Any]) -> dict[str, Any]:
    before = molecule_snapshot(source)
    count_status = "ok"
    count_error_type = None
    count_error_message = None
    count = None
    try:
        count = GetStereoisomerCount(source, enumeration_options(profile))
    except Exception as error:  # noqa: BLE001 - exception identity is oracle data
        count_status = "error"
        count_error_type = type(error).__name__
        count_error_message = str(error)

    outputs = []
    enumeration_status = "ok"
    enumeration_error_type = None
    enumeration_error_message = None
    try:
        for output_index, isomer in enumerate(
            EnumerateStereoisomers(source, enumeration_options(profile))
        ):
            outputs.append(
                {
                    "output_index": output_index,
                    "state": molecule_snapshot(isomer),
                }
            )
    except Exception as error:  # noqa: BLE001 - exception identity is oracle data
        enumeration_status = "error"
        enumeration_error_type = type(error).__name__
        enumeration_error_message = str(error)
    return {
        "profile": profile_snapshot(profile),
        "source_before": before,
        "count": {
            "status": count_status,
            "value": count,
            "error_type": count_error_type,
            "error_message": count_error_message,
        },
        "enumeration": {
            "status": enumeration_status,
            "error_type": enumeration_error_type,
            "error_message": enumeration_error_message,
            "outputs": outputs,
        },
        "source_after": molecule_snapshot(source),
    }


def record_for_case(
    case: dict[str, Any], profiles: dict[str, dict[str, Any]]
) -> dict[str, Any]:
    input_kind, materialized_input = materialize_input(case["input"])
    source = {
        "kind": input_kind,
        "value": materialized_input,
        "declared": case["input"],
        "mutations": case.get("mutations", []),
        "upstream_sources": case["source"],
        "branches": case["branches"],
    }
    try:
        mol = parse_input(case["input"])
        if mol is None:
            return {
                "schema_version": 1,
                "case_id": case["id"],
                "source": source,
                "parse_status": "none",
                "parse_error_type": None,
                "parse_error_message": None,
                "initial_state": None,
                "potential_stereo_runs": [],
                "enumeration_runs": [],
                "final_source_state": None,
            }
        apply_mutations(mol, case.get("mutations", []))
    except Exception as error:  # noqa: BLE001 - exception identity is oracle data
        return {
            "schema_version": 1,
            "case_id": case["id"],
            "source": source,
            "parse_status": "error",
            "parse_error_type": type(error).__name__,
            "parse_error_message": str(error),
            "initial_state": None,
            "potential_stereo_runs": [],
            "enumeration_runs": [],
            "final_source_state": None,
        }

    initial = molecule_snapshot(mol)
    potential_runs = [
        potential_stereo_run(mol, profile)
        for profile in case.get("potential_stereo_profiles", [])
    ]
    enumeration_runs = [
        enumeration_run(mol, profiles[profile_id])
        for profile_id in case.get("profiles", [])
    ]
    return {
        "schema_version": 1,
        "case_id": case["id"],
        "source": source,
        "parse_status": "ok",
        "parse_error_type": None,
        "parse_error_message": None,
        "initial_state": initial,
        "potential_stereo_runs": potential_runs,
        "enumeration_runs": enumeration_runs,
        "final_source_state": molecule_snapshot(mol),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    document = json.loads(args.input.read_text(encoding="utf-8"))
    if document.get("schema_version") != 1:
        raise SystemExit("focused stereoisomer fixture must have schema_version 1")
    if not isinstance(document.get("option_profiles"), list) or not isinstance(
        document.get("cases"), list
    ):
        raise SystemExit("focused stereoisomer fixture requires option_profiles and cases")
    profiles = {profile["id"]: profile for profile in document["option_profiles"]}
    if len(profiles) != len(document["option_profiles"]):
        raise SystemExit("focused stereoisomer fixture profile IDs must be unique")
    records = [record_for_case(case, profiles) for case in document["cases"]]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True, separators=(",", ":")))
            handle.write("\n")
    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
