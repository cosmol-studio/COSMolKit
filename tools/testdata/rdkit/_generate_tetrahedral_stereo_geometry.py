#!/usr/bin/env python3
"""Generate RDKit ETKDG geometry golden data for tetrahedral stereo checks."""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
from importlib.metadata import version
from pathlib import Path
from typing import Any, Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

EXPECTED_RDKIT_VERSION = "2026.3.1"
ASSIGN_FIXTURE = Path(
    "testdata/stereo/fixtures/assign_atom_chiral_tags_from_structure_cases.json"
)
ASSIGN_OUTPUT_NAME = "assign_atom_chiral_tags_from_structure.jsonl"
NON_TETRAHEDRAL_ENV = "RDK_ENABLE_NONTETRAHEDRAL_STEREO"


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


def build_record(smiles: str) -> dict[str, object] | None:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return None

    chiral_centers = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if str(atom.GetChiralTag()) in {"CHI_TETRAHEDRAL_CW", "CHI_TETRAHEDRAL_CCW"}
    ]
    if not chiral_centers:
        return None

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    res = AllChem.EmbedMolecule(mol, params)
    if res != 0:
        return {
            "smiles": smiles,
            "rdkit_ok": False,
            "centers": chiral_centers,
            "positions": None,
            "error": "EmbedMolecule failed",
        }

    conf = mol.GetConformer()
    return {
        "smiles": smiles,
        "rdkit_ok": True,
        "centers": chiral_centers,
        "positions": conf.GetPositions().tolist(),
        "error": None,
    }


def deep_merge(base: dict[str, Any], update: dict[str, Any]) -> dict[str, Any]:
    result = dict(base)
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def decode_number(value: int | float | str) -> float:
    if isinstance(value, (int, float)):
        return float(value)
    tokens = {
        "negative_zero": -0.0,
        "min_subnormal": float.fromhex("0x0.0000000000001p-1022"),
        "max_finite": sys.float_info.max,
        "nextafter_0.1_down": math.nextafter(0.1, -math.inf),
        "nextafter_0.1_up": math.nextafter(0.1, math.inf),
        "nextafter_negative_0.1_down": math.nextafter(-0.1, -math.inf),
        "nextafter_negative_0.1_up": math.nextafter(-0.1, math.inf),
        "nextafter_zero_tolerance_down": math.nextafter(1.0e-16, -math.inf),
        "zero_tolerance": 1.0e-16,
        "nextafter_zero_tolerance_up": math.nextafter(1.0e-16, math.inf),
        "cos_100_degrees": math.cos(100.0 * math.pi / 180.0),
        "sin_100_degrees": math.sin(100.0 * math.pi / 180.0),
        "nan": math.nan,
        "positive_infinity": math.inf,
        "negative_infinity": -math.inf,
    }
    if value not in tokens:
        raise ValueError(f"unknown numeric token: {value!r}")
    return tokens[value]


def decode_point(values: list[int | float | str]) -> list[float]:
    if len(values) != 3:
        raise ValueError(f"coordinate must contain exactly three values: {values!r}")
    return [decode_number(value) for value in values]


def set_property(owner: Any, name: str, value: Any) -> None:
    if isinstance(value, bool):
        owner.SetBoolProp(name, value)
    elif isinstance(value, int):
        owner.SetIntProp(name, value)
    elif isinstance(value, float):
        owner.SetDoubleProp(name, value)
    else:
        owner.SetProp(name, str(value))


def configure_atom(atom: Chem.Atom, config: dict[str, Any]) -> None:
    atom.SetFormalCharge(int(config.get("formal_charge", 0)))
    atom.SetNumExplicitHs(int(config.get("explicit_hydrogens", 0)))
    tag_name = config.get("chiral_tag", "CHI_UNSPECIFIED")
    atom.SetChiralTag(getattr(Chem.ChiralType, tag_name))
    permutation = config.get("chiral_permutation")
    if permutation is not None:
        atom.SetUnsignedProp("_chiralPermutation", int(permutation))
    for name, value in config.get("props", {}).items():
        set_property(atom, name, value)


def configure_bond(bond: Chem.Bond, config: dict[str, Any]) -> None:
    direction = config.get("direction", "NONE")
    bond.SetBondDir(getattr(Chem.BondDir, direction))
    unknown_stereo = config.get("unknown_stereo")
    if unknown_stereo is not None:
        bond.SetIntProp("_UnknownStereo", int(unknown_stereo))
    for name, value in config.get("props", {}).items():
        set_property(bond, name, value)


def add_bond(mol: Chem.RWMol, config: dict[str, Any]) -> None:
    begin = int(config["begin"])
    end = int(config["end"])
    mol.AddBond(begin, end, getattr(Chem.BondType, config.get("type", "SINGLE")))
    configure_bond(mol.GetBondBetweenAtoms(begin, end), config)


def build_explicit_molecule(case: dict[str, Any]) -> Chem.RWMol:
    mol = Chem.RWMol()
    for config in case["atoms"]:
        atom = Chem.Atom(int(config["atomic_number"]))
        configure_atom(atom, config)
        mol.AddAtom(atom)
    for config in case["bonds"]:
        add_bond(mol, config)
    return mol


def build_single_center_molecule(case: dict[str, Any]) -> Chem.RWMol:
    mol = Chem.RWMol()
    center = case["center"]
    center_atom = Chem.Atom(int(center["atomic_number"]))
    configure_atom(center_atom, center)
    mol.AddAtom(center_atom)
    ligand_atomic_numbers = case["ligand_atomic_numbers"]
    for ligand_index, _ in enumerate(case["ligands"]):
        mol.AddAtom(Chem.Atom(int(ligand_atomic_numbers[ligand_index])))

    overrides = {int(item["ligand"]): item for item in case.get("bond_overrides", [])}
    for ligand_index, _ in enumerate(case["ligands"]):
        config = {
            "begin": 0,
            "end": ligand_index + 1,
            "type": case["bond_type"],
            "direction": case["bond_direction"],
            "unknown_stereo": case["unknown_stereo"],
        }
        config.update(overrides.get(ligand_index, {}))
        config.pop("ligand", None)
        if config.pop("reverse", False):
            config["begin"], config["end"] = config["end"], config["begin"]
        add_bond(mol, config)
    return mol


def add_conformers(mol: Chem.RWMol, case: dict[str, Any]) -> None:
    for conformer_config in case.get("conformers", []):
        if "coordinates" in conformer_config:
            coordinates = conformer_config["coordinates"]
        else:
            center = conformer_config.get("center_position", case["center_position"])
            ligands = conformer_config.get("ligands", case["ligands"])
            coordinates = [center, *ligands]
        if len(coordinates) != mol.GetNumAtoms():
            raise ValueError(
                f"{case['case_id']}: conformer row count {len(coordinates)} "
                f"does not match atom count {mol.GetNumAtoms()}"
            )
        conformer = Chem.Conformer(mol.GetNumAtoms())
        conformer.SetId(int(conformer_config["id"]))
        conformer.Set3D(bool(conformer_config["is_3d"]))
        for atom_index, coordinate in enumerate(coordinates):
            x, y, z = decode_point(coordinate)
            conformer.SetAtomPosition(atom_index, Point3D(x, y, z))
        mol.AddConformer(conformer, assignId=False)


def build_assign_molecule(defaults: dict[str, Any], raw_case: dict[str, Any]) -> tuple[Chem.RWMol, dict[str, Any]]:
    case = deep_merge(defaults, raw_case)
    if case["builder"] == "explicit":
        mol = build_explicit_molecule(case)
    elif case["builder"] == "single_center":
        mol = build_single_center_molecule(case)
    else:
        raise ValueError(f"unsupported fixture builder: {case['builder']!r}")
    for name, value in case.get("molecule_props", {}).items():
        set_property(mol, name, value)
    add_conformers(mol, case)
    mol.UpdatePropertyCache(strict=False)
    return mol, case


def json_value(value: Any) -> Any:
    if isinstance(value, float):
        return value.hex()
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    if isinstance(value, (tuple, list)):
        return [json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): json_value(item) for key, item in sorted(value.items())}
    return str(value)


def props(owner: Any) -> dict[str, Any]:
    return json_value(owner.GetPropsAsDict(includePrivate=True, includeComputed=False))


def snapshot_molecule(mol: Chem.Mol) -> dict[str, Any]:
    atoms = []
    for atom in mol.GetAtoms():
        atoms.append(
            {
                "index": atom.GetIdx(),
                "atomic_number": atom.GetAtomicNum(),
                "formal_charge": atom.GetFormalCharge(),
                "explicit_hydrogens": atom.GetNumExplicitHs(),
                "implicit_hydrogens": atom.GetNumImplicitHs(),
                "chiral_tag": str(atom.GetChiralTag()),
                "chiral_permutation": props(atom).get("_chiralPermutation"),
                "non_explicit_3d_chirality": props(atom).get(
                    "_NonExplicit3DChirality"
                ),
                "props": props(atom),
            }
        )
    bonds = []
    for bond in mol.GetBonds():
        bonds.append(
            {
                "index": bond.GetIdx(),
                "begin": bond.GetBeginAtomIdx(),
                "end": bond.GetEndAtomIdx(),
                "type": str(bond.GetBondType()),
                "direction": str(bond.GetBondDir()),
                "unknown_stereo": props(bond).get("_UnknownStereo"),
                "props": props(bond),
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
                "props": props(conformer),
            }
        )
    return {
        "atom_count": mol.GetNumAtoms(),
        "bond_count": mol.GetNumBonds(),
        "atoms": atoms,
        "bonds": bonds,
        "conformers": conformers,
        "molecule_props": props(mol),
        "stereochem_done": props(mol).get("_StereochemDone"),
    }


def selected_conformer_id(mol: Chem.Mol, conf_id: int) -> int | None:
    conformers = list(mol.GetConformers())
    if not conformers:
        return None
    if conf_id < 0:
        return conformers[0].GetId()
    for conformer in conformers:
        if conformer.GetId() == conf_id:
            return conf_id
    return None


def run_assign_case(payload: dict[str, Any]) -> dict[str, Any]:
    mol, case = build_assign_molecule(payload["defaults"], payload["case"])
    conf_id = int(case["conf_id"])
    before = snapshot_molecule(mol)
    selected_id = selected_conformer_id(mol, conf_id)
    status = "ok"
    error_type = None
    error_text = None
    try:
        Chem.AssignAtomChiralTagsFromStructure(
            mol,
            confId=conf_id,
            replaceExistingTags=bool(case["replace_existing_tags"]),
        )
    except Exception as error:
        status = "error"
        error_type = type(error).__name__
        error_text = str(error)
    return {
        "case_id": case["case_id"],
        "selection_reason": case["selection_reason"],
        "environment": case["environment"],
        "conf_id": conf_id,
        "replace_existing_tags": bool(case["replace_existing_tags"]),
        "status": status,
        "error_type": error_type,
        "error_text": error_text,
        "selected_conformer_id": selected_id,
        "before": before,
        "after": snapshot_molecule(mol),
    }


def octahedral_case(raw: dict[str, Any]) -> dict[str, Any]:
    axes = {
        "px": [1.0, 0.0, 0.0],
        "nx": [-1.0, 0.0, 0.0],
        "py": [0.0, 1.0, 0.0],
        "ny": [0.0, -1.0, 0.0],
        "pz": [0.0, 0.0, 1.0],
        "nz": [0.0, 0.0, -1.0],
    }
    return {
        "case_id": raw["case_id"],
        "selection_reason": (
            f"OctahedralPermFrom3D {raw['branch']} with {raw['volume_sign']} VOLTEST"
        ),
        "center": {"atomic_number": 15},
        "ligands": [axes[name] for name in raw["neighbor_order"]],
        "conformers": [{"id": 0, "is_3d": True}],
    }


def generate_assign_records(output: Path) -> None:
    fixture = json.loads(ASSIGN_FIXTURE.read_text(encoding="utf-8"))
    if fixture.get("schema_version") != 1:
        raise ValueError(f"unsupported assign fixture schema: {fixture.get('schema_version')!r}")
    cases = list(fixture["cases"])
    cases.extend(octahedral_case(case) for case in fixture["octahedral_switch_cases"])
    records = []
    for case in cases:
        merged = deep_merge(fixture["defaults"], case)
        environment = os.environ.copy()
        env_config = merged["environment"]
        if env_config["mode"] == "unset":
            environment.pop(NON_TETRAHEDRAL_ENV, None)
        elif env_config["mode"] == "set":
            environment[NON_TETRAHEDRAL_ENV] = env_config["value"]
        else:
            raise ValueError(f"unsupported environment mode: {env_config['mode']!r}")
        result = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--assign-case-worker"],
            input=json.dumps({"defaults": fixture["defaults"], "case": case}),
            cwd=Path(__file__).resolve().parents[3],
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"assign oracle worker failed for {case['case_id']}:\n{result.stderr}"
            )
        lines = [line for line in result.stdout.splitlines() if line.strip()]
        if len(lines) != 1:
            raise RuntimeError(
                f"assign oracle worker emitted {len(lines)} rows for {case['case_id']}"
            )
        record = json.loads(lines[0])
        if record.get("case_id") != case["case_id"]:
            raise RuntimeError(f"assign oracle worker case mismatch for {case['case_id']}")
        records.append(record)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=True, allow_nan=False))
            handle.write("\n")
    print(f"Wrote {len(records)} records to {output}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assign-case-worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("testdata/smiles/corpus/smiles_small.smi"),
        help="input SMILES file (default: testdata/smiles/corpus/smiles_small.smi)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("testdata/stereo/expected/rdkit/smiles_small/tetrahedral_stereo_geometry.jsonl"),
        help="output JSONL path (default: testdata/stereo/expected/rdkit/smiles_small/tetrahedral_stereo_geometry.jsonl)",
    )
    args = parser.parse_args()

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    if args.assign_case_worker:
        payload = json.loads(sys.stdin.read())
        print(json.dumps(run_assign_case(payload), ensure_ascii=True, allow_nan=False))
        return
    if args.output.name == ASSIGN_OUTPUT_NAME:
        generate_assign_records(args.output)
        return
    records = [record for smiles in iter_smiles(args.input) if (record := build_record(smiles))]
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w", encoding="utf-8") as f:
        for record in records:
            f.write(json.dumps(record, ensure_ascii=True))
            f.write("\n")

    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
