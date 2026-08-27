"""Shared pinned-RDKit tautomer oracle used by fixtures and corpus audits."""

from __future__ import annotations

import csv
import gzip
import json
from importlib.metadata import version
from pathlib import Path
from typing import Any, Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize


EXPECTED_RDKIT_VERSION = "2026.3.1"
EXPECTED_RDKIT_SOURCE_REVISION = "351f8f378f8ad6bbd517980c38896e66bf907af8c"
PROFILE_PATH = Path(__file__).with_name("tautomer_profile.json")
RDLogger.DisableLog("rdApp.*")


def assert_rdkit_version() -> None:
    actual = version("rdkit")
    if actual != EXPECTED_RDKIT_VERSION:
        raise SystemExit(
            f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, got {actual}"
        )


def load_profile(path: Path = PROFILE_PATH) -> dict[str, Any]:
    profile = json.loads(path.read_text(encoding="utf-8"))
    if profile.get("schema_version") != 1:
        raise ValueError("tautomer profile schema_version must be 1")
    if profile.get("rdkit_version") != EXPECTED_RDKIT_VERSION:
        raise ValueError("tautomer profile has the wrong RDKit version")
    if profile.get("source_revision") != EXPECTED_RDKIT_SOURCE_REVISION:
        raise ValueError("tautomer profile has the wrong RDKit source revision")
    branches = profile.get("branches")
    if not isinstance(branches, list) or not branches:
        raise ValueError("tautomer profile must contain branches")
    names = [branch.get("name") for branch in branches]
    if any(not isinstance(name, str) or not name for name in names):
        raise ValueError("tautomer branch names must be non-empty strings")
    if len(names) != len(set(names)):
        raise ValueError("tautomer branch names must be unique")
    return profile


def iter_input_records(path: Path) -> Iterable[dict[str, Any]]:
    if path.suffix == ".json":
        document = json.loads(path.read_text(encoding="utf-8"))
        for row, case in enumerate(document["cases"]):
            yield {"row": row, "sanitize": True, "remove_hs": True, **case}
        return
    if path.name.endswith(".csv.gz"):
        with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
            for row, columns in enumerate(csv.reader(stream)):
                if len(columns) != 2:
                    raise ValueError(f"{path}:{row + 1}: expected exactly two columns")
                yield {
                    "row": row,
                    "case_id": f"{path.stem}:{row}",
                    "smiles": columns[0],
                    "expected_canonical_smiles": columns[1],
                    "sanitize": True,
                    "remove_hs": True,
                    "source": str(path),
                }
        return
    row = 0
    for raw in path.read_text(encoding="utf-8").splitlines():
        smiles = raw.strip()
        if not smiles or smiles.startswith("#"):
            continue
        yield {
            "row": row,
            "case_id": f"{path.stem}:{row}",
            "smiles": smiles,
            "sanitize": True,
            "remove_hs": True,
            "source": str(path),
        }
        row += 1


def error_record(error: BaseException) -> dict[str, str]:
    return {"type": type(error).__name__, "message": str(error)}


def molecule_state(molecule: Chem.Mol) -> dict[str, Any]:
    atoms = []
    for atom in molecule.GetAtoms():
        atoms.append(
            {
                "atomic_number": atom.GetAtomicNum(),
                "formal_charge": atom.GetFormalCharge(),
                "explicit_hydrogens": atom.GetNumExplicitHs(),
                "no_implicit": atom.GetNoImplicit(),
                "isotope": atom.GetIsotope(),
                "radical_electrons": atom.GetNumRadicalElectrons(),
                "aromatic": atom.GetIsAromatic(),
                "chiral_tag": str(atom.GetChiralTag()),
                "hybridization": str(atom.GetHybridization()),
                "cip_code": atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None,
            }
        )
    bonds = []
    stereo_names = {
        "STEREONONE": "NONE",
        "STEREOANY": "ANY",
        "STEREOZ": "Z",
        "STEREOE": "E",
        "STEREOCIS": "CIS",
        "STEREOTRANS": "TRANS",
        "STEREOATROPCW": "ATROP_CW",
        "STEREOATROPCCW": "ATROP_CCW",
    }
    for bond in molecule.GetBonds():
        bonds.append(
            {
                "begin": bond.GetBeginAtomIdx(),
                "end": bond.GetEndAtomIdx(),
                "bond_type": str(bond.GetBondType()),
                "aromatic": bond.GetIsAromatic(),
                "conjugated": bond.GetIsConjugated(),
                "direction": str(bond.GetBondDir()),
                "stereo": stereo_names[str(bond.GetStereo())],
                "stereo_atoms": list(bond.GetStereoAtoms()),
            }
        )
    return {
        "isomeric_smiles": Chem.MolToSmiles(molecule, isomericSmiles=True),
        "atoms": atoms,
        "bonds": bonds,
    }


def make_enumerator(branch: dict[str, Any]) -> Any:
    if branch["catalog"] == "v1":
        enumerator = rdMolStandardize.GetV1TautomerEnumerator()
        enumerator.SetMaxTautomers(branch["max_tautomers"])
        enumerator.SetMaxTransforms(branch["max_transforms"])
        enumerator.SetRemoveSp3Stereo(branch["remove_sp3_stereo"])
        enumerator.SetRemoveBondStereo(branch["remove_bond_stereo"])
        enumerator.SetReassignStereo(branch["reassign_stereo"])
        return enumerator
    parameters = rdMolStandardize.CleanupParameters()
    parameters.maxTautomers = branch["max_tautomers"]
    parameters.maxTransforms = branch["max_transforms"]
    parameters.tautomerRemoveSp3Stereo = branch["remove_sp3_stereo"]
    parameters.tautomerRemoveBondStereo = branch["remove_bond_stereo"]
    parameters.tautomerRemoveIsotopicHs = branch["remove_isotopic_hydrogens"]
    parameters.tautomerReassignStereo = branch["reassign_stereo"]
    return rdMolStandardize.TautomerEnumerator(parameters)


def score_state(molecule: Chem.Mol, enumerator: Any) -> dict[str, int]:
    ring = rdMolStandardize.ScoreRings(molecule)
    substructure = rdMolStandardize.ScoreSubstructs(molecule)
    hetero_hydrogen = rdMolStandardize.ScoreHeteroHs(molecule)
    return {
        "ring": ring,
        "substructure": substructure,
        "hetero_hydrogen": hetero_hydrogen,
        "total": enumerator.ScoreTautomer(molecule),
    }


def parse_molecule(record: dict[str, Any]) -> Chem.Mol | None:
    parameters = Chem.SmilesParserParams()
    parameters.sanitize = bool(record.get("sanitize", True))
    parameters.removeHs = bool(record.get("remove_hs", True))
    return Chem.MolFromSmiles(record["smiles"], parameters)


def enumerate_branch(molecule: Chem.Mol, branch: dict[str, Any]) -> dict[str, Any]:
    enumerator = make_enumerator(branch)
    try:
        result = enumerator.Enumerate(molecule)
        tautomers = list(result)
        canonical = enumerator.PickCanonical(result)
        return {
            "ok": True,
            "error": None,
            "status": str(result.status),
            "ordered_smiles": list(result.smiles),
            "modified_atoms": list(result.modifiedAtoms),
            "modified_bonds": list(result.modifiedBonds),
            "molecule_states": [molecule_state(value) for value in tautomers],
            "scores": [score_state(value, enumerator) for value in tautomers],
            "canonical_smiles": Chem.MolToSmiles(canonical, isomericSmiles=True),
            "canonical_state": molecule_state(canonical),
        }
    except Exception as error:  # noqa: BLE001
        return {
            "ok": False,
            "error": error_record(error),
            "status": None,
            "ordered_smiles": [],
            "modified_atoms": [],
            "modified_bonds": [],
            "molecule_states": [],
            "scores": [],
            "canonical_smiles": None,
            "canonical_state": None,
        }


def canonicalize_branch(molecule: Chem.Mol, branch: dict[str, Any]) -> dict[str, Any]:
    enumerator = make_enumerator(branch)
    try:
        canonical = enumerator.Canonicalize(molecule)
        return {
            "ok": True,
            "error": None,
            "canonical_smiles": Chem.MolToSmiles(
                canonical, isomericSmiles=True
            ),
            "canonical_state": molecule_state(canonical),
            "canonical_score": score_state(canonical, enumerator),
        }
    except Exception as error:  # noqa: BLE001
        return {
            "ok": False,
            "error": error_record(error),
            "canonical_smiles": None,
            "canonical_state": None,
            "canonical_score": None,
        }


def endpoint_input_observation(
    smiles: str, profile: dict[str, Any], branch_names: list[str]
) -> dict[str, Any]:
    record = {
        "smiles": smiles,
        "parse": {"ok": True, "error": None},
        "branches": {},
    }
    try:
        molecule = Chem.MolFromSmiles(smiles)
    except Exception as error:  # noqa: BLE001
        return {
            **record,
            "parse": {"ok": False, "error": error_record(error)},
        }
    if molecule is None:
        return {
            **record,
            "parse": {
                "ok": False,
                "error": {
                    "type": "NullMolecule",
                    "message": "MolFromSmiles returned None",
                },
            },
        }
    by_name = {branch["name"]: branch for branch in profile["branches"]}
    return {
        **record,
        "branches": {
            name: {
                "parameters": by_name[name],
                **canonicalize_branch(molecule, by_name[name]),
            }
            for name in branch_names
        },
    }


def build_record(
    record: dict[str, Any], profile: dict[str, Any], branch_names: list[str]
) -> dict[str, Any]:
    base = {
        "schema_version": 1,
        "row": int(record["row"]),
        "case_id": str(record["case_id"]),
        "smiles": str(record["smiles"]),
        "sanitize": bool(record.get("sanitize", True)),
        "remove_hs": bool(record.get("remove_hs", True)),
        "source": str(record.get("source", "")),
        "expected_canonical_smiles": record.get("expected_canonical_smiles"),
    }
    try:
        molecule = parse_molecule(record)
    except Exception as error:  # noqa: BLE001
        return {**base, "parse": {"ok": False, "error": error_record(error)}, "branches": {}}
    if molecule is None:
        return {
            **base,
            "parse": {"ok": False, "error": {"type": "NullMolecule", "message": "MolFromSmiles returned None"}},
            "branches": {},
        }
    by_name = {branch["name"]: branch for branch in profile["branches"]}
    result = {
        **base,
        "parse": {"ok": True, "error": None},
        "branches": {
            name: {"parameters": by_name[name], **enumerate_branch(molecule, by_name[name])}
            for name in branch_names
        },
    }
    paired_smiles = record.get("expected_canonical_smiles")
    if paired_smiles is not None:
        result["input_tautomer"] = endpoint_input_observation(
            str(paired_smiles), profile, branch_names
        )
        permutation = Chem.MolToRandomSmilesVect(
            molecule,
            1,
            randomSeed=int(record["row"]) + 1,
            isomericSmiles=True,
            kekuleSmiles=False,
            allBondsExplicit=False,
            allHsExplicit=False,
        )[0]
        result["atom_order_permutation"] = endpoint_input_observation(
            permutation, profile, branch_names
        )
    return result
