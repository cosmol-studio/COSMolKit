#!/usr/bin/env python3
"""Generate RDKit conformer-generation golden data."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDistGeom, rdGeometry


EXPECTED_RDKIT_VERSION = "2026.3.1"
REPO_ROOT = Path(__file__).resolve().parents[3]
INVENTORY_PATH = REPO_ROOT / "testdata/conformer/fixtures/rdkit_inventory.jsonl"


def assert_rdkit_version() -> None:
    actual = version("rdkit")
    assert actual == EXPECTED_RDKIT_VERSION, (
        f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, got {actual}"
    )


def load_fixture_inventory() -> dict[str, Path]:
    mapping: dict[str, Path] = {}
    for line in INVENTORY_PATH.read_text(encoding="utf-8").splitlines():
        record = json.loads(line)
        fixture = record.get("fixture")
        if fixture is None:
            continue
        mapping[record["source"]] = REPO_ROOT / "testdata/conformer/fixtures" / fixture
    return mapping


FIXTURE_INVENTORY = load_fixture_inventory()


def load_fixture_mol(source: str) -> Chem.Mol:
    path = FIXTURE_INVENTORY.get(source, REPO_ROOT / source)
    text = path.read_text(encoding="utf-8")
    mol = Chem.MolFromMolBlock(text, sanitize=True, removeHs=False)
    if mol is None:
        raise ValueError(f"RDKit failed to parse fixture {source}")
    return mol


def load_smiles_mol(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        raise ValueError(f"RDKit failed to parse SMILES {smiles}")
    return mol


def preset_default() -> rdDistGeom.EmbedParameters:
    return rdDistGeom.EmbedParameters()


def preset_kdg() -> rdDistGeom.EmbedParameters:
    return rdDistGeom.KDG()


def preset_etdg() -> rdDistGeom.EmbedParameters:
    return rdDistGeom.ETDG()


def preset_etdg_v2() -> rdDistGeom.EmbedParameters:
    return rdDistGeom.ETDGv2()


def preset_etkdg() -> rdDistGeom.EmbedParameters:
    return rdDistGeom.ETKDG()


def preset_etkdg_v2() -> rdDistGeom.EmbedParameters:
    return rdDistGeom.ETKDGv2()


def preset_etkdg_v3() -> rdDistGeom.EmbedParameters:
    return rdDistGeom.ETKDGv3()


def preset_sr_etkdg_v3() -> rdDistGeom.EmbedParameters:
    return rdDistGeom.srETKDGv3()


def coords_for_conformer(mol: Chem.Mol, conf_id: int) -> list[list[float]]:
    conf = mol.GetConformer(conf_id)
    return [
        [float(conf.GetAtomPosition(idx)[axis]) for axis in range(3)]
        for idx in range(mol.GetNumAtoms())
    ]


def all_conformer_coords(mol: Chem.Mol, conf_ids: list[int]) -> list[list[list[float]]]:
    return [coords_for_conformer(mol, conf_id) for conf_id in conf_ids]


def failure_counts(params: rdDistGeom.EmbedParameters) -> list[int] | None:
    if not hasattr(params, "GetFailureCounts"):
        return None
    return [int(value) for value in params.GetFailureCounts()]


def case_record(case: dict[str, object]) -> dict[str, object]:
    params = case["preset"]()
    for key, value in case.get("attrs", {}).items():
        setattr(params, key, value)
    if "coord_map" in case:
        params.SetCoordMap(
            {
                int(atom_idx): rdGeometry.Point3D(*xyz)
                for atom_idx, xyz in case["coord_map"].items()
            }
        )
    if "cpci" in case:
        params.SetCPCI(
            {
                tuple(int(part) for part in key.split("-")): float(value)
                for key, value in case["cpci"].items()
            }
        )
    param_hook = case.get("param_hook")
    if param_hook is not None:
        param_hook(params)

    loader = case["loader"]
    mol = loader(case["source"])

    record: dict[str, object] = {
        "case_id": case["case_id"],
        "mode": case["mode"],
        "source_kind": case["source_kind"],
        "source": case["source"],
        "preset": case["preset_name"],
        "attrs": case.get("attrs", {}),
        "rdkit_ok": True,
        "status": None,
        "ids": None,
        "conformers": None,
        "failure_counts": None,
        "error": None,
    }

    try:
        if case["mode"] == "single":
            status = int(rdDistGeom.EmbedMolecule(mol, params))
            record["status"] = status
            if status >= 0:
                record["conformers"] = [coords_for_conformer(mol, status)]
            else:
                record["conformers"] = []
        elif case["mode"] == "multi":
            ids = [int(value) for value in rdDistGeom.EmbedMultipleConfs(mol, case["num_confs"], params)]
            record["ids"] = ids
            record["conformers"] = all_conformer_coords(mol, ids)
        else:
            raise ValueError(f"unsupported mode {case['mode']}")
        record["failure_counts"] = failure_counts(params)
    except Exception as err:
        record["rdkit_ok"] = False
        record["error"] = str(err)

    return record


CASES: list[dict[str, object]] = [
    {
        "case_id": "single_dg_simple_torsion",
        "mode": "single",
        "source_kind": "fixture_mol",
        "source": "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.dg.mol",
        "loader": load_fixture_mol,
        "preset_name": "DG",
        "preset": preset_default,
        "attrs": {"randomSeed": 42, "numThreads": 1},
    },
    {
        "case_id": "single_kdg_simple_torsion",
        "mode": "single",
        "source_kind": "fixture_mol",
        "source": "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.kdg.mol",
        "loader": load_fixture_mol,
        "preset_name": "KDG",
        "preset": preset_kdg,
        "attrs": {"randomSeed": 42, "numThreads": 1},
    },
    {
        "case_id": "single_etdg_simple_torsion",
        "mode": "single",
        "source_kind": "fixture_mol",
        "source": "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etdg.mol",
        "loader": load_fixture_mol,
        "preset_name": "ETDG",
        "preset": preset_etdg,
        "attrs": {"randomSeed": 42, "numThreads": 1},
    },
    {
        "case_id": "single_etdgv2_torsion",
        "mode": "single",
        "source_kind": "fixture_mol",
        "source": "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/torsion.etkdg.v2.mol",
        "loader": load_fixture_mol,
        "preset_name": "ETDGv2",
        "preset": preset_etdg_v2,
        "attrs": {"randomSeed": 42, "numThreads": 1},
    },
    {
        "case_id": "single_etkdg_simple_torsion",
        "mode": "single",
        "source_kind": "fixture_mol",
        "source": "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etkdg.mol",
        "loader": load_fixture_mol,
        "preset_name": "ETKDG",
        "preset": preset_etkdg,
        "attrs": {"randomSeed": 42, "numThreads": 1},
    },
    {
        "case_id": "single_etkdgv2_torsion",
        "mode": "single",
        "source_kind": "fixture_mol",
        "source": "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/torsion.etkdg.v2.mol",
        "loader": load_fixture_mol,
        "preset_name": "ETKDGv2",
        "preset": preset_etkdg_v2,
        "attrs": {"randomSeed": 42, "numThreads": 1},
    },
    {
        "case_id": "single_etkdgv3_macrocycle",
        "mode": "single",
        "source_kind": "fixture_mol",
        "source": "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.macrocycle.etkdgv3.mol",
        "loader": load_fixture_mol,
        "preset_name": "ETKDGv3",
        "preset": preset_etkdg_v3,
        "attrs": {"randomSeed": 42, "numThreads": 1},
    },
    {
        "case_id": "single_sretkdgv3_smallring",
        "mode": "single",
        "source_kind": "fixture_mol",
        "source": "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.smallring.etkdgv3.mol",
        "loader": load_fixture_mol,
        "preset_name": "srETKDGv3",
        "preset": preset_sr_etkdg_v3,
        "attrs": {"randomSeed": 42, "numThreads": 1},
    },
    {
        "case_id": "single_chirality_failure_fixture",
        "mode": "single",
        "source_kind": "fixture_mol",
        "source": "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/chirality_failure_test.mol",
        "loader": load_fixture_mol,
        "preset_name": "ETKDGv3",
        "preset": preset_etkdg_v3,
        "attrs": {"randomSeed": 61453, "numThreads": 1, "trackFailures": True, "maxIterations": 50},
    },
    {
        "case_id": "single_coordmap_randomcoords",
        "mode": "single",
        "source_kind": "fixture_mol",
        "source": "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etkdg.mol",
        "loader": load_fixture_mol,
        "preset_name": "ETKDGv3",
        "preset": preset_etkdg_v3,
        "attrs": {"randomSeed": 12648430, "numThreads": 1, "useRandomCoords": True},
        "coord_map": {
            0: [0.0, 0.0, 0.0],
            1: [0.0, 0.0, 1.5],
            2: [0.0, 1.5, 1.5],
        },
    },
    {
        "case_id": "single_cpci_etkdgv3",
        "mode": "single",
        "source_kind": "fixture_mol",
        "source": "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etkdg.mol",
        "loader": load_fixture_mol,
        "preset_name": "ETKDGv3",
        "preset": preset_etkdg_v3,
        "attrs": {"randomSeed": 12648430, "numThreads": 1},
        "cpci": {
            "0-3": 0.5,
            "1-4": -0.25,
        },
    },
    {
        "case_id": "multi_embed_fragments_separately",
        "mode": "multi",
        "num_confs": 3,
        "source_kind": "smiles",
        "source": "CC.CC",
        "loader": load_smiles_mol,
        "preset_name": "ETKDG",
        "preset": preset_etkdg,
        "attrs": {
            "randomSeed": 61453,
            "numThreads": 1,
            "embedFragmentsSeparately": True,
        },
    },
    {
        "case_id": "multi_etkdg_seeded",
        "mode": "multi",
        "num_confs": 10,
        "source_kind": "fixture_mol",
        "source": "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etkdg.mol",
        "loader": load_fixture_mol,
        "preset_name": "ETKDG",
        "preset": preset_etkdg,
        "attrs": {"randomSeed": 61453, "numThreads": 1, "trackFailures": True, "timeout": 1},
    },
    {
        "case_id": "multi_etkdg_pruned_symmetry",
        "mode": "multi",
        "num_confs": 10,
        "source_kind": "fixture_mol",
        "source": "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etkdg.mol",
        "loader": load_fixture_mol,
        "preset_name": "ETKDG",
        "preset": preset_etkdg,
        "attrs": {
            "randomSeed": 61453,
            "numThreads": 1,
            "pruneRmsThresh": 0.5,
            "useSymmetryForPruning": True,
        },
    },
    {
        "case_id": "multi_etkdg_sequential_seeds",
        "mode": "multi",
        "num_confs": 4,
        "source_kind": "fixture_mol",
        "source": "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etkdg.mol",
        "loader": load_fixture_mol,
        "preset_name": "ETKDG",
        "preset": preset_etkdg,
        "attrs": {
            "randomSeed": 61453,
            "numThreads": 1,
            "enableSequentialRandomSeeds": True,
        },
    },
    {
        "case_id": "multi_force_trans_amides_true",
        "mode": "multi",
        "num_confs": 10,
        "source_kind": "smiles",
        "source": "CC(=O)NC",
        "loader": load_smiles_mol,
        "preset_name": "DG",
        "preset": preset_default,
        "attrs": {
            "randomSeed": 61453,
            "numThreads": 1,
            "forceTransAmides": True,
            "useExpTorsionAnglePrefs": False,
            "useBasicKnowledge": True,
        },
    },
    {
        "case_id": "multi_force_trans_amides_false",
        "mode": "multi",
        "num_confs": 10,
        "source_kind": "smiles",
        "source": "CC(=O)NC",
        "loader": load_smiles_mol,
        "preset_name": "DG",
        "preset": preset_default,
        "attrs": {
            "randomSeed": 61453,
            "numThreads": 1,
            "forceTransAmides": False,
            "useExpTorsionAnglePrefs": False,
            "useBasicKnowledge": True,
        },
    },
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("testdata/conformer/expected/rdkit/smiles_small/conformer_generation.jsonl"),
        help="output JSONL path (default: testdata/conformer/expected/rdkit/smiles_small/conformer_generation.jsonl)",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="unused compatibility flag accepted for generate_all.py",
    )
    args = parser.parse_args()

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    records = [case_record(case) for case in CASES]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=True))
            handle.write("\n")
    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
