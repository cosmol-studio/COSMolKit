"""Original 19 conformer cases through the public Python API.

Set COSMOLKIT_CONFORMER_ORACLE_19 to the pinned JSONL oracle when it is
not installed under testdata/conformer/expected/rdkit/smiles_small.
Fixture cases require the public SDF reader with the original IO48 flags.
"""
import hashlib
import json
import os
import tempfile
from pathlib import Path

import cosmolkit as ck
import pytest

WORKSPACE = Path(__file__).resolve().parents[2]
REFERENCE = Path(os.environ.get(
    "COSMOLKIT_CONFORMER_ORACLE_19",
    WORKSPACE / "testdata/conformer/expected/rdkit/smiles_small/conformer_generation.jsonl",
))
if not REFERENCE.is_file():
    raise FileNotFoundError(
        f"Missing original19 oracle: {REFERENCE}; set COSMOLKIT_CONFORMER_ORACLE_19 "
        "to the pinned original JSONL file. Do not regenerate the reference."
    )
OUTPUT = (Path(os.environ["COSMOLKIT_CONFORMER_NATIVE_OUTPUT"])
          if "COSMOLKIT_CONFORMER_NATIVE_OUTPUT" in os.environ
          else Path(tempfile.mkdtemp(prefix="cosmolkit-conformer19-")))
RAW = REFERENCE.read_bytes()
assert hashlib.sha256(RAW).hexdigest() == "710255e3ace6d263cbbe47e7e68dd71db2cb5ca8986456c78c042bbd9dc6c352"
RECORDS = [json.loads(line) for line in RAW.splitlines()]
assert len(RECORDS) == 19
# Original fixture identities, independent of generated audit artifacts.
FIXTURES = {
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/Issue3238580.1.mol": {
        "fixture": "rdkit/test_data/Issue3238580.1.mol",
        "sha256": "a2f85be9181d1c86c728eca3dc0686b4e1b78fd5f1dbce5e118ca2dc3fd25117",
        "bytes": 1204
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/Issue3238580.2.mol": {
        "fixture": "rdkit/test_data/Issue3238580.2.mol",
        "sha256": "a2efdac34209553a196d568ede259b4cba24137d792740aa86ac8e0caf49a301",
        "bytes": 1486
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/Issue3238580.3.mol": {
        "fixture": "rdkit/test_data/Issue3238580.3.mol",
        "sha256": "2e19b95df2811e52f50f480fe93d2a4630588aa89aece311565999a99a69bcb9",
        "bytes": 1486
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/Issue3483968.mol": {
        "fixture": "rdkit/test_data/Issue3483968.mol",
        "sha256": "7e3025087d79ec44520ffd8beb437aaf4f63a9c43c2a3212040c47d680bcb8e2",
        "bytes": 3431
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/atropisomers.sdf": {
        "fixture": "rdkit/test_data/atropisomers.sdf",
        "sha256": "2a9a8af87fbebeaaa5a959b776755d1bc21991ec856d252bb5f3275377ad9c2d",
        "bytes": 8444
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/chirality_failure_test.mol": {
        "fixture": "rdkit/test_data/chirality_failure_test.mol",
        "sha256": "91c26205a7648da730c5e677bd6911402b78f3a2c0a5552e1264b56ebc01c327",
        "bytes": 2160
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/cis_trans_cases.csv": {
        "fixture": "rdkit/test_data/cis_trans_cases.csv",
        "sha256": "89780520fe7ba0750def25db0c9d168cb304d592c1b5188c3952d5275165d1d4",
        "bytes": 9029
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/combi_coords.sdf": {
        "fixture": "rdkit/test_data/combi_coords.sdf",
        "sha256": "a951d20bbe09bd87c9d602ea53b6380ab1bb3998f5908061ea53cecb38e8d070",
        "bytes": 208327
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/constrain1.sdf": {
        "fixture": "rdkit/test_data/constrain1.sdf",
        "sha256": "33166694858fea066fbe8450cb450503c9d0e3eb0979cbe3ab03f8e7592244b3",
        "bytes": 2087
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/embedDistOpti.sdf": {
        "fixture": "rdkit/test_data/embedDistOpti.sdf",
        "sha256": "0c110aa347055bb0e6d4fcfec2074685f2fdf7bfb14afe54f1ecca7547472059",
        "bytes": 225825
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/initCoords.random.sdf": {
        "fixture": "rdkit/test_data/initCoords.random.sdf",
        "sha256": "e2ee890eede4a7d15ceed72952adf74c9d78f77682ff6c377f9dca1773d3a4f4",
        "bytes": 15596
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/initCoords.sdf": {
        "fixture": "rdkit/test_data/initCoords.sdf",
        "sha256": "a89538f7aa29d7881da8c5dd29d147572bdf071b4060a97c749a907d5d7e2952",
        "bytes": 7379
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.dg.mol": {
        "fixture": "rdkit/test_data/simple_torsion.dg.mol",
        "sha256": "859e5f22dc1f2a3f090114e4e9fb1282000a3fe2034c7fd2774c411fb4e29d66",
        "bytes": 1055
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etdg.mol": {
        "fixture": "rdkit/test_data/simple_torsion.etdg.mol",
        "sha256": "40c67906d3e421703beba6fcf43c556e308ef3739a9f29a81ee37db773941684",
        "bytes": 1055
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etkdg.mol": {
        "fixture": "rdkit/test_data/simple_torsion.etkdg.mol",
        "sha256": "40c67906d3e421703beba6fcf43c556e308ef3739a9f29a81ee37db773941684",
        "bytes": 1055
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.etkdg.new.mol": {
        "fixture": "rdkit/test_data/simple_torsion.etkdg.new.mol",
        "sha256": "40c67906d3e421703beba6fcf43c556e308ef3739a9f29a81ee37db773941684",
        "bytes": 1055
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.kdg.mol": {
        "fixture": "rdkit/test_data/simple_torsion.kdg.mol",
        "sha256": "859e5f22dc1f2a3f090114e4e9fb1282000a3fe2034c7fd2774c411fb4e29d66",
        "bytes": 1055
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.macrocycle.etkdg.mol": {
        "fixture": "rdkit/test_data/simple_torsion.macrocycle.etkdg.mol",
        "sha256": "de17114b8d8de6637b5a5b9eb60186d625ac1d1c49530d6d3b66d0efae9cbed7",
        "bytes": 2645
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.macrocycle.etkdgv3.mol": {
        "fixture": "rdkit/test_data/simple_torsion.macrocycle.etkdgv3.mol",
        "sha256": "b04a72fc4b5a5b7bf7d1112c48fda4fdbfa7797b2fbca2b9d0d7cf4801ef1eee",
        "bytes": 2728
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.macrocycle1.etkdg.mol": {
        "fixture": "rdkit/test_data/simple_torsion.macrocycle1.etkdg.mol",
        "sha256": "64d009826a3e6086a01976682acf0a04855f212b611986cfa452fe1dd9b3638a",
        "bytes": 2645
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/simple_torsion.smallring.etkdgv3.mol": {
        "fixture": "rdkit/test_data/simple_torsion.smallring.etkdgv3.mol",
        "sha256": "85a2e16f2cca6d7efc0d0af50f08416fb03fc96ac8eac18f024d105c37100ae3",
        "bytes": 1566
    },
    "third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/torsion.etkdg.v2.mol": {
        "fixture": "rdkit/test_data/torsion.etkdg.v2.mol",
        "sha256": "9c2f58d523a454d673e3b26dab353d6b631c1e500004b82b3088d79873f52087",
        "bytes": 1483
    }
}
FACTORIES = {"DG": "dg", "KDG": "kdg", "ETDG": "etdg", "ETDGv2": "etdg_v2", "ETKDG": "etkdg", "ETKDGv2": "etkdg_v2", "ETKDGv3": "etkdg_v3", "srETKDGv3": "sr_etkdg_v3"}
COUNTS = {"multi_embed_fragments_separately": 3, "multi_etkdg_seeded": 10, "multi_etkdg_pruned_symmetry": 10, "multi_etkdg_sequential_seeds": 4, "multi_force_trans_amides_true": 10, "multi_force_trans_amides_false": 10}


def molecule(record):
    kind, source = record["source_kind"], record["source"]
    if kind == "fixture_mol":
        identity = FIXTURES[source]
        raw = (WORKSPACE / "testdata/conformer/fixtures" / identity["fixture"]).read_bytes()
        assert len(raw) == identity["bytes"]
        assert hashlib.sha256(raw).hexdigest() == identity["sha256"]
        params = ck.SdfReadParams(sanitize=True, remove_hydrogens=False, process_property_lists=False)
        return ck.Molecule.from_sdf_with_params(raw.decode(), params)
    result = ck.Molecule.from_smiles(source)
    if kind == "smiles_with_hydrogens":
        return result.with_hydrogens()
    assert kind == "smiles"
    return result


def parameters(record):
    params = getattr(ck.EmbedParams, FACTORIES[record["preset"]])()
    name = record["case_id"]
    params = configured(params, num_threads=1)
    if name in {"single_dg_simple_torsion", "single_kdg_simple_torsion", "single_etdg_simple_torsion", "single_etdgv2_torsion", "single_etkdg_simple_torsion", "single_etkdgv2_torsion", "single_etkdgv3_macrocycle", "single_sretkdgv3_smallring"}:
        params = configured(params, random_seed=42)
    elif name == "single_chirality_failure_fixture":
        params = configured(params, random_seed=0xF00D)
        params = configured(params, track_failures=True)
        params = configured(params, max_iterations=50)
    elif name == "single_coordmap_randomcoords":
        params = configured(params, random_seed=0xC0FFEE)
        params = configured(params, use_random_coords=True)
        params = configured(params, coord_map={0: [0.0, 0.0, 0.0], 1: [0.0, 0.0, 1.5], 2: [0.0, 1.5, 1.5]})
    elif name == "single_cpci_etkdgv3":
        params = configured(params, random_seed=0xC0FFEE)
        params = configured(params, cpci={(0, 3): 0.5, (1, 4): -0.25})
    elif name in {"single_etkdgv3_x0_ring_connectivity_first", "single_etkdgv3_x0_ring_connectivity_second"}:
        params = configured(params, max_iterations=3)
        params = configured(params, random_seed=61453)
        params = configured(params, timeout=0)
    elif name in COUNTS:
        params = configured(params, random_seed=0xF00D)
        if name == "multi_embed_fragments_separately":
            params = configured(params, embed_fragments_separately=True)
        elif name == "multi_etkdg_seeded":
            params = configured(params, track_failures=True)
            params = configured(params, timeout=1)
        elif name == "multi_etkdg_pruned_symmetry":
            params = configured(params, prune_rms_thresh=0.5)
            params = configured(params, use_symmetry_for_pruning=True)
        elif name == "multi_etkdg_sequential_seeds":
            params = configured(params, enable_sequential_random_seeds=True)
        else:
            params = configured(params, force_trans_amides=name == "multi_force_trans_amides_true")
            params = configured(params, use_exp_torsion_angle_prefs=False)
            params = configured(params, use_basic_knowledge=True)
    else:
        raise AssertionError(f"unknown original case {name}")
    return params


@pytest.mark.parametrize("record", RECORDS, ids=[record["case_id"] for record in RECORDS])
def test_original19_native(record):
    actual = {"case_id": record["case_id"], "matched": False, "stage": "source"}
    try:
        assert record["rdkit_ok"], record.get("error")
        source = molecule(record)
        params = parameters(record)
        actual["stage"] = "generate"
        if record["mode"] == "single":
            result = source.with_3d_conformer_result(params)
            actual["status"] = result.conf_id()
            assert actual["status"] == record["status"]
        else:
            assert record["mode"] == "multi"
            result = source.with_3d_conformers_result(COUNTS[record["case_id"]], params)
            actual["ids"] = result.conf_ids()
            assert actual["ids"] == record["ids"]
        actual["failure_counts"] = result.params().failures
        assert actual["failure_counts"] == (record.get("failure_counts") or [])
        actual["conformers"] = [conf.coordinates() for conf in result.molecule().conformers_3d()]
        expected = record["conformers"]
        assert len(actual["conformers"]) == len(expected)
        for conf, reference in zip(actual["conformers"], expected):
            assert len(conf) == len(reference)
            for atom, reference_atom in zip(conf, reference):
                assert len(atom) == 3
                for axis in range(3):
                    assert abs(atom[axis] - reference_atom[axis]) <= 1.0e-6
        actual["matched"] = True
        actual["stage"] = "complete"
    except Exception as error:
        actual["error"] = {"type": type(error).__name__, "message": str(error)}
        raise
    finally:
        OUTPUT.mkdir(parents=True, exist_ok=True)
        path = OUTPUT / (record["case_id"] + ".json")
        assert not path.exists(), "actual outputs require a new per-run directory"
        path.write_text(json.dumps(actual, indent=2) + "\n")


def configured(params, **changes):
    values = {field: getattr(params, field) for field in ['max_iterations', 'num_threads', 'random_seed', 'clear_confs', 'use_random_coords', 'box_size_mult', 'rand_neg_eig', 'num_zero_fail', 'coord_map', 'optimizer_force_tol', 'ignore_smoothing_failures', 'enforce_chirality', 'use_exp_torsion_angle_prefs', 'use_basic_knowledge', 'verbose', 'basin_thresh', 'prune_rms_thresh', 'only_heavy_atoms_for_rms', 'et_version', 'embed_fragments_separately', 'use_small_ring_torsions', 'use_macrocycle_torsions', 'use_macrocycle14config', 'timeout', 'cpci', 'force_trans_amides', 'use_symmetry_for_pruning', 'bounds_mat_force_scaling', 'track_failures', 'enable_sequential_random_seeds', 'symmetrize_conjugated_terminal_groups_for_pruning']}
    values.update(changes)
    return ck.EmbedParams(**values)
