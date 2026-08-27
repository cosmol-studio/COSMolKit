from __future__ import annotations

import json
from pathlib import Path

import _generate_tautomer_catalog_golden as catalog_generator
import _generate_tautomer_golden as tautomer_generator
import _tautomer_oracle as oracle
from _expected_schema import SCHEMAS, validate_jsonl_output
from generate_all import GENERATOR_SPECS, PROFILE_INPUTS, REPO_ROOT


def test_tautomer_generators_and_profiles_are_registered() -> None:
    specs = {item.output: item for item in GENERATOR_SPECS if item.domain == "tautomer"}
    assert set(specs) == {"tautomer.jsonl", "tautomer_catalog.jsonl"}
    assert specs["tautomer.jsonl"].generator_dependencies == (
        "_tautomer_oracle.py",
        "tautomer_profile.json",
    )
    assert specs["tautomer.jsonl"].deterministic_shards == 16
    assert specs["tautomer_catalog.jsonl"].profiles == frozenset(
        {"tautomer_focused"}
    )
    assert {
        "tautomer_focused",
        "tautomer_pcs_1k",
        "tautomer_pcs_100k",
        "smiles_5000",
    } <= set(PROFILE_INPUTS)
    assert {"tautomer.jsonl", "tautomer_catalog.jsonl"} <= set(SCHEMAS)


def test_profile_pins_sources_branches_and_exact_comparison_fields() -> None:
    profile = oracle.load_profile()
    assert profile["source_revision"] == "351f8f378f8ad6bbd517980c38896e66bf907af8c"
    assert profile["comparison_fields"] == [
        "parse",
        "enumeration_outcome",
        "ordered_smiles",
        "molecule_states",
        "status",
        "modified_atoms",
        "modified_bonds",
        "scores",
        "canonical_smiles",
        "canonical_state",
    ]
    for source in profile["source_paths"]:
        assert (REPO_ROOT / "third_party/rdkit" / source).is_file()
    names = [branch["name"] for branch in profile["branches"]]
    assert names == [
        "default",
        "v1",
        "max_tautomers_1",
        "max_transforms_1",
        "retain_sp3_stereo",
        "retain_bond_stereo",
        "retain_isotopic_hydrogens",
        "no_reassign_stereo",
    ]
    assert profile["corpus_branches"] == ["default", "v1"]
    assert profile["chembl_branches"] == names


def test_catalog_generator_reads_every_pinned_tuple_with_provenance() -> None:
    records = {
        name: catalog_generator.parse_catalog(name, path)
        for name, path in catalog_generator.CATALOGS.items()
    }
    assert {name: len(values) for name, values in records.items()} == {
        "current": 37,
        "v1": 36,
    }
    assert records["current"][0] == {
        "schema_version": 1,
        "catalog": "current",
        "index": 0,
        "source_path": "Code/GraphMol/MolStandardize/TautomerCatalog/tautomerTransforms.in",
        "source_line": 18,
        "name": "1,3 (thio)keto/enol f",
        "smarts": "[CX4!H0R{0-2}]-[C;z{1-2}]=[O,S,Se,Te;X1]",
        "bonds": "",
        "charges": "",
    }
    assert all(record["index"] == index for values in records.values() for index, record in enumerate(values))


def test_oracle_records_are_deterministic_complete_and_schema_valid(tmp_path: Path) -> None:
    oracle.assert_rdkit_version()
    profile = oracle.load_profile()
    inputs = [
        {
            "row": 0,
            "case_id": "keto-enol",
            "smiles": "CC(C)=O",
            "sanitize": True,
            "remove_hs": True,
            "source": "unit",
        },
        {
            "row": 1,
            "case_id": "invalid",
            "smiles": "not-a-smiles",
            "sanitize": True,
            "remove_hs": True,
            "source": "unit",
        },
    ]
    first = [oracle.build_record(record, profile, ["default", "v1"]) for record in inputs]
    second = [oracle.build_record(record, profile, ["default", "v1"]) for record in inputs]
    assert first == second
    assert first[0]["parse"] == {"ok": True, "error": None}
    assert list(first[0]["branches"]) == ["default", "v1"]
    assert "input_tautomer" not in first[0]
    assert "atom_order_permutation" not in first[0]
    for branch in first[0]["branches"].values():
        assert branch["ok"]
        assert branch["ordered_smiles"] == sorted(branch["ordered_smiles"])
        assert len(branch["molecule_states"]) == len(branch["scores"])
        assert branch["canonical_smiles"] == branch["canonical_state"]["isomeric_smiles"]
    assert not first[1]["parse"]["ok"]
    assert first[1]["branches"] == {}

    output = tmp_path / "tautomer.jsonl"
    output.write_text(
        "".join(json.dumps(record, sort_keys=True) + "\n" for record in first),
        encoding="utf-8",
    )
    _, count = validate_jsonl_output(output, output.name)
    assert count == 2


def test_pcs_adapter_preserves_both_source_columns_and_row_identity() -> None:
    path = PROFILE_INPUTS["tautomer_pcs_1k"]
    assert path == REPO_ROOT / "testdata/tautomer/corpus/rdkit/1kPCS_tautomer.csv.gz"
    records = oracle.iter_input_records(path)
    first = next(records)
    assert first["row"] == 0
    assert first["smiles"] == "CCOC(=O)c1cc(-c2ccccn2)c(=O)n2ccccc12"
    assert first["expected_canonical_smiles"] == first["smiles"]
    assert first["source"] == str(path)


def test_pcs_oracle_records_paired_input_and_atom_order_observations() -> None:
    oracle.assert_rdkit_version()
    profile = oracle.load_profile()
    source = next(oracle.iter_input_records(PROFILE_INPUTS["tautomer_pcs_1k"]))
    record = oracle.build_record(source, profile, ["default", "v1"])
    assert record["input_tautomer"]["smiles"] == source["expected_canonical_smiles"]
    assert record["atom_order_permutation"]["smiles"]
    for observation in (
        record["input_tautomer"],
        record["atom_order_permutation"],
    ):
        assert observation["parse"] == {"ok": True, "error": None}
        assert list(observation["branches"]) == ["default", "v1"]
        for branch in observation["branches"].values():
            assert branch["ok"]
            assert branch["canonical_smiles"] == branch["canonical_state"]["isomeric_smiles"]
            assert branch["canonical_score"]["total"] == sum(
                branch["canonical_score"][name]
                for name in ("ring", "substructure", "hetero_hydrogen")
            )


def test_parallel_tautomer_generation_preserves_exact_input_order(tmp_path: Path) -> None:
    oracle.assert_rdkit_version()
    profile = oracle.load_profile()
    records = [
        {
            "row": index,
            "case_id": f"parallel-{index}",
            "smiles": smiles,
            "sanitize": True,
            "remove_hs": True,
            "source": "unit",
        }
        for index, smiles in enumerate(("CC(C)=O", "c1ccccc1O", "CCN", "O=CNC"))
    ]
    serial = tmp_path / "serial.jsonl"
    parallel = tmp_path / "parallel.jsonl"
    tautomer_generator.generate_records(
        records, profile, ["default", "v1"], serial, shards=3, jobs=1
    )
    tautomer_generator.generate_records(
        records, profile, ["default", "v1"], parallel, shards=3, jobs=2
    )
    assert parallel.read_bytes() == serial.read_bytes()
