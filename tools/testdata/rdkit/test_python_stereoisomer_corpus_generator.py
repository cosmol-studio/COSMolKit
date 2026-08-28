from __future__ import annotations

import hashlib
import json
from pathlib import Path

from rdkit import Chem

import _generate_python_stereoisomer_corpus_golden as generator
from _expected_schema import SCHEMAS, validate_jsonl_output
from generate_all import GENERATOR_SPECS, validate_cached_domain


def serialize_records(records: list[dict[str, object]]) -> bytes:
    return b"".join(
        (
            json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n"
        ).encode()
        for record in records
    )


def test_corpus_schema_generator_and_profiles_are_registered_once() -> None:
    specs = [
        item
        for item in GENERATOR_SPECS
        if item.output == "python_stereoisomer_corpus.jsonl"
    ]
    assert len(specs) == 1
    spec = specs[0]
    assert spec.script == "_generate_python_stereoisomer_corpus_golden.py"
    assert spec.domain == "stereo"
    assert spec.suites == frozenset(
        {"all", "stereo", "python-stereoisomer-corpus"}
    )
    assert spec.profiles == frozenset(
        {"python_stereoisomer_small", "python_stereoisomer_5000"}
    )
    assert spec.output_schema_version == 2
    assert "python_stereoisomer_corpus.jsonl" in SCHEMAS


def test_stereo_generators_select_only_compatible_profile_schemas() -> None:
    stereo_specs = [item for item in GENERATOR_SPECS if item.domain == "stereo"]
    assert all(item.profiles is not None for item in stereo_specs)

    expected_outputs = {
        "python_stereoisomer_focused": {"python_stereoisomer.jsonl"},
        "python_stereoisomer_small": {"python_stereoisomer_corpus.jsonl"},
        "python_stereoisomer_5000": {"python_stereoisomer_corpus.jsonl"},
    }
    for profile, outputs in expected_outputs.items():
        assert {
            item.output
            for item in stereo_specs
            if profile in item.profiles and "stereo" in item.suites
        } == outputs


def test_profile_order_and_options_match_the_bounded_parity_contract() -> None:
    assert generator.POTENTIAL_PROFILES == (
        {"id": "preserve_possible", "clean_it": False, "flag_possible": True},
        {"id": "clean_possible", "clean_it": True, "flag_possible": True},
    )
    assert [profile["id"] for profile in generator.ENUMERATION_PROFILES] == [
        "default_bounded",
        "all_assigned_bounded",
        "non_unique_bounded",
        "seeded_three",
    ]
    assert generator.ENUMERATION_PROFILES == (
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


def test_generation_is_complete_deterministic_and_preserves_input_order(
    tmp_path: Path,
) -> None:
    source = tmp_path / "input.smi"
    source.write_text("# retained comment\nCC\n\nCC(F)=CC(Cl)C\n", encoding="utf-8")
    rows = list(generator.iter_smiles(source))
    assert rows == [(0, "CC"), (1, "CC(F)=CC(Cl)C")]

    first = [generator.build_record(row, smiles) for row, smiles in rows]
    second = [generator.build_record(row, smiles) for row, smiles in rows]
    assert first == second
    assert [record["row"] for record in first] == [0, 1]
    assert [record["smiles"] for record in first] == ["CC", "CC(F)=CC(Cl)C"]
    for record in first:
        assert record["parse_status"] == "ok"
        assert [
            run["profile"]["id"] for run in record["potential_stereo"]
        ] == ["preserve_possible", "clean_possible"]
        assert [run["profile"]["id"] for run in record["enumeration"]] == [
            "default_bounded",
            "all_assigned_bounded",
            "non_unique_bounded",
            "seeded_three",
        ]
        assert all(
            isinstance(run["theoretical_count"], str)
            and run["theoretical_count"].isdigit()
            for run in record["enumeration"]
        )

    first_bytes = serialize_records(first)
    assert first_bytes == serialize_records(second)
    output = tmp_path / "python_stereoisomer_corpus.jsonl"
    output.write_bytes(first_bytes)
    checksum, records_written = validate_jsonl_output(output, output.name)
    assert records_written == len(rows)
    assert checksum == hashlib.sha256(first_bytes).hexdigest()


def test_enumeration_direction_gauge_is_stable_under_component_wide_flip() -> None:
    molecule = Chem.MolFromSmiles("F/C=C/C=C/Cl")
    assert molecule is not None
    Chem.SetDoubleBondNeighborDirections(molecule)
    flipped = Chem.Mol(molecule)
    for bond in flipped.GetBonds():
        if str(bond.GetBondDir()) == "ENDDOWNRIGHT":
            bond.SetBondDir(Chem.BondDir.ENDUPRIGHT)
        elif str(bond.GetBondDir()) == "ENDUPRIGHT":
            bond.SetBondDir(Chem.BondDir.ENDDOWNRIGHT)

    assert generator.molecule_stereo_state(molecule) != generator.molecule_stereo_state(
        flipped
    )
    assert generator.molecule_stereo_state(
        molecule, normalize_enumeration_directions=True
    ) == generator.molecule_stereo_state(
        flipped, normalize_enumeration_directions=True
    )


def test_cached_output_rejects_partial_or_changed_corpus(tmp_path: Path) -> None:
    output_name = "python_stereoisomer_corpus.jsonl"
    output = tmp_path / output_name
    record = generator.build_record(0, "CC")
    output.write_bytes(serialize_records([record]))
    checksum, records = validate_jsonl_output(output, output_name)
    top_identity = {"schema_version": 1, "domain": "stereo"}
    identity = {
        "path": output_name,
        "output_schema_version": 2,
        "generator": [],
        "inputs": [],
        "options": {"profile": "test"},
        "platform": {"system": "test", "machine": "test"},
    }
    manifest = {
        **top_identity,
        "outputs": [{**identity, "sha256": checksum, "records": records}],
    }
    (tmp_path / "manifest.json").write_text(
        json.dumps(manifest), encoding="utf-8"
    )
    assert validate_cached_domain(tmp_path, top_identity, [identity])[0]

    output.write_text("", encoding="utf-8")
    valid, reason = validate_cached_domain(tmp_path, top_identity, [identity])
    assert not valid
    assert "checksum" in reason or "record count" in reason

    output.write_bytes(serialize_records([{**record, "smiles": "CCC"}]))
    valid, reason = validate_cached_domain(tmp_path, top_identity, [identity])
    assert not valid
    assert "checksum" in reason
