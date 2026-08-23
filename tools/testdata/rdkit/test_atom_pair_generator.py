from __future__ import annotations

import hashlib
import json
from pathlib import Path

import _generate_atom_pair_fingerprint_golden as atom_pair_generator
from _expected_schema import SCHEMAS, validate_jsonl_output
from generate_all import (
    GENERATOR_SPECS,
    REPO_ROOT,
    validate_cached_domain,
)


def serialize_records(records: list[dict[str, object]]) -> bytes:
    return b"".join(
        (json.dumps(record, ensure_ascii=True, sort_keys=True) + "\n").encode()
        for record in records
    )


def test_atom_pair_schema_and_generator_identity_are_registered() -> None:
    specs = [item for item in GENERATOR_SPECS if item.output == "atom_pair_fingerprint.jsonl"]
    assert len(specs) == 1
    spec = specs[0]
    assert spec.script == "_generate_atom_pair_fingerprint_golden.py"
    assert spec.domain == "fingerprint"
    assert "atom-pair" in spec.suites
    assert spec.generator_dependencies == ("atom_pair_fingerprint_profile.json",)
    assert "atom_pair_fingerprint.jsonl" in SCHEMAS


def test_profile_pins_source_version_options_and_exact_comparison_fields() -> None:
    profile = atom_pair_generator.load_profile()
    assert profile["schema_version"] == 1
    assert profile["rdkit_version"] == atom_pair_generator.EXPECTED_RDKIT_VERSION
    assert profile["comparison_fields"] == {
        "sparse_count": ["length", "nonzero_elements"],
        "count": ["length", "nonzero_elements"],
        "sparse_bit": ["length", "on_bits"],
        "explicit_bit": ["length", "on_bits"],
        "additional_output": [
            "atom_counts",
            "atom_to_bits",
            "bit_info_map",
            "atoms_per_bit",
        ],
    }
    for source_path in profile["source_paths"]:
        assert (REPO_ROOT / "third_party/rdkit" / source_path).is_file()

    branches = profile["branches"]
    names = [branch["name"] for branch in branches]
    assert names == [
        "default",
        "chiral",
        "count_simulation_false",
        "custom_bounds_collision",
        "distance_2_4",
        "from_first",
        "ignore_first",
        "custom_atom_index_plus_11",
        "num_bits_per_feature_2",
        "additional_output",
    ]
    required = {
        "name",
        "minDistance",
        "maxDistance",
        "includeChirality",
        "use2D",
        "countSimulation",
        "countBounds",
        "fpSize",
        "numBitsPerFeature",
    }
    assert all(required <= set(branch) for branch in branches)
    by_name = {branch["name"]: branch for branch in branches}
    assert by_name["custom_bounds_collision"]["countBounds"] == [1, 3, 5]
    assert by_name["custom_bounds_collision"]["fpSize"] == 257
    assert by_name["distance_2_4"]["minDistance"] == 2
    assert by_name["distance_2_4"]["maxDistance"] == 4
    assert by_name["from_first"]["fromAtoms"] == "first"
    assert by_name["ignore_first"]["ignoreAtoms"] == "first"
    assert (
        by_name["custom_atom_index_plus_11"]["customAtomInvariants"]
        == "index_plus_11"
    )
    assert by_name["num_bits_per_feature_2"]["numBitsPerFeature"] == 2
    assert by_name["additional_output"]["additionalOutput"] is True


def test_generation_is_complete_deterministic_and_preserves_input_order(
    tmp_path: Path,
) -> None:
    atom_pair_generator.assert_rdkit_version()
    profile = atom_pair_generator.load_profile()
    smiles = ["CCCO", "C[C@H](O)F"]
    first = atom_pair_generator.build_records(smiles, profile)
    second = atom_pair_generator.build_records(smiles, profile)
    assert first == second
    assert [record["row"] for record in first] == [0, 1]
    assert [record["smiles"] for record in first] == smiles

    branch_names = [branch["name"] for branch in profile["branches"]]
    for record in first:
        assert list(record["branches"]) == branch_names
        for branch in profile["branches"]:
            generated = record["branches"][branch["name"]]
            assert generated["parameters"] == branch
            if branch["name"] == "num_bits_per_feature_2":
                # RDKit's sparse AtomPair count vector has a fixed 2^23
                # length. The source-derived second feature bit is not
                # folded on this path, so SparseIntVect rejects it while the
                # other three result forms remain valid.
                assert not generated["sparse_count"]["ok"]
                assert generated["sparse_count"]["error"].startswith("IndexError:")
                assert all(
                    generated[name]["ok"]
                    for name in ("count", "sparse_bit", "explicit_bit")
                )
            else:
                assert all(
                    generated[name]["ok"]
                    for name in (
                        "sparse_count",
                        "count",
                        "sparse_bit",
                        "explicit_bit",
                    )
                )

    first_bytes = serialize_records(first)
    second_bytes = serialize_records(second)
    assert first_bytes == second_bytes
    output = tmp_path / "atom_pair_fingerprint.jsonl"
    output.write_bytes(first_bytes)
    checksum, records = validate_jsonl_output(output, output.name)
    assert records == len(smiles)
    assert checksum == hashlib.sha256(first_bytes).hexdigest()


def test_resolved_atom_arguments_are_explicit_and_source_meaningful() -> None:
    mol = atom_pair_generator.Chem.MolFromSmiles("CCCO")
    assert mol is not None
    profile = atom_pair_generator.load_profile()
    by_name = {branch["name"]: branch for branch in profile["branches"]}
    assert atom_pair_generator.resolved_arguments(mol, by_name["default"]) == {}
    assert atom_pair_generator.resolved_arguments(mol, by_name["from_first"]) == {
        "fromAtoms": [0]
    }
    assert atom_pair_generator.resolved_arguments(mol, by_name["ignore_first"]) == {
        "ignoreAtoms": [0]
    }
    assert atom_pair_generator.resolved_arguments(
        mol, by_name["custom_atom_index_plus_11"]
    ) == {"customAtomInvariants": [11, 12, 13, 14]}


def test_cached_domain_rejects_missing_or_incomplete_atom_pair_output(
    tmp_path: Path,
) -> None:
    output_name = "atom_pair_fingerprint.jsonl"
    output = tmp_path / output_name
    record = {
        "row": 0,
        "smiles": "C",
        "rdkit_ok": True,
        "branches": {},
        "error": None,
    }
    output.write_text(json.dumps(record, sort_keys=True) + "\n", encoding="utf-8")
    checksum, records = validate_jsonl_output(output, output_name)
    top_identity = {"schema_version": 1, "domain": "fingerprint"}
    identity = {
        "path": output_name,
        "output_schema_version": 1,
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

    output.unlink()
    valid, reason = validate_cached_domain(tmp_path, top_identity, [identity])
    assert not valid
    assert "cannot validate" in reason
