import json
import cosmolkit
import pytest


def _molecules():
    return [
        cosmolkit.Molecule.from_smiles("CCCCO"),
        cosmolkit.Molecule.from_smiles("CCCCC"),
        cosmolkit.Molecule.from_smiles("CC(C)CC"),
    ]


def test_topological_torsion_generator_defaults_and_scalar_vector_forms():
    molecule = cosmolkit.Molecule.from_smiles("CCCCO")
    generator = cosmolkit.TopologicalTorsionFingerprintGenerator()
    options = generator.settings()

    assert options.include_chirality is False
    assert options.torsion_atom_count == 4
    assert options.count_simulation is True
    assert options.count_bounds == [1, 2, 4, 8]
    assert options.fp_size == 2048
    assert options.bits_per_feature == 1
    assert options.only_shortest_paths is False

    sparse_count = molecule.fingerprint_topological_torsion_sparse_count_with_generator(generator)
    sparse_bit = molecule.fingerprint_topological_torsion_sparse_with_generator(generator)
    count = molecule.fingerprint_topological_torsion_count_with_generator(generator)
    bit = molecule.fingerprint_topological_torsion_with_generator(generator)

    assert sum(sparse_count.nonzero_elements().values()) == 2
    assert len(sparse_bit.on_bits()) == 2
    assert count.length() == 2048
    assert sum(count.nonzero_elements().values()) == 2
    assert bit.n_bits() == 2048
    assert len(bit.on_bits()) == 2


def test_topological_torsion_bulk_forms_preserve_order_and_match_scalar_calls():
    molecules = _molecules()
    generator = cosmolkit.TopologicalTorsionFingerprintGenerator(params=cosmolkit.TopologicalTorsionParams(fp_size=512))

    sparse_counts = generator.sparse_counts(molecules, num_threads=2)
    sparse_bits = generator.sparse_fingerprints(molecules, num_threads=2)
    counts = generator.counts(molecules, num_threads=2)
    bits = generator.fingerprints(molecules, num_threads=2)

    assert len(sparse_counts) == len(molecules)
    assert len(sparse_bits) == len(molecules)
    assert len(counts) == len(molecules)
    assert len(bits) == len(molecules)
    for index, molecule in enumerate(molecules):
        assert sparse_counts[index].nonzero_elements() == molecule.fingerprint_topological_torsion_sparse_count_with_generator(
            generator
        ).nonzero_elements()
        assert sparse_bits[index].on_bits() == molecule.fingerprint_topological_torsion_sparse_with_generator(
            generator
        ).on_bits()
        assert counts[index].nonzero_elements() == molecule.fingerprint_topological_torsion_count_with_generator(
            generator
        ).nonzero_elements()
        assert bits[index].on_bits() == molecule.fingerprint_topological_torsion_with_generator(generator).on_bits()


def test_topological_torsion_options_are_live_and_mutate_the_generator():
    molecule = cosmolkit.Molecule.from_smiles("CCCCCC")
    generator = cosmolkit.TopologicalTorsionFingerprintGenerator()
    options = generator.settings()

    baseline = molecule.fingerprint_topological_torsion_with_generator(generator)
    options.fp_size = 257
    options.count_simulation = False
    options.include_chirality = True
    options.bits_per_feature = 2
    options.torsion_atom_count = 3
    options.only_shortest_paths = True
    options.set_count_bounds([1, 3, 5])

    changed = molecule.fingerprint_topological_torsion_with_generator(generator)
    assert options.fp_size == 257
    assert options.count_simulation is False
    assert options.include_chirality is True
    assert options.bits_per_feature == 2
    assert options.torsion_atom_count == 3
    assert options.only_shortest_paths is True
    assert options.count_bounds == [1, 3, 5]
    assert changed.n_bits() == 257
    assert changed.on_bits()
    assert (baseline.n_bits(), baseline.on_bits()) != (changed.n_bits(), changed.on_bits())


def test_none_and_empty_atom_selections_remain_distinct():
    molecule = cosmolkit.Molecule.from_smiles("CCCCO")
    generator = cosmolkit.TopologicalTorsionFingerprintGenerator()

    default = molecule.fingerprint_topological_torsion_sparse_count_with_generator(generator).nonzero_elements()
    assert default
    absent = molecule.fingerprint_topological_torsion_sparse_count_with_generator(
        generator, params=cosmolkit.TopologicalTorsionCallParams(from_atoms=None)
    ).nonzero_elements()
    empty = molecule.fingerprint_topological_torsion_sparse_count_with_generator(
        generator, params=cosmolkit.TopologicalTorsionCallParams(from_atoms=[])
    ).nonzero_elements()
    assert absent == default
    assert empty == {}
    unfiltered = molecule.fingerprint_topological_torsion_sparse_count_with_generator(
        generator, params=cosmolkit.TopologicalTorsionCallParams(ignore_atoms=[])
    ).nonzero_elements()
    assert unfiltered == default
    with pytest.raises(cosmolkit.TopologicalTorsionReadError, match="bad atom invariants size"):
        molecule.fingerprint_topological_torsion_sparse_count_with_generator(
            generator, params=cosmolkit.TopologicalTorsionCallParams(custom_atom_invariants=[])
        )

    legacy_default = molecule.fingerprint_topological_torsion_sparse_count_legacy().nonzero_elements()
    legacy_empty = molecule.fingerprint_topological_torsion_sparse_count_legacy_with_params(
        cosmolkit.LegacyTopologicalTorsionParams(from_atoms=[])
    ).nonzero_elements()
    assert legacy_default == default
    assert legacy_empty == {}


def test_custom_atom_invariants_and_selection_route_through_the_shared_core():
    molecule = cosmolkit.Molecule.from_smiles("CCCCC")
    generator = cosmolkit.TopologicalTorsionFingerprintGenerator(params=cosmolkit.TopologicalTorsionParams(count_simulation=False))

    default = molecule.fingerprint_topological_torsion_sparse_count_with_generator(generator).nonzero_elements()
    custom = molecule.fingerprint_topological_torsion_sparse_count_with_generator(generator, params=cosmolkit.TopologicalTorsionCallParams(custom_atom_invariants=[10, 20, 30, 40, 50]
    )).nonzero_elements()
    rooted = molecule.fingerprint_topological_torsion_sparse_count_with_generator(generator, params=cosmolkit.TopologicalTorsionCallParams(from_atoms=[0]
    )).nonzero_elements()
    ignored = molecule.fingerprint_topological_torsion_sparse_count_with_generator(generator, params=cosmolkit.TopologicalTorsionCallParams(ignore_atoms=[0]
    )).nonzero_elements()

    assert custom != default
    assert sum(rooted.values()) == 1
    assert sum(ignored.values()) == 1


def test_additional_output_allocations_are_populated():
    generator = cosmolkit.TopologicalTorsionFingerprintGenerator(params=cosmolkit.TopologicalTorsionParams(count_simulation=False))
    output = cosmolkit.FingerprintAdditionalOutput()
    output.allocate_atom_to_bits()
    output.allocate_bit_info_map()
    output.allocate_bit_paths()
    output.allocate_atom_counts()
    output.allocate_atoms_per_bit()

    first_molecule = cosmolkit.Molecule.from_smiles("CCCCO")
    first = first_molecule.fingerprint_topological_torsion_sparse_count_with_generator(generator, output=output
    )
    atom_to_bits = output.atom_to_bits()
    atom_counts = output.atom_counts()
    bit_paths = output.bit_paths()
    atoms_per_bit = output.atoms_per_bit()
    assert atom_to_bits is not None
    assert atom_counts is not None
    assert bit_paths is not None
    assert atoms_per_bit is not None
    assert len(atom_to_bits) == first_molecule.num_atoms()
    assert len(atom_counts) == first_molecule.num_atoms()
    assert output.bit_info_map() == {}
    assert set(bit_paths) == set(first.nonzero_elements())
    assert atoms_per_bit == bit_paths


def test_generator_json_round_trip_preserves_options_and_fingerprints():
    molecule = cosmolkit.Molecule.from_smiles("CC[C@H](F)Cl")
    generator = cosmolkit.TopologicalTorsionFingerprintGenerator(params=cosmolkit.TopologicalTorsionParams(
        include_chirality=True,
        torsion_atom_count=3,
        count_simulation=False,
        count_bounds=[1, 3, 7],
        fp_size=777,
    ))
    payload = generator.to_json()
    assert isinstance(json.loads(payload), dict)

    restored = cosmolkit.TopologicalTorsionFingerprintGenerator.from_json(payload)
    restored_by_function = cosmolkit.TopologicalTorsionFingerprintGenerator.from_json(payload)
    assert restored.settings().include_chirality is True
    assert restored.settings().torsion_atom_count == 3
    assert restored.settings().count_simulation is False
    assert restored.settings().count_bounds == [1, 3, 7]
    assert restored.settings().fp_size == 777
    expected = molecule.fingerprint_topological_torsion_count_with_generator(generator).nonzero_elements()
    assert molecule.fingerprint_topological_torsion_count_with_generator(restored).nonzero_elements() == expected
    assert molecule.fingerprint_topological_torsion_count_with_generator(restored_by_function).nonzero_elements() == expected


def test_legacy_wrappers_are_thin_deterministic_adapters_to_modern_generation():
    molecule = cosmolkit.Molecule.from_smiles("CCCCO")

    unfolded = molecule.fingerprint_topological_torsion_sparse_count_legacy()
    modern_generator = cosmolkit.TopologicalTorsionFingerprintGenerator(params=cosmolkit.TopologicalTorsionParams(
        count_simulation=False
    ))
    modern_sparse = molecule.fingerprint_topological_torsion_sparse_count_with_generator(modern_generator)
    assert unfolded.nonzero_elements() == modern_sparse.nonzero_elements()

    hashed = molecule.fingerprint_topological_torsion_count_legacy_with_params(
        cosmolkit.LegacyTopologicalTorsionParams(fp_size=1000)
    )
    modern_generator = cosmolkit.TopologicalTorsionFingerprintGenerator(params=cosmolkit.TopologicalTorsionParams(
        count_simulation=False, fp_size=1000
    ))
    modern_count = molecule.fingerprint_topological_torsion_count_with_generator(modern_generator)
    assert hashed.nonzero_elements() == modern_count.nonzero_elements()
    assert sorted(hashed.nonzero_elements()) == [24, 288]

    legacy_bit = molecule.fingerprint_topological_torsion_legacy_with_params(
        cosmolkit.LegacyTopologicalTorsionParams(fp_size=2048)
    )
    assert legacy_bit.n_bits() == 2048
    assert legacy_bit.on_bits()
    assert legacy_bit.on_bits() == molecule.fingerprint_topological_torsion_legacy_with_params(
        cosmolkit.LegacyTopologicalTorsionParams(fp_size=2048)
    ).on_bits()


def test_atom_pair_and_torsion_helper_utilities_share_the_same_codes():
    parameters = cosmolkit.AtomPairsParameters
    assert parameters.version() == "1.1.0"
    assert parameters.code_size() == (
        parameters.num_type_bits()
        + parameters.num_pi_bits()
        + parameters.num_branch_bits()
    )
    assert parameters.max_path_length() == (1 << parameters.num_path_bits()) - 1
    assert 6 in parameters.atom_types()

    molecule = cosmolkit.Molecule.from_smiles("CCCC")
    atom_code = molecule.with_atom_pair_atom_code(0).code
    atom_explanation = cosmolkit.AtomCodeExplanation.from_code(atom_code)
    assert atom_explanation.symbol() == "C"

    score = molecule.topological_torsion_path_score([0, 1, 2, 3], 4)
    explanation = cosmolkit.explain_path_score(score, 4)
    assert len(explanation) == 4
    assert all(entry[0] == "C" for entry in explanation)
    assert molecule.topological_torsion_ids() == [score]


def test_python_surface_returns_typed_errors_for_invalid_inputs():
    molecule = cosmolkit.Molecule.from_smiles("CCCCO")
    generator = cosmolkit.TopologicalTorsionFingerprintGenerator()

    # RDKit's factory accepts 8; construction does not calculate the
    # unfolded size (whose shift would be undefined for this setting).
    eight = cosmolkit.TopologicalTorsionFingerprintGenerator(
        params=cosmolkit.TopologicalTorsionParams(torsion_atom_count=8)
    )
    assert eight.settings().torsion_atom_count == 8
    assert molecule.fingerprint_topological_torsion_count_with_generator(eight).nonzero_elements() == {}
    zero_size = cosmolkit.TopologicalTorsionFingerprintGenerator(params=cosmolkit.TopologicalTorsionParams(fp_size=0))
    with pytest.raises(ValueError, match="fingerprint size"):
        _ = molecule.fingerprint_topological_torsion_with_generator(zero_size)
    with pytest.raises(ValueError, match="atom invariants"):
        _ = molecule.fingerprint_topological_torsion_count_with_generator(generator, params=cosmolkit.TopologicalTorsionCallParams(custom_atom_invariants=[1]))
    with pytest.raises(ValueError):
        _ = cosmolkit.TopologicalTorsionFingerprintGenerator.from_json("not json")
    with pytest.raises(IndexError, match="out of range"):
        _ = molecule.topological_torsion_path_score([0, 1, 2, 99], 4)
    with pytest.raises(cosmolkit.TopologicalTorsionPathScoreError, match="size must be greater than zero") as caught:
        _ = molecule.topological_torsion_path_score([], 0)
    assert caught.value.kind == "ZeroSize"


def test_live_invalid_options_never_panic_or_return_out_of_range_vectors():
    molecule = cosmolkit.Molecule.from_smiles("CCCC")
    generator = cosmolkit.TopologicalTorsionFingerprintGenerator()
    options = generator.settings()

    options.set_count_bounds([])
    assert molecule.fingerprint_topological_torsion_sparse_count_with_generator(generator).nonzero_elements()
    with pytest.raises(cosmolkit.TopologicalTorsionReadError, match="effectiveSize / countBounds.size") as caught:
        _ = molecule.fingerprint_topological_torsion_sparse_with_generator(generator)
    assert caught.value.domain == "fingerprints" and caught.value.kind == "Generator"
    assert isinstance(caught.value.__cause__, ValueError)
    assert "C++-undefined operation" in str(caught.value.__cause__)
    with pytest.raises(ValueError, match="Count bounds are empty"):
        _ = molecule.fingerprint_topological_torsion_with_generator(generator)

    options.count_simulation = False
    options.fp_size = 0
    assert molecule.fingerprint_topological_torsion_sparse_count_with_generator(generator).nonzero_elements()
    assert molecule.fingerprint_topological_torsion_sparse_with_generator(generator).on_bits()
    with pytest.raises(cosmolkit.TopologicalTorsionReadError, match="outside vector length 0") as caught:
        _ = molecule.fingerprint_topological_torsion_count_with_generator(generator)
    assert caught.value.kind == "Generator"
    assert isinstance(caught.value.__cause__, ValueError)
    with pytest.raises(cosmolkit.TopologicalTorsionReadError, match="outside vector length 0") as caught:
        _ = molecule.fingerprint_topological_torsion_with_generator(generator)
    assert caught.value.kind == "Generator"
    assert isinstance(caught.value.__cause__, ValueError)

    options.fp_size = 2048
    options.bits_per_feature = 0
    zero_width = molecule.fingerprint_topological_torsion_count_with_generator(generator).nonzero_elements()
    options.bits_per_feature = 1
    assert molecule.fingerprint_topological_torsion_count_with_generator(generator).nonzero_elements() == zero_width


def test_repeated_python_calls_are_deterministic_over_the_rust_core():
    molecule = cosmolkit.Molecule.from_smiles("CC[C@H](F)Cl")
    generator = cosmolkit.TopologicalTorsionFingerprintGenerator(params=cosmolkit.TopologicalTorsionParams(
        include_chirality=True, fp_size=4096
    ))

    first = molecule.fingerprint_topological_torsion_count_with_generator(generator).nonzero_elements()
    second = molecule.fingerprint_topological_torsion_count_with_generator(generator).nonzero_elements()
    repeated = cosmolkit.TopologicalTorsionFingerprintGenerator(params=cosmolkit.TopologicalTorsionParams(
        include_chirality=True, fp_size=4096
    ))
    third = molecule.fingerprint_topological_torsion_count_with_generator(repeated).nonzero_elements()
    assert first == second == third
