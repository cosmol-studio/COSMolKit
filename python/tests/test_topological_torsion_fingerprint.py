import json
from typing import cast

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
    generator = cosmolkit.get_topological_torsion_generator()
    options = generator.get_options()

    assert options.include_chirality is False
    assert options.torsion_atom_count == 4
    assert options.count_simulation is True
    assert options.count_bounds == [1, 2, 4, 8]
    assert options.fp_size == 2048
    assert options.num_bits_per_feature == 1
    assert options.only_shortest_paths is False

    sparse_count = generator.get_sparse_count_fingerprint(molecule)
    sparse_bit = generator.get_sparse_fingerprint(molecule)
    count = generator.get_count_fingerprint(molecule)
    bit = generator.get_fingerprint(molecule)

    assert sum(sparse_count.nonzero_elements().values()) == 2
    assert len(sparse_bit.on_bits()) == 2
    assert count.size() == 2048
    assert sum(count.nonzero_elements().values()) == 2
    assert bit.n_bits() == 2048
    assert len(bit.on_bits()) == 2


def test_topological_torsion_bulk_forms_preserve_order_and_match_scalar_calls():
    molecules = _molecules()
    generator = cosmolkit.get_topological_torsion_generator(fp_size=512)

    sparse_counts = generator.get_sparse_count_fingerprints(molecules, num_threads=2)
    sparse_bits = generator.get_sparse_fingerprints(molecules, num_threads=2)
    counts = generator.get_count_fingerprints(molecules, num_threads=2)
    bits = generator.get_fingerprints(molecules, num_threads=2)

    assert len(sparse_counts) == len(molecules)
    assert len(sparse_bits) == len(molecules)
    assert len(counts) == len(molecules)
    assert len(bits) == len(molecules)
    for index, molecule in enumerate(molecules):
        assert sparse_counts[index].nonzero_elements() == generator.get_sparse_count_fingerprint(
            molecule
        ).nonzero_elements()
        assert sparse_bits[index].on_bits() == generator.get_sparse_fingerprint(
            molecule
        ).on_bits()
        assert counts[index].nonzero_elements() == generator.get_count_fingerprint(
            molecule
        ).nonzero_elements()
        assert bits[index].on_bits() == generator.get_fingerprint(molecule).on_bits()


def test_topological_torsion_options_are_live_and_mutate_the_generator():
    molecule = cosmolkit.Molecule.from_smiles("CCCCCC")
    generator = cosmolkit.get_topological_torsion_generator()
    options = generator.get_options()

    baseline = generator.get_fingerprint(molecule)
    options.fp_size = 257
    options.count_simulation = False
    options.include_chirality = True
    options.num_bits_per_feature = 2
    options.torsion_atom_count = 3
    options.only_shortest_paths = True
    options.set_count_bounds([1, 3, 5])

    changed = generator.get_fingerprint(molecule)
    assert options.fp_size == 257
    assert options.count_simulation is False
    assert options.include_chirality is True
    assert options.num_bits_per_feature == 2
    assert options.torsion_atom_count == 3
    assert options.only_shortest_paths is True
    assert options.count_bounds == [1, 3, 5]
    assert changed.n_bits() == 257
    assert changed.on_bits()
    assert (baseline.n_bits(), baseline.on_bits()) != (changed.n_bits(), changed.on_bits())


def test_empty_python_lists_are_converted_to_absent_optional_arguments():
    molecule = cosmolkit.Molecule.from_smiles("CCCCO")
    generator = cosmolkit.get_topological_torsion_generator()

    default = generator.get_sparse_count_fingerprint(molecule).nonzero_elements()
    empty = generator.get_sparse_count_fingerprint(
        molecule,
        from_atoms=[],
        ignore_atoms=[],
        custom_atom_invariants=[],
        custom_bond_invariants=[],
    ).nonzero_elements()
    assert empty == default

    legacy_default = cosmolkit.get_topological_torsion_fingerprint(
        molecule
    ).nonzero_elements()
    legacy_empty = cosmolkit.get_topological_torsion_fingerprint(
        molecule,
        from_atoms=[],
        ignore_atoms=[],
        atom_invariants=[],
    ).nonzero_elements()
    assert legacy_empty == legacy_default


def test_custom_atom_invariants_and_selection_route_through_the_shared_core():
    molecule = cosmolkit.Molecule.from_smiles("CCCCC")
    generator = cosmolkit.get_topological_torsion_generator(count_simulation=False)

    default = generator.get_sparse_count_fingerprint(molecule).nonzero_elements()
    custom = generator.get_sparse_count_fingerprint(
        molecule, custom_atom_invariants=[10, 20, 30, 40, 50]
    ).nonzero_elements()
    rooted = generator.get_sparse_count_fingerprint(
        molecule, from_atoms=[0]
    ).nonzero_elements()
    ignored = generator.get_sparse_count_fingerprint(
        molecule, ignore_atoms=[0]
    ).nonzero_elements()

    assert custom != default
    assert sum(rooted.values()) == 1
    assert sum(ignored.values()) == 1


def test_additional_output_allocations_are_populated():
    generator = cosmolkit.get_topological_torsion_generator(count_simulation=False)
    output = cosmolkit.AdditionalOutput()
    output.allocate_atom_to_bits()
    output.allocate_bit_info_map()
    output.allocate_bit_paths()
    output.allocate_atom_counts()
    output.allocate_atoms_per_bit()

    first_molecule = cosmolkit.Molecule.from_smiles("CCCCO")
    first = generator.get_sparse_count_fingerprint(
        first_molecule, additional_output=output
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
    generator = cosmolkit.get_topological_torsion_generator(
        include_chirality=True,
        torsion_atom_count=3,
        count_simulation=False,
        count_bounds=[1, 3, 7],
        fp_size=777,
    )
    payload = generator.to_json()
    assert isinstance(json.loads(payload), dict)

    restored = cosmolkit.TopologicalTorsionFingerprintGenerator.from_json(payload)
    restored_by_function = cosmolkit.topological_torsion_generator_from_json(payload)
    assert restored.get_options().include_chirality is True
    assert restored.get_options().torsion_atom_count == 3
    assert restored.get_options().count_simulation is False
    assert restored.get_options().count_bounds == [1, 3, 7]
    assert restored.get_options().fp_size == 777
    expected = generator.get_count_fingerprint(molecule).nonzero_elements()
    assert restored.get_count_fingerprint(molecule).nonzero_elements() == expected
    assert restored_by_function.get_count_fingerprint(molecule).nonzero_elements() == expected


def test_legacy_wrappers_are_thin_deterministic_adapters_to_modern_generation():
    molecule = cosmolkit.Molecule.from_smiles("CCCCO")

    unfolded = cosmolkit.get_topological_torsion_fingerprint(molecule)
    modern_sparse = cosmolkit.get_topological_torsion_generator(
        count_simulation=False
    ).get_sparse_count_fingerprint(molecule)
    assert unfolded.nonzero_elements() == modern_sparse.nonzero_elements()

    hashed = cosmolkit.get_hashed_topological_torsion_fingerprint(
        molecule, n_bits=1000
    )
    modern_count = cosmolkit.get_topological_torsion_generator(
        count_simulation=False, fp_size=1000
    ).get_count_fingerprint(molecule)
    assert hashed.nonzero_elements() == modern_count.nonzero_elements()
    assert sorted(hashed.nonzero_elements()) == [24, 288]

    legacy_bit = cosmolkit.get_hashed_topological_torsion_fingerprint_as_bit_vect(
        molecule, n_bits=2048
    )
    assert legacy_bit.n_bits() == 2048
    assert legacy_bit.on_bits()
    assert legacy_bit.on_bits() == cosmolkit.get_hashed_topological_torsion_fingerprint_as_bit_vect(
        molecule, n_bits=2048
    ).on_bits()


def test_atom_pair_and_torsion_helper_utilities_share_the_same_codes():
    parameters = cosmolkit.AtomPairsParameters
    assert parameters.VERSION == "1.1.0"
    assert parameters.CODE_SIZE == (
        parameters.NUM_TYPE_BITS
        + parameters.NUM_PI_BITS
        + parameters.NUM_BRANCH_BITS
    )
    assert parameters.MAX_PATH_LENGTH == (1 << parameters.NUM_PATH_BITS) - 1
    assert 6 in parameters.atom_types

    molecule = cosmolkit.Molecule.from_smiles("CCCC")
    atom_code = cosmolkit.get_atom_code(molecule, 0)
    atom_explanation = cast(tuple[str, int, int], cosmolkit.explain_atom_code(atom_code))
    assert atom_explanation[0] == "C"

    score = cosmolkit.py_score_path(molecule, [0, 1, 2, 3], 4)
    explanation = cosmolkit.explain_path_score(score, 4)
    assert len(explanation) == 4
    assert all(entry[0] == "C" for entry in explanation)
    assert cosmolkit.get_topological_torsion_fingerprint_as_ids(molecule) == [score]


def test_python_surface_returns_typed_errors_for_invalid_inputs():
    molecule = cosmolkit.Molecule.from_smiles("CCCCO")
    generator = cosmolkit.get_topological_torsion_generator()

    with pytest.raises(ValueError, match="topological torsion"):
        _ = cosmolkit.get_topological_torsion_generator(torsion_atom_count=8)
    zero_size = cosmolkit.get_topological_torsion_generator(fp_size=0)
    with pytest.raises(ValueError, match="fingerprint size"):
        _ = zero_size.get_fingerprint(molecule)
    with pytest.raises(ValueError, match="atom invariants"):
        _ = generator.get_count_fingerprint(molecule, custom_atom_invariants=[1])
    with pytest.raises(ValueError):
        _ = cosmolkit.topological_torsion_generator_from_json("not json")
    with pytest.raises(IndexError, match="out of range"):
        _ = cosmolkit.py_score_path(molecule, [0, 1, 2, 99], 4)
    with pytest.raises(ValueError, match="size must be greater than zero"):
        _ = cosmolkit.py_score_path(molecule, [], 0)


def test_live_invalid_options_never_panic_or_return_out_of_range_vectors():
    molecule = cosmolkit.Molecule.from_smiles("CCCC")
    generator = cosmolkit.get_topological_torsion_generator()
    options = generator.get_options()

    options.set_count_bounds([])
    assert generator.get_sparse_count_fingerprint(molecule).nonzero_elements()
    with pytest.raises(ValueError, match="Count bounds are empty"):
        _ = generator.get_sparse_fingerprint(molecule)
    with pytest.raises(ValueError, match="Count bounds are empty"):
        _ = generator.get_fingerprint(molecule)

    options.count_simulation = False
    options.fp_size = 0
    assert generator.get_sparse_count_fingerprint(molecule).nonzero_elements()
    assert generator.get_sparse_fingerprint(molecule).on_bits()
    with pytest.raises(ValueError, match="fingerprint size"):
        _ = generator.get_count_fingerprint(molecule)
    with pytest.raises(ValueError, match="fingerprint size"):
        _ = generator.get_fingerprint(molecule)

    options.fp_size = 2048
    options.num_bits_per_feature = 0
    zero_width = generator.get_count_fingerprint(molecule).nonzero_elements()
    options.num_bits_per_feature = 1
    assert generator.get_count_fingerprint(molecule).nonzero_elements() == zero_width


def test_repeated_python_calls_are_deterministic_over_the_rust_core():
    molecule = cosmolkit.Molecule.from_smiles("CC[C@H](F)Cl")
    generator = cosmolkit.get_topological_torsion_generator(
        include_chirality=True, fp_size=4096
    )

    first = generator.get_count_fingerprint(molecule).nonzero_elements()
    second = generator.get_count_fingerprint(molecule).nonzero_elements()
    third = cosmolkit.get_topological_torsion_generator(
        include_chirality=True, fp_size=4096
    ).get_count_fingerprint(molecule).nonzero_elements()
    assert first == second == third
