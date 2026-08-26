from concurrent.futures import ThreadPoolExecutor
from typing import TypeVar

import cosmolkit
import pytest


SMILES = ["C", "CCO", "c1ccccc1O", "C[C@H](O)F"]
T = TypeVar("T")


def present_values(values: list[T | None]) -> list[T]:
    assert all(value is not None for value in values)
    return [value for value in values if value is not None]


def fingerprint_key(value) -> tuple[int, list[int]]:
    return value.n_bits(), value.on_bits()


def every_fingerprint_family_key(molecule, torsion_generator):
    return {
        "morgan": fingerprint_key(molecule.fingerprint_morgan()),
        "maccs": fingerprint_key(molecule.maccs_fingerprint()),
        "avalon": fingerprint_key(molecule.avalon_fingerprint()),
        "rdk": fingerprint_key(molecule.topological_fingerprint()),
        "atom_pair": fingerprint_key(molecule.fingerprint_atom_pair()),
        "topological_torsion": fingerprint_key(
            torsion_generator.get_fingerprint(molecule)
        ),
        "layered": fingerprint_key(molecule.fingerprint_layered()),
    }


def test_layered_scalar_defaults_types_counts_masks_roots_and_immutability():
    molecule = cosmolkit.Molecule.from_smiles("CCO")
    before = molecule.to_smiles()

    default = molecule.fingerprint_layered()
    assert default.n_bits() == 2048
    assert default.on_bits() == [
        92,
        360,
        596,
        610,
        611,
        674,
        867,
        1044,
        1111,
        1783,
        1784,
    ]
    assert molecule.to_smiles() == before

    topology_mask = molecule.fingerprint_layered(layers=0x01)
    assert topology_mask.on_bits() == [674, 867]
    masked = molecule.fingerprint_layered_with_output(
        atom_counts=[10, 20, 30], set_only_bits=topology_mask
    )
    assert masked.fingerprint().on_bits() == [674, 867]
    assert masked.atom_counts() == [12, 23, 32]

    no_counts = molecule.fingerprint_layered_with_output()
    assert no_counts.atom_counts() is None
    assert no_counts.fingerprint().on_bits() == default.on_bits()

    assert molecule.fingerprint_layered(from_atoms=[]).on_bits() == []
    assert molecule.fingerprint_layered(
        branched_paths=False, from_atoms=[0]
    ).on_bits() == [360, 596, 610, 611, 674, 867, 1044, 1111, 1783, 1784]
    assert molecule.to_smiles() == before


def test_layered_scalar_preserves_source_errors_and_unsigned_layer_bits():
    molecule = cosmolkit.Molecule.from_smiles("CCO")

    assert molecule.fingerprint_layered(
        layers=0xFFFF_FFC0, atom_counts=[5, 6, 7]
    ).on_bits() == []
    with pytest.raises(ValueError, match="minPath==0"):
        molecule.fingerprint_layered(min_path=0)
    with pytest.raises(ValueError, match="maxPath<minPath"):
        molecule.fingerprint_layered(min_path=3, max_path=2)
    with pytest.raises(ValueError, match="bad atomCounts size"):
        molecule.fingerprint_layered(atom_counts=[0, 0])
    with pytest.raises(ValueError, match="bad setOnlyBits size"):
        molecule.fingerprint_layered(
            fp_size=64, set_only_bits=molecule.fingerprint_layered()
        )
    with pytest.raises(ValueError, match="fromAtoms contains atom index out of range"):
        molecule.fingerprint_layered(from_atoms=[3])


def test_layered_batch_matches_ordered_scalar_calls_and_keeps_invalid_positions():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(SMILES)
    values = present_values(
        batch.fingerprint_layered_list(
            layers=0x07,
            fp_size=256,
            branched_paths=False,
            from_atoms=[0],
            n_jobs=4,
            progress_bar=False,
        )
    )
    assert [value.on_bits() for value in values] == [
        cosmolkit.Molecule.from_smiles(smiles)
        .fingerprint_layered(
            layers=0x07,
            fp_size=256,
            branched_paths=False,
            from_atoms=[0],
        )
        .on_bits()
        for smiles in SMILES
    ]

    outputs = present_values(
        batch.fingerprint_layered_with_output_list(
            layers=0x01,
            atom_counts=[0, 0, 0, 0, 0, 0, 0, 0],
            n_jobs=3,
            progress_bar=False,
        )
    )
    assert [value.fingerprint().on_bits() for value in outputs] == [
        cosmolkit.Molecule.from_smiles(smiles)
        .fingerprint_layered(layers=0x01)
        .on_bits()
        for smiles in SMILES
    ]
    assert all(value.atom_counts() is not None for value in outputs)

    kept = cosmolkit.MoleculeBatch.from_smiles_list(
        ["CC", "not-a-smiles", "CCC"], errors="keep"
    )
    kept_values = kept.fingerprint_layered_list(n_jobs=2, progress_bar=False)
    assert [value is not None for value in kept_values] == [True, False, True]
    assert [(error.index(), error.operation()) for error in kept.errors()] == [
        (1, "batch.from_smiles_list")
    ]


def test_layered_batch_thread_counts_repeats_errors_and_concurrency_are_stable():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(SMILES)
    baseline = [
        value.on_bits()
        for value in present_values(
            batch.fingerprint_layered_list(n_jobs=1, progress_bar=False)
        )
    ]
    for n_jobs in (1, 2, 4):
        for _ in range(3):
            current = present_values(
                batch.fingerprint_layered_list(n_jobs=n_jobs, progress_bar=False)
            )
            assert [value.on_bits() for value in current] == baseline

    with pytest.raises(ValueError, match="n_jobs"):
        batch.fingerprint_layered_list(n_jobs=0)
    with pytest.raises(
        cosmolkit.BatchValidationError, match="batch.layered_fingerprint"
    ):
        batch.fingerprint_layered_list(from_atoms=[99], n_jobs=2, progress_bar=False)

    def layered_call() -> list[list[int]]:
        return [
            value.on_bits()
            for value in present_values(
                batch.fingerprint_layered_list(n_jobs=4, progress_bar=False)
            )
        ]

    def morgan_call() -> list[list[int]]:
        return [
            value.on_bits()
            for value in present_values(
                batch.fingerprint_morgan_list(n_jobs=3, progress_bar=False)
            )
        ]

    with ThreadPoolExecutor(max_workers=4) as executor:
        layered_futures = [executor.submit(layered_call) for _ in range(3)]
        morgan_future = executor.submit(morgan_call)
        assert [future.result() for future in layered_futures] == [
            baseline,
            baseline,
            baseline,
        ]
        assert len(morgan_future.result()) == len(SMILES)


def test_layered_calls_compose_with_every_fingerprint_family_without_crosstalk():
    molecule = cosmolkit.Molecule.from_smiles("CC[C@H](F)Cl.c1ccncc1O")
    torsion_generator = cosmolkit.get_topological_torsion_generator()
    before_smiles = molecule.to_smiles()

    family_baseline = every_fingerprint_family_key(molecule, torsion_generator)
    topology_baseline = fingerprint_key(
        molecule.fingerprint_layered(
            layers=0x01,
            min_path=1,
            max_path=4,
            fp_size=257,
            atom_counts=[7] * molecule.num_atoms(),
            branched_paths=False,
            from_atoms=[0, 6],
        )
    )
    ring_baseline = fingerprint_key(
        molecule.fingerprint_layered(
            layers=0x38, min_path=2, max_path=6, fp_size=509
        )
    )

    def complete_call():
        families = every_fingerprint_family_key(molecule, torsion_generator)
        topology = fingerprint_key(
            molecule.fingerprint_layered(
                layers=0x01,
                min_path=1,
                max_path=4,
                fp_size=257,
                atom_counts=[7] * molecule.num_atoms(),
                branched_paths=False,
                from_atoms=[0, 6],
            )
        )
        ring = fingerprint_key(
            molecule.fingerprint_layered(
                layers=0x38, min_path=2, max_path=6, fp_size=509
            )
        )
        return families, topology, ring

    expected = family_baseline, topology_baseline, ring_baseline
    for _ in range(4):
        assert complete_call() == expected

    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = [executor.submit(complete_call) for _ in range(12)]
        assert [future.result() for future in futures] == [expected] * 12

    assert every_fingerprint_family_key(molecule, torsion_generator) == family_baseline
    assert molecule.to_smiles() == before_smiles
