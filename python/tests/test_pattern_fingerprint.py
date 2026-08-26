from concurrent.futures import ThreadPoolExecutor

import cosmolkit
import pytest


SMILES = ["C", "CC", "CCO", "c1ccccc1O", "C[C@H](O)F"]


def present_bits(values: list[cosmolkit.Fingerprint | None]) -> list[list[int]]:
    assert all(value is not None for value in values)
    return [value.on_bits() for value in values if value is not None]


def test_pattern_scalar_defaults_options_type_and_exact_bits():
    molecule = cosmolkit.Molecule.from_smiles("CC")

    default = molecule.pattern_fingerprint()
    explicit = molecule.pattern_fingerprint(n_bits=2048, tautomeric=False)
    tautomeric = molecule.pattern_fingerprint(n_bits=2048, tautomeric=True)

    assert isinstance(default, cosmolkit.Fingerprint)
    assert default.n_bits() == 2048
    assert default.on_bits() == [429, 778, 1022, 1061, 1236, 1295]
    assert explicit.on_bits() == default.on_bits()
    assert tautomeric.on_bits() == [429, 776, 778, 1022, 1061, 1236, 1295]


def test_pattern_scalar_preserves_input_and_reports_argument_errors():
    molecule = cosmolkit.Molecule.from_smiles("c1ccccc1O")
    before = molecule.to_smiles()

    first = molecule.pattern_fingerprint(n_bits=127, tautomeric=True)
    second = molecule.pattern_fingerprint(n_bits=127, tautomeric=True)

    assert first.on_bits() == second.on_bits()
    assert molecule.to_smiles() == before
    with pytest.raises(ValueError, match="fingerprint requires n_bits > 0"):
        molecule.pattern_fingerprint(n_bits=0)
    with pytest.raises(TypeError):
        molecule.pattern_fingerprint(n_bits="2048")  # pyright: ignore[reportArgumentType]
    with pytest.raises(TypeError):
        molecule.pattern_fingerprint(tautomeric="yes")  # pyright: ignore[reportArgumentType]


def test_pattern_batch_matches_scalar_order_and_keeps_invalid_positions():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(
        ["CC", "not-a-smiles", "c1ccccc1O"], errors="keep"
    )
    values = batch.pattern_fingerprint_list(
        n_bits=257, tautomeric=True, n_jobs=2, progress_bar=False
    )

    assert [value is not None for value in values] == [True, False, True]
    assert [value.on_bits() if value is not None else None for value in values] == [
        cosmolkit.Molecule.from_smiles("CC")
        .pattern_fingerprint(n_bits=257, tautomeric=True)
        .on_bits(),
        None,
        cosmolkit.Molecule.from_smiles("c1ccccc1O")
        .pattern_fingerprint(n_bits=257, tautomeric=True)
        .on_bits(),
    ]
    with pytest.raises(cosmolkit.BatchValidationError, match="batch.pattern_fingerprint"):
        batch.pattern_fingerprint_list(n_bits=0, n_jobs=2, progress_bar=False)
    with pytest.raises(ValueError, match="n_jobs"):
        batch.pattern_fingerprint_list(n_jobs=0)


def test_pattern_batch_is_deterministic_across_workers_and_concurrent_calls():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(SMILES).with_progress_bar(False)
    expected = present_bits(
        batch.pattern_fingerprint_list(n_bits=511, tautomeric=True, n_jobs=1)
    )

    for n_jobs in (1, 2, 4):
        assert present_bits(
            batch.pattern_fingerprint_list(
                n_bits=511, tautomeric=True, n_jobs=n_jobs
            )
        ) == expected

    def pattern_call() -> list[list[int]]:
        return present_bits(
            batch.pattern_fingerprint_list(
                n_bits=511, tautomeric=True, n_jobs=2, progress_bar=False
            )
        )

    with ThreadPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(lambda _: pattern_call(), range(8)))
    assert results == [expected] * 8


def fingerprint_family_snapshot(molecule: cosmolkit.Molecule) -> dict[str, list[int] | str]:
    torsion = cosmolkit.get_topological_torsion_generator()
    return {
        "smiles": molecule.to_smiles(),
        "avalon": molecule.avalon_fingerprint().on_bits(),
        "morgan": molecule.fingerprint_morgan().on_bits(),
        "atom_pair": molecule.fingerprint_atom_pair().on_bits(),
        "pattern": molecule.pattern_fingerprint().on_bits(),
        "tautomeric_pattern": molecule.pattern_fingerprint(
            n_bits=257, tautomeric=True
        ).on_bits(),
        "topological": molecule.topological_fingerprint().on_bits(),
        "maccs": molecule.maccs_fingerprint().on_bits(),
        "topological_torsion": torsion.get_fingerprint(molecule).on_bits(),
    }


def test_pattern_composes_with_every_fingerprint_family_on_one_shared_molecule():
    molecule = cosmolkit.Molecule.from_smiles("CCOc1ccc(C(=O)N[C@@H](C)F)cc1")
    expected = fingerprint_family_snapshot(molecule)

    for _ in range(4):
        assert molecule.pattern_fingerprint(n_bits=257, tautomeric=True).on_bits() == expected[
            "tautomeric_pattern"
        ]
        assert molecule.fingerprint_atom_pair().on_bits() == expected["atom_pair"]
        assert molecule.maccs_fingerprint().on_bits() == expected["maccs"]
        assert molecule.fingerprint_morgan().on_bits() == expected["morgan"]
        assert molecule.topological_fingerprint().on_bits() == expected["topological"]
        assert molecule.avalon_fingerprint().on_bits() == expected["avalon"]
        assert (
            cosmolkit.get_topological_torsion_generator()
            .get_fingerprint(molecule)
            .on_bits()
            == expected["topological_torsion"]
        )
        assert molecule.pattern_fingerprint().on_bits() == expected["pattern"]
        assert fingerprint_family_snapshot(molecule) == expected

    with ThreadPoolExecutor(max_workers=8) as executor:
        results = list(executor.map(lambda _: fingerprint_family_snapshot(molecule), range(24)))
    assert results == [expected] * 24
    assert molecule.to_smiles() == expected["smiles"]
