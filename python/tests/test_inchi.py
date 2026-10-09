from pathlib import Path

import pytest

import cosmolkit


METHANE_INCHI = "InChI=1S/CH4/h1H4"
METHANE_INCHI_KEY = "VNWKTOKETHGBQD-UHFFFAOYSA-N"
STEREO_ISOTOPE_INCHI = "InChI=1S/CHBrClF/c2-1(3)4/t1-/m0/s1/i1+1"
REPO_ROOT = Path(__file__).resolve().parents[2]


def test_inchi_python_four_entry_points_match_exact_methane_results() -> None:
    source = cosmolkit.Molecule.from_smiles("C")
    before = source.to_smiles()

    inchi = source.to_inchi()
    assert inchi == METHANE_INCHI
    assert source.to_inchi_key() == METHANE_INCHI_KEY
    assert cosmolkit.inchi_to_key(inchi) == METHANE_INCHI_KEY

    parsed = cosmolkit.Molecule.from_inchi(inchi, sanitize=False, remove_hs=False)
    assert parsed is not None
    assert parsed.num_atoms() == 1
    assert parsed.num_bonds() == 0
    assert parsed.atoms()[0].atomic_number() == 6
    assert parsed.atoms()[0].explicit_hydrogens() == 4

    assert source.to_smiles() == before


def test_inchi_python_matches_source_stereo_cleanup_for_isotopic_center() -> None:
    parsed = cosmolkit.Molecule.from_inchi(
        STEREO_ISOTOPE_INCHI, sanitize=False, remove_hs=False
    )
    assert parsed is not None
    assert parsed.num_atoms() == 4
    assert parsed.num_bonds() == 3
    carbon = parsed.atoms()[0]
    assert carbon.atomic_number() == 6
    assert carbon.isotope() == 13
    assert carbon.chiral_tag_name() == "CHI_UNSPECIFIED"


def test_inchi_python_preserves_relative_and_racemic_stereo_options() -> None:
    source = cosmolkit.Molecule.from_smiles("F[C@H](Cl)Br")

    assert source.to_inchi(options="-SRel") == (
        "InChI=1/CHBrClF/c2-1(3)4/h1H/t1-/s2"
    )
    assert source.to_inchi(options="-SRac") == (
        "InChI=1/CHBrClF/c2-1(3)4/h1H/t1-/s3"
    )


def test_inchi_python_preserves_cationic_aromatic_nitrogen_charge() -> None:
    source = cosmolkit.Molecule.from_smiles("C[n+]1ccccc1")
    inchi = source.to_inchi()

    parsed = cosmolkit.Molecule.from_inchi(inchi, sanitize=False, remove_hs=False)
    assert parsed is not None
    nitrogen = next(atom for atom in parsed.atoms() if atom.atomic_number() == 7)
    assert nitrogen.formal_charge() == 1


def test_from_inchi_covers_all_sanitize_remove_hs_valence_branches() -> None:
    inchi = "InChI=1S/HNO3/c2-1(3)4/h(H,2,3,4)"

    for sanitize in (False, True):
        for remove_hs in (False, True):
            parsed = cosmolkit.Molecule.from_inchi(
                inchi, sanitize=sanitize, remove_hs=remove_hs
            )
            assert parsed is not None
            nitrogen = next(atom for atom in parsed.atoms() if atom.atomic_number() == 7)
            cached_nitrogen = parsed.atom_metadata(recalculate=False)[nitrogen.id()]
            expected_valence = 4 if sanitize else 5
            assert cached_nitrogen.explicit_valence() == expected_valence
            assert cached_nitrogen.total_valence() == expected_valence


def test_inchi_python_invalid_key_input_has_canonical_structured_error() -> None:
    # The canonical Rust facade returns Result, projected as InchiError;
    # rejected input is no longer transported as warning + None.
    with pytest.raises(cosmolkit.InchiError) as captured:
        cosmolkit.inchi_to_key("")
    error = captured.value
    assert error.domain == "inchi"
    assert error.operation == "inchi_to_key"
    assert error.kind == cosmolkit.InchiErrorKind.InvalidInput
    assert error.detail == "the InChI engine returned no identifier"


def test_inchi_python_sanitize_rejection_has_canonical_structured_error() -> None:
    inchi = (
        "InChI=1S/C8H16O6S2/c9-5-8(14-16(11,12)13)7(10)6-15-3-1-2-4-15/"
        "h7-10H,1-6H2/t7-,8+/m0/s1"
    )

    with pytest.raises(cosmolkit.InchiError) as captured:
        cosmolkit.Molecule.from_inchi(inchi)
    assert captured.value.domain == "inchi"
    assert captured.value.kind == cosmolkit.InchiErrorKind.Toolkit
    assert captured.value.operation == "mol_from_inchi"
    assert "hydrogen-removal sanitize failed" in captured.value.detail
    assert "Explicit valence for atom # 15 S, 7" in captured.value.detail


def test_inchi_python_supported_substance_group_input_matches_rdkit() -> None:
    fixture = (
        REPO_ROOT
        / "testdata/rdkit_builtin/fixtures/Code/GraphMol/FileParsers/Issue3432136_1.mol"
    )
    molecule = cosmolkit.Molecule.from_mol(
        fixture.read_text(), sanitize=False, remove_hs=False
    )

    # Frozen from pyproject.toml's RDKit 2026.3.1 for this exact MOL fixture,
    # with sanitize=False/removeHs=False. This input is no longer unsupported.
    assert molecule.to_inchi() == "InChI=1S/C5H12/c1-4-5(2)3/h5H,4H2,1-3H3"


def test_inchi_python_error_kind_vocabulary_includes_allocation_failure() -> None:
    # The facade uses one InchiError plus typed kinds, not Python exception
    # subclasses. Constructing a fake exception never exercised allocation.
    names = {"AllocationFailed", "UnsupportedState", "InvalidInput", "InvalidSourceOutput", "SanitizeFailed", "Toolkit", "SourcePort"}
    assert {name for name in vars(cosmolkit.InchiErrorKind) if name[0].isupper()} == names
    kinds = [getattr(cosmolkit.InchiErrorKind, name) for name in sorted(names)]
    assert all(isinstance(kind, cosmolkit.InchiErrorKind) for kind in kinds)
    assert all(left != right for i, left in enumerate(kinds) for right in kinds[i + 1:])


def test_inchi_python_surface_uses_molecule_methods_and_project_naming() -> None:
    molecule = cosmolkit.Molecule.from_smiles("C")

    assert callable(molecule.to_inchi)
    assert callable(molecule.to_inchi_key)
    assert callable(cosmolkit.Molecule.from_inchi)
    assert callable(cosmolkit.inchi_to_key)
    assert not hasattr(cosmolkit, "Chem")
    assert not hasattr(cosmolkit, "InchiToInchiKey")
