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
    assert parsed.atoms()[0].atomic_num() == 6
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
    assert carbon.atomic_num() == 6
    assert carbon.isotope() == 13
    assert carbon.chiral_tag_name() == "CHI_UNSPECIFIED"


def test_inchi_python_preserves_relative_and_racemic_stereo_options() -> None:
    source = cosmolkit.Molecule.from_smiles("F[C@H](Cl)Br")

    assert source.to_inchi("-SRel") == (
        "InChI=1/CHBrClF/c2-1(3)4/h1H/t1-/s2"
    )
    assert source.to_inchi("-SRac") == (
        "InChI=1/CHBrClF/c2-1(3)4/h1H/t1-/s3"
    )


def test_inchi_python_preserves_cationic_aromatic_nitrogen_charge() -> None:
    source = cosmolkit.Molecule.from_smiles("C[n+]1ccccc1")
    inchi = source.to_inchi()

    parsed = cosmolkit.Molecule.from_inchi(inchi, sanitize=False, remove_hs=False)
    assert parsed is not None
    nitrogen = next(atom for atom in parsed.atoms() if atom.atomic_num() == 7)
    assert nitrogen.formal_charge() == 1


def test_from_inchi_covers_all_sanitize_remove_hs_valence_branches() -> None:
    inchi = "InChI=1S/HNO3/c2-1(3)4/h(H,2,3,4)"

    for sanitize in (False, True):
        for remove_hs in (False, True):
            parsed = cosmolkit.Molecule.from_inchi(
                inchi, sanitize=sanitize, remove_hs=remove_hs
            )
            assert parsed is not None
            nitrogen = next(atom for atom in parsed.atoms() if atom.atomic_num() == 7)
            expected_valence = 4 if sanitize else 5
            assert nitrogen.explicit_valence() == expected_valence
            assert nitrogen.total_valence() == expected_valence


def test_inchi_python_exposes_source_diagnostic_as_structured_warning() -> None:
    with pytest.warns(cosmolkit.InchiDiagnosticWarning) as records:
        key = cosmolkit.inchi_to_key("")

    assert key is None
    assert len(records) == 1
    diagnostic = records[0].message
    assert isinstance(diagnostic, cosmolkit.InchiDiagnosticWarning)
    assert diagnostic.level == "error"
    assert diagnostic.message == "Invalid InChI prefix in generating InChI Key\n"


def test_inchi_python_returns_none_for_rdkit_mol_sanitize_exception() -> None:
    inchi = (
        "InChI=1S/C8H16O6S2/c9-5-8(14-16(11,12)13)7(10)6-15-3-1-2-4-15/"
        "h7-10H,1-6H2/t7-,8+/m0/s1"
    )

    assert cosmolkit.Molecule.from_inchi(inchi) is None


def test_inchi_python_rejects_unsupported_molecule_state_structurally() -> None:
    fixture = (
        REPO_ROOT
        / "testdata/rdkit_builtin/fixtures/Code/GraphMol/FileParsers/Issue3432136_1.mol"
    )
    molecule = cosmolkit.Molecule.read_mol_from_str(
        fixture.read_text(), sanitize=False, remove_hs=False
    )

    with pytest.raises(cosmolkit.InchiUnsupportedStateError) as captured:
        _ = molecule.to_inchi()

    error = captured.value
    assert isinstance(error, cosmolkit.InchiError)
    assert error.operation == "mol_to_inchi"
    assert error.kind == "unsupported_state"
    assert error.detail == "InChI bridge does not model substance groups"


def test_inchi_python_allocation_error_contract_is_deterministic_and_structured() -> None:
    error = cosmolkit.InchiAllocationError(
        "mol_to_inchi failed (AllocationFailed): AllocationFailed",
        "mol_to_inchi",
        "allocation_failed",
        "AllocationFailed",
    )

    assert isinstance(error, cosmolkit.InchiError)
    assert error.operation == "mol_to_inchi"
    assert error.kind == "allocation_failed"
    assert error.detail == "AllocationFailed"
    assert str(error) == "mol_to_inchi failed (AllocationFailed): AllocationFailed"


def test_inchi_python_surface_uses_molecule_methods_and_project_naming() -> None:
    molecule = cosmolkit.Molecule.from_smiles("C")

    assert callable(molecule.to_inchi)
    assert callable(molecule.to_inchi_key)
    assert callable(cosmolkit.Molecule.from_inchi)
    assert callable(cosmolkit.inchi_to_key)
    assert not hasattr(cosmolkit, "Chem")
    assert not hasattr(cosmolkit, "InchiToInchiKey")
