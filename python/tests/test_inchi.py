from pathlib import Path

import pytest

import cosmolkit


METHANE_INCHI = "InChI=1S/CH4/h1H4"
METHANE_INCHI_KEY = "VNWKTOKETHGBQD-UHFFFAOYSA-N"
STEREO_ISOTOPE_INCHI = "InChI=1S/CHBrClF/c2-1(3)4/t1-/m0/s1/i1+1"
REPO_ROOT = Path(__file__).resolve().parents[2]


def test_inchi_python_four_scalar_apis_match_exact_methane_results() -> None:
    source = cosmolkit.Molecule.from_smiles("C")
    before = source.to_smiles()

    inchi = cosmolkit.Chem.MolToInchi(source)
    assert inchi == METHANE_INCHI
    assert cosmolkit.Chem.MolToInchiKey(source) == METHANE_INCHI_KEY
    assert cosmolkit.InchiToInchiKey(inchi) == METHANE_INCHI_KEY

    parsed = cosmolkit.Chem.MolFromInchi(inchi, False, False)
    assert parsed is not None
    assert parsed.num_atoms() == 1
    assert parsed.num_bonds() == 0
    assert parsed.atoms()[0].atomic_num() == 6
    assert parsed.atoms()[0].explicit_hydrogens() == 4

    assert source.to_smiles() == before


def test_inchi_python_preserves_source_defined_stereo_and_isotope() -> None:
    parsed = cosmolkit.Chem.MolFromInchi(STEREO_ISOTOPE_INCHI, False, False)
    assert parsed is not None
    assert parsed.num_atoms() == 4
    assert parsed.num_bonds() == 3
    carbon = parsed.atoms()[0]
    assert carbon.atomic_num() == 6
    assert carbon.isotope() == 13
    assert carbon.chiral_tag_name() in {
        "CHI_TETRAHEDRAL_CW",
        "CHI_TETRAHEDRAL_CCW",
    }


def test_inchi_python_exposes_source_diagnostic_as_structured_warning() -> None:
    with pytest.warns(cosmolkit.InchiDiagnosticWarning) as records:
        key = cosmolkit.InchiToInchiKey("")

    assert key == ""
    assert len(records) == 1
    diagnostic = records[0].message
    assert isinstance(diagnostic, cosmolkit.InchiDiagnosticWarning)
    assert diagnostic.level == "error"
    assert diagnostic.message == "Invalid InChI prefix in generating InChI Key\n"


def test_inchi_python_rejects_unsupported_molecule_state_structurally() -> None:
    fixture = (
        REPO_ROOT
        / "testdata/rdkit_builtin/fixtures/Code/GraphMol/FileParsers/Issue3432136_1.mol"
    )
    molecule = cosmolkit.Molecule.read_mol_from_str(
        fixture.read_text(), sanitize=False, remove_hs=False
    )

    with pytest.raises(cosmolkit.InchiUnsupportedStateError) as captured:
        cosmolkit.Chem.MolToInchi(molecule)

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
