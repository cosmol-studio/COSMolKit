"""Exercise canonical parameters across the Python/Rust transport boundary."""
import cosmolkit as ck
import pytest


def test_hydrogen_parameters_and_in_place_result():
    molecule = ck.Molecule.from_smiles("CO")
    expanded = molecule.with_hydrogens_with_params(ck.AddHsParams(only_on_atoms=[1]))
    assert molecule.num_atoms() == 2
    assert expanded.num_atoms() == 3
    assert molecule.add_hydrogens_with_params_(ck.AddHsParams(only_on_atoms=[1])) is None
    assert molecule.num_atoms() == 3
    assert molecule.without_hydrogens_with_params(ck.RemoveHsParams()).num_atoms() == 2


def test_matrix_parameters_reach_the_owner():
    molecule = ck.Molecule.from_smiles("C=C")
    ordinary = molecule.distance_matrix()
    weighted = molecule.distance_matrix_with_params(ck.DistanceMatrixParams(use_bond_order=True))
    assert ordinary.dimension() == 2
    assert ordinary.get(0, 1) == 1.0
    assert weighted.get(0, 1) == 0.5
    assert weighted.values() == [0.0, 0.5, 0.5, 0.0]
    assert weighted.get(2, 0) is None


def test_fragment_selection_and_seed_transport():
    molecule = ck.Molecule.from_smiles("CCO")
    assert molecule.to_fragment_smiles([0, 1]) == "CC"
    assert molecule.to_fragment_smiles_with_params(ck.FragmentSmilesWriteParams([1, 2])) == "CO"
    assert molecule.to_cx_smiles_with_params(ck.CxSmilesWriteParams(fields=ck.CxSmilesFields.NONE)) == "CCO"
    assert molecule.to_random_smiles(5, 42) == molecule.to_random_smiles_with_params(5, 42, ck.RandomSmilesWriteParams())


def test_failed_coordinate_change_is_atomic():
    molecule = ck.Molecule.from_smiles("CC")
    editor = molecule.to_builder()
    editor.add_3d_conformer([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    molecule = editor.build()
    before = molecule.distance_matrix_3d().values()
    with pytest.raises(ck.OperationError):
        molecule.set_atom_position_(99, [1.0, 2.0, 3.0])
    assert molecule.distance_matrix_3d().values() == before
    moved = molecule.with_atom_position_with_params(1, [3.0, 0.0, 0.0], ck.AtomPositionParams(conformer_id=0))
    assert molecule.distance_matrix_3d().get(0, 1) == 2.0
    assert moved.distance_matrix_3d().get(0, 1) == 3.0


def test_sanitize_flags_preserve_invalid_bits_error():
    flags = ck.SanitizeOperations.PROPERTIES | ck.SanitizeOperations.KEKULIZE
    assert flags.contains(ck.SanitizeOperations.PROPERTIES)
    assert not flags.is_empty()
    with pytest.raises(ck.SanitizeError) as error:
        ck.SanitizeOperations.from_bits(1 << 31)
    assert error.value.domain == "sanitize"
    assert error.value.kind == "InvalidOperations"
    assert error.value.bits == 1 << 31
    assert error.value.unknown_bits == 1 << 31
    assert error.value.__cause__ is None


def test_explicit_valence_query_and_in_place_assignment():
    molecule = ck.Molecule.from_smiles("CC")
    assert molecule.assign_valence_with_params_(ck.ValenceParams(strict=True)) is None
    assert not molecule.has_valence_violation(0)
    with pytest.raises(ck.ValenceError):
        molecule.has_valence_violation(100)


def test_missing_matrix_conformer_preserves_requested_id():
    molecule = ck.Molecule.from_smiles("CC")
    with pytest.raises(ck.MatrixError) as error:
        molecule.distance_matrix_3d_with_params(ck.DistanceMatrix3dParams(conformer_id=99))
    assert error.value.kind == "ConformerNotFound"
    assert error.value.conformer_id == 99


@pytest.mark.parametrize(
    "smiles, stage, cause_type, kind, payload",
    [
        ("CN(C)(C)C", ck.SanitizeStage.Properties, ck.ValenceError, "InvalidValence",
         {"atom": 1, "atomic_number": 7, "formal_charge": 0, "phase": "Explicit", "calculated": 4}),
        ("c1cccc1", ck.SanitizeStage.Kekulize, ck.KekulizeError, "NotKekulizable",
         {"problem_atoms": [0, 1, 2, 3, 4]}),
    ],
)
def test_chemistry_problem_preserves_typed_stage_and_domain_cause(smiles, stage, cause_type, kind, payload):
    molecule = ck.Molecule.from_smiles_with_params(
        smiles, ck.SmilesParseParams(sanitize=False, remove_hydrogens=False)
    )
    report = molecule.detect_chemistry_problems()
    assert len(report.problems) == 1
    problem = report.problems[0]
    assert problem.operation == stage
    assert isinstance(problem.error, ck.ChemistryProblemError)
    cause = problem.error.__cause__
    assert isinstance(cause, cause_type)
    assert cause.kind == kind
    for field, expected in payload.items():
        assert getattr(cause, field) == expected
