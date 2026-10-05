"""Pinned source conditions discovered during independent whole-block review."""
import ast
import math
from pathlib import Path

import cosmolkit as ck
import pytest


@pytest.mark.parametrize("operation", ["missing_conformer", "no_match", "best_map", "write"])
def test_wrapper_weight_length_error_precedes_domain_errors(operation):
    builder = ck.Molecule.from_smiles("CC").to_builder()
    builder.add_3d_conformer([[0., 0., 0.], [1., 0., 0.]])
    molecule = builder.build()
    mapping = [ck.AlignmentAtomMap(0, 0), ck.AlignmentAtomMap(1, 1)]
    before = molecule.to_binary()
    expected_type = ck.OperationError if operation == "write" else ck.AlignmentError
    with pytest.raises(expected_type) as error:
        if operation == "best_map":
            molecule.best_alignment_to_with_params(molecule, ck.BestAlignmentParameters(
                probe_conformer_id=999, atom_maps=[mapping], weights=[1.],
            ))
        elif operation == "no_match":
            molecule.alignment_transform_to(ck.Molecule.from_smiles("OO"), ck.AlignmentParameters(weights=[1.]))
        else:
            params = ck.AlignmentParameters(probe_conformer_id=999, atom_map=mapping, weights=[1.])
            if operation == "write":
                molecule.align_to_with_params_(molecule, params)
            else:
                molecule.alignment_transform_to_with_params(molecule, params)
    cause = error.value.__cause__ if operation == "write" else error.value
    assert isinstance(cause, ck.AlignmentError)
    assert cause.kind == "WeightCountMismatch"
    assert (cause.map_len, cause.weight_len) == (2, 1)
    assert molecule.to_binary() == before


@pytest.mark.parametrize("operation", [
    "alignment_transform", "best_alignment", "coordinate_rmsd", "align_to",
    "all_conformers", "align_conformers", "alignment_empty_map", "align_empty_map",
])
def test_empty_optional_sequences_follow_source_wrapper_without_changing_fields(operation):
    builder = ck.Molecule.from_smiles("CC").to_builder()
    for _ in range(2):
        builder.add_3d_conformer([[0., 0., 0.], [1., 0., 0.]])
    molecule = builder.build()
    if operation in ["alignment_transform", "align_to", "alignment_empty_map", "align_empty_map"]:
        params = ck.AlignmentParameters(**({"atom_map": []} if "empty_map" in operation else {"weights": []}))
        if "empty_map" in operation:
            assert params.atom_map == []
        else:
            assert params.weights == []
        if operation in ["alignment_transform", "alignment_empty_map"]:
            assert molecule.alignment_transform_to(molecule, params).rmsd() == 0.
            assert molecule.alignment_transform_to_with_params(molecule, params).rmsd() == 0.
        else:
            assert molecule.with_alignment_to(molecule, params)[1].rmsd() == 0.
            assert molecule.with_alignment_to_with_params(molecule, params)[1].rmsd() == 0.
            assert molecule.align_to_(molecule, params).rmsd() == 0.
            assert molecule.align_to_with_params_(molecule, params).rmsd() == 0.
    elif operation == "best_alignment":
        params = ck.BestAlignmentParameters(weights=[])
        assert params.weights == []
        assert molecule.best_alignment_to(molecule, params).rmsd() == 0.
        assert molecule.best_alignment_to_with_params(molecule, params).rmsd() == 0.
        assert molecule.best_rmsd_to(molecule, params) == 0.
        assert molecule.best_rmsd_to_with_params(molecule, params) == 0.
    elif operation == "coordinate_rmsd":
        params = ck.CoordinateRmsdParameters(weights=[])
        assert params.weights == []
        assert molecule.coordinate_rmsd_to(molecule, params) == 0.
        assert molecule.coordinate_rmsd_to_with_params(molecule, params) == 0.
    elif operation == "all_conformers":
        params = ck.AllConformerRmsdParameters(weights=[])
        assert params.weights == []
        assert [row.rmsd() for row in molecule.all_conformer_best_rmsds(params)] == [0.]
        assert [row.rmsd() for row in molecule.all_conformer_best_rmsds_with_params(params)] == [0.]
    else:
        params = ck.ConformerAlignmentParameters(weights=[])
        assert params.weights == []
        assert molecule.with_aligned_conformers(params)[1].rmsds() == [0.]
        assert molecule.with_aligned_conformers_with_params(params)[1].rmsds() == [0.]
        assert molecule.align_conformers_(params).rmsds() == [0.]
        assert molecule.align_conformers_with_params_(params).rmsds() == [0.]


@pytest.mark.parametrize("count", [1, 2])
def test_empty_atom_indices_select_all_like_pinned_source_wrapper(count):
    builder = ck.Molecule.from_smiles("CC").to_builder()
    builder.add_3d_conformer([[0., 0., 0.], [1., 0., 0.]])
    if count == 2:
        builder.add_3d_conformer([[2., 1., 0.], [3., 1., 0.]])
    source = builder.build()
    before = source.to_binary()
    params = ck.ConformerAlignmentParameters(atom_indices=[])
    aligned, report = source.with_aligned_conformers(params)
    assert report.rmsds() == [0.] * (count - 1)
    expected = [[[0., 0., 0.], [1., 0., 0.]]] * count
    assert [c.coordinates() for c in aligned.conformers_3d()] == expected
    assert source.to_binary() == before
    assert source.align_conformers_with_params_(params).rmsds() == report.rmsds()
    assert [c.coordinates() for c in source.conformers_3d()] == expected
    with pytest.raises(ck.OperationError) as error:
        source.align_conformers_(ck.ConformerAlignmentParameters(atom_indices=[2]))
    assert isinstance(error.value.__cause__, ck.AlignmentError)
    assert error.value.__cause__.kind == "ProbeAtomOutOfRange"


def test_positive_infinite_weight_preserves_pinned_jacobi_nan_branch_rotation():
    builder = ck.Molecule.from_smiles("CC").to_builder()
    builder.add_3d_conformer([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    molecule = builder.build()
    params = ck.AlignmentParameters(
        atom_map=[ck.AlignmentAtomMap(0, 0), ck.AlignmentAtomMap(1, 1)],
        weights=[math.inf, 1.0],
    )
    result = molecule.alignment_transform_to(molecule, params)
    assert math.isnan(result.rmsd())
    matrix = result.transform().matrix()
    # AlignPoints.cpp:195 rotates only when fabs(b)>0, false for NaN.
    assert [row[:3] for row in matrix[:3]] == [
        [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]
    ]
    assert all(math.isnan(matrix[row][3]) for row in range(3))
    assert matrix[3] == [0.0, 0.0, 0.0, 1.0]


def test_registered_alignment_error_has_one_source_typed_standard_stub_declaration():
    path = Path(__file__).parents[1] / "cosmolkit.pyi"
    tree = ast.parse(path.read_text())
    matches = [node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == "AlignmentError"]
    assert len(matches) == 1
    declaration = matches[0]
    assert [ast.unparse(base) for base in declaration.bases] == ["builtins.ValueError"]
    fields = {node.target.id: ast.unparse(node.annotation) for node in declaration.body if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name)}
    required = {"domain": "builtins.str", "kind": "builtins.str", "id": "builtins.int", "index": "builtins.int", "atom_count": "builtins.int", "map_len": "builtins.int", "weight_len": "builtins.int", "message": "builtins.str"}
    assert {key: fields.get(key) for key in required} == required
    exported = next(ast.literal_eval(node.value) for node in tree.body if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id == "__all__" for target in node.targets))
    assert exported.count("AlignmentError") == 1
