from pathlib import Path

import cosmolkit
import numpy as np
import pytest


def test_with_2d_coords_returns_new_molecule_without_mutating_input():
    mol = cosmolkit.Molecule.from_smiles("CCO")

    mol_2d = mol.with_2d_coords()

    assert mol is not mol_2d
    assert mol.num_conformers() == 0
    assert not mol.has_2d_coords()
    assert mol_2d.num_conformers() == 0
    assert mol_2d.has_2d_coords()


def test_from_smiles_sanitize_flag_and_molecule_sanitize_are_value_style():
    raw = cosmolkit.Molecule.from_smiles("CN(=O)=O", sanitize=False)

    sanitized = raw.sanitize()

    assert [atom.formal_charge() for atom in raw.atoms()] == [0, 0, 0, 0]
    assert [atom.formal_charge() for atom in sanitized.atoms()] == [0, 1, -1, 0]
    assert raw is not sanitized


def test_sdf_sanitize_false_is_not_silently_ignored():
    sdf = """ethane
     COSMolKit      2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
$$$$
"""
    with pytest.raises(ValueError, match="sanitize=False is not implemented for SDF"):
        cosmolkit.Molecule.read_sdf_record_from_str(sdf, sanitize=False)


def test_sanitize_strict_false_is_not_silently_ignored():
    mol = cosmolkit.Molecule.from_smiles("CCO")

    with pytest.raises(ValueError, match="strict=False sanitization is not implemented"):
        mol.sanitize(strict=False)


def test_molecule_batch_sanitize_flag_and_transform_are_not_noops():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(
        ["CN(=O)=O"],
        sanitize=False,
        errors="raise",
        n_jobs=1,
    )

    with pytest.raises(cosmolkit.BatchValidationError, match="to_smiles"):
        batch.to_smiles_list(canonical=False)

    sanitized = batch.sanitize(errors="raise", n_jobs=1).to_smiles_list()

    assert sanitized == ["C[N+](=O)[O-]"]


def test_structural_array_access_returns_numpy_arrays():
    mol = cosmolkit.Molecule.from_smiles("CCO").with_2d_coords()

    coords = mol.coords_2d()
    bounds = mol.dg_bounds_matrix()

    assert isinstance(coords, np.ndarray)
    assert coords.shape == (3, 3)
    assert isinstance(bounds, np.ndarray)
    assert bounds.shape == (3, 3)


def test_structural_array_access_supports_numpy_operations():
    mol = cosmolkit.Molecule.from_smiles("CCO").with_2d_coords()

    coords = mol.coords_2d()
    bounds = mol.dg_bounds_matrix()

    centered = coords - coords.mean(axis=0)
    assert centered.shape == coords.shape
    assert np.allclose(centered.mean(axis=0), np.zeros(3))
    assert np.asarray(coords) is coords
    assert np.asarray(bounds) is bounds
    assert np.isclose(bounds[0, 1], bounds[0][1])
    assert np.all(np.diag(bounds) == 0.0)


def test_structural_array_access_can_bridge_to_torch_if_installed():
    torch = pytest.importorskip("torch")
    mol = cosmolkit.Molecule.from_smiles("CCO").with_2d_coords()

    tensor = torch.from_numpy(mol.coords_2d())

    assert tuple(tensor.shape) == (3, 3)
    assert str(tensor.dtype) == "torch.float64"


def test_write_png_auto_prepares_drawing_without_mutating_input(tmp_path: Path):
    mol = cosmolkit.Molecule.from_smiles("CCO")
    output = tmp_path / "ethanol.png"

    mol.write_png(str(output))


    assert output.exists()
    assert output.stat().st_size > 0
    assert mol.num_conformers() == 0


def test_read_sdf_from_str_helpers_return_molecules():
    sdf = """ethane
     COSMolKit      2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
$$$$
"""
    one = cosmolkit.Molecule.read_sdf_record_from_str(sdf)
    many = cosmolkit.Molecule.read_sdf_records_from_str(sdf)

    assert len(one) == 2
    assert len(many) == 1
    assert len(many[0]) == 2


def test_read_sdf_coordinate_dim_can_be_forced():
    sdf = """flat
     COSMolKit      2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""
    mol_2d = cosmolkit.Molecule.read_sdf_record_from_str(sdf, coordinate_dim="2d")
    mol_3d = cosmolkit.Molecule.read_sdf_record_from_str(sdf, coordinate_dim="3d")

    assert mol_2d.coords_2d().shape == (1, 3)
    assert np.allclose(mol_3d.coords_3d(), np.array([[0.0, 0.0, 0.0]]))


def test_molecule_batch_processes_in_order_and_preserves_single_molecule_inputs():
    mol = cosmolkit.Molecule.from_smiles("CCO")
    batch = cosmolkit.MoleculeBatch.from_smiles_list(
        ["CCO", "c1ccccc1", "CC(=O)O"], errors="keep", n_jobs=2
    )

    prepared = batch.add_hydrogens(errors="keep", n_jobs=2).compute_2d_coords(
        errors="keep", n_jobs=2
    )

    assert len(batch) == 3
    assert batch.valid_mask() == [True, True, True]
    assert prepared.to_smiles_list() == ["[H]OC([H])([H])C([H])([H])[H]", "[H]c1c([H])c([H])c([H])c([H])c1[H]", "[H]OC(=O)C([H])([H])[H]"]
    assert mol.to_smiles() == "CCO"


def test_molecule_to_smiles_exposes_writer_options():
    benzene = cosmolkit.Molecule.from_smiles("c1ccccc1")
    ethanol = cosmolkit.Molecule.from_smiles("CCO")
    mapped = cosmolkit.Molecule.from_smiles("[CH3:7][OH:2]")

    assert benzene.to_smiles(kekule=True) == "C1=CC=CC=C1"
    assert ethanol.to_smiles(all_bonds_explicit=True) == "C-C-O"
    assert ethanol.to_smiles(all_hs_explicit=True) == "[CH3][CH2][OH]"
    assert ethanol.to_smiles(canonical=False, rooted_at_atom=2) == "OCC"
    assert mapped.to_smiles(canonical=False) == "[CH3:7][OH:2]"
    assert mapped.to_smiles(ignore_atom_map_numbers=True) == "OC"


def test_molecule_batch_parallel_jobs_are_value_style_and_inherited():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(["CCO", "CC"], errors="keep")

    configured = batch.with_parallel_jobs(2)
    prepared = configured.add_hydrogens(errors="keep").compute_2d_coords(errors="keep")

    assert batch.parallel_jobs() is None
    assert configured.parallel_jobs() == 2
    assert prepared.parallel_jobs() == 2
    assert prepared.to_smiles_list(n_jobs=1)[0].startswith("[H]O")
    assert configured.with_parallel_jobs(None).parallel_jobs() is None
    with pytest.raises(ValueError, match="n_jobs must be >= 1"):
        batch.with_parallel_jobs(0)


def test_molecule_batch_to_list_indexing_and_iteration_return_molecules_or_none():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(["CCO", "C1CC", "CC"], errors="keep")

    values = batch.to_list()

    assert [value.to_smiles() if value is not None else None for value in values] == [
        "CCO",
        None,
        "CC",
    ]
    assert batch[0].to_smiles() == "CCO"
    assert batch[-1].to_smiles() == "CC"
    assert batch[1] is None
    assert [value.to_smiles() if value is not None else None for value in batch] == [
        "CCO",
        None,
        "CC",
    ]
    with pytest.raises(IndexError):
        _ = batch[3]


def test_molecule_batch_getitem_supports_slices_masks_and_index_lists():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(
        ["CCO", "C1CC", "CC", "N"],
        errors="keep",
    ).with_parallel_jobs(2)

    tail = batch[2:]
    assert isinstance(tail, cosmolkit.MoleculeBatch)
    assert tail.parallel_jobs() == 2
    assert tail.to_smiles_list() == ["CC", "N"]

    reversed_batch = batch[::-1]
    assert reversed_batch.to_smiles_list() == ["N", "CC", None, "CCO"]

    valid = batch[batch.valid_mask()]
    assert valid.to_smiles_list() == ["CCO", "CC", "N"]

    selected = batch[[3, 0, -2]]
    assert selected.to_smiles_list() == ["N", "CCO", "CC"]

    assert len(batch[[]]) == 0
    with pytest.raises(IndexError, match="boolean mask length"):
        _ = batch[[True, False]]
    with pytest.raises(IndexError):
        _ = batch[[99]]
    with pytest.raises(TypeError):
        _ = batch[True]


def test_molecule_batch_dg_bounds_matrix_list_returns_numpy_arrays_or_none():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(["CCO", "C1CC", "CC"], errors="keep")

    matrices = batch.dg_bounds_matrix_list(n_jobs=2)

    assert isinstance(matrices[0], np.ndarray)
    assert matrices[0].shape == (3, 3)
    assert matrices[1] is None
    assert isinstance(matrices[2], np.ndarray)
    assert matrices[2].shape == (2, 2)


def test_molecule_batch_keeps_errors_and_filters_valid_records(tmp_path: Path):
    batch = cosmolkit.MoleculeBatch.from_smiles_list(["CCO", "C1CC"], errors="keep", n_jobs=2)

    assert len(batch) == 2
    assert batch.valid_mask() == [True, False]
    assert batch.invalid_mask() == [False, True]
    assert batch.valid_count() == 1
    assert batch.invalid_count() == 1
    assert batch.errors()[0].index() == 1
    assert batch.filter_valid().to_smiles_list() == ["CCO"]

    report = batch.compute_2d_coords(errors="skip", n_jobs=2).to_sdf(
        str(tmp_path / "valid.sdf"), errors="skip", n_jobs=2
    )
    assert report.success() == 1
    assert report.failed() == 0
    assert (tmp_path / "valid.sdf").exists()


def test_molecule_batch_exports_use_custom_filenames(tmp_path: Path):
    batch = (
        cosmolkit.MoleculeBatch.from_smiles_list(["CCO", "C1CC", "CC"], errors="keep")
        .with_parallel_jobs(2)
        .compute_2d_coords(errors="keep")
    )

    image_report = batch.to_images(
        str(tmp_path / "images"),
        format="svg",
        errors="skip",
        filenames=["ethanol", "bad.svg", None],
    )
    assert image_report.success() == 2
    assert (tmp_path / "images" / "ethanol.svg").exists()
    assert (tmp_path / "images" / "000002.svg").exists()

    sdf_report = batch.to_sdf_files(
        str(tmp_path / "sdf"),
        format="v2000",
        errors="skip",
        filenames=["ethanol", "bad.sdf", None],
    )
    assert sdf_report.success() == 2
    assert (tmp_path / "sdf" / "ethanol.sdf").exists()
    assert (tmp_path / "sdf" / "000002.sdf").exists()

    with pytest.raises(cosmolkit.BatchValidationError, match="FilenameError"):
        batch.to_images(str(tmp_path / "bad"), format="svg", filenames=["../x", None, None])


def test_molecule_batch_parallel_smiles_writer_options():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(
        ["[*:1]C", "[13CH3:7][C@H](F)Cl"],
        errors="keep",
        n_jobs=2,
    )

    assert batch.to_smiles_list(n_jobs=2) == ["C[*:1]", "F[C@H](Cl)[13CH3:7]"]
    assert batch.to_smiles_list(canonical=False, n_jobs=2) == [
        "[*:1]C",
        "[13CH3:7][C@H](F)Cl",
    ]
    assert batch.to_smiles_list(ignore_atom_map_numbers=True, n_jobs=2) == [
        "C[*]",
        "F[C@H](Cl)[13CH3]",
    ]
    assert batch.to_smiles_list(all_bonds_explicit=True, n_jobs=2)[0] == "C-[*:1]"
    assert batch.to_smiles_list(rooted_at_atom=0, n_jobs=2) == [
        "[*:1]C",
        "[13CH3:7][C@H](F)Cl",
    ]


def test_molecule_batch_raise_aggregates_errors():
    with pytest.raises(cosmolkit.BatchValidationError) as excinfo:
        cosmolkit.MoleculeBatch.from_smiles_list(["C1CC", "N1"], errors="raise", n_jobs=2)

    assert "batch validation failed" in str(excinfo.value)
    errors = excinfo.value.errors()
    assert [error.error_type() for error in errors] == [
        cosmolkit.BatchErrorType.SMILES_PARSE,
        cosmolkit.BatchErrorType.SMILES_PARSE,
    ]
    assert errors[0].error_type_code() == int(cosmolkit.BatchErrorType.SMILES_PARSE)
    assert errors[0].error_type_name() == "SmilesParseError"


def test_batch_errors_expose_intenum_types_and_mode_enum_is_accepted():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(
        ["CCO", "C1CC"],
        errors=cosmolkit.BatchErrorMode.KEEP,
        n_jobs=2,
    )

    errors = batch.errors()
    assert len(errors) == 1
    assert errors[0].error_type() == cosmolkit.BatchErrorType.SMILES_PARSE
    assert errors[0].error_type_name() == "SmilesParseError"
    assert errors[0].as_dict()[3] == ("error_type", "SmilesParseError")
    assert cosmolkit.BATCH_ERROR_TYPE_MAP["SmilesParseError"] == (
        cosmolkit.BatchErrorType.SMILES_PARSE
    )
    assert cosmolkit.BATCH_ERROR_MODE_MAP["keep"] == cosmolkit.BatchErrorMode.KEEP

    filtered = batch.compute_2d_coords(errors=cosmolkit.BatchErrorMode.SKIP)
    assert len(filtered) == 1
