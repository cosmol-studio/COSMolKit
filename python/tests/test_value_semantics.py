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


def test_molecule_batch_raise_aggregates_errors():
    with pytest.raises(cosmolkit.BatchValidationError) as excinfo:
        cosmolkit.MoleculeBatch.from_smiles_list(["C1CC", "N1"], errors="raise", n_jobs=2)

    assert "batch validation failed" in str(excinfo.value)
