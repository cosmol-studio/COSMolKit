from pathlib import Path

import cosmolkit
import numpy as np
import pytest


def assert_coordinate_rows_match_atoms(mol: cosmolkit.Molecule) -> None:
    atom_count = len(mol)
    if mol.has_2d_coordinates():
        assert mol.coordinates_2d().shape == (atom_count, 3)
    for conformer_index in range(mol.num_conformers()):
        assert mol.coordinates_3d(conformer_index).shape == (atom_count, 3)


def require_valid_molecules(batch: cosmolkit.MoleculeBatch) -> list[cosmolkit.Molecule]:
    molecules = list(batch)
    assert all(mol is not None for mol in molecules)
    return [mol for mol in molecules if mol is not None]


def test_with_2d_coordinates_returns_new_molecule_without_mutating_input():
    mol = cosmolkit.Molecule.from_smiles("CCO")

    assert not hasattr(mol, "with_2d_coords")
    assert not hasattr(mol, "has_2d_coords")
    assert not hasattr(mol, "coords_2d")
    assert not hasattr(mol, "coords_3d")
    mol_2d = mol.with_2d_coordinates()

    assert mol is not mol_2d
    assert mol.num_conformers() == 0
    assert not mol.has_2d_coordinates()
    assert mol_2d.num_conformers() == 0
    assert mol_2d.has_2d_coordinates()


def test_from_smiles_sanitize_flag_and_molecule_sanitize_are_value_style():
    raw = cosmolkit.Molecule.from_smiles("CN(=O)=O", sanitize=False)

    sanitized = raw.sanitize()

    assert [atom.formal_charge() for atom in raw.atoms()] == [0, 0, 0, 0]
    assert [atom.formal_charge() for atom in sanitized.atoms()] == [0, 1, -1, 0]
    assert raw is not sanitized


def test_in_place_molecule_methods_mutate_receiver_and_return_none():
    mol = cosmolkit.Molecule.from_smiles("CCO")

    result = mol.add_hydrogens_()

    assert result is None
    assert mol.num_atoms() == 9
    assert mol.num_bonds() == 8

    assert mol.compute_2d_coordinates_() is None
    assert mol.has_2d_coordinates()

    assert mol.remove_hydrogens_(sanitize=False) is None
    assert mol.num_atoms() == 3
    assert mol.num_bonds() == 2


def test_in_place_molecule_methods_preserve_shared_python_values():
    mol = cosmolkit.Molecule.from_smiles("CCO")
    shared = mol.with_2d_coordinates()

    mol.add_hydrogens_()

    assert shared.num_atoms() == 3
    assert shared.num_bonds() == 2
    assert mol.num_atoms() == 9
    assert mol.num_bonds() == 8


def test_in_place_sanitize_and_kekulize_methods_match_value_methods():
    raw = cosmolkit.Molecule.from_smiles("CN(=O)=O", sanitize=False)
    expected = raw.sanitize()

    raw.sanitize_()

    assert raw.to_smiles() == expected.to_smiles()

    benzene = cosmolkit.Molecule.from_smiles("c1ccccc1")
    expected_kekule = benzene.with_kekulized_bonds(clear_aromatic_flags=False)

    benzene.kekulize_(clear_aromatic_flags=False)

    assert benzene.to_smiles(kekule=True) == expected_kekule.to_smiles(kekule=True)


def test_sdf_sanitize_false_is_supported_and_not_silently_sanitized():
    sdf = """ethane
     COSMolKit      2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
$$$$
"""
    mol = cosmolkit.Molecule.read_sdf_from_str(sdf, sanitize=False)

    assert mol.num_atoms() == 2
    assert mol.num_bonds() == 1


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
    mol = cosmolkit.Molecule.from_smiles("CCO").with_2d_coordinates()

    coords = mol.coordinates_2d()
    bounds = mol.dg_bounds_matrix()

    assert isinstance(coords, np.ndarray)
    assert coords.shape == (3, 3)
    assert isinstance(bounds, np.ndarray)
    assert bounds.shape == (3, 3)


def test_sdf_string_export_is_explicit_about_2d_or_3d_coordinates():
    mol = cosmolkit.Molecule.from_smiles("CCO")

    sdf_text = mol.to_2d_sdf_string(format="v2000")
    assert "2D" in sdf_text.splitlines()[1]
    assert not mol.has_2d_coordinates()
    assert "V2000" in mol.to_2d_sdf_string(format="v2000", include_stereo=False)
    aromatic_sdf = cosmolkit.Molecule.from_smiles("c1ccccc1").to_2d_sdf_string(
        format="v2000", kekulize=False
    )
    assert "  4  0" in aromatic_sdf

    with pytest.raises(ValueError, match="3D coordinates are required"):
        mol.to_3d_sdf_string(format="v2000")

    methane_3d = """methane_3d
  COSMolKit      3D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6291    0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6291   -0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6291    0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6291   -0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
M  END
$$$$
"""
    mol_3d = cosmolkit.Molecule.read_sdf_from_str(methane_3d, coordinate_dim="3d")
    assert "3D" in mol_3d.to_3d_sdf_string(format="v2000").splitlines()[1]
    assert "2D" in mol_3d.to_2d_sdf_string(format="v2000").splitlines()[1]
    assert mol_3d.num_conformers() == 1
    assert not mol_3d.has_2d_coordinates()


def test_molfile_read_matches_single_record_without_sdf_separator(tmp_path: Path):
    mol = cosmolkit.Molecule.from_smiles("CCO").with_2d_coordinates()
    mol_text = mol.to_2d_sdf_string(format="v2000").replace("$$$$\n", "")

    parsed = cosmolkit.Molecule.read_mol_from_str(mol_text, coordinate_dim="2d")
    assert parsed.to_smiles() == "CCO"
    assert parsed.has_2d_coordinates()

    path = tmp_path / "ethanol.mol"
    path.write_text(mol_text, encoding="utf-8")
    from_file = cosmolkit.Molecule.read_mol(str(path), coordinate_dim="2d")
    assert from_file.to_smiles() == "CCO"
    assert from_file.has_2d_coordinates()


def test_structural_array_access_supports_numpy_operations():
    mol = cosmolkit.Molecule.from_smiles("CCO").with_2d_coordinates()

    coords = mol.coordinates_2d()
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
    mol = cosmolkit.Molecule.from_smiles("CCO").with_2d_coordinates()

    tensor = torch.from_numpy(mol.coordinates_2d())

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
    one = cosmolkit.Molecule.read_sdf_from_str(sdf)
    many = cosmolkit.MoleculeBatch.read_sdf_records_from_str(sdf)

    assert len(one) == 2
    assert len(many) == 1
    first_many = many[0]
    assert first_many is not None
    assert len(first_many) == 2


def test_read_sdf_coordinate_dim_can_be_forced():
    sdf = """flat
     COSMolKit      2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""
    mol_2d = cosmolkit.Molecule.read_sdf_from_str(sdf, coordinate_dim="2d")
    mol_3d = cosmolkit.Molecule.read_sdf_from_str(sdf, coordinate_dim="3d")

    assert mol_2d.coordinates_2d().shape == (1, 3)
    assert np.allclose(mol_3d.coordinates_3d(), np.array([[0.0, 0.0, 0.0]]))


def test_sdf_coordinate_dim_is_applied_to_file_dataset_reader_and_batch(tmp_path: Path):
    sdf = """flat
     COSMolKit      2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""
    path = tmp_path / "flat.sdf"
    path.write_text(sdf)

    single = cosmolkit.Molecule.read_sdf(str(path), coordinate_dim="3d")
    from_file = cosmolkit.MoleculeBatch.read_sdf(
        str(path), errors="raise", n_jobs=1, coordinate_dim="3d"
    )[0]
    from_file_progress = cosmolkit.MoleculeBatch.read_sdf(
        str(path), errors="raise", progress_bar=True, coordinate_dim="3d"
    )[0]
    from_text_batch = cosmolkit.MoleculeBatch.read_sdf_records_from_str(
        sdf, errors="raise", n_jobs=1, coordinate_dim="3d"
    )[0]
    dataset = cosmolkit.SdfDataset.open(str(path), coordinate_dim="3d")
    from_dataset_index = dataset[0].molecule()
    from_dataset_iter = next(iter(dataset)).molecule()
    from_dataset_batch = next(dataset.batches(size=1, errors="raise"))[0]
    from_reader_batch = next(
        cosmolkit.SdfReader.open(str(path), coordinate_dim="3d").batches(
            size=1, errors="raise"
        )
    )[0]

    molecules = [
        single,
        from_file,
        from_file_progress,
        from_text_batch,
        from_dataset_index,
        from_dataset_iter,
        from_dataset_batch,
        from_reader_batch,
    ]
    for mol in molecules:
        assert mol is not None
        assert mol.num_conformers() == 1
        assert np.allclose(mol.coordinates_3d(), np.array([[0.0, 0.0, 0.0]]))


def test_molecule_batch_read_sdf_reads_file_records(tmp_path: Path):
    sdf = """ethane
     COSMolKit      2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
$$$$
oxygen
     COSMolKit      2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""
    path = tmp_path / "mols.sdf"
    path.write_text(sdf)

    batch = cosmolkit.MoleculeBatch.read_sdf(str(path), errors="raise", n_jobs=1)

    assert len(batch) == 2
    assert batch.valid_mask() == [True, True]
    first = batch[0]
    second = batch[1]
    assert first is not None
    assert second is not None
    assert len(first) == 2
    assert len(second) == 1


def test_molecule_batch_read_sdf_progress_bar_reads_file_records(tmp_path: Path):
    sdf = """carbon
     COSMolKit      2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
oxygen
     COSMolKit      2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""
    path = tmp_path / "progress.sdf"
    path.write_text(sdf)

    batch = cosmolkit.MoleculeBatch.read_sdf(
        str(path), errors="raise", progress_bar=True
    )

    assert len(batch) == 2
    assert batch.valid_mask() == [True, True]


def test_sdf_dataset_supports_indexing_iteration_and_batches(tmp_path: Path):
    sdf = """carbon
     COSMolKit      2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
>  <supplier_id>
D008

$$$$
oxygen
     COSMolKit      2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
nitrogen
     COSMolKit      2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""
    path = tmp_path / "indexed.sdf"
    path.write_text(sdf)

    ds = cosmolkit.SdfDataset.open(str(path))

    assert len(ds) == 3
    assert ds.metadata(0).title() == "carbon"
    assert ds[1].title() == "oxygen"
    assert ds[-1].title() == "nitrogen"
    assert ds[0].data_field("supplier_id") == "D008"
    head = ds[:2]
    assert len(head) == 2
    assert head.valid_mask() == [True, True]
    masked = ds[[True, False, True]]
    assert len(masked) == 2
    assert [mol.to_smiles() for mol in require_valid_molecules(masked)] == ["C", "N"]
    selected = ds[[2, 0]]
    assert [mol.to_smiles() for mol in require_valid_molecules(selected)] == ["N", "C"]
    with pytest.raises(TypeError, match="scalar boolean"):
        ds[True]
    with pytest.raises(TypeError, match="must not mix bool and int"):
        ds[[True, 1]]
    with pytest.raises(IndexError, match="boolean mask length"):
        ds[[True, False]]
    assert [record.title() for record in ds] == ["carbon", "oxygen", "nitrogen"]

    batches = list(ds.batches(size=2, errors="raise"))
    assert [len(batch) for batch in batches] == [2, 1]
    assert batches[0].valid_mask() == [True, True]
    assert batches[1].valid_mask() == [True]

    reader_batches = list(cosmolkit.SdfReader.open(str(path)).batches(size=2))
    assert [len(batch) for batch in reader_batches] == [2, 1]


def test_protein_from_pdb_str_projects_protein_chains():
    pdb = """\
ATOM      1  N   ALA A   1      11.104  13.207   9.900  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.210  13.912  10.555  1.00 20.00           C  
ATOM      3  C   ALA A   1      13.470  13.079  10.413  1.00 20.00           C  
HETATM    4  O   HOH A   2      18.000  10.000   8.000  1.00 10.00           O  
HETATM    5  C1  LIG B   1      18.500  11.000   8.500  1.00 10.00           C  
"""
    protein = cosmolkit.Protein.from_pdb_str(pdb)

    assert protein.num_chains() == 1
    assert protein.num_residues() == 1
    assert protein.num_atoms() == 3
    chain = protein[0]
    assert chain.kind() == "Protein"
    assert [residue.name() for residue in chain.residues()] == ["ALA"]
    assert [residue.kind() for residue in chain.residues()] == ["AminoAcid"]
    residue = chain.residues()[0]
    assert residue.code() == cosmolkit.ResidueCode.ALA
    assert residue.info().kind() == cosmolkit.ResidueInfoKind.AA
    assert residue.info().is_amino_acid()
    assert residue.one_letter_code() == "A"
    assert residue.fasta_code() == "A"
    assert cosmolkit.residue_code_from_name("TRY") == cosmolkit.ResidueCode.TRP
    assert cosmolkit.find_tabulated_residue_idx("h2o") == 154
    assert cosmolkit.get_residue_info(154).code() == cosmolkit.ResidueCode.HOH
    assert cosmolkit.find_tabulated_residue("MSE").fasta_code() == "X"
    assert cosmolkit.expand_one_letter("m", cosmolkit.ResidueInfoKind.AA) == "MET"
    assert cosmolkit.expand_protein_one_letter("m") == "MET"
    assert cosmolkit.expand_one_letter_sequence(
        "ACD(MSE)", cosmolkit.ResidueInfoKind.AA
    ) == ["ALA", "CYS", "ASP", "MSE"]
    assert cosmolkit.expand_protein_one_letter_string("m") == ["MET"]


def test_protein_from_pdb_str_does_not_fail_on_supported_metal_hetatm():
    pdb = """\
ATOM      1  N   ALA A   1      11.104  13.207   9.900  1.00 20.00           N
HETATM    2 HG    HG     2      -2.213  10.563  24.265  1.00 32.73          HG
HETATM    3 CD    CD     3      -3.467  18.396  77.649  0.50 39.48          CD
"""
    protein = cosmolkit.Protein.from_pdb_str(pdb)

    assert protein.num_atoms() == 1


def test_molecule_batch_read_sdf_file_error_modes(tmp_path: Path):
    sdf = """ok
     COSMolKit      2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
bad
     COSMolKit      2D

bad counts
M  END
$$$$
"""
    path = tmp_path / "with_error.sdf"
    path.write_text(sdf)

    kept = cosmolkit.MoleculeBatch.read_sdf(str(path), errors="keep", n_jobs=1)
    assert len(kept) == 2
    assert kept.valid_mask() == [True, False]
    assert kept.errors()[0].index() == 1
    assert kept.errors()[0].operation() == "read_sdf"

    skipped = cosmolkit.MoleculeBatch.read_sdf(str(path), errors="keep", n_jobs=1).filter_valid()
    assert len(skipped) == 1
    assert skipped.valid_mask() == [True]

    with pytest.raises(cosmolkit.BatchValidationError):
        cosmolkit.MoleculeBatch.read_sdf(str(path), errors="raise", n_jobs=1)


def test_without_hydrogens_filters_coordinate_rows_for_removed_atoms():
    sdf = """interleaved_h
  COSMolKit  3D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  3  4  1  0
M  END
$$$$
"""
    mol_3d = cosmolkit.Molecule.read_sdf_from_str(sdf, coordinate_dim="3d")
    removed_3d = mol_3d.without_hydrogens(sanitize=False)

    assert len(removed_3d) == 2
    assert removed_3d.coordinates_3d().shape == (2, 3)
    assert np.allclose(removed_3d.coordinates_3d(), np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]))

    mol_2d = cosmolkit.Molecule.read_sdf_from_str(sdf, coordinate_dim="2d")
    removed_2d = mol_2d.without_hydrogens(sanitize=False)

    assert len(removed_2d) == 2
    assert removed_2d.coordinates_2d().shape == (2, 3)
    assert np.allclose(removed_2d.coordinates_2d(), np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]))

    both = mol_3d.with_2d_coordinates().without_hydrogens(sanitize=False)
    assert len(both) == 2
    assert both.coordinates_2d().shape == (2, 3)
    assert both.coordinates_3d().shape == (2, 3)


def test_public_transforms_keep_coordinate_rows_aligned_with_atoms():
    base_2d = cosmolkit.Molecule.from_smiles("CCO").with_2d_coordinates()
    assert_coordinate_rows_match_atoms(base_2d)
    assert_coordinate_rows_match_atoms(base_2d.with_hydrogens())
    assert_coordinate_rows_match_atoms(base_2d.sanitize())

    benzene = cosmolkit.Molecule.from_smiles("c1ccccc1").with_2d_coordinates()
    assert_coordinate_rows_match_atoms(benzene.with_kekulized_bonds(clear_aromatic_flags=False))

    sdf = """interleaved_h
  COSMolKit  3D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  3  4  1  0
M  END
$$$$
"""
    base_3d = cosmolkit.Molecule.read_sdf_from_str(sdf, coordinate_dim="3d")
    assert_coordinate_rows_match_atoms(base_3d)
    assert_coordinate_rows_match_atoms(base_3d.without_hydrogens(sanitize=False))

    base_both = base_3d.with_2d_coordinates()
    assert_coordinate_rows_match_atoms(base_both)
    assert_coordinate_rows_match_atoms(base_both.without_hydrogens(sanitize=False))

    batch = cosmolkit.MoleculeBatch.read_sdf_records_from_str(
        sdf, coordinate_dim="3d", errors="raise", n_jobs=2
    )
    batch_removed = batch.without_hydrogens(errors="raise", n_jobs=2)
    first_removed = batch_removed[0]
    assert first_removed is not None
    assert_coordinate_rows_match_atoms(first_removed)


def test_molecule_batch_processes_in_order_and_preserves_single_molecule_inputs():
    mol = cosmolkit.Molecule.from_smiles("CCO")
    batch = cosmolkit.MoleculeBatch.from_smiles_list(
        ["CCO", "c1ccccc1", "CC(=O)O"], errors="keep", n_jobs=2
    )

    prepared = batch.with_hydrogens(errors="keep", n_jobs=2).with_2d_coordinates(
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
    assert mapped.to_smiles(canonical=False, ignore_atom_map_numbers=True) == "CO"
    assert mapped.to_smiles(ignore_atom_map_numbers=True) == "[CH3:7][OH:2]"


def test_molecule_batch_parallel_jobs_are_value_style_and_inherited():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(["CCO", "CC"], errors="keep")

    configured = batch.with_parallel_jobs(2)
    prepared = configured.with_hydrogens(errors="keep").with_2d_coordinates(errors="keep")

    assert batch.parallel_jobs() is None
    assert configured.parallel_jobs() == 2
    assert prepared.parallel_jobs() == 2
    first_smiles = prepared.to_smiles_list(n_jobs=1)[0]
    assert first_smiles is not None
    assert first_smiles.startswith("[H]O")
    assert configured.with_parallel_jobs(None).parallel_jobs() is None
    with pytest.raises(ValueError, match="n_jobs must be >= 1"):
        batch.with_parallel_jobs(0)


def test_molecule_batch_progress_bar_is_value_style_and_overridable():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(["CCO", "CC"], errors="keep")

    configured = batch.with_progress_bar(True)
    prepared = configured.with_2d_coordinates(errors="keep", progress_bar=False)

    assert batch.progress_bar() is None
    assert configured.progress_bar() is True
    assert prepared.progress_bar() is True
    assert configured.with_progress_bar(None).progress_bar() is None
    assert prepared.to_smiles_list(progress_bar=False) == ["CCO", "CC"]


def test_molecule_batch_to_list_indexing_and_iteration_return_molecules_or_none():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(["CCO", "C1CC", "CC"], errors="keep")

    values = batch.to_list()

    assert [value.to_smiles() if value is not None else None for value in values] == [
        "CCO",
        None,
        "CC",
    ]
    first_batch = batch[0]
    assert first_batch is not None
    assert first_batch.to_smiles() == "CCO"
    last_batch = batch[-1]
    assert last_batch is not None
    assert last_batch.to_smiles() == "CC"
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

    report = batch.with_2d_coordinates(errors="keep", n_jobs=2).to_sdf(
        str(tmp_path / "valid.sdf"), errors="keep", n_jobs=2
    )
    assert report.success() == 1
    assert report.failed() == 0
    assert (tmp_path / "valid.sdf").exists()


def test_molecule_batch_exports_use_custom_filenames(tmp_path: Path):
    batch = (
        cosmolkit.MoleculeBatch.from_smiles_list(["CCO", "C1CC", "CC"], errors="keep")
        .with_parallel_jobs(2)
        .with_2d_coordinates(errors="keep")
    )

    image_report = batch.to_images(
        str(tmp_path / "images"),
        format="svg",
        errors="keep",
        filenames=["ethanol", "bad.svg", None],
    )
    assert image_report.success() == 2
    assert (tmp_path / "images" / "ethanol.svg").exists()
    assert (tmp_path / "images" / "mol_2.svg").exists()

    sdf_report = batch.to_sdf_files(
        str(tmp_path / "sdf"),
        format="v2000",
        errors="keep",
        filenames=["ethanol", "bad.sdf", None],
    )
    assert sdf_report.success() == 2
    assert (tmp_path / "sdf" / "ethanol.sdf").exists()
    assert (tmp_path / "sdf" / "mol_2.sdf").exists()

    with pytest.raises(cosmolkit.BatchValidationError, match="invalid filename"):
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
    assert batch.to_smiles_list(
        canonical=False, ignore_atom_map_numbers=True, n_jobs=2
    ) == [
        "*C",
        "[13CH3][C@H](F)Cl",
    ]
    assert batch.to_smiles_list(ignore_atom_map_numbers=True, n_jobs=2) == [
        "[*:1]C",
        "[13CH3:7][C@H](F)Cl",
    ]
    assert batch.to_smiles_list(all_bonds_explicit=True, n_jobs=2)[0] == "C-[*:1]"
    assert batch.to_smiles_list(rooted_at_atom=0, n_jobs=2) == [
        "[*:1]C",
        "[13CH3:7][C@H](F)Cl",
    ]
    first_valid = batch.filter_valid()[0]
    assert first_valid is not None
    assert first_valid.to_smiles(rooted_at_atom=0) == "[*:1]C"


def test_molecule_batch_raise_aggregates_errors():
    with pytest.raises(cosmolkit.BatchValidationError) as excinfo:
        cosmolkit.MoleculeBatch.from_smiles_list(["C1CC", "N1"], errors="raise", n_jobs=2)

    assert "batch validation failed" in str(excinfo.value)
    errors = excinfo.value.errors()
    assert [error.operation() for error in errors] == [
        "batch.from_smiles_list",
        "batch.from_smiles_list",
    ]
    assert all(error.message() for error in errors)


def test_batch_errors_expose_intenum_types_and_mode_enum_is_accepted():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(
        ["CCO", "C1CC"],
        errors=cosmolkit.BatchErrorMode.KEEP,
        n_jobs=2,
    )

    errors = batch.errors()
    assert len(errors) == 1
    assert errors[0].operation() == "batch.from_smiles_list"
    assert errors[0].message() == "unclosed ring"
    assert errors[0].as_dict() == [
        ("index", "1"),
        ("operation", "batch.from_smiles_list"),
        ("message", errors[0].message()),
    ]
    assert cosmolkit.BATCH_ERROR_MODE_MAP["keep"] == cosmolkit.BatchErrorMode.KEEP

    filtered = batch.filter_valid().with_2d_coordinates(errors=cosmolkit.BatchErrorMode.KEEP)
    assert len(filtered) == 1
