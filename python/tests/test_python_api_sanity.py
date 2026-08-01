from __future__ import annotations

import pickle
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pytest

import cosmolkit


def atom_signature(mol: cosmolkit.Molecule) -> list[tuple[int, int, int, str, int | None]]:
    return [
        (
            atom.idx(),
            atom.atomic_num(),
            atom.formal_charge(),
            atom.chiral_tag().name,
            atom.isotope(),
        )
        for atom in mol.atoms()
    ]


def bond_signature(
    mol: cosmolkit.Molecule,
) -> list[tuple[int, int, int, str, str, str, bool]]:
    return [
        (
            bond.idx(),
            bond.begin_atom_idx(),
            bond.end_atom_idx(),
            bond.bond_type().name,
            bond.bond_dir().name,
            bond.stereo().name,
            bond.is_aromatic(),
        )
        for bond in mol.bonds()
    ]


def assert_coord_rows_match_atoms(mol: cosmolkit.Molecule) -> None:
    atom_count = len(mol)
    if mol.has_2d_coordinates():
        coords2d = mol.coordinates_2d()
        assert coords2d.shape == (atom_count, 3), coords2d.shape
        assert np.allclose(coords2d[:, 2], 0.0)
    for conformer_index in range(mol.num_conformers()):
        coords3d = mol.coordinates_3d(conformer_index)
        assert coords3d.shape == (atom_count, 3), coords3d.shape


def regression_fixture_text(name: str) -> str:
    return (
        Path(__file__).resolve().parents[2]
        / "testdata"
        / "stereo"
        / "fixtures"
        / name
    ).read_text(encoding="utf-8")


KEKULE_BENZENE_MOL = """kekule_benzene
  COSMolKit      2D

  6  6  0  0  0  0  0  0  0  0999 V2000
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  1  2  0
M  END
"""


def assert_kekule_benzene_unsanitized(mol: cosmolkit.Molecule) -> None:
    assert len(mol) == 6
    assert [bond.bond_type() for bond in mol.bonds()] == [
        cosmolkit.BondOrder.SINGLE,
        cosmolkit.BondOrder.DOUBLE,
        cosmolkit.BondOrder.SINGLE,
        cosmolkit.BondOrder.DOUBLE,
        cosmolkit.BondOrder.SINGLE,
        cosmolkit.BondOrder.DOUBLE,
    ]
    assert not any(atom.is_aromatic() for atom in mol.atoms())
    assert not any(bond.is_aromatic() for bond in mol.bonds())


def test_smiles_and_value_semantics_preserve_expected_graph_features():
    base = cosmolkit.Molecule.from_smiles("F[C@H](Cl)[13CH3:7]")
    assert len(base) == 4
    assert base.to_smiles(canonical=False) == "F[C@H](Cl)[13CH3:7]"
    assert base.to_smiles(isomeric_smiles=False) == "FC(Cl)[CH3:7]"
    assert base.to_smiles(canonical=False, ignore_atom_map_numbers=True) == "F[C@H](Cl)[13CH3]"
    assert any(atom.chiral_tag() != cosmolkit.ChiralTag.CHI_UNSPECIFIED for atom in base.atoms())
    assert atom_signature(base)[-1][-1] == 13

    with_h = base.with_hydrogens()
    assert with_h is not base
    assert len(with_h) > len(base)
    assert with_h.to_smiles(canonical=False).count("[H]") > 0
    assert len(base) == 4

    removed = with_h.without_hydrogens(sanitize=False)
    assert removed is not with_h
    assert len(removed) == len(base)
    assert removed.to_smiles(canonical=False) == base.to_smiles(canonical=False)
    assert atom_signature(removed) == atom_signature(base)
    assert bond_signature(removed) == bond_signature(base)


def test_molecule_supports_python_pickle_roundtrip():
    mol = (
        cosmolkit.Molecule.from_smiles("F[C@H](Cl)[13CH3:7]")
        .with_hydrogens()
        .with_2d_coordinates()
    )

    restored = pickle.loads(pickle.dumps(mol))

    assert restored is not mol
    assert restored.to_smiles(canonical=False) == mol.to_smiles(canonical=False)
    assert atom_signature(restored) == atom_signature(mol)
    assert bond_signature(restored) == bond_signature(mol)
    assert restored.has_2d_coordinates()
    assert np.allclose(restored.coordinates_2d(), mol.coordinates_2d())


def test_molecule_pickle_rebuild_rejects_invalid_state():
    with pytest.raises(ValueError, match="unsupported Molecule pickle schema"):
        cosmolkit._rebuild_molecule_from_pickle(
            {
                "kind": "cosmolkit.Molecule",
                "pickle_schema": 999,
                "core_format": "cosmolkit-molecule-archive",
                "payload": b"",
            }
        )


def test_coordinate_and_sdf_roundtrip_behaviors_are_consistent():
    base = cosmolkit.Molecule.from_smiles("CCO")
    mol2d = base.with_2d_coordinates()
    assert not base.has_2d_coordinates()
    assert mol2d.has_2d_coordinates()
    assert mol2d.num_conformers() == 0
    assert_coord_rows_match_atoms(mol2d)

    sdf2d = mol2d.to_2d_sdf_string(format="v2000", include_stereo=True, kekulize=True)
    assert "V2000" in sdf2d
    assert "2D" in sdf2d.splitlines()[1]
    restored2d = cosmolkit.Molecule.read_sdf_from_str(sdf2d, coordinate_dim="2d")
    assert restored2d.to_smiles() == "CCO"
    assert_coord_rows_match_atoms(restored2d)

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
    mol3d = cosmolkit.Molecule.read_sdf_from_str(methane_3d, coordinate_dim="3d")
    assert mol3d.num_conformers() == 1
    assert not mol3d.has_2d_coordinates()
    assert_coord_rows_match_atoms(mol3d)

    sdf3d = mol3d.to_3d_sdf_string(format="v2000")
    assert "3D" in sdf3d.splitlines()[1]
    restored3d = cosmolkit.Molecule.read_sdf_from_str(sdf3d, coordinate_dim="3d")
    assert restored3d.num_conformers() == 1
    assert np.allclose(restored3d.coordinates_3d(), mol3d.coordinates_3d())

    both = mol3d.with_2d_coordinates()
    assert both.has_2d_coordinates()
    assert both.num_conformers() == 1
    assert_coord_rows_match_atoms(both)

    removed = both.without_hydrogens(sanitize=False)
    assert len(removed) == 1
    assert removed.has_2d_coordinates()
    assert removed.num_conformers() == 1
    assert removed.coordinates_2d().shape == (1, 3)
    assert removed.coordinates_3d().shape == (1, 3)


def test_read_mol_from_str_accepts_sanitize_false():
    mol = cosmolkit.Molecule.read_mol_from_str(
        KEKULE_BENZENE_MOL,
        coordinate_dim="2d",
        sanitize=False,
    )

    assert_kekule_benzene_unsanitized(mol)


def test_read_mol_accepts_sanitize_false(tmp_path: Path):
    path = tmp_path / "kekule_benzene.mol"
    path.write_text(KEKULE_BENZENE_MOL, encoding="utf-8")

    mol = cosmolkit.Molecule.read_mol(
        str(path),
        coordinate_dim="2d",
        sanitize=False,
    )

    assert_kekule_benzene_unsanitized(mol)


def test_read_sdf_from_str_accepts_sanitize_false():
    mol = cosmolkit.Molecule.read_sdf_from_str(
        f"{KEKULE_BENZENE_MOL}$$$$\n",
        coordinate_dim="2d",
        sanitize=False,
    )

    assert_kekule_benzene_unsanitized(mol)


def test_read_sdf_accepts_sanitize_false(tmp_path: Path):
    path = tmp_path / "kekule_benzene.sdf"
    path.write_text(f"{KEKULE_BENZENE_MOL}$$$$\n", encoding="utf-8")

    mol = cosmolkit.Molecule.read_sdf(
        str(path),
        coordinate_dim="2d",
        sanitize=False,
    )

    assert_kekule_benzene_unsanitized(mol)


def test_1aid_sdf_reader_preserves_rdkit_ring_stereo_after_remove_hs():
    mol = cosmolkit.Molecule.read_sdf_from_str(
        regression_fixture_text("1aid_ligand.sdf"),
        coordinate_dim="3d",
        sanitize=True,
        remove_hs=True,
    )

    tags = [atom.chiral_tag() for atom in mol.atoms()]
    assert tags[0] == cosmolkit.ChiralTag.CHI_TETRAHEDRAL_CW
    assert tags[3] == cosmolkit.ChiralTag.CHI_TETRAHEDRAL_CW
    centers = {center: ligands for center, ligands in mol.tetrahedral_stereo()}
    assert {0, 3}.issubset(centers)


def test_forcefield_wrappers_optimize_existing_3d_conformer_by_value():
    ethanol_3d = """ethanol_3d
  COSMolKit      3D

  9  8  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1000    1.2000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6000    0.9000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6000   -0.9000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    1.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.9000   -0.9000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.7000    0.0000    1.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.9000    1.2000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  1  4  1  0
  1  5  1  0
  1  6  1  0
  2  7  1  0
  2  8  1  0
  3  9  1  0
M  END
"""
    mol = cosmolkit.Molecule.read_mol_from_str(ethanol_3d, coordinate_dim="3d")
    original_coords = mol.coordinates_3d().copy()

    assert cosmolkit.uff_has_all_molecule_params(mol)
    assert mol.has_uff_params()

    result = mol.with_uff_optimized(max_iters=200)
    optimized = result.molecule()

    assert result.needs_more() is False
    assert result.energy() >= 0.0
    assert optimized is not mol
    assert optimized.num_conformers() == 1
    assert np.allclose(mol.coordinates_3d(), original_coords)
    assert not np.allclose(optimized.coordinates_3d(), original_coords)

    confs = cosmolkit.uff_optimize_molecule_confs(mol, max_iters=50)
    assert len(confs.conformer_results()) == mol.num_conformers()
    assert confs.molecule().num_conformers() == mol.num_conformers()

    mmff_available = cosmolkit.mmff_has_all_molecule_params(mol)
    assert mol.has_mmff_params() == mmff_available
    if mmff_available:
        mmff = mol.with_mmff_optimized(max_iters=50)
        assert mmff.molecule().num_conformers() == 1
        assert isinstance(mmff.needs_more(), bool)


def test_conformer_generation_python_api_exposes_native_embedding_and_parameters():
    base = cosmolkit.Molecule.from_smiles("CC(=O)NC").with_hydrogens()
    params = cosmolkit.EmbedParameters.etkdg_v3()
    params.random_seed = 0xF00D
    params.num_threads = 1
    params.max_iterations = 50
    params.track_failures = True

    result = base.with_3d_conformer_result(params)
    embedded = result.molecule()

    assert embedded is not base
    assert base.num_conformers() == 0
    assert embedded.num_conformers() == 1
    assert embedded.coordinates_3d().shape == (len(embedded), 3)
    assert result.conf_id() == 0
    assert result.ok() is True
    assert result.params().failures == params.failures
    assert params.failures and sum(params.failures) == 0
    assert params.et_version == 2
    assert params.use_exp_torsion_angle_prefs is True
    assert params.use_basic_knowledge is True
    assert params.use_macrocycle_torsions is True
    assert params.use_small_ring_torsions is False

    multi_params = cosmolkit.EmbedParameters.etkdg()
    multi_params.random_seed = 123
    multi_params.num_threads = 1
    multi_params.prune_rms_thresh = -1.0
    multi_params.enable_sequential_random_seeds = True
    multi_result = base.with_3d_conformers_result(3, multi_params)
    multi = multi_result.molecule()
    assert multi.num_conformers() == 3
    assert multi_result.conf_ids() == [0, 1, 2]
    assert multi_result.generated_count() == 3
    assert multi_result.requested_num_confs() == 3
    for conformer_index in range(multi.num_conformers()):
        assert multi.coordinates_3d(conformer_index).shape == (len(multi), 3)

    mutable = cosmolkit.Molecule.from_smiles("CCO").with_hydrogens()
    mutable_params = cosmolkit.EmbedParameters.kdg()
    mutable_params.random_seed = 77
    mutable.embed_3d_conformer_(mutable_params)
    assert mutable.num_conformers() == 1

    json_params = cosmolkit.EmbedParameters.dg()
    json_params.update_from_json(
        '{"randomSeed": 17, "useRandomCoords": true, "boxSizeMult": 3.5, "forceTransAmides": false, "trackFailures": true}'
    )
    assert json_params.random_seed == 17
    assert json_params.use_random_coords is True
    assert json_params.box_size_mult == 3.5
    assert json_params.force_trans_amides is False
    assert json_params.track_failures is True
    assert '"randomSeed":"17"' in json_params.to_json()

    mapped_params = cosmolkit.EmbedParameters.etkdg_v3()
    mapped_params.random_seed = 0xC0FFEE
    mapped_params.num_threads = 1
    mapped_params.use_random_coords = True
    mapped_params.coord_map = {
        0: (0.0, 0.0, 0.0),
        1: (0.0, 0.0, 1.5),
        2: (0.0, 1.5, 1.5),
    }
    mapped_params.cpci = {(0, 3): 0.5, (1, 4): -0.25}

    assert mapped_params.coord_map == {
        0: (0.0, 0.0, 0.0),
        1: (0.0, 0.0, 1.5),
        2: (0.0, 1.5, 1.5),
    }
    assert mapped_params.cpci == {(0, 3): 0.5, (1, 4): -0.25}

    mapped = base.with_3d_conformer(mapped_params)
    assert np.allclose(mapped.coordinates_3d()[0], [0.0, 0.0, 0.0])
    assert np.allclose(mapped.coordinates_3d()[1], [0.0, 0.0, 1.5])
    assert np.allclose(mapped.coordinates_3d()[2], [0.0, 1.5, 1.5])


def test_conformer_generation_failure_tracking_and_forcefield_post_optimization():
    fixture = (
        Path(__file__).resolve().parents[2]
        / "testdata/rdkit_builtin/fixtures/Code/GraphMol/DistGeomHelpers/chirality_failure_test.mol"
    )
    chiral = cosmolkit.Molecule.read_mol(str(fixture), coordinate_dim="auto", sanitize=True)
    chiral_params = cosmolkit.EmbedParameters.etkdg_v3()
    chiral_params.random_seed = 0xF00D
    chiral_params.num_threads = 1
    chiral_params.max_iterations = 50
    chiral_params.track_failures = True

    failed_result = chiral.with_3d_conformer_result(chiral_params)
    failed = failed_result.molecule()

    assert failed.num_conformers() == 0
    assert failed_result.conf_id() == -1
    assert failed_result.ok() is False
    assert failed_result.params().failures == chiral_params.failures
    assert chiral_params.failures
    assert sum(chiral_params.failures) > 0

    mol = cosmolkit.Molecule.from_smiles("CCO").with_hydrogens()
    embed_params = cosmolkit.EmbedParameters.etkdg_v3()
    embed_params.random_seed = 61453
    embed_params.num_threads = 1
    embedded = mol.with_3d_conformer(embed_params)
    coords_before = embedded.coordinates_3d().copy()

    uff = embedded.with_uff_optimized(max_iters=100)
    assert uff.molecule().num_conformers() == 1
    assert uff.energy() >= 0.0
    assert isinstance(uff.needs_more(), bool)
    assert uff.status_code() in (0, 1)
    assert np.allclose(embedded.coordinates_3d(), coords_before)
    assert not np.allclose(uff.molecule().coordinates_3d(), coords_before)

    if embedded.has_mmff_params():
        mmff = embedded.with_mmff_optimized(max_iters=50)
        assert mmff.molecule().num_conformers() == 1
        assert isinstance(mmff.needs_more(), bool)
        assert mmff.status_code() in (-1, 0, 1)


def test_editing_commit_boundary_matches_sanitize_behavior():
    invalid_editor = cosmolkit.Molecule.from_smiles("CC").edit()
    oxygen_a = invalid_editor.add_atom("O")
    oxygen_b = invalid_editor.add_atom("O")
    invalid_editor.add_bond(1, oxygen_a, order="double")
    invalid_editor.add_bond(1, oxygen_b, order="double")

    with pytest.raises(ValueError, match="sanitize failed"):
        _ = invalid_editor.commit()

    edited = invalid_editor.commit(sanitize=False)
    assert len(edited) == 4
    assert [bond.bond_type() for bond in edited.bonds()][-2:] == [
        cosmolkit.BondOrder.DOUBLE,
        cosmolkit.BondOrder.DOUBLE,
    ]

    valid_editor = cosmolkit.Molecule.from_smiles("CC").edit()
    oxygen = valid_editor.add_atom("O")
    valid_editor.add_bond(1, oxygen, order="single")
    valid = valid_editor.commit()
    assert valid.to_smiles(canonical=False) == "CCO"

    metal_editor = cosmolkit.Molecule.from_smiles("C").edit()
    hg = metal_editor.add_atom("Hg")
    metal = metal_editor.commit(sanitize=False)
    assert metal.atoms()[hg].atomic_num() == 80


def test_read_mol_stops_at_m_end_and_ignores_trailing_sdf_text(tmp_path: Path):
    mol2d = cosmolkit.Molecule.from_smiles("CCO").with_2d_coordinates()
    mol_path = tmp_path / "ethanol.mol"
    sdf_text = mol2d.to_2d_sdf_string(format="v2000").replace(
        "$$$$\n",
        ">  <supplier_id>\nD008\n\n$$$$\n",
    )
    _ = mol_path.write_text(sdf_text, encoding="utf-8")

    from_text = cosmolkit.Molecule.read_mol_from_str(sdf_text, coordinate_dim="2d")
    from_file = cosmolkit.Molecule.read_mol(str(mol_path), coordinate_dim="2d")

    assert from_text.to_smiles() == "CCO"
    assert from_file.to_smiles() == "CCO"
    assert from_text.has_2d_coordinates()
    assert from_file.has_2d_coordinates()
    assert len(from_file) == len(mol2d)
    assert_coord_rows_match_atoms(from_file)


def test_molfile_atomic_symbol_normalizes_uppercase_second_letter():
    mol_text = """brtest
  COSMolKit      2D

  3  2  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 BR  0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
M  END
"""
    mol = cosmolkit.Molecule.read_mol_from_str(mol_text, coordinate_dim="2d")

    assert [atom.atomic_num() for atom in mol.atoms()] == [6, 35, 35]


def test_molfile_invalid_mrv_sma_rejects_record():
    mol_text = """invalid_mrv_sma
  COSMolKit      2D

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  MRV SMA   1 MyDogHasFleas
M  END
"""
    with pytest.raises(ValueError, match="Cannot parse smarts"):
        _ = cosmolkit.Molecule.read_mol_from_str(mol_text, coordinate_dim="2d")


def test_read_mol2_from_str_and_file(tmp_path: Path):
    mol2_text = """@<TRIPOS>MOLECULE
ethanol
9 8
SMALL
NO_CHARGES
@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 C.3
2 C2 1.5 0.0 0.0 C.3
3 O1 3.0 0.0 0.0 O.3
4 H1 -0.5 0.9 0.0 H
5 H2 -0.5 -0.9 0.0 H
6 H3 0.0 0.0 1.0 H
7 H4 1.5 0.9 0.0 H
8 H5 1.5 0.0 1.0 H
9 H6 3.5 0.0 0.8 H
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
3 1 4 1
4 1 5 1
5 1 6 1
6 2 7 1
7 2 8 1
8 3 9 1
"""
    from_text = cosmolkit.Molecule.read_mol2_from_str(mol2_text)
    assert from_text.to_smiles() == "CCO"
    assert from_text.num_conformers() == 1

    path = tmp_path / "ethanol.mol2"
    _ = path.write_text(mol2_text, encoding="utf-8")
    from_file = cosmolkit.Molecule.read_mol2(str(path))
    assert from_file.to_smiles() == "CCO"
    assert from_file.num_conformers() == 1

    with pytest.raises(ValueError, match="unsupported MOL2 variant"):
        _ = cosmolkit.Molecule.read_mol2_from_str(mol2_text, variant="tripos")


def test_1aid_mol2_reader_preserves_rdkit_ring_stereo_after_remove_hs():
    mol = cosmolkit.Molecule.read_mol2_from_str(
        regression_fixture_text("1aid_ligand.mol2"),
        sanitize=True,
        remove_hs=True,
    )

    tags = [atom.chiral_tag() for atom in mol.atoms()]
    assert tags[0] == cosmolkit.ChiralTag.CHI_TETRAHEDRAL_CW
    assert tags[3] == cosmolkit.ChiralTag.CHI_TETRAHEDRAL_CW
    centers = {center: ligands for center, ligands in mol.tetrahedral_stereo()}
    assert {0, 3}.issubset(centers)


def test_fingerprint_and_stereo_outputs_are_structurally_reasonable():
    chiral = cosmolkit.Molecule.from_smiles("F[C@H](Cl)Br")
    stereo = chiral.tetrahedral_stereo()
    assert len(stereo) == 1
    center, ligands = stereo[0]
    assert center == 1
    assert ligands == [0, 2, 3, None]

    opposite = cosmolkit.Molecule.from_smiles("F[C@@H](Cl)Br")
    opposite_stereo = opposite.tetrahedral_stereo()
    assert opposite_stereo == [(1, [0, 3, 2, None])]

    fp = chiral.fingerprint_morgan(radius=2, n_bits=256)
    same = cosmolkit.Molecule.from_smiles("F[C@H](Cl)Br").fingerprint_morgan(radius=2, n_bits=256)
    other = cosmolkit.Molecule.from_smiles("CCO").fingerprint_morgan(radius=2, n_bits=256)
    assert fp.tanimoto(same) == 1.0
    assert 0.0 <= fp.tanimoto(other) < 1.0

    result = chiral.fingerprint_morgan_with_output(radius=2, n_bits=256)
    additional = result.additional_output()
    assert result.fingerprint().n_bits() == 256
    assert len(additional.atom_counts()) == len(chiral)
    assert isinstance(additional.bit_info_map(), dict)
    with pytest.raises(ValueError, match="avalon_fingerprint.*not implemented"):
        chiral.avalon_fingerprint(n_bits=256)
    with pytest.raises(ValueError, match="topological_fingerprint.*not implemented"):
        chiral.topological_fingerprint(n_bits=256)
    assert chiral.maccs_fingerprint().n_bits() == 166
    chiral.perceive_stereochemistry()


def test_smarts_parser_is_exposed_as_python_metadata():
    query = cosmolkit.parse_smarts("[#6]-O")

    assert query.num_atoms() == 2
    assert query.num_bonds() == 1
    assert query.ring_closures() == []
    assert "SmartsMolecule" in repr(query)

    with pytest.raises(ValueError, match="SMARTS|smarts|bracket"):
        _ = cosmolkit.parse_smarts("[#6")


def test_fragment_hash_pickle_and_scaffold_bindings_are_available():
    disconnected = cosmolkit.Molecule.from_smiles("CC.O")
    fragments = disconnected.fragments()
    assert len(fragments) == 2
    assert sorted(fragment.to_smiles() for fragment in fragments) == ["CC", "O"]
    assert disconnected.largest_fragment().to_smiles() == "CC"

    aromatic = cosmolkit.Molecule.from_smiles("c1ccccc1CCO")
    assert isinstance(aromatic.hash(), int)
    assert isinstance(aromatic.hash_with_ranks([0] * len(aromatic)), int)
    assert aromatic.murcko_scaffold().num_atoms() > 0
    assert aromatic.net_scaffold().num_atoms() > 0

    payload = aromatic.mol_to_binary()
    restored_method = cosmolkit.Molecule.mol_from_binary(payload)
    restored_fn = cosmolkit.mol_from_binary(payload)
    assert restored_method.to_smiles() == aromatic.to_smiles()
    assert restored_fn.to_smiles() == aromatic.to_smiles()
    assert cosmolkit.mol_to_binary(aromatic) == payload
    assert cosmolkit.version() == cosmolkit.__version__


def test_draw_and_substructure_bindings_are_available():
    mol = cosmolkit.Molecule.from_smiles("c1ccccc1O").with_2d_coordinates()
    png = mol.to_png(width=200, height=150)
    assert bytes(png).startswith(b"\x89PNG\r\n\x1a\n")

    pdb_block = mol.to_pdb_block()
    assert isinstance(pdb_block, str)
    assert "HETATM" in pdb_block or "ATOM" in pdb_block

    pdb_mol = cosmolkit.Molecule.from_pdb_block(
        """\
HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 10.00           C  
HETATM    2  O1  LIG A   1       1.200   0.000   0.000  1.00 10.00           O  
""",
        sanitize=False,
        remove_hs=False,
        proximity_bonding=True,
    )
    assert pdb_mol.num_atoms() == 2
    assert pdb_mol.num_bonds() == 1

    metal_pdb_mol = cosmolkit.Molecule.from_pdb_block(
        """\
HETATM    1 HG    HG     1      -2.213  10.563  24.265  1.00 32.73          HG
HETATM    2 CD    CD     1      -3.467  18.396  77.649  0.50 39.48          CD
""",
        sanitize=False,
        remove_hs=False,
        proximity_bonding=False,
    )
    assert [atom.atomic_num() for atom in metal_pdb_mol.atoms()] == [80, 48]

    mmcif_mol = cosmolkit.Molecule.from_mmcif_block(
        """\
data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 C C1 . LIG A 1 0.000 0.000 0.000
HETATM 2 O O1 . LIG A 1 1.200 0.000 0.000
""",
        sanitize=False,
        remove_hs=False,
        proximity_bonding=True,
    )
    assert mmcif_mol.num_atoms() == 2
    assert mmcif_mol.num_bonds() == 1

    xyz_mol = cosmolkit.Molecule.from_xyz_block(
        """\
3
water
O 0.000 0.000 0.000
H 0.758 0.000 0.504
H -0.758 0.000 0.504
"""
    )
    assert xyz_mol.num_atoms() == 3
    assert xyz_mol.num_bonds() == 0
    assert xyz_mol.num_conformers() == 1
    assert np.allclose(
        xyz_mol.coordinates_3d(),
        np.array(
            [
                [0.0, 0.0, 0.0],
                [0.758, 0.0, 0.504],
                [-0.758, 0.0, 0.504],
            ]
        ),
    )

    query = cosmolkit.Molecule.from_smiles("c1ccccc1")
    assert cosmolkit.has_substruct_match(mol, query) is True
    first = cosmolkit.get_substruct_match(mol, query)
    assert first is not None
    assert len(first.atom_mapping()) == len(query)
    matches = cosmolkit.get_substruct_matches(mol, query)
    assert len(matches) >= 1
    matches_limited = cosmolkit.get_substruct_matches_with_params(
        mol, query, max_matches=1, uniquify=True
    )
    assert len(matches_limited) == 1


def test_batch_api_combinations_preserve_order_shapes_and_record_alignment():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(
        ["CCO", "c1ccccc1", "bad", "F[C@H](Cl)Br"],
        sanitize=True,
        errors=cosmolkit.BatchErrorMode.KEEP,
        n_jobs=2,
    ).with_parallel_jobs(2).with_progress_bar(False)

    assert len(batch) == 4
    assert batch.parallel_jobs() == 2
    assert batch.progress_bar() is False
    assert batch.valid_mask() == [True, True, False, True]
    assert batch.valid_count() == 3
    assert batch.invalid_count() == 1
    assert batch.errors()[0].operation() == "batch.from_smiles_list"

    selected = batch[[0, 3]]
    assert len(selected) == 2
    sliced = batch[1:]
    assert len(sliced) == 3
    masked = batch[batch.valid_mask()]
    assert len(masked) == 3

    prepared = batch.with_hydrogens(errors="keep").with_2d_coordinates(errors="keep")
    assert prepared.valid_mask() == [True, True, False, True]
    mols = prepared.to_list()
    assert mols[2] is None
    assert all(m is None or m.has_2d_coordinates() for m in mols)
    for mol in mols:
        if mol is not None:
            assert_coord_rows_match_atoms(mol)

    smiles = prepared.to_smiles_list(canonical=False)
    assert smiles[0] is not None and "[H]" in smiles[0]
    assert smiles[2] is None

    dg = prepared.dg_bounds_matrix_list()
    assert mols[0] is not None
    assert dg[0] is not None and dg[0].shape[0] == len(mols[0])
    assert dg[2] is None

    fps = prepared.fingerprint_morgan_list(n_bits=128)
    assert fps[0] is not None
    assert fps[2] is None

    with TemporaryDirectory() as td:
        td_path = Path(td)
        img_report = prepared.to_images(
            str(td_path / "images"),
            format="svg",
            errors="keep",
            filenames=["ethanol.svg", "benzene.svg", None, "chiral.svg"],
        )
        sdf_report = prepared.to_sdf_files(
            str(td_path / "sdf"),
            format="v2000",
            errors="keep",
            filenames=["ethanol.sdf", "benzene.sdf", None, "chiral.sdf"],
        )
        assert img_report.total() == 4
        assert img_report.success() == 3
        assert sdf_report.success() == 3
        assert (td_path / "images" / "ethanol.svg").exists()
        assert (td_path / "images" / "benzene.svg").exists()
        assert (td_path / "sdf" / "chiral.sdf").exists()


def test_from_rdkit_roundtrip_keeps_expected_stereo_when_rdkit_is_available():
    rdkit_chem = pytest.importorskip("rdkit.Chem")
    rd_mol = rdkit_chem.MolFromSmiles("C/C=C/C")
    assert rd_mol is not None
    bridged = cosmolkit.Molecule.from_rdkit(rd_mol)
    assert bridged.to_smiles(canonical=False) == "C/C=C/C"
    assert bridged.bonds()[1].stereo() == cosmolkit.BondStereo.E
