from __future__ import annotations

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
    if mol.has_2d_coords():
        coords2d = mol.coords_2d()
        assert coords2d.shape == (atom_count, 3), coords2d.shape
        assert np.allclose(coords2d[:, 2], 0.0)
    for conformer_index in range(mol.num_conformers()):
        coords3d = mol.coords_3d(conformer_index)
        assert coords3d.shape == (atom_count, 3), coords3d.shape


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


def test_coordinate_and_sdf_roundtrip_behaviors_are_consistent():
    base = cosmolkit.Molecule.from_smiles("CCO")
    mol2d = base.with_2d_coords()
    assert not base.has_2d_coords()
    assert mol2d.has_2d_coords()
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
    assert not mol3d.has_2d_coords()
    assert_coord_rows_match_atoms(mol3d)

    sdf3d = mol3d.to_3d_sdf_string(format="v2000")
    assert "3D" in sdf3d.splitlines()[1]
    restored3d = cosmolkit.Molecule.read_sdf_from_str(sdf3d, coordinate_dim="3d")
    assert restored3d.num_conformers() == 1
    assert np.allclose(restored3d.coords_3d(), mol3d.coords_3d())

    both = mol3d.with_2d_coords()
    assert both.has_2d_coords()
    assert both.num_conformers() == 1
    assert_coord_rows_match_atoms(both)

    removed = both.without_hydrogens(sanitize=False)
    assert len(removed) == 1
    assert removed.has_2d_coords()
    assert removed.num_conformers() == 1
    assert removed.coords_2d().shape == (1, 3)
    assert removed.coords_3d().shape == (1, 3)


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


def test_read_mol_rejects_sdf_record_separator_but_accepts_plain_mol(tmp_path: Path):
    mol2d = cosmolkit.Molecule.from_smiles("CCO").with_2d_coords()
    mol_path = tmp_path / "ethanol.mol"
    mol2d.write_sdf(str(mol_path), format="v2000")

    with pytest.raises(ValueError, match="Extra non-molfile content after M  END"):
        cosmolkit.Molecule.read_mol(str(mol_path), coordinate_dim="2d")

    mol_text = mol_path.read_text(encoding="utf-8").replace("$$$$\n", "")
    mol_path.write_text(mol_text, encoding="utf-8")
    reread = cosmolkit.Molecule.read_mol(str(mol_path), coordinate_dim="2d")
    assert reread.has_2d_coords()
    assert len(reread) == len(mol2d)
    assert_coord_rows_match_atoms(reread)


def test_fingerprint_and_stereo_outputs_are_structurally_reasonable():
    chiral = cosmolkit.Molecule.from_smiles("F[C@H](Cl)Br")
    stereo = chiral.tetrahedral_stereo()
    assert len(stereo) == 1
    center, ligands = stereo[0]
    assert center == 1
    assert len(ligands) == 4

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
    assert chiral.avalon_fingerprint(n_bits=256).n_bits() == 256
    assert chiral.topological_fingerprint(n_bits=256).n_bits() == 256
    assert chiral.maccs_fingerprint().n_bits() == 166
    chiral.perceive_stereochemistry()


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
    mol = cosmolkit.Molecule.from_smiles("c1ccccc1O").with_2d_coords()
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

    prepared = batch.add_hydrogens(errors="keep").compute_2d_coords(errors="keep")
    assert prepared.valid_mask() == [True, True, False, True]
    mols = prepared.to_list()
    assert mols[2] is None
    assert all(m is None or m.has_2d_coords() for m in mols)
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
