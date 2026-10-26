"""Small ordinary regressions for structural text IO and binding call forms."""
import pytest
import cosmolkit as ck

PDB = "HETATM    1  C1  LIG A   7       1.000   2.000   3.000  0.50 12.00           C  \nHETATM    2  H1  LIG A   7       1.800   2.000   3.000  1.00  4.00           H  \nCONECT    1    2\nEND\n"
MMCIF = "data_lig\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\nHETATM 1 C C1 . LIG A A 7 1.000 2.000 3.000 0.50 12.00\n"


def positions(text: str) -> list[tuple[float, ...]]:
    return [tuple(float(row[a:b]) for a, b in ((30, 38), (38, 46), (46, 54)))
            for row in text.splitlines() if row.startswith(("ATOM  ", "HETATM"))]


def test_structural_readers_default_params_kwargs_and_conveniences():
    params = ck.BioMoleculeParams(remove_hs=False, proximity_bonding=False)
    for read in (ck.Molecule.from_pdb, ck.mol_from_pdb):
        assert read(PDB).num_atoms == 1
        assert read(PDB, params).num_atoms == 2
        assert read(PDB, remove_hs=False, proximity_bonding=False).num_bonds == 1
        assert positions(read(PDB).to_pdb()) == [(1.0, 2.0, 3.0)]
    assert ck.Molecule.from_pdb_with_params(PDB, params).num_atoms == 2
    for read in (ck.Molecule.from_mmcif, ck.mol_from_mmcif):
        assert positions(read(MMCIF).to_pdb()) == [(1.0, 2.0, 3.0)]
        assert read(MMCIF, params).num_atoms == 1
        assert read(MMCIF, proximity_bonding=False).num_bonds == 0
    assert ck.Molecule.from_mmcif_with_params(MMCIF, params).num_atoms == 1


def test_pdb_unique_explicit_all_and_atomic_writable_options():
    source = ck.mol_from_smiles("C")
    assert positions(source.to_pdb()) == [(0.0, 0.0, 0.0)]
    builder = source.to_builder()
    builder.set_2d_coordinates([[1.25, -2.5]])
    builder.add_3d_conformer([[4.0, 5.0, 6.0]])
    molecule = builder.build()
    with pytest.raises(ck.MolecularIoError) as caught:
        molecule.to_pdb()
    cause = caught.value.__cause__
    assert isinstance(cause, ck.PdbWriteError)
    assert (cause.variant.name, cause.two_d, cause.three_d) == ("AmbiguousCoordinates", 1, 1)
    assert positions(molecule.to_pdb(coordinate_selection=ck.MolCoordinateSelection.TwoD(0))) == [(1.25, -2.5, 0.0)]
    params = ck.PdbWriteParams(coordinate_selection=ck.MolCoordinateSelection.ThreeD(0))
    assert positions(molecule.to_pdb(params)) == [(4.0, 5.0, 6.0)]
    params.coordinate_selection = None
    params.all_conformers = True
    assert positions(molecule.to_pdb_with_params(params)) == [(1.25, -2.5, 0.0), (4.0, 5.0, 6.0)]
    assert molecule.to_pdb(all_conformers=True).count("ENDMDL\n") == 2
    params.flavor = ck.PdbWriteFlags.MASTER
    assert "MASTER" in molecule.to_pdb(params)
    with pytest.raises(TypeError):
        setattr(params, "all_conformers", "yes")
    assert params.all_conformers is True
    assert positions(source.to_pdb()) == [(0.0, 0.0, 0.0)]
    coordinates = molecule.coordinates_2d
    assert coordinates is not None
    assert coordinates.tolist() == [[1.25, -2.5]]
    with pytest.raises(ck.MolecularIoError) as caught:
        molecule.to_pdb(coordinate_selection=ck.MolCoordinateSelection.ThreeD(91))
    cause = caught.value.__cause__
    assert isinstance(cause, ck.PdbWriteError)
    assert (cause.variant.name, cause.id) == ("MissingCoordinate", 91)
    with pytest.raises(ck.MolecularIoError) as caught:
        molecule.to_pdb(all_conformers=True, coordinate_selection=ck.MolCoordinateSelection.TwoD(0))
    cause = caught.value.__cause__
    assert isinstance(cause, ck.PdbWriteError)
    assert cause.variant.name == "ConflictingCoordinateSelection"


def test_pdb_aromatic_copy_and_typed_reader_failure():
    molecule = ck.mol_from_smiles("c1ccccc1")
    before = molecule.to_smiles()
    assert "CONECT    1    2    6    6\n" in molecule.to_pdb()
    assert molecule.to_smiles() == before
    with pytest.raises(ck.MolecularIoError) as caught:
        ck.mol_from_mmcif("not a CIF block")
    assert isinstance(caught.value.__cause__, ck.BioMmcifReadError)
    assert caught.value.variant.name == "MmcifRead"
    for name in ("from_pdb_block", "from_mmcif_block", "to_pdb_block"):
        assert not hasattr(ck.Molecule, name)
