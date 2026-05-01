from pathlib import Path

import cosmolkit


def test_with_2d_coords_returns_new_molecule_without_mutating_input():
    mol = cosmolkit.Molecule.from_smiles("CCO")

    mol_2d = mol.with_2d_coords()

    assert mol is not mol_2d
    assert mol.num_conformers() == 0
    assert mol_2d.num_conformers() == 1


def test_write_png_auto_prepares_drawing_without_mutating_input(tmp_path: Path):
    mol = cosmolkit.Molecule.from_smiles("CCO")
    output = tmp_path / "ethanol.png"

    mol.write_png(str(output))


    assert output.exists()
    assert output.stat().st_size > 0
    assert mol.num_conformers() == 0
