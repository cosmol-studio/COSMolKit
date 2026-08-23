from pathlib import Path

import pytest

import cosmolkit


PDB = """\
HEADER    PYTHON MMCIF WRITER TEST               23-AUG-26   0XYZ
ATOM      1  N   ALA A   1      11.104  13.207   9.900  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.210  13.912  10.555  1.00 21.00           C
TER       3      ALA A   1
END
"""


def test_biostructure_mmcif_string_writer_defaults_reparse_without_mutation() -> None:
    structure = cosmolkit.BioStructure.from_pdb_str(PDB)
    before = (
        repr(structure),
        structure.name(),
        structure.num_models(),
        structure.num_chains(),
        structure.num_residues(),
        structure.num_atoms(),
        [atom.position() for atom in structure.atoms()],
    )

    output = structure.to_mmcif()
    reparsed = cosmolkit.BioStructure.from_mmcif_str(output, "python-roundtrip.cif")

    assert output.startswith("data_")
    assert "_atom_site.Cartn_x" in output
    assert reparsed.num_models() == structure.num_models()
    assert reparsed.num_chains() == structure.num_chains()
    assert reparsed.num_residues() == structure.num_residues()
    assert reparsed.num_atoms() == structure.num_atoms()
    assert [atom.position() for atom in reparsed.atoms()] == [
        atom.position() for atom in structure.atoms()
    ]
    assert before == (
        repr(structure),
        structure.name(),
        structure.num_models(),
        structure.num_chains(),
        structure.num_residues(),
        structure.num_atoms(),
        [atom.position() for atom in structure.atoms()],
    )


def test_biostructure_mmcif_options_and_file_writer_share_serialization(
    tmp_path: Path,
) -> None:
    structure = cosmolkit.BioStructure.from_pdb_str(PDB)
    groups = cosmolkit.MmcifOutputGroups(False)
    groups.atoms = True
    groups.block_name = True
    groups.group_pdb = True
    options = cosmolkit.MmcifWriteOptions()
    options.groups = groups
    options.compact = True

    expected = structure.to_mmcif(options)
    output_path = tmp_path / "selected.cif"
    structure.write_mmcif(str(output_path), options)

    assert output_path.read_text(encoding="utf-8") == expected
    assert "_atom_site." in expected
    assert "_entry." not in expected
    assert "_entity." not in expected
    assert options.groups.atoms is True
    assert options.groups.entry is False
    assert cosmolkit.MmcifWriteOptions().groups.auth_all is False


def test_biostructure_mmcif_writer_errors_and_api_ownership(tmp_path: Path) -> None:
    structure = cosmolkit.BioStructure.from_pdb_str(PDB)

    with pytest.raises(OSError):
        structure.write_mmcif(str(tmp_path / "missing" / "output.cif"))

    assert not hasattr(cosmolkit.Protein, "to_mmcif")
    assert not hasattr(cosmolkit.Protein, "write_mmcif")
    assert not hasattr(cosmolkit.Molecule, "to_mmcif")
    assert not hasattr(cosmolkit.Molecule, "write_mmcif")
