"""Small BIO capability regressions against the released 0.3.0 interface.

Names follow the current public contract (from_pdb/read and *Ref). Expected
hierarchy, residue metadata and coordinates were observed with cosmolkit 0.3.0;
the tests require neither that installation nor an external reference corpus.
"""

from pathlib import Path
import runpy

import pytest

import cosmolkit as ck


PDB = """\
ATOM      1  N   MET A   1      11.104  13.207   9.900  1.00 20.00           N
ATOM      2  CA  MET A   1      12.210  13.912  10.555  1.00 20.00           C
ATOM      3  C   MET A   1      13.470  13.079  10.413  1.00 20.00           C
ATOM      4  N   GLY A   2      14.530  13.650  10.980  1.00 20.00           N
ATOM      5  CA  GLY A   2      15.790  12.920  10.910  1.00 20.00           C
HETATM    6  O   HOH A   3      18.000  10.000   8.000  1.00 10.00           O
HETATM    7  C1  LIG B   1      18.500  11.000   8.500  1.00 10.00           C
"""


def test_complete_structure_and_protein_projection_retain_baseline_capabilities(tmp_path):
    structure = ck.BioStructure.from_pdb(PDB)
    assert (structure.num_models(), structure.num_chains(), structure.num_residues(), structure.num_atoms()) == (1, 2, 4, 7)
    protein = structure.protein()
    assert (protein.num_models(), protein.num_chains(), protein.num_residues(), protein.num_atoms()) == (1, 1, 2, 5)
    chain = protein[0]
    assert isinstance(chain, ck.ProteinChainRef)
    assert (chain.id(), chain.kind()) == (0, "Protein")
    residues = chain.residues()
    assert all(isinstance(residue, ck.ProteinResidueRef) for residue in residues)
    assert [
        (r.id(), r.name(), int(r.code()), r.one_letter_code(), r.fasta_code(),
         r.canonical_one_letter_code(), r.is_standard(), r.is_modified_amino_acid(),
         int(r.parent_standard_code()))
        for r in residues
    ] == [
        (0, "MET", 16, "M", "M", "M", True, False, 16),
        (1, "GLY", 11, "G", "G", "G", True, False, 11),
    ]
    assert [
        (a.id(), a.name(), a.atomic_number(), a.position())
        for r in residues for a in r.atoms()
    ] == [
        (0, "N", 7, (11.104, 13.207, 9.9)),
        (1, "CA", 6, (12.21, 13.912, 10.555)),
        (2, "C", 6, (13.47, 13.079, 10.413)),
        (3, "N", 7, (14.53, 13.65, 10.98)),
        (4, "CA", 6, (15.79, 12.92, 10.91)),
    ]
    assert isinstance(protein.atoms()[0], ck.ProteinAtomRef)
    # File ingress and structural output must survive the API renaming too.
    path = tmp_path / "structure.pdb"
    path.write_text(PDB)
    assert ck.BioStructure.read(path).num_atoms() == 7
    assert ck.Protein.read(path).num_atoms() == 5
    mmcif = structure.to_mmcif()
    reparsed = ck.BioStructure.from_mmcif(mmcif)
    assert [a.position() for a in reparsed.atoms()] == [a.position() for a in structure.atoms()]
    assert reparsed.protein().num_atoms() == 5
    assert structure.to_molecule(sanitize=False, remove_hs=False).num_atoms() == 7
    assert structure.num_atoms() == 7


@pytest.mark.parametrize("name", ["protein_from_pdb", "protein_contact_summary", "structure_models"])
def test_bio_examples_execute_complete_workflows(name):
    runpy.run_path(str(Path(__file__).resolve().parents[1] / "examples" / f"{name}.py"), run_name="__main__")
