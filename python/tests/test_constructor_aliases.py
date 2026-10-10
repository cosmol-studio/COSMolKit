"""Constructor conveniences reference existing factories, not another API."""

import pytest
from rdkit import Chem

import cosmolkit as ck


ALIASES = [
    ("mol_from_smiles", "Molecule", "from_smiles"),
    ("mol_from_mol", "Molecule", "from_mol"),
    ("mol_from_sdf", "Molecule", "from_sdf"),
    ("mol_from_mol2", "Molecule", "from_mol2"),
    ("mol_from_xyz_block", "Molecule", "from_xyz_block"),
    ("mol_from_inchi", "Molecule", "from_inchi"),
    ("mol_from_binary", "Molecule", "from_binary"),
    ("mol_from_rdkit", "Molecule", "from_rdkit"),
    ("mols_from_smiles_list", "MoleculeBatch", "from_smiles_list"),
    ("mols_from_sdf_records", "MoleculeBatch", "from_sdf_records"),
    ("mols_from_records", "MoleculeBatch", "from_records"),
    ("reaction_from_smirks", "Reaction", "from_smirks"),
    ("bio_from_pdb", "BioStructure", "from_pdb"),
    ("bio_from_mmcif", "BioStructure", "from_mmcif"),
    ("protein_from_pdb", "Protein", "from_pdb"),
    ("protein_from_mmcif", "Protein", "from_mmcif"),
]

PDB = "ATOM      1  N   ALA A   1      11.104  13.207   9.900  1.00 20.00           N\nEND\n"
MOL2 = """@<TRIPOS>MOLECULE
methane
1 0 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
1 C1 0.0 0.0 0.0 C.3 1 MOL 0.0
@<TRIPOS>BOND
"""


@pytest.mark.parametrize("name,owner,method", ALIASES)
def test_constructor_shortcut_is_the_existing_callable(name, owner, method):
    original = getattr(getattr(ck, owner), method)
    alias = getattr(ck, name)
    assert alias == original
    assert alias.__doc__ == original.__doc__


def test_smiles_shortcut_inherits_bytes_params_and_keywords():
    assert ck.mol_from_smiles(b"CCO").to_smiles() == "CCO"
    params = ck.SmilesParseParams(sanitize=False)
    for mol in [
        ck.mol_from_smiles("CN(=O)=O", params),
        ck.mol_from_smiles("CN(=O)=O", sanitize=False),
    ]:
        assert mol.to_smiles(canonical=False) == "CN(=O)=O"
    with pytest.raises(TypeError):
        ck.mol_from_smiles("CCO", params, sanitize=True)
    with pytest.raises(TypeError):
        ck.mol_from_smiles("CCO", unknown_option=True)


def test_smiles_shortcut_preserves_error_type_message_and_context():
    errors = []
    for factory in [ck.Molecule.from_smiles, ck.mol_from_smiles]:
        with pytest.raises(Exception) as error:
            factory("C1CC")
        errors.append((type(error.value), error.value.args,
                       getattr(error.value, "domain", None), getattr(error.value, "kind", None)))
    assert errors[0] == errors[1]


def test_molecular_format_shortcuts_use_original_factories():
    source = ck.Molecule.from_smiles("C").with_2d_coordinates()
    inputs = {
        "mol_from_mol": source.to_mol(),
        "mol_from_sdf": source.to_sdf(),
        "mol_from_mol2": MOL2,
        "mol_from_xyz_block": "1\ncarbon\nC 0.0 0.0 0.0\n",
        "mol_from_inchi": "InChI=1S/CH4/h1H4",
        "mol_from_binary": source.to_binary(),
        "mol_from_rdkit": Chem.MolFromSmiles("C"),
    }
    methods = {name: method for name, _, method in ALIASES}
    for name, data in inputs.items():
        actual = getattr(ck, name)(data)
        expected = getattr(ck.Molecule, methods[name])(data)
        assert isinstance(actual, ck.Molecule)
        assert actual.to_smiles() == expected.to_smiles()
        assert actual.num_atoms() == expected.num_atoms()


def test_batch_shortcuts_inherit_execution_configuration_and_record_state():
    inputs = ["CN(=O)=O", "C1CC"]
    options = dict(sanitize=False, errors="keep", n_jobs=1)
    actual = ck.mols_from_smiles_list(inputs, **options)
    expected = ck.MoleculeBatch.from_smiles_list(inputs, **options)
    assert actual.to_smiles_list(canonical=False) == expected.to_smiles_list(canonical=False)
    assert len(actual.errors()) == len(expected.errors()) == 1
    mol = ck.mol_from_smiles("C").with_2d_coordinates()
    sdf = mol.to_sdf()
    assert ck.mols_from_sdf_records(sdf, errors="keep", n_jobs=1).to_smiles_list() == ["C"]
    assert ck.mols_from_records([mol], ck.BatchErrorMode.KEEP).to_smiles_list() == ["C"]


def test_reaction_shortcut_inherits_configuration():
    text = "[C:1]>>[C:1]"
    params = ck.ReactionParseParams(use_smiles=True)
    for reaction in [ck.reaction_from_smirks(text, params),
                     ck.reaction_from_smirks(text, use_smiles=True)]:
        assert isinstance(reaction, ck.Reaction)


@pytest.mark.parametrize("owner,pdb_name,cif_name", [
    (ck.BioStructure, "bio_from_pdb", "bio_from_mmcif"),
    (ck.Protein, "protein_from_pdb", "protein_from_mmcif"),
])
def test_bio_shortcuts_retain_types_coordinates_and_configuration(owner, pdb_name, cif_name):
    actual = getattr(ck, pdb_name)(PDB, ignore_ter=True)
    expected = owner.from_pdb(PDB, ck.BioPdbReadParams(ignore_ter=True))
    assert isinstance(actual, owner)
    assert actual.num_atoms() == expected.num_atoms() == 1
    assert [atom.position() for atom in actual.atoms()] == [atom.position() for atom in expected.atoms()]
    # Use the complete structure writer for both structural and Protein ingress.
    cif = ck.BioStructure.from_pdb(PDB).to_mmcif()
    assert isinstance(getattr(ck, cif_name)(cif), owner)
    assert getattr(ck, cif_name)(cif).num_atoms() == owner.from_mmcif(cif).num_atoms()
