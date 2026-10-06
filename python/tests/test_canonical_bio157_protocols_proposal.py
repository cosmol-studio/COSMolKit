"""Fixed canonical BIO protocol proposals; p1/ROOT acceptance is pending."""
import gc
import pickle
from types import MappingProxyType
import pytest
import cosmolkit as ck

PDB=("ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \n"
     "HETATM    2  O   HOH B   2       4.000   5.000   6.000  1.00 20.00           O  \nEND\n")

@pytest.mark.parametrize("name,alias",[("TRY","TRP"),("WAT","HOH"),("H2O","HOH"),("+A","DA"),("+C","DC"),("+G","DG"),("+I","DI"),("+T","DT"),("+U","DU"),("+N","DN")])
def test_residue_map_source_aliases(name,alias):
    assert isinstance(ck.RESIDUE_CODE_MAP,MappingProxyType)
    assert ck.RESIDUE_CODE_MAP[name] is getattr(ck.ResidueCode,alias)
    with pytest.raises(TypeError): ck.RESIDUE_CODE_MAP[name]=ck.ResidueCode.ALA

@pytest.mark.parametrize("kind",list(ck.ResidueInfoKind))
def test_residue_kind_map_pickle(kind):
    assert ck.RESIDUE_INFO_KIND_MAP[kind.name] is kind
    assert pickle.loads(pickle.dumps(kind)) is kind
    with pytest.raises(TypeError): ck.RESIDUE_INFO_KIND_MAP[kind.name]=kind

@pytest.mark.parametrize("name",["ALA","GLY","MSE","HOH","DA","A","UNK"])
def test_residue_owned_results(name):
    info=ck.find_residue_info(name)
    assert info.code()==ck.residue_code(name)
    assert ck.residue_info(ck.find_residue_info_index(name)).code()==info.code()
    assert pickle.loads(pickle.dumps(info.code())) is info.code()
    assert isinstance(repr(info),str)
    with pytest.raises(AttributeError): info.fake_field=1

def test_expansion_and_checked_error_boundaries():
    assert ck.expand_one_letter("A",ck.ResidueInfoKind.AA)=="ALA"
    assert ck.expand_one_letter_sequence("AG",ck.ResidueInfoKind.AA)==["ALA","GLY"]
    with pytest.raises(ValueError):ck.expand_one_letter("AA",ck.ResidueInfoKind.AA)
    with pytest.raises(ValueError):ck.expand_one_letter("A",99)
    with pytest.raises(IndexError):ck.residue_info(1_000_000)

@pytest.mark.parametrize("name", ["ALA", "ala", "not_a_residue"])
def test_residue_identity_preserves_input_and_returns_owned_info(name):
    identity = ck.ResidueIdentity.new(name)
    assert identity.name() == name
    assert identity.code() is ck.residue_code(name)
    info = identity.info()
    assert info.code() is identity.code()
    assert identity.is_tabulated() == info.found()
    del identity
    gc.collect()
    assert info.code() is ck.residue_code(name)

def test_checked_residue_info_returns_optional_value():
    index = ck.find_residue_info_index("ALA")
    assert ck.residue_info_checked(index).code() is ck.ResidueCode.ALA
    assert ck.residue_info_checked(1_000_000) is None

def test_structure_hierarchy_snapshot_cow_and_lifetime():
    structure=ck.BioStructure.from_pdb(PDB)
    assert (structure.num_models(),structure.num_chains(),structure.num_residues(),structure.num_atoms())==(1,2,2,2)
    model=structure[-1];chain=model.chains()[0];residue=chain.residues()[0];atom=residue.atoms()[0]
    original=atom.position()
    structure.translate_((2.,0.,0.))
    assert atom.position()==original
    assert structure.atoms()[0].position()==(3.,2.,3.)
    with pytest.raises(IndexError):structure[1]
    with pytest.raises(IndexError):structure[-2]
    with pytest.raises(AttributeError):atom.fake_field=1
    del structure,model,chain,residue;gc.collect()
    assert atom.position()==original
    assert atom.source().serial()==1

def test_protein_filtered_views_and_value_transform():
    protein=ck.Protein.from_pdb(PDB)
    assert (protein.num_models(),protein.num_chains(),protein.num_residues(),protein.num_atoms())==(1,1,1,1)
    chain=protein[-1];residue=chain.residues()[0];atom=residue.atoms()[0]
    assert residue.code() is ck.ResidueCode.ALA
    assert residue.is_standard()
    changed=protein.with_translated_coordinates((1.,2.,3.))
    assert atom.position()==(1.,2.,3.)
    assert changed.atoms()[0].position()==(2.,4.,6.)
    protein.translate_((3.,0.,0.))
    assert atom.position()==(1.,2.,3.)
    with pytest.raises(IndexError):protein[1]
    del protein,chain,residue;gc.collect()
    assert atom.name()=="CA"

def test_selection_projection_retains_source_order():
    structure=ck.BioStructure.from_pdb(PDB)
    selection=ck.BioSelection.from_cid("/1/A/1")
    assert structure.selected_atom_ids(selection)==[0]
    selected=structure.with_selection(selection)
    assert selected.num_atoms()==1
    assert structure.num_atoms()==2
    assert selected.atoms()[0].source().serial()==1

def test_structure_optional_queries_and_atom_id_boundaries():
    structure = ck.BioStructure.from_pdb(PDB)
    assert structure.has_origx() is False
    assert structure.origx().approx(structure.origx(), 0.0)
    assert structure.ncs_oper_identity_id() is None
    assert structure.resolution() == 0.0
    assert structure.ter_status() == 0
    assert structure.atom_position(0) == list(structure.atoms()[0].position())
    assert structure.atom_position(2**32 - 1) is None
    atoms = structure.residue_atoms(0)
    assert [atom.source().serial() for atom in atoms] == [1]
    assert structure.residue_atoms(2**32 - 1) is None
    assert structure.find_entity("missing") is None
    assert structure.find_entity_of_subchain("missing") is None
    structure.translate_((2.0, 0.0, 0.0))
    del structure
    gc.collect()
    assert atoms[0].position() == (1.0, 2.0, 3.0)

def test_crystal_and_transform_snapshot_lifetimes():
    structure = ck.BioStructure.from_pdb(
        "CRYST1   10.000   11.000   12.000  90.00  90.00  90.00 P 1           1\n" + PDB
    )
    crystal = structure.crystal()
    transform = structure.origx()
    del structure
    gc.collect()
    assert crystal.space_group_number() == 1
    assert transform.approx(transform, 0.0)
    assert not transform.approx(transform, float("nan"))

def test_coordinate_block_snapshot_survives_mutation_and_owner_drop():
    structure = ck.BioStructure.from_pdb(PDB)
    assert structure.validate() is None
    coordinates = structure.coordinates()
    original = coordinates.positions
    assert len(coordinates) == structure.num_atoms()
    assert original == [list(atom.position()) for atom in structure.atoms()]
    original[0][0] = -100.0
    assert coordinates.positions[0] == [1.0, 2.0, 3.0]
    structure.translate_((2.0, 0.0, 0.0))
    assert structure.validate() is None
    assert structure.coordinates().positions[0] == [3.0, 2.0, 3.0]
    del structure
    gc.collect()
    assert coordinates.positions[0] == [1.0, 2.0, 3.0]
    with pytest.raises(AttributeError):
        coordinates.positions = []

@pytest.mark.parametrize("raw", [b"", b"CA", b" CA "])
def test_atom_name_preserves_ascii_bytes(raw):
    name = ck.AtomName.from_ascii(raw)
    assert name.as_bytes() == raw
    assert name.as_str() == raw.decode("ascii")
    with pytest.raises(AttributeError):
        name.fake_field = 1

@pytest.mark.parametrize("raw", [b"ABCDE", b"\xff"])
def test_atom_name_rejects_invalid_bytes(raw):
    assert ck.AtomName.from_ascii(raw) is None

def test_atom_lookup_element_and_snapshot_boundaries():
    structure = ck.BioStructure.from_pdb(PDB)
    name = ck.AtomName.from_ascii(b" CA ")
    atom_id, atom = structure.find_atom(0, name, ck.AltLocRequest.Any, ck.Element.C)
    assert atom_id == 0
    assert atom.source().serial() == 1
    assert structure.find_atom(0, name, ck.AltLocRequest.Any, ck.Element.O) is None
    assert structure.find_atom(0, ck.AtomName.from_ascii(b"CA"), ck.AltLocRequest.Any, None) is None
    assert structure.find_atom(2**32 - 1, name, ck.AltLocRequest.Any, None) is None
    assert structure.atom_by_altloc(0, name, None)[0] == atom_id
    structure.translate_((2.0, 0.0, 0.0))
    del structure
    gc.collect()
    assert atom.position() == (1.0, 2.0, 3.0)

def test_atom_lookup_exact_altloc_and_blank_compatibility():
    alternate_a = PDB.splitlines()[0][:16] + "A" + PDB.splitlines()[0][17:]
    alternate_b = PDB.splitlines()[0][:16] + "B" + PDB.splitlines()[0][17:]
    structure = ck.BioStructure.from_pdb(alternate_a + "\n" + alternate_b + "\nEND\n")
    name = ck.AtomName.from_ascii(b" CA ")
    label = ck.AltLocLabel(ord("B"))
    assert label.value() == ord("B")
    assert structure.find_atom(0, name, ck.AltLocRequest.Any, None)[0] == 0
    assert structure.find_atom(0, name, ck.AltLocRequest.Exact(label), None)[0] == 1
    assert structure.atom_by_altloc(0, name, label)[0] == 1
    assert structure.find_atom(0, name, ck.AltLocRequest.Exact(None), None) is None
    blank = ck.BioStructure.from_pdb(PDB)
    assert blank.find_atom(0, name, ck.AltLocRequest.Exact(label), None)[0] == 0
    with pytest.raises(ck.BioStructureError) as missing:
        blank.atom_by_altloc(0, name, label)
    assert missing.value.domain == "bio"
    assert missing.value.kind == "AtomNotFound"

def test_atom_lookup_errors_retain_payload_and_structure():
    structure = ck.BioStructure.from_pdb(PDB)
    original = structure.coordinates().positions
    with pytest.raises(ck.BioStructureError) as invalid:
        structure.atom_by_altloc(2**32 - 1, ck.AtomName.from_ascii(b" CA "), None)
    error = invalid.value
    assert (error.domain, error.kind) == ("bio", "RowReferenceOutOfBounds")
    assert (error.table, error.index, error.table_len) == ("residues", 2**32 - 1, 2)
    with pytest.raises(ck.BioStructureError) as missing:
        structure.atom_by_altloc(0, ck.AtomName.from_ascii(b" N  "), None)
    assert (missing.value.domain, missing.value.kind) == ("bio", "AtomNotFound")
    assert structure.coordinates().positions == original
    assert structure.validate() is None
