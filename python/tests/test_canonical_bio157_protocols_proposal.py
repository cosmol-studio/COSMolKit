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
