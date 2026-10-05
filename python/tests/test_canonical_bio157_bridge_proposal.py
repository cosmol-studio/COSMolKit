"""Source-backed BIO molecular-conversion proposals, pending p1/ROOT review."""
import itertools
import pytest
import cosmolkit as ck

def atom(serial,name,element,x,y=0.,z=0.,alt=" ",res="LIG"):
    return f"HETATM{serial:5d} {name:4s}{alt}{res:3s} A   1    {x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{20.:6.2f}          {element:>2s}  \n"

ETHANE=atom(1," C1 ","C",0.)+atom(2," C2 ","C",1.4)+atom(3," H1 ","H",-1.)+"CONECT    1    2    3\nEND\n"

@pytest.mark.parametrize("sanitize,remove_hs,proximity",list(itertools.product((False,True),repeat=3)))
def test_pipeline_controls_and_checked_constructor(sanitize,remove_hs,proximity):
    source=ck.BioStructure.from_pdb(ETHANE)
    params=ck.BioMoleculeParams(sanitize=sanitize,remove_hs=remove_hs,proximity_bonding=proximity)
    molecule=source.to_molecule_with_params(params)
    assert isinstance(molecule,ck.Molecule)
    assert molecule.num_atoms()==(2 if sanitize and remove_hs else 3)
    assert molecule.num_bonds()==(1 if sanitize and remove_hs else 2)
    assert source.num_atoms()==3
    explicit=source.to_molecule(sanitize=sanitize,remove_hs=remove_hs,proximity_bonding=proximity)
    assert explicit.num_atoms()==molecule.num_atoms()
    assert explicit.num_bonds()==molecule.num_bonds()

def test_default_bridge_delegation():
    source=ck.BioStructure.from_pdb(ETHANE)
    assert source.to_molecule().num_atoms()==2
    assert source.to_molecule_with_params(ck.BioMoleculeParams()).num_atoms()==2

@pytest.mark.parametrize("flavor,expected",[(0,1),(1,5)])
def test_source_atom_filters_and_original_coordinate_bytes(flavor,expected):
    text=(atom(1," C1 ","C",0.)+atom(2," C2 ","C",10.,alt="B")+
          atom(3," C3 ","C",20.,res="DUM")+atom(4," Q1 ","C",30.)+
          atom(5," C5 ","C",9999.,9999.,9999.)+"END\n")
    source=ck.BioStructure.from_pdb(text)
    assert source.num_atoms()==5
    result=source.to_molecule(sanitize=False,remove_hs=False,proximity_bonding=False,flavor=flavor)
    assert result.num_atoms()==expected
    assert result.num_bonds()==0

def test_repeated_conect_targets_reuse_existing_bond_order_owner():
    text=atom(1," C1 ","C",0.)+atom(2," C2 ","C",1.4)+"CONECT    1    2    2\nEND\n"
    molecule=ck.BioStructure.from_pdb(text).to_molecule(remove_hs=False,proximity_bonding=False)
    assert molecule.to_smiles()=="C=C"

def test_all_hierarchy_models_become_one_detached_candidate():
    text="MODEL        1\n"+atom(1," C1 ","C",0.)+"ENDMDL\nMODEL        2\n"+atom(2," C2 ","C",10.)+"ENDMDL\nEND\n"
    source=ck.BioStructure.from_pdb(text)
    result=source.to_molecule(sanitize=False,remove_hs=False,proximity_bonding=False)
    assert source.num_models()==2
    assert result.num_atoms()==2

def test_empty_structure_checked_conversion():
    source=ck.BioStructure.from_pdb("END\n")
    assert source.to_molecule().num_atoms()==0

def test_real_sanitize_failure_preserves_source_and_domain_error():
    text="".join(atom(i," C  ","C",10.*i) for i in range(1,7))+"CONECT    1    2    3    4    5\nCONECT    1    6\nEND\n"
    source=ck.BioStructure.from_pdb(text)
    before=[row.position() for row in source.atoms()]
    with pytest.raises(ck.BioMoleculeError) as observed:
        source.to_molecule(proximity_bonding=False,remove_hs=False)
    assert observed.value.domain=="bio"
    assert observed.value.kind=="Conversion"
    assert observed.value.__cause__ is not None
    assert [row.position() for row in source.atoms()]==before
