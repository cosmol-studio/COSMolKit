from __future__ import annotations

import enum
import pytest
import cosmolkit as ck


@pytest.mark.parametrize("name,count", [("BondOrder",22),("ChiralTag",9),("BondDirection",7),("BondStereo",8)])
def test_complete_dynamic_vocabulary_transport(name,count):
    cls=getattr(ck,name)
    assert issubclass(cls,enum.IntEnum)
    assert [int(v) for v in cls]==list(range(count))
    for member in cls:
        assert cls(int(member)) is member
        assert getattr(cls,member.name) is member
    with pytest.raises(ValueError):cls(count)
    with pytest.raises(ValueError):cls(-1)


def test_atom_scalar_and_core_metadata_records_have_one_canonical_context():
    molecule=ck.Molecule.from_smiles("CCO")
    atoms=molecule.atoms()
    metadata=molecule.atom_metadata()
    assert len(molecule)==3
    assert molecule.num_bonds()==2
    assert [a.id() for a in atoms]==[0,1,2]
    assert [a.atomic_number() for a in atoms]==[6,6,8]
    assert [a.element() for a in atoms]==[ck.Element.C,ck.Element.C,ck.Element.O]
    assert [a.degree() for a in atoms]==[1,2,1]
    assert [a.explicit_hydrogens() for a in atoms]==[0,0,0]
    assert [a.implicit_hydrogens() for a in atoms]==[3,2,1]
    assert [a.total_hydrogens() for a in atoms]==[3,2,1]
    assert [a.explicit_valence() for a in atoms]==[1,2,1]
    assert [a.total_valence() for a in atoms]==[4,4,2]
    assert [a.formal_charge() for a in atoms]==[0,0,0]
    assert [a.isotope() for a in atoms]==[None,None,None]
    assert [a.atom_map() for a in atoms]==[None,None,None]
    assert [a.radical_electrons() for a in atoms]==[0,0,0]
    assert len(atoms) == len(metadata)
    for a,m in zip(atoms,metadata):
        for field in ("degree","explicit_valence","implicit_hydrogens","total_hydrogens","total_valence"):
            assert getattr(a,field)()==getattr(m,field)()
        assert a.no_implicit() is False
        assert a.is_aromatic() is False
        assert a.chiral_tag()==ck.ChiralTag.CHI_UNSPECIFIED
        assert a.chiral_tag_code()==0
        assert a.chiral_tag_name()=="CHI_UNSPECIFIED"
        assert a.cip_neighbor_order() is None
        assert a.cip_descriptor() is None
        assert f"id={a.id()}" in repr(a)
    assert repr(molecule)=="Molecule(num_atoms=3, num_bonds=2)"
    assert len(ck.Molecule.new())==0
    assert ck.Molecule.new().atoms()==[]
    assert ck.Molecule.new().bonds()==[]
    assert ck.Molecule.new().atom_metadata()==[]


def test_atom_isotope_map_aromatic_charge_explicit_hydrogen_noimplicit():
    # Both oversized source-number spellings were actually rejected by the
    # pinned source-built native reference; preserve them as error conditions.
    with pytest.raises(ck.SmilesError):ck.Molecule.from_smiles("[13CH3:4294967295]")
    with pytest.raises(ck.SmilesError):ck.Molecule.from_smiles("[13CH3:2147483647]")
    atom=ck.Molecule.from_smiles("[13CH3:7]").atoms()[0]
    assert atom.isotope()==13
    assert atom.atom_map()==7
    # Canonical detached UInt storage is independent of SMILES's number lexer.
    builder=ck.MoleculeBuilder.new()
    index=builder.add_atom(ck.AtomSpec(ck.Element.C).with_atom_map(4294967295))
    assert builder.build().atoms()[index].atom_map()==4294967295
    assert atom.explicit_hydrogens()==3
    assert atom.no_implicit() is True
    assert atom.radical_electrons()==1
    assert atom.degree()==0
    assert atom.total_hydrogens()==3
    positive=ck.Molecule.from_smiles("[NH4+]").atoms()[0]
    assert positive.formal_charge()==1
    assert positive.explicit_hydrogens()==4
    assert positive.implicit_hydrogens()==0
    assert positive.total_hydrogens()==4
    assert all(a.is_aromatic() for a in ck.Molecule.from_smiles("c1ccccc1").atoms())


def test_bond_scalar_stereo_and_dynamic_vocabulary_records():
    bonds=ck.Molecule.from_smiles("F/C=C/F").bonds()
    assert [b.id() for b in bonds]==[0,1,2]
    assert [(b.begin(),b.end()) for b in bonds]==[(0,1),(1,2),(2,3)]
    assert [b.order() for b in bonds]==[ck.BondOrder.SINGLE,ck.BondOrder.DOUBLE,ck.BondOrder.SINGLE]
    for b in bonds:
        assert b.order_code()==int(b.order())
        assert b.order_name()==b.order().name
        assert b.direction_code()==int(b.direction())
        assert b.direction_name()==b.direction().name
        assert b.stereo_code()==int(b.stereo())
        assert b.stereo_name()==b.stereo().name
        assert b.is_aromatic() is False
        assert b.cip_neighbor_order() is None
        assert b.cip_descriptor() is None
        assert f"id={b.id()}" in repr(b)
    assert bonds[1].stereo() in (ck.BondStereo.STEREOE,ck.BondStereo.STEREOTRANS)
    assert bonds[1].stereo_atoms()==[0,3]
    assert bonds[0].stereo_atoms() is None
    assert all(b.order()==ck.BondOrder.AROMATIC and b.is_aromatic() for b in ck.Molecule.from_smiles("c1ccccc1").bonds())


def test_read_values_and_cip_options_are_immutable():
    molecule=ck.Molecule.from_smiles("CC")
    for value in (molecule.atoms()[0],molecule.bonds()[0],molecule.atom_metadata()[0]):
        with pytest.raises(AttributeError):value.id=8
    options=ck.CipLabelOptions()
    assert options.atoms is None and options.bonds is None
    assert options.max_recursive_iterations==0
    options=ck.CipLabelOptions(atoms=[],bonds=[1,0,1],max_recursive_iterations=3)
    assert options.atoms==[] and options.bonds==[1,0,1]
    assert options.max_recursive_iterations==3
    copy=options.bonds
    copy.append(2)
    assert options.bonds==[1,0,1]
    with pytest.raises(AttributeError):options.atoms=[0]
    with pytest.raises(OverflowError):ck.CipLabelOptions(atoms=[-1])
    with pytest.raises(OverflowError):ck.CipLabelOptions(max_recursive_iterations=-1)


def test_sanitize_value_and_generated_inplace_are_equivalent_and_atomic():
    molecule=ck.Molecule.from_smiles("CCO")
    before=molecule.to_smiles()
    sanitized=molecule.sanitize()
    assert sanitized is not molecule
    assert sanitized.to_smiles()==before
    assert molecule.to_smiles()==before
    assert molecule.sanitize_() is None
    assert molecule.to_smiles()==before
