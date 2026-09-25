use cosmolkit_model::{
    Atom, AtomId, AtomPdbResidueInfo, AtomPropertyError, AtomSpec, ChiralTag, Element,
    Hybridization,
};

fn full_residue() -> AtomPdbResidueInfo {
    AtomPdbResidueInfo::new(" CA ", 42, "ALA", 7, "A", true)
        .with_alt_loc("B")
        .with_insertion_code("I")
        .with_occupancy(0.75)
        .with_temp_factor(12.5)
        .with_secondary_structure(3)
        .with_segment_number(9)
        .with_monomer_class("LGRP")
}

fn full_spec() -> AtomSpec {
    AtomSpec::new(Element::C)
        .with_element(Element::N)
        .with_formal_charge(-1)
        .with_explicit_hydrogens(2)
        .with_chiral_tag(ChiralTag::TetrahedralCw)
        .with_chiral_permutation(4)
        .with_unknown_stereo(true)
        .with_mol_parity(1)
        .with_mol_inversion_flag(2)
        .with_implicit_hydrogen(true)
        .with_tracked_isotopic_hydrogens(vec![2, 3])
        .with_aromatic(true)
        .with_isotope(15)
        .with_atom_map(8)
        .with_no_implicit(true)
        .with_radical_electrons(1)
        .with_hybridization(Hybridization::Sp2)
        .with_pdb_residue_info(full_residue())
        .with_prop("ordinary", "kept")
        .unwrap()
        .with_computed_prop("computed", "cached")
        .unwrap()
}

#[test]
fn atom_id_is_a_stable_ordered_display_value() {
    let first = AtomId::new(3);
    let second = AtomId::new(8);
    assert_eq!(first.index(), 3);
    assert!(first < second);
    assert_eq!(first.to_string(), "3");
}

#[test]
fn zero_isotope_is_absent_at_every_atom_construction_and_update_boundary() {
    // RDKit Atom::getIsotope() uses zero for an unspecified isotope. The
    // detached model projects that sentinel to None, including explicit zero.
    let unspecified = AtomSpec::new(Element::C);
    let zero = AtomSpec::new(Element::C).with_isotope(13).with_isotope(0);
    assert_eq!(zero.isotope(), None);
    assert_eq!(zero, unspecified);
    let mut atom = Atom::from_spec(AtomId::new(0), zero);
    assert_eq!(atom.isotope(), None);
    for isotope in [1, 13, u16::MAX] {
        atom.set_isotope(Some(isotope));
        assert_eq!(atom.isotope(), Some(isotope));
        atom.set_isotope(Some(0));
        assert_eq!(atom.isotope(), None);
        assert_eq!(atom, Atom::from_spec(AtomId::new(0), unspecified.clone()));
    }
    atom.set_isotope(None);
    assert_eq!(atom.isotope(), None);
}

#[test]
fn pdb_residue_info_covers_source_defaults_and_all_fields() {
    let default = AtomPdbResidueInfo::default();
    assert_eq!(default.atom_name(), "");
    assert_eq!(default.serial_number(), 0);
    assert_eq!(default.alt_loc(), "");
    assert_eq!(default.residue_name(), "");
    assert_eq!(default.residue_number(), 0);
    assert_eq!(default.chain_id(), "");
    assert_eq!(default.insertion_code(), "");
    assert_eq!(default.occupancy(), 1.0);
    assert_eq!(default.temp_factor(), 0.0);
    assert!(!default.is_hetero_atom());
    assert_eq!(default.secondary_structure(), 0);
    assert_eq!(default.segment_number(), 0);
    assert_eq!(default.monomer_class(), "");

    let full = full_residue();
    assert_eq!(full.atom_name(), " CA ");
    assert_eq!(full.serial_number(), 42);
    assert_eq!(full.alt_loc(), "B");
    assert_eq!(full.residue_name(), "ALA");
    assert_eq!(full.residue_number(), 7);
    assert_eq!(full.chain_id(), "A");
    assert_eq!(full.insertion_code(), "I");
    assert_eq!(full.occupancy(), 0.75);
    assert_eq!(full.temp_factor(), 12.5);
    assert!(full.is_hetero_atom());
    assert_eq!(full.secondary_structure(), 3);
    assert_eq!(full.segment_number(), 9);
    assert_eq!(full.monomer_class(), "LGRP");
    assert_eq!(full, full.clone());

    let positive_zero = AtomPdbResidueInfo::default().with_temp_factor(0.0);
    let negative_zero = AtomPdbResidueInfo::default().with_temp_factor(-0.0);
    assert_ne!(positive_zero, negative_zero);
}

#[test]
fn atom_spec_covers_every_default_builder_getter_and_optional_reset() {
    let default = AtomSpec::new(Element::C);
    assert_eq!(default.element(), Element::C);
    assert_eq!(default.formal_charge(), 0);
    assert_eq!(default.explicit_hydrogens(), 0);
    assert_eq!(default.chiral_tag(), ChiralTag::Unspecified);
    assert_eq!(default.chiral_permutation(), None);
    assert!(!default.unknown_stereo());
    assert_eq!(default.mol_parity(), None);
    assert_eq!(default.mol_inversion_flag(), None);
    assert!(!default.implicit_hydrogen());
    assert!(default.tracked_isotopic_hydrogens().is_empty());
    assert!(!default.is_aromatic());
    assert_eq!(default.isotope(), None);
    assert_eq!(default.atom_map(), None);
    assert!(!default.no_implicit());
    assert_eq!(default.radical_electrons(), 0);
    assert_eq!(default.hybridization(), Hybridization::Unspecified);
    assert!(default.props().is_empty());
    assert!(default.computed_prop_names().is_empty());
    assert_eq!(default.pdb_residue_info(), None);

    let full = full_spec();
    assert_eq!(full.element(), Element::N);
    assert_eq!(full.formal_charge(), -1);
    assert_eq!(full.explicit_hydrogens(), 2);
    assert_eq!(full.chiral_tag(), ChiralTag::TetrahedralCw);
    assert_eq!(full.chiral_permutation(), Some(4));
    assert!(full.unknown_stereo());
    assert_eq!(full.mol_parity(), Some(1));
    assert_eq!(full.mol_inversion_flag(), Some(2));
    assert!(full.implicit_hydrogen());
    assert_eq!(full.tracked_isotopic_hydrogens(), &[2, 3]);
    assert!(full.is_aromatic());
    assert_eq!(full.isotope(), Some(15));
    assert_eq!(full.atom_map(), Some(8));
    assert!(full.no_implicit());
    assert_eq!(full.radical_electrons(), 1);
    assert_eq!(full.hybridization(), Hybridization::Sp2);
    assert_eq!(full.prop("ordinary"), Some("kept"));
    assert!(full.is_prop_computed("computed"));
    assert_eq!(full.pdb_residue_info(), Some(&full_residue()));

    let reset = full
        .without_chiral_permutation()
        .without_mol_parity()
        .without_mol_inversion_flag()
        .without_tracked_isotopic_hydrogens()
        .without_isotope()
        .without_atom_map()
        .without_pdb_residue_info();
    assert_eq!(reset.chiral_permutation(), None);
    assert_eq!(reset.mol_parity(), None);
    assert_eq!(reset.mol_inversion_flag(), None);
    assert!(reset.tracked_isotopic_hydrogens().is_empty());
    assert_eq!(reset.isotope(), None);
    assert_eq!(reset.atom_map(), None);
    assert_eq!(reset.pdb_residue_info(), None);
}

#[test]
fn atom_from_spec_preserves_every_fact_and_detached_setters_cover_both_states() {
    let mut atom = Atom::from_spec(AtomId::new(5), full_spec());
    assert_eq!(atom.id(), AtomId::new(5));
    assert_eq!(atom.element(), Element::N);
    assert_eq!(atom.atomic_number(), 7);
    assert_eq!(atom.formal_charge(), -1);
    assert_eq!(atom.explicit_hydrogens(), 2);
    assert_eq!(atom.chiral_tag(), ChiralTag::TetrahedralCw);
    assert_eq!(atom.chiral_permutation(), Some(4));
    assert!(atom.unknown_stereo());
    assert_eq!(atom.mol_parity(), Some(1));
    assert_eq!(atom.mol_inversion_flag(), Some(2));
    assert!(atom.implicit_hydrogen());
    assert_eq!(atom.tracked_isotopic_hydrogens(), &[2, 3]);
    assert!(atom.is_aromatic());
    assert_eq!(atom.isotope(), Some(15));
    assert_eq!(atom.atom_map(), Some(8));
    assert!(atom.no_implicit());
    assert_eq!(atom.radical_electrons(), 1);
    assert_eq!(atom.hybridization(), Hybridization::Sp2);
    assert_eq!(atom.prop("ordinary"), Some("kept"));
    assert!(atom.is_prop_computed("computed"));
    assert_eq!(atom.pdb_residue_info(), Some(&full_residue()));

    atom = atom.with_id(AtomId::new(1));
    atom.set_element(Element::O);
    atom.set_formal_charge(1);
    atom.set_explicit_hydrogens(0);
    atom.set_chiral_tag(ChiralTag::Unspecified);
    atom.set_chiral_permutation(None);
    atom.set_unknown_stereo(false);
    atom.set_mol_parity(None);
    atom.set_mol_inversion_flag(None);
    atom.set_implicit_hydrogen(false);
    atom.set_tracked_isotopic_hydrogens(Vec::new());
    atom.set_aromatic(false);
    atom.set_isotope(None);
    atom.set_atom_map(None);
    atom.set_no_implicit(false);
    atom.set_radical_electrons(0);
    atom.set_hybridization(Hybridization::Unspecified);
    atom.set_pdb_residue_info(None);

    assert_eq!(atom.id(), AtomId::new(1));
    assert_eq!(atom.element(), Element::O);
    assert_eq!(atom.formal_charge(), 1);
    assert_eq!(atom.explicit_hydrogens(), 0);
    assert_eq!(atom.chiral_tag(), ChiralTag::Unspecified);
    assert_eq!(atom.chiral_permutation(), None);
    assert!(!atom.unknown_stereo());
    assert_eq!(atom.mol_parity(), None);
    assert_eq!(atom.mol_inversion_flag(), None);
    assert!(!atom.implicit_hydrogen());
    assert!(atom.tracked_isotopic_hydrogens().is_empty());
    assert!(!atom.is_aromatic());
    assert_eq!(atom.isotope(), None);
    assert_eq!(atom.atom_map(), None);
    assert!(!atom.no_implicit());
    assert_eq!(atom.radical_electrons(), 0);
    assert_eq!(atom.hybridization(), Hybridization::Unspecified);
    assert_eq!(atom.pdb_residue_info(), None);
}

#[test]
fn checked_atom_properties_cover_empty_overwrite_membership_and_clear() {
    assert_eq!(
        AtomSpec::new(Element::C).with_prop("", "value"),
        Err(AtomPropertyError::EmptyKey)
    );
    assert_eq!(
        AtomSpec::new(Element::C).with_computed_prop("", "value"),
        Err(AtomPropertyError::EmptyKey)
    );

    let spec = AtomSpec::new(Element::C)
        .with_computed_prop("cache", "first")
        .unwrap()
        .with_computed_prop("cache", "second")
        .unwrap()
        .with_prop("ordinary", "first")
        .unwrap()
        .with_prop("ordinary", "second")
        .unwrap();
    assert_eq!(spec.prop("cache"), Some("second"));
    assert_eq!(spec.prop("ordinary"), Some("second"));
    assert_eq!(spec.computed_prop_names().len(), 1);

    let mut atom = Atom::from_spec(AtomId::new(0), spec);
    assert_eq!(atom.set_prop("", "value"), Err(AtomPropertyError::EmptyKey));
    assert_eq!(
        atom.set_computed_prop("", "value"),
        Err(AtomPropertyError::EmptyKey)
    );
    atom.set_prop("cache", "ordinary overwrite").unwrap();
    assert!(atom.is_prop_computed("cache"));
    atom.set_computed_prop("cache", "computed overwrite")
        .unwrap();
    assert_eq!(atom.computed_prop_names().len(), 1);

    atom.clear_prop("missing");
    atom.clear_prop("cache");
    assert_eq!(atom.prop("cache"), None);
    assert!(!atom.is_prop_computed("cache"));

    atom.set_computed_prop("temporary", "gone").unwrap();
    atom.clear_computed_props();
    assert_eq!(atom.prop("temporary"), None);
    assert_eq!(atom.prop("ordinary"), Some("second"));
    assert!(atom.computed_prop_names().is_empty());
}
