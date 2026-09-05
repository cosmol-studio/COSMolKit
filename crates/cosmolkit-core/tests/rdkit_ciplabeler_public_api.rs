use std::str::FromStr;

use cosmolkit_core::{
    AtomSpec, BondOrder, BondSpec, CipDescriptor, CipDescriptorError, Element, Molecule, MoleculeBuilder,
};

const DESCRIPTORS: [(CipDescriptor, &str); 10] = [
    (CipDescriptor::R, "R"),
    (CipDescriptor::S, "S"),
    (CipDescriptor::LowerR, "r"),
    (CipDescriptor::LowerS, "s"),
    (CipDescriptor::E, "E"),
    (CipDescriptor::Z, "Z"),
    (CipDescriptor::M, "M"),
    (CipDescriptor::P, "P"),
    (CipDescriptor::LowerM, "m"),
    (CipDescriptor::LowerP, "p"),
];

fn molecule_with_atom_descriptor(value: Option<&str>) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let atom = value.map_or_else(
        || AtomSpec::new(Element::C),
        |value| AtomSpec::new(Element::C).with_prop("_CIPCode", value),
    );
    builder.add_atom(atom);
    builder.build().expect("single atom must build")
}

fn molecule_with_bond_descriptor(value: Option<&str>) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let begin = builder.add_atom(AtomSpec::new(Element::C));
    let end = builder.add_atom(AtomSpec::new(Element::C));
    let bond = value.map_or_else(
        || BondSpec::new(begin, end, BondOrder::Single),
        |value| BondSpec::new(begin, end, BondOrder::Single).with_prop("_CIPCode", value),
    );
    builder.add_bond(bond).expect("single bond must build");
    builder.build().expect("two-atom molecule must build")
}

#[test]
fn typed_modern_cip_descriptors_roundtrip_stable_source_spellings() {
    for (descriptor, spelling) in DESCRIPTORS {
        assert_eq!(descriptor.as_str(), spelling);
        assert_eq!(descriptor.to_string(), spelling);
        assert_eq!(CipDescriptor::from_str(spelling), Ok(descriptor));
    }
}

#[test]
fn atom_and_bond_queries_return_every_supported_modern_descriptor() {
    for (descriptor, spelling) in DESCRIPTORS {
        let atom_molecule = molecule_with_atom_descriptor(Some(spelling));
        assert_eq!(atom_molecule.atoms()[0].cip_descriptor(), Ok(Some(descriptor)));

        let bond_molecule = molecule_with_bond_descriptor(Some(spelling));
        assert_eq!(bond_molecule.bonds()[0].cip_descriptor(), Ok(Some(descriptor)));
    }
}

#[test]
fn atom_and_bond_queries_distinguish_absent_from_malformed_state() {
    let atom_molecule = molecule_with_atom_descriptor(None);
    assert_eq!(atom_molecule.atoms()[0].cip_descriptor(), Ok(None));

    let bond_molecule = molecule_with_bond_descriptor(None);
    assert_eq!(bond_molecule.bonds()[0].cip_descriptor(), Ok(None));

    let malformed = "not-a-modern-cip-descriptor";
    let expected = Err(CipDescriptorError::InvalidStoredDescriptor {
        value: malformed.to_owned(),
    });
    assert_eq!(
        molecule_with_atom_descriptor(Some(malformed)).atoms()[0].cip_descriptor(),
        expected
    );
    assert_eq!(
        molecule_with_bond_descriptor(Some(malformed)).bonds()[0].cip_descriptor(),
        expected
    );
}

#[test]
fn unsupported_non_tetrahedral_codes_are_not_exposed_as_assigned_descriptors() {
    for value in ["SP_4", "TBPY_5", "OC_6"] {
        assert_eq!(
            CipDescriptor::from_str(value),
            Err(CipDescriptorError::InvalidStoredDescriptor {
                value: value.to_owned(),
            })
        );
        assert_eq!(
            molecule_with_atom_descriptor(Some(value)).atoms()[0].cip_descriptor(),
            Err(CipDescriptorError::InvalidStoredDescriptor {
                value: value.to_owned(),
            })
        );
        assert_eq!(
            molecule_with_bond_descriptor(Some(value)).bonds()[0].cip_descriptor(),
            Err(CipDescriptorError::InvalidStoredDescriptor {
                value: value.to_owned(),
            })
        );
    }
}
