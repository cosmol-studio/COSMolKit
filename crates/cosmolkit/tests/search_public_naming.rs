#![cfg(feature = "cap-search")]

use cosmolkit::{
    BINDING_CONTRACT, BindingDefault, BindingKind, BindingOwner, QueryGraph, SmartsParseError,
    SmartsParseParams,
};

#[test]
fn flat_functions_and_bound_query_factories_share_the_parser() {
    let _: fn(&str) -> Result<QueryGraph, SmartsParseError> = cosmolkit::parse_smarts;
    let _: fn(&str) -> Result<QueryGraph, SmartsParseError> = cosmolkit::search::from_smarts;
    for (text, atoms, bonds) in [("C(N)O", 3, 2), ("C1CC1", 3, 3), ("N<-C", 2, 1)] {
        let flat = cosmolkit::parse_smarts(text).unwrap();
        let factory = cosmolkit::search::from_smarts(text).unwrap();
        for query in [&flat, &factory] {
            assert_eq!((query.num_atoms(), query.num_bonds()), (atoms, bonds));
        }
        assert_eq!(
            cosmolkit::write_smarts(&flat, &Default::default()).unwrap(),
            cosmolkit::write_smarts(&factory, &Default::default()).unwrap()
        );
    }
    let params = SmartsParseParams {
        merge_hs: true,
        ..Default::default()
    };
    for query in [
        cosmolkit::parse_smarts_with_params("C[H]", &params).unwrap(),
        cosmolkit::search::from_smarts_with_params("C[H]", &params).unwrap(),
    ] {
        assert_eq!((query.num_atoms(), query.num_bonds()), (1, 0));
    }
    assert!(params.merge_hs);
    for text in ["C(", "[#6", "C |notCX| name"] {
        let flat = cosmolkit::parse_smarts(text).unwrap_err();
        let factory = cosmolkit::search::from_smarts(text).unwrap_err();
        assert_eq!(format!("{flat:?}"), format!("{factory:?}"));
    }
}

#[test]
fn registry_records_both_static_factories_and_the_flat_functions() {
    for (id, name, parameter_names) in [
        ("QueryGraph.from_smarts", "from_smarts", &["text"][..]),
        (
            "QueryGraph.from_smarts_with_params",
            "from_smarts_with_params",
            &["text", "params"][..],
        ),
    ] {
        let row = BINDING_CONTRACT
            .iter()
            .find(|row| row.semantic_id == id)
            .unwrap();
        assert_eq!(row.owner, BindingOwner::Type);
        assert_eq!(row.python_name, name);
        assert_eq!(row.feature, "cap-search");
        let callable = row.callable.unwrap();
        assert_eq!(callable.kind, BindingKind::Static);
        assert_eq!(callable.receiver, None);
        assert_eq!(
            callable
                .parameters
                .iter()
                .map(|p| p.name)
                .collect::<Vec<_>>(),
            parameter_names
        );
        assert!(
            callable
                .parameters
                .iter()
                .all(|p| p.default == BindingDefault::Required)
        );
    }
    for name in [
        "parse_smarts",
        "parse_smarts_with_params",
        "compile_query",
        "write_smarts",
        "write_cx_smarts",
    ] {
        let row = BINDING_CONTRACT
            .iter()
            .find(|row| row.semantic_id == format!("search.{name}"))
            .unwrap();
        assert_eq!(row.owner, BindingOwner::Module);
        assert_eq!(row.python_name, name);
        let path: String = row
            .rust_path
            .chars()
            .filter(|c| !c.is_whitespace())
            .collect();
        assert_eq!(path, format!("crate::{name}"));
    }
}

#[cfg(feature = "cap-smiles")]
#[test]
fn concrete_smarts_preserves_source_atom_attributes_and_rooting() {
    // SmartsWrite.cpp getNonQueryAtomSmarts/getNonQueryBondSmarts:
    // organic carriers use atomic numbers, explicit H appears only with
    // chirality, and traversal uses the requested original atom index.
    for (input, expected) in [
        ("CCO", "[#6]-[#6]-[#8]"),
        ("c1ccccc1", "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"),
        ("[13CH3][NH3+:7]", "[13#6]-[#7+:7]"),
        ("[Mg+2].[Cl-]", "[Mg+2].[#17-]"),
        ("", ""),
    ] {
        let molecule = cosmolkit::Molecule::from_smiles(input).unwrap();
        let before = molecule.to_smiles().unwrap();
        assert_eq!(
            molecule.to_smarts().unwrap().as_bytes(),
            expected.as_bytes(),
            "{input}"
        );
        assert_eq!(
            molecule.to_cx_smarts().unwrap().as_bytes(),
            expected.as_bytes(),
            "{input}"
        );
        assert_eq!(molecule.to_smiles().unwrap(), before);
    }
    let molecule = cosmolkit::Molecule::from_smiles("CCO").unwrap();
    let params = cosmolkit::SmartsWriteParams {
        rooted_at_atom: Some(2),
        ..Default::default()
    };
    assert_eq!(
        molecule.to_smarts_with_params(&params).unwrap().as_bytes(),
        b"[#8]-[#6]-[#6]"
    );
    let params = cosmolkit::SmartsWriteParams {
        rooted_at_atom: Some(3),
        ..Default::default()
    };
    assert!(matches!(
        molecule.to_smarts_with_params(&params),
        Err(cosmolkit::SmartsWriteError::RootedAtomOutOfRange { atom: 3 })
    ));
}

#[cfg(feature = "cap-smiles")]
#[test]
fn concrete_smarts_keeps_nonquery_mapping_and_cx_dative_source_branches() {
    let molecule = cosmolkit::Molecule::from_smiles("[NH3:7]->[Cu+2]").unwrap();
    let params = cosmolkit::SmartsWriteParams {
        include_atom_maps: false,
        ..Default::default()
    };
    // Native getNonQueryAtomSmarts reads molAtomMapNumber unconditionally;
    // MolToCXSmarts disables the dative token and emits its CX extension.
    assert_eq!(
        molecule.to_smarts_with_params(&params).unwrap().as_bytes(),
        b"[#7:7]->[Cu+2]"
    );
    assert_eq!(
        molecule
            .to_cx_smarts_with_params(&params)
            .unwrap()
            .as_bytes(),
        b"[#7:7]-[Cu+2] |C:0.0|"
    );
    let params = cosmolkit::SmartsWriteParams {
        rooted_at_atom: Some(1),
        ..Default::default()
    };
    assert_eq!(
        molecule.to_smarts_with_params(&params).unwrap().as_bytes(),
        b"[Cu+2]<-[#7:7]"
    );
}

#[test]
fn concrete_smarts_registry_exposes_read_only_parameterized_methods() {
    for name in [
        "to_smarts",
        "to_smarts_with_params",
        "to_cx_smarts",
        "to_cx_smarts_with_params",
    ] {
        let row = BINDING_CONTRACT
            .iter()
            .find(|row| row.semantic_id == format!("Molecule.{name}"))
            .unwrap();
        assert_eq!(row.owner, BindingOwner::Molecule);
        assert_eq!(row.python_name, name);
        assert_eq!(row.feature, "cap-search");
        let callable = row.callable.unwrap();
        assert_eq!(callable.kind, BindingKind::Instance);
        assert_eq!(
            callable.parameters.len(),
            usize::from(name.ends_with("_with_params"))
        );
    }
}

#[cfg(feature = "cap-smiles")]
#[test]
fn concrete_smarts_rooted_tetrahedral_order_follows_source_permutation() {
    let molecule = cosmolkit::Molecule::from_smiles("F[C@H](Cl)Br").unwrap();
    // Canon.cpp uses the incoming bond, ring closures, and outgoing bonds as
    // trueOrder. Roots F/Cl/Br give [0,1,2]/[1,0,2]/[2,0,1]; starting at the
    // center adds the source explicit-H inversion. These are fixed source
    // permutation observations, not strings generated by another writer.
    for (root, expected) in [
        (0, "[#9]-[#6@H](-[#17])-[#35]"),
        (1, "[#6@@H](-[#9])(-[#17])-[#35]"),
        (2, "[#17]-[#6@@H](-[#9])-[#35]"),
        (3, "[#35]-[#6@H](-[#9])-[#17]"),
    ] {
        let params = cosmolkit::SmartsWriteParams {
            rooted_at_atom: Some(root),
            ..Default::default()
        };
        assert_eq!(
            molecule.to_smarts_with_params(&params).unwrap().as_bytes(),
            expected.as_bytes(),
            "root={root}"
        );
    }
    let params = cosmolkit::SmartsWriteParams {
        isomeric_smiles: false,
        rooted_at_atom: Some(1),
        ..Default::default()
    };
    assert_eq!(
        molecule.to_smarts_with_params(&params).unwrap().as_bytes(),
        b"[#6](-[#9])(-[#17])-[#35]"
    );
}

#[test]
fn concrete_smarts_keeps_raw_symbol_bytes_and_signed_map_read_errors() {
    use cosmolkit::{
        Atom, AtomId, AtomSpec, CoordinateBlock, Element, MoleculeBuilder, MoleculeProperties,
        PropertyValue, TopologyBlock,
    };
    let mut atom = Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C));
    atom.set_prop(
        "smilesSymbol",
        PropertyValue::String(vec![b'X', 0xff].into()),
    )
    .unwrap();
    atom.set_prop("molAtomMapNumber", PropertyValue::Int(-7))
        .unwrap();
    let topology =
        TopologyBlock::try_from_parts(vec![atom.clone()], vec![], vec![], vec![]).unwrap();
    let molecule = MoleculeBuilder::from_parts(
        topology,
        CoordinateBlock::default(),
        MoleculeProperties::default(),
    )
    .build()
    .unwrap();
    assert_eq!(
        molecule.to_smarts().unwrap().as_bytes(),
        &[b'[', b'X', 0xff, b':', b'-', b'7', b']']
    );
    atom.set_prop("molAtomMapNumber", PropertyValue::UInt(2_147_483_648))
        .unwrap();
    let topology = TopologyBlock::try_from_parts(vec![atom], vec![], vec![], vec![]).unwrap();
    let molecule = MoleculeBuilder::from_parts(
        topology,
        CoordinateBlock::default(),
        MoleculeProperties::default(),
    )
    .build()
    .unwrap();
    assert!(matches!(
        molecule.to_smarts(),
        Err(cosmolkit::SmartsWriteError::CxAtomPropertyInt {
            property: "molAtomMapNumber",
            source: cosmolkit_core::PropertyIntReadError::UnsignedOverflow {
                value: 2_147_483_648
            },
            ..
        })
    ));
}

#[test]
fn concrete_cx_smarts_preserves_first_source_xyz_even_when_not_three_dimensional() {
    use cosmolkit::{
        Atom, AtomId, AtomSpec, Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension,
        Element, MoleculeBuilder, MoleculeProperties, TopologyBlock,
    };
    let topology = TopologyBlock::try_from_parts(
        vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
        vec![],
        vec![],
        vec![],
    )
    .unwrap();
    let coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(8, vec![[9.0, 8.0]])],
        conformers_3d: vec![Conformer3D::new(7, vec![[1.0, 2.0, 3.0]], false)],
        source_coordinate_dim: None,
        source_conformer_order: Some(vec![CoordinateDimension::ThreeD, CoordinateDimension::TwoD]),
    };
    let molecule =
        MoleculeBuilder::from_parts(topology, coordinates, MoleculeProperties::default())
            .build()
            .unwrap();
    assert_eq!(
        molecule.to_cx_smarts().unwrap().as_bytes(),
        b"[#6] |(1,2,)|"
    );
    assert_eq!(
        molecule.conformers_3d()[0].coordinates(),
        &[[1.0, 2.0, 3.0]]
    );
    assert!(!molecule.conformers_3d()[0].is_3d());
}
