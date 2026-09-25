#[path = "../src/stereo.rs"]
mod stereo;

use cosmolkit_smiles::{
    SmilesParseError, SmilesWriteParams, parse_smiles, write_smiles, write_smiles_with_params,
};
use cosmolkit_types::{BondOrder, ChiralTag};

#[test]
fn swap_counting_and_tetrahedral_helpers_cover_source_boundaries() {
    assert_eq!(
        stereo::count_swaps_to_interconvert::<u8>(&[], vec![]),
        Some(0)
    );
    assert_eq!(stereo::count_swaps_to_interconvert(&[1], vec![1]), Some(0));
    assert_eq!(
        stereo::count_swaps_to_interconvert(&[1, 2, 3, 4], vec![4, 3, 2, 1]),
        Some(2)
    );
    assert_eq!(
        stereo::count_swaps_to_interconvert(&[1, 1, 2], vec![1, 2, 1]),
        Some(1)
    );
    assert_eq!(stereo::count_swaps_to_interconvert(&[1], vec![]), None);
    assert_eq!(
        stereo::count_swaps_to_interconvert(&[1, 2], vec![1, 3]),
        None
    );

    assert_eq!(
        stereo::invert_tetrahedral_tag(ChiralTag::TetrahedralCw),
        ChiralTag::TetrahedralCcw
    );
    assert_eq!(
        stereo::invert_tetrahedral_tag(ChiralTag::TetrahedralCcw),
        ChiralTag::TetrahedralCw
    );
    for tag in [
        ChiralTag::Unspecified,
        ChiralTag::Other,
        ChiralTag::Tetrahedral,
        ChiralTag::Allene,
        ChiralTag::SquarePlanar,
        ChiralTag::TrigonalBipyramidal,
        ChiralTag::Octahedral,
    ] {
        assert_eq!(stereo::invert_tetrahedral_tag(tag), tag);
    }

    assert!(stereo::atom_has_fourth_valence(1, false));
    assert!(stereo::atom_has_fourth_valence(0, true));
    assert!(!stereo::atom_has_fourth_valence(0, false));
    assert!(!stereo::atom_has_fourth_valence(2, false));

    assert!(stereo::chiral_atom_needs_tag_inversion(
        3, 1, true, true, 0, true
    ));
    assert!(stereo::chiral_atom_needs_tag_inversion(
        3, 0, false, false, 1, false
    ));
    for arguments in [
        (2, 1, true, false, 1, false),
        (4, 1, true, false, 1, false),
        (3, 2, true, false, 0, false),
        (3, 0, false, true, 1, false),
        (3, 0, false, false, 0, false),
        (3, 0, false, false, 2, false),
        (3, 0, false, false, 1, true),
    ] {
        assert!(!stereo::chiral_atom_needs_tag_inversion(
            arguments.0,
            arguments.1,
            arguments.2,
            arguments.3,
            arguments.4,
            arguments.5,
        ));
    }
}

#[test]
fn bond_order_projection_covers_the_complete_vocabulary() {
    for (order, expected) in [
        (BondOrder::Single, 1.0),
        (BondOrder::Double, 2.0),
        (BondOrder::Triple, 3.0),
        (BondOrder::Quadruple, 4.0),
        (BondOrder::Quintuple, 5.0),
        (BondOrder::Hextuple, 6.0),
        (BondOrder::Aromatic, 1.5),
        (BondOrder::OneAndHalf, 1.5),
        (BondOrder::TwoAndHalf, 2.5),
        (BondOrder::ThreeAndHalf, 3.5),
        (BondOrder::FourAndHalf, 4.5),
        (BondOrder::FiveAndHalf, 5.5),
        (BondOrder::Dative, 1.0),
        (BondOrder::DativeOne, 1.0),
        (BondOrder::DativeLeft, 1.0),
        (BondOrder::DativeRight, 1.0),
        (BondOrder::Ionic, 0.0),
        (BondOrder::Hydrogen, 0.0),
        (BondOrder::ThreeCenter, 0.0),
        (BondOrder::Other, 0.0),
        (BondOrder::Zero, 0.0),
        (BondOrder::Unspecified, 0.0),
    ] {
        assert_eq!(stereo::bond_order_as_double(order), expected, "{order:?}");
    }
}

#[test]
fn implicit_nontetrahedral_neighbors_follow_source_positions() {
    for (tag, maximum) in [
        (ChiralTag::SquarePlanar, 4),
        (ChiralTag::TrigonalBipyramidal, 5),
        (ChiralTag::Octahedral, 6),
    ] {
        assert_eq!(stereo::nontetrahedral_max_neighbors(tag), Some(maximum));
    }
    for tag in [
        ChiralTag::Unspecified,
        ChiralTag::TetrahedralCw,
        ChiralTag::TetrahedralCcw,
        ChiralTag::Other,
        ChiralTag::Tetrahedral,
        ChiralTag::Allene,
    ] {
        assert_eq!(stereo::nontetrahedral_max_neighbors(tag), None);
    }

    let mut first = vec![Some(cosmolkit_model::BondId::new(0))];
    stereo::insert_implicit_nontetrahedral_neighbors(&mut first, ChiralTag::SquarePlanar, true);
    assert_eq!(
        first,
        [None, None, None, Some(cosmolkit_model::BondId::new(0))]
    );

    let mut later = vec![
        Some(cosmolkit_model::BondId::new(0)),
        Some(cosmolkit_model::BondId::new(1)),
    ];
    stereo::insert_implicit_nontetrahedral_neighbors(&mut later, ChiralTag::SquarePlanar, false);
    assert_eq!(
        later,
        [
            Some(cosmolkit_model::BondId::new(0)),
            None,
            None,
            Some(cosmolkit_model::BondId::new(1)),
        ]
    );

    let mut complete = vec![None; 4];
    stereo::insert_implicit_nontetrahedral_neighbors(&mut complete, ChiralTag::SquarePlanar, false);
    assert_eq!(complete, vec![None; 4]);
}

#[test]
fn parser_and_writer_preserve_tetrahedral_storage_and_traversal_order() {
    for (input, atom_index, stored, noncanonical) in [
        (
            "[C@H](F)(Cl)Br",
            0,
            ChiralTag::TetrahedralCw,
            "[C@H](F)(Cl)Br",
        ),
        (
            "[C@@H](F)(Cl)Br",
            0,
            ChiralTag::TetrahedralCcw,
            "[C@@H](F)(Cl)Br",
        ),
        ("F[C@H](Cl)Br", 1, ChiralTag::TetrahedralCcw, "F[C@H](Cl)Br"),
        (
            "F[C@]1(CCO1)Br",
            1,
            ChiralTag::TetrahedralCcw,
            "F[C@@]1(Br)CCO1",
        ),
    ] {
        let record = parse_smiles(input, &Default::default()).unwrap();
        assert_eq!(
            record.topology.atoms[atom_index].chiral_tag(),
            stored,
            "{input}"
        );
        assert_eq!(
            write_smiles_with_params(
                &record,
                &SmilesWriteParams {
                    canonical: false,
                    ..Default::default()
                },
            )
            .unwrap(),
            noncanonical,
            "{input}"
        );
    }
}

#[test]
fn parser_and_writer_roundtrip_sp_tb_oh_typed_permutations() {
    for (input, atom_index, stored, canonical) in [
        ("[Pt@SP1](F)(Cl)(Br)I", 0, 1, "[F][Pt@SP1]([Cl])([Br])[I]"),
        ("F[Pt@SP1](Cl)Br", 1, 2, "[F][Pt@SP1]([Cl])[Br]"),
        ("[P@TB1](F)(Cl)(Br)(I)N", 0, 1, "N[P@TB8](F)(Cl)(Br)I"),
        ("F[P@TB1](Cl)(Br)I", 1, 3, "F[P@TB1](Cl)(Br)I"),
        (
            "[Co@OH1](F)(Cl)(Br)(I)(N)O",
            0,
            1,
            "[NH2][Co@OH9]([OH])([F])([Cl])([Br])[I]",
        ),
        (
            "F[Co@OH1](Cl)(Br)(I)N",
            1,
            3,
            "[NH2][Co@OH24]([F])([Cl])([Br])[I]",
        ),
    ] {
        let record = parse_smiles(input, &Default::default()).unwrap();
        assert_eq!(
            record.topology.atoms[atom_index].chiral_permutation(),
            Some(stored)
        );
        let output = write_smiles(&record).unwrap();
        assert_eq!(output, canonical, "{input}");
        let reparsed = parse_smiles(&output, &Default::default()).unwrap();
        let center = reparsed
            .topology
            .atoms
            .iter()
            .find(|atom| atom.chiral_tag() == record.topology.atoms[atom_index].chiral_tag())
            .expect("stereo center");
        assert!(center.chiral_permutation().is_some());
    }
}

#[test]
fn parser_rejects_invalid_source_permutations_without_partial_output() {
    for input in ["[C@TH0]", "[C@TH3]", "[C@SP4]", "[C@TB21]", "[C@OH31]"] {
        assert!(matches!(
            parse_smiles(input, &Default::default()),
            Err(SmilesParseError::Atom { .. })
        ));
    }
}
