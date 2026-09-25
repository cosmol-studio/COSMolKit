use cosmolkit_bio::{
    AltLocLabel, AtomAddress, AtomName, AtomSourceIds, BioAsu, BioAtomId, BioAtomRow,
    BioBasicRefinementInfo, BioCalcFlag, BioChainId, BioChainRow, BioCisPep, BioConnection,
    BioConnectionKind, BioCoordinateBlock, BioCoordinateFormat, BioDiffractionInfo, BioEntityId,
    BioExperimentInfo, BioExperimentalCrystalInfo, BioHelix, BioHelixClass, BioMetadata, BioModRes,
    BioModelId, BioModelRow, BioRefinementInfo, BioRefinementRestraint, BioReflectionsInfo,
    BioResidueId, BioResidueRow, BioRowSpan, BioSheet, BioSiftsUnpResidue,
    BioSoftwareClassification, BioSoftwareItem, BioStrand, BioStructure, BioStructureError,
    BioStructureParts, BioStructureSourceState, BioTlsGroup, BioTlsSelection, BioTransform,
    ChainKind, ChainSourceIds, EntityKind, PdbAtomSerial, PdbChainId, PdbSeqId, ResidueAddress,
    ResidueInfoKind, ResidueName, ResidueSourceIds,
};
use cosmolkit_types::Element;
use std::collections::BTreeMap;

fn atom_row_with_isotope(name: [u8; 4], isotope_mass_number: Option<u16>) -> BioAtomRow {
    BioAtomRow::new(
        BioResidueId::new(7),
        AtomName::from_ascii(name).unwrap(),
        Element::H,
        isotope_mass_number,
        Some(AltLocLabel::new(b'B')),
        -2,
        BioCalcFlag::Calculated,
        0.625,
        17.25,
        [1.0, -2.0, 3.0, -4.0, 5.0, -6.0],
        12,
        0.375,
        AtomSourceIds::new(Some(PdbAtomSerial::new(481))),
    )
}

#[test]
fn bio_atom_row_distinguishes_unspecified_hydrogen_and_deuterium_state() {
    let ordinary_hydrogen = atom_row_with_isotope(*b" D1 ", None);
    let deuterium = atom_row_with_isotope(*b" D1 ", Some(2));

    assert_eq!(ordinary_hydrogen.element(), Element::H);
    assert_eq!(deuterium.element(), Element::H);
    assert_eq!(ordinary_hydrogen.isotope_mass_number(), None);
    assert_eq!(deuterium.isotope_mass_number(), Some(2));
    assert_ne!(ordinary_hydrogen, deuterium);
    assert_eq!(deuterium.clone(), deuterium);

    let independently_named_deuterium = atom_row_with_isotope(*b" CA ", Some(2));
    assert_eq!(independently_named_deuterium.element(), Element::H);
    assert_eq!(independently_named_deuterium.isotope_mass_number(), Some(2));
    assert_ne!(independently_named_deuterium.name(), deuterium.name());

    assert_eq!(deuterium.residue_id(), BioResidueId::new(7));
    assert_eq!(deuterium.altloc(), Some(AltLocLabel::new(b'B')));
    assert_eq!(deuterium.formal_charge(), -2);
    assert_eq!(deuterium.calc_flag(), BioCalcFlag::Calculated);
    assert_eq!(deuterium.occupancy(), 0.625);
    assert_eq!(deuterium.b_iso(), 17.25);
    assert_eq!(deuterium.anisou(), &[1.0, -2.0, 3.0, -4.0, 5.0, -6.0]);
    assert_eq!(deuterium.tls_group_id(), 12);
    assert_eq!(deuterium.fraction(), 0.375);
    assert_eq!(
        deuterium.source().serial().map(PdbAtomSerial::value),
        Some(481)
    );
}

fn address(
    sequence_number: Option<i32>,
    insertion_code: Option<u8>,
    segment: &[u8],
    name: &[u8],
) -> ResidueAddress {
    ResidueAddress::new(
        sequence_number,
        insertion_code,
        segment,
        ResidueName::from_ascii(name).unwrap(),
    )
    .unwrap()
}

fn atom_address(
    chain: &[u8],
    sequence_number: Option<i32>,
    insertion_code: Option<u8>,
    segment: &[u8],
    residue_name: &[u8],
    logical_atom_name: &str,
    altloc: Option<u8>,
) -> AtomAddress {
    AtomAddress::new(
        PdbChainId::from_ascii(chain).unwrap(),
        address(sequence_number, insertion_code, segment, residue_name),
        logical_atom_name,
        altloc,
    )
}

#[test]
fn residue_address_matches_gemmi_sequence_segment_and_name_rules() {
    let source = address(Some(42), Some(b'A'), b"SEG1", b"ALA");
    let insertion_case_variant = address(Some(42), Some(b'a'), b"SEG1", b"ALA");
    assert!(source.matches(&insertion_case_variant));
    assert_eq!(source, insertion_case_variant);

    assert!(!source.matches(&address(Some(43), Some(b'A'), b"SEG1", b"ALA")));
    assert!(!source.matches(&address(Some(42), Some(b'B'), b"SEG1", b"ALA")));

    let other_segment = address(Some(42), Some(b'A'), b"SEG2", b"ALA");
    assert!(!source.matches(&other_segment));
    assert!(source.matches_without_segment(&other_segment));

    let other_name = address(Some(42), Some(b'A'), b"SEG1", b"GLY");
    assert!(!source.matches(&other_name));
    assert!(!source.matches_without_segment(&other_name));

    // Gemmi's bit-mask comparison is not a general Unicode/ASCII lowercase
    // operation: it clears exactly bit 0x20 after XOR, including punctuation.
    assert!(
        address(Some(42), Some(b'@'), b"SEG1", b"ALA").matches(&address(
            Some(42),
            Some(b'`'),
            b"SEG1",
            b"ALA"
        ))
    );
}

#[test]
fn residue_address_canonicalizes_gemmi_absence_and_keeps_exact_text() {
    let absent = address(None, None, b"", b"ALA");
    let source_sentinels = address(Some(i32::MIN), Some(b' '), b"", b"ALA");
    assert!(absent.matches(&source_sentinels));
    assert_eq!(absent.sequence_number(), None);
    assert_eq!(source_sentinels.insertion_code(), None);
    assert_eq!(absent.segment(), "");
    assert_eq!(absent.name().as_str(), "ALA");

    let exact_text = address(Some(-7), None, b"A b!", b"MSE ");
    assert_eq!(exact_text.sequence_number(), Some(-7));
    assert_eq!(exact_text.segment(), "A b!");
    assert_eq!(exact_text.name().as_str(), "MSE ");

    assert!(
        ResidueAddress::new(
            Some(1),
            None,
            b"ABCDE",
            ResidueName::from_ascii(b"ALA").unwrap()
        )
        .is_none()
    );
    assert!(
        ResidueAddress::new(
            Some(1),
            None,
            &[0xc3, 0xa9],
            ResidueName::from_ascii(b"ALA").unwrap()
        )
        .is_none()
    );
}

#[test]
fn atom_address_matches_gemmi_chain_residue_name_and_altloc_fields() {
    let base = atom_address(b"A", Some(42), Some(b'A'), b"SEG1", b"ALA", "CA", None);
    assert_eq!(base.chain_name().as_str(), "A");
    assert_eq!(
        base.residue(),
        address(Some(42), Some(b'A'), b"SEG1", b"ALA")
    );
    assert_eq!(base.logical_atom_name(), "CA");
    assert_eq!(base.altloc(), 0);

    let insertion_case_variant =
        atom_address(b"A", Some(42), Some(b'a'), b"SEG1", b"ALA", "CA", Some(0));
    assert_eq!(base, insertion_case_variant);
    assert_eq!(
        base,
        atom_address(
            b"A",
            Some(42),
            Some(b'A'),
            b"SEG1",
            b"ALA",
            "CA",
            Some(b'\0')
        )
    );

    assert_ne!(
        base,
        atom_address(b"a", Some(42), Some(b'A'), b"SEG1", b"ALA", "CA", None)
    );
    assert_ne!(
        base,
        atom_address(b"A", Some(43), Some(b'A'), b"SEG1", b"ALA", "CA", None)
    );
    assert_ne!(
        base,
        atom_address(b"A", Some(42), Some(b'B'), b"SEG1", b"ALA", "CA", None)
    );
    assert_ne!(
        base,
        atom_address(b"A", Some(42), Some(b'A'), b"SEG2", b"ALA", "CA", None)
    );
    assert_ne!(
        base,
        atom_address(b"A", Some(42), Some(b'A'), b"SEG1", b"GLY", "CA", None)
    );
    assert_ne!(
        base,
        atom_address(
            b"A",
            Some(42),
            Some(b'A'),
            b"SEG1",
            b"ALA",
            "CA",
            Some(b' ')
        )
    );
    assert_ne!(
        atom_address(
            b"A",
            Some(42),
            Some(b'A'),
            b"SEG1",
            b"ALA",
            "CA",
            Some(b'A')
        ),
        atom_address(
            b"A",
            Some(42),
            Some(b'A'),
            b"SEG1",
            b"ALA",
            "CA",
            Some(b'a')
        )
    );
}

#[test]
fn atom_address_retains_logical_atom_name_without_four_column_padding() {
    let short = atom_address(b"A", Some(1), None, b"", b"ALA", "CA", None);
    let padded = atom_address(b"A", Some(1), None, b"", b"ALA", " CA ", None);
    let long = atom_address(
        b"A",
        Some(1),
        None,
        b"",
        b"ALA",
        "logical atom name beyond four bytes",
        None,
    );

    assert_ne!(short, padded);
    assert_eq!(short.logical_atom_name(), "CA");
    assert_eq!(padded.logical_atom_name(), " CA ");
    assert_eq!(
        long.logical_atom_name(),
        "logical atom name beyond four bytes"
    );
}

#[test]
fn bio_connection_defaults_match_gemmi_member_initializers() {
    assert_eq!(BioConnectionKind::Covale as u8, 0);
    assert_eq!(BioConnectionKind::Disulf as u8, 1);
    assert_eq!(BioConnectionKind::Hydrog as u8, 2);
    assert_eq!(BioConnectionKind::MetalC as u8, 3);
    assert_eq!(BioConnectionKind::Unknown as u8, 4);
    assert_eq!(BioConnectionKind::default(), BioConnectionKind::Unknown);

    assert_eq!(BioAsu::Same as u8, 0);
    assert_eq!(BioAsu::Different as u8, 1);
    assert_eq!(BioAsu::Any as u8, 2);
    assert_eq!(BioAsu::default(), BioAsu::Any);

    let connection = BioConnection::default();
    assert_eq!(connection.name, "");
    assert_eq!(connection.link_id, "");
    assert_eq!(connection.kind, BioConnectionKind::Unknown);
    assert_eq!(connection.asu, BioAsu::Any);
    assert_eq!(connection.reported_distance.to_bits(), 0.0_f64.to_bits());
    assert_eq!(connection.reported_sym, [0_i16; 4]);

    for partner in [&connection.partner1, &connection.partner2] {
        assert_eq!(partner.chain_name().as_str(), "");
        let residue = partner.residue();
        assert_eq!(residue.sequence_number(), None);
        assert_eq!(residue.insertion_code(), None);
        assert_eq!(residue.segment(), "");
        assert_eq!(residue.name().as_str(), "");
        assert_eq!(partner.logical_atom_name(), "");
        assert_eq!(partner.altloc(), 0);
    }
}

#[test]
fn bio_connection_retains_both_partners_and_reported_fields() {
    let partner1 = atom_address(
        b"A",
        Some(12),
        Some(b'A'),
        b"SEG1",
        b"CYS",
        "SG",
        Some(b'B'),
    );
    let partner2 = atom_address(b"B", Some(31), None, b"", b"CYS", "SG", None);
    let connection = BioConnection {
        name: "disulf1".to_owned(),
        link_id: "SSBOND".to_owned(),
        kind: BioConnectionKind::Disulf,
        asu: BioAsu::Different,
        partner1: partner1.clone(),
        partner2: partner2.clone(),
        reported_distance: 2.032,
        reported_sym: [1, -2, 3, -4],
    };

    assert_eq!(connection.name, "disulf1");
    assert_eq!(connection.link_id, "SSBOND");
    assert_eq!(connection.kind, BioConnectionKind::Disulf);
    assert_eq!(connection.asu, BioAsu::Different);
    assert_eq!(connection.partner1, partner1);
    assert_eq!(connection.partner2, partner2);
    assert_eq!(connection.reported_distance.to_bits(), 2.032_f64.to_bits());
    assert_eq!(connection.reported_sym, [1, -2, 3, -4]);
}

#[test]
fn bio_cis_pep_defaults_match_gemmi_member_initializers() {
    let cis_peptide = BioCisPep::default();

    assert_eq!(cis_peptide.partner_c, AtomAddress::default());
    assert_eq!(cis_peptide.partner_n, AtomAddress::default());
    assert_eq!(cis_peptide.model_num, 0);
    assert_eq!(cis_peptide.only_altloc, 0);
    assert!(cis_peptide.reported_angle.is_nan());
}

#[test]
fn bio_cis_pep_retains_two_addresses_model_altloc_and_reported_angle() {
    let partner_c = atom_address(b"A", Some(7), None, b"", b"PRO", "C", Some(b'A'));
    let partner_n = atom_address(b"A", Some(8), None, b"", b"PRO", "N", Some(b'A'));
    let cis_peptide = BioCisPep {
        partner_c: partner_c.clone(),
        partner_n: partner_n.clone(),
        model_num: 2,
        only_altloc: b'A',
        reported_angle: 178.25,
    };

    assert_eq!(cis_peptide.partner_c, partner_c);
    assert_eq!(cis_peptide.partner_n, partner_n);
    assert_eq!(cis_peptide.model_num, 2);
    assert_eq!(cis_peptide.only_altloc, b'A');
    assert_eq!(cis_peptide.reported_angle.to_bits(), 178.25_f64.to_bits());
}

#[test]
fn bio_mod_res_default_matches_gemmi_member_defaults() {
    let mod_res = BioModRes::default();

    assert_eq!(mod_res.chain_name.as_str(), "");
    assert_eq!(mod_res.res_id.sequence_number(), None);
    assert_eq!(mod_res.res_id.insertion_code(), None);
    assert_eq!(mod_res.res_id.segment(), "");
    assert_eq!(mod_res.res_id.name().as_str(), "");
    assert_eq!(mod_res.parent_comp_id, "");
    assert_eq!(mod_res.mod_id, "");
    assert_eq!(mod_res.details, "");
}

#[test]
fn bio_mod_res_retains_source_strings_and_approved_residue_name_boundary() {
    let chain_name = PdbChainId::from_ascii(b"ABCD").unwrap();
    let res_id = address(Some(55), Some(b'B'), b"SEG1", b"MOD4");
    let mod_res = BioModRes {
        chain_name,
        res_id,
        parent_comp_id: "PARENT-COMPONENT-LONGER-THAN-FOUR".to_owned(),
        mod_id: "Refmac extension identifier".to_owned(),
        details: "retained source details text".to_owned(),
    };

    assert_eq!(mod_res.chain_name.as_str(), "ABCD");
    assert_eq!(mod_res.res_id.sequence_number(), Some(55));
    assert_eq!(mod_res.res_id.insertion_code(), Some(b'B'));
    assert_eq!(mod_res.res_id.segment(), "SEG1");
    assert_eq!(mod_res.res_id.name().as_str(), "MOD4");
    assert_eq!(mod_res.parent_comp_id, "PARENT-COMPONENT-LONGER-THAN-FOUR");
    assert_eq!(mod_res.mod_id, "Refmac extension identifier");
    assert_eq!(mod_res.details, "retained source details text");

    assert!(PdbChainId::from_ascii(b"ABCDE").is_none());
    assert!(ResidueName::from_ascii(b"ABCDE").is_none());
    assert!(ResidueName::from_ascii(&[0xc3, 0xa9]).is_none());
}

#[test]
fn bio_helix_defaults_match_gemmi_member_initializers_and_unresolved_addresses() {
    let helix = BioHelix::default();

    assert_eq!(helix.start.chain_name().as_str(), "");
    assert_eq!(helix.start.residue().sequence_number(), None);
    assert_eq!(helix.start.residue().insertion_code(), None);
    assert_eq!(helix.start.residue().segment(), "");
    assert_eq!(helix.start.residue().name().as_str(), "");
    assert_eq!(helix.start.logical_atom_name(), "");
    assert_eq!(helix.start.altloc(), 0);

    assert_eq!(helix.end.chain_name().as_str(), "");
    assert_eq!(helix.end.residue().sequence_number(), None);
    assert_eq!(helix.end.residue().insertion_code(), None);
    assert_eq!(helix.end.residue().segment(), "");
    assert_eq!(helix.end.residue().name().as_str(), "");
    assert_eq!(helix.end.logical_atom_name(), "");
    assert_eq!(helix.end.altloc(), 0);

    assert_eq!(helix.pdb_helix_class, BioHelixClass::UnknownHelix);
    assert_eq!(helix.length, -1);
}

#[test]
fn bio_helix_integer_class_setter_matches_gemmi_range_and_ordinals() {
    let mut helix = BioHelix::default();
    let expected = [
        BioHelixClass::RAlpha,
        BioHelixClass::ROmega,
        BioHelixClass::RPi,
        BioHelixClass::RGamma,
        BioHelixClass::R310,
        BioHelixClass::LAlpha,
        BioHelixClass::LOmega,
        BioHelixClass::LGamma,
        BioHelixClass::Helix27,
        BioHelixClass::HelixPolyProlineNone,
    ];

    helix.set_helix_class_as_int(0);
    assert_eq!(helix.pdb_helix_class, BioHelixClass::UnknownHelix);

    for (source_code, expected_class) in (1..=10).zip(expected) {
        helix.set_helix_class_as_int(source_code);
        assert_eq!(helix.pdb_helix_class, expected_class);
        assert_eq!(helix.pdb_helix_class as i32, source_code);
    }

    for out_of_range in [i32::MIN, -1, 0, 11, i32::MAX] {
        helix.set_helix_class_as_int(out_of_range);
        assert_eq!(helix.pdb_helix_class, BioHelixClass::HelixPolyProlineNone);
    }
}

#[test]
fn bio_sheet_defaults_and_id_constructor_match_gemmi() {
    let empty = BioSheet::default();
    assert_eq!(empty.name, "");
    assert!(empty.strands.is_empty());

    let named = BioSheet::new("SHEET-A");
    assert_eq!(named.name, "SHEET-A");
    assert!(named.strands.is_empty());
}

#[test]
fn bio_strand_preserves_ordered_endpoint_and_hbond_addresses_without_coordinates() {
    let start = atom_address(b"A", Some(11), Some(b'A'), b"", b"GLY", "CA", Some(b'B'));
    let end = atom_address(b"B", Some(12), None, b"", b"GLY", "C", None);
    let hbond_atom2 = atom_address(b"B", Some(9), None, b"", b"GLY", "O", None);
    let hbond_atom1 = atom_address(b"A", Some(10), None, b"", b"GLY", "N", None);
    let first = BioStrand::new(
        start.clone(),
        end.clone(),
        hbond_atom2.clone(),
        hbond_atom1.clone(),
        -1,
        "strand-one".to_owned(),
    );

    // These valid source addresses intentionally refer to no structure or
    // coordinates; Gemmi's value type stores addresses without resolving them.
    let absent_start = atom_address(b"NONE", Some(900), None, b"", b"UNK", "Q1", None);
    let absent_end = atom_address(b"NONE", Some(901), None, b"", b"UNK", "Q2", None);
    let second = BioStrand::new(
        absent_start.clone(),
        absent_end.clone(),
        AtomAddress::default(),
        AtomAddress::default(),
        4,
        "strand-two".to_owned(),
    );

    let mut sheet = BioSheet::new("SHEET-A");
    sheet.strands.push(first);
    sheet.strands.push(second);

    assert_eq!(sheet.strands.len(), 2);
    assert_eq!(sheet.strands[0].name, "strand-one");
    assert_eq!(sheet.strands[0].start, start);
    assert_eq!(sheet.strands[0].end, end);
    assert_eq!(sheet.strands[0].hbond_atom2, hbond_atom2);
    assert_eq!(sheet.strands[0].hbond_atom1, hbond_atom1);
    assert_eq!(sheet.strands[0].sense, -1);

    assert_eq!(sheet.strands[1].name, "strand-two");
    assert_eq!(sheet.strands[1].start, absent_start);
    assert_eq!(sheet.strands[1].end, absent_end);
    assert_eq!(sheet.strands[1].hbond_atom2, AtomAddress::default());
    assert_eq!(sheet.strands[1].hbond_atom1, AtomAddress::default());
    assert_eq!(sheet.strands[1].sense, 4);
}

#[test]
fn bio_software_item_defaults_match_gemmi_member_initializers() {
    let item = BioSoftwareItem::default();

    assert_eq!(item.name, "");
    assert_eq!(item.version, "");
    assert_eq!(item.date, "");
    assert_eq!(item.description, "");
    assert_eq!(item.contact_author, "");
    assert_eq!(item.contact_author_email, "");
    assert_eq!(item.classification, BioSoftwareClassification::Unspecified);
}

#[test]
fn bio_software_classification_preserves_gemmi_order_and_ordinals() {
    let values = [
        BioSoftwareClassification::DataCollection,
        BioSoftwareClassification::DataExtraction,
        BioSoftwareClassification::DataProcessing,
        BioSoftwareClassification::DataReduction,
        BioSoftwareClassification::DataScaling,
        BioSoftwareClassification::ModelBuilding,
        BioSoftwareClassification::Phasing,
        BioSoftwareClassification::Refinement,
        BioSoftwareClassification::Unspecified,
    ];

    assert_eq!(
        values.map(|classification| classification as i32),
        [0, 1, 2, 3, 4, 5, 6, 7, 8]
    );
}

#[test]
fn bio_software_item_preserves_every_source_string_field() {
    let item = BioSoftwareItem {
        name: "tool-name".to_owned(),
        version: "v-1".to_owned(),
        date: "2026-09-24".to_owned(),
        description: "desc".to_owned(),
        contact_author: "author".to_owned(),
        contact_author_email: "author@example.invalid".to_owned(),
        classification: BioSoftwareClassification::Refinement,
    };

    assert_eq!(
        [
            item.name.as_str(),
            item.version.as_str(),
            item.date.as_str(),
            item.description.as_str(),
            item.contact_author.as_str(),
            item.contact_author_email.as_str(),
        ],
        [
            "tool-name",
            "v-1",
            "2026-09-24",
            "desc",
            "author",
            "author@example.invalid",
        ]
    );
    assert_eq!(item.classification as i32, 7);
}

#[test]
fn bio_reflections_info_defaults_all_source_missing_statistics_to_nan() {
    let info = BioReflectionsInfo::default();

    assert!(info.resolution_high.is_nan());
    assert!(info.resolution_low.is_nan());
    assert!(info.completeness.is_nan());
    assert!(info.redundancy.is_nan());
    assert!(info.r_merge.is_nan());
    assert!(info.r_sym.is_nan());
    assert!(info.mean_i_over_sigma.is_nan());
}

#[test]
fn bio_reflections_info_preserves_each_statistic_in_source_order() {
    let info = BioReflectionsInfo {
        resolution_high: 9.0,
        resolution_low: 1.0,
        completeness: 0.875,
        redundancy: 4.0,
        r_merge: 0.125,
        r_sym: 0.25,
        mean_i_over_sigma: 32.0,
    };

    assert_eq!(
        [
            info.resolution_high,
            info.resolution_low,
            info.completeness,
            info.redundancy,
            info.r_merge,
            info.r_sym,
            info.mean_i_over_sigma,
        ],
        [9.0, 1.0, 0.875, 4.0, 0.125, 0.25, 32.0]
    );
}

#[test]
fn bio_basic_refinement_info_defaults_match_gemmi_member_initializers() {
    let info = BioBasicRefinementInfo::default();

    assert!(info.resolution_high.is_nan());
    assert!(info.resolution_low.is_nan());
    assert!(info.completeness.is_nan());
    assert_eq!(info.reflection_count, -1);
    assert_eq!(info.work_set_count, -1);
    assert_eq!(info.rfree_set_count, -1);
    assert!(info.r_all.is_nan());
    assert!(info.r_work.is_nan());
    assert!(info.r_free.is_nan());
    assert!(info.cc_fo_fc_work.is_nan());
    assert!(info.cc_fo_fc_free.is_nan());
    assert!(info.fsc_work.is_nan());
    assert!(info.fsc_free.is_nan());
    assert!(info.cc_intensity_work.is_nan());
    assert!(info.cc_intensity_free.is_nan());
}

#[test]
fn bio_basic_refinement_info_preserves_total_and_per_bin_values() {
    let total = BioBasicRefinementInfo {
        resolution_high: 9.5,
        resolution_low: 0.5,
        completeness: 0.8,
        reflection_count: 17,
        work_set_count: 13,
        rfree_set_count: 4,
        r_all: 0.11,
        r_work: 0.12,
        r_free: 0.13,
        cc_fo_fc_work: 0.91,
        cc_fo_fc_free: 0.92,
        fsc_work: 0.93,
        fsc_free: 0.94,
        cc_intensity_work: 0.95,
        cc_intensity_free: 0.96,
    };
    let per_bin: Vec<BioBasicRefinementInfo> = vec![total];

    let expected = BioBasicRefinementInfo {
        resolution_high: 9.5,
        resolution_low: 0.5,
        completeness: 0.8,
        reflection_count: 17,
        work_set_count: 13,
        rfree_set_count: 4,
        r_all: 0.11,
        r_work: 0.12,
        r_free: 0.13,
        cc_fo_fc_work: 0.91,
        cc_fo_fc_free: 0.92,
        fsc_work: 0.93,
        fsc_free: 0.94,
        cc_intensity_work: 0.95,
        cc_intensity_free: 0.96,
    };
    assert_eq!(total, expected);
    assert_eq!(per_bin, vec![expected]);
}

#[test]
fn bio_refinement_restraint_defaults_and_name_constructor_match_gemmi() {
    let default = BioRefinementRestraint::default();

    assert!(default.name.is_empty());
    assert_eq!(default.count, -1);
    assert!(default.weight.is_nan());
    assert!(default.function.is_empty());
    assert!(default.dev_ideal.is_nan());

    let named = BioRefinementRestraint::new("bond length");
    assert_eq!(named.name, "bond length");
    assert_eq!(named.count, -1);
    assert!(named.weight.is_nan());
    assert!(named.function.is_empty());
    assert!(named.dev_ideal.is_nan());
}

#[test]
fn bio_refinement_restraint_preserves_every_populated_source_field() {
    let restraint = BioRefinementRestraint {
        name: "angle".to_owned(),
        count: 12,
        weight: 0.75,
        function: "Harmonic".to_owned(),
        dev_ideal: 1.25,
    };

    assert_eq!(
        restraint,
        BioRefinementRestraint {
            name: "angle".to_owned(),
            count: 12,
            weight: 0.75,
            function: "Harmonic".to_owned(),
            dev_ideal: 1.25,
        }
    );
}

#[test]
fn bio_refinement_info_defaults_match_gemmi_member_initializers() {
    let info = BioRefinementInfo::default();

    assert!(info.basic.resolution_high.is_nan());
    assert!(info.basic.resolution_low.is_nan());
    assert!(info.basic.completeness.is_nan());
    assert_eq!(info.basic.reflection_count, -1);
    assert_eq!(info.basic.work_set_count, -1);
    assert_eq!(info.basic.rfree_set_count, -1);
    assert!(info.basic.r_all.is_nan());
    assert!(info.basic.r_work.is_nan());
    assert!(info.basic.r_free.is_nan());
    assert!(info.basic.cc_fo_fc_work.is_nan());
    assert!(info.basic.cc_fo_fc_free.is_nan());
    assert!(info.basic.fsc_work.is_nan());
    assert!(info.basic.fsc_free.is_nan());
    assert!(info.basic.cc_intensity_work.is_nan());
    assert!(info.basic.cc_intensity_free.is_nan());

    assert!(info.id.is_empty());
    assert!(info.cross_validation_method.is_empty());
    assert!(info.rfree_selection_method.is_empty());
    assert_eq!(info.bin_count, -1);
    assert!(info.bins.is_empty());
    assert!(info.mean_b.is_nan());
    assert!(info.aniso_b.iter().all(|value| value.is_nan()));
    assert!(info.luzzati_error.is_nan());
    assert!(info.dpi_blow_r.is_nan());
    assert!(info.dpi_blow_rfree.is_nan());
    assert!(info.dpi_cruickshank_r.is_nan());
    assert!(info.dpi_cruickshank_rfree.is_nan());
    assert!(info.restr_stats.is_empty());
    assert!(info.tls_groups.is_empty());
    assert!(info.remarks.is_empty());
}

#[test]
fn bio_refinement_info_preserves_all_fields_and_collection_order() {
    let basic = BioBasicRefinementInfo {
        resolution_high: 2.1,
        resolution_low: 40.0,
        completeness: 0.98,
        reflection_count: 812,
        work_set_count: 770,
        rfree_set_count: 42,
        r_all: 0.16,
        r_work: 0.15,
        r_free: 0.19,
        cc_fo_fc_work: 0.95,
        cc_fo_fc_free: 0.89,
        fsc_work: 0.92,
        fsc_free: 0.87,
        cc_intensity_work: 0.91,
        cc_intensity_free: 0.86,
    };
    let first_bin = BioBasicRefinementInfo {
        resolution_high: 2.1,
        resolution_low: 3.0,
        completeness: 0.99,
        reflection_count: 311,
        work_set_count: 295,
        rfree_set_count: 16,
        r_all: 0.14,
        r_work: 0.13,
        r_free: 0.18,
        cc_fo_fc_work: 0.96,
        cc_fo_fc_free: 0.9,
        fsc_work: 0.93,
        fsc_free: 0.88,
        cc_intensity_work: 0.92,
        cc_intensity_free: 0.87,
    };
    let second_bin = BioBasicRefinementInfo {
        resolution_high: 3.0,
        resolution_low: 40.0,
        completeness: 0.95,
        reflection_count: 501,
        work_set_count: 475,
        rfree_set_count: 26,
        r_all: 0.18,
        r_work: 0.17,
        r_free: 0.21,
        cc_fo_fc_work: 0.93,
        cc_fo_fc_free: 0.85,
        fsc_work: 0.9,
        fsc_free: 0.84,
        cc_intensity_work: 0.89,
        cc_intensity_free: 0.83,
    };
    let first_restraint = BioRefinementRestraint {
        name: "bond length".to_owned(),
        count: 121,
        weight: 0.8,
        function: "Harmonic".to_owned(),
        dev_ideal: 0.02,
    };
    let second_restraint = BioRefinementRestraint {
        name: "angle".to_owned(),
        count: 87,
        weight: 1.25,
        function: "Cosine".to_owned(),
        dev_ideal: 1.5,
    };
    let tls = BioTlsGroup {
        num_id: 4,
        id: "TLS-4".to_owned(),
        selections: vec![BioTlsSelection {
            chain: PdbChainId::from_ascii(b"AB").unwrap(),
            res_begin: PdbSeqId::new(8, Some(b'A')),
            res_end: PdbSeqId::new(29, Some(b'B')),
            details: "source-ordered selection".to_owned(),
        }],
        origin: [1.25, -2.5, 3.75],
        t: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        l: [7.0, 8.0, 9.0, 10.0, 11.0, 12.0],
        s: [[13.0, 14.0, 15.0], [16.0, 17.0, 18.0], [19.0, 20.0, 21.0]],
    };
    let info = BioRefinementInfo {
        basic,
        id: "refinement-2".to_owned(),
        cross_validation_method: "cross-validation method".to_owned(),
        rfree_selection_method: "free-set method".to_owned(),
        bin_count: 2,
        bins: vec![first_bin, second_bin],
        mean_b: 24.5,
        aniso_b: [31.0, 32.0, 33.0, 34.0, 35.0, 36.0],
        luzzati_error: 0.17,
        dpi_blow_r: 0.021,
        dpi_blow_rfree: 0.032,
        dpi_cruickshank_r: 0.043,
        dpi_cruickshank_rfree: 0.054,
        restr_stats: vec![first_restraint.clone(), second_restraint.clone()],
        tls_groups: vec![tls.clone()],
        remarks: "REMARK 3 refinement details".to_owned(),
    };

    assert_eq!(info.basic, basic);
    assert_eq!(info.id, "refinement-2");
    assert_eq!(info.cross_validation_method, "cross-validation method");
    assert_eq!(info.rfree_selection_method, "free-set method");
    assert_eq!(info.bin_count, 2);
    assert_eq!(info.bins, [first_bin, second_bin]);
    assert_eq!(info.mean_b, 24.5);
    assert_eq!(info.aniso_b, [31.0, 32.0, 33.0, 34.0, 35.0, 36.0]);
    assert_eq!(info.luzzati_error, 0.17);
    assert_eq!(info.dpi_blow_r, 0.021);
    assert_eq!(info.dpi_blow_rfree, 0.032);
    assert_eq!(info.dpi_cruickshank_r, 0.043);
    assert_eq!(info.dpi_cruickshank_rfree, 0.054);
    assert_eq!(info.restr_stats, [first_restraint, second_restraint]);
    assert_eq!(info.tls_groups, [tls]);
    assert_eq!(info.remarks, "REMARK 3 refinement details");
}

#[test]
fn bio_metadata_defaults_match_gemmi_member_initializers() {
    let metadata = BioMetadata::default();

    assert!(metadata.authors.is_empty());
    assert!(metadata.experiments.is_empty());
    assert!(metadata.crystals.is_empty());
    assert!(metadata.refinement.is_empty());
    assert!(metadata.software.is_empty());
    assert!(metadata.solved_by.is_empty());
    assert!(metadata.starting_model.is_empty());
    assert!(metadata.remark_300_detail.is_empty());
}

#[test]
fn bio_metadata_preserves_all_source_ordered_members() {
    let metadata = BioMetadata {
        authors: vec!["Author B".to_owned(), "Author A".to_owned()],
        experiments: vec![
            BioExperimentInfo {
                method: "X-ray diffraction".to_owned(),
                number_of_crystals: 1,
                unique_reflections: 120,
                reflections: BioReflectionsInfo {
                    resolution_high: 2.1,
                    resolution_low: 38.0,
                    completeness: 0.97,
                    redundancy: 3.2,
                    r_merge: 0.08,
                    r_sym: 0.09,
                    mean_i_over_sigma: 15.0,
                },
                b_wilson: 19.5,
                shells: Vec::new(),
                diffraction_ids: vec!["dif-A".to_owned()],
            },
            BioExperimentInfo {
                method: "Electron microscopy".to_owned(),
                number_of_crystals: 2,
                unique_reflections: 240,
                reflections: BioReflectionsInfo {
                    resolution_high: 3.0,
                    resolution_low: 52.0,
                    completeness: 0.88,
                    redundancy: 2.4,
                    r_merge: 0.12,
                    r_sym: 0.13,
                    mean_i_over_sigma: 9.0,
                },
                b_wilson: 27.0,
                shells: Vec::new(),
                diffraction_ids: vec!["dif-B".to_owned()],
            },
        ],
        crystals: vec![
            BioExperimentalCrystalInfo {
                id: "crystal-B".to_owned(),
                description: "second source crystal".to_owned(),
                ph: 7.25,
                ph_range: "7.0-7.5".to_owned(),
                diffractions: Vec::new(),
            },
            BioExperimentalCrystalInfo {
                id: "crystal-A".to_owned(),
                description: "first source crystal".to_owned(),
                ph: 5.5,
                ph_range: "5.0-6.0".to_owned(),
                diffractions: Vec::new(),
            },
        ],
        refinement: vec![
            BioRefinementInfo {
                id: "refinement-B".to_owned(),
                mean_b: 25.0,
                ..BioRefinementInfo::default()
            },
            BioRefinementInfo {
                id: "refinement-A".to_owned(),
                mean_b: 17.5,
                ..BioRefinementInfo::default()
            },
        ],
        software: vec![
            BioSoftwareItem {
                name: "refiner B".to_owned(),
                version: "2.0".to_owned(),
                classification: BioSoftwareClassification::Refinement,
                ..BioSoftwareItem::default()
            },
            BioSoftwareItem {
                name: "extractor A".to_owned(),
                version: "1.0".to_owned(),
                classification: BioSoftwareClassification::DataExtraction,
                ..BioSoftwareItem::default()
            },
        ],
        solved_by: "molecular replacement".to_owned(),
        starting_model: "template-42".to_owned(),
        remark_300_detail: "biological assembly details".to_owned(),
    };

    assert_eq!(metadata.authors, ["Author B", "Author A"]);
    assert_eq!(
        metadata
            .experiments
            .iter()
            .map(|experiment| experiment.method.as_str())
            .collect::<Vec<_>>(),
        ["X-ray diffraction", "Electron microscopy"]
    );
    assert_eq!(metadata.experiments[0].reflections.r_merge, 0.08);
    assert_eq!(metadata.experiments[1].diffraction_ids, ["dif-B"]);
    assert_eq!(
        metadata
            .crystals
            .iter()
            .map(|crystal| crystal.id.as_str())
            .collect::<Vec<_>>(),
        ["crystal-B", "crystal-A"]
    );
    assert_eq!(metadata.crystals[0].ph, 7.25);
    assert_eq!(
        metadata
            .refinement
            .iter()
            .map(|refinement| refinement.id.as_str())
            .collect::<Vec<_>>(),
        ["refinement-B", "refinement-A"]
    );
    assert_eq!(metadata.refinement[1].mean_b, 17.5);
    assert_eq!(
        metadata
            .software
            .iter()
            .map(|software| software.name.as_str())
            .collect::<Vec<_>>(),
        ["refiner B", "extractor A"]
    );
    assert_eq!(
        metadata.software[0].classification,
        BioSoftwareClassification::Refinement
    );
    assert_eq!(metadata.solved_by, "molecular replacement");
    assert_eq!(metadata.starting_model, "template-42");
    assert_eq!(metadata.remark_300_detail, "biological assembly details");
}

#[test]
fn bio_structure_source_state_defaults_match_gemmi_structure_and_origx() {
    let state = BioStructureSourceState::default();

    assert!(state.name.is_empty());
    assert_eq!(state.resolution.to_bits(), 0.0_f64.to_bits());
    assert!(state.conect_map.is_empty());
    assert!(!state.has_d_fraction);
    assert_eq!(state.non_ascii_line, 0);
    assert_eq!(state.ter_status, 0);
    assert!(!state.has_origx);
    assert!(state.info.is_empty());
    assert!(state.raw_remarks.is_empty());
    assert_eq!(
        state.origx.matrix(),
        &[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    );
    assert_eq!(state.origx.translation(), &[0.0, 0.0, 0.0]);
}

#[test]
fn bio_structure_source_state_preserves_info_keys_remarks_and_conect_multiplicity() {
    let state = BioStructureSourceState {
        conect_map: BTreeMap::from([(42, vec![7, 7]), (7, vec![12, 12, 20])]),
        info: BTreeMap::from([
            ("_vendor.MixedCase ".to_owned(), "raw value ".to_owned()),
            ("_entry.id".to_owned(), "1abc".to_owned()),
            ("_cell.Z_PDB".to_owned(), "123".to_owned()),
        ]),
        raw_remarks: vec![
            "REMARK 2 RESOLUTION.    1.80 ANGSTROMS.".to_owned(),
            "REMARK 300  biological assembly details  ".to_owned(),
        ],
        ..BioStructureSourceState::default()
    };

    assert_eq!(
        state.info.keys().map(String::as_str).collect::<Vec<_>>(),
        ["_cell.Z_PDB", "_entry.id", "_vendor.MixedCase "]
    );
    assert_eq!(
        state.info.get("_entry.id").map(String::as_str),
        Some("1abc")
    );
    assert_eq!(
        state.info.get("_vendor.MixedCase ").map(String::as_str),
        Some("raw value ")
    );
    assert_eq!(
        state.raw_remarks,
        [
            "REMARK 2 RESOLUTION.    1.80 ANGSTROMS.",
            "REMARK 300  biological assembly details  ",
        ]
    );
    assert_eq!(
        state.conect_map,
        BTreeMap::from([(7, vec![12, 12, 20]), (42, vec![7, 7])])
    );
}

fn empty_bio_structure_parts() -> BioStructureParts {
    BioStructureParts {
        input_format: BioCoordinateFormat::Unknown,
        models: Vec::new(),
        chains: Vec::new(),
        residues: Vec::new(),
        atoms: Vec::new(),
        entities: Vec::new(),
        connections: Vec::new(),
        cispeps: Vec::new(),
        mod_residues: Vec::new(),
        helices: Vec::new(),
        sheets: Vec::new(),
        metadata: BioMetadata::default(),
        source_state: BioStructureSourceState::default(),
        coordinates: BioCoordinateBlock::default(),
        crystal: None,
        ncs_operators: Vec::new(),
        assemblies: Vec::new(),
    }
}

fn one_finite_coordinate_parts() -> BioStructureParts {
    let mut parts = empty_bio_structure_parts();
    parts.models.push(BioModelRow::new(
        BioRowSpan::<BioChainId>::new(0, 1).unwrap(),
        Some(1),
    ));
    parts.chains.push(BioChainRow::new(
        BioModelId::new(0),
        None,
        BioRowSpan::<BioResidueId>::new(0, 1).unwrap(),
        ChainKind::Protein,
        ChainSourceIds::new(
            Some(PdbChainId::from_ascii(b"A").unwrap()),
            Some("A".to_owned()),
        ),
    ));
    parts.residues.push(BioResidueRow::new(
        BioChainId::new(0),
        BioRowSpan::<BioAtomId>::new(0, 1).unwrap(),
        ResidueName::from_ascii(b"ALA").unwrap(),
        ResidueInfoKind::Aa,
        EntityKind::Polymer,
        None,
        None,
        ResidueSourceIds::new(None, None, None, None, None).unwrap(),
        BioSiftsUnpResidue::default(),
    ));
    parts.atoms.push(BioAtomRow::new(
        BioResidueId::new(0),
        AtomName::from_ascii(*b" CA ").unwrap(),
        Element::C,
        None,
        None,
        0,
        BioCalcFlag::NotSet,
        1.0,
        20.0,
        [0.0; 6],
        -1,
        0.0,
        AtomSourceIds::new(None),
    ));
    parts.coordinates = BioCoordinateBlock::new(vec![[1.25, -2.5, 4.0]]);
    parts
}

#[test]
fn bio_structure_round_trips_source_metadata_and_clone_is_independent() {
    let partner1 = atom_address(b"AA", Some(12), Some(b'A'), b"SEG1", b"CYS", "SG", None);
    let partner2 = atom_address(
        b"BB",
        Some(19),
        Some(b'B'),
        b"SEG2",
        b"CYS",
        "SG",
        Some(b'A'),
    );
    let connection = BioConnection {
        name: "disulfide-12-19".to_owned(),
        link_id: "SSBOND".to_owned(),
        kind: BioConnectionKind::Disulf,
        asu: BioAsu::Different,
        partner1: partner1.clone(),
        partner2: partner2.clone(),
        reported_distance: 2.03,
        reported_sym: [1, -2, 3, -4],
    };
    let cispep = BioCisPep {
        partner_c: partner1.clone(),
        partner_n: partner2.clone(),
        model_num: 2,
        only_altloc: b'B',
        reported_angle: f64::NAN,
    };
    let mod_residue = BioModRes {
        chain_name: PdbChainId::from_ascii(b"AA").unwrap(),
        res_id: address(Some(12), Some(b'A'), b"SEG1", b"MSE"),
        parent_comp_id: "MET".to_owned(),
        mod_id: "MSE".to_owned(),
        details: "selenomethionine".to_owned(),
    };
    let helix = BioHelix {
        start: partner1.clone(),
        end: partner2.clone(),
        pdb_helix_class: BioHelixClass::RAlpha,
        length: 8,
    };
    let sheet = BioSheet {
        name: "sheet-z".to_owned(),
        strands: vec![BioStrand::new(
            partner1.clone(),
            partner2.clone(),
            partner2.clone(),
            partner1.clone(),
            -1,
            "strand-2".to_owned(),
        )],
    };

    let experiment = BioExperimentInfo {
        method: "X-ray diffraction".to_owned(),
        number_of_crystals: 2,
        unique_reflections: 712,
        reflections: BioReflectionsInfo {
            resolution_high: 1.25,
            resolution_low: 37.0,
            completeness: 0.975,
            redundancy: 4.5,
            r_merge: 0.045,
            r_sym: 0.051,
            mean_i_over_sigma: 21.0,
        },
        b_wilson: 18.75,
        shells: vec![BioReflectionsInfo {
            resolution_high: 2.0,
            resolution_low: 37.0,
            completeness: 0.91,
            redundancy: 3.0,
            r_merge: 0.08,
            r_sym: 0.09,
            mean_i_over_sigma: 12.0,
        }],
        diffraction_ids: vec!["diffrn-B".to_owned(), "diffrn-A".to_owned()],
    };
    let mut crystal = BioExperimentalCrystalInfo::default();
    crystal.id = "crystal-3".to_owned();
    crystal.description = "orthorhombic crystal".to_owned();
    crystal.ph = 6.8;
    crystal.ph_range = "6.6-7.0".to_owned();
    let mut refinement = BioRefinementInfo::default();
    refinement.id = "refine-2".to_owned();
    refinement.remarks = "source refinement remark".to_owned();
    let software = BioSoftwareItem {
        name: "refiner".to_owned(),
        version: "2.4".to_owned(),
        date: "2026-09-24".to_owned(),
        description: "refinement program".to_owned(),
        contact_author: "A. Author".to_owned(),
        contact_author_email: "author@example.invalid".to_owned(),
        classification: BioSoftwareClassification::Refinement,
    };
    let metadata = BioMetadata {
        authors: vec!["Author Z".to_owned(), "Author A".to_owned()],
        experiments: vec![experiment.clone()],
        crystals: vec![crystal],
        refinement: vec![refinement],
        software: vec![software.clone()],
        solved_by: "molecular replacement".to_owned(),
        starting_model: "template-42".to_owned(),
        remark_300_detail: "assembly annotation".to_owned(),
    };
    let source_state = BioStructureSourceState {
        name: "entry-z".to_owned(),
        resolution: 1.84,
        conect_map: BTreeMap::from([(17, vec![9, 9, 12]), (9, vec![17])]),
        has_d_fraction: true,
        non_ascii_line: 203,
        ter_status: b'y',
        has_origx: true,
        origx: BioTransform::new(
            [[1.0, 0.25, 0.0], [0.0, 1.0, -0.5], [0.125, 0.0, 1.0]],
            [4.0, -3.0, 2.0],
        ),
        info: BTreeMap::from([
            ("_entry.id".to_owned(), "1xyz".to_owned()),
            ("_vendor.Flag ".to_owned(), " raw value ".to_owned()),
        ]),
        raw_remarks: vec!["REMARK 2 RESOLUTION. 1.84 ANGSTROMS.".to_owned()],
    };

    let mut parts = empty_bio_structure_parts();
    parts.connections.push(connection.clone());
    parts.cispeps.push(cispep);
    parts.mod_residues.push(mod_residue.clone());
    parts.helices.push(helix.clone());
    parts.sheets.push(sheet.clone());
    parts.metadata = metadata;
    parts.source_state = source_state.clone();

    let structure = BioStructure::from_parts(parts).unwrap();
    let returned = structure.clone().into_parts();
    assert_eq!(returned.connections, [connection]);
    assert_eq!(returned.cispeps.len(), 1);
    assert_eq!(returned.cispeps[0].partner_c, partner1);
    assert_eq!(returned.cispeps[0].partner_n, partner2);
    assert_eq!(returned.cispeps[0].model_num, 2);
    assert_eq!(returned.cispeps[0].only_altloc, b'B');
    assert!(returned.cispeps[0].reported_angle.is_nan());
    assert_eq!(returned.mod_residues, [mod_residue]);
    assert_eq!(returned.helices, [helix]);
    assert_eq!(returned.sheets, [sheet]);
    assert_eq!(returned.metadata.authors, ["Author Z", "Author A"]);
    assert_eq!(returned.metadata.experiments, [experiment]);
    assert_eq!(returned.metadata.crystals[0].id, "crystal-3");
    assert_eq!(
        returned.metadata.crystals[0].description,
        "orthorhombic crystal"
    );
    assert_eq!(returned.metadata.crystals[0].ph, 6.8);
    assert_eq!(returned.metadata.crystals[0].ph_range, "6.6-7.0");
    assert_eq!(returned.metadata.refinement[0].id, "refine-2");
    assert_eq!(
        returned.metadata.refinement[0].remarks,
        "source refinement remark"
    );
    assert!(returned.metadata.refinement[0].mean_b.is_nan());
    assert_eq!(returned.metadata.software, [software]);
    assert_eq!(returned.metadata.solved_by, "molecular replacement");
    assert_eq!(returned.metadata.starting_model, "template-42");
    assert_eq!(returned.metadata.remark_300_detail, "assembly annotation");
    assert_eq!(returned.source_state, source_state);

    let mut clone_parts = structure.clone().into_parts();
    clone_parts.connections[0].name.push_str("-clone");
    clone_parts.metadata.authors[0].push_str("-clone");
    clone_parts.source_state.raw_remarks[0].push_str("-clone");
    clone_parts.helices[0].length = 99;
    clone_parts.sheets[0].strands[0].sense = 1;
    assert_eq!(structure.connections()[0].name, "disulfide-12-19");
    assert_eq!(structure.metadata().authors[0], "Author Z");
    assert_eq!(
        structure.source_state().raw_remarks[0],
        "REMARK 2 RESOLUTION. 1.84 ANGSTROMS."
    );
    assert_eq!(structure.helices()[0].length, 8);
    assert_eq!(structure.sheets()[0].strands[0].sense, -1);
}

#[test]
fn bio_structure_metadata_defaults_and_nan_sentinels_survive_transport() {
    let empty = BioStructure::from_parts(empty_bio_structure_parts()).unwrap();
    assert!(empty.connections().is_empty());
    assert!(empty.cispeps().is_empty());
    assert!(empty.mod_residues().is_empty());
    assert!(empty.helices().is_empty());
    assert!(empty.sheets().is_empty());
    assert_eq!(empty.metadata(), &BioMetadata::default());
    assert_eq!(empty.source_state(), &BioStructureSourceState::default());

    let mut parts = empty_bio_structure_parts();
    parts.cispeps.push(BioCisPep::default());
    parts
        .metadata
        .experiments
        .push(BioExperimentInfo::default());
    parts
        .metadata
        .crystals
        .push(BioExperimentalCrystalInfo::default());
    parts.metadata.refinement.push(BioRefinementInfo::default());
    parts.metadata.software.push(BioSoftwareItem::default());
    let structure = BioStructure::from_parts(parts).unwrap();
    let returned = structure.into_parts();

    assert!(returned.cispeps[0].reported_angle.is_nan());
    assert!(returned.metadata.experiments[0].b_wilson.is_nan());
    assert!(
        returned.metadata.experiments[0]
            .reflections
            .resolution_high
            .is_nan()
    );
    assert!(returned.metadata.crystals[0].ph.is_nan());
    assert!(returned.metadata.refinement[0].mean_b.is_nan());
    assert!(
        returned.metadata.refinement[0]
            .aniso_b
            .iter()
            .all(|component| component.is_nan())
    );
    assert_eq!(
        returned.metadata.software[0].classification,
        BioSoftwareClassification::Unspecified
    );
}

#[test]
fn bio_structure_validation_keeps_unresolved_source_addresses_and_nan_sentinels() {
    let absent_partner_1 = atom_address(
        b"X1",
        Some(902),
        Some(b'Z'),
        b"NOPE",
        b"ZZZ",
        "Q1",
        Some(b'A'),
    );
    let absent_partner_2 = atom_address(b"X2", Some(-711), None, b"MISS", b"XXX", "Q2", None);
    let mut parts = one_finite_coordinate_parts();
    parts.connections.push(BioConnection {
        name: "unresolved-link".to_owned(),
        link_id: String::new(),
        kind: BioConnectionKind::Unknown,
        asu: BioAsu::Any,
        partner1: absent_partner_1.clone(),
        partner2: absent_partner_2.clone(),
        reported_distance: 2.5,
        reported_sym: [0; 4],
    });
    parts.cispeps.push(BioCisPep {
        partner_c: absent_partner_1.clone(),
        partner_n: absent_partner_2.clone(),
        model_num: 9,
        only_altloc: 0,
        reported_angle: f64::NAN,
    });
    parts.mod_residues.push(BioModRes {
        chain_name: PdbChainId::from_ascii(b"X3").unwrap(),
        res_id: address(Some(1188), None, b"", b"YYY"),
        parent_comp_id: "YYY".to_owned(),
        mod_id: "UNRESOLVED".to_owned(),
        details: String::new(),
    });
    parts.helices.push(BioHelix {
        start: absent_partner_1.clone(),
        end: absent_partner_2.clone(),
        pdb_helix_class: BioHelixClass::UnknownHelix,
        length: -1,
    });
    parts.sheets.push(BioSheet {
        name: "unresolved-sheet".to_owned(),
        strands: vec![BioStrand::new(
            absent_partner_1.clone(),
            absent_partner_2.clone(),
            absent_partner_2.clone(),
            absent_partner_1.clone(),
            0,
            String::new(),
        )],
    });
    parts.metadata.refinement.push(BioRefinementInfo::default());

    assert_eq!(BioStructure::validate_parts(&parts), Ok(()));
    let structure = BioStructure::from_parts(parts).unwrap();
    assert_eq!(structure.validate(), Ok(()));
    assert_eq!(structure.coordinates().positions(), &[[1.25, -2.5, 4.0]]);
    assert!(structure.cispeps()[0].reported_angle.is_nan());
    assert!(structure.metadata().refinement[0].mean_b.is_nan());
    assert_eq!(structure.connections()[0].partner1, absent_partner_1);
    assert_eq!(structure.connections()[0].partner2, absent_partner_2);
    assert_eq!(
        structure.mod_residues()[0].res_id.sequence_number(),
        Some(1188)
    );

    let mut invalid_span = structure.into_parts();
    invalid_span.models[0] =
        BioModelRow::new(BioRowSpan::<BioChainId>::new(1, 1).unwrap(), Some(1));
    assert_eq!(
        BioStructure::validate_parts(&invalid_span),
        Err(BioStructureError::NonContiguousSpan {
            table: "models->chains",
            expected_start: 0,
            actual_start: 1,
        })
    );
}

#[test]
fn bio_experiment_info_defaults_match_gemmi_member_initializers() {
    let info = BioExperimentInfo::default();

    assert!(info.method.is_empty());
    assert_eq!(info.number_of_crystals, -1);
    assert_eq!(info.unique_reflections, -1);
    assert!(info.reflections.resolution_high.is_nan());
    assert!(info.reflections.resolution_low.is_nan());
    assert!(info.reflections.completeness.is_nan());
    assert!(info.reflections.redundancy.is_nan());
    assert!(info.reflections.r_merge.is_nan());
    assert!(info.reflections.r_sym.is_nan());
    assert!(info.reflections.mean_i_over_sigma.is_nan());
    assert!(info.b_wilson.is_nan());
    assert!(info.shells.is_empty());
    assert!(info.diffraction_ids.is_empty());
}

#[test]
fn bio_experiment_info_preserves_source_ordered_fields_shells_and_diffraction_ids() {
    let info = BioExperimentInfo {
        method: "X-RAY DIFFRACTION".to_owned(),
        number_of_crystals: 3,
        unique_reflections: 401,
        reflections: BioReflectionsInfo {
            resolution_high: 9.0,
            resolution_low: 1.0,
            completeness: 0.875,
            redundancy: 4.0,
            r_merge: 0.125,
            r_sym: 0.25,
            mean_i_over_sigma: 32.0,
        },
        b_wilson: 18.5,
        shells: vec![
            BioReflectionsInfo {
                resolution_high: 2.0,
                resolution_low: 9.0,
                completeness: 0.75,
                redundancy: 3.0,
                r_merge: 0.1,
                r_sym: 0.2,
                mean_i_over_sigma: 12.0,
            },
            BioReflectionsInfo {
                resolution_high: 1.0,
                resolution_low: 3.0,
                completeness: 0.95,
                redundancy: 5.0,
                r_merge: 0.03,
                r_sym: 0.04,
                mean_i_over_sigma: 40.0,
            },
        ],
        diffraction_ids: vec!["diffrn-z".to_owned(), "diffrn-a".to_owned()],
    };

    assert_eq!(info.method, "X-RAY DIFFRACTION");
    assert_eq!(info.number_of_crystals, 3);
    assert_eq!(info.unique_reflections, 401);
    assert_eq!(
        [
            info.reflections.resolution_high,
            info.reflections.resolution_low,
            info.reflections.completeness,
            info.reflections.redundancy,
            info.reflections.r_merge,
            info.reflections.r_sym,
            info.reflections.mean_i_over_sigma,
            info.b_wilson,
        ],
        [9.0, 1.0, 0.875, 4.0, 0.125, 0.25, 32.0, 18.5]
    );
    assert_eq!(
        info.shells,
        [
            BioReflectionsInfo {
                resolution_high: 2.0,
                resolution_low: 9.0,
                completeness: 0.75,
                redundancy: 3.0,
                r_merge: 0.1,
                r_sym: 0.2,
                mean_i_over_sigma: 12.0,
            },
            BioReflectionsInfo {
                resolution_high: 1.0,
                resolution_low: 3.0,
                completeness: 0.95,
                redundancy: 5.0,
                r_merge: 0.03,
                r_sym: 0.04,
                mean_i_over_sigma: 40.0,
            },
        ]
    );
    assert_eq!(info.diffraction_ids, ["diffrn-z", "diffrn-a"]);
}

#[test]
fn bio_diffraction_info_defaults_match_gemmi_member_initializers() {
    let info = BioDiffractionInfo::default();

    assert!(info.id.is_empty());
    assert!(info.temperature.is_nan());
    assert!(info.source.is_empty());
    assert!(info.source_type.is_empty());
    assert!(info.synchrotron.is_empty());
    assert!(info.beamline.is_empty());
    assert!(info.wavelengths.is_empty());
    assert!(info.scattering_type.is_empty());
    assert_eq!(info.mono_or_laue, 0);
    assert!(info.monochromator.is_empty());
    assert!(info.collection_date.is_empty());
    assert!(info.optics.is_empty());
    assert!(info.detector.is_empty());
    assert!(info.detector_make.is_empty());
}

#[test]
fn bio_diffraction_info_preserves_source_ordered_fields_and_byte_code() {
    let info = BioDiffractionInfo {
        id: "diffrn-9".to_owned(),
        temperature: 298.15,
        source: "synchrotron source".to_owned(),
        source_type: "synchrotron".to_owned(),
        synchrotron: "facility".to_owned(),
        beamline: "BL-7".to_owned(),
        wavelengths: "0.9795,1.0000".to_owned(),
        scattering_type: "x-ray".to_owned(),
        mono_or_laue: b'M',
        monochromator: "Si(111)".to_owned(),
        collection_date: "2026-09-24".to_owned(),
        optics: "focused".to_owned(),
        detector: "pixel array".to_owned(),
        detector_make: "detector model".to_owned(),
    };

    assert_eq!(
        [
            info.id.as_str(),
            info.source.as_str(),
            info.source_type.as_str(),
            info.synchrotron.as_str(),
            info.beamline.as_str(),
            info.wavelengths.as_str(),
            info.scattering_type.as_str(),
            info.monochromator.as_str(),
            info.collection_date.as_str(),
            info.optics.as_str(),
            info.detector.as_str(),
            info.detector_make.as_str(),
        ],
        [
            "diffrn-9",
            "synchrotron source",
            "synchrotron",
            "facility",
            "BL-7",
            "0.9795,1.0000",
            "x-ray",
            "Si(111)",
            "2026-09-24",
            "focused",
            "pixel array",
            "detector model",
        ]
    );
    assert_eq!(info.temperature, 298.15);
    assert_eq!(info.mono_or_laue, b'M');
}

#[test]
fn bio_experimental_crystal_info_defaults_match_gemmi_member_initializers() {
    let info = BioExperimentalCrystalInfo::default();

    assert!(info.id.is_empty());
    assert!(info.description.is_empty());
    assert!(info.ph.is_nan());
    assert!(info.ph_range.is_empty());
    assert!(info.diffractions.is_empty());
}

#[test]
fn bio_experimental_crystal_info_preserves_fields_and_ordered_diffractions() {
    let diffraction = |id: &str, wavelengths: &str| {
        let mut info = BioDiffractionInfo::default();
        info.id = id.to_owned();
        info.wavelengths = wavelengths.to_owned();
        info
    };
    let info = BioExperimentalCrystalInfo {
        id: "xtal-9".to_owned(),
        description: "growth batch 9".to_owned(),
        ph: 5.7,
        ph_range: "4.8-6.1".to_owned(),
        diffractions: vec![
            diffraction("diffrn-9", "0.9795"),
            diffraction("diffrn-a", "1.0000"),
        ],
    };

    assert_eq!(info.id, "xtal-9");
    assert_eq!(info.description, "growth batch 9");
    assert_eq!(info.ph, 5.7);
    assert_eq!(info.ph_range, "4.8-6.1");
    assert_eq!(
        info.diffractions
            .iter()
            .map(|diffraction| (diffraction.id.as_str(), diffraction.wavelengths.as_str()))
            .collect::<Vec<_>>(),
        [("diffrn-9", "0.9795"), ("diffrn-a", "1.0000")]
    );
}

#[test]
fn bio_tls_selection_defaults_match_gemmi_member_initializers() {
    let selection = BioTlsSelection::default();

    assert_eq!(selection.chain.as_str(), "");
    assert_eq!(selection.res_begin.seq_num(), i32::MIN);
    assert_eq!(selection.res_begin.ins_code(), None);
    assert_eq!(selection.res_end.seq_num(), i32::MIN);
    assert_eq!(selection.res_end.ins_code(), None);
    assert_eq!(selection.details, "");
}

#[test]
fn bio_tls_selection_preserves_chain_range_insertions_and_details() {
    let selection = BioTlsSelection {
        chain: PdbChainId::from_ascii(b"ABCD").unwrap(),
        res_begin: PdbSeqId::new(-7, Some(b'B')),
        res_end: PdbSeqId::new(23, None),
        details: "selection details".to_owned(),
    };

    assert_eq!(selection.chain.as_str(), "ABCD");
    assert_eq!(selection.res_begin.seq_num(), -7);
    assert_eq!(selection.res_begin.ins_code(), Some(b'B'));
    assert_eq!(selection.res_end.seq_num(), 23);
    assert_eq!(selection.res_end.ins_code(), None);
    assert_eq!(selection.details, "selection details");
}

#[test]
fn bio_tls_group_defaults_match_gemmi_member_initializers() {
    let group = BioTlsGroup::default();

    assert_eq!(group.num_id, -1);
    assert_eq!(group.id, "");
    assert!(group.selections.is_empty());
    assert_eq!(group.origin, [0.0; 3]);
    assert!(group.t.iter().all(|component| component.is_nan()));
    assert!(group.l.iter().all(|component| component.is_nan()));
    assert!(
        group
            .s
            .iter()
            .flat_map(|row| row.iter())
            .all(|component| component.is_nan())
    );
}

#[test]
fn bio_tls_group_preserves_ids_selections_and_all_tensor_components() {
    let first_selection = BioTlsSelection {
        chain: PdbChainId::from_ascii(b"ABCD").unwrap(),
        res_begin: PdbSeqId::new(-7, Some(b'B')),
        res_end: PdbSeqId::new(23, None),
        details: "sel-one".to_owned(),
    };
    let second_selection = BioTlsSelection {
        chain: PdbChainId::from_ascii(b"Z").unwrap(),
        res_begin: PdbSeqId::new(24, None),
        res_end: PdbSeqId::new(31, Some(b'C')),
        details: "sel-two".to_owned(),
    };
    let group = BioTlsGroup {
        num_id: 7,
        id: "tls-7".to_owned(),
        selections: vec![first_selection.clone(), second_selection.clone()],
        origin: [1.25, -2.5, 3.75],
        t: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        l: [-1.0, -2.0, -3.0, -4.0, -5.0, -6.0],
        s: [[11.0, 12.0, 13.0], [21.0, 22.0, 23.0], [31.0, 32.0, 33.0]],
    };

    assert_eq!(group.num_id, 7);
    assert_eq!(group.id, "tls-7");
    assert_eq!(group.selections, [first_selection, second_selection]);
    assert_eq!(group.origin, [1.25, -2.5, 3.75]);
    assert_eq!(group.t, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    assert_eq!(group.l, [-1.0, -2.0, -3.0, -4.0, -5.0, -6.0]);
    assert_eq!(
        group.s,
        [[11.0, 12.0, 13.0], [21.0, 22.0, 23.0], [31.0, 32.0, 33.0]]
    );
}
