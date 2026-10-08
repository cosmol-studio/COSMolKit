use cosmolkit_core::KekulizeError;
use cosmolkit_model::{
    AtomId, BondId, Conformer2D, Conformer3D, CoordinateBlock, MoleculeProperties, PropertyValue,
    SGroupConnection, SGroupData, StereoGroup, StereoGroupKind, SubstanceGroup, SubstanceGroupId,
    SubstanceGroupKind, TopologyBlock, set_stereo_group_write_id,
};
use cosmolkit_smiles::{
    CxCoordinateSelection, CxSmilesFields, CxSmilesWriteParams, SmilesParseError, SmilesRecord,
    SmilesWriteParams, parse_smiles, parse_smiles_complete_source, write_cx_smiles_with_params,
    write_smiles_with_params,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag};

fn kekule_cx_params(fields: CxSmilesFields) -> CxSmilesWriteParams {
    cx_params(fields, true, true)
}

fn cx_params(
    fields: CxSmilesFields,
    do_isomeric_smiles: bool,
    do_kekule: bool,
) -> CxSmilesWriteParams {
    CxSmilesWriteParams {
        smiles: SmilesWriteParams {
            canonical: true,
            do_isomeric_smiles,
            clean_stereo: true,
            do_kekule,
            ..Default::default()
        },
        fields,
        ..Default::default()
    }
}

fn cx_coordinate_params(selection: CxCoordinateSelection) -> CxSmilesWriteParams {
    CxSmilesWriteParams {
        coordinate_selection: selection,
        ..cx_params(CxSmilesFields::COORDS, true, false)
    }
}

fn cx_atom_label_params() -> CxSmilesWriteParams {
    CxSmilesWriteParams {
        smiles: SmilesWriteParams {
            canonical: false,
            ..Default::default()
        },
        fields: CxSmilesFields::ATOM_LABELS,
        ..Default::default()
    }
}

fn cx_atom_value_params() -> CxSmilesWriteParams {
    CxSmilesWriteParams {
        smiles: SmilesWriteParams {
            canonical: false,
            ..Default::default()
        },
        fields: CxSmilesFields::MOLFILE_VALUES,
        ..Default::default()
    }
}

fn cx_radical_params() -> CxSmilesWriteParams {
    CxSmilesWriteParams {
        smiles: SmilesWriteParams {
            canonical: false,
            ..Default::default()
        },
        fields: CxSmilesFields::RADICALS,
        ..Default::default()
    }
}

fn cx_atom_prop_params() -> CxSmilesWriteParams {
    CxSmilesWriteParams {
        smiles: SmilesWriteParams {
            canonical: false,
            ..Default::default()
        },
        fields: CxSmilesFields::ATOM_PROPS,
        ..Default::default()
    }
}

fn cx_data_sgroup_params(canonical: bool) -> CxSmilesWriteParams {
    CxSmilesWriteParams {
        smiles: SmilesWriteParams {
            canonical,
            do_isomeric_smiles: true,
            clean_stereo: false,
            do_kekule: false,
            ..Default::default()
        },
        fields: CxSmilesFields::SGROUPS,
        ..Default::default()
    }
}

fn cx_polymer_params(canonical: bool) -> CxSmilesWriteParams {
    CxSmilesWriteParams {
        smiles: SmilesWriteParams {
            canonical,
            do_isomeric_smiles: true,
            clean_stereo: false,
            do_kekule: false,
            ..Default::default()
        },
        fields: CxSmilesFields::POLYMER,
        ..Default::default()
    }
}

// These original fixtures compare UTF-8 source strings. Decode only at the
// test observation boundary: invalid bytes fail instead of being substituted.
fn fixture_writer_text(text: cosmolkit_model::PropertyText) -> String {
    String::from_utf8(text.into_bytes()).expect("original writer fixture UTF-8 bytes")
}

#[test]
fn cx_writer_kekulizes_before_wedge_work_on_a_private_copy() {
    // Pinned RDKit 2026.03.1 direct MolToCXSmiles profile: sanitized
    // c1ccccc1, bond 0 set to BEGINWEDGE, canonical/isomeric/cleanStereo and
    // doKekule enabled, CX_NONE fields, default RestoreBondDirOptionClear.
    // The source output is C1=CC=CC=C1 and the input bond remains wedged.
    let mut record = parse_smiles("c1ccccc1", &Default::default()).expect("parse benzene");
    record.topology.bonds[0].set_direction(BondDirection::BeginWedge);
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &kekule_cx_params(CxSmilesFields::NONE))
            .map(fixture_writer_text)
            .unwrap(),
        "C1=CC=CC=C1"
    );
    assert_eq!(record, before, "CX writing must not mutate its caller");
}

#[test]
fn cx_writer_kekulize_error_preserves_typed_core_source_before_stereo_work() {
    // Pinned RDKit 2026.03.1 direct MolToCXSmiles first raises
    // AtomKekulizeException("non-ring atom 0 marked aromatic") for this input,
    // before the writer reaches its malformed _ringStereoAtoms string.
    let mut record = parse_smiles("C", &Default::default()).expect("parse carbon");
    record.topology.atoms[0].set_aromatic(true);
    record.topology.atoms[0]
        .set_prop("_ringStereoAtoms", "not-an-int-vector")
        .expect("set source malformed property");
    let before = record.clone();

    let error = write_cx_smiles_with_params(&record, &kekule_cx_params(CxSmilesFields::NONE))
        .map(fixture_writer_text)
        .expect_err("CX Kekulize runs before writer stereo/property processing");
    assert!(matches!(
        &error,
        SmilesParseError::WriterKekulize(KekulizeError::AromaticAtomOutsideRing { atom })
            if atom.index() == 0
    ));
    assert!(
        std::error::Error::source(&error)
            .and_then(|source| source.downcast_ref::<KekulizeError>())
            .is_some()
    );
    assert_eq!(record, before, "failed CX writing must preserve its caller");
}

#[test]
fn cx_writer_returns_empty_before_extensions_for_an_empty_molecule() {
    // Pinned RDKit 2026.03.1 MolToCXSmiles(Chem.Mol(), doKekule=true,
    // CX_ALL) returns the empty string before extension emission.
    let record = SmilesRecord {
        topology: TopologyBlock::default(),
        coordinates: CoordinateBlock::default(),
        properties: MoleculeProperties::default(),
    };
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &kekule_cx_params(CxSmilesFields::ALL))
            .map(fixture_writer_text)
            .unwrap(),
        ""
    );
    assert_eq!(record, before, "empty CX writing must preserve its caller");
}

#[test]
fn cx_writer_default_clear_preserves_wiggly_cfg_after_pre_kekulization() {
    // Pinned RDKit 2026.03.1 direct MolToCXSmiles profile: sanitized c1ccccc1,
    // bond 0 _MolFileBondCfg=2, canonical/isomeric/cleanStereo and doKekule
    // enabled, CX_ALL fields, default RestoreBondDirOptionClear. The wrapper
    // first kekulizes its private copy, then preserves the wiggly cfg for CX
    // emission; the exact output is C1=CC=CC=C1 |w:2.1|.
    let mut record = parse_smiles("c1ccccc1", &Default::default()).expect("parse benzene");
    record.topology.bonds[0]
        .set_prop("_MolFileBondCfg", "2")
        .expect("set source wiggly-bond config");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_params(CxSmilesFields::ALL, true, true))
            .map(fixture_writer_text)
            .unwrap(),
        "C1=CC=CC=C1 |w:2.1|"
    );
    assert_eq!(record, before, "CX preprocessing must preserve its caller");
}

#[test]
fn cx_writer_default_clear_removes_nonwiggly_bond_config_and_direction() {
    // Pinned RDKit 2026.03.1 direct MolToCXSmiles profile: sanitized OC(C)Cl,
    // bond 1 _MolFileBondCfg=1 plus BEGINWEDGE, canonical/isomeric/cleanStereo,
    // doKekule disabled, CX_ALL fields and default Clear. Non-wiggly cfg and
    // direction are both removed before extension emission; output is CC(O)Cl.
    let mut record = parse_smiles("OC(C)Cl", &Default::default()).expect("parse chlorohydrin");
    record.topology.bonds[1]
        .set_prop("_MolFileBondCfg", "1")
        .expect("set source non-wiggly bond config");
    record.topology.bonds[1].set_direction(BondDirection::BeginWedge);
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_params(CxSmilesFields::ALL, true, false))
            .map(fixture_writer_text)
            .unwrap(),
        "CC(O)Cl"
    );
    assert_eq!(record, before, "CX preprocessing must preserve its caller");
}

#[test]
fn cx_writer_nonisomeric_profile_masks_wiggly_bond_config() {
    // Pinned RDKit 2026.03.1 direct MolToCXSmiles profile: sanitized OC(C)Cl,
    // bond 1 _MolFileBondCfg=2, canonical/nonisomeric/cleanStereo, doKekule
    // disabled, CX_ALL fields and default Clear. The wrapper retains cfg=2
    // during preparation, then removes CX_BOND_CFG for nonisomeric output;
    // the exact output is CC(O)Cl.
    let mut record = parse_smiles("OC(C)Cl", &Default::default()).expect("parse chlorohydrin");
    record.topology.bonds[1]
        .set_prop("_MolFileBondCfg", "2")
        .expect("set source wiggly-bond config");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_params(CxSmilesFields::ALL, false, false))
            .map(fixture_writer_text)
            .unwrap(),
        "CC(O)Cl".to_owned()
    );
    assert_eq!(record, before, "CX preprocessing must preserve its caller");
}

#[test]
fn cx_writer_bond_config_serializes_prepared_wedges_wiggly_indices_and_field_gate() {
    // Pinned RDKit 2026.03.1 direct MolToCXSmiles profile: sanitized records,
    // canonical/isomeric/cleanStereo, doKekule=false, and the named CX fields.
    // Coordinate wedging produces the exact wU/wD blocks below. The OC wiggly
    // fixture proves that both the begin-atom index and bond index use final
    // output order, while disabling CX_BOND_CFG omits the block entirely.
    for (input, expected) in [
        (
            "CC(O)Cl |(-3.9163,5.4767,;-3.9163,3.9367,;-2.5826,3.1667,;-5.25,3.1667,),wU:1.0|",
            "C[C@@H](O)Cl |(-3.9163,5.4767,;-3.9163,3.9367,;-2.5826,3.1667,;-5.25,3.1667,),wU:1.0|",
        ),
        (
            "CC(O)Cl |(-3.9163,5.4767,;-3.9163,3.9367,;-2.5826,3.1667,;-5.25,3.1667,),wD:1.0|",
            "C[C@H](O)Cl |(-3.9163,5.4767,;-3.9163,3.9367,;-2.5826,3.1667,;-5.25,3.1667,),wD:1.0|",
        ),
    ] {
        let record = parse_smiles(input, &Default::default()).expect("parse prepared wedge row");
        let before = record.clone();
        assert_eq!(
            write_cx_smiles_with_params(
                &record,
                &cx_params(
                    CxSmilesFields::COORDS | CxSmilesFields::BOND_CFG,
                    true,
                    false
                ),
            )
            .map(fixture_writer_text)
            .unwrap(),
            expected.to_owned()
        );
        assert_eq!(record, before, "prepared wedge writing must be immutable");
    }

    let wiggly = parse_smiles("OC |w:0.0|", &Default::default()).expect("parse wiggly row");
    let before = wiggly.clone();
    assert_eq!(
        write_cx_smiles_with_params(&wiggly, &cx_params(CxSmilesFields::BOND_CFG, true, false),)
            .map(fixture_writer_text)
            .unwrap(),
        "CO |w:1.0|".to_owned()
    );
    assert_eq!(
        write_cx_smiles_with_params(&wiggly, &cx_params(CxSmilesFields::NONE, true, false))
            .map(fixture_writer_text)
            .unwrap(),
        "CO".to_owned()
    );
    assert_eq!(wiggly, before, "bond-config field gating must be immutable");
}

#[test]
fn cx_writer_coordinate_and_hydrogen_bonds_use_directed_begin_and_exact_type() {
    // Pinned RDKit 2026.03.1 MolToCXSmiles on sanitized NOC after replacing
    // source bond 1 with the named exact BondType. Canonical/isomeric output,
    // doKekule=false and CX_COORDINATE_BONDS|CX_HYDROGEN_BONDS produce these
    // exact rows. DATIVEONE is distinct from DATIVE and emits no C: section.
    let fields = CxSmilesFields::COORDINATE_BONDS | CxSmilesFields::HYDROGEN_BONDS;
    for (order, expected) in [
        (BondOrder::Dative, "C[OH]N |C:1.0|"),
        (BondOrder::Hydrogen, "CON |H:1.0|"),
        (BondOrder::DativeOne, "C~ON"),
    ] {
        let mut record = parse_smiles("NOC", &Default::default()).expect("parse NOC fixture");
        record.topology.bonds[1].set_order(order);
        let before = record.clone();

        assert_eq!(
            write_cx_smiles_with_params(&record, &cx_params(fields, true, false))
                .map(fixture_writer_text)
                .unwrap(),
            expected,
            "source bond type {order:?}"
        );
        assert_eq!(record, before, "directed bond writing must be immutable");
    }
}

#[test]
fn cx_writer_zero_bonds_use_final_bond_positions_in_source_order() {
    // Pinned RDKit 2026.03.1 MolToCXSmiles on sanitized NOCF after replacing
    // source bonds 0 and 2 with exact ZERO bonds. Canonical/isomeric output,
    // doKekule=false and CX_ZERO_BONDS emit the two final bond positions in
    // traversal order, retaining the nonzero bond between them.
    let mut record = parse_smiles("NOCF", &Default::default()).expect("parse NOCF fixture");
    record.topology.bonds[0].set_order(BondOrder::Zero);
    record.topology.bonds[2].set_order(BondOrder::Zero);
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_params(CxSmilesFields::ZERO_BONDS, true, false),)
            .map(fixture_writer_text)
            .unwrap(),
        "N~OC~F |Z:0,2|".to_owned()
    );
    assert_eq!(record, before, "zero-bond writing must be immutable");
}

fn ring_stereo_params(canonical: bool) -> CxSmilesWriteParams {
    CxSmilesWriteParams {
        smiles: SmilesWriteParams {
            canonical,
            do_isomeric_smiles: true,
            clean_stereo: false,
            do_kekule: false,
            ..Default::default()
        },
        fields: CxSmilesFields::BOND_CFG,
        ..Default::default()
    }
}

#[test]
fn cx_writer_ring_stereo_filters_nonring_small_ring_and_non_double_bonds() {
    // Pinned RDKit 2026.03.1 direct MolToCXSmiles with canonical=false,
    // isomeric=true, cleanStereo=false, doKekule=false and CX_BOND_CFG.
    // STEREOANY is ignored for a non-ring double bond, a seven-membered ring
    // double bond, and a ten-membered ring single bond, respectively.
    for (input, bond_index, expected) in [
        ("CC=C", 1, "CC=C"),
        ("C1=CCCCCC1", 0, "C1=CCCCCC1"),
        ("C1CCCCCCCCC1", 0, "C1CCCCCCCCC1"),
    ] {
        let mut record = parse_smiles(input, &Default::default())
            .unwrap_or_else(|error| panic!("parse ring eligibility row {input:?}: {error}"));
        record.topology.bonds[bond_index]
            .set_stereo(BondStereo::Any)
            .expect("set source STEREOANY state");
        let before = record.clone();

        assert_eq!(
            write_cx_smiles_with_params(&record, &ring_stereo_params(false))
                .map(fixture_writer_text)
                .unwrap(),
            expected,
            "ring eligibility row {input:?}"
        );
        assert_eq!(record, before, "eligibility writing must be immutable");
    }
}

#[test]
fn cx_writer_ring_stereo_concatenates_cis_trans_and_unknown_source_blocks() {
    // Pinned RDKit 2026.03.1 under the same noncanonical profile. The three
    // ten-membered rings carry CIS, TRANS and ANY in source component order.
    // get_ringbond_cistrans_block returns c + t + ctu without separators
    // between categories, while each category retains final bond positions.
    let mut record = parse_smiles(
        "C1=CCCCCCCCC1.C1=CCCCCCCCC1.C1=CCCCCCCCC1",
        &Default::default(),
    )
    .expect("parse three eligible source rings");
    record.topology.bonds[0].set_stereo_atoms(Some([AtomId::new(9), AtomId::new(2)]));
    record.topology.bonds[0]
        .set_stereo(BondStereo::Cis)
        .expect("set source CIS state");
    record.topology.bonds[9].set_stereo_atoms(Some([AtomId::new(19), AtomId::new(12)]));
    record.topology.bonds[9]
        .set_stereo(BondStereo::Trans)
        .expect("set source TRANS state");
    record.topology.bonds[18]
        .set_stereo(BondStereo::Any)
        .expect("set source ANY state");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &ring_stereo_params(false))
            .map(fixture_writer_text)
            .unwrap(),
        "C1=C\\CCCCCCCC/1.C1=C/CCCCCCCC/1.C1=CCCCCCCCC1 |c:0t:10ctu:20|"
    );
    assert_eq!(record, before, "ring stereo writing must be immutable");
}

#[test]
fn cx_writer_ring_stereo_preserves_boolean_atom_order_subscripts() {
    // Pinned RDKit 2026.03.1 canonical output for the branched ten-membered
    // ring. Both stereochemical reference pairs enter the degree>2 neighbor
    // loops. The source expression atomOrder[nbrIdx < o1] performs a boolean
    // subscript; replacing it with an output-position comparison emits c:1
    // for these TRANS rows instead of the pinned t:1.
    for (references, expected) in [
        (
            [AtomId::new(1), AtomId::new(4)],
            "F/C1=C(\\Cl)CCCCCCCC1 |t:1|",
        ),
        (
            [AtomId::new(11), AtomId::new(3)],
            "F/C1=C(\\Cl)CCCCCCCC1 |t:1|",
        ),
    ] {
        let mut record = parse_smiles("C1(Cl)=C(F)CCCCCCCC1", &Default::default())
            .expect("parse branched eligible ring");
        record.topology.bonds[1].set_stereo_atoms(Some(references));
        record.topology.bonds[1]
            .set_stereo(BondStereo::Trans)
            .expect("set source TRANS state");
        let before = record.clone();

        assert_eq!(
            write_cx_smiles_with_params(&record, &ring_stereo_params(true))
                .map(fixture_writer_text)
                .unwrap(),
            expected,
            "stereo references {references:?}"
        );
        assert_eq!(record, before, "neighbor-order writing must be immutable");
    }

    let mut cis =
        parse_smiles("C1(Cl)=C(F)CCCCCCCC1", &Default::default()).expect("parse branched CIS ring");
    cis.topology.bonds[1].set_stereo_atoms(Some([AtomId::new(1), AtomId::new(4)]));
    cis.topology.bonds[1]
        .set_stereo(BondStereo::Cis)
        .expect("set source CIS state");
    assert_eq!(
        write_cx_smiles_with_params(&cis, &ring_stereo_params(true))
            .map(fixture_writer_text)
            .unwrap(),
        "F/C1=C(/Cl)CCCCCCCC1 |c:1|"
    );
}

#[test]
fn cx_writer_clean_stereo_guard_and_group_cleanup_match_legacy_oracle() {
    // Pinned RDKit 2026.03.1 legacy profile. The two components have source
    // tetrahedral centers at atom rows 1 and 6; atom 6 is then made
    // unspecified while both remain in AND read group 7. With cleanStereo
    // false, the outer cleanup is skipped and output is `&1:1,6`. With
    // cleanStereo true, cleanup keeps only the still-specified atom and emits
    // `&1:6`. A present `_StereochemDone=0` skips assignment by property
    // presence but still reaches outer cleanup; absent/present marker results
    // therefore match for each cleanStereo setting.
    let parsed = parse_smiles("F[C@](Cl)(Br)I.F[C@](Cl)(Br)I", &Default::default())
        .expect("parse two source stereocenters");
    for (clean_stereo, done_marker, expected) in [
        (false, false, "FC(Cl)(Br)I.F[C@](Cl)(Br)I |&1:1,6|"),
        (false, true, "FC(Cl)(Br)I.F[C@](Cl)(Br)I |&1:1,6|"),
        (true, false, "FC(Cl)(Br)I.F[C@](Cl)(Br)I |&1:6|"),
        (true, true, "FC(Cl)(Br)I.F[C@](Cl)(Br)I |&1:6|"),
    ] {
        let mut record = parsed.clone();
        record.topology.atoms[6].set_chiral_tag(ChiralTag::Unspecified);
        record.topology.stereo_groups = vec![
            StereoGroup::new(
                StereoGroupKind::And,
                vec![AtomId::new(1), AtomId::new(6)],
                Vec::new(),
            )
            .with_id(7),
        ];
        if done_marker {
            record
                .properties
                .set_prop("_StereochemDone", "0")
                .expect("set present zero-like marker");
        }
        let before = record.clone();
        let params = CxSmilesWriteParams {
            smiles: SmilesWriteParams {
                canonical: true,
                do_isomeric_smiles: true,
                clean_stereo,
                ..Default::default()
            },
            fields: CxSmilesFields::ENHANCED_STEREO,
            ..Default::default()
        };

        assert_eq!(
            write_cx_smiles_with_params(&record, &params)
                .map(fixture_writer_text)
                .unwrap(),
            expected,
            "cleanStereo={clean_stereo}, done marker present={done_marker}"
        );
        assert_eq!(record, before, "CX writing must preserve its caller");
    }
}

#[test]
fn cx_writer_stereo_group_ids_keep_type_namespaces_holes_and_tie_order() {
    // Pinned RDKit 2026.03.1 direct MolToCXSmiles with canonical=false,
    // isomeric=true, cleanStereo=false, doKekule=false and
    // CX_ENHANCEDSTEREO. Sorting occurs by type, existing write ID and mapped
    // atom indices before assignment. OR ID 5 is duplicated, so its lower
    // mapped group keeps 5 and the later group receives 4. Missing OR and AND
    // IDs independently skip their respective reserved ID 1 and receive 2.
    let component = "F[C@](Cl)(Br)I";
    let input = std::iter::repeat_n(component, 8)
        .collect::<Vec<_>>()
        .join(".");
    let mut record = parse_smiles(&input, &Default::default()).expect("parse eight stereocenters");
    record.topology.stereo_groups = [
        (StereoGroupKind::Absolute, 31, 77, 1),
        (StereoGroupKind::Or, 7, 5, 6),
        (StereoGroupKind::Or, 11, 5, 11),
        (StereoGroupKind::Or, 13, 1, 16),
        (StereoGroupKind::Or, 17, 3, 21),
        (StereoGroupKind::Or, 19, 0, 26),
        (StereoGroupKind::And, 23, 1, 31),
        (StereoGroupKind::And, 29, 0, 36),
    ]
    .into_iter()
    .map(|(kind, read_id, write_id, atom)| {
        let mut group =
            StereoGroup::new(kind, vec![AtomId::new(atom)], Vec::new()).with_id(read_id);
        set_stereo_group_write_id(&mut group, write_id);
        group
    })
    .collect();
    let before = record.clone();
    let params = CxSmilesWriteParams {
        smiles: SmilesWriteParams {
            canonical: false,
            do_isomeric_smiles: true,
            clean_stereo: false,
            do_kekule: false,
            ..Default::default()
        },
        fields: CxSmilesFields::ENHANCED_STEREO,
        ..Default::default()
    };

    assert_eq!(
        write_cx_smiles_with_params(&record, &params)
            .map(fixture_writer_text)
            .unwrap(),
        concat!(
            "F[C@](Cl)(Br)I.F[C@](Cl)(Br)I.F[C@](Cl)(Br)I.F[C@](Cl)(Br)I.",
            "F[C@](Cl)(Br)I.F[C@](Cl)(Br)I.F[C@](Cl)(Br)I.F[C@](Cl)(Br)I ",
            "|a:1,o2:26,o1:16,o3:21,o5:6,o4:11,&2:36,&1:31|"
        )
        .to_owned()
    );
    assert_eq!(record, before, "group ID assignment must be immutable");
}

#[test]
fn cx_writer_enhanced_stereo_maps_atom_and_atrop_bond_members_in_group_order() {
    // Pinned RDKit 2026.03.1 direct MolToCXSmiles on the same biphenyl state:
    // canonical/isomeric, cleanStereo=false, doKekule=false, and
    // CX_BOND_ATROPISOMER|CX_ENHANCEDSTEREO. The input group order is
    // AND, OR, ABS; source sorting emits ABS, OR, AND. The bond-only OR group
    // receives mapped axis endpoint 3 from the shared atrop wedge map, while
    // direct atoms 0 and 11 map to output positions 10 and 9.
    let mut record = parse_smiles_complete_source("c1ccccc1-c1ccccc1", &Default::default())
        .expect("parse biphenyl enhanced-stereo fixture");
    let axis_id = record
        .topology
        .bonds
        .iter()
        .position(|bond| {
            (bond.begin() == AtomId::new(5) && bond.end() == AtomId::new(6))
                || (bond.begin() == AtomId::new(6) && bond.end() == AtomId::new(5))
        })
        .map(BondId::new)
        .expect("find central biphenyl bond");
    let axis = &mut record.topology.bonds[axis_id.index()];
    axis.set_stereo_atoms(Some([AtomId::new(4), AtomId::new(7)]));
    axis.set_stereo(BondStereo::AtropCw)
        .expect("set source atrop stereo");
    record.topology.stereo_groups = vec![
        StereoGroup::new(StereoGroupKind::And, vec![AtomId::new(11)], Vec::new()).with_id(29),
        StereoGroup::new(StereoGroupKind::Or, Vec::new(), vec![axis_id]).with_id(19),
        StereoGroup::new(StereoGroupKind::Absolute, vec![AtomId::new(0)], Vec::new()).with_id(31),
    ];
    let before = record.clone();
    let params = CxSmilesWriteParams {
        smiles: SmilesWriteParams {
            canonical: true,
            do_isomeric_smiles: true,
            clean_stereo: false,
            do_kekule: false,
            ..Default::default()
        },
        fields: CxSmilesFields::BOND_ATROPISOMER | CxSmilesFields::ENHANCED_STEREO,
        ..Default::default()
    };

    assert_eq!(
        write_cx_smiles_with_params(&record, &params)
            .map(fixture_writer_text)
            .unwrap(),
        "c1ccc(-c2ccccc2)cc1 |wU:3.2,a:10,o1:3,&1:9|"
    );
    assert_eq!(
        record, before,
        "enhanced stereo rendering must be immutable"
    );
}

#[test]
fn cx_writer_enhanced_stereo_preserves_overlapping_groups_with_atrop_wedge_map() {
    // Pinned RDKit 2026.03.1 direct MolToCXSmiles on biphenyl with the
    // canonical/isomeric, cleanStereo=false, doKekule=false profile and
    // CX_BOND_ATROPISOMER|CX_ENHANCEDSTEREO. AND/OR/ABS atom groups overlap
    // at original atoms 0, 6, and 11; a separate AND group covers atom 3.
    // The bond-only OR group receives mapped endpoint 4 through the same
    // Atropisomer wedge map used by the bond extension.
    let mut record = parse_smiles_complete_source("c1ccccc1-c1ccccc1", &Default::default())
        .expect("parse biphenyl overlapping enhanced-stereo fixture");
    let axis_id = record
        .topology
        .bonds
        .iter()
        .position(|bond| {
            (bond.begin() == AtomId::new(5) && bond.end() == AtomId::new(6))
                || (bond.begin() == AtomId::new(6) && bond.end() == AtomId::new(5))
        })
        .map(BondId::new)
        .expect("find central biphenyl bond");
    let axis = &mut record.topology.bonds[axis_id.index()];
    axis.set_stereo_atoms(Some([AtomId::new(4), AtomId::new(7)]));
    axis.set_stereo(BondStereo::AtropCw)
        .expect("set source atrop stereo");
    record.topology.stereo_groups = [
        (
            StereoGroupKind::And,
            vec![AtomId::new(0), AtomId::new(6)],
            Vec::new(),
            3,
        ),
        (
            StereoGroupKind::Or,
            vec![AtomId::new(6), AtomId::new(11)],
            Vec::new(),
            5,
        ),
        (StereoGroupKind::Or, Vec::new(), vec![axis_id], 7),
        (
            StereoGroupKind::Absolute,
            vec![AtomId::new(0), AtomId::new(11)],
            Vec::new(),
            9,
        ),
        (StereoGroupKind::And, vec![AtomId::new(3)], Vec::new(), 11),
    ]
    .into_iter()
    .map(|(kind, atoms, bonds, read_id)| StereoGroup::new(kind, atoms, bonds).with_id(read_id))
    .collect();
    let before = record.clone();
    let params = CxSmilesWriteParams {
        smiles: SmilesWriteParams {
            canonical: true,
            do_isomeric_smiles: true,
            clean_stereo: false,
            do_kekule: false,
            ..Default::default()
        },
        fields: CxSmilesFields::BOND_ATROPISOMER | CxSmilesFields::ENHANCED_STEREO,
        ..Default::default()
    };

    assert_eq!(
        write_cx_smiles_with_params(&record, &params)
            .map(fixture_writer_text)
            .unwrap(),
        "c1ccc(-c2ccccc2)cc1 |wD:4.4,a:9,10,o1:3,10,o2:4,&1:3,9,&2:6|"
    );
    assert_eq!(
        record, before,
        "overlapping enhanced-stereo rendering must preserve its caller"
    );
}

#[test]
fn cx_writer_nonisomeric_mask_removes_stereo_extensions_but_keeps_base_write() {
    // Pinned RDKit 2026.03.1 direct MolToCXSmiles profile. The wrapper removes
    // exactly enhanced stereo and bond configuration fields for
    // doIsomericSmiles=false, while retaining the base output and other masks.
    let mut record = parse_smiles("F[C@](Cl)(Br)I.F[C@](Cl)(Br)I", &Default::default())
        .expect("parse two source stereocenters");
    record.topology.stereo_groups = vec![
        StereoGroup::new(
            StereoGroupKind::And,
            vec![AtomId::new(1), AtomId::new(6)],
            Vec::new(),
        )
        .with_id(7),
    ];
    let before = record.clone();
    let params = CxSmilesWriteParams {
        smiles: SmilesWriteParams {
            canonical: true,
            do_isomeric_smiles: false,
            clean_stereo: true,
            ..Default::default()
        },
        fields: CxSmilesFields::ENHANCED_STEREO | CxSmilesFields::BOND_CFG,
        ..Default::default()
    };

    assert_eq!(
        write_cx_smiles_with_params(&record, &params)
            .map(fixture_writer_text)
            .unwrap(),
        "FC(Cl)(Br)I.FC(Cl)(Br)I"
    );
    assert_eq!(
        record, before,
        "nonisomeric CX writing must preserve its caller"
    );
}

#[test]
fn cx_writer_clean_stereo_matches_pinned_raw_large_ring_cis_output() {
    // Pinned RDKit 2026.03.1 direct MolToCXSmiles legacy profile: raw
    // sanitize=false/removeHs=false input has pending |c:5|, Cis on bond 5
    // with references [4,7], no directions on bonds 4/6, canonical/isomeric/
    // cleanStereo=true, CX_BOND_CFG and default RestoreBondDirOptionClear.
    // The source emits `C1=CCCCCCCCC1`; CX cistrans is removed by clean legacy
    // assignment, not by a serializer special case.
    let parser = cosmolkit_smiles::SmilesParseParams {
        sanitize: false,
        remove_hydrogens: false,
        ..Default::default()
    };
    let record = parse_smiles("C1CCCCC=CCCC1 |c:5|", &parser).expect("parse raw ring CXSMILES");
    assert_eq!(
        record.properties.prop("_needsDetectBondStereo"),
        Some(&PropertyValue::Int(1))
    );
    assert_eq!(record.properties.prop("_StereochemDone"), None);
    assert_eq!(
        record.topology.bonds[5].stereo(),
        cosmolkit_types::BondStereo::Cis
    );
    assert_eq!(
        record.topology.bonds[5].stereo_atoms(),
        Some([AtomId::new(4), AtomId::new(7)])
    );
    assert_eq!(record.topology.bonds[4].direction(), BondDirection::None);
    assert_eq!(record.topology.bonds[6].direction(), BondDirection::None);
    let before = record.clone();
    let params = CxSmilesWriteParams {
        smiles: SmilesWriteParams {
            canonical: true,
            do_isomeric_smiles: true,
            clean_stereo: true,
            ..Default::default()
        },
        fields: CxSmilesFields::BOND_CFG,
        ..Default::default()
    };

    assert_eq!(
        write_cx_smiles_with_params(&record, &params)
            .map(fixture_writer_text)
            .unwrap(),
        "C1=CCCCCCCCC1"
    );
    assert_eq!(record, before, "raw CX writing must preserve its caller");
}

#[test]
fn cx_writer_multiple_component_coordinates_follow_base_output_maps() {
    // Pinned RDKit 2026.03.1 direct MolToCXSmiles profile: canonical writing
    // reorders OC.CO to CO.CO and remaps coordinates using the same flattened
    // atom output order: input rows 1,0,2,3.
    let record = parse_smiles("OC.CO |(0,0,;1,0,;2,0,;3,0,)|", &Default::default())
        .expect("parse two components with 2D coordinates");
    let before = record.clone();
    let params = cx_params(CxSmilesFields::COORDS, true, false);

    assert_eq!(
        write_cx_smiles_with_params(&record, &params)
            .map(fixture_writer_text)
            .unwrap(),
        "CO.CO |(1,0,;0,0,;2,0,;3,0,)|"
    );
    assert_eq!(
        record, before,
        "coordinate-map emission must preserve its caller"
    );
}

#[test]
fn cx_writer_emits_selected_coordinate_dimension_and_rejects_invalid_rows() {
    // Pinned RDKit 2026.03.1 direct MolToCXSmiles get_coords_block: XY is
    // always emitted in atomOrder, Z is emitted only when conf.is3D(), and a
    // formatted zero Z value has no token.
    let mut three_d = parse_smiles("CC", &Default::default()).expect("parse ethane");
    three_d.coordinates = CoordinateBlock {
        conformers_3d: vec![Conformer3D::new(
            17,
            vec![[1.0, 2.0, 9.0], [3.0, 4.0, 0.0]],
            true,
        )],
        ..Default::default()
    };
    let three_d_before = three_d.clone();
    assert_eq!(
        write_cx_smiles_with_params(&three_d, &cx_coordinate_params(CxCoordinateSelection::Auto))
            .map(fixture_writer_text)
            .unwrap(),
        "CC |(1,2,9;3,4,)|".to_owned()
    );
    assert_eq!(three_d, three_d_before);

    let mut xy_only_3d = parse_smiles("CC", &Default::default()).expect("parse ethane");
    xy_only_3d.coordinates = CoordinateBlock {
        conformers_3d: vec![Conformer3D::new(
            18,
            vec![[1.0, 2.0, 9.0], [3.0, 4.0, 8.0]],
            false,
        )],
        ..Default::default()
    };
    let xy_only_before = xy_only_3d.clone();
    assert_eq!(
        write_cx_smiles_with_params(
            &xy_only_3d,
            &cx_coordinate_params(CxCoordinateSelection::Auto)
        )
        .map(fixture_writer_text)
        .unwrap(),
        "CC |(1,2,;3,4,)|".to_owned()
    );
    assert_eq!(xy_only_3d, xy_only_before);

    // Detached CoordinateBlock validation reports malformed source rows before
    // coordinate indexing, retaining the dimension, stored ID and row counts.
    let mut short_2d = parse_smiles("CC", &Default::default()).expect("parse ethane");
    short_2d.coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(19, vec![[1.0, 2.0]])],
        ..Default::default()
    };
    let short_2d_before = short_2d.clone();
    assert_eq!(
        write_cx_smiles_with_params(
            &short_2d,
            &cx_coordinate_params(CxCoordinateSelection::Auto)
        )
        .map(fixture_writer_text)
        .unwrap_err(),
        SmilesParseError::Model("2D conformer 19 has 1 coordinate rows, expected 2".to_owned())
    );
    assert_eq!(short_2d, short_2d_before);

    let mut long_3d = parse_smiles("CC", &Default::default()).expect("parse ethane");
    long_3d.coordinates = CoordinateBlock {
        conformers_3d: vec![Conformer3D::new(
            20,
            vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
            true,
        )],
        ..Default::default()
    };
    let long_3d_before = long_3d.clone();
    assert_eq!(
        write_cx_smiles_with_params(&long_3d, &cx_coordinate_params(CxCoordinateSelection::Auto))
            .map(fixture_writer_text)
            .unwrap_err(),
        SmilesParseError::Model("3D conformer 20 has 3 coordinate rows, expected 2".to_owned())
    );
    assert_eq!(long_3d, long_3d_before);
}

#[test]
fn cx_writer_coordinate_numbers_follow_percent_g_precision_rounding_and_sign_rules() {
    // Pinned RDKit 2026.03.1 CXSmilesOps.cpp get_coords_block uses boost
    // format "%g": six significant digits, fixed/scientific exponent choice,
    // ties-to-even rounding, and zero_small_vals before formatting.
    let cases: [([[f64; 3]; 2], &str); 2] = [
        (
            [
                [0.00123456789, 1.23456789, 100_000.0],
                [1_000_000.0, 999_999.9, -0.0],
            ],
            "CC |(0.00123457,1.23457,100000;1e+06,1e+06,)|",
        ),
        (
            [
                [12_344.249, 12_344.25, 12_344.251],
                [-12_344.249, -12_344.25, -12_344.251],
            ],
            "CC |(12344.2,12344.2,12344.3;-12344.2,-12344.2,-12344.3)|",
        ),
    ];

    for (coordinates, expected) in cases {
        let mut record = parse_smiles("CC", &Default::default()).expect("parse ethane");
        record.coordinates = CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(31, coordinates.to_vec(), true)],
            ..Default::default()
        };
        let before = record.clone();
        assert_eq!(
            write_cx_smiles_with_params(
                &record,
                &cx_coordinate_params(CxCoordinateSelection::Auto)
            )
            .map(fixture_writer_text)
            .unwrap(),
            expected.to_owned()
        );
        assert_eq!(
            record, before,
            "numeric CX rendering must preserve its caller"
        );
    }
}

#[test]
fn cx_writer_zero_small_uses_strict_threshold_for_both_signs() {
    // Pinned RDKit 2026.03.1 CXSmilesOps.cpp::zero_small_vals maps magnitudes
    // below 1e-4 to positive zero, while equality and values above it survive.
    let mut record = parse_smiles("CC", &Default::default()).expect("parse ethane");
    record.coordinates = CoordinateBlock {
        conformers_3d: vec![Conformer3D::new(
            32,
            vec![
                [0.0000999, 0.0001, 0.0001001],
                [-0.0000999, -0.0001, -0.0001001],
            ],
            true,
        )],
        ..Default::default()
    };
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_coordinate_params(CxCoordinateSelection::Auto))
            .map(fixture_writer_text)
            .unwrap(),
        "CC |(0,0.0001,0.0001001;0,-0.0001,-0.0001001)|".to_owned()
    );
    assert_eq!(
        record, before,
        "near-zero rendering must preserve its caller"
    );
}

#[test]
fn cx_writer_bad_bond_config_propagates_source_typed_error_immutably() {
    // Pinned RDKit 2026.03.1 getPropIfPresent<unsigned int> raises
    // RuntimeError("bad any_cast") for this string-valued source property
    // during the default Clear direction phase, before the post-base stage.
    let mut record = parse_smiles("CC", &Default::default()).expect("parse ethane");
    record.topology.bonds[0]
        .set_prop("_MolFileBondCfg", "bad")
        .expect("set malformed source-typed property");
    let before = record.clone();
    let params = cx_params(CxSmilesFields::BOND_CFG, true, false);

    let error = write_cx_smiles_with_params(&record, &params)
        .map(fixture_writer_text)
        .expect_err("source unsigned-int property read must fail visibly");
    assert!(matches!(
        error,
        SmilesParseError::WriterNumeric(cosmolkit_core::PropertyUIntReadError::Lexical { value, .. }) if value.as_bytes() == b"bad"
    ));
    assert_eq!(record, before, "failed CX writing must preserve its caller");
}

#[test]
fn cx_writer_each_source_field_mask_dispatches_its_block() {
    // Pinned RDKit 2026.03.1 direct MolToCXSmiles outputs with canonical,
    // isomeric, cleanStereo=true, doKekule=false and exactly the named field
    // bit. Inputs are identical CX records; each expected string includes the
    // source base output and canonical atom/bond map.
    let cases = [
        (
            "*C |$star;carbon$|",
            CxSmilesFields::ATOM_LABELS,
            "*C |$star;carbon$|",
        ),
        (
            "OC |$_AV:o;c$|",
            CxSmilesFields::MOLFILE_VALUES,
            "CO |$_AV:c;o$|",
        ),
        (
            "OC |(0,0,;1,0,)|",
            CxSmilesFields::COORDS,
            "CO |(1,0,;0,0,)|",
        ),
        ("OC |^1:0|", CxSmilesFields::RADICALS, "C[O] |^1:1|"),
        (
            "OC |atomProp:0.k.v|",
            CxSmilesFields::ATOM_PROPS,
            "CO |atomProp:1.k.v|",
        ),
        (
            "OC1CCC(F)C1 |LN:1:1.3.2.6|",
            CxSmilesFields::LINKNODES,
            "OC1CCC(F)C1 |LN:1:1.3.2.6|",
        ),
        (
            "F[C@H](Cl)Br.O[C@@H](N)I |&7:1,o2:5|",
            CxSmilesFields::ENHANCED_STEREO,
            "F[C@H](Cl)Br.N[C@H](O)I |o1:5,&1:1|",
        ),
        (
            "OCC |SgD:0,2:foo:bar::::|",
            CxSmilesFields::SGROUPS,
            "CCO |SgD:2,0:foo:bar::::|",
        ),
        (
            "OCC |Sg:n:0,1:lab:ht:::|",
            CxSmilesFields::POLYMER,
            "CCO |Sg:n:2,1:lab:ht:::|",
        ),
        ("CC |w:1.0|", CxSmilesFields::BOND_CFG, "CC |w:1.0|"),
        (
            "NOC |C:1.0|",
            CxSmilesFields::COORDINATE_BONDS,
            "C[OH]N |C:1.1|",
        ),
        ("NOC |H:2.1|", CxSmilesFields::HYDROGEN_BONDS, "CON |H:0.0|"),
        ("NOC |Z:0|", CxSmilesFields::ZERO_BONDS, "CO~N |Z:1|"),
    ];

    for (input, fields, expected) in cases {
        let record = parse_smiles(input, &Default::default())
            .unwrap_or_else(|error| panic!("parse source fixture {input:?}: {error}"));
        let before = record.clone();
        assert_eq!(
            write_cx_smiles_with_params(&record, &cx_params(fields, true, false))
                .map(fixture_writer_text)
                .unwrap(),
            expected,
            "CX field mask {} for {input:?}",
            fields.bits()
        );
        assert_eq!(record, before, "field dispatch must preserve {input:?}");
    }

    // Pinned source fixture is sanitized biphenyl with the central single bond
    // marked STEREOATROPCW and stereo references [4, 7]. Its canonical
    // CX_BOND_ATROPISOMER output is `c1ccc(-c2ccccc2)cc1 |wU:3.10|`.
    let mut atrop = parse_smiles_complete_source("c1ccccc1-c1ccccc1", &Default::default())
        .expect("parse biphenyl atropisomer fixture");
    let axis = atrop
        .topology
        .bonds
        .iter_mut()
        .find(|bond| {
            (bond.begin() == AtomId::new(5) && bond.end() == AtomId::new(6))
                || (bond.begin() == AtomId::new(6) && bond.end() == AtomId::new(5))
        })
        .expect("find pinned central biphenyl bond");
    axis.set_stereo(cosmolkit_types::BondStereo::AtropCw)
        .expect("set source atrop stereo");
    axis.set_stereo_atoms(Some([AtomId::new(4), AtomId::new(7)]));
    let before = atrop.clone();

    assert_eq!(
        write_cx_smiles_with_params(
            &atrop,
            &cx_params(CxSmilesFields::BOND_ATROPISOMER, true, false),
        )
        .map(fixture_writer_text)
        .unwrap(),
        "c1ccc(-c2ccccc2)cc1 |wU:3.10|".to_owned()
    );
    assert_eq!(
        atrop, before,
        "atrop field dispatch must preserve its input"
    );
}

#[test]
fn cx_writer_link_nodes_apply_source_validation_and_degree_shapes() {
    // Pinned RDKit 2026.03.1 direct MolToCXSmiles profile: canonical=false,
    // isomeric=true, cleanStereo=false, doKekule=false and CX_LINKNODES.
    // The actual types.h constant is "_molLinkNodes". Keep the original
    // unrecognized "_MolFileLinkNodes" inputs and test their source no-op,
    // then exercise every same payload under the actual source key.
    let params = CxSmilesWriteParams {
        smiles: SmilesWriteParams {
            canonical: false,
            do_isomeric_smiles: true,
            clean_stereo: false,
            do_kekule: false,
            ..Default::default()
        },
        fields: CxSmilesFields::LINKNODES,
        ..Default::default()
    };
    let cases = [
        ("degree two", "CCC", "1 3 2 2 1 2 3", "CCC |LN:1:1.3|"),
        (
            "degree three",
            "CC(C)C",
            "1 4 2 2 1 2 3",
            "CC(C)C |LN:1:1.4.0.2|",
        ),
        (
            "signed unsigned lexical conversion",
            "CCC",
            "-2 4294967295 2 2 1 2 3",
            "CCC |LN:1:4294967294.4294967295|",
        ),
        ("malformed token", "CCC", "1 x 2 2 1 2 3", "CCC"),
        ("fewer than five values", "CCC", "1 3", "CCC"),
        ("missing bond values", "CCC", "1 3 2 2 1", "CCC"),
        ("zero minimum", "CCC", "0 3 2 2 1 2 3", "CCC"),
        ("descending counts", "CCC", "3 1 2 2 1 2 3", "CCC"),
        (
            "unsupported bond count",
            "CC(C)C",
            "1 4 3 2 1 2 3 2 4",
            "CC(C)C",
        ),
        ("different centers", "CC(C)C", "1 4 2 2 1 3 4", "CC(C)C"),
        ("missing declared bond", "CCCC", "1 4 2 2 4 2 3", "CCCC"),
        (
            "first missing bond short circuits later range",
            "CCCC",
            "1 4 2 2 4 2 9",
            "CCCC",
        ),
    ];

    for (name, smiles, link_nodes, expected) in cases {
        let mut record = parse_smiles(smiles, &Default::default())
            .unwrap_or_else(|error| panic!("parse {name} fixture: {error}"));
        record
            .properties
            .set_prop("molFileLinkNodes", link_nodes)
            .expect("set source link-node property");
        let original_before = record.clone();
        let original_output = write_cx_smiles_with_params(&record, &params).unwrap();
        let absent = parse_smiles(smiles, &Default::default()).unwrap();
        assert_eq!(
            original_output,
            write_cx_smiles_with_params(&absent, &params).unwrap(),
            "unrecognized original property key: {name}"
        );
        assert_eq!(record, original_before, "original input preserved: {name}");
        record
            .properties
            .set_prop("_molLinkNodes", link_nodes)
            .expect("set actual source key with same original payload");
        let before = record.clone();
        assert_eq!(
            write_cx_smiles_with_params(&record, &params)
                .map(fixture_writer_text)
                .unwrap(),
            expected,
            "{name}"
        );
        assert_eq!(record, before, "{name} must preserve its caller");
    }
}

#[test]
fn cx_writer_link_nodes_propagate_source_range_failures() {
    // ROMol::getBondBetweenAtoms checks the first endpoint and center before
    // testing edge existence, then reaches the second pair only when the
    // first bond exists. The detached writer projects that range failure into
    // the existing typed CX error without mutating the caller.
    let params = CxSmilesWriteParams {
        smiles: SmilesWriteParams {
            canonical: false,
            clean_stereo: false,
            do_kekule: false,
            ..Default::default()
        },
        fields: CxSmilesFields::LINKNODES,
        ..Default::default()
    };
    for (name, smiles, link_nodes, expected_message) in [
        (
            "second endpoint out of range",
            "CCCC",
            "1 4 2 2 1 2 9",
            "link-node atom index 8 is out of range for 4 atoms",
        ),
        (
            "zero endpoint wraps before range check",
            "CCC",
            "1 3 2 2 0 2 3",
            "link-node atom index 4294967295 is out of range for 3 atoms",
        ),
    ] {
        let mut record = parse_smiles(smiles, &Default::default())
            .unwrap_or_else(|error| panic!("parse {name} fixture: {error}"));
        record
            .properties
            .set_prop("molFileLinkNodes", link_nodes)
            .expect("set source link-node property");
        let original_before = record.clone();
        let original_output = write_cx_smiles_with_params(&record, &params).unwrap();
        let absent = parse_smiles(smiles, &Default::default()).unwrap();
        assert_eq!(
            original_output,
            write_cx_smiles_with_params(&absent, &params).unwrap(),
            "unrecognized original property key: {name}"
        );
        assert_eq!(record, original_before, "original input preserved: {name}");
        record
            .properties
            .set_prop("_molLinkNodes", link_nodes)
            .expect("set actual source key with same original payload");
        let before = record.clone();
        assert_eq!(
            write_cx_smiles_with_params(&record, &params).map(fixture_writer_text),
            Err(SmilesParseError::Cx(expected_message.to_owned())),
            "{name}"
        );
        assert_eq!(record, before, "{name} must preserve its caller");
    }
}

#[test]
fn cx_writer_link_nodes_preserve_direct_canonical_atom_order_mapping() {
    // The pinned serializer computes revOrder but does not use it. These two
    // source records therefore index canonical atomOrder directly by their
    // original center and endpoint indices, while retaining record order.
    let record = parse_smiles("FC1CCC(O)C1 |LN:1:1.3.2.6,4:1.4.3.6|", &Default::default())
        .expect("parse canonical link-node mapping fixture");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_params(CxSmilesFields::LINKNODES, true, false),)
            .map(fixture_writer_text)
            .unwrap(),
        "OC1CCC(F)C1 |LN:4:1.3.3.6,1:1.4.2.6|".to_owned()
    );
    assert_eq!(record, before, "link-node writing must preserve its caller");
}

#[test]
fn cx_writer_data_sgroups_preserve_typed_members_fields_and_group_order() {
    // Pinned RDKit 2026.03.1 direct canonical/isomeric, cleanStereo=false,
    // doKekule=false, CX_SGROUPS output. DAT group order and each group's
    // member order are retained while members map through canonical revOrder.
    let mut record = parse_smiles("OCC", &Default::default()).expect("parse DAT fixture");
    let mut first = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
        .with_atoms(vec![AtomId::new(2), AtomId::new(0)])
        .with_data(SGroupData {
            field_name: Some("first".into()),
            field_info: Some("i,n".into()),
            query_op: Some("q:o".into()),
            values: vec!["v1".into()],
            ..SGroupData::default()
        });
    first.set_prop("FIELDTAG", "t|g");
    let second = SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
        .with_atoms(vec![AtomId::new(1)])
        .with_data(SGroupData {
            field_name: Some("second".into()),
            values: vec!["v2".into()],
            ..SGroupData::default()
        });
    record.topology.substance_groups = vec![first, second];
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_data_sgroup_params(true))
            .map(fixture_writer_text)
            .unwrap(),
        "CCO |SgD:0,2:first:v1:q:o:i,n:t|g:,SgD:1:second:v2::::|"
    );
    assert_eq!(record, before, "DAT writing must preserve its caller");
}

#[test]
fn cx_writer_data_sgroups_emit_all_typed_values_and_empty_members() {
    // A pinned V3000 DAT group with DATAFIELDS ["first", "second"] emits
    // both values in vector order. A source DAT group with no atoms retains
    // get_sgroup_data_block's surprising seekp result instead of disappearing.
    let mut values = parse_smiles("CO", &Default::default()).expect("parse multi-value DAT");
    values.topology.substance_groups = vec![
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_atoms(vec![AtomId::new(1), AtomId::new(0)])
            .with_data(SGroupData {
                field_name: Some("FIELD".into()),
                field_info: Some("INFO".into()),
                query_op: Some("=".into()),
                values: vec!["first".into(), "second".into()],
                ..SGroupData::default()
            }),
    ];
    let before_values = values.clone();
    assert_eq!(
        write_cx_smiles_with_params(&values, &cx_data_sgroup_params(false))
            .map(fixture_writer_text)
            .unwrap(),
        "CO |SgD:1,0:FIELD:first,second:=:INFO::|"
    );
    assert_eq!(values, before_values, "multi-value DAT must be immutable");

    let mut empty = parse_smiles("CC", &Default::default()).expect("parse empty-member DAT");
    empty.topology.substance_groups = vec![
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data).with_data(
            SGroupData {
                field_name: Some("FIELD".into()),
                values: vec!["VALUE".into()],
                ..SGroupData::default()
            },
        ),
    ];
    let before_empty = empty.clone();
    assert_eq!(
        write_cx_smiles_with_params(&empty, &cx_data_sgroup_params(false))
            .map(fixture_writer_text)
            .unwrap(),
        "CC |SgD:FIELD:VALUE::::|"
    );
    assert_eq!(empty, before_empty, "empty-member DAT must be immutable");
}

#[test]
fn cx_writer_data_sgroups_do_not_requote_source_literal_field_bytes() {
    // read_text_to decodes numeric entities, while get_sgroup_data_block
    // appends the resulting bytes directly and never calls quote_string.
    // The deliberately ambiguous output below is the pinned source result.
    let record = parse_smiles(
        "CC |SgD:0,1:a&#46;b:c&#59;d:q&#58;o:i&#44;n:t&#124;g:|",
        &Default::default(),
    )
    .expect("parse escaped DAT fixture");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_data_sgroup_params(false))
            .map(fixture_writer_text)
            .unwrap(),
        "CC |SgD:0,1:a.b:c;d:q:o:i,n:t|g:|"
    );
    assert_eq!(record, before, "literal DAT writing must be immutable");
}

#[test]
fn cx_writer_polymer_sgroups_map_all_supported_types_in_source_order() {
    // Pinned RDKit 2026.03.1 direct noncanonical/isomeric,
    // cleanStereo=false, doKekule=false, CX_POLYMER output. reverseTypemap
    // keeps the lexicographically first code for each TYPE, so plain COP and
    // COP/ALT both emit `alt`; canonical MIX and COM emit `mix` and `c`.
    let mut record = parse_smiles("C", &Default::default()).expect("parse polymer fixture");
    let kinds = [
        (SubstanceGroupKind::StructuralRepeatUnit, None, "n"),
        (SubstanceGroupKind::Monomer, None, "mon"),
        (SubstanceGroupKind::Mer, None, "mer"),
        (SubstanceGroupKind::Copolymer, None, "alt"),
        (SubstanceGroupKind::Copolymer, Some("ALT"), "alt"),
        (SubstanceGroupKind::Copolymer, Some("RAN"), "ran"),
        (SubstanceGroupKind::Copolymer, Some("BLO"), "blk"),
        (SubstanceGroupKind::Crosslink, None, "xl"),
        (SubstanceGroupKind::Modification, None, "mod"),
        (SubstanceGroupKind::Mixture, None, "mix"),
        (SubstanceGroupKind::Formulation, None, "f"),
        (SubstanceGroupKind::AnyPolymer, None, "any"),
        (SubstanceGroupKind::Generic("GEN".into()), None, "gen"),
        (SubstanceGroupKind::MixtureComponent, None, "c"),
        (SubstanceGroupKind::Graft, None, "grf"),
    ];
    record.topology.substance_groups = kinds
        .into_iter()
        .enumerate()
        .map(|(index, (kind, subtype, _))| {
            let group = SubstanceGroup::new(SubstanceGroupId::new(index), kind)
                .with_atoms(vec![AtomId::new(0)])
                .with_label(format!("L{index}"))
                .with_connection(SGroupConnection::Either);
            match subtype {
                Some(subtype) => group.with_subtype(subtype),
                None => group,
            }
        })
        .collect();
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_polymer_params(false))
            .map(fixture_writer_text)
            .unwrap(),
        "C |Sg:n:0:L0:eu:::,Sg:mon:0:L1:eu:::,Sg:mer:0:L2:eu:::,Sg:alt:0:L3:eu:::,Sg:alt:0:L4:eu:::,Sg:ran:0:L5:eu:::,Sg:blk:0:L6:eu:::,Sg:xl:0:L7:eu:::,Sg:mod:0:L8:eu:::,Sg:mix:0:L9:eu:::,Sg:f:0:L10:eu:::,Sg:any:0:L11:eu:::,Sg:gen:0:L12:eu:::,Sg:c:0:L13:eu:::,Sg:grf:0:L14:eu:::|"
    );
    assert_eq!(record, before, "polymer type writing must be immutable");
}

#[test]
fn cx_writer_polymer_sgroups_preserve_source_field_and_connection_order() {
    // Pinned profile as above. Members, groups and labels retain source order;
    // CONNECT is lowercased in place after the label field.
    let mut record = parse_smiles("OC", &Default::default()).expect("parse connectivity fixture");
    record.topology.substance_groups = vec![
        SubstanceGroup::new(
            SubstanceGroupId::new(0),
            SubstanceGroupKind::StructuralRepeatUnit,
        )
        .with_atoms(vec![AtomId::new(1), AtomId::new(0)])
        .with_label("head")
        .with_connection(SGroupConnection::HeadToHead),
        SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Monomer)
            .with_atoms(vec![AtomId::new(0)])
            .with_label("tail")
            .with_connection(SGroupConnection::HeadToTail),
        SubstanceGroup::new(SubstanceGroupId::new(2), SubstanceGroupKind::Mer)
            .with_atoms(vec![AtomId::new(1)])
            .with_label("either")
            .with_connection(SGroupConnection::Either),
        SubstanceGroup::new(SubstanceGroupId::new(3), SubstanceGroupKind::Graft)
            .with_atoms(vec![AtomId::new(0)])
            .with_label("unknown")
            .with_connection(SGroupConnection::Unknown("MiXeD".into())),
    ];
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_polymer_params(false))
            .map(fixture_writer_text)
            .unwrap(),
        "OC |Sg:n:1,0:head:hh:::,Sg:mon:0:tail:ht:::,Sg:mer:1:either:eu:::,Sg:grf:0:unknown:mixed:::|"
    );
    assert_eq!(record, before, "polymer field writing must be immutable");
}

#[test]
fn cx_writer_polymer_crossings_use_direct_bond_order_index_and_pair_order() {
    // Pinned RDKit 2026.03.1 canonical/isomeric, cleanStereo=false,
    // doKekule=false, CX_POLYMER output. The canonical bond-order vector is
    // the non-self-inverse cycle [4,0,1,2,3]. get_sgroup_polymer_block uses
    // bondOrder[v] directly and emits odd XBCORR entries in pair order.
    let record = parse_smiles(
        "C(COCN)C |Sg:n:0,1,2,3,4,5::ht:0,1,2,3,4:4,3,2,1,0:|",
        &Default::default(),
    )
    .expect("parse non-self-inverse polymer crossing fixture");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_polymer_params(true))
            .map(fixture_writer_text)
            .unwrap(),
        "CCCOCN |Sg:n:1,2,3,4,5,0::ht:4,0,1,2,3:3,2,1,0,4:|"
    );
    assert_eq!(record, before, "crossing-bond writing must be immutable");
}

#[test]
fn cx_writer_polymer_crossings_apply_source_length_thresholds() {
    // Pinned source basic/minimal row: finalized XBHEAD [0] and XBCORR [0,3]
    // do not pass the strict size>1 and size>2 writer guards. Construct the
    // equivalent typed owner state directly because star_e is query-only and
    // deliberately rejected by the protected concrete-molecule lowerer.
    let mut record = parse_smiles("CCCO", &Default::default())
        .expect("parse concrete minimal polymer crossing fixture");
    record.topology.substance_groups = vec![
        SubstanceGroup::new(
            SubstanceGroupId::new(0),
            SubstanceGroupKind::StructuralRepeatUnit,
        )
        .with_atoms(vec![
            AtomId::new(0),
            AtomId::new(1),
            AtomId::new(2),
            AtomId::new(3),
        ])
        .with_connection(SGroupConnection::Either)
        .with_head_crossing_bonds(vec![BondId::new(0)])
        .with_crossing_bond_correspondence(vec![BondId::new(0), BondId::new(2)]),
    ];
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_polymer_params(false))
            .map(fixture_writer_text)
            .unwrap(),
        "CCCO |Sg:n:0,1,2,3::eu:::|"
    );
    assert_eq!(record, before, "threshold writing must be immutable");
}

#[test]
fn cx_writer_sgroup_hierarchy_maps_selected_parent_children_in_source_order() {
    // Pinned RDKit 2026.03.1 direct MolToCXSmiles profile: noncanonical,
    // isomeric, cleanStereo=false, doKekule=false. Output indices are assigned
    // to DAT rows first and polymer rows second. The SUP row is excluded, and
    // child indices retain SGroup source-row order across the two families.
    let mut record =
        parse_smiles("CCC", &Default::default()).expect("parse concrete mixed hierarchy fixture");
    record.topology.substance_groups = vec![
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_atoms(vec![AtomId::new(0)])
            .with_parent(SubstanceGroupId::new(3))
            .with_data(SGroupData {
                field_name: Some("A".into()),
                ..Default::default()
            })
            .with_prop("index", "40")
            .unwrap(),
        SubstanceGroup::new(
            SubstanceGroupId::new(1),
            SubstanceGroupKind::StructuralRepeatUnit,
        )
        .with_atoms(vec![AtomId::new(1)])
        .with_parent(SubstanceGroupId::new(3))
        .with_connection(SGroupConnection::HeadToTail)
        .with_prop("index", "7")
        .unwrap(),
        SubstanceGroup::new(SubstanceGroupId::new(2), SubstanceGroupKind::Superatom)
            .with_atoms(vec![AtomId::new(2)])
            .with_parent(SubstanceGroupId::new(3))
            .with_prop("index", "8")
            .unwrap(),
        SubstanceGroup::new(
            SubstanceGroupId::new(3),
            SubstanceGroupKind::StructuralRepeatUnit,
        )
        .with_atoms(vec![AtomId::new(2)])
        .with_connection(SGroupConnection::HeadToTail)
        .with_prop("index", "900")
        .unwrap(),
        SubstanceGroup::new(SubstanceGroupId::new(4), SubstanceGroupKind::Data)
            .with_atoms(vec![AtomId::new(1)])
            .with_parent(SubstanceGroupId::new(3))
            .with_data(SGroupData {
                field_name: Some("C".into()),
                ..Default::default()
            })
            .with_prop("index", "42")
            .unwrap(),
    ];
    let before = record.clone();

    for (params, expected) in [
        (
            CxSmilesWriteParams {
                fields: CxSmilesFields::SGROUPS | CxSmilesFields::POLYMER,
                ..cx_data_sgroup_params(false)
            },
            "CCC |SgD:0:A:::::,SgD:1:C:::::,Sg:n:1::ht:::,Sg:n:2::ht:::,SgH:3:0.2.1|",
        ),
        (
            cx_data_sgroup_params(false),
            "CCC |SgD:0:A:::::,SgD:1:C:::::|",
        ),
        (
            cx_polymer_params(false),
            "CCC |Sg:n:1::ht:::,Sg:n:2::ht:::,SgH:1:0|",
        ),
    ] {
        assert_eq!(
            write_cx_smiles_with_params(&record, &params)
                .map(fixture_writer_text)
                .unwrap(),
            expected
        );
        assert_eq!(record, before, "hierarchy writing must be immutable");
    }
}

#[test]
fn cx_writer_emits_source_field_order_and_comma_pipe_separators() {
    // Pinned RDKit 2026.03.1 direct output for coordinates, labels, molfile
    // values and atom props: each nonempty block follows source order with one
    // comma, and the combined extension has exactly one opening/closing pipe.
    let record = parse_smiles(
        "OC |(0,0,;1,0,),$oxygen;carbon$,$_AV:o;c$,atomProp:0.k.v|",
        &Default::default(),
    )
    .expect("parse multi-field CX fixture");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(
            &record,
            &cx_params(
                CxSmilesFields::COORDS
                    | CxSmilesFields::ATOM_LABELS
                    | CxSmilesFields::MOLFILE_VALUES
                    | CxSmilesFields::ATOM_PROPS,
                true,
                false,
            ),
        )
        .map(fixture_writer_text)
        .unwrap(),
        "CO |(1,0,;0,0,),$carbon;oxygen$,$_AV:c;o$,atomProp:1.k.v|".to_owned()
    );
    assert_eq!(
        record, before,
        "combined field emission must preserve input"
    );

    // The pinned hierarchy helper is reached for either SGroups or POLYMER
    // selection and emits only the selected family of parent/child groups.
    for (input, fields, expected) in [
        (
            "CC |SgD:0:PARENT:p::::,SgD:1:CHILD:c::::,SgH:0:1|",
            CxSmilesFields::SGROUPS,
            "CC |SgD:0:PARENT:p::::,SgD:1:CHILD:c::::,SgH:0:1|",
        ),
        (
            "CCC |Sg:n:0,1::ht:::,Sg:n:1,2::ht:::,SgH:0:1|",
            CxSmilesFields::POLYMER,
            "CCC |Sg:n:0,1::ht:::,Sg:n:1,2::ht:::,SgH:0:1|",
        ),
    ] {
        let record = parse_smiles(input, &Default::default())
            .unwrap_or_else(|error| panic!("parse hierarchy fixture {input:?}: {error}"));
        assert_eq!(
            write_cx_smiles_with_params(&record, &cx_params(fields, true, false))
                .map(fixture_writer_text)
                .unwrap(),
            expected,
            "hierarchy dispatch for {}",
            fields.bits()
        );
    }
}

#[test]
fn cx_writer_nonisomeric_mask_removes_only_stereo_field_bits() {
    // Pinned RDKit 2026.03.1 direct MolToCXSmiles profile. Source input has
    // enhanced stereo, a wiggly bond, atom labels and an atom property.
    // Nonisomeric output clears CX_ENHANCEDSTEREO and CX_BOND_CFG while the
    // labels and atom property remain emitted under the same CX_ALL mask.
    let record = parse_smiles(
        "F[C@H](Cl)Br.CC |&7:1,w:4.3,atomProp:0.keep.yes,$f;chiral;cl;br;c;c$|",
        &Default::default(),
    )
    .expect("parse isomeric field-mask fixture");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_params(CxSmilesFields::ALL, true, false))
            .map(fixture_writer_text)
            .unwrap(),
        "CC.F[C@H](Cl)Br |$c;c;f;chiral;cl;br$,atomProp:2.keep.yes,w:0.0,&1:3|"
    );
    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_params(CxSmilesFields::ALL, false, false))
            .map(fixture_writer_text)
            .unwrap(),
        "CC.FC(Cl)Br |$c;c;f;chiral;cl;br$,atomProp:2.keep.yes|".to_owned()
    );
    assert_eq!(
        record, before,
        "isomeric mask projection must preserve input"
    );
}

#[test]
fn cx_writer_atom_labels_follow_pseudo_attachment_and_fallback_order() {
    // Pinned RDKit 2026.03.1 get_atomlabel_block accepts only Pol/Mod as
    // dummyLabel pseudoatoms, then checks attachment points and atomLabel.
    let params = cx_atom_label_params();
    let cases: [(Option<&str>, Option<&str>, Option<&str>, Option<&str>, &str); 8] = [
        (Some("A"), None, None, None, "*"),
        (Some("Pol"), None, None, None, "* |$Pol_p$|"),
        (Some("Mod"), None, None, None, "* |$Mod_p$|"),
        (
            Some("A"),
            Some("generic"),
            Some("1"),
            Some("ordinary"),
            "* |$generic_p$|",
        ),
        (Some("A"), None, None, Some("ordinary"), "* |$ordinary$|"),
        (Some("A"), None, Some("1"), Some("ordinary"), "* |$_AP1$|"),
        (Some("A"), None, Some("2"), Some("ordinary"), "* |$_AP2$|"),
        (
            Some("A"),
            None,
            Some("3"),
            Some("ordinary"),
            "* |$ordinary$|",
        ),
    ];

    for (dummy, generic, attachment, ordinary, expected) in cases {
        let mut record = parse_smiles("*", &Default::default()).expect("parse dummy atom");
        if let Some(value) = dummy {
            record.topology.atoms[0]
                .set_prop("dummyLabel", value)
                .expect("set dummy label");
        }
        if let Some(value) = generic {
            record.topology.atoms[0]
                .set_prop("_QueryAtomGenericLabel", value)
                .expect("set generic query label");
        }
        if let Some(value) = attachment {
            record.topology.atoms[0]
                .set_prop("_fromAttchpt", value)
                .expect("set attachment point");
        }
        if let Some(value) = ordinary {
            record.topology.atoms[0]
                .set_prop("atomLabel", value)
                .expect("set ordinary atom label");
        }
        let before = record.clone();

        let original_expected = if attachment.is_some() && generic.is_none() {
            if ordinary.is_some() {
                "* |$ordinary$|"
            } else {
                "*"
            }
        } else {
            expected
        };
        assert_eq!(
            write_cx_smiles_with_params(&record, &params)
                .map(fixture_writer_text)
                .unwrap(),
            expected,
            "dummy={dummy:?}, generic={generic:?}, attachment={attachment:?}, ordinary={ordinary:?}"
        );
        assert_eq!(record, before, "CX label writing must preserve its caller");
        if let Some(value) = attachment {
            let mut source_record = record.clone();
            source_record.topology.atoms[0]
                .set_prop("_fromAttchpt", value)
                .expect("same original attachment payload at actual source key");
            let source_before = source_record.clone();
            assert_eq!(
                write_cx_smiles_with_params(&source_record, &params).unwrap(),
                expected.into()
            );
            assert_eq!(
                source_record, source_before,
                "source attachment fixture remains immutable"
            );
        }
    }
}

#[test]
fn cx_writer_atom_label_holes_trailing_slots_and_text_follow_source() {
    // Pinned RDKit quote_string returns text unchanged. Empty mapped positions
    // remain semicolon slots, including the final empty atom slot.
    let params = cx_atom_label_params();
    let mut record = parse_smiles("CCCC", &Default::default()).expect("parse butane");
    record.topology.atoms[0]
        .set_prop("atomLabel", "left.with.dot")
        .expect("set dotted source label");
    record.topology.atoms[2]
        .set_prop("atomLabel", "tail")
        .expect("set trailing source label");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &params)
            .map(fixture_writer_text)
            .unwrap(),
        "CCCC |$left.with.dot;;tail;$|"
    );
    assert_eq!(record, before, "CX label slots must preserve their caller");

    let mut leading_hole = parse_smiles("CCCC", &Default::default()).expect("parse butane");
    leading_hole.topology.atoms[1]
        .set_prop("atomLabel", "middle")
        .expect("set middle source label");
    assert_eq!(
        write_cx_smiles_with_params(&leading_hole, &params)
            .map(fixture_writer_text)
            .unwrap(),
        "CCCC |$;middle;;$|"
    );
}

#[test]
fn cx_writer_atom_label_clears_semicolon_only_source_block() {
    // The pinned final find_if_not clears a block made only from separators,
    // including a literal semicolon returned unchanged by quote_string.
    let mut record = parse_smiles("C", &Default::default()).expect("parse carbon");
    record.topology.atoms[0]
        .set_prop("atomLabel", ";")
        .expect("set semicolon label");
    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_atom_label_params())
            .map(fixture_writer_text)
            .unwrap(),
        "C"
    );
}

#[test]
fn cx_writer_mixed_dimension_selection_is_explicit_for_both_record_orders() {
    // Pinned RDKit CX parsing assigns IDs in source record order. Its implicit
    // first-conformer outputs are retained here as reference evidence:
    // 2D-first -> `CC |(0,0,;1,0,)|`; 3D-first -> `CC |(0,0,1;1,0,2)|`.
    // CK-COORD-002 deliberately makes Auto ambiguous for both inputs; selecting
    // either stored dimension and ID preserves that conformer's source output.
    let cases = [
        (
            "CC |(0,0;1,0)(0,0,1;1,0,2)|",
            0,
            1,
            true,
            "CC |(0,0,;1,0,)|",
        ),
        (
            "CC |(0,0,1;1,0,2)(0,0;1,0)|",
            1,
            0,
            false,
            "CC |(0,0,1;1,0,2)|",
        ),
    ];

    for (input, expected_2d_id, expected_3d_id, pinned_default_is_2d, pinned_default_output) in
        cases
    {
        let record = parse_smiles(input, &Default::default())
            .unwrap_or_else(|error| panic!("parse mixed-dimension input {input:?}: {error}"));
        assert_eq!(record.coordinates.conformers_2d.len(), 0, "{input}");
        assert_eq!(record.coordinates.conformers_3d.len(), 2, "{input}");
        let false_3d = record
            .coordinates
            .conformers_3d
            .iter()
            .find(|c| !c.is_3d())
            .unwrap();
        let true_3d = record
            .coordinates
            .conformers_3d
            .iter()
            .find(|c| c.is_3d())
            .unwrap();
        let two_d_id = false_3d.id();
        let three_d_id = true_3d.id();
        assert_eq!(false_3d.coordinates(), &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        assert_eq!(true_3d.coordinates(), &[[0.0, 0.0, 1.0], [1.0, 0.0, 2.0]]);
        assert_eq!(
            two_d_id, expected_2d_id,
            "source coordinate record ID: {input}"
        );
        assert_eq!(
            three_d_id, expected_3d_id,
            "source coordinate record ID: {input}"
        );
        let before = record.clone();

        assert_eq!(
            write_cx_smiles_with_params(
                &record,
                &cx_coordinate_params(CxCoordinateSelection::Auto)
            )
            .map(fixture_writer_text),
            Err(SmilesParseError::AmbiguousCoordinateSelection {
                two_d_count: 0,
                three_d_count: 2
            })
        );
        // Native Point3D with is3D=false remains a ThreeD storage set.
        assert!(matches!(write_cx_smiles_with_params(&record,
            &cx_coordinate_params(CxCoordinateSelection::TwoD { id: two_d_id })),
            Err(SmilesParseError::MissingCoordinateSelection {
                selection: CxCoordinateSelection::TwoD { id }
            }) if id == two_d_id));
        let explicit_2d = write_cx_smiles_with_params(
            &record,
            &cx_coordinate_params(CxCoordinateSelection::ThreeD { id: two_d_id }),
        )
        .map(fixture_writer_text)
        .unwrap();
        let explicit_3d = write_cx_smiles_with_params(
            &record,
            &cx_coordinate_params(CxCoordinateSelection::ThreeD { id: three_d_id }),
        )
        .map(fixture_writer_text)
        .unwrap();

        assert_eq!(explicit_2d, "CC |(0,0,;1,0,)|".to_owned());
        assert_eq!(explicit_3d, "CC |(0,0,1;1,0,2)|".to_owned());
        assert_eq!(
            if pinned_default_is_2d {
                std::str::from_utf8(explicit_2d.as_bytes())
                    .expect("original explicit coordinate fixture is UTF-8")
            } else {
                std::str::from_utf8(explicit_3d.as_bytes())
                    .expect("original explicit coordinate fixture is UTF-8")
            },
            pinned_default_output,
            "the pinned first-insertion output remains represented by explicit selection"
        );
        assert_eq!(record, before, "coordinate selection must preserve {input}");
        // A caller-owned typed 2D set is distinct from native Point3D storage
        // whose is3D flag is false. Exercise the explicit typed route as well.
        let mut typed_two_d = record.clone();
        typed_two_d.coordinates = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(two_d_id, vec![[0.0, 0.0], [1.0, 0.0]])],
            ..Default::default()
        };
        let typed_before = typed_two_d.clone();
        assert_eq!(
            write_cx_smiles_with_params(
                &typed_two_d,
                &cx_coordinate_params(CxCoordinateSelection::TwoD { id: two_d_id }),
            )
            .map(fixture_writer_text)
            .unwrap(),
            explicit_2d,
        );
        assert_eq!(
            typed_two_d, typed_before,
            "explicit typed 2D fixture preserved"
        );
    }
}

#[test]
fn cx_writer_selection_uses_dimension_scoped_ids_and_rejects_nonunique_or_missing_sets() {
    let empty = parse_smiles("CC", &Default::default()).expect("parse coordinate-free ethane");
    let empty_before = empty.clone();
    assert_eq!(
        write_cx_smiles_with_params(&empty, &cx_coordinate_params(CxCoordinateSelection::Auto))
            .map(fixture_writer_text)
            .unwrap(),
        "CC".to_owned()
    );
    assert!(matches!(
        write_cx_smiles_with_params(
            &empty,
            &cx_coordinate_params(CxCoordinateSelection::TwoD { id: 99 })
        )
        .map(fixture_writer_text),
        Err(SmilesParseError::MissingCoordinateSelection {
            selection: CxCoordinateSelection::TwoD { id: 99 }
        })
    ));
    assert_eq!(empty, empty_before);

    let sole = parse_smiles("CC |(4,5;6,7)|", &Default::default())
        .expect("parse one stored 2D coordinate set");
    let sole_before = sole.clone();
    assert_eq!(
        write_cx_smiles_with_params(&sole, &cx_coordinate_params(CxCoordinateSelection::Auto))
            .map(fixture_writer_text)
            .unwrap(),
        "CC |(4,5,;6,7,)|".to_owned()
    );
    assert_eq!(sole, sole_before);

    let mut sole_3d = parse_smiles("CC", &Default::default()).expect("parse ethane for 3D");
    sole_3d.coordinates = CoordinateBlock {
        conformers_3d: vec![Conformer3D::new(
            23,
            vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
            true,
        )],
        ..Default::default()
    };
    let sole_3d_before = sole_3d.clone();
    assert_eq!(
        write_cx_smiles_with_params(&sole_3d, &cx_coordinate_params(CxCoordinateSelection::Auto))
            .map(fixture_writer_text)
            .unwrap(),
        "CC |(1,2,3;4,5,6)|".to_owned()
    );
    assert_eq!(
        write_cx_smiles_with_params(
            &sole_3d,
            &cx_coordinate_params(CxCoordinateSelection::ThreeD { id: 23 })
        )
        .map(fixture_writer_text)
        .unwrap(),
        "CC |(1,2,3;4,5,6)|".to_owned()
    );
    assert!(matches!(
        write_cx_smiles_with_params(
            &sole_3d,
            &cx_coordinate_params(CxCoordinateSelection::TwoD { id: 23 })
        )
        .map(fixture_writer_text),
        Err(SmilesParseError::MissingCoordinateSelection {
            selection: CxCoordinateSelection::TwoD { id: 23 }
        })
    ));
    assert!(matches!(
        write_cx_smiles_with_params(
            &sole_3d,
            &cx_coordinate_params(CxCoordinateSelection::ThreeD { id: 24 })
        )
        .map(fixture_writer_text),
        Err(SmilesParseError::MissingCoordinateSelection {
            selection: CxCoordinateSelection::ThreeD { id: 24 }
        })
    ));
    assert_eq!(sole_3d, sole_3d_before);

    let mut multiple_2d = parse_smiles("CC", &Default::default()).expect("parse ethane");
    multiple_2d.coordinates.conformers_2d = vec![
        Conformer2D::new(7, vec![[0.0, 0.0], [1.0, 0.0]]),
        Conformer2D::new(21, vec![[2.0, 0.0], [3.0, 0.0]]),
    ];
    let multiple_before = multiple_2d.clone();
    assert_eq!(
        write_cx_smiles_with_params(
            &multiple_2d,
            &cx_coordinate_params(CxCoordinateSelection::Auto)
        )
        .map(fixture_writer_text),
        Err(SmilesParseError::AmbiguousCoordinateSelection {
            two_d_count: 2,
            three_d_count: 0
        })
    );
    assert_eq!(
        write_cx_smiles_with_params(
            &multiple_2d,
            &cx_coordinate_params(CxCoordinateSelection::TwoD { id: 21 })
        )
        .map(fixture_writer_text)
        .unwrap(),
        "CC |(2,0,;3,0,)|".to_owned()
    );
    assert!(matches!(
        write_cx_smiles_with_params(
            &multiple_2d,
            &cx_coordinate_params(CxCoordinateSelection::TwoD { id: 1 })
        )
        .map(fixture_writer_text),
        Err(SmilesParseError::MissingCoordinateSelection {
            selection: CxCoordinateSelection::TwoD { id: 1 }
        })
    ));
    assert_eq!(multiple_2d, multiple_before);

    let mut same_id_both_dimensions =
        parse_smiles("CC", &Default::default()).expect("parse ethane for dimension-scoped IDs");
    same_id_both_dimensions.coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(5, vec![[10.0, 20.0], [30.0, 40.0]])],
        conformers_3d: vec![Conformer3D::new(
            5,
            vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
            true,
        )],
        ..Default::default()
    };
    let same_id_before = same_id_both_dimensions.clone();
    assert!(matches!(
        write_cx_smiles_with_params(
            &same_id_both_dimensions,
            &cx_coordinate_params(CxCoordinateSelection::Auto)
        )
        .map(fixture_writer_text),
        Err(SmilesParseError::AmbiguousCoordinateSelection {
            two_d_count: 1,
            three_d_count: 1
        })
    ));
    assert_eq!(
        write_cx_smiles_with_params(
            &same_id_both_dimensions,
            &cx_coordinate_params(CxCoordinateSelection::TwoD { id: 5 })
        )
        .map(fixture_writer_text)
        .unwrap(),
        "CC |(10,20,;30,40,)|".to_owned()
    );
    assert_eq!(
        write_cx_smiles_with_params(
            &same_id_both_dimensions,
            &cx_coordinate_params(CxCoordinateSelection::ThreeD { id: 5 })
        )
        .map(fixture_writer_text)
        .unwrap(),
        "CC |(1,2,3;4,5,6)|".to_owned()
    );
    assert_eq!(same_id_both_dimensions, same_id_before);
}

#[test]
fn cx_writer_atom_values_preserve_mapped_holes_and_literal_text() {
    // Pinned get_value_block copies each present value literally and retains
    // one semicolon slot for every atom in the base writer's output order.
    let mut record = parse_smiles("CCCC", &Default::default()).expect("parse butane");
    record.topology.atoms[0]
        .set_prop("molFileValue", "left.with.dot")
        .expect("set first molfile value");
    record.topology.atoms[2]
        .set_prop("molFileValue", "tail")
        .expect("set third molfile value");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_atom_value_params())
            .map(fixture_writer_text)
            .unwrap(),
        "CCCC |$_AV:left.with.dot;;tail;$|"
    );
    assert_eq!(record, before, "value output must preserve its input");
}

#[test]
fn cx_writer_atom_values_preserve_leading_interior_and_trailing_holes() {
    // The pinned block has an empty value slot for every atom lacking the
    // property, including leading, interior and trailing positions.
    let mut record = parse_smiles("CCCC", &Default::default()).expect("parse butane");
    record.topology.atoms[1]
        .set_prop("molFileValue", "middle")
        .expect("set second molfile value");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_atom_value_params())
            .map(fixture_writer_text)
            .unwrap(),
        "CCCC |$_AV:;middle;;$|"
    );
    assert_eq!(record, before, "value output must preserve its input");
}

#[test]
fn cx_writer_atom_values_emit_present_empty_values_and_all_empty_slots() {
    // getCXExtensions gates the block on property presence, while
    // get_value_block preserves every slot even when all present text is empty.
    let mut single = parse_smiles("C", &Default::default()).expect("parse carbon");
    single.topology.atoms[0]
        .set_prop("molFileValue", "")
        .expect("set present empty molfile value");
    let single_before = single.clone();
    assert_eq!(
        write_cx_smiles_with_params(&single, &cx_atom_value_params())
            .map(fixture_writer_text)
            .unwrap(),
        "C |$_AV:$|"
    );
    assert_eq!(single, single_before, "empty value output preserves input");

    let mut multiple = parse_smiles("CCC", &Default::default()).expect("parse propane");
    multiple.topology.atoms[1]
        .set_prop("molFileValue", "")
        .expect("set present empty molfile value");
    let multiple_before = multiple.clone();
    assert_eq!(
        write_cx_smiles_with_params(&multiple, &cx_atom_value_params())
            .map(fixture_writer_text)
            .unwrap(),
        "CCC |$_AV:;;$|"
    );
    assert_eq!(multiple, multiple_before, "empty slots preserve input");
}

#[test]
fn cx_writer_radicals_group_supported_counts_and_keep_output_atom_order() {
    // Pinned get_radical_block orders groups by electron count and positions
    // within each group by the base writer's atom output order.
    let mut record = parse_smiles("CCCCC", &Default::default()).expect("parse pentane");
    for (atom, count) in [2, 1, 2, 3, 1].into_iter().enumerate() {
        record.topology.atoms[atom].set_radical_electrons(count);
    }
    let before = record.clone();
    let params = cx_radical_params();
    let base = write_smiles_with_params(&record, &params.smiles)
        .map(fixture_writer_text)
        .unwrap();

    assert_eq!(
        write_cx_smiles_with_params(&record, &params)
            .map(fixture_writer_text)
            .unwrap(),
        format!("{base} |^1:1,4,^2:0,2,^5:3|")
    );
    assert_eq!(record, before, "radical extension writing preserves input");
}

#[test]
fn cx_writer_radicals_preserve_source_default_for_unsupported_counts() {
    // Pinned source warns for count 4 but still appends its output index with
    // no marker; the CX caller serializes this nonempty source block.
    let mut record = parse_smiles("CCC", &Default::default()).expect("parse propane");
    record.topology.atoms[0].set_radical_electrons(4);
    record.topology.atoms[1].set_radical_electrons(1);
    record.topology.atoms[2].set_radical_electrons(2);
    let before = record.clone();
    let params = cx_radical_params();
    let base = write_smiles_with_params(&record, &params.smiles)
        .map(fixture_writer_text)
        .unwrap();

    assert_eq!(
        write_cx_smiles_with_params(&record, &params)
            .map(fixture_writer_text)
            .unwrap(),
        format!("{base} |^1:1,^2:2,0|")
    );
    assert_eq!(record, before, "radical extension writing preserves input");
}

#[test]
fn cx_writer_radicals_emit_no_block_when_all_electron_counts_are_zero() {
    let record = parse_smiles("CCC", &Default::default()).expect("parse propane");
    let params = cx_radical_params();
    let base = write_smiles_with_params(&record, &params.smiles)
        .map(fixture_writer_text)
        .unwrap();

    assert_eq!(
        write_cx_smiles_with_params(&record, &params)
            .map(fixture_writer_text)
            .unwrap(),
        base
    );
}

#[test]
fn cx_writer_atom_property_escapes_each_period_with_exact_source_token() {
    let mut record = parse_smiles("C", &Default::default()).expect("parse carbon");
    record.topology.atoms[0]
        .set_prop("a.b", "c.d")
        .expect("set dotted atom property");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_atom_prop_params())
            .map(fixture_writer_text)
            .unwrap(),
        "C |atomProp:0.a&#46;b.c&#46;d|"
    );
    assert_eq!(record, before, "property escaping must preserve input");
}

#[test]
fn cx_writer_atom_property_preserves_empty_values() {
    let mut record = parse_smiles("C", &Default::default()).expect("parse carbon");
    record.topology.atoms[0]
        .set_prop("empty", "")
        .expect("set empty atom property");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_atom_prop_params())
            .map(fixture_writer_text)
            .unwrap(),
        "C |atomProp:0.empty.|"
    );
    assert_eq!(
        record, before,
        "empty property serialization preserves input"
    );
}

#[test]
fn cx_writer_atom_property_copies_non_ascii_text_without_normalization() {
    let mut record = parse_smiles("C", &Default::default()).expect("parse carbon");
    record.topology.atoms[0]
        .set_prop("clé", "λ雪")
        .expect("set non-ASCII atom property");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_atom_prop_params())
            .map(fixture_writer_text)
            .unwrap(),
        "C |atomProp:0.clé.λ雪|"
    );
    assert_eq!(
        record, before,
        "non-ASCII property serialization preserves input"
    );
}

#[test]
fn cx_writer_atom_properties_use_exact_pseudoatom_exclusions() {
    // Pinned get_atom_props_block skips dummyLabel only for attachment points,
    // the separate literal '*', and the exact {Pol, Mod} pseudoatom set.
    // These values were incorrectly suppressed by the former local 8-item set.
    let params = cx_atom_prop_params();
    for (label, expected) in [
        ("Het", "* |atomProp:0.dummyLabel.Het|"),
        ("Any", "* |atomProp:0.dummyLabel.Any|"),
        ("A", "* |atomProp:0.dummyLabel.A|"),
        ("Q", "* |atomProp:0.dummyLabel.Q|"),
        ("X", "* |atomProp:0.dummyLabel.X|"),
        ("Pol", "*"),
        ("Mod", "*"),
        ("*", "*"),
    ] {
        let mut record = parse_smiles("*", &Default::default()).expect("parse dummy atom");
        record.topology.atoms[0]
            .set_prop("dummyLabel", label)
            .expect("set dummy label");
        let before = record.clone();
        assert_eq!(
            write_cx_smiles_with_params(&record, &params)
                .map(fixture_writer_text)
                .unwrap(),
            expected,
            "dummyLabel={label:?}"
        );
        assert_eq!(
            record, before,
            "writing dummyLabel={label:?} preserves input"
        );
    }

    let mut attachment = parse_smiles("*", &Default::default()).expect("parse attachment");
    attachment.topology.atoms[0]
        .set_prop("dummyLabel", "Het")
        .expect("set dummy label");
    attachment.topology.atoms[0]
        .set_prop("_fromAttachPoint", "1")
        .expect("set attachment marker");
    let before = attachment.clone();
    assert_eq!(
        write_cx_smiles_with_params(&attachment, &params).unwrap(),
        "* |atomProp:0.dummyLabel.Het|".into(),
        "original unrecognized attachment key is ignored"
    );
    assert_eq!(attachment, before);
    attachment.topology.atoms[0]
        .set_prop("_fromAttchpt", "1")
        .expect("same original payload at source key");
    let source_before = attachment.clone();
    assert_eq!(
        write_cx_smiles_with_params(&attachment, &params)
            .map(fixture_writer_text)
            .unwrap(),
        "*",
        "attachment-point dummy labels stay excluded"
    );
    assert_eq!(attachment, source_before);
}

#[test]
fn cx_writer_atom_properties_exclude_skipped_private_and_computed_properties() {
    let mut record = parse_smiles("*", &Default::default()).expect("parse dummy atom");
    for (name, value) in [
        ("atomLabel", "label"),
        ("molFileValue", "value"),
        ("molParity", "1"),
        ("molAtomMapNumber", "7"),
        ("molStereoCare", "1"),
        ("molRxnExachg", "1"),
        ("molInversionFlag", "1"),
        ("_private", "private"),
    ] {
        record.topology.atoms[0]
            .set_prop(name, value)
            .unwrap_or_else(|error| panic!("set {name}: {error}"));
    }
    record.topology.atoms[0]
        .set_computed_prop("computed", "computed")
        .expect("set computed property");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_atom_prop_params())
            .map(fixture_writer_text)
            .unwrap(),
        "*"
    );
}

#[test]
fn cx_writer_atom_properties_preserve_source_insertion_order() {
    // Pinned Dict::setVal appends new keys and getPropList returns that order.
    // Keep the original executable counterexample and exact pinned output.
    let record = parse_smiles("C |atomProp:0.z.last:0.a.first|", &Default::default())
        .expect("parse ordered atom properties");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_atom_prop_params())
            .map(fixture_writer_text)
            .unwrap(),
        "C |atomProp:0.z.last:0.a.first|"
    );
    assert_eq!(record, before, "property writing preserves input");
}

#[test]
fn cx_writer_atom_properties_overwrite_without_moving_the_key() {
    let record = parse_smiles(
        "C |atomProp:0.z.first:0.a.middle:0.z.last|",
        &Default::default(),
    )
    .expect("parse duplicate ordered atom properties");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_atom_prop_params())
            .map(fixture_writer_text)
            .unwrap(),
        "C |atomProp:0.z.last:0.a.middle|"
    );
    assert_eq!(record, before, "property overwrite writing preserves input");
}

#[test]
fn cx_writer_atom_properties_clear_then_reinsert_appends_the_key() {
    let mut record = parse_smiles("C |atomProp:0.z.last:0.a.first|", &Default::default())
        .expect("parse ordered atom properties");
    record.topology.atoms[0].clear_prop("z");
    record.topology.atoms[0]
        .set_prop("z", "reinserted")
        .expect("reinsert atom property");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_atom_prop_params())
            .map(fixture_writer_text)
            .unwrap(),
        "C |atomProp:0.a.first:0.z.reinserted|"
    );
    assert_eq!(
        record, before,
        "property reinsertion writing preserves input"
    );
}

#[test]
fn cx_writer_atom_properties_filter_interleaved_keys_without_reordering() {
    let mut record = parse_smiles("C |atomProp:0.z.last|", &Default::default())
        .expect("parse initial atom property");
    let atom = &mut record.topology.atoms[0];
    atom.set_prop("_private", "hidden").expect("set private");
    atom.set_computed_prop("computed", "hidden")
        .expect("set computed");
    atom.set_prop("molParity", "1").expect("set skipped");
    atom.set_prop("a", "first").expect("set ordinary");
    atom.set_prop("molAtomMapNumber", "7")
        .expect("set second skipped");
    atom.set_prop("m", "middle").expect("set second ordinary");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_atom_prop_params())
            .map(fixture_writer_text)
            .unwrap(),
        "C |atomProp:0.z.last:0.a.first:0.m.middle|"
    );
    assert_eq!(record, before, "filtered property writing preserves input");
}

#[test]
fn cx_writer_atom_properties_project_typed_values_after_source_filters() {
    let mut record = parse_smiles("C", &Default::default()).expect("parse carbon");
    let atom = &mut record.topology.atoms[0];
    atom.set_prop("Text", PropertyValue::String("a.b".into()))
        .expect("set String property");
    atom.set_prop("Count", PropertyValue::Int(7))
        .expect("set Int property");
    atom.set_prop("Real", PropertyValue::Double(-0.0))
        .expect("set Double property");
    atom.set_prop("Active", PropertyValue::Bool(true))
        .expect("set Bool property");
    atom.set_prop("Count", PropertyValue::Int(7))
        .expect("overwrite Int property in place");
    atom.set_prop("atomLabel", PropertyValue::Int(9))
        .expect("set typed source-skipped property");
    atom.set_prop("_private", PropertyValue::Double(2.5))
        .expect("set typed private property");
    atom.set_computed_prop("computed", PropertyValue::Bool(false))
        .expect("set typed computed property");
    let before = record.clone();

    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_atom_prop_params())
            .map(fixture_writer_text)
            .unwrap(),
        "C |atomProp:0.Text.a&#46;b:0.Count.7:0.Real.-0:0.Active.1|"
    );
    assert_eq!(record, before, "typed property writing preserves input");
}

#[test]
fn cx_labels_and_values_use_dict_string_conversion_for_scalar_tags() {
    // RDKit 2026.03.1: MolToCXSmiles after SetIntProp(atomLabel, 9)
    // and SetBoolProp(molFileValue, true), independently checked at integration.
    let mut record = parse_smiles("C", &Default::default()).unwrap();
    record.topology.atoms[0]
        .set_prop("atomLabel", 9_i32)
        .unwrap();
    record.topology.atoms[0]
        .set_prop("molFileValue", true)
        .unwrap();
    let before = record.clone();
    assert_eq!(
        write_cx_smiles_with_params(&record, &Default::default())
            .map(fixture_writer_text)
            .unwrap(),
        "C |$9$,$_AV:1$|"
    );
    assert_eq!(record, before);
}

#[test]
fn cx_writer_atom_properties_survive_clone_and_topology_remap() {
    let source = parse_smiles(
        "C.N |atomProp:0.z.carbon:0.a.first:1.y.nitrogen:1.b.second|",
        &Default::default(),
    )
    .expect("parse disconnected ordered properties");
    let source_before = source.clone();
    let mut remapped = source.clone();
    let (topology, _) = remapped
        .topology
        .reordered_atoms(&[AtomId::new(1), AtomId::new(0)])
        .expect("reverse detached atom rows");
    remapped.topology = topology;
    let remapped_before = remapped.clone();

    assert_eq!(
        write_cx_smiles_with_params(&remapped, &cx_atom_prop_params())
            .map(fixture_writer_text)
            .unwrap(),
        "N.C |atomProp:0.y.nitrogen:0.b.second:1.z.carbon:1.a.first|"
    );
    assert_eq!(
        source, source_before,
        "topology remap preserves source input"
    );
    assert_eq!(
        remapped, remapped_before,
        "property writing preserves remapped input"
    );
}

#[test]
fn cpp_identifier_spellings_are_ordinary_property_keys_without_aliases() {
    // RDGeneral/types.h: the literal keys differ from these C++ identifiers.
    // Preserve the original mis-spelled inputs as explicit negative coverage.
    let mut record = parse_smiles("*", &Default::default()).unwrap();
    let atom = &mut record.topology.atoms[0];
    atom.set_prop("dummyLabel", "Het").unwrap();
    atom.set_prop("atomLabel", "ordinary").unwrap();
    atom.set_prop("_fromAttachPoint", "1").unwrap();
    atom.set_prop("molRxnExactChange", "1").unwrap();
    record
        .properties
        .set_prop("_MolFileLinkNodes", "1 3 2 1 2 1 3")
        .unwrap();
    let before = record.clone();
    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_atom_label_params())
            .map(fixture_writer_text)
            .unwrap(),
        "* |$ordinary$|"
    );
    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_atom_prop_params())
            .map(fixture_writer_text)
            .unwrap(),
        "* |atomProp:0.dummyLabel.Het:0.molRxnExactChange.1|"
    );
    assert_eq!(
        write_cx_smiles_with_params(&record, &cx_params(CxSmilesFields::LINKNODES, true, false))
            .map(fixture_writer_text)
            .unwrap(),
        "*"
    );
    assert_eq!(record, before);
}
