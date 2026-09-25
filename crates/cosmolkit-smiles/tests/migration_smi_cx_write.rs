use cosmolkit_core::KekulizeError;
use cosmolkit_model::{
    AtomId, CoordinateBlock, MoleculeProperties, StereoGroup, StereoGroupKind, TopologyBlock,
};
use cosmolkit_smiles::{
    CxSmilesFields, CxSmilesWriteParams, SmilesParseError, SmilesRecord, SmilesWriteParams,
    parse_smiles, write_cx_smiles_with_params,
};
use cosmolkit_types::{BondDirection, ChiralTag};

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
    }
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
        write_cx_smiles_with_params(&record, &kekule_cx_params(CxSmilesFields::NONE)).unwrap(),
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
        write_cx_smiles_with_params(&record, &kekule_cx_params(CxSmilesFields::ALL)).unwrap(),
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
        write_cx_smiles_with_params(&record, &cx_params(CxSmilesFields::ALL, true, true)).unwrap(),
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
        write_cx_smiles_with_params(&record, &cx_params(CxSmilesFields::ALL, true, false)).unwrap(),
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
            .unwrap(),
        "CC(O)Cl"
    );
    assert_eq!(record, before, "CX preprocessing must preserve its caller");
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
        };

        assert_eq!(
            write_cx_smiles_with_params(&record, &params).unwrap(),
            expected,
            "cleanStereo={clean_stereo}, done marker present={done_marker}"
        );
        assert_eq!(record, before, "CX writing must preserve its caller");
    }
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
    };

    assert_eq!(
        write_cx_smiles_with_params(&record, &params).unwrap(),
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
    assert_eq!(record.properties.prop("_needsDetectBondStereo"), Some("1"));
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
    };

    assert_eq!(
        write_cx_smiles_with_params(&record, &params).unwrap(),
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
        write_cx_smiles_with_params(&record, &params).unwrap(),
        "CO.CO |(1,0,;0,0,;2,0,;3,0,)|"
    );
    assert_eq!(
        record, before,
        "coordinate-map emission must preserve its caller"
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
        .expect_err("source unsigned-int property read must fail visibly");
    assert!(matches!(
        error,
        SmilesParseError::WriterStereo(message) if message.contains("bad_any_cast")
    ));
    assert_eq!(record, before, "failed CX writing must preserve its caller");
}
