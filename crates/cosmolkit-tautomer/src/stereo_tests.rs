use super::*;
pub(crate) fn fixture_from_smiles(
    text: &str,
) -> Result<TautomerRecord, Box<dyn std::error::Error>> {
    let params = cosmolkit_smiles::SmilesParseParams::default();
    let parsed = cosmolkit_smiles::parse_smiles(text, &params)?;
    let result = remove_hydrogens_with_params(
        parsed.topology,
        parsed.coordinates,
        parsed.properties,
        &RemoveHsParams {
            update_explicit_count: true,
            sanitize: true,
            ..Default::default()
        },
    )?;
    let mut valence = result.final_valence;
    let mut rings = result.final_rings;
    let record = cosmolkit_smiles::finalize_smiles_stereo(
        cosmolkit_smiles::SmilesRecord {
            topology: result.topology,
            coordinates: result.coordinates,
            properties: result.properties,
        },
        &params,
        &mut valence,
        &mut rings,
    )?;
    Ok(prepared(TautomerRecordView {
        topology: &record.topology,
        coordinates: &record.coordinates,
        properties: &record.properties,
        valence: valence.as_ref(),
        rings: rings.as_ref(),
    })?)
}
fn string_prop<'a>(atom: &'a cosmolkit_model::Atom, key: &str) -> Option<&'a str> {
    match atom.prop(key) {
        Some(cosmolkit_model::PropertyValue::String(value)) => {
            Some(std::str::from_utf8(value.as_bytes()).expect("original fixed ASCII CIP fixture"))
        }
        _ => None,
    }
}
fn set_tautomer_stereo_and_isotopic_hydrogens(
    source: &TautomerRecord,
    tautomer: &mut TautomerRecord,
    atoms: &BTreeSet<AtomId>,
    bonds: &BTreeSet<BondId>,
    params: TautomerParams,
) -> Result<bool, TautomerRunError> {
    super::set_tautomer_stereo_and_isotopic_hydrogens(
        source,
        tautomer,
        atoms,
        bonds,
        params,
        &CoordinateBlock::default(),
    )
}
fn marked_atoms(indices: impl IntoIterator<Item = usize>) -> BTreeSet<AtomId> {
    indices.into_iter().map(AtomId::new).collect()
}

fn marked_bonds(indices: impl IntoIterator<Item = usize>) -> BTreeSet<BondId> {
    indices.into_iter().map(BondId::new).collect()
}

#[test]
fn stereo_and_isotopic_hydrogens_sp2_and_remove_sp3_clear_chiral_and_cip_state() {
    let mut source = fixture_from_smiles("[C@H](F)(Cl)Br").expect("chiral source");
    source.topology.atoms[0].set_prop("_CIPCode", "R").unwrap();
    let source_before = source.clone();

    for (hybridization, remove_sp3_stereo) in
        [(Hybridization::Sp2, false), (Hybridization::Sp3, true)]
    {
        let mut tautomer = source.clone();
        tautomer.topology.atoms[0].set_hybridization(hybridization);
        let options = TautomerParams::default()
            .with_remove_sp3_stereo(remove_sp3_stereo)
            .with_reassign_stereo(false);
        let changed = set_tautomer_stereo_and_isotopic_hydrogens(
            &source,
            &mut tautomer,
            &marked_atoms([0]),
            &BTreeSet::new(),
            options,
        )
        .expect("apply atom stereo transition");

        assert!(changed);
        assert_eq!(
            tautomer.topology.atoms.as_slice()[0].chiral_tag(),
            ChiralTag::Unspecified
        );
        assert_eq!(string_prop(&tautomer.topology.atoms[0], "_CIPCode"), None);
        assert_eq!(
            tautomer.properties.prop("_StereochemDone"),
            Some(&cosmolkit_model::PropertyValue::Int(1))
        );
        assert_eq!(source, source_before);
    }
}

#[test]
fn stereo_and_isotopic_hydrogens_sp3_restore_copies_source_tag_and_present_cip_only() {
    let mut source = fixture_from_smiles("[C@H](F)(Cl)Br").expect("chiral source");
    source.topology.atoms[0].set_prop("_CIPCode", "S").unwrap();
    let mut tautomer = source.clone();
    let atom = &mut tautomer.topology.atoms[0];
    atom.set_hybridization(Hybridization::Sp3);
    atom.set_chiral_tag(match source.topology.atoms.as_slice()[0].chiral_tag() {
        ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
        _ => ChiralTag::TetrahedralCw,
    });
    atom.set_prop("_CIPCode", "R").unwrap();

    let changed = set_tautomer_stereo_and_isotopic_hydrogens(
        &source,
        &mut tautomer,
        &marked_atoms([0]),
        &BTreeSet::new(),
        TautomerParams::default()
            .with_remove_sp3_stereo(false)
            .with_reassign_stereo(false),
    )
    .expect("restore source stereo");

    assert!(changed);
    assert_eq!(
        tautomer.topology.atoms.as_slice()[0].chiral_tag(),
        source.topology.atoms.as_slice()[0].chiral_tag()
    );
    assert_eq!(
        string_prop(&tautomer.topology.atoms[0], "_CIPCode"),
        Some("S")
    );

    source.topology.atoms[0].clear_prop("_CIPCode");
    tautomer.topology.atoms[0]
        .set_prop("_CIPCode", "R")
        .unwrap();
    set_tautomer_stereo_and_isotopic_hydrogens(
        &source,
        &mut tautomer,
        &marked_atoms([0]),
        &BTreeSet::new(),
        TautomerParams::default()
            .with_remove_sp3_stereo(false)
            .with_reassign_stereo(false),
    )
    .expect("source-absent CIP branch");
    assert_eq!(
        string_prop(&tautomer.topology.atoms[0], "_CIPCode"),
        Some("R")
    );
}

#[test]
fn stereo_and_isotopic_hydrogens_remove_or_zero_total_h_clear_tracked_isotopes() {
    let atom = AtomId::new(0);
    let topology = cosmolkit_model::TopologyBlock::try_from_parts(
        vec![cosmolkit_model::Atom::from_spec(
            atom,
            cosmolkit_model::AtomSpec::new(cosmolkit_model::Element::C)
                .with_no_implicit(true)
                .with_explicit_hydrogens(1)
                .with_tracked_isotopic_hydrogens(vec![2, 3]),
        )],
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    let source = prepared(TautomerRecordView {
        topology: &topology,
        coordinates: &CoordinateBlock::default(),
        properties: &MoleculeProperties::default(),
        valence: None,
        rings: None,
    })
    .expect("isotopic-H source");

    let mut retained = source.clone();
    set_tautomer_stereo_and_isotopic_hydrogens(
        &source,
        &mut retained,
        &marked_atoms([atom.index()]),
        &BTreeSet::new(),
        TautomerParams::default()
            .with_remove_isotopic_hydrogens(false)
            .with_reassign_stereo(false),
    )
    .expect("retain tracked isotopes");
    assert_eq!(
        retained.topology.atoms.as_slice()[0].tracked_isotopic_hydrogens(),
        &[2, 3]
    );

    let mut removed_by_option = source.clone();
    set_tautomer_stereo_and_isotopic_hydrogens(
        &source,
        &mut removed_by_option,
        &marked_atoms([0]),
        &BTreeSet::new(),
        TautomerParams::default().with_reassign_stereo(false),
    )
    .expect("remove tracked isotopes by option");
    assert!(
        removed_by_option.topology.atoms.as_slice()[0]
            .tracked_isotopic_hydrogens()
            .is_empty()
    );

    let mut removed_at_zero_h = source.clone();
    removed_at_zero_h.topology.atoms[0].set_explicit_hydrogens(0);
    set_tautomer_stereo_and_isotopic_hydrogens(
        &source,
        &mut removed_at_zero_h,
        &marked_atoms([0]),
        &BTreeSet::new(),
        TautomerParams::default()
            .with_remove_isotopic_hydrogens(false)
            .with_reassign_stereo(false),
    )
    .expect("remove tracked isotopes at zero total H");
    assert!(
        removed_at_zero_h.topology.atoms.as_slice()[0]
            .tracked_isotopic_hydrogens()
            .is_empty()
    );
}

#[test]
fn stereo_and_isotopic_hydrogens_non_ring_double_bond_removal_sets_any_and_clears_dirs() {
    let source = fixture_from_smiles("F/C=C/Cl").expect("E/Z source");
    let source_before = source.clone();
    let double_bond = source
        .topology
        .bonds
        .as_slice()
        .iter()
        .find(|bond| bond.order() == BondOrder::Double)
        .expect("double bond")
        .id();
    assert!(is_stereo_beyond_any(
        source.topology.bonds.as_slice()[double_bond.index()].stereo()
    ));
    let directional = source
        .topology
        .bonds
        .as_slice()
        .iter()
        .filter(|bond| {
            matches!(
                bond.direction(),
                BondDirection::EndDownRight | BondDirection::EndUpRight
            )
        })
        .map(cosmolkit_model::Bond::id)
        .collect::<Vec<_>>();
    assert!(!directional.is_empty());
    let mut tautomer = source.clone();

    set_tautomer_stereo_and_isotopic_hydrogens(
        &source,
        &mut tautomer,
        &BTreeSet::new(),
        &marked_bonds([double_bond.index()]),
        TautomerParams::default().with_reassign_stereo(false),
    )
    .expect("remove non-ring bond stereo");

    assert_eq!(
        tautomer.topology.bonds.as_slice()[double_bond.index()].stereo(),
        BondStereo::Any
    );
    assert_eq!(
        tautomer.topology.bonds.as_slice()[double_bond.index()].stereo_atoms(),
        None
    );
    for bond_id in directional {
        assert_eq!(
            tautomer.topology.bonds.as_slice()[bond_id.index()].direction(),
            BondDirection::None
        );
    }
    assert_eq!(source, source_before);
}

#[test]
fn stereo_and_isotopic_hydrogens_preservation_restores_bond_stereo_atoms_and_dirs() {
    let source = fixture_from_smiles("F/C=C/Cl").expect("E/Z source");
    let double_bond = source
        .topology
        .bonds
        .as_slice()
        .iter()
        .find(|bond| bond.order() == BondOrder::Double)
        .expect("double bond")
        .id();
    let mut tautomer = source.clone();
    tautomer.topology.bonds[double_bond.index()].set_stereo_atoms(None);
    tautomer.topology.bonds[double_bond.index()]
        .set_stereo(BondStereo::Any)
        .unwrap();
    for bond in &mut tautomer.topology.bonds {
        if matches!(
            bond.direction(),
            BondDirection::EndDownRight | BondDirection::EndUpRight
        ) {
            bond.set_direction(BondDirection::None);
        }
    }

    let changed = set_tautomer_stereo_and_isotopic_hydrogens(
        &source,
        &mut tautomer,
        &BTreeSet::new(),
        &marked_bonds([double_bond.index()]),
        TautomerParams::default()
            .with_remove_bond_stereo(false)
            .with_reassign_stereo(false),
    )
    .expect("restore source bond stereo");

    assert!(changed);
    assert_eq!(
        tautomer.topology.bonds.as_slice()[double_bond.index()].stereo(),
        source.topology.bonds.as_slice()[double_bond.index()].stereo()
    );
    assert_eq!(
        tautomer.topology.bonds.as_slice()[double_bond.index()].stereo_atoms(),
        source.topology.bonds.as_slice()[double_bond.index()].stereo_atoms()
    );
    for (actual, expected) in tautomer
        .topology
        .bonds
        .as_slice()
        .iter()
        .zip(source.topology.bonds.as_slice())
    {
        assert_eq!(actual.direction(), expected.direction());
    }
}

#[test]
fn stereo_and_isotopic_hydrogens_ring_and_non_double_bonds_use_none() {
    let ring_source = fixture_from_smiles("C1=CCCCC1").expect("ring alkene");
    let ring_bond = ring_source
        .topology
        .bonds
        .as_slice()
        .iter()
        .find(|bond| bond.order() == BondOrder::Double)
        .expect("ring double bond")
        .id();
    let mut ring_tautomer = ring_source.clone();
    ring_tautomer.topology.bonds[ring_bond.index()]
        .set_stereo(BondStereo::Any)
        .unwrap();
    set_tautomer_stereo_and_isotopic_hydrogens(
        &ring_source,
        &mut ring_tautomer,
        &BTreeSet::new(),
        &marked_bonds([ring_bond.index()]),
        TautomerParams::default().with_reassign_stereo(false),
    )
    .expect("ring fallback");
    assert_eq!(
        ring_tautomer.topology.bonds.as_slice()[ring_bond.index()].stereo(),
        BondStereo::None
    );

    let source = fixture_from_smiles("CC").expect("single bond");
    let mut tautomer = source.clone();
    tautomer.topology.bonds[0]
        .set_stereo(BondStereo::Any)
        .unwrap();
    set_tautomer_stereo_and_isotopic_hydrogens(
        &source,
        &mut tautomer,
        &BTreeSet::new(),
        &marked_bonds([0]),
        TautomerParams::default().with_reassign_stereo(false),
    )
    .expect("single bond cleanup");
    assert_eq!(
        tautomer.topology.bonds.as_slice()[0].stereo(),
        BondStereo::None
    );
    assert_eq!(tautomer.topology.bonds.as_slice()[0].stereo_atoms(), None);
}

#[test]
fn stereo_and_isotopic_hydrogens_reassignment_reapplies_any_contract() {
    let source = fixture_from_smiles("F/C=C/Cl").expect("E/Z source");
    let double_bond = source
        .topology
        .bonds
        .as_slice()
        .iter()
        .find(|bond| bond.order() == BondOrder::Double)
        .expect("double bond")
        .id();
    let mut tautomer = source.clone();

    set_tautomer_stereo_and_isotopic_hydrogens(
        &source,
        &mut tautomer,
        &BTreeSet::new(),
        &marked_bonds([double_bond.index()]),
        TautomerParams::default(),
    )
    .expect("reassign and reapply explicit undefined stereo");

    assert_eq!(
        tautomer.topology.bonds.as_slice()[double_bond.index()].stereo(),
        BondStereo::Any
    );
    assert_eq!(
        tautomer.topology.bonds.as_slice()[double_bond.index()].stereo_atoms(),
        None
    );
    assert_eq!(
        tautomer.properties.prop("_StereochemDone"),
        Some(&cosmolkit_model::PropertyValue::Int(1))
    );
    assert!(
        tautomer
            .properties
            .is_prop_computed("_StereochemDone")
            .unwrap()
    );
}

#[test]
fn stereo_and_isotopic_hydrogens_executes_every_option_combination_without_source_mutation() {
    let source = fixture_from_smiles("F/C=C/[C@H](Cl)Br").expect("combined source");
    let source_before = source.clone();
    let double_bond = source
        .topology
        .bonds
        .as_slice()
        .iter()
        .find(|bond| bond.order() == BondOrder::Double)
        .expect("double bond")
        .id();
    let chiral_atom = source
        .topology
        .atoms
        .as_slice()
        .iter()
        .find(|atom| atom.chiral_tag() != ChiralTag::Unspecified)
        .expect("chiral atom")
        .id();

    for bits in 0_u8..16 {
        let options = TautomerParams::default()
            .with_remove_sp3_stereo(bits & 1 != 0)
            .with_remove_bond_stereo(bits & 2 != 0)
            .with_remove_isotopic_hydrogens(bits & 4 != 0)
            .with_reassign_stereo(bits & 8 != 0);
        let mut tautomer = source.clone();
        set_tautomer_stereo_and_isotopic_hydrogens(
            &source,
            &mut tautomer,
            &marked_atoms([chiral_atom.index()]),
            &marked_bonds([double_bond.index()]),
            options,
        )
        .unwrap_or_else(|error| panic!("option mask {bits:#06b}: {error}"));
        assert_eq!(source, source_before, "option mask {bits:#06b}");
        if options.remove_bond_stereo() {
            assert_eq!(
                tautomer.topology.bonds.as_slice()[double_bond.index()].stereo(),
                BondStereo::Any,
                "option mask {bits:#06b}"
            );
        }
        if !options.reassign_stereo() {
            assert_eq!(
                tautomer.properties.prop("_StereochemDone"),
                Some(&cosmolkit_model::PropertyValue::Int(1))
            );
            assert_eq!(
                tautomer
                    .properties
                    .is_prop_computed("_StereochemDone")
                    .unwrap(),
                source
                    .properties
                    .is_prop_computed("_StereochemDone")
                    .unwrap()
            );
        }
    }
}

#[test]
fn source_without_stereo_atom_pair_retains_candidate_pair() {
    let mut source = fixture_from_smiles("F/C=C/Cl").unwrap();
    let bond = source
        .topology
        .bonds
        .iter()
        .find(|b| b.order() == BondOrder::Double)
        .unwrap()
        .id();
    let mut candidate = source.clone();
    let pair = candidate.topology.bonds[bond.index()].stereo_atoms();
    assert!(pair.is_some());
    source.topology.bonds[bond.index()].set_stereo_atoms(None);
    source.topology.bonds[bond.index()]
        .set_stereo(BondStereo::None)
        .unwrap();
    let before = source.clone();
    let changed = set_tautomer_stereo_and_isotopic_hydrogens(
        &source,
        &mut candidate,
        &BTreeSet::new(),
        &marked_bonds([bond.index()]),
        TautomerParams::default()
            .with_remove_bond_stereo(false)
            .with_reassign_stereo(false),
    )
    .unwrap();
    assert!(changed);
    assert_eq!(candidate.topology.bonds[bond.index()].stereo_atoms(), pair);
    assert_eq!(
        candidate.topology.bonds[bond.index()].stereo(),
        BondStereo::None
    );
    assert_eq!(source, before);
}
#[test]
fn cip_restoration_uses_source_string_projection() {
    let mut source = fixture_from_smiles("[C@H](F)(Cl)Br").unwrap();
    source.topology.atoms[0]
        .set_prop("_CIPCode", 17_i32)
        .unwrap();
    let mut candidate = source.clone();
    candidate.topology.atoms[0].set_hybridization(Hybridization::Sp3);
    let before = source.clone();
    set_tautomer_stereo_and_isotopic_hydrogens(
        &source,
        &mut candidate,
        &marked_atoms([0]),
        &BTreeSet::new(),
        TautomerParams::default()
            .with_remove_sp3_stereo(false)
            .with_reassign_stereo(false),
    )
    .unwrap();
    assert_eq!(
        string_prop(&candidate.topology.atoms[0], "_CIPCode"),
        Some("17")
    );
    assert_eq!(source, before);
}
#[test]
fn ring_detection_retains_source_fast_ring_update_without_stereo_reassignment() {
    let source = fixture_from_smiles("C1=CCCCC1").unwrap();
    let mut candidate = source.clone();
    candidate.rings = RingInfo::new(
        RingFindType::OtherOrUnknown,
        candidate.topology.atoms.len(),
        candidate.topology.bonds.len(),
    );
    let before = source.clone();
    let bond = source
        .topology
        .bonds
        .iter()
        .find(|b| b.order() == BondOrder::Double)
        .unwrap()
        .id();
    set_tautomer_stereo_and_isotopic_hydrogens(
        &source,
        &mut candidate,
        &BTreeSet::new(),
        &marked_bonds([bond.index()]),
        TautomerParams::default().with_reassign_stereo(false),
    )
    .unwrap();
    assert!(candidate.rings.is_find_fast_or_better());
    assert!(!candidate.rings.is_symm_sssr());
    assert_eq!(
        candidate.topology.bonds[bond.index()].stereo(),
        BondStereo::None
    );
    assert_eq!(source, before);
}

#[test]
fn single_query_endpoint_transfers_hydrogen_in_source_order() {
    let source = fixture_from_smiles("[CH4]").unwrap();
    let candidate = kekulized(&source).unwrap();
    let catalog = crate::TautomerCatalog::from_data(&[("same endpoint", "[CH4]", "", "")]).unwrap();
    let transform = &catalog.transforms()[0];
    let matches = transform_matches(&candidate, &CoordinateBlock::default(), transform).unwrap();
    assert_eq!(matches.len(), 1);
    let before = source.clone();
    let attempt = apply_tautomer_transform_match(
        &source,
        &candidate,
        &CoordinateBlock::default(),
        transform,
        &matches[0],
        &BTreeSet::new(),
        &BTreeSet::new(),
        &|_| false,
        TautomerParams::default(),
    )
    .unwrap();
    let TautomerExpansionAttempt::Product(product) = attempt else {
        panic!("same endpoint must produce a retained candidate: {attempt:?}")
    };
    assert_eq!(product.tautomer.topology.atoms[0].explicit_hydrogens(), 4);
    assert_eq!(product.canonical_smiles.as_bytes(), b"C");
    assert_eq!(source, before);
}

#[test]
fn source_property_failure_cip_guard_retains_prefix_and_absent_branch() {
    let source = fixture_from_smiles("[C@H](F)(Cl)Br").unwrap();
    let source_before = source.clone();
    for cip_present in [false, true] {
        let mut tautomer = source.clone();
        tautomer.topology.atoms[0].clear_prop("_CIPCode").unwrap();
        if cip_present {
            tautomer.topology.atoms[0]
                .set_prop("_CIPCode", "R")
                .unwrap();
        }
        tautomer.topology.atoms[0]
            .set_prop("__computedProps", cosmolkit_model::PropertyValue::Int(7))
            .unwrap();
        let before = tautomer.clone();
        let options = TautomerParams::default().with_reassign_stereo(false);
        let result = set_tautomer_stereo_and_isotopic_hydrogens(
            &source,
            &mut tautomer,
            &marked_atoms([0]),
            &BTreeSet::new(),
            options,
        );
        if cip_present {
            assert!(matches!(
                result,
                Err(TautomerRunError::AtomProperty(
                    cosmolkit_model::AtomPropertyError::ComputedListKind(_)
                ))
            ));
            let mut expected_prefix = before;
            expected_prefix.topology.atoms[0].set_chiral_tag(ChiralTag::Unspecified);
            assert_eq!(tautomer, expected_prefix);
        } else {
            assert!(result.unwrap());
            assert_eq!(tautomer.topology.atoms[0].prop("_CIPCode"), None);
            assert_eq!(
                tautomer.topology.atoms[0].chiral_tag(),
                ChiralTag::Unspecified
            );
            assert_eq!(
                tautomer.topology.atoms[0].prop("__computedProps"),
                Some(&cosmolkit_model::PropertyValue::Int(7))
            );
            assert_eq!(
                tautomer.properties.prop("_StereochemDone"),
                Some(&cosmolkit_model::PropertyValue::Int(1))
            );
        }
        assert_eq!(source, source_before);
    }
}
