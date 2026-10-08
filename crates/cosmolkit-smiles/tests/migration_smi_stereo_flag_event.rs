use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, CoordinateBlock, MoleculeProperties,
    PropertyValue, TopologyBlock,
};
use cosmolkit_smiles::{
    SmilesParseParams, SmilesRecord, SmilesStereoError, finalize_smiles_stereo,
};
use cosmolkit_types::{BondDirection, BondOrder, Element};

fn record() -> SmilesRecord {
    let atoms = (0..4)
        .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
        .collect();
    let bonds = vec![
        Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                .with_direction(BondDirection::Unknown),
        ),
        Bond::from_spec(
            BondId::new(1),
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Double),
        ),
        Bond::from_spec(
            BondId::new(2),
            BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Single),
        ),
    ];
    let mut properties = MoleculeProperties::default();
    properties.set_prop("head", "kept").unwrap();
    properties
        .set_computed_prop("_needsDetectBondStereo", 0_i32)
        .unwrap();
    properties.set_prop("tail", "kept").unwrap();
    SmilesRecord {
        topology: TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap(),
        coordinates: CoordinateBlock::default(),
        properties,
    }
}

#[test]
fn source_flag_effect_is_cleared_after_successful_direction_processing() {
    // Source SmilesParse.cpp direction reconstruction precedes clearProp;
    // ordinary setProp retains prior computed membership until that clear.
    let result = finalize_smiles_stereo(
        record(),
        &SmilesParseParams {
            sanitize: true,
            remove_hydrogens: false,
            ..Default::default()
        },
        &mut None,
        &mut None,
    )
    .unwrap();
    assert!(result.properties.prop("_needsDetectBondStereo").is_none());
    assert_eq!(
        result.properties.prop("__computedProps"),
        Some(&PropertyValue::StringVector(vec!["_StereochemDone".into()]))
    );
    assert_eq!(
        result.properties.prop("head"),
        Some(&PropertyValue::from("kept"))
    );
    assert_eq!(
        result.properties.prop("tail"),
        Some(&PropertyValue::from("kept"))
    );
}

#[test]
fn source_computed_list_cast_failure_is_propagated_by_actual_smiles_finalizer() {
    // RDProps setProp(computed=false) does not cast __computedProps; clearProp
    // always casts it. findSSSR clears extraRings first, before direction
    // processing or later stereo assignment. Preserve the original cast cause.
    let mut input = record();
    input.properties.set_prop("__computedProps", 3_i32).unwrap();
    let before = input.clone();
    let result = finalize_smiles_stereo(
        input,
        &SmilesParseParams {
            sanitize: true,
            remove_hydrogens: false,
            ..Default::default()
        },
        &mut None,
        &mut None,
    );
    assert!(matches!(
        result,
        Err(SmilesStereoError::Rings(
            cosmolkit_core::RingFindingError::MoleculeProperty(
                cosmolkit_model::MoleculePropertyError::ComputedListKind(_)
            )
        ))
    ));
    assert_eq!(
        before.properties.prop("_needsDetectBondStereo"),
        Some(&PropertyValue::Int(0))
    );
    assert_eq!(
        before.properties.prop("__computedProps"),
        Some(&PropertyValue::Int(3))
    );
}
