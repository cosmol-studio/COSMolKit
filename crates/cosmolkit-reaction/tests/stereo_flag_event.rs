use cosmolkit_core::{ValenceModel, assign_valence_with_options_for_topology, fast_find_rings};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, CoordinateBlock, MoleculeProperties,
    PropertyValue, TopologyBlock,
};
use cosmolkit_reaction::{
    Reaction, ReactionInput, ReactionRunParams, ReactionValidationParams, initialize_reaction,
    parse_smirks, run_reactants,
};
use cosmolkit_types::{BondDirection, BondOrder, Element};

#[test]
fn actual_reaction_product_retains_executed_source_stereo_flag() {
    // ReactionRunner.cpp generateOneProductSet dispatches direction generation
    // iff any input SINGLE has direction>NONE. Chirality.cpp setStereoForBond
    // writes ordinary Int1 after the successful squiggle stereo write.
    // This exercises returned ReactionProduct.properties, not a second store.
    let reaction = parse_smirks("[C:1]-[C:2]=[C:3]-[C:4]>>[C:1]-[C:2]=[C:3]-[C:4]").unwrap();
    // Source fully mapped SINGLE bonds keep product template direction; they
    // do not inherit the reactant's direction. Vary both native guard inputs.
    for (direction, template_direction, executed_write) in [
        (BondDirection::None, BondDirection::None, false),
        (BondDirection::None, BondDirection::Unknown, false),
        (BondDirection::Unknown, BondDirection::None, false),
        (BondDirection::Unknown, BondDirection::Unknown, true),
    ] {
        let mut products = reaction.product_templates().to_vec();
        products[0].bonds_mut()[0]
            .bond_mut()
            .set_direction(template_direction);
        let current =
            Reaction::from_templates(reaction.reactant_templates().to_vec(), products, vec![])
                .unwrap();
        let mut current =
            initialize_reaction(&current, &ReactionValidationParams::default()).unwrap();
        let topology = TopologyBlock::try_from_parts(
            (0..4)
                .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            vec![
                Bond::from_spec(
                    BondId::new(0),
                    BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                        .with_direction(direction),
                ),
                Bond::from_spec(
                    BondId::new(1),
                    BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Double),
                ),
                Bond::from_spec(
                    BondId::new(2),
                    BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Single),
                ),
            ],
            vec![],
            vec![],
        )
        .unwrap();
        let rings = fast_find_rings(&topology).unwrap();
        let valence =
            assign_valence_with_options_for_topology(&topology, ValenceModel::RdkitLike, false)
                .unwrap();
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let outputs = run_reactants(
            &mut current,
            &[ReactionInput {
                topology: &topology,
                coordinates: &coordinates,
                properties: &properties,
                rings: Some(&rings),
                valence: Some(&valence),
            }],
            &ReactionRunParams {
                max_products: 1,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(outputs.len(), 1);
        assert_eq!(outputs[0].len(), 1);
        let product = &outputs[0][0];
        if executed_write {
            assert_eq!(
                product.properties.prop("_needsDetectBondStereo"),
                Some(&PropertyValue::Int(1))
            );
        } else {
            assert!(product.properties.prop("_needsDetectBondStereo").is_none());
        }
        assert!(properties.prop("_needsDetectBondStereo").is_none());
    }
}
