#![cfg(all(feature = "cap-reaction", feature = "cap-smiles"))]
use cosmolkit::*;

#[test]
#[cfg(feature = "cap-sanitize")]
fn sanitized_reaction_products_retain_source_ring_state_for_nitrogen_stereo() {
    // RDKit 2026.03.1, ReactionFromSmarts -> RunReactants -> SanitizeMol
    // -> MolToSmiles. Raw and sanitized products retain the same atom tag;
    // the existing SymmSSSR cache controls nitrogen's legacy stereo legality.
    let source = Molecule::from_smiles("C1C[C@H]2CC[C@H]2C1").unwrap();
    let before = source.clone();
    let mut reaction = Reaction::from_smirks("[C;H1:1]>>[N:1]").unwrap();
    let products = source
        .reaction_products_from_inputs(&mut reaction, &[&source], &ReactionRunParams::default())
        .unwrap();
    assert_eq!(products.len(), 2);
    for (group, expected) in products.iter().zip(["C1C[C@@H]2CCN2C1", "C1C[C@H]2CCN2C1"]) {
        assert_eq!(group.len(), 1);
        let raw = &group[0];
        assert_eq!(raw.to_smiles().unwrap().as_bytes(), expected.as_bytes());
        let raw_before = raw.clone();
        let sanitized = raw.sanitize().unwrap();
        let sanitized_before = sanitized.clone();
        assert_eq!(
            sanitized.to_smiles().unwrap().as_bytes(),
            expected.as_bytes()
        );
        assert_eq!(
            sanitized
                .without_hydrogens()
                .unwrap()
                .to_smiles()
                .unwrap()
                .as_bytes(),
            expected.as_bytes()
        );
        assert_eq!(raw.topology(), raw_before.topology());
        assert_eq!(raw.properties(), raw_before.properties());
        assert_eq!(sanitized.topology(), sanitized_before.topology());
        assert_eq!(sanitized.properties(), sanitized_before.properties());
    }
    assert_eq!(source.topology(), before.topology());
    assert_eq!(source.properties(), before.properties());
}

#[test]
fn reaction_writer_clears_unpaired_template_bond_directions() {
    // Pinned RDKit ReactionToSmarts delegates each role to MolToSmarts.
    for (input, expected) in [
        ("[C:1]>>[C:1]/C=C", "[C:1]>>[C:1]C=C"),
        (r"[C:1]>>[C:1]\C=C", "[C:1]>>[C:1]C=C"),
        (
            "[C:1]/[C:2]=[C:3]>>[C:1]/[C:2]=[C:3]",
            "[C:1][C:2]=[C:3]>>[C:1][C:2]=[C:3]",
        ),
        (
            "[C:1]/[C:2]=[C:3]>>[C:1][C:2]=[C:3]",
            "[C:1][C:2]=[C:3]>>[C:1][C:2]=[C:3]",
        ),
    ] {
        let reaction = Reaction::from_smirks(input).unwrap();
        let reactant_before = reaction.reactant_template(0).unwrap().clone();
        let product_before = reaction.product_template(0).unwrap().clone();
        let nonisomeric = ReactionWriteParams {
            isomeric_smiles: false,
            ..Default::default()
        };
        for output in [
            reaction.to_smirks(),
            reaction.to_smirks_with_params(&nonisomeric),
            reaction.to_cx_smirks(),
            reaction.to_cx_smirks_with_params(&nonisomeric),
        ] {
            assert_eq!(output.unwrap().as_bytes(), expected.as_bytes(), "{input}");
        }
        assert_eq!(reaction.reactant_template(0).unwrap(), &reactant_before);
        assert_eq!(reaction.product_template(0).unwrap(), &product_before);
    }
}

#[test]
fn reaction_values_preserve_templates_settings_and_initialization() {
    let rxn = Reaction::from_smirks("[C:1]>O>[N:1]").unwrap();
    assert_eq!(
        (
            rxn.num_reactant_templates(),
            rxn.num_agent_templates(),
            rxn.num_product_templates()
        ),
        (1, 1, 1)
    );
    assert!(!rxn.is_initialized());
    assert!(
        rxn.validate_with_params(&ReactionValidationParams::new(true))
            .unwrap()
            .is_valid()
    );
    assert!(rxn.with_initialized().unwrap().is_initialized());
    assert!(!rxn.is_initialized());
    assert_eq!(rxn.without_agents().removed_templates().len(), 1);
    assert_eq!(rxn.without_agents().reaction().num_agent_templates(), 0);
    assert_eq!(rxn.num_agent_templates(), 1);
    let text = rxn.to_smirks().unwrap();
    let roundtrip = parse_smirks(std::str::from_utf8(text.as_bytes()).unwrap()).unwrap();
    assert_eq!(roundtrip.num_agent_templates(), 1);
    assert!(matches!(
        rxn.reactant_template(9),
        Err(ReactionModelError::TemplateIndex {
            index: 9,
            count: 1,
            ..
        })
    ));
    assert!(matches!(
        Reaction::from_smirks("CC"),
        Err(ReactionParseError::Separators { count: 0 })
    ));
}

#[test]
fn product_groups_initialize_same_reaction_and_keep_sources() {
    let source = Molecule::from_smiles("C").unwrap();
    let peer = source.clone();
    let mut rxn = Reaction::from_smirks("[C:1]>>[N:1]").unwrap();
    let sets = source.reaction_products(&mut rxn, 0).unwrap();
    assert!(rxn.is_initialized());
    assert_eq!(sets.len(), 1);
    assert_eq!(sets[0].len(), 1);
    assert_eq!(sets[0][0].topology().atoms[0].element(), Element::N);
    assert_eq!(source.topology(), peer.topology());
    assert_eq!(source.coordinates_2d(), peer.coordinates_2d());
    assert_eq!(source.conformers_3d(), peer.conformers_3d());
    assert_eq!(source.properties(), peer.properties());
    let mut pair = Reaction::from_smirks("[C:1].[C:2]>>[C:1].[C:2]").unwrap();
    let sets = source
        .reaction_products_from_inputs(
            &mut pair,
            &[&source, &source],
            &ReactionRunParams::default(),
        )
        .unwrap();
    assert_eq!(sets.len(), 1);
    assert_eq!(sets[0].len(), 2);
    assert_eq!(source.topology(), peer.topology());
    assert!(matches!(
        source.reaction_products_from_inputs(&mut pair, &[&source], &ReactionRunParams::default()),
        Err(OperationError::ReactionRun(
            ReactionRunError::ReactantArity {
                expected: 2,
                actual: 1
            }
        ))
    ));
}

#[test]
fn apply_value_and_inplace_preserve_atomicity_and_source_bool() {
    let mut source = Molecule::from_smiles("C").unwrap();
    let peer = source.clone();
    let mut rxn = Reaction::from_smirks("[C:1]>>[N:1]").unwrap();
    let result = source.apply_reaction(&mut rxn).unwrap();
    assert!(result.changed());
    assert_eq!(result.molecule().topology().atoms[0].element(), Element::N);
    assert_eq!(source.topology(), peer.topology());
    assert!(source.apply_reaction_(&mut rxn).unwrap());
    assert_eq!(source.topology().atoms[0].element(), Element::N);
    assert_eq!(peer.topology().atoms[0].element(), Element::C);
    let before = source.clone();
    let mut growth = Reaction::from_smirks("[N:1]>>[N:1]C").unwrap();
    assert!(matches!(
        source.apply_reaction_(&mut growth),
        Err(OperationError::ReactionApply(
            ReactionApplyError::AddsProductAtom { .. }
        ))
    ));
    assert_eq!(source.topology(), before.topology());
    assert_eq!(source.coordinates_2d(), before.coordinates_2d());
    assert_eq!(source.conformers_3d(), before.conformers_3d());
    assert_eq!(source.properties(), before.properties());
}

#[test]
fn reaction_parameter_properties_and_contracts_are_registered() {
    assert_eq!(ReactionParseParams::default().allow_cxsmiles(), true);
    assert_eq!(ReactionRunParams::default().max_products(), 1000);
    assert_eq!(ReactionCoordinateSelection::three_d(42).id(), Some(42));
    assert!(ReactionCoordinateSelection::three_d(42).is_3d());
    assert_eq!(ReactionCoordinateSelection::auto().id(), None);
    for name in [
        "Molecule.reaction_products",
        "Molecule.reaction_products_with_params",
        "Molecule.reaction_products_from_inputs",
        "Reaction.run",
        "Molecule.apply_reaction",
        "Molecule.apply_reaction_",
        "Molecule.apply_reaction_with_params",
        "Molecule.apply_reaction_with_params_",
    ] {
        assert!(
            BINDING_CONTRACT
                .iter()
                .any(|entry| entry.semantic_id == name),
            "{name}"
        );
    }
}
