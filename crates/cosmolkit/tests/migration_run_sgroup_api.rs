use cosmolkit::{
    AtomId, AtomSpec, BINDING_CONTRACT, BindingExposure, BindingItem, BindingOwner, BindingParity,
    BindingSupport, BindingTypeRole, BondId, BondOrder, BondSpec, Element, MoleculeBuilder,
    OperationError, SGroupBracket, SGroupCState, SGroupDisplay, StateModel, SubstanceGroup,
    SubstanceGroupId, SubstanceGroupKind,
};

#[test]
fn canonical_sgroup_values_and_builder_accessors_are_public_with_exact_signatures() {
    let _: fn(usize) -> SubstanceGroupId = SubstanceGroupId::new;
    let _: fn(&SubstanceGroupId) -> usize = SubstanceGroupId::index;
    let _: fn(SubstanceGroupId, SubstanceGroupKind) -> SubstanceGroup = SubstanceGroup::new;
    let _: fn([[f64; 3]; 3]) -> SGroupBracket = SGroupBracket::new;
    let _: for<'a> fn(&'a SGroupBracket) -> &'a [[f64; 3]; 3] = SGroupBracket::points;
    let _: fn(BondId, [f64; 3]) -> SGroupCState = SGroupCState::new;
    let _: fn(&SGroupCState) -> BondId = SGroupCState::bond;
    let _: for<'a> fn(&'a SGroupCState) -> &'a [f64; 3] = SGroupCState::vector;
    let _: for<'a> fn(&'a SGroupDisplay) -> &'a [SGroupBracket] = SGroupDisplay::brackets;
    let _: for<'a> fn(&'a SubstanceGroup) -> &'a [BondId] = SubstanceGroup::head_crossing_bonds;
    let _: for<'a> fn(&'a SubstanceGroup) -> &'a [BondId] =
        SubstanceGroup::crossing_bond_correspondence;
    let _: for<'a> fn(&'a SubstanceGroup) -> Option<&'a SGroupDisplay> = SubstanceGroup::display;
    let _: for<'a> fn(&'a SubstanceGroup) -> &'a [SGroupCState] = SubstanceGroup::cstates;
    let _: for<'a> fn(&'a MoleculeBuilder) -> &'a [SubstanceGroup] =
        MoleculeBuilder::substance_groups;
    let _: fn(&mut MoleculeBuilder, SubstanceGroup) -> Result<SubstanceGroupId, OperationError> =
        MoleculeBuilder::add_substance_group;
}

#[test]
fn canonical_sgroup_xyz_and_typed_bond_references_are_observable_through_builder() {
    let points = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
    let bracket = SGroupBracket::new(points);
    let cstate = SGroupCState::new(BondId::new(0), [-1.0, -2.0, -3.0]);
    let group = SubstanceGroup::new(SubstanceGroupId::new(91), SubstanceGroupKind::Superatom)
        .with_atoms(vec![AtomId::new(0), AtomId::new(1)])
        .with_bonds(vec![BondId::new(0)])
        .with_display(SGroupDisplay {
            brackets: vec![bracket],
            ..Default::default()
        })
        .with_cstates(vec![cstate])
        .with_head_crossing_bonds(vec![BondId::new(0), BondId::new(0)])
        .with_crossing_bond_correspondence(vec![BondId::new(0), BondId::new(0), BondId::new(0)]);

    let mut builder = MoleculeBuilder::new();
    let carbon = builder.add_atom(AtomSpec::new(Element::C));
    let oxygen = builder.add_atom(AtomSpec::new(Element::O));
    builder
        .add_bond(BondSpec::new(carbon, oxygen, BondOrder::Single))
        .expect("valid bond");
    let id = builder
        .add_substance_group(group)
        .expect("canonical SGroup is valid");
    assert_eq!(id.index(), 0);

    let observed = &builder.substance_groups()[0];
    assert_eq!(observed.display().unwrap().brackets()[0].points(), &points);
    assert_eq!(observed.cstates()[0].bond(), BondId::new(0));
    assert_eq!(observed.cstates()[0].vector(), &[-1.0, -2.0, -3.0]);
    assert_eq!(
        observed.head_crossing_bonds(),
        &[BondId::new(0), BondId::new(0)]
    );
    assert_eq!(
        observed.crossing_bond_correspondence(),
        &[BondId::new(0), BondId::new(0), BondId::new(0)]
    );
    for raw_sidecar in ["XBHEAD", "XBCORR", "_headCrossings", "_tailCrossings"] {
        assert_eq!(observed.props().get(raw_sidecar), None, "{raw_sidecar}");
    }
}

#[test]
fn binding_contract_has_only_canonical_sgroup_values_constructors_and_readers() {
    let expected = [
        "types.SubstanceGroupId",
        "types.SubstanceGroupKind",
        "types.SubstanceGroup",
        "types.SGroupBracket",
        "types.SGroupCState",
        "types.SGroupDisplay",
        "SubstanceGroupId.new",
        "SubstanceGroup.new",
        "SubstanceGroupId.index",
        "SGroupBracket.new",
        "SGroupCState.new",
        "SGroupBracket.points",
        "SGroupCState.bond",
        "SGroupCState.vector",
        "SGroupDisplay.brackets",
        "SubstanceGroup.id",
        "SubstanceGroup.kind",
        "SubstanceGroup.display",
        "SubstanceGroup.cstates",
        "SubstanceGroup.head_crossing_bonds",
        "SubstanceGroup.crossing_bond_correspondence",
    ];
    let rows = BINDING_CONTRACT
        .iter()
        .filter(|row| expected.contains(&row.semantic_id))
        .collect::<Vec<_>>();
    assert_eq!(
        rows.iter().map(|row| row.semantic_id).collect::<Vec<_>>(),
        expected
    );
    assert!(rows.iter().all(|row| row.owner == BindingOwner::Type));
    assert!(
        rows.iter()
            .all(|row| row.exposure == BindingExposure::Public)
    );
    assert!(rows.iter().all(|row| row.feature == "runtime"));
    assert!(
        rows.iter()
            .all(|row| row.support == BindingSupport::Supported)
    );
    assert!(
        rows.iter()
            .all(|row| row.parity == BindingParity::NotApplicable)
    );
    assert!(rows[..6].iter().all(|row| {
        row.item == BindingItem::Type && row.type_role == Some(BindingTypeRole::Value)
    }));
    for index in [6, 7, 9, 10] {
        let row = rows[index];
        assert!(
            row.item == BindingItem::Callable
                && row.callable.unwrap().state_model == StateModel::ValueReturning
                && row.callable.unwrap().operation_semantic_id.is_none()
        );
    }
    for index in std::iter::once(8).chain(11..rows.len()) {
        let row = rows[index];
        assert!(
            row.item == BindingItem::Callable
                && row.callable.unwrap().state_model == StateModel::ReadOnly
                && row.callable.unwrap().operation_semantic_id.is_none()
        );
    }
    assert_eq!(rows[19].javascript_name, "headCrossingBonds");
    assert_eq!(rows[20].javascript_name, "crossingBondCorrespondence");

    assert!(BINDING_CONTRACT.iter().all(|row| {
        !["XBHEAD", "XBCORR", "_headCrossings", "_tailCrossings"]
            .iter()
            .any(|legacy| row.semantic_id.contains(legacy))
    }));
}
