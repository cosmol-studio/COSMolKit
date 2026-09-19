use cosmolkit::{
    Atom, AtomId, AtomSpec, BINDING_CONTRACT, BindingExposure, BindingItem, BindingOwner,
    BindingParity, BindingSupport, BindingTypeRole, StateModel, TemplateAttachment,
    TemplateAttachmentOrder, TemplateAttachmentOrderError,
};

#[test]
fn canonical_template_attachment_types_and_accessors_are_reexported() {
    let _: fn(&TemplateAttachment) -> AtomId = TemplateAttachment::target;
    let _: for<'a> fn(&'a TemplateAttachment) -> &'a str = TemplateAttachment::label;
    let _: for<'a> fn(&'a TemplateAttachmentOrder) -> &'a [TemplateAttachment] =
        TemplateAttachmentOrder::entries;
    let _: for<'a> fn(&'a Atom) -> Option<&'a TemplateAttachmentOrder> =
        Atom::template_attachment_order;
    let _: Option<TemplateAttachmentOrderError> = None;
    let _: Option<AtomSpec> = None;
}

#[test]
fn binding_contract_has_exact_template_attachment_value_and_accessor_rows() {
    let expected = [
        "types.TemplateAttachment",
        "types.TemplateAttachmentOrder",
        "types.TemplateAttachmentOrderError",
        "TemplateAttachment.target",
        "TemplateAttachment.label",
        "TemplateAttachmentOrder.entries",
        "Atom.template_attachment_order",
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
    assert_eq!(rows[0].item, BindingItem::Type);
    assert_eq!(rows[0].type_role, Some(BindingTypeRole::Value));
    assert_eq!(rows[1].item, BindingItem::Type);
    assert_eq!(rows[1].type_role, Some(BindingTypeRole::Value));
    assert_eq!(rows[2].item, BindingItem::Type);
    assert_eq!(rows[2].type_role, Some(BindingTypeRole::Error));
    for row in &rows[3..] {
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.callable.unwrap().state_model, StateModel::ReadOnly);
        assert!(row.callable.unwrap().parameters.is_empty());
        assert!(row.callable.unwrap().operation_semantic_id.is_none());
    }
    assert_eq!(rows[4].rust_path, "crate :: TemplateAttachment :: label");
    assert_eq!(rows[6].javascript_name, "templateAttachmentOrder");
}
