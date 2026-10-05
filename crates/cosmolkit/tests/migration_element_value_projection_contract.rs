//! Exact public value-callable contracts and result field schema witnesses.
use cosmolkit as ck;
const _: fn(u8) -> Option<ck::Element> = ck::Element::from_atomic_number;
const _: for<'a> fn(&'a str) -> Option<ck::Element> = ck::Element::from_symbol;
const _: fn(ck::Element) -> u8 = ck::Element::atomic_number;
const _: fn(ck::Element) -> &'static str = ck::Element::symbol;
#[cfg(feature = "cap-valence")]
const _: fn(ck::Element) -> ck::ElementInfo = ck::element_info;
// Result schema witnesses refer to actual fields, not nonexistent Rust getters.
const _: for<'a> fn(&'a ck::ElementInfo) -> ck::Element = |info| info.element;
const _: for<'a> fn(&'a ck::ElementInfo) -> &'static str = |info| info.symbol;
const _: for<'a> fn(&'a ck::ElementInfo) -> u8 = |info| info.atomic_number;
const _: for<'a> fn(&'a ck::ElementInfo) -> u8 = |info| info.period;
const _: for<'a> fn(&'a ck::ElementInfo) -> i32 = |info| info.outer_electrons;
const _: for<'a> fn(&'a ck::ElementInfo) -> &'static [i32] = |info| info.valences;
const _: for<'a> fn(&'a ck::ElementInfo) -> f64 = |info| info.rb0;
const _: for<'a> fn(&'a ck::ElementInfo) -> f64 = |info| info.atomic_weight;

fn contract(
    method: &str,
    python: &str,
    javascript: &str,
    kind: ck::BindingKind,
    receiver: Option<ck::BindingReceiver>,
    parameters: &[(&str, &str)],
    output: &str,
) {
    let matches: Vec<_> = ck::BINDING_CONTRACT
        .iter()
        .filter(|entry| entry.semantic_id == format!("Element.{method}"))
        .collect();
    assert_eq!(matches.len(), 1);
    let row = matches[0];
    assert_eq!(row.item, ck::BindingItem::Callable);
    assert_eq!(row.owner, ck::BindingOwner::Type);
    assert_eq!(
        row.rust_path.replace(' ', ""),
        format!("crate::Element::{method}")
    );
    assert_eq!(row.python_name, python);
    assert_eq!(row.javascript_name, javascript);
    assert_eq!(row.feature, "metadata");
    assert_eq!(row.status, ck::FunctionStatus::Experimental);
    assert_eq!(row.type_role, None);
    let call = row.callable.unwrap();
    assert_eq!(call.kind, kind);
    assert_eq!(call.receiver, receiver);
    assert_eq!(call.output_type.replace(' ', ""), output);
    assert_eq!(call.error_type, None);
    assert_eq!(call.operation_semantic_id, None);
    assert_eq!(call.state_model, ck::StateModel::ValueReturning);
    assert_eq!(call.parameters.len(), parameters.len());
    for (actual, &(name, ty)) in call.parameters.iter().zip(parameters) {
        assert_eq!(actual.name, name);
        assert_eq!(actual.type_name.replace(' ', ""), ty);
        assert_eq!(actual.default, ck::BindingDefault::Required);
    }
}

#[test]
fn number_constructor_is_an_exact_static_optional_value() {
    contract(
        "from_atomic_number",
        "from_atomic_number",
        "fromAtomicNumber",
        ck::BindingKind::Static,
        None,
        &[("atomic_number", "u8")],
        "Option<crate::Element>",
    );
    for number in 0..=118 {
        assert_eq!(
            ck::Element::from_atomic_number(number)
                .unwrap()
                .atomic_number(),
            number
        );
    }
    for number in 119..=255 {
        assert!(ck::Element::from_atomic_number(number).is_none());
    }
}

#[test]
fn symbol_constructor_is_an_exact_static_optional_value() {
    contract(
        "from_symbol",
        "from_symbol",
        "fromSymbol",
        ck::BindingKind::Static,
        None,
        &[("symbol", "&str")],
        "Option<crate::Element>",
    );
    for element in ck::Element::iter_with_dummy() {
        assert_eq!(ck::Element::from_symbol(element.symbol()), Some(element));
    }
    assert_eq!(ck::Element::from_symbol("Uut"), Some(ck::Element::NH));
    assert_eq!(ck::Element::from_symbol("Uup"), Some(ck::Element::MC));
    for symbol in ["", "c", " C", "C ", "Uuo"] {
        assert!(ck::Element::from_symbol(symbol).is_none());
    }
}

#[test]
fn atomic_number_accessor_has_an_owned_copy_receiver() {
    contract(
        "atomic_number",
        "atomic_number",
        "atomicNumber",
        ck::BindingKind::Instance,
        Some(ck::BindingReceiver::Owned),
        &[],
        "u8",
    );
    let element = ck::Element::C;
    for _ in 0..3 {
        assert_eq!(element.atomic_number(), 6);
    }
}

#[test]
fn symbol_accessor_has_an_owned_copy_receiver() {
    contract(
        "symbol",
        "symbol",
        "symbol",
        ck::BindingKind::Instance,
        Some(ck::BindingReceiver::Owned),
        &[],
        "&'staticstr",
    );
    let element = ck::Element::C;
    for _ in 0..3 {
        assert_eq!(element.symbol(), "C");
    }
}

#[cfg(feature = "cap-valence")]
#[test]
fn result_fields_and_module_callable_preserve_existing_owner_contract() {
    let row = ck::BINDING_CONTRACT
        .iter()
        .find(|row| row.semantic_id == "module.element_info")
        .unwrap();
    assert_eq!(row.owner, ck::BindingOwner::Module);
    assert_eq!(row.item, ck::BindingItem::Callable);
    assert_eq!(row.feature, "cap-valence");
    assert_eq!(row.status, ck::FunctionStatus::Experimental);
    assert_eq!(row.python_name, "element_info");
    assert_eq!(row.javascript_name, "elementInfo");
    assert_eq!(row.rust_path.replace(' ', ""), "crate::element_info");
    let call = row.callable.unwrap();
    assert_eq!(call.kind, ck::BindingKind::Module);
    assert_eq!(call.receiver, None);
    assert_eq!(call.state_model, ck::StateModel::ReadOnly);
    assert_eq!(call.output_type.replace(' ', ""), "crate::ElementInfo");
    assert_eq!(call.error_type, None);
    assert_eq!(call.operation_semantic_id, None);
    assert_eq!(call.parameters.len(), 1);
    assert_eq!(call.parameters[0].name, "element");
    assert_eq!(
        call.parameters[0].type_name.replace(' ', ""),
        "crate::Element"
    );
    assert_eq!(call.parameters[0].default, ck::BindingDefault::Required);
    assert!(
        ck::BINDING_CONTRACT
            .iter()
            .all(|row| !row.semantic_id.starts_with("ElementInfo."))
    );
    let ck::ElementInfo {
        element,
        symbol,
        atomic_number,
        period,
        outer_electrons,
        valences,
        rb0,
        atomic_weight,
    } = ck::element_info(ck::Element::C);
    assert_eq!(
        (element, symbol, atomic_number, period, outer_electrons),
        (ck::Element::C, "C", 6, 2, 4)
    );
    assert_eq!(valences, &[4]);
    assert_eq!(rb0.to_bits(), 0.77_f64.to_bits());
    assert_eq!(atomic_weight.to_bits(), 12.011_f64.to_bits());
}
