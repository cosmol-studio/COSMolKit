#[path = "../src/status.rs"]
mod status;

#[path = "../src/binding.rs"]
mod binding;

use binding::expand_binding_contract;
use quote::quote;

struct Molecule {
    value: String,
}

struct Protein;
struct BioStructure;

impl Protein {
    fn into_bio_structure(self) -> BioStructure {
        BioStructure
    }
}

#[derive(Clone, Copy)]
struct BorrowedView<'a> {
    value: &'a str,
}

impl<'a> BorrowedView<'a> {
    fn value(self) -> &'a str {
        self.value
    }
}

impl Molecule {
    fn property<'a>(&'a self, key: &str) -> Option<&'a str> {
        (key == "value").then_some(self.value.as_str())
    }
}

enum BindingItem {
    Callable,
    Type,
}
enum BindingOwner {
    Molecule,
    Module,
    Type,
}
enum BindingKind {
    Instance,
    Static,
    Module,
}
enum BindingDefault {
    Required,
    Value(&'static str),
}
struct BindingParameterContract {
    name: &'static str,
    type_name: &'static str,
    default: BindingDefault,
}
enum BindingReceiver {
    Shared,
    Mutable,
    Owned,
}
struct BindingCallableContract {
    kind: BindingKind,
    receiver: Option<BindingReceiver>,
    parameters: &'static [BindingParameterContract],
    output_type: &'static str,
    error_type: Option<&'static str>,
    state_model: StateModel,
    operation_semantic_id: Option<&'static str>,
}
enum BindingTypeRole {
    Value,
    Parameter,
    Result,
    Error,
}
enum FunctionStatus {
    Parity {
        reference: &'static str,
    },
    ParityWithDifferences {
        reference: &'static str,
        explanation: &'static str,
    },
    Native,
    Experimental,
}
enum StateModel {
    ValueReturning,
    InPlace,
    ReadOnly,
}
enum BindingPropertyAccess {
    Getter,
    Setter,
}
struct BindingContractEntry {
    semantic_id: &'static str,
    item: BindingItem,
    owner: BindingOwner,
    rust_path: &'static str,
    python_name: &'static str,
    python_property: Option<BindingPropertyAccess>,
    javascript_name: &'static str,
    feature: &'static str,
    required_capabilities: &'static [&'static str],
    status: FunctionStatus,
    callable: Option<BindingCallableContract>,
    type_role: Option<BindingTypeRole>,
}

cosmolkit_macros::binding_contract! {
    static HRTB_BINDINGS = [
        {
            semantic_id: "Molecule.property",
            item: callable,
            owner: molecule,
            rust: crate::Molecule::property,
            python: "property",
            javascript: "property",
            feature: "runtime",
            kind: instance,
            parameters: [{ name: key, type: &str, default: required }],
            output: Option<&str>,
            error: none,
            state: read_only,
            operation: none,
            signature: for<'a, 'b> fn(&'a crate::Molecule, &'b str) -> Option<&'a str>,
        }
    ];
}

cosmolkit_macros::binding_contract! {
    static GENERIC_RECEIVER_BINDINGS = [
        {
            semantic_id: "BorrowedView.value",
            item: callable,
            owner: type_,
            rust: crate::BorrowedView::value,
            python: "value",
            javascript: "value",
            feature: "test",
            kind: instance, receiver: owned,
            parameters: [],
            output: &str,
            error: none,
            state: value_returning,
            operation: none,
            signature: for<'a> fn(crate::BorrowedView<'a>) -> &'a str,
        }
    ];
}

cosmolkit_macros::binding_contract! {
    static OWNED_PROTEIN_BINDINGS = [
        {
            semantic_id: "Protein.into_bio_structure",
            item: callable,
            owner: type_,
            rust: crate::Protein::into_bio_structure,
            python: "into_bio_structure",
            javascript: "intoBioStructure",
            feature: "cap-bio",
            kind: instance, receiver: owned,
            parameters: [],
            output: crate::BioStructure,
            error: none,
            state: value_returning,
            operation: none,
            signature: fn(crate::Protein) -> crate::BioStructure,
        }
    ];
}

fn compact(tokens: proc_macro2::TokenStream) -> String {
    tokens
        .to_string()
        .chars()
        .filter(|character| !character.is_whitespace())
        .collect()
}

fn callable_entry() -> &'static str {
    r#"#[cfg(feature = "cap-descriptors")]
    {
        semantic_id: "Molecule.molecular_weight",
        item: callable,
        owner: molecule,
        rust: crate::Molecule::molecular_weight,
        python: "molecular_weight",
        javascript: "molecularWeight",
        feature: "cap-descriptors",
        kind: instance,
        parameters: [],
        output: f64,
        error: crate::OperationError,
        state: read_only,
        operation: none,
        signature: fn(&crate::Molecule) -> Result<f64, crate::OperationError>,
    }"#
}

fn registry_with(entries: &str) -> proc_macro2::TokenStream {
    format!("pub static API = [{entries}];")
        .parse()
        .expect("tokens")
}

trait TestTokens {
    fn into_tokens(self) -> proc_macro2::TokenStream;
}

impl TestTokens for proc_macro2::TokenStream {
    fn into_tokens(self) -> proc_macro2::TokenStream {
        self
    }
}

impl TestTokens for &str {
    fn into_tokens(self) -> proc_macro2::TokenStream {
        self.parse().expect("tokens")
    }
}

impl TestTokens for String {
    fn into_tokens(self) -> proc_macro2::TokenStream {
        self.parse().expect("tokens")
    }
}

fn error_for(source: impl TestTokens) -> String {
    let tokens = source.into_tokens();
    match expand_binding_contract(tokens) {
        Ok(_) => panic!("fixture unexpectedly expanded"),
        Err(error) => error.to_string(),
    }
}

#[test]
fn empty_registry_emits_one_empty_canonical_slice() {
    let generated = compact(expand_binding_contract(quote!(pub static API = [];)).unwrap());
    assert_eq!(
        generated,
        "pubstaticAPI:&[crate::BindingContractEntry]=&[];"
    );
}

#[test]
fn mixed_callable_and_type_entries_preserve_order_and_entry_local_cfg() {
    let type_entry = r#"#[cfg(feature = "cap-search")]
    {
        semantic_id: "types.QueryGraph", item: type, owner: type_,
        rust: crate::QueryGraph, python: "QueryGraph", javascript: "QueryGraph",
        feature: "cap-search", role: value,
    }"#;
    let generated = compact(
        expand_binding_contract(registry_with(&format!(
            "{},{}",
            callable_entry(),
            type_entry
        )))
        .unwrap(),
    );
    let callable = generated.find("Molecule.molecular_weight").unwrap();
    let value_type = generated.find("types.QueryGraph").unwrap();
    assert!(callable < value_type);
    assert!(generated.contains("#[cfg(feature=\"cap-descriptors\")]crate::BindingContractEntry"));
    assert!(generated.contains("#[cfg(feature=\"cap-search\")]crate::BindingContractEntry"));
    assert!(generated.contains("item:crate::BindingItem::Callable"));
    assert!(generated.contains("item:crate::BindingItem::Type"));
    assert!(generated.contains("type_role:Some(crate::BindingTypeRole::Value)"));
}

#[test]
fn every_entry_uses_the_full_signature_without_an_exposure_switch() {
    let generated = compact(expand_binding_contract(registry_with(callable_entry())).unwrap());
    assert!(generated.contains(
        "const__BINDING_ASSERT_API_0:fn(&crate::Molecule)->Result<f64,crate::OperationError>=crate::Molecule::molecular_weight;"
    ));
}

#[test]
fn linked_molecule_status_is_generated_from_the_operation_not_redeclared() {
    let source = callable_entry().replace("operation: none", "operation: \"molecular_weight\"");
    let generated = compact(expand_binding_contract(registry_with(&source)).unwrap());
    assert!(generated.contains("status:crate::Molecule::__FUNCTION_STATUS_MOLECULAR_WEIGHT"));
    for status in [
        "experimental",
        "native",
        "parity(\"RDKit\")",
        "parity_with_differences(\"RDKit\", \"approved conditions\")",
    ] {
        let repeated = source.replace(
            "feature: \"cap-descriptors\",",
            &format!("feature: \"cap-descriptors\", status: {status},"),
        );
        assert!(error_for(registry_with(&repeated)).contains("status is inherited"));
    }
}

#[test]
fn higher_ranked_borrowed_instance_output_compiles_with_exact_signature() {
    let molecule = Molecule {
        value: "kept".to_owned(),
    };
    let property: for<'a, 'b> fn(&'a Molecule, &'b str) -> Option<&'a str> = Molecule::property;
    assert_eq!(property(&molecule, "value"), Some("kept"));
    assert_eq!(HRTB_BINDINGS.len(), 1);
    assert_eq!(HRTB_BINDINGS[0].semantic_id, "Molecule.property");
}

#[test]
fn lifetime_parameterized_type_owned_receiver_keeps_exact_public_signature() {
    fn assert_value_signature<'a>() {
        let _: fn(BorrowedView<'a>) -> &'a str = BorrowedView::value;
    }

    let view = BorrowedView { value: "kept" };
    assert_value_signature();
    assert_eq!(BorrowedView::value(view), "kept");
    assert_eq!(GENERIC_RECEIVER_BINDINGS.len(), 1);
    assert_eq!(
        GENERIC_RECEIVER_BINDINGS[0].semantic_id,
        "BorrowedView.value"
    );
}

#[test]
fn lifetime_parameterized_type_owned_receiver_still_rejects_wrong_base_type() {
    let invalid = r#"{
        semantic_id: "BorrowedView.value", item: callable, owner: type_,
        rust: crate::BorrowedView::value, python: "value", javascript: "value",
        feature: "test", kind: instance, receiver: owned, parameters: [], output: &str,
        error: none, state: value_returning, operation: none,
        signature: for<'a> fn(crate::Molecule) -> &'a str,
    }"#;
    assert!(
        error_for(registry_with(invalid))
            .contains("binding signature has the wrong instance receiver")
    );
}

#[test]
fn unsafe_abi_and_variadic_binding_signatures_remain_rejected() {
    let base = callable_entry();
    for invalid in [
        base.replace("signature: fn", "signature: unsafe fn"),
        base.replace("signature: fn", "signature: extern \"C\" fn"),
        base.replace("fn(&crate::Molecule)", "fn(&crate::Molecule, ...)"),
    ] {
        assert!(
            error_for(registry_with(&invalid))
                .contains("safe non-variadic Rust fn with the default ABI")
        );
    }
}

#[test]
fn instance_static_and_module_signatures_cover_parameters_defaults_and_direct_outputs() {
    let constructor = r#"{
        semantic_id: "Molecule.from_smiles", item: callable, owner: molecule,
        rust: crate::Molecule::from_smiles, python: "from_smiles", javascript: "fromSmiles",
        feature: "cap-smiles", kind: static_,
        parameters: [{ name: text, type: &str, default: required }],
        output: crate::Molecule, error: crate::OperationError, state: value_returning,
        operation: none,
        signature: fn(&str) -> Result<crate::Molecule, crate::OperationError>,
    }"#;
    let configured = r#"{
        semantic_id: "Molecule.to_smiles_with_params", item: callable, owner: molecule,
        rust: crate::Molecule::to_smiles_with_params,
        python: "to_smiles_with_params", javascript: "toSmilesWithParams",
        feature: "cap-smiles", kind: instance,
        parameters: [
            { name: params, type: &crate::SmilesWriteParams, default: required },
            { name: include_cx, type: bool, default: false },
        ],
        output: String, error: crate::OperationError, state: read_only, operation: none,
        signature: fn(&crate::Molecule, &crate::SmilesWriteParams, bool)
            -> Result<String, crate::OperationError>,
    }"#;
    let module = r#"{
        semantic_id: "module.version", item: callable, owner: module,
        rust: crate::version, python: "version", javascript: "version",
        feature: "metadata", kind: module, parameters: [],
        output: &'static str, error: none, state: read_only, operation: none,
        signature: fn() -> &'static str,
    }"#;
    let generated = compact(
        expand_binding_contract(registry_with(&format!(
            "{constructor},{configured},{module}"
        )))
        .unwrap(),
    );
    assert!(generated.contains("kind:crate::BindingKind::Static"));
    assert!(generated.contains("kind:crate::BindingKind::Instance"));
    assert!(generated.contains("kind:crate::BindingKind::Module"));
    assert!(generated.contains("default:crate::BindingDefault::Required"));
    assert!(generated.contains("default:crate::BindingDefault::Value(stringify!(false))"));
    assert!(generated.contains("error_type:None"));
}

#[test]
fn all_type_roles_and_four_function_statuses_are_structured() {
    let roles = ["value", "parameter", "result", "error"];
    let statuses = [
        "parity(\"RDKit\")",
        "parity_with_differences(\"Gemmi\", \"Approved difference under explicitly documented conditions\")",
        "native",
        "experimental",
    ];
    let mut entries = Vec::new();
    for role in roles {
        for status in statuses {
            let index = entries.len();
            entries.push(format!("{{semantic_id:\"types.Type{index}\",item:type,owner:type_,rust:crate::Type{index},python:\"Type{index}\",javascript:\"Type{index}\",feature:\"types\",status:{status},role:{role}}}"));
        }
    }
    assert_eq!(entries.len(), 16);
    let entries = entries.join(",");
    let generated = compact(expand_binding_contract(registry_with(&entries)).unwrap());
    for variant in ["Value", "Parameter", "Result", "Error"] {
        assert!(generated.contains(&format!("BindingTypeRole::{variant}")));
    }
    for variant in ["Parity", "ParityWithDifferences", "Native", "Experimental"] {
        assert!(generated.contains(&format!("FunctionStatus::{variant}")));
    }
    assert!(generated.contains("reference:\"RDKit\""));
    assert!(generated.contains("reference:\"Gemmi\""));
    assert_eq!(generated.matches("assert_public_type::<").count(), 16);
}

#[test]
fn unknown_missing_duplicate_and_wrong_item_fields_fail_closed() {
    let base = callable_entry();
    assert!(
        error_for(registry_with(
            &base.replace("item: callable,", "mystery: callable,")
        ))
        .contains("unknown binding entry field `mystery`")
    );
    assert!(
        error_for(registry_with(
            &base.replace("feature: \"cap-descriptors\",", "")
        ))
        .contains("missing `feature`")
    );
    assert!(
        error_for(registry_with(
            &base.replace("item: callable,", "item: callable, item: callable,")
        ))
        .contains("duplicate binding field `item`")
    );
    let type_with_kind = r#"pub static API = [{semantic_id:"types.QueryGraph",item:type,owner:type_,rust:crate::QueryGraph,python:"QueryGraph",javascript:"QueryGraph",feature:"cap-search",role:value,kind:instance}];"#;
    assert!(error_for(type_with_kind).contains("type binding entry cannot declare `kind`"));
}

#[test]
fn duplicate_id_and_each_owner_scoped_projection_fail_closed() {
    let first = callable_entry();
    let second = first.replace("molecular_weight", "exact_molecular_weight");
    let duplicate_id = format!("pub static API = [{first},{first}];");
    assert!(error_for(duplicate_id).contains("duplicate binding semantic_id"));

    let duplicate_rust = second.replace(
        "crate::Molecule::exact_molecular_weight",
        "crate::Molecule::molecular_weight",
    );
    assert!(
        error_for(format!("pub static API = [{first},{duplicate_rust}];"))
            .contains("duplicate Rust binding projection")
    );
    let duplicate_python = second.replace(
        "python: \"exact_molecular_weight\"",
        "python: \"molecular_weight\"",
    );
    assert!(
        error_for(format!("pub static API = [{first},{duplicate_python}];"))
            .contains("Python callable name must equal")
    );
    let duplicate_javascript = second.replace(
        "javascript: \"exactMolecularWeight\"",
        "javascript: \"molecularWeight\"",
    );
    assert!(
        error_for(format!(
            "pub static API = [{first},{duplicate_javascript}];"
        ))
        .contains("JavaScript callable name must")
    );
}

#[test]
fn receiver_state_and_trailing_underscore_contracts_fail_closed() {
    let base = callable_entry();
    assert!(
        error_for(registry_with(
            &base.replace("owner: molecule", "owner: module")
        ))
        .contains("invalid callable owner/kind")
    );
    assert!(
        error_for(registry_with(
            &base.replace("state: read_only", "state: in_place")
        ))
        .contains("only in-place callable names")
    );
    let misplaced_suffix = base.replace("molecular_weight", "molecular_weight_");
    assert!(error_for(registry_with(&misplaced_suffix)).contains("only in-place callable names"));
    let bad_receiver = base.replace("fn(&crate::Molecule)", "fn(&mut crate::Molecule)");
    assert!(error_for(registry_with(&bad_receiver)).contains("wrong instance receiver"));
}

#[test]
fn type_owned_constructor_names_are_scoped_to_the_actual_receiver_type() {
    let constructor = |type_name: &str| {
        format!(
            r#"{{
                semantic_id: "{type_name}.new", item: callable, owner: type_,
                rust: crate::{type_name}::new, python: "new", javascript: "new",
                feature: "test", kind: static_, parameters: [],
                output: crate::{type_name}, error: none, state: value_returning,
                operation: none, signature: fn() -> crate::{type_name}
            }}"#
        )
    };
    let input = format!(
        "pub static API = [{},{}];",
        constructor("FirstValue"),
        constructor("SecondValue")
    );
    assert!(expand_binding_contract(input.parse().unwrap()).is_ok());
}

#[test]
fn type_declaration_projection_collision_is_rejected_in_the_export_scope() {
    let first = r#"{
        semantic_id:"first.Value",item:type,owner:type_,rust:crate::first::Value,
        python:"Value",javascript:"Value",feature:"test",role:value
    }"#;
    let second = r#"{
        semantic_id:"second.Value",item:type,owner:type_,rust:crate::second::Value,
        python:"Value",javascript:"Value",feature:"test",role:value
    }"#;
    let error = error_for(format!("pub static API = [{first},{second}];"));
    assert!(
        error.contains("duplicate Python binding projection")
            || error.contains("duplicate JavaScript binding projection"),
        "unexpected error: {error}"
    );
}

#[test]
fn same_receiver_callable_projection_collision_is_rejected() {
    let constructor = r#"{
        semantic_id:"FirstValue.new",item:callable,owner:type_,
        rust:crate::FirstValue::new,python:"new",javascript:"new",
        feature:"test",kind:static_,parameters:[],output:crate::FirstValue,error:none,
        state:value_returning,operation:none,signature:fn()->crate::FirstValue
    }"#;
    let error = error_for(format!("pub static API = [{constructor},{constructor}];"));
    assert!(
        error.contains("duplicate binding semantic_id")
            || error.contains("duplicate Rust binding projection")
            || error.contains("duplicate Python binding projection")
            || error.contains("duplicate JavaScript binding projection"),
        "unexpected error: {error}"
    );
}

#[test]
fn parameter_default_and_full_signature_disagreements_fail_closed() {
    let parameterized = r#"pub static API = [{
        semantic_id:"Molecule.configured",item:callable,owner:molecule,
        rust:crate::Molecule::configured,python:"configured",javascript:"configured",
        feature:"test",kind:instance,parameters:[
            {name:first,type:usize,default:1usize},
            {name:second,type:bool,default:required}
        ],output:usize,error:none,state:read_only,operation:none,
        signature:fn(&crate::Molecule,usize,bool)->usize
    }];"#;
    assert!(error_for(parameterized).contains("required parameter cannot follow"));
    let wrong_count =
        callable_entry().replace("fn(&crate::Molecule)", "fn(&crate::Molecule, bool)");
    assert!(error_for(registry_with(&wrong_count)).contains("argument count disagrees"));
    let wrong_return = callable_entry().replace(
        "Result<f64, crate::OperationError>",
        "Result<String, crate::OperationError>",
    );
    assert!(error_for(registry_with(&wrong_return)).contains("return type disagrees"));
}

#[test]
fn operation_identity_canonical_names_and_projections_fail_closed() {
    let operation = callable_entry().replace("operation: none", "operation: \"different\"");
    assert!(error_for(registry_with(&operation)).contains("operation id must equal"));
    let legacy = callable_entry().replace("molecular_weight", "calc_molecular_weight");
    assert!(error_for(registry_with(&legacy)).contains("non-canonical public prefix"));
    let python = callable_entry().replace("python: \"molecular_weight\"", "python: \"weight\"");
    assert!(error_for(registry_with(&python)).contains("Python callable name must equal"));
    let javascript = callable_entry().replace(
        "javascript: \"molecularWeight\"",
        "javascript: \"molecular_weight\"",
    );
    assert!(error_for(registry_with(&javascript)).contains("JavaScript callable name must"));
}

#[test]
fn python_properties_and_constructor_are_explicit_validated_projections() {
    let getter = r#"{
        semantic_id:"Params.limit", item:callable, owner:type_,
        rust:crate::Params::limit, python:"limit", python_property:getter,
        javascript:"limit", feature:"test", kind:instance, receiver:shared,
        parameters:[], output:usize, error:none, state:read_only, operation:none,
        signature:fn(&crate::Params)->usize,
    }"#;
    let setter = r#"{
        semantic_id:"Params.set_limit", item:callable, owner:type_,
        rust:crate::Params::set_limit, python:"limit", python_property:setter,
        javascript:"setLimit", feature:"test", kind:instance, receiver:mutable,
        parameters:[{name:value,type:usize,default:required}],
        output:(), error:none, state:in_place, operation:none,
        signature:fn(&mut crate::Params,usize)->(),
    }"#;
    let constructor = r#"{
        semantic_id:"Params.new", item:callable, owner:type_,
        rust:crate::Params::new, python:"__new__", javascript:"new",
        feature:"test", kind:static_, parameters:[], output:crate::Params,
        error:none, state:value_returning, operation:none, signature:fn()->crate::Params,
    }"#;
    let valid = format!("{getter},{setter},{constructor}");
    let generated = compact(expand_binding_contract(registry_with(&valid)).unwrap());
    assert!(generated.contains("Some(crate::BindingPropertyAccess::Getter)"));
    assert!(generated.contains("Some(crate::BindingPropertyAccess::Setter)"));
    assert!(
        error_for(registry_with(
            &setter.replace("python_property:setter", "python_property:getter")
        ))
        .contains("getter requires")
    );
    assert!(
        error_for(registry_with(
            &getter.replace("python_property:getter", "python_property:setter")
        ))
        .contains("setter requires")
    );
    assert!(
        error_for(registry_with(
            &constructor.replace("python:\"__new__\"", "python:\"__init__\"")
        ))
        .contains("Python callable name must equal")
    );
}

#[test]
fn cfg_features_and_invalid_function_statuses_fail_closed() {
    let cfg = callable_entry().replace("feature = \"cap-descriptors\"", "feature = \"other\"");
    assert!(error_for(registry_with(&cfg)).contains("cfg feature disagrees"));
    let compound = callable_entry().replace(
        "feature = \"cap-descriptors\"",
        "all(feature = \"cap-descriptors\")",
    );
    assert!(error_for(registry_with(&compound)).contains("compound or non-feature"));
    for invalid in [
        "supported",
        "supported_with_rdkit_parity",
        "unsupported",
        "preserved_only",
    ] {
        let source = callable_entry().replace(
            "feature: \"cap-descriptors\",",
            &format!("feature: \"cap-descriptors\", status: {invalid},"),
        );
        assert!(error_for(registry_with(&source)).contains("status must be"));
    }
    for (status, expected) in [
        ("parity(\"\")", "reference must not be empty"),
        ("parity(\"  \")", "reference must not be empty"),
        (
            "parity_with_differences(\"Gemmi\", \"\")",
            "explanation must not be empty",
        ),
        (
            "parity_with_differences(\"\", \"reason\")",
            "reference must not be empty",
        ),
        (
            "parity(\"RDKit\", \"extra\")",
            "unexpected function status arguments",
        ),
    ] {
        let source = callable_entry().replace(
            "feature: \"cap-descriptors\",",
            &format!("feature: \"cap-descriptors\", status: {status},"),
        );
        assert!(error_for(registry_with(&source)).contains(expected));
    }
    for field in [
        "exposure: registered",
        "exposure: public",
        "support: experimental",
        "parity: required_now",
    ] {
        let source = callable_entry().replace(
            "feature: \"cap-descriptors\",",
            &format!("feature: \"cap-descriptors\", {field},"),
        );
        assert!(error_for(registry_with(&source)).contains("unknown binding entry field"));
    }
}

#[test]
fn generated_contract_contains_no_runtime_registry_wrapper_or_domain_implementation() {
    let generated = compact(expand_binding_contract(registry_with(callable_entry())).unwrap());
    for forbidden in [
        "structMolecule",
        "structOpParts",
        "MOLECULE_OPS",
        "SUPPORT_MATRIX",
        "finish_in_place",
        "cosmolkit_core",
        "hydrogens::",
        "compat",
    ] {
        assert!(!generated.contains(forbidden), "unexpected `{forbidden}`");
    }
}

#[test]
fn type_owned_callables_cover_static_borrowed_mutable_and_consuming_receivers() {
    let entries = r#"
    {
        semantic_id:"MoleculeBuilder.new",item:callable,owner:type_,
        rust:crate::MoleculeBuilder::new,python:"new",javascript:"new",
        feature:"runtime",kind:static_,parameters:[],output:crate::MoleculeBuilder,error:none,
        state:value_returning,operation:none,
        signature:fn()->crate::MoleculeBuilder
    },
    {
        semantic_id:"MoleculeBuilder.len",item:callable,owner:type_,
        rust:crate::MoleculeBuilder::len,python:"len",javascript:"len",
        feature:"runtime",kind:instance,parameters:[],output:usize,error:none,state:read_only,
        operation:none,signature:fn(&crate::MoleculeBuilder)->usize
    },
    {
        semantic_id:"MoleculeBuilder.push",item:callable,owner:type_,
        rust:crate::MoleculeBuilder::push,python:"push",javascript:"push",
        feature:"runtime",kind:instance,parameters:[{name:value,type:usize,default:required}],
        output:(),error:none,state:in_place,operation:none,
        signature:fn(&mut crate::MoleculeBuilder,usize)->()
    },
    {
        semantic_id:"MoleculeBuilder.build",item:callable,owner:type_,
        rust:crate::MoleculeBuilder::build,python:"build",javascript:"build",
        feature:"runtime",kind: instance, receiver: owned,parameters:[],output:crate::Molecule,
        error:crate::OperationError,state:value_returning,operation:none,
        signature:fn(crate::MoleculeBuilder)->Result<crate::Molecule,crate::OperationError>
    }"#;
    let generated = compact(expand_binding_contract(registry_with(entries)).unwrap());
    assert!(generated.contains("fn()->crate::MoleculeBuilder=crate::MoleculeBuilder::new"));
    assert!(generated.contains("fn(&crate::MoleculeBuilder)->usize=crate::MoleculeBuilder::len"));
    assert!(
        generated.contains("fn(&mutcrate::MoleculeBuilder,usize)->()=crate::MoleculeBuilder::push")
    );
    assert!(generated.contains(
        "fn(crate::MoleculeBuilder)->Result<crate::Molecule,crate::OperationError>=crate::MoleculeBuilder::build"
    ));
}

#[test]
fn invalid_type_owned_kind_and_receiver_shapes_fail_closed() {
    let valid = r#"{
        semantic_id:"MoleculeBuilder.push",item:callable,owner:type_,
        rust:crate::MoleculeBuilder::push,python:"push",javascript:"push",
        feature:"runtime",kind:instance,parameters:[{name:value,type:usize,default:required}],
        output:(),error:none,state:in_place,operation:none,
        signature:fn(&mut crate::MoleculeBuilder,usize)->()
    }"#;
    let module_kind = valid.replace("kind:instance", "kind:module");
    assert!(error_for(registry_with(&module_kind)).contains("owner/kind combination"));

    let static_in_place = valid.replace("kind:instance", "kind:static_").replace(
        "signature:fn(&mut crate::MoleculeBuilder,usize)",
        "signature:fn(usize)",
    );
    assert!(error_for(registry_with(&static_in_place)).contains("requires an instance receiver"));

    let immutable_receiver = valid.replace(
        "signature:fn(&mut crate::MoleculeBuilder,usize)",
        "signature:fn(&crate::MoleculeBuilder,usize)",
    );
    assert!(error_for(registry_with(&immutable_receiver)).contains("wrong instance receiver"));
}

#[test]
fn in_place_javascript_projection_strips_only_the_rust_python_suffix() {
    let entry = r#"{
        semantic_id:"Molecule.set_atom_position_",item:callable,owner:molecule,
        rust:crate::Molecule::set_atom_position_,python:"set_atom_position_",
        javascript:"setAtomPosition",feature:"cap-transforms",kind:instance,
        parameters:[{name:atom,type:usize,default:required}],output:(),
        error:crate::OperationError,state:in_place,operation:"set_atom_position_",
        signature:fn(&mut crate::Molecule,usize)->Result<(),crate::OperationError>
    }"#;
    expand_binding_contract(registry_with(entry)).expect("canonical in-place projection");

    let suffixed_javascript = entry.replace(
        "javascript:\"setAtomPosition\"",
        "javascript:\"setAtomPosition_\"",
    );
    assert!(
        error_for(registry_with(&suffixed_javascript))
            .contains("JavaScript callable name must be `setAtomPosition`")
    );
}

#[test]
fn current_cosmolkit_registry_validates_every_cfg_gated_projection() {
    let file = syn::parse_file(include_str!(
        "../../cosmolkit/src/binding_contract/registry.rs"
    ))
    .expect("current binding registry parses as Rust");
    let tokens = file
        .items
        .into_iter()
        .find_map(|item| match item {
            syn::Item::Macro(item) if item.mac.path.is_ident("binding_contract") => {
                Some(item.mac.tokens)
            }
            _ => None,
        })
        .expect("one binding_contract invocation");
    expand_binding_contract(tokens).expect("all current projections satisfy canonical naming");
}

#[test]
fn bio_objects_require_inplace_suffix_and_mutable_receiver() {
    for target in ["BioStructure", "Protein"] {
        let entry = format!(
            r#"{{
            semantic_id:"{target}.translate_",item:callable,owner:type_,
            rust:crate::{target}::translate_,python:"translate_",javascript:"translate",
            feature:"cap-bio",kind:instance,parameters:[],output:(),error:none,state:in_place,operation:none,
            signature:fn(&mut crate::{target})->()
        }}"#
        );
        expand_binding_contract(registry_with(&entry)).unwrap();
        assert!(
            error_for(registry_with(&entry.replace("translate_", "translate")))
                .contains("only in-place callable names")
        );
        assert!(
            error_for(registry_with(&entry.replace("&mut", "&")))
                .contains("wrong instance receiver")
        );
        assert!(
            error_for(registry_with(
                &entry.replace("state:in_place", "state:read_only")
            ))
            .contains("only in-place callable names")
        );
    }
}

#[test]
fn receiver_capability_is_independent_of_business_type() {
    for target in ["BioStructure", "Protein", "OtherValue"] {
        let entry = format!(
            r#"{{
            semantic_id:"{target}.convert", item:callable, owner:type_,
            rust:crate::{target}::convert, python:"convert", javascript:"convert",
            feature:"cap-bio", kind:instance, receiver:owned, parameters:[], output:crate::{target}, error:none,
            state:value_returning, operation:none, signature:fn(crate::{target})->crate::{target}
        }}"#
        );
        assert!(expand_binding_contract(registry_with(&entry)).is_ok());
        for invalid in [
            entry.replace("receiver:owned", "receiver:shared"),
            entry.replace("receiver:owned", "receiver:mutable"),
            entry.replace("receiver:owned,", ""),
            entry.replace("receiver:owned", "receiver:unknown"),
            entry.replace("kind:instance", "kind:static_"),
            entry.replace("state:value_returning", "state:read_only"),
            entry.replace("state:value_returning", "state:in_place"),
            entry.replace("operation:none", "operation:\"convert\""),
            entry.replace("signature:fn(", "signature:fn(&"),
        ] {
            assert!(
                expand_binding_contract(registry_with(&invalid)).is_err(),
                "{invalid}"
            );
        }
        let borrowed = entry
            .replace("receiver:owned", "receiver:shared")
            .replace("signature:fn(", "signature:fn(&");
        assert!(expand_binding_contract(registry_with(&borrowed)).is_ok());
    }
    assert_eq!(OWNED_PROTEIN_BINDINGS.len(), 1);
}

#[test]
fn value_and_inplace_same_verb_keep_distinct_canonical_javascript_names() {
    let entries = r#"{
        semantic_id:"Molecule.sanitize",item:callable,owner:molecule,
        rust:crate::Molecule::sanitize,python:"sanitize",javascript:"sanitize",
        feature:"cap-sanitize",kind:instance,parameters:[],output:crate::Molecule,
        error:crate::OperationError,state:value_returning,operation:"sanitize",
        signature:fn(&crate::Molecule)->Result<crate::Molecule,crate::OperationError>
    },{
        semantic_id:"Molecule.sanitize_",item:callable,owner:molecule,
        rust:crate::Molecule::sanitize_,python:"sanitize_",javascript:"sanitize_",
        feature:"cap-sanitize",kind:instance,parameters:[],output:(),
        error:crate::OperationError,state:in_place,operation:"sanitize_",
        signature:fn(&mut crate::Molecule)->Result<(),crate::OperationError>
    }"#;
    expand_binding_contract(registry_with(entries)).expect("distinct paired projections");
    let ambiguous = entries.replace(r#"javascript:"sanitize_""#, r#"javascript:"sanitize""#);
    assert!(
        error_for(registry_with(&ambiguous))
            .contains("JavaScript callable name must be `sanitize_`")
    );
    let wrong = entries.replace(
        r#"javascript:"sanitize_""#,
        r#"javascript:"sanitizeInPlace""#,
    );
    assert!(
        error_for(registry_with(&wrong)).contains("JavaScript callable name must be `sanitize_`")
    );
}
