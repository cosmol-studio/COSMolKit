//! Private parser and generator for the canonical binding-contract registry.
//!
//! This module emits metadata and compile-time type assertions only. It does
//! not generate public API methods, runtime behavior, or domain algorithms.

use std::collections::HashSet;

use quote::{ToTokens, format_ident, quote};
use syn::{
    Attribute, Expr, Ident, LitStr, Meta, Path, ReturnType, Token, Type, TypeFnPtr, Visibility,
    braced, bracketed,
    ext::IdentExt,
    parse::{Parse, ParseStream},
    punctuated::Punctuated,
    visit_mut::{self, VisitMut},
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ItemClass {
    Callable,
    Type,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Owner {
    Molecule,
    Module,
    Type,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Exposure {
    Registered,
    Public,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum CallableKind {
    Instance,
    Static,
    Module,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum StateModel {
    ValueReturning,
    InPlace,
    ReadOnly,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum TypeRole {
    Value,
    Parameter,
    Result,
    Error,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Support {
    Unsupported,
    PreservedOnly,
    Experimental,
    Supported,
    SupportedWithRdkitParity,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Parity {
    NotApplicable,
    RequiredWhenSupported,
    RequiredNow,
}

#[derive(Clone)]
enum ParameterDefault {
    Required,
    Value(Expr),
}

#[derive(Clone)]
struct BindingParameter {
    name: Ident,
    ty: Type,
    default: ParameterDefault,
}

impl Parse for BindingParameter {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        let content;
        braced!(content in input);
        let mut name = None;
        let mut ty = None;
        let mut default = None;
        while !content.is_empty() {
            let key = content.call(Ident::parse_any)?;
            content.parse::<Token![:]>()?;
            match key.to_string().as_str() {
                "name" => set_once(&mut name, content.parse()?, &key)?,
                "type" => set_once(&mut ty, content.parse()?, &key)?,
                "default" => {
                    let value = if content.peek(Ident) {
                        let fork = content.fork();
                        let candidate: Ident = fork.parse()?;
                        if candidate == "required" && (fork.is_empty() || fork.peek(Token![,])) {
                            content.parse::<Ident>()?;
                            ParameterDefault::Required
                        } else {
                            ParameterDefault::Value(content.parse()?)
                        }
                    } else {
                        ParameterDefault::Value(content.parse()?)
                    };
                    set_once(&mut default, value, &key)?;
                }
                other => return Err(unknown_field(&key, "parameter", other)),
            }
            consume_comma(&content)?;
        }
        Ok(Self {
            name: required(name, "parameter.name")?,
            ty: required(ty, "parameter.type")?,
            default: required(default, "parameter.default")?,
        })
    }
}

#[derive(Clone)]
enum ErrorType {
    None,
    Typed(Type),
}
#[derive(Clone)]
enum OperationLink {
    None,
    Id(LitStr),
}

struct BindingRegistry {
    visibility: Visibility,
    name: Ident,
    entries: Vec<BindingEntry>,
}

struct BindingEntry {
    cfg_attrs: Vec<Attribute>,
    semantic_id: LitStr,
    item: ItemClass,
    owner: Owner,
    rust: Path,
    python: LitStr,
    javascript: LitStr,
    feature: LitStr,
    exposure: Exposure,
    support: Support,
    parity: Parity,
    callable: Option<CallablePayload>,
    type_role: Option<TypeRole>,
}

struct CallablePayload {
    kind: CallableKind,
    parameters: Vec<BindingParameter>,
    output: Type,
    error: ErrorType,
    state: StateModel,
    operation: OperationLink,
    signature: TypeFnPtr,
}

#[derive(Default)]
struct BindingEntryDraft {
    semantic_id: Option<LitStr>,
    item: Option<Ident>,
    owner: Option<Ident>,
    rust: Option<Path>,
    python: Option<LitStr>,
    javascript: Option<LitStr>,
    feature: Option<LitStr>,
    exposure: Option<Ident>,
    support: Option<Ident>,
    parity: Option<Ident>,
    kind: Option<Ident>,
    parameters: Option<Vec<BindingParameter>>,
    output: Option<Type>,
    error: Option<ErrorType>,
    state: Option<Ident>,
    operation: Option<OperationLink>,
    signature: Option<TypeFnPtr>,
    role: Option<Ident>,
}

impl Parse for BindingRegistry {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        let visibility = input.parse()?;
        input.parse::<Token![static]>()?;
        let name = input.parse()?;
        input.parse::<Token![=]>()?;
        let content;
        bracketed!(content in input);
        let mut entries = Vec::new();
        while !content.is_empty() {
            let cfg_attrs = content.call(Attribute::parse_outer)?;
            let entry_content;
            braced!(entry_content in content);
            entries.push(parse_binding_entry(cfg_attrs, &entry_content)?);
            consume_comma(&content)?;
        }
        if input.peek(Token![;]) {
            input.parse::<Token![;]>()?;
        }
        if !input.is_empty() {
            return Err(input.error("unexpected tokens after binding registry"));
        }
        validate_registry(&entries)?;
        Ok(Self {
            visibility,
            name,
            entries,
        })
    }
}

fn parse_binding_entry(
    cfg_attrs: Vec<Attribute>,
    input: ParseStream<'_>,
) -> syn::Result<BindingEntry> {
    let mut draft = BindingEntryDraft::default();
    while !input.is_empty() {
        let key = input.call(Ident::parse_any)?;
        input.parse::<Token![:]>()?;
        match key.to_string().as_str() {
            "semantic_id" => set_once(&mut draft.semantic_id, input.parse()?, &key)?,
            "item" => set_once(&mut draft.item, input.call(Ident::parse_any)?, &key)?,
            "owner" => set_once(&mut draft.owner, input.parse()?, &key)?,
            "rust" => set_once(&mut draft.rust, input.parse()?, &key)?,
            "python" => set_once(&mut draft.python, input.parse()?, &key)?,
            "javascript" => set_once(&mut draft.javascript, input.parse()?, &key)?,
            "feature" => set_once(&mut draft.feature, input.parse()?, &key)?,
            "exposure" => set_once(&mut draft.exposure, input.parse()?, &key)?,
            "support" => set_once(&mut draft.support, input.parse()?, &key)?,
            "parity" => set_once(&mut draft.parity, input.parse()?, &key)?,
            "kind" => set_once(&mut draft.kind, input.parse()?, &key)?,
            "parameters" => {
                let values;
                bracketed!(values in input);
                let parameters =
                    Punctuated::<BindingParameter, Token![,]>::parse_terminated(&values)?
                        .into_iter()
                        .collect();
                set_once(&mut draft.parameters, parameters, &key)?;
            }
            "output" => set_once(&mut draft.output, input.parse()?, &key)?,
            "error" => {
                let value = parse_none_or_type(input)?;
                set_once(&mut draft.error, value, &key)?;
            }
            "state" => set_once(&mut draft.state, input.parse()?, &key)?,
            "operation" => {
                let value = parse_none_or_string(input)?;
                set_once(&mut draft.operation, value, &key)?;
            }
            "signature" => set_once(&mut draft.signature, input.parse()?, &key)?,
            "role" => set_once(&mut draft.role, input.parse()?, &key)?,
            other => return Err(unknown_field(&key, "binding entry", other)),
        }
        consume_comma(input)?;
    }
    finish_entry(cfg_attrs, draft)
}

fn finish_entry(cfg_attrs: Vec<Attribute>, draft: BindingEntryDraft) -> syn::Result<BindingEntry> {
    let semantic_id = required(draft.semantic_id, "semantic_id")?;
    let item_ident = required(draft.item, "item")?;
    let item = parse_item(&item_ident)?;
    let owner_ident = required(draft.owner, "owner")?;
    let owner = parse_owner(&owner_ident)?;
    let rust = required(draft.rust, "rust")?;
    let python = required(draft.python, "python")?;
    let javascript = required(draft.javascript, "javascript")?;
    let feature = required(draft.feature, "feature")?;
    let exposure_ident = required(draft.exposure, "exposure")?;
    let exposure = parse_exposure(&exposure_ident)?;
    let support_ident = required(draft.support, "support")?;
    let support = parse_support(&support_ident)?;
    let parity_ident = required(draft.parity, "parity")?;
    let parity = parse_parity(&parity_ident)?;
    require_nonempty(&semantic_id, "semantic_id")?;
    require_nonempty(&python, "python")?;
    require_nonempty(&javascript, "javascript")?;
    require_nonempty(&feature, "feature")?;
    validate_cfg(&cfg_attrs, &feature)?;
    validate_support_parity(support, parity, &parity_ident)?;

    let (callable, type_role) = match item {
        ItemClass::Callable => {
            if draft.role.is_some() {
                return Err(syn::Error::new_spanned(
                    draft.role,
                    "callable binding entry cannot declare `role`",
                ));
            }
            let kind_ident = required(draft.kind, "kind")?;
            let state_ident = required(draft.state, "state")?;
            let payload = CallablePayload {
                kind: parse_kind(&kind_ident)?,
                parameters: required(draft.parameters, "parameters")?,
                output: required(draft.output, "output")?,
                error: required(draft.error, "error")?,
                state: parse_state(&state_ident)?,
                operation: required(draft.operation, "operation")?,
                signature: required(draft.signature, "signature")?,
            };
            (Some(payload), None)
        }
        ItemClass::Type => {
            reject_present(draft.kind, "kind", "type", &item_ident)?;
            reject_present(draft.parameters, "parameters", "type", &item_ident)?;
            reject_present(draft.output, "output", "type", &item_ident)?;
            reject_present(draft.error, "error", "type", &item_ident)?;
            reject_present(draft.state, "state", "type", &item_ident)?;
            reject_present(draft.operation, "operation", "type", &item_ident)?;
            reject_present(draft.signature, "signature", "type", &item_ident)?;
            if owner != Owner::Type {
                return Err(syn::Error::new_spanned(
                    owner_ident,
                    "type binding entry requires `owner: type_`",
                ));
            }
            let role_ident = required(draft.role, "role")?;
            (None, Some(parse_type_role(&role_ident)?))
        }
    };
    Ok(BindingEntry {
        cfg_attrs,
        semantic_id,
        item,
        owner,
        rust,
        python,
        javascript,
        feature,
        exposure,
        support,
        parity,
        callable,
        type_role,
    })
}

fn validate_registry(entries: &[BindingEntry]) -> syn::Result<()> {
    let mut ids = HashSet::new();
    let mut rust = HashSet::new();
    for entry in entries {
        let owner = owner_key(entry.owner);
        insert_unique(
            &mut ids,
            entry.semantic_id.value(),
            &entry.semantic_id,
            "duplicate binding semantic_id",
        )?;
        insert_unique(
            &mut rust,
            format!("{owner}:{}", tokens(&entry.rust)),
            &entry.rust,
            "duplicate Rust binding projection in one owner scope",
        )?;
    }
    for entry in entries {
        if let Some(payload) = &entry.callable {
            validate_callable(
                &entry.semantic_id,
                entry.owner,
                &entry.rust,
                &entry.python,
                &entry.javascript,
                payload,
            )?;
        } else {
            validate_type_names(
                &entry.semantic_id,
                &entry.rust,
                &entry.python,
                &entry.javascript,
            )?;
        }
    }
    let mut python = HashSet::new();
    let mut javascript = HashSet::new();
    for entry in entries {
        let owner = projection_owner_key(entry);
        insert_unique(
            &mut python,
            format!("{owner}:{}", entry.python.value()),
            &entry.python,
            "duplicate Python binding projection in one owner scope",
        )?;
        insert_unique(
            &mut javascript,
            format!("{owner}:{}", entry.javascript.value()),
            &entry.javascript,
            "duplicate JavaScript binding projection in one owner scope",
        )?;
    }
    Ok(())
}

fn validate_callable(
    semantic_id: &LitStr,
    owner: Owner,
    rust: &Path,
    python: &LitStr,
    javascript: &LitStr,
    payload: &CallablePayload,
) -> syn::Result<()> {
    match (owner, payload.kind) {
        (Owner::Molecule, CallableKind::Instance | CallableKind::Static)
        | (Owner::Type, CallableKind::Instance | CallableKind::Static)
        | (Owner::Module, CallableKind::Module) => {}
        _ => {
            return Err(syn::Error::new_spanned(
                rust,
                "invalid callable owner/kind combination",
            ));
        }
    }
    if payload.state == StateModel::InPlace && payload.kind != CallableKind::Instance {
        return Err(syn::Error::new_spanned(
            rust,
            "in-place binding requires an instance receiver",
        ));
    }
    let rust_name = rust_last_name(rust)?;
    let logical_name = logical_name(semantic_id)?;
    if rust_name != logical_name {
        return Err(syn::Error::new_spanned(
            rust,
            "Rust callable name must equal the semantic_id final component",
        ));
    }
    for forbidden in ["calc_", "mol_to_", "mol_from_", "get_"] {
        if rust_name.starts_with(forbidden) {
            return Err(syn::Error::new_spanned(
                rust,
                format!("non-canonical public prefix `{forbidden}` is forbidden"),
            ));
        }
    }
    if python.value() != rust_name {
        return Err(syn::Error::new_spanned(
            python,
            "Python callable name must equal the canonical Rust name",
        ));
    }
    let requires_trailing_underscore =
        owner == Owner::Molecule && payload.state == StateModel::InPlace;
    if rust_name.ends_with('_') != requires_trailing_underscore {
        return Err(syn::Error::new_spanned(
            rust,
            "only in-place callable names may end in `_`",
        ));
    }
    let canonical_js = snake_to_camel(rust_name.trim_end_matches('_'));
    if javascript.value() != canonical_js {
        return Err(syn::Error::new_spanned(
            javascript,
            format!("JavaScript callable name must be `{canonical_js}`"),
        ));
    }
    let mut names = HashSet::new();
    let mut saw_default = false;
    for parameter in &payload.parameters {
        if !names.insert(parameter.name.to_string()) {
            return Err(syn::Error::new_spanned(
                &parameter.name,
                "duplicate binding parameter name",
            ));
        }
        match parameter.default {
            ParameterDefault::Required if saw_default => {
                return Err(syn::Error::new_spanned(
                    &parameter.name,
                    "required parameter cannot follow a defaulted parameter",
                ));
            }
            ParameterDefault::Required => {}
            ParameterDefault::Value(_) => saw_default = true,
        }
    }
    if let OperationLink::Id(operation) = &payload.operation {
        require_nonempty(operation, "operation")?;
        if operation.value() != logical_name {
            return Err(syn::Error::new_spanned(
                operation,
                "operation id must equal the callable semantic_id final component",
            ));
        }
    }
    validate_signature(owner, rust, payload)
}

fn validate_signature(owner: Owner, rust: &Path, payload: &CallablePayload) -> syn::Result<()> {
    if payload.signature.unsafety.is_some()
        || payload.signature.abi.is_some()
        || payload.signature.variadic.is_some()
    {
        return Err(syn::Error::new_spanned(
            &payload.signature,
            "binding signature must be a safe non-variadic Rust fn with the default ABI",
        ));
    }
    let inputs: Vec<&Type> = payload.signature.inputs.iter().map(|arg| &arg.ty).collect();
    let receiver_count = usize::from(payload.kind == CallableKind::Instance);
    if inputs.len() != receiver_count + payload.parameters.len() {
        return Err(syn::Error::new_spanned(
            &payload.signature,
            "binding signature argument count disagrees with receiver and parameters",
        ));
    }
    if payload.kind == CallableKind::Instance {
        let expected: Type = match owner {
            Owner::Molecule if payload.state == StateModel::InPlace => {
                syn::parse_quote!(&mut crate::Molecule)
            }
            Owner::Molecule => syn::parse_quote!(&crate::Molecule),
            Owner::Type => {
                let mut receiver_path = rust.clone();
                receiver_path.segments.pop();
                receiver_path.segments.pop_punct();
                if receiver_path.segments.is_empty() {
                    return Err(syn::Error::new_spanned(
                        rust,
                        "type-owned callable path must include its receiver type",
                    ));
                }
                match payload.state {
                    StateModel::InPlace => syn::parse_quote!(&mut #receiver_path),
                    StateModel::ReadOnly => syn::parse_quote!(&#receiver_path),
                    StateModel::ValueReturning => syn::parse_quote!(#receiver_path),
                }
            }
            Owner::Module => {
                return Err(syn::Error::new_spanned(
                    rust,
                    "module callable cannot have an instance receiver",
                ));
            }
        };
        if signature_metadata_tokens(inputs[0]) != signature_metadata_tokens(&expected) {
            return Err(syn::Error::new_spanned(
                inputs[0],
                "binding signature has the wrong instance receiver",
            ));
        }
    } else if owner == Owner::Molecule && payload.state == StateModel::InPlace {
        return Err(syn::Error::new_spanned(
            &payload.signature,
            "static Molecule callable cannot be in-place",
        ));
    }
    for (actual, parameter) in inputs[receiver_count..].iter().zip(&payload.parameters) {
        if signature_metadata_tokens(actual) != signature_metadata_tokens(&parameter.ty) {
            return Err(syn::Error::new_spanned(
                *actual,
                format!(
                    "signature type disagrees with parameter `{}`",
                    parameter.name
                ),
            ));
        }
    }
    let actual_output = match &payload.signature.output {
        ReturnType::Default => syn::parse_quote!(()),
        ReturnType::Type(_, ty) => ty.clone(),
    };
    let expected_output: Type = match &payload.error {
        ErrorType::None => payload.output.clone(),
        ErrorType::Typed(error) => {
            let output = &payload.output;
            syn::parse_quote!(Result<#output, #error>)
        }
    };
    if signature_metadata_tokens(&actual_output) != signature_metadata_tokens(&expected_output) {
        return Err(syn::Error::new_spanned(
            &payload.signature.output,
            "binding signature return type disagrees with output/error metadata",
        ));
    }
    Ok(())
}

/// Binding metadata records public value types, while an explicit bare-fn
/// binder records which borrowed input owns a borrowed output. Erase only
/// reference lifetime spellings for metadata agreement; the generated const
/// assertion retains the complete declared signature and asks rustc to verify
/// the exact higher-ranked relationship against the public item.
fn signature_metadata_tokens(ty: &Type) -> String {
    struct ElideReferenceLifetimes;

    impl VisitMut for ElideReferenceLifetimes {
        fn visit_type_reference_mut(&mut self, reference: &mut syn::TypeReference) {
            reference.lifetime = None;
            visit_mut::visit_type_reference_mut(self, reference);
        }
    }

    let mut normalized = ty.clone();
    ElideReferenceLifetimes.visit_type_mut(&mut normalized);
    tokens(&normalized)
}

fn validate_type_names(
    semantic_id: &LitStr,
    rust: &Path,
    python: &LitStr,
    javascript: &LitStr,
) -> syn::Result<()> {
    let rust_name = rust_last_name(rust)?;
    if logical_name(semantic_id)? != rust_name {
        return Err(syn::Error::new_spanned(
            rust,
            "Rust type name must equal the semantic_id final component",
        ));
    }
    if python.value() != rust_name || javascript.value() != rust_name {
        return Err(syn::Error::new_spanned(
            rust,
            "type projections must preserve the canonical Rust type name",
        ));
    }
    Ok(())
}

fn validate_cfg(attrs: &[Attribute], feature: &LitStr) -> syn::Result<()> {
    if attrs.len() > 1 {
        return Err(syn::Error::new_spanned(
            &attrs[1],
            "binding entry accepts at most one cfg attribute",
        ));
    }
    let Some(attribute) = attrs.first() else {
        return Ok(());
    };
    if !attribute.path().is_ident("cfg") {
        return Err(syn::Error::new_spanned(
            attribute,
            "binding entry accepts only cfg(feature = \"...\")",
        ));
    }
    let Meta::List(list) = &attribute.meta else {
        return Err(syn::Error::new_spanned(
            attribute,
            "binding cfg must be cfg(feature = \"...\")",
        ));
    };
    let values = list.parse_args_with(Punctuated::<Meta, Token![,]>::parse_terminated)?;
    let Some(Meta::NameValue(value)) = values.first() else {
        return Err(syn::Error::new_spanned(
            list,
            "compound or non-feature binding cfg is forbidden",
        ));
    };
    if values.len() != 1 || !value.path.is_ident("feature") {
        return Err(syn::Error::new_spanned(
            list,
            "compound or non-feature binding cfg is forbidden",
        ));
    }
    let Expr::Lit(expression) = &value.value else {
        return Err(syn::Error::new_spanned(
            &value.value,
            "binding cfg feature must be a string literal",
        ));
    };
    let syn::Lit::Str(cfg_feature) = &expression.lit else {
        return Err(syn::Error::new_spanned(
            &expression.lit,
            "binding cfg feature must be a string literal",
        ));
    };
    if cfg_feature.value() != feature.value() {
        return Err(syn::Error::new_spanned(
            cfg_feature,
            "binding cfg feature disagrees with the feature field",
        ));
    }
    Ok(())
}

fn validate_support_parity(support: Support, parity: Parity, span: &Ident) -> syn::Result<()> {
    if (parity == Parity::RequiredNow) != (support == Support::SupportedWithRdkitParity) {
        Err(syn::Error::new_spanned(
            span,
            "required_now and supported_with_rdkit_parity must be declared together",
        ))
    } else {
        Ok(())
    }
}

pub(crate) fn expand_binding_contract(
    input: proc_macro2::TokenStream,
) -> syn::Result<proc_macro2::TokenStream> {
    expand_registry(syn::parse2::<BindingRegistry>(input)?)
}

fn expand_registry(registry: BindingRegistry) -> syn::Result<proc_macro2::TokenStream> {
    let BindingRegistry {
        visibility,
        name,
        entries,
    } = registry;
    let mut values = Vec::with_capacity(entries.len());
    let mut assertions = Vec::new();
    for (index, entry) in entries.iter().enumerate() {
        let cfg = &entry.cfg_attrs;
        let semantic_id = &entry.semantic_id;
        let item = item_tokens(entry.item);
        let owner = owner_tokens(entry.owner);
        let rust = &entry.rust;
        let python = &entry.python;
        let javascript = &entry.javascript;
        let feature = &entry.feature;
        let exposure = exposure_tokens(entry.exposure);
        let support = support_tokens(entry.support);
        let parity = parity_tokens(entry.parity);
        let (callable, role) = match (&entry.callable, entry.type_role) {
            (Some(payload), None) => {
                let kind = kind_tokens(payload.kind);
                let state = state_tokens(payload.state);
                let parameters = payload.parameters.iter().map(|parameter| {
                    let name = parameter.name.to_string();
                    let ty = &parameter.ty;
                    let default = match &parameter.default {
                        ParameterDefault::Required => quote!(crate::BindingDefault::Required),
                        ParameterDefault::Value(value) => quote!(crate::BindingDefault::Value(stringify!(#value))),
                    };
                    quote!(crate::BindingParameterContract { name: #name, type_name: stringify!(#ty), default: #default })
                });
                let output = &payload.output;
                let error = match &payload.error {
                    ErrorType::None => quote!(None),
                    ErrorType::Typed(error) => quote!(Some(stringify!(#error))),
                };
                let operation = match &payload.operation {
                    OperationLink::None => quote!(None),
                    OperationLink::Id(id) => quote!(Some(#id)),
                };
                (
                    quote!(Some(crate::BindingCallableContract {
                        kind: #kind, parameters: &[#(#parameters),*], output_type: stringify!(#output),
                        error_type: #error, state_model: #state, operation_semantic_id: #operation,
                    })),
                    quote!(None),
                )
            }
            (None, Some(role)) => {
                let role = type_role_tokens(role);
                (quote!(None), quote!(Some(#role)))
            }
            _ => {
                return Err(syn::Error::new_spanned(
                    semantic_id,
                    "validated binding entry lost its payload",
                ));
            }
        };
        values.push(quote! {
            #(#cfg)*
            crate::BindingContractEntry {
                semantic_id: #semantic_id, item: #item, owner: #owner,
                rust_path: stringify!(#rust), python_name: #python, javascript_name: #javascript,
                feature: #feature, exposure: #exposure, support: #support, parity: #parity,
                callable: #callable, type_role: #role,
            }
        });
        if entry.exposure == Exposure::Public {
            let assertion = format_ident!("__BINDING_ASSERT_{}_{}", name, index);
            if let Some(payload) = &entry.callable {
                let signature = &payload.signature;
                assertions.push(quote! { #(#cfg)* const #assertion: #signature = #rust; });
            } else {
                assertions.push(quote! { #(#cfg)* const #assertion: fn() = || {
                    fn assert_public_type<T>() {} assert_public_type::<#rust>();
                }; });
            }
        }
    }
    Ok(quote! {
        #visibility static #name: &[crate::BindingContractEntry] = &[#(#values),*];
        #(#assertions)*
    })
}

fn parse_none_or_type(input: ParseStream<'_>) -> syn::Result<ErrorType> {
    if input.peek(Ident) {
        let fork = input.fork();
        let candidate: Ident = fork.parse()?;
        if candidate == "none" && (fork.is_empty() || fork.peek(Token![,])) {
            input.parse::<Ident>()?;
            return Ok(ErrorType::None);
        }
    }
    Ok(ErrorType::Typed(input.parse()?))
}

fn parse_none_or_string(input: ParseStream<'_>) -> syn::Result<OperationLink> {
    if input.peek(Ident) {
        let fork = input.fork();
        let candidate: Ident = fork.parse()?;
        if candidate == "none" && (fork.is_empty() || fork.peek(Token![,])) {
            input.parse::<Ident>()?;
            return Ok(OperationLink::None);
        }
    }
    Ok(OperationLink::Id(input.parse()?))
}

fn parse_item(v: &Ident) -> syn::Result<ItemClass> {
    parse_enum(
        v,
        &["callable", "type"],
        |n| {
            if n == "callable" {
                ItemClass::Callable
            } else {
                ItemClass::Type
            }
        },
        "item must be `callable` or `type`",
    )
}
fn parse_owner(v: &Ident) -> syn::Result<Owner> {
    parse_enum(
        v,
        &["molecule", "module", "type_"],
        |n| match n {
            "molecule" => Owner::Molecule,
            "module" => Owner::Module,
            _ => Owner::Type,
        },
        "owner must be `molecule`, `module`, or `type_`",
    )
}
fn parse_exposure(v: &Ident) -> syn::Result<Exposure> {
    parse_enum(
        v,
        &["registered", "public"],
        |n| {
            if n == "registered" {
                Exposure::Registered
            } else {
                Exposure::Public
            }
        },
        "exposure must be `registered` or `public`",
    )
}
fn parse_kind(v: &Ident) -> syn::Result<CallableKind> {
    parse_enum(
        v,
        &["instance", "static_", "module"],
        |n| match n {
            "instance" => CallableKind::Instance,
            "static_" => CallableKind::Static,
            _ => CallableKind::Module,
        },
        "kind must be `instance`, `static_`, or `module`",
    )
}
fn parse_state(v: &Ident) -> syn::Result<StateModel> {
    parse_enum(
        v,
        &["value_returning", "in_place", "read_only"],
        |n| match n {
            "value_returning" => StateModel::ValueReturning,
            "in_place" => StateModel::InPlace,
            _ => StateModel::ReadOnly,
        },
        "state must be `value_returning`, `in_place`, or `read_only`",
    )
}
fn parse_type_role(v: &Ident) -> syn::Result<TypeRole> {
    parse_enum(
        v,
        &["value", "parameter", "result", "error"],
        |n| match n {
            "value" => TypeRole::Value,
            "parameter" => TypeRole::Parameter,
            "result" => TypeRole::Result,
            _ => TypeRole::Error,
        },
        "role must be `value`, `parameter`, `result`, or `error`",
    )
}
fn parse_support(v: &Ident) -> syn::Result<Support> {
    parse_enum(
        v,
        &[
            "unsupported",
            "preserved_only",
            "experimental",
            "supported",
            "supported_with_rdkit_parity",
        ],
        |n| match n {
            "unsupported" => Support::Unsupported,
            "preserved_only" => Support::PreservedOnly,
            "experimental" => Support::Experimental,
            "supported" => Support::Supported,
            _ => Support::SupportedWithRdkitParity,
        },
        "unsupported binding support status",
    )
}
fn parse_parity(v: &Ident) -> syn::Result<Parity> {
    parse_enum(
        v,
        &["not_applicable", "required_when_supported", "required_now"],
        |n| match n {
            "not_applicable" => Parity::NotApplicable,
            "required_when_supported" => Parity::RequiredWhenSupported,
            _ => Parity::RequiredNow,
        },
        "unsupported binding parity policy",
    )
}

fn parse_enum<T>(
    value: &Ident,
    allowed: &[&str],
    convert: impl FnOnce(&str) -> T,
    message: &str,
) -> syn::Result<T> {
    let name = value.to_string();
    if allowed.contains(&name.as_str()) {
        Ok(convert(&name))
    } else {
        Err(syn::Error::new_spanned(value, message))
    }
}

fn item_tokens(v: ItemClass) -> proc_macro2::TokenStream {
    match v {
        ItemClass::Callable => quote!(crate::BindingItem::Callable),
        ItemClass::Type => quote!(crate::BindingItem::Type),
    }
}
fn owner_tokens(v: Owner) -> proc_macro2::TokenStream {
    match v {
        Owner::Molecule => quote!(crate::BindingOwner::Molecule),
        Owner::Module => quote!(crate::BindingOwner::Module),
        Owner::Type => quote!(crate::BindingOwner::Type),
    }
}
fn exposure_tokens(v: Exposure) -> proc_macro2::TokenStream {
    match v {
        Exposure::Registered => quote!(crate::BindingExposure::Registered),
        Exposure::Public => quote!(crate::BindingExposure::Public),
    }
}
fn kind_tokens(v: CallableKind) -> proc_macro2::TokenStream {
    match v {
        CallableKind::Instance => quote!(crate::BindingKind::Instance),
        CallableKind::Static => quote!(crate::BindingKind::Static),
        CallableKind::Module => quote!(crate::BindingKind::Module),
    }
}
fn state_tokens(v: StateModel) -> proc_macro2::TokenStream {
    match v {
        StateModel::ValueReturning => quote!(crate::StateModel::ValueReturning),
        StateModel::InPlace => quote!(crate::StateModel::InPlace),
        StateModel::ReadOnly => quote!(crate::StateModel::ReadOnly),
    }
}
fn type_role_tokens(v: TypeRole) -> proc_macro2::TokenStream {
    match v {
        TypeRole::Value => quote!(crate::BindingTypeRole::Value),
        TypeRole::Parameter => quote!(crate::BindingTypeRole::Parameter),
        TypeRole::Result => quote!(crate::BindingTypeRole::Result),
        TypeRole::Error => quote!(crate::BindingTypeRole::Error),
    }
}
fn support_tokens(v: Support) -> proc_macro2::TokenStream {
    match v {
        Support::Unsupported => quote!(crate::BindingSupport::Unsupported),
        Support::PreservedOnly => quote!(crate::BindingSupport::PreservedOnly),
        Support::Experimental => quote!(crate::BindingSupport::Experimental),
        Support::Supported => quote!(crate::BindingSupport::Supported),
        Support::SupportedWithRdkitParity => {
            quote!(crate::BindingSupport::SupportedWithRdkitParity)
        }
    }
}
fn parity_tokens(v: Parity) -> proc_macro2::TokenStream {
    match v {
        Parity::NotApplicable => quote!(crate::BindingParity::NotApplicable),
        Parity::RequiredWhenSupported => quote!(crate::BindingParity::RequiredWhenSupported),
        Parity::RequiredNow => quote!(crate::BindingParity::RequiredNow),
    }
}

fn set_once<T>(slot: &mut Option<T>, value: T, key: &Ident) -> syn::Result<()> {
    if slot.replace(value).is_some() {
        Err(syn::Error::new_spanned(
            key,
            format!("duplicate binding field `{key}`"),
        ))
    } else {
        Ok(())
    }
}
fn required<T>(value: Option<T>, field: &str) -> syn::Result<T> {
    value.ok_or_else(|| {
        syn::Error::new(
            proc_macro2::Span::call_site(),
            format!("binding entry is missing `{field}`"),
        )
    })
}
fn reject_present<T>(
    value: Option<T>,
    field: &str,
    item: &str,
    span: &impl ToTokens,
) -> syn::Result<()> {
    if value.is_some() {
        Err(syn::Error::new_spanned(
            span,
            format!("{item} binding entry cannot declare `{field}`"),
        ))
    } else {
        Ok(())
    }
}
fn require_nonempty(value: &LitStr, field: &str) -> syn::Result<()> {
    if value.value().is_empty() {
        Err(syn::Error::new_spanned(
            value,
            format!("binding `{field}` cannot be empty"),
        ))
    } else {
        Ok(())
    }
}
fn insert_unique<T: ToTokens>(
    set: &mut HashSet<String>,
    key: String,
    span: &T,
    message: &str,
) -> syn::Result<()> {
    if set.insert(key) {
        Ok(())
    } else {
        Err(syn::Error::new_spanned(span, message))
    }
}
fn consume_comma(input: ParseStream<'_>) -> syn::Result<()> {
    if input.peek(Token![,]) {
        input.parse::<Token![,]>()?;
    }
    Ok(())
}
fn unknown_field(key: &Ident, context: &str, field: &str) -> syn::Error {
    syn::Error::new_spanned(key, format!("unknown {context} field `{field}`"))
}
fn rust_last_name(path: &Path) -> syn::Result<String> {
    path.segments
        .last()
        .map(|s| s.ident.to_string())
        .ok_or_else(|| syn::Error::new_spanned(path, "Rust binding path cannot be empty"))
}
fn logical_name(id: &LitStr) -> syn::Result<String> {
    id.value()
        .rsplit('.')
        .next()
        .filter(|n| !n.is_empty())
        .map(ToOwned::to_owned)
        .ok_or_else(|| syn::Error::new_spanned(id, "semantic_id has no logical name"))
}
fn snake_to_camel(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    let mut upper = false;
    for c in value.chars() {
        if c == '_' {
            upper = true;
        } else if upper {
            out.extend(c.to_uppercase());
            upper = false;
        } else {
            out.push(c);
        }
    }
    out
}
fn owner_key(owner: Owner) -> &'static str {
    match owner {
        Owner::Molecule => "molecule",
        Owner::Module => "module",
        Owner::Type => "type",
    }
}
fn projection_owner_key(entry: &BindingEntry) -> String {
    match entry.owner {
        Owner::Molecule | Owner::Module => owner_key(entry.owner).to_owned(),
        Owner::Type => {
            if entry.callable.is_none() {
                // Type declarations share the language-level exported type
                // namespace; distinct Rust paths do not create distinct
                // Python or JavaScript declaration scopes.
                return owner_key(Owner::Type).to_owned();
            }
            let mut owner = entry.rust.clone();
            owner.segments.pop();
            owner.segments.pop_punct();
            format!("type:{}", tokens(&owner))
        }
    }
}
fn tokens(value: &impl ToTokens) -> String {
    value.to_token_stream().to_string()
}
