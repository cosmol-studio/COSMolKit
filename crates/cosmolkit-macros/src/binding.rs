//! Private parser and generator for the canonical binding-contract registry.
//!
//! This module emits metadata and compile-time type assertions only. It does
//! not generate public API methods, runtime behavior, or domain algorithms.

use std::collections::{HashMap, HashSet};

use crate::status::FunctionStatus;

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
enum CallableKind {
    Instance,
    Static,
    Module,
    Constructor,
}
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Receiver {
    Shared,
    Mutable,
    Owned,
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
enum PythonProperty {
    Getter,
    Setter,
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

struct BindingProperty {
    name: Ident,
    rust: Path,
    signature: TypeFnPtr,
}
impl Parse for BindingProperty {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        let content;
        braced!(content in input);
        let mut name = None;
        let mut rust = None;
        let mut signature = None;
        while !content.is_empty() {
            let key = content.call(Ident::parse_any)?;
            content.parse::<Token![:]>()?;
            match key.to_string().as_str() {
                "name" => set_once(&mut name, content.parse()?, &key)?,
                "rust" => set_once(&mut rust, content.parse()?, &key)?,
                "signature" => set_once(&mut signature, content.parse()?, &key)?,
                other => return Err(unknown_field(&key, "read-only property", other)),
            }
            consume_comma(&content)?;
        }
        Ok(Self {
            name: required(name, "property.name")?,
            rust: required(rust, "property.rust")?,
            signature: required(signature, "property.signature")?,
        })
    }
}
struct BindingKeywordProjection {
    parameters: LitStr,
    target: LitStr,
}
impl Parse for BindingKeywordProjection {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        let content;
        braced!(content in input);
        let mut parameters = None;
        let mut target = None;
        while !content.is_empty() {
            let key = content.call(Ident::parse_any)?;
            content.parse::<Token![:]>()?;
            match key.to_string().as_str() {
                "parameters" => set_once(&mut parameters, content.parse()?, &key)?,
                "target" => set_once(&mut target, content.parse()?, &key)?,
                other => return Err(unknown_field(&key, "Python keyword projection", other)),
            }
            consume_comma(&content)?;
        }
        Ok(Self {
            parameters: required(parameters, "python_keywords.parameters")?,
            target: required(target, "python_keywords.target")?,
        })
    }
}

struct BindingEntry {
    cfg_attrs: Vec<Attribute>,
    semantic_id: LitStr,
    item: ItemClass,
    owner: Owner,
    rust: Path,
    python: LitStr,
    python_property: Option<PythonProperty>,
    javascript: LitStr,
    feature: LitStr,
    requires: Vec<LitStr>,
    status: Option<FunctionStatus>,
    callable: Option<CallablePayload>,
    type_role: Option<TypeRole>,
    properties: Vec<BindingProperty>,
    python_keywords: Option<BindingKeywordProjection>,
    python_adapters: Vec<PythonAdapter>,
}

/// Language-object ingress composes registered Rust APIs; it is not a fake
/// Rust callable accepting a Python object or a second chemistry registry.
struct PythonAdapter {
    name: Ident,
    targets: Vec<LitStr>,
}
impl Parse for PythonAdapter {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        let content;
        braced!(content in input);
        let mut name = None;
        let mut targets = None;
        while !content.is_empty() {
            let key = content.call(Ident::parse_any)?;
            content.parse::<Token![:]>()?;
            match key.to_string().as_str() {
                "name" => set_once(&mut name, content.parse()?, &key)?,
                "targets" => {
                    let values;
                    bracketed!(values in content);
                    set_once(
                        &mut targets,
                        Punctuated::<LitStr, Token![,]>::parse_terminated(&values)?
                            .into_iter()
                            .collect(),
                        &key,
                    )?;
                }
                other => return Err(unknown_field(&key, "Python adapter", other)),
            }
            consume_comma(&content)?;
        }
        Ok(Self {
            name: required(name, "adapter.name")?,
            targets: required(targets, "adapter.targets")?,
        })
    }
}

struct CallablePayload {
    kind: CallableKind,
    receiver: Option<Receiver>,
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
    python_property: Option<Ident>,
    javascript: Option<LitStr>,
    feature: Option<LitStr>,
    requires: Option<Vec<LitStr>>,
    status: Option<FunctionStatus>,
    kind: Option<Ident>,
    receiver: Option<Ident>,
    parameters: Option<Vec<BindingParameter>>,
    output: Option<Type>,
    error: Option<ErrorType>,
    state: Option<Ident>,
    operation: Option<OperationLink>,
    signature: Option<TypeFnPtr>,
    role: Option<Ident>,
    properties: Option<Vec<BindingProperty>>,
    python_keywords: Option<BindingKeywordProjection>,
    python_adapters: Option<Vec<PythonAdapter>>,
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
            "python_property" => set_once(&mut draft.python_property, input.parse()?, &key)?,
            "javascript" => set_once(&mut draft.javascript, input.parse()?, &key)?,
            "feature" => set_once(&mut draft.feature, input.parse()?, &key)?,
            "requires" => {
                let values;
                bracketed!(values in input);
                let required = Punctuated::<LitStr, Token![,]>::parse_terminated(&values)?
                    .into_iter()
                    .collect();
                set_once(&mut draft.requires, required, &key)?;
            }
            "status" => set_once(&mut draft.status, input.parse()?, &key)?,
            "kind" => set_once(&mut draft.kind, input.parse()?, &key)?,
            "receiver" => set_once(&mut draft.receiver, input.parse()?, &key)?,
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
            "properties" => {
                let values;
                bracketed!(values in input);
                set_once(
                    &mut draft.properties,
                    Punctuated::<BindingProperty, Token![,]>::parse_terminated(&values)?
                        .into_iter()
                        .collect(),
                    &key,
                )?;
            }
            "python_keywords" => set_once(&mut draft.python_keywords, input.parse()?, &key)?,
            "python_adapters" => {
                let values;
                bracketed!(values in input);
                set_once(
                    &mut draft.python_adapters,
                    Punctuated::<PythonAdapter, Token![,]>::parse_terminated(&values)?
                        .into_iter()
                        .collect(),
                    &key,
                )?;
            }
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
    let python_property = draft
        .python_property
        .map(|value| match value.to_string().as_str() {
            "getter" => Ok(PythonProperty::Getter),
            "setter" => Ok(PythonProperty::Setter),
            _ => Err(syn::Error::new_spanned(
                value,
                "python_property must be getter or setter",
            )),
        })
        .transpose()?;
    if python_property.is_some() && item != ItemClass::Callable {
        return Err(syn::Error::new_spanned(
            &python,
            "python_property requires a callable accessor",
        ));
    }
    let javascript = required(draft.javascript, "javascript")?;
    let feature = required(draft.feature, "feature")?;
    require_nonempty(&semantic_id, "semantic_id")?;
    require_nonempty(&python, "python")?;
    require_nonempty(&javascript, "javascript")?;
    require_nonempty(&feature, "feature")?;
    validate_cfg(&cfg_attrs, &feature)?;
    let requires = match draft.requires {
        None => Vec::new(),
        Some(values) => {
            validate_required_capabilities(&feature, &values)?;
            values
        }
    };

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
            let kind = parse_kind(&kind_ident)?;
            let state = parse_state(&state_ident)?;
            let receiver = match draft.receiver {
                Some(value) => Some(match value.to_string().as_str() {
                    "shared" => Receiver::Shared,
                    "mutable" => Receiver::Mutable,
                    "owned" => Receiver::Owned,
                    _ => {
                        return Err(syn::Error::new_spanned(
                            value,
                            "receiver must be shared, mutable or owned",
                        ));
                    }
                }),
                None if kind == CallableKind::Instance => Some(if state == StateModel::InPlace {
                    Receiver::Mutable
                } else {
                    Receiver::Shared
                }),
                None => None,
            };
            let payload = CallablePayload {
                receiver,
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
            reject_present(draft.receiver, "receiver", "type", &item_ident)?;
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
        python_property,
        javascript,
        feature,
        requires,
        status: draft.status,
        callable,
        type_role,
        properties: draft.properties.unwrap_or_default(),
        python_keywords: draft.python_keywords,
        python_adapters: draft.python_adapters.unwrap_or_default(),
    })
}

fn validate_registry(entries: &[BindingEntry]) -> syn::Result<()> {
    for entry in entries {
        if !entry.python_adapters.is_empty() && entry.item != ItemClass::Type {
            return Err(syn::Error::new_spanned(
                &entry.semantic_id,
                "Python adapters require a type declaration",
            ));
        }
        let mut adapter_names = HashSet::new();
        let type_name = rust_last_name(&entry.rust)?;
        for adapter in &entry.python_adapters {
            let name = adapter.name.to_string();
            if !adapter_names.insert(name.clone()) {
                return Err(syn::Error::new_spanned(
                    &adapter.name,
                    "duplicate Python adapter",
                ));
            }
            if entries.iter().any(|candidate| {
                candidate.item == ItemClass::Callable
                    && candidate.semantic_id.value() == format!("{}.{}", type_name, name)
            }) {
                return Err(syn::Error::new_spanned(
                    &adapter.name,
                    "Python adapter duplicates a canonical Rust callable",
                ));
            }
            if adapter.targets.is_empty() {
                return Err(syn::Error::new_spanned(
                    &adapter.name,
                    "Python adapter requires registered Rust targets",
                ));
            }
            let mut targets = HashSet::new();
            for target in &adapter.targets {
                if !targets.insert(target.value()) {
                    return Err(syn::Error::new_spanned(
                        target,
                        "duplicate Python adapter target",
                    ));
                }
                if !entries.iter().any(|candidate| {
                    candidate.item == ItemClass::Callable
                        && candidate.semantic_id.value() == target.value()
                }) {
                    return Err(syn::Error::new_spanned(
                        target,
                        "Python adapter target is not a registered Rust callable",
                    ));
                }
            }
        }
        if !entry.properties.is_empty() && entry.item != ItemClass::Type {
            return Err(syn::Error::new_spanned(
                &entry.semantic_id,
                "properties require a type declaration",
            ));
        }
        let mut property_names = HashSet::new();
        for property in &entry.properties {
            if !property_names.insert(property.name.to_string()) {
                return Err(syn::Error::new_spanned(
                    &property.name,
                    "duplicate read-only property",
                ));
            }
            if property.signature.inputs.len() != 1
                || property.signature.unsafety.is_some()
                || property.signature.variadic.is_some()
            {
                return Err(syn::Error::new_spanned(
                    &property.signature,
                    "property requires a safe unary getter signature",
                ));
            }
            let syn::Type::Reference(receiver) = &property.signature.inputs[0].ty else {
                return Err(syn::Error::new_spanned(
                    &property.signature,
                    "property requires a shared receiver",
                ));
            };
            let expected = &entry.rust;
            let element = &receiver.elem;
            if receiver.mutability.is_some()
                || quote!(#expected).to_string() != quote!(#element).to_string()
            {
                return Err(syn::Error::new_spanned(
                    &property.signature,
                    "property receiver must match its declared type",
                ));
            }
        }
        if let Some(projection) = &entry.python_keywords {
            if entry.item != ItemClass::Callable {
                return Err(syn::Error::new_spanned(
                    &entry.semantic_id,
                    "python_keywords require a callable",
                ));
            }
            let constructor = entries
                .iter()
                .find(|candidate| candidate.semantic_id.value() == projection.parameters.value())
                .ok_or_else(|| {
                    syn::Error::new_spanned(
                        &projection.parameters,
                        "keyword parameter constructor is not registered",
                    )
                })?;
            let Some(constructor) = &constructor.callable else {
                return Err(syn::Error::new_spanned(
                    &projection.parameters,
                    "keyword parameters must reference a callable constructor",
                ));
            };
            if !matches!(
                constructor.kind,
                CallableKind::Static | CallableKind::Constructor
            ) || constructor
                .parameters
                .iter()
                .any(|parameter| matches!(parameter.default, ParameterDefault::Required))
            {
                return Err(syn::Error::new_spanned(
                    &projection.parameters,
                    "keyword constructor must be static with explicit defaults",
                ));
            }
            let target = entries
                .iter()
                .find(|candidate| candidate.semantic_id.value() == projection.target.value())
                .ok_or_else(|| {
                    syn::Error::new_spanned(&projection.target, "keyword target is not registered")
                })?;
            let Some(target) = &target.callable else {
                return Err(syn::Error::new_spanned(
                    &projection.target,
                    "keyword target must be callable",
                ));
            };
            if target.parameters.len() != 1 {
                return Err(syn::Error::new_spanned(
                    &projection.target,
                    "keyword target must take one immutable parameter object",
                ));
            }
            let syn::Type::Reference(reference) = &target.parameters[0].ty else {
                return Err(syn::Error::new_spanned(
                    &projection.target,
                    "keyword target must borrow parameters",
                ));
            };
            let element = &reference.elem;
            let output = &constructor.output;
            if reference.mutability.is_some()
                || quote!(#element).to_string() != quote!(#output).to_string()
            {
                return Err(syn::Error::new_spanned(
                    &projection.target,
                    "keyword constructor output must match borrowed target parameter",
                ));
            }
        }
    }
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
            let rust_name = rust_last_name(&entry.rust)?;
            // Value and in-place forms can share their domain verb (sanitize /
            // sanitize_). Preserve the mutation suffix only for such a pair;
            // stripping it would give two different behaviors one JS name.
            let retain_in_place_suffix = payload.state == StateModel::InPlace
                && rust_name.ends_with('_')
                && entries.iter().any(|other| {
                    projection_owner_key(other) == projection_owner_key(entry)
                        && other
                            .callable
                            .as_ref()
                            .is_some_and(|callable| callable.state == StateModel::ValueReturning)
                        && rust_last_name(&other.rust)
                            .is_ok_and(|name| name == rust_name.trim_end_matches('_'))
                });
            validate_callable(
                &entry.semantic_id,
                entry.owner,
                &entry.rust,
                &entry.python,
                &entry.javascript,
                payload,
                retain_in_place_suffix,
                entry.python_property,
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
    let mut python = HashMap::new();
    let mut javascript = HashSet::new();
    for entry in entries {
        let owner = projection_owner_key(entry);
        let key = format!("{owner}:{}", entry.python.value());
        if let Some(previous) = python.insert(key, entry.python_property) {
            if !matches!(
                (previous, entry.python_property),
                (Some(PythonProperty::Getter), Some(PythonProperty::Setter))
                    | (Some(PythonProperty::Setter), Some(PythonProperty::Getter))
            ) {
                return Err(syn::Error::new_spanned(
                    &entry.python,
                    "duplicate Python binding projection in one owner scope",
                ));
            }
        }
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
    retain_in_place_suffix: bool,
    python_property: Option<PythonProperty>,
) -> syn::Result<()> {
    if (payload.kind == CallableKind::Instance) != payload.receiver.is_some() {
        return Err(syn::Error::new_spanned(
            rust,
            "receiver is only valid on an instance callable",
        ));
    }
    if let Some(receiver) = payload.receiver {
        let valid = match receiver {
            Receiver::Shared => payload.state != StateModel::InPlace,
            Receiver::Mutable => payload.state == StateModel::InPlace,
            Receiver::Owned => {
                payload.state == StateModel::ValueReturning
                    && matches!(payload.operation, OperationLink::None)
            }
        };
        if !valid {
            return Err(syn::Error::new_spanned(
                rust,
                "receiver conflicts with state or operation lifecycle",
            ));
        }
    }
    match (owner, payload.kind) {
        (Owner::Molecule, CallableKind::Instance | CallableKind::Static)
        | (
            Owner::Type,
            CallableKind::Instance | CallableKind::Static | CallableKind::Constructor,
        )
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
    if payload.kind == CallableKind::Constructor
        && (rust_name != "new"
            || python.value() != "__new__"
            || payload.state != StateModel::ValueReturning)
    {
        return Err(syn::Error::new_spanned(
            python,
            "constructor projection requires canonical Rust new, Python __new__, and value_returning",
        ));
    }
    let expected_python = match python_property {
        Some(PythonProperty::Getter) => {
            if payload.kind != CallableKind::Instance
                || payload.receiver != Some(Receiver::Shared)
                || payload.state != StateModel::ReadOnly
                || !payload.parameters.is_empty()
            {
                return Err(syn::Error::new_spanned(
                    python,
                    "Python getter requires a shared, zero-argument read-only accessor",
                ));
            }
            rust_name.as_str()
        }
        Some(PythonProperty::Setter) => {
            if payload.kind != CallableKind::Instance
                || payload.receiver != Some(Receiver::Mutable)
                || payload.state != StateModel::InPlace
                || payload.parameters.len() != 1
                || !matches!(&payload.output, Type::Tuple(tuple) if tuple.elems.is_empty())
            {
                return Err(syn::Error::new_spanned(
                    python,
                    "Python setter requires a mutable, single-argument unit accessor",
                ));
            }
            rust_name.strip_prefix("set_").ok_or_else(|| {
                syn::Error::new_spanned(python, "Python setter requires a Rust set_ accessor")
            })?
        }
        None if owner == Owner::Type
            && matches!(
                payload.kind,
                CallableKind::Static | CallableKind::Constructor
            )
            && rust_name == "new"
            && python.value() == "__new__" =>
        {
            "__new__"
        }
        None => rust_name.as_str(),
    };
    if python.value() != expected_python {
        return Err(syn::Error::new_spanned(
            python,
            "Python callable name must equal the canonical Rust name",
        ));
    }
    let structural_object = owner == Owner::Type
        && rust.segments.iter().rev().nth(1).is_some_and(|segment| {
            segment.ident == "BioStructure"
                || segment.ident == "Protein"
                || segment.ident == "MolecularForceField"
        });
    let requires_trailing_underscore =
        (owner == Owner::Molecule || structural_object) && payload.state == StateModel::InPlace;
    if rust_name.ends_with('_') != requires_trailing_underscore {
        return Err(syn::Error::new_spanned(
            rust,
            "only in-place callable names may end in `_`",
        ));
    }
    let mut canonical_js = snake_to_camel(rust_name.trim_end_matches('_'));
    if retain_in_place_suffix {
        canonical_js.push('_');
    }
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
        let receiver_path: Path = match owner {
            Owner::Molecule => syn::parse_quote!(crate::Molecule),
            Owner::Type => {
                let mut path = rust.clone();
                path.segments.pop();
                path.segments.pop_punct();
                if path.segments.is_empty() {
                    return Err(syn::Error::new_spanned(
                        rust,
                        "type-owned callable needs a receiver type",
                    ));
                }
                path
            }
            Owner::Module => {
                return Err(syn::Error::new_spanned(
                    rust,
                    "module callable cannot have an instance receiver",
                ));
            }
        };
        let expected: Type = match payload.receiver.expect("validated instance receiver") {
            Receiver::Shared => syn::parse_quote!(&#receiver_path),
            Receiver::Mutable => syn::parse_quote!(&mut #receiver_path),
            Receiver::Owned => syn::parse_quote!(#receiver_path),
        };
        if receiver_metadata_tokens(inputs[0]) != receiver_metadata_tokens(&expected) {
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
/// binder records which borrowed input owns a borrowed output. Erase lifetime
/// spellings (including lifetime arguments on borrowed view types) only for
/// metadata agreement; the generated const assertion retains the complete
/// declared signature and asks rustc to verify the exact higher-ranked
/// relationship against the public item. Type and const arguments remain
/// exact.
fn signature_metadata_tokens(ty: &Type) -> String {
    struct ElideReferenceLifetimes;

    impl VisitMut for ElideReferenceLifetimes {
        fn visit_type_reference_mut(&mut self, reference: &mut syn::TypeReference) {
            reference.lifetime = None;
            visit_mut::visit_type_reference_mut(self, reference);
        }

        fn visit_lifetime_mut(&mut self, lifetime: &mut syn::Lifetime) {
            *lifetime = syn::parse_quote!('_);
        }
    }

    let mut normalized = ty.clone();
    ElideReferenceLifetimes.visit_type_mut(&mut normalized);
    tokens(&normalized)
}

/// A type-owned method path names its receiver as `path::Type::method` and
/// therefore cannot carry the receiver type's generic arguments. Compare the
/// receiver's ownership/reference shape and canonical base path here; keep
/// the complete declared type for the generated public signature assertion.
fn receiver_metadata_tokens(ty: &Type) -> String {
    fn erase_terminal_path_arguments(ty: &mut Type) {
        match ty {
            Type::Reference(reference) => erase_terminal_path_arguments(&mut reference.elem),
            Type::Paren(paren) => erase_terminal_path_arguments(&mut paren.elem),
            Type::Group(group) => erase_terminal_path_arguments(&mut group.elem),
            Type::Path(path) if path.qself.is_none() => {
                if let Some(segment) = path.path.segments.last_mut() {
                    segment.arguments = syn::PathArguments::None;
                }
            }
            _ => {}
        }
    }

    let mut normalized = ty.clone();
    erase_terminal_path_arguments(&mut normalized);
    signature_metadata_tokens(&normalized)
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

fn validate_required_capabilities(owner: &LitStr, requires: &[LitStr]) -> syn::Result<()> {
    fn is_capability(value: &str) -> bool {
        value.strip_prefix("cap-").is_some_and(|tail| {
            !tail.is_empty()
                && tail.split('-').all(|part| {
                    !part.is_empty()
                        && part
                            .bytes()
                            .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit())
                })
        })
    }
    if requires.is_empty() {
        return Err(syn::Error::new_spanned(
            owner,
            "requires must contain at least one capability",
        ));
    }
    if !is_capability(&owner.value()) {
        return Err(syn::Error::new_spanned(
            owner,
            "requires owner must use a cap- capability spelling",
        ));
    }
    let mut seen = HashSet::new();
    for requirement in requires {
        let value = requirement.value();
        if !is_capability(&value) {
            return Err(syn::Error::new_spanned(
                requirement,
                "requires must use nonempty cap- capability spelling",
            ));
        }
        if value == owner.value() {
            return Err(syn::Error::new_spanned(
                requirement,
                "requires must not repeat its owner capability",
            ));
        }
        if !seen.insert(value) {
            return Err(syn::Error::new_spanned(
                requirement,
                "duplicate required capability",
            ));
        }
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
    // Native archives are deliberately absent from wasm32, including under
    // `full`. This one platform boundary is not a second capability selector.
    let native_archive: Attribute = syn::parse_quote!(
        #[cfg(all(feature = "cap-serialization", not(target_arch = "wasm32")))]
    );
    if feature.value() == "cap-serialization"
        && attribute.to_token_stream().to_string() == native_archive.to_token_stream().to_string()
    {
        return Ok(());
    }
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
    let mut property_values = Vec::new();
    let mut keyword_values = Vec::new();
    let mut adapter_values = Vec::new();
    for (index, entry) in entries.iter().enumerate() {
        // The effective gate is generated once from the owning declaration;
        // registry rows and every signature/type assertion use the same value.
        let requires = &entry.requires;
        let cfg = if requires.is_empty() {
            entry.cfg_attrs.clone()
        } else {
            let owner = &entry.feature;
            let mut attrs = entry.cfg_attrs.clone();
            attrs.push(syn::parse_quote!(#[cfg(all(feature = #owner, #(feature = #requires),*))]));
            attrs
        };
        let semantic_id = &entry.semantic_id;
        for adapter in &entry.python_adapters {
            let name = adapter.name.to_string();
            let targets = &adapter.targets;
            let mut adapter_cfg = cfg.clone();
            for target in targets {
                let target = entries
                    .iter()
                    .find(|row| row.semantic_id.value() == target.value())
                    .expect("validated adapter target");
                adapter_cfg.extend(target.cfg_attrs.clone());
                if !target.requires.is_empty() {
                    let requires = &target.requires;
                    adapter_cfg.push(syn::parse_quote!(#[cfg(all(#(feature = #requires),*))]));
                }
            }
            adapter_values.push(
                quote! { #(#adapter_cfg)* crate::BindingPythonAdapterContract {
                    type_semantic_id: #semantic_id, name: #name, targets: &[#(#targets),*],
                }},
            );
        }
        for property in &entry.properties {
            let property_name = property.name.to_string();
            let getter = &property.rust;
            let signature = &property.signature;
            let output = match &signature.output {
                ReturnType::Type(_, ty) => quote!(#ty),
                ReturnType::Default => quote!(()),
            };
            property_values.push(quote! { #(#cfg)* crate::BindingPropertyContract {
                type_semantic_id: #semantic_id, name: #property_name, rust_path: stringify!(#getter), output_type: stringify!(#output),
            }});
            assertions
                .push(quote! { #(#cfg)* const _: fn() = || { let _: #signature = #getter; }; });
        }
        if let Some(projection) = &entry.python_keywords {
            let parameters = &projection.parameters;
            let target = &projection.target;
            keyword_values.push(quote! { #(#cfg)* crate::BindingKeywordContract {
                semantic_id: #semantic_id, parameters_semantic_id: #parameters, target_semantic_id: #target,
            }});
        }
        let item = item_tokens(entry.item);
        let owner = owner_tokens(entry.owner);
        let rust = &entry.rust;
        let python = &entry.python;
        let python_property = match entry.python_property {
            None => quote!(None),
            Some(PythonProperty::Getter) => quote!(Some(crate::BindingPropertyAccess::Getter)),
            Some(PythonProperty::Setter) => quote!(Some(crate::BindingPropertyAccess::Setter)),
        };
        let javascript = &entry.javascript;
        let feature = &entry.feature;
        let status = match (&entry.callable, &entry.status) {
            (Some(payload), status)
                if entry.owner == Owner::Molecule
                    && matches!(payload.operation, OperationLink::Id(_)) =>
            {
                if status.is_some() {
                    return Err(syn::Error::new_spanned(
                        semantic_id,
                        "linked molecule status is inherited from its operation declaration",
                    ));
                }
                let OperationLink::Id(operation) = &payload.operation else {
                    unreachable!()
                };
                let constant = format_ident!(
                    "__FUNCTION_STATUS_{}",
                    operation.value().to_ascii_uppercase()
                );
                quote!(crate::Molecule::#constant)
            }
            (_, Some(status)) => status.tokens(),
            (_, None) => quote!(crate::FunctionStatus::Experimental),
        };
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
                let receiver = match payload.receiver {
                    Some(Receiver::Shared) => quote!(Some(crate::BindingReceiver::Shared)),
                    Some(Receiver::Mutable) => quote!(Some(crate::BindingReceiver::Mutable)),
                    Some(Receiver::Owned) => quote!(Some(crate::BindingReceiver::Owned)),
                    None => quote!(None),
                };
                (
                    quote!(Some(crate::BindingCallableContract {
                        kind: #kind, receiver: #receiver, parameters: &[#(#parameters),*], output_type: stringify!(#output),
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
                rust_path: stringify!(#rust), python_name: #python, python_property: #python_property, javascript_name: #javascript,
                feature: #feature, required_capabilities: &[#(#requires),*], status: #status,
                callable: #callable, type_role: #role,
            }
        });
        {
            let assertion = format_ident!("__BINDING_ASSERT_{}_{}", name, index);
            if let Some(payload) = &entry.callable {
                let signature = &payload.signature;
                if let Some(bound_lifetimes) = &signature.lifetimes {
                    let parameters = &bound_lifetimes.lifetimes;
                    let mut instantiated_signature = signature.clone();
                    instantiated_signature.lifetimes = None;
                    assertions.push(quote! {
                        #(#cfg)* const #assertion: fn() = || {
                            fn assert_signature<#parameters>() {
                                let _: #instantiated_signature = #rust;
                            }
                        };
                    });
                } else {
                    assertions.push(quote! { #(#cfg)* const #assertion: #signature = #rust; });
                }
            } else {
                assertions.push(quote! { #(#cfg)* const #assertion: fn() = || {
                    fn assert_public_type<T>() {} assert_public_type::<#rust>();
                }; });
            }
        }
    }
    let properties_name = format_ident!("{}_PROPERTIES", name);
    let keywords_name = format_ident!("{}_KEYWORDS", name);
    let adapters_name = format_ident!("{}_PYTHON_ADAPTERS", name);
    let adapters = if adapter_values.is_empty() {
        quote! {}
    } else {
        quote! { #visibility static #adapters_name: &[crate::BindingPythonAdapterContract] = &[#(#adapter_values),*]; }
    };
    let properties = if property_values.is_empty() {
        quote! {}
    } else {
        quote! { #visibility static #properties_name: &[crate::BindingPropertyContract] = &[#(#property_values),*]; }
    };
    let keywords = if keyword_values.is_empty() {
        quote! {}
    } else {
        quote! { #visibility static #keywords_name: &[crate::BindingKeywordContract] = &[#(#keyword_values),*]; }
    };
    Ok(quote! {
        #properties
        #keywords
        #adapters
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

fn parse_kind(v: &Ident) -> syn::Result<CallableKind> {
    parse_enum(
        v,
        &["instance", "static_", "module", "constructor"],
        |n| match n {
            "instance" => CallableKind::Instance,
            "static_" => CallableKind::Static,
            "constructor" => CallableKind::Constructor,
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

fn kind_tokens(v: CallableKind) -> proc_macro2::TokenStream {
    match v {
        CallableKind::Instance => quote!(crate::BindingKind::Instance),
        CallableKind::Static | CallableKind::Constructor => quote!(crate::BindingKind::Static),
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
