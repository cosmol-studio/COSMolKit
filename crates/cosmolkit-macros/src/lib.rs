use proc_macro::{Ident as MacroIdent, Literal, Spacing, TokenStream, TokenTree};
use quote::{format_ident, quote};
use std::path::PathBuf;
use syn::{
    Expr, FnArg, Ident, ItemFn, LitBool, LitStr, Pat, PatType, Path, Token, braced, bracketed,
    parenthesized, parse::Parse, parse::ParseStream, parse_macro_input, parse_quote,
    punctuated::Punctuated, spanned::Spanned,
};

struct MolOpBodyAttr {
    op_name: syn::Ident,
    _comma: Token![,],
    parts_name: syn::Ident,
}

struct MoleculeOpsInput {
    ops: Vec<OpEntry>,
}

struct OpEntry {
    name: Ident,
    params: Punctuated<PatType, Token![,]>,
    fields: OpFields,
}

#[derive(Default)]
struct OpFields {
    method: Option<Ident>,
    impl_fn: Option<Ident>,
    domain: Option<Ident>,
    kind: Option<Ident>,
    topology_edit: Option<Ident>,
    access: Option<AccessFields>,
    may_mutate: Vec<Ident>,
    auto_remap: Vec<Ident>,
    derived_effects: Option<DerivedEffectFields>,
    semantic_preconditions: Vec<Ident>,
    requires_mapping: Option<Ident>,
    allows_noop: Option<LitBool>,
    feature: Option<Path>,
    parity: Option<Ident>,
    io_roundtrip: Option<LitBool>,
    invariant_profile: Option<LitStr>,
    parity_profile: Option<LitStr>,
    default_method: Option<Ident>,
    default_args: Vec<Expr>,
    inplace: Option<LitBool>,
    inplace_method: Option<Ident>,
    default_inplace_method: Option<Ident>,
}

#[derive(Default)]
struct AccessFields {
    read: Vec<Ident>,
    write: Vec<Ident>,
}

#[derive(Default)]
struct DerivedEffectFields {
    recompute: Vec<Ident>,
    preserve: Vec<Ident>,
    invalidate: Vec<Ident>,
}

impl Parse for MolOpBodyAttr {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        Ok(Self {
            op_name: input.parse()?,
            _comma: input.parse()?,
            parts_name: input.parse()?,
        })
    }
}

impl Parse for MoleculeOpsInput {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        let mut ops = Vec::new();
        while !input.is_empty() {
            ops.push(input.parse()?);
        }
        Ok(Self { ops })
    }
}

impl Parse for OpEntry {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        let first: Ident = input.parse()?;
        let name = if first == "op" { input.parse()? } else { first };

        let params = if input.peek(syn::token::Paren) {
            let content;
            parenthesized!(content in input);
            content.parse_terminated(PatType::parse, Token![,])?
        } else {
            Punctuated::new()
        };

        let content;
        braced!(content in input);
        let mut fields = OpFields::default();
        while !content.is_empty() {
            let key: Ident = content.parse()?;
            content.parse::<Token![:]>()?;
            match key.to_string().as_str() {
                "method" => fields.method = Some(content.parse()?),
                "impl_fn" => fields.impl_fn = Some(content.parse()?),
                "domain" => fields.domain = Some(content.parse()?),
                "kind" => fields.kind = Some(content.parse()?),
                "topology_edit" => fields.topology_edit = Some(content.parse()?),
                "access" => fields.access = Some(parse_access_fields(&content)?),
                "may_mutate" => fields.may_mutate = parse_ident_list(&content)?,
                "auto_remap" => fields.auto_remap = parse_ident_list(&content)?,
                "derived_effects" => {
                    fields.derived_effects = Some(parse_derived_effect_fields(&content)?)
                }
                "semantic_preconditions" => {
                    fields.semantic_preconditions = parse_ident_list(&content)?
                }
                "must_handle" => {
                    return Err(syn::Error::new(
                        key.span(),
                        "use `derived_effects: { require_handle: [...], recompute: [...] }`; `must_handle` is a derived compatibility view",
                    ));
                }
                "needs_update" => {
                    return Err(syn::Error::new(
                        key.span(),
                        "use `derived_effects: { invalidate: [...], recompute: [...] }`; `needs_update` is a derived compatibility view",
                    ));
                }
                "invalidates_old" => {
                    return Err(syn::Error::new(
                        key.span(),
                        "use `derived_effects.invalidate` to declare stale cache states that must be cleared or updated",
                    ));
                }
                "invalidates" => {
                    return Err(syn::Error::new(
                        key.span(),
                        "use `derived_effects.invalidate` to declare stale cache states that must be cleared or updated",
                    ));
                }
                "requires_mapping" => fields.requires_mapping = Some(content.parse()?),
                "report" => {
                    let _ = parse_ident_list(&content)?;
                }
                "allows_noop" => fields.allows_noop = Some(content.parse()?),
                "feature" => fields.feature = Some(content.parse()?),
                "parity" => fields.parity = Some(content.parse()?),
                "io_roundtrip" => fields.io_roundtrip = Some(content.parse()?),
                "invariant_profile" => fields.invariant_profile = Some(content.parse()?),
                "parity_profile" => fields.parity_profile = Some(content.parse()?),
                "default_method" => fields.default_method = Some(content.parse()?),
                "default_args" => fields.default_args = parse_expr_list(&content)?,
                "inplace" => fields.inplace = Some(content.parse()?),
                "inplace_method" => fields.inplace_method = Some(content.parse()?),
                "default_inplace_method" => fields.default_inplace_method = Some(content.parse()?),
                "rdkit_parity" => {
                    return Err(syn::Error::new(
                        key.span(),
                        "use `parity: required_when_supported | required_now | not_applicable`; `rdkit_parity` is ambiguous",
                    ));
                }
                other => {
                    return Err(syn::Error::new(
                        key.span(),
                        format!("unknown molecule_ops! field `{other}`"),
                    ));
                }
            }
            if content.peek(Token![,]) {
                content.parse::<Token![,]>()?;
            }
        }

        Ok(Self {
            name,
            params,
            fields,
        })
    }
}

fn parse_ident_list(input: ParseStream<'_>) -> syn::Result<Vec<Ident>> {
    let content;
    bracketed!(content in input);
    let items = content.parse_terminated(Ident::parse, Token![,])?;
    Ok(items.into_iter().collect())
}

fn parse_access_fields(input: ParseStream<'_>) -> syn::Result<AccessFields> {
    let content;
    braced!(content in input);
    let mut fields = AccessFields::default();
    while !content.is_empty() {
        let key: Ident = content.parse()?;
        content.parse::<Token![:]>()?;
        match key.to_string().as_str() {
            "read" => fields.read = parse_ident_list(&content)?,
            "write" => fields.write = parse_ident_list(&content)?,
            other => {
                return Err(syn::Error::new(
                    key.span(),
                    format!("unknown access field `{other}`"),
                ));
            }
        }
        if content.peek(Token![,]) {
            content.parse::<Token![,]>()?;
        }
    }
    Ok(fields)
}

fn parse_derived_effect_fields(input: ParseStream<'_>) -> syn::Result<DerivedEffectFields> {
    let content;
    braced!(content in input);
    let mut fields = DerivedEffectFields::default();
    while !content.is_empty() {
        let key: Ident = content.parse()?;
        content.parse::<Token![:]>()?;
        match key.to_string().as_str() {
            "recompute" => fields.recompute = parse_ident_list(&content)?,
            "preserve" => fields.preserve = parse_ident_list(&content)?,
            "invalidate" => fields.invalidate = parse_ident_list(&content)?,
            "requires" => {
                return Err(syn::Error::new(
                    key.span(),
                    "`derived_effects.requires` was removed; recompute from topology or declare preserve/invalidate as appropriate",
                ));
            }
            "unsupported" => {
                return Err(syn::Error::new(
                    key.span(),
                    "`derived_effects.unsupported` was removed; unsupported behavior should be expressed by returning a structured error",
                ));
            }
            "require_handle" => {
                return Err(syn::Error::new(
                    key.span(),
                    "`require_handle` was removed with `derived_effects.requires`; use recompute/preserve/invalidate only",
                ));
            }
            other => {
                return Err(syn::Error::new(
                    key.span(),
                    format!("unknown derived_effects field `{other}`"),
                ));
            }
        }
        if content.peek(Token![,]) {
            content.parse::<Token![,]>()?;
        }
    }
    Ok(fields)
}

fn parse_expr_list(input: ParseStream<'_>) -> syn::Result<Vec<Expr>> {
    let content;
    bracketed!(content in input);
    let items = content.parse_terminated(Expr::parse, Token![,])?;
    Ok(items.into_iter().collect())
}

#[proc_macro_attribute]
pub fn mol_op_body(attr: TokenStream, item: TokenStream) -> TokenStream {
    let attr = parse_macro_input!(attr as MolOpBodyAttr);
    let mut func = parse_macro_input!(item as ItemFn);

    for input in &func.sig.inputs {
        if matches!(input, FnArg::Receiver(_)) {
            return syn::Error::new(
                input.span(),
                "mol_op_body operation bodies must not receive self",
            )
            .to_compile_error()
            .into();
        }
    }

    let _op_name = attr.op_name;
    let parts_name = attr.parts_name;
    let operation_params = std::mem::take(&mut func.sig.inputs);
    func.sig
        .inputs
        .push(parse_quote!(#parts_name: &mut crate::ops::OpParts<'_>));
    func.sig.inputs.extend(operation_params);

    quote!(#func).into()
}

#[proc_macro_attribute]
pub fn bio_op_body(attr: TokenStream, item: TokenStream) -> TokenStream {
    let attr = parse_macro_input!(attr as MolOpBodyAttr);
    let mut func = parse_macro_input!(item as ItemFn);

    for input in &func.sig.inputs {
        if matches!(input, FnArg::Receiver(_)) {
            return syn::Error::new(
                input.span(),
                "bio_op_body operation bodies must not receive self",
            )
            .to_compile_error()
            .into();
        }
    }

    let _op_name = attr.op_name;
    let parts_name = attr.parts_name;
    let operation_params = std::mem::take(&mut func.sig.inputs);
    func.sig
        .inputs
        .push(parse_quote!(#parts_name: &mut crate::bio_ops::BioOpParts<'_>));
    func.sig.inputs.extend(operation_params);

    quote!(#func).into()
}

#[proc_macro]
pub fn molecule_ops(input: TokenStream) -> TokenStream {
    let registry = parse_macro_input!(input as MoleculeOpsInput);
    expand_molecule_ops(registry)
}

fn expand_molecule_ops(registry: MoleculeOpsInput) -> TokenStream {
    let mut spec_defs = Vec::new();
    let mut op_refs = Vec::new();
    let mut support_entries = Vec::new();
    let mut invariant_entries = Vec::new();
    let mut parity_entries = Vec::new();
    let mut methods = Vec::new();

    for op in registry.ops {
        let expanded = match expand_op(op) {
            Ok(expanded) => expanded,
            Err(err) => return err.to_compile_error().into(),
        };
        spec_defs.push(expanded.spec_def);
        op_refs.push(expanded.op_ref);
        support_entries.push(expanded.support_entry);
        invariant_entries.push(expanded.invariant_entry);
        if let Some(entry) = expanded.parity_entry {
            parity_entries.push(entry);
        }
        methods.push(expanded.method);
        if let Some(method) = expanded.default_method {
            methods.push(method);
        }
        if let Some(method) = expanded.inplace_method {
            methods.push(method);
        }
        if let Some(method) = expanded.default_inplace_method {
            methods.push(method);
        }
    }

    quote! {
        #(#spec_defs)*

        pub const MOLECULE_OPS: &[&crate::ops::MoleculeOpSpec] = &[
            #(#op_refs,)*
        ];

        pub const SUPPORT_MATRIX: &[crate::ops::SupportMatrixEntry] = &[
            #(#support_entries,)*
        ];

        pub const OPERATION_INVARIANT_MATRIX: &[crate::ops::OperationInvariantEntry] = &[
            #(#invariant_entries,)*
        ];

        pub const PARITY_MATRIX: &[crate::ops::ParityMatrixEntry] = &[
            #(#parity_entries,)*
        ];

        impl crate::molecule::Molecule {
            #(#methods)*
        }
    }
    .into()
}

struct ExpandedOp {
    spec_def: proc_macro2::TokenStream,
    op_ref: proc_macro2::TokenStream,
    support_entry: proc_macro2::TokenStream,
    invariant_entry: proc_macro2::TokenStream,
    parity_entry: Option<proc_macro2::TokenStream>,
    method: proc_macro2::TokenStream,
    default_method: Option<proc_macro2::TokenStream>,
    inplace_method: Option<proc_macro2::TokenStream>,
    default_inplace_method: Option<proc_macro2::TokenStream>,
}

fn expand_op(op: OpEntry) -> syn::Result<ExpandedOp> {
    let method = required_field(op.fields.method, &op.name, "method")?;
    let impl_fn = required_field(op.fields.impl_fn, &op.name, "impl_fn")?;
    let kind = required_field(op.fields.kind, &op.name, "kind")?;
    let feature = required_field(op.fields.feature, &op.name, "feature")?;
    let parity = required_field(op.fields.parity, &op.name, "parity")?;
    let invariant_profile =
        required_field(op.fields.invariant_profile, &op.name, "invariant_profile")?;
    let domain = op
        .fields
        .domain
        .map_or_else(|| ident_with_span("topology", op.name.span()), Ok)?;
    let topology_edit = op
        .fields
        .topology_edit
        .map_or_else(|| ident_with_span("none", op.name.span()), Ok)?;
    let access = required_field(op.fields.access, &op.name, "access")?;
    let derived_effects = required_field(op.fields.derived_effects, &op.name, "derived_effects")?;
    let requires_mapping = op
        .fields
        .requires_mapping
        .map_or_else(|| ident_with_span("none", op.name.span()), Ok)?;
    let allows_noop = op
        .fields
        .allows_noop
        .map_or_else(|| parse_quote!(true), |value| value);
    let io_roundtrip = op
        .fields
        .io_roundtrip
        .map_or_else(|| parse_quote!(false), |value| value);
    let inplace_enabled = op.fields.inplace.as_ref().is_some_and(|value| value.value);
    let inplace_method_ident = op.fields.inplace_method.clone();
    let default_inplace_method_ident = op.fields.default_inplace_method.clone();
    let default_method_ident = op.fields.default_method.clone();
    let default_args_for_inplace = op.fields.default_args.clone();

    let spec_ident = format_ident!("{}_SPEC", op.name.to_string().to_ascii_uppercase());
    let method_name = method.to_string();
    let impl_fn_name = impl_fn.to_string();
    let domain_expr = domain_expr(&domain)?;
    let kind_expr = kind_expr(&kind)?;
    let topology_edit_expr = topology_edit_expr(&topology_edit)?;
    let access_expr = access_expr(&access)?;
    let may_mutate_expr = block_set_expr(&op.fields.may_mutate)?;
    let auto_remap_expr = block_set_expr(&op.fields.auto_remap)?;
    let semantic_preconditions_expr =
        semantic_precondition_set_expr(&op.fields.semantic_preconditions)?;
    let recompute_expr = derived_state_expr(&derived_effects.recompute)?;
    let preserve_expr = derived_state_expr(&derived_effects.preserve)?;
    let invalidate_expr = derived_state_expr(&derived_effects.invalidate)?;
    let mapping_expr = mapping_expr(&requires_mapping)?;
    let parity_expr = parity_expr(&parity)?;
    let call_args = param_idents(&op.params)?;
    let params = op.params;

    let spec_def = quote! {
        pub const #spec_ident: crate::ops::MoleculeOpSpec = crate::ops::MoleculeOpSpec {
            method: #method_name,
            impl_fn: #impl_fn_name,
            domain: #domain_expr,
            kind: #kind_expr,
            topology_edit: #topology_edit_expr,
            access: #access_expr,
            may_mutate: #may_mutate_expr,
            auto_remap: #auto_remap_expr,
            derived_effects: crate::ops::DerivedEffects::new(
                #recompute_expr,
                #preserve_expr,
                #invalidate_expr,
            ),
            semantic_preconditions: #semantic_preconditions_expr,
            requires_mapping: #mapping_expr,
            allows_noop: #allows_noop,
            support: crate::#feature.status,
            parity: #parity_expr,
            io_roundtrip: #io_roundtrip,
        };
    };

    let op_ref = quote! { &#spec_ident };
    let support_entry = quote! {
        crate::ops::SupportMatrixEntry {
            feature: &crate::#feature,
            operation: Some(&#spec_ident),
        }
    };
    let invariant_entry = quote! {
        crate::ops::OperationInvariantEntry::for_operation(&#spec_ident, #invariant_profile)
    };

    let parity_entry = if matches!(
        parity.to_string().as_str(),
        "not_applicable" | "NotApplicable"
    ) {
        None
    } else {
        let profile = required_field(op.fields.parity_profile, &op.name, "parity_profile")?;
        Some(quote! {
            crate::ops::ParityMatrixEntry {
                operation: &#spec_ident,
                feature: &crate::#feature,
                profile: #profile,
                rdkit_version: None,
            }
        })
    };

    let method_def = quote! {
        pub fn #method(&self, #params) -> Result<crate::molecule::Molecule, crate::ops::OperationError> {
            if matches!(
                crate::ops::#spec_ident.support,
                crate::SupportStatus::Unsupported { .. }
            ) {
                return Err(crate::ops::OperationError::UnsupportedFeature {
                    operation: &crate::ops::#spec_ident,
                    source: crate::UnsupportedFeatureError::from_spec(&crate::#feature),
                });
            }
            let mut parts = crate::ops::OpParts::new(self, &crate::ops::#spec_ident)?;
            let outcome = #impl_fn(&mut parts, #(#call_args),*)?;
            parts.finish(outcome)
        }
    };

    let inplace_method_name = inplace_method_ident
        .clone()
        .unwrap_or_else(|| format_ident!("{}_", method));
    let inplace_method = if inplace_enabled {
        Some(quote! {
            pub fn #inplace_method_name(&mut self, #params) -> Result<(), crate::ops::OperationError> {
                if matches!(
                    crate::ops::#spec_ident.support,
                    crate::SupportStatus::Unsupported { .. }
                ) {
                    return Err(crate::ops::OperationError::UnsupportedFeature {
                        operation: &crate::ops::#spec_ident,
                        source: crate::UnsupportedFeatureError::from_spec(&crate::#feature),
                    });
                }
                let mut parts = crate::ops::OpParts::new_in_place(self, &crate::ops::#spec_ident)?;
                let outcome = match #impl_fn(&mut parts, #(#call_args),*) {
                    Ok(outcome) => outcome,
                    Err(error) => {
                        parts.abort_in_place();
                        return Err(error);
                    }
                };
                parts.finish_in_place(outcome)
            }
        })
    } else {
        None
    };

    let default_method = if let Some(default_method) = op.fields.default_method {
        let default_args = op.fields.default_args;
        Some(quote! {
            pub fn #default_method(&self) -> Result<crate::molecule::Molecule, crate::ops::OperationError> {
                self.#method(#(#default_args),*)
            }
        })
    } else {
        None
    };

    let default_inplace_method = if inplace_enabled
        && let Some(default_method) = default_method_ident
    {
        let default_inplace_method =
            default_inplace_method_ident.unwrap_or_else(|| format_ident!("{}_", default_method));
        Some(quote! {
            pub fn #default_inplace_method(&mut self) -> Result<(), crate::ops::OperationError> {
                self.#inplace_method_name(#(#default_args_for_inplace),*)
            }
        })
    } else {
        None
    };

    Ok(ExpandedOp {
        spec_def,
        op_ref,
        support_entry,
        invariant_entry,
        parity_entry,
        method: method_def,
        default_method,
        inplace_method,
        default_inplace_method,
    })
}

fn required_field<T>(value: Option<T>, op: &Ident, field: &str) -> syn::Result<T> {
    value.ok_or_else(|| {
        syn::Error::new(
            op.span(),
            format!("operation `{op}` is missing required field `{field}`"),
        )
    })
}

fn ident_with_span(value: &str, span: proc_macro2::Span) -> syn::Result<Ident> {
    Ok(Ident::new(value, span))
}

fn param_idents(params: &Punctuated<PatType, Token![,]>) -> syn::Result<Vec<Ident>> {
    params
        .iter()
        .map(|param| {
            let Pat::Ident(pat_ident) = param.pat.as_ref() else {
                return Err(syn::Error::new(
                    param.pat.span(),
                    "operation parameters must use simple identifier patterns",
                ));
            };
            Ok(pat_ident.ident.clone())
        })
        .collect()
}

fn domain_expr(ident: &Ident) -> syn::Result<proc_macro2::TokenStream> {
    match ident.to_string().as_str() {
        "topology" | "Topology" => Ok(quote!(crate::ops::OperationDomain::Topology)),
        "coordinate" | "Coordinate" => Ok(quote!(crate::ops::OperationDomain::Coordinate)),
        other => Err(syn::Error::new(
            ident.span(),
            format!("unknown operation domain `{other}`"),
        )),
    }
}

fn kind_expr(ident: &Ident) -> syn::Result<proc_macro2::TokenStream> {
    match ident.to_string().as_str() {
        "strong" | "Strong" => Ok(quote!(crate::ops::MoleculeOpKind::Strong)),
        "weak" | "Weak" => Ok(quote!(crate::ops::MoleculeOpKind::Weak)),
        other => Err(syn::Error::new(
            ident.span(),
            format!("unknown operation kind `{other}`"),
        )),
    }
}

fn topology_edit_expr(ident: &Ident) -> syn::Result<proc_macro2::TokenStream> {
    match ident.to_string().as_str() {
        "none" | "None" => Ok(quote!(crate::ops::TopologyEditKind::None)),
        "local" | "Local" => Ok(quote!(crate::ops::TopologyEditKind::Local)),
        "compacting" | "Compacting" => Ok(quote!(crate::ops::TopologyEditKind::Compacting)),
        "appending" | "Appending" | "append" | "Append" => {
            Ok(quote!(crate::ops::TopologyEditKind::Appending))
        }
        "renumbering" | "Renumbering" => Ok(quote!(crate::ops::TopologyEditKind::Renumbering)),
        "merge" | "Merge" => Ok(quote!(crate::ops::TopologyEditKind::Merge)),
        other => Err(syn::Error::new(
            ident.span(),
            format!("unknown topology edit kind `{other}`"),
        )),
    }
}

fn mapping_expr(ident: &Ident) -> syn::Result<proc_macro2::TokenStream> {
    match ident.to_string().as_str() {
        "none" | "None" => Ok(quote!(crate::ops::MappingRequirement::None)),
        "identity" | "Identity" => Ok(quote!(crate::ops::MappingRequirement::Identity)),
        "required" | "Required" => Ok(quote!(crate::ops::MappingRequirement::Required)),
        other => Err(syn::Error::new(
            ident.span(),
            format!("unknown mapping requirement `{other}`"),
        )),
    }
}

fn parity_expr(ident: &Ident) -> syn::Result<proc_macro2::TokenStream> {
    match ident.to_string().as_str() {
        "not_applicable" | "NotApplicable" => Ok(quote!(crate::ops::ParityPolicy::NotApplicable)),
        "required_when_supported" | "RequiredWhenSupported" => {
            Ok(quote!(crate::ops::ParityPolicy::RequiredWhenSupported))
        }
        "required_now" | "RequiredNow" => Ok(quote!(crate::ops::ParityPolicy::RequiredNow)),
        other => Err(syn::Error::new(
            ident.span(),
            format!("unknown parity policy `{other}`"),
        )),
    }
}

fn block_set_expr(items: &[Ident]) -> syn::Result<proc_macro2::TokenStream> {
    union_expr(
        items,
        quote!(crate::ops::BlockSet::NONE),
        |item| match item.to_string().as_str() {
            "topology" | "Topology" => Ok(quote!(crate::ops::BlockSet::TOPOLOGY)),
            "coordinates" | "Coordinates" | "conformers" | "Conformers" => {
                Ok(quote!(crate::ops::BlockSet::COORDINATES))
            }
            "properties" | "Properties" | "props" | "Props" => {
                Ok(quote!(crate::ops::BlockSet::PROPERTIES))
            }
            "derived_cache" | "DerivedCache" => Ok(quote!(crate::ops::BlockSet::DERIVED_CACHE)),
            other => Err(syn::Error::new(
                item.span(),
                format!("unknown mutable block `{other}`"),
            )),
        },
    )
}

fn access_expr(access: &AccessFields) -> syn::Result<proc_macro2::TokenStream> {
    let read = block_set_expr(&access.read)?;
    let write = block_set_expr(&access.write)?;
    Ok(quote!(crate::ops::BlockAccess::new(#read, #write)))
}

fn derived_state_expr(items: &[Ident]) -> syn::Result<proc_macro2::TokenStream> {
    union_expr(
        items,
        quote!(crate::DerivedState::NONE),
        |item| match item.to_string().as_str() {
            "rings" | "Rings" => Ok(quote!(crate::DerivedState::RINGS)),
            "ring_families" | "RingFamilies" => Ok(quote!(crate::DerivedState::RING_FAMILIES)),
            "valence" | "Valence" => Ok(quote!(crate::DerivedState::VALENCE)),
            "aromaticity" | "Aromaticity" => Ok(quote!(crate::DerivedState::AROMATICITY)),
            "stereo" | "Stereo" => Ok(quote!(crate::DerivedState::STEREO)),
            "coordinates" | "Coordinates" => Ok(quote!(crate::DerivedState::COORDINATES)),
            "drawing" | "Drawing" => Ok(quote!(crate::DerivedState::DRAWING)),
            "fingerprint" | "Fingerprint" => Ok(quote!(crate::DerivedState::FINGERPRINT)),
            other => Err(syn::Error::new(
                item.span(),
                format!("unknown derived state `{other}`"),
            )),
        },
    )
}

fn semantic_precondition_set_expr(items: &[Ident]) -> syn::Result<proc_macro2::TokenStream> {
    union_expr(
        items,
        quote!(crate::ops::SemanticPreconditionSet::NONE),
        |item| match item.to_string().as_str() {
            "trusted_bond_topology" | "TrustedBondTopology" => Ok(quote!(
                crate::ops::SemanticPreconditionSet::TRUSTED_BOND_TOPOLOGY
            )),
            "hydrogen_ownership_represented" | "HydrogenOwnershipRepresented" => Ok(quote!(
                crate::ops::SemanticPreconditionSet::HYDROGEN_OWNERSHIP_REPRESENTED
            )),
            "valence_computable" | "ValenceComputable" => Ok(quote!(
                crate::ops::SemanticPreconditionSet::VALENCE_COMPUTABLE
            )),
            other => Err(syn::Error::new(
                item.span(),
                format!("unknown semantic precondition `{other}`"),
            )),
        },
    )
}

fn union_expr<F>(
    items: &[Ident],
    none: proc_macro2::TokenStream,
    mut map: F,
) -> syn::Result<proc_macro2::TokenStream>
where
    F: FnMut(&Ident) -> syn::Result<proc_macro2::TokenStream>,
{
    let mut expr = none;
    for item in items {
        let next = map(item)?;
        expr = quote!(#expr.union(#next));
    }
    Ok(expr)
}

#[proc_macro]
pub fn rdkit_uff_params(input: TokenStream) -> TokenStream {
    let path = manifest_relative_path(input);
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", path.display()));
    let mut out = String::from("&[");
    for (line_no, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        assert!(
            cols.len() == 3,
            "{}:{} expected 3 tab-separated columns",
            path.display(),
            line_no + 1
        );
        let key = escape_str(cols[0]);
        let r1: f64 = cols[1]
            .parse()
            .unwrap_or_else(|err| panic!("{}:{} invalid r1: {err}", path.display(), line_no + 1));
        let xi: f64 = cols[2]
            .parse()
            .unwrap_or_else(|err| panic!("{}:{} invalid Xi: {err}", path.display(), line_no + 1));
        out.push_str(&format!("(\"{key}\",{r1:?}_f64,{xi:?}_f64),"));
    }
    out.push(']');
    out.parse().expect("generated UFF params tokens")
}

#[proc_macro]
pub fn rdkit_uff_param_match(input: TokenStream) -> TokenStream {
    let (path, label_ident) = manifest_relative_path_and_ident(input);
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", path.display()));
    let mut out = format!("match {label_ident} {{");
    for (line_no, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        assert!(
            cols.len() == 3,
            "{}:{} expected 3 tab-separated columns",
            path.display(),
            line_no + 1
        );
        let key = escape_str(cols[0]);
        let r1: f64 = cols[1]
            .parse()
            .unwrap_or_else(|err| panic!("{}:{} invalid r1: {err}", path.display(), line_no + 1));
        let xi: f64 = cols[2]
            .parse()
            .unwrap_or_else(|err| panic!("{}:{} invalid Xi: {err}", path.display(), line_no + 1));
        out.push_str(&format!(
            "\"{key}\"=>Some(UffAtomicParams{{r1:{r1:?}_f64,gmp_xi:{xi:?}_f64}}),"
        ));
    }
    out.push_str("_=>None}");
    out.parse().expect("generated UFF param match tokens")
}

#[proc_macro]
pub fn rdkit_periodic_rvdw(input: TokenStream) -> TokenStream {
    let path = manifest_relative_path(input);
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", path.display()));
    let mut values = [1.7_f64; 119];
    for (line_no, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        assert!(
            cols.len() == 2,
            "{}:{} expected 2 tab-separated columns",
            path.display(),
            line_no + 1
        );
        let atomic_num: u8 = cols[0].parse().unwrap_or_else(|err| {
            panic!(
                "{}:{} invalid atomic number: {err}",
                path.display(),
                line_no + 1
            )
        });
        let rvdw: f64 = cols[1]
            .parse()
            .unwrap_or_else(|err| panic!("{}:{} invalid rvdw: {err}", path.display(), line_no + 1));
        values[usize::from(atomic_num)] = rvdw;
    }
    let mut out = String::from("[");
    for value in values {
        out.push_str(&format!("{value:?}_f64,"));
    }
    out.push(']');
    out.parse().expect("generated periodic rvdw tokens")
}

fn manifest_relative_path(input: TokenStream) -> PathBuf {
    let Some(TokenTree::Literal(lit)) = input.into_iter().next() else {
        panic!("expected a single string literal path");
    };
    let relative = parse_string_literal(lit);
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR must be set");
    PathBuf::from(manifest_dir).join(relative)
}

fn manifest_relative_path_and_ident(input: TokenStream) -> (PathBuf, MacroIdent) {
    let mut iter = input.into_iter();
    let Some(TokenTree::Literal(lit)) = iter.next() else {
        panic!("expected a string literal path followed by an identifier");
    };
    match iter.next() {
        Some(TokenTree::Punct(punct))
            if punct.as_char() == ',' && punct.spacing() == Spacing::Alone => {}
        _ => panic!("expected comma after path"),
    }
    let Some(TokenTree::Ident(ident)) = iter.next() else {
        panic!("expected identifier after comma");
    };
    if iter.next().is_some() {
        panic!("unexpected tokens after identifier");
    }
    let relative = parse_string_literal(lit);
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR must be set");
    (PathBuf::from(manifest_dir).join(relative), ident)
}

fn parse_string_literal(lit: Literal) -> String {
    let source = lit.to_string();
    source
        .strip_prefix('"')
        .and_then(|value| value.strip_suffix('"'))
        .expect("expected a plain string literal path")
        .to_string()
}

fn escape_str(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
        .replace('\r', "\\r")
        .replace('\t', "\\t")
}

// ---------------------------------------------------------------------------
// bio_structure_ops! macro
// ---------------------------------------------------------------------------

struct BioOpsInput {
    ops: Vec<BioOpEntry>,
}

struct BioOpEntry {
    name: Ident,
    params: Punctuated<PatType, Token![,]>,
    fields: BioOpFields,
}

#[derive(Default)]
struct BioOpFields {
    method: Option<Ident>,
    impl_fn: Option<Ident>,
    domain: Option<Ident>,
    kind: Option<Ident>,
    edit_kind: Option<Ident>,
    may_mutate: Vec<Ident>,
    auto_remap: Vec<Ident>,
    must_handle: Vec<Ident>,
    needs_update: Vec<Ident>,
    requires_mapping: Option<Ident>,
    allows_noop: Option<LitBool>,
    parity: Option<Ident>,
    io_roundtrip: Option<LitBool>,
    invariant_profile: Option<LitStr>,
    parity_profile: Option<LitStr>,
    feature: Option<Path>,
}

impl Parse for BioOpsInput {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        let mut ops = Vec::new();
        while !input.is_empty() {
            ops.push(input.parse()?);
        }
        Ok(Self { ops })
    }
}

impl Parse for BioOpEntry {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        let first: Ident = input.parse()?;
        let name = if first == "op" { input.parse()? } else { first };

        let params = if input.peek(syn::token::Paren) {
            let content;
            parenthesized!(content in input);
            content.parse_terminated(PatType::parse, Token![,])?
        } else {
            Punctuated::new()
        };

        let content;
        braced!(content in input);
        let mut fields = BioOpFields::default();
        while !content.is_empty() {
            let key: Ident = content.parse()?;
            content.parse::<Token![:]>()?;
            match key.to_string().as_str() {
                "method" => fields.method = Some(content.parse()?),
                "impl_fn" => fields.impl_fn = Some(content.parse()?),
                "domain" => fields.domain = Some(content.parse()?),
                "kind" => fields.kind = Some(content.parse()?),
                "edit_kind" => fields.edit_kind = Some(content.parse()?),
                "may_mutate" => fields.may_mutate = parse_ident_list(&content)?,
                "auto_remap" => fields.auto_remap = parse_ident_list(&content)?,
                "must_handle" => fields.must_handle = parse_ident_list(&content)?,
                "needs_update" => fields.needs_update = parse_ident_list(&content)?,
                "requires_mapping" => fields.requires_mapping = Some(content.parse()?),
                "allows_noop" => fields.allows_noop = Some(content.parse()?),
                "parity" => fields.parity = Some(content.parse()?),
                "io_roundtrip" => fields.io_roundtrip = Some(content.parse()?),
                "invariant_profile" => fields.invariant_profile = Some(content.parse()?),
                "parity_profile" => fields.parity_profile = Some(content.parse()?),
                "feature" => fields.feature = Some(content.parse()?),
                other => {
                    return Err(syn::Error::new(
                        key.span(),
                        format!("unknown bio_structure_ops! field `{other}`"),
                    ));
                }
            }
            if content.peek(Token![,]) {
                content.parse::<Token![,]>()?;
            }
        }
        Ok(Self {
            name,
            params,
            fields,
        })
    }
}

#[proc_macro]
pub fn bio_structure_ops(input: TokenStream) -> TokenStream {
    let registry = parse_macro_input!(input as BioOpsInput);
    expand_bio_structure_ops(registry)
}

fn expand_bio_structure_ops(registry: BioOpsInput) -> TokenStream {
    let mut spec_defs = Vec::new();
    let mut op_refs = Vec::new();
    let mut support_entries = Vec::new();
    let mut invariant_entries = Vec::new();
    let mut parity_entries = Vec::new();
    let mut methods = Vec::new();

    for op in registry.ops {
        let expanded = match expand_bio_op(op) {
            Ok(e) => e,
            Err(err) => return err.to_compile_error().into(),
        };
        spec_defs.push(expanded.spec_def);
        op_refs.push(expanded.op_ref);
        support_entries.push(expanded.support_entry);
        invariant_entries.push(expanded.invariant_entry);
        if let Some(e) = expanded.parity_entry {
            parity_entries.push(e);
        }
        methods.push(expanded.method);
    }

    quote! {
        #(#spec_defs)*

        pub const BIO_STRUCTURE_OPS: &[&crate::bio_ops::BioStructureOpSpec] = &[
            #(#op_refs,)*
        ];

        pub const BIO_SUPPORT_MATRIX: &[crate::bio_ops::BioSupportMatrixEntry] = &[
            #(#support_entries,)*
        ];

        pub const BIO_OPERATION_INVARIANT_MATRIX: &[crate::bio_ops::BioOperationInvariantEntry] = &[
            #(#invariant_entries,)*
        ];

        pub const BIO_PARITY_MATRIX: &[crate::bio_ops::BioParityMatrixEntry] = &[
            #(#parity_entries,)*
        ];

        impl crate::bio::BioStructure {
            #(#methods)*
        }
    }
    .into()
}

struct ExpandedBioOp {
    spec_def: proc_macro2::TokenStream,
    op_ref: proc_macro2::TokenStream,
    support_entry: proc_macro2::TokenStream,
    invariant_entry: proc_macro2::TokenStream,
    parity_entry: Option<proc_macro2::TokenStream>,
    method: proc_macro2::TokenStream,
}

fn expand_bio_op(op: BioOpEntry) -> syn::Result<ExpandedBioOp> {
    let method = required_field(op.fields.method, &op.name, "method")?;
    let impl_fn = required_field(op.fields.impl_fn, &op.name, "impl_fn")?;
    let kind = required_field(op.fields.kind, &op.name, "kind")?;
    let feature = required_field(op.fields.feature, &op.name, "feature")?;
    let parity = required_field(op.fields.parity, &op.name, "parity")?;
    let invariant_profile =
        required_field(op.fields.invariant_profile, &op.name, "invariant_profile")?;

    let domain = op
        .fields
        .domain
        .map_or_else(|| ident_with_span("hierarchy", op.name.span()), Ok)?;
    let edit_kind = op
        .fields
        .edit_kind
        .map_or_else(|| ident_with_span("none", op.name.span()), Ok)?;
    let requires_mapping = op
        .fields
        .requires_mapping
        .map_or_else(|| ident_with_span("none", op.name.span()), Ok)?;
    let allows_noop = op
        .fields
        .allows_noop
        .map_or_else(|| parse_quote!(true), |v| v);
    let io_roundtrip = op
        .fields
        .io_roundtrip
        .map_or_else(|| parse_quote!(false), |v| v);

    let spec_ident = format_ident!("BIO_{}_SPEC", op.name.to_string().to_ascii_uppercase());
    let method_name = method.to_string();
    let impl_fn_name = impl_fn.to_string();
    let domain_expr = bio_domain_expr(&domain)?;
    let kind_expr = bio_kind_expr(&kind)?;
    let edit_kind_expr = bio_edit_kind_expr(&edit_kind)?;
    let may_mutate_expr = bio_block_set_expr(&op.fields.may_mutate)?;
    let auto_remap_expr = bio_block_set_expr(&op.fields.auto_remap)?;
    let must_handle_expr = bio_state_set_expr(&op.fields.must_handle)?;
    let needs_update_expr = bio_derived_state_expr(&op.fields.needs_update)?;
    let mapping_expr = bio_mapping_expr(&requires_mapping)?;
    let parity_expr = bio_parity_expr(&parity)?;
    let call_args = param_idents(&op.params)?;
    let params = op.params;

    let support_expr = quote!(crate::#feature.status);

    let spec_def = quote! {
        pub const #spec_ident: crate::bio_ops::BioStructureOpSpec = crate::bio_ops::BioStructureOpSpec {
            method: #method_name,
            impl_fn: #impl_fn_name,
            domain: #domain_expr,
            kind: #kind_expr,
            edit_kind: #edit_kind_expr,
            may_mutate: #may_mutate_expr,
            auto_remap: #auto_remap_expr,
            must_handle: #must_handle_expr,
            needs_update: #needs_update_expr,
            requires_mapping: #mapping_expr,
            allows_noop: #allows_noop,
            support: #support_expr,
            parity: #parity_expr,
            io_roundtrip: #io_roundtrip,
        };
    };

    let op_ref = quote! { &#spec_ident };
    let support_entry = quote! {
        crate::bio_ops::BioSupportMatrixEntry {
            feature: &crate::#feature,
            operation: &#spec_ident,
        }
    };
    let invariant_entry = quote! {
        crate::bio_ops::BioOperationInvariantEntry {
            operation: &#spec_ident,
            profile: #invariant_profile,
        }
    };

    // GemmiWhenApplicable / BiopythonWhenApplicable / PdbSpecRequired mean
    // "parity required once support is promoted" — no entry until then.
    // Only RequiredNow demands an immediate parity_profile and matrix entry.
    let parity_entry = match parity.to_string().as_str() {
        "not_applicable"
        | "NotApplicable"
        | "gemmi_when_applicable"
        | "GemmiWhenApplicable"
        | "biopython_when_applicable"
        | "BiopythonWhenApplicable"
        | "pdb_spec_required"
        | "PdbSpecRequired" => None,
        "required_now" | "RequiredNow" => {
            let profile = required_field(op.fields.parity_profile, &op.name, "parity_profile")?;
            Some(quote! {
                crate::bio_ops::BioParityMatrixEntry {
                    operation: &#spec_ident,
                    profile: #profile,
                }
            })
        }
        other => {
            return Err(syn::Error::new(
                parity.span(),
                format!("unknown bio parity policy `{other}`"),
            ));
        }
    };

    let method_def = quote! {
        pub fn #method(&self, #params) -> Result<crate::bio::BioStructure, crate::bio_ops::BioOperationError> {
            if let crate::SupportStatus::Unsupported { reason } = crate::bio_ops::#spec_ident.support {
                return Err(crate::bio_ops::BioOperationError::Unsupported {
                    operation: &crate::bio_ops::#spec_ident,
                    reason,
                });
            }
            let mut parts = crate::bio_ops::BioOpParts::new(self, &crate::bio_ops::#spec_ident);
            let outcome = #impl_fn(&mut parts, #(#call_args),*)?;
            parts.finish(outcome)
        }
    };

    Ok(ExpandedBioOp {
        spec_def,
        op_ref,
        support_entry,
        invariant_entry,
        parity_entry,
        method: method_def,
    })
}

fn bio_domain_expr(ident: &Ident) -> syn::Result<proc_macro2::TokenStream> {
    match ident.to_string().as_str() {
        "selection" | "Selection" => Ok(quote!(crate::bio_ops::BioOpDomain::Selection)),
        "hierarchy" | "Hierarchy" => Ok(quote!(crate::bio_ops::BioOpDomain::Hierarchy)),
        "coordinate" | "Coordinate" => Ok(quote!(crate::bio_ops::BioOpDomain::Coordinate)),
        "assembly" | "Assembly" => Ok(quote!(crate::bio_ops::BioOpDomain::Assembly)),
        "annotation" | "Annotation" => Ok(quote!(crate::bio_ops::BioOpDomain::Annotation)),
        "bonding" | "Bonding" => Ok(quote!(crate::bio_ops::BioOpDomain::Bonding)),
        "polymer" | "Polymer" => Ok(quote!(crate::bio_ops::BioOpDomain::Polymer)),
        "chemistry_bridge" | "ChemistryBridge" => {
            Ok(quote!(crate::bio_ops::BioOpDomain::ChemistryBridge))
        }
        other => Err(syn::Error::new(
            ident.span(),
            format!("unknown bio op domain `{other}`"),
        )),
    }
}

fn bio_kind_expr(ident: &Ident) -> syn::Result<proc_macro2::TokenStream> {
    match ident.to_string().as_str() {
        "strong" | "Strong" => Ok(quote!(crate::bio_ops::BioOpKind::Strong)),
        "weak" | "Weak" => Ok(quote!(crate::bio_ops::BioOpKind::Weak)),
        other => Err(syn::Error::new(
            ident.span(),
            format!("unknown bio op kind `{other}`"),
        )),
    }
}

fn bio_edit_kind_expr(ident: &Ident) -> syn::Result<proc_macro2::TokenStream> {
    match ident.to_string().as_str() {
        "none" | "None" => Ok(quote!(crate::bio_ops::BioEditKind::None)),
        "local" | "Local" => Ok(quote!(crate::bio_ops::BioEditKind::Local)),
        "compacting" | "Compacting" => Ok(quote!(crate::bio_ops::BioEditKind::Compacting)),
        "expanding" | "Expanding" => Ok(quote!(crate::bio_ops::BioEditKind::Expanding)),
        "renumbering" | "Renumbering" => Ok(quote!(crate::bio_ops::BioEditKind::Renumbering)),
        "splitting" | "Splitting" => Ok(quote!(crate::bio_ops::BioEditKind::Splitting)),
        "merging" | "Merging" => Ok(quote!(crate::bio_ops::BioEditKind::Merging)),
        "transforming" | "Transforming" => Ok(quote!(crate::bio_ops::BioEditKind::Transforming)),
        other => Err(syn::Error::new(
            ident.span(),
            format!("unknown bio edit kind `{other}`"),
        )),
    }
}

fn bio_mapping_expr(ident: &Ident) -> syn::Result<proc_macro2::TokenStream> {
    match ident.to_string().as_str() {
        "none" | "None" => Ok(quote!(crate::bio_ops::MappingRequirement::None)),
        "identity" | "Identity" => Ok(quote!(crate::bio_ops::MappingRequirement::Identity)),
        "required" | "Required" => Ok(quote!(crate::bio_ops::MappingRequirement::Required)),
        other => Err(syn::Error::new(
            ident.span(),
            format!("unknown mapping requirement `{other}`"),
        )),
    }
}

fn bio_parity_expr(ident: &Ident) -> syn::Result<proc_macro2::TokenStream> {
    match ident.to_string().as_str() {
        "not_applicable" | "NotApplicable" => {
            Ok(quote!(crate::bio_ops::BioParityPolicy::NotApplicable))
        }
        "gemmi_when_applicable" | "GemmiWhenApplicable" => {
            Ok(quote!(crate::bio_ops::BioParityPolicy::GemmiWhenApplicable))
        }
        "biopython_when_applicable" | "BiopythonWhenApplicable" => Ok(quote!(
            crate::bio_ops::BioParityPolicy::BiopythonWhenApplicable
        )),
        "pdb_spec_required" | "PdbSpecRequired" => {
            Ok(quote!(crate::bio_ops::BioParityPolicy::PdbSpecRequired))
        }
        "required_now" | "RequiredNow" => Ok(quote!(crate::bio_ops::BioParityPolicy::RequiredNow)),
        other => Err(syn::Error::new(
            ident.span(),
            format!("unknown bio parity policy `{other}`"),
        )),
    }
}

fn bio_block_set_expr(items: &[Ident]) -> syn::Result<proc_macro2::TokenStream> {
    union_expr(
        items,
        quote!(crate::bio_ops::BioBlockSet::NONE),
        |item| match item.to_string().as_str() {
            "atoms" | "Atoms" => Ok(quote!(crate::bio_ops::BioBlockSet::ATOMS)),
            "residues" | "Residues" => Ok(quote!(crate::bio_ops::BioBlockSet::RESIDUES)),
            "chains" | "Chains" => Ok(quote!(crate::bio_ops::BioBlockSet::CHAINS)),
            "entities" | "Entities" => Ok(quote!(crate::bio_ops::BioBlockSet::ENTITIES)),
            "models" | "Models" => Ok(quote!(crate::bio_ops::BioBlockSet::MODELS)),
            "coordinates" | "Coordinates" => Ok(quote!(crate::bio_ops::BioBlockSet::COORDINATES)),
            "bonds" | "Bonds" => Ok(quote!(crate::bio_ops::BioBlockSet::BONDS)),
            "assemblies" | "Assemblies" => Ok(quote!(crate::bio_ops::BioBlockSet::ASSEMBLIES)),
            "annotations" | "Annotations" => Ok(quote!(crate::bio_ops::BioBlockSet::ANNOTATIONS)),
            "derived_cache" | "DerivedCache" => {
                Ok(quote!(crate::bio_ops::BioBlockSet::DERIVED_CACHE))
            }
            "properties" | "Properties" => Ok(quote!(crate::bio_ops::BioBlockSet::PROPERTIES)),
            other => Err(syn::Error::new(
                item.span(),
                format!("unknown bio block `{other}`"),
            )),
        },
    )
}

fn bio_state_set_expr(items: &[Ident]) -> syn::Result<proc_macro2::TokenStream> {
    union_expr(
        items,
        quote!(crate::bio_ops::BioStateSet::NONE),
        |item| match item.to_string().as_str() {
            "hierarchy" | "Hierarchy" => Ok(quote!(crate::bio_ops::BioStateSet::HIERARCHY)),
            "residue_spans" | "ResidueSpans" => {
                Ok(quote!(crate::bio_ops::BioStateSet::RESIDUE_SPANS))
            }
            "chain_spans" | "ChainSpans" => Ok(quote!(crate::bio_ops::BioStateSet::CHAIN_SPANS)),
            "model_spans" | "ModelSpans" => Ok(quote!(crate::bio_ops::BioStateSet::MODEL_SPANS)),
            "coordinate_alignment" | "CoordinateAlignment" => {
                Ok(quote!(crate::bio_ops::BioStateSet::COORDINATE_ALIGNMENT))
            }
            "entity_mapping" | "EntityMapping" => {
                Ok(quote!(crate::bio_ops::BioStateSet::ENTITY_MAPPING))
            }
            "altloc_groups" | "AltlocGroups" => {
                Ok(quote!(crate::bio_ops::BioStateSet::ALTLOC_GROUPS))
            }
            "assembly_references" | "AssemblyReferences" => {
                Ok(quote!(crate::bio_ops::BioStateSet::ASSEMBLY_REFERENCES))
            }
            "bond_references" | "BondReferences" => {
                Ok(quote!(crate::bio_ops::BioStateSet::BOND_REFERENCES))
            }
            "selection_provenance" | "SelectionProvenance" => {
                Ok(quote!(crate::bio_ops::BioStateSet::SELECTION_PROVENANCE))
            }
            "polymer_annotation" | "PolymerAnnotation" => {
                Ok(quote!(crate::bio_ops::BioStateSet::POLYMER_ANNOTATION))
            }
            "secondary_structure" | "SecondaryStructure" => {
                Ok(quote!(crate::bio_ops::BioStateSet::SECONDARY_STRUCTURE))
            }
            other => Err(syn::Error::new(
                item.span(),
                format!("unknown bio state `{other}`"),
            )),
        },
    )
}

fn bio_derived_state_expr(items: &[Ident]) -> syn::Result<proc_macro2::TokenStream> {
    union_expr(
        items,
        quote!(crate::bio_ops::BioDerivedState::NONE),
        |item| match item.to_string().as_str() {
            "atom_index" | "AtomIndex" => Ok(quote!(crate::bio_ops::BioDerivedState::ATOM_INDEX)),
            "residue_index" | "ResidueIndex" => {
                Ok(quote!(crate::bio_ops::BioDerivedState::RESIDUE_INDEX))
            }
            "chain_index" | "ChainIndex" => {
                Ok(quote!(crate::bio_ops::BioDerivedState::CHAIN_INDEX))
            }
            "entity_index" | "EntityIndex" => {
                Ok(quote!(crate::bio_ops::BioDerivedState::ENTITY_INDEX))
            }
            "sequence_cache" | "SequenceCache" => {
                Ok(quote!(crate::bio_ops::BioDerivedState::SEQUENCE_CACHE))
            }
            "polymer_cache" | "PolymerCache" => {
                Ok(quote!(crate::bio_ops::BioDerivedState::POLYMER_CACHE))
            }
            "altloc_cache" | "AltlocCache" => {
                Ok(quote!(crate::bio_ops::BioDerivedState::ALTLOC_CACHE))
            }
            "assembly_cache" | "AssemblyCache" => {
                Ok(quote!(crate::bio_ops::BioDerivedState::ASSEMBLY_CACHE))
            }
            "bond_cache" | "BondCache" => Ok(quote!(crate::bio_ops::BioDerivedState::BOND_CACHE)),
            "backbone_geometry" | "BackboneGeometry" => {
                Ok(quote!(crate::bio_ops::BioDerivedState::BACKBONE_GEOMETRY))
            }
            "sidechain_geometry" | "SidechainGeometry" => {
                Ok(quote!(crate::bio_ops::BioDerivedState::SIDECHAIN_GEOMETRY))
            }
            "nucleic_geometry" | "NucleicGeometry" => {
                Ok(quote!(crate::bio_ops::BioDerivedState::NUCLEIC_GEOMETRY))
            }
            "secondary_structure" | "SecondaryStructure" => {
                Ok(quote!(crate::bio_ops::BioDerivedState::SECONDARY_STRUCTURE))
            }
            "contact_map" | "ContactMap" => {
                Ok(quote!(crate::bio_ops::BioDerivedState::CONTACT_MAP))
            }
            "graph_cache" | "GraphCache" => {
                Ok(quote!(crate::bio_ops::BioDerivedState::GRAPH_CACHE))
            }
            other => Err(syn::Error::new(
                item.span(),
                format!("unknown bio derived state `{other}`"),
            )),
        },
    )
}
