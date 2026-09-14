//! Private parsers and validators for operation declarations.
//!
//! This is native proc-macro infrastructure. Runtime values and domain
//! behavior remain owned by the crate in which generated tokens expand.

use std::collections::HashSet;

use quote::format_ident;
use syn::{
    Attribute, Expr, Ident, LitBool, LitStr, Pat, PatType, Path, Token, Type, braced, bracketed,
    parenthesized, parse::Parse, parse::ParseStream,
};

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) enum MoleculeBlock {
    Topology,
    Coordinates,
    Properties,
    DerivedCache,
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) enum DerivedState {
    Rings,
    RingFamilies,
    Valence,
    Aromaticity,
    Stereo,
    Coordinates,
    Drawing,
    Fingerprint,
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) enum SemanticPrecondition {
    TrustedBondTopology,
    HydrogenOwnershipRepresented,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum MoleculeOutput {
    Single,
    Multiple,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum MoleculeDomain {
    Topology,
    Coordinate,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum OperationKind {
    Weak,
    Strong,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum TopologyEditKind {
    None,
    Local,
    Compacting,
    Expanding,
    Reordering,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum MappingRequirement {
    None,
    Identity,
    Required,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum MoleculeParity {
    NotApplicable,
    RequiredWhenSupported,
    RequiredNow,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum CipStatePolicy {
    Preserve,
    Clear,
    Recompute,
    TautomerSourceTransition,
}

#[derive(Clone, Debug)]
pub(crate) struct AccessFields {
    pub(crate) read: Vec<MoleculeBlock>,
    pub(crate) write: Vec<MoleculeBlock>,
}

#[derive(Clone, Debug)]
pub(crate) struct DerivedEffectFields {
    pub(crate) recompute: Vec<DerivedState>,
    pub(crate) preserve: Vec<DerivedState>,
    pub(crate) invalidate: Vec<DerivedState>,
    pub(crate) operation_defined: Vec<DerivedState>,
}

#[derive(Clone)]
pub(crate) struct MoleculeFields {
    pub(crate) method: Ident,
    pub(crate) docs: Option<LitStr>,
    pub(crate) impl_fn: Path,
    pub(crate) output: MoleculeOutput,
    pub(crate) result_type: Option<Type>,
    pub(crate) assemble_fn: Option<Path>,
    pub(crate) domain: MoleculeDomain,
    pub(crate) kind: OperationKind,
    pub(crate) topology_edit: TopologyEditKind,
    pub(crate) access: AccessFields,
    pub(crate) may_mutate: Vec<MoleculeBlock>,
    pub(crate) auto_remap: Vec<MoleculeBlock>,
    pub(crate) derived_effects: DerivedEffectFields,
    pub(crate) cip_state: CipStatePolicy,
    pub(crate) semantic_preconditions: Vec<SemanticPrecondition>,
    pub(crate) requires_mapping: MappingRequirement,
    pub(crate) feature: Path,
    pub(crate) parity: MoleculeParity,
    pub(crate) io_roundtrip: bool,
    pub(crate) invariant_profile: LitStr,
    pub(crate) parity_profile: Option<LitStr>,
    pub(crate) default_method: Option<Ident>,
    pub(crate) default_args: Vec<Expr>,
    pub(crate) inplace: bool,
    pub(crate) inplace_method: Option<Ident>,
    pub(crate) inplace_docs: Option<LitStr>,
    pub(crate) default_inplace_method: Option<Ident>,
}

#[derive(Clone)]
pub(crate) struct MoleculeOperation {
    pub(crate) cfg_attrs: Vec<Attribute>,
    pub(crate) name: Ident,
    pub(crate) params: Vec<PatType>,
    pub(crate) fields: MoleculeFields,
}

#[derive(Clone)]
pub(crate) struct MoleculeRegistry {
    pub(crate) operations: Vec<MoleculeOperation>,
}

#[derive(Default)]
struct RawAccessFields {
    read: Option<Vec<Ident>>,
    write: Option<Vec<Ident>>,
}

#[derive(Default)]
struct RawDerivedEffectFields {
    recompute: Option<Vec<Ident>>,
    preserve: Option<Vec<Ident>>,
    invalidate: Option<Vec<Ident>>,
    operation_defined: Option<Vec<Ident>>,
}

#[derive(Default)]
struct RawMoleculeFields {
    method: Option<Ident>,
    docs: Option<LitStr>,
    impl_fn: Option<Path>,
    output: Option<Ident>,
    result_type: Option<Type>,
    assemble_fn: Option<Path>,
    domain: Option<Ident>,
    kind: Option<Ident>,
    topology_edit: Option<Ident>,
    access: Option<RawAccessFields>,
    may_mutate: Option<Vec<Ident>>,
    auto_remap: Option<Vec<Ident>>,
    derived_effects: Option<RawDerivedEffectFields>,
    cip_state: Option<Ident>,
    semantic_preconditions: Option<Vec<Ident>>,
    requires_mapping: Option<Ident>,
    feature: Option<Path>,
    parity: Option<Ident>,
    io_roundtrip: Option<LitBool>,
    invariant_profile: Option<LitStr>,
    parity_profile: Option<LitStr>,
    default_method: Option<Ident>,
    default_args: Option<Vec<Expr>>,
    inplace: Option<LitBool>,
    inplace_method: Option<Ident>,
    inplace_docs: Option<LitStr>,
    default_inplace_method: Option<Ident>,
}

impl Parse for MoleculeRegistry {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        let mut operations = Vec::new();
        while !input.is_empty() {
            operations.push(input.parse()?);
        }
        validate_molecule_registry(&operations)?;
        Ok(Self { operations })
    }
}

impl Parse for MoleculeOperation {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        let cfg_attrs = parse_cfg_attributes(input)?;
        let first: Ident = input.parse()?;
        let name = if first == "op" { input.parse()? } else { first };
        let params = if input.peek(syn::token::Paren) {
            let content;
            parenthesized!(content in input);
            content
                .parse_terminated(PatType::parse, Token![,])?
                .into_iter()
                .collect()
        } else {
            Vec::new()
        };
        validate_parameters(&name, &params)?;

        let content;
        braced!(content in input);
        let mut raw = RawMoleculeFields::default();
        let mut seen = HashSet::new();
        while !content.is_empty() {
            let key: Ident = content.parse()?;
            reject_duplicate_field(&mut seen, &key, "molecule_ops")?;
            content.parse::<Token![:]>()?;
            match key.to_string().as_str() {
                "method" => raw.method = Some(content.parse()?),
                "docs" => raw.docs = Some(content.parse()?),
                "impl_fn" => raw.impl_fn = Some(content.parse()?),
                "output" => raw.output = Some(content.parse()?),
                "result_type" => raw.result_type = Some(content.parse()?),
                "assemble_fn" => raw.assemble_fn = Some(content.parse()?),
                "domain" => raw.domain = Some(content.parse()?),
                "kind" => raw.kind = Some(content.parse()?),
                "topology_edit" => raw.topology_edit = Some(content.parse()?),
                "access" => raw.access = Some(parse_access_fields(&content)?),
                "may_mutate" => raw.may_mutate = Some(parse_ident_list(&content, "may_mutate")?),
                "auto_remap" => raw.auto_remap = Some(parse_ident_list(&content, "auto_remap")?),
                "derived_effects" => {
                    raw.derived_effects = Some(parse_derived_effect_fields(&content)?)
                }
                "cip_state" => raw.cip_state = Some(content.parse()?),
                "semantic_preconditions" => {
                    raw.semantic_preconditions =
                        Some(parse_ident_list(&content, "semantic_preconditions")?)
                }
                "requires_mapping" => raw.requires_mapping = Some(content.parse()?),
                "feature" => raw.feature = Some(content.parse()?),
                "parity" => raw.parity = Some(content.parse()?),
                "io_roundtrip" => raw.io_roundtrip = Some(content.parse()?),
                "invariant_profile" => raw.invariant_profile = Some(content.parse()?),
                "parity_profile" => raw.parity_profile = Some(content.parse()?),
                "default_method" => raw.default_method = Some(content.parse()?),
                "default_args" => raw.default_args = Some(parse_expr_list(&content)?),
                "inplace" => raw.inplace = Some(content.parse()?),
                "inplace_method" => raw.inplace_method = Some(content.parse()?),
                "inplace_docs" => raw.inplace_docs = Some(content.parse()?),
                "default_inplace_method" => raw.default_inplace_method = Some(content.parse()?),
                "must_handle" => {
                    return Err(syn::Error::new(
                        key.span(),
                        "'must_handle' is BioStructure state and was removed from molecule_ops; use the current derived_effects categories",
                    ));
                }
                "needs_update" => {
                    return Err(syn::Error::new(
                        key.span(),
                        "use derived_effects invalidate/recompute; 'needs_update' is not a molecule declaration field",
                    ));
                }
                "invalidates_old" | "invalidates" => {
                    return Err(syn::Error::new(
                        key.span(),
                        "use derived_effects.invalidate to declare stale derived states",
                    ));
                }
                "rdkit_parity" => {
                    return Err(syn::Error::new(
                        key.span(),
                        "use the structured parity policy; 'rdkit_parity' is ambiguous",
                    ));
                }
                "report" => {
                    return Err(syn::Error::new(
                        key.span(),
                        "'report' was a discarded legacy field and is not canonical",
                    ));
                }
                other => {
                    return Err(syn::Error::new(
                        key.span(),
                        format!("unknown molecule_ops field '{other}'"),
                    ));
                }
            }
            parse_optional_comma(&content)?;
        }

        let fields = finish_molecule_fields(&name, params.len(), raw)?;
        Ok(Self {
            cfg_attrs,
            name,
            params,
            fields,
        })
    }
}

fn finish_molecule_fields(
    operation: &Ident,
    parameter_count: usize,
    raw: RawMoleculeFields,
) -> syn::Result<MoleculeFields> {
    let method = required(raw.method, operation, "method")?;
    reject_value_suffix(&method, "method")?;
    let impl_fn = required(raw.impl_fn, operation, "impl_fn")?;
    let output = parse_molecule_output(raw.output.as_ref())?;
    let domain = parse_molecule_domain(raw.domain.as_ref())?;
    let kind = parse_operation_kind(required(raw.kind, operation, "kind")?)?;
    let topology_edit = parse_topology_edit(raw.topology_edit.as_ref())?;
    let access = finish_access(required(raw.access, operation, "access")?)?;
    let may_mutate = parse_molecule_blocks(raw.may_mutate.clone().unwrap_or_default())?;
    let auto_remap = parse_molecule_blocks(raw.auto_remap.clone().unwrap_or_default())?;
    let derived_effects =
        finish_derived_effects(required(raw.derived_effects, operation, "derived_effects")?)?;
    let cip_ident = required(raw.cip_state, operation, "cip_state")?;
    let cip_state = parse_cip_state(&cip_ident)?;
    let semantic_preconditions =
        parse_semantic_preconditions(raw.semantic_preconditions.clone().unwrap_or_default())?;
    let requires_mapping = parse_mapping_requirement(raw.requires_mapping.as_ref())?;
    let feature = required(raw.feature.clone(), operation, "feature")?;
    let parity_ident = required(raw.parity.clone(), operation, "parity")?;
    let parity = parse_molecule_parity(&parity_ident)?;
    let invariant_profile = required(
        raw.invariant_profile.clone(),
        operation,
        "invariant_profile",
    )?;
    let io_roundtrip = raw.io_roundtrip.as_ref().is_some_and(LitBool::value);
    let inplace = raw.inplace.as_ref().is_some_and(LitBool::value);

    validate_molecule_relationships(
        operation,
        &method,
        parameter_count,
        output,
        raw.result_type.as_ref(),
        raw.assemble_fn.as_ref(),
        kind,
        topology_edit,
        &access,
        &may_mutate,
        &auto_remap,
        &derived_effects,
        cip_state,
        requires_mapping,
        parity,
        raw.parity_profile.as_ref(),
        raw.default_method.as_ref(),
        raw.default_args.as_deref().unwrap_or_default(),
        inplace,
        raw.inplace_method.as_ref(),
        raw.inplace_docs.as_ref(),
        raw.default_inplace_method.as_ref(),
    )?;

    let inplace_method = inplace.then(|| {
        raw.inplace_method
            .clone()
            .unwrap_or_else(|| format_ident!("{}_", method))
    });
    let default_inplace_method = if inplace {
        raw.default_method.as_ref().map(|default_method| {
            raw.default_inplace_method
                .clone()
                .unwrap_or_else(|| format_ident!("{}_", default_method))
        })
    } else {
        None
    };

    Ok(MoleculeFields {
        method,
        docs: raw.docs,
        impl_fn,
        output,
        result_type: raw.result_type,
        assemble_fn: raw.assemble_fn,
        domain,
        kind,
        topology_edit,
        access,
        may_mutate,
        auto_remap,
        derived_effects,
        cip_state,
        semantic_preconditions,
        requires_mapping,
        feature,
        parity,
        io_roundtrip,
        invariant_profile,
        parity_profile: raw.parity_profile,
        default_method: raw.default_method,
        default_args: raw.default_args.unwrap_or_default(),
        inplace,
        inplace_method,
        inplace_docs: raw.inplace_docs,
        default_inplace_method,
    })
}

#[allow(clippy::too_many_arguments)]
fn validate_molecule_relationships(
    operation: &Ident,
    method: &Ident,
    parameter_count: usize,
    output: MoleculeOutput,
    result_type: Option<&Type>,
    assemble_fn: Option<&Path>,
    kind: OperationKind,
    topology_edit: TopologyEditKind,
    access: &AccessFields,
    may_mutate: &[MoleculeBlock],
    auto_remap: &[MoleculeBlock],
    derived_effects: &DerivedEffectFields,
    cip_state: CipStatePolicy,
    requires_mapping: MappingRequirement,
    parity: MoleculeParity,
    parity_profile: Option<&LitStr>,
    default_method: Option<&Ident>,
    default_args: &[Expr],
    inplace: bool,
    inplace_method: Option<&Ident>,
    inplace_docs: Option<&LitStr>,
    default_inplace_method: Option<&Ident>,
) -> syn::Result<()> {
    for block in may_mutate {
        if !access.write.contains(block) {
            return Err(syn::Error::new(
                operation.span(),
                "every may_mutate block must also appear in access.write",
            ));
        }
    }
    for block in auto_remap {
        if !access.write.contains(block) || !may_mutate.contains(block) {
            return Err(syn::Error::new(
                operation.span(),
                "every auto_remap block must appear in access.write and may_mutate",
            ));
        }
    }
    let has_derived_effect = !derived_effects.recompute.is_empty()
        || !derived_effects.preserve.is_empty()
        || !derived_effects.invalidate.is_empty()
        || !derived_effects.operation_defined.is_empty();
    if has_derived_effect && !access.write.contains(&MoleculeBlock::DerivedCache) {
        return Err(syn::Error::new(
            operation.span(),
            "declared derived effects require derived_cache write access",
        ));
    }
    validate_operation_defined(operation, derived_effects)?;
    validate_cip_transition(operation, method, output, access, cip_state)?;

    let changes_indices = matches!(
        topology_edit,
        TopologyEditKind::Compacting | TopologyEditKind::Expanding | TopologyEditKind::Reordering
    );
    if changes_indices
        && (kind != OperationKind::Strong || requires_mapping != MappingRequirement::Required)
    {
        return Err(syn::Error::new(
            operation.span(),
            "index-changing topology edits must be strong and require a mapping",
        ));
    }
    if kind == OperationKind::Strong && !changes_indices {
        return Err(syn::Error::new(
            operation.span(),
            "strong molecule operations must declare an index-changing topology edit",
        ));
    }
    if kind == OperationKind::Weak && changes_indices {
        return Err(syn::Error::new(
            operation.span(),
            "weak molecule operations cannot declare an index-changing topology edit",
        ));
    }

    if output == MoleculeOutput::Multiple && inplace {
        return Err(syn::Error::new(
            operation.span(),
            "multiple-output molecule operations cannot generate an in-place wrapper",
        ));
    }
    match (output, result_type, assemble_fn) {
        (MoleculeOutput::Single, _, Some(path)) => {
            return Err(syn::Error::new_spanned(
                path,
                "assemble_fn is only valid for multiple-output molecule operations",
            ));
        }
        (MoleculeOutput::Single, Some(result), None) if inplace => {
            return Err(syn::Error::new_spanned(
                result,
                "pending result_type cannot generate an in-place wrapper",
            ));
        }
        (MoleculeOutput::Single, Some(result), None) => {
            if !matches!(result, Type::Path(path) if path.qself.is_none()
                && path.path.segments.iter().all(|segment| matches!(segment.arguments, syn::PathArguments::None)))
            {
                return Err(syn::Error::new_spanned(
                    result,
                    "pending result_type must be a plain type path without generic arguments",
                ));
            }
        }
        (MoleculeOutput::Multiple, Some(_), None) => {
            return Err(syn::Error::new(
                operation.span(),
                "a multiple-output result_type requires assemble_fn",
            ));
        }
        (MoleculeOutput::Multiple, None, Some(path)) => {
            return Err(syn::Error::new_spanned(
                path,
                "assemble_fn requires result_type",
            ));
        }
        _ => {}
    }

    if !inplace {
        if let Some(value) = inplace_method {
            return Err(syn::Error::new_spanned(
                value,
                "inplace_method requires inplace: true",
            ));
        }
        if let Some(value) = inplace_docs {
            return Err(syn::Error::new_spanned(
                value,
                "inplace_docs requires inplace: true",
            ));
        }
        if let Some(value) = default_inplace_method {
            return Err(syn::Error::new_spanned(
                value,
                "default_inplace_method requires inplace: true",
            ));
        }
    }
    if let Some(value) = inplace_method {
        require_inplace_suffix(value, "inplace_method")?;
    }
    if let Some(value) = default_method {
        reject_value_suffix(value, "default_method")?;
        // Default arguments replace a trailing parameter suffix. Any leading
        // parameters remain explicit on the generated short method. This
        // keeps one operation/spec for APIs such as `(atom, position,
        // params)` -> `(atom, position)` instead of manufacturing a second
        // operation row for the convenience projection.
        if default_args.len() > parameter_count {
            return Err(syn::Error::new_spanned(
                value,
                format!(
                    "default_args has {} entries but the operation has only {parameter_count} parameters",
                    default_args.len()
                ),
            ));
        }
    } else if !default_args.is_empty() {
        return Err(syn::Error::new(
            operation.span(),
            "default_args requires default_method",
        ));
    }
    if let Some(value) = default_inplace_method {
        if default_method.is_none() {
            return Err(syn::Error::new_spanned(
                value,
                "default_inplace_method requires default_method",
            ));
        }
        require_inplace_suffix(value, "default_inplace_method")?;
    }

    match (parity, parity_profile) {
        (MoleculeParity::NotApplicable, Some(profile)) => {
            return Err(syn::Error::new_spanned(
                profile,
                "parity: not_applicable must not declare parity_profile",
            ));
        }
        (MoleculeParity::RequiredWhenSupported | MoleculeParity::RequiredNow, None) => {
            return Err(syn::Error::new(
                operation.span(),
                "this parity policy requires parity_profile",
            ));
        }
        _ => {}
    }
    Ok(())
}

fn validate_operation_defined(operation: &Ident, effects: &DerivedEffectFields) -> syn::Result<()> {
    if effects.operation_defined.is_empty() {
        return Ok(());
    }
    let allowed_operation = matches!(
        operation.to_string().as_str(),
        "without_hydrogens" | "without_hydrogens_with_params"
    );
    let valence_only = effects.operation_defined == [DerivedState::Valence];
    if allowed_operation && valence_only {
        return Ok(());
    }
    Err(syn::Error::new(
        operation.span(),
        "operation_defined is permitted only for valence in the hydrogen-removal family; widening it requires explicit human-author approval",
    ))
}

fn validate_cip_transition(
    operation: &Ident,
    method: &Ident,
    output: MoleculeOutput,
    access: &AccessFields,
    cip_state: CipStatePolicy,
) -> syn::Result<()> {
    if cip_state != CipStatePolicy::TautomerSourceTransition {
        return Ok(());
    }
    let valid = operation == "enumerate_tautomers_with_options"
        && method == "enumerate_tautomers_with_options"
        && output == MoleculeOutput::Multiple
        && access.write.contains(&MoleculeBlock::Topology)
        && access.write.contains(&MoleculeBlock::Properties);
    if valid {
        Ok(())
    } else {
        Err(syn::Error::new(
            operation.span(),
            "tautomer_source_transition is permitted only for the multiple-output enumerate_tautomers_with_options operation with topology and properties write access",
        ))
    }
}

fn validate_molecule_registry(operations: &[MoleculeOperation]) -> syn::Result<()> {
    let mut operation_names = HashSet::new();
    let mut method_names = HashSet::new();
    for operation in operations {
        let operation_name = operation.name.to_string();
        if !operation_names.insert(operation_name.clone()) {
            return Err(syn::Error::new_spanned(
                &operation.name,
                format!("duplicate molecule operation '{operation_name}'"),
            ));
        }
        for method in operation.generated_methods() {
            let name = method.to_string();
            if !method_names.insert(name.clone()) {
                return Err(syn::Error::new_spanned(
                    method,
                    format!("duplicate generated molecule method '{name}'"),
                ));
            }
        }
    }
    Ok(())
}

impl MoleculeOperation {
    fn generated_methods(&self) -> Vec<&Ident> {
        let mut methods = vec![&self.fields.method];
        if let Some(method) = self.fields.default_method.as_ref() {
            methods.push(method);
        }
        if let Some(method) = self.fields.inplace_method.as_ref() {
            methods.push(method);
        }
        if let Some(method) = self.fields.default_inplace_method.as_ref() {
            methods.push(method);
        }
        methods
    }
}

fn parse_access_fields(input: ParseStream<'_>) -> syn::Result<RawAccessFields> {
    let content;
    braced!(content in input);
    let mut fields = RawAccessFields::default();
    let mut seen = HashSet::new();
    while !content.is_empty() {
        let key: Ident = content.parse()?;
        reject_duplicate_field(&mut seen, &key, "access")?;
        content.parse::<Token![:]>()?;
        match key.to_string().as_str() {
            "read" => fields.read = Some(parse_ident_list(&content, "access.read")?),
            "write" => fields.write = Some(parse_ident_list(&content, "access.write")?),
            other => {
                return Err(syn::Error::new(
                    key.span(),
                    format!("unknown access field '{other}'"),
                ));
            }
        }
        parse_optional_comma(&content)?;
    }
    Ok(fields)
}

fn finish_access(raw: RawAccessFields) -> syn::Result<AccessFields> {
    let read = parse_molecule_blocks(raw.read.unwrap_or_default())?;
    let write = parse_molecule_blocks(raw.write.unwrap_or_default())?;
    if let Some(block) = read.iter().find(|block| write.contains(block)) {
        return Err(syn::Error::new(
            proc_macro2::Span::call_site(),
            format!("block '{block:?}' cannot appear in access.read and access.write"),
        ));
    }
    Ok(AccessFields { read, write })
}

fn parse_derived_effect_fields(input: ParseStream<'_>) -> syn::Result<RawDerivedEffectFields> {
    let content;
    braced!(content in input);
    let mut fields = RawDerivedEffectFields::default();
    let mut seen = HashSet::new();
    while !content.is_empty() {
        let key: Ident = content.parse()?;
        reject_duplicate_field(&mut seen, &key, "derived_effects")?;
        content.parse::<Token![:]>()?;
        match key.to_string().as_str() {
            "recompute" => {
                fields.recompute = Some(parse_ident_list(&content, "derived_effects.recompute")?)
            }
            "preserve" => {
                fields.preserve = Some(parse_ident_list(&content, "derived_effects.preserve")?)
            }
            "invalidate" => {
                fields.invalidate = Some(parse_ident_list(&content, "derived_effects.invalidate")?)
            }
            "operation_defined" => {
                fields.operation_defined = Some(parse_ident_list(
                    &content,
                    "derived_effects.operation_defined",
                )?)
            }
            "requires" => {
                return Err(syn::Error::new(
                    key.span(),
                    "derived_effects.requires was removed; use recompute, preserve, or invalidate",
                ));
            }
            "unsupported" => {
                return Err(syn::Error::new(
                    key.span(),
                    "derived_effects.unsupported was removed; unsupported capabilities return structured errors",
                ));
            }
            "require_handle" => {
                return Err(syn::Error::new(
                    key.span(),
                    "require_handle was removed with derived_effects.requires",
                ));
            }
            other => {
                return Err(syn::Error::new(
                    key.span(),
                    format!("unknown derived_effects field '{other}'"),
                ));
            }
        }
        parse_optional_comma(&content)?;
    }
    Ok(fields)
}

fn finish_derived_effects(raw: RawDerivedEffectFields) -> syn::Result<DerivedEffectFields> {
    let recompute = parse_derived_states(raw.recompute.unwrap_or_default())?;
    let preserve = parse_derived_states(raw.preserve.unwrap_or_default())?;
    let invalidate = parse_derived_states(raw.invalidate.unwrap_or_default())?;
    let operation_defined = parse_derived_states(raw.operation_defined.unwrap_or_default())?;
    let groups = [&recompute, &preserve, &invalidate, &operation_defined];
    let mut seen = HashSet::new();
    for group in groups {
        for state in group {
            if !seen.insert(*state) {
                return Err(syn::Error::new(
                    proc_macro2::Span::call_site(),
                    format!("derived state '{state:?}' appears in more than one effect category"),
                ));
            }
        }
    }
    Ok(DerivedEffectFields {
        recompute,
        preserve,
        invalidate,
        operation_defined,
    })
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) enum BioBlock {
    Atoms,
    Residues,
    Chains,
    Entities,
    Models,
    Coordinates,
    Bonds,
    Assemblies,
    Annotations,
    DerivedCache,
    Properties,
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) enum BioState {
    Hierarchy,
    ResidueSpans,
    ChainSpans,
    ModelSpans,
    CoordinateAlignment,
    EntityMapping,
    AltlocGroups,
    AssemblyReferences,
    BondReferences,
    SelectionProvenance,
    PolymerAnnotation,
    SecondaryStructure,
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) enum BioDerivedState {
    AtomIndex,
    ResidueIndex,
    ChainIndex,
    EntityIndex,
    SequenceCache,
    PolymerCache,
    AltlocCache,
    AssemblyCache,
    BondCache,
    BackboneGeometry,
    SidechainGeometry,
    NucleicGeometry,
    SecondaryStructure,
    ContactMap,
    GraphCache,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum BioDomain {
    Selection,
    Hierarchy,
    Coordinate,
    Assembly,
    Annotation,
    Bonding,
    Polymer,
    ChemistryBridge,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum BioEditKind {
    None,
    Local,
    Compacting,
    Expanding,
    Renumbering,
    Splitting,
    Merging,
    Transforming,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum BioParity {
    NotApplicable,
    GemmiWhenApplicable,
    BiopythonWhenApplicable,
    PdbSpecRequired,
    RequiredNow,
}

#[derive(Clone)]
pub(crate) struct BioFields {
    pub(crate) method: Ident,
    pub(crate) impl_fn: Path,
    pub(crate) domain: BioDomain,
    pub(crate) kind: OperationKind,
    pub(crate) edit_kind: BioEditKind,
    pub(crate) may_mutate: Vec<BioBlock>,
    pub(crate) auto_remap: Vec<BioBlock>,
    pub(crate) must_handle: Vec<BioState>,
    pub(crate) needs_update: Vec<BioDerivedState>,
    pub(crate) requires_mapping: MappingRequirement,
    pub(crate) feature: Path,
    pub(crate) parity: BioParity,
    pub(crate) io_roundtrip: bool,
    pub(crate) invariant_profile: LitStr,
    pub(crate) parity_profile: Option<LitStr>,
}

#[derive(Clone)]
pub(crate) struct BioOperation {
    pub(crate) cfg_attrs: Vec<Attribute>,
    pub(crate) name: Ident,
    pub(crate) params: Vec<PatType>,
    pub(crate) fields: BioFields,
}

#[derive(Clone)]
pub(crate) struct BioRegistry {
    pub(crate) operations: Vec<BioOperation>,
}

#[derive(Default)]
struct RawBioFields {
    method: Option<Ident>,
    impl_fn: Option<Path>,
    domain: Option<Ident>,
    kind: Option<Ident>,
    edit_kind: Option<Ident>,
    may_mutate: Option<Vec<Ident>>,
    auto_remap: Option<Vec<Ident>>,
    must_handle: Option<Vec<Ident>>,
    needs_update: Option<Vec<Ident>>,
    requires_mapping: Option<Ident>,
    feature: Option<Path>,
    parity: Option<Ident>,
    io_roundtrip: Option<LitBool>,
    invariant_profile: Option<LitStr>,
    parity_profile: Option<LitStr>,
}

impl Parse for BioRegistry {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        let mut operations = Vec::new();
        while !input.is_empty() {
            operations.push(input.parse()?);
        }
        validate_bio_registry(&operations)?;
        Ok(Self { operations })
    }
}

impl Parse for BioOperation {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        let cfg_attrs = parse_cfg_attributes(input)?;
        let first: Ident = input.parse()?;
        let name = if first == "op" { input.parse()? } else { first };
        let params = if input.peek(syn::token::Paren) {
            let content;
            parenthesized!(content in input);
            content
                .parse_terminated(PatType::parse, Token![,])?
                .into_iter()
                .collect()
        } else {
            Vec::new()
        };
        validate_parameters(&name, &params)?;

        let content;
        braced!(content in input);
        let mut raw = RawBioFields::default();
        let mut seen = HashSet::new();
        while !content.is_empty() {
            let key: Ident = content.parse()?;
            reject_duplicate_field(&mut seen, &key, "bio_structure_ops")?;
            content.parse::<Token![:]>()?;
            match key.to_string().as_str() {
                "method" => raw.method = Some(content.parse()?),
                "impl_fn" => raw.impl_fn = Some(content.parse()?),
                "domain" => raw.domain = Some(content.parse()?),
                "kind" => raw.kind = Some(content.parse()?),
                "edit_kind" => raw.edit_kind = Some(content.parse()?),
                "may_mutate" => raw.may_mutate = Some(parse_ident_list(&content, "may_mutate")?),
                "auto_remap" => raw.auto_remap = Some(parse_ident_list(&content, "auto_remap")?),
                "must_handle" => raw.must_handle = Some(parse_ident_list(&content, "must_handle")?),
                "needs_update" => {
                    raw.needs_update = Some(parse_ident_list(&content, "needs_update")?)
                }
                "requires_mapping" => raw.requires_mapping = Some(content.parse()?),
                "feature" => raw.feature = Some(content.parse()?),
                "parity" => raw.parity = Some(content.parse()?),
                "io_roundtrip" => raw.io_roundtrip = Some(content.parse()?),
                "invariant_profile" => raw.invariant_profile = Some(content.parse()?),
                "parity_profile" => raw.parity_profile = Some(content.parse()?),
                other => {
                    return Err(syn::Error::new(
                        key.span(),
                        format!("unknown bio_structure_ops field '{other}'"),
                    ));
                }
            }
            parse_optional_comma(&content)?;
        }

        let method = required(raw.method, &name, "method")?;
        reject_value_suffix(&method, "method")?;
        let impl_fn = required(raw.impl_fn, &name, "impl_fn")?;
        let domain = parse_bio_domain(raw.domain.as_ref())?;
        let kind = parse_operation_kind(required(raw.kind, &name, "kind")?)?;
        let edit_kind = parse_bio_edit_kind(raw.edit_kind.as_ref())?;
        let may_mutate = parse_bio_blocks(raw.may_mutate.unwrap_or_default())?;
        let auto_remap = parse_bio_blocks(raw.auto_remap.unwrap_or_default())?;
        let must_handle = parse_bio_states(raw.must_handle.unwrap_or_default())?;
        let needs_update = parse_bio_derived_states(raw.needs_update.unwrap_or_default())?;
        let requires_mapping = parse_mapping_requirement(raw.requires_mapping.as_ref())?;
        let feature = required(raw.feature, &name, "feature")?;
        let parity_ident = required(raw.parity, &name, "parity")?;
        let parity = parse_bio_parity(&parity_ident)?;
        let invariant_profile = required(raw.invariant_profile, &name, "invariant_profile")?;
        validate_bio_relationships(
            &name,
            kind,
            edit_kind,
            &may_mutate,
            &auto_remap,
            requires_mapping,
            parity,
            raw.parity_profile.as_ref(),
        )?;

        Ok(Self {
            cfg_attrs,
            name,
            params,
            fields: BioFields {
                method,
                impl_fn,
                domain,
                kind,
                edit_kind,
                may_mutate,
                auto_remap,
                must_handle,
                needs_update,
                requires_mapping,
                feature,
                parity,
                io_roundtrip: raw.io_roundtrip.as_ref().is_some_and(LitBool::value),
                invariant_profile,
                parity_profile: raw.parity_profile,
            },
        })
    }
}

fn parse_cfg_attributes(input: ParseStream<'_>) -> syn::Result<Vec<Attribute>> {
    let attributes = input.call(Attribute::parse_outer)?;
    for attribute in &attributes {
        if !attribute.path().is_ident("cfg") {
            return Err(syn::Error::new_spanned(
                attribute,
                "operation declarations accept only outer #[cfg(...)] attributes",
            ));
        }
    }
    Ok(attributes)
}

fn validate_bio_relationships(
    operation: &Ident,
    kind: OperationKind,
    edit_kind: BioEditKind,
    may_mutate: &[BioBlock],
    auto_remap: &[BioBlock],
    mapping: MappingRequirement,
    parity: BioParity,
    parity_profile: Option<&LitStr>,
) -> syn::Result<()> {
    for block in auto_remap {
        if !may_mutate.contains(block) {
            return Err(syn::Error::new(
                operation.span(),
                "every Bio auto_remap block must also appear in may_mutate",
            ));
        }
    }
    if matches!(
        edit_kind,
        BioEditKind::Expanding | BioEditKind::Splitting | BioEditKind::Merging
    ) {
        return Err(syn::Error::new(
            operation.span(),
            "Bio expanding, splitting, and merging declarations require the unresolved source-indexed/per-output mapping contract",
        ));
    }
    let changes_indices = matches!(
        edit_kind,
        BioEditKind::Compacting | BioEditKind::Renumbering
    );
    if changes_indices && (kind != OperationKind::Strong || mapping != MappingRequirement::Required)
    {
        return Err(syn::Error::new(
            operation.span(),
            "index-changing Bio operations must be strong and require a mapping",
        ));
    }
    if kind == OperationKind::Strong && !changes_indices {
        return Err(syn::Error::new(
            operation.span(),
            "strong Bio operations currently require compacting or renumbering edits",
        ));
    }
    if kind == OperationKind::Weak && changes_indices {
        return Err(syn::Error::new(
            operation.span(),
            "weak Bio operations cannot compact or renumber hierarchy rows",
        ));
    }
    match (parity, parity_profile) {
        (BioParity::RequiredNow, None) => {
            return Err(syn::Error::new(
                operation.span(),
                "parity: required_now requires parity_profile",
            ));
        }
        (BioParity::RequiredNow, Some(_)) => {}
        (_, Some(profile)) => {
            return Err(syn::Error::new_spanned(
                profile,
                "only parity: required_now may declare a current Bio parity_profile",
            ));
        }
        (_, None) => {}
    }
    Ok(())
}

fn validate_bio_registry(operations: &[BioOperation]) -> syn::Result<()> {
    let mut operation_names = HashSet::new();
    let mut method_names = HashSet::new();
    for operation in operations {
        let operation_name = operation.name.to_string();
        if !operation_names.insert(operation_name.clone()) {
            return Err(syn::Error::new_spanned(
                &operation.name,
                format!("duplicate Bio operation '{operation_name}'"),
            ));
        }
        let method_name = operation.fields.method.to_string();
        if !method_names.insert(method_name.clone()) {
            return Err(syn::Error::new_spanned(
                &operation.fields.method,
                format!("duplicate generated Bio method '{method_name}'"),
            ));
        }
    }
    Ok(())
}

fn validate_parameters(operation: &Ident, params: &[PatType]) -> syn::Result<()> {
    let mut names = HashSet::new();
    for parameter in params {
        let Pat::Ident(pattern) = parameter.pat.as_ref() else {
            return Err(syn::Error::new_spanned(
                &parameter.pat,
                format!("operation '{operation}' parameters must use simple identifiers"),
            ));
        };
        if pattern.subpat.is_some() {
            return Err(syn::Error::new_spanned(
                pattern,
                format!("operation '{operation}' parameters must not use subpatterns"),
            ));
        }
        let name = pattern.ident.to_string();
        if !names.insert(name.clone()) {
            return Err(syn::Error::new_spanned(
                &pattern.ident,
                format!("operation '{operation}' has duplicate parameter '{name}'"),
            ));
        }
    }
    Ok(())
}

fn parse_ident_list(input: ParseStream<'_>, field: &str) -> syn::Result<Vec<Ident>> {
    let content;
    bracketed!(content in input);
    let values = content.parse_terminated(Ident::parse, Token![,])?;
    let mut seen = HashSet::new();
    let mut result = Vec::new();
    for value in values {
        let name = value.to_string();
        if !seen.insert(name.clone()) {
            return Err(syn::Error::new_spanned(
                value,
                format!("duplicate '{name}' in '{field}'"),
            ));
        }
        result.push(value);
    }
    Ok(result)
}

fn parse_expr_list(input: ParseStream<'_>) -> syn::Result<Vec<Expr>> {
    let content;
    bracketed!(content in input);
    Ok(content
        .parse_terminated(Expr::parse, Token![,])?
        .into_iter()
        .collect())
}

fn parse_optional_comma(input: ParseStream<'_>) -> syn::Result<()> {
    if input.peek(Token![,]) {
        input.parse::<Token![,]>()?;
    }
    Ok(())
}

fn reject_duplicate_field(seen: &mut HashSet<String>, key: &Ident, owner: &str) -> syn::Result<()> {
    let value = key.to_string();
    if seen.insert(value.clone()) {
        Ok(())
    } else {
        Err(syn::Error::new_spanned(
            key,
            format!("duplicate {owner} field '{value}'"),
        ))
    }
}

fn required<T>(value: Option<T>, operation: &Ident, field: &str) -> syn::Result<T> {
    value.ok_or_else(|| {
        syn::Error::new(
            operation.span(),
            format!("operation '{operation}' is missing '{field}'"),
        )
    })
}

fn reject_value_suffix(method: &Ident, field: &str) -> syn::Result<()> {
    if method.to_string().ends_with('_') {
        Err(syn::Error::new_spanned(
            method,
            format!("value-style '{field}' must not end with '_'"),
        ))
    } else {
        Ok(())
    }
}

fn require_inplace_suffix(method: &Ident, field: &str) -> syn::Result<()> {
    if method.to_string().ends_with('_') {
        Ok(())
    } else {
        Err(syn::Error::new_spanned(
            method,
            format!("in-place '{field}' must end with '_'"),
        ))
    }
}

fn parse_molecule_blocks(values: Vec<Ident>) -> syn::Result<Vec<MoleculeBlock>> {
    values
        .into_iter()
        .map(|value| match value.to_string().as_str() {
            "topology" => Ok(MoleculeBlock::Topology),
            "coordinates" => Ok(MoleculeBlock::Coordinates),
            "properties" => Ok(MoleculeBlock::Properties),
            "derived_cache" => Ok(MoleculeBlock::DerivedCache),
            other => Err(syn::Error::new_spanned(
                value,
                format!("unknown molecule block '{other}'"),
            )),
        })
        .collect()
}

fn parse_derived_states(values: Vec<Ident>) -> syn::Result<Vec<DerivedState>> {
    values
        .into_iter()
        .map(|value| match value.to_string().as_str() {
            "rings" => Ok(DerivedState::Rings),
            "ring_families" => Ok(DerivedState::RingFamilies),
            "valence" => Ok(DerivedState::Valence),
            "aromaticity" => Ok(DerivedState::Aromaticity),
            "stereo" => Ok(DerivedState::Stereo),
            "coordinates" => Ok(DerivedState::Coordinates),
            "drawing" => Ok(DerivedState::Drawing),
            "fingerprint" => Ok(DerivedState::Fingerprint),
            other => Err(syn::Error::new_spanned(
                value,
                format!("unknown derived state '{other}'"),
            )),
        })
        .collect()
}

fn parse_semantic_preconditions(values: Vec<Ident>) -> syn::Result<Vec<SemanticPrecondition>> {
    values
        .into_iter()
        .map(|value| match value.to_string().as_str() {
            "trusted_bond_topology" => Ok(SemanticPrecondition::TrustedBondTopology),
            "hydrogen_ownership_represented" => {
                Ok(SemanticPrecondition::HydrogenOwnershipRepresented)
            }
            other => Err(syn::Error::new_spanned(
                value,
                format!("unknown semantic precondition '{other}'"),
            )),
        })
        .collect()
}

fn parse_molecule_output(value: Option<&Ident>) -> syn::Result<MoleculeOutput> {
    match value.map(ToString::to_string).as_deref() {
        None | Some("single") => Ok(MoleculeOutput::Single),
        Some("multiple") => Ok(MoleculeOutput::Multiple),
        Some(other) => Err(syn::Error::new_spanned(
            value.expect("present value"),
            format!("unknown molecule output '{other}'"),
        )),
    }
}

fn parse_molecule_domain(value: Option<&Ident>) -> syn::Result<MoleculeDomain> {
    match value.map(ToString::to_string).as_deref() {
        None | Some("topology") => Ok(MoleculeDomain::Topology),
        Some("coordinate") => Ok(MoleculeDomain::Coordinate),
        Some(other) => Err(syn::Error::new_spanned(
            value.expect("present value"),
            format!("unknown operation domain '{other}'"),
        )),
    }
}

fn parse_operation_kind(value: Ident) -> syn::Result<OperationKind> {
    match value.to_string().as_str() {
        "weak" => Ok(OperationKind::Weak),
        "strong" => Ok(OperationKind::Strong),
        other => Err(syn::Error::new_spanned(
            value,
            format!("unknown operation kind '{other}'"),
        )),
    }
}

fn parse_topology_edit(value: Option<&Ident>) -> syn::Result<TopologyEditKind> {
    match value.map(ToString::to_string).as_deref() {
        None | Some("none") => Ok(TopologyEditKind::None),
        Some("local") => Ok(TopologyEditKind::Local),
        Some("compacting") => Ok(TopologyEditKind::Compacting),
        Some("expanding") => Ok(TopologyEditKind::Expanding),
        Some("reordering") => Ok(TopologyEditKind::Reordering),
        Some(other) => Err(syn::Error::new_spanned(
            value.expect("present value"),
            format!("unknown topology edit kind '{other}'"),
        )),
    }
}

fn parse_mapping_requirement(value: Option<&Ident>) -> syn::Result<MappingRequirement> {
    match value.map(ToString::to_string).as_deref() {
        None | Some("none") => Ok(MappingRequirement::None),
        Some("identity") => Ok(MappingRequirement::Identity),
        Some("required") => Ok(MappingRequirement::Required),
        Some(other) => Err(syn::Error::new_spanned(
            value.expect("present value"),
            format!("unknown mapping requirement '{other}'"),
        )),
    }
}

fn parse_molecule_parity(value: &Ident) -> syn::Result<MoleculeParity> {
    match value.to_string().as_str() {
        "not_applicable" => Ok(MoleculeParity::NotApplicable),
        "required_when_supported" => Ok(MoleculeParity::RequiredWhenSupported),
        "required_now" => Ok(MoleculeParity::RequiredNow),
        other => Err(syn::Error::new_spanned(
            value,
            format!("unknown molecule parity policy '{other}'"),
        )),
    }
}

fn parse_cip_state(value: &Ident) -> syn::Result<CipStatePolicy> {
    match value.to_string().as_str() {
        "preserve" => Ok(CipStatePolicy::Preserve),
        "clear" => Ok(CipStatePolicy::Clear),
        "recompute" => Ok(CipStatePolicy::Recompute),
        "tautomer_source_transition" => Ok(CipStatePolicy::TautomerSourceTransition),
        other => Err(syn::Error::new_spanned(
            value,
            format!("unknown CIP state policy '{other}'"),
        )),
    }
}

fn parse_bio_blocks(values: Vec<Ident>) -> syn::Result<Vec<BioBlock>> {
    values
        .into_iter()
        .map(|value| match value.to_string().as_str() {
            "atoms" => Ok(BioBlock::Atoms),
            "residues" => Ok(BioBlock::Residues),
            "chains" => Ok(BioBlock::Chains),
            "entities" => Ok(BioBlock::Entities),
            "models" => Ok(BioBlock::Models),
            "coordinates" => Ok(BioBlock::Coordinates),
            "bonds" => Ok(BioBlock::Bonds),
            "assemblies" => Ok(BioBlock::Assemblies),
            "annotations" => Ok(BioBlock::Annotations),
            "derived_cache" => Ok(BioBlock::DerivedCache),
            "properties" => Ok(BioBlock::Properties),
            other => Err(syn::Error::new_spanned(
                value,
                format!("unknown Bio block '{other}'"),
            )),
        })
        .collect()
}

fn parse_bio_states(values: Vec<Ident>) -> syn::Result<Vec<BioState>> {
    values
        .into_iter()
        .map(|value| match value.to_string().as_str() {
            "hierarchy" => Ok(BioState::Hierarchy),
            "residue_spans" => Ok(BioState::ResidueSpans),
            "chain_spans" => Ok(BioState::ChainSpans),
            "model_spans" => Ok(BioState::ModelSpans),
            "coordinate_alignment" => Ok(BioState::CoordinateAlignment),
            "entity_mapping" => Ok(BioState::EntityMapping),
            "altloc_groups" => Ok(BioState::AltlocGroups),
            "assembly_references" => Ok(BioState::AssemblyReferences),
            "bond_references" => Ok(BioState::BondReferences),
            "selection_provenance" => Ok(BioState::SelectionProvenance),
            "polymer_annotation" => Ok(BioState::PolymerAnnotation),
            "secondary_structure" => Ok(BioState::SecondaryStructure),
            other => Err(syn::Error::new_spanned(
                value,
                format!("unknown Bio handled state '{other}'"),
            )),
        })
        .collect()
}

fn parse_bio_derived_states(values: Vec<Ident>) -> syn::Result<Vec<BioDerivedState>> {
    values
        .into_iter()
        .map(|value| match value.to_string().as_str() {
            "atom_index" => Ok(BioDerivedState::AtomIndex),
            "residue_index" => Ok(BioDerivedState::ResidueIndex),
            "chain_index" => Ok(BioDerivedState::ChainIndex),
            "entity_index" => Ok(BioDerivedState::EntityIndex),
            "sequence_cache" => Ok(BioDerivedState::SequenceCache),
            "polymer_cache" => Ok(BioDerivedState::PolymerCache),
            "altloc_cache" => Ok(BioDerivedState::AltlocCache),
            "assembly_cache" => Ok(BioDerivedState::AssemblyCache),
            "bond_cache" => Ok(BioDerivedState::BondCache),
            "backbone_geometry" => Ok(BioDerivedState::BackboneGeometry),
            "sidechain_geometry" => Ok(BioDerivedState::SidechainGeometry),
            "nucleic_geometry" => Ok(BioDerivedState::NucleicGeometry),
            "secondary_structure" => Ok(BioDerivedState::SecondaryStructure),
            "contact_map" => Ok(BioDerivedState::ContactMap),
            "graph_cache" => Ok(BioDerivedState::GraphCache),
            other => Err(syn::Error::new_spanned(
                value,
                format!("unknown Bio derived state '{other}'"),
            )),
        })
        .collect()
}

fn parse_bio_domain(value: Option<&Ident>) -> syn::Result<BioDomain> {
    match value.map(ToString::to_string).as_deref() {
        None | Some("hierarchy") => Ok(BioDomain::Hierarchy),
        Some("selection") => Ok(BioDomain::Selection),
        Some("coordinate") => Ok(BioDomain::Coordinate),
        Some("assembly") => Ok(BioDomain::Assembly),
        Some("annotation") => Ok(BioDomain::Annotation),
        Some("bonding") => Ok(BioDomain::Bonding),
        Some("polymer") => Ok(BioDomain::Polymer),
        Some("chemistry_bridge") => Ok(BioDomain::ChemistryBridge),
        Some(other) => Err(syn::Error::new_spanned(
            value.expect("present value"),
            format!("unknown Bio operation domain '{other}'"),
        )),
    }
}

fn parse_bio_edit_kind(value: Option<&Ident>) -> syn::Result<BioEditKind> {
    match value.map(ToString::to_string).as_deref() {
        None | Some("none") => Ok(BioEditKind::None),
        Some("local") => Ok(BioEditKind::Local),
        Some("compacting") => Ok(BioEditKind::Compacting),
        Some("expanding") => Ok(BioEditKind::Expanding),
        Some("renumbering") => Ok(BioEditKind::Renumbering),
        Some("splitting") => Ok(BioEditKind::Splitting),
        Some("merging") => Ok(BioEditKind::Merging),
        Some("transforming") => Ok(BioEditKind::Transforming),
        Some(other) => Err(syn::Error::new_spanned(
            value.expect("present value"),
            format!("unknown Bio edit kind '{other}'"),
        )),
    }
}

fn parse_bio_parity(value: &Ident) -> syn::Result<BioParity> {
    match value.to_string().as_str() {
        "not_applicable" => Ok(BioParity::NotApplicable),
        "gemmi_when_applicable" => Ok(BioParity::GemmiWhenApplicable),
        "biopython_when_applicable" => Ok(BioParity::BiopythonWhenApplicable),
        "pdb_spec_required" => Ok(BioParity::PdbSpecRequired),
        "required_now" => Ok(BioParity::RequiredNow),
        other => Err(syn::Error::new_spanned(
            value,
            format!("unknown Bio parity policy '{other}'"),
        )),
    }
}
