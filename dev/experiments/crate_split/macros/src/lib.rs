use proc_macro::TokenStream;
use quote::{ToTokens, format_ident, quote};
use syn::{
    FnArg, Ident, ItemFn, LitBool, LitStr, Pat, PatType, Path, Token, Type, braced, bracketed,
    parenthesized, parse::Parse, parse::ParseStream, parse_macro_input, parse_quote,
    punctuated::Punctuated, spanned::Spanned,
};

struct BodyAttribute {
    operation: Ident,
    _comma: Token![,],
    parts: Ident,
}

impl Parse for BodyAttribute {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        Ok(Self {
            operation: input.parse()?,
            _comma: input.parse()?,
            parts: input.parse()?,
        })
    }
}

#[proc_macro_attribute]
pub fn mol_op_body(attribute: TokenStream, item: TokenStream) -> TokenStream {
    expand_body(attribute, item, false)
}

#[proc_macro_attribute]
pub fn mol_multi_op_body(attribute: TokenStream, item: TokenStream) -> TokenStream {
    expand_body(attribute, item, true)
}

fn expand_body(attribute: TokenStream, item: TokenStream, multiple_output: bool) -> TokenStream {
    let attribute = parse_macro_input!(attribute as BodyAttribute);
    let mut function = parse_macro_input!(item as ItemFn);
    for input in &function.sig.inputs {
        if matches!(input, FnArg::Receiver(_)) {
            return syn::Error::new(input.span(), "operation bodies must not receive self")
                .to_compile_error()
                .into();
        }
    }
    let operation = attribute.operation;
    let parts = attribute.parts;
    let parameters = std::mem::take(&mut function.sig.inputs);
    let access = context_name(&operation, "Access");
    let context_type: Type = if multiple_output {
        parse_quote!(&mut crate::MultiOutputOpParts<'_, crate::#access>)
    } else {
        parse_quote!(&mut crate::OpParts<'_, crate::#access>)
    };
    function
        .sig
        .inputs
        .push(parse_quote!(#parts: #context_type));
    function.sig.inputs.extend(parameters);
    quote!(#function).into()
}

struct Registry {
    operations: Vec<Operation>,
}

struct Operation {
    name: Ident,
    parameters: Punctuated<PatType, Token![,]>,
    fields: Fields,
}

#[derive(Default)]
struct Fields {
    method: Option<Ident>,
    impl_fn: Option<Ident>,
    output: Option<Ident>,
    result_type: Option<Type>,
    domain: Option<Ident>,
    kind: Option<Ident>,
    topology_edit: Option<Ident>,
    access: Option<Access>,
    may_mutate: Vec<Ident>,
    auto_remap: Vec<Ident>,
    derived_effects: Option<Effects>,
    cip_state: Option<Ident>,
    semantic_preconditions: Vec<Ident>,
    requires_mapping: Option<Ident>,
    feature: Option<Path>,
    parity: Option<Ident>,
    parity_profile: Option<LitStr>,
    io_roundtrip: Option<LitBool>,
    invariant_profile: Option<LitStr>,
    inplace: Option<LitBool>,
    inplace_method: Option<Ident>,
}

#[derive(Default)]
struct Access {
    read: Vec<Ident>,
    write: Vec<Ident>,
}

#[derive(Default)]
struct Effects {
    recompute: Vec<Ident>,
    preserve: Vec<Ident>,
    invalidate: Vec<Ident>,
    operation_defined: Vec<Ident>,
}

impl Parse for Registry {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        let mut operations = Vec::new();
        while !input.is_empty() {
            operations.push(input.parse()?);
        }
        Ok(Self { operations })
    }
}

impl Parse for Operation {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        let marker: Ident = input.parse()?;
        if marker != "op" {
            return Err(syn::Error::new(marker.span(), "expected `op`"));
        }
        let name = input.parse()?;
        let parameters_content;
        parenthesized!(parameters_content in input);
        let parameters = parameters_content.parse_terminated(PatType::parse, Token![,])?;
        let content;
        braced!(content in input);
        let mut fields = Fields::default();
        while !content.is_empty() {
            let key: Ident = content.parse()?;
            content.parse::<Token![:]>()?;
            match key.to_string().as_str() {
                "method" => fields.method = Some(content.parse()?),
                "impl_fn" => fields.impl_fn = Some(content.parse()?),
                "output" => fields.output = Some(content.parse()?),
                "result_type" => fields.result_type = Some(content.parse()?),
                "domain" => fields.domain = Some(content.parse()?),
                "kind" => fields.kind = Some(content.parse()?),
                "topology_edit" => fields.topology_edit = Some(content.parse()?),
                "access" => fields.access = Some(parse_access(&content)?),
                "may_mutate" => fields.may_mutate = parse_ident_list(&content)?,
                "auto_remap" => fields.auto_remap = parse_ident_list(&content)?,
                "derived_effects" => fields.derived_effects = Some(parse_effects(&content)?),
                "cip_state" => fields.cip_state = Some(content.parse()?),
                "semantic_preconditions" => {
                    fields.semantic_preconditions = parse_ident_list(&content)?
                }
                "requires_mapping" => fields.requires_mapping = Some(content.parse()?),
                "feature" => fields.feature = Some(content.parse()?),
                "parity" => fields.parity = Some(content.parse()?),
                "parity_profile" => fields.parity_profile = Some(content.parse()?),
                "io_roundtrip" => fields.io_roundtrip = Some(content.parse()?),
                "invariant_profile" => fields.invariant_profile = Some(content.parse()?),
                "inplace" => fields.inplace = Some(content.parse()?),
                "inplace_method" => fields.inplace_method = Some(content.parse()?),
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
            parameters,
            fields,
        })
    }
}

fn parse_ident_list(input: ParseStream<'_>) -> syn::Result<Vec<Ident>> {
    let content;
    bracketed!(content in input);
    Ok(content
        .parse_terminated(Ident::parse, Token![,])?
        .into_iter()
        .collect())
}

fn parse_access(input: ParseStream<'_>) -> syn::Result<Access> {
    let content;
    braced!(content in input);
    let mut access = Access::default();
    while !content.is_empty() {
        let key: Ident = content.parse()?;
        content.parse::<Token![:]>()?;
        match key.to_string().as_str() {
            "read" => access.read = parse_ident_list(&content)?,
            "write" => access.write = parse_ident_list(&content)?,
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
    Ok(access)
}

fn parse_effects(input: ParseStream<'_>) -> syn::Result<Effects> {
    let content;
    braced!(content in input);
    let mut effects = Effects::default();
    while !content.is_empty() {
        let key: Ident = content.parse()?;
        content.parse::<Token![:]>()?;
        match key.to_string().as_str() {
            "recompute" => effects.recompute = parse_ident_list(&content)?,
            "preserve" => effects.preserve = parse_ident_list(&content)?,
            "invalidate" => effects.invalidate = parse_ident_list(&content)?,
            "operation_defined" => effects.operation_defined = parse_ident_list(&content)?,
            "requires" | "unsupported" | "require_handle" => {
                return Err(syn::Error::new(
                    key.span(),
                    "removed derived-effect field is not accepted",
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
    Ok(effects)
}

fn required<T>(value: Option<T>, operation: &Ident, field: &str) -> syn::Result<T> {
    value.ok_or_else(|| {
        syn::Error::new(
            operation.span(),
            format!("operation `{operation}` is missing `{field}`"),
        )
    })
}

fn strings(values: &[Ident]) -> proc_macro2::TokenStream {
    let values = values.iter().map(ToString::to_string);
    quote!(&[#(#values),*])
}

fn block_set(values: &[Ident]) -> syn::Result<proc_macro2::TokenStream> {
    let mut value = quote!(crate::BlockSet::NONE);
    for block in values {
        let block = match block.to_string().as_str() {
            "topology" => quote!(crate::BlockSet::TOPOLOGY),
            "coordinates" => quote!(crate::BlockSet::COORDINATES),
            "properties" => quote!(crate::BlockSet::PROPERTIES),
            other => {
                return Err(syn::Error::new(
                    block.span(),
                    format!("unknown molecule block `{other}`"),
                ));
            }
        };
        value = quote!(#value.union(#block));
    }
    Ok(value)
}

fn reject_duplicates(left: &[Ident], right: &[Ident], label: &str) -> syn::Result<()> {
    if let Some(overlap) = left
        .iter()
        .find(|item| right.iter().any(|candidate| candidate == *item))
    {
        return Err(syn::Error::new(
            overlap.span(),
            format!("`{label}` categories overlap at `{overlap}`"),
        ));
    }
    Ok(())
}

fn same_ident_set(left: &[Ident], right: &[Ident]) -> bool {
    left.len() == right.len()
        && left
            .iter()
            .all(|item| right.iter().any(|candidate| candidate == item))
}

fn access_mode(block: &str, access: &Access) -> proc_macro2::TokenStream {
    if access.write.iter().any(|candidate| candidate == block) {
        quote!(crate::AccessMode::Write)
    } else if access.read.iter().any(|candidate| candidate == block) {
        quote!(crate::AccessMode::Read)
    } else {
        quote!(crate::AccessMode::None)
    }
}

fn parameter_names(parameters: &Punctuated<PatType, Token![,]>) -> syn::Result<Vec<Ident>> {
    parameters
        .iter()
        .map(|parameter| match parameter.pat.as_ref() {
            Pat::Ident(pattern) => Ok(pattern.ident.clone()),
            pattern => Err(syn::Error::new(
                pattern.span(),
                "operation parameters must be simple identifiers",
            )),
        })
        .collect()
}

fn context_name(operation: &Ident, suffix: &str) -> Ident {
    let mut value = String::new();
    for part in operation.to_string().split('_') {
        let mut chars = part.chars();
        if let Some(first) = chars.next() {
            value.extend(first.to_uppercase());
            value.extend(chars);
        }
    }
    value.push_str(suffix);
    format_ident!("{value}")
}

#[proc_macro]
pub fn molecule_ops(input: TokenStream) -> TokenStream {
    let registry = parse_macro_input!(input as Registry);
    match expand_registry(registry) {
        Ok(output) => output.into(),
        Err(error) => error.to_compile_error().into(),
    }
}

fn expand_registry(registry: Registry) -> syn::Result<proc_macro2::TokenStream> {
    let mut specs = Vec::new();
    let mut spec_names = Vec::new();
    let mut contexts = Vec::new();
    let mut methods = Vec::new();

    for operation in registry.operations {
        let name = operation.name;
        let fields = operation.fields;
        let method = required(fields.method, &name, "method")?;
        let implementation = required(fields.impl_fn, &name, "impl_fn")?;
        let output = required(fields.output, &name, "output")?;
        let domain = required(fields.domain, &name, "domain")?;
        let kind = required(fields.kind, &name, "kind")?;
        let topology_edit = required(fields.topology_edit, &name, "topology_edit")?;
        let access = required(fields.access, &name, "access")?;
        let effects = required(fields.derived_effects, &name, "derived_effects")?;
        let cip_state = required(fields.cip_state, &name, "cip_state")?;
        let mapping = required(fields.requires_mapping, &name, "requires_mapping")?;
        let feature = required(fields.feature, &name, "feature")?;
        let parity = required(fields.parity, &name, "parity")?;
        let invariant = required(fields.invariant_profile, &name, "invariant_profile")?;
        let io_roundtrip = required(fields.io_roundtrip, &name, "io_roundtrip")?;
        reject_duplicates(&access.read, &access.write, "access read/write")?;
        reject_duplicates(&effects.recompute, &effects.preserve, "derived effects")?;
        reject_duplicates(&effects.recompute, &effects.invalidate, "derived effects")?;
        reject_duplicates(
            &effects.recompute,
            &effects.operation_defined,
            "derived effects",
        )?;
        reject_duplicates(&effects.preserve, &effects.invalidate, "derived effects")?;
        reject_duplicates(
            &effects.preserve,
            &effects.operation_defined,
            "derived effects",
        )?;
        reject_duplicates(
            &effects.invalidate,
            &effects.operation_defined,
            "derived effects",
        )?;
        let multiple = output == "multiple";
        if multiple && fields.inplace.as_ref().is_some_and(|value| value.value) {
            return Err(syn::Error::new(
                output.span(),
                "multiple-output operations cannot generate in-place wrappers",
            ));
        }
        if parity != "not_applicable" && fields.parity_profile.is_none() {
            return Err(syn::Error::new(
                parity.span(),
                "parity-enabled operation requires parity_profile",
            ));
        }
        if !effects.operation_defined.is_empty()
            && !(name == "remove_hydrogens"
                && effects.operation_defined.len() == 1
                && effects.operation_defined[0] == "valence")
        {
            return Err(syn::Error::new(
                name.span(),
                "operation_defined is restricted to remove_hydrogens valence",
            ));
        }

        let spec = format_ident!("{}_SPEC", name.to_string().to_ascii_uppercase());
        let method_name = method.to_string();
        let implementation_name = implementation.to_string();
        let domain = domain.to_string();
        let kind = kind.to_string();
        let topology_edit = topology_edit.to_string();
        let feature_name = feature.to_token_stream().to_string();
        let feature = quote!(stringify!(#feature));
        let parity_name = parity.to_string();
        let parity_profile = fields
            .parity_profile
            .map_or_else(|| quote!(None), |value| quote!(Some(#value)));
        let topology_access = access_mode("topology", &access);
        let coordinate_access = access_mode("coordinates", &access);
        let property_access = access_mode("properties", &access);
        let may_mutate = block_set(&fields.may_mutate)?;
        let auto_remap = block_set(&fields.auto_remap)?;
        if !same_ident_set(&fields.may_mutate, &access.write) {
            return Err(syn::Error::new(
                name.span(),
                "may_mutate must exactly match access.write in this experiment",
            ));
        }
        let recompute = strings(&effects.recompute);
        let preserve = strings(&effects.preserve);
        let invalidate = strings(&effects.invalidate);
        let operation_defined = strings(&effects.operation_defined);
        let preconditions = strings(&fields.semantic_preconditions);
        let output_value = if multiple {
            quote!(crate::OperationOutput::Multiple)
        } else {
            quote!(crate::OperationOutput::Single)
        };
        let result_type = fields.result_type.map_or_else(
            || {
                if multiple {
                    quote!(Some("Vec<Molecule>"))
                } else {
                    quote!(None)
                }
            },
            |result| quote!(Some(stringify!(#result))),
        );
        let mapping_value = if mapping == "required" {
            quote!(crate::MappingRequirement::Required)
        } else if mapping == "identity" {
            quote!(crate::MappingRequirement::Identity)
        } else {
            quote!(crate::MappingRequirement::None)
        };
        let requires_mapping = mapping != "none";
        let cip_value = match cip_state.to_string().as_str() {
            "clear_computed" => quote!(crate::CipStatePolicy::ClearComputed),
            "assign" => quote!(crate::CipStatePolicy::Assign),
            "tautomer_source_transition" => {
                quote!(crate::CipStatePolicy::TautomerSourceTransition)
            }
            _ => quote!(crate::CipStatePolicy::Preserve),
        };

        specs.push(quote! {
            pub const #spec: crate::MoleculeOpSpec = crate::MoleculeOpSpec {
                method: #method_name,
                impl_fn: #implementation_name,
                domain: #domain,
                kind: #kind,
                output: #output_value,
                result_type: #result_type,
                access: crate::BlockAccess {
                    topology: #topology_access,
                    coordinates: #coordinate_access,
                    properties: #property_access,
                },
                may_mutate: #may_mutate,
                topology_edit: Some(#topology_edit),
                auto_remap: #auto_remap,
                derived_effects: crate::DerivedEffects {
                    recompute: #recompute,
                    preserve: #preserve,
                    invalidate: #invalidate,
                    operation_defined: #operation_defined,
                },
                semantic_preconditions: #preconditions,
                requires_mapping: #requires_mapping,
                support: "implemented",
                parity: #parity_name,
                parity_profile: #parity_profile,
                io_roundtrip: if #io_roundtrip { "required" } else { "not_applicable" },
                invariant_profile: #invariant,
                feature: #feature,
                mapping_requirement: #mapping_value,
                cip_state: #cip_value,
            };
        });
        spec_names.push(spec.clone());

        let access_marker = context_name(&name, "Access");
        let can_read_topology = access
            .read
            .iter()
            .chain(access.write.iter())
            .any(|block| block == "topology");
        let can_read_coordinates = access
            .read
            .iter()
            .chain(access.write.iter())
            .any(|block| block == "coordinates");
        let can_read_properties = access
            .read
            .iter()
            .chain(access.write.iter())
            .any(|block| block == "properties");
        let can_write_topology = access.write.iter().any(|block| block == "topology");
        let can_write_coordinates = access.write.iter().any(|block| block == "coordinates");
        let can_write_properties = access.write.iter().any(|block| block == "properties");
        let writes_all = can_write_topology && can_write_coordinates && can_write_properties;

        let topology_read = access.read.iter().any(|block| block == "topology").then(|| {
            quote! {
                pub(crate) fn topology(&self) -> Result<&crate::TopologyBlock, crate::OperationError> {
                    crate::OpParts::read_topology_runtime(self)
                }
            }
        });
        let coordinates_read = access.read.iter().any(|block| block == "coordinates").then(|| {
            quote! {
                pub(crate) fn coordinates(&self) -> Result<&crate::CoordinateBlock, crate::OperationError> {
                    crate::OpParts::read_coordinates_runtime(self)
                }
            }
        });
        let properties_read = access.read.iter().any(|block| block == "properties").then(|| {
            quote! {
                pub(crate) fn properties(&self) -> Result<&crate::MoleculeProperties, crate::OperationError> {
                    crate::OpParts::read_properties_runtime(self)
                }
            }
        });
        let multi_topology_read = can_read_topology.then(|| {
            quote! {
                pub(crate) fn topology(&self) -> Result<&crate::TopologyBlock, crate::OperationError> {
                    crate::MultiOutputOpParts::source_topology_runtime(self)
                }
            }
        });
        let multi_coordinates_read = can_read_coordinates.then(|| {
            quote! {
                pub(crate) fn coordinates(&self) -> Result<&crate::CoordinateBlock, crate::OperationError> {
                    crate::MultiOutputOpParts::source_coordinates_runtime(self)
                }
            }
        });
        let multi_properties_read = can_read_properties.then(|| {
            quote! {
                pub(crate) fn properties(&self) -> Result<&crate::MoleculeProperties, crate::OperationError> {
                    crate::MultiOutputOpParts::source_properties_runtime(self)
                }
            }
        });
        let topology_write = can_write_topology.then(|| {
            quote! {
                pub(crate) fn begin_topology_mut(&mut self) -> Result<crate::TopologyBlock, crate::OperationError> {
                    crate::OpParts::checkout_topology_mut(self)
                }
                pub(crate) fn commit_topology(&mut self, topology: crate::TopologyBlock) -> Result<(), crate::OperationError> {
                    crate::OpParts::install_topology(self, topology)
                }
            }
        });
        let coordinates_write = can_write_coordinates.then(|| {
            quote! {
                pub(crate) fn begin_coordinates_mut(&mut self) -> Result<crate::CoordinateBlock, crate::OperationError> {
                    crate::OpParts::checkout_coordinates_mut(self)
                }
                pub(crate) fn commit_coordinates(&mut self, coordinates: crate::CoordinateBlock) -> Result<(), crate::OperationError> {
                    crate::OpParts::install_coordinates(self, coordinates)
                }
            }
        });
        let properties_write = can_write_properties.then(|| {
            quote! {
                pub(crate) fn begin_properties_mut(&mut self) -> Result<crate::MoleculeProperties, crate::OperationError> {
                    crate::OpParts::checkout_properties_mut(self)
                }
                pub(crate) fn commit_properties(&mut self, properties: crate::MoleculeProperties) -> Result<(), crate::OperationError> {
                    crate::OpParts::install_properties(self, properties)
                }
            }
        });
        let all_writable = writes_all.then(|| {
            quote! {
                pub(crate) fn extract_all_writable(&mut self) -> Result<(crate::TopologyBlock, crate::CoordinateBlock, crate::MoleculeProperties), crate::OperationError> {
                    crate::OpParts::checkout_all_writable(self)
                }
                pub(crate) fn commit_all_writable(&mut self, topology: crate::TopologyBlock, coordinates: crate::CoordinateBlock, properties: crate::MoleculeProperties) -> Result<(), crate::OperationError> {
                    crate::OpParts::install_all_writable(self, topology, coordinates, properties)
                }
            }
        });
        let capability_definition = if multiple {
            quote! {
                pub(crate) struct #access_marker;
                impl<'a> crate::MultiOutputOpParts<'a, #access_marker> {
                    #multi_topology_read
                    #multi_coordinates_read
                    #multi_properties_read
                    pub(crate) fn emit_all(&mut self, candidates: Vec<(crate::TopologyBlock, crate::CoordinateBlock, crate::MoleculeProperties)>) -> Result<(), crate::OperationError> {
                        crate::MultiOutputOpParts::emit_all_runtime(self, candidates)
                    }
                    pub(crate) fn finish(self) -> Result<Vec<crate::Molecule>, crate::OperationError> {
                        crate::MultiOutputOpParts::finish_runtime(self)
                    }
                }
            }
        } else {
            quote! {
                pub(crate) struct #access_marker;
                impl<'a> crate::OpParts<'a, #access_marker> {
                    #topology_read
                    #coordinates_read
                    #properties_read
                    #topology_write
                    #coordinates_write
                    #properties_write
                    #all_writable
                }
            }
        };
        contexts.push(capability_definition);

        let parameters = operation.parameters;
        let arguments = parameter_names(&parameters)?;
        let method_definition = if multiple {
            quote! {
                #[cfg(feature = #feature_name)]
                pub fn #method(&self, #parameters) -> Result<Vec<crate::Molecule>, crate::OperationError> {
                    let mut context = crate::MultiOutputOpParts::<crate::#access_marker>::new(self, &#spec)?;
                    #implementation(&mut context, #(#arguments),*)?;
                    context.finish()
                }
            }
        } else {
            quote! {
                #[cfg(feature = #feature_name)]
                pub fn #method(&self, #parameters) -> Result<crate::Molecule, crate::OperationError> {
                    let mut context = crate::OpParts::<crate::#access_marker>::new(self, &#spec)?;
                    #implementation(&mut context, #(#arguments),*)?;
                    context.finish()
                }
            }
        };
        methods.push(method_definition);

        if fields.inplace.as_ref().is_some_and(|value| value.value) {
            let inplace_method = required(fields.inplace_method, &name, "inplace_method")?;
            methods.push(quote! {
                #[cfg(feature = #feature_name)]
                pub fn #inplace_method(&mut self, #parameters) -> Result<(), crate::OperationError> {
                    let mut context = crate::OpParts::<crate::#access_marker>::new_in_place(self, &#spec)?;
                    if let Err(error) = #implementation(&mut context, #(#arguments),*) {
                        context.abort_in_place();
                        return Err(error);
                    }
                    context.finish_in_place()
                }
            });
        }
    }

    let support_entries = spec_names
        .iter()
        .map(|spec| quote!((#spec.method, #spec.support)));
    let invariant_entries = spec_names
        .iter()
        .map(|spec| quote!((#spec.method, #spec.invariant_profile)));
    let parity_entries = spec_names
        .iter()
        .map(|spec| quote!((#spec.method, #spec.parity, #spec.parity_profile)));

    Ok(quote! {
        #(#specs)*
        mod __operation_capabilities {
            #(#contexts)*
        }
        pub(crate) use __operation_capabilities::*;
        pub static MOLECULE_OPS: &[crate::MoleculeOpSpec] = &[#(#spec_names),*];
        pub static SUPPORT_MATRIX: &[(&str, &str)] = &[#(#support_entries),*];
        pub static OPERATION_INVARIANT_MATRIX: &[(&str, &str)] = &[#(#invariant_entries),*];
        pub static PARITY_MATRIX: &[(&str, &str, Option<&str>)] = &[#(#parity_entries),*];

        impl crate::Molecule {
            #(#methods)*
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use quote::ToTokens;
    use std::collections::BTreeSet;
    use syn::{ImplItem, Item};

    fn context_methods(expanded: proc_macro2::TokenStream, capability: &str) -> BTreeSet<String> {
        let file = syn::parse2::<syn::File>(expanded).expect("generated registry must parse");
        fn collect(items: &[Item], capability: &str, methods: &mut BTreeSet<String>) {
            for item in items {
                match item {
                    Item::Impl(item_impl) => {
                        if item_impl
                            .self_ty
                            .to_token_stream()
                            .to_string()
                            .contains(capability)
                        {
                            methods.extend(item_impl.items.iter().filter_map(|item| match item {
                                ImplItem::Fn(function) => Some(function.sig.ident.to_string()),
                                _ => None,
                            }));
                        }
                    }
                    Item::Mod(item_mod) => {
                        if let Some((_, nested)) = &item_mod.content {
                            collect(nested, capability, methods);
                        }
                    }
                    _ => {}
                }
            }
        }

        let mut methods = BTreeSet::new();
        collect(&file.items, capability, &mut methods);
        methods
    }

    #[test]
    fn registry_generates_compile_time_context_from_declared_access() {
        let registry = syn::parse_str::<Registry>(
            r#"
            op conformer() {
                method: with_conformer,
                impl_fn: conformer_impl,
                output: single,
                domain: coordinate,
                kind: coordinate,
                topology_edit: none,
                access: { read: [topology, properties], write: [coordinates] },
                may_mutate: [coordinates],
                auto_remap: [],
                derived_effects: {
                    recompute: [],
                    preserve: [topology],
                    invalidate: [coordinate_dependent_stereo],
                },
                cip_state: preserve,
                semantic_preconditions: [],
                requires_mapping: none,
                feature: CONFORMER_FEATURE,
                parity: required_now,
                parity_profile: "experiment",
                io_roundtrip: false,
                invariant_profile: "coordinates",
                inplace: false,
            }
            "#,
        )
        .expect("registry input must parse");

        let expanded = expand_registry(registry).expect("registry must expand");
        let expanded_text = expanded.to_string();
        assert!(expanded_text.contains("OpParts < 'a , ConformerAccess >"));
        assert!(!expanded_text.contains("ConformerContext"));
        let methods = context_methods(expanded, "ConformerAccess");
        for expected in [
            "topology",
            "properties",
            "begin_coordinates_mut",
            "commit_coordinates",
        ] {
            assert!(methods.contains(expected), "missing capability {expected}");
        }
        for forbidden in [
            "coordinates",
            "begin_topology_mut",
            "commit_topology",
            "begin_properties_mut",
            "commit_properties",
            "record_topology_mapping",
            "finish",
            "new",
        ] {
            assert!(
                !methods.contains(forbidden),
                "undeclared capability {forbidden} was generated"
            );
        }
    }
}
