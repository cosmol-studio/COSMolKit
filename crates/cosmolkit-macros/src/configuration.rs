//! Native declarations supply field names and Rust-resolved types, never a second schema.
use proc_macro2::TokenStream;
use quote::{format_ident, quote};
use syn::{Attribute, Data, DeriveInput, Fields, Ident, Path, Token, Type, Visibility, braced};

pub(crate) fn derive(input: DeriveInput) -> syn::Result<TokenStream> {
    let mut native_storage = false;
    let mut required = false;
    let mut factory = false;
    let mut stateful = false;
    let mut accessors = false;
    for attribute in &input.attrs {
        if attribute.path().is_ident("configuration") {
            attribute.parse_nested_meta(|meta| {
                if meta.path.is_ident("scalar") {
                    native_storage = true;
                    Ok(())
                } else if meta.path.is_ident("required") {
                    required = true;
                    Ok(())
                } else if meta.path.is_ident("factory") {
                    factory = true;
                    Ok(())
                } else if meta.path.is_ident("stateful") {
                    stateful = true;
                    Ok(())
                } else if meta.path.is_ident("accessors") {
                    accessors = true;
                    Ok(())
                } else {
                    Err(meta.error("expected scalar storage for immutable input fields"))
                }
            })?;
        }
    }
    if !input.generics.params.is_empty() || input.generics.where_clause.is_some() {
        return Err(syn::Error::new_spanned(
            input.generics,
            "Configuration requires a concrete input record",
        ));
    }
    let Data::Struct(data) = input.data else {
        return Err(syn::Error::new_spanned(
            input.ident,
            "Configuration requires a named-field input record",
        ));
    };
    let Fields::Named(fields) = data.fields else {
        return Err(syn::Error::new_spanned(
            input.ident,
            "Configuration requires named public input fields",
        ));
    };
    if fields.named.is_empty() {
        return Err(syn::Error::new_spanned(
            input.ident,
            "Configuration requires at least one input field",
        ));
    }
    let metadata = format_ident!("__configuration_{}", input.ident);
    let exporter = format_ident!("__ck_configuration_{}", input.ident);
    let mut names = Vec::new();
    let mut aliases = Vec::new();
    let mut types = Vec::new();
    let mut attributes = Vec::new();
    let mut flattened_exports = Vec::new();
    for (index, field) in fields.named.into_iter().enumerate() {
        let policy = field.attrs.iter().find_map(|a| {
            a.path()
                .is_ident("configuration")
                .then(|| a.meta.require_list().ok().map(|m| m.tokens.to_string()))
                .flatten()
                .filter(|value| matches!(value.as_str(), "state" | "adapter"))
        });
        if let Some(policy) = policy {
            if !stateful {
                return Err(syn::Error::new_spanned(
                    field,
                    "state/adapter fields require an explicit stateful projection",
                ));
            }
            // Runtime state is not editable input. Adapter inputs retain their
            // declared type and marker for language-specific ownership adapters.
            if policy == "state" {
                continue;
            }
        }
        if !matches!(field.vis, Visibility::Public(_)) && !accessors {
            return Err(syn::Error::new_spanned(
                field,
                "Configuration cannot project private or stateful fields",
            ));
        }
        if !field.attrs.is_empty() {
            // Conditional fields and custom construction policies need an explicit protocol.
            for attribute in &field.attrs {
                if !attribute.path().is_ident("doc") && !attribute.path().is_ident("configuration")
                {
                    return Err(syn::Error::new_spanned(
                        attribute,
                        "Configuration field attributes require an explicit supported projection protocol",
                    ));
                }
            }
        }
        names.push(field.ident.unwrap());
        let mut attrs = field.attrs;
        if attrs.iter().any(|a| {
            a.path().is_ident("configuration")
                && a.meta
                    .require_list()
                    .is_ok_and(|m| m.tokens.to_string() == "flatten")
        }) {
            let Type::Path(ty) = &field.ty else {
                return Err(syn::Error::new_spanned(
                    &field.ty,
                    "flatten requires a native configuration type",
                ));
            };
            let mut exporter = ty.path.clone();
            let last = exporter.segments.last_mut().unwrap();
            last.ident = format_ident!("__ck_configuration_{}", last.ident);
            if exporter.segments.len() == 1 {
                exporter = syn::parse_quote!(crate::#exporter);
            }
            let alias = format_ident!("Flatten{index}");
            flattened_exports.push(quote!(pub use #exporter as #alias;));
        }
        if required {
            attrs.push(syn::parse_quote!(#[configuration(required)]));
        }
        attributes.push(attrs);
        aliases.push(format_ident!("Field{index}"));
        types.push(field.ty);
    }
    let marker = Ident::new(
        if stateful || accessors {
            "adapter"
        } else {
            match (native_storage, factory) {
                (false, false) => "fields",
                (true, false) => "native",
                (false, true) => "fields_factory",
                (true, true) => "native_factory",
            }
        },
        input.ident.span(),
    );
    Ok(quote! {
        #[doc(hidden)]
        #[allow(non_snake_case)]
        pub mod #metadata {
            #[allow(unused_imports)]
            use super::*;
            #(pub type #aliases = #types;)*
            #(#flattened_exports)*
        }
        #[doc(hidden)]
        #[macro_export]
        macro_rules! #exporter {
            ($consumer:path, $native:path, using $projection:path) => {
                $consumer! { @#marker $native { #(#(#attributes)* #names: $crate::#metadata::#aliases),* } methods {} using $projection }
            };
            ($consumer:path, $native:path, {$($method:item)*}, {$($context:tt)*}) => {
                $consumer! { @resume {$($context)*} @#marker $native { #(#(#attributes)* #names: $crate::#metadata::#aliases),* } methods {$($method)*} }
            };
            ($consumer:path, $native:path, {$($method:item)*}) => {
                $consumer! { @#marker $native { #(#(#attributes)* #names: $crate::#metadata::#aliases),* } methods {$($method)*} }
            };
            ($consumer:path, $native:path) => {
                $crate::#exporter!($consumer, $native, {});
            };
        }
    })
}

struct Field {
    attrs: Vec<Attribute>,
    name: Ident,
    ty: Type,
    children: Option<Vec<Field>>,
}
impl syn::parse::Parse for Field {
    fn parse(input: syn::parse::ParseStream<'_>) -> syn::Result<Self> {
        let attrs = input.call(Attribute::parse_outer)?;
        let name = input.parse()?;
        input.parse::<Token![:]>()?;
        let ty = input.parse()?;
        let children = if input.peek(Token![=>]) {
            input.parse::<Token![=>]>()?;
            let content;
            braced!(content in input);
            Some(
                content
                    .parse_terminated(Field::parse, Token![,])?
                    .into_iter()
                    .collect(),
            )
        } else {
            None
        };
        Ok(Self {
            attrs,
            name,
            ty,
            children,
        })
    }
}
impl quote::ToTokens for Field {
    fn to_tokens(&self, tokens: &mut TokenStream) {
        let Self {
            attrs,
            name,
            ty,
            children,
        } = self;
        let nested = children.as_ref().map(|fields| quote!(=> {#(#fields),*}));
        tokens.extend(quote!(#(#attrs)* #name: #ty #nested));
    }
}
struct Projection {
    marker: Ident,
    native: Path,
    fields: Vec<Field>,
    methods: TokenStream,
    consumer: Option<Path>,
}
impl syn::parse::Parse for Projection {
    fn parse(input: syn::parse::ParseStream<'_>) -> syn::Result<Self> {
        input.parse::<Token![@]>()?;
        let marker: Ident = input.parse()?;
        if ![
            "fields",
            "native",
            "fields_factory",
            "native_factory",
            "adapter",
        ]
        .iter()
        .any(|name| marker == name)
        {
            return Err(syn::Error::new_spanned(
                marker,
                "expected native Configuration field metadata",
            ));
        }
        let native = input.parse()?;
        let body;
        braced!(body in input);
        let fields = body
            .parse_terminated(Field::parse, Token![,])?
            .into_iter()
            .collect();
        let methods = parse_methods(input)?;
        let consumer = if input.is_empty() {
            None
        } else {
            let marker: Ident = input.parse()?;
            if marker != "using" {
                return Err(syn::Error::new_spanned(
                    marker,
                    "expected explicit state adapter",
                ));
            }
            Some(input.parse()?)
        };
        Ok(Self {
            marker,
            native,
            fields,
            methods,
            consumer,
        })
    }
}

fn pending_flatten(fields: &mut [Field]) -> Option<&mut Field> {
    for field in fields {
        if field.children.is_some() {
            if let Some(pending) = pending_flatten(field.children.as_mut().unwrap()) {
                return Some(pending);
            }
        } else if field.attrs.iter().any(|a| {
            a.path().is_ident("configuration")
                && a.meta
                    .require_list()
                    .is_ok_and(|m| m.tokens.to_string() == "flatten")
        }) {
            return Some(field);
        }
    }
    None
}

struct Resume {
    parent: Projection,
    child: Projection,
}
impl syn::parse::Parse for Resume {
    fn parse(input: syn::parse::ParseStream<'_>) -> syn::Result<Self> {
        input.parse::<Token![@]>()?;
        let marker: Ident = input.parse()?;
        if marker != "resume" {
            return Err(syn::Error::new_spanned(marker, "expected resume"));
        }
        let content;
        braced!(content in input);
        Ok(Self {
            parent: content.parse()?,
            child: input.parse()?,
        })
    }
}

fn collect_leaves<'a>(
    fields: &'a [Field],
    prefix: TokenStream,
    leaves: &mut Vec<(&'a Field, TokenStream)>,
) {
    for field in fields {
        let name = &field.name;
        let path = if prefix.is_empty() {
            quote!(#name)
        } else {
            quote!(#prefix.#name)
        };
        if let Some(children) = &field.children {
            collect_leaves(children, path, leaves);
        } else {
            leaves.push((field, path));
        }
    }
}

fn record_literal(
    fields: &[Field],
    values: &std::collections::HashMap<String, TokenStream>,
) -> TokenStream {
    let assignments = fields.iter().map(|field| {
        let name = &field.name;
        let value = if let Some(children) = &field.children {
            let ty = &field.ty;
            let nested = record_literal(children, values);
            quote!({ use #ty as NativeGroup; NativeGroup { #nested } })
        } else {
            values[&name.to_string()].clone()
        };
        quote!(#name: #value)
    });
    quote!(#(#assignments,)*)
}

fn parse_methods(input: syn::parse::ParseStream<'_>) -> syn::Result<TokenStream> {
    if input.is_empty() {
        return Ok(TokenStream::new());
    }
    if input.peek(Token![,]) {
        input.parse::<Token![,]>()?;
    }
    let name: Ident = input.parse()?;
    if name != "methods" {
        return Err(syn::Error::new_spanned(name, "expected methods"));
    }
    let content;
    braced!(content in input);
    content.parse()
}

struct Request {
    native: Path,
    methods: TokenStream,
    consumer: Option<Path>,
}
impl syn::parse::Parse for Request {
    fn parse(input: syn::parse::ParseStream<'_>) -> syn::Result<Self> {
        let native = input.parse()?;
        let ahead = input.fork();
        let custom = if ahead.peek(Token![,]) {
            ahead.parse::<Token![,]>()?;
            ahead.peek(Ident) && ahead.parse::<Ident>()? == "using"
        } else {
            false
        };
        if custom {
            input.parse::<Token![,]>()?;
            input.parse::<Ident>()?;
            Ok(Self {
                native,
                methods: TokenStream::new(),
                consumer: Some(input.parse()?),
            })
        } else {
            Ok(Self {
                native,
                methods: parse_methods(input)?,
                consumer: None,
            })
        }
    }
}

fn emit_projection(mut projection: Projection) -> syn::Result<TokenStream> {
    if let Some(field) = pending_flatten(&mut projection.fields) {
        let Type::Path(ty) = &field.ty else {
            return Err(syn::Error::new_spanned(
                &field.ty,
                "expected native metadata alias",
            ));
        };
        let mut exporter = ty.path.clone();
        let last = exporter.segments.last_mut().unwrap();
        last.ident = format_ident!(
            "Flatten{}",
            last.ident
                .to_string()
                .strip_prefix("Field")
                .ok_or_else(|| syn::Error::new_spanned(&field.ty, "expected native field alias"))?
        );
        let child = field.ty.clone();
        let native = &projection.native;
        let marker = &projection.marker;
        let fields = &projection.fields;
        let methods = &projection.methods;
        let consumer = projection.consumer.as_ref().map(|path| quote!(using #path));
        return Ok(
            quote! { #exporter!(::cosmolkit::project_configuration, #child, {}, { @#marker #native {#(#fields),*} methods {#methods} #consumer }); },
        );
    }
    let Projection {
        marker,
        native,
        fields,
        methods,
        consumer,
    } = projection;
    if let Some(consumer) = consumer {
        let mut leaves = Vec::new();
        collect_leaves(&fields, TokenStream::new(), &mut leaves);
        // Explicit Python ownership adapters generate their ordinary properties
        // from this callback and implement marked fields themselves. Keep that
        // protocol unchanged; the raw metadata export includes adapter fields.
        leaves.retain(|(field, _)| {
            !field.attrs.iter().any(|attribute| {
                attribute.path().is_ident("configuration")
                    && attribute
                        .meta
                        .require_list()
                        .is_ok_and(|meta| meta.tokens.to_string() == "adapter")
            })
        });
        let names = leaves.iter().map(|(field, _)| &field.name);
        let types: Vec<_> = leaves.iter().map(|(field, _)| &field.ty).collect();
        return Ok(quote! {
            #consumer! {#(#names: #types),*}
        });
    }
    if marker == "adapter" {
        return Err(syn::Error::new_spanned(
            native,
            "stateful/accessor configuration requires its explicit ownership adapter",
        ));
    }
    let name = &native.segments.last().unwrap().ident;
    let factory = marker == "fields_factory" || marker == "native_factory";
    let mut arguments = Vec::new();
    let mut signatures = Vec::new();
    let mut inputs = std::collections::HashMap::new();
    let mut factory_inputs = Vec::new();
    let mut accessors = Vec::new();
    let mut required_arguments = Vec::new();
    let mut required_signatures = Vec::new();
    let mut leaves = Vec::new();
    collect_leaves(&fields, TokenStream::new(), &mut leaves);
    for (definition, native_field) in leaves {
        let Field {
            attrs,
            name: field,
            ty,
            ..
        } = definition;
        if inputs.contains_key(&field.to_string()) {
            return Err(syn::Error::new_spanned(
                field,
                "flattened configuration field names must be unique",
            ));
        }
        let mut required = false;
        let mut none_default = false;
        let mut docs = Vec::new();
        for attr in attrs {
            if attr.path().is_ident("doc") {
                docs.push(attr);
            } else if attr.path().is_ident("configuration") {
                attr.parse_nested_meta(|meta| {
                    if meta.path.is_ident("required") {
                        required = true;
                        Ok(())
                    } else if meta.path.is_ident("none_default") {
                        none_default = true;
                        Ok(())
                    } else {
                        Err(meta.error("expected required or none_default"))
                    }
                })?;
            }
        }
        if required && none_default {
            return Err(syn::Error::new_spanned(
                field,
                "required input cannot have a None default",
            ));
        }
        let adapter = quote!(<#ty as crate::configuration::FieldProjection>::Adapter);
        let input = if none_default {
            quote!(crate::parameters::FieldInput<#adapter>)
        } else {
            quote!(crate::parameters::ContainerInput<#adapter>)
        };
        let value = if none_default {
            quote!(value.or_default(|| <#native>::default().#native_field))
        } else {
            quote!(value.0)
        };
        let ctor_value = if none_default {
            quote!(#field.or_default(|| <#native>::default().#native_field))
        } else {
            quote!(#field.0)
        };
        if required {
            required_arguments.push(quote!(#field: #input));
            required_signatures.push(quote!(#field));
        } else {
            arguments.push(quote!(#field: #input));
            let default = if none_default {
                quote!(crate::parameters::FieldInput::<#adapter>::Default)
            } else {
                quote!(crate::parameters::ContainerInput::<#adapter>(<#native>::default().#native_field))
            };
            signatures.push(quote!(#field=#default));
        }
        inputs.insert(field.to_string(), ctor_value.clone());
        factory_inputs.push(ctor_value);
        let setter = format_ident!("set_{}", field);
        accessors.push(quote! {
            #(#docs)*
            #[getter]
            fn #field(slf: pyo3::PyRef<'_, Self>) -> crate::parameters::FieldOutput<#adapter, Self> {
                crate::parameters::FieldOutput::new(
                    slf.into(), |native: &#native| &native.#native_field,
                    |native: &mut #native| &mut native.#native_field,
                )
            }
            #[setter]
            fn #setter(slf: &pyo3::Bound<'_, Self>, value: #input) -> pyo3::PyResult<()> {
                let value = #value;
                slf.try_borrow_mut()?.inner.write(slf.py(), |native| native.#native_field = value)
            }
        });
    }
    let inputs = record_literal(&fields, &inputs);
    let factory_method = if factory {
        quote! {
            #[staticmethod]
            #[pyo3(name = "new", signature = (#(#required_signatures,)* #(#signatures,)*))]
            fn registered_new(py: pyo3::Python<'_>, #(#required_arguments,)* #(#arguments,)*) -> pyo3::PyResult<Self> {
                Self::from_rust(py, <#native>::new(#(#factory_inputs),*))
            }
        }
    } else {
        quote!()
    };
    Ok(quote! {
        #[pyo3::pyclass(module = "cosmolkit")]
        pub struct #name { pub(crate) inner: crate::parameters::NativeField<#native> }
        impl #name {
            pub(crate) fn from_rust(_py: pyo3::Python<'_>, value: #native) -> pyo3::PyResult<Self> { Ok(Self { inner: crate::parameters::NativeField::Owned(value) }) }
            pub(crate) fn to_rust(&self, py: pyo3::Python<'_>) -> pyo3::PyResult<#native> { self.inner.read(py, Clone::clone) }
        }
        impl crate::parameters::NativeOwner for #name {
            type Value = #native;
            fn read_native(&self, py: pyo3::Python<'_>, visit: &mut dyn FnMut(&Self::Value)) -> pyo3::PyResult<()> { self.inner.read(py, visit) }
            fn write_native(&mut self, py: pyo3::Python<'_>, visit: &mut dyn FnMut(&mut Self::Value)) -> pyo3::PyResult<()> { self.inner.write(py, visit) }
        }
        impl crate::parameters::ConfigurationView for #name {
            fn from_field(inner: crate::parameters::NativeField<#native>) -> Self { Self { inner } }
        }
        impl crate::parameters::Configuration for #name {
            type Native = #native;
            fn from_native(py: pyo3::Python<'_>, value: Self::Native) -> pyo3::PyResult<Self> { Self::from_rust(py, value) }
            fn to_native(&self, py: pyo3::Python<'_>) -> pyo3::PyResult<Self::Native> { self.to_rust(py) }
        }
        impl crate::configuration::FieldProjection for #native {
            type Adapter = crate::parameters::LiveConfigurationField<#name>;
        }
        #[pyo3::pymethods]
        impl #name {
            #[new]
            #[pyo3(signature = (*, #(#required_signatures,)* #(#signatures,)*))]
            fn new(py: pyo3::Python<'_>, #(#required_arguments,)* #(#arguments,)*) -> pyo3::PyResult<Self> {
                use #native as Native;
                Self::from_rust(py, Native { #inputs })
            }
            #(#accessors)*
            fn __traverse__(&self, visit: pyo3::gc::PyVisit<'_>) -> Result<(), pyo3::PyTraverseError> {
                self.inner.traverse(visit)
            }
            #factory_method
            #methods
        }
    })
}

pub(crate) fn project(input: TokenStream) -> syn::Result<TokenStream> {
    if input
        .clone()
        .into_iter()
        .next()
        .is_some_and(|token| token.to_string() == "@")
    {
        if input
            .clone()
            .into_iter()
            .nth(1)
            .is_some_and(|token| token.to_string() == "resume")
        {
            let Resume { mut parent, child } = syn::parse2(input)?;
            let field = pending_flatten(&mut parent.fields).ok_or_else(|| {
                syn::Error::new_spanned(&parent.native, "no pending native flatten field")
            })?;
            field.children = Some(child.fields);
            return emit_projection(parent);
        }
        return emit_projection(syn::parse2(input)?);
    }
    let Request {
        native,
        methods,
        consumer,
    } = syn::parse2(input)?;
    if native
        .segments
        .iter()
        .any(|segment| !segment.arguments.is_empty())
    {
        return Err(syn::Error::new_spanned(
            native,
            "project_configuration requires a concrete native type path",
        ));
    }
    let mut exporter = native.clone();
    let last = exporter.segments.last_mut().unwrap();
    last.ident = format_ident!("__ck_configuration_{}", last.ident);
    Ok(if let Some(consumer) = consumer {
        quote! { #exporter!(::cosmolkit::project_configuration, #native, using #consumer); }
    } else {
        quote! { #exporter!(::cosmolkit::project_configuration, #native, {#methods}); }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adapter_metadata_keeps_types_without_rebuilding_state() {
        let expanded = derive(syn::parse_quote! {
            #[configuration(stateful)]
            struct Options {
                pub limit: usize,
                #[configuration(adapter)]
                pub callback: Option<Callback>,
                #[configuration(state)]
                pub cache: Cache,
            }
        })
        .unwrap()
        .to_string();
        assert!(expanded.contains("type Field1 = Option < Callback >"));
        assert!(expanded.contains("configuration (adapter)"));
        assert!(!expanded.contains("Field2"));
        assert!(!expanded.contains("cache :"));
        let projected = emit_projection(syn::parse_quote! {
            @adapter Options {
                limit: usize,
                #[configuration(adapter)]
                callback: Option<Callback>
            } methods {} using scalar_properties
        })
        .unwrap()
        .to_string();
        assert_eq!(projected, "scalar_properties ! { limit : usize }");
    }

    #[test]
    fn rejects_non_input_shapes_instead_of_guessing_a_projection() {
        for input in [
            syn::parse_quote!(
                struct Private {
                    state: u32,
                }
            ),
            syn::parse_quote!(
                struct Generic<T> {
                    pub value: T,
                }
            ),
            syn::parse_quote!(
                struct Tuple(pub u32);
            ),
            syn::parse_quote!(
                enum Choice {
                    First,
                    Second,
                }
            ),
            syn::parse_quote!(
                struct Conditional {
                    #[cfg(feature = "extra")]
                    pub value: u32,
                }
            ),
        ] {
            assert!(derive(input).is_err());
        }
    }
}
