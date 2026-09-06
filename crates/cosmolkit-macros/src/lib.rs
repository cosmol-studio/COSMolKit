//! Operation-body macros for the new `cosmolkit` runtime.
//!
//! This crate intentionally has a much smaller scope than the legacy
//! `cosmolkit-macros` crate. It only projects a registry operation name into
//! the corresponding compile-time capability context. Registry generation and
//! source-backed operation implementation will migrate independently.

use proc_macro::TokenStream;
use quote::{format_ident, quote};
use syn::{
    FnArg, Ident, ItemFn, LitStr, Token, Type, Visibility, braced, bracketed, parse::Parse,
    parse::ParseStream, parse_macro_input, parse_quote,
};

struct BodyAttribute {
    operation: Ident,
    _comma: Token![,],
    context: Ident,
}

impl Parse for BodyAttribute {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        Ok(Self {
            operation: input.parse()?,
            _comma: input.parse()?,
            context: input.parse()?,
        })
    }
}

fn access_marker(operation: &Ident) -> Ident {
    let mut name = String::new();
    for component in operation.to_string().split('_') {
        let mut chars = component.chars();
        if let Some(first) = chars.next() {
            name.extend(first.to_uppercase());
            name.extend(chars);
        }
    }
    format_ident!("{name}Access")
}

fn expand(attribute: TokenStream, item: TokenStream, multiple: bool) -> TokenStream {
    let BodyAttribute {
        operation, context, ..
    } = parse_macro_input!(attribute as BodyAttribute);
    let mut function = parse_macro_input!(item as ItemFn);

    for input in &function.sig.inputs {
        if matches!(input, FnArg::Receiver(_)) {
            return syn::Error::new_spanned(input, "operation bodies must not receive self")
                .to_compile_error()
                .into();
        }
    }

    let access = access_marker(&operation);
    let context_type: Type = if multiple {
        parse_quote!(&mut crate::MultiOutputOpParts<'_, crate::#access>)
    } else {
        parse_quote!(&mut crate::OpParts<'_, crate::#access>)
    };
    let parameters = std::mem::take(&mut function.sig.inputs);
    function
        .sig
        .inputs
        .push(parse_quote!(#context: #context_type));
    function.sig.inputs.extend(parameters);

    quote!(#function).into()
}

/// Injects a single-output operation context into an operation body.
///
/// ```ignore
/// #[mol_op_body(remove_hydrogens, context)]
/// fn remove_hydrogens_impl() -> Result<(), OperationError> {
///     let _ = context;
///     Ok(())
/// }
/// ```
#[proc_macro_attribute]
pub fn mol_op_body(attribute: TokenStream, item: TokenStream) -> TokenStream {
    expand(attribute, item, false)
}

/// Injects a multiple-output operation context into an operation body.
#[proc_macro_attribute]
pub fn mol_multi_op_body(attribute: TokenStream, item: TokenStream) -> TokenStream {
    expand(attribute, item, true)
}

struct BindingRegistry {
    visibility: Visibility,
    name: Ident,
    entries: Vec<BindingEntry>,
}

struct BindingEntry {
    semantic_id: LitStr,
    owner: Ident,
    kind: Ident,
    rust_method: Ident,
    python_name: LitStr,
    javascript_name: LitStr,
    feature: LitStr,
    return_kind: Ident,
    state_model: Ident,
    signature: Ident,
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
            let entry_content;
            braced!(entry_content in content);
            entries.push(parse_binding_entry(&entry_content)?);
            if content.peek(Token![,]) {
                content.parse::<Token![,]>()?;
            }
        }
        if input.peek(Token![;]) {
            input.parse::<Token![;]>()?;
        }
        Ok(Self {
            visibility,
            name,
            entries,
        })
    }
}

fn parse_binding_entry(input: ParseStream<'_>) -> syn::Result<BindingEntry> {
    let mut semantic_id = None;
    let mut owner = None;
    let mut kind = None;
    let mut rust_method = None;
    let mut python_name = None;
    let mut javascript_name = None;
    let mut feature = None;
    let mut return_kind = None;
    let mut state_model = None;
    let mut signature = None;
    while !input.is_empty() {
        let key: Ident = input.parse()?;
        input.parse::<Token![:]>()?;
        match key.to_string().as_str() {
            "semantic_id" => semantic_id = Some(input.parse()?),
            "owner" => owner = Some(input.parse()?),
            "kind" => kind = Some(input.parse()?),
            "rust" => rust_method = Some(input.parse()?),
            "python" => python_name = Some(input.parse()?),
            "javascript" => javascript_name = Some(input.parse()?),
            "feature" => feature = Some(input.parse()?),
            "returns" => return_kind = Some(input.parse()?),
            "state" => state_model = Some(input.parse()?),
            "signature" => signature = Some(input.parse()?),
            other => {
                return Err(syn::Error::new(
                    key.span(),
                    format!("unknown binding field `{other}`"),
                ));
            }
        }
        if input.peek(Token![,]) {
            input.parse::<Token![,]>()?;
        }
    }
    Ok(BindingEntry {
        semantic_id: required(semantic_id, "semantic_id")?,
        owner: required(owner, "owner")?,
        kind: required(kind, "kind")?,
        rust_method: required(rust_method, "rust")?,
        python_name: required(python_name, "python")?,
        javascript_name: required(javascript_name, "javascript")?,
        feature: required(feature, "feature")?,
        return_kind: required(return_kind, "returns")?,
        state_model: required(state_model, "state")?,
        signature: required(signature, "signature")?,
    })
}

fn required<T>(value: Option<T>, field: &str) -> syn::Result<T> {
    value.ok_or_else(|| {
        syn::Error::new(
            proc_macro2::Span::call_site(),
            format!("binding entry is missing `{field}`"),
        )
    })
}

fn owner(value: &Ident) -> syn::Result<proc_macro2::TokenStream> {
    match value.to_string().as_str() {
        "molecule" => Ok(quote!(crate::BindingOwner::Molecule)),
        "module" => Ok(quote!(crate::BindingOwner::Module)),
        _ => Err(syn::Error::new_spanned(
            value,
            "owner must be `molecule` or `module`",
        )),
    }
}

fn kind(value: &Ident) -> syn::Result<proc_macro2::TokenStream> {
    match value.to_string().as_str() {
        "instance" => Ok(quote!(crate::BindingKind::Instance)),
        "static" => Ok(quote!(crate::BindingKind::Static)),
        "module" => Ok(quote!(crate::BindingKind::Module)),
        _ => Err(syn::Error::new_spanned(
            value,
            "kind must be `instance`, `static`, or `module`",
        )),
    }
}

fn return_kind(value: &Ident) -> syn::Result<proc_macro2::TokenStream> {
    let variant = match value.to_string().as_str() {
        "result_molecule" => quote!(ResultMolecule),
        "result_unit" => quote!(ResultUnit),
        "result_bool" => quote!(ResultBool),
        "result_f64" => quote!(ResultF64),
        "result_string" => quote!(ResultString),
        "string" => quote!(String),
        "unit" => quote!(Unit),
        _ => {
            return Err(syn::Error::new_spanned(
                value,
                "unsupported binding return kind",
            ));
        }
    };
    Ok(quote!(crate::ReturnKind::#variant))
}

fn state_model(value: &Ident) -> syn::Result<proc_macro2::TokenStream> {
    let variant = match value.to_string().as_str() {
        "value_returning" => quote!(ValueReturning),
        "in_place" => quote!(InPlace),
        "read_only" => quote!(ReadOnly),
        _ => {
            return Err(syn::Error::new_spanned(
                value,
                "unsupported binding state model",
            ));
        }
    };
    Ok(quote!(crate::StateModel::#variant))
}

fn signature_assertion(
    index: usize,
    method: &Ident,
    signature: &Ident,
) -> syn::Result<proc_macro2::TokenStream> {
    let assertion = format_ident!("__BINDING_SIGNATURE_ASSERT_{index}");
    let value = match signature.to_string().as_str() {
        "value_molecule" => {
            quote!(const #assertion: fn(&crate::Molecule) -> Result<crate::Molecule, crate::OperationError> = crate::Molecule::#method;)
        }
        "in_place_unit" => {
            quote!(const #assertion: fn(&mut crate::Molecule) -> Result<(), crate::OperationError> = crate::Molecule::#method;)
        }
        "read_only_f64" => {
            quote!(const #assertion: fn(&crate::Molecule) -> Result<f64, crate::OperationError> = crate::Molecule::#method;)
        }
        "read_only_string" => {
            quote!(const #assertion: fn(&crate::Molecule) -> Result<String, crate::OperationError> = crate::Molecule::#method;)
        }
        _ => {
            return Err(syn::Error::new_spanned(
                signature,
                "unsupported binding Rust signature",
            ));
        }
    };
    Ok(value)
}

fn expand_binding_registry(input: BindingRegistry) -> syn::Result<proc_macro2::TokenStream> {
    let BindingRegistry {
        visibility,
        name,
        entries,
    } = input;
    let mut values = Vec::with_capacity(entries.len());
    let mut assertions = Vec::with_capacity(entries.len());
    for (index, entry) in entries.iter().enumerate() {
        let BindingEntry {
            semantic_id,
            owner: owner_name,
            kind: kind_name,
            rust_method,
            python_name,
            javascript_name,
            feature,
            return_kind: return_name,
            state_model: state_name,
            signature,
        } = entry;
        let owner = owner(owner_name)?;
        let kind = kind(kind_name)?;
        let return_kind = return_kind(return_name)?;
        let state_model = state_model(state_name)?;
        let assertion = signature_assertion(index, rust_method, signature)?;
        values.push(quote! {
            crate::BindingContractEntry {
                semantic_id: #semantic_id,
                owner: #owner,
                kind: #kind,
                rust_name: stringify!(#rust_method),
                python_name: #python_name,
                javascript_name: #javascript_name,
                feature: #feature,
                return_kind: #return_kind,
                state_model: #state_model,
            }
        });
        assertions.push(assertion);
    }
    Ok(quote! {
        #visibility static #name: &[crate::BindingContractEntry] = &[#(#values),*];
        #(#assertions)*
    })
}

/// Declares the canonical cross-language naming and type contract.
#[proc_macro]
pub fn binding_contract(input: TokenStream) -> TokenStream {
    match syn::parse::<BindingRegistry>(input).and_then(expand_binding_registry) {
        Ok(output) => output.into(),
        Err(error) => error.to_compile_error().into(),
    }
}

#[cfg(test)]
mod tests {
    use super::access_marker;
    use syn::parse_quote;

    #[test]
    fn operation_names_map_to_capability_markers() {
        assert_eq!(
            access_marker(&parse_quote!(remove_hydrogens)),
            "RemoveHydrogensAccess"
        );
        assert_eq!(
            access_marker(&parse_quote!(with_3d_conformer)),
            "With3dConformerAccess"
        );
    }
}
