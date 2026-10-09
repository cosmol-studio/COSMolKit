//! Compile-time operation and binding declarations for COSMolKit.

#[allow(dead_code)]
mod binding;
mod bio;
#[allow(dead_code)]
mod declaration;
mod javascript_options;
#[allow(dead_code)]
mod matrices;
#[allow(dead_code)]
mod projection;
mod python_configuration;
mod python_enum;
mod result;
mod status;
#[allow(dead_code)]
mod wrappers;

use proc_macro::TokenStream;
use quote::quote;

/// Derive a plain-options constructor and TypeScript options interface from
/// an existing wasm-bindgen constructor. Native constructor validation wins.
#[proc_macro_attribute]
pub fn javascript_options(attribute: TokenStream, item: TokenStream) -> TokenStream {
    let result = syn::parse::<syn::LitStr>(attribute)
        .and_then(|name| syn::parse(item).and_then(|item| javascript_options::expand(item, name)));
    match result {
        Ok(tokens) => tokens.into(),
        Err(error) => error.to_compile_error().into(),
    }
}

/// Accept a declared Python enum member or its exact snake-case string at
/// every native extraction boundary. Outputs remain the original enum class.
#[proc_macro_attribute]
pub fn python_enum(attribute: TokenStream, item: TokenStream) -> TokenStream {
    let existing_methods = attribute.to_string() == "existing_methods";
    if !attribute.is_empty() && !existing_methods {
        return syn::Error::new(
            proc_macro2::Span::call_site(),
            "expected no arguments or existing_methods",
        )
        .to_compile_error()
        .into();
    }
    match syn::parse(item).and_then(|item| python_enum::expand(item, existing_methods)) {
        Ok(tokens) => tokens.into(),
        Err(error) => error.to_compile_error().into(),
    }
}

/// Generate configuration setters from the binding's actual constructor.
/// Each setter validates a replacement through that constructor before commit.
/// `fieldwise` retains preset-only state in records with public Rust fields.
#[proc_macro_attribute]
pub fn python_configuration(attribute: TokenStream, item: TokenStream) -> TokenStream {
    let fieldwise = attribute.to_string() == "fieldwise";
    if !attribute.is_empty() && !fieldwise {
        return syn::Error::new(
            proc_macro2::Span::call_site(),
            "expected no arguments or fieldwise",
        )
        .to_compile_error()
        .into();
    }
    match syn::parse(item).and_then(|item| python_configuration::expand(item, fieldwise)) {
        Ok(tokens) => tokens.into(),
        Err(error) => error.to_compile_error().into(),
    }
}

/// Converts one marked pending molecule field after runtime finalization.
/// This derive is for registry-owned result types in the runtime crate.
#[proc_macro_derive(MoleculeResult, attributes(pending_molecule))]
pub fn molecule_result(input: TokenStream) -> TokenStream {
    match syn::parse(input).and_then(result::expand) {
        Ok(tokens) => tokens.into(),
        Err(error) => error.to_compile_error().into(),
    }
}

/// Injects a single-output operation context into an operation body.
///
/// ```ignore
/// #[mol_op_body(remove_hs, context)]
/// fn remove_hydrogens_impl() -> Result<(), OperationError> {
///     let _ = context;
///     Ok(())
/// }
/// ```
#[proc_macro_attribute]
pub fn mol_op_body(attribute: TokenStream, item: TokenStream) -> TokenStream {
    projection::expand_body(attribute, item, projection::BodyClass::MoleculeSingle)
}

/// Injects a multiple-output operation context into an operation body.
#[proc_macro_attribute]
pub fn mol_multi_op_body(attribute: TokenStream, item: TokenStream) -> TokenStream {
    projection::expand_body(attribute, item, projection::BodyClass::MoleculeMultiple)
}

/// Declares Molecule operation specifications, matrices, and access markers.
#[proc_macro]
pub fn molecule_ops(input: TokenStream) -> TokenStream {
    match syn::parse::<declaration::MoleculeRegistry>(input).and_then(|registry| {
        let markers = projection::expand_molecule_access_markers(&registry)?;
        let matrices = matrices::expand_molecule_matrices(&registry)?;
        let wrappers = wrappers::expand_molecule_wrappers(&registry)?;
        Ok(quote!(#markers #matrices #wrappers))
    }) {
        Ok(output) => output.into(),
        Err(error) => error.to_compile_error().into(),
    }
}

/// Generates lightweight BIO value/in-place pairs, field access and metadata.
#[proc_macro]
pub fn bio_structure_ops(input: TokenStream) -> TokenStream {
    match syn::parse::<bio::Registry>(input).map(bio::expand) {
        Ok(output) => output.into(),
        Err(error) => error.to_compile_error().into(),
    }
}

/// Declares the canonical cross-language naming and type contract.
#[proc_macro]
pub fn binding_contract(input: TokenStream) -> TokenStream {
    match binding::expand_binding_contract(input.into()) {
        Ok(output) => output.into(),
        Err(error) => error.to_compile_error().into(),
    }
}
