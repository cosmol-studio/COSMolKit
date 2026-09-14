//! Compile-time operation and binding declarations for COSMolKit.

#[allow(dead_code)]
mod binding;
#[allow(dead_code)]
mod declaration;
#[allow(dead_code)]
mod matrices;
#[allow(dead_code)]
mod projection;
mod result;
#[allow(dead_code)]
mod wrappers;

use proc_macro::TokenStream;
use quote::quote;

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
/// #[mol_op_body(remove_hydrogens, context)]
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

/// Injects a BioStructure operation context into an operation body.
#[proc_macro_attribute]
pub fn bio_op_body(attribute: TokenStream, item: TokenStream) -> TokenStream {
    projection::expand_body(attribute, item, projection::BodyClass::Bio)
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

/// Declares BioStructure operation specifications, matrices, and access markers.
#[proc_macro]
pub fn bio_structure_ops(input: TokenStream) -> TokenStream {
    match syn::parse::<declaration::BioRegistry>(input).and_then(|registry| {
        let markers = projection::expand_bio_access_markers(&registry)?;
        let matrices = matrices::expand_bio_matrices(&registry)?;
        let wrappers = wrappers::expand_bio_wrappers(&registry)?;
        Ok(quote!(#markers #matrices #wrappers))
    }) {
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
