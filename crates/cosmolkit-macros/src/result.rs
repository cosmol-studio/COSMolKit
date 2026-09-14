//! Derives structural result conversion; no domain callback or runtime state.

use quote::{format_ident, quote};
use syn::{Data, DeriveInput, Fields, GenericParam, PathArguments, Type};

pub(crate) fn expand(input: DeriveInput) -> syn::Result<proc_macro2::TokenStream> {
    let name = &input.ident;
    let fail = |message: &str| syn::Error::new_spanned(name, message);
    if input.generics.params.len() != 1 || input.generics.where_clause.is_some() {
        return Err(fail(
            "MoleculeResult requires exactly one unbounded type parameter",
        ));
    }
    let Some(GenericParam::Type(parameter)) = input.generics.params.first() else {
        return Err(fail("MoleculeResult requires a molecule type parameter"));
    };
    if !parameter.bounds.is_empty() || parameter.default.is_none() {
        return Err(fail(
            "MoleculeResult requires an unbounded parameter defaulting to Molecule",
        ));
    }
    let molecule_type = &parameter.default.as_ref().unwrap().1;
    let Type::Path(default_path) = molecule_type else {
        return Err(fail("MoleculeResult default must be Molecule"));
    };
    if default_path
        .path
        .segments
        .last()
        .is_none_or(|s| s.ident != "Molecule")
    {
        return Err(fail("MoleculeResult default must be Molecule"));
    }
    let Data::Struct(data) = &input.data else {
        return Err(fail("MoleculeResult requires a named-field struct"));
    };
    let Fields::Named(fields) = &data.fields else {
        return Err(fail("MoleculeResult requires a named-field struct"));
    };
    let marked: Vec<_> = fields
        .named
        .iter()
        .filter(|f| {
            f.attrs
                .iter()
                .any(|a| a.path().is_ident("pending_molecule"))
        })
        .collect();
    if marked.len() != 1 {
        return Err(fail(
            "MoleculeResult requires exactly one #[pending_molecule] field",
        ));
    }
    let generic = &parameter.ident;
    let is_parameter = |ty: &Type| matches!(ty, Type::Path(path) if path.qself.is_none() && path.path.is_ident(generic));
    let field = marked[0];
    let optional = match &field.ty {
        ty if is_parameter(ty) => false,
        Type::Path(path)
            if path.qself.is_none()
                && path.path.segments.len() == 1
                && path.path.segments[0].ident == "Option" =>
        {
            let PathArguments::AngleBracketed(args) = &path.path.segments[0].arguments else {
                return Err(fail("pending field must be M or Option<M>"));
            };
            if args.args.len() != 1
                || !matches!(args.args.first(), Some(syn::GenericArgument::Type(ty)) if is_parameter(ty))
            {
                return Err(fail("pending field must be M or Option<M>"));
            }
            true
        }
        _ => return Err(fail("pending field must be M or Option<M>")),
    };
    let names: Vec<_> = fields
        .named
        .iter()
        .map(|f| f.ident.as_ref().unwrap())
        .collect();
    let locals: Vec<_> = (0..names.len())
        .map(|index| format_ident!("__cosmolkit_result_field_{index}"))
        .collect();
    let pending_index = names
        .iter()
        .position(|name| Some(*name) == field.ident.as_ref())
        .unwrap();
    let pending = &locals[pending_index];
    let conversion = if optional {
        quote!(#pending.map(|pending| finalizer.resolve(pending)).transpose()?)
    } else {
        quote!(finalizer.resolve(#pending)?)
    };
    let assignments = names
        .iter()
        .zip(&locals)
        .enumerate()
        .map(|(index, (field, local))| {
            if index == pending_index {
                quote!(#field: #conversion)
            } else {
                quote!(#field: #local)
            }
        });
    Ok(quote! {
        impl<Access> crate::PendingResult<Access> for #name<crate::PendingMolecule<Access>> {
            type Finished = #name<#molecule_type>;
            fn resolve_pending(
                self,
                finalizer: &mut crate::ResultFinalizer<'_, Access>,
            ) -> Result<Self::Finished, crate::OperationError> {
                let #name { #(#names: #locals),* } = self;
                Ok(#name { #(#assignments),* })
            }
        }
    })
}
