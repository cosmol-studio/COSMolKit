//! Private declaration-driven operation wrapper generation.
//!
//! Generated tokens name runtime-owned values at the call site. This module
//! owns no molecule state, operation runtime, or domain behavior.

use quote::{format_ident, quote};
use syn::{Ident, Pat};

use crate::declaration::{
    BioOperation, BioRegistry, MoleculeOperation, MoleculeOutput, MoleculeRegistry,
};

pub(crate) fn expand_molecule_wrappers(
    registry: &MoleculeRegistry,
) -> syn::Result<proc_macro2::TokenStream> {
    let operations = registry
        .operations
        .iter()
        .map(expand_molecule_operation)
        .collect::<syn::Result<Vec<_>>>()?;

    Ok(quote! {
        impl crate::Molecule {
            #(#operations)*
        }
    })
}

pub(crate) fn expand_bio_wrappers(registry: &BioRegistry) -> syn::Result<proc_macro2::TokenStream> {
    let operations = registry
        .operations
        .iter()
        .map(expand_bio_operation)
        .collect::<syn::Result<Vec<_>>>()?;

    Ok(quote! {
        impl crate::BioStructure {
            #(#operations)*
        }
    })
}

fn expand_molecule_operation(
    operation: &MoleculeOperation,
) -> syn::Result<proc_macro2::TokenStream> {
    let cfg = &operation.cfg_attrs;
    let name = &operation.name;
    let fields = &operation.fields;
    let method = &fields.method;
    let params = &operation.params;
    let call_args = call_arguments(operation)?;
    let impl_fn = &fields.impl_fn;
    let feature = &fields.feature;
    let spec = format_ident!("{}_SPEC", name.to_string().to_ascii_uppercase());
    let docs = fields.docs.as_ref().map(|text| quote!(#[doc = #text]));

    let support_check = quote! {
        if matches!(
            #spec.support,
            crate::SupportStatus::Unsupported { .. }
        ) {
            return Err(crate::ops::OperationError::UnsupportedFeature {
                operation: &#spec,
                source: crate::UnsupportedFeatureError::from_spec(&#feature),
            });
        }
    };

    let primary = match (
        fields.output,
        fields.result_type.as_ref(),
        fields.assemble_fn.as_ref(),
    ) {
        (MoleculeOutput::Single, None, None) => quote! {
            #(#cfg)*
            #docs
            pub fn #method(&self, #(#params),*) -> Result<crate::Molecule, crate::ops::OperationError> {
                #support_check
                let mut parts = crate::OpParts::new(self, &#spec)?;
                #impl_fn(&mut parts, #(#call_args),*)?;
                parts.finish()
            }
        },
        (MoleculeOutput::Single, Some(result), None) => quote! {
            #(#cfg)*
            #docs
            pub fn #method(&self, #(#params),*) -> Result<(crate::Molecule, #result), crate::ops::OperationError> {
                #support_check
                let mut parts = crate::OpParts::new(self, &#spec)?;
                let result = #impl_fn(&mut parts, #(#call_args),*)?;
                let molecule = parts.finish()?;
                Ok((molecule, result))
            }
        },
        (MoleculeOutput::Multiple, None, None) => quote! {
            #(#cfg)*
            #docs
            pub fn #method(&self, #(#params),*) -> Result<Vec<crate::Molecule>, crate::ops::OperationError> {
                #support_check
                let mut parts = crate::MultiOutputOpParts::new(self, &#spec)?;
                #impl_fn(&mut parts, #(#call_args),*)?;
                parts.finish()
            }
        },
        (MoleculeOutput::Multiple, Some(result), Some(assemble)) => quote! {
            #(#cfg)*
            #docs
            pub fn #method(&self, #(#params),*) -> Result<#result, crate::ops::OperationError> {
                #support_check
                let mut parts = crate::MultiOutputOpParts::new(self, &#spec)?;
                let metadata = #impl_fn(&mut parts, #(#call_args),*)?;
                let molecules = parts.finish()?;
                #assemble(molecules, metadata)
            }
        },
        _ => {
            return Err(syn::Error::new_spanned(
                name,
                "validated molecule result/assembler relationship became inconsistent",
            ));
        }
    };

    let inplace = if fields.inplace {
        let inplace_method = fields.inplace_method.as_ref().ok_or_else(|| {
            syn::Error::new_spanned(name, "validated in-place method name is missing")
        })?;
        let inplace_docs = fields
            .inplace_docs
            .as_ref()
            .map(|text| quote!(#[doc = #text]));
        if let Some(result) = fields.result_type.as_ref() {
            quote! {
                #(#cfg)*
                #inplace_docs
                pub fn #inplace_method(&mut self, #(#params),*) -> Result<#result, crate::ops::OperationError> {
                    #support_check
                    let mut parts = crate::OpParts::new_in_place(self, &#spec)?;
                    let result = match #impl_fn(&mut parts, #(#call_args),*) {
                        Ok(result) => result,
                        Err(error) => {
                            parts.abort_in_place();
                            return Err(error);
                        }
                    };
                    parts.finish_in_place()?;
                    Ok(result)
                }
            }
        } else {
            quote! {
                #(#cfg)*
                #inplace_docs
                pub fn #inplace_method(&mut self, #(#params),*) -> Result<(), crate::ops::OperationError> {
                    #support_check
                    let mut parts = crate::OpParts::new_in_place(self, &#spec)?;
                    if let Err(error) = #impl_fn(&mut parts, #(#call_args),*) {
                        parts.abort_in_place();
                        return Err(error);
                    }
                    parts.finish_in_place()
                }
            }
        }
    } else {
        quote!()
    };

    let default = if let Some(default_method) = fields.default_method.as_ref() {
        let default_args = &fields.default_args;
        let forwarded_count = params.len() - default_args.len();
        let forwarded_params = &params[..forwarded_count];
        let forwarded_args = &call_args[..forwarded_count];
        let forwarded_signature = if forwarded_params.is_empty() {
            quote!()
        } else {
            quote!(, #(#forwarded_params),*)
        };
        let return_type = molecule_value_return_type(fields);
        quote! {
            #(#cfg)*
            pub fn #default_method(&self #forwarded_signature) -> Result<#return_type, crate::ops::OperationError> {
                self.#method(#(#forwarded_args,)* #(#default_args),*)
            }
        }
    } else {
        quote!()
    };

    let default_inplace = match (
        fields.default_inplace_method.as_ref(),
        fields.inplace_method.as_ref(),
    ) {
        (Some(default_method), Some(inplace_method)) => {
            let default_args = &fields.default_args;
            let forwarded_count = params.len() - default_args.len();
            let forwarded_params = &params[..forwarded_count];
            let forwarded_args = &call_args[..forwarded_count];
            let forwarded_signature = if forwarded_params.is_empty() {
                quote!()
            } else {
                quote!(, #(#forwarded_params),*)
            };
            let return_type = fields
                .result_type
                .as_ref()
                .map_or_else(|| quote!(()), |result| quote!(#result));
            quote! {
                #(#cfg)*
                pub fn #default_method(&mut self #forwarded_signature) -> Result<#return_type, crate::ops::OperationError> {
                    self.#inplace_method(#(#forwarded_args,)* #(#default_args),*)
                }
            }
        }
        (None, _) => quote!(),
        (Some(_), None) => {
            return Err(syn::Error::new_spanned(
                name,
                "validated default in-place method lost its primary in-place method",
            ));
        }
    };

    Ok(quote! {
        #primary
        #default
        #inplace
        #default_inplace
    })
}

fn molecule_value_return_type(
    fields: &crate::declaration::MoleculeFields,
) -> proc_macro2::TokenStream {
    match (fields.output, fields.result_type.as_ref()) {
        (MoleculeOutput::Single, None) => quote!(crate::Molecule),
        (MoleculeOutput::Single, Some(result)) => quote!((crate::Molecule, #result)),
        (MoleculeOutput::Multiple, None) => quote!(Vec<crate::Molecule>),
        (MoleculeOutput::Multiple, Some(result)) => quote!(#result),
    }
}

fn expand_bio_operation(operation: &BioOperation) -> syn::Result<proc_macro2::TokenStream> {
    let cfg = &operation.cfg_attrs;
    let name = &operation.name;
    let method = &operation.fields.method;
    let params = &operation.params;
    let call_args = call_arguments_from_params(name, params)?;
    let impl_fn = &operation.fields.impl_fn;
    let feature = &operation.fields.feature;
    let spec = format_ident!("BIO_{}_SPEC", name.to_string().to_ascii_uppercase());

    Ok(quote! {
        #(#cfg)*
        pub fn #method(&self, #(#params),*) -> Result<crate::BioStructure, crate::bio_ops::BioOperationError> {
            if let crate::SupportStatus::Unsupported { reason } = #spec.support {
                return Err(crate::bio_ops::BioOperationError::Unsupported {
                    operation: &#spec,
                    reason,
                });
            }
            let _feature = &#feature;
            let mut parts = crate::BioOpParts::new(self, &#spec);
            #impl_fn(&mut parts, #(#call_args),*)?;
            parts.finish()
        }
    })
}

fn call_arguments(operation: &MoleculeOperation) -> syn::Result<Vec<&Ident>> {
    call_arguments_from_params(&operation.name, &operation.params)
}

fn call_arguments_from_params<'a>(
    operation: &Ident,
    params: &'a [syn::PatType],
) -> syn::Result<Vec<&'a Ident>> {
    params
        .iter()
        .map(|parameter| match parameter.pat.as_ref() {
            Pat::Ident(binding) if binding.subpat.is_none() => Ok(&binding.ident),
            pattern => Err(syn::Error::new_spanned(
                pattern,
                format!("operation `{operation}` has a non-forwardable parameter pattern"),
            )),
        })
        .collect()
}
