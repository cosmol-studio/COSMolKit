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
        .filter(|operation| operation.fields.receiver_type.is_none())
        .map(expand_molecule_operation)
        .collect::<syn::Result<Vec<_>>>()?;
    let reconstruction_operations = registry
        .operations
        .iter()
        .filter(|operation| operation.fields.receiver_type.is_some())
        .map(expand_molecule_operation)
        .collect::<syn::Result<Vec<_>>>()?;

    Ok(quote! {
        impl crate::Molecule {
            #(#operations)*
        }
        #(#reconstruction_operations)*
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
    let visibility = &fields.method_visibility;
    let error_type = &fields.error_type;
    let params = &operation.params;
    let call_args = call_arguments(operation)?;
    let impl_fn = &fields.impl_fn;
    let spec = format_ident!("{}_SPEC", name.to_string().to_ascii_uppercase());
    let docs = fields.docs.as_ref().map(|text| quote!(#[doc = #text]));

    if let Some(receiver_type) = &fields.receiver_type {
        let result = fields
            .result_type
            .as_ref()
            .expect("validated reconstruction result");
        let assemble = fields
            .assemble_fn
            .as_ref()
            .expect("validated reconstruction assembler");
        return Ok(quote! {
            #(#cfg)*
            impl #receiver_type {
                #docs
                #visibility fn #method(&mut self, #(#params),*) -> Result<#result, #error_type> {
                    let mut parts = crate::MultiOutputOpParts::new_reconstruction(&#spec)?;
                    let metadata = #impl_fn(&mut parts, self, #(#call_args),*)?;
                    let molecules = parts.finish()?;
                    #assemble(molecules, metadata).map_err(<#error_type as ::core::convert::From<crate::ops::OperationError>>::from)
                }
            }
        });
    }

    let detached_result = fields
        .report_type
        .as_ref()
        .or(fields.inplace_result_type.as_ref());
    let primary = if let Some(result) = detached_result {
        let return_type = molecule_value_return_type(fields);
        let committed_result = if let Some(public_result) = fields.report_result_type.as_ref() {
            quote!(<#public_result as ::core::convert::From<(crate::Molecule, #result)>>::from((molecule, result)))
        } else if fields.report_type.is_some() {
            quote!((molecule, result))
        } else {
            quote!({
                let _ = result;
                molecule
            })
        };
        quote! {
            #(#cfg)*
            #docs
            #visibility fn #method(&self, #(#params),*) -> Result<#return_type, #error_type> {
                let mut parts = crate::OpParts::new(self, &#spec)?;
                let result: #result = #impl_fn(&mut parts, #(#call_args),*)?;
                let molecule = parts.finish()?;
                Ok(#committed_result)
            }
        }
    } else {
        match (
            fields.output,
            fields.result_type.as_ref(),
            fields.assemble_fn.as_ref(),
        ) {
            (MoleculeOutput::Single, None, None) => quote! {
                #(#cfg)*
                #docs
                #visibility fn #method(&self, #(#params),*) -> Result<crate::Molecule, #error_type> {

                    let mut parts = crate::OpParts::new(self, &#spec)?;
                    #impl_fn(&mut parts, #(#call_args),*)?;
                    parts.finish().map_err(<#error_type as ::core::convert::From<crate::ops::OperationError>>::from)
                }
            },
            (MoleculeOutput::Single, Some(result), None) => {
                let marker = crate::projection::access_marker(name)?;
                quote! {
                    #(#cfg)*
                    #docs
                    #visibility fn #method(&self, #(#params),*) -> Result<#result, #error_type> {

                        let mut parts = crate::OpParts::new(self, &#spec)?;
                        let pending: #result<crate::PendingMolecule<#marker>> = #impl_fn(&mut parts, #(#call_args),*)?;
                        parts.finish_result(pending).map_err(<#error_type as ::core::convert::From<crate::ops::OperationError>>::from)
                    }
                }
            }
            (MoleculeOutput::LazyMultiple, None, None) => quote! {
                #(#cfg)*
                #docs
                #visibility fn #method(&self, #(#params),*) -> Result<crate::StereoisomerIterator, #error_type> {
                    let mut parts = crate::MultiOutputOpParts::new(self, &#spec)?;
                    #impl_fn(&mut parts, #(#call_args),*)?;
                    parts.finish_lazy().map_err(<#error_type as ::core::convert::From<crate::ops::OperationError>>::from)
                }
            },
            (MoleculeOutput::Multiple, None, None) => quote! {
                #(#cfg)*
                #docs
                #visibility fn #method(&self, #(#params),*) -> Result<Vec<crate::Molecule>, #error_type> {

                    let mut parts = crate::MultiOutputOpParts::new(self, &#spec)?;
                    #impl_fn(&mut parts, #(#call_args),*)?;
                    parts.finish().map_err(<#error_type as ::core::convert::From<crate::ops::OperationError>>::from)
                }
            },
            (MoleculeOutput::Multiple, Some(result), Some(assemble)) => quote! {
                #(#cfg)*
                #docs
                #visibility fn #method(&self, #(#params),*) -> Result<#result, #error_type> {

                    let mut parts = crate::MultiOutputOpParts::new(self, &#spec)?;
                    let metadata = #impl_fn(&mut parts, #(#call_args),*)?;
                    let molecules = parts.finish()?;
                    #assemble(molecules, metadata).map_err(<#error_type as ::core::convert::From<crate::ops::OperationError>>::from)
                }
            },
            _ => {
                return Err(syn::Error::new_spanned(
                    name,
                    "validated molecule result/assembler relationship became inconsistent",
                ));
            }
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
        if let Some(result) = detached_result {
            let return_type = fields
                .inplace_result_type
                .as_ref()
                .or(fields.report_result_type.as_ref())
                .unwrap_or(result);
            let committed_return = if let Some(inplace_result) = fields.inplace_result_type.as_ref()
            {
                // D2: named value report and scalar in-place report share the
                // same body and checked finish. No snapshot is required for
                // this projection; From is checked by the Rust type system.
                quote! {
                    parts.finish_in_place()?;
                    Ok(<#inplace_result as ::core::convert::From<#result>>::from(result))
                }
            } else if let Some(public_result) = fields.report_result_type.as_ref() {
                // Borrowing the transaction ends at finish_in_place. Snapshot only
                // the finalized live value, whose blocks remain Arc-shared.
                quote! {
                    parts.finish_in_place()?;
                    Ok(<#public_result as ::core::convert::From<(crate::Molecule, #result)>>::from((self.clone(), result)))
                }
            } else {
                quote! {
                    parts.finish_in_place()?;
                    Ok(result)
                }
            };
            quote! {
                #(#cfg)*
                #inplace_docs
                #visibility fn #inplace_method(&mut self, #(#params),*) -> Result<#return_type, #error_type> {

                    let mut parts = crate::OpParts::new_in_place(self, &#spec)?;
                    let result = match #impl_fn(&mut parts, #(#call_args),*) {
                        Ok(result) => result,
                        Err(error) => {
                            parts.abort_in_place();
                            return Err(<#error_type as ::core::convert::From<crate::ops::OperationError>>::from(error));
                        }
                    };
                    #committed_return
                }
            }
        } else {
            quote! {
                #(#cfg)*
                #inplace_docs
                #visibility fn #inplace_method(&mut self, #(#params),*) -> Result<(), #error_type> {

                    let mut parts = crate::OpParts::new_in_place(self, &#spec)?;
                    if let Err(error) = #impl_fn(&mut parts, #(#call_args),*) {
                        parts.abort_in_place();
                        return Err(<#error_type as ::core::convert::From<crate::ops::OperationError>>::from(error));
                    }
                    parts.finish_in_place().map_err(<#error_type as ::core::convert::From<crate::ops::OperationError>>::from)
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
            #visibility fn #default_method(&self #forwarded_signature) -> Result<#return_type, #error_type> {
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
                .inplace_result_type
                .as_ref()
                .or(fields.report_result_type.as_ref())
                .or(detached_result)
                .map_or_else(|| quote!(()), |result| quote!(#result));
            quote! {
                #(#cfg)*
                #visibility fn #default_method(&mut self #forwarded_signature) -> Result<#return_type, #error_type> {
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

    let statuses = [
        Some(&fields.method),
        fields.default_method.as_ref(),
        fields.inplace_method.as_ref(),
        fields.default_inplace_method.as_ref(),
    ]
    .into_iter()
    .flatten()
    .map(|method| {
        let constant = format_ident!(
            "__FUNCTION_STATUS_{}",
            method.to_string().to_ascii_uppercase()
        );
        quote! { #(#cfg)* pub(crate) const #constant: crate::FunctionStatus = #spec.status; }
    });
    Ok(quote! {
        #(#statuses)*
        #primary
        #inplace
        #default
        #default_inplace
    })
}

fn molecule_value_return_type(
    fields: &crate::declaration::MoleculeFields,
) -> proc_macro2::TokenStream {
    if let Some(result) = fields.report_result_type.as_ref() {
        return quote!(#result);
    }
    if let Some(report) = fields.report_type.as_ref() {
        if let Some(result) = &fields.report_result_type {
            return quote!(#result);
        }
        return quote!((crate::Molecule, #report));
    }
    match (
        fields.output,
        fields.result_type.as_ref(),
        fields.assemble_fn.as_ref(),
    ) {
        (MoleculeOutput::Single, None, None) => quote!(crate::Molecule),
        (MoleculeOutput::Single, Some(result), None) => quote!(#result),
        (MoleculeOutput::Multiple, None, None) => quote!(Vec<crate::Molecule>),
        (MoleculeOutput::LazyMultiple, None, None) => quote!(crate::StereoisomerIterator),
        (MoleculeOutput::Multiple, Some(result), Some(_)) => quote!(#result),
        _ => unreachable!("molecule result/assembler shape was validated before wrapper expansion"),
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
