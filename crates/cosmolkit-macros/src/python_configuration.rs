//! Language configuration glue only: reuse the existing constructor's validation.
use proc_macro2::TokenStream;
use quote::{format_ident, quote};
use syn::{FnArg, ImplItem, ItemImpl, Pat, Type};

fn setter_type(ty: &Type) -> Type {
    match ty {
        Type::Reference(reference) => match &*reference.elem {
            Type::Path(path)
                if path.path.segments.len() == 1
                    && !matches!(
                        path.path.segments[0].ident.to_string().as_str(),
                        "str" | "Bound"
                    ) =>
            {
                let target = &reference.elem;
                syn::parse_quote!(pyo3::Bound<'_, #target>)
            }
            _ => ty.clone(),
        },
        Type::Path(path) => {
            let mut path = path.clone();
            if let Some(segment) = path.path.segments.last_mut()
                && segment.ident == "Option"
                && let syn::PathArguments::AngleBracketed(arguments) = &mut segment.arguments
            {
                for argument in &mut arguments.args {
                    if let syn::GenericArgument::Type(ty) = argument {
                        *ty = setter_type(ty);
                    }
                }
            }
            Type::Path(path)
        }
        _ => ty.clone(),
    }
}

pub(crate) fn expand(mut implementation: ItemImpl, fieldwise: bool) -> syn::Result<TokenStream> {
    let constructor = implementation
        .items
        .iter()
        .find_map(|item| match item {
            ImplItem::Fn(method) if method.attrs.iter().any(|a| a.path().is_ident("new")) => {
                Some(method)
            }
            _ => None,
        })
        .ok_or_else(|| {
            syn::Error::new_spanned(
                &implementation,
                "configuration requires a #[new] constructor",
            )
        })?;
    let fields = constructor
        .sig
        .inputs
        .iter()
        .filter_map(|argument| match argument {
            FnArg::Typed(argument) => match (&*argument.pat, &*argument.ty) {
                (_, Type::Path(ty))
                    if ty
                        .path
                        .segments
                        .last()
                        .is_some_and(|part| part.ident == "Python") =>
                {
                    None
                }
                (Pat::Ident(name), _) => Some((
                    name.ident.clone(),
                    argument.ty.clone(),
                    argument.attrs.clone(),
                )),
                _ => None,
            },
            _ => None,
        })
        .collect::<Vec<_>>();
    let names = fields
        .iter()
        .map(|(name, _, _)| name.to_string())
        .collect::<Vec<_>>();
    for (field, ty, attributes) in fields {
        let ty = setter_type(&ty);
        let setter = format_ident!("set_{field}");
        if implementation
            .items
            .iter()
            .any(|item| matches!(item, ImplItem::Fn(method) if method.sig.ident == setter))
        {
            continue;
        }
        let name = field.to_string();
        let commit = if fieldwise {
            quote!(std::mem::swap(&mut current.inner.#field, &mut replacement.inner.#field);)
        } else {
            quote!(std::mem::swap(&mut *current, &mut *replacement);)
        };
        implementation.items.push(syn::parse2(quote! {
            #[setter(#field)]
            fn #setter(slf: &pyo3::Bound<'_, Self>, #(#attributes)* value: #ty) -> pyo3::PyResult<()> {
                use pyo3::prelude::*;
                let kwargs = pyo3::types::PyDict::new(slf.py());
                for key in [#(#names),*] {
                    if key != #name {
                        kwargs.set_item(key, slf.getattr(key)?)?;
                    }
                }
                kwargs.set_item(#name, value)?;
                // No mutation until both conversion and constructor validation succeed.
                let candidate = slf.get_type().call((), Some(&kwargs))?;
                let candidate = candidate.cast::<Self>()?;
                let mut replacement = candidate.try_borrow_mut()?;
                let mut current = slf.try_borrow_mut()?;
                #commit
                Ok(())
            }
        })?);
    }
    refresh(implementation)
}

pub(crate) fn refresh(mut implementation: ItemImpl) -> syn::Result<TokenStream> {
    implementation.attrs.insert(
        0,
        syn::parse_quote! {
            #[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
        },
    );
    // The Python projection keeps a native child object for a live nested
    // configuration view. Refresh only that detached value after validated
    // parent assignment; never expose a mutable Rust/runtime reference.
    implementation.items.push(syn::parse_quote! {
        #[gen_stub(skip)]
        fn _configuration_replace(
            slf: &pyo3::Bound<'_, Self>,
            value: &pyo3::Bound<'_, Self>,
        ) -> pyo3::PyResult<()> {
            use pyo3::prelude::*;
            if slf.is(value) {
                return Ok(());
            }
            std::mem::swap(&mut *slf.try_borrow_mut()?, &mut *value.try_borrow_mut()?);
            Ok(())
        }
    });
    Ok(quote!(#implementation))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constructor_fields_generate_typed_setters_without_replacing_existing_ones() {
        let implementation: ItemImpl = syn::parse_quote! {
            impl Params {
                #[new]
                fn new(py: Python<'_>, count: usize, flag: bool) -> Self { unimplemented!() }
                #[setter]
                fn set_flag(&mut self, value: bool) {}
            }
        };
        let output = expand(implementation, false).unwrap().to_string();
        assert_eq!(output.matches("fn set_count").count(), 1);
        assert_eq!(output.matches("fn set_flag").count(), 1);
        assert!(!output.contains("fn set_py"));
        assert!(output.contains("value : usize"));
        assert_eq!(output.matches("fn _configuration_replace").count(), 1);
        assert!(output.find("get_type").unwrap() < output.find("current =").unwrap());
    }

    #[test]
    fn missing_constructor_is_rejected() {
        assert!(expand(syn::parse_quote!(impl Params {}), false).is_err());
    }

    #[test]
    fn fieldwise_commit_keeps_state_outside_the_constructor() {
        let output = expand(
            syn::parse_quote! {
                impl Params {
                    #[new]
                    fn new(count: usize) -> Self { unimplemented!() }
                }
            },
            true,
        )
        .unwrap()
        .to_string();
        assert!(output.contains("current . inner . count"));
        assert!(!output.contains("& mut * current"));
    }
}
