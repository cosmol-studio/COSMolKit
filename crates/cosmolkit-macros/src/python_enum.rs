//! Enum extraction is language projection, not a chemistry conversion.
use proc_macro2::TokenStream;
use quote::quote;
use syn::{Fields, ItemEnum, Meta};

fn spelling(name: &str) -> String {
    let chars: Vec<_> = name.chars().collect();
    let mut result = String::new();
    for (index, ch) in chars.iter().copied().enumerate() {
        if index > 0 {
            let previous = chars[index - 1];
            let boundary = (ch.is_ascii_uppercase()
                && (previous.is_ascii_lowercase()
                    || previous.is_ascii_uppercase()
                        && chars.get(index + 1).is_some_and(char::is_ascii_lowercase)))
                || ch.is_ascii_digit() && previous.is_ascii_lowercase();
            if boundary {
                result.push('_');
            }
        }
        result.push(ch.to_ascii_lowercase());
    }
    result
}

pub(crate) fn expand(mut item: ItemEnum, existing_methods: bool) -> syn::Result<TokenStream> {
    for variant in &item.variants {
        if !matches!(variant.fields, Fields::Unit) {
            return Err(syn::Error::new_spanned(
                variant,
                "Python enum inputs require unit variants",
            ));
        }
    }
    let attribute = item
        .attrs
        .iter_mut()
        .find(|attr| attr.path().is_ident("pyclass"))
        .ok_or_else(|| syn::Error::new_spanned(&item.ident, "python_enum requires pyclass"))?;
    let arguments = attribute
        .parse_args_with(syn::punctuated::Punctuated::<Meta, syn::Token![,]>::parse_terminated)?;
    let arguments: Vec<_> = arguments
        .into_iter()
        .filter(|arg| {
            !arg.path().is_ident("from_py_object") && !arg.path().is_ident("skip_from_py_object")
        })
        .collect();
    *attribute = syn::parse_quote!(#[pyclass(#(#arguments,)* skip_from_py_object)]);
    let name = &item.ident;
    let variants: Vec<_> = item.variants.iter().map(|v| &v.ident).collect();
    let strings: Vec<_> = variants.iter().map(|v| spelling(&v.to_string())).collect();
    let expected = strings.join(", ");
    let methods = (!existing_methods).then(|| {
        quote! {
            #[pyo3::pymethods]
            impl #name {
                #[classattr]
                fn _enum_string_values() -> Vec<(&'static str, Self)> {
                    Self::enum_string_values()
                }
            }
        }
    });
    Ok(quote! {
        #item

        impl<'a, 'py> pyo3::FromPyObject<'a, 'py> for #name {
            type Error = pyo3::PyErr;
            fn extract(value: pyo3::Borrowed<'a, 'py, pyo3::PyAny>) -> Result<Self, pyo3::PyErr> {
                use pyo3::prelude::*;
                if value.is_instance_of::<Self>() {
                    return Ok(value.extract::<pyo3::PyRef<'_, Self>>()?.clone());
                }
                if value.is_instance_of::<pyo3::types::PyString>() {
                    let text = value.extract::<std::borrow::Cow<'_, str>>()?;
                    return match text.as_ref() {
                        #(#strings => Ok(Self::#variants),)*
                        _ => Err(pyo3::exceptions::PyValueError::new_err(format!(
                            "invalid {} string {:?}; expected one of: {}", stringify!(#name), text, #expected
                        ))),
                    };
                }
                Err(pyo3::exceptions::PyTypeError::new_err(concat!(
                    "expected ", stringify!(#name), " or str"
                )))
            }
        }

        impl #name {
            // Generated from the same variants as native extraction. Stub
            // generation and contract checks consume this, never a second list.
            fn enum_string_values() -> Vec<(&'static str, Self)> {
                vec![#((#strings, Self::#variants)),*]
            }
        }
        #methods
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_spelling_includes_dimension_suffixes_and_acronyms() {
        for (member, text) in [
            ("Preserve", "preserve"),
            ("Require3D", "require_3d"),
            ("TwoD", "two_d"),
            ("V2000", "v2000"),
            ("NonStrict", "non_strict"),
        ] {
            assert_eq!(spelling(member), text);
        }
    }

    #[test]
    fn replaces_clone_only_extraction_and_preserves_class_options() {
        let output = expand(
            syn::parse_quote! {
                #[pyclass(module = "cosmolkit", frozen, from_py_object)]
                #[derive(Clone)]
                enum Mode { First, Second }
            },
            false,
        )
        .unwrap()
        .to_string();
        assert!(output.contains("skip_from_py_object"));
        assert!(output.contains("frozen"));
        assert!(output.contains("PyValueError"));
        assert!(output.contains("PyTypeError"));
        assert!(output.contains("_enum_string_values"));
    }

    #[test]
    fn rejects_non_enum_vocabulary() {
        assert!(
            expand(
                syn::parse_quote! { #[pyclass] enum Mode { Value(u32) } },
                false
            )
            .is_err()
        );
    }
}
