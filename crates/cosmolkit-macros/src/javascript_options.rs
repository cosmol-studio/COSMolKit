//! Project constructor fields without a parallel configuration inventory.
use proc_macro2::TokenStream;
use quote::{format_ident, quote};
use syn::{FnArg, ImplItem, ItemImpl, LitStr, Meta, Pat};

pub(crate) fn expand(item: ItemImpl, interface: LitStr) -> syn::Result<TokenStream> {
    let constructor = item
        .items
        .iter()
        .find_map(|member| match member {
            ImplItem::Fn(method) if method.sig.ident == "new" => Some(method),
            _ => None,
        })
        .ok_or_else(|| {
            syn::Error::new_spanned(&item, "options require an existing new constructor")
        })?;
    let mut fields = Vec::new();
    let mut declaration = format!("export interface {} {{\n", interface.value());
    for input in &constructor.sig.inputs {
        let FnArg::Typed(argument) = input else {
            return Err(syn::Error::new_spanned(
                input,
                "constructor inputs must be named",
            ));
        };
        let Pat::Ident(binding) = &*argument.pat else {
            return Err(syn::Error::new_spanned(
                input,
                "constructor inputs must be named",
            ));
        };
        let mut name = String::new();
        let mut upper = false;
        for ch in binding.ident.to_string().chars() {
            if ch == '_' {
                upper = true;
            } else if upper {
                name.push(ch.to_ascii_uppercase());
                upper = false;
            } else {
                name.push(ch);
            }
        }
        let mut ty = None;
        for attr in &argument.attrs {
            if !attr.path().is_ident("wasm_bindgen") {
                continue;
            }
            let options = attr.parse_args_with(
                syn::punctuated::Punctuated::<Meta, syn::Token![,]>::parse_terminated,
            )?;
            for option in options {
                if let Meta::NameValue(value) = option
                    && value.path.is_ident("unchecked_optional_param_type")
                    && let syn::Expr::Lit(literal) = value.value
                    && let syn::Lit::Str(text) = literal.lit
                {
                    ty = Some(text.value());
                }
            }
        }
        let ty = ty.ok_or_else(|| {
            syn::Error::new_spanned(
                input,
                "options require the constructor's explicit optional JavaScript type",
            )
        })?;
        declaration.push_str(&format!("  {name}?: {ty};\n"));
        fields.push(name);
    }
    declaration.push_str("}\n");
    let ty = &item.self_ty;
    let constant = format_ident!("{}_OPTIONS_TYPESCRIPT", interface.value().to_uppercase());
    Ok(quote! {
        #item
        impl #ty {
            pub(crate) fn from_js_options(value: &wasm_bindgen::JsValue) -> Result<Self, wasm_bindgen::JsValue> {
                if !value.is_object() || value.is_null() || js_sys::Array::is_array(value) {
                    return Err(js_sys::TypeError::new("expected a configuration options object").into());
                }
                let object = js_sys::Object::from(value.clone());
                let prototype = js_sys::Object::get_prototype_of(&object);
                if !prototype.is_null() && prototype != js_sys::Object::get_prototype_of(&js_sys::Object::new()) {
                    return Err(js_sys::TypeError::new("expected a plain configuration options object").into());
                }
                for key in js_sys::Object::keys(&object).iter() {
                    let name = key.as_string().ok_or_else(|| js_sys::TypeError::new("invalid option name"))?;
                    if ![#(#fields),*].contains(&name.as_str()) {
                        return Err(js_sys::TypeError::new(&format!("unknown configuration option: {name}")).into());
                    }
                }
                Self::new(#(js_sys::Reflect::get(value, &wasm_bindgen::JsValue::from_str(#fields))?),*)
            }
        }
        #[wasm_bindgen::prelude::wasm_bindgen(typescript_custom_section)]
        const #constant: &'static str = #declaration;
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn options_reuse_constructor_validation_and_declared_field_types() {
        let output = expand(syn::parse_quote! {
            impl Params {
                fn new(#[wasm_bindgen(unchecked_optional_param_type="boolean")] use_chirality: JsValue) -> Result<Self, JsValue> { unimplemented!() }
            }
        }, syn::parse_quote!("ParamsOptions")).unwrap().to_string();
        assert!(output.contains("useChirality?: boolean"));
        assert!(output.contains("Self :: new"));
        assert!(output.contains("unknown configuration option"));
    }
}
