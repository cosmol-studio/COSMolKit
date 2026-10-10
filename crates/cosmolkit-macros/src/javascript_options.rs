//! Project constructor fields without a parallel configuration inventory.
use proc_macro2::TokenStream;
use quote::{format_ident, quote};
use syn::{FnArg, ImplItem, ItemImpl, LitStr, Meta, Pat};

pub(crate) struct Options {
    interface: LitStr,
    clone_snapshot: bool,
}

impl syn::parse::Parse for Options {
    fn parse(input: syn::parse::ParseStream) -> syn::Result<Self> {
        let interface = input.parse()?;
        let clone_snapshot = if input.is_empty() {
            false
        } else {
            input.parse::<syn::Token![,]>()?;
            let option: syn::Ident = input.parse()?;
            if option != "clone_snapshot" {
                return Err(syn::Error::new_spanned(option, "expected clone_snapshot"));
            }
            true
        };
        Ok(Self {
            interface,
            clone_snapshot,
        })
    }
}

#[cfg(test)]
fn expand(item: ItemImpl, interface: LitStr) -> syn::Result<TokenStream> {
    expand_with_options(
        item,
        Options {
            interface,
            clone_snapshot: false,
        },
    )
}

pub(crate) fn expand_with_options(
    mut item: ItemImpl,
    options: Options,
) -> syn::Result<TokenStream> {
    let Options {
        interface,
        clone_snapshot,
    } = options;
    let constructor = item
        .items
        .iter()
        .find_map(|member| match member {
            ImplItem::Fn(method)
                if method.attrs.iter().any(|attr| {
                    attr.path().is_ident("wasm_bindgen") && attr.parse_args_with(
                    syn::punctuated::Punctuated::<Meta, syn::Token![,]>::parse_terminated,
                ).is_ok_and(|options| options.iter().any(|option| matches!(option,
                    Meta::Path(path) if path.is_ident("constructor"))))
                }) =>
            {
                Some(method)
            }
            _ => None,
        })
        .or_else(|| {
            item.items.iter().find_map(|member| match member {
                ImplItem::Fn(method) if method.sig.ident == "new" => Some(method),
                _ => None,
            })
        })
        .ok_or_else(|| {
            syn::Error::new_spanned(&item, "options require an existing new constructor")
        })?;
    let constructor_name = constructor.sig.ident.clone();
    let mut fields = Vec::new();
    let mut inputs = Vec::new();
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
        inputs.push((binding.ident.clone(), name.clone(), ty));
        fields.push(name);
    }
    declaration.push_str("}\n");
    // Rebuild through the existing constructor before replacing the value.
    // This preserves every other configurable field and gives setters exactly
    // the constructor's validation and failure atomicity, including nested
    // configuration objects. Never generate setters for calculated outputs.
    for (field, name, javascript_type) in &inputs {
        let setter = format_ident!("set_{}", field);
        if let Some(ImplItem::Fn(method)) = item
            .items
            .iter_mut()
            .find(|member| matches!(member, ImplItem::Fn(method) if method.sig.ident == setter))
        {
            // A handwritten validating setter still projects the constructor's
            // declared type. JsValue is transport, never a public `any` escape.
            for argument in &mut method.sig.inputs {
                if let FnArg::Typed(argument) = argument
                    && matches!(&*argument.ty, syn::Type::Path(path) if path.path.segments.last().is_some_and(|segment| segment.ident == "JsValue"))
                    && !argument
                        .attrs
                        .iter()
                        .any(|attr| attr.path().is_ident("wasm_bindgen"))
                {
                    argument.attrs.push(
                        syn::parse_quote!(#[wasm_bindgen(unchecked_param_type = #javascript_type)]),
                    );
                }
            }
            continue;
        }
        let arguments = inputs.iter().map(|(other, _, _)| {
            if other == field {
                return Ok(quote! { value });
            }
            let getter = item.items.iter().find_map(|member| match member {
                ImplItem::Fn(method) if method.sig.ident == *other && method.sig.receiver().is_some() => Some(method),
                _ => None,
            }).ok_or_else(|| syn::Error::new_spanned(other, "writable configuration requires its constructor field's getter"))?;
            let fallible = matches!(&getter.sig.output,
                syn::ReturnType::Type(_, ty) if matches!(&**ty,
                    syn::Type::Path(path) if path.path.segments.last().is_some_and(|segment| segment.ident == "Result")));
            Ok(if fallible {
                quote! { wasm_bindgen::JsValue::from(self.#other()?) }
            } else {
                quote! { wasm_bindgen::JsValue::from(self.#other()) }
            })
        }).collect::<syn::Result<Vec<_>>>()?;
        item.items.push(syn::parse2(quote! {
            #[wasm_bindgen(setter, js_name = #name)]
            pub fn #setter(&mut self,
                #[wasm_bindgen(unchecked_param_type = #javascript_type)]
                value: wasm_bindgen::JsValue,
            ) -> Result<(), wasm_bindgen::JsValue> {
                let replacement = Self::#constructor_name(#(#arguments),*)?;
                *self = replacement;
                Ok(())
            }
        })?);
    }
    let ty = &item.self_ty;
    let constant = format_ident!("{}_OPTIONS_TYPESCRIPT", interface.value().to_uppercase());
    let visitor = format_ident!("visit_{}_configuration", interface.value());
    let visitor_name = visitor.to_string();
    let visitor_js = format!("export function {visitor_name}(value, visit) {{ visit(value); }}");
    let defaults = inputs
        .iter()
        .map(|_| quote! { wasm_bindgen::JsValue::UNDEFINED });
    let getters = inputs
        .iter()
        .map(|(field, _, _)| {
            let fallible = item.items.iter().any(|member| {
                matches!(member,
            ImplItem::Fn(method) if method.sig.ident == *field && matches!(&method.sig.output,
                syn::ReturnType::Type(_, ty) if matches!(&**ty, syn::Type::Path(path)
                    if path.path.segments.last().is_some_and(|segment| segment.ident == "Result"))))
            });
            if fallible {
                quote! { wasm_bindgen::JsValue::from(parameters.#field()?) }
            } else {
                quote! { wasm_bindgen::JsValue::from(parameters.#field()) }
            }
        })
        .collect::<Vec<_>>();
    let snapshot = if clone_snapshot {
        // Some configurations contain preset state beyond constructor fields.
        // Their explicit Clone implementation preserves that state as well.
        quote! { Ok::<Self, wasm_bindgen::JsValue>(parameters.clone()) }
    } else {
        quote! { Self::#constructor_name(#(#getters),*) }
    };
    Ok(quote! {
        #item
        #[wasm_bindgen::prelude::wasm_bindgen(inline_js = #visitor_js)]
        extern "C" {
            #[wasm_bindgen(catch, js_name = #visitor_name)]
            fn #visitor(value: &wasm_bindgen::JsValue, visit: &mut dyn FnMut(&#ty)) -> Result<(), wasm_bindgen::JsValue>;
        }
        impl #ty {
            pub(crate) fn from_js_configuration(value: &wasm_bindgen::JsValue) -> Result<Self, wasm_bindgen::JsValue> {
                if value.is_undefined() {
                    return Self::#constructor_name(#(#defaults),*);
                }
                let mut snapshot = None;
                if #visitor(value, &mut |parameters: &Self| {
                    snapshot = Some((|| #snapshot)());
                }).is_ok() {
                    return snapshot.ok_or_else(|| js_sys::TypeError::new("invalid configuration instance"))?;
                }
                Self::from_js_options(value)
            }
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
                Self::#constructor_name(#(js_sys::Reflect::get(value, &wasm_bindgen::JsValue::from_str(#fields))?),*)
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
    fn handwritten_setter_keeps_constructor_type_instead_of_any() {
        let output = expand(syn::parse_quote! {
            impl Params {
                fn new(#[wasm_bindgen(unchecked_optional_param_type="boolean")] flag: JsValue) -> Result<Self, JsValue> { unimplemented!() }
                fn flag(&self) -> bool { false }
                #[wasm_bindgen(setter)]
                fn set_flag(&mut self, value: JsValue) -> Result<(), JsValue> { unimplemented!() }
            }
        }, syn::parse_quote!("ParamsOptions")).unwrap().to_string();
        assert!(output.contains("unchecked_param_type = \"boolean\""));
        assert_eq!(output.matches("fn set_flag").count(), 1);
    }
    #[test]
    fn actual_constructor_and_snapshot_preserve_factory_configuration() {
        let output = expand_with_options(syn::parse_quote! {
            impl Params {
                #[wasm_bindgen(constructor)]
                fn create(#[wasm_bindgen(unchecked_optional_param_type="number")] count: JsValue) -> Result<Self, JsValue> { unimplemented!() }
                fn count(&self) -> u32 { 1 }
            }
        }, Options { interface: syn::parse_quote!("ParamsOptions"), clone_snapshot: true }).unwrap().to_string();
        assert!(output.contains("Self :: create"));
        assert!(output.contains("parameters . clone ()"));
        assert!(!output.contains("Self :: new"));
    }
    #[test]
    fn options_reuse_constructor_validation_and_declared_field_types() {
        let output = expand(syn::parse_quote! {
            impl Params {
                fn new(#[wasm_bindgen(unchecked_optional_param_type="boolean")] use_chirality: JsValue) -> Result<Self, JsValue> { unimplemented!() }
                fn use_chirality(&self) -> bool { false }
            }
        }, syn::parse_quote!("ParamsOptions")).unwrap().to_string();
        assert!(output.contains("useChirality?: boolean"));
        assert!(output.contains("Self :: new"));
        assert!(output.contains("unknown configuration option"));
        assert!(output.contains("set_use_chirality"));
    }

    #[test]
    fn setter_preserves_other_fields_and_validates_before_assignment() {
        let output = expand(
            syn::parse_quote! {
                impl Params {
                    fn new(
                        #[wasm_bindgen(unchecked_optional_param_type="boolean")] flag: JsValue,
                        #[wasm_bindgen(unchecked_optional_param_type="number")] count: JsValue,
                    ) -> Result<Self, JsValue> { unimplemented!() }
                    fn flag(&self) -> bool { false }
                    fn count(&self) -> u32 { 1 }
                }
            },
            syn::parse_quote!("ParamsOptions"),
        )
        .unwrap()
        .to_string();
        assert!(output.contains("JsValue :: from (self . count ())"));
        assert!(output.contains("let replacement = Self :: new"));
        assert!(output.contains("* self = replacement"));
    }
}
