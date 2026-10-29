//! Declaration generation only; runtime adaptation stays in the binding crate.
mod type_union;
use heck::ToLowerCamelCase;
use proc_macro::TokenStream;
use quote::{format_ident, quote};
use syn::{
    Field, Ident, Path, Token, Type, braced,
    parse::{Parse, ParseStream},
    punctuated::Punctuated,
};

/// A JS union whose members follow the same cfg gates as native error branches.
#[proc_macro]
pub fn type_union(input: TokenStream) -> TokenStream {
    type_union::expand(input)
}

struct Declaration {
    docs: Vec<syn::Attribute>,
    name: Ident,
    native: Path,
    fields: Punctuated<Field, Token![,]>,
    factory: bool,
}
impl Parse for Declaration {
    fn parse(input: ParseStream) -> syn::Result<Self> {
        let docs = input.call(syn::Attribute::parse_outer)?;
        let name = input.parse()?;
        input.parse::<Token![,]>()?;
        let native = input.parse()?;
        let content;
        braced!(content in input);
        let fields = content.parse_terminated(Field::parse_named, Token![,])?;
        let factory = if input.peek(Ident) {
            let marker: Ident = input.parse()?;
            if marker != "factory" {
                return Err(syn::Error::new_spanned(marker, "expected factory"));
            }
            true
        } else {
            false
        };
        Ok(Self {
            docs,
            name,
            native,
            fields,
            factory,
        })
    }
}

/// Input fields come from the canonical Configuration declaration.
#[proc_macro]
pub fn parameters(input: TokenStream) -> TokenStream {
    let Declaration {
        docs,
        name,
        native,
        fields,
        factory,
    } = syn::parse_macro_input!(input as Declaration);
    let factory_marker = if factory { quote!(factory) } else { quote!() };
    for field in &fields {
        if field.attrs.iter().any(|a| {
            a.path().is_ident("configuration")
                && a.parse_args::<Ident>().is_ok_and(|v| v == "flatten")
        }) {
            let mut field_type = &field.ty;
            while let Type::Group(group) = field_type {
                field_type = &group.elem;
            }
            let Type::Path(ty) = field_type else {
                return syn::Error::new_spanned(
                    &field.ty,
                    "flatten requires native configuration metadata",
                )
                .into_compile_error()
                .into();
            };
            let mut exporter = ty.path.clone();
            let last = exporter.segments.last_mut().unwrap();
            let original = last.ident.to_string();
            let Some(index) = original.strip_prefix("Field") else {
                return syn::Error::new_spanned(
                    &field.ty,
                    "flatten requires a canonical field metadata alias",
                )
                .into_compile_error()
                .into();
            };
            last.ident = format_ident!("Flatten{index}");
            let field_name = field.ident.as_ref().unwrap();
            let child = &ty.path;
            return quote! {
                #exporter!(crate::project_parameters,#child,{}, {
                    @flatten #(#docs)* #name, #native {#fields} #factory_marker, #field_name
                });
            }
            .into();
        }
    }
    let mut required = Vec::new();
    for field in &fields {
        let mut is_required = false;
        for attribute in &field.attrs {
            if attribute.path().is_ident("configuration") {
                let flag = match attribute.parse_args::<Ident>() {
                    Ok(flag) => flag,
                    Err(error) => return error.into_compile_error().into(),
                };
                if flag == "required" {
                    is_required = true;
                } else if flag != "none_default" && flag != "adapter" {
                    return syn::Error::new_spanned(
                        flag,
                        "configuration policy needs an explicit projection adapter",
                    )
                    .into_compile_error()
                    .into();
                }
            }
        }
        required.push(is_required);
    }
    let options = format_ident!("{}Options", name);
    let input_name = format_ident!("{}Input", name);
    let value_name = format_ident!("{}Value", name);
    let name_string = name.to_string();
    let input_type = format!("{name} | {options}");
    let names: Vec<_> = fields.iter().map(|f| f.ident.as_ref().unwrap()).collect();
    let field_docs: Vec<Vec<_>> = fields
        .iter()
        .map(|f| {
            f.attrs
                .iter()
                .filter(|a| a.path().is_ident("doc"))
                .collect()
        })
        .collect();
    let types: Vec<_> = fields.iter().map(|f| &f.ty).collect();
    let native_paths: Vec<_> = fields
        .iter()
        .map(|f| {
            match f
                .attrs
                .iter()
                .find(|a| a.path().is_ident("projection_path"))
            {
                Some(a) => a.meta.require_list().unwrap().tokens.clone(),
                None => {
                    let ident = f.ident.as_ref().unwrap();
                    quote!(#ident)
                }
            }
        })
        .collect();
    let js_names: Vec<_> = names
        .iter()
        .map(|f| f.to_string().to_lower_camel_case())
        .collect();
    let setters: Vec<_> = names.iter().map(|f| format_ident!("set_{}", f)).collect();
    // Explicit binding declarations select native accessor methods when the
    // canonical configuration deliberately keeps storage private.
    let accessor_reads: Vec<_> = fields
        .iter()
        .map(|field| {
            field
                .attrs
                .iter()
                .find(|a| a.path().is_ident("projection_read"))
                .map(|a| a.meta.require_list().unwrap().tokens.clone())
        })
        .collect();
    let assignments: Vec<_> = fields
        .iter()
        .zip(&native_paths)
        .map(|(field, path)| {
            match field
                .attrs
                .iter()
                .find(|a| a.path().is_ident("projection_write"))
            {
                Some(a) => {
                    let setter = &a.meta.require_list().unwrap().tokens;
                    quote!(parent.#setter(value))
                }
                None => quote!(parent.#path=value),
            }
        })
        .collect();
    let getters: Vec<_> = types.iter().zip(&native_paths).zip(&accessor_reads).map(|((ty,path),read)| {
        match read {
            Some(read) => quote!(<#ty as crate::fields::Field>::get(crate::storage::NativeField::owned(self.inner.read(|parent|parent.#read)?))),
            None => quote!(<#ty as crate::fields::Field>::get(self.inner.child(|value|Some(&value.#path),|value|Some(&mut value.#path)))),
        }
    }).collect();
    // Native metadata uses type aliases. Resolve the public type through the
    // generated getter, rather than asking Tsify to interpret an opaque alias.
    let public_types: Vec<_> = js_names
        .iter()
        .map(|field| format!("{name}[\"{field}\"]"))
        .collect();
    let option_annotations: Vec<_> = public_types
        .iter()
        .zip(&required)
        .map(|(ty, required)| {
            if *required {
                quote!(#[tsify(type = #ty)])
            } else {
                quote!(#[tsify(optional, type = #ty)])
            }
        })
        .collect();
    let has_required = required.iter().any(|value| *value);
    let constructor_input = if has_required {
        quote!(tsify::Ts<#options>)
    } else {
        quote!(Option<tsify::Ts<#options>>)
    };
    let normalize_options = if has_required {
        quote!(let options = Some(options);)
    } else {
        quote!()
    };
    let options_constructor = if has_required {
        quote!(Self::new(tsify::Ts::new_unchecked(raw.clone())))
    } else {
        quote!(Self::new(Some(tsify::Ts::new_unchecked(raw.clone()))))
    };
    let omitted_input = if has_required {
        quote!(Err(js_sys::TypeError::new(
            "required configuration fields are missing"
        )
        .into()))
    } else {
        quote!(visit(&<#native>::default()))
    };
    let finish_native = if factory {
        quote!(let native = <#native>::new(#(native.#native_paths),*);)
    } else {
        quote!()
    };
    let native_options = if !factory && !has_required {
        quote! {
            pub(crate) fn options_on(mut native: #native, raw: Option<&wasm_bindgen::JsValue>) -> Result<#native, wasm_bindgen::JsValue> {
                if let Some(raw) = raw {
                    crate::input::check_options(raw, &[#(#js_names),*])?;
                    #(
                        let value = js_sys::Reflect::get(raw, &wasm_bindgen::JsValue::from_str(#js_names))?;
                        if !value.is_undefined() {
                            let value = <#types as crate::fields::Field>::decode(&value)?;
                            let parent = &mut native;
                            #assignments;
                        }
                    )*
                }
                Ok(native)
            }
        }
    } else {
        quote!()
    };
    let construction = if !factory && !has_required {
        quote! {
            let raw = options.map(|value|value.js_value());
            let native = Self::options_on(<#native>::default(), raw.as_ref())?;
        }
    } else if !required.is_empty() && required.iter().all(|value| *value) {
        let native_value = if factory {
            quote!(<#native>::new(#(#names),*))
        } else {
            quote!(#native { #(#names),* })
        };
        quote! {
            let raw=options.js_value();
            crate::input::check_options(&raw, &[#(#js_names),*])?;
            #(
                let value=js_sys::Reflect::get(&raw,&wasm_bindgen::JsValue::from_str(#js_names))?;
                if value.is_undefined(){return Err(js_sys::TypeError::new(concat!("missing required option: ",#js_names)).into());}
                let #names=<#types as crate::fields::Field>::decode(&value)?;
            )*
            let native=#native_value;
        }
    } else {
        quote! {
            #normalize_options
            let mut native = <#native>::default();
            if let Some(options) = options {
                let raw = options.js_value();
                crate::input::check_options(&raw, &[#(#js_names),*])?;
                #(
                    let value = js_sys::Reflect::get(&raw, &wasm_bindgen::JsValue::from_str(#js_names))?;
                    if !value.is_undefined() {
                        let value = <#types as crate::fields::Field>::decode(&value)?;
                        let parent = &mut native;
                        #assignments;
                    } else if #required {
                        return Err(js_sys::TypeError::new(concat!("missing required option: ", #js_names)).into());
                    }
                )*
            }
            #finish_native
        }
    };
    quote! {
        #[derive(tsify::Tsify)]
        pub struct #options {
            #(
                #(#field_docs)*
                #[serde(rename = #js_names)]
                #option_annotations
                pub #names: #types,
            )*
        }
        #[wasm_bindgen::prelude::wasm_bindgen]
        #(#docs)*
        pub struct #name { inner: crate::storage::NativeField<#native> }
        #[wasm_bindgen::prelude::wasm_bindgen]
        extern "C" {
            #[wasm_bindgen(typescript_type = #input_type)]
            pub type #input_name;
            #[wasm_bindgen(typescript_type = #name_string)]
            pub type #value_name;
        }
        const _: () = {
            #[wasm_bindgen::prelude::wasm_bindgen]
            extern "C" {
                #[wasm_bindgen(catch, js_namespace = Reflect, js_name = apply)]
                fn borrow_configuration(
                    visit: &mut dyn FnMut(&#name),
                    receiver: &wasm_bindgen::JsValue,
                    arguments: &js_sys::Array,
                ) -> Result<wasm_bindgen::JsValue, wasm_bindgen::JsValue>;
            }
            impl #name {
                #native_options
                pub(crate) fn with_input<R>(value: Option<&#input_name>, visit: impl FnOnce(&#native) -> Result<R, wasm_bindgen::JsValue>) -> Result<R, wasm_bindgen::JsValue> {
                    let Some(value) = value else { return #omitted_input };
                    let raw: &wasm_bindgen::JsValue = value.as_ref();
                    if crate::input::is_plain_options(raw) {
                        return #options_constructor?.inner.read(visit)?;
                    }
                    // Official wasm-bindgen callback ABI borrows the class. It
                    // does not run the owning TryFromJsValue conversion.
                    Self::with_instance(raw, visit)
                }
                fn from_instance(raw: &wasm_bindgen::JsValue) -> Result<#native, wasm_bindgen::JsValue> {
                    Self::with_instance(raw, |native| Ok(native.clone()))
                }
                fn with_instance<R>(raw: &wasm_bindgen::JsValue, visit: impl FnOnce(&#native) -> Result<R, wasm_bindgen::JsValue>) -> Result<R, wasm_bindgen::JsValue> {
                    crate::input::with_instance(raw, #name_string,
                        |visit, receiver, args| borrow_configuration(&mut |params| visit(params), receiver, args), |params| {
                        params.inner.read(visit)?
                    })
                }
            }
        };
        impl crate::fields::Field for #native {
            type Input = #value_name;
            type Output = #name;
            fn input(value: Self::Input) -> Result<Self, wasm_bindgen::JsValue> {
                Self::decode(value.as_ref())
            }
            fn decode(value: &wasm_bindgen::JsValue) -> Result<Self, wasm_bindgen::JsValue> {
                #name::from_instance(value)
            }
            fn get(inner: crate::storage::NativeField<Self>) -> Result<Self::Output, wasm_bindgen::JsValue> {
                Ok(#name { inner })
            }
        }
        #[wasm_bindgen::prelude::wasm_bindgen]
        impl #name {
            #[wasm_bindgen(constructor)]
            pub fn new(options: #constructor_input) -> Result<Self, wasm_bindgen::JsValue> {
                #construction
                Ok(Self { inner: crate::storage::NativeField::owned(native) })
            }
            pub fn copy(&self) -> Result<Self, wasm_bindgen::JsValue> {
                Ok(Self { inner: crate::storage::NativeField::owned(self.inner.read(Clone::clone)?) })
            }
            #(
                #(#field_docs)*
                #[wasm_bindgen(getter, js_name = #js_names)]
                pub fn #names(&self) -> Result<<#types as crate::fields::Field>::Output, wasm_bindgen::JsValue> {
                    #getters
                }
                #[wasm_bindgen(setter, js_name = #js_names)]
                pub fn #setters(&self, value: <#types as crate::fields::Field>::Input) -> Result<(), wasm_bindgen::JsValue> {
                    let value = <#types as crate::fields::Field>::input(value)?;
                    self.inner.write(|parent| { #assignments; })
                }
            )*
        }
    }.into()
}

struct FlattenDeclaration {
    parent: Declaration,
    field: Ident,
    children: Punctuated<Field, Token![,]>,
}
impl Parse for FlattenDeclaration {
    fn parse(input: ParseStream) -> syn::Result<Self> {
        let parent = input.parse()?;
        input.parse::<Token![,]>()?;
        let field = input.parse()?;
        input.parse::<Token![,]>()?;
        let children;
        braced!(children in input);
        Ok(Self {
            parent,
            field,
            children: children.parse_terminated(Field::parse_named, Token![,])?,
        })
    }
}

/// Resume canonical nested metadata expansion; no field/type/default list is
/// maintained here. Native paths remain visible in generated accessors.
#[proc_macro]
pub fn flatten_parameters(input: TokenStream) -> TokenStream {
    let FlattenDeclaration {
        parent:
            Declaration {
                docs,
                name,
                native,
                fields,
                factory,
            },
        field,
        children,
    } = syn::parse_macro_input!(input as FlattenDeclaration);
    let mut expanded = Vec::new();
    for definition in fields {
        if definition.ident.as_ref() != Some(&field) {
            expanded.push(definition);
            continue;
        }
        let prefix = definition
            .attrs
            .iter()
            .find(|a| a.path().is_ident("projection_path"))
            .map(|a| a.meta.require_list().unwrap().tokens.clone())
            .unwrap_or_else(|| quote!(#field));
        let accessor = definition
            .attrs
            .iter()
            .any(|a| a.path().is_ident("projection_accessors"));
        for mut child in children.clone() {
            let ident = child.ident.as_ref().unwrap();
            if accessor {
                let setter = format_ident!("set_{}", ident);
                child
                    .attrs
                    .push(syn::parse_quote!(#[projection_read(#ident())]));
                child
                    .attrs
                    .push(syn::parse_quote!(#[projection_write(#setter)]));
            } else {
                child
                    .attrs
                    .push(syn::parse_quote!(#[projection_path(#prefix.#ident)]));
            }
            expanded.push(child);
        }
    }
    let mut names = std::collections::HashSet::new();
    for field in &expanded {
        if !names.insert(field.ident.as_ref().unwrap().to_string()) {
            return syn::Error::new_spanned(
                field,
                "flattened configuration fields must have unique names",
            )
            .into_compile_error()
            .into();
        }
    }
    let factory_marker = if factory { quote!(factory) } else { quote!() };
    quote!(cosmolkit_tsify_macros::parameters!(#(#docs)* #name,#native{#(#expanded),*} #factory_marker);)
        .into()
}
