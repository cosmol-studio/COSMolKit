//! Feature-local declaration strings, without postprocessing generated .d.ts.
use proc_macro::TokenStream;
use quote::{format_ident, quote};
use syn::{
    Attribute, Ident, LitStr, Token,
    parse::{Parse, ParseStream},
};

struct Union {
    name: Ident,
    parts: Vec<(Vec<Attribute>, LitStr)>,
}

impl Parse for Union {
    fn parse(input: ParseStream) -> syn::Result<Self> {
        let name = input.parse()?;
        input.parse::<Token![;]>()?;
        let mut parts = Vec::new();
        while !input.is_empty() {
            let attrs = input.call(Attribute::parse_outer)?;
            for attr in &attrs {
                if !attr.path().is_ident("cfg") {
                    return Err(syn::Error::new_spanned(attr, "only cfg is supported"));
                }
            }
            parts.push((attrs, input.parse()?));
            input.parse::<Token![;]>()?;
        }
        Ok(Self { name, parts })
    }
}

pub fn expand(input: TokenStream) -> TokenStream {
    let Union { name, parts } = syn::parse_macro_input!(input as Union);
    let public_name = name.to_string();
    let section = format_ident!("__TS_UNION_{}", public_name.to_uppercase());
    let members = parts.iter().map(|(attrs, text)| {
        quote! {
            #(#attrs)* concat!(" | ", #text),
        }
    });
    quote! {
        #[wasm_bindgen::prelude::wasm_bindgen]
        extern "C" {
            #[wasm_bindgen(typescript_type = #public_name)]
            #[derive(Clone)]
            pub type #name;
        }
        #[wasm_bindgen::prelude::wasm_bindgen(typescript_custom_section)]
        const #section: &str = {
            const PARTS: &[&str] = &[
                concat!("export type ", #public_name, " = never"),
                #(#members)*
                ";",
            ];
            const LENGTH: usize = {
                let mut length = 0;
                let mut i = 0;
                while i < PARTS.len() { length += PARTS[i].len(); i += 1; }
                length
            };
            const BYTES: [u8; LENGTH] = {
                let mut bytes = [0; LENGTH];
                let mut offset = 0;
                let mut i = 0;
                while i < PARTS.len() {
                    let part = PARTS[i].as_bytes();
                    let mut j = 0;
                    while j < part.len() {
                        bytes[offset] = part[j];
                        offset += 1; j += 1;
                    }
                    i += 1;
                }
                bytes
            };
            // All fragments are UTF-8 string literals; concatenation preserves UTF-8.
            match str::from_utf8(&BYTES) {
                Ok(text) => text,
                Err(_) => panic!("invalid declaration text"),
            }
        };
    }
    .into()
}
