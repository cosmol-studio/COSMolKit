use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::fs;
use std::path::PathBuf;

use quote::{ToTokens, format_ident, quote};
use syn::visit_mut::{self, VisitMut};
use syn::{Attribute, Fields, File, Item, ItemStruct, Type, TypeBareFn, TypePath, parse_quote};

struct OwnedTypeRewrite;

impl OwnedTypeRewrite {
    fn replacement_path(path: &TypePath) -> Option<Type> {
        let ident = path.path.segments.last()?.ident.to_string();
        let replacement = match ident.as_str() {
            "c_schar" | "c_char" => parse_quote!(i8),
            "c_uchar" => parse_quote!(u8),
            "c_short" => parse_quote!(i16),
            "c_ushort" => parse_quote!(u16),
            "c_int" => parse_quote!(i32),
            "c_uint" => parse_quote!(u32),
            "c_long" | "c_longlong" | "__clock_t" | "__off_t" | "__off64_t" => {
                parse_quote!(i64)
            }
            "c_ulong" | "c_ulonglong" => parse_quote!(u64),
            "c_float" => parse_quote!(f32),
            "c_double" => parse_quote!(f64),
            "c_void" => parse_quote!(SourceVoid),
            "_IO_FILE" => parse_quote!(SourceFile),
            "__builtin_va_list" | "__gnuc_va_list" => parse_quote!(SourceVaList),
            _ => return None,
        };
        Some(replacement)
    }

    fn rewrite_bare_function(function: &mut TypeBareFn) {
        function.unsafety = None;
        function.abi = None;
    }
}

impl VisitMut for OwnedTypeRewrite {
    fn visit_type_mut(&mut self, node: &mut Type) {
        match node {
            Type::Ptr(pointer) => {
                let element = (*pointer.elem).clone();
                let mut replacement: Type = if pointer.mutability.is_some() {
                    parse_quote!(SourceMutPointer<#element>)
                } else {
                    parse_quote!(SourceConstPointer<#element>)
                };
                self.visit_type_mut(&mut replacement);
                *node = replacement;
            }
            Type::Path(path) => {
                if let Some(mut replacement) = Self::replacement_path(path) {
                    self.visit_type_mut(&mut replacement);
                    *node = replacement;
                } else {
                    visit_mut::visit_type_path_mut(self, path);
                }
            }
            Type::BareFn(function) => {
                Self::rewrite_bare_function(function);
                visit_mut::visit_type_bare_fn_mut(self, function);
            }
            _ => visit_mut::visit_type_mut(self, node),
        }
    }
}

fn item_name(item: &Item) -> Option<String> {
    match item {
        Item::Const(item) => Some(item.ident.to_string()),
        Item::Struct(item) => Some(item.ident.to_string()),
        Item::Type(item) => Some(item.ident.to_string()),
        _ => None,
    }
}

fn skipped_name(name: &str) -> bool {
    name == "__BindgenUnionField"
        || name == "tagSplitLong"
        || name == "BnsAltPath"
        || name == "__builtin_va_list"
        || name == "__gnuc_va_list"
        || name == "__va_list_tag"
        || name.starts_with("_IO_")
        || name.starts_with("__off")
        || name == "__clock_t"
}

fn owned_derives() -> Vec<Attribute> {
    vec![parse_quote!(#[derive(Clone, Debug, PartialEq)])]
}

fn normalize_struct(item: &mut ItemStruct) {
    item.attrs = owned_derives();
}

fn default_expression(field_type: &Type) -> proc_macro2::TokenStream {
    match field_type {
        Type::Array(array) => {
            let element = default_expression(&array.elem);
            quote!(::std::array::from_fn(|_| #element))
        }
        _ => quote!(::std::default::Default::default()),
    }
}

fn default_implementation(item: &ItemStruct) -> Item {
    let name = &item.ident;
    let constructor = match &item.fields {
        Fields::Named(fields) => {
            let initializers = fields.named.iter().map(|field| {
                let field_name = field.ident.as_ref().expect("named field");
                let value = default_expression(&field.ty);
                quote!(#field_name: #value)
            });
            quote!(Self { #(#initializers),* })
        }
        Fields::Unnamed(fields) => {
            let initializers = fields
                .unnamed
                .iter()
                .map(|field| default_expression(&field.ty));
            quote!(Self(#(#initializers),*))
        }
        Fields::Unit => quote!(Self),
    };
    syn::parse2(quote! {
        impl ::std::default::Default for #name {
            fn default() -> Self {
                #constructor
            }
        }
    })
    .expect("generate Default implementation")
}

fn transform(mut file: File) -> File {
    file.items.retain(|item| {
        !matches!(item, Item::ForeignMod(_) | Item::Impl(_))
            && item_name(item).is_none_or(|name| !skipped_name(&name))
    });

    let mut rewrite = OwnedTypeRewrite;
    let mut transformed_items = Vec::new();
    for mut item in file.items {
        match item {
            Item::Struct(ref mut item) => normalize_struct(item),
            Item::Type(ref mut item) => item.attrs.clear(),
            Item::Const(ref mut item) => item.attrs.clear(),
            _ => {}
        }
        rewrite.visit_item_mut(&mut item);
        let default = match &item {
            Item::Struct(item) => Some(default_implementation(item)),
            _ => None,
        };
        transformed_items.push(item);
        if let Some(default) = default {
            transformed_items.push(default);
        }
    }
    file.items = transformed_items;
    file
}

fn deduplicate(items: Vec<Item>) -> Vec<Item> {
    let mut named = BTreeMap::<String, String>::new();
    let mut unnamed = BTreeSet::<String>::new();
    let mut deduplicated = Vec::new();
    for item in items {
        let rendered = item.to_token_stream().to_string();
        if let Some(name) = item_name(&item) {
            if let Some(previous) = named.get(&name) {
                assert_eq!(
                    previous, &rendered,
                    "conflicting configured declarations named {name}"
                );
                continue;
            }
            named.insert(name, rendered);
        } else if !unnamed.insert(rendered) {
            continue;
        }
        deduplicated.push(item);
    }
    deduplicated
}

fn main() {
    let mut paths: Vec<PathBuf> = env::args_os().skip(1).map(PathBuf::from).collect();
    assert!(
        paths.len() >= 2,
        "expected one or more inputs and one output"
    );
    let output = paths.pop().expect("output Rust path");

    let mut inputs = paths.into_iter();
    let header_input = inputs.next().expect("configured header input");
    let header_source = fs::read_to_string(&header_input).expect("read header bindings");
    let header_file = syn::parse_file(&header_source).expect("parse header bindings");
    let mut items = deduplicate(transform(header_file).items);

    for input in inputs {
        let source = fs::read_to_string(&input).expect("read bindgen Rust input");
        let parsed = syn::parse_file(&source).expect("parse bindgen Rust input");
        let local_items = deduplicate(transform(parsed).items);
        if local_items.is_empty() {
            continue;
        }
        let stem = input.file_stem().expect("input filename").to_string_lossy();
        let stem = stem
            .split_once('_')
            .map_or(stem.as_ref(), |(_, source_stem)| source_stem);
        let module = format_ident!("local_{stem}");
        items.push(parse_quote! {
            pub(crate) mod #module {
                use super::*;
                #(#local_items)*
            }
        });
    }
    let transformed = File {
        shebang: None,
        attrs: Vec::new(),
        items,
    };
    let body = prettyplease::unparse(&transformed);
    let generated = format!(
        "// Generated from IUPAC InChI v1.07.5 at commit\n\
         // 11a87982bb518f57ac013f0b258c283655e1ea1d.\n\
         // Input: dev/inchi_type_inventory_wrapper.h plus all 60 production C files.\n\
         // Configuration: COMPILE_ANSI_ONLY, TARGET_API_LIB, GCC/Clang LP64.\n\
         // Generator declaration baseline: bindgen-cli 0.72.1.\n\
         // The checked-in result is pure Rust and contains no FFI declarations.\n\
         use super::*;\n\n{body}"
    );
    fs::write(output, generated).expect("write owned Rust source types");

    let rendered = transformed.to_token_stream().to_string();
    assert!(!rendered.contains("* mut"));
    assert!(!rendered.contains("* const"));
    assert!(!rendered.contains("extern \"C\""));
}
