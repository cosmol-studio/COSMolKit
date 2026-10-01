//! Lightweight BIO transactions: field references, COW and one body per pair.
use proc_macro2::TokenStream;
use quote::{format_ident, quote};
use std::collections::BTreeSet;
use syn::{
    Ident, Path, Token, Type, braced, bracketed, parenthesized,
    parse::{Parse, ParseStream},
    punctuated::Punctuated,
};

pub struct Registry(Vec<Operation>);
struct Field {
    name: Ident,
    ty: Type,
}
impl Parse for Field {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        let name = input.parse()?;
        input.parse::<Token![:]>()?;
        Ok(Self {
            name,
            ty: input.parse()?,
        })
    }
}
struct Operation {
    name: Ident,
    args: Vec<Field>,
    targets: Vec<Ident>,
    value: Ident,
    inplace: Ident,
    body: Path,
    access: Ident,
    read: Vec<Field>,
    write: Vec<Field>,
    replace: Vec<Field>,
}
fn key(input: ParseStream<'_>, expected: &str) -> syn::Result<()> {
    let found: Ident = input.parse()?;
    if found != expected {
        return Err(syn::Error::new(
            found.span(),
            format!("expected {expected}"),
        ));
    }
    input.parse::<Token![:]>()?;
    Ok(())
}
fn fields(input: ParseStream<'_>) -> syn::Result<Vec<Field>> {
    let inner;
    bracketed!(inner in input);
    let fields = Punctuated::<Field, Token![,]>::parse_terminated(&inner)?;
    input.parse::<Token![,]>()?;
    Ok(fields.into_iter().collect())
}
impl Parse for Registry {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        let mut ops = Vec::new();
        let mut methods = BTreeSet::new();
        let mut ids = BTreeSet::new();
        let mut accesses = BTreeSet::new();
        while !input.is_empty() {
            let name: Ident = input.parse()?;
            let args;
            parenthesized!(args in input);
            let args = Punctuated::<Field, Token![,]>::parse_terminated(&args)?
                .into_iter()
                .collect::<Vec<_>>();
            let inner;
            braced!(inner in input);
            key(&inner, "targets")?;
            let targets;
            bracketed!(targets in inner);
            let targets = Punctuated::<Ident, Token![,]>::parse_terminated(&targets)?
                .into_iter()
                .collect::<Vec<_>>();
            inner.parse::<Token![,]>()?;
            key(&inner, "value")?;
            let value: Ident = inner.parse()?;
            inner.parse::<Token![,]>()?;
            key(&inner, "inplace")?;
            let inplace: Ident = inner.parse()?;
            inner.parse::<Token![,]>()?;
            key(&inner, "body")?;
            let body = inner.parse()?;
            inner.parse::<Token![,]>()?;
            key(&inner, "access")?;
            let access: Ident = inner.parse()?;
            inner.parse::<Token![,]>()?;
            key(&inner, "read")?;
            let read = fields(&inner)?;
            key(&inner, "write")?;
            let write = fields(&inner)?;
            let replace = if inner.is_empty() {
                Vec::new()
            } else {
                key(&inner, "replace")?;
                fields(&inner)?
            };
            if !inner.is_empty() {
                return Err(inner.error("unexpected BIO operation field"));
            }
            if targets.is_empty() || read.is_empty() && write.is_empty() && replace.is_empty() {
                return Err(syn::Error::new(
                    name.span(),
                    "BIO operation requires targets and field access",
                ));
            }
            if !ids.insert(name.to_string()) || !accesses.insert(access.to_string()) {
                return Err(syn::Error::new(
                    name.span(),
                    "duplicate BIO operation or access type",
                ));
            }
            if value.to_string().ends_with('_') || !inplace.to_string().ends_with('_') {
                return Err(syn::Error::new(
                    name.span(),
                    "only the in-place method must end in underscore",
                ));
            }
            let mut names = BTreeSet::new();
            for f in read.iter().chain(write.iter()).chain(replace.iter()) {
                if !names.insert(f.name.to_string()) {
                    return Err(syn::Error::new(
                        f.name.span(),
                        "duplicate or overlapping BIO field access",
                    ));
                }
            }
            let mut arg_names = BTreeSet::new();
            for arg in &args {
                if !arg_names.insert(arg.name.to_string())
                    || ["self", "working", "data", "access"]
                        .contains(&arg.name.to_string().as_str())
                {
                    return Err(syn::Error::new(
                        arg.name.span(),
                        "duplicate or reserved argument name",
                    ));
                }
            }
            for target in &targets {
                if target != "BioStructure" && target != "Protein" {
                    return Err(syn::Error::new(
                        target.span(),
                        "BIO target must be BioStructure or Protein",
                    ));
                }
                for method in [&value, &inplace] {
                    if !methods.insert((target.to_string(), method.to_string())) {
                        return Err(syn::Error::new(
                            method.span(),
                            "duplicate BIO target method",
                        ));
                    }
                }
            }
            ops.push(Operation {
                name,
                args,
                targets,
                value,
                inplace,
                body,
                access,
                read,
                write,
                replace,
            });
        }
        Ok(Self(ops))
    }
}

pub fn expand(registry: Registry) -> TokenStream {
    let mut definitions = Vec::new();
    let mut entries = Vec::new();
    for op in registry.0 {
        let Operation {
            name,
            args,
            targets,
            value,
            inplace,
            body,
            access,
            read,
            write,
            replace,
        } = op;
        let rn: Vec<_> = read.iter().map(|f| &f.name).collect();
        let rt: Vec<_> = read.iter().map(|f| &f.ty).collect();
        let wn: Vec<_> = write.iter().map(|f| &f.name).collect();
        let wt: Vec<_> = write.iter().map(|f| &f.ty).collect();
        let xn: Vec<_> = replace.iter().map(|f| &f.name).collect();
        let xt: Vec<_> = replace.iter().map(|f| &f.ty).collect();
        let replacement = format_ident!("{}Replacement", access);
        let an: Vec<_> = args.iter().map(|f| &f.name).collect();
        let at: Vec<_> = args.iter().map(|f| &f.ty).collect();
        definitions.push(quote! {
            pub(super) struct #access<'a> {
                #(pub(super) #rn: &'a #rt,)*
                #(pub(super) #wn: &'a mut #wt,)*
                #(pub(super) #xn: &'a std::sync::Arc<#xt>,)*
            }
        });
        if !replace.is_empty() {
            definitions.push(quote! {
                pub(super) struct #replacement {
                    #(pub(super) #xn: std::sync::Arc<#xt>,)*
                }
            });
        }
        // Replacement permissions produce a typed output, not mutable access
        // to old blocks. Installation stays here, within the storage module;
        // no body receives a setter or a complete detached input object.
        let invoke = quote! {
            let data = working.operation_data_mut();
            let access = #access {
                #(#rn: &data.#rn,)*
                #(#wn: std::sync::Arc::make_mut(&mut data.#wn),)*
                #(#xn: &data.#xn,)*
            };
            #body(access, #(#an),*)?
        };
        let execute = if replace.is_empty() {
            quote! { { #invoke; } }
        } else {
            quote! {
                let replacement: #replacement = { #invoke };
                {
                    let data = working.operation_data_mut();
                    #(data.#xn = replacement.#xn;)*
                }
            }
        };
        for target in targets {
            definitions.push(quote! {
                impl #target {
                    /// Return a new value; unchanged blocks remain shared.
                    pub fn #value(&self, #(#an: #at),*) -> Result<Self, crate::BioOperationError> {
                        let mut working = self.clone();
                        #execute
                        working.validate_operation()?;
                        Ok(working)
                    }
                    /// Replace this value only after the common body and validation succeed.
                    pub fn #inplace(&mut self, #(#an: #at),*) -> Result<(), crate::BioOperationError> {
                        let result = self.#value(#(#an),*)?;
                        *self = result;
                        Ok(())
                    }
                }
            });
            entries.push(quote! { BioOperationSpec {
                id: stringify!(#name), target: stringify!(#target),
                value_method: stringify!(#value), inplace_method: stringify!(#inplace),
                read: &[#(stringify!(#rn)),*], write: &[#(stringify!(#wn)),*],
                replace: &[#(stringify!(#xn)),*],
            } });
        }
    }
    quote! {
        #(#definitions)*
        pub(super) const BIO_STRUCTURE_OPS: &[BioOperationSpec] = &[#(#entries),*];
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const VALID: &str = "translate(offset: [f64; 3]) { targets: [BioStructure, Protein], value: with_translation, inplace: translate_, body: super::translate, access: TranslationAccess, read: [], write: [coordinates: Coordinates], }";
    #[test]
    fn valid_pair_has_one_body_call_per_target() {
        let output = expand(syn::parse_str(VALID).unwrap()).to_string();
        assert_eq!(output.matches("super :: translate").count(), 2);
        assert_eq!(output.matches("Arc :: make_mut").count(), 2);
        assert_eq!(output.matches("validate_operation").count(), 2);
    }
    #[test]
    fn invalid_permissions_names_and_duplicate_targets_are_rejected() {
        for bad in [
            VALID.replace("read: []", "read: [coordinates: Coordinates]"),
            VALID.replace("translate_,", "translate,"),
            VALID.replace("BioStructure, Protein", "BioStructure, BioStructure"),
            VALID.replace("BioStructure, Protein", "Molecule"),
            VALID.replace("offset:", "working:"),
        ] {
            assert!(syn::parse_str::<Registry>(&bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn replacement_is_declared_generic_and_never_materializes_old_blocks() {
        let declaration = VALID.replace(
            "write: [coordinates: Coordinates],",
            "write: [], replace: [coordinates: Coordinates, atoms: Vec<Atom>],",
        );
        let output = expand(syn::parse_str(&declaration).unwrap()).to_string();
        assert!(output.contains("TranslationAccessReplacement"));
        assert!(output.contains("coordinates : & 'a std :: sync :: Arc < Coordinates >"));
        assert!(output.contains("atoms : std :: sync :: Arc < Vec < Atom > >"));
        assert_eq!(
            output
                .matches("data . coordinates = replacement . coordinates")
                .count(),
            2
        );
        assert_eq!(
            output.matches("data . atoms = replacement . atoms").count(),
            2
        );
        assert_eq!(output.matches("Arc :: make_mut").count(), 0);
        assert_eq!(output.matches("super :: translate").count(), 2);
        assert_eq!(output.matches("validate_operation").count(), 2);
    }

    #[test]
    fn replacement_permissions_must_be_disjoint_and_nonduplicated() {
        let good = VALID.replace(
            "write: [coordinates: Coordinates],",
            "write: [], replace: [coordinates: Coordinates],",
        );
        for bad in [
            good.replace("read: []", "read: [coordinates: Coordinates]"),
            good.replace("write: []", "write: [coordinates: Coordinates]"),
            good.replace(
                "replace: [coordinates: Coordinates]",
                "replace: [coordinates: Coordinates, coordinates: Coordinates]",
            ),
        ] {
            assert!(syn::parse_str::<Registry>(&bad).is_err(), "{bad}");
        }
    }
}
