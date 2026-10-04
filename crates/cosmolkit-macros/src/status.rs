//! Shared parsing for the four behavior commitments. No execution authority.

use quote::quote;
use syn::{
    Ident, LitStr, Token, parenthesized,
    parse::{Parse, ParseStream},
};

#[derive(Clone)]
pub(crate) enum FunctionStatus {
    Parity(LitStr),
    ParityWithDifferences(LitStr, LitStr),
    Native,
    Experimental,
}

impl Default for FunctionStatus {
    fn default() -> Self {
        Self::Experimental
    }
}

impl Parse for FunctionStatus {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        let name: Ident = input.parse()?;
        match name.to_string().as_str() {
            "experimental" => Ok(Self::Experimental),
            "native" => Ok(Self::Native),
            "parity" | "parity_with_differences" => {
                let content;
                parenthesized!(content in input);
                let reference: LitStr = content.parse()?;
                if reference.value().trim().is_empty() {
                    return Err(syn::Error::new_spanned(
                        reference,
                        "parity reference must not be empty",
                    ));
                }
                let status = if name == "parity" {
                    Self::Parity(reference)
                } else {
                    content.parse::<Token![,]>()?;
                    let explanation: LitStr = content.parse()?;
                    if explanation.value().trim().is_empty() {
                        return Err(syn::Error::new_spanned(
                            explanation,
                            "approved difference explanation must not be empty",
                        ));
                    }
                    Self::ParityWithDifferences(reference, explanation)
                };
                if !content.is_empty() {
                    return Err(content.error("unexpected function status arguments"));
                }
                Ok(status)
            }
            _ => Err(syn::Error::new_spanned(
                name,
                "status must be parity, parity_with_differences, native or experimental",
            )),
        }
    }
}

impl FunctionStatus {
    pub(crate) fn tokens(&self) -> proc_macro2::TokenStream {
        match self {
            Self::Parity(reference) => {
                quote!(crate::FunctionStatus::Parity { reference: #reference })
            }
            Self::ParityWithDifferences(reference, explanation) => {
                quote!(crate::FunctionStatus::ParityWithDifferences { reference: #reference, explanation: #explanation })
            }
            Self::Native => quote!(crate::FunctionStatus::Native),
            Self::Experimental => quote!(crate::FunctionStatus::Experimental),
        }
    }
}
