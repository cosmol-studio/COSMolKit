//! Private operation-body capability projection.
//!
//! Generated tokens name runtime-owned types at the call site. This module
//! neither defines runtime state nor depends on a runtime or domain crate.

use std::collections::HashMap;

use proc_macro::TokenStream;
use quote::{format_ident, quote};
use syn::{FnArg, Ident, ItemFn, Pat, Token, Type, parse::Parse, parse::ParseStream, parse_quote};

use crate::declaration::{
    BioBlock, BioDerivedState, BioOperation, BioRegistry, BioState, MappingRequirement,
    MoleculeBlock, MoleculeOperation, MoleculeOutput, MoleculeRegistry, TopologyEditKind,
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum BodyClass {
    MoleculeSingle,
    MoleculeMultiple,
    Bio,
}

impl BodyClass {
    const fn discriminator(self) -> usize {
        match self {
            Self::MoleculeSingle => 0,
            Self::MoleculeMultiple => 1,
            Self::Bio => 2,
        }
    }

    const fn macro_name(self) -> &'static str {
        match self {
            Self::MoleculeSingle => "mol_op_body",
            Self::MoleculeMultiple => "mol_multi_op_body",
            Self::Bio => "bio_op_body",
        }
    }
}

pub(crate) struct BodyAttribute {
    pub(crate) operation: Ident,
    _comma: Token![,],
    pub(crate) context: Ident,
}

impl Parse for BodyAttribute {
    fn parse(input: ParseStream<'_>) -> syn::Result<Self> {
        Ok(Self {
            operation: input.parse()?,
            _comma: input.parse()?,
            context: input.parse()?,
        })
    }
}

pub(crate) fn access_marker(operation: &Ident) -> syn::Result<Ident> {
    let operation_name = operation.to_string();
    if operation_name.starts_with("r#") {
        return Err(syn::Error::new_spanned(
            operation,
            "operation names used by body macros must not be raw identifiers",
        ));
    }

    let mut marker_name = String::new();
    for component in operation_name.split('_') {
        if component.is_empty() {
            return Err(syn::Error::new_spanned(
                operation,
                "operation names used by body macros must not contain empty snake-case components",
            ));
        }
        let mut chars = component.chars();
        let first = chars.next().expect("nonempty component checked above");
        marker_name.extend(first.to_uppercase());
        marker_name.extend(chars);
    }
    marker_name.push_str("Access");
    Ok(format_ident!("{marker_name}", span = operation.span()))
}

fn context_type(class: BodyClass, marker: &Ident) -> syn::Result<Type> {
    match class {
        BodyClass::MoleculeSingle => syn::parse2(quote!(&mut crate::OpParts<'_, crate::#marker>)),
        BodyClass::MoleculeMultiple => {
            syn::parse2(quote!(&mut crate::MultiOutputOpParts<'_, crate::#marker>))
        }
        BodyClass::Bio => syn::parse2(quote!(&mut crate::BioOpParts<'_, crate::#marker>)),
    }
}

fn reject_receiver_or_collision(
    function: &ItemFn,
    context: &Ident,
    class: BodyClass,
) -> syn::Result<()> {
    for input in &function.sig.inputs {
        match input {
            FnArg::Receiver(_) => {
                return Err(syn::Error::new_spanned(
                    input,
                    format!(
                        "{} operation bodies must not receive self",
                        class.macro_name()
                    ),
                ));
            }
            FnArg::Typed(argument) => {
                if let Pat::Ident(binding) = argument.pat.as_ref()
                    && binding.ident == *context
                {
                    return Err(syn::Error::new_spanned(
                        &binding.ident,
                        format!(
                            "{} context name `{context}` conflicts with an existing function parameter",
                            class.macro_name()
                        ),
                    ));
                }
            }
        }
    }
    Ok(())
}

pub(crate) fn expand_body_tokens(
    attribute: proc_macro2::TokenStream,
    item: proc_macro2::TokenStream,
    class: BodyClass,
) -> syn::Result<proc_macro2::TokenStream> {
    let BodyAttribute {
        operation, context, ..
    } = syn::parse2(attribute)?;
    let mut function: ItemFn = syn::parse2(item)?;
    reject_receiver_or_collision(&function, &context, class)?;

    let marker = access_marker(&operation)?;
    let context_type = context_type(class, &marker)?;
    let original_parameters = std::mem::take(&mut function.sig.inputs);
    function
        .sig
        .inputs
        .push(parse_quote!(#context: #context_type));
    function.sig.inputs.extend(original_parameters);

    let expected_class = class.discriminator();
    Ok(quote! {
        const _: [(); #expected_class] = [(); crate::#marker::__COSMOLKIT_BODY_CLASS];
        #function
    })
}

pub(crate) fn expand_body(
    attribute: TokenStream,
    item: TokenStream,
    class: BodyClass,
) -> TokenStream {
    match expand_body_tokens(attribute.into(), item.into(), class) {
        Ok(output) => output.into(),
        Err(error) => error.to_compile_error().into(),
    }
}

fn validate_unique_markers<'a>(operations: impl IntoIterator<Item = &'a Ident>) -> syn::Result<()> {
    let mut seen = HashMap::<String, Ident>::new();
    for operation in operations {
        let marker = access_marker(operation)?;
        let marker_name = marker.to_string();
        if let Some(previous) = seen.insert(marker_name.clone(), operation.clone()) {
            let mut error = syn::Error::new_spanned(
                operation,
                format!(
                    "operation `{operation}` generates access marker `{marker_name}`, which conflicts with operation `{previous}`"
                ),
            );
            error.combine(syn::Error::new_spanned(
                previous,
                format!("first operation generating `{marker_name}` is here"),
            ));
            return Err(error);
        }
    }
    Ok(())
}

const fn molecule_block_bit(block: MoleculeBlock) -> u64 {
    match block {
        MoleculeBlock::Topology => 1 << 0,
        MoleculeBlock::Coordinates => 1 << 1,
        MoleculeBlock::Properties => 1 << 2,
        MoleculeBlock::DerivedCache => 1 << 3,
    }
}

fn molecule_block_mask(blocks: &[MoleculeBlock]) -> u64 {
    blocks
        .iter()
        .fold(0, |mask, block| mask | molecule_block_bit(*block))
}

const fn bio_block_bit(block: BioBlock) -> u64 {
    match block {
        BioBlock::Atoms => 1 << 0,
        BioBlock::Residues => 1 << 1,
        BioBlock::Chains => 1 << 2,
        BioBlock::Entities => 1 << 3,
        BioBlock::Models => 1 << 4,
        BioBlock::Coordinates => 1 << 5,
        BioBlock::Bonds => 1 << 6,
        BioBlock::Assemblies => 1 << 7,
        BioBlock::Annotations => 1 << 8,
        BioBlock::DerivedCache => 1 << 9,
        BioBlock::Properties => 1 << 10,
    }
}

fn bio_block_mask(blocks: &[BioBlock]) -> u64 {
    blocks
        .iter()
        .fold(0, |mask, block| mask | bio_block_bit(*block))
}

const fn bio_state_bit(state: BioState) -> u64 {
    match state {
        BioState::Hierarchy => 1 << 0,
        BioState::ResidueSpans => 1 << 1,
        BioState::ChainSpans => 1 << 2,
        BioState::ModelSpans => 1 << 3,
        BioState::CoordinateAlignment => 1 << 4,
        BioState::EntityMapping => 1 << 5,
        BioState::AltlocGroups => 1 << 6,
        BioState::AssemblyReferences => 1 << 7,
        BioState::BondReferences => 1 << 8,
        BioState::SelectionProvenance => 1 << 9,
        BioState::PolymerAnnotation => 1 << 10,
        BioState::SecondaryStructure => 1 << 11,
    }
}

fn bio_state_mask(states: &[BioState]) -> u64 {
    states
        .iter()
        .fold(0, |mask, state| mask | bio_state_bit(*state))
}

const fn bio_derived_state_bit(state: BioDerivedState) -> u64 {
    match state {
        BioDerivedState::AtomIndex => 1 << 0,
        BioDerivedState::ResidueIndex => 1 << 1,
        BioDerivedState::ChainIndex => 1 << 2,
        BioDerivedState::EntityIndex => 1 << 3,
        BioDerivedState::SequenceCache => 1 << 4,
        BioDerivedState::PolymerCache => 1 << 5,
        BioDerivedState::AltlocCache => 1 << 6,
        BioDerivedState::AssemblyCache => 1 << 7,
        BioDerivedState::BondCache => 1 << 8,
        BioDerivedState::BackboneGeometry => 1 << 9,
        BioDerivedState::SidechainGeometry => 1 << 10,
        BioDerivedState::NucleicGeometry => 1 << 11,
        BioDerivedState::SecondaryStructure => 1 << 12,
        BioDerivedState::ContactMap => 1 << 13,
        BioDerivedState::GraphCache => 1 << 14,
    }
}

fn bio_derived_state_mask(states: &[BioDerivedState]) -> u64 {
    states
        .iter()
        .fold(0, |mask, state| mask | bio_derived_state_bit(*state))
}

fn expand_molecule_marker(operation: &MoleculeOperation) -> syn::Result<proc_macro2::TokenStream> {
    let cfg_attrs = &operation.cfg_attrs;
    let marker = access_marker(&operation.name)?;
    let class = match operation.fields.output {
        MoleculeOutput::Single => BodyClass::MoleculeSingle,
        MoleculeOutput::Multiple => BodyClass::MoleculeMultiple,
    }
    .discriminator();
    let read = molecule_block_mask(&operation.fields.access.read);
    let write = molecule_block_mask(&operation.fields.access.write);
    let mut methods = Vec::new();

    for block in &operation.fields.access.read {
        methods.push(match block {
            MoleculeBlock::Topology => quote! {
                pub(crate) fn topology(&self) -> Result<&cosmolkit_model::TopologyBlock, crate::OperationError> {
                    self.read_topology_runtime()
                }
            },
            MoleculeBlock::Coordinates => quote! {
                pub(crate) fn coordinates(&self) -> Result<&cosmolkit_model::CoordinateBlock, crate::OperationError> {
                    self.read_coordinates_runtime()
                }
            },
            MoleculeBlock::Properties => quote! {
                pub(crate) fn properties(&self) -> Result<&cosmolkit_model::MoleculeProperties, crate::OperationError> {
                    self.read_properties_runtime()
                }
            },
            MoleculeBlock::DerivedCache => quote! {
                pub(crate) fn derived_cache(&self) -> Result<&crate::molecule::DerivedCacheBlock, crate::OperationError> {
                    self.read_derived_cache_runtime()
                }
            },
        });
    }

    for block in &operation.fields.access.write {
        methods.extend(match block {
            MoleculeBlock::Topology => vec![
                quote! {
                    pub(crate) fn checkout_topology(&mut self) -> Result<cosmolkit_model::TopologyBlock, crate::OperationError> {
                        self.checkout_topology_runtime()
                    }
                },
                quote! {
                    pub(crate) fn install_topology(&mut self, value: cosmolkit_model::TopologyBlock) -> Result<(), crate::OperationError> {
                        self.install_topology_runtime(value)
                    }
                },
            ],
            MoleculeBlock::Coordinates => vec![
                quote! {
                    pub(crate) fn checkout_coordinates(&mut self) -> Result<cosmolkit_model::CoordinateBlock, crate::OperationError> {
                        self.checkout_coordinates_runtime()
                    }
                },
                quote! {
                    pub(crate) fn install_coordinates(&mut self, value: cosmolkit_model::CoordinateBlock) -> Result<(), crate::OperationError> {
                        self.install_coordinates_runtime(value)
                    }
                },
            ],
            MoleculeBlock::Properties => vec![
                quote! {
                    pub(crate) fn checkout_properties(&mut self) -> Result<cosmolkit_model::MoleculeProperties, crate::OperationError> {
                        self.checkout_properties_runtime()
                    }
                },
                quote! {
                    pub(crate) fn install_properties(&mut self, value: cosmolkit_model::MoleculeProperties) -> Result<(), crate::OperationError> {
                        self.install_properties_runtime(value)
                    }
                },
            ],
            MoleculeBlock::DerivedCache => vec![
                quote! {
                    pub(crate) fn checkout_derived_cache(&mut self) -> Result<crate::molecule::DerivedCacheBlock, crate::OperationError> {
                        self.checkout_derived_cache_runtime()
                    }
                },
                quote! {
                    pub(crate) fn install_derived_cache(&mut self, value: crate::molecule::DerivedCacheBlock) -> Result<(), crate::OperationError> {
                        self.install_derived_cache_runtime(value)
                    }
                },
            ],
        });
    }

    if operation.fields.topology_edit != TopologyEditKind::None {
        methods.push(quote! {
            pub(crate) fn record_topology_edit(&mut self, edit: crate::TopologyEditKind) -> Result<(), crate::OperationError> {
                self.record_topology_edit_runtime(edit)
            }
        });
    }
    if operation.fields.requires_mapping != MappingRequirement::None {
        methods.push(quote! {
            pub(crate) fn record_topology_mapping(&mut self, mapping: cosmolkit_model::TopologyMapping) -> Result<(), crate::OperationError> {
                self.record_topology_mapping_runtime(mapping)
            }
        });
    }
    if !operation.fields.auto_remap.is_empty() {
        methods.push(quote! {
            pub(crate) fn apply_runtime_remap(&mut self) -> Result<(), crate::OperationError> {
                self.apply_runtime_remap_runtime()
            }
        });
    }
    if !operation.fields.derived_effects.recompute.is_empty()
        || !operation
            .fields
            .derived_effects
            .operation_defined
            .is_empty()
    {
        methods.push(quote! {
            pub(crate) fn mark_cache_updated(&mut self, states: crate::DerivedState) -> Result<(), crate::OperationError> {
                self.mark_cache_updated_runtime(states)
            }
        });
    }
    if !operation.fields.derived_effects.recompute.is_empty()
        || !operation.fields.derived_effects.invalidate.is_empty()
        || !operation
            .fields
            .derived_effects
            .operation_defined
            .is_empty()
    {
        methods.push(quote! {
            pub(crate) fn clear_cache(&mut self, states: crate::DerivedState) -> Result<(), crate::OperationError> {
                self.clear_cache_runtime(states)
            }
        });
    }
    if !operation.fields.derived_effects.preserve.is_empty() {
        methods.push(quote! {
            pub(crate) fn prove_preserved(
                &mut self,
                states: crate::DerivedState,
                proof: crate::PreservationProof,
            ) -> Result<(), crate::OperationError> {
                self.prove_preserved_runtime(states, proof)
            }
        });
    }
    methods.push(quote! {
        pub(crate) fn apply_cip_policy(&mut self) -> Result<(), crate::OperationError> {
            self.apply_cip_policy_runtime()
        }
    });

    let capability_impl = match operation.fields.output {
        MoleculeOutput::Single => quote! {
            #(#cfg_attrs)*
            impl<'a> crate::OpParts<'a, #marker> {
                #(#methods)*
            }
        },
        MoleculeOutput::Multiple => {
            let read_methods = operation.fields.access.read.iter().filter_map(|block| match block {
                MoleculeBlock::Topology => Some(quote! {
                    pub(crate) fn topology(&self) -> Result<&cosmolkit_model::TopologyBlock, crate::OperationError> {
                        self.source_topology_runtime()
                    }
                }),
                MoleculeBlock::Coordinates => Some(quote! {
                    pub(crate) fn coordinates(&self) -> Result<&cosmolkit_model::CoordinateBlock, crate::OperationError> {
                        self.source_coordinates_runtime()
                    }
                }),
                MoleculeBlock::Properties => Some(quote! {
                    pub(crate) fn properties(&self) -> Result<&cosmolkit_model::MoleculeProperties, crate::OperationError> {
                        self.source_properties_runtime()
                    }
                }),
                MoleculeBlock::DerivedCache => None,
            });
            quote! {
                #(#cfg_attrs)*
                impl<'a> crate::MultiOutputOpParts<'a, #marker> {
                    #(#read_methods)*

                    pub(crate) fn emit_all(
                        &mut self,
                        candidates: Vec<(
                            cosmolkit_model::TopologyBlock,
                            cosmolkit_model::CoordinateBlock,
                            cosmolkit_model::MoleculeProperties,
                        )>,
                    ) -> Result<(), crate::OperationError> {
                        self.emit_all_runtime(candidates)
                    }
                }
            }
        }
    };
    Ok(quote! {
        #(#cfg_attrs)*
        pub(crate) struct #marker;
        #(#cfg_attrs)*
        impl #marker {
            pub(crate) const __COSMOLKIT_BODY_CLASS: usize = #class;
            pub(crate) const __COSMOLKIT_ACCESS_READ: u64 = #read;
            pub(crate) const __COSMOLKIT_ACCESS_WRITE: u64 = #write;
        }
        #capability_impl
    })
}

fn expand_bio_marker(operation: &BioOperation) -> syn::Result<proc_macro2::TokenStream> {
    let cfg_attrs = &operation.cfg_attrs;
    let marker = access_marker(&operation.name)?;
    let class = BodyClass::Bio.discriminator();
    let may_mutate = bio_block_mask(&operation.fields.may_mutate);
    let auto_remap = bio_block_mask(&operation.fields.auto_remap);
    let must_handle = bio_state_mask(&operation.fields.must_handle);
    let needs_update = bio_derived_state_mask(&operation.fields.needs_update);
    Ok(quote! {
        #(#cfg_attrs)*
        pub(crate) struct #marker;
        #(#cfg_attrs)*
        impl #marker {
            pub(crate) const __COSMOLKIT_BODY_CLASS: usize = #class;
            pub(crate) const __COSMOLKIT_MAY_MUTATE: u64 = #may_mutate;
            pub(crate) const __COSMOLKIT_AUTO_REMAP: u64 = #auto_remap;
            pub(crate) const __COSMOLKIT_MUST_HANDLE: u64 = #must_handle;
            pub(crate) const __COSMOLKIT_NEEDS_UPDATE: u64 = #needs_update;
        }
    })
}

pub(crate) fn expand_molecule_access_markers(
    registry: &MoleculeRegistry,
) -> syn::Result<proc_macro2::TokenStream> {
    validate_unique_markers(registry.operations.iter().map(|operation| &operation.name))?;
    let markers = registry
        .operations
        .iter()
        .map(expand_molecule_marker)
        .collect::<syn::Result<Vec<_>>>()?;
    Ok(quote!(#(#markers)*))
}

pub(crate) fn expand_bio_access_markers(
    registry: &BioRegistry,
) -> syn::Result<proc_macro2::TokenStream> {
    validate_unique_markers(registry.operations.iter().map(|operation| &operation.name))?;
    let markers = registry
        .operations
        .iter()
        .map(expand_bio_marker)
        .collect::<syn::Result<Vec<_>>>()?;
    Ok(quote!(#(#markers)*))
}
