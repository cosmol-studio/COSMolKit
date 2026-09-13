//! Private single-source operation specification and matrix generation.
//!
//! Generated tokens name runtime-owned values at the call site. This module
//! owns no runtime state and performs no operation or domain behavior.

use quote::{format_ident, quote};

use crate::declaration::{
    AccessFields, BioBlock, BioDerivedState, BioDomain, BioEditKind, BioOperation, BioParity,
    BioRegistry, BioState, CipStatePolicy, DerivedState, MappingRequirement, MoleculeBlock,
    MoleculeDomain, MoleculeOperation, MoleculeOutput, MoleculeParity, MoleculeRegistry,
    OperationKind, SemanticPrecondition, TopologyEditKind,
};

pub(crate) fn expand_molecule_matrices(
    registry: &MoleculeRegistry,
) -> syn::Result<proc_macro2::TokenStream> {
    let rows = registry
        .operations
        .iter()
        .map(expand_molecule_operation)
        .collect::<syn::Result<Vec<_>>>()?;
    let specs = rows.iter().map(|row| &row.spec);
    let operation_rows = rows.iter().map(|row| &row.operation_row);
    let support_rows = rows.iter().map(|row| &row.support_row);
    let invariant_rows = rows.iter().map(|row| &row.invariant_row);
    let parity_rows = rows.iter().filter_map(|row| row.parity_row.as_ref());

    Ok(quote! {
        #(#specs)*

        pub const MOLECULE_OPS: &[&crate::ops::MoleculeOpSpec] = &[
            #(#operation_rows,)*
        ];
        pub const SUPPORT_MATRIX: &[crate::ops::SupportMatrixEntry] = &[
            #(#support_rows,)*
        ];
        pub const OPERATION_INVARIANT_MATRIX: &[crate::ops::OperationInvariantEntry] = &[
            #(#invariant_rows,)*
        ];
        pub const PARITY_MATRIX: &[crate::ops::ParityMatrixEntry] = &[
            #(#parity_rows,)*
        ];
    })
}

pub(crate) fn expand_bio_matrices(registry: &BioRegistry) -> syn::Result<proc_macro2::TokenStream> {
    let rows = registry
        .operations
        .iter()
        .map(expand_bio_operation)
        .collect::<syn::Result<Vec<_>>>()?;
    let specs = rows.iter().map(|row| &row.spec);
    let operation_rows = rows.iter().map(|row| &row.operation_row);
    let support_rows = rows.iter().map(|row| &row.support_row);
    let invariant_rows = rows.iter().map(|row| &row.invariant_row);
    let parity_rows = rows.iter().filter_map(|row| row.parity_row.as_ref());

    Ok(quote! {
        #(#specs)*

        pub const BIO_STRUCTURE_OPS: &[&crate::bio_ops::BioStructureOpSpec] = &[
            #(#operation_rows,)*
        ];
        pub const BIO_SUPPORT_MATRIX: &[crate::bio_ops::BioSupportMatrixEntry] = &[
            #(#support_rows,)*
        ];
        pub const BIO_OPERATION_INVARIANT_MATRIX: &[crate::bio_ops::BioOperationInvariantEntry] = &[
            #(#invariant_rows,)*
        ];
        pub const BIO_PARITY_MATRIX: &[crate::bio_ops::BioParityMatrixEntry] = &[
            #(#parity_rows,)*
        ];
    })
}

struct MatrixRows {
    spec: proc_macro2::TokenStream,
    operation_row: proc_macro2::TokenStream,
    support_row: proc_macro2::TokenStream,
    invariant_row: proc_macro2::TokenStream,
    parity_row: Option<proc_macro2::TokenStream>,
}

fn expand_molecule_operation(operation: &MoleculeOperation) -> syn::Result<MatrixRows> {
    let cfg = &operation.cfg_attrs;
    let name = &operation.name;
    let fields = &operation.fields;
    let spec_ident = format_ident!("{}_SPEC", name.to_string().to_ascii_uppercase());
    let method = fields.method.to_string();
    let impl_fn = fields.impl_fn.to_token_stream().to_string();
    let output = molecule_output(fields.output);
    let result_type = molecule_result_type(operation);
    let domain = molecule_domain(fields.domain);
    let kind = operation_kind(fields.kind);
    let topology_edit = topology_edit(fields.topology_edit);
    let access = molecule_access(&fields.access);
    let may_mutate = molecule_blocks(&fields.may_mutate);
    let auto_remap = molecule_blocks(&fields.auto_remap);
    let recompute = derived_states(&fields.derived_effects.recompute);
    let preserve = derived_states(&fields.derived_effects.preserve);
    let invalidate = derived_states(&fields.derived_effects.invalidate);
    let operation_defined = derived_states(&fields.derived_effects.operation_defined);
    let cip_state = cip_state(fields.cip_state);
    let preconditions = semantic_preconditions(&fields.semantic_preconditions);
    let mapping = mapping_requirement(fields.requires_mapping, false);
    let feature = &fields.feature;
    let parity = molecule_parity(fields.parity);
    let io_roundtrip = fields.io_roundtrip;
    let invariant_profile = &fields.invariant_profile;

    let spec = quote! {
        #(#cfg)*
        pub const #spec_ident: crate::ops::MoleculeOpSpec = crate::ops::MoleculeOpSpec {
            method: #method,
            impl_fn: #impl_fn,
            output: #output,
            result_type: #result_type,
            domain: #domain,
            kind: #kind,
            topology_edit: #topology_edit,
            access: #access,
            may_mutate: #may_mutate,
            auto_remap: #auto_remap,
            derived_effects: crate::ops::DerivedEffects::new(
                #recompute,
                #preserve,
                #invalidate,
                #operation_defined,
            ),
            cip_state: #cip_state,
            semantic_preconditions: #preconditions,
            requires_mapping: #mapping,
            support: #feature.status,
            parity: #parity,
            io_roundtrip: #io_roundtrip,
        };
    };
    let operation_row = quote! { #(#cfg)* &#spec_ident };
    let support_row = quote! {
        #(#cfg)*
        crate::ops::SupportMatrixEntry {
            feature: &#feature,
            operation: Some(&#spec_ident),
        }
    };
    let invariant_row = quote! {
        #(#cfg)*
        crate::ops::OperationInvariantEntry::for_operation(&#spec_ident, #invariant_profile)
    };
    let parity_row = match fields.parity {
        MoleculeParity::NotApplicable => None,
        MoleculeParity::RequiredWhenSupported | MoleculeParity::RequiredNow => {
            let profile = fields.parity_profile.as_ref().ok_or_else(|| {
                syn::Error::new_spanned(name, "molecule parity matrix row requires parity_profile")
            })?;
            Some(quote! {
                #(#cfg)*
                crate::ops::ParityMatrixEntry {
                    operation: &#spec_ident,
                    feature: &#feature,
                    profile: #profile,
                    rdkit_version: None,
                }
            })
        }
    };

    Ok(MatrixRows {
        spec,
        operation_row,
        support_row,
        invariant_row,
        parity_row,
    })
}

fn expand_bio_operation(operation: &BioOperation) -> syn::Result<MatrixRows> {
    let cfg = &operation.cfg_attrs;
    let name = &operation.name;
    let fields = &operation.fields;
    let spec_ident = format_ident!("BIO_{}_SPEC", name.to_string().to_ascii_uppercase());
    let method = fields.method.to_string();
    let impl_fn = fields.impl_fn.to_token_stream().to_string();
    let domain = bio_domain(fields.domain);
    let kind = bio_kind(fields.kind);
    let edit_kind = bio_edit_kind(fields.edit_kind);
    let may_mutate = bio_blocks(&fields.may_mutate);
    let auto_remap = bio_blocks(&fields.auto_remap);
    let must_handle = bio_states(&fields.must_handle);
    let needs_update = bio_derived_states(&fields.needs_update);
    let mapping = mapping_requirement(fields.requires_mapping, true);
    let feature = &fields.feature;
    let parity = bio_parity(fields.parity);
    let io_roundtrip = fields.io_roundtrip;
    let invariant_profile = &fields.invariant_profile;

    let spec = quote! {
        #(#cfg)*
        pub const #spec_ident: crate::bio_ops::BioStructureOpSpec = crate::bio_ops::BioStructureOpSpec {
            method: #method,
            impl_fn: #impl_fn,
            domain: #domain,
            kind: #kind,
            edit_kind: #edit_kind,
            may_mutate: #may_mutate,
            auto_remap: #auto_remap,
            must_handle: #must_handle,
            needs_update: #needs_update,
            requires_mapping: #mapping,
            support: #feature.status,
            parity: #parity,
            io_roundtrip: #io_roundtrip,
        };
    };
    let operation_row = quote! { #(#cfg)* &#spec_ident };
    let support_row = quote! {
        #(#cfg)*
        crate::bio_ops::BioSupportMatrixEntry {
            feature: &#feature,
            operation: &#spec_ident,
        }
    };
    let invariant_row = quote! {
        #(#cfg)*
        crate::bio_ops::BioOperationInvariantEntry {
            operation: &#spec_ident,
            profile: #invariant_profile,
        }
    };
    let parity_row = if fields.parity == BioParity::RequiredNow {
        let profile = fields.parity_profile.as_ref().ok_or_else(|| {
            syn::Error::new_spanned(name, "Bio parity matrix row requires parity_profile")
        })?;
        Some(quote! {
            #(#cfg)*
            crate::bio_ops::BioParityMatrixEntry {
                operation: &#spec_ident,
                profile: #profile,
            }
        })
    } else {
        None
    };

    Ok(MatrixRows {
        spec,
        operation_row,
        support_row,
        invariant_row,
        parity_row,
    })
}

fn molecule_result_type(operation: &MoleculeOperation) -> proc_macro2::TokenStream {
    let fields = &operation.fields;
    match (
        fields.output,
        fields.result_type.as_ref(),
        fields.assemble_fn.as_ref(),
    ) {
        (MoleculeOutput::Single, None, _) => quote!("Molecule"),
        (MoleculeOutput::Single, Some(result), None) => {
            quote!(stringify!((crate::Molecule, #result)))
        }
        (MoleculeOutput::Single, Some(result), Some(_)) => quote!(stringify!(#result)),
        (MoleculeOutput::Multiple, None, _) => quote!("Vec<Molecule>"),
        (MoleculeOutput::Multiple, Some(result), _) => quote!(stringify!(#result)),
    }
}

fn molecule_output(value: MoleculeOutput) -> proc_macro2::TokenStream {
    match value {
        MoleculeOutput::Single => quote!(crate::ops::MoleculeOpOutput::Single),
        MoleculeOutput::Multiple => quote!(crate::ops::MoleculeOpOutput::Multiple),
    }
}

fn molecule_domain(value: MoleculeDomain) -> proc_macro2::TokenStream {
    match value {
        MoleculeDomain::Topology => quote!(crate::ops::OperationDomain::Topology),
        MoleculeDomain::Coordinate => quote!(crate::ops::OperationDomain::Coordinate),
    }
}

fn operation_kind(value: OperationKind) -> proc_macro2::TokenStream {
    match value {
        OperationKind::Weak => quote!(crate::ops::MoleculeOpKind::Weak),
        OperationKind::Strong => quote!(crate::ops::MoleculeOpKind::Strong),
    }
}

fn topology_edit(value: TopologyEditKind) -> proc_macro2::TokenStream {
    match value {
        TopologyEditKind::None => quote!(crate::ops::TopologyEditKind::None),
        TopologyEditKind::Local => quote!(crate::ops::TopologyEditKind::Local),
        TopologyEditKind::Compacting => quote!(crate::ops::TopologyEditKind::Compacting),
        TopologyEditKind::Expanding => quote!(crate::ops::TopologyEditKind::Appending),
        TopologyEditKind::Reordering => quote!(crate::ops::TopologyEditKind::Renumbering),
    }
}

fn molecule_access(access: &AccessFields) -> proc_macro2::TokenStream {
    let read = molecule_blocks(&access.read);
    let write = molecule_blocks(&access.write);
    quote!(crate::ops::BlockAccess::new(#read, #write))
}

fn molecule_blocks(values: &[MoleculeBlock]) -> proc_macro2::TokenStream {
    union(
        values.iter().map(|value| match value {
            MoleculeBlock::Topology => quote!(crate::ops::BlockSet::TOPOLOGY),
            MoleculeBlock::Coordinates => quote!(crate::ops::BlockSet::COORDINATES),
            MoleculeBlock::Properties => quote!(crate::ops::BlockSet::PROPERTIES),
            MoleculeBlock::DerivedCache => quote!(crate::ops::BlockSet::DERIVED_CACHE),
        }),
        quote!(crate::ops::BlockSet::NONE),
    )
}

fn derived_states(values: &[DerivedState]) -> proc_macro2::TokenStream {
    union(
        values.iter().map(|value| match value {
            DerivedState::Rings => quote!(crate::DerivedState::RINGS),
            DerivedState::RingFamilies => quote!(crate::DerivedState::RING_FAMILIES),
            DerivedState::Valence => quote!(crate::DerivedState::VALENCE),
            DerivedState::Aromaticity => quote!(crate::DerivedState::AROMATICITY),
            DerivedState::Stereo => quote!(crate::DerivedState::STEREO),
            DerivedState::Coordinates => quote!(crate::DerivedState::COORDINATES),
            DerivedState::Drawing => quote!(crate::DerivedState::DRAWING),
            DerivedState::Fingerprint => quote!(crate::DerivedState::FINGERPRINT),
        }),
        quote!(crate::DerivedState::NONE),
    )
}

fn semantic_preconditions(values: &[SemanticPrecondition]) -> proc_macro2::TokenStream {
    union(
        values.iter().map(|value| match value {
            SemanticPrecondition::TrustedBondTopology => {
                quote!(crate::ops::SemanticPreconditionSet::TRUSTED_BOND_TOPOLOGY)
            }
            SemanticPrecondition::HydrogenOwnershipRepresented => {
                quote!(crate::ops::SemanticPreconditionSet::HYDROGEN_OWNERSHIP_REPRESENTED)
            }
        }),
        quote!(crate::ops::SemanticPreconditionSet::NONE),
    )
}

fn cip_state(value: CipStatePolicy) -> proc_macro2::TokenStream {
    match value {
        CipStatePolicy::Preserve => quote!(crate::ops::CipStatePolicy::Preserve),
        CipStatePolicy::Clear => quote!(crate::ops::CipStatePolicy::ClearComputed),
        CipStatePolicy::Recompute => quote!(crate::ops::CipStatePolicy::Assign),
        CipStatePolicy::TautomerSourceTransition => {
            quote!(crate::ops::CipStatePolicy::TautomerSourceTransition)
        }
    }
}

fn molecule_parity(value: MoleculeParity) -> proc_macro2::TokenStream {
    match value {
        MoleculeParity::NotApplicable => quote!(crate::ops::ParityPolicy::NotApplicable),
        MoleculeParity::RequiredWhenSupported => {
            quote!(crate::ops::ParityPolicy::RequiredWhenSupported)
        }
        MoleculeParity::RequiredNow => quote!(crate::ops::ParityPolicy::RequiredNow),
    }
}

fn bio_domain(value: BioDomain) -> proc_macro2::TokenStream {
    match value {
        BioDomain::Selection => quote!(crate::bio_ops::BioOpDomain::Selection),
        BioDomain::Hierarchy => quote!(crate::bio_ops::BioOpDomain::Hierarchy),
        BioDomain::Coordinate => quote!(crate::bio_ops::BioOpDomain::Coordinate),
        BioDomain::Assembly => quote!(crate::bio_ops::BioOpDomain::Assembly),
        BioDomain::Annotation => quote!(crate::bio_ops::BioOpDomain::Annotation),
        BioDomain::Bonding => quote!(crate::bio_ops::BioOpDomain::Bonding),
        BioDomain::Polymer => quote!(crate::bio_ops::BioOpDomain::Polymer),
        BioDomain::ChemistryBridge => quote!(crate::bio_ops::BioOpDomain::ChemistryBridge),
    }
}

fn bio_kind(value: OperationKind) -> proc_macro2::TokenStream {
    match value {
        OperationKind::Weak => quote!(crate::bio_ops::BioOpKind::Weak),
        OperationKind::Strong => quote!(crate::bio_ops::BioOpKind::Strong),
    }
}

fn bio_edit_kind(value: BioEditKind) -> proc_macro2::TokenStream {
    match value {
        BioEditKind::None => quote!(crate::bio_ops::BioEditKind::None),
        BioEditKind::Local => quote!(crate::bio_ops::BioEditKind::Local),
        BioEditKind::Compacting => quote!(crate::bio_ops::BioEditKind::Compacting),
        BioEditKind::Expanding => quote!(crate::bio_ops::BioEditKind::Expanding),
        BioEditKind::Renumbering => quote!(crate::bio_ops::BioEditKind::Renumbering),
        BioEditKind::Splitting => quote!(crate::bio_ops::BioEditKind::Splitting),
        BioEditKind::Merging => quote!(crate::bio_ops::BioEditKind::Merging),
        BioEditKind::Transforming => quote!(crate::bio_ops::BioEditKind::Transforming),
    }
}

fn mapping_requirement(value: MappingRequirement, bio: bool) -> proc_macro2::TokenStream {
    let variant = match value {
        MappingRequirement::None => quote!(None),
        MappingRequirement::Identity => quote!(Identity),
        MappingRequirement::Required => quote!(Required),
    };
    if bio {
        quote!(crate::bio_ops::MappingRequirement::#variant)
    } else {
        quote!(crate::ops::MappingRequirement::#variant)
    }
}

fn bio_parity(value: BioParity) -> proc_macro2::TokenStream {
    match value {
        BioParity::NotApplicable => quote!(crate::bio_ops::BioParityPolicy::NotApplicable),
        BioParity::GemmiWhenApplicable => {
            quote!(crate::bio_ops::BioParityPolicy::GemmiWhenApplicable)
        }
        BioParity::BiopythonWhenApplicable => {
            quote!(crate::bio_ops::BioParityPolicy::BiopythonWhenApplicable)
        }
        BioParity::PdbSpecRequired => quote!(crate::bio_ops::BioParityPolicy::PdbSpecRequired),
        BioParity::RequiredNow => quote!(crate::bio_ops::BioParityPolicy::RequiredNow),
    }
}

fn bio_blocks(values: &[BioBlock]) -> proc_macro2::TokenStream {
    union(
        values.iter().map(|value| match value {
            BioBlock::Atoms => quote!(crate::bio_ops::BioBlockSet::ATOMS),
            BioBlock::Residues => quote!(crate::bio_ops::BioBlockSet::RESIDUES),
            BioBlock::Chains => quote!(crate::bio_ops::BioBlockSet::CHAINS),
            BioBlock::Entities => quote!(crate::bio_ops::BioBlockSet::ENTITIES),
            BioBlock::Models => quote!(crate::bio_ops::BioBlockSet::MODELS),
            BioBlock::Coordinates => quote!(crate::bio_ops::BioBlockSet::COORDINATES),
            BioBlock::Bonds => quote!(crate::bio_ops::BioBlockSet::BONDS),
            BioBlock::Assemblies => quote!(crate::bio_ops::BioBlockSet::ASSEMBLIES),
            BioBlock::Annotations => quote!(crate::bio_ops::BioBlockSet::ANNOTATIONS),
            BioBlock::DerivedCache => quote!(crate::bio_ops::BioBlockSet::DERIVED_CACHE),
            BioBlock::Properties => quote!(crate::bio_ops::BioBlockSet::PROPERTIES),
        }),
        quote!(crate::bio_ops::BioBlockSet::NONE),
    )
}

fn bio_states(values: &[BioState]) -> proc_macro2::TokenStream {
    union(
        values.iter().map(|value| match value {
            BioState::Hierarchy => quote!(crate::bio_ops::BioStateSet::HIERARCHY),
            BioState::ResidueSpans => quote!(crate::bio_ops::BioStateSet::RESIDUE_SPANS),
            BioState::ChainSpans => quote!(crate::bio_ops::BioStateSet::CHAIN_SPANS),
            BioState::ModelSpans => quote!(crate::bio_ops::BioStateSet::MODEL_SPANS),
            BioState::CoordinateAlignment => {
                quote!(crate::bio_ops::BioStateSet::COORDINATE_ALIGNMENT)
            }
            BioState::EntityMapping => quote!(crate::bio_ops::BioStateSet::ENTITY_MAPPING),
            BioState::AltlocGroups => quote!(crate::bio_ops::BioStateSet::ALTLOC_GROUPS),
            BioState::AssemblyReferences => {
                quote!(crate::bio_ops::BioStateSet::ASSEMBLY_REFERENCES)
            }
            BioState::BondReferences => quote!(crate::bio_ops::BioStateSet::BOND_REFERENCES),
            BioState::SelectionProvenance => {
                quote!(crate::bio_ops::BioStateSet::SELECTION_PROVENANCE)
            }
            BioState::PolymerAnnotation => quote!(crate::bio_ops::BioStateSet::POLYMER_ANNOTATION),
            BioState::SecondaryStructure => {
                quote!(crate::bio_ops::BioStateSet::SECONDARY_STRUCTURE)
            }
        }),
        quote!(crate::bio_ops::BioStateSet::NONE),
    )
}

fn bio_derived_states(values: &[BioDerivedState]) -> proc_macro2::TokenStream {
    union(
        values.iter().map(|value| match value {
            BioDerivedState::AtomIndex => quote!(crate::bio_ops::BioDerivedState::ATOM_INDEX),
            BioDerivedState::ResidueIndex => quote!(crate::bio_ops::BioDerivedState::RESIDUE_INDEX),
            BioDerivedState::ChainIndex => quote!(crate::bio_ops::BioDerivedState::CHAIN_INDEX),
            BioDerivedState::EntityIndex => quote!(crate::bio_ops::BioDerivedState::ENTITY_INDEX),
            BioDerivedState::SequenceCache => {
                quote!(crate::bio_ops::BioDerivedState::SEQUENCE_CACHE)
            }
            BioDerivedState::PolymerCache => quote!(crate::bio_ops::BioDerivedState::POLYMER_CACHE),
            BioDerivedState::AltlocCache => quote!(crate::bio_ops::BioDerivedState::ALTLOC_CACHE),
            BioDerivedState::AssemblyCache => {
                quote!(crate::bio_ops::BioDerivedState::ASSEMBLY_CACHE)
            }
            BioDerivedState::BondCache => quote!(crate::bio_ops::BioDerivedState::BOND_CACHE),
            BioDerivedState::BackboneGeometry => {
                quote!(crate::bio_ops::BioDerivedState::BACKBONE_GEOMETRY)
            }
            BioDerivedState::SidechainGeometry => {
                quote!(crate::bio_ops::BioDerivedState::SIDECHAIN_GEOMETRY)
            }
            BioDerivedState::NucleicGeometry => {
                quote!(crate::bio_ops::BioDerivedState::NUCLEIC_GEOMETRY)
            }
            BioDerivedState::SecondaryStructure => {
                quote!(crate::bio_ops::BioDerivedState::SECONDARY_STRUCTURE)
            }
            BioDerivedState::ContactMap => quote!(crate::bio_ops::BioDerivedState::CONTACT_MAP),
            BioDerivedState::GraphCache => quote!(crate::bio_ops::BioDerivedState::GRAPH_CACHE),
        }),
        quote!(crate::bio_ops::BioDerivedState::NONE),
    )
}

fn union(
    values: impl IntoIterator<Item = proc_macro2::TokenStream>,
    none: proc_macro2::TokenStream,
) -> proc_macro2::TokenStream {
    values.into_iter().fold(
        none,
        |accumulator, value| quote!(#accumulator.union(#value)),
    )
}

trait ToTokenString {
    fn to_token_stream(&self) -> proc_macro2::TokenStream;
}

impl ToTokenString for syn::Path {
    fn to_token_stream(&self) -> proc_macro2::TokenStream {
        quote!(#self)
    }
}
