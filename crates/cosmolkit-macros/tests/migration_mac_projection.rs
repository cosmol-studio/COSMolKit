extern crate proc_macro2 as proc_macro;

#[path = "../src/declaration.rs"]
mod declaration;
#[path = "../src/projection.rs"]
mod projection;

use std::fs;
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

use declaration::{BioRegistry, MoleculeRegistry};
use proc_macro2::TokenStream;
use projection::expand_molecule_access_markers;
use projection::{BodyClass, access_marker, expand_bio_access_markers, expand_body_tokens};
use quote::quote;
use syn::{Item, ItemFn, parse_quote};

fn compact(tokens: &TokenStream) -> String {
    tokens
        .to_string()
        .chars()
        .filter(|ch| !ch.is_whitespace())
        .collect()
}

fn expand(attribute: TokenStream, item: TokenStream, class: BodyClass) -> TokenStream {
    expand_body_tokens(attribute, item, class)
        .unwrap_or_else(|error| panic!("expected body expansion to succeed: {error}"))
}

fn expanded_function(tokens: &TokenStream) -> ItemFn {
    let file = syn::parse_file(&tokens.to_string())
        .unwrap_or_else(|error| panic!("expected expanded tokens to parse: {error}"));
    file.items
        .into_iter()
        .find_map(|item| match item {
            Item::Fn(function) => Some(function),
            _ => None,
        })
        .expect("expanded output must contain the operation body")
}

fn molecule_operation(name: &str, method: &str, output: &str, read: &str, write: &str) -> String {
    format!(
        r#"
        op {name} {{
            method: {method},
            impl_fn: crate::{name}_impl,
            output: {output},
            kind: weak,
            access: {{ read: [{read}], write: [{write}] }},
            derived_effects: {{
                recompute: [], preserve: [], invalidate: [], operation_defined: [],
            }},
            cip_state: preserve,
            feature: crate::FEATURE,
            parity: not_applicable,
            invariant_profile: "projection",
        }}
        "#
    )
}

fn parse_molecule(source: &str) -> MoleculeRegistry {
    syn::parse_str(source)
        .unwrap_or_else(|error| panic!("expected molecule declaration to parse: {error}"))
}

fn bio_operation(name: &str, method: &str, may_mutate: &str, auto_remap: &str) -> String {
    format!(
        r#"
        op {name} {{
            method: {method},
            impl_fn: crate::{name}_impl,
            domain: selection,
            kind: strong,
            edit_kind: compacting,
            may_mutate: [{may_mutate}],
            auto_remap: [{auto_remap}],
            must_handle: [
                hierarchy, residue_spans, chain_spans, model_spans,
                coordinate_alignment, entity_mapping, altloc_groups,
                assembly_references, bond_references, selection_provenance,
                polymer_annotation, secondary_structure,
            ],
            needs_update: [
                atom_index, residue_index, chain_index, entity_index,
                sequence_cache, polymer_cache, altloc_cache, assembly_cache,
                bond_cache, backbone_geometry, sidechain_geometry,
                nucleic_geometry, secondary_structure, contact_map, graph_cache,
            ],
            requires_mapping: required,
            feature: crate::BIO_FEATURE,
            parity: not_applicable,
            invariant_profile: "bio-projection",
        }}
        "#
    )
}

fn parse_bio(source: &str) -> BioRegistry {
    syn::parse_str(source)
        .unwrap_or_else(|error| panic!("expected Bio declaration to parse: {error}"))
}

static RUSTC_CASE: AtomicU64 = AtomicU64::new(0);

fn rustc(source: &str) -> Output {
    let case = RUSTC_CASE.fetch_add(1, Ordering::Relaxed);
    let directory = std::env::temp_dir().join(format!(
        "cosmolkit_mac_projection_{}_{}",
        std::process::id(),
        case
    ));
    fs::create_dir(&directory).expect("create isolated rustc fixture directory");
    let source_path = directory.join("fixture.rs");
    fs::write(&source_path, source).expect("write isolated rustc fixture");
    let compiler = std::env::var_os("RUSTC").unwrap_or_else(|| "rustc".into());
    let output = Command::new(compiler)
        .args(["--edition=2024", "--crate-type=lib"])
        .arg(&source_path)
        .arg("--out-dir")
        .arg(&directory)
        .output()
        .expect("run rustc for projection fixture");
    fs::remove_dir_all(&directory).expect("remove isolated rustc fixture directory");
    output
}

fn rustc_diagnostic(source: &str) -> String {
    let output = rustc(source);
    assert!(!output.status.success(), "fixture unexpectedly compiled");
    String::from_utf8_lossy(&output.stderr).into_owned()
}

#[test]
fn each_body_class_injects_its_exact_private_context_first() {
    let cases = [
        (
            BodyClass::MoleculeSingle,
            quote!(inspect, parts),
            quote!(
                fn body(value: usize) {}
            ),
            "parts:&mutcrate::OpParts<'_,crate::InspectAccess>",
            "[();0usize]",
        ),
        (
            BodyClass::MoleculeMultiple,
            quote!(enumerate, outputs),
            quote!(
                fn body(value: usize) {}
            ),
            "outputs:&mutcrate::MultiOutputOpParts<'_,crate::EnumerateAccess>",
            "[();1usize]",
        ),
        (
            BodyClass::Bio,
            quote!(remove_waters, bio),
            quote!(
                fn body(value: usize) {}
            ),
            "bio:&mutcrate::BioOpParts<'_,crate::RemoveWatersAccess>",
            "[();2usize]",
        ),
    ];

    for (class, attribute, item, context, assertion) in cases {
        let output = compact(&expand(attribute, item, class));
        assert!(output.contains(context), "missing context in {output}");
        assert!(
            output.contains(assertion),
            "missing class assertion in {output}"
        );
        assert!(output.contains("value:usize"));
        assert!(output.find(context).unwrap() < output.find("value:usize").unwrap());
    }
}

#[test]
fn zero_one_and_multiple_original_parameters_keep_order_and_types() {
    let zero = expanded_function(&expand(
        quote!(inspect, parts),
        quote!(
            fn zero() {}
        ),
        BodyClass::MoleculeSingle,
    ));
    assert_eq!(zero.sig.inputs.len(), 1);

    let one = expanded_function(&expand(
        quote!(inspect, parts),
        quote!(
            fn one(value: Option<&'static str>) {}
        ),
        BodyClass::MoleculeSingle,
    ));
    assert_eq!(one.sig.inputs.len(), 2);
    let one_input = one.sig.inputs.iter().nth(1).expect("original parameter");
    assert_eq!(compact(&quote!(#one_input)), "value:Option<&'staticstr>");

    let multiple = expanded_function(&expand(
        quote!(inspect, parts),
        quote!(
            fn many<'a, T>(first: &'a [T], second: (usize, bool), third: Option<T>) {}
        ),
        BodyClass::MoleculeSingle,
    ));
    let inputs = multiple.sig.inputs.iter().skip(1).collect::<Vec<_>>();
    assert_eq!(
        compact(&quote!(#(#inputs),*)),
        "first:&'a[T],second:(usize,bool),third:Option<T>"
    );
}

#[test]
fn function_shape_attributes_qualifiers_generics_where_return_and_body_are_preserved() {
    let original: ItemFn = parse_quote! {
        #[inline(always)]
        pub(crate) const unsafe extern "C" fn shaped<'a, T: Copy>(value: &'a T) -> usize
        where
            T: Clone,
        {
            let _ = value;
            7
        }
    };
    let expanded = expanded_function(&expand(
        quote!(inspect, parts),
        quote!(#original),
        BodyClass::MoleculeSingle,
    ));
    let expanded_attrs = &expanded.attrs;
    let original_attrs = &original.attrs;
    assert_eq!(
        compact(&quote!(#(#expanded_attrs)*)),
        compact(&quote!(#(#original_attrs)*))
    );
    let expanded_vis = &expanded.vis;
    let original_vis = &original.vis;
    assert_eq!(
        compact(&quote!(#expanded_vis)),
        compact(&quote!(#original_vis))
    );
    let expanded_constness = &expanded.sig.constness;
    let original_constness = &original.sig.constness;
    assert_eq!(
        compact(&quote!(#expanded_constness)),
        compact(&quote!(#original_constness))
    );
    let expanded_safety = &expanded.sig.safety;
    let original_safety = &original.sig.safety;
    assert_eq!(
        compact(&quote!(#expanded_safety)),
        compact(&quote!(#original_safety))
    );
    let expanded_abi = &expanded.sig.abi;
    let original_abi = &original.sig.abi;
    assert_eq!(
        compact(&quote!(#expanded_abi)),
        compact(&quote!(#original_abi))
    );
    assert_eq!(expanded.sig.ident, original.sig.ident);
    let expanded_generics = &expanded.sig.generics;
    let original_generics = &original.sig.generics;
    assert_eq!(
        compact(&quote!(#expanded_generics)),
        compact(&quote!(#original_generics))
    );
    let expanded_output = &expanded.sig.output;
    let original_output = &original.sig.output;
    assert_eq!(
        compact(&quote!(#expanded_output)),
        compact(&quote!(#original_output))
    );
    let expanded_block = &expanded.block;
    let original_block = &original.block;
    assert_eq!(
        compact(&quote!(#expanded_block)),
        compact(&quote!(#original_block))
    );

    let asynchronous = expanded_function(&expand(
        quote!(inspect, parts),
        quote!(
            pub async fn asynchronous(value: usize) -> usize {
                value
            }
        ),
        BodyClass::MoleculeSingle,
    ));
    assert!(asynchronous.sig.asyncness.is_some());
}

#[test]
fn marker_spelling_handles_components_and_rejects_raw_or_empty_components() {
    assert_eq!(
        access_marker(&parse_quote!(inspect)).unwrap(),
        "InspectAccess"
    );
    assert_eq!(
        access_marker(&parse_quote!(with_hydrogens)).unwrap(),
        "WithHydrogensAccess"
    );
    assert_eq!(
        access_marker(&parse_quote!(with_3d_conformer)).unwrap(),
        "With3dConformerAccess"
    );
    assert!(
        access_marker(&syn::parse_str("r#type").unwrap())
            .unwrap_err()
            .to_string()
            .contains("raw identifiers")
    );
    assert!(
        access_marker(&parse_quote!(foo__bar))
            .unwrap_err()
            .to_string()
            .contains("empty snake-case components")
    );
}

#[test]
fn declaration_rejects_distinct_names_that_generate_the_same_marker() {
    let source = format!(
        "{}{}",
        molecule_operation("foo_bar", "foo_bar", "single", "", ""),
        molecule_operation("foo_Bar", "foo_bar_upper", "single", "", "")
    );
    let registry = parse_molecule(&source);
    let error = expand_molecule_access_markers(&registry)
        .unwrap_err()
        .to_string();
    assert!(error.contains("FooBarAccess"));
    assert!(error.contains("foo_bar"));
    assert!(error.contains("foo_Bar"));
}

#[test]
fn molecule_projection_covers_every_block_in_read_and_write_without_inference() {
    let source = format!(
        "{}{}{}",
        molecule_operation(
            "read_all",
            "read_all",
            "single",
            "topology, coordinates, properties, derived_cache",
            "",
        ),
        molecule_operation(
            "write_all",
            "write_all",
            "multiple",
            "",
            "topology, coordinates, properties, derived_cache",
        ),
        molecule_operation("narrow", "narrow", "single", "topology", "coordinates"),
    );
    let tokens = compact(&expand_molecule_access_markers(&parse_molecule(&source)).unwrap());
    assert!(tokens.contains("ReadAllAccess"));
    assert!(tokens.contains("__COSMOLKIT_BODY_CLASS:usize=0usize"));
    assert!(tokens.contains("__COSMOLKIT_ACCESS_READ:u64=15u64"));
    assert!(tokens.contains("__COSMOLKIT_ACCESS_WRITE:u64=15u64"));
    assert!(tokens.contains("WriteAllAccess"));
    assert!(tokens.contains("__COSMOLKIT_BODY_CLASS:usize=1usize"));
    let mut narrow_impls = tokens.split("implNarrowAccess");
    let _before_narrow_impl = narrow_impls.next().expect("token prefix");
    let narrow = narrow_impls.next().expect("narrow marker impl");
    assert!(
        narrow_impls.next().is_none(),
        "narrow marker must have exactly one impl in {tokens}"
    );
    assert!(narrow.contains("__COSMOLKIT_ACCESS_READ:u64=1u64"));
    assert!(narrow.contains("__COSMOLKIT_ACCESS_WRITE:u64=2u64"));
}

#[test]
fn bio_projection_keeps_all_four_permission_vocabularies_distinct() {
    let all_blocks = "atoms, residues, chains, entities, models, coordinates, bonds, assemblies, annotations, derived_cache, properties";
    let tokens = compact(
        &expand_bio_access_markers(&parse_bio(&bio_operation(
            "remove_waters",
            "without_waters",
            all_blocks,
            all_blocks,
        )))
        .unwrap(),
    );
    assert!(tokens.contains("pub(crate)structRemoveWatersAccess;"));
    assert!(tokens.contains("__COSMOLKIT_BODY_CLASS:usize=2usize"));
    assert!(tokens.contains("__COSMOLKIT_MAY_MUTATE:u64=2047u64"));
    assert!(tokens.contains("__COSMOLKIT_AUTO_REMAP:u64=2047u64"));
    assert!(tokens.contains("__COSMOLKIT_MUST_HANDLE:u64=4095u64"));
    assert!(tokens.contains("__COSMOLKIT_NEEDS_UPDATE:u64=32767u64"));
    assert!(!tokens.contains("ACCESS_READ"));
    assert!(!tokens.contains("ACCESS_WRITE"));
}

#[test]
fn malformed_attributes_fail_closed() {
    let malformed = [
        quote!(),
        quote!(inspect),
        quote!(inspect,),
        quote!(, parts),
        quote!(crate::inspect, parts),
        quote!("inspect", parts),
        quote!(inspect, parts, extra),
        quote!(inspect, , parts),
    ];
    for attribute in malformed {
        assert!(
            expand_body_tokens(
                attribute.clone(),
                quote!(
                    fn body() {}
                ),
                BodyClass::MoleculeSingle,
            )
            .is_err(),
            "malformed attribute unexpectedly succeeded: {attribute}"
        );
    }
}

#[test]
fn receiver_diagnostics_name_the_invoked_macro() {
    let cases = [
        (BodyClass::MoleculeSingle, "mol_op_body"),
        (BodyClass::MoleculeMultiple, "mol_multi_op_body"),
        (BodyClass::Bio, "bio_op_body"),
    ];
    for (class, macro_name) in cases {
        let error = expand_body_tokens(
            quote!(inspect, parts),
            quote!(
                fn body(&self) {}
            ),
            class,
        )
        .unwrap_err()
        .to_string();
        assert!(error.contains(macro_name));
        assert!(error.contains("must not receive self"));
    }
}

#[test]
fn direct_context_binding_collision_is_rejected_for_each_body_class() {
    for class in [
        BodyClass::MoleculeSingle,
        BodyClass::MoleculeMultiple,
        BodyClass::Bio,
    ] {
        let error = expand_body_tokens(
            quote!(inspect, parts),
            quote!(
                fn body(parts: usize) {}
            ),
            class,
        )
        .unwrap_err()
        .to_string();
        assert!(error.contains("context name `parts`"));
        assert!(error.contains("conflicts with an existing function parameter"));
    }
}

#[test]
fn unknown_operation_identity_is_a_compile_failure() {
    let body = expand(
        quote!(unknown_operation, parts),
        quote!(
            fn body() {}
        ),
        BodyClass::MoleculeSingle,
    );
    let diagnostic = rustc_diagnostic(&format!(
        "#![allow(dead_code)] struct OpParts<'a, T>(&'a mut T); {body}"
    ));
    assert!(diagnostic.contains("UnknownOperationAccess"));
    assert!(diagnostic.contains("cannot find type") || diagnostic.contains("cannot find"));
}

#[test]
fn every_cross_class_mismatch_is_a_compile_failure() {
    let cases = [
        (
            BodyClass::MoleculeMultiple,
            "struct MultiOutputOpParts<'a, T>(&'a mut T);",
            0,
        ),
        (BodyClass::Bio, "struct BioOpParts<'a, T>(&'a mut T);", 0),
        (
            BodyClass::MoleculeSingle,
            "struct OpParts<'a, T>(&'a mut T);",
            2,
        ),
    ];
    for (class, context_type, declared_class) in cases {
        let body = expand(
            quote!(inspect, parts),
            quote!(
                fn body() {}
            ),
            class,
        );
        let diagnostic = rustc_diagnostic(&format!(
            "#![allow(dead_code)] {context_type} struct InspectAccess; impl InspectAccess {{ const __COSMOLKIT_BODY_CLASS: usize = {declared_class}; }} {body}"
        ));
        assert!(diagnostic.contains("mismatched types"));
    }
}

#[test]
fn generated_marker_is_private_fieldless_zero_sized_and_couples_to_a_valid_body() {
    let registry = parse_molecule(&molecule_operation("inspect", "inspect", "single", "", ""));
    let marker = expand_molecule_access_markers(&registry).unwrap();
    let body = expand(
        quote!(inspect, parts),
        quote!(
            fn body(value: usize) {
                let _ = value;
            }
        ),
        BodyClass::MoleculeSingle,
    );
    let marker_text = compact(&marker);
    assert!(marker_text.contains("pub(crate)structInspectAccess;"));
    assert!(!marker_text.contains("pubstructInspectAccess"));
    let output = rustc(&format!(
        "#![allow(dead_code)] struct OperationError; struct PreservationProof; struct OpParts<'a, T>(&'a mut T); impl<'a, T> OpParts<'a, T> {{ fn apply_cip_policy_runtime(&mut self) -> Result<(), OperationError> {{ Ok(()) }} }} {marker} const _: [(); 0] = [(); std::mem::size_of::<InspectAccess>()]; {body}"
    ));
    assert!(
        output.status.success(),
        "valid declaration/body coupling failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn operation_callsite_fixture(body: TokenStream) -> String {
    let registry = parse_molecule(&molecule_operation(
        "inspect",
        "inspect",
        "single",
        "topology",
        "coordinates",
    ));
    let marker = expand_molecule_access_markers(&registry).unwrap();
    let body = expand(quote!(inspect, parts), body, BodyClass::MoleculeSingle);
    format!(
        r#"
        #![allow(dead_code)]
        mod cosmolkit_model {{
            #[derive(Clone, Default)] pub(crate) struct TopologyBlock;
            #[derive(Clone, Default)] pub(crate) struct CoordinateBlock;
            #[derive(Clone, Default)] pub(crate) struct MoleculeProperties;
            #[derive(Clone, Default)] pub(crate) struct TopologyMapping;
        }}
        mod molecule {{ #[derive(Clone, Default)] pub(crate) struct DerivedCacheBlock; }}
        #[derive(Debug)] struct OperationError;
        struct PreservationProof;
        struct DerivedState;
        struct TopologyEditKind;
        struct MultiOutputOpParts<'a, T>(&'a mut T);
        struct OpParts<'a, T>(std::marker::PhantomData<(&'a (), T)>);
        impl<'a, T> OpParts<'a, T> {{
            fn read_topology_runtime(&self) -> Result<&cosmolkit_model::TopologyBlock, OperationError> {{
                static VALUE: cosmolkit_model::TopologyBlock = cosmolkit_model::TopologyBlock;
                Ok(&VALUE)
            }}
            fn checkout_coordinates_runtime(&mut self) -> Result<cosmolkit_model::CoordinateBlock, OperationError> {{
                Ok(cosmolkit_model::CoordinateBlock)
            }}
            fn install_coordinates_runtime(&mut self, _: cosmolkit_model::CoordinateBlock) -> Result<(), OperationError> {{ Ok(()) }}
            fn apply_cip_policy_runtime(&mut self) -> Result<(), OperationError> {{ Ok(()) }}
        }}
        {marker}
        {body}
        "#
    )
}

#[test]
fn declared_operation_body_methods_compile_at_the_real_generated_callsite() {
    let output = rustc(&operation_callsite_fixture(quote! {
        fn inspect_impl() -> Result<(), OperationError> {
            let _ = parts.topology()?;
            let coordinates = parts.checkout_coordinates()?;
            parts.install_coordinates(coordinates)?;
            parts.apply_cip_policy()
        }
    }));
    assert!(
        output.status.success(),
        "declared operation capabilities failed to compile: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn undeclared_operation_body_method_is_a_compile_failure() {
    let diagnostic = rustc_diagnostic(&operation_callsite_fixture(quote! {
        fn inspect_impl() -> Result<(), OperationError> {
            let _ = parts.properties()?;
            Ok(())
        }
    }));
    assert!(
        diagnostic.contains("no method named `properties`"),
        "{diagnostic}"
    );
}

#[test]
fn projection_contains_only_declaration_derived_block_methods() {
    let registry = parse_molecule(&molecule_operation(
        "inspect",
        "inspect",
        "single",
        "topology",
        "coordinates",
    ));
    let tokens = compact(&expand_molecule_access_markers(&registry).unwrap());
    for required in [
        "fntopology",
        "fncheckout_coordinates",
        "fninstall_coordinates",
        "fnapply_cip_policy",
    ] {
        assert!(
            tokens.contains(required),
            "missing generated method {required}"
        );
    }
    for forbidden in ["fncoordinates", "fnproperties", "fncheckout_topology"] {
        assert!(
            !tokens.contains(forbidden),
            "generated undeclared method {forbidden}"
        );
    }
}
