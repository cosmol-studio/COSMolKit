#[path = "../src/declaration.rs"]
mod declaration;
#[path = "../src/matrices.rs"]
mod matrices;

use declaration::{BioRegistry, MoleculeRegistry};
use matrices::{expand_bio_matrices, expand_molecule_matrices};
use proc_macro2::TokenStream;

fn compact(tokens: &TokenStream) -> String {
    tokens
        .to_string()
        .chars()
        .filter(|character| !character.is_whitespace())
        .collect()
}

fn molecule(source: &str) -> String {
    let registry = syn::parse_str::<MoleculeRegistry>(source)
        .unwrap_or_else(|error| panic!("expected Molecule registry: {error}"));
    compact(
        &expand_molecule_matrices(&registry)
            .unwrap_or_else(|error| panic!("expected Molecule matrices: {error}")),
    )
}

fn bio(source: &str) -> String {
    let registry = syn::parse_str::<BioRegistry>(source)
        .unwrap_or_else(|error| panic!("expected Bio registry: {error}"));
    compact(
        &expand_bio_matrices(&registry)
            .unwrap_or_else(|error| panic!("expected Bio matrices: {error}")),
    )
}

fn molecule_operation(name: &str, extra: &str) -> String {
    format!(
        r#"
        op {name} {{
            method: {name},
            impl_fn: crate::{name}_impl,
            kind: weak,
            access: {{ read: [], write: [] }},
            derived_effects: {{
                recompute: [], preserve: [], invalidate: [], operation_defined: [],
            }},
            cip_state: preserve,
            feature: crate::FEATURE,
            parity: not_applicable,
            invariant_profile: "empty",
            {extra}
        }}
        "#
    )
}

fn bio_operation(name: &str, domain: &str, extra: &str) -> String {
    format!(
        r#"
        op {name} {{
            method: {name},
            impl_fn: crate::{name}_impl,
            domain: {domain},
            kind: weak,
            edit_kind: local,
            may_mutate: [],
            auto_remap: [],
            must_handle: [],
            needs_update: [],
            feature: crate::BIO_FEATURE,
            parity: not_applicable,
            invariant_profile: "empty",
            {extra}
        }}
        "#
    )
}

#[test]
fn empty_registries_emit_exactly_the_eight_typed_tables() {
    let molecule = molecule("");
    for table in [
        "MOLECULE_OPS",
        "SUPPORT_MATRIX",
        "OPERATION_INVARIANT_MATRIX",
        "PARITY_MATRIX",
    ] {
        assert_eq!(molecule.matches(table).count(), 1, "{table}: {molecule}");
    }

    let bio = bio("");
    for table in [
        "BIO_STRUCTURE_OPS",
        "BIO_SUPPORT_MATRIX",
        "BIO_OPERATION_INVARIANT_MATRIX",
        "BIO_PARITY_MATRIX",
    ] {
        assert_eq!(bio.matches(table).count(), 1, "{table}: {bio}");
    }
}

#[test]
fn molecule_rows_preserve_every_block_effect_precondition_and_profile() {
    let output = molecule(
        r#"
        #[cfg(feature = "matrix-a")]
        op inspect {
            method: inspect,
            impl_fn: crate::inspect_impl,
            domain: coordinate,
            kind: weak,
            topology_edit: local,
            access: {
                read: [topology, properties],
                write: [coordinates, derived_cache],
            },
            may_mutate: [coordinates, derived_cache],
            auto_remap: [coordinates],
            derived_effects: {
                recompute: [rings, ring_families],
                preserve: [valence, aromaticity],
                invalidate: [stereo, coordinates, drawing, fingerprint],
                operation_defined: [],
            },
            cip_state: recompute,
            semantic_preconditions: [trusted_bond_topology, hydrogen_ownership_represented],
            requires_mapping: identity,
            feature: crate::INSPECT_FEATURE,
            parity: required_now,
            io_roundtrip: true,
            invariant_profile: "inspect-invariants",
            parity_profile: "inspect-parity",
        }
        "#,
    );

    for expected in [
        "pubconstINSPECT_SPEC:crate::ops::MoleculeOpSpec",
        "method:\"inspect\"",
        "impl_fn:\"crate::inspect_impl\"",
        "output:crate::ops::MoleculeOpOutput::Single",
        "result_type:\"Molecule\"",
        "domain:crate::ops::OperationDomain::Coordinate",
        "kind:crate::ops::MoleculeOpKind::Weak",
        "topology_edit:crate::ops::TopologyEditKind::Local",
        "crate::ops::BlockSet::TOPOLOGY",
        "crate::ops::BlockSet::COORDINATES",
        "crate::ops::BlockSet::PROPERTIES",
        "crate::ops::BlockSet::DERIVED_CACHE",
        "crate::DerivedState::RINGS",
        "crate::DerivedState::RING_FAMILIES",
        "crate::DerivedState::VALENCE",
        "crate::DerivedState::AROMATICITY",
        "crate::DerivedState::STEREO",
        "crate::DerivedState::COORDINATES",
        "crate::DerivedState::DRAWING",
        "crate::DerivedState::FINGERPRINT",
        "crate::ops::CipStatePolicy::Assign",
        "crate::ops::SemanticPreconditionSet::TRUSTED_BOND_TOPOLOGY",
        "crate::ops::SemanticPreconditionSet::HYDROGEN_OWNERSHIP_REPRESENTED",
        "requires_mapping:crate::ops::MappingRequirement::Identity",
        "support:crate::INSPECT_FEATURE.status",
        "parity:crate::ops::ParityPolicy::RequiredNow",
        "io_roundtrip:true",
    ] {
        assert!(output.contains(expected), "missing {expected} in {output}");
    }
    let invariant_start = output
        .find("pubconstOPERATION_INVARIANT_MATRIX")
        .expect("Molecule invariant matrix must be generated");
    let parity_start = output
        .find("pubconstPARITY_MATRIX")
        .expect("Molecule parity matrix must be generated");
    let invariant_matrix = &output[invariant_start..parity_start];
    assert_eq!(
        invariant_matrix
            .matches(
                "crate::ops::OperationInvariantEntry::for_operation(&INSPECT_SPEC,\"inspect-invariants\")",
            )
            .count(),
        1,
        "invariant row must bind INSPECT_SPEC and its invariant profile exactly once: {invariant_matrix}"
    );
    assert!(
        !invariant_matrix.contains("inspect-parity"),
        "parity profile must not leak into invariant row: {invariant_matrix}"
    );
    let parity_matrix = &output[parity_start..];
    assert!(
        parity_matrix.contains("profile:\"inspect-parity\""),
        "parity row must retain its distinct profile: {parity_matrix}"
    );
    assert_eq!(output.matches("#[cfg(feature=\"matrix-a\")]").count(), 5);
    assert_eq!(output.matches("&INSPECT_SPEC").count(), 4);
    assert_eq!(output.matches("&crate::INSPECT_FEATURE").count(), 2);
    assert!(!output.contains("implcrate::Molecule"));
    assert!(!output.contains("fninspect"));
}

#[test]
fn molecule_output_result_and_topology_edit_branches_are_exact() {
    let single_typed = molecule(&molecule_operation("typed", "result_type: crate::Report,"));
    assert!(
        single_typed.contains("result_type:stringify!((crate::Molecule,crate::Report))"),
        "{single_typed}"
    );

    let multiple = molecule(&molecule_operation("many", "output: multiple,"));
    assert!(multiple.contains("output:crate::ops::MoleculeOpOutput::Multiple"));
    assert!(multiple.contains("result_type:\"Vec<Molecule>\""));

    let typed_multiple = molecule(&molecule_operation(
        "many_typed",
        "output: multiple, result_type: crate::Report, assemble_fn: crate::assemble,",
    ));
    assert!(typed_multiple.contains("result_type:stringify!(crate::Report)"));

    for (name, edit, generated) in [
        ("compact", "compacting", "Compacting"),
        ("expand", "expanding", "Appending"),
        ("reorder", "reordering", "Renumbering"),
    ] {
        let source = format!(
            r#"
            op {name} {{
                method: {name}, impl_fn: crate::{name}_impl,
                kind: strong, topology_edit: {edit},
                access: {{ read: [], write: [topology] }},
                may_mutate: [topology],
                derived_effects: {{
                    recompute: [], preserve: [], invalidate: [], operation_defined: [],
                }},
                cip_state: clear, requires_mapping: required,
                feature: crate::FEATURE, parity: required_when_supported,
                invariant_profile: "strong", parity_profile: "strong-parity",
            }}
            "#
        );
        let output = molecule(&source);
        assert!(
            output.contains(&format!("TopologyEditKind::{generated}")),
            "{output}"
        );
        assert!(output.contains("CipStatePolicy::ClearComputed"));
        assert!(output.contains("MappingRequirement::Required"));
        assert!(output.contains("ParityPolicy::RequiredWhenSupported"));
    }
}

#[test]
fn molecule_tautomer_transition_and_operation_defined_valence_are_not_rewritten() {
    let tautomer = molecule(
        r#"
        op enumerate_tautomers_with_options(options: crate::Options) {
            method: enumerate_tautomers_with_options,
            impl_fn: crate::enumerate_tautomers_impl,
            output: multiple,
            kind: weak,
            access: { read: [], write: [topology, properties, derived_cache] },
            may_mutate: [topology, properties, derived_cache],
            derived_effects: {
                recompute: [], preserve: [], invalidate: [], operation_defined: [],
            },
            cip_state: tautomer_source_transition,
            feature: crate::TAUTOMER_FEATURE,
            parity: not_applicable,
            invariant_profile: "tautomer",
        }
        "#,
    );
    assert!(tautomer.contains("CipStatePolicy::TautomerSourceTransition"));
    assert!(!tautomer.contains("ParityMatrixEntry{"));

    let hydrogens = molecule(
        r#"
        op without_hydrogens {
            method: without_hydrogens,
            impl_fn: crate::without_hydrogens_impl,
            kind: weak,
            access: { read: [], write: [derived_cache] },
            may_mutate: [derived_cache],
            derived_effects: {
                recompute: [], preserve: [], invalidate: [], operation_defined: [valence],
            },
            cip_state: preserve,
            feature: crate::HYDROGENS_FEATURE,
            parity: not_applicable,
            invariant_profile: "hydrogens",
        }
        "#,
    );
    assert!(hydrogens.contains("crate::DerivedState::VALENCE"));
}

#[test]
fn molecule_parity_rows_follow_policy_without_fabricating_profiles() {
    let none = molecule(&molecule_operation("native", ""));
    assert_eq!(none.matches("ParityMatrixEntry{").count(), 0);

    for (name, parity, generated) in [
        (
            "planned",
            "required_when_supported",
            "RequiredWhenSupported",
        ),
        ("current", "required_now", "RequiredNow"),
    ] {
        let source = molecule_operation(
            name,
            &format!("parity_profile: \"{name}-profile\", parity: {parity},"),
        )
        .replacen("parity: not_applicable,", "", 1);
        let output = molecule(&source);
        assert_eq!(output.matches("ParityMatrixEntry{").count(), 1);
        assert!(output.contains(&format!("ParityPolicy::{generated}")));
        assert!(output.contains(&format!("profile:\"{name}-profile\"")));
    }
}

#[test]
fn bio_rows_preserve_all_blocks_states_derived_states_and_required_parity() {
    let output = bio(r#"
        #[cfg(feature = "bio-matrix")]
        op remove_waters {
            method: without_waters,
            impl_fn: crate::remove_waters_impl,
            domain: selection,
            kind: strong,
            edit_kind: compacting,
            may_mutate: [
                atoms, residues, chains, entities, models, coordinates, bonds,
                assemblies, annotations, derived_cache, properties,
            ],
            auto_remap: [coordinates, bonds, assemblies, annotations, properties],
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
            feature: crate::BIO_SELECTION_FEATURE,
            parity: required_now,
            io_roundtrip: true,
            invariant_profile: "bio-compact",
            parity_profile: "gemmi-remove-waters",
        }
        "#);

    for expected in [
        "pubconstBIO_REMOVE_WATERS_SPEC:crate::bio_ops::BioStructureOpSpec",
        "method:\"without_waters\"",
        "impl_fn:\"crate::remove_waters_impl\"",
        "domain:crate::bio_ops::BioOpDomain::Selection",
        "kind:crate::bio_ops::BioOpKind::Strong",
        "edit_kind:crate::bio_ops::BioEditKind::Compacting",
        "requires_mapping:crate::bio_ops::MappingRequirement::Required",
        "support:crate::BIO_SELECTION_FEATURE.status",
        "parity:crate::bio_ops::BioParityPolicy::RequiredNow",
        "profile:\"bio-compact\"",
        "profile:\"gemmi-remove-waters\"",
    ] {
        assert!(output.contains(expected), "missing {expected} in {output}");
    }
    for block in [
        "ATOMS",
        "RESIDUES",
        "CHAINS",
        "ENTITIES",
        "MODELS",
        "COORDINATES",
        "BONDS",
        "ASSEMBLIES",
        "ANNOTATIONS",
        "DERIVED_CACHE",
        "PROPERTIES",
    ] {
        assert!(output.contains(&format!("BioBlockSet::{block}")), "{block}");
    }
    for state in [
        "HIERARCHY",
        "RESIDUE_SPANS",
        "CHAIN_SPANS",
        "MODEL_SPANS",
        "COORDINATE_ALIGNMENT",
        "ENTITY_MAPPING",
        "ALTLOC_GROUPS",
        "ASSEMBLY_REFERENCES",
        "BOND_REFERENCES",
        "SELECTION_PROVENANCE",
        "POLYMER_ANNOTATION",
        "SECONDARY_STRUCTURE",
    ] {
        assert!(output.contains(&format!("BioStateSet::{state}")), "{state}");
    }
    for state in [
        "ATOM_INDEX",
        "RESIDUE_INDEX",
        "CHAIN_INDEX",
        "ENTITY_INDEX",
        "SEQUENCE_CACHE",
        "POLYMER_CACHE",
        "ALTLOC_CACHE",
        "ASSEMBLY_CACHE",
        "BOND_CACHE",
        "BACKBONE_GEOMETRY",
        "SIDECHAIN_GEOMETRY",
        "NUCLEIC_GEOMETRY",
        "SECONDARY_STRUCTURE",
        "CONTACT_MAP",
        "GRAPH_CACHE",
    ] {
        assert!(
            output.contains(&format!("BioDerivedState::{state}")),
            "{state}"
        );
    }
    assert_eq!(output.matches("#[cfg(feature=\"bio-matrix\")]").count(), 5);
    assert_eq!(output.matches("&BIO_REMOVE_WATERS_SPEC").count(), 4);
    assert!(!output.contains("implcrate::bio::BioStructure"));
}

#[test]
fn bio_domain_mapping_and_parity_omission_branches_are_exact() {
    for (index, (domain, generated)) in [
        ("selection", "Selection"),
        ("hierarchy", "Hierarchy"),
        ("coordinate", "Coordinate"),
        ("assembly", "Assembly"),
        ("annotation", "Annotation"),
        ("bonding", "Bonding"),
        ("polymer", "Polymer"),
        ("chemistry_bridge", "ChemistryBridge"),
    ]
    .into_iter()
    .enumerate()
    {
        let output = bio(&bio_operation(&format!("domain_{index}"), domain, ""));
        assert!(output.contains(&format!("BioOpDomain::{generated}")));
        assert!(output.contains("MappingRequirement::None"));
    }

    for parity in [
        "not_applicable",
        "gemmi_when_applicable",
        "biopython_when_applicable",
        "pdb_spec_required",
    ] {
        let source = bio_operation("parity_case", "hierarchy", &format!("parity: {parity},"))
            .replacen("parity: not_applicable,", "", 1);
        let output = bio(&source);
        assert_eq!(
            output.matches("BioParityMatrixEntry{").count(),
            0,
            "{output}"
        );
    }
}

#[test]
fn bio_strong_renumbering_and_identity_mapping_are_distinct() {
    let output = bio(r#"
        op renumber {
            method: renumber, impl_fn: crate::renumber_impl,
            domain: hierarchy, kind: strong, edit_kind: renumbering,
            may_mutate: [atoms], auto_remap: [],
            must_handle: [hierarchy], needs_update: [atom_index],
            requires_mapping: required,
            feature: crate::BIO_FEATURE,
            parity: not_applicable,
            invariant_profile: "renumber",
        }
        "#);
    assert!(output.contains("BioEditKind::Renumbering"));
    assert!(output.contains("MappingRequirement::Required"));

    let identity = bio(&bio_operation(
        "identity_edit",
        "hierarchy",
        "requires_mapping: identity,",
    ));
    assert!(identity.contains("MappingRequirement::Identity"));
}

#[test]
fn operation_cfg_is_explicit_and_non_cfg_attributes_are_rejected() {
    let output = molecule(&format!(
        "#[cfg(any(feature = \"a\", feature = \"b\"))] {}",
        molecule_operation("conditional", "")
    ));
    assert_eq!(
        output
            .matches("#[cfg(any(feature=\"a\",feature=\"b\"))]")
            .count(),
        4
    );

    let error = syn::parse_str::<MoleculeRegistry>(&format!(
        "#[allow(dead_code)] {}",
        molecule_operation("invalid_attribute", "")
    ))
    .err()
    .expect("non-cfg operation attribute must fail")
    .to_string();
    assert_eq!(
        error,
        "operation declarations accept only outer #[cfg(...)] attributes"
    );
}

#[test]
fn operation_order_and_feature_identity_are_shared_across_all_rows() {
    let source = format!(
        "{}{}",
        molecule_operation("first", ""),
        molecule_operation("second", "")
    );
    let output = molecule(&source);
    for table in [
        "MOLECULE_OPS",
        "SUPPORT_MATRIX",
        "OPERATION_INVARIANT_MATRIX",
    ] {
        let table = output.split(table).nth(1).expect("table suffix");
        assert!(table.find("FIRST_SPEC").unwrap() < table.find("SECOND_SPEC").unwrap());
    }
    assert_eq!(output.matches("feature:&crate::FEATURE").count(), 2);
    assert_eq!(output.matches("support:crate::FEATURE.status").count(), 2);
}
