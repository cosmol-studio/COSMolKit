#[path = "../src/declaration.rs"]
mod declaration;

use declaration::{
    BioBlock, BioDerivedState, BioDomain, BioEditKind, BioParity, BioRegistry, BioState,
    CipStatePolicy, DerivedState, MappingRequirement, MoleculeBlock, MoleculeDomain,
    MoleculeOutput, MoleculeParity, MoleculeRegistry, OperationKind, SemanticPrecondition,
    TopologyEditKind,
};

fn molecule_source() -> String {
    r#"
        op inspect(value: usize, enabled: bool) {
            method: inspect_with_params,
            docs: "inspect a detached molecule",
            impl_fn: crate::inspect_impl,
            output: single,
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
            invariant_profile: "coordinate-write",
            parity_profile: "inspect-rdkit",
            default_method: inspect,
            default_args: [0usize, false],
            inplace: true,
            inplace_method: inspect_with_params_,
            inplace_docs: "mutate explicitly",
            default_inplace_method: inspect_,
        }
    "#
    .to_owned()
}

fn bio_source() -> String {
    r#"
        op remove_waters(selector: &str) {
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
            invariant_profile: "bio-compacting",
            parity_profile: "gemmi-remove-waters",
        }
    "#
    .to_owned()
}

fn parse_molecule(source: &str) -> MoleculeRegistry {
    match syn::parse_str(source) {
        Ok(value) => value,
        Err(error) => panic!("expected molecule declaration to parse: {error}"),
    }
}

fn parse_bio(source: &str) -> BioRegistry {
    match syn::parse_str(source) {
        Ok(value) => value,
        Err(error) => panic!("expected Bio declaration to parse: {error}"),
    }
}

fn molecule_error(source: &str) -> String {
    match syn::parse_str::<MoleculeRegistry>(source) {
        Ok(_) => panic!("expected molecule declaration to fail"),
        Err(error) => error.to_string(),
    }
}

fn bio_error(source: &str) -> String {
    match syn::parse_str::<BioRegistry>(source) {
        Ok(_) => panic!("expected Bio declaration to fail"),
        Err(error) => error.to_string(),
    }
}

fn replace(source: &str, from: &str, to: &str) -> String {
    assert!(
        source.contains(from),
        "fixture is missing replacement: {from}"
    );
    source.replacen(from, to, 1)
}

fn without_inplace_fixture(source: &str) -> String {
    let source = replace(source, "inplace: true", "inplace: false");
    let source = replace(
        &source,
        "            inplace_method: inspect_with_params_,\n",
        "",
    );
    let source = replace(
        &source,
        "            inplace_docs: \"mutate explicitly\",\n",
        "",
    );
    replace(
        &source,
        "            default_inplace_method: inspect_,\n",
        "",
    )
}

#[test]
fn molecule_full_declaration_preserves_types_order_and_every_set_value() {
    let registry = parse_molecule(&molecule_source());
    let operation = &registry.operations[0];
    let fields = &operation.fields;

    assert_eq!(operation.name.to_string(), "inspect");
    assert_eq!(operation.params.len(), 2);
    assert_eq!(fields.output, MoleculeOutput::Single);
    assert_eq!(fields.domain, MoleculeDomain::Coordinate);
    assert_eq!(fields.kind, OperationKind::Weak);
    assert_eq!(fields.topology_edit, TopologyEditKind::Local);
    assert_eq!(
        fields.access.read,
        [MoleculeBlock::Topology, MoleculeBlock::Properties]
    );
    assert_eq!(
        fields.access.write,
        [MoleculeBlock::Coordinates, MoleculeBlock::DerivedCache]
    );
    assert_eq!(
        fields.derived_effects.recompute,
        [DerivedState::Rings, DerivedState::RingFamilies]
    );
    assert_eq!(
        fields.derived_effects.preserve,
        [DerivedState::Valence, DerivedState::Aromaticity]
    );
    assert_eq!(
        fields.derived_effects.invalidate,
        [
            DerivedState::Stereo,
            DerivedState::Coordinates,
            DerivedState::Drawing,
            DerivedState::Fingerprint,
        ]
    );
    assert_eq!(
        fields.semantic_preconditions,
        [
            SemanticPrecondition::TrustedBondTopology,
            SemanticPrecondition::HydrogenOwnershipRepresented,
        ]
    );
    assert_eq!(fields.cip_state, CipStatePolicy::Recompute);
    assert_eq!(fields.requires_mapping, MappingRequirement::Identity);
    assert_eq!(fields.parity, MoleculeParity::RequiredNow);
    assert!(fields.io_roundtrip);
    assert!(fields.inplace);
    assert_eq!(
        fields.inplace_method.as_ref().unwrap(),
        "inspect_with_params_"
    );
    assert_eq!(fields.default_inplace_method.as_ref().unwrap(), "inspect_");
    assert_eq!(fields.default_args.len(), 2);
}

#[test]
fn molecule_defaults_and_optional_op_introducer_are_supported() {
    let source = r#"
        inspect {
            method: inspect,
            impl_fn: crate::inspect_impl,
            kind: weak,
            access: { read: [], write: [] },
            derived_effects: {
                recompute: [], preserve: [], invalidate: [], operation_defined: [],
            },
            cip_state: preserve,
            feature: crate::INSPECT_FEATURE,
            parity: not_applicable,
            invariant_profile: "read-only",
        }
    "#;
    let registry = parse_molecule(source);
    let fields = &registry.operations[0].fields;
    assert_eq!(fields.output, MoleculeOutput::Single);
    assert_eq!(fields.domain, MoleculeDomain::Topology);
    assert_eq!(fields.topology_edit, TopologyEditKind::None);
    assert_eq!(fields.requires_mapping, MappingRequirement::None);
    assert!(!fields.io_roundtrip);
    assert!(!fields.inplace);
}

#[test]
fn molecule_strong_edit_and_mapping_branches_parse() {
    for edit in ["compacting", "expanding", "reordering"] {
        let source = format!(
            r#"
                op reshape {{
                    method: reshape,
                    impl_fn: crate::reshape_impl,
                    kind: strong,
                    topology_edit: {edit},
                    access: {{ read: [], write: [topology] }},
                    may_mutate: [topology],
                    derived_effects: {{
                        recompute: [], preserve: [], invalidate: [], operation_defined: [],
                    }},
                    cip_state: clear,
                    requires_mapping: required,
                    feature: crate::RESHAPE_FEATURE,
                    parity: required_when_supported,
                    invariant_profile: "strong",
                    parity_profile: "reshape-rdkit",
                }}
            "#
        );
        let fields = &parse_molecule(&source).operations[0].fields;
        assert_eq!(fields.kind, OperationKind::Strong);
        assert_eq!(fields.requires_mapping, MappingRequirement::Required);
    }
}

#[test]
fn molecule_multiple_typed_result_and_tautomer_transition_parse() {
    let source = r#"
        op enumerate_tautomers_with_options(options: crate::Options) {
            method: enumerate_tautomers_with_options,
            impl_fn: crate::enumerate_impl,
            output: multiple,
            result_type: crate::TautomerResult,
            assemble_fn: crate::assemble_tautomers,
            kind: weak,
            access: { read: [coordinates], write: [topology, properties] },
            may_mutate: [topology, properties],
            derived_effects: {
                recompute: [], preserve: [], invalidate: [], operation_defined: [],
            },
            cip_state: tautomer_source_transition,
            feature: crate::TAUTOMER_FEATURE,
            parity: required_when_supported,
            invariant_profile: "multiple",
            parity_profile: "tautomer-rdkit",
        }
    "#;
    let fields = &parse_molecule(source).operations[0].fields;
    assert_eq!(fields.output, MoleculeOutput::Multiple);
    assert!(fields.result_type.is_some());
    assert!(fields.assemble_fn.is_some());
    assert_eq!(fields.cip_state, CipStatePolicy::TautomerSourceTransition);
}

#[test]
fn molecule_hydrogen_operation_defined_allowlist_is_exact() {
    let source = r#"
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
            feature: crate::HYDROGEN_FEATURE,
            parity: required_when_supported,
            invariant_profile: "hydrogen",
            parity_profile: "hydrogen-rdkit",
        }
    "#;
    assert_eq!(
        parse_molecule(source).operations[0]
            .fields
            .derived_effects
            .operation_defined,
        [DerivedState::Valence]
    );
    let bad = replace(source, "op without_hydrogens", "op sanitize");
    assert!(molecule_error(&bad).contains("explicit human-author approval"));
}

#[test]
fn molecule_missing_unknown_and_removed_fields_have_focused_diagnostics() {
    let missing = replace(
        &molecule_source(),
        "            method: inspect_with_params,\n",
        "",
    );
    assert!(molecule_error(&missing).contains("missing 'method'"));
    let unknown = replace(
        &molecule_source(),
        "            docs: \"inspect a detached molecule\",",
        "            mystery: true,",
    );
    assert!(molecule_error(&unknown).contains("unknown molecule_ops field 'mystery'"));

    for (field, diagnostic) in [
        ("must_handle", "removed from molecule_ops"),
        ("needs_update", "not a molecule declaration field"),
        ("invalidates_old", "derived_effects.invalidate"),
        ("invalidates", "derived_effects.invalidate"),
        ("rdkit_parity", "structured parity policy"),
        ("report", "discarded legacy field"),
    ] {
        let source = replace(
            &molecule_source(),
            "            docs: \"inspect a detached molecule\",",
            &format!("            {field}: [],"),
        );
        assert!(
            molecule_error(&source).contains(diagnostic),
            "field {field}"
        );
    }
}

#[test]
fn molecule_rejects_duplicate_fields_list_members_and_access_overlap() {
    let duplicate_field = replace(
        &molecule_source(),
        "            method: inspect_with_params,",
        "            method: inspect_with_params,\n            method: inspect_again,",
    );
    assert!(molecule_error(&duplicate_field).contains("duplicate molecule_ops field 'method'"));

    let duplicate_list = replace(
        &molecule_source(),
        "read: [topology, properties]",
        "read: [topology, topology]",
    );
    assert!(molecule_error(&duplicate_list).contains("duplicate 'topology' in 'access.read'"));

    let overlap = replace(
        &molecule_source(),
        "write: [coordinates, derived_cache]",
        "write: [topology, coordinates, derived_cache]",
    );
    assert!(molecule_error(&overlap).contains("cannot appear in access.read and access.write"));
}

#[test]
fn molecule_rejects_duplicate_nested_fields_and_derived_state_overlap() {
    let duplicate_nested = replace(
        &molecule_source(),
        "read: [topology, properties],",
        "read: [topology, properties],\n                read: [],",
    );
    assert!(molecule_error(&duplicate_nested).contains("duplicate access field 'read'"));

    let overlap = replace(
        &molecule_source(),
        "preserve: [valence, aromaticity]",
        "preserve: [rings, valence]",
    );
    assert!(molecule_error(&overlap).contains("appears in more than one effect category"));
}

#[test]
fn molecule_rejects_removed_derived_effect_fields() {
    for (field, diagnostic) in [
        ("requires", "derived_effects.requires was removed"),
        ("unsupported", "structured errors"),
        ("require_handle", "removed with derived_effects.requires"),
    ] {
        let source = replace(
            &molecule_source(),
            "operation_defined: [],",
            &format!("{field}: [],"),
        );
        assert!(
            molecule_error(&source).contains(diagnostic),
            "field {field}"
        );
    }
}

#[test]
fn molecule_rejects_unknown_typed_vocabulary() {
    for (from, to, diagnostic) in [
        (
            "output: single",
            "output: stream",
            "unknown molecule output",
        ),
        (
            "domain: coordinate",
            "domain: drawing",
            "unknown operation domain",
        ),
        ("kind: weak", "kind: maybe", "unknown operation kind"),
        (
            "topology_edit: local",
            "topology_edit: shuffle",
            "unknown topology edit",
        ),
        (
            "coordinates, derived_cache",
            "unknown_block, derived_cache",
            "unknown molecule block",
        ),
        (
            "recompute: [rings, ring_families]",
            "recompute: [mystery]",
            "unknown derived state",
        ),
        (
            "cip_state: recompute",
            "cip_state: retain",
            "unknown CIP state policy",
        ),
        (
            "requires_mapping: identity",
            "requires_mapping: optional",
            "unknown mapping requirement",
        ),
        (
            "parity: required_now",
            "parity: perhaps",
            "unknown molecule parity policy",
        ),
        (
            "trusted_bond_topology, hydrogen_ownership_represented",
            "unknown_precondition",
            "unknown semantic precondition",
        ),
    ] {
        assert!(molecule_error(&replace(&molecule_source(), from, to)).contains(diagnostic));
    }
}

#[test]
fn molecule_permission_relationships_are_enforced() {
    let mutation = replace(
        &molecule_source(),
        "may_mutate: [coordinates, derived_cache]",
        "may_mutate: [topology]",
    );
    assert!(molecule_error(&mutation).contains("may_mutate block"));

    let remap = replace(
        &molecule_source(),
        "auto_remap: [coordinates]",
        "auto_remap: [properties]",
    );
    assert!(molecule_error(&remap).contains("auto_remap block"));

    let cache = replace(
        &molecule_source(),
        "write: [coordinates, derived_cache]",
        "write: [coordinates]",
    );
    let cache = replace(
        &cache,
        "may_mutate: [coordinates, derived_cache]",
        "may_mutate: [coordinates]",
    );
    assert!(molecule_error(&cache).contains("derived effects require derived_cache write access"));
}

#[test]
fn molecule_edit_mapping_relationships_are_enforced() {
    let weak_index_change = replace(
        &molecule_source(),
        "topology_edit: local",
        "topology_edit: compacting",
    );
    assert!(molecule_error(&weak_index_change).contains("must be strong and require a mapping"));

    let strong_local = replace(&molecule_source(), "kind: weak", "kind: strong");
    assert!(molecule_error(&strong_local).contains("strong molecule operations"));
}

#[test]
fn molecule_output_result_and_assembler_relationships_are_enforced() {
    let multiple_inplace = replace(&molecule_source(), "output: single", "output: multiple");
    assert!(molecule_error(&multiple_inplace).contains("cannot generate an in-place wrapper"));

    let result_without_assembler = replace(
        &molecule_source(),
        "            output: single,",
        "            output: multiple,\n            result_type: crate::ResultType,",
    );
    let result_without_assembler = without_inplace_fixture(&result_without_assembler);
    assert!(molecule_error(&result_without_assembler).contains("requires assemble_fn"));

    let assembler_without_result = replace(
        &molecule_source(),
        "            output: single,",
        "            output: multiple,\n            assemble_fn: crate::assemble,",
    );
    let assembler_without_result = without_inplace_fixture(&assembler_without_result);
    assert!(molecule_error(&assembler_without_result).contains("requires result_type"));

    let single_assembler = replace(
        &molecule_source(),
        "            output: single,",
        "            output: single,\n            assemble_fn: crate::assemble,",
    );
    assert!(molecule_error(&single_assembler).contains("only for multiple-output"));
}

#[test]
fn molecule_default_and_inplace_relationships_are_enforced() {
    let no_inplace = replace(&molecule_source(), "inplace: true", "inplace: false");
    assert!(molecule_error(&no_inplace).contains("inplace_method requires inplace: true"));

    let bad_inplace_name = replace(
        &molecule_source(),
        "inplace_method: inspect_with_params_",
        "inplace_method: inspect_with_params",
    );
    assert!(molecule_error(&bad_inplace_name).contains("must end with '_'"));

    let bad_default_arity = replace(
        &molecule_source(),
        "default_args: [0usize, false]",
        "default_args: [0usize, false, true]",
    );
    assert!(molecule_error(&bad_default_arity).contains("has 3 entries"));

    let no_default = replace(&molecule_source(), "default_method: inspect,\n", "");
    assert!(molecule_error(&no_default).contains("default_args requires default_method"));
}

#[test]
fn molecule_method_suffix_and_parity_rules_are_enforced() {
    let value_suffix = replace(
        &molecule_source(),
        "method: inspect_with_params",
        "method: inspect_with_params_",
    );
    assert!(molecule_error(&value_suffix).contains("value-style 'method' must not end"));

    let no_profile = replace(
        &molecule_source(),
        "            parity_profile: \"inspect-rdkit\",\n",
        "",
    );
    assert!(molecule_error(&no_profile).contains("requires parity_profile"));

    let not_applicable = replace(
        &molecule_source(),
        "parity: required_now",
        "parity: not_applicable",
    );
    assert!(molecule_error(&not_applicable).contains("must not declare parity_profile"));
}

#[test]
fn molecule_tautomer_transition_guardrail_rejects_each_wrong_shape() {
    let good = r#"
        enumerate_tautomers_with_options {
            method: enumerate_tautomers_with_options,
            impl_fn: crate::enumerate,
            output: multiple,
            kind: weak,
            access: { read: [], write: [topology, properties] },
            may_mutate: [topology, properties],
            derived_effects: { recompute: [], preserve: [], invalidate: [], operation_defined: [] },
            cip_state: tautomer_source_transition,
            feature: crate::TAUTOMER_FEATURE,
            parity: not_applicable,
            invariant_profile: "tautomer",
        }
    "#;
    parse_molecule(good);
    let topology_only = replace(good, "write: [topology, properties]", "write: [topology]");
    let topology_only = replace(
        &topology_only,
        "may_mutate: [topology, properties]",
        "may_mutate: [topology]",
    );
    for bad in [
        replace(good, "enumerate_tautomers_with_options {", "tautomers {"),
        replace(
            good,
            "method: enumerate_tautomers_with_options",
            "method: tautomers",
        ),
        replace(good, "output: multiple", "output: single"),
        topology_only,
    ] {
        assert!(molecule_error(&bad).contains("tautomer_source_transition is permitted only"));
    }
}

#[test]
fn molecule_registry_and_parameter_identities_are_unique() {
    let one = r#"
        op first(value: usize) {
            method: first,
            impl_fn: crate::first,
            kind: weak,
            access: { read: [], write: [] },
            derived_effects: { recompute: [], preserve: [], invalidate: [], operation_defined: [] },
            cip_state: preserve,
            feature: crate::FEATURE,
            parity: not_applicable,
            invariant_profile: "read",
        }
    "#;
    assert!(
        molecule_error(&format!("{one}{one}")).contains("duplicate molecule operation 'first'")
    );
    let method_collision = format!("{one}{}", one.replace("op first", "op second"));
    assert!(
        molecule_error(&method_collision).contains("duplicate generated molecule method 'first'")
    );
    let duplicate_parameter = one.replace("value: usize", "value: usize, value: bool");
    assert!(molecule_error(&duplicate_parameter).contains("duplicate parameter 'value'"));
    let destructured = one.replace("value: usize", "(value, other): (usize, usize)");
    assert!(molecule_error(&destructured).contains("must use simple identifiers"));
}

#[test]
fn bio_full_declaration_preserves_every_typed_set_value() {
    let registry = parse_bio(&bio_source());
    let operation = &registry.operations[0];
    let fields = &operation.fields;
    assert_eq!(operation.params.len(), 1);
    assert_eq!(fields.domain, BioDomain::Selection);
    assert_eq!(fields.kind, OperationKind::Strong);
    assert_eq!(fields.edit_kind, BioEditKind::Compacting);
    assert_eq!(fields.may_mutate.len(), 11);
    assert_eq!(fields.auto_remap.len(), 5);
    assert_eq!(fields.must_handle.len(), 12);
    assert_eq!(fields.needs_update.len(), 15);
    assert!(fields.may_mutate.contains(&BioBlock::Atoms));
    assert!(fields.must_handle.contains(&BioState::SecondaryStructure));
    assert!(fields.needs_update.contains(&BioDerivedState::GraphCache));
    assert_eq!(fields.requires_mapping, MappingRequirement::Required);
    assert_eq!(fields.parity, BioParity::RequiredNow);
    assert!(fields.io_roundtrip);
}

#[test]
fn bio_defaults_optional_op_and_all_domains_parse() {
    let base = r#"
        inspect {
            method: inspect,
            impl_fn: crate::inspect,
            kind: weak,
            feature: crate::FEATURE,
            parity: not_applicable,
            invariant_profile: "read",
        }
    "#;
    let fields = &parse_bio(base).operations[0].fields;
    assert_eq!(fields.domain, BioDomain::Hierarchy);
    assert_eq!(fields.edit_kind, BioEditKind::None);
    assert_eq!(fields.requires_mapping, MappingRequirement::None);

    for domain in [
        "selection",
        "hierarchy",
        "coordinate",
        "assembly",
        "annotation",
        "bonding",
        "polymer",
        "chemistry_bridge",
    ] {
        parse_bio(&replace(
            base,
            "kind: weak,",
            &format!("domain: {domain},\n            kind: weak,"),
        ));
    }
}

#[test]
fn bio_edit_and_mapping_branches_are_exact() {
    let weak = r#"
        op edit {
            method: edit,
            impl_fn: crate::edit,
            kind: weak,
            edit_kind: local,
            feature: crate::FEATURE,
            parity: not_applicable,
            invariant_profile: "local",
        }
    "#;
    for edit in ["none", "local", "transforming"] {
        parse_bio(&replace(
            weak,
            "edit_kind: local",
            &format!("edit_kind: {edit}"),
        ));
    }
    for edit in ["compacting", "renumbering"] {
        let source = replace(weak, "kind: weak", "kind: strong");
        let source = replace(&source, "edit_kind: local", &format!("edit_kind: {edit}"));
        let source = replace(
            &source,
            "feature: crate::FEATURE,",
            "requires_mapping: required,\n            feature: crate::FEATURE,",
        );
        parse_bio(&source);
    }
    for edit in ["expanding", "splitting", "merging"] {
        let source = replace(weak, "edit_kind: local", &format!("edit_kind: {edit}"));
        assert!(
            bio_error(&source).contains("unresolved source-indexed/per-output mapping contract")
        );
    }
}

#[test]
fn bio_parity_branches_are_exact() {
    let base = r#"
        inspect {
            method: inspect,
            impl_fn: crate::inspect,
            kind: weak,
            feature: crate::FEATURE,
            parity: not_applicable,
            invariant_profile: "read",
        }
    "#;
    for parity in [
        "not_applicable",
        "gemmi_when_applicable",
        "biopython_when_applicable",
        "pdb_spec_required",
    ] {
        parse_bio(&replace(
            base,
            "parity: not_applicable",
            &format!("parity: {parity}"),
        ));
    }
    let required = replace(
        base,
        "parity: not_applicable,",
        "parity: required_now,\n            parity_profile: \"bio-current\",",
    );
    parse_bio(&required);
    let missing = replace(
        &required,
        "            parity_profile: \"bio-current\",\n",
        "",
    );
    assert!(bio_error(&missing).contains("required_now requires parity_profile"));
    let future_with_profile = replace(
        base,
        "parity: not_applicable,",
        "parity: gemmi_when_applicable,\n            parity_profile: \"not-current\",",
    );
    assert!(bio_error(&future_with_profile).contains("only parity: required_now"));
}

#[test]
fn bio_rejects_unknown_duplicate_and_inconsistent_fields() {
    let unknown = replace(&bio_source(), "domain: selection", "domain: mystery");
    assert!(bio_error(&unknown).contains("unknown Bio operation domain"));

    let duplicate = replace(
        &bio_source(),
        "method: without_waters,",
        "method: without_waters,\n            method: again,",
    );
    assert!(bio_error(&duplicate).contains("duplicate bio_structure_ops field 'method'"));

    let duplicate_member = replace(
        &bio_source(),
        "auto_remap: [coordinates, bonds, assemblies, annotations, properties]",
        "auto_remap: [coordinates, coordinates]",
    );
    assert!(bio_error(&duplicate_member).contains("duplicate 'coordinates' in 'auto_remap'"));

    let invalid_remap = replace(
        &bio_source(),
        "auto_remap: [coordinates, bonds, assemblies, annotations, properties]",
        "auto_remap: [coordinates, bonds, assemblies, annotations, properties]",
    );
    let invalid_remap = replace(
        &invalid_remap,
        "models, coordinates, bonds",
        "models, bonds",
    );
    assert!(bio_error(&invalid_remap).contains("Bio auto_remap block"));
}

#[test]
fn bio_rejects_unknown_blocks_states_and_policies() {
    for (from, to, diagnostic) in [
        ("atoms, residues", "mystery, residues", "unknown Bio block"),
        (
            "hierarchy, residue_spans",
            "mystery, residue_spans",
            "unknown Bio handled state",
        ),
        (
            "atom_index, residue_index",
            "mystery, residue_index",
            "unknown Bio derived state",
        ),
        (
            "edit_kind: compacting",
            "edit_kind: mystery",
            "unknown Bio edit kind",
        ),
        (
            "requires_mapping: required",
            "requires_mapping: mystery",
            "unknown mapping requirement",
        ),
        (
            "parity: required_now",
            "parity: mystery",
            "unknown Bio parity policy",
        ),
    ] {
        assert!(bio_error(&replace(&bio_source(), from, to)).contains(diagnostic));
    }
}

#[test]
fn bio_registry_and_parameter_identities_are_unique() {
    let one = r#"
        op first(value: usize) {
            method: first,
            impl_fn: crate::first,
            kind: weak,
            feature: crate::FEATURE,
            parity: not_applicable,
            invariant_profile: "read",
        }
    "#;
    assert!(bio_error(&format!("{one}{one}")).contains("duplicate Bio operation 'first'"));
    let method_collision = format!("{one}{}", one.replace("op first", "op second"));
    assert!(bio_error(&method_collision).contains("duplicate generated Bio method 'first'"));
    let duplicate_parameter = one.replace("value: usize", "value: usize, value: bool");
    assert!(bio_error(&duplicate_parameter).contains("duplicate parameter 'value'"));
}

#[test]
fn historical_declaration_shape_matrix_remains_parseable() {
    let mut source = String::new();
    for index in 0..24 {
        let (domain, cip) = if index % 2 == 0 {
            ("topology", "preserve")
        } else {
            ("coordinate", "clear")
        };
        source.push_str(&format!(
            r#"
                op historical_{index}(value: usize) {{
                    method: historical_{index},
                    impl_fn: crate::historical_{index},
                    domain: {domain},
                    kind: weak,
                    topology_edit: local,
                    access: {{ read: [topology], write: [] }},
                    derived_effects: {{
                        recompute: [], preserve: [], invalidate: [], operation_defined: [],
                    }},
                    cip_state: {cip},
                    feature: crate::HISTORICAL_FEATURE,
                    parity: not_applicable,
                    invariant_profile: "historical-shape",
                }}
            "#
        ));
    }
    assert_eq!(parse_molecule(&source).operations.len(), 24);
    parse_bio(&bio_source());
}

#[test]
fn current_cosmolkit_operations_have_disjoint_access_for_every_cfg_gate() {
    let file = syn::parse_file(include_str!("../../cosmolkit/src/ops/registry.rs"))
        .expect("current operation registry parses as Rust");
    let tokens = file
        .items
        .into_iter()
        .find_map(|item| match item {
            syn::Item::Macro(item) if item.mac.path.is_ident("molecule_ops") => {
                Some(item.mac.tokens)
            }
            _ => None,
        })
        .expect("one molecule_ops invocation");
    let registry: MoleculeRegistry =
        syn::parse2(tokens).expect("all current operation declarations satisfy access validation");

    for operation in &registry.operations {
        for block in &operation.fields.access.read {
            assert!(
                !operation.fields.access.write.contains(block),
                "{} exposes {block:?} through both read and write",
                operation.name
            );
        }
    }

    let potential = registry
        .operations
        .iter()
        .find(|operation| operation.name == "potential_stereo")
        .expect("potential_stereo remains registered");
    assert!(potential.fields.access.read.is_empty());
    assert_eq!(
        potential.fields.access.write,
        vec![MoleculeBlock::Topology, MoleculeBlock::DerivedCache]
    );
}
