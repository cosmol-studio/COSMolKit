use std::collections::HashSet;

use cosmolkit::{
    BINDING_CONTRACT, BindingCallableContract, BindingDefault, BindingExposure, BindingItem,
    BindingKind, BindingOwner, BindingParameterContract, BindingParity, BindingSupport,
    BindingTypeRole, StateModel,
};

fn entry(semantic_id: &str) -> &'static cosmolkit::BindingContractEntry {
    BINDING_CONTRACT
        .iter()
        .find(|entry| entry.semantic_id == semantic_id)
        .unwrap_or_else(|| panic!("missing binding contract entry {semantic_id}"))
}

#[cfg(feature = "io")]
#[test]
fn sdf_reader_contracts_distinguish_concrete_and_preserving_results() {
    for (id, output) in [
        ("Molecule.from_sdf", "crate::Molecule"),
        ("Molecule.from_sdf_with_params", "crate::Molecule"),
        ("SdfRecord.from_sdf", "crate::SdfRecord"),
        ("SdfRecord.from_sdf_with_params", "crate::SdfRecord"),
    ] {
        let contract = entry(id);
        assert_eq!(contract.feature, "io");
        assert_eq!(contract.exposure, BindingExposure::Public);
        assert_eq!(contract.support, BindingSupport::Supported);
        assert_eq!(contract.parity, BindingParity::RequiredWhenSupported);
        let callable = contract.callable.unwrap();
        assert_eq!(callable.output_type.replace(' ', ""), output);
        assert_eq!(
            callable.error_type.unwrap().replace(' ', ""),
            "crate::SdfError"
        );
        assert_eq!(callable.state_model, StateModel::ValueReturning);
        assert_eq!(callable.kind, BindingKind::Static);
    }
    for (id, role) in [
        ("types.SdfRecord", BindingTypeRole::Result),
        ("types.SdfGraph", BindingTypeRole::Result),
        ("types.SdfCoordinateMode", BindingTypeRole::Parameter),
        ("types.SdfReadParams", BindingTypeRole::Parameter),
        ("types.SdfError", BindingTypeRole::Error),
    ] {
        let contract = entry(id);
        assert_eq!(contract.type_role, Some(role));
        assert_eq!(contract.feature, "io");
        assert_eq!(contract.exposure, BindingExposure::Public);
        assert_eq!(contract.support, BindingSupport::Supported);
        assert_eq!(contract.parity, BindingParity::RequiredWhenSupported);
    }
    for (id, output, error) in [
        ("SdfRecord.graph", "&crate::SdfGraph", None),
        (
            "SdfRecord.molecule",
            "&crate::Molecule",
            Some("crate::SdfError"),
        ),
        (
            "SdfRecord.query_graph",
            "&crate::QueryGraph",
            Some("crate::SdfError"),
        ),
        ("SdfRecord.data_fields", "&[(String,String)]", None),
        ("SdfRecord.properties", "&crate::MoleculeProperties", None),
        (
            "SdfRecord.substance_groups",
            "&[crate::SubstanceGroup]",
            None,
        ),
        (
            "SdfRecord.source_coordinate_dim",
            "Option<crate::CoordinateDimension>",
            None,
        ),
    ] {
        let contract = entry(id);
        assert_eq!(contract.owner, BindingOwner::Type);
        assert_eq!(contract.exposure, BindingExposure::Public);
        assert_eq!(contract.support, BindingSupport::Supported);
        assert_eq!(contract.parity, BindingParity::RequiredWhenSupported);
        let callable = contract.callable.unwrap();
        assert_eq!(callable.kind, BindingKind::Instance);
        assert_eq!(callable.state_model, StateModel::ReadOnly);
        assert_eq!(callable.operation_semantic_id, None);
        assert_eq!(callable.output_type.replace(' ', ""), output);
        assert_eq!(
            callable.error_type.map(|name| name.replace(' ', "")),
            error.map(str::to_owned)
        );
    }
    for id in [
        "types.MolBlockReadParams",
        "types.MolBlockError",
        "Molecule.from_molblock",
        "Molecule.from_molblock_with_params",
    ] {
        assert!(
            BINDING_CONTRACT.iter().all(|entry| entry.semantic_id != id),
            "unimplemented interface must not appear in the public registry: {id}",
        );
    }
}

#[test]
fn canonical_registry_preserves_order_and_feature_local_subsets() {
    let mut expected = vec![
        "types.Molecule",
        "types.MoleculeBuilder",
        "Molecule.new",
        "Molecule.from_parts",
        "Molecule.to_builder",
        "Molecule.num_atoms",
        "Molecule.num_bonds",
        "Molecule.atoms",
        "Molecule.bonds",
        "Molecule.atom",
        "Molecule.bond",
        "Molecule.topology",
        "Molecule.coordinates_2d",
        "Molecule.conformers_3d",
        "Molecule.properties",
        "Molecule.property",
        "MoleculeBuilder.new",
        "MoleculeBuilder.from_parts",
        "MoleculeBuilder.build",
        "MoleculeBuilder.add_atom",
        "MoleculeBuilder.add_bond",
        "MoleculeBuilder.set_atom_formal_charge",
        "MoleculeBuilder.set_bond_order",
        "MoleculeBuilder.remove_bond_between_atoms",
        "MoleculeBuilder.degree",
        "MoleculeBuilder.neighbor_bonds",
        "MoleculeBuilder.bond_between_atoms",
        "MoleculeBuilder.atoms",
        "MoleculeBuilder.bonds",
        "MoleculeBuilder.substance_groups",
        "MoleculeBuilder.stereo_groups",
        "MoleculeBuilder.coordinates",
        "MoleculeBuilder.properties",
        "MoleculeBuilder.set_2d_coordinates",
        "MoleculeBuilder.add_2d_conformer",
        "MoleculeBuilder.add_3d_conformer",
        "MoleculeBuilder.add_substance_group",
        "MoleculeBuilder.add_stereo_group",
        "MoleculeBuilder.with_name",
        "MoleculeBuilder.with_property",
        "MoleculeBuilder.with_sdf_data_field",
        "MoleculeBuilder.with_properties",
        "types.SubstanceGroupId",
        "types.SubstanceGroupKind",
        "types.SubstanceGroup",
        "types.SGroupBracket",
        "types.SGroupCState",
        "types.SGroupDisplay",
        "SubstanceGroupId.new",
        "SubstanceGroup.new",
        "SubstanceGroupId.index",
        "SGroupBracket.new",
        "SGroupCState.new",
        "SGroupBracket.points",
        "SGroupCState.bond",
        "SGroupCState.vector",
        "SGroupDisplay.brackets",
        "SubstanceGroup.id",
        "SubstanceGroup.kind",
        "SubstanceGroup.display",
        "SubstanceGroup.cstates",
        "SubstanceGroup.head_crossing_bonds",
        "SubstanceGroup.crossing_bond_correspondence",
        "types.OperationError",
        "types.TemplateAttachment",
        "types.TemplateAttachmentOrder",
        "types.TemplateAttachmentOrderError",
        "TemplateAttachment.target",
        "TemplateAttachment.label",
        "TemplateAttachmentOrder.entries",
        "Atom.template_attachment_order",
        "types.SupportStatus",
        "types.ParityPolicy",
        "types.FeatureSpec",
        "types.MoleculeOpSpec",
        "types.SupportMatrixEntry",
        "types.OperationInvariantEntry",
        "types.ParityMatrixEntry",
        "types.FeatureSpecIter",
        "module.feature_specs",
        "module.feature_spec",
        "module.operation_specs",
        "module.operation_spec",
        "module.support_matrix",
        "module.operation_invariant_matrix",
        "module.parity_matrix",
        "module.operation_invariant",
        "module.operation_parity",
        "module.version",
    ];
    if cfg!(feature = "matrices") {
        expected.extend([
            "types.DenseMatrix",
            "types.DistanceMatrixParams",
            "types.DistanceMatrix3dParams",
            "types.MatrixError",
            "Molecule.distance_matrix",
            "Molecule.distance_matrix_with_params",
            "Molecule.distance_matrix_3d",
            "Molecule.distance_matrix_3d_with_params",
        ]);
    }
    if cfg!(feature = "depict") {
        expected.extend([
            "types.Coordinate2DParams",
            "types.Coordinate2DError",
            "types.Coordinate2DTemplateError",
            "types.Coordinate2DLayoutError",
            "Molecule.with_2d_coordinates",
            "Molecule.with_2d_coordinates_with_params",
        ]);
    }
    if cfg!(feature = "transforms") {
        expected.extend([
            "types.AtomPositionParams",
            "types.TransformError",
            "Molecule.with_atom_position",
            "Molecule.with_atom_position_with_params",
            "Molecule.set_atom_position_",
            "Molecule.set_atom_position_with_params_",
        ]);
    }
    if cfg!(feature = "sanitize") {
        expected.extend([
            "types.SanitizeOperations",
            "types.SanitizeStage",
            "types.SanitizeParams",
            "types.SanitizeError",
            "types.ChemistryProblemError",
            "types.ChemistryProblem",
            "types.ChemistryProblemReport",
            "Molecule.sanitize",
            "Molecule.sanitize_with_params",
            "Molecule.detect_chemistry_problems",
            "Molecule.detect_chemistry_problems_with_params",
        ]);
    }
    if cfg!(feature = "kekulize") {
        expected.extend([
            "types.KekulizeParams",
            "types.KekulizeError",
            "Molecule.with_kekulized_bonds",
            "Molecule.with_kekulized_bonds_with_params",
            "Molecule.kekulize_bonds_",
            "Molecule.kekulize_bonds_with_params_",
        ]);
    }
    if cfg!(feature = "aromaticity") {
        expected.extend([
            "types.AromaticityModel",
            "types.AromaticityParams",
            "types.AromaticityError",
            "Molecule.with_assigned_aromaticity",
            "Molecule.with_assigned_aromaticity_with_params",
            "Molecule.assign_aromaticity_",
            "Molecule.assign_aromaticity_with_params_",
        ]);
    }
    if cfg!(feature = "valence") {
        expected.extend([
            "types.ValenceModel",
            "types.ValenceParams",
            "types.ValenceError",
            "Molecule.with_assigned_valence",
            "Molecule.with_assigned_valence_with_params",
            "Molecule.assign_valence_",
            "Molecule.assign_valence_with_params_",
            "Molecule.has_valence_violation",
        ]);
    }
    if cfg!(feature = "radicals") {
        expected.extend([
            "Molecule.with_assigned_radicals",
            "Molecule.assign_radicals_",
        ]);
    }
    if cfg!(feature = "rings") {
        expected.extend([
            "types.RingSearchParams",
            "Molecule.with_assigned_rings",
            "Molecule.assign_rings_",
            "Molecule.with_assigned_ring_families",
            "Molecule.with_assigned_ring_families_with_params",
            "Molecule.assign_ring_families_",
            "Molecule.assign_ring_families_with_params_",
        ]);
    }
    if cfg!(feature = "stereo") {
        expected.extend([
            "types.StructureTagParams",
            "types.StereoError",
            "Molecule.with_chiral_tags_from_structure",
            "Molecule.with_chiral_tags_from_structure_with_params",
            "Molecule.assign_chiral_tags_from_structure_",
            "Molecule.assign_chiral_tags_from_structure_with_params_",
            "types.PotentialStereoParams",
            "types.PotentialStereoType",
            "types.PotentialStereoSpecified",
            "types.PotentialStereoDescriptor",
            "types.PotentialStereoCenter",
            "types.PotentialStereoInfo",
            "types.RingStereoRelation",
            "types.PotentialStereoResult",
            "types.PotentialStereoError",
            "Molecule.potential_stereo",
            "Molecule.potential_stereo_with_params",
            "types.CipDescriptor",
            "types.CipDescriptorError",
            "types.CipLabelOptions",
            "types.CipLabelerError",
            "Molecule.with_cip_labels",
            "Molecule.with_cip_labels_with_options",
            "Molecule.assign_cip_labels_",
            "Molecule.assign_cip_labels_with_options_",
        ]);
    }
    if cfg!(feature = "descriptors") {
        expected.extend([
            "Molecule.molecular_weight",
            "Molecule.exact_molecular_weight",
            "Molecule.molecular_formula",
        ]);
    }
    if cfg!(feature = "hydrogens") {
        expected.extend([
            "types.AddHsParams",
            "types.HydrogenError",
            "types.RemoveHsParams",
            "Molecule.with_hydrogens",
            "Molecule.with_hydrogens_with_params",
            "Molecule.add_hydrogens_",
            "Molecule.add_hydrogens_with_params_",
            "Molecule.without_hydrogens",
            "Molecule.without_hydrogens_with_params",
            "Molecule.remove_hydrogens_",
            "Molecule.remove_hydrogens_with_params_",
        ]);
    }
    if cfg!(feature = "bio") {
        expected.extend([
            "types.ResidueInfoKind",
            "types.ResidueCode",
            "types.ResidueInfo",
            "types.PdbAtomSerial",
            "types.PdbChainId",
            "types.PdbSeqId",
            "types.AtomName",
            "types.ResidueName",
            "types.AltLocLabel",
            "types.AtomSourceIds",
            "types.ResidueSourceIds",
            "types.ChainSourceIds",
            "types.EntitySourceIds",
            "types.ResidueSequenceError",
            "module.residue_info",
            "module.residue_info_checked",
            "module.find_residue_info_index",
            "module.find_residue_info",
            "module.residue_code",
            "module.expand_one_letter",
            "module.expand_one_letter_sequence",
            "types.BioStructure",
            "types.BioStructureParts",
            "types.BioStructureError",
            "types.BioCoordinateFormat",
            "types.BioCalcFlag",
            "types.EntityKind",
            "types.PolymerKind",
            "types.ResidueKind",
            "types.ChainKind",
            "types.BioRowSpan",
            "types.BioAtomId",
            "types.BioResidueId",
            "types.BioChainId",
            "types.BioEntityId",
            "types.BioModelId",
            "types.BioAssemblyId",
            "types.BioAltLocGroupId",
            "types.BioAtomRow",
            "types.BioResidueRow",
            "types.BioChainRow",
            "types.BioEntityRow",
            "types.BioModelRow",
            "types.BioCoordinateBlock",
            "types.BioTransform",
            "types.BioCrystalCell",
            "types.BioCrystalInfo",
            "types.BioNcsOperator",
            "types.BioAssemblyOperator",
            "types.BioAssemblyGenerator",
            "types.BioAssemblySpecialKind",
            "types.BioAssembly",
            "types.AltLocRequest",
            "types.BioEntityDbRef",
            "types.BioSiftsUnpResidue",
            "types.Protein",
            "types.ProteinProjectionError",
            "types.ProteinChainRef",
            "types.ProteinResidueRef",
            "types.ProteinAtomRef",
            "BioStructure.protein",
            "Protein.num_models",
            "Protein.num_chains",
            "Protein.num_residues",
            "Protein.num_atoms",
            "Protein.chains",
            "Protein.chain",
            "Protein.residues",
            "Protein.atoms",
            "ProteinChainRef.id",
            "ProteinChainRef.row",
            "ProteinChainRef.kind",
            "ProteinChainRef.source",
            "ProteinChainRef.residues",
            "ProteinChainRef.atoms",
            "ProteinResidueRef.id",
            "ProteinResidueRef.row",
            "ProteinResidueRef.name",
            "ProteinResidueRef.kind",
            "ProteinResidueRef.info",
            "ProteinResidueRef.code",
            "ProteinResidueRef.one_letter_code",
            "ProteinResidueRef.fasta_code",
            "ProteinResidueRef.is_standard",
            "ProteinResidueRef.chain",
            "ProteinResidueRef.atoms",
            "ProteinAtomRef.id",
            "ProteinAtomRef.row",
            "ProteinAtomRef.name",
            "ProteinAtomRef.element",
            "ProteinAtomRef.altloc",
            "ProteinAtomRef.residue",
            "ProteinAtomRef.position",
        ]);
    }
    if cfg!(feature = "io") {
        expected.extend([
            "types.SdfRecord",
            "types.SdfGraph",
            "SdfRecord.graph",
            "SdfRecord.molecule",
            "SdfRecord.query_graph",
            "SdfRecord.data_fields",
            "SdfRecord.properties",
            "SdfRecord.substance_groups",
            "SdfRecord.source_coordinate_dim",
            "SdfRecord.from_sdf",
            "SdfRecord.from_sdf_with_params",
            "types.SdfCoordinateMode",
            "types.SdfReadParams",
            "types.SdfError",
            "Molecule.from_sdf",
            "Molecule.from_sdf_with_params",
        ]);
    }
    if cfg!(feature = "smiles") {
        expected.extend([
            "types.SmilesParseParams",
            "types.SmilesError",
            "types.SmilesStereoError",
            "Molecule.from_smiles",
            "Molecule.from_smiles_with_params",
        ]);
    }

    if cfg!(feature = "fingerprints") {
        // Exact newly registered public sparse-count surface, in registry order.
        expected.splice(
            0..0,
            [
                "types.SparseCountFingerprint",
                "types.SparseCountFingerprint32",
                "types.FingerprintError",
                "SparseCountFingerprint.new",
                "SparseCountFingerprint.length",
                "SparseCountFingerprint.value",
                "SparseCountFingerprint.set_value",
                "SparseCountFingerprint.nonzero_elements",
                "SparseCountFingerprint.total_value",
                "SparseCountFingerprint.fuzzy_and",
                "SparseCountFingerprint.fuzzy_or",
                "SparseCountFingerprint.with_added",
                "SparseCountFingerprint.with_subtracted",
                "SparseCountFingerprint.with_added_scalar",
                "SparseCountFingerprint.with_subtracted_scalar",
                "SparseCountFingerprint.with_multiplied_scalar",
                "SparseCountFingerprint.with_divided_scalar",
                "SparseCountFingerprint32.new",
                "SparseCountFingerprint32.length",
                "SparseCountFingerprint32.value",
                "SparseCountFingerprint32.set_value",
                "SparseCountFingerprint32.nonzero_elements",
                "SparseCountFingerprint32.total_value",
                "SparseCountFingerprint32.fuzzy_and",
                "SparseCountFingerprint32.fuzzy_or",
                "SparseCountFingerprint32.with_added",
                "SparseCountFingerprint32.with_subtracted",
                "SparseCountFingerprint32.with_added_scalar",
                "SparseCountFingerprint32.with_subtracted_scalar",
                "SparseCountFingerprint32.with_multiplied_scalar",
                "SparseCountFingerprint32.with_divided_scalar",
            ],
        );
    }
    assert_eq!(
        BINDING_CONTRACT
            .iter()
            .map(|entry| entry.semantic_id)
            .collect::<Vec<_>>(),
        expected
    );
    assert_eq!(
        BINDING_CONTRACT
            .iter()
            .map(|entry| entry.semantic_id)
            .collect::<HashSet<_>>()
            .len(),
        BINDING_CONTRACT.len(),
        "semantic identities must be unique"
    );
}

#[test]
fn always_present_entries_have_exact_type_and_module_payloads() {
    let molecule = entry("types.Molecule");
    assert_eq!(molecule.item, BindingItem::Type);
    assert_eq!(molecule.owner, BindingOwner::Type);
    assert!(molecule.rust_path.ends_with("Molecule"));
    assert_eq!(molecule.python_name, "Molecule");
    assert_eq!(molecule.javascript_name, "Molecule");
    assert_eq!(molecule.feature, "runtime");
    assert_eq!(molecule.exposure, BindingExposure::Public);
    assert_eq!(molecule.support, BindingSupport::Supported);
    assert_eq!(molecule.parity, BindingParity::NotApplicable);
    assert_eq!(molecule.callable, None);
    assert_eq!(molecule.type_role, Some(BindingTypeRole::Value));

    let error = entry("types.OperationError");
    assert_eq!(error.item, BindingItem::Type);
    assert_eq!(error.owner, BindingOwner::Type);
    assert!(error.rust_path.ends_with("OperationError"));
    assert_eq!(error.type_role, Some(BindingTypeRole::Error));
    assert_eq!(error.callable, None);

    let version = entry("module.version");
    assert_eq!(version.item, BindingItem::Callable);
    assert_eq!(version.owner, BindingOwner::Module);
    assert_eq!(version.rust_path.replace(' ', ""), "crate::version");
    assert_eq!(version.python_name, "version");
    assert_eq!(version.javascript_name, "version");
    assert_eq!(version.feature, "metadata");
    assert_eq!(version.exposure, BindingExposure::Public);
    assert_eq!(version.support, BindingSupport::Supported);
    assert_eq!(version.parity, BindingParity::NotApplicable);
    assert_eq!(version.type_role, None);
    let callable = version.callable.expect("version callable metadata");
    assert_eq!(callable.kind, BindingKind::Module);
    assert!(callable.parameters.is_empty());
    assert_eq!(callable.output_type.replace(' ', ""), "&'staticstr");
    assert_eq!(callable.error_type, None);
    assert_eq!(callable.state_model, StateModel::ReadOnly);
    assert_eq!(callable.operation_semantic_id, None);
}

#[cfg(feature = "descriptors")]
#[test]
fn descriptor_entries_are_exact_canonical_read_only_methods() {
    for (semantic_id, rust_name, javascript, output) in [
        (
            "Molecule.molecular_weight",
            "molecular_weight",
            "molecularWeight",
            "f64",
        ),
        (
            "Molecule.exact_molecular_weight",
            "exact_molecular_weight",
            "exactMolecularWeight",
            "f64",
        ),
        (
            "Molecule.molecular_formula",
            "molecular_formula",
            "molecularFormula",
            "String",
        ),
    ] {
        let entry = entry(semantic_id);
        assert_eq!(entry.item, BindingItem::Callable);
        assert_eq!(entry.owner, BindingOwner::Molecule);
        assert!(entry.rust_path.replace(' ', "").ends_with(rust_name));
        assert_eq!(entry.python_name, rust_name);
        assert_eq!(entry.javascript_name, javascript);
        assert_eq!(entry.feature, "descriptors");
        assert_eq!(entry.exposure, BindingExposure::Public);
        assert_eq!(entry.support, BindingSupport::Experimental);
        assert_eq!(entry.parity, BindingParity::RequiredWhenSupported);
        assert_eq!(entry.type_role, None);
        let callable = entry.callable.expect("descriptor callable metadata");
        assert_eq!(callable.kind, BindingKind::Instance);
        assert!(callable.parameters.is_empty());
        assert_eq!(callable.output_type, output);
        assert_eq!(
            callable.error_type.map(|name| name.replace(' ', "")),
            Some("crate::OperationError".to_owned())
        );
        assert_eq!(callable.state_model, StateModel::ReadOnly);
        assert_eq!(callable.operation_semantic_id, None);
    }
}

#[cfg(feature = "hydrogens")]
#[test]
fn hydrogen_entries_bind_value_and_in_place_names_with_parity_support() {
    for (semantic_id, rust_name, javascript, output, state) in [
        (
            "Molecule.with_hydrogens",
            "with_hydrogens",
            "withHydrogens",
            "crate::Molecule",
            StateModel::ValueReturning,
        ),
        (
            "Molecule.add_hydrogens_",
            "add_hydrogens_",
            "addHydrogens",
            "()",
            StateModel::InPlace,
        ),
        (
            "Molecule.without_hydrogens",
            "without_hydrogens",
            "withoutHydrogens",
            "crate::Molecule",
            StateModel::ValueReturning,
        ),
        (
            "Molecule.remove_hydrogens_",
            "remove_hydrogens_",
            "removeHydrogens",
            "()",
            StateModel::InPlace,
        ),
    ] {
        let entry = entry(semantic_id);
        assert_eq!(entry.owner, BindingOwner::Molecule);
        assert!(entry.rust_path.replace(' ', "").ends_with(rust_name));
        assert_eq!(entry.python_name, rust_name);
        assert_eq!(entry.javascript_name, javascript);
        assert_eq!(entry.feature, "hydrogens");
        assert_eq!(entry.exposure, BindingExposure::Public);
        assert_eq!(entry.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(entry.parity, BindingParity::RequiredNow);
        let callable = entry.callable.expect("hydrogen callable metadata");
        assert_eq!(callable.kind, BindingKind::Instance);
        assert!(callable.parameters.is_empty());
        assert_eq!(callable.output_type.replace(' ', ""), output);
        assert_eq!(
            callable.error_type.map(|name| name.replace(' ', "")),
            Some("crate::OperationError".to_owned())
        );
        assert_eq!(callable.state_model, state);
        assert_eq!(callable.operation_semantic_id, Some(rust_name));
    }
}

#[test]
fn schema_represents_every_closed_enum_and_parameter_payload_branch() {
    static PARAMETERS: [BindingParameterContract; 2] = [
        BindingParameterContract {
            name: "required_input",
            type_name: "&str",
            default: BindingDefault::Required,
        },
        BindingParameterContract {
            name: "configured_input",
            type_name: "bool",
            default: BindingDefault::Value("false"),
        },
    ];
    let callable = BindingCallableContract {
        kind: BindingKind::Static,
        parameters: &PARAMETERS,
        output_type: "u32",
        error_type: Some("OperationError"),
        state_model: StateModel::ValueReturning,
        operation_semantic_id: Some("from_example"),
    };
    assert_eq!(callable.kind, BindingKind::Static);
    assert_eq!(callable.parameters[0].default, BindingDefault::Required);
    assert_eq!(
        callable.parameters[1].default,
        BindingDefault::Value("false")
    );
    assert_eq!(callable.error_type, Some("OperationError"));
    assert_eq!(callable.state_model, StateModel::ValueReturning);

    assert_eq!(
        [
            BindingOwner::Molecule,
            BindingOwner::Module,
            BindingOwner::Type
        ]
        .len(),
        3
    );
    assert_eq!(
        [
            BindingKind::Instance,
            BindingKind::Static,
            BindingKind::Module
        ]
        .len(),
        3
    );
    assert_eq!(
        [BindingExposure::Registered, BindingExposure::Public].len(),
        2
    );
    assert_eq!(
        [
            BindingTypeRole::Value,
            BindingTypeRole::Parameter,
            BindingTypeRole::Result,
            BindingTypeRole::Error
        ]
        .len(),
        4
    );
    assert_eq!(
        [
            BindingSupport::Unsupported,
            BindingSupport::PreservedOnly,
            BindingSupport::Experimental,
            BindingSupport::Supported,
            BindingSupport::SupportedWithRdkitParity
        ]
        .len(),
        5
    );
    assert_eq!(
        [
            BindingParity::NotApplicable,
            BindingParity::RequiredWhenSupported,
            BindingParity::RequiredNow
        ]
        .len(),
        3
    );
    assert_eq!(
        [
            StateModel::ValueReturning,
            StateModel::InPlace,
            StateModel::ReadOnly
        ]
        .len(),
        3
    );
}

#[test]
fn registry_source_has_one_declaration_and_no_legacy_schema_or_cfg_products() {
    let registry = include_str!("../src/binding_contract/registry.rs");
    let types = include_str!("../src/binding_contract/types.rs");

    assert_eq!(registry.matches("binding_contract!").count(), 1);
    assert_eq!(registry.matches("pub static BINDING_CONTRACT").count(), 1);
    assert!(!registry.contains("cfg(all"));
    assert!(!registry.contains("cfg(not"));
    assert!(!registry.contains("ReturnKind"));
    assert!(!types.contains("ReturnKind"));
    assert!(!types.contains("type Binding"));
}
