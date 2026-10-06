#![allow(unexpected_cfgs)]

use std::collections::HashSet;
#[cfg(feature = "cap-forcefields")]
use std::path::PathBuf;
#[cfg(feature = "cap-forcefields")]
use std::process::{Command, Output};

use cosmolkit::{
    BINDING_CONTRACT, BindingCallableContract, BindingDefault, BindingItem, BindingKind,
    BindingOwner, BindingParameterContract, BindingTypeRole, FunctionStatus, StateModel,
};

fn entry(semantic_id: &str) -> &'static cosmolkit::BindingContractEntry {
    BINDING_CONTRACT
        .iter()
        .find(|entry| entry.semantic_id == semantic_id)
        .unwrap_or_else(|| panic!("missing binding contract entry {semantic_id}"))
}

fn canonical_svg_identity_status() -> FunctionStatus {
    FunctionStatus::ParityWithDifferences {
        reference: "RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8 MolDraw2DSVG",
        explanation: "ROOT-SVG-CANONICAL-METADATA-20261005: public SVG declares ck=https://kit.cosmol.org/ instead of the pinned source renderer identity; all other drawing bytes retain their source comparison.",
    }
}

#[cfg(cosmolkit_uff_param_api_probe)]
mod uff_param_api_compile_probe {
    #[cfg(cosmolkit_uff_param_api_case = "query_available")]
    pub fn query_available(
        molecule: &cosmolkit::Molecule,
    ) -> Result<bool, cosmolkit::UffParameterQueryError> {
        let _kind_method: fn(&cosmolkit::UffParameterError) -> cosmolkit::UffParameterErrorKind =
            cosmolkit::UffParameterError::kind;
        let _kind_values = [
            cosmolkit::UffParameterErrorKind::Preparation,
            cosmolkit::UffParameterErrorKind::ParameterTable,
            cosmolkit::UffParameterErrorKind::Typing,
        ];
        let _query_method: fn(
            &cosmolkit::Molecule,
        ) -> Result<bool, cosmolkit::UffParameterQueryError> =
            cosmolkit::Molecule::uff_has_all_molecule_params;
        molecule.uff_has_all_molecule_params()
    }

    #[cfg(cosmolkit_uff_param_api_case = "base_api_available")]
    pub fn base_api_available() {
        let _molecule = cosmolkit::Molecule::new();
    }

    #[cfg(cosmolkit_uff_param_api_case = "query_unavailable")]
    pub fn query_unavailable(molecule: &cosmolkit::Molecule) {
        use cosmolkit::{UffParameterError, UffParameterErrorKind, UffParameterQueryError};

        let _: Option<UffParameterError> = None;
        let _: Option<UffParameterErrorKind> = None;
        let _: Option<UffParameterQueryError> = None;
        let _ = molecule.uff_has_all_molecule_params();
    }

    #[cfg(cosmolkit_uff_param_api_case = "forcefields_only_exclusions")]
    pub fn forcefields_only_exclusions(molecule: &cosmolkit::Molecule) {
        use cosmolkit::{
            ForceFieldError, ForceFieldOptions, RingSearchParams, ValenceParams,
            mmff_has_all_molecule_params, mmff_optimize, uff_has_all_molecule_params,
        };

        let _ = (
            ForceFieldError::Unsupported,
            ForceFieldOptions::default(),
            RingSearchParams::default(),
            ValenceParams::default(),
            mmff_has_all_molecule_params,
            mmff_optimize,
            uff_has_all_molecule_params,
        );
        let _ = molecule.with_assigned_valence();
        let _ = molecule.with_assigned_rings();
        let _ = cosmolkit::forcefields::UffParameterQueryError::Cache;
    }
}

#[cfg(feature = "cap-forcefields")]
fn uff_param_p10_compile_check(case: &str, features: Option<&str>) -> Output {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let workspace = manifest_dir
        .parent()
        .and_then(|path| path.parent())
        .expect("cosmolkit crate must be nested under the workspace crates directory");
    let target = workspace.join("target/uff-param-api-compile");
    let inherited = std::env::var("RUSTFLAGS").unwrap_or_default();
    let rustflags = format!(
        "{inherited} --cfg cosmolkit_uff_param_api_probe --cfg=cosmolkit_uff_param_api_case=\"{case}\""
    );
    let mut command = Command::new(std::env::var("CARGO").unwrap_or_else(|_| "cargo".into()));
    command
        .current_dir(workspace)
        .env("CARGO_TARGET_DIR", target)
        .env("CARGO_INCREMENTAL", "0")
        .env("RUSTFLAGS", rustflags)
        .args([
            "check",
            "--quiet",
            "-p",
            "cosmolkit",
            "--test",
            "migration_run_schema",
            "--no-default-features",
        ]);
    if let Some(features) = features {
        command.args(["--features", features]);
    }
    command
        .output()
        .expect("run the external cosmolkit feature-isolation compile proof")
}

#[cfg(feature = "cap-io")]
#[test]
fn sdf_reader_contracts_distinguish_concrete_and_preserving_results() {
    for (id, output) in [
        ("Molecule.from_sdf", "crate::Molecule"),
        ("Molecule.from_sdf_with_params", "crate::Molecule"),
        ("SdfRecord.from_sdf", "crate::SdfRecord"),
        ("SdfRecord.from_sdf_with_params", "crate::SdfRecord"),
    ] {
        let contract = entry(id);
        assert_eq!(contract.feature, "cap-io");
        assert_eq!(contract.status, FunctionStatus::Experimental);
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
        assert_eq!(contract.feature, "cap-io");
        assert_eq!(contract.status, FunctionStatus::Experimental);
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
        assert_eq!(contract.status, FunctionStatus::Experimental);
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
        "types.TopologyBlock",
        "types.CoordinateBlock",
        "types.Conformer2D",
        "types.Conformer3D",
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
        "types.PropertyValue",
        "types.PropertyValueKind",
        "types.PropertyValueError",
        "PropertyValueError.expected",
        "PropertyValueError.actual",
        "PropertyValue.kind",
        "PropertyValue.as_string",
        "PropertyValue.as_int",
        "PropertyValue.as_uint",
        "PropertyValue.as_int_vector",
        "PropertyValue.as_double",
        "PropertyValue.as_bool",
        "types.TemplateAttachment",
        "types.TemplateAttachmentOrder",
        "types.TemplateAttachmentOrderError",
        "TemplateAttachment.target",
        "TemplateAttachment.label",
        "TemplateAttachmentOrder.entries",
        "Atom.template_attachment_order",
        "types.FunctionStatus",
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
    if cfg!(feature = "cap-forcefields") {
        let operation_error = expected
            .iter()
            .position(|semantic_id| *semantic_id == "types.OperationError")
            .expect("the existing OperationError registry fixture is present")
            + 1;
        expected.splice(
            operation_error..operation_error,
            [
                "types.UffParameterQueryError",
                "types.UffParameterError",
                "types.UffParameterErrorKind",
                "UffParameterError.kind",
            ],
        );
    }
    if cfg!(feature = "cap-matrices") {
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
    if cfg!(feature = "cap-depict") {
        expected.extend([
            "types.Coordinate2DParams",
            "types.Coordinate2DError",
            "types.Coordinate2DTemplateError",
            "types.Coordinate2DLayoutError",
        ]);
    }
    expected.push("Molecule.has_2d_coordinates");
    if cfg!(feature = "cap-depict") {
        expected.extend([
            "Molecule.with_2d_coordinates",
            "Molecule.with_2d_coordinates_with_params",
            "types.DrawingError",
            "Molecule.to_svg",
            "Molecule.to_png",
            "Molecule.compute_2d_coordinates_",
            "Molecule.compute_2d_coordinates_with_params_",
            "types.DrawingWriteError",
            "Molecule.write_svg",
            "Molecule.write_png",
        ]);
    }
    if cfg!(feature = "cap-transforms") {
        expected.extend([
            "types.AtomPositionParams",
            "types.TransformError",
            "Molecule.with_atom_position",
            "Molecule.with_atom_position_with_params",
            "Molecule.set_atom_position_",
            "Molecule.set_atom_position_with_params_",
        ]);
    }
    if cfg!(feature = "cap-sanitize") {
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
    if cfg!(feature = "cap-kekulize") {
        expected.extend([
            "types.KekulizeParams",
            "types.KekulizeError",
            "Molecule.with_kekulized_bonds",
            "Molecule.with_kekulized_bonds_with_params",
            "Molecule.kekulize_bonds_",
            "Molecule.kekulize_bonds_with_params_",
        ]);
    }
    if cfg!(feature = "cap-aromaticity") {
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
    if cfg!(feature = "cap-valence") {
        expected.extend(["types.ValenceModel", "types.ValenceParams"]);
    }
    // Shared model vocabulary exists independently of valence computation.
    expected.push("types.ValenceError");
    if cfg!(feature = "cap-valence") {
        expected.extend([
            "Molecule.with_assigned_valence",
            "Molecule.with_assigned_valence_with_params",
            "Molecule.assign_valence_",
            "Molecule.assign_valence_with_params_",
            "Molecule.has_valence_violation",
        ]);
    }
    if cfg!(feature = "cap-radicals") {
        expected.extend([
            "Molecule.with_assigned_radicals",
            "Molecule.assign_radicals_",
        ]);
    }
    if cfg!(feature = "cap-rings") {
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
    if cfg!(feature = "cap-stereo") {
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
    if cfg!(feature = "cap-descriptors") {
        expected.extend([
            "Molecule.hall_kier_alpha",
            "Molecule.hall_kier_alpha_with_contributions",
            "Molecule.kappa_1",
            "Molecule.kappa_2",
            "Molecule.kappa_3",
            "Molecule.phi",
            "Molecule.mqns",
            "Molecule.chi_0_v",
            "Molecule.chi_1_v",
            "Molecule.chi_2_v",
            "Molecule.chi_3_v",
            "Molecule.chi_4_v",
            "Molecule.chi_n_v",
            "Molecule.chi_0_n",
            "Molecule.chi_1_n",
            "Molecule.chi_2_n",
            "Molecule.chi_3_n",
            "Molecule.chi_4_n",
            "Molecule.chi_n_n",
            "Molecule.molecular_weight",
            "Molecule.exact_molecular_weight",
            "Molecule.molecular_formula",
            "types.DescriptorError",
            "types.DescriptorReadError",
            "Molecule.num_heavy_atoms",
            "Molecule.total_atom_count",
            "Molecule.num_rings",
            "Molecule.num_heterocycles",
            "Molecule.num_heteroatoms",
            "Molecule.num_hba",
            "Molecule.num_hbd",
            "Molecule.num_aromatic_rings",
            "Molecule.num_saturated_rings",
            "Molecule.num_aliphatic_rings",
            "Molecule.num_aromatic_heterocycles",
            "Molecule.num_aromatic_carbocycles",
            "Molecule.num_aliphatic_heterocycles",
            "Molecule.num_aliphatic_carbocycles",
            "Molecule.num_saturated_heterocycles",
            "Molecule.num_saturated_carbocycles",
            "Molecule.lipinski_hba",
            "Molecule.lipinski_hbd",
            "Molecule.fraction_csp3",
        ]);
    }
    if cfg!(feature = "cap-forcefields") {
        expected.extend([
            "Molecule.uff_has_all_molecule_params",
            "types.UffOptimizationParams",
            "types.UffOptimizationResult",
            "types.UffOptimizationError",
            "types.UffOptimizationErrorKind",
            "types.UffConformerOptimizationParams",
            "types.UffConformerOptimizationResult",
            "types.UffConformerResult",
            "UffOptimizationError.kind",
            "Molecule.with_uff_optimized",
            "Molecule.with_uff_optimized_with_params",
            "Molecule.with_uff_optimized_confs",
            "Molecule.with_uff_optimized_confs_with_params",
        ]);
    }
    if cfg!(feature = "cap-hydrogens") {
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
    if cfg!(feature = "cap-bio") {
        expected.extend([
            "types.ResidueInfoKind",
            "ResidueInfoKind.name",
            "types.ResidueCode",
            "types.ResidueCodeParseError",
            "ResidueCodeParseError.input",
            "types.ResidueIdentity",
            "ResidueIdentity.new",
            "ResidueIdentity.name",
            "ResidueIdentity.code",
            "ResidueIdentity.info",
            "ResidueIdentity.is_tabulated",
            "types.ResidueInfo",
            "ResidueInfo.canonical_one_letter_code",
            "ResidueInfo.parent_standard_code",
            "ResidueInfo.is_modified_amino_acid",
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
            "types.ProteinSelectionSummary",
            "types.ProteinChainIter",
            "types.ProteinResidueIter",
            "types.ProteinAtomIter",
            "types.ProteinChainRef",
            "types.ProteinResidueRef",
            "types.ProteinAtomRef",
            "BioStructure.from_parts",
            "BioStructure.validate_parts",
            "BioStructure.validate",
            "BioStructure.into_parts",
            "BioStructure.input_format",
            "BioStructure.models",
            "BioStructure.chains",
            "BioStructure.residues",
            "BioStructure.atoms",
            "BioStructure.entities",
            "BioStructure.connections",
            "BioStructure.cispeps",
            "BioStructure.mod_residues",
            "BioStructure.helices",
            "BioStructure.sheets",
            "BioStructure.ncs_operators",
            "BioStructure.assemblies",
            "BioStructure.metadata",
            "BioStructure.source_state",
            "BioStructure.coordinates",
            "BioStructure.crystal",
            "BioStructure.find_entity",
            "BioStructure.find_entity_of_subchain",
            "BioStructure.find_atom",
            "BioStructure.atom_by_altloc",
            "types.BioStrand",
            "types.BioSoftwareItem",
            "types.BioReflectionsInfo",
            "types.BioBasicRefinementInfo",
            "types.BioRefinementRestraint",
            "types.BioExperimentInfo",
            "types.BioDiffractionInfo",
            "types.BioExperimentalCrystalInfo",
            "types.BioTlsSelection",
            "types.BioTlsGroup",
            "types.BioRefinementInfo",
            "types.BioConnectionKind",
            "types.BioAsu",
            "types.BioHelixClass",
            "types.BioSoftwareClassification",
            "types.AtomAddress",
            "types.ResidueAddress",
            "types.BioConnection",
            "types.BioCisPep",
            "types.BioModRes",
            "types.BioHelix",
            "types.BioSheet",
            "types.BioMetadata",
            "types.BioStructureSourceState",
            "BioStructure.num_models",
            "BioStructure.name",
            "BioStructure.has_origx",
            "BioStructure.origx",
            "BioStructure.ncs_oper_identity_id",
            "BioStructure.resolution",
            "BioStructure.ter_status",
            "BioStructure.num_chains",
            "BioStructure.num_residues",
            "BioStructure.num_atoms",
            "BioStructure.num_entities",
            "BioStructure.atom_position",
            "BioStructure.residue_atoms",
            "BioStructure.protein",
            "Protein.input_format",
            "Protein.as_bio_structure",
            "Protein.into_bio_structure",
            "Protein.selection_summary",
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
    expected.push("types.CoordinateDimension");
    if cfg!(feature = "cap-io") {
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
    if cfg!(feature = "cap-smiles") {
        expected.extend([
            "types.SmilesParseParams",
            "types.SmilesError",
            "types.SmilesStereoError",
            "Molecule.from_smiles",
            "Molecule.from_smiles_with_params",
        ]);
    }

    if cfg!(feature = "cap-fingerprints") {
        // Exact Morgan configuration and FingerprintAdditionalOutput surface precedes
        // the existing sparse-count entries in registry order.
        expected.splice(
            0..0,
            [
                "types.MorganParams",
                "types.MorganInvariants",
                "types.MorganFingerprintParams",
                "types.MorganReadError",
                "types.FingerprintAdditionalOutput",
                "types.Fingerprint",
                "types.SparseBitFingerprint",
                "types.SparseCountFingerprint",
                "types.SparseCountFingerprint32",
                "types.FingerprintError",
                "Molecule.morgan_sparse_count_fingerprint",
                "Molecule.morgan_sparse_count_fingerprint_with_params",
                "Molecule.morgan_sparse_fingerprint",
                "Molecule.morgan_sparse_fingerprint_with_params",
                "Molecule.morgan_count_fingerprint",
                "Molecule.morgan_count_fingerprint_with_params",
                "Molecule.morgan_fingerprint",
                "Molecule.morgan_fingerprint_with_params",
                "FingerprintAdditionalOutput.new",
                "FingerprintAdditionalOutput.default",
                "FingerprintAdditionalOutput.allocate_atom_counts",
                "FingerprintAdditionalOutput.allocate_atom_to_bits",
                "FingerprintAdditionalOutput.allocate_bit_info_map",
                "FingerprintAdditionalOutput.allocate_bit_paths",
                "FingerprintAdditionalOutput.allocate_atoms_per_bit",
                "FingerprintAdditionalOutput.atom_counts",
                "FingerprintAdditionalOutput.atom_to_bits",
                "FingerprintAdditionalOutput.bit_info_map",
                "FingerprintAdditionalOutput.bit_paths",
                "FingerprintAdditionalOutput.atoms_per_bit",
                "Fingerprint.n_bits",
                "Fingerprint.on_bits",
                "SparseBitFingerprint.n_bits",
                "SparseBitFingerprint.on_bits",
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
    if cfg!(feature = "cap-bio") {
        // Approved associated readers and lightweight operations precede the
        // existing registry. Preserve exact ordering and feature-local coverage.
        expected.splice(
            0..0,
            [
                "types.BioOperationError",
                "BioStructure.with_translated_coordinates",
                "BioStructure.translate_",
                "Protein.with_translated_coordinates",
                "Protein.translate_",
                "BioStructure.with_selection",
                "BioStructure.retain_selection_",
                "Protein.with_selection",
                "Protein.retain_selection_",
                "types.BioSelectionCopyError",
                "types.BioSelectionCopyCause",
                "types.BioRowTraverseError",
                "types.BioRowModelError",
                "types.BioRowChainError",
                "types.BioPdbReadParams",
                "types.BioPdbReadError",
                "types.BioPdbReadStage",
                "types.BioMmcifReadError",
                "types.BioMmcifReadStage",
                "types.BioReadParams",
                "types.BioReadError",
                "types.ProteinReadError",
                "types.BioMmcifWriteParams",
                "types.BioMmcifWriteError",
                "BioStructure.to_mmcif_with_params",
                "BioStructure.to_mmcif",
                "BioStructure.write_mmcif_with_params",
                "BioStructure.write_mmcif",
                "types.BioPdbWriteParams",
                "types.BioPdbWriteError",
                "BioStructure.to_pdb_with_params",
                "BioStructure.to_pdb",
                "BioStructure.write_pdb_with_params",
                "BioStructure.write_pdb",
                "BioStructure.from_text_with_params",
                "BioStructure.from_text",
                "BioStructure.read_with_format",
                "BioStructure.read",
                "Protein.from_text_with_params",
                "Protein.from_text",
                "Protein.read_with_format",
                "Protein.read",
                "BioStructure.from_pdb",
                "BioStructure.from_pdb_with_params",
                "BioStructure.from_mmcif",
                "Protein.from_pdb",
                "Protein.from_pdb_with_params",
                "Protein.from_mmcif",
                "BioPdbReadError.stage",
                "BioPdbReadError.line_number",
                "BioPdbReadError.record_tag",
                "BioMmcifReadError.stage",
            ],
        );
    }
    if cfg!(feature = "cap-smiles") {
        expected.splice(
            0..0,
            [
                "types.SmilesWriteParams",
                "types.CxSmilesWriteParams",
                "types.CxSmilesFields",
                "types.CxCoordinateSelection",
                "types.RandomSmilesWriteParams",
                "types.FragmentSmilesWriteParams",
                "types.FragmentCxSmilesWriteParams",
                "types.SmilesWriteError",
                "Molecule.to_smiles",
                "Molecule.to_smiles_with_params",
                "Molecule.to_cx_smiles",
                "Molecule.to_cx_smiles_with_params",
                "Molecule.to_fragment_smiles",
                "Molecule.to_fragment_smiles_with_params",
                "Molecule.to_fragment_cx_smiles",
                "Molecule.to_fragment_cx_smiles_with_params",
                "Molecule.to_random_smiles",
                "Molecule.to_random_smiles_with_params",
                "types.BioSelection",
                "BioSelection.from_cid",
                "Protein.selected_atom_ids",
                "BioStructure.selected_atom_ids",
                "BioSelection.to_cid",
                "types.BioSelectionParseError",
                "types.BioSelectionMatchError",
            ],
        );
    }
    if cfg!(feature = "cap-valence") {
        expected.insert(0, "module.element_info");
    }
    expected.splice(
        0..0,
        [
            "types.Element",
            "Element.from_atomic_number",
            "Element.from_symbol",
            "Element.atomic_number",
            "Element.symbol",
            "types.ElementInfo",
        ],
    );
    // Complete public blocks are static independent fixtures. Keep every original
    // identity, ordering assertion and feature-local consumer unchanged.
    // Source evidence: binding_contract/registry.rs complete UFF/MMFF, TAU,
    // SEARCH, FP, descriptors and canonical atom/bond/builder declarations.
    {
        let mut additions = Vec::new();
        if cfg!(feature = "cap-forcefields") {
            additions.extend([
                "types.UffEvaluationParams",
                "types.UffEnergyGradient",
                "Molecule.uff_energy_gradient",
                "Molecule.uff_energy_gradient_with_params",
                "UffEnergyGradient.energy",
                "UffEnergyGradient.gradient",
                "UffOptimizationResult.molecule",
                "UffOptimizationResult.status_code",
                "UffOptimizationResult.needs_more",
                "UffOptimizationResult.energy",
                "UffConformerResult.conformer_id",
                "UffConformerResult.status_code",
                "UffConformerResult.needs_more",
                "UffConformerResult.energy",
                "UffConformerOptimizationResult.molecule",
                "UffConformerOptimizationResult.conformer_results",
            ]);
        }
        if !additions.is_empty() {
            let insertion = expected
                .iter()
                .position(|id| *id == "types.Element")
                .expect("original fixture anchor types.Element");
            expected.splice(insertion..insertion, additions);
        }
    }
    {
        let mut additions = Vec::new();
        additions.extend(["types.QueryGraph"]);
        if cfg!(feature = "cap-forcefields") {
            additions.extend([
                "MmffAtomProperties.atom_type",
                "MmffAtomProperties.formal_charge",
                "MmffAtomProperties.partial_charge",
                "types.MmffEvaluationParams",
                "types.MmffEnergyGradient",
                "Molecule.mmff_energy_gradient",
                "Molecule.mmff_energy_gradient_with_params",
                "MmffEnergyGradient.energy",
                "MmffEnergyGradient.gradient",
                "MmffOptimizeMoleculeResult.molecule",
                "MmffOptimizeMoleculeResult.needs_more",
                "MmffOptimizeMoleculeResult.status_code",
                "MmffOptimizeMoleculeConfsResult.molecule",
                "MmffOptimizeMoleculeConfsResult.conformer_results",
                "MmffOptimizeMoleculeConfResult.needs_more",
                "MmffOptimizeMoleculeConfResult.status_code",
                "MmffOptimizeMoleculeConfResult.energy",
                "MmffProperties.is_valid",
                "MmffProperties.variant",
                "MmffProperties.atoms",
                "MmffProperties.atom_type",
                "MmffProperties.formal_charge",
                "MmffProperties.partial_charge",
                "types.MmffOptimizationParams",
                "types.MmffConformerOptimizationParams",
                "types.MmffOptimizeMoleculeResult",
                "types.MmffOptimizeMoleculeConfResult",
                "types.MmffOptimizeMoleculeConfsResult",
                "types.MmffOptimizationError",
                "Molecule.with_mmff_optimized",
                "Molecule.with_mmff_optimized_with_params",
                "Molecule.with_mmff_optimized_confs",
                "Molecule.with_mmff_optimized_confs_with_params",
            ]);
        }
        if cfg!(feature = "cap-tautomer") {
            additions.extend([
                "types.TautomerParams",
                "types.TautomerScoreParams",
                "types.TautomerScoreTerm",
                "types.TautomerScore",
                "types.TautomerEnumeration",
                "types.TautomerEnumerationStatus",
                "types.TautomerRunError",
                "types.TautomerCatalogError",
                "types.TautomerMoleculeView",
                "types.TautomerProgress",
                "Molecule.enumerate_tautomers",
                "Molecule.enumerate_tautomers_with_params",
                "Molecule.canonical_tautomer",
                "Molecule.canonical_tautomer_with_params",
                "Molecule.tautomer_score",
                "Molecule.tautomer_score_with_params",
            ]);
        }
        if cfg!(feature = "cap-search") {
            additions.extend([
                "types.SmartsParseParams",
                "types.SmartsParseError",
                "types.SubstructMatchParams",
                "types.SubstructMatchError",
                "types.MatchResult",
                "types.CompiledQuery",
                "types.QueryCompileError",
                "types.MatchError",
                "types.SmartsWriteParams",
                "types.SmartsWriteError",
                "search.parse_smarts",
                "search.parse_smarts_with_params",
                "QueryGraph.from_smarts",
                "QueryGraph.from_smarts_with_params",
                "search.compile_query",
                "search.write_smarts",
                "search.write_cx_smarts",
                "Molecule.substruct_match",
                "Molecule.substruct_matches",
                "Molecule.has_substruct_match",
                "Molecule.substruct_matches_with_params",
                "Molecule.substruct_matches_compiled",
            ]);
        }
        let insertion = expected
            .iter()
            .position(|id| *id == "types.ElementInfo")
            .unwrap()
            + 1
            + usize::from(cfg!(feature = "cap-valence"));
        expected.splice(insertion..insertion, additions);
    }
    {
        let mut additions = Vec::new();
        if cfg!(feature = "cap-fingerprints") {
            additions.extend([
                "types.AtomPairParams",
                "types.AtomPairFingerprintParams",
                "types.AtomPairAtomInvariantsGenerator",
                "types.AtomPairReadError",
                "Molecule.atom_pair_fingerprint",
                "Molecule.atom_pair_fingerprint_with_params",
                "Molecule.atom_pair_sparse_fingerprint",
                "Molecule.atom_pair_sparse_fingerprint_with_params",
                "Molecule.atom_pair_count_fingerprint",
                "Molecule.atom_pair_count_fingerprint_with_params",
                "Molecule.atom_pair_sparse_count_fingerprint",
                "Molecule.atom_pair_sparse_count_fingerprint_with_params",
                "AtomPairAtomInvariantsGenerator.info_string",
                "AtomPairAtomInvariantsGenerator.to_json",
                "types.TopologicalTorsionParams",
                "types.TopologicalTorsionFingerprintParams",
                "types.TopologicalTorsionReadError",
                "Molecule.topological_torsion_fingerprint",
                "Molecule.topological_torsion_fingerprint_with_params",
                "Molecule.topological_torsion_sparse_fingerprint",
                "Molecule.topological_torsion_sparse_fingerprint_with_params",
                "Molecule.topological_torsion_count_fingerprint",
                "Molecule.topological_torsion_count_fingerprint_with_params",
                "Molecule.topological_torsion_sparse_count_fingerprint",
                "Molecule.topological_torsion_sparse_count_fingerprint_with_params",
            ]);
        }
        if !additions.is_empty() {
            let insertion = expected
                .iter()
                .position(|id| *id == "types.MorganParams")
                .expect("original fixture anchor types.MorganParams");
            expected.splice(insertion..insertion, additions);
        }
    }
    {
        let mut additions = vec!["types.AtomMetadata"];
        if cfg!(feature = "cap-valence") {
            additions.push("Molecule.atom_metadata");
        }
        let insertion = expected
            .iter()
            .position(|id| *id == "types.Molecule")
            .expect("original fixture anchor types.Molecule");
        expected.splice(insertion..insertion, additions);
    }
    {
        let mut additions = Vec::new();
        if cfg!(feature = "cap-descriptors") {
            additions.extend(["Molecule.chi_0", "Molecule.chi_1"]);
        }
        if !additions.is_empty() {
            let insertion = expected
                .iter()
                .position(|id| *id == "Molecule.hall_kier_alpha")
                .expect("original fixture anchor Molecule.hall_kier_alpha");
            expected.splice(insertion..insertion, additions);
        }
    }
    {
        let mut additions = Vec::new();
        if cfg!(feature = "cap-descriptors") {
            additions.extend([
                "types.RotatableBondsOptions",
                "Molecule.num_amide_bonds",
                "Molecule.num_spiro_atoms",
                "Molecule.num_bridgehead_atoms",
                "Molecule.num_atom_stereo_centers",
                "Molecule.num_unspecified_atom_stereo_centers",
                "Molecule.num_rotatable_bonds",
                "Molecule.num_rotatable_bonds_with_params",
                "Molecule.molecular_weight_with_params",
                "Molecule.exact_molecular_weight_with_params",
                "Molecule.molecular_formula_with_params",
            ]);
        }
        if !additions.is_empty() {
            let insertion = expected
                .iter()
                .position(|id| *id == "types.DescriptorError")
                .expect("original fixture anchor types.DescriptorError");
            expected.splice(insertion..insertion, additions);
        }
    }
    {
        let mut additions = Vec::new();
        if cfg!(feature = "cap-forcefields") {
            additions.extend([
                "types.MmffProperties",
                "types.MmffPropertiesParams",
                "types.MmffAtomProperties",
                "types.MmffVariant",
                "types.MmffMolPropertiesError",
                "Molecule.mmff_has_all_molecule_params",
                "Molecule.mmff_properties",
                "Molecule.mmff_properties_with_params",
            ]);
        }
        if !additions.is_empty() {
            let insertion = expected
                .iter()
                .position(|id| *id == "Molecule.uff_has_all_molecule_params")
                .expect("original fixture anchor Molecule.uff_has_all_molecule_params");
            expected.splice(insertion..insertion, additions);
        }
    }
    {
        let mut additions = Vec::new();
        if cfg!(feature = "cap-descriptors") {
            additions.extend([
                "types.CrippenTotals",
                "types.LabuteAsaContributions",
                "Molecule.crippen_descriptors",
                "Molecule.labute_asa",
                "Molecule.labute_asa_contributions",
                "Molecule.tpsa",
                "Molecule.slogp_vsa",
                "Molecule.smr_vsa",
                "Molecule.slogp_vsa_1",
                "Molecule.slogp_vsa_2",
                "Molecule.slogp_vsa_3",
                "Molecule.slogp_vsa_4",
                "Molecule.slogp_vsa_5",
                "Molecule.slogp_vsa_6",
                "Molecule.slogp_vsa_7",
                "Molecule.slogp_vsa_8",
                "Molecule.slogp_vsa_9",
                "Molecule.slogp_vsa_10",
                "Molecule.slogp_vsa_11",
                "Molecule.slogp_vsa_12",
                "Molecule.smr_vsa_1",
                "Molecule.smr_vsa_2",
                "Molecule.smr_vsa_3",
                "Molecule.smr_vsa_4",
                "Molecule.smr_vsa_5",
                "Molecule.smr_vsa_6",
                "Molecule.smr_vsa_7",
                "Molecule.smr_vsa_8",
                "Molecule.smr_vsa_9",
                "Molecule.smr_vsa_10",
                "Molecule.crippen_descriptors_with_params",
                "Molecule.labute_asa_with_params",
                "Molecule.labute_asa_contributions_with_params",
                "Molecule.tpsa_with_params",
                "Molecule.slogp_vsa_with_params",
                "Molecule.smr_vsa_with_params",
                "Molecule.qed",
                "Molecule.chi_0_v_with_params",
                "Molecule.chi_1_v_with_params",
                "Molecule.chi_2_v_with_params",
                "Molecule.chi_3_v_with_params",
                "Molecule.chi_4_v_with_params",
                "Molecule.chi_n_v_with_params",
                "Molecule.chi_0_n_with_params",
                "Molecule.chi_1_n_with_params",
                "Molecule.chi_2_n_with_params",
                "Molecule.chi_3_n_with_params",
                "Molecule.chi_4_n_with_params",
                "Molecule.chi_n_n_with_params",
            ]);
        }
        additions.extend([
            "types.Atom",
            "types.Bond",
            "types.BondOrder",
            "types.BondDirection",
            "types.BondStereo",
            "types.ChiralTag",
            "types.AtomSpec",
            "types.BondSpec",
        ]);
        additions.extend([
            "Atom.id",
            "Atom.element",
            "Atom.atomic_number",
            "Atom.formal_charge",
            "Atom.chiral_tag",
            "Atom.chiral_tag_code",
            "Atom.chiral_tag_name",
            "Atom.isotope",
            "Atom.atom_map",
            "Atom.is_aromatic",
            "Atom.explicit_hydrogens",
            "Atom.no_implicit",
            "Atom.radical_electrons",
            "Bond.id",
            "Bond.begin",
            "Bond.end",
            "Bond.order",
            "Bond.order_code",
            "Bond.order_name",
            "Bond.direction",
            "Bond.direction_code",
            "Bond.direction_name",
            "Bond.stereo",
            "Bond.stereo_code",
            "Bond.stereo_name",
            "Bond.stereo_atoms",
            "Bond.is_aromatic",
            "Atom.cip_descriptor",
            "Atom.cip_neighbor_order",
            "Atom.cip_rank",
            "Bond.cip_descriptor",
            "Bond.cip_neighbor_order",
        ]);
        if cfg!(feature = "cap-stereo") {
            additions.extend(["Molecule.cip_computed"]);
        }
        additions.extend([
            "AtomSpec.new",
            "BondSpec.new",
            "AtomSpec.with_formal_charge",
            "AtomSpec.with_explicit_hydrogens",
            "AtomSpec.with_atom_map",
            "AtomSpec.with_isotope",
            "AtomSpec.with_no_implicit",
        ]);
        if cfg!(feature = "cap-sanitize") {
            additions.extend(["Molecule.sanitize_", "Molecule.sanitize_with_params_"]);
        }
        expected.extend(additions);
    }
    // Complete source-audited FP v7 public closure follows the existing registry.
    if cfg!(feature = "cap-fingerprints") {
        expected.extend([
            "types.LegacyTopologicalTorsionParams",
            "types.TopologicalTorsionFingerprintGenerator",
            "types.TopologicalTorsionSettings",
            "types.TopologicalTorsionCallParams",
            "TopologicalTorsionFingerprintGenerator.new",
            "TopologicalTorsionFingerprintGenerator.from_json",
            "TopologicalTorsionFingerprintGenerator.settings",
            "TopologicalTorsionFingerprintGenerator.info_string",
            "TopologicalTorsionFingerprintGenerator.to_json",
            "TopologicalTorsionFingerprintGenerator.fingerprints",
            "TopologicalTorsionFingerprintGenerator.sparse_fingerprints",
            "TopologicalTorsionFingerprintGenerator.counts",
            "TopologicalTorsionFingerprintGenerator.sparse_counts",
            "TopologicalTorsionSettings.torsion_atom_count",
            "TopologicalTorsionSettings.set_torsion_atom_count",
            "TopologicalTorsionSettings.only_shortest_paths",
            "TopologicalTorsionSettings.set_only_shortest_paths",
            "TopologicalTorsionSettings.include_chirality",
            "TopologicalTorsionSettings.set_include_chirality",
            "TopologicalTorsionSettings.count_simulation",
            "TopologicalTorsionSettings.set_count_simulation",
            "TopologicalTorsionSettings.fp_size",
            "TopologicalTorsionSettings.set_fp_size",
            "TopologicalTorsionSettings.bits_per_feature",
            "TopologicalTorsionSettings.set_bits_per_feature",
            "TopologicalTorsionSettings.count_bounds",
            "TopologicalTorsionSettings.set_count_bounds",
            "TopologicalTorsionSettings.params",
            "Molecule.topological_torsion_fingerprint_with_generator",
            "Molecule.topological_torsion_sparse_fingerprint_with_generator",
            "Molecule.topological_torsion_count_fingerprint_with_generator",
            "Molecule.topological_torsion_sparse_count_fingerprint_with_generator",
            "LegacyTopologicalTorsionParams.new",
            "Molecule.legacy_topological_torsion_sparse_count_fingerprint",
            "Molecule.legacy_topological_torsion_sparse_count_fingerprint_with_params",
            "Molecule.legacy_topological_torsion_count_fingerprint",
            "Molecule.legacy_topological_torsion_count_fingerprint_with_params",
            "Molecule.legacy_topological_torsion_fingerprint",
            "Molecule.legacy_topological_torsion_fingerprint_with_params",
        ]);
    }
    // Complete canonical feature blocks: literal identities preserve the original
    // ordering/uniqueness checks and independent feature-local expectations.
    if cfg!(feature = "cap-fingerprints") {
        expected.insert(0, "types.FingerprintPreparationError");
    }
    if cfg!(feature = "cap-fingerprints") {
        let position = expected
            .iter()
            .position(|id| *id == "types.Element")
            .expect("the existing types.Element fixture is present");
        expected.splice(
            position..position,
            [
                "types.AtomPairAtomCodeResult",
                "Molecule.with_atom_pair_atom_code",
            ],
        );
    }
    if cfg!(feature = "cap-tautomer") {
        let position = expected
            .iter()
            .position(|id| *id == "types.TautomerParams")
            .expect("the existing types.TautomerParams fixture is present");
        expected.splice(
            position..position,
            [
                "default_tautomer_score_terms",
                "TautomerParams.max_tautomers",
                "TautomerParams.set_max_tautomers",
                "TautomerParams.with_max_tautomers",
                "TautomerParams.max_transforms",
                "TautomerParams.set_max_transforms",
                "TautomerParams.with_max_transforms",
                "TautomerParams.remove_sp3_stereo",
                "TautomerParams.set_remove_sp3_stereo",
                "TautomerParams.with_remove_sp3_stereo",
                "TautomerParams.remove_bond_stereo",
                "TautomerParams.set_remove_bond_stereo",
                "TautomerParams.with_remove_bond_stereo",
                "TautomerParams.remove_isotopic_hydrogens",
                "TautomerParams.set_remove_isotopic_hydrogens",
                "TautomerParams.with_remove_isotopic_hydrogens",
                "TautomerParams.reassign_stereo",
                "TautomerParams.set_reassign_stereo",
                "TautomerParams.with_reassign_stereo",
                "TautomerParams.v1",
                "TautomerParams.from_transform_data",
                "TautomerParams.from_transform_file",
                "TautomerParams.transform_count",
                "TautomerParams.callback",
                "TautomerParams.set_callback",
                "TautomerParams.scorer",
                "TautomerParams.set_scorer",
                "TautomerEnumeration.len",
                "TautomerEnumeration.is_empty",
                "TautomerEnumeration.status",
                "TautomerEnumeration.modified_atoms",
                "TautomerEnumeration.modified_bonds",
                "TautomerEnumeration.canonical_smiles",
                "TautomerEnumeration.get",
                "TautomerEnumeration.iter",
                "TautomerEnumeration.entries",
                "tautomer.canonical_tautomer_from_molecules",
                "tautomer.canonical_tautomer_from_molecules_with_params",
                "TautomerMoleculeView.to_owned",
                "TautomerMoleculeView.num_atoms",
                "TautomerMoleculeView.num_bonds",
                "TautomerMoleculeView.atom_metadata",
                "TautomerMoleculeView.atom_degree",
                "TautomerMoleculeView.atoms",
                "TautomerMoleculeView.bonds",
                "TautomerMoleculeView.properties",
                "TautomerMoleculeView.atom",
                "TautomerMoleculeView.bond",
                "TautomerMoleculeView.to_smiles",
                "TautomerMoleculeView.tautomer_score",
                "TautomerProgress.to_owned",
                "TautomerProgress.len",
                "TautomerProgress.is_empty",
                "TautomerProgress.status",
                "TautomerProgress.num_transforms",
                "TautomerProgress.modified_atoms",
                "TautomerProgress.modified_bonds",
                "TautomerProgress.entries",
                "TautomerScoreTerm.new",
                "TautomerScoreTerm.name",
                "TautomerScoreTerm.smarts",
                "TautomerScoreTerm.score",
                "TautomerScore.ring",
                "TautomerScore.substructure",
                "TautomerScore.hetero_hydrogen",
                "TautomerScore.total",
                "TautomerEnumeration.canonical_tautomer",
                "TautomerEnumeration.canonical_tautomer_with_params",
            ],
        );
    }
    if cfg!(all(feature = "cap-io", feature = "cap-bio")) {
        let position = expected
            .iter()
            .position(|id| *id == "types.BioPdbReadError")
            .expect("the existing types.BioPdbReadError fixture is present");
        expected.splice(
            position..position,
            [
                "types.BioMoleculeParams",
                "types.BioMoleculeError",
                "types.BioMoleculeConversionError",
                "BioStructure.to_molecule_with_params",
                "BioStructure.to_molecule",
                "Protein.to_molecule_with_params",
                "Protein.to_molecule",
            ],
        );
    }
    if cfg!(feature = "cap-bio") {
        let position = expected
            .iter()
            .position(|id| *id == "types.BioMmcifWriteParams")
            .expect("the existing types.BioMmcifWriteParams fixture is present");
        expected.splice(
            position..position,
            ["BioCrystalInfo.space_group_number", "BioTransform.approx"],
        );
    }
    if cfg!(feature = "cap-fingerprints") {
        let position = expected
            .iter()
            .position(|id| *id == "types.TopologicalTorsionParams")
            .expect("the existing types.TopologicalTorsionParams fixture is present");
        expected.splice(
            position..position,
            [
                "types.AtomCodeExplanation",
                "errors.AtomCodeExplanationError",
            ],
        );
    }
    let position = expected
        .iter()
        .position(|id| *id == "types.PropertyValue")
        .expect("the existing types.PropertyValue fixture is present");
    expected.splice(
        position..position,
        [
            "types.MoleculeProperties",
            "types.SdfPropertyList",
            "types.SdfPropertyListTarget",
            "MoleculeProperties.name",
            "MoleculeProperties.sdf_data_fields",
            "MoleculeProperties.sdf_property_lists",
            "MoleculeProperties.props",
            "MoleculeProperties.prop",
            "MoleculeProperties.is_prop_computed",
            "MoleculeProperties.computed_prop_names",
            "SdfPropertyList.target",
            "SdfPropertyList.name",
            "SdfPropertyList.values",
        ],
    );
    let position = expected
        .iter()
        .position(|id| *id == "types.BondOrder")
        .expect("the existing types.BondOrder fixture is present");
    expected.splice(position..position, ["types.Hybridization"]);
    let position = expected
        .iter()
        .position(|id| *id == "types.BondOrder")
        .expect("the existing types.BondOrder fixture is present");
    expected.splice(position..position, ["Atom.hybridization"]);
    if cfg!(feature = "cap-fingerprints") {
        expected.extend([
            "Molecule.topological_torsion_ids",
            "Molecule.topological_torsion_ids_with_params",
            "AtomCodeExplanation.from_code",
            "AtomCodeExplanation.symbol",
            "AtomCodeExplanation.branch_count",
            "AtomCodeExplanation.pi_electrons",
            "AtomCodeExplanation.chirality",
        ]);
    }
    if cfg!(feature = "cap-fingerprints") {
        let position = expected
            .iter()
            .position(|id| *id == "types.MorganParams")
            .expect("the existing types.MorganParams fixture is present");
        expected.splice(
            position..position,
            [
                "types.FingerprintJsonError",
                "MorganParams.info_string",
                "MorganParams.to_json",
                "MorganParams.with_json",
                "AtomPairParams.info_string",
                "AtomPairParams.to_json",
                "AtomPairParams.with_json",
                "TopologicalTorsionParams.info_string",
                "TopologicalTorsionParams.to_json",
                "TopologicalTorsionParams.with_json",
            ],
        );
    }
    if cfg!(feature = "cap-fingerprints") {
        let position = expected
            .iter()
            .position(|id| *id == "types.FingerprintPreparationError")
            .expect("the existing fingerprint preparation error fixture is present");
        expected.splice(
            position..position,
            [
                "types.MorganAtomInvariantsGenerator",
                "types.MorganBondInvariantsGenerator",
                "types.MorganFingerprintGenerator",
                "types.MorganSettings",
                "types.MorganCallParams",
                "MorganAtomInvariantsGenerator.connectivity",
                "MorganAtomInvariantsGenerator.features",
                "MorganAtomInvariantsGenerator.atom_pair",
                "MorganBondInvariantsGenerator.new",
                "MorganBondInvariantsGenerator.use_bond_types",
                "MorganBondInvariantsGenerator.include_chirality",
                "MorganCallParams.new",
                "MorganFingerprintGenerator.new",
                "MorganFingerprintGenerator.from_json",
                "MorganFingerprintGenerator.settings",
                "MorganFingerprintGenerator.info_string",
                "MorganFingerprintGenerator.to_json",
                "MorganFingerprintGenerator.fingerprints",
                "Molecule.morgan_fingerprint_with_generator",
                "MorganFingerprintGenerator.counts",
                "Molecule.morgan_count_fingerprint_with_generator",
                "MorganFingerprintGenerator.sparse_fingerprints",
                "Molecule.morgan_sparse_fingerprint_with_generator",
                "MorganFingerprintGenerator.sparse_counts",
                "Molecule.morgan_sparse_count_fingerprint_with_generator",
                "MorganSettings.radius",
                "MorganSettings.set_radius",
                "MorganSettings.only_nonzero_invariants",
                "MorganSettings.set_only_nonzero_invariants",
                "MorganSettings.include_redundant_environments",
                "MorganSettings.set_include_redundant_environments",
                "MorganSettings.include_chirality",
                "MorganSettings.set_include_chirality",
                "MorganSettings.count_simulation",
                "MorganSettings.set_count_simulation",
                "MorganSettings.fp_size",
                "MorganSettings.set_fp_size",
                "MorganSettings.bits_per_feature",
                "MorganSettings.set_bits_per_feature",
                "MorganSettings.count_bounds",
                "MorganSettings.set_count_bounds",
                "MorganSettings.params",
            ],
        );
    }
    if cfg!(feature = "cap-serialization") {
        expected.extend([
            "types.PickleError",
            "Molecule.to_binary",
            "Molecule.from_binary",
        ]);
    }
    if cfg!(feature = "cap-hashing") {
        expected.extend([
            "types.CipRankError",
            "types.MoleculeHashError",
            "Molecule.molecular_hash",
            "Molecule.molecular_hash_with_ranks",
        ]);
    }
    // Element metadata and QueryGraph are the canonical stable prefix.
    // Move only these literal fixtures; retain all other original relative order.
    let mut prefix = vec![
        "types.Element",
        "Element.from_atomic_number",
        "Element.from_symbol",
        "Element.atomic_number",
        "Element.symbol",
        "types.ElementInfo",
    ];
    if cfg!(feature = "cap-valence") {
        prefix.push("module.element_info");
    }
    prefix.push("types.QueryGraph");
    for id in &prefix {
        let position = expected
            .iter()
            .position(|existing| existing == id)
            .expect("original canonical prefix fixture is present");
        expected.remove(position);
    }
    if cfg!(feature = "cap-fingerprints") {
        prefix.extend([
            "errors.TopologicalTorsionPathScoreError",
            "Molecule.topological_torsion_path_score",
            "explain_path_score",
            "types.AtomPairsParameters",
            "AtomPairsParameters.version",
            "AtomPairsParameters.num_type_bits",
            "AtomPairsParameters.num_pi_bits",
            "AtomPairsParameters.num_branch_bits",
            "AtomPairsParameters.num_chiral_bits",
            "AtomPairsParameters.code_size",
            "AtomPairsParameters.num_path_bits",
            "AtomPairsParameters.max_path_length",
            "AtomPairsParameters.num_atom_pair_fingerprint_bits",
            "AtomPairsParameters.atom_types",
        ]);
    }
    expected.splice(0..0, prefix);
    expected.extend(["types.LigandRef", "types.TetrahedralStereo"]);
    if cfg!(feature = "cap-stereo") {
        expected.extend([
            "types.StereoReadError",
            "Molecule.tetrahedral_stereo",
            "Molecule.perceive_stereochemistry",
            "Molecule.find_chiral_centers",
        ]);
    }
    if cfg!(feature = "cap-alignment") {
        expected.extend([
            "types.AlignmentAtomMap",
            "types.AlignmentParameters",
            "types.BestAlignmentParameters",
            "types.CoordinateRmsdParameters",
            "types.AllConformerRmsdParameters",
            "types.ConformerAlignmentParameters",
            "types.AlignmentResult",
            "types.AlignmentTransform",
            "types.ConformerRmsd",
            "types.ConformerAlignmentReport",
            "types.AlignmentError",
            "Molecule.alignment_transform_to",
            "Molecule.alignment_transform_to_with_params",
            "Molecule.best_alignment_to",
            "Molecule.best_alignment_to_with_params",
            "Molecule.best_rmsd_to",
            "Molecule.best_rmsd_to_with_params",
            "Molecule.coordinate_rmsd_to",
            "Molecule.coordinate_rmsd_to_with_params",
            "Molecule.all_conformer_best_rmsds",
            "Molecule.all_conformer_best_rmsds_with_params",
            "Molecule.with_alignment_to",
            "Molecule.with_alignment_to_with_params",
            "Molecule.align_to_",
            "Molecule.align_to_with_params_",
            "Molecule.with_aligned_conformers",
            "Molecule.with_aligned_conformers_with_params",
            "Molecule.align_conformers_",
            "Molecule.align_conformers_with_params_",
            "AlignmentResult.rmsd",
            "AlignmentResult.transform",
            "AlignmentResult.atom_map",
            "AlignmentTransform.matrix",
            "ConformerRmsd.rmsd",
            "ConformerRmsd.probe_conformer_id",
            "ConformerRmsd.reference_conformer_id",
            "ConformerAlignmentReport.rmsds",
            "AlignmentAtomMap.new",
            "AlignmentParameters.new",
            "BestAlignmentParameters.new",
            "CoordinateRmsdParameters.new",
            "AllConformerRmsdParameters.new",
            "ConformerAlignmentParameters.new",
        ]);
    }
    // Exact accepted MAIN fingerprint tail and new coordinate/conformer declarations.
    if cfg!(feature = "cap-fingerprints") {
        expected.extend([
            "types.MaccsFingerprintParams",
            "types.MaccsFingerprintError",
            "Molecule.maccs_fingerprint",
            "Molecule.maccs_fingerprint_raw",
            "Molecule.maccs_fingerprint_with_params",
            "types.LayeredFingerprintParams",
            "types.LayeredFingerprintLayers",
            "types.LayeredFingerprintResult",
            "types.LayeredFingerprintError",
            "Molecule.layered_fingerprint",
            "Molecule.layered_fingerprint_with_params",
            "Molecule.layered_fingerprint_with_output",
            "Molecule.layered_fingerprint_with_output_with_params",
            "layered_query_fingerprint_with_params",
            "layered_query_fingerprint_with_output_with_params",
            "LayeredFingerprintResult.fingerprint",
            "LayeredFingerprintResult.atom_counts",
            "LayeredFingerprintLayers.bits",
            "LayeredFingerprintLayers.from_bits_retain",
            "Fingerprint.from_on_bits",
            "Fingerprint.tanimoto",
        ]);
    }
    if cfg!(feature = "cap-transforms") {
        expected.extend([
            "types.CoordinateZPolicy",
            "types.Coordinate2DInputParams",
            "types.Coordinate3DInputParams",
            "types.Replace3DCoordinatesParams",
            "types.CoordinateInputError",
            "types.Coordinate3DReadError",
            "Molecule.coordinates_3d",
            "CoordinateZPolicy.from_name",
            "Molecule.with_2d_coordinate_block",
            "Molecule.with_2d_coordinate_block_with_params",
            "Molecule.set_2d_coordinates_",
            "Molecule.set_2d_coordinates_with_params_",
            "Molecule.with_3d_coordinates",
            "Molecule.with_3d_coordinates_with_params",
            "Molecule.set_3d_coordinates_",
            "Molecule.set_3d_coordinates_with_params_",
            "Molecule.with_added_3d_conformer",
            "Molecule.with_added_3d_conformer_with_params",
            "Molecule.add_3d_conformer_",
            "Molecule.add_3d_conformer_with_params_",
            "Molecule.with_only_3d_conformer",
            "Molecule.with_only_3d_conformer_with_params",
            "Molecule.set_only_3d_conformer_",
            "Molecule.set_only_3d_conformer_with_params_",
            "Molecule.with_cleared_3d_conformers",
            "Molecule.clear_3d_conformers_",
        ]);
    }
    if cfg!(feature = "cap-conformer") {
        expected.extend([
            "types.EmbedMoleculeResult",
            "EmbedMoleculeResult.molecule",
            "EmbedMoleculeResult.params",
            "types.EmbedMultipleConfsResult",
            "EmbedMultipleConfsResult.molecule",
            "EmbedMultipleConfsResult.params",
            "EmbedMoleculeResult.conf_id",
            "EmbedMoleculeResult.ok",
            "EmbedMultipleConfsResult.conf_ids",
            "EmbedMultipleConfsResult.requested_num_confs",
            "EmbedMultipleConfsResult.generated_count",
            "Molecule.with_3d_conformer",
            "Molecule.with_3d_conformer_with_params",
            "Molecule.embed_3d_conformer_",
            "Molecule.embed_3d_conformer_with_params_",
            "Molecule.with_3d_conformer_result",
            "Molecule.with_3d_conformer_result_with_params",
            "Molecule.embed_3d_conformer_result_",
            "Molecule.embed_3d_conformer_result_with_params_",
            "Molecule.with_3d_conformers",
            "Molecule.with_3d_conformers_with_params",
            "Molecule.embed_3d_conformers_",
            "Molecule.embed_3d_conformers_with_params_",
            "Molecule.with_3d_conformers_result",
            "Molecule.with_3d_conformers_result_with_params",
            "Molecule.embed_3d_conformers_result_",
            "Molecule.embed_3d_conformers_result_with_params_",
            "types.EmbedParams",
            "EmbedParams.max_iterations",
            "EmbedParams.num_threads",
            "EmbedParams.random_seed",
            "EmbedParams.clear_confs",
            "EmbedParams.use_random_coords",
            "EmbedParams.box_size_mult",
            "EmbedParams.rand_neg_eig",
            "EmbedParams.num_zero_fail",
            "EmbedParams.coord_map",
            "EmbedParams.optimizer_force_tol",
            "EmbedParams.ignore_smoothing_failures",
            "EmbedParams.enforce_chirality",
            "EmbedParams.use_exp_torsion_angle_prefs",
            "EmbedParams.use_basic_knowledge",
            "EmbedParams.verbose",
            "EmbedParams.basin_thresh",
            "EmbedParams.prune_rms_thresh",
            "EmbedParams.only_heavy_atoms_for_rms",
            "EmbedParams.et_version",
            "EmbedParams.embed_fragments_separately",
            "EmbedParams.use_small_ring_torsions",
            "EmbedParams.use_macrocycle_torsions",
            "EmbedParams.use_macrocycle14config",
            "EmbedParams.timeout",
            "EmbedParams.cpci",
            "EmbedParams.force_trans_amides",
            "EmbedParams.use_symmetry_for_pruning",
            "EmbedParams.bounds_mat_force_scaling",
            "EmbedParams.track_failures",
            "EmbedParams.failures",
            "EmbedParams.enable_sequential_random_seeds",
            "EmbedParams.symmetrize_conjugated_terminal_groups_for_pruning",
            "EmbedParams.new",
            "EmbedParams.dg",
            "EmbedParams.kdg",
            "EmbedParams.etdg",
            "EmbedParams.etdg_v2",
            "EmbedParams.etkdg",
            "EmbedParams.etkdg_v2",
            "EmbedParams.etkdg_v3",
            "EmbedParams.sr_etkdg_v3",
            "EmbedParams.to_json",
            "EmbedParams.with_json",
            "Molecule.num_3d_conformers",
            "Molecule.dg_bounds_matrix",
        ]);
    }
    #[cfg(feature = "cap-fingerprints")]
    {
        expected.extend([
            "types.PatternFingerprintParams",
            "types.PatternFingerprintError",
            "Molecule.pattern_fingerprint",
            "Molecule.pattern_fingerprint_with_params",
            "pattern_query_fingerprint",
            "pattern_query_fingerprint_with_params",
            "types.TopologicalFingerprintParams",
            "types.TopologicalFingerprintOutputRequest",
            "types.TopologicalFingerprintOutput",
            "types.TopologicalFingerprintResult",
            "types.TopologicalFingerprintError",
            "Molecule.topological_fingerprint",
            "Molecule.topological_fingerprint_with_params",
            "Molecule.topological_fingerprint_with_output",
            "Molecule.topological_fingerprint_with_output_with_params",
            "topological_query_fingerprint_with_params",
            "topological_query_fingerprint_with_output_with_params",
            "TopologicalFingerprintResult.fingerprint",
            "TopologicalFingerprintResult.atom_bits",
            "TopologicalFingerprintResult.bit_info",
        ]);
    }
    // ROOT-authorized BATCH registry fixture; fixed declaration order and feature gates.
    // Delivery proposal: exact pinned current IO prerequisite declaration order.
    if cfg!(feature = "cap-io") {
        expected.extend([
            "SdfRecord.from_query_graph",
            "SdfRecord.to_mol",
            "SdfRecord.to_mol_with_params",
            "SdfRecord.to_sdf",
            "SdfRecord.to_sdf_with_params",
            "types.PropertyStringError",
            "PropertyStringError.kind",
            "Molecule.atom_property_string",
            "Molecule.bond_property_string",
            "types.MolecularIoError",
            "types.Mol2ReadParams",
            "types.Mol2Type",
            "types.Mol2ReadError",
            "types.Mol2PostError",
            "types.XyzReadError",
            "types.XyzWriteError",
            "types.XyzWriteParams",
            "types.SdfDataset",
            "types.SdfDatasetIterator",
            "types.SdfRecordMetadata",
            "types.SdfRecordStream",
            "Molecule.from_xyz_block",
            "Molecule.read_xyz",
            "Molecule.from_mol2",
            "Molecule.read_mol2",
            "Molecule.read_mol",
            "Molecule.read_sdf",
            "Molecule.from_mol2_with_params",
            "Molecule.read_mol2_with_params",
            "Molecule.read_mol_with_params",
            "Molecule.read_sdf_with_params",
            "Molecule.from_mol",
            "Molecule.from_mol_with_params",
            "Molecule.to_xyz",
            "Molecule.to_xyz_with_params",
            "Molecule.write_xyz",
            "Molecule.write_xyz_with_params",
            "SdfDataset.open",
            "SdfDataset.open_with_params",
            "SdfRecordStream.open",
            "SdfRecordStream.open_with_params",
            "SdfDataset.len",
            "SdfDataset.is_empty",
            "SdfDataset.path",
            "SdfDataset.iter",
            "SdfDataset.metadata",
            "SdfDataset.record",
            "SdfDataset.record_with_params",
            "SdfDataset.record_text",
            "SdfRecordStream.next_record",
            "SdfRecordStream.is_end",
            "SdfRecordStream.records_consumed",
            "SdfRecordStream.bytes_consumed",
            "SdfRecordStream.lines_consumed",
            "SdfRecordMetadata.index",
            "SdfRecordMetadata.byte_offset",
            "SdfRecordMetadata.byte_len",
            "SdfRecordMetadata.byte_range",
            "SdfRecordMetadata.line_range",
            "SdfRecordMetadata.title",
            "SdfRecord.title",
            "SdfRecord.index",
            "SdfRecord.data_field",
        ]);
    }
    if cfg!(feature = "cap-io") && cfg!(feature = "cap-batch") {
        expected.extend([
            "types.SdfBatchIterator",
            "types.SdfReaderBatchIterator",
            "MoleculeBatch.from_sdf_records",
            "MoleculeBatch.from_sdf_records_with_params",
            "MoleculeBatch.read_sdf",
            "MoleculeBatch.read_sdf_with_params",
            "MoleculeBatch.from_dataset_indices",
            "MoleculeBatch.get",
            "MoleculeBatch.records",
            "SdfDataset.batches",
            "SdfBatchIterator.next_batch",
            "SdfReaderBatchIterator.next_batch",
        ]);
    }
    if cfg!(feature = "cap-io") {
        expected.extend([
            "types.SdfReader",
            "SdfReader.open",
            "SdfReader.open_with_params",
            "SdfReader.path",
            "SdfReader.params",
        ]);
    }
    if cfg!(feature = "cap-io") && cfg!(feature = "cap-batch") {
        expected.extend(["SdfReader.batches"]);
    }
    if cfg!(feature = "cap-io") {
        expected.extend([
            "types.MolWriteError",
            "types.MolBlockWriteParams",
            "types.MolCoordinateSelection",
            "types.SdfFormat",
            "Molecule.to_mol",
            "Molecule.to_mol_with_params",
            "Molecule.to_sdf",
            "Molecule.to_sdf_with_params",
            "Molecule.to_sdf_2d",
            "Molecule.to_sdf_2d_with_params",
            "Molecule.to_sdf_3d",
            "Molecule.to_sdf_3d_with_params",
            "Molecule.write_mol",
            "Molecule.write_mol_with_params",
            "Molecule.write_sdf",
            "Molecule.write_sdf_with_params",
            "Molecule.write_sdf_files",
            "Molecule.write_sdf_files_with_params",
        ]);
    }
    if cfg!(feature = "cap-io") && cfg!(feature = "cap-batch") {
        expected.extend([
            "types.BatchExportParams",
            "MoleculeBatch.to_sdf",
            "MoleculeBatch.to_sdf_files",
            "MoleculeBatch.to_sdf_with_params",
            "MoleculeBatch.to_sdf_files_with_params",
            "BatchExportReport.total",
            "BatchExportReport.success",
            "BatchExportReport.failed",
            "BatchExportReport.errors",
            "SdfRecordStream.batches",
        ]);
    }
    if cfg!(feature = "cap-batch") {
        expected.extend([
            "types.MoleculeBatch",
            "types.BatchRecord",
            "types.BatchError",
            "types.BatchErrorMode",
            "types.BatchValidationError",
            "types.BatchParams",
            "types.BatchExportReport",
            "MoleculeBatch.len",
            "MoleculeBatch.is_empty",
            "MoleculeBatch.valid_mask",
            "MoleculeBatch.invalid_mask",
            "MoleculeBatch.valid_count",
            "MoleculeBatch.invalid_count",
            "MoleculeBatch.errors",
            "MoleculeBatch.parallel_jobs",
            "MoleculeBatch.progress_bar",
            "MoleculeBatch.to_list",
            "MoleculeBatch.with_valid_records",
            "MoleculeBatch.with_parallel_jobs",
            "MoleculeBatch.with_progress_bar",
            "MoleculeBatch.from_records",
        ]);
        if cfg!(feature = "cap-smiles") {
            expected.extend([
                "MoleculeBatch.from_smiles_list",
                "MoleculeBatch.from_smiles_list_with_params",
            ]);
        }
        if cfg!(feature = "cap-sanitize") {
            expected.extend([
                "MoleculeBatch.sanitize",
                "MoleculeBatch.sanitize_with_params",
            ]);
        }
        if cfg!(feature = "cap-hydrogens") {
            expected.extend([
                "MoleculeBatch.with_hydrogens",
                "MoleculeBatch.with_hydrogens_with_params",
                "MoleculeBatch.without_hydrogens",
                "MoleculeBatch.without_hydrogens_with_params",
            ]);
        }
        if cfg!(feature = "cap-kekulize") {
            expected.extend([
                "MoleculeBatch.with_kekulized_bonds",
                "MoleculeBatch.with_kekulized_bonds_with_params",
            ]);
        }
        if cfg!(feature = "cap-depict") {
            expected.extend([
                "MoleculeBatch.with_2d_coordinates",
                "MoleculeBatch.with_2d_coordinates_with_params",
            ]);
        }
        expected.extend(["types.BatchQueryParams"]);
        if cfg!(feature = "cap-smiles") {
            expected.extend([
                "MoleculeBatch.to_smiles_list",
                "MoleculeBatch.to_smiles_list_with_params",
            ]);
        }
        if cfg!(feature = "cap-conformer") {
            expected.extend([
                "MoleculeBatch.dg_bounds_matrix_list",
                "MoleculeBatch.dg_bounds_matrix_list_with_params",
            ]);
        }
        if cfg!(feature = "cap-depict") {
            expected.extend([
                "MoleculeBatch.to_svg_list",
                "MoleculeBatch.to_svg_list_with_params",
                "types.BatchImageParams",
                "types.BatchImageError",
                "MoleculeBatch.to_images",
                "MoleculeBatch.to_images_with_params",
            ]);
        }
        if cfg!(feature = "cap-fingerprints") {
            expected.extend([
                "MoleculeBatch.fingerprint_atom_pair_list",
                "MoleculeBatch.fingerprint_atom_pair_list_with_params",
                "MoleculeBatch.fingerprint_atom_pair_sparse_count_list",
                "MoleculeBatch.fingerprint_atom_pair_sparse_count_list_with_params",
                "MoleculeBatch.fingerprint_atom_pair_count_list",
                "MoleculeBatch.fingerprint_atom_pair_count_list_with_params",
                "MoleculeBatch.fingerprint_atom_pair_sparse_bits_list",
                "MoleculeBatch.fingerprint_atom_pair_sparse_bits_list_with_params",
                "MoleculeBatch.fingerprint_layered_list",
                "MoleculeBatch.fingerprint_layered_list_with_params",
                "MoleculeBatch.fingerprint_layered_with_output_list",
                "MoleculeBatch.fingerprint_layered_with_output_list_with_params",
                "MoleculeBatch.pattern_fingerprint_list",
                "MoleculeBatch.pattern_fingerprint_list_with_params",
                "MoleculeBatch.fingerprint_morgan_list",
                "MoleculeBatch.fingerprint_morgan_list_with_params",
                "types.BatchFingerprintOutput",
                "MoleculeBatch.fingerprint_atom_pair_with_output_list",
                "MoleculeBatch.fingerprint_atom_pair_with_output_list_with_params",
                "MoleculeBatch.fingerprint_morgan_with_output_list",
                "MoleculeBatch.fingerprint_morgan_with_output_list_with_params",
                "MoleculeBatch.fingerprint_morgan_list_with_generator_params",
                "MoleculeBatch.fingerprint_morgan_with_output_list_with_generator_params",
                "types.BatchFingerprintAdditionalOutput",
                "types.BatchFingerprintOutputError",
                "BatchFingerprintAdditionalOutput.atom_counts",
                "BatchFingerprintAdditionalOutput.atom_to_bits",
                "BatchFingerprintAdditionalOutput.bit_info_map",
                "BatchFingerprintAdditionalOutput.bit_paths",
                "BatchFingerprintAdditionalOutput.atoms_per_bit",
                "BatchFingerprintOutput.fingerprint",
                "BatchFingerprintOutput.additional_output",
            ]);
        }
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

#[cfg(feature = "cap-forcefields")]
#[test]
fn uff_public_registry_has_exact_value_result_contract() {
    for (id, role) in [
        ("types.UffOptimizationParams", BindingTypeRole::Parameter),
        ("types.UffOptimizationResult", BindingTypeRole::Result),
        ("types.UffOptimizationError", BindingTypeRole::Error),
        ("types.UffOptimizationErrorKind", BindingTypeRole::Value),
    ] {
        let row = entry(id);
        assert_eq!(row.status, FunctionStatus::Experimental);
        assert_eq!(row.feature, "cap-forcefields");
        assert_eq!(row.type_role, Some(role));
    }
    for (method, javascript, parameter_count) in [
        ("with_uff_optimized", "withUffOptimized", 0),
        (
            "with_uff_optimized_with_params",
            "withUffOptimizedWithParams",
            1,
        ),
    ] {
        let row = entry(&format!("Molecule.{method}"));
        assert_eq!(row.status, FunctionStatus::Experimental);
        assert_eq!(row.python_name, method);
        assert_eq!(row.javascript_name, javascript);
        let callable = row.callable.unwrap();
        assert_eq!(callable.kind, BindingKind::Instance);
        assert_eq!(callable.state_model, StateModel::ValueReturning);
        assert_eq!(callable.parameters.len(), parameter_count);
        assert_eq!(
            callable.output_type.replace(' ', ""),
            "crate::UffOptimizationResult"
        );
        assert_eq!(
            callable.error_type.map(|name| name.replace(' ', "")),
            Some("crate::OperationError".into())
        );
        assert_eq!(callable.operation_semantic_id, Some(method));
    }
    let _: fn(
        &cosmolkit::Molecule,
    ) -> Result<cosmolkit::UffOptimizationResult, cosmolkit::OperationError> =
        cosmolkit::Molecule::with_uff_optimized;
    let _: fn(
        &cosmolkit::Molecule,
        &cosmolkit::UffOptimizationParams,
    ) -> Result<cosmolkit::UffOptimizationResult, cosmolkit::OperationError> =
        cosmolkit::Molecule::with_uff_optimized_with_params;
    let _: fn(&cosmolkit::UffOptimizationError) -> cosmolkit::UffOptimizationErrorKind =
        cosmolkit::UffOptimizationError::kind;
    let params = cosmolkit::UffOptimizationParams::default();
    assert_eq!(params.max_iterations, 1000);
    assert_eq!(params.vdw_threshold.to_bits(), 10.0_f64.to_bits());
    assert!(params.ignore_interfragment_interactions);
    assert_eq!(params.conformer_id, None);
}

#[cfg(feature = "cap-forcefields")]
#[test]
fn uff_param_p10_registry_has_exact_read_only_query_contract() {
    let forcefields_ids = BINDING_CONTRACT
        .iter()
        .filter(|row| row.feature == "cap-forcefields")
        .map(|row| row.semantic_id)
        .collect::<Vec<_>>();
    assert_eq!(
        forcefields_ids,
        vec![
            "types.UffEvaluationParams",
            "types.UffEnergyGradient",
            "Molecule.uff_energy_gradient",
            "Molecule.uff_energy_gradient_with_params",
            "UffEnergyGradient.energy",
            "UffEnergyGradient.gradient",
            "UffOptimizationResult.molecule",
            "UffOptimizationResult.status_code",
            "UffOptimizationResult.needs_more",
            "UffOptimizationResult.energy",
            "UffConformerResult.conformer_id",
            "UffConformerResult.status_code",
            "UffConformerResult.needs_more",
            "UffConformerResult.energy",
            "UffConformerOptimizationResult.molecule",
            "UffConformerOptimizationResult.conformer_results",
            "MmffAtomProperties.atom_type",
            "MmffAtomProperties.formal_charge",
            "MmffAtomProperties.partial_charge",
            "types.MmffEvaluationParams",
            "types.MmffEnergyGradient",
            "Molecule.mmff_energy_gradient",
            "Molecule.mmff_energy_gradient_with_params",
            "MmffEnergyGradient.energy",
            "MmffEnergyGradient.gradient",
            "MmffOptimizeMoleculeResult.molecule",
            "MmffOptimizeMoleculeResult.needs_more",
            "MmffOptimizeMoleculeResult.status_code",
            "MmffOptimizeMoleculeConfsResult.molecule",
            "MmffOptimizeMoleculeConfsResult.conformer_results",
            "MmffOptimizeMoleculeConfResult.needs_more",
            "MmffOptimizeMoleculeConfResult.status_code",
            "MmffOptimizeMoleculeConfResult.energy",
            "MmffProperties.is_valid",
            "MmffProperties.variant",
            "MmffProperties.atoms",
            "MmffProperties.atom_type",
            "MmffProperties.formal_charge",
            "MmffProperties.partial_charge",
            "types.MmffOptimizationParams",
            "types.MmffConformerOptimizationParams",
            "types.MmffOptimizeMoleculeResult",
            "types.MmffOptimizeMoleculeConfResult",
            "types.MmffOptimizeMoleculeConfsResult",
            "types.MmffOptimizationError",
            "Molecule.with_mmff_optimized",
            "Molecule.with_mmff_optimized_with_params",
            "Molecule.with_mmff_optimized_confs",
            "Molecule.with_mmff_optimized_confs_with_params",
            "types.UffParameterQueryError",
            "types.UffParameterError",
            "types.UffParameterErrorKind",
            "UffParameterError.kind",
            "types.MmffProperties",
            "types.MmffPropertiesParams",
            "types.MmffAtomProperties",
            "types.MmffVariant",
            "types.MmffMolPropertiesError",
            "Molecule.mmff_has_all_molecule_params",
            "Molecule.mmff_properties",
            "Molecule.mmff_properties_with_params",
            "Molecule.uff_has_all_molecule_params",
            "types.UffOptimizationParams",
            "types.UffOptimizationResult",
            "types.UffOptimizationError",
            "types.UffOptimizationErrorKind",
            "types.UffConformerOptimizationParams",
            "types.UffConformerOptimizationResult",
            "types.UffConformerResult",
            "UffOptimizationError.kind",
            "Molecule.with_uff_optimized",
            "Molecule.with_uff_optimized_with_params",
            "Molecule.with_uff_optimized_confs",
            "Molecule.with_uff_optimized_confs_with_params",
        ],
        "complete forcefields exposes the exact UFF and MMFF identities"
    );

    for (semantic_id, rust_name, python_name, javascript_name, role) in [
        (
            "types.UffParameterQueryError",
            "UffParameterQueryError",
            "UffParameterQueryError",
            "UffParameterQueryError",
            BindingTypeRole::Error,
        ),
        (
            "types.UffParameterError",
            "UffParameterError",
            "UffParameterError",
            "UffParameterError",
            BindingTypeRole::Error,
        ),
        (
            "types.UffParameterErrorKind",
            "UffParameterErrorKind",
            "UffParameterErrorKind",
            "UffParameterErrorKind",
            BindingTypeRole::Value,
        ),
    ] {
        let row = entry(semantic_id);
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.owner, BindingOwner::Type);
        assert!(row.rust_path.ends_with(rust_name));
        assert_eq!(row.python_name, python_name);
        assert_eq!(row.javascript_name, javascript_name);
        assert_eq!(row.feature, "cap-forcefields");
        assert_eq!(row.status, FunctionStatus::Experimental);
        assert_eq!(row.callable, None);
        assert_eq!(row.type_role, Some(role));
    }

    let kind = entry("UffParameterError.kind");
    assert_eq!(kind.item, BindingItem::Callable);
    assert_eq!(kind.owner, BindingOwner::Type);
    assert_eq!(
        kind.rust_path.replace(' ', ""),
        "crate::UffParameterError::kind"
    );
    assert_eq!(kind.python_name, "kind");
    assert_eq!(kind.javascript_name, "kind");
    assert_eq!(kind.feature, "cap-forcefields");
    assert_eq!(kind.status, FunctionStatus::Experimental);
    assert_eq!(kind.type_role, None);
    let kind_callable = kind.callable.expect("kind accessor callable metadata");
    assert_eq!(kind_callable.kind, BindingKind::Instance);
    assert!(kind_callable.parameters.is_empty());
    assert_eq!(
        kind_callable.output_type.replace(' ', ""),
        "crate::UffParameterErrorKind"
    );
    assert_eq!(kind_callable.error_type, None);
    assert_eq!(kind_callable.state_model, StateModel::ReadOnly);
    assert_eq!(kind_callable.operation_semantic_id, None);
    let _kind_signature: fn(&cosmolkit::UffParameterError) -> cosmolkit::UffParameterErrorKind =
        cosmolkit::UffParameterError::kind;
    let kind_name = |value| match value {
        cosmolkit::UffParameterErrorKind::Preparation => "Preparation",
        cosmolkit::UffParameterErrorKind::ParameterTable => "ParameterTable",
        cosmolkit::UffParameterErrorKind::Typing => "Typing",
    };
    assert_eq!(
        [
            kind_name(cosmolkit::UffParameterErrorKind::Preparation),
            kind_name(cosmolkit::UffParameterErrorKind::ParameterTable),
            kind_name(cosmolkit::UffParameterErrorKind::Typing),
        ],
        ["Preparation", "ParameterTable", "Typing"]
    );

    let query = entry("Molecule.uff_has_all_molecule_params");
    assert_eq!(query.item, BindingItem::Callable);
    assert_eq!(query.owner, BindingOwner::Molecule);
    assert_eq!(
        query.rust_path.replace(' ', ""),
        "crate::Molecule::uff_has_all_molecule_params"
    );
    assert_eq!(query.python_name, "uff_has_all_molecule_params");
    assert_eq!(query.javascript_name, "uffHasAllMoleculeParams");
    assert_eq!(query.feature, "cap-forcefields");
    assert_eq!(query.status, FunctionStatus::Experimental);
    assert_eq!(query.type_role, None);
    let query_callable = query.callable.expect("query callable metadata");
    assert_eq!(query_callable.kind, BindingKind::Instance);
    assert!(query_callable.parameters.is_empty());
    assert_eq!(query_callable.output_type, "bool");
    assert_eq!(
        query_callable.error_type.map(|name| name.replace(' ', "")),
        Some("crate::UffParameterQueryError".to_owned())
    );
    assert_eq!(query_callable.state_model, StateModel::ReadOnly);
    assert_eq!(query_callable.operation_semantic_id, None);
    let _query_signature: fn(
        &cosmolkit::Molecule,
    ) -> Result<bool, cosmolkit::UffParameterQueryError> =
        cosmolkit::Molecule::uff_has_all_molecule_params;
}

#[cfg(feature = "cap-forcefields")]
#[test]
fn uff_param_p10_external_feature_isolation() {
    let query_available = uff_param_p10_compile_check("query_available", Some("forcefields"));
    assert!(
        query_available.status.success(),
        "forcefields-only public query consumer failed to compile:\n{}",
        String::from_utf8_lossy(&query_available.stderr)
    );

    let query_capability_only =
        uff_param_p10_compile_check("query_available", Some("cap-forcefields"));
    assert!(
        query_capability_only.status.success(),
        "cap-forcefields-only public query consumer failed to compile:\n{}",
        String::from_utf8_lossy(&query_capability_only.stderr)
    );

    let base_without_forcefields = uff_param_p10_compile_check("base_api_available", None);
    assert!(
        base_without_forcefields.status.success(),
        "base Molecule positive control failed without forcefields:\n{}",
        String::from_utf8_lossy(&base_without_forcefields.stderr)
    );

    let query_without_forcefields = uff_param_p10_compile_check("query_unavailable", None);
    assert!(!query_without_forcefields.status.success());
    let query_errors = String::from_utf8_lossy(&query_without_forcefields.stderr);
    assert!(query_errors.contains("error[E0432]"), "{query_errors}");
    for absent_item in [
        "UffParameterError",
        "UffParameterErrorKind",
        "UffParameterQueryError",
    ] {
        assert!(
            query_errors.contains(absent_item),
            "missing {absent_item}: {query_errors}"
        );
    }
    assert!(
        query_errors.contains("no method named `uff_has_all_molecule_params`"),
        "query method absence was not the intended compile failure:\n{query_errors}"
    );

    let forcefields_only_exclusions =
        uff_param_p10_compile_check("forcefields_only_exclusions", Some("cap-forcefields"));
    assert!(!forcefields_only_exclusions.status.success());
    let excluded_errors = String::from_utf8_lossy(&forcefields_only_exclusions.stderr);
    assert!(
        excluded_errors.contains("error[E0432]"),
        "{excluded_errors}"
    );
    assert!(
        excluded_errors.contains("error[E0599]"),
        "{excluded_errors}"
    );
    assert!(
        excluded_errors.contains("error[E0603]"),
        "{excluded_errors}"
    );
    for absent_item in [
        "ValenceParams",
        "RingSearchParams",
        "ForceFieldError",
        "ForceFieldOptions",
        "mmff_has_all_molecule_params",
        "mmff_optimize",
        "uff_has_all_molecule_params",
        "with_assigned_valence",
        "with_assigned_rings",
        "forcefields",
    ] {
        assert!(
            excluded_errors.contains(absent_item),
            "missing {absent_item}: {excluded_errors}"
        );
    }
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
    assert_eq!(molecule.status, FunctionStatus::Experimental);
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
    assert_eq!(version.status, FunctionStatus::Experimental);
    assert_eq!(version.type_role, None);
    let callable = version.callable.expect("version callable metadata");
    assert_eq!(callable.kind, BindingKind::Module);
    assert!(callable.parameters.is_empty());
    assert_eq!(callable.output_type.replace(' ', ""), "&'staticstr");
    assert_eq!(callable.error_type, None);
    assert_eq!(callable.state_model, StateModel::ReadOnly);
    assert_eq!(callable.operation_semantic_id, None);
}

#[cfg(feature = "cap-depict")]
#[test]
fn drawing_entries_have_exact_experimental_shared_query_contracts() {
    let error = entry("types.DrawingError");
    assert_eq!(error.item, BindingItem::Type);
    assert_eq!(error.owner, BindingOwner::Type);
    assert_eq!(error.rust_path.replace(' ', ""), "crate::DrawingError");
    assert_eq!(error.python_name, "DrawingError");
    assert_eq!(error.javascript_name, "DrawingError");
    assert_eq!(error.feature, "cap-depict");
    assert_eq!(error.status, FunctionStatus::Experimental);
    assert_eq!(error.type_role, Some(BindingTypeRole::Error));
    assert_eq!(error.callable, None);
    for (name, javascript, output) in [
        ("to_svg", "toSvg", "String"),
        ("to_png", "toPng", "Vec<u8>"),
    ] {
        let row = entry(&format!("Molecule.{name}"));
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Molecule);
        assert_eq!(
            row.rust_path.replace(' ', ""),
            format!("crate::Molecule::{name}")
        );
        assert_eq!(row.python_name, name);
        assert_eq!(row.javascript_name, javascript);
        assert_eq!(row.feature, "cap-depict");
        assert_eq!(
            row.status,
            if name == "to_svg" {
                canonical_svg_identity_status()
            } else {
                FunctionStatus::Experimental
            }
        );
        let callable = row.callable.unwrap();
        assert_eq!(callable.kind, BindingKind::Instance);
        assert_eq!(callable.receiver, Some(cosmolkit::BindingReceiver::Shared));
        assert_eq!(callable.state_model, StateModel::ReadOnly);
        assert_eq!(callable.operation_semantic_id, None);
        assert_eq!(callable.output_type.replace(' ', ""), output);
        assert_eq!(
            callable.error_type.unwrap().replace(' ', ""),
            "crate::DrawingError"
        );
        assert_eq!(callable.parameters.len(), 2);
        for (parameter, expected_name) in callable.parameters.iter().zip(["width", "height"]) {
            assert_eq!(parameter.name, expected_name);
            assert_eq!(parameter.type_name, "u32");
            assert_eq!(parameter.default, BindingDefault::Required);
        }
    }
}

#[cfg(feature = "cap-descriptors")]
#[test]
fn descriptor_entries_are_exact_canonical_read_only_methods() {
    for (semantic_id, rust_name, javascript, output, error) in [
        (
            "Molecule.molecular_weight",
            "molecular_weight",
            "molecularWeight",
            "f64",
            "crate::DescriptorReadError",
        ),
        (
            "Molecule.exact_molecular_weight",
            "exact_molecular_weight",
            "exactMolecularWeight",
            "f64",
            "crate::DescriptorReadError",
        ),
        (
            "Molecule.molecular_formula",
            "molecular_formula",
            "molecularFormula",
            "String",
            "crate::DescriptorReadError",
        ),
        (
            "Molecule.num_heavy_atoms",
            "num_heavy_atoms",
            "numHeavyAtoms",
            "u32",
            "crate::DescriptorReadError",
        ),
        (
            "Molecule.total_atom_count",
            "total_atom_count",
            "totalAtomCount",
            "u32",
            "crate::DescriptorReadError",
        ),
        (
            "Molecule.lipinski_hba",
            "lipinski_hba",
            "lipinskiHba",
            "u32",
            "crate::DescriptorReadError",
        ),
        (
            "Molecule.lipinski_hbd",
            "lipinski_hbd",
            "lipinskiHbd",
            "u32",
            "crate::DescriptorReadError",
        ),
        (
            "Molecule.fraction_csp3",
            "fraction_csp3",
            "fractionCsp3",
            "f64",
            "crate::DescriptorReadError",
        ),
    ] {
        let entry = entry(semantic_id);
        assert_eq!(entry.item, BindingItem::Callable);
        assert_eq!(entry.owner, BindingOwner::Molecule);
        assert!(entry.rust_path.replace(' ', "").ends_with(rust_name));
        assert_eq!(entry.python_name, rust_name);
        assert_eq!(entry.javascript_name, javascript);
        assert_eq!(entry.feature, "cap-descriptors");
        assert_eq!(entry.status, FunctionStatus::Experimental);
        assert_eq!(entry.type_role, None);
        let callable = entry.callable.expect("descriptor callable metadata");
        assert_eq!(callable.kind, BindingKind::Instance);
        assert!(callable.parameters.is_empty());
        assert_eq!(callable.output_type, output);
        assert_eq!(
            callable.error_type.map(|name| name.replace(' ', "")),
            Some(error.replace(' ', ""))
        );
        assert_eq!(callable.state_model, StateModel::ReadOnly);
        assert_eq!(callable.operation_semantic_id, None);
    }
}

#[cfg(feature = "cap-descriptors")]
#[test]
fn descriptor_error_type_entries_carry_error_roles_and_exact_names() {
    for (semantic_id, rust_type, projection) in [
        (
            "types.DescriptorError",
            "DescriptorError",
            "DescriptorError",
        ),
        (
            "types.DescriptorReadError",
            "DescriptorReadError",
            "DescriptorReadError",
        ),
    ] {
        let entry = entry(semantic_id);
        assert_eq!(entry.item, BindingItem::Type);
        assert_eq!(entry.owner, BindingOwner::Type);
        assert!(entry.rust_path.replace(' ', "").ends_with(rust_type));
        assert_eq!(entry.python_name, projection);
        assert_eq!(entry.javascript_name, projection);
        assert_eq!(entry.feature, "cap-descriptors");
        assert_eq!(entry.status, FunctionStatus::Experimental);
        assert_eq!(entry.type_role, Some(BindingTypeRole::Error));
        assert!(entry.callable.is_none());
    }
}

#[cfg(feature = "cap-hydrogens")]
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
        assert_eq!(entry.feature, "cap-hydrogens");
        assert_eq!(entry.status, FunctionStatus::Experimental);
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
        receiver: None,
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
            FunctionStatus::Parity { reference: "RDKit" },
            FunctionStatus::ParityWithDifferences {
                reference: "Gemmi",
                explanation: "Approved difference for the documented boundary"
            },
            FunctionStatus::Native,
            FunctionStatus::Experimental
        ]
        .len(),
        4
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
    for old_field in ["exposure:", "support:", "parity:"] {
        assert!(!registry.contains(old_field));
    }
}

#[test]
fn status_commitments_are_per_function_and_shared_with_registered_operations() {
    let mut parity = Vec::new();
    for contract in BINDING_CONTRACT {
        let fuzzy = matches!(
            contract.semantic_id,
            "SparseCountFingerprint.fuzzy_and"
                | "SparseCountFingerprint.fuzzy_or"
                | "SparseCountFingerprint32.fuzzy_and"
                | "SparseCountFingerprint32.fuzzy_or"
        );
        let native_coordinates = matches!(
            contract.semantic_id,
            "types.CoordinateZPolicy"
                | "types.Coordinate2DInputParams"
                | "types.Coordinate3DInputParams"
                | "types.Replace3DCoordinatesParams"
                | "types.CoordinateInputError"
                | "types.Coordinate3DReadError"
                | "Molecule.coordinates_3d"
                | "CoordinateZPolicy.from_name"
                | "Molecule.with_2d_coordinate_block"
                | "Molecule.with_2d_coordinate_block_with_params"
                | "Molecule.set_2d_coordinates_"
                | "Molecule.set_2d_coordinates_with_params_"
                | "Molecule.with_3d_coordinates"
                | "Molecule.with_3d_coordinates_with_params"
                | "Molecule.set_3d_coordinates_"
                | "Molecule.set_3d_coordinates_with_params_"
                | "Molecule.with_added_3d_conformer"
                | "Molecule.with_added_3d_conformer_with_params"
                | "Molecule.add_3d_conformer_"
                | "Molecule.add_3d_conformer_with_params_"
                | "Molecule.with_only_3d_conformer"
                | "Molecule.with_only_3d_conformer_with_params"
                | "Molecule.set_only_3d_conformer_"
                | "Molecule.set_only_3d_conformer_with_params_"
                | "Molecule.with_cleared_3d_conformers"
                | "Molecule.clear_3d_conformers_"
        );
        // These exact project-native BATCH statuses remain declared in the source registry.
        let native_batch = matches!(
            contract.semantic_id,
            "types.MoleculeBatch"
                | "types.BatchRecord"
                | "types.BatchError"
                | "types.BatchErrorMode"
                | "types.BatchValidationError"
                | "types.BatchParams"
                | "types.BatchExportReport"
                | "MoleculeBatch.len"
                | "MoleculeBatch.is_empty"
                | "MoleculeBatch.valid_mask"
                | "MoleculeBatch.invalid_mask"
                | "MoleculeBatch.valid_count"
                | "MoleculeBatch.invalid_count"
                | "MoleculeBatch.errors"
                | "MoleculeBatch.parallel_jobs"
                | "MoleculeBatch.progress_bar"
                | "MoleculeBatch.to_list"
                | "MoleculeBatch.with_valid_records"
                | "MoleculeBatch.with_parallel_jobs"
                | "MoleculeBatch.with_progress_bar"
                | "MoleculeBatch.from_records"
                | "MoleculeBatch.from_smiles_list"
                | "MoleculeBatch.from_smiles_list_with_params"
                | "MoleculeBatch.sanitize"
                | "MoleculeBatch.sanitize_with_params"
                | "MoleculeBatch.with_hydrogens"
                | "MoleculeBatch.with_hydrogens_with_params"
                | "MoleculeBatch.without_hydrogens"
                | "MoleculeBatch.without_hydrogens_with_params"
                | "MoleculeBatch.with_kekulized_bonds"
                | "MoleculeBatch.with_kekulized_bonds_with_params"
                | "MoleculeBatch.with_2d_coordinates"
                | "MoleculeBatch.with_2d_coordinates_with_params"
                | "types.BatchQueryParams"
                | "MoleculeBatch.to_smiles_list"
                | "MoleculeBatch.to_smiles_list_with_params"
                | "MoleculeBatch.dg_bounds_matrix_list"
                | "MoleculeBatch.dg_bounds_matrix_list_with_params"
                | "MoleculeBatch.to_svg_list"
                | "MoleculeBatch.to_svg_list_with_params"
                | "types.BatchImageParams"
                | "types.BatchImageError"
                | "MoleculeBatch.to_images"
                | "MoleculeBatch.to_images_with_params"
                | "MoleculeBatch.fingerprint_atom_pair_list"
                | "MoleculeBatch.fingerprint_atom_pair_list_with_params"
                | "MoleculeBatch.fingerprint_atom_pair_sparse_count_list"
                | "MoleculeBatch.fingerprint_atom_pair_sparse_count_list_with_params"
                | "MoleculeBatch.fingerprint_atom_pair_count_list"
                | "MoleculeBatch.fingerprint_atom_pair_count_list_with_params"
                | "MoleculeBatch.fingerprint_atom_pair_sparse_bits_list"
                | "MoleculeBatch.fingerprint_atom_pair_sparse_bits_list_with_params"
                | "MoleculeBatch.fingerprint_layered_list"
                | "MoleculeBatch.fingerprint_layered_list_with_params"
                | "MoleculeBatch.fingerprint_layered_with_output_list"
                | "MoleculeBatch.fingerprint_layered_with_output_list_with_params"
                | "MoleculeBatch.pattern_fingerprint_list"
                | "MoleculeBatch.pattern_fingerprint_list_with_params"
                | "MoleculeBatch.fingerprint_morgan_list"
                | "MoleculeBatch.fingerprint_morgan_list_with_params"
                | "types.BatchFingerprintOutput"
                | "MoleculeBatch.fingerprint_atom_pair_with_output_list"
                | "MoleculeBatch.fingerprint_atom_pair_with_output_list_with_params"
                | "MoleculeBatch.fingerprint_morgan_with_output_list"
                | "MoleculeBatch.fingerprint_morgan_with_output_list_with_params"
                | "MoleculeBatch.fingerprint_morgan_list_with_generator_params"
                | "MoleculeBatch.fingerprint_morgan_with_output_list_with_generator_params"
                | "types.BatchFingerprintAdditionalOutput"
                | "types.BatchFingerprintOutputError"
                | "BatchFingerprintAdditionalOutput.atom_counts"
                | "BatchFingerprintAdditionalOutput.atom_to_bits"
                | "BatchFingerprintAdditionalOutput.bit_info_map"
                | "BatchFingerprintAdditionalOutput.bit_paths"
                | "BatchFingerprintAdditionalOutput.atoms_per_bit"
                | "BatchFingerprintOutput.fingerprint"
                | "BatchFingerprintOutput.additional_output"
        );
        let expected = if native_coordinates || native_batch {
            FunctionStatus::Native
        } else if fuzzy {
            FunctionStatus::Parity { reference: "RDKit" }
        } else if matches!(
            contract.semantic_id,
            "Molecule.to_svg" | "Molecule.write_svg"
        ) {
            canonical_svg_identity_status()
        } else if matches!(
            contract.semantic_id,
            "types.LigandRef"
                | "types.TetrahedralStereo"
                | "types.StereoReadError"
                | "Molecule.tetrahedral_stereo"
                | "Molecule.perceive_stereochemistry"
                | "Molecule.find_chiral_centers"
        ) {
            FunctionStatus::Native
        } else {
            FunctionStatus::Experimental
        };
        assert_eq!(contract.status, expected, "{}", contract.semantic_id);
        if fuzzy {
            parity.push(contract.semantic_id);
        }
    }
    assert_eq!(
        parity.len(),
        if cfg!(feature = "cap-fingerprints") {
            4
        } else {
            0
        }
    );
    for operation in cosmolkit::operation_specs() {
        let binding = entry(&format!("Molecule.{}", operation.method));
        assert_eq!(operation.status, binding.status);
    }
}
