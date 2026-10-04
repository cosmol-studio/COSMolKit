use std::collections::HashSet;

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
            "Molecule.with_2d_coordinates",
            "Molecule.with_2d_coordinates_with_params",
            "types.DrawingError",
            "Molecule.to_svg",
            "Molecule.to_png",
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
            "Molecule.molecular_weight",
            "Molecule.exact_molecular_weight",
            "Molecule.molecular_formula",
            "types.DescriptorError",
            "types.DescriptorReadError",
            "Molecule.num_heavy_atoms",
            "Molecule.total_atom_count",
            "Molecule.num_rings",
            "Molecule.num_heterocycles",
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
        assert_eq!(row.status, FunctionStatus::Experimental);
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
            "crate::OperationError",
        ),
        (
            "Molecule.exact_molecular_weight",
            "exact_molecular_weight",
            "exactMolecularWeight",
            "f64",
            "crate::OperationError",
        ),
        (
            "Molecule.molecular_formula",
            "molecular_formula",
            "molecularFormula",
            "String",
            "crate::OperationError",
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
        let expected = if fuzzy {
            FunctionStatus::Parity { reference: "RDKit" }
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
