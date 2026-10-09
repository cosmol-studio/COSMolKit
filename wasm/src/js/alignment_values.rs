//! Alignment result and error transport, without algorithm implementations.
#[cfg(feature = "cap-alignment")]
use crate::alignment_parameters::AlignmentAtomMap;
use cosmolkit_wasm::rust as ck;
use js_sys::{Array, Error, Reflect};
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;

#[cfg(feature = "cap-alignment")]
#[wasm_bindgen]
pub struct AlignmentTransform {
    pub(crate) inner: ck::AlignmentTransform,
}

#[cfg(feature = "cap-alignment")]
#[wasm_bindgen]
impl AlignmentTransform {
    #[wasm_bindgen(unchecked_return_type = "number[][]")]
    pub fn matrix(&self) -> JsValue {
        self.inner
            .matrix()
            .iter()
            .map(|row| row.iter().map(|n| JsValue::from_f64(*n)).collect::<Array>())
            .collect::<Array>()
            .into()
    }
}

#[cfg(feature = "cap-alignment")]
#[wasm_bindgen]
pub struct AlignmentResult {
    pub(crate) inner: ck::AlignmentResult,
}

#[cfg(feature = "cap-alignment")]
#[wasm_bindgen]
impl AlignmentResult {
    pub fn rmsd(&self) -> f64 {
        self.inner.rmsd()
    }
    pub fn transform(&self) -> AlignmentTransform {
        AlignmentTransform {
            inner: *self.inner.transform(),
        }
    }
    #[wasm_bindgen(js_name = atomMap)]
    pub fn atom_map(&self) -> Vec<AlignmentAtomMap> {
        self.inner
            .atom_map()
            .iter()
            .map(|inner| AlignmentAtomMap { inner: *inner })
            .collect()
    }
}

#[cfg(feature = "cap-alignment")]
#[wasm_bindgen]
pub struct ConformerRmsd {
    pub(crate) inner: ck::ConformerRmsd,
}

#[cfg(feature = "cap-alignment")]
#[wasm_bindgen]
impl ConformerRmsd {
    pub fn rmsd(&self) -> f64 {
        self.inner.rmsd()
    }
    #[wasm_bindgen(js_name = probeConformerId)]
    pub fn probe_conformer_id(&self) -> usize {
        self.inner.probe_conformer_id()
    }
    #[wasm_bindgen(js_name = referenceConformerId)]
    pub fn reference_conformer_id(&self) -> usize {
        self.inner.reference_conformer_id()
    }
}

#[cfg(feature = "cap-alignment")]
#[wasm_bindgen]
pub struct ConformerAlignmentReport {
    pub(crate) inner: ck::ConformerAlignmentReport,
}

#[cfg(feature = "cap-alignment")]
#[wasm_bindgen]
impl ConformerAlignmentReport {
    pub fn rmsds(&self) -> Vec<f64> {
        self.inner.rmsds().to_vec()
    }
}

#[cfg(feature = "cap-alignment")]
#[wasm_bindgen]
pub struct AlignmentError {
    pub(crate) inner: ck::AlignmentError,
}

#[cfg(feature = "cap-alignment")]
#[wasm_bindgen]
impl AlignmentError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "alignment".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        alignment_kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter)]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
}

#[cfg(feature = "cap-alignment")]
fn alignment_kind(source: &ck::AlignmentError) -> &'static str {
    use ck::AlignmentError as E;
    match source {
        E::InvalidCoordinates(..) => "InvalidCoordinates",
        E::Matching(..) => "Matching",
        E::QueryGraph(..) => "QueryGraph",
        E::QueryParse(..) => "QueryParse",
        E::ThreadSelection(..) => "ThreadSelection",
        E::NoConformers => "NoConformers",
        E::ConformerNotFound { .. } => "ConformerNotFound",
        E::EmptyAtomMap => "EmptyAtomMap",
        E::ProbeAtomOutOfRange { .. } => "ProbeAtomOutOfRange",
        E::ReferenceAtomOutOfRange { .. } => "ReferenceAtomOutOfRange",
        E::WeightCountMismatch { .. } => "WeightCountMismatch",
        E::NonPositiveWeight { .. } => "NonPositiveWeight",
        E::NoSubstructureMatch => "NoSubstructureMatch",
        E::TerminalGroupSymmetrization { .. } => "TerminalGroupSymmetrization",
        E::NumericalPrecondition { .. } => "NumericalPrecondition",
        E::WorkerTerminated => "WorkerTerminated",
    }
}

pub(crate) fn set(error: &JsValue, name: &str, value: JsValue) -> Result<(), JsValue> {
    if !Reflect::set(error, &name.into(), &value)? {
        return Err(js_sys::TypeError::new("cannot set binding error detail").into());
    }
    Ok(())
}

#[cfg(feature = "cap-alignment")]
pub(crate) fn alignment_error(source: &ck::AlignmentError) -> Result<JsValue, JsValue> {
    use ck::AlignmentError as E;
    let error = Error::new(&source.to_string());
    error.set_name("AlignmentError");
    let error: JsValue = error.into();
    set(&error, "domain", "alignment".into())?;
    set(&error, "kind", alignment_kind(source).into())?;
    set(
        &error,
        "detail",
        AlignmentError {
            inner: source.clone(),
        }
        .into(),
    )?;
    match source {
        E::ConformerNotFound { id } => set(&error, "id", (*id).into())?,
        E::ProbeAtomOutOfRange { index, atom_count }
        | E::ReferenceAtomOutOfRange { index, atom_count } => {
            set(&error, "index", JsValue::from_f64(*index as f64))?;
            set(&error, "atomCount", JsValue::from_f64(*atom_count as f64))?;
        }
        E::WeightCountMismatch {
            map_len,
            weight_len,
        } => {
            set(&error, "mapLen", JsValue::from_f64(*map_len as f64))?;
            set(&error, "weightLen", JsValue::from_f64(*weight_len as f64))?;
        }
        E::NonPositiveWeight { index } => set(&error, "index", JsValue::from_f64(*index as f64))?,
        E::TerminalGroupSymmetrization { message } | E::NumericalPrecondition { message } => {
            set(&error, "message", (*message).into())?
        }
        _ => {}
    }
    if let Some(cause) = source.source() {
        set(&error, "cause", source_error(cause)?)?;
    }
    Ok(error)
}

pub(crate) fn source_error(source: &(dyn RustError + 'static)) -> Result<JsValue, JsValue> {
    #[cfg(feature = "cap-search")]
    if let Some(e) = source.downcast_ref::<ck::SmartsParseError>() {
        return crate::query_construction::parse_error(e);
    }
    #[cfg(feature = "cap-reaction")]
    if let Some(e) = source.downcast_ref::<ck::ReactionModelError>() {
        return crate::reaction_errors::model_error(e);
    }
    #[cfg(feature = "cap-reaction")]
    if let Some(e) = source.downcast_ref::<ck::ReactionParseError>() {
        return crate::reaction_errors::parse_error(e);
    }
    #[cfg(feature = "cap-reaction")]
    if let Some(e) = source.downcast_ref::<ck::ReactionRunError>() {
        return crate::reaction_errors::run_error(e);
    }
    #[cfg(feature = "cap-reaction")]
    if let Some(e) = source.downcast_ref::<ck::ReactionApplyError>() {
        return crate::reaction_errors::apply_error(e);
    }
    #[cfg(feature = "cap-reaction")]
    if let Some(e) = source.downcast_ref::<ck::ReactionProductError>() {
        return crate::reaction_errors::product_error(e);
    }
    #[cfg(feature = "cap-reaction")]
    if let Some(e) = source.downcast_ref::<ck::ReactionWriteError>() {
        return crate::reaction_errors::write_error(e);
    }
    #[cfg(feature = "cap-reaction")]
    if let Some(e) = source.downcast_ref::<ck::ReactionValidationError>() {
        return crate::reaction_errors::validation_error(e);
    }
    #[cfg(feature = "cap-reaction")]
    if let Some(e) = source.downcast_ref::<ck::ReactionInitializationError>() {
        return crate::reaction_errors::initialization_error(e);
    }
    #[cfg(feature = "cap-search")]
    if let Some(e) = source.downcast_ref::<ck::SubstructMatchError>() {
        return crate::search_errors::substruct_error(e);
    }
    #[cfg(feature = "cap-search")]
    if let Some(e) = source.downcast_ref::<ck::MatchError>() {
        return crate::search_errors::match_error(e);
    }
    #[cfg(feature = "cap-search")]
    if let Some(e) = source.downcast_ref::<ck::QueryCompileError>() {
        return crate::search_errors::compile_error(e);
    }
    #[cfg(feature = "cap-search")]
    if let Some(e) = source.downcast_ref::<ck::SmartsWriteError>() {
        return crate::search_errors::write_error(e);
    }

    if let Some(e) = source.downcast_ref::<ck::ValenceError>() {
        return crate::valence_errors::valence_error(e);
    }
    if let Some(e) = source.downcast_ref::<ck::ChemistryProblemError>() {
        return crate::sanitize::problem_error(e);
    }

    if let Some(e) = source.downcast_ref::<ck::MatrixError>() {
        return crate::matrix_errors::matrix_error(e);
    }

    if let Some(e) = source.downcast_ref::<ck::MolWriteError>() {
        return crate::mol_write_errors::mol_write_error(e);
    }

    if let Some(source) = source.downcast_ref::<ck::XyzReadError>() {
        return crate::xyz_mol2_errors::xyz_read_error(source);
    }
    if let Some(source) = source.downcast_ref::<ck::XyzWriteError>() {
        return crate::xyz_mol2_errors::xyz_write_error(source);
    }
    if let Some(source) = source.downcast_ref::<ck::Mol2ReadError>() {
        return crate::xyz_mol2_errors::mol2_read_error(source);
    }
    if let Some(source) = source.downcast_ref::<ck::Mol2PostError>() {
        return crate::xyz_mol2_errors::mol2_post_error(source);
    }
    if let Some(source) = source.downcast_ref::<ck::PropertyStringError>() {
        return crate::property_strings::string_error(source);
    }

    if let Some(source) = source.downcast_ref::<ck::SdfError>() {
        return crate::io_errors::sdf_error(source);
    }
    if let Some(source) = source.downcast_ref::<ck::MolecularIoError>() {
        return crate::io_errors::io_error(source);
    }
    if let Some(source) = source.downcast_ref::<ck::PropertyValueError>() {
        return crate::property_values::property_error(source);
    }

    #[cfg(feature = "cap-bio")]
    if let Some(source) = source.downcast_ref::<ck::BioMoleculeError>() {
        return crate::bio_conversion::molecule_error(source);
    }
    #[cfg(feature = "cap-bio")]
    if let Some(source) = source.downcast_ref::<ck::BioMoleculeConversionError>() {
        return crate::bio_conversion::conversion_error(source);
    }

    #[cfg(feature = "cap-fingerprints")]
    if let Some(source) = source.downcast_ref::<ck::TopologicalTorsionPathScoreError>() {
        return crate::path_codes::score_error(source);
    }
    #[cfg(feature = "cap-fingerprints")]
    if let Some(source) = source.downcast_ref::<ck::AtomCodeExplanationError>() {
        return crate::path_codes::explanation_error(source);
    }

    #[cfg(feature = "cap-conformer")]
    if let Some(source) = source.downcast_ref::<ck::ConformerError>() {
        return crate::conformer_errors::params_error(source);
    }
    #[cfg(feature = "cap-conformer")]
    if let Some(source) = source.downcast_ref::<ck::ConformerRunError>() {
        return crate::conformer_errors::run_error(source);
    }
    #[cfg(feature = "cap-bio")]
    if let Some(source) = source.downcast_ref::<ck::ProteinReadError>() {
        return crate::bio_protein::read_error(source);
    }
    #[cfg(feature = "cap-bio")]
    if let Some(source) = source.downcast_ref::<ck::BioPdbWriteError>() {
        return crate::bio_write_errors::pdb_error(source);
    }
    #[cfg(feature = "cap-bio")]
    if let Some(source) = source.downcast_ref::<ck::BioMmcifWriteError>() {
        return crate::bio_write_errors::mmcif_error(source);
    }
    #[cfg(feature = "cap-bio")]
    if let Some(source) = source.downcast_ref::<ck::BioSelectionParseError>() {
        return crate::bio_selection::parse_error(source);
    }
    #[cfg(feature = "cap-bio")]
    if let Some(source) = source.downcast_ref::<ck::BioSelectionMatchError>() {
        return crate::bio_selection::match_error(source);
    }
    #[cfg(feature = "cap-bio")]
    if let Some(source) = source.downcast_ref::<ck::BioSelectionCopyError>() {
        return crate::bio_selection::copy_error(source);
    }
    #[cfg(feature = "cap-bio")]
    if let Some(source) = source.downcast_ref::<ck::BioSelectionCopyCause>() {
        return crate::bio_selection::copy_cause_error(source);
    }
    #[cfg(feature = "cap-bio")]
    if let Some(source) = source.downcast_ref::<ck::BioOperationError>() {
        return crate::bio_selection::operation_error(source);
    }
    #[cfg(feature = "cap-bio")]
    if let Some(source) = source.downcast_ref::<ck::ProteinProjectionError>() {
        return crate::bio_selection::protein_projection_error(source);
    }
    #[cfg(feature = "cap-bio")]
    if let Some(source) = source.downcast_ref::<ck::BioStructureError>() {
        return crate::bio_hierarchy_errors::structure_error(source);
    }
    #[cfg(feature = "cap-bio")]
    if let Some(source) = source.downcast_ref::<ck::BioRowModelError>() {
        return crate::bio_hierarchy_errors::row_model_error(source);
    }
    #[cfg(feature = "cap-bio")]
    if let Some(source) = source.downcast_ref::<ck::BioRowChainError>() {
        return crate::bio_hierarchy_errors::row_chain_error(source);
    }
    #[cfg(feature = "cap-bio")]
    if let Some(source) = source.downcast_ref::<ck::BioRowTraverseError>() {
        return crate::bio_hierarchy_errors::row_traverse_error(source);
    }
    #[cfg(feature = "cap-bio")]
    if let Some(source) = source.downcast_ref::<ck::BioReadError>() {
        return crate::bio_read_errors::read_error(source);
    }
    #[cfg(feature = "cap-bio")]
    if let Some(source) = source.downcast_ref::<ck::BioPdbReadError>() {
        return crate::bio_read_errors::pdb_error(source);
    }
    #[cfg(feature = "cap-bio")]
    if let Some(source) = source.downcast_ref::<ck::BioMmcifReadError>() {
        return crate::bio_read_errors::mmcif_error(source);
    }
    #[cfg(feature = "cap-bio")]
    if let Some(source) = source.downcast_ref::<ck::ResidueSequenceError>() {
        return crate::bio_residue::sequence_error(source);
    }
    #[cfg(feature = "cap-fingerprints")]
    if let Some(source) = source.downcast_ref::<ck::MorganReadError>() {
        return crate::fingerprint_source_errors::morgan_error(source);
    }
    #[cfg(feature = "cap-forcefields")]
    if let Some(source) = source.downcast_ref::<ck::UffOptimizationError>() {
        return crate::uff::uff_optimization_error(source);
    }
    #[cfg(feature = "cap-forcefields")]
    if let Some(source) = source.downcast_ref::<ck::MmffOptimizationError>() {
        return crate::mmff::mmff_optimization_error(source);
    }
    #[cfg(feature = "cap-fingerprints")]
    if let Some(source) = source.downcast_ref::<ck::MoleculeHashError>() {
        return crate::hashing::hash_error(source);
    }
    #[cfg(feature = "cap-fingerprints")]
    if let Some(source) = source.downcast_ref::<ck::CipRankError>() {
        return crate::hashing::cip_error(source);
    }
    #[cfg(feature = "cap-forcefields")]
    if let Some(source) = source.downcast_ref::<ck::UffParameterQueryError>() {
        return crate::forcefield_properties::uff_query_error(source);
    }
    #[cfg(feature = "cap-forcefields")]
    if let Some(source) = source.downcast_ref::<ck::UffParameterError>() {
        return crate::forcefield_properties::uff_parameter_error(source);
    }
    #[cfg(feature = "cap-forcefields")]
    if let Some(source) = source.downcast_ref::<ck::MmffMolPropertiesError>() {
        return crate::forcefield_properties::mmff_properties_error(source);
    }
    #[cfg(feature = "cap-fingerprints")]
    if let Some(source) = source.downcast_ref::<ck::TopologicalFingerprintError>() {
        return crate::topological::topological_error(source);
    }
    #[cfg(feature = "cap-fingerprints")]
    if let Some(source) = source.downcast_ref::<ck::MaccsFingerprintError>() {
        return crate::maccs::maccs_error(source);
    }
    #[cfg(feature = "cap-fingerprints")]
    if let Some(source) = source.downcast_ref::<ck::TopologicalTorsionReadError>() {
        return crate::fingerprint_source_errors::topological_torsion_error(source);
    }
    #[cfg(feature = "cap-fingerprints")]
    if let Some(source) = source.downcast_ref::<ck::PatternFingerprintError>() {
        return crate::layered_pattern_errors::pattern_error(source);
    }
    #[cfg(feature = "cap-fingerprints")]
    if let Some(source) = source.downcast_ref::<ck::LayeredFingerprintError>() {
        return crate::layered_pattern_errors::layered_error(source);
    }
    #[cfg(all(feature = "cap-fingerprints", feature = "cap-batch"))]
    if let Some(source) = source.downcast_ref::<ck::BatchFingerprintOutputError>() {
        return crate::fingerprint_source_errors::output_error(source);
    }
    #[cfg(feature = "cap-fingerprints")]
    if let Some(source) = source.downcast_ref::<ck::FingerprintJsonError>() {
        return crate::fingerprint_source_errors::json_error(source);
    }
    #[cfg(feature = "cap-fingerprints")]
    if let Some(source) = source.downcast_ref::<ck::FingerprintPreparationError>() {
        return crate::fingerprint_source_errors::preparation_error(source);
    }
    #[cfg(feature = "cap-fingerprints")]
    if let Some(source) = source.downcast_ref::<ck::AtomPairReadError>() {
        return crate::fingerprint_source_errors::atom_pair_error(source);
    }
    #[cfg(feature = "cap-fingerprints")]
    if let Some(source) = source.downcast_ref::<ck::FingerprintError>() {
        return crate::fingerprint_errors::fingerprint_error(source);
    }
    #[cfg(all(feature = "cap-batch", feature = "cap-depict"))]
    if let Some(source) = source.downcast_ref::<ck::BatchImageError>() {
        return crate::image_errors::batch_image_error(source);
    }
    #[cfg(feature = "cap-descriptors")]
    if let Some(source) = source.downcast_ref::<ck::DescriptorReadError>() {
        return crate::descriptor_errors::descriptor_read_error(source);
    }
    #[cfg(feature = "cap-descriptors")]
    if let Some(source) = source.downcast_ref::<ck::DescriptorError>() {
        return crate::descriptor_errors::descriptor_error(source);
    }
    #[cfg(feature = "cap-depict")]
    if let Some(source) = source.downcast_ref::<ck::DrawingWriteError>() {
        return crate::image_errors::drawing_write_error(source);
    }
    if let Some(source) = source.downcast_ref::<std::io::Error>() {
        return crate::io_errors::filesystem_error(source);
    }

    #[cfg(feature = "cap-depict")]
    if let Some(source) = source.downcast_ref::<ck::DrawingError>() {
        return crate::drawing_errors::drawing_error(source);
    }

    if let Some(source) = source.downcast_ref::<ck::SanitizeError>() {
        return crate::transform_errors::sanitize_error(source);
    }
    if let Some(source) = source.downcast_ref::<ck::HydrogenError>() {
        return crate::transform_errors::hydrogen_error(source);
    }
    if let Some(source) = source.downcast_ref::<ck::KekulizeError>() {
        return crate::transform_errors::kekulize_error(source);
    }
    #[cfg(feature = "cap-depict")]
    if let Some(source) = source.downcast_ref::<ck::Coordinate2DTemplateError>() {
        return crate::depict_errors::template_error(source);
    }
    #[cfg(feature = "cap-depict")]
    if let Some(source) = source.downcast_ref::<ck::Coordinate2DLayoutError>() {
        return crate::depict_errors::layout_error(source);
    }
    #[cfg(feature = "cap-depict")]
    if let Some(source) = source.downcast_ref::<ck::Coordinate2DError>() {
        return crate::transform_errors::coordinate_2d_error(source);
    }
    if let Some(source) = source.downcast_ref::<ck::SmilesError>() {
        return crate::smiles_boundary::smiles_error(source);
    }
    #[cfg(feature = "cap-batch")]
    if let Some(source) = source.downcast_ref::<ck::BatchValidationError>() {
        return crate::batch_boundary::batch_validation_error(source);
    }
    #[cfg(feature = "cap-batch")]
    if let Some(source) = source.downcast_ref::<ck::BatchError>() {
        return crate::batch_boundary::batch_record_error(source);
    }
    if let Some(source) = source.downcast_ref::<ck::SmilesWriteError>() {
        return crate::smiles_boundary::smiles_write_error(source);
    }
    if let Some(source) = source.downcast_ref::<ck::AromaticityError>() {
        return crate::aromaticity_boundary::aromaticity_error(source);
    }
    #[cfg(feature = "cap-alignment")]
    if let Some(source) = source.downcast_ref::<ck::AlignmentError>() {
        return alignment_error(source);
    }
    if let Some(source) = source.downcast_ref::<ck::OperationError>() {
        return operation_error(source);
    }
    // Python's canonical mapper also keeps unknown source types as native
    // errors with their original text and recursive source chain.
    let error: JsValue = Error::new(&source.to_string()).into();
    if let Some(cause) = source.source() {
        set(&error, "cause", source_error(cause)?)?;
    }
    Ok(error)
}

pub(crate) fn operation_error(source: &ck::OperationError) -> Result<JsValue, JsValue> {
    use ck::OperationError as E;
    let kind = match source {
        E::AtomPropertyIndex { .. } => "AtomPropertyIndex",
        E::ReservedAtomPropertyKey { .. } => "ReservedAtomPropertyKey",
        #[cfg(feature = "cap-reaction")]
        E::ReactionRun(..) => "ReactionRun",
        #[cfg(feature = "cap-reaction")]
        E::ReactionApply(..) => "ReactionApply",
        #[cfg(feature = "full")]
        E::Enumeration(..) => "Enumeration",
        E::AtomProperty(..) => "AtomProperty",
        E::BondProperty(..) => "BondProperty",
        E::InvalidReconstructionOrigin { .. } => "InvalidReconstructionOrigin",
        #[cfg(feature = "cap-alignment")]
        E::Alignment(..) => "Alignment",
        #[cfg(feature = "cap-conformer")]
        E::Conformer(..) => "Conformer",
        #[cfg(feature = "full")]
        E::Tautomer(..) => "Tautomer",
        E::UnsupportedFeature { .. } => "UnsupportedFeature",
        E::Unsupported { .. } => "Unsupported",
        #[cfg(feature = "cap-transforms")]
        E::Fragments(..) => "Fragments",
        E::EmptyFragments => "EmptyFragments",
        E::OutputMismatch { .. } => "OutputMismatch",
        E::AccessDenied { .. } => "AccessDenied",
        E::BlockCheckedOut { .. } => "BlockCheckedOut",
        E::BlockNotCheckedOut { .. } => "BlockNotCheckedOut",
        E::IncompleteCommit { .. } => "IncompleteCommit",
        E::TopologyEditContract { .. } => "TopologyEditContract",
        E::MappingContract { .. } => "MappingContract",
        E::InvalidTopologyMapping { .. } => "InvalidTopologyMapping",
        E::AutoRemapContract { .. } => "AutoRemapContract",
        E::OperationContract { .. } => "OperationContract",
        E::SemanticPreconditionContract { .. } => "SemanticPreconditionContract",
        E::CoordinateAppendRequiresValues { .. } => "CoordinateAppendRequiresValues",
        E::DerivedEffectContract { .. } => "DerivedEffectContract",
        E::CipStateContract { .. } => "CipStateContract",
        E::InvalidTopology(..) => "InvalidTopology",
        #[cfg(feature = "cap-fingerprints")]
        E::AtomCode(..) => "AtomCode",
        E::InvalidTopologyEdit(..) => "InvalidTopologyEdit",
        E::InvalidCoordinates(..) => "InvalidCoordinates",
        E::InvalidProperty(..) => "InvalidProperty",
        E::InvalidPropertyList { .. } => "InvalidPropertyList",
        E::InvalidDerivedCache { .. } => "InvalidDerivedCache",
        E::Valence(..) => "Valence",
        E::Radical(..) => "Radical",
        E::Rings(..) => "Rings",
        E::PotentialStereo(..) => "PotentialStereo",
        E::Stereo(..) => "Stereo",
        E::CipLabeler(..) => "CipLabeler",
        E::Transform(..) => "Transform",
        E::CoordinateInput(..) => "CoordinateInput",
        #[cfg(feature = "cap-depict")]
        E::Coordinate2D(..) => "Coordinate2D",
        #[cfg(feature = "cap-forcefields")]
        E::UffOptimization(..) => "UffOptimization",
        #[cfg(feature = "cap-forcefields")]
        E::MmffOptimization(..) => "MmffOptimization",
        E::Kekulize(..) => "Kekulize",
        E::Aromaticity(..) => "Aromaticity",
        E::Sanitize(..) => "Sanitize",
        E::Hydrogen(..) => "Hydrogen",
        E::InvalidAlgorithmResult { .. } => "InvalidAlgorithmResult",
        E::Algorithm { .. } => "Algorithm",
        #[cfg(feature = "cap-hashing")]
        E::Scaffold(..) => "Scaffold",
    };
    let error = Error::new(&source.to_string());
    error.set_name("OperationError");
    let error: JsValue = error.into();
    set(&error, "domain", "operation".into())?;
    set(&error, "kind", kind.into())?;
    match source {
        E::AtomPropertyIndex { atom, atom_count } => {
            set(&error,"atomIndex",(atom.index() as f64).into())?;
            set(&error,"atomCount",(*atom_count as f64).into())?;
        }
        E::ReservedAtomPropertyKey { key } => set(&error,"key",key.as_str().into())?,
        _ => {},
    }
    match source {
        #[cfg(feature = "cap-reaction")]
        E::ReactionRun(cause) => set(&error, "cause", crate::reaction_errors::run_error(cause)?)?,
        #[cfg(feature = "cap-reaction")]
        E::ReactionApply(cause) => {
            set(&error, "cause", crate::reaction_errors::apply_error(cause)?)?
        }
        _ => {
            if let Some(cause) = source.source() {
                set(&error, "cause", source_error(cause)?)?;
            }
        }
    }
    Ok(error)
}
