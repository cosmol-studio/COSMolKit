//! Structured errors copied from public facade values, with real causes.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use js_sys::Array;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct SanitizeError {
    inner: ck::SanitizeError,
}
#[wasm_bindgen]
impl SanitizeError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "sanitize".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        sanitize_error_kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number | null")]
    pub fn bits(&self) -> JsValue {
        match self.inner {
            ck::SanitizeError::InvalidOperations { bits, .. } => bits.into(),
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=unknownBits,unchecked_return_type="number | null")]
    pub fn unknown_bits(&self) -> JsValue {
        match self.inner {
            ck::SanitizeError::InvalidOperations { unknown_bits, .. } => unknown_bits.into(),
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter, unchecked_return_type = "SanitizeStage | null")]
    pub fn stage(&self) -> JsValue {
        use ck::SanitizeError as E;
        match &self.inner {
            E::InvalidTopology { stage, .. }
            | E::InvalidQueryState { stage, .. }
            | E::Cleanup { stage, .. }
            | E::Properties { stage, .. }
            | E::Rings { stage, .. }
            | E::Kekulize { stage, .. }
            | E::Radicals { stage, .. }
            | E::Aromaticity { stage, .. }
            | E::Conjugation { stage, .. }
            | E::Hybridization { stage, .. }
            | E::Atropisomers { stage, .. }
            | E::Chirality { stage, .. }
            | E::AdjustHs { stage, .. } => {
                JsValue::from(crate::transform_parameters::stage(*stage) as u32)
            }
            E::InvalidOperations { .. }
            | E::MoleculeProperty(..)
            | E::BondProperty(..)
            | E::AtomProperty(..) => JsValue::NULL,
        }
    }
}
fn sanitize_error_kind(source: &ck::SanitizeError) -> &'static str {
    use ck::SanitizeError as E;
    match source {
        E::MoleculeProperty(..) => "MoleculeProperty",
        E::BondProperty(..) => "BondProperty",
        E::AtomProperty(..) => "AtomProperty",
        E::InvalidOperations { .. } => "InvalidOperations",
        E::InvalidTopology { .. } => "InvalidTopology",
        E::InvalidQueryState { .. } => "InvalidQueryState",
        E::Cleanup { .. } => "Cleanup",
        E::Properties { .. } => "Properties",
        E::Rings { .. } => "Rings",
        E::Kekulize { .. } => "Kekulize",
        E::Radicals { .. } => "Radicals",
        E::Aromaticity { .. } => "Aromaticity",
        E::Conjugation { .. } => "Conjugation",
        E::Hybridization { .. } => "Hybridization",
        E::Atropisomers { .. } => "Atropisomers",
        E::Chirality { .. } => "Chirality",
        E::AdjustHs { .. } => "AdjustHs",
    }
}
pub(crate) fn sanitize_error(source: &ck::SanitizeError) -> Result<JsValue, JsValue> {
    use ck::SanitizeError as E;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("SanitizeError");
    let error: JsValue = error.into();
    set(&error, "domain", "sanitize".into())?;
    set(&error, "kind", sanitize_error_kind(source).into())?;
    set(
        &error,
        "detail",
        SanitizeError {
            inner: source.clone(),
        }
        .into(),
    )?;
    match source {
        E::MoleculeProperty(..) | E::BondProperty(..) | E::AtomProperty(..) => {}
        E::InvalidOperations { bits, unknown_bits } => {
            set(&error, "bits", JsValue::from_f64(*bits as f64))?;
            set(
                &error,
                "unknownBits",
                JsValue::from_f64(*unknown_bits as f64),
            )?;
        }
        E::InvalidTopology { stage, source } => {
            set(
                &error,
                "stage",
                JsValue::from(crate::transform_parameters::stage(*stage) as u32),
            )?;
            let _ = source; // Retained through Error::source below.
        }
        E::InvalidQueryState { stage, source } => {
            set(
                &error,
                "stage",
                JsValue::from(crate::transform_parameters::stage(*stage) as u32),
            )?;
            let _ = source; // Retained through Error::source below.
        }
        E::Cleanup { stage, source } => {
            set(
                &error,
                "stage",
                JsValue::from(crate::transform_parameters::stage(*stage) as u32),
            )?;
            let _ = source; // Retained through Error::source below.
        }
        E::Properties { stage, source } => {
            set(
                &error,
                "stage",
                JsValue::from(crate::transform_parameters::stage(*stage) as u32),
            )?;
            let _ = source; // Retained through Error::source below.
        }
        E::Rings { stage, source } => {
            set(
                &error,
                "stage",
                JsValue::from(crate::transform_parameters::stage(*stage) as u32),
            )?;
            let _ = source; // Retained through Error::source below.
        }
        E::Kekulize { stage, source } => {
            set(
                &error,
                "stage",
                JsValue::from(crate::transform_parameters::stage(*stage) as u32),
            )?;
            let _ = source; // Retained through Error::source below.
        }
        E::Radicals { stage, source } => {
            set(
                &error,
                "stage",
                JsValue::from(crate::transform_parameters::stage(*stage) as u32),
            )?;
            let _ = source; // Retained through Error::source below.
        }
        E::Aromaticity { stage, source } => {
            set(
                &error,
                "stage",
                JsValue::from(crate::transform_parameters::stage(*stage) as u32),
            )?;
            let _ = source; // Retained through Error::source below.
        }
        E::Conjugation { stage, source } => {
            set(
                &error,
                "stage",
                JsValue::from(crate::transform_parameters::stage(*stage) as u32),
            )?;
            let _ = source; // Retained through Error::source below.
        }
        E::Hybridization { stage, source } => {
            set(
                &error,
                "stage",
                JsValue::from(crate::transform_parameters::stage(*stage) as u32),
            )?;
            let _ = source; // Retained through Error::source below.
        }
        E::Atropisomers { stage, source } => {
            set(
                &error,
                "stage",
                JsValue::from(crate::transform_parameters::stage(*stage) as u32),
            )?;
            let _ = source; // Retained through Error::source below.
        }
        E::Chirality { stage, source } => {
            set(
                &error,
                "stage",
                JsValue::from(crate::transform_parameters::stage(*stage) as u32),
            )?;
            let _ = source; // Retained through Error::source below.
        }
        E::AdjustHs { stage, source } => {
            set(
                &error,
                "stage",
                JsValue::from(crate::transform_parameters::stage(*stage) as u32),
            )?;
            let _ = source; // Retained through Error::source below.
        }
    }
    if let Some(cause) = source.source() {
        set(&error, "cause", source_error(cause)?)?;
    }
    Ok(error)
}
#[wasm_bindgen]
pub struct HydrogenError {
    inner: ck::HydrogenError,
}
#[wasm_bindgen]
impl HydrogenError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "hydrogens".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        hydrogen_error_kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
}
fn hydrogen_error_kind(source: &ck::HydrogenError) -> &'static str {
    use ck::HydrogenError as E;
    match source {
        E::StereoOrder(..) => "StereoOrder",
        E::BondProperty(..) => "BondProperty",
        E::AtomProperty(..) => "AtomProperty",
        E::MissingSourceState { .. } => "MissingSourceState",
        E::InvalidTopology(..) => "InvalidTopology",
        E::InvalidCoordinates(..) => "InvalidCoordinates",
        E::TopologyEdit(..) => "TopologyEdit",
        E::InvalidMapping(..) => "InvalidMapping",
        E::InvalidQueryState(..) => "InvalidQueryState",
        E::InvalidProperty(..) => "InvalidProperty",
        E::Valence(..) => "Valence",
        E::Sanitize(..) => "Sanitize",
        E::InvalidPropertyList { .. } => "InvalidPropertyList",
        E::Unsupported { .. } => "Unsupported",
        E::OnlyOnAtomOutOfRange { .. } => "OnlyOnAtomOutOfRange",
        E::InvalidAdditionPlan { .. } => "InvalidAdditionPlan",
        E::InvalidRemovalCandidate { .. } => "InvalidRemovalCandidate",
        E::ValenceAssignmentLength { .. } => "ValenceAssignmentLength",
        E::ExplicitHydrogenOverflow { .. } => "ExplicitHydrogenOverflow",
        E::InvalidStereoTransition { .. } => "InvalidStereoTransition",
        E::CoordinatePlacement { .. } => "CoordinatePlacement",
    }
}
pub(crate) fn hydrogen_error(source: &ck::HydrogenError) -> Result<JsValue, JsValue> {
    use ck::HydrogenError as E;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("HydrogenError");
    let error: JsValue = error.into();
    set(&error, "domain", "hydrogens".into())?;
    set(&error, "kind", hydrogen_error_kind(source).into())?;
    set(
        &error,
        "detail",
        HydrogenError {
            inner: source.clone(),
        }
        .into(),
    )?;
    match source {
        E::StereoOrder(..)
        | E::BondProperty(..)
        | E::AtomProperty(..)
        | E::MissingSourceState { .. } => {}
        E::InvalidTopology(_) => {}
        E::InvalidCoordinates(_) => {}
        E::TopologyEdit(_) => {}
        E::InvalidMapping(_) => {}
        E::InvalidQueryState(_) => {}
        E::InvalidProperty(_) => {}
        E::Valence(_) => {}
        E::Sanitize(_) => {}
        E::InvalidPropertyList {
            target,
            name,
            expected_rows,
            actual_rows,
        } => {
            set(
                &error,
                "target",
                JsValue::from_str(match target {
                    ck::SdfPropertyListTarget::Atom => "atom",
                    ck::SdfPropertyListTarget::Bond => "bond",
                }),
            )?;
            set(&error, "name", crate::host_values::text(name)?.into())?;
            set(
                &error,
                "expectedRows",
                JsValue::from_f64(*expected_rows as f64),
            )?;
            set(&error, "actualRows", JsValue::from_f64(*actual_rows as f64))?;
        }
        E::Unsupported { operation, reason } => {
            set(&error, "operation", JsValue::from_str(operation))?;
            set(&error, "reason", JsValue::from_str(reason))?;
        }
        E::OnlyOnAtomOutOfRange { atom, atom_count } => {
            set(&error, "atom", JsValue::from_f64(atom.index() as f64))?;
            set(&error, "atomCount", JsValue::from_f64(*atom_count as f64))?;
        }
        E::InvalidAdditionPlan { addition, reason } => {
            set(
                &error,
                "addition",
                addition.map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64)),
            )?;
            set(&error, "reason", JsValue::from_str(reason))?;
        }
        E::InvalidRemovalCandidate {
            position,
            atom,
            reason,
        } => {
            set(&error, "position", JsValue::from_f64(*position as f64))?;
            set(&error, "atom", JsValue::from_f64(atom.index() as f64))?;
            set(&error, "reason", JsValue::from_str(reason))?;
        }
        E::ValenceAssignmentLength {
            field,
            expected,
            actual,
        } => {
            set(&error, "field", JsValue::from_str(field))?;
            set(&error, "expected", JsValue::from_f64(*expected as f64))?;
            set(&error, "actual", JsValue::from_f64(*actual as f64))?;
        }
        E::ExplicitHydrogenOverflow { atom, current } => {
            set(&error, "atom", JsValue::from_f64(atom.index() as f64))?;
            set(&error, "current", JsValue::from_f64(*current as f64))?;
        }
        E::InvalidStereoTransition { bond, reason } => {
            set(&error, "bond", JsValue::from_f64(bond.index() as f64))?;
            set(&error, "reason", JsValue::from_str(reason))?;
        }
        E::CoordinatePlacement {
            addition,
            conformer,
            dimension,
            reason,
        } => {
            set(&error, "addition", JsValue::from_f64(*addition as f64))?;
            set(&error, "conformer", JsValue::from_f64(*conformer as f64))?;
            set(&error, "dimension", JsValue::from_str(dimension))?;
            set(&error, "reason", JsValue::from_str(reason))?;
        }
    }
    if let Some(cause) = source.source() {
        set(&error, "cause", source_error(cause)?)?;
    }
    Ok(error)
}
#[wasm_bindgen]
pub struct KekulizeError {
    inner: ck::KekulizeError,
}
#[wasm_bindgen]
impl KekulizeError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "kekulize".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        kekulize_error_kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
    #[wasm_bindgen(getter,js_name=expected,unchecked_return_type="number | null")]
    pub fn payload_0(&self) -> JsValue {
        use ck::KekulizeError as E;
        match &self.inner {
            E::AtomSelectionLength { expected, .. }
            | E::BondSelectionLength { expected, .. }
            | E::MatchingStateLength { expected, .. } => JsValue::from(*expected as u32),
            _ => JsValue::NULL,
        }
    }

    #[wasm_bindgen(getter,js_name=actual,unchecked_return_type="number | null")]
    pub fn payload_1(&self) -> JsValue {
        use ck::KekulizeError as E;
        match &self.inner {
            E::AtomSelectionLength { actual, .. }
            | E::BondSelectionLength { actual, .. }
            | E::MatchingStateLength { actual, .. } => JsValue::from(*actual as u32),
            _ => JsValue::NULL,
        }
    }

    #[wasm_bindgen(getter,js_name=atom,unchecked_return_type="number | null")]
    pub fn payload_2(&self) -> JsValue {
        use ck::KekulizeError as E;
        match &self.inner {
            E::CandidateAtomOutOfRange { atom, .. }
            | E::DuplicateCandidateAtom { atom }
            | E::DuplicateDoneAtom { atom }
            | E::MissingBacktrackAnchor { atom }
            | E::AromaticAtomOutsideRing { atom }
            | E::PostconditionValenceMismatch { atom, .. }
            | E::IntegerOverflow { atom, .. } => JsValue::from(atom.index() as u32),
            _ => JsValue::NULL,
        }
    }

    #[wasm_bindgen(getter,js_name=atomCount,unchecked_return_type="number | null")]
    pub fn payload_3(&self) -> JsValue {
        use ck::KekulizeError as E;
        match &self.inner {
            E::CandidateAtomOutOfRange { atom_count, .. } => JsValue::from(*atom_count as u32),
            _ => JsValue::NULL,
        }
    }

    #[wasm_bindgen(getter,js_name=field,unchecked_return_type="string | null")]
    pub fn payload_4(&self) -> JsValue {
        use ck::KekulizeError as E;
        match &self.inner {
            E::MatchingStateLength { field, .. } | E::IntegerOverflow { field, .. } => {
                JsValue::from(*field)
            }
            _ => JsValue::NULL,
        }
    }

    #[wasm_bindgen(getter,js_name=begin,unchecked_return_type="number | null")]
    pub fn payload_5(&self) -> JsValue {
        use ck::KekulizeError as E;
        match &self.inner {
            E::MissingCandidateBond { begin, .. } => JsValue::from(begin.index() as u32),
            _ => JsValue::NULL,
        }
    }

    #[wasm_bindgen(getter,js_name=end,unchecked_return_type="number | null")]
    pub fn payload_6(&self) -> JsValue {
        use ck::KekulizeError as E;
        match &self.inner {
            E::MissingCandidateBond { end, .. } => JsValue::from(end.index() as u32),
            _ => JsValue::NULL,
        }
    }

    #[wasm_bindgen(getter,js_name=questions,unchecked_return_type="number | null")]
    pub fn payload_7(&self) -> JsValue {
        use ck::KekulizeError as E;
        match &self.inner {
            E::QuestionSubsetOverflow { questions, .. } => JsValue::from(*questions as u32),
            _ => JsValue::NULL,
        }
    }

    #[wasm_bindgen(getter,js_name=bitWidth,unchecked_return_type="number | null")]
    pub fn payload_8(&self) -> JsValue {
        use ck::KekulizeError as E;
        match &self.inner {
            E::QuestionSubsetOverflow { bit_width, .. } => JsValue::from(*bit_width),
            _ => JsValue::NULL,
        }
    }

    #[wasm_bindgen(getter,js_name=problemAtoms,unchecked_return_type="number[] | null")]
    pub fn payload_9(&self) -> JsValue {
        use ck::KekulizeError as E;
        match &self.inner {
            E::NotKekulizable { problem_atoms } => problem_atoms
                .iter()
                .map(|id| JsValue::from(id.index() as u32))
                .collect::<Array>()
                .into(),
            _ => JsValue::NULL,
        }
    }

    #[wasm_bindgen(getter,js_name=before,unchecked_return_type="number | null")]
    pub fn payload_10(&self) -> JsValue {
        use ck::KekulizeError as E;
        match &self.inner {
            E::PostconditionValenceMismatch { before, .. } => JsValue::from(*before),
            _ => JsValue::NULL,
        }
    }

    #[wasm_bindgen(getter,js_name=after,unchecked_return_type="number | null")]
    pub fn payload_11(&self) -> JsValue {
        use ck::KekulizeError as E;
        match &self.inner {
            E::PostconditionValenceMismatch { after, .. } => JsValue::from(*after),
            _ => JsValue::NULL,
        }
    }

    #[wasm_bindgen(getter,js_name=bond,unchecked_return_type="number | null")]
    pub fn payload_12(&self) -> JsValue {
        use ck::KekulizeError as E;
        match &self.inner {
            E::UnsupportedQueryState { bond, .. } => JsValue::from(bond.index() as u32),
            _ => JsValue::NULL,
        }
    }

    #[wasm_bindgen(getter,js_name=reason,unchecked_return_type="string | null")]
    pub fn payload_13(&self) -> JsValue {
        use ck::KekulizeError as E;
        match &self.inner {
            E::UnsupportedQueryState { detail, .. } => JsValue::from(*detail),
            _ => JsValue::NULL,
        }
    }
}
fn kekulize_error_kind(source: &ck::KekulizeError) -> &'static str {
    use ck::KekulizeError as E;
    match source {
        E::InvalidTopology(..) => "InvalidTopology",
        E::InvalidQueryState(..) => "InvalidQueryState",
        E::RingFinding(..) => "RingFinding",
        E::Valence(..) => "Valence",
        E::CanonicalRank(..) => "CanonicalRank",
        E::AtomSelectionLength { .. } => "AtomSelectionLength",
        E::BondSelectionLength { .. } => "BondSelectionLength",
        E::CandidateAtomOutOfRange { .. } => "CandidateAtomOutOfRange",
        E::DuplicateCandidateAtom { .. } => "DuplicateCandidateAtom",
        E::DuplicateDoneAtom { .. } => "DuplicateDoneAtom",
        E::MissingBacktrackAnchor { .. } => "MissingBacktrackAnchor",
        E::AromaticAtomOutsideRing { .. } => "AromaticAtomOutsideRing",
        E::MatchingStateLength { .. } => "MatchingStateLength",
        E::MissingCandidateBond { .. } => "MissingCandidateBond",
        E::QuestionSubsetOverflow { .. } => "QuestionSubsetOverflow",
        E::NotKekulizable { .. } => "NotKekulizable",
        E::PostconditionValenceMismatch { .. } => "PostconditionValenceMismatch",
        E::UnsupportedQueryState { .. } => "UnsupportedQueryState",
        E::IntegerOverflow { .. } => "IntegerOverflow",
    }
}
pub(crate) fn kekulize_error(source: &ck::KekulizeError) -> Result<JsValue, JsValue> {
    use ck::KekulizeError as E;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("KekulizeError");
    let error: JsValue = error.into();
    set(&error, "domain", "kekulize".into())?;
    set(&error, "kind", kekulize_error_kind(source).into())?;
    set(
        &error,
        "detail",
        KekulizeError {
            inner: source.clone(),
        }
        .into(),
    )?;
    match source {
        E::InvalidTopology(_) => {}
        E::InvalidQueryState(_) => {}
        E::RingFinding(_) => {}
        E::Valence(_) => {}
        E::CanonicalRank(_) => {}
        E::AtomSelectionLength { expected, actual } => {
            set(&error, "expected", JsValue::from_f64(*expected as f64))?;
            set(&error, "actual", JsValue::from_f64(*actual as f64))?;
        }
        E::BondSelectionLength { expected, actual } => {
            set(&error, "expected", JsValue::from_f64(*expected as f64))?;
            set(&error, "actual", JsValue::from_f64(*actual as f64))?;
        }
        E::CandidateAtomOutOfRange { atom, atom_count } => {
            set(&error, "atom", JsValue::from_f64(atom.index() as f64))?;
            set(&error, "atomCount", JsValue::from_f64(*atom_count as f64))?;
        }
        E::DuplicateCandidateAtom { atom } => {
            set(&error, "atom", JsValue::from_f64(atom.index() as f64))?;
        }
        E::DuplicateDoneAtom { atom } => {
            set(&error, "atom", JsValue::from_f64(atom.index() as f64))?;
        }
        E::MissingBacktrackAnchor { atom } => {
            set(&error, "atom", JsValue::from_f64(atom.index() as f64))?;
        }
        E::AromaticAtomOutsideRing { atom } => {
            set(&error, "atom", JsValue::from_f64(atom.index() as f64))?;
        }
        E::MatchingStateLength {
            field,
            expected,
            actual,
        } => {
            set(&error, "field", JsValue::from_str(field))?;
            set(&error, "expected", JsValue::from_f64(*expected as f64))?;
            set(&error, "actual", JsValue::from_f64(*actual as f64))?;
        }
        E::MissingCandidateBond { begin, end } => {
            set(&error, "begin", JsValue::from_f64(begin.index() as f64))?;
            set(&error, "end", JsValue::from_f64(end.index() as f64))?;
        }
        E::QuestionSubsetOverflow {
            questions,
            bit_width,
        } => {
            set(&error, "questions", JsValue::from_f64(*questions as f64))?;
            set(&error, "bitWidth", JsValue::from_f64(*bit_width as f64))?;
        }
        E::NotKekulizable { problem_atoms } => {
            set(
                &error,
                "problemAtoms",
                problem_atoms
                    .iter()
                    .map(|id| JsValue::from_f64(id.index() as f64))
                    .collect::<Array>()
                    .into(),
            )?;
        }
        E::PostconditionValenceMismatch {
            atom,
            before,
            after,
        } => {
            set(&error, "atom", JsValue::from_f64(atom.index() as f64))?;
            set(&error, "before", JsValue::from_f64(*before as f64))?;
            set(&error, "after", JsValue::from_f64(*after as f64))?;
        }
        E::UnsupportedQueryState { bond, detail } => {
            set(&error, "bond", JsValue::from_f64(bond.index() as f64))?;
            set(&error, "reason", JsValue::from_str(detail))?;
        }
        E::IntegerOverflow { atom, field } => {
            set(&error, "atom", JsValue::from_f64(atom.index() as f64))?;
            set(&error, "field", JsValue::from_str(field))?;
        }
    }
    if let Some(cause) = source.source() {
        set(&error, "cause", source_error(cause)?)?;
    }
    Ok(error)
}
#[wasm_bindgen]
pub struct Coordinate2DError {
    inner: ck::Coordinate2DError,
}
#[wasm_bindgen]
impl Coordinate2DError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "depict".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        coordinate_2d_error_kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
}
fn coordinate_2d_error_kind(source: &ck::Coordinate2DError) -> &'static str {
    use ck::Coordinate2DError as E;
    match source {
        E::InvalidTopology(..) => "InvalidTopology",
        E::PropertyCache(..) => "PropertyCache",
        E::RingFinding(..) => "RingFinding",
        E::StereoAssignment(..) => "StereoAssignment",
        E::TemplateLoading(..) => "TemplateLoading",
        E::Fragment(..) => "Fragment",
        E::CoordinateValidation(..) => "CoordinateValidation",
        E::ConformerIdOverflow { .. } => "ConformerIdOverflow",
        E::CoordGenUnavailable => "CoordGenUnavailable",
    }
}
pub(crate) fn coordinate_2d_error(source: &ck::Coordinate2DError) -> Result<JsValue, JsValue> {
    use ck::Coordinate2DError as E;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("Coordinate2DError");
    let error: JsValue = error.into();
    set(&error, "domain", "depict".into())?;
    set(&error, "kind", coordinate_2d_error_kind(source).into())?;
    set(
        &error,
        "detail",
        Coordinate2DError {
            inner: source.clone(),
        }
        .into(),
    )?;
    match source {
        E::InvalidTopology(_) => {}
        E::PropertyCache(_) => {}
        E::RingFinding(_) => {}
        E::StereoAssignment(_) => {}
        E::TemplateLoading(_) => {}
        E::Fragment(_) => {}
        E::CoordinateValidation(_) => {}
        E::ConformerIdOverflow { max_id } => {
            set(&error, "maxId", JsValue::from_f64(*max_id as f64))?;
        }
        E::CoordGenUnavailable => {}
    }
    if let Some(cause) = source.source() {
        set(&error, "cause", source_error(cause)?)?;
    }
    Ok(error)
}
