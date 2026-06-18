use super::{
    ConfSeqBaseConformerError, ConfSeqDecodeError, ConfSeqDecodeOptions, ConfSeqTemplateBackend,
    build_confseq_base_template, build_distance_geometry_template, decode_from_template,
    parse_confseq,
};

/// Coarse phase where a ConfSeq decode attempt failed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConfSeqDiagnosticPhase {
    Parse,
    DistanceGeometryTemplate,
    BaseTemplate,
    Decode,
}

/// Diagnostic outcome for one ConfSeq candidate.
#[derive(Debug, Clone, PartialEq)]
pub struct ConfSeqDiagnostic {
    pub phase: Option<ConfSeqDiagnosticPhase>,
    pub error: Option<ConfSeqDecodeError>,
    pub base_error: Option<ConfSeqBaseConformerError>,
    pub parsed: bool,
    pub distance_geometry_template_built: bool,
    pub base_template_built: bool,
    pub distance_geometry_decoded: bool,
    pub base_decoded: bool,
}

impl ConfSeqDiagnostic {
    fn success() -> Self {
        Self {
            phase: None,
            error: None,
            base_error: None,
            parsed: true,
            distance_geometry_template_built: true,
            base_template_built: true,
            distance_geometry_decoded: true,
            base_decoded: true,
        }
    }

    fn failed(
        phase: ConfSeqDiagnosticPhase,
        error: ConfSeqDecodeError,
        parsed: bool,
        distance_geometry_template_built: bool,
        base_template_built: bool,
        distance_geometry_decoded: bool,
        base_decoded: bool,
    ) -> Self {
        let base_error = match &error {
            ConfSeqDecodeError::BaseConformer(error) => Some(error.clone()),
            _ => None,
        };
        Self {
            phase: Some(phase),
            error: Some(error),
            base_error,
            parsed,
            distance_geometry_template_built,
            base_template_built,
            distance_geometry_decoded,
            base_decoded,
        }
    }
}

/// Run both ConfSeq template backends far enough to classify failures.
///
/// This is intended for corpus development and benchmark scripts. It does not
/// align conformers or compute RMSD; callers can do that from the normal decode
/// API after selecting successful candidates.
pub fn diagnose_confseq_candidate(
    in_smiles: &str,
    confseq: &str,
    options: &ConfSeqDecodeOptions,
) -> ConfSeqDiagnostic {
    let parsed = match parse_confseq(in_smiles, confseq) {
        Ok(parsed) => parsed,
        Err(error) => {
            return ConfSeqDiagnostic::failed(
                ConfSeqDiagnosticPhase::Parse,
                error,
                false,
                false,
                false,
                false,
                false,
            );
        }
    };

    let mut dg_options = options.clone();
    dg_options.template_backend = ConfSeqTemplateBackend::DistanceGeometry;
    let dg_template = match build_distance_geometry_template(
        &parsed.stripped_smiles,
        &parsed.chiral_tags_by_atom,
        &dg_options,
    ) {
        Ok(template) => template,
        Err(error) => {
            return ConfSeqDiagnostic::failed(
                ConfSeqDiagnosticPhase::DistanceGeometryTemplate,
                error,
                true,
                false,
                false,
                false,
                false,
            );
        }
    };

    let dg_decoded = decode_from_template(&dg_template, &parsed, &dg_options).is_ok();

    let base_template =
        match build_confseq_base_template(&parsed.stripped_smiles, &parsed.chiral_tags_by_atom) {
            Ok(template) => template,
            Err(error) => {
                return ConfSeqDiagnostic::failed(
                    ConfSeqDiagnosticPhase::BaseTemplate,
                    error,
                    true,
                    true,
                    false,
                    dg_decoded,
                    false,
                );
            }
        };

    let mut base_options = options.clone();
    base_options.template_backend = ConfSeqTemplateBackend::BaseConformer;
    if let Err(error) = decode_from_template(&base_template, &parsed, &base_options) {
        return ConfSeqDiagnostic::failed(
            ConfSeqDiagnosticPhase::Decode,
            error,
            true,
            true,
            true,
            dg_decoded,
            false,
        );
    }

    let mut diagnostic = ConfSeqDiagnostic::success();
    diagnostic.distance_geometry_decoded = dg_decoded;
    diagnostic
}

pub fn diagnose_confseq_batch(
    in_smiles: &[String],
    confseq: &[String],
    options: &ConfSeqDecodeOptions,
) -> Result<Vec<ConfSeqDiagnostic>, ConfSeqDecodeError> {
    if in_smiles.len() != confseq.len() {
        return Err(ConfSeqDecodeError::BatchLengthMismatch {
            in_smiles: in_smiles.len(),
            confseq: confseq.len(),
        });
    }
    if matches!(options.num_threads, Some(0)) {
        return Err(ConfSeqDecodeError::InvalidThreadCount);
    }
    Ok(in_smiles
        .iter()
        .zip(confseq)
        .map(|(smiles, confseq)| diagnose_confseq_candidate(smiles, confseq, options))
        .collect())
}
