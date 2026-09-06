use super::{
    ConfSeqDecodeError, ConfSeqDecodeOptions, ConfSeqFastGeometryError, ConfSeqTemplateBackend,
    build_confseq_base_template, build_distance_geometry_template, decode_from_template,
    parse_confseq, prepare_p_chiral_embedding_molecule,
};
use crate::{ChiralTag, Molecule};
use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq)]
pub struct ParsedConfSeqDiagnostics {
    pub stripped_smiles: String,
    pub chiral_tags_by_atom: HashMap<usize, ChiralTag>,
    pub raw_dihedral_literals_by_pair: HashMap<(usize, usize), String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DihedralPairDiagnostics {
    pub parsed_pairs: Vec<(usize, usize)>,
    pub template_pairs: Vec<(usize, usize)>,
    pub missing_from_tokens: Vec<(usize, usize)>,
    pub extra_token_pairs: Vec<(usize, usize)>,
    pub raw_dihedral_literals_by_pair: HashMap<(usize, usize), String>,
}

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
    pub base_error: Option<ConfSeqFastGeometryError>,
    pub parsed: bool,
    pub distance_geometry_template_built: bool,
    pub base_template_built: bool,
    pub distance_geometry_decoded: bool,
    pub base_decoded: bool,
}

pub fn parse_confseq_for_diagnostics(
    in_smiles: &str,
    confseq: &str,
) -> Result<ParsedConfSeqDiagnostics, ConfSeqDecodeError> {
    let parsed = parse_confseq(in_smiles, confseq)?;
    Ok(ParsedConfSeqDiagnostics {
        stripped_smiles: parsed.stripped_smiles,
        chiral_tags_by_atom: parsed.chiral_tags_by_atom,
        raw_dihedral_literals_by_pair: parsed.raw_dihedral_literals_by_pair,
    })
}

pub fn prepare_p_chiral_for_diagnostics(
    smiles: &str,
    chiral_tags_by_atom: &HashMap<usize, ChiralTag>,
) -> Result<Molecule, ConfSeqDecodeError> {
    let molecule = Molecule::from_smiles(smiles)
        .map_err(|err| ConfSeqDecodeError::SmilesParse(err.to_string()))?;
    prepare_p_chiral_embedding_molecule(molecule, chiral_tags_by_atom)
}

pub fn distance_geometry_dihedral_pair_diagnostics(
    in_smiles: &str,
    confseq: &str,
    options: &ConfSeqDecodeOptions,
) -> Result<DihedralPairDiagnostics, ConfSeqDecodeError> {
    let parsed = parse_confseq(in_smiles, confseq)?;
    let mut dg_options = options.clone();
    dg_options.template_backend = ConfSeqTemplateBackend::DistanceGeometry;
    let template = build_distance_geometry_template(
        &parsed.stripped_smiles,
        &parsed.chiral_tags_by_atom,
        &dg_options,
    )?;
    let mut parsed_pairs = parsed
        .dihedral_angles_by_pair
        .keys()
        .copied()
        .collect::<Vec<_>>();
    let mut template_pairs = template
        .dihedrals_by_pair
        .keys()
        .copied()
        .collect::<Vec<_>>();
    parsed_pairs.sort_unstable();
    template_pairs.sort_unstable();
    let parsed_set = parsed_pairs
        .iter()
        .copied()
        .collect::<std::collections::HashSet<_>>();
    let template_set = template_pairs
        .iter()
        .copied()
        .collect::<std::collections::HashSet<_>>();
    let mut missing_from_tokens = template_pairs
        .iter()
        .copied()
        .filter(|pair| !parsed_set.contains(pair))
        .collect::<Vec<_>>();
    let mut extra_token_pairs = parsed_pairs
        .iter()
        .copied()
        .filter(|pair| !template_set.contains(pair))
        .collect::<Vec<_>>();
    missing_from_tokens.sort_unstable();
    extra_token_pairs.sort_unstable();
    Ok(DihedralPairDiagnostics {
        parsed_pairs,
        template_pairs,
        missing_from_tokens,
        extra_token_pairs,
        raw_dihedral_literals_by_pair: parsed.raw_dihedral_literals_by_pair,
    })
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
            ConfSeqDecodeError::FastGeometry(error) => Some(error.clone()),
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
    base_options.template_backend = ConfSeqTemplateBackend::FastGeometry;
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
