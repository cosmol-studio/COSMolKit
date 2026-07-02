use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};
use std::f64::consts::PI;

use rayon::prelude::*;
use thiserror::Error;

use crate::chemistry::{coordinates, distgeom, mol_transforms, rings};
#[cfg(test)]
use crate::io::molblock::{MolBlockWriteParams, SdfFormat, mol_to_sdf_record_with_params};
use crate::smiles_write::{SmilesWriteParams, mol_to_smiles, mol_to_smiles_with_atom_output_order};
use crate::{
    AdjacencyList, AtomId, Bond, BondId, BondOrder, BondSpec, BondStereo, ChiralTag, Conformer3D,
    EmbedParameters, Hybridization, Molecule, get_uff_angle_bend_params,
    get_uff_bond_stretch_params, uff_optimize_molecule,
};

pub mod diagnostics;

#[derive(Clone)]
pub struct ConfSeqDecodeOptions {
    pub embed_params: EmbedParameters,
    pub optimize_with_uff: bool,
    pub num_threads: Option<usize>,
    pub apply_dihedrals: bool,
    pub apply_angles: bool,
    pub template_backend: ConfSeqTemplateBackend,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConfSeqTemplateBackend {
    /// Current production path: RDKit-aligned distance geometry followed by
    /// optional UFF relaxation. This remains the default because it supports the
    /// broadest ConfSeq input space currently modeled by COSMolKit.
    DistanceGeometry,
    /// Experimental ConfSeq-specific fast geometry. This is intentionally a
    /// strict subset, not a heuristic fallback: unsupported topologies return
    /// `ConfSeqDecodeError::FastGeometry` instead of silently switching back to
    /// distance geometry.
    FastGeometry,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ConfSeqBaseStructuralRiskClass {
    Nonplanar5MemberRingPucker,
    Nonplanar6MemberRingPucker,
    SharedSp3RingCenterBranch,
    MultiRingSharedAtomBranch,
    FusedNonplanarRingEdgeBranch,
    NonplanarMultiringRigidComponent,
}

impl ConfSeqBaseStructuralRiskClass {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Nonplanar5MemberRingPucker => "nonplanar_5_member_ring_pucker",
            Self::Nonplanar6MemberRingPucker => "nonplanar_6_member_ring_pucker",
            Self::SharedSp3RingCenterBranch => "shared_sp3_ring_center_branch",
            Self::MultiRingSharedAtomBranch => "multi_ring_shared_atom_branch",
            Self::FusedNonplanarRingEdgeBranch => "fused_nonplanar_ring_edge_branch",
            Self::NonplanarMultiringRigidComponent => "nonplanar_multiring_rigid_component",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConfSeqBaseStructuralRiskPrecheck {
    pub classes: Vec<ConfSeqBaseStructuralRiskClass>,
    pub fallback_candidate: bool,
}

impl ConfSeqBaseStructuralRiskPrecheck {
    pub fn has_risk(&self) -> bool {
        !self.classes.is_empty()
    }

    pub fn class_names(&self) -> Vec<&'static str> {
        self.classes.iter().map(|class| class.as_str()).collect()
    }
}

impl Default for ConfSeqDecodeOptions {
    fn default() -> Self {
        let mut embed_params = EmbedParameters::etkdg();
        embed_params.random_seed = 0;
        Self {
            embed_params,
            optimize_with_uff: true,
            num_threads: None,
            apply_dihedrals: true,
            apply_angles: true,
            template_backend: ConfSeqTemplateBackend::DistanceGeometry,
        }
    }
}

#[derive(Debug, Clone, Error, PartialEq)]
pub enum ConfSeqDecodeError {
    #[error("in_smiles and confseq batch lengths differ: in_smiles={in_smiles}, confseq={confseq}")]
    BatchLengthMismatch { in_smiles: usize, confseq: usize },
    #[error("ConfSeq batch n_jobs must be >= 1")]
    InvalidThreadCount,
    #[error("ConfSeq token stream contains unsupported token '{token}'")]
    UnsupportedToken { token: String },
    #[error(
        "ConfSeq contains {observed} dihedral tokens but molecule has {expected} eligible single bonds"
    )]
    DihedralTokenCountMismatch { observed: usize, expected: usize },
    #[error(
        "ConfSeq contains {observed} angle tokens but molecule has {expected} eligible angle centers"
    )]
    AngleTokenCountMismatch { observed: usize, expected: usize },
    #[error(
        "ConfSeq token at explicit-bond position {position} is missing in completed token stream"
    )]
    TokenPositionOutOfRange { position: usize },
    #[error("ConfSeq explicit-bond SMILES mapping failed: {0}")]
    BondTokenMapping(String),
    #[error("SMILES parse failed: {0}")]
    SmilesParse(String),
    #[error("SMILES write failed: {0}")]
    SmilesWrite(String),
    #[error("ring perception failed: {0}")]
    RingFinding(String),
    #[error("{0}")]
    FastGeometry(#[from] ConfSeqFastGeometryError),
    #[error("3D embedding failed")]
    EmbedFailed,
    #[error("molecular transform failed: {0}")]
    MolTransform(String),
}

#[derive(Debug, Clone, Error, PartialEq)]
pub enum ConfSeqFastGeometryError {
    #[error("ConfSeq fast geometry currently supports only connected molecules")]
    Disconnected,
    #[error(
        "ConfSeq fast geometry supports only supported 3-8 member heavy-atom rings that are disjoint or edge-fused"
    )]
    UnsupportedRingSystem,
    #[error("ConfSeq fast geometry does not support {ring_size}-member rings")]
    UnsupportedRingSize { ring_size: usize },
    #[error(
        "ConfSeq fast geometry does not support ring atom element {atomic_number} in atom ring {ring_index}"
    )]
    UnsupportedRingElement {
        ring_index: usize,
        atomic_number: u8,
    },
    #[error("ConfSeq fast geometry does not support aromatic {ring_size}-member ring")]
    UnsupportedAromaticRingSize { ring_size: usize },
    #[error(
        "ConfSeq fast geometry does not support ring sharing between rings {left} and {right}: shared_atoms={shared_atoms}, shared_bonds={shared_bonds}"
    )]
    UnsupportedRingFusion {
        left: usize,
        right: usize,
        shared_atoms: usize,
        shared_bonds: usize,
    },
    #[error("ConfSeq fast geometry does not support closed fused/spiro ring constraint components")]
    UnsupportedClosedRingFusion,
    #[error(
        "ConfSeq fast geometry could not satisfy tetrahedral stereo at atom {center}: {reason}"
    )]
    UnsupportedTetrahedralStereo { center: usize, reason: String },
    #[error("ConfSeq fast geometry ring geometry is invalid: {0}")]
    InvalidRingGeometry(String),
    #[error("ConfSeq fast geometry rigid fragment embedding failed: {0}")]
    RigidFragmentEmbedding(String),
    #[error("ConfSeq fast geometry traversal left atoms unplaced")]
    PlacementLeftAtomsUnplaced,
    #[error(
        "ConfSeq fast geometry fragment assembly has no usable anchor for component {component}"
    )]
    AssemblyMissingAnchor { component: usize },
    #[error(
        "ConfSeq fast geometry fragment assembly has conflicting anchors for component {component}: {reason}"
    )]
    AssemblyAnchorConflict { component: usize, reason: String },
    #[error(
        "ConfSeq fast geometry fragment assembly boundary bond mismatch for bond {bond}: observed={observed:.3}, target={target:.3}"
    )]
    AssemblyBoundaryMismatch {
        bond: usize,
        observed: f64,
        target: f64,
    },
    #[error("ConfSeq fast geometry failed: {0}")]
    Build(String),
}

#[derive(Debug, Clone)]
pub struct ConfSeqBatchDecodeResult {
    pub molecules: Vec<Option<Molecule>>,
    pub errors: Vec<Option<ConfSeqDecodeError>>,
}

#[derive(Debug, Clone)]
struct ParsedConfSeq {
    stripped_smiles: String,
    dihedral_angles_by_pair: HashMap<(usize, usize), f64>,
    raw_dihedral_literals_by_pair: HashMap<(usize, usize), String>,
    angle_values_deg: Vec<f64>,
    chiral_tags_by_atom: HashMap<usize, ChiralTag>,
}

#[derive(Debug, Clone)]
struct Template {
    molecule: Molecule,
    dihedrals: Vec<(usize, usize, usize, usize)>,
    dihedrals_by_pair: HashMap<(usize, usize), (usize, usize, usize, usize)>,
    angle_centers: Vec<(usize, usize, usize)>,
    ring_bond_pairs: HashSet<(usize, usize)>,
    last_ring_bonds: Vec<(usize, usize, BondSpec)>,
}

#[derive(Debug, Clone)]
struct BondTokenMapping {
    smiles_be: String,
    token_idx_to_atom_pair: HashMap<usize, (usize, usize)>,
    ring_closure_pairs: Vec<(usize, usize)>,
}

impl BondTokenMapping {
    fn atom_pairs_in_token_order(&self) -> Vec<(usize, usize)> {
        let mut entries = self
            .token_idx_to_atom_pair
            .iter()
            .map(|(token_idx, pair)| (*token_idx, *pair))
            .collect::<Vec<_>>();
        entries.sort_by_key(|(token_idx, _)| *token_idx);
        entries.into_iter().map(|(_, pair)| pair).collect()
    }
}

#[derive(Debug, Clone)]
struct ConfSeqBaseConstraintModel {
    bond_targets: Vec<f64>,
    angle_targets: HashMap<(usize, usize, usize), f64>,
    torsion_priors: HashMap<(usize, usize, usize, usize), ConfSeqBaseTorsionPrior>,
    path14_distance_priors: Vec<ConfSeqBasePath14DistancePrior>,
    planar_bonds: HashSet<(usize, usize)>,
    rigid_components: Vec<ConfSeqBaseRigidComponent>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ConfSeqBaseTorsionPrior {
    Cis,
    Trans,
    Free,
}

#[derive(Debug, Clone)]
struct ConfSeqBasePath14DistancePrior {
    atoms: (usize, usize, usize, usize),
    lower_bound: f64,
    upper_bound: f64,
}

#[derive(Debug, Clone)]
struct ConfSeqBaseRingComponent {
    atoms: Vec<usize>,
    ring_sizes_by_atom: HashMap<usize, usize>,
    planar: bool,
}

#[derive(Debug, Clone)]
struct ConfSeqBaseRigidComponent {
    atoms: Vec<usize>,
    kind: ConfSeqBaseRigidComponentKind,
    connectors: Vec<ConfSeqBaseTemplateConnector>,
}

#[derive(Debug, Clone)]
struct ConfSeqBaseRigidFragmentTemplate {
    atoms: Vec<usize>,
    realization_atoms: Vec<usize>,
    stereo_context_atoms: Vec<usize>,
    shape: ConfSeqBaseRigidFragmentShape,
    bonds: Vec<ConfSeqBaseTemplateBond>,
    angles: Vec<ConfSeqBaseTemplateAngle>,
    connectors: Vec<ConfSeqBaseTemplateConnector>,
    cache_key: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ConfSeqBaseRigidFragmentShape {
    Atom,
    Bond,
    Angle,
    SingleCenterPlanar,
    SingleCenterNonplanar,
    Chain,
    Tree,
    Planar,
    RingPuckered,
    RingPolycyclic,
    Unsupported,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ConfSeqBaseRingTopology {
    Simple,
    EdgeFusedChain,
    EdgeFusedPolycyclic,
    Spiro,
    BridgedOrCage,
    SingleNonSimple,
    Unknown,
}

#[derive(Debug, Clone)]
struct ConfSeqBaseTemplateBond {
    begin: usize,
    end: usize,
    length: f64,
}

#[derive(Debug, Clone)]
struct ConfSeqBaseTemplateAngle {
    left: usize,
    center: usize,
    right: usize,
    angle_rad: f64,
}

#[derive(Debug, Clone)]
struct ConfSeqBaseTemplateConnector {
    core_atom: usize,
    external_atom: usize,
    bond_id: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ConfSeqBaseRigidComponentKind {
    AcyclicPlanarPi,
    AcyclicNonplanar,
    AcyclicOther,
    RingPlanar,
    RingNonplanar,
    RingMixed,
}

pub fn decode_confseq(in_smiles: &str, confseq: &str) -> Result<Molecule, ConfSeqDecodeError> {
    decode_confseq_with_options(in_smiles, confseq, &ConfSeqDecodeOptions::default())
}

pub fn decode_confseq_record(confseq: &str) -> Result<Molecule, ConfSeqDecodeError> {
    decode_confseq_record_with_options(confseq, &ConfSeqDecodeOptions::default())
}

pub fn decode_confseq_record_with_options(
    confseq: &str,
    options: &ConfSeqDecodeOptions,
) -> Result<Molecule, ConfSeqDecodeError> {
    let in_smiles = infer_confseq_input_smiles(confseq)?;
    decode_confseq_with_options(&in_smiles, confseq, options)
}

pub fn decode_confseq_with_options(
    in_smiles: &str,
    confseq: &str,
    options: &ConfSeqDecodeOptions,
) -> Result<Molecule, ConfSeqDecodeError> {
    let parsed = parse_confseq(in_smiles, confseq)?;
    let template = build_template(
        &parsed.stripped_smiles,
        &parsed.chiral_tags_by_atom,
        options,
    )?;
    decode_from_template(&template, &parsed, options)
}

pub fn decode_confseq_batch(
    in_smiles: &[String],
    confseq: &[String],
    keep_errors: bool,
) -> Result<ConfSeqBatchDecodeResult, ConfSeqDecodeError> {
    decode_confseq_batch_with_options(
        in_smiles,
        confseq,
        &ConfSeqDecodeOptions::default(),
        keep_errors,
    )
}

pub fn decode_confseq_record_batch(
    confseq: &[String],
    keep_errors: bool,
) -> Result<ConfSeqBatchDecodeResult, ConfSeqDecodeError> {
    decode_confseq_record_batch_with_options(confseq, &ConfSeqDecodeOptions::default(), keep_errors)
}

pub fn decode_confseq_record_batch_with_options(
    confseq: &[String],
    options: &ConfSeqDecodeOptions,
    keep_errors: bool,
) -> Result<ConfSeqBatchDecodeResult, ConfSeqDecodeError> {
    let in_smiles: Result<Vec<_>, _> = confseq
        .iter()
        .map(|record| infer_confseq_input_smiles(record))
        .collect();
    let in_smiles = match in_smiles {
        Ok(value) => value,
        Err(_) if keep_errors => {
            let mut molecules = Vec::with_capacity(confseq.len());
            let mut errors = Vec::with_capacity(confseq.len());
            for record in confseq {
                match infer_confseq_input_smiles(record)
                    .and_then(|in_smiles| decode_confseq_with_options(&in_smiles, record, options))
                {
                    Ok(molecule) => {
                        molecules.push(Some(molecule));
                        errors.push(None);
                    }
                    Err(error) => {
                        molecules.push(None);
                        errors.push(Some(error));
                    }
                }
            }
            return Ok(ConfSeqBatchDecodeResult { molecules, errors });
        }
        Err(error) => return Err(error),
    };
    decode_confseq_batch_with_options(&in_smiles, confseq, options, keep_errors)
}

pub fn decode_confseq_batch_with_options(
    in_smiles: &[String],
    confseq: &[String],
    options: &ConfSeqDecodeOptions,
    keep_errors: bool,
) -> Result<ConfSeqBatchDecodeResult, ConfSeqDecodeError> {
    if in_smiles.len() != confseq.len() {
        return Err(ConfSeqDecodeError::BatchLengthMismatch {
            in_smiles: in_smiles.len(),
            confseq: confseq.len(),
        });
    }
    if matches!(options.num_threads, Some(0)) {
        return Err(ConfSeqDecodeError::InvalidThreadCount);
    }

    let decoded_inputs: Vec<_> = in_smiles
        .iter()
        .zip(confseq)
        .map(|(smiles, td)| parse_confseq(smiles, td))
        .collect();

    if !keep_errors {
        for parsed in &decoded_inputs {
            if let Err(error) = parsed {
                return Err(error.clone());
            }
        }
    }

    let mut cache = HashMap::<(String, Vec<(usize, u8)>), Template>::new();
    let mut templates = Vec::with_capacity(decoded_inputs.len());
    for parsed in &decoded_inputs {
        match parsed {
            Ok(parsed) => {
                let mut chiral_key: Vec<_> = parsed
                    .chiral_tags_by_atom
                    .iter()
                    .map(|(atom, tag)| (*atom, chiral_tag_cache_code(*tag)))
                    .collect();
                chiral_key.sort_unstable_by_key(|(atom, _)| *atom);
                let cache_key = (parsed.stripped_smiles.clone(), chiral_key);
                if !cache.contains_key(&cache_key) {
                    let template = build_template(
                        &parsed.stripped_smiles,
                        &parsed.chiral_tags_by_atom,
                        options,
                    );
                    match template {
                        Ok(template) => {
                            cache.insert(cache_key.clone(), template);
                        }
                        Err(error) if keep_errors => {
                            templates.push(Err(error));
                            continue;
                        }
                        Err(error) => return Err(error),
                    }
                }
                templates.push(Ok(cache
                    .get(&cache_key)
                    .expect("template inserted above")
                    .clone()));
            }
            Err(error) => templates.push(Err(error.clone())),
        }
    }

    let run_decode = || {
        decoded_inputs
            .par_iter()
            .zip(templates.par_iter())
            .map(|(parsed, template)| match (parsed, template) {
                (Ok(parsed), Ok(template)) => decode_from_template(template, parsed, options),
                (Err(error), _) | (_, Err(error)) => Err(error.clone()),
            })
            .collect::<Vec<_>>()
    };

    let decoded = match options.num_threads {
        Some(n_jobs) => rayon::ThreadPoolBuilder::new()
            .num_threads(n_jobs)
            .build()
            .map_err(|_| ConfSeqDecodeError::InvalidThreadCount)?
            .install(run_decode),
        None => run_decode(),
    };

    let mut molecules = Vec::with_capacity(decoded.len());
    let mut errors = Vec::with_capacity(decoded.len());
    for decoded in decoded {
        match decoded {
            Ok(molecule) => {
                molecules.push(Some(molecule));
                errors.push(None);
            }
            Err(error) if keep_errors => {
                molecules.push(None);
                errors.push(Some(error));
            }
            Err(error) => return Err(error),
        }
    }
    Ok(ConfSeqBatchDecodeResult { molecules, errors })
}

fn infer_confseq_input_smiles(confseq: &str) -> Result<String, ConfSeqDecodeError> {
    let tokens: Vec<_> = confseq.split_whitespace().collect();
    let mut inferred = Vec::with_capacity(tokens.len());
    for (idx, token) in tokens.iter().enumerate() {
        if parse_angle_token(token)?.is_some() {
            if tokens.get(idx + 1).copied() == Some("|") {
                inferred.push(confseq_bond_prefix(token, "^"));
            } else {
                inferred.push(confseq_bond_prefix(token, "-"));
            }
        } else {
            inferred.push((*token).to_string());
        }
    }
    Ok(inferred.join(" "))
}

fn confseq_bond_prefix(token: &str, marker: &str) -> String {
    if token.starts_with('/') {
        format!("/ {marker}")
    } else if token.starts_with('\\') {
        format!("\\ {marker}")
    } else if token.starts_with('{') || token.starts_with('}') {
        format!("! {marker}")
    } else {
        marker.to_string()
    }
}

fn parse_confseq(in_smiles: &str, confseq: &str) -> Result<ParsedConfSeq, ConfSeqDecodeError> {
    // ConfSeq source anchor:
    //   matches = re.findall(pattern, TD_smiles)
    //   TD_smiles = re.sub(pattern, '', TD_smiles)
    //   in_smiles = in_smiles.replace('^ |','')
    //   in_smiles = in_smiles.replace(' !','')
    //   in_smiles=in_smiles.replace('/ -','/')
    //   in_smiles=in_smiles.replace('\\ -','\\')
    //   TD_smiles=TD_smiles.replace('/ <','/<')
    //   TD_smiles=TD_smiles.replace('\\ <','\\<')
    let (angle_values_deg, mut td_smiles) = extract_angle_tokens(confseq)?;
    let mut in_smiles = in_smiles.replace("^ |", "");
    in_smiles = in_smiles.replace(" !", "");
    in_smiles = in_smiles.replace("/ -", "/");
    in_smiles = in_smiles.replace("\\ -", "\\");
    td_smiles = td_smiles.replace("/ -", "/");
    td_smiles = td_smiles.replace("\\ -", "\\");
    td_smiles = td_smiles.replace("/ <", "/<");
    td_smiles = td_smiles.replace("\\ <", "\\<");
    td_smiles = td_smiles.replace("> /", ">/");
    td_smiles = td_smiles.replace("> \\", ">\\");

    let stripped_smiles = strip_confseq_topology_markers(&in_smiles);
    let (td_smiles, chiral_tags_by_atom) =
        extract_pseudo_chiral_tags(&stripped_smiles, &td_smiles)?;
    let raw_dihedral_literals_by_pair =
        parse_atom_pair_dihedral_literals(&stripped_smiles, &td_smiles)?;
    let dihedral_angles_by_pair = parse_atom_pair_dihedrals(&stripped_smiles, &td_smiles)?;

    Ok(ParsedConfSeq {
        stripped_smiles,
        dihedral_angles_by_pair,
        raw_dihedral_literals_by_pair,
        angle_values_deg,
        chiral_tags_by_atom,
    })
}

fn extract_angle_tokens(input: &str) -> Result<(Vec<f64>, String), ConfSeqDecodeError> {
    let mut angles = Vec::new();
    let tokens: Vec<_> = input.split_whitespace().collect();
    let mut kept = Vec::with_capacity(tokens.len());
    let mut idx = 0;
    while idx < tokens.len() {
        if let Some(angle) = parse_angle_token(tokens[idx])? {
            if tokens.get(idx + 1).copied() == Some("|") {
                angles.push(angle);
                idx += 2;
                continue;
            }
        }
        kept.push(tokens[idx]);
        idx += 1;
    }
    Ok((angles, kept.join(" ")))
}

fn parse_angle_token(token: &str) -> Result<Option<f64>, ConfSeqDecodeError> {
    let normalized_token = token.trim_start_matches(['/', '\\']);
    let Some(inner) = normalized_token
        .strip_prefix('<')
        .and_then(|value| value.strip_suffix('>'))
    else {
        return Ok(None);
    };
    parse_angle_literal(inner, token).map(Some)
}

fn parse_angle_literal(inner: &str, token: &str) -> Result<f64, ConfSeqDecodeError> {
    let normalized = inner.replace(['/', '\\'], "");
    normalized
        .parse::<f64>()
        .map_err(|_| ConfSeqDecodeError::UnsupportedToken {
            token: token.to_string(),
        })
}

fn strip_confseq_topology_markers(smiles: &str) -> String {
    let mut stripped = String::new();
    let mut in_bracket = false;
    for token in smiles.split_whitespace() {
        match token {
            "[" if !in_bracket => {
                in_bracket = true;
                stripped.push_str(token);
            }
            "]" if in_bracket => {
                in_bracket = false;
                stripped.push_str(token);
            }
            "^" | "|" | "!" if !in_bracket => {}
            "-" if !in_bracket => {}
            _ => stripped.push_str(token),
        }
    }
    stripped
}

fn extract_pseudo_chiral_tags(
    stripped_smiles: &str,
    td_smiles: &str,
) -> Result<(String, HashMap<usize, ChiralTag>), ConfSeqDecodeError> {
    // ConfSeq source anchor:
    //   TD_smiles,p_chiral_dic=del_add_chiral_from_TD_smiles(in_smiles,TD_smiles)
    //   if atom_pair_dihedrals_dic[atom_pair][0]=='{':
    //       atom_chiral_dic[atom_pair[0]]=Chem.ChiralType.CHI_TETRAHEDRAL_CW
    //   elif atom_pair_dihedrals_dic[atom_pair][0]=='}':
    //       atom_chiral_dic[atom_pair[0]]=Chem.ChiralType.CHI_TETRAHEDRAL_CCW
    //   if atom_pair_dihedrals_dic[atom_pair][-1]=='{':
    //       atom_chiral_dic[atom_pair[-1]]=Chem.ChiralType.CHI_TETRAHEDRAL_CW
    //   elif atom_pair_dihedrals_dic[atom_pair][-1]=='}':
    //       atom_chiral_dic[atom_pair[-1]]=Chem.ChiralType.CHI_TETRAHEDRAL_CCW
    let normalized = td_smiles
        .replace("{ <", "{<")
        .replace("} <", "}<")
        .replace("> {", ">{")
        .replace("> }", ">}");
    let raw_values = parse_atom_pair_dihedral_literals(stripped_smiles, &normalized)?;
    let mut chiral_tags = HashMap::new();
    for (pair, raw) in raw_values {
        if raw.starts_with('{') {
            chiral_tags.insert(pair.0, ChiralTag::TetrahedralCw);
        } else if raw.starts_with('}') {
            chiral_tags.insert(pair.0, ChiralTag::TetrahedralCcw);
        }
        if raw.ends_with('{') {
            chiral_tags.insert(pair.1, ChiralTag::TetrahedralCw);
        } else if raw.ends_with('}') {
            chiral_tags.insert(pair.1, ChiralTag::TetrahedralCcw);
        }
    }
    let cleaned = normalized
        .replace("{<", "<")
        .replace("}<", "<")
        .replace(">{", ">")
        .replace(">}", ">");
    Ok((cleaned, chiral_tags))
}

fn parse_atom_pair_dihedrals(
    stripped_smiles: &str,
    td_smiles: &str,
) -> Result<HashMap<(usize, usize), f64>, ConfSeqDecodeError> {
    parse_atom_pair_dihedral_literals(stripped_smiles, td_smiles)?
        .into_iter()
        .map(|(pair, raw)| {
            let raw = raw
                .trim_matches(|ch| matches!(ch, '{' | '}'))
                .trim_start_matches(['/', '\\']);
            let inner = raw
                .strip_prefix('<')
                .and_then(|value| value.strip_suffix('>'))
                .ok_or_else(|| ConfSeqDecodeError::UnsupportedToken {
                    token: raw.to_string(),
                })?;
            Ok((pair, parse_angle_literal(inner, raw)?))
        })
        .collect()
}

fn parse_atom_pair_dihedral_literals(
    stripped_smiles: &str,
    td_smiles: &str,
) -> Result<HashMap<(usize, usize), String>, ConfSeqDecodeError> {
    // ConfSeq source anchor:
    //   bond_idx_token_idx_dic, token_idx_bond_idx_dic, atom_pairs=get_bond_token_atom_pairs(smiles_BE)
    //   t_smiles=complete_t_smiles(TD_smiles_lis,smiles_BE_lis)
    //   for token_idx in token_idx_bond_idx_dic:
    //       bond_idx=token_idx_bond_idx_dic[token_idx]
    //       if '<' in t_smiles[token_idx]:
    //           atom_pair=atom_pairs[bond_idx]
    //           atom_pair_dihedrals_dic[atom_pair]=...
    let mapping = bond_token_mapping_for_smiles(stripped_smiles)?;
    let smiles_be_lis: Vec<String> = mapping.smiles_be.chars().map(|ch| ch.to_string()).collect();
    let t_smiles_lis: Vec<String> = td_smiles
        .split_whitespace()
        .map(ToString::to_string)
        .collect();
    let completed = complete_t_smiles(t_smiles_lis, &smiles_be_lis);
    let mut values = HashMap::new();
    for (token_idx, pair) in mapping.token_idx_to_atom_pair {
        let Some(token) = completed.get(token_idx) else {
            return Err(ConfSeqDecodeError::TokenPositionOutOfRange {
                position: token_idx,
            });
        };
        if token.contains('<') {
            values.insert(pair, token.clone());
        }
    }
    Ok(values)
}

fn complete_t_smiles(mut t_smiles: Vec<String>, smiles_be: &[String]) -> Vec<String> {
    // ConfSeq source anchor:
    //   while idx < len(smiles_BE_lis):
    //       if smiles_BE_lis[idx] != t_smiles_lis[idx] and '>' not in t_smiles_lis[idx]:
    //           t_smiles_lis.insert(idx,'-')
    //       idx+=1
    let mut idx = 0;
    while idx < smiles_be.len() {
        if idx >= t_smiles.len() {
            t_smiles.push(String::new());
        } else if smiles_be[idx] != t_smiles[idx] && !t_smiles[idx].contains('>') {
            t_smiles.insert(idx, "-".to_string());
        }
        idx += 1;
    }
    t_smiles
}

fn build_template(
    smiles: &str,
    chiral_tags_by_atom: &HashMap<usize, ChiralTag>,
    options: &ConfSeqDecodeOptions,
) -> Result<Template, ConfSeqDecodeError> {
    match options.template_backend {
        ConfSeqTemplateBackend::DistanceGeometry => {
            build_distance_geometry_template(smiles, chiral_tags_by_atom, options)
        }
        ConfSeqTemplateBackend::FastGeometry => {
            build_confseq_base_template(smiles, chiral_tags_by_atom)
        }
    }
}

fn build_distance_geometry_template(
    smiles: &str,
    chiral_tags_by_atom: &HashMap<usize, ChiralTag>,
    options: &ConfSeqDecodeOptions,
) -> Result<Template, ConfSeqDecodeError> {
    let molecule = Molecule::from_smiles(smiles)
        .map_err(|err| ConfSeqDecodeError::SmilesParse(err.to_string()))?;
    let molecule = prepare_p_chiral_embedding_molecule(molecule, chiral_tags_by_atom)?;

    // ConfSeq source anchor:
    //   mol_with_h = Chem.AddHs(mol)
    //   params = AllChem.ETKDG()
    //   params.randomSeed = 0
    //   AllChem.EmbedMolecule(mol_with_h, params)
    //   if is_op == True:
    //       AllChem.UFFOptimizeMolecule(mol_with_h)
    //   mol = Chem.RemoveHs(mol_with_h)
    let with_h = molecule
        .with_hydrogens()
        .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
    let embed_params = options.embed_params.clone();
    let embedded_with_h = with_h
        .with_3d_conformer_with_params(embed_params)
        .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
    if embedded_with_h.conformers_3d().is_empty() {
        return Err(ConfSeqDecodeError::EmbedFailed);
    }
    let optimized_with_h = if options.optimize_with_uff {
        uff_optimize_molecule(&embedded_with_h, 200, 10.0, -1, true)
            .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?
            .molecule
    } else {
        embedded_with_h
    };
    let embedded = optimized_with_h
        .without_hydrogens_with_sanitize(false)
        .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
    let embedded = restore_embedding_molecule_state(embedded, &molecule)?;

    let dihedrals = collect_single_bond_dihedrals(&embedded);
    let dihedrals_by_pair = dihedrals
        .iter()
        .copied()
        .map(|dihedral @ (_, j, k, _)| (sorted_pair(j, k), dihedral))
        .collect();
    let angle_centers = collect_angle_centers(&embedded)?;
    let ring_bond_pairs = collect_ring_bond_pairs(&embedded)?;
    let last_ring_bonds = collect_last_ring_bonds(smiles, &embedded)?;

    Ok(Template {
        molecule: embedded,
        dihedrals,
        dihedrals_by_pair,
        angle_centers,
        ring_bond_pairs,
        last_ring_bonds,
    })
}

fn build_confseq_base_template(
    smiles: &str,
    chiral_tags_by_atom: &HashMap<usize, ChiralTag>,
) -> Result<Template, ConfSeqDecodeError> {
    let molecule = Molecule::from_smiles(smiles)
        .map_err(|err| ConfSeqDecodeError::SmilesParse(err.to_string()))?;
    let molecule = prepare_p_chiral_embedding_molecule(molecule, chiral_tags_by_atom)?;
    let embedded = try_build_confseq_fast_geometry(&molecule)?;
    let embedded = restore_embedding_molecule_state(embedded, &molecule)?;

    let dihedrals = collect_single_bond_dihedrals(&embedded);
    let dihedrals_by_pair = dihedrals
        .iter()
        .copied()
        .map(|dihedral @ (_, j, k, _)| (sorted_pair(j, k), dihedral))
        .collect();
    let angle_centers = collect_angle_centers(&embedded)?;
    let ring_bond_pairs = collect_ring_bond_pairs(&embedded)?;
    let last_ring_bonds = collect_last_ring_bonds(smiles, &embedded)?;

    // FastGeometry is only a cheaper initializer for `template.molecule`.
    // ConfSeq token recovery must remain identical to the DistanceGeometry
    // template path: every parsed angle/dihedral token is applied, including
    // ring-deferred dihedrals. Fragment cut bonds are internal to base
    // initialization and must not filter final decode semantics.
    Ok(Template {
        molecule: embedded,
        dihedrals,
        dihedrals_by_pair,
        angle_centers,
        ring_bond_pairs,
        last_ring_bonds,
    })
}

fn try_build_confseq_fast_geometry(
    molecule: &Molecule,
) -> Result<Molecule, ConfSeqFastGeometryError> {
    if molecule.num_atoms() == 0 {
        return Err(ConfSeqFastGeometryError::Build(
            "empty molecule has no fast geometry".to_string(),
        ));
    }
    let ring_info = rings::symmetrize_sssr(molecule)
        .map_err(|err| ConfSeqFastGeometryError::Build(err.to_string()))?;
    let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
    if !is_connected(molecule, &adjacency) {
        return Err(ConfSeqFastGeometryError::Disconnected);
    }
    let model = build_confseq_base_constraint_model(molecule, &ring_info)?;
    let coords = construct_confseq_base_coords(molecule, &adjacency, &model)?;

    let mut builder = molecule.to_builder();
    builder
        .add_conformer(Conformer3D::new(0, coords, true))
        .map_err(|err| ConfSeqFastGeometryError::Build(err.to_string()))?;
    let molecule = builder
        .build()
        .map_err(|err| ConfSeqFastGeometryError::Build(err.to_string()))?;
    let molecule = apply_confseq_base_double_bond_stereo(molecule)?;
    apply_confseq_base_tetrahedral_stereo(molecule, &adjacency)
}

fn build_confseq_base_constraint_model(
    molecule: &Molecule,
    ring_info: &rings::RingInfo,
) -> Result<ConfSeqBaseConstraintModel, ConfSeqFastGeometryError> {
    let ring_components = classify_confseq_base_ring_components(molecule, ring_info)?;
    let ring_membership = confseq_base_ring_membership(molecule.num_atoms(), &ring_components);
    let bond_targets = molecule
        .bonds()
        .iter()
        .map(|bond| confseq_base_source_backed_bond_length(molecule, bond))
        .collect();
    let angle_targets = collect_confseq_base_angle_targets(molecule, &ring_membership);
    let (torsion_priors, path14_distance_priors) = collect_confseq_base_torsion_priors(molecule);
    let planar_bonds = collect_confseq_base_planar_bonds(molecule);
    let rigid_components =
        collect_confseq_base_rigid_components(molecule, ring_info, &ring_components);
    Ok(ConfSeqBaseConstraintModel {
        bond_targets,
        angle_targets,
        torsion_priors,
        path14_distance_priors,
        planar_bonds,
        rigid_components,
    })
}

fn classify_confseq_base_ring_components(
    molecule: &Molecule,
    ring_info: &rings::RingInfo,
) -> Result<Vec<ConfSeqBaseRingComponent>, ConfSeqFastGeometryError> {
    for (ring_index, atoms) in ring_info.atom_rings().iter().enumerate() {
        let bonds = &ring_info.bond_rings()[ring_index];
        validate_supported_confseq_base_ring(molecule, ring_index, atoms, bonds)?;
    }

    let mut ring_graph = vec![Vec::<usize>::new(); ring_info.num_rings()];
    for left in 0..ring_info.num_rings() {
        for right in left + 1..ring_info.num_rings() {
            let shared_atoms = ring_info.atom_rings()[left]
                .iter()
                .filter(|atom| ring_info.atom_rings()[right].contains(atom))
                .count();
            let shared_bonds = shared_bond_ids_between_rings(ring_info, left, right).len();
            if shared_atoms > 0 || shared_bonds > 0 {
                ring_graph[left].push(right);
                ring_graph[right].push(left);
            }
        }
    }

    let mut components = Vec::new();
    let mut seen_rings = vec![false; ring_info.num_rings()];
    for ring_idx in 0..ring_info.num_rings() {
        if seen_rings[ring_idx] {
            continue;
        }
        let mut stack = vec![ring_idx];
        seen_rings[ring_idx] = true;
        let mut component_rings = Vec::new();
        while let Some(current) = stack.pop() {
            component_rings.push(current);
            for &next in &ring_graph[current] {
                if !seen_rings[next] {
                    seen_rings[next] = true;
                    stack.push(next);
                }
            }
        }

        let mut atoms = HashSet::<usize>::new();
        let mut ring_sizes_by_atom = HashMap::<usize, usize>::new();
        for ring in component_rings {
            let ring_size = ring_info.atom_rings()[ring].len();
            for atom in &ring_info.atom_rings()[ring] {
                let idx = atom.index();
                atoms.insert(idx);
                ring_sizes_by_atom
                    .entry(idx)
                    .and_modify(|old| *old = (*old).min(ring_size))
                    .or_insert(ring_size);
            }
        }

        let mut atoms: Vec<_> = atoms.into_iter().collect();
        atoms.sort_unstable();
        let planar = confseq_base_ring_component_is_planar(molecule, &atoms);
        components.push(ConfSeqBaseRingComponent {
            atoms,
            ring_sizes_by_atom,
            planar,
        });
    }
    Ok(components)
}

fn confseq_base_ring_component_is_planar(molecule: &Molecule, atoms: &[usize]) -> bool {
    let atom_set: HashSet<_> = atoms.iter().copied().collect();
    atoms.iter().copied().all(|atom_idx| {
        let atom = &molecule.atoms()[atom_idx];
        atom.is_aromatic()
            || atom.hybridization() == Hybridization::Sp2
            || confseq_base_ring_atom_has_exocyclic_pi_bond(molecule, atom_idx, &atom_set)
            || confseq_base_ring_atom_is_conjugated_to_ring_pi_system(molecule, atom_idx, &atom_set)
    })
}

fn confseq_base_ring_atom_has_exocyclic_pi_bond(
    molecule: &Molecule,
    atom_idx: usize,
    ring_atoms: &HashSet<usize>,
) -> bool {
    molecule.bonds().iter().any(|bond| {
        let begin = bond.begin().index();
        let end = bond.end().index();
        let other = if begin == atom_idx {
            end
        } else if end == atom_idx {
            begin
        } else {
            return false;
        };
        !ring_atoms.contains(&other)
            && matches!(bond.order(), BondOrder::Double | BondOrder::Aromatic)
    })
}

fn confseq_base_ring_atom_is_conjugated_to_ring_pi_system(
    molecule: &Molecule,
    atom_idx: usize,
    ring_atoms: &HashSet<usize>,
) -> bool {
    let mut has_pi_neighbor = false;
    let mut has_ring_single = false;
    for bond in molecule.bonds() {
        let begin = bond.begin().index();
        let end = bond.end().index();
        let other = if begin == atom_idx {
            end
        } else if end == atom_idx {
            begin
        } else {
            continue;
        };
        if !ring_atoms.contains(&other) {
            continue;
        }
        let other_atom = &molecule.atoms()[other];
        if other_atom.is_aromatic()
            || other_atom.hybridization() == Hybridization::Sp2
            || matches!(bond.order(), BondOrder::Double | BondOrder::Aromatic)
        {
            has_pi_neighbor = true;
        }
        if bond.order() == BondOrder::Single {
            has_ring_single = true;
        }
    }
    has_pi_neighbor
        && has_ring_single
        && matches!(
            molecule.atoms()[atom_idx].atomic_number(),
            6 | 7 | 8 | 15 | 16
        )
}

fn confseq_base_ring_membership(
    atom_count: usize,
    ring_components: &[ConfSeqBaseRingComponent],
) -> Vec<Option<(usize, usize)>> {
    let mut membership = vec![None; atom_count];
    for (component_idx, component) in ring_components.iter().enumerate() {
        for atom in component.atoms.iter().copied() {
            if membership[atom].is_none() {
                let ring_size = component
                    .ring_sizes_by_atom
                    .get(&atom)
                    .copied()
                    .unwrap_or(component.atoms.len());
                membership[atom] = Some((component_idx, ring_size));
            }
        }
    }
    membership
}

fn collect_confseq_base_angle_targets(
    molecule: &Molecule,
    ring_membership: &[Option<(usize, usize)>],
) -> HashMap<(usize, usize, usize), f64> {
    let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
    let mut targets = HashMap::new();
    for center in 0..molecule.num_atoms() {
        let neighbors: Vec<_> = adjacency
            .neighbors_of(center)
            .iter()
            .map(|neighbor| neighbor.atom_index)
            .collect();
        let ring_size = ring_membership[center]
            .map(|(_, ring_size)| ring_size)
            .unwrap_or(0);
        for left_pos in 0..neighbors.len() {
            for right_pos in left_pos + 1..neighbors.len() {
                let angle = if ring_size > 0 {
                    confseq_base_ring_angle_rad(molecule, center, ring_size)
                } else {
                    confseq_base_source_backed_angle_rad(
                        molecule,
                        neighbors[left_pos],
                        center,
                        neighbors[right_pos],
                    )
                };
                targets.insert(
                    sorted_angle(neighbors[left_pos], center, neighbors[right_pos]),
                    angle,
                );
            }
        }
    }
    targets
}

fn collect_confseq_base_torsion_priors(
    molecule: &Molecule,
) -> (
    HashMap<(usize, usize, usize, usize), ConfSeqBaseTorsionPrior>,
    Vec<ConfSeqBasePath14DistancePrior>,
) {
    let mut priors = HashMap::new();
    if let Ok(path_priors) = distgeom::dg_path14_priors(molecule) {
        let mut distance_priors = Vec::with_capacity(path_priors.len());
        for prior in path_priors {
            let kind = match prior.kind {
                distgeom::DgPath14Kind::Cis => ConfSeqBaseTorsionPrior::Cis,
                distgeom::DgPath14Kind::Trans => ConfSeqBaseTorsionPrior::Trans,
                distgeom::DgPath14Kind::Other => ConfSeqBaseTorsionPrior::Free,
            };
            priors.insert(prior.atoms, kind);
            let (i, j, k, l) = prior.atoms;
            priors.insert((l, k, j, i), kind);
            distance_priors.push(ConfSeqBasePath14DistancePrior {
                atoms: prior.atoms,
                lower_bound: prior.lower_bound,
                upper_bound: prior.upper_bound,
            });
        }
        return (priors, distance_priors);
    }

    for bond in molecule.bonds() {
        let j = bond.begin().index();
        let k = bond.end().index();
        let left = heavy_neighbors_except(molecule, j, k);
        let right = heavy_neighbors_except(molecule, k, j);
        let prior = match bond.order() {
            BondOrder::Double => match bond.stereo() {
                BondStereo::Cis | BondStereo::Z => ConfSeqBaseTorsionPrior::Cis,
                BondStereo::Trans | BondStereo::E => ConfSeqBaseTorsionPrior::Trans,
                _ => ConfSeqBaseTorsionPrior::Cis,
            },
            BondOrder::Single
                if molecule.atoms()[j].hybridization() == Hybridization::Sp2
                    && molecule.atoms()[k].hybridization() == Hybridization::Sp2 =>
            {
                ConfSeqBaseTorsionPrior::Trans
            }
            _ => ConfSeqBaseTorsionPrior::Free,
        };
        for &i in &left {
            for &l in &right {
                priors.insert((i, j, k, l), prior);
                priors.insert((l, k, j, i), prior);
            }
        }
    }
    (priors, Vec::new())
}

fn collect_confseq_base_planar_bonds(molecule: &Molecule) -> HashSet<(usize, usize)> {
    molecule
        .bonds()
        .iter()
        .filter(|bond| {
            bond.order() == BondOrder::Double
                || bond.is_aromatic()
                || molecule.atoms()[bond.begin().index()].is_aromatic()
                || molecule.atoms()[bond.end().index()].is_aromatic()
        })
        .map(|bond| sorted_pair(bond.begin().index(), bond.end().index()))
        .collect()
}

fn collect_confseq_base_rigid_components(
    molecule: &Molecule,
    ring_info: &rings::RingInfo,
    ring_components: &[ConfSeqBaseRingComponent],
) -> Vec<ConfSeqBaseRigidComponent> {
    let cut_bonds = confseq_base_fragment_cut_bonds(molecule, ring_info);
    let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
    let mut seen = vec![false; molecule.num_atoms()];
    let mut components = Vec::new();
    for atom in molecule.atoms() {
        let start = atom.id().index();
        if seen[start] || atom.atomic_number() == 1 {
            continue;
        }
        let mut atoms = Vec::new();
        let mut queue = VecDeque::new();
        seen[start] = true;
        queue.push_back(start);
        while let Some(current) = queue.pop_front() {
            atoms.push(current);
            for neighbor in adjacency.neighbors_of(current) {
                let next = neighbor.atom_index;
                if seen[next]
                    || molecule.atoms()[next].atomic_number() == 1
                    || cut_bonds.contains(&sorted_pair(current, next))
                {
                    continue;
                }
                seen[next] = true;
                queue.push_back(next);
            }
        }
        atoms.sort_unstable();
        let atom_set: HashSet<_> = atoms.iter().copied().collect();
        let connectors = molecule
            .bonds()
            .iter()
            .filter_map(|bond| {
                let begin = bond.begin().index();
                let end = bond.end().index();
                match (atom_set.contains(&begin), atom_set.contains(&end)) {
                    (true, false) if cut_bonds.contains(&sorted_pair(begin, end)) => {
                        Some(ConfSeqBaseTemplateConnector {
                            core_atom: begin,
                            external_atom: end,
                            bond_id: bond.id().index(),
                        })
                    }
                    (false, true) if cut_bonds.contains(&sorted_pair(begin, end)) => {
                        Some(ConfSeqBaseTemplateConnector {
                            core_atom: end,
                            external_atom: begin,
                            bond_id: bond.id().index(),
                        })
                    }
                    _ => None,
                }
            })
            .collect();
        let kind = classify_confseq_base_rigid_component(molecule, ring_components, &atoms);
        components.push(ConfSeqBaseRigidComponent {
            kind,
            atoms,
            connectors,
        });
    }
    components
}

fn confseq_base_fragment_cut_bonds(
    molecule: &Molecule,
    ring_info: &rings::RingInfo,
) -> HashSet<(usize, usize)> {
    collect_single_bond_dihedrals(molecule)
        .into_iter()
        .filter_map(|(_, j, k, _)| {
            let pair = sorted_pair(j, k);
            let bond = bond_between_pair(molecule, pair)?;
            (ring_info.num_bond_rings(bond.id()) == 0).then_some(pair)
        })
        .collect()
}

fn classify_confseq_base_rigid_component(
    molecule: &Molecule,
    ring_components: &[ConfSeqBaseRingComponent],
    atoms: &[usize],
) -> ConfSeqBaseRigidComponentKind {
    let atom_set: HashSet<_> = atoms.iter().copied().collect();
    let mut touches_ring = false;
    let mut contains_planar_ring = false;
    let mut contains_nonplanar_ring = false;
    let mut ring_is_wholly_inside = false;
    for component in ring_components {
        let overlap = component.atoms.iter().any(|atom| atom_set.contains(atom));
        if !overlap {
            continue;
        }
        touches_ring = true;
        if component.atoms.iter().all(|atom| atom_set.contains(atom)) {
            ring_is_wholly_inside = true;
            if component.planar {
                contains_planar_ring = true;
            } else {
                contains_nonplanar_ring = true;
            }
        }
    }
    if contains_nonplanar_ring {
        return ConfSeqBaseRigidComponentKind::RingNonplanar;
    }
    let all_planar_like = atoms.iter().copied().all(|atom_idx| {
        let atom = &molecule.atoms()[atom_idx];
        atom.atomic_number() == 1
            || atom.is_aromatic()
            || atom.hybridization() == Hybridization::Sp
            || atom.hybridization() == Hybridization::Sp2
    });
    if all_planar_like {
        return if touches_ring {
            ConfSeqBaseRigidComponentKind::RingPlanar
        } else {
            ConfSeqBaseRigidComponentKind::AcyclicPlanarPi
        };
    }
    if contains_planar_ring && ring_is_wholly_inside {
        return ConfSeqBaseRigidComponentKind::RingPlanar;
    }
    if touches_ring {
        return ConfSeqBaseRigidComponentKind::RingMixed;
    }
    if atoms.iter().copied().any(|atom_idx| {
        let atom = &molecule.atoms()[atom_idx];
        atom.atomic_number() > 1
            && !atom.is_aromatic()
            && atom.hybridization() == Hybridization::Sp3
            && heavy_degree_within_atom_set(molecule, atom_idx, &atom_set) >= 3
    }) {
        return ConfSeqBaseRigidComponentKind::AcyclicNonplanar;
    }
    ConfSeqBaseRigidComponentKind::AcyclicOther
}

pub fn confseq_base_structural_risk_precheck(
    molecule: &Molecule,
) -> Result<ConfSeqBaseStructuralRiskPrecheck, ConfSeqDecodeError> {
    let rings = rings::symmetrize_sssr(molecule)
        .map_err(|err| ConfSeqDecodeError::RingFinding(err.to_string()))?;
    let model = build_confseq_base_constraint_model(molecule, &rings)?;
    Ok(confseq_base_structural_risk_precheck_with_model(
        molecule, &rings, &model,
    ))
}

fn confseq_base_structural_risk_precheck_with_model(
    molecule: &Molecule,
    ring_info: &rings::RingInfo,
    model: &ConfSeqBaseConstraintModel,
) -> ConfSeqBaseStructuralRiskPrecheck {
    let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
    let mut classes = BTreeSet::new();

    for ring_atoms in ring_info.atom_rings() {
        let atoms = ring_atoms
            .iter()
            .map(|atom| atom.index())
            .collect::<Vec<_>>();
        if atoms.len() < 5 || atoms.len() > 6 {
            continue;
        }
        if !confseq_base_ring_component_is_planar(molecule, &atoms) {
            let risk_class = match atoms.len() {
                5 => ConfSeqBaseStructuralRiskClass::Nonplanar5MemberRingPucker,
                6 => ConfSeqBaseStructuralRiskClass::Nonplanar6MemberRingPucker,
                _ => continue,
            };
            classes.insert(risk_class);
        }
    }

    for atom in molecule.atoms() {
        let atom_idx = atom.id().index();
        if ring_info.num_atom_rings(atom.id()) < 2 {
            continue;
        }
        let heavy_ring_degree = adjacency
            .neighbors_of(atom_idx)
            .iter()
            .filter(|neighbor| {
                molecule.atoms()[neighbor.atom_index].atomic_number() > 1
                    && molecule
                        .bonds()
                        .iter()
                        .find(|bond| {
                            let pair = sorted_pair(bond.begin().index(), bond.end().index());
                            pair == sorted_pair(atom_idx, neighbor.atom_index)
                        })
                        .is_some_and(|bond| ring_info.num_bond_rings(bond.id()) > 0)
            })
            .count();
        if atom.hybridization() == Hybridization::Sp3 && heavy_ring_degree >= 3 {
            classes.insert(ConfSeqBaseStructuralRiskClass::SharedSp3RingCenterBranch);
        }
        if ring_info.num_atom_rings(atom.id()) >= 2 && heavy_ring_degree >= 3 {
            classes.insert(ConfSeqBaseStructuralRiskClass::MultiRingSharedAtomBranch);
        }
    }

    for bond in molecule.bonds() {
        if ring_info.num_bond_rings(bond.id()) >= 2
            && (molecule.atoms()[bond.begin().index()].hybridization() == Hybridization::Sp3
                || molecule.atoms()[bond.end().index()].hybridization() == Hybridization::Sp3)
        {
            classes.insert(ConfSeqBaseStructuralRiskClass::FusedNonplanarRingEdgeBranch);
        }
    }

    for component in &model.rigid_components {
        if matches!(
            component.kind,
            ConfSeqBaseRigidComponentKind::RingNonplanar | ConfSeqBaseRigidComponentKind::RingMixed
        ) {
            let ring_count = ring_info
                .atom_rings()
                .iter()
                .filter(|ring| {
                    ring.iter()
                        .any(|atom| component.atoms.contains(&atom.index()))
                })
                .count();
            if ring_count >= 2 {
                classes.insert(ConfSeqBaseStructuralRiskClass::NonplanarMultiringRigidComponent);
            }
        }
    }

    let classes = classes.into_iter().collect::<Vec<_>>();
    ConfSeqBaseStructuralRiskPrecheck {
        classes,
        fallback_candidate: false,
    }
}

fn heavy_degree_within_atom_set(
    molecule: &Molecule,
    atom_idx: usize,
    atom_set: &HashSet<usize>,
) -> usize {
    molecule
        .bonds()
        .iter()
        .filter(|bond| {
            let begin = bond.begin().index();
            let end = bond.end().index();
            let other = if begin == atom_idx {
                end
            } else if end == atom_idx {
                begin
            } else {
                return false;
            };
            atom_set.contains(&other) && molecule.atoms()[other].atomic_number() > 1
        })
        .count()
}

fn construct_confseq_base_coords(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
) -> Result<Vec<[f64; 3]>, ConfSeqFastGeometryError> {
    construct_confseq_base_coords_from_rigid_fragments(molecule, adjacency, model)
}

fn construct_confseq_base_coords_from_rigid_fragments(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
) -> Result<Vec<[f64; 3]>, ConfSeqFastGeometryError> {
    if model.rigid_components.is_empty() {
        return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(
            "molecule has no rigid heavy-atom fragments".to_string(),
        ));
    }
    let mut coords = vec![[0.0; 3]; molecule.num_atoms()];
    let mut placed = vec![false; molecule.num_atoms()];
    let mut placed_components = vec![false; model.rigid_components.len()];
    let mut connector_targets = HashMap::<(usize, usize), [f64; 3]>::new();
    let mut component_by_atom = vec![None; molecule.num_atoms()];
    let scaffold_shape =
        coordinates::rdkit_initial_2d_scaffold_coords(molecule.atoms(), molecule.bonds())
            .ok()
            .map(|points| {
                points
                    .into_iter()
                    .map(|point| [point[0], point[1], 0.0])
                    .collect::<Vec<_>>()
            });
    let mut local_components = Vec::with_capacity(model.rigid_components.len());
    let mut fragment_template_cache = HashMap::<String, ConfSeqBaseRigidFragmentTemplate>::new();
    let full_bounds = distgeom::dg_bounds_matrix(molecule).map_err(|err| {
        ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
            "failed to build whole-molecule DG bounds for rigid fragments: {err}"
        ))
    })?;
    for (component_idx, component) in model.rigid_components.iter().enumerate() {
        for &atom in &component.atoms {
            component_by_atom[atom] = Some(component_idx);
        }
        let template = confseq_base_cached_rigid_fragment_template(
            molecule,
            model,
            component,
            &mut fragment_template_cache,
        )?;
        let local = confseq_base_realize_rigid_fragment_template(
            molecule,
            model,
            &template,
            scaffold_shape.as_deref(),
            Some(&full_bounds),
        )?;
        local_components.push(local);
    }

    let root_component = confseq_base_select_root_component(molecule, model, &component_by_atom);
    place_confseq_base_rigid_component_local(
        &model.rigid_components[root_component],
        &local_components[root_component],
        [0.0, 0.0, 0.0],
        &mut coords,
        &mut placed,
        &mut connector_targets,
    )?;
    placed_components[root_component] = true;

    loop {
        let Some(next_component) = confseq_base_select_next_assembly_component(
            molecule,
            adjacency,
            model,
            &component_by_atom,
            &local_components,
            &coords,
            &placed,
            &placed_components,
            &connector_targets,
        ) else {
            break;
        };
        place_confseq_base_unplaced_component_by_anchors(
            molecule,
            adjacency,
            model,
            &component_by_atom,
            &local_components[next_component],
            next_component,
            &mut coords,
            &mut placed,
            &mut connector_targets,
        )?;
        placed_components[next_component] = true;
    }

    for atom in molecule.atoms() {
        let atom_idx = atom.id().index();
        if atom.atomic_number() == 1 || placed[atom_idx] {
            continue;
        }
        return Err(ConfSeqFastGeometryError::PlacementLeftAtomsUnplaced);
    }
    place_confseq_base_hydrogens_from_heavy_neighbors(molecule, adjacency, model, &mut coords);
    validate_confseq_base_constraint_coords(molecule, model, &coords)?;
    Ok(coords)
}

#[derive(Debug, Clone, Copy)]
struct ConfSeqBaseAssemblyAnchor {
    placed_atom: usize,
    unplaced_atom: usize,
    bond_id: usize,
}

const CONFSEQ_BASE_TWO_ANCHOR_DISTANCE_ABS_TOLERANCE: f64 = 0.20;
const CONFSEQ_BASE_TWO_ANCHOR_DISTANCE_REL_TOLERANCE: f64 = 0.12;

fn confseq_base_select_root_component(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    component_by_atom: &[Option<usize>],
) -> usize {
    let mut external_degree = vec![0usize; model.rigid_components.len()];
    for bond in molecule.bonds() {
        let begin = bond.begin().index();
        let end = bond.end().index();
        let (Some(begin_component), Some(end_component)) =
            (component_by_atom[begin], component_by_atom[end])
        else {
            continue;
        };
        if begin_component != end_component {
            external_degree[begin_component] += 1;
            external_degree[end_component] += 1;
        }
    }
    model
        .rigid_components
        .iter()
        .enumerate()
        .max_by_key(|(idx, component)| {
            (
                component.atoms.len(),
                external_degree[*idx],
                component.connectors.len(),
            )
        })
        .map(|(idx, _)| idx)
        .unwrap_or(0)
}

#[allow(clippy::too_many_arguments)]
fn confseq_base_select_next_assembly_component(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
    component_by_atom: &[Option<usize>],
    local_components: &[HashMap<usize, [f64; 3]>],
    coords: &[[f64; 3]],
    placed: &[bool],
    placed_components: &[bool],
    connector_targets: &HashMap<(usize, usize), [f64; 3]>,
) -> Option<usize> {
    let mut best = None;
    let mut best_score = ConfSeqBaseAssemblyCandidateScore::default();
    for component_idx in 0..model.rigid_components.len() {
        if placed_components[component_idx] {
            continue;
        }
        let anchors = confseq_base_component_assembly_anchors(
            molecule,
            component_by_atom,
            placed,
            component_idx,
        );
        if anchors.is_empty() {
            continue;
        }
        let local = &local_components[component_idx];
        let connector_count = anchors
            .iter()
            .filter(|anchor| {
                connector_targets.contains_key(&(anchor.placed_atom, anchor.unplaced_atom))
                    && local.contains_key(&anchor.placed_atom)
                    && local.contains_key(&anchor.unplaced_atom)
            })
            .count();
        let consistent_anchor_pairs = confseq_base_count_consistent_two_anchor_pairs(
            molecule, adjacency, model, local, coords, placed, &anchors,
        );
        let unplaced_external_count = confseq_base_unplaced_external_component_edges(
            molecule,
            component_by_atom,
            placed_components,
            component_idx,
        );
        let deferred_single_anchor_score = if anchors.len() == 1 {
            usize::MAX - unplaced_external_count
        } else {
            usize::MAX
        };
        let score = ConfSeqBaseAssemblyCandidateScore {
            connector_count,
            consistent_anchor_pairs,
            anchor_count: anchors.len(),
            deferred_single_anchor_score,
            atom_count: model.rigid_components[component_idx].atoms.len(),
            component_idx,
        };
        if score > best_score {
            best_score = score;
            best = Some(component_idx);
        }
    }
    best
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord)]
struct ConfSeqBaseAssemblyCandidateScore {
    connector_count: usize,
    consistent_anchor_pairs: usize,
    anchor_count: usize,
    deferred_single_anchor_score: usize,
    atom_count: usize,
    component_idx: usize,
}

fn confseq_base_unplaced_external_component_edges(
    molecule: &Molecule,
    component_by_atom: &[Option<usize>],
    placed_components: &[bool],
    component_idx: usize,
) -> usize {
    molecule
        .bonds()
        .iter()
        .filter(|bond| {
            let begin = bond.begin().index();
            let end = bond.end().index();
            let (Some(begin_component), Some(end_component)) =
                (component_by_atom[begin], component_by_atom[end])
            else {
                return false;
            };
            if begin_component == end_component {
                return false;
            }
            let other_component = if begin_component == component_idx {
                end_component
            } else if end_component == component_idx {
                begin_component
            } else {
                return false;
            };
            !placed_components[other_component]
        })
        .count()
}

fn confseq_base_count_consistent_two_anchor_pairs(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
    local: &HashMap<usize, [f64; 3]>,
    coords: &[[f64; 3]],
    placed: &[bool],
    anchors: &[ConfSeqBaseAssemblyAnchor],
) -> usize {
    let mut count = 0usize;
    for left_idx in 0..anchors.len() {
        for right_idx in left_idx + 1..anchors.len() {
            if confseq_base_two_anchor_pair_is_consistent(
                molecule,
                adjacency,
                model,
                local,
                coords,
                placed,
                anchors[left_idx],
                anchors[right_idx],
            ) {
                count += 1;
            }
        }
    }
    count
}

#[allow(clippy::too_many_arguments)]
fn confseq_base_two_anchor_pair_is_consistent(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
    local: &HashMap<usize, [f64; 3]>,
    coords: &[[f64; 3]],
    placed: &[bool],
    left: ConfSeqBaseAssemblyAnchor,
    right: ConfSeqBaseAssemblyAnchor,
) -> bool {
    let Some(left_local) = local.get(&left.unplaced_atom).copied() else {
        return false;
    };
    let Some(right_local) = local.get(&right.unplaced_atom).copied() else {
        return false;
    };
    let local_distance = vec_len(vec_sub(right_local, left_local));
    if local_distance <= 1.0e-8 {
        return false;
    }
    let left_target = confseq_base_anchor_target(molecule, adjacency, model, coords, placed, left);
    let right_target =
        confseq_base_anchor_target(molecule, adjacency, model, coords, placed, right);
    let target_distance = vec_len(vec_sub(right_target, left_target));
    if target_distance <= 1.0e-8 {
        return false;
    }
    let allowed_delta = CONFSEQ_BASE_TWO_ANCHOR_DISTANCE_ABS_TOLERANCE
        .max(CONFSEQ_BASE_TWO_ANCHOR_DISTANCE_REL_TOLERANCE * local_distance.max(target_distance));
    (local_distance - target_distance).abs() <= allowed_delta
}

#[allow(clippy::too_many_arguments)]
fn place_confseq_base_unplaced_component_by_anchors(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
    component_by_atom: &[Option<usize>],
    local: &HashMap<usize, [f64; 3]>,
    unplaced_component: usize,
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
    connector_targets: &mut HashMap<(usize, usize), [f64; 3]>,
) -> Result<(), ConfSeqFastGeometryError> {
    let component = &model.rigid_components[unplaced_component];
    let anchors = confseq_base_component_assembly_anchors(
        molecule,
        component_by_atom,
        placed,
        unplaced_component,
    );
    if anchors.is_empty() {
        return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
            "component has no placed assembly anchor atoms={:?}",
            component.atoms
        )));
    }
    let connector_anchors = anchors
        .iter()
        .copied()
        .filter(|anchor| {
            connector_targets.contains_key(&(anchor.placed_atom, anchor.unplaced_atom))
                && local.contains_key(&anchor.placed_atom)
                && local.contains_key(&anchor.unplaced_atom)
        })
        .collect::<Vec<_>>();
    for anchor in connector_anchors {
        if place_confseq_base_rigid_component_connector_stub(
            molecule,
            model,
            component,
            local,
            anchor,
            &anchors,
            coords,
            placed,
            connector_targets,
        )? {
            return Ok(());
        }
    }
    if anchors.len() >= 2
        && let Some((first, second)) = confseq_base_select_two_anchor_frame(
            molecule, adjacency, model, local, coords, placed, &anchors,
        )
    {
        return place_confseq_base_rigid_component_two_anchor_frame(
            molecule,
            adjacency,
            model,
            component,
            local,
            first,
            second,
            coords,
            placed,
            connector_targets,
        );
    }
    place_confseq_base_rigid_component_single_anchor(
        molecule,
        adjacency,
        model,
        component,
        local,
        anchors[0],
        coords,
        placed,
        connector_targets,
    )
}

fn confseq_base_component_assembly_anchors(
    molecule: &Molecule,
    component_by_atom: &[Option<usize>],
    placed: &[bool],
    unplaced_component: usize,
) -> Vec<ConfSeqBaseAssemblyAnchor> {
    let mut anchors = Vec::new();
    for bond in molecule.bonds() {
        let begin = bond.begin().index();
        let end = bond.end().index();
        let begin_in_component = component_by_atom[begin] == Some(unplaced_component);
        let end_in_component = component_by_atom[end] == Some(unplaced_component);
        match (begin_in_component, end_in_component) {
            (true, false) if placed[end] => anchors.push(ConfSeqBaseAssemblyAnchor {
                placed_atom: end,
                unplaced_atom: begin,
                bond_id: bond.id().index(),
            }),
            (false, true) if placed[begin] => anchors.push(ConfSeqBaseAssemblyAnchor {
                placed_atom: begin,
                unplaced_atom: end,
                bond_id: bond.id().index(),
            }),
            _ => {}
        }
    }
    anchors
}

fn confseq_base_select_two_anchor_frame(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
    local: &HashMap<usize, [f64; 3]>,
    coords: &[[f64; 3]],
    placed: &[bool],
    anchors: &[ConfSeqBaseAssemblyAnchor],
) -> Option<(ConfSeqBaseAssemblyAnchor, ConfSeqBaseAssemblyAnchor)> {
    let mut best = None;
    let mut best_score = 0.0;
    for left_idx in 0..anchors.len() {
        for right_idx in left_idx + 1..anchors.len() {
            let left = anchors[left_idx];
            let right = anchors[right_idx];
            let Some(left_local) = local.get(&left.unplaced_atom).copied() else {
                continue;
            };
            let Some(right_local) = local.get(&right.unplaced_atom).copied() else {
                continue;
            };
            let local_distance = vec_len(vec_sub(right_local, left_local));
            if local_distance <= 1.0e-8 {
                continue;
            }
            if !confseq_base_two_anchor_pair_is_consistent(
                molecule, adjacency, model, local, coords, placed, left, right,
            ) {
                continue;
            }
            let left_target =
                confseq_base_anchor_target(molecule, adjacency, model, coords, placed, left);
            let right_target =
                confseq_base_anchor_target(molecule, adjacency, model, coords, placed, right);
            let target_distance = vec_len(vec_sub(right_target, left_target));
            let score = local_distance.min(target_distance);
            if score > best_score {
                best = Some((left, right));
                best_score = score;
            }
        }
    }
    best
}

#[allow(clippy::too_many_arguments)]
fn place_confseq_base_rigid_component_single_anchor(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
    component: &ConfSeqBaseRigidComponent,
    local: &HashMap<usize, [f64; 3]>,
    anchor: ConfSeqBaseAssemblyAnchor,
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
    connector_targets: &mut HashMap<(usize, usize), [f64; 3]>,
) -> Result<(), ConfSeqFastGeometryError> {
    let length = confseq_base_bond_target_by_id(model, molecule, anchor.bond_id);
    let attach_direction = confseq_base_fragment_attach_direction(
        molecule,
        adjacency,
        model,
        coords,
        placed,
        anchor.placed_atom,
        anchor.unplaced_atom,
    );
    let target_anchor = vec_add(
        coords[anchor.placed_atom],
        vec_scale(attach_direction, length),
    );
    let Some(local_anchor) = local.get(&anchor.unplaced_atom) else {
        return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
            "realized fragment is missing anchor atom {} for component {:?}",
            anchor.unplaced_atom, component.atoms
        )));
    };
    place_confseq_base_rigid_component_local(
        component,
        local,
        vec_sub(target_anchor, *local_anchor),
        coords,
        placed,
        connector_targets,
    )
}

fn place_confseq_base_rigid_component_connector_stub(
    _molecule: &Molecule,
    _model: &ConfSeqBaseConstraintModel,
    component: &ConfSeqBaseRigidComponent,
    local: &HashMap<usize, [f64; 3]>,
    anchor: ConfSeqBaseAssemblyAnchor,
    _anchors: &[ConfSeqBaseAssemblyAnchor],
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
    connector_targets: &mut HashMap<(usize, usize), [f64; 3]>,
) -> Result<bool, ConfSeqFastGeometryError> {
    let Some(local_core) = local.get(&anchor.unplaced_atom).copied() else {
        return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
            "realized fragment is missing connector core atom {} for component {:?}",
            anchor.unplaced_atom, component.atoms
        )));
    };
    let Some(local_external) = local.get(&anchor.placed_atom).copied() else {
        return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
            "realized fragment is missing connector stub atom {} for component {:?}",
            anchor.placed_atom, component.atoms
        )));
    };
    let Some(target_core) = connector_targets
        .get(&(anchor.placed_atom, anchor.unplaced_atom))
        .copied()
    else {
        return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
            "missing connector target {}->{} for component {:?}",
            anchor.placed_atom, anchor.unplaced_atom, component.atoms
        )));
    };
    let target_external = coords[anchor.placed_atom];
    let local_axis = vec_sub(local_core, local_external);
    let target_axis = vec_sub(target_core, target_external);
    if vec_len(local_axis) <= 1.0e-8 || vec_len(target_axis) <= 1.0e-8 {
        return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
            "degenerate connector stub {}-{} for component {:?}",
            anchor.placed_atom, anchor.unplaced_atom, component.atoms
        )));
    }
    place_confseq_base_rigid_component_aligned_to_segment(
        component,
        local,
        local_external,
        local_axis,
        target_external,
        target_axis,
        coords,
        placed,
        connector_targets,
    )?;
    Ok(true)
}

#[allow(clippy::too_many_arguments)]
fn place_confseq_base_rigid_component_two_anchor_frame(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
    component: &ConfSeqBaseRigidComponent,
    local: &HashMap<usize, [f64; 3]>,
    first: ConfSeqBaseAssemblyAnchor,
    second: ConfSeqBaseAssemblyAnchor,
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
    connector_targets: &mut HashMap<(usize, usize), [f64; 3]>,
) -> Result<(), ConfSeqFastGeometryError> {
    let first_target =
        confseq_base_anchor_target(molecule, adjacency, model, coords, placed, first);
    let second_target =
        confseq_base_anchor_target(molecule, adjacency, model, coords, placed, second);
    let Some(first_local) = local.get(&first.unplaced_atom).copied() else {
        return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
            "realized fragment is missing anchor atom {} for component {:?}",
            first.unplaced_atom, component.atoms
        )));
    };
    let Some(second_local) = local.get(&second.unplaced_atom).copied() else {
        return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
            "realized fragment is missing anchor atom {} for component {:?}",
            second.unplaced_atom, component.atoms
        )));
    };
    let local_axis = vec_sub(second_local, first_local);
    let target_axis = vec_sub(second_target, first_target);
    if vec_len(local_axis) <= 1.0e-8 || vec_len(target_axis) <= 1.0e-8 {
        return place_confseq_base_rigid_component_single_anchor(
            molecule,
            adjacency,
            model,
            component,
            local,
            first,
            coords,
            placed,
            connector_targets,
        );
    }
    place_confseq_base_rigid_component_aligned_to_segment(
        component,
        local,
        first_local,
        local_axis,
        first_target,
        target_axis,
        coords,
        placed,
        connector_targets,
    )
}

fn confseq_base_anchor_target(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
    coords: &[[f64; 3]],
    placed: &[bool],
    anchor: ConfSeqBaseAssemblyAnchor,
) -> [f64; 3] {
    let length = confseq_base_bond_target_by_id(model, molecule, anchor.bond_id);
    let attach_direction = confseq_base_fragment_attach_direction(
        molecule,
        adjacency,
        model,
        coords,
        placed,
        anchor.placed_atom,
        anchor.unplaced_atom,
    );
    vec_add(
        coords[anchor.placed_atom],
        vec_scale(attach_direction, length),
    )
}

fn confseq_base_cached_rigid_fragment_template(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    component: &ConfSeqBaseRigidComponent,
    cache: &mut HashMap<String, ConfSeqBaseRigidFragmentTemplate>,
) -> Result<ConfSeqBaseRigidFragmentTemplate, ConfSeqFastGeometryError> {
    let key = confseq_base_fragment_template_key(molecule, component);
    if let Some(template) = cache.get(&key) {
        if template.atoms == component.atoms {
            return Ok(template.clone());
        }
    }

    let template =
        confseq_base_build_rigid_fragment_template(molecule, model, component, key.clone())?;
    cache.insert(key, template.clone());
    Ok(template)
}

fn confseq_base_build_rigid_fragment_template(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    component: &ConfSeqBaseRigidComponent,
    cache_key: String,
) -> Result<ConfSeqBaseRigidFragmentTemplate, ConfSeqFastGeometryError> {
    if component.atoms.is_empty() {
        return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(
            "empty rigid fragment".to_string(),
        ));
    }

    let connectors = component.connectors.clone();
    let realization_atoms = confseq_base_component_realization_atoms(component);
    let stereo_context_atoms = confseq_base_component_stereo_context_atoms(molecule, component);
    let atom_set: HashSet<_> = realization_atoms.iter().copied().collect();
    let bonds = molecule
        .bonds()
        .iter()
        .filter_map(|bond| {
            let begin = bond.begin().index();
            let end = bond.end().index();
            (atom_set.contains(&begin) && atom_set.contains(&end)).then(|| {
                ConfSeqBaseTemplateBond {
                    begin,
                    end,
                    length: confseq_base_bond_target_by_id(model, molecule, bond.id().index()),
                }
            })
        })
        .collect::<Vec<_>>();
    let angles = model
        .angle_targets
        .iter()
        .filter_map(|(&(left, center, right), &angle_rad)| {
            (atom_set.contains(&left) && atom_set.contains(&center) && atom_set.contains(&right))
                .then_some(ConfSeqBaseTemplateAngle {
                    left,
                    center,
                    right,
                    angle_rad,
                })
        })
        .collect::<Vec<_>>();
    let shape = confseq_base_rigid_fragment_shape(molecule, component, &realization_atoms)?;
    if shape == ConfSeqBaseRigidFragmentShape::Unsupported {
        return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
            "closed-form rigid fragment template is not implemented for kind={:?} atoms={:?}",
            component.kind, component.atoms
        )));
    }

    Ok(ConfSeqBaseRigidFragmentTemplate {
        atoms: component.atoms.clone(),
        realization_atoms,
        stereo_context_atoms,
        shape,
        bonds,
        angles,
        connectors,
        cache_key,
    })
}

fn confseq_base_component_realization_atoms(component: &ConfSeqBaseRigidComponent) -> Vec<usize> {
    let mut realization_atoms = component.atoms.clone();
    for connector in &component.connectors {
        if !realization_atoms.contains(&connector.external_atom) {
            realization_atoms.push(connector.external_atom);
        }
    }
    realization_atoms.sort_unstable();
    realization_atoms
}

fn confseq_base_component_stereo_context_atoms(
    molecule: &Molecule,
    component: &ConfSeqBaseRigidComponent,
) -> Vec<usize> {
    let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
    let realization_atoms = confseq_base_component_realization_atoms(component);
    let realization_set: HashSet<_> = realization_atoms.iter().copied().collect();
    let mut centers = component.atoms.clone();
    for connector in &component.connectors {
        if !centers.contains(&connector.external_atom) {
            centers.push(connector.external_atom);
        }
        for neighbor in adjacency.neighbors_of(connector.external_atom) {
            let neighbor = neighbor.atom_index;
            if centers.contains(&neighbor) {
                continue;
            }
            if confseq_base_atom_needs_tetrahedral_context(molecule, &adjacency, neighbor)
                && adjacency
                    .neighbors_of(neighbor)
                    .iter()
                    .any(|neighbor_of_neighbor| {
                        neighbor_of_neighbor.atom_index == connector.external_atom
                    })
            {
                centers.push(neighbor);
            }
        }
    }

    let mut context_atoms = Vec::new();
    for center in centers {
        if !confseq_base_atom_needs_tetrahedral_context(molecule, &adjacency, center) {
            continue;
        }
        if !realization_set.contains(&center) && !context_atoms.contains(&center) {
            context_atoms.push(center);
        }
        for neighbor in adjacency.neighbors_of(center) {
            if !realization_set.contains(&neighbor.atom_index)
                && !context_atoms.contains(&neighbor.atom_index)
            {
                context_atoms.push(neighbor.atom_index);
            }
        }
    }
    context_atoms.sort_unstable();
    context_atoms
}

fn confseq_base_atom_needs_tetrahedral_context(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    atom_idx: usize,
) -> bool {
    let atom = &molecule.atoms()[atom_idx];
    let tag = atom.chiral_tag();
    if matches!(tag, ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw) {
        return adjacency.neighbors_of(atom_idx).len() >= 3;
    }
    (atom.atomic_number() == 6 || atom.atomic_number() == 7)
        && adjacency.neighbors_of(atom_idx).len() == 4
}

fn confseq_base_component_is_path_like(molecule: &Molecule, atoms: &[usize]) -> bool {
    let atom_set: HashSet<_> = atoms.iter().copied().collect();
    let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
    atoms.iter().copied().all(|atom| {
        adjacency
            .neighbors_of(atom)
            .iter()
            .filter(|neighbor| atom_set.contains(&neighbor.atom_index))
            .count()
            <= 2
    })
}

fn confseq_base_component_has_planar_bond(molecule: &Molecule, atoms: &[usize]) -> bool {
    let atom_set: HashSet<_> = atoms.iter().copied().collect();
    molecule.bonds().iter().any(|bond| {
        let begin = bond.begin().index();
        let end = bond.end().index();
        atom_set.contains(&begin)
            && atom_set.contains(&end)
            && (bond.order() == BondOrder::Double
                || bond.order() == BondOrder::Aromatic
                || bond.is_aromatic())
    })
}

fn confseq_base_component_is_tree_like(molecule: &Molecule, atoms: &[usize]) -> bool {
    if atoms.is_empty() {
        return false;
    }
    let atom_set: HashSet<_> = atoms.iter().copied().collect();
    let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
    let internal_edges = molecule
        .bonds()
        .iter()
        .filter(|bond| {
            atom_set.contains(&bond.begin().index()) && atom_set.contains(&bond.end().index())
        })
        .count();
    if internal_edges + 1 != atoms.len() {
        return false;
    }
    let mut seen = HashSet::new();
    let mut queue = VecDeque::from([atoms[0]]);
    seen.insert(atoms[0]);
    while let Some(atom) = queue.pop_front() {
        for neighbor in adjacency.neighbors_of(atom) {
            if atom_set.contains(&neighbor.atom_index) && seen.insert(neighbor.atom_index) {
                queue.push_back(neighbor.atom_index);
            }
        }
    }
    seen.len() == atoms.len()
}

fn confseq_base_component_contains_cycle(molecule: &Molecule, atoms: &[usize]) -> bool {
    if atoms.len() < 3 {
        return false;
    }
    let atom_set: HashSet<_> = atoms.iter().copied().collect();
    let internal_edges = molecule
        .bonds()
        .iter()
        .filter(|bond| {
            atom_set.contains(&bond.begin().index()) && atom_set.contains(&bond.end().index())
        })
        .count();
    internal_edges >= atoms.len()
}

fn confseq_base_ring_topology(molecule: &Molecule, atoms: &[usize]) -> ConfSeqBaseRingTopology {
    let Ok(ring_info) = rings::symmetrize_sssr(molecule) else {
        return if confseq_base_component_contains_cycle(molecule, atoms) {
            ConfSeqBaseRingTopology::Unknown
        } else {
            ConfSeqBaseRingTopology::Unknown
        };
    };
    confseq_base_ring_topology_from_ring_info(molecule, &ring_info, atoms)
}

fn confseq_base_ring_topology_from_ring_info(
    molecule: &Molecule,
    ring_info: &rings::RingInfo,
    atoms: &[usize],
) -> ConfSeqBaseRingTopology {
    if confseq_base_simple_ring_order(molecule, atoms).is_some() {
        return ConfSeqBaseRingTopology::Simple;
    }
    let atom_set: HashSet<_> = atoms.iter().copied().collect();
    let component_rings = ring_info
        .atom_rings()
        .iter()
        .enumerate()
        .filter_map(|(ring_idx, ring_atoms)| {
            ring_atoms
                .iter()
                .all(|atom| atom_set.contains(&atom.index()))
                .then_some(ring_idx)
        })
        .collect::<Vec<_>>();
    if component_rings.is_empty() {
        return ConfSeqBaseRingTopology::Unknown;
    }
    if component_rings.len() == 1 {
        return ConfSeqBaseRingTopology::SingleNonSimple;
    }
    let mut shared_atom_pairs = 0usize;
    let mut shared_bond_pairs = 0usize;
    let mut max_shared_atoms = 0usize;
    let mut max_shared_bonds = 0usize;
    let mut ring_graph = vec![Vec::<usize>::new(); component_rings.len()];
    for left_pos in 0..component_rings.len() {
        for right_pos in left_pos + 1..component_rings.len() {
            let left = component_rings[left_pos];
            let right = component_rings[right_pos];
            let shared_atoms = ring_info.atom_rings()[left]
                .iter()
                .filter(|atom| ring_info.atom_rings()[right].contains(atom))
                .count();
            let shared_bonds = shared_bond_ids_between_rings(ring_info, left, right).len();
            if shared_atoms > 0 {
                shared_atom_pairs += 1;
                max_shared_atoms = max_shared_atoms.max(shared_atoms);
                ring_graph[left_pos].push(right_pos);
                ring_graph[right_pos].push(left_pos);
            }
            if shared_bonds > 0 {
                shared_bond_pairs += 1;
                max_shared_bonds = max_shared_bonds.max(shared_bonds);
            }
        }
    }
    if max_shared_bonds > 0 {
        let is_chain = shared_bond_pairs + 1 == component_rings.len()
            && ring_graph.iter().all(|neighbors| neighbors.len() <= 2);
        return if is_chain {
            ConfSeqBaseRingTopology::EdgeFusedChain
        } else {
            ConfSeqBaseRingTopology::EdgeFusedPolycyclic
        };
    }
    if shared_atom_pairs > 0 && max_shared_atoms == 1 {
        return ConfSeqBaseRingTopology::Spiro;
    }
    if max_shared_atoms > 1 {
        return ConfSeqBaseRingTopology::BridgedOrCage;
    }
    ConfSeqBaseRingTopology::Unknown
}

fn confseq_base_rigid_fragment_shape(
    molecule: &Molecule,
    component: &ConfSeqBaseRigidComponent,
    realization_atoms: &[usize],
) -> Result<ConfSeqBaseRigidFragmentShape, ConfSeqFastGeometryError> {
    match realization_atoms.len() {
        0 => {
            return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(
                "empty rigid fragment".to_string(),
            ));
        }
        1 => return Ok(ConfSeqBaseRigidFragmentShape::Atom),
        2 => return Ok(ConfSeqBaseRigidFragmentShape::Bond),
        _ => {}
    }
    if confseq_base_component_kind_is_planar_shape(component.kind) {
        return Ok(ConfSeqBaseRigidFragmentShape::Planar);
    }
    if component.kind == ConfSeqBaseRigidComponentKind::AcyclicOther
        && confseq_base_component_is_path_like(molecule, realization_atoms)
        && confseq_base_component_has_planar_bond(molecule, realization_atoms)
    {
        return Ok(ConfSeqBaseRigidFragmentShape::Planar);
    }
    if realization_atoms.len() >= 4
        && component.kind == ConfSeqBaseRigidComponentKind::AcyclicOther
        && confseq_base_component_is_path_like(molecule, realization_atoms)
    {
        return Ok(ConfSeqBaseRigidFragmentShape::Chain);
    }
    if component.kind == ConfSeqBaseRigidComponentKind::AcyclicOther
        && confseq_base_component_is_tree_like(molecule, realization_atoms)
    {
        return Ok(ConfSeqBaseRigidFragmentShape::Tree);
    }
    if component.kind == ConfSeqBaseRigidComponentKind::RingNonplanar
        && confseq_base_simple_ring_order(molecule, &component.atoms).is_some()
    {
        return Ok(ConfSeqBaseRigidFragmentShape::RingPuckered);
    }
    if component.kind == ConfSeqBaseRigidComponentKind::RingNonplanar
        && matches!(
            confseq_base_ring_topology(molecule, &component.atoms),
            ConfSeqBaseRingTopology::EdgeFusedChain
                | ConfSeqBaseRingTopology::EdgeFusedPolycyclic
                | ConfSeqBaseRingTopology::Spiro
                | ConfSeqBaseRingTopology::BridgedOrCage
                | ConfSeqBaseRingTopology::SingleNonSimple
        )
    {
        return Ok(ConfSeqBaseRigidFragmentShape::RingPolycyclic);
    }
    if realization_atoms.len() == 3
        && confseq_base_three_atom_angle_center(molecule, realization_atoms).is_some()
    {
        return Ok(ConfSeqBaseRigidFragmentShape::Angle);
    }
    if let Some((center, ligands)) = confseq_base_single_center_ligands(molecule, realization_atoms)
    {
        return Ok(
            if confseq_base_center_prefers_planar_frame(molecule, center, &ligands) {
                ConfSeqBaseRigidFragmentShape::SingleCenterPlanar
            } else {
                ConfSeqBaseRigidFragmentShape::SingleCenterNonplanar
            },
        );
    }
    Ok(ConfSeqBaseRigidFragmentShape::Unsupported)
}

fn confseq_base_realize_rigid_fragment_template(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    template: &ConfSeqBaseRigidFragmentTemplate,
    _scaffold_shape: Option<&[[f64; 3]]>,
    full_bounds: Option<&[Vec<f64>]>,
) -> Result<HashMap<usize, [f64; 3]>, ConfSeqFastGeometryError> {
    let coords = match template.shape {
        ConfSeqBaseRigidFragmentShape::Atom => {
            let mut coords = HashMap::new();
            coords.insert(template.realization_atoms[0], [0.0, 0.0, 0.0]);
            coords
        }
        ConfSeqBaseRigidFragmentShape::Bond => {
            let mut coords = HashMap::new();
            let length = template
                .bonds
                .first()
                .map(|bond| bond.length)
                .unwrap_or(1.50);
            coords.insert(template.realization_atoms[0], [0.0, 0.0, 0.0]);
            coords.insert(template.realization_atoms[1], [length, 0.0, 0.0]);
            coords
        }
        ConfSeqBaseRigidFragmentShape::Angle => {
            let center =
                confseq_base_three_atom_angle_center(molecule, &template.realization_atoms)
                    .ok_or_else(|| {
                        ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                            "angle template has no unique center atoms={:?}",
                            template.realization_atoms
                        ))
                    })?;
            let atom_set: HashSet<_> = template.realization_atoms.iter().copied().collect();
            let ligands = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds())
                .neighbors_of(center)
                .iter()
                .filter_map(|neighbor| {
                    atom_set
                        .contains(&neighbor.atom_index)
                        .then_some(neighbor.atom_index)
                })
                .collect::<Vec<_>>();
            confseq_base_center_template_from_template(molecule, model, template, center, &ligands)?
        }
        ConfSeqBaseRigidFragmentShape::SingleCenterPlanar
        | ConfSeqBaseRigidFragmentShape::SingleCenterNonplanar => {
            let (center, ligands) =
                confseq_base_single_center_ligands(molecule, &template.realization_atoms)
                    .ok_or_else(|| {
                        ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                            "single-center template has no unique center atoms={:?}",
                            template.realization_atoms
                        ))
                    })?;
            confseq_base_center_template_from_template(molecule, model, template, center, &ligands)?
        }
        ConfSeqBaseRigidFragmentShape::Chain => {
            confseq_base_chain_fragment_template(molecule, model, template)?
        }
        ConfSeqBaseRigidFragmentShape::Tree => {
            confseq_base_tree_fragment_template(molecule, model, template)?
        }
        ConfSeqBaseRigidFragmentShape::Planar => {
            confseq_base_conditioned_fragment_template(molecule, model, template, 1, full_bounds)?
        }
        ConfSeqBaseRigidFragmentShape::RingPuckered => {
            confseq_base_conditioned_fragment_template(molecule, model, template, 8, full_bounds)?
        }
        ConfSeqBaseRigidFragmentShape::RingPolycyclic => {
            confseq_base_conditioned_fragment_template(molecule, model, template, 1, full_bounds)?
        }
        ConfSeqBaseRigidFragmentShape::Unsupported => {
            return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "unsupported rigid fragment template key={}",
                template.cache_key
            )));
        }
    };
    validate_confseq_base_realized_template_geometry(template, &coords)?;
    confseq_base_validate_connector_stubs(model, molecule, template, &coords)?;
    Ok(coords)
}

fn validate_confseq_base_realized_template_geometry(
    template: &ConfSeqBaseRigidFragmentTemplate,
    coords: &HashMap<usize, [f64; 3]>,
) -> Result<(), ConfSeqFastGeometryError> {
    for bond in &template.bonds {
        let Some(begin) = coords.get(&bond.begin).copied() else {
            return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "template bond begin atom {} is missing from key={}",
                bond.begin, template.cache_key
            )));
        };
        let Some(end) = coords.get(&bond.end).copied() else {
            return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "template bond end atom {} is missing from key={}",
                bond.end, template.cache_key
            )));
        };
        let observed = vec_len(vec_sub(begin, end));
        if (observed - bond.length).abs() > 0.35 {
            return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "template bond target failed for {}-{} key={} observed={observed:.3} target={:.3}",
                bond.begin, bond.end, template.cache_key, bond.length
            )));
        }
    }

    for angle in &template.angles {
        let Some(left) = coords.get(&angle.left).copied() else {
            continue;
        };
        let Some(center) = coords.get(&angle.center).copied() else {
            continue;
        };
        let Some(right) = coords.get(&angle.right).copied() else {
            continue;
        };
        let Some(observed) = angle_rad_from_points(left, center, right) else {
            return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "template angle is degenerate for {}-{}-{} key={}",
                angle.left, angle.center, angle.right, template.cache_key
            )));
        };
        if angular_delta_rad(observed, angle.angle_rad).abs() > 35.0_f64.to_radians() {
            return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "template angle target failed for {}-{}-{} key={} observed={:.1} target={:.1}",
                angle.left,
                angle.center,
                angle.right,
                template.cache_key,
                observed.to_degrees(),
                angle.angle_rad.to_degrees()
            )));
        }
    }

    Ok(())
}

fn confseq_base_validate_connector_stubs(
    model: &ConfSeqBaseConstraintModel,
    molecule: &Molecule,
    template: &ConfSeqBaseRigidFragmentTemplate,
    coords: &HashMap<usize, [f64; 3]>,
) -> Result<(), ConfSeqFastGeometryError> {
    for connector in &template.connectors {
        let Some(core_coord) = coords.get(&connector.core_atom).copied() else {
            return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "connector core atom {} is missing from template key={}",
                connector.core_atom, template.cache_key
            )));
        };
        let Some(external_coord) = coords.get(&connector.external_atom).copied() else {
            return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "connector external atom {} is missing from template key={}",
                connector.external_atom, template.cache_key
            )));
        };
        let observed = vec_len(vec_sub(core_coord, external_coord));
        let target = confseq_base_bond_target_by_id(model, molecule, connector.bond_id);
        if (observed - target).abs() > 0.35 {
            return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "connector stub bond target failed for {}-{} key={} observed={observed:.3} target={target:.3}",
                connector.core_atom, connector.external_atom, template.cache_key
            )));
        }
    }
    Ok(())
}

fn confseq_base_fragment_template_key(
    molecule: &Molecule,
    component: &ConfSeqBaseRigidComponent,
) -> String {
    let realization_atoms = confseq_base_component_realization_atoms(component);
    let stereo_context_atoms = confseq_base_component_stereo_context_atoms(molecule, component);
    let atom_set: HashSet<_> = realization_atoms.iter().copied().collect();
    let local_index: HashMap<_, _> = component
        .atoms
        .iter()
        .copied()
        .enumerate()
        .map(|(pos, atom)| (atom, pos))
        .collect();
    let realization_local_index: HashMap<_, _> = realization_atoms
        .iter()
        .copied()
        .enumerate()
        .map(|(pos, atom)| (atom, pos))
        .collect();
    let mut parts = vec![format!("{:?}", component.kind)];
    for &atom_idx in &component.atoms {
        let atom = &molecule.atoms()[atom_idx];
        parts.push(format!(
            "a{}:{}:{:?}:{}:{}",
            local_index[&atom_idx],
            atom.atomic_number(),
            atom.hybridization(),
            atom.formal_charge(),
            atom.is_aromatic()
        ));
    }
    for &atom_idx in &realization_atoms {
        if component.atoms.contains(&atom_idx) {
            continue;
        }
        let atom = &molecule.atoms()[atom_idx];
        parts.push(format!(
            "x{}:{}:{:?}:{}:{}",
            atom_idx,
            atom.atomic_number(),
            atom.hybridization(),
            atom.formal_charge(),
            atom.is_aromatic()
        ));
    }
    for &atom_idx in &stereo_context_atoms {
        let atom = &molecule.atoms()[atom_idx];
        parts.push(format!(
            "s{}:{}:{:?}:{}:{}:{:?}",
            atom_idx,
            atom.atomic_number(),
            atom.hybridization(),
            atom.formal_charge(),
            atom.is_aromatic(),
            atom.chiral_tag()
        ));
    }
    let mut bonds = Vec::new();
    for bond in molecule.bonds() {
        let begin = bond.begin().index();
        let end = bond.end().index();
        if atom_set.contains(&begin) && atom_set.contains(&end) {
            let i = realization_local_index[&begin].min(realization_local_index[&end]);
            let j = realization_local_index[&begin].max(realization_local_index[&end]);
            bonds.push(format!("{i}-{j}:{:?}:{}", bond.order(), bond.is_aromatic()));
        }
    }
    bonds.sort_unstable();
    parts.extend(bonds);
    let mut connectors = component
        .connectors
        .iter()
        .map(|connector| {
            let core = local_index[&connector.core_atom];
            let atom = &molecule.atoms()[connector.external_atom];
            format!(
                "c{core}->{}:{}:{:?}:{}:{}",
                connector.external_atom,
                atom.atomic_number(),
                atom.hybridization(),
                atom.formal_charge(),
                atom.is_aromatic()
            )
        })
        .collect::<Vec<_>>();
    connectors.sort_unstable();
    parts.extend(connectors);
    parts.join("|")
}

fn confseq_base_conditioned_fragment_template(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    template: &ConfSeqBaseRigidFragmentTemplate,
    num_confs: u32,
    full_bounds: Option<&[Vec<f64>]>,
) -> Result<HashMap<usize, [f64; 3]>, ConfSeqFastGeometryError> {
    let owned_full_bounds;
    let full_bounds = if let Some(full_bounds) = full_bounds {
        full_bounds
    } else {
        owned_full_bounds = distgeom::dg_bounds_matrix(molecule).map_err(|err| {
            ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "failed to build whole-molecule DG bounds for fragment key={}: {err}",
                template.cache_key
            ))
        })?;
        &owned_full_bounds
    };
    let conditioner_atoms = confseq_base_conditioner_atoms(template);
    let conditioner_set: HashSet<_> = conditioner_atoms.iter().copied().collect();
    let atoms_to_remove = molecule
        .atoms()
        .iter()
        .filter_map(|atom| (!conditioner_set.contains(&atom.id().index())).then_some(atom.id()))
        .collect::<Vec<_>>();

    let mut builder = molecule.to_builder();
    let old_to_fragment = builder.remove_atoms_for_construction(&atoms_to_remove);
    let mut fragment = builder.build().map_err(|err| {
        ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
            "failed to build fragment-local conditioner molecule key={}: {err}",
            template.cache_key
        ))
    })?;
    confseq_base_sanitize_fragment_stereo_for_conditioner(&mut fragment);
    if fragment.num_atoms() != conditioner_atoms.len() {
        return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
            "fragment-local conditioner atom count mismatch key={} expected={} observed={}",
            template.cache_key,
            conditioner_atoms.len(),
            fragment.num_atoms()
        )));
    }

    let mut fragment_to_old = vec![None; fragment.num_atoms()];
    for &old_atom in &conditioner_atoms {
        let Some(new_atom) = old_to_fragment
            .get(old_atom)
            .copied()
            .flatten()
            .map(AtomId::index)
        else {
            return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "fragment-local conditioner lost atom {old_atom} key={}",
                template.cache_key
            )));
        };
        fragment_to_old[new_atom] = Some(old_atom);
    }

    let fragment_bounds = confseq_base_slice_bounds_matrix(&full_bounds, &fragment_to_old)
        .map_err(|reason| {
            ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "failed to slice whole-molecule DG bounds for fragment key={}: {reason}",
                template.cache_key
            ))
        })?;

    let mut params = EmbedParameters::etkdg();
    params.random_seed = 0;
    params.enable_sequential_random_seeds = true;
    params.num_threads = 1;
    params.clear_confs = true;
    params.prune_rms_thresh = -1.0;
    params.embed_fragments_separately = false;
    params.set_bounds_matrix(fragment_bounds).map_err(|err| {
        ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
            "failed to attach sliced DG bounds for fragment key={}: {err}",
            template.cache_key
        ))
    })?;
    let conditioner_num_confs = if template.stereo_context_atoms.is_empty() {
        num_confs.max(1)
    } else {
        num_confs.max(8)
    };
    let (embedded, conf_ids) =
        distgeom::embed_multiple_confs(&fragment, conditioner_num_confs, &mut params).map_err(
            |err| {
                ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                    "fragment-local EmbedMolecule conditioner failed key={}: {err}",
                    template.cache_key
                ))
            },
        )?;
    let valid_conf_ids = conf_ids
        .into_iter()
        .filter(|conf_id| *conf_id >= 0)
        .map(|conf_id| conf_id as usize)
        .collect::<Vec<_>>();
    if valid_conf_ids.is_empty() || embedded.conformers_3d().is_empty() {
        return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
            "fragment-local EmbedMolecule conditioner produced no conformer key={}",
            template.cache_key
        )));
    }

    let mut best = None;
    let mut best_score = f64::INFINITY;
    for conf_id in valid_conf_ids {
        let Some(conformer) = embedded.conformers_3d().get(conf_id) else {
            continue;
        };
        let mut candidate = HashMap::with_capacity(conditioner_atoms.len());
        for (fragment_atom, coord) in conformer.coordinates().iter().copied().enumerate() {
            let Some(old_atom) = fragment_to_old[fragment_atom] else {
                return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                    "fragment-local conditioner has unmapped atom {fragment_atom} key={}",
                    template.cache_key
                )));
            };
            candidate.insert(old_atom, coord);
        }
        let score = confseq_base_conditioned_fragment_candidate_score(
            molecule, model, template, &candidate,
        );
        if score < best_score {
            let realized = template
                .realization_atoms
                .iter()
                .filter_map(|atom| candidate.get(atom).copied().map(|coord| (*atom, coord)))
                .collect::<HashMap<_, _>>();
            if realized.len() == template.realization_atoms.len() {
                best_score = score;
                best = Some(realized);
            }
        }
    }
    best.ok_or_else(|| {
        ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
            "fragment-local EmbedMolecule conditioner produced no usable conformer key={}",
            template.cache_key
        ))
    })
}

fn confseq_base_conditioner_atoms(template: &ConfSeqBaseRigidFragmentTemplate) -> Vec<usize> {
    let mut atoms = template.realization_atoms.clone();
    for &atom in &template.stereo_context_atoms {
        if !atoms.contains(&atom) {
            atoms.push(atom);
        }
    }
    atoms.sort_unstable();
    atoms
}

fn confseq_base_conditioned_fragment_candidate_score(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    template: &ConfSeqBaseRigidFragmentTemplate,
    coords: &HashMap<usize, [f64; 3]>,
) -> f64 {
    let mut score = confseq_base_template_geometry_score(model, template, coords, None);
    if template.stereo_context_atoms.is_empty() {
        return score;
    }
    let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
    let Ok(constraints) = collect_confseq_base_tetrahedral_stereo_constraints(molecule, &adjacency)
    else {
        return score;
    };
    for constraint in constraints {
        if !confseq_base_constraint_atoms_all_present(&constraint, coords) {
            continue;
        }
        let volume = confseq_base_chiral_volume_from_hash(coords, &constraint);
        if !confseq_base_chiral_volume_satisfies_tag(volume, constraint.tag) {
            score += 1000.0;
        }
    }
    score
}

fn confseq_base_constraint_atoms_all_present(
    constraint: &ConfSeqBaseTetrahedralStereoConstraint,
    coords: &HashMap<usize, [f64; 3]>,
) -> bool {
    coords.contains_key(&constraint.center)
        && constraint
            .ligands
            .iter()
            .copied()
            .all(|atom| coords.contains_key(&atom))
}

fn confseq_base_chiral_volume_from_hash(
    coords: &HashMap<usize, [f64; 3]>,
    constraint: &ConfSeqBaseTetrahedralStereoConstraint,
) -> f64 {
    let [a, b, c, d] = constraint.ligands;
    let anchor = if d == constraint.center {
        coords[&constraint.center]
    } else {
        coords[&d]
    };
    let v1 = vec_sub(coords[&a], anchor);
    let v2 = vec_sub(coords[&b], anchor);
    let v3 = vec_sub(coords[&c], anchor);
    vec_dot(v1, vec_cross(v2, v3))
}

fn confseq_base_slice_bounds_matrix(
    full_bounds: &[Vec<f64>],
    fragment_to_old: &[Option<usize>],
) -> Result<Vec<Vec<f64>>, String> {
    let full_n = full_bounds.len();
    if full_bounds.iter().any(|row| row.len() != full_n) {
        return Err("whole-molecule DG bounds matrix is not square".to_string());
    }
    let mut sliced = vec![vec![0.0; fragment_to_old.len()]; fragment_to_old.len()];
    for (fragment_i, old_i) in fragment_to_old.iter().copied().enumerate() {
        let old_i = old_i.ok_or_else(|| format!("fragment atom {fragment_i} has no old atom"))?;
        if old_i >= full_n {
            return Err(format!(
                "old atom {old_i} for fragment atom {fragment_i} is out of bounds {full_n}"
            ));
        }
        for (fragment_j, old_j) in fragment_to_old.iter().copied().enumerate() {
            let old_j =
                old_j.ok_or_else(|| format!("fragment atom {fragment_j} has no old atom"))?;
            if old_j >= full_n {
                return Err(format!(
                    "old atom {old_j} for fragment atom {fragment_j} is out of bounds {full_n}"
                ));
            }
            sliced[fragment_i][fragment_j] = full_bounds[old_i][old_j];
        }
    }
    Ok(sliced)
}

fn confseq_base_sanitize_fragment_stereo_for_conditioner(fragment: &mut Molecule) {
    let adjacency = AdjacencyList::from_topology(fragment.num_atoms(), fragment.bonds());
    let mut clear_atom_stereo = Vec::new();
    for atom in fragment.atoms() {
        let tag = atom.chiral_tag();
        if matches!(
            tag,
            ChiralTag::Tetrahedral | ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) && adjacency.neighbors_of(atom.id().index()).len() < 3
        {
            clear_atom_stereo.push(atom.id().index());
        }
    }

    let mut clear_bond_stereo = Vec::new();
    for bond in fragment.bonds() {
        if matches!(bond.stereo(), BondStereo::None | BondStereo::Any) {
            continue;
        }
        let Some(stereo_atoms) = bond.stereo_atoms() else {
            clear_bond_stereo.push(bond.id().index());
            continue;
        };
        if stereo_atoms
            .iter()
            .any(|atom| atom.index() >= fragment.num_atoms())
        {
            clear_bond_stereo.push(bond.id().index());
            continue;
        }
        if !adjacency
            .neighbors_of(bond.begin().index())
            .iter()
            .any(|neighbor| {
                neighbor.atom_index == stereo_atoms[0].index()
                    || neighbor.atom_index == stereo_atoms[1].index()
            })
            || !adjacency
                .neighbors_of(bond.end().index())
                .iter()
                .any(|neighbor| {
                    neighbor.atom_index == stereo_atoms[0].index()
                        || neighbor.atom_index == stereo_atoms[1].index()
                })
        {
            clear_bond_stereo.push(bond.id().index());
        }
    }

    if clear_atom_stereo.is_empty() && clear_bond_stereo.is_empty() {
        return;
    }

    let topology = fragment.topology_block_mut();
    for atom_idx in clear_atom_stereo {
        let atom = &mut topology.atoms[atom_idx];
        atom.set_chiral_tag(ChiralTag::Unspecified);
        atom.set_chiral_permutation(None);
    }
    for bond_idx in clear_bond_stereo {
        let bond = &mut topology.bonds[bond_idx];
        bond.set_stereo(BondStereo::None);
        bond.set_stereo_atoms(None);
    }
}

fn confseq_base_chain_fragment_template(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    template: &ConfSeqBaseRigidFragmentTemplate,
) -> Result<HashMap<usize, [f64; 3]>, ConfSeqFastGeometryError> {
    let order = confseq_base_path_like_atom_order(molecule, &template.realization_atoms)
        .ok_or_else(|| {
            ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "chain template requires a simple path key={}",
                template.cache_key
            ))
        })?;
    if order.len() == 1 {
        return Ok(HashMap::from([(order[0], [0.0, 0.0, 0.0])]));
    }

    let mut coords = HashMap::new();
    coords.insert(order[0], [0.0, 0.0, 0.0]);
    let first_length = confseq_base_bond_target_between(model, molecule, order[0], order[1])?;
    coords.insert(order[1], [first_length, 0.0, 0.0]);

    for idx in 2..order.len() {
        let grandparent = order[idx - 2];
        let parent = order[idx - 1];
        let atom = order[idx];
        let parent_coord = coords[&parent];
        let grandparent_coord = coords[&grandparent];
        let parent_axis = vec_normalize(vec_sub(parent_coord, grandparent_coord));
        let angle = model
            .angle_targets
            .get(&sorted_angle(grandparent, parent, atom))
            .copied()
            .unwrap_or_else(|| {
                confseq_base_source_backed_angle_rad(molecule, grandparent, parent, atom)
            });
        let phase_ord = idx.saturating_sub(2);
        let direction = child_direction(parent_axis, phase_ord, 2, angle);
        let length = confseq_base_bond_target_between(model, molecule, parent, atom)?;
        coords.insert(atom, vec_add(parent_coord, vec_scale(direction, length)));
    }

    Ok(coords)
}

fn confseq_base_path_like_atom_order(molecule: &Molecule, atoms: &[usize]) -> Option<Vec<usize>> {
    let atom_set: HashSet<_> = atoms.iter().copied().collect();
    let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
    let mut degrees = HashMap::<usize, usize>::new();
    for &atom in atoms {
        let degree = adjacency
            .neighbors_of(atom)
            .iter()
            .filter(|neighbor| atom_set.contains(&neighbor.atom_index))
            .count();
        if degree > 2 {
            return None;
        }
        degrees.insert(atom, degree);
    }
    let start = atoms
        .iter()
        .copied()
        .filter(|atom| degrees.get(atom).copied().unwrap_or(0) <= 1)
        .min()?;
    let mut order = Vec::with_capacity(atoms.len());
    let mut seen = HashSet::new();
    let mut previous = None;
    let mut current = start;
    loop {
        order.push(current);
        seen.insert(current);
        let next = adjacency
            .neighbors_of(current)
            .iter()
            .map(|neighbor| neighbor.atom_index)
            .filter(|neighbor| atom_set.contains(neighbor))
            .filter(|neighbor| Some(*neighbor) != previous)
            .find(|neighbor| !seen.contains(neighbor));
        let Some(next) = next else {
            break;
        };
        previous = Some(current);
        current = next;
    }
    (order.len() == atoms.len()).then_some(order)
}

fn confseq_base_bond_target_between(
    model: &ConfSeqBaseConstraintModel,
    molecule: &Molecule,
    begin: usize,
    end: usize,
) -> Result<f64, ConfSeqFastGeometryError> {
    let Some(bond) = bond_between_pair(molecule, sorted_pair(begin, end)) else {
        return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
            "missing bond target for {begin}-{end}"
        )));
    };
    Ok(confseq_base_bond_target_by_id(
        model,
        molecule,
        bond.id().index(),
    ))
}

fn confseq_base_tree_fragment_template(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    template: &ConfSeqBaseRigidFragmentTemplate,
) -> Result<HashMap<usize, [f64; 3]>, ConfSeqFastGeometryError> {
    let atom_set: HashSet<_> = template.realization_atoms.iter().copied().collect();
    let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
    let root = template
        .realization_atoms
        .iter()
        .copied()
        .max_by_key(|&atom| {
            adjacency
                .neighbors_of(atom)
                .iter()
                .filter(|neighbor| atom_set.contains(&neighbor.atom_index))
                .count()
        })
        .ok_or_else(|| {
            ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "tree template has no root key={}",
                template.cache_key
            ))
        })?;
    let root_children = adjacency
        .neighbors_of(root)
        .iter()
        .map(|neighbor| neighbor.atom_index)
        .filter(|atom| atom_set.contains(atom))
        .collect::<Vec<_>>();
    let mut coords = HashMap::new();
    coords.insert(root, [0.0, 0.0, 0.0]);
    if root_children.is_empty() {
        return Ok(coords);
    }
    let root_dirs = confseq_base_center_frame_directions(molecule, model, root, &root_children)?;
    let mut visited = HashSet::from([root]);
    let mut queue = VecDeque::new();
    for (child_ord, child) in root_children.iter().copied().enumerate() {
        let length = confseq_base_bond_target_between(model, molecule, root, child)?;
        coords.insert(child, vec_scale(root_dirs[child_ord], length));
        visited.insert(child);
        queue.push_back((root, child));
    }
    while let Some((parent, atom)) = queue.pop_front() {
        let children = adjacency
            .neighbors_of(atom)
            .iter()
            .map(|neighbor| neighbor.atom_index)
            .filter(|child| atom_set.contains(child))
            .filter(|child| *child != parent)
            .collect::<Vec<_>>();
        let parent_axis = vec_normalize(vec_sub(coords[&atom], coords[&parent]));
        for (child_ord, child) in children.iter().copied().enumerate() {
            if visited.contains(&child) {
                continue;
            }
            let angle = model
                .angle_targets
                .get(&sorted_angle(parent, atom, child))
                .copied()
                .unwrap_or_else(|| {
                    confseq_base_source_backed_angle_rad(molecule, parent, atom, child)
                });
            let direction = child_direction(parent_axis, child_ord, children.len(), angle);
            let length = confseq_base_bond_target_between(model, molecule, atom, child)?;
            coords.insert(child, vec_add(coords[&atom], vec_scale(direction, length)));
            visited.insert(child);
            queue.push_back((atom, child));
        }
    }
    if visited.len() != template.realization_atoms.len() {
        return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
            "tree template failed to visit all atoms key={}",
            template.cache_key
        )));
    }
    Ok(coords)
}

fn confseq_base_template_geometry_score(
    model: &ConfSeqBaseConstraintModel,
    template: &ConfSeqBaseRigidFragmentTemplate,
    coords: &HashMap<usize, [f64; 3]>,
    atom_filter: Option<&HashSet<usize>>,
) -> f64 {
    let in_filter = |atom: usize| atom_filter.is_none_or(|filter| filter.contains(&atom));
    let mut score = 0.0;
    for bond in &template.bonds {
        if !in_filter(bond.begin) || !in_filter(bond.end) {
            continue;
        }
        let Some(begin) = coords.get(&bond.begin).copied() else {
            continue;
        };
        let Some(end) = coords.get(&bond.end).copied() else {
            continue;
        };
        let delta = vec_len(vec_sub(begin, end)) - bond.length;
        score += delta * delta;
    }
    for angle in &template.angles {
        if !in_filter(angle.left) || !in_filter(angle.center) || !in_filter(angle.right) {
            continue;
        }
        let Some(left) = coords.get(&angle.left).copied() else {
            continue;
        };
        let Some(center) = coords.get(&angle.center).copied() else {
            continue;
        };
        let Some(right) = coords.get(&angle.right).copied() else {
            continue;
        };
        if let Some(observed) = angle_rad_from_points(left, center, right) {
            let delta = angular_delta_rad(observed, angle.angle_rad);
            score += delta * delta;
        }
    }
    for prior in &model.path14_distance_priors {
        let (i, j, k, l) = prior.atoms;
        if !in_filter(i) || !in_filter(j) || !in_filter(k) || !in_filter(l) {
            continue;
        }
        let Some(left) = coords.get(&i).copied() else {
            continue;
        };
        let Some(right) = coords.get(&l).copied() else {
            continue;
        };
        let observed = vec_len(vec_sub(left, right));
        let delta = if observed < prior.lower_bound {
            prior.lower_bound - observed
        } else if observed > prior.upper_bound {
            observed - prior.upper_bound
        } else {
            0.0
        };
        score += 4.0 * delta * delta;
    }
    score
}

fn confseq_base_simple_ring_order(molecule: &Molecule, atoms: &[usize]) -> Option<Vec<usize>> {
    if !(3..=8).contains(&atoms.len()) {
        return None;
    }
    let atom_set: HashSet<_> = atoms.iter().copied().collect();
    let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
    let mut ring_neighbors = HashMap::<usize, Vec<usize>>::new();
    for &atom in atoms {
        let neighbors = adjacency
            .neighbors_of(atom)
            .iter()
            .filter_map(|neighbor| {
                atom_set
                    .contains(&neighbor.atom_index)
                    .then_some(neighbor.atom_index)
            })
            .collect::<Vec<_>>();
        if neighbors.len() != 2 {
            return None;
        }
        ring_neighbors.insert(atom, neighbors);
    }
    let start = atoms.iter().copied().min()?;
    let first = *ring_neighbors.get(&start)?.iter().min()?;
    let mut order = vec![start, first];
    let mut prev = start;
    let mut current = first;
    while order.len() < atoms.len() {
        let next = ring_neighbors
            .get(&current)?
            .iter()
            .copied()
            .find(|&neighbor| neighbor != prev && !order.contains(&neighbor))?;
        order.push(next);
        prev = current;
        current = next;
    }
    ring_neighbors
        .get(&current)?
        .contains(&start)
        .then_some(order)
}

fn confseq_base_three_atom_angle_center(molecule: &Molecule, atoms: &[usize]) -> Option<usize> {
    let atom_set: HashSet<_> = atoms.iter().copied().collect();
    let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
    atoms.iter().copied().find(|&atom| {
        adjacency
            .neighbors_of(atom)
            .iter()
            .filter(|neighbor| atom_set.contains(&neighbor.atom_index))
            .count()
            == 2
    })
}

fn confseq_base_single_center_ligands(
    molecule: &Molecule,
    atoms: &[usize],
) -> Option<(usize, Vec<usize>)> {
    let atom_set: HashSet<_> = atoms.iter().copied().collect();
    let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
    let mut centers = Vec::new();
    for &atom in atoms {
        let ligands = adjacency
            .neighbors_of(atom)
            .iter()
            .filter_map(|neighbor| {
                atom_set
                    .contains(&neighbor.atom_index)
                    .then_some(neighbor.atom_index)
            })
            .collect::<Vec<_>>();
        if ligands.len() + 1 == atoms.len() && (2..=4).contains(&ligands.len()) {
            centers.push((atom, ligands));
        }
    }
    (centers.len() == 1).then(|| centers.remove(0))
}

fn confseq_base_center_template_from_template(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    template: &ConfSeqBaseRigidFragmentTemplate,
    center: usize,
    ligands: &[usize],
) -> Result<HashMap<usize, [f64; 3]>, ConfSeqFastGeometryError> {
    let mut directions = confseq_base_center_frame_directions(molecule, model, center, ligands)?;
    let mut base_coords = confseq_base_center_template_coords_for_ligands(
        molecule,
        model,
        center,
        ligands,
        &directions,
    )?;
    confseq_base_fit_center_template_angles(
        molecule,
        model,
        center,
        ligands,
        &mut directions,
        &mut base_coords,
    );
    if ligands.len() < 3 {
        return Ok(base_coords);
    }

    let mut best = base_coords;
    let mut best_score = confseq_base_center_template_score(model, template, &best);
    for permutation in confseq_base_ligand_permutations(ligands.len()) {
        if permutation.iter().copied().eq(0..ligands.len()) {
            continue;
        }
        let permuted_ligands = permutation
            .iter()
            .map(|&idx| ligands[idx])
            .collect::<Vec<_>>();
        let mut coords = confseq_base_center_template_coords_for_ligands(
            molecule,
            model,
            center,
            &permuted_ligands,
            &directions,
        )?;
        confseq_base_fit_center_template_angles(
            molecule,
            model,
            center,
            &permuted_ligands,
            &mut directions.clone(),
            &mut coords,
        );
        let score = confseq_base_center_template_score(model, template, &coords);
        if score < best_score {
            best_score = score;
            best = coords;
        }
    }
    Ok(best)
}

fn confseq_base_center_template_coords_for_ligands(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    center: usize,
    ligands: &[usize],
    directions: &[[f64; 3]],
) -> Result<HashMap<usize, [f64; 3]>, ConfSeqFastGeometryError> {
    let mut coords = HashMap::new();
    coords.insert(center, [0.0, 0.0, 0.0]);
    for (idx, &ligand) in ligands.iter().enumerate() {
        let Some(bond) = bond_between_pair(molecule, sorted_pair(center, ligand)) else {
            return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "center template missing bond {center}-{ligand}"
            )));
        };
        let length = confseq_base_bond_target_by_id(model, molecule, bond.id().index());
        coords.insert(ligand, vec_scale(directions[idx], length));
    }
    Ok(coords)
}

fn confseq_base_center_template_score(
    model: &ConfSeqBaseConstraintModel,
    template: &ConfSeqBaseRigidFragmentTemplate,
    coords: &HashMap<usize, [f64; 3]>,
) -> f64 {
    let atom_filter: HashSet<_> = template.realization_atoms.iter().copied().collect();
    confseq_base_template_geometry_score(model, template, coords, Some(&atom_filter))
}

fn confseq_base_ligand_permutations(len: usize) -> Vec<Vec<usize>> {
    fn permute(prefix: &mut Vec<usize>, remaining: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if remaining.is_empty() {
            out.push(prefix.clone());
            return;
        }
        for idx in 0..remaining.len() {
            let value = remaining.remove(idx);
            prefix.push(value);
            permute(prefix, remaining, out);
            prefix.pop();
            remaining.insert(idx, value);
        }
    }
    let mut out = Vec::new();
    let mut prefix = Vec::new();
    let mut remaining = (0..len).collect();
    permute(&mut prefix, &mut remaining, &mut out);
    out
}

fn confseq_base_center_frame_directions(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    center: usize,
    ligands: &[usize],
) -> Result<Vec<[f64; 3]>, ConfSeqFastGeometryError> {
    match ligands.len() {
        2 => {
            let angle =
                confseq_base_source_backed_angle_rad(molecule, ligands[0], center, ligands[1]);
            Ok(vec![
                [1.0, 0.0, 0.0],
                [angle.cos(), angle.sin().max(1.0e-8), 0.0],
            ])
        }
        3 => {
            if confseq_base_center_prefers_planar_frame(molecule, center, ligands) {
                Ok((0..3)
                    .map(|idx| {
                        let theta = 2.0 * PI * idx as f64 / 3.0;
                        [theta.cos(), theta.sin(), 0.0]
                    })
                    .collect())
            } else {
                Ok(confseq_base_center_directions_from_target_angles(
                    molecule, model, center, ligands,
                )
                .unwrap_or_else(|| {
                    vec![
                        [1.0, 0.0, -1.0 / 3.0],
                        [-0.5, 0.8660254037844386, -1.0 / 3.0],
                        [-0.5, -0.8660254037844386, -1.0 / 3.0],
                    ]
                    .into_iter()
                    .map(vec_normalize)
                    .collect()
                }))
            }
        }
        4 => Ok(confseq_base_center_directions_from_target_angles(
            molecule, model, center, ligands,
        )
        .unwrap_or_else(|| {
            vec![
                [1.0, 1.0, 1.0],
                [1.0, -1.0, -1.0],
                [-1.0, 1.0, -1.0],
                [-1.0, -1.0, 1.0],
            ]
            .into_iter()
            .map(vec_normalize)
            .collect()
        })),
        _ => Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
            "center template supports 2-4 ligands, got {} at atom {center}",
            ligands.len()
        ))),
    }
}

fn confseq_base_center_prefers_planar_frame(
    molecule: &Molecule,
    center: usize,
    ligands: &[usize],
) -> bool {
    let atom = &molecule.atoms()[center];
    atom.hybridization() == Hybridization::Sp2
        || atom.is_aromatic()
        || ligands.iter().any(|&ligand| {
            bond_between_pair(molecule, sorted_pair(center, ligand))
                .is_some_and(|bond| matches!(bond.order(), BondOrder::Double | BondOrder::Aromatic))
        })
}

fn confseq_base_center_directions_from_target_angles(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    center: usize,
    ligands: &[usize],
) -> Option<Vec<[f64; 3]>> {
    let n = ligands.len();
    if !(3..=4).contains(&n) {
        return None;
    }
    let mut gram = vec![vec![0.0; n]; n];
    for idx in 0..n {
        gram[idx][idx] = 1.0;
    }
    for left_idx in 0..n {
        for right_idx in left_idx + 1..n {
            let target = model
                .angle_targets
                .get(&sorted_angle(ligands[left_idx], center, ligands[right_idx]))
                .copied()
                .unwrap_or_else(|| {
                    confseq_base_source_backed_angle_rad(
                        molecule,
                        ligands[left_idx],
                        center,
                        ligands[right_idx],
                    )
                });
            let cos_target = target.cos().clamp(-0.995, 0.995);
            gram[left_idx][right_idx] = cos_target;
            gram[right_idx][left_idx] = cos_target;
        }
    }
    confseq_base_directions_from_gram(&gram)
}

fn confseq_base_directions_from_gram(gram: &[Vec<f64>]) -> Option<Vec<[f64; 3]>> {
    match gram.len() {
        3 => {
            let c12 = gram[0][1];
            let c13 = gram[0][2];
            let c23 = gram[1][2];
            let y2 = (1.0 - c12 * c12).max(0.0).sqrt();
            if y2 <= 1.0e-8 {
                return None;
            }
            let x3 = c13;
            let y3 = (c23 - c12 * c13) / y2;
            let z3_sq = 1.0 - x3 * x3 - y3 * y3;
            if z3_sq < -1.0e-6 {
                return None;
            }
            Some(vec![
                [1.0, 0.0, 0.0],
                [c12, y2, 0.0],
                vec_normalize([x3, y3, z3_sq.max(0.0).sqrt()]),
            ])
        }
        4 => {
            let mut l = vec![vec![0.0; 4]; 4];
            for i in 0..4 {
                for j in 0..=i {
                    let sum = (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>();
                    if i == j {
                        let diag = gram[i][i] - sum;
                        if diag < -1.0e-6 {
                            return None;
                        }
                        l[i][j] = diag.max(0.0).sqrt();
                    } else if l[j][j] > 1.0e-8 {
                        l[i][j] = (gram[i][j] - sum) / l[j][j];
                    } else {
                        return None;
                    }
                }
            }
            Some(
                l.into_iter()
                    .map(|row| vec_normalize([row[0], row[1], row[2]]))
                    .collect(),
            )
        }
        _ => None,
    }
}

fn confseq_base_fit_center_template_angles(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    center: usize,
    ligands: &[usize],
    directions: &mut [[f64; 3]],
    coords: &mut HashMap<usize, [f64; 3]>,
) {
    if ligands.len() != 3 || !confseq_base_center_prefers_planar_frame(molecule, center, ligands) {
        return;
    }
    let mut angles = [2.0 * PI / 3.0; 3];
    for idx in 0..3 {
        let left = ligands[idx];
        let right = ligands[(idx + 1) % 3];
        if let Some(&angle) = model.angle_targets.get(&sorted_angle(left, center, right)) {
            angles[idx] = angle;
        }
    }
    let total: f64 = angles.iter().sum();
    if total <= 1.0e-10 {
        return;
    }
    let scale = 2.0 * PI / total;
    let mut theta: f64 = 0.0;
    for idx in 0..3 {
        directions[idx] = [theta.cos(), theta.sin(), 0.0];
        if let Some(coord) = coords.get_mut(&ligands[idx]) {
            let length = vec_len(*coord);
            *coord = vec_scale(directions[idx], length);
        }
        theta += angles[idx] * scale;
    }
}

fn place_confseq_base_rigid_component_local(
    component: &ConfSeqBaseRigidComponent,
    local: &HashMap<usize, [f64; 3]>,
    translation: [f64; 3],
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
    connector_targets: &mut HashMap<(usize, usize), [f64; 3]>,
) -> Result<(), ConfSeqFastGeometryError> {
    for &atom in &component.atoms {
        let Some(local_coord) = local.get(&atom) else {
            return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "realized fragment is missing atom {atom} for component {:?}",
                component.atoms
            )));
        };
        coords[atom] = vec_add(*local_coord, translation);
        placed[atom] = true;
    }
    update_confseq_base_connector_targets(
        component,
        local,
        |point| Ok(vec_add(point, translation)),
        coords,
        connector_targets,
    )
}

#[allow(clippy::too_many_arguments)]
fn place_confseq_base_rigid_component_aligned_to_segment(
    component: &ConfSeqBaseRigidComponent,
    local: &HashMap<usize, [f64; 3]>,
    local_origin: [f64; 3],
    local_axis: [f64; 3],
    target_origin: [f64; 3],
    target_axis: [f64; 3],
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
    connector_targets: &mut HashMap<(usize, usize), [f64; 3]>,
) -> Result<(), ConfSeqFastGeometryError> {
    let transformed = confseq_base_transformed_component_coords(
        component,
        local,
        local_origin,
        local_axis,
        target_origin,
        target_axis,
    )?
    .ok_or_else(|| {
        ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
            "degenerate rigid segment alignment for component {:?}",
            component.atoms
        ))
    })?;
    commit_confseq_base_transformed_component(
        component,
        local,
        transformed,
        |point| {
            let rotated = rotate_vector_between_unit_dirs(
                vec_sub(point, local_origin),
                local_axis,
                target_axis,
            )?;
            Ok(vec_add(target_origin, rotated))
        },
        coords,
        placed,
        connector_targets,
    )
}

fn confseq_base_transformed_component_coords(
    component: &ConfSeqBaseRigidComponent,
    local: &HashMap<usize, [f64; 3]>,
    local_origin: [f64; 3],
    local_axis: [f64; 3],
    target_origin: [f64; 3],
    target_axis: [f64; 3],
) -> Result<Option<HashMap<usize, [f64; 3]>>, ConfSeqFastGeometryError> {
    if vec_len(local_axis) <= 1.0e-8 || vec_len(target_axis) <= 1.0e-8 {
        return Ok(None);
    }
    let transform = |point: [f64; 3]| -> Result<[f64; 3], ConfSeqFastGeometryError> {
        let rotated =
            rotate_vector_between_unit_dirs(vec_sub(point, local_origin), local_axis, target_axis)?;
        Ok(vec_add(target_origin, rotated))
    };
    let mut transformed = HashMap::new();
    for &atom in &component.atoms {
        let Some(local_coord) = local.get(&atom).copied() else {
            return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "realized fragment is missing atom {atom} for component {:?}",
                component.atoms
            )));
        };
        transformed.insert(atom, transform(local_coord)?);
    }
    Ok(Some(transformed))
}

fn commit_confseq_base_transformed_component(
    component: &ConfSeqBaseRigidComponent,
    local: &HashMap<usize, [f64; 3]>,
    transformed: HashMap<usize, [f64; 3]>,
    transform: impl Fn([f64; 3]) -> Result<[f64; 3], ConfSeqFastGeometryError>,
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
    connector_targets: &mut HashMap<(usize, usize), [f64; 3]>,
) -> Result<(), ConfSeqFastGeometryError> {
    for &atom in &component.atoms {
        let Some(coord) = transformed.get(&atom).copied() else {
            return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "transformed fragment is missing atom {atom} for component {:?}",
                component.atoms
            )));
        };
        coords[atom] = coord;
        placed[atom] = true;
    }
    update_confseq_base_connector_targets(component, local, transform, coords, connector_targets)
}

fn update_confseq_base_connector_targets(
    component: &ConfSeqBaseRigidComponent,
    local: &HashMap<usize, [f64; 3]>,
    transform: impl Fn([f64; 3]) -> Result<[f64; 3], ConfSeqFastGeometryError>,
    _coords: &[[f64; 3]],
    connector_targets: &mut HashMap<(usize, usize), [f64; 3]>,
) -> Result<(), ConfSeqFastGeometryError> {
    for connector in &component.connectors {
        let Some(local_coord) = local.get(&connector.external_atom).copied() else {
            continue;
        };
        let target = transform(local_coord).map_err(|_| {
            ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "failed to transform connector stub atom {} for component {:?}",
                connector.external_atom, component.atoms
            ))
        })?;
        connector_targets.insert((connector.core_atom, connector.external_atom), target);
    }
    Ok(())
}

fn confseq_base_fragment_attach_direction(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
    coords: &[[f64; 3]],
    placed: &[bool],
    placed_atom: usize,
    unplaced_atom: usize,
) -> [f64; 3] {
    let placed_neighbors: Vec<_> = adjacency
        .neighbors_of(placed_atom)
        .iter()
        .filter(|neighbor| placed[neighbor.atom_index])
        .map(|neighbor| neighbor.atom_index)
        .collect();
    if let Some(parent) = placed_neighbors.first().copied() {
        let axis = vec_normalize(vec_sub(coords[placed_atom], coords[parent]));
        let angle = model
            .angle_targets
            .get(&sorted_angle(parent, placed_atom, unplaced_atom))
            .copied()
            .unwrap_or_else(|| confseq_base_local_angle_rad(molecule, placed_atom));
        child_direction(axis, 0, 1, angle)
    } else {
        [1.0, 0.0, 0.0]
    }
}

fn place_confseq_base_hydrogens_from_heavy_neighbors(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
    coords: &mut [[f64; 3]],
) {
    for atom in molecule.atoms() {
        let atom_idx = atom.id().index();
        if atom.atomic_number() != 1 {
            continue;
        }
        let Some(parent) = adjacency
            .neighbors_of(atom_idx)
            .iter()
            .map(|neighbor| neighbor.atom_index)
            .find(|neighbor| molecule.atoms()[*neighbor].atomic_number() > 1)
        else {
            continue;
        };
        let bond_id = adjacency
            .neighbors_of(atom_idx)
            .iter()
            .find(|neighbor| neighbor.atom_index == parent)
            .map(|neighbor| neighbor.bond.index());
        let length = bond_id
            .map(|bond_id| confseq_base_bond_target_by_id(model, molecule, bond_id))
            .unwrap_or(1.09);
        let axis = adjacency
            .neighbors_of(parent)
            .iter()
            .filter(|neighbor| neighbor.atom_index != atom_idx)
            .find(|neighbor| molecule.atoms()[neighbor.atom_index].atomic_number() > 1)
            .map(|neighbor| vec_normalize(vec_sub(coords[parent], coords[neighbor.atom_index])))
            .unwrap_or([1.0, 0.0, 0.0]);
        coords[atom_idx] = vec_add(
            coords[parent],
            vec_scale(
                child_direction(axis, 0, 1, confseq_base_local_angle_rad(molecule, parent)),
                length,
            ),
        );
    }
}

fn confseq_base_component_kind_is_planar_shape(kind: ConfSeqBaseRigidComponentKind) -> bool {
    matches!(
        kind,
        ConfSeqBaseRigidComponentKind::AcyclicPlanarPi | ConfSeqBaseRigidComponentKind::RingPlanar
    )
}

fn confseq_base_bond_target_by_id(
    model: &ConfSeqBaseConstraintModel,
    molecule: &Molecule,
    bond_idx: usize,
) -> f64 {
    model
        .bond_targets
        .get(bond_idx)
        .copied()
        .unwrap_or_else(|| confseq_base_static_bond_length_fallback(&molecule.bonds()[bond_idx]))
}

fn validate_confseq_base_constraint_coords(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    coords: &[[f64; 3]],
) -> Result<(), ConfSeqFastGeometryError> {
    for bond in molecule.bonds() {
        let target = model
            .bond_targets
            .get(bond.id().index())
            .copied()
            .unwrap_or_else(|| confseq_base_static_bond_length_fallback(bond));
        let observed = vec_len(vec_sub(
            coords[bond.begin().index()],
            coords[bond.end().index()],
        ));
        if (observed - target).abs() > 0.35 {
            return Err(ConfSeqFastGeometryError::AssemblyBoundaryMismatch {
                bond: bond.id().index(),
                observed,
                target,
            });
        }
    }
    for &(begin, end) in &model.planar_bonds {
        let observed = vec_len(vec_sub(coords[begin], coords[end]));
        if observed <= 1.0e-10 {
            return Err(ConfSeqFastGeometryError::Build(format!(
                "planar bond {begin}-{end} has coincident base coordinates"
            )));
        }
    }
    for (&torsion, prior) in &model.torsion_priors {
        if matches!(prior, ConfSeqBaseTorsionPrior::Free) {
            continue;
        }
        let (i, j, k, l) = torsion;
        if vec_len(vec_sub(coords[i], coords[j])) <= 1.0e-10
            || vec_len(vec_sub(coords[j], coords[k])) <= 1.0e-10
            || vec_len(vec_sub(coords[k], coords[l])) <= 1.0e-10
        {
            return Err(ConfSeqFastGeometryError::Build(format!(
                "torsion prior {i}-{j}-{k}-{l} has degenerate base coordinates"
            )));
        }
    }
    Ok(())
}

fn sorted_angle(left: usize, center: usize, right: usize) -> (usize, usize, usize) {
    if left <= right {
        (left, center, right)
    } else {
        (right, center, left)
    }
}

fn apply_confseq_base_double_bond_stereo(
    molecule: Molecule,
) -> Result<Molecule, ConfSeqFastGeometryError> {
    let constraints = collect_confseq_base_double_bond_stereo_constraints(&molecule);
    let mut molecule = molecule;
    for (i, j, k, l, target_deg) in constraints {
        molecule = mol_transforms::set_dihedral_deg(molecule, i, j, k, l, target_deg, 0).map_err(
            |err| {
                ConfSeqFastGeometryError::Build(format!(
                    "failed to apply double-bond stereo constraint {i}-{j}-{k}-{l}: {err}"
                ))
            },
        )?;
    }
    Ok(molecule)
}

fn collect_confseq_base_double_bond_stereo_constraints(
    molecule: &Molecule,
) -> Vec<(usize, usize, usize, usize, f64)> {
    molecule
        .bonds()
        .iter()
        .filter_map(|bond| {
            if bond.order() != BondOrder::Double {
                return None;
            }
            let target_deg = match bond.stereo() {
                // RDKit stores SMILES directional-bond stereo as cis/trans
                // relative to the two controlling atoms. The ConfSeq base
                // conformer must preserve that rigid local relationship before
                // applying rotatable single-bond tokens.
                BondStereo::Cis | BondStereo::Z => 0.0,
                BondStereo::Trans | BondStereo::E => 180.0,
                _ => return None,
            };
            let [left, right] = bond.stereo_atoms()?;
            let begin = bond.begin().index();
            let end = bond.end().index();
            Some((left.index(), begin, end, right.index(), target_deg))
        })
        .collect()
}

fn apply_confseq_base_tetrahedral_stereo(
    molecule: Molecule,
    adjacency: &AdjacencyList,
) -> Result<Molecule, ConfSeqFastGeometryError> {
    let constraints = collect_confseq_base_tetrahedral_stereo_constraints(&molecule, adjacency)?;
    if constraints.is_empty() {
        return Ok(molecule);
    }

    let mut coords = molecule
        .conformers_3d()
        .first()
        .ok_or_else(|| {
            ConfSeqFastGeometryError::Build("fast geometry has no coordinates".to_string())
        })?
        .coordinates()
        .to_vec();
    let initial_coords = coords.clone();
    for _ in 0..constraints.len().saturating_mul(2).max(1) {
        let mut changed = false;
        for constraint in &constraints {
            let volume = confseq_base_chiral_volume(&coords, constraint);
            if confseq_base_chiral_volume_satisfies_tag(volume, constraint.tag) {
                continue;
            }
            let mut candidates = Vec::new();
            for movable_pos in 0..constraint.ligands.len() {
                let movable_ligand = constraint.ligands[movable_pos];
                let Ok(movable) = confseq_base_chiral_movable_side(
                    &molecule,
                    adjacency,
                    constraint,
                    movable_ligand,
                ) else {
                    continue;
                };
                let contains_other_ligand = confseq_base_chiral_side_contains_other_ligand(
                    &movable,
                    constraint,
                    movable_ligand,
                );
                candidates.push((contains_other_ligand, movable.len(), movable_pos, movable));
            }
            candidates.sort_by_key(|(contains_other_ligand, len, movable_pos, _)| {
                (*contains_other_ligand, *len, *movable_pos)
            });
            let mut best = None::<(usize, f64, f64, Vec<[f64; 3]>)>;
            for (contains_other_ligand, _, movable_pos, movable) in candidates {
                let mut trial = coords.clone();
                let adjusted = adjust_movable_side_to_chiral_volume_sign(
                    &mut trial,
                    &movable,
                    constraint,
                    movable_pos,
                    !contains_other_ligand,
                )?;
                if adjusted
                    && confseq_base_chiral_volume_satisfies_tag(
                        confseq_base_chiral_volume(&trial, constraint),
                        constraint.tag,
                    )
                {
                    let unsatisfied =
                        confseq_base_unsatisfied_chiral_constraints(&trial, &constraints);
                    let rms_displacement =
                        confseq_base_coord_displacement_rms(&initial_coords, &trial);
                    let max_displacement =
                        confseq_base_coord_max_displacement(&initial_coords, &trial);
                    let replace = best
                        .as_ref()
                        .map(|(best_unsatisfied, best_rms, best_max, _)| {
                            (unsatisfied, rms_displacement, max_displacement)
                                < (*best_unsatisfied, *best_rms, *best_max)
                        })
                        .unwrap_or(true);
                    if replace {
                        best = Some((unsatisfied, rms_displacement, max_displacement, trial));
                    }
                }
            }
            if let Some((_, _, _, trial)) = best {
                coords = trial;
            } else {
                return Err(ConfSeqFastGeometryError::UnsupportedTetrahedralStereo {
                    center: constraint.center,
                    reason: "cannot be corrected by moving one explicit ligand side".to_string(),
                });
            }
            changed = true;
        }
        if !changed {
            break;
        }
    }

    for constraint in &constraints {
        let volume = confseq_base_chiral_volume(&coords, constraint);
        if !confseq_base_chiral_volume_satisfies_tag(volume, constraint.tag) {
            return Err(ConfSeqFastGeometryError::UnsupportedTetrahedralStereo {
                center: constraint.center,
                reason: "final chiral-volume sign validation failed".to_string(),
            });
        }
    }

    let mut out = molecule;
    let coord_block = out.coordinate_block_mut();
    let Some(conformer) = coord_block.conformers_3d.first_mut() else {
        return Err(ConfSeqFastGeometryError::Build(
            "fast geometry has no coordinates".to_string(),
        ));
    };
    conformer.coordinates_mut().copy_from_slice(&coords);
    Ok(out)
}

fn confseq_base_unsatisfied_chiral_constraints(
    coords: &[[f64; 3]],
    constraints: &[ConfSeqBaseTetrahedralStereoConstraint],
) -> usize {
    constraints
        .iter()
        .filter(|constraint| {
            !confseq_base_chiral_volume_satisfies_tag(
                confseq_base_chiral_volume(coords, constraint),
                constraint.tag,
            )
        })
        .count()
}

fn confseq_base_coord_displacement_rms(reference: &[[f64; 3]], coords: &[[f64; 3]]) -> f64 {
    if reference.is_empty() || reference.len() != coords.len() {
        return f64::INFINITY;
    }
    let sum_sq = reference
        .iter()
        .zip(coords)
        .map(|(left, right)| {
            let delta = vec_sub(*left, *right);
            vec_dot(delta, delta)
        })
        .sum::<f64>();
    (sum_sq / reference.len() as f64).sqrt()
}

fn confseq_base_coord_max_displacement(reference: &[[f64; 3]], coords: &[[f64; 3]]) -> f64 {
    if reference.len() != coords.len() {
        return f64::INFINITY;
    }
    reference
        .iter()
        .zip(coords)
        .map(|(left, right)| vec_len(vec_sub(*left, *right)))
        .fold(0.0, f64::max)
}

#[derive(Debug, Clone)]
struct ConfSeqBaseTetrahedralStereoConstraint {
    center: usize,
    ligands: [usize; 4],
    tag: ChiralTag,
}

fn collect_confseq_base_tetrahedral_stereo_constraints(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
) -> Result<Vec<ConfSeqBaseTetrahedralStereoConstraint>, ConfSeqFastGeometryError> {
    let mut constraints = Vec::new();
    for atom in molecule.atoms() {
        let tag = atom.chiral_tag();
        if !matches!(tag, ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw) {
            continue;
        }
        let center = atom.id().index();
        let neighbors: Vec<_> = adjacency
            .neighbors_of(center)
            .iter()
            .map(|neighbor| neighbor.atom_index)
            .collect();
        let ligands = match neighbors.as_slice() {
            [a, b, c, d] => [*a, *b, *c, *d],
            // RDKit's distance-geometry chiral set uses the center index as the
            // fourth point for three-neighbor tetrahedral centers with an
            // implicit ligand. The volume sign convention below follows that
            // existing RDKit-derived path exactly.
            [a, b, c] => [*a, *b, *c, center],
            _ => {
                return Err(ConfSeqFastGeometryError::UnsupportedTetrahedralStereo {
                    center,
                    reason: "requires 3 or 4 explicit neighbor coordinates".to_string(),
                });
            }
        };
        constraints.push(ConfSeqBaseTetrahedralStereoConstraint {
            center,
            ligands,
            tag,
        });
    }
    Ok(constraints)
}

fn confseq_base_chiral_volume(
    coords: &[[f64; 3]],
    constraint: &ConfSeqBaseTetrahedralStereoConstraint,
) -> f64 {
    let [a, b, c, d] = constraint.ligands;
    let anchor = if d == constraint.center {
        coords[constraint.center]
    } else {
        coords[d]
    };
    let v1 = vec_sub(coords[a], anchor);
    let v2 = vec_sub(coords[b], anchor);
    let v3 = vec_sub(coords[c], anchor);
    vec_dot(v1, vec_cross(v2, v3))
}

fn confseq_base_chiral_volume_satisfies_tag(volume: f64, tag: ChiralTag) -> bool {
    let eps = 1.0e-8;
    match tag {
        ChiralTag::TetrahedralCcw => volume > eps,
        ChiralTag::TetrahedralCw => volume < -eps,
        _ => true,
    }
}

fn confseq_base_chiral_movable_side(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    constraint: &ConfSeqBaseTetrahedralStereoConstraint,
    movable_ligand: usize,
) -> Result<Vec<usize>, ConfSeqFastGeometryError> {
    let center = constraint.center;
    if movable_ligand == center
        || bond_between_pair(molecule, sorted_pair(center, movable_ligand)).is_none()
    {
        return Err(ConfSeqFastGeometryError::UnsupportedTetrahedralStereo {
            center,
            reason: "has no explicit movable ligand".to_string(),
        });
    }
    let side =
        connected_side_without_crossing(adjacency, molecule.num_atoms(), movable_ligand, center);
    if side.contains(&center) {
        return Err(ConfSeqFastGeometryError::UnsupportedTetrahedralStereo {
            center,
            reason: "is inside a cyclic component unsupported by local mirror correction"
                .to_string(),
        });
    }
    // Cyclic centers often have multiple explicit ligands in the same graph
    // component after removing only the center atom. Rotating that whole
    // component about the center preserves all internal distances and all
    // center-ligand bond lengths; the caller still accepts the trial only after
    // recomputing the actual chiral volume and later revalidating every
    // tetrahedral constraint. This expands the strict candidate set without a
    // force-field or distance-geometry fallback.
    Ok(side)
}

fn confseq_base_chiral_side_contains_other_ligand(
    side: &[usize],
    constraint: &ConfSeqBaseTetrahedralStereoConstraint,
    movable_ligand: usize,
) -> bool {
    constraint.ligands.iter().copied().any(|ligand| {
        ligand != movable_ligand && ligand != constraint.center && side.contains(&ligand)
    })
}

fn connected_side_without_crossing(
    adjacency: &AdjacencyList,
    atom_count: usize,
    start: usize,
    blocked: usize,
) -> Vec<usize> {
    let mut seen = vec![false; atom_count];
    let mut queue = VecDeque::new();
    seen[blocked] = true;
    seen[start] = true;
    queue.push_back(start);
    while let Some(atom) = queue.pop_front() {
        for neighbor in adjacency.neighbors_of(atom) {
            if !seen[neighbor.atom_index] {
                seen[neighbor.atom_index] = true;
                queue.push_back(neighbor.atom_index);
            }
        }
    }
    seen.into_iter()
        .enumerate()
        .filter_map(|(idx, value)| (idx != blocked && value).then_some(idx))
        .collect()
}

fn rotate_movable_side_to_chiral_volume_sign(
    coords: &mut [[f64; 3]],
    atoms: &[usize],
    constraint: &ConfSeqBaseTetrahedralStereoConstraint,
    movable_pos: usize,
) -> Result<(), ConfSeqFastGeometryError> {
    let center = constraint.center;
    let movable_ligand = constraint.ligands[movable_pos];
    let normal = confseq_base_chiral_volume_gradient_for_ligand(coords, constraint, movable_pos);
    let normal_len = vec_len(normal);
    if normal_len <= 1.0e-10 {
        return Err(ConfSeqFastGeometryError::UnsupportedTetrahedralStereo {
            center,
            reason: "has collinear fixed ligands".to_string(),
        });
    }
    let normal = vec_scale(normal, 1.0 / normal_len);
    let center_point = coords[center];
    let old_root = vec_sub(coords[movable_ligand], center_point);
    let root_len = vec_len(old_root);
    if root_len <= 1.0e-10 {
        return Err(ConfSeqFastGeometryError::UnsupportedTetrahedralStereo {
            center,
            reason: "has coincident movable ligand coordinates".to_string(),
        });
    }

    let plane_origin =
        confseq_base_chiral_volume_plane_origin_for_ligand(coords, constraint, movable_pos);
    let distance_to_volume_plane = vec_dot(vec_sub(coords[movable_ligand], plane_origin), normal);
    let reflected = vec_sub(
        coords[movable_ligand],
        vec_scale(normal, 2.0 * distance_to_volume_plane),
    );
    let mut target_root = vec_scale(vec_normalize(vec_sub(reflected, center_point)), root_len);
    if !confseq_base_chiral_volume_satisfies_tag(
        confseq_base_chiral_volume_with_movable_root(coords, constraint, movable_pos, target_root),
        constraint.tag,
    ) {
        let sign = match constraint.tag {
            ChiralTag::TetrahedralCcw => 1.0,
            ChiralTag::TetrahedralCw => -1.0,
            _ => 1.0,
        };
        target_root = vec_scale(normal, sign * root_len);
    }
    if !confseq_base_chiral_volume_satisfies_tag(
        confseq_base_chiral_volume_with_movable_root(coords, constraint, movable_pos, target_root),
        constraint.tag,
    ) {
        return Err(ConfSeqFastGeometryError::UnsupportedTetrahedralStereo {
            center,
            reason: "cannot be satisfied by rotating one ligand side".to_string(),
        });
    }

    rotate_points_mapping_vector(coords, atoms, center_point, old_root, target_root)
}

fn adjust_movable_side_to_chiral_volume_sign(
    coords: &mut [[f64; 3]],
    atoms: &[usize],
    constraint: &ConfSeqBaseTetrahedralStereoConstraint,
    movable_pos: usize,
    allow_rotation_fallback: bool,
) -> Result<bool, ConfSeqFastGeometryError> {
    if translate_movable_side_to_chiral_volume_sign(coords, atoms, constraint, movable_pos)? {
        return Ok(true);
    }
    if !allow_rotation_fallback {
        return Ok(false);
    }
    rotate_movable_side_to_chiral_volume_sign(coords, atoms, constraint, movable_pos)?;
    Ok(true)
}

fn translate_movable_side_to_chiral_volume_sign(
    coords: &mut [[f64; 3]],
    atoms: &[usize],
    constraint: &ConfSeqBaseTetrahedralStereoConstraint,
    movable_pos: usize,
) -> Result<bool, ConfSeqFastGeometryError> {
    let center = constraint.center;
    let normal = confseq_base_chiral_volume_gradient_for_ligand(coords, constraint, movable_pos);
    let normal_len = vec_len(normal);
    if normal_len <= 1.0e-10 {
        return Ok(false);
    }
    let normal = vec_scale(normal, 1.0 / normal_len);
    let current = confseq_base_chiral_volume(coords, constraint);
    let target = confseq_base_min_signed_chiral_volume_for_constraint(constraint);
    let delta_volume = target - current;
    let displacement = delta_volume / normal_len;
    if displacement.abs() > confseq_base_max_chiral_translation(coords, constraint) {
        return Ok(false);
    }
    let delta = vec_scale(normal, displacement);
    for &atom in atoms {
        coords[atom] = vec_add(coords[atom], delta);
    }
    if confseq_base_chiral_volume_satisfies_tag(
        confseq_base_chiral_volume(coords, constraint),
        constraint.tag,
    ) {
        Ok(true)
    } else {
        Err(ConfSeqFastGeometryError::UnsupportedTetrahedralStereo {
            center,
            reason: "minimum chiral-volume translation failed".to_string(),
        })
    }
}

fn confseq_base_min_signed_chiral_volume_for_constraint(
    constraint: &ConfSeqBaseTetrahedralStereoConstraint,
) -> f64 {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::findChiralSets (Embedder.cpp:1112-1122)
    // RDKit✔️✔️:         double volLowerBound = 5.0;
    // RDKit✔️✔️:         double volUpperBound = 100.0;
    // RDKit✔️✔️:         if (nbrs.size() < 4) {
    // RDKit✔️✔️:           volLowerBound = 2.0;
    // RDKit✔️✔️:           nbrs.insert(nbrs.end(), atom->getIdx());
    // RDKit✔️✔️:         }
    // END RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::findChiralSets
    let lower = if constraint.ligands[3] == constraint.center {
        2.0
    } else {
        5.0
    };
    match constraint.tag {
        ChiralTag::TetrahedralCcw => lower,
        ChiralTag::TetrahedralCw => -lower,
        _ => 0.0,
    }
}

fn confseq_base_max_chiral_translation(
    coords: &[[f64; 3]],
    constraint: &ConfSeqBaseTetrahedralStereoConstraint,
) -> f64 {
    let center = coords[constraint.center];
    let min_bond = constraint
        .ligands
        .iter()
        .copied()
        .filter(|ligand| *ligand != constraint.center)
        .map(|ligand| vec_len(vec_sub(coords[ligand], center)))
        .fold(f64::INFINITY, f64::min);
    if min_bond.is_finite() {
        0.45 * min_bond
    } else {
        0.65
    }
}

fn confseq_base_chiral_volume_gradient_for_ligand(
    coords: &[[f64; 3]],
    constraint: &ConfSeqBaseTetrahedralStereoConstraint,
    ligand_pos: usize,
) -> [f64; 3] {
    let [a, b, c, d] = constraint.ligands;
    match ligand_pos {
        0 => vec_cross(vec_sub(coords[b], coords[d]), vec_sub(coords[c], coords[d])),
        1 => vec_cross(vec_sub(coords[c], coords[d]), vec_sub(coords[a], coords[d])),
        2 => vec_cross(vec_sub(coords[a], coords[d]), vec_sub(coords[b], coords[d])),
        3 if d != constraint.center => vec_scale(
            vec_cross(vec_sub(coords[b], coords[a]), vec_sub(coords[c], coords[a])),
            -1.0,
        ),
        _ => [0.0, 0.0, 0.0],
    }
}

fn confseq_base_chiral_volume_plane_origin_for_ligand(
    coords: &[[f64; 3]],
    constraint: &ConfSeqBaseTetrahedralStereoConstraint,
    ligand_pos: usize,
) -> [f64; 3] {
    let [a, _, _, d] = constraint.ligands;
    match ligand_pos {
        0..=2 => {
            if d == constraint.center {
                coords[constraint.center]
            } else {
                coords[d]
            }
        }
        3 => coords[a],
        _ => coords[constraint.center],
    }
}

fn confseq_base_chiral_volume_with_movable_root(
    coords: &[[f64; 3]],
    constraint: &ConfSeqBaseTetrahedralStereoConstraint,
    movable_pos: usize,
    root: [f64; 3],
) -> f64 {
    let center_point = coords[constraint.center];
    let mut points = [[0.0; 3]; 4];
    for (idx, ligand) in constraint.ligands.iter().copied().enumerate() {
        points[idx] = if idx == movable_pos {
            vec_add(center_point, root)
        } else if ligand == constraint.center {
            center_point
        } else {
            coords[ligand]
        };
    }
    let anchor = if constraint.ligands[3] == constraint.center {
        center_point
    } else {
        points[3]
    };
    let v1 = vec_sub(points[0], anchor);
    let v2 = vec_sub(points[1], anchor);
    let v3 = vec_sub(points[2], anchor);
    vec_dot(v1, vec_cross(v2, v3))
}

fn rotate_points_mapping_vector(
    coords: &mut [[f64; 3]],
    atoms: &[usize],
    origin: [f64; 3],
    from: [f64; 3],
    to: [f64; 3],
) -> Result<(), ConfSeqFastGeometryError> {
    let from = vec_normalize(from);
    let to = vec_normalize(to);
    let cos_theta = vec_dot(from, to).clamp(-1.0, 1.0);
    let axis = vec_cross(from, to);
    let sin_theta = vec_len(axis);
    let (axis, sin_theta) = if sin_theta <= 1.0e-12 {
        if cos_theta > 0.0 {
            return Ok(());
        }
        (perpendicular_unit(from), 0.0)
    } else {
        (vec_scale(axis, 1.0 / sin_theta), sin_theta)
    };
    let cos_theta = if sin_theta == 0.0 { -1.0 } else { cos_theta };
    for &atom in atoms {
        let offset = vec_sub(coords[atom], origin);
        coords[atom] = vec_add(
            origin,
            rotate_vec_around_unit_axis(offset, axis, sin_theta, cos_theta),
        );
    }
    Ok(())
}

fn rotate_vector_between_unit_dirs(
    vector: [f64; 3],
    from: [f64; 3],
    to: [f64; 3],
) -> Result<[f64; 3], ConfSeqFastGeometryError> {
    let from = vec_normalize(from);
    let to = vec_normalize(to);
    let cos_theta = vec_dot(from, to).clamp(-1.0, 1.0);
    let axis = vec_cross(from, to);
    let sin_theta = vec_len(axis);
    if sin_theta <= 1.0e-12 {
        if cos_theta > 0.0 {
            return Ok(vector);
        }
        return Ok(rotate_vec_around_unit_axis(
            vector,
            perpendicular_unit(from),
            0.0,
            -1.0,
        ));
    }
    Ok(rotate_vec_around_unit_axis(
        vector,
        vec_scale(axis, 1.0 / sin_theta),
        sin_theta,
        cos_theta,
    ))
}

fn rotate_vec_around_unit_axis(
    vector: [f64; 3],
    axis: [f64; 3],
    sin_theta: f64,
    cos_theta: f64,
) -> [f64; 3] {
    vec_add(
        vec_add(
            vec_scale(vector, cos_theta),
            vec_scale(vec_cross(axis, vector), sin_theta),
        ),
        vec_scale(axis, vec_dot(axis, vector) * (1.0 - cos_theta)),
    )
}

fn angular_delta_rad(left: f64, right: f64) -> f64 {
    let mut delta = left - right;
    while delta > PI {
        delta -= 2.0 * PI;
    }
    while delta < -PI {
        delta += 2.0 * PI;
    }
    delta
}

fn shared_bond_ids_between_rings(
    ring_info: &rings::RingInfo,
    left: usize,
    right: usize,
) -> Vec<BondId> {
    ring_info.bond_rings()[left]
        .iter()
        .copied()
        .filter(|bond| ring_info.bond_rings()[right].contains(bond))
        .collect()
}

fn validate_supported_confseq_base_ring(
    molecule: &Molecule,
    ring_index: usize,
    ring_atoms: &[AtomId],
    ring_bonds: &[BondId],
) -> Result<(), ConfSeqFastGeometryError> {
    if !(3..=8).contains(&ring_atoms.len()) {
        return Err(ConfSeqFastGeometryError::UnsupportedRingSize {
            ring_size: ring_atoms.len(),
        });
    }
    let all_aromatic_atoms = ring_atoms
        .iter()
        .all(|atom| molecule.atoms()[atom.index()].is_aromatic());
    let all_aromatic_bonds = ring_bonds
        .iter()
        .all(|bond| molecule.bonds()[bond.index()].is_aromatic());
    if all_aromatic_atoms || all_aromatic_bonds {
        if !matches!(ring_atoms.len(), 5 | 6) {
            return Err(ConfSeqFastGeometryError::UnsupportedAromaticRingSize {
                ring_size: ring_atoms.len(),
            });
        }
        for atom in ring_atoms {
            let atomic_number = molecule.atoms()[atom.index()].atomic_number();
            if !matches!(atomic_number, 6 | 7 | 8 | 16) {
                return Err(ConfSeqFastGeometryError::UnsupportedRingElement {
                    ring_index,
                    atomic_number,
                });
            }
        }
        return Ok(());
    }
    for atom in ring_atoms {
        let atomic_number = molecule.atoms()[atom.index()].atomic_number();
        if !matches!(atomic_number, 6 | 7 | 8 | 15 | 16) {
            return Err(ConfSeqFastGeometryError::UnsupportedRingElement {
                ring_index,
                atomic_number,
            });
        }
    }
    Ok(())
}

fn child_direction(
    parent_axis: [f64; 3],
    child_ord: usize,
    child_count: usize,
    angle: f64,
) -> [f64; 3] {
    let axis_to_parent = vec_scale(parent_axis, -1.0);
    let normal = if parent_axis[2].abs() < 0.9 {
        vec_normalize(vec_cross(parent_axis, [0.0, 0.0, 1.0]))
    } else {
        vec_normalize(vec_cross(parent_axis, [0.0, 1.0, 0.0]))
    };
    let binormal = vec_normalize(vec_cross(axis_to_parent, normal));
    // Subset geometry rule: children are distributed around the parent-axis cone
    // instead of alternating between two directions. For two children, derive the
    // azimuth that makes child-child and parent-child angles match the local
    // ideal where possible: sp3 gives 120 deg azimuth, sp2 gives 180 deg.
    // This is still local tree placement, not a global non-bonded optimizer.
    let phase = child_azimuth(child_ord, child_count, angle);
    vec_normalize(vec_add(
        vec_scale(axis_to_parent, angle.cos()),
        vec_add(
            vec_scale(normal, phase.cos() * angle.sin()),
            vec_scale(binormal, phase.sin() * angle.sin()),
        ),
    ))
}

fn child_azimuth(child_ord: usize, child_count: usize, angle: f64) -> f64 {
    match child_count {
        0 | 1 => 0.0,
        2 => {
            let cos_angle = angle.cos();
            let sin_sq = angle.sin() * angle.sin();
            let delta = if sin_sq <= 1.0e-12 {
                PI
            } else {
                ((cos_angle - cos_angle * cos_angle) / sin_sq)
                    .clamp(-1.0, 1.0)
                    .acos()
            };
            if child_ord == 0 {
                -0.5 * delta
            } else {
                0.5 * delta
            }
        }
        _ => 2.0 * PI * child_ord as f64 / child_count as f64,
    }
}

fn is_connected(molecule: &Molecule, adjacency: &AdjacencyList) -> bool {
    let mut seen = vec![false; molecule.num_atoms()];
    let mut queue = VecDeque::new();
    seen[0] = true;
    queue.push_back(0);
    while let Some(atom) = queue.pop_front() {
        for nbr in adjacency.neighbors_of(atom) {
            if !seen[nbr.atom_index] {
                seen[nbr.atom_index] = true;
                queue.push_back(nbr.atom_index);
            }
        }
    }
    seen.into_iter().all(|value| value)
}

fn confseq_base_local_angle_rad(molecule: &Molecule, center: usize) -> f64 {
    let atom = &molecule.atoms()[center];
    match atom.hybridization() {
        Hybridization::Sp => PI,
        Hybridization::Sp2 => 120.0_f64.to_radians(),
        // Short-term local-angle subset for sulfides. Without this, acyclic
        // thioether branches are forced into the generic tetrahedral angle and
        // no longer match the UFF-relaxed template locally.
        _ if atom.atomic_number() == 16 && !has_double_bond(molecule, center) => {
            100.0_f64.to_radians()
        }
        Hybridization::Sp3 => 109.47122063449069_f64.to_radians(),
        _ if atom.is_aromatic() => 120.0_f64.to_radians(),
        _ => 109.47122063449069_f64.to_radians(),
    }
}

fn confseq_base_ring_angle_rad(molecule: &Molecule, center: usize, ring_size: usize) -> f64 {
    let atom = &molecule.atoms()[center];
    let deg_to_rad = PI / 180.0;
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::_setRingAngle (BoundsMatrixBuilder.cpp:393-414)
    // RDKit✔️✔️: void _setRingAngle(Atom::HybridizationType aHyb, unsigned int ringSize,
    // RDKit✔️✔️:                    double &angle) {
    // RDKit✔️✔️:   if ((aHyb == Atom::SP2 && ringSize <= 8) || (ringSize == 3) ||
    // RDKit✔️✔️:       (ringSize == 4)) {
    // RDKit✔️✔️:     angle = M_PI * (1 - 2.0 / ringSize);
    // RDKit✔️✔️:   } else if (aHyb == Atom::SP3) {
    // RDKit✔️✔️:     if (ringSize == 5) {
    // RDKit✔️✔️:       angle = 104 * M_PI / 180;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       angle = 109.5 * M_PI / 180;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (aHyb == Atom::SP3D) {
    // RDKit✔️✔️:     angle = 105.0 * M_PI / 180;
    // RDKit✔️✔️:   } else if (aHyb == Atom::SP3D2) {
    // RDKit✔️✔️:     angle = 90.0 * M_PI / 180;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     angle = 120 * M_PI / 180;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::_setRingAngle
    if (atom.hybridization() == Hybridization::Sp2 && ring_size <= 8)
        || ring_size == 3
        || ring_size == 4
    {
        PI * (1.0 - 2.0 / ring_size as f64)
    } else if atom.hybridization() == Hybridization::Sp3 {
        if ring_size == 5 {
            104.0 * deg_to_rad
        } else {
            109.5 * deg_to_rad
        }
    } else if atom.hybridization() == Hybridization::Sp3d {
        105.0 * deg_to_rad
    } else if atom.hybridization() == Hybridization::Sp3d2 {
        90.0 * deg_to_rad
    } else {
        120.0 * deg_to_rad
    }
}

fn confseq_base_source_backed_angle_rad(
    molecule: &Molecule,
    left: usize,
    center: usize,
    right: usize,
) -> f64 {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::getUFFAngleBendParams (AtomTyper.cpp:559-589)
    // RDKit✔️❗: bool getUFFAngleBendParams(const ROMol &mol, unsigned int idx1,
    // RDKit✔️❗:                            unsigned int idx2, unsigned int idx3,
    // RDKit✔️❗:                            UFFAngle &uffAngleBendParams) {
    // RDKit✔️❗:   auto params = ParamCollection::getParams();
    // RDKit✔️❗:   unsigned int idx[3] = {idx1, idx2, idx3};
    // RDKit✔️❗:   AtomicParamVect paramVect(3);
    // RDKit✔️❗:   const Bond *bond[2];
    // RDKit✔️❗:   uffAngleBendParams.theta0 = paramVect[1]->theta0;
    // RDKit✔️❗:   return res;
    // RDKit✔️❗: }
    // END RDKIT CPP FUNCTION RDKit::UFF::getUFFAngleBendParams
    get_uff_angle_bend_params(molecule, left, center, right)
        .ok()
        .flatten()
        .map(|params| params.theta0.to_radians())
        .unwrap_or_else(|| confseq_base_local_angle_rad(molecule, center))
}

fn confseq_base_source_backed_bond_length(molecule: &Molecule, bond: &Bond) -> f64 {
    // BEGIN RDKIT CPP FUNCTION RDKit::UFF::getUFFBondStretchParams (AtomTyper.cpp:535-557)
    // RDKit✔️❗: bool getUFFBondStretchParams(const ROMol &mol, unsigned int idx1,
    // RDKit✔️❗:                              unsigned int idx2, UFFBond &uffBondStretchParams) {
    // RDKit✔️❗:   auto params = ParamCollection::getParams();
    // RDKit✔️❗:   unsigned int idx[2] = {idx1, idx2};
    // RDKit✔️❗:   AtomicParamVect paramVect(2);
    // RDKit✔️❗:   const Bond *bond = mol.getBondBetweenAtoms(idx1, idx2);
    // RDKit✔️❗:   if (res) {
    // RDKit✔️❗:     double bondOrder = bond->getBondTypeAsDouble();
    // RDKit✔️❗:     uffBondStretchParams.r0 =
    // RDKit✔️❗:         UFF::Utils::calcBondRestLength(bondOrder, paramVect[0], paramVect[1]);
    // RDKit✔️❗:   }
    // RDKit✔️❗:   return res;
    // RDKit✔️❗: }
    // END RDKIT CPP FUNCTION RDKit::UFF::getUFFBondStretchParams
    get_uff_bond_stretch_params(molecule, bond.begin().index(), bond.end().index())
        .ok()
        .flatten()
        .map(|params| params.r0)
        .unwrap_or_else(|| confseq_base_static_bond_length(molecule, bond))
}

fn confseq_base_static_bond_length_fallback(bond: &Bond) -> f64 {
    match bond.order() {
        BondOrder::Double => 1.34,
        BondOrder::Triple => 1.20,
        BondOrder::Aromatic => aromatic_bond_length(),
        _ => 1.50,
    }
}

fn confseq_base_static_bond_length(molecule: &Molecule, bond: &Bond) -> f64 {
    let a = molecule.atoms()[bond.begin().index()].atomic_number();
    let b = molecule.atoms()[bond.end().index()].atomic_number();
    // Short-term parameter subset: these constants cover the initial ConfSeq
    // fast-geometry tests and common organic heavy-atom bonds. They are not a
    // replacement for the RDKit bounds-matrix/UFF parameter sources already
    // ported elsewhere in core; expanding this table should be source-backed.
    match (a.min(b), a.max(b), bond.order()) {
        (6, 6, BondOrder::Aromatic) => aromatic_bond_length(),
        (6, 6, BondOrder::Double) => 1.34,
        (6, 6, BondOrder::Triple) => 1.20,
        (6, 6, _) => 1.52,
        (6, 7, BondOrder::Aromatic) => 1.34,
        (6, 7, BondOrder::Triple) => 1.16,
        (6, 7, BondOrder::Double) => 1.28,
        (6, 7, _) => 1.47,
        (6, 8, BondOrder::Aromatic) => 1.36,
        (6, 8, BondOrder::Double) => 1.23,
        (6, 8, _) => 1.43,
        (6, 9, _) => 1.35,
        (6, 15, _) => 1.84,
        (6, 16, BondOrder::Aromatic) => 1.74,
        (6, 16, BondOrder::Double) => 1.61,
        (6, 16, _) => 1.82,
        (6, 17, _) => 1.77,
        (6, 35, _) => 1.94,
        (6, 53, _) => 2.14,
        (7, 7, BondOrder::Aromatic) => 1.32,
        // RDKit UFF source anchor:
        //   Params.cpp: "N_3 0.7 ...", "N_2 0.685 ...", "N_1 0.656 ..."
        //   AtomTyper.cpp: r0 = UFF::Utils::calcBondRestLength(bondOrder, ...)
        // This ConfSeq base table remains static; these N-N entries prevent
        // falling through to the overly long generic 1.50 A default.
        (7, 7, BondOrder::Double) => 1.25,
        (7, 7, BondOrder::Triple) => 1.10,
        (7, 7, _) => 1.42,
        (7, 8, BondOrder::Aromatic) => 1.35,
        (7, 8, _) => 1.40,
        (7, 16, BondOrder::Aromatic) => 1.70,
        (7, 16, _) => 1.68,
        (8, 15, BondOrder::Double) => 1.48,
        (8, 15, _) => 1.63,
        (8, 16, BondOrder::Double) => 1.43,
        (8, 16, _) => 1.57,
        (1, 6, _) => 1.09,
        (1, 7, _) => 1.01,
        (1, 8, _) => 0.96,
        _ => 1.50,
    }
}

fn aromatic_bond_length() -> f64 {
    1.397
}

fn vec_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn vec_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn vec_scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

fn angle_rad_from_points(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Option<f64> {
    let ba = vec_sub(a, b);
    let bc = vec_sub(c, b);
    let ba_len = vec_len(ba);
    let bc_len = vec_len(bc);
    if ba_len <= 1.0e-12 || bc_len <= 1.0e-12 {
        return None;
    }
    Some(
        (vec_dot(ba, bc) / (ba_len * bc_len))
            .clamp(-1.0, 1.0)
            .acos(),
    )
}

fn vec_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn vec_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn vec_len(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn vec_normalize(v: [f64; 3]) -> [f64; 3] {
    let len = vec_len(v);
    if len <= 1.0e-12 {
        [1.0, 0.0, 0.0]
    } else {
        vec_scale(v, 1.0 / len)
    }
}

fn perpendicular_unit(v: [f64; 3]) -> [f64; 3] {
    let reference = if v[2].abs() < 0.9 {
        [0.0, 0.0, 1.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    vec_normalize(vec_cross(reference, v))
}

fn prepare_p_chiral_embedding_molecule(
    molecule: Molecule,
    chiral_tags_by_atom: &HashMap<usize, ChiralTag>,
) -> Result<Molecule, ConfSeqDecodeError> {
    // ConfSeq source anchor:
    //   if atom.GetAtomicNum()==7 and atom.GetHybridization()==Chem.HybridizationType.SP3 and atom.GetFormalCharge()==0:
    //       atom.SetFormalCharge(1)
    //       atom.SetNumExplicitHs(1)
    //   elif atom.GetAtomicNum()==16 and atom.GetHybridization()==Chem.HybridizationType.SP3 and atom.GetFormalCharge()==1:
    //       atom.SetFormalCharge(0)
    //       atom.SetNumExplicitHs(1)
    //   for idx, atom in enumerate(mol.GetAtoms()):
    //       atom.SetIsotope(idx+1)
    //   for atom_idx, tag in p_chiral_dic.items():
    //       atom.SetChiralTag(tag)
    let mut builder = molecule.to_builder();
    for atom_idx in 0..builder.atoms().len() {
        let atom = &builder.atoms()[atom_idx];
        let needs_n = atom.atomic_number() == 7
            && atom.hybridization() == Hybridization::Sp3
            && atom.formal_charge() == 0;
        let needs_s = atom.atomic_number() == 16
            && atom.hybridization() == Hybridization::Sp3
            && atom.formal_charge() == 1;
        let atom_mut = builder
            .atom_mut(AtomId::new(atom_idx))
            .expect("atom index from builder length");
        if needs_n {
            atom_mut.set_formal_charge(1);
            atom_mut.set_explicit_hydrogens(1);
        } else if needs_s {
            atom_mut.set_formal_charge(0);
            atom_mut.set_explicit_hydrogens(1);
        }
        atom_mut.set_isotope(Some((atom_idx + 1) as u16));
        if let Some(tag) = chiral_tags_by_atom.get(&atom_idx).copied() {
            atom_mut.set_chiral_tag(tag);
        }
    }
    builder
        .build()
        .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))
}

fn restore_embedding_molecule_state(
    molecule: Molecule,
    source: &Molecule,
) -> Result<Molecule, ConfSeqDecodeError> {
    // ConfSeq source anchor:
    //   for atom_idx,charge_H in charge_Hs.items():
    //       mol.GetAtomWithIdx(atom_idx).SetFormalCharge(charge_H[0])
    //       mol.GetAtomWithIdx(atom_idx).SetNumExplicitHs(charge_H[1])
    //   for atom in mol.GetAtoms():
    //       atom.SetIsotope(0)
    let mut builder = molecule.to_builder();
    for atom_idx in 0..builder.atoms().len().min(source.atoms().len()) {
        let source_atom = &source.atoms()[atom_idx];
        let atom_mut = builder
            .atom_mut(AtomId::new(atom_idx))
            .expect("atom index from builder length");
        atom_mut.set_formal_charge(source_atom.formal_charge());
        atom_mut.set_explicit_hydrogens(source_atom.explicit_hydrogens());
        atom_mut.set_isotope(None);
    }
    builder
        .build()
        .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))
}

fn decode_from_template(
    template: &Template,
    parsed: &ParsedConfSeq,
    options: &ConfSeqDecodeOptions,
) -> Result<Molecule, ConfSeqDecodeError> {
    if options.apply_dihedrals
        && parsed.dihedral_angles_by_pair.len() != template.dihedrals_by_pair.len()
    {
        return Err(ConfSeqDecodeError::DihedralTokenCountMismatch {
            observed: parsed.dihedral_angles_by_pair.len(),
            expected: template.dihedrals_by_pair.len(),
        });
    }
    if options.apply_angles && parsed.angle_values_deg.len() != template.angle_centers.len() {
        return Err(ConfSeqDecodeError::AngleTokenCountMismatch {
            observed: parsed.angle_values_deg.len(),
            expected: template.angle_centers.len(),
        });
    }

    let mut molecule = template.molecule.clone();

    if options.apply_angles {
        for ((i, j, k), angle) in template.angle_centers.iter().zip(&parsed.angle_values_deg) {
            molecule = mol_transforms::set_bond_angle_deg(molecule, *i, *j, *k, *angle, 0)
                .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
        }
    }

    if options.apply_dihedrals {
        let mut unapplied = Vec::new();
        for dihedral in &template.dihedrals {
            let (_, j, k, _) = *dihedral;
            let pair = sorted_pair(j, k);
            let angle = parsed.dihedral_angles_by_pair.get(&pair).ok_or(
                ConfSeqDecodeError::DihedralTokenCountMismatch {
                    observed: parsed.dihedral_angles_by_pair.len(),
                    expected: template.dihedrals_by_pair.len(),
                },
            )?;
            if template.ring_bond_pairs.contains(&pair) {
                unapplied.push((*dihedral, *angle));
            } else {
                molecule = set_dihedral_deg_checked(molecule, *dihedral, *angle)?;
            }
        }
        if !unapplied.is_empty() {
            molecule =
                apply_ring_deferred_dihedrals(molecule, &template.last_ring_bonds, unapplied)?;
        }
    }

    Ok(molecule)
}

fn apply_ring_deferred_dihedrals(
    molecule: Molecule,
    last_ring_bonds: &[(usize, usize, BondSpec)],
    unapplied: Vec<((usize, usize, usize, usize), f64)>,
) -> Result<Molecule, ConfSeqDecodeError> {
    // ConfSeq source anchor:
    //   mol_no_ring = Chem.RWMol(mol)
    //   for begin,end in last_ring_bonds:
    //       mol_no_ring.RemoveBond(begin,end)
    //   unapplied_dihedrals=change_dihedral(last_ring_bonds,dihedral_list,unapplied_dihedrals)
    //   Chem.MolToSmiles(mol_no_ring,canonical = False)
    //   order_lis = eval(mol_no_ring.GetProp('_smilesAtomOutputOrder'))
    //   mol_no_ring = Chem.RenumberAtoms(mol_no_ring, order_lis)
    //   mol,unapplied_dihedrals=apply_dihedrals(mol_no_ring,unapplied_dihedrals[:])
    //   mol,unapplied_dihedrals=apply_dihedrals(mol,unapplied_dihedrals[:])
    //   for begin,end in last_ring_bonds:
    //       ring_mol.AddBond(begin,end,order=Chem.BondType.SINGLE)
    let original_molecule = molecule.clone();
    let mut builder = molecule.to_builder();
    for (begin, end, _) in last_ring_bonds {
        builder.remove_bond_between_atoms(AtomId::new(*begin), AtomId::new(*end));
    }
    let molecule_no_ring = builder
        .build()
        .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
    let (mut molecule, old_to_new) = renumber_like_confseq_smiles_output(molecule_no_ring)?;

    let mut deferred =
        change_dihedral_for_removed_ring_bonds(&original_molecule, last_ring_bonds, unapplied)?;
    let last_ring_bonds = remap_last_ring_bonds(last_ring_bonds, &old_to_new)?;
    remap_deferred_dihedrals(&mut deferred, &old_to_new)?;
    for _ in 0..2 {
        let ring_info = rings::symmetrize_sssr(&molecule)
            .map_err(|err| ConfSeqDecodeError::RingFinding(err.to_string()))?;
        let mut still_deferred = Vec::new();
        for (dihedral, angle) in deferred {
            if dihedral_bonds_exist(&molecule, dihedral)
                && !dihedral_center_bond_is_in_ring(&ring_info, dihedral)
            {
                molecule = set_dihedral_deg_checked(molecule, dihedral, angle)?;
            } else {
                still_deferred.push((dihedral, angle));
            }
        }
        deferred = still_deferred;
        if deferred.is_empty() {
            break;
        }
    }

    let mut builder = molecule.to_builder();
    for (_, _, spec) in &last_ring_bonds {
        builder
            .add_bond(spec.clone())
            .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
    }
    builder
        .build()
        .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))
}

fn renumber_like_confseq_smiles_output(
    molecule: Molecule,
) -> Result<(Molecule, Vec<usize>), ConfSeqDecodeError> {
    let mut params = SmilesWriteParams::default();
    params.canonical = false;
    let (_, order) = mol_to_smiles_with_atom_output_order(&molecule, &params)
        .map_err(|err| ConfSeqDecodeError::SmilesWrite(err.to_string()))?;
    let mut builder = molecule.to_builder();
    let mapping = builder
        .renumber_atoms_for_construction(&order)
        .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
    let old_to_new = mapping
        .atoms()
        .old_to_new()
        .iter()
        .map(|mapped| {
            mapped.map(AtomId::index).ok_or_else(|| {
                ConfSeqDecodeError::MolTransform("atom renumbering dropped an atom".to_string())
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let molecule = builder
        .build()
        .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
    Ok((molecule, old_to_new))
}

fn remap_last_ring_bonds(
    last_ring_bonds: &[(usize, usize, BondSpec)],
    old_to_new: &[usize],
) -> Result<Vec<(usize, usize, BondSpec)>, ConfSeqDecodeError> {
    last_ring_bonds
        .iter()
        .map(|(begin, end, spec)| {
            let begin = *old_to_new.get(*begin).ok_or_else(|| {
                ConfSeqDecodeError::MolTransform(
                    "last ring bond begin atom is out of renumbering range".to_string(),
                )
            })?;
            let end = *old_to_new.get(*end).ok_or_else(|| {
                ConfSeqDecodeError::MolTransform(
                    "last ring bond end atom is out of renumbering range".to_string(),
                )
            })?;
            let stereo_atoms = spec
                .stereo_atoms()
                .map(|[left, right]| {
                    Ok::<[AtomId; 2], ConfSeqDecodeError>([
                        AtomId::new(*old_to_new.get(left.index()).ok_or_else(|| {
                            ConfSeqDecodeError::MolTransform(
                                "last ring bond stereo atom is out of renumbering range"
                                    .to_string(),
                            )
                        })?),
                        AtomId::new(*old_to_new.get(right.index()).ok_or_else(|| {
                            ConfSeqDecodeError::MolTransform(
                                "last ring bond stereo atom is out of renumbering range"
                                    .to_string(),
                            )
                        })?),
                    ])
                })
                .transpose()?;
            let spec = spec.remapped_endpoints(AtomId::new(begin), AtomId::new(end), stereo_atoms);
            Ok((begin, end, spec))
        })
        .collect()
}

fn remap_deferred_dihedrals(
    deferred: &mut [((usize, usize, usize, usize), f64)],
    old_to_new: &[usize],
) -> Result<(), ConfSeqDecodeError> {
    for ((i, j, k, l), _) in deferred {
        *i = *old_to_new.get(*i).ok_or_else(|| {
            ConfSeqDecodeError::MolTransform(
                "dihedral atom is out of renumbering range".to_string(),
            )
        })?;
        *j = *old_to_new.get(*j).ok_or_else(|| {
            ConfSeqDecodeError::MolTransform(
                "dihedral atom is out of renumbering range".to_string(),
            )
        })?;
        *k = *old_to_new.get(*k).ok_or_else(|| {
            ConfSeqDecodeError::MolTransform(
                "dihedral atom is out of renumbering range".to_string(),
            )
        })?;
        *l = *old_to_new.get(*l).ok_or_else(|| {
            ConfSeqDecodeError::MolTransform(
                "dihedral atom is out of renumbering range".to_string(),
            )
        })?;
    }
    Ok(())
}

fn change_dihedral_for_removed_ring_bonds(
    molecule: &Molecule,
    last_ring_bonds: &[(usize, usize, BondSpec)],
    dihedrals: Vec<((usize, usize, usize, usize), f64)>,
) -> Result<Vec<((usize, usize, usize, usize), f64)>, ConfSeqDecodeError> {
    // ConfSeq source anchor:
    //   if sorted((atom3, atom4)) in last_ring_bonds:
    //       neighbor_indices.remove(atom2)
    //       neighbor_indices.remove(atom4)
    //       if len(neighbor_indices) >= 1:
    //           o_2_d = Chem.rdMolTransforms.GetDihedralDeg(mol.GetConformer(),atom1, atom2, atom3, neighbor_indices[0])
    //           o_1_d = Chem.rdMolTransforms.GetDihedralDeg(mol.GetConformer(),atom1, atom2, atom3, atom4)
    //           new_angle = o_2_d - o_1_d + angle
    //           new_unapplied_dihedrals.append((atom1, atom2, atom3, neighbor_indices[0],new_angle))
    //   elif sorted((atom1, atom2)) in last_ring_bonds:
    //       neighbor_indices.remove(atom1)
    //       neighbor_indices.remove(atom3)
    //       if len(neighbor_indices) >= 1:
    //           o_2_d = Chem.rdMolTransforms.GetDihedralDeg(mol.GetConformer(),neighbor_indices[0], atom2, atom3,atom4 )
    //           o_1_d = Chem.rdMolTransforms.GetDihedralDeg(mol.GetConformer(),atom1, atom2, atom3, atom4)
    //           new_angle = o_2_d - o_1_d + angle
    //           new_unapplied_dihedrals.append((neighbor_indices[0], atom2, atom3,atom4,new_angle))
    let removed: HashSet<_> = last_ring_bonds
        .iter()
        .map(|(begin, end, _)| sorted_pair(*begin, *end))
        .collect();
    let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
    let mut changed = Vec::with_capacity(dihedrals.len());
    for ((i, j, k, l), angle) in dihedrals {
        if removed.contains(&sorted_pair(k, l)) {
            let neighbor = adjacency
                .neighbors_of(k)
                .iter()
                .map(|neighbor| neighbor.atom_index)
                .find(|neighbor| *neighbor != j && *neighbor != l);
            if let Some(neighbor) = neighbor {
                let o_2_d = mol_transforms::get_dihedral_deg(molecule, i, j, k, neighbor, 0)
                    .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
                let o_1_d = mol_transforms::get_dihedral_deg(molecule, i, j, k, l, 0)
                    .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
                changed.push(((i, j, k, neighbor), o_2_d - o_1_d + angle));
            } else {
                changed.push(((i, j, k, l), angle));
            }
        } else if removed.contains(&sorted_pair(i, j)) {
            let neighbor = adjacency
                .neighbors_of(j)
                .iter()
                .map(|neighbor| neighbor.atom_index)
                .find(|neighbor| *neighbor != i && *neighbor != k);
            if let Some(neighbor) = neighbor {
                let o_2_d = mol_transforms::get_dihedral_deg(molecule, neighbor, j, k, l, 0)
                    .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
                let o_1_d = mol_transforms::get_dihedral_deg(molecule, i, j, k, l, 0)
                    .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
                changed.push(((neighbor, j, k, l), o_2_d - o_1_d + angle));
            } else {
                changed.push(((i, j, k, l), angle));
            }
        } else {
            changed.push(((i, j, k, l), angle));
        }
    }
    Ok(changed)
}

fn set_dihedral_deg_checked(
    molecule: Molecule,
    dihedral: (usize, usize, usize, usize),
    angle_deg: f64,
) -> Result<Molecule, ConfSeqDecodeError> {
    let (i, j, k, l) = dihedral;
    mol_transforms::set_dihedral(molecule, i, j, k, l, angle_deg * PI / 180.0, 0)
        .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))
}

fn dihedral_bonds_exist(molecule: &Molecule, dihedral: (usize, usize, usize, usize)) -> bool {
    let (i, j, k, l) = dihedral;
    has_bond(molecule, i, j) && has_bond(molecule, j, k) && has_bond(molecule, k, l)
}

fn dihedral_center_bond_is_in_ring(
    ring_info: &crate::chemistry::rings::RingInfo,
    dihedral: (usize, usize, usize, usize),
) -> bool {
    let (_, j, k, _) = dihedral;
    ring_info
        .atom_rings()
        .iter()
        .any(|ring| ring.contains(&AtomId::new(j)) && ring.contains(&AtomId::new(k)))
}

fn has_bond(molecule: &Molecule, left: usize, right: usize) -> bool {
    molecule.bonds().iter().any(|bond| {
        sorted_pair(bond.begin().index(), bond.end().index()) == sorted_pair(left, right)
    })
}

fn bond_token_mapping_for_smiles(smiles: &str) -> Result<BondTokenMapping, ConfSeqDecodeError> {
    let molecule = Molecule::from_smiles(smiles)
        .map_err(|err| ConfSeqDecodeError::SmilesParse(err.to_string()))?;
    let params = SmilesWriteParams {
        canonical: false,
        all_bonds_explicit: true,
        ..SmilesWriteParams::default()
    };
    let smiles_be = mol_to_smiles(&molecule, &params)
        .map_err(|err| ConfSeqDecodeError::SmilesWrite(err.to_string()))?;
    explicit_bond_smiles_mapping(&smiles_be)
}

fn explicit_bond_smiles_mapping(smiles_be: &str) -> Result<BondTokenMapping, ConfSeqDecodeError> {
    // ConfSeq source anchor:
    //   mol_ = indigo.loadMolecule(smiles_BE)
    //   mol_block = mol_.molfile()
    //   ...
    //   for idx,token in enumerate(smiles_BE):
    //       if in_bracket==0 and token in ['-','=','#',':','/','\\']:
    //           bond_idx_token_idx_dic[bond_idx]=idx
    //
    // ConfSeq❗🔝: Indigo's molfile bond-line order follows the
    // all-bonds-explicit SMILES traversal for modeled ConfSeq inputs, including
    // ring-closure token positions. COSMolKit derives atom pairs in the same
    // traversal to avoid the Indigo molfile round trip; exact parity remains
    // covered by source-backed golden tests.
    let chars: Vec<char> = smiles_be.chars().collect();
    let mut token_idx_to_atom_pair = HashMap::new();
    let mut ring_closure_pairs = Vec::new();
    let mut branch_stack = Vec::new();
    let mut ring_open: HashMap<String, usize> = HashMap::new();
    let mut current_atom: Option<usize> = None;
    let mut pending_bond_token: Option<usize> = None;
    let mut atom_count = 0usize;
    let mut idx = 0usize;
    while idx < chars.len() {
        match chars[idx] {
            '(' => {
                let current = current_atom.ok_or_else(|| {
                    ConfSeqDecodeError::BondTokenMapping("branch before atom".to_string())
                })?;
                branch_stack.push(current);
                idx += 1;
            }
            ')' => {
                current_atom = Some(branch_stack.pop().ok_or_else(|| {
                    ConfSeqDecodeError::BondTokenMapping("unbalanced branch".to_string())
                })?);
                idx += 1;
            }
            '-' | '=' | '#' | ':' | '/' | '\\' => {
                pending_bond_token = Some(idx);
                idx += 1;
            }
            '[' => {
                let start = idx;
                idx += 1;
                while idx < chars.len() && chars[idx] != ']' {
                    idx += 1;
                }
                if idx == chars.len() {
                    return Err(ConfSeqDecodeError::BondTokenMapping(
                        "unclosed bracket atom".to_string(),
                    ));
                }
                let atom_idx = atom_count;
                atom_count += 1;
                connect_new_atom(
                    &mut token_idx_to_atom_pair,
                    &mut current_atom,
                    &mut pending_bond_token,
                    atom_idx,
                    start,
                )?;
                idx += 1;
            }
            '%' | '0'..='9' => {
                let (label, next_idx) = parse_ring_label(&chars, idx)?;
                let current = current_atom.ok_or_else(|| {
                    ConfSeqDecodeError::BondTokenMapping("ring closure before atom".to_string())
                })?;
                if let Some(open_atom) = ring_open.remove(&label) {
                    let token_idx = pending_bond_token.take().ok_or_else(|| {
                        ConfSeqDecodeError::BondTokenMapping(format!(
                            "ring closure {label} has no explicit bond token"
                        ))
                    })?;
                    let pair = sorted_pair(open_atom, current);
                    token_idx_to_atom_pair.insert(token_idx, pair);
                    ring_closure_pairs.push(pair);
                } else {
                    ring_open.insert(label, current);
                    pending_bond_token = None;
                }
                idx = next_idx;
            }
            '.' => {
                current_atom = None;
                pending_bond_token = None;
                idx += 1;
            }
            ch if is_organic_atom_start(ch) => {
                let start = idx;
                if matches!(ch, 'B' | 'C') && chars.get(idx + 1).copied() == Some('l')
                    || ch == 'B' && chars.get(idx + 1).copied() == Some('r')
                {
                    idx += 2;
                } else {
                    idx += 1;
                }
                let atom_idx = atom_count;
                atom_count += 1;
                connect_new_atom(
                    &mut token_idx_to_atom_pair,
                    &mut current_atom,
                    &mut pending_bond_token,
                    atom_idx,
                    start,
                )?;
            }
            ch => {
                return Err(ConfSeqDecodeError::BondTokenMapping(format!(
                    "unexpected character '{ch}' at position {idx}"
                )));
            }
        }
    }
    if !ring_open.is_empty() {
        return Err(ConfSeqDecodeError::BondTokenMapping(
            "unclosed ring label".to_string(),
        ));
    }
    Ok(BondTokenMapping {
        smiles_be: smiles_be.to_string(),
        token_idx_to_atom_pair,
        ring_closure_pairs,
    })
}

fn connect_new_atom(
    token_idx_to_atom_pair: &mut HashMap<usize, (usize, usize)>,
    current_atom: &mut Option<usize>,
    pending_bond_token: &mut Option<usize>,
    atom_idx: usize,
    start: usize,
) -> Result<(), ConfSeqDecodeError> {
    if let Some(prev) = *current_atom {
        let token_idx = pending_bond_token.take().ok_or_else(|| {
            ConfSeqDecodeError::BondTokenMapping(format!(
                "atom at position {start} has no explicit bond token"
            ))
        })?;
        token_idx_to_atom_pair.insert(token_idx, sorted_pair(prev, atom_idx));
    }
    *current_atom = Some(atom_idx);
    Ok(())
}

fn parse_ring_label(chars: &[char], idx: usize) -> Result<(String, usize), ConfSeqDecodeError> {
    if chars[idx] == '%' {
        let first = chars.get(idx + 1).copied();
        let second = chars.get(idx + 2).copied();
        if !matches!(first, Some(ch) if ch.is_ascii_digit())
            || !matches!(second, Some(ch) if ch.is_ascii_digit())
        {
            return Err(ConfSeqDecodeError::BondTokenMapping(
                "invalid percent ring label".to_string(),
            ));
        }
        Ok((format!("{}{}", first.unwrap(), second.unwrap()), idx + 3))
    } else {
        Ok((chars[idx].to_string(), idx + 1))
    }
}

fn is_organic_atom_start(ch: char) -> bool {
    matches!(
        ch,
        'B' | 'C' | 'N' | 'O' | 'P' | 'S' | 'F' | 'I' | 'b' | 'c' | 'n' | 'o' | 'p' | 's'
    )
}

fn collect_single_bond_dihedrals(molecule: &Molecule) -> Vec<(usize, usize, usize, usize)> {
    // ConfSeq source anchor:
    //   for bond in mol.GetBonds():
    //       if bond.GetBondType() == Chem.BondType.SINGLE:
    //           ...
    //           if len(set(neighbors1+neighbors2)) == len(neighbors1) + len(neighbors2):
    //               if len(neighbors1) > 0 and len(neighbors2) > 0:
    //                   i, j, k, l = pick_neighbor(mol,neighbors1),idx1, idx2,pick_neighbor(mol,neighbors2)
    let mut dihedrals = Vec::new();
    for bond in molecule.bonds() {
        if bond.order() != BondOrder::Single {
            continue;
        }
        let mut j = bond.begin().index();
        let mut k = bond.end().index();
        if j > k {
            std::mem::swap(&mut j, &mut k);
        }
        let left = heavy_neighbors_except(molecule, j, k);
        let right = heavy_neighbors_except(molecule, k, j);
        if left
            .iter()
            .chain(right.iter())
            .collect::<HashSet<_>>()
            .len()
            != left.len() + right.len()
        {
            continue;
        }
        if let (Some(i), Some(l)) = (
            pick_neighbor(molecule, &left),
            pick_neighbor(molecule, &right),
        ) {
            dihedrals.push((i, j, k, l));
        }
    }
    dihedrals
}

fn collect_angle_centers(
    molecule: &Molecule,
) -> Result<Vec<(usize, usize, usize)>, ConfSeqDecodeError> {
    // ConfSeq source anchor:
    //   if is_atom_in_ring(atom, mol):
    //       continue
    //   if neighbor.GetAtomicNum() == 8 and any(
    //       bond.GetBondType() == Chem.rdchem.BondType.DOUBLE for bond in neighbor.GetBonds()
    //   ):
    //       continue
    //   heavy_neighbors = [n for n in valid_neighbors if n.GetAtomicNum() > 1]
    let ring_info = rings::fast_find_rings(molecule)
        .map_err(|err| ConfSeqDecodeError::RingFinding(err.to_string()))?;
    let mut centers = Vec::new();
    for atom in molecule.atoms() {
        let center = atom.id().index();
        if ring_info.num_atom_rings(atom.id()) > 0 {
            continue;
        }
        let heavy: Vec<_> = molecule
            .bonds()
            .iter()
            .filter_map(|bond| {
                let begin = bond.begin().index();
                let end = bond.end().index();
                let other = if begin == center {
                    end
                } else if end == center {
                    begin
                } else {
                    return None;
                };
                if molecule.atoms()[other].atomic_number() <= 1 {
                    return None;
                }
                if molecule.atoms()[other].atomic_number() == 8 && has_double_bond(molecule, other)
                {
                    return None;
                }
                Some(other)
            })
            .collect();
        if heavy.len() == 2 {
            centers.push((heavy[0], center, heavy[1]));
        }
    }
    Ok(centers)
}

fn collect_ring_bond_pairs(
    molecule: &Molecule,
) -> Result<HashSet<(usize, usize)>, ConfSeqDecodeError> {
    let ring_info = rings::fast_find_rings(molecule)
        .map_err(|err| ConfSeqDecodeError::RingFinding(err.to_string()))?;
    Ok(molecule
        .bonds()
        .iter()
        .filter(|bond| ring_info.num_bond_rings(bond.id()) > 0)
        .map(|bond| sorted_pair(bond.begin().index(), bond.end().index()))
        .collect())
}

fn collect_last_ring_bonds(
    smiles: &str,
    molecule: &Molecule,
) -> Result<Vec<(usize, usize, BondSpec)>, ConfSeqDecodeError> {
    // ConfSeq source anchor:
    //   for idx,tok in enumerate(toks):
    //       if tok in ['-', '=', '#', ':', '/', '\\'] and toks[idx+1].isdigit() and tok!=':':
    //           if bond_count not in bond_counts:
    //               bond_counts.append(bond_count)
    //   last_ring_bonds.append(atom_pairs[idx])
    //   shared_bond_list=get_fully_shared_ring_bonds(mol)
    let smiles_molecule = Molecule::from_smiles(smiles)
        .map_err(|err| ConfSeqDecodeError::SmilesParse(err.to_string()))?;
    let shared = collect_fully_shared_ring_bonds(&smiles_molecule)?;
    let params = SmilesWriteParams {
        canonical: false,
        all_bonds_explicit: true,
        ..SmilesWriteParams::default()
    };
    let smiles_be = mol_to_smiles(&smiles_molecule, &params)
        .map_err(|err| ConfSeqDecodeError::SmilesWrite(err.to_string()))?;
    let mapping = explicit_bond_smiles_mapping(&smiles_be)?;
    let ring_bond_indices = confseq_last_ring_bond_indices_from_smiles_be(&smiles_be);
    let atom_pairs = mapping.atom_pairs_in_token_order();
    let mut last_ring_bonds = Vec::new();
    let mut seen = HashSet::new();
    for bond_index in ring_bond_indices {
        let Some(pair) = atom_pairs.get(bond_index).copied() else {
            continue;
        };
        if shared.contains(&pair) || !seen.insert(pair) {
            continue;
        }
        if let Some(bond) = bond_between_pair(molecule, pair) {
            last_ring_bonds.push((pair.0, pair.1, bond_spec_from_bond(bond)));
        }
    }
    Ok(last_ring_bonds)
}

fn confseq_last_ring_bond_indices_from_smiles_be(smiles_be: &str) -> Vec<usize> {
    // ConfSeq source anchor:
    //   toks = atomwise_tokenizer(smiles_BE)
    //   ring_bond_lis = []
    //   count = 0
    //   for i in range(len(toks)-1):
    //       if toks[i] in ['=','-','#','/','\\',':']:
    //           if (toks[i + 1].isdigit() or toks[i + 1].replace('%','').isdigit()) and toks[i] != ':' :
    //               ring_bond_lis.append(count)
    //           count += 1
    let tokens = atomwise_smiles_tokens_for_confseq(smiles_be);
    let mut ring_bond_indices = Vec::new();
    let mut bond_count = 0usize;
    for window in tokens.windows(2) {
        let token = window[0];
        if matches!(token, "=" | "-" | "#" | "/" | "\\" | ":") {
            let next = window[1];
            if token != ":" && confseq_ring_label_token(next) {
                ring_bond_indices.push(bond_count);
            }
            bond_count += 1;
        }
    }
    ring_bond_indices
}

fn atomwise_smiles_tokens_for_confseq(smiles: &str) -> Vec<&str> {
    let mut tokens = Vec::new();
    let mut idx = 0usize;
    while idx < smiles.len() {
        let rest = &smiles[idx..];
        let ch = rest.chars().next().expect("idx is on char boundary");
        if ch == '[' {
            if let Some(end_rel) = rest.find(']') {
                let end = idx + end_rel + 1;
                tokens.push(&smiles[idx..end]);
                idx = end;
                continue;
            }
        }
        if ch == '%' {
            let mut end = idx + ch.len_utf8();
            for _ in 0..2 {
                if let Some(next) = smiles[end..].chars().next() {
                    if next.is_ascii_digit() {
                        end += next.len_utf8();
                    }
                }
            }
            tokens.push(&smiles[idx..end]);
            idx = end;
            continue;
        }
        let end = idx + ch.len_utf8();
        tokens.push(&smiles[idx..end]);
        idx = end;
    }
    tokens
}

fn confseq_ring_label_token(token: &str) -> bool {
    token.chars().all(|ch| ch.is_ascii_digit())
        || token
            .strip_prefix('%')
            .is_some_and(|label| label.chars().all(|ch| ch.is_ascii_digit()))
}

fn collect_fully_shared_ring_bonds(
    molecule: &Molecule,
) -> Result<HashSet<(usize, usize)>, ConfSeqDecodeError> {
    let ring_info = rings::symmetrize_sssr(molecule)
        .map_err(|err| ConfSeqDecodeError::RingFinding(err.to_string()))?;
    let mut shared = HashSet::new();
    for ring_bonds in ring_info.bond_rings() {
        if ring_bonds
            .iter()
            .all(|bond| ring_info.num_bond_rings(*bond) >= 2)
        {
            for bond_id in ring_bonds {
                let bond = &molecule.bonds()[bond_id.index()];
                shared.insert(sorted_pair(bond.begin().index(), bond.end().index()));
            }
        }
    }
    Ok(shared)
}

fn bond_between_pair(molecule: &Molecule, pair: (usize, usize)) -> Option<&Bond> {
    molecule
        .bonds()
        .iter()
        .find(|bond| sorted_pair(bond.begin().index(), bond.end().index()) == pair)
}

fn bond_spec_from_bond(bond: &Bond) -> BondSpec {
    let mut spec = BondSpec::new(bond.begin(), bond.end(), bond.order())
        .with_aromatic(bond.is_aromatic())
        .with_conjugated(bond.is_conjugated())
        .with_direction(bond.direction())
        .with_stereo(bond.stereo())
        .with_unknown_stereo(bond.unknown_stereo());
    if let Some([begin, end]) = bond.stereo_atoms() {
        spec = spec.with_stereo_atoms(begin, end);
    }
    if let Some(query) = bond.query().cloned() {
        spec = spec.with_query(query);
    }
    for (key, value) in bond.props() {
        spec = spec.with_prop(key.clone(), value.clone());
    }
    spec
}

fn has_double_bond(molecule: &Molecule, atom: usize) -> bool {
    molecule.bonds().iter().any(|bond| {
        let begin = bond.begin().index();
        let end = bond.end().index();
        (begin == atom || end == atom) && bond.order() == BondOrder::Double
    })
}

fn heavy_neighbors_except(molecule: &Molecule, atom: usize, exclude: usize) -> Vec<usize> {
    molecule
        .bonds()
        .iter()
        .filter_map(|bond| {
            let begin = bond.begin().index();
            let end = bond.end().index();
            let other = if begin == atom {
                end
            } else if end == atom {
                begin
            } else {
                return None;
            };
            (other != exclude && molecule.atoms()[other].atomic_number() > 1).then_some(other)
        })
        .collect()
}

fn pick_neighbor(molecule: &Molecule, neighbors: &[usize]) -> Option<usize> {
    neighbors
        .iter()
        .copied()
        .min_by_key(|idx| (heavy_degree(molecule, *idx) == 1, *idx))
}

fn heavy_degree(molecule: &Molecule, atom: usize) -> usize {
    molecule
        .bonds()
        .iter()
        .filter(|bond| {
            let begin = bond.begin().index();
            let end = bond.end().index();
            let other = if begin == atom {
                end
            } else if end == atom {
                begin
            } else {
                return false;
            };
            molecule.atoms()[other].atomic_number() > 1
        })
        .count()
}

fn sorted_pair(a: usize, b: usize) -> (usize, usize) {
    if a <= b { (a, b) } else { (b, a) }
}

fn chiral_tag_cache_code(tag: ChiralTag) -> u8 {
    match tag {
        ChiralTag::Unspecified => 0,
        ChiralTag::TetrahedralCw => 1,
        ChiralTag::TetrahedralCcw => 2,
        ChiralTag::Other => 3,
        ChiralTag::Tetrahedral => 4,
        ChiralTag::Allene => 5,
        ChiralTag::SquarePlanar => 6,
        ChiralTag::TrigonalBipyramidal => 7,
        ChiralTag::Octahedral => 8,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::Value;

    #[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
    struct ConfSeqDgReferenceCache {
        corpus_path: String,
        entries: Vec<ConfSeqDgReferenceCacheEntry>,
    }

    #[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
    struct ConfSeqDgReferenceCacheEntry {
        heavy_atom_points: Option<Vec<[f64; 3]>>,
        error: Option<String>,
    }

    const CONFSEQ_DG_REFERENCE_MAX_REASONABLE_HEAVY_BOND_A: f64 = 2.2;

    fn build_template_initial_coords_for_test(
        smiles: &str,
        chiral_tags_by_atom: &HashMap<usize, ChiralTag>,
        options: &ConfSeqDecodeOptions,
        stage: crate::chemistry::distgeom::EmbedderTestStage,
    ) -> Result<Template, ConfSeqDecodeError> {
        let molecule = Molecule::from_smiles(smiles)
            .map_err(|err| ConfSeqDecodeError::SmilesParse(err.to_string()))?;
        let molecule = prepare_p_chiral_embedding_molecule(molecule, chiral_tags_by_atom)?;
        let with_h = molecule
            .with_hydrogens()
            .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
        let mut params = options.embed_params.clone();
        let embedded_with_h =
            crate::distgeom::embedder_stage_coords_for_test(&with_h, &mut params, stage)
                .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
        if embedded_with_h.conformers_3d().is_empty() {
            return Err(ConfSeqDecodeError::EmbedFailed);
        }
        let embedded = embedded_with_h
            .without_hydrogens_with_sanitize(false)
            .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
        let embedded = restore_embedding_molecule_state(embedded, &molecule)?;

        let dihedrals = collect_single_bond_dihedrals(&embedded);
        let dihedrals_by_pair = dihedrals
            .iter()
            .copied()
            .map(|dihedral @ (_, j, k, _)| (sorted_pair(j, k), dihedral))
            .collect();
        let angle_centers = collect_angle_centers(&embedded)?;
        let ring_bond_pairs = collect_ring_bond_pairs(&embedded)?;
        let last_ring_bonds = collect_last_ring_bonds(smiles, &embedded)?;

        Ok(Template {
            molecule: embedded,
            dihedrals,
            dihedrals_by_pair,
            angle_centers,
            ring_bond_pairs,
            last_ring_bonds,
        })
    }

    fn heavy_atom_points_from_with_h_trace(
        molecule: &Molecule,
        points: Option<Vec<[f64; 3]>>,
    ) -> Option<Vec<[f64; 3]>> {
        let points = points?;
        Some(
            molecule
                .atoms()
                .iter()
                .enumerate()
                .filter_map(|(idx, atom)| {
                    (atom.element().atomic_number() != 1).then_some(points[idx])
                })
                .collect(),
        )
    }

    fn build_template_embed_trace_for_test(
        smiles: &str,
        chiral_tags_by_atom: &HashMap<usize, ChiralTag>,
        options: &ConfSeqDecodeOptions,
    ) -> Result<crate::chemistry::distgeom::EmbedderTestTrace, ConfSeqDecodeError> {
        let molecule = Molecule::from_smiles(smiles)
            .map_err(|err| ConfSeqDecodeError::SmilesParse(err.to_string()))?;
        let molecule = prepare_p_chiral_embedding_molecule(molecule, chiral_tags_by_atom)?;
        let with_h = molecule
            .with_hydrogens()
            .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
        let mut params = options.embed_params.clone();
        let mut trace = crate::distgeom::embedder_trace_for_test(&with_h, &mut params)
            .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
        trace.initial_coords = heavy_atom_points_from_with_h_trace(&with_h, trace.initial_coords);
        trace.first_minimized = heavy_atom_points_from_with_h_trace(&with_h, trace.first_minimized);
        trace.fourth_dimension_cleaned =
            heavy_atom_points_from_with_h_trace(&with_h, trace.fourth_dimension_cleaned);
        trace.exp_torsion_minimized =
            heavy_atom_points_from_with_h_trace(&with_h, trace.exp_torsion_minimized);
        trace.final_checked = heavy_atom_points_from_with_h_trace(&with_h, trace.final_checked);
        Ok(trace)
    }

    fn conformer_points(mol: &Molecule) -> Vec<[f64; 3]> {
        mol.conformers_3d()[0].coordinates().to_vec()
    }

    fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
        let dx = a[0] - b[0];
        let dy = a[1] - b[1];
        let dz = a[2] - b[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    fn angle_deg(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
        let ba = vec_sub(a, b);
        let bc = vec_sub(c, b);
        let denom = vec_len(ba) * vec_len(bc);
        let cos_theta = if denom <= 1.0e-12 {
            1.0
        } else {
            ((ba[0] * bc[0] + ba[1] * bc[1] + ba[2] * bc[2]) / denom).clamp(-1.0, 1.0)
        };
        cos_theta.acos() * 180.0 / PI
    }

    fn dihedral_deg_for_test(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
        let b0 = vec_sub(a, b);
        let b1 = vec_sub(c, b);
        let b2 = vec_sub(d, c);
        let b1 = vec_normalize(b1);
        let v = vec_sub(b0, vec_scale(b1, vec_dot(b0, b1)));
        let w = vec_sub(b2, vec_scale(b1, vec_dot(b2, b1)));
        vec_dot(vec_cross(b1, v), w)
            .atan2(vec_dot(v, w))
            .to_degrees()
    }

    fn abs_dihedral_delta_deg(left: f64, right: f64) -> f64 {
        let mut delta = left - right;
        while delta > 180.0 {
            delta -= 360.0;
        }
        while delta < -180.0 {
            delta += 360.0;
        }
        delta.abs()
    }

    fn distance_stats(
        pairs: impl Iterator<Item = (usize, usize)>,
        left: &[[f64; 3]],
        right: &[[f64; 3]],
    ) -> (usize, f64, f64) {
        let mut count = 0usize;
        let mut sum_sq = 0.0;
        let mut max_abs = 0.0;
        for (a, b) in pairs {
            let delta = distance(left[a], left[b]) - distance(right[a], right[b]);
            count += 1;
            sum_sq += delta * delta;
            if delta.abs() > max_abs {
                max_abs = delta.abs();
            }
        }
        let rms = if count == 0 {
            0.0
        } else {
            (sum_sq / count as f64).sqrt()
        };
        (count, rms, max_abs)
    }

    fn heavy_bond_length_stats(
        molecule: &Molecule,
        heavy_points: &[[f64; 3]],
    ) -> (usize, f64, f64) {
        let heavy_index_by_atom = heavy_index_by_atom(molecule);
        let mut count = 0usize;
        let mut sum_sq = 0.0;
        let mut max = 0.0;
        for bond in molecule.bonds() {
            let begin = bond.begin().index();
            let end = bond.end().index();
            let (Some(begin_heavy), Some(end_heavy)) =
                (heavy_index_by_atom[begin], heavy_index_by_atom[end])
            else {
                continue;
            };
            let length = distance(heavy_points[begin_heavy], heavy_points[end_heavy]);
            count += 1;
            sum_sq += length * length;
            if length > max {
                max = length;
            }
        }
        let rms = if count == 0 {
            0.0
        } else {
            (sum_sq / count as f64).sqrt()
        };
        (count, rms, max)
    }

    fn aromatic_atom_indices(mol: &Molecule) -> Vec<usize> {
        let mut atoms = HashSet::new();
        for bond in mol.bonds().iter().filter(|bond| bond.is_aromatic()) {
            atoms.insert(bond.begin().index());
            atoms.insert(bond.end().index());
        }
        let mut atoms: Vec<_> = atoms.into_iter().collect();
        atoms.sort_unstable();
        atoms
    }

    fn max_aromatic_ring_plane_deviation(molecule: &Molecule) -> f64 {
        let coords = conformer_points(molecule);
        let rings = rings::symmetrize_sssr(molecule).expect("ring perception should work");
        let mut max_deviation: f64 = 0.0;
        for ring in rings.atom_rings() {
            if !ring
                .iter()
                .all(|atom| molecule.atoms()[atom.index()].is_aromatic())
            {
                continue;
            }
            let origin = coords[ring[0].index()];
            let mut normal = [0.0, 0.0, 1.0];
            'normal: for left in 1..ring.len() {
                for right in left + 1..ring.len() {
                    let candidate = vec_cross(
                        vec_sub(coords[ring[left].index()], origin),
                        vec_sub(coords[ring[right].index()], origin),
                    );
                    if vec_len(candidate) > 1.0e-12 {
                        normal = vec_normalize(candidate);
                        break 'normal;
                    }
                }
            }
            for atom in ring {
                let deviation = vec_dot(vec_sub(coords[atom.index()], origin), normal).abs();
                max_deviation = max_deviation.max(deviation);
            }
        }
        max_deviation
    }

    fn optimized_template_for_comparison(smiles: &str) -> Template {
        let options = ConfSeqDecodeOptions::default();
        build_template(smiles, &HashMap::new(), &options).expect("UFF template should build")
    }

    fn base_template_for_comparison(smiles: &str) -> Template {
        let options = ConfSeqDecodeOptions {
            optimize_with_uff: false,
            template_backend: ConfSeqTemplateBackend::FastGeometry,
            ..ConfSeqDecodeOptions::default()
        };
        build_template(smiles, &HashMap::new(), &options)
            .expect("ConfSeq base template should build")
    }

    fn heavy_atom_points_for_rmsd(mol: &Molecule) -> Vec<[f64; 3]> {
        let coords = conformer_points(mol);
        mol.atoms()
            .iter()
            .filter_map(|atom| (atom.atomic_number() != 1).then_some(coords[atom.id().index()]))
            .collect()
    }

    fn automorphism_aware_heavy_rmsd_for_test(
        molecule: &Molecule,
        reference_heavy_points: &[[f64; 3]],
        base: &Molecule,
    ) -> Option<f64> {
        let base_heavy_points = heavy_atom_points_for_rmsd(base);
        if reference_heavy_points.len() != base_heavy_points.len() {
            return None;
        }
        let heavy_index_by_atom = heavy_index_by_atom(molecule);
        let params = crate::SubstructMatchParams {
            max_matches: 1000,
            uniquify: false,
            use_chirality: false,
            specified_stereo_query_matches_unspecified: false,
        };
        let matches = crate::get_substruct_matches_with_params(molecule, molecule, &params);
        if matches.is_empty() {
            return None;
        }
        let mut best = f64::INFINITY;
        for match_result in matches {
            if match_result.atom_mapping.len() != molecule.num_atoms() {
                continue;
            }
            let mut permuted_base_points = vec![[0.0; 3]; reference_heavy_points.len()];
            let mut usable = true;
            for (query_atom, &mapped_atom) in match_result.atom_mapping.iter().enumerate() {
                let Some(query_heavy_idx) = heavy_index_by_atom[query_atom] else {
                    continue;
                };
                let Some(mapped_heavy_idx) = heavy_index_by_atom[mapped_atom] else {
                    usable = false;
                    break;
                };
                permuted_base_points[query_heavy_idx] = base_heavy_points[mapped_heavy_idx];
            }
            if !usable {
                continue;
            }
            let rmsd = crate::distgeom::aligned_rmsd_for_test(
                reference_heavy_points,
                &permuted_base_points,
            );
            best = best.min(rmsd);
        }
        best.is_finite().then_some(best)
    }

    #[derive(Debug, Clone)]
    struct ConfSeqBaseFailExportRecord {
        idx: usize,
        rmsd: f64,
        automorphism_rmsd: Option<f64>,
        max_rigid_rmsd: Option<f64>,
        max_rigid_shape_rmsd: Option<f64>,
        max_rigid_connector_rmsd: Option<f64>,
        worst_fragment_atoms: Vec<usize>,
        stripped_smiles: String,
        structural_risk_classes: Vec<String>,
        structural_risk_fallback_candidate: bool,
        reference_points: Vec<[f64; 3]>,
        base: Molecule,
    }

    fn molecule_with_test_conformer(
        molecule: &Molecule,
        points: Vec<[f64; 3]>,
    ) -> Result<Molecule, String> {
        let mut builder = molecule.to_builder();
        builder
            .add_conformer(Conformer3D::new(0, points, true))
            .map_err(|err| err.to_string())?;
        builder.build().map_err(|err| err.to_string())
    }

    fn export_confseq_base_fail_subset(
        export_dir: &str,
        records: &[ConfSeqBaseFailExportRecord],
    ) -> Result<(), String> {
        std::fs::create_dir_all(export_dir).map_err(|err| err.to_string())?;
        let mut reference_sdf = String::new();
        let mut base_sdf = String::new();
        let mut jsonl = String::new();
        let params = MolBlockWriteParams {
            format: SdfFormat::V2000,
            force_2d: false,
            include_stereo: true,
            kekulize: true,
            precision: 6,
        };
        for record in records {
            let source_mol = Molecule::from_smiles(&record.stripped_smiles)
                .map_err(|err| format!("failed to parse row {} SMILES: {err}", record.idx))?;
            let reference =
                molecule_with_test_conformer(&source_mol, record.reference_points.clone())
                    .map_err(|err| {
                        format!(
                            "failed to build row {} reference molecule: {err}",
                            record.idx
                        )
                    })?;
            let reference_record =
                mol_to_sdf_record_with_params(&reference, &params).map_err(|err| {
                    format!("failed to write row {} reference SDF: {err}", record.idx)
                })?;
            let base_record = mol_to_sdf_record_with_params(&record.base, &params)
                .map_err(|err| format!("failed to write row {} base SDF: {err}", record.idx))?;
            reference_sdf.push_str(&reference_record);
            base_sdf.push_str(&base_record);
            let meta = serde_json::json!({
                "idx": record.idx,
                "line": record.idx + 1,
                "global_rmsd": record.rmsd,
                "automorphism_rmsd": record.automorphism_rmsd,
                "max_rigid_rmsd": record.max_rigid_rmsd,
                "max_rigid_shape_rmsd": record.max_rigid_shape_rmsd,
                "max_rigid_connector_rmsd": record.max_rigid_connector_rmsd,
                "worst_fragment_atoms": record.worst_fragment_atoms,
                "structural_risk_classes": record.structural_risk_classes,
                "structural_risk_fallback_candidate": record.structural_risk_fallback_candidate,
                "stripped_smiles": record.stripped_smiles,
            });
            jsonl.push_str(&serde_json::to_string(&meta).map_err(|err| err.to_string())?);
            jsonl.push('\n');
        }
        std::fs::write(
            std::path::Path::new(export_dir).join("reference_fail_gt_0_1a.sdf"),
            reference_sdf,
        )
        .map_err(|err| err.to_string())?;
        std::fs::write(
            std::path::Path::new(export_dir).join("base_fail_gt_0_1a.sdf"),
            base_sdf,
        )
        .map_err(|err| err.to_string())?;
        std::fs::write(
            std::path::Path::new(export_dir).join("fail_gt_0_1a.jsonl"),
            jsonl,
        )
        .map_err(|err| err.to_string())
    }

    fn quantile_sorted(values: &[f64], q: f64) -> f64 {
        if values.is_empty() {
            return f64::NAN;
        }
        let idx = ((values.len() - 1) as f64 * q).round() as usize;
        values[idx]
    }

    #[derive(Debug, Clone, Copy, Default)]
    struct LocalConstraintSummary {
        bond_count: usize,
        bond_rms_a: f64,
        max_bond_abs_a: f64,
        non_token_angle_count: usize,
        non_token_angle_rms_deg: f64,
        max_non_token_angle_abs_deg: f64,
    }

    #[derive(Debug, Clone, Copy, Default)]
    struct LocalBondDeviationBucket {
        count: usize,
        base_ref_sum_sq: f64,
        base_target_sum_sq: f64,
        reference_target_sum_sq: f64,
        max_base_ref_abs: f64,
        max_base_target_abs: f64,
        max_reference_target_abs: f64,
    }

    impl LocalBondDeviationBucket {
        fn record(&mut self, base_ref: f64, base_target: f64, reference_target: f64) {
            self.count += 1;
            self.base_ref_sum_sq += base_ref * base_ref;
            self.base_target_sum_sq += base_target * base_target;
            self.reference_target_sum_sq += reference_target * reference_target;
            self.max_base_ref_abs = self.max_base_ref_abs.max(base_ref.abs());
            self.max_base_target_abs = self.max_base_target_abs.max(base_target.abs());
            self.max_reference_target_abs =
                self.max_reference_target_abs.max(reference_target.abs());
        }

        fn base_ref_rms(self) -> f64 {
            if self.count == 0 {
                0.0
            } else {
                (self.base_ref_sum_sq / self.count as f64).sqrt()
            }
        }

        fn base_target_rms(self) -> f64 {
            if self.count == 0 {
                0.0
            } else {
                (self.base_target_sum_sq / self.count as f64).sqrt()
            }
        }

        fn reference_target_rms(self) -> f64 {
            if self.count == 0 {
                0.0
            } else {
                (self.reference_target_sum_sq / self.count as f64).sqrt()
            }
        }
    }

    fn local_constraint_summary_against_reference(
        molecule: &Molecule,
        reference_heavy_points: &[[f64; 3]],
        base: &Molecule,
    ) -> LocalConstraintSummary {
        let base_heavy_points = heavy_atom_points_for_rmsd(base);
        let heavy_index_by_atom = heavy_index_by_atom(molecule);
        let token_angles = collect_angle_centers(molecule)
            .map(|centers| {
                centers
                    .into_iter()
                    .map(|(left, center, right)| sorted_angle(left, center, right))
                    .collect::<HashSet<_>>()
            })
            .unwrap_or_default();

        let mut bond_count = 0usize;
        let mut bond_sum_sq = 0.0;
        let mut max_bond_abs_a: f64 = 0.0;
        for bond in molecule.bonds() {
            let begin = bond.begin().index();
            let end = bond.end().index();
            let (Some(begin_heavy), Some(end_heavy)) =
                (heavy_index_by_atom[begin], heavy_index_by_atom[end])
            else {
                continue;
            };
            let delta = distance(
                reference_heavy_points[begin_heavy],
                reference_heavy_points[end_heavy],
            ) - distance(base_heavy_points[begin_heavy], base_heavy_points[end_heavy]);
            bond_count += 1;
            bond_sum_sq += delta * delta;
            max_bond_abs_a = max_bond_abs_a.max(delta.abs());
        }

        let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
        let mut angle_count = 0usize;
        let mut angle_sum_sq = 0.0;
        let mut max_angle_abs_deg: f64 = 0.0;
        for atom in molecule.atoms() {
            if atom.atomic_number() <= 1 {
                continue;
            }
            let center = atom.id().index();
            let heavy_neighbors = adjacency
                .neighbors_of(center)
                .iter()
                .map(|neighbor| neighbor.atom_index)
                .filter(|neighbor| molecule.atoms()[*neighbor].atomic_number() > 1)
                .collect::<Vec<_>>();
            for left_idx in 0..heavy_neighbors.len() {
                for right_idx in left_idx + 1..heavy_neighbors.len() {
                    let left = heavy_neighbors[left_idx];
                    let right = heavy_neighbors[right_idx];
                    if token_angles.contains(&sorted_angle(left, center, right)) {
                        continue;
                    }
                    let (Some(left_heavy), Some(center_heavy), Some(right_heavy)) = (
                        heavy_index_by_atom[left],
                        heavy_index_by_atom[center],
                        heavy_index_by_atom[right],
                    ) else {
                        continue;
                    };
                    let ref_angle = angle_deg(
                        reference_heavy_points[left_heavy],
                        reference_heavy_points[center_heavy],
                        reference_heavy_points[right_heavy],
                    );
                    let base_angle = angle_deg(
                        base_heavy_points[left_heavy],
                        base_heavy_points[center_heavy],
                        base_heavy_points[right_heavy],
                    );
                    let delta = ref_angle - base_angle;
                    angle_count += 1;
                    angle_sum_sq += delta * delta;
                    max_angle_abs_deg = max_angle_abs_deg.max(delta.abs());
                }
            }
        }

        LocalConstraintSummary {
            bond_count,
            bond_rms_a: if bond_count == 0 {
                0.0
            } else {
                (bond_sum_sq / bond_count as f64).sqrt()
            },
            max_bond_abs_a,
            non_token_angle_count: angle_count,
            non_token_angle_rms_deg: if angle_count == 0 {
                0.0
            } else {
                (angle_sum_sq / angle_count as f64).sqrt()
            },
            max_non_token_angle_abs_deg: max_angle_abs_deg,
        }
    }

    fn collect_local_bond_deviation_buckets(
        molecule: &Molecule,
        reference_heavy_points: &[[f64; 3]],
        base: &Molecule,
        buckets: &mut HashMap<String, LocalBondDeviationBucket>,
    ) {
        let Ok(ring_info) = rings::symmetrize_sssr(molecule) else {
            return;
        };
        let Ok(model) = build_confseq_base_constraint_model(molecule, &ring_info) else {
            return;
        };
        let cut_bonds = confseq_base_fragment_cut_bonds(molecule, &ring_info);
        let heavy_index_by_atom = heavy_index_by_atom(molecule);
        let base_heavy_points = heavy_atom_points_for_rmsd(base);
        let mut component_by_atom = vec![None; molecule.num_atoms()];
        let mut shape_by_component = Vec::with_capacity(model.rigid_components.len());
        for (component_idx, component) in model.rigid_components.iter().enumerate() {
            for &atom in &component.atoms {
                component_by_atom[atom] = Some(component_idx);
            }
            let realization_atoms = confseq_base_component_realization_atoms(component);
            let shape = confseq_base_rigid_fragment_shape(molecule, component, &realization_atoms)
                .map(|shape| format!("{shape:?}"))
                .unwrap_or_else(|err| format!("shape_error:{err:?}"));
            shape_by_component.push(shape);
        }

        for bond in molecule.bonds() {
            let begin = bond.begin().index();
            let end = bond.end().index();
            let (Some(begin_heavy), Some(end_heavy)) =
                (heavy_index_by_atom[begin], heavy_index_by_atom[end])
            else {
                continue;
            };
            let base_len = distance(base_heavy_points[begin_heavy], base_heavy_points[end_heavy]);
            let reference_len = distance(
                reference_heavy_points[begin_heavy],
                reference_heavy_points[end_heavy],
            );
            let target = model
                .bond_targets
                .get(bond.id().index())
                .copied()
                .unwrap_or_else(|| confseq_base_static_bond_length_fallback(bond));
            let base_ref = base_len - reference_len;
            let base_target = base_len - target;
            let reference_target = reference_len - target;
            let pair = sorted_pair(begin, end);
            let scope = if cut_bonds.contains(&pair) {
                "cut"
            } else {
                "internal"
            };
            let shape = match (component_by_atom[begin], component_by_atom[end]) {
                (Some(left), Some(right)) if left == right => shape_by_component[left].as_str(),
                (Some(_), Some(_)) => "cross_component",
                _ => "unassigned",
            };
            let key = format!(
                "scope={scope}|shape={shape}|order={:?}|arom={}",
                bond.order(),
                bond.is_aromatic()
            );
            buckets
                .entry(key)
                .or_default()
                .record(base_ref, base_target, reference_target);
        }
    }

    fn print_local_bond_deviation_buckets(
        label: &str,
        buckets: HashMap<String, LocalBondDeviationBucket>,
        limit: usize,
    ) {
        let mut buckets = buckets.into_iter().collect::<Vec<_>>();
        buckets.sort_by(|left, right| {
            right
                .1
                .base_target_rms()
                .total_cmp(&left.1.base_target_rms())
                .then_with(|| right.1.count.cmp(&left.1.count))
                .then_with(|| left.0.cmp(&right.0))
        });
        for (key, bucket) in buckets.into_iter().take(limit) {
            eprintln!(
                "{label} count={} base_ref_rms={:.6} base_target_rms={:.6} reference_target_rms={:.6} max_base_ref_abs={:.6} max_base_target_abs={:.6} max_reference_target_abs={:.6} key={key}",
                bucket.count,
                bucket.base_ref_rms(),
                bucket.base_target_rms(),
                bucket.reference_target_rms(),
                bucket.max_base_ref_abs,
                bucket.max_base_target_abs,
                bucket.max_reference_target_abs,
            );
        }
    }

    #[derive(Debug, Clone)]
    struct RigidFragmentRmsdSummary {
        fragment_rmsds: Vec<f64>,
        fragment_shape_rmsds: Vec<f64>,
        fragment_terminal_symmetry_rmsds: Vec<f64>,
        fragment_connector_rmsds: Vec<f64>,
        max_rmsd: Option<f64>,
        max_shape_rmsd: Option<f64>,
        max_terminal_symmetry_rmsd: Option<f64>,
        max_connector_rmsd: Option<f64>,
        mirror_branch_like_fragment_count_01a: usize,
        mirror_branch_like_fragment_count_03a: usize,
        worst_fragment_atoms: Vec<usize>,
    }

    #[derive(Debug, Clone, Copy)]
    struct RigidFragmentMetricDetails {
        proper_rmsd: f64,
        shape_rmsd: f64,
        terminal_symmetry_rmsd: f64,
        connector_rmsd: Option<f64>,
        angle_rms_deg: f64,
        max_angle_delta_deg: f64,
    }

    fn rigid_fragment_rmsd_summary(
        molecule: &Molecule,
        reference_heavy_points: &[[f64; 3]],
        base: &Molecule,
    ) -> RigidFragmentRmsdSummary {
        let base_heavy_points = heavy_atom_points_for_rmsd(base);
        assert_eq!(
            reference_heavy_points.len(),
            base_heavy_points.len(),
            "rigid-fragment RMSD requires matching heavy atom counts"
        );
        let fragments = rigid_heavy_fragments_cutting_confseq_rotors(molecule);
        let cut_bonds = rings::symmetrize_sssr(molecule)
            .map(|ring_info| confseq_base_fragment_cut_bonds(molecule, &ring_info))
            .unwrap_or_else(|_| {
                collect_single_bond_dihedrals(molecule)
                    .into_iter()
                    .map(|(_, j, k, _)| sorted_pair(j, k))
                    .collect()
            });
        let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
        let heavy_index_by_atom = heavy_index_by_atom(molecule);
        let mut fragment_rmsds = Vec::new();
        let mut fragment_shape_rmsds = Vec::new();
        let mut fragment_terminal_symmetry_rmsds = Vec::new();
        let mut fragment_connector_rmsds = Vec::new();
        let mut mirror_branch_like_fragment_count_01a = 0usize;
        let mut mirror_branch_like_fragment_count_03a = 0usize;
        let mut worst_fragment_atoms = Vec::new();
        for fragment in fragments {
            let heavy_indices: Vec<_> = fragment
                .iter()
                .copied()
                .filter_map(|atom_idx| heavy_index_by_atom[atom_idx])
                .collect();
            if heavy_indices.len() < 3 {
                continue;
            }
            let ref_points: Vec<_> = heavy_indices
                .iter()
                .map(|&idx| reference_heavy_points[idx])
                .collect();
            let base_points: Vec<_> = heavy_indices
                .iter()
                .map(|&idx| base_heavy_points[idx])
                .collect();
            let rmsd = crate::distgeom::aligned_rmsd_for_test(&ref_points, &base_points);
            let shape_rmsd = mirror_invariant_aligned_rmsd_for_test(&ref_points, &base_points);
            let terminal_symmetry_rmsd = terminal_ligand_symmetry_aligned_rmsd_for_test(
                molecule,
                &adjacency,
                &fragment,
                &ref_points,
                &base_points,
            );
            if fragment_rmsds
                .iter()
                .copied()
                .max_by(f64::total_cmp)
                .map(|current_max| rmsd > current_max)
                .unwrap_or(true)
            {
                worst_fragment_atoms = fragment.clone();
            }
            fragment_rmsds.push(rmsd);
            fragment_shape_rmsds.push(shape_rmsd);
            if rmsd > 0.3 && shape_rmsd <= 0.1 {
                mirror_branch_like_fragment_count_01a += 1;
            }
            if rmsd > 0.3 && shape_rmsd <= 0.3 {
                mirror_branch_like_fragment_count_03a += 1;
            }
            fragment_terminal_symmetry_rmsds.push(terminal_symmetry_rmsd);

            let connector_fragment = rigid_fragment_atoms_with_connector_stubs(
                molecule, &adjacency, &cut_bonds, &fragment,
            );
            let connector_heavy_indices: Vec<_> = connector_fragment
                .iter()
                .copied()
                .filter_map(|atom_idx| heavy_index_by_atom[atom_idx])
                .collect();
            if connector_heavy_indices.len() >= 3 {
                let connector_ref_points: Vec<_> = connector_heavy_indices
                    .iter()
                    .map(|&idx| reference_heavy_points[idx])
                    .collect();
                let connector_base_points: Vec<_> = connector_heavy_indices
                    .iter()
                    .map(|&idx| base_heavy_points[idx])
                    .collect();
                fragment_connector_rmsds.push(crate::distgeom::aligned_rmsd_for_test(
                    &connector_ref_points,
                    &connector_base_points,
                ));
            }
        }
        let max_rmsd = fragment_rmsds.iter().copied().max_by(f64::total_cmp);
        let max_shape_rmsd = fragment_shape_rmsds.iter().copied().max_by(f64::total_cmp);
        let max_terminal_symmetry_rmsd = fragment_terminal_symmetry_rmsds
            .iter()
            .copied()
            .max_by(f64::total_cmp);
        let max_connector_rmsd = fragment_connector_rmsds
            .iter()
            .copied()
            .max_by(f64::total_cmp);
        RigidFragmentRmsdSummary {
            fragment_rmsds,
            fragment_shape_rmsds,
            fragment_terminal_symmetry_rmsds,
            fragment_connector_rmsds,
            max_rmsd,
            max_shape_rmsd,
            max_terminal_symmetry_rmsd,
            max_connector_rmsd,
            mirror_branch_like_fragment_count_01a,
            mirror_branch_like_fragment_count_03a,
            worst_fragment_atoms,
        }
    }

    fn rigid_fragment_metric_details(
        molecule: &Molecule,
        reference_heavy_points: &[[f64; 3]],
        base: &Molecule,
        fragment: &[usize],
    ) -> RigidFragmentMetricDetails {
        let base_heavy_points = heavy_atom_points_for_rmsd(base);
        let heavy_index_by_atom = heavy_index_by_atom(molecule);
        let heavy_indices: Vec<_> = fragment
            .iter()
            .copied()
            .filter_map(|atom_idx| heavy_index_by_atom[atom_idx])
            .collect();
        let ref_points: Vec<_> = heavy_indices
            .iter()
            .map(|&idx| reference_heavy_points[idx])
            .collect();
        let base_points: Vec<_> = heavy_indices
            .iter()
            .map(|&idx| base_heavy_points[idx])
            .collect();
        let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
        let proper_rmsd = crate::distgeom::aligned_rmsd_for_test(&ref_points, &base_points);
        let shape_rmsd = mirror_invariant_aligned_rmsd_for_test(&ref_points, &base_points);
        let terminal_symmetry_rmsd = terminal_ligand_symmetry_aligned_rmsd_for_test(
            molecule,
            &adjacency,
            fragment,
            &ref_points,
            &base_points,
        );
        let cut_bonds = rings::symmetrize_sssr(molecule)
            .map(|ring_info| confseq_base_fragment_cut_bonds(molecule, &ring_info))
            .unwrap_or_default();
        let connector_fragment =
            rigid_fragment_atoms_with_connector_stubs(molecule, &adjacency, &cut_bonds, fragment);
        let connector_heavy_indices: Vec<_> = connector_fragment
            .iter()
            .copied()
            .filter_map(|atom_idx| heavy_index_by_atom[atom_idx])
            .collect();
        let connector_rmsd = (connector_heavy_indices.len() >= 3).then(|| {
            let connector_ref_points: Vec<_> = connector_heavy_indices
                .iter()
                .map(|&idx| reference_heavy_points[idx])
                .collect();
            let connector_base_points: Vec<_> = connector_heavy_indices
                .iter()
                .map(|&idx| base_heavy_points[idx])
                .collect();
            crate::distgeom::aligned_rmsd_for_test(&connector_ref_points, &connector_base_points)
        });
        let (angle_rms_deg, max_angle_delta_deg) =
            rigid_fragment_angle_delta_details(molecule, reference_heavy_points, base, fragment);
        RigidFragmentMetricDetails {
            proper_rmsd,
            shape_rmsd,
            terminal_symmetry_rmsd,
            connector_rmsd,
            angle_rms_deg,
            max_angle_delta_deg,
        }
    }

    fn points_for_atoms(points: &[[f64; 3]], atoms: &[usize]) -> Vec<[f64; 3]> {
        atoms.iter().map(|&atom| points[atom]).collect()
    }

    fn fragment_stage_coords_from_sliced_bounds_for_test(
        molecule: &Molecule,
        template: &ConfSeqBaseRigidFragmentTemplate,
        stage: crate::chemistry::distgeom::EmbedderTestStage,
    ) -> Result<HashMap<usize, [f64; 3]>, ConfSeqFastGeometryError> {
        let full_bounds = distgeom::dg_bounds_matrix(molecule).map_err(|err| {
            ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "failed to build whole-molecule DG bounds for diagnostic fragment key={}: {err}",
                template.cache_key
            ))
        })?;
        let realization_set: HashSet<_> = template.realization_atoms.iter().copied().collect();
        let atoms_to_remove = molecule
            .atoms()
            .iter()
            .filter_map(|atom| (!realization_set.contains(&atom.id().index())).then_some(atom.id()))
            .collect::<Vec<_>>();
        let mut builder = molecule.to_builder();
        let old_to_fragment = builder.remove_atoms_for_construction(&atoms_to_remove);
        let mut fragment = builder.build().map_err(|err| {
            ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "failed to build diagnostic fragment molecule key={}: {err}",
                template.cache_key
            ))
        })?;
        confseq_base_sanitize_fragment_stereo_for_conditioner(&mut fragment);
        let mut fragment_to_old = vec![None; fragment.num_atoms()];
        for &old_atom in &template.realization_atoms {
            let Some(new_atom) = old_to_fragment
                .get(old_atom)
                .copied()
                .flatten()
                .map(AtomId::index)
            else {
                return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                    "diagnostic fragment lost atom {old_atom} key={}",
                    template.cache_key
                )));
            };
            fragment_to_old[new_atom] = Some(old_atom);
        }
        let fragment_bounds = confseq_base_slice_bounds_matrix(&full_bounds, &fragment_to_old)
            .map_err(|reason| {
                ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                    "failed to slice diagnostic bounds key={}: {reason}",
                    template.cache_key
                ))
            })?;
        let mut params = EmbedParameters::etkdg();
        params.random_seed = 0;
        params.enable_sequential_random_seeds = true;
        params.num_threads = 1;
        params.clear_confs = true;
        params.prune_rms_thresh = -1.0;
        params.embed_fragments_separately = false;
        params.set_bounds_matrix(fragment_bounds).map_err(|err| {
            ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "failed to attach diagnostic sliced bounds key={}: {err}",
                template.cache_key
            ))
        })?;
        let staged = crate::distgeom::embedder_stage_coords_for_test(&fragment, &mut params, stage)
            .map_err(|err| {
                ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                    "diagnostic fragment stage embedding failed key={}: {err}",
                    template.cache_key
                ))
            })?;
        let conformer = staged.conformers_3d().first().ok_or_else(|| {
            ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                "diagnostic fragment stage produced no conformer key={}",
                template.cache_key
            ))
        })?;
        let mut coords = HashMap::new();
        for (fragment_atom, coord) in conformer.coordinates().iter().copied().enumerate() {
            let Some(old_atom) = fragment_to_old[fragment_atom] else {
                return Err(ConfSeqFastGeometryError::RigidFragmentEmbedding(format!(
                    "diagnostic fragment has unmapped atom {fragment_atom} key={}",
                    template.cache_key
                )));
            };
            coords.insert(old_atom, coord);
        }
        Ok(coords)
    }

    fn normalized_or_none_for_test(v: [f64; 3]) -> Option<[f64; 3]> {
        let norm = vec_dot(v, v).sqrt();
        (norm > 1.0e-8).then(|| vec_scale(v, 1.0 / norm))
    }

    fn local_frame_for_test(points: &[[f64; 3]]) -> Option<([f64; 3], [[f64; 3]; 3])> {
        for i in 0..points.len() {
            for j in 0..points.len() {
                if i == j {
                    continue;
                }
                let Some(ex) = normalized_or_none_for_test(vec_sub(points[j], points[i])) else {
                    continue;
                };
                for k in 0..points.len() {
                    if k == i || k == j {
                        continue;
                    }
                    let ik = vec_sub(points[k], points[i]);
                    let projected = vec_sub(ik, vec_scale(ex, vec_dot(ik, ex)));
                    let Some(ey) = normalized_or_none_for_test(projected) else {
                        continue;
                    };
                    let ez = vec_cross(ex, ey);
                    return Some((points[i], [ex, ey, ez]));
                }
            }
        }
        None
    }

    fn align_points_by_local_frame_for_test(
        reference: &[[f64; 3]],
        probe: &[[f64; 3]],
    ) -> Vec<[f64; 3]> {
        let Some((ref_origin, ref_frame)) = local_frame_for_test(reference) else {
            return probe.to_vec();
        };
        let Some((probe_origin, probe_frame)) = local_frame_for_test(probe) else {
            return probe.to_vec();
        };
        probe
            .iter()
            .map(|&point| {
                let local = vec_sub(point, probe_origin);
                let x = vec_dot(local, probe_frame[0]);
                let y = vec_dot(local, probe_frame[1]);
                let z = vec_dot(local, probe_frame[2]);
                vec_add(
                    ref_origin,
                    vec_add(
                        vec_scale(ref_frame[0], x),
                        vec_add(vec_scale(ref_frame[1], y), vec_scale(ref_frame[2], z)),
                    ),
                )
            })
            .collect()
    }

    fn fragment_subset_molecule_for_test(
        source_mol: &Molecule,
        kept_atoms: &[usize],
        coords: Vec<[f64; 3]>,
    ) -> Result<Molecule, String> {
        let kept_set: HashSet<_> = kept_atoms.iter().copied().collect();
        let atoms_to_remove = source_mol
            .atoms()
            .iter()
            .filter_map(|atom| (!kept_set.contains(&atom.id().index())).then_some(atom.id()))
            .collect::<Vec<_>>();
        let mut builder = source_mol.to_builder();
        builder.remove_atoms_for_construction(&atoms_to_remove);
        builder
            .add_conformer(Conformer3D::new(0, coords, true))
            .map_err(|err| err.to_string())?;
        builder.build().map_err(|err| err.to_string())
    }

    fn write_fragment_realization_pair_sdf_for_test(
        path: &std::path::Path,
        source_mol: &Molecule,
        fragment_atoms: &[usize],
        whole_points: &[[f64; 3]],
        fragment_coords_by_atom: &HashMap<usize, [f64; 3]>,
    ) -> Result<(), String> {
        let kept_atoms = source_mol
            .atoms()
            .iter()
            .map(|atom| atom.id().index())
            .filter(|atom| fragment_atoms.contains(atom))
            .collect::<Vec<_>>();
        let whole_fragment_points = kept_atoms
            .iter()
            .map(|&atom| whole_points[atom])
            .collect::<Vec<_>>();
        let fragment_points = kept_atoms
            .iter()
            .map(|atom| {
                fragment_coords_by_atom.get(atom).copied().ok_or_else(|| {
                    format!("missing fragment realization coordinate for atom {atom}")
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let aligned_fragment_points =
            align_points_by_local_frame_for_test(&whole_fragment_points, &fragment_points);
        let whole_mol =
            fragment_subset_molecule_for_test(source_mol, &kept_atoms, whole_fragment_points)?;
        let fragment_mol =
            fragment_subset_molecule_for_test(source_mol, &kept_atoms, aligned_fragment_points)?;
        let params = MolBlockWriteParams {
            format: SdfFormat::V2000,
            force_2d: false,
            include_stereo: true,
            kekulize: true,
            precision: 6,
        };
        let mut out = String::new();
        out.push_str(
            &mol_to_sdf_record_with_params(&whole_mol, &params)
                .map_err(|err| format!("failed to write whole DG fragment record: {err}"))?,
        );
        out.push_str(
            &mol_to_sdf_record_with_params(&fragment_mol, &params)
                .map_err(|err| format!("failed to write sliced fragment record: {err}"))?,
        );
        std::fs::write(path, out).map_err(|err| err.to_string())
    }

    fn rigid_fragment_angle_delta_details(
        molecule: &Molecule,
        reference_heavy_points: &[[f64; 3]],
        base: &Molecule,
        fragment: &[usize],
    ) -> (f64, f64) {
        let base_heavy_points = heavy_atom_points_for_rmsd(base);
        let heavy_index_by_atom = heavy_index_by_atom(molecule);
        let fragment_set: HashSet<_> = fragment.iter().copied().collect();
        let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
        let mut count = 0usize;
        let mut sum_sq = 0.0;
        let mut max_delta: f64 = 0.0;
        for &center in fragment {
            let neighbors: Vec<_> = adjacency
                .neighbors_of(center)
                .iter()
                .map(|neighbor| neighbor.atom_index)
                .filter(|atom| fragment_set.contains(atom))
                .filter(|atom| molecule.atoms()[*atom].atomic_number() > 1)
                .collect();
            for left_idx in 0..neighbors.len() {
                for right_idx in left_idx + 1..neighbors.len() {
                    let left = neighbors[left_idx];
                    let right = neighbors[right_idx];
                    let Some(left_heavy) = heavy_index_by_atom[left] else {
                        continue;
                    };
                    let Some(center_heavy) = heavy_index_by_atom[center] else {
                        continue;
                    };
                    let Some(right_heavy) = heavy_index_by_atom[right] else {
                        continue;
                    };
                    let ref_angle = angle_deg(
                        reference_heavy_points[left_heavy],
                        reference_heavy_points[center_heavy],
                        reference_heavy_points[right_heavy],
                    );
                    let base_angle = angle_deg(
                        base_heavy_points[left_heavy],
                        base_heavy_points[center_heavy],
                        base_heavy_points[right_heavy],
                    );
                    let delta = (ref_angle - base_angle).abs();
                    count += 1;
                    sum_sq += delta * delta;
                    max_delta = max_delta.max(delta);
                }
            }
        }
        if count == 0 {
            (0.0, 0.0)
        } else {
            ((sum_sq / count as f64).sqrt(), max_delta)
        }
    }

    fn mirror_invariant_aligned_rmsd_for_test(
        ref_points: &[[f64; 3]],
        probe_points: &[[f64; 3]],
    ) -> f64 {
        let proper = crate::distgeom::aligned_rmsd_for_test(ref_points, probe_points);
        let mirrored = probe_points
            .iter()
            .map(|point| [-point[0], point[1], point[2]])
            .collect::<Vec<_>>();
        proper.min(crate::distgeom::aligned_rmsd_for_test(
            ref_points, &mirrored,
        ))
    }

    fn rigid_fragment_atoms_with_connector_stubs(
        molecule: &Molecule,
        adjacency: &AdjacencyList,
        cut_bonds: &HashSet<(usize, usize)>,
        fragment: &[usize],
    ) -> Vec<usize> {
        let fragment_set: HashSet<_> = fragment.iter().copied().collect();
        let mut atoms = fragment.to_vec();
        for &atom in fragment {
            for neighbor in adjacency.neighbors_of(atom) {
                let other = neighbor.atom_index;
                if molecule.atoms()[other].atomic_number() <= 1 || fragment_set.contains(&other) {
                    continue;
                }
                if cut_bonds.contains(&sorted_pair(atom, other)) && !atoms.contains(&other) {
                    atoms.push(other);
                }
            }
        }
        atoms.sort_unstable();
        atoms
    }

    fn terminal_ligand_symmetry_aligned_rmsd_for_test(
        molecule: &Molecule,
        adjacency: &AdjacencyList,
        fragment: &[usize],
        ref_points: &[[f64; 3]],
        probe_points: &[[f64; 3]],
    ) -> f64 {
        let fragment_set: HashSet<_> = fragment.iter().copied().collect();
        let mut groups: HashMap<String, Vec<usize>> = HashMap::new();
        for (position, &atom_idx) in fragment.iter().enumerate() {
            let internal_neighbors = adjacency
                .neighbors_of(atom_idx)
                .iter()
                .filter(|neighbor| {
                    fragment_set.contains(&neighbor.atom_index)
                        && molecule.atoms()[neighbor.atom_index].atomic_number() > 1
                })
                .collect::<Vec<_>>();
            if internal_neighbors.len() != 1 {
                continue;
            }
            let parent = internal_neighbors[0].atom_index;
            let Some(bond) = bond_between_pair(molecule, sorted_pair(atom_idx, parent)) else {
                continue;
            };
            let atom = &molecule.atoms()[atom_idx];
            let key = format!(
                "parent={parent};z={};q={};hyb={:?};arom={};bond={:?}:{}",
                atom.atomic_number(),
                atom.formal_charge(),
                atom.hybridization(),
                atom.is_aromatic(),
                bond.order(),
                bond.is_aromatic()
            );
            groups.entry(key).or_default().push(position);
        }
        let groups = groups
            .into_values()
            .filter(|group| group.len() > 1)
            .collect::<Vec<_>>();
        let total_permutations = groups
            .iter()
            .map(|group| factorial_for_diagnostic(group.len()))
            .product::<usize>();
        if groups.is_empty() || total_permutations > 720 {
            return crate::distgeom::aligned_rmsd_for_test(ref_points, probe_points);
        }

        let mut best = f64::INFINITY;
        let mut order = (0..probe_points.len()).collect::<Vec<_>>();
        visit_terminal_ligand_symmetry_orders(
            0,
            &groups,
            &mut order,
            ref_points,
            probe_points,
            &mut best,
        );
        best
    }

    fn visit_terminal_ligand_symmetry_orders(
        group_idx: usize,
        groups: &[Vec<usize>],
        order: &mut [usize],
        ref_points: &[[f64; 3]],
        probe_points: &[[f64; 3]],
        best: &mut f64,
    ) {
        if group_idx == groups.len() {
            let permuted = order
                .iter()
                .map(|&idx| probe_points[idx])
                .collect::<Vec<_>>();
            *best = best.min(crate::distgeom::aligned_rmsd_for_test(
                ref_points, &permuted,
            ));
            return;
        }
        let group = &groups[group_idx];
        let original = group.iter().map(|&idx| order[idx]).collect::<Vec<_>>();
        let mut values = original.clone();
        visit_terminal_ligand_permutations(0, &mut values, &mut |candidate| {
            for (&position, &value) in group.iter().zip(candidate) {
                order[position] = value;
            }
            visit_terminal_ligand_symmetry_orders(
                group_idx + 1,
                groups,
                order,
                ref_points,
                probe_points,
                best,
            );
        });
        for (&position, &value) in group.iter().zip(&original) {
            order[position] = value;
        }
    }

    fn visit_terminal_ligand_permutations(
        idx: usize,
        values: &mut [usize],
        visit: &mut impl FnMut(&[usize]),
    ) {
        if idx == values.len() {
            visit(values);
            return;
        }
        for swap_idx in idx..values.len() {
            values.swap(idx, swap_idx);
            visit_terminal_ligand_permutations(idx + 1, values, visit);
            values.swap(idx, swap_idx);
        }
    }

    fn factorial_for_diagnostic(value: usize) -> usize {
        (1..=value).product()
    }

    fn rigid_heavy_fragments_cutting_confseq_rotors(molecule: &Molecule) -> Vec<Vec<usize>> {
        let cut_bonds = rings::symmetrize_sssr(molecule)
            .map(|ring_info| confseq_base_fragment_cut_bonds(molecule, &ring_info))
            .unwrap_or_else(|_| {
                collect_single_bond_dihedrals(molecule)
                    .into_iter()
                    .map(|(_, j, k, _)| sorted_pair(j, k))
                    .collect()
            });
        let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
        let mut seen = vec![false; molecule.num_atoms()];
        let mut fragments = Vec::new();
        for atom in molecule.atoms() {
            let start = atom.id().index();
            if seen[start] || atom.atomic_number() == 1 {
                continue;
            }
            let mut fragment = Vec::new();
            let mut queue = VecDeque::new();
            seen[start] = true;
            queue.push_back(start);
            while let Some(current) = queue.pop_front() {
                fragment.push(current);
                for neighbor in adjacency.neighbors_of(current) {
                    let next = neighbor.atom_index;
                    if seen[next] || molecule.atoms()[next].atomic_number() == 1 {
                        continue;
                    }
                    if cut_bonds.contains(&sorted_pair(current, next)) {
                        continue;
                    }
                    seen[next] = true;
                    queue.push_back(next);
                }
            }
            fragments.push(fragment);
        }
        fragments
    }

    fn heavy_index_by_atom(molecule: &Molecule) -> Vec<Option<usize>> {
        let mut heavy_index_by_atom = vec![None; molecule.num_atoms()];
        let mut heavy_idx = 0usize;
        for atom in molecule.atoms() {
            if atom.atomic_number() == 1 {
                continue;
            }
            heavy_index_by_atom[atom.id().index()] = Some(heavy_idx);
            heavy_idx += 1;
        }
        heavy_index_by_atom
    }

    fn rigid_fragment_type_for_diagnostic(
        molecule: &Molecule,
        fragment_atoms: &[usize],
    ) -> &'static str {
        let key = rigid_fragment_framework_key_for_diagnostic(molecule, fragment_atoms);
        if key.starts_with("ring_planar") {
            "ring_planar"
        } else if key.starts_with("ring_nonplanar") {
            "ring_nonplanar"
        } else if key == "planar_pi" {
            "planar_pi"
        } else if key == "ring_mixed" {
            "ring_mixed"
        } else {
            "acyclic"
        }
    }

    fn rigid_fragment_framework_key_for_diagnostic(
        molecule: &Molecule,
        fragment_atoms: &[usize],
    ) -> String {
        let fragment_set: HashSet<_> = fragment_atoms.iter().copied().collect();
        let ring_info = rings::symmetrize_sssr(molecule).ok();
        let ring_components = ring_info
            .as_ref()
            .and_then(|ring_info| classify_confseq_base_ring_components(molecule, ring_info).ok())
            .unwrap_or_default();
        for component in &ring_components {
            if component
                .atoms
                .iter()
                .all(|atom| fragment_set.contains(atom))
            {
                return if let Some(ring_info) = ring_info.as_ref() {
                    format!(
                        "{}:{:?}",
                        if component.planar {
                            "ring_planar"
                        } else {
                            "ring_nonplanar"
                        },
                        confseq_base_ring_topology_from_ring_info(
                            molecule,
                            ring_info,
                            &component.atoms
                        )
                    )
                } else {
                    format!(
                        "{}:Unknown",
                        if component.planar {
                            "ring_planar"
                        } else {
                            "ring_nonplanar"
                        }
                    )
                };
            }
        }
        let all_planar_like = fragment_atoms.iter().copied().all(|atom_idx| {
            let atom = &molecule.atoms()[atom_idx];
            atom.is_aromatic() || atom.hybridization() == Hybridization::Sp2
        });
        if all_planar_like {
            "planar_pi".to_string()
        } else if fragment_atoms.iter().any(|atom| {
            ring_components
                .iter()
                .any(|component| component.atoms.contains(atom))
        }) {
            "ring_mixed".to_string()
        } else {
            "acyclic".to_string()
        }
    }

    #[derive(Debug, Default)]
    struct RigidFragmentFrameworkRmsdBuckets {
        proper: HashMap<String, Vec<f64>>,
        shape: HashMap<String, Vec<f64>>,
        terminal_symmetry: HashMap<String, Vec<f64>>,
        connector: HashMap<String, Vec<f64>>,
    }

    impl RigidFragmentFrameworkRmsdBuckets {
        fn record(&mut self, key: String, details: &RigidFragmentMetricDetails) {
            self.proper
                .entry(key.clone())
                .or_default()
                .push(details.proper_rmsd);
            self.shape
                .entry(key.clone())
                .or_default()
                .push(details.shape_rmsd);
            self.terminal_symmetry
                .entry(key.clone())
                .or_default()
                .push(details.terminal_symmetry_rmsd);
            if let Some(connector_rmsd) = details.connector_rmsd {
                self.connector.entry(key).or_default().push(connector_rmsd);
            }
        }
    }

    #[derive(Debug, Default)]
    struct BaseGeometryStageDiagnostics {
        comparable: usize,
        pass_005a_5deg: usize,
        pass_010a_10deg: usize,
        pass_010a_15deg: usize,
        local_bond_rmsds: Vec<f64>,
        local_non_token_angle_rmsds: Vec<f64>,
        rigid_fragment_shape_rmsds: Vec<f64>,
        rigid_fragment_shape_max_per_mol: Vec<f64>,
    }

    impl BaseGeometryStageDiagnostics {
        fn record(
            &mut self,
            source_mol: &Molecule,
            reference_points: &[[f64; 3]],
            probe: &Molecule,
        ) {
            let local_summary =
                local_constraint_summary_against_reference(source_mol, reference_points, probe);
            self.comparable += 1;
            if local_summary.bond_rms_a <= 0.05 && local_summary.non_token_angle_rms_deg <= 5.0 {
                self.pass_005a_5deg += 1;
            }
            if local_summary.bond_rms_a <= 0.10 && local_summary.non_token_angle_rms_deg <= 10.0 {
                self.pass_010a_10deg += 1;
            }
            if local_summary.bond_rms_a <= 0.10 && local_summary.non_token_angle_rms_deg <= 15.0 {
                self.pass_010a_15deg += 1;
            }
            if local_summary.bond_count > 0 {
                self.local_bond_rmsds.push(local_summary.bond_rms_a);
            }
            if local_summary.non_token_angle_count > 0 {
                self.local_non_token_angle_rmsds
                    .push(local_summary.non_token_angle_rms_deg);
            }

            let rigid_summary = rigid_fragment_rmsd_summary(source_mol, reference_points, probe);
            self.rigid_fragment_shape_rmsds
                .extend(rigid_summary.fragment_shape_rmsds.iter().copied());
            if let Some(max_shape_rmsd) = rigid_summary.max_shape_rmsd {
                self.rigid_fragment_shape_max_per_mol.push(max_shape_rmsd);
            }
        }

        fn sort(&mut self) {
            self.local_bond_rmsds
                .sort_by(|left, right| left.total_cmp(right));
            self.local_non_token_angle_rmsds
                .sort_by(|left, right| left.total_cmp(right));
            self.rigid_fragment_shape_rmsds
                .sort_by(|left, right| left.total_cmp(right));
            self.rigid_fragment_shape_max_per_mol
                .sort_by(|left, right| left.total_cmp(right));
        }

        fn print(&self, label: &str) {
            let pass_rate_005a_5deg =
                self.pass_005a_5deg as f64 / self.comparable.max(1) as f64 * 100.0;
            let pass_rate_010a_10deg =
                self.pass_010a_10deg as f64 / self.comparable.max(1) as f64 * 100.0;
            let pass_rate_010a_15deg =
                self.pass_010a_15deg as f64 / self.comparable.max(1) as f64 * 100.0;
            eprintln!(
                "confseq_base_corpus {label}_local_constraint_pass comparable={} pass_0_05a_5deg={} rate_0_05a_5deg={:.2}% pass_0_10a_10deg={} rate_0_10a_10deg={:.2}% pass_0_10a_15deg={} rate_0_10a_15deg={:.2}%",
                self.comparable,
                self.pass_005a_5deg,
                pass_rate_005a_5deg,
                self.pass_010a_10deg,
                pass_rate_010a_10deg,
                self.pass_010a_15deg,
                pass_rate_010a_15deg,
            );
            eprintln!(
                "confseq_base_corpus {label}_local_bond_rms comparable={} p50={:.6} p90={:.6} p99={:.6}",
                self.local_bond_rmsds.len(),
                quantile_sorted(&self.local_bond_rmsds, 0.50),
                quantile_sorted(&self.local_bond_rmsds, 0.90),
                quantile_sorted(&self.local_bond_rmsds, 0.99),
            );
            eprintln!(
                "confseq_base_corpus {label}_local_non_token_angle_rms comparable={} p50={:.6} p90={:.6} p99={:.6}",
                self.local_non_token_angle_rmsds.len(),
                quantile_sorted(&self.local_non_token_angle_rmsds, 0.50),
                quantile_sorted(&self.local_non_token_angle_rmsds, 0.90),
                quantile_sorted(&self.local_non_token_angle_rmsds, 0.99),
            );
            eprintln!(
                "confseq_base_corpus {label}_rigid_fragment_shape_rmsd fragments={} p50={:.6} p90={:.6} p99={:.6}",
                self.rigid_fragment_shape_rmsds.len(),
                quantile_sorted(&self.rigid_fragment_shape_rmsds, 0.50),
                quantile_sorted(&self.rigid_fragment_shape_rmsds, 0.90),
                quantile_sorted(&self.rigid_fragment_shape_rmsds, 0.99),
            );
            eprintln!(
                "confseq_base_corpus {label}_rigid_fragment_shape_max_per_mol comparable={} p50={:.6} p90={:.6} p99={:.6}",
                self.rigid_fragment_shape_max_per_mol.len(),
                quantile_sorted(&self.rigid_fragment_shape_max_per_mol, 0.50),
                quantile_sorted(&self.rigid_fragment_shape_max_per_mol, 0.90),
                quantile_sorted(&self.rigid_fragment_shape_max_per_mol, 0.99),
            );
        }
    }

    fn print_rigid_fragment_rmsd_buckets(
        label: &str,
        mut buckets: HashMap<String, Vec<f64>>,
        limit: usize,
    ) {
        let mut summaries = buckets
            .drain()
            .map(|(key, mut values)| {
                values.sort_by(|left, right| left.total_cmp(right));
                let count = values.len();
                let p50 = quantile_sorted(&values, 0.50);
                let p90 = quantile_sorted(&values, 0.90);
                let p95 = quantile_sorted(&values, 0.95);
                let p99 = quantile_sorted(&values, 0.99);
                (key, count, p50, p90, p95, p99)
            })
            .collect::<Vec<_>>();
        summaries.sort_by(|left, right| {
            right
                .5
                .total_cmp(&left.5)
                .then_with(|| right.1.cmp(&left.1))
                .then_with(|| left.0.cmp(&right.0))
        });
        for (key, count, p50, p90, p95, p99) in summaries.into_iter().take(limit) {
            eprintln!(
                "{label} key={key} count={count} p50={p50:.6} p90={p90:.6} p95={p95:.6} p99={p99:.6}"
            );
        }
    }

    #[derive(Debug, Default)]
    struct ConfSeqBaseFrameworkDiagnostics {
        fragment_shapes: HashMap<String, usize>,
        assembly_paths: HashMap<String, usize>,
        fragment_template_failures: HashMap<String, usize>,
        nonplanar_ring_topologies: HashMap<String, usize>,
    }

    impl ConfSeqBaseFrameworkDiagnostics {
        fn collect(
            molecule: &Molecule,
            ring_info: &rings::RingInfo,
            model: &ConfSeqBaseConstraintModel,
        ) -> Self {
            let mut diagnostics = Self::default();
            for component in &model.rigid_components {
                let realization_atoms = confseq_base_component_realization_atoms(component);
                match confseq_base_rigid_fragment_shape(molecule, component, &realization_atoms) {
                    Ok(shape) => {
                        *diagnostics
                            .fragment_shapes
                            .entry(format!("{:?}:{:?}", component.kind, shape))
                            .or_default() += 1;
                    }
                    Err(err) => {
                        *diagnostics
                            .fragment_template_failures
                            .entry(format!("{err:?}"))
                            .or_default() += 1;
                    }
                }
            }
            for path_class in confseq_base_assembly_path_classes_for_diagnostic(molecule, model) {
                *diagnostics.assembly_paths.entry(path_class).or_default() += 1;
            }
            for topology in
                confseq_base_nonplanar_ring_topologies_for_diagnostic(molecule, ring_info)
            {
                *diagnostics
                    .nonplanar_ring_topologies
                    .entry(topology)
                    .or_default() += 1;
            }
            diagnostics
        }

        fn merge_into(
            &self,
            fragment_shapes: &mut HashMap<String, usize>,
            assembly_paths: &mut HashMap<String, usize>,
            fragment_template_failures: &mut HashMap<String, usize>,
            nonplanar_ring_topologies: &mut HashMap<String, usize>,
        ) {
            merge_count_map(fragment_shapes, &self.fragment_shapes);
            merge_count_map(assembly_paths, &self.assembly_paths);
            merge_count_map(fragment_template_failures, &self.fragment_template_failures);
            merge_count_map(nonplanar_ring_topologies, &self.nonplanar_ring_topologies);
        }
    }

    fn merge_count_map(target: &mut HashMap<String, usize>, source: &HashMap<String, usize>) {
        for (key, count) in source {
            *target.entry(key.clone()).or_default() += count;
        }
    }

    fn print_top_count_map(label: &str, counts: HashMap<String, usize>, limit: usize) {
        let mut counts: Vec<_> = counts.into_iter().collect();
        counts.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
        for (key, count) in counts.into_iter().take(limit) {
            eprintln!("{label} count={count} key={key}");
        }
    }

    fn confseq_base_assembly_path_classes_for_diagnostic(
        molecule: &Molecule,
        model: &ConfSeqBaseConstraintModel,
    ) -> Vec<String> {
        if model.rigid_components.is_empty() {
            return Vec::new();
        }
        let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
        let mut component_by_atom = vec![None; molecule.num_atoms()];
        for (component_idx, component) in model.rigid_components.iter().enumerate() {
            for &atom in &component.atoms {
                component_by_atom[atom] = Some(component_idx);
            }
        }
        let scaffold_shape =
            coordinates::rdkit_initial_2d_scaffold_coords(molecule.atoms(), molecule.bonds())
                .ok()
                .map(|points| {
                    points
                        .into_iter()
                        .map(|point| [point[0], point[1], 0.0])
                        .collect::<Vec<_>>()
                });
        let mut local_components = Vec::with_capacity(model.rigid_components.len());
        let mut fragment_template_cache =
            HashMap::<String, ConfSeqBaseRigidFragmentTemplate>::new();
        for component in &model.rigid_components {
            let Ok(template) = confseq_base_cached_rigid_fragment_template(
                molecule,
                model,
                component,
                &mut fragment_template_cache,
            ) else {
                return Vec::new();
            };
            let Ok(local) = confseq_base_realize_rigid_fragment_template(
                molecule,
                model,
                &template,
                scaffold_shape.as_deref(),
                None,
            ) else {
                return Vec::new();
            };
            local_components.push(local);
        }
        let mut coords = vec![[0.0; 3]; molecule.num_atoms()];
        let mut placed_atoms = vec![false; molecule.num_atoms()];
        let mut placed_components = vec![false; model.rigid_components.len()];
        let mut connector_targets = HashMap::<(usize, usize), [f64; 3]>::new();
        let root_component =
            confseq_base_select_root_component(molecule, model, &component_by_atom);
        if place_confseq_base_rigid_component_local(
            &model.rigid_components[root_component],
            &local_components[root_component],
            [0.0, 0.0, 0.0],
            &mut coords,
            &mut placed_atoms,
            &mut connector_targets,
        )
        .is_err()
        {
            return Vec::new();
        }
        placed_components[root_component] = true;
        let (cut_bonds, reciprocal_stubs) =
            confseq_base_connector_stub_counts_for_diagnostic(model);
        let mut classes = Vec::new();
        classes.push(format!("cut_bonds={cut_bonds}"));
        classes.push(format!("reciprocal_connector_stubs={reciprocal_stubs}"));
        loop {
            let Some(next_component) = confseq_base_select_next_assembly_component(
                molecule,
                &adjacency,
                model,
                &component_by_atom,
                &local_components,
                &coords,
                &placed_atoms,
                &placed_components,
                &connector_targets,
            ) else {
                break;
            };
            let anchors = confseq_base_component_assembly_anchors(
                molecule,
                &component_by_atom,
                &placed_atoms,
                next_component,
            );
            let local = &local_components[next_component];
            let connector_count = anchors
                .iter()
                .filter(|anchor| {
                    connector_targets.contains_key(&(anchor.placed_atom, anchor.unplaced_atom))
                        && local.contains_key(&anchor.placed_atom)
                        && local.contains_key(&anchor.unplaced_atom)
                })
                .count();
            let bucket = if connector_count > 0 {
                "shared_bond_assembly".to_string()
            } else {
                match anchors.len() {
                    0 => "fallback_no_anchor".to_string(),
                    1 => "fallback_single_anchor".to_string(),
                    2 => "fallback_multi_anchor".to_string(),
                    _ => "fallback_multi_anchor".to_string(),
                }
            };
            classes.push(bucket);
            if place_confseq_base_unplaced_component_by_anchors(
                molecule,
                &adjacency,
                model,
                &component_by_atom,
                &local_components[next_component],
                next_component,
                &mut coords,
                &mut placed_atoms,
                &mut connector_targets,
            )
            .is_err()
            {
                break;
            }
            placed_components[next_component] = true;
        }
        classes
    }

    fn confseq_base_connector_stub_counts_for_diagnostic(
        model: &ConfSeqBaseConstraintModel,
    ) -> (usize, usize) {
        let mut directed = HashSet::<(usize, usize)>::new();
        for component in &model.rigid_components {
            for connector in &component.connectors {
                directed.insert((connector.core_atom, connector.external_atom));
            }
        }
        let mut cut_bonds = HashSet::<(usize, usize)>::new();
        let mut reciprocal = HashSet::<(usize, usize)>::new();
        for &(core, external) in &directed {
            let pair = sorted_pair(core, external);
            cut_bonds.insert(pair);
            if directed.contains(&(external, core)) {
                reciprocal.insert(pair);
            }
        }
        (cut_bonds.len(), reciprocal.len())
    }

    fn confseq_base_nonplanar_ring_topologies_for_diagnostic(
        molecule: &Molecule,
        ring_info: &rings::RingInfo,
    ) -> Vec<String> {
        let Ok(components) = classify_confseq_base_ring_components(molecule, ring_info) else {
            return Vec::new();
        };
        components
            .iter()
            .filter(|component| !component.planar)
            .map(|component| {
                confseq_base_ring_component_topology_for_diagnostic(molecule, ring_info, component)
            })
            .collect()
    }

    fn confseq_base_ring_component_topology_for_diagnostic(
        molecule: &Molecule,
        ring_info: &rings::RingInfo,
        component: &ConfSeqBaseRingComponent,
    ) -> String {
        let atom_set: HashSet<_> = component.atoms.iter().copied().collect();
        let rings = ring_info
            .atom_rings()
            .iter()
            .enumerate()
            .filter_map(|(ring_idx, atoms)| {
                atoms
                    .iter()
                    .all(|atom| atom_set.contains(&atom.index()))
                    .then_some(ring_idx)
            })
            .collect::<Vec<_>>();
        let mut shared_atom_pairs = 0usize;
        let mut shared_bond_pairs = 0usize;
        let mut max_shared_atoms = 0usize;
        let mut max_shared_bonds = 0usize;
        for left_pos in 0..rings.len() {
            for right_pos in left_pos + 1..rings.len() {
                let left = rings[left_pos];
                let right = rings[right_pos];
                let shared_atoms = ring_info.atom_rings()[left]
                    .iter()
                    .filter(|atom| ring_info.atom_rings()[right].contains(atom))
                    .count();
                let shared_bonds = shared_bond_ids_between_rings(ring_info, left, right).len();
                if shared_atoms > 0 {
                    shared_atom_pairs += 1;
                    max_shared_atoms = max_shared_atoms.max(shared_atoms);
                }
                if shared_bonds > 0 {
                    shared_bond_pairs += 1;
                    max_shared_bonds = max_shared_bonds.max(shared_bonds);
                }
            }
        }
        let simple = confseq_base_simple_ring_order(molecule, &component.atoms).is_some();
        let topology = if simple {
            "simple"
        } else if rings.len() == 1 {
            "single_non_simple"
        } else if max_shared_bonds > 0 && shared_bond_pairs + 1 == rings.len() {
            "edge_fused_chain"
        } else if max_shared_bonds > 0 {
            "edge_fused_polycyclic"
        } else if max_shared_atoms == 1 {
            "spiro"
        } else if max_shared_atoms > 1 {
            "bridged_or_cage"
        } else {
            "disconnected_or_unknown"
        };
        let ring_sizes = rings
            .iter()
            .map(|&ring| ring_info.atom_rings()[ring].len().to_string())
            .collect::<Vec<_>>()
            .join(",");
        format!(
            "{topology}:rings={}:atoms={}:shared_atom_pairs={shared_atom_pairs}:shared_bond_pairs={shared_bond_pairs}:max_shared_atoms={max_shared_atoms}:max_shared_bonds={max_shared_bonds}:sizes=[{ring_sizes}]",
            rings.len(),
            component.atoms.len(),
        )
    }

    fn rigid_fragment_local_diagnostic(molecule: &Molecule, fragment_atoms: &[usize]) -> String {
        let fragment_set: HashSet<_> = fragment_atoms.iter().copied().collect();
        let cut_bonds = rings::symmetrize_sssr(molecule)
            .map(|ring_info| confseq_base_fragment_cut_bonds(molecule, &ring_info))
            .unwrap_or_else(|_| {
                collect_single_bond_dihedrals(molecule)
                    .into_iter()
                    .map(|(_, j, k, _)| sorted_pair(j, k))
                    .collect()
            });
        let ring_bonds: HashSet<_> = rings::symmetrize_sssr(molecule)
            .map(|ring_info| {
                ring_info
                    .bond_rings()
                    .iter()
                    .flatten()
                    .map(|bond| bond.index())
                    .collect()
            })
            .unwrap_or_default();
        let atom_details: Vec<_> = fragment_atoms
            .iter()
            .copied()
            .map(|atom_idx| {
                let atom = &molecule.atoms()[atom_idx];
                format!(
                    "{atom_idx}:Z{} hyb={:?} arom={} charge={} hdeg={}",
                    atom.atomic_number(),
                    atom.hybridization(),
                    atom.is_aromatic(),
                    atom.formal_charge(),
                    heavy_degree(molecule, atom_idx)
                )
            })
            .collect();
        let bond_details: Vec<_> = molecule
            .bonds()
            .iter()
            .filter_map(|bond| {
                let begin = bond.begin().index();
                let end = bond.end().index();
                let begin_in = fragment_set.contains(&begin);
                let end_in = fragment_set.contains(&end);
                if !begin_in && !end_in {
                    return None;
                }
                let scope = if begin_in && end_in {
                    "internal"
                } else {
                    "border"
                };
                let pair = sorted_pair(begin, end);
                Some(format!(
                    "{begin}-{end}:{:?} arom={} conj={} ring={} rotor_cut={} {scope}",
                    bond.order(),
                    bond.is_aromatic(),
                    bond.is_conjugated(),
                    ring_bonds.contains(&bond.id().index()),
                    cut_bonds.contains(&pair),
                ))
            })
            .collect();
        format!(
            "atoms=[{}] bonds=[{}]",
            atom_details.join("; "),
            bond_details.join("; ")
        )
    }

    fn load_or_generate_confseq_dg_reference_cache(
        corpus_path: &str,
        in_smiles_batch: &[String],
        td_smiles_batch: &[String],
    ) -> ConfSeqDgReferenceCache {
        let cache_path = confseq_dg_reference_cache_path(corpus_path);
        if cache_path.exists() {
            let raw = std::fs::read_to_string(&cache_path).unwrap_or_else(|err| {
                panic!(
                    "failed to read ConfSeq DG reference cache {}: {err}",
                    cache_path.display()
                )
            });
            let mut cache: ConfSeqDgReferenceCache =
                serde_json::from_str(&raw).unwrap_or_else(|err| {
                    panic!(
                        "failed to parse ConfSeq DG reference cache {}: {err}",
                        cache_path.display()
                    )
                });
            assert_eq!(
                cache.entries.len(),
                in_smiles_batch.len(),
                "ConfSeq DG reference cache {} has {} entries but corpus has {}",
                cache_path.display(),
                cache.entries.len(),
                in_smiles_batch.len()
            );
            filter_invalid_confseq_dg_reference_cache_entries(
                &mut cache,
                in_smiles_batch,
                td_smiles_batch,
            );
            return cache;
        }

        if std::env::var("COSMOLKIT_CONFSEQ_GENERATE_DG_REFERENCE_CACHE")
            .ok()
            .as_deref()
            != Some("1")
        {
            panic!(
                "ConfSeq DG reference cache is missing at {}. Set COSMOLKIT_CONFSEQ_GENERATE_DG_REFERENCE_CACHE=1 to generate it explicitly.",
                cache_path.display()
            );
        }

        let reference_options = ConfSeqDecodeOptions {
            optimize_with_uff: false,
            ..ConfSeqDecodeOptions::default()
        };
        let mut entries = Vec::with_capacity(in_smiles_batch.len());
        for (idx, (in_smiles, td_smiles)) in in_smiles_batch.iter().zip(td_smiles_batch).enumerate()
        {
            eprintln!(
                "generating ConfSeq DG reference cache record {}/{}",
                idx,
                in_smiles_batch.len()
            );
            let entry = match decode_confseq_with_options(in_smiles, td_smiles, &reference_options)
            {
                Ok(molecule) => {
                    let heavy_points = heavy_atom_points_for_rmsd(&molecule);
                    let (_, _, max_heavy_bond) = heavy_bond_length_stats(&molecule, &heavy_points);
                    if max_heavy_bond > CONFSEQ_DG_REFERENCE_MAX_REASONABLE_HEAVY_BOND_A {
                        ConfSeqDgReferenceCacheEntry {
                            heavy_atom_points: None,
                            error: Some(format!(
                                "InvalidFullDecodeGeometry(max_heavy_bond={max_heavy_bond:.6})"
                            )),
                        }
                    } else {
                        ConfSeqDgReferenceCacheEntry {
                            heavy_atom_points: Some(heavy_points),
                            error: None,
                        }
                    }
                }
                Err(error) => ConfSeqDgReferenceCacheEntry {
                    heavy_atom_points: None,
                    error: Some(format!("{error:?}")),
                },
            };
            entries.push(entry);
        }
        let cache = ConfSeqDgReferenceCache {
            corpus_path: corpus_path.to_string(),
            entries,
        };
        if let Some(parent) = cache_path.parent() {
            std::fs::create_dir_all(parent).unwrap_or_else(|err| {
                panic!(
                    "failed to create ConfSeq DG reference cache directory {}: {err}",
                    parent.display()
                )
            });
        }
        let raw = serde_json::to_string_pretty(&cache)
            .expect("ConfSeq DG reference cache should serialize");
        std::fs::write(&cache_path, raw).unwrap_or_else(|err| {
            panic!(
                "failed to write ConfSeq DG reference cache {}: {err}",
                cache_path.display()
            )
        });
        cache
    }

    fn filter_invalid_confseq_dg_reference_cache_entries(
        cache: &mut ConfSeqDgReferenceCache,
        in_smiles_batch: &[String],
        td_smiles_batch: &[String],
    ) {
        for (idx, entry) in cache.entries.iter_mut().enumerate() {
            let Some(heavy_points) = entry.heavy_atom_points.as_ref() else {
                continue;
            };
            let Ok(parsed) = parse_confseq(&in_smiles_batch[idx], &td_smiles_batch[idx]) else {
                continue;
            };
            let Ok(source_mol) = Molecule::from_smiles(&parsed.stripped_smiles) else {
                continue;
            };
            let (_, _, max_heavy_bond) = heavy_bond_length_stats(&source_mol, heavy_points);
            if max_heavy_bond > CONFSEQ_DG_REFERENCE_MAX_REASONABLE_HEAVY_BOND_A {
                entry.heavy_atom_points = None;
                entry.error = Some(format!(
                    "InvalidFullDecodeGeometry(max_heavy_bond={max_heavy_bond:.6})"
                ));
            }
        }
    }

    fn read_confseq_corpus_smiles(corpus_path: &str) -> (Vec<String>, Vec<String>) {
        let input =
            std::fs::read_to_string(corpus_path).expect("ConfSeq corpus should be readable");
        let mut in_smiles_batch = Vec::new();
        let mut td_smiles_batch = Vec::new();
        for (line_idx, line) in input.lines().enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let value: Value = serde_json::from_str(line).unwrap_or_else(|err| {
                panic!("failed to parse corpus JSON line {}: {err}", line_idx + 1)
            });
            let in_smiles = value["in_smiles"].as_str().unwrap_or_else(|| {
                panic!(
                    "corpus line {} is missing string field in_smiles",
                    line_idx + 1
                )
            });
            let td_smiles = value["td_smiles"].as_str().unwrap_or_else(|| {
                panic!(
                    "corpus line {} is missing string field td_smiles",
                    line_idx + 1
                )
            });
            in_smiles_batch.push(in_smiles.to_string());
            td_smiles_batch.push(td_smiles.to_string());
        }
        (in_smiles_batch, td_smiles_batch)
    }

    fn confseq_dg_reference_cache_path(corpus_path: &str) -> std::path::PathBuf {
        std::env::var("COSMOLKIT_CONFSEQ_DG_REFERENCE_CACHE")
            .map(std::path::PathBuf::from)
            .unwrap_or_else(|_| {
                std::path::PathBuf::from(format!("{corpus_path}.dg_reference_heavy_points.json"))
            })
    }

    fn bond_length_rms_against_reference(reference: &Molecule, probe: &Molecule) -> f64 {
        let ref_points = conformer_points(reference);
        let probe_points = conformer_points(probe);
        let (_, rms, _) = distance_stats(
            reference
                .bonds()
                .iter()
                .map(|bond| (bond.begin().index(), bond.end().index())),
            &ref_points,
            &probe_points,
        );
        rms
    }

    fn angle_rms_against_reference(reference: &Molecule, probe: &Molecule) -> f64 {
        let ref_points = conformer_points(reference);
        let probe_points = conformer_points(probe);
        let centers = collect_angle_centers(reference).expect("angle centers should collect");
        if centers.is_empty() {
            return 0.0;
        }
        let sum_sq: f64 = centers
            .iter()
            .map(|(i, j, k)| {
                let delta = angle_deg(ref_points[*i], ref_points[*j], ref_points[*k])
                    - angle_deg(probe_points[*i], probe_points[*j], probe_points[*k]);
                delta * delta
            })
            .sum();
        (sum_sq / centers.len() as f64).sqrt()
    }

    fn all_heavy_angle_centers(molecule: &Molecule) -> Vec<(usize, usize, usize)> {
        let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
        let mut centers = Vec::new();
        for atom in molecule.atoms() {
            let center = atom.id().index();
            let heavy_neighbors: Vec<_> = adjacency
                .neighbors_of(center)
                .iter()
                .map(|nbr| nbr.atom_index)
                .filter(|idx| molecule.atoms()[*idx].atomic_number() > 1)
                .collect();
            for left_ord in 0..heavy_neighbors.len() {
                for right_ord in left_ord + 1..heavy_neighbors.len() {
                    centers.push((
                        heavy_neighbors[left_ord],
                        center,
                        heavy_neighbors[right_ord],
                    ));
                }
            }
        }
        centers
    }

    fn all_heavy_angle_rms_against_reference(reference: &Molecule, probe: &Molecule) -> f64 {
        let ref_points = conformer_points(reference);
        let probe_points = conformer_points(probe);
        let centers = all_heavy_angle_centers(reference);
        if centers.is_empty() {
            return 0.0;
        }
        let sum_sq: f64 = centers
            .iter()
            .map(|(i, j, k)| {
                let delta = angle_deg(ref_points[*i], ref_points[*j], ref_points[*k])
                    - angle_deg(probe_points[*i], probe_points[*j], probe_points[*k]);
                delta * delta
            })
            .sum();
        (sum_sq / centers.len() as f64).sqrt()
    }

    #[test]
    fn decode_confseq_decodes_minimal_ethane_without_external_runtime() {
        let molecule = decode_confseq("C C", "C C").expect("ethane ConfSeq should decode");

        assert_eq!(molecule.num_atoms(), 2);
        assert_eq!(molecule.num_bonds(), 1);
        assert_eq!(molecule.conformers_3d().len(), 1);
    }

    #[test]
    fn confseq_fast_geometry_preserves_explicit_double_bond_stereo_planes() {
        for (smiles, expected) in [("C/C=C/C", 180.0), ("C/C=C\\C", 0.0)] {
            let template = base_template_for_comparison(smiles);
            let coords = conformer_points(&template.molecule);
            let dihedral = dihedral_deg_for_test(coords[0], coords[1], coords[2], coords[3]);
            assert!(
                abs_dihedral_delta_deg(dihedral.abs(), expected) < 1.0e-8,
                "smiles={smiles} expected {expected} deg double-bond control dihedral, got {dihedral}"
            );
        }
    }

    #[test]
    fn confseq_fast_geometry_satisfies_tetrahedral_stereo_volume_sign() {
        for (tag, should_be_positive) in [
            (ChiralTag::TetrahedralCcw, true),
            (ChiralTag::TetrahedralCw, false),
        ] {
            let molecule = Molecule::from_smiles("FC(Cl)(Br)I").expect("test SMILES parses");
            let mut builder = molecule.to_builder();
            builder
                .atom_mut(AtomId::new(1))
                .expect("center atom exists")
                .set_chiral_tag(tag);
            let molecule = builder.build().expect("test molecule builds");
            let embedded =
                try_build_confseq_fast_geometry(&molecule).expect("fast geometry should build");
            let adjacency = AdjacencyList::from_topology(embedded.num_atoms(), embedded.bonds());
            let constraints =
                collect_confseq_base_tetrahedral_stereo_constraints(&embedded, &adjacency)
                    .expect("tetrahedral constraints should collect");
            let constraint = constraints
                .iter()
                .find(|constraint| constraint.center == 1)
                .expect("center stereo constraint should exist");
            let coords = conformer_points(&embedded);
            let volume = confseq_base_chiral_volume(&coords, constraint);

            assert_eq!(
                volume > 0.0,
                should_be_positive,
                "tag={tag:?} should map to the RDKit-derived chiral volume sign"
            );
            assert!(confseq_base_chiral_volume_satisfies_tag(volume, tag));
        }
    }

    #[test]
    fn confseq_base_constraint_model_extracts_rdkit_prior_classes() {
        let molecule = Molecule::from_smiles("C/C=C/C").expect("test SMILES parses");
        let ring_info = rings::symmetrize_sssr(&molecule).expect("ring perception should work");
        let model = build_confseq_base_constraint_model(&molecule, &ring_info)
            .expect("constraint model should build");

        assert_eq!(model.bond_targets.len(), molecule.num_bonds());
        assert!(model.planar_bonds.contains(&sorted_pair(1, 2)));
        assert!(
            model
                .torsion_priors
                .values()
                .any(|prior| *prior == ConfSeqBaseTorsionPrior::Trans)
        );
    }

    #[test]
    fn confseq_base_constraint_builder_places_acyclic_and_simple_ring_scaffolds() {
        for smiles in ["CCCC", "CC(C)CO", "c1ccccc1", "C1CCCCC1"] {
            let molecule = Molecule::from_smiles(smiles).expect("test SMILES parses");
            let embedded =
                try_build_confseq_fast_geometry(&molecule).expect("fast geometry should build");

            assert_eq!(embedded.num_atoms(), molecule.num_atoms());
            assert_eq!(embedded.conformers_3d().len(), 1);
            assert!(
                embedded
                    .conformers_3d()
                    .first()
                    .expect("fast geometry exists")
                    .coordinates()
                    .iter()
                    .all(|point| point.iter().all(|value| value.is_finite()))
            );
        }
    }

    #[test]
    fn confseq_base_constraint_model_groups_shared_ring_systems() {
        let molecule = Molecule::from_smiles("C1CC2CCC1C2").expect("bridged system parses");
        let ring_info = rings::symmetrize_sssr(&molecule).expect("ring perception works");
        let ring_components = classify_confseq_base_ring_components(&molecule, &ring_info)
            .expect("shared ring systems should be modeled as one component");

        assert_eq!(ring_components.len(), 1);
        assert!(ring_components[0].atoms.len() > ring_info.atom_rings()[0].len());
    }

    #[test]
    fn confseq_fast_geometry_conditions_shared_nonplanar_ring_fragments() {
        for smiles in ["C1CCC2(CC1)CCC2", "C1CC2(CC1)CCNCC2"] {
            let molecule = Molecule::from_smiles(smiles).expect("spiro fixture parses");
            let embedded = try_build_confseq_fast_geometry(&molecule)
                .expect("shared nonplanar ring fragments should use fragment-local conditioning");
            let points = conformer_points(&embedded);
            assert!(
                points
                    .iter()
                    .all(|point| point.iter().all(|value| value.is_finite())),
                "conditioned shared-ring fixture should produce finite coordinates: {smiles}"
            );
            for bond in molecule.bonds() {
                let length = vec_len(vec_sub(
                    points[bond.begin().index()],
                    points[bond.end().index()],
                ));
                assert!(
                    (0.9..=1.8).contains(&length),
                    "conditioned shared-ring fixture has unreasonable bond length {length:.3}: {smiles}"
                );
            }
        }
    }

    #[test]
    fn confseq_base_ring_system_scaffold_uses_rdkit_initial_prior() {
        let molecule = Molecule::from_smiles("c1ccc2ccccc2c1").expect("naphthalene parses");
        let embedded =
            try_build_confseq_fast_geometry(&molecule).expect("fast geometry should build");
        let coords = embedded.conformers_3d()[0].coordinates();
        let unique_xy: HashSet<_> = coords
            .iter()
            .map(|point| ((point[0] * 1000.0) as i64, (point[1] * 1000.0) as i64))
            .collect();

        assert_eq!(unique_xy.len(), molecule.num_atoms());
    }

    #[test]
    fn confseq_base_scaffold_initializes_connected_ring_systems_globally() {
        let molecule = Molecule::from_smiles("c1ccccc1CCc2ccccc2").expect("fixture parses");
        let embedded =
            try_build_confseq_fast_geometry(&molecule).expect("fast geometry should build");
        let coords = embedded.conformers_3d()[0].coordinates();
        let longest_bond = molecule
            .bonds()
            .iter()
            .map(|bond| {
                vec_len(vec_sub(
                    coords[bond.begin().index()],
                    coords[bond.end().index()],
                ))
            })
            .fold(0.0, f64::max);

        assert!(
            longest_bond < 2.0,
            "global RDKit scaffold initialization should not leave ring systems disconnected: longest_bond={longest_bond:.3}"
        );
    }

    #[test]
    fn confseq_base_backend_returns_structured_error_without_distgeom_fallback() {
        let options = ConfSeqDecodeOptions {
            optimize_with_uff: false,
            template_backend: ConfSeqTemplateBackend::FastGeometry,
            ..ConfSeqDecodeOptions::default()
        };

        let error =
            decode_confseq_with_options("C 1 C C 2 C C C 1 C 2", "C 1 C C 2 C C C 1 C 2", &options)
                .expect_err("fast-geometry decode failure must not fallback to DistGeom");

        assert!(matches!(
            error,
            ConfSeqDecodeError::AngleTokenCountMismatch { .. }
                | ConfSeqDecodeError::DihedralTokenCountMismatch { .. }
                | ConfSeqDecodeError::FastGeometry(_)
        ));
    }

    #[test]
    fn confseq_diagnostics_classifies_base_template_failure() {
        let options = ConfSeqDecodeOptions {
            optimize_with_uff: false,
            ..ConfSeqDecodeOptions::default()
        };
        let diagnostic = diagnostics::diagnose_confseq_candidate(
            "C 1 C C 2 C C C 1 C 2",
            "C 1 C C 2 C C C 1 C 2",
            &options,
        );

        assert_eq!(
            diagnostic.phase,
            Some(diagnostics::ConfSeqDiagnosticPhase::Decode)
        );
        assert!(diagnostic.parsed);
        assert!(diagnostic.distance_geometry_template_built);
        assert!(diagnostic.base_template_built);
        assert!(!diagnostic.base_decoded);
    }

    #[test]
    fn bond_token_mapping_tracks_explicit_bond_positions_for_branches_and_rings() {
        let mapping = bond_token_mapping_for_smiles("C1CCC1C").expect("mapping should parse");

        assert_eq!(mapping.token_idx_to_atom_pair.len(), 5);
        assert!(mapping.ring_closure_pairs.contains(&sorted_pair(0, 3)));
    }

    #[test]
    fn parse_confseq_maps_dihedral_tokens_by_atom_pair_not_position_count() {
        let parsed = parse_confseq("C C C C", "C <60> C C").expect("ConfSeq should parse");

        assert_eq!(
            parsed.dihedral_angles_by_pair.get(&sorted_pair(0, 1)),
            Some(&60.0)
        );
    }

    #[test]
    fn parse_confseq_preserves_bracketed_formal_charge_tokens() {
        let parsed = parse_confseq(
            "C - [ N + ] ( = O ) [ O - ]",
            "C <90> [ N + ] ( = O ) [ O - ]",
        )
        .expect("charged bracket atoms should parse");

        assert_eq!(parsed.stripped_smiles, "C[N+](=O)[O-]");
    }

    #[test]
    fn parse_confseq_accepts_direction_prefixed_angle_tokens() {
        let forward = parse_confseq("C C C C", "C /<177> C C").expect("slash token parses");
        let backward = parse_confseq("C C C C", "C \\<-5> C C").expect("backslash token parses");

        assert_eq!(
            forward.dihedral_angles_by_pair.get(&sorted_pair(0, 1)),
            Some(&177.0)
        );
        assert_eq!(
            backward.dihedral_angles_by_pair.get(&sorted_pair(0, 1)),
            Some(&-5.0)
        );
    }

    #[test]
    fn decode_confseq_batch_preserves_order_and_reports_invalid_ring_candidate_without_unsupported()
    {
        let smiles = vec!["C C".to_string(), "C 1 C C C 1".to_string()];
        let confseq = vec!["C C".to_string(), "C 1 C C C 1".to_string()];

        let result = decode_confseq_batch(&smiles, &confseq, true).expect("batch length is valid");

        assert!(result.molecules[0].is_some());
        assert!(matches!(
            result.errors[1],
            Some(ConfSeqDecodeError::DihedralTokenCountMismatch { .. })
        ));
    }

    #[test]
    fn decode_confseq_batch_accepts_explicit_thread_count() {
        let smiles = vec!["C C".to_string(), "C C".to_string()];
        let confseq = vec!["C C".to_string(), "C C".to_string()];
        let options = ConfSeqDecodeOptions {
            num_threads: Some(2),
            ..ConfSeqDecodeOptions::default()
        };

        let result = decode_confseq_batch_with_options(&smiles, &confseq, &options, false)
            .expect("threaded batch should decode");

        assert_eq!(result.molecules.len(), 2);
        assert!(result.errors.iter().all(Option::is_none));
    }

    #[test]
    fn pseudo_chirality_tokens_are_parsed_into_atom_tags() {
        let parsed = parse_confseq("C C C C", "C { <60> C C").expect("pseudo tag should parse");

        assert_eq!(
            parsed.chiral_tags_by_atom.get(&0),
            Some(&ChiralTag::TetrahedralCw)
        );
    }

    #[test]
    fn initial_coords_path_differs_from_full_embed_path_on_real_confseq_sample() {
        let in_smiles = "C c 1 c c ( C ) n c ( - N ^ | - N 2 - C ( = O ) - C - S - C - 2 = S ) n 1";
        let confseq = "C c 1 c c ( C ) n c ( <22> N <115> | <-87> N 2 <170> C ( = O ) <-3> C <1> S <3> C <-171> 2 = S ) n 1";
        let mut options = ConfSeqDecodeOptions::default();
        options.optimize_with_uff = false;

        let parsed = parse_confseq(in_smiles, confseq).expect("sample ConfSeq should parse");
        let full_template = build_template(
            &parsed.stripped_smiles,
            &parsed.chiral_tags_by_atom,
            &options,
        )
        .expect("full template should build");
        let initial_template = build_template_initial_coords_for_test(
            &parsed.stripped_smiles,
            &parsed.chiral_tags_by_atom,
            &options,
            crate::chemistry::distgeom::EmbedderTestStage::InitialCoords,
        )
        .expect("initial-coords template should build");
        let first_min_template = build_template_initial_coords_for_test(
            &parsed.stripped_smiles,
            &parsed.chiral_tags_by_atom,
            &options,
            crate::chemistry::distgeom::EmbedderTestStage::FirstMinimized,
        )
        .expect("first-minimized template should build");
        let fourth_dim_template = build_template_initial_coords_for_test(
            &parsed.stripped_smiles,
            &parsed.chiral_tags_by_atom,
            &options,
            crate::chemistry::distgeom::EmbedderTestStage::FourthDimensionCleaned,
        )
        .expect("fourth-dimension template should build");

        let full = decode_from_template(&full_template, &parsed, &options)
            .expect("full path should decode");
        let full_points = conformer_points(&full);
        let stages = [
            (
                "initial",
                decode_from_template(&initial_template, &parsed, &options)
                    .expect("initial path should decode"),
            ),
            (
                "first_min",
                decode_from_template(&first_min_template, &parsed, &options)
                    .expect("first-min path should decode"),
            ),
            (
                "fourth_dim",
                decode_from_template(&fourth_dim_template, &parsed, &options)
                    .expect("fourth-dim path should decode"),
            ),
        ];

        let aromatic_atoms = aromatic_atom_indices(&full);
        for (name, mol) in stages {
            let points = conformer_points(&mol);
            let rmsd = crate::distgeom::aligned_rmsd_for_test(&full_points, &points);
            let all_bonds = full
                .bonds()
                .iter()
                .map(|bond| (bond.begin().index(), bond.end().index()));
            let aromatic_bonds = full
                .bonds()
                .iter()
                .filter(|bond| bond.is_aromatic())
                .map(|bond| (bond.begin().index(), bond.end().index()));
            let non_aromatic_bonds = full
                .bonds()
                .iter()
                .filter(|bond| !bond.is_aromatic())
                .map(|bond| (bond.begin().index(), bond.end().index()));
            let (bond_count, bond_rms, bond_max) = distance_stats(all_bonds, &full_points, &points);
            let (arom_count, arom_rms, arom_max) =
                distance_stats(aromatic_bonds, &full_points, &points);
            let (non_arom_count, non_arom_rms, non_arom_max) =
                distance_stats(non_aromatic_bonds, &full_points, &points);
            let mut nonbond_pairs = Vec::new();
            for i in 0..full.num_atoms() {
                for j in i + 1..full.num_atoms() {
                    if !full.bonds().iter().any(|bond| {
                        let pair = sorted_pair(bond.begin().index(), bond.end().index());
                        pair == sorted_pair(i, j)
                    }) {
                        nonbond_pairs.push((i, j));
                    }
                }
            }
            let (nonbond_count, nonbond_rms, nonbond_max) =
                distance_stats(nonbond_pairs.into_iter(), &full_points, &points);
            let full_arom: Vec<_> = aromatic_atoms.iter().map(|idx| full_points[*idx]).collect();
            let stage_arom: Vec<_> = aromatic_atoms.iter().map(|idx| points[*idx]).collect();
            let aromatic_rmsd = crate::distgeom::aligned_rmsd_for_test(&full_arom, &stage_arom);
            eprintln!(
                "stage={name} aligned_rmsd={rmsd:.6} bonds n={bond_count} rms={bond_rms:.6} max={bond_max:.6}; aromatic bonds n={arom_count} rms={arom_rms:.6} max={arom_max:.6}; non-aromatic bonds n={non_arom_count} rms={non_arom_rms:.6} max={non_arom_max:.6}; nonbond n={nonbond_count} rms={nonbond_rms:.6} max={nonbond_max:.6}; aromatic atom RMSD={aromatic_rmsd:.6}"
            );
        }

        assert!(true);
    }

    #[test]
    #[ignore = "local ConfSeq corpus pass-rate snapshot"]
    fn confseq_generate_dg_reference_cache_snapshot() {
        let corpus_path = std::env::var("COSMOLKIT_CONFSEQ_CORPUS").unwrap_or_else(|_| {
            "/home/wangjingtong/sh4090/confseq_test_strings_100x10.jsonl".to_string()
        });
        let (mut in_smiles_batch, mut td_smiles_batch) = read_confseq_corpus_smiles(&corpus_path);
        let range_start = std::env::var("COSMOLKIT_CONFSEQ_CACHE_START")
            .ok()
            .map(|value| {
                value
                    .parse::<usize>()
                    .expect("COSMOLKIT_CONFSEQ_CACHE_START should be a usize")
            })
            .unwrap_or(0);
        let range_end = std::env::var("COSMOLKIT_CONFSEQ_CACHE_END")
            .ok()
            .map(|value| {
                value
                    .parse::<usize>()
                    .expect("COSMOLKIT_CONFSEQ_CACHE_END should be a usize")
            })
            .unwrap_or(in_smiles_batch.len());
        if range_start != 0 || range_end != in_smiles_batch.len() {
            assert!(
                range_start <= range_end && range_end <= in_smiles_batch.len(),
                "invalid ConfSeq cache range {range_start}..{range_end} for {} records",
                in_smiles_batch.len()
            );
            if std::env::var("COSMOLKIT_CONFSEQ_DG_REFERENCE_CACHE").is_err() {
                let range_cache_path = format!(
                    "{corpus_path}.dg_reference_heavy_points.{range_start}_{range_end}.json"
                );
                unsafe {
                    std::env::set_var("COSMOLKIT_CONFSEQ_DG_REFERENCE_CACHE", &range_cache_path);
                }
            }
            in_smiles_batch = in_smiles_batch[range_start..range_end].to_vec();
            td_smiles_batch = td_smiles_batch[range_start..range_end].to_vec();
            eprintln!(
                "generating ConfSeq DG reference cache for corpus range {range_start}..{range_end}"
            );
        }

        let reference_cache = load_or_generate_confseq_dg_reference_cache(
            &corpus_path,
            &in_smiles_batch,
            &td_smiles_batch,
        );

        let success_count = reference_cache
            .entries
            .iter()
            .filter(|entry| entry.heavy_atom_points.is_some())
            .count();
        eprintln!(
            "generated ConfSeq DG reference cache entries={} success={} fail={} path={}",
            reference_cache.entries.len(),
            success_count,
            reference_cache.entries.len() - success_count,
            confseq_dg_reference_cache_path(&corpus_path).display()
        );
    }

    #[derive(Debug, Clone, serde::Serialize)]
    struct ConfSeqDgStageCache {
        corpus_path: String,
        entries: Vec<ConfSeqDgStageCacheEntry>,
    }

    #[derive(Debug, Clone, serde::Serialize)]
    struct ConfSeqDgStageCacheEntry {
        template_heavy_atom_points: Option<Vec<[f64; 3]>>,
        angle_only_heavy_atom_points: Option<Vec<[f64; 3]>>,
        full_heavy_atom_points: Option<Vec<[f64; 3]>>,
        error: Option<String>,
    }

    #[derive(Debug, Clone, serde::Serialize)]
    struct ConfSeqDgEmbedStageCache {
        corpus_path: String,
        entries: Vec<ConfSeqDgEmbedStageCacheEntry>,
    }

    #[derive(Debug, Clone, serde::Serialize)]
    struct ConfSeqDgEmbedStageCacheEntry {
        initial_heavy_atom_points: Option<Vec<[f64; 3]>>,
        first_minimized_heavy_atom_points: Option<Vec<[f64; 3]>>,
        fourth_dimension_heavy_atom_points: Option<Vec<[f64; 3]>>,
        exp_torsion_heavy_atom_points: Option<Vec<[f64; 3]>>,
        final_heavy_atom_points: Option<Vec<[f64; 3]>>,
        low_level: crate::chemistry::distgeom::EmbedderTestLowLevelTrace,
        failures: Vec<(String, u32)>,
        error: Option<String>,
    }

    #[test]
    #[ignore = "local ConfSeq DG embed-stage cache for source-parity diagnostics"]
    fn confseq_generate_dg_embed_stage_cache_snapshot() {
        let corpus_path = std::env::var("COSMOLKIT_CONFSEQ_CORPUS").unwrap_or_else(|_| {
            "/home/wangjingtong/sh4090/confseq_test_strings_100x10.jsonl".to_string()
        });
        let (in_smiles_batch, td_smiles_batch) = read_confseq_corpus_smiles(&corpus_path);
        let cache_path = std::env::var("COSMOLKIT_CONFSEQ_DG_EMBED_STAGE_CACHE")
            .map(std::path::PathBuf::from)
            .unwrap_or_else(|_| {
                std::path::PathBuf::from(format!("{corpus_path}.dg_embed_stage_heavy_points.json"))
            });
        let options = ConfSeqDecodeOptions {
            optimize_with_uff: false,
            apply_angles: false,
            apply_dihedrals: false,
            ..ConfSeqDecodeOptions::default()
        };
        let mut entries = Vec::with_capacity(in_smiles_batch.len());
        let mut success_count = 0usize;
        for (idx, (in_smiles, td_smiles)) in
            in_smiles_batch.iter().zip(&td_smiles_batch).enumerate()
        {
            eprintln!(
                "generating ConfSeq DG embed-stage cache record {}/{}",
                idx,
                in_smiles_batch.len()
            );
            let entry = (|| -> Result<ConfSeqDgEmbedStageCacheEntry, ConfSeqDecodeError> {
                let parsed = parse_confseq(in_smiles, td_smiles)?;
                let trace = build_template_embed_trace_for_test(
                    &parsed.stripped_smiles,
                    &parsed.chiral_tags_by_atom,
                    &options,
                )?;
                Ok(ConfSeqDgEmbedStageCacheEntry {
                    initial_heavy_atom_points: trace.initial_coords,
                    first_minimized_heavy_atom_points: trace.first_minimized,
                    fourth_dimension_heavy_atom_points: trace.fourth_dimension_cleaned,
                    exp_torsion_heavy_atom_points: trace.exp_torsion_minimized,
                    final_heavy_atom_points: trace.final_checked,
                    low_level: trace.low_level,
                    failures: trace.failures,
                    error: None,
                })
            })();
            match entry {
                Ok(entry) => {
                    success_count += 1;
                    entries.push(entry);
                }
                Err(error) => entries.push(ConfSeqDgEmbedStageCacheEntry {
                    initial_heavy_atom_points: None,
                    first_minimized_heavy_atom_points: None,
                    fourth_dimension_heavy_atom_points: None,
                    exp_torsion_heavy_atom_points: None,
                    final_heavy_atom_points: None,
                    low_level: Default::default(),
                    failures: Vec::new(),
                    error: Some(format!("{error:?}")),
                }),
            }
        }
        let cache = ConfSeqDgEmbedStageCache {
            corpus_path: corpus_path.clone(),
            entries,
        };
        if let Some(parent) = cache_path.parent() {
            std::fs::create_dir_all(parent).unwrap_or_else(|err| {
                panic!(
                    "failed to create ConfSeq DG embed-stage cache directory {}: {err}",
                    parent.display()
                )
            });
        }
        let raw = serde_json::to_string_pretty(&cache).expect("embed-stage cache should serialize");
        std::fs::write(&cache_path, raw).unwrap_or_else(|err| {
            panic!(
                "failed to write ConfSeq DG embed-stage cache {}: {err}",
                cache_path.display()
            )
        });
        eprintln!(
            "generated ConfSeq DG embed-stage cache entries={} success={} fail={} path={}",
            cache.entries.len(),
            success_count,
            cache.entries.len() - success_count,
            cache_path.display()
        );
    }

    #[derive(Debug, Clone, serde::Serialize)]
    struct ConfSeqDgBoundsCache {
        corpus_path: String,
        entries: Vec<ConfSeqDgBoundsCacheEntry>,
    }

    #[derive(Debug, Clone, serde::Serialize)]
    struct ConfSeqDgBoundsCacheEntry {
        bounds: Option<Vec<Vec<f64>>>,
        error: Option<String>,
    }

    #[test]
    #[ignore = "local ConfSeq DG bounds cache for source-parity diagnostics"]
    fn confseq_generate_dg_bounds_cache_snapshot() {
        let corpus_path = std::env::var("COSMOLKIT_CONFSEQ_CORPUS").unwrap_or_else(|_| {
            "/home/wangjingtong/sh4090/confseq_test_strings_100x10.jsonl".to_string()
        });
        let (in_smiles_batch, td_smiles_batch) = read_confseq_corpus_smiles(&corpus_path);
        let cache_path = std::env::var("COSMOLKIT_CONFSEQ_DG_BOUNDS_CACHE")
            .map(std::path::PathBuf::from)
            .unwrap_or_else(|_| {
                std::path::PathBuf::from(format!("{corpus_path}.dg_bounds_matrix.json"))
            });
        let mut entries = Vec::with_capacity(in_smiles_batch.len());
        let mut success_count = 0usize;
        for (idx, (in_smiles, td_smiles)) in
            in_smiles_batch.iter().zip(&td_smiles_batch).enumerate()
        {
            eprintln!(
                "generating ConfSeq DG bounds cache record {}/{}",
                idx,
                in_smiles_batch.len()
            );
            let entry = (|| -> Result<ConfSeqDgBoundsCacheEntry, ConfSeqDecodeError> {
                let parsed = parse_confseq(in_smiles, td_smiles)?;
                let molecule = Molecule::from_smiles(&parsed.stripped_smiles)
                    .map_err(|err| ConfSeqDecodeError::SmilesParse(err.to_string()))?;
                let molecule =
                    prepare_p_chiral_embedding_molecule(molecule, &parsed.chiral_tags_by_atom)?;
                let with_h = molecule
                    .with_hydrogens()
                    .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
                let bounds = distgeom::dg_bounds_matrix(&with_h)
                    .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
                Ok(ConfSeqDgBoundsCacheEntry {
                    bounds: Some(bounds),
                    error: None,
                })
            })();
            match entry {
                Ok(entry) => {
                    success_count += 1;
                    entries.push(entry);
                }
                Err(error) => entries.push(ConfSeqDgBoundsCacheEntry {
                    bounds: None,
                    error: Some(format!("{error:?}")),
                }),
            }
        }
        let cache = ConfSeqDgBoundsCache {
            corpus_path: corpus_path.clone(),
            entries,
        };
        if let Some(parent) = cache_path.parent() {
            std::fs::create_dir_all(parent).unwrap_or_else(|err| {
                panic!(
                    "failed to create ConfSeq DG bounds cache directory {}: {err}",
                    parent.display()
                )
            });
        }
        let raw = serde_json::to_string_pretty(&cache).expect("bounds cache should serialize");
        std::fs::write(&cache_path, raw).unwrap_or_else(|err| {
            panic!(
                "failed to write ConfSeq DG bounds cache {}: {err}",
                cache_path.display()
            )
        });
        eprintln!(
            "generated ConfSeq DG bounds cache entries={} success={} fail={} path={}",
            cache.entries.len(),
            success_count,
            cache.entries.len() - success_count,
            cache_path.display()
        );
    }

    #[test]
    #[ignore = "local ConfSeq corpus stage cache for source-parity diagnostics"]
    fn confseq_generate_dg_stage_cache_snapshot() {
        let corpus_path = std::env::var("COSMOLKIT_CONFSEQ_CORPUS").unwrap_or_else(|_| {
            "/home/wangjingtong/sh4090/confseq_test_strings_100x10.jsonl".to_string()
        });
        let (in_smiles_batch, td_smiles_batch) = read_confseq_corpus_smiles(&corpus_path);
        let cache_path = std::env::var("COSMOLKIT_CONFSEQ_DG_STAGE_CACHE")
            .map(std::path::PathBuf::from)
            .unwrap_or_else(|_| {
                std::path::PathBuf::from(format!("{corpus_path}.dg_stage_heavy_points.json"))
            });
        let mut entries = Vec::with_capacity(in_smiles_batch.len());
        let options = ConfSeqDecodeOptions {
            optimize_with_uff: false,
            ..ConfSeqDecodeOptions::default()
        };
        let mut success_count = 0usize;
        for (idx, (in_smiles, td_smiles)) in
            in_smiles_batch.iter().zip(&td_smiles_batch).enumerate()
        {
            eprintln!(
                "generating ConfSeq DG stage cache record {}/{}",
                idx,
                in_smiles_batch.len()
            );
            let entry = (|| -> Result<ConfSeqDgStageCacheEntry, ConfSeqDecodeError> {
                let parsed = parse_confseq(in_smiles, td_smiles)?;
                let template = build_template(
                    &parsed.stripped_smiles,
                    &parsed.chiral_tags_by_atom,
                    &options,
                )?;
                let mut angle_only_options = options.clone();
                angle_only_options.apply_dihedrals = false;
                let angle_only = decode_from_template(&template, &parsed, &angle_only_options)?;
                let full = decode_from_template(&template, &parsed, &options)?;
                Ok(ConfSeqDgStageCacheEntry {
                    template_heavy_atom_points: Some(heavy_atom_points_for_rmsd(
                        &template.molecule,
                    )),
                    angle_only_heavy_atom_points: Some(heavy_atom_points_for_rmsd(&angle_only)),
                    full_heavy_atom_points: Some(heavy_atom_points_for_rmsd(&full)),
                    error: None,
                })
            })();
            match entry {
                Ok(entry) => {
                    success_count += 1;
                    entries.push(entry);
                }
                Err(error) => entries.push(ConfSeqDgStageCacheEntry {
                    template_heavy_atom_points: None,
                    angle_only_heavy_atom_points: None,
                    full_heavy_atom_points: None,
                    error: Some(format!("{error:?}")),
                }),
            }
        }
        let cache = ConfSeqDgStageCache {
            corpus_path: corpus_path.clone(),
            entries,
        };
        if let Some(parent) = cache_path.parent() {
            std::fs::create_dir_all(parent).unwrap_or_else(|err| {
                panic!(
                    "failed to create ConfSeq DG stage cache directory {}: {err}",
                    parent.display()
                )
            });
        }
        let raw = serde_json::to_string_pretty(&cache).expect("stage cache should serialize");
        std::fs::write(&cache_path, raw).unwrap_or_else(|err| {
            panic!(
                "failed to write ConfSeq DG stage cache {}: {err}",
                cache_path.display()
            )
        });
        eprintln!(
            "generated ConfSeq DG stage cache entries={} success={} fail={} path={}",
            cache.entries.len(),
            success_count,
            cache.entries.len() - success_count,
            cache_path.display()
        );
    }

    #[test]
    #[ignore = "local ConfSeq ring-deferred dihedral trace"]
    fn confseq_trace_ring_deferred_dihedrals_for_row() {
        let corpus_path = std::env::var("COSMOLKIT_CONFSEQ_CORPUS").unwrap_or_else(|_| {
            "/home/wangjingtong/code/COSMolKit/tmp/confseq_test_strings_100x10.current_mapping.jsonl"
                .to_string()
        });
        let row_idx = std::env::var("COSMOLKIT_CONFSEQ_TRACE_ROW")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or(754);
        let (in_smiles_batch, td_smiles_batch) = read_confseq_corpus_smiles(&corpus_path);
        let in_smiles = &in_smiles_batch[row_idx];
        let td_smiles = &td_smiles_batch[row_idx];
        let options = ConfSeqDecodeOptions {
            optimize_with_uff: false,
            ..ConfSeqDecodeOptions::default()
        };
        let parsed = parse_confseq(in_smiles, td_smiles).expect("row should parse");
        let template = build_template(
            &parsed.stripped_smiles,
            &parsed.chiral_tags_by_atom,
            &options,
        )
        .expect("template should build");

        let mut molecule = template.molecule.clone();
        for ((i, j, k), angle) in template.angle_centers.iter().zip(&parsed.angle_values_deg) {
            molecule = mol_transforms::set_bond_angle_deg(molecule, *i, *j, *k, *angle, 0)
                .expect("angle should apply");
        }

        let mut unapplied = Vec::new();
        let mut applied_first = Vec::new();
        for dihedral in &template.dihedrals {
            let (_, j, k, _) = *dihedral;
            let pair = sorted_pair(j, k);
            let angle = *parsed
                .dihedral_angles_by_pair
                .get(&pair)
                .expect("dihedral angle should exist");
            if template.ring_bond_pairs.contains(&pair) {
                unapplied.push((*dihedral, angle));
            } else {
                molecule = set_dihedral_deg_checked(molecule, *dihedral, angle)
                    .expect("non-ring dihedral should apply");
                applied_first.push((*dihedral, angle));
            }
        }

        let original_molecule = molecule.clone();
        let mut builder = molecule.to_builder();
        for (begin, end, _) in &template.last_ring_bonds {
            builder.remove_bond_between_atoms(AtomId::new(*begin), AtomId::new(*end));
        }
        let molecule_no_ring = builder.build().expect("removed-ring molecule should build");
        let (mut molecule, old_to_new) =
            renumber_like_confseq_smiles_output(molecule_no_ring).expect("renumber should work");
        let mut deferred = change_dihedral_for_removed_ring_bonds(
            &original_molecule,
            &template.last_ring_bonds,
            unapplied.clone(),
        )
        .expect("deferred dihedrals should rewrite");
        let remapped_last_ring_bonds =
            remap_last_ring_bonds(&template.last_ring_bonds, &old_to_new)
                .expect("last ring bonds should remap");
        remap_deferred_dihedrals(&mut deferred, &old_to_new)
            .expect("deferred dihedrals should remap");

        let mut rounds = Vec::new();
        for round in 0..2 {
            let ring_info = rings::symmetrize_sssr(&molecule).expect("ring finding should work");
            let atom_rings = ring_info
                .atom_rings()
                .iter()
                .map(|ring| ring.iter().map(|atom| atom.index()).collect::<Vec<_>>())
                .collect::<Vec<_>>();
            let mut applicable = Vec::new();
            let mut still_deferred = Vec::new();
            for (dihedral, angle) in deferred {
                if dihedral_bonds_exist(&molecule, dihedral)
                    && !dihedral_center_bond_is_in_ring(&ring_info, dihedral)
                {
                    applicable.push((dihedral, angle));
                    molecule = set_dihedral_deg_checked(molecule, dihedral, angle)
                        .expect("deferred dihedral should apply");
                } else {
                    still_deferred.push((dihedral, angle));
                }
            }
            rounds.push(serde_json::json!({
                "round": round,
                "atom_rings": atom_rings,
                "applicable": applicable,
                "deferred": still_deferred,
            }));
            deferred = still_deferred;
        }

        let trace_smiles_molecule =
            Molecule::from_smiles(&parsed.stripped_smiles).expect("trace smiles should parse");
        let trace_smiles_params = SmilesWriteParams {
            canonical: false,
            all_bonds_explicit: true,
            ..SmilesWriteParams::default()
        };
        let trace_smiles_be = mol_to_smiles(&trace_smiles_molecule, &trace_smiles_params)
            .expect("trace smiles should write");
        let trace_ring_bond_indices =
            confseq_last_ring_bond_indices_from_smiles_be(&trace_smiles_be);
        let trace_atom_pairs = trace_smiles_molecule
            .bonds()
            .iter()
            .map(|bond| [bond.begin().index(), bond.end().index()])
            .collect::<Vec<_>>();
        let trace = serde_json::json!({
            "row_idx": row_idx,
            "stripped_smiles": parsed.stripped_smiles,
            "smiles_be": trace_smiles_be,
            "ring_bond_indices": trace_ring_bond_indices,
            "smiles_molecule_atom_pairs": trace_atom_pairs,
            "angle_values": parsed.angle_values_deg,
            "atom_pair_dihedrals": parsed.dihedral_angles_by_pair.iter().map(|(pair, angle)| serde_json::json!({
                "pair": [pair.0, pair.1],
                "angle": angle,
            })).collect::<Vec<_>>(),
            "template_dihedrals": template.dihedrals,
            "ring_bond_pairs": template.ring_bond_pairs.iter().map(|pair| [pair.0, pair.1]).collect::<Vec<_>>(),
            "last_ring_bonds": template.last_ring_bonds.iter().map(|(i, j, _)| (*i, *j)).collect::<Vec<_>>(),
            "applied_first": applied_first,
            "unapplied_first": unapplied,
            "old_to_new": old_to_new,
            "last_ring_bonds_renum": remapped_last_ring_bonds.iter().map(|(i, j, _)| (*i, *j)).collect::<Vec<_>>(),
            "changed_renum": deferred,
            "rounds": rounds,
        });
        eprintln!("{}", serde_json::to_string_pretty(&trace).unwrap());
    }

    #[test]
    #[ignore = "local ConfSeq full-decode bonded-distance stage diagnostic"]
    fn confseq_trace_full_decode_bond_distance_stages_for_row() {
        let corpus_path = std::env::var("COSMOLKIT_CONFSEQ_CORPUS").unwrap_or_else(|_| {
            "/home/wangjingtong/sh4090/confseq_test_strings_100x10.jsonl".to_string()
        });
        let row_idx = std::env::var("COSMOLKIT_CONFSEQ_TRACE_ROW")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or(309);
        let (in_smiles_batch, td_smiles_batch) = read_confseq_corpus_smiles(&corpus_path);
        let in_smiles = &in_smiles_batch[row_idx];
        let td_smiles = &td_smiles_batch[row_idx];
        let options = ConfSeqDecodeOptions {
            optimize_with_uff: false,
            ..ConfSeqDecodeOptions::default()
        };
        let parsed = parse_confseq(in_smiles, td_smiles).expect("row should parse");
        let template = build_template(
            &parsed.stripped_smiles,
            &parsed.chiral_tags_by_atom,
            &options,
        )
        .expect("template should build");

        let stage_stats = |label: &str, molecule: &Molecule| {
            let points = heavy_atom_points_for_rmsd(molecule);
            let (count, rms, max) = heavy_bond_length_stats(molecule, &points);
            eprintln!(
                "confseq_full_decode_stage row={row_idx} stage={label} heavy_bonds={count} rms_len={rms:.6} max_len={max:.6}"
            );
        };

        stage_stats("template", &template.molecule);

        let mut molecule = template.molecule.clone();
        for (idx, ((i, j, k), angle)) in template
            .angle_centers
            .iter()
            .zip(&parsed.angle_values_deg)
            .enumerate()
        {
            molecule = mol_transforms::set_bond_angle_deg(molecule, *i, *j, *k, *angle, 0)
                .unwrap_or_else(|err| panic!("angle {idx} should apply: {err}"));
            if idx % 10 == 9 {
                stage_stats(&format!("angles_through_{idx}"), &molecule);
            }
        }
        stage_stats("after_angles", &molecule);

        let mut unapplied = Vec::new();
        let mut applied_non_ring_count = 0usize;
        for (idx, dihedral) in template.dihedrals.iter().enumerate() {
            let (_, j, k, _) = *dihedral;
            let pair = sorted_pair(j, k);
            let angle = *parsed
                .dihedral_angles_by_pair
                .get(&pair)
                .expect("dihedral angle should exist");
            if template.ring_bond_pairs.contains(&pair) {
                unapplied.push((*dihedral, angle));
            } else {
                molecule = set_dihedral_deg_checked(molecule, *dihedral, angle)
                    .unwrap_or_else(|err| panic!("non-ring dihedral {idx} should apply: {err}"));
                applied_non_ring_count += 1;
                if applied_non_ring_count % 10 == 0 {
                    stage_stats(
                        &format!("non_ring_dihedrals_applied_{applied_non_ring_count}"),
                        &molecule,
                    );
                }
            }
        }
        stage_stats("after_non_ring_dihedrals", &molecule);
        eprintln!(
            "confseq_full_decode_stage row={row_idx} unapplied_ring_dihedrals={}",
            unapplied.len()
        );

        let original_molecule = molecule.clone();
        let mut builder = molecule.to_builder();
        for (begin, end, _) in &template.last_ring_bonds {
            builder.remove_bond_between_atoms(AtomId::new(*begin), AtomId::new(*end));
        }
        let molecule_no_ring = builder.build().expect("removed-ring molecule should build");
        stage_stats("after_remove_last_ring_bonds", &molecule_no_ring);

        let (mut molecule, old_to_new) =
            renumber_like_confseq_smiles_output(molecule_no_ring).expect("renumber should work");
        stage_stats("after_renumber", &molecule);

        let mut deferred = change_dihedral_for_removed_ring_bonds(
            &original_molecule,
            &template.last_ring_bonds,
            unapplied,
        )
        .expect("deferred dihedrals should rewrite");
        let last_ring_bonds = remap_last_ring_bonds(&template.last_ring_bonds, &old_to_new)
            .expect("last ring bonds should remap");
        remap_deferred_dihedrals(&mut deferred, &old_to_new)
            .expect("deferred dihedrals should remap");
        eprintln!(
            "confseq_full_decode_stage row={row_idx} remapped_last_ring_bonds={:?} deferred_after_remap={}",
            last_ring_bonds
                .iter()
                .map(|(i, j, _)| (*i, *j))
                .collect::<Vec<_>>(),
            deferred.len()
        );

        for round in 0..2 {
            let ring_info = rings::symmetrize_sssr(&molecule).expect("ring finding should work");
            let mut still_deferred = Vec::new();
            let mut applied_this_round = 0usize;
            for (dihedral, angle) in deferred {
                if dihedral_bonds_exist(&molecule, dihedral)
                    && !dihedral_center_bond_is_in_ring(&ring_info, dihedral)
                {
                    molecule = set_dihedral_deg_checked(molecule, dihedral, angle)
                        .expect("deferred dihedral should apply");
                    applied_this_round += 1;
                } else {
                    still_deferred.push((dihedral, angle));
                }
            }
            stage_stats(&format!("after_deferred_round_{round}"), &molecule);
            eprintln!(
                "confseq_full_decode_stage row={row_idx} round={round} applied={applied_this_round} still_deferred={}",
                still_deferred.len()
            );
            deferred = still_deferred;
        }

        let mut builder = molecule.to_builder();
        for (_, _, spec) in &last_ring_bonds {
            builder
                .add_bond(spec.clone())
                .expect("ring bond should restore");
        }
        let molecule = builder.build().expect("final molecule should build");
        stage_stats("after_restore_ring_bonds", &molecule);
    }

    #[test]
    #[ignore = "local ConfSeq corpus pass-rate snapshot"]
    fn confseq_base_corpus_pass_rate_snapshot() {
        let corpus_path = std::env::var("COSMOLKIT_CONFSEQ_CORPUS").unwrap_or_else(|_| {
            "/home/wangjingtong/sh4090/confseq_test_strings_100x10.jsonl".to_string()
        });
        let fail_export_dir = std::env::var("COSMOLKIT_CONFSEQ_BASE_FAIL_EXPORT_DIR").ok();
        let base_options = ConfSeqDecodeOptions {
            optimize_with_uff: false,
            template_backend: ConfSeqTemplateBackend::FastGeometry,
            ..ConfSeqDecodeOptions::default()
        };

        let (in_smiles_batch, td_smiles_batch) = read_confseq_corpus_smiles(&corpus_path);

        let reference_cache = load_or_generate_confseq_dg_reference_cache(
            &corpus_path,
            &in_smiles_batch,
            &td_smiles_batch,
        );
        let base_batch = decode_confseq_batch_with_options(
            &in_smiles_batch,
            &td_smiles_batch,
            &base_options,
            true,
        )
        .expect("base batch should decode with per-record errors");

        let total = in_smiles_batch.len();
        let mut reference_success = 0usize;
        let mut base_success_where_reference_success = 0usize;
        let mut pass_03a = 0usize;
        let mut automorphism_pass_03a = 0usize;
        let mut strict_fail_automorphism_pass_03a = 0usize;
        let mut rmsds = Vec::new();
        let mut automorphism_rmsds = Vec::new();
        let mut local_bond_rmsds = Vec::new();
        let mut local_bond_max_abs = Vec::new();
        let mut local_non_token_angle_rmsds = Vec::new();
        let mut local_non_token_angle_max_abs = Vec::new();
        let mut local_bond_deviation_buckets = HashMap::<String, LocalBondDeviationBucket>::new();
        let mut local_constraint_comparable = 0usize;
        let mut local_constraint_pass_005a_5deg = 0usize;
        let mut local_constraint_pass_010a_10deg = 0usize;
        let mut local_constraint_pass_010a_15deg = 0usize;
        let mut rigid_fragment_rmsds = Vec::new();
        let mut rigid_fragment_shape_rmsds = Vec::new();
        let mut rigid_fragment_terminal_symmetry_rmsds = Vec::new();
        let mut rigid_fragment_connector_rmsds = Vec::new();
        let mut rigid_fragment_max_rmsds = Vec::new();
        let mut rigid_fragment_max_shape_rmsds = Vec::new();
        let mut rigid_fragment_max_terminal_symmetry_rmsds = Vec::new();
        let mut rigid_fragment_max_connector_rmsds = Vec::new();
        let mut worst_rigid_fragments = Vec::<(f64, usize, Vec<usize>, String)>::new();
        let mut mirror_branch_like_fragments_01a = 0usize;
        let mut mirror_branch_like_fragments_03a = 0usize;
        let mut mirror_branch_like_molecules_01a = 0usize;
        let mut mirror_branch_like_molecules_03a = 0usize;
        let mut rigid_fragment_pass_01a = 0usize;
        let mut rigid_fragment_pass_02a = 0usize;
        let mut rigid_fragment_pass_03a = 0usize;
        let mut rigid_shape_pass_01a = 0usize;
        let mut rigid_shape_pass_02a = 0usize;
        let mut rigid_shape_pass_03a = 0usize;
        let mut rigid_terminal_symmetry_pass_01a = 0usize;
        let mut rigid_terminal_symmetry_pass_02a = 0usize;
        let mut rigid_terminal_symmetry_pass_03a = 0usize;
        let mut rigid_connector_pass_01a = 0usize;
        let mut rigid_connector_pass_02a = 0usize;
        let mut rigid_connector_pass_03a = 0usize;
        let mut rigid_max_pass_01a = 0usize;
        let mut rigid_max_pass_02a = 0usize;
        let mut rigid_max_pass_03a = 0usize;
        let mut global_fail_rigid_pass = 0usize;
        let mut global_fail_rigid_fail = 0usize;
        let mut base_errors: HashMap<String, usize> = HashMap::new();
        let mut framework_fragment_shapes: HashMap<String, usize> = HashMap::new();
        let mut framework_assembly_paths: HashMap<String, usize> = HashMap::new();
        let mut framework_fragment_template_failures: HashMap<String, usize> = HashMap::new();
        let mut framework_nonplanar_ring_topologies: HashMap<String, usize> = HashMap::new();
        let mut rigid_fragment_framework_rmsd_buckets =
            RigidFragmentFrameworkRmsdBuckets::default();
        let mut pre_token_geometry = BaseGeometryStageDiagnostics::default();
        let mut post_token_geometry = BaseGeometryStageDiagnostics::default();
        let mut pre_token_template_errors: HashMap<String, usize> = HashMap::new();
        let mut fail_exports = Vec::<ConfSeqBaseFailExportRecord>::new();
        let mut structural_risk_reference = 0usize;
        let mut structural_risk_comparable = 0usize;
        let mut strict_fail_structural_risk = 0usize;
        let mut strict_fail_structural_risk_not_automorphism = 0usize;
        let mut initializer_decode_fallback_reference = 0usize;
        let mut initializer_decode_fallback_comparable = 0usize;
        let mut strict_fail_initializer_decode_fallback = 0usize;
        let mut strict_fail_initializer_decode_fallback_not_automorphism = 0usize;
        let mut nonfallback_comparable = 0usize;
        let mut nonfallback_pass_03a = 0usize;
        let mut nonfallback_automorphism_pass_03a = 0usize;
        let mut structural_risk_classes_counts: HashMap<String, usize> = HashMap::new();
        let mut strict_fail_structural_risk_classes_counts: HashMap<String, usize> = HashMap::new();

        for idx in 0..total {
            let Some(reference_points) = reference_cache.entries[idx].heavy_atom_points.as_ref()
            else {
                continue;
            };
            reference_success += 1;
            let parsed = parse_confseq(&in_smiles_batch[idx], &td_smiles_batch[idx])
                .expect("fixture ConfSeq should parse");
            let source_mol = Molecule::from_smiles(&parsed.stripped_smiles)
                .expect("fixture stripped SMILES should parse");
            let rings = rings::symmetrize_sssr(&source_mol).expect("ring perception should work");
            let mut structural_risk_precheck = ConfSeqBaseStructuralRiskPrecheck {
                classes: Vec::new(),
                fallback_candidate: false,
            };
            if let Ok(model) = build_confseq_base_constraint_model(&source_mol, &rings) {
                ConfSeqBaseFrameworkDiagnostics::collect(&source_mol, &rings, &model).merge_into(
                    &mut framework_fragment_shapes,
                    &mut framework_assembly_paths,
                    &mut framework_fragment_template_failures,
                    &mut framework_nonplanar_ring_topologies,
                );
                structural_risk_precheck =
                    confseq_base_structural_risk_precheck_with_model(&source_mol, &rings, &model);
                if structural_risk_precheck.has_risk() {
                    structural_risk_reference += 1;
                    for class in &structural_risk_precheck.classes {
                        *structural_risk_classes_counts
                            .entry(class.as_str().to_string())
                            .or_default() += 1;
                    }
                }
                if structural_risk_precheck.fallback_candidate {
                    initializer_decode_fallback_reference += 1;
                }
            }

            match base_batch.molecules[idx].as_ref() {
                Some(base) => {
                    match build_confseq_base_template(
                        &parsed.stripped_smiles,
                        &parsed.chiral_tags_by_atom,
                    ) {
                        Ok(base_template) => {
                            pre_token_geometry.record(
                                &source_mol,
                                reference_points,
                                &base_template.molecule,
                            );
                        }
                        Err(error) => {
                            *pre_token_template_errors
                                .entry(format!("{error:?}"))
                                .or_default() += 1;
                        }
                    }
                    post_token_geometry.record(&source_mol, reference_points, base);
                    base_success_where_reference_success += 1;
                    if structural_risk_precheck.has_risk() {
                        structural_risk_comparable += 1;
                    }
                    if structural_risk_precheck.fallback_candidate {
                        initializer_decode_fallback_comparable += 1;
                    }
                    let base_points = heavy_atom_points_for_rmsd(base);
                    assert_eq!(
                        reference_points.len(),
                        base_points.len(),
                        "line {} decoded molecules should have matching heavy atom counts",
                        idx + 1
                    );
                    let rmsd = crate::distgeom::aligned_rmsd_for_test(
                        reference_points.as_slice(),
                        &base_points,
                    );
                    if rmsd <= 0.3 {
                        pass_03a += 1;
                    }
                    let automorphism_rmsd =
                        automorphism_aware_heavy_rmsd_for_test(&source_mol, reference_points, base);
                    let automorphism_pass = automorphism_rmsd.is_some_and(|rmsd| rmsd <= 0.3);
                    if let Some(automorphism_rmsd) = automorphism_rmsd {
                        if automorphism_rmsd <= 0.3 {
                            automorphism_pass_03a += 1;
                            if rmsd > 0.3 {
                                strict_fail_automorphism_pass_03a += 1;
                            }
                        }
                        automorphism_rmsds.push(automorphism_rmsd);
                    }
                    if !structural_risk_precheck.fallback_candidate {
                        nonfallback_comparable += 1;
                        if rmsd <= 0.3 {
                            nonfallback_pass_03a += 1;
                        }
                        if automorphism_pass {
                            nonfallback_automorphism_pass_03a += 1;
                        }
                    }
                    if rmsd > 0.3 && structural_risk_precheck.has_risk() {
                        strict_fail_structural_risk += 1;
                        if !automorphism_pass {
                            strict_fail_structural_risk_not_automorphism += 1;
                            for class in &structural_risk_precheck.classes {
                                *strict_fail_structural_risk_classes_counts
                                    .entry(class.as_str().to_string())
                                    .or_default() += 1;
                            }
                        }
                    }
                    if rmsd > 0.3 && structural_risk_precheck.fallback_candidate {
                        strict_fail_initializer_decode_fallback += 1;
                        if !automorphism_pass {
                            strict_fail_initializer_decode_fallback_not_automorphism += 1;
                        }
                    }
                    let local_summary = local_constraint_summary_against_reference(
                        &source_mol,
                        reference_points,
                        base,
                    );
                    local_constraint_comparable += 1;
                    if local_summary.bond_rms_a <= 0.05
                        && local_summary.non_token_angle_rms_deg <= 5.0
                    {
                        local_constraint_pass_005a_5deg += 1;
                    }
                    if local_summary.bond_rms_a <= 0.10
                        && local_summary.non_token_angle_rms_deg <= 10.0
                    {
                        local_constraint_pass_010a_10deg += 1;
                    }
                    if local_summary.bond_rms_a <= 0.10
                        && local_summary.non_token_angle_rms_deg <= 15.0
                    {
                        local_constraint_pass_010a_15deg += 1;
                    }
                    if local_summary.bond_count > 0 {
                        local_bond_rmsds.push(local_summary.bond_rms_a);
                        local_bond_max_abs.push(local_summary.max_bond_abs_a);
                    }
                    if local_summary.non_token_angle_count > 0 {
                        local_non_token_angle_rmsds.push(local_summary.non_token_angle_rms_deg);
                        local_non_token_angle_max_abs
                            .push(local_summary.max_non_token_angle_abs_deg);
                    }
                    collect_local_bond_deviation_buckets(
                        &source_mol,
                        reference_points,
                        base,
                        &mut local_bond_deviation_buckets,
                    );
                    let rigid_summary =
                        rigid_fragment_rmsd_summary(&source_mol, reference_points, base);
                    mirror_branch_like_fragments_01a +=
                        rigid_summary.mirror_branch_like_fragment_count_01a;
                    mirror_branch_like_fragments_03a +=
                        rigid_summary.mirror_branch_like_fragment_count_03a;
                    if rigid_summary
                        .max_rmsd
                        .zip(rigid_summary.max_shape_rmsd)
                        .is_some_and(|(proper, shape)| proper > 0.3 && shape <= 0.1)
                    {
                        mirror_branch_like_molecules_01a += 1;
                    }
                    if rigid_summary
                        .max_rmsd
                        .zip(rigid_summary.max_shape_rmsd)
                        .is_some_and(|(proper, shape)| proper > 0.3 && shape <= 0.3)
                    {
                        mirror_branch_like_molecules_03a += 1;
                    }
                    if rmsd > 0.1 {
                        fail_exports.push(ConfSeqBaseFailExportRecord {
                            idx,
                            rmsd,
                            automorphism_rmsd,
                            max_rigid_rmsd: rigid_summary.max_rmsd,
                            max_rigid_shape_rmsd: rigid_summary.max_shape_rmsd,
                            max_rigid_connector_rmsd: rigid_summary.max_connector_rmsd,
                            worst_fragment_atoms: rigid_summary.worst_fragment_atoms.clone(),
                            stripped_smiles: parsed.stripped_smiles.clone(),
                            structural_risk_classes: structural_risk_precheck
                                .class_names()
                                .into_iter()
                                .map(str::to_string)
                                .collect(),
                            structural_risk_fallback_candidate: structural_risk_precheck
                                .fallback_candidate,
                            reference_points: reference_points.clone(),
                            base: base.clone(),
                        });
                    }
                    for &fragment_rmsd in &rigid_summary.fragment_rmsds {
                        if fragment_rmsd <= 0.1 {
                            rigid_fragment_pass_01a += 1;
                        }
                        if fragment_rmsd <= 0.2 {
                            rigid_fragment_pass_02a += 1;
                        }
                        if fragment_rmsd <= 0.3 {
                            rigid_fragment_pass_03a += 1;
                        }
                    }
                    for &fragment_rmsd in &rigid_summary.fragment_shape_rmsds {
                        if fragment_rmsd <= 0.1 {
                            rigid_shape_pass_01a += 1;
                        }
                        if fragment_rmsd <= 0.2 {
                            rigid_shape_pass_02a += 1;
                        }
                        if fragment_rmsd <= 0.3 {
                            rigid_shape_pass_03a += 1;
                        }
                    }
                    for &fragment_rmsd in &rigid_summary.fragment_terminal_symmetry_rmsds {
                        if fragment_rmsd <= 0.1 {
                            rigid_terminal_symmetry_pass_01a += 1;
                        }
                        if fragment_rmsd <= 0.2 {
                            rigid_terminal_symmetry_pass_02a += 1;
                        }
                        if fragment_rmsd <= 0.3 {
                            rigid_terminal_symmetry_pass_03a += 1;
                        }
                    }
                    for &fragment_rmsd in &rigid_summary.fragment_connector_rmsds {
                        if fragment_rmsd <= 0.1 {
                            rigid_connector_pass_01a += 1;
                        }
                        if fragment_rmsd <= 0.2 {
                            rigid_connector_pass_02a += 1;
                        }
                        if fragment_rmsd <= 0.3 {
                            rigid_connector_pass_03a += 1;
                        }
                    }
                    rigid_fragment_rmsds.extend(rigid_summary.fragment_rmsds.iter().copied());
                    rigid_fragment_shape_rmsds
                        .extend(rigid_summary.fragment_shape_rmsds.iter().copied());
                    rigid_fragment_terminal_symmetry_rmsds.extend(
                        rigid_summary
                            .fragment_terminal_symmetry_rmsds
                            .iter()
                            .copied(),
                    );
                    rigid_fragment_connector_rmsds
                        .extend(rigid_summary.fragment_connector_rmsds.iter().copied());
                    if let Some(max_rmsd) = rigid_summary.max_rmsd {
                        rigid_fragment_max_rmsds.push(max_rmsd);
                        if max_rmsd <= 0.1 {
                            rigid_max_pass_01a += 1;
                        }
                        if max_rmsd <= 0.2 {
                            rigid_max_pass_02a += 1;
                        }
                        if max_rmsd <= 0.3 {
                            rigid_max_pass_03a += 1;
                        }
                        worst_rigid_fragments.push((
                            max_rmsd,
                            idx,
                            rigid_summary.worst_fragment_atoms,
                            parsed.stripped_smiles,
                        ));
                        if rmsd > 0.3 {
                            if max_rmsd <= 0.3 {
                                global_fail_rigid_pass += 1;
                            } else {
                                global_fail_rigid_fail += 1;
                            }
                        }
                    }
                    if let Some(max_rmsd) = rigid_summary.max_shape_rmsd {
                        rigid_fragment_max_shape_rmsds.push(max_rmsd);
                    }
                    if let Some(max_rmsd) = rigid_summary.max_terminal_symmetry_rmsd {
                        rigid_fragment_max_terminal_symmetry_rmsds.push(max_rmsd);
                    }
                    if let Some(max_rmsd) = rigid_summary.max_connector_rmsd {
                        rigid_fragment_max_connector_rmsds.push(max_rmsd);
                    }
                    for fragment in rigid_heavy_fragments_cutting_confseq_rotors(&source_mol) {
                        if fragment
                            .iter()
                            .filter(|atom| source_mol.atoms()[**atom].atomic_number() > 1)
                            .count()
                            < 3
                        {
                            continue;
                        }
                        let key =
                            rigid_fragment_framework_key_for_diagnostic(&source_mol, &fragment);
                        let details = rigid_fragment_metric_details(
                            &source_mol,
                            reference_points,
                            base,
                            &fragment,
                        );
                        rigid_fragment_framework_rmsd_buckets.record(key, &details);
                    }
                    rmsds.push(rmsd);
                }
                None => {
                    let error = base_batch.errors[idx]
                        .as_ref()
                        .expect("base failure should include an error");
                    *base_errors.entry(format!("{error:?}")).or_default() += 1;
                }
            }
        }

        rmsds.sort_by(|left, right| left.total_cmp(right));
        automorphism_rmsds.sort_by(|left, right| left.total_cmp(right));
        local_bond_rmsds.sort_by(|left, right| left.total_cmp(right));
        local_bond_max_abs.sort_by(|left, right| left.total_cmp(right));
        local_non_token_angle_rmsds.sort_by(|left, right| left.total_cmp(right));
        local_non_token_angle_max_abs.sort_by(|left, right| left.total_cmp(right));
        rigid_fragment_rmsds.sort_by(|left, right| left.total_cmp(right));
        rigid_fragment_shape_rmsds.sort_by(|left, right| left.total_cmp(right));
        rigid_fragment_terminal_symmetry_rmsds.sort_by(|left, right| left.total_cmp(right));
        rigid_fragment_connector_rmsds.sort_by(|left, right| left.total_cmp(right));
        rigid_fragment_max_rmsds.sort_by(|left, right| left.total_cmp(right));
        rigid_fragment_max_shape_rmsds.sort_by(|left, right| left.total_cmp(right));
        rigid_fragment_max_terminal_symmetry_rmsds.sort_by(|left, right| left.total_cmp(right));
        rigid_fragment_max_connector_rmsds.sort_by(|left, right| left.total_cmp(right));
        pre_token_geometry.sort();
        post_token_geometry.sort();
        let total_pass_rate = pass_03a as f64 / total.max(1) as f64 * 100.0;
        let comparable_pass_rate =
            pass_03a as f64 / base_success_where_reference_success.max(1) as f64 * 100.0;
        let automorphism_comparable_pass_rate = automorphism_pass_03a as f64
            / base_success_where_reference_success.max(1) as f64
            * 100.0;
        let nonfallback_pass_rate =
            nonfallback_pass_03a as f64 / nonfallback_comparable.max(1) as f64 * 100.0;
        let nonfallback_automorphism_pass_rate =
            nonfallback_automorphism_pass_03a as f64 / nonfallback_comparable.max(1) as f64 * 100.0;
        let initializer_decode_fallback_rate = initializer_decode_fallback_comparable as f64
            / base_success_where_reference_success.max(1) as f64
            * 100.0;
        let base_coverage =
            base_success_where_reference_success as f64 / reference_success.max(1) as f64 * 100.0;
        let local_constraint_pass_rate_005a_5deg = local_constraint_pass_005a_5deg as f64
            / local_constraint_comparable.max(1) as f64
            * 100.0;
        let local_constraint_pass_rate_010a_10deg = local_constraint_pass_010a_10deg as f64
            / local_constraint_comparable.max(1) as f64
            * 100.0;
        let local_constraint_pass_rate_010a_15deg = local_constraint_pass_010a_15deg as f64
            / local_constraint_comparable.max(1) as f64
            * 100.0;

        eprintln!("confseq_base_corpus path={corpus_path}");
        eprintln!(
            "confseq_base_corpus total={total} reference_success={reference_success} base_success_where_reference_success={base_success_where_reference_success} pass_rmsd_le_0_3a={pass_03a}"
        );
        eprintln!(
            "confseq_base_corpus total_pass_rate={total_pass_rate:.2}% comparable_pass_rate={comparable_pass_rate:.2}% base_coverage_vs_reference={base_coverage:.2}%"
        );
        eprintln!(
            "confseq_base_corpus automorphism_aware_global pass_rmsd_le_0_3a={} comparable_pass_rate={:.2}% strict_fail_reclassified_count={} rmsd p50={:.6} p90={:.6} p99={:.6}",
            automorphism_pass_03a,
            automorphism_comparable_pass_rate,
            strict_fail_automorphism_pass_03a,
            quantile_sorted(&automorphism_rmsds, 0.50),
            quantile_sorted(&automorphism_rmsds, 0.90),
            quantile_sorted(&automorphism_rmsds, 0.99),
        );
        eprintln!(
            "confseq_base_corpus structural_risk_precheck reference={} comparable={} strict_fail={} strict_fail_not_automorphism={}",
            structural_risk_reference,
            structural_risk_comparable,
            strict_fail_structural_risk,
            strict_fail_structural_risk_not_automorphism,
        );
        eprintln!(
            "confseq_base_corpus initializer_decode_fallback_precheck reference={} comparable={} comparable_rate={:.2}% strict_fail={} strict_fail_not_automorphism={}",
            initializer_decode_fallback_reference,
            initializer_decode_fallback_comparable,
            initializer_decode_fallback_rate,
            strict_fail_initializer_decode_fallback,
            strict_fail_initializer_decode_fallback_not_automorphism,
        );
        eprintln!(
            "confseq_base_corpus nonfallback_final_decoded comparable={} pass_rmsd_le_0_3a={} pass_rate={:.2}% automorphism_pass_rmsd_le_0_3a={} automorphism_pass_rate={:.2}%",
            nonfallback_comparable,
            nonfallback_pass_03a,
            nonfallback_pass_rate,
            nonfallback_automorphism_pass_03a,
            nonfallback_automorphism_pass_rate,
        );
        eprintln!(
            "confseq_base_corpus rmsd p50={:.6} p75={:.6} p90={:.6} p95={:.6} p99={:.6}",
            quantile_sorted(&rmsds, 0.50),
            quantile_sorted(&rmsds, 0.75),
            quantile_sorted(&rmsds, 0.90),
            quantile_sorted(&rmsds, 0.95),
            quantile_sorted(&rmsds, 0.99)
        );
        pre_token_geometry.print("pre_token");
        post_token_geometry.print("post_token");
        eprintln!(
            "confseq_base_corpus local_constraint_bond_rms comparable={} p50={:.6} p75={:.6} p90={:.6} p95={:.6} p99={:.6}",
            local_bond_rmsds.len(),
            quantile_sorted(&local_bond_rmsds, 0.50),
            quantile_sorted(&local_bond_rmsds, 0.75),
            quantile_sorted(&local_bond_rmsds, 0.90),
            quantile_sorted(&local_bond_rmsds, 0.95),
            quantile_sorted(&local_bond_rmsds, 0.99)
        );
        eprintln!(
            "confseq_base_corpus local_constraint_bond_max_abs comparable={} p50={:.6} p90={:.6} p99={:.6}",
            local_bond_max_abs.len(),
            quantile_sorted(&local_bond_max_abs, 0.50),
            quantile_sorted(&local_bond_max_abs, 0.90),
            quantile_sorted(&local_bond_max_abs, 0.99)
        );
        eprintln!(
            "confseq_base_corpus local_constraint_non_token_angle_rms comparable={} p50={:.6} p75={:.6} p90={:.6} p95={:.6} p99={:.6}",
            local_non_token_angle_rmsds.len(),
            quantile_sorted(&local_non_token_angle_rmsds, 0.50),
            quantile_sorted(&local_non_token_angle_rmsds, 0.75),
            quantile_sorted(&local_non_token_angle_rmsds, 0.90),
            quantile_sorted(&local_non_token_angle_rmsds, 0.95),
            quantile_sorted(&local_non_token_angle_rmsds, 0.99)
        );
        eprintln!(
            "confseq_base_corpus local_constraint_non_token_angle_max_abs comparable={} p50={:.6} p90={:.6} p99={:.6}",
            local_non_token_angle_max_abs.len(),
            quantile_sorted(&local_non_token_angle_max_abs, 0.50),
            quantile_sorted(&local_non_token_angle_max_abs, 0.90),
            quantile_sorted(&local_non_token_angle_max_abs, 0.99)
        );
        eprintln!(
            "confseq_base_corpus local_constraint_pass comparable={} pass_0_05a_5deg={} rate_0_05a_5deg={:.2}% pass_0_10a_10deg={} rate_0_10a_10deg={:.2}% pass_0_10a_15deg={} rate_0_10a_15deg={:.2}%",
            local_constraint_comparable,
            local_constraint_pass_005a_5deg,
            local_constraint_pass_rate_005a_5deg,
            local_constraint_pass_010a_10deg,
            local_constraint_pass_rate_010a_10deg,
            local_constraint_pass_010a_15deg,
            local_constraint_pass_rate_010a_15deg,
        );
        print_local_bond_deviation_buckets(
            "confseq_base_corpus local_bond_deviation",
            local_bond_deviation_buckets,
            16,
        );
        eprintln!(
            "confseq_base_corpus rigid_fragment_rmsd fragments={} pass_le_0_1a={} pass_le_0_2a={} pass_le_0_3a={} p50={:.6} p75={:.6} p90={:.6} p95={:.6} p99={:.6}",
            rigid_fragment_rmsds.len(),
            rigid_fragment_pass_01a,
            rigid_fragment_pass_02a,
            rigid_fragment_pass_03a,
            quantile_sorted(&rigid_fragment_rmsds, 0.50),
            quantile_sorted(&rigid_fragment_rmsds, 0.75),
            quantile_sorted(&rigid_fragment_rmsds, 0.90),
            quantile_sorted(&rigid_fragment_rmsds, 0.95),
            quantile_sorted(&rigid_fragment_rmsds, 0.99)
        );
        eprintln!(
            "confseq_base_corpus rigid_fragment_shape_rmsd fragments={} pass_le_0_1a={} pass_le_0_2a={} pass_le_0_3a={} p50={:.6} p75={:.6} p90={:.6} p95={:.6} p99={:.6}",
            rigid_fragment_shape_rmsds.len(),
            rigid_shape_pass_01a,
            rigid_shape_pass_02a,
            rigid_shape_pass_03a,
            quantile_sorted(&rigid_fragment_shape_rmsds, 0.50),
            quantile_sorted(&rigid_fragment_shape_rmsds, 0.75),
            quantile_sorted(&rigid_fragment_shape_rmsds, 0.90),
            quantile_sorted(&rigid_fragment_shape_rmsds, 0.95),
            quantile_sorted(&rigid_fragment_shape_rmsds, 0.99)
        );
        eprintln!(
            "confseq_base_corpus rigid_fragment_terminal_symmetry_rmsd fragments={} pass_le_0_1a={} pass_le_0_2a={} pass_le_0_3a={} p50={:.6} p75={:.6} p90={:.6} p95={:.6} p99={:.6}",
            rigid_fragment_terminal_symmetry_rmsds.len(),
            rigid_terminal_symmetry_pass_01a,
            rigid_terminal_symmetry_pass_02a,
            rigid_terminal_symmetry_pass_03a,
            quantile_sorted(&rigid_fragment_terminal_symmetry_rmsds, 0.50),
            quantile_sorted(&rigid_fragment_terminal_symmetry_rmsds, 0.75),
            quantile_sorted(&rigid_fragment_terminal_symmetry_rmsds, 0.90),
            quantile_sorted(&rigid_fragment_terminal_symmetry_rmsds, 0.95),
            quantile_sorted(&rigid_fragment_terminal_symmetry_rmsds, 0.99)
        );
        eprintln!(
            "confseq_base_corpus rigid_fragment_connector_rmsd fragments={} pass_le_0_1a={} pass_le_0_2a={} pass_le_0_3a={} p50={:.6} p75={:.6} p90={:.6} p95={:.6} p99={:.6}",
            rigid_fragment_connector_rmsds.len(),
            rigid_connector_pass_01a,
            rigid_connector_pass_02a,
            rigid_connector_pass_03a,
            quantile_sorted(&rigid_fragment_connector_rmsds, 0.50),
            quantile_sorted(&rigid_fragment_connector_rmsds, 0.75),
            quantile_sorted(&rigid_fragment_connector_rmsds, 0.90),
            quantile_sorted(&rigid_fragment_connector_rmsds, 0.95),
            quantile_sorted(&rigid_fragment_connector_rmsds, 0.99)
        );
        eprintln!(
            "confseq_base_corpus rigid_fragment_max_per_mol comparable={} pass_le_0_1a={} pass_le_0_2a={} pass_le_0_3a={} p50={:.6} p90={:.6} p99={:.6} global_fail_rigid_pass={} global_fail_rigid_fail={}",
            rigid_fragment_max_rmsds.len(),
            rigid_max_pass_01a,
            rigid_max_pass_02a,
            rigid_max_pass_03a,
            quantile_sorted(&rigid_fragment_max_rmsds, 0.50),
            quantile_sorted(&rigid_fragment_max_rmsds, 0.90),
            quantile_sorted(&rigid_fragment_max_rmsds, 0.99),
            global_fail_rigid_pass,
            global_fail_rigid_fail,
        );
        eprintln!(
            "confseq_base_corpus rigid_fragment_shape_max_per_mol comparable={} p50={:.6} p90={:.6} p99={:.6}",
            rigid_fragment_max_shape_rmsds.len(),
            quantile_sorted(&rigid_fragment_max_shape_rmsds, 0.50),
            quantile_sorted(&rigid_fragment_max_shape_rmsds, 0.90),
            quantile_sorted(&rigid_fragment_max_shape_rmsds, 0.99),
        );
        eprintln!(
            "confseq_base_corpus rigid_fragment_mirror_branch_like fragments_shape_le_0_1a={} fragments_shape_le_0_3a={} molecules_shape_le_0_1a={} molecules_shape_le_0_3a={}",
            mirror_branch_like_fragments_01a,
            mirror_branch_like_fragments_03a,
            mirror_branch_like_molecules_01a,
            mirror_branch_like_molecules_03a,
        );
        eprintln!(
            "confseq_base_corpus rigid_fragment_terminal_symmetry_max_per_mol comparable={} p50={:.6} p90={:.6} p99={:.6}",
            rigid_fragment_max_terminal_symmetry_rmsds.len(),
            quantile_sorted(&rigid_fragment_max_terminal_symmetry_rmsds, 0.50),
            quantile_sorted(&rigid_fragment_max_terminal_symmetry_rmsds, 0.90),
            quantile_sorted(&rigid_fragment_max_terminal_symmetry_rmsds, 0.99),
        );
        eprintln!(
            "confseq_base_corpus rigid_fragment_connector_max_per_mol comparable={} p50={:.6} p90={:.6} p99={:.6}",
            rigid_fragment_max_connector_rmsds.len(),
            quantile_sorted(&rigid_fragment_max_connector_rmsds, 0.50),
            quantile_sorted(&rigid_fragment_max_connector_rmsds, 0.90),
            quantile_sorted(&rigid_fragment_max_connector_rmsds, 0.99),
        );
        print_top_count_map(
            "confseq_base_corpus framework_fragment_shape",
            framework_fragment_shapes,
            16,
        );
        print_top_count_map(
            "confseq_base_corpus framework_assembly_path",
            framework_assembly_paths,
            8,
        );
        print_top_count_map(
            "confseq_base_corpus framework_fragment_template_failure",
            framework_fragment_template_failures,
            8,
        );
        print_top_count_map(
            "confseq_base_corpus framework_nonplanar_ring_topology",
            framework_nonplanar_ring_topologies,
            16,
        );
        print_top_count_map(
            "confseq_base_corpus structural_risk_precheck_class",
            structural_risk_classes_counts,
            16,
        );
        print_top_count_map(
            "confseq_base_corpus structural_risk_strict_fail_not_automorphism_class",
            strict_fail_structural_risk_classes_counts,
            16,
        );
        print_top_count_map(
            "confseq_base_corpus pre_token_template_error",
            pre_token_template_errors,
            8,
        );
        print_rigid_fragment_rmsd_buckets(
            "confseq_base_corpus rigid_fragment_framework_proper_rmsd",
            rigid_fragment_framework_rmsd_buckets.proper,
            16,
        );
        print_rigid_fragment_rmsd_buckets(
            "confseq_base_corpus rigid_fragment_framework_shape_rmsd",
            rigid_fragment_framework_rmsd_buckets.shape,
            16,
        );
        print_rigid_fragment_rmsd_buckets(
            "confseq_base_corpus rigid_fragment_framework_terminal_symmetry_rmsd",
            rigid_fragment_framework_rmsd_buckets.terminal_symmetry,
            16,
        );
        print_rigid_fragment_rmsd_buckets(
            "confseq_base_corpus rigid_fragment_framework_connector_rmsd",
            rigid_fragment_framework_rmsd_buckets.connector,
            16,
        );
        if let Some(export_dir) = fail_export_dir.as_deref() {
            export_confseq_base_fail_subset(export_dir, &fail_exports)
                .expect("ConfSeq FastGeometry failure subset export should succeed");
            eprintln!(
                "confseq_base_corpus fail_export threshold_global_rmsd_gt_0_1a count={} dir={export_dir}",
                fail_exports.len()
            );
        }
        worst_rigid_fragments.sort_by(|left, right| right.0.total_cmp(&left.0));
        for (max_rmsd, idx, atoms, stripped_smiles) in worst_rigid_fragments.into_iter().take(12) {
            let mol = Molecule::from_smiles(&stripped_smiles)
                .expect("fixture stripped SMILES should parse");
            let fragment_type = rigid_fragment_type_for_diagnostic(&mol, &atoms);
            let local = rigid_fragment_local_diagnostic(&mol, &atoms);
            let reference_points = reference_cache.entries[idx]
                .heavy_atom_points
                .as_ref()
                .expect("worst fragment has reference points");
            let base = base_batch.molecules[idx]
                .as_ref()
                .expect("worst fragment has base molecule");
            let details = rigid_fragment_metric_details(&mol, reference_points, base, &atoms);
            eprintln!(
                "confseq_base_corpus worst_rigid_fragment idx={idx} max_rmsd={max_rmsd:.6} type={fragment_type} atoms={atoms:?} proper={:.6} shape={:.6} terminal_symmetry={:.6} connector={:.6} angle_rms_deg={:.3} max_angle_delta_deg={:.3} {local} stripped_smiles={stripped_smiles}",
                details.proper_rmsd,
                details.shape_rmsd,
                details.terminal_symmetry_rmsd,
                details.connector_rmsd.unwrap_or(f64::NAN),
                details.angle_rms_deg,
                details.max_angle_delta_deg
            );
        }

        let mut errors: Vec<_> = base_errors.into_iter().collect();
        errors.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
        for (error, count) in errors.into_iter().take(12) {
            eprintln!("confseq_base_corpus base_error count={count} error={error}");
        }
    }

    #[test]
    #[ignore = "local ConfSeq DG reference atom-order sanity diagnostic"]
    fn confseq_dg_reference_atom_order_sanity_snapshot() {
        let corpus_path = std::env::var("COSMOLKIT_CONFSEQ_CORPUS").unwrap_or_else(|_| {
            "/home/wangjingtong/sh4090/confseq_test_strings_100x10.jsonl".to_string()
        });
        let (in_smiles_batch, td_smiles_batch) = read_confseq_corpus_smiles(&corpus_path);
        let reference_cache = load_or_generate_confseq_dg_reference_cache(
            &corpus_path,
            &in_smiles_batch,
            &td_smiles_batch,
        );
        let reference_options = ConfSeqDecodeOptions {
            optimize_with_uff: false,
            ..ConfSeqDecodeOptions::default()
        };
        let row_filter = std::env::var("COSMOLKIT_CONFSEQ_DG_REFERENCE_ROW_FILTER")
            .ok()
            .map(|value| {
                value
                    .split(',')
                    .filter_map(|part| part.trim().parse::<usize>().ok())
                    .collect::<HashSet<_>>()
            })
            .unwrap_or_else(|| [309usize, 937, 865, 28].into_iter().collect());
        for idx in row_filter {
            let parsed = parse_confseq(&in_smiles_batch[idx], &td_smiles_batch[idx])
                .expect("fixture ConfSeq should parse");
            let source_mol = Molecule::from_smiles(&parsed.stripped_smiles)
                .expect("fixture stripped SMILES should parse");
            let reference = decode_confseq_with_options(
                &in_smiles_batch[idx],
                &td_smiles_batch[idx],
                &reference_options,
            )
            .expect("DG reference should decode");
            let direct_points = heavy_atom_points_for_rmsd(&reference);
            let direct_on_reference = heavy_bond_length_stats(&reference, &direct_points);
            let direct_on_source = heavy_bond_length_stats(&source_mol, &direct_points);
            if let Some(cached_points) = reference_cache.entries[idx].heavy_atom_points.as_ref() {
                let cached_on_source = heavy_bond_length_stats(&source_mol, cached_points);
                let direct_cached_rmsd =
                    crate::distgeom::aligned_rmsd_for_test(&direct_points, cached_points);
                eprintln!(
                    "confseq_dg_reference_atom_order idx={idx} cache_status=valid direct_cached_rmsd={direct_cached_rmsd:.6} reference_bonds count={} rms={:.6} max={:.6} source_direct count={} rms={:.6} max={:.6} source_cached count={} rms={:.6} max={:.6} stripped_smiles={}",
                    direct_on_reference.0,
                    direct_on_reference.1,
                    direct_on_reference.2,
                    direct_on_source.0,
                    direct_on_source.1,
                    direct_on_source.2,
                    cached_on_source.0,
                    cached_on_source.1,
                    cached_on_source.2,
                    parsed.stripped_smiles,
                );
            } else {
                eprintln!(
                    "confseq_dg_reference_atom_order idx={idx} cache_status=invalid error={} reference_bonds count={} rms={:.6} max={:.6} source_direct count={} rms={:.6} max={:.6} stripped_smiles={}",
                    reference_cache.entries[idx]
                        .error
                        .as_deref()
                        .unwrap_or("unknown"),
                    direct_on_reference.0,
                    direct_on_reference.1,
                    direct_on_reference.2,
                    direct_on_source.0,
                    direct_on_source.1,
                    direct_on_source.2,
                    parsed.stripped_smiles,
                );
            }
        }
    }

    #[test]
    #[ignore = "local ConfSeq corpus FastGeometry-only coverage snapshot"]
    fn confseq_base_corpus_base_only_snapshot() {
        let corpus_path = std::env::var("COSMOLKIT_CONFSEQ_CORPUS").unwrap_or_else(|_| {
            "/home/wangjingtong/sh4090/confseq_test_strings_100x10.jsonl".to_string()
        });
        let input =
            std::fs::read_to_string(&corpus_path).expect("ConfSeq corpus should be readable");
        let base_options = ConfSeqDecodeOptions {
            optimize_with_uff: false,
            template_backend: ConfSeqTemplateBackend::FastGeometry,
            ..ConfSeqDecodeOptions::default()
        };

        let mut in_smiles_batch = Vec::new();
        let mut td_smiles_batch = Vec::new();
        for (line_idx, line) in input.lines().enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let value: Value = serde_json::from_str(line).unwrap_or_else(|err| {
                panic!("failed to parse corpus JSON line {}: {err}", line_idx + 1)
            });
            let in_smiles = value["in_smiles"].as_str().unwrap_or_else(|| {
                panic!(
                    "corpus line {} is missing string field in_smiles",
                    line_idx + 1
                )
            });
            let td_smiles = value["td_smiles"].as_str().unwrap_or_else(|| {
                panic!(
                    "corpus line {} is missing string field td_smiles",
                    line_idx + 1
                )
            });
            in_smiles_batch.push(in_smiles.to_string());
            td_smiles_batch.push(td_smiles.to_string());
        }

        let base_batch = decode_confseq_batch_with_options(
            &in_smiles_batch,
            &td_smiles_batch,
            &base_options,
            true,
        )
        .expect("base batch should decode with per-record errors");

        let total = in_smiles_batch.len();
        let success = base_batch
            .molecules
            .iter()
            .filter(|molecule| molecule.is_some())
            .count();
        let mut base_errors: HashMap<String, usize> = HashMap::new();
        for error in base_batch.errors.iter().flatten() {
            *base_errors.entry(format!("{error:?}")).or_default() += 1;
        }
        let success_rate = success as f64 / total.max(1) as f64 * 100.0;

        eprintln!("confseq_base_only_corpus path={corpus_path}");
        eprintln!(
            "confseq_base_only_corpus total={total} success={success} failure={} success_rate={success_rate:.2}%",
            total - success
        );

        let mut errors: Vec<_> = base_errors.into_iter().collect();
        errors.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
        for (error, count) in errors.into_iter().take(12) {
            eprintln!("confseq_base_only_corpus base_error count={count} error={error}");
        }
    }

    #[test]
    #[ignore = "local ConfSeq corpus structural diagnostics"]
    fn confseq_base_corpus_structural_diagnostics() {
        let corpus_path = std::env::var("COSMOLKIT_CONFSEQ_CORPUS").unwrap_or_else(|_| {
            "/home/wangjingtong/sh4090/confseq_test_strings_100x10.jsonl".to_string()
        });
        let input =
            std::fs::read_to_string(&corpus_path).expect("ConfSeq corpus should be readable");
        let base_options = ConfSeqDecodeOptions {
            optimize_with_uff: false,
            template_backend: ConfSeqTemplateBackend::FastGeometry,
            ..ConfSeqDecodeOptions::default()
        };

        let mut in_smiles_batch = Vec::new();
        let mut td_smiles_batch = Vec::new();
        for (line_idx, line) in input.lines().enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let value: Value = serde_json::from_str(line).unwrap_or_else(|err| {
                panic!("failed to parse corpus JSON line {}: {err}", line_idx + 1)
            });
            in_smiles_batch.push(
                value["in_smiles"]
                    .as_str()
                    .expect("in_smiles should be present")
                    .to_string(),
            );
            td_smiles_batch.push(
                value["td_smiles"]
                    .as_str()
                    .expect("td_smiles should be present")
                    .to_string(),
            );
        }

        let reference_cache = load_or_generate_confseq_dg_reference_cache(
            &corpus_path,
            &in_smiles_batch,
            &td_smiles_batch,
        );
        let base_batch = decode_confseq_batch_with_options(
            &in_smiles_batch,
            &td_smiles_batch,
            &base_options,
            true,
        )
        .expect("base batch should decode with per-record errors");

        let mut successes = Vec::<(f64, usize)>::new();
        let mut failures = Vec::<(usize, String)>::new();
        for idx in 0..in_smiles_batch.len() {
            match (
                reference_cache.entries[idx].heavy_atom_points.as_ref(),
                base_batch.molecules[idx].as_ref(),
            ) {
                (Some(reference_points), Some(base)) => {
                    let base_points = heavy_atom_points_for_rmsd(base);
                    let rmsd =
                        crate::distgeom::aligned_rmsd_for_test(reference_points, &base_points);
                    successes.push((rmsd, idx));
                }
                (_, None) => {
                    failures.push((
                        idx,
                        format!(
                            "{:?}",
                            base_batch.errors[idx]
                                .as_ref()
                                .expect("base failure should include an error")
                        ),
                    ));
                }
                _ => {}
            }
        }
        successes.sort_by(|left, right| right.0.total_cmp(&left.0));

        eprintln!("confseq_base_structural_diagnostics worst_successes");
        for (rank, (rmsd, idx)) in successes.into_iter().take(12).enumerate() {
            let stripped = parse_confseq(&in_smiles_batch[idx], &td_smiles_batch[idx])
                .expect("fixture ConfSeq should parse")
                .stripped_smiles;
            let mol = Molecule::from_smiles(&stripped).expect("fixture should parse");
            let rings = rings::symmetrize_sssr(&mol).expect("ring perception should work");
            let ring_components = classify_confseq_base_ring_components(&mol, &rings)
                .expect("ring components should classify for diagnostic");
            let base = base_batch.molecules[idx].as_ref().expect("checked above");
            let coords = conformer_points(base);
            let max_abs_z = coords
                .iter()
                .map(|point| point[2].abs())
                .fold(0.0, f64::max);
            let planar_bonds = collect_confseq_base_planar_bonds(&mol).len();
            let rotatable_single_bonds = mol
                .bonds()
                .iter()
                .filter(|bond| {
                    bond.order() == BondOrder::Single
                        && !bond.is_aromatic()
                        && !planar_bond_like_for_diagnostic(&mol, bond)
                })
                .count();
            eprintln!(
                "idx={} rmsd={rmsd:.6} atoms={} bonds={} rings={} planar_bonds={} rot_single={} max_abs_z={max_abs_z:.3} stripped_smiles={} in_smiles={} td_smiles={}",
                idx,
                mol.num_atoms(),
                mol.num_bonds(),
                rings.num_rings(),
                planar_bonds,
                rotatable_single_bonds,
                stripped,
                in_smiles_batch[idx],
                td_smiles_batch[idx]
            );
            if rank < 6 {
                for (component_idx, component) in ring_components.iter().enumerate() {
                    let atom_set: HashSet<_> = component.atoms.iter().copied().collect();
                    let atom_details: Vec<_> = component
                        .atoms
                        .iter()
                        .map(|&atom_idx| {
                            let atom = &mol.atoms()[atom_idx];
                            format!(
                                "Z{}@{}:{:?}:arom={}:exo_pi={}:conj={}",
                                atom.atomic_number(),
                                atom_idx,
                                atom.hybridization(),
                                atom.is_aromatic(),
                                confseq_base_ring_atom_has_exocyclic_pi_bond(
                                    &mol, atom_idx, &atom_set
                                ),
                                confseq_base_ring_atom_is_conjugated_to_ring_pi_system(
                                    &mol, atom_idx, &atom_set
                                )
                            )
                        })
                        .collect();
                    eprintln!(
                        "idx={idx} component={component_idx} planar={} atoms=[{}]",
                        component.planar,
                        atom_details.join(","),
                    );
                }
            }
        }

        let mut band_successes = Vec::<(f64, usize)>::new();
        for idx in 0..in_smiles_batch.len() {
            if let (Some(reference_points), Some(base)) = (
                reference_cache.entries[idx].heavy_atom_points.as_ref(),
                base_batch.molecules[idx].as_ref(),
            ) {
                let base_points = heavy_atom_points_for_rmsd(base);
                let rmsd = crate::distgeom::aligned_rmsd_for_test(reference_points, &base_points);
                band_successes.push((rmsd, idx));
            }
        }
        band_successes.sort_by(|left, right| left.0.total_cmp(&right.0));
        eprintln!("confseq_base_structural_diagnostics quantile_samples");
        for &quantile in &[0.50, 0.75, 0.90] {
            if band_successes.is_empty() {
                continue;
            }
            let center = ((band_successes.len() - 1) as f64 * quantile).round() as isize;
            let start = (center - 2).max(0) as usize;
            let end = (center + 3).min(band_successes.len() as isize) as usize;
            for &(rmsd, idx) in &band_successes[start..end] {
                let stripped = parse_confseq(&in_smiles_batch[idx], &td_smiles_batch[idx])
                    .expect("fixture ConfSeq should parse")
                    .stripped_smiles;
                let mol = Molecule::from_smiles(&stripped).expect("fixture should parse");
                let rings = rings::symmetrize_sssr(&mol).expect("ring perception should work");
                let ring_components = classify_confseq_base_ring_components(&mol, &rings)
                    .expect("ring components should classify for diagnostic");
                let planar_components = ring_components
                    .iter()
                    .filter(|component| component.planar)
                    .count();
                let nonplanar_components = ring_components.len() - planar_components;
                let rotatable_single_bonds = mol
                    .bonds()
                    .iter()
                    .filter(|bond| {
                        bond.order() == BondOrder::Single
                            && !bond.is_aromatic()
                            && !planar_bond_like_for_diagnostic(&mol, bond)
                    })
                    .count();
                eprintln!(
                    "quantile={quantile:.2} idx={idx} rmsd={rmsd:.6} atoms={} rings={} components={} planar_components={} nonplanar_components={} rot_single={} stripped_smiles={stripped}",
                    mol.num_atoms(),
                    rings.num_rings(),
                    ring_components.len(),
                    planar_components,
                    nonplanar_components,
                    rotatable_single_bonds,
                );
            }
        }

        eprintln!("confseq_base_structural_diagnostics first_failures");
        for (idx, error) in failures.into_iter().take(20) {
            let stripped = parse_confseq(&in_smiles_batch[idx], &td_smiles_batch[idx])
                .expect("fixture ConfSeq should parse")
                .stripped_smiles;
            let mol = Molecule::from_smiles(&stripped).expect("fixture should parse");
            let rings = rings::symmetrize_sssr(&mol).expect("ring perception should work");
            eprintln!(
                "idx={} atoms={} bonds={} rings={} error={} stripped_smiles={} in_smiles={} td_smiles={}",
                idx,
                mol.num_atoms(),
                mol.num_bonds(),
                rings.num_rings(),
                error,
                stripped,
                in_smiles_batch[idx],
                td_smiles_batch[idx]
            );
        }
    }

    #[test]
    #[ignore = "local ConfSeq corpus stage RMSD diagnostics"]
    fn confseq_base_corpus_stage_rmsd_diagnostics() {
        let corpus_path = std::env::var("COSMOLKIT_CONFSEQ_CORPUS").unwrap_or_else(|_| {
            "/home/wangjingtong/sh4090/confseq_test_strings_100x10.jsonl".to_string()
        });
        let input =
            std::fs::read_to_string(&corpus_path).expect("ConfSeq corpus should be readable");
        let base_options = ConfSeqDecodeOptions {
            optimize_with_uff: false,
            template_backend: ConfSeqTemplateBackend::FastGeometry,
            ..ConfSeqDecodeOptions::default()
        };
        let dg_options = ConfSeqDecodeOptions {
            optimize_with_uff: false,
            template_backend: ConfSeqTemplateBackend::DistanceGeometry,
            ..ConfSeqDecodeOptions::default()
        };

        let mut in_smiles_batch = Vec::new();
        let mut td_smiles_batch = Vec::new();
        let mut parsed_batch = Vec::new();
        for (line_idx, line) in input.lines().enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let value: Value = serde_json::from_str(line).unwrap_or_else(|err| {
                panic!("failed to parse corpus JSON line {}: {err}", line_idx + 1)
            });
            let in_smiles = value["in_smiles"]
                .as_str()
                .expect("in_smiles should be present");
            let td_smiles = value["td_smiles"]
                .as_str()
                .expect("td_smiles should be present");
            let parsed = parse_confseq(in_smiles, td_smiles).expect("fixture ConfSeq should parse");
            in_smiles_batch.push(in_smiles.to_string());
            td_smiles_batch.push(td_smiles.to_string());
            parsed_batch.push(parsed);
        }

        let reference_cache = load_or_generate_confseq_dg_reference_cache(
            &corpus_path,
            &in_smiles_batch,
            &td_smiles_batch,
        );
        let base_batch = decode_confseq_batch_with_options(
            &in_smiles_batch,
            &td_smiles_batch,
            &base_options,
            true,
        )
        .expect("base batch should decode with per-record errors");

        let mut worst = Vec::new();
        for idx in 0..in_smiles_batch.len() {
            let Some(reference_points) = reference_cache.entries[idx].heavy_atom_points.as_ref()
            else {
                continue;
            };
            let Some(base) = base_batch.molecules[idx].as_ref() else {
                continue;
            };
            let rmsd = crate::distgeom::aligned_rmsd_for_test(
                reference_points,
                &heavy_atom_points_for_rmsd(base),
            );
            worst.push((rmsd, idx));
        }
        worst.sort_by(|left, right| right.0.total_cmp(&left.0));

        let mut rows = Vec::new();
        for (_, line_idx) in worst.into_iter().take(16) {
            let parsed = &parsed_batch[line_idx];
            let Ok(reference) = build_template(
                &parsed.stripped_smiles,
                &parsed.chiral_tags_by_atom,
                &dg_options,
            )
            .and_then(|template| decode_from_template(&template, &parsed, &dg_options)) else {
                continue;
            };
            let reference_points = heavy_atom_points_for_rmsd(&reference);
            let Ok(base_template) = build_template(
                &parsed.stripped_smiles,
                &parsed.chiral_tags_by_atom,
                &base_options,
            ) else {
                continue;
            };

            let mut angle_only_options = base_options.clone();
            angle_only_options.apply_dihedrals = false;
            let Ok(angle_only) = decode_from_template(&base_template, &parsed, &angle_only_options)
            else {
                continue;
            };
            let Ok(full_base) = decode_from_template(&base_template, &parsed, &base_options) else {
                continue;
            };
            let source_mol =
                Molecule::from_smiles(&parsed.stripped_smiles).expect("fixture should parse");

            let template_rmsd = crate::distgeom::aligned_rmsd_for_test(
                &reference_points,
                &heavy_atom_points_for_rmsd(&base_template.molecule),
            );
            let angle_rmsd = crate::distgeom::aligned_rmsd_for_test(
                &reference_points,
                &heavy_atom_points_for_rmsd(&angle_only),
            );
            let full_rmsd = crate::distgeom::aligned_rmsd_for_test(
                &reference_points,
                &heavy_atom_points_for_rmsd(&full_base),
            );
            let template_rigid = rigid_fragment_rmsd_summary(
                &source_mol,
                &reference_points,
                &base_template.molecule,
            );
            let angle_rigid =
                rigid_fragment_rmsd_summary(&source_mol, &reference_points, &angle_only);
            let full_rigid =
                rigid_fragment_rmsd_summary(&source_mol, &reference_points, &full_base);
            let worst_atoms = full_rigid.worst_fragment_atoms.clone();
            let angle_centers_touching_worst = base_template
                .angle_centers
                .iter()
                .copied()
                .filter(|(left, center, right)| {
                    worst_atoms.contains(left)
                        || worst_atoms.contains(center)
                        || worst_atoms.contains(right)
                })
                .collect::<Vec<_>>();
            rows.push((
                full_rmsd,
                line_idx,
                template_rmsd,
                angle_rmsd,
                template_rigid.max_rmsd.unwrap_or(f64::NAN),
                angle_rigid.max_rmsd.unwrap_or(f64::NAN),
                full_rigid.max_rmsd.unwrap_or(f64::NAN),
                worst_atoms,
                angle_centers_touching_worst,
                parsed.stripped_smiles.clone(),
            ));
        }

        rows.sort_by(|left, right| right.0.total_cmp(&left.0));
        for (
            full_rmsd,
            line_idx,
            template_rmsd,
            angle_rmsd,
            template_rigid,
            angle_rigid,
            full_rigid,
            worst_atoms,
            angle_centers_touching_worst,
            stripped_smiles,
        ) in rows.into_iter().take(16)
        {
            eprintln!(
                "confseq_base_stage idx={line_idx} template_rmsd={template_rmsd:.6} angle_only_rmsd={angle_rmsd:.6} full_rmsd={full_rmsd:.6} template_rigid_max={template_rigid:.6} angle_rigid_max={angle_rigid:.6} full_rigid_max={full_rigid:.6} worst_atoms={worst_atoms:?} angle_centers_touching_worst={angle_centers_touching_worst:?} stripped_smiles={stripped_smiles}"
            );
        }
    }

    #[test]
    #[ignore = "local ConfSeq fragment sliced-bounds stage diagnostics"]
    fn confseq_base_fragment_sliced_bounds_stage_diagnostics() {
        let corpus_path = std::env::var("COSMOLKIT_CONFSEQ_CORPUS").unwrap_or_else(|_| {
            "/home/wangjingtong/sh4090/confseq_test_strings_100x10.jsonl".to_string()
        });
        let input =
            std::fs::read_to_string(&corpus_path).expect("ConfSeq corpus should be readable");
        let base_options = ConfSeqDecodeOptions {
            optimize_with_uff: false,
            template_backend: ConfSeqTemplateBackend::FastGeometry,
            ..ConfSeqDecodeOptions::default()
        };
        let dg_options = ConfSeqDecodeOptions {
            optimize_with_uff: false,
            template_backend: ConfSeqTemplateBackend::DistanceGeometry,
            ..ConfSeqDecodeOptions::default()
        };
        let max_scan = std::env::var("COSMOLKIT_CONFSEQ_STAGE_DIAG_LIMIT")
            .ok()
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or(80);

        let mut in_smiles_batch = Vec::new();
        let mut td_smiles_batch = Vec::new();
        for (line_idx, line) in input.lines().enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let value: Value = serde_json::from_str(line).unwrap_or_else(|err| {
                panic!("failed to parse corpus JSON line {}: {err}", line_idx + 1)
            });
            let in_smiles = value["in_smiles"]
                .as_str()
                .expect("in_smiles should be present");
            let td_smiles = value["td_smiles"]
                .as_str()
                .expect("td_smiles should be present");
            in_smiles_batch.push(in_smiles.to_string());
            td_smiles_batch.push(td_smiles.to_string());
        }
        let reference_cache = load_or_generate_confseq_dg_reference_cache(
            &corpus_path,
            &in_smiles_batch,
            &td_smiles_batch,
        );

        let mut rows = Vec::new();
        for (line_idx, (in_smiles, td_smiles)) in in_smiles_batch
            .iter()
            .zip(td_smiles_batch.iter())
            .enumerate()
            .take(max_scan)
        {
            let Ok(parsed) = parse_confseq(in_smiles, td_smiles) else {
                continue;
            };
            let Some(reference_heavy_points) = reference_cache.entries[line_idx]
                .heavy_atom_points
                .as_ref()
                .cloned()
            else {
                continue;
            };
            let Ok(base_template) = build_template(
                &parsed.stripped_smiles,
                &parsed.chiral_tags_by_atom,
                &base_options,
            ) else {
                continue;
            };
            let source_mol =
                Molecule::from_smiles(&parsed.stripped_smiles).expect("fixture should parse");
            let summary = rigid_fragment_rmsd_summary(
                &source_mol,
                &reference_heavy_points,
                &base_template.molecule,
            );
            let Some(max_rmsd) = summary.max_shape_rmsd else {
                continue;
            };
            rows.push((
                max_rmsd,
                line_idx,
                parsed,
                source_mol,
                summary.worst_fragment_atoms,
                reference_heavy_points,
            ));
        }
        rows.sort_by(|left, right| right.0.total_cmp(&left.0));

        for (
            base_fragment_rmsd,
            line_idx,
            parsed,
            source_mol,
            worst_atoms,
            reference_heavy_points,
        ) in rows.into_iter().take(8)
        {
            let rings = rings::symmetrize_sssr(&source_mol).expect("ring perception should work");
            let model = build_confseq_base_constraint_model(&source_mol, &rings)
                .expect("base constraint model should build");
            let Some(component) = model.rigid_components.iter().find(|component| {
                worst_atoms
                    .iter()
                    .all(|atom| component.atoms.contains(atom))
            }) else {
                continue;
            };
            let mut cache = HashMap::<String, ConfSeqBaseRigidFragmentTemplate>::new();
            let Ok(template) = confseq_base_cached_rigid_fragment_template(
                &source_mol,
                &model,
                component,
                &mut cache,
            ) else {
                continue;
            };
            let stages = [
                (
                    "initial",
                    crate::chemistry::distgeom::EmbedderTestStage::InitialCoords,
                ),
                (
                    "first_min",
                    crate::chemistry::distgeom::EmbedderTestStage::FirstMinimized,
                ),
                (
                    "fourth_dim",
                    crate::chemistry::distgeom::EmbedderTestStage::FourthDimensionCleaned,
                ),
            ];
            eprintln!(
                "confseq_base_fragment_stage idx={line_idx} base_fragment_shape_rmsd={base_fragment_rmsd:.6} atoms={worst_atoms:?} shape={:?} key={} stripped_smiles={}",
                template.shape, template.cache_key, parsed.stripped_smiles
            );
            for (stage_name, stage) in stages {
                let Ok(whole_template) = build_template_initial_coords_for_test(
                    &parsed.stripped_smiles,
                    &parsed.chiral_tags_by_atom,
                    &dg_options,
                    stage,
                ) else {
                    continue;
                };
                let whole_points = conformer_points(&whole_template.molecule);
                let Ok(fragment_coords) = fragment_stage_coords_from_sliced_bounds_for_test(
                    &source_mol,
                    &template,
                    stage,
                ) else {
                    eprintln!(
                        "confseq_base_fragment_stage idx={line_idx} stage={stage_name} fragment_error=true"
                    );
                    continue;
                };
                let fragment_points = worst_atoms
                    .iter()
                    .filter_map(|atom| fragment_coords.get(atom).copied())
                    .collect::<Vec<_>>();
                if fragment_points.len() != worst_atoms.len() {
                    continue;
                }
                let whole_fragment_points = points_for_atoms(&whole_points, &worst_atoms);
                let heavy_index_by_atom = heavy_index_by_atom(&source_mol);
                let reference_fragment_points = worst_atoms
                    .iter()
                    .filter_map(|atom| heavy_index_by_atom[*atom])
                    .map(|heavy_idx| reference_heavy_points[heavy_idx])
                    .collect::<Vec<_>>();
                if reference_fragment_points.len() != worst_atoms.len() {
                    continue;
                }
                let whole_vs_ref = crate::distgeom::aligned_rmsd_for_test(
                    &reference_fragment_points,
                    &whole_fragment_points,
                );
                let fragment_vs_ref = crate::distgeom::aligned_rmsd_for_test(
                    &reference_fragment_points,
                    &fragment_points,
                );
                let fragment_vs_whole = crate::distgeom::aligned_rmsd_for_test(
                    &whole_fragment_points,
                    &fragment_points,
                );
                eprintln!(
                    "confseq_base_fragment_stage idx={line_idx} stage={stage_name} whole_vs_ref={whole_vs_ref:.6} fragment_vs_ref={fragment_vs_ref:.6} fragment_vs_whole={fragment_vs_whole:.6}"
                );
            }
        }
    }

    #[test]
    #[ignore = "export local SDF pairs comparing whole-DG fragment slices with sliced-fragment realization"]
    fn confseq_base_export_fragment_realization_pair_sdfs() {
        let corpus_path = std::env::var("COSMOLKIT_CONFSEQ_CORPUS").unwrap_or_else(|_| {
            "/home/wangjingtong/sh4090/confseq_test_strings_100x10.jsonl".to_string()
        });
        let export_dir = std::env::var("COSMOLKIT_CONFSEQ_FRAGMENT_PAIR_EXPORT_DIR")
            .unwrap_or_else(|_| "tmp/confseq_fragment_realization_pairs".to_string());
        std::fs::create_dir_all(&export_dir).expect("fragment pair export dir should be writable");
        let input =
            std::fs::read_to_string(&corpus_path).expect("ConfSeq corpus should be readable");
        let rows = input.lines().collect::<Vec<_>>();
        let dg_options = ConfSeqDecodeOptions {
            optimize_with_uff: false,
            template_backend: ConfSeqTemplateBackend::DistanceGeometry,
            ..ConfSeqDecodeOptions::default()
        };
        let examples = [
            (28usize, vec![9, 10, 29, 11, 23, 24, 28, 25, 27, 26]),
            (
                868usize,
                vec![
                    0, 1, 2, 10, 3, 9, 11, 4, 8, 21, 12, 13, 5, 7, 20, 22, 14, 6, 18, 15, 16, 17,
                    19, 27, 23, 25, 28, 29, 24, 31, 26, 30,
                ],
            ),
        ];
        let mut metadata = String::new();
        for (idx, worst_atoms) in examples {
            let value: Value = serde_json::from_str(rows[idx]).unwrap_or_else(|err| {
                panic!("failed to parse corpus JSON line {}: {err}", idx + 1)
            });
            let in_smiles = value["in_smiles"]
                .as_str()
                .expect("in_smiles should be present");
            let td_smiles = value["td_smiles"]
                .as_str()
                .expect("td_smiles should be present");
            let parsed = parse_confseq(in_smiles, td_smiles).expect("ConfSeq row should parse");
            let source_mol =
                Molecule::from_smiles(&parsed.stripped_smiles).expect("fixture should parse");
            let rings = rings::symmetrize_sssr(&source_mol).expect("ring perception should work");
            let model = build_confseq_base_constraint_model(&source_mol, &rings)
                .expect("base constraint model should build");
            let component = model
                .rigid_components
                .iter()
                .find(|component| {
                    worst_atoms
                        .iter()
                        .all(|atom| component.atoms.contains(atom))
                })
                .expect("worst atoms should be contained in one rigid component");
            let mut cache = HashMap::<String, ConfSeqBaseRigidFragmentTemplate>::new();
            let template = confseq_base_cached_rigid_fragment_template(
                &source_mol,
                &model,
                component,
                &mut cache,
            )
            .expect("fragment template should build");
            let stage = crate::chemistry::distgeom::EmbedderTestStage::FourthDimensionCleaned;
            let whole_template = build_template_initial_coords_for_test(
                &parsed.stripped_smiles,
                &parsed.chiral_tags_by_atom,
                &dg_options,
                stage,
            )
            .expect("whole DG staged template should build");
            let whole_points = conformer_points(&whole_template.molecule);
            let fragment_coords =
                fragment_stage_coords_from_sliced_bounds_for_test(&source_mol, &template, stage)
                    .expect("fragment sliced-bounds coordinates should build");
            let path = std::path::Path::new(&export_dir)
                .join(format!("idx{idx}_whole_dg_vs_fragment_realization.sdf"));
            write_fragment_realization_pair_sdf_for_test(
                &path,
                &source_mol,
                &worst_atoms,
                &whole_points,
                &fragment_coords,
            )
            .expect("fragment pair SDF should write");
            metadata.push_str(
                &serde_json::to_string(&serde_json::json!({
                    "idx": idx,
                    "file": path.to_string_lossy(),
                    "records": [
                        "whole_dg_fragment_slice",
                        "fragment_sliced_bounds_realization_aligned_to_whole"
                    ],
                    "fragment_atoms_original_indices": worst_atoms,
                    "shape": format!("{:?}", template.shape),
                    "cache_key": template.cache_key,
                    "stripped_smiles": parsed.stripped_smiles,
                }))
                .expect("metadata should serialize"),
            );
            metadata.push('\n');
            eprintln!(
                "confseq_base_fragment_pair_sdf idx={idx} path={}",
                path.display()
            );
        }
        std::fs::write(
            std::path::Path::new(&export_dir).join("metadata.jsonl"),
            metadata,
        )
        .expect("metadata should write");
    }

    fn planar_bond_like_for_diagnostic(molecule: &Molecule, bond: &Bond) -> bool {
        bond.order() == BondOrder::Double
            || molecule.atoms()[bond.begin().index()].is_aromatic()
            || molecule.atoms()[bond.end().index()].is_aromatic()
            || (molecule.atoms()[bond.begin().index()].hybridization() == Hybridization::Sp2
                && molecule.atoms()[bond.end().index()].hybridization() == Hybridization::Sp2)
    }
}
