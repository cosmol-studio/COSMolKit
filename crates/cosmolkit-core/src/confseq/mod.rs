use std::collections::{HashMap, HashSet, VecDeque};
use std::f64::consts::PI;

use rayon::prelude::*;
use thiserror::Error;

use crate::chemistry::{coordinates, distgeom, mol_transforms, rings};
use crate::smiles_write::{SmilesWriteParams, mol_to_smiles};
use crate::{
    AdjacencyList, AtomId, Bond, BondId, BondOrder, BondSpec, BondStereo, ChiralTag, Conformer3D,
    EmbedParameters, Hybridization, Molecule, get_uff_bond_stretch_params, uff_optimize_molecule,
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
    /// Experimental ConfSeq-specific base conformer. This is intentionally a
    /// strict subset, not a heuristic fallback: unsupported topologies return
    /// `ConfSeqDecodeError::BaseConformer` instead of silently switching back to
    /// distance geometry.
    BaseConformer,
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
    BaseConformer(#[from] ConfSeqBaseConformerError),
    #[error("3D embedding failed")]
    EmbedFailed,
    #[error("molecular transform failed: {0}")]
    MolTransform(String),
}

#[derive(Debug, Clone, Error, PartialEq)]
pub enum ConfSeqBaseConformerError {
    #[error("ConfSeq base conformer currently supports only connected molecules")]
    Disconnected,
    #[error(
        "ConfSeq base conformer supports only supported 3-8 member heavy-atom rings that are disjoint or edge-fused"
    )]
    UnsupportedRingSystem,
    #[error("ConfSeq base conformer does not support {ring_size}-member rings")]
    UnsupportedRingSize { ring_size: usize },
    #[error(
        "ConfSeq base conformer does not support ring atom element {atomic_number} in atom ring {ring_index}"
    )]
    UnsupportedRingElement {
        ring_index: usize,
        atomic_number: u8,
    },
    #[error("ConfSeq base conformer does not support aromatic {ring_size}-member ring")]
    UnsupportedAromaticRingSize { ring_size: usize },
    #[error(
        "ConfSeq base conformer does not support ring sharing between rings {left} and {right}: shared_atoms={shared_atoms}, shared_bonds={shared_bonds}"
    )]
    UnsupportedRingFusion {
        left: usize,
        right: usize,
        shared_atoms: usize,
        shared_bonds: usize,
    },
    #[error(
        "ConfSeq base conformer does not support closed fused/spiro ring constraint components"
    )]
    UnsupportedClosedRingFusion,
    #[error(
        "ConfSeq base conformer could not satisfy tetrahedral stereo at atom {center}: {reason}"
    )]
    UnsupportedTetrahedralStereo { center: usize, reason: String },
    #[error("ConfSeq base conformer ring geometry is invalid: {0}")]
    InvalidRingGeometry(String),
    #[error("ConfSeq base conformer traversal left atoms unplaced")]
    PlacementLeftAtomsUnplaced,
    #[error("ConfSeq base conformer failed: {0}")]
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
    angle_values_deg: Vec<f64>,
    chiral_tags_by_atom: HashMap<usize, ChiralTag>,
}

#[derive(Debug, Clone)]
struct Template {
    molecule: Molecule,
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

#[derive(Debug, Clone)]
struct ConfSeqBaseConstraintModel {
    bond_targets: Vec<f64>,
    angle_targets: HashMap<(usize, usize, usize), f64>,
    torsion_priors: HashMap<(usize, usize, usize, usize), ConfSeqBaseTorsionPrior>,
    path14_distance_priors: Vec<ConfSeqBasePath14DistancePrior>,
    planar_bonds: HashSet<(usize, usize)>,
    ring_components: Vec<ConfSeqBaseRingComponent>,
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
    rings: Vec<Vec<usize>>,
    ring_sizes_by_atom: HashMap<usize, usize>,
    planar: bool,
}

pub fn decode_confseq(in_smiles: &str, confseq: &str) -> Result<Molecule, ConfSeqDecodeError> {
    decode_confseq_with_options(in_smiles, confseq, &ConfSeqDecodeOptions::default())
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
    let dihedral_angles_by_pair = parse_atom_pair_dihedrals(&stripped_smiles, &td_smiles)?;

    Ok(ParsedConfSeq {
        stripped_smiles,
        dihedral_angles_by_pair,
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
        ConfSeqTemplateBackend::BaseConformer => {
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
        .into_iter()
        .map(|dihedral @ (_, j, k, _)| (sorted_pair(j, k), dihedral))
        .collect();
    let angle_centers = collect_angle_centers(&embedded)?;
    let ring_bond_pairs = collect_ring_bond_pairs(&embedded)?;
    let last_ring_bonds = collect_last_ring_bonds(smiles, &embedded)?;

    Ok(Template {
        molecule: embedded,
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
    let embedded = try_build_confseq_base_conformer(&molecule)?;
    let embedded = restore_embedding_molecule_state(embedded, &molecule)?;

    let dihedrals = collect_single_bond_dihedrals(&embedded);
    let dihedrals_by_pair = dihedrals
        .into_iter()
        .map(|dihedral @ (_, j, k, _)| (sorted_pair(j, k), dihedral))
        .collect();
    let angle_centers = collect_angle_centers(&embedded)?;
    let ring_bond_pairs = collect_ring_bond_pairs(&embedded)?;
    let last_ring_bonds = collect_last_ring_bonds(smiles, &embedded)?;

    Ok(Template {
        molecule: embedded,
        dihedrals_by_pair,
        angle_centers,
        ring_bond_pairs,
        last_ring_bonds,
    })
}

fn try_build_confseq_base_conformer(
    molecule: &Molecule,
) -> Result<Molecule, ConfSeqBaseConformerError> {
    if molecule.num_atoms() == 0 {
        return Err(ConfSeqBaseConformerError::Build(
            "empty molecule has no base conformer".to_string(),
        ));
    }
    let ring_info = rings::symmetrize_sssr(molecule)
        .map_err(|err| ConfSeqBaseConformerError::Build(err.to_string()))?;
    let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
    if !is_connected(molecule, &adjacency) {
        return Err(ConfSeqBaseConformerError::Disconnected);
    }
    let model = build_confseq_base_constraint_model(molecule, &ring_info)?;
    let coords = construct_confseq_base_coords(molecule, &adjacency, &model)?;

    let mut builder = molecule.to_builder();
    builder
        .add_conformer(Conformer3D::new(0, coords, true))
        .map_err(|err| ConfSeqBaseConformerError::Build(err.to_string()))?;
    let molecule = builder
        .build()
        .map_err(|err| ConfSeqBaseConformerError::Build(err.to_string()))?;
    let molecule = apply_confseq_base_double_bond_stereo(molecule)?;
    apply_confseq_base_tetrahedral_stereo(molecule, &adjacency)
}

fn build_confseq_base_constraint_model(
    molecule: &Molecule,
    ring_info: &rings::RingInfo,
) -> Result<ConfSeqBaseConstraintModel, ConfSeqBaseConformerError> {
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
    Ok(ConfSeqBaseConstraintModel {
        bond_targets,
        angle_targets,
        torsion_priors,
        path14_distance_priors,
        planar_bonds,
        ring_components,
    })
}

fn classify_confseq_base_ring_components(
    molecule: &Molecule,
    ring_info: &rings::RingInfo,
) -> Result<Vec<ConfSeqBaseRingComponent>, ConfSeqBaseConformerError> {
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
        let mut rings = Vec::new();
        for ring in component_rings {
            let ring_size = ring_info.atom_rings()[ring].len();
            rings.push(
                ring_info.atom_rings()[ring]
                    .iter()
                    .map(|atom| atom.index())
                    .collect(),
            );
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
            rings,
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
        let angle = if ring_size > 0 {
            confseq_base_ring_angle_rad(molecule, center, ring_size)
        } else {
            confseq_base_local_angle_rad(molecule, center)
        };
        for left_pos in 0..neighbors.len() {
            for right_pos in left_pos + 1..neighbors.len() {
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

fn construct_confseq_base_coords(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
) -> Result<Vec<[f64; 3]>, ConfSeqBaseConformerError> {
    let mut candidates = Vec::<Vec<[f64; 3]>>::new();
    if let Ok(coords) =
        construct_confseq_base_coords_from_rdkit_scaffold(molecule, adjacency, model)
    {
        candidates.push(coords);
    }
    if model
        .ring_components
        .iter()
        .any(|component| !component.planar)
    {
        if let Ok(coords) =
            construct_confseq_base_coords_from_local_3d_prior(molecule, adjacency, model)
        {
            candidates.push(coords);
        }
    }
    let seeded_candidates = candidates.clone();
    for coords in &seeded_candidates {
        candidates.extend(confseq_base_component_twist_candidates(
            molecule, adjacency, model, coords,
        ));
    }
    if candidates.is_empty() {
        return construct_confseq_base_coords_from_local_3d_prior(molecule, adjacency, model);
    }

    let chiral_constraints =
        collect_confseq_base_tetrahedral_stereo_constraints(molecule, adjacency)?;
    candidates
        .into_iter()
        .min_by(|left, right| {
            let left_score =
                confseq_base_candidate_score(molecule, model, &chiral_constraints, left);
            let right_score =
                confseq_base_candidate_score(molecule, model, &chiral_constraints, right);
            left_score.total_cmp(&right_score)
        })
        .ok_or(ConfSeqBaseConformerError::PlacementLeftAtomsUnplaced)
}

fn confseq_base_component_twist_candidates(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
    coords: &[[f64; 3]],
) -> Vec<Vec<[f64; 3]>> {
    let component_by_atom = confseq_base_ring_component_by_atom(molecule.num_atoms(), model);
    let mut out = Vec::new();
    let seed_path14_penalty = confseq_base_path14_distance_penalty(model, coords);
    for bond in molecule.bonds() {
        if bond.order() != BondOrder::Single || bond.is_aromatic() {
            continue;
        }
        let a = bond.begin().index();
        let b = bond.end().index();
        let a_component = component_by_atom[a];
        let b_component = component_by_atom[b];
        if a_component == b_component {
            continue;
        }
        let rotates_component = match (a_component, b_component) {
            (Some(_), Some(right)) => Some((a, b, model.ring_components[right].atoms.clone())),
            (Some(_), None) => Some((
                a,
                b,
                confseq_base_side_atoms_across_bond(adjacency, molecule.num_atoms(), a, b),
            )),
            (None, Some(right)) => Some((a, b, model.ring_components[right].atoms.clone())),
            (None, None) => None,
        };
        let Some((anchor, pivot, rotate_atoms)) = rotates_component else {
            continue;
        };
        if rotate_atoms.is_empty() || !rotate_atoms.contains(&pivot) {
            continue;
        }
        for angle in [
            60.0_f64.to_radians(),
            90.0_f64.to_radians(),
            120.0_f64.to_radians(),
            -60.0_f64.to_radians(),
            -90.0_f64.to_radians(),
            -120.0_f64.to_radians(),
        ] {
            let candidate =
                confseq_base_rotate_atoms_around_bond(coords, anchor, pivot, &rotate_atoms, angle);
            let candidate_path14_penalty = confseq_base_path14_distance_penalty(model, &candidate);
            if candidate_path14_penalty + 1.0e-6 < seed_path14_penalty * 0.75 {
                out.push(candidate);
            }
        }
    }
    out
}

fn confseq_base_path14_distance_penalty(
    model: &ConfSeqBaseConstraintModel,
    coords: &[[f64; 3]],
) -> f64 {
    let mut penalty = 0.0;
    for prior in &model.path14_distance_priors {
        let (i, _, _, l) = prior.atoms;
        let observed = vec_len(vec_sub(coords[i], coords[l]));
        let deficit = if observed < prior.lower_bound {
            prior.lower_bound - observed
        } else if observed > prior.upper_bound {
            observed - prior.upper_bound
        } else {
            0.0
        };
        penalty += deficit * deficit * 1.5;
    }
    penalty
}

fn confseq_base_ring_component_by_atom(
    atom_count: usize,
    model: &ConfSeqBaseConstraintModel,
) -> Vec<Option<usize>> {
    let mut component_by_atom = vec![None; atom_count];
    for (component_idx, component) in model.ring_components.iter().enumerate() {
        for &atom in &component.atoms {
            component_by_atom[atom] = Some(component_idx);
        }
    }
    component_by_atom
}

fn confseq_base_side_atoms_across_bond(
    adjacency: &AdjacencyList,
    atom_count: usize,
    anchor: usize,
    pivot: usize,
) -> Vec<usize> {
    let mut atoms = Vec::new();
    let mut seen = vec![false; atom_count];
    let mut queue = VecDeque::new();
    seen[anchor] = true;
    seen[pivot] = true;
    queue.push_back(pivot);
    while let Some(atom) = queue.pop_front() {
        atoms.push(atom);
        for neighbor in adjacency.neighbors_of(atom) {
            if !seen[neighbor.atom_index] {
                seen[neighbor.atom_index] = true;
                queue.push_back(neighbor.atom_index);
            }
        }
    }
    atoms
}

fn confseq_base_rotate_atoms_around_bond(
    coords: &[[f64; 3]],
    anchor: usize,
    pivot: usize,
    rotate_atoms: &[usize],
    angle: f64,
) -> Vec<[f64; 3]> {
    let mut rotated = coords.to_vec();
    let axis = vec_normalize(vec_sub(coords[pivot], coords[anchor]));
    let sin_theta = angle.sin();
    let cos_theta = angle.cos();
    for &atom in rotate_atoms {
        if atom == anchor {
            continue;
        }
        let offset = vec_sub(coords[atom], coords[pivot]);
        rotated[atom] = vec_add(
            coords[pivot],
            rotate_vec_around_unit_axis(offset, axis, sin_theta, cos_theta),
        );
    }
    rotated
}

fn confseq_base_candidate_score(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    chiral_constraints: &[ConfSeqBaseTetrahedralStereoConstraint],
    coords: &[[f64; 3]],
) -> f64 {
    let mut score = 0.0;
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
        let delta = observed - target;
        score += delta * delta * 8.0;
    }
    for (&(left, center, right), &target) in &model.angle_targets {
        let observed =
            angle_rad_from_points(coords[left], coords[center], coords[right]).unwrap_or(PI);
        let delta = angular_delta_rad(observed, target);
        score += delta * delta;
    }
    for constraint in chiral_constraints {
        let volume = confseq_base_chiral_volume(coords, constraint);
        let target = confseq_base_min_signed_chiral_volume_for_constraint(constraint);
        let deficit = match constraint.tag {
            ChiralTag::TetrahedralCcw => (target - volume).max(0.0),
            ChiralTag::TetrahedralCw => (volume - target).max(0.0),
            _ => 0.0,
        };
        score += deficit * deficit * 0.04;
    }
    for prior in &model.path14_distance_priors {
        let (i, _, _, l) = prior.atoms;
        let observed = vec_len(vec_sub(coords[i], coords[l]));
        let deficit = if observed < prior.lower_bound {
            prior.lower_bound - observed
        } else if observed > prior.upper_bound {
            observed - prior.upper_bound
        } else {
            0.0
        };
        score += deficit * deficit * 1.5;
    }
    score
}

fn construct_confseq_base_coords_from_local_3d_prior(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
) -> Result<Vec<[f64; 3]>, ConfSeqBaseConformerError> {
    let mut coords = vec![[0.0; 3]; molecule.num_atoms()];
    let mut placed = vec![false; molecule.num_atoms()];
    let ring_membership =
        confseq_base_ring_membership(molecule.num_atoms(), &model.ring_components);
    let scaffold = if model.ring_components.is_empty() {
        None
    } else {
        Some(
            coordinates::rdkit_initial_2d_scaffold_coords(molecule.atoms(), molecule.bonds())
                .map_err(|err| {
                    ConfSeqBaseConformerError::Build(format!(
                        "RDKit initial scaffold construction failed for ring components: {err}"
                    ))
                })?,
        )
    };

    if let Some(component) = model.ring_components.first() {
        if component.planar {
            place_confseq_constraint_ring_component(
                component,
                scaffold.as_deref().expect("ring scaffold should exist"),
                &mut coords,
                &mut placed,
            );
        } else {
            place_confseq_nonplanar_ring_component(
                molecule,
                adjacency,
                model,
                component,
                None,
                &mut coords,
                &mut placed,
            )?;
        }
    }

    if model.ring_components.is_empty() {
        place_confseq_constraint_tree(
            molecule,
            adjacency,
            model,
            &ring_membership,
            0,
            None,
            [0.0, 0.0, 0.0],
            None,
            &mut coords,
            &mut placed,
        )?;
    } else {
        propagate_confseq_constraint_from_placed_atoms(
            molecule,
            adjacency,
            model,
            &ring_membership,
            scaffold.as_deref(),
            &mut coords,
            &mut placed,
        )?;
    }

    if placed.iter().any(|value| !*value) {
        return Err(ConfSeqBaseConformerError::PlacementLeftAtomsUnplaced);
    }
    validate_confseq_base_constraint_coords(molecule, model, &coords)?;
    Ok(coords)
}

fn construct_confseq_base_coords_from_rdkit_scaffold(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
) -> Result<Vec<[f64; 3]>, ConfSeqBaseConformerError> {
    let scaffold =
        coordinates::rdkit_initial_2d_scaffold_coords(molecule.atoms(), molecule.bonds()).map_err(
            |err| {
                ConfSeqBaseConformerError::Build(format!(
                    "RDKit initial scaffold construction failed for base conformer: {err}"
                ))
            },
        )?;
    let mut scaled_coords: Vec<_> = scaffold
        .into_iter()
        .map(|point| [point[0], point[1], 0.0])
        .collect();
    rescale_confseq_base_scaffold_bonds(molecule, model, &mut scaled_coords);
    let coords = if validate_confseq_base_constraint_coords(molecule, model, &scaled_coords).is_ok()
    {
        scaled_coords
    } else {
        let scaffold_direction_coords = construct_confseq_base_coords_from_scaffold_directions(
            molecule,
            adjacency,
            model,
            &scaled_coords,
        );
        if validate_confseq_base_constraint_coords(molecule, model, &scaffold_direction_coords)
            .is_ok()
        {
            scaffold_direction_coords
        } else if !model.ring_components.is_empty() {
            construct_confseq_base_coords_from_local_3d_prior(molecule, adjacency, model)?
        } else {
            scaffold_direction_coords
        }
    };
    validate_confseq_base_constraint_coords(molecule, model, &coords)?;
    Ok(coords)
}

fn rescale_confseq_base_scaffold_bonds(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    coords: &mut [[f64; 3]],
) {
    let mut observed_sum = 0.0;
    let mut target_sum = 0.0;
    let mut count = 0usize;
    for bond in molecule.bonds() {
        let observed = vec_len(vec_sub(
            coords[bond.begin().index()],
            coords[bond.end().index()],
        ));
        if observed <= 1.0e-10 {
            continue;
        }
        let target = model
            .bond_targets
            .get(bond.id().index())
            .copied()
            .unwrap_or_else(|| confseq_base_static_bond_length_fallback(bond));
        observed_sum += observed;
        target_sum += target;
        count += 1;
    }
    if count == 0 || observed_sum <= 1.0e-10 {
        return;
    }
    let scale = target_sum / observed_sum;
    let center = confseq_base_coord_centroid(coords);
    for coord in coords {
        *coord = vec_add(center, vec_scale(vec_sub(*coord, center), scale));
    }
}

fn construct_confseq_base_coords_from_scaffold_directions(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
    scaffold: &[[f64; 3]],
) -> Vec<[f64; 3]> {
    let mut coords = vec![[0.0; 3]; molecule.num_atoms()];
    let mut placed = vec![false; molecule.num_atoms()];
    let mut queue = VecDeque::new();
    placed[0] = true;
    queue.push_back(0);

    while let Some(parent) = queue.pop_front() {
        for neighbor in adjacency.neighbors_of(parent) {
            let child = neighbor.atom_index;
            if placed[child] {
                continue;
            }
            let direction = vec_normalize(vec_sub(scaffold[child], scaffold[parent]));
            let length = model
                .bond_targets
                .get(neighbor.bond.index())
                .copied()
                .unwrap_or_else(|| {
                    confseq_base_static_bond_length_fallback(
                        &molecule.bonds()[neighbor.bond.index()],
                    )
                });
            coords[child] = vec_add(coords[parent], vec_scale(direction, length));
            placed[child] = true;
            queue.push_back(child);
        }
    }

    coords
}

fn place_confseq_constraint_ring_component(
    component: &ConfSeqBaseRingComponent,
    scaffold: &[[f64; 2]],
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
) {
    for &atom in &component.atoms {
        let point = scaffold[atom];
        coords[atom] = [point[0], point[1], 0.0];
        placed[atom] = true;
    }
}

#[allow(clippy::too_many_arguments)]
fn place_confseq_constraint_tree(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
    ring_membership: &[Option<(usize, usize)>],
    atom_idx: usize,
    parent: Option<usize>,
    point: [f64; 3],
    incoming_axis: Option<[f64; 3]>,
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
) -> Result<(), ConfSeqBaseConformerError> {
    if placed[atom_idx] {
        return Ok(());
    }
    coords[atom_idx] = point;
    placed[atom_idx] = true;

    let axis = incoming_axis.unwrap_or([1.0, 0.0, 0.0]);
    let children: Vec<_> = adjacency
        .neighbors_of(atom_idx)
        .iter()
        .filter(|neighbor| Some(neighbor.atom_index) != parent)
        .filter(|neighbor| !placed[neighbor.atom_index])
        .collect();
    let child_count = children.len();
    for (child_ord, neighbor) in children.into_iter().enumerate() {
        let length = model
            .bond_targets
            .get(neighbor.bond.index())
            .copied()
            .unwrap_or_else(|| {
                confseq_base_static_bond_length_fallback(&molecule.bonds()[neighbor.bond.index()])
            });
        let angle = parent
            .and_then(|parent| {
                model
                    .angle_targets
                    .get(&sorted_angle(parent, atom_idx, neighbor.atom_index))
                    .copied()
            })
            .unwrap_or_else(|| confseq_base_local_angle_rad(molecule, atom_idx));
        let dir = confseq_constraint_child_direction(
            molecule,
            adjacency,
            model,
            coords,
            placed,
            atom_idx,
            parent,
            neighbor.atom_index,
            axis,
            child_ord,
            child_count,
            angle,
            length,
        );
        let child_point = vec_add(coords[atom_idx], vec_scale(dir, length));
        place_confseq_constraint_tree(
            molecule,
            adjacency,
            model,
            ring_membership,
            neighbor.atom_index,
            Some(atom_idx),
            child_point,
            Some(dir),
            coords,
            placed,
        )?;
    }
    Ok(())
}

fn propagate_confseq_constraint_from_placed_atoms(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
    ring_membership: &[Option<(usize, usize)>],
    scaffold: Option<&[[f64; 2]]>,
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
) -> Result<(), ConfSeqBaseConformerError> {
    loop {
        let mut progressed = false;
        for atom_idx in 0..molecule.num_atoms() {
            if !placed[atom_idx] {
                continue;
            }
            let base_axis = adjacency
                .neighbors_of(atom_idx)
                .iter()
                .find(|neighbor| placed[neighbor.atom_index])
                .map(|neighbor| {
                    vec_normalize(vec_sub(coords[atom_idx], coords[neighbor.atom_index]))
                })
                .unwrap_or([1.0, 0.0, 0.0]);
            let children: Vec<_> = adjacency
                .neighbors_of(atom_idx)
                .iter()
                .filter(|neighbor| !placed[neighbor.atom_index])
                .collect();
            for (child_ord, neighbor) in children.iter().enumerate() {
                let length = model
                    .bond_targets
                    .get(neighbor.bond.index())
                    .copied()
                    .unwrap_or_else(|| {
                        confseq_base_static_bond_length_fallback(
                            &molecule.bonds()[neighbor.bond.index()],
                        )
                    });
                let angle = confseq_base_local_angle_rad(molecule, atom_idx);
                let parent = adjacency
                    .neighbors_of(atom_idx)
                    .iter()
                    .find(|candidate| placed[candidate.atom_index])
                    .map(|candidate| candidate.atom_index);
                let angle = parent
                    .and_then(|parent| {
                        model
                            .angle_targets
                            .get(&sorted_angle(parent, atom_idx, neighbor.atom_index))
                            .copied()
                    })
                    .unwrap_or(angle);
                let dir = confseq_constraint_child_direction(
                    molecule,
                    adjacency,
                    model,
                    coords,
                    placed,
                    atom_idx,
                    parent,
                    neighbor.atom_index,
                    base_axis,
                    child_ord,
                    children.len(),
                    angle,
                    length,
                );
                let placed_neighbors: Vec<_> = adjacency
                    .neighbors_of(neighbor.atom_index)
                    .iter()
                    .filter(|candidate| placed[candidate.atom_index])
                    .map(|candidate| candidate.atom_index)
                    .collect();
                if ring_membership[neighbor.atom_index].is_none() && placed_neighbors.len() >= 2 {
                    if let Some(point) = place_atom_from_two_placed_neighbors(
                        molecule,
                        model,
                        coords,
                        neighbor.atom_index,
                        placed_neighbors[0],
                        placed_neighbors[1],
                        dir,
                    ) {
                        coords[neighbor.atom_index] = point;
                        placed[neighbor.atom_index] = true;
                        progressed = true;
                        continue;
                    }
                }
                if let Some((component_idx, _)) = ring_membership[neighbor.atom_index] {
                    if !model.ring_components[component_idx].planar {
                        if model.ring_components[component_idx]
                            .atoms
                            .iter()
                            .any(|atom| placed[*atom])
                        {
                            continue;
                        }
                        place_confseq_nonplanar_ring_component(
                            molecule,
                            adjacency,
                            model,
                            &model.ring_components[component_idx],
                            Some(ConfSeqNonplanarRingAttachment {
                                anchor: atom_idx,
                                attachment: neighbor.atom_index,
                                bond_length: length,
                                direction: dir,
                            }),
                            coords,
                            placed,
                        )?;
                        progressed = true;
                        continue;
                    }
                    if model.ring_components[component_idx]
                        .atoms
                        .iter()
                        .any(|atom| placed[*atom])
                    {
                        continue;
                    }
                    place_confseq_constraint_ring_component_from_attachment(
                        &model.ring_components[component_idx],
                        scaffold.ok_or_else(|| {
                            ConfSeqBaseConformerError::Build(
                                "ring component placement requires RDKit scaffold coordinates"
                                    .to_string(),
                            )
                        })?,
                        atom_idx,
                        neighbor.atom_index,
                        length,
                        dir,
                        coords,
                        placed,
                    )?;
                } else {
                    coords[neighbor.atom_index] = vec_add(coords[atom_idx], vec_scale(dir, length));
                    placed[neighbor.atom_index] = true;
                }
                progressed = true;
            }
        }
        if !progressed {
            break;
        }
    }
    Ok(())
}

#[derive(Debug, Clone, Copy)]
struct ConfSeqNonplanarRingAttachment {
    anchor: usize,
    attachment: usize,
    bond_length: f64,
    direction: [f64; 3],
}

fn place_confseq_nonplanar_ring_component(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
    component: &ConfSeqBaseRingComponent,
    attachment: Option<ConfSeqNonplanarRingAttachment>,
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
) -> Result<(), ConfSeqBaseConformerError> {
    let component_atoms: HashSet<_> = component.atoms.iter().copied().collect();
    let seed = if let Some(attachment) = attachment {
        coords[attachment.attachment] = vec_add(
            coords[attachment.anchor],
            vec_scale(attachment.direction, attachment.bond_length),
        );
        placed[attachment.attachment] = true;
        attachment.attachment
    } else {
        let seed = component.atoms[0];
        coords[seed] = [0.0, 0.0, 0.0];
        placed[seed] = true;
        seed
    };

    initialize_confseq_nonplanar_ring_path(
        molecule, adjacency, model, component, seed, attachment, coords, placed,
    )?;

    for _ in 0..component.atoms.len().saturating_mul(3).max(1) {
        let mut progressed = false;
        for ring in &component.rings {
            if let Some(segment) =
                confseq_ring_unplaced_segment_between_placed_anchors(ring, placed)
            {
                place_confseq_nonplanar_ring_segment(
                    molecule, model, ring, &segment, coords, placed,
                )?;
                progressed = true;
            }
        }
        if progressed {
            continue;
        }
        for &atom in &component.atoms {
            if placed[atom] {
                continue;
            }
            let placed_neighbors: Vec<_> = adjacency
                .neighbors_of(atom)
                .iter()
                .filter(|candidate| {
                    component_atoms.contains(&candidate.atom_index) && placed[candidate.atom_index]
                })
                .map(|candidate| candidate.atom_index)
                .collect();
            if placed_neighbors.len() < 2 {
                continue;
            }
            let preferred_dir = confseq_preferred_direction_from_two_neighbors(
                coords[placed_neighbors[0]],
                coords[placed_neighbors[1]],
            );
            if let Some(point) = place_atom_from_two_placed_neighbors(
                molecule,
                model,
                coords,
                atom,
                placed_neighbors[0],
                placed_neighbors[1],
                preferred_dir,
            ) {
                coords[atom] = point;
                placed[atom] = true;
                progressed = true;
            }
        }
        if progressed {
            continue;
        }
        for &parent in &component.atoms {
            if !placed[parent] {
                continue;
            }
            let parent_axis = adjacency
                .neighbors_of(parent)
                .iter()
                .find(|neighbor| placed[neighbor.atom_index])
                .map(|neighbor| vec_normalize(vec_sub(coords[parent], coords[neighbor.atom_index])))
                .unwrap_or([1.0, 0.0, 0.0]);
            let children: Vec<_> = adjacency
                .neighbors_of(parent)
                .iter()
                .filter(|neighbor| component_atoms.contains(&neighbor.atom_index))
                .filter(|neighbor| !placed[neighbor.atom_index])
                .collect();
            for (child_ord, neighbor) in children.iter().enumerate() {
                let child = neighbor.atom_index;
                let length = confseq_base_bond_target_by_id(model, molecule, neighbor.bond.index());
                let angle_parent = adjacency
                    .neighbors_of(parent)
                    .iter()
                    .find(|candidate| placed[candidate.atom_index])
                    .map(|candidate| candidate.atom_index);
                let angle = angle_parent
                    .and_then(|angle_parent| {
                        model
                            .angle_targets
                            .get(&sorted_angle(angle_parent, parent, child))
                            .copied()
                    })
                    .unwrap_or_else(|| {
                        confseq_base_ring_angle_rad(
                            molecule,
                            parent,
                            component
                                .ring_sizes_by_atom
                                .get(&parent)
                                .copied()
                                .unwrap_or(component.atoms.len()),
                        )
                    });
                let dir = confseq_constraint_child_direction(
                    molecule,
                    adjacency,
                    model,
                    coords,
                    placed,
                    parent,
                    angle_parent,
                    child,
                    parent_axis,
                    child_ord,
                    children.len(),
                    angle,
                    length,
                );
                let placed_neighbors: Vec<_> = adjacency
                    .neighbors_of(child)
                    .iter()
                    .filter(|candidate| placed[candidate.atom_index])
                    .map(|candidate| candidate.atom_index)
                    .collect();
                coords[child] = if placed_neighbors.len() >= 2 {
                    place_atom_from_two_placed_neighbors(
                        molecule,
                        model,
                        coords,
                        child,
                        placed_neighbors[0],
                        placed_neighbors[1],
                        dir,
                    )
                    .unwrap_or_else(|| vec_add(coords[parent], vec_scale(dir, length)))
                } else {
                    vec_add(coords[parent], vec_scale(dir, length))
                };
                placed[child] = true;
                progressed = true;
            }
        }
        if !progressed {
            break;
        }
    }

    if component.atoms.iter().any(|atom| !placed[*atom]) {
        return Err(ConfSeqBaseConformerError::PlacementLeftAtomsUnplaced);
    }
    Ok(())
}

fn confseq_ring_unplaced_segment_between_placed_anchors(
    ring: &[usize],
    placed: &[bool],
) -> Option<Vec<usize>> {
    if ring.len() < 3 {
        return None;
    }
    let mut best = None::<Vec<usize>>;
    for start_pos in 0..ring.len() {
        if !placed[ring[start_pos]] {
            continue;
        }
        for distance in 2..ring.len() {
            let end_pos = (start_pos + distance) % ring.len();
            if !placed[ring[end_pos]] {
                continue;
            }
            let mut segment = Vec::with_capacity(distance + 1);
            let mut has_unplaced = false;
            for offset in 0..=distance {
                let atom = ring[(start_pos + offset) % ring.len()];
                has_unplaced |= !placed[atom];
                segment.push(atom);
            }
            if has_unplaced
                && best
                    .as_ref()
                    .map(|current| segment.len() < current.len())
                    .unwrap_or(true)
            {
                best = Some(segment);
            }
        }
    }
    best
}

fn place_confseq_nonplanar_ring_segment(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    ring: &[usize],
    segment: &[usize],
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
) -> Result<(), ConfSeqBaseConformerError> {
    if segment.len() < 3 {
        return Ok(());
    }
    let oriented_ring = confseq_oriented_ring_from_segment(ring, segment)?;
    let local = closed_confseq_nonplanar_ring_coords(molecule, model, &oriented_ring)?;
    let start = segment[0];
    let end = *segment.last().expect("checked above");
    let Some(local_start_pos) = oriented_ring.iter().position(|atom| *atom == start) else {
        return Err(ConfSeqBaseConformerError::Build(format!(
            "nonplanar ring segment start atom {start} is absent from ring"
        )));
    };
    let Some(local_end_pos) = oriented_ring.iter().position(|atom| *atom == end) else {
        return Err(ConfSeqBaseConformerError::Build(format!(
            "nonplanar ring segment end atom {end} is absent from ring"
        )));
    };
    let local_anchor = vec_sub(local[local_end_pos], local[local_start_pos]);
    let target_anchor = vec_sub(coords[end], coords[start]);
    let local_len = vec_len(local_anchor);
    let target_len = vec_len(target_anchor);
    if local_len <= 1.0e-10 || target_len <= 1.0e-10 {
        return Ok(());
    }
    let scale = target_len / local_len;

    for &atom in segment.iter().skip(1).take(segment.len() - 2) {
        if placed[atom] {
            continue;
        }
        let Some(local_pos) = oriented_ring
            .iter()
            .position(|candidate| *candidate == atom)
        else {
            return Err(ConfSeqBaseConformerError::Build(format!(
                "nonplanar ring segment atom {atom} is absent from oriented ring"
            )));
        };
        let offset = vec_sub(local[local_pos], local[local_start_pos]);
        let offset = rotate_vector_between_unit_dirs(offset, local_anchor, target_anchor)?;
        coords[atom] = vec_add(coords[start], vec_scale(offset, scale));
        placed[atom] = true;
    }
    Ok(())
}

fn confseq_oriented_ring_from_segment(
    ring: &[usize],
    segment: &[usize],
) -> Result<Vec<usize>, ConfSeqBaseConformerError> {
    if ring.is_empty() || segment.len() < 2 {
        return Ok(ring.to_vec());
    }
    let start = segment[0];
    let next = segment[1];
    let Some(start_pos) = ring.iter().position(|atom| *atom == start) else {
        return Err(ConfSeqBaseConformerError::Build(format!(
            "nonplanar ring segment start atom {start} is absent from ring"
        )));
    };
    let forward_next = ring[(start_pos + 1) % ring.len()];
    let reverse_next = ring[(start_pos + ring.len() - 1) % ring.len()];
    let forward = if forward_next == next {
        true
    } else if reverse_next == next {
        false
    } else {
        return Err(ConfSeqBaseConformerError::Build(format!(
            "nonplanar ring segment {start}-{next} is not consecutive in ring"
        )));
    };
    let mut out = Vec::with_capacity(ring.len());
    for offset in 0..ring.len() {
        let pos = if forward {
            (start_pos + offset) % ring.len()
        } else {
            (start_pos + ring.len() - offset) % ring.len()
        };
        out.push(ring[pos]);
    }
    Ok(out)
}

fn confseq_preferred_direction_from_two_neighbors(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    let midpoint = vec_scale(vec_add(a, b), 0.5);
    let span = vec_sub(b, a);
    let normal = perpendicular_unit(span);
    vec_normalize(vec_add(midpoint, normal))
}

#[allow(clippy::too_many_arguments)]
fn initialize_confseq_nonplanar_ring_path(
    molecule: &Molecule,
    _adjacency: &AdjacencyList,
    model: &ConfSeqBaseConstraintModel,
    component: &ConfSeqBaseRingComponent,
    seed: usize,
    attachment: Option<ConfSeqNonplanarRingAttachment>,
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
) -> Result<(), ConfSeqBaseConformerError> {
    let Some(ring) = component
        .rings
        .iter()
        .find(|ring| ring.contains(&seed))
        .or_else(|| component.rings.first())
    else {
        return Ok(());
    };
    let seed_pos = ring.iter().position(|atom| *atom == seed).unwrap_or(0);
    let mut path = Vec::with_capacity(ring.len());
    for offset in 0..ring.len() {
        path.push(ring[(seed_pos + offset) % ring.len()]);
    }

    if path.len() < 3 {
        return Ok(());
    }

    let local = closed_confseq_nonplanar_ring_coords(molecule, model, &path)?;
    let first = path[0];
    let second = path[1];
    let target_origin = coords[first];
    let local_seed_to_second = vec_normalize(vec_sub(local[1], local[0]));
    let target_seed_to_second = attachment
        .map(|attachment| vec_scale(attachment.direction, -1.0))
        .unwrap_or([1.0, 0.0, 0.0]);
    for (pos, &atom) in path.iter().enumerate() {
        if placed[atom] && atom != first {
            continue;
        }
        let offset = vec_sub(local[pos], local[0]);
        let offset =
            rotate_vector_between_unit_dirs(offset, local_seed_to_second, target_seed_to_second)?;
        coords[atom] = vec_add(target_origin, offset);
        placed[atom] = true;
    }
    placed[second] = true;
    Ok(())
}

fn closed_confseq_nonplanar_ring_coords(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    ring: &[usize],
) -> Result<Vec<[f64; 3]>, ConfSeqBaseConformerError> {
    let n = ring.len();
    let mut edge_lengths = Vec::with_capacity(n);
    for pos in 0..n {
        let a = ring[pos];
        let b = ring[(pos + 1) % n];
        let Some(bond) = bond_between_pair(molecule, sorted_pair(a, b)) else {
            return Err(ConfSeqBaseConformerError::Build(format!(
                "nonplanar closed ring has no bond {a}-{b}"
            )));
        };
        edge_lengths.push(confseq_base_bond_target_by_id(
            model,
            molecule,
            bond.id().index(),
        ));
    }

    let amplitude = confseq_nonplanar_ring_pucker_amplitude(molecule, ring);
    let z_values: Vec<_> = (0..n)
        .map(|pos| {
            if n % 2 == 0 {
                if pos % 2 == 0 { amplitude } else { -amplitude }
            } else {
                amplitude * (2.0 * PI * pos as f64 / n as f64).sin()
            }
        })
        .collect();
    let mut xy_lengths = Vec::with_capacity(n);
    for pos in 0..n {
        let dz = z_values[(pos + 1) % n] - z_values[pos];
        xy_lengths.push(
            (edge_lengths[pos] * edge_lengths[pos] - dz * dz)
                .max(0.25)
                .sqrt(),
        );
    }
    let radius = closed_polygon_radius_for_chords(&xy_lengths);
    let mut angles = Vec::with_capacity(n);
    angles.push(0.0);
    let mut theta = 0.0;
    for chord in xy_lengths.iter().take(n - 1) {
        theta += 2.0 * (chord / (2.0 * radius)).clamp(-1.0, 1.0).asin();
        angles.push(theta);
    }

    let mut coords: Vec<_> = angles
        .into_iter()
        .zip(z_values)
        .map(|(angle, z)| [radius * angle.cos(), radius * angle.sin(), z])
        .collect();
    let centroid = confseq_base_coord_centroid(&coords);
    for coord in &mut coords {
        *coord = vec_sub(*coord, centroid);
    }
    Ok(coords)
}

fn closed_polygon_radius_for_chords(chords: &[f64]) -> f64 {
    let max_chord = chords.iter().copied().fold(0.0, f64::max);
    let mut low = (max_chord * 0.5).max(1.0e-6);
    let mut high = low;
    while chords
        .iter()
        .map(|chord| 2.0 * (chord / (2.0 * high)).clamp(-1.0, 1.0).asin())
        .sum::<f64>()
        > 2.0 * PI
    {
        high *= 2.0;
    }
    for _ in 0..32 {
        let mid = 0.5 * (low + high);
        let sum = chords
            .iter()
            .map(|chord| 2.0 * (chord / (2.0 * mid)).clamp(-1.0, 1.0).asin())
            .sum::<f64>();
        if sum > 2.0 * PI {
            low = mid;
        } else {
            high = mid;
        }
    }
    high
}

fn confseq_nonplanar_ring_pucker_amplitude(molecule: &Molecule, ring: &[usize]) -> f64 {
    if ring.iter().all(|&atom| {
        molecule.atoms()[atom].is_aromatic()
            || molecule.atoms()[atom].hybridization() == Hybridization::Sp2
    }) {
        return 0.0;
    }
    match ring.len() {
        3 | 4 => 0.12,
        5 => 0.22,
        6 => 0.32,
        _ => 0.24,
    }
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

fn place_atom_from_two_placed_neighbors(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    coords: &[[f64; 3]],
    atom: usize,
    neighbor_a: usize,
    neighbor_b: usize,
    preferred_dir: [f64; 3],
) -> Option<[f64; 3]> {
    let bond_a = bond_between_pair(molecule, sorted_pair(atom, neighbor_a))?;
    let bond_b = bond_between_pair(molecule, sorted_pair(atom, neighbor_b))?;
    let radius_a = model
        .bond_targets
        .get(bond_a.id().index())
        .copied()
        .unwrap_or_else(|| confseq_base_static_bond_length_fallback(bond_a));
    let radius_b = model
        .bond_targets
        .get(bond_b.id().index())
        .copied()
        .unwrap_or_else(|| confseq_base_static_bond_length_fallback(bond_b));
    sphere_sphere_intersection_point(
        coords[neighbor_a],
        radius_a,
        coords[neighbor_b],
        radius_b,
        preferred_dir,
    )
}

fn sphere_sphere_intersection_point(
    center_a: [f64; 3],
    radius_a: f64,
    center_b: [f64; 3],
    radius_b: f64,
    preferred_dir: [f64; 3],
) -> Option<[f64; 3]> {
    let axis = vec_sub(center_b, center_a);
    let distance = vec_len(axis);
    if distance <= 1.0e-10 {
        return None;
    }
    let unit = vec_scale(axis, 1.0 / distance);
    let x = (radius_a * radius_a - radius_b * radius_b + distance * distance) / (2.0 * distance);
    let height_sq = radius_a * radius_a - x * x;
    if height_sq < -1.0e-8 {
        return None;
    }
    let base = vec_add(center_a, vec_scale(unit, x));
    let height = height_sq.max(0.0).sqrt();
    let side = vec_normalize(vec_sub(
        preferred_dir,
        vec_scale(unit, vec_dot(preferred_dir, unit)),
    ));
    let side = if vec_len(side) <= 1.0e-10 {
        perpendicular_unit(unit)
    } else {
        side
    };
    Some(vec_add(base, vec_scale(side, height)))
}

#[allow(clippy::too_many_arguments)]
fn place_confseq_constraint_ring_component_from_attachment(
    component: &ConfSeqBaseRingComponent,
    scaffold: &[[f64; 2]],
    anchor: usize,
    attachment: usize,
    bond_length: f64,
    attachment_dir: [f64; 3],
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
) -> Result<(), ConfSeqBaseConformerError> {
    let scaffold_anchor = [scaffold[anchor][0], scaffold[anchor][1], 0.0];
    let scaffold_attachment = [scaffold[attachment][0], scaffold[attachment][1], 0.0];
    let scaffold_dir = vec_normalize(vec_sub(scaffold_attachment, scaffold_anchor));
    let target_attachment = vec_add(coords[anchor], vec_scale(attachment_dir, bond_length));
    for &atom in &component.atoms {
        let point = [scaffold[atom][0], scaffold[atom][1], 0.0];
        let offset = vec_sub(point, scaffold_attachment);
        let rotated = rotate_vector_between_unit_dirs(offset, scaffold_dir, attachment_dir)?;
        coords[atom] = vec_add(target_attachment, rotated);
        placed[atom] = true;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn confseq_constraint_child_direction(
    _molecule: &Molecule,
    _adjacency: &AdjacencyList,
    _model: &ConfSeqBaseConstraintModel,
    _coords: &[[f64; 3]],
    placed: &[bool],
    center: usize,
    parent: Option<usize>,
    _child: usize,
    parent_axis: [f64; 3],
    child_ord: usize,
    child_count: usize,
    angle: f64,
    _length: f64,
) -> [f64; 3] {
    let fallback = child_direction(parent_axis, child_ord, child_count, angle);
    let Some(parent) = parent else {
        return fallback;
    };
    if !placed[parent] || !placed[center] {
        return fallback;
    }
    fallback
}

fn validate_confseq_base_constraint_coords(
    molecule: &Molecule,
    model: &ConfSeqBaseConstraintModel,
    coords: &[[f64; 3]],
) -> Result<(), ConfSeqBaseConformerError> {
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
            return Err(ConfSeqBaseConformerError::Build(format!(
                "base constraint bond target failed for bond {}: observed={observed:.3}, target={target:.3}",
                bond.id().index()
            )));
        }
    }
    for &(begin, end) in &model.planar_bonds {
        let observed = vec_len(vec_sub(coords[begin], coords[end]));
        if observed <= 1.0e-10 {
            return Err(ConfSeqBaseConformerError::Build(format!(
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
            return Err(ConfSeqBaseConformerError::Build(format!(
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
) -> Result<Molecule, ConfSeqBaseConformerError> {
    let constraints = collect_confseq_base_double_bond_stereo_constraints(&molecule);
    let mut molecule = molecule;
    for (i, j, k, l, target_deg) in constraints {
        molecule = mol_transforms::set_dihedral_deg(molecule, i, j, k, l, target_deg, 0).map_err(
            |err| {
                ConfSeqBaseConformerError::Build(format!(
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
) -> Result<Molecule, ConfSeqBaseConformerError> {
    let constraints = collect_confseq_base_tetrahedral_stereo_constraints(&molecule, adjacency)?;
    if constraints.is_empty() {
        return Ok(molecule);
    }

    let mut coords = molecule
        .conformers_3d()
        .first()
        .ok_or_else(|| {
            ConfSeqBaseConformerError::Build("base conformer has no coordinates".to_string())
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
                return Err(ConfSeqBaseConformerError::UnsupportedTetrahedralStereo {
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
            return Err(ConfSeqBaseConformerError::UnsupportedTetrahedralStereo {
                center: constraint.center,
                reason: "final chiral-volume sign validation failed".to_string(),
            });
        }
    }

    let mut out = molecule;
    let coord_block = out.coordinate_block_mut();
    let Some(conformer) = coord_block.conformers_3d.first_mut() else {
        return Err(ConfSeqBaseConformerError::Build(
            "base conformer has no coordinates".to_string(),
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
) -> Result<Vec<ConfSeqBaseTetrahedralStereoConstraint>, ConfSeqBaseConformerError> {
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
                return Err(ConfSeqBaseConformerError::UnsupportedTetrahedralStereo {
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
) -> Result<Vec<usize>, ConfSeqBaseConformerError> {
    let center = constraint.center;
    if movable_ligand == center
        || bond_between_pair(molecule, sorted_pair(center, movable_ligand)).is_none()
    {
        return Err(ConfSeqBaseConformerError::UnsupportedTetrahedralStereo {
            center,
            reason: "has no explicit movable ligand".to_string(),
        });
    }
    let side =
        connected_side_without_crossing(adjacency, molecule.num_atoms(), movable_ligand, center);
    if side.contains(&center) {
        return Err(ConfSeqBaseConformerError::UnsupportedTetrahedralStereo {
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
) -> Result<(), ConfSeqBaseConformerError> {
    let center = constraint.center;
    let movable_ligand = constraint.ligands[movable_pos];
    let normal = confseq_base_chiral_volume_gradient_for_ligand(coords, constraint, movable_pos);
    let normal_len = vec_len(normal);
    if normal_len <= 1.0e-10 {
        return Err(ConfSeqBaseConformerError::UnsupportedTetrahedralStereo {
            center,
            reason: "has collinear fixed ligands".to_string(),
        });
    }
    let normal = vec_scale(normal, 1.0 / normal_len);
    let center_point = coords[center];
    let old_root = vec_sub(coords[movable_ligand], center_point);
    let root_len = vec_len(old_root);
    if root_len <= 1.0e-10 {
        return Err(ConfSeqBaseConformerError::UnsupportedTetrahedralStereo {
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
        return Err(ConfSeqBaseConformerError::UnsupportedTetrahedralStereo {
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
) -> Result<bool, ConfSeqBaseConformerError> {
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
) -> Result<bool, ConfSeqBaseConformerError> {
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
        Err(ConfSeqBaseConformerError::UnsupportedTetrahedralStereo {
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
) -> Result<(), ConfSeqBaseConformerError> {
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
) -> Result<[f64; 3], ConfSeqBaseConformerError> {
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

fn dihedral_rad_from_points(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let b0 = vec_sub(a, b);
    let b1 = vec_sub(c, b);
    let b2 = vec_sub(d, c);
    let b1 = vec_normalize(b1);
    let v = vec_sub(b0, vec_scale(b1, vec_dot(b0, b1)));
    let w = vec_sub(b2, vec_scale(b1, vec_dot(b2, b1)));
    vec_dot(vec_cross(b1, v), w).atan2(vec_dot(v, w))
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
) -> Result<(), ConfSeqBaseConformerError> {
    if !(3..=8).contains(&ring_atoms.len()) {
        return Err(ConfSeqBaseConformerError::UnsupportedRingSize {
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
            return Err(ConfSeqBaseConformerError::UnsupportedAromaticRingSize {
                ring_size: ring_atoms.len(),
            });
        }
        for atom in ring_atoms {
            let atomic_number = molecule.atoms()[atom.index()].atomic_number();
            if !matches!(atomic_number, 6 | 7 | 8 | 16) {
                return Err(ConfSeqBaseConformerError::UnsupportedRingElement {
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
            return Err(ConfSeqBaseConformerError::UnsupportedRingElement {
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
    // base-conformer tests and common organic heavy-atom bonds. They are not a
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

fn confseq_base_coord_centroid(coords: &[[f64; 3]]) -> [f64; 3] {
    if coords.is_empty() {
        return [0.0, 0.0, 0.0];
    }
    let sum = coords
        .iter()
        .fold([0.0, 0.0, 0.0], |acc, point| vec_add(acc, *point));
    vec_scale(sum, 1.0 / coords.len() as f64)
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
        for (pair, dihedral) in &template.dihedrals_by_pair {
            let angle = parsed.dihedral_angles_by_pair.get(pair).ok_or(
                ConfSeqDecodeError::DihedralTokenCountMismatch {
                    observed: parsed.dihedral_angles_by_pair.len(),
                    expected: template.dihedrals_by_pair.len(),
                },
            )?;
            if template.ring_bond_pairs.contains(pair) {
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
    //   mol,unapplied_dihedrals=apply_dihedrals(mol_no_ring,unapplied_dihedrals[:])
    //   mol,unapplied_dihedrals=apply_dihedrals(mol,unapplied_dihedrals[:])
    //   for begin,end in last_ring_bonds:
    //       ring_mol.AddBond(begin,end,order=Chem.BondType.SINGLE)
    //
    // COSMolKit keeps atom indices stable while removing/re-adding bonds, so the RDKit
    // MolToSmiles-based renumbering step is not needed here.
    let mut builder = molecule.to_builder();
    for (begin, end, _) in last_ring_bonds {
        builder.remove_bond_between_atoms(AtomId::new(*begin), AtomId::new(*end));
    }
    let mut molecule = builder
        .build()
        .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;

    let mut deferred = change_dihedral_for_removed_ring_bonds(last_ring_bonds, unapplied);
    for _ in 0..2 {
        let mut still_deferred = Vec::new();
        for (dihedral, angle) in deferred {
            if dihedral_bonds_exist(&molecule, dihedral) {
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
    for (_, _, spec) in last_ring_bonds {
        builder
            .add_bond(spec.clone())
            .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))?;
    }
    builder
        .build()
        .map_err(|err| ConfSeqDecodeError::MolTransform(err.to_string()))
}

fn change_dihedral_for_removed_ring_bonds(
    last_ring_bonds: &[(usize, usize, BondSpec)],
    dihedrals: Vec<((usize, usize, usize, usize), f64)>,
) -> Vec<((usize, usize, usize, usize), f64)> {
    // ConfSeq source anchor:
    //   if sorted([dihedral[0],dihedral[1]]) in last_ring_bonds:
    //       dihedral=(dihedral[1],dihedral[2],dihedral[3],dihedral[0])
    //   elif sorted([dihedral[2],dihedral[3]]) in last_ring_bonds:
    //       dihedral=(dihedral[2],dihedral[1],dihedral[0],dihedral[3])
    let removed: HashSet<_> = last_ring_bonds
        .iter()
        .map(|(begin, end, _)| sorted_pair(*begin, *end))
        .collect();
    dihedrals
        .into_iter()
        .map(|((i, j, k, l), angle)| {
            if removed.contains(&sorted_pair(i, j)) {
                ((j, k, l, i), angle)
            } else if removed.contains(&sorted_pair(k, l)) {
                ((k, j, i, l), angle)
            } else {
                ((i, j, k, l), angle)
            }
        })
        .collect()
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
    // ConfSeq✔️🔝: COSMolKit generates the same explicit-bond token positions while
    // avoiding Indigo molfile round-trip allocation; atom-pair and token-index maps
    // are produced in one linear scan over all-bonds-explicit SMILES.
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
    let mapping = bond_token_mapping_for_smiles(smiles)?;
    let shared = collect_fully_shared_ring_bonds(molecule)?;
    let mut last_ring_bonds = Vec::new();
    let mut seen = HashSet::new();
    for pair in mapping.ring_closure_pairs {
        if shared.contains(&pair) || !seen.insert(pair) {
            continue;
        }
        if let Some(bond) = bond_between_pair(molecule, pair) {
            last_ring_bonds.push((pair.0, pair.1, bond_spec_from_bond(bond)));
        }
    }
    Ok(last_ring_bonds)
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
            .into_iter()
            .map(|dihedral @ (_, j, k, _)| (sorted_pair(j, k), dihedral))
            .collect();
        let angle_centers = collect_angle_centers(&embedded)?;
        let ring_bond_pairs = collect_ring_bond_pairs(&embedded)?;
        let last_ring_bonds = collect_last_ring_bonds(smiles, &embedded)?;

        Ok(Template {
            molecule: embedded,
            dihedrals_by_pair,
            angle_centers,
            ring_bond_pairs,
            last_ring_bonds,
        })
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
            template_backend: ConfSeqTemplateBackend::BaseConformer,
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

    fn quantile_sorted(values: &[f64], q: f64) -> f64 {
        if values.is_empty() {
            return f64::NAN;
        }
        let idx = ((values.len() - 1) as f64 * q).round() as usize;
        values[idx]
    }

    #[derive(Debug, Clone)]
    struct RigidFragmentRmsdSummary {
        fragment_rmsds: Vec<f64>,
        max_rmsd: Option<f64>,
        worst_fragment_atoms: Vec<usize>,
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
        let heavy_index_by_atom = heavy_index_by_atom(molecule);
        let mut fragment_rmsds = Vec::new();
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
            if fragment_rmsds
                .iter()
                .copied()
                .max_by(f64::total_cmp)
                .map(|current_max| rmsd > current_max)
                .unwrap_or(true)
            {
                worst_fragment_atoms = fragment;
            }
            fragment_rmsds.push(rmsd);
        }
        let max_rmsd = fragment_rmsds.iter().copied().max_by(f64::total_cmp);
        RigidFragmentRmsdSummary {
            fragment_rmsds,
            max_rmsd,
            worst_fragment_atoms,
        }
    }

    fn rigid_heavy_fragments_cutting_confseq_rotors(molecule: &Molecule) -> Vec<Vec<usize>> {
        let cut_bonds: HashSet<_> = collect_single_bond_dihedrals(molecule)
            .into_iter()
            .map(|(_, j, k, _)| sorted_pair(j, k))
            .collect();
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
        model: &ConfSeqBaseConstraintModel,
        fragment_atoms: &[usize],
    ) -> &'static str {
        let fragment_set: HashSet<_> = fragment_atoms.iter().copied().collect();
        for component in &model.ring_components {
            if component
                .atoms
                .iter()
                .all(|atom| fragment_set.contains(atom))
            {
                return if component.planar {
                    "ring_planar"
                } else {
                    "ring_nonplanar"
                };
            }
        }
        let all_planar_like = fragment_atoms.iter().copied().all(|atom_idx| {
            let atom = &molecule.atoms()[atom_idx];
            atom.is_aromatic() || atom.hybridization() == Hybridization::Sp2
        });
        if all_planar_like {
            "planar_pi"
        } else if fragment_atoms.iter().any(|atom| {
            model
                .ring_components
                .iter()
                .any(|component| component.atoms.contains(atom))
        }) {
            "ring_mixed"
        } else {
            "acyclic"
        }
    }

    fn rigid_fragment_local_diagnostic(molecule: &Molecule, fragment_atoms: &[usize]) -> String {
        let fragment_set: HashSet<_> = fragment_atoms.iter().copied().collect();
        let cut_bonds: HashSet<_> = collect_single_bond_dihedrals(molecule)
            .into_iter()
            .map(|(_, j, k, _)| sorted_pair(j, k))
            .collect();
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
            let cache: ConfSeqDgReferenceCache = serde_json::from_str(&raw).unwrap_or_else(|err| {
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
        let reference_batch = decode_confseq_batch_with_options(
            in_smiles_batch,
            td_smiles_batch,
            &reference_options,
            true,
        )
        .expect("reference batch should decode with per-record errors");

        let entries = reference_batch
            .molecules
            .iter()
            .zip(reference_batch.errors.iter())
            .map(|(molecule, error)| ConfSeqDgReferenceCacheEntry {
                heavy_atom_points: molecule.as_ref().map(heavy_atom_points_for_rmsd),
                error: error.as_ref().map(|error| format!("{error:?}")),
            })
            .collect();
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
    fn confseq_base_conformer_preserves_explicit_double_bond_stereo_planes() {
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
    fn confseq_base_conformer_satisfies_tetrahedral_stereo_volume_sign() {
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
                try_build_confseq_base_conformer(&molecule).expect("base conformer should build");
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
                try_build_confseq_base_conformer(&molecule).expect("base conformer should build");

            assert_eq!(embedded.num_atoms(), molecule.num_atoms());
            assert_eq!(embedded.conformers_3d().len(), 1);
            assert!(
                embedded
                    .conformers_3d()
                    .first()
                    .expect("base conformer exists")
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
        let model = build_confseq_base_constraint_model(&molecule, &ring_info)
            .expect("shared ring systems should be modeled as one component");

        assert_eq!(model.ring_components.len(), 1);
        assert!(model.ring_components[0].atoms.len() > ring_info.atom_rings()[0].len());
    }

    #[test]
    fn confseq_base_conformer_initializes_shared_ring_system_scaffolds() {
        for smiles in ["C1CCC2(CC1)CCC2", "C1CC2(CC1)CCNCC2"] {
            let molecule = Molecule::from_smiles(smiles).expect("spiro fixture parses");
            let embedded =
                try_build_confseq_base_conformer(&molecule).expect("base conformer should build");
            let coords = embedded.conformers_3d()[0].coordinates();

            assert!(
                coords
                    .iter()
                    .all(|point| point.iter().all(|value| value.is_finite())),
                "shared ring fixture should produce finite base coordinates: {smiles}"
            );
        }
    }

    #[test]
    fn confseq_base_ring_system_scaffold_uses_rdkit_initial_prior() {
        let molecule = Molecule::from_smiles("c1ccc2ccccc2c1").expect("naphthalene parses");
        let embedded =
            try_build_confseq_base_conformer(&molecule).expect("base conformer should build");
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
            try_build_confseq_base_conformer(&molecule).expect("base conformer should build");
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
            template_backend: ConfSeqTemplateBackend::BaseConformer,
            ..ConfSeqDecodeOptions::default()
        };

        let error =
            decode_confseq_with_options("C 1 C C 2 C C C 1 C 2", "C 1 C C 2 C C C 1 C 2", &options)
                .expect_err("base-conformer decode failure must not fallback to DistGeom");

        assert!(matches!(
            error,
            ConfSeqDecodeError::AngleTokenCountMismatch { .. }
                | ConfSeqDecodeError::DihedralTokenCountMismatch { .. }
                | ConfSeqDecodeError::BaseConformer(_)
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
    fn confseq_base_corpus_pass_rate_snapshot() {
        let corpus_path = std::env::var("COSMOLKIT_CONFSEQ_CORPUS").unwrap_or_else(|_| {
            "/home/wangjingtong/sh4090/confseq_test_strings_100x10.jsonl".to_string()
        });
        let input =
            std::fs::read_to_string(&corpus_path).expect("ConfSeq corpus should be readable");
        let base_options = ConfSeqDecodeOptions {
            optimize_with_uff: false,
            template_backend: ConfSeqTemplateBackend::BaseConformer,
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
        let mut rmsds = Vec::new();
        let mut rigid_fragment_rmsds = Vec::new();
        let mut rigid_fragment_max_rmsds = Vec::new();
        let mut worst_rigid_fragments = Vec::<(f64, usize, Vec<usize>, String)>::new();
        let mut rigid_fragment_pass_01a = 0usize;
        let mut rigid_fragment_pass_02a = 0usize;
        let mut rigid_fragment_pass_03a = 0usize;
        let mut rigid_max_pass_01a = 0usize;
        let mut rigid_max_pass_02a = 0usize;
        let mut rigid_max_pass_03a = 0usize;
        let mut global_fail_rigid_pass = 0usize;
        let mut global_fail_rigid_fail = 0usize;
        let mut base_errors: HashMap<String, usize> = HashMap::new();

        for idx in 0..total {
            let Some(reference_points) = reference_cache.entries[idx].heavy_atom_points.as_ref()
            else {
                continue;
            };
            reference_success += 1;

            match base_batch.molecules[idx].as_ref() {
                Some(base) => {
                    base_success_where_reference_success += 1;
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
                    let parsed = parse_confseq(&in_smiles_batch[idx], &td_smiles_batch[idx])
                        .expect("fixture ConfSeq should parse");
                    let source_mol = Molecule::from_smiles(&parsed.stripped_smiles)
                        .expect("fixture stripped SMILES should parse");
                    let rigid_summary =
                        rigid_fragment_rmsd_summary(&source_mol, reference_points, base);
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
                    rigid_fragment_rmsds.extend(rigid_summary.fragment_rmsds.iter().copied());
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
        rigid_fragment_rmsds.sort_by(|left, right| left.total_cmp(right));
        rigid_fragment_max_rmsds.sort_by(|left, right| left.total_cmp(right));
        let total_pass_rate = pass_03a as f64 / total.max(1) as f64 * 100.0;
        let comparable_pass_rate =
            pass_03a as f64 / base_success_where_reference_success.max(1) as f64 * 100.0;
        let base_coverage =
            base_success_where_reference_success as f64 / reference_success.max(1) as f64 * 100.0;

        eprintln!("confseq_base_corpus path={corpus_path}");
        eprintln!(
            "confseq_base_corpus total={total} reference_success={reference_success} base_success_where_reference_success={base_success_where_reference_success} pass_rmsd_le_0_3a={pass_03a}"
        );
        eprintln!(
            "confseq_base_corpus total_pass_rate={total_pass_rate:.2}% comparable_pass_rate={comparable_pass_rate:.2}% base_coverage_vs_reference={base_coverage:.2}%"
        );
        eprintln!(
            "confseq_base_corpus rmsd p50={:.6} p75={:.6} p90={:.6} p95={:.6} p99={:.6}",
            quantile_sorted(&rmsds, 0.50),
            quantile_sorted(&rmsds, 0.75),
            quantile_sorted(&rmsds, 0.90),
            quantile_sorted(&rmsds, 0.95),
            quantile_sorted(&rmsds, 0.99)
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
        worst_rigid_fragments.sort_by(|left, right| right.0.total_cmp(&left.0));
        for (max_rmsd, idx, atoms, stripped_smiles) in worst_rigid_fragments.into_iter().take(12) {
            let mol = Molecule::from_smiles(&stripped_smiles)
                .expect("fixture stripped SMILES should parse");
            let rings = rings::symmetrize_sssr(&mol).expect("ring perception should work");
            let model = build_confseq_base_constraint_model(&mol, &rings)
                .expect("base constraint model should build for diagnostic");
            let fragment_type = rigid_fragment_type_for_diagnostic(&mol, &model, &atoms);
            let local = rigid_fragment_local_diagnostic(&mol, &atoms);
            eprintln!(
                "confseq_base_corpus worst_rigid_fragment idx={idx} max_rmsd={max_rmsd:.6} type={fragment_type} atoms={atoms:?} {local} stripped_smiles={stripped_smiles}"
            );
        }

        let mut errors: Vec<_> = base_errors.into_iter().collect();
        errors.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
        for (error, count) in errors.into_iter().take(12) {
            eprintln!("confseq_base_corpus base_error count={count} error={error}");
        }
    }

    #[test]
    #[ignore = "local ConfSeq corpus BaseConformer-only coverage snapshot"]
    fn confseq_base_corpus_base_only_snapshot() {
        let corpus_path = std::env::var("COSMOLKIT_CONFSEQ_CORPUS").unwrap_or_else(|_| {
            "/home/wangjingtong/sh4090/confseq_test_strings_100x10.jsonl".to_string()
        });
        let input =
            std::fs::read_to_string(&corpus_path).expect("ConfSeq corpus should be readable");
        let base_options = ConfSeqDecodeOptions {
            optimize_with_uff: false,
            template_backend: ConfSeqTemplateBackend::BaseConformer,
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
            template_backend: ConfSeqTemplateBackend::BaseConformer,
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
            let model = build_confseq_base_constraint_model(&mol, &rings)
                .expect("base constraint model should build for diagnostic");
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
                for (component_idx, component) in model.ring_components.iter().enumerate() {
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
                        "idx={idx} component={component_idx} planar={} atoms=[{}] rings={:?}",
                        component.planar,
                        atom_details.join(","),
                        component.rings
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
                let model = build_confseq_base_constraint_model(&mol, &rings)
                    .expect("base constraint model should build for diagnostic");
                let planar_components = model
                    .ring_components
                    .iter()
                    .filter(|component| component.planar)
                    .count();
                let nonplanar_components = model.ring_components.len() - planar_components;
                let base = base_batch.molecules[idx].as_ref().expect("checked above");
                let coords = conformer_points(base);
                let path14_penalty = confseq_base_path14_distance_penalty(&model, &coords);
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
                    "quantile={quantile:.2} idx={idx} rmsd={rmsd:.6} atoms={} rings={} components={} planar_components={} nonplanar_components={} rot_single={} path14_penalty={path14_penalty:.6} stripped_smiles={stripped}",
                    mol.num_atoms(),
                    rings.num_rings(),
                    model.ring_components.len(),
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
            template_backend: ConfSeqTemplateBackend::BaseConformer,
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
            rows.push((
                full_rmsd,
                line_idx,
                template_rmsd,
                angle_rmsd,
                parsed.stripped_smiles.clone(),
            ));
        }

        rows.sort_by(|left, right| right.0.total_cmp(&left.0));
        for (full_rmsd, line_idx, template_rmsd, angle_rmsd, stripped_smiles) in
            rows.into_iter().take(16)
        {
            eprintln!(
                "confseq_base_stage idx={line_idx} template_rmsd={template_rmsd:.6} angle_only_rmsd={angle_rmsd:.6} full_rmsd={full_rmsd:.6} stripped_smiles={stripped_smiles}"
            );
        }
    }

    #[test]
    #[ignore = "local ConfSeq base candidate diagnostics"]
    fn confseq_base_candidate_diagnostics() {
        let smiles = std::env::var("COSMOLKIT_CONFSEQ_BASE_DIAG_SMILES")
            .unwrap_or_else(|_| "Cc1cc(C)nc(NN2C(=O)CSC2=S)n1".to_string());
        let mol = Molecule::from_smiles(&smiles).expect("diagnostic SMILES should parse");
        let rings = rings::symmetrize_sssr(&mol).expect("ring perception should work");
        let adjacency = AdjacencyList::from_topology(mol.num_atoms(), mol.bonds());
        let model = build_confseq_base_constraint_model(&mol, &rings)
            .expect("base constraint model should build");
        let chiral_constraints =
            collect_confseq_base_tetrahedral_stereo_constraints(&mol, &adjacency)
                .expect("chiral constraints should collect");

        let candidates = [
            (
                "rdkit_scaffold",
                construct_confseq_base_coords_from_rdkit_scaffold(&mol, &adjacency, &model),
            ),
            (
                "local_3d_prior",
                construct_confseq_base_coords_from_local_3d_prior(&mol, &adjacency, &model),
            ),
        ];
        eprintln!("confseq_base_candidate_diag smiles={smiles}");
        for (name, candidate) in candidates {
            match candidate {
                Ok(coords) => {
                    let total =
                        confseq_base_candidate_score(&mol, &model, &chiral_constraints, &coords);
                    let (violations, penalty, max_deficit) =
                        confseq_base_path14_distance_score_details(&model, &coords);
                    let max_abs_z = coords
                        .iter()
                        .map(|point| point[2].abs())
                        .fold(0.0, f64::max);
                    eprintln!(
                        "confseq_base_candidate_diag candidate={name} total_score={total:.6} path14_violations={violations} path14_penalty={penalty:.6} max_path14_deficit={max_deficit:.6} max_abs_z={max_abs_z:.3}"
                    );
                }
                Err(err) => {
                    eprintln!("confseq_base_candidate_diag candidate={name} error={err:?}");
                }
            }
        }
        for (component_idx, component) in model.ring_components.iter().enumerate() {
            eprintln!(
                "confseq_base_candidate_diag component={component_idx} planar={} atoms={:?} rings={:?}",
                component.planar, component.atoms, component.rings
            );
        }
    }

    fn confseq_base_path14_distance_score_details(
        model: &ConfSeqBaseConstraintModel,
        coords: &[[f64; 3]],
    ) -> (usize, f64, f64) {
        let mut violations = 0usize;
        let mut penalty = 0.0;
        let mut max_deficit: f64 = 0.0;
        for prior in &model.path14_distance_priors {
            let (i, _, _, l) = prior.atoms;
            let observed = vec_len(vec_sub(coords[i], coords[l]));
            let deficit = if observed < prior.lower_bound {
                prior.lower_bound - observed
            } else if observed > prior.upper_bound {
                observed - prior.upper_bound
            } else {
                0.0
            };
            if deficit > 0.0 {
                violations += 1;
                penalty += deficit * deficit * 1.5;
                max_deficit = max_deficit.max(deficit);
            }
        }
        (violations, penalty, max_deficit)
    }

    fn planar_bond_like_for_diagnostic(molecule: &Molecule, bond: &Bond) -> bool {
        bond.order() == BondOrder::Double
            || molecule.atoms()[bond.begin().index()].is_aromatic()
            || molecule.atoms()[bond.end().index()].is_aromatic()
            || (molecule.atoms()[bond.begin().index()].hybridization() == Hybridization::Sp2
                && molecule.atoms()[bond.end().index()].hybridization() == Hybridization::Sp2)
    }
}
