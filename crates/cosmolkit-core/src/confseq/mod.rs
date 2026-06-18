use std::collections::{HashMap, HashSet, VecDeque};
use std::f64::consts::PI;

use rayon::prelude::*;
use thiserror::Error;

use crate::chemistry::{mol_transforms, rings};
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
    let geometry = ConfSeqBaseGeometry::new(molecule, &ring_info);

    let coords = if ring_info.num_rings() == 0 {
        build_acyclic_confseq_base_coords(molecule, &adjacency, &geometry)
    } else {
        build_supported_aromatic_ring_system_confseq_base_coords(
            molecule, &adjacency, &ring_info, &geometry,
        )?
    };

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
    for _ in 0..constraints.len().saturating_mul(2).max(1) {
        let mut changed = false;
        for constraint in &constraints {
            let volume = confseq_base_chiral_volume(&coords, constraint);
            if confseq_base_chiral_volume_satisfies_tag(volume, constraint.tag) {
                continue;
            }
            let mut corrected = false;
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
            // Prefer independent ligand-side moves before whole cyclic-side
            // moves. Both candidate classes preserve topology distances and are
            // accepted only after explicit volume revalidation; the ordering is
            // deterministic minimal perturbation when several corrections work.
            candidates.sort_by_key(|(contains_other_ligand, len, movable_pos, _)| {
                (*contains_other_ligand, *len, *movable_pos)
            });
            for (_, _, movable_pos, movable) in candidates {
                let mut trial = coords.clone();
                rotate_movable_side_to_chiral_volume_sign(
                    &mut trial,
                    &movable,
                    constraint,
                    movable_pos,
                )?;
                if confseq_base_chiral_volume_satisfies_tag(
                    confseq_base_chiral_volume(&trial, constraint),
                    constraint.tag,
                ) {
                    coords = trial;
                    corrected = true;
                    break;
                }
            }
            if !corrected {
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

fn build_acyclic_confseq_base_coords(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    geometry: &ConfSeqBaseGeometry,
) -> Vec<[f64; 3]> {
    // Current subset: connected acyclic heavy-atom graphs are placed as a
    // deterministic tree. This deliberately solves only the ConfSeq base-frame
    // requirement: reasonable bond lengths and local angles before token-driven
    // dihedral/angle application. It is not a conformer search and does not
    // model non-bonded interactions.
    let mut coords = vec![[0.0; 3]; molecule.num_atoms()];
    let mut placed = vec![false; molecule.num_atoms()];
    place_confseq_tree(
        molecule,
        adjacency,
        0,
        None,
        [0.0, 0.0, 0.0],
        &mut coords,
        &mut placed,
        geometry,
    );
    coords
}

fn build_supported_aromatic_ring_system_confseq_base_coords(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    ring_info: &rings::RingInfo,
    geometry: &ConfSeqBaseGeometry,
) -> Result<Vec<[f64; 3]>, ConfSeqBaseConformerError> {
    let atom_to_ring = confseq_base_atom_to_primary_ring(molecule, ring_info)?;
    let ring_bonds = confseq_base_ring_bond_set(ring_info)?;
    let mut coords = vec![[0.0; 3]; molecule.num_atoms()];
    let mut placed = vec![false; molecule.num_atoms()];
    let mut ring_placed = vec![false; ring_info.num_rings()];

    match atom_to_ring[0] {
        Some(ring_idx) => {
            place_confseq_ring_component(
                molecule,
                ring_info,
                ring_idx,
                None,
                &mut coords,
                &mut placed,
                &mut ring_placed,
                geometry,
            )?;
            propagate_confseq_base_from_placed_rings(
                molecule,
                adjacency,
                ring_info,
                &atom_to_ring,
                &ring_bonds,
                &mut coords,
                &mut placed,
                &mut ring_placed,
                geometry,
            )?;
        }
        None => {
            place_confseq_tree_with_ring_blocks(
                molecule,
                adjacency,
                ring_info,
                &atom_to_ring,
                &ring_bonds,
                0,
                None,
                [0.0, 0.0, 0.0],
                &mut coords,
                &mut placed,
                &mut ring_placed,
                geometry,
            )?;
        }
    }

    if placed.iter().any(|value| !*value) {
        return Err(ConfSeqBaseConformerError::PlacementLeftAtomsUnplaced);
    }

    Ok(coords)
}

fn confseq_base_atom_to_primary_ring(
    molecule: &Molecule,
    ring_info: &rings::RingInfo,
) -> Result<Vec<Option<usize>>, ConfSeqBaseConformerError> {
    let mut atom_to_ring = vec![None; molecule.num_atoms()];
    for (ring_idx, ring) in ring_info.atom_rings().iter().enumerate() {
        validate_supported_confseq_base_ring(
            molecule,
            ring_idx,
            ring,
            &ring_info.bond_rings()[ring_idx],
        )?;
        for atom in ring {
            if atom_to_ring[atom.index()].is_none() {
                atom_to_ring[atom.index()] = Some(ring_idx);
            }
        }
    }
    validate_confseq_base_ring_sharing(molecule, ring_info)?;
    Ok(atom_to_ring)
}

fn confseq_base_ring_bond_set(
    ring_info: &rings::RingInfo,
) -> Result<HashSet<usize>, ConfSeqBaseConformerError> {
    let mut ring_bonds = HashSet::new();
    for ring in ring_info.bond_rings() {
        for bond in ring {
            ring_bonds.insert(bond.index());
        }
    }
    Ok(ring_bonds)
}

fn validate_confseq_base_ring_sharing(
    _molecule: &Molecule,
    ring_info: &rings::RingInfo,
) -> Result<(), ConfSeqBaseConformerError> {
    let mut fusion_graph = vec![Vec::new(); ring_info.num_rings()];
    for left in 0..ring_info.num_rings() {
        for right in left + 1..ring_info.num_rings() {
            let shared_atoms = ring_info.atom_rings()[left]
                .iter()
                .filter(|atom| ring_info.atom_rings()[right].contains(atom))
                .count();
            let shared_bonds = shared_bond_ids_between_rings(ring_info, left, right).len();
            if shared_atoms == 0 && shared_bonds == 0 {
                continue;
            }
            let edge_fused = shared_atoms == 2 && shared_bonds == 1;
            let spiro_fused = shared_atoms == 1 && shared_bonds == 0;
            if !edge_fused && !spiro_fused {
                return Err(ConfSeqBaseConformerError::UnsupportedRingFusion {
                    left,
                    right,
                    shared_atoms,
                    shared_bonds,
                });
            }
            fusion_graph[left].push(right);
            fusion_graph[right].push(left);
        }
    }
    if !confseq_base_fusion_graph_is_forest(&fusion_graph) {
        return Err(ConfSeqBaseConformerError::UnsupportedClosedRingFusion);
    }
    Ok(())
}

fn confseq_base_fusion_graph_is_forest(fusion_graph: &[Vec<usize>]) -> bool {
    // Current fused-ring subset is deliberately a forest of edge fusions and
    // single-atom spiro joins. These can be propagated from one already-placed
    // ring constraint at a time. Closed pericondensed or polyspiro graphs remain
    // rejected until solved as one global constrained component.
    fn visit(node: usize, parent: Option<usize>, graph: &[Vec<usize>], seen: &mut [bool]) -> bool {
        seen[node] = true;
        for &next in &graph[node] {
            if Some(next) == parent {
                continue;
            }
            if seen[next] || !visit(next, Some(node), graph, seen) {
                return false;
            }
        }
        true
    }

    let mut seen = vec![false; fusion_graph.len()];
    for node in 0..fusion_graph.len() {
        if !seen[node] && !visit(node, None, fusion_graph, &mut seen) {
            return false;
        }
    }
    true
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

fn shared_atom_ids_between_rings(
    ring_info: &rings::RingInfo,
    left: usize,
    right: usize,
) -> Vec<AtomId> {
    ring_info.atom_rings()[left]
        .iter()
        .copied()
        .filter(|atom| ring_info.atom_rings()[right].contains(atom))
        .collect()
}

fn place_confseq_ring_component(
    molecule: &Molecule,
    ring_info: &rings::RingInfo,
    ring_idx: usize,
    anchor: Option<(usize, [f64; 3], [f64; 3])>,
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
    ring_placed: &mut [bool],
    geometry: &ConfSeqBaseGeometry,
) -> Result<(), ConfSeqBaseConformerError> {
    place_confseq_ring_block(
        molecule,
        ring_info,
        ring_idx,
        anchor,
        coords,
        placed,
        ring_placed,
        geometry,
    )?;

    loop {
        let mut progressed = false;
        for candidate in 0..ring_info.num_rings() {
            if ring_placed[candidate] {
                continue;
            }
            let Some(anchor_ring) = (0..ring_info.num_rings()).find(|placed_ring| {
                ring_placed[*placed_ring]
                    && shared_bond_ids_between_rings(ring_info, *placed_ring, candidate).len() == 1
                    && shared_atom_ids_between_rings(ring_info, *placed_ring, candidate).len() == 2
            }) else {
                let Some(anchor_ring) = (0..ring_info.num_rings()).find(|placed_ring| {
                    ring_placed[*placed_ring]
                        && shared_bond_ids_between_rings(ring_info, *placed_ring, candidate)
                            .is_empty()
                        && shared_atom_ids_between_rings(ring_info, *placed_ring, candidate).len()
                            == 1
                }) else {
                    continue;
                };
                place_confseq_spiro_ring_on_shared_atom(
                    molecule,
                    ring_info,
                    anchor_ring,
                    candidate,
                    coords,
                    placed,
                    ring_placed,
                    geometry,
                )?;
                progressed = true;
                continue;
            };
            place_confseq_fused_ring_on_shared_bond(
                molecule,
                ring_info,
                anchor_ring,
                candidate,
                coords,
                placed,
                ring_placed,
                geometry,
            )?;
            progressed = true;
        }
        if !progressed {
            break;
        }
    }
    Ok(())
}

fn place_confseq_ring_block(
    molecule: &Molecule,
    ring_info: &rings::RingInfo,
    ring_idx: usize,
    anchor: Option<(usize, [f64; 3], [f64; 3])>,
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
    ring_placed: &mut [bool],
    geometry: &ConfSeqBaseGeometry,
) -> Result<(), ConfSeqBaseConformerError> {
    if ring_placed[ring_idx] {
        return Ok(());
    }
    let ring = &ring_info.atom_rings()[ring_idx];
    // ConfSeq-base ring block: each supported 3-8 member heavy-atom ring is
    // initially placed as a deterministic cyclic polygon using local bond
    // lengths. This is only a scaffold for subsequent ConfSeq angle/dihedral
    // application; it is not a low-energy ring-conformation search. Spiro,
    // bridged, and closed pericondensed systems remain unsupported because
    // they require solving shared constraints as a global component.
    let ring_points = confseq_base_ring_local_points(molecule, ring, geometry)?;

    match anchor {
        Some((anchor_atom, anchor_point, radial_direction)) => {
            let Some(anchor_pos) = ring.iter().position(|atom| atom.index() == anchor_atom) else {
                return Err(ConfSeqBaseConformerError::Build(
                    "ring anchor atom is not a member of the target ring".to_string(),
                ));
            };
            let local_anchor = ring_points[anchor_pos];
            let local_u = confseq_base_ring_substituent_direction_from_points(
                molecule,
                ring,
                anchor_atom,
                &ring_points,
            );
            let local_v = perpendicular_unit(local_u);
            let local_w = vec_normalize(vec_cross(local_u, local_v));
            let world_u = vec_normalize(radial_direction);
            let world_v = perpendicular_unit(world_u);
            let world_w = vec_normalize(vec_cross(world_u, world_v));
            for (atom, point) in ring.iter().zip(ring_points) {
                let rel = vec_sub(point, local_anchor);
                let x = vec_dot(rel, local_u);
                let y = vec_dot(rel, local_v);
                let z = vec_dot(rel, local_w);
                coords[atom.index()] = vec_add(
                    anchor_point,
                    vec_add(
                        vec_scale(world_u, x),
                        vec_add(vec_scale(world_v, y), vec_scale(world_w, z)),
                    ),
                );
                placed[atom.index()] = true;
            }
        }
        None => {
            for (atom, point) in ring.iter().zip(ring_points) {
                coords[atom.index()] = point;
                placed[atom.index()] = true;
            }
        }
    }
    ring_placed[ring_idx] = true;
    Ok(())
}

fn place_confseq_spiro_ring_on_shared_atom(
    molecule: &Molecule,
    ring_info: &rings::RingInfo,
    placed_ring_idx: usize,
    candidate_ring_idx: usize,
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
    ring_placed: &mut [bool],
    geometry: &ConfSeqBaseGeometry,
) -> Result<(), ConfSeqBaseConformerError> {
    let shared = shared_atom_ids_between_rings(ring_info, placed_ring_idx, candidate_ring_idx);
    if shared.len() != 1
        || !shared_bond_ids_between_rings(ring_info, placed_ring_idx, candidate_ring_idx).is_empty()
    {
        return Err(ConfSeqBaseConformerError::UnsupportedRingFusion {
            left: placed_ring_idx,
            right: candidate_ring_idx,
            shared_atoms: shared.len(),
            shared_bonds: shared_bond_ids_between_rings(
                ring_info,
                placed_ring_idx,
                candidate_ring_idx,
            )
            .len(),
        });
    }
    let anchor = shared[0].index();
    if !placed[anchor] {
        return Err(ConfSeqBaseConformerError::Build(
            "shared spiro atom is not placed".to_string(),
        ));
    }
    let normal = confseq_base_ring_plane_normal(ring_info, placed_ring_idx, coords);
    place_confseq_ring_block(
        molecule,
        ring_info,
        candidate_ring_idx,
        Some((anchor, coords[anchor], normal)),
        coords,
        placed,
        ring_placed,
        geometry,
    )
}

fn place_confseq_fused_ring_on_shared_bond(
    molecule: &Molecule,
    ring_info: &rings::RingInfo,
    placed_ring_idx: usize,
    candidate_ring_idx: usize,
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
    ring_placed: &mut [bool],
    geometry: &ConfSeqBaseGeometry,
) -> Result<(), ConfSeqBaseConformerError> {
    let shared = shared_bond_ids_between_rings(ring_info, placed_ring_idx, candidate_ring_idx);
    if shared.len() != 1 {
        return Err(ConfSeqBaseConformerError::UnsupportedRingFusion {
            left: placed_ring_idx,
            right: candidate_ring_idx,
            shared_atoms: shared_atom_ids_between_rings(
                ring_info,
                placed_ring_idx,
                candidate_ring_idx,
            )
            .len(),
            shared_bonds: shared.len(),
        });
    }
    let shared_bond = &molecule.bonds()[shared[0].index()];
    let begin = shared_bond.begin().index();
    let end = shared_bond.end().index();
    if !placed[begin] || !placed[end] {
        return Err(ConfSeqBaseConformerError::Build(
            "shared fused-ring bond is not placed".to_string(),
        ));
    }

    let ring = &ring_info.atom_rings()[candidate_ring_idx];
    let ring_points = confseq_base_ring_local_points(molecule, ring, geometry)?;
    let Some(begin_pos) = ring.iter().position(|atom| atom.index() == begin) else {
        return Err(ConfSeqBaseConformerError::Build(
            "shared fused-ring begin atom is not in candidate ring".to_string(),
        ));
    };
    let Some(end_pos) = ring.iter().position(|atom| atom.index() == end) else {
        return Err(ConfSeqBaseConformerError::Build(
            "shared fused-ring end atom is not in candidate ring".to_string(),
        ));
    };

    let local_begin = ring_points[begin_pos];
    let local_end = ring_points[end_pos];
    let local_u = vec_normalize(vec_sub(local_end, local_begin));
    let local_centroid = centroid_3d(&ring_points);
    let local_centroid_from_bond = vec_sub(
        local_centroid,
        vec_scale(vec_add(local_begin, local_end), 0.5),
    );
    let local_v = vec_normalize(vec_sub(
        local_centroid_from_bond,
        vec_scale(local_u, vec_dot(local_centroid_from_bond, local_u)),
    ));
    let local_v = if vec_len(local_v) <= 1.0e-12 {
        perpendicular_unit(local_u)
    } else {
        local_v
    };
    let local_w = vec_normalize(vec_cross(local_u, local_v));
    let local_centroid_side = vec_dot(vec_sub(local_centroid, local_begin), local_v);
    let world_begin = coords[begin];
    let world_end = coords[end];
    let world_u = vec_normalize(vec_sub(world_end, world_begin));
    let normal = confseq_base_ring_plane_normal(ring_info, placed_ring_idx, coords);
    let mut world_v = vec_normalize(vec_cross(normal, world_u));
    let world_w = vec_normalize(vec_cross(world_u, world_v));
    let existing_centroid = confseq_base_ring_centroid(ring_info, placed_ring_idx, coords);
    let midpoint = vec_scale(vec_add(world_begin, world_end), 0.5);
    let existing_centroid_side = vec_dot(vec_sub(existing_centroid, midpoint), world_v);
    if existing_centroid_side * local_centroid_side > 0.0 {
        world_v = vec_scale(world_v, -1.0);
    }

    for (atom, point) in ring.iter().zip(ring_points) {
        let rel = vec_sub(point, local_begin);
        let x = vec_dot(rel, local_u);
        let y = vec_dot(rel, local_v);
        let z = vec_dot(rel, local_w);
        coords[atom.index()] = vec_add(
            world_begin,
            vec_add(
                vec_scale(world_u, x),
                vec_add(vec_scale(world_v, y), vec_scale(world_w, z)),
            ),
        );
        placed[atom.index()] = true;
    }
    ring_placed[candidate_ring_idx] = true;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn propagate_confseq_base_from_placed_rings(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    ring_info: &rings::RingInfo,
    atom_to_ring: &[Option<usize>],
    ring_bonds: &HashSet<usize>,
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
    ring_placed: &mut [bool],
    geometry: &ConfSeqBaseGeometry,
) -> Result<(), ConfSeqBaseConformerError> {
    for ring_idx in 0..ring_info.num_rings() {
        if !ring_placed[ring_idx] {
            continue;
        }
        for atom in &ring_info.atom_rings()[ring_idx] {
            propagate_confseq_base_from_atom(
                molecule,
                adjacency,
                ring_info,
                atom_to_ring,
                ring_bonds,
                atom.index(),
                None,
                coords,
                placed,
                ring_placed,
                geometry,
            )?;
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn place_confseq_tree_with_ring_blocks(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    ring_info: &rings::RingInfo,
    atom_to_ring: &[Option<usize>],
    ring_bonds: &HashSet<usize>,
    atom_idx: usize,
    parent: Option<usize>,
    point: [f64; 3],
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
    ring_placed: &mut [bool],
    geometry: &ConfSeqBaseGeometry,
) -> Result<(), ConfSeqBaseConformerError> {
    coords[atom_idx] = point;
    placed[atom_idx] = true;
    propagate_confseq_base_from_atom(
        molecule,
        adjacency,
        ring_info,
        atom_to_ring,
        ring_bonds,
        atom_idx,
        parent,
        coords,
        placed,
        ring_placed,
        geometry,
    )
}

#[allow(clippy::too_many_arguments)]
fn propagate_confseq_base_from_atom(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    ring_info: &rings::RingInfo,
    atom_to_ring: &[Option<usize>],
    ring_bonds: &HashSet<usize>,
    atom_idx: usize,
    parent: Option<usize>,
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
    ring_placed: &mut [bool],
    geometry: &ConfSeqBaseGeometry,
) -> Result<(), ConfSeqBaseConformerError> {
    let base_axis = parent
        .map(|parent| vec_normalize(vec_sub(coords[atom_idx], coords[parent])))
        .unwrap_or([1.0, 0.0, 0.0]);
    let children: Vec<_> = adjacency
        .neighbors_of(atom_idx)
        .iter()
        .filter(|nbr| Some(nbr.atom_index) != parent)
        .filter(|nbr| !ring_bonds.contains(&nbr.bond.index()))
        .filter(|nbr| !placed[nbr.atom_index])
        .collect();
    let child_count = children.len();
    for (child_ord, nbr) in children.into_iter().enumerate() {
        let bond = &molecule.bonds()[nbr.bond.index()];
        let length = geometry.bond_length(bond);
        let dir = if let Some(ring_idx) = atom_to_ring[atom_idx] {
            let directions = confseq_base_ring_substituent_directions(
                molecule,
                ring_info,
                ring_idx,
                atom_idx,
                coords,
                child_count,
            );
            directions.get(child_ord).copied().unwrap_or_else(|| {
                confseq_base_ring_substituent_direction(
                    molecule, ring_info, ring_idx, atom_idx, coords,
                )
            })
        } else {
            child_direction(
                base_axis,
                child_ord,
                child_count,
                confseq_base_local_angle_rad(molecule, atom_idx),
            )
        };
        let child_point = vec_add(coords[atom_idx], vec_scale(dir, length));
        if let Some(child_ring_idx) = atom_to_ring[nbr.atom_index] {
            place_confseq_ring_component(
                molecule,
                ring_info,
                child_ring_idx,
                Some((nbr.atom_index, child_point, vec_scale(dir, -1.0))),
                coords,
                placed,
                ring_placed,
                geometry,
            )?;
            propagate_confseq_base_from_placed_rings(
                molecule,
                adjacency,
                ring_info,
                atom_to_ring,
                ring_bonds,
                coords,
                placed,
                ring_placed,
                geometry,
            )?;
        } else {
            place_confseq_tree_with_ring_blocks(
                molecule,
                adjacency,
                ring_info,
                atom_to_ring,
                ring_bonds,
                nbr.atom_index,
                Some(atom_idx),
                child_point,
                coords,
                placed,
                ring_placed,
                geometry,
            )?;
        }
    }
    Ok(())
}

fn confseq_base_ring_substituent_direction(
    molecule: &Molecule,
    ring_info: &rings::RingInfo,
    ring_idx: usize,
    atom_idx: usize,
    coords: &[[f64; 3]],
) -> [f64; 3] {
    let ring = &ring_info.atom_rings()[ring_idx];
    let points: Vec<_> = ring.iter().map(|atom| coords[atom.index()]).collect();
    confseq_base_ring_substituent_direction_from_points(molecule, ring, atom_idx, &points)
}

fn confseq_base_ring_substituent_directions(
    molecule: &Molecule,
    ring_info: &rings::RingInfo,
    ring_idx: usize,
    atom_idx: usize,
    coords: &[[f64; 3]],
    child_count: usize,
) -> Vec<[f64; 3]> {
    let ring = &ring_info.atom_rings()[ring_idx];
    let points: Vec<_> = ring.iter().map(|atom| coords[atom.index()]).collect();
    confseq_base_ring_substituent_directions_from_points(
        molecule,
        ring,
        atom_idx,
        &points,
        child_count,
    )
}

fn confseq_base_ring_substituent_direction_from_points(
    molecule: &Molecule,
    ring: &[AtomId],
    atom_idx: usize,
    points: &[[f64; 3]],
) -> [f64; 3] {
    let Some(pos) = ring.iter().position(|atom| atom.index() == atom_idx) else {
        return [1.0, 0.0, 0.0];
    };
    if ring.len() < 3 || points.len() != ring.len() {
        return [1.0, 0.0, 0.0];
    }

    let center = points[pos];
    let prev = points[(pos + ring.len() - 1) % ring.len()];
    let next = points[(pos + 1) % ring.len()];
    let to_prev = vec_normalize(vec_sub(prev, center));
    let to_next = vec_normalize(vec_sub(next, center));
    let bisector = vec_add(to_prev, to_next);
    if vec_len(bisector) <= 1.0e-12 {
        let centroid = centroid_3d(points);
        return vec_normalize(vec_sub(center, centroid));
    }

    let atom = &molecule.atoms()[atom_idx];
    if atom.is_aromatic() || atom.hybridization() == Hybridization::Sp2 {
        let centroid = centroid_3d(points);
        let first = vec_normalize(bisector);
        let second = vec_scale(first, -1.0);
        let radial = vec_normalize(vec_sub(center, centroid));
        return if vec_dot(first, radial) >= vec_dot(second, radial) {
            first
        } else {
            second
        };
    }

    // Choose the external unit vector whose angles to both ring neighbors match
    // the center's local valence angle when geometrically possible. This keeps
    // ring substituents constrained by bond lengths and local angles without
    // invoking a force field or conformer search.
    let internal_angle = angle_rad(prev, center, next);
    let half_cos = (0.5 * internal_angle).cos();
    let target_cos = confseq_base_local_angle_rad(molecule, atom_idx).cos();
    if half_cos.abs() <= 1.0e-12 {
        let centroid = centroid_3d(points);
        return vec_normalize(vec_sub(center, centroid));
    }
    let in_bisector_scale = (target_cos / half_cos).clamp(-1.0, 1.0);
    let in_plane = vec_scale(vec_normalize(bisector), in_bisector_scale);
    let normal = vec_cross(to_prev, to_next);
    let out_of_plane = (1.0 - in_bisector_scale * in_bisector_scale)
        .max(0.0)
        .sqrt();
    if vec_len(normal) <= 1.0e-12 || out_of_plane <= 1.0e-12 {
        return vec_normalize(in_plane);
    }

    let normal = vec_normalize(normal);
    let first = vec_normalize(vec_add(in_plane, vec_scale(normal, out_of_plane)));
    let second = vec_normalize(vec_add(in_plane, vec_scale(normal, -out_of_plane)));
    let centroid = centroid_3d(points);
    let radial = vec_normalize(vec_sub(center, centroid));
    if vec_dot(first, radial) >= vec_dot(second, radial) {
        first
    } else {
        second
    }
}

fn confseq_base_ring_substituent_directions_from_points(
    molecule: &Molecule,
    ring: &[AtomId],
    atom_idx: usize,
    points: &[[f64; 3]],
    child_count: usize,
) -> Vec<[f64; 3]> {
    if child_count == 0 {
        return Vec::new();
    }
    if child_count == 1 {
        return vec![confseq_base_ring_substituent_direction_from_points(
            molecule, ring, atom_idx, points,
        )];
    }
    let repeated_single = || {
        vec![
            confseq_base_ring_substituent_direction_from_points(molecule, ring, atom_idx, points);
            child_count
        ]
    };
    if child_count != 2 {
        return repeated_single();
    }

    let Some(pos) = ring.iter().position(|atom| atom.index() == atom_idx) else {
        return repeated_single();
    };
    if ring.len() < 3 || points.len() != ring.len() {
        return repeated_single();
    }

    let center = points[pos];
    let prev = points[(pos + ring.len() - 1) % ring.len()];
    let next = points[(pos + 1) % ring.len()];
    let to_prev = vec_normalize(vec_sub(prev, center));
    let to_next = vec_normalize(vec_sub(next, center));
    let bisector = vec_add(to_prev, to_next);
    if vec_len(bisector) <= 1.0e-12 {
        return repeated_single();
    }

    let atom = &molecule.atoms()[atom_idx];
    if atom.is_aromatic() || atom.hybridization() == Hybridization::Sp2 {
        return repeated_single();
    }

    let internal_angle = angle_rad(prev, center, next);
    let half_cos = (0.5 * internal_angle).cos();
    let target_cos = confseq_base_local_angle_rad(molecule, atom_idx).cos();
    if half_cos.abs() <= 1.0e-12 {
        return repeated_single();
    }
    let in_bisector_scale = (target_cos / half_cos).clamp(-1.0, 1.0);
    let in_plane = vec_scale(vec_normalize(bisector), in_bisector_scale);
    let normal = vec_cross(to_prev, to_next);
    let out_of_plane = (1.0 - in_bisector_scale * in_bisector_scale)
        .max(0.0)
        .sqrt();
    if vec_len(normal) <= 1.0e-12 || out_of_plane <= 1.0e-12 {
        return repeated_single();
    }

    let normal = vec_normalize(normal);
    let first = vec_normalize(vec_add(in_plane, vec_scale(normal, out_of_plane)));
    let second = vec_normalize(vec_add(in_plane, vec_scale(normal, -out_of_plane)));
    // Analytic paired solution for two external substituents on an sp3 ring
    // atom with two ring neighbors: both directions satisfy the same target
    // angle to each ring bond and occupy opposite sides of the local ring
    // plane. This avoids the previous invalid collapse where both children
    // reused the single-substituent direction.
    vec![first, second]
}

fn confseq_base_ring_centroid(
    ring_info: &rings::RingInfo,
    ring_idx: usize,
    coords: &[[f64; 3]],
) -> [f64; 3] {
    let ring = &ring_info.atom_rings()[ring_idx];
    let sum = ring
        .iter()
        .map(|atom| coords[atom.index()])
        .fold([0.0; 3], vec_add);
    vec_scale(sum, 1.0 / ring.len() as f64)
}

fn confseq_base_ring_plane_normal(
    ring_info: &rings::RingInfo,
    ring_idx: usize,
    coords: &[[f64; 3]],
) -> [f64; 3] {
    let ring = &ring_info.atom_rings()[ring_idx];
    if ring.len() < 3 {
        return [0.0, 0.0, 1.0];
    }
    let origin = coords[ring[0].index()];
    for left in 1..ring.len() {
        for right in left + 1..ring.len() {
            let normal = vec_cross(
                vec_sub(coords[ring[left].index()], origin),
                vec_sub(coords[ring[right].index()], origin),
            );
            if vec_len(normal) > 1.0e-12 {
                return vec_normalize(normal);
            }
        }
    }
    [0.0, 0.0, 1.0]
}

fn confseq_base_ring_side_lengths(
    molecule: &Molecule,
    ring: &[AtomId],
    geometry: &ConfSeqBaseGeometry,
) -> Result<Vec<f64>, ConfSeqBaseConformerError> {
    let mut lengths = Vec::with_capacity(ring.len());
    for pos in 0..ring.len() {
        let begin = ring[pos].index();
        let end = ring[(pos + 1) % ring.len()].index();
        let Some(bond) = bond_between_pair(molecule, sorted_pair(begin, end)) else {
            return Err(ConfSeqBaseConformerError::Build(
                "ring atom order does not map to adjacent ring bonds".to_string(),
            ));
        };
        lengths.push(geometry.bond_length(bond));
    }
    Ok(lengths)
}

fn confseq_base_ring_local_points(
    molecule: &Molecule,
    ring: &[AtomId],
    geometry: &ConfSeqBaseGeometry,
) -> Result<Vec<[f64; 3]>, ConfSeqBaseConformerError> {
    let side_lengths = confseq_base_ring_side_lengths(molecule, ring, geometry)?;
    let points = if is_confseq_base_saturated_five_membered_ring(molecule, ring)? {
        saturated_five_membered_envelope_points(molecule, ring, &side_lengths)?
    } else if is_confseq_base_puckerable_six_membered_ring(molecule, ring)? {
        saturated_six_membered_chair_points(molecule, ring, &side_lengths)?
    } else {
        planar_cyclic_polygon_points(&side_lengths).map(|points| {
            points
                .into_iter()
                .map(|point| [point[0], point[1], 0.0])
                .collect()
        })?
    };
    Ok(orient_ring_local_points_for_embedded_tetrahedral_stereo(
        molecule, ring, points,
    ))
}

fn orient_ring_local_points_for_embedded_tetrahedral_stereo(
    molecule: &Molecule,
    ring: &[AtomId],
    points: Vec<[f64; 3]>,
) -> Vec<[f64; 3]> {
    if ring.len() < 3 || points.iter().all(|point| point[2].abs() <= 1.0e-12) {
        return points;
    }
    let constraints = collect_ring_local_tetrahedral_constraints(molecule, ring);
    if constraints.is_empty() {
        return points;
    }

    let mirrored: Vec<_> = points
        .iter()
        .map(|point| [point[0], point[1], -point[2]])
        .collect();
    let original_unsatisfied =
        count_unsatisfied_ring_local_tetrahedral_constraints(&points, &constraints);
    let mirrored_unsatisfied =
        count_unsatisfied_ring_local_tetrahedral_constraints(&mirrored, &constraints);
    if mirrored_unsatisfied < original_unsatisfied {
        mirrored
    } else {
        points
    }
}

#[derive(Debug, Clone)]
struct RingLocalTetrahedralConstraint {
    ligands: [RingLocalTetrahedralPoint; 4],
    tag: ChiralTag,
}

#[derive(Debug, Clone, Copy)]
enum RingLocalTetrahedralPoint {
    Ring(usize),
    Substituent(usize),
}

fn collect_ring_local_tetrahedral_constraints(
    molecule: &Molecule,
    ring: &[AtomId],
) -> Vec<RingLocalTetrahedralConstraint> {
    let mut constraints = Vec::new();
    for (center_pos, atom) in ring.iter().enumerate() {
        let center = atom.index();
        let tag = molecule.atoms()[center].chiral_tag();
        if !matches!(tag, ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw) {
            continue;
        }
        let mut ligand_points = Vec::new();
        for bond in molecule.bonds() {
            let begin = bond.begin().index();
            let end = bond.end().index();
            let neighbor = if begin == center {
                end
            } else if end == center {
                begin
            } else {
                continue;
            };
            if let Some(pos) = ring.iter().position(|atom| atom.index() == neighbor) {
                ligand_points.push(RingLocalTetrahedralPoint::Ring(pos));
            } else {
                ligand_points.push(RingLocalTetrahedralPoint::Substituent(center_pos));
            };
        }
        let ligands = match ligand_points.as_slice() {
            [a, b, c] => [*a, *b, *c, RingLocalTetrahedralPoint::Ring(center_pos)],
            [a, b, c, d] => [*a, *b, *c, *d],
            _ => continue,
        };
        constraints.push(RingLocalTetrahedralConstraint { ligands, tag });
    }
    constraints
}

fn count_unsatisfied_ring_local_tetrahedral_constraints(
    points: &[[f64; 3]],
    constraints: &[RingLocalTetrahedralConstraint],
) -> usize {
    constraints
        .iter()
        .filter(|constraint| {
            !confseq_base_chiral_volume_satisfies_tag(
                ring_local_tetrahedral_volume(points, constraint),
                constraint.tag,
            )
        })
        .count()
}

fn ring_local_tetrahedral_volume(
    points: &[[f64; 3]],
    constraint: &RingLocalTetrahedralConstraint,
) -> f64 {
    let [a, b, c, d] = constraint.ligands;
    let anchor = ring_local_tetrahedral_point(points, d);
    let v1 = vec_sub(ring_local_tetrahedral_point(points, a), anchor);
    let v2 = vec_sub(ring_local_tetrahedral_point(points, b), anchor);
    let v3 = vec_sub(ring_local_tetrahedral_point(points, c), anchor);
    vec_dot(v1, vec_cross(v2, v3))
}

fn ring_local_tetrahedral_point(points: &[[f64; 3]], point: RingLocalTetrahedralPoint) -> [f64; 3] {
    match point {
        RingLocalTetrahedralPoint::Ring(pos) => points[pos],
        RingLocalTetrahedralPoint::Substituent(pos) => {
            let center = points[pos];
            let prev = points[(pos + points.len() - 1) % points.len()];
            let next = points[(pos + 1) % points.len()];
            let to_prev = vec_normalize(vec_sub(prev, center));
            let to_next = vec_normalize(vec_sub(next, center));
            let bisector = vec_add(to_prev, to_next);
            let direction = if vec_len(bisector) <= 1.0e-12 {
                vec_normalize(vec_sub(center, centroid_3d(points)))
            } else {
                vec_scale(vec_normalize(bisector), -1.0)
            };
            vec_add(center, direction)
        }
    }
}

fn is_confseq_base_saturated_five_membered_ring(
    molecule: &Molecule,
    ring: &[AtomId],
) -> Result<bool, ConfSeqBaseConformerError> {
    if ring.len() != 5 {
        return Ok(false);
    }
    for atom in ring {
        let atom = &molecule.atoms()[atom.index()];
        if atom.is_aromatic()
            || atom.hybridization() != Hybridization::Sp3
            || !matches!(atom.atomic_number(), 6 | 7 | 8 | 15 | 16)
        {
            return Ok(false);
        }
    }
    for pos in 0..ring.len() {
        let begin = ring[pos].index();
        let end = ring[(pos + 1) % ring.len()].index();
        let Some(bond) = bond_between_pair(molecule, sorted_pair(begin, end)) else {
            return Err(ConfSeqBaseConformerError::Build(
                "ring atom order does not map to adjacent ring bonds".to_string(),
            ));
        };
        if bond.order() != BondOrder::Single || bond.is_aromatic() {
            return Ok(false);
        }
    }
    Ok(true)
}

fn is_confseq_base_puckerable_six_membered_ring(
    molecule: &Molecule,
    ring: &[AtomId],
) -> Result<bool, ConfSeqBaseConformerError> {
    if ring.len() != 6 {
        return Ok(false);
    }
    let has_sp3_center = ring
        .iter()
        .any(|atom| molecule.atoms()[atom.index()].hybridization() == Hybridization::Sp3);
    if !has_sp3_center {
        return Ok(false);
    }
    for atom in ring {
        let atom = &molecule.atoms()[atom.index()];
        if atom.is_aromatic() || !matches!(atom.atomic_number(), 6 | 7 | 8 | 15 | 16) {
            return Ok(false);
        }
    }
    Ok(true)
}

fn saturated_five_membered_envelope_points(
    molecule: &Molecule,
    ring: &[AtomId],
    side_lengths: &[f64],
) -> Result<Vec<[f64; 3]>, ConfSeqBaseConformerError> {
    debug_assert_eq!(ring.len(), 5);
    debug_assert_eq!(side_lengths.len(), 5);
    let mut candidate_envelopes: Vec<_> = ring
        .iter()
        .enumerate()
        .filter_map(|(pos, atom)| {
            matches!(
                molecule.atoms()[atom.index()].chiral_tag(),
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
            )
            .then_some(pos)
        })
        .collect();
    if candidate_envelopes.is_empty() {
        candidate_envelopes.extend(0..ring.len());
    }

    let mut best_points = None;
    let mut best_score = f64::INFINITY;
    for envelope_pos in candidate_envelopes {
        let Some((points, score)) =
            optimized_five_membered_envelope_for_atom(molecule, ring, side_lengths, envelope_pos)
        else {
            continue;
        };
        if score < best_score {
            best_score = score;
            best_points = Some(points);
        }
    }
    best_points.ok_or_else(|| {
        ConfSeqBaseConformerError::Build(
            "saturated five-membered ring cannot form an envelope base polygon".to_string(),
        )
    })
}

fn optimized_five_membered_envelope_for_atom(
    molecule: &Molecule,
    ring: &[AtomId],
    side_lengths: &[f64],
    envelope_pos: usize,
) -> Option<(Vec<[f64; 3]>, f64)> {
    let prev_edge = (envelope_pos + side_lengths.len() - 1) % side_lengths.len();
    let next_edge = envelope_pos;
    let upper = 0.5 * side_lengths[prev_edge].min(side_lengths[next_edge]) * (1.0 - 1.0e-10);
    if !upper.is_finite() || upper <= 0.0 {
        return None;
    }
    let mut best_height = None;
    let mut best_score = f64::INFINITY;
    for sample in 0..=96 {
        let height = upper * sample as f64 / 96.0;
        let score = five_membered_envelope_local_angle_score(
            molecule,
            ring,
            side_lengths,
            envelope_pos,
            height,
        );
        if score.is_finite() && score < best_score {
            best_score = score;
            best_height = Some(height);
        }
    }
    let mut center = best_height?;
    let mut step = upper / 96.0;
    for _ in 0..48 {
        let left = (center - step).max(0.0);
        let right = (center + step).min(upper);
        let left_mid = left + (right - left) / 3.0;
        let right_mid = right - (right - left) / 3.0;
        let left_score = five_membered_envelope_local_angle_score(
            molecule,
            ring,
            side_lengths,
            envelope_pos,
            left_mid,
        );
        let right_score = five_membered_envelope_local_angle_score(
            molecule,
            ring,
            side_lengths,
            envelope_pos,
            right_mid,
        );
        center = if left_score <= right_score {
            left_mid
        } else {
            right_mid
        };
        step *= 0.5;
    }
    let points = five_membered_envelope_points_for_height(side_lengths, envelope_pos, center)?;
    let score = five_membered_envelope_local_angle_score(
        molecule,
        ring,
        side_lengths,
        envelope_pos,
        center,
    );
    Some((points, score))
}

fn five_membered_envelope_local_angle_score(
    molecule: &Molecule,
    ring: &[AtomId],
    side_lengths: &[f64],
    envelope_pos: usize,
    height: f64,
) -> f64 {
    let Some(points) = five_membered_envelope_points_for_height(side_lengths, envelope_pos, height)
    else {
        return f64::INFINITY;
    };
    let mut score = 0.0;
    for idx in 0..5 {
        let prev = (idx + 4) % 5;
        let next = (idx + 1) % 5;
        let observed = angle_rad(points[prev], points[idx], points[next]);
        let target = confseq_base_ring_angle_rad(molecule, ring[idx].index(), ring.len());
        let delta = observed - target;
        score += delta * delta;
    }
    score
}

fn five_membered_envelope_points_for_height(
    side_lengths: &[f64],
    envelope_pos: usize,
    height: f64,
) -> Option<Vec<[f64; 3]>> {
    debug_assert_eq!(side_lengths.len(), 5);
    let mut xy_lengths = side_lengths.to_vec();
    let prev_edge = (envelope_pos + side_lengths.len() - 1) % side_lengths.len();
    let next_edge = envelope_pos;
    for edge in [prev_edge, next_edge] {
        let xy_sq = side_lengths[edge] * side_lengths[edge] - height * height;
        if xy_sq <= 0.0 {
            return None;
        }
        xy_lengths[edge] = xy_sq.sqrt();
    }
    let xy_points = planar_cyclic_polygon_points(&xy_lengths).ok()?;
    Some(
        xy_points
            .into_iter()
            .enumerate()
            .map(|(idx, point)| {
                [
                    point[0],
                    point[1],
                    if idx == envelope_pos { height } else { 0.0 },
                ]
            })
            .collect(),
    )
}

fn saturated_six_membered_chair_points(
    molecule: &Molecule,
    ring: &[AtomId],
    side_lengths: &[f64],
) -> Result<Vec<[f64; 3]>, ConfSeqBaseConformerError> {
    // Analytic chair subset for equal-edge saturated six-membered rings. For
    // edge length L, choosing in-plane radius a=L*sqrt(8/9) and alternating
    // z=+/-L/6 makes adjacent bond lengths L and ring bond angles tetrahedral.
    if side_lengths
        .iter()
        .all(|length| (*length - side_lengths[0]).abs() <= 1.0e-12)
        && ring
            .iter()
            .all(|atom| molecule.atoms()[atom.index()].atomic_number() == 6)
    {
        return Ok(equal_edge_saturated_six_membered_chair_points(
            side_lengths[0],
        ));
    }

    // Mixed saturated extension: keep every ring edge length exact and choose a
    // single chair-like puckering height by deterministic local-angle least
    // squares. This is a constrained base scaffold, not a force field or
    // conformer search; non-bonded terms and global energy are intentionally
    // absent.
    optimized_puckered_six_membered_ring_points(molecule, ring, side_lengths)
}

fn equal_edge_saturated_six_membered_chair_points(length: f64) -> Vec<[f64; 3]> {
    let radius = length * (8.0_f64 / 9.0).sqrt();
    let z = length / 6.0;
    (0..6)
        .map(|idx| {
            let theta = idx as f64 * PI / 3.0;
            [
                radius * theta.cos(),
                radius * theta.sin(),
                if idx % 2 == 0 { z } else { -z },
            ]
        })
        .collect()
}

fn optimized_puckered_six_membered_ring_points(
    molecule: &Molecule,
    ring: &[AtomId],
    side_lengths: &[f64],
) -> Result<Vec<[f64; 3]>, ConfSeqBaseConformerError> {
    debug_assert_eq!(ring.len(), 6);
    debug_assert_eq!(side_lengths.len(), 6);
    let min_side = side_lengths.iter().copied().fold(f64::INFINITY, f64::min);
    if !min_side.is_finite() || min_side <= 0.0 {
        return Err(ConfSeqBaseConformerError::Build(
            "saturated six-membered ring has invalid side lengths".to_string(),
        ));
    }
    let upper = 0.5 * min_side * (1.0 - 1.0e-10);
    let mut best_height = None;
    let mut best_score = f64::INFINITY;
    for sample in 0..=96 {
        let height = upper * sample as f64 / 96.0;
        let score =
            puckered_six_membered_ring_local_angle_score(molecule, ring, side_lengths, height);
        if score.is_finite() && score < best_score {
            best_score = score;
            best_height = Some(height);
        }
    }
    let Some(mut center) = best_height else {
        return Err(ConfSeqBaseConformerError::Build(
            "saturated six-membered ring cannot form a puckered base polygon".to_string(),
        ));
    };
    let mut step = upper / 96.0;
    for _ in 0..48 {
        let left = (center - step).max(0.0);
        let right = (center + step).min(upper);
        let left_mid = left + (right - left) / 3.0;
        let right_mid = right - (right - left) / 3.0;
        let left_score =
            puckered_six_membered_ring_local_angle_score(molecule, ring, side_lengths, left_mid);
        let right_score =
            puckered_six_membered_ring_local_angle_score(molecule, ring, side_lengths, right_mid);
        center = if left_score <= right_score {
            left_mid
        } else {
            right_mid
        };
        step *= 0.5;
    }
    puckered_six_membered_ring_points_for_height(side_lengths, center).ok_or_else(|| {
        ConfSeqBaseConformerError::Build(
            "saturated six-membered ring optimized height is invalid".to_string(),
        )
    })
}

fn puckered_six_membered_ring_local_angle_score(
    molecule: &Molecule,
    ring: &[AtomId],
    side_lengths: &[f64],
    height: f64,
) -> f64 {
    let Some(points) = puckered_six_membered_ring_points_for_height(side_lengths, height) else {
        return f64::INFINITY;
    };
    let mut score = 0.0;
    for idx in 0..6 {
        let prev = (idx + 5) % 6;
        let next = (idx + 1) % 6;
        let observed = angle_rad(points[prev], points[idx], points[next]);
        let target = confseq_base_ring_angle_rad(molecule, ring[idx].index(), ring.len());
        let delta = observed - target;
        score += delta * delta;
    }
    score
}

fn puckered_six_membered_ring_points_for_height(
    side_lengths: &[f64],
    height: f64,
) -> Option<Vec<[f64; 3]>> {
    let z_delta = 2.0 * height;
    let mut xy_lengths = Vec::with_capacity(side_lengths.len());
    for side in side_lengths {
        let xy_sq = side * side - z_delta * z_delta;
        if xy_sq <= 0.0 {
            return None;
        }
        xy_lengths.push(xy_sq.sqrt());
    }
    let xy_points = planar_cyclic_polygon_points(&xy_lengths).ok()?;
    Some(
        xy_points
            .into_iter()
            .enumerate()
            .map(|(idx, point)| {
                [
                    point[0],
                    point[1],
                    if idx % 2 == 0 { height } else { -height },
                ]
            })
            .collect(),
    )
}

fn planar_cyclic_polygon_points(
    side_lengths: &[f64],
) -> Result<Vec<[f64; 2]>, ConfSeqBaseConformerError> {
    let max_side = side_lengths.iter().copied().fold(0.0, f64::max);
    let total: f64 = side_lengths.iter().sum();
    if side_lengths.len() < 3 || max_side <= 0.0 || max_side * 2.0 >= total {
        return Err(ConfSeqBaseConformerError::Build(
            "ring side lengths cannot form a polygon".to_string(),
        ));
    }

    let mut low = max_side * 0.5 * (1.0 + 1.0e-12);
    let mut high = low * 2.0;
    while cyclic_polygon_angle_sum(side_lengths, high) > 2.0 * PI {
        high *= 2.0;
    }
    if cyclic_polygon_angle_sum(side_lengths, low) < 2.0 * PI {
        return Err(ConfSeqBaseConformerError::Build(
            "ring side lengths cannot form a cyclic polygon".to_string(),
        ));
    }
    for _ in 0..80 {
        let mid = 0.5 * (low + high);
        if cyclic_polygon_angle_sum(side_lengths, mid) > 2.0 * PI {
            low = mid;
        } else {
            high = mid;
        }
    }

    let radius = 0.5 * (low + high);
    let mut theta = 0.0_f64;
    let mut points = Vec::with_capacity(side_lengths.len());
    for side in side_lengths {
        points.push([radius * theta.cos(), radius * theta.sin()]);
        theta += 2.0 * (side / (2.0 * radius)).clamp(-1.0, 1.0).asin();
    }
    Ok(points)
}

fn cyclic_polygon_angle_sum(side_lengths: &[f64], radius: f64) -> f64 {
    side_lengths
        .iter()
        .map(|side| 2.0 * (side / (2.0 * radius)).clamp(-1.0, 1.0).asin())
        .sum()
}

fn centroid_3d(points: &[[f64; 3]]) -> [f64; 3] {
    let mut sum = [0.0, 0.0, 0.0];
    for point in points {
        sum = vec_add(sum, *point);
    }
    vec_scale(sum, 1.0 / points.len() as f64)
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

fn place_confseq_tree(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    atom_idx: usize,
    parent: Option<usize>,
    point: [f64; 3],
    coords: &mut [[f64; 3]],
    placed: &mut [bool],
    geometry: &ConfSeqBaseGeometry,
) {
    coords[atom_idx] = point;
    placed[atom_idx] = true;

    let base_axis = parent
        .map(|parent| vec_normalize(vec_sub(coords[atom_idx], coords[parent])))
        .unwrap_or([1.0, 0.0, 0.0]);
    let children: Vec<_> = adjacency
        .neighbors_of(atom_idx)
        .iter()
        .filter(|nbr| Some(nbr.atom_index) != parent && !placed[nbr.atom_index])
        .collect();
    let child_count = children.len();
    for (child_ord, nbr) in children.into_iter().enumerate() {
        let bond = &molecule.bonds()[nbr.bond.index()];
        let length = geometry.bond_length(bond);
        let angle = confseq_base_local_angle_rad(molecule, atom_idx);
        let dir = child_direction(base_axis, child_ord, child_count, angle);
        let child_point = vec_add(coords[atom_idx], vec_scale(dir, length));
        place_confseq_tree(
            molecule,
            adjacency,
            nbr.atom_index,
            Some(atom_idx),
            child_point,
            coords,
            placed,
            geometry,
        );
    }
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

struct ConfSeqBaseGeometry {
    bond_lengths: Vec<f64>,
}

impl ConfSeqBaseGeometry {
    fn new(molecule: &Molecule, ring_info: &rings::RingInfo) -> Self {
        let ring_bonds: HashSet<_> = ring_info
            .bond_rings()
            .iter()
            .flat_map(|ring| ring.iter().map(|bond| bond.index()))
            .collect();
        let bond_lengths = molecule
            .bonds()
            .iter()
            .map(|bond| {
                // Non-aromatic ring blocks are constrained base scaffolds; keep
                // their static edge lengths to avoid hidden ring-shape drift.
                // Aromatic ring and non-ring bonds use the source-backed UFF
                // rest length already ported in core.
                if ring_bonds.contains(&bond.id().index()) && !bond.is_aromatic() {
                    confseq_base_static_bond_length(molecule, bond)
                } else {
                    confseq_base_source_backed_bond_length(molecule, bond)
                }
            })
            .collect();
        Self { bond_lengths }
    }

    fn bond_length(&self, bond: &Bond) -> f64 {
        self.bond_lengths
            .get(bond.id().index())
            .copied()
            .unwrap_or_else(|| confseq_base_static_bond_length_fallback(bond))
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

fn confseq_base_bond_length(molecule: &Molecule, bond: &Bond) -> f64 {
    confseq_base_source_backed_bond_length(molecule, bond)
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

fn angle_rad(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ba = vec_sub(a, b);
    let bc = vec_sub(c, b);
    let denom = vec_len(ba) * vec_len(bc);
    if denom <= 1.0e-12 {
        return 0.0;
    }
    (vec_dot(ba, bc) / denom).clamp(-1.0, 1.0).acos()
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
    fn confseq_base_conformer_corrects_ring_bound_tetrahedral_stereo_with_cyclic_side_fallback() {
        let smiles =
            "Cc1ccc(c(c1F)C(=O)N2C[C@H](C3(C2)CC[NH2+]CC3)C(=O)N[C@@H]4CCC[C@@H]([C@@H]4C)C)F";
        let molecule = Molecule::from_smiles(smiles).expect("test SMILES parses");
        let embedded =
            try_build_confseq_base_conformer(&molecule).expect("base conformer should build");
        let coords = conformer_points(&embedded);
        let adjacency = AdjacencyList::from_topology(embedded.num_atoms(), embedded.bonds());
        let constraints =
            collect_confseq_base_tetrahedral_stereo_constraints(&embedded, &adjacency)
                .expect("tetrahedral constraints should collect");

        assert!(
            constraints.iter().any(
                |constraint| constraint.ligands.iter().copied().any(|ligand| {
                    ligand != constraint.center
                        && confseq_base_chiral_movable_side(
                            &embedded, &adjacency, constraint, ligand,
                        )
                        .is_ok_and(|side| {
                            confseq_base_chiral_side_contains_other_ligand(
                                &side, constraint, ligand,
                            )
                        })
                })
            ),
            "fixture should exercise the cyclic-side stereo correction candidate"
        );
        for constraint in &constraints {
            let volume = confseq_base_chiral_volume(&coords, constraint);
            assert!(
                confseq_base_chiral_volume_satisfies_tag(volume, constraint.tag),
                "tetrahedral stereo at atom {} should satisfy {:?}, volume={volume}",
                constraint.center,
                constraint.tag
            );
        }
    }

    #[test]
    fn confseq_base_conformer_separates_two_substituents_on_ring_bound_sp3_center() {
        let smiles = "C[C@H](c1[n-]nnn1)[NH2+][C@@H]2CS(=O)(=O)c3c2cccc3";
        let molecule = Molecule::from_smiles(smiles).expect("test SMILES parses");
        let embedded =
            try_build_confseq_base_conformer(&molecule).expect("base conformer should build");
        let coords = conformer_points(&embedded);

        let sulfur = embedded
            .atoms()
            .iter()
            .find(|atom| {
                atom.atomic_number() == 16 && heavy_degree(&embedded, atom.id().index()) == 4
            })
            .expect("fixture should contain a tetravalent ring-bound sulfur")
            .id()
            .index();
        let oxygen_neighbors: Vec<_> = embedded
            .bonds()
            .iter()
            .filter_map(|bond| {
                let begin = bond.begin().index();
                let end = bond.end().index();
                let other = if begin == sulfur {
                    end
                } else if end == sulfur {
                    begin
                } else {
                    return None;
                };
                (embedded.atoms()[other].atomic_number() == 8).then_some(other)
            })
            .collect();

        assert_eq!(oxygen_neighbors.len(), 2);
        let angle = angle_deg(
            coords[oxygen_neighbors[0]],
            coords[sulfur],
            coords[oxygen_neighbors[1]],
        );
        assert!(
            angle > 30.0,
            "two non-ring substituents on ring-bound sp3 sulfur must not collapse, angle={angle:.6}"
        );
    }

    #[test]
    fn confseq_base_conformer_uses_analytic_chair_for_saturated_carbon_six_rings() {
        let template = base_template_for_comparison("C1CCCCC1");
        let coords = conformer_points(&template.molecule);

        for idx in 0..6 {
            let next = (idx + 1) % 6;
            let bond = bond_between_pair(&template.molecule, sorted_pair(idx, next))
                .expect("cyclohexane edge should exist");
            let expected = confseq_base_static_bond_length(&template.molecule, bond);
            assert!(
                (distance(coords[idx], coords[next]) - expected).abs() < 1.0e-10,
                "cyclohexane chair edge {idx}-{next} should use the base ring C-C single length"
            );
        }
        for idx in 0..6 {
            let prev = (idx + 5) % 6;
            let next = (idx + 1) % 6;
            let angle = angle_deg(coords[prev], coords[idx], coords[next]);
            assert!(
                (angle - 109.47122063449069).abs() < 1.0e-10,
                "cyclohexane chair angle at {idx} should be tetrahedral, got {angle}"
            );
        }
    }

    #[test]
    fn confseq_base_conformer_uses_puckered_six_rings_for_saturated_heterocycles() {
        for smiles in ["C1CCNCC1", "C1COCCN1", "C1CCSCC1"] {
            let template = base_template_for_comparison(smiles);
            let coords = conformer_points(&template.molecule);
            let rings = rings::symmetrize_sssr(&template.molecule)
                .expect("heterocycle ring perception should work");
            let ring = &rings.atom_rings()[0];

            for pos in 0..ring.len() {
                let begin = ring[pos].index();
                let end = ring[(pos + 1) % ring.len()].index();
                let bond = bond_between_pair(&template.molecule, sorted_pair(begin, end))
                    .expect("ring edge should exist");
                let expected = confseq_base_static_bond_length(&template.molecule, bond);
                assert!(
                    (distance(coords[begin], coords[end]) - expected).abs() < 1.0e-8,
                    "saturated heterocycle {smiles} should preserve ring edge {begin}-{end}"
                );
            }

            let max_abs_z = ring
                .iter()
                .map(|atom| coords[atom.index()][2].abs())
                .fold(0.0, f64::max);
            assert!(
                max_abs_z > 0.05,
                "saturated heterocycle {smiles} should not collapse to a planar ring"
            );
        }
    }

    #[test]
    fn confseq_base_conformer_preserves_saturated_five_ring_edges_without_forcefield() {
        let template = base_template_for_comparison("C1CCCC1");
        let coords = conformer_points(&template.molecule);

        let max_abs_z = coords
            .iter()
            .map(|point| point[2].abs())
            .fold(0.0, f64::max);
        assert!(max_abs_z.is_finite());

        for idx in 0..5 {
            let next = (idx + 1) % 5;
            let bond = bond_between_pair(&template.molecule, sorted_pair(idx, next))
                .expect("cyclopentane edge should exist");
            let expected = confseq_base_static_bond_length(&template.molecule, bond);
            assert!(
                (distance(coords[idx], coords[next]) - expected).abs() < 1.0e-10,
                "cyclopentane base edge {idx}-{next} should use the base ring C-C single length"
            );
        }
    }

    #[test]
    fn confseq_base_ring_substituents_satisfy_local_center_angles() {
        for smiles in ["Cc1ccccc1", "CC1CCNCC1", "C1CCN(CC1)C"] {
            let template = base_template_for_comparison(smiles);
            let molecule = &template.molecule;
            let coords = conformer_points(molecule);
            let adjacency = AdjacencyList::from_topology(molecule.num_atoms(), molecule.bonds());
            let rings = rings::symmetrize_sssr(molecule).expect("ring perception should work");

            let mut checked = 0;
            for ring in rings.atom_rings() {
                let ring_atoms: HashSet<_> = ring.iter().map(|atom| atom.index()).collect();
                for atom in ring {
                    let center = atom.index();
                    let Some(outside) = adjacency
                        .neighbors_of(center)
                        .iter()
                        .map(|nbr| nbr.atom_index)
                        .find(|idx| !ring_atoms.contains(idx))
                    else {
                        continue;
                    };
                    let ring_neighbors: Vec<_> = adjacency
                        .neighbors_of(center)
                        .iter()
                        .map(|nbr| nbr.atom_index)
                        .filter(|idx| ring_atoms.contains(idx))
                        .collect();
                    assert_eq!(ring_neighbors.len(), 2);
                    let target = confseq_base_local_angle_rad(molecule, center).to_degrees();
                    for ring_neighbor in ring_neighbors {
                        let observed =
                            angle_deg(coords[ring_neighbor], coords[center], coords[outside]);
                        assert!(
                            (observed - target).abs() < 1.0e-8,
                            "smiles={smiles} center={center} ring_neighbor={ring_neighbor} outside={outside} expected {target}, got {observed}"
                        );
                    }
                    checked += 1;
                }
            }
            assert!(
                checked > 0,
                "test molecule {smiles} should have a ring substituent"
            );
        }
    }

    #[test]
    fn confseq_base_conformer_builds_acyclic_small_molecules_with_local_geometry() {
        for smiles in [
            "CCCC",
            "CCO",
            "CC(C)C",
            "CC(C)(C)C",
            "CC(C)CO",
            "CC=O",
            "CC#N",
            "C=CC",
            "C#CC",
            "CC(=O)N",
            "COC(=O)C",
            "CCF",
            "CCCl",
            "CCBr",
            "CCI",
            "CCSC",
            "CS(=O)C",
            "CS(=O)(=O)C",
            "CP(=O)(O)O",
            "C(F)(Cl)Br",
            "CC(=O)OC",
        ] {
            let uff_template = optimized_template_for_comparison(smiles);
            let base_template = base_template_for_comparison(smiles);
            let bond_rms =
                bond_length_rms_against_reference(&uff_template.molecule, &base_template.molecule);
            let angle_rms =
                angle_rms_against_reference(&uff_template.molecule, &base_template.molecule);
            let all_heavy_angle_rms = all_heavy_angle_rms_against_reference(
                &uff_template.molecule,
                &base_template.molecule,
            );

            eprintln!(
                "confseq_base acyclic smiles={smiles} bond_rms={bond_rms:.6} angle_rms={angle_rms:.6} all_heavy_angle_rms={all_heavy_angle_rms:.6}"
            );
            assert!(
                bond_rms < 0.08,
                "acyclic base bond lengths should stay near UFF-relaxed path"
            );
            assert!(
                angle_rms < 8.0,
                "acyclic base local angles should stay near UFF-relaxed path"
            );
            assert!(
                all_heavy_angle_rms < 8.0,
                "acyclic base all-heavy local angles should stay near UFF-relaxed path"
            );
        }
    }

    #[test]
    fn confseq_base_conformer_builds_single_aromatic_ring_as_rigid_planar_polygon() {
        for smiles in [
            "c1ccccc1",
            "Cc1ccccc1",
            "CCc1ccccc1",
            "Cc1ccc(C)cc1",
            "c1ccncc1",
            "Cc1ccncc1",
            "c1cnccn1",
            "c1ccoc1",
            "Cc1ccoc1",
            "c1ccsc1",
            "c1cc[nH]c1",
            "c1ncc[nH]1",
        ] {
            let uff_template = optimized_template_for_comparison(smiles);
            let base_template = base_template_for_comparison(smiles);
            let bond_rms =
                bond_length_rms_against_reference(&uff_template.molecule, &base_template.molecule);
            let all_heavy_angle_rms = all_heavy_angle_rms_against_reference(
                &uff_template.molecule,
                &base_template.molecule,
            );
            let aromatic_atoms = aromatic_atom_indices(&base_template.molecule);
            let coords = conformer_points(&base_template.molecule);
            let max_abs_z = aromatic_atoms
                .iter()
                .map(|idx| coords[*idx][2].abs())
                .fold(0.0, f64::max);

            eprintln!(
                "confseq_base aromatic_ring smiles={smiles} bond_rms={bond_rms:.6} all_heavy_angle_rms={all_heavy_angle_rms:.6} max_ring_abs_z={max_abs_z:.6}"
            );
            assert!(bond_rms < 0.10);
            assert!(all_heavy_angle_rms < 12.0);
            assert!(max_abs_z < 1.0e-9);
        }
    }

    #[test]
    fn confseq_base_conformer_builds_disjoint_aromatic_ring_blocks() {
        for smiles in [
            "c1ccccc1-c2ccccc2",
            "c1ccccc1Oc2ccccc2",
            "c1ccccc1CCc2ccccc2",
            "Cc1ccc(Oc2ccncc2)cc1",
            "c1ccoc1-c2ccsc2",
        ] {
            let uff_template = optimized_template_for_comparison(smiles);
            let base_template = base_template_for_comparison(smiles);
            let bond_rms =
                bond_length_rms_against_reference(&uff_template.molecule, &base_template.molecule);
            let all_heavy_angle_rms = all_heavy_angle_rms_against_reference(
                &uff_template.molecule,
                &base_template.molecule,
            );
            let coords = conformer_points(&base_template.molecule);
            let rings = rings::fast_find_rings(&base_template.molecule)
                .expect("base template ring perception should work");

            for ring in rings.atom_rings() {
                let max_abs_z = ring
                    .iter()
                    .map(|atom| coords[atom.index()][2].abs())
                    .fold(0.0, f64::max);
                assert!(
                    max_abs_z < 1.0e-9,
                    "each aromatic ring block stays planar in the base frame"
                );
            }

            eprintln!(
                "confseq_base disjoint_aromatic smiles={smiles} bond_rms={bond_rms:.6} all_heavy_angle_rms={all_heavy_angle_rms:.6}"
            );
            assert!(bond_rms < 0.12);
            assert!(all_heavy_angle_rms < 18.0);
        }
    }

    #[test]
    fn confseq_base_conformer_builds_edge_fused_aromatic_ring_components() {
        for smiles in [
            "c1ccc2ccccc2c1",
            "c1ccc2ncccc2c1",
            "c1ccc2occc2c1",
            "c1ccc2cc3ccccc3cc2c1",
            "c1ccc2c(c1)ccc3ccccc23",
            "c1ccc2c(c1)[nH]c3ccccc23",
        ] {
            let uff_template = optimized_template_for_comparison(smiles);
            let base_template = base_template_for_comparison(smiles);
            let bond_rms =
                bond_length_rms_against_reference(&uff_template.molecule, &base_template.molecule);
            let all_heavy_angle_rms = all_heavy_angle_rms_against_reference(
                &uff_template.molecule,
                &base_template.molecule,
            );
            let aromatic_atoms = aromatic_atom_indices(&base_template.molecule);
            let coords = conformer_points(&base_template.molecule);
            let max_abs_z = aromatic_atoms
                .iter()
                .map(|idx| coords[*idx][2].abs())
                .fold(0.0, f64::max);

            eprintln!(
                "confseq_base fused_aromatic smiles={smiles} bond_rms={bond_rms:.6} all_heavy_angle_rms={all_heavy_angle_rms:.6} max_abs_z={max_abs_z:.6}"
            );
            assert!(bond_rms < 0.12);
            assert!(all_heavy_angle_rms < 18.0);
            assert!(max_abs_z < 1.0e-9);
        }
    }

    #[test]
    fn confseq_base_conformer_builds_linked_and_substituted_fused_aromatic_components() {
        for smiles in [
            "Cc1ccc2ccccc2c1",
            "CCc1ccc2ccccc2c1",
            "COc1ccc2ccccc2c1",
            "c1ccc2ccccc2c1-c3ccccc3",
            "c1ccc2ccccc2c1CCc3ccccc3",
            "c1ccc2ccccc2c1-c3ccc4ccccc4c3",
        ] {
            let uff_template = optimized_template_for_comparison(smiles);
            let base_template = base_template_for_comparison(smiles);
            let bond_rms =
                bond_length_rms_against_reference(&uff_template.molecule, &base_template.molecule);
            let all_heavy_angle_rms = all_heavy_angle_rms_against_reference(
                &uff_template.molecule,
                &base_template.molecule,
            );
            let aromatic_atoms = aromatic_atom_indices(&base_template.molecule);
            let coords = conformer_points(&base_template.molecule);
            let max_abs_z = aromatic_atoms
                .iter()
                .map(|idx| coords[*idx][2].abs())
                .fold(0.0, f64::max);

            eprintln!(
                "confseq_base linked_fused smiles={smiles} bond_rms={bond_rms:.6} all_heavy_angle_rms={all_heavy_angle_rms:.6} max_abs_z={max_abs_z:.6}"
            );
            assert!(bond_rms < 0.12);
            assert!(all_heavy_angle_rms < 18.0);
            assert!(max_abs_z < 1.0e-9);
        }
    }

    #[test]
    fn confseq_base_conformer_builds_branched_aromatic_linker_components() {
        for smiles in [
            "c1ccccc1Oc2ccccc2",
            "c1ccccc1Sc2ccccc2",
            "c1ccccc1S(=O)(=O)c2ccccc2",
            "c1ccc2ccccc2c1S(=O)(=O)c3ccccc3",
            "c1ccccc1N(c2ccccc2)c3ccccc3",
            "c1ccc2ccccc2c1N(c3ccccc3)c4ccccc4",
        ] {
            let uff_template = optimized_template_for_comparison(smiles);
            let base_template = base_template_for_comparison(smiles);
            let bond_rms =
                bond_length_rms_against_reference(&uff_template.molecule, &base_template.molecule);
            let all_heavy_angle_rms = all_heavy_angle_rms_against_reference(
                &uff_template.molecule,
                &base_template.molecule,
            );
            let max_ring_plane_deviation =
                max_aromatic_ring_plane_deviation(&base_template.molecule);

            eprintln!(
                "confseq_base branched_aromatic smiles={smiles} bond_rms={bond_rms:.6} all_heavy_angle_rms={all_heavy_angle_rms:.6} max_ring_plane_deviation={max_ring_plane_deviation:.6}"
            );
            assert!(bond_rms < 0.12);
            assert!(all_heavy_angle_rms < 18.0);
            assert!(max_ring_plane_deviation < 1.0e-9);
        }
    }

    #[test]
    fn confseq_base_conformer_rejects_ring_systems_outside_current_subset() {
        let molecule = Molecule::from_smiles("C1CC2CCC1C2").expect("bridged system parses");
        let error = try_build_confseq_base_conformer(&molecule)
            .expect_err("bridged ring systems are not in the current ConfSeq-base subset");

        assert!(matches!(
            error,
            ConfSeqBaseConformerError::UnsupportedRingFusion { .. }
                | ConfSeqBaseConformerError::UnsupportedClosedRingFusion
        ));
    }

    #[test]
    fn confseq_base_conformer_builds_spiro_ring_forest() {
        for smiles in ["C1CCC2(CC1)CCC2", "C1CC2(CC1)CCNCC2"] {
            let base_template = base_template_for_comparison(smiles);
            let coords = conformer_points(&base_template.molecule);
            let bond_rms = bond_length_rms_against_reference(
                &optimized_template_for_comparison(smiles).molecule,
                &base_template.molecule,
            );
            let max_origin_distance = coords
                .iter()
                .map(|point| vec_len(*point))
                .fold(0.0, f64::max);

            eprintln!(
                "confseq_base spiro smiles={smiles} bond_rms={bond_rms:.6} max_origin_distance={max_origin_distance:.6}"
            );
            assert!(bond_rms < 0.20);
            assert!(max_origin_distance.is_finite());
        }
    }

    #[test]
    fn confseq_base_fusion_graph_guard_rejects_closed_fusion_components() {
        assert!(confseq_base_fusion_graph_is_forest(&[
            vec![1],
            vec![0, 2],
            vec![1]
        ]));
        assert!(!confseq_base_fusion_graph_is_forest(&[
            vec![1, 2],
            vec![0, 2],
            vec![0, 1]
        ]));
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
                .expect_err("unsupported base-conformer ring must not fallback to DistGeom");

        assert!(matches!(
            error,
            ConfSeqDecodeError::BaseConformer(
                ConfSeqBaseConformerError::UnsupportedRingFusion { .. }
                    | ConfSeqBaseConformerError::UnsupportedClosedRingFusion
            )
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
            Some(diagnostics::ConfSeqDiagnosticPhase::BaseTemplate)
        );
        assert!(diagnostic.parsed);
        assert!(diagnostic.distance_geometry_template_built);
        assert!(!diagnostic.base_template_built);
        assert!(matches!(
            diagnostic.base_error,
            Some(
                ConfSeqBaseConformerError::UnsupportedRingFusion { .. }
                    | ConfSeqBaseConformerError::UnsupportedClosedRingFusion
            )
        ));
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
}
