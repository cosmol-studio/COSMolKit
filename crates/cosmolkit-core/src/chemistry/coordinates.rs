//! 2D coordinate generation ported from the RDKit depiction algorithm.
//!
//! Source reproduction protocol: dev/source_reproduction_protocol.md
//!
//! Marker convention:
//!   RDKit✔️✔️: one-to-one port
//!   RDKit✔️❌: adapted for new-core API, same algorithm
//!   RDKit❗✔️: algorithm-equivalent via different architecture
//!
//! Source layout (`molblock.rs`, line numbers):
//!
//!   Types + helpers (~33 – 55)   norm, normalize, rotate, cos/sin/sqrt
//!   TreeEmbeddedAtom             ~825
//!   RdkitHybridization           ~758
//!   atom_depict_rank             ~838
//!   rdkit_ring_radius            ~847
//!   n_outer_electrons_rdkit      ~852
//!   rdkit_num_bonds_plus_lone_pairs     ~856
//!   rdkit_hybridizations_for_depict     ~893
//!   compute_sub_angle            ~1764
//!   rotation_dir                 ~1899
//!   canonicalize_component       ~769
//!   rdkit_rank_atoms_by_rank     ~1705
//!   rdkit_set_nbr_order          ~1723
//!   rdkit_compute_initial_coords_strict ~503
//!   compute_2d_coords            ~446
//!
//! RDKit marker convention is defined in `dev/source_reproduction_protocol.md`.
//! Local shorthand:
//!   RDKit❗✔️: ported from compute_2d_coords_impl (intentionally different architecture)
//!   RDKit✔️✔️: one-to-one port
//!   RDKit✔️❌: adapted for new-core API, same algorithm
//!   RDKit❗❌: not yet ported (unimplemented stub)

use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};
use std::f64::consts::PI;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{OnceLock, RwLock};

unsafe extern "C" {
    #[link_name = "cos"]
    fn c_libm_cos(x: f64) -> f64;
    #[link_name = "sin"]
    fn c_libm_sin(x: f64) -> f64;
    #[link_name = "acos"]
    fn c_libm_acos(x: f64) -> f64;
    #[link_name = "sqrt"]
    fn c_libm_sqrt(x: f64) -> f64;
}

use crate::AtomId;
use crate::ChiralTag;
use crate::Molecule;
use crate::atom::Atom;
use crate::bond::{Bond, BondDirection, BondOrder, BondStereo};
use crate::builder::MoleculeBuilder;
use crate::search::smarts_parse::parse_smarts;
use crate::search::substruct::get_substruct_match;
use crate::valence::periodic_table_outer_electrons;
use crate::{AtomQueryPredicate, QueryNode};
use crate::{AtomSpec, BondSpec, CoordinateDimension, Element};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub(crate) enum Coordinate2DError {
    #[error("{0}")]
    InvalidInput(&'static str),
    #[error("{0}")]
    UnsupportedFeature(String),
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RdkitTemplateLineParts {
    smarts_body: String,
    cx_block: Option<String>,
    trailing_name: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct TemplateGraphBond {
    begin_atom_idx: usize,
    end_atom_idx: usize,
    query: crate::QueryNode<crate::BondQueryPredicate>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RdkitTemplateGraphModel {
    atom_queries: Vec<crate::QueryNode<crate::AtomQueryPredicate>>,
    bonds: Vec<TemplateGraphBond>,
}

#[derive(Debug, Clone, PartialEq)]
struct RdkitTemplateRuntimeModel {
    graph: RdkitTemplateGraphModel,
    coords_2d: Option<Vec<[f64; 2]>>,
    conformer_3d: Option<Vec<[f64; 3]>>,
    source_coordinate_dim: Option<CoordinateDimension>,
    fragment_count: usize,
    bond_ring_counts: Vec<usize>,
}

type RdkitTemplateBuckets = HashMap<usize, Vec<RdkitTemplateRuntimeModel>>;

const RDKIT_TEMPLATE_SMARTS_SOURCE: &str = include_str!("../data/rdkit_depictor_template_smarts.h");

#[derive(Debug, Default)]
struct RdkitCoordinateTemplateRegistry {
    templates: RdkitTemplateBuckets,
}

type RdkitIntPoint2DMap = BTreeMap<usize, (f64, f64)>;
const RDKIT_FORMER_NBR_INDICES_PROP: &str = "__formerNbrIndices";
const RDKIT_FORMER_IDX_PROP: &str = "__formerIdx";
const RDKIT_COLLISION_THRES: f64 = 0.70;
const RDKIT_HETEROATOM_COLL_SCALE: f64 = 1.3;
const RDKIT_BOND_THRES: f64 = 0.50;
const RDKIT_MAX_COLL_ITERS: usize = 15;
const RDKIT_NUM_BONDS_FLIPS: usize = 3;
const RDKIT_ANGLE_OPEN: f64 = 0.1222;
const RDKIT_RANDOM_MODULUS: u64 = 2_147_483_647;
const RDKIT_RANDOM_MULTIPLIER: u64 = 48_271;

#[derive(Debug, Clone)]
struct RdkitMinStdRand {
    state: u64,
}

impl Default for RdkitMinStdRand {
    fn default() -> Self {
        Self { state: 42 }
    }
}

impl RdkitMinStdRand {
    fn new(seed: i32) -> Self {
        if seed > 0 {
            Self { state: seed as u64 }
        } else {
            Self::default()
        }
    }

    fn next_raw(&mut self) -> u32 {
        self.state = (self.state * RDKIT_RANDOM_MULTIPLIER) % RDKIT_RANDOM_MODULUS;
        self.state as u32
    }

    fn next_bounded(&mut self, upper_inclusive: usize) -> usize {
        if upper_inclusive == 0 {
            return 0;
        }
        (self.next_raw() as usize) % (upper_inclusive + 1)
    }

    fn next_unit_f64(&mut self) -> f64 {
        self.next_raw() as f64 / RDKIT_RANDOM_MODULUS as f64
    }
}

fn parse_rdkit_template_line(line: &str) -> Result<RdkitTemplateLineParts, Coordinate2DError> {
    // RDKit source: SmilesParse.cpp / CXSmilesOps.cpp template-loading call chain
    // RDKit✔️✔️:   std::string lsmarts, name, cxPart;
    // RDKit✔️✔️:   preprocessSmiles(smarts, params, lsmarts, name, cxPart);
    // RDKit✔️✔️:   auto res = toMol(labelRecursivePatterns(lsmarts), smarts_parse, lsmarts);
    // RDKit✔️✔️:   handleCXPartAndName(res.get(), params, cxPart, name);
    // RDKit✔️✔️:   if (!res || cxPart.empty()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (params.allowCXSMILES) {
    // RDKit✔️✔️:     if (*pos == '|') {
    // RDKit✔️✔️:       SmilesParseOps::parseCXExtensions(*res, cxPart, pos);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    //
    // RDKit✔️❌: We port only the line-splitting contract needed before the
    // later template graph/coordinate model exists. Unsupported suffix forms
    // fail closed instead of being guessed into partial template state.
    let trimmed = line.trim();
    if trimmed.is_empty() {
        return Err(Coordinate2DError::UnsupportedFeature(
            "RDKit template line parsing requires a non-empty SMARTS line".to_string(),
        ));
    }

    let Some(split_idx) = trimmed.find([' ', '\t']) else {
        return Ok(RdkitTemplateLineParts {
            smarts_body: trimmed.to_string(),
            cx_block: None,
            trailing_name: None,
        });
    };

    let smarts_body = trimmed[..split_idx].to_string();
    let suffix = trimmed[split_idx..].trim();
    if smarts_body.is_empty() {
        return Err(Coordinate2DError::UnsupportedFeature(
            "RDKit template line parsing produced an empty SMARTS body".to_string(),
        ));
    }
    if suffix.is_empty() {
        return Ok(RdkitTemplateLineParts {
            smarts_body,
            cx_block: None,
            trailing_name: None,
        });
    }
    if !suffix.starts_with('|') {
        return Err(Coordinate2DError::UnsupportedFeature(
            "RDKit template line parsing does not support template suffix text without a CX block"
                .to_string(),
        ));
    }

    let Some(cx_end) = suffix[1..].find('|').map(|idx| idx + 1) else {
        return Err(Coordinate2DError::UnsupportedFeature(
            "RDKit template line parsing requires a closing '|' for the CX block".to_string(),
        ));
    };
    let cx_block = suffix[..=cx_end].to_string();
    let trailing_name = suffix[cx_end + 1..].trim();

    Ok(RdkitTemplateLineParts {
        smarts_body,
        cx_block: Some(cx_block),
        trailing_name: (!trailing_name.is_empty()).then(|| trailing_name.to_string()),
    })
}

fn parse_rdkit_template_graph_model(
    smarts: &str,
) -> Result<RdkitTemplateGraphModel, Coordinate2DError> {
    let parsed = parse_smarts(smarts).map_err(|error| {
        Coordinate2DError::UnsupportedFeature(format!(
            "RDKit template SMARTS graph parsing failed: {error}"
        ))
    })?;
    let bonds = expand_template_smarts_bonds(smarts)?;
    let atom_count_from_topology = bonds
        .iter()
        .flat_map(|bond| [bond.begin_atom_idx, bond.end_atom_idx])
        .max()
        .map_or(parsed.num_atoms(), |max_idx| max_idx + 1);
    if atom_count_from_topology != parsed.num_atoms() {
        return Err(Coordinate2DError::UnsupportedFeature(format!(
            "RDKit template SMARTS graph atom count mismatch: query parser={}, topology parser={atom_count_from_topology}",
            parsed.num_atoms(),
        )));
    }
    Ok(RdkitTemplateGraphModel {
        atom_queries: parsed.atom_queries,
        bonds,
    })
}

fn build_rdkit_template_runtime_model(
    line: &str,
) -> Result<RdkitTemplateRuntimeModel, Coordinate2DError> {
    let parts = parse_rdkit_template_line(line)?;
    let graph = parse_rdkit_template_graph_model(&parts.smarts_body)?;
    let (coords_2d, conformer_3d, source_coordinate_dim) = match parts.cx_block.as_deref() {
        Some(cx_block) => parse_template_cx_coordinate_block(cx_block, graph.atom_queries.len())?,
        None => (None, None, None),
    };
    let topology_molecule = build_template_topology_probe_molecule(&graph)?;
    let fragment_count = crate::fragment::get_num_fragments(&topology_molecule);
    let rings = crate::symmetrize_sssr(&topology_molecule).map_err(|error| {
        Coordinate2DError::UnsupportedFeature(format!(
            "RDKit template runtime ring initialization failed: {error}"
        ))
    })?;
    let bond_ring_counts = topology_molecule
        .bonds()
        .iter()
        .map(|bond| rings.num_bond_rings(bond.id()))
        .collect();
    Ok(RdkitTemplateRuntimeModel {
        graph,
        coords_2d,
        conformer_3d,
        source_coordinate_dim,
        fragment_count,
        bond_ring_counts,
    })
}

fn assert_valid_rdkit_template(
    runtime: &RdkitTemplateRuntimeModel,
    smiles: &str,
) -> Result<(), Coordinate2DError> {
    // RDKit✔️✔️:   // template must have 2D coordinates
    // RDKit✔️✔️:   if (mol.getNumConformers() == 0) {
    // RDKit✔️✔️:     std::string msg = "Template missing coordinates: " + smiles;
    // RDKit✔️✔️:     throw RDDepict::DepictException(msg);
    // RDKit✔️✔️:   }
    if runtime.coords_2d.is_none() && runtime.conformer_3d.is_none() {
        return Err(Coordinate2DError::UnsupportedFeature(format!(
            "Template missing coordinates: {smiles}"
        )));
    }
    // RDKit✔️✔️:   if (mol.getConformer().is3D()) {
    // RDKit✔️✔️:     std::string msg =
    // RDKit✔️✔️:         "Template has 3D coordinates, 2D coordinates required: " + smiles;
    // RDKit✔️✔️:     throw RDDepict::DepictException(msg);
    // RDKit✔️✔️:   }
    if runtime.source_coordinate_dim == Some(CoordinateDimension::ThreeD) {
        return Err(Coordinate2DError::UnsupportedFeature(format!(
            "Template has 3D coordinates, 2D coordinates required: {smiles}"
        )));
    }

    // RDKit✔️✔️:   // Make sure this template is a single ring system (spiro'd ring systems are
    // RDKit✔️✔️:   // OK). We can check this by ensuring that every bond is in a ring and that
    // RDKit✔️✔️:   // there is only one connected component in the molecular graph.
    // RDKit✔️✔️:   if (RDKit::MolOps::getMolFrags(mol).size() != 1) {
    // RDKit✔️✔️:     std::string msg =
    // RDKit✔️✔️:         "Template consists of multiple fragments, single fragment required: " +
    // RDKit✔️✔️:         smiles;
    // RDKit✔️✔️:     throw RDDepict::DepictException(msg);
    // RDKit✔️✔️:   }
    if runtime.fragment_count != 1 {
        return Err(Coordinate2DError::UnsupportedFeature(format!(
            "Template consists of multiple fragments, single fragment required: {smiles}"
        )));
    }

    // RDKit✔️✔️:   if (mol.getNumAtoms() == 1) {
    // RDKit✔️✔️:     std::string msg = "Template is not a ring system: " + smiles;
    // RDKit✔️✔️:     throw RDDepict::DepictException(msg);
    // RDKit✔️✔️:   }
    if runtime.graph.atom_queries.len() == 1 {
        return Err(Coordinate2DError::UnsupportedFeature(format!(
            "Template is not a ring system: {smiles}"
        )));
    }
    // RDKit✔️✔️:   auto ri = mol.getRingInfo();
    // RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
    // RDKit✔️✔️:     if (!ri->numBondRings(i)) {
    // RDKit✔️✔️:       std::string msg = "Template is not a ring system: " + smiles;
    // RDKit✔️✔️:       throw RDDepict::DepictException(msg);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if runtime.bond_ring_counts.iter().any(|count| *count == 0) {
        return Err(Coordinate2DError::UnsupportedFeature(format!(
            "Template is not a ring system: {smiles}"
        )));
    }
    Ok(())
}

fn for_each_rdkit_template_file_line<F>(
    template_path: &str,
    mut visit: F,
) -> Result<(), Coordinate2DError>
where
    F: FnMut(&str) -> Result<(), Coordinate2DError>,
{
    // RDKit✔️✔️: void CoordinateTemplates::loadTemplatesFromPath(
    // RDKit✔️✔️:     const std::string &templatePath,
    // RDKit✔️✔️:     std::unordered_map<unsigned int, std::vector<std::shared_ptr<RDKit::ROMol>>>
    // RDKit✔️✔️:         &templates) {
    // RDKit✔️✔️:   std::ifstream cxsmiles(templatePath);
    // RDKit✔️✔️:   if (!cxsmiles) {
    // RDKit✔️✔️:     std::string msg = "Could not open file " + templatePath;
    // RDKit✔️✔️:     throw RDDepict::DepictException(msg);
    // RDKit✔️✔️:   }
    let file = File::open(template_path).map_err(|_| {
        Coordinate2DError::UnsupportedFeature(format!("Could not open file {template_path}"))
    })?;

    // RDKit✔️✔️:   // Try loading templates in from this directory, if unsuccessful, keep current
    // RDKit✔️✔️:   // templates
    // RDKit✔️✔️:   std::string line;
    // RDKit✔️✔️:   while (std::getline(cxsmiles, line)) {
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let line = line.map_err(|error| {
            Coordinate2DError::UnsupportedFeature(format!(
                "RDKit template file read failed for {template_path}: {error}"
            ))
        })?;
        visit(&line)?;
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   cxsmiles.close();
    Ok(())
}

fn parse_rdkit_template_line_for_loading(
    template_path: &str,
    line: &str,
) -> Result<RdkitTemplateRuntimeModel, Coordinate2DError> {
    build_rdkit_template_runtime_model(line).map_err(|_| {
        Coordinate2DError::UnsupportedFeature(format!(
            "Could not load templates from {template_path}: Invalid smarts"
        ))
    })
}

fn load_rdkit_templates_from_path(
    template_path: &str,
    templates: &mut RdkitTemplateBuckets,
) -> Result<(), Coordinate2DError> {
    // RDKit✔️✔️: void CoordinateTemplates::loadTemplatesFromPath(
    // RDKit✔️✔️:     const std::string &templatePath,
    // RDKit✔️✔️:     std::unordered_map<unsigned int, std::vector<std::shared_ptr<RDKit::ROMol>>>
    // RDKit✔️✔️:         &templates) {
    for_each_rdkit_template_file_line(template_path, |line| {
        // RDKit✔️✔️:     RDKit::ROMol *mol_ptr = RDKit::SmartsToMol(line);
        // RDKit✔️✔️:     if (!mol_ptr) {
        // RDKit✔️✔️:       std::string msg =
        // RDKit✔️✔️:           "Could not load templates from " + templatePath + ": Invalid smarts";
        // RDKit✔️✔️:       cxsmiles.close();
        // RDKit✔️✔️:       throw RDDepict::DepictException(msg);
        // RDKit✔️✔️:     }
        let runtime = parse_rdkit_template_line_for_loading(template_path, line)?;
        // RDKit✔️❗:     std::shared_ptr<RDKit::ROMol> mol(mol_ptr);
        // RDKit✔️✔️:     // Initialize ring info using symmetrizeSSSR to match depictor ring counting
        // RDKit✔️✔️:     RDKit::VECT_INT_VECT arings;
        // RDKit✔️✔️:     RDKit::MolOps::symmetrizeSSSR(*mol, arings);
        // RDKit✔️✔️:     try {
        // RDKit✔️✔️:       assertValidTemplate(*mol, line);
        // RDKit✔️✔️:     } catch (RDDepict::DepictException &e) {
        // RDKit✔️✔️:       cxsmiles.close();
        // RDKit✔️✔️:       throw e;
        // RDKit✔️✔️:     }
        assert_valid_rdkit_template(&runtime, line)?;
        // RDKit✔️✔️:     templates[mol->getNumAtoms()].push_back(mol);
        templates
            .entry(runtime.graph.atom_queries.len())
            .or_default()
            .push(runtime);
        Ok(())
    })?;
    // RDKit✔️✔️:   cxsmiles.close();
    // RDKit✔️✔️: }
    Ok(())
}

static RDKIT_COORDINATE_TEMPLATE_REGISTRY: OnceLock<RwLock<RdkitCoordinateTemplateRegistry>> =
    OnceLock::new();

fn rdkit_coordinate_template_registry() -> &'static RwLock<RdkitCoordinateTemplateRegistry> {
    // RDKit✔️✔️:   static CoordinateTemplates template_mols;
    // RDKit✔️✔️:   return template_mols;
    RDKIT_COORDINATE_TEMPLATE_REGISTRY.get_or_init(|| {
        // RDKit✔️✔️: private:
        // RDKit✔️✔️:   CoordinateTemplates() { loadDefaultTemplates(); }
        RwLock::new(build_default_rdkit_coordinate_template_registry())
    })
}

fn rdkit_has_template_of_size(atom_count: usize) -> bool {
    // RDKit✔️✔️: bool hasTemplateOfSize(unsigned int atomCount) {
    // RDKit✔️✔️:   if (m_templates.find(atomCount) != m_templates.end()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    rdkit_coordinate_template_registry()
        .read()
        .expect("template registry lock poisoned")
        .templates
        .contains_key(&atom_count)
}

fn rdkit_matching_templates(atom_count: usize) -> Vec<RdkitTemplateRuntimeModel> {
    // RDKit✔️✔️: const std::vector<std::shared_ptr<RDKit::ROMol>> &getMatchingTemplates(
    // RDKit✔️✔️:     unsigned int atomCount) {
    // RDKit✔️✔️:   return m_templates[atomCount];
    // RDKit✔️✔️: }
    rdkit_coordinate_template_registry()
        .read()
        .expect("template registry lock poisoned")
        .templates
        .get(&atom_count)
        .cloned()
        .unwrap_or_default()
}

fn rdkit_default_template_smarts() -> impl Iterator<Item = &'static str> {
    RDKIT_TEMPLATE_SMARTS_SOURCE.lines().filter_map(|line| {
        let trimmed = line.trim();
        if !trimmed.starts_with('"') {
            return None;
        }
        let content = trimmed.strip_prefix('"')?;
        let end = content.rfind('"')?;
        Some(&content[..end])
    })
}

fn load_default_rdkit_ring_system_templates() {
    // RDKit✔️✔️: void loadDefaultTemplates() {
    // RDKit✔️✔️:   clearTemplates();
    // RDKit✔️✔️:   // load default templates into m_templates map by atom count
    // RDKit✔️✔️:   for (const auto &smarts : TEMPLATE_SMARTS) {
    // RDKit✔️✔️:     std::shared_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts));
    // RDKit✔️✔️:     if (!mol) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // Initialize ring info using symmetrizeSSSR to match depictor ring counting
    // RDKit✔️✔️:     RDKit::VECT_INT_VECT arings;
    // RDKit✔️✔️:     RDKit::MolOps::symmetrizeSSSR(*mol, arings);
    // RDKit✔️✔️:     m_templates[mol->getNumAtoms()].push_back(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let mut registry = rdkit_coordinate_template_registry()
        .write()
        .expect("template registry lock poisoned");
    registry.templates.clear();
    for smarts in rdkit_default_template_smarts() {
        let Ok(runtime) = build_rdkit_template_runtime_model(smarts) else {
            continue;
        };
        registry
            .templates
            .entry(runtime.graph.atom_queries.len())
            .or_default()
            .push(runtime);
    }
}

fn build_default_rdkit_coordinate_template_registry() -> RdkitCoordinateTemplateRegistry {
    let mut templates = RdkitTemplateBuckets::new();
    for smarts in rdkit_default_template_smarts() {
        let Ok(runtime) = build_rdkit_template_runtime_model(smarts) else {
            continue;
        };
        templates
            .entry(runtime.graph.atom_queries.len())
            .or_default()
            .push(runtime);
    }
    RdkitCoordinateTemplateRegistry { templates }
}

fn set_rdkit_ring_system_templates(template_path: &str) -> Result<(), Coordinate2DError> {
    // RDKit✔️✔️: void CoordinateTemplates::setRingSystemTemplates(
    // RDKit✔️✔️:     const std::string &templatePath) {
    // RDKit✔️✔️:   // Try loading templates in from this directory, if unsuccessful, keep current
    // RDKit✔️✔️:   // templates
    // RDKit✔️✔️:   std::unordered_map<unsigned int, std::vector<std::shared_ptr<RDKit::ROMol>>>
    // RDKit✔️✔️:       templates;
    // RDKit✔️✔️:   loadTemplatesFromPath(templatePath, templates);
    let mut templates = RdkitTemplateBuckets::new();
    load_rdkit_templates_from_path(template_path, &mut templates)?;
    // RDKit✔️✔️:   clearTemplates();
    // RDKit✔️✔️:   m_templates = std::move(templates);
    // RDKit✔️✔️: }
    rdkit_coordinate_template_registry()
        .write()
        .expect("template registry lock poisoned")
        .templates = templates;
    Ok(())
}

fn add_rdkit_ring_system_templates(template_path: &str) -> Result<(), Coordinate2DError> {
    // RDKit✔️✔️: void CoordinateTemplates::addRingSystemTemplates(
    // RDKit✔️✔️:     const std::string &templatePath) {
    // RDKit✔️✔️:   // Try loading templates in from this directory, if unsuccessful, keep current
    // RDKit✔️✔️:   // templates
    // RDKit✔️✔️:   std::unordered_map<unsigned int, std::vector<std::shared_ptr<RDKit::ROMol>>>
    // RDKit✔️✔️:       templates;
    // RDKit✔️✔️:   loadTemplatesFromPath(templatePath, templates);
    let mut templates = RdkitTemplateBuckets::new();
    load_rdkit_templates_from_path(template_path, &mut templates)?;
    // RDKit✔️✔️:   for (auto &kv : templates) {
    // RDKit✔️✔️:     m_templates[kv.first].insert(m_templates[kv.first].begin(),
    // RDKit✔️✔️:                                  kv.second.begin(), kv.second.end());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let mut registry = rdkit_coordinate_template_registry()
        .write()
        .expect("template registry lock poisoned");
    for (atom_count, mut incoming) in templates {
        let bucket = registry.templates.entry(atom_count).or_default();
        incoming.append(bucket);
        *bucket = incoming;
    }
    Ok(())
}

fn rdkit_load_default_ring_system_templates() {
    // RDKit✔️✔️: void loadDefaultRingSystemTemplates() {
    // RDKit✔️✔️:   CoordinateTemplates &coordinate_templates =
    // RDKit✔️✔️:       CoordinateTemplates::getRingSystemTemplates();
    // RDKit✔️✔️:   coordinate_templates.loadDefaultTemplates();
    // RDKit✔️✔️: }
    load_default_rdkit_ring_system_templates();
}

fn parse_template_cx_coordinate_block(
    cx_block: &str,
    atom_count: usize,
) -> Result<
    (
        Option<Vec<[f64; 2]>>,
        Option<Vec<[f64; 3]>>,
        Option<CoordinateDimension>,
    ),
    Coordinate2DError,
> {
    if !cx_block.starts_with("|(") || !cx_block.ends_with(")|") {
        return Err(Coordinate2DError::UnsupportedFeature(
            "RDKit template runtime model currently supports only a leading CX coordinate block"
                .to_string(),
        ));
    }
    let payload = &cx_block[2..cx_block.len() - 2];
    let mut coords_3d = vec![[0.0_f64, 0.0_f64, 0.0_f64]; atom_count];
    let mut saw_z_token = false;
    for (atom_idx, token) in payload.split(';').enumerate() {
        if atom_idx >= atom_count || token.is_empty() {
            continue;
        }
        let mut pieces = token.split(',');
        if let Some(x) = pieces.next()
            && !x.is_empty()
        {
            coords_3d[atom_idx][0] = x.parse::<f64>().map_err(|_| {
                Coordinate2DError::UnsupportedFeature(
                    "RDKit template CX coordinate parsing failed for x".to_string(),
                )
            })?;
        }
        if let Some(y) = pieces.next()
            && !y.is_empty()
        {
            coords_3d[atom_idx][1] = y.parse::<f64>().map_err(|_| {
                Coordinate2DError::UnsupportedFeature(
                    "RDKit template CX coordinate parsing failed for y".to_string(),
                )
            })?;
        }
        if let Some(z) = pieces.next()
            && !z.is_empty()
        {
            coords_3d[atom_idx][2] = z.parse::<f64>().map_err(|_| {
                Coordinate2DError::UnsupportedFeature(
                    "RDKit template CX coordinate parsing failed for z".to_string(),
                )
            })?;
            saw_z_token = true;
        }
    }
    let is_3d = saw_z_token && coords_3d.iter().any(|coord| coord[2].abs() > 1e-3);
    if is_3d {
        Ok((None, Some(coords_3d), Some(CoordinateDimension::ThreeD)))
    } else {
        let coords_2d = coords_3d.iter().map(|coord| [coord[0], coord[1]]).collect();
        Ok((Some(coords_2d), None, Some(CoordinateDimension::TwoD)))
    }
}

fn build_template_topology_probe_molecule(
    graph: &RdkitTemplateGraphModel,
) -> Result<Molecule, Coordinate2DError> {
    let mut builder = MoleculeBuilder::new();
    let atom_ids: Vec<_> = (0..graph.atom_queries.len())
        .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
        .collect();
    for bond in &graph.bonds {
        builder
            .add_bond(BondSpec::new(
                atom_ids[bond.begin_atom_idx],
                atom_ids[bond.end_atom_idx],
                bond_order_for_template_probe(&bond.query),
            ))
            .map_err(|error| {
                Coordinate2DError::UnsupportedFeature(format!(
                    "RDKit template topology probe build failed: {error}"
                ))
            })?;
    }
    builder.build().map_err(|error| {
        Coordinate2DError::UnsupportedFeature(format!(
            "RDKit template topology probe molecule build failed: {error}"
        ))
    })
}

fn bond_order_for_template_probe(query: &crate::QueryNode<crate::BondQueryPredicate>) -> BondOrder {
    match query {
        crate::QueryNode::Predicate(crate::BondQueryPredicate::Order(order)) => *order,
        crate::QueryNode::Predicate(crate::BondQueryPredicate::IsAromatic(true)) => {
            BondOrder::Aromatic
        }
        _ => BondOrder::Single,
    }
}

fn expand_template_smarts_bonds(smarts: &str) -> Result<Vec<TemplateGraphBond>, Coordinate2DError> {
    #[derive(Clone)]
    struct ParserState<'a> {
        chars: &'a [char],
        pos: usize,
        next_atom_idx: usize,
        bonds: Vec<TemplateGraphBond>,
        ring_open: BTreeMap<u8, (usize, crate::QueryNode<crate::BondQueryPredicate>)>,
    }

    fn bond_query_for_char(ch: char) -> crate::QueryNode<crate::BondQueryPredicate> {
        match ch {
            '-' => crate::QueryNode::Predicate(crate::BondQueryPredicate::Order(BondOrder::Single)),
            '=' => crate::QueryNode::Predicate(crate::BondQueryPredicate::Order(BondOrder::Double)),
            '#' => crate::QueryNode::Predicate(crate::BondQueryPredicate::Order(BondOrder::Triple)),
            ':' => crate::QueryNode::Predicate(crate::BondQueryPredicate::IsAromatic(true)),
            '~' => crate::QueryNode::Predicate(crate::BondQueryPredicate::Any),
            // RDKit source: SMARTS `/` and `\` preserve directional stereo on
            // the query bond but do not contribute an extra graph-level bond
            // predicate during default substructure matching.
            '/' | '\\' => crate::QueryNode::Predicate(crate::BondQueryPredicate::Any),
            _ => crate::QueryNode::Predicate(crate::BondQueryPredicate::Any),
        }
    }

    fn is_bond_char(ch: char) -> bool {
        matches!(ch, '-' | '=' | '#' | ':' | '~' | '/' | '\\')
    }

    fn is_organic_subset_symbol(chars: &[char], pos: usize) -> Option<usize> {
        let tail = chars.get(pos..)?;
        let two_char = if tail.len() >= 2 {
            match (tail[0], tail[1]) {
                ('C', 'l') | ('B', 'r') | ('S', 'i') | ('A', 's') | ('S', 'e') | ('T', 'e') => {
                    Some(2)
                }
                _ => None,
            }
        } else {
            None
        };
        two_char.or_else(|| {
            matches!(
                tail[0],
                '*' | 'B'
                    | 'C'
                    | 'N'
                    | 'O'
                    | 'S'
                    | 'P'
                    | 'F'
                    | 'I'
                    | 'H'
                    | 'R'
                    | 'X'
                    | 'D'
                    | 'v'
                    | 'V'
                    | 'r'
                    | 'u'
                    | 'A'
                    | 'T'
                    | 'Z'
                    | 'K'
                    | 'W'
                    | 'U'
                    | 'Y'
                    | 'G'
                    | 'L'
                    | 'J'
                    | 'E'
                    | 'M'
                    | 'Q'
                    | 'c'
                    | 'n'
                    | 'o'
                    | 's'
                    | 'p'
                    | 'a'
                    | 'b'
            )
            .then_some(1)
        })
    }

    fn parse_ring_number(chars: &[char], pos: &mut usize) -> Result<u8, Coordinate2DError> {
        if chars.get(*pos) == Some(&'%') {
            let d1 = chars.get(*pos + 1).copied();
            let d2 = chars.get(*pos + 2).copied();
            match (d1, d2) {
                (Some(c1), Some(c2)) if c1.is_ascii_digit() && c2.is_ascii_digit() => {
                    *pos += 3;
                    Ok(((c1 as u8 - b'0') * 10) + (c2 as u8 - b'0'))
                }
                _ => Err(Coordinate2DError::UnsupportedFeature(
                    "RDKit template SMARTS topology parser expected two digits after '%'"
                        .to_string(),
                )),
            }
        } else if let Some(ch) = chars.get(*pos).copied() {
            if ch.is_ascii_digit() {
                *pos += 1;
                Ok(ch as u8 - b'0')
            } else {
                Err(Coordinate2DError::UnsupportedFeature(
                    "RDKit template SMARTS topology parser expected a ring-closure digit"
                        .to_string(),
                ))
            }
        } else {
            Err(Coordinate2DError::UnsupportedFeature(
                "RDKit template SMARTS topology parser reached end while reading ring closure"
                    .to_string(),
            ))
        }
    }

    fn skip_atom(chars: &[char], pos: &mut usize) -> Result<(), Coordinate2DError> {
        match chars.get(*pos).copied() {
            Some('[') => {
                let mut depth = 1usize;
                *pos += 1;
                while let Some(ch) = chars.get(*pos).copied() {
                    *pos += 1;
                    match ch {
                        '[' => depth += 1,
                        ']' => {
                            depth -= 1;
                            if depth == 0 {
                                return Ok(());
                            }
                        }
                        _ => {}
                    }
                }
                Err(Coordinate2DError::UnsupportedFeature(
                    "RDKit template SMARTS topology parser found an unclosed bracket atom"
                        .to_string(),
                ))
            }
            Some(_) => {
                let Some(consumed) = is_organic_subset_symbol(chars, *pos) else {
                    return Err(Coordinate2DError::UnsupportedFeature(
                        "RDKit template SMARTS topology parser encountered an unsupported atom token"
                            .to_string(),
                    ));
                };
                *pos += consumed;
                Ok(())
            }
            None => Err(Coordinate2DError::UnsupportedFeature(
                "RDKit template SMARTS topology parser expected an atom token".to_string(),
            )),
        }
    }

    fn parse_atom(state: &mut ParserState<'_>) -> Result<usize, Coordinate2DError> {
        skip_atom(state.chars, &mut state.pos)?;
        let atom_idx = state.next_atom_idx;
        state.next_atom_idx += 1;
        Ok(atom_idx)
    }

    fn add_bond(
        state: &mut ParserState<'_>,
        begin_atom_idx: usize,
        end_atom_idx: usize,
        query: crate::QueryNode<crate::BondQueryPredicate>,
    ) {
        state.bonds.push(TemplateGraphBond {
            begin_atom_idx,
            end_atom_idx,
            query,
        });
    }

    fn parse_chain(
        state: &mut ParserState<'_>,
        mut current_atom_idx: usize,
    ) -> Result<usize, Coordinate2DError> {
        while state.pos < state.chars.len() {
            match state.chars[state.pos] {
                ')' => break,
                '(' => {
                    state.pos += 1;
                    let _ = parse_chain(state, current_atom_idx)?;
                    if state.chars.get(state.pos) != Some(&')') {
                        return Err(Coordinate2DError::UnsupportedFeature(
                            "RDKit template SMARTS topology parser expected ')' to close a branch"
                                .to_string(),
                        ));
                    }
                    state.pos += 1;
                }
                '.' => {
                    state.pos += 1;
                    current_atom_idx = parse_atom(state)?;
                }
                _ => {
                    let mut query = crate::QueryNode::Predicate(crate::BondQueryPredicate::Any);
                    if is_bond_char(state.chars[state.pos]) {
                        query = bond_query_for_char(state.chars[state.pos]);
                        state.pos += 1;
                    }
                    if state.pos >= state.chars.len() {
                        return Err(Coordinate2DError::UnsupportedFeature(
                            "RDKit template SMARTS topology parser ended after a bond token"
                                .to_string(),
                        ));
                    }
                    if state.chars[state.pos].is_ascii_digit() || state.chars[state.pos] == '%' {
                        let ring_number = parse_ring_number(state.chars, &mut state.pos)?;
                        if let Some((open_atom_idx, open_query)) =
                            state.ring_open.remove(&ring_number)
                        {
                            if open_query
                                != crate::QueryNode::Predicate(crate::BondQueryPredicate::Any)
                                && query
                                    != crate::QueryNode::Predicate(crate::BondQueryPredicate::Any)
                                && open_query != query
                            {
                                return Err(Coordinate2DError::UnsupportedFeature(
                                    "RDKit template SMARTS topology parser does not support mismatched explicit ring-closure bond queries"
                                        .to_string(),
                                ));
                            }
                            let resolved_query = if query
                                == crate::QueryNode::Predicate(crate::BondQueryPredicate::Any)
                            {
                                open_query
                            } else {
                                query
                            };
                            add_bond(state, open_atom_idx, current_atom_idx, resolved_query);
                        } else {
                            state
                                .ring_open
                                .insert(ring_number, (current_atom_idx, query));
                        }
                    } else {
                        let next_atom_idx = parse_atom(state)?;
                        add_bond(state, current_atom_idx, next_atom_idx, query);
                        current_atom_idx = next_atom_idx;
                    }
                }
            }
        }
        Ok(current_atom_idx)
    }

    let chars: Vec<char> = smarts.chars().collect();
    if chars.is_empty() {
        return Ok(Vec::new());
    }
    let mut state = ParserState {
        chars: &chars,
        pos: 0,
        next_atom_idx: 0,
        bonds: Vec::new(),
        ring_open: BTreeMap::new(),
    };
    let first_atom_idx = parse_atom(&mut state)?;
    let _ = parse_chain(&mut state, first_atom_idx)?;
    if !state.ring_open.is_empty() {
        return Err(Coordinate2DError::UnsupportedFeature(
            "RDKit template SMARTS topology parser found an unbalanced ring closure".to_string(),
        ));
    }
    Ok(state.bonds)
}

static PREFER_COORD_GEN: AtomicBool = AtomicBool::new(false);

pub(crate) fn set_prefer_coord_gen(value: bool) {
    PREFER_COORD_GEN.store(value, Ordering::Relaxed);
}

pub(crate) fn prefer_coord_gen() -> bool {
    PREFER_COORD_GEN.load(Ordering::Relaxed)
}

pub(crate) const fn is_coordgen_support_available() -> bool {
    false
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Compute2DCoordParameters<'a> {
    pub(crate) coord_map: Option<&'a BTreeMap<usize, [f64; 2]>>,
    pub(crate) canon_orient: bool,
    pub(crate) clear_confs: bool,
    pub(crate) n_flips_per_sample: u32,
    pub(crate) n_samples: u32,
    pub(crate) sample_seed: i32,
    pub(crate) permute_deg4_nodes: bool,
    pub(crate) force_rdkit: bool,
    pub(crate) use_ring_templates: bool,
}

impl Default for Compute2DCoordParameters<'_> {
    fn default() -> Self {
        Self {
            coord_map: None,
            canon_orient: false,
            clear_confs: true,
            n_flips_per_sample: 0,
            n_samples: 0,
            sample_seed: 0,
            permute_deg4_nodes: false,
            force_rdkit: false,
            use_ring_templates: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct With2DCoordinatesParams {
    pub canon_orient: bool,
    pub clear_confs: bool,
    pub n_flips_per_sample: u32,
    pub n_samples: u32,
    pub sample_seed: i32,
    pub permute_deg4_nodes: bool,
    pub force_rdkit: bool,
    pub use_ring_templates: bool,
}

impl Default for With2DCoordinatesParams {
    fn default() -> Self {
        Self {
            canon_orient: true,
            clear_confs: true,
            n_flips_per_sample: 0,
            n_samples: 0,
            sample_seed: 0,
            permute_deg4_nodes: false,
            force_rdkit: false,
            use_ring_templates: false,
        }
    }
}

impl With2DCoordinatesParams {
    #[must_use]
    pub(crate) fn as_compute_params(&self) -> Compute2DCoordParameters<'static> {
        Compute2DCoordParameters {
            coord_map: None,
            canon_orient: self.canon_orient,
            clear_confs: self.clear_confs,
            n_flips_per_sample: self.n_flips_per_sample,
            n_samples: self.n_samples,
            sample_seed: self.sample_seed,
            permute_deg4_nodes: self.permute_deg4_nodes,
            force_rdkit: self.force_rdkit,
            use_ring_templates: self.use_ring_templates,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Compute2DCoordsMimicDistMatParameters {
    pub(crate) canon_orient: bool,
    pub(crate) clear_confs: bool,
    pub(crate) weight_dist_mat: f64,
    pub(crate) n_flips_per_sample: u32,
    pub(crate) n_samples: u32,
    pub(crate) sample_seed: i32,
    pub(crate) permute_deg4_nodes: bool,
    pub(crate) force_rdkit: bool,
}

impl Default for Compute2DCoordsMimicDistMatParameters {
    fn default() -> Self {
        Self {
            canon_orient: true,
            clear_confs: true,
            weight_dist_mat: 0.5,
            n_flips_per_sample: 3,
            n_samples: 100,
            sample_seed: 25,
            permute_deg4_nodes: true,
            force_rdkit: false,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct ConstrainedDepictionParams {
    pub(crate) accept_failure: bool,
    pub(crate) force_rdkit: bool,
    pub(crate) allow_rgroups: bool,
    pub(crate) align_only: bool,
    pub(crate) adjust_molblock_wedging: bool,
    pub(crate) existing_conf_id: isize,
    pub(crate) use_ring_templates: bool,
}

impl Default for ConstrainedDepictionParams {
    fn default() -> Self {
        Self {
            accept_failure: false,
            force_rdkit: false,
            allow_rgroups: false,
            align_only: false,
            adjust_molblock_wedging: true,
            existing_conf_id: -1,
            use_ring_templates: false,
        }
    }
}

// ── RDKit-compatible math wrappers ──────────────────────────────────────────
// These force a specific math-lib code path that RDKit's C++ expects.
// ── Basic math wrappers (RDKit Numerics, Geometry/point.h) ─────────────────

/// RDKit❗✔️: cos wrapper — matches RDKit's cos usage in depiction geometry
fn rdkit_cos(x: f64) -> f64 {
    // RDKit depiction code calls <cmath> directly; use the same C libm entry
    // points here instead of Rust intrinsics to reduce last-bit drift.
    unsafe { c_libm_cos(x) }
}

/// RDKit❗✔️: sin wrapper — matches RDKit's sin usage in depiction geometry
fn rdkit_sin(x: f64) -> f64 {
    unsafe { c_libm_sin(x) }
}

/// RDKit❗✔️: acos wrapper — matches RDKit's acos usage in depiction geometry
fn rdkit_acos(x: f64) -> f64 {
    unsafe { c_libm_acos(x) }
}

/// RDKit❗✔️: sqrt wrapper — matches RDKit's sqrt usage in depiction geometry
fn rdkit_sqrt(x: f64) -> f64 {
    unsafe { c_libm_sqrt(x) }
}

/// RDKit❗✔️: 2D vector norm — RDKit Point2D::length() equivalent
fn norm(v: (f64, f64)) -> f64 {
    rdkit_sqrt(v.0 * v.0 + v.1 * v.1)
}

const RDKIT_ZERO_TOLERANCE: f64 = 1.0e-16;

/// RDKit✔️✔️: 2D vector normalization — RDKit Point2D::normalize() equivalent
fn normalize(v: (f64, f64)) -> (f64, f64) {
    let n = norm(v);
    if n < RDKIT_ZERO_TOLERANCE {
        panic!("Cannot normalize a zero length vector");
    }
    (v.0 / n, v.1 / n)
}

/// RDKit❗✔️: 2D vector rotation — RDKit Transform2D::Rotate equivalent
fn rotate(v: (f64, f64), angle: f64) -> (f64, f64) {
    let c = rdkit_cos(angle);
    let s = rdkit_sin(angle);
    (v.0 * c - v.1 * s, v.0 * s + v.1 * c)
}

/// RDKit✔️✔️: rotate point around center — RDKit Transform2D rotation
fn rotate_around(p: (f64, f64), center: (f64, f64), angle: f64) -> (f64, f64) {
    let trans = transform2d_set_transform_center_angle(center, angle);
    transform2d_point(p, trans)
}

/// RDKit✔️✔️: computeAngle() geometry helper used by EmbeddedFrag
fn compute_angle(center: (f64, f64), p1: (f64, f64), p2: (f64, f64)) -> f64 {
    let t1 = normalize((p1.0 - center.0, p1.1 - center.1));
    let t2 = normalize((p2.0 - center.0, p2.1 - center.1));
    let mut dot_prod = t1.0 * t2.0 + t1.1 * t2.1;
    if dot_prod < -1.0 {
        dot_prod = -1.0;
    } else if dot_prod > 1.0 {
        dot_prod = 1.0;
    }
    rdkit_acos(dot_prod)
}

/// RDKit❗✔️: computeNormal() geometry helper used by EmbeddedFrag
fn compute_normal(center: (f64, f64), nbr: (f64, f64)) -> (f64, f64) {
    normalize((center.1 - nbr.1, nbr.0 - center.0))
}

/// RDKit❗✔️: rotation direction — RDKit Transform2D orientation detection
fn rotation_dir(center: (f64, f64), loc1: (f64, f64), loc2: (f64, f64), rem_angle: f64) -> i32 {
    let pt1 = (loc1.0 - center.0, loc1.1 - center.1);
    let pt2 = (loc2.0 - center.0, loc2.1 - center.1);
    let cross = (pt1.0 * pt2.1 - pt1.1 * pt2.0) * (PI - rem_angle);
    if cross >= 0.0 { -1 } else { 1 }
}

// ── 2D transform helpers (RDKit EmbeddedFrag.cpp) ────────────────────────────

#[allow(dead_code)]
/// RDKit❗✔️: 3x3 matrix multiplication — RDKit Transform2D::Mul3 equivalent
fn transform2d_mul3(lhs: [f64; 9], rhs: [f64; 9]) -> [f64; 9] {
    let mut out = [0.0; 9];
    for i in 0..3 {
        let id_a = i * 3;
        for j in 0..3 {
            let id_c = id_a + j;
            out[id_c] = 0.0;
            for k in 0..3 {
                out[id_c] += lhs[id_a + k] * rhs[k * 3 + j];
            }
        }
    }
    out
}

#[derive(Clone, Copy)]
struct RdkitTransform2D {
    data: [f64; 9],
}

impl RdkitTransform2D {
    fn identity() -> Self {
        let mut out = Self { data: [0.0; 9] };
        out.data[0] = 1.0;
        out.data[4] = 1.0;
        out.data[8] = 1.0;
        out
    }

    fn assign(&mut self, other: Self) {
        self.data = other.data;
    }

    fn mul_assign(&mut self, rhs: Self) {
        self.data = transform2d_mul3(self.data, rhs.data);
    }

    fn set_translation(&mut self, pt: (f64, f64)) {
        let mut i = 2usize;
        self.data[i] = pt.0;
        i += 3;
        self.data[i] = pt.1;
        i += 3;
        self.data[i] = 1.0;
    }

    fn transform_point(&self, pt: (f64, f64)) -> (f64, f64) {
        (
            self.data[0] * pt.0 + self.data[1] * pt.1 + self.data[2],
            self.data[3] * pt.0 + self.data[4] * pt.1 + self.data[5],
        )
    }

    fn set_transform_center_angle(pt: (f64, f64), angle: f64) -> Self {
        let mut this = Self::identity();
        let mut trans1 = Self::identity();
        trans1.set_translation((-pt.0, -pt.1));
        this.data[0] = rdkit_cos(angle);
        this.data[1] = -rdkit_sin(angle);
        this.data[3] = rdkit_sin(angle);
        this.data[4] = rdkit_cos(angle);
        this.mul_assign(trans1);

        let mut trans2 = Self::identity();
        trans2.set_translation(pt);
        trans2.mul_assign(this);
        this.assign(trans2);
        this
    }

    fn set_transform_two_point(
        ref1: (f64, f64),
        ref2: (f64, f64),
        pt1: (f64, f64),
        pt2: (f64, f64),
    ) -> Self {
        let rvec = (ref2.0 - ref1.0, ref2.1 - ref1.1);
        let pvec = (pt2.0 - pt1.0, pt2.1 - pt1.1);
        let dp = rvec.0 * pvec.0 + rvec.1 * pvec.1;
        let lp = norm(rvec) * norm(pvec);
        if lp <= 0.0 {
            return Self::identity();
        }

        let mut cval = dp / lp;
        if cval < -1.0 {
            cval = -1.0;
        } else if cval > 1.0 {
            cval = 1.0;
        }
        let mut ang = rdkit_acos(cval);
        let cross = pvec.0 * rvec.1 - pvec.1 * rvec.0;
        if cross < 0.0 {
            ang *= -1.0;
        }

        let mut this = Self::identity();
        this.data[0] = rdkit_cos(ang);
        this.data[1] = -rdkit_sin(ang);
        this.data[3] = rdkit_sin(ang);
        this.data[4] = rdkit_cos(ang);

        let npt1 = this.transform_point(pt1);
        this.data[2] = ref1.0 - npt1.0;
        this.data[5] = ref1.1 - npt1.1;
        this
    }

    fn to_affine(self) -> [f64; 6] {
        [
            self.data[0],
            self.data[1],
            self.data[2],
            self.data[3],
            self.data[4],
            self.data[5],
        ]
    }
}

#[allow(dead_code)]
/// RDKit❗✔️: convert 3×3 matrix to affine 2×6 — RDKit Transform2D::toAffine
fn transform2d_to_affine(data: [f64; 9]) -> [f64; 6] {
    [data[0], data[1], data[2], data[3], data[4], data[5]]
}

#[allow(dead_code)]
/// RDKit❗✔️: set rotation+translation from center+angle — RDKit Transform2D::setTransform
fn transform2d_set_transform_center_angle(pt: (f64, f64), angle: f64) -> [f64; 6] {
    RdkitTransform2D::set_transform_center_angle(pt, angle).to_affine()
}

/// RDKit❗✔️: apply affine transform to point — RDKit Transform2D::transformPoint
fn transform2d_point(pt: (f64, f64), data: [f64; 6]) -> (f64, f64) {
    (
        data[0] * pt.0 + data[1] * pt.1 + data[2],
        data[3] * pt.0 + data[4] * pt.1 + data[5],
    )
}

/// RDKit❗✔️: two-point transform — RDKit Transform2D::setTransform mapping p1→p1'
fn transform2d_set_transform_two_point(
    ref1: (f64, f64),
    ref2: (f64, f64),
    pt1: (f64, f64),
    pt2: (f64, f64),
) -> [f64; 6] {
    RdkitTransform2D::set_transform_two_point(ref1, ref2, pt1, pt2).to_affine()
}

fn rdkit_embed_ring(ring: &[usize]) -> RdkitIntPoint2DMap {
    // RDKit✔️✔️: RDGeom::INT_POINT2D_MAP embedRing(const RDKit::INT_VECT &ring) {
    // RDKit✔️✔️:   // The process here is very straight forward
    // RDKit✔️✔️:   // we take the center of the ring to lies at the origin put the first
    // RDKit✔️✔️:   // point at the origin and then sweep
    // RDKit✔️✔️:   // anticlockwise so by an angle A = 360/n for the next point
    // RDKit✔️✔️:   // the length of the arm (l) we want to sweep is easy to compute given the
    // RDKit✔️✔️:   // bond length (b) we want to use for each bond in the ring (for now
    // RDKit✔️✔️:   // we will assume that this bond length is the same for all bonds in the ring
    // RDKit✔️✔️:   //  l = b/sqrt(2*(1 - cos(A))
    // RDKit✔️✔️:   // the above formula derives from the triangle formula, where side 'c' is
    // RDKit✔️✔️:   // given
    // RDKit✔️✔️:   // in terms of sides 'a' and 'b' as
    // RDKit✔️✔️:   // c = a^2 + b^2 - 2.a.b.cos(A)
    // RDKit✔️✔️:   // where A is the angle between a and b
    // RDKit✔️✔️:   // compute the sweep angle
    // RDKit✔️✔️:   unsigned int na = ring.size();
    let na = ring.len();
    // RDKit✔️✔️:   double ang = 2 * M_PI / na;
    let ang = 2.0 * PI / na as f64;
    // RDKit✔️✔️:   // compute the arm length
    // RDKit✔️✔️:   double al = BOND_LEN / (sqrt(2 * (1 - cos(ang))));
    let al = 1.5 / rdkit_sqrt(2.0 * (1.0 - rdkit_cos(ang)));
    // RDKit✔️✔️:   RDGeom::INT_POINT2D_MAP res;
    let mut res = RdkitIntPoint2DMap::new();
    // RDKit✔️✔️:   for (unsigned int i = 0; i < na; ++i) {
    // RDKit✔️✔️:     auto x = al * cos(i * ang);
    // RDKit✔️✔️:     auto y = al * sin(i * ang);
    // RDKit✔️✔️:     RDGeom::Point2D loc(x, y);
    // RDKit✔️✔️:     res[ring[i]] = loc;
    // RDKit✔️✔️:   }
    for (i, atom_idx) in ring.iter().copied().enumerate() {
        let x = al * rdkit_cos(i as f64 * ang);
        let y = al * rdkit_sin(i as f64 * ang);
        res.insert(atom_idx, (x, y));
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    res
}

fn rdkit_transform_points(nring_cor: &mut RdkitIntPoint2DMap, trans: [f64; 6]) {
    // RDKit✔️✔️: void transformPoints(RDGeom::INT_POINT2D_MAP &nringCor,
    // RDKit✔️✔️:                      const RDGeom::Transform2D &trans) {
    // RDKit✔️✔️:   std::for_each(nringCor.begin(), nringCor.end(),
    // RDKit✔️✔️:                 [&trans](auto &elem) { trans.TransformPoint(elem.second); });
    // RDKit✔️✔️: }
    for point in nring_cor.values_mut() {
        *point = transform2d_point(*point, trans);
    }
}

fn rdkit_compute_bisect_point(
    rcr: (f64, f64),
    ang: f64,
    nb1: (f64, f64),
    nb2: (f64, f64),
) -> (f64, f64) {
    // RDKit✔️✔️: RDGeom::Point2D computeBisectPoint(const RDGeom::Point2D &rcr, double ang,
    // RDKit✔️✔️:                                    const RDGeom::Point2D &nb1,
    // RDKit✔️✔️:                                    const RDGeom::Point2D &nb2) {
    // RDKit✔️✔️:   RDGeom::Point2D cloc = nb1;
    let mut cloc = nb1;
    // RDKit✔️✔️:   cloc += nb2;
    cloc.0 += nb2.0;
    cloc.1 += nb2.1;
    // RDKit✔️✔️:   cloc *= 0.5;
    cloc.0 *= 0.5;
    cloc.1 *= 0.5;
    // RDKit✔️✔️:   if (ang > M_PI) {
    if ang > PI {
        // RDKit✔️✔️:     // invert the cloc
        // RDKit✔️✔️:     cloc -= rcr;
        cloc.0 -= rcr.0;
        cloc.1 -= rcr.1;
        // RDKit✔️✔️:     cloc *= -1.0;
        cloc.0 *= -1.0;
        cloc.1 *= -1.0;
        // RDKit✔️✔️:     cloc += rcr;
        cloc.0 += rcr.0;
        cloc.1 += rcr.1;
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return cloc;
    // RDKit✔️✔️: }
    cloc
}

fn rdkit_reflect_point(point: (f64, f64), loc1: (f64, f64), loc2: (f64, f64)) -> (f64, f64) {
    // RDKit✔️✔️: RDGeom::Point2D reflectPoint(const RDGeom::Point2D &point,
    // RDKit✔️✔️:                              const RDGeom::Point2D &loc1,
    // RDKit✔️✔️:                              const RDGeom::Point2D &loc2) {
    // RDKit✔️✔️:   RDGeom::Point2D org(0.0, 0.0);
    let org = (0.0, 0.0);
    // RDKit✔️✔️:   RDGeom::Point2D xaxis(1.0, 0.0);
    let xaxis = (1.0, 0.0);
    // RDKit✔️✔️:   RDGeom::Point2D cent = (loc1 + loc2);
    let mut cent = (loc1.0 + loc2.0, loc1.1 + loc2.1);
    // RDKit✔️✔️:   cent *= 0.5;
    cent.0 *= 0.5;
    cent.1 *= 0.5;

    // RDKit✔️✔️:   RDGeom::Transform2D trans;
    // RDKit✔️✔️:   trans.SetTransform(org, xaxis, cent, loc1);
    let trans = transform2d_set_transform_two_point(org, xaxis, cent, loc1);
    // RDKit✔️✔️:   /// reverse transform
    // RDKit✔️✔️:   RDGeom::Transform2D itrans;
    // RDKit✔️✔️:   itrans.SetTransform(cent, loc1, org, xaxis);
    let itrans = transform2d_set_transform_two_point(cent, loc1, org, xaxis);

    // RDKit✔️✔️:   RDGeom::INT_POINT2D_MAP_I nci;
    // RDKit✔️✔️:   RDGeom::Point2D res;
    // RDKit✔️✔️:   res = point;
    // RDKit✔️✔️:   trans.TransformPoint(res);
    let mut res = transform2d_point(point, trans);
    // RDKit✔️✔️:   res.y = -res.y;
    res.1 = -res.1;
    // RDKit✔️✔️:   itrans.TransformPoint(res);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    transform2d_point(res, itrans)
}

fn rdkit_reflect_points(coord_map: &mut RdkitIntPoint2DMap, loc1: (f64, f64), loc2: (f64, f64)) {
    // RDKit✔️✔️: void reflectPoints(RDGeom::INT_POINT2D_MAP &coordMap,
    // RDKit✔️✔️:                    const RDGeom::Point2D &loc1, const RDGeom::Point2D &loc2) {
    // RDKit✔️✔️:   std::for_each(coordMap.begin(), coordMap.end(), [&loc1, &loc2](auto &elem) {
    // RDKit✔️✔️:     reflectPoint(elem.second, loc1, loc2);
    // RDKit✔️✔️:   });
    // RDKit✔️✔️: }
    for point in coord_map.values_mut() {
        *point = rdkit_reflect_point(*point, loc1, loc2);
    }
}

// ── Types ────────────────────────────────────────────────────────────────────

// BEGIN RDKIT CPP TYPE: ConjugHybrid.h HybridizationType
// RDKit✔️✔️: enum HybridizationType { UNSPECIFIED, S, SP, SP2, SP3, SP3D, SP3D2 };
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum RdkitHybridization {
    Unspecified,
    S,
    Sp,
    Sp2,
    Sp3,
    Sp3d,
    Sp3d2,
}

// BEGIN RDKIT CPP STRUCT: EmbeddedFrag.h EmbeddedAtom
// RDKit✔️✔️: struct EmbeddedAtom { Point2D loc; Point2D normal; bool ccw; ... }
#[derive(Clone, Debug)]
struct TreeEmbeddedAtom {
    loc: (f64, f64),
    normal: (f64, f64),
    ccw: bool,
    cis_trans_nbr: Option<usize>,
    angle: f64,
    nbr1: Option<usize>,
    nbr2: Option<usize>,
    rot_dir: i32,
    pending: Vec<usize>,
    d_density: f64,
    df_fixed: bool,
}

// BEGIN RDKIT CPP STRUCT: EmbeddedFrag.h EmbeddedFrag
// RDKit✔️❌: class EmbeddedFrag { INT_EATOM_MAP d_eatoms; INT_LIST d_attachPts; ... }
#[derive(Clone, Debug, Default)]
struct RdkitEmbeddedFrag {
    eatoms: BTreeMap<usize, TreeEmbeddedAtom>,
    attach_pts: VecDeque<usize>,
    done: bool,
}

impl RdkitEmbeddedFrag {
    fn size(&self) -> usize {
        self.eatoms.len()
    }

    fn is_done(&self) -> bool {
        self.done
    }

    fn mark_done(&mut self) {
        self.done = true;
    }

    fn get_embedded_atoms(&self) -> &BTreeMap<usize, TreeEmbeddedAtom> {
        &self.eatoms
    }

    // RDKit✔️✔️: EmbeddedFrag::EmbeddedFrag(unsigned int aid, const RDKit::ROMol *mol) { ... }
    fn from_single_atom(
        aid: usize,
        atoms: &[Atom],
        bonds: &[Bond],
        adjacency: &[Vec<usize>],
        degree: &[usize],
        cip_ranks: &[u32],
    ) -> Self {
        let mut frag = Self::default();
        frag.eatoms.insert(
            aid,
            TreeEmbeddedAtom {
                loc: (0.0, 0.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        frag.done = false;
        frag.update_new_neighs(aid, atoms, bonds, adjacency, degree, cip_ranks);
        frag
    }

    // RDKit✔️✔️: EmbeddedFrag::EmbeddedFrag(const RDKit::ROMol *mol,
    // RDKit✔️✔️:                            const RDGeom::INT_POINT2D_MAP &coordMap) { ... }
    fn from_coord_map(
        coord_map: &BTreeMap<usize, (f64, f64)>,
        atoms: &[Atom],
        bonds: &[Bond],
        adjacency: &[Vec<usize>],
        degree: &[usize],
        cip_ranks: &[u32],
        atom_ring_counts: &[usize],
    ) -> Self {
        let mut frag = Self::default();
        for (&aid, &loc) in coord_map {
            frag.eatoms.insert(
                aid,
                TreeEmbeddedAtom {
                    loc,
                    normal: (0.0, 0.0),
                    ccw: true,
                    cis_trans_nbr: None,
                    angle: -1.0,
                    nbr1: None,
                    nbr2: None,
                    rot_dir: 0,
                    pending: Vec::new(),
                    d_density: -1.0,
                    df_fixed: true,
                },
            );
            frag.done = false;
        }
        frag.setup_new_neighs(atoms, bonds, adjacency, degree, cip_ranks);
        frag.setup_attachment_points(adjacency, atom_ring_counts);
        frag
    }

    // RDKit❗❌: EmbeddedFrag::EmbeddedFrag(const RDKit::ROMol *mol,
    // RDKit❗❌:                            const RDKit::VECT_INT_VECT &fusedRings,
    // RDKit❗❌:                            bool useRingTemplates) { ... }
    fn from_fused_ring_atom_order(
        ring_atom_coords: &BTreeMap<usize, (f64, f64)>,
        atoms: &[Atom],
        bonds: &[Bond],
        adjacency: &[Vec<usize>],
        degree: &[usize],
        cip_ranks: &[u32],
        atom_ring_counts: &[usize],
    ) -> Self {
        let mut frag = Self::from_coord_map(
            ring_atom_coords,
            atoms,
            bonds,
            adjacency,
            degree,
            cip_ranks,
            atom_ring_counts,
        );
        frag.done = false;
        frag
    }

    // RDKit✔️✔️: EmbeddedFrag::EmbeddedFrag(const RDKit::Bond *dblBond) { ... }
    fn from_double_bond(bond: &Bond) -> Option<Self> {
        if bond.order() != BondOrder::Double {
            return None;
        }
        let stereo_atoms = bond.stereo_atoms()?;
        let stereo = bond.stereo();
        if matches!(stereo, BondStereo::Any | BondStereo::None) {
            return None;
        }
        let begin = bond.begin().index();
        let end = bond.end().index();
        let mut frag = Self::default();
        frag.eatoms.insert(
            begin,
            TreeEmbeddedAtom {
                loc: (0.0, 0.0),
                normal: (0.0, -1.0),
                ccw: false,
                cis_trans_nbr: Some(stereo_atoms[0].index()),
                angle: -1.0,
                nbr1: Some(end),
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        let (end_normal, end_ccw) = if matches!(stereo, BondStereo::Z | BondStereo::Cis) {
            ((0.0, -1.0), true)
        } else {
            ((0.0, 1.0), false)
        };
        frag.eatoms.insert(
            end,
            TreeEmbeddedAtom {
                loc: (BOND_LEN, 0.0),
                normal: end_normal,
                ccw: end_ccw,
                cis_trans_nbr: Some(stereo_atoms[1].index()),
                angle: -1.0,
                nbr1: Some(begin),
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        frag.done = false;
        Some(frag)
    }

    // RDKit✔️✔️: int EmbeddedFrag::findNumNeigh(const RDGeom::Point2D &pt, double radius) { ... }
    fn find_num_neigh(&self, pt: (f64, f64), radius: f64) -> i32 {
        self.eatoms
            .values()
            .filter(|st| norm((st.loc.0 - pt.0, st.loc.1 - pt.1)) < radius)
            .count() as i32
    }

    // RDKit✔️✔️: void EmbeddedFrag::computeNbrsAndAng(unsigned int aid,
    // RDKit✔️✔️:                                      const RDKit::INT_VECT &doneNbrs) { ... }
    fn compute_nbrs_and_ang(
        &mut self,
        aid: usize,
        done_nbrs: &[usize],
        atom_ring_counts: &[usize],
    ) {
        // RDKit✔️✔️:   auto winner = anglePairs.back();
        // RDKit✔️✔️:   for (auto pr : boost::adaptors::reverse(anglePairs)) {
        // RDKit✔️✔️:     if ((dp_mol->getRingInfo()->numAtomRings(pr.second.first) <= 1) &&
        // RDKit✔️✔️:         (dp_mol->getRingInfo()->numAtomRings(pr.second.second) <= 1)) {
        // RDKit✔️✔️:       winner = pr;
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        if done_nbrs.len() < 3 {
            return;
        }
        let mut angle_pairs = Vec::<(f64, (usize, usize))>::new();
        // RDKit✔️✔️:   for (auto nbi1 = doneNbrs.begin(); nbi1 != doneNbrs.end(); ++nbi1) {
        // RDKit✔️✔️:     auto nbi3 = nbi1;
        // RDKit✔️✔️:     for (auto nbi2 = nbi3++; nbi2 != doneNbrs.end(); ++nbi2) {
        // RDKit✔️✔️:       ang = computeAngle(d_eatoms[aid].loc, d_eatoms[*nbi1].loc,
        // RDKit✔️✔️:                          d_eatoms[*nbi2].loc);
        // RDKit✔️✔️:       auto nbrPair = std::make_pair((*nbi1), (*nbi2));
        // RDKit✔️✔️:       anglePairs.emplace_back(ang, nbrPair);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        for i in 0..done_nbrs.len() {
            for j in i..done_nbrs.len() {
                let a = done_nbrs[i];
                let b = done_nbrs[j];
                angle_pairs.push((
                    compute_angle(
                        self.eatoms[&aid].loc,
                        self.eatoms[&a].loc,
                        self.eatoms[&b].loc,
                    ),
                    (a, b),
                ));
            }
        }
        angle_pairs.sort_by(|x, y| x.0.total_cmp(&y.0));
        let Some(mut winner) = angle_pairs.last().copied() else {
            return;
        };
        for pair in angle_pairs.iter().rev() {
            let (a, b) = pair.1;
            if atom_ring_counts.get(a).copied().unwrap_or(0) <= 1
                && atom_ring_counts.get(b).copied().unwrap_or(0) <= 1
            {
                winner = *pair;
                break;
            }
        }
        let (w_ang, (wnb1, wnb2)) = winner;
        let mut nb1 = wnb1;
        let mut nb2 = wnb2;
        for &(_, (a, b)) in angle_pairs.iter() {
            if wnb1 == a {
                nb2 = wnb1;
                nb1 = b;
                break;
            } else if wnb1 == b {
                nb2 = wnb1;
                nb1 = a;
                break;
            } else if wnb2 == a {
                nb2 = wnb2;
                nb1 = b;
                break;
            } else if wnb2 == b {
                nb2 = wnb2;
                nb1 = a;
                break;
            }
        }
        let rot_dir = rotation_dir(
            self.eatoms[&aid].loc,
            self.eatoms[&nb1].loc,
            self.eatoms[&nb2].loc,
            w_ang,
        );
        if let Some(st) = self.eatoms.get_mut(&aid) {
            st.rot_dir = rot_dir;
            st.nbr1 = Some(nb1);
            st.nbr2 = Some(nb2);
            st.angle = 2.0 * PI - w_ang;
        }
    }

    // RDKit✔️✔️: void EmbeddedFrag::updateNewNeighs(unsigned int aid) { ... }
    fn update_new_neighs(
        &mut self,
        aid: usize,
        atoms: &[Atom],
        bonds: &[Bond],
        adjacency: &[Vec<usize>],
        degree: &[usize],
        cip_ranks: &[u32],
    ) {
        let mut heavy = Vec::new();
        let mut hydrogens = Vec::new();
        for &nb in &adjacency[aid] {
            if self.eatoms.contains_key(&nb) {
                continue;
            }
            if atoms[nb].atomic_number() != 1 {
                heavy.push(nb);
            } else {
                hydrogens.push(nb);
            }
        }
        heavy.extend(hydrogens);
        if !heavy.is_empty() && (degree[aid] < 4 || heavy.len() < 3) {
            heavy = rdkit_rank_atoms_by_rank(atoms, &heavy, degree, true);
        } else if degree[aid] >= 4 && heavy.len() >= 3 {
            heavy = rdkit_set_nbr_order(aid, &heavy, atoms, bonds, adjacency, degree);
        }
        if let Some(st) = self.eatoms.get_mut(&aid) {
            st.pending = heavy;
            if !st.pending.is_empty() && !self.attach_pts.iter().any(|&x| x == aid) {
                self.attach_pts.push_back(aid);
            }
        }
    }

    // RDKit✔️✔️: void EmbeddedFrag::setupNewNeighs() { ... }
    fn setup_new_neighs(
        &mut self,
        atoms: &[Atom],
        bonds: &[Bond],
        adjacency: &[Vec<usize>],
        degree: &[usize],
        cip_ranks: &[u32],
    ) {
        self.attach_pts.clear();
        let atom_ids: Vec<usize> = self.eatoms.keys().copied().collect();
        for aid in atom_ids {
            self.update_new_neighs(aid, atoms, bonds, adjacency, degree, cip_ranks);
        }
        let mut ranked = self.attach_pts.iter().copied().collect::<Vec<_>>();
        rdkit_rank_atoms_by_rank_into(atoms, &mut ranked, degree, true);
        self.attach_pts = ranked.into();
    }

    // RDKit✔️✔️: int EmbeddedFrag::findNeighbor(unsigned int aid) { ... }
    fn find_neighbor(&self, aid: usize, adjacency: &[Vec<usize>]) -> Option<usize> {
        adjacency[aid]
            .iter()
            .copied()
            .find(|nb| self.eatoms.contains_key(nb))
    }

    // RDKit✔️✔️: void EmbeddedFrag::setupAttachmentPoints() { ... }
    fn setup_attachment_points(&mut self, adjacency: &[Vec<usize>], atom_ring_counts: &[usize]) {
        let attach_ids: Vec<usize> = self.attach_pts.iter().copied().collect();
        for dai in attach_ids {
            let enbrs = self
                .eatoms
                .get(&dai)
                .map(|st| st.pending.clone())
                .unwrap_or_default();
            let done_nbrs: Vec<usize> = adjacency[dai]
                .iter()
                .copied()
                .filter(|nb| !enbrs.contains(nb) && self.eatoms.contains_key(nb))
                .collect();
            if done_nbrs.is_empty() {
                if let Some(st) = self.eatoms.get_mut(&dai) {
                    st.normal = (1.0, 0.0);
                    st.angle = -1.0;
                }
            } else if done_nbrs.len() == 1 {
                let nbid = done_nbrs[0];
                let nb_loc = self.eatoms[&nbid].loc;
                if let Some(st) = self.eatoms.get_mut(&dai) {
                    st.nbr1 = Some(nbid);
                    st.normal = compute_normal(st.loc, nb_loc);
                }
            } else if done_nbrs.len() == 2 {
                let nb1 = done_nbrs[0];
                let nb2 = done_nbrs[1];
                let loc = self.eatoms[&dai].loc;
                let ang = compute_angle(loc, self.eatoms[&nb1].loc, self.eatoms[&nb2].loc);
                if let Some(st) = self.eatoms.get_mut(&dai) {
                    st.nbr1 = Some(nb1);
                    st.nbr2 = Some(nb2);
                    st.angle = ang;
                }
            } else {
                self.compute_nbrs_and_ang(dai, &done_nbrs, atom_ring_counts);
            }
        }
    }

    // RDKit✔️✔️: void EmbeddedFrag::addAtomToAtomWithAng(unsigned int aid,
    // RDKit✔️✔️:                                         unsigned int toAid) { ... }
    fn add_atom_to_atom_with_ang(&mut self, aid: usize, to_aid: usize) -> Option<()> {
        let ref_state = self.eatoms.get(&to_aid)?.clone();
        let nnbr = ref_state.pending.len() as f64;
        let rem_angle = 2.0 * PI - ref_state.angle;
        let mut curr_angle = rem_angle / (1.0 + nnbr);
        if let Some(st) = self.eatoms.get_mut(&to_aid) {
            st.angle += curr_angle;
        }

        let nb1 = self.eatoms.get(&ref_state.nbr1?)?.loc;
        let nb2 = self.eatoms.get(&ref_state.nbr2?)?.loc;
        if self.eatoms.get(&to_aid)?.rot_dir == 0 {
            let rd = rotation_dir(ref_state.loc, nb1, nb2, rem_angle);
            if let Some(st) = self.eatoms.get_mut(&to_aid) {
                st.rot_dir = rd;
            }
        }
        curr_angle *= self.eatoms.get(&to_aid)?.rot_dir as f64;
        let mut curr_loc = rotate_around(nb2, ref_state.loc, curr_angle);
        if rem_angle.abs() - PI < 1e-3 {
            let curr_loc2 = rotate_around(nb2, ref_state.loc, -curr_angle);
            if self.find_num_neigh(curr_loc, 0.5) > self.find_num_neigh(curr_loc2, 0.5) {
                curr_loc = curr_loc2;
            }
        }
        if let Some(st) = self.eatoms.get_mut(&to_aid) {
            st.nbr2 = Some(aid);
        }

        let tpt = (curr_loc.0 - ref_state.loc.0, curr_loc.1 - ref_state.loc.1);
        let probe_normal = (-tpt.1, tpt.0);
        let tp1 = (curr_loc.0 + probe_normal.0, curr_loc.1 + probe_normal.1);
        let tp2 = (curr_loc.0 - probe_normal.0, curr_loc.1 - probe_normal.1);
        let nccw = self.find_num_neigh(tp1, 2.5);
        let ncw = self.find_num_neigh(tp2, 2.5);
        let mut normal = normalize(probe_normal);
        let (ccw, out_normal) = if nccw < ncw {
            (false, normal)
        } else {
            normal = (-normal.0, -normal.1);
            (true, normal)
        };
        self.eatoms.insert(
            aid,
            TreeEmbeddedAtom {
                loc: curr_loc,
                normal: out_normal,
                ccw,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: Some(to_aid),
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        Some(())
    }

    // RDKit✔️✔️: void EmbeddedFrag::addAtomToAtomWithNoAng(unsigned int aid,
    // RDKit✔️✔️:                                           unsigned int toAid) { ... }
    fn add_atom_to_atom_with_no_ang(
        &mut self,
        aid: usize,
        to_aid: usize,
        atoms: &[Atom],
        bonds: &[Bond],
        degree: &[usize],
    ) -> Option<()> {
        let ref_state = self.eatoms.get(&to_aid)?.clone();
        let mut ref_atom_ccw = ref_state.ccw;
        let mut curr_loc = ref_state.normal;
        if ref_state.cis_trans_nbr.is_some_and(|ct| ct != aid) {
            ref_atom_ccw = !ref_atom_ccw;
            curr_loc = (-curr_loc.0, -curr_loc.1);
        }

        let hybridizations = rdkit_hybridizations_for_depict(atoms, bonds, degree).ok()?;
        let mut angle = compute_sub_angle(degree[to_aid], hybridizations[to_aid]);
        let mut flip_norm = false;
        if self.eatoms.get(&to_aid)?.nbr1.is_some() {
            if let Some(st) = self.eatoms.get_mut(&to_aid) {
                st.angle = angle;
                st.nbr2 = Some(aid);
            }
        } else {
            let norm = self.eatoms.get(&to_aid)?.normal;
            let rot = rotate(norm, angle);
            if let Some(st) = self.eatoms.get_mut(&to_aid) {
                st.normal = rot;
                st.nbr1 = Some(aid);
            }
            flip_norm = true;
        }

        angle -= PI / 2.0;
        if !ref_atom_ccw {
            angle *= -1.0;
        }
        curr_loc = rotate(curr_loc, angle);
        curr_loc.0 *= BOND_LEN;
        curr_loc.1 *= BOND_LEN;
        curr_loc.0 += ref_state.loc.0;
        curr_loc.1 += ref_state.loc.1;

        let tpt = (ref_state.loc.0 - curr_loc.0, ref_state.loc.1 - curr_loc.1);
        let mut normal = (-tpt.1, tpt.0);
        if ref_atom_ccw ^ flip_norm {
            normal = (-normal.0, -normal.1);
        }
        normal = normalize(normal);
        self.eatoms.insert(
            aid,
            TreeEmbeddedAtom {
                loc: curr_loc,
                normal,
                ccw: (!ref_atom_ccw) ^ flip_norm,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: Some(to_aid),
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        Some(())
    }

    // RDKit✔️✔️: void EmbeddedFrag::addNonRingAtom(unsigned int aid,
    // RDKit✔️✔️:                                    unsigned int toAid) { ... }
    fn add_non_ring_atom(
        &mut self,
        aid: usize,
        to_aid: usize,
        atoms: &[Atom],
        bonds: &[Bond],
        adjacency: &[Vec<usize>],
        degree: &[usize],
        cip_ranks: &[u32],
    ) -> Option<()> {
        if self.eatoms.contains_key(&aid) || !self.eatoms.contains_key(&to_aid) {
            return None;
        }
        if self.eatoms.get(&to_aid)?.angle > 0.0 {
            self.add_atom_to_atom_with_ang(aid, to_aid)?;
        } else {
            self.add_atom_to_atom_with_no_ang(aid, to_aid, atoms, bonds, degree)?;
        }
        if let Some(st) = self.eatoms.get_mut(&to_aid) {
            st.pending.retain(|&x| x != aid);
        }
        self.update_new_neighs(aid, atoms, bonds, adjacency, degree, cip_ranks);
        Some(())
    }

    // RDKit✔️✔️: void EmbeddedFrag::expandEfrag(RDKit::INT_LIST &nratms,
    // RDKit✔️✔️:                               std::list<EmbeddedFrag> &efrags) { ... }
    fn expand_efrag(
        &mut self,
        nratms: &mut Vec<usize>,
        efrags: &mut Vec<Self>,
        atoms: &[Atom],
        bonds: &[Bond],
        adjacency: &[Vec<usize>],
        degree: &[usize],
        cip_ranks: &[u32],
    ) -> Option<()> {
        self.merge_frags_with_common(efrags, atoms, bonds, adjacency, degree, cip_ranks)?;
        while !self.attach_pts.is_empty() {
            let aid = self.attach_pts.front().copied()?;
            let nbrs = self.eatoms.get(&aid)?.pending.clone();
            if nbrs.is_empty() {
                return None;
            }
            for nbri in nbrs {
                if let Some(pos) = nratms.iter().position(|&x| x == nbri) {
                    self.add_non_ring_atom(nbri, aid, atoms, bonds, adjacency, degree, cip_ranks)?;
                    nratms.remove(pos);
                } else {
                    let found = efrags
                        .iter()
                        .position(|frag| !frag.done && frag.eatoms.contains_key(&nbri));
                    if let Some(idx) = found {
                        let mut other = efrags.remove(idx);
                        self.merge_no_common(
                            &mut other, aid, nbri, atoms, bonds, adjacency, degree, cip_ranks,
                        )?;
                        let remove_attach = self
                            .eatoms
                            .get(&nbri)
                            .is_some_and(|st| st.pending.is_empty())
                            && self.attach_pts.iter().any(|&x| x == nbri);
                        if remove_attach {
                            self.attach_pts.retain(|&x| x != nbri);
                        }
                    }
                }
            }
            self.attach_pts.pop_front();
            if let Some(st) = self.eatoms.get_mut(&aid) {
                st.pending.clear();
            }
            self.merge_frags_with_common(efrags, atoms, bonds, adjacency, degree, cip_ranks)?;
        }
        Some(())
    }

    // RDKit✔️❌: bool EmbeddedFrag::matchToTemplate(const RDKit::INT_VECT &ringSystemAtoms,
    // RDKit✔️❌:                                    unsigned int ring_count) { ... }
    fn match_to_template(
        &mut self,
        atoms: &[Atom],
        bonds: &[Bond],
        adjacency: &[Vec<usize>],
        degree: &[usize],
        atom_ring_counts: &[usize],
        ring_system_atoms: &[usize],
        ring_count: usize,
    ) -> Option<bool> {
        // RDKit✔️✔️:   if (!coordinate_templates.hasTemplateOfSize(ringSystemAtoms.size())) {
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        if !rdkit_has_template_of_size(ring_system_atoms.len()) {
            return Some(false);
        }

        // RDKit✔️❌:   RDKit::RWMol rs_mol(*dp_mol, true);
        // RDKit✔️❌:   ... mark non-ring-system atoms as dummy, count kept bonds ...
        let rs_mol = build_rdkit_ring_system_molecule(atoms, bonds, ring_system_atoms).ok()?;
        let num_bonds = rs_mol.num_bonds();
        let rs_degree_counts = rdkit_template_degree_counts(&rs_mol, None);

        // RDKit✔️✔️:   for (const auto &mol : coordinate_templates.getMatchingTemplates(ringSystemAtoms.size())) {
        for template in rdkit_matching_templates(ring_system_atoms.len()) {
            if template.graph.bonds.len() != num_bonds {
                continue;
            }
            let template_mol = build_rdkit_template_query_molecule(&template).ok()?;
            if crate::symmetrize_sssr(&template_mol)
                .ok()
                .map(|rings| rings.num_rings())
                .unwrap_or(usize::MAX)
                != ring_count
            {
                continue;
            }
            let template_degree_counts = rdkit_template_degree_counts(&template_mol, None);
            if rs_degree_counts != template_degree_counts {
                continue;
            }
            let Some(match_result) = get_substruct_match(&rs_mol, &template_mol) else {
                continue;
            };
            if !rdkit_check_stereo_chemistry(&rs_mol, &template, &match_result.atom_mapping) {
                continue;
            }
            for (template_aidx, &rs_local_aidx) in match_result.atom_mapping.iter().enumerate() {
                let Some(&rs_aidx) = ring_system_atoms.get(rs_local_aidx) else {
                    return None;
                };
                if let Some(coord) = template
                    .coords_2d
                    .as_ref()
                    .and_then(|coords| coords.get(template_aidx))
                {
                    self.eatoms.insert(
                        rs_aidx,
                        TreeEmbeddedAtom {
                            loc: (coord[0], coord[1]),
                            normal: (0.0, 0.0),
                            ccw: true,
                            cis_trans_nbr: None,
                            angle: -1.0,
                            nbr1: None,
                            nbr2: None,
                            rot_dir: 0,
                            pending: Vec::new(),
                            d_density: -1.0,
                            df_fixed: true,
                        },
                    );
                }
            }
            self.setup_new_neighs(atoms, bonds, adjacency, degree, &[]);
            self.setup_attachment_points(adjacency, atom_ring_counts);
            return Some(true);
        }
        Some(false)
    }

    // RDKit✔️✔️: static void mirrorTransRingAtoms(const RDKit::ROMol &mol,
    // RDKit✔️✔️:                                  const RDKit::INT_VECT &ring,
    // RDKit✔️✔️:                                  RDGeom::INT_POINT2D_MAP &coords) { ... }
    fn mirror_trans_ring_atoms(
        bonds: &[Bond],
        ring: &[usize],
        coords: &mut BTreeMap<usize, (f64, f64)>,
    ) {
        // RDKit✔️✔️:   RDKit::INT_VECT transRingAtoms;
        // RDKit✔️✔️:   for (size_t i = 0; i < ring.size(); ++i) {
        for i in 0..ring.len() {
            // RDKit✔️✔️:     const auto atom1 = ring[i];
            // RDKit✔️✔️:     const auto atom2 = ring[(i + 1) % ring.size()];
            let atom1 = ring[i];
            let atom2 = ring[(i + 1) % ring.len()];
            let Some(bond) = bonds.iter().find(|bond| {
                let begin = bond.begin().index();
                let end = bond.end().index();
                (begin == atom1 && end == atom2) || (begin == atom2 && end == atom1)
            }) else {
                continue;
            };
            // RDKit✔️✔️:     const auto bond = mol.getBondBetweenAtoms(atom1, atom2);
            // RDKit✔️✔️:     if (bond->getBondType() != RDKit::Bond::DOUBLE) {
            // RDKit✔️✔️:       continue;
            // RDKit✔️✔️:     }
            if bond.order() != BondOrder::Double {
                continue;
            }
            // RDKit✔️✔️:     const auto stype = bond->getStereo();
            // RDKit✔️✔️:     if (stype <= RDKit::Bond::STEREOANY) {
            // RDKit✔️✔️:       continue;
            // RDKit✔️✔️:     }
            let stype = bond.stereo();
            if matches!(stype, BondStereo::None | BondStereo::Any) {
                continue;
            }
            // RDKit✔️✔️:     const auto &neighbors = bond->getStereoAtoms();
            // RDKit✔️✔️:     if (neighbors.size() != 2) {
            // RDKit✔️✔️:       continue;
            // RDKit✔️✔️:     }
            let Some(neighbors) = bond.stereo_atoms() else {
                continue;
            };
            // RDKit✔️✔️:     const auto leftIsIn =
            // RDKit✔️✔️:         std::find(ring.begin(), ring.end(), neighbors[0]) != ring.end();
            // RDKit✔️✔️:     const auto rightIsIn =
            // RDKit✔️✔️:         std::find(ring.begin(), ring.end(), neighbors[1]) != ring.end();
            let left_is_in = ring.contains(&neighbors[0].index());
            let right_is_in = ring.contains(&neighbors[1].index());
            // RDKit✔️✔️:     bool isTrans = false;
            // RDKit✔️✔️:     if (stype == RDKit::Bond::STEREOTRANS || stype == RDKit::Bond::STEREOE) {
            // RDKit✔️✔️:       if (leftIsIn == rightIsIn) {
            // RDKit✔️✔️:         isTrans = true;
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:     } else if (leftIsIn != rightIsIn) {
            // RDKit✔️✔️:       isTrans = true;
            // RDKit✔️✔️:     }
            let is_trans = if matches!(stype, BondStereo::Trans | BondStereo::E) {
                left_is_in == right_is_in
            } else {
                left_is_in != right_is_in
            };
            // RDKit✔️✔️:     if (!isTrans) {
            // RDKit✔️✔️:       continue;
            // RDKit✔️✔️:     }
            if !is_trans {
                continue;
            }
            // RDKit✔️✔️:     const auto left = ring[(i + ring.size() - 1) % ring.size()];
            // RDKit✔️✔️:     const auto right = atom2;
            let left = ring[(i + ring.len() - 1) % ring.len()];
            let right = atom2;
            // RDKit✔️✔️:     const auto last = coords[left];
            // RDKit✔️✔️:     const auto ref = coords[right];
            // RDKit✔️✔️:     const auto interest = coords[atom1];
            let last = coords[&left];
            let ref_pt = coords[&right];
            let interest = coords[&atom1];
            // RDKit✔️✔️:     const auto d = last - ref;
            let d = (last.0 - ref_pt.0, last.1 - ref_pt.1);
            let dot = d.0 * d.0 + d.1 * d.1;
            // RDKit✔️✔️:     const double a = (d.x * d.x - d.y * d.y) / d.dotProduct(d);
            // RDKit✔️✔️:     const double b = 2 * d.x * d.y / d.dotProduct(d);
            let a = (d.0 * d.0 - d.1 * d.1) / dot;
            let b = 2.0 * d.0 * d.1 / dot;
            // RDKit✔️✔️:     const double x =
            // RDKit✔️✔️:         a * (interest.x - ref.x) + b * (interest.y - ref.y) + ref.x;
            // RDKit✔️✔️:     const double y =
            // RDKit✔️✔️:         b * (interest.x - ref.x) - a * (interest.y - ref.y) + ref.y;
            let x = a * (interest.0 - ref_pt.0) + b * (interest.1 - ref_pt.1) + ref_pt.0;
            let y = b * (interest.0 - ref_pt.0) - a * (interest.1 - ref_pt.1) + ref_pt.1;
            // RDKit✔️✔️:     coords[atom1] = RDGeom::Point2D(x, y);
            // RDKit✔️✔️:   }
            coords.insert(atom1, (x, y));
        }
        // RDKit✔️✔️: }
    }

    // RDKit✔️✔️: void EmbeddedFrag::initFromRingCoords(const RDKit::INT_VECT &ring,
    // RDKit✔️✔️:                                      const RDGeom::INT_POINT2D_MAP &nringMap) { ... }
    fn init_from_ring_coords(&mut self, ring: &[usize], nring_map: &BTreeMap<usize, (f64, f64)>) {
        // RDKit✔️✔️:   double largestAngle = M_PI * (1 - (2.0 / ring.size()));
        let largest_angle = PI * (1.0 - (2.0 / ring.len() as f64));
        // RDKit✔️✔️:   auto prev = ring.back();
        let mut prev = *ring.last().expect("ring must be non-empty");
        // RDKit✔️✔️:   unsigned int cnt = 0;
        let mut cnt = 0usize;
        // RDKit✔️✔️:   for (auto ai : ring) {
        for &ai in ring {
            // RDKit✔️✔️:     EmbeddedAtom eatm;
            // RDKit✔️✔️:     eatm.loc = nringMap.at(ai);
            // RDKit✔️✔️:     eatm.aid = ai;
            // RDKit✔️✔️:     eatm.angle = largestAngle;
            // RDKit✔️✔️:     eatm.nbr1 = prev;
            let eatm = TreeEmbeddedAtom {
                loc: nring_map[&ai],
                normal: (0.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: largest_angle,
                nbr1: Some(prev),
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            };
            // RDKit✔️✔️:     if (cnt) {
            // RDKit✔️✔️:       d_eatoms[prev].nbr2 = ai;
            // RDKit✔️✔️:     }
            if cnt > 0 {
                if let Some(prev_atom) = self.eatoms.get_mut(&prev) {
                    prev_atom.nbr2 = Some(ai);
                }
            }
            // RDKit✔️✔️:     d_eatoms[ai] = eatm;
            // RDKit✔️✔️:     prev = ai;
            // RDKit✔️✔️:     cnt++;
            // RDKit✔️✔️:   }
            self.eatoms.insert(ai, eatm);
            prev = ai;
            cnt += 1;
        }
        // RDKit✔️✔️:   d_eatoms[prev].nbr2 = ring.front();
        // RDKit✔️✔️: }
        if let Some(prev_atom) = self.eatoms.get_mut(&prev) {
            prev_atom.nbr2 = ring.first().copied();
        }
    }

    // RDKit✔️✔️: void EmbeddedFrag::mergeRing(const EmbeddedFrag &embRing, unsigned int nCommon,
    // RDKit✔️✔️:                             const RDKit::INT_VECT &pinAtoms) { ... }
    fn merge_ring(&mut self, emb_ring: &Self, n_common: usize, pin_atoms: &[usize]) {
        // RDKit✔️✔️:   const auto &oatoms = embRing.GetEmbeddedAtoms();
        // RDKit✔️✔️:   for (const auto &ori : oatoms) {
        for (&aid, ori) in &emb_ring.eatoms {
            // RDKit✔️✔️:     auto aid = ori.first;
            // RDKit✔️✔️:     if (d_eatoms.find(aid) == d_eatoms.end()) {
            // RDKit✔️✔️:       d_eatoms[aid] = ori.second;
            // RDKit✔️✔️:     } else {
            if !self.eatoms.contains_key(&aid) {
                self.eatoms.insert(aid, ori.clone());
            } else {
                // RDKit✔️✔️:       if (nCommon <= 2) {
                // RDKit✔️✔️:         if (std::find(pinAtoms.begin(), pinAtoms.end(), aid) !=
                // RDKit✔️✔️:             pinAtoms.end()) {
                if n_common <= 2 && pin_atoms.contains(&aid) {
                    // RDKit✔️✔️:           d_eatoms[aid].angle += ori.second.angle;
                    let st = self.eatoms.get_mut(&aid).expect("embedded atom exists");
                    st.angle += ori.angle;
                    // RDKit✔️✔️:           if (d_eatoms[aid].nbr1 == ori.second.nbr1) {
                    // RDKit✔️✔️:             d_eatoms[aid].nbr1 = ori.second.nbr2;
                    // RDKit✔️✔️:           } else if (d_eatoms[aid].nbr1 == ori.second.nbr2) {
                    // RDKit✔️✔️:             d_eatoms[aid].nbr1 = ori.second.nbr1;
                    // RDKit✔️✔️:           } else if (d_eatoms[aid].nbr2 == ori.second.nbr1) {
                    // RDKit✔️✔️:             d_eatoms[aid].nbr2 = ori.second.nbr2;
                    // RDKit✔️✔️:           } else if (d_eatoms[aid].nbr2 == ori.second.nbr2) {
                    // RDKit✔️✔️:             d_eatoms[aid].nbr2 = ori.second.nbr1;
                    // RDKit✔️✔️:           }
                    if st.nbr1 == ori.nbr1 {
                        st.nbr1 = ori.nbr2;
                    } else if st.nbr1 == ori.nbr2 {
                        st.nbr1 = ori.nbr1;
                    } else if st.nbr2 == ori.nbr1 {
                        st.nbr2 = ori.nbr2;
                    } else if st.nbr2 == ori.nbr2 {
                        st.nbr2 = ori.nbr1;
                    }
                }
            }
        }
        // RDKit✔️✔️: }
    }

    fn transform(&mut self, trans: [f64; 6]) {
        for st in self.eatoms.values_mut() {
            let loc = st.loc;
            st.loc = transform2d_point(loc, trans);
            let temp = transform2d_point((loc.0 + st.normal.0, loc.1 + st.normal.1), trans);
            st.normal = (temp.0 - st.loc.0, temp.1 - st.loc.1);
        }
    }

    fn translate(&mut self, shift: (f64, f64)) {
        for st in self.eatoms.values_mut() {
            st.loc.0 += shift.0;
            st.loc.1 += shift.1;
        }
    }

    fn canonicalize_orientation(&mut self) {
        // BEGIN RDKIT CPP FUNCTION: EmbeddedFrag::canonicalizeOrientation() (EmbeddedFrag.cpp)
        // RDKit✔️✔️: void EmbeddedFrag::canonicalizeOrientation() {
        // RDKit✔️✔️:   if (d_eatoms.size() <= 1) {
        // RDKit✔️✔️:     return;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   RDGeom::Point2D cent(0.0, 0.0);
        // RDKit✔️✔️:   for (const auto &elem : d_eatoms) {
        // RDKit✔️✔️:     cent += elem.second.loc;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   cent *= (1.0 / d_eatoms.size());
        // RDKit✔️✔️:   double xx = 0.0;
        // RDKit✔️✔️:   double xy = 0.0;
        // RDKit✔️✔️:   double yy = 0.0;
        // RDKit✔️✔️:   for (auto &elem : d_eatoms) {
        // RDKit✔️✔️:     elem.second.loc -= cent;
        // RDKit✔️✔️:     xx += (elem.second.loc.x) * (elem.second.loc.x);
        // RDKit✔️✔️:     xy += (elem.second.loc.x) * (elem.second.loc.y);
        // RDKit✔️✔️:     yy += (elem.second.loc.y) * (elem.second.loc.y);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   RDGeom::Point2D eig1, eig2;
        // RDKit✔️✔️:   auto d = (xx - yy) * (xx - yy) + 4 * xy * xy;
        // RDKit✔️✔️:   d = sqrt(d);
        // RDKit✔️✔️:   eig1.x = 2 * xy;
        // RDKit✔️✔️:   eig1.y = (yy - xx) + d;
        // RDKit✔️✔️:   if (eig1.length() <= 1e-4) {
        // RDKit✔️✔️:     return;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   auto eVal1 = (xx + yy + d) / 2;
        // RDKit✔️✔️:   eig1.normalize();
        // RDKit✔️✔️:   eig2.x = 2 * xy;
        // RDKit✔️✔️:   eig2.y = (yy - xx) - d;
        // RDKit✔️✔️:   auto eVal2 = (xx + yy - d) / 2;
        // RDKit✔️✔️:   if (eig2.length() > 1e-4) {
        // RDKit✔️✔️:     eig2.normalize();
        // RDKit✔️✔️:     if (eVal2 > eVal1) {
        // RDKit✔️✔️:       std::swap(eig1, eig2);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   RDGeom::Transform2D trans;
        // RDKit✔️✔️:   trans.setVal(0, 0, eig1.x);
        // RDKit✔️✔️:   trans.setVal(1, 0, -eig1.y);
        // RDKit✔️✔️:   trans.setVal(0, 1, eig1.y);
        // RDKit✔️✔️:   trans.setVal(1, 1, eig1.x);
        // RDKit✔️✔️:   this->Transform(trans);
        // RDKit✔️✔️: }
        if self.eatoms.len() <= 1 {
            return;
        }
        let mut cent = (0.0f64, 0.0f64);
        for st in self.eatoms.values() {
            cent.0 += st.loc.0;
            cent.1 += st.loc.1;
        }
        let scale = 1.0 / self.eatoms.len() as f64;
        cent.0 *= scale;
        cent.1 *= scale;

        let (mut xx, mut xy, mut yy) = (0.0f64, 0.0f64, 0.0f64);
        for st in self.eatoms.values_mut() {
            st.loc.0 -= cent.0;
            st.loc.1 -= cent.1;
            xx += st.loc.0 * st.loc.0;
            xy += st.loc.0 * st.loc.1;
            yy += st.loc.1 * st.loc.1;
        }

        let d = ((xx - yy) * (xx - yy) + 4.0 * xy * xy).sqrt();
        let mut eig1 = (2.0 * xy, (yy - xx) + d);
        let eig1_len = norm(eig1);
        if eig1_len <= 1.0e-4 {
            return;
        }
        let e_val1 = (xx + yy + d) / 2.0;
        eig1 = (eig1.0 / eig1_len, eig1.1 / eig1_len);

        let mut eig2 = (2.0 * xy, (yy - xx) - d);
        let e_val2 = (xx + yy - d) / 2.0;
        let eig2_len = norm(eig2);
        if eig2_len > 1.0e-4 {
            eig2 = (eig2.0 / eig2_len, eig2.1 / eig2_len);
            if e_val2 > e_val1 {
                std::mem::swap(&mut eig1, &mut eig2);
            }
        }
        if debug_depict_row_active(58) {
            eprintln!(
                "COSMOL_CANON centroid=({:.17},{:.17}) bits=({:#018x},{:#018x}) xx={:.17} xy={:.17} yy={:.17} d={:.17} eig1=({:.17},{:.17}) eig2=({:.17},{:.17}) e1={:.17} e2={:.17}",
                cent.0,
                cent.1,
                cent.0.to_bits(),
                cent.1.to_bits(),
                xx,
                xy,
                yy,
                d,
                eig1.0,
                eig1.1,
                eig2.0,
                eig2.1,
                e_val1,
                e_val2
            );
            for (aid, st) in self.eatoms.iter().take(6) {
                eprintln!(
                    "COSMOL_CANON_PRE atom={} loc=({:.17},{:.17}) bits=({:#018x},{:#018x})",
                    aid,
                    st.loc.0,
                    st.loc.1,
                    st.loc.0.to_bits(),
                    st.loc.1.to_bits()
                );
            }
        }
        let trans = [eig1.0, eig1.1, 0.0, -eig1.1, eig1.0, 0.0];
        self.transform(trans);
        if debug_depict_row_active(58) {
            eprintln!(
                "COSMOL_CANON_TRANS data=[{:.17},{:.17},{:.17},{:.17},{:.17},{:.17}] bits=[{:#018x},{:#018x},{:#018x},{:#018x},{:#018x},{:#018x}]",
                trans[0],
                trans[1],
                trans[2],
                trans[3],
                trans[4],
                trans[5],
                trans[0].to_bits(),
                trans[1].to_bits(),
                trans[2].to_bits(),
                trans[3].to_bits(),
                trans[4].to_bits(),
                trans[5].to_bits()
            );
            for (aid, st) in self.eatoms.iter().take(6) {
                eprintln!(
                    "COSMOL_CANON_POST atom={} loc=({:.17},{:.17}) bits=({:#018x},{:#018x})",
                    aid,
                    st.loc.0,
                    st.loc.1,
                    st.loc.0.to_bits(),
                    st.loc.1.to_bits()
                );
            }
        }
    }

    fn compute_box(&self) -> Option<(f64, f64, f64, f64)> {
        let mut px = -1.0e8f64;
        let mut nx = 1.0e8f64;
        let mut py = -1.0e8f64;
        let mut ny = 1.0e8f64;
        for st in self.eatoms.values() {
            px = px.max(st.loc.0);
            nx = nx.min(st.loc.0);
            py = py.max(st.loc.1);
            ny = ny.min(st.loc.1);
        }
        if self.eatoms.is_empty() {
            None
        } else {
            Some((px, -nx, py, -ny))
        }
    }

    fn reflect(&mut self, loc1: (f64, f64), loc2: (f64, f64)) {
        for st in self.eatoms.values_mut() {
            let temp = (st.loc.0 + st.normal.0, st.loc.1 + st.normal.1);
            st.loc = rdkit_reflect_point(st.loc, loc1, loc2);
            let temp = rdkit_reflect_point(temp, loc1, loc2);
            st.normal = (temp.0 - st.loc.0, temp.1 - st.loc.1);
            st.ccw = !st.ccw;
        }
    }

    // RDKit✔️❌: RDGeom::Transform2D EmbeddedFrag::computeOneAtomTrans(
    // RDKit✔️❌:     unsigned int commAid, const EmbeddedFrag &other) { ... }
    fn compute_one_atom_trans(&self, comm_aid: usize, other: &Self) -> Option<[f64; 6]> {
        let rcr = self.eatoms.get(&comm_aid)?.loc;
        let oeatm = other.eatoms.get(&comm_aid)?;
        let ccr = oeatm.loc;
        let onb1 = oeatm.nbr1?;
        let onb2 = oeatm.nbr2?;
        let onb1_loc = other.eatoms.get(&onb1)?.loc;
        let onb2_loc = other.eatoms.get(&onb2)?.loc;
        let mid_pt = (
            (onb1_loc.0 + onb2_loc.0) * 0.5,
            (onb1_loc.1 + onb2_loc.1) * 0.5,
        );
        let nb1 = self.eatoms.get(&comm_aid)?.nbr1?;
        let nb2 = self.eatoms.get(&comm_aid)?.nbr2?;
        let nbp1 = self.eatoms.get(&nb1)?.loc;
        let nbp2 = self.eatoms.get(&nb2)?.loc;
        let ang = self.eatoms.get(&comm_aid)?.angle;
        let largest_angle = 2.0 * PI - ang;
        let bpt = rdkit_compute_bisect_point(rcr, largest_angle, nbp1, nbp2);
        Some(transform2d_set_transform_two_point(rcr, bpt, ccr, mid_pt))
    }

    // RDKit✔️✔️: RDGeom::Transform2D EmbeddedFrag::computeTwoAtomTrans(
    // RDKit✔️✔️:     unsigned int aid1, unsigned int aid2,
    // RDKit✔️✔️:     const RDGeom::INT_POINT2D_MAP &nringCor) { ... }
    fn compute_two_atom_trans(
        &self,
        aid1: usize,
        aid2: usize,
        nring_cor: &BTreeMap<usize, (f64, f64)>,
    ) -> Option<[f64; 6]> {
        let loc1 = *nring_cor.get(&aid1)?;
        let loc2 = *nring_cor.get(&aid2)?;
        let ref1 = self.eatoms.get(&aid1)?.loc;
        let ref2 = self.eatoms.get(&aid2)?.loc;
        Some(transform2d_set_transform_two_point(ref1, ref2, loc1, loc2))
    }

    // RDKit✔️❌: void EmbeddedFrag::reflectIfNecessaryDensity(EmbeddedFrag &embFrag,
    // RDKit✔️❌:                                             unsigned int aid1,
    // RDKit✔️❌:                                             unsigned int aid2) { ... }
    fn reflect_if_necessary_density(&self, emb_frag: &mut Self, aid1: usize, aid2: usize) {
        let pin1 = self.eatoms[&aid1].loc;
        let pin2 = self.eatoms[&aid2].loc;
        let mut density_normal = 0.0f64;
        let mut density_reflect = 0.0f64;
        for (&oa, ost) in &emb_frag.eatoms {
            if self.eatoms.contains_key(&oa) {
                continue;
            }
            let loc = ost.loc;
            let rloc = rdkit_reflect_point(loc, pin1, pin2);
            for tst in self.eatoms.values() {
                let d = norm((tst.loc.0 - loc.0, tst.loc.1 - loc.1));
                let rd = norm((tst.loc.0 - rloc.0, tst.loc.1 - rloc.1));
                density_normal += if d > 1.0e-3 { 1.0 / d } else { 1000.0 };
                density_reflect += if rd > 1.0e-3 { 1.0 / rd } else { 1000.0 };
            }
        }
        if density_normal - density_reflect > 1.0e-4 {
            emb_frag.reflect(pin1, pin2);
        }
    }

    // RDKit✔️✔️: void EmbeddedFrag::computeDistMat(DOUBLE_SMART_PTR &dmat) { ... }
    fn compute_dist_mat(&self, num_atoms: usize) -> Vec<f64> {
        let dsize = num_atoms.saturating_mul(num_atoms.saturating_sub(1)) / 2;
        let mut dmat = vec![-1.0; dsize];
        for (&ai0, sti) in &self.eatoms {
            for (&aj0, stj) in self.eatoms.range(..ai0) {
                let (mut ai, mut aj) = (ai0, aj0);
                if ai < aj {
                    std::mem::swap(&mut ai, &mut aj);
                }
                let idx = (ai * (ai - 1) / 2) + aj;
                dmat[idx] = norm((stj.loc.0 - sti.loc.0, stj.loc.1 - sti.loc.1));
            }
        }
        dmat
    }

    // RDKit✔️✔️: double EmbeddedFrag::mimicDistMatAndDensityCostFunc(
    // RDKit✔️✔️:     const DOUBLE_SMART_PTR *dmat, double mimicDmatWt) { ... }
    fn mimic_dist_mat_and_density_cost_func(
        &self,
        num_atoms: usize,
        dmat: Option<&[f64]>,
        mimic_dmat_wt: f64,
    ) -> f64 {
        if num_atoms < 2 {
            return 0.0;
        }
        let dsize = num_atoms * (num_atoms - 1) / 2;
        let dmat_2d = self.compute_dist_mat(num_atoms);
        let mut res1 = 0.0;
        let mut res2 = 0.0;
        for i in 0..dsize {
            let d = dmat_2d[i];
            let d2 = d * d;
            if d2 > 1.0e-3 {
                res1 += 1.0 / d2;
            } else {
                res1 += 1000.0;
            }
            if let Some(src) = dmat {
                if src.get(i).copied().unwrap_or(-1.0) >= 0.0 {
                    let dd = d - src[i];
                    res2 += dd * dd;
                }
            }
        }
        let wt = mimic_dmat_wt.clamp(0.0, 1.0);
        ((1.0 - wt) * res1) + (wt * res2)
    }

    // RDKit✔️✔️: void EmbeddedFrag::permuteBonds(unsigned int aid, unsigned int aid1,
    // RDKit✔️✔️:                                 unsigned int aid2) { ... }
    fn permute_bonds(&mut self, aid: usize, aid1: usize, aid2: usize, adjacency: &[Vec<usize>]) {
        let rl1 = self.eatoms[&aid].loc;
        let rl2 = (
            (self.eatoms[&aid1].loc.0 + self.eatoms[&aid2].loc.0) * 0.5,
            (self.eatoms[&aid1].loc.1 + self.eatoms[&aid2].loc.1) * 0.5,
        );
        let mut frag_a = Vec::new();
        let mut frag_b = Vec::new();
        recurse_atom_one_side(aid1, aid, adjacency, &mut frag_a);
        recurse_atom_one_side(aid2, aid, adjacency, &mut frag_b);
        for fi in frag_a {
            if let Some(st) = self.eatoms.get_mut(&fi) {
                let temp = (st.loc.0 + st.normal.0, st.loc.1 + st.normal.1);
                st.loc = rdkit_reflect_point(st.loc, rl1, rl2);
                let temp = rdkit_reflect_point(temp, rl1, rl2);
                st.normal = (temp.0 - st.loc.0, temp.1 - st.loc.1);
                st.ccw = !st.ccw;
            }
        }
        for fi in frag_b {
            if let Some(st) = self.eatoms.get_mut(&fi) {
                let temp = (st.loc.0 + st.normal.0, st.loc.1 + st.normal.1);
                st.loc = rdkit_reflect_point(st.loc, rl1, rl2);
                let temp = rdkit_reflect_point(temp, rl1, rl2);
                st.normal = (temp.0 - st.loc.0, temp.1 - st.loc.1);
                st.ccw = !st.ccw;
            }
        }
    }

    // RDKit✔️✔️: void EmbeddedFrag::randomSampleFlipsAndPermutations(
    // RDKit✔️✔️:     unsigned int nBondsPerSample, unsigned int nSamples, int seed,
    // RDKit✔️✔️:     const DOUBLE_SMART_PTR *dmat, double mimicDmatWt, bool permuteDeg4Nodes) { ... }
    fn random_sample_flips_and_permutations(
        &mut self,
        atoms: &[Atom],
        bonds: &[Bond],
        comp: &[usize],
        adjacency: &[Vec<usize>],
        n_bonds_per_sample: usize,
        n_samples: usize,
        seed: i32,
        dmat: Option<&[f64]>,
        mimic_dmat_wt: f64,
        permute_deg4_nodes: bool,
    ) {
        let rot_bonds = rdkit_get_all_rotatable_bonds(bonds, |bid| {
            is_ring_bond_in_component(bid, bonds, comp, adjacency)
        });
        let nb = rot_bonds.len();

        let mut nd4 = 0usize;
        let mut deg4nodes = Vec::<usize>::new();
        let mut deg4_nbr_bids = Vec::<Vec<usize>>::new();
        let mut deg4_nbr_aids = Vec::<Vec<usize>>::new();
        if permute_deg4_nodes {
            for caid in 0..atoms.len() {
                if adjacency[caid].len() == 4
                    && !adjacency[caid].iter().any(|&nbr| {
                        bond_between_idx(bonds, caid, nbr).is_some_and(|bid| {
                            is_ring_bond_in_component(bid, bonds, comp, adjacency)
                        })
                    })
                {
                    let (aids, bids) = rdkit_get_nbr_atom_and_bond_ids(caid, bonds, adjacency);
                    let all_in = aids
                        .iter()
                        .all(|aid| self.eatoms.get(aid).is_some_and(|st| !st.df_fixed));
                    if all_in {
                        deg4nodes.push(caid);
                        deg4_nbr_bids.push(bids);
                        deg4_nbr_aids.push(aids);
                    }
                }
            }
            nd4 = deg4nodes.len();
        }

        let nt = nb + nd4;
        if nt == 0 {
            return;
        }
        let n_per_sample = nt.min(n_bonds_per_sample);
        let mut rng = RdkitMinStdRand::new(seed);

        let mut best_crd_map = BTreeMap::<usize, (f64, f64)>::new();
        let mut best_dens =
            self.mimic_dist_mat_and_density_cost_func(atoms.len(), dmat, mimic_dmat_wt);
        for (&aid, st) in &self.eatoms {
            best_crd_map.insert(aid, st.loc);
        }

        for _ in 0..n_samples {
            for _ in 0..n_per_sample {
                let ri = rng.next_bounded(nt - 1);
                if ri < nb {
                    self.flip_about_bond(
                        rot_bonds[ri],
                        bonds,
                        adjacency,
                        |bid| is_ring_bond_in_component(bid, bonds, comp, adjacency),
                        true,
                    );
                } else {
                    let d4i = ri - nb;
                    let ai = deg4nodes[d4i];
                    let nbr_locs: Vec<(f64, f64)> = deg4_nbr_aids[d4i]
                        .iter()
                        .map(|aid| self.eatoms[aid].loc)
                        .collect();
                    let bnd_pairs = rdkit_find_bond_pairs_to_permute_deg4(
                        self.eatoms[&ai].loc,
                        &deg4_nbr_bids[d4i],
                        &nbr_locs,
                    );
                    let fbi = if rng.next_unit_f64() > 0.5 { 1 } else { 0 };
                    let (bid1, bid2) = bnd_pairs[fbi];
                    let aid1 = if bonds[bid1].begin().index() == ai {
                        bonds[bid1].end().index()
                    } else {
                        bonds[bid1].begin().index()
                    };
                    let aid2 = if bonds[bid2].begin().index() == ai {
                        bonds[bid2].end().index()
                    } else {
                        bonds[bid2].begin().index()
                    };
                    self.permute_bonds(ai, aid1, aid2, adjacency);
                }
            }
            let density =
                self.mimic_dist_mat_and_density_cost_func(atoms.len(), dmat, mimic_dmat_wt);
            if best_dens - density > 1.0e-4 {
                best_dens = density;
                for (&aid, st) in &self.eatoms {
                    best_crd_map.insert(aid, st.loc);
                }
            }
        }
        for (&aid, loc) in &best_crd_map {
            self.eatoms.get_mut(&aid).expect("embedded atom").loc = *loc;
        }
    }

    // RDKit✔️✔️: std::vector<PAIR_I_I> EmbeddedFrag::findCollisions(const double *dmat,
    // RDKit✔️✔️:                                                    bool includeBonds) { ... }
    fn find_collisions(
        &mut self,
        atoms: &[Atom],
        bonds: &[Bond],
        adjacency: &[Vec<usize>],
        dmat: &[f64],
        include_bonds: bool,
    ) -> Vec<(usize, usize)> {
        let mut res = Vec::new();
        for st in self.eatoms.values_mut() {
            st.d_density = 0.0;
        }
        let col_thres2 = RDKIT_COLLISION_THRES * RDKIT_COLLISION_THRES;
        let atom_ids: Vec<usize> = self.eatoms.keys().copied().collect();
        for &ai in &atom_ids {
            let atom_type_factor1 = if atoms[ai].atomic_number() != 6 {
                RDKIT_HETEROATOM_COLL_SCALE
            } else {
                1.0
            };
            for &aj in atom_ids.iter().filter(|&&aj| aj < ai) {
                let atom_type_factor2 = if atoms[aj].atomic_number() != 6 {
                    RDKIT_HETEROATOM_COLL_SCALE
                } else {
                    1.0
                };
                let pti = self.eatoms[&ai].loc;
                let ptj = self.eatoms[&aj].loc;
                let dx = ptj.0 - pti.0;
                let dy = ptj.1 - pti.1;
                let mut d2 = dx * dx + dy * dy;
                let add = if d2 > 1.0e-3 { 1.0 / d2 } else { 1000.0 };
                if let Some(st) = self.eatoms.get_mut(&ai) {
                    st.d_density += add;
                }
                if let Some(st) = self.eatoms.get_mut(&aj) {
                    st.d_density += add;
                }
                d2 /= atom_type_factor1 * atom_type_factor2;
                if d2 < col_thres2 {
                    res.push((ai, aj));
                }
            }
        }
        if include_bonds {
            let bond_thres2 = RDKIT_BOND_THRES * RDKIT_BOND_THRES;
            for b1 in bonds {
                let bid1 = b1.id().index();
                let beg1 = b1.begin().index();
                let end1 = b1.end().index();
                if !(self.eatoms.contains_key(&beg1) && self.eatoms.contains_key(&end1)) {
                    continue;
                }
                let pbeg1 = self.eatoms[&beg1].loc;
                let pend1 = self.eatoms[&end1].loc;
                let v1 = (pend1.0 - pbeg1.0, pend1.1 - pbeg1.1);
                let avg1 = ((pend1.0 + pbeg1.0) * 0.5, (pend1.1 + pbeg1.1) * 0.5);
                for b2 in bonds.iter().filter(|b2| b2.id().index() > bid1) {
                    let beg2 = b2.begin().index();
                    let end2 = b2.end().index();
                    if !(self.eatoms.contains_key(&beg2) && self.eatoms.contains_key(&end2)) {
                        continue;
                    }
                    let pbeg2 = self.eatoms[&beg2].loc;
                    let pend2 = self.eatoms[&end2].loc;
                    let avg2 = (
                        ((pend2.0 + pbeg2.0) * 0.5) - avg1.0,
                        ((pend2.1 + pbeg2.1) * 0.5) - avg1.1,
                    );
                    let avg2_len_sq = avg2.0 * avg2.0 + avg2.1 * avg2.1;
                    if avg2_len_sq < 0.5 && avg2_len_sq < bond_thres2 {
                        let v2 = (pbeg2.0 - pbeg1.0, pbeg2.1 - pbeg1.1);
                        let v3 = (pend2.0 - pbeg1.0, pend2.1 - pbeg1.1);
                        let val_prod = cross_val(v1, v2) * cross_val(v1, v3);
                        if val_prod < -1.0e-6 {
                            res.push(find_closest_pair(beg1, end1, beg2, end2, atoms.len(), dmat));
                        }
                    }
                }
            }
        }
        let _ = adjacency;
        res
    }

    // RDKit✔️✔️: double EmbeddedFrag::totalDensity() { ... }
    fn total_density(&self) -> f64 {
        self.eatoms.values().map(|st| st.d_density).sum()
    }

    // RDKit✔️✔️: void EmbeddedFrag::flipAboutBond(unsigned int bondId, bool flipEnd) { ... }
    fn flip_about_bond(
        &mut self,
        bond_id: usize,
        bonds: &[Bond],
        adjacency: &[Vec<usize>],
        is_ring_bond: impl Fn(usize) -> bool,
        flip_end: bool,
    ) {
        if is_ring_bond(bond_id) {
            return;
        }
        let bond = &bonds[bond_id];
        let mut beg_aid = bond.begin().index();
        let mut end_aid = bond.end().index();
        if !flip_end {
            std::mem::swap(&mut beg_aid, &mut end_aid);
        }
        let mut beg_loc = self.eatoms[&beg_aid].loc;
        let mut end_loc = self.eatoms[&end_aid].loc;
        let mut end_side_aids = Vec::new();
        recurse_atom_one_side(end_aid, beg_aid, adjacency, &mut end_side_aids);
        let n_atoms_fixed = self.eatoms.values().filter(|st| st.df_fixed).count();
        let n_end_atoms_fixed = if n_atoms_fixed > 0 {
            end_side_aids
                .iter()
                .filter(|aid| self.eatoms.get(aid).is_some_and(|st| st.df_fixed))
                .count()
        } else {
            0
        };
        let mut end_side_flip = true;
        if n_end_atoms_fixed > 0 {
            return;
        } else {
            let nats = self.eatoms.len();
            let n_end_side = end_side_aids.len();
            if (nats - n_end_side) < n_end_side {
                end_side_flip = false;
            }
        }
        for (&aid, st) in &mut self.eatoms {
            let on_end_side = end_side_aids.contains(&aid);
            if end_side_flip ^ !on_end_side {
                let temp = (st.loc.0 + st.normal.0, st.loc.1 + st.normal.1);
                st.loc = rdkit_reflect_point(st.loc, beg_loc, end_loc);
                let temp = rdkit_reflect_point(temp, beg_loc, end_loc);
                st.normal = (temp.0 - st.loc.0, temp.1 - st.loc.1);
                st.ccw = !st.ccw;
                if aid == beg_aid {
                    beg_loc = st.loc;
                }
                if aid == end_aid {
                    end_loc = st.loc;
                }
            }
        }
    }

    // RDKit✔️✔️: void EmbeddedFrag::openAngles(const double *dmat, unsigned int aid1,
    // RDKit✔️✔️:                               unsigned int aid2) { ... }
    fn open_angles(
        &mut self,
        adjacency: &[Vec<usize>],
        num_atoms: usize,
        dmat: &[f64],
        aid1: usize,
        aid2: usize,
    ) {
        let deg1 = adjacency[aid1].len();
        let deg2 = adjacency[aid2].len();
        let fixed1 = self.eatoms[&aid1].df_fixed;
        let fixed2 = self.eatoms[&aid2].df_fixed;
        if (deg1 > 1 || fixed1) && (deg2 > 1 || fixed2) {
            return;
        }
        let (aid_a, aid_b, kind) = if (deg1 == 1 && !fixed1) && (deg2 == 1 && !fixed2) {
            (
                find_deg1_neighbor(adjacency, aid1),
                find_deg1_neighbor(adjacency, aid2),
                1,
            )
        } else if (deg1 == 1 && !fixed1) && (deg2 > 1 || fixed2) {
            let aid_a = find_deg1_neighbor(adjacency, aid1);
            (
                aid_a,
                aid_a.map(|a| find_closest_neighbor(adjacency, dmat, num_atoms, a, aid2)),
                2,
            )
        } else {
            let aid_b = find_deg1_neighbor(adjacency, aid2);
            (
                aid_b.map(|b| find_closest_neighbor(adjacency, dmat, num_atoms, b, aid1)),
                aid_b,
                3,
            )
        };
        let (Some(aid_a), Some(aid_b)) = (aid_a, aid_b) else {
            return;
        };
        let v2 = (
            self.eatoms[&aid1].loc.0 - self.eatoms[&aid_a].loc.0,
            self.eatoms[&aid1].loc.1 - self.eatoms[&aid_a].loc.1,
        );
        let v1 = (
            self.eatoms[&aid_b].loc.0 - self.eatoms[&aid_a].loc.0,
            self.eatoms[&aid_b].loc.1 - self.eatoms[&aid_a].loc.1,
        );
        let cross = cross_val(v1, v2);
        match kind {
            1 => {
                let mut angle = RDKIT_ANGLE_OPEN;
                if cross < 0.0 {
                    angle *= -1.0;
                }
                let p1 = rotate_around(self.eatoms[&aid1].loc, self.eatoms[&aid_a].loc, angle);
                let p2 = rotate_around(self.eatoms[&aid2].loc, self.eatoms[&aid_b].loc, -angle);
                self.eatoms.get_mut(&aid1).expect("embedded atom").loc = p1;
                self.eatoms.get_mut(&aid2).expect("embedded atom").loc = p2;
            }
            2 => {
                let mut angle = 2.0 * RDKIT_ANGLE_OPEN;
                if cross < 0.0 {
                    angle *= -1.0;
                }
                let p1 = rotate_around(self.eatoms[&aid1].loc, self.eatoms[&aid_a].loc, angle);
                self.eatoms.get_mut(&aid1).expect("embedded atom").loc = p1;
            }
            3 => {
                let mut angle = -2.0 * RDKIT_ANGLE_OPEN;
                if cross < 0.0 {
                    angle *= -1.0;
                }
                let p2 = rotate_around(self.eatoms[&aid2].loc, self.eatoms[&aid_b].loc, angle);
                self.eatoms.get_mut(&aid2).expect("embedded atom").loc = p2;
            }
            _ => {}
        }
    }

    // RDKit✔️✔️: void EmbeddedFrag::removeCollisionsBondFlip() { ... }
    fn remove_collisions_bond_flip(
        &mut self,
        atoms: &[Atom],
        bonds: &[Bond],
        comp: &[usize],
        adjacency: &[Vec<usize>],
    ) {
        let dmat = component_graph_distance_matrix(atoms.len(), comp, adjacency);
        let mut colls = self.find_collisions(atoms, bonds, adjacency, &dmat, true);
        let mut done_bonds = BTreeMap::<usize, usize>::new();
        let mut iter = 0usize;
        while iter < RDKIT_MAX_COLL_ITERS && !colls.is_empty() {
            let ncols = colls.len();
            let c_aids = colls[0];
            let rot_bonds = rdkit_get_rotatable_bonds_between(
                c_aids.0,
                c_aids.1,
                atoms.len(),
                bonds,
                adjacency,
                |bid| is_ring_bond_in_component(bid, bonds, comp, adjacency),
            );
            let prev_density = self.total_density();
            for ri in rot_bonds {
                let done = *done_bonds.get(&ri).unwrap_or(&0);
                if done >= RDKIT_NUM_BONDS_FLIPS {
                    continue;
                }
                done_bonds.insert(ri, done + 1);
                self.flip_about_bond(
                    ri,
                    bonds,
                    adjacency,
                    |bid| is_ring_bond_in_component(bid, bonds, comp, adjacency),
                    true,
                );
                colls = self.find_collisions(atoms, bonds, adjacency, &dmat, true);
                let new_density = self.total_density();
                if colls.len() < ncols {
                    done_bonds.insert(ri, RDKIT_NUM_BONDS_FLIPS);
                    break;
                } else if colls.len() == ncols && new_density < prev_density {
                    break;
                } else {
                    self.flip_about_bond(
                        ri,
                        bonds,
                        adjacency,
                        |bid| is_ring_bond_in_component(bid, bonds, comp, adjacency),
                        true,
                    );
                    colls = self.find_collisions(atoms, bonds, adjacency, &dmat, true);
                    self.flip_about_bond(
                        ri,
                        bonds,
                        adjacency,
                        |bid| is_ring_bond_in_component(bid, bonds, comp, adjacency),
                        false,
                    );
                    colls = self.find_collisions(atoms, bonds, adjacency, &dmat, true);
                    let new_density = self.total_density();
                    if colls.len() < ncols {
                        done_bonds.insert(ri, RDKIT_NUM_BONDS_FLIPS);
                        break;
                    } else if colls.len() == ncols && new_density < prev_density {
                        break;
                    } else {
                        self.flip_about_bond(
                            ri,
                            bonds,
                            adjacency,
                            |bid| is_ring_bond_in_component(bid, bonds, comp, adjacency),
                            false,
                        );
                        colls = self.find_collisions(atoms, bonds, adjacency, &dmat, true);
                    }
                }
            }
            iter += 1;
        }
    }

    // RDKit✔️✔️: void EmbeddedFrag::removeCollisionsOpenAngles() { ... }
    fn remove_collisions_open_angles(
        &mut self,
        atoms: &[Atom],
        bonds: &[Bond],
        comp: &[usize],
        adjacency: &[Vec<usize>],
    ) {
        let dmat = component_graph_distance_matrix(atoms.len(), comp, adjacency);
        for (aid1, aid2) in self.find_collisions(
            atoms,
            bonds,
            adjacency,
            &self.compute_dist_mat(atoms.len()),
            false,
        ) {
            self.open_angles(adjacency, atoms.len(), &dmat, aid1, aid2);
        }
    }

    // RDKit✔️✔️: void EmbeddedFrag::removeCollisionsShortenBonds() { ... }
    fn remove_collisions_shorten_bonds(
        &mut self,
        atoms: &[Atom],
        bonds: &[Bond],
        comp: &[usize],
        adjacency: &[Vec<usize>],
    ) {
        let coord_dmat = self.compute_dist_mat(atoms.len());
        let mut colls = self.find_collisions(atoms, bonds, adjacency, &coord_dmat, false);
        let mut ncols = colls.len();
        let mut iter = 0usize;
        while ncols > 0 && iter < RDKIT_MAX_COLL_ITERS {
            let (mut aid1, mut aid2) = colls[0];
            let mut fixed1 = self.eatoms[&aid1].df_fixed;
            let mut fixed2 = self.eatoms[&aid2].df_fixed;
            if fixed1 && fixed2 {
                colls.remove(0);
                ncols = colls.len();
                iter += 1;
                continue;
            }
            let mut deg1 = adjacency[aid1].len();
            let mut deg2 = adjacency[aid2].len();
            if fixed1 || (deg2 > deg1 && !fixed2) {
                std::mem::swap(&mut deg1, &mut deg2);
                std::mem::swap(&mut aid1, &mut aid2);
                std::mem::swap(&mut fixed1, &mut fixed2);
            }
            if let Some(mut path) = shortest_path_in_component(aid1, aid2, comp, adjacency) {
                if path.first().copied() == Some(aid1) {
                    path.remove(0);
                }
                let n_open = any_non_ring_bonds_on_path(aid1, &path, bonds, comp, adjacency);
                if n_open > 0 {
                    if deg1 == 1
                        && let Some(aid_a) = find_deg1_neighbor(adjacency, aid1)
                    {
                        let mut loc = (
                            self.eatoms[&aid1].loc.0 - self.eatoms[&aid_a].loc.0,
                            self.eatoms[&aid1].loc.1 - self.eatoms[&aid_a].loc.1,
                        );
                        loc.0 *= 0.9;
                        loc.1 *= 0.9;
                        if norm(loc) > 0.75 {
                            self.eatoms.get_mut(&aid1).expect("embedded atom").loc = (
                                self.eatoms[&aid_a].loc.0 + loc.0,
                                self.eatoms[&aid_a].loc.1 + loc.1,
                            );
                        }
                    }
                    if deg2 == 1
                        && !fixed2
                        && let Some(aid_a) = find_deg1_neighbor(adjacency, aid2)
                    {
                        let mut loc = (
                            self.eatoms[&aid2].loc.0 - self.eatoms[&aid_a].loc.0,
                            self.eatoms[&aid2].loc.1 - self.eatoms[&aid_a].loc.1,
                        );
                        loc.0 *= 0.9;
                        loc.1 *= 0.9;
                        if norm(loc) > 0.75 {
                            self.eatoms.get_mut(&aid2).expect("embedded atom").loc = (
                                self.eatoms[&aid_a].loc.0 + loc.0,
                                self.eatoms[&aid_a].loc.1 + loc.1,
                            );
                        }
                    }
                } else {
                    let mut r_path = Vec::<usize>::new();
                    let mut nbr_map = BTreeMap::<usize, Vec<usize>>::new();
                    recurse_deg_two_ring_atoms_component(
                        aid1,
                        bonds,
                        comp,
                        adjacency,
                        &mut r_path,
                        &mut nbr_map,
                    );
                    if r_path.is_empty() {
                        recurse_deg_two_ring_atoms_component(
                            aid2,
                            bonds,
                            comp,
                            adjacency,
                            &mut r_path,
                            &mut nbr_map,
                        );
                    }
                    let mut move_map = BTreeMap::<usize, (f64, f64)>::new();
                    for &rpi in &r_path {
                        if self.eatoms[&rpi].df_fixed {
                            continue;
                        }
                        let nbrs = &nbr_map[&rpi];
                        let mut mv = (
                            (self.eatoms[&nbrs[0]].loc.0 + self.eatoms[&nbrs[1]].loc.0) * 0.5
                                - self.eatoms[&rpi].loc.0,
                            (self.eatoms[&nbrs[0]].loc.1 + self.eatoms[&nbrs[1]].loc.1) * 0.5
                                - self.eatoms[&rpi].loc.1,
                        );
                        let len = norm(mv);
                        if len > 1.0e-12 {
                            mv.0 = (mv.0 / len) * RDKIT_COLLISION_THRES;
                            mv.1 = (mv.1 / len) * RDKIT_COLLISION_THRES;
                            move_map.insert(rpi, mv);
                        }
                    }
                    for rpi in r_path {
                        if let Some(mv) = move_map.get(&rpi) {
                            let cur = self.eatoms[&rpi].loc;
                            self.eatoms.get_mut(&rpi).expect("embedded atom").loc =
                                (cur.0 + mv.0, cur.1 + mv.1);
                        }
                    }
                }
                colls = self.find_collisions(atoms, bonds, adjacency, &coord_dmat, false);
            } else {
                colls.remove(0);
            }
            ncols = colls.len();
            iter += 1;
        }
    }

    fn merge_no_common(
        &mut self,
        emb_obj: &mut Self,
        to_aid: usize,
        nbr_aid: usize,
        atoms: &[Atom],
        bonds: &[Bond],
        adjacency: &[Vec<usize>],
        degree: &[usize],
        cip_ranks: &[u32],
    ) -> Option<()> {
        self.add_non_ring_atom(nbr_aid, to_aid, atoms, bonds, adjacency, degree, cip_ranks)?;
        emb_obj.add_non_ring_atom(to_aid, nbr_aid, atoms, bonds, adjacency, degree, cip_ranks)?;
        self.merge_with_common(
            emb_obj,
            vec![to_aid, nbr_aid],
            atoms,
            bonds,
            adjacency,
            degree,
            cip_ranks,
        )
    }

    fn find_common_atoms(&self, other: &Self) -> Vec<usize> {
        let mut common = Vec::new();
        for aid in self.eatoms.keys() {
            if other.eatoms.contains_key(aid) {
                common.push(*aid);
            }
        }
        common
    }

    fn merge_with_common(
        &mut self,
        other: &mut Self,
        mut common: Vec<usize>,
        atoms: &[Atom],
        bonds: &[Bond],
        adjacency: &[Vec<usize>],
        degree: &[usize],
        cip_ranks: &[u32],
    ) -> Option<()> {
        let mut ct_case = 0u8;
        if common.len() == 1 {
            let comm_aid = common[0];
            let other_atom;
            if self
                .eatoms
                .get(&comm_aid)
                .and_then(|st| st.cis_trans_nbr)
                .is_some()
            {
                ct_case = 2;
                other_atom = self.eatoms.get(&comm_aid)?.nbr1;
                if let Some(aid) = other_atom {
                    other.add_non_ring_atom(
                        aid, comm_aid, atoms, bonds, adjacency, degree, cip_ranks,
                    )?;
                }
            } else if other
                .eatoms
                .get(&comm_aid)
                .and_then(|st| st.cis_trans_nbr)
                .is_some()
            {
                ct_case = 1;
                other_atom = other.eatoms.get(&comm_aid)?.nbr1;
                if let Some(aid) = other_atom {
                    self.add_non_ring_atom(
                        aid, comm_aid, atoms, bonds, adjacency, degree, cip_ranks,
                    )?;
                }
            } else {
                other_atom = self.eatoms.get(&comm_aid)?.nbr1;
                if let Some(aid) = other_atom {
                    other.add_non_ring_atom(
                        aid, comm_aid, atoms, bonds, adjacency, degree, cip_ranks,
                    )?;
                }
            }
            if let Some(aid) = other_atom {
                common.push(aid);
            }
        }

        let trans = if common.len() == 1 {
            self.compute_one_atom_trans(common[0], other)?
        } else {
            let aid1 = common[0];
            let aid2 = common[1];
            let ref1 = self.eatoms.get(&aid1)?.loc;
            let ref2 = self.eatoms.get(&aid2)?.loc;
            let oth1 = other.eatoms.get(&aid1)?.loc;
            let oth2 = other.eatoms.get(&aid2)?.loc;
            transform2d_set_transform_two_point(ref1, ref2, oth1, oth2)
        };
        other.transform(trans);

        if common.len() >= 2 {
            let aid1 = common[0];
            let aid2 = common[1];
            if ct_case > 0 {
                let p1_loc = self.eatoms.get(&aid1)?.loc;
                let (p1_norm, ring_atom) = if ct_case == 1 {
                    let aid1_state = other.eatoms.get(&aid1)?;
                    (aid1_state.normal, aid1_state.cis_trans_nbr?)
                } else {
                    let aid1_state = self.eatoms.get(&aid1)?;
                    (aid1_state.normal, aid1_state.cis_trans_nbr?)
                };
                let r_atm_loc = if ct_case == 1 {
                    self.eatoms.get(&ring_atom)?.loc
                } else {
                    other.eatoms.get(&ring_atom)?.loc
                };
                let r_rel = (r_atm_loc.0 - p1_loc.0, r_atm_loc.1 - p1_loc.1);
                let dot = r_rel.0 * p1_norm.0 + r_rel.1 * p1_norm.1;
                if dot < 0.0 {
                    let p2_loc = self.eatoms.get(&aid2)?.loc;
                    other.reflect(p1_loc, p2_loc);
                }
            } else if common.len() == 2 {
                self.reflect_if_necessary_density(other, aid1, aid2);
            } else {
                let pt1 = self.eatoms.get(&aid1)?.loc;
                let pt2 = self.eatoms.get(&aid2)?.loc;
                let normal = (-(pt2.1 - pt1.1), pt2.0 - pt1.0);
                let oth3 = other.eatoms.get(&common[2])?.loc;
                let pt3 = self.eatoms.get(&common[2])?.loc;
                let dot1 = normal.0 * (pt3.0 - pt1.0) + normal.1 * (pt3.1 - pt1.1);
                let dot2 = normal.0 * (oth3.0 - pt1.0) + normal.1 * (oth3.1 - pt1.1);
                if dot1 * dot2 < 0.0 {
                    other.reflect(pt1, pt2);
                }
            }
        }
        for (&aid, ost) in &other.eatoms {
            if !common.contains(&aid) {
                self.eatoms.insert(aid, ost.clone());
                if !ost.pending.is_empty() && !self.attach_pts.iter().any(|&x| x == aid) {
                    self.attach_pts.push_back(aid);
                }
            } else if let Some(mst) = self.eatoms.get_mut(&aid) {
                if ost.cis_trans_nbr.is_some() {
                    mst.cis_trans_nbr = ost.cis_trans_nbr;
                    mst.normal = ost.normal;
                    mst.ccw = ost.ccw;
                }
                if ost.angle > 0.0 {
                    mst.angle = ost.angle;
                    mst.nbr1 = ost.nbr1;
                    mst.nbr2 = ost.nbr2;
                }
            }
        }
        for aid in common {
            if self.eatoms.contains_key(&aid) {
                self.update_new_neighs(aid, atoms, bonds, adjacency, degree, cip_ranks);
            }
        }
        other.done = true;
        Some(())
    }

    fn merge_frags_with_common(
        &mut self,
        efrags: &mut Vec<Self>,
        atoms: &[Atom],
        bonds: &[Bond],
        adjacency: &[Vec<usize>],
        degree: &[usize],
        cip_ranks: &[u32],
    ) -> Option<()> {
        loop {
            let found = efrags.iter().enumerate().find_map(|(idx, frag)| {
                if frag.done {
                    return None;
                }
                let common = self.find_common_atoms(frag);
                (!common.is_empty()).then_some((idx, common))
            });
            let Some((idx, common)) = found else { break };
            let common_for_cleanup = common.clone();
            let mut other = efrags.remove(idx);
            self.merge_with_common(
                &mut other, common, atoms, bonds, adjacency, degree, cip_ranks,
            )?;
            for cai in common_for_cleanup {
                let remove_attach = self
                    .eatoms
                    .get(&cai)
                    .is_some_and(|st| st.pending.is_empty())
                    && self.attach_pts.iter().any(|&x| x == cai);
                if remove_attach {
                    self.attach_pts.retain(|&x| x != cai);
                }
            }
        }
        Some(())
    }

    // RDKit✔️❌: void EmbeddedFrag::embedFusedRings(const RDKit::VECT_INT_VECT &fusedRings,
    // RDKit✔️❌:                                    bool useRingTemplates) { ... }
    fn embed_fused_rings(
        &mut self,
        atoms: &[Atom],
        bonds: &[Bond],
        adjacency: &[Vec<usize>],
        degree: &[usize],
        cip_ranks: &[u32],
        atom_ring_counts: &[usize],
        fused_rings: &[Vec<usize>],
        use_ring_templates: bool,
    ) -> Option<()> {
        let mut funion: Vec<usize> = fused_rings
            .iter()
            .flatten()
            .copied()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        if use_ring_templates
            && (fused_rings.len() > 1 || (fused_rings.len() == 1 && fused_rings[0].len() > 8))
        {
            if self.match_to_template(
                atoms,
                bonds,
                adjacency,
                degree,
                atom_ring_counts,
                &funion,
                fused_rings.len(),
            )? {
                return Some(());
            }
        }

        let mut coords = Vec::with_capacity(fused_rings.len());
        for ring in fused_rings {
            let mut ring_coords = rdkit_embed_ring(ring);
            Self::mirror_trans_ring_atoms(bonds, ring, &mut ring_coords);
            coords.push(ring_coords);
        }

        let mut done_rings = Vec::<usize>::new();
        if use_ring_templates {
            let (core_rings, core_ring_ids) = rdkit_find_core_rings(fused_rings, bonds);
            if core_rings.len() > 1 && core_rings.len() < fused_rings.len() {
                let core_union: Vec<usize> = core_rings
                    .iter()
                    .flatten()
                    .copied()
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .collect();
                if self.match_to_template(
                    atoms,
                    bonds,
                    adjacency,
                    degree,
                    atom_ring_counts,
                    &core_union,
                    core_rings.len(),
                )? {
                    done_rings = core_ring_ids;
                }
            }
        }

        if done_rings.is_empty() {
            let first_ring_id = rdkit_pick_first_ring_to_embed(degree, fused_rings);
            self.init_from_ring_coords(&fused_rings[first_ring_id], &coords[first_ring_id]);
            done_rings.push(first_ring_id);
        }
        funion = fused_rings
            .iter()
            .flatten()
            .copied()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        while self.eatoms.len() < funion.len() {
            let (next_id, common_atom_ids) =
                rdkit_find_next_ring_to_embed(&done_rings, fused_rings)?;
            let mut emb_ring = Self::default();
            emb_ring.init_from_ring_coords(&fused_rings[next_id], &coords[next_id]);
            let mut pin_atoms = Vec::new();
            if common_atom_ids.len() == 1 {
                let aid = common_atom_ids[0];
                let trans = self.compute_one_atom_trans(aid, &emb_ring)?;
                emb_ring.transform(trans);
                pin_atoms.push(aid);
            } else {
                let aid1 = *common_atom_ids.first()?;
                let aid2 = *common_atom_ids.last()?;
                pin_atoms.push(aid1);
                pin_atoms.push(aid2);
                let trans = self.compute_two_atom_trans(aid1, aid2, &coords[next_id])?;
                emb_ring.transform(trans);
                self.reflect_if_necessary_density(&mut emb_ring, aid1, aid2);
            }
            self.merge_ring(&emb_ring, common_atom_ids.len(), &pin_atoms);
            done_rings.push(next_id);
        }
        // RDKit✔️✔️: constructor path for fused rings returns to embedFusedSystems(),
        // RDKit✔️✔️: which then calls setupNewNeighs() only. setupAttachmentPoints()
        // RDKit✔️✔️: is intentionally not run here because ring traversal order has
        // RDKit✔️✔️: already initialized nbr1/nbr2/angle for ring atoms.
        self.setup_new_neighs(atoms, bonds, adjacency, degree, cip_ranks);
        Some(())
    }
}

fn build_rdkit_template_query_molecule(
    template: &RdkitTemplateRuntimeModel,
) -> Result<Molecule, Coordinate2DError> {
    let mut builder = MoleculeBuilder::new();
    let atom_ids: Vec<_> = template
        .graph
        .atom_queries
        .iter()
        .map(|query| {
            builder.add_atom(
                AtomSpec::new(Element::from_atomic_number(0).unwrap()).with_query(query.clone()),
            )
        })
        .collect();
    for bond in &template.graph.bonds {
        builder
            .add_bond(
                BondSpec::new(
                    atom_ids[bond.begin_atom_idx],
                    atom_ids[bond.end_atom_idx],
                    bond_order_for_template_probe(&bond.query),
                )
                .with_query(bond.query.clone()),
            )
            .map_err(|error| {
                Coordinate2DError::UnsupportedFeature(format!(
                    "RDKit template query build failed: {error}"
                ))
            })?;
    }
    if let Some(coords) = &template.coords_2d {
        builder
            .set_2d_coordinates(coords.clone())
            .map_err(|error| {
                Coordinate2DError::UnsupportedFeature(format!(
                    "RDKit template query coordinate build failed: {error}"
                ))
            })?;
    }
    builder.build().map_err(|error| {
        Coordinate2DError::UnsupportedFeature(format!(
            "RDKit template query molecule build failed: {error}"
        ))
    })
}

fn build_rdkit_ring_system_molecule(
    atoms: &[Atom],
    bonds: &[Bond],
    ring_system_atoms: &[usize],
) -> Result<Molecule, Coordinate2DError> {
    let ring_set: BTreeSet<usize> = ring_system_atoms.iter().copied().collect();
    let mut builder = MoleculeBuilder::new();
    let mut atom_map = BTreeMap::<usize, crate::AtomId>::new();
    for &aid in ring_system_atoms {
        atom_map.insert(aid, builder.add_atom(rdkit_atom_to_spec(&atoms[aid])));
    }
    for bond in bonds {
        let begin = bond.begin().index();
        let end = bond.end().index();
        if !ring_set.contains(&begin) || !ring_set.contains(&end) {
            continue;
        }
        let mut spec = BondSpec::new(atom_map[&begin], atom_map[&end], bond.order())
            .with_stereo(bond.stereo())
            .with_aromatic(bond.is_aromatic())
            .with_conjugated(bond.is_conjugated())
            .with_direction(bond.direction())
            .with_unknown_stereo(bond.unknown_stereo());
        if let Some([a0, a1]) = bond.stereo_atoms() {
            if let (Some(&mapped0), Some(&mapped1)) =
                (atom_map.get(&a0.index()), atom_map.get(&a1.index()))
            {
                spec = spec.with_stereo_atoms(mapped0, mapped1);
            }
        }
        if let Some(query) = bond.query() {
            spec = spec.with_query(query.clone());
        }
        builder.add_bond(spec).map_err(|error| {
            Coordinate2DError::UnsupportedFeature(format!(
                "RDKit ring-system molecule build failed: {error}"
            ))
        })?;
    }
    builder.build().map_err(|error| {
        Coordinate2DError::UnsupportedFeature(format!(
            "RDKit ring-system molecule finalization failed: {error}"
        ))
    })
}

fn build_rdkit_molecule_from_slices(
    atoms: &[Atom],
    bonds: &[Bond],
) -> Result<Molecule, Coordinate2DError> {
    let mut builder = MoleculeBuilder::new();
    let mut atom_map = Vec::with_capacity(atoms.len());
    for atom in atoms {
        atom_map.push(builder.add_atom(rdkit_atom_to_spec(atom)));
    }
    for bond in bonds {
        let begin = bond.begin().index();
        let end = bond.end().index();
        let mut spec = BondSpec::new(atom_map[begin], atom_map[end], bond.order())
            .with_stereo(bond.stereo())
            .with_aromatic(bond.is_aromatic())
            .with_conjugated(bond.is_conjugated())
            .with_direction(bond.direction())
            .with_unknown_stereo(bond.unknown_stereo());
        if let Some([a0, a1]) = bond.stereo_atoms() {
            spec = spec.with_stereo_atoms(atom_map[a0.index()], atom_map[a1.index()]);
        }
        if let Some(query) = bond.query() {
            spec = spec.with_query(query.clone());
        }
        builder.add_bond(spec).map_err(|error| {
            Coordinate2DError::UnsupportedFeature(format!(
                "RDKit whole-molecule build failed: {error}"
            ))
        })?;
    }
    builder.build().map_err(|error| {
        Coordinate2DError::UnsupportedFeature(format!(
            "RDKit whole-molecule finalization failed: {error}"
        ))
    })
}

fn build_rdkit_depict_molecule_from_slices(
    atoms: &[Atom],
    bonds: &[Bond],
) -> Result<Molecule, Coordinate2DError> {
    let mut molecule = build_rdkit_molecule_from_slices(atoms, bonds)?;
    let valence =
        crate::assign_valence_with_options(&molecule, crate::ValenceModel::RdkitLike, false)
        .map_err(|error| {
            Coordinate2DError::UnsupportedFeature(format!(
                "RDKit updatePropertyCache(false) equivalent failed during 2D coordinate generation: {error}"
            ))
        })?;
    molecule.derived_cache_mut().valence = Some(valence);
    let rings = crate::symmetrize_sssr(&molecule).map_err(|error| {
        Coordinate2DError::UnsupportedFeature(format!(
            "RDKit GetSymmSSSR equivalent failed during 2D coordinate generation: {error}"
        ))
    })?;
    molecule.derived_cache_mut().rings = Some(rings);
    Ok(molecule)
}

fn rdkit_atom_to_spec(atom: &Atom) -> AtomSpec {
    let mut spec = AtomSpec::new(atom.element())
        .with_formal_charge(atom.formal_charge())
        .with_explicit_hydrogens(atom.explicit_hydrogens())
        .with_chiral_tag(atom.chiral_tag())
        .with_aromatic(atom.is_aromatic())
        .with_no_implicit(atom.no_implicit())
        .with_radical_electrons(atom.radical_electrons())
        .with_hybridization(atom.hybridization());
    if let Some(isotope) = atom.isotope() {
        spec = spec.with_isotope(isotope);
    }
    if let Some(atom_map) = atom.atom_map() {
        spec = spec.with_atom_map(atom_map);
    }
    if let Some(chiral_perm) = atom.chiral_permutation() {
        spec = spec.with_chiral_permutation(chiral_perm);
    }
    if atom.unknown_stereo() {
        spec = spec.with_unknown_stereo(true);
    }
    if let Some(mol_parity) = atom.mol_parity() {
        spec = spec.with_mol_parity(mol_parity);
    }
    if let Some(mol_inv) = atom.mol_inversion_flag() {
        spec = spec.with_mol_inversion_flag(mol_inv);
    }
    if atom.implicit_hydrogen() {
        spec = spec.with_implicit_hydrogen(true);
    }
    let tracked = atom.tracked_isotopic_hydrogens().to_vec();
    if !tracked.is_empty() {
        spec = spec.with_tracked_isotopic_hydrogens(tracked);
    }
    if let Some(query) = atom.query() {
        spec = spec.with_query(query.clone());
    }
    for (key, value) in atom.props() {
        spec = spec.with_prop(key.clone(), value.clone());
    }
    if let Some(pdb_info) = atom.pdb_residue_info() {
        spec = spec.with_pdb_residue_info(pdb_info.clone());
    }
    spec
}

fn rdkit_template_degree_counts(mol: &Molecule, dummy_atomic_num: Option<u8>) -> [i32; 5] {
    let adjacency = crate::AdjacencyList::from_topology(mol.num_atoms(), mol.bonds());
    let mut counts = [0i32; 5];
    for atom in mol.atoms() {
        if dummy_atomic_num.is_some_and(|anum| atom.atomic_number() == anum) {
            continue;
        }
        let mut degree = 0usize;
        for nbr in adjacency.neighbors_of(atom.id().index()) {
            if dummy_atomic_num
                .is_some_and(|anum| mol.atoms()[nbr.atom_index].atomic_number() == anum)
            {
                continue;
            }
            degree += 1;
            if degree == 4 {
                break;
            }
        }
        counts[degree] += 1;
    }
    counts
}

fn rdkit_check_stereo_chemistry(
    mol: &Molecule,
    template: &RdkitTemplateRuntimeModel,
    atom_mapping: &[usize],
) -> bool {
    // RDKit✔️❌: static bool checkStereoChemistry(const RDKit::ROMol &mol,
    // RDKit✔️❌:                                  const RDKit::ROMol &template_mol,
    // RDKit✔️❌:                                  RDKit::MatchVectType match) { ... }
    let Some(coords) = &template.coords_2d else {
        return false;
    };
    let inv_map: BTreeMap<usize, usize> = atom_mapping
        .iter()
        .enumerate()
        .map(|(template_idx, &mol_idx)| (mol_idx, template_idx))
        .collect();
    let adjacency = crate::AdjacencyList::from_topology(mol.num_atoms(), mol.bonds());
    for bond in mol.bonds() {
        if bond.order() != BondOrder::Double
            || matches!(bond.stereo(), BondStereo::Any | BondStereo::None)
        {
            continue;
        }
        let Some(stereo_atoms) = bond.stereo_atoms() else {
            continue;
        };
        let atom1_neighbor1 = stereo_atoms[0].index();
        let atom2_neighbor1 = stereo_atoms[1].index();
        let atom1 = bond.begin().index();
        let atom2 = bond.end().index();

        let atom1_neighbor2 = if adjacency.neighbors_of(bond.begin().index()).len() > 2 {
            adjacency
                .neighbors_of(bond.begin().index())
                .iter()
                .map(|a| a.atom_index)
                .find(|&nbr| nbr != atom1_neighbor1 && nbr != atom2)
        } else {
            None
        };
        let atom2_neighbor2 = if adjacency.neighbors_of(bond.end().index()).len() > 2 {
            adjacency
                .neighbors_of(bond.end().index())
                .iter()
                .map(|a| a.atom_index)
                .find(|&nbr| nbr != atom2_neighbor1 && nbr != atom1)
        } else {
            None
        };

        let template_atom1 = inv_map.get(&atom1).copied();
        let template_atom2 = inv_map.get(&atom2).copied();
        let mut template_atom1_neighbor1 = inv_map.get(&atom1_neighbor1).copied();
        let mut template_atom2_neighbor1 = inv_map.get(&atom2_neighbor1).copied();
        let template_atom1_neighbor2 = atom1_neighbor2.and_then(|idx| inv_map.get(&idx).copied());
        let template_atom2_neighbor2 = atom2_neighbor2.and_then(|idx| inv_map.get(&idx).copied());

        let mut swap_stereo = false;
        if template_atom1_neighbor1.is_none() {
            template_atom1_neighbor1 = template_atom1_neighbor2;
            swap_stereo = !swap_stereo;
        }
        if template_atom2_neighbor1.is_none() {
            template_atom2_neighbor1 = template_atom2_neighbor2;
            swap_stereo = !swap_stereo;
        }
        let (
            Some(template_atom1),
            Some(template_atom2),
            Some(template_atom1_neighbor1),
            Some(template_atom2_neighbor1),
        ) = (
            template_atom1,
            template_atom2,
            template_atom1_neighbor1,
            template_atom2_neighbor1,
        )
        else {
            return false;
        };
        let atom1_loc = coords[template_atom1];
        let atom2_loc = coords[template_atom2];
        let atom1_neighbor_loc = coords[template_atom1_neighbor1];
        let atom2_neighbor_loc = coords[template_atom2_neighbor1];
        let v12 = (
            atom1_neighbor_loc[0] - atom1_loc[0],
            atom1_neighbor_loc[1] - atom1_loc[1],
        );
        let v42 = (
            atom2_neighbor_loc[0] - atom1_loc[0],
            atom2_neighbor_loc[1] - atom1_loc[1],
        );
        let v32 = (atom2_loc[0] - atom1_loc[0], atom2_loc[1] - atom1_loc[1]);
        let cross1 = v32.0 * v12.1 - v32.1 * v12.0;
        let cross2 = v32.0 * v42.1 - v32.1 * v42.0;
        let mut is_cis = cross1 * cross2 > 0.0;
        if swap_stereo {
            is_cis = !is_cis;
        }
        let bond_is_cis = matches!(bond.stereo(), BondStereo::Z | BondStereo::Cis);
        if is_cis != bond_is_cis {
            return false;
        }
    }
    true
}

// ── Pure functions ───────────────────────────────────────────────────────────

// BEGIN RDKIT CPP FUNCTION: atomDepictRank (DepictUtils.h)
// RDKit✔️✔️: unsigned int atomDepictRank(unsigned int atomicNum, unsigned int degree)
fn atom_depict_rank(atomic_num: u8, degree: usize) -> usize {
    let anum = if atomic_num == 1 {
        1000usize
    } else {
        atomic_num as usize
    };
    100 * anum + degree
}

// BEGIN RDKIT CPP FUNCTION: ringRadius (EmbeddedFrag.cpp)
// RDKit✔️✔️: ringRadius(unsigned ringSize, double bondLen)
fn rdkit_ring_radius(ring_size: usize, bond_len: f64) -> f64 {
    let ang = 2.0 * PI / ring_size as f64;
    bond_len / rdkit_sqrt(2.0 * (1.0 - rdkit_cos(ang)))
}

// BEGIN RDKIT CPP FUNCTION: numBondsPlusLonePairs (ConjugHybrid.cpp)
// RDKit✔️❌: adapted for new-core API
fn rdkit_num_bonds_plus_lone_pairs(
    atomic_num: u8,
    graph_degree: usize,
    explicit_hydrogens: u8,
    explicit_valence: i32,
    implicit_hydrogens: i32,
    radical_electrons: u8,
    formal_charge: i8,
    zero_or_outgoing_dative_bonds: usize,
) -> Option<i32> {
    let deg = graph_degree as i32 + i32::from(explicit_hydrogens) + implicit_hydrogens
        - zero_or_outgoing_dative_bonds as i32;

    if atomic_num <= 1 {
        return Some(deg);
    }
    let nouter = n_outer_electrons_rdkit(atomic_num)?;
    let total_valence = explicit_valence + implicit_hydrogens;
    let charge = formal_charge as i32;
    let free_electrons = nouter - (total_valence + charge);
    if total_valence + nouter - charge < 8 {
        let radicals = radical_electrons as i32;
        Some(deg + (free_electrons - radicals) / 2 + radicals)
    } else {
        Some(deg + free_electrons / 2)
    }
}

// BEGIN RDKIT CPP FUNCTION: ComputeSubAngle (DepictUtils.h)
// RDKit✔️✔️: double computeSubAngle(unsigned int degree, HybridizationType type)
fn compute_sub_angle(degree: usize, hybridization: RdkitHybridization) -> f64 {
    match hybridization {
        RdkitHybridization::Unspecified | RdkitHybridization::Sp3 => {
            if degree == 4 {
                PI / 2.0
            } else {
                2.0 * PI / 3.0
            }
        }
        RdkitHybridization::Sp2 => 2.0 * PI / 3.0,
        RdkitHybridization::S
        | RdkitHybridization::Sp
        | RdkitHybridization::Sp3d
        | RdkitHybridization::Sp3d2 => 2.0 * PI / degree as f64,
    }
}

// BEGIN RDKIT CPP FUNCTION (standalone): canonicalizeComponent
// RDKit✔️✔️: (RDKit DepictUtils.cpp)
fn canonicalize_component(component: &mut Vec<(usize, (f64, f64))>) {
    if component.len() <= 1 {
        return;
    }
    let n = component.len() as f64;
    let (mut cx, mut cy) = (0.0f64, 0.0f64);
    for &(_, (x, y)) in component.iter() {
        cx += x;
        cy += y;
    }
    let inv_n = 1.0 / n;
    cx *= inv_n;
    cy *= inv_n;

    let (mut xx, mut xy, mut yy) = (0.0f64, 0.0f64, 0.0f64);
    let mut centered: Vec<(f64, f64)> = Vec::with_capacity(component.len());
    for &(_, (x, y)) in component.iter() {
        let px = x - cx;
        let py = y - cy;
        centered.push((px, py));
        xx += px * px;
        xy += px * py;
        yy += py * py;
    }

    let d = ((xx - yy) * (xx - yy) + 4.0 * xy * xy).sqrt();
    let mut eig1 = (2.0 * xy, (yy - xx) + d);
    let eig1_len = norm(eig1);
    if eig1_len <= 1e-4 {
        for (i, (_, pos)) in component.iter_mut().enumerate() {
            *pos = centered[i];
        }
        return;
    }
    let e_val1 = (xx + yy + d) / 2.0;
    eig1 = (eig1.0 / eig1_len, eig1.1 / eig1_len);

    let mut eig2 = (2.0 * xy, (yy - xx) - d);
    let e_val2 = (xx + yy - d) / 2.0;
    let eig2_len = norm(eig2);
    if eig2_len > 1e-4 {
        eig2 = (eig2.0 / eig2_len, eig2.1 / eig2_len);
        if e_val2 > e_val1 {
            std::mem::swap(&mut eig1, &mut eig2);
        }
    }

    for (i, (_, pos)) in component.iter_mut().enumerate() {
        let (px, py) = centered[i];
        // trans = [[eig1.x, eig1.y],[-eig1.y,eig1.x]]
        let rx = px * eig1.0 + py * eig1.1;
        let ry = -px * eig1.1 + py * eig1.0;
        *pos = (rx, ry);
    }
}

fn component_box_rdkit(component: &[(usize, (f64, f64))]) -> (f64, f64, f64, f64) {
    let mut px = -1.0e8f64;
    let mut nx = 1.0e8f64;
    let mut py = -1.0e8f64;
    let mut ny = 1.0e8f64;
    for &(_, (x, y)) in component {
        px = px.max(x);
        nx = nx.min(x);
        py = py.max(y);
        ny = ny.min(y);
    }
    if component.is_empty() {
        (0.0, 0.0, 0.0, 0.0)
    } else {
        (px, -nx, py, -ny)
    }
}

// ── Outer electrons ─────────────────────────────────────────────────────────

/// RDKit❗✔️: outer electron count from atomic number — RDKit PeriodicTable::getNouter
fn n_outer_electrons_rdkit(atomic_num: u8) -> Option<i32> {
    periodic_table_outer_electrons(atomic_num).ok()
}

// ── Hybridization for depiction ─────────────────────────────────────────────

// BEGIN RDKIT CPP FUNCTION: getHybridization (RDDepictor.cpp / ConjugHybrid.cpp)
// RDKit✔️❌: adapted for new-core API
fn rdkit_hybridizations_for_depict(
    atoms: &[Atom],
    bonds: &[Bond],
    degree: &[usize],
) -> Result<Vec<RdkitHybridization>, Coordinate2DError> {
    // RDKit✔️❌: void setHybridization(ROMol &mol) {
    // RDKit✔️❌:   for (auto atom : mol.atoms()) {
    // RDKit✔️❌:     if (atom->getAtomicNum() == 0) {
    // RDKit✔️❌:       atom->setHybridization(Atom::UNSPECIFIED);
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       switch (atom->getChiralTag()) {
    // RDKit✔️❌:         case Atom::ChiralType::CHI_TETRAHEDRAL:
    // RDKit✔️❌:         case Atom::ChiralType::CHI_TETRAHEDRAL_CW:
    // RDKit✔️❌:         case Atom::ChiralType::CHI_TETRAHEDRAL_CCW:
    // RDKit✔️❌:           if (atom->getTotalDegree() == 4) {
    // RDKit✔️❌:             atom->setHybridization(Atom::HybridizationType::SP3);
    // RDKit✔️❌:             continue;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case Atom::ChiralType::CHI_SQUAREPLANAR:
    // RDKit✔️❌:           if (atom->getTotalDegree() <= 4 && atom->getTotalDegree() >= 2) {
    // RDKit✔️❌:             atom->setHybridization(Atom::HybridizationType::SP2D);
    // RDKit✔️❌:             continue;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
    // RDKit✔️❌:           if (atom->getTotalDegree() <= 5 && atom->getTotalDegree() >= 2) {
    // RDKit✔️❌:             atom->setHybridization(Atom::HybridizationType::SP3D);
    // RDKit✔️❌:             continue;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case Atom::ChiralType::CHI_OCTAHEDRAL:
    // RDKit✔️❌:           if (atom->getTotalDegree() <= 6 && atom->getTotalDegree() >= 2) {
    // RDKit✔️❌:             atom->setHybridization(Atom::HybridizationType::SP3D2);
    // RDKit✔️❌:             continue;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         default:
    // RDKit✔️❌:           break;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       int norbs;
    // RDKit✔️❌:       if (atom->getAtomicNum() < 89) {
    // RDKit✔️❌:         norbs = numBondsPlusLonePairs(atom);
    // RDKit✔️❌:       } else {
    // RDKit✔️❌:         norbs = atom->getTotalDegree();
    // RDKit✔️❌:       }
    // RDKit✔️❌:       switch (norbs) {
    // RDKit✔️❌:         case 0:
    // RDKit✔️❌:         case 1:
    // RDKit✔️❌:           atom->setHybridization(Atom::S);
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case 2:
    // RDKit✔️❌:           atom->setHybridization(Atom::SP);
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case 3:
    // RDKit✔️❌:           atom->setHybridization(Atom::SP2);
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case 4:
    // RDKit✔️❌:           if (atom->getTotalDegree() > 3 ||
    // RDKit✔️❌:               !MolOps::atomHasConjugatedBond(atom)) {
    // RDKit✔️❌:             atom->setHybridization(Atom::SP3);
    // RDKit✔️❌:           } else {
    // RDKit✔️❌:             atom->setHybridization(Atom::SP2);
    // RDKit✔️❌:           }
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case 5:
    // RDKit✔️❌:           atom->setHybridization(Atom::SP3D);
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case 6:
    // RDKit✔️❌:           atom->setHybridization(Atom::SP3D2);
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         default:
    // RDKit✔️❌:           atom->setHybridization(Atom::UNSPECIFIED);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    let depict_mol = build_rdkit_depict_molecule_from_slices(atoms, bonds)?;
    let assignment =
        crate::assign_valence(&depict_mol, crate::ValenceModel::RdkitLike).map_err(|error| {
            Coordinate2DError::UnsupportedFeature(format!(
                "RDKit hybridization valence assignment failed: {error}"
            ))
        })?;
    let radicals = crate::assign_radicals(&depict_mol).map_err(|error| {
        Coordinate2DError::UnsupportedFeature(format!(
            "RDKit hybridization radical assignment failed: {error}"
        ))
    })?;

    let mut out = Vec::with_capacity(atoms.len());
    for (idx, atom) in atoms.iter().enumerate() {
        if atom.atomic_number() == 0 {
            out.push(RdkitHybridization::Unspecified);
            continue;
        }

        let total_degree = degree[idx] as i32
            + i32::from(atom.explicit_hydrogens())
            + assignment.implicit_hydrogens[idx];
        match atom.chiral_tag() {
            ChiralTag::Tetrahedral | ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
                if total_degree == 4 =>
            {
                out.push(RdkitHybridization::Sp3);
                continue;
            }
            ChiralTag::TrigonalBipyramidal if (2..=5).contains(&total_degree) => {
                out.push(RdkitHybridization::Sp3d);
                continue;
            }
            ChiralTag::Octahedral if (2..=6).contains(&total_degree) => {
                out.push(RdkitHybridization::Sp3d2);
                continue;
            }
            _ => {}
        }

        let zero_or_outgoing_dative_bonds = bonds
            .iter()
            .filter(|bond| {
                let touches_atom = bond.begin().index() == idx || bond.end().index() == idx;
                touches_atom
                    && (matches!(bond.order(), BondOrder::Zero)
                        || (matches!(bond.order(), BondOrder::Dative) && bond.end().index() != idx))
            })
            .count();

        let norbs = if atom.atomic_number() < 89 {
            rdkit_num_bonds_plus_lone_pairs(
                atom.atomic_number(),
                degree[idx],
                atom.explicit_hydrogens(),
                assignment.explicit_valence[idx] as i32,
                assignment.implicit_hydrogens[idx],
                radicals[idx],
                atom.formal_charge(),
                zero_or_outgoing_dative_bonds,
            )
            .ok_or_else(|| {
                Coordinate2DError::UnsupportedFeature(format!(
                    "RDKit hybridization outer-electron lookup failed for atom {idx}"
                ))
            })?
        } else {
            degree[idx] as i32
                + atom.explicit_hydrogens() as i32
                + assignment.implicit_hydrogens[idx]
                - zero_or_outgoing_dative_bonds as i32
        };
        // Fallback: if norbs is 0 or negative (lookup failure), use Unspecified
        if norbs <= 0 {
            out.push(RdkitHybridization::Unspecified);
            continue;
        }

        let has_conjugated_bond = bonds.iter().any(|bond| {
            (bond.begin().index() == idx || bond.end().index() == idx) && bond.is_conjugated()
        });

        out.push(match norbs {
            0 | 1 => RdkitHybridization::S,
            2 => RdkitHybridization::Sp,
            3 => RdkitHybridization::Sp2,
            4 => {
                if total_degree > 3 || !has_conjugated_bond {
                    RdkitHybridization::Sp3
                } else {
                    RdkitHybridization::Sp2
                }
            }
            5 => RdkitHybridization::Sp3d,
            6 => RdkitHybridization::Sp3d2,
            _ => RdkitHybridization::Unspecified,
        });
    }
    Ok(out)
}

// ── Ranking helpers ─────────────────────────────────────────────────────────

// BEGIN RDKIT CPP FUNCTION: rankAtomsByRank (DepictUtils.cpp)
// RDKit✔️✔️: template <class T>
// RDKit✔️✔️: T rankAtomsByRank(const RDKit::ROMol &mol, const T &commAtms, bool ascending) {
// RDKit✔️✔️:   const auto natms = commAtms.size();
// RDKit✔️✔️:   INT_PAIR_VECT rankAid;
// RDKit✔️✔️:   rankAid.reserve(natms);
// RDKit✔️✔️:   for (const auto aid : commAtms) {
// RDKit✔️✔️:     unsigned int rank = aid;
// RDKit✔️✔️:     const auto at = mol.getAtomWithIdx(aid);
// RDKit✔️✔️:     if (!at->getPropIfPresent(RDKit::common_properties::_CIPRank, rank)) {
// RDKit✔️✔️:       if (at->getPropIfPresent(RDKit::common_properties::_ChiralAtomRank, rank)) {
// RDKit✔️✔️:         rank = mol.getNumAtoms() - rank;
// RDKit✔️✔️:       }
// RDKit✔️✔️:       rank += mol.getNumAtoms() * getAtomDepictRank(at);
// RDKit✔️✔️:     }
// RDKit✔️✔️:     rankAid.emplace_back(rank, aid);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (ascending) {
// RDKit✔️✔️:     std::stable_sort(rankAid.begin(), rankAid.end(),
// RDKit✔️✔️:                      [](const auto &e1, const auto &e2) { return e1 < e2; });
// RDKit✔️✔️:   } else {
// RDKit✔️✔️:     std::stable_sort(rankAid.begin(), rankAid.end(),
// RDKit✔️✔️:                      [](const auto &e1, const auto &e2) { return e1 > e2; });
// RDKit✔️✔️:   }
// RDKit✔️✔️:   T res;
// RDKit✔️✔️:   std::for_each(rankAid.begin(), rankAid.end(),
// RDKit✔️✔️:                 [&res](const auto &elem) { res.push_back(elem.second); });
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
fn rdkit_rank_atoms_by_rank_into(
    atom_slice: &[Atom],
    order: &mut [usize],
    degree: &[usize],
    ascending: bool,
) {
    let mol_atom_count = atom_slice.len() as u32;
    let mut rank_aid: Vec<(u32, usize)> = Vec::with_capacity(order.len());
    for &aid in order.iter() {
        let mut rank = aid as u32;
        let atom = &atom_slice[aid];
        if let Some(cip_rank) = atom
            .prop("_CIPRank")
            .and_then(|value| value.parse::<u32>().ok())
        {
            rank = cip_rank;
        } else if let Some(chiral_rank) = atom
            .prop("_ChiralAtomRank")
            .and_then(|value| value.parse::<u32>().ok())
        {
            rank = mol_atom_count - chiral_rank;
        } else {
            rank += mol_atom_count * atom_depict_rank(atom.atomic_number(), degree[aid]) as u32;
        }
        rank_aid.push((rank, aid));
    }
    if ascending {
        rank_aid.sort_by(|e1, e2| e1.cmp(e2));
    } else {
        rank_aid.sort_by(|e1, e2| e2.cmp(e1));
    }
    for (dst, (_, aid)) in order.iter_mut().zip(rank_aid.into_iter()) {
        *dst = aid;
    }
}

fn rdkit_rank_atoms_by_rank(
    atom_slice: &[Atom],
    order: &[usize],
    degree: &[usize],
    ascending: bool,
) -> Vec<usize> {
    let mut res = order.to_vec();
    rdkit_rank_atoms_by_rank_into(atom_slice, &mut res, degree, ascending);
    res
}

// BEGIN RDKIT CPP FUNCTION: setNbrOrder (DepictUtils.cpp)
// RDKit✔️✔️: RDKit::INT_VECT setNbrOrder(unsigned int aid, const RDKit::INT_VECT &nbrs,
// RDKit✔️✔️:                                   const RDKit::ROMol &mol) {
// RDKit✔️✔️:   RDKit::INT_VECT thold = nbrs;
// RDKit✔️✔️:   int ref = -1;
// RDKit✔️✔️:   for (auto n : mol.atomNeighbors(mol.getAtomWithIdx(aid))) {
// RDKit✔️✔️:     if (std::find(nbrs.begin(), nbrs.end(), n->getIdx()) == nbrs.end()) {
// RDKit✔️✔️:       ref = n->getIdx();
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (ref > -1) {
// RDKit✔️✔️:     thold.push_back(ref);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (thold.size() <= 3) {
// RDKit✔️✔️:     return rankAtomsByRank(mol, nbrs);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   thold = rankAtomsByRank(mol, thold);
// RDKit✔️✔️:   unsigned int ln = thold.size();
// RDKit✔️✔️:   std::swap(thold[ln - 3], thold[ln - 2]);
// RDKit✔️✔️:   if (ref > -1) {
// RDKit✔️✔️:     auto cptr = std::find(thold.begin(), thold.end(), ref);
// RDKit✔️✔️:     PRECONDITION(cptr != thold.end(), "");
// RDKit✔️✔️:     RDKit::INT_VECT res;
// RDKit✔️✔️:     cptr++;
// RDKit✔️✔️:     while (cptr != thold.end()) {
// RDKit✔️✔️:       res.push_back(*cptr);
// RDKit✔️✔️:       cptr++;
// RDKit✔️✔️:     }
// RDKit✔️✔️:     cptr = thold.begin();
// RDKit✔️✔️:     while ((*cptr) != ref) {
// RDKit✔️✔️:       res.push_back(*cptr);
// RDKit✔️✔️:       cptr++;
// RDKit✔️✔️:     }
// RDKit✔️✔️:     return res;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return thold;
// RDKit✔️✔️: }
fn rdkit_set_nbr_order(
    aid: usize,
    nbrs: &[usize],
    atoms: &[Atom],
    _bonds: &[Bond],
    adjacency: &[Vec<usize>],
    degree: &[usize],
) -> Vec<usize> {
    let mut thold = nbrs.to_vec();
    let mut ref_atom: Option<usize> = None;
    for &nb in &adjacency[aid] {
        if !nbrs.contains(&nb) {
            ref_atom = Some(nb);
        }
    }
    if let Some(r) = ref_atom {
        thold.push(r);
    }
    if thold.len() <= 3 {
        return rdkit_rank_atoms_by_rank(atoms, nbrs, degree, true);
    }

    thold = rdkit_rank_atoms_by_rank(atoms, &thold, degree, true);
    let ln = thold.len();
    thold.swap(ln - 3, ln - 2);

    if let Some(r) = ref_atom {
        if let Some(pos) = thold.iter().position(|&a| a == r) {
            let mut out = Vec::with_capacity(nbrs.len());
            out.extend_from_slice(&thold[pos + 1..]);
            out.extend_from_slice(&thold[..pos]);
            return out;
        }
    }
    thold
}

// BEGIN RDKIT CPP FUNCTION: getNbrAtomAndBondIds (DepictUtils.cpp)
// RDKit✔️✔️: void getNbrAtomAndBondIds(unsigned int aid, const RDKit::ROMol *mol,
// RDKit✔️✔️:                           RDKit::INT_VECT &aids, RDKit::INT_VECT &bids) {
// RDKit✔️✔️:   CHECK_INVARIANT(mol, "");
// RDKit✔️✔️:   unsigned int na = mol->getNumAtoms();
// RDKit✔️✔️:   URANGE_CHECK(aid, na);
// RDKit✔️✔️:   for (auto nbr : mol->atomNeighbors(mol->getAtomWithIdx(aid))) {
// RDKit✔️✔️:     auto bi = mol->getBondBetweenAtoms(aid, nbr->getIdx())->getIdx();
// RDKit✔️✔️:     aids.push_back(nbr->getIdx());
// RDKit✔️✔️:     bids.push_back(bi);
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
fn rdkit_get_nbr_atom_and_bond_ids(
    aid: usize,
    bonds: &[Bond],
    adjacency: &[Vec<usize>],
) -> (Vec<usize>, Vec<usize>) {
    let mut aids = Vec::new();
    let mut bids = Vec::new();
    for &nbr in &adjacency[aid] {
        if let Some((bid, _)) = bonds.iter().enumerate().find(|(_, bond)| {
            let begin = bond.begin().index();
            let end = bond.end().index();
            (begin == aid && end == nbr) || (begin == nbr && end == aid)
        }) {
            aids.push(nbr);
            bids.push(bid);
        }
    }
    (aids, bids)
}

// BEGIN RDKIT CPP FUNCTION: findBondsPairsToPermuteDeg4 (DepictUtils.cpp)
// RDKit✔️✔️: INT_PAIR_VECT findBondsPairsToPermuteDeg4(const RDGeom::Point2D &center,
// RDKit✔️✔️:                                           const RDKit::INT_VECT &nbrBids,
// RDKit✔️✔️:                                           const VECT_C_POINT &nbrLocs) {
// RDKit✔️✔️:   INT_PAIR_VECT res;
// RDKit✔️✔️:   CHECK_INVARIANT(nbrBids.size() == 4, "");
// RDKit✔️✔️:   CHECK_INVARIANT(nbrLocs.size() == 4, "");
// RDKit✔️✔️:   std::vector<RDGeom::Point2D> nbrPts;
// RDKit✔️✔️:   nbrPts.reserve(nbrLocs.size());
// RDKit✔️✔️:   for (const auto &nloc : nbrLocs) {
// RDKit✔️✔️:     RDGeom::Point2D v = (*nloc) - center;
// RDKit✔️✔️:     nbrPts.push_back(v);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   double dp1 = nbrPts[0].dotProduct(nbrPts[1]);
// RDKit✔️✔️:   if (fabs(dp1) < 1.e-3) {
// RDKit✔️✔️:     INT_PAIR p1(nbrBids[0], nbrBids[1]);
// RDKit✔️✔️:     res.push_back(p1);
// RDKit✔️✔️:     double dp2 = nbrPts[0].dotProduct(nbrPts[2]);
// RDKit✔️✔️:     if (fabs(dp2) < 1.e-3) {
// RDKit✔️✔️:       INT_PAIR p2(nbrBids[0], nbrBids[2]);
// RDKit✔️✔️:       res.push_back(p2);
// RDKit✔️✔️:     } else {
// RDKit✔️✔️:       INT_PAIR p2(nbrBids[0], nbrBids[3]);
// RDKit✔️✔️:       res.push_back(p2);
// RDKit✔️✔️:     }
// RDKit✔️✔️:     return res;
// RDKit✔️✔️:   } else {
// RDKit✔️✔️:     INT_PAIR p1(nbrBids[0], nbrBids[2]);
// RDKit✔️✔️:     res.push_back(p1);
// RDKit✔️✔️:     INT_PAIR p2(nbrBids[0], nbrBids[3]);
// RDKit✔️✔️:     res.push_back(p2);
// RDKit✔️✔️:     return res;
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
fn rdkit_find_bond_pairs_to_permute_deg4(
    center: (f64, f64),
    nbr_bids: &[usize],
    nbr_locs: &[(f64, f64)],
) -> Vec<(usize, usize)> {
    assert_eq!(nbr_bids.len(), 4);
    assert_eq!(nbr_locs.len(), 4);

    let nbr_pts: Vec<(f64, f64)> = nbr_locs
        .iter()
        .map(|&(x, y)| (x - center.0, y - center.1))
        .collect();

    let dp1 = nbr_pts[0].0 * nbr_pts[1].0 + nbr_pts[0].1 * nbr_pts[1].1;
    if dp1.abs() < 1.0e-3 {
        let mut res = vec![(nbr_bids[0], nbr_bids[1])];
        let dp2 = nbr_pts[0].0 * nbr_pts[2].0 + nbr_pts[0].1 * nbr_pts[2].1;
        if dp2.abs() < 1.0e-3 {
            res.push((nbr_bids[0], nbr_bids[2]));
        } else {
            res.push((nbr_bids[0], nbr_bids[3]));
        }
        res
    } else {
        vec![(nbr_bids[0], nbr_bids[2]), (nbr_bids[0], nbr_bids[3])]
    }
}

// BEGIN RDKIT CPP FUNCTION: pickFirstRingToEmbed (DepictUtils.cpp)
// RDKit✔️✔️: int pickFirstRingToEmbed(const RDKit::ROMol &mol,
// RDKit✔️✔️:                          const RDKit::VECT_INT_VECT &fusedRings) {
// RDKit✔️✔️:   int res = -1;
// RDKit✔️✔️:   unsigned int maxSize = 0;
// RDKit✔️✔️:   int subs, minsubs = static_cast<int>(1e8);
// RDKit✔️✔️:   int cnt = 0;
// RDKit✔️✔️:   for (const auto &fusedRing : fusedRings) {
// RDKit✔️✔️:     subs = 0;
// RDKit✔️✔️:     for (auto rii : fusedRing) {
// RDKit✔️✔️:       if (mol.getAtomWithIdx(rii)->getDegree() > 2) {
// RDKit✔️✔️:         ++subs;
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:     if (subs < minsubs) {
// RDKit✔️✔️:       res = cnt;
// RDKit✔️✔️:       minsubs = subs;
// RDKit✔️✔️:       maxSize = fusedRing.size();
// RDKit✔️✔️:     } else if (subs == minsubs) {
// RDKit✔️✔️:       if (fusedRing.size() > maxSize) {
// RDKit✔️✔️:         res = cnt;
// RDKit✔️✔️:         maxSize = fusedRing.size();
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:     cnt++;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
fn rdkit_pick_first_ring_to_embed(degree: &[usize], fused_rings: &[Vec<usize>]) -> usize {
    let mut res = usize::MAX;
    let mut max_size = 0usize;
    let mut min_subs = 100_000_000usize;
    let mut cnt = 0usize;
    for fused_ring in fused_rings {
        let mut subs = 0usize;
        for &rii in fused_ring {
            if degree[rii] > 2 {
                subs += 1;
            }
        }
        if subs < min_subs {
            res = cnt;
            min_subs = subs;
            max_size = fused_ring.len();
        } else if subs == min_subs && fused_ring.len() > max_size {
            res = cnt;
            max_size = fused_ring.len();
        }
        cnt += 1;
    }
    res
}

// RDKit✔️✔️: RDKit::VECT_INT_VECT findCoreRings(const RDKit::VECT_INT_VECT &fusedRings,
// RDKit✔️✔️:                                    RDKit::INT_VECT &coreRingsIds,
// RDKit✔️✔️:                                    const RDKit::ROMol &mol) {
// RDKit✔️✔️:   // simplify the fused rings to a set of core rings by iteratively removing
// RDKit✔️✔️:   // rings that share only one or two consecutive atoms. These
// RDKit✔️✔️:   // are trivial to embed after the core rings have been embedded and will make
// RDKit✔️✔️:   // template matching more powerful since it will not be affected by the side
// RDKit✔️✔️:   // rings
// RDKit✔️✔️:
// RDKit✔️✔️:   boost::dynamic_bitset<> removedRings(fusedRings.size());
// RDKit✔️✔️:   bool removedARing = false;
// RDKit✔️✔️:   do {
// RDKit✔️✔️:     removedARing = false;
// RDKit✔️✔️:     for (unsigned int currRingId = 0; currRingId < fusedRings.size();
// RDKit✔️✔️:          currRingId++) {
// RDKit✔️✔️:       if (removedRings[currRingId] || removedARing) {
// RDKit✔️✔️:         continue;
// RDKit✔️✔️:       }
// RDKit✔️✔️:       auto nIntersectingAtoms = 0u;
// RDKit✔️✔️:       int aid1 = -1;
// RDKit✔️✔️:       int aid2 = -1;
// RDKit✔️✔️:       for (unsigned int otherRingId = 0; otherRingId < fusedRings.size();
// RDKit✔️✔️:            otherRingId++) {
// RDKit✔️✔️:         if (currRingId == otherRingId || removedRings[otherRingId]) {
// RDKit✔️✔️:           continue;
// RDKit✔️✔️:         }
// RDKit✔️✔️:         RDKit::INT_VECT commmonAtoms;
// RDKit✔️✔️:         RDKit::Intersect(fusedRings[currRingId], fusedRings[otherRingId],
// RDKit✔️✔️:                          commmonAtoms);
// RDKit✔️✔️:         for (auto rii : commmonAtoms) {
// RDKit✔️✔️:           if (rii != aid1 && rii != aid2) {
// RDKit✔️✔️:             ++nIntersectingAtoms;
// RDKit✔️✔️:             if (aid1 == -1) {
// RDKit✔️✔️:               aid1 = rii;
// RDKit✔️✔️:             } else {
// RDKit✔️✔️:               aid2 = rii;
// RDKit✔️✔️:             }
// RDKit✔️✔️:             if (nIntersectingAtoms == 2) {
// RDKit✔️✔️:               break;
// RDKit✔️✔️:             }
// RDKit✔️✔️:           }
// RDKit✔️✔️:         }
// RDKit✔️✔️:       }
// RDKit✔️✔️:       // note that the set of rings is not SSSR because we use symmetrizeSSSR,
// RDKit✔️✔️:       // so we cannot force a check for only one fused ring. Instead we make
// RDKit✔️✔️:       // sure that this ring shares only one atom or one bond (two consecutive
// RDKit✔️✔️:       // atoms)
// RDKit✔️✔️:       if (nIntersectingAtoms == 1 ||
// RDKit✔️✔️:           (nIntersectingAtoms == 2 &&
// RDKit✔️✔️:            mol.getBondBetweenAtoms(aid1, aid2) != nullptr)) {
// RDKit✔️✔️:         removedRings[currRingId] = true;
// RDKit✔️✔️:         removedARing = true;
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   } while (removedARing);
// RDKit✔️✔️:   RDKit::VECT_INT_VECT res;
// RDKit✔️✔️:   for (unsigned int currRingId = 0; currRingId < fusedRings.size();
// RDKit✔️✔️:        currRingId++) {
// RDKit✔️✔️:     if (!removedRings[currRingId]) {
// RDKit✔️✔️:       res.push_back(fusedRings[currRingId]);
// RDKit✔️✔️:       coreRingsIds.push_back(currRingId);
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
fn rdkit_find_core_rings(
    fused_rings: &[Vec<usize>],
    bonds: &[Bond],
) -> (Vec<Vec<usize>>, Vec<usize>) {
    let mut removed_rings = vec![false; fused_rings.len()];
    let mut removed_a_ring;
    loop {
        removed_a_ring = false;
        for curr_ring_id in 0..fused_rings.len() {
            if removed_rings[curr_ring_id] || removed_a_ring {
                continue;
            }
            let mut n_intersecting_atoms = 0usize;
            let mut aid1: Option<usize> = None;
            let mut aid2: Option<usize> = None;
            for other_ring_id in 0..fused_rings.len() {
                if curr_ring_id == other_ring_id || removed_rings[other_ring_id] {
                    continue;
                }
                for &rii in &fused_rings[curr_ring_id] {
                    if fused_rings[other_ring_id].contains(&rii)
                        && Some(rii) != aid1
                        && Some(rii) != aid2
                    {
                        n_intersecting_atoms += 1;
                        if aid1.is_none() {
                            aid1 = Some(rii);
                        } else {
                            aid2 = Some(rii);
                        }
                        if n_intersecting_atoms == 2 {
                            break;
                        }
                    }
                }
            }
            if n_intersecting_atoms == 1
                || (n_intersecting_atoms == 2
                    && aid1
                        .zip(aid2)
                        .and_then(|(a, b)| bond_between_idx(bonds, a, b))
                        .is_some())
            {
                removed_rings[curr_ring_id] = true;
                removed_a_ring = true;
            }
        }
        if !removed_a_ring {
            break;
        }
    }
    let mut res = Vec::new();
    let mut core_ring_ids = Vec::new();
    for (curr_ring_id, ring) in fused_rings.iter().enumerate() {
        if !removed_rings[curr_ring_id] {
            res.push(ring.clone());
            core_ring_ids.push(curr_ring_id);
        }
    }
    (res, core_ring_ids)
}

// RDKit✔️✔️: RDKit::INT_VECT findNextRingToEmbed(const RDKit::INT_VECT &doneRings,
// RDKit✔️✔️:                                     const RDKit::VECT_INT_VECT &fusedRings,
// RDKit✔️✔️:                                     int &nextId) {
// RDKit✔️✔️:   // REVIEW: We are changing this after Issue166
// RDKit✔️✔️:   // Originally the ring that have maximum number of atoms in common with the
// RDKit✔️✔️:   // atoms
// RDKit✔️✔️:   // that have already been embedded will be the ring that will get embedded.
// RDKit✔️✔️:   // But
// RDKit✔️✔️:   // if we can find a ring with two atoms in common with the embedded atoms, we
// RDKit✔️✔️:   // will choose that first before systems with more than 2 atoms in common.
// RDKit✔️✔️:   // Cases with two atoms in common are in general flat systems to start with
// RDKit✔️✔️:   // and can be embedded cleanly. when there are more than 2 atoms in common,
// RDKit✔️✔️:   // these are most likely bridged systems, which are screwed up anyway, might
// RDKit✔️✔️:   // as well screw them up later if we do not have a system with two rings in
// RDKit✔️✔️:   // common then we will return the ring with max, common atoms
// RDKit✔️✔️:   PRECONDITION(doneRings.size() > 0, "");
// RDKit✔️✔️:   PRECONDITION(fusedRings.size() > 1, "");
// RDKit✔️✔️:
// RDKit✔️✔️:   RDKit::INT_VECT commonAtoms, res, doneAtoms, notDone;
// RDKit✔️✔️:   for (int i = 0; i < rdcast<int>(fusedRings.size()); i++) {
// RDKit✔️✔️:     if (std::find(doneRings.begin(), doneRings.end(), i) == doneRings.end()) {
// RDKit✔️✔️:       notDone.push_back(i);
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   RDKit::Union(fusedRings, doneAtoms, &notDone);
// RDKit✔️✔️:
// RDKit✔️✔️:   int maxCommonAtoms = 0;
// RDKit✔️✔️:
// RDKit✔️✔️:   int currRingId = 0;
// RDKit✔️✔️:   for (const auto &fusedRing : fusedRings) {
// RDKit✔️✔️:     if (std::find(doneRings.begin(), doneRings.end(), currRingId) !=
// RDKit✔️✔️:         doneRings.end()) {
// RDKit✔️✔️:       currRingId++;
// RDKit✔️✔️:       continue;
// RDKit✔️✔️:     }
// RDKit✔️✔️:     commonAtoms.clear();
// RDKit✔️✔️:     int numCommonAtoms = 0;
// RDKit✔️✔️:     for (auto rii : fusedRing) {
// RDKit✔️✔️:       if (std::find(doneAtoms.begin(), doneAtoms.end(), (rii)) !=
// RDKit✔️✔️:           doneAtoms.end()) {
// RDKit✔️✔️:         commonAtoms.push_back(rii);
// RDKit✔️✔️:         numCommonAtoms++;
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:     if (numCommonAtoms == 2) {
// RDKit✔️✔️:       // if we found a ring with two atoms in common get out
// RDKit✔️✔️:       nextId = currRingId;
// RDKit✔️✔️:       return commonAtoms;  // FIX: this causes the rendering to be non-canonical
// RDKit✔️✔️:     }
// RDKit✔️✔️:     if (numCommonAtoms > maxCommonAtoms) {
// RDKit✔️✔️:       maxCommonAtoms = numCommonAtoms;
// RDKit✔️✔️:       nextId = currRingId;
// RDKit✔️✔️:       res = commonAtoms;
// RDKit✔️✔️:     }
// RDKit✔️✔️:     ++currRingId;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   // here is an additional constrain we will put on the common atoms it is quite
// RDKit✔️✔️:   // likely that the common atoms form a chain (it is possible we can construct
// RDKit✔️✔️:   // some weird cases where this does not hold true - but for now we will assume
// RDKit✔️✔️:   // this is true. However the IDs in the res may not be in the order of going
// RDKit✔️✔️:   // from one end of the chain to the other -
// RDKit✔️✔️:   //  here is an example C1CCC(CC12)CCC2
// RDKit✔️✔️:   // - two rings here with three atoms in common
// RDKit✔️✔️:   // let ring1:(0,1,2,3,4,5) be a ring that is already embedded, then let
// RDKit✔️✔️:   // ring2:(4,3,6,7,8,5) be the ring that we found to be the next ring we should
// RDKit✔️✔️:   // embed. The commonAtoms are (4,3,5) - note that they will be in this order
// RDKit✔️✔️:   // since the rings are always traversed in order. Now we would like these
// RDKit✔️✔️:   // common atoms to be returned in the order (5,4,3) - then we have a
// RDKit✔️✔️:   // continuous chain, we can do this by simply looking at the original ring
// RDKit✔️✔️:   // order (4,3,6,7,8,5) and observing that 5 need to come to the front
// RDKit✔️✔️:
// RDKit✔️✔️:   // find out how many atoms from the end we need to move to the front
// RDKit✔️✔️:   unsigned int cmnLst = 0;
// RDKit✔️✔️:   unsigned int nCmn = res.size();
// RDKit✔️✔️:   for (unsigned int i = 0; i < nCmn; i++) {
// RDKit✔️✔️:     if (res[i] == fusedRings[nextId][i]) {
// RDKit✔️✔️:       cmnLst++;
// RDKit✔️✔️:     } else {
// RDKit✔️✔️:       break;
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   // now do the moving if we have to
// RDKit✔️✔️:   if ((cmnLst > 0) && (cmnLst < res.size())) {
// RDKit✔️✔️:     RDKit::INT_VECT tempV = res;
// RDKit✔️✔️:
// RDKit✔️✔️:     for (unsigned int i = cmnLst; i < nCmn; i++) {
// RDKit✔️✔️:       res[i - cmnLst] = tempV[i];
// RDKit✔️✔️:     }
// RDKit✔️✔️:     unsigned int nMov = nCmn - cmnLst;
// RDKit✔️✔️:     for (unsigned int i = 0; i < cmnLst; i++) {
// RDKit✔️✔️:       res[nMov + i] = tempV[i];
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   POSTCONDITION(res.size() > 0, "");
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
fn rdkit_find_next_ring_to_embed(
    done_rings: &[usize],
    fused_rings: &[Vec<usize>],
) -> Option<(usize, Vec<usize>)> {
    let mut not_done = Vec::new();
    for i in 0..fused_rings.len() {
        if !done_rings.contains(&i) {
            not_done.push(i);
        }
    }

    let mut done_atoms = Vec::new();
    for (ring_id, ring) in fused_rings.iter().enumerate() {
        if not_done.contains(&ring_id) {
            continue;
        }
        for &atom_id in ring {
            if !done_atoms.contains(&atom_id) {
                done_atoms.push(atom_id);
            }
        }
    }

    let mut res = Vec::new();
    let mut next_id = usize::MAX;
    let mut max_common_atoms = 0usize;

    for (curr_ring_id, fused_ring) in fused_rings.iter().enumerate() {
        if done_rings.contains(&curr_ring_id) {
            continue;
        }
        let mut common_atoms = Vec::new();
        let mut num_common_atoms = 0usize;
        for &rii in fused_ring {
            if done_atoms.contains(&rii) {
                common_atoms.push(rii);
                num_common_atoms += 1;
            }
        }
        if num_common_atoms == 2 {
            return Some((curr_ring_id, common_atoms));
        }
        if num_common_atoms > max_common_atoms {
            max_common_atoms = num_common_atoms;
            next_id = curr_ring_id;
            res = common_atoms;
        }
    }

    if res.is_empty() || next_id == usize::MAX {
        return None;
    }

    let mut cmn_lst = 0usize;
    let n_cmn = res.len();
    for i in 0..n_cmn {
        if res[i] == fused_rings[next_id][i] {
            cmn_lst += 1;
        } else {
            break;
        }
    }
    if cmn_lst > 0 && cmn_lst < res.len() {
        let temp_v = res.clone();
        for i in cmn_lst..n_cmn {
            res[i - cmn_lst] = temp_v[i];
        }
        let n_mov = n_cmn - cmn_lst;
        for i in 0..cmn_lst {
            res[n_mov + i] = temp_v[i];
        }
    }

    Some((next_id, res))
}

// RDKit✔️✔️: RDKit::INT_VECT getAllRotatableBonds(const RDKit::ROMol &mol) {
// RDKit✔️✔️:   RDKit::INT_VECT res;
// RDKit✔️✔️:   for (const auto bond : mol.bonds()) {
// RDKit✔️✔️:     int bid = bond->getIdx();
// RDKit✔️✔️:     if ((bond->getStereo() <= RDKit::Bond::STEREOANY) &&
// RDKit✔️✔️:         (!(mol.getRingInfo()->numBondRings(bid)))) {
// RDKit✔️✔️:       res.push_back(bid);
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
fn rdkit_get_all_rotatable_bonds(
    bonds: &[Bond],
    is_ring_bond: impl Fn(usize) -> bool,
) -> Vec<usize> {
    let mut res = Vec::new();
    for bond in bonds {
        let bid = bond.id().index();
        if matches!(bond.stereo(), BondStereo::None | BondStereo::Any) && !is_ring_bond(bid) {
            res.push(bid);
        }
    }
    res
}

// RDKit✔️✔️: RDKit::INT_VECT getRotatableBonds(const RDKit::ROMol &mol, unsigned int aid1,
// RDKit✔️✔️:                                   unsigned int aid2) {
// RDKit✔️✔️:   PRECONDITION(aid1 < mol.getNumAtoms(), "");
// RDKit✔️✔️:   PRECONDITION(aid2 < mol.getNumAtoms(), "");
// RDKit✔️✔️:
// RDKit✔️✔️:   RDKit::INT_LIST path = RDKit::MolOps::getShortestPath(mol, aid1, aid2);
// RDKit✔️✔️:   RDKit::INT_VECT res;
// RDKit✔️✔️:   if (path.size() >= 4) {
// RDKit✔️✔️:     // remove the first atom (aid1) and last atom (aid2)
// RDKit✔️✔️:     CHECK_INVARIANT(static_cast<unsigned int>(path.front()) == aid1,
// RDKit✔️✔️:                     "bad first element");
// RDKit✔️✔️:     path.pop_front();
// RDKit✔️✔️:     CHECK_INVARIANT(static_cast<unsigned int>(path.back()) == aid2,
// RDKit✔️✔️:                     "bad last element");
// RDKit✔️✔️:     path.pop_back();
// RDKit✔️✔️:
// RDKit✔️✔️:     auto pid = path.front();
// RDKit✔️✔️:     for (auto aid : path) {
// RDKit✔️✔️:       if (aid == pid) {
// RDKit✔️✔️:         continue;
// RDKit✔️✔️:       }
// RDKit✔️✔️:       const RDKit::Bond *bond = mol.getBondBetweenAtoms(pid, aid);
// RDKit✔️✔️:       int bid = bond->getIdx();
// RDKit✔️✔️:       if ((bond->getStereo() <= RDKit::Bond::STEREOANY) &&
// RDKit✔️✔️:           (!(mol.getRingInfo()->numBondRings(bid)))) {
// RDKit✔️✔️:         res.push_back(bid);
// RDKit✔️✔️:       }
// RDKit✔️✔️:       pid = aid;
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
fn rdkit_get_rotatable_bonds_between(
    aid1: usize,
    aid2: usize,
    num_atoms: usize,
    bonds: &[Bond],
    adjacency: &[Vec<usize>],
    is_ring_bond: impl Fn(usize) -> bool,
) -> Vec<usize> {
    if aid1 >= num_atoms || aid2 >= num_atoms {
        return Vec::new();
    }

    let mut q = VecDeque::<usize>::new();
    let mut prev = BTreeMap::<usize, usize>::new();
    q.push_back(aid1);
    prev.insert(aid1, aid1);
    while let Some(u) = q.pop_front() {
        if u == aid2 {
            break;
        }
        for &v in &adjacency[u] {
            if let std::collections::btree_map::Entry::Vacant(e) = prev.entry(v) {
                e.insert(u);
                q.push_back(v);
            }
        }
    }
    if !prev.contains_key(&aid2) {
        return Vec::new();
    }

    let mut path = vec![aid2];
    let mut cur = aid2;
    while cur != aid1 {
        cur = prev[&cur];
        path.push(cur);
    }
    path.reverse();

    let mut res = Vec::new();
    if path.len() >= 4 {
        if path.first().copied() != Some(aid1) || path.last().copied() != Some(aid2) {
            return Vec::new();
        }
        path.remove(0);
        path.pop();
        let mut pid = path[0];
        for aid in path {
            if aid == pid {
                continue;
            }
            if let Some(bid) = bond_between_idx(bonds, pid, aid) {
                let bond = &bonds[bid];
                if matches!(bond.stereo(), BondStereo::None | BondStereo::Any) && !is_ring_bond(bid)
                {
                    res.push(bid);
                }
            }
            pid = aid;
        }
    }
    res
}

// RDKit✔️✔️: void _recurseAtomOneSide(unsigned int endAid, unsigned int begAid,
// RDKit✔️✔️:                          const RDKit::ROMol *mol, RDKit::INT_VECT &flipAids) { ... }
fn recurse_atom_one_side(
    end_aid: usize,
    beg_aid: usize,
    adjacency: &[Vec<usize>],
    flip_aids: &mut Vec<usize>,
) {
    flip_aids.push(end_aid);
    for &nbr in &adjacency[end_aid] {
        if nbr != beg_aid && !flip_aids.contains(&nbr) {
            recurse_atom_one_side(nbr, beg_aid, adjacency, flip_aids);
        }
    }
}

// RDKit✔️✔️: double _crossVal(const RDGeom::Point2D &v1, const RDGeom::Point2D &v2) { ... }
fn cross_val(v1: (f64, f64), v2: (f64, f64)) -> f64 {
    v1.0 * v2.1 - v2.0 * v1.1
}

// RDKit✔️✔️: PAIR_I_I _findClosestPair(unsigned int beg1, unsigned int end1,
// RDKit✔️✔️:                           unsigned int beg2, unsigned int end2,
// RDKit✔️✔️:                           const RDKit::ROMol &mol, const double *dmat) { ... }
fn find_closest_pair(
    beg1: usize,
    end1: usize,
    beg2: usize,
    end2: usize,
    num_atoms: usize,
    dmat: &[f64],
) -> (usize, usize) {
    let candidates = [
        (dmat[beg1 * num_atoms + beg2], (beg1, beg2)),
        (dmat[beg1 * num_atoms + end2], (beg1, end2)),
        (dmat[end1 * num_atoms + beg2], (end1, beg2)),
        (dmat[end1 * num_atoms + end2], (end1, end2)),
    ];
    let mut best = candidates[0];
    for candidate in candidates.iter().skip(1) {
        if candidate.0 < best.0 {
            best = *candidate;
        }
    }
    best.1
}

// RDKit✔️✔️: unsigned int _findDeg1Neighbor(const RDKit::ROMol *mol, unsigned int aid) { ... }
fn find_deg1_neighbor(adjacency: &[Vec<usize>], aid: usize) -> Option<usize> {
    (adjacency[aid].len() == 1).then_some(adjacency[aid][0])
}

// RDKit✔️✔️: unsigned int _findClosestNeighbor(const RDKit::ROMol *mol, const double *dmat,
// RDKit✔️✔️:                                   unsigned int aid1, unsigned int aid2) { ... }
fn find_closest_neighbor(
    adjacency: &[Vec<usize>],
    dmat: &[f64],
    num_atoms: usize,
    aid1: usize,
    aid2: usize,
) -> usize {
    let base = aid1 * num_atoms;
    let mut res = 0usize;
    let mut mdist = 1.0e8;
    for &nbr in &adjacency[aid2] {
        let d = dmat[base + nbr];
        if d < mdist {
            mdist = d;
            res = nbr;
        }
    }
    res
}

fn component_graph_distance_matrix(
    num_atoms: usize,
    comp: &[usize],
    adjacency: &[Vec<usize>],
) -> Vec<f64> {
    let comp_set: BTreeSet<usize> = comp.iter().copied().collect();
    let mut out = vec![f64::INFINITY; num_atoms * num_atoms];
    for &start in comp {
        let mut q = VecDeque::<usize>::new();
        out[start * num_atoms + start] = 0.0;
        q.push_back(start);
        while let Some(cur) = q.pop_front() {
            let next_dist = out[start * num_atoms + cur] + 1.0;
            for &nb in &adjacency[cur] {
                if !comp_set.contains(&nb) {
                    continue;
                }
                let slot = &mut out[start * num_atoms + nb];
                if slot.is_infinite() {
                    *slot = next_dist;
                    q.push_back(nb);
                }
            }
        }
    }
    out
}

fn shortest_path_in_component(
    src: usize,
    dst: usize,
    comp: &[usize],
    adjacency: &[Vec<usize>],
) -> Option<Vec<usize>> {
    let comp_set: BTreeSet<usize> = comp.iter().copied().collect();
    let mut q = VecDeque::<usize>::new();
    let mut prev = BTreeMap::<usize, usize>::new();
    q.push_back(src);
    prev.insert(src, src);
    while let Some(u) = q.pop_front() {
        if u == dst {
            break;
        }
        for &v in &adjacency[u] {
            if !comp_set.contains(&v) {
                continue;
            }
            if let std::collections::btree_map::Entry::Vacant(e) = prev.entry(v) {
                e.insert(u);
                q.push_back(v);
            }
        }
    }
    if !prev.contains_key(&dst) {
        return None;
    }
    let mut path = vec![dst];
    let mut cur = dst;
    while cur != src {
        cur = prev[&cur];
        path.push(cur);
    }
    path.reverse();
    Some(path)
}

fn is_ring_bond_in_component(
    bond_idx: usize,
    bonds: &[Bond],
    comp: &[usize],
    adjacency: &[Vec<usize>],
) -> bool {
    let comp_set: BTreeSet<usize> = comp.iter().copied().collect();
    let b = &bonds[bond_idx];
    let u = b.begin().index();
    let v = b.end().index();
    let mut stack = vec![u];
    let mut seen = BTreeSet::<usize>::new();
    seen.insert(u);
    while let Some(cur) = stack.pop() {
        for &nb in &adjacency[cur] {
            if !comp_set.contains(&nb) {
                continue;
            }
            if (cur == u && nb == v) || (cur == v && nb == u) {
                continue;
            }
            if nb == v {
                return true;
            }
            if seen.insert(nb) {
                stack.push(nb);
            }
        }
    }
    false
}

fn any_non_ring_bonds_on_path(
    aid: usize,
    path: &[usize],
    bonds: &[Bond],
    comp: &[usize],
    adjacency: &[Vec<usize>],
) -> usize {
    let mut prev = aid;
    let mut n_open = 0usize;
    for &pi in path {
        if let Some(bid) = bond_between_idx(bonds, prev, pi)
            && !is_ring_bond_in_component(bid, bonds, comp, adjacency)
        {
            n_open += 1;
        }
        prev = pi;
    }
    n_open
}

fn recurse_deg_two_ring_atoms_component(
    aid: usize,
    bonds: &[Bond],
    comp: &[usize],
    adjacency: &[Vec<usize>],
    r_path: &mut Vec<usize>,
    nbr_map: &mut BTreeMap<usize, Vec<usize>>,
) {
    let nbrs: Vec<usize> = adjacency[aid]
        .iter()
        .copied()
        .filter(|&nbr| {
            comp.contains(&nbr)
                && bond_between_idx(bonds, aid, nbr)
                    .is_some_and(|bid| is_ring_bond_in_component(bid, bonds, comp, adjacency))
        })
        .collect();
    if nbrs.len() != 2 {
        return;
    }
    r_path.push(aid);
    nbr_map.insert(aid, nbrs.clone());
    for nbr in nbrs {
        if !r_path.contains(&nbr) {
            recurse_deg_two_ring_atoms_component(nbr, bonds, comp, adjacency, r_path, nbr_map);
        }
    }
}

fn rdkit_query_contains_atomic_number(
    query: &QueryNode<AtomQueryPredicate>,
    atomic_number: u8,
) -> bool {
    match query {
        QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(n)) => *n == atomic_number,
        QueryNode::Predicate(_) => false,
        QueryNode::And(children) | QueryNode::Or(children) => children
            .iter()
            .any(|child| rdkit_query_contains_atomic_number(child, atomic_number)),
        QueryNode::Not(child) => rdkit_query_contains_atomic_number(child, atomic_number),
    }
}

fn rdkit_parse_former_nbr_indices_prop(value: &str) -> Vec<usize> {
    value
        .split(',')
        .filter(|part| !part.is_empty())
        .filter_map(|part| part.parse::<usize>().ok())
        .collect()
}

// RDKit❗✔️: bool isAtomTerminalRGroupOrQueryHydrogen(const Atom *atom) {
// RDKit❗✔️:   return (atom->getDegree() == 1 && isAtomDummy(atom)) ||
// RDKit❗✔️:          (atom->hasQuery() &&
// RDKit❗✔️:           describeQuery(atom).find("AtomAtomicNum 1 = val") !=
// RDKit❗✔️:               std::string::npos);
// RDKit❗✔️: }
fn rdkit_is_atom_terminal_rgroup_or_query_hydrogen(
    atom_idx: usize,
    atoms: &[Atom],
    adjacency: &crate::AdjacencyList,
) -> bool {
    let atom = &atoms[atom_idx];
    (adjacency.neighbors_of(atom_idx).len() == 1 && atom.atomic_number() == 0)
        || atom
            .query()
            .is_some_and(|query| rdkit_query_contains_atomic_number(query, 1))
}

// RDKit❗✔️: bool hasTerminalRGroupOrQueryHydrogen(const RDKit::ROMol &mol) {
// RDKit❗✔️:   // we do not need the allowRGroups logic if there are no
// RDKit❗✔️:   // terminal dummy atoms
// RDKit❗✔️:   auto atoms = mol.atoms();
// RDKit❗✔️:   return std::any_of(atoms.begin(), atoms.end(),
// RDKit❗✔️:                      RDKit::isAtomTerminalRGroupOrQueryHydrogen);
// RDKit❗✔️: }
fn rdkit_has_terminal_rgroup_or_query_hydrogen(mol: &Molecule) -> bool {
    let adjacency = &mol.topology_block().adjacency;
    mol.atoms().iter().any(|atom| {
        rdkit_is_atom_terminal_rgroup_or_query_hydrogen(atom.id().index(), mol.atoms(), adjacency)
    })
}

// RDKit❗✔️: The `adjustQueryProperties(...)` prelude is ported only for the
// RDKit flags used here: converting eligible single bonds into a
// `single|aromatic` query before terminal R-group reduction.
fn rdkit_adjust_query_properties_for_depiction_template(mol: &mut Molecule) {
    let adjacency = mol.topology_block().adjacency.clone();
    let atoms = mol.atoms().to_vec();
    for bond in &mut mol.topology_block_mut().bonds {
        let begin = bond.begin().index();
        let end = bond.end().index();
        if bond.query().is_some() || bond.order() != BondOrder::Single {
            continue;
        }
        let begin_is_aromatic = atoms[begin].is_aromatic();
        let end_is_aromatic = atoms[end].is_aromatic();
        let replace = if begin_is_aromatic ^ end_is_aromatic {
            (begin_is_aromatic && adjacency.neighbors_of(end).len() == 1)
                || (end_is_aromatic && adjacency.neighbors_of(begin).len() == 1)
        } else {
            begin_is_aromatic && end_is_aromatic
        };
        if replace {
            bond.set_query(Some(QueryNode::predicate(
                crate::BondQueryPredicate::OrderIn(vec![BondOrder::Single, BondOrder::Aromatic]),
            )));
        }
    }
}

// RDKit❗✔️: std::unique_ptr<RDKit::RWMol> prepareTemplateForRGroups(
// RDKit❗✔️:     RDKit::RWMol &templateMol) {
// RDKit❗✔️:   auto queryParams = RDKit::MolOps::AdjustQueryParameters::noAdjustments();
// RDKit❗✔️:   queryParams.adjustSingleBondsToDegreeOneNeighbors = true;
// RDKit❗✔️:   queryParams.adjustSingleBondsBetweenAromaticAtoms = true;
// RDKit❗✔️:   RDKit::MolOps::adjustQueryProperties(templateMol, &queryParams);
// RDKit❗✔️:   std::map<unsigned int, unsigned int> removedIdxToNbrIdx;
// RDKit❗✔️:   std::unique_ptr<RDKit::RWMol> reducedTemplateMol;
// RDKit❗✔️:   for (const auto &bond : templateMol.bonds()) {
// RDKit❗✔️:     int atomIdxToRemove = -1;
// RDKit❗✔️:     int nbrIdx = -1;
// RDKit❗✔️:     auto beginAtom = bond->getBeginAtom();
// RDKit❗✔️:     auto endAtom = bond->getEndAtom();
// RDKit❗✔️:     if (RDKit::isAtomTerminalRGroupOrQueryHydrogen(beginAtom) &&
// RDKit❗✔️:         endAtom->hasQuery()) {
// RDKit❗✔️:       atomIdxToRemove = beginAtom->getIdx();
// RDKit❗✔️:       nbrIdx = endAtom->getIdx();
// RDKit❗✔️:     } else if (RDKit::isAtomTerminalRGroupOrQueryHydrogen(endAtom) &&
// RDKit❗✔️:                beginAtom->hasQuery()) {
// RDKit❗✔️:       atomIdxToRemove = endAtom->getIdx();
// RDKit❗✔️:       nbrIdx = beginAtom->getIdx();
// RDKit❗✔️:     }
// RDKit❗✔️:     if (atomIdxToRemove != -1) {
// RDKit❗✔️:       removedIdxToNbrIdx[atomIdxToRemove] = nbrIdx;
// RDKit❗✔️:     }
// RDKit❗✔️:   }
// RDKit❗✔️:   if (!removedIdxToNbrIdx.empty()) {
// RDKit❗✔️:     reducedTemplateMol.reset(new RDKit::RWMol(templateMol));
// RDKit❗✔️:     for (auto reducedTemplateAtom : reducedTemplateMol->atoms()) {
// RDKit❗✔️:       auto formerIdx = reducedTemplateAtom->getIdx();
// RDKit❗✔️:       reducedTemplateAtom->setProp(FORMER_IDX, formerIdx);
// RDKit❗✔️:       auto it = removedIdxToNbrIdx.find(formerIdx);
// RDKit❗✔️:       if (it != removedIdxToNbrIdx.end()) {
// RDKit❗✔️:         auto otherAtom = reducedTemplateMol->getAtomWithIdx(it->second);
// RDKit❗✔️:         std::vector<unsigned int> formerNbrIndices;
// RDKit❗✔️:         otherAtom->getPropIfPresent(FORMER_NBR_INDICES, formerNbrIndices);
// RDKit❗✔️:         formerNbrIndices.push_back(formerIdx);
// RDKit❗✔️:         otherAtom->setProp(FORMER_NBR_INDICES, formerNbrIndices);
// RDKit❗✔️:       }
// RDKit❗✔️:     }
// RDKit❗✔️:     reducedTemplateMol->beginBatchEdit();
// RDKit❗✔️:     for (const auto &pair : removedIdxToNbrIdx) {
// RDKit❗✔️:       reducedTemplateMol->removeAtom(
// RDKit❗✔️:           reducedTemplateMol->getAtomWithIdx(pair.first));
// RDKit❗✔️:     }
// RDKit❗✔️:     reducedTemplateMol->commitBatchEdit();
// RDKit❗✔️:   }
// RDKit❗✔️:   return reducedTemplateMol;
// RDKit❗✔️: }
fn rdkit_prepare_template_for_rgroups(template_mol: &Molecule) -> Option<Molecule> {
    let mut prepared = template_mol.clone();
    rdkit_adjust_query_properties_for_depiction_template(&mut prepared);

    let atoms = prepared.atoms().to_vec();
    let adjacency = prepared.topology_block().adjacency.clone();
    let mut removed_idx_to_nbr_idx = BTreeMap::<usize, usize>::new();
    for bond in prepared.bonds() {
        let begin = bond.begin().index();
        let end = bond.end().index();
        let begin_atom = &atoms[begin];
        let end_atom = &atoms[end];
        if rdkit_is_atom_terminal_rgroup_or_query_hydrogen(begin, &atoms, &adjacency)
            && end_atom.query().is_some()
        {
            removed_idx_to_nbr_idx.insert(begin, end);
        } else if rdkit_is_atom_terminal_rgroup_or_query_hydrogen(end, &atoms, &adjacency)
            && begin_atom.query().is_some()
        {
            removed_idx_to_nbr_idx.insert(end, begin);
        }
    }
    if removed_idx_to_nbr_idx.is_empty() {
        return None;
    }

    let atom_count = prepared.num_atoms();
    let mut neighbor_props = vec![Vec::<usize>::new(); atom_count];
    for (&removed_idx, &nbr_idx) in &removed_idx_to_nbr_idx {
        neighbor_props[nbr_idx].push(removed_idx);
    }
    for atom in &mut prepared.topology_block_mut().atoms {
        let former_idx = atom.id().index();
        atom.set_prop(RDKIT_FORMER_IDX_PROP, former_idx.to_string());
        if !neighbor_props[former_idx].is_empty() {
            let encoded = neighbor_props[former_idx]
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",");
            atom.set_prop(RDKIT_FORMER_NBR_INDICES_PROP, encoded);
        }
    }
    let atoms_to_remove = removed_idx_to_nbr_idx
        .keys()
        .copied()
        .map(crate::AtomId::new)
        .collect::<Vec<_>>();
    prepared
        .topology_block_mut()
        .remove_atoms_with_mapping(&atoms_to_remove);
    Some(prepared)
}

// RDKit❗✔️: void reducedToFullMatches(const RDKit::RWMol &reducedQuery,
// RDKit❗✔️:                           const RDKit::RWMol &molHs,
// RDKit❗✔️:                           std::vector<RDKit::MatchVectType> &matches) {
// RDKit❗✔️:   boost::dynamic_bitset<> molHsMatches(molHs.getNumAtoms());
// RDKit❗✔️:   for (auto &match : matches) {
// RDKit❗✔️:     molHsMatches.reset();
// RDKit❗✔️:     for (const auto &pair : match) {
// RDKit❗✔️:       molHsMatches.set(pair.second);
// RDKit❗✔️:     }
// RDKit❗✔️:     RDKit::MatchVectType newMatch;
// RDKit❗✔️:     for (auto pairIt = match.begin(); pairIt != match.end(); ++pairIt) {
// RDKit❗✔️:       const auto reducedQueryAtom = reducedQuery.getAtomWithIdx(pairIt->first);
// RDKit❗✔️:       const auto molAtom = molHs.getAtomWithIdx(pairIt->second);
// RDKit❗✔️:       unsigned int formerIdx;
// RDKit❗✔️:       reducedQueryAtom->getProp(FORMER_IDX, formerIdx);
// RDKit❗✔️:       pairIt->first = formerIdx;
// RDKit❗✔️:       std::vector<unsigned int> formerNbrIndices;
// RDKit❗✔️:       reducedQueryAtom->getPropIfPresent(FORMER_NBR_INDICES, formerNbrIndices);
// RDKit❗✔️:       for (const auto &molNbr : molHs.atomNeighbors(molAtom)) {
// RDKit❗✔️:         if (formerNbrIndices.empty()) {
// RDKit❗✔️:           break;
// RDKit❗✔️:         }
// RDKit❗✔️:         auto molNbrIdx = molNbr->getIdx();
// RDKit❗✔️:         if (!molHsMatches.test(molNbrIdx)) {
// RDKit❗✔️:           auto formerNbrIdx = formerNbrIndices.back();
// RDKit❗✔️:           formerNbrIndices.pop_back();
// RDKit❗✔️:           newMatch.emplace_back(formerNbrIdx, molNbrIdx);
// RDKit❗✔️:         }
// RDKit❗✔️:       }
// RDKit❗✔️:     }
// RDKit❗✔️:     auto matchSize = match.size();
// RDKit❗✔️:     match.resize(matchSize + newMatch.size());
// RDKit❗✔️:     std::move(newMatch.begin(), newMatch.end(), match.begin() + matchSize);
// RDKit❗✔️:   }
// RDKit❗✔️: }
fn rdkit_reduced_to_full_matches(
    reduced_query: &Molecule,
    mol_hs: &Molecule,
    matches: &mut Vec<Vec<(usize, usize)>>,
) {
    let adjacency = &mol_hs.topology_block().adjacency;
    let mut mol_hs_matches = vec![false; mol_hs.num_atoms()];
    for match_vec in matches.iter_mut() {
        mol_hs_matches.fill(false);
        for &(_, mol_idx) in match_vec.iter() {
            if let Some(slot) = mol_hs_matches.get_mut(mol_idx) {
                *slot = true;
            }
        }
        let mut new_match = Vec::new();
        for pair in match_vec.iter_mut() {
            let reduced_query_atom = &reduced_query.atoms()[pair.0];
            let former_idx = reduced_query_atom
                .prop(RDKIT_FORMER_IDX_PROP)
                .and_then(|value| value.parse::<usize>().ok())
                .unwrap_or(pair.0);
            let mol_atom_idx = pair.1;
            pair.0 = former_idx;
            let mut former_nbr_indices = reduced_query_atom
                .prop(RDKIT_FORMER_NBR_INDICES_PROP)
                .map(rdkit_parse_former_nbr_indices_prop)
                .unwrap_or_default();
            for mol_nbr in adjacency.neighbors_of(mol_atom_idx) {
                if former_nbr_indices.is_empty() {
                    break;
                }
                let mol_nbr_idx = mol_nbr.atom_index;
                if !mol_hs_matches[mol_nbr_idx] {
                    let former_nbr_idx = former_nbr_indices.pop().expect("checked non-empty");
                    new_match.push((former_nbr_idx, mol_nbr_idx));
                }
            }
        }
        match_vec.extend(new_match);
    }
}

fn rdkit_select_3d_conformer_index(
    mol: &Molecule,
    conf_id: isize,
) -> Result<usize, Coordinate2DError> {
    let conformers = mol.conformers_3d();
    if conformers.is_empty() {
        return Err(Coordinate2DError::InvalidInput(
            "constrained depiction requires an available 3D conformer",
        ));
    }
    if conf_id == -1 {
        return Ok(0);
    }
    let requested = usize::try_from(conf_id).map_err(|_| {
        Coordinate2DError::InvalidInput("constrained depiction conformer id must be >= -1")
    })?;
    conformers
        .iter()
        .position(|conformer| conformer.id() == requested)
        .ok_or(Coordinate2DError::InvalidInput(
            "constrained depiction conformer id is out of range",
        ))
}

// RDKit✔️✔️: bool isAtomTerminalRGroupOrQueryHydrogen(const Atom *atom) {
// RDKit✔️✔️:   return (atom->getDegree() == 1 && isAtomDummy(atom)) ||
// RDKit✔️✔️:          (atom->hasQuery() &&
// RDKit✔️✔️:           describeQuery(atom).find("AtomAtomicNum 1 = val") !=
// RDKit✔️✔️:               std::string::npos);
// RDKit✔️✔️: }
// RDKit✔️❌: const MatchVectType &getMostSubstitutedCoreMatch(
// RDKit✔️❌:     const ROMol &mol, const ROMol &core,
// RDKit✔️❌:     const std::vector<MatchVectType> &matches) {
// RDKit✔️❌:   detail::ScoreMatchesByDegreeOfCoreSubstitution matchScorer(mol, core,
// RDKit✔️❌:                                                              matches);
// RDKit✔️❌:   return matchScorer.getMostSubstitutedCoreMatch();
// RDKit✔️❌: }
// RDKit✔️❌: std::vector<MatchVectType> sortMatchesByDegreeOfCoreSubstitution(
// RDKit✔️❌:     const ROMol &mol, const ROMol &core,
// RDKit✔️❌:     const std::vector<MatchVectType> &matches) {
// RDKit✔️❌:   detail::ScoreMatchesByDegreeOfCoreSubstitution matchScorer(mol, core,
// RDKit✔️❌:                                                              matches);
// RDKit✔️❌:   return matchScorer.sortMatchesByDegreeOfCoreSubstitution();
// RDKit✔️❌: }
fn rdkit_score_match_by_degree_of_core_substitution(
    mol: &Molecule,
    query: &Molecule,
    match_vec: &[(usize, usize)],
) -> f64 {
    let na = mol.num_atoms();
    let sum_indices = (na * (na + 1) / 2) as f64;
    let mut penalty = 0.0;
    let mut i = 0.0;
    for &(query_idx, mol_idx) in match_vec {
        i += mol_idx as f64;
        let query_atom = &query.atoms()[query_idx];
        let mol_atom = &mol.atoms()[mol_idx];
        if mol_atom.atomic_number() == 1
            && rdkit_is_atom_terminal_rgroup_or_query_hydrogen(
                query_idx,
                query.atoms(),
                &query.topology_block().adjacency,
            )
        {
            penalty += 1.0;
        }
    }
    penalty + i / sum_indices
}

fn rdkit_get_most_substituted_core_match(
    mol: &Molecule,
    query: &Molecule,
    matches: &[Vec<(usize, usize)>],
) -> Option<Vec<(usize, usize)>> {
    matches
        .iter()
        .min_by(|left, right| {
            rdkit_score_match_by_degree_of_core_substitution(mol, query, left).total_cmp(
                &rdkit_score_match_by_degree_of_core_substitution(mol, query, right),
            )
        })
        .cloned()
}

fn rdkit_sort_matches_by_degree_of_core_substitution(
    mol: &Molecule,
    query: &Molecule,
    matches: &[Vec<(usize, usize)>],
) -> Vec<Vec<(usize, usize)>> {
    let mut indexed = matches
        .iter()
        .cloned()
        .map(|match_vec| {
            let score = rdkit_score_match_by_degree_of_core_substitution(mol, query, &match_vec);
            (score, match_vec)
        })
        .collect::<Vec<_>>();
    indexed.sort_by(|left, right| left.0.total_cmp(&right.0));
    indexed
        .into_iter()
        .map(|(_, match_vec)| match_vec)
        .collect()
}

// RDKit❗✔️: void invertMolBlockWedgingInfo(ROMol &mol) {
// RDKit❗✔️:   for (auto b : mol.bonds()) {
// RDKit❗✔️:     int bond_dir = -1;
// RDKit❗✔️:     if (b->getPropIfPresent<int>(common_properties::_MolFileBondStereo,
// RDKit❗✔️:                                  bond_dir)) {
// RDKit❗✔️:       if (bond_dir == 1) {
// RDKit❗✔️:         b->setProp<int>(common_properties::_MolFileBondStereo, 6);
// RDKit❗✔️:       } else if (bond_dir == 6) {
// RDKit❗✔️:         b->setProp<int>(common_properties::_MolFileBondStereo, 1);
// RDKit❗✔️:       }
// RDKit❗✔️:     }
// RDKit❗✔️:     int cfg = -1;
// RDKit❗✔️:     if (b->getPropIfPresent<int>(common_properties::_MolFileBondCfg, cfg)) {
// RDKit❗✔️:       if (cfg == 1) {
// RDKit❗✔️:         b->setProp<int>(common_properties::_MolFileBondCfg, 3);
// RDKit❗✔️:       } else if (cfg == 3) {
// RDKit❗✔️:         b->setProp<int>(common_properties::_MolFileBondCfg, 1);
// RDKit❗✔️:       }
// RDKit❗✔️:     }
// RDKit❗✔️:   }
// RDKit❗✔️: }
fn rdkit_invert_molblock_wedging_info(mol: &mut Molecule) {
    for bond in &mut mol.topology_block_mut().bonds {
        if let Some(bond_dir) = bond.prop("_MolFileBondStereo").map(str::to_string) {
            match bond_dir.as_str() {
                "1" => bond.set_prop("_MolFileBondStereo", "6"),
                "6" => bond.set_prop("_MolFileBondStereo", "1"),
                _ => {}
            }
        }
        if let Some(cfg) = bond.prop("_MolFileBondCfg").map(str::to_string) {
            match cfg.as_str() {
                "1" => bond.set_prop("_MolFileBondCfg", "3"),
                "3" => bond.set_prop("_MolFileBondCfg", "1"),
                _ => {}
            }
        }
    }
}

// RDKit✔️✔️: void clearMolBlockWedgingInfo(ROMol &mol) {
// RDKit✔️✔️:   for (auto b : mol.bonds()) {
// RDKit✔️✔️:     b->clearProp(common_properties::_MolFileBondStereo);
// RDKit✔️✔️:     b->clearProp(common_properties::_MolFileBondCfg);
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
fn rdkit_clear_molblock_wedging_info(mol: &mut Molecule) {
    for bond in &mut mol.topology_block_mut().bonds {
        bond.clear_prop("_MolFileBondStereo");
        bond.clear_prop("_MolFileBondCfg");
    }
}

fn rdkit_select_2d_conformer_index(
    mol: &Molecule,
    conf_id: isize,
) -> Result<usize, Coordinate2DError> {
    let conformers = mol.conformers_2d();
    if conformers.is_empty() {
        return Err(Coordinate2DError::InvalidInput(
            "constrained depiction requires an available 2D conformer",
        ));
    }
    if conf_id == -1 {
        return Ok(0);
    }
    let requested = usize::try_from(conf_id).map_err(|_| {
        Coordinate2DError::InvalidInput("constrained depiction conformer id must be >= -1")
    })?;
    conformers
        .iter()
        .position(|conformer| conformer.id() == requested)
        .ok_or(Coordinate2DError::InvalidInput(
            "constrained depiction conformer id is out of range",
        ))
}

// RDKit✔️❌: void removeAllConformersButOne(RDKit::ROMol &mol, int confId) {
// RDKit✔️❌:   std::vector<int> conformerIndicesToRemove;
// RDKit✔️❌:   for (auto confIt = mol.beginConformers(); confIt != mol.endConformers();
// RDKit✔️❌:        ++confIt) {
// RDKit✔️❌:     int i = (*confIt)->getId();
// RDKit✔️❌:     if ((confId != -1 && i == confId) ||
// RDKit✔️❌:         (confId == -1 && confIt == mol.beginConformers())) {
// RDKit✔️❌:       continue;
// RDKit✔️❌:     }
// RDKit✔️❌:     conformerIndicesToRemove.push_back(i);
// RDKit✔️❌:   }
// RDKit✔️❌:   for (auto i : conformerIndicesToRemove) {
// RDKit✔️❌:     mol.removeConformer(i);
// RDKit✔️❌:   }
// RDKit✔️❌:   CHECK_INVARIANT(mol.getNumConformers() == 1, "");
// RDKit✔️❌:   mol.getConformer().setId(0);
// RDKit✔️❌: }
fn rdkit_remove_all_2d_conformers_but_one(
    mol: &mut Molecule,
    conf_id: isize,
) -> Result<(), Coordinate2DError> {
    let keep_index = rdkit_select_2d_conformer_index(mol, conf_id)?;
    let survivor = {
        let conformers = &mol.coordinate_block().conformers_2d;
        conformers
            .get(keep_index)
            .cloned()
            .ok_or(Coordinate2DError::InvalidInput(
                "constrained depiction conformer id is out of range",
            ))?
    };
    let coordinate_block = mol.coordinate_block_mut();
    coordinate_block.conformers_2d.clear();
    coordinate_block.conformers_2d.push(survivor.with_id(0));
    Ok(())
}

const RDKIT_ALIGN_POINTS_TOLERANCE: f64 = 1.0e-6;
const RDKIT_ALIGN_POINTS_MAX_ITERATIONS: usize = 50;

fn rdkit_transform3d_identity() -> [[f64; 4]; 4] {
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
}

fn rdkit_transform3d_mul(lhs: &[[f64; 4]; 4], rhs: &[[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let mut out = [[0.0; 4]; 4];
    for row in 0..4 {
        for col in 0..4 {
            out[row][col] = (0..4).map(|k| lhs[row][k] * rhs[k][col]).sum();
        }
    }
    out
}

fn rdkit_transform3d_transform_point(trans: &[[f64; 4]; 4], point: [f64; 3]) -> [f64; 3] {
    [
        trans[0][0] * point[0] + trans[0][1] * point[1] + trans[0][2] * point[2] + trans[0][3],
        trans[1][0] * point[0] + trans[1][1] * point[1] + trans[1][2] * point[2] + trans[1][3],
        trans[2][0] * point[0] + trans[2][1] * point[1] + trans[2][2] * point[2] + trans[2][3],
    ]
}

fn rdkit_transform3d_set_translation(trans: &mut [[f64; 4]; 4], move_vec: [f64; 3]) {
    trans[0][3] = move_vec[0];
    trans[1][3] = move_vec[1];
    trans[2][3] = move_vec[2];
}

fn rdkit_transform3d_set_rotation_from_quaternion(trans: &mut [[f64; 4]; 4], quaternion: [f64; 4]) {
    let q0 = quaternion[0];
    let q1 = quaternion[1];
    let q2 = quaternion[2];
    let q3 = quaternion[3];
    let n = q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3;
    let s = if n > 0.0 { 2.0 / n } else { 0.0 };
    let x = q1 * s;
    let y = q2 * s;
    let z = q3 * s;
    let wx = q0 * x;
    let wy = q0 * y;
    let wz = q0 * z;
    let xx = q1 * x;
    let xy = q1 * y;
    let xz = q1 * z;
    let yy = q2 * y;
    let yz = q2 * z;
    let zz = q3 * z;

    *trans = rdkit_transform3d_identity();
    trans[0][0] = 1.0 - (yy + zz);
    trans[0][1] = xy - wz;
    trans[0][2] = xz + wy;
    trans[1][0] = xy + wz;
    trans[1][1] = 1.0 - (xx + zz);
    trans[1][2] = yz - wx;
    trans[2][0] = xz - wy;
    trans[2][1] = yz + wx;
    trans[2][2] = 1.0 - (xx + yy);
}

fn rdkit_transform3d_reflect(trans: &mut [[f64; 4]; 4]) {
    for row in trans.iter_mut().take(3) {
        for cell in row.iter_mut().take(3) {
            *cell = -*cell;
        }
    }
}

fn rdkit_weighted_sum_of_points(points: &[[f64; 3]]) -> [f64; 3] {
    points.iter().fold([0.0, 0.0, 0.0], |mut acc, point| {
        acc[0] += point[0];
        acc[1] += point[1];
        acc[2] += point[2];
        acc
    })
}

fn rdkit_weighted_sum_of_len_sq(points: &[[f64; 3]]) -> f64 {
    points.iter().fold(0.0, |acc, point| {
        acc + point[0] * point[0] + point[1] * point[1] + point[2] * point[2]
    })
}

fn rdkit_compute_covariance_mat(
    ref_points: &[[f64; 3]],
    probe_points: &[[f64; 3]],
) -> [[f64; 3]; 3] {
    let mut cov_mat = [[0.0; 3]; 3];
    for (rpt, ppt) in ref_points.iter().zip(probe_points.iter()) {
        cov_mat[0][0] += ppt[0] * rpt[0];
        cov_mat[0][1] += ppt[0] * rpt[1];
        cov_mat[0][2] += ppt[0] * rpt[2];
        cov_mat[1][0] += ppt[1] * rpt[0];
        cov_mat[1][1] += ppt[1] * rpt[1];
        cov_mat[1][2] += ppt[1] * rpt[2];
        cov_mat[2][0] += ppt[2] * rpt[0];
        cov_mat[2][1] += ppt[2] * rpt[1];
        cov_mat[2][2] += ppt[2] * rpt[2];
    }
    cov_mat
}

fn rdkit_reflect_covariance_mat(cov_mat: &mut [[f64; 3]; 3]) {
    for row in cov_mat.iter_mut() {
        for cell in row.iter_mut() {
            *cell = -*cell;
        }
    }
}

fn rdkit_convert_cov_mat_to_quad(
    cov_mat: &[[f64; 3]; 3],
    rpt_sum: [f64; 3],
    ppt_sum: [f64; 3],
    wts_sum: f64,
) -> [[f64; 4]; 4] {
    let px_rx = cov_mat[0][0] - (ppt_sum[0] / wts_sum) * rpt_sum[0];
    let px_ry = cov_mat[0][1] - (ppt_sum[0] / wts_sum) * rpt_sum[1];
    let px_rz = cov_mat[0][2] - (ppt_sum[0] / wts_sum) * rpt_sum[2];
    let py_rx = cov_mat[1][0] - (ppt_sum[1] / wts_sum) * rpt_sum[0];
    let py_ry = cov_mat[1][1] - (ppt_sum[1] / wts_sum) * rpt_sum[1];
    let py_rz = cov_mat[1][2] - (ppt_sum[1] / wts_sum) * rpt_sum[2];
    let pz_rx = cov_mat[2][0] - (ppt_sum[2] / wts_sum) * rpt_sum[0];
    let pz_ry = cov_mat[2][1] - (ppt_sum[2] / wts_sum) * rpt_sum[1];
    let pz_rz = cov_mat[2][2] - (ppt_sum[2] / wts_sum) * rpt_sum[2];

    let mut quad = [[0.0; 4]; 4];
    quad[0][0] = -2.0 * (px_rx + py_ry + pz_rz);
    quad[1][1] = -2.0 * (px_rx - py_ry - pz_rz);
    quad[2][2] = -2.0 * (py_ry - pz_rz - px_rx);
    quad[3][3] = -2.0 * (pz_rz - px_rx - py_ry);
    quad[0][1] = 2.0 * (py_rz - pz_ry);
    quad[1][0] = quad[0][1];
    quad[0][2] = 2.0 * (pz_rx - px_rz);
    quad[2][0] = quad[0][2];
    quad[0][3] = 2.0 * (px_ry - py_rx);
    quad[3][0] = quad[0][3];
    quad[1][2] = -2.0 * (px_ry + py_rx);
    quad[2][1] = quad[1][2];
    quad[1][3] = -2.0 * (pz_rx + px_rz);
    quad[3][1] = quad[1][3];
    quad[2][3] = -2.0 * (py_rz + pz_ry);
    quad[3][2] = quad[2][3];
    quad
}

fn rdkit_align_points_jacobi(
    mut quad: [[f64; 4]; 4],
    max_iter: usize,
) -> ([f64; 4], [[f64; 4]; 4]) {
    let mut eigen_vecs = [[0.0; 4]; 4];
    let mut eigen_vals = [0.0; 4];
    for j in 0..4 {
        eigen_vecs[j][j] = 1.0;
        eigen_vals[j] = quad[j][j];
    }

    for _ in 0..max_iter {
        let mut diag_norm = 0.0;
        let mut off_diag_norm = 0.0;
        for j in 0..4 {
            diag_norm += eigen_vals[j].abs();
            for row in quad.iter().take(j) {
                off_diag_norm += row[j].abs();
            }
        }
        if diag_norm.abs() > 1.0e-16 && (off_diag_norm / diag_norm) <= RDKIT_ALIGN_POINTS_TOLERANCE
        {
            break;
        }
        for j in 1..4 {
            for i in 0..j {
                let b = quad[i][j];
                if b.abs() <= 0.0 {
                    continue;
                }
                let dma = eigen_vals[j] - eigen_vals[i];
                let t = if (dma.abs() + b.abs()) <= dma.abs() {
                    b / dma
                } else {
                    let q = 0.5 * dma / b;
                    let mut t = 1.0 / (q.abs() + rdkit_sqrt(1.0 + q * q));
                    if q < 0.0 {
                        t = -t;
                    }
                    t
                };
                let c = 1.0 / rdkit_sqrt(t * t + 1.0);
                let s = t * c;
                quad[i][j] = 0.0;
                for k in 0..i {
                    let atemp = c * quad[k][i] - s * quad[k][j];
                    quad[k][j] = s * quad[k][i] + c * quad[k][j];
                    quad[k][i] = atemp;
                }
                for k in (i + 1)..j {
                    let atemp = c * quad[i][k] - s * quad[k][j];
                    quad[k][j] = s * quad[i][k] + c * quad[k][j];
                    quad[i][k] = atemp;
                }
                for k in (j + 1)..4 {
                    let atemp = c * quad[i][k] - s * quad[j][k];
                    quad[j][k] = s * quad[i][k] + c * quad[j][k];
                    quad[i][k] = atemp;
                }
                for row in &mut eigen_vecs {
                    let vtemp = c * row[i] - s * row[j];
                    row[j] = s * row[i] + c * row[j];
                    row[i] = vtemp;
                }
                let dtemp = c * c * eigen_vals[i] + s * s * eigen_vals[j] - 2.0 * c * s * b;
                eigen_vals[j] = s * s * eigen_vals[i] + c * c * eigen_vals[j] + 2.0 * c * s * b;
                eigen_vals[i] = dtemp;
            }
        }
    }

    for j in 0..3 {
        let mut k = j;
        let mut dtemp = eigen_vals[k];
        for (i, val) in eigen_vals.iter().enumerate().skip(j + 1) {
            if *val < dtemp {
                k = i;
                dtemp = eigen_vals[k];
            }
        }
        if k > j {
            eigen_vals[k] = eigen_vals[j];
            eigen_vals[j] = dtemp;
            for row in &mut eigen_vecs {
                row.swap(k, j);
            }
        }
    }

    (eigen_vals, eigen_vecs)
}

fn rdkit_align_points(
    ref_points: &[[f64; 3]],
    probe_points: &[[f64; 3]],
    reflect: bool,
    max_iterations: usize,
) -> Result<(f64, [[f64; 4]; 4]), Coordinate2DError> {
    if ref_points.len() != probe_points.len() {
        return Err(Coordinate2DError::InvalidInput(
            "alignment requires matching point counts",
        ));
    }
    let npt = ref_points.len();
    if npt == 0 {
        return Err(Coordinate2DError::InvalidInput(
            "alignment requires at least one point",
        ));
    }

    let mut trans = rdkit_transform3d_identity();
    let wts_sum = npt as f64;
    let mut rpt_sum = rdkit_weighted_sum_of_points(ref_points);
    let ppt_sum = rdkit_weighted_sum_of_points(probe_points);
    let rpt_sum_len_sq = rdkit_weighted_sum_of_len_sq(ref_points);
    let ppt_sum_len_sq = rdkit_weighted_sum_of_len_sq(probe_points);

    let mut cov_mat = rdkit_compute_covariance_mat(ref_points, probe_points);
    if reflect {
        rpt_sum = [-rpt_sum[0], -rpt_sum[1], -rpt_sum[2]];
        rdkit_reflect_covariance_mat(&mut cov_mat);
    }
    let quad = rdkit_convert_cov_mat_to_quad(&cov_mat, rpt_sum, ppt_sum, wts_sum);
    let (eigen_vals, eigen_vecs) = rdkit_align_points_jacobi(quad, max_iterations);
    let quaternion = [
        eigen_vecs[0][0],
        eigen_vecs[1][0],
        eigen_vecs[2][0],
        eigen_vecs[3][0],
    ];
    rdkit_transform3d_set_rotation_from_quaternion(&mut trans, quaternion);
    if reflect {
        rdkit_transform3d_reflect(&mut trans);
    }
    let mut ssr = eigen_vals[0]
        - ((ppt_sum[0] * ppt_sum[0] + ppt_sum[1] * ppt_sum[1] + ppt_sum[2] * ppt_sum[2])
            + (rpt_sum[0] * rpt_sum[0] + rpt_sum[1] * rpt_sum[1] + rpt_sum[2] * rpt_sum[2]))
            / wts_sum
        + rpt_sum_len_sq
        + ppt_sum_len_sq;
    if ssr < 0.0 && ssr.abs() < RDKIT_ALIGN_POINTS_TOLERANCE {
        ssr = 0.0;
    }
    if reflect {
        rpt_sum = [-rpt_sum[0], -rpt_sum[1], -rpt_sum[2]];
    }
    let transformed_ppt_sum = rdkit_transform3d_transform_point(&trans, ppt_sum);
    let move_vec = [
        (rpt_sum[0] - transformed_ppt_sum[0]) / wts_sum,
        (rpt_sum[1] - transformed_ppt_sum[1]) / wts_sum,
        (rpt_sum[2] - transformed_ppt_sum[2]) / wts_sum,
    ];
    rdkit_transform3d_set_translation(&mut trans, move_vec);
    Ok((ssr, trans))
}

fn rdkit_apply_transform_to_2d_conformer(
    mol: &mut Molecule,
    conf_id: isize,
    trans: &[[f64; 4]; 4],
) -> Result<(), Coordinate2DError> {
    let conf_index = rdkit_select_2d_conformer_index(mol, conf_id)?;
    let coords = mol.coordinate_block_mut().conformers_2d[conf_index].coordinates_mut();
    for point in coords.iter_mut() {
        let transformed = rdkit_transform3d_transform_point(trans, [point[0], point[1], 0.0]);
        point[0] = transformed[0];
        point[1] = transformed[1];
    }
    Ok(())
}

fn rdkit_get_alignment_transform(
    probe_mol: &Molecule,
    reference_mol: &Molecule,
    probe_conf_id: isize,
    reference_conf_id: isize,
    atom_map: &[(usize, usize)],
) -> Result<(f64, [[f64; 4]; 4]), Coordinate2DError> {
    let probe_conf_index = rdkit_select_2d_conformer_index(probe_mol, probe_conf_id)?;
    let reference_conf_index = rdkit_select_2d_conformer_index(reference_mol, reference_conf_id)?;
    let probe_coords = probe_mol.conformers_2d()[probe_conf_index].coordinates();
    let reference_coords = reference_mol.conformers_2d()[reference_conf_index].coordinates();

    let mut ref_points = Vec::with_capacity(atom_map.len());
    let mut probe_points = Vec::with_capacity(atom_map.len());
    for &(probe_atom_idx, reference_atom_idx) in atom_map {
        let probe = probe_coords
            .get(probe_atom_idx)
            .ok_or(Coordinate2DError::InvalidInput(
                "probe atom index in alignment atom map is out of range",
            ))?;
        let reference =
            reference_coords
                .get(reference_atom_idx)
                .ok_or(Coordinate2DError::InvalidInput(
                    "reference atom index in alignment atom map is out of range",
                ))?;
        probe_points.push([probe[0], probe[1], 0.0]);
        ref_points.push([reference[0], reference[1], 0.0]);
    }
    let (ssr, trans) = rdkit_align_points(
        &ref_points,
        &probe_points,
        false,
        RDKIT_ALIGN_POINTS_MAX_ITERATIONS,
    )?;
    Ok((rdkit_sqrt(ssr / probe_points.len() as f64), trans))
}

fn rdkit_get_best_alignment_transform(
    probe_mol: &Molecule,
    reference_mol: &Molecule,
    probe_conf_id: isize,
    reference_conf_id: isize,
    matches: &[Vec<(usize, usize)>],
    max_matches: usize,
) -> Result<([[f64; 4]; 4], Vec<(usize, usize)>), Coordinate2DError> {
    let mut best_rmsd = f64::INFINITY;
    let mut best_trans = rdkit_transform3d_identity();
    let mut best_match = None;
    for match_vec in matches.iter().take(max_matches) {
        let (rmsd, trans) = rdkit_get_alignment_transform(
            probe_mol,
            reference_mol,
            probe_conf_id,
            reference_conf_id,
            match_vec,
        )?;
        if rmsd < best_rmsd {
            best_rmsd = rmsd;
            best_trans = trans;
            best_match = Some(match_vec.clone());
        }
    }
    best_match
        .map(|match_vec| (best_trans, match_vec))
        .ok_or(Coordinate2DError::InvalidInput(
            "alignment requires at least one substructure match",
        ))
}

fn rdkit_substruct_matches_unordered(mol: &Molecule, query: &Molecule) -> Vec<Vec<(usize, usize)>> {
    let params = crate::SubstructMatchParams {
        max_matches: 1000,
        uniquify: false,
        use_chirality: false,
        specified_stereo_query_matches_unspecified: false,
    };
    crate::get_substruct_matches_with_params(mol, query, &params)
        .into_iter()
        .map(|match_result| {
            match_result
                .atom_mapping
                .into_iter()
                .enumerate()
                .collect::<Vec<_>>()
        })
        .collect()
}

fn rdkit_add_hs_copy(mol: &Molecule) -> Result<Molecule, Coordinate2DError> {
    mol.clone().with_hydrogens().map_err(|error| match error {
        crate::OperationError::InvalidInput { .. }
        | crate::OperationError::Chemistry { .. }
        | crate::OperationError::Unsupported { .. } => Coordinate2DError::UnsupportedFeature(
            format!("RDKit addHs equivalent failed during constrained depiction: {error}"),
        ),
        crate::OperationError::UnsupportedFeature { source, .. } => {
            Coordinate2DError::UnsupportedFeature(format!(
                "RDKit addHs equivalent failed during constrained depiction: {source}"
            ))
        }
        other => Coordinate2DError::UnsupportedFeature(format!(
            "RDKit addHs equivalent failed during constrained depiction: {other}"
        )),
    })
}

fn rdkit_append_2d_conformer(
    mol: &mut Molecule,
    coords: Vec<[f64; 2]>,
) -> Result<usize, Coordinate2DError> {
    if coords.len() != mol.num_atoms() {
        return Err(Coordinate2DError::InvalidInput(
            "2D conformer row count must match atom count",
        ));
    }
    let coordinate_block = mol.coordinate_block_mut();
    let next_id = coordinate_block
        .conformers_2d
        .iter()
        .map(crate::Conformer2D::id)
        .max()
        .map_or(0, |max_id| max_id + 1);
    coordinate_block
        .conformers_2d
        .push(crate::Conformer2D::new(next_id, coords));
    Ok(next_id)
}

fn rdkit_set_single_2d_conformer(
    mol: &mut Molecule,
    coords: Vec<[f64; 2]>,
) -> Result<usize, Coordinate2DError> {
    if coords.len() != mol.num_atoms() {
        return Err(Coordinate2DError::InvalidInput(
            "2D conformer row count must match atom count",
        ));
    }
    let coordinate_block = mol.coordinate_block_mut();
    coordinate_block.conformers_2d.clear();
    coordinate_block
        .conformers_2d
        .push(crate::Conformer2D::new(0, coords));
    Ok(0)
}

fn rdkit_compute_2d_coords_for_mol(
    mol: &mut Molecule,
    coord_map: Option<&BTreeMap<usize, [f64; 2]>>,
    canon_orient: bool,
    clear_confs: bool,
    force_rdkit: bool,
    use_ring_templates: bool,
) -> Result<usize, Coordinate2DError> {
    let coords = compute_2d_coords_with_options(
        mol.atoms(),
        mol.bonds(),
        coord_map,
        canon_orient,
        true,
        0,
        0,
        0,
        false,
        force_rdkit,
        use_ring_templates,
    )?;
    if clear_confs {
        rdkit_set_single_2d_conformer(mol, coords)
    } else {
        rdkit_append_2d_conformer(mol, coords)
    }
}

// RDKit✔️❌: void generateDepictionMatching2DStructure(
// RDKit✔️❌:     RDKit::ROMol &mol, const RDKit::ROMol &reference,
// RDKit✔️❌:     const RDKit::MatchVectType &refMatchVect, int confId,
// RDKit✔️❌:     const ConstrainedDepictionParams &p) {
// RDKit✔️❌:   if (refMatchVect.size() > reference.getNumAtoms()) {
// RDKit✔️❌:     throw DepictException(
// RDKit✔️❌:         "When a refMatchVect is provided, it must have size "
// RDKit✔️❌:         "<= number of atoms in the reference");
// RDKit✔️❌:   }
// RDKit✔️❌:   for (const auto &mv : refMatchVect) {
// RDKit✔️❌:     if (mv.first > static_cast<int>(reference.getNumAtoms())) {
// RDKit✔️❌:       throw DepictException(
// RDKit✔️❌:           "Reference atom index in refMatchVect out of range");
// RDKit✔️❌:     }
// RDKit✔️❌:     if (mv.second > static_cast<int>(mol.getNumAtoms())) {
// RDKit✔️❌:       throw DepictException("Molecule atom index in refMatchVect out of range");
// RDKit✔️❌:     }
// RDKit✔️❌:   }
// RDKit✔️❌:   bool hasExistingCoords = mol.getNumConformers() > 0;
// RDKit✔️❌:   bool shouldClearWedgingInfo = p.adjustMolBlockWedging && !hasExistingCoords;
// RDKit✔️❌:   bool shouldInvertWedgingIfRequired = false;
// RDKit✔️❌:   RDGeom::Transform3D trans;
// RDKit✔️❌:   if (p.alignOnly) { ... } else { ... }
// RDKit✔️❌:   if (shouldClearWedgingInfo) {
// RDKit✔️❌:     RDKit::Chirality::clearMolBlockWedgingInfo(mol);
// RDKit✔️❌:   } else if (shouldInvertWedgingIfRequired) {
// RDKit✔️❌:     invertWedgingIfMolHasFlipped(mol, trans);
// RDKit✔️❌:   }
// RDKit✔️❌: }
pub(crate) fn generate_depiction_matching_2d_structure_with_ref_match(
    mol: &mut Molecule,
    reference: &Molecule,
    ref_match_vect: &[(usize, usize)],
    conf_id: isize,
    params: &ConstrainedDepictionParams,
) -> Result<(), Coordinate2DError> {
    if ref_match_vect.len() > reference.num_atoms() {
        return Err(Coordinate2DError::InvalidInput(
            "When a refMatchVect is provided, it must have size <= number of atoms in the reference",
        ));
    }
    for &(reference_atom_idx, mol_atom_idx) in ref_match_vect {
        if reference_atom_idx >= reference.num_atoms() {
            return Err(Coordinate2DError::InvalidInput(
                "Reference atom index in refMatchVect out of range",
            ));
        }
        if mol_atom_idx >= mol.num_atoms() {
            return Err(Coordinate2DError::InvalidInput(
                "Molecule atom index in refMatchVect out of range",
            ));
        }
    }
    let has_existing_coords = !mol.conformers_2d().is_empty();
    let mut should_clear_wedging_info = params.adjust_molblock_wedging && !has_existing_coords;
    let mut should_invert_wedging_if_required = false;
    let mut trans = rdkit_transform3d_identity();

    if params.align_only {
        if !has_existing_coords {
            let _ = rdkit_compute_2d_coords_for_mol(
                mol,
                None,
                false,
                true,
                params.force_rdkit,
                params.use_ring_templates,
            )?;
        }
        let atom_map = ref_match_vect
            .iter()
            .map(|&(reference_atom_idx, mol_atom_idx)| (mol_atom_idx, reference_atom_idx))
            .collect::<Vec<_>>();
        let (_, alignment_trans) = rdkit_get_alignment_transform(
            mol,
            reference,
            params.existing_conf_id,
            conf_id,
            &atom_map,
        )?;
        trans = alignment_trans;
        rdkit_apply_transform_to_2d_conformer(mol, params.existing_conf_id, &trans)?;
        rdkit_remove_all_2d_conformers_but_one(mol, params.existing_conf_id)?;
        if !should_clear_wedging_info {
            should_invert_wedging_if_required = params.adjust_molblock_wedging;
        }
    } else {
        let reference_conf_index = rdkit_select_2d_conformer_index(reference, conf_id)?;
        let reference_coords = reference.conformers_2d()[reference_conf_index].coordinates();
        let mut coord_map = BTreeMap::new();
        for &(reference_atom_idx, mol_atom_idx) in ref_match_vect {
            let point = reference_coords[reference_atom_idx];
            coord_map.insert(mol_atom_idx, point);
        }
        let new_conf_id = rdkit_compute_2d_coords_for_mol(
            mol,
            Some(&coord_map),
            false,
            !(params.adjust_molblock_wedging && has_existing_coords),
            params.force_rdkit,
            params.use_ring_templates,
        )?;
        if params.adjust_molblock_wedging {
            const RMSD_THRESHOLD: f64 = 1.0e-2;
            const MSD_THRESHOLD: f64 = RMSD_THRESHOLD * RMSD_THRESHOLD;
            if !should_clear_wedging_info {
                let mut mol_matching_indices = vec![false; mol.num_atoms()];
                for &(_, mol_atom_idx) in ref_match_vect {
                    mol_matching_indices[mol_atom_idx] = true;
                }
                should_clear_wedging_info = mol.bonds().iter().any(|bond| {
                    (bond.prop("_MolFileBondStereo").is_some()
                        || bond.prop("_MolFileBondCfg").is_some())
                        && (!mol_matching_indices[bond.begin().index()]
                            || !mol_matching_indices[bond.end().index()])
                });
            }
            if !should_clear_wedging_info {
                let new_conf_index = rdkit_select_2d_conformer_index(mol, new_conf_id as isize)?;
                let mol_pos = mol.conformers_2d()[new_conf_index].coordinates();
                should_clear_wedging_info =
                    ref_match_vect
                        .iter()
                        .any(|&(reference_atom_idx, mol_atom_idx)| {
                            let dx =
                                mol_pos[mol_atom_idx][0] - reference_coords[reference_atom_idx][0];
                            let dy =
                                mol_pos[mol_atom_idx][1] - reference_coords[reference_atom_idx][1];
                            dx * dx + dy * dy > MSD_THRESHOLD
                        });
            }
            if !should_clear_wedging_info {
                let identity_match = ref_match_vect
                    .iter()
                    .map(|&(_, mol_atom_idx)| (mol_atom_idx, mol_atom_idx))
                    .collect::<Vec<_>>();
                let (rmsd, alignment_trans) = rdkit_get_alignment_transform(
                    mol,
                    mol,
                    new_conf_id as isize,
                    params.existing_conf_id,
                    &identity_match,
                )?;
                trans = alignment_trans;
                if rmsd > RMSD_THRESHOLD {
                    should_clear_wedging_info = true;
                } else {
                    should_invert_wedging_if_required = true;
                }
            }
        }
        if has_existing_coords {
            rdkit_remove_all_2d_conformers_but_one(mol, new_conf_id as isize)?;
        }
    }
    if should_clear_wedging_info {
        rdkit_clear_molblock_wedging_info(mol);
    } else if should_invert_wedging_if_required {
        rdkit_invert_wedging_if_mol_has_flipped(mol, &trans);
    }
    Ok(())
}

// RDKit✔️❌: void generateDepictionMatching2DStructure(
// RDKit✔️❌:     RDKit::ROMol &mol, const RDKit::ROMol &reference,
// RDKit✔️❌:     const RDKit::MatchVectType &refMatchVect, int confId, bool forceRDKit) {
// RDKit✔️❌:   ConstrainedDepictionParams p;
// RDKit✔️❌:   p.forceRDKit = forceRDKit;
// RDKit✔️❌:   generateDepictionMatching2DStructure(mol, reference, refMatchVect, confId, p);
// RDKit✔️❌: }
pub(crate) fn generate_depiction_matching_2d_structure_with_ref_match_force_rdkit(
    mol: &mut Molecule,
    reference: &Molecule,
    ref_match_vect: &[(usize, usize)],
    conf_id: isize,
    force_rdkit: bool,
) -> Result<(), Coordinate2DError> {
    let params = ConstrainedDepictionParams {
        force_rdkit,
        ..ConstrainedDepictionParams::default()
    };
    generate_depiction_matching_2d_structure_with_ref_match(
        mol,
        reference,
        ref_match_vect,
        conf_id,
        &params,
    )
}

// RDKit✔️❌: RDKit::MatchVectType generateDepictionMatching2DStructure(
// RDKit✔️❌:     RDKit::ROMol &mol, const RDKit::ROMol &reference, int confId,
// RDKit✔️❌:     const RDKit::ROMol *referencePattern,
// RDKit✔️❌:     const ConstrainedDepictionParams &params) { ... }
pub(crate) fn generate_depiction_matching_2d_structure(
    mol: &mut Molecule,
    reference: &Molecule,
    conf_id: isize,
    reference_pattern: Option<&Molecule>,
    params: &ConstrainedDepictionParams,
) -> Result<Vec<(usize, usize)>, Coordinate2DError> {
    let mut reference_hs = None;
    let mut mol_hs = None;
    let mut query_adj = None;
    let mut match_vect = Vec::<(usize, usize)>::new();
    let mut pattern_to_ref_matches = Vec::<Vec<(usize, usize)>>::new();
    let mut pattern_to_ref_match = Vec::<(usize, usize)>::new();
    let query = reference_pattern.unwrap_or(reference);
    let mut pattern_to_ref_mapping = vec![-1isize; query.num_atoms()];
    let mut p = *params;
    p.allow_rgroups = p.allow_rgroups && rdkit_has_terminal_rgroup_or_query_hydrogen(query);

    let mut reduced_query = None;
    let mut prb_mol = mol.clone();
    let mut ref_mol = query.clone();
    if p.allow_rgroups {
        let mol_hs_value = rdkit_add_hs_copy(mol)?;
        let mut query_adj_value = query.clone();
        reduced_query = rdkit_prepare_template_for_rgroups(&query_adj_value);
        prb_mol = mol_hs_value.clone();
        ref_mol = reduced_query
            .clone()
            .unwrap_or_else(|| query_adj_value.clone());
        mol_hs = Some(mol_hs_value);
        query_adj = Some(query_adj_value);
    }

    if let Some(reference_pattern) = reference_pattern {
        if p.allow_rgroups {
            let reference_hs_value = rdkit_add_hs_copy(reference)?;
            let query_adj_ref = query_adj.as_ref().expect("allow_rgroups implies query_adj");
            pattern_to_ref_matches =
                rdkit_substruct_matches_unordered(&reference_hs_value, &ref_mol);
            if let Some(reduced_query) = reduced_query.as_ref() {
                rdkit_reduced_to_full_matches(
                    reduced_query,
                    &reference_hs_value,
                    &mut pattern_to_ref_matches,
                );
            }
            if !pattern_to_ref_matches.is_empty() {
                pattern_to_ref_match = rdkit_get_most_substituted_core_match(
                    &reference_hs_value,
                    query_adj_ref,
                    &pattern_to_ref_matches,
                )
                .unwrap_or_default();
            }
            reference_hs = Some(reference_hs_value);
        } else if let Some(match_result) = crate::get_substruct_match(reference, reference_pattern)
        {
            pattern_to_ref_match = match_result.atom_mapping.into_iter().enumerate().collect();
        }
        if pattern_to_ref_match.is_empty() {
            return Err(Coordinate2DError::InvalidInput(
                "Reference pattern does not map to reference.",
            ));
        }
        let num_ref_atoms = reference.num_atoms();
        for &(pattern_idx, reference_idx) in &pattern_to_ref_match {
            if p.allow_rgroups && reference_idx >= num_ref_atoms {
                continue;
            }
            if pattern_idx >= pattern_to_ref_mapping.len() {
                return Err(Coordinate2DError::InvalidInput(
                    "reference pattern atom index out of range while building constrained depiction mapping",
                ));
            }
            pattern_to_ref_mapping[pattern_idx] = reference_idx as isize;
        }
    } else {
        for (idx, slot) in pattern_to_ref_mapping.iter_mut().enumerate() {
            *slot = idx as isize;
        }
    }

    if p.align_only {
        let mut matches = rdkit_substruct_matches_unordered(&prb_mol, &ref_mol);
        if !matches.is_empty() {
            if p.allow_rgroups {
                if let Some(reduced_query) = reduced_query.as_ref() {
                    rdkit_reduced_to_full_matches(
                        reduced_query,
                        mol_hs.as_ref().expect("allow_rgroups implies mol_hs"),
                        &mut matches,
                    );
                }
                let query_adj_ref = query_adj.as_ref().expect("allow_rgroups implies query_adj");
                matches = rdkit_sort_matches_by_degree_of_core_substitution(
                    &prb_mol,
                    query_adj_ref,
                    &matches,
                );
                let mut max_matched_heavies = -1isize;
                let mut max_pruned_match_size = -1isize;
                let mut pruned_matches = Vec::new();
                let num_mol_atoms = mol.num_atoms();
                for match_vec in matches {
                    let mut n_matched_heavies = 0isize;
                    let mut pruned_match = Vec::new();
                    for &(pattern_idx, mol_idx) in &match_vec {
                        let ref_atom = &query_adj_ref.atoms()[pattern_idx];
                        if rdkit_is_atom_terminal_rgroup_or_query_hydrogen(
                            pattern_idx,
                            query_adj_ref.atoms(),
                            &query_adj_ref.topology_block().adjacency,
                        ) {
                            if mol_idx >= num_mol_atoms {
                                continue;
                            }
                            n_matched_heavies += 1;
                        }
                        let ref_idx = pattern_to_ref_mapping[pattern_idx];
                        if ref_idx == -1 {
                            continue;
                        }
                        pruned_match.push((mol_idx, ref_idx as usize));
                    }
                    if n_matched_heavies < max_matched_heavies {
                        break;
                    }
                    max_matched_heavies = n_matched_heavies;
                    let pruned_match_size = pruned_match.len() as isize;
                    if pruned_match_size > max_pruned_match_size {
                        max_pruned_match_size = pruned_match_size;
                        pruned_matches.clear();
                    }
                    if pruned_match_size == max_pruned_match_size {
                        pruned_matches.push(pruned_match);
                    }
                }
                matches = pruned_matches;
            } else {
                for match_vec in &mut matches {
                    for pair in match_vec.iter_mut() {
                        let ref_idx = pattern_to_ref_mapping[pair.0];
                        pair.0 = pair.1;
                        pair.1 = ref_idx as usize;
                    }
                }
            }
            if mol.conformers_2d().is_empty() {
                let _ = rdkit_compute_2d_coords_for_mol(
                    mol,
                    None,
                    false,
                    true,
                    p.force_rdkit,
                    p.use_ring_templates,
                )?;
                if p.adjust_molblock_wedging {
                    rdkit_clear_molblock_wedging_info(mol);
                    p.adjust_molblock_wedging = false;
                }
            }
            const MAX_MATCHES: usize = 1000;
            let (trans, best_match) = rdkit_get_best_alignment_transform(
                mol,
                reference,
                p.existing_conf_id,
                conf_id,
                &matches,
                MAX_MATCHES,
            )?;
            match_vect = best_match
                .into_iter()
                .map(|(mol_idx, reference_idx)| (reference_idx, mol_idx))
                .collect();
            rdkit_apply_transform_to_2d_conformer(mol, p.existing_conf_id, &trans)?;
            rdkit_remove_all_2d_conformers_but_one(mol, p.existing_conf_id)?;
            if p.adjust_molblock_wedging {
                rdkit_invert_wedging_if_mol_has_flipped(mol, &trans);
            }
        }
    } else {
        if p.allow_rgroups {
            let mut matches = rdkit_substruct_matches_unordered(&prb_mol, &ref_mol);
            if !matches.is_empty() {
                if let Some(reduced_query) = reduced_query.as_ref() {
                    rdkit_reduced_to_full_matches(
                        reduced_query,
                        mol_hs.as_ref().expect("allow_rgroups implies mol_hs"),
                        &mut matches,
                    );
                }
                let query_adj_ref = query_adj.as_ref().expect("allow_rgroups implies query_adj");
                let num_mol_atoms = mol.num_atoms();
                if let Some(best_match) =
                    rdkit_get_most_substituted_core_match(&prb_mol, query_adj_ref, &matches)
                {
                    for (pattern_idx, mol_idx) in best_match {
                        if mol_idx < num_mol_atoms && pattern_to_ref_mapping[pattern_idx] != -1 {
                            match_vect.push((pattern_idx, mol_idx));
                        }
                    }
                }
            }
        } else if let Some(match_result) = crate::get_substruct_match(&prb_mol, &ref_mol) {
            match_vect = match_result.atom_mapping.into_iter().enumerate().collect();
        }
        if !match_vect.is_empty() {
            for pair in &mut match_vect {
                pair.0 = pattern_to_ref_mapping[pair.0] as usize;
            }
            generate_depiction_matching_2d_structure_with_ref_match(
                mol,
                reference,
                &match_vect,
                conf_id,
                &p,
            )?;
        }
    }

    if match_vect.is_empty() {
        if p.accept_failure {
            let _ = rdkit_compute_2d_coords_for_mol(
                mol,
                None,
                false,
                true,
                p.force_rdkit,
                p.use_ring_templates,
            )?;
            if p.adjust_molblock_wedging {
                rdkit_clear_molblock_wedging_info(mol);
            }
        } else {
            return Err(Coordinate2DError::InvalidInput(
                "Substructure match with reference not found.",
            ));
        }
    }
    Ok(match_vect)
}

// RDKit✔️❌: RDKit::MatchVectType generateDepictionMatching2DStructure(
// RDKit✔️❌:     RDKit::ROMol &mol, const RDKit::ROMol &reference, int confId,
// RDKit✔️❌:     const RDKit::ROMol *referencePattern, bool acceptFailure, bool forceRDKit,
// RDKit✔️❌:     bool allowOptionalAttachments) { ... }
pub(crate) fn generate_depiction_matching_2d_structure_simple(
    mol: &mut Molecule,
    reference: &Molecule,
    conf_id: isize,
    reference_pattern: Option<&Molecule>,
    accept_failure: bool,
    force_rdkit: bool,
    allow_optional_attachments: bool,
) -> Result<Vec<(usize, usize)>, Coordinate2DError> {
    let params = ConstrainedDepictionParams {
        accept_failure,
        force_rdkit,
        allow_rgroups: allow_optional_attachments,
        ..ConstrainedDepictionParams::default()
    };
    generate_depiction_matching_2d_structure(mol, reference, conf_id, reference_pattern, &params)
}

// RDKit✔️❌: void generateDepictionMatching3DStructure(RDKit::ROMol &mol,
// RDKit✔️❌:                                           const RDKit::ROMol &reference,
// RDKit✔️❌:                                           int confId,
// RDKit✔️❌:                                           RDKit::ROMol *referencePattern,
// RDKit✔️❌:                                           bool acceptFailure, bool forceRDKit) { ... }
pub(crate) fn generate_depiction_matching_3d_structure(
    mol: &mut Molecule,
    reference: &Molecule,
    conf_id: isize,
    reference_pattern: Option<&Molecule>,
    accept_failure: bool,
    force_rdkit: bool,
) -> Result<(), Coordinate2DError> {
    let num_ats = mol.num_atoms();
    if reference_pattern.is_none() && reference.num_atoms() < num_ats {
        if accept_failure {
            let _ = rdkit_compute_2d_coords_for_mol(mol, None, false, true, force_rdkit, false)?;
            return Ok(());
        }
        return Err(Coordinate2DError::InvalidInput(
            "Reference molecule not compatible with target molecule.",
        ));
    }

    let mut mol_to_ref = vec![-1isize; num_ats];
    if let Some(reference_pattern) = reference_pattern.filter(|pattern| pattern.num_atoms() > 0) {
        let mol_match_vect = crate::get_substruct_match(mol, reference_pattern);
        let ref_match_vect = crate::get_substruct_match(reference, reference_pattern);
        let (Some(mol_match_vect), Some(ref_match_vect)) = (mol_match_vect, ref_match_vect) else {
            if accept_failure {
                let _ =
                    rdkit_compute_2d_coords_for_mol(mol, None, false, true, force_rdkit, false)?;
                return Ok(());
            }
            return Err(Coordinate2DError::InvalidInput(
                "Reference pattern didn't match molecule or reference.",
            ));
        };
        for i in 0..mol_match_vect.atom_mapping.len() {
            mol_to_ref[mol_match_vect.atom_mapping[i]] = ref_match_vect.atom_mapping[i] as isize;
        }
    } else {
        for (idx, slot) in mol_to_ref.iter_mut().enumerate() {
            *slot = idx as isize;
        }
    }

    let reference_conf_index = rdkit_select_3d_conformer_index(reference, conf_id)?;
    let conf = reference.conformers_3d()[reference_conf_index].coordinates();
    let mut dmat = vec![-1.0; num_ats * (num_ats - 1) / 2];
    for i in 0..num_ats {
        if mol_to_ref[i] == -1 {
            continue;
        }
        let cds_i = conf[i];
        for j in (i + 1)..num_ats {
            if mol_to_ref[j] == -1 {
                continue;
            }
            let cds_j = conf[mol_to_ref[j] as usize];
            let dx = cds_i[0] - cds_j[0];
            let dy = cds_i[1] - cds_j[1];
            let dz = cds_i[2] - cds_j[2];
            dmat[(j * (j - 1) / 2) + i] = rdkit_sqrt(dx * dx + dy * dy + dz * dz);
        }
    }
    let coords = compute_2d_coords_mimic_distmat_with_params(
        mol.atoms(),
        mol.bonds(),
        Some(&dmat),
        &Compute2DCoordsMimicDistMatParameters {
            canon_orient: false,
            clear_confs: true,
            weight_dist_mat: 0.5,
            n_flips_per_sample: 3,
            n_samples: 100,
            sample_seed: 25,
            permute_deg4_nodes: true,
            force_rdkit,
        },
    )?;
    let _ = rdkit_set_single_2d_conformer(mol, coords)?;
    Ok(())
}

// RDKit❗✔️: bool invertWedgingIfMolHasFlipped(RDKit::ROMol &mol,
// RDKit❗✔️:                                   const RDGeom::Transform3D &trans) {
// RDKit❗✔️:   constexpr double FLIP_THRESHOLD = -0.99;
// RDKit❗✔️:   auto zRot = trans.getVal(2, 2);
// RDKit❗✔️:   bool shouldFlip = zRot < FLIP_THRESHOLD;
// RDKit❗✔️:   if (shouldFlip) {
// RDKit❗✔️:     RDKit::Chirality::invertMolBlockWedgingInfo(mol);
// RDKit❗✔️:   }
// RDKit❗✔️:   return shouldFlip;
// RDKit❗✔️: }
fn rdkit_invert_wedging_if_mol_has_flipped(mol: &mut Molecule, trans: &[[f64; 4]; 4]) -> bool {
    const FLIP_THRESHOLD: f64 = -0.99;
    let z_rot = trans[2][2];
    let should_flip = z_rot < FLIP_THRESHOLD;
    if should_flip {
        rdkit_invert_molblock_wedging_info(mol);
    }
    should_flip
}

// ── Main orchestrator ────────────────────────────────────────────────────────

// BEGIN RDKIT CPP FUNCTION: computeInitialCoords (RDDepictor.cpp)
// RDKit✔️❗: ported with cip_ranks parameter and explicit fragment return.
fn rdkit_compute_initial_efrags_strict(
    atoms: &[Atom],
    bonds: &[Bond],
    cip_ranks: &[u32],
    coord_map: Option<&BTreeMap<usize, [f64; 2]>>,
    use_ring_templates: bool,
) -> Result<Vec<RdkitEmbeddedFrag>, Coordinate2DError> {
    let n = atoms.len();
    let mut adjacency = vec![Vec::<usize>::new(); n];
    let mut degree = vec![0usize; n];
    for b in bonds {
        adjacency[b.begin().index()].push(b.end().index());
        adjacency[b.end().index()].push(b.begin().index());
        degree[b.begin().index()] += 1;
        degree[b.end().index()] += 1;
    }
    let atom_ranks: Vec<i32> = (0..n)
        .map(|i| atom_depict_rank(atoms[i].atomic_number(), degree[i]) as i32)
        .collect();
    let stereo_mol = build_rdkit_molecule_from_slices(atoms, bonds)?;
    // BEGIN RDKIT CPP FUNCTION: computeInitialCoords ring finding (RDDepictor.cpp)
    // RDKit✔️✔️:   bool includeDativeBonds = true;
    // RDKit✔️✔️:   RDKit::MolOps::symmetrizeSSSR(mol, arings, includeDativeBonds);
    let ring_info =
        crate::rings::symmetrize_sssr_with_options(&stereo_mol, true, false).map_err(|error| {
            Coordinate2DError::UnsupportedFeature(format!(
                "RDKit symmetrizeSSSR equivalent failed during 2D coordinate generation: {error}"
            ))
        })?;
    let arings: Vec<Vec<usize>> = ring_info
        .atom_rings()
        .iter()
        .map(|ring| ring.iter().map(|aid| aid.index()).collect())
        .collect();
    let atom_ring_counts: Vec<usize> = (0..n)
        .map(|aid| ring_info.num_atom_rings(AtomId::new(aid)))
        .collect();
    let ring_bond_ids: BTreeSet<usize> = ring_info
        .bond_rings()
        .iter()
        .flat_map(|ring| ring.iter().map(|bid| bid.index()))
        .collect();

    let mut efrags = Vec::new();
    let coord_map_tuples = coord_map.map(|coord_map| {
        coord_map
            .iter()
            .map(|(&aid, &coord)| (aid, (coord[0], coord[1])))
            .collect::<BTreeMap<usize, (f64, f64)>>()
    });

    let pre_spec = coord_map_tuples
        .as_ref()
        .is_some_and(|coord_map| coord_map.len() > 1);
    if let Some(coord_map) = coord_map_tuples
        .as_ref()
        .filter(|coord_map| coord_map.len() > 1)
    {
        efrags.push(RdkitEmbeddedFrag::from_coord_map(
            coord_map,
            atoms,
            bonds,
            &adjacency,
            &degree,
            cip_ranks,
            &atom_ring_counts,
        ));
    }

    if !arings.is_empty() {
        embed_fused_systems(
            atoms,
            bonds,
            &adjacency,
            &degree,
            cip_ranks,
            &atom_ring_counts,
            &arings,
            &mut efrags,
            coord_map_tuples.as_ref(),
            use_ring_templates,
        );
    }

    for seed in embed_nontetrahedral_stereo(atoms, bonds, &stereo_mol, &atom_ranks) {
        let coord_map = seed.into_iter().collect::<BTreeMap<usize, (f64, f64)>>();
        efrags.push(RdkitEmbeddedFrag::from_coord_map(
            &coord_map,
            atoms,
            bonds,
            &adjacency,
            &degree,
            cip_ranks,
            &atom_ring_counts,
        ));
    }

    embed_cis_trans_systems(
        atoms,
        bonds,
        &adjacency,
        &degree,
        cip_ranks,
        &ring_bond_ids,
        &mut efrags,
    );

    let mut nratms = get_non_embedded_atoms(n, &efrags);
    let mut mri = if pre_spec {
        Some(0usize)
    } else {
        find_largest_frag(&efrags)
    };

    while mri.is_some() || !nratms.is_empty() {
        if mri.is_none() {
            let mut best_pos = None;
            let mut best_rank = i32::MAX;
            for (pos, &aid) in nratms.iter().enumerate() {
                let mut rank = atom_ranks[aid];
                rank *= i32::try_from(n).unwrap_or(i32::MAX);
                rank += i32::try_from(aid).unwrap_or(i32::MAX);
                if rank < best_rank {
                    best_rank = rank;
                    best_pos = Some(pos);
                }
            }
            let Some(best_pos) = best_pos else { break };
            let aid = nratms.remove(best_pos);
            efrags.push(RdkitEmbeddedFrag::from_single_atom(
                aid, atoms, bonds, &adjacency, &degree, cip_ranks,
            ));
            mri = Some(efrags.len() - 1);
        }
        let Some(idx) = mri else { break };
        efrags[idx].mark_done();
        let mut curr = std::mem::take(&mut efrags[idx]);
        curr.expand_efrag(
            &mut nratms,
            &mut efrags,
            atoms,
            bonds,
            &adjacency,
            &degree,
            cip_ranks,
        );
        // RDKit✔️❌: std::list iterators remain valid across erases of other
        // RDKit✔️❌: fragments during expandEfrag(). The Vec-backed port keeps a
        // RDKit✔️❌: placeholder at the current slot and restores the expanded
        // RDKit✔️❌: fragment into that placeholder even if earlier erases shift
        // RDKit✔️❌: its index left.
        if let Some(slot) = efrags.iter().position(|frag| frag.eatoms.is_empty()) {
            efrags[slot] = curr;
        } else {
            efrags.push(curr);
        }
        mri = find_largest_frag(&efrags);
    }

    Ok(efrags)
}

// BEGIN RDKIT CPP FUNCTION: copyCoordinate (RDDepictor.cpp)
// RDKit✔️✔️: copy embedded fragment coordinates into a flat conformer-like array.
fn copy_coordinate_from_efrags(num_atoms: usize, efrags: &[RdkitEmbeddedFrag]) -> Vec<[f64; 2]> {
    let mut out = vec![[0.0f64, 0.0f64]; num_atoms];
    for efrag in efrags {
        for (&aid, state) in efrag.get_embedded_atoms() {
            out[aid] = [state.loc.0, state.loc.1];
        }
    }
    out
}

fn efrag_component_atom_ids(efrag: &RdkitEmbeddedFrag) -> Vec<usize> {
    let mut comp: Vec<usize> = efrag.get_embedded_atoms().keys().copied().collect();
    comp.sort_unstable();
    comp
}

fn debug_depict_row_active(row: usize) -> bool {
    std::env::var("COSMOLKIT_DEBUG_DEPICT_ROW")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        == Some(row)
}

fn debug_print_efrag_stage(label: &str, efrags: &[RdkitEmbeddedFrag]) {
    if !debug_depict_row_active(58) {
        return;
    }
    eprintln!("COSMOL_STAGE label={} count={}", label, efrags.len());
    for (frag_idx, frag) in efrags.iter().enumerate() {
        let mut atom_ids: Vec<usize> = frag.get_embedded_atoms().keys().copied().collect();
        atom_ids.sort_unstable();
        eprintln!(
            "COSMOL_STAGE frag={} size={} atoms={:?}",
            frag_idx,
            atom_ids.len(),
            atom_ids
        );
        for atom_id in atom_ids {
            let state = &frag.get_embedded_atoms()[&atom_id];
            eprintln!(
                "COSMOL_STAGE atom={} loc=({:.17},{:.17}) bits=({:#018x},{:#018x})",
                atom_id,
                state.loc.0,
                state.loc.1,
                state.loc.0.to_bits(),
                state.loc.1.to_bits()
            );
        }
    }
}

// ── Public API (staged entry point) ─────────────────────────────────────────

fn require_default_compute_2d_coord_parameters(
    params: &Compute2DCoordParameters<'_>,
) -> Result<(), Coordinate2DError> {
    if !params.clear_confs {
        return Err(Coordinate2DError::UnsupportedFeature(
            "Compute2DCoordParameters.clear_confs=false is not yet ported in the active Rust compute2DCoords path".to_string(),
        ));
    }
    Ok(())
}

// BEGIN RDKIT CPP FUNCTION: Chirality::legacyStereoPerception keepGoing gate (Chirality.cpp)
// RDKit✔️✔️: bool hasStereoAtoms = false;
// RDKit✔️✔️: for (auto atom : mol.atoms()) {
// RDKit✔️✔️:   if (!hasStereoAtoms && atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
// RDKit✔️✔️:       atom->getChiralTag() != Atom::CHI_OTHER) {
// RDKit✔️✔️:     hasStereoAtoms = true;
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// RDKit✔️✔️: bool hasStereoBonds = false;
// RDKit✔️✔️: for (auto bond : mol.bonds()) {
// RDKit✔️✔️:   if (!hasStereoBonds && bond->getBondType() == Bond::DOUBLE) {
// RDKit✔️✔️:     for (auto nbond : mol.atomBonds(bond->getBeginAtom())) {
// RDKit✔️✔️:       if (nbond->getBondDir() == Bond::ENDDOWNRIGHT ||
// RDKit✔️✔️:           nbond->getBondDir() == Bond::ENDUPRIGHT) {
// RDKit✔️✔️:         hasStereoBonds = true;
// RDKit✔️✔️:         break;
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:     if (!hasStereoBonds) {
// RDKit✔️✔️:       for (auto nbond : mol.atomBonds(bond->getEndAtom())) {
// RDKit✔️✔️:         if (nbond->getBondDir() == Bond::ENDDOWNRIGHT ||
// RDKit✔️✔️:             nbond->getBondDir() == Bond::ENDUPRIGHT) {
// RDKit✔️✔️:           hasStereoBonds = true;
// RDKit✔️✔️:           break;
// RDKit✔️✔️:         }
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// RDKit✔️✔️: bool keepGoing = hasStereoAtoms | hasStereoBonds;
// RDKit✔️✔️: if (!keepGoing) {
// RDKit✔️✔️:   keepGoing = flagPossibleStereoCenters &&
// RDKit✔️✔️:               (hasPotentialStereoAtoms || hasPotentialStereoBonds);
// RDKit✔️✔️: }
fn rdkit_depict_has_stereo_bond_dirs(depict_mol: &Molecule) -> bool {
    depict_mol.bonds().iter().any(|bond| {
        if bond.order() != BondOrder::Double {
            return false;
        }
        let begin = bond.begin().index();
        let end = bond.end().index();
        depict_mol.bonds().iter().any(|nbond| {
            (nbond.begin().index() == begin
                || nbond.end().index() == begin
                || nbond.begin().index() == end
                || nbond.end().index() == end)
                && matches!(
                    nbond.direction(),
                    BondDirection::EndDownRight | BondDirection::EndUpRight
                )
        })
    })
}

fn rdkit_depict_has_potential_stereo_bonds(
    depict_mol: &Molecule,
) -> Result<bool, Coordinate2DError> {
    for bond in depict_mol.bonds() {
        if bond.order() != BondOrder::Double {
            continue;
        }
        if crate::stereo::should_detect_double_bond_stereo(depict_mol, bond.id()).map_err(
            |error| {
                Coordinate2DError::UnsupportedFeature(format!(
                    "RDKit shouldDetectDoubleBondStereo equivalent failed during 2D coordinate stereochemistry assignment: {error}"
                ))
            },
        )? {
            return Ok(true);
        }
    }
    Ok(false)
}

fn rdkit_depict_needs_atom_chiral_code_ranks(depict_mol: &Molecule) -> bool {
    depict_mol.atoms().iter().any(|atom| {
        !matches!(atom.chiral_tag(), ChiralTag::Unspecified | ChiralTag::Other)
            && atom.prop("_CIPCode").is_none()
    })
}

fn rdkit_depict_needs_bond_stereo_ranks(depict_mol: &Molecule) -> bool {
    depict_mol
        .bonds()
        .iter()
        .any(|bond| bond.order() == BondOrder::Double && bond.stereo() == BondStereo::None)
}

fn rdkit_depict_ordering_cip_ranks(
    depict_mol: &mut Molecule,
) -> Result<Vec<u32>, Coordinate2DError> {
    // BEGIN RDKIT CPP FUNCTION: MolOps::assignStereochemistry(mol, false)
    // RDKit✔️✔️: void assignStereochemistry(ROMol &mol, bool cleanIt, bool force,
    // RDKit✔️✔️:                            bool flagPossibleStereoCenters) {
    // RDKit✔️✔️:   if (!force && mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (mol.needsUpdatePropertyCache()) {
    // RDKit✔️✔️:     mol.updatePropertyCache(false);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   Chirality::legacyStereoPerception(mol, cleanIt, flagPossibleStereoCenters);
    // RDKit✔️✔️:   mol.setProp(common_properties::_StereochemDone, 1, true);
    // RDKit✔️✔️: }
    //
    // BEGIN RDKIT CPP FUNCTION: legacyStereoPerception(mol, false, false)
    // RDKit✔️✔️: bool hasStereoAtoms = false;
    // RDKit✔️✔️: bool hasPotentialStereoAtoms = false;
    // RDKit✔️✔️: for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:   if (!hasStereoAtoms && atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit✔️✔️:       atom->getChiralTag() != Atom::CHI_OTHER) {
    // RDKit✔️✔️:     hasStereoAtoms = true;
    // RDKit✔️✔️:   } else if (!hasPotentialStereoAtoms) {
    // RDKit✔️✔️:     UINT_VECT ranks;
    // RDKit✔️✔️:     Chirality::INT_PAIR_VECT nbrs;
    // RDKit✔️✔️:     hasPotentialStereoAtoms =
    // RDKit✔️✔️:         isAtomPotentialChiralCenter(atom, mol, ranks, nbrs).first;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: bool hasStereoBonds = false;
    // RDKit✔️✔️: bool hasPotentialStereoBonds = false;
    // RDKit✔️✔️: for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:   if (!hasStereoBonds && bond->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:     bool isSpecified = false;
    // RDKit✔️✔️:     ...
    // RDKit✔️✔️:     if (!hasPotentialStereoBonds && !isSpecified &&
    // RDKit✔️✔️:         shouldDetectDoubleBondStereo(bond)) {
    // RDKit✔️✔️:       hasPotentialStereoBonds = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!cleanIt && hasStereoBonds && hasPotentialStereoBonds) {
    // RDKit✔️✔️:     break;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: UINT_VECT atomRanks;
    // RDKit✔️✔️: bool keepGoing = hasStereoAtoms | hasStereoBonds;
    // RDKit✔️✔️: if (!keepGoing) {
    // RDKit✔️✔️:   keepGoing = flagPossibleStereoCenters &&
    // RDKit✔️✔️:               (hasPotentialStereoAtoms || hasPotentialStereoBonds);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: while (keepGoing) {
    // RDKit✔️✔️:   if (hasStereoAtoms || hasPotentialStereoAtoms) {
    // RDKit✔️✔️:     std::tie(hasStereoAtoms, changedStereoAtoms) =
    // RDKit✔️✔️:         Chirality::assignAtomChiralCodes(mol, atomRanks,
    // RDKit✔️✔️:                                        flagPossibleStereoCenters);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (hasStereoBonds || hasPotentialStereoBonds) {
    // RDKit✔️✔️:     std::tie(hasStereoBonds, changedStereoBonds) =
    // RDKit✔️✔️:         Chirality::assignBondStereoCodes(mol, atomRanks);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   keepGoing = (hasStereoAtoms || hasStereoBonds) &&
    // RDKit✔️✔️:               (changedStereoAtoms || changedStereoBonds);
    // RDKit✔️✔️:   if (keepGoing) {
    // RDKit✔️✔️:     Chirality::rerankAtoms(mol, atomRanks);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: legacyStereoPerception(mol, false, false)
    let existing_cip_ranks = depict_mol
        .atoms()
        .iter()
        .map(|atom| {
            atom.prop("_CIPRank")
                .and_then(|value| value.parse::<u32>().ok())
        })
        .collect::<Option<Vec<_>>>();
    // RDKit✔️✔️: if (!force && mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit✔️✔️:   return;
    // RDKit✔️✔️: }
    if depict_mol.prop("_StereochemDone").is_some() {
        return Ok(existing_cip_ranks.unwrap_or_default());
    }
    if let Some(ranks) = existing_cip_ranks {
        return Ok(ranks);
    }
    let needs_atom_ranks = rdkit_depict_needs_atom_chiral_code_ranks(depict_mol);
    let has_stereo_none_double = depict_mol
        .bonds()
        .iter()
        .any(|bond| bond.order() == BondOrder::Double && bond.stereo() == BondStereo::None);
    let runs_bond_stereo_assignment = rdkit_depict_has_stereo_bond_dirs(depict_mol)
        || rdkit_depict_has_potential_stereo_bonds(depict_mol)?;
    let needs_bond_ranks = runs_bond_stereo_assignment && has_stereo_none_double;
    if needs_atom_ranks || needs_bond_ranks {
        let ranks = crate::stereo::assign_atom_cip_ranks_in_place(depict_mol).map_err(|error| {
            Coordinate2DError::UnsupportedFeature(format!(
                "RDKit assignAtomCIPRanks equivalent failed during 2D coordinate stereochemistry assignment: {error}"
            ))
        })?;
        return Ok(ranks);
    }
    Ok(Vec::new())
}

/// Compute RDKit-compatible 2D coordinates for a molecule using the currently
/// ported default RDKit control path.
fn compute_2d_coords_default_path(
    atoms: &[Atom],
    bonds: &[Bond],
) -> Result<Vec<[f64; 2]>, Coordinate2DError> {
    let n = atoms.len();
    if n == 0 {
        return Err(Coordinate2DError::InvalidInput(
            "empty molecule has no 2D coordinates",
        ));
    }
    let mut depict_mol = build_rdkit_depict_molecule_from_slices(atoms, bonds)?;
    let cip_ranks = rdkit_depict_ordering_cip_ranks(&mut depict_mol)?;
    let depict_atoms = depict_mol.atoms();
    let depict_bonds = depict_mol.bonds();

    let mut efrags =
        rdkit_compute_initial_efrags_strict(depict_atoms, depict_bonds, &cip_ranks, None, false)?;
    shift_coords(&mut efrags);
    Ok(copy_coordinate_from_efrags(n, &efrags))
}

/// Compute RDKit-compatible 2D coordinates for a molecule.
///
/// This is the legacy narrow entry point used by current call sites.
pub(crate) fn compute_2d_coords(
    atoms: &[Atom],
    bonds: &[Bond],
) -> Result<Vec<[f64; 2]>, Coordinate2DError> {
    compute_2d_coords_with_params(atoms, bonds, &Compute2DCoordParameters::default())
}

/// Return the RDKit-depictor initial embedded-fragment scaffold before final
/// collision cleanup and sampling. ConfSeq FastGeometry uses this as a cheap,
/// source-backed ring-system prior instead of running distance geometry.
pub(crate) fn rdkit_initial_2d_scaffold_coords(
    atoms: &[Atom],
    bonds: &[Bond],
) -> Result<Vec<[f64; 2]>, Coordinate2DError> {
    let n = atoms.len();
    if n == 0 {
        return Err(Coordinate2DError::InvalidInput(
            "empty molecule has no 2D coordinates",
        ));
    }
    let mut depict_mol = build_rdkit_depict_molecule_from_slices(atoms, bonds)?;
    let cip_ranks = rdkit_depict_ordering_cip_ranks(&mut depict_mol)?;
    let depict_atoms = depict_mol.atoms();
    let depict_bonds = depict_mol.bonds();
    let mut efrags =
        rdkit_compute_initial_efrags_strict(depict_atoms, depict_bonds, &cip_ranks, None, false)?;
    shift_coords(&mut efrags);
    Ok(copy_coordinate_from_efrags(n, &efrags))
}

// BEGIN RDKIT CPP FUNCTION: compute2DCoords overload wrapper (RDDepictor.cpp)
// RDKit✔️✔️: unsigned int compute2DCoords(RDKit::ROMol &mol,
// RDKit✔️✔️:                              const RDGeom::INT_POINT2D_MAP *coordMap,
// RDKit✔️✔️:                              bool canonOrient, bool clearConfs,
// RDKit✔️✔️:                              unsigned int nFlipsPerSample,
// RDKit✔️✔️:                              unsigned int nSamples, int sampleSeed,
// RDKit✔️✔️:                              bool permuteDeg4Nodes, bool forceRDKit,
// RDKit✔️✔️:                              bool useRingTemplates) {
// RDKit✔️✔️:   Compute2DCoordParameters params;
// RDKit✔️✔️:   params.coordMap = coordMap;
// RDKit✔️✔️:   params.canonOrient = canonOrient;
// RDKit✔️✔️:   params.clearConfs = clearConfs;
// RDKit✔️✔️:   params.nFlipsPerSample = nFlipsPerSample;
// RDKit✔️✔️:   params.nSamples = nSamples;
// RDKit✔️✔️:   params.sampleSeed = sampleSeed;
// RDKit✔️✔️:   params.permuteDeg4Nodes = permuteDeg4Nodes;
// RDKit✔️✔️:   params.forceRDKit = forceRDKit;
// RDKit✔️✔️:   params.useRingTemplates = useRingTemplates;
// RDKit✔️✔️:   return compute2DCoords(mol, params);
// RDKit✔️✔️: }
pub(crate) fn compute_2d_coords_with_options(
    atoms: &[Atom],
    bonds: &[Bond],
    coord_map: Option<&BTreeMap<usize, [f64; 2]>>,
    canon_orient: bool,
    clear_confs: bool,
    n_flips_per_sample: u32,
    n_samples: u32,
    sample_seed: i32,
    permute_deg4_nodes: bool,
    force_rdkit: bool,
    use_ring_templates: bool,
) -> Result<Vec<[f64; 2]>, Coordinate2DError> {
    let params = Compute2DCoordParameters {
        coord_map,
        canon_orient,
        clear_confs,
        n_flips_per_sample,
        n_samples,
        sample_seed,
        permute_deg4_nodes,
        force_rdkit,
        use_ring_templates,
    };
    compute_2d_coords_with_params(atoms, bonds, &params)
}

// BEGIN RDKIT CPP FUNCTION: compute2DCoords params overload (RDDepictor.h / RDDepictor.cpp)
// RDKit✔️✔️: struct RDKIT_DEPICTOR_EXPORT Compute2DCoordParameters {
// RDKit✔️✔️:   const RDGeom::INT_POINT2D_MAP *coordMap = nullptr;
// RDKit✔️✔️:   bool canonOrient = false;
// RDKit✔️✔️:   bool clearConfs = true;
// RDKit✔️✔️:   unsigned int nFlipsPerSample = 0;
// RDKit✔️✔️:   unsigned int nSamples = 0;
// RDKit✔️✔️:   int sampleSeed = 0;
// RDKit✔️✔️:   bool permuteDeg4Nodes = false;
// RDKit✔️✔️:   bool forceRDKit = false;
// RDKit✔️✔️:   bool useRingTemplates = false;
// RDKit✔️✔️: };
// RDKit❗✔️: unsigned int compute2DCoords(RDKit::ROMol &mol,
// RDKit❗✔️:                              const Compute2DCoordParameters &params);
fn compute_2d_coords_with_molecule_state(
    atoms: &[Atom],
    bonds: &[Bond],
    stereochemistry_done: bool,
    params: &Compute2DCoordParameters<'_>,
) -> Result<Vec<[f64; 2]>, Coordinate2DError> {
    // BEGIN RDKIT CPP FUNCTION: compute2DCoords params overload (RDDepictor.cpp CoordGen routing)
    // RDKit✔️✔️: #ifdef RDK_BUILD_COORDGEN_SUPPORT
    // RDKit✔️✔️:   if (!params.forceRDKit && preferCoordGen) {
    // RDKit✔️✔️:     RDKit::CoordGen::CoordGenParams coordgen_params;
    // RDKit✔️✔️:     if (params.coordMap) {
    // RDKit✔️✔️:       coordgen_params.coordMap = *params.coordMap;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     auto cid = RDKit::CoordGen::addCoords(mol, &coordgen_params);
    // RDKit✔️✔️:     return cid;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️: #endif
    if !params.force_rdkit && prefer_coord_gen() {
        return Err(Coordinate2DError::UnsupportedFeature(
            "preferCoordGen=true requires CoordGen support, which is not enabled in this Rust runtime".to_string(),
        ));
    }
    require_default_compute_2d_coord_parameters(params)?;
    let n = atoms.len();
    if n == 0 {
        return Err(Coordinate2DError::InvalidInput(
            "empty molecule has no 2D coordinates",
        ));
    }
    let mut depict_mol = build_rdkit_depict_molecule_from_slices(atoms, bonds)?;
    if stereochemistry_done {
        depict_mol.properties_mut().set_prop("_StereochemDone", "1");
    }
    let cip_ranks = rdkit_depict_ordering_cip_ranks(&mut depict_mol)?;
    let depict_atoms = depict_mol.atoms();
    let depict_bonds = depict_mol.bonds();
    let mut adjacency = vec![Vec::<usize>::new(); n];
    for bond in depict_bonds {
        adjacency[bond.begin().index()].push(bond.end().index());
        adjacency[bond.end().index()].push(bond.begin().index());
    }
    let mut efrags = rdkit_compute_initial_efrags_strict(
        depict_atoms,
        depict_bonds,
        &cip_ranks,
        params.coord_map,
        params.use_ring_templates,
    )?;
    debug_print_efrag_stage("initial", &efrags);

    for frag in &mut efrags {
        let comp = efrag_component_atom_ids(frag);
        if params.n_samples > 0 && params.n_flips_per_sample > 0 {
            frag.random_sample_flips_and_permutations(
                depict_atoms,
                depict_bonds,
                &comp,
                &adjacency,
                params.n_flips_per_sample as usize,
                params.n_samples as usize,
                params.sample_seed,
                None,
                0.0,
                params.permute_deg4_nodes,
            );
        } else {
            frag.remove_collisions_bond_flip(depict_atoms, depict_bonds, &comp, &adjacency);
        }
    }
    for frag in &mut efrags {
        let comp = efrag_component_atom_ids(frag);
        frag.remove_collisions_open_angles(depict_atoms, depict_bonds, &comp, &adjacency);
    }
    for frag in &mut efrags {
        let comp = efrag_component_atom_ids(frag);
        frag.remove_collisions_shorten_bonds(depict_atoms, depict_bonds, &comp, &adjacency);
    }
    debug_print_efrag_stage("post_cleanup", &efrags);
    if params
        .coord_map
        .is_none_or(|coord_map| coord_map.is_empty())
        && params.canon_orient
    {
        for frag in &mut efrags {
            frag.canonicalize_orientation();
        }
        debug_print_efrag_stage("canonicalize", &efrags);
    }
    shift_coords(&mut efrags);
    debug_print_efrag_stage("shift", &efrags);
    let mut coords = copy_coordinate_from_efrags(n, &efrags);
    if let Some(coord_map) = params.coord_map.filter(|coord_map| coord_map.len() == 1) {
        let (&ref_idx, &ref_pos) = coord_map.iter().next().expect("single-entry coordMap");
        let shift = [
            ref_pos[0] - coords[ref_idx][0],
            ref_pos[1] - coords[ref_idx][1],
        ];
        for coord in &mut coords {
            coord[0] += shift[0];
            coord[1] += shift[1];
        }
    }
    Ok(coords)
}

pub(crate) fn compute_2d_coords_with_params(
    atoms: &[Atom],
    bonds: &[Bond],
    params: &Compute2DCoordParameters<'_>,
) -> Result<Vec<[f64; 2]>, Coordinate2DError> {
    compute_2d_coords_with_molecule_state(atoms, bonds, false, params)
}

pub(crate) fn compute_2d_coords_with_properties_and_params(
    atoms: &[Atom],
    bonds: &[Bond],
    properties: &crate::MoleculeProperties,
    params: &Compute2DCoordParameters<'_>,
) -> Result<Vec<[f64; 2]>, Coordinate2DError> {
    compute_2d_coords_with_molecule_state(
        atoms,
        bonds,
        properties.prop("_StereochemDone").is_some(),
        params,
    )
}

// BEGIN RDKIT CPP FUNCTION: compute2DCoordsMimicDistMat (RDDepictor.h / RDDepictor.cpp)
// RDKit✔️✔️: RDKIT_DEPICTOR_EXPORT unsigned int compute2DCoordsMimicDistMat(
// RDKit✔️✔️:     RDKit::ROMol &mol, const DOUBLE_SMART_PTR *dmat = nullptr,
// RDKit✔️✔️:     bool canonOrient = true, bool clearConfs = true, double weightDistMat = 0.5,
// RDKit✔️✔️:     unsigned int nFlipsPerSample = 3, unsigned int nSamples = 100,
// RDKit✔️✔️:     int sampleSeed = 25, bool permuteDeg4Nodes = true, bool forceRDKit = false);
// RDKit✔️❗: unsigned int compute2DCoordsMimicDistMat(
// RDKit✔️❗:     RDKit::ROMol &mol, const DOUBLE_SMART_PTR *dmat, bool canonOrient,
// RDKit✔️❗:     bool clearConfs, double weightDistMat, unsigned int nFlipsPerSample,
// RDKit✔️❗:     unsigned int nSamples, int sampleSeed, bool permuteDeg4Nodes, bool) {
// RDKit✔️❗:   std::list<EmbeddedFrag> efrags;
// RDKit✔️❗:   computeInitialCoords(mol, nullptr, efrags, false);
// RDKit✔️❗:   for (auto &eri : efrags) {
// RDKit✔️❗:     eri.randomSampleFlipsAndPermutations(nFlipsPerSample, nSamples, sampleSeed,
// RDKit✔️❗:                                          dmat, weightDistMat, permuteDeg4Nodes);
// RDKit✔️❗:   }
// RDKit✔️❗:   if (canonOrient && efrags.size()) {
// RDKit✔️❗:     for (auto &eri : efrags) {
// RDKit✔️❗:       eri.canonicalizeOrientation();
// RDKit✔️❗:     }
// RDKit✔️❗:   }
// RDKit✔️❗:   DepictorLocal::_shiftCoords(efrags);
// RDKit✔️❗:   return copyCoordinate(mol, efrags, clearConfs);
// RDKit✔️❗: }
pub(crate) fn compute_2d_coords_mimic_distmat_with_params(
    atoms: &[Atom],
    bonds: &[Bond],
    dmat: Option<&[f64]>,
    params: &Compute2DCoordsMimicDistMatParameters,
) -> Result<Vec<[f64; 2]>, Coordinate2DError> {
    if !params.force_rdkit && prefer_coord_gen() {
        return Err(Coordinate2DError::UnsupportedFeature(
            "preferCoordGen=true requires CoordGen support, which is not enabled in this Rust runtime".to_string(),
        ));
    }
    if !params.clear_confs {
        return Err(Coordinate2DError::UnsupportedFeature(
            "compute2DCoordsMimicDistMat clear_confs=false is not yet ported in the active Rust path".to_string(),
        ));
    }
    let n = atoms.len();
    if n == 0 {
        return Err(Coordinate2DError::InvalidInput(
            "empty molecule has no 2D coordinates",
        ));
    }

    let cip_ranks = vec![0u32; n];
    let mut adjacency = vec![Vec::<usize>::new(); n];
    for bond in bonds {
        adjacency[bond.begin().index()].push(bond.end().index());
        adjacency[bond.end().index()].push(bond.begin().index());
    }
    let mut efrags = rdkit_compute_initial_efrags_strict(atoms, bonds, &cip_ranks, None, false)?;
    for frag in &mut efrags {
        let comp = efrag_component_atom_ids(frag);
        frag.random_sample_flips_and_permutations(
            atoms,
            bonds,
            &comp,
            &adjacency,
            params.n_flips_per_sample as usize,
            params.n_samples as usize,
            params.sample_seed,
            dmat,
            params.weight_dist_mat,
            params.permute_deg4_nodes,
        );
    }
    if params.canon_orient && !efrags.is_empty() {
        for frag in &mut efrags {
            frag.canonicalize_orientation();
        }
    }
    shift_coords(&mut efrags);
    Ok(copy_coordinate_from_efrags(n, &efrags))
}

// BEGIN RDKIT CPP FUNCTION: Add2DCoordsToMol (Basement/Depictor.cpp)
// RDKit✔️✔️: int Add2DCoordsToMol(ROMol &mol, bool useDLL) {
// RDKit✔️✔️:   int res;
// RDKit✔️✔️: #ifdef WIN32_DLLBUILD
// RDKit✔️✔️:   if (useDLL) {
// RDKit✔️✔️:     res = Add2DCoordsToMolDLL(mol);
// RDKit✔️✔️:   } else {
// RDKit✔️✔️:     res = 0;
// RDKit✔️✔️:   }
// RDKit✔️✔️: #else
// RDKit✔️✔️:   res = 0;
// RDKit✔️✔️: #endif
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
pub(crate) fn add_2d_coords_to_molecule(_molecule: &mut Molecule, _use_dll: bool) -> bool {
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::{AtomSpec, Element};
    use crate::bond::BondSpec;
    use crate::builder::MoleculeBuilder;
    use std::sync::{Mutex, OnceLock};
    use tempfile::NamedTempFile;

    #[test]
    fn test_single_atom() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C));
        let mol = builder.build().unwrap();
        let coords = compute_2d_coords(mol.atoms(), mol.bonds()).unwrap();
        assert_eq!(coords.len(), 1);
        assert_eq!(coords[0], [0.0, 0.0]);
    }

    #[test]
    fn test_two_atoms() {
        let mut builder = MoleculeBuilder::new();
        let c1 = builder.add_atom(AtomSpec::new(Element::C));
        let c2 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(c1, c2, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();
        let coords = compute_2d_coords(mol.atoms(), mol.bonds()).unwrap();
        assert_eq!(coords.len(), 2);
        assert_eq!(coords[0], [0.0, 0.0]);
        let dx = coords[1][0] - coords[0][0];
        let dy = coords[1][1] - coords[0][1];
        assert!(((dx * dx + dy * dy).sqrt() - BOND_LEN).abs() < 1e-6);
    }

    #[test]
    fn test_three_linear_chain() {
        let mut builder = MoleculeBuilder::new();
        let c1 = builder.add_atom(AtomSpec::new(Element::C));
        let c2 = builder.add_atom(AtomSpec::new(Element::C));
        let c3 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(c1, c2, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c2, c3, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();
        let coords = compute_2d_coords(mol.atoms(), mol.bonds()).unwrap();
        assert_eq!(coords.len(), 3);
        assert_eq!(coords[0], [0.0, 0.0]);
        let dx01 = coords[1][0] - coords[0][0];
        let dy01 = coords[1][1] - coords[0][1];
        let dx12 = coords[2][0] - coords[1][0];
        let dy12 = coords[2][1] - coords[1][1];
        assert!(((dx01 * dx01 + dy01 * dy01).sqrt() - BOND_LEN).abs() < 1e-6);
        assert!(((dx12 * dx12 + dy12 * dy12).sqrt() - BOND_LEN).abs() < 1e-6);
        assert_ne!(coords[1], coords[2]);
        assert_ne!(coords[0], coords[2]);
    }

    #[test]
    fn test_atom_depict_rank() {
        assert_eq!(atom_depict_rank(6, 2), 602); // C, deg 2
        assert_eq!(atom_depict_rank(1, 1), 100001); // H, deg 1 (special: anum=1000)
        assert_eq!(atom_depict_rank(8, 1), 801); // O, deg 1
    }

    #[test]
    fn test_rdkit_ring_radius() {
        let r = rdkit_ring_radius(6, 1.5);
        assert!((r - 1.5).abs() < 0.01); // hexagon radius ~= bond length
    }

    #[test]
    fn test_compute_sub_angle() {
        let a = compute_sub_angle(4, RdkitHybridization::Sp3);
        assert!((a - PI / 2.0).abs() < 1e-6);
        let a = compute_sub_angle(3, RdkitHybridization::Sp2);
        assert!((a - 2.0 * PI / 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_canonicalize_component_single() {
        let mut comp = vec![(0usize, (1.0, 2.0))];
        canonicalize_component(&mut comp);
        assert_eq!(comp[0].1, (1.0, 2.0));
    }

    #[test]
    fn test_canonicalize_component_two() {
        let mut comp = vec![(0usize, (1.0, 0.0)), (1usize, (-1.0, 0.0))];
        canonicalize_component(&mut comp);
        // Should preserve relative positions, just rotate to principal axes
        assert_eq!(comp.len(), 2);
    }

    #[test]
    fn compute2d_params_defaults_match_rdkit_header_defaults() {
        let params = Compute2DCoordParameters::default();
        assert!(params.coord_map.is_none());
        assert!(!params.canon_orient);
        assert!(params.clear_confs);
        assert_eq!(params.n_flips_per_sample, 0);
        assert_eq!(params.n_samples, 0);
        assert_eq!(params.sample_seed, 0);
        assert!(!params.permute_deg4_nodes);
        assert!(!params.force_rdkit);
        assert!(!params.use_ring_templates);
    }

    #[test]
    fn compute2d_mimic_distmat_params_defaults_match_rdkit_header_defaults() {
        let params = Compute2DCoordsMimicDistMatParameters::default();
        assert!(params.canon_orient);
        assert!(params.clear_confs);
        assert!((params.weight_dist_mat - 0.5).abs() < 1.0e-12);
        assert_eq!(params.n_flips_per_sample, 3);
        assert_eq!(params.n_samples, 100);
        assert_eq!(params.sample_seed, 25);
        assert!(params.permute_deg4_nodes);
        assert!(!params.force_rdkit);
    }

    #[test]
    fn compute2d_default_params_and_options_entrypoints_produce_same_coords() {
        let mut builder = MoleculeBuilder::new();
        let c1 = builder.add_atom(AtomSpec::new(Element::C));
        let c2 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(c1, c2, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();

        let from_params = compute_2d_coords_with_params(
            mol.atoms(),
            mol.bonds(),
            &Compute2DCoordParameters::default(),
        )
        .unwrap();
        let from_options = compute_2d_coords_with_options(
            mol.atoms(),
            mol.bonds(),
            None,
            false,
            true,
            0,
            0,
            0,
            false,
            false,
            false,
        )
        .unwrap();

        assert_eq!(from_params, from_options);
    }

    #[test]
    fn compute2d_options_rejects_clear_confs_false_until_molecule_writeback_is_ported() {
        let mut builder = MoleculeBuilder::new();
        let c1 = builder.add_atom(AtomSpec::new(Element::C));
        let c2 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(c1, c2, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();

        let err = compute_2d_coords_with_options(
            mol.atoms(),
            mol.bonds(),
            None,
            false,
            false,
            0,
            0,
            0,
            false,
            false,
            false,
        )
        .unwrap_err();

        match err {
            Coordinate2DError::UnsupportedFeature(message) => {
                assert!(message.contains("clear_confs=false"));
            }
            other => panic!("expected UnsupportedFeature, got {other:?}"),
        }
    }

    #[test]
    fn compute2d_mimic_distmat_rejects_clear_confs_false_until_molecule_writeback_is_ported() {
        let mut builder = MoleculeBuilder::new();
        let c1 = builder.add_atom(AtomSpec::new(Element::C));
        let c2 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(c1, c2, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();

        let err = compute_2d_coords_mimic_distmat_with_params(
            mol.atoms(),
            mol.bonds(),
            None,
            &Compute2DCoordsMimicDistMatParameters {
                clear_confs: false,
                ..Compute2DCoordsMimicDistMatParameters::default()
            },
        )
        .unwrap_err();

        match err {
            Coordinate2DError::UnsupportedFeature(message) => {
                assert!(message.contains("clear_confs=false"));
            }
            other => panic!("expected UnsupportedFeature, got {other:?}"),
        }
    }

    #[test]
    fn compute2d_single_atom_coord_map_translates_final_layout_to_reference_atom() {
        let _lock = prefer_coordgen_test_lock().lock().unwrap();
        let _guard = PreferCoordGenGuard::capture();
        set_prefer_coord_gen(false);
        let mut builder = MoleculeBuilder::new();
        let c1 = builder.add_atom(AtomSpec::new(Element::C));
        let c2 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(c1, c2, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();

        let baseline = compute_2d_coords_with_options(
            mol.atoms(),
            mol.bonds(),
            None,
            false,
            true,
            0,
            0,
            0,
            false,
            false,
            false,
        )
        .unwrap();

        let mut coord_map = BTreeMap::new();
        coord_map.insert(0usize, [5.0, -2.0]);
        let coords = compute_2d_coords_with_options(
            mol.atoms(),
            mol.bonds(),
            Some(&coord_map),
            false,
            true,
            0,
            0,
            0,
            false,
            false,
            false,
        )
        .unwrap();

        assert_eq!(coords[0], [5.0, -2.0]);
        let baseline_delta = [
            baseline[1][0] - baseline[0][0],
            baseline[1][1] - baseline[0][1],
        ];
        let shifted_delta = [coords[1][0] - coords[0][0], coords[1][1] - coords[0][1]];
        assert!((shifted_delta[0] - baseline_delta[0]).abs() < 1.0e-12);
        assert!((shifted_delta[1] - baseline_delta[1]).abs() < 1.0e-12);
    }

    #[test]
    fn compute2d_canon_orient_runs_without_prespecified_coords() {
        let _lock = prefer_coordgen_test_lock().lock().unwrap();
        let _guard = PreferCoordGenGuard::capture();
        set_prefer_coord_gen(false);

        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        let a2 = builder.add_atom(AtomSpec::new(Element::C));
        let a3 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(a1, a2, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(a1, a3, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();

        let coords_without = compute_2d_coords_with_options(
            mol.atoms(),
            mol.bonds(),
            None,
            false,
            true,
            0,
            0,
            0,
            false,
            false,
            false,
        )
        .unwrap();

        let coords_with = compute_2d_coords_with_options(
            mol.atoms(),
            mol.bonds(),
            None,
            true,
            true,
            0,
            0,
            0,
            false,
            false,
            false,
        )
        .unwrap();

        assert_eq!(coords_with.len(), 4);
        assert!(
            coords_with
                .iter()
                .all(|coord| coord[0].is_finite() && coord[1].is_finite())
        );
        assert_ne!(coords_with, coords_without);
    }

    #[test]
    fn compute2d_prespec_fragment_order_does_not_depend_on_component_index_zip() {
        let _lock = prefer_coordgen_test_lock().lock().unwrap();
        let _guard = PreferCoordGenGuard::capture();
        set_prefer_coord_gen(false);
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        let a2 = builder.add_atom(AtomSpec::new(Element::C));
        let a3 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(a2, a3, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();

        let mut coord_map = BTreeMap::new();
        coord_map.insert(2usize, [10.0, 10.0]);
        coord_map.insert(3usize, [11.5, 10.0]);

        let coords = compute_2d_coords_with_options(
            mol.atoms(),
            mol.bonds(),
            Some(&coord_map),
            false,
            true,
            0,
            0,
            0,
            false,
            false,
            false,
        )
        .unwrap();

        assert_eq!(coords[2], [10.0, 10.0]);
        assert_eq!(coords[3], [11.5, 10.0]);
        assert!(coords[0][0].is_finite() && coords[0][1].is_finite());
        assert!(coords[1][0].is_finite() && coords[1][1].is_finite());
    }

    #[test]
    fn compute2d_options_and_params_match_for_non_default_supported_flags() {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..4)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect();
        for (a, b) in [(0, 1), (1, 2), (2, 3)] {
            builder
                .add_bond(BondSpec::new(atoms[a], atoms[b], BondOrder::Single))
                .unwrap();
        }
        let mol = builder.build().unwrap();

        let params = Compute2DCoordParameters {
            coord_map: Some(&BTreeMap::from([(0usize, [1.5, -0.5])])),
            canon_orient: true,
            clear_confs: true,
            n_flips_per_sample: 1,
            n_samples: 1,
            sample_seed: 7,
            permute_deg4_nodes: true,
            force_rdkit: true,
            use_ring_templates: true,
        };
        let from_params = compute_2d_coords_with_params(mol.atoms(), mol.bonds(), &params).unwrap();
        let from_options = compute_2d_coords_with_options(
            mol.atoms(),
            mol.bonds(),
            params.coord_map,
            params.canon_orient,
            params.clear_confs,
            params.n_flips_per_sample,
            params.n_samples,
            params.sample_seed,
            params.permute_deg4_nodes,
            params.force_rdkit,
            params.use_ring_templates,
        )
        .unwrap();

        assert_eq!(from_options, from_params);
    }

    #[test]
    fn compute2d_options_and_params_match_for_ring_template_flag_on_ring_system() {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect();
        for (a, b) in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)] {
            builder
                .add_bond(BondSpec::new(atoms[a], atoms[b], BondOrder::Single))
                .unwrap();
        }
        let mol = builder.build().unwrap();

        let params = Compute2DCoordParameters {
            force_rdkit: true,
            use_ring_templates: true,
            ..Compute2DCoordParameters::default()
        };
        let from_params = compute_2d_coords_with_params(mol.atoms(), mol.bonds(), &params).unwrap();
        let from_options = compute_2d_coords_with_options(
            mol.atoms(),
            mol.bonds(),
            None,
            false,
            true,
            0,
            0,
            0,
            false,
            true,
            true,
        )
        .unwrap();

        assert_eq!(from_options, from_params);
    }

    #[test]
    fn compute2d_mimic_distmat_ignores_negative_entries_and_matches_missing_reference() {
        let _lock = prefer_coordgen_test_lock().lock().unwrap();
        let _guard = PreferCoordGenGuard::capture();
        set_prefer_coord_gen(false);
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        let a2 = builder.add_atom(AtomSpec::new(Element::C));
        let a3 = builder.add_atom(AtomSpec::new(Element::C));
        for (a, b) in [(a0, a1), (a1, a2), (a2, a3)] {
            builder
                .add_bond(BondSpec::new(a, b, BondOrder::Single))
                .unwrap();
        }
        let mol = builder.build().unwrap();
        let params = Compute2DCoordsMimicDistMatParameters {
            weight_dist_mat: 1.0,
            n_flips_per_sample: 1,
            n_samples: 8,
            sample_seed: 17,
            permute_deg4_nodes: false,
            ..Compute2DCoordsMimicDistMatParameters::default()
        };

        let without_ref =
            compute_2d_coords_mimic_distmat_with_params(mol.atoms(), mol.bonds(), None, &params)
                .unwrap();
        let negative_ref = vec![-1.0; mol.num_atoms() * (mol.num_atoms() - 1) / 2];
        let with_negative = compute_2d_coords_mimic_distmat_with_params(
            mol.atoms(),
            mol.bonds(),
            Some(&negative_ref),
            &params,
        )
        .unwrap();

        assert_eq!(without_ref, with_negative);
    }

    #[test]
    fn compute2d_mimic_distmat_weight_changes_sample_selection() {
        let _lock = prefer_coordgen_test_lock().lock().unwrap();
        let _guard = PreferCoordGenGuard::capture();
        set_prefer_coord_gen(false);
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..4)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect();
        for (a, b) in [(0, 1), (1, 2), (2, 3)] {
            builder
                .add_bond(BondSpec::new(atoms[a], atoms[b], BondOrder::Single))
                .unwrap();
        }
        let mol = builder.build().unwrap();
        let dmat = vec![1.5, 3.0, 0.5, 1.5, 3.0, 1.5];

        let density_only = compute_2d_coords_mimic_distmat_with_params(
            mol.atoms(),
            mol.bonds(),
            Some(&dmat),
            &Compute2DCoordsMimicDistMatParameters {
                weight_dist_mat: 0.0,
                n_flips_per_sample: 1,
                n_samples: 12,
                sample_seed: 19,
                permute_deg4_nodes: false,
                ..Compute2DCoordsMimicDistMatParameters::default()
            },
        )
        .unwrap();
        let mimic_only = compute_2d_coords_mimic_distmat_with_params(
            mol.atoms(),
            mol.bonds(),
            Some(&dmat),
            &Compute2DCoordsMimicDistMatParameters {
                weight_dist_mat: 1.0,
                n_flips_per_sample: 1,
                n_samples: 12,
                sample_seed: 19,
                permute_deg4_nodes: false,
                ..Compute2DCoordsMimicDistMatParameters::default()
            },
        )
        .unwrap();

        assert_ne!(density_only, mimic_only);
    }

    #[test]
    fn compute2d_mimic_distmat_is_deterministic_for_same_positive_seed() {
        let _lock = prefer_coordgen_test_lock().lock().unwrap();
        let _guard = PreferCoordGenGuard::capture();
        set_prefer_coord_gen(false);
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..5)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect();
        for (a, b) in [(0, 1), (1, 2), (2, 3), (3, 4)] {
            builder
                .add_bond(BondSpec::new(atoms[a], atoms[b], BondOrder::Single))
                .unwrap();
        }
        let mol = builder.build().unwrap();
        let dmat = vec![1.5, 3.0, 4.5, 6.0, 1.5, 3.0, 4.5, 1.5, 3.0, 1.5];
        let params = Compute2DCoordsMimicDistMatParameters {
            weight_dist_mat: 0.75,
            n_flips_per_sample: 1,
            n_samples: 16,
            sample_seed: 23,
            permute_deg4_nodes: false,
            ..Compute2DCoordsMimicDistMatParameters::default()
        };

        let coords1 = compute_2d_coords_mimic_distmat_with_params(
            mol.atoms(),
            mol.bonds(),
            Some(&dmat),
            &params,
        )
        .unwrap();
        let coords2 = compute_2d_coords_mimic_distmat_with_params(
            mol.atoms(),
            mol.bonds(),
            Some(&dmat),
            &params,
        )
        .unwrap();

        assert_eq!(coords1, coords2);
    }

    #[test]
    fn prefer_coordgen_defaults_false_like_rdkit_global() {
        set_prefer_coord_gen(false);
        assert!(!prefer_coord_gen());
        assert!(!is_coordgen_support_available());
    }

    #[test]
    fn compute2d_prefers_coordgen_only_when_not_forced_to_rdkit() {
        let _lock = prefer_coordgen_test_lock().lock().unwrap();
        let _guard = PreferCoordGenGuard::capture();
        let mut builder = MoleculeBuilder::new();
        let c1 = builder.add_atom(AtomSpec::new(Element::C));
        let c2 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(c1, c2, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();

        set_prefer_coord_gen(true);
        let err = compute_2d_coords_with_params(
            mol.atoms(),
            mol.bonds(),
            &Compute2DCoordParameters::default(),
        )
        .unwrap_err();
        match err {
            Coordinate2DError::UnsupportedFeature(message) => {
                assert!(message.contains("CoordGen"));
            }
            other => panic!("expected UnsupportedFeature, got {other:?}"),
        }

        let coords = compute_2d_coords_with_params(
            mol.atoms(),
            mol.bonds(),
            &Compute2DCoordParameters {
                force_rdkit: true,
                ..Compute2DCoordParameters::default()
            },
        )
        .unwrap();
        assert_eq!(coords.len(), 2);
        set_prefer_coord_gen(false);
    }

    #[test]
    fn add2dcoords_to_mol_wrapper_returns_false_without_windows_dll_path() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C));
        let mut mol = builder.build().unwrap();
        assert!(!add_2d_coords_to_molecule(&mut mol, true));
        assert!(!add_2d_coords_to_molecule(&mut mol, false));
    }

    #[test]
    fn embed_ring_places_regular_polygon_with_rdkit_bond_length() {
        let coords = rdkit_embed_ring(&[10, 11, 12, 13, 14, 15]);
        assert_eq!(coords.len(), 6);
        for point in coords.values() {
            assert!((norm(*point) - 1.5).abs() < 1e-8);
        }
        let ordered = [10usize, 11, 12, 13, 14, 15]
            .into_iter()
            .map(|idx| coords[&idx])
            .collect::<Vec<_>>();
        for window in ordered.windows(2) {
            let delta = (window[0].0 - window[1].0, window[0].1 - window[1].1);
            assert!((norm(delta) - 1.5).abs() < 1e-8);
        }
        let wrap = (ordered[0].0 - ordered[5].0, ordered[0].1 - ordered[5].1);
        assert!((norm(wrap) - 1.5).abs() < 1e-8);
    }

    #[test]
    fn transform_points_applies_two_point_affine_mapping() {
        let mut points = BTreeMap::from([(0usize, (0.0, 0.0)), (1usize, (1.0, 0.0))]);
        let trans =
            transform2d_set_transform_two_point((2.0, 3.0), (2.0, 4.0), (0.0, 0.0), (1.0, 0.0));
        rdkit_transform_points(&mut points, trans);
        assert!((points[&0].0 - 2.0).abs() < 1e-8);
        assert!((points[&0].1 - 3.0).abs() < 1e-8);
        assert!((points[&1].0 - 2.0).abs() < 1e-8);
        assert!((points[&1].1 - 4.0).abs() < 1e-8);
    }

    #[test]
    fn rank_atoms_uses_cip_ranks_before_depict_rank() {
        let mut builder = MoleculeBuilder::new();
        let c = builder.add_atom(AtomSpec::new(Element::C).with_prop("_CIPRank", "30"));
        let o = builder.add_atom(AtomSpec::new(Element::O).with_prop("_CIPRank", "10"));
        let n = builder.add_atom(AtomSpec::new(Element::N).with_prop("_CIPRank", "20"));
        builder
            .add_bond(BondSpec::new(c, o, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c, n, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();
        let degree = vec![2usize, 1, 1];
        let order = rdkit_rank_atoms_by_rank(mol.atoms(), &[0, 1, 2], &degree, true);
        assert_eq!(order, vec![1, 2, 0]);
    }

    #[test]
    fn rank_atoms_uses_chiral_atom_rank_when_cip_rank_prop_is_absent() {
        let mut builder = MoleculeBuilder::new();
        let c = builder.add_atom(AtomSpec::new(Element::C).with_prop("_ChiralAtomRank", "0"));
        let h = builder.add_atom(AtomSpec::new(Element::H).with_prop("_ChiralAtomRank", "1"));
        builder
            .add_bond(BondSpec::new(c, h, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();
        let degree = vec![1usize, 1];
        let order = rdkit_rank_atoms_by_rank(mol.atoms(), &[0, 1], &degree, true);
        assert_eq!(order, vec![1, 0]);
    }

    #[test]
    fn rank_atoms_falls_back_to_depict_rank_and_supports_descending() {
        let mut builder = MoleculeBuilder::new();
        let c = builder.add_atom(AtomSpec::new(Element::C));
        let o = builder.add_atom(AtomSpec::new(Element::O));
        let h = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(BondSpec::new(c, o, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c, h, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();
        let degree = vec![2usize, 1, 1];
        let ascending = rdkit_rank_atoms_by_rank(mol.atoms(), &[0, 1, 2], &degree, true);
        let descending = rdkit_rank_atoms_by_rank(mol.atoms(), &[0, 1, 2], &degree, false);
        assert_eq!(ascending, vec![0, 1, 2]);
        assert_eq!(descending, vec![2, 1, 0]);
    }

    #[test]
    fn set_nbr_order_returns_ranked_neighbors_when_no_reference_neighbor_exists() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C));
        let a = builder.add_atom(AtomSpec::new(Element::O));
        let b = builder.add_atom(AtomSpec::new(Element::N));
        builder
            .add_bond(BondSpec::new(center, a, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, b, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();
        let adjacency = vec![vec![1, 2], vec![0], vec![0]];
        let degree = vec![2usize, 1, 1];
        let ordered =
            rdkit_set_nbr_order(0, &[2, 1], mol.atoms(), mol.bonds(), &adjacency, &degree);
        assert_eq!(ordered, vec![2, 1]);
    }

    #[test]
    fn set_nbr_order_rotates_around_reference_neighbor_for_degree_four_center() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C));
        let a = builder.add_atom(AtomSpec::new(Element::O));
        let b = builder.add_atom(AtomSpec::new(Element::N));
        let c = builder.add_atom(AtomSpec::new(Element::F));
        let d = builder.add_atom(AtomSpec::new(Element::CL));
        for nbr in [a, b, c, d] {
            builder
                .add_bond(BondSpec::new(center, nbr, BondOrder::Single))
                .unwrap();
        }
        let mol = builder.build().unwrap();
        let adjacency = vec![vec![1, 2, 3, 4], vec![0], vec![0], vec![0], vec![0]];
        let degree = vec![4usize, 1, 1, 1, 1];
        let ordered =
            rdkit_set_nbr_order(0, &[1, 2, 3], mol.atoms(), mol.bonds(), &adjacency, &degree);
        assert_eq!(ordered, vec![2, 3, 1]);
    }

    #[test]
    fn nbr_atom_bond_ids_returns_neighbor_and_bond_indices_in_adjacency_order() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::O));
        let a2 = builder.add_atom(AtomSpec::new(Element::N));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(a0, a2, BondOrder::Double))
            .unwrap();
        let mol = builder.build().unwrap();
        let adjacency = vec![vec![1, 2], vec![0], vec![0]];
        let (aids, bids) = rdkit_get_nbr_atom_and_bond_ids(0, mol.bonds(), &adjacency);
        assert_eq!(aids, vec![1, 2]);
        assert_eq!(bids, vec![0, 1]);
    }

    #[test]
    fn permute_deg4_returns_perpendicular_pairs_when_first_two_neighbors_are_orthogonal() {
        let pairs = rdkit_find_bond_pairs_to_permute_deg4(
            (0.0, 0.0),
            &[10, 11, 12, 13],
            &[(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)],
        );
        assert_eq!(pairs, vec![(10, 11), (10, 13)]);
    }

    #[test]
    fn permute_deg4_returns_second_and_third_pairs_when_first_two_neighbors_are_opposite() {
        let pairs = rdkit_find_bond_pairs_to_permute_deg4(
            (0.0, 0.0),
            &[20, 21, 22, 23],
            &[(1.0, 0.0), (-1.0, 0.0), (0.0, 1.0), (0.0, -1.0)],
        );
        assert_eq!(pairs, vec![(20, 22), (20, 23)]);
    }

    #[test]
    fn pick_first_ring_prefers_fewer_substituents() {
        let degree = vec![3usize, 2, 2, 2, 2, 2, 2, 2];
        let fused_rings = vec![vec![0, 1, 2, 3, 4, 5], vec![1, 2, 6, 7]];
        let picked = rdkit_pick_first_ring_to_embed(&degree, &fused_rings);
        assert_eq!(picked, 1);
    }

    #[test]
    fn pick_first_ring_breaks_equal_substituent_ties_by_larger_ring() {
        let degree = vec![2usize, 2, 2, 2, 2, 2, 2];
        let fused_rings = vec![vec![0, 1, 2, 3], vec![1, 2, 4, 5, 6]];
        let picked = rdkit_pick_first_ring_to_embed(&degree, &fused_rings);
        assert_eq!(picked, 1);
    }

    #[test]
    fn find_core_rings_iteratively_removes_singly_fused_side_rings_in_iteration_order() {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..14)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect();
        for (a, b) in [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (2, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 3),
            (7, 10),
            (10, 11),
            (11, 12),
            (12, 13),
            (13, 6),
        ] {
            builder
                .add_bond(BondSpec::new(atoms[a], atoms[b], BondOrder::Single))
                .unwrap();
        }
        let mol = builder.build().unwrap();
        let fused_rings = vec![
            vec![0, 1, 2, 3, 4, 5],
            vec![2, 3, 6, 7, 8, 9],
            vec![6, 7, 10, 11, 12, 13],
        ];
        let (core_rings, core_ring_ids) = rdkit_find_core_rings(&fused_rings, mol.bonds());
        assert_eq!(core_rings, vec![vec![6, 7, 10, 11, 12, 13]]);
        assert_eq!(core_ring_ids, vec![2]);
    }

    #[test]
    fn find_core_rings_stops_shared_atom_scan_after_first_two_hits() {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..7)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect();
        for (a, b) in [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 0),
            (4, 5),
            (5, 6),
            (6, 2),
        ] {
            builder
                .add_bond(BondSpec::new(atoms[a], atoms[b], BondOrder::Single))
                .unwrap();
        }
        let mol = builder.build().unwrap();
        let fused_rings = vec![vec![0, 1, 2, 3, 4], vec![2, 3, 4, 5, 6]];
        let (core_rings, core_ring_ids) = rdkit_find_core_rings(&fused_rings, mol.bonds());
        assert_eq!(core_rings, vec![vec![2, 3, 4, 5, 6]]);
        assert_eq!(core_ring_ids, vec![1]);
    }

    #[test]
    fn find_next_ring_prefers_two_common_atoms_before_larger_overlap() {
        let done_rings = vec![0usize];
        let fused_rings = vec![
            vec![0, 1, 2, 3, 4, 5],
            vec![2, 3, 6, 7, 8, 9],
            vec![1, 2, 3, 10, 11, 12],
        ];
        let (next_id, common_atoms) =
            rdkit_find_next_ring_to_embed(&done_rings, &fused_rings).unwrap();
        assert_eq!(next_id, 1);
        assert_eq!(common_atoms, vec![2, 3]);
    }

    #[test]
    fn find_next_ring_rotates_common_atom_chain_to_ring_order_endpoints() {
        let done_rings = vec![0usize];
        let fused_rings = vec![vec![0, 1, 2, 3, 4, 5], vec![4, 3, 6, 7, 8, 5]];
        let (next_id, common_atoms) =
            rdkit_find_next_ring_to_embed(&done_rings, &fused_rings).unwrap();
        assert_eq!(next_id, 1);
        assert_eq!(common_atoms, vec![5, 4, 3]);
    }

    #[test]
    fn get_all_rotatable_bonds_filters_out_ring_and_stereo_bonds() {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..8)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect();
        for (idx, (a, b)) in [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 6), (6, 7)]
            .into_iter()
            .enumerate()
        {
            let spec = if idx == 4 {
                BondSpec::new(atoms[a], atoms[b], BondOrder::Single)
                    .with_stereo_atoms(atoms[3], atoms[6])
                    .with_stereo(BondStereo::Cis)
            } else {
                BondSpec::new(atoms[a], atoms[b], BondOrder::Single)
            };
            builder.add_bond(spec).unwrap();
        }
        let mol = builder.build().unwrap();
        let ring_bond_ids = BTreeSet::from([0usize, 1, 2]);
        let rotatable =
            rdkit_get_all_rotatable_bonds(mol.bonds(), |bid| ring_bond_ids.contains(&bid));
        assert_eq!(rotatable, vec![3, 5, 6]);
    }

    #[test]
    fn get_rotatable_bonds_uses_shortest_path_without_endpoint_bonds() {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..5)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect();
        for (a, b) in [(0, 1), (1, 2), (2, 3), (3, 4)] {
            builder
                .add_bond(BondSpec::new(atoms[a], atoms[b], BondOrder::Single))
                .unwrap();
        }
        let mol = builder.build().unwrap();
        let adjacency = vec![vec![1], vec![0, 2], vec![1, 3], vec![2, 4], vec![3]];
        let rotatable = rdkit_get_rotatable_bonds_between(
            0,
            4,
            mol.num_atoms(),
            mol.bonds(),
            &adjacency,
            |_| false,
        );
        assert_eq!(rotatable, vec![1, 2]);
    }

    #[test]
    fn has_terminal_rgroup_or_query_hydrogen_detects_terminal_dummy_and_atomic_number_query_h() {
        let mut builder = MoleculeBuilder::new();
        let q = builder.add_atom(
            AtomSpec::new(Element::C).with_query(QueryNode::predicate(AtomQueryPredicate::Any)),
        );
        let dummy = builder.add_atom(AtomSpec::new(Element::from_atomic_number(0).unwrap()));
        let query_h = builder.add_atom(
            AtomSpec::new(Element::from_atomic_number(0).unwrap())
                .with_query(QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1))),
        );
        builder
            .add_bond(BondSpec::new(q, dummy, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(q, query_h, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();
        assert!(rdkit_has_terminal_rgroup_or_query_hydrogen(&mol));
    }

    #[test]
    fn prepare_template_for_rgroups_removes_terminal_dummy_and_records_former_indices() {
        let mut builder = MoleculeBuilder::new();
        let query_center = builder.add_atom(
            AtomSpec::new(Element::C).with_query(QueryNode::predicate(AtomQueryPredicate::Any)),
        );
        let dummy = builder.add_atom(AtomSpec::new(Element::from_atomic_number(0).unwrap()));
        builder
            .add_bond(BondSpec::new(query_center, dummy, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();

        let reduced = rdkit_prepare_template_for_rgroups(&mol).unwrap();
        assert_eq!(reduced.num_atoms(), 1);
        assert_eq!(reduced.atoms()[0].prop(RDKIT_FORMER_IDX_PROP), Some("0"));
        assert_eq!(
            reduced.atoms()[0].prop(RDKIT_FORMER_NBR_INDICES_PROP),
            Some("1")
        );
    }

    #[test]
    fn prepare_template_for_rgroups_preserves_query_hydrogen_neighbor_metadata() {
        let mut builder = MoleculeBuilder::new();
        let query_center = builder.add_atom(
            AtomSpec::new(Element::C)
                .with_aromatic(true)
                .with_query(QueryNode::predicate(AtomQueryPredicate::Any)),
        );
        let query_h = builder.add_atom(
            AtomSpec::new(Element::from_atomic_number(0).unwrap())
                .with_query(QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1))),
        );
        let carbon = builder.add_atom(AtomSpec::new(Element::C).with_aromatic(true));
        builder
            .add_bond(BondSpec::new(query_center, query_h, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(query_center, carbon, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();

        let reduced = rdkit_prepare_template_for_rgroups(&mol).unwrap();
        assert_eq!(reduced.num_atoms(), 2);
        assert_eq!(reduced.atoms()[0].prop(RDKIT_FORMER_IDX_PROP), Some("0"));
        assert_eq!(
            reduced.atoms()[0].prop(RDKIT_FORMER_NBR_INDICES_PROP),
            Some("1")
        );
        assert_eq!(reduced.atoms()[1].prop(RDKIT_FORMER_IDX_PROP), Some("2"));
        let query = reduced.bonds()[0].query();
        assert_eq!(
            query,
            Some(&QueryNode::predicate(crate::BondQueryPredicate::OrderIn(
                vec![BondOrder::Single, BondOrder::Aromatic,]
            )))
        );
    }

    #[test]
    fn reduced_to_full_matches_restores_former_indices_and_appends_unmatched_neighbors() {
        let mut reduced_query_builder = MoleculeBuilder::new();
        let rq0 = reduced_query_builder.add_atom(
            AtomSpec::new(Element::C)
                .with_query(QueryNode::predicate(AtomQueryPredicate::Any))
                .with_prop(RDKIT_FORMER_IDX_PROP, "0")
                .with_prop(RDKIT_FORMER_NBR_INDICES_PROP, "2,3"),
        );
        let rq1 = reduced_query_builder.add_atom(
            AtomSpec::new(Element::N)
                .with_query(QueryNode::predicate(AtomQueryPredicate::Any))
                .with_prop(RDKIT_FORMER_IDX_PROP, "1"),
        );
        reduced_query_builder
            .add_bond(BondSpec::new(rq0, rq1, BondOrder::Single))
            .unwrap();
        let reduced_query = reduced_query_builder.build().unwrap();

        let mut mol_builder = MoleculeBuilder::new();
        let m0 = mol_builder.add_atom(AtomSpec::new(Element::C));
        let m1 = mol_builder.add_atom(AtomSpec::new(Element::N));
        let m2 = mol_builder.add_atom(AtomSpec::new(Element::O));
        let m3 = mol_builder.add_atom(AtomSpec::new(Element::F));
        mol_builder
            .add_bond(BondSpec::new(m0, m1, BondOrder::Single))
            .unwrap();
        mol_builder
            .add_bond(BondSpec::new(m0, m2, BondOrder::Single))
            .unwrap();
        mol_builder
            .add_bond(BondSpec::new(m0, m3, BondOrder::Single))
            .unwrap();
        let mol_hs = mol_builder.build().unwrap();

        let mut matches = vec![vec![(0usize, 0usize), (1usize, 1usize)]];
        rdkit_reduced_to_full_matches(&reduced_query, &mol_hs, &mut matches);
        assert_eq!(matches, vec![vec![(0, 0), (1, 1), (3, 2), (2, 3)]]);
    }

    #[test]
    fn invert_wedging_if_mol_has_flipped_respects_threshold_and_inverts_props() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        let a2 = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(
                BondSpec::new(a0, a1, BondOrder::Single)
                    .with_prop("_MolFileBondCfg", "1")
                    .with_prop("_MolFileBondStereo", "1"),
            )
            .unwrap();
        builder
            .add_bond(
                BondSpec::new(a1, a2, BondOrder::Single)
                    .with_prop("_MolFileBondCfg", "2")
                    .with_prop("_MolFileBondStereo", "4"),
            )
            .unwrap();
        let mut mol = builder.build().unwrap();

        let not_flipped = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, -0.99, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        assert!(!rdkit_invert_wedging_if_mol_has_flipped(
            &mut mol,
            &not_flipped
        ));
        assert_eq!(mol.bonds()[0].prop("_MolFileBondCfg"), Some("1"));
        assert_eq!(mol.bonds()[0].prop("_MolFileBondStereo"), Some("1"));
        assert_eq!(mol.bonds()[1].prop("_MolFileBondCfg"), Some("2"));
        assert_eq!(mol.bonds()[1].prop("_MolFileBondStereo"), Some("4"));

        let flipped = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, -0.991, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        assert!(rdkit_invert_wedging_if_mol_has_flipped(&mut mol, &flipped));
        assert_eq!(mol.bonds()[0].prop("_MolFileBondCfg"), Some("3"));
        assert_eq!(mol.bonds()[0].prop("_MolFileBondStereo"), Some("6"));
        assert_eq!(mol.bonds()[1].prop("_MolFileBondCfg"), Some("2"));
        assert_eq!(mol.bonds()[1].prop("_MolFileBondStereo"), Some("4"));
    }

    #[test]
    fn clear_molblock_wedging_info_removes_stereo_and_cfg_props() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(
                BondSpec::new(a0, a1, BondOrder::Single)
                    .with_prop("_MolFileBondCfg", "1")
                    .with_prop("_MolFileBondStereo", "6"),
            )
            .unwrap();
        let mut mol = builder.build().unwrap();

        rdkit_clear_molblock_wedging_info(&mut mol);

        assert_eq!(mol.bonds()[0].prop("_MolFileBondCfg"), None);
        assert_eq!(mol.bonds()[0].prop("_MolFileBondStereo"), None);
    }

    #[test]
    fn remove_all_2d_conformers_but_one_keeps_requested_id_and_resets_to_zero() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C));
        builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(BondSpec::new(
                crate::AtomId::new(0),
                crate::AtomId::new(1),
                BondOrder::Single,
            ))
            .unwrap();
        builder
            .set_2d_coordinates(vec![[0.0, 0.0], [1.0, 0.0]])
            .unwrap();
        builder
            .add_2d_conformer(vec![[0.0, 1.0], [1.0, 1.0]])
            .unwrap();
        builder
            .add_2d_conformer(vec![[0.0, 2.0], [1.0, 2.0]])
            .unwrap();
        let mut mol = builder.build().unwrap();

        rdkit_remove_all_2d_conformers_but_one(&mut mol, 2).unwrap();

        assert_eq!(mol.conformers_2d().len(), 1);
        assert_eq!(mol.conformers_2d()[0].id(), 0);
        assert_eq!(
            mol.conformers_2d()[0].coordinates(),
            &[[0.0, 2.0], [1.0, 2.0]]
        );
    }

    #[test]
    fn remove_all_2d_conformers_but_one_minus_one_keeps_first_conformer() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C));
        builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(BondSpec::new(
                crate::AtomId::new(0),
                crate::AtomId::new(1),
                BondOrder::Single,
            ))
            .unwrap();
        builder
            .set_2d_coordinates(vec![[0.0, 0.0], [1.0, 0.0]])
            .unwrap();
        builder
            .add_2d_conformer(vec![[0.0, 1.0], [1.0, 1.0]])
            .unwrap();
        let mut mol = builder.build().unwrap();

        rdkit_remove_all_2d_conformers_but_one(&mut mol, -1).unwrap();

        assert_eq!(mol.conformers_2d().len(), 1);
        assert_eq!(mol.conformers_2d()[0].id(), 0);
        assert_eq!(
            mol.conformers_2d()[0].coordinates(),
            &[[0.0, 0.0], [1.0, 0.0]]
        );
    }

    #[test]
    fn constrained_depiction_align_only_aligns_existing_conformer_and_keeps_one() {
        let mut ref_builder = MoleculeBuilder::new();
        let r0 = ref_builder.add_atom(AtomSpec::new(Element::C));
        let r1 = ref_builder.add_atom(AtomSpec::new(Element::C));
        ref_builder
            .add_bond(BondSpec::new(r0, r1, BondOrder::Single))
            .unwrap();
        ref_builder
            .set_2d_coordinates(vec![[0.0, 0.0], [2.0, 0.0]])
            .unwrap();
        let reference = ref_builder.build().unwrap();

        let mut mol_builder = MoleculeBuilder::new();
        let m0 = mol_builder.add_atom(AtomSpec::new(Element::C));
        let m1 = mol_builder.add_atom(AtomSpec::new(Element::C));
        mol_builder
            .add_bond(BondSpec::new(m0, m1, BondOrder::Single))
            .unwrap();
        mol_builder
            .set_2d_coordinates(vec![[0.0, 0.0], [0.0, 2.0]])
            .unwrap();
        mol_builder
            .add_2d_conformer(vec![[10.0, 10.0], [11.0, 10.0]])
            .unwrap();
        let mut mol = mol_builder.build().unwrap();

        let params = ConstrainedDepictionParams {
            align_only: true,
            existing_conf_id: 0,
            ..ConstrainedDepictionParams::default()
        };
        generate_depiction_matching_2d_structure_with_ref_match(
            &mut mol,
            &reference,
            &[(0, 0), (1, 1)],
            0,
            &params,
        )
        .unwrap();

        assert_eq!(mol.conformers_2d().len(), 1);
        assert_eq!(mol.conformers_2d()[0].id(), 0);
        let coords = mol.conformers_2d()[0].coordinates();
        assert!((coords[0][0] - 0.0).abs() < 1.0e-6);
        assert!((coords[0][1] - 0.0).abs() < 1.0e-6);
        assert!((coords[1][0] - 2.0).abs() < 1.0e-6);
        assert!((coords[1][1] - 0.0).abs() < 1.0e-6);
    }

    #[test]
    fn constrained_depiction_hard_constraint_uses_reference_core_coords_and_keeps_new_conformer() {
        let _lock = prefer_coordgen_test_lock().lock().unwrap();
        let _guard = PreferCoordGenGuard::capture();
        set_prefer_coord_gen(false);
        let mut ref_builder = MoleculeBuilder::new();
        let r0 = ref_builder.add_atom(AtomSpec::new(Element::C));
        let r1 = ref_builder.add_atom(AtomSpec::new(Element::C));
        ref_builder
            .add_bond(BondSpec::new(r0, r1, BondOrder::Single))
            .unwrap();
        ref_builder
            .set_2d_coordinates(vec![[1.0, 1.0], [3.0, 1.0]])
            .unwrap();
        let reference = ref_builder.build().unwrap();

        let mut mol_builder = MoleculeBuilder::new();
        let m0 = mol_builder.add_atom(AtomSpec::new(Element::C));
        let m1 = mol_builder.add_atom(AtomSpec::new(Element::C));
        let m2 = mol_builder.add_atom(AtomSpec::new(Element::O));
        mol_builder
            .add_bond(BondSpec::new(m0, m1, BondOrder::Single))
            .unwrap();
        mol_builder
            .add_bond(BondSpec::new(m1, m2, BondOrder::Single))
            .unwrap();
        mol_builder
            .set_2d_coordinates(vec![[9.0, 9.0], [10.0, 9.0], [11.0, 9.0]])
            .unwrap();
        let mut mol = mol_builder.build().unwrap();

        generate_depiction_matching_2d_structure_with_ref_match(
            &mut mol,
            &reference,
            &[(0, 0), (1, 1)],
            0,
            &ConstrainedDepictionParams::default(),
        )
        .unwrap();

        assert_eq!(mol.conformers_2d().len(), 1);
        let coords = mol.conformers_2d()[0].coordinates();
        assert!((coords[0][0] - 1.0).abs() < 1.0e-6);
        assert!((coords[0][1] - 1.0).abs() < 1.0e-6);
        assert!((coords[1][0] - 3.0).abs() < 1.0e-6);
        assert!((coords[1][1] - 1.0).abs() < 1.0e-6);
    }

    #[test]
    fn constrained_depiction_outer_accept_failure_generates_unconstrained_coords() {
        let _lock = prefer_coordgen_test_lock().lock().unwrap();
        let _guard = PreferCoordGenGuard::capture();
        set_prefer_coord_gen(false);
        let mut ref_builder = MoleculeBuilder::new();
        let r0 = ref_builder.add_atom(AtomSpec::new(Element::C));
        ref_builder.set_2d_coordinates(vec![[0.0, 0.0]]).unwrap();
        let reference = ref_builder.build().unwrap();

        let mut mol_builder = MoleculeBuilder::new();
        let m0 = mol_builder.add_atom(AtomSpec::new(Element::O));
        let m1 = mol_builder.add_atom(AtomSpec::new(Element::O));
        mol_builder
            .add_bond(BondSpec::new(m0, m1, BondOrder::Single))
            .unwrap();
        let mut mol = mol_builder.build().unwrap();

        let match_vect = generate_depiction_matching_2d_structure(
            &mut mol,
            &reference,
            0,
            None,
            &ConstrainedDepictionParams {
                accept_failure: true,
                ..ConstrainedDepictionParams::default()
            },
        )
        .unwrap();

        assert!(match_vect.is_empty());
        assert_eq!(mol.conformers_2d().len(), 1);
        assert_eq!(mol.conformers_2d()[0].coordinates().len(), 2);
    }

    #[test]
    fn constrained_depiction_outer_reference_pattern_returns_reference_atom_mapping() {
        let _lock = prefer_coordgen_test_lock().lock().unwrap();
        let _guard = PreferCoordGenGuard::capture();
        set_prefer_coord_gen(false);
        let mut ref_builder = MoleculeBuilder::new();
        let r0 = ref_builder.add_atom(AtomSpec::new(Element::C));
        let r1 = ref_builder.add_atom(AtomSpec::new(Element::C));
        let r2 = ref_builder.add_atom(AtomSpec::new(Element::O));
        ref_builder
            .add_bond(BondSpec::new(r0, r1, BondOrder::Single))
            .unwrap();
        ref_builder
            .add_bond(BondSpec::new(r1, r2, BondOrder::Single))
            .unwrap();
        ref_builder
            .set_2d_coordinates(vec![[0.0, 0.0], [1.5, 0.0], [3.0, 0.0]])
            .unwrap();
        let reference = ref_builder.build().unwrap();

        let mut pattern_builder = MoleculeBuilder::new();
        let p0 = pattern_builder.add_atom(AtomSpec::new(Element::C));
        let p1 = pattern_builder.add_atom(AtomSpec::new(Element::O));
        pattern_builder
            .add_bond(BondSpec::new(p0, p1, BondOrder::Single))
            .unwrap();
        let reference_pattern = pattern_builder.build().unwrap();

        let mut mol_builder = MoleculeBuilder::new();
        let m0 = mol_builder.add_atom(AtomSpec::new(Element::C));
        let m1 = mol_builder.add_atom(AtomSpec::new(Element::C));
        let m2 = mol_builder.add_atom(AtomSpec::new(Element::O));
        mol_builder
            .add_bond(BondSpec::new(m0, m1, BondOrder::Single))
            .unwrap();
        mol_builder
            .add_bond(BondSpec::new(m1, m2, BondOrder::Single))
            .unwrap();
        let mut mol = mol_builder.build().unwrap();

        let match_vect = generate_depiction_matching_2d_structure(
            &mut mol,
            &reference,
            0,
            Some(&reference_pattern),
            &ConstrainedDepictionParams::default(),
        )
        .unwrap();

        assert_eq!(match_vect, vec![(1, 1), (2, 2)]);
        assert_eq!(mol.conformers_2d().len(), 1);
    }

    #[test]
    fn constrained_depiction_matching_3d_structure_accept_failure_falls_back_to_2d_coords() {
        let _lock = prefer_coordgen_test_lock().lock().unwrap();
        let _guard = PreferCoordGenGuard::capture();
        set_prefer_coord_gen(false);
        let mut ref_builder = MoleculeBuilder::new();
        ref_builder.add_atom(AtomSpec::new(Element::C));
        let reference = ref_builder.build().unwrap();

        let mut mol_builder = MoleculeBuilder::new();
        let m0 = mol_builder.add_atom(AtomSpec::new(Element::C));
        let m1 = mol_builder.add_atom(AtomSpec::new(Element::C));
        mol_builder
            .add_bond(BondSpec::new(m0, m1, BondOrder::Single))
            .unwrap();
        let mut mol = mol_builder.build().unwrap();

        generate_depiction_matching_3d_structure(&mut mol, &reference, -1, None, true, false)
            .unwrap();

        assert_eq!(mol.conformers_2d().len(), 1);
        assert_eq!(mol.conformers_2d()[0].coordinates().len(), 2);
    }

    #[test]
    fn constrained_depiction_matching_3d_structure_matches_mimic_distmat_path() {
        let _lock = prefer_coordgen_test_lock().lock().unwrap();
        let _guard = PreferCoordGenGuard::capture();
        set_prefer_coord_gen(false);
        let mut ref_builder = MoleculeBuilder::new();
        let r0 = ref_builder.add_atom(AtomSpec::new(Element::C));
        let r1 = ref_builder.add_atom(AtomSpec::new(Element::C));
        ref_builder
            .add_bond(BondSpec::new(r0, r1, BondOrder::Single))
            .unwrap();
        ref_builder
            .add_3d_conformer(vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
            .unwrap();
        let reference = ref_builder.build().unwrap();

        let mut mol_builder = MoleculeBuilder::new();
        let m0 = mol_builder.add_atom(AtomSpec::new(Element::C));
        let m1 = mol_builder.add_atom(AtomSpec::new(Element::C));
        mol_builder
            .add_bond(BondSpec::new(m0, m1, BondOrder::Single))
            .unwrap();
        let mut mol = mol_builder.build().unwrap();

        let expected = compute_2d_coords_mimic_distmat_with_params(
            mol.atoms(),
            mol.bonds(),
            Some(&[2.0]),
            &Compute2DCoordsMimicDistMatParameters {
                canon_orient: false,
                clear_confs: true,
                weight_dist_mat: 0.5,
                n_flips_per_sample: 3,
                n_samples: 100,
                sample_seed: 25,
                permute_deg4_nodes: true,
                force_rdkit: false,
            },
        )
        .unwrap();

        generate_depiction_matching_3d_structure(&mut mol, &reference, 0, None, false, false)
            .unwrap();

        assert_eq!(mol.conformers_2d().len(), 1);
        assert_eq!(mol.conformers_2d()[0].coordinates(), expected.as_slice());
    }

    #[test]
    fn straighten_depiction_minimize_rotation_keeps_smallest_adjustment() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .unwrap();
        builder
            .add_2d_conformer(vec![[0.0, 0.0], [1.0, 0.1]])
            .unwrap();
        let mut mol = builder.build().unwrap();

        straighten_depiction(&mut mol, 0, true).unwrap();

        let coords = mol.conformers_2d()[0].coordinates();
        let dx = coords[1][0] - coords[0][0];
        let dy = coords[1][1] - coords[0][1];
        let theta = dy.atan2(dx).to_degrees();
        assert!(theta.abs() < 1.0);
    }

    #[test]
    fn normalize_depiction_canonicalize_zero_centers_without_rotating() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .unwrap();
        builder
            .add_2d_conformer(vec![[1.0, 1.0], [2.0, 1.0]])
            .unwrap();
        let mut mol = builder.build().unwrap();

        let scale = normalize_depiction(&mut mol, 0, 0, 1.0).unwrap();

        assert!((scale - 1.0).abs() < 1.0e-12);
        let coords = mol.conformers_2d()[0].coordinates();
        assert!((coords[0][0] + 0.5).abs() < 1.0e-8);
        assert!((coords[1][0] - 0.5).abs() < 1.0e-8);
        assert!(coords[0][1].abs() < 1.0e-8);
        assert!(coords[1][1].abs() < 1.0e-8);
    }

    #[test]
    fn normalize_depiction_negative_canonicalize_rotates_ninety_degrees() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .unwrap();
        builder
            .add_2d_conformer(vec![[0.0, 0.0], [2.0, 0.0]])
            .unwrap();
        let mut mol = builder.build().unwrap();

        normalize_depiction(&mut mol, 0, -1, 1.0).unwrap();

        let coords = mol.conformers_2d()[0].coordinates();
        assert!(coords[0][0].abs() < 1.0e-8);
        assert!(coords[1][0].abs() < 1.0e-8);
        assert!((coords[0][1] + 1.0).abs() < 1.0e-8);
        assert!((coords[1][1] - 1.0).abs() < 1.0e-8);
    }

    #[test]
    fn normalize_depiction_negative_scale_uses_most_common_bond_length() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        let a2 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(a1, a2, BondOrder::Single))
            .unwrap();
        builder
            .add_2d_conformer(vec![[0.0, 0.0], [2.0, 0.0], [4.0, 0.0]])
            .unwrap();
        let mut mol = builder.build().unwrap();

        let scale = normalize_depiction(&mut mol, 0, 0, -1.0).unwrap();

        assert!((scale - 0.75).abs() < 1.0e-12);
        let coords = mol.conformers_2d()[0].coordinates();
        let dx = coords[1][0] - coords[0][0];
        assert!((dx.abs() - 1.5).abs() < 1.0e-8);
    }

    #[test]
    fn embedded_frag_ctor_single_atom_updates_pending_neighbors() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::O));
        let a2 = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(a0, a2, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();
        let adjacency = vec![vec![1, 2], vec![0], vec![0]];
        let degree = vec![2usize, 1, 1];
        let frag = RdkitEmbeddedFrag::from_single_atom(
            0,
            mol.atoms(),
            mol.bonds(),
            &adjacency,
            &degree,
            &[],
        );
        assert_eq!(frag.attach_pts.iter().copied().collect::<Vec<_>>(), vec![0]);
        assert_eq!(frag.eatoms[&0].pending, vec![1, 2]);
    }

    #[test]
    fn embedded_frag_coord_map_setup_populates_attachment_points_and_normals() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        let a2 = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(a1, a2, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();
        let adjacency = vec![vec![1], vec![0, 2], vec![1]];
        let degree = vec![1usize, 2, 1];
        let atom_ring_counts = vec![0usize; mol.num_atoms()];
        let coords = BTreeMap::from([(0usize, (0.0, 0.0)), (1usize, (1.5, 0.0))]);
        let frag = RdkitEmbeddedFrag::from_coord_map(
            &coords,
            mol.atoms(),
            mol.bonds(),
            &adjacency,
            &degree,
            &[],
            &atom_ring_counts,
        );
        assert_eq!(frag.attach_pts.iter().copied().collect::<Vec<_>>(), vec![1]);
        assert_eq!(frag.eatoms[&1].pending, vec![2]);
        assert_eq!(frag.eatoms[&1].nbr1, Some(0));
        assert!((frag.eatoms[&1].normal.0 - 0.0).abs() < 1e-8);
        assert!((frag.eatoms[&1].normal.1 + 1.0).abs() < 1e-8);
    }

    #[test]
    fn embedded_frag_double_bond_ctor_seeds_cis_trans_fragment() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        let a2 = builder.add_atom(AtomSpec::new(Element::F));
        let a3 = builder.add_atom(AtomSpec::new(Element::CL));
        builder
            .add_bond(
                BondSpec::new(a0, a1, BondOrder::Double)
                    .with_stereo(BondStereo::Cis)
                    .with_stereo_atoms(a2, a3),
            )
            .unwrap();
        let mol = builder.build().unwrap();
        let frag = RdkitEmbeddedFrag::from_double_bond(&mol.bonds()[0]).unwrap();
        assert_eq!(frag.eatoms[&0].nbr1, Some(1));
        assert_eq!(frag.eatoms[&1].nbr1, Some(0));
        assert_eq!(frag.eatoms[&0].cis_trans_nbr, Some(2));
        assert_eq!(frag.eatoms[&1].cis_trans_nbr, Some(3));
        assert_eq!(frag.eatoms[&1].normal, (0.0, -1.0));
        assert!(frag.eatoms[&1].ccw);
    }

    #[test]
    fn compute_nbrs_and_attachment_points_set_angle_neighbors_and_find_neighbor() {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..4)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect();
        for (a, b) in [(0, 1), (0, 2), (0, 3)] {
            builder
                .add_bond(BondSpec::new(atoms[a], atoms[b], BondOrder::Single))
                .unwrap();
        }
        let mol = builder.build().unwrap();
        let adjacency = vec![vec![1, 2, 3], vec![0], vec![0], vec![0]];
        let degree = vec![3usize, 1, 1, 1];
        let atom_ring_counts = vec![0usize; mol.num_atoms()];
        let mut frag = RdkitEmbeddedFrag::from_coord_map(
            &BTreeMap::from([
                (0usize, (0.0, 0.0)),
                (1usize, (1.0, 0.0)),
                (2usize, (0.0, 1.0)),
            ]),
            mol.atoms(),
            mol.bonds(),
            &adjacency,
            &degree,
            &[],
            &atom_ring_counts,
        );
        frag.setup_attachment_points(&adjacency, &atom_ring_counts);
        assert_eq!(frag.find_neighbor(3, &adjacency), Some(0));
        assert_eq!(frag.eatoms[&0].pending, vec![3]);
        assert!(frag.find_num_neigh((0.1, 0.1), 1.1) >= 2);
        assert!(frag.eatoms[&0].nbr1.is_some());
        assert!(frag.eatoms[&0].angle > 0.0);
    }

    #[test]
    fn stereo_template_match_accepts_matching_cis_double_bond_layout() {
        let mut mol_builder = MoleculeBuilder::new();
        let a0 = mol_builder.add_atom(AtomSpec::new(Element::C));
        let a1 = mol_builder.add_atom(AtomSpec::new(Element::C));
        let a2 = mol_builder.add_atom(AtomSpec::new(Element::F));
        let a3 = mol_builder.add_atom(AtomSpec::new(Element::CL));
        mol_builder
            .add_bond(
                BondSpec::new(a0, a1, BondOrder::Double)
                    .with_stereo(BondStereo::Cis)
                    .with_stereo_atoms(a2, a3),
            )
            .unwrap();
        mol_builder
            .add_bond(BondSpec::new(a0, a2, BondOrder::Single))
            .unwrap();
        mol_builder
            .add_bond(BondSpec::new(a1, a3, BondOrder::Single))
            .unwrap();
        let mol = mol_builder.build().unwrap();

        let template =
            build_rdkit_template_runtime_model("[#6](/[F])=[#6](\\Cl) |(0,0,;1,0,;0,1,;1,1,)|")
                .unwrap();
        let template_mol = build_rdkit_template_query_molecule(&template).unwrap();
        let match_result = get_substruct_match(&mol, &template_mol).expect("template should match");
        assert!(rdkit_check_stereo_chemistry(
            &mol,
            &template,
            &match_result.atom_mapping
        ));
    }

    #[test]
    fn stereo_template_match_rejects_mismatched_double_bond_layout() {
        let mut mol_builder = MoleculeBuilder::new();
        let a0 = mol_builder.add_atom(AtomSpec::new(Element::C));
        let a1 = mol_builder.add_atom(AtomSpec::new(Element::C));
        let a2 = mol_builder.add_atom(AtomSpec::new(Element::F));
        let a3 = mol_builder.add_atom(AtomSpec::new(Element::CL));
        mol_builder
            .add_bond(
                BondSpec::new(a0, a1, BondOrder::Double)
                    .with_stereo(BondStereo::Trans)
                    .with_stereo_atoms(a2, a3),
            )
            .unwrap();
        mol_builder
            .add_bond(BondSpec::new(a0, a2, BondOrder::Single))
            .unwrap();
        mol_builder
            .add_bond(BondSpec::new(a1, a3, BondOrder::Single))
            .unwrap();
        let mol = mol_builder.build().unwrap();

        let template =
            build_rdkit_template_runtime_model("[#6](/[F])=[#6](\\Cl) |(0,0,;1,0,;0,1,;1,1,)|")
                .unwrap();
        let template_mol = build_rdkit_template_query_molecule(&template).unwrap();
        let match_result =
            get_substruct_match(&mol, &template_mol).expect("template should graph-match");
        assert!(!rdkit_check_stereo_chemistry(
            &mol,
            &template,
            &match_result.atom_mapping
        ));
    }

    #[test]
    fn match_to_template_returns_false_when_size_bucket_is_missing() {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..3)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect();
        for (a, b) in [(0, 1), (1, 2), (2, 0)] {
            builder
                .add_bond(BondSpec::new(atoms[a], atoms[b], BondOrder::Single))
                .unwrap();
        }
        let mol = builder.build().unwrap();
        let adjacency = vec![vec![1, 2], vec![0, 2], vec![1, 0]];
        let degree = vec![2usize, 2, 2];
        let atom_ring_counts = vec![1usize; mol.num_atoms()];
        let mut frag = RdkitEmbeddedFrag::default();
        let result = frag
            .match_to_template(
                mol.atoms(),
                mol.bonds(),
                &adjacency,
                &degree,
                &atom_ring_counts,
                &[0, 1, 2],
                1,
            )
            .expect("helper should not error");
        let bucket_exists = rdkit_has_template_of_size(3);
        assert_eq!(result, bucket_exists && !frag.eatoms.is_empty());
    }

    #[test]
    fn mirror_trans_ring_atoms_reflects_trans_ring_atom_into_ring() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        let a2 = builder.add_atom(AtomSpec::new(Element::C));
        let a3 = builder.add_atom(AtomSpec::new(Element::C));
        let a4 = builder.add_atom(AtomSpec::new(Element::F));
        let a5 = builder.add_atom(AtomSpec::new(Element::CL));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(
                BondSpec::new(a1, a2, BondOrder::Double)
                    .with_stereo(BondStereo::Trans)
                    .with_stereo_atoms(a0, a3),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(a2, a3, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(a3, a0, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(a1, a4, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(a2, a5, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();

        let ring = vec![0usize, 1, 2, 3];
        let mut coords = BTreeMap::from([
            (0usize, (0.0, 0.0)),
            (1usize, (1.0, 1.0)),
            (2usize, (2.0, 0.0)),
            (3usize, (1.0, 0.0)),
        ]);
        RdkitEmbeddedFrag::mirror_trans_ring_atoms(mol.bonds(), &ring, &mut coords);
        assert!((coords[&1].0 - 1.0).abs() < 1e-8);
        assert!((coords[&1].1 + 1.0).abs() < 1e-8);
    }

    #[test]
    fn init_ring_coords_sets_prev_and_next_neighbors() {
        let mut frag = RdkitEmbeddedFrag::default();
        let ring = vec![10usize, 11, 12, 13];
        let coords = rdkit_embed_ring(&ring);
        frag.init_from_ring_coords(&ring, &coords);
        let expected = PI * (1.0 - (2.0 / ring.len() as f64));
        assert_eq!(frag.eatoms[&10].nbr1, Some(13));
        assert_eq!(frag.eatoms[&10].nbr2, Some(11));
        assert_eq!(frag.eatoms[&11].nbr1, Some(10));
        assert_eq!(frag.eatoms[&11].nbr2, Some(12));
        assert!((frag.eatoms[&12].angle - expected).abs() < 1e-8);
        assert_eq!(frag.eatoms[&13].nbr2, Some(10));
    }

    #[test]
    fn merge_ring_updates_pin_atom_neighbor_and_angle() {
        let mut master = RdkitEmbeddedFrag::default();
        master.eatoms.insert(
            0,
            TreeEmbeddedAtom {
                loc: (0.0, 0.0),
                normal: (0.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: 1.0,
                nbr1: Some(3),
                nbr2: Some(1),
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        master.eatoms.insert(
            1,
            TreeEmbeddedAtom {
                loc: (1.0, 0.0),
                normal: (0.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: 1.5,
                nbr1: Some(0),
                nbr2: Some(2),
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );

        let mut emb_ring = RdkitEmbeddedFrag::default();
        emb_ring.eatoms.insert(
            0,
            TreeEmbeddedAtom {
                loc: (0.0, 0.0),
                normal: (0.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: 2.0,
                nbr1: Some(3),
                nbr2: Some(4),
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        emb_ring.eatoms.insert(
            4,
            TreeEmbeddedAtom {
                loc: (-1.0, 0.0),
                normal: (0.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: 2.5,
                nbr1: Some(0),
                nbr2: Some(5),
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );

        master.merge_ring(&emb_ring, 2, &[0]);
        assert!((master.eatoms[&0].angle - 3.0).abs() < 1e-8);
        assert_eq!(master.eatoms[&0].nbr1, Some(4));
        assert!(master.eatoms.contains_key(&4));
    }

    #[test]
    fn embed_fused_rings_seeds_first_ring_when_templates_disabled() {
        let _lock = template_registry_test_lock()
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        let _guard = TemplateRegistrySnapshotGuard::capture();

        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect();
        for (a, b) in [(0, 1), (1, 2), (2, 3), (3, 0), (2, 4), (4, 5), (5, 3)] {
            builder
                .add_bond(BondSpec::new(atoms[a], atoms[b], BondOrder::Single))
                .unwrap();
        }
        let mol = builder.build().unwrap();
        let adjacency = crate::AdjacencyList::from_topology(mol.num_atoms(), mol.bonds());
        let adjacency: Vec<Vec<usize>> = (0..mol.num_atoms())
            .map(|i| {
                adjacency
                    .neighbors_of(i)
                    .iter()
                    .map(|n| n.atom_index)
                    .collect()
            })
            .collect();
        let degree: Vec<usize> = adjacency.iter().map(Vec::len).collect();
        let atom_ring_counts = vec![1usize; mol.num_atoms()];
        let fused_rings = vec![vec![0usize, 1, 2, 3], vec![2usize, 4, 5, 3]];

        let mut frag = RdkitEmbeddedFrag::default();
        frag.embed_fused_rings(
            mol.atoms(),
            mol.bonds(),
            &adjacency,
            &degree,
            &[],
            &atom_ring_counts,
            &fused_rings,
            false,
        )
        .expect("embedding should succeed");

        let expected_first = rdkit_pick_first_ring_to_embed(&degree, &fused_rings);
        for aid in &fused_rings[expected_first] {
            assert!(frag.eatoms.contains_key(aid));
        }
        assert_eq!(frag.eatoms.len(), 6);
    }

    #[test]
    fn embed_fused_rings_whole_system_template_populates_exact_template_coords() {
        let _lock = template_registry_test_lock()
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        let _guard = TemplateRegistrySnapshotGuard::capture();

        let file = NamedTempFile::new().unwrap();
        std::fs::write(
            file.path(),
            "[*]1[*][*][*]2[*][*][*]12 |(0,0,;1,0,;2,0,;2,1,;1.5,2,;0.5,2,;0,1,)|\n",
        )
        .unwrap();
        set_rdkit_ring_system_templates(file.path().to_str().unwrap()).unwrap();

        let mut builder = MoleculeBuilder::new();
        let atoms = vec![
            builder.add_atom(AtomSpec::new(Element::C)),
            builder.add_atom(AtomSpec::new(Element::O)),
            builder.add_atom(AtomSpec::new(Element::C)),
            builder.add_atom(AtomSpec::new(Element::C)),
            builder.add_atom(AtomSpec::new(Element::C)),
            builder.add_atom(AtomSpec::new(Element::C)),
            builder.add_atom(AtomSpec::new(Element::C)),
        ];
        for (a, b) in [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 6),
            (6, 0),
            (3, 4),
            (4, 5),
            (5, 6),
        ] {
            builder
                .add_bond(BondSpec::new(atoms[a], atoms[b], BondOrder::Single))
                .unwrap();
        }
        let mol = builder.build().unwrap();
        let adjacency = crate::AdjacencyList::from_topology(mol.num_atoms(), mol.bonds());
        let adjacency: Vec<Vec<usize>> = (0..mol.num_atoms())
            .map(|i| {
                adjacency
                    .neighbors_of(i)
                    .iter()
                    .map(|n| n.atom_index)
                    .collect()
            })
            .collect();
        let degree: Vec<usize> = adjacency.iter().map(Vec::len).collect();
        let atom_ring_counts = vec![1usize; mol.num_atoms()];
        let fused_rings = vec![vec![0usize, 1, 2, 3, 6], vec![3usize, 4, 5, 6]];

        let mut frag = RdkitEmbeddedFrag::default();
        frag.embed_fused_rings(
            mol.atoms(),
            mol.bonds(),
            &adjacency,
            &degree,
            &[],
            &atom_ring_counts,
            &fused_rings,
            true,
        )
        .expect("template embedding should succeed");

        let template = rdkit_matching_templates(7)
            .pop()
            .expect("template bucket should exist");
        let coords = template
            .coords_2d
            .expect("template should carry coordinates");
        let mut got_coords: Vec<(i64, i64)> = frag
            .eatoms
            .values()
            .map(|st| {
                (
                    (st.loc.0 * 1_000_000.0).round() as i64,
                    (st.loc.1 * 1_000_000.0).round() as i64,
                )
            })
            .collect();
        let mut expected_coords: Vec<(i64, i64)> = coords
            .iter()
            .map(|coord| {
                (
                    (coord[0] * 1_000_000.0).round() as i64,
                    (coord[1] * 1_000_000.0).round() as i64,
                )
            })
            .collect();
        got_coords.sort_unstable();
        expected_coords.sort_unstable();
        assert_eq!(got_coords, expected_coords);
    }

    #[test]
    fn embed_fused_rings_core_template_fallback_populates_core_coords_before_side_ring_attachment()
    {
        let _lock = template_registry_test_lock()
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        let _guard = TemplateRegistrySnapshotGuard::capture();
        let mol = Molecule::from_smiles("C1C2CC3C1CC1(C2)NC31").unwrap();
        let adjacency = crate::AdjacencyList::from_topology(mol.num_atoms(), mol.bonds());
        let adjacency: Vec<Vec<usize>> = (0..mol.num_atoms())
            .map(|i| {
                adjacency
                    .neighbors_of(i)
                    .iter()
                    .map(|n| n.atom_index)
                    .collect()
            })
            .collect();
        let degree: Vec<usize> = adjacency.iter().map(Vec::len).collect();
        let ring_info = crate::symmetrize_sssr(&mol).unwrap();
        let atom_ring_counts: Vec<usize> = (0..mol.num_atoms())
            .map(|aid| ring_info.num_atom_rings(AtomId::new(aid)))
            .collect();
        let fused_rings: Vec<Vec<usize>> = ring_info
            .atom_rings()
            .iter()
            .map(|ring| ring.iter().map(|aid| aid.index()).collect())
            .collect();
        let (core_rings, core_ring_ids) = rdkit_find_core_rings(&fused_rings, mol.bonds());
        assert!(
            core_rings.len() > 1 && core_rings.len() < fused_rings.len(),
            "fixture must exercise core-template fallback: fused_rings={fused_rings:?} core_rings={core_rings:?}"
        );
        let core_union: Vec<usize> = core_rings
            .iter()
            .flatten()
            .copied()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let core_mol =
            build_rdkit_ring_system_molecule(mol.atoms(), mol.bonds(), &core_union).unwrap();
        let template = build_rdkit_template_runtime_model(
            "C1C2CC3CC1CC3C2 |(-7.01,3.13,;-7.71,4.35,;-7.01,5.56,;-5.61,5.56,;-4.91,4.35,;-5.61,3.13,;-4.28,3.57,;-4.28,5.13,;-6.34,4.05,)|",
        )
        .unwrap();
        rdkit_coordinate_template_registry()
            .write()
            .expect("template registry lock poisoned")
            .templates = HashMap::from([(core_union.len(), vec![template.clone()])]);
        let template_mol = build_rdkit_template_query_molecule(&template).unwrap();
        assert_eq!(core_ring_ids.len(), core_rings.len());
        assert_eq!(core_mol.num_bonds(), template_mol.num_bonds());
        assert_eq!(
            crate::symmetrize_sssr(&core_mol).unwrap().num_rings(),
            core_rings.len()
        );
        assert_eq!(
            crate::symmetrize_sssr(&template_mol).unwrap().num_rings(),
            core_rings.len()
        );
        assert_eq!(
            rdkit_template_degree_counts(&core_mol, None),
            rdkit_template_degree_counts(&template_mol, None)
        );
        assert!(
            get_substruct_match(&core_mol, &template_mol).is_some(),
            "core ring union should substructure-match the explicit core template"
        );

        let mut core_frag = RdkitEmbeddedFrag::default();
        assert!(
            core_frag
                .match_to_template(
                    mol.atoms(),
                    mol.bonds(),
                    &adjacency,
                    &degree,
                    &atom_ring_counts,
                    &core_union,
                    core_rings.len(),
                )
                .expect("core match helper should not error"),
            "core ring union should template-match before side-ring attachment"
        );

        let mut frag = RdkitEmbeddedFrag::default();
        frag.embed_fused_rings(
            mol.atoms(),
            mol.bonds(),
            &adjacency,
            &degree,
            &[],
            &atom_ring_counts,
            &fused_rings,
            true,
        )
        .expect("core-template fallback embedding should succeed");

        assert_eq!(frag.eatoms.len(), mol.num_atoms());
        let template = rdkit_matching_templates(core_union.len())
            .pop()
            .expect("template bucket should exist");
        let coords = template
            .coords_2d
            .expect("template should carry coordinates");
        let mut got_core_coords: Vec<(i64, i64)> = core_union
            .iter()
            .map(|aid| {
                let st = frag.eatoms.get(aid).expect("core atom should be embedded");
                (
                    (st.loc.0 * 1_000_000.0).round() as i64,
                    (st.loc.1 * 1_000_000.0).round() as i64,
                )
            })
            .collect();
        let mut expected_coords: Vec<(i64, i64)> = coords
            .iter()
            .map(|coord| {
                (
                    (coord[0] * 1_000_000.0).round() as i64,
                    (coord[1] * 1_000_000.0).round() as i64,
                )
            })
            .collect();
        got_core_coords.sort_unstable();
        expected_coords.sort_unstable();
        assert_eq!(got_core_coords, expected_coords);
    }

    #[test]
    fn bisect_point_flips_midpoint_for_obtuse_angles() {
        let acute = rdkit_compute_bisect_point((0.0, 0.0), PI / 2.0, (1.0, 0.0), (0.0, 1.0));
        assert!((acute.0 - 0.5).abs() < 1e-8);
        assert!((acute.1 - 0.5).abs() < 1e-8);
        let obtuse = rdkit_compute_bisect_point((0.0, 0.0), 1.5 * PI, (1.0, 0.0), (0.0, 1.0));
        assert!((obtuse.0 + 0.5).abs() < 1e-8);
        assert!((obtuse.1 + 0.5).abs() < 1e-8);
    }

    #[test]
    fn reflect_point_mirrors_across_line() {
        let reflected = rdkit_reflect_point((1.0, 2.0), (0.0, 0.0), (0.0, 3.0));
        assert!((reflected.0 + 1.0).abs() < 1e-8);
        assert!((reflected.1 - 2.0).abs() < 1e-8);
    }

    #[test]
    fn reflect_points_mirrors_all_points_in_map() {
        let mut points = BTreeMap::from([(0usize, (1.0, 2.0)), (1usize, (-2.0, -1.0))]);
        rdkit_reflect_points(&mut points, (0.0, 0.0), (0.0, 3.0));
        assert!((points[&0].0 + 1.0).abs() < 1e-8);
        assert!((points[&0].1 - 2.0).abs() < 1e-8);
        assert!((points[&1].0 - 2.0).abs() < 1e-8);
        assert!((points[&1].1 + 1.0).abs() < 1e-8);
    }

    #[test]
    fn compute_one_atom_trans_aligns_shared_atom_and_neighbor_midpoint_bisector() {
        let mut master = RdkitEmbeddedFrag::default();
        master.eatoms.insert(
            0,
            TreeEmbeddedAtom {
                loc: (0.0, 0.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: 1.5 * PI,
                nbr1: Some(1),
                nbr2: Some(2),
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        master.eatoms.insert(
            1,
            TreeEmbeddedAtom {
                loc: (1.0, 0.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        master.eatoms.insert(
            2,
            TreeEmbeddedAtom {
                loc: (0.0, 1.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );

        let mut other = RdkitEmbeddedFrag::default();
        other.eatoms.insert(
            0,
            TreeEmbeddedAtom {
                loc: (3.0, 2.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: 1.5 * PI,
                nbr1: Some(3),
                nbr2: Some(4),
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        other.eatoms.insert(
            3,
            TreeEmbeddedAtom {
                loc: (4.0, 2.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        other.eatoms.insert(
            4,
            TreeEmbeddedAtom {
                loc: (3.0, 3.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );

        let trans = master.compute_one_atom_trans(0, &other).unwrap();
        other.transform(trans);
        let shared = other.eatoms[&0].loc;
        assert!((shared.0 - 0.0).abs() < 1e-8);
        assert!((shared.1 - 0.0).abs() < 1e-8);

        let original_mid = ((4.0 + 3.0) * 0.5, (2.0 + 3.0) * 0.5);
        let n3 = other.eatoms[&3].loc;
        let n4 = other.eatoms[&4].loc;
        let mid = ((n3.0 + n4.0) * 0.5, (n3.1 + n4.1) * 0.5);
        let original_dx = original_mid.0 - 3.0f64;
        let original_dy = original_mid.1 - 2.0f64;
        let original_dist = (original_dx * original_dx + original_dy * original_dy).sqrt();
        let new_dist = (mid.0 * mid.0 + mid.1 * mid.1).sqrt();
        assert!((new_dist - original_dist).abs() < 1e-8);
    }

    #[test]
    fn add_atom_with_no_ang_first_seed_rotates_normal_and_places_new_atom() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();
        let adjacency = vec![vec![1], vec![0]];
        let degree = vec![1usize, 1];

        let mut frag = RdkitEmbeddedFrag::from_single_atom(
            0,
            mol.atoms(),
            mol.bonds(),
            &adjacency,
            &degree,
            &[],
        );
        frag.add_atom_to_atom_with_no_ang(1, 0, mol.atoms(), mol.bonds(), &degree)
            .expect("first no-angle attachment should succeed");

        assert_eq!(frag.eatoms[&0].nbr1, Some(1));
        assert!(frag.eatoms[&0].angle < 0.0);
        assert_eq!(frag.eatoms[&1].nbr1, Some(0));
        let loc = frag.eatoms[&1].loc;
        assert!((norm(loc) - BOND_LEN).abs() < 1e-8);
        assert!(norm(frag.eatoms[&1].normal) > 0.999999);
    }

    #[test]
    fn add_atom_with_ang_updates_angle_and_sets_new_embedded_atom() {
        let mut frag = RdkitEmbeddedFrag::default();
        frag.eatoms.insert(
            0,
            TreeEmbeddedAtom {
                loc: (0.0, 0.0),
                normal: (0.0, 1.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: PI,
                nbr1: Some(1),
                nbr2: Some(2),
                rot_dir: 0,
                pending: vec![3],
                d_density: -1.0,
                df_fixed: false,
            },
        );
        frag.eatoms.insert(
            1,
            TreeEmbeddedAtom {
                loc: (-1.5, 0.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        frag.eatoms.insert(
            2,
            TreeEmbeddedAtom {
                loc: (1.5, 0.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );

        frag.add_atom_to_atom_with_ang(3, 0)
            .expect("angle-based attachment should succeed");

        assert!(frag.eatoms[&0].angle > PI);
        assert_eq!(frag.eatoms[&0].nbr2, Some(3));
        assert_eq!(frag.eatoms[&3].nbr1, Some(0));
        assert!(norm(frag.eatoms[&3].normal) > 0.999999);
    }

    #[test]
    fn add_non_ring_atom_removes_pending_neighbor_and_updates_new_atom_pending() {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..3)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect();
        builder
            .add_bond(BondSpec::new(atoms[0], atoms[1], BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(atoms[1], atoms[2], BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();
        let adjacency = vec![vec![1], vec![0, 2], vec![1]];
        let degree = vec![1usize, 2, 1];

        let mut frag = RdkitEmbeddedFrag::from_single_atom(
            1,
            mol.atoms(),
            mol.bonds(),
            &adjacency,
            &degree,
            &[],
        );
        assert_eq!(frag.eatoms[&1].pending, vec![0, 2]);
        frag.add_non_ring_atom(0, 1, mol.atoms(), mol.bonds(), &adjacency, &degree, &[])
            .expect("non-ring attachment should succeed");

        assert!(!frag.eatoms[&1].pending.contains(&0));
        assert_eq!(frag.eatoms[&0].nbr1, Some(1));
        assert!(frag.attach_pts.iter().any(|&aid| aid == 0 || aid == 1));
    }

    #[test]
    fn expand_efrag_grows_single_atom_seed_across_pending_chain() {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..3)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect();
        builder
            .add_bond(BondSpec::new(atoms[0], atoms[1], BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(atoms[1], atoms[2], BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();
        let adjacency = vec![vec![1], vec![0, 2], vec![1]];
        let degree = vec![1usize, 2, 1];

        let mut frag = RdkitEmbeddedFrag::from_single_atom(
            0,
            mol.atoms(),
            mol.bonds(),
            &adjacency,
            &degree,
            &[],
        );
        let mut nratms = vec![1usize, 2];
        let mut efrags = Vec::<RdkitEmbeddedFrag>::new();
        frag.expand_efrag(
            &mut nratms,
            &mut efrags,
            mol.atoms(),
            mol.bonds(),
            &adjacency,
            &degree,
            &[],
        )
        .expect("expansion should succeed");

        assert!(nratms.is_empty());
        assert!(efrags.is_empty());
        assert_eq!(frag.eatoms.len(), 3);
        assert!(frag.attach_pts.is_empty());
        assert!(frag.eatoms[&0].pending.is_empty());
        assert!(frag.eatoms[&1].pending.is_empty());
        assert!(frag.eatoms[&2].pending.is_empty());
    }

    #[test]
    fn find_collisions_sets_density_and_reports_close_pair() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();
        let adjacency = vec![vec![1], vec![0]];

        let mut frag = RdkitEmbeddedFrag::default();
        frag.eatoms.insert(
            0,
            TreeEmbeddedAtom {
                loc: (0.0, 0.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        frag.eatoms.insert(
            1,
            TreeEmbeddedAtom {
                loc: (0.1, 0.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );

        let dmat = frag.compute_dist_mat(mol.num_atoms());
        let collisions = frag.find_collisions(mol.atoms(), mol.bonds(), &adjacency, &dmat, false);
        assert_eq!(collisions, vec![(1, 0)]);
        assert!(frag.eatoms[&0].d_density > 0.0);
        assert!(frag.eatoms[&1].d_density > 0.0);
        assert!(frag.total_density() > 0.0);
    }

    #[test]
    fn mimic_distmat_clamps_weight_and_uses_density_when_reference_absent() {
        let mut frag = RdkitEmbeddedFrag::default();
        frag.eatoms.insert(
            0,
            TreeEmbeddedAtom {
                loc: (0.0, 0.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        frag.eatoms.insert(
            1,
            TreeEmbeddedAtom {
                loc: (1.0, 0.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        let no_ref = frag.mimic_dist_mat_and_density_cost_func(2, None, -1.0);
        assert!((no_ref - 1.0).abs() < 1.0e-8);

        let ref_dmat = vec![2.0];
        let mimic_only = frag.mimic_dist_mat_and_density_cost_func(2, Some(&ref_dmat), 2.0);
        assert!((mimic_only - 1.0).abs() < 1.0e-8);
    }

    #[test]
    fn permute_bonds_reflects_both_selected_neighbor_fragments_about_bisector() {
        let adjacency = vec![vec![1, 2, 3, 4], vec![0], vec![0], vec![0], vec![0]];
        let mut frag = RdkitEmbeddedFrag::default();
        for (aid, loc) in [
            (0usize, (0.0, 0.0)),
            (1usize, (-1.0, 0.0)),
            (2usize, (0.0, 1.0)),
            (3usize, (1.0, 0.0)),
            (4usize, (0.0, -1.0)),
        ] {
            frag.eatoms.insert(
                aid,
                TreeEmbeddedAtom {
                    loc,
                    normal: (1.0, 0.0),
                    ccw: true,
                    cis_trans_nbr: None,
                    angle: -1.0,
                    nbr1: None,
                    nbr2: None,
                    rot_dir: 0,
                    pending: Vec::new(),
                    d_density: -1.0,
                    df_fixed: false,
                },
            );
        }

        frag.permute_bonds(0, 1, 2, &adjacency);
        assert!((frag.eatoms[&1].loc.0 - 0.0).abs() < 1.0e-8);
        assert!((frag.eatoms[&1].loc.1 - 1.0).abs() < 1.0e-8);
        assert!((frag.eatoms[&2].loc.0 + 1.0).abs() < 1.0e-8);
        assert!((frag.eatoms[&2].loc.1 - 0.0).abs() < 1.0e-8);
        assert_eq!(frag.eatoms[&3].loc, (1.0, 0.0));
        assert_eq!(frag.eatoms[&4].loc, (0.0, -1.0));
    }

    #[test]
    fn random_sample_flips_and_permutations_keeps_best_sample_when_flip_is_worse() {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..4)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect();
        for (a, b) in [(0, 1), (1, 2), (2, 3)] {
            builder
                .add_bond(BondSpec::new(atoms[a], atoms[b], BondOrder::Single))
                .unwrap();
        }
        let mol = builder.build().unwrap();
        let adjacency = vec![vec![1], vec![0, 2], vec![1, 3], vec![2]];
        let comp = vec![0usize, 1, 2, 3];

        let mut frag = RdkitEmbeddedFrag::default();
        for (aid, loc) in [
            (0usize, (0.0, 0.0)),
            (1usize, (1.0, 0.0)),
            (2usize, (2.0, 0.0)),
            (3usize, (3.0, 1.0)),
        ] {
            frag.eatoms.insert(
                aid,
                TreeEmbeddedAtom {
                    loc,
                    normal: (1.0, 0.0),
                    ccw: true,
                    cis_trans_nbr: None,
                    angle: -1.0,
                    nbr1: None,
                    nbr2: None,
                    rot_dir: 0,
                    pending: Vec::new(),
                    d_density: -1.0,
                    df_fixed: false,
                },
            );
        }

        let before = frag.eatoms[&3].loc;
        frag.random_sample_flips_and_permutations(
            mol.atoms(),
            mol.bonds(),
            &comp,
            &adjacency,
            1,
            1,
            7,
            None,
            0.0,
            false,
        );
        assert_eq!(frag.eatoms[&3].loc, before);
    }

    #[test]
    fn square_planar_nontetrahedral_seed_is_wired_into_compute_2d_pipeline() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(
            AtomSpec::new(Element::PT)
                .with_chiral_tag(ChiralTag::SquarePlanar)
                .with_chiral_permutation(1),
        );
        let a1 = builder.add_atom(AtomSpec::new(Element::CL));
        let a2 = builder.add_atom(AtomSpec::new(Element::CL));
        let a3 = builder.add_atom(AtomSpec::new(Element::CL));
        let a4 = builder.add_atom(AtomSpec::new(Element::CL));
        for nbr in [a1, a2, a3, a4] {
            builder
                .add_bond(BondSpec::new(center, nbr, BondOrder::Single))
                .unwrap();
        }
        let mol = builder.build().unwrap();
        let coords = compute_2d_coords(mol.atoms(), mol.bonds()).unwrap();

        assert!((coords[center.index()][0]).abs() < 1.0e-8);
        assert!((coords[center.index()][1]).abs() < 1.0e-8);
        let expected = ISQRT2 * BOND_LEN;
        assert!((coords[a1.index()][0] - expected).abs() < 1.0e-6);
        assert!((coords[a1.index()][1] - expected).abs() < 1.0e-6);
    }

    #[test]
    fn trigonal_bipyramidal_nontetrahedral_seed_is_wired_into_compute_2d_pipeline() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(
            AtomSpec::new(Element::P)
                .with_chiral_tag(ChiralTag::TrigonalBipyramidal)
                .with_chiral_permutation(1),
        );
        let n1 = builder.add_atom(AtomSpec::new(Element::F));
        let n2 = builder.add_atom(AtomSpec::new(Element::CL));
        let n3 = builder.add_atom(AtomSpec::new(Element::BR));
        let n4 = builder.add_atom(AtomSpec::new(Element::I));
        let n5 = builder.add_atom(AtomSpec::new(Element::O));
        for nbr in [n1, n2, n3, n4, n5] {
            builder
                .add_bond(BondSpec::new(center, nbr, BondOrder::Single))
                .unwrap();
        }
        let mol = builder.build().unwrap();
        let coords = compute_2d_coords(mol.atoms(), mol.bonds()).unwrap();
        let atom_ranks: Vec<i32> = (0..mol.num_atoms())
            .map(|i| {
                let degree = mol
                    .bonds()
                    .iter()
                    .filter(|bond| bond.begin().index() == i || bond.end().index() == i)
                    .count();
                atom_depict_rank(mol.atoms()[i].atomic_number(), degree) as i32
            })
            .collect();
        let expected = embed_tbp(center.index(), mol.atoms(), mol.bonds(), &mol, &atom_ranks);

        let mut actual_pairs = Vec::new();
        let mut expected_pairs = Vec::new();
        let ids = [
            center.index(),
            n1.index(),
            n2.index(),
            n3.index(),
            n4.index(),
            n5.index(),
        ];
        for i in 0..ids.len() {
            for j in 0..i {
                let ai = ids[i];
                let aj = ids[j];
                let adx = coords[ai][0] - coords[aj][0];
                let ady = coords[ai][1] - coords[aj][1];
                actual_pairs.push(((adx * adx + ady * ady).sqrt() * 1_000_000.0).round() as i64);
                let ex = expected[&ai].0 - expected[&aj].0;
                let ey = expected[&ai].1 - expected[&aj].1;
                expected_pairs.push(((ex * ex + ey * ey).sqrt() * 1_000_000.0).round() as i64);
            }
        }
        actual_pairs.sort_unstable();
        expected_pairs.sort_unstable();
        assert_eq!(actual_pairs, expected_pairs);
    }

    #[test]
    fn octahedral_nontetrahedral_seed_is_wired_into_compute_2d_pipeline() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(
            AtomSpec::new(Element::CO)
                .with_chiral_tag(ChiralTag::Octahedral)
                .with_chiral_permutation(1),
        );
        let n1 = builder.add_atom(AtomSpec::new(Element::F));
        let n2 = builder.add_atom(AtomSpec::new(Element::CL));
        let n3 = builder.add_atom(AtomSpec::new(Element::BR));
        let n4 = builder.add_atom(AtomSpec::new(Element::I));
        let n5 = builder.add_atom(AtomSpec::new(Element::N));
        let n6 = builder.add_atom(AtomSpec::new(Element::O));
        for nbr in [n1, n2, n3, n4, n5, n6] {
            builder
                .add_bond(BondSpec::new(center, nbr, BondOrder::Single))
                .unwrap();
        }
        let mol = builder.build().unwrap();
        let coords = compute_2d_coords(mol.atoms(), mol.bonds()).unwrap();

        assert!((coords[center.index()][0]).abs() < 1.0e-8);
        assert!((coords[center.index()][1]).abs() < 1.0e-8);
        let axis_like = [
            n1.index(),
            n2.index(),
            n3.index(),
            n4.index(),
            n5.index(),
            n6.index(),
        ]
        .iter()
        .filter(|&&idx| coords[idx][0].abs() < 1.0e-6 || coords[idx][1].abs() < 1.0e-6)
        .count();
        assert!(axis_like >= 2);
    }

    #[test]
    fn flip_about_bond_reflects_smaller_side_of_non_ring_bond() {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..4)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect();
        builder
            .add_bond(BondSpec::new(atoms[0], atoms[1], BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(atoms[1], atoms[2], BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(atoms[2], atoms[3], BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();
        let adjacency = vec![vec![1], vec![0, 2], vec![1, 3], vec![2]];

        let mut frag = RdkitEmbeddedFrag::default();
        frag.eatoms.insert(
            0,
            TreeEmbeddedAtom {
                loc: (0.0, 0.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        frag.eatoms.insert(
            1,
            TreeEmbeddedAtom {
                loc: (1.0, 0.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        frag.eatoms.insert(
            2,
            TreeEmbeddedAtom {
                loc: (2.0, 0.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        frag.eatoms.insert(
            3,
            TreeEmbeddedAtom {
                loc: (3.0, 1.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );

        frag.flip_about_bond(1, mol.bonds(), &adjacency, |_| false, true);
        assert!((frag.eatoms[&2].loc.0 - 2.0).abs() < 1e-8);
        assert!((frag.eatoms[&2].loc.1 - 0.0).abs() < 1e-8);
        assert!((frag.eatoms[&3].loc.0 - 3.0).abs() < 1e-8);
        assert!((frag.eatoms[&3].loc.1 + 1.0).abs() < 1e-8);
    }

    #[test]
    fn remove_collisions_open_angles_rotates_degree_one_collision_pair() {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..4)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect();
        for (a, b) in [(0, 1), (1, 2), (2, 3)] {
            builder
                .add_bond(BondSpec::new(atoms[a], atoms[b], BondOrder::Single))
                .unwrap();
        }
        let mol = builder.build().unwrap();
        let adjacency = vec![vec![1], vec![0, 2], vec![1, 3], vec![2]];
        let comp = vec![0usize, 1, 2, 3];

        let mut frag = RdkitEmbeddedFrag::default();
        for (aid, loc) in [
            (0usize, (0.0, 0.2)),
            (1usize, (0.0, 0.0)),
            (2usize, (1.0, 0.0)),
            (3usize, (1.0, 0.2)),
        ] {
            frag.eatoms.insert(
                aid,
                TreeEmbeddedAtom {
                    loc,
                    normal: (1.0, 0.0),
                    ccw: true,
                    cis_trans_nbr: None,
                    angle: -1.0,
                    nbr1: None,
                    nbr2: None,
                    rot_dir: 0,
                    pending: Vec::new(),
                    d_density: -1.0,
                    df_fixed: false,
                },
            );
        }

        let before0 = frag.eatoms[&0].loc;
        let before3 = frag.eatoms[&3].loc;
        frag.remove_collisions_open_angles(mol.atoms(), mol.bonds(), &comp, &adjacency);
        assert_ne!(frag.eatoms[&0].loc, before0);
        assert_ne!(frag.eatoms[&3].loc, before3);
    }

    #[test]
    fn remove_collisions_shorten_bonds_moves_terminal_atom_inward_on_open_path() {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..2)
            .map(|_| builder.add_atom(AtomSpec::new(Element::O)))
            .collect();
        builder
            .add_bond(BondSpec::new(atoms[0], atoms[1], BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();
        let adjacency = vec![vec![1], vec![0]];
        let comp = vec![0usize, 1];

        let mut frag = RdkitEmbeddedFrag::default();
        frag.eatoms.insert(
            0,
            TreeEmbeddedAtom {
                loc: (0.0, 0.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        frag.eatoms.insert(
            1,
            TreeEmbeddedAtom {
                loc: (0.9, 0.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );

        frag.remove_collisions_shorten_bonds(mol.atoms(), mol.bonds(), &comp, &adjacency);
        assert!((frag.eatoms[&1].loc.0 - 0.81).abs() < 1e-8);
    }

    #[test]
    fn compute_two_atom_trans_maps_two_shared_atoms_onto_reference_pair() {
        let mut master = RdkitEmbeddedFrag::default();
        master.eatoms.insert(
            0,
            TreeEmbeddedAtom {
                loc: (0.0, 0.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        master.eatoms.insert(
            1,
            TreeEmbeddedAtom {
                loc: (2.0, 0.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );

        let mut other = RdkitEmbeddedFrag::default();
        other.eatoms.insert(
            0,
            TreeEmbeddedAtom {
                loc: (1.0, 1.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );
        other.eatoms.insert(
            1,
            TreeEmbeddedAtom {
                loc: (1.0, 3.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: -1.0,
                df_fixed: false,
            },
        );

        let coords = BTreeMap::from([
            (0usize, other.eatoms[&0].loc),
            (1usize, other.eatoms[&1].loc),
        ]);
        let trans = master.compute_two_atom_trans(0, 1, &coords).unwrap();
        other.transform(trans);
        assert!((other.eatoms[&0].loc.0 - 0.0).abs() < 1e-8);
        assert!((other.eatoms[&0].loc.1 - 0.0).abs() < 1e-8);
        assert!((other.eatoms[&1].loc.0 - 2.0).abs() < 1e-8);
        assert!((other.eatoms[&1].loc.1 - 0.0).abs() < 1e-8);
    }

    #[test]
    fn canonicalize_component_centers_fragment_and_aligns_major_axis() {
        let mut comp = vec![
            (0usize, (1.0, 1.0)),
            (1usize, (2.0, 2.0)),
            (2usize, (3.0, 3.0)),
        ];
        canonicalize_component(&mut comp);
        let sum_x: f64 = comp.iter().map(|(_, p)| p.0).sum();
        let sum_y: f64 = comp.iter().map(|(_, p)| p.1).sum();
        assert!(sum_x.abs() < 1e-8);
        assert!(sum_y.abs() < 1e-8);
        let sum_abs_y: f64 = comp.iter().map(|(_, p)| p.1.abs()).sum();
        assert!(sum_abs_y < 1e-8);
    }

    #[test]
    fn compute2d_disconnected_components_use_component_boxes_for_spacing() {
        let _lock = prefer_coordgen_test_lock().lock().unwrap();
        let _guard = PreferCoordGenGuard::capture();
        set_prefer_coord_gen(false);
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        let a2 = builder.add_atom(AtomSpec::new(Element::C));
        let a3 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(a2, a3, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();
        let coords = compute_2d_coords(mol.atoms(), mol.bonds()).unwrap();

        let first = vec![
            (0usize, (coords[0][0], coords[0][1])),
            (1usize, (coords[1][0], coords[1][1])),
        ];
        let second = vec![
            (2usize, (coords[2][0], coords[2][1])),
            (3usize, (coords[3][0], coords[3][1])),
        ];
        let first_box = component_box_rdkit(&first);
        let first_x_extent = first_box.0 + first_box.1;
        let first_y_extent = first_box.2 + first_box.3;

        if first_x_extent > first_y_extent {
            let first_max_y = first
                .iter()
                .map(|(_, p)| p.1)
                .fold(f64::NEG_INFINITY, f64::max);
            let second_min_y = second
                .iter()
                .map(|(_, p)| p.1)
                .fold(f64::INFINITY, f64::min);
            assert!((second_min_y - first_max_y - 1.0).abs() < 1e-8);
        } else {
            let first_max_x = first
                .iter()
                .map(|(_, p)| p.0)
                .fold(f64::NEG_INFINITY, f64::max);
            let second_min_x = second
                .iter()
                .map(|(_, p)| p.0)
                .fold(f64::INFINITY, f64::min);
            assert!((second_min_x - first_max_x - 1.0).abs() < 1e-8);
        }
    }

    #[test]
    fn component_box_rdkit_matches_embedded_frag_box_sign_convention() {
        let comp = vec![
            (0usize, (-2.0, 1.0)),
            (1usize, (3.0, -4.0)),
            (2usize, (1.0, 2.0)),
        ];
        assert_eq!(component_box_rdkit(&comp), (3.0, 2.0, 2.0, 4.0));
    }

    #[test]
    fn embed_fused_systems_groups_fused_rings_and_respects_coordmap_template_gate() {
        let mol = Molecule::from_smiles("C1C2CC3C1CC1(C2)NC31").unwrap();
        let atoms = mol.atoms();
        let bonds = mol.bonds();
        let n = atoms.len();
        let mut adjacency = vec![Vec::<usize>::new(); n];
        let mut degree = vec![0usize; n];
        for bond in bonds {
            adjacency[bond.begin().index()].push(bond.end().index());
            adjacency[bond.end().index()].push(bond.begin().index());
            degree[bond.begin().index()] += 1;
            degree[bond.end().index()] += 1;
        }
        let cip_ranks = vec![0u32; n];
        let fused_rings: Vec<Vec<usize>> = crate::symmetrize_sssr(&mol)
            .unwrap()
            .atom_rings()
            .iter()
            .map(|ring| ring.iter().map(|aid| aid.index()).collect())
            .collect();
        let mut coord_map = BTreeMap::new();
        coord_map.insert(0usize, (0.0, 0.0));
        coord_map.insert(1usize, (1.0, 0.0));
        let mut efrags = Vec::new();
        let atom_ring_counts = vec![1usize; atoms.len()];
        embed_fused_systems(
            atoms,
            bonds,
            &adjacency,
            &degree,
            &cip_ranks,
            &atom_ring_counts,
            &fused_rings,
            &mut efrags,
            Some(&coord_map),
            true,
        );
        assert_eq!(efrags.len(), 1);
        let embedded: BTreeSet<_> = efrags[0].eatoms.keys().copied().collect();
        let expected: BTreeSet<_> = fused_rings.iter().flatten().copied().collect();
        assert_eq!(embedded, expected);
    }

    #[test]
    fn embed_cis_trans_systems_adds_only_non_ring_stereo_double_bonds() {
        let mol = Molecule::from_smiles("F/C=C/F").unwrap();
        let atoms = mol.atoms();
        let bonds = mol.bonds();
        let n = atoms.len();
        let mut adjacency = vec![Vec::<usize>::new(); n];
        let mut degree = vec![0usize; n];
        for bond in bonds {
            adjacency[bond.begin().index()].push(bond.end().index());
            adjacency[bond.end().index()].push(bond.begin().index());
            degree[bond.begin().index()] += 1;
            degree[bond.end().index()] += 1;
        }
        let cip_ranks = vec![0u32; n];
        let mut efrags = Vec::new();
        embed_cis_trans_systems(
            atoms,
            bonds,
            &adjacency,
            &degree,
            &cip_ranks,
            &BTreeSet::new(),
            &mut efrags,
        );
        assert_eq!(efrags.len(), 1);
        assert_eq!(efrags[0].eatoms.len(), 2);
        let coords: Vec<_> = efrags[0].eatoms.values().map(|st| st.loc).collect();
        assert!(coords.contains(&(0.0, 0.0)));
        assert!(coords.contains(&(BOND_LEN, 0.0)));
    }

    #[test]
    fn get_non_embedded_atoms_returns_only_uncovered_atom_indices() {
        let frag_a = RdkitEmbeddedFrag {
            eatoms: BTreeMap::from([
                (
                    0usize,
                    TreeEmbeddedAtom {
                        loc: (0.0, 0.0),
                        normal: (1.0, 0.0),
                        ccw: true,
                        cis_trans_nbr: None,
                        angle: -1.0,
                        nbr1: None,
                        nbr2: None,
                        rot_dir: 0,
                        pending: Vec::new(),
                        d_density: 0.0,
                        df_fixed: false,
                    },
                ),
                (
                    2usize,
                    TreeEmbeddedAtom {
                        loc: (1.0, 0.0),
                        normal: (1.0, 0.0),
                        ccw: true,
                        cis_trans_nbr: None,
                        angle: -1.0,
                        nbr1: None,
                        nbr2: None,
                        rot_dir: 0,
                        pending: Vec::new(),
                        d_density: 0.0,
                        df_fixed: false,
                    },
                ),
            ]),
            attach_pts: VecDeque::new(),
            done: false,
        };
        let frag_b = RdkitEmbeddedFrag {
            eatoms: BTreeMap::from([(
                4usize,
                TreeEmbeddedAtom {
                    loc: (2.0, 0.0),
                    normal: (1.0, 0.0),
                    ccw: true,
                    cis_trans_nbr: None,
                    angle: -1.0,
                    nbr1: None,
                    nbr2: None,
                    rot_dir: 0,
                    pending: Vec::new(),
                    d_density: 0.0,
                    df_fixed: false,
                },
            )]),
            attach_pts: VecDeque::new(),
            done: false,
        };
        assert_eq!(get_non_embedded_atoms(6, &[frag_a, frag_b]), vec![1, 3, 5]);
    }

    #[test]
    fn find_largest_frag_skips_done_fragments() {
        let mut small = RdkitEmbeddedFrag::default();
        small.eatoms.insert(
            0,
            TreeEmbeddedAtom {
                loc: (0.0, 0.0),
                normal: (1.0, 0.0),
                ccw: true,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: None,
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
                d_density: 0.0,
                df_fixed: false,
            },
        );
        let mut large_done = small.clone();
        large_done.eatoms.insert(
            1,
            TreeEmbeddedAtom {
                loc: (1.0, 0.0),
                ..small.eatoms[&0].clone()
            },
        );
        large_done.mark_done();
        let mut medium = small.clone();
        medium.eatoms.insert(
            2,
            TreeEmbeddedAtom {
                loc: (2.0, 0.0),
                ..small.eatoms[&0].clone()
            },
        );
        assert_eq!(find_largest_frag(&[small, large_done, medium]), Some(2));
    }

    #[test]
    fn shift_coords_translates_second_fragment_by_rdkit_box_rule() {
        let mk_atom = |loc| TreeEmbeddedAtom {
            loc,
            normal: (1.0, 0.0),
            ccw: true,
            cis_trans_nbr: None,
            angle: -1.0,
            nbr1: None,
            nbr2: None,
            rot_dir: 0,
            pending: Vec::new(),
            d_density: 0.0,
            df_fixed: false,
        };
        let mut efrags = vec![
            RdkitEmbeddedFrag {
                eatoms: BTreeMap::from([
                    (0usize, mk_atom((0.0, 0.0))),
                    (1usize, mk_atom((1.0, 0.0))),
                ]),
                attach_pts: VecDeque::new(),
                done: false,
            },
            RdkitEmbeddedFrag {
                eatoms: BTreeMap::from([
                    (2usize, mk_atom((0.0, 0.0))),
                    (3usize, mk_atom((1.0, 0.0))),
                ]),
                attach_pts: VecDeque::new(),
                done: false,
            },
        ];
        shift_coords(&mut efrags);
        let shifted = &efrags[1].eatoms;
        assert_eq!(shifted[&2].loc, (0.0, 1.0));
        assert_eq!(shifted[&3].loc, (1.0, 1.0));
    }

    #[test]
    fn copy_sign_and_thetabin_match_rdkit_helper_shape() {
        assert_eq!(copy_sign(2.5, -0.2, 0.1), -2.5);
        assert_eq!(copy_sign(2.5, -0.05, 0.1), 2.5);
        let bin = ThetaBin {
            d_theta_avg: 1.25,
            theta_values: vec![1.0, 1.5],
        };
        assert_eq!(bin.d_theta_avg, 1.25);
        assert_eq!(bin.theta_values, vec![1.0, 1.5]);
    }

    #[test]
    fn compute_initial_coords_prespec_coord_map_seeds_first_fragment() {
        let mol = Molecule::from_smiles("CC.C").unwrap();
        let mut coord_map = BTreeMap::new();
        coord_map.insert(0usize, [2.0, 3.0]);
        coord_map.insert(1usize, [3.5, 3.0]);

        let efrags = rdkit_compute_initial_efrags_strict(
            mol.atoms(),
            mol.bonds(),
            &vec![0; mol.num_atoms()],
            Some(&coord_map),
            false,
        )
        .unwrap();

        assert!(!efrags.is_empty());
        assert_eq!(efrags[0].get_embedded_atoms()[&0].loc, (2.0, 3.0));
        assert_eq!(efrags[0].get_embedded_atoms()[&1].loc, (3.5, 3.0));
    }

    #[test]
    fn compute_initial_coords_single_coord_map_entry_does_not_create_prespec_fragment() {
        let mol = Molecule::from_smiles("CC").unwrap();
        let mut coord_map = BTreeMap::new();
        coord_map.insert(0usize, [8.0, 9.0]);

        let efrags = rdkit_compute_initial_efrags_strict(
            mol.atoms(),
            mol.bonds(),
            &vec![0; mol.num_atoms()],
            Some(&coord_map),
            false,
        )
        .unwrap();

        assert!(!efrags.is_empty());
        assert_ne!(efrags[0].get_embedded_atoms()[&0].loc, (8.0, 9.0));
    }

    #[test]
    fn compute_initial_coords_embeds_fused_and_cis_trans_seed_fragments_before_expansion() {
        let fused = Molecule::from_smiles("C1CCC2CCCCC2C1").unwrap();
        let fused_efrags = rdkit_compute_initial_efrags_strict(
            fused.atoms(),
            fused.bonds(),
            &vec![0; fused.num_atoms()],
            None,
            false,
        )
        .unwrap();
        assert!(fused_efrags.iter().any(|frag| frag.size() >= 2));

        let cis = Molecule::from_smiles("F/C=C/F").unwrap();
        let cis_efrags = rdkit_compute_initial_efrags_strict(
            cis.atoms(),
            cis.bonds(),
            &vec![0; cis.num_atoms()],
            None,
            false,
        )
        .unwrap();
        assert!(cis_efrags.iter().any(|frag| frag.size() >= 2));
    }

    #[test]
    fn copy_coordinate_from_efrags_writes_embedded_positions_by_atom_index() {
        let mk_atom = |loc| TreeEmbeddedAtom {
            loc,
            normal: (1.0, 0.0),
            ccw: true,
            cis_trans_nbr: None,
            angle: -1.0,
            nbr1: None,
            nbr2: None,
            rot_dir: 0,
            pending: Vec::new(),
            d_density: 0.0,
            df_fixed: false,
        };
        let efrags = vec![
            RdkitEmbeddedFrag {
                eatoms: BTreeMap::from([
                    (1usize, mk_atom((1.0, 2.0))),
                    (3usize, mk_atom((3.0, 4.0))),
                ]),
                attach_pts: VecDeque::new(),
                done: false,
            },
            RdkitEmbeddedFrag {
                eatoms: BTreeMap::from([(0usize, mk_atom((-1.0, -2.0)))]),
                attach_pts: VecDeque::new(),
                done: false,
            },
        ];

        let coords = copy_coordinate_from_efrags(4, &efrags);
        assert_eq!(coords[0], [-1.0, -2.0]);
        assert_eq!(coords[1], [1.0, 2.0]);
        assert_eq!(coords[2], [0.0, 0.0]);
        assert_eq!(coords[3], [3.0, 4.0]);
    }

    #[test]
    fn template_line_parse_splits_smarts_and_cx_block() {
        let parts = parse_rdkit_template_line("C1CC1 |(0,0,;1,0,;0,1,)|").unwrap();
        assert_eq!(parts.smarts_body, "C1CC1");
        assert_eq!(parts.cx_block.as_deref(), Some("|(0,0,;1,0,;0,1,)|"));
        assert_eq!(parts.trailing_name, None);
    }

    #[test]
    fn template_line_parse_preserves_trailing_name_after_cx_block() {
        let parts = parse_rdkit_template_line("C1CC1 |(0,0,;1,0,;0,1,)| cyclopropane").unwrap();
        assert_eq!(parts.smarts_body, "C1CC1");
        assert_eq!(parts.cx_block.as_deref(), Some("|(0,0,;1,0,;0,1,)|"));
        assert_eq!(parts.trailing_name.as_deref(), Some("cyclopropane"));
    }

    #[test]
    fn template_line_parse_rejects_suffix_without_cx_prefix() {
        let err = parse_rdkit_template_line("C1CC1 trailing-text").unwrap_err();
        match err {
            Coordinate2DError::UnsupportedFeature(message) => {
                assert!(message.contains("suffix text without a CX block"));
            }
            other => panic!("expected UnsupportedFeature, got {other:?}"),
        }
    }

    #[test]
    fn template_line_parse_rejects_unclosed_cx_block() {
        let err = parse_rdkit_template_line("C1CC1 |(0,0,;1,0,;0,1,)").unwrap_err();
        match err {
            Coordinate2DError::UnsupportedFeature(message) => {
                assert!(message.contains("closing '|'"));
            }
            other => panic!("expected UnsupportedFeature, got {other:?}"),
        }
    }

    #[test]
    fn template_line_parse_rejects_empty_lines() {
        let err = parse_rdkit_template_line("   ").unwrap_err();
        match err {
            Coordinate2DError::UnsupportedFeature(message) => {
                assert!(message.contains("non-empty SMARTS line"));
            }
            other => panic!("expected UnsupportedFeature, got {other:?}"),
        }
    }

    #[test]
    fn template_graph_model_expands_ring_closure_into_explicit_bond() {
        let model = parse_rdkit_template_graph_model("C1CC1").unwrap();
        assert_eq!(model.atom_queries.len(), 3);
        assert_eq!(model.bonds.len(), 3);
        assert_eq!(
            model
                .bonds
                .iter()
                .map(|bond| (bond.begin_atom_idx, bond.end_atom_idx))
                .collect::<Vec<_>>(),
            vec![(0, 1), (1, 2), (0, 2)]
        );
    }

    #[test]
    fn template_graph_model_preserves_branch_connectivity() {
        let model = parse_rdkit_template_graph_model("C(O)N").unwrap();
        assert_eq!(model.atom_queries.len(), 3);
        assert_eq!(model.bonds.len(), 2);
        assert_eq!(
            model
                .bonds
                .iter()
                .map(|bond| (bond.begin_atom_idx, bond.end_atom_idx))
                .collect::<Vec<_>>(),
            vec![(0, 1), (0, 2)]
        );
    }

    #[test]
    fn template_graph_model_preserves_explicit_bond_query_kinds() {
        let model = parse_rdkit_template_graph_model("C=C").unwrap();
        assert_eq!(model.bonds.len(), 1);
        assert_eq!(
            model.bonds[0].query,
            crate::QueryNode::Predicate(crate::BondQueryPredicate::Order(BondOrder::Double))
        );
    }

    #[test]
    fn template_graph_model_rejects_unbalanced_ring_closure() {
        let err = parse_rdkit_template_graph_model("C1CC").unwrap_err();
        match err {
            Coordinate2DError::UnsupportedFeature(message) => {
                assert!(message.contains("unbalanced ring closure"));
            }
            other => panic!("expected UnsupportedFeature, got {other:?}"),
        }
    }

    #[test]
    fn template_runtime_model_preserves_2d_coords_and_ring_counts() {
        let runtime = build_rdkit_template_runtime_model("C1CC1 |(0,0,;1,0,;0,1,)|").unwrap();
        assert_eq!(
            runtime.source_coordinate_dim,
            Some(CoordinateDimension::TwoD)
        );
        assert_eq!(
            runtime.coords_2d,
            Some(vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        );
        assert!(runtime.conformer_3d.is_none());
        assert_eq!(runtime.fragment_count, 1);
        assert_eq!(runtime.bond_ring_counts, vec![1, 1, 1]);
    }

    #[test]
    fn template_runtime_model_marks_nonzero_z_coords_as_3d() {
        let runtime = build_rdkit_template_runtime_model("C1CC1 |(0,0,1,;1,0,0,;0,1,0,)|").unwrap();
        assert_eq!(
            runtime.source_coordinate_dim,
            Some(CoordinateDimension::ThreeD)
        );
        assert!(runtime.coords_2d.is_none());
        assert_eq!(
            runtime.conformer_3d,
            Some(vec![[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        );
    }

    #[test]
    fn template_runtime_model_tracks_multifragment_topology() {
        let runtime = build_rdkit_template_runtime_model("C.C |(0,0,;1,0,)|").unwrap();
        assert_eq!(runtime.fragment_count, 2);
        assert_eq!(runtime.bond_ring_counts, Vec::<usize>::new());
    }

    #[test]
    fn template_validation_accepts_single_fragment_ring_with_2d_coords() {
        let runtime = build_rdkit_template_runtime_model("C1CC1 |(0,0,;1,0,;0,1,)|").unwrap();
        assert!(assert_valid_rdkit_template(&runtime, "C1CC1 |(0,0,;1,0,;0,1,)|").is_ok());
    }

    #[test]
    fn template_validation_rejects_missing_coordinates() {
        let runtime = build_rdkit_template_runtime_model("C1CC1").unwrap();
        let err = assert_valid_rdkit_template(&runtime, "C1CC1").unwrap_err();
        match err {
            Coordinate2DError::UnsupportedFeature(message) => {
                assert_eq!(message, "Template missing coordinates: C1CC1");
            }
            other => panic!("expected UnsupportedFeature, got {other:?}"),
        }
    }

    #[test]
    fn template_validation_rejects_3d_coordinates() {
        let runtime = build_rdkit_template_runtime_model("C1CC1 |(0,0,1,;1,0,0,;0,1,0,)|").unwrap();
        let err = assert_valid_rdkit_template(&runtime, "C1CC1 3d").unwrap_err();
        match err {
            Coordinate2DError::UnsupportedFeature(message) => {
                assert_eq!(
                    message,
                    "Template has 3D coordinates, 2D coordinates required: C1CC1 3d"
                );
            }
            other => panic!("expected UnsupportedFeature, got {other:?}"),
        }
    }

    #[test]
    fn template_validation_rejects_multiple_fragments() {
        let runtime = build_rdkit_template_runtime_model("C.C |(0,0,;1,0,)|").unwrap();
        let err = assert_valid_rdkit_template(&runtime, "C.C").unwrap_err();
        match err {
            Coordinate2DError::UnsupportedFeature(message) => {
                assert_eq!(
                    message,
                    "Template consists of multiple fragments, single fragment required: C.C"
                );
            }
            other => panic!("expected UnsupportedFeature, got {other:?}"),
        }
    }

    #[test]
    fn template_validation_rejects_single_atom_template() {
        let runtime = build_rdkit_template_runtime_model("C |(0,0,)|").unwrap();
        let err = assert_valid_rdkit_template(&runtime, "C").unwrap_err();
        match err {
            Coordinate2DError::UnsupportedFeature(message) => {
                assert_eq!(message, "Template is not a ring system: C");
            }
            other => panic!("expected UnsupportedFeature, got {other:?}"),
        }
    }

    #[test]
    fn template_validation_rejects_non_ring_bonds() {
        let runtime = build_rdkit_template_runtime_model("CC |(0,0,;1,0,)|").unwrap();
        let err = assert_valid_rdkit_template(&runtime, "CC").unwrap_err();
        match err {
            Coordinate2DError::UnsupportedFeature(message) => {
                assert_eq!(message, "Template is not a ring system: CC");
            }
            other => panic!("expected UnsupportedFeature, got {other:?}"),
        }
    }

    #[test]
    fn template_loading_rejects_missing_file() {
        let mut templates = RdkitTemplateBuckets::new();
        let err = load_rdkit_templates_from_path(
            "/definitely/not/present/rdkit_templates.cxsmiles",
            &mut templates,
        )
        .unwrap_err();
        match err {
            Coordinate2DError::UnsupportedFeature(message) => {
                assert_eq!(
                    message,
                    "Could not open file /definitely/not/present/rdkit_templates.cxsmiles"
                );
            }
            other => panic!("expected UnsupportedFeature, got {other:?}"),
        }
    }

    #[test]
    fn template_loading_rejects_invalid_smarts() {
        let file = NamedTempFile::new().unwrap();
        std::fs::write(file.path(), "C1CC\n").unwrap();
        let mut templates = RdkitTemplateBuckets::new();
        let err = load_rdkit_templates_from_path(file.path().to_str().unwrap(), &mut templates)
            .unwrap_err();
        match err {
            Coordinate2DError::UnsupportedFeature(message) => {
                assert_eq!(
                    message,
                    format!(
                        "Could not load templates from {}: Invalid smarts",
                        file.path().display()
                    )
                );
            }
            other => panic!("expected UnsupportedFeature, got {other:?}"),
        }
    }

    #[test]
    fn template_loading_forwards_template_validation_failure() {
        let file = NamedTempFile::new().unwrap();
        std::fs::write(file.path(), "CC |(0,0,;1,0,)|\n").unwrap();
        let mut templates = RdkitTemplateBuckets::new();
        let err = load_rdkit_templates_from_path(file.path().to_str().unwrap(), &mut templates)
            .unwrap_err();
        match err {
            Coordinate2DError::UnsupportedFeature(message) => {
                assert_eq!(message, "Template is not a ring system: CC |(0,0,;1,0,)|");
            }
            other => panic!("expected UnsupportedFeature, got {other:?}"),
        }
    }

    #[test]
    fn template_loading_preserves_duplicate_insertion_order_within_size_bucket() {
        let file = NamedTempFile::new().unwrap();
        std::fs::write(
            file.path(),
            concat!("C1CC1 |(0,0,;1,0,;0,1,)|\n", "C1CC1 |(5,5,;6,5,;5,6,)|\n"),
        )
        .unwrap();
        let mut templates = RdkitTemplateBuckets::new();
        load_rdkit_templates_from_path(file.path().to_str().unwrap(), &mut templates).unwrap();
        let bucket = templates.get(&3).unwrap();
        assert_eq!(bucket.len(), 2);
        assert_eq!(
            bucket[0].coords_2d,
            Some(vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        );
        assert_eq!(
            bucket[1].coords_2d,
            Some(vec![[5.0, 5.0], [6.0, 5.0], [5.0, 6.0]])
        );
    }

    struct TemplateRegistrySnapshotGuard {
        snapshot: RdkitTemplateBuckets,
    }

    impl TemplateRegistrySnapshotGuard {
        fn capture() -> Self {
            let snapshot = rdkit_coordinate_template_registry()
                .read()
                .expect("template registry lock poisoned")
                .templates
                .clone();
            Self { snapshot }
        }
    }

    impl Drop for TemplateRegistrySnapshotGuard {
        fn drop(&mut self) {
            rdkit_coordinate_template_registry()
                .write()
                .expect("template registry lock poisoned")
                .templates = self.snapshot.clone();
        }
    }

    fn template_registry_test_lock() -> &'static Mutex<()> {
        static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
        LOCK.get_or_init(|| Mutex::new(()))
    }

    struct PreferCoordGenGuard {
        snapshot: bool,
    }

    impl PreferCoordGenGuard {
        fn capture() -> Self {
            Self {
                snapshot: prefer_coord_gen(),
            }
        }
    }

    impl Drop for PreferCoordGenGuard {
        fn drop(&mut self) {
            set_prefer_coord_gen(self.snapshot);
        }
    }

    fn prefer_coordgen_test_lock() -> &'static Mutex<()> {
        static LOCK: OnceLock<Mutex<()>> = OnceLock::new();
        LOCK.get_or_init(|| Mutex::new(()))
    }

    #[test]
    fn ring_system_templates_default_registry_contains_templates() {
        let _lock = template_registry_test_lock().lock().unwrap();
        let _guard = TemplateRegistrySnapshotGuard::capture();
        rdkit_load_default_ring_system_templates();
        let expected = build_default_rdkit_coordinate_template_registry();
        let actual = rdkit_coordinate_template_registry()
            .read()
            .expect("template registry lock poisoned")
            .templates
            .clone();
        assert_eq!(actual, expected.templates);
    }

    #[test]
    fn ring_system_templates_set_replaces_registry_on_success() {
        let _lock = template_registry_test_lock().lock().unwrap();
        let _guard = TemplateRegistrySnapshotGuard::capture();
        let file = NamedTempFile::new().unwrap();
        std::fs::write(file.path(), "C1CC1 |(0,0,;1,0,;0,1,)|\n").unwrap();
        set_rdkit_ring_system_templates(file.path().to_str().unwrap()).unwrap();
        assert_eq!(rdkit_matching_templates(3).len(), 1);
        assert!(rdkit_matching_templates(4).is_empty());
    }

    #[test]
    fn ring_system_templates_set_failure_preserves_current_registry() {
        let _lock = template_registry_test_lock().lock().unwrap();
        let _guard = TemplateRegistrySnapshotGuard::capture();
        let seed = NamedTempFile::new().unwrap();
        std::fs::write(seed.path(), "C1CC1 |(0,0,;1,0,;0,1,)|\n").unwrap();
        set_rdkit_ring_system_templates(seed.path().to_str().unwrap()).unwrap();
        let before = rdkit_matching_templates(3);

        let invalid = NamedTempFile::new().unwrap();
        std::fs::write(invalid.path(), "C1CC\n").unwrap();
        let err = set_rdkit_ring_system_templates(invalid.path().to_str().unwrap()).unwrap_err();
        match err {
            Coordinate2DError::UnsupportedFeature(message) => {
                assert_eq!(
                    message,
                    format!(
                        "Could not load templates from {}: Invalid smarts",
                        invalid.path().display()
                    )
                );
            }
            other => panic!("expected UnsupportedFeature, got {other:?}"),
        }
        assert_eq!(rdkit_matching_templates(3), before);
    }

    #[test]
    fn ring_system_templates_add_prepends_new_templates_within_bucket() {
        let _lock = template_registry_test_lock().lock().unwrap();
        let _guard = TemplateRegistrySnapshotGuard::capture();
        let seed = NamedTempFile::new().unwrap();
        std::fs::write(seed.path(), "C1CC1 |(0,0,;1,0,;0,1,)|\n").unwrap();
        set_rdkit_ring_system_templates(seed.path().to_str().unwrap()).unwrap();

        let added = NamedTempFile::new().unwrap();
        std::fs::write(added.path(), "C1CC1 |(5,5,;6,5,;5,6,)|\n").unwrap();
        add_rdkit_ring_system_templates(added.path().to_str().unwrap()).unwrap();

        let bucket = rdkit_matching_templates(3);
        assert_eq!(bucket.len(), 2);
        assert_eq!(
            bucket[0].coords_2d,
            Some(vec![[5.0, 5.0], [6.0, 5.0], [5.0, 6.0]])
        );
        assert_eq!(
            bucket[1].coords_2d,
            Some(vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        );
    }

    #[test]
    fn ring_system_templates_load_default_resets_registry() {
        let _lock = template_registry_test_lock().lock().unwrap();
        let _guard = TemplateRegistrySnapshotGuard::capture();
        let file = NamedTempFile::new().unwrap();
        std::fs::write(file.path(), "C1CC1 |(0,0,;1,0,;0,1,)|\n").unwrap();
        set_rdkit_ring_system_templates(file.path().to_str().unwrap()).unwrap();
        assert_eq!(rdkit_matching_templates(3).len(), 1);

        rdkit_load_default_ring_system_templates();
        let expected = build_default_rdkit_coordinate_template_registry();
        let actual = rdkit_coordinate_template_registry()
            .read()
            .expect("template registry lock poisoned")
            .templates
            .clone();
        assert_eq!(actual, expected.templates);
    }
}

// ── Non-tetrahedral stereo embedding (RDDepictor.cpp:46-243) ────────────────
//
// Ported from RDKit RDDepictor.cpp (lines 46-243)
// Source reproduction protocol: dev/source_reproduction_protocol.md

/// RDKit bond length constant.
const BOND_LEN: f64 = 1.5;
/// RDKit 1/sqrt(2) constant for square planar layout.
// Preserve the source's rounded depiction constant instead of substituting a
// mathematically related standard-library value.
#[allow(clippy::approx_constant)]
const ISQRT2: f64 = 0.707_107;
/// RDKit sqrt(3)/2 constant for trigonal bipyramidal layout.
const SQRT3_2: f64 = 0.866_025;

// BEGIN RDKIT CPP FUNCTION: getRankedAtomNeighbors (RDDepictor.cpp:46-58)
// RDKit✔️✔️: std::vector<const RDKit::Atom *> getRankedAtomNeighbors(
// RDKit✔️✔️:     const RDKit::ROMol &mol, const RDKit::Atom *atom,
// RDKit✔️✔️:     const std::vector<int> &atomRanks)
/// Get neighbors of an atom sorted by atom rank (ascending).
fn get_ranked_atom_neighbors(atom_idx: usize, atom_ranks: &[i32], bonds: &[Bond]) -> Vec<usize> {
    let mut nbrs: Vec<usize> = bonds
        .iter()
        .filter_map(|bond| {
            if bond.begin().index() == atom_idx {
                Some(bond.end().index())
            } else if bond.end().index() == atom_idx {
                Some(bond.begin().index())
            } else {
                None
            }
        })
        .collect();
    nbrs.sort_by_key(|&nbr| atom_ranks[nbr]);
    nbrs
}

// BEGIN RDKIT CPP FUNCTION: embedSquarePlanar (RDDepictor.cpp:60-96)
// RDKit✔️✔️: void embedSquarePlanar(const RDKit::ROMol &mol, const RDKit::Atom *atom,
// RDKit✔️✔️:                        std::list<EmbeddedFrag> &efrags,
// RDKit✔️✔️:                        const std::vector<int> &atomRanks)
/// Embed a square planar center with its ligands. Produces a coordinate map.
fn embed_square_planar(
    atom_idx: usize,
    atoms: &[Atom],
    bonds: &[Bond],
    atom_ranks: &[i32],
) -> HashMap<usize, (f64, f64)> {
    let tag = atoms[atom_idx].chiral_tag();
    if tag != ChiralTag::SquarePlanar {
        return HashMap::new();
    }

    let nbrs = get_ranked_atom_neighbors(atom_idx, atom_ranks, bonds);
    if nbrs.is_empty() {
        return HashMap::new();
    }

    // RDKit ideal points for square planar: corners of a rotated square
    let ideal_points = [
        (ISQRT2 * BOND_LEN, ISQRT2 * BOND_LEN),
        (ISQRT2 * BOND_LEN, -ISQRT2 * BOND_LEN),
        (-ISQRT2 * BOND_LEN, -ISQRT2 * BOND_LEN),
        (-ISQRT2 * BOND_LEN, ISQRT2 * BOND_LEN),
    ];

    let mut coord_map = HashMap::new();
    coord_map.insert(atom_idx, (0.0, 0.0));
    coord_map.insert(nbrs[0], ideal_points[0]);

    let mut q2_full = false;
    for &nbr in &nbrs {
        if nbr == nbrs[0] {
            continue;
        }
        // Use ideal angle lookup
        let angle = crate::stereo::get_ideal_angle_between_ligands_by_slice(
            atom_idx, nbrs[0], nbr, atoms, bonds,
        );
        if (angle - 180.0).abs() < 0.1 {
            coord_map.insert(nbr, ideal_points[2]);
        } else {
            if !q2_full {
                coord_map.insert(nbr, ideal_points[1]);
                q2_full = true;
            } else {
                coord_map.insert(nbr, ideal_points[3]);
            }
        }
    }
    coord_map
}

// BEGIN RDKIT CPP FUNCTION: embedTBP (RDDepictor.cpp:98-134)
// RDKit✔️✔️: void embedTBP(const RDKit::ROMol &mol, const RDKit::Atom *atom,
// RDKit✔️✔️:                std::list<EmbeddedFrag> &efrags,
// RDKit✔️✔️:                const std::vector<int> &atomRanks)
/// Embed a trigonal bipyramidal center with its ligands.
fn embed_tbp(
    atom_idx: usize,
    atoms: &[Atom],
    bonds: &[Bond],
    _mol: &Molecule,
    atom_ranks: &[i32],
) -> HashMap<usize, (f64, f64)> {
    let tag = atoms[atom_idx].chiral_tag();
    if tag != ChiralTag::TrigonalBipyramidal {
        return HashMap::new();
    }

    let nbrs = get_ranked_atom_neighbors(atom_idx, atom_ranks, bonds);
    if nbrs.is_empty() {
        return HashMap::new();
    }

    // RDKit ideal points for TBP
    let ideal_points = [
        (0.0, BOND_LEN),                        // axial
        (0.0, -BOND_LEN),                       // axial
        (-SQRT3_2 * BOND_LEN, BOND_LEN / 2.0),  // equatorial
        (-SQRT3_2 * BOND_LEN, -BOND_LEN / 2.0), // equatorial
        (BOND_LEN, 0.0),                        // equatorial
    ];

    let mut coord_map = HashMap::new();
    coord_map.insert(atom_idx, (0.0, 0.0));

    let axial1 =
        crate::stereo::get_trigonal_bipyramidal_axial_atom_by_slice(atom_idx, 1, atoms, bonds);
    let axial2 =
        crate::stereo::get_trigonal_bipyramidal_axial_atom_by_slice(atom_idx, -1, atoms, bonds);

    if let Some(a1) = axial1 {
        coord_map.insert(a1, ideal_points[0]);
    }
    if let Some(a2) = axial2 {
        coord_map.insert(a2, ideal_points[1]);
    }

    let mut which_eq = 2usize;
    for &nbr in &nbrs {
        if Some(nbr) != axial1 && Some(nbr) != axial2 {
            coord_map.insert(nbr, ideal_points[which_eq]);
            which_eq += 1;
        }
    }
    coord_map
}

// BEGIN RDKIT CPP FUNCTION: embedOctahedral (RDDepictor.cpp:136-211)
// RDKit✔️✔️: void embedOctahedral(const RDKit::ROMol &mol, const RDKit::Atom *atom,
// RDKit✔️✔️:                     std::list<EmbeddedFrag> &efrags,
// RDKit✔️✔️:                     const std::vector<int> &atomRanks)
/// Embed an octahedral center with its ligands.
fn embed_octahedral(
    atom_idx: usize,
    atoms: &[Atom],
    bonds: &[Bond],
    _mol: &Molecule,
    atom_ranks: &[i32],
) -> HashMap<usize, (f64, f64)> {
    let tag = atoms[atom_idx].chiral_tag();
    if tag != ChiralTag::Octahedral {
        return HashMap::new();
    }

    let nbrs = get_ranked_atom_neighbors(atom_idx, atom_ranks, bonds);
    if nbrs.is_empty() {
        return HashMap::new();
    }

    // RDKit ideal points for octahedral
    let ideal_points = [
        (0.0, BOND_LEN),                        // axial
        (0.0, -BOND_LEN),                       // axial
        (SQRT3_2 * BOND_LEN, BOND_LEN / 2.0),   // equatorial
        (SQRT3_2 * BOND_LEN, -BOND_LEN / 2.0),  // equatorial
        (-SQRT3_2 * BOND_LEN, -BOND_LEN / 2.0), // equatorial
        (-SQRT3_2 * BOND_LEN, BOND_LEN / 2.0),  // equatorial
    ];

    let mut coord_map = HashMap::new();
    coord_map.insert(atom_idx, (0.0, 0.0));

    // Determine axial atoms by finding a pair with 180° angle, or fallback to first all-90
    let mut axial1: Option<usize> = None;
    let mut axial2: Option<usize> = None;
    for i in 0..nbrs.len() {
        let mut all90 = true;
        for j in (i + 1)..nbrs.len() {
            let angle = crate::stereo::get_ideal_angle_between_ligands_by_slice(
                atom_idx, nbrs[i], nbrs[j], atoms, bonds,
            );
            if (angle - 180.0).abs() < 0.1 {
                axial1 = Some(nbrs[i]);
                axial2 = Some(nbrs[j]);
                all90 = false;
                break;
            } else if (angle - 90.0).abs() > 0.1 {
                all90 = false;
            }
        }
        if all90 {
            axial1 = Some(nbrs[i]);
        }
        if axial1.is_some() {
            break;
        }
    }

    if let Some(a1) = axial1 {
        coord_map.insert(a1, ideal_points[0]);
    }
    if let Some(a2) = axial2 {
        coord_map.insert(a2, ideal_points[1]);
    }

    let mut ref_eq_atom1: Option<usize> = None;
    let mut ref_eq_atom2: Option<usize> = None;
    for &nbr in &nbrs {
        if Some(nbr) == axial1 || Some(nbr) == axial2 {
            continue;
        }
        if ref_eq_atom1.is_none() {
            ref_eq_atom1 = Some(nbr);
            coord_map.insert(nbr, ideal_points[2]);
            let across =
                crate::stereo::get_chiral_across_atom_by_atom_by_slice(atom_idx, nbr, atoms, bonds);
            if let Some(aa) = across {
                ref_eq_atom2 = Some(aa);
                coord_map.insert(aa, ideal_points[4]);
            }
        } else if Some(nbr) != ref_eq_atom1 && Some(nbr) != ref_eq_atom2 {
            coord_map.insert(nbr, ideal_points[3]);
            let across2 =
                crate::stereo::get_chiral_across_atom_by_atom_by_slice(atom_idx, nbr, atoms, bonds);
            if let Some(aa) = across2 {
                coord_map.insert(aa, ideal_points[5]);
            }
            break;
        }
    }
    coord_map
}

// BEGIN RDKIT CPP FUNCTION: embedNontetrahedralStereo (RDDepictor.cpp:213-243)
// RDKit✔️✔️: void embedNontetrahedralStereo(const RDKit::ROMol &mol,
// RDKit✔️✔️:                                std::list<EmbeddedFrag> &efrags,
// RDKit✔️✔️:                                const std::vector<int> &atomRanks)
/// Embed all non-tetrahedral stereo centers in a molecule.
/// Returns a list of coordinate maps, one per non-tetrahedral center.
fn embed_nontetrahedral_stereo(
    atoms: &[Atom],
    bonds: &[Bond],
    mol: &Molecule,
    atom_ranks: &[i32],
) -> Vec<HashMap<usize, (f64, f64)>> {
    let mut results = Vec::new();
    for atom in atoms {
        let tag = atom.chiral_tag();
        match tag {
            ChiralTag::SquarePlanar => {
                let cm = embed_square_planar(atom.id().index(), atoms, bonds, atom_ranks);
                if !cm.is_empty() {
                    results.push(cm);
                }
            }
            ChiralTag::TrigonalBipyramidal => {
                let cm = embed_tbp(atom.id().index(), atoms, bonds, mol, atom_ranks);
                if !cm.is_empty() {
                    results.push(cm);
                }
            }
            ChiralTag::Octahedral => {
                let cm = embed_octahedral(atom.id().index(), atoms, bonds, mol, atom_ranks);
                if !cm.is_empty() {
                    results.push(cm);
                }
            }
            _ => {}
        }
    }
    results
}
// END RDKIT CPP FUNCTION: embedNontetrahedralStereo
//
// Embed helper: find bond between two atoms by index.
fn bond_between_idx(bonds: &[Bond], a: usize, b: usize) -> Option<usize> {
    bonds.iter().find_map(|bond| {
        if (bond.begin().index() == a && bond.end().index() == b)
            || (bond.begin().index() == b && bond.end().index() == a)
        {
            Some(bond.id().index())
        } else {
            None
        }
    })
}

// BEGIN RDKIT CPP FUNCTION: embedFusedSystems (RDDepictor.cpp:246-296)
// RDKit✔️✔️: void embedFusedSystems(const RDKit::ROMol &mol,
// RDKit✔️✔️:                        const RDKit::VECT_INT_VECT &arings,
// RDKit✔️✔️:                        std::list<EmbeddedFrag> &efrags,
// RDKit✔️✔️:                        const RDGeom::INT_POINT2D_MAP *coordMap,
// RDKit✔️✔️:                        bool useRingTemplates) {
fn embed_fused_systems(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &[Vec<usize>],
    degree: &[usize],
    cip_ranks: &[u32],
    atom_ring_counts: &[usize],
    arings: &[Vec<usize>],
    efrags: &mut Vec<RdkitEmbeddedFrag>,
    coord_map: Option<&BTreeMap<usize, (f64, f64)>>,
    use_ring_templates: bool,
) {
    let mut ring_neighbors = vec![Vec::<usize>::new(); arings.len()];
    for i in 0..arings.len() {
        for j in i + 1..arings.len() {
            if arings[i].iter().any(|aid| arings[j].contains(aid)) {
                ring_neighbors[i].push(j);
                ring_neighbors[j].push(i);
            }
        }
    }

    fn pick_fused(curr: usize, neigh: &[Vec<usize>], done: &mut [bool], out: &mut Vec<usize>) {
        done[curr] = true;
        out.push(curr);
        for &nb in &neigh[curr] {
            if !done[nb] {
                pick_fused(nb, neigh, done, out);
            }
        }
    }

    let mut fus_done = vec![false; arings.len()];
    let mut curr = 0usize;
    while curr < arings.len() {
        let mut fused = Vec::new();
        pick_fused(curr, &ring_neighbors, &mut fus_done, &mut fused);
        let frings: Vec<Vec<usize>> = fused.iter().map(|&rid| arings[rid].clone()).collect();

        let mut allow_ring_templates = use_ring_templates;
        if use_ring_templates {
            if let Some(coord_map) = coord_map {
                let mut count = 0usize;
                for ring in &frings {
                    for &aid in ring {
                        if coord_map.contains_key(&aid) {
                            count += 1;
                        }
                    }
                }
                allow_ring_templates = count < 2;
            }
        }

        let mut efrag = RdkitEmbeddedFrag::default();
        let _ = efrag.embed_fused_rings(
            atoms,
            bonds,
            adjacency,
            degree,
            cip_ranks,
            atom_ring_counts,
            &frings,
            allow_ring_templates,
        );
        efrag.setup_new_neighs(atoms, bonds, adjacency, degree, cip_ranks);
        efrags.push(efrag);

        let mut next = None;
        for (idx, done) in fus_done.iter().enumerate() {
            if !*done {
                next = Some(idx);
                break;
            }
        }
        let Some(next_idx) = next else { break };
        curr = next_idx;
    }
}

// BEGIN RDKIT CPP FUNCTION: embedCisTransSystems (RDDepictor.cpp:298-319)
// RDKit✔️✔️: void embedCisTransSystems(const RDKit::ROMol &mol,
// RDKit✔️✔️:                           std::list<EmbeddedFrag> &efrags) {
fn embed_cis_trans_systems(
    atoms: &[Atom],
    bonds: &[Bond],
    adjacency: &[Vec<usize>],
    degree: &[usize],
    cip_ranks: &[u32],
    ring_bond_ids: &BTreeSet<usize>,
    efrags: &mut Vec<RdkitEmbeddedFrag>,
) {
    for bond in bonds {
        if bond.order() != BondOrder::Double {
            continue;
        }
        if !matches!(
            bond.stereo(),
            BondStereo::Cis | BondStereo::Trans | BondStereo::E | BondStereo::Z
        ) {
            continue;
        }
        if ring_bond_ids.contains(&bond.id().index()) {
            continue;
        }
        if bond.stereo_atoms().is_none_or(|atoms| atoms.len() != 2) {
            continue;
        }
        if let Some(mut efrag) = RdkitEmbeddedFrag::from_double_bond(bond) {
            efrag.setup_new_neighs(atoms, bonds, adjacency, degree, cip_ranks);
            efrags.push(efrag);
        }
    }
}

// BEGIN RDKIT CPP FUNCTION: getNonEmbeddedAtoms / _findLargestFrag / _shiftCoords /
// BEGIN RDKIT CPP FUNCTION: copySign / ThetaBin (RDDepictor.cpp:321-401)
// RDKit✔️✔️: helpers matching DepictorLocal utility surface.
fn get_non_embedded_atoms(num_atoms: usize, efrags: &[RdkitEmbeddedFrag]) -> Vec<usize> {
    let mut done = vec![false; num_atoms];
    for efrag in efrags {
        for &aid in efrag.get_embedded_atoms().keys() {
            done[aid] = true;
        }
    }
    (0..num_atoms).filter(|&aid| !done[aid]).collect()
}

fn find_largest_frag(efrags: &[RdkitEmbeddedFrag]) -> Option<usize> {
    let mut best = None;
    let mut best_size = 0usize;
    for (idx, efrag) in efrags.iter().enumerate() {
        if !efrag.is_done() && efrag.size() > best_size {
            best_size = efrag.size();
            best = Some(idx);
        }
    }
    best
}

fn shift_coords(efrags: &mut [RdkitEmbeddedFrag]) {
    if efrags.is_empty() {
        return;
    }
    let Some((mut xmax, xmin, mut ymax, ymin)) = efrags[0].compute_box() else {
        return;
    };
    if std::env::var("COSMOLKIT_DEBUG_DEPICT_ROW").ok().as_deref() == Some("58") {
        eprintln!(
            "COSMOL_SHIFT frag=0 box=({:.17},{:.17},{:.17},{:.17}) bits=({:#018x},{:#018x},{:#018x},{:#018x})",
            xmax,
            xmin,
            ymax,
            ymin,
            xmax.to_bits(),
            xmin.to_bits(),
            ymax.to_bits(),
            ymin.to_bits()
        );
    }
    for efrag in efrags.iter_mut().skip(1) {
        let Some((xp, xn, yp, yn)) = efrag.compute_box() else {
            continue;
        };
        let mut shift = (0.0, 0.0);
        let xshift = !(xmax + xmin > ymax + ymin);
        if xshift {
            shift.0 = xmax + xn + 1.0;
            xmax += xp + xn + 1.0;
        } else {
            shift.1 = ymax + yn + 1.0;
            ymax += yp + yn + 1.0;
        }
        if std::env::var("COSMOLKIT_DEBUG_DEPICT_ROW").ok().as_deref() == Some("58") {
            eprintln!(
                "COSMOL_SHIFT frag=next xshift={} box=({:.17},{:.17},{:.17},{:.17}) shift=({:.17},{:.17}) bits=({:#018x},{:#018x},{:#018x},{:#018x}) ({:#018x},{:#018x})",
                xshift,
                xp,
                xn,
                yp,
                yn,
                shift.0,
                shift.1,
                xp.to_bits(),
                xn.to_bits(),
                yp.to_bits(),
                yn.to_bits(),
                shift.0.to_bits(),
                shift.1.to_bits()
            );
        }
        efrag.translate(shift);
    }
}

fn copy_sign(to: f64, from: f64, tol: f64) -> f64 {
    if from < -tol { -to.abs() } else { to.abs() }
}

#[derive(Clone, Debug, Default)]
struct ThetaBin {
    d_theta_avg: f64,
    theta_values: Vec<f64>,
}

// RDKit✔️❌: RDGeom::Transform3D *computeCanonicalTransform(const Conformer &conf,
// RDKit✔️❌:                                                const RDGeom::Point3D *center,
// RDKit✔️❌:                                                bool normalizeCovar,
// RDKit✔️❌:                                                bool ignoreHs,
// RDKit✔️❌:                                                double *eigenValues) { ... }
// RDKit✔️❌: 2D planar adaptation of MolTransforms canonical transform for
// RDKit✔️❌: depiction conformers stored as [x,y] rows.
fn rdkit_compute_canonical_transform_for_2d_coords(
    coords: &[[f64; 2]],
    center: Option<[f64; 2]>,
) -> [[f64; 4]; 4] {
    let mut trans = rdkit_transform3d_identity();
    let origin = if let Some(center) = center {
        center
    } else if coords.is_empty() {
        [0.0, 0.0]
    } else {
        let mut sum = [0.0, 0.0];
        for point in coords {
            sum[0] += point[0];
            sum[1] += point[1];
        }
        [sum[0] / coords.len() as f64, sum[1] / coords.len() as f64]
    };
    if coords.len() > 1 {
        let mut xx = 0.0;
        let mut xy = 0.0;
        let mut yy = 0.0;
        for point in coords {
            let dx = point[0] - origin[0];
            let dy = point[1] - origin[1];
            xx += dx * dx;
            xy += dx * dy;
            yy += dy * dy;
        }
        let d = ((xx - yy) * (xx - yy) + 4.0 * xy * xy).sqrt();
        let mut eig1 = (2.0 * xy, (yy - xx) + d);
        let eig1_len = norm(eig1);
        if eig1_len > 1.0e-4 {
            eig1 = (eig1.0 / eig1_len, eig1.1 / eig1_len);
        } else {
            eig1 = (1.0, 0.0);
        }
        let eig2 = (-eig1.1, eig1.0);
        trans[0][0] = eig1.0;
        trans[0][1] = eig1.1;
        trans[1][0] = eig2.0;
        trans[1][1] = eig2.1;
        trans[2][2] = 1.0;
    }
    let neg_origin = [-origin[0], -origin[1], 0.0];
    let transformed_origin = rdkit_transform3d_transform_point(&trans, neg_origin);
    rdkit_transform3d_set_translation(&mut trans, transformed_origin);
    trans
}

// RDKit✔️❌: void straightenDepiction(RDKit::ROMol &mol, int confId, bool minimizeRotation) {
// RDKit✔️❌:   if (!mol.getNumBonds()) {
// RDKit✔️❌:     return;
// RDKit✔️❌:   }
// RDKit✔️❌:   constexpr double RAD2DEG = 180. / M_PI;
// RDKit✔️❌:   constexpr double DEG2RAD = M_PI / 180.;
// RDKit✔️❌:   constexpr double ALMOST_ZERO = 1.e-5;
// RDKit✔️❌:   constexpr double INCR_DEG = 30.;
// RDKit✔️❌:   constexpr double HALF_INCR_DEG = 0.5 * INCR_DEG;
// RDKit✔️❌:   constexpr double QUARTER_INCR_DEG = 0.25 * INCR_DEG;
// RDKit✔️❌:   ... theta binning and rotation selection ...
pub(crate) fn straighten_depiction(
    mol: &mut Molecule,
    conf_id: isize,
    minimize_rotation: bool,
) -> Result<(), Coordinate2DError> {
    if mol.num_bonds() == 0 {
        return Ok(());
    }
    const RAD2DEG: f64 = 180.0 / PI;
    const DEG2RAD: f64 = PI / 180.0;
    const ALMOST_ZERO: f64 = 1.0e-5;
    const INCR_DEG: f64 = 30.0;
    const HALF_INCR_DEG: f64 = 0.5 * INCR_DEG;
    const QUARTER_INCR_DEG: f64 = 0.25 * INCR_DEG;

    let conf_index = rdkit_select_2d_conformer_index(mol, conf_id)?;
    let coords_snapshot = mol.conformers_2d()[conf_index].coordinates().to_vec();
    let mut theta_bins = std::collections::HashMap::<i32, ThetaBin>::new();
    for bond in mol.bonds() {
        let bi = bond.begin().index();
        let ei = bond.end().index();
        let mut bv = [
            coords_snapshot[bi][0] - coords_snapshot[ei][0],
            coords_snapshot[bi][1] - coords_snapshot[ei][1],
        ];
        bv[0] = if bv[0] < 0.0 {
            bv[0].min(-ALMOST_ZERO)
        } else {
            bv[0].max(ALMOST_ZERO)
        };
        let theta = RAD2DEG * (bv[1] / bv[0]).atan();
        let mut d_theta = (-theta).rem_euclid(INCR_DEG);
        if d_theta > HALF_INCR_DEG {
            d_theta -= INCR_DEG;
        }
        let theta_key = (d_theta + copy_sign(0.5, d_theta, ALMOST_ZERO)) as i32;
        let theta_bin = theta_bins.entry(theta_key).or_default();
        theta_bin.d_theta_avg += d_theta;
        theta_bin.theta_values.push(theta);
    }
    if theta_bins.is_empty() {
        return Ok(());
    }
    let mut d_theta_smallest = f64::MAX;
    for theta_bin in theta_bins.values_mut() {
        theta_bin.d_theta_avg /= theta_bin.theta_values.len() as f64;
        if theta_bin.d_theta_avg.abs() < d_theta_smallest.abs() {
            d_theta_smallest = theta_bin.d_theta_avg;
        }
    }
    let min_rotation_bin = theta_bins
        .values()
        .max_by(|a, b| {
            a.theta_values
                .len()
                .cmp(&b.theta_values.len())
                .then_with(|| b.d_theta_avg.abs().total_cmp(&a.d_theta_avg.abs()))
        })
        .expect("theta_bins is non-empty");
    let mut d_theta_min = min_rotation_bin.d_theta_avg;
    if !minimize_rotation {
        let mut count_60_vs_30 = [0u32, 0u32];
        for theta in &min_rotation_bin.theta_values {
            let abs_theta = (*theta + d_theta_min).abs();
            if abs_theta < ALMOST_ZERO {
                continue;
            }
            let idx = (((abs_theta + 0.5) / INCR_DEG) as usize) % 2;
            count_60_vs_30[idx] += 1;
        }
        if count_60_vs_30[0] > count_60_vs_30[1] {
            d_theta_min -= copy_sign(INCR_DEG, d_theta_min, ALMOST_ZERO);
        }
    } else if d_theta_smallest.abs() < ALMOST_ZERO
        || (d_theta_smallest.abs() < d_theta_min.abs() && d_theta_min.abs() > QUARTER_INCR_DEG)
    {
        d_theta_min = d_theta_smallest;
    }
    if d_theta_min.abs() > ALMOST_ZERO {
        let trans = transform2d_set_transform_center_angle((0.0, 0.0), d_theta_min * DEG2RAD);
        let coords = mol.coordinate_block_mut().conformers_2d[conf_index].coordinates_mut();
        for point in coords.iter_mut() {
            let rotated = transform2d_point((point[0], point[1]), trans);
            point[0] = rotated.0;
            point[1] = rotated.1;
        }
    }
    Ok(())
}

// RDKit✔️❌: double normalizeDepiction(RDKit::ROMol &mol, int confId, int canonicalize,
// RDKit✔️❌:                           double scaleFactor) {
// RDKit✔️❌:   constexpr double SCALE_FACTOR_THRESHOLD = 1.e-5;
// RDKit✔️❌:   if (!mol.getNumBonds()) {
// RDKit✔️❌:     return -1.;
// RDKit✔️❌:   }
// RDKit✔️❌:   ... most-common-bond-length scaling and canonical transform ...
pub(crate) fn normalize_depiction(
    mol: &mut Molecule,
    conf_id: isize,
    canonicalize: i32,
    mut scale_factor: f64,
) -> Result<f64, Coordinate2DError> {
    const SCALE_FACTOR_THRESHOLD: f64 = 1.0e-5;
    const RDKIT_BOND_LEN: f64 = 1.5;
    if mol.num_bonds() == 0 {
        return Ok(-1.0);
    }
    let conf_index = rdkit_select_2d_conformer_index(mol, conf_id)?;
    if scale_factor < 0.0 {
        let coords = mol.conformers_2d()[conf_index].coordinates();
        let mut most_common_bond_length_int = -1i32;
        let mut max_count = 0u32;
        let mut binned_bond_lengths = std::collections::HashMap::<i32, u32>::new();
        for bond in mol.bonds() {
            let begin = bond.begin().index();
            let end = bond.end().index();
            let dx = coords[begin][0] - coords[end][0];
            let dy = coords[begin][1] - coords[end][1];
            let bond_length = ((rdkit_sqrt(dx * dx + dy * dy) * 10.0) + 0.5) as i32;
            let count = binned_bond_lengths.entry(bond_length).or_insert(0);
            *count += 1;
            if *count > max_count {
                max_count = *count;
                most_common_bond_length_int = bond_length;
            }
        }
        if most_common_bond_length_int > 0 {
            let most_common_bond_length = most_common_bond_length_int as f64 * 0.1;
            scale_factor = RDKIT_BOND_LEN / most_common_bond_length;
        }
    }

    let coords_snapshot = mol.conformers_2d()[conf_index].coordinates().to_vec();
    let mut canon_trans = if canonicalize != 0 {
        Some(rdkit_compute_canonical_transform_for_2d_coords(
            &coords_snapshot,
            None,
        ))
    } else {
        let mut trans = rdkit_transform3d_identity();
        let mut centroid = [0.0, 0.0];
        for point in &coords_snapshot {
            centroid[0] += point[0];
            centroid[1] += point[1];
        }
        centroid[0] /= coords_snapshot.len() as f64;
        centroid[1] /= coords_snapshot.len() as f64;
        rdkit_transform3d_set_translation(&mut trans, [-centroid[0], -centroid[1], 0.0]);
        Some(trans)
    };
    if canonicalize < 0 {
        let rotate90 = transform2d_set_transform_center_angle((0.0, 0.0), PI / 2.0);
        let mut rotate90_3d = rdkit_transform3d_identity();
        rotate90_3d[0][0] = rotate90[0];
        rotate90_3d[0][1] = rotate90[1];
        rotate90_3d[0][3] = rotate90[2];
        rotate90_3d[1][0] = rotate90[3];
        rotate90_3d[1][1] = rotate90[4];
        rotate90_3d[1][3] = rotate90[5];
        if let Some(trans) = canon_trans.as_mut() {
            *trans = rdkit_transform3d_mul(&rotate90_3d, trans);
        }
    }

    let is_scale_factor_sane = scale_factor > SCALE_FACTOR_THRESHOLD;
    let coords = mol.coordinate_block_mut().conformers_2d[conf_index].coordinates_mut();
    if is_scale_factor_sane && (scale_factor - 1.0).abs() > SCALE_FACTOR_THRESHOLD {
        let mut trans = rdkit_transform3d_identity();
        trans[0][0] = scale_factor;
        trans[1][1] = scale_factor;
        if let Some(canon_trans) = canon_trans {
            trans = rdkit_transform3d_mul(&trans, &canon_trans);
        }
        for point in coords.iter_mut() {
            let transformed = rdkit_transform3d_transform_point(&trans, [point[0], point[1], 0.0]);
            point[0] = transformed[0];
            point[1] = transformed[1];
        }
    } else if let Some(canon_trans) = canon_trans {
        for point in coords.iter_mut() {
            let transformed =
                rdkit_transform3d_transform_point(&canon_trans, [point[0], point[1], 0.0]);
            point[0] = transformed[0];
            point[1] = transformed[1];
        }
    }
    if !is_scale_factor_sane {
        scale_factor = -1.0;
    }
    Ok(scale_factor)
}
// ── End non-tetrahedral embed functions ────────────────────────────────────
