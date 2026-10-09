//! Complete source-backed CrystalFF preferences over detached topology and QueryGraph.
use cosmolkit_core::{RingInfo, ValenceAssignment};
use cosmolkit_model::{
    BondId, BondOrder, BondQueryPredicate, CoordinateBlock, Hybridization, QueryGraph, QueryNode,
    TopologyBlock,
};
use cosmolkit_search::{
    SearchTarget, SmartsParseError, SmartsParseParams, SubstructMatchParams,
    build_prepared_query_match_context, parse_smarts,
    try_get_substruct_matches_with_params_and_context,
};
use std::{
    collections::BTreeMap,
    env,
    sync::{Arc, Mutex, OnceLock},
    time::Instant,
};
const MIN_MACROCYCLE_SIZE: usize = 9;

const TORSION_PREFERENCES_V1: &str = include_str!("params/v1.txt");
const TORSION_PREFERENCES_V2: &str = include_str!("params/v2.txt");
const TORSION_PREFERENCES_SMALL_RINGS: &str = include_str!("params/smallrings.txt");
const TORSION_PREFERENCES_MACROCYCLES: &str = include_str!("params/macrocycles.txt");

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum CrystalffTorsionPreferencesError {
    #[error("RDKit CrystalFF ETversion must be 1 or 2, got {version}")]
    InvalidExperimentalTorsionVersion { version: u32 },
    #[error("RDKit CrystalFF torsion-preference line is missing SMARTS pattern: {line}")]
    MissingSmartsPattern { line: String },
    #[error(
        "RDKit CrystalFF torsion-preference line ended before all sign/force-constant pairs were parsed: {line}"
    )]
    IncompleteParameterLine { line: String },
    #[error(
        "RDKit CrystalFF torsion-preference token '{token}' is not a valid integer in line: {line}"
    )]
    InvalidIntegerToken { token: String, line: String },
    #[error(
        "RDKit CrystalFF torsion-preference token '{token}' is not a valid floating-point value in line: {line}"
    )]
    InvalidFloatToken { token: String, line: String },
    #[error(
        "RDKit CrystalFF parameter source {constant_name} has unexpected C++ string-constant format"
    )]
    InvalidParameterSourceFormat { constant_name: &'static str },
    #[error(
        "RDKit CrystalFF parameter source {constant_name} contains unsupported escape sequence \\{escape}"
    )]
    UnsupportedEscapeSequence {
        constant_name: &'static str,
        escape: char,
    },
    #[error(
        "RDKit CrystalFF parameter source {constant_name} ended inside a quoted string literal"
    )]
    UnterminatedQuotedLiteral { constant_name: &'static str },
    #[error(
        "RDKit CrystalFF parameter source {constant_name} ended with a trailing backslash escape"
    )]
    TrailingEscape { constant_name: &'static str },
    #[error("RDKit CrystalFF molecule has no atoms")]
    EmptyMolecule,
    #[error(transparent)]
    PreparedTarget(#[from] cosmolkit_search::QueryMatchContextError),
    #[error(transparent)]
    Match(#[from] cosmolkit_search::SubstructMatchError),
    #[error("RDKit CrystalFF torsion SMARTS query build failed for '{smarts}': {reason}")]
    QueryBuildFailed { smarts: String, reason: String },
    #[error("RDKit CrystalFF requires initialized ring information: {reason}")]
    RingInfoUnavailable { reason: String },
    #[error(
        "RDKit CrystalFF match for SMARTS '{smarts}' did not contain the central bond between atoms {aid2} and {aid3}"
    )]
    MissingCentralBond {
        smarts: String,
        aid2: usize,
        aid3: usize,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub struct CrystalFFDetails {
    pub exp_torsion_atoms: Vec<Vec<i32>>,
    pub exp_torsion_angles: Vec<(Vec<i32>, Vec<f64>)>,
    pub improper_atoms: Vec<Vec<i32>>,
    pub bonds: Vec<(i32, i32)>,
    pub angles: Vec<Vec<i32>>,
    pub atom_nums: Vec<i32>,
    pub bounds_mat_force_scaling: f64,
    pub constrained_atoms: Vec<bool>,
}

impl Default for CrystalFFDetails {
    fn default() -> Self {
        Self {
            exp_torsion_atoms: Vec::new(),
            exp_torsion_angles: Vec::new(),
            improper_atoms: Vec::new(),
            bonds: Vec::new(),
            angles: Vec::new(),
            atom_nums: Vec::new(),
            bounds_mat_force_scaling: 1.0,
            constrained_atoms: Vec::new(),
        }
    }
}

pub type CrystalffTorsionBondMatch = (usize, Vec<usize>, ExpTorsionAngle);

#[derive(Debug, Clone)]
pub(crate) struct ExpTorsionAngle {
    torsion_idx: usize,
    smarts: String,
    force_constants: Vec<f64>,
    signs: Vec<i32>,
    query_molecule: Result<QueryGraph, SmartsParseError>,
    idx: [usize; 4],
}

impl ExpTorsionAngle {
    #[must_use]
    pub(crate) const fn torsion_idx(&self) -> usize {
        self.torsion_idx
    }

    #[must_use]
    pub(crate) fn smarts(&self) -> &str {
        &self.smarts
    }

    #[must_use]
    pub(crate) fn force_constants(&self) -> &[f64] {
        &self.force_constants
    }

    #[must_use]
    pub(crate) fn signs(&self) -> &[i32] {
        &self.signs
    }

    #[must_use]
    pub(crate) fn query_molecule(&self) -> Result<&QueryGraph, &SmartsParseError> {
        self.query_molecule.as_ref()
    }

    #[must_use]
    pub(crate) const fn idx(&self) -> [usize; 4] {
        self.idx
    }
}

#[derive(Debug, Clone)]
pub(crate) struct ExpTorsionAngleCollection {
    param_data: String,
    params: Vec<ExpTorsionAngle>,
}

impl ExpTorsionAngleCollection {
    #[must_use]
    pub(crate) fn param_data(&self) -> &str {
        &self.param_data
    }

    #[must_use]
    pub(crate) fn params(&self) -> &[ExpTorsionAngle] {
        &self.params
    }

    pub(crate) fn new(param_data: &str) -> Result<Self, CrystalffTorsionPreferencesError> {
        // BEGIN RDKIT CPP CONSTRUCTOR ForceFields::CrystalFF::ExpTorsionAngleCollection::ExpTorsionAngleCollection (TorsionPreferences.cpp:102-141)
        // RDKit✔️✔️: ExpTorsionAngleCollection::ExpTorsionAngleCollection(
        // RDKit✔️✔️:     const std::string &paramData) {
        // RDKit✔️✔️:   boost::char_separator<char> tabSep(" ", "", boost::drop_empty_tokens);
        // RDKit✔️✔️:   std::istringstream inStream(paramData);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   std::string inLine = RDKit::getLine(inStream);
        // RDKit✔️✔️:   unsigned int torsionIdx = 0;
        // RDKit✔️✔️:   while (!inStream.eof()) {
        // RDKit✔️✔️:     if (inLine[0] != '#') {
        // RDKit✔️✔️:       ExpTorsionAngle angle;
        // RDKit✔️✔️:       tokenizer tokens(inLine, tabSep);
        // RDKit✔️✔️:       tokenizer::iterator token = tokens.begin();
        // RDKit✔️✔️:       angle.smarts = *token;
        // RDKit✔️✔️:       angle.torsionIdx = torsionIdx++;
        // RDKit✔️✔️:       ++token;
        // RDKit✔️✔️:       for (unsigned int i = 0; i < 12; i += 2) {
        // RDKit✔️✔️:         angle.signs.push_back(boost::lexical_cast<int>(*token));
        // RDKit✔️✔️:         ++token;
        // RDKit✔️✔️:         angle.V.push_back(boost::lexical_cast<double>(*token));
        // RDKit✔️✔️:         ++token;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       angle.dp_pattern.reset(SmartsToMol(angle.smarts));
        // RDKit✔️✔️:       // get the atom indices for atom 1, 2, 3, 4 in the pattern
        // RDKit✔️✔️:       for (unsigned int i = 0; i < (angle.dp_pattern.get())->getNumAtoms();
        // RDKit✔️✔️:            ++i) {
        // RDKit✔️✔️:         Atom const *atom = (angle.dp_pattern.get())->getAtomWithIdx(i);
        // RDKit✔️✔️:         int num;
        // RDKit✔️✔️:         if (atom->getPropIfPresent("molAtomMapNumber", num)) {
        // RDKit✔️✔️:           if (num > 0 && num < 5) {
        // RDKit✔️✔️:             angle.idx[num - 1] = i;
        // RDKit✔️✔️:           }
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       d_params.push_back(std::move(angle));
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     inLine = RDKit::getLine(inStream);
        // RDKit✔️✔️:   }  // while loop
        // RDKit✔️✔️:   // std::cerr << "Exp. torsion angles = " << d_params.size() << " "
        // RDKit✔️✔️:   //    << d_params[d_params.size()-1].smarts << std::endl;
        // RDKit✔️✔️: }
        let mut params = Vec::new();

        for line in param_data.lines() {
            if line.starts_with('#') {
                continue;
            }
            let angle = parse_exp_torsion_angle_line(line, params.len())?;
            params.push(angle);
        }

        Ok(Self {
            param_data: param_data.to_owned(),
            params,
        })
    }

    pub(crate) fn get_params(
        version: u32,
        use_small_ring_torsions: bool,
        use_macrocycle_torsions: bool,
        param_data: &str,
    ) -> Result<Arc<Self>, CrystalffTorsionPreferencesError> {
        // BEGIN RDKIT CPP METHOD ForceFields::CrystalFF::ExpTorsionAngleCollection::getParams (TorsionPreferences.cpp:72-100)
        // RDKit✔️✔️: const ExpTorsionAngleCollection *ExpTorsionAngleCollection::getParams(
        // RDKit✔️✔️:     unsigned int version, bool useSmallRingTorsions, bool useMacrocycleTorsions,
        // RDKit✔️✔️:     const std::string &paramData) {
        // RDKit✔️✔️:   std::string params;
        // RDKit✔️✔️:   if (paramData.empty()) {
        // RDKit✔️✔️:     switch (version) {
        // RDKit✔️✔️:       case 1:
        // RDKit✔️✔️:         params = torsionPreferencesV1;
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       case 2:
        // RDKit✔️✔️:         params = torsionPreferencesV2;
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       default:
        // RDKit✔️✔️:         throw ValueErrorException("ETversion must be 1 or 2.");
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     params = paramData;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (useSmallRingTorsions) {
        // RDKit✔️✔️:     params += torsionPreferencesSmallRings;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (useMacrocycleTorsions) {
        // RDKit✔️✔️:     params += torsionPreferencesMacrocycles;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return &(param_flyweight(params).get());
        // RDKit✔️✔️: }
        // BEGIN RDKIT CPP TYPEDEF ForceFields::CrystalFF::param_flyweight (TorsionPreferences.cpp:59-63)
        // RDKit✔️✔️: typedef boost::flyweight<
        // RDKit✔️✔️:     boost::flyweights::key_value<std::string, ExpTorsionAngleCollection>,
        // RDKit✔️✔️:     boost::flyweights::no_tracking>
        // RDKit✔️✔️:     param_flyweight;
        // END RDKIT CPP TYPEDEF ForceFields::CrystalFF::param_flyweight
        let mut params = if param_data.is_empty() {
            match version {
                1 => torsion_preferences_v1()?,
                2 => torsion_preferences_v2()?,
                _ => {
                    return Err(
                        CrystalffTorsionPreferencesError::InvalidExperimentalTorsionVersion {
                            version,
                        },
                    );
                }
            }
            .to_owned()
        } else {
            param_data.to_owned()
        };

        if use_small_ring_torsions {
            params.push_str(torsion_preferences_small_rings()?);
        }

        if use_macrocycle_torsions {
            params.push_str(torsion_preferences_macrocycles()?);
        }

        crystalff_param_flyweight(&params)
    }
}

fn crystalff_param_flyweight(
    params: &str,
) -> Result<Arc<ExpTorsionAngleCollection>, CrystalffTorsionPreferencesError> {
    static PARAM_FLYWEIGHT: OnceLock<Mutex<BTreeMap<String, Arc<ExpTorsionAngleCollection>>>> =
        OnceLock::new();

    let cache = PARAM_FLYWEIGHT.get_or_init(|| Mutex::new(BTreeMap::new()));
    {
        let cache_guard = cache.lock().expect("torsion parameter cache poisoned");
        if let Some(collection) = cache_guard.get(params) {
            return Ok(Arc::clone(collection));
        }
    }

    let collection = Arc::new(ExpTorsionAngleCollection::new(params)?);
    let mut cache_guard = cache.lock().expect("torsion parameter cache poisoned");
    Ok(Arc::clone(
        cache_guard
            .entry(params.to_owned())
            .or_insert_with(|| Arc::clone(&collection)),
    ))
}

pub(crate) fn get_experimental_torsions(
    mol: &TopologyBlock,
    ring_info: &RingInfo,
    valence: &ValenceAssignment,
    details: &mut CrystalFFDetails,
    torsion_bonds: &mut Vec<CrystalffTorsionBondMatch>,
    use_exp_torsions: bool,
    use_small_ring_torsions: bool,
    use_macrocycle_torsions: bool,
    use_basic_knowledge: bool,
    version: u32,
    verbose: bool,
) -> Result<(), CrystalffTorsionPreferencesError> {
    // BEGIN RDKIT CPP FUNCTION ForceFields::CrystalFF::getExperimentalTorsions (TorsionPreferences.cpp:143-286)
    // RDKit✔️❗: void getExperimentalTorsions(
    // RDKit✔️❗:     const RDKit::ROMol &mol, CrystalFFDetails &details,
    // RDKit✔️❗:     std::vector<std::tuple<unsigned int, std::vector<unsigned int>,
    // RDKit✔️❗:                            const ExpTorsionAngle *>> &torsionBonds,
    // RDKit✔️❗:     bool useExpTorsions, bool useSmallRingTorsions, bool useMacrocycleTorsions,
    // RDKit✔️❗:     bool useBasicKnowledge, unsigned int version, bool verbose) {
    torsion_bonds.clear();
    let nb = mol.bonds.len();
    let na = mol.atoms.len();
    // RDKit✔️❗:   torsionBonds.clear();
    // RDKit✔️❗:   unsigned int nb = mol.getNumBonds();
    // RDKit✔️❗:   unsigned int na = mol.getNumAtoms();
    if na == 0 {
        // RDKit✔️❗:   if (!na) {
        // RDKit✔️❗:     throw ValueErrorException("molecule has no atoms");
        // RDKit✔️❗:   }
        return Err(CrystalffTorsionPreferencesError::EmptyMolecule);
    }

    // RDKit✔️❗:   // check that vectors are empty
    // RDKit✔️❗:   details.expTorsionAtoms.clear();
    // RDKit✔️❗:   details.expTorsionAngles.clear();
    // RDKit✔️❗:   details.improperAtoms.clear();
    details.exp_torsion_atoms.clear();
    details.exp_torsion_angles.clear();
    details.improper_atoms.clear();

    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(
        mol,
        &coordinates,
        &mol.stereo_groups,
        Some(ring_info),
        Some(valence),
    );
    // Reuse validated borrowed topology, rings and valence for every torsion SMARTS.
    // Behavior pending paired oracle; one O(V+E) validation, no cloned graph/cache.
    let query_context = build_prepared_query_match_context(mol, ring_info, valence)
        .map_err(CrystalffTorsionPreferencesError::PreparedTarget)?;
    let excluded_bonds = compute_excluded_bridged_bonds(nb, &ring_info);
    let mut done_bonds = vec![false; nb];

    if use_exp_torsions {
        let params = ExpTorsionAngleCollection::get_params(
            version,
            use_small_ring_torsions,
            use_macrocycle_torsions,
            "",
        )?;
        for param in params.params() {
            let query_molecule = param.query_molecule().map_err(|reason| {
                CrystalffTorsionPreferencesError::QueryBuildFailed {
                    smarts: param.smarts().to_owned(),
                    reason: reason.to_string(),
                }
            })?;
            let trace_row34_timing =
                env::var("COSMOLKIT_ROW34_TORSION_TIMING").ok().as_deref() == Some("1") && na == 44;
            let trace_row64_timing =
                env::var("RDKIT_ROW64_TRACE").ok().as_deref() == Some("1") && na == 106;
            if env::var("RDKIT_ROW61_TRACE").ok().as_deref() == Some("1")
                && na == 26
                && param.torsion_idx() == 349
            {
                println!(
                    "row61_torsion349_begin smarts={} query_atoms={} query_bonds={}",
                    param.smarts(),
                    query_molecule.num_atoms(),
                    query_molecule.num_bonds()
                );
            }
            let match_start = (trace_row34_timing || trace_row64_timing).then(Instant::now);
            let matches = try_get_substruct_matches_with_params_and_context(
                &target,
                query_molecule,
                &SubstructMatchParams {
                    max_matches: 1000,
                    uniquify: false,
                    use_chirality: false,
                    specified_stereo_query_matches_unspecified: false,
                    ..Default::default()
                },
                &query_context,
            )
            .map_err(CrystalffTorsionPreferencesError::Match)?;
            if env::var("RDKIT_ROW61_TRACE").ok().as_deref() == Some("1")
                && na == 26
                && param.torsion_idx() == 349
            {
                for (match_idx, match_result) in matches.iter().enumerate() {
                    println!(
                        "row61_torsion349_raw_match idx={match_idx} atom_mapping={:?}",
                        match_result.atom_mapping
                    );
                }
            }
            if let Some(start) = match_start {
                let elapsed = start.elapsed().as_secs_f64();
                if trace_row34_timing && elapsed >= 0.01 {
                    println!(
                        "row34_torsion_timing smarts={} elapsed={elapsed:.6} matches={}",
                        param.smarts(),
                        matches.len()
                    );
                } else if trace_row64_timing {
                    println!(
                        "row64_torsion_timing smarts={} elapsed={elapsed:.6} matches={}",
                        param.smarts(),
                        matches.len()
                    );
                }
            }
            for match_result in matches {
                let aid1 = match_result.atom_mapping[param.idx()[0]];
                let aid2 = match_result.atom_mapping[param.idx()[1]];
                let aid3 = match_result.atom_mapping[param.idx()[2]];
                let aid4 = match_result.atom_mapping[param.idx()[3]];
                let Some(bond) = bond_between_atoms(mol, aid2, aid3) else {
                    return Err(CrystalffTorsionPreferencesError::MissingCentralBond {
                        smarts: param.smarts().to_owned(),
                        aid2,
                        aid3,
                    });
                };
                let bid2 = bond.id().index();

                if excluded_bonds[bid2] || ring_info.num_bond_rings(BondId::new(bid2)) > 3 {
                    done_bonds[bid2] = true;
                }
                if !done_bonds[bid2] {
                    if !details.constrained_atoms.is_empty()
                        && details
                            .constrained_atoms
                            .get(aid1)
                            .copied()
                            .unwrap_or(false)
                        && details
                            .constrained_atoms
                            .get(aid2)
                            .copied()
                            .unwrap_or(false)
                        && details
                            .constrained_atoms
                            .get(aid3)
                            .copied()
                            .unwrap_or(false)
                        && details
                            .constrained_atoms
                            .get(aid4)
                            .copied()
                            .unwrap_or(false)
                    {
                        continue;
                    }
                    if env::var("COSMOLKIT_ROW16_TORSION_TRACE").ok().as_deref() == Some("1")
                        && na == 4
                    {
                        println!(
                            "row16_torsion_match bond_idx={bid2} atoms={:?} torsion_idx={} smarts={} signs={:?} force_constants={:?}",
                            vec![aid1, aid2, aid3, aid4],
                            param.torsion_idx(),
                            param.smarts(),
                            param.signs(),
                            param.force_constants()
                        );
                    }
                    if env::var("RDKIT_ROW34_TRACE").ok().as_deref() == Some("1") && na == 44 {
                        println!(
                            "row34_torsion_match bond_idx={bid2} atoms={:?} torsion_idx={} smarts={} signs={:?} force_constants={:?} excluded={} already_done={}",
                            vec![aid1, aid2, aid3, aid4],
                            param.torsion_idx(),
                            param.smarts(),
                            param.signs(),
                            param.force_constants(),
                            excluded_bonds[bid2],
                            done_bonds[bid2]
                        );
                    }
                    if env::var("COSMOLKIT_ROW20_TORSION_MATCH_TRACE")
                        .ok()
                        .as_deref()
                        == Some("1")
                        && na == 13
                    {
                        println!(
                            "row20_torsion_match bond_idx={bid2} atoms={:?} torsion_idx={} smarts={} signs={:?} force_constants={:?} excluded={} already_done={}",
                            vec![aid1, aid2, aid3, aid4],
                            param.torsion_idx(),
                            param.smarts(),
                            param.signs(),
                            param.force_constants(),
                            excluded_bonds[bid2],
                            done_bonds[bid2]
                        );
                    }
                    if env::var("RDKIT_ROW57_TRACE").ok().as_deref() == Some("1") && na == 27 {
                        println!(
                            "row57_torsion_match bond_idx={bid2} atoms={:?} torsion_idx={} smarts={} signs={:?} force_constants={:?} excluded={} already_done={}",
                            vec![aid1, aid2, aid3, aid4],
                            param.torsion_idx(),
                            param.smarts(),
                            param.signs(),
                            param.force_constants(),
                            excluded_bonds[bid2],
                            done_bonds[bid2]
                        );
                    }
                    if env::var("RDKIT_ROW61_TRACE").ok().as_deref() == Some("1") && na == 26 {
                        println!(
                            "row61_torsion_match bond_idx={bid2} atoms={:?} torsion_idx={} smarts={} signs={:?} force_constants={:?} excluded={} already_done={}",
                            vec![aid1, aid2, aid3, aid4],
                            param.torsion_idx(),
                            param.smarts(),
                            param.signs(),
                            param.force_constants(),
                            excluded_bonds[bid2],
                            done_bonds[bid2]
                        );
                    }
                    torsion_bonds.push((bid2, vec![aid1, aid2, aid3, aid4], param.clone()));
                    done_bonds[bid2] = true;
                    details.exp_torsion_atoms.push(vec![
                        i32::try_from(aid1).expect("atom index fits i32"),
                        i32::try_from(aid2).expect("atom index fits i32"),
                        i32::try_from(aid3).expect("atom index fits i32"),
                        i32::try_from(aid4).expect("atom index fits i32"),
                    ]);
                    details
                        .exp_torsion_angles
                        .push((param.signs().to_vec(), param.force_constants().to_vec()));
                    let _ = verbose;
                }
            }
        }
    }

    if use_basic_knowledge {
        let mut done_atoms = vec![false; na];
        for aid2 in 0..na {
            if done_atoms[aid2] {
                continue;
            }
            let mut atoms = vec![-1; 4];
            atoms[1] = i32::try_from(aid2).expect("atom index fits i32");
            let atom2 = &mol.atoms[aid2];
            let at2_atomic_num = atom2.atomic_number();
            if matches!(at2_atomic_num, 6 | 7 | 8)
                && atom2.hybridization() == Hybridization::Sp2
                && mol.adjacency.neighbors_of(aid2).len() == 3
            {
                let mut i = 0usize;
                let mut is_bound_to_sp2o = 0i32;
                for neighbor in mol.adjacency.neighbors_of(aid2) {
                    let atom_x = &mol.atoms[neighbor.atom_index];
                    atoms[i] = i32::try_from(atom_x.id().index()).expect("atom index fits i32");
                    if is_bound_to_sp2o == 0 {
                        is_bound_to_sp2o = i32::from(
                            at2_atomic_num == 6
                                && atom_x.atomic_number() == 8
                                && atom_x.hybridization() == Hybridization::Sp2,
                        );
                    }
                    if i == 0 {
                        i += 1;
                    }
                    i += 1;
                }
                atoms.push(i32::from(at2_atomic_num));
                atoms.push(is_bound_to_sp2o);
                details.improper_atoms.push(atoms);
            }
        }

        for atom_ring in ring_info.atom_rings() {
            let r_size = atom_ring.len();
            if !(4..=6).contains(&r_size) {
                continue;
            }
            for i in 0..r_size {
                let aid1 = atom_ring[i].index();
                let aid2 = atom_ring[(i + 1) % r_size].index();
                let aid3 = atom_ring[(i + 2) % r_size].index();
                let aid4 = atom_ring[(i + 3) % r_size].index();
                let bid2 = bond_between_atoms(mol, aid2, aid3)
                    .expect("ring atom sequence must have a central bond")
                    .id()
                    .index();
                if !done_bonds[bid2]
                    && mol.atoms[aid1].hybridization() == Hybridization::Sp2
                    && mol.atoms[aid2].hybridization() == Hybridization::Sp2
                    && mol.atoms[aid3].hybridization() == Hybridization::Sp2
                    && mol.atoms[aid4].hybridization() == Hybridization::Sp2
                {
                    done_bonds[bid2] = true;
                    details.exp_torsion_atoms.push(vec![
                        i32::try_from(aid1).expect("atom index fits i32"),
                        i32::try_from(aid2).expect("atom index fits i32"),
                        i32::try_from(aid3).expect("atom index fits i32"),
                        i32::try_from(aid4).expect("atom index fits i32"),
                    ]);
                    let mut signs = vec![1; 6];
                    signs[1] = -1;
                    let mut fconsts = vec![0.0; 6];
                    fconsts[1] = 100.0;
                    details.exp_torsion_angles.push((signs, fconsts));
                }
            }
        }
    }

    if env::var("RDKIT_ROW61_TRACE").ok().as_deref() == Some("1") && na == 26 {
        println!(
            "row61_details exp_torsion_atoms={:?} exp_torsion_angles={:?} improper_atoms={:?}",
            details.exp_torsion_atoms, details.exp_torsion_angles, details.improper_atoms
        );
    }

    Ok(())
}

pub fn get_experimental_torsions_without_bonds(
    mol: &TopologyBlock,
    ring_info: &RingInfo,
    valence: &ValenceAssignment,
    details: &mut CrystalFFDetails,
    use_exp_torsions: bool,
    use_small_ring_torsions: bool,
    use_macrocycle_torsions: bool,
    use_basic_knowledge: bool,
    version: u32,
    verbose: bool,
) -> Result<(), CrystalffTorsionPreferencesError> {
    // BEGIN RDKIT CPP FUNCTION ForceFields::CrystalFF::getExperimentalTorsions overload (TorsionPreferences.cpp:288-297)
    // RDKit✔️✔️: void getExperimentalTorsions(const RDKit::ROMol &mol, CrystalFFDetails &details,
    // RDKit✔️✔️:                              bool useExpTorsions, bool useSmallRingTorsions,
    // RDKit✔️✔️:                              bool useMacrocycleTorsions, bool useBasicKnowledge,
    // RDKit✔️✔️:                              unsigned int version, bool verbose) {
    // RDKit✔️✔️:   std::vector<std::tuple<unsigned int, std::vector<unsigned int>,
    // RDKit✔️✔️:                          const ExpTorsionAngle *>>
    // RDKit✔️✔️:       torsionBonds;
    // RDKit✔️✔️:   getExperimentalTorsions(mol, details, torsionBonds, useExpTorsions,
    // RDKit✔️✔️:                           useSmallRingTorsions, useMacrocycleTorsions,
    // RDKit✔️✔️:                           useBasicKnowledge, version, verbose);
    // RDKit✔️✔️: }
    let mut torsion_bonds = Vec::new();
    get_experimental_torsions(
        mol,
        ring_info,
        valence,
        details,
        &mut torsion_bonds,
        use_exp_torsions,
        use_small_ring_torsions,
        use_macrocycle_torsions,
        use_basic_knowledge,
        version,
        verbose,
    )
}

fn torsion_preferences_v1() -> Result<&'static str, CrystalffTorsionPreferencesError> {
    // BEGIN RECOVERY GEO-13 SOURCE torsion_preferences_v1
    // RDKit✔️✔️: const std::string torsionPreferencesV1 =
    // RDKit✔️✔️:     "[O:1]=[C:2]!@;-[O:3]~[CH0:4] -1 78.2 1 0.0 1 0.0 1 0.0 1 0.0 1 0.0\n"
    // RDKit✔️✔️:     "[O:1]=[C:2]([N])!@;-[O:3]~[C:4] -1 79.1 1 0.0 1 0.0 1 0.0 1 0.0 1 0.0\n"
    // RDKit✔️✔️:     "[O:1]=[C:2]!@;-[O:3]~[C:4] -1 100.0 1 0.0 1 0.0 1 0.0 1 0.0 1 0.0\n"
    // END RECOVERY GEO-13 SOURCE torsion_preferences_v1

    Ok(TORSION_PREFERENCES_V1)
}

fn torsion_preferences_v2() -> Result<&'static str, CrystalffTorsionPreferencesError> {
    Ok(TORSION_PREFERENCES_V2)
}

fn torsion_preferences_small_rings() -> Result<&'static str, CrystalffTorsionPreferencesError> {
    Ok(TORSION_PREFERENCES_SMALL_RINGS)
}

fn torsion_preferences_macrocycles() -> Result<&'static str, CrystalffTorsionPreferencesError> {
    Ok(TORSION_PREFERENCES_MACROCYCLES)
}

fn parse_exp_torsion_angle_line(
    line: &str,
    torsion_idx: usize,
) -> Result<ExpTorsionAngle, CrystalffTorsionPreferencesError> {
    let mut tokens = line.split_whitespace();
    let smarts = tokens
        .next()
        .ok_or_else(|| CrystalffTorsionPreferencesError::MissingSmartsPattern {
            line: line.to_owned(),
        })?
        .to_owned();
    let mut signs = Vec::with_capacity(6);
    let mut force_constants = Vec::with_capacity(6);

    for _ in 0..6 {
        let sign_token = tokens.next().ok_or_else(|| {
            CrystalffTorsionPreferencesError::IncompleteParameterLine {
                line: line.to_owned(),
            }
        })?;
        let sign = sign_token.parse::<i32>().map_err(|_| {
            CrystalffTorsionPreferencesError::InvalidIntegerToken {
                token: sign_token.to_owned(),
                line: line.to_owned(),
            }
        })?;
        signs.push(sign);

        let force_token = tokens.next().ok_or_else(|| {
            CrystalffTorsionPreferencesError::IncompleteParameterLine {
                line: line.to_owned(),
            }
        })?;
        let force_constant = force_token.parse::<f64>().map_err(|_| {
            CrystalffTorsionPreferencesError::InvalidFloatToken {
                token: force_token.to_owned(),
                line: line.to_owned(),
            }
        })?;
        force_constants.push(force_constant);
    }

    // RDKit✔️✔️:       angle.dp_pattern.reset(SmartsToMol(angle.smarts));
    // RDKit✔️✔️:       // get the atom indices for atom 1, 2, 3, 4 in the pattern
    // RDKit✔️✔️:       for (unsigned int i = 0; i < (angle.dp_pattern.get())->getNumAtoms();
    // RDKit✔️✔️:            ++i) {
    // RDKit✔️✔️:         Atom const *atom = (angle.dp_pattern.get())->getAtomWithIdx(i);
    // RDKit✔️✔️:         int num;
    // RDKit✔️✔️:         if (atom->getPropIfPresent("molAtomMapNumber", num)) {
    // RDKit✔️✔️:           if (num > 0 && num < 5) {
    // RDKit✔️✔️:             angle.idx[num - 1] = i;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // Local complexity review: CrystalFF performs one canonical SMARTS
    // compilation and one linear atom-map scan. The retained Result preserves
    // the collection's deferred parse-error behavior without reparsing.
    let query_molecule = parse_smarts(&smarts, &SmartsParseParams::default());
    let idx = query_molecule
        .as_ref()
        .map_or([0; 4], map_pattern_atom_indices);

    Ok(ExpTorsionAngle {
        torsion_idx,
        smarts,
        force_constants,
        signs,
        query_molecule,
        idx,
    })
}

fn bond_between_atoms(mol: &TopologyBlock, a: usize, b: usize) -> Option<&cosmolkit_model::Bond> {
    mol.bonds.iter().find(|bond| {
        let begin = bond.begin().index();
        let end = bond.end().index();
        (begin == a && end == b) || (begin == b && end == a)
    })
}

fn compute_excluded_bridged_bonds(nb: usize, ring_info: &RingInfo) -> Vec<bool> {
    let bond_rings = ring_info.bond_rings();
    let mut excluded_bonds = vec![false; nb];
    for (ri_idx, rii) in bond_rings.iter().enumerate() {
        let mut rs1 = vec![false; nb];
        for bond_id in rii {
            rs1[bond_id.index()] = true;
        }
        for rjj in bond_rings.iter().skip(ri_idx + 1) {
            if rii.len() >= MIN_MACROCYCLE_SIZE && rjj.len() >= MIN_MACROCYCLE_SIZE {
                continue;
            }
            let mut n_in_common = 0usize;
            for bond_id in rjj {
                if rs1[bond_id.index()] {
                    n_in_common += 1;
                    if n_in_common > 1 {
                        break;
                    }
                }
            }
            if n_in_common > 1 {
                if rii.len() < MIN_MACROCYCLE_SIZE {
                    for bond_id in rii {
                        excluded_bonds[bond_id.index()] = true;
                    }
                }
                if rjj.len() < MIN_MACROCYCLE_SIZE {
                    for bond_id in rjj {
                        excluded_bonds[bond_id.index()] = true;
                    }
                }
            }
        }
    }
    excluded_bonds
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct CrystalffTemplateBond {
    begin_atom_idx: usize,
    end_atom_idx: usize,
    query: QueryNode<BondQueryPredicate>,
}

fn expand_crystalff_smarts_bonds(
    query_molecule: &QueryGraph,
) -> Result<Vec<CrystalffTemplateBond>, String> {
    // Local complexity review: this is one O(E) pass over the canonical query
    // graph with one clone per bond query. It performs no SMARTS tokenization,
    // graph reconstruction, fallback, or consumer-local predicate synthesis.
    query_molecule
        .bonds()
        .iter()
        .enumerate()
        .map(|(bond_idx, bond)| {
            let query = bond.predicate().clone();
            Ok(CrystalffTemplateBond {
                begin_atom_idx: bond.begin().index(),
                end_atom_idx: bond.end().index(),
                query,
            })
        })
        .collect()
}

fn map_pattern_atom_indices(query_molecule: &QueryGraph) -> [usize; 4] {
    // RDKit✔️✔️:       for (unsigned int i = 0; i < (angle.dp_pattern.get())->getNumAtoms();
    // RDKit✔️✔️:            ++i) {
    // RDKit✔️✔️:         Atom const *atom = (angle.dp_pattern.get())->getAtomWithIdx(i);
    // RDKit✔️✔️:         int num;
    // RDKit✔️✔️:         if (atom->getPropIfPresent("molAtomMapNumber", num)) {
    // RDKit✔️✔️:           if (num > 0 && num < 5) {
    // RDKit✔️✔️:             angle.idx[num - 1] = i;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // Local complexity review: one allocation-free O(V) pass over the
    // canonical compiled query atoms matches RDKit's dp_pattern traversal.
    let mut idx = [0; 4];
    for (atom_idx, atom) in query_molecule.atoms().iter().enumerate() {
        if let Some(map_number) = atom.atom_map()
            && (1..=4).contains(&map_number)
        {
            idx[usize::try_from(map_number - 1).expect("atom-map index fits usize")] = atom_idx;
        }
    }
    idx
}

#[cfg(test)]
mod tests {
    #[test]
    fn recovery_geo13_v1_duplicate_removal_preserves_every_remaining_ordered_parameter() {
        let data = torsion_preferences_v1().unwrap();
        let first = data.lines().next().unwrap();
        assert_eq!(
            first,
            "[O:1]=[C:2]!@;-[O:3]~[CH0:4] -1 78.2 1 0.0 1 0.0 1 0.0 1 0.0 1 0.0"
        );
        assert_eq!(data.lines().filter(|line| *line == first).count(), 1);
        // Reinsert exactly the removed .1 entry, rather than reordering or
        // deduplicating arbitrary SMARTS. The complete vendored file is also
        // checked byte-for-byte against the official .6 file in GEO-13.json.
        let previous_data = format!("{first}\n{data}");
        let previous = ExpTorsionAngleCollection::new(&previous_data).unwrap();
        let current = ExpTorsionAngleCollection::get_params(1, false, false, "").unwrap();
        // C++ excludes 15 commented-out rules: the active table shrinks
        // from 376 to 375, despite 391/390 newline literals in raw text.
        assert_eq!(previous.params().len(), 376);
        assert_eq!(current.params().len(), 375);
        for (index, angle) in current.params().iter().enumerate() {
            let old_index = if index == 0 { 0 } else { index + 1 };
            let old = &previous.params()[old_index];
            assert_eq!(angle.torsion_idx(), index);
            assert_eq!(old.torsion_idx(), old_index);
            assert_eq!(angle.smarts(), old.smarts());
            assert_eq!(angle.idx(), old.idx());
            assert_eq!(angle.signs(), old.signs());
            assert_eq!(angle.force_constants(), old.force_constants());
            for cos_phi in [-1.0, -0.25, 0.0, 0.5, 1.0] {
                let energy = |parameter: &super::ExpTorsionAngle| {
                    crate::crystalff::torsion::calc_torsion_energy_m6(
                        parameter.force_constants(),
                        parameter.signs(),
                        cos_phi,
                    )
                };
                assert_eq!(energy(angle), energy(old));
            }
        }
    }

    #[test]
    fn recovery_geo13_v1_duplicate_removal_keeps_first_match_and_shifts_later_indices() {
        // Cases derived from the first three source SMARTS; both official
        // endpoints confirm atom tuples, first-match selection and parameters.
        for (smiles, expected_indices, expected_atoms, expected_barriers) in [
            (
                "CC(=O)OC(C)(C)C",
                vec![0, 41],
                vec![vec![2, 1, 3, 4], vec![5, 4, 3, 1]],
                vec![78.2, 0.0],
            ),
            ("NC(=O)OC", vec![1], vec![vec![2, 1, 3, 4]], vec![79.1]),
            ("CC(=O)OC", vec![2], vec![vec![2, 1, 3, 4]], vec![100.0]),
        ] {
            let mol = fixture_smiles(smiles).unwrap();
            let mut quiet_records = None;
            for verbose in [false, true] {
                let mut details = CrystalFFDetails::default();
                let mut bonds = Vec::new();
                get_experimental_torsions(
                    &mol,
                    &mut details,
                    &mut bonds,
                    true,
                    false,
                    false,
                    false,
                    1,
                    verbose,
                )
                .unwrap();
                let indices: Vec<_> = bonds.iter().map(|entry| entry.2.torsion_idx()).collect();
                let atoms: Vec<_> = bonds.iter().map(|entry| entry.1.clone()).collect();
                assert_eq!(indices, expected_indices, "{smiles}, verbose={verbose}");
                assert_eq!(atoms, expected_atoms);
                assert_eq!(
                    details.exp_torsion_atoms,
                    expected_atoms
                        .iter()
                        .map(|atoms| atoms.iter().map(|&atom| atom as i32).collect::<Vec<_>>())
                        .collect::<Vec<_>>()
                );
                assert_eq!(
                    bonds.iter().map(|entry| entry.0).collect::<Vec<_>>(),
                    if expected_indices.len() == 2 {
                        vec![2, 3]
                    } else {
                        vec![2]
                    }
                );
                for (index, (_, _, angle)) in bonds.iter().enumerate() {
                    assert_eq!(angle.force_constants()[0], expected_barriers[index]);
                    assert_eq!(
                        details.exp_torsion_angles[index],
                        (angle.signs().to_vec(), angle.force_constants().to_vec())
                    );
                    if index == 0 {
                        assert_eq!(angle.signs(), &[-1, 1, 1, 1, 1, 1]);
                        assert_eq!(&angle.force_constants()[1..], &[0.0; 5]);
                        for cos_phi in [-1.0, 0.0, 1.0] {
                            let energy = crate::crystalff::torsion::calc_torsion_energy_m6(
                                angle.force_constants(),
                                angle.signs(),
                                cos_phi,
                            );
                            assert_eq!(energy, expected_barriers[index] * (1.0 - cos_phi));
                        }
                    }
                }
                let records = (indices, atoms, details.exp_torsion_angles);
                if let Some(quiet) = &quiet_records {
                    assert_eq!(&records, quiet);
                } else {
                    quiet_records = Some(records);
                }
            }
        }
    }

    #[test]
    fn recovery_geo13_v1_duplicate_removal_preserves_optional_table_order_and_offsets() {
        let base = torsion_preferences_v1().unwrap();
        for (small, macrocycles) in [(false, false), (true, false), (false, true), (true, true)] {
            let mut expected = base.to_owned();
            if small {
                expected.push_str(torsion_preferences_small_rings().unwrap());
            }
            if macrocycles {
                expected.push_str(torsion_preferences_macrocycles().unwrap());
            }
            let params = ExpTorsionAngleCollection::get_params(1, small, macrocycles, "").unwrap();
            assert_eq!(params.param_data(), expected);
            assert_eq!(params.params().len(), expected.lines().count());
            for (index, angle) in params.params().iter().enumerate() {
                assert_eq!(angle.torsion_idx(), index);
                assert_eq!(
                    angle.smarts(),
                    expected
                        .lines()
                        .nth(index)
                        .unwrap()
                        .split_whitespace()
                        .next()
                        .unwrap()
                );
            }
            if small || macrocycles {
                let first_appended = if small {
                    torsion_preferences_small_rings().unwrap()
                } else {
                    torsion_preferences_macrocycles().unwrap()
                };
                assert_eq!(
                    params.params()[375].smarts(),
                    first_appended.split_whitespace().next().unwrap()
                );
            }
        }
    }

    use super::*;
    use cosmolkit_search::{atom_matches_query, bond_matches_query, compile_query_fixture};
    fn fixture_smiles(text: &str) -> Result<TopologyBlock, String> {
        let record =
            cosmolkit_smiles::parse_smiles(text, &Default::default()).map_err(|e| e.to_string())?;
        // Same default sanitization as the original Molecule::from_smiles fixture.
        cosmolkit_core::sanitize_topology(&record.topology, &Default::default())
            .map(|r| r.topology)
            .map_err(|e| e.to_string())
    }
    fn fixture_hydrogens(topology: TopologyBlock) -> TopologyBlock {
        cosmolkit_core::add_hydrogens_impl(topology, CoordinateBlock::default(), Default::default())
            .unwrap()
            .topology
    }
    fn fixture_molblock(text: &str) -> TopologyBlock {
        let (topology, _, _) = cosmolkit_io::sdf::read_v2000_detached(text).unwrap();
        cosmolkit_core::sanitize_topology(&topology, &Default::default())
            .unwrap()
            .topology
    }
    fn test_target(topology: &TopologyBlock) -> SearchTarget<'_> {
        static EMPTY: std::sync::OnceLock<CoordinateBlock> = std::sync::OnceLock::new();
        SearchTarget::new(
            topology,
            EMPTY.get_or_init(CoordinateBlock::default),
            &topology.stereo_groups,
            None,
            None,
        )
    }
    // Fixture forwarding performs chemistry through existing owners, not a test algorithm.
    fn get_substruct_matches_with_params(
        mol: &TopologyBlock,
        query: &QueryGraph,
        params: &SubstructMatchParams,
    ) -> Vec<cosmolkit_search::SubstructMatchResult> {
        // Original Molecule::from_smiles fixtures supplied their sanitized
        // cached state. Detached topology fixtures must explicitly retain that
        // state at SearchTarget; Search correctly does not recompute absent
        // caches when a reached getter requires them. These are fixture-only
        // calls to the canonical chemistry owners, not matcher fallbacks.
        let rings = cosmolkit_core::symmetrized_sssr(mol, &Default::default()).unwrap();
        let valence = cosmolkit_core::assign_valence_with_options_for_topology(
            mol,
            cosmolkit_core::ValenceModel::RdkitLike,
            false,
        )
        .unwrap();
        let empty_coordinates = CoordinateBlock::default();
        let target = SearchTarget::new(
            mol,
            &empty_coordinates,
            &mol.stereo_groups,
            Some(&rings),
            Some(&valence),
        );
        cosmolkit_search::try_get_substruct_matches_with_params(&target, query, params).unwrap()
    }
    fn get_experimental_torsions(
        mol: &TopologyBlock,
        details: &mut CrystalFFDetails,
        bonds: &mut Vec<CrystalffTorsionBondMatch>,
        exp: bool,
        small: bool,
        macrocycle: bool,
        basic: bool,
        version: u32,
        verbose: bool,
    ) -> Result<(), CrystalffTorsionPreferencesError> {
        // Preserve the source empty-molecule error before fixture preparation.
        if mol.atoms.is_empty() {
            return Err(CrystalffTorsionPreferencesError::EmptyMolecule);
        }
        let rings = cosmolkit_core::symmetrized_sssr(mol, &Default::default()).unwrap();
        let valence = cosmolkit_core::assign_valence_for_topology(
            mol,
            cosmolkit_core::ValenceModel::RdkitLike,
        )
        .unwrap();
        super::get_experimental_torsions(
            mol, &rings, &valence, details, bonds, exp, small, macrocycle, basic, version, verbose,
        )
    }
    fn get_experimental_torsions_without_bonds(
        mol: &TopologyBlock,
        details: &mut CrystalFFDetails,
        exp: bool,
        small: bool,
        macrocycle: bool,
        basic: bool,
        version: u32,
        verbose: bool,
    ) -> Result<(), CrystalffTorsionPreferencesError> {
        get_experimental_torsions(
            mol,
            details,
            &mut Vec::new(),
            exp,
            small,
            macrocycle,
            basic,
            version,
            verbose,
        )
    }
    const VALID_BASE_PARAM_LINE: &str = "[C:1][C:2][C:3][C:4] 1 0 -1 0 1 0 -1 0 1 0 -1 0\n";
    const VALID_OVERRIDE_PARAM_LINE: &str =
        "[N:1][C:2][O:3][S:4] -1 1.1 1 1.2 -1 1.3 1 1.4 -1 1.5 1 1.6\n";

    fn conformer_fixture_path(relative_path: &str) -> std::path::PathBuf {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../..")
            .join("testdata/conformer/fixtures")
            .join(relative_path)
    }

    #[test]
    fn crystalff_exptorsionanglecollection_get_params_selects_version_1_table() {
        let params = ExpTorsionAngleCollection::get_params(1, false, false, "").unwrap();

        assert_eq!(params.param_data(), torsion_preferences_v1().unwrap());
        assert_ne!(params.param_data(), torsion_preferences_v2().unwrap());
    }

    #[test]
    fn crystalff_exptorsionanglecollection_get_params_selects_version_2_table() {
        let params = ExpTorsionAngleCollection::get_params(2, false, false, "").unwrap();

        assert_eq!(params.param_data(), torsion_preferences_v2().unwrap());
        assert_ne!(params.param_data(), torsion_preferences_v1().unwrap());
    }

    #[test]
    fn crystalff_exptorsionanglecollection_get_params_rejects_invalid_version() {
        let err = ExpTorsionAngleCollection::get_params(3, false, false, "")
            .expect_err("expected invalid-version failure");

        assert_eq!(
            err,
            CrystalffTorsionPreferencesError::InvalidExperimentalTorsionVersion { version: 3 }
        );
    }

    #[test]
    fn crystalff_exptorsionanglecollection_get_params_uses_explicit_param_override() {
        let custom = VALID_OVERRIDE_PARAM_LINE;

        let params = ExpTorsionAngleCollection::get_params(2, false, false, custom).unwrap();

        assert_eq!(params.param_data(), custom);
    }

    #[test]
    fn crystalff_exptorsionanglecollection_get_params_appends_small_ring_table() {
        let base = VALID_BASE_PARAM_LINE;

        let params = ExpTorsionAngleCollection::get_params(2, true, false, base).unwrap();

        assert_eq!(
            params.param_data(),
            format!("{base}{}", torsion_preferences_small_rings().unwrap())
        );
    }

    #[test]
    fn crystalff_exptorsionanglecollection_get_params_appends_macrocycle_table() {
        let base = VALID_BASE_PARAM_LINE;

        let params = ExpTorsionAngleCollection::get_params(2, false, true, base).unwrap();

        assert_eq!(
            params.param_data(),
            format!("{base}{}", torsion_preferences_macrocycles().unwrap())
        );
    }

    #[test]
    fn crystalff_exptorsionanglecollection_get_params_appends_both_optional_tables_in_source_order()
    {
        let base = VALID_BASE_PARAM_LINE;

        let params = ExpTorsionAngleCollection::get_params(2, true, true, base).unwrap();

        assert_eq!(
            params.param_data(),
            format!(
                "{base}{}{}",
                torsion_preferences_small_rings().unwrap(),
                torsion_preferences_macrocycles().unwrap()
            )
        );
    }

    #[test]
    fn crystalff_exptorsionanglecollection_constructor_skips_comment_lines() {
        let params = ExpTorsionAngleCollection::new(
            "# comment\n\
             [C:1][C:2][C:3][C:4] 1 0.1 -1 0.2 1 0.3 -1 0.4 1 0.5 -1 0.6\n",
        )
        .unwrap();

        assert_eq!(params.params().len(), 1);
        assert_eq!(params.params()[0].torsion_idx(), 0);
    }

    #[test]
    fn crystalff_exptorsionanglecollection_constructor_parses_tokens_and_assigns_torsion_indices() {
        let params = ExpTorsionAngleCollection::new(
            "[C:1][C:2][C:3][C:4] 1 0.1 -1 0.2 1 0.3 -1 0.4 1 0.5 -1 0.6\n\
             [N:1][C:2][O:3][S:4] -1 1.1 1 1.2 -1 1.3 1 1.4 -1 1.5 1 1.6\n",
        )
        .unwrap();

        assert_eq!(params.params().len(), 2);
        assert_eq!(params.params()[0].torsion_idx(), 0);
        assert_eq!(params.params()[1].torsion_idx(), 1);
        assert_eq!(params.params()[0].signs(), &[1, -1, 1, -1, 1, -1]);
        assert_eq!(
            params.params()[1].force_constants(),
            &[1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
        );
    }

    #[test]
    fn crystalff_query_default_bonds_do_not_match_azide_double_bond_chain() {
        let mol = fixture_hydrogens(fixture_smiles("[N-]=[N+]=N").expect("parse azide"));
        let query = compile_query_fixture("[*:1][X3,X2:2]=[X3,X2:3][*:4]").expect("build query");

        let matches = get_substruct_matches_with_params(
            &mol,
            &query,
            &SubstructMatchParams {
                max_matches: 1000,
                uniquify: true,
                use_chirality: false,
                specified_stereo_query_matches_unspecified: false,
                ..Default::default()
            },
        );

        assert!(
            matches.is_empty(),
            "implicit SMARTS bonds must remain single-or-aromatic, not match azide double bonds: {matches:?}"
        );
    }

    #[test]
    fn smarts_consumer_crystalff_unspecified_bond() {
        let query = parse_smarts("CC", &SmartsParseParams::default()).unwrap();
        let expanded = super::expand_crystalff_smarts_bonds(&query).unwrap();

        assert_eq!(expanded.len(), 1);
        assert_eq!(
            expanded[0].query,
            QueryNode::Predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Aromatic,
            ]))
        );
    }

    #[test]
    fn smarts_consumer_crystalff_bonds() {
        let query = parse_smarts("[#6]1~[#7]=[#8]-1", &SmartsParseParams::default()).unwrap();
        let expanded = super::expand_crystalff_smarts_bonds(&query).unwrap();

        assert_eq!(expanded.len(), query.num_bonds());
        assert_eq!(
            expanded
                .iter()
                .map(|bond| (bond.begin_atom_idx, bond.end_atom_idx, &bond.query))
                .collect::<Vec<_>>(),
            query
                .bonds()
                .iter()
                .map(|bond| (bond.begin().index(), bond.end().index(), bond.predicate()))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn smarts_consumer_crystalff_compile() {
        let params = ExpTorsionAngleCollection::new("CCCC 1 0 1 0 1 0 1 0 1 0 1 0\n").unwrap();

        let angle = &params.params()[0];
        assert_eq!(angle.smarts(), "CCCC");
        let query = angle.query_molecule().expect("canonical compiled query");
        assert_eq!(query.num_atoms(), 4);
        assert_eq!(query.num_bonds(), 3);
    }

    #[test]
    fn crystalff_exptorsionanglecollection_constructor_retains_smarts_parse_errors() {
        let params = ExpTorsionAngleCollection::new("[C:1 1 0 1 0 1 0 1 0 1 0 1 0\n").unwrap();

        let angle = &params.params()[0];
        assert_eq!(angle.smarts(), "[C:1");
        assert!(angle.query_molecule().is_err());
        assert_eq!(angle.idx(), [0; 4]);
    }

    #[test]
    fn smarts_consumer_crystalff_maps() {
        let params =
            ExpTorsionAngleCollection::new("[C:4][C:2][C:1][C:3] 1 0 1 0 1 0 1 0 1 0 1 0\n")
                .unwrap();

        assert_eq!(params.params()[0].idx(), [2, 1, 3, 0]);
    }

    #[test]
    fn crystalff_exptorsionanglecollection_constructor_ignores_out_of_range_atom_maps() {
        let params =
            ExpTorsionAngleCollection::new("[C:7][C:2][C:0][C:4] 1 0 1 0 1 0 1 0 1 0 1 0\n")
                .unwrap();

        assert_eq!(params.params()[0].idx(), [0, 1, 0, 3]);
    }

    #[test]
    fn crystalff_exptorsionanglecollection_constructor_maps_modeled_query_atoms() {
        let params =
            ExpTorsionAngleCollection::new("[C:1][$([C]):2][C:3][C:4] 1 0 1 0 1 0 1 0 1 0 1 0\n")
                .unwrap();
        let query_atom_count = params.params()[0]
            .query_molecule()
            .expect("query molecule")
            .num_atoms();

        assert_eq!(query_atom_count, 4);
        assert!(
            params.params()[0]
                .idx()
                .iter()
                .all(|idx| *idx < query_atom_count)
        );
    }

    #[test]
    fn crystalff_get_experimental_torsions_matches_linear_alkane_smarts() {
        let mol = fixture_smiles("CCCCC").expect("pentane");
        let mut details = CrystalFFDetails::default();
        let mut torsion_bonds = Vec::new();

        get_experimental_torsions(
            &mol,
            &mut details,
            &mut torsion_bonds,
            true,
            false,
            false,
            false,
            2,
            false,
        )
        .expect("pentane torsion preferences");

        assert_eq!(details.exp_torsion_atoms.len(), 2);
        assert_eq!(details.exp_torsion_angles.len(), 2);
        assert_eq!(torsion_bonds.len(), 2);
        assert_eq!(torsion_bonds[0].0, 1);
        assert_eq!(torsion_bonds[1].0, 2);
        assert_eq!(details.exp_torsion_atoms[0], vec![0, 1, 2, 3]);
        assert_eq!(details.exp_torsion_atoms[1], vec![1, 2, 3, 4]);
        assert_eq!(details.exp_torsion_angles[0].0.len(), 6);
        assert_eq!(details.exp_torsion_angles[0].1.len(), 6);
    }

    #[test]
    fn crystalff_x0_torsion_rule_rejects_ring_connected_central_atoms() {
        let cases = [
            (
                "O=C1Nc2ccc(Cl)cc2C(c2ccccc2Cl)=NC1O",
                vec![12, 11, 10, 18],
                vec![9, 10, 18, 19],
                13,
            ),
            (
                "CN1C2CCC1CC1(CN=C(c3cn(C)c4ccccc34)O1)C2",
                vec![12, 11, 10, 9],
                vec![8, 9, 10, 11],
                11,
            ),
        ];

        for (smiles, rejected_atoms, retained_atoms, expected_count) in cases {
            let mol = fixture_hydrogens(fixture_smiles(smiles).expect("parse"));
            let mut details = CrystalFFDetails::default();

            get_experimental_torsions_without_bonds(
                &mol,
                &mut details,
                true,
                false,
                true,
                true,
                2,
                false,
            )
            .expect("ETKDGv3 torsion preferences");

            assert_eq!(details.exp_torsion_atoms.len(), expected_count, "{smiles}");
            assert!(
                !details.exp_torsion_atoms.contains(&rejected_atoms),
                "x0 rule must reject the ring-connected central atom in {smiles}"
            );
            assert!(
                details.exp_torsion_atoms.contains(&retained_atoms),
                "{smiles}"
            );
        }
    }

    #[test]
    fn crystalff_get_experimental_torsions_applies_small_ring_and_macrocycle_tables() {
        let small_ring = fixture_smiles("C1COCC1").expect("tetrahydrofuran-like ring");
        let macrocycle = fixture_smiles("C1COCCCCCCC1").expect("macrocycle");
        let mut details = CrystalFFDetails::default();

        get_experimental_torsions_without_bonds(
            &small_ring,
            &mut details,
            true,
            false,
            false,
            false,
            1,
            false,
        )
        .expect("small ring without optional table");
        let small_ring_base_count = details.exp_torsion_atoms.len();
        assert_eq!(
            details.exp_torsion_atoms.len(),
            details.exp_torsion_angles.len()
        );

        get_experimental_torsions_without_bonds(
            &small_ring,
            &mut details,
            true,
            true,
            false,
            false,
            1,
            false,
        )
        .expect("small ring with optional table");
        assert!(!details.exp_torsion_atoms.is_empty());
        assert!(details.exp_torsion_atoms.len() >= small_ring_base_count);
        assert_eq!(
            details.exp_torsion_atoms.len(),
            details.exp_torsion_angles.len()
        );

        get_experimental_torsions_without_bonds(
            &macrocycle,
            &mut details,
            true,
            false,
            false,
            false,
            1,
            false,
        )
        .expect("macrocycle without optional table");
        let macrocycle_base_count = details.exp_torsion_atoms.len();
        assert_eq!(
            details.exp_torsion_atoms.len(),
            details.exp_torsion_angles.len()
        );

        get_experimental_torsions_without_bonds(
            &macrocycle,
            &mut details,
            true,
            false,
            true,
            false,
            1,
            false,
        )
        .expect("macrocycle with optional table");
        assert!(!details.exp_torsion_atoms.is_empty());
        assert!(details.exp_torsion_atoms.len() >= macrocycle_base_count);
        assert_eq!(
            details.exp_torsion_atoms.len(),
            details.exp_torsion_angles.len()
        );
    }

    #[test]
    fn crystalff_get_experimental_torsions_excludes_bridged_small_ring_bonds() {
        let mol = fixture_smiles("O[C@H]1C[C@H]2CC[C@]1(C)C2(C)C").expect("bridged ring");
        let mut details = CrystalFFDetails::default();

        get_experimental_torsions_without_bonds(
            &mol,
            &mut details,
            true,
            true,
            false,
            false,
            1,
            false,
        )
        .expect("bridged small ring path");

        assert_eq!(details.exp_torsion_atoms.len(), 0);
        assert_eq!(details.exp_torsion_angles.len(), 0);
    }

    #[test]
    fn crystalff_get_experimental_torsions_skips_fully_constrained_matches() {
        let mol = fixture_smiles("CCCC").expect("butane");
        let mut details = CrystalFFDetails {
            constrained_atoms: vec![true; mol.atoms.len()],
            ..CrystalFFDetails::default()
        };
        let mut torsion_bonds = Vec::new();

        get_experimental_torsions(
            &mol,
            &mut details,
            &mut torsion_bonds,
            true,
            false,
            false,
            false,
            1,
            false,
        )
        .expect("constrained torsion path");

        assert!(
            torsion_bonds.is_empty(),
            "unexpected torsion matches: {torsion_bonds:?}"
        );
        assert!(details.exp_torsion_atoms.is_empty());
        assert!(details.exp_torsion_angles.is_empty());
    }

    #[test]
    fn crystalff_get_experimental_torsions_adds_basic_knowledge_improper_terms() {
        let mol = fixture_smiles("CC(=O)O").expect("acetic acid");
        let mut details = CrystalFFDetails::default();

        get_experimental_torsions_without_bonds(
            &mol,
            &mut details,
            false,
            false,
            false,
            true,
            1,
            false,
        )
        .expect("basic knowledge improper terms");

        assert_eq!(details.improper_atoms.len(), 1);
        let improper = &details.improper_atoms[0];
        assert_eq!(improper.len(), 6);
        assert_eq!(improper[1], 1);
        assert_eq!(improper[4], 6);
        assert_eq!(improper[5], 1);
        assert!(improper[..4].contains(&0));
        assert!(improper[..4].contains(&2));
        assert!(improper[..4].contains(&3));
    }

    #[test]
    fn crystalff_get_experimental_torsions_injects_flat_ring_torsions_from_basic_knowledge() {
        let mol = fixture_smiles("c1ccccc1").expect("benzene");
        let mut details = CrystalFFDetails::default();

        get_experimental_torsions_without_bonds(
            &mol,
            &mut details,
            false,
            false,
            false,
            true,
            1,
            false,
        )
        .expect("flat ring torsion injection");

        assert_eq!(details.improper_atoms.len(), 0);
        assert_eq!(details.exp_torsion_atoms.len(), 6);
        assert_eq!(details.exp_torsion_angles.len(), 6);
        for (signs, fconsts) in &details.exp_torsion_angles {
            assert_eq!(signs, &vec![1, -1, 1, 1, 1, 1]);
            assert_eq!(fconsts, &vec![0.0, 100.0, 0.0, 0.0, 0.0, 0.0]);
        }
    }

    #[test]
    fn crystalff_get_experimental_torsions_rejects_empty_molecule() {
        let mol = TopologyBlock::default();
        let mut details = CrystalFFDetails::default();

        let err = get_experimental_torsions_without_bonds(
            &mol,
            &mut details,
            true,
            false,
            false,
            false,
            1,
            false,
        )
        .expect_err("empty molecule should fail");

        assert_eq!(err, CrystalffTorsionPreferencesError::EmptyMolecule);
    }

    #[test]
    fn crystalff_experimental_torsions_full_overload_covers_tables_and_policies() {
        let pentane = fixture_smiles("CCCCC").expect("pentane");
        let mut details = CrystalFFDetails {
            exp_torsion_atoms: vec![vec![99]],
            exp_torsion_angles: vec![(vec![99], vec![99.0])],
            improper_atoms: vec![vec![99]],
            ..CrystalFFDetails::default()
        };
        let mut torsion_bonds = vec![(
            99,
            vec![99],
            ExpTorsionAngleCollection::new("[C:1][C:2][C:3][C:4] 1 0 1 0 1 0 1 0 1 0 1 0\n")
                .unwrap()
                .params()[0]
                .clone(),
        )];

        get_experimental_torsions(
            &pentane,
            &mut details,
            &mut torsion_bonds,
            true,
            false,
            false,
            false,
            2,
            false,
        )
        .expect("pentane SMARTS torsion preferences");

        assert_eq!(
            details.exp_torsion_atoms,
            vec![vec![0, 1, 2, 3], vec![1, 2, 3, 4]]
        );
        assert_eq!(details.exp_torsion_angles.len(), 2);
        assert!(details.improper_atoms.is_empty());
        assert_eq!(
            torsion_bonds
                .iter()
                .map(|entry| entry.0)
                .collect::<Vec<_>>(),
            vec![1, 2]
        );

        let small_ring = fixture_smiles("C1COCC1").expect("tetrahydrofuran-like ring");
        get_experimental_torsions_without_bonds(
            &small_ring,
            &mut details,
            true,
            true,
            false,
            false,
            1,
            false,
        )
        .expect("small-ring optional torsion table");
        assert!(!details.exp_torsion_atoms.is_empty());
        assert_eq!(
            details.exp_torsion_atoms.len(),
            details.exp_torsion_angles.len()
        );

        let macrocycle = fixture_smiles("C1COCCCCCCC1").expect("macrocycle");
        get_experimental_torsions_without_bonds(
            &macrocycle,
            &mut details,
            true,
            false,
            true,
            false,
            1,
            false,
        )
        .expect("macrocycle optional torsion table");
        assert!(!details.exp_torsion_atoms.is_empty());
        assert_eq!(
            details.exp_torsion_atoms.len(),
            details.exp_torsion_angles.len()
        );

        let bridged = fixture_smiles("O[C@H]1C[C@H]2CC[C@]1(C)C2(C)C").expect("bridged ring");
        get_experimental_torsions_without_bonds(
            &bridged,
            &mut details,
            true,
            true,
            false,
            false,
            1,
            false,
        )
        .expect("bridged-ring exclusion");
        assert!(details.exp_torsion_atoms.is_empty());
        assert!(details.exp_torsion_angles.is_empty());

        let butane = fixture_smiles("CCCC").expect("butane");
        details.constrained_atoms = vec![true; butane.atoms.len()];
        get_experimental_torsions(
            &butane,
            &mut details,
            &mut torsion_bonds,
            true,
            false,
            false,
            false,
            1,
            false,
        )
        .expect("fully constrained torsion skip");
        assert!(
            torsion_bonds.is_empty(),
            "unexpected torsion matches: {torsion_bonds:?}"
        );
        assert!(details.exp_torsion_atoms.is_empty());

        let err = get_experimental_torsions_without_bonds(
            &butane,
            &mut details,
            true,
            false,
            false,
            false,
            3,
            false,
        )
        .expect_err("unsupported ETversion should fail");
        assert_eq!(
            err,
            CrystalffTorsionPreferencesError::InvalidExperimentalTorsionVersion { version: 3 }
        );
    }

    #[test]
    fn crystalff_rdkit_simple_torsion_fixture_matches_rdkit_get_experimental_torsions() {
        let path = conformer_fixture_path("rdkit/test_data/simple_torsion.etkdg.mol");
        let text = std::fs::read_to_string(&path).expect("read RDKit simple_torsion.etkdg fixture");
        let mol = fixture_molblock(&text);

        let mut details = CrystalFFDetails::default();
        let mut torsion_bonds = Vec::new();
        get_experimental_torsions(
            &mol,
            &mut details,
            &mut torsion_bonds,
            true,
            false,
            false,
            false,
            2,
            false,
        )
        .expect("simple_torsion.etkdg torsion preferences");

        assert_eq!(torsion_bonds.len(), 1);
        assert_eq!(torsion_bonds[0].0, 1);
        assert_eq!(torsion_bonds[0].1, vec![0, 1, 2, 3]);
        assert_eq!(torsion_bonds[0].2.torsion_idx(), 229);
        assert_eq!(
            torsion_bonds[0].2.smarts(),
            "[!#1:1][CX4H2:2]!@;-[CX4H2:3][!#1:4]"
        );
        assert_eq!(details.exp_torsion_atoms, vec![vec![0, 1, 2, 3]]);
        assert_eq!(
            details.exp_torsion_angles,
            vec![(vec![1, 1, 1, 1, 1, 1], vec![0.0, 0.0, 4.0, 0.0, 0.0, 0.0])]
        );
    }

    #[test]
    fn crystalff_rdkit_smallring_fixture_matches_rdkit_sr_etkdgv3_experimental_torsions() {
        let path = conformer_fixture_path("rdkit/test_data/simple_torsion.smallring.etkdgv3.mol");
        let text =
            std::fs::read_to_string(&path).expect("read RDKit simple_torsion.smallring fixture");
        let mol = fixture_molblock(&text);

        let mut details = CrystalFFDetails::default();
        let mut torsion_bonds = Vec::new();
        get_experimental_torsions(
            &mol,
            &mut details,
            &mut torsion_bonds,
            true,
            true,
            false,
            false,
            2,
            false,
        )
        .expect("simple_torsion.smallring torsion preferences");

        let actual: Vec<(usize, Vec<usize>, usize, String)> = torsion_bonds
            .iter()
            .map(|(bond_idx, atom_indices, angle)| {
                (
                    *bond_idx,
                    atom_indices.clone(),
                    angle.torsion_idx(),
                    angle.smarts().to_owned(),
                )
            })
            .collect();
        let expected = vec![
            (
                1,
                vec![0, 1, 2, 3],
                392,
                "[!#1;r{5-8}:1]@[CX4;r{5-8}:2]@;-[CX4;r{5-8}:3]@[!#1;r{5-8}:4]".to_owned(),
            ),
            (
                4,
                vec![0, 5, 4, 3],
                392,
                "[!#1;r{5-8}:1]@[CX4;r{5-8}:2]@;-[CX4;r{5-8}:3]@[!#1;r{5-8}:4]".to_owned(),
            ),
            (
                5,
                vec![1, 0, 5, 4],
                392,
                "[!#1;r{5-8}:1]@[CX4;r{5-8}:2]@;-[CX4;r{5-8}:3]@[!#1;r{5-8}:4]".to_owned(),
            ),
            (
                2,
                vec![1, 2, 3, 4],
                392,
                "[!#1;r{5-8}:1]@[CX4;r{5-8}:2]@;-[CX4;r{5-8}:3]@[!#1;r{5-8}:4]".to_owned(),
            ),
            (
                0,
                vec![2, 1, 0, 5],
                392,
                "[!#1;r{5-8}:1]@[CX4;r{5-8}:2]@;-[CX4;r{5-8}:3]@[!#1;r{5-8}:4]".to_owned(),
            ),
            (
                3,
                vec![2, 3, 4, 5],
                392,
                "[!#1;r{5-8}:1]@[CX4;r{5-8}:2]@;-[CX4;r{5-8}:3]@[!#1;r{5-8}:4]".to_owned(),
            ),
        ];
        assert_eq!(actual, expected);
        assert_eq!(
            details.exp_torsion_atoms,
            vec![
                vec![0, 1, 2, 3],
                vec![0, 5, 4, 3],
                vec![1, 0, 5, 4],
                vec![1, 2, 3, 4],
                vec![2, 1, 0, 5],
                vec![2, 3, 4, 5],
            ]
        );
        assert_eq!(
            details.exp_torsion_angles,
            vec![
                (vec![1, 1, 1, 1, 1, 1], vec![0.0, 0.0, 30.0, 0.0, 0.0, 0.0]),
                (vec![1, 1, 1, 1, 1, 1], vec![0.0, 0.0, 30.0, 0.0, 0.0, 0.0]),
                (vec![1, 1, 1, 1, 1, 1], vec![0.0, 0.0, 30.0, 0.0, 0.0, 0.0]),
                (vec![1, 1, 1, 1, 1, 1], vec![0.0, 0.0, 30.0, 0.0, 0.0, 0.0]),
                (vec![1, 1, 1, 1, 1, 1], vec![0.0, 0.0, 30.0, 0.0, 0.0, 0.0]),
                (vec![1, 1, 1, 1, 1, 1], vec![0.0, 0.0, 30.0, 0.0, 0.0, 0.0]),
            ]
        );
    }

    #[test]
    fn crystalff_simple_torsion_query_matches_fixture_like_rdkit() {
        let path = conformer_fixture_path("rdkit/test_data/simple_torsion.etkdg.mol");
        let text = std::fs::read_to_string(&path).expect("read RDKit simple_torsion.etkdg fixture");
        let mol = fixture_molblock(&text);
        let query = compile_query_fixture("[!#1:1][CX4H2:2]!@;-[CX4H2:3][!#1:4]")
            .expect("build CrystalFF query");
        let matches = get_substruct_matches_with_params(
            &mol,
            &query,
            &SubstructMatchParams {
                max_matches: 1000,
                uniquify: false,
                use_chirality: false,
                specified_stereo_query_matches_unspecified: false,
                ..Default::default()
            },
        );

        assert_eq!(
            matches
                .iter()
                .map(|m| m.atom_mapping.clone())
                .collect::<Vec<_>>(),
            vec![vec![0, 1, 2, 3], vec![3, 2, 1, 0]]
        );
    }

    #[test]
    fn crystalff_aliphatic_carbon_query_does_not_match_row34_aromatic_carbon() {
        let mol =
            fixture_hydrogens(fixture_smiles("CCCC1CCC(c2ccc(OCC)cc2)CC1").expect("parse row34"));
        let query = compile_query_fixture("[C:1][CX4:2]!@;-[CX3:3][c:4]").expect("build query");

        let matches = get_substruct_matches_with_params(
            &mol,
            &query,
            &SubstructMatchParams {
                max_matches: 1000,
                uniquify: false,
                use_chirality: false,
                specified_stereo_query_matches_unspecified: false,
                ..Default::default()
            },
        );

        assert!(
            !matches.iter().any(|m| m.atom_mapping == vec![5, 6, 7, 8]),
            "aliphatic [CX3] must not match aromatic carbon in row34: {matches:?}"
        );
    }

    #[test]
    fn crystalff_carboxylate_query_does_not_match_neutral_row20_carboxylic_acid() {
        let mol = fixture_hydrogens(fixture_smiles("N[C@@H](C)C(=O)O").expect("parse row20"));
        let query = compile_query_fixture("[O:1]=[C:2]([O-])!@;-[CX4H1:3][H:4]")
            .expect("build carboxylate query");

        let matches = get_substruct_matches_with_params(
            &mol,
            &query,
            &SubstructMatchParams {
                max_matches: 1000,
                uniquify: false,
                use_chirality: false,
                specified_stereo_query_matches_unspecified: false,
                ..Default::default()
            },
        );

        assert!(
            matches.is_empty(),
            "carboxylate SMARTS must not match neutral row20 acid: {matches:?}"
        );
    }

    #[test]
    fn crystalff_recursive_carbonyl_query_matches_row57_hydrazide_torsion_like_rdkit() {
        let mol = fixture_hydrogens(fixture_smiles("CCCCCC(=O)NNC(N)=S").expect("parse row57"));
        let query = compile_query_fixture("[$(C=O):1][NX3:2]!@;-[!#1:3][!#1:4]")
            .expect("build recursive carbonyl query");

        let matches = get_substruct_matches_with_params(
            &mol,
            &query,
            &SubstructMatchParams {
                max_matches: 1000,
                uniquify: false,
                use_chirality: false,
                specified_stereo_query_matches_unspecified: false,
                ..Default::default()
            },
        );

        assert!(
            matches.iter().any(|m| m.atom_mapping == vec![5, 7, 8, 9]),
            "recursive carbonyl SMARTS must match row57 torsion [5,7,8,9]: {matches:?}"
        );
    }

    #[test]
    #[ignore = "debug helper for row34 torsion SMARTS parity investigation"]
    fn debug_row34_hydrogen_torsion_query_matches() {
        let smarts = "[a:1][c:2]!@;-[CX4H1:3][H:4]";
        let mol =
            fixture_hydrogens(fixture_smiles("CCCC1CCC(c2ccc(OCC)cc2)CC1").expect("parse row34"));
        let query = compile_query_fixture(smarts).expect("build CrystalFF query");

        let matches = get_substruct_matches_with_params(
            &mol,
            &query,
            &SubstructMatchParams {
                max_matches: 1000,
                uniquify: false,
                use_chirality: false,
                specified_stereo_query_matches_unspecified: false,
                ..Default::default()
            },
        );

        println!("row34_debug_smarts={smarts}");
        println!(
            "row34_debug_matches={:?}",
            matches
                .iter()
                .map(|m| m.atom_mapping.clone())
                .collect::<Vec<_>>()
        );

        for (query_atom_idx, query_atom) in query.atoms().iter().enumerate() {
            let matching_atoms: Vec<_> = mol
                .atoms
                .iter()
                .enumerate()
                .filter_map(|(mol_atom_idx, mol_atom)| {
                    let query_node = query_atom.predicate();
                    Some(query_node).and_then(|query_node| {
                        atom_matches_query(mol_atom, query_node, &test_target(&mol))
                            .unwrap()
                            .then_some(mol_atom_idx)
                    })
                })
                .collect();
            println!(
                "row34_query_atom_matches query_atom_idx={query_atom_idx} query={:?} matches={matching_atoms:?}",
                query_atom.predicate()
            );
        }

        for (query_bond_idx, query_bond) in query.bonds().iter().enumerate() {
            let matching_bonds: Vec<_> = mol
                .bonds
                .iter()
                .enumerate()
                .filter_map(|(mol_bond_idx, mol_bond)| {
                    let query_node = query_bond.predicate();
                    Some(query_node).and_then(|query_node| {
                        bond_matches_query(mol_bond, query_node, &test_target(&mol))
                            .unwrap()
                            .then_some(mol_bond_idx)
                    })
                })
                .collect();
            println!(
                "row34_query_bond_matches query_bond_idx={query_bond_idx} query={:?} matches={matching_bonds:?}",
                query_bond.predicate()
            );
        }
    }

    #[test]
    #[ignore = "debug helper for row20 carboxylate torsion SMARTS parity investigation"]
    fn debug_row20_carboxylate_torsion_query_matches() {
        let smarts = "[O:1]=[C:2]([O-])!@;-[CX4H1:3][H:4]";
        let mol = fixture_hydrogens(fixture_smiles("N[C@@H](C)C(=O)O").expect("parse row20"));
        let query = compile_query_fixture(smarts).expect("build CrystalFF query");

        let matches = get_substruct_matches_with_params(
            &mol,
            &query,
            &SubstructMatchParams {
                max_matches: 1000,
                uniquify: false,
                use_chirality: false,
                specified_stereo_query_matches_unspecified: false,
                ..Default::default()
            },
        );

        println!("row20_debug_smarts={smarts}");
        println!(
            "row20_debug_matches={:?}",
            matches
                .iter()
                .map(|m| m.atom_mapping.clone())
                .collect::<Vec<_>>()
        );

        for (query_atom_idx, query_atom) in query.atoms().iter().enumerate() {
            let matching_atoms: Vec<_> = mol
                .atoms
                .iter()
                .enumerate()
                .filter_map(|(mol_atom_idx, mol_atom)| {
                    let query_node = query_atom.predicate();
                    Some(query_node).and_then(|query_node| {
                        atom_matches_query(mol_atom, query_node, &test_target(&mol))
                            .unwrap()
                            .then_some(mol_atom_idx)
                    })
                })
                .collect();
            println!(
                "row20_query_atom_matches query_atom_idx={query_atom_idx} query={:?} matches={matching_atoms:?}",
                query_atom.predicate()
            );
        }

        for (query_bond_idx, query_bond) in query.bonds().iter().enumerate() {
            let matching_bonds: Vec<_> = mol
                .bonds
                .iter()
                .enumerate()
                .filter_map(|(mol_bond_idx, mol_bond)| {
                    let query_node = query_bond.predicate();
                    Some(query_node).and_then(|query_node| {
                        bond_matches_query(mol_bond, query_node, &test_target(&mol))
                            .unwrap()
                            .then_some(mol_bond_idx)
                    })
                })
                .collect();
            println!(
                "row20_query_bond_matches query_bond_idx={query_bond_idx} query={:?} matches={matching_bonds:?}",
                query_bond.predicate()
            );
        }

        for (mol_atom_idx, atom) in mol.atoms.iter().enumerate() {
            println!(
                "row20_mol_atom idx={mol_atom_idx} atomic_num={} formal_charge={} explicit_hs={} aromatic={} hyb={:?} atom_map={:?}",
                atom.atomic_number(),
                atom.formal_charge(),
                atom.explicit_hydrogens(),
                atom.is_aromatic(),
                atom.hybridization(),
                atom.atom_map()
            );
        }
    }

    #[test]
    #[ignore = "debug helper for row89 recursive carbonyl torsion SMARTS parity investigation"]
    fn debug_row89_recursive_carbonyl_torsion_query_matches() {
        let smarts = "[$(C=O):1][NX3:2]!@;-[!#1:3][!#1:4]";
        let mol = fixture_hydrogens(
            fixture_smiles("COC(=O)c4ccccc4(NC(=O)n3c1ccccc1sc2ccccc23)").expect("parse row89"),
        );
        let query = compile_query_fixture(smarts).expect("build CrystalFF query");

        let matches = get_substruct_matches_with_params(
            &mol,
            &query,
            &SubstructMatchParams {
                max_matches: 1000,
                uniquify: false,
                use_chirality: false,
                specified_stereo_query_matches_unspecified: false,
                ..Default::default()
            },
        );

        println!("row89_debug_smarts={smarts}");
        println!(
            "row89_debug_query_idx={:?}",
            super::map_pattern_atom_indices(&query)
        );
        println!(
            "row89_debug_matches={:?}",
            matches
                .iter()
                .map(|m| m.atom_mapping.clone())
                .collect::<Vec<_>>()
        );

        for (query_atom_idx, query_atom) in query.atoms().iter().enumerate() {
            let matching_atoms: Vec<_> = mol
                .atoms
                .iter()
                .enumerate()
                .filter_map(|(mol_atom_idx, mol_atom)| {
                    let query_node = query_atom.predicate();
                    Some(query_node).and_then(|query_node| {
                        atom_matches_query(mol_atom, query_node, &test_target(&mol))
                            .unwrap()
                            .then_some(mol_atom_idx)
                    })
                })
                .collect();
            println!(
                "row89_query_atom_matches query_atom_idx={query_atom_idx} atom_map={:?} query={:?} matches={matching_atoms:?}",
                query_atom.atom_map(),
                query_atom.predicate()
            );
        }

        for (query_bond_idx, query_bond) in query.bonds().iter().enumerate() {
            let matching_bonds: Vec<_> = mol
                .bonds
                .iter()
                .enumerate()
                .filter_map(|(mol_bond_idx, mol_bond)| {
                    let query_node = query_bond.predicate();
                    Some(query_node).and_then(|query_node| {
                        bond_matches_query(mol_bond, query_node, &test_target(&mol))
                            .unwrap()
                            .then_some(mol_bond_idx)
                    })
                })
                .collect();
            println!(
                "row89_query_bond_matches query_bond_idx={query_bond_idx} query={:?} matches={matching_bonds:?}",
                query_bond.predicate()
            );
        }
    }

    #[test]
    #[ignore = "debug helper for row61 torsion_idx=349 substructure parity investigation"]
    fn debug_row61_torsion349_query_matches() {
        let smarts = "[cH0:1][c:2]([cH0])!@;-[CX3:3]=[O:4]";
        let mol = fixture_hydrogens(
            fixture_smiles("COC(=O)c1cccc([N+](=O)[O-])c1C(=O)OC").expect("parse row61"),
        );
        let query = compile_query_fixture(smarts).expect("build CrystalFF query");

        let matches = get_substruct_matches_with_params(
            &mol,
            &query,
            &SubstructMatchParams {
                max_matches: 1000,
                uniquify: false,
                use_chirality: false,
                specified_stereo_query_matches_unspecified: false,
                ..Default::default()
            },
        );

        println!("row61_debug_smarts={smarts}");
        println!(
            "row61_debug_query_idx={:?}",
            super::map_pattern_atom_indices(&query)
        );
        println!(
            "row61_debug_matches={:?}",
            matches
                .iter()
                .map(|m| m.atom_mapping.clone())
                .collect::<Vec<_>>()
        );

        for (query_atom_idx, query_atom) in query.atoms().iter().enumerate() {
            let matching_atoms: Vec<_> = mol
                .atoms
                .iter()
                .enumerate()
                .filter_map(|(mol_atom_idx, mol_atom)| {
                    let query_node = query_atom.predicate();
                    Some(query_node).and_then(|query_node| {
                        atom_matches_query(mol_atom, query_node, &test_target(&mol))
                            .unwrap()
                            .then_some(mol_atom_idx)
                    })
                })
                .collect();
            println!(
                "row61_query_atom_matches query_atom_idx={query_atom_idx} atom_map={:?} query={:?} matches={matching_atoms:?}",
                query_atom.atom_map(),
                query_atom.predicate()
            );
        }

        for (query_bond_idx, query_bond) in query.bonds().iter().enumerate() {
            let matching_bonds: Vec<_> = mol
                .bonds
                .iter()
                .enumerate()
                .filter_map(|(mol_bond_idx, mol_bond)| {
                    let query_node = query_bond.predicate();
                    Some(query_node).and_then(|query_node| {
                        bond_matches_query(mol_bond, query_node, &test_target(&mol))
                            .unwrap()
                            .then_some(mol_bond_idx)
                    })
                })
                .collect();
            println!(
                "row61_query_bond_matches query_bond_idx={query_bond_idx} query={:?} matches={matching_bonds:?}",
                query_bond.predicate()
            );
        }
    }

    #[test]
    fn crystalff_etv2_fixed_fixture_first_boundary_diagnostic() {
        let path = conformer_fixture_path("rdkit/test_data/torsion.etkdg.v2.mol");
        let mol = fixture_molblock(&std::fs::read_to_string(path).unwrap());
        let rings = cosmolkit_core::symmetrized_sssr(&mol, &Default::default()).unwrap();
        let valence = cosmolkit_core::assign_valence_for_topology(
            &mol,
            cosmolkit_core::ValenceModel::RdkitLike,
        )
        .unwrap();
        for (i, atom) in mol.atoms.iter().enumerate() {
            println!(
                "etv2_atom id={i} z={} arom={} hybrid={:?} h={} valence={}",
                atom.atomic_number(),
                atom.is_aromatic(),
                atom.hybridization(),
                i32::from(atom.explicit_hydrogens()) + valence.implicit_hydrogens[i],
                valence.explicit_valence[i]
            );
        }
        println!(
            "etv2_bonds {:?}",
            mol.bonds
                .iter()
                .map(|b| (
                    b.id().index(),
                    b.begin().index(),
                    b.end().index(),
                    b.order(),
                    b.is_aromatic()
                ))
                .collect::<Vec<_>>()
        );
        println!("etv2_rings {:?}", rings);
        let mut details = CrystalFFDetails::default();
        let mut bonds = Vec::new();
        get_experimental_torsions(
            &mol,
            &mut details,
            &mut bonds,
            true,
            false,
            false,
            false,
            2,
            false,
        )
        .unwrap();
        for (bond, atoms, params) in &bonds {
            println!(
                "etv2_torsion bond={bond} atoms={atoms:?} idx={} smarts={} signs={:?} force={:?}",
                params.torsion_idx(),
                params.smarts(),
                params.signs(),
                params.force_constants()
            );
        }
        assert_eq!(
            bonds.len(),
            1,
            "fixed upstream GetExperimentalTorsions boundary"
        );
        assert_eq!(
            (bonds[0].0, bonds[0].1.clone(), bonds[0].2.torsion_idx()),
            (6, vec![0, 6, 7, 8], 39)
        );
        assert_eq!(
            details.exp_torsion_angles,
            vec![(vec![1, -1, 1, 1, 1, 1], vec![0.0, 7.2, 0.0, 0.0, 0.0, 0.0])]
        );
    }
}
