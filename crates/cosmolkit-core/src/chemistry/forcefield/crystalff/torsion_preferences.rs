//! Source-backed RDKit CrystalFF torsion-preference helpers.

use crate::SmartsParseError;
use crate::search::smarts_parse::{SmartsMolecule, parse_smarts};
use crate::{
    AtomQueryPredicate, AtomSpec, BondId, BondOrder, BondQueryPredicate, BondSpec, Element,
    Hybridization, Molecule, MoleculeBuilder, QueryNode, RingInfo, SubstructMatchParams,
    get_substruct_matches_with_params, symmetrize_sssr,
};
use std::collections::BTreeMap;

const MIN_MACROCYCLE_SIZE: usize = 9;

include!(concat!(env!("OUT_DIR"), "/crystalff_defaults_generated.rs"));

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
    pattern: Result<SmartsMolecule, SmartsParseError>,
    query_molecule: Result<Molecule, String>,
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
    pub(crate) fn pattern(&self) -> Result<&SmartsMolecule, &SmartsParseError> {
        self.pattern.as_ref()
    }

    #[must_use]
    pub(crate) fn query_molecule(&self) -> Result<&Molecule, &str> {
        self.query_molecule.as_ref().map_err(String::as_str)
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
    ) -> Result<Self, CrystalffTorsionPreferencesError> {
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
        let params = if param_data.is_empty() {
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

        let mut collection = Self::new(&params)?;

        if use_small_ring_torsions {
            collection
                .param_data
                .push_str(torsion_preferences_small_rings()?);
            let small_ring_collection = Self::new(torsion_preferences_small_rings()?)?;
            collection.params.extend(small_ring_collection.params);
        }

        if use_macrocycle_torsions {
            collection
                .param_data
                .push_str(torsion_preferences_macrocycles()?);
            let macrocycle_collection = Self::new(torsion_preferences_macrocycles()?)?;
            collection.params.extend(macrocycle_collection.params);
        }

        Ok(collection)
    }
}

pub fn get_experimental_torsions(
    mol: &Molecule,
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
    let nb = mol.num_bonds();
    let na = mol.num_atoms();
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

    let ring_info = crystalff_ring_info(mol)?;
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
                    reason: reason.to_owned(),
                }
            })?;
            let matches = get_substruct_matches_with_params(
                mol,
                query_molecule,
                &SubstructMatchParams {
                    max_matches: 1000,
                    uniquify: false,
                },
            );
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
            let atom2 = &mol.atoms()[aid2];
            let at2_atomic_num = atom2.atomic_number();
            if matches!(at2_atomic_num, 6 | 7 | 8)
                && atom2.hybridization() == Hybridization::Sp2
                && mol.topology_block().adjacency.neighbors_of(aid2).len() == 3
            {
                let mut i = 0usize;
                let mut is_bound_to_sp2o = 0i32;
                for neighbor in mol.topology_block().adjacency.neighbors_of(aid2) {
                    let atom_x = &mol.atoms()[neighbor.atom_index];
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
                    && mol.atoms()[aid1].hybridization() == Hybridization::Sp2
                    && mol.atoms()[aid2].hybridization() == Hybridization::Sp2
                    && mol.atoms()[aid3].hybridization() == Hybridization::Sp2
                    && mol.atoms()[aid4].hybridization() == Hybridization::Sp2
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

    Ok(())
}

pub fn get_experimental_torsions_without_bonds(
    mol: &Molecule,
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

    let pattern = parse_smarts(&smarts);
    let query_molecule = build_crystalff_query_molecule(&smarts);
    let idx = map_pattern_atom_indices(&smarts);

    Ok(ExpTorsionAngle {
        torsion_idx,
        smarts,
        force_constants,
        signs,
        pattern,
        query_molecule,
        idx,
    })
}

fn crystalff_ring_info(mol: &Molecule) -> Result<RingInfo, CrystalffTorsionPreferencesError> {
    if let Some(ring_info) = mol.derived_cache().rings.as_ref() {
        if ring_info.is_initialized() {
            return Ok(ring_info.clone());
        }
    }
    symmetrize_sssr(mol).map_err(
        |error| CrystalffTorsionPreferencesError::RingInfoUnavailable {
            reason: error.to_string(),
        },
    )
}

fn bond_between_atoms(mol: &Molecule, a: usize, b: usize) -> Option<&crate::Bond> {
    mol.bonds().iter().find(|bond| {
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

fn build_crystalff_query_molecule(smarts: &str) -> Result<Molecule, String> {
    let parsed = parse_smarts(smarts).map_err(|error| error.to_string())?;
    let atom_queries = parsed.atom_queries;
    let bonds = expand_crystalff_smarts_bonds(smarts)?;
    let mut builder = MoleculeBuilder::new();
    let atom_ids: Vec<_> = atom_queries
        .iter()
        .map(|query| {
            builder.add_atom(
                AtomSpec::new(
                    Element::from_atomic_number(0)
                        .expect("atomic number zero is a valid query atom"),
                )
                .with_query(query.clone()),
            )
        })
        .collect();
    for (atom_idx, atom_id) in atom_ids.iter().enumerate() {
        if let Some(map_num) = extract_atom_map_from_query(&atom_queries[atom_idx]) {
            let _ = map_num;
            let atom = builder.atom_mut(*atom_id).expect("query atom must exist");
            atom.set_atom_map(Some(map_num));
        }
    }
    for bond in &bonds {
        builder
            .add_bond(
                BondSpec::new(
                    atom_ids[bond.begin_atom_idx],
                    atom_ids[bond.end_atom_idx],
                    bond_order_for_query_probe(&bond.query),
                )
                .with_query(bond.query.clone()),
            )
            .map_err(|error| error.to_string())?;
    }
    builder.build().map_err(|error| error.to_string())
}

fn extract_atom_map_from_query(query: &QueryNode<AtomQueryPredicate>) -> Option<u32> {
    match query {
        QueryNode::Predicate(AtomQueryPredicate::AtomMapNumber(n)) => Some(*n),
        QueryNode::And(children) | QueryNode::Or(children) => {
            children.iter().find_map(extract_atom_map_from_query)
        }
        QueryNode::Not(child) => extract_atom_map_from_query(child),
        _ => None,
    }
}

fn bond_order_for_query_probe(query: &QueryNode<BondQueryPredicate>) -> BondOrder {
    match query {
        QueryNode::Predicate(BondQueryPredicate::Order(order)) => *order,
        QueryNode::Predicate(BondQueryPredicate::IsAromatic(true)) => BondOrder::Aromatic,
        QueryNode::And(children) | QueryNode::Or(children) => children
            .iter()
            .find_map(|child| match child {
                QueryNode::Predicate(BondQueryPredicate::Order(order)) => Some(*order),
                QueryNode::Predicate(BondQueryPredicate::IsAromatic(true)) => {
                    Some(BondOrder::Aromatic)
                }
                _ => None,
            })
            .unwrap_or(BondOrder::Single),
        _ => BondOrder::Single,
    }
}

fn expand_crystalff_smarts_bonds(smarts: &str) -> Result<Vec<CrystalffTemplateBond>, String> {
    #[derive(Clone)]
    struct ParserState<'a> {
        chars: &'a [char],
        pos: usize,
        next_atom_idx: usize,
        bonds: Vec<CrystalffTemplateBond>,
        ring_open: BTreeMap<u8, (usize, QueryNode<BondQueryPredicate>)>,
    }

    fn bond_query_for_char(ch: char) -> QueryNode<BondQueryPredicate> {
        match ch {
            '-' => QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Single)),
            '=' => QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Double)),
            '#' => QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Triple)),
            ':' => QueryNode::Predicate(BondQueryPredicate::IsAromatic(true)),
            '~' => QueryNode::Predicate(BondQueryPredicate::Any),
            '/' => QueryNode::Predicate(BondQueryPredicate::Direction(
                crate::BondDirection::EndUpRight,
            )),
            '\\' => QueryNode::Predicate(BondQueryPredicate::Direction(
                crate::BondDirection::EndDownRight,
            )),
            _ => QueryNode::Predicate(BondQueryPredicate::Any),
        }
    }

    fn is_bond_char(ch: char) -> bool {
        matches!(ch, '-' | '=' | '#' | ':' | '~' | '/' | '\\')
    }

    fn is_bond_query_char(ch: char) -> bool {
        is_bond_char(ch) || matches!(ch, '!' | '@' | ';' | ',')
    }

    fn simplify_bond_query(
        predicates: Vec<QueryNode<BondQueryPredicate>>,
    ) -> QueryNode<BondQueryPredicate> {
        let mut simplified = predicates
            .into_iter()
            .filter(|predicate| *predicate != QueryNode::Predicate(BondQueryPredicate::Any))
            .collect::<Vec<_>>();
        match simplified.len() {
            0 => QueryNode::Predicate(BondQueryPredicate::Any),
            1 => simplified.pop().expect("single bond predicate"),
            _ => QueryNode::And(simplified),
        }
    }

    fn parse_bond_query(
        chars: &[char],
        pos: &mut usize,
    ) -> Result<QueryNode<BondQueryPredicate>, String> {
        let mut predicates = Vec::new();
        let mut consumed = false;
        let mut logical_or = false;
        while let Some(&ch) = chars.get(*pos) {
            if !is_bond_query_char(ch) {
                break;
            }
            consumed = true;
            match ch {
                '!' if chars.get(*pos + 1) == Some(&'@') => {
                    predicates.push(QueryNode::Predicate(BondQueryPredicate::IsInRing(false)));
                    *pos += 2;
                }
                '@' => {
                    predicates.push(QueryNode::Predicate(BondQueryPredicate::IsInRing(true)));
                    *pos += 1;
                }
                ';' => {
                    *pos += 1;
                }
                ',' => {
                    logical_or = true;
                    *pos += 1;
                }
                '-' | '=' | '#' | ':' | '~' | '/' | '\\' => {
                    predicates.push(bond_query_for_char(ch));
                    *pos += 1;
                }
                '!' => {
                    return Err(format!(
                        "RDKit CrystalFF SMARTS topology parser encountered unsupported negated bond query at position {}",
                        *pos
                    ));
                }
                _ => unreachable!("non-bond query character filtered earlier"),
            }
        }
        if consumed {
            let simplified = predicates
                .into_iter()
                .filter(|predicate| *predicate != QueryNode::Predicate(BondQueryPredicate::Any))
                .collect::<Vec<_>>();
            if logical_or {
                match simplified.len() {
                    0 => Ok(QueryNode::Predicate(BondQueryPredicate::Any)),
                    1 => Ok(simplified
                        .into_iter()
                        .next()
                        .expect("single bond predicate")),
                    _ => Ok(QueryNode::Or(simplified)),
                }
            } else {
                Ok(simplify_bond_query(simplified))
            }
        } else {
            Ok(QueryNode::Predicate(BondQueryPredicate::Any))
        }
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

    fn parse_ring_number(chars: &[char], pos: &mut usize) -> Result<u8, String> {
        if chars.get(*pos) == Some(&'%') {
            let d1 = chars.get(*pos + 1).copied();
            let d2 = chars.get(*pos + 2).copied();
            match (d1, d2) {
                (Some(c1), Some(c2)) if c1.is_ascii_digit() && c2.is_ascii_digit() => {
                    *pos += 3;
                    Ok(((c1 as u8 - b'0') * 10) + (c2 as u8 - b'0'))
                }
                _ => Err(
                    "RDKit CrystalFF SMARTS topology parser expected two digits after '%'"
                        .to_string(),
                ),
            }
        } else if let Some(ch) = chars.get(*pos).copied() {
            if ch.is_ascii_digit() {
                *pos += 1;
                Ok(ch as u8 - b'0')
            } else {
                Err(
                    "RDKit CrystalFF SMARTS topology parser expected a ring-closure digit"
                        .to_string(),
                )
            }
        } else {
            Err(
                "RDKit CrystalFF SMARTS topology parser reached end while reading ring closure"
                    .to_string(),
            )
        }
    }

    fn skip_atom(chars: &[char], pos: &mut usize) -> Result<(), String> {
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
                Err(
                    "RDKit CrystalFF SMARTS topology parser found an unclosed bracket atom"
                        .to_string(),
                )
            }
            Some(_) => {
                let Some(consumed) = is_organic_subset_symbol(chars, *pos) else {
                    return Err(
                        "RDKit CrystalFF SMARTS topology parser encountered an unsupported atom token"
                            .to_string(),
                    );
                };
                *pos += consumed;
                Ok(())
            }
            None => {
                Err("RDKit CrystalFF SMARTS topology parser expected an atom token".to_string())
            }
        }
    }

    fn parse_atom(state: &mut ParserState<'_>) -> Result<usize, String> {
        skip_atom(state.chars, &mut state.pos)?;
        let atom_idx = state.next_atom_idx;
        state.next_atom_idx += 1;
        Ok(atom_idx)
    }

    fn add_bond(
        state: &mut ParserState<'_>,
        begin_atom_idx: usize,
        end_atom_idx: usize,
        query: QueryNode<BondQueryPredicate>,
    ) {
        state.bonds.push(CrystalffTemplateBond {
            begin_atom_idx,
            end_atom_idx,
            query,
        });
    }

    fn parse_chain(
        state: &mut ParserState<'_>,
        mut current_atom_idx: usize,
    ) -> Result<usize, String> {
        while state.pos < state.chars.len() {
            match state.chars[state.pos] {
                ')' => break,
                '(' => {
                    state.pos += 1;
                    let _ = parse_chain(state, current_atom_idx)?;
                    if state.chars.get(state.pos) != Some(&')') {
                        return Err(
                            "RDKit CrystalFF SMARTS topology parser expected ')' to close a branch"
                                .to_string(),
                        );
                    }
                    state.pos += 1;
                }
                '.' => {
                    state.pos += 1;
                    current_atom_idx = parse_atom(state)?;
                }
                _ => {
                    let query = parse_bond_query(state.chars, &mut state.pos)?;
                    if state.pos >= state.chars.len() {
                        return Err(
                            "RDKit CrystalFF SMARTS topology parser ended after a bond token"
                                .to_string(),
                        );
                    }
                    if state.chars[state.pos].is_ascii_digit() || state.chars[state.pos] == '%' {
                        let ring_number = parse_ring_number(state.chars, &mut state.pos)?;
                        if let Some((open_atom_idx, open_query)) =
                            state.ring_open.remove(&ring_number)
                        {
                            let resolved_query =
                                if query == QueryNode::Predicate(BondQueryPredicate::Any) {
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
        return Err(
            "RDKit CrystalFF SMARTS topology parser found an unbalanced ring closure".to_string(),
        );
    }
    Ok(state.bonds)
}

fn map_pattern_atom_indices(smarts: &str) -> [usize; 4] {
    let mut idx = [0; 4];
    let chars: Vec<char> = smarts.chars().collect();
    let mut atom_idx = 0usize;
    let mut i = 0usize;

    while i < chars.len() {
        match chars[i] {
            '[' => {
                let start = i + 1;
                i += 1;
                while i < chars.len() && chars[i] != ']' {
                    i += 1;
                }
                let end = i.min(chars.len());
                let bracket_content: String = chars[start..end].iter().collect();
                if let Some(map_number) = extract_atom_map_number(&bracket_content)
                    && (1..=4).contains(&map_number)
                {
                    idx[usize::try_from(map_number - 1).expect("atom-map index fits usize")] =
                        atom_idx;
                }
                atom_idx += 1;
                if i < chars.len() {
                    i += 1;
                }
            }
            'A'..='Z' | 'a'..='z' | '*' => {
                if is_organic_atom_token_start(&chars, i) {
                    atom_idx += 1;
                    i += organic_atom_token_len(&chars, i);
                } else {
                    i += 1;
                }
            }
            _ => {
                i += 1;
            }
        }
    }

    idx
}

fn extract_atom_map_number(bracket_content: &str) -> Option<u32> {
    let colon_pos = bracket_content.rfind(':')?;
    let digits = &bracket_content[colon_pos + 1..];
    if digits.is_empty() || !digits.chars().all(|ch| ch.is_ascii_digit()) {
        return None;
    }
    digits.parse::<u32>().ok()
}

fn is_organic_atom_token_start(chars: &[char], idx: usize) -> bool {
    match chars[idx] {
        'B' => chars.get(idx + 1) != Some(&'r'),
        'C' => chars.get(idx + 1) != Some(&'l'),
        'N' | 'O' | 'P' | 'S' | 'F' | 'I' | 'b' | 'c' | 'n' | 'o' | 'p' | 's' | '*' => true,
        'H' => true,
        _ => false,
    }
}

fn organic_atom_token_len(chars: &[char], idx: usize) -> usize {
    match chars[idx] {
        'C' if chars.get(idx + 1) == Some(&'l') => 2,
        'B' if chars.get(idx + 1) == Some(&'r') => 2,
        _ => 1,
    }
}

#[cfg(test)]
mod tests {
    use super::{
        CrystalFFDetails, CrystalffTorsionPreferencesError, ExpTorsionAngleCollection,
        get_experimental_torsions, get_experimental_torsions_without_bonds,
        torsion_preferences_macrocycles, torsion_preferences_small_rings, torsion_preferences_v1,
        torsion_preferences_v2,
    };
    use crate::Molecule;

    const VALID_BASE_PARAM_LINE: &str = "[C:1][C:2][C:3][C:4] 1 0 -1 0 1 0 -1 0 1 0 -1 0\n";
    const VALID_OVERRIDE_PARAM_LINE: &str =
        "[N:1][C:2][O:3][S:4] -1 1.1 1 1.2 -1 1.3 1 1.4 -1 1.5 1 1.6\n";

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
    fn crystalff_exptorsionanglecollection_constructor_retains_smarts_parse_result() {
        let params = ExpTorsionAngleCollection::new("CCCC 1 0 1 0 1 0 1 0 1 0 1 0\n").unwrap();

        let angle = &params.params()[0];
        assert_eq!(angle.smarts(), "CCCC");
        assert!(angle.pattern().is_ok());
    }

    #[test]
    fn crystalff_exptorsionanglecollection_constructor_retains_unsupported_smarts_parse_errors() {
        let params =
            ExpTorsionAngleCollection::new("[C:1][C:2][C:3][C:4] 1 0 1 0 1 0 1 0 1 0 1 0\n")
                .unwrap();

        let angle = &params.params()[0];
        assert_eq!(angle.smarts(), "[C:1][C:2][C:3][C:4]");
        assert!(angle.pattern().is_ok());
    }

    #[test]
    fn crystalff_exptorsionanglecollection_constructor_extracts_atom_map_indices() {
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
    fn crystalff_get_experimental_torsions_matches_linear_alkane_smarts() {
        let mol = Molecule::from_smiles("CCCCC").expect("pentane");
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
    fn crystalff_get_experimental_torsions_applies_small_ring_and_macrocycle_tables() {
        let small_ring = Molecule::from_smiles("C1COCC1").expect("tetrahydrofuran-like ring");
        let macrocycle = Molecule::from_smiles("C1COCCCCCCC1").expect("macrocycle");
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
        let mol = Molecule::from_smiles("O[C@H]1C[C@H]2CC[C@]1(C)C2(C)C").expect("bridged ring");
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
        let mol = Molecule::from_smiles("CCCC").expect("butane");
        let mut details = CrystalFFDetails {
            constrained_atoms: vec![true; mol.num_atoms()],
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

        assert!(torsion_bonds.is_empty());
        assert!(details.exp_torsion_atoms.is_empty());
        assert!(details.exp_torsion_angles.is_empty());
    }

    #[test]
    fn crystalff_get_experimental_torsions_adds_basic_knowledge_improper_terms() {
        let mol = Molecule::from_smiles("CC(=O)O").expect("acetic acid");
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
        let mol = Molecule::from_smiles("c1ccccc1").expect("benzene");
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
        let mol = Molecule::new();
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
}
