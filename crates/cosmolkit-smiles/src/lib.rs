//! Detached SMILES parser and writer.
//!
//! The parser owns SMILES syntax and lowers directly into model blocks.  It
//! does not construct a live `Molecule`; the facade is responsible for
//! installing the returned values into its runtime state.

use std::collections::{BTreeMap, HashMap, HashSet};

use cosmolkit_cx::parse_cx_extensions_with_atom_window;
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondDirection, BondId, BondSpec, CoordinateBlock,
    MoleculeProperties, PropertyValue, TopologyBlock,
};
use cosmolkit_types::{BondOrder, ChiralTag, Element};

mod canonical_rank;
mod cx_lowering;
mod cx_writer;
mod finalize_stereo;
mod fragment;
mod stereo;
mod writer;

#[doc(hidden)]
pub use cx_writer::{CoordinateSource, select_cx_coordinates};
pub use cx_writer::{
    CxCoordinateSelection, CxSmilesFields, CxSmilesWriteParams, write_cx_smiles,
    write_cx_smiles_with_params,
};
pub use finalize_stereo::{SmilesStereoError, finalize_smiles_stereo};
pub use fragment::FragmentWriteInputError;
#[doc(hidden)]
pub use writer::{SmartsTraversalError, prepare_smarts_serialization_topology};

pub use writer::in_organic_subset;
#[doc(hidden)]
pub use writer::{AtomColor, DfsBondSymbols, MolStackElem, dfs_build_query_stack};
pub use writer::{
    RandomSmilesWriteParams, SmilesWriteOutput, SmilesWriteParams, write_fragment_cx_smiles,
    write_fragment_smiles_output, write_random_smiles_vector, write_smiles,
    write_smiles_with_params, write_smiles_with_random,
};

const CXSMILES_BOND_IDX_PROP: &str = "_cxsmilesBondIdx";

#[derive(Debug, Clone, PartialEq)]
pub struct SmilesRecord {
    pub topology: TopologyBlock,
    pub coordinates: CoordinateBlock,
    pub properties: MoleculeProperties,
}

/// Borrowed detached input for serialization. No live molecule or commit authority.
#[derive(Clone, Copy)]
pub struct SmilesRecordView<'a> {
    pub topology: &'a TopologyBlock,
    pub coordinates: &'a CoordinateBlock,
    pub properties: &'a MoleculeProperties,
    /// Existing source ring state; absence requests the source fallback.
    pub rings: Option<&'a cosmolkit_core::RingInfo>,
}

impl<'a> From<&'a SmilesRecord> for SmilesRecordView<'a> {
    fn from(record: &'a SmilesRecord) -> Self {
        Self {
            topology: &record.topology,
            coordinates: &record.coordinates,
            properties: &record.properties,
            rings: None,
        }
    }
}

impl<'a> From<&SmilesRecordView<'a>> for SmilesRecordView<'a> {
    fn from(record: &SmilesRecordView<'a>) -> Self {
        *record
    }
}

impl SmilesRecordView<'_> {
    /// Materialize only at the source algorithm's owned working-copy boundary.
    fn to_owned_record(self) -> SmilesRecord {
        SmilesRecord {
            topology: self.topology.clone(),
            coordinates: self.coordinates.clone(),
            properties: self.properties.clone(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum SmilesParseError {
    #[error("{0}")]
    StereoGroup(#[from] cosmolkit_model::StereoGroupError),

    #[error("detached writer topology edit failed: {0}")]
    TopologyEdit(#[source] cosmolkit_model::TopologyEditError),

    #[error("CX link-node atom {atom} is outside {atom_count} source atoms")]
    CxLinkAtomOutOfRange { atom: usize, atom_count: usize },
    #[error("CX link-node order index {index} is outside {count} entries")]
    CxLinkOrderOutOfRange { index: usize, count: usize },

    #[error("CX zero-bond source bond {bond} is outside {bond_count} bonds")]
    CxZeroBondOutOfRange { bond: BondId, bond_count: usize },

    #[error("CX typed-bond source bond {bond} is outside {bond_count} bonds")]
    CxTypedBondOutOfRange { bond: BondId, bond_count: usize },

    #[error("CX atom-property source atom {atom} is outside {atom_count} atoms")]
    CxAtomPropertyAtomOutOfRange { atom: AtomId, atom_count: usize },

    #[error("CX conformer atom {atom} is outside {atom_count} coordinate rows")]
    CxCoordinateAtomOutOfRange { atom: AtomId, atom_count: usize },

    #[error("source canonicalization invariant failed: {0}")]
    WriterCanonicalInvariant(&'static str),
    #[error("source signed property read failed: {0}")]
    WriterInt(#[source] cosmolkit_core::PropertyIntReadError),
    #[error("source size_t property read failed: {0}")]
    WriterULong(#[source] cosmolkit_core::PropertyULongReadError),
    #[error("source ring preparation failed: {0}")]
    WriterRings(#[source] cosmolkit_core::RingFindingError),
    #[error("source stereo perception failed: {0}")]
    WriterLegacyStereo(#[source] cosmolkit_core::LegacyStereoError),
    #[error("source potential stereo getter failed: {0}")]
    WriterPotentialStereo(#[source] cosmolkit_core::PotentialStereoError),
    #[error("source non-tetrahedral permutation failed: {0}")]
    WriterParserStereoOrder(#[source] cosmolkit_core::parser_stereo_order::ParserStereoOrderError),

    #[error("Too many rings open at once. SMILES cannot be generated.")]
    TraversalTooManyOpenRings,
    #[error("source traversal {state} index {index} is outside {count} entries")]
    TraversalStateIndex {
        state: &'static str,
        index: usize,
        count: usize,
    },

    #[error("source stereo order failed: {0}")]
    WriterStereoOrder(#[source] cosmolkit_core::StereoOrderError),
    #[error("CX property value has the wrong source type: {0}")]
    WriterPropertyKind(#[source] cosmolkit_model::PropertyValueError),
    #[error("invalid CX coordinate state: {0}")]
    Coordinates(#[source] cosmolkit_model::CoordinateValidationError),
    #[error("source unsigned property read failed: {0}")]
    WriterNumeric(#[source] cosmolkit_core::PropertyUIntReadError),
    #[error("CX required property read failed at atom {atom:?}: {source}")]
    WriterRequiredProperty {
        atom: AtomId,
        #[source]
        source: cosmolkit_core::RequiredPropertyStringError,
    },
    #[error("property listing failed on atom {atom}: {source}")]
    WriterPropertyList {
        atom: cosmolkit_model::AtomId,
        #[source]
        source: cosmolkit_model::AtomPropertyError,
    },
    #[error("atom property operation failed: {0}")]
    AtomProperty(#[from] cosmolkit_model::AtomPropertyError),
    #[error("bond property operation failed: {0}")]
    BondProperty(#[from] cosmolkit_model::BondValueError),
    #[error("molecule property operation failed: {0}")]
    MoleculeProperty(#[from] cosmolkit_model::MoleculePropertyError),
    #[error("property value kind mismatch: {0}")]
    PropertyKind(#[from] cosmolkit_model::PropertyValueError),
    #[error("parser carrier finalization failed: {0}")]
    ParserCarrier(#[from] cosmolkit_core::parser_helpers::ParserCarrierError),
    #[error("unsupported SMILES token '{token}' at byte {offset}")]
    Unsupported { token: char, offset: usize },
    #[error("invalid SMILES syntax at byte {offset}: {message}")]
    Syntax { offset: usize, message: String },
    #[error("invalid atom at byte {offset}: {text}")]
    Atom { offset: usize, text: String },
    #[error("unclosed ring index {index}")]
    UnclosedRing { index: u32 },
    #[error("invalid CX extension: {0}")]
    Cx(String),
    #[error("unsupported CX record: {0}")]
    UnsupportedCx(&'static str),
    #[error("unsupported SMILES writer branch: {0}")]
    UnsupportedWriter(&'static str),
    #[error("canonical SMILES ranking failed: {0}")]
    CanonicalRank(String),
    #[error("SMILES writer valence preparation failed: {0}")]
    WriterValence(String),
    #[error("SMILES writer kekulization failed: {0}")]
    WriterKekulize(#[source] cosmolkit_core::KekulizeError),
    #[error("SMILES writer property conversion failed: {0}")]
    WriterProperty(#[source] cosmolkit_core::PropertyStringError),
    #[error("SMILES writer stereochemistry preparation failed: {0}")]
    WriterStereo(String),
    #[error("SMILES sanitize stage failed: {0}")]
    ParserSanitize(#[source] cosmolkit_core::SanitizeError),
    #[error("SMILES hydrogen removal stage failed: {0}")]
    ParserRemoveHydrogens(#[source] cosmolkit_core::HydrogenError),
    #[error("SMILES final stereo stage failed: {0}")]
    ParserStereo(#[source] finalize_stereo::SmilesStereoError),
    #[error("SMILES source atropisomer stage failed: {0}")]
    ParserAtropisomer(#[source] cosmolkit_core::AtropisomerError),
    #[error("SMILES writer stereo-group inversion failed: {0}")]
    WriterStereoBond(#[source] cosmolkit_model::BondValueError),
    #[error("CX coordinate selection is ambiguous ({two_d_count} 2D and {three_d_count} 3D sets)")]
    AmbiguousCoordinateSelection {
        two_d_count: usize,
        three_d_count: usize,
    },
    #[error("requested CX coordinate set {selection:?} is unavailable")]
    MissingCoordinateSelection { selection: CxCoordinateSelection },
    #[error("root atom index {atom_index} is out of range for {atom_count} atoms")]
    WriterRootAtomOutOfRange {
        atom_index: usize,
        atom_count: usize,
    },
    #[error("invalid detached model: {0}")]
    Model(String),
    #[error("SMILES replacement key must not be empty")]
    EmptyReplacementKey,
    #[error("SMILES replacements do not converge; cyclic key '{key}' remains active")]
    ReplacementCycle { key: String },
}

impl From<cosmolkit_model::TopologyEditError> for SmilesParseError {
    fn from(error: cosmolkit_model::TopologyEditError) -> Self {
        // Rust-only typed transport: source clearComputedProps throws the
        // underlying atom/bond property error. The detached edit wrapper must
        // retain that same category and reason, rather than flatten to text.
        match error {
            cosmolkit_model::TopologyEditError::AtomProperty(source) => Self::AtomProperty(source),
            cosmolkit_model::TopologyEditError::InvalidBond(source) => Self::BondProperty(source),
            cosmolkit_model::TopologyEditError::MoleculeProperty(source) => {
                Self::MoleculeProperty(source)
            }
            other => Self::TopologyEdit(other),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmilesParseParams {
    pub sanitize: bool,
    pub allow_cxsmiles: bool,
    pub strict_cxsmiles: bool,
    pub parse_name: bool,
    pub remove_hs: bool,
    pub skip_cleanup: bool,
    pub debug_parse: bool,
    pub replacements: BTreeMap<String, String>,
}

impl Default for SmilesParseParams {
    fn default() -> Self {
        // BEGIN RDKIT CPP TYPE v2::SmilesParse::SmilesParserParams
        // RDKit✔️✔️:   bool sanitize = true;
        // RDKit✔️✔️:   bool allowCXSMILES = true;
        // RDKit✔️✔️:   bool strictCXSMILES = true;
        // RDKit✔️✔️:   bool parseName = true;
        // RDKit✔️✔️:   bool removeHs = true;
        // RDKit✔️✔️:   bool skipCleanup = false;
        // RDKit✔️✔️:   bool debugParse = false;
        // RDKit✔️✔️:   std::map<std::string, std::string> replacements;
        // END RDKIT CPP TYPE v2::SmilesParse::SmilesParserParams
        Self {
            sanitize: true,
            allow_cxsmiles: true,
            strict_cxsmiles: true,
            parse_name: true,
            remove_hs: true,
            skip_cleanup: false,
            debug_parse: false,
            replacements: BTreeMap::new(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PreprocessedSmiles {
    smiles: String,
    name: String,
    cx_part: String,
}

fn replacement_key_is_cyclic(
    start: &str,
    current: &str,
    replacements: &BTreeMap<String, String>,
    visiting: &mut HashSet<String>,
) -> bool {
    if !visiting.insert(current.to_owned()) {
        return current == start;
    }
    let cyclic = replacements.get(current).is_some_and(|replacement| {
        replacements.keys().any(|next| {
            replacement.contains(next)
                && (next == start || replacement_key_is_cyclic(start, next, replacements, visiting))
        })
    });
    visiting.remove(current);
    cyclic
}

fn apply_replacements(
    input: &str,
    replacements: &BTreeMap<String, String>,
) -> Result<String, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION preprocessSmiles (replacement loop)
    // RDKit✔️✔️:   if (!params.replacements.empty()) {
    // RDKit✔️✔️:     std::string smi = lsmiles;
    // RDKit✔️✔️:     for (auto loopAgain = true; loopAgain;) {
    // RDKit✔️✔️:       loopAgain = false;
    // RDKit✔️✔️:       for (const auto &pr : params.replacements) {
    // RDKit✔️✔️:         if (smi.find(pr.first) != std::string::npos) {
    // RDKit✔️✔️:           loopAgain = true;
    // RDKit✔️✔️:           boost::replace_all(smi, pr.first, pr.second);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     lsmiles = smi;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION preprocessSmiles (replacement loop)
    // BTreeMap preserves std::map's key ordering. RDKit assigns termination
    // responsibility to callers; the Rust boundary rejects an active cyclic
    // rewrite explicitly instead of hanging forever. Each convergent pass is
    // otherwise the same ordered replace-all loop and remains O(p * r * n).
    if replacements.keys().any(String::is_empty) {
        return Err(SmilesParseError::EmptyReplacementKey);
    }
    if replacements.is_empty() {
        return Ok(input.to_owned());
    }
    let cyclic_keys = replacements
        .keys()
        .filter(|key| replacement_key_is_cyclic(key, key, replacements, &mut HashSet::new()))
        .cloned()
        .collect::<Vec<_>>();
    let mut smiles = input.to_owned();
    let mut seen = HashSet::new();
    loop {
        if !seen.insert(smiles.clone()) {
            let key = replacements
                .keys()
                .find(|key| smiles.contains(key.as_str()))
                .cloned()
                .unwrap_or_default();
            return Err(SmilesParseError::ReplacementCycle { key });
        }
        let mut changed = false;
        for (from, to) in replacements {
            if smiles.contains(from) {
                smiles = smiles.replace(from, to);
                changed = true;
            }
        }
        if !changed {
            return Ok(smiles);
        }
        if let Some(key) = cyclic_keys.iter().find(|key| smiles.contains(key.as_str())) {
            return Err(SmilesParseError::ReplacementCycle { key: key.clone() });
        }
    }
}

fn preprocess_smiles(
    input: &str,
    params: &SmilesParseParams,
) -> Result<PreprocessedSmiles, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION preprocessSmiles (name/CX partition)
    // RDKit✔️✔️:   cxPart = "";
    // RDKit✔️✔️:   name = "";
    // RDKit✔️✔️:   if (params.parseName && !params.allowCXSMILES) {
    // RDKit✔️✔️:     size_t sidx = smiles.find_first_of(" \t");
    // RDKit✔️✔️:     if (sidx != std::string::npos && sidx != 0) {
    // RDKit✔️✔️:       lsmiles = smiles.substr(0, sidx);
    // RDKit✔️✔️:       name = boost::trim_copy(smiles.substr(sidx, smiles.size() - sidx));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (params.allowCXSMILES) {
    // RDKit✔️✔️:     size_t sidx = smiles.find_first_of(" \t");
    // RDKit✔️✔️:     if (sidx != std::string::npos && sidx != 0) {
    // RDKit✔️✔️:       lsmiles = smiles.substr(0, sidx);
    // RDKit✔️✔️:       cxPart = boost::trim_copy(smiles.substr(sidx, smiles.size() - sidx));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (lsmiles.empty()) {
    // RDKit✔️✔️:     lsmiles = smiles;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION preprocessSmiles (name/CX partition)
    let mut smiles = String::new();
    let mut name = String::new();
    let mut cx_part = String::new();
    if params.parse_name && !params.allow_cxsmiles {
        if let Some(split) = input.find([' ', '\t'])
            && split != 0
        {
            smiles = input[..split].to_owned();
            name = input[split..].trim().to_owned();
        }
    } else if params.allow_cxsmiles
        && let Some(split) = input.find([' ', '\t'])
        && split != 0
    {
        smiles = input[..split].to_owned();
        cx_part = input[split..].trim().to_owned();
    }
    if smiles.is_empty() {
        smiles = input.to_owned();
    }
    smiles = apply_replacements(&smiles, &params.replacements)?;
    Ok(PreprocessedSmiles {
        smiles,
        name,
        cx_part,
    })
}

fn atom_error(text: &str, offset: usize) -> SmilesParseError {
    SmilesParseError::Atom {
        offset,
        text: text.into(),
    }
}

fn parse_source_number(
    text: &str,
    cursor: &mut usize,
    offset: usize,
) -> Result<u32, SmilesParseError> {
    // RDKit✔️✔️: nonzero_number:  NONZERO_DIGIT_TOKEN
    // RDKit✔️✔️: | nonzero_number digit {
    // RDKit✔️✔️:   if($1 >= std::numeric_limits<std::int32_t>::max()/10 ||
    // RDKit✔️✔️:      $1*10 >= std::numeric_limits<std::int32_t>::max()-$2 ){
    // RDKit✔️✔️:      yyerror(input,molList,branchPoints,scanner,start_token,current_token_position,"number too large");
    // RDKit✔️✔️:      yyErrorCleanup(molList);
    // RDKit✔️✔️:      YYABORT;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   $$ = $1*10 + $2;
    // RDKit✔️✔️:   }
    // Source rejects the entire max/10 prefix, even when the next digit could
    // fit signed32. One bounded multiply/check per digit preserves its cost.
    let bytes = text.as_bytes();
    let Some(first) = bytes.get(*cursor).copied().filter(u8::is_ascii_digit) else {
        return Err(atom_error(text, offset));
    };
    *cursor += 1;
    if first == b'0' {
        return Ok(0);
    }
    let mut value = u32::from(first - b'0');
    while let Some(digit) = bytes.get(*cursor).copied().filter(u8::is_ascii_digit) {
        let digit = u32::from(digit - b'0');
        let max = i32::MAX as u32;
        if value >= max / 10 || value * 10 >= max - digit {
            return Err(atom_error(text, offset));
        }
        value = value * 10 + digit;
        *cursor += 1;
    }
    Ok(value)
}

fn parse_bracket_element(
    text: &str,
    cursor: &mut usize,
    offset: usize,
) -> Result<(Element, bool), SmilesParseError> {
    if text.as_bytes().get(*cursor) == Some(&b'*') {
        *cursor += 1;
        return Ok((Element::DUMMY, false));
    }
    if text.as_bytes().get(*cursor) == Some(&b'#') {
        *cursor += 1;
        let atomic_number = parse_source_number(text, cursor, offset)?;
        let atomic_number = u8::try_from(atomic_number).map_err(|_| atom_error(text, offset))?;
        return Element::from_atomic_number(atomic_number)
            .map(|element| (element, false))
            .ok_or_else(|| atom_error(text, offset));
    }

    for width in [3, 2, 1] {
        let Some(token) = text.get(*cursor..*cursor + width) else {
            continue;
        };
        let aromatic = matches!(
            token,
            "b" | "c" | "n" | "o" | "p" | "s" | "as" | "se" | "te"
        );
        let element = match token {
            "b" => Some(Element::B),
            "c" => Some(Element::C),
            "n" => Some(Element::N),
            "o" => Some(Element::O),
            "p" => Some(Element::P),
            "s" => Some(Element::S),
            "as" => Some(Element::AS),
            "se" => Some(Element::SE),
            // RDKit✔️✔️: <IN_ATOM_STATE>te   {	yylval->atom = new Atom( 52 );
            // RDKit✔️✔️: 			yylval->atom->setIsAromatic(true);
            // RDKit✔️✔️: 				return AROMATIC_ATOM_TOKEN;
            // RDKit✔️✔️: 			}
            // Behavior: this token is legal only inside brackets, matching
            // IN_ATOM_STATE; parse_simple_atom keeps the bare token invalid.
            // Complexity: one fixed token arm, no allocation or graph scan.
            "te" => Some(Element::TE),
            _ => Element::from_symbol(token),
        };
        if let Some(element) = element {
            *cursor += width;
            return Ok((element, aromatic));
        }
    }
    Err(atom_error(text, offset))
}

fn normalize_chiral_spec(
    mut spec: AtomSpec,
    text: &str,
    offset: usize,
) -> Result<AtomSpec, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION CheckChiralitySpecifications
    // RDKit✔️✔️: static const std::map<int, int> permutationLimits = {
    // RDKit✔️✔️:     {RDKit::Atom::ChiralType::CHI_TETRAHEDRAL, 2},
    // RDKit✔️✔️:     {RDKit::Atom::ChiralType::CHI_ALLENE, 2},
    // RDKit✔️✔️:     {RDKit::Atom::ChiralType::CHI_SQUAREPLANAR, 3},
    // RDKit✔️✔️:     {RDKit::Atom::ChiralType::CHI_OCTAHEDRAL, 30},
    // RDKit✔️✔️:     {RDKit::Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL, 20}};
    // RDKit✔️✔️: if (!checkChiralPermutation(atom->getChiralTag(), permutation)) {
    // RDKit✔️✔️:   throw SmilesParseException(error);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (atom->getChiralTag() == RDKit::Atom::ChiralType::CHI_TETRAHEDRAL) {
    // RDKit✔️✔️:   if (permutation == 0 || permutation == 1) {
    // RDKit✔️✔️:     atom->setChiralTag(RDKit::Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    // RDKit✔️✔️:     atom->clearProp(common_properties::_chiralPermutation);
    // RDKit✔️✔️:   } else if (permutation == 2) {
    // RDKit✔️✔️:     atom->setChiralTag(RDKit::Atom::ChiralType::CHI_TETRAHEDRAL_CW);
    // RDKit✔️✔️:     atom->clearProp(common_properties::_chiralPermutation);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION CheckChiralitySpecifications
    let Some(permutation) = spec.chiral_permutation() else {
        return Ok(spec);
    };
    let limit = match spec.chiral_tag() {
        ChiralTag::Tetrahedral | ChiralTag::Allene => 2,
        ChiralTag::SquarePlanar => 3,
        ChiralTag::TrigonalBipyramidal => 20,
        ChiralTag::Octahedral => 30,
        _ => return Ok(spec),
    };
    if permutation > limit {
        return Err(atom_error(text, offset));
    }
    if spec.chiral_tag() == ChiralTag::Tetrahedral {
        spec = spec
            .with_chiral_tag(if permutation < 2 {
                ChiralTag::TetrahedralCcw
            } else {
                ChiralTag::TetrahedralCw
            })
            .without_chiral_permutation();
    }
    Ok(spec)
}

fn parse_atom(text: &str, offset: usize) -> Result<AtomSpec, SmilesParseError> {
    // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy bracket atomd/charge_element/h_element/chiral_element
    // RDKit✔️✔️: | ATOM_OPEN_TOKEN charge_element COLON_TOKEN number ATOM_CLOSE_TOKEN
    // RDKit✔️✔️: {
    // RDKit✔️✔️:   $$ = $2;
    // RDKit✔️✔️:   $$->setNoImplicit(true);
    // RDKit✔️✔️:   $$->setProp(RDKit::common_properties::molAtomMapNumber,$4);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: | ATOM_OPEN_TOKEN charge_element ATOM_CLOSE_TOKEN
    // RDKit✔️✔️: {
    // RDKit✔️✔️:   $$ = $2;
    // RDKit✔️✔️:   $2->setNoImplicit(true);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: charge_element: h_element
    // RDKit✔️✔️: | h_element PLUS_TOKEN { $1->setFormalCharge(1); }
    // RDKit✔️✔️: | h_element PLUS_TOKEN PLUS_TOKEN { $1->setFormalCharge(2); }
    // RDKit✔️✔️: | h_element PLUS_TOKEN number { $1->setFormalCharge($3); }
    // RDKit✔️✔️: | h_element MINUS_TOKEN { $1->setFormalCharge(-1); }
    // RDKit✔️✔️: | h_element MINUS_TOKEN MINUS_TOKEN { $1->setFormalCharge(-2); }
    // RDKit✔️✔️: | h_element MINUS_TOKEN number { $1->setFormalCharge(-$3); }
    // RDKit✔️✔️: chiral_element: element
    // RDKit✔️✔️: | element AT_TOKEN { $1->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW); }
    // RDKit✔️✔️: | element AT_TOKEN AT_TOKEN { $1->setChiralTag(Atom::CHI_TETRAHEDRAL_CW); }
    // END RDKIT CPP GRAMMAR ACTION smiles.yy bracket atomd/charge_element/h_element/chiral_element
    // Both parsers consume the bracket payload once from left to right. The
    // Rust cursor uses no token allocation and has O(n) time/O(1) extra space,
    // matching the generated lexer/parser complexity for this modeled grammar.
    let mut cursor = 0;
    let isotope = if text.as_bytes().first().is_some_and(u8::is_ascii_digit) {
        Some(parse_source_number(text, &mut cursor, offset)?)
    } else {
        None
    };
    let (element, aromatic) = parse_bracket_element(text, &mut cursor, offset)?;
    let mut spec = AtomSpec::new(element)
        .with_aromatic(aromatic)
        .with_no_implicit(true);
    if let Some(isotope) = isotope {
        spec = spec.with_isotope(u16::try_from(isotope).map_err(|_| atom_error(text, offset))?);
    }

    if text.as_bytes().get(cursor) == Some(&b'@') {
        cursor += 1;
        if text.as_bytes().get(cursor) == Some(&b'@') {
            cursor += 1;
            spec = spec.with_chiral_tag(ChiralTag::TetrahedralCw);
        } else {
            let spaces = text[cursor..]
                .bytes()
                .take_while(|byte| *byte == b' ')
                .count();
            let class_start = cursor + spaces;
            let class = [
                ("TH", ChiralTag::Tetrahedral),
                ("AL", ChiralTag::Allene),
                ("SP", ChiralTag::SquarePlanar),
                ("TB", ChiralTag::TrigonalBipyramidal),
                ("OH", ChiralTag::Octahedral),
            ]
            .into_iter()
            .find(|(name, _)| text[class_start..].starts_with(name));
            if let Some((name, tag)) = class {
                cursor = class_start + name.len();
                let permutation = if text.as_bytes().get(cursor).is_some_and(u8::is_ascii_digit) {
                    parse_source_number(text, &mut cursor, offset)?
                } else {
                    0
                };
                if permutation == 0
                    && text
                        .as_bytes()
                        .get(cursor.saturating_sub(1))
                        .is_some_and(u8::is_ascii_digit)
                {
                    return Err(atom_error(text, offset));
                }
                spec = spec
                    .with_chiral_tag(tag)
                    .with_chiral_permutation(permutation);
            } else {
                spec = spec.with_chiral_tag(ChiralTag::TetrahedralCcw);
            }
        }
    }

    if text.as_bytes().get(cursor) == Some(&b'H') {
        cursor += 1;
        let count = if text.as_bytes().get(cursor).is_some_and(u8::is_ascii_digit) {
            parse_source_number(text, &mut cursor, offset)?
        } else {
            1
        };
        spec = spec
            .with_explicit_hydrogens(u8::try_from(count).map_err(|_| atom_error(text, offset))?);
    }

    match text.as_bytes().get(cursor) {
        Some(b'+') => {
            cursor += 1;
            let charge = if text.as_bytes().get(cursor) == Some(&b'+') {
                cursor += 1;
                2
            } else if text.as_bytes().get(cursor).is_some_and(u8::is_ascii_digit) {
                parse_source_number(text, &mut cursor, offset)?
            } else {
                1
            };
            spec = spec
                .with_formal_charge(i8::try_from(charge).map_err(|_| atom_error(text, offset))?);
        }
        Some(b'-') => {
            cursor += 1;
            let magnitude = if text.as_bytes().get(cursor) == Some(&b'-') {
                cursor += 1;
                2
            } else if text.as_bytes().get(cursor).is_some_and(u8::is_ascii_digit) {
                parse_source_number(text, &mut cursor, offset)?
            } else {
                1
            };
            let magnitude = i16::try_from(magnitude).map_err(|_| atom_error(text, offset))?;
            spec = spec.with_formal_charge(
                i8::try_from(-magnitude).map_err(|_| atom_error(text, offset))?,
            );
        }
        _ => {}
    }

    if text.as_bytes().get(cursor) == Some(&b':') {
        cursor += 1;
        spec = spec.with_atom_map(parse_source_number(text, &mut cursor, offset)?);
    }
    if cursor != text.len() {
        return Err(atom_error(text, offset));
    }
    normalize_chiral_spec(spec, text, offset)
}

fn parse_simple_atom(input: &str, offset: usize) -> Result<(AtomSpec, usize), SmilesParseError> {
    // BEGIN RDKIT CPP LEXER RULES smiles.ll simple atoms
    // RDKit✔️✔️: B  { yylval->atom = new Atom(5);return ORGANIC_ATOM_TOKEN; }
    // RDKit✔️✔️: C  { yylval->atom = new Atom(6);return ORGANIC_ATOM_TOKEN; }
    // RDKit✔️✔️: N  { yylval->atom = new Atom(7);return ORGANIC_ATOM_TOKEN; }
    // RDKit✔️✔️: O  { yylval->atom = new Atom(8);return ORGANIC_ATOM_TOKEN; }
    // RDKit✔️✔️: P  { yylval->atom = new Atom(15);return ORGANIC_ATOM_TOKEN; }
    // RDKit✔️✔️: S  { yylval->atom = new Atom(16);return ORGANIC_ATOM_TOKEN; }
    // RDKit✔️✔️: F  { yylval->atom = new Atom(9);return ORGANIC_ATOM_TOKEN; }
    // RDKit✔️✔️: Cl { yylval->atom = new Atom(17);return ORGANIC_ATOM_TOKEN; }
    // RDKit✔️✔️: Br { yylval->atom = new Atom(35);return ORGANIC_ATOM_TOKEN; }
    // RDKit✔️✔️: I  { yylval->atom = new Atom(53);return ORGANIC_ATOM_TOKEN; }
    // RDKit✔️✔️: b		    {	yylval->atom = new Atom ( 5 );
    // RDKit✔️✔️: 			yylval->atom->setIsAromatic(true);
    // RDKit✔️✔️: 				return AROMATIC_ATOM_TOKEN;
    // RDKit✔️✔️: 			}
    // RDKit✔️✔️: c		    {	yylval->atom = new Atom ( 6 );
    // RDKit✔️✔️: 			yylval->atom->setIsAromatic(true);
    // RDKit✔️✔️: 				return AROMATIC_ATOM_TOKEN;
    // RDKit✔️✔️: 			}
    // RDKit✔️✔️: n		    {	yylval->atom = new Atom( 7 );
    // RDKit✔️✔️: 			yylval->atom->setIsAromatic(true);
    // RDKit✔️✔️: 				return AROMATIC_ATOM_TOKEN;
    // RDKit✔️✔️: 			}
    // RDKit✔️✔️: o		    {	yylval->atom = new Atom( 8 );
    // RDKit✔️✔️: 			yylval->atom->setIsAromatic(true);
    // RDKit✔️✔️: 				return AROMATIC_ATOM_TOKEN;
    // RDKit✔️✔️: 			}
    // RDKit✔️✔️: p		    {	yylval->atom = new Atom( 15 );
    // RDKit✔️✔️: 			yylval->atom->setIsAromatic(true);
    // RDKit✔️✔️: 				return AROMATIC_ATOM_TOKEN;
    // RDKit✔️✔️: 			}
    // RDKit✔️✔️: s		    {	yylval->atom = new Atom( 16 );
    // RDKit✔️✔️: 			yylval->atom->setIsAromatic(true);
    // RDKit✔️✔️: 				return AROMATIC_ATOM_TOKEN;
    // RDKit✔️✔️: 			}
    // END RDKIT CPP LEXER RULES smiles.ll simple atoms
    let (element, aromatic, consumed) = if input.starts_with("Cl") {
        (Element::CL, false, 2)
    } else if input.starts_with("Br") {
        (Element::BR, false, 2)
    } else {
        match input.as_bytes().first().copied() {
            Some(b'B') => (Element::B, false, 1),
            Some(b'C') => (Element::C, false, 1),
            Some(b'N') => (Element::N, false, 1),
            Some(b'O') => (Element::O, false, 1),
            Some(b'P') => (Element::P, false, 1),
            Some(b'S') => (Element::S, false, 1),
            Some(b'F') => (Element::F, false, 1),
            Some(b'I') => (Element::I, false, 1),
            Some(b'b') => (Element::B, true, 1),
            Some(b'c') => (Element::C, true, 1),
            Some(b'n') => (Element::N, true, 1),
            Some(b'o') => (Element::O, true, 1),
            Some(b'p') => (Element::P, true, 1),
            Some(b's') => (Element::S, true, 1),
            Some(b'*') => (Element::DUMMY, false, 1),
            _ => return Err(atom_error(input, offset)),
        }
    };
    Ok((AtomSpec::new(element).with_aromatic(aromatic), consumed))
}

fn get_unspecified_bond_type(atom1_aromatic: bool, atom2_aromatic: bool) -> BondOrder {
    // BEGIN RDKIT CPP FUNCTION GetUnspecifiedBondType
    // RDKit✔️✔️: Bond::BondType GetUnspecifiedBondType(const RWMol *mol, const Atom *atom1,
    // RDKit✔️✔️:                                       const Atom *atom2) {
    // RDKit✔️✔️:   PRECONDITION(mol, "no molecule");
    // RDKit✔️✔️:   PRECONDITION(atom1, "no atom1");
    // RDKit✔️✔️:   PRECONDITION(atom2, "no atom2");
    // RDKit✔️✔️:   Bond::BondType res;
    // RDKit✔️✔️:   if (atom1->getIsAromatic() && atom2->getIsAromatic()) {
    // RDKit✔️✔️:     res = Bond::AROMATIC;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res = Bond::SINGLE;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION GetUnspecifiedBondType
    if atom1_aromatic && atom2_aromatic {
        BondOrder::Aromatic
    } else {
        BondOrder::Single
    }
}

fn resolved_bond_spec(
    pending: Option<BondOrder>,
    atoms: &[Atom],
    begin: AtomId,
    end: AtomId,
    direction: BondDirection,
    query: Option<cosmolkit_model::QueryNode<cosmolkit_model::BondQueryPredicate>>,
) -> BondSpec {
    // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy explicit dative bond orientation
    // RDKit✔️✔️:   if( $2->getBondType() == Bond::DATIVER ){
    // RDKit✔️✔️:     $2->setBeginAtomIdx(atomIdx1);
    // RDKit✔️✔️:     $2->setEndAtomIdx(atomIdx2);
    // RDKit✔️✔️:     $2->setBondType(Bond::DATIVE);
    // RDKit✔️✔️:   }else if ( $2->getBondType() == Bond::DATIVEL ){
    // RDKit✔️✔️:     $2->setBeginAtomIdx(atomIdx2);
    // RDKit✔️✔️:     $2->setEndAtomIdx(atomIdx1);
    // RDKit✔️✔️:     $2->setBondType(Bond::DATIVE);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     $2->setBeginAtomIdx(atomIdx1);
    // RDKit✔️✔️:     $2->setEndAtomIdx(atomIdx2);
    // RDKit✔️✔️:   }
    // END RDKIT CPP GRAMMAR ACTION smiles.yy explicit dative bond orientation
    let order = pending.unwrap_or_else(|| {
        get_unspecified_bond_type(
            atoms[begin.index()].is_aromatic(),
            atoms[end.index()].is_aromatic(),
        )
    });
    let spec = match order {
        BondOrder::DativeLeft => BondSpec::new(end, begin, BondOrder::Dative),
        BondOrder::DativeRight => BondSpec::new(begin, end, BondOrder::Dative),
        order => BondSpec::new(begin, end, order).with_aromatic(order == BondOrder::Aromatic),
    };
    let spec = spec.with_direction(direction);
    match query {
        Some(query) => spec.with_query(query),
        None => spec,
    }
}

fn bond_order(symbol: char) -> Result<BondOrder, SmilesParseError> {
    match symbol {
        '-' => Ok(BondOrder::Single),
        '=' => Ok(BondOrder::Double),
        '#' => Ok(BondOrder::Triple),
        '$' => Ok(BondOrder::Quadruple),
        ':' => Ok(BondOrder::Aromatic),
        _ => Err(SmilesParseError::Unsupported {
            token: symbol,
            offset: 0,
        }),
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RingPartial {
    atom: AtomId,
    order: Option<BondOrder>,
    query: Option<cosmolkit_model::QueryNode<cosmolkit_model::BondQueryPredicate>>,
    direction: BondDirection,
    offset: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PendingRingClosure {
    ring: u32,
    opening: RingPartial,
    closing: RingPartial,
    cx_bond_index: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RingClosureRecord {
    ring: u32,
    bond: Option<BondId>,
}

fn push_smiles_bond(bonds: &mut Vec<Bond>, spec: BondSpec, cx_bond_index: Option<u32>) {
    // BEGIN COMPLETE RDKit .6 CHEM26 ring_implicit
    // RDKit❗✔️: | mol ring_number {
    // RDKit❗✔️:   RWMol * mp = (*molList)[$$];
    // RDKit❗✔️:   Atom *atom=mp->getActiveAtom();
    // RDKit❗✔️:   mp->setAtomBookmark(atom,$2);
    // RDKit❗✔️:
    // RDKit❗✔️:   Bond *newB = mp->createPartialBond(atom->getIdx(),
    // RDKit❗✔️: 				     Bond::UNSPECIFIED);
    // RDKit❗✔️:   mp->setBondBookmark(newB,$2);
    // RDKit❗✔️:   newB->setProp(RDKit::common_properties::_unspecifiedOrder,1);
    // RDKit❗✔️:   if(!(mp->getAllBondsWithBookmark($2).size()%2)){
    // RDKit❗✔️:     newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);
    // RDKit❗✔️:
    // RDKit❗✔️:   INT_VECT tmp;
    // RDKit❗✔️:   atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
    // RDKit❗✔️:   tmp.push_back(-($2+1));
    // RDKit❗✔️:   atom->setProp(RDKit::common_properties::_RingClosures,tmp);
    // RDKit❗✔️: }
    // END COMPLETE RDKit .6 CHEM26 ring_implicit
    // BEGIN COMPLETE RDKit .6 CHEM26 ring_explicit
    // RDKit❗✔️: | mol BOND_TOKEN ring_number {
    // RDKit❗✔️:   RWMol * mp = (*molList)[$$];
    // RDKit❗✔️:   Atom *atom=mp->getActiveAtom();
    // RDKit❗✔️:   Bond *newB = mp->createPartialBond(atom->getIdx(),
    // RDKit❗✔️: 				     $2->getBondType());
    // RDKit❗✔️:   if($2->hasProp(RDKit::common_properties::_unspecifiedOrder)){
    // RDKit❗✔️:     newB->setProp(RDKit::common_properties::_unspecifiedOrder,1);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   newB->setBondDir($2->getBondDir());
    // RDKit❗✔️:   mp->setAtomBookmark(atom,$3);
    // RDKit❗✔️:   mp->setBondBookmark(newB,$3);
    // RDKit❗✔️:   if(!(mp->getAllBondsWithBookmark($3).size()%2)){
    // RDKit❗✔️:     newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);
    // RDKit❗✔️:
    // RDKit❗✔️:   INT_VECT tmp;
    // RDKit❗✔️:   atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
    // RDKit❗✔️:   tmp.push_back(-($3+1));
    // RDKit❗✔️:   atom->setProp(RDKit::common_properties::_RingClosures,tmp);
    // RDKit❗✔️:   delete $2;
    // RDKit❗✔️: }
    // END COMPLETE RDKit .6 CHEM26 ring_explicit
    // BEGIN COMPLETE RDKit .6 CHEM26 ring_single
    // RDKit❗✔️: | mol MINUS_TOKEN ring_number {
    // RDKit❗✔️:   RWMol * mp = (*molList)[$$];
    // RDKit❗✔️:   Atom *atom=mp->getActiveAtom();
    // RDKit❗✔️:   Bond *newB = mp->createPartialBond(atom->getIdx(),
    // RDKit❗✔️: 				     Bond::SINGLE);
    // RDKit❗✔️:   mp->setAtomBookmark(atom,$3);
    // RDKit❗✔️:   mp->setBondBookmark(newB,$3);
    // RDKit❗✔️:   if(!(mp->getAllBondsWithBookmark($3).size()%2)){
    // RDKit❗✔️:     newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);
    // RDKit❗✔️:
    // RDKit❗✔️:   INT_VECT tmp;
    // RDKit❗✔️:   atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
    // RDKit❗✔️:   tmp.push_back(-($3+1));
    // RDKit❗✔️:   atom->setProp(RDKit::common_properties::_RingClosures,tmp);
    // RDKit❗✔️: }
    // END COMPLETE RDKit .6 CHEM26 ring_single
    // One private storage helper is shared by the native grammar reductions.
    // Ordinary productions reserve the logical counter but pass no property;
    // delayed ring closures alone pass their previously reserved parse slot.
    // Existing detached bond construction and validation remain unchanged.
    let spec = if let Some(source_index) = cx_bond_index {
        spec.with_prop(CXSMILES_BOND_IDX_PROP, PropertyValue::UInt(source_index))
            .expect("the internal CXSMILES bond-index property key is non-empty")
    } else {
        spec
    };
    let index = bonds.len();
    bonds.push(Bond::from_spec(BondId::new(index), spec));
}

fn opposite_bond_direction(direction: BondDirection) -> BondDirection {
    match direction {
        BondDirection::EndDownRight => BondDirection::EndUpRight,
        BondDirection::EndUpRight => BondDirection::EndDownRight,
        other => other,
    }
}

fn merged_ring_direction(target: &RingPartial, source: &RingPartial) -> BondDirection {
    // BEGIN RDKIT CPP FUNCTION swapBondDirIfNeeded
    // RDKit✔️✔️: void swapBondDirIfNeeded(Bond *bond1, const Bond *bond2) {
    // RDKit✔️✔️:   if (bond1->getBondDir() == Bond::NONE && bond2->getBondDir() != Bond::NONE) {
    // RDKit✔️✔️:     bond1->setBondDir(bond2->getBondDir());
    // RDKit✔️✔️:     if (bond1->getBeginAtom() != bond2->getBeginAtom()) {
    // RDKit✔️✔️:       switch (bond1->getBondDir()) {
    // RDKit✔️✔️:         case Bond::ENDDOWNRIGHT:
    // RDKit✔️✔️:           bond1->setBondDir(Bond::ENDUPRIGHT);
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case Bond::ENDUPRIGHT:
    // RDKit✔️✔️:           bond1->setBondDir(Bond::ENDDOWNRIGHT);
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION swapBondDirIfNeeded
    if target.direction != BondDirection::None || source.direction == BondDirection::None {
        return target.direction;
    }
    if target.atom == source.atom {
        source.direction
    } else {
        opposite_bond_direction(source.direction)
    }
}

fn parse_ring_number(text: &str, cursor: &mut usize) -> Result<u32, SmilesParseError> {
    // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy ring_number
    // RDKit✔️✔️: ring_number:  digit
    // RDKit✔️✔️: | PERCENT_TOKEN NONZERO_DIGIT_TOKEN digit { $$ = $2*10+$3; }
    // RDKit✔️✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit GROUP_CLOSE_TOKEN { $$ = $3; }
    // RDKit✔️✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit GROUP_CLOSE_TOKEN { $$ = $3*10+$4; }
    // RDKit✔️✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit GROUP_CLOSE_TOKEN { $$ = $3*100+$4*10+$5; }
    // RDKit✔️✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit GROUP_CLOSE_TOKEN { $$ = $3*1000+$4*100+$5*10+$6; }
    // RDKit✔️✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit digit GROUP_CLOSE_TOKEN { $$ = $3*10000+$4*1000+$5*100+$6*10+$7; }
    // END RDKIT CPP GRAMMAR ACTION smiles.yy ring_number
    let offset = *cursor;
    let bytes = text.as_bytes();
    match bytes.get(*cursor).copied() {
        Some(digit @ b'0'..=b'9') => {
            *cursor += 1;
            Ok(u32::from(digit - b'0'))
        }
        Some(b'%') => {
            *cursor += 1;
            if bytes.get(*cursor) == Some(&b'(') {
                *cursor += 1;
                let digit_start = *cursor;
                while bytes.get(*cursor).is_some_and(u8::is_ascii_digit) {
                    *cursor += 1;
                }
                let digits = &text[digit_start..*cursor];
                if digits.is_empty() || digits.len() > 5 || bytes.get(*cursor) != Some(&b')') {
                    return Err(SmilesParseError::Syntax {
                        offset,
                        message: "ring index in `%(...)` requires one to five digits".into(),
                    });
                }
                *cursor += 1;
                digits.parse().map_err(|_| SmilesParseError::Syntax {
                    offset,
                    message: "ring index is too large".into(),
                })
            } else {
                let Some(tens @ b'1'..=b'9') = bytes.get(*cursor).copied() else {
                    return Err(SmilesParseError::Syntax {
                        offset,
                        message: "two-digit ring index after `%` must start with 1-9".into(),
                    });
                };
                let Some(ones @ b'0'..=b'9') = bytes.get(*cursor + 1).copied() else {
                    return Err(SmilesParseError::Syntax {
                        offset,
                        message: "ring index requires two digits after `%`".into(),
                    });
                };
                *cursor += 2;
                Ok(u32::from(tens - b'0') * 10 + u32::from(ones - b'0'))
            }
        }
        _ => Err(SmilesParseError::Syntax {
            offset,
            message: "expected ring index".into(),
        }),
    }
}

fn check_ring_closure_branch_status(atoms: &mut [Atom], degrees: &[usize], atom: AtomId) {
    // BEGIN RDKIT CPP FUNCTION CheckRingClosureBranchStatus
    // RDKit✔️✔️: void CheckRingClosureBranchStatus(RDKit::Atom *atom, RDKit::RWMol *mp) {
    // RDKit✔️✔️:   // github #786 and #1652: if the ring closure comes after a branch,
    // RDKit✔️✔️:   // the stereochem is wrong.
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
    // RDKit✔️✔️:   PRECONDITION(mp, "bad mol");
    // RDKit✔️✔️:   if (atom->getIdx() != mp->getNumAtoms(true) - 1 &&
    // RDKit✔️✔️:       (atom->getDegree() == 1 ||
    // RDKit✔️✔️:        (atom->getDegree() == 2 && atom->getIdx() != 0) ||
    // RDKit✔️✔️:        (atom->getDegree() == 3 && atom->getIdx() == 0)) &&
    // RDKit✔️✔️:       (atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit✔️✔️:        atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW)) {
    // RDKit✔️✔️:     atom->invertChirality();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION CheckRingClosureBranchStatus
    let degree = degrees[atom.index()];
    let should_invert = atom.index() != atoms.len().saturating_sub(1)
        && (degree == 1
            || (degree == 2 && atom.index() != 0)
            || (degree == 3 && atom.index() == 0));
    if should_invert {
        let inverted = match atoms[atom.index()].chiral_tag() {
            ChiralTag::TetrahedralCw => Some(ChiralTag::TetrahedralCcw),
            ChiralTag::TetrahedralCcw => Some(ChiralTag::TetrahedralCw),
            _ => None,
        };
        if let Some(inverted) = inverted {
            atoms[atom.index()].set_chiral_tag(inverted);
        }
    }
}

fn close_ring_closures(
    atoms: &[Atom],
    bonds: &mut Vec<Bond>,
    closures: &mut [PendingRingClosure],
    ring_closures_by_atom: &mut [Vec<RingClosureRecord>],
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION CloseMolRings (pair validation and bond selection)
    // RDKit✔️✔️:       while (atomIt != atomsEnd) {
    // RDKit✔️✔️:         Atom *atom1 = *atomIt;
    // RDKit✔️✔️:         ++atomIt;
    // RDKit✔️✔️:         if (!toleratePartials && atomIt == atomsEnd) {
    // RDKit✔️✔️:           ReportParseError("unclosed ring");
    // RDKit✔️✔️:         } else if (atomIt != atomsEnd && *atomIt == atom1) {
    // RDKit✔️✔️:           // make sure we don't try to connect an atom to itself
    // RDKit✔️✔️:           // this was github #1925
    // RDKit✔️✔️:           auto fmt =
    // RDKit✔️✔️:               boost::format{
    // RDKit✔️✔️:                   "duplicated ring closure %1% bonds atom %2% to itself"} %
    // RDKit✔️✔️:               bookmark.first % atom1->getIdx();
    // RDKit✔️✔️:           std::string msg = fmt.str();
    // RDKit✔️✔️:           ReportParseError(msg.c_str(), true);
    // RDKit✔️✔️:         } else if (mol->getBondBetweenAtoms(atom1->getIdx(),
    // RDKit✔️✔️:                                             (*atomIt)->getIdx()) != nullptr) {
    // RDKit✔️✔️:           auto fmt =
    // RDKit✔️✔️:               boost::format{
    // RDKit✔️✔️:                   "ring closure %1% duplicates bond between atom %2% and atom "
    // RDKit✔️✔️:                   "%3%"} %
    // RDKit✔️✔️:               bookmark.first % atom1->getIdx() % (*atomIt)->getIdx();
    // RDKit✔️✔️:           std::string msg = fmt.str();
    // RDKit✔️✔️:           ReportParseError(msg.c_str(), true);
    // RDKit✔️✔️:         } else if (atomIt != atomsEnd) {
    // RDKit✔️✔️:           // we actually found an atom, so connect it to the first
    // RDKit✔️✔️:           Atom *atom2 = *atomIt;
    // RDKit✔️✔️:           ++atomIt;
    // RDKit✔️✔️:
    // RDKit✔️✔️:           // We're guaranteed two partial bonds, one for each time
    // RDKit✔️✔️:           // the ring index was used.  We give the first specification
    // RDKit✔️✔️:           // priority.
    // RDKit✔️✔️:           if (!bond1->hasProp(common_properties::_unspecifiedOrder)) {
    // RDKit✔️✔️:             matchedBond = bond1;
    // RDKit✔️✔️:             matchedBond->setEndAtomIdx(atom2->getIdx());
    // RDKit✔️✔️:           } else {
    // RDKit✔️✔️:             matchedBond = bond2;
    // RDKit✔️✔️:             matchedBond->setEndAtomIdx(atom1->getIdx());
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           if (matchedBond->getBondType() == Bond::UNSPECIFIED &&
    // RDKit✔️✔️:               !matchedBond->hasQuery()) {
    // RDKit✔️✔️:             Bond::BondType bondT = GetUnspecifiedBondType(mol, atom1, atom2);
    // RDKit✔️✔️:             matchedBond->setBondType(bondT);
    // RDKit✔️✔️:           }
    // END RDKIT CPP FUNCTION CloseMolRings (pair validation and bond selection)
    closures.sort_by_key(|closure| closure.ring);
    let mut bond_pairs = bonds
        .iter()
        .map(|bond| {
            let endpoints = (bond.begin().index(), bond.end().index());
            if endpoints.0 <= endpoints.1 {
                endpoints
            } else {
                (endpoints.1, endpoints.0)
            }
        })
        .collect::<std::collections::HashSet<_>>();
    for closure in closures {
        let atom1 = closure.opening.atom;
        let atom2 = closure.closing.atom;
        if atom1 == atom2 {
            return Err(SmilesParseError::Syntax {
                offset: closure.closing.offset,
                message: format!(
                    "duplicated ring closure {} bonds atom {} to itself",
                    closure.ring,
                    atom1.index()
                ),
            });
        }
        let pair = if atom1.index() <= atom2.index() {
            (atom1.index(), atom2.index())
        } else {
            (atom2.index(), atom1.index())
        };
        if !bond_pairs.insert(pair) {
            return Err(SmilesParseError::Syntax {
                offset: closure.closing.offset,
                message: format!(
                    "ring closure {} duplicates bond between atom {} and atom {}",
                    closure.ring,
                    atom1.index(),
                    atom2.index()
                ),
            });
        }

        let (selected, other) = if closure.opening.order.is_some() {
            (&mut closure.opening, &mut closure.closing)
        } else {
            (&mut closure.closing, &mut closure.opening)
        };
        let order = selected.order.unwrap_or_else(|| {
            get_unspecified_bond_type(
                atoms[atom1.index()].is_aromatic(),
                atoms[atom2.index()].is_aromatic(),
            )
        });
        // RDKit✔️✔️:             if (matchedBond->getBondType() == Bond::DATIVEL) {
        // RDKit✔️✔️:               matchedBond->setBeginAtomIdx(atom2->getIdx());
        // RDKit✔️✔️:               matchedBond->setEndAtomIdx(atom1->getIdx());
        // RDKit✔️✔️:               matchedBond->setBondType(Bond::DATIVE);
        // RDKit✔️✔️:             } else if (matchedBond->getBondType() == Bond::DATIVER) {
        // RDKit✔️✔️:               matchedBond->setEndAtomIdx(atom2->getIdx());
        // RDKit✔️✔️:               matchedBond->setBondType(Bond::DATIVE);
        // RDKit✔️✔️:             } else {
        // RDKit✔️✔️:               matchedBond->setEndAtomIdx(atom2->getIdx());
        // RDKit✔️✔️:             }
        // The same source block appears for the retained closing partial with
        // atom1/atom2 exchanged. `selected` and `other` encode that exchange.
        let (begin, end, order) = match order {
            BondOrder::DativeLeft => (other.atom, selected.atom, BondOrder::Dative),
            BondOrder::DativeRight => (selected.atom, other.atom, BondOrder::Dative),
            order => (selected.atom, other.atom, order),
        };
        let direction = merged_ring_direction(&selected, &other);
        let bond = BondId::new(bonds.len());
        let spec = BondSpec::new(begin, end, order)
            .with_aromatic(order == BondOrder::Aromatic)
            .with_direction(direction);
        // The source retains the selected partial Bond object itself, including
        // its query. Transfer that identity rather than infer it from order.
        let spec = match selected.query.take() {
            Some(query) => spec.with_query(query),
            None => spec,
        };
        // Drop the unused source partial query as closeMolRings deletes it.
        let _ = other.query.take();
        push_smiles_bond(bonds, spec, Some(closure.cx_bond_index));
        // RDKit✔️✔️:             *closurePos = bondIdx - 1;
        // Each grammar occurrence left a ring-number placeholder at the atom.
        // Replace the first still-unresolved occurrence exactly where it was
        // parsed so GetBondOrdering retains SMILES token order.
        for atom in [atom1, atom2] {
            let record = ring_closures_by_atom[atom.index()]
                .iter_mut()
                .find(|record| record.ring == closure.ring && record.bond.is_none())
                .ok_or_else(|| {
                    SmilesParseError::Model(format!(
                        "ring closure {} is missing its atom-order placeholder",
                        closure.ring
                    ))
                })?;
            record.bond = Some(bond);
        }
    }
    Ok(())
}

fn adjust_atom_chirality_flags(
    atoms: &mut [Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    ring_closures_by_atom: &[Vec<RingClosureRecord>],
    smiles_start_atoms: &[bool],
) -> Result<(), SmilesParseError> {
    let rings = ring_closures_by_atom
        .iter()
        .map(|records| {
            records
                .iter()
                .map(|record| {
                    record.bond.ok_or_else(|| {
                        SmilesParseError::Model(format!(
                            "ring closure {} remained unresolved during chirality adjustment",
                            record.ring
                        ))
                    })
                })
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()?;
    let assignments = cosmolkit_core::parser_helpers::parser_chirality_assignments(
        atoms,
        bonds,
        |index| {
            adjacency
                .neighbors_of(index)
                .iter()
                .map(|neighbor| (neighbor.atom_index, neighbor.bond))
        },
        &rings,
        smiles_start_atoms,
    )
    .map_err(|error| match error {
        cosmolkit_core::parser_helpers::ParserCarrierError::PropertyString(source) => {
            SmilesParseError::WriterProperty(source)
        }
        cosmolkit_core::parser_helpers::ParserCarrierError::Model(message) => {
            SmilesParseError::Model(message)
        }
        cosmolkit_core::parser_helpers::ParserCarrierError::AtomProperty(source) => {
            SmilesParseError::AtomProperty(source)
        }
    })?;
    for (atom, (tag, permutation)) in atoms.iter_mut().zip(assignments) {
        atom.set_chiral_tag(tag);
        atom.set_chiral_permutation(permutation);
    }
    Ok(())
}

/// Complete delayed component cleanup after reaction-wide CX annotations.
#[doc(hidden)]
pub fn cleanup_after_parsing(record: &mut SmilesRecord) -> Result<(), SmilesParseError> {
    cleanup_after_parsing_impl(record, true)
}

fn cleanup_after_parsing_impl(
    record: &mut SmilesRecord,
    cleanup_nontetrahedral: bool,
) -> Result<(), SmilesParseError> {
    // RDKit✔️❌: void CleanupAfterParsing(RWMol *mol) {
    // RDKit✔️❌:   PRECONDITION(mol, "no molecule");
    // RDKit✔️❌:   for (auto atom : mol->atoms()) {
    // RDKit✔️❌:     atom->clearProp(common_properties::_RingClosures);
    // RDKit✔️❌:     atom->clearProp(common_properties::_SmilesStart);
    // RDKit✔️❌:     std::string label;
    // RDKit✔️❌:     if (atom->getAtomicNum() == 0 &&
    // RDKit✔️❌:         atom->getPropIfPresent(common_properties::atomLabel, label)) {
    // RDKit✔️❌:       // marvinsketch can output higher labels than _AP1 and _AP2, but they
    // RDKit✔️❌:       // aren't part of the MOL file spec so we don't treat them as attachment
    // RDKit✔️❌:       // points
    // RDKit✔️❌:       if (label == "_AP1") {
    // RDKit✔️❌:         atom->setProp(common_properties::_fromAttachPoint, 1);
    // RDKit✔️❌:       } else if (label == "_AP2") {
    // RDKit✔️❌:         atom->setProp(common_properties::_fromAttachPoint, 2);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   for (auto bond : mol->bonds()) {
    // RDKit✔️❌:     bond->clearProp(common_properties::_unspecifiedOrder);
    // RDKit✔️❌:     bond->clearProp("_cxsmilesBondIdx");
    // RDKit✔️❌:   }
    // RDKit✔️❌:   for (auto sg : RDKit::getSubstanceGroups(*mol)) {
    // RDKit✔️❌:     sg.clearProp("_cxsmilesindex");
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (!Chirality::getAllowNontetrahedralChirality()) {
    // RDKit✔️❌:     bool needWarn = false;
    // RDKit✔️❌:     for (auto atom : mol->atoms()) {
    // RDKit✔️❌:       if (atom->hasProp(common_properties::_chiralPermutation)) {
    // RDKit✔️❌:         needWarn = true;
    // RDKit✔️❌:         atom->clearProp(common_properties::_chiralPermutation);
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (atom->getChiralTag() > Atom::ChiralType::CHI_OTHER) {
    // RDKit✔️❌:         needWarn = true;
    // RDKit✔️❌:         atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (needWarn) {
    // RDKit✔️❌:       BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:           << "ignoring non-tetrahedral stereo specification since setAllowNontetrahedralChirality() is false."
    // RDKit✔️❌:           << std::endl;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // Destination composition preserves all reached error/effect ordering.
    // Atom algorithms remain in CORE; copying each SGroup matches native
    // by-value iteration, with the known canonical tree/order allocation cost.

    cosmolkit_core::parser_helpers::cleanup_parser_atoms(&mut record.topology.atoms)?;
    for bond in &mut record.topology.bonds {
        bond.clear_prop("_unspecifiedOrder")?;
        bond.clear_prop(CXSMILES_BOND_IDX_PROP)?;
    }
    cosmolkit_core::parser_helpers::cleanup_parser_substance_groups(
        &record.topology.substance_groups,
    )?;
    if cleanup_nontetrahedral {
        cosmolkit_core::parser_helpers::cleanup_parser_nontetrahedral_atoms(
            &mut record.topology.atoms,
        )?;
    }
    Ok(())
}

/// Parse SMILES into detached topology, coordinate, and property blocks.
pub fn parse_smiles(
    input: &str,
    params: &SmilesParseParams,
) -> Result<SmilesRecord, SmilesParseError> {
    // MAIN's detached parsing seam: the live facade owns the subsequent
    // sanitize/removeHs/finalize sequence and installs its derived carriers.
    parse_smiles_stages(input, params, false)
}

/// Complete source parser composition for detached reaction components.
/// Unlike the live-facade parsing seam, this includes requested final stages.
#[doc(hidden)]
pub fn parse_smiles_complete_source(
    input: &str,
    params: &SmilesParseParams,
) -> Result<SmilesRecord, SmilesParseError> {
    parse_smiles_stages(input, params, true)
}

fn parse_smiles_stages(
    input: &str,
    params: &SmilesParseParams,
    complete: bool,
) -> Result<SmilesRecord, SmilesParseError> {
    // BEGIN COMPLETE RDKit .6 CHEM26 ordinary_implicit
    // RDKit❗✔️: | mol atomd       {
    // RDKit❗✔️:   RWMol *mp = (*molList)[$$];
    // RDKit❗✔️:   Atom *a1 = mp->getActiveAtom();
    // RDKit❗✔️:   int atomIdx1=a1->getIdx();
    // RDKit❗✔️:   int atomIdx2=mp->addAtom($2,true,true);
    // RDKit❗✔️:   mp->addBond(atomIdx1,atomIdx2,
    // RDKit❗✔️: 	      SmilesParseOps::GetUnspecifiedBondType(mp,a1,mp->getAtomWithIdx(atomIdx2)));
    // RDKit❗✔️:   ++numBondsParsed;
    // RDKit❗✔️:   //delete $2;
    // RDKit❗✔️: }
    // END COMPLETE RDKit .6 CHEM26 ordinary_implicit
    // BEGIN COMPLETE RDKit .6 CHEM26 ordinary_explicit
    // RDKit❗✔️: | mol BOND_TOKEN atomd  {
    // RDKit❗✔️:   RWMol *mp = (*molList)[$$];
    // RDKit❗✔️:   int atomIdx1 = mp->getActiveAtom()->getIdx();
    // RDKit❗✔️:   int atomIdx2 = mp->addAtom($3,true,true);
    // RDKit❗✔️:   if( $2->getBondType() == Bond::DATIVER ){
    // RDKit❗✔️:     $2->setBeginAtomIdx(atomIdx1);
    // RDKit❗✔️:     $2->setEndAtomIdx(atomIdx2);
    // RDKit❗✔️:     $2->setBondType(Bond::DATIVE);
    // RDKit❗✔️:   }else if ( $2->getBondType() == Bond::DATIVEL ){
    // RDKit❗✔️:     $2->setBeginAtomIdx(atomIdx2);
    // RDKit❗✔️:     $2->setEndAtomIdx(atomIdx1);
    // RDKit❗✔️:     $2->setBondType(Bond::DATIVE);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     $2->setBeginAtomIdx(atomIdx1);
    // RDKit❗✔️:     $2->setEndAtomIdx(atomIdx2);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   ++numBondsParsed;
    // RDKit❗✔️:   mp->addBond($2,true);
    // RDKit❗✔️:   //delete $3;
    // RDKit❗✔️: }
    // END COMPLETE RDKit .6 CHEM26 ordinary_explicit
    // BEGIN COMPLETE RDKit .6 CHEM26 ordinary_single
    // RDKit❗✔️: | mol MINUS_TOKEN atomd {
    // RDKit❗✔️:   RWMol *mp = (*molList)[$$];
    // RDKit❗✔️:   int atomIdx1 = mp->getActiveAtom()->getIdx();
    // RDKit❗✔️:   int atomIdx2 = mp->addAtom($3,true,true);
    // RDKit❗✔️:   mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
    // RDKit❗✔️:   ++numBondsParsed;
    // RDKit❗✔️:   //delete $3;
    // RDKit❗✔️: }
    // END COMPLETE RDKit .6 CHEM26 ordinary_single
    // BEGIN COMPLETE RDKit .6 CHEM26 branch_implicit
    // RDKit❗✔️: | mol branch_open_token atomd {
    // RDKit❗✔️:   RWMol *mp = (*molList)[$$];
    // RDKit❗✔️:   Atom *a1 = mp->getActiveAtom();
    // RDKit❗✔️:   int atomIdx1=a1->getIdx();
    // RDKit❗✔️:   int atomIdx2=mp->addAtom($3,true,true);
    // RDKit❗✔️:   mp->addBond(atomIdx1,atomIdx2,
    // RDKit❗✔️: 	      SmilesParseOps::GetUnspecifiedBondType(mp,a1,mp->getAtomWithIdx(atomIdx2)));
    // RDKit❗✔️:   ++numBondsParsed;
    // RDKit❗✔️:   branchPoints.push_back({atomIdx1, $2});
    // RDKit❗✔️: }
    // END COMPLETE RDKit .6 CHEM26 branch_implicit
    // BEGIN COMPLETE RDKit .6 CHEM26 branch_explicit
    // RDKit❗✔️: | mol branch_open_token BOND_TOKEN atomd  {
    // RDKit❗✔️:   RWMol *mp = (*molList)[$$];
    // RDKit❗✔️:   int atomIdx1 = mp->getActiveAtom()->getIdx();
    // RDKit❗✔️:   int atomIdx2 = mp->addAtom($4,true,true);
    // RDKit❗✔️:   if( $3->getBondType() == Bond::DATIVER ){
    // RDKit❗✔️:     $3->setBeginAtomIdx(atomIdx1);
    // RDKit❗✔️:     $3->setEndAtomIdx(atomIdx2);
    // RDKit❗✔️:     $3->setBondType(Bond::DATIVE);
    // RDKit❗✔️:   }else if ( $3->getBondType() == Bond::DATIVEL ){
    // RDKit❗✔️:     $3->setBeginAtomIdx(atomIdx2);
    // RDKit❗✔️:     $3->setEndAtomIdx(atomIdx1);
    // RDKit❗✔️:     $3->setBondType(Bond::DATIVE);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     $3->setBeginAtomIdx(atomIdx1);
    // RDKit❗✔️:     $3->setEndAtomIdx(atomIdx2);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   ++numBondsParsed;
    // RDKit❗✔️:   mp->addBond($3,true);
    // RDKit❗✔️:
    // RDKit❗✔️:   branchPoints.push_back({atomIdx1, $2});
    // RDKit❗✔️: }
    // END COMPLETE RDKit .6 CHEM26 branch_explicit
    // BEGIN COMPLETE RDKit .6 CHEM26 branch_single
    // RDKit❗✔️: | mol branch_open_token MINUS_TOKEN atomd {
    // RDKit❗✔️:   RWMol *mp = (*molList)[$$];
    // RDKit❗✔️:   int atomIdx1 = mp->getActiveAtom()->getIdx();
    // RDKit❗✔️:   int atomIdx2 = mp->addAtom($4,true,true);
    // RDKit❗✔️:   mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
    // RDKit❗✔️:   ++numBondsParsed;
    // RDKit❗✔️:   branchPoints.push_back({atomIdx1, $2});
    // RDKit❗✔️: }
    // END COMPLETE RDKit .6 CHEM26 branch_single
    // These six source reductions coalesce into the two native atom-token
    // branches below. Both reserve exactly one unsigned parse-order slot before
    // storing each ordinary edge; ring reservations remain in the ring branch.

    // BEGIN COMPLETE PINNED SF253 MolFromSmiles
    // RDKit❗❌: std::unique_ptr<RWMol> MolFromSmiles(const std::string &smiles,
    // RDKit❗❌:                                      const SmilesParserParams &params) {
    // RDKit❗❌:   // Calling MolFromSmiles in a multithreaded context is generally safe *unless*
    // RDKit❗❌:   // the value of debugParse is different for different threads. The if
    // RDKit❗❌:   // statement below avoids a TSAN warning in the case where multiple threads
    // RDKit❗❌:   // all use the same value for debugParse.
    // RDKit❗❌:   if (yysmiles_debug != params.debugParse) {
    // RDKit❗❌:     yysmiles_debug = params.debugParse;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::string lsmiles, name, cxPart;
    // RDKit❗❌:   preprocessSmiles(smiles, params, lsmiles, name, cxPart);
    // RDKit❗❌:   // strip any leading/trailing whitespace:
    // RDKit❗❌:   // boost::trim_if(smi,boost::is_any_of(" \t\r\n"));
    // RDKit❗❌:   auto res = toMol(lsmiles, smiles_parse, lsmiles);
    // RDKit❗❌:   if (!res) {
    // RDKit❗❌:     return res;
    // RDKit❗❌:   }
    // RDKit❗❌:   handleCXPartAndName(res.get(), params, cxPart, name);
    // RDKit❗❌:
    // RDKit❗❌:   // get a conformer
    // RDKit❗❌:   const Conformer *conf = nullptr, *conf3d = nullptr;
    // RDKit❗❌:   if (res && res->getNumConformers() > 0) {
    // RDKit❗❌:     for (unsigned int confId = 0; confId < res->getNumConformers(); ++confId) {
    // RDKit❗❌:       auto *testConf = &res->getConformer(confId);
    // RDKit❗❌:       if (!testConf->is3D()) {
    // RDKit❗❌:         if (conf == nullptr) {  // only take the first 2d conf
    // RDKit❗❌:           conf = testConf;
    // RDKit❗❌:         }
    // RDKit❗❌:       } else {
    // RDKit❗❌:         if (conf3d == nullptr) {  // only take the first 3d conf
    // RDKit❗❌:           conf3d = testConf;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       if (conf != nullptr && conf3d != nullptr) {
    // RDKit❗❌:         break;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (res->hasProp(SmilesParseOps::detail::_needsDetectAtomStereo)) {
    // RDKit❗❌:     // we encountered a wedged bond in the CXSMILES,
    // RDKit❗❌:     // these need to be handled the same way they were in mol files
    // RDKit❗❌:     res->clearProp(SmilesParseOps::detail::_needsDetectAtomStereo);
    // RDKit❗❌:
    // RDKit❗❌:     if (conf) {
    // RDKit❗❌:       MolOps::assignChiralTypesFromBondDirs(*res, conf->getId());
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // if we read a 3D conformer, set the stereo:
    // RDKit❗❌:   // if (res->getNumConformers() && res->getConformer().is3D()) {
    // RDKit❗❌:   if (!conf && conf3d) {
    // RDKit❗❌:     res->updatePropertyCache(false);
    // RDKit❗❌:     MolOps::assignChiralTypesFrom3D(*res, conf3d->getId(), true);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (conf) {
    // RDKit❗❌:     Atropisomers::detectAtropisomerChirality(*res, conf);
    // RDKit❗❌:   } else if (conf3d) {
    // RDKit❗❌:     Atropisomers::detectAtropisomerChirality(*res, conf3d);
    // RDKit❗❌:   } else {
    // RDKit❗❌:     Atropisomers::detectAtropisomerChirality(*res, nullptr);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (res && (params.sanitize || params.removeHs)) {
    // RDKit❗❌:     if (params.removeHs) {
    // RDKit❗❌:       MolOps::RemoveHsParameters rhp;
    // RDKit❗❌:       rhp.updateExplicitCount = true;
    // RDKit❗❌:       MolOps::removeHs(*res, rhp, params.sanitize);
    // RDKit❗❌:     } else if (params.sanitize) {
    // RDKit❗❌:       MolOps::sanitizeMol(*res);
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     if (res->hasProp(SmilesParseOps::detail::_needsDetectBondStereo)) {
    // RDKit❗❌:       // we encountered either wiggly bond in the CXSMILES,
    // RDKit❗❌:       // these need to be handled the same way they were in mol files
    // RDKit❗❌:       if (conf || conf3d) {
    // RDKit❗❌:         MolOps::clearSingleBondDirFlags(*res);
    // RDKit❗❌:       }
    // RDKit❗❌:       MolOps::setDoubleBondNeighborDirections(*res, conf ? conf : conf3d);
    // RDKit❗❌:     }
    // RDKit❗❌:     res->clearProp(SmilesParseOps::detail::_needsDetectBondStereo);
    // RDKit❗❌:     // figure out stereochemistry:
    // RDKit❗❌:     bool cleanIt = true, force = true, flagPossible = true;
    // RDKit❗❌:     MolOps::assignStereochemistry(*res, cleanIt, force, flagPossible);
    // RDKit❗❌:   } else {
    // RDKit❗❌:     //  we still need to do something about double bond stereochemistry
    // RDKit❗❌:     //  (was github issue 337)
    // RDKit❗❌:     //  now that atom stereochem has been perceived, the wedging
    // RDKit❗❌:     //  information is no longer needed, so we clear
    // RDKit❗❌:     //  single bond dir flags:
    // RDKit❗❌:     MolOps::clearSingleBondDirFlags(*res, true);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (res && res->hasProp(common_properties::_NeedsQueryScan)) {
    // RDKit❗❌:     res->clearProp(common_properties::_NeedsQueryScan);
    // RDKit❗❌:     if (!params.sanitize) {
    // RDKit❗❌:       // we know that this can be the ring bond query, do ring perception if we
    // RDKit❗❌:       // need to:
    // RDKit❗❌:       MolOps::fastFindRings(*res);
    // RDKit❗❌:     }
    // RDKit❗❌:     QueryOps::completeMolQueries(res.get(), 0xDEADBEEF);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (res) {
    // RDKit❗❌:     if (!params.skipCleanup) {
    // RDKit❗❌:       SmilesParseOps::CleanupAfterParsing(res.get());
    // RDKit❗❌:     }
    // RDKit❗❌:     if (!name.empty()) {
    // RDKit❗❌:       res->setProp(common_properties::_Name, name);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // END COMPLETE PINNED SF253 MolFromSmiles
    // Canonical CORE owners execute the source sanitize/H-removal branches.
    // Atrop candidate ordering remains user-authorized deferred equivalence.
    // Concrete SMILES query-only capabilities remain explicitly unsupported.
    // Local cost: sanitizer detached topology copy and geometry XY lift are
    // known extra allocation; no duplicated chemistry algorithm is introduced.
    // BEGIN RDKIT CPP FUNCTION MolFromSmiles (debug selection)
    // RDKit✔️✔️:   if (yysmiles_debug != params.debugParse) {
    // RDKit✔️✔️:     yysmiles_debug = params.debugParse;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION MolFromSmiles (debug selection)
    // This hand-written parser has no process-global generated-parser trace
    // switch. The option is therefore chemistry-neutral and thread-safe.
    let _debug_parse = params.debug_parse;
    let preprocessed = preprocess_smiles(input, params)?;
    let graph_text = preprocessed.smiles.trim();
    let bytes = graph_text.as_bytes();
    let mut atoms = Vec::<Atom>::new();
    let mut degrees = Vec::<usize>::new();
    let mut bonds = Vec::<Bond>::new();
    let mut rings = HashMap::<u32, RingPartial>::new();
    let mut ring_closures = Vec::<PendingRingClosure>::new();
    let mut ring_closures_by_atom = Vec::<Vec<RingClosureRecord>>::new();
    let mut smiles_start_atoms = Vec::<bool>::new();
    let mut branches = Vec::<AtomId>::new();
    let mut branch_needs_atom = false;
    let mut current = None::<AtomId>;
    let mut pending = None;
    let mut pending_query = None;
    let mut pending_direction = BondDirection::None;
    let mut next_cx_bond_index = 0_u32;
    let mut index = 0;
    while index < bytes.len() {
        // RDKit✔️✔️: mol: atomd {
        // RDKit✔️✔️: | mol BOND_TOKEN atomd  {
        // RDKit✔️✔️: | mol MINUS_TOKEN atomd {
        // RDKit✔️✔️: | mol SEPARATOR_TOKEN atomd {
        // RDKit✔️✔️: | mol BOND_TOKEN ring_number {
        // RDKit✔️✔️: | mol branch_open_token BOND_TOKEN atomd  {
        // RDKit✔️✔️: | mol branch_open_token MINUS_TOKEN atomd {
        // These smiles.yy productions consume a bond exactly once, followed
        // by an atom or (outside a new branch) a ring number. A separator and
        // a branch opening require an atom, not another separator or branch.
        // Constant-time grammar state; no rescans or token allocations.
        let token = bytes[index] as char;
        let atom_token = token.is_ascii_alphabetic() || matches!(token, '[' | '*');
        let ring_token = token.is_ascii_digit() || token == '%';
        let bond_token = matches!(token, '-' | '=' | '#' | ':' | '$' | '~' | '/' | '\\')
            || graph_text[index..].starts_with("<-");
        let has_pending_bond = pending.is_some() || pending_direction != BondDirection::None;
        let missing_operand = if has_pending_bond {
            !atom_token && !(ring_token && !branch_needs_atom)
        } else if current.is_none() {
            !atom_token
        } else if branch_needs_atom {
            !atom_token && !bond_token
        } else {
            false
        };
        if missing_operand {
            return Err(SmilesParseError::Syntax {
                offset: index,
                message: "expected atom or ring operand".into(),
            });
        }
        // BEGIN RDKIT CPP LEXER RULES smiles.ll dative bonds
        // RDKit✔️✔️: \-\>  { yylval->bond = new Bond(Bond::DATIVER);
        // RDKit✔️✔️:       return BOND_TOKEN; }
        // RDKit✔️✔️: \<\-  { yylval->bond = new Bond(Bond::DATIVEL);
        // RDKit✔️✔️:       return BOND_TOKEN; }
        // END RDKIT CPP LEXER RULES smiles.ll dative bonds
        if graph_text[index..].starts_with("->") {
            pending = Some(BondOrder::DativeRight);
            pending_query = None;
            index += 2;
            continue;
        }
        if graph_text[index..].starts_with("<-") {
            pending = Some(BondOrder::DativeLeft);
            pending_query = None;
            index += 2;
            continue;
        }
        match bytes[index] as char {
            '-' | '=' | '#' | ':' | '$' => {
                pending = Some(bond_order(bytes[index] as char)?);
                pending_query = None;
                index += 1;
            }
            '~' => {
                // RDKit✔️✔️: \~	{ yylval->bond = new QueryBond();
                // RDKit✔️✔️: 	  yylval->bond->setQuery(makeBondNullQuery());
                // RDKit✔️✔️: 	  return BOND_TOKEN;  }
                pending = Some(BondOrder::Unspecified);
                pending_query = Some(cosmolkit_model::QueryNode::predicate(
                    cosmolkit_model::BondQueryPredicate::Any,
                ));
                index += 1;
            }
            '/' => {
                pending_direction = BondDirection::EndUpRight;
                index += 1;
            }
            '\\' => {
                // RDKit✔️✔️: [\\]{1,2}    { yylval->bond = new Bond(Bond::UNSPECIFIED);
                // RDKit✔️✔️: 	yylval->bond->setProp(RDKit::common_properties::_unspecifiedOrder,1);
                // RDKit✔️✔️: 	yylval->bond->setBondDir(Bond::ENDDOWNRIGHT);
                // RDKit✔️✔️: 	return BOND_TOKEN;  }
                // Flex consumes one or two backslashes as ONE token, not two
                // consecutive bonds. A third backslash remains invalid.
                pending_direction = BondDirection::EndDownRight;
                index += if bytes.get(index + 1) == Some(&b'\\') {
                    2
                } else {
                    1
                };
            }
            '(' => {
                branches.push(current.ok_or_else(|| SmilesParseError::Syntax {
                    offset: index,
                    message: "branch has no preceding atom".into(),
                })?);
                branch_needs_atom = true;
                index += 1;
            }
            ')' => {
                current = Some(branches.pop().ok_or_else(|| SmilesParseError::Syntax {
                    offset: index,
                    message: "unmatched branch close".into(),
                })?);
                index += 1;
            }
            '.' => {
                current = None;
                index += 1;
            }
            '0'..='9' | '%' => {
                let ring_offset = index;
                let ring = parse_ring_number(graph_text, &mut index)?;
                let atom = current.ok_or_else(|| SmilesParseError::Syntax {
                    offset: ring_offset,
                    message: "ring index has no preceding atom".into(),
                })?;
                check_ring_closure_branch_status(&mut atoms, &degrees, atom);
                ring_closures_by_atom[atom.index()].push(RingClosureRecord { ring, bond: None });
                let partial = RingPartial {
                    atom,
                    order: pending,
                    query: pending_query.take(),
                    direction: pending_direction,
                    offset: ring_offset,
                };
                if let Some(opening) = rings.remove(&ring) {
                    ring_closures.push(PendingRingClosure {
                        ring,
                        opening,
                        closing: partial,
                        cx_bond_index: take_smiles_bond_source_index(&mut next_cx_bond_index),
                    });
                } else {
                    rings.insert(ring, partial);
                }
                pending = None;
                pending_direction = BondDirection::None;
            }
            '[' => {
                let start = index;
                let end = graph_text[index + 1..]
                    .find(']')
                    .map(|value| index + 1 + value)
                    .ok_or_else(|| SmilesParseError::Syntax {
                        offset: index,
                        message: "unclosed bracket atom".into(),
                    })?;
                let spec = parse_atom(&graph_text[index + 1..end], start)?;
                let atom = AtomId::new(atoms.len());
                smiles_start_atoms.push(current.is_none());
                atoms.push(Atom::from_spec(atom, spec));
                degrees.push(0);
                ring_closures_by_atom.push(Vec::new());
                if let Some(previous) = current {
                    let spec = resolved_bond_spec(
                        pending,
                        &atoms,
                        previous,
                        atom,
                        pending_direction,
                        pending_query.take(),
                    );
                    let _ = take_smiles_bond_source_index(&mut next_cx_bond_index);
                    push_smiles_bond(&mut bonds, spec, None);
                    degrees[previous.index()] += 1;
                    degrees[atom.index()] += 1;
                }
                current = Some(atom);
                branch_needs_atom = false;
                pending = None;
                pending_query = None;
                pending_direction = BondDirection::None;
                index = end + 1;
            }
            token if token.is_ascii_alphabetic() || token == '*' => {
                let start = index;
                let (spec, consumed) = parse_simple_atom(&graph_text[start..], start)?;
                index += consumed;
                let atom = AtomId::new(atoms.len());
                smiles_start_atoms.push(current.is_none());
                atoms.push(Atom::from_spec(atom, spec));
                degrees.push(0);
                ring_closures_by_atom.push(Vec::new());
                if let Some(previous) = current {
                    let spec = resolved_bond_spec(
                        pending,
                        &atoms,
                        previous,
                        atom,
                        pending_direction,
                        pending_query.take(),
                    );
                    let _ = take_smiles_bond_source_index(&mut next_cx_bond_index);
                    push_smiles_bond(&mut bonds, spec, None);
                    degrees[previous.index()] += 1;
                    degrees[atom.index()] += 1;
                }
                current = Some(atom);
                branch_needs_atom = false;
                pending = None;
                pending_query = None;
                pending_direction = BondDirection::None;
            }
            token => {
                return Err(SmilesParseError::Unsupported {
                    token,
                    offset: index,
                });
            }
        }
    }
    // RDKit✔️✔️: | meta_start error EOS_TOKEN{
    // RDKit✔️✔️:   yyerrok;
    // RDKit✔️✔️:   yyErrorCleanup(molList);
    // RDKit✔️✔️:   YYABORT;
    // A dangling bond or separator cannot reduce to the grammar's mol.
    if pending.is_some()
        || pending_direction != BondDirection::None
        || (!bytes.is_empty() && current.is_none())
    {
        return Err(SmilesParseError::Syntax {
            offset: graph_text.len(),
            message: "missing final atom or ring operand".into(),
        });
    }
    if let Some((&index, _)) = rings.iter().next() {
        return Err(SmilesParseError::UnclosedRing { index });
    }
    if !branches.is_empty() {
        return Err(SmilesParseError::Syntax {
            offset: graph_text.len(),
            message: "unclosed branch".into(),
        });
    }
    close_ring_closures(
        &atoms,
        &mut bonds,
        &mut ring_closures,
        &mut ring_closures_by_atom,
    )?;
    let adjacency = AdjacencyList::from_topology(atoms.len(), &bonds);
    adjust_atom_chirality_flags(
        &mut atoms,
        &bonds,
        &adjacency,
        &ring_closures_by_atom,
        &smiles_start_atoms,
    )?;
    let topology = TopologyBlock {
        adjacency,
        atoms,
        bonds,
        ..TopologyBlock::default()
    };
    topology.validate().map_err(|error| match error {
        cosmolkit_model::TopologyValidationError::StereoGroup(cause) => {
            SmilesParseError::StereoGroup(cause)
        }
        other => SmilesParseError::Model(other.to_string()),
    })?;
    let mut record = SmilesRecord {
        topology,
        coordinates: CoordinateBlock::default(),
        properties: MoleculeProperties::default(),
    };
    let mut name = preprocessed.name;
    let cx = preprocessed.cx_part;
    if !cx.is_empty() {
        // BEGIN RDKIT CPP FUNCTION handleCXPartAndName
        // RDKit✔️✔️:   std::string::const_iterator pos = cxPart.cbegin();
        // RDKit✔️✔️:   bool cxfailed = false;
        // RDKit✔️✔️:   if (params.allowCXSMILES) {
        // RDKit✔️✔️:     if (*pos == '|') {
        // RDKit✔️✔️:       try {
        // RDKit✔️✔️:         SmilesParseOps::parseCXExtensions(*res, cxPart, pos);
        // RDKit✔️✔️:       } catch (...) {
        // RDKit✔️✔️:         cxfailed = true;
        // RDKit✔️✔️:         if (params.strictCXSMILES) {
        // RDKit✔️✔️:           throw;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       res->setProp("_CXSMILES_Data", std::string(cxPart.cbegin(), pos));
        // RDKit✔️✔️:     } else if (params.strictCXSMILES && !params.parseName &&
        // RDKit✔️✔️:                pos != cxPart.cend()) {
        // RDKit✔️✔️:       throw RDKit::SmilesParseException(
        // RDKit✔️✔️:           "CXSMILES extension does not start with | and parseName=false");
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (!cxfailed && params.parseName && pos != cxPart.end()) {
        // RDKit✔️✔️:     std::string nmpart(pos, cxPart.cend());
        // RDKit✔️✔️:     name = boost::trim_copy(nmpart);
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION handleCXPartAndName
        if params.allow_cxsmiles && cx.starts_with('|') {
            // RDKit❗❌: New duplicate failures preserve source-completed CX
            // records. An extra syntax pass isolates them without extending
            // unchanged malformed-CX recovery; see the private adapter's cost.
            let atom_count = u32::try_from(record.topology.atoms.len()).map_err(|_| {
                SmilesParseError::Cx("CX atom count exceeds source unsigned32 domain".into())
            })?;
            let progress =
                cosmolkit_cx::parse_cx_extensions_progress_with_atom_window(&cx, 0, atom_count);
            if cx_lowering::has_enhanced_membership_collision(
                &progress,
                record.topology.atoms.len(),
            ) {
                let mut cursor = 0;
                let applied = cx_lowering::apply_enhanced_duplicate_progress(
                    &mut record,
                    &progress,
                    &mut cursor,
                );
                if let Err(error) = applied {
                    if params.strict_cxsmiles {
                        return Err(error);
                    }
                }
                record
                    .properties
                    .set_prop("_CXSMILES_Data", &cx[..cursor])?;
                // Both new rejection forms suppress name parsing in lax mode.
            } else {
                match parse_cx_extensions_with_atom_window(&cx, 0, record.topology.atoms.len()) {
                    Ok(parsed) => {
                        match cx_lowering::apply_cx_to_smiles_record(&mut record, &parsed) {
                            Ok(()) => {
                                record
                                    .properties
                                    .set_prop("_CXSMILES_Data", &cx[..parsed.consumed()])
                                    .map_err(|error| SmilesParseError::Model(error.to_string()))?;
                                if params.parse_name
                                    && let Some(parsed_name) = cx
                                        .get(parsed.consumed()..)
                                        .map(str::trim)
                                        .filter(|name| !name.is_empty())
                                {
                                    name = parsed_name.to_owned();
                                }
                            }
                            Err(error) if params.strict_cxsmiles => return Err(error),
                            Err(_) => record
                                .properties
                                .set_prop("_CXSMILES_Data", "")
                                .expect("the internal CXSMILES data property key is non-empty"),
                        }
                    }
                    Err(error) if params.strict_cxsmiles => {
                        return Err(SmilesParseError::Cx(error.to_string()));
                    }
                    Err(_) => record
                        .properties
                        .set_prop("_CXSMILES_Data", "")
                        .expect("the internal CXSMILES data property key is non-empty"),
                }
            }
        } else if params.allow_cxsmiles && params.strict_cxsmiles && !params.parse_name {
            return Err(SmilesParseError::Cx(
                "CXSMILES extension does not start with | and parseName=false".into(),
            ));
        } else if params.parse_name {
            name = cx.trim().to_owned();
        }
    }
    apply_parser_wedge_stereo(&mut record)?;
    apply_parser_3d_stereo(&mut record)?;
    apply_parser_atrop_stereo(&mut record)?;
    if complete {
        // Preserve the actual source-selected conformer identity across removeHs;
        // the hydrogen owner remaps coordinate rows without changing these IDs.
        let (conf, conf3d) = parser_conformers(&record.coordinates)?;
        let selected_conformer = conf.or(conf3d).map(|row| match row {
            CoordinateSource::TwoD(row) => row.id(),
            CoordinateSource::ThreeD(row) => row.id(),
        });
        let mut final_valence = None;
        let mut final_rings = None;
        if params.remove_hs {
            let removed = cosmolkit_core::remove_hydrogens_with_params(
                record.topology,
                record.coordinates,
                record.properties,
                &cosmolkit_core::RemoveHsParams {
                    update_explicit_count: true,
                    sanitize: params.sanitize,
                    ..Default::default()
                },
            )
            .map_err(SmilesParseError::ParserRemoveHydrogens)?;
            record = SmilesRecord {
                topology: removed.topology,
                coordinates: removed.coordinates,
                properties: removed.properties,
            };
            final_valence = removed.final_valence;
            final_rings = removed.final_rings;
        } else if params.sanitize {
            let sanitized = cosmolkit_core::sanitize_topology(
                &record.topology,
                &cosmolkit_core::SanitizeParams::default(),
            )
            .map_err(SmilesParseError::ParserSanitize)?;
            record.properties.clear_computed_props()?;
            if let Some(count) = sanitized.aromatic_ring_count {
                // Native numArom has the source setAromaticity integral value.
                record
                    .properties
                    .set_computed_prop("numArom", count as i32)?;
            }
            record.topology = sanitized.topology;
            final_valence = sanitized.final_valence;
            final_rings = sanitized.final_rings;
        }
        record = finalize_stereo::finalize_smiles_stereo_with_conformer(
            record,
            params,
            &mut final_valence,
            &mut final_rings,
            selected_conformer,
        )
        .map_err(SmilesParseError::ParserStereo)?;
        if let Some(valence) = final_valence {
            for (index, atom) in record.topology.atoms.iter_mut().enumerate() {
                atom.set_source_valence_facts(cosmolkit_model::SourceAtomValenceFacts {
                    explicit_valence: valence.explicit_valence[index] as i8,
                    implicit_valence: valence.implicit_hydrogens[index] as i8,
                });
            }
        }
        // Query-only CX atom records are rejected by the existing representation
        // boundary. No fake query completion is performed on concrete atoms.
        if record.properties.prop("_NeedsQueryScan").is_some() {
            return Err(SmilesParseError::UnsupportedCx(
                "CX atom-query completion requires the canonical QueryGraph",
            ));
        }
    }
    if !params.skip_cleanup {
        cleanup_after_parsing_impl(&mut record, complete)?;
    }
    if !name.is_empty() {
        record.properties = record.properties.with_name(&name);
    }
    record.topology.validate().map_err(|error| match error {
        cosmolkit_model::TopologyValidationError::StereoGroup(cause) => {
            SmilesParseError::StereoGroup(cause)
        }
        other => SmilesParseError::Model(other.to_string()),
    })?;
    record
        .coordinates
        .validate_for_atom_count(record.topology.atoms.len())
        .map_err(|error| SmilesParseError::Model(error.to_string()))?;
    Ok(record)
}

fn parser_conformers(
    coordinates: &CoordinateBlock,
) -> Result<(Option<CoordinateSource<'_>>, Option<CoordinateSource<'_>>), SmilesParseError> {
    // BEGIN RDKIT COMPLETE SOURCE MolFromSmiles conformer selection stage
    // RDKit✔️✔️:   // get a conformer
    // RDKit✔️✔️:   const Conformer *conf = nullptr, *conf3d = nullptr;
    // RDKit✔️✔️:   if (res && res->getNumConformers() > 0) {
    // RDKit✔️✔️:     for (unsigned int confId = 0; confId < res->getNumConformers(); ++confId) {
    // RDKit✔️✔️:       auto *testConf = &res->getConformer(confId);
    // RDKit✔️✔️:       if (!testConf->is3D()) {
    // RDKit✔️✔️:         if (conf == nullptr) {  // only take the first 2d conf
    // RDKit✔️✔️:           conf = testConf;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         if (conf3d == nullptr) {  // only take the first 3d conf
    // RDKit✔️✔️:           conf3d = testConf;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (conf != nullptr && conf3d != nullptr) {
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // END RDKIT COMPLETE SOURCE MolFromSmiles conformer selection stage
    // The source scans explicit numeric IDs 0..getNumConformers, not array
    // positions or default selection. Fresh CX producers assign these IDs.
    // is3D is independent of XYZ storage. No ID sorting, dimension guess or
    // copying; one source-shaped ID lookup per iteration with early exit.
    // O(C^2) worst case like the source's repeated linear ID getter; O(1)
    // extra storage, borrowed original rows. Empty collection returns no pair.
    let count = (coordinates.conformers_2d.len() + coordinates.conformers_3d.len()) as u32;
    let mut two_d = None;
    let mut three_d = None;
    for id in 0..count {
        let row = cx_writer::source_conformer_by_id(coordinates, id as i32)?;
        let is_3d = match row {
            CoordinateSource::TwoD(_) => false,
            CoordinateSource::ThreeD(row) => row.is_3d(),
        };
        if !is_3d {
            if two_d.is_none() {
                two_d = Some(row);
            }
        } else if three_d.is_none() {
            three_d = Some(row);
        }
        if two_d.is_some() && three_d.is_some() {
            break;
        }
    }
    Ok((two_d, three_d))
}

fn apply_parser_wedge_stereo(record: &mut SmilesRecord) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT COMPLETE SOURCE MolFromSmiles pending wedge stage
    // RDKit✔️✔️:     // we encountered a wedged bond in the CXSMILES,
    // RDKit✔️✔️:     // these need to be handled the same way they were in mol files
    // RDKit✔️✔️:     res->clearProp(SmilesParseOps::detail::_needsDetectAtomStereo);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (conf) {
    // RDKit✔️✔️:       MolOps::assignChiralTypesFromBondDirs(*res, conf->getId());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT COMPLETE SOURCE MolFromSmiles pending wedge stage
    // Consume only after handleCXPartAndName, clear pending before the core
    // assignment/error, and use the exact selected false-is3D conformer. The
    // normal native Point3D-backed 2D path is borrowed with no allocation.
    // Legacy detached XY rows need the existing exact XYZ(z=0) lift: O(V)
    // allocation versus the native borrowed Point3D row, a known worse cost.
    let (two_d, _) = parser_conformers(&record.coordinates)?;
    if record.properties.prop("_needsDetectAtomStereo").is_none() {
        return Ok(());
    }
    record.properties.clear_prop("_needsDetectAtomStereo")?;
    match two_d {
        Some(CoordinateSource::ThreeD(row)) => {
            cosmolkit_core::assign_chiral_types_from_bond_dirs(&mut record.topology, row, false)
                .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
        }
        Some(CoordinateSource::TwoD(row)) => {
            let xyz = cosmolkit_model::Conformer3D::new(
                row.id(),
                row.coordinates()
                    .iter()
                    .map(|xy| [xy[0], xy[1], 0.0])
                    .collect(),
                false,
            );
            cosmolkit_core::assign_chiral_types_from_bond_dirs(&mut record.topology, &xyz, false)
                .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
        }
        None => {}
    }
    Ok(())
}

fn apply_parser_3d_stereo(record: &mut SmilesRecord) -> Result<(), SmilesParseError> {
    // BEGIN COMPLETE PINNED MolFromSmiles source3D stage
    // RDKit✔️❌:   // if we read a 3D conformer, set the stereo:
    // RDKit✔️❌:   // if (res->getNumConformers() && res->getConformer().is3D()) {
    // RDKit✔️❌:   if (!conf && conf3d) {
    // RDKit✔️❌:     res->updatePropertyCache(false);
    // RDKit✔️❌:     MolOps::assignChiralTypesFrom3D(*res, conf3d->getId(), true);
    // RDKit✔️❌:   }
    // END COMPLETE PINNED MolFromSmiles source3D stage
    // Use the same genuine numeric-ID/flag selection. A false-is3D Point3D
    // row suppresses this stage even when a true 3D row is also present.
    // The one CORE cache kernel computes native signed8 cache values, then
    // the owned parser result retains those exact facts before geometry.
    // Errors propagate; this owned parse result is discarded on error, as
    // the source unique_ptr result is. No live molecule authority is exposed.
    // Cost: readonly cache result has two O(V) scalar buffers; the existing
    // detached structure transform owns a topology copy versus native in-place
    // mutation. Repeated borrowed selection adds O(C^2). Known worse cost.
    let (two_d, three_d) = parser_conformers(&record.coordinates)?;
    if two_d.is_some() {
        return Ok(());
    }
    let Some(CoordinateSource::ThreeD(conformer)) = three_d else {
        return Ok(());
    };
    let conformer_id = conformer.id() as i32;
    let valence = cosmolkit_core::assign_valence_with_options_for_topology(
        &record.topology,
        cosmolkit_core::ValenceModel::RdkitLike,
        false,
    )
    .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
    for (index, atom) in record.topology.atoms.iter_mut().enumerate() {
        atom.set_source_valence_facts(cosmolkit_model::SourceAtomValenceFacts {
            explicit_valence: valence.explicit_valence[index] as i8,
            implicit_valence: valence.implicit_hydrogens[index] as i8,
        });
    }
    // Native assignChiralTypesFrom3D clears the marker before perception.
    record.properties.clear_prop("_StereochemDone")?;
    let assignment = cosmolkit_core::assign_chiral_tags_from_structure(
        &record.topology,
        &record.coordinates,
        &valence,
        &cosmolkit_core::StructureTagParams {
            conformer_id,
            replace_existing_tags: true,
        },
    )
    .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?;
    record.topology = assignment.topology;
    Ok(())
}

fn apply_parser_atrop_stereo(record: &mut SmilesRecord) -> Result<(), SmilesParseError> {
    // BEGIN COMPLETE PINNED MolFromSmiles sourceatrop stage
    // RDKit❗❌:   if (conf) {
    // RDKit❗❌:     Atropisomers::detectAtropisomerChirality(*res, conf);
    // RDKit❗❌:   } else if (conf3d) {
    // RDKit❗❌:     Atropisomers::detectAtropisomerChirality(*res, conf3d);
    // RDKit❗❌:   } else {
    // RDKit❗❌:     Atropisomers::detectAtropisomerChirality(*res, nullptr);
    // RDKit❗❌:   }
    // END COMPLETE PINNED MolFromSmiles sourceatrop stage
    // Borrow the genuine first false-is3D row, else true-is3D row, else null.
    // CORE owns all cache/conjugation/hybridization/detection algorithms.
    // This source-owned parse result applies explicit detached writes; it is
    // discarded on error. The remaining full-function comparison includes
    // native pointer-set traversal/diagnostic order (currently unresolved).
    let (two_d, three_d) = parser_conformers(&record.coordinates)?;
    let conformer = two_d.or(three_d).map(|row| match row {
        CoordinateSource::TwoD(row) => cosmolkit_core::AtropisomerConformer::TwoD(row),
        CoordinateSource::ThreeD(row) => cosmolkit_core::AtropisomerConformer::ThreeD(row),
    });
    let assignment = cosmolkit_core::detect_atropisomer_chirality(&record.topology, conformer)
        .map_err(SmilesParseError::ParserAtropisomer)?;
    for (id, facts) in assignment.atom_valence_updates {
        record.topology.atoms[id.index()].set_source_valence_facts(facts);
    }
    if let Some(flags) = assignment.conjugated_bonds {
        for (bond, flag) in record.topology.bonds.iter_mut().zip(flags) {
            bond.set_conjugated(flag);
        }
    }
    if let Some(hybridization) = assignment.hybridization {
        for (atom, value) in record.topology.atoms.iter_mut().zip(hybridization.values) {
            atom.set_hybridization(value);
        }
    }
    for update in assignment.bond_updates {
        record.topology.bonds[update.bond.index()].set_stereo(update.stereo)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_model::{BondStereo, PropertyValue, StereoGroupKind, SubstanceGroupKind};

    #[test]
    fn bond_and_separator_tokens_require_the_source_grammar_operands() {
        // Pinned smiles.yy: mol BOND_TOKEN atomd / mol BOND_TOKEN ring_number /
        // mol SEPARATOR_TOKEN atomd. Neither repeated nor dangling tokens form mol.
        for smiles in [
            "C==C", "C#", "C..C", "C--C", "C/#C", "C=", "C/", "C~", "C->", "C.", ".C", "=CC",
            "C()", "C(.C)", "C=(O)", "C(1)CC1", "C(=)O", "C\\\\\\C", "C//C", "C/\\C",
        ] {
            assert!(
                matches!(
                    parse_smiles(smiles, &Default::default()),
                    Err(SmilesParseError::Syntax { .. })
                ),
                "must reject {smiles:?}"
            );
        }
        for smiles in [
            "",
            "C=C",
            "C#N",
            "C.C",
            "C(=O)O",
            "C-1CCCCC1",
            "C1.C1",
            "C(Cl)(F)Br",
            "C/C=C\\C",
            "C/C=C\\\\C",
            "C\\\\1CCCCC1",
            "N->[Cu]",
            "C~N",
        ] {
            parse_smiles(smiles, &Default::default())
                .unwrap_or_else(|error| panic!("must accept {smiles:?}: {error}"));
        }
    }

    fn string_property(value: Option<&PropertyValue>) -> Option<&str> {
        match value {
            Some(PropertyValue::String(value)) => Some(fixed_property_text(value)),
            _ => None,
        }
    }

    #[test]
    fn tilde_query_identity_survives_chain_branch_ring_selection_and_serialization() {
        use cosmolkit_model::{BondQueryPredicate, QueryNode};
        for (smiles, query_count) in [
            ("C~N", 1),
            ("C(~N)O", 1),
            ("C~1CC1", 1),
            ("C1CC~1", 1),
            ("C~1CC=1", 1),
            ("C=1CC~1", 0),
        ] {
            let record = parse_smiles(smiles, &Default::default()).unwrap();
            record.topology.validate().unwrap();
            assert_eq!(
                record
                    .topology
                    .bonds
                    .iter()
                    .filter(|bond| bond.query().is_some())
                    .count(),
                query_count,
                "{smiles}"
            );
            for bond in &record.topology.bonds {
                if let Some(query) = bond.query() {
                    assert_eq!(query, &QueryNode::predicate(BondQueryPredicate::Any));
                }
            }
            let copied = record.clone();
            assert_eq!(copied.topology, record.topology);
            let text = write_smiles(&copied).unwrap();
            let roundtrip = parse_smiles(fixed_property_text(&text), &Default::default()).unwrap();
            assert_eq!(
                roundtrip
                    .topology
                    .bonds
                    .iter()
                    .filter(|bond| bond.query().is_some())
                    .count(),
                query_count,
                "{smiles} -> {text:?}"
            );
        }
    }
    #[test]
    fn parses_and_writes_detached_ethanol() {
        let record = parse_smiles("CCO", &Default::default()).expect("parse");
        assert_eq!(record.topology.atoms.len(), 3);
        assert_eq!(
            (write_smiles(&record).expect("write")).as_bytes(),
            ("CCO").as_bytes()
        );
        assert_eq!(
            parse_smiles("C%12CC%12", &Default::default())
                .unwrap()
                .topology
                .bonds
                .len(),
            3
        );
    }

    #[test]
    fn parses_bracket_attributes_and_cx_coordinates() {
        let record = parse_smiles(
            "[13CH3+]CO |(0,0,0;1,0,0;2,0,0)| ethanol",
            &Default::default(),
        )
        .expect("parse");
        assert_eq!(record.topology.atoms[0].isotope(), Some(13));
        assert_eq!(record.topology.atoms[0].formal_charge(), 1);
        // Source: `parse_coords` sets `is3D = true` for a third token but then
        // `conf->set3D(is3D && hasNonZeroZCoords(*conf))`, and
        // `get_coords_block` emits Z only for `conf.is3D()`. An all-zero Z
        // column sets is3D=false while retaining the source Point3D rows.
        assert_eq!(record.coordinates.conformers_3d.len(), 1);
        assert_eq!(record.coordinates.conformers_2d.len(), 0);
        assert!(!record.coordinates.conformers_3d[0].is_3d());
        assert_eq!(
            record.coordinates.conformers_3d[0].coordinates(),
            &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
        );
        assert_eq!(
            record.properties.name().map(|value| value.as_bytes()),
            Some("ethanol".as_bytes())
        );
    }

    #[test]
    fn implicit_bonds_follow_rdkit_aromatic_endpoint_rules() {
        let record = parse_smiles("cc.cC.Cc.c-c.[c][n]", &Default::default()).unwrap();
        let orders = record
            .topology
            .bonds
            .iter()
            .map(Bond::order)
            .collect::<Vec<_>>();
        assert_eq!(
            orders,
            [
                BondOrder::Aromatic,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Single,
                BondOrder::Aromatic,
            ]
        );
        assert!(
            record
                .topology
                .atoms
                .iter()
                .take(2)
                .all(|atom| !atom.no_implicit())
        );
    }

    #[test]
    fn preprocesses_plain_and_cx_names_like_rdkit() {
        let plain = parse_smiles("CC   ethanol sample", &Default::default()).unwrap();
        assert_eq!(
            plain.properties.name().map(|value| value.as_bytes()),
            Some("ethanol sample".as_bytes())
        );

        let cx = parse_smiles("CC |$foo;bar$| named sample", &Default::default()).unwrap();
        assert_eq!(
            cx.properties.name().map(|value| value.as_bytes()),
            Some("named sample".as_bytes())
        );
        assert_eq!(
            string_property(cx.topology.atoms[0].prop("atomLabel")),
            Some("foo")
        );

        let no_cx = parse_smiles(
            "CC |not cx when disabled|",
            &SmilesParseParams {
                allow_cxsmiles: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            no_cx.properties.name().map(|value| value.as_bytes()),
            Some("|not cx when disabled|".as_bytes())
        );
    }

    #[test]
    fn parse_name_false_preserves_rdkit_strict_cx_errors() {
        let error = parse_smiles(
            "CC ethanol",
            &SmilesParseParams {
                parse_name: false,
                ..Default::default()
            },
        )
        .unwrap_err();
        assert!(matches!(error, SmilesParseError::Cx(_)));

        let error = parse_smiles(
            "CC ethanol",
            &SmilesParseParams {
                allow_cxsmiles: false,
                parse_name: false,
                ..Default::default()
            },
        )
        .unwrap_err();
        assert!(matches!(
            error,
            SmilesParseError::Unsupported { token: ' ', .. }
        ));
    }

    #[test]
    fn ring_closures_preserve_first_specification_orientation_and_direction() {
        let first_double = parse_smiles("C=1CCCCC-1", &Default::default()).unwrap();
        let closure = &first_double.topology.bonds[5];
        assert_eq!(closure.order(), BondOrder::Double);
        assert_eq!(
            (closure.begin(), closure.end()),
            (AtomId::new(0), AtomId::new(5))
        );

        let first_single = parse_smiles("C-1CCCCC=1", &Default::default()).unwrap();
        let closure = &first_single.topology.bonds[5];
        assert_eq!(closure.order(), BondOrder::Single);
        assert_eq!(
            (closure.begin(), closure.end()),
            (AtomId::new(0), AtomId::new(5))
        );

        let opening_direction = parse_smiles("C/1CCCCC1", &Default::default()).unwrap();
        let closure = &opening_direction.topology.bonds[5];
        assert_eq!(
            (closure.begin(), closure.end()),
            (AtomId::new(5), AtomId::new(0))
        );
        assert_eq!(closure.direction(), BondDirection::EndDownRight);

        let closing_direction = parse_smiles("C1CCCCC/1", &Default::default()).unwrap();
        let closure = &closing_direction.topology.bonds[5];
        assert_eq!(
            (closure.begin(), closure.end()),
            (AtomId::new(5), AtomId::new(0))
        );
        assert_eq!(closure.direction(), BondDirection::EndUpRight);

        let aromatic = parse_smiles("c1ccccc1", &Default::default()).unwrap();
        let closure = &aromatic.topology.bonds[5];
        assert_eq!(closure.order(), BondOrder::Aromatic);
        assert!(closure.is_aromatic());
        assert_eq!(
            (closure.begin(), closure.end()),
            (AtomId::new(5), AtomId::new(0))
        );
    }

    #[test]
    fn dative_and_quadruple_ring_closures_preserve_source_orientation() {
        for (input, begin, end) in [
            ("N->1CCCCC1", 0, 5),
            ("N<-1CCCCC1", 5, 0),
            ("N1CCCCC->1", 5, 0),
            ("N1CCCCC<-1", 0, 5),
        ] {
            let record = parse_smiles(input, &Default::default()).unwrap();
            let closure = &record.topology.bonds[5];
            assert_eq!(closure.order(), BondOrder::Dative, "{input}");
            assert_eq!(
                (closure.begin(), closure.end()),
                (AtomId::new(begin), AtomId::new(end)),
                "{input}"
            );
        }

        for (input, begin, end) in [("C$1CCCCC1", 0, 5), ("C1CCCCC$1", 5, 0)] {
            let record = parse_smiles(input, &Default::default()).unwrap();
            let closure = &record.topology.bonds[5];
            assert_eq!(closure.order(), BondOrder::Quadruple, "{input}");
            assert_eq!(
                (closure.begin(), closure.end()),
                (AtomId::new(begin), AtomId::new(end)),
                "{input}"
            );
        }

        let right = parse_smiles("N->C", &Default::default()).unwrap();
        assert_eq!(right.topology.bonds[0].order(), BondOrder::Dative);
        assert_eq!(
            (
                right.topology.bonds[0].begin(),
                right.topology.bonds[0].end()
            ),
            (AtomId::new(0), AtomId::new(1))
        );
        let left = parse_smiles("N<-C", &Default::default()).unwrap();
        assert_eq!(left.topology.bonds[0].order(), BondOrder::Dative);
        assert_eq!(
            (left.topology.bonds[0].begin(), left.topology.bonds[0].end()),
            (AtomId::new(1), AtomId::new(0))
        );
    }

    #[test]
    fn ring_closures_support_extended_labels_and_source_ordering() {
        for input in ["C%12CC%12", "C%(0)CC%(0)", "C%(123)CC%(123)"] {
            let record = parse_smiles(input, &Default::default()).unwrap();
            assert_eq!(record.topology.bonds.len(), 3, "{input}");
            assert_eq!(
                (
                    record.topology.bonds[2].begin(),
                    record.topology.bonds[2].end()
                ),
                (AtomId::new(2), AtomId::new(0)),
                "{input}"
            );
        }

        let repeated = parse_smiles("C1CC11CC1", &Default::default()).unwrap();
        assert_eq!(repeated.topology.bonds.len(), 6);
        assert_eq!(
            (
                repeated.topology.bonds[4].begin(),
                repeated.topology.bonds[4].end()
            ),
            (AtomId::new(2), AtomId::new(0))
        );
        assert_eq!(
            (
                repeated.topology.bonds[5].begin(),
                repeated.topology.bonds[5].end()
            ),
            (AtomId::new(4), AtomId::new(2))
        );

        for input in ["C%01CC%01", "C%()CC%()", "C%(123456)CC%(123456)"] {
            assert!(
                matches!(
                    parse_smiles(input, &Default::default()),
                    Err(SmilesParseError::Syntax { .. })
                ),
                "{input}"
            );
        }
    }

    #[test]
    fn ring_closures_reject_self_and_existing_bond_pairs() {
        for input in ["C11", "C1C1"] {
            let error = parse_smiles(input, &Default::default()).unwrap_err();
            assert!(
                matches!(error, SmilesParseError::Syntax { .. }),
                "{error:?}"
            );
        }
    }

    #[test]
    fn ring_closure_branch_status_and_cx_indices_follow_source_order() {
        let branch_after = parse_smiles("F[C@](Cl)1CCCC1", &Default::default()).unwrap();
        assert_eq!(
            branch_after.topology.atoms[1].chiral_tag(),
            ChiralTag::TetrahedralCw
        );
        let branch_before = parse_smiles("F[C@]1(Cl)CCCC1", &Default::default()).unwrap();
        assert_eq!(
            branch_before.topology.atoms[1].chiral_tag(),
            ChiralTag::TetrahedralCcw
        );

        let cx = parse_smiles("C1CC1CC |Z:2|", &Default::default()).unwrap();
        assert_eq!(cx.topology.bonds.len(), 5);
        assert_eq!(cx.topology.bonds[4].order(), BondOrder::Zero);
        assert_eq!(
            (cx.topology.bonds[4].begin(), cx.topology.bonds[4].end()),
            (AtomId::new(2), AtomId::new(0))
        );
    }

    #[test]
    fn bracket_aromatic_tellurium_preserves_source_atom_state() {
        let record = parse_smiles("c1cc[te]c1", &Default::default()).unwrap();
        assert_eq!(record.topology.atoms.len(), 5);
        assert_eq!(record.topology.bonds.len(), 5);
        let tellurium = &record.topology.atoms[3];
        assert_eq!(tellurium.element(), Element::TE);
        assert_eq!(tellurium.atomic_number(), 52);
        assert!(tellurium.is_aromatic());
        assert!(tellurium.no_implicit());
        assert_eq!(tellurium.explicit_hydrogens(), 0);
        assert!(
            record
                .topology
                .bonds
                .iter()
                .all(|bond| bond.order() == BondOrder::Aromatic)
        );

        let decorated = parse_smiles("[125teH+:7]", &Default::default()).unwrap();
        let atom = &decorated.topology.atoms[0];
        assert_eq!(atom.element(), Element::TE);
        assert!(atom.is_aromatic());
        assert_eq!(atom.isotope(), Some(125));
        assert_eq!(atom.explicit_hydrogens(), 1);
        assert_eq!(atom.formal_charge(), 1);
        assert_eq!(atom.atom_map(), Some(7));
        assert!(
            !parse_smiles("[Te]", &Default::default())
                .unwrap()
                .topology
                .atoms[0]
                .is_aromatic()
        );
        assert!(parse_smiles("c1cctec1", &Default::default()).is_err());
    }

    #[test]
    fn bracket_parser_preserves_rdkit_atom_fields_and_repeated_charges() {
        let record = parse_smiles(
            "[13C@TH2H3-:5].[2HH1-].[Pt++].[Cl--].[#6H2-:7].[seH]",
            &Default::default(),
        )
        .unwrap();
        let atoms = &record.topology.atoms;

        assert_eq!(atoms[0].isotope(), Some(13));
        assert_eq!(atoms[0].chiral_tag(), ChiralTag::TetrahedralCw);
        assert_eq!(atoms[0].chiral_permutation(), None);
        assert_eq!(atoms[0].explicit_hydrogens(), 3);
        assert_eq!(atoms[0].formal_charge(), -1);
        assert_eq!(atoms[0].atom_map(), Some(5));
        assert!(atoms[0].no_implicit());

        assert_eq!(atoms[1].element(), Element::H);
        assert_eq!(atoms[1].isotope(), Some(2));
        assert_eq!(atoms[1].explicit_hydrogens(), 1);
        assert_eq!(atoms[1].formal_charge(), -1);
        assert_eq!(atoms[2].element(), Element::PT);
        assert_eq!(atoms[2].formal_charge(), 2);
        assert_eq!(atoms[3].element(), Element::CL);
        assert_eq!(atoms[3].formal_charge(), -2);
        assert_eq!(atoms[4].element(), Element::C);
        assert_eq!(atoms[4].explicit_hydrogens(), 2);
        assert_eq!(atoms[4].atom_map(), Some(7));
        assert_eq!(atoms[5].element(), Element::SE);
        assert!(atoms[5].is_aromatic());
        assert_eq!(atoms[5].explicit_hydrogens(), 1);
        assert!(atoms.iter().all(Atom::no_implicit));
    }

    #[test]
    fn parser_adjusts_tetrahedral_tags_to_rdkit_storage_order() {
        for (input, atom_index, expected) in [
            ("[C@H](F)(Cl)Br", 0, ChiralTag::TetrahedralCw),
            ("[C@@H](F)(Cl)Br", 0, ChiralTag::TetrahedralCcw),
            ("F[C@H](Cl)Br", 1, ChiralTag::TetrahedralCcw),
            ("Br[C@@H](Cl)F", 1, ChiralTag::TetrahedralCw),
            ("N[C@](F)(Cl)Br", 1, ChiralTag::TetrahedralCcw),
            ("F[C@]1(Br)CCO1", 1, ChiralTag::TetrahedralCcw),
            ("F[C@]1(CCO1)Br", 1, ChiralTag::TetrahedralCcw),
        ] {
            let record = parse_smiles(input, &Default::default()).unwrap();
            assert_eq!(
                record.topology.atoms[atom_index].chiral_tag(),
                expected,
                "{input}"
            );
        }
    }

    #[test]
    fn bracket_parser_preserves_non_tetrahedral_chiral_classes() {
        let record =
            parse_smiles("[C@AL1].[C@SP3].[C@TB20].[C@OH30]", &Default::default()).unwrap();
        let atoms = &record.topology.atoms;
        assert_eq!(atoms[0].chiral_tag(), ChiralTag::Allene);
        assert_eq!(atoms[0].chiral_permutation(), Some(1));
        assert_eq!(atoms[1].chiral_tag(), ChiralTag::SquarePlanar);
        assert_eq!(atoms[1].chiral_permutation(), Some(3));
        assert_eq!(atoms[2].chiral_tag(), ChiralTag::TrigonalBipyramidal);
        assert_eq!(atoms[2].chiral_permutation(), Some(20));
        assert_eq!(atoms[3].chiral_tag(), ChiralTag::Octahedral);
        assert_eq!(atoms[3].chiral_permutation(), Some(30));
    }

    #[test]
    fn parser_converts_nontetrahedral_permutations_to_rdkit_storage_order() {
        for (input, atom_index, expected_permutation) in [
            ("[Pt@SP1](F)(Cl)(Br)I", 0, 1),
            ("I[Pt@SP1](Br)(Cl)F", 1, 1),
            ("[Pt@SP1](F)(Cl)Br", 0, 1),
            ("F[Pt@SP1](Cl)Br", 1, 2),
            ("[P@TB1](F)(Cl)(Br)(I)N", 0, 1),
            ("N[P@TB1](I)(Br)(Cl)F", 1, 1),
            ("[P@TB1](F)(Cl)(Br)I", 0, 18),
            ("F[P@TB1](Cl)(Br)I", 1, 3),
            ("[Co@OH1](F)(Cl)(Br)(I)(N)O", 0, 1),
            ("O[Co@OH1](N)(I)(Br)(Cl)F", 1, 1),
            ("[Co@OH1](F)(Cl)(Br)(I)N", 0, 22),
            ("F[Co@OH1](Cl)(Br)(I)N", 1, 3),
        ] {
            let record = parse_smiles(input, &Default::default()).unwrap();
            assert_eq!(
                record.topology.atoms[atom_index].chiral_permutation(),
                Some(expected_permutation),
                "{input}"
            );
        }
    }

    #[test]
    fn bracket_parser_rejects_source_invalid_charge_and_chiral_forms() {
        for input in ["[C@TH0]", "[C@TH3]", "[C@SP4]", "[C+++]"] {
            assert!(
                matches!(
                    parse_smiles(input, &Default::default()),
                    Err(SmilesParseError::Atom { .. })
                ),
                "{input}"
            );
        }
    }

    #[test]
    fn lowers_cx_coordinate_hydrogen_zero_bonds_and_radicals_like_rdkit() {
        let coordinate = parse_smiles("NO |C:1.0|", &Default::default()).unwrap();
        assert_eq!(coordinate.topology.bonds[0].order(), BondOrder::Dative);
        assert_eq!(coordinate.topology.bonds[0].begin(), AtomId::new(1));
        assert_eq!(coordinate.topology.bonds[0].end(), AtomId::new(0));

        let hydrogen = parse_smiles("NO |H:1.0|", &Default::default()).unwrap();
        assert_eq!(hydrogen.topology.bonds[0].order(), BondOrder::Hydrogen);
        assert_eq!(hydrogen.topology.bonds[0].begin(), AtomId::new(1));

        let zero = parse_smiles("CC~CC |Z:1|", &Default::default()).unwrap();
        assert_eq!(zero.topology.bonds[1].order(), BondOrder::Zero);

        let radicals = parse_smiles("CCC |^1:0,^4:1,^7:2|", &Default::default()).unwrap();
        assert_eq!(radicals.topology.atoms[0].radical_electrons(), 1);
        assert_eq!(radicals.topology.atoms[1].radical_electrons(), 2);
        assert_eq!(radicals.topology.atoms[2].radical_electrons(), 3);
    }

    #[test]
    fn lowers_cx_enhanced_wedge_and_double_bond_stereo_like_rdkit() {
        let enhanced = parse_smiles("CCC |o1:1,o1:2|", &Default::default()).unwrap();
        assert_eq!(enhanced.topology.stereo_groups.len(), 1);
        assert_eq!(
            enhanced.topology.stereo_groups[0].kind(),
            StereoGroupKind::Or
        );
        assert_eq!(enhanced.topology.stereo_groups[0].id(), Some(1));
        assert_eq!(
            enhanced.topology.stereo_groups[0].atoms(),
            &[AtomId::new(1), AtomId::new(2)]
        );

        let wedge = parse_smiles("CC |wU:1.0|", &Default::default()).unwrap();
        assert_eq!(wedge.topology.bonds[0].begin(), AtomId::new(1));
        assert_eq!(wedge.topology.bonds[0].end(), AtomId::new(0));
        assert_eq!(
            wedge.topology.bonds[0].prop("_MolFileBondCfg"),
            Some(&cosmolkit_model::PropertyValue::UInt(1))
        );
        assert_eq!(wedge.properties.prop("_needsDetectAtomStereo"), None);
        assert_eq!(wedge.topology.atoms[1].chiral_tag(), ChiralTag::Unspecified);

        let double = parse_smiles("CC=CC |c:1|", &Default::default()).unwrap();
        assert_eq!(double.topology.bonds[1].stereo(), BondStereo::Cis);
        assert_eq!(
            double.topology.bonds[1].stereo_atoms(),
            Some([AtomId::new(0), AtomId::new(3)])
        );
        assert_eq!(
            double.properties.prop("_needsDetectBondStereo"),
            Some(&cosmolkit_model::PropertyValue::Int(1))
        );
    }

    #[test]
    fn lowers_cx_linknodes_sgroups_and_variable_attachments_like_rdkit() {
        let link = parse_smiles("C1CC1 |LN:1:1.3|", &Default::default()).unwrap();
        assert_eq!(
            link.properties.prop("_molLinkNodes"),
            Some(&cosmolkit_model::PropertyValue::String(
                "1 3 2 2 1 2 3".into()
            ))
        );

        let data = parse_smiles("CCO |SgD:2,1:FIELD:info::::|", &Default::default()).unwrap();
        assert_eq!(data.topology.substance_groups.len(), 1);
        let data_group = &data.topology.substance_groups[0];
        assert_eq!(data_group.kind(), &SubstanceGroupKind::Data);
        assert_eq!(data_group.atoms(), &[AtomId::new(2), AtomId::new(1)]);
        assert_eq!(
            data_group
                .props()
                .get("FIELDNAME".as_bytes())
                .map(|value| cosmolkit_core::property_value_to_string(value).unwrap()),
            Some(cosmolkit_model::PropertyText::from("FIELD"))
        );
        assert_eq!(data_group.data_fields(), &["info".into()]);

        let polymer = parse_smiles("CC |Sg:n:0::ht|", &Default::default()).unwrap();
        assert_eq!(polymer.topology.substance_groups.len(), 1);
        assert_eq!(
            polymer.topology.substance_groups[0].kind(),
            &SubstanceGroupKind::StructuralRepeatUnit
        );

        let attachment = parse_smiles("CO*.C1=CC=NC=C1 |m:2:3.5.4|", &Default::default()).unwrap();
        assert_eq!(
            string_property(attachment.topology.bonds[1].prop("_MolFileBondEndPts")),
            Some("(3 4 6 5)")
        );
        assert_eq!(
            string_property(attachment.topology.bonds[1].prop("_MolFileBondAttach")),
            Some("ANY")
        );
    }

    #[test]
    fn query_only_cx_records_are_atomic_strict_errors_and_nonstrict_recovery() {
        let strict = parse_smiles("CC |rb:0:0|", &Default::default()).unwrap_err();
        assert!(
            matches!(strict, SmilesParseError::UnsupportedCx(_)),
            "{strict:?}"
        );

        let non_strict = parse_smiles(
            "CC |rb:0:0| ethanol",
            &SmilesParseParams {
                strict_cxsmiles: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(non_strict.topology.atoms.len(), 2);
        assert_eq!(
            non_strict.properties.prop("_CXSMILES_Data"),
            Some(&cosmolkit_model::PropertyValue::String("".into()))
        );
        assert_eq!(non_strict.properties.name(), None);
    }

    #[test]
    fn cx_lowering_cleans_parser_only_bond_and_sgroup_indices() {
        let record = parse_smiles("CC |Sg:n:0::ht|", &Default::default()).unwrap();
        assert!(
            record
                .topology
                .bonds
                .iter()
                .all(|bond| bond.prop(CXSMILES_BOND_IDX_PROP).is_none())
        );
        assert!(record.topology.substance_groups.iter().all(|group| {
            group
                .props()
                .get("_cxsmilesindex".as_bytes())
                .map(|value| cosmolkit_core::property_value_to_string(value).unwrap())
                == Some(cosmolkit_model::PropertyText::from("0"))
        }));
    }
}

fn take_smiles_bond_source_index(counter: &mut u32) -> u32 {
    // RDKit✔️✔️:   ++numBondsParsed;
    // Unsigned32 increment wraps modulo2^32. Ordinary callers discard the old
    // value; ring reservation callers retain it for the delayed closing edge.
    // One copy/add, O(1) scratch; this helper does not assign any property.
    let previous = *counter;
    *counter = counter.wrapping_add(1);
    previous
}

#[cfg(test)]
mod uint_complete_source_condition_cells {
    use super::*;
    // FROZEN UINT CONDITION: COUNTER_0
    #[test]
    fn uint_cell_counter_0_lib() {
        let mut counter = 0_u32;
        let previous = take_smiles_bond_source_index(&mut counter);
        assert_eq!(previous, 0_u32);
        assert_eq!(counter, 1_u32);
        let mut bonds = vec![];
        push_smiles_bond(
            &mut bonds,
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            Some(previous),
        );
        assert_eq!(
            bonds[0].prop("_cxsmilesBondIdx"),
            Some(&cosmolkit_model::PropertyValue::UInt(0_u32))
        );
    }
    // FROZEN UINT CONDITION: COUNTER_1
    #[test]
    fn uint_cell_counter_1_lib() {
        let mut counter = 1_u32;
        let previous = take_smiles_bond_source_index(&mut counter);
        assert_eq!(previous, 1_u32);
        assert_eq!(counter, 2_u32);
        let mut bonds = vec![];
        push_smiles_bond(
            &mut bonds,
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            Some(previous),
        );
        assert_eq!(
            bonds[0].prop("_cxsmilesBondIdx"),
            Some(&cosmolkit_model::PropertyValue::UInt(1_u32))
        );
    }
    // FROZEN UINT CONDITION: COUNTER_2147483646
    #[test]
    fn uint_cell_counter_2147483646_lib() {
        let mut counter = 2147483646_u32;
        let previous = take_smiles_bond_source_index(&mut counter);
        assert_eq!(previous, 2147483646_u32);
        assert_eq!(counter, 2147483647_u32);
        let mut bonds = vec![];
        push_smiles_bond(
            &mut bonds,
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            Some(previous),
        );
        assert_eq!(
            bonds[0].prop("_cxsmilesBondIdx"),
            Some(&cosmolkit_model::PropertyValue::UInt(2147483646_u32))
        );
    }
    // FROZEN UINT CONDITION: COUNTER_2147483647
    #[test]
    fn uint_cell_counter_2147483647_lib() {
        let mut counter = 2147483647_u32;
        let previous = take_smiles_bond_source_index(&mut counter);
        assert_eq!(previous, 2147483647_u32);
        assert_eq!(counter, 2147483648_u32);
        let mut bonds = vec![];
        push_smiles_bond(
            &mut bonds,
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            Some(previous),
        );
        assert_eq!(
            bonds[0].prop("_cxsmilesBondIdx"),
            Some(&cosmolkit_model::PropertyValue::UInt(2147483647_u32))
        );
    }
    // FROZEN UINT CONDITION: COUNTER_2147483648
    #[test]
    fn uint_cell_counter_2147483648_lib() {
        let mut counter = 2147483648_u32;
        let previous = take_smiles_bond_source_index(&mut counter);
        assert_eq!(previous, 2147483648_u32);
        assert_eq!(counter, 2147483649_u32);
        let mut bonds = vec![];
        push_smiles_bond(
            &mut bonds,
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            Some(previous),
        );
        assert_eq!(
            bonds[0].prop("_cxsmilesBondIdx"),
            Some(&cosmolkit_model::PropertyValue::UInt(2147483648_u32))
        );
    }
    // FROZEN UINT CONDITION: COUNTER_4294967295
    #[test]
    fn uint_cell_counter_4294967295_lib() {
        let mut counter = 4294967295_u32;
        let previous = take_smiles_bond_source_index(&mut counter);
        assert_eq!(previous, 4294967295_u32);
        assert_eq!(counter, 0_u32);
        let mut bonds = vec![];
        push_smiles_bond(
            &mut bonds,
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            Some(previous),
        );
        assert_eq!(
            bonds[0].prop("_cxsmilesBondIdx"),
            Some(&cosmolkit_model::PropertyValue::UInt(4294967295_u32))
        );
    }
}

#[doc(hidden)]
pub use cx_writer::select_cx_coordinates_from_sets;

#[doc(hidden)]
pub use cx_writer::format_cx_coordinate;

#[doc(hidden)]
pub use cx_writer::assign_stereo_group_ids;

#[cfg(test)]
fn fixed_property_text(value: &cosmolkit_model::PropertyText) -> &str {
    std::str::from_utf8(value.as_bytes()).expect("original fixed fixture text is UTF8")
}

#[cfg(test)]
mod source_parser_conformer_tests {
    use super::parse_smiles_complete_source as parse_smiles;
    use super::*;
    use cosmolkit_model::{Conformer3D, CoordinateDimension};

    #[test]
    fn parser_uses_numeric_source_ids_and_explicit_is3d_not_default_front_or_geometry() {
        let rows = CoordinateBlock {
            conformers_3d: vec![
                Conformer3D::new(1, vec![[1.0, 2.0, 3.0]], false),
                Conformer3D::new(0, vec![[4.0, 5.0, 6.0]], false),
                Conformer3D::new(2, vec![[0.0, 0.0, 0.0]], true),
            ],
            source_conformer_order: Some(vec![CoordinateDimension::ThreeD; 3]),
            ..Default::default()
        };
        let before = rows.clone();
        let (two_d, three_d) = parser_conformers(&rows).unwrap();
        match two_d.unwrap() {
            CoordinateSource::ThreeD(row) => assert!(std::ptr::eq(row, &rows.conformers_3d[1])),
            _ => panic!("native 2D flag lives in XYZ storage"),
        }
        match three_d.unwrap() {
            CoordinateSource::ThreeD(row) => assert!(std::ptr::eq(row, &rows.conformers_3d[2])),
            _ => panic!("explicit source 3D flag must remain independent of geometry"),
        }
        match cx_writer::source_conformer_by_id(&rows, -1).unwrap() {
            CoordinateSource::ThreeD(row) => assert!(std::ptr::eq(row, &rows.conformers_3d[0])),
            _ => panic!("native getter uses actual insertion front"),
        }
        assert!(matches!(
            select_cx_coordinates(&rows, CxCoordinateSelection::Auto),
            Err(SmilesParseError::AmbiguousCoordinateSelection {
                two_d_count: 0,
                three_d_count: 3,
            })
        ));
        assert_eq!(rows, before);
    }

    #[test]
    fn parser_preserves_source_missing_numeric_id_error_and_empty_case() {
        let empty = CoordinateBlock::default();
        let (two_d, three_d) = parser_conformers(&empty).unwrap();
        assert!(two_d.is_none() && three_d.is_none());
        let rows = CoordinateBlock {
            conformers_3d: vec![
                Conformer3D::new(0, vec![[0.0, 0.0, 0.0]], false),
                Conformer3D::new(19, vec![[1.0, 2.0, 3.0]], false),
            ],
            ..Default::default()
        };
        let before = rows.clone();
        assert!(
            matches!(parser_conformers(&rows), Err(SmilesParseError::Model(message)) if message == "Can't find conformation with ID: 1")
        );
        assert_eq!(rows, before);
    }

    #[test]
    fn parser_wedge_uses_first_2d_flag_even_after_a_3d_conformer() {
        let two_d = "(-3.9163,5.4767,;-3.9163,3.9367,;-2.5826,3.1667,;-5.25,3.1667,)";
        let three_d = "(-3.9163,5.4767,1;-3.9163,3.9367,1;-2.5826,3.1667,1;-5.25,3.1667,1)";
        for blocks in [format!("{three_d}{two_d}"), format!("{two_d}{three_d}")] {
            for (wedge, expected) in [
                ("wU", ChiralTag::TetrahedralCw),
                ("wD", ChiralTag::TetrahedralCcw),
            ] {
                // Original exact native fixture geometry and opposite wedge directions.
                let input = format!("CC(O)Cl |{blocks},{wedge}:1.0|");
                let record = parse_smiles(
                    &input,
                    &SmilesParseParams {
                        sanitize: false,
                        remove_hs: false,
                        ..Default::default()
                    },
                )
                .unwrap();
                assert_eq!(record.topology.atoms[1].chiral_tag(), expected, "{input}");
                assert_eq!(record.topology.atoms[1].explicit_hydrogens(), 1);
                assert_eq!(record.properties.prop("_needsDetectAtomStereo"), None);
                assert_eq!(record.coordinates.conformers_3d.len(), 2);
                assert!(record.coordinates.conformers_2d.is_empty());
                assert_eq!(
                    record.coordinates.conformers_3d[0].is_3d(),
                    blocks.starts_with(three_d)
                );
                assert_eq!(
                    record.coordinates.conformers_3d[1].is_3d(),
                    !blocks.starts_with(three_d)
                );
            }
        }
    }
}

#[cfg(test)]
mod source_parser_3d_tests {
    use super::parse_smiles_complete_source as parse_smiles;
    use super::*;
    use cosmolkit_model::{PropertyValue, SourceAtomValenceFacts};
    const XYZ: &str = "(1,0,0;0,0,0;0,1,0;0,0,1)";
    const XY: &str = "(1,0,0;0,0,0;0,1,0;0,0,0)";
    fn params() -> SmilesParseParams {
        SmilesParseParams {
            sanitize: false,
            remove_hs: false,
            ..Default::default()
        }
    }
    #[test]
    fn true_3d_updates_all_cache_rows_then_emits_native_tetrahedral_int_property() {
        let record = parse_smiles(&format!("FC(Cl)Br |{XYZ}|"), &params()).unwrap();
        assert_eq!(
            record.topology.atoms[1].chiral_tag(),
            ChiralTag::TetrahedralCcw
        );
        assert_eq!(
            record.topology.atoms[1].prop("_NonExplicit3DChirality"),
            Some(&PropertyValue::Int(1))
        );
        assert_eq!(record.topology.atoms[1].explicit_hydrogens(), 0);
        for (index, atom) in record.topology.atoms.iter().enumerate() {
            assert_eq!(
                atom.source_valence_facts(),
                SourceAtomValenceFacts {
                    explicit_valence: if index == 1 { 3 } else { 1 },
                    implicit_valence: if index == 1 { 1 } else { 0 },
                }
            );
        }
        assert!(record.coordinates.conformers_3d[0].is_3d());
        assert_eq!(
            record.coordinates.conformers_3d[0].coordinates(),
            &[
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0]
            ]
        );
    }
    #[test]
    fn any_false_is3d_conformer_suppresses_3d_stage_in_both_source_orders() {
        for blocks in [format!("{XYZ}{XY}"), format!("{XY}{XYZ}")] {
            let mut record = parse_smiles(&format!("FC(Cl)Br |{blocks}|"), &params()).unwrap();
            for atom in &record.topology.atoms {
                assert_eq!(atom.chiral_tag(), ChiralTag::Unspecified);
                assert_eq!(
                    atom.source_valence_facts(),
                    SourceAtomValenceFacts::UNINITIALIZED
                );
            }
            record.properties.set_prop("_StereochemDone", "0").unwrap();
            let before = record.clone();
            apply_parser_3d_stereo(&mut record).unwrap();
            assert_eq!(record, before);
            assert_eq!(
                record.coordinates.conformers_3d[0].is_3d(),
                blocks.starts_with(XYZ)
            );
            assert_eq!(
                record.coordinates.conformers_3d[1].is_3d(),
                !blocks.starts_with(XYZ)
            );
            assert!(record.coordinates.conformers_2d.is_empty());
        }
    }
    #[test]
    fn true_3d_clears_done_by_presence_and_replaces_existing_tag_without_new_nonexplicit_marker() {
        let mut record = parse_smiles(&format!("FC(Cl)Br |{XYZ}|"), &params()).unwrap();
        record.properties.set_prop("_StereochemDone", "0").unwrap();
        record.topology.atoms[1].set_chiral_tag(ChiralTag::TetrahedralCw);
        record.topology.atoms[1]
            .clear_prop("_NonExplicit3DChirality")
            .unwrap();
        for atom in &mut record.topology.atoms {
            atom.set_source_valence_facts(SourceAtomValenceFacts::UNINITIALIZED);
        }
        let coordinates_before = record.coordinates.clone();
        apply_parser_3d_stereo(&mut record).unwrap();
        assert_eq!(record.properties.prop("_StereochemDone"), None);
        assert_eq!(
            record.topology.atoms[1].chiral_tag(),
            ChiralTag::TetrahedralCcw
        );
        assert_eq!(
            record.topology.atoms[1].prop("_NonExplicit3DChirality"),
            None
        );
        assert_eq!(
            record.topology.atoms[1].source_valence_facts(),
            SourceAtomValenceFacts {
                explicit_valence: 3,
                implicit_valence: 1
            }
        );
        assert_eq!(record.coordinates, coordinates_before);
    }
}

#[cfg(test)]
mod source_parser_atrop_tests {
    use super::parse_smiles_complete_source as parse_smiles;
    use super::*;
    use cosmolkit_model::{BondStereo, SourceAtomValenceFacts};
    use cosmolkit_types::Hybridization;
    fn params() -> SmilesParseParams {
        SmilesParseParams {
            sanitize: false,
            remove_hs: false,
            ..Default::default()
        }
    }
    const TWO_D: &str = "(0,1,0;0,0,0;0,-1,0;1,0,0;1,1,0;1,-1,0)";
    const THREE_D: &str = "(0,1,0;0,0,0;0,-1,0;1,0,0;1,0,1;1,0,-1)";
    fn check_native_global_prelude(record: &SmilesRecord) {
        let cache = [(1, 3), (4, 0), (2, 0), (4, 0), (2, 0), (1, 3)];
        let hybs = [
            Hybridization::Sp3,
            Hybridization::Sp2,
            Hybridization::Sp2,
            Hybridization::Sp2,
            Hybridization::Sp2,
            Hybridization::Sp3,
        ];
        for (i, atom) in record.topology.atoms.iter().enumerate() {
            assert_eq!(
                atom.source_valence_facts(),
                SourceAtomValenceFacts {
                    explicit_valence: cache[i].0,
                    implicit_valence: cache[i].1
                }
            );
            assert_eq!(atom.hybridization(), hybs[i]);
        }
        assert_eq!(
            record
                .topology
                .bonds
                .iter()
                .map(|b| b.is_conjugated())
                .collect::<Vec<_>>(),
            vec![false, true, true, true, false]
        );
    }
    #[test]
    fn no_conformer_native_wedges_trigger_global_cache_conjugation_hybridization_and_atrop_tag() {
        let record = parse_smiles("CC(=O)C(=O)C |wU:1.0,3.4|", &params()).unwrap();
        check_native_global_prelude(&record);
        assert_eq!(record.topology.bonds[2].stereo(), BondStereo::AtropCcw);
        assert!(
            record.coordinates.conformers_2d.is_empty()
                && record.coordinates.conformers_3d.is_empty()
        );
        assert_eq!(record.properties.prop("_needsDetectAtomStereo"), None);
    }
    #[test]
    fn native_degree_failure_keeps_only_endpoint_cache_effects_and_skips_hybridization() {
        let record = parse_smiles("CCCC |wU:1.0,2.2|", &params()).unwrap();
        for (i, atom) in record.topology.atoms.iter().enumerate() {
            let expected = if i == 1 || i == 2 {
                SourceAtomValenceFacts {
                    explicit_valence: 2,
                    implicit_valence: 2,
                }
            } else {
                SourceAtomValenceFacts::UNINITIALIZED
            };
            assert_eq!(atom.source_valence_facts(), expected);
            assert_eq!(atom.hybridization(), Hybridization::Unspecified);
        }
        assert!(
            record
                .topology
                .bonds
                .iter()
                .all(|b| b.stereo() == BondStereo::None && !b.is_conjugated())
        );
    }
    #[test]
    fn source_atrop_uses_first_false_is3d_row_before_any_true_row_in_both_orders() {
        for blocks in [format!("{THREE_D}{TWO_D}"), format!("{TWO_D}{THREE_D}")] {
            let record =
                parse_smiles(&format!("CC(=O)C(=O)C |{blocks},wU:1.0,3.4|"), &params()).unwrap();
            check_native_global_prelude(&record);
            assert_eq!(record.topology.bonds[2].stereo(), BondStereo::AtropCcw);
            assert_eq!(
                record.coordinates.conformers_3d[0].is_3d(),
                blocks.starts_with(THREE_D)
            );
            assert_eq!(
                record.coordinates.conformers_3d[1].is_3d(),
                !blocks.starts_with(THREE_D)
            );
        }
        let record =
            parse_smiles(&format!("CC(=O)C(=O)C |{THREE_D},wU:1.0,3.4|"), &params()).unwrap();
        check_native_global_prelude(&record);
        assert_eq!(record.topology.bonds[2].stereo(), BondStereo::AtropCw);
        let small = "(0,0.0001,0;0,0,0;0,-0.0001,0;1,0,0;1,0,0.0001;1,0,-0.0001)";
        let record =
            parse_smiles(&format!("CC(=O)C(=O)C |{small},wU:1.0,3.4|"), &params()).unwrap();
        check_native_global_prelude(&record);
        assert_eq!(record.topology.bonds[2].stereo(), BondStereo::None);
    }
    #[test]
    fn native_normalization_failure_propagates_the_original_core_error_type() {
        let huge = "1.7976931348623157e308";
        let coords = format!("(0,1,0;0,0,0;0,-1,0;{huge},0,0;{huge},0,1;{huge},0,-1)");
        assert!(
            matches!(parse_smiles(&format!("CC(=O)C(=O)C |{coords},wU:1.0,3.4|"),&params()),
            Err(SmilesParseError::ParserAtropisomer(cosmolkit_core::AtropisomerError::Normalization(
                cosmolkit_core::StereoError::ZeroLengthVector {center,neighbor}
            ))) if center==AtomId::new(1) && neighbor==AtomId::new(3))
        );
    }
}

#[cfg(test)]
mod source_full_parser_composition_tests {
    use super::parse_smiles_complete_source as parse_smiles;
    use super::*;
    use cosmolkit_model::{Element, PropertyValue};

    #[test]
    fn source_flag_matrix_removes_hydrogens_and_finishes_stereo_at_the_actual_parser_boundary() {
        for sanitize in [false, true] {
            for remove_hs in [false, true] {
                let params = SmilesParseParams {
                    sanitize,
                    remove_hs,
                    ..Default::default()
                };
                let record = parse_smiles("[H]C sample", &params).unwrap();
                assert_eq!(record.topology.atoms.len(), if remove_hs { 1 } else { 2 });
                assert_eq!(record.properties.name(), Some(&"sample".into()));
                assert_eq!(
                    record.properties.prop("_StereochemDone"),
                    if sanitize || remove_hs {
                        Some(&PropertyValue::Int(1))
                    } else {
                        None
                    }
                );
                if sanitize || remove_hs {
                    assert!(
                        record
                            .properties
                            .is_prop_computed("_StereochemDone")
                            .unwrap()
                    );
                }
                record.topology.validate().unwrap();
                record
                    .coordinates
                    .validate_for_atom_count(record.topology.atoms.len())
                    .unwrap();
            }
        }
    }

    #[test]
    fn source_sanitize_failures_propagate_from_the_reached_canonical_owner() {
        let raw = SmilesParseParams {
            sanitize: false,
            remove_hs: false,
            ..Default::default()
        };
        assert_eq!(
            parse_smiles("C(C)(C)(C)(C)C", &raw)
                .unwrap()
                .topology
                .atoms
                .len(),
            6
        );
        let sanitize = SmilesParseParams {
            sanitize: true,
            remove_hs: false,
            ..Default::default()
        };
        assert!(matches!(
            parse_smiles("C(C)(C)(C)(C)C", &sanitize),
            Err(SmilesParseError::ParserSanitize(
                cosmolkit_core::SanitizeError::Properties { .. }
            ))
        ));
        assert!(matches!(
            parse_smiles("C(C)(C)(C)(C)C", &Default::default()),
            Err(SmilesParseError::ParserRemoveHydrogens(
                cosmolkit_core::HydrogenError::Sanitize(
                    cosmolkit_core::SanitizeError::Properties { .. }
                )
            ))
        ));
    }

    #[test]
    fn source_default_hydrogen_parameters_keep_isotopic_atoms() {
        let ordinary = parse_smiles("[H]C", &Default::default()).unwrap();
        let isotope = parse_smiles("[2H]C", &Default::default()).unwrap();
        assert_eq!(ordinary.topology.atoms.len(), 1);
        assert_eq!(ordinary.topology.atoms[0].element(), Element::C);
        assert_eq!(isotope.topology.atoms.len(), 2);
        assert_eq!(isotope.topology.atoms[0].element(), Element::H);
        assert_eq!(isotope.topology.atoms[0].isotope(), Some(2));
    }

    #[test]
    fn source_remove_hydrogens_remaps_cx_coordinates_before_stereo_finalization() {
        let record = parse_smiles("[H]C |(0,0,;1,0,)|", &Default::default()).unwrap();
        assert_eq!(record.topology.atoms.len(), 1);
        assert_eq!(record.coordinates.conformers_3d.len(), 1);
        assert_eq!(record.coordinates.conformers_3d[0].id(), 0);
        assert!(!record.coordinates.conformers_3d[0].is_3d());
        assert_eq!(
            record.coordinates.conformers_3d[0].coordinates(),
            &[[1.0, 0.0, 0.0]]
        );
    }

    #[test]
    fn source_aromatic_count_keeps_integral_type_on_both_sanitize_dispatches() {
        for remove_hs in [false, true] {
            let record = parse_smiles(
                "c1ccccc1",
                &SmilesParseParams {
                    remove_hs,
                    ..Default::default()
                },
            )
            .unwrap();
            assert_eq!(
                record.properties.prop("numArom"),
                Some(&PropertyValue::Int(1))
            );
            assert!(record.properties.is_prop_computed("numArom").unwrap());
            assert!(record.topology.atoms.iter().all(|atom| atom.is_aromatic()));
        }
    }
}

#[doc(hidden)]
pub use writer::canonicalize_fragment_source;

#[doc(hidden)]
pub use writer::canonicalize_fragment_from_bond_mask_source;

#[doc(hidden)]
pub use writer::canonicalize_query_fragment_source;

#[doc(hidden)]
pub use cx_writer::{write_cx_coordinates_from_source, zero_small_cx_coordinate};

#[doc(hidden)]
pub use cx_writer::quote_cx_atom_property;

#[doc(hidden)]
pub use cx_writer::write_query_cx_atom_properties_source;

#[doc(hidden)]
pub use cx_writer::write_cx_coord_or_hydrogen_bonds_source;

#[doc(hidden)]
pub use cx_writer::write_cx_zero_bonds_source;

#[doc(hidden)]
pub use cx_writer::{emit_cx_link_node_warning_source, write_query_cx_link_nodes_source};

#[doc(hidden)]
pub use cx_writer::append_cx_extension_source;

#[cfg(test)]
mod recovery_chem26 {
    use super::*;
    fn params() -> SmilesParseParams {
        SmilesParseParams {
            sanitize: false,
            remove_hs: false,
            skip_cleanup: true,
            ..Default::default()
        }
    }
    #[test]
    fn six_reductions_keep_only_reserved_ring_properties_and_map_every_slot() {
        for (text, slots, properties) in [
            ("CCC", vec![0, 1], vec![]),
            ("C-C", vec![0], vec![]),
            ("C=C", vec![0], vec![]),
            ("C(-C)C", vec![0, 1], vec![]),
            ("C(=C)C", vec![0, 1], vec![]),
            ("CC(C)C", vec![0, 1, 2], vec![]),
            ("C.C", vec![], vec![]),
            ("C1CC1.CC", vec![0, 1, 3, 2], vec![(3, 2)]),
            ("C1CC1CC", vec![0, 1, 4, 2, 3], vec![(4, 2)]),
            (
                "C1CC1C2CC2C",
                vec![0, 1, 6, 2, 3, 4, 7, 5],
                vec![(6, 2), (7, 6)],
            ),
        ] {
            let record = parse_smiles(text, &params()).unwrap();
            assert_eq!(record.topology.bonds.len(), slots.len(), "{text}");
            let actual = record
                .topology
                .bonds
                .iter()
                .filter_map(|b| {
                    b.prop(CXSMILES_BOND_IDX_PROP).map(|v| {
                        (
                            b.id().index(),
                            cosmolkit_core::property_value_to_uint(v).unwrap(),
                        )
                    })
                })
                .collect::<Vec<_>>();
            assert_eq!(actual, properties, "{text}");
            let _ = slots;
        }
    }
    #[test]
    fn all_cx_bond_families_resolve_ordinary_and_ring_parse_slots() {
        for (raw, physical, atom, other) in [(0, 0, 0, 1), (2, 4, 0, 2), (3, 2, 3, 2), (4, 3, 4, 3)]
        {
            for (kind, order) in [("C", BondOrder::Dative), ("H", BondOrder::Hydrogen)] {
                let r = parse_smiles(&format!("C1CC1CC |{kind}:{atom}.{raw}|"), &params()).unwrap();
                let b = &r.topology.bonds[physical];
                assert_eq!(b.order(), order);
                assert_eq!(
                    (b.begin(), b.end()),
                    (AtomId::new(atom), AtomId::new(other))
                );
            }
            let r = parse_smiles(&format!("C1CC1CC |Z:{raw}|"), &params()).unwrap();
            assert_eq!(r.topology.bonds[physical].order(), BondOrder::Zero);
            for (kind, cfg) in [("w", 2), ("wU", 1), ("wD", 3)] {
                let r = parse_smiles(&format!("C1CC1CC |{kind}:{atom}.{raw}|"), &params()).unwrap();
                assert_eq!(
                    r.topology.bonds[physical].prop("_MolFileBondCfg"),
                    Some(&PropertyValue::UInt(cfg))
                );
            }
        }
        for (text, raw, physical) in [("C1=CC1CC", 0, 0), ("C=1CC1CC", 2, 4), ("C1CC1=CC", 3, 2)] {
            for (kind, expected) in [
                ("c", cosmolkit_types::BondStereo::Cis),
                ("t", cosmolkit_types::BondStereo::Trans),
                ("ctu", cosmolkit_types::BondStereo::Any),
            ] {
                // Real shared lowerer before the separate parser stereo finalizer.
                let mut r = parse_smiles(text, &params()).unwrap();
                let cx = cosmolkit_cx::parse_cx_extensions(format!("|{kind}:{raw}|")).unwrap();
                cx_lowering::apply_cx_to_smiles_record(&mut r, &cx).unwrap();
                assert_eq!(r.topology.bonds[physical].stereo(), expected);
            }
        }
    }
    #[test]
    fn outer_raw_bond_guard_stays_before_fallback_and_endpoint_checks() {
        let mut r = parse_smiles("C1CC1CC", &params()).unwrap();
        let before = r.clone();
        for text in [
            "|C:0.5|",
            "|H:0.5|",
            "|Z:5|",
            "|wU:0.5|",
            "|c:5|",
            "|C:99.0|",
            "|wD:99.0|",
        ] {
            let cx = cosmolkit_cx::parse_cx_extensions(text).unwrap();
            cx_lowering::apply_cx_to_smiles_record(&mut r, &cx).unwrap();
            assert_eq!(r, before, "{text}");
        }
        let cx = cosmolkit_cx::parse_cx_extensions("|Z:5,3|").unwrap();
        cx_lowering::apply_cx_to_smiles_record(&mut r, &cx).unwrap();
        assert_eq!(r.topology.bonds[2].order(), BondOrder::Zero);
        assert_eq!(r.topology.bonds[4].order(), BondOrder::Single);
    }
}
