//! Detached SMILES parser and writer.
//!
//! The parser owns SMILES syntax and lowers directly into model blocks.  It
//! does not construct a live `Molecule`; the facade is responsible for
//! installing the returned values into its runtime state.

use std::collections::{BTreeMap, HashMap, HashSet};

use cosmolkit_cx::parse_cx_extensions;
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondDirection, BondId, BondSpec, CoordinateBlock,
    MoleculeProperties, TopologyBlock,
};
use cosmolkit_types::{BondOrder, ChiralTag, Element};

mod canonical_rank;
mod cx_lowering;
mod cx_writer;
mod stereo;
mod writer;

pub use cx_writer::{
    CxSmilesFields, CxSmilesWriteParams, write_cx_smiles, write_cx_smiles_with_params,
};
pub use writer::{SmilesWriteParams, write_smiles, write_smiles_with_params};

const CXSMILES_BOND_IDX_PROP: &str = "_cxsmilesBondIdx";

#[derive(Debug, Clone, PartialEq)]
pub struct SmilesRecord {
    pub topology: TopologyBlock,
    pub coordinates: CoordinateBlock,
    pub properties: MoleculeProperties,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SmilesParseError {
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
    #[error("SMILES writer stereochemistry preparation failed: {0}")]
    WriterStereo(String),
    #[error("invalid detached model: {0}")]
    Model(String),
    #[error("SMILES replacement key must not be empty")]
    EmptyReplacementKey,
    #[error("SMILES replacements do not converge; cyclic key '{key}' remains active")]
    ReplacementCycle { key: String },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmilesParseParams {
    pub sanitize: bool,
    pub allow_cxsmiles: bool,
    pub strict_cxsmiles: bool,
    pub parse_name: bool,
    pub remove_hydrogens: bool,
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
            remove_hydrogens: true,
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
        value = value
            .checked_mul(10)
            .and_then(|value| value.checked_add(u32::from(digit - b'0')))
            .filter(|value| *value <= i32::MAX as u32)
            .ok_or_else(|| atom_error(text, offset))?;
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
        let aromatic = matches!(token, "b" | "c" | "n" | "o" | "p" | "s" | "as" | "se");
        let element = match token {
            "b" => Some(Element::B),
            "c" => Some(Element::C),
            "n" => Some(Element::N),
            "o" => Some(Element::O),
            "p" => Some(Element::P),
            "s" => Some(Element::S),
            "as" => Some(Element::AS),
            "se" => Some(Element::SE),
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
    spec.with_direction(direction)
}

fn bond_order(symbol: char) -> Result<BondOrder, SmilesParseError> {
    match symbol {
        '-' => Ok(BondOrder::Single),
        '=' => Ok(BondOrder::Double),
        '#' => Ok(BondOrder::Triple),
        '$' => Ok(BondOrder::Quadruple),
        ':' => Ok(BondOrder::Aromatic),
        '~' => Ok(BondOrder::Unspecified),
        _ => Err(SmilesParseError::Unsupported {
            token: symbol,
            offset: 0,
        }),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RingPartial {
    atom: AtomId,
    order: Option<BondOrder>,
    direction: BondDirection,
    offset: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PendingRingClosure {
    ring: u32,
    opening: RingPartial,
    closing: RingPartial,
    cx_bond_index: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RingClosureRecord {
    ring: u32,
    bond: Option<BondId>,
}

fn push_smiles_bond(bonds: &mut Vec<Bond>, spec: BondSpec, cx_bond_index: usize) {
    // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy _cxsmilesBondIdx assignment
    // RDKit✔️✔️: res->setProp("_cxsmilesBondIdx", numBondsParsed++);
    // END RDKIT CPP GRAMMAR ACTION smiles.yy _cxsmilesBondIdx assignment
    let index = bonds.len();
    bonds.push(Bond::from_spec(
        BondId::new(index),
        spec.with_prop(CXSMILES_BOND_IDX_PROP, cx_bond_index.to_string())
            .expect("the internal CXSMILES bond-index property key is non-empty"),
    ));
}

fn opposite_bond_direction(direction: BondDirection) -> BondDirection {
    match direction {
        BondDirection::EndDownRight => BondDirection::EndUpRight,
        BondDirection::EndUpRight => BondDirection::EndDownRight,
        other => other,
    }
}

fn merged_ring_direction(target: RingPartial, source: RingPartial) -> BondDirection {
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
            (closure.opening, closure.closing)
        } else {
            (closure.closing, closure.opening)
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
        let direction = merged_ring_direction(selected, other);
        let bond = BondId::new(bonds.len());
        push_smiles_bond(
            bonds,
            BondSpec::new(begin, end, order)
                .with_aromatic(order == BondOrder::Aromatic)
                .with_direction(direction),
            closure.cx_bond_index,
        );
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

fn get_bond_ordering(
    atom: AtomId,
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    ring_records: &[RingClosureRecord],
) -> Result<(Vec<BondId>, usize), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION GetBondOrdering
    // RDKit✔️✔️: unsigned int GetBondOrdering(INT_LIST &bondOrdering, const RDKit::RWMol *mol,
    // RDKit✔️✔️:                              const RDKit::Atom *atom) {
    // RDKit✔️✔️:   INT_VECT ringClosures;
    // RDKit✔️✔️:   atom->getPropIfPresent(common_properties::_RingClosures, ringClosures);
    // RDKit✔️✔️:   std::list<SIZET_PAIR> neighbors;
    // RDKit✔️✔️:   neighbors.emplace_back(atom->getIdx(), -1);
    // RDKit✔️✔️:   for (auto nbrIdx : boost::make_iterator_range(mol->getAtomNeighbors(atom))) {
    // RDKit✔️✔️:     const Bond *nbrBond = mol->getBondBetweenAtoms(atom->getIdx(), nbrIdx);
    // RDKit✔️✔️:     if (std::find(ringClosures.begin(), ringClosures.end(),
    // RDKit✔️✔️:                   static_cast<int>(nbrBond->getIdx())) == ringClosures.end()) {
    // RDKit✔️✔️:       neighbors.emplace_back(nbrIdx, nbrBond->getIdx());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   neighbors.sort();
    // RDKit✔️✔️:   auto selfPos = neighbors.begin();
    // RDKit✔️✔️:   if (selfPos->first != atom->getIdx()) {
    // RDKit✔️✔️:     ++selfPos;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   CHECK_INVARIANT(selfPos->first == atom->getIdx(), "weird atom ordering");
    // RDKit✔️✔️:   for (auto neighborIt = neighbors.begin(); neighborIt != neighbors.end();
    // RDKit✔️✔️:        ++neighborIt) {
    // RDKit✔️✔️:     if (neighborIt != selfPos) {
    // RDKit✔️✔️:       bondOrdering.push_back(rdcast<int>(neighborIt->second));
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       bondOrdering.insert(bondOrdering.end(), ringClosures.begin(),
    // RDKit✔️✔️:                           ringClosures.end());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return ringClosures.size();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION GetBondOrdering
    let ring_closures = ring_records
        .iter()
        .map(|record| {
            record.bond.ok_or_else(|| {
                SmilesParseError::Model(format!(
                    "ring closure {} remained unresolved during chirality adjustment",
                    record.ring
                ))
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut neighbors = vec![(atom.index(), None)];
    for neighbor in adjacency.neighbors_of(atom.index()) {
        if !ring_closures.contains(&neighbor.bond) {
            neighbors.push((neighbor.atom_index, Some(neighbor.bond)));
        }
    }
    neighbors.sort_by_key(|(neighbor, _)| *neighbor);
    let self_position = neighbors
        .iter()
        .position(|(neighbor, bond)| *neighbor == atom.index() && bond.is_none())
        .ok_or_else(|| {
            SmilesParseError::Model("SMILES atom is absent from bond ordering".into())
        })?;
    let mut ordering = Vec::with_capacity(adjacency.neighbors_of(atom.index()).len());
    for (position, (_, bond)) in neighbors.into_iter().enumerate() {
        if position == self_position {
            ordering.extend(ring_closures.iter().copied());
        } else {
            ordering.push(bond.expect("only the self sentinel has no bond"));
        }
    }
    // Validate that every generated ID still refers to an incident bond. This
    // turns the source invariant into a structured detached-model error.
    if ordering.iter().any(|bond| {
        let bond = &bonds[bond.index()];
        bond.begin() != atom && bond.end() != atom
    }) {
        return Err(SmilesParseError::Model(
            "SMILES bond ordering contains a non-incident bond".into(),
        ));
    }
    Ok((ordering, ring_closures.len()))
}

fn adjust_atom_chirality_flags(
    atoms: &mut [Atom],
    bonds: &[Bond],
    adjacency: &AdjacencyList,
    ring_closures_by_atom: &[Vec<RingClosureRecord>],
    smiles_start_atoms: &[bool],
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION AdjustAtomChiralityFlags
    // RDKit✔️✔️: void AdjustAtomChiralityFlags(RWMol *mol) {
    // RDKit✔️✔️:   for (auto atom : mol->atoms()) {
    // RDKit✔️✔️:     Atom::ChiralType chiralType = atom->getChiralTag();
    // RDKit✔️✔️:     if (chiralType == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit✔️✔️:         chiralType == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit✔️✔️:       INT_LIST bondOrdering;
    // RDKit✔️✔️:       unsigned int numClosures = GetBondOrdering(bondOrdering, mol, atom);
    // RDKit✔️✔️:       int nSwaps = atom->getPerturbationOrder(bondOrdering);
    // RDKit✔️✔️:       if (Canon::chiralAtomNeedsTagInversion(
    // RDKit✔️✔️:               *mol, atom, atom->hasProp(common_properties::_SmilesStart),
    // RDKit✔️✔️:               numClosures)) {
    // RDKit✔️✔️:         ++nSwaps;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (nSwaps % 2) {
    // RDKit✔️✔️:         atom->invertChirality();
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (chiralType == Atom::CHI_SQUAREPLANAR ||
    // RDKit✔️✔️:                chiralType == Atom::CHI_TRIGONALBIPYRAMIDAL ||
    // RDKit✔️✔️:                chiralType == Atom::CHI_OCTAHEDRAL) {
    // RDKit✔️✔️:       INT_LIST bonds;
    // RDKit✔️✔️:       GetBondOrdering(bonds, mol, atom);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       unsigned int ref_max = Chirality::getMaxNbors(chiralType);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // insert (-1) for hydrogens or missing ligands, where these are placed
    // RDKit✔️✔️:       // depends on if it is the first atom or not
    // RDKit✔️✔️:       if (bonds.size() < ref_max) {
    // RDKit✔️✔️:         if (atom->hasProp(common_properties::_SmilesStart)) {
    // RDKit✔️✔️:           bonds.insert(bonds.begin(), ref_max - bonds.size(), -1);
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           bonds.insert(++bonds.begin(), ref_max - bonds.size(), -1);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       atom->setProp(common_properties::_chiralPermutation,
    // RDKit✔️✔️:                     Chirality::getChiralPermutation(atom, bonds, true));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION AdjustAtomChiralityFlags
    let mut invert = vec![false; atoms.len()];
    let mut nontetrahedral_permutations = vec![None; atoms.len()];
    for atom_index in 0..atoms.len() {
        let atom = &atoms[atom_index];
        let atom_id = AtomId::new(atom_index);
        match atom.chiral_tag() {
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw => {
                let (ordering, num_closures) = get_bond_ordering(
                    atom_id,
                    bonds,
                    adjacency,
                    &ring_closures_by_atom[atom_index],
                )?;
                let storage_order = adjacency
                    .neighbors_of(atom_index)
                    .iter()
                    .map(|neighbor| neighbor.bond)
                    .collect::<Vec<_>>();
                let mut swaps = stereo::count_swaps_to_interconvert(&ordering, storage_order)
                    .ok_or_else(|| {
                        SmilesParseError::Model(
                            "SMILES and storage bond orderings are not permutations".into(),
                        )
                    })?;
                let unsaturated = adjacency.neighbors_of(atom_index).iter().any(|neighbor| {
                    stereo::bond_order_as_double(bonds[neighbor.bond.index()].order()) > 1.0
                });
                if stereo::chiral_atom_needs_tag_inversion(
                    adjacency.neighbors_of(atom_index).len(),
                    atom.explicit_hydrogens(),
                    smiles_start_atoms[atom_index],
                    stereo::atom_has_fourth_valence(atom.explicit_hydrogens(), false),
                    num_closures,
                    unsaturated,
                ) {
                    swaps += 1;
                }
                invert[atom_index] = swaps % 2 == 1;
            }
            ChiralTag::SquarePlanar | ChiralTag::TrigonalBipyramidal | ChiralTag::Octahedral => {
                let (ordering, _) = get_bond_ordering(
                    atom_id,
                    bonds,
                    adjacency,
                    &ring_closures_by_atom[atom_index],
                )?;
                let mut probe = ordering.into_iter().map(Some).collect::<Vec<_>>();
                stereo::insert_implicit_nontetrahedral_neighbors(
                    &mut probe,
                    atom.chiral_tag(),
                    smiles_start_atoms[atom_index],
                );
                let incident = adjacency
                    .neighbors_of(atom_index)
                    .iter()
                    .map(|neighbor| neighbor.bond)
                    .collect::<Vec<_>>();
                nontetrahedral_permutations[atom_index] = Some(
                    stereo::nontetrahedral_chiral_permutation(
                        atom.chiral_permutation().unwrap_or(0),
                        atom.chiral_tag(),
                        bonds.len(),
                        &incident,
                        &probe,
                        true,
                    )
                    .map_err(|error| SmilesParseError::Model(error.to_string()))?,
                );
            }
            _ => {}
        }
    }
    for ((atom, invert), permutation) in atoms
        .iter_mut()
        .zip(invert)
        .zip(nontetrahedral_permutations)
    {
        if invert {
            atom.set_chiral_tag(stereo::invert_tetrahedral_tag(atom.chiral_tag()));
        }
        if let Some(permutation) = permutation {
            atom.set_chiral_permutation(Some(permutation));
        }
    }
    Ok(())
}

fn cleanup_after_parsing(record: &mut SmilesRecord) {
    // BEGIN RDKIT CPP FUNCTION CleanupAfterParsing
    // RDKit✔️✔️: for (auto atom : mol->atoms()) {
    // RDKit✔️✔️:   atom->clearProp(common_properties::_RingClosures);
    // RDKit✔️✔️:   atom->clearProp(common_properties::_SmilesStart);
    // RDKit✔️✔️:   std::string label;
    // RDKit✔️✔️:   if (atom->getAtomicNum() == 0 &&
    // RDKit✔️✔️:       atom->getPropIfPresent(common_properties::atomLabel, label)) {
    // RDKit✔️✔️:     if (label == "_AP1") {
    // RDKit✔️✔️:       atom->setProp(common_properties::_fromAttachPoint, 1);
    // RDKit✔️✔️:     } else if (label == "_AP2") {
    // RDKit✔️✔️:       atom->setProp(common_properties::_fromAttachPoint, 2);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: for (auto bond : mol->bonds()) {
    // RDKit✔️✔️:   bond->clearProp(common_properties::_unspecifiedOrder);
    // RDKit✔️✔️:   bond->clearProp("_cxsmilesBondIdx");
    // RDKit✔️✔️: }
    // RDKit✔️✔️: for (auto sg : RDKit::getSubstanceGroups(*mol)) {
    // RDKit✔️✔️:   sg.clearProp("_cxsmilesindex");
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION CleanupAfterParsing
    for atom in &mut record.topology.atoms {
        atom.clear_prop("_RingClosures");
        atom.clear_prop("_SmilesStart");
        if atom.atomic_number() == 0 {
            match atom.prop("atomLabel") {
                Some("_AP1") => atom
                    .set_prop("_fromAttachPoint", "1")
                    .expect("the internal attachment-point property key is non-empty"),
                Some("_AP2") => atom
                    .set_prop("_fromAttachPoint", "2")
                    .expect("the internal attachment-point property key is non-empty"),
                _ => {}
            }
        }
    }
    for bond in &mut record.topology.bonds {
        bond.clear_prop("_unspecifiedOrder");
        bond.clear_prop(CXSMILES_BOND_IDX_PROP);
    }
    for group in &mut record.topology.substance_groups {
        group.clear_prop("_cxsmilesindex");
    }
}

/// Parse SMILES into detached topology, coordinate, and property blocks.
pub fn parse_smiles(
    input: &str,
    params: &SmilesParseParams,
) -> Result<SmilesRecord, SmilesParseError> {
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
    let mut current = None::<AtomId>;
    let mut pending = None;
    let mut pending_direction = BondDirection::None;
    let mut next_cx_bond_index = 0;
    let mut index = 0;
    while index < bytes.len() {
        // BEGIN RDKIT CPP LEXER RULES smiles.ll dative bonds
        // RDKit✔️✔️: \-\>  { yylval->bond = new Bond(Bond::DATIVER);
        // RDKit✔️✔️:       return BOND_TOKEN; }
        // RDKit✔️✔️: \<\-  { yylval->bond = new Bond(Bond::DATIVEL);
        // RDKit✔️✔️:       return BOND_TOKEN; }
        // END RDKIT CPP LEXER RULES smiles.ll dative bonds
        if graph_text[index..].starts_with("->") {
            pending = Some(BondOrder::DativeRight);
            index += 2;
            continue;
        }
        if graph_text[index..].starts_with("<-") {
            pending = Some(BondOrder::DativeLeft);
            index += 2;
            continue;
        }
        match bytes[index] as char {
            '-' | '=' | '#' | ':' | '~' | '$' => {
                pending = Some(bond_order(bytes[index] as char)?);
                index += 1;
            }
            '/' => {
                pending_direction = BondDirection::EndUpRight;
                index += 1;
            }
            '\\' => {
                pending_direction = BondDirection::EndDownRight;
                index += 1;
            }
            '(' => {
                branches.push(current.ok_or_else(|| SmilesParseError::Syntax {
                    offset: index,
                    message: "branch has no preceding atom".into(),
                })?);
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
                    direction: pending_direction,
                    offset: ring_offset,
                };
                if let Some(opening) = rings.remove(&ring) {
                    ring_closures.push(PendingRingClosure {
                        ring,
                        opening,
                        closing: partial,
                        cx_bond_index: next_cx_bond_index,
                    });
                    next_cx_bond_index += 1;
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
                    let spec =
                        resolved_bond_spec(pending, &atoms, previous, atom, pending_direction);
                    push_smiles_bond(&mut bonds, spec, next_cx_bond_index);
                    degrees[previous.index()] += 1;
                    degrees[atom.index()] += 1;
                    next_cx_bond_index += 1;
                }
                current = Some(atom);
                pending = None;
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
                    let spec =
                        resolved_bond_spec(pending, &atoms, previous, atom, pending_direction);
                    push_smiles_bond(&mut bonds, spec, next_cx_bond_index);
                    degrees[previous.index()] += 1;
                    degrees[atom.index()] += 1;
                    next_cx_bond_index += 1;
                }
                current = Some(atom);
                pending = None;
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
    topology
        .validate()
        .map_err(|error| SmilesParseError::Model(error.to_string()))?;
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
            match parse_cx_extensions(&cx) {
                Ok(parsed) => match cx_lowering::apply_cx_to_smiles_record(&mut record, &parsed) {
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
                },
                Err(error) if params.strict_cxsmiles => {
                    return Err(SmilesParseError::Cx(error.to_string()));
                }
                Err(_) => record
                    .properties
                    .set_prop("_CXSMILES_Data", "")
                    .expect("the internal CXSMILES data property key is non-empty"),
            }
        } else if params.allow_cxsmiles && params.strict_cxsmiles && !params.parse_name {
            return Err(SmilesParseError::Cx(
                "CXSMILES extension does not start with | and parseName=false".into(),
            ));
        } else if params.parse_name {
            name = cx.trim().to_owned();
        }
    }
    if !name.is_empty() {
        record.properties = record.properties.with_name(&name);
    }
    // BEGIN RDKIT CPP FUNCTION MolFromSmiles (cleanup gate)
    // RDKit✔️✔️:   if (res) {
    // RDKit✔️✔️:     if (!params.skipCleanup) {
    // RDKit✔️✔️:       SmilesParseOps::CleanupAfterParsing(res.get());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!name.empty()) {
    // RDKit✔️✔️:       res->setProp(common_properties::_Name, name);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION MolFromSmiles (cleanup gate)
    if !params.skip_cleanup {
        cleanup_after_parsing(&mut record);
    }
    record
        .topology
        .validate()
        .map_err(|error| SmilesParseError::Model(error.to_string()))?;
    record
        .coordinates
        .validate_for_atom_count(record.topology.atoms.len())
        .map_err(|error| SmilesParseError::Model(error.to_string()))?;
    Ok(record)
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_model::{BondStereo, StereoGroupKind, SubstanceGroupKind};

    #[test]
    fn parses_and_writes_detached_ethanol() {
        let record = parse_smiles("CCO", &Default::default()).expect("parse");
        assert_eq!(record.topology.atoms.len(), 3);
        assert_eq!(write_smiles(&record).expect("write"), "CCO");
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
        // column is therefore a 2D conformer, not 3D.
        assert_eq!(record.coordinates.conformers_3d.len(), 0);
        assert_eq!(record.coordinates.conformers_2d.len(), 1);
        assert_eq!(record.properties.name(), Some("ethanol"));
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
        assert_eq!(plain.properties.name(), Some("ethanol sample"));

        let cx = parse_smiles("CC |$foo;bar$| named sample", &Default::default()).unwrap();
        assert_eq!(cx.properties.name(), Some("named sample"));
        assert_eq!(cx.topology.atoms[0].prop("atomLabel"), Some("foo"));

        let no_cx = parse_smiles(
            "CC |not cx when disabled|",
            &SmilesParseParams {
                allow_cxsmiles: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(no_cx.properties.name(), Some("|not cx when disabled|"));
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
        assert_eq!(wedge.topology.bonds[0].prop("_MolFileBondCfg"), Some("1"));
        assert_eq!(wedge.properties.prop("_needsDetectAtomStereo"), None);
        assert_eq!(wedge.topology.atoms[1].chiral_tag(), ChiralTag::Unspecified);

        let double = parse_smiles("CC=CC |c:1|", &Default::default()).unwrap();
        assert_eq!(double.topology.bonds[1].stereo(), BondStereo::Cis);
        assert_eq!(
            double.topology.bonds[1].stereo_atoms(),
            Some([AtomId::new(0), AtomId::new(3)])
        );
        assert_eq!(double.properties.prop("_needsDetectBondStereo"), Some("1"));
    }

    #[test]
    fn lowers_cx_linknodes_sgroups_and_variable_attachments_like_rdkit() {
        let link = parse_smiles("C1CC1 |LN:1:1.3|", &Default::default()).unwrap();
        assert_eq!(
            link.properties.prop("_MolFileLinkNodes"),
            Some("1 3 2 2 1 2 3")
        );

        let data = parse_smiles("CCO |SgD:2,1:FIELD:info::::|", &Default::default()).unwrap();
        assert_eq!(data.topology.substance_groups.len(), 1);
        let data_group = &data.topology.substance_groups[0];
        assert_eq!(data_group.kind(), &SubstanceGroupKind::Data);
        assert_eq!(data_group.atoms(), &[AtomId::new(2), AtomId::new(1)]);
        assert_eq!(
            data_group.props().get("FIELDNAME").map(String::as_str),
            Some("FIELD")
        );
        assert_eq!(data_group.data_fields(), &["info".to_owned()]);

        let polymer = parse_smiles("CC |Sg:n:0::ht|", &Default::default()).unwrap();
        assert_eq!(polymer.topology.substance_groups.len(), 1);
        assert_eq!(
            polymer.topology.substance_groups[0].kind(),
            &SubstanceGroupKind::StructuralRepeatUnit
        );

        let attachment = parse_smiles("CO*.C1=CC=NC=C1 |m:2:3.5.4|", &Default::default()).unwrap();
        assert_eq!(
            attachment.topology.bonds[1].prop("_MolFileBondEndPts"),
            Some("(3 4 6 5)")
        );
        assert_eq!(
            attachment.topology.bonds[1].prop("_MolFileBondAttach"),
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
        assert_eq!(non_strict.properties.prop("_CXSMILES_Data"), Some(""));
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
        assert!(
            record
                .topology
                .substance_groups
                .iter()
                .all(|group| { group.props().get("_cxsmilesindex").is_none() })
        );
    }
}
