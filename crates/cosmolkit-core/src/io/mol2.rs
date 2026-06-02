//! RDKit-compatible MOL2 reader types.
//!
//! Parser behavior is ported step by step from RDKit
//! `Code/GraphMol/FileParsers/Mol2FileParser.cpp` and declarations in
//! `Code/GraphMol/FileParsers/FileParsers.h`.

use std::{fs, path::Path};

use crate::{
    AdjacencyList, AtomId, AtomQueryPredicate, AtomSpec, BondOrder, BondSpec, Conformer3D, Element,
    Molecule, MoleculeBuilder, QueryNode, UnsupportedFeatureError, rdkit_valence_list,
    rings::find_sssr_from_parts,
    valence::{bond_valence_contrib, periodic_table_outer_electrons},
};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum Mol2ReadError {
    #[error(transparent)]
    UnsupportedFeature(#[from] UnsupportedFeatureError),
    #[error("MOL2 I/O failed: {0}")]
    Io(String),
    #[error("MOL2 parse failed: {0}")]
    Parse(String),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mol2Type {
    // BEGIN RDKIT CPP DECLARATION third_party/rdkit/Code/GraphMol/FileParsers/FileParsers.h :: Mol2Type
    // RDKit✔️✔️: typedef enum {
    // RDKit✔️✔️:   CORINA = 0  //!< supports output from Corina and some dbtranslate output
    // RDKit✔️✔️: } Mol2Type;
    // END RDKIT CPP DECLARATION third_party/rdkit/Code/GraphMol/FileParsers/FileParsers.h :: Mol2Type
    Corina,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Mol2ReadParams {
    pub sanitize: bool,
    pub remove_hs: bool,
    pub variant: Mol2Type,
    pub cleanup_substructures: bool,
}

impl Default for Mol2ReadParams {
    fn default() -> Self {
        // BEGIN RDKIT CPP DECLARATION third_party/rdkit/Code/GraphMol/FileParsers/FileParsers.h :: Mol2ParserParams
        // RDKit✔️✔️: struct Mol2ParserParams {
        // RDKit✔️✔️:   bool sanitize = true; /**< sanitize the molecule after building it */
        // RDKit✔️✔️:   bool removeHs = true; /**< remove Hs after constructing the molecule */
        // RDKit✔️✔️:   Mol2Type variant = Mol2Type::CORINA; /**< the atom type definitions to use */
        // RDKit✔️✔️:   bool cleanupSubstructures =
        // RDKit✔️✔️:       true; /**< toggles recognition and cleanup of common substructures */
        // RDKit✔️✔️: };
        // END RDKIT CPP DECLARATION third_party/rdkit/Code/GraphMol/FileParsers/FileParsers.h :: Mol2ParserParams
        Self {
            sanitize: true,
            remove_hs: true,
            variant: Mol2Type::Corina,
            cleanup_substructures: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Mol2Record {
    pub molecule: Molecule,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Mol2SectionOffsets {
    molecule_start: usize,
    atom_start: usize,
    bond_start: Option<usize>,
    charge_start: Option<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Mol2MoleculeHeader {
    name: String,
    n_atoms: u32,
    n_bonds: u32,
    charge_type: String,
}

#[derive(Debug, Clone, PartialEq)]
struct Mol2AtomLine {
    atom_spec: AtomSpec,
    atom_name: String,
    position: [f64; 3],
    sybyl_atom_type: String,
    sybyl_symbol: String,
    no_implicit: bool,
    partial_charge: Option<String>,
}

#[derive(Debug, Clone)]
struct Mol2AtomBlock {
    builder: MoleculeBuilder,
    idx_corresp: Vec<Option<AtomId>>,
    has_h_atoms: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Mol2BondLine {
    bond_spec: BondSpec,
    sybyl_bond_type: String,
    zero_based_begin_idx: u32,
    zero_based_end_idx: u32,
}

#[derive(Debug, Clone)]
struct Mol2BondBlock {
    builder: MoleculeBuilder,
    n_bad_bonds: u32,
}

#[derive(Debug, Clone)]
struct Mol2UnityAtomAttrBlock {
    builder: MoleculeBuilder,
    next_offset: usize,
}

fn rdkit_get_line_at_with_eof<'a>(input: &'a str, offset: &mut usize) -> (&'a str, bool) {
    // BEGIN RDKIT CPP HELPER third_party/rdkit/Code/RDGeneral/StreamOps.h :: getLine
    // RDKit✔️✔️: inline std::string getLine(std::istream *inStream) {
    // RDKit✔️✔️:   std::string res;
    // RDKit✔️✔️:   std::getline(*inStream, res);
    // RDKit✔️✔️:   if (!res.empty() && (res.back() == '\r')) {
    // RDKit✔️✔️:     res.resize(res.length() - 1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP HELPER third_party/rdkit/Code/RDGeneral/StreamOps.h :: getLine
    if *offset >= input.len() {
        return ("", true);
    }

    let remaining = &input[*offset..];
    let (line, next_offset, eof_after_read) = match remaining.find('\n') {
        Some(newline_idx) => (&remaining[..newline_idx], *offset + newline_idx + 1, false),
        None => (remaining, input.len(), true),
    };
    *offset = next_offset;
    (line.strip_suffix('\r').unwrap_or(line), eof_after_read)
}

fn rdkit_get_line_at<'a>(input: &'a str, offset: &mut usize) -> &'a str {
    rdkit_get_line_at_with_eof(input, offset).0
}

fn rdkit_get_line_at_checked<'a>(
    input: &'a str,
    offset: &mut usize,
) -> Result<&'a str, Mol2ReadError> {
    let (line, eof_after_read) = rdkit_get_line_at_with_eof(input, offset);
    if eof_after_read {
        return Err(Mol2ReadError::Parse("premature EOF".to_string()));
    }
    Ok(line)
}

fn parse_rdkit_unsigned_int_token(token: &str) -> Option<u32> {
    let (negative, digits) = match token.as_bytes().first() {
        Some(b'+') => (false, &token[1..]),
        Some(b'-') => (true, &token[1..]),
        _ => (false, token),
    };
    if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
        return None;
    }

    let value = digits.parse::<u64>().ok()?;
    if value > u32::MAX as u64 {
        return None;
    }
    Some(if negative {
        0u32.wrapping_sub(value as u32)
    } else {
        value as u32
    })
}

fn parse_mol2_molecule_header_like_rdkit(
    input: &str,
    molecule_start: usize,
) -> Result<Mol2MoleculeHeader, Mol2ReadError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2DataStream molecule header
    // RDKit✔️✔️: inStream.seekg(molStart, std::ios::beg);
    // RDKit✔️✔️: tempStr = getLine(inStream);
    // RDKit✔️✔️: auto res = std::make_unique<RWMol>();
    // RDKit✔️✔️: boost::trim_right(tempStr);
    // RDKit✔️✔️: res->setProp(common_properties::_Name, tempStr);
    // RDKit✔️✔️:
    // RDKit✔️✔️: tempStr = getLine(inStream);
    // RDKit✔️✔️: tokenizer tokens(tempStr, sep);
    // RDKit✔️✔️: if (tokens.begin() == tokens.end()) {
    // RDKit✔️✔️:   throw FileParseException("Empty counts line");
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: unsigned int nAtoms = 0, nBonds = 0;
    // RDKit✔️✔️: tokenizer::const_iterator itemIt = tokens.begin();
    // RDKit✔️✔️: // counts line, this is where we really get started
    // RDKit✔️✔️: try {
    // RDKit✔️✔️:   nAtoms = boost::lexical_cast<unsigned int>(*itemIt);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   ++itemIt;
    // RDKit✔️✔️:   if (itemIt != tokens.end()) {
    // RDKit✔️✔️:     nBonds = boost::lexical_cast<unsigned int>(*itemIt);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: } catch (boost::bad_lexical_cast &) {
    // RDKit✔️✔️:   std::ostringstream errout;
    // RDKit✔️✔️:   errout << "Cannot convert " << *itemIt << " to unsigned int";
    // RDKit✔️✔️:   throw FileParseException(errout.str());
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: if (nAtoms == 0) {
    // RDKit✔️✔️:   throw FileParseException("molecule has no atoms");
    // RDKit✔️✔️: }
    // RDKit✔️✔️: tempStr = getLine(inStream);  // mol_type - ignore
    // RDKit✔️✔️: tempStr = getLine(inStream);
    // RDKit✔️✔️: boost::trim(tempStr);
    // RDKit✔️✔️: res->setProp("_TriposChargeType", tempStr);
    // RDKit✔️✔️: // stop here since we don't support anything else from the MOLECULE block
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2DataStream molecule header
    let mut offset = molecule_start;
    let name = rdkit_get_line_at(input, &mut offset)
        .trim_end_matches(|ch: char| ch.is_ascii_whitespace())
        .to_string();

    let counts_line = rdkit_get_line_at(input, &mut offset);
    let mut tokens = counts_line
        .split([' ', '\t', '\n'])
        .filter(|token| !token.is_empty());
    let atom_token = tokens
        .next()
        .ok_or_else(|| Mol2ReadError::Parse("Empty counts line".to_string()))?;
    let n_atoms = parse_rdkit_unsigned_int_token(atom_token).ok_or_else(|| {
        Mol2ReadError::Parse(format!("Cannot convert {atom_token} to unsigned int"))
    })?;

    let n_bonds = if let Some(bond_token) = tokens.next() {
        parse_rdkit_unsigned_int_token(bond_token).ok_or_else(|| {
            Mol2ReadError::Parse(format!("Cannot convert {bond_token} to unsigned int"))
        })?
    } else {
        0
    };

    if n_atoms == 0 {
        return Err(Mol2ReadError::Parse("molecule has no atoms".to_string()));
    }

    let _mol_type = rdkit_get_line_at(input, &mut offset);
    let charge_type = rdkit_get_line_at(input, &mut offset)
        .trim_matches(|ch: char| ch.is_ascii_whitespace())
        .to_string();

    Ok(Mol2MoleculeHeader {
        name,
        n_atoms,
        n_bonds,
        charge_type,
    })
}

fn rdkit_atomic_number_for_symbol_like_rdkit(symbol: &str) -> Result<u8, Mol2ReadError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/PeriodicTable.h :: getAtomicNumber
    // RDKit✔️✔️: int getAtomicNumber(const std::string &elementSymbol) const {
    // RDKit✔️✔️:   // this little optimization actually makes a measurable difference
    // RDKit✔️✔️:   // in molecule-construction time
    // RDKit✔️✔️:   int anum = -1;
    // RDKit✔️✔️:   if (elementSymbol == "C") {
    // RDKit✔️✔️:     anum = 6;
    // RDKit✔️✔️:   } else if (elementSymbol == "N") {
    // RDKit✔️✔️:     anum = 7;
    // RDKit✔️✔️:   } else if (elementSymbol == "O") {
    // RDKit✔️✔️:     anum = 8;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     STR_UINT_MAP::const_iterator iter = byname.find(elementSymbol);
    // RDKit✔️✔️:     if (iter != byname.end()) {
    // RDKit✔️✔️:       anum = iter->second;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   POSTCONDITION(anum > -1, "Element '" + elementSymbol + "' not found");
    // RDKit✔️✔️:   return anum;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/PeriodicTable.h :: getAtomicNumber
    let atomic_number = match symbol {
        "*" => 0,
        "H" => 1,
        "He" => 2,
        "Li" => 3,
        "Be" => 4,
        "B" => 5,
        "C" => 6,
        "N" => 7,
        "O" => 8,
        "F" => 9,
        "Ne" => 10,
        "Na" => 11,
        "Mg" => 12,
        "Al" => 13,
        "Si" => 14,
        "P" => 15,
        "S" => 16,
        "Cl" => 17,
        "Ar" => 18,
        "K" => 19,
        "Ca" => 20,
        "Sc" => 21,
        "Ti" => 22,
        "V" => 23,
        "Cr" => 24,
        "Mn" => 25,
        "Fe" => 26,
        "Co" => 27,
        "Ni" => 28,
        "Cu" => 29,
        "Zn" => 30,
        "Ga" => 31,
        "Ge" => 32,
        "As" => 33,
        "Se" => 34,
        "Br" => 35,
        "Kr" => 36,
        "Rb" => 37,
        "Sr" => 38,
        "Y" => 39,
        "Zr" => 40,
        "Nb" => 41,
        "Mo" => 42,
        "Tc" => 43,
        "Ru" => 44,
        "Rh" => 45,
        "Pd" => 46,
        "Ag" => 47,
        "Cd" => 48,
        "In" => 49,
        "Sn" => 50,
        "Sb" => 51,
        "Te" => 52,
        "I" => 53,
        "Xe" => 54,
        "Cs" => 55,
        "Ba" => 56,
        "La" => 57,
        "Ce" => 58,
        "Pr" => 59,
        "Nd" => 60,
        "Pm" => 61,
        "Sm" => 62,
        "Eu" => 63,
        "Gd" => 64,
        "Tb" => 65,
        "Dy" => 66,
        "Ho" => 67,
        "Er" => 68,
        "Tm" => 69,
        "Yb" => 70,
        "Lu" => 71,
        "Hf" => 72,
        "Ta" => 73,
        "W" => 74,
        "Re" => 75,
        "Os" => 76,
        "Ir" => 77,
        "Pt" => 78,
        "Au" => 79,
        "Hg" => 80,
        "Tl" => 81,
        "Pb" => 82,
        "Bi" => 83,
        "Po" => 84,
        "At" => 85,
        "Rn" => 86,
        "Fr" => 87,
        "Ra" => 88,
        "Ac" => 89,
        "Th" => 90,
        "Pa" => 91,
        "U" => 92,
        "Np" => 93,
        "Pu" => 94,
        "Am" => 95,
        "Cm" => 96,
        "Bk" => 97,
        "Cf" => 98,
        "Es" => 99,
        "Fm" => 100,
        "Md" => 101,
        "No" => 102,
        "Lr" => 103,
        "Rf" => 104,
        "Db" => 105,
        "Sg" => 106,
        "Bh" => 107,
        "Hs" => 108,
        "Mt" => 109,
        "Ds" => 110,
        "Rg" => 111,
        "Cn" => 112,
        "Nh" | "Uut" => 113,
        "Fl" => 114,
        "Mc" | "Uup" => 115,
        "Lv" => 116,
        "Ts" => 117,
        "Og" => 118,
        _ => {
            return Err(Mol2ReadError::Parse(format!(
                "Element '{symbol}' not found"
            )));
        }
    };
    Ok(atomic_number)
}

fn mol2_atom_spec_from_sybyl_symbol_like_rdkit(
    symbol: &str,
) -> Result<Option<AtomSpec>, Mol2ReadError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: ParseMol2FileAtomLine atom type handling
    // RDKit✔️✔️:   // bad symbols:
    // RDKit✔️✔️:   // LP is not an atom so remove it ...
    // RDKit✔️✔️:   if (symb == "LP") {
    // RDKit✔️✔️:     delete res;
    // RDKit✔️✔️:     return nullptr;
    // RDKit✔️✔️:   } else if (symb == "ANY" || symb == "Du") {
    // RDKit✔️✔️:     // queryAtoms
    // RDKit✔️✔️:     // according to the SYBYL spec, these match anything
    // RDKit✔️✔️:     auto *query = new QueryAtom(0);
    // RDKit✔️✔️:     query->setQuery(makeAtomNullQuery());
    // RDKit✔️✔️:     delete res;
    // RDKit✔️✔️:     res = query;
    // RDKit✔️✔️:   } else if (symb == "HEV") {
    // RDKit✔️✔️:     auto *query = new QueryAtom(1);
    // RDKit✔️✔️:     query->getQuery()->setNegation(true);
    // RDKit✔️✔️:     delete res;
    // RDKit✔️✔️:     res = query;
    // RDKit✔️✔️:   } else if (symb == "HET") {
    // RDKit✔️✔️:     // Tripos: N,O,P,S
    // RDKit✔️✔️:     auto *query = new QueryAtom(7);
    // RDKit✔️✔️:     query->expandQuery(makeAtomNumQuery(8), Queries::COMPOSITE_OR, true);
    // RDKit✔️✔️:     query->expandQuery(makeAtomNumQuery(15), Queries::COMPOSITE_OR, true);
    // RDKit✔️✔️:     query->expandQuery(makeAtomNumQuery(16), Queries::COMPOSITE_OR, true);
    // RDKit✔️✔️:     delete res;
    // RDKit✔️✔️:     res = query;
    // RDKit✔️✔️:   } else if (symb == "HAL") {
    // RDKit✔️✔️:     // Tripos: F,Cl,Br,I
    // RDKit✔️✔️:     auto *query = new QueryAtom(9);
    // RDKit✔️✔️:     query->expandQuery(makeAtomNumQuery(17), Queries::COMPOSITE_OR, true);
    // RDKit✔️✔️:     query->expandQuery(makeAtomNumQuery(35), Queries::COMPOSITE_OR, true);
    // RDKit✔️✔️:     query->expandQuery(makeAtomNumQuery(53), Queries::COMPOSITE_OR, true);
    // RDKit✔️✔️:     delete res;
    // RDKit✔️✔️:     res = query;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res->setAtomicNum(PeriodicTable::getTable()->getAtomicNumber(symb));
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: ParseMol2FileAtomLine atom type handling
    let spec = match symbol {
        "LP" => return Ok(None),
        "ANY" | "Du" => {
            AtomSpec::new(Element::DUMMY).with_query(QueryNode::predicate(AtomQueryPredicate::Any))
        }
        "HEV" => AtomSpec::new(Element::H).with_query(QueryNode::not(QueryNode::predicate(
            AtomQueryPredicate::AtomicNumber(1),
        ))),
        "HET" => AtomSpec::new(Element::N).with_query(QueryNode::predicate(
            AtomQueryPredicate::AtomicNumberIn(vec![7, 8, 15, 16]),
        )),
        "HAL" => AtomSpec::new(Element::F).with_query(QueryNode::predicate(
            AtomQueryPredicate::AtomicNumberIn(vec![9, 17, 35, 53]),
        )),
        _ => {
            let atomic_number = rdkit_atomic_number_for_symbol_like_rdkit(symbol)?;
            AtomSpec::new(Element::from_atomic_number(atomic_number).expect("u8 atomic number"))
        }
    };
    Ok(Some(spec))
}

fn parse_mol2_atom_line_like_rdkit(atom_line: &str) -> Result<Option<Mol2AtomLine>, Mol2ReadError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: ParseMol2FileAtomLine base fields
    // RDKit✔️✔️: Atom *ParseMol2FileAtomLine(const std::string atomLine, RDGeom::Point3D &pos) {
    // RDKit✔️✔️:   typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
    // RDKit✔️✔️:   boost::char_separator<char> sep(" \t\n");
    // RDKit✔️✔️:   std::string tAN, tAT;
    // RDKit✔️✔️:   tokenizer tokens(atomLine, sep);
    // RDKit✔️✔️:   tokenizer::const_iterator itemIt = tokens.begin();
    // RDKit✔️✔️:   if (itemIt == tokens.end()) {
    // RDKit✔️✔️:     throw FileParseException("no info in mol2 atom line");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto *res = new Atom();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // skip TriposAtomId
    // RDKit✔️✔️:   ++itemIt;
    // RDKit✔️✔️:   if (itemIt == tokens.end()) {
    // RDKit✔️✔️:     delete res;
    // RDKit✔️✔️:     throw FileParseException("premature end of mol2 atom line");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // the sybyl atom name does not necessarily make sense - into atom property
    // RDKit✔️✔️:   tAN = *itemIt;
    // RDKit✔️✔️:   ++itemIt;
    // RDKit✔️✔️:   if (itemIt == tokens.end()) {
    // RDKit✔️✔️:     delete res;
    // RDKit✔️✔️:     throw FileParseException("premature end of mol2 atom line");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     pos.x = boost::lexical_cast<double>(*itemIt);
    // RDKit✔️✔️:     ++itemIt;
    // RDKit✔️✔️:     if (itemIt == tokens.end()) {
    // RDKit✔️✔️:       delete res;
    // RDKit✔️✔️:       throw FileParseException("premature end of mol2 atom line");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     pos.y = boost::lexical_cast<double>(*itemIt);
    // RDKit✔️✔️:     ++itemIt;
    // RDKit✔️✔️:     if (itemIt == tokens.end()) {
    // RDKit✔️✔️:       delete res;
    // RDKit✔️✔️:       throw FileParseException("premature end of mol2 atom line");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     pos.z = boost::lexical_cast<double>(*itemIt);
    // RDKit✔️✔️:     ++itemIt;
    // RDKit✔️✔️:     if (itemIt == tokens.end()) {
    // RDKit✔️✔️:       delete res;
    // RDKit✔️✔️:       throw FileParseException("premature end of mol2 atom line");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️✔️:     delete res;
    // RDKit✔️✔️:     throw FileParseException("Cannot process mol2 coordinates.");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // now it becomes interesting - this is the SYBYL atom type. I put this into
    // RDKit✔️✔️:   // an atom property and deduce the
    // RDKit✔️✔️:   // symbol from that (everything before the dot)
    // RDKit✔️✔️:   tAT = *itemIt;
    // RDKit✔️✔️:   std::string symb = (*itemIt).substr(0, (*itemIt).find('.'));
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // now assign the properties
    // RDKit✔️✔️:   res->setProp("_TriposAtomName", tAN);  // maybe remove that since it's
    // RDKit✔️✔️:                                            // useless?
    // RDKit✔️✔️:   res->setProp(common_properties::_TriposAtomType, tAT);
    // RDKit✔️✔️:   // no implicit hydrogens for mol2 files
    // RDKit✔️✔️:   res->setNoImplicit(true);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // up to here the fields must be written - check if we have more
    // RDKit✔️✔️:   // next comes the subst_id, subst_name - skip those
    // RDKit✔️✔️:   ++itemIt;
    // RDKit✔️✔️:   if (itemIt == tokens.end()) {
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++itemIt;
    // RDKit✔️✔️:   if (itemIt == tokens.end()) {
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++itemIt;
    // RDKit✔️✔️:   // the Partial charge in the file
    // RDKit✔️✔️:   if (itemIt != tokens.end()) {
    // RDKit✔️✔️:     res->setProp("_TriposPartialCharge", *itemIt);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // we skip the status bit ...
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: ParseMol2FileAtomLine base fields
    let mut tokens = atom_line
        .split([' ', '\t', '\n'])
        .filter(|token| !token.is_empty());
    let _tripos_atom_id = tokens
        .next()
        .ok_or_else(|| Mol2ReadError::Parse("no info in mol2 atom line".to_string()))?;

    let atom_name = tokens
        .next()
        .ok_or_else(|| Mol2ReadError::Parse("premature end of mol2 atom line".to_string()))?
        .to_string();

    let x_token = tokens
        .next()
        .ok_or_else(|| Mol2ReadError::Parse("premature end of mol2 atom line".to_string()))?;
    let x = x_token
        .parse::<f64>()
        .map_err(|_| Mol2ReadError::Parse("Cannot process mol2 coordinates.".to_string()))?;
    let y_token = tokens
        .next()
        .ok_or_else(|| Mol2ReadError::Parse("premature end of mol2 atom line".to_string()))?;
    let y = y_token
        .parse::<f64>()
        .map_err(|_| Mol2ReadError::Parse("Cannot process mol2 coordinates.".to_string()))?;
    let z_token = tokens
        .next()
        .ok_or_else(|| Mol2ReadError::Parse("premature end of mol2 atom line".to_string()))?;
    let z = z_token
        .parse::<f64>()
        .map_err(|_| Mol2ReadError::Parse("Cannot process mol2 coordinates.".to_string()))?;

    let sybyl_atom_type = tokens
        .next()
        .ok_or_else(|| Mol2ReadError::Parse("premature end of mol2 atom line".to_string()))?
        .to_string();
    let sybyl_symbol = sybyl_atom_type
        .split_once('.')
        .map_or(sybyl_atom_type.as_str(), |(symbol, _)| symbol)
        .to_string();

    let _subst_id = tokens.next();
    let _subst_name = tokens.next();
    let partial_charge = tokens.next().map(str::to_string);

    let Some(atom_spec) = mol2_atom_spec_from_sybyl_symbol_like_rdkit(&sybyl_symbol)? else {
        return Ok(None);
    };
    let mut atom_spec = atom_spec
        .with_prop("_TriposAtomName", atom_name.clone())
        .with_prop("_TriposAtomType", sybyl_atom_type.clone())
        .with_no_implicit(true);
    if let Some(partial_charge) = &partial_charge {
        atom_spec = atom_spec.with_prop("_TriposPartialCharge", partial_charge.clone());
    }

    Ok(Some(Mol2AtomLine {
        atom_spec,
        atom_name,
        position: [x, y, z],
        sybyl_atom_type,
        sybyl_symbol,
        no_implicit: true,
        partial_charge,
    }))
}

fn parse_mol2_atom_block_like_rdkit(
    input: &str,
    atom_start: usize,
    n_atoms: u32,
) -> Result<Mol2AtomBlock, Mol2ReadError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: ParseMol2AtomBlock
    // RDKit✔️✔️: void ParseMol2AtomBlock(std::istream *inStream, RWMol *res, unsigned int nAtoms,
    // RDKit✔️✔️:                         INT_VECT &idxCorresp) {
    // RDKit✔️✔️:   PRECONDITION(inStream, "inStream not valid");
    // RDKit✔️✔️:   PRECONDITION(!inStream->eof(), "inStream is at eof");
    // RDKit✔️✔️:   PRECONDITION(res, "RWMol not valid");
    // RDKit✔️✔️:   PRECONDITION(idxCorresp.size() == nAtoms, "vector size mismatch");
    // RDKit✔️✔️:   unsigned int nLP = 0;
    // RDKit✔️✔️:   bool hasHAtoms = false;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<RDGeom::Point3D> threeDPs;
    // RDKit✔️✔️:   threeDPs.reserve(nAtoms);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit✔️✔️:     std::string tempStr = getLine(inStream);
    // RDKit✔️✔️:     if (inStream->eof()) {
    // RDKit✔️✔️:       throw FileParseException("premature EOF");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     RDGeom::Point3D pos;
    // RDKit✔️✔️:     Atom *atom = ParseMol2FileAtomLine(tempStr, pos);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // if atom is NULL then we hit LP
    // RDKit✔️✔️:     if (atom) {
    // RDKit✔️✔️:       int aid = res->addAtom(atom, false, true);
    // RDKit✔️✔️:       idxCorresp[i] = aid;
    // RDKit✔️✔️:       threeDPs.push_back(pos);
    // RDKit✔️✔️:       if (atom->getSymbol() == "H") {
    // RDKit✔️✔️:         hasHAtoms = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       ++nLP;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // mol2 files need to have hydrogen atoms otherwise formal charge estimation
    // RDKit✔️✔️:   // will be problematic
    // RDKit✔️✔️:   if (!hasHAtoms) {
    // RDKit✔️✔️:     std::string nm;
    // RDKit✔️✔️:     res->getProp(common_properties::_Name, nm);
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog) << nm
    // RDKit✔️✔️:                             << ": Warning - no explicit hydrogens in "
    // RDKit✔️✔️:                                "mol2 file but needed for formal charge "
    // RDKit✔️✔️:                                "estimation."
    // RDKit✔️✔️:                             << std::endl;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // create conformer based on 3DPoints and add to RWMol
    // RDKit✔️✔️:   auto *conf = new Conformer(nAtoms - nLP);
    // RDKit✔️✔️:   std::vector<RDGeom::Point3D>::const_iterator threeDPsIt = threeDPs.begin();
    // RDKit✔️✔️:   for (unsigned int i = 0; i < threeDPs.size(); ++i) {
    // RDKit✔️✔️:     conf->setAtomPos(i, *threeDPsIt);
    // RDKit✔️✔️:     ++threeDPsIt;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res->addConformer(conf, true);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   POSTCONDITION(res->getNumAtoms() == (nAtoms - nLP),
    // RDKit✔️✔️:                 "Wrong number of atoms in molecule");
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: ParseMol2AtomBlock
    let n_atoms_usize = usize::try_from(n_atoms)
        .map_err(|_| Mol2ReadError::Parse("atom count does not fit platform usize".to_string()))?;
    let mut builder = MoleculeBuilder::new();
    let mut idx_corresp = vec![None; n_atoms_usize];
    let mut has_h_atoms = false;
    let mut three_dps = Vec::with_capacity(n_atoms_usize);
    let mut offset = atom_start;

    for slot in idx_corresp.iter_mut().take(n_atoms_usize) {
        let temp_str = rdkit_get_line_at_checked(input, &mut offset)?;
        if let Some(atom_line) = parse_mol2_atom_line_like_rdkit(temp_str)? {
            let is_hydrogen = atom_line.atom_spec.element() == Element::H;
            let atom_id = builder.add_atom(atom_line.atom_spec);
            *slot = Some(atom_id);
            three_dps.push(atom_line.position);
            if is_hydrogen {
                has_h_atoms = true;
            }
        }
    }

    builder
        .add_conformer(Conformer3D::new(0, three_dps, true))
        .map_err(|err| Mol2ReadError::Parse(err.to_string()))?;

    let actual_atoms = builder.atoms().len();
    let expected_atoms = idx_corresp.iter().filter(|atom| atom.is_some()).count();
    if actual_atoms != expected_atoms {
        return Err(Mol2ReadError::Parse(
            "Wrong number of atoms in molecule".to_string(),
        ));
    }

    Ok(Mol2AtomBlock {
        builder,
        idx_corresp,
        has_h_atoms,
    })
}

fn parse_mol2_bond_line_like_rdkit(
    bond_line: &str,
    idx_corresp: &[Option<AtomId>],
) -> Result<Option<Mol2BondLine>, Mol2ReadError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: ParseMol2FileBondLine
    // RDKit✔️✔️: Bond *ParseMol2FileBondLine(const std::string bondLine,
    // RDKit✔️✔️:                             const INT_VECT &idxCorresp) {
    // RDKit✔️✔️:   unsigned int idx1, idx2;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
    // RDKit✔️✔️:   boost::char_separator<char> sep(" \t\n");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   tokenizer tokens(bondLine, sep);
    // RDKit✔️✔️:   tokenizer::const_iterator itemIt = tokens.begin();
    // RDKit✔️✔️:   if (itemIt == tokens.end()) {
    // RDKit✔️✔️:     throw FileParseException("no info in mol2 bond line");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     // tripos bond id skip
    // RDKit✔️✔️:     ++itemIt;
    // RDKit✔️✔️:     if (itemIt == tokens.end()) {
    // RDKit✔️✔️:       throw FileParseException("no info in mol2 bond line");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     idx1 = boost::lexical_cast<unsigned int>(*itemIt);
    // RDKit✔️✔️:     ++itemIt;
    // RDKit✔️✔️:     if (itemIt == tokens.end()) {
    // RDKit✔️✔️:       throw FileParseException("no info in mol2 bond line");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     idx2 = boost::lexical_cast<unsigned int>(*itemIt);
    // RDKit✔️✔️:     ++itemIt;
    // RDKit✔️✔️:     if (itemIt == tokens.end()) {
    // RDKit✔️✔️:       throw FileParseException("no info in mol2 bond line");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️✔️:     throw FileParseException("Cannot process mol2 bonds.");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // adjust the numbering
    // RDKit✔️✔️:   idx1--;
    // RDKit✔️✔️:   idx2--;
    // RDKit✔️✔️:
    // RDKit❗✔️:   // if either of both ends of the bond is not an atom in the mol - return NULL
    // RDKit❗✔️:   if (!(idx1 < idxCorresp.size() || idx2 < idxCorresp.size())) {
    // RDKit❗✔️:     throw FileParseException("index mismatch");
    // RDKit❗✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (idxCorresp[idx1] < 0 || idxCorresp[idx2] < 0) {
    // RDKit✔️✔️:     return nullptr;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // lexical casts of bond types is not really useful in this case since we
    // RDKit✔️✔️:   // could also have strings as values ...
    // RDKit✔️✔️:   // use if/else if/end statements for now - maybe change to mapping later - no
    // RDKit✔️✔️:   // idea if it's worth it ...?
    // RDKit✔️✔️:   Bond::BondType type;
    // RDKit✔️✔️:   std::string tBType = *itemIt;
    // RDKit✔️✔️:   if (tBType == "1" || tBType == "am") {
    // RDKit✔️✔️:     type = Bond::SINGLE;
    // RDKit✔️✔️:   } else if (tBType == "2") {
    // RDKit✔️✔️:     type = Bond::DOUBLE;
    // RDKit✔️✔️:   } else if (tBType == "3") {
    // RDKit✔️✔️:     type = Bond::TRIPLE;
    // RDKit✔️✔️:   } else if (tBType == "ar") {
    // RDKit✔️✔️:     type = Bond::AROMATIC;
    // RDKit✔️✔️:   } else if (tBType == "du" || tBType == "un") {
    // RDKit✔️✔️:     type = Bond::UNSPECIFIED;  // have to come back here - see if there is any
    // RDKit✔️✔️:                                // comment in th documentation ...
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     // this happens only if some weird thing is in the file or if we encounter a
    // RDKit✔️✔️:     // "nc" - not connected ...
    // RDKit✔️✔️:     // but why would anyone specify a bond which is not connected?
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog) << "Warning - unsupported bond type: " << tBType
    // RDKit✔️✔️:                             << " ignored!" << std::endl;
    // RDKit✔️✔️:     return nullptr;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto *res = new Bond(type);
    // RDKit✔️✔️:   res->setBeginAtomIdx(idxCorresp[idx1]);
    // RDKit✔️✔️:   res->setEndAtomIdx(idxCorresp[idx2]);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: ParseMol2FileBondLine
    let mut tokens = bond_line
        .split([' ', '\t', '\n'])
        .filter(|token| !token.is_empty());
    let _tripos_bond_id = tokens
        .next()
        .ok_or_else(|| Mol2ReadError::Parse("no info in mol2 bond line".to_string()))?;

    let idx1_token = tokens
        .next()
        .ok_or_else(|| Mol2ReadError::Parse("no info in mol2 bond line".to_string()))?;
    let idx1 = parse_rdkit_unsigned_int_token(idx1_token)
        .ok_or_else(|| Mol2ReadError::Parse("Cannot process mol2 bonds.".to_string()))?;

    let idx2_token = tokens
        .next()
        .ok_or_else(|| Mol2ReadError::Parse("no info in mol2 bond line".to_string()))?;
    let idx2 = parse_rdkit_unsigned_int_token(idx2_token)
        .ok_or_else(|| Mol2ReadError::Parse("Cannot process mol2 bonds.".to_string()))?;

    let sybyl_bond_type = tokens
        .next()
        .ok_or_else(|| Mol2ReadError::Parse("no info in mol2 bond line".to_string()))?
        .to_string();

    let idx1 = idx1.wrapping_sub(1);
    let idx2 = idx2.wrapping_sub(1);
    let idx1_usize =
        usize::try_from(idx1).map_err(|_| Mol2ReadError::Parse("index mismatch".to_string()))?;
    let idx2_usize =
        usize::try_from(idx2).map_err(|_| Mol2ReadError::Parse("index mismatch".to_string()))?;

    let Some(begin) = idx_corresp.get(idx1_usize).copied().flatten() else {
        if idx1_usize >= idx_corresp.len() {
            return Err(Mol2ReadError::Parse("index mismatch".to_string()));
        }
        return Ok(None);
    };
    let Some(end) = idx_corresp.get(idx2_usize).copied().flatten() else {
        if idx2_usize >= idx_corresp.len() {
            return Err(Mol2ReadError::Parse("index mismatch".to_string()));
        }
        return Ok(None);
    };

    let order = match sybyl_bond_type.as_str() {
        "1" | "am" => BondOrder::Single,
        "2" => BondOrder::Double,
        "3" => BondOrder::Triple,
        "ar" => BondOrder::Aromatic,
        "du" | "un" => BondOrder::Unspecified,
        _ => return Ok(None),
    };

    Ok(Some(Mol2BondLine {
        bond_spec: BondSpec::new(begin, end, order),
        sybyl_bond_type,
        zero_based_begin_idx: idx1,
        zero_based_end_idx: idx2,
    }))
}

fn parse_mol2_bond_block_like_rdkit(
    input: &str,
    bond_start: usize,
    n_bonds: u32,
    atom_block: Mol2AtomBlock,
) -> Result<Mol2BondBlock, Mol2ReadError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: ParseMol2BondBlock
    // RDKit✔️✔️: void ParseMol2BondBlock(std::istream *inStream, RWMol *res, unsigned int nBonds,
    // RDKit✔️✔️:                         const INT_VECT &idxCorresp) {
    // RDKit✔️✔️:   PRECONDITION(inStream, "inStream not valid");
    // RDKit✔️✔️:   PRECONDITION(!inStream->eof(), "inStream is at eof");
    // RDKit✔️✔️:   PRECONDITION(res, "RWMol not valid");
    // RDKit✔️✔️:   unsigned int nBadBonds = 0;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nBonds; ++i) {
    // RDKit✔️✔️:     std::string tempStr = getLine(inStream);
    // RDKit✔️✔️:     if (inStream->eof()) {
    // RDKit✔️✔️:       throw FileParseException("premature EOF");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     Bond *bond = ParseMol2FileBondLine(tempStr, idxCorresp);
    // RDKit✔️✔️:     // if something weird happened there will be no bond for that line
    // RDKit✔️✔️:     if (bond) {
    // RDKit✔️✔️:       // if we got an aromatic bond set the flag on the bond and the connected
    // RDKit✔️✔️:       // atoms
    // RDKit✔️✔️:       if (bond->getBondType() == Bond::AROMATIC) {
    // RDKit✔️✔️:         bond->setIsAromatic(true);
    // RDKit✔️✔️:         res->getAtomWithIdx(bond->getBeginAtomIdx())->setIsAromatic(true);
    // RDKit✔️✔️:         res->getAtomWithIdx(bond->getEndAtomIdx())->setIsAromatic(true);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       res->addBond(bond, true);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       nBadBonds++;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   POSTCONDITION(res->getNumBonds() == (nBonds - nBadBonds),
    // RDKit✔️✔️:                 "Wrong number of atoms in molecule");
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: ParseMol2BondBlock
    let mut builder = atom_block.builder;
    let mut n_bad_bonds = 0_u32;
    let mut offset = bond_start;

    for _ in 0..n_bonds {
        let temp_str = rdkit_get_line_at_checked(input, &mut offset)?;
        let Some(mut bond_line) =
            parse_mol2_bond_line_like_rdkit(temp_str, &atom_block.idx_corresp)?
        else {
            n_bad_bonds = n_bad_bonds.wrapping_add(1);
            continue;
        };

        if bond_line.bond_spec.order() == BondOrder::Aromatic {
            let begin = bond_line.bond_spec.begin();
            let end = bond_line.bond_spec.end();
            bond_line.bond_spec = bond_line.bond_spec.with_aromatic(true);
            let Some(begin_atom) = builder.atom_mut(begin) else {
                return Err(Mol2ReadError::Parse("index mismatch".to_string()));
            };
            begin_atom.set_aromatic(true);
            let Some(end_atom) = builder.atom_mut(end) else {
                return Err(Mol2ReadError::Parse("index mismatch".to_string()));
            };
            end_atom.set_aromatic(true);
        }

        builder
            .add_bond(bond_line.bond_spec)
            .map_err(|err| Mol2ReadError::Parse(err.to_string()))?;
    }

    let expected_bonds = usize::try_from(n_bonds.wrapping_sub(n_bad_bonds))
        .map_err(|_| Mol2ReadError::Parse("bond count does not fit platform usize".to_string()))?;
    if builder.bonds().len() != expected_bonds {
        return Err(Mol2ReadError::Parse(
            "Wrong number of atoms in molecule".to_string(),
        ));
    }

    Ok(Mol2BondBlock {
        builder,
        n_bad_bonds,
    })
}

fn read_formal_charges_from_attr_like_rdkit(
    input: &str,
    charge_start: usize,
    bond_block: Mol2BondBlock,
) -> Result<Mol2UnityAtomAttrBlock, Mol2ReadError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: readFormalChargesFromAttr
    // RDKit✔️✔️: void readFormalChargesFromAttr(std::istream *inStream, RWMol *res) {
    // RDKit✔️✔️:   PRECONDITION(inStream, "inStream not valid");
    // RDKit✔️✔️:   PRECONDITION(!inStream->eof(), "inStream is at eof");
    // RDKit✔️✔️:   PRECONDITION(res, "RWMol not valid");
    // RDKit✔️✔️:   typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
    // RDKit✔️✔️:   boost::char_separator<char> sep(" \t\n");
    // RDKit✔️✔️:   bool readNextAtomAttribs = true;
    // RDKit✔️✔️:   unsigned int atomIdx = 0, noAtomAttr = 0;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // std::streampos stPos = inStream->tellg();
    // RDKit✔️✔️:   std::string tempStr = getLine(inStream);
    // RDKit✔️✔️:   // there needs to be at least one entry
    // RDKit✔️✔️:   if (inStream->eof()) {
    // RDKit✔️✔️:     throw FileParseException("premature EOF in readFormalCharges");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   while (readNextAtomAttribs) {
    // RDKit✔️✔️:     tokenizer tokens(tempStr, sep);
    // RDKit✔️✔️:     tokenizer::const_iterator itemIt = tokens.begin();
    // RDKit✔️✔️:     try {
    // RDKit✔️✔️:       atomIdx = boost::lexical_cast<unsigned int>(*itemIt);
    // RDKit✔️✔️:       ++itemIt;
    // RDKit✔️✔️:       noAtomAttr = boost::lexical_cast<unsigned int>(*itemIt);
    // RDKit✔️✔️:     } catch (boost::bad_lexical_cast &) {
    // RDKit✔️✔️:       throw FileParseException("Cannot process mol2 UnityAtomAttr.");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (unsigned int i = 0; i < noAtomAttr; ++i) {
    // RDKit✔️✔️:       std::string tempStr = getLine(inStream);
    // RDKit✔️✔️:       if (inStream->eof()) {
    // RDKit✔️✔️:         throw FileParseException("premature EOF in readFormalCharges");
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       tokenizer atmTokens(tempStr, sep);
    // RDKit✔️✔️:       itemIt = atmTokens.begin();
    // RDKit✔️✔️:       if (*itemIt == "AtomExpr") {
    // RDKit✔️✔️:         int formCharge = 0;
    // RDKit✔️✔️:         ++itemIt;
    // RDKit✔️✔️:         // FIX:what if an atom has multiple properties? Seems like they might be
    // RDKit✔️✔️:         // separated
    // RDKit✔️✔️:         // with ";" ... need to look at that in more detail!
    // RDKit✔️✔️:         if ((*itemIt).find("=") == std::string::npos) {
    // RDKit✔️✔️:           try {
    // RDKit✔️✔️:             formCharge = boost::lexical_cast<int>(*itemIt);
    // RDKit✔️✔️:           } catch (boost::bad_lexical_cast &) {
    // RDKit✔️✔️:             readNextAtomAttribs = false;
    // RDKit✔️✔️:             throw FileParseException("Cannot process mol2 formal charge.");
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           // assign the charge
    // RDKit✔️✔️:           res->getAtomWithIdx(atomIdx - 1)->setFormalCharge(formCharge);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }  // endfor
    // RDKit✔️✔️:     // check if we have finished reading UNITY_ATOM_ATTR
    // RDKit✔️✔️:     if (inStream->eof()) {
    // RDKit✔️✔️:       readNextAtomAttribs = false;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       tempStr = getLine(inStream);
    // RDKit✔️✔️:       if (tempStr == "" || tempStr[0] == '@' || tempStr[0] == '#') {
    // RDKit✔️✔️:         readNextAtomAttribs = false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: readFormalChargesFromAttr
    let mut builder = bond_block.builder;
    let mut offset = charge_start;
    let mut temp_str = rdkit_get_line_at_checked(input, &mut offset)?;

    loop {
        let mut tokens = temp_str
            .split([' ', '\t', '\n'])
            .filter(|token| !token.is_empty());
        let atom_idx_token = tokens.next().ok_or_else(|| {
            Mol2ReadError::Parse("Cannot process mol2 UnityAtomAttr.".to_string())
        })?;
        let atom_idx = parse_rdkit_unsigned_int_token(atom_idx_token).ok_or_else(|| {
            Mol2ReadError::Parse("Cannot process mol2 UnityAtomAttr.".to_string())
        })?;
        let no_atom_attr_token = tokens.next().ok_or_else(|| {
            Mol2ReadError::Parse("Cannot process mol2 UnityAtomAttr.".to_string())
        })?;
        let no_atom_attr = parse_rdkit_unsigned_int_token(no_atom_attr_token).ok_or_else(|| {
            Mol2ReadError::Parse("Cannot process mol2 UnityAtomAttr.".to_string())
        })?;

        for _ in 0..no_atom_attr {
            let attr_line = rdkit_get_line_at_checked(input, &mut offset)?;
            let mut attr_tokens = attr_line
                .split([' ', '\t', '\n'])
                .filter(|token| !token.is_empty());
            let Some(attr_name) = attr_tokens.next() else {
                continue;
            };
            if attr_name == "AtomExpr" {
                let Some(charge_token) = attr_tokens.next() else {
                    return Err(Mol2ReadError::Parse(
                        "Cannot process mol2 formal charge.".to_string(),
                    ));
                };
                if !charge_token.contains('=') {
                    let form_charge = charge_token.parse::<i32>().map_err(|_| {
                        Mol2ReadError::Parse("Cannot process mol2 formal charge.".to_string())
                    })?;
                    let form_charge = i8::try_from(form_charge).map_err(|_| {
                        Mol2ReadError::Parse("Cannot process mol2 formal charge.".to_string())
                    })?;
                    let atom_idx = atom_idx.wrapping_sub(1);
                    let atom_idx = usize::try_from(atom_idx)
                        .map_err(|_| Mol2ReadError::Parse("index mismatch".to_string()))?;
                    let Some(atom) = builder.atom_mut(AtomId::new(atom_idx)) else {
                        return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                    };
                    atom.set_formal_charge(form_charge);
                }
            }
        }

        if offset >= input.len() {
            break;
        }
        temp_str = rdkit_get_line_at(input, &mut offset);
        if temp_str.is_empty() || temp_str.starts_with('@') || temp_str.starts_with('#') {
            break;
        }
    }

    Ok(Mol2UnityAtomAttrBlock {
        builder,
        next_offset: offset,
    })
}

fn fix_nitro_substructure_and_charge_like_rdkit(
    builder: &mut MoleculeBuilder,
    atom_id: AtomId,
) -> Result<(), Mol2ReadError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: fixNitroSubstructureAndCharge
    // RDKit✔️✔️: void fixNitroSubstructureAndCharge(RWMol *res, unsigned int atIdx) {
    // RDKit✔️✔️:   unsigned int noODblNeighbors = 0;
    // RDKit✔️✔️:   ROMol::ADJ_ITER nbrIdxIt, nbrEndIdxIt;
    // RDKit✔️✔️:   unsigned int toModIdx = 0;
    // RDKit✔️✔️:   boost::tie(nbrIdxIt, nbrEndIdxIt) =
    // RDKit✔️✔️:       res->getAtomNeighbors(res->getAtomWithIdx(atIdx));
    // RDKit✔️✔️:   while (nbrIdxIt != nbrEndIdxIt) {
    // RDKit✔️✔️:     Bond *curBond = res->getBondBetweenAtoms(atIdx, *nbrIdxIt);
    // RDKit✔️✔️:     if (res->getAtomWithIdx(*nbrIdxIt)->getAtomicNum() == 8 &&
    // RDKit✔️✔️:         curBond->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:       ++noODblNeighbors;
    // RDKit✔️✔️:       toModIdx = *nbrIdxIt;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++nbrIdxIt;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (noODblNeighbors == 2) {
    // RDKit✔️✔️:     res->getBondBetweenAtoms(atIdx, toModIdx)->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:     res->getAtomWithIdx(atIdx)->setFormalCharge(1);
    // RDKit✔️✔️:     res->getAtomWithIdx(toModIdx)->setFormalCharge(-1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: fixNitroSubstructureAndCharge
    if atom_id.index() >= builder.atoms().len() {
        return Err(Mol2ReadError::Parse("index mismatch".to_string()));
    }

    let mut no_o_double_neighbors = 0_u32;
    let mut to_mod_idx = None;
    let mut to_mod_bond = None;

    for bond_id in builder.neighbor_bonds(atom_id).iter().copied() {
        let bond = builder
            .bond(bond_id)
            .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_string()))?;
        let neighbor_id = if bond.begin() == atom_id {
            bond.end()
        } else {
            bond.begin()
        };
        let neighbor = builder
            .atoms()
            .get(neighbor_id.index())
            .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_string()))?;
        if neighbor.atomic_number() == 8 && bond.order() == BondOrder::Double {
            no_o_double_neighbors = no_o_double_neighbors.wrapping_add(1);
            to_mod_idx = Some(neighbor_id);
            to_mod_bond = Some(bond_id);
        }
    }

    if no_o_double_neighbors == 2 {
        let to_mod_bond =
            to_mod_bond.ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_string()))?;
        let to_mod_idx =
            to_mod_idx.ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_string()))?;
        builder
            .set_bond_order(to_mod_bond, BondOrder::Single)
            .map_err(|err| Mol2ReadError::Parse(err.to_string()))?;
        let Some(atom) = builder.atom_mut(atom_id) else {
            return Err(Mol2ReadError::Parse("index mismatch".to_string()));
        };
        atom.set_formal_charge(1);
        let Some(atom) = builder.atom_mut(to_mod_idx) else {
            return Err(Mol2ReadError::Parse("index mismatch".to_string()));
        };
        atom.set_formal_charge(-1);
    }

    Ok(())
}

fn guess_formal_charges_like_rdkit(builder: &mut MoleculeBuilder) -> Result<(), Mol2ReadError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: guessFormalCharges
    // RDKit✔️✔️: void guessFormalCharges(RWMol *res) {
    // RDKit✔️✔️:   // FIX: this whole thing has problems with positively charged pyridines et al.
    // RDKit✔️✔️:   for (RWMol::AtomIterator atomIt = res->beginAtoms();
    // RDKit✔️✔️:        atomIt != res->endAtoms(); ++atomIt) {
    // RDKit✔️✔️:     Atom *at = (*atomIt);
    // RDKit✔️✔️:     // assign only if no formal charge set on that atom and atom is not carbon
    // RDKit✔️✔️:     // (the latter
    // RDKit✔️✔️:     // might needs changing later on - let's see) and not for query atoms (dummy
    // RDKit✔️✔️:     // etc.)
    // RDKit✔️✔️:     // since this happens only during cleanup of bad substructures
    // RDKit✔️✔️:     // FIX: check nitro compounds et al.
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (at->getFormalCharge() == 0 && at->getSymbol() != "C" &&
    // RDKit✔️✔️:         !(at->hasQuery())) {
    // RDKit✔️✔️:       // have to calculate the explicit valence without call to
    // RDKit✔️✔️:       // calcExplicitValence since
    // RDKit✔️✔️:       // this will barf when I have e.g. an uncharged 4 valent nitrogen ...
    // RDKit✔️✔️:       int noAromBonds = 0;
    // RDKit✔️✔️:       double accum = 0;  // FIX: could this give non int values ?
    // RDKit✔️✔️:       for (const auto bnd : res->atomBonds(at)) {
    // RDKit✔️✔️:         accum += bnd->getValenceContrib(at);
    // RDKit✔️✔️:         if (bnd->getBondType() == Bond::AROMATIC) {
    // RDKit✔️✔️:           ++noAromBonds;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // Assumption: if there is an aromatic bridge atom the accum will be 4.5
    // RDKit✔️✔️:       //(three aromatic bonds), e.g. naphthalenes. However those are not charged
    // RDKit✔️✔️:       // so we
    // RDKit✔️✔️:       // can stop here
    // RDKit✔️✔️:       // FIX: a structure such as c12cccc[n+]2cccn1 will not work if we just
    // RDKit✔️✔️:       // continue
    // RDKit✔️✔️:       // for noAromBonds>2
    // RDKit✔️✔️:       // special case checking - if atom c then go on
    // RDKit✔️✔️:       if (noAromBonds > 2 && at->getSymbol() == "C") {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // dbtranslate problems - for mols with no UNITY_ATOM_ATTR set - we won't
    // RDKit✔️✔️:       // be able to guess
    // RDKit✔️✔️:       // if this mol needs to be guessed or not ... hence, we don't guess on
    // RDKit✔️✔️:       // atoms with
    // RDKit✔️✔️:       // ar specification and not atom type X.ar and in 5-membered ring (there
    // RDKit✔️✔️:       // is stuff like
    // RDKit✔️✔️:       // c1cc[o+]cc1 that should be charged
    // RDKit✔️✔️:       // e.g. 5ring with N.pl3 as NH atom or other atoms without ar
    // RDKit✔️✔️:       // specification in aromatic ring
    // RDKit✔️✔️:       // FIX: do we need make sure this only happens for atoms in ring?
    // RDKit✔️✔️:       if (!res->getRingInfo()->isSssrOrBetter()) {
    // RDKit✔️✔️:         MolOps::findSSSR(*res);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       auto tATT = at->getProp<std::string>(common_properties::_TriposAtomType);
    // RDKit✔️✔️:       if (at->getIsAromatic() && tATT.find("ar") == std::string::npos &&
    // RDKit✔️✔️:           res->getRingInfo()->isAtomInRingOfSize(at->getIdx(), 5)) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // for dbtranslate a problem will also occur for N.ar with 3 aromatic
    // RDKit✔️✔️:       // bonds and no charge
    // RDKit✔️✔️:       // assigned for these we don't assign charges (corina will create
    // RDKit✔️✔️:       // kekulized input for this
    // RDKit✔️✔️:       //(at least in most cases) - anyway, throw a warning!
    // RDKit✔️✔️:       if (noAromBonds == 3 && tATT == "N.ar") {
    // RDKit✔️✔️:         auto nm = res->getProp<std::string>(common_properties::_Name);
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << nm
    // RDKit✔️✔️:             << ": warning - aromatic N with 3 aromatic bonds - "
    // RDKit✔️✔️:                "skipping charge guess for this atom"
    // RDKit✔️✔️:             << std::endl;
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // sometimes things like benzimidazoles can have only one bond of the
    // RDKit✔️✔️:       // imidazole ring as aromatic and the other one as a single bond ...
    // RDKit✔️✔️:       // catch that this way - see also the trick from GL
    // RDKit✔️✔️:       auto expVal = static_cast<int>(std::round(accum + 0.1));
    // RDKit✔️✔️:       const auto &valens =
    // RDKit✔️✔️:           PeriodicTable::getTable()->getValenceList(at->getAtomicNum());
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // check default valence and compare to expVal - chg
    // RDKit✔️✔️:       // the hypothesis is that we prefer positively charged atoms over
    // RDKit✔️✔️:       // negatively charged ones
    // RDKit✔️✔️:       // for multi default valence atoms (e.g. CS(O)(O) should end up being
    // RDKit✔️✔️:       // C[S+]([O-])[O-] rather
    // RDKit✔️✔️:       // than C[S-][O-][O-] but that might change based no different examples
    // RDKit✔️✔️:       int nElectrons =
    // RDKit✔️✔️:           PeriodicTable::getTable()->getNouterElecs(at->getAtomicNum());
    // RDKit✔️✔️:       int assignChg;
    // RDKit✔️✔️:       if (nElectrons >= 4) {
    // RDKit✔️✔️:         assignChg = expVal - valens.front();
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         assignChg = valens.front() - expVal;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (assignChg > 0 && nElectrons >= 4) {
    // RDKit✔️✔️:         for (auto vi : valens) {
    // RDKit✔️✔️:           // Since we do this only for nocharged atoms we can get away without
    // RDKit✔️✔️:           // including the charge into this apart from that we do not assign
    // RDKit✔️✔️:           // charges higher than +/- 1 for atoms with multiple valence states
    // RDKit✔️✔️:           // otherwise the early break would have to go away which in turn would
    // RDKit✔️✔️:           // result in things like [S+4] for sulfonamides
    // RDKit✔️✔️:           assignChg = expVal - vi;
    // RDKit✔️✔️:           if (vi <= expVal && assignChg < 2) {
    // RDKit✔️✔️:             break;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (assignChg) {
    // RDKit✔️✔️:         // no aromatic atom will get aa abs(charge) > 1
    // RDKit✔️✔️:         if (at->getIsAromatic() && abs(assignChg) > 1) {
    // RDKit✔️✔️:           at->setFormalCharge((assignChg > 0) -
    // RDKit✔️✔️:                               (assignChg < 0));  // this results in -1 or +1
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           at->setFormalCharge(assignChg);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         // corina will create strange nitro groups which look like N(=O)(=O)
    // RDKit✔️✔️:         // which in turn will result in
    // RDKit✔️✔️:         //[N+2](=O)(=O) this needs to be fixed at this stage since otherwise we
    // RDKit✔️✔️:         // would have to check on all
    // RDKit✔️✔️:         // N.pl3 during the cleanup substructures step
    // RDKit✔️✔️:         if (assignChg == 2 && expVal == 5 && at->getSymbol() == "N") {
    // RDKit✔️✔️:           fixNitroSubstructureAndCharge(res, at->getIdx());
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         // what if expVal != at->calcExplicitValence();?
    // RDKit✔️✔️:         // cannot imagine a case where that will happen now.
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: guessFormalCharges
    let mut ring_info = None;

    for idx in 0..builder.atoms().len() {
        let atom_id = AtomId::new(idx);
        let atom = &builder.atoms()[idx];
        if atom.formal_charge() != 0 || atom.atomic_number() == 6 || atom.query().is_some() {
            continue;
        }

        let mut no_aromatic_bonds = 0_i32;
        let mut accum = 0.0_f64;
        for bond_id in builder.neighbor_bonds(atom_id).iter().copied() {
            let bond = builder
                .bond(bond_id)
                .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_string()))?;
            accum += bond_valence_contrib(bond, atom_id)
                .map_err(|err| Mol2ReadError::Parse(err.to_string()))?;
            if bond.order() == BondOrder::Aromatic {
                no_aromatic_bonds += 1;
            }
        }

        if no_aromatic_bonds > 2 && builder.atoms()[idx].atomic_number() == 6 {
            continue;
        }

        if ring_info.is_none() {
            let adjacency =
                AdjacencyList::try_from_topology(builder.atoms().len(), builder.bonds())
                    .map_err(|err| Mol2ReadError::Parse(err.to_string()))?;
            ring_info = Some(
                find_sssr_from_parts(builder.atoms().len(), builder.bonds(), &adjacency)
                    .map_err(|err| Mol2ReadError::Parse(err.to_string()))?,
            );
        }
        let tripos_atom_type = builder.atoms()[idx]
            .prop("_TriposAtomType")
            .ok_or_else(|| Mol2ReadError::Parse("Missing _TriposAtomType".to_string()))?
            .to_string();
        if builder.atoms()[idx].is_aromatic()
            && !tripos_atom_type.contains("ar")
            && ring_info
                .as_ref()
                .is_some_and(|info| info.is_atom_in_ring_of_size(atom_id, 5))
        {
            continue;
        }

        if no_aromatic_bonds == 3 && tripos_atom_type == "N.ar" {
            continue;
        }

        let explicit_valence = (accum + 0.1).round() as i32;
        let valences = rdkit_valence_list(builder.atoms()[idx].atomic_number())
            .map_err(|err| Mol2ReadError::Parse(err.to_string()))?
            .ok_or_else(|| Mol2ReadError::Parse("Missing RDKit valence list".to_string()))?;
        let n_electrons = periodic_table_outer_electrons(builder.atoms()[idx].atomic_number())
            .map_err(|err| Mol2ReadError::Parse(err.to_string()))?;
        let mut assign_charge = if n_electrons >= 4 {
            explicit_valence - valences[0]
        } else {
            valences[0] - explicit_valence
        };
        if assign_charge > 0 && n_electrons >= 4 {
            for &valence in valences {
                assign_charge = explicit_valence - valence;
                if valence <= explicit_valence && assign_charge < 2 {
                    break;
                }
            }
        }
        if assign_charge != 0 {
            let formal_charge = if builder.atoms()[idx].is_aromatic()
                && assign_charge.unsigned_abs() > 1
            {
                (assign_charge > 0) as i8 - (assign_charge < 0) as i8
            } else {
                i8::try_from(assign_charge)
                    .map_err(|_| Mol2ReadError::Parse("formal charge out of range".to_string()))?
            };
            let Some(atom) = builder.atom_mut(atom_id) else {
                return Err(Mol2ReadError::Parse("index mismatch".to_string()));
            };
            atom.set_formal_charge(formal_charge);
            if assign_charge == 2
                && explicit_valence == 5
                && builder.atoms()[idx].atomic_number() == 7
            {
                fix_nitro_substructure_and_charge_like_rdkit(builder, atom_id)?;
            }
        }
    }

    Ok(())
}

fn check_no_h_neighbors_n_oxide_like_rdkit(
    builder: &MoleculeBuilder,
    atom_id: AtomId,
    to_mod_idx: &mut Option<AtomId>,
) -> Result<u32, Mol2ReadError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: chkNoHNeighbNOx
    // RDKit✔️✔️: unsigned int chkNoHNeighbNOx(RWMol *res, ROMol::ADJ_ITER atIdxIt,
    // RDKit✔️✔️:                              int &toModIdx) {
    // RDKit✔️✔️:   Atom *at = res->getAtomWithIdx(*atIdxIt);
    // RDKit✔️✔️:   unsigned int noHNbrs = 0;
    // RDKit✔️✔️:   ROMol::ADJ_ITER nbrIdxIt, nbrEndIdxIt;
    // RDKit✔️✔️:   boost::tie(nbrIdxIt, nbrEndIdxIt) = res->getAtomNeighbors(at);
    // RDKit✔️✔️:   while (nbrIdxIt != nbrEndIdxIt) {
    // RDKit✔️✔️:     if (res->getAtomWithIdx(*nbrIdxIt)->getAtomicNum() == 1) {
    // RDKit✔️✔️:       ++noHNbrs;
    // RDKit✔️✔️:     } else if (res->getAtomWithIdx(*nbrIdxIt)->getAtomicNum() == 8 &&
    // RDKit✔️✔️:                res->getAtomDegree(res->getAtomWithIdx(*nbrIdxIt)) == 1) {
    // RDKit✔️✔️:       // this is a N in an N-oxide constellation
    // RDKit✔️✔️:       // we can do the above if clause since mol2 have explicit hydrogens
    // RDKit✔️✔️:       toModIdx = *atIdxIt;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++nbrIdxIt;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return noHNbrs;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: chkNoHNeighbNOx
    if atom_id.index() >= builder.atoms().len() {
        return Err(Mol2ReadError::Parse("index mismatch".to_string()));
    }

    let mut no_h_neighbors = 0_u32;

    for bond_id in builder.neighbor_bonds(atom_id).iter().copied() {
        let bond = builder
            .bond(bond_id)
            .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_string()))?;
        let neighbor_id = if bond.begin() == atom_id {
            bond.end()
        } else {
            bond.begin()
        };
        let neighbor = builder
            .atoms()
            .get(neighbor_id.index())
            .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_string()))?;
        if neighbor.atomic_number() == 1 {
            no_h_neighbors = no_h_neighbors.wrapping_add(1);
        } else if neighbor.atomic_number() == 8 && builder.degree(neighbor_id) == 1 {
            *to_mod_idx = Some(atom_id);
        }
    }

    Ok(no_h_neighbors)
}

fn clean_up_mol2_substructures_like_rdkit(
    builder: &mut MoleculeBuilder,
) -> Result<bool, Mol2ReadError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: cleanUpMol2Substructures
    // RDKit✔️✔️: bool cleanUpMol2Substructures(RWMol *res) {
    // RDKit✔️✔️:   // NOTE: check the nitro fix in guess formal charges!
    // RDKit✔️✔️:   boost::dynamic_bitset<> isFixed(res->getNumAtoms());
    // RDKit✔️✔️:   for (auto at : res->atoms()) {
    // RDKit✔️✔️:     unsigned int idx = at->getIdx();
    // RDKit✔️✔️:     // make sure we haven't finished this atom already
    // RDKit✔️✔️:     if (isFixed[idx]) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     auto tAT = at->getProp<std::string>(common_properties::_TriposAtomType);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (tAT == "N.4") {
    // RDKit✔️✔️:       at->setFormalCharge(1);
    // RDKit✔️✔️:     } else if (tAT == "O.co2") {
    // RDKit✔️✔️:       // negatively charged carboxylates with O.co2
    // RDKit✔️✔️:       // according to Tripos, those should only appear in carboxylates and
    // RDKit✔️✔️:       // phosphates,
    // RDKit✔️✔️:       if (at->getDegree() != 1) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "Warning - O.co2 with degree >1." << std::endl;
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       auto nbrs = res->atomNeighbors(at);
    // RDKit✔️✔️:       // this should return only the C.2
    // RDKit✔️✔️:       auto nbr = *nbrs.begin();
    // RDKit✔️✔️:       auto tATT = nbr->getProp<std::string>(common_properties::_TriposAtomType);
    // RDKit✔️✔️:       if (tATT == "P.3") {
    // RDKit✔️✔️:         // special case for phosphates
    // RDKit✔️✔️:         // we keep the first bond to O.co2 as double and make the rest single
    // RDKit✔️✔️:         Bond *b = res->getBondBetweenAtoms(idx, nbr->getIdx());
    // RDKit✔️✔️:         b->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:         b->setIsAromatic(false);
    // RDKit✔️✔️:         at->setIsAromatic(false);
    // RDKit✔️✔️:         isFixed[idx] = 1;
    // RDKit✔️✔️:         for (auto onbr : res->atomNeighbors(nbr)) {
    // RDKit✔️✔️:           if (onbr->getAtomicNum() == 8 && !isFixed[onbr->getIdx()] &&
    // RDKit✔️✔️:               onbr->getProp<std::string>(common_properties::_TriposAtomType) ==
    // RDKit✔️✔️:                   "O.co2") {
    // RDKit✔️✔️:             Bond *ob = res->getBondBetweenAtoms(nbr->getIdx(), onbr->getIdx());
    // RDKit✔️✔️:             ob->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:             ob->setIsAromatic(false);
    // RDKit✔️✔️:             onbr->setFormalCharge(-1);
    // RDKit✔️✔️:             onbr->setIsAromatic(false);
    // RDKit✔️✔️:             isFixed[onbr->getIdx()] = 1;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         nbr->setIsAromatic(false);
    // RDKit✔️✔️:         isFixed[nbr->getIdx()] = 1;
    // RDKit✔️✔️:
    // RDKit✔️✔️:       } else if (tATT == "C.2" || tATT == "S.o2") {
    // RDKit✔️✔️:         // carboxylates and sulfonates
    // RDKit✔️✔️:         // this should return only the bond between C.2 and O.co2
    // RDKit✔️✔️:         Bond *b = res->getBondBetweenAtoms(idx, nbr->getIdx());
    // RDKit✔️✔️:         if (!isFixed[nbr->getIdx()]) {
    // RDKit✔️✔️:           // the first occurrence is negatively charged and has a single bond
    // RDKit✔️✔️:           b->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:           b->setIsAromatic(false);
    // RDKit✔️✔️:           at->setFormalCharge(-1);
    // RDKit✔️✔️:           at->setIsAromatic(false);
    // RDKit✔️✔️:           nbr->setIsAromatic(false);
    // RDKit✔️✔️:           isFixed[idx] = 1;
    // RDKit✔️✔️:           isFixed[nbr->getIdx()] = 1;
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           // the other occurrences are not charged and have a double bond
    // RDKit✔️✔️:           b->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:           b->setIsAromatic(false);
    // RDKit✔️✔️:           at->setIsAromatic(false);
    // RDKit✔️✔️:           isFixed[idx] = 1;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         std::string nm;
    // RDKit✔️✔️:         res->getProp(common_properties::_Name, nm);
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << nm << ": warning - O.co2 with non C.2 or S.o2 neighbor."
    // RDKit✔️✔️:             << std::endl;
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (tAT == "C.cat") {
    // RDKit✔️✔️:       // positively charged guanidinium groups with C.cat
    // RDKit✔️✔️:       // according to Tripos these should only appear in guanidinium groups
    // RDKit❌❌:       // for the structural fix - the last nitrogen with the least number of
    // RDKit❌❌:       // heavy atoms will get the double bond and the positive charge.
    // RDKit❌❌:       // remember : this is not canonical!
    // RDKit✔️✔️:       // first - set the C.cat as fixed
    // RDKit✔️✔️:       isFixed[idx] = 1;
    // RDKit✔️✔️:       ROMol::ADJ_ITER nbrIdxIt, endNbrsIdxIt, tmpIdxIt;
    // RDKit❌❌:       unsigned int lowestDeg = 100;
    // RDKit✔️✔️:       boost::tie(nbrIdxIt, endNbrsIdxIt) = res->getAtomNeighbors(at);
    // RDKit✔️✔️:       // one problem of programs like Corina is, that they will create also
    // RDKit✔️✔️:       // C.cat
    // RDKit✔️✔️:       // for groups that are not guanidinium. We cannot fix all, but the charged
    // RDKit✔️✔️:       // amidine
    // RDKit✔️✔️:       // in a ring is taken care of too.
    // RDKit✔️✔️:       tmpIdxIt = nbrIdxIt;
    // RDKit✔️✔️:       // declare and initialise toModIdx
    // RDKit✔️✔️:       int toModIdx = -1;
    // RDKit✔️✔️:       unsigned int noNNeighbors = 0;
    // RDKit✔️✔️:       while (tmpIdxIt != endNbrsIdxIt) {
    // RDKit✔️✔️:         if (res->getAtomWithIdx(*tmpIdxIt)->getSymbol() == "N") {
    // RDKit✔️✔️:           ++noNNeighbors;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         ++tmpIdxIt;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (noNNeighbors < 2 || noNNeighbors > 3) {
    // RDKit✔️✔️:         std::string nm;
    // RDKit✔️✔️:         res->getProp(common_properties::_Name, nm);
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << nm << ": Error - C.Cat with bad number of N neighbors."
    // RDKit✔️✔️:             << std::endl;
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       } else if (noNNeighbors == 2) {
    // RDKit✔️✔️:         // the idea is that we assign the positive charge according to the
    // RDKit✔️✔️:         // following precedence:
    // RDKit✔️✔️:         // 1. is part of N-oxide
    // RDKit✔️✔️:         // 2. atom with highest number of hydrogen atoms
    // RDKit✔️✔️:         // 3. atom in ring
    // RDKit✔️✔️:         // 4. random
    // RDKit✔️✔️:         // first we identify the N atoms
    // RDKit✔️✔️:         ROMol::ADJ_ITER idxIt1 = nbrIdxIt, idxIt2 = nbrIdxIt;
    // RDKit✔️✔️:         bool firstIdent = false;
    // RDKit✔️✔️:         while (nbrIdxIt != endNbrsIdxIt) {
    // RDKit✔️✔️:           if (res->getAtomWithIdx(*nbrIdxIt)->getSymbol() == "N") {
    // RDKit✔️✔️:             // fix the bond to one - only the modified N will have a double bond
    // RDKit✔️✔️:             // to C.cat
    // RDKit✔️✔️:             res->getBondBetweenAtoms(idx, *nbrIdxIt)->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:             res->getBondBetweenAtoms(idx, *nbrIdxIt)->setIsAromatic(false);
    // RDKit✔️✔️:             res->getAtomWithIdx(*nbrIdxIt)->setIsAromatic(false);
    // RDKit✔️✔️:             // FIX: what is happening if we hit an atom that was fixed before -
    // RDKit✔️✔️:             // probably nothing.
    // RDKit✔️✔️:             // since I cannot think of a case where this is a problem - throw a
    // RDKit✔️✔️:             // warning
    // RDKit✔️✔️:             if (isFixed[*nbrIdxIt]) {
    // RDKit✔️✔️:               std::string nm;
    // RDKit✔️✔️:               res->getProp(common_properties::_Name, nm);
    // RDKit✔️✔️:               BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:                   << nm << ": warning - charged amidine and isFixed atom."
    // RDKit✔️✔️:                   << std::endl;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:             isFixed[*nbrIdxIt] = 1;
    // RDKit✔️✔️:             if (firstIdent) {
    // RDKit✔️✔️:               idxIt2 = nbrIdxIt;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               idxIt1 = nbrIdxIt;
    // RDKit✔️✔️:               firstIdent = true;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           ++nbrIdxIt;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         // now that we know which are the relevant atoms we check the above
    // RDKit✔️✔️:         // features
    // RDKit✔️✔️:         // is part of N-oxide?
    // RDKit✔️✔️:         // number of hydrogens on each neighbour
    // RDKit✔️✔️:         unsigned int noHNbrs1 = chkNoHNeighbNOx(res, idxIt1, toModIdx);
    // RDKit✔️✔️:         unsigned int noHNbrs2 = chkNoHNeighbNOx(res, idxIt2, toModIdx);
    // RDKit✔️✔️:         if (toModIdx < 0) {
    // RDKit✔️✔️:           // no N-oxide
    // RDKit✔️✔️:           if (noHNbrs1 != noHNbrs2) {
    // RDKit✔️✔️:             if (noHNbrs1 > noHNbrs2) {
    // RDKit✔️✔️:               toModIdx = *idxIt1;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               toModIdx = *idxIt2;  // this is random if both have the same
    // RDKit✔️✔️:                                      // number of atoms
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else {
    // RDKit✔️✔️:             // perceive the rings
    // RDKit✔️✔️:             if (!res->getRingInfo()->isSssrOrBetter()) {
    // RDKit✔️✔️:               MolOps::findSSSR(*res);
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:             // then we check if both atoms are in a ring
    // RDKit✔️✔️:             unsigned int rIdx1 = res->getRingInfo()->numAtomRings((*idxIt1));
    // RDKit✔️✔️:             unsigned int rIdx2 = res->getRingInfo()->numAtomRings((*idxIt2));
    // RDKit✔️✔️:             if (rIdx1 > rIdx2) {
    // RDKit✔️✔️:               toModIdx = *idxIt1;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               toModIdx = *idxIt2;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         res->getBondBetweenAtoms(idx, toModIdx)->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:         res->getBondBetweenAtoms(idx, toModIdx)->setIsAromatic(false);
    // RDKit✔️✔️:         res->getAtomWithIdx(toModIdx)->setFormalCharge(1);
    // RDKit✔️✔️:         at->setIsAromatic(false);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         while (nbrIdxIt != endNbrsIdxIt) {
    // RDKit✔️✔️:           if (!isFixed[*nbrIdxIt]) {
    // RDKit✔️✔️:             // we get in here if this N.pl3 was not seen / fixed before
    // RDKit✔️✔️:             Atom *nbr = res->getAtomWithIdx(*nbrIdxIt);
    // RDKit✔️✔️:             // get the number of heavy atoms connected to this atom
    // RDKit✔️✔️:             ROMol::ADJ_ITER nbrNbrIdxIt, nbrEndNbrsIdxIt;
    // RDKit✔️✔️:             unsigned int hvyAtDeg = 0;
    // RDKit✔️✔️:             boost::tie(nbrNbrIdxIt, nbrEndNbrsIdxIt) =
    // RDKit✔️✔️:                 res->getAtomNeighbors(nbr);
    // RDKit✔️✔️:             while (nbrNbrIdxIt != nbrEndNbrsIdxIt) {
    // RDKit✔️✔️:               if (res->getAtomWithIdx(*nbrNbrIdxIt)->getAtomicNum() > 1) {
    // RDKit✔️✔️:                 std::string nbrAT;
    // RDKit✔️✔️:                 res->getAtomWithIdx(*nbrNbrIdxIt)
    // RDKit✔️✔️:                     ->getProp(common_properties::_TriposAtomType, nbrAT);
    // RDKit✔️✔️:                 if (nbrAT == "C.cat") {
    // RDKit✔️✔️:                   hvyAtDeg += 2;  // that way we reduce the risk of ionising the
    // RDKit✔️✔️:                                   // N attached to another C.cat ...
    // RDKit✔️✔️:                 } else {
    // RDKit✔️✔️:                   ++hvyAtDeg;
    // RDKit✔️✔️:                 }
    // RDKit✔️✔️:               }
    // RDKit✔️✔️:               ++nbrNbrIdxIt;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:             // now check for lowest heavy atom degree
    // RDKit✔️✔️:             if (hvyAtDeg < lowestDeg) {
    // RDKit✔️✔️:               toModIdx = *nbrIdxIt;
    // RDKit✔️✔️:               lowestDeg = hvyAtDeg;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:             // modify the bond between C.Cat and the N.pl3
    // RDKit✔️✔️:             Bond *b = res->getBondBetweenAtoms(idx, *nbrIdxIt);
    // RDKit✔️✔️:             b->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:             b->setIsAromatic(false);
    // RDKit✔️✔️:             nbr->setIsAromatic(false);
    // RDKit✔️✔️:             // set N.pl3 as fixed
    // RDKit✔️✔️:             isFixed[*nbrIdxIt] = 1;
    // RDKit✔️✔️:           } else {
    // RDKit✔️✔️:             // the N is already fixed - since we don't touch this atom make the
    // RDKit✔️✔️:             // bond to single
    // RDKit✔️✔️:             // FIX: check on 3-way symmetric guanidinium mol -
    // RDKit✔️✔️:             //     this could produce a only single bonded C.cat for bad H mols
    // RDKit✔️✔️:             res->getBondBetweenAtoms(idx, *nbrIdxIt)->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:             res->getBondBetweenAtoms(idx, *nbrIdxIt)->setIsAromatic(false);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           ++nbrIdxIt;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         // now modify the respective N and the C.cat
    // RDKit✔️✔️:         Bond *b = res->getBondBetweenAtoms(idx, toModIdx);
    // RDKit✔️✔️:         b->setBondType(Bond::DOUBLE);
    // RDKit✔️✔️:         b->setIsAromatic(false);
    // RDKit✔️✔️:         res->getAtomWithIdx(toModIdx)->setFormalCharge(1);
    // RDKit✔️✔️:         at->setIsAromatic(false);
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: cleanUpMol2Substructures
    let mut is_fixed = vec![false; builder.atoms().len()];

    for idx in 0..builder.atoms().len() {
        if is_fixed[idx] {
            continue;
        }
        let atom_id = AtomId::new(idx);
        let tripos_atom_type = builder.atoms()[idx]
            .prop("_TriposAtomType")
            .ok_or_else(|| Mol2ReadError::Parse("Missing _TriposAtomType".to_string()))?
            .to_string();

        if tripos_atom_type == "N.4" {
            let Some(atom) = builder.atom_mut(atom_id) else {
                return Err(Mol2ReadError::Parse("index mismatch".to_string()));
            };
            atom.set_formal_charge(1);
        } else if tripos_atom_type == "O.co2" {
            if builder.degree(atom_id) != 1 {
                return Ok(false);
            }

            let bond_id = *builder
                .neighbor_bonds(atom_id)
                .first()
                .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_string()))?;
            let bond = builder
                .bond(bond_id)
                .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_string()))?;
            let nbr_id = if bond.begin() == atom_id {
                bond.end()
            } else {
                bond.begin()
            };
            let nbr_type = builder.atoms()[nbr_id.index()]
                .prop("_TriposAtomType")
                .ok_or_else(|| Mol2ReadError::Parse("Missing _TriposAtomType".to_string()))?
                .to_string();

            if nbr_type == "P.3" {
                let Some(bond) = builder.bond_mut(bond_id) else {
                    return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                };
                bond.set_order(BondOrder::Double);
                bond.set_aromatic(false);
                let Some(atom) = builder.atom_mut(atom_id) else {
                    return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                };
                atom.set_aromatic(false);
                is_fixed[idx] = true;

                let onbr_ids = builder.neighbor_bonds(nbr_id).to_vec();
                for onbr_bond_id in onbr_ids {
                    let onbr_bond = builder
                        .bond(onbr_bond_id)
                        .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_string()))?;
                    let onbr_id = if onbr_bond.begin() == nbr_id {
                        onbr_bond.end()
                    } else {
                        onbr_bond.begin()
                    };
                    if builder.atoms()[onbr_id.index()].atomic_number() == 8
                        && !is_fixed[onbr_id.index()]
                        && builder.atoms()[onbr_id.index()].prop("_TriposAtomType") == Some("O.co2")
                    {
                        let Some(bond) = builder.bond_mut(onbr_bond_id) else {
                            return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                        };
                        bond.set_order(BondOrder::Single);
                        bond.set_aromatic(false);
                        let Some(atom) = builder.atom_mut(onbr_id) else {
                            return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                        };
                        atom.set_formal_charge(-1);
                        atom.set_aromatic(false);
                        is_fixed[onbr_id.index()] = true;
                    }
                }

                let Some(atom) = builder.atom_mut(nbr_id) else {
                    return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                };
                atom.set_aromatic(false);
                is_fixed[nbr_id.index()] = true;
            } else if nbr_type == "C.2" || nbr_type == "S.o2" {
                if !is_fixed[nbr_id.index()] {
                    let Some(bond) = builder.bond_mut(bond_id) else {
                        return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                    };
                    bond.set_order(BondOrder::Single);
                    bond.set_aromatic(false);
                    let Some(atom) = builder.atom_mut(atom_id) else {
                        return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                    };
                    atom.set_formal_charge(-1);
                    atom.set_aromatic(false);
                    let Some(atom) = builder.atom_mut(nbr_id) else {
                        return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                    };
                    atom.set_aromatic(false);
                    is_fixed[idx] = true;
                    is_fixed[nbr_id.index()] = true;
                } else {
                    let Some(bond) = builder.bond_mut(bond_id) else {
                        return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                    };
                    bond.set_order(BondOrder::Double);
                    bond.set_aromatic(false);
                    let Some(atom) = builder.atom_mut(atom_id) else {
                        return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                    };
                    atom.set_aromatic(false);
                    is_fixed[idx] = true;
                }
            } else {
                return Ok(false);
            }
        } else if tripos_atom_type == "C.cat" {
            is_fixed[idx] = true;

            let mut neighbor_pairs = Vec::new();
            let mut nitrogen_neighbors = Vec::new();
            for bond_id in builder.neighbor_bonds(atom_id).iter().copied() {
                let bond = builder
                    .bond(bond_id)
                    .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_string()))?;
                let neighbor_id = if bond.begin() == atom_id {
                    bond.end()
                } else {
                    bond.begin()
                };
                neighbor_pairs.push((neighbor_id, bond_id));
                if builder.atoms()[neighbor_id.index()].atomic_number() == 7 {
                    nitrogen_neighbors.push((neighbor_id, bond_id));
                }
            }

            if nitrogen_neighbors.len() < 2 || nitrogen_neighbors.len() > 3 {
                return Ok(false);
            }
            if nitrogen_neighbors.len() == 3 {
                let mut lowest_degree = 100_u32;
                let mut to_mod_idx = None;
                let mut to_mod_bond = None;
                for (neighbor_id, bond_id) in neighbor_pairs.iter().copied() {
                    if !is_fixed[neighbor_id.index()] {
                        let mut heavy_atom_degree = 0_u32;
                        for neighbor_bond_id in builder.neighbor_bonds(neighbor_id).iter().copied()
                        {
                            let neighbor_bond =
                                builder.bond(neighbor_bond_id).ok_or_else(|| {
                                    Mol2ReadError::Parse("index mismatch".to_string())
                                })?;
                            let neighbor_neighbor_id = if neighbor_bond.begin() == neighbor_id {
                                neighbor_bond.end()
                            } else {
                                neighbor_bond.begin()
                            };
                            if builder.atoms()[neighbor_neighbor_id.index()].atomic_number() > 1 {
                                if builder.atoms()[neighbor_neighbor_id.index()]
                                    .prop("_TriposAtomType")
                                    == Some("C.cat")
                                {
                                    heavy_atom_degree = heavy_atom_degree.wrapping_add(2);
                                } else {
                                    heavy_atom_degree = heavy_atom_degree.wrapping_add(1);
                                }
                            }
                        }
                        if heavy_atom_degree < lowest_degree {
                            to_mod_idx = Some(neighbor_id);
                            to_mod_bond = Some(bond_id);
                            lowest_degree = heavy_atom_degree;
                        }
                        let Some(bond) = builder.bond_mut(bond_id) else {
                            return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                        };
                        bond.set_order(BondOrder::Single);
                        bond.set_aromatic(false);
                        let Some(atom) = builder.atom_mut(neighbor_id) else {
                            return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                        };
                        atom.set_aromatic(false);
                        is_fixed[neighbor_id.index()] = true;
                    } else {
                        let Some(bond) = builder.bond_mut(bond_id) else {
                            return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                        };
                        bond.set_order(BondOrder::Single);
                        bond.set_aromatic(false);
                    }
                }
                let to_mod_idx =
                    to_mod_idx.ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_string()))?;
                let to_mod_bond = to_mod_bond
                    .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_string()))?;
                let Some(bond) = builder.bond_mut(to_mod_bond) else {
                    return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                };
                bond.set_order(BondOrder::Double);
                bond.set_aromatic(false);
                let Some(atom) = builder.atom_mut(to_mod_idx) else {
                    return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                };
                atom.set_formal_charge(1);
                let Some(atom) = builder.atom_mut(atom_id) else {
                    return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                };
                atom.set_aromatic(false);
            } else {
                let (first_n_id, first_bond_id) = nitrogen_neighbors[0];
                let (second_n_id, second_bond_id) = nitrogen_neighbors[1];
                for (nitrogen_id, bond_id) in nitrogen_neighbors.iter().copied() {
                    let Some(bond) = builder.bond_mut(bond_id) else {
                        return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                    };
                    bond.set_order(BondOrder::Single);
                    bond.set_aromatic(false);
                    let Some(atom) = builder.atom_mut(nitrogen_id) else {
                        return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                    };
                    atom.set_aromatic(false);
                    is_fixed[nitrogen_id.index()] = true;
                }

                let mut to_mod_idx = None;
                let no_h_neighbors_1 =
                    check_no_h_neighbors_n_oxide_like_rdkit(builder, first_n_id, &mut to_mod_idx)?;
                let no_h_neighbors_2 =
                    check_no_h_neighbors_n_oxide_like_rdkit(builder, second_n_id, &mut to_mod_idx)?;
                let to_mod_id = if let Some(to_mod_idx) = to_mod_idx {
                    to_mod_idx
                } else if no_h_neighbors_1 != no_h_neighbors_2 {
                    if no_h_neighbors_1 > no_h_neighbors_2 {
                        first_n_id
                    } else {
                        second_n_id
                    }
                } else {
                    let adjacency =
                        AdjacencyList::try_from_topology(builder.atoms().len(), builder.bonds())
                            .map_err(|err| Mol2ReadError::Parse(err.to_string()))?;
                    let ring_info =
                        find_sssr_from_parts(builder.atoms().len(), builder.bonds(), &adjacency)
                            .map_err(|err| Mol2ReadError::Parse(err.to_string()))?;
                    let first_ring_count = ring_info.num_atom_rings(first_n_id);
                    let second_ring_count = ring_info.num_atom_rings(second_n_id);
                    if first_ring_count > second_ring_count {
                        first_n_id
                    } else {
                        second_n_id
                    }
                };
                let to_mod_bond_id = if to_mod_id == first_n_id {
                    first_bond_id
                } else if to_mod_id == second_n_id {
                    second_bond_id
                } else {
                    return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                };
                let Some(bond) = builder.bond_mut(to_mod_bond_id) else {
                    return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                };
                bond.set_order(BondOrder::Double);
                bond.set_aromatic(false);
                let Some(atom) = builder.atom_mut(to_mod_id) else {
                    return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                };
                atom.set_formal_charge(1);
                let Some(atom) = builder.atom_mut(atom_id) else {
                    return Err(Mol2ReadError::Parse("index mismatch".to_string()));
                };
                atom.set_aromatic(false);
            }
        }
    }

    Ok(true)
}

fn scan_mol2_sections_like_rdkit(input: &str) -> Result<Mol2SectionOffsets, Mol2ReadError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2DataStream section scan
    // RDKit✔️✔️: std::string tempStr, lineBeg;
    // RDKit✔️✔️: typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
    // RDKit✔️✔️: boost::char_separator<char> sep(" \t\n");
    // RDKit✔️✔️: Utils::LocaleSwitcher ls;
    // RDKit✔️✔️: std::streampos molStart = 0, atomStart = 0, bondStart = 0, chargeStart = 0;
    // RDKit✔️✔️: while (!inStream.eof() && !inStream.fail()) {
    // RDKit✔️✔️:   tempStr = getLine(inStream);
    // RDKit✔️✔️:   if (inStream.eof()) {
    // RDKit✔️✔️:     break;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (tempStr != "" && tempStr[0] == '@') {
    // RDKit✔️✔️:     tokenizer tokens(tempStr, sep);
    // RDKit✔️✔️:     std::string firstToken = *tokens.begin();
    // RDKit✔️✔️:     if (firstToken == "@<TRIPOS>MOLECULE") {
    // RDKit✔️✔️:       if (!molStart) {
    // RDKit✔️✔️:         molStart = inStream.tellg();
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (firstToken == "@<TRIPOS>ATOM") {
    // RDKit✔️✔️:       atomStart = inStream.tellg();
    // RDKit✔️✔️:     } else if (firstToken == "@<TRIPOS>BOND") {
    // RDKit✔️✔️:       bondStart = inStream.tellg();
    // RDKit✔️✔️:     } else if (firstToken == "@<TRIPOS>UNITY_ATOM_ATTR") {
    // RDKit✔️✔️:       chargeStart = inStream.tellg();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (!molStart) {
    // RDKit✔️✔️:   throw FileParseException("No MOLECULE block found in Mol2 data");
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (!atomStart) {
    // RDKit✔️✔️:   throw FileParseException("No ATOM block found in Mol2 data");
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2DataStream section scan
    let mut molecule_start = None;
    let mut atom_start = None;
    let mut bond_start = None;
    let mut charge_start = None;
    let mut offset = 0usize;

    for raw_line in input.split_inclusive('\n') {
        let end_offset = offset + raw_line.len();
        let temp_str = raw_line.trim_end_matches(['\r', '\n']);

        if temp_str.starts_with('@') {
            let first_token = temp_str.split_whitespace().next().unwrap_or("");
            match first_token {
                "@<TRIPOS>MOLECULE" => {
                    if molecule_start.is_none() {
                        molecule_start = Some(end_offset);
                    } else {
                        break;
                    }
                }
                "@<TRIPOS>ATOM" => atom_start = Some(end_offset),
                "@<TRIPOS>BOND" => bond_start = Some(end_offset),
                "@<TRIPOS>UNITY_ATOM_ATTR" => charge_start = Some(end_offset),
                _ => {}
            }
        }

        offset = end_offset;
    }

    let molecule_start = molecule_start
        .ok_or_else(|| Mol2ReadError::Parse("No MOLECULE block found in Mol2 data".to_string()))?;
    let atom_start = atom_start
        .ok_or_else(|| Mol2ReadError::Parse("No ATOM block found in Mol2 data".to_string()))?;

    Ok(Mol2SectionOffsets {
        molecule_start,
        atom_start,
        bond_start,
        charge_start,
    })
}

pub fn mol_from_mol2_data_stream_like_rdkit(
    input: &str,
    params: Mol2ReadParams,
) -> Result<Option<Mol2Record>, Mol2ReadError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2DataStream unsanitized assembly
    // RDKit✔️✔️:   if (inStream.eof()) {
    // RDKit✔️✔️:     inStream.clear();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   inStream.seekg(molStart, std::ios::beg);
    // RDKit✔️✔️:   tempStr = getLine(inStream);
    // RDKit✔️✔️:   auto res = std::make_unique<RWMol>();
    // RDKit✔️✔️:   boost::trim_right(tempStr);
    // RDKit✔️✔️:   res->setProp(common_properties::_Name, tempStr);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   tempStr = getLine(inStream);
    // RDKit✔️✔️:   tokenizer tokens(tempStr, sep);
    // RDKit✔️✔️:   if (tokens.begin() == tokens.end()) {
    // RDKit✔️✔️:     throw FileParseException("Empty counts line");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int nAtoms = 0, nBonds = 0;
    // RDKit✔️✔️:   tokenizer::const_iterator itemIt = tokens.begin();
    // RDKit✔️✔️:   // counts line, this is where we really get started
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     nAtoms = boost::lexical_cast<unsigned int>(*itemIt);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     ++itemIt;
    // RDKit✔️✔️:     if (itemIt != tokens.end()) {
    // RDKit✔️✔️:       nBonds = boost::lexical_cast<unsigned int>(*itemIt);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️✔️:     std::ostringstream errout;
    // RDKit✔️✔️:     errout << "Cannot convert " << *itemIt << " to unsigned int";
    // RDKit✔️✔️:     throw FileParseException(errout.str());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (nAtoms == 0) {
    // RDKit✔️✔️:     throw FileParseException("molecule has no atoms");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   tempStr = getLine(inStream);  // mol_type - ignore
    // RDKit✔️✔️:   tempStr = getLine(inStream);
    // RDKit✔️✔️:   boost::trim(tempStr);
    // RDKit✔️✔️:   res->setProp("_TriposChargeType", tempStr);
    // RDKit✔️✔️:   // stop here since we don't support anything else from the MOLECULE block
    // RDKit✔️✔️:   INT_VECT idxCorresp(nAtoms, -1);
    // RDKit✔️✔️:   inStream.seekg(atomStart, std::ios::beg);
    // RDKit✔️✔️:   ParseMol2AtomBlock(&inStream, res.get(), nAtoms, idxCorresp);
    // RDKit✔️✔️:   if (nBonds) {
    // RDKit✔️✔️:     // stop here since we don't support anything else from the MOLECULE block
    // RDKit✔️✔️:     inStream.seekg(bondStart, std::ios::beg);
    // RDKit✔️✔️:     ParseMol2BondBlock(&inStream, res.get(), nBonds, idxCorresp);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!chargeStart) {
    // RDKit✔️✔️:     bool molFixed;
    // RDKit✔️✔️:     if (params.cleanupSubstructures) {
    // RDKit✔️✔️:       molFixed = cleanUpMol2Substructures(res.get());
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       molFixed = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (!molFixed) {
    // RDKit✔️✔️:       return nullptr;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // mol2 format does not support formal charge information, hence we need to
    // RDKit✔️✔️:     // guess it based on default and explicit valences
    // RDKit✔️✔️:     guessFormalCharges(res.get());
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     inStream.seekg(chargeStart, std::ios::beg);
    // RDKit✔️✔️:     readFormalChargesFromAttr(&inStream, res.get());
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2DataStream unsanitized assembly
    let offsets = scan_mol2_sections_like_rdkit(input)?;
    let header = parse_mol2_molecule_header_like_rdkit(input, offsets.molecule_start)?;
    let atom_block = parse_mol2_atom_block_like_rdkit(input, offsets.atom_start, header.n_atoms)?;

    let mut bond_block = if header.n_bonds == 0 {
        Mol2BondBlock {
            builder: atom_block.builder,
            n_bad_bonds: 0,
        }
    } else {
        let bond_start = offsets
            .bond_start
            .ok_or_else(|| Mol2ReadError::Parse("No BOND block found".to_string()))?;
        parse_mol2_bond_block_like_rdkit(input, bond_start, header.n_bonds, atom_block)?
    };

    bond_block.builder = bond_block
        .builder
        .with_name(header.name)
        .with_property("_TriposChargeType", header.charge_type);

    let builder = if let Some(charge_start) = offsets.charge_start {
        read_formal_charges_from_attr_like_rdkit(input, charge_start, bond_block)?.builder
    } else {
        if params.cleanup_substructures
            && !clean_up_mol2_substructures_like_rdkit(&mut bond_block.builder)?
        {
            return Ok(None);
        }
        guess_formal_charges_like_rdkit(&mut bond_block.builder)?;
        bond_block.builder
    };

    let molecule = builder
        .build()
        .map_err(|err| Mol2ReadError::Parse(err.to_string()))?;
    Ok(Some(Mol2Record { molecule }))
}

fn finish_mol2_read_like_rdkit(
    mut record: Mol2Record,
    params: Mol2ReadParams,
) -> Result<Mol2Record, Mol2ReadError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2DataStream finalization order
    // RDKit✔️✔️:   // set chirality prior to sanitization since it happens from 3D and it's not
    // RDKit✔️✔️:   // possible anymore once the hydrogens are removed
    // RDKit✔️✔️:   // FIX: for now this is only for the first conformer - need to be changed once
    // RDKit✔️✔️:   // we use multiconformer files
    // RDKit✔️✔️:   MolOps::assignChiralTypesFrom3D(*res);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (res && params.sanitize) {
    // RDKit✔️✔️:     MolOps::cleanUp(*res);
    // RDKit✔️✔️:     ...
    // RDKit✔️✔️:     res->updatePropertyCache(false);
    // RDKit✔️✔️:     MolOps::assignStereochemistry(*res, true, true);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2DataStream finalization order
    assign_chiral_types_from_3d_for_mol2_like_rdkit(&mut record.molecule);
    record.molecule = sanitize_mol2_molecule_like_rdkit(record.molecule, params)?;
    Ok(record)
}

pub fn mol_from_mol2_block_like_rdkit(
    mol_block: &str,
    params: Mol2ReadParams,
) -> Result<Option<Mol2Record>, Mol2ReadError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2Block
    // RDKit✔️✔️: std::unique_ptr<RWMol> MolFromMol2Block(const std::string &molBlock,
    // RDKit✔️✔️:                                         const Mol2ParserParams &params) {
    // RDKit✔️✔️:   std::istringstream inStream(molBlock);
    // RDKit✔️✔️:   return MolFromMol2DataStream(inStream, params);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2Block
    mol_from_mol2_data_stream_like_rdkit(mol_block, params)?
        .map(|record| finish_mol2_read_like_rdkit(record, params))
        .transpose()
}

pub fn mol_from_mol2_file_like_rdkit(
    path: impl AsRef<Path>,
    params: Mol2ReadParams,
) -> Result<Option<Mol2Record>, Mol2ReadError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2File
    // RDKit✔️✔️: std::ifstream inStream(fName.c_str(), std::ios_base::binary);
    // RDKit✔️✔️: if (!inStream || (inStream.bad())) {
    // RDKit✔️✔️:   std::ostringstream errout;
    // RDKit✔️✔️:   errout << "Bad input file " << fName;
    // RDKit✔️✔️:   throw BadFileException(errout.str());
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (!inStream.eof()) {
    // RDKit✔️✔️:   return MolFromMol2DataStream(inStream, params);
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   return nullptr;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2File
    let path = path.as_ref();
    let text = fs::read_to_string(path)
        .map_err(|_| Mol2ReadError::Io(format!("Bad input file {}", path.display())))?;
    if text.is_empty() {
        Ok(None)
    } else {
        mol_from_mol2_block_like_rdkit(&text, params)
    }
}

pub fn read_mol2_from_str(s: &str) -> Result<Option<Mol2Record>, Mol2ReadError> {
    read_mol2_from_str_with_params(s, Mol2ReadParams::default())
}

pub fn read_mol2_from_str_with_params(
    s: &str,
    params: Mol2ReadParams,
) -> Result<Option<Mol2Record>, Mol2ReadError> {
    // BEGIN RDKIT CPP INLINE third_party/rdkit/Code/GraphMol/FileParsers/FileParsers.h :: Mol2BlockToMol
    // RDKit✔️✔️: inline RWMol *Mol2BlockToMol(const std::string &molBlock, bool sanitize = true,
    // RDKit✔️✔️:                              bool removeHs = true,
    // RDKit✔️✔️:                              Mol2Type variant = Mol2Type::CORINA,
    // RDKit✔️✔️:                              bool cleanupSubstructures = true) {
    // RDKit✔️✔️:   v2::FileParsers::Mol2ParserParams ps;
    // RDKit✔️✔️:   ps.sanitize = sanitize;
    // RDKit✔️✔️:   ps.removeHs = removeHs;
    // RDKit✔️✔️:   ps.variant = variant;
    // RDKit✔️✔️:   ps.cleanupSubstructures = cleanupSubstructures;
    // RDKit✔️✔️:   return v2::FileParsers::MolFromMol2Block(molBlock, ps).release();
    // RDKit✔️✔️: }
    // END RDKIT CPP INLINE third_party/rdkit/Code/GraphMol/FileParsers/FileParsers.h :: Mol2BlockToMol
    mol_from_mol2_block_like_rdkit(s, params)
}

pub fn read_mol2_file(path: impl AsRef<Path>) -> Result<Option<Mol2Record>, Mol2ReadError> {
    read_mol2_file_with_params(path, Mol2ReadParams::default())
}

pub fn read_mol2_file_with_params(
    path: impl AsRef<Path>,
    params: Mol2ReadParams,
) -> Result<Option<Mol2Record>, Mol2ReadError> {
    // BEGIN RDKIT CPP INLINE third_party/rdkit/Code/GraphMol/FileParsers/FileParsers.h :: Mol2FileToMol
    // RDKit✔️✔️: inline RWMol *Mol2FileToMol(const std::string &fName, bool sanitize = true,
    // RDKit✔️✔️:                             bool removeHs = true,
    // RDKit✔️✔️:                             Mol2Type variant = Mol2Type::CORINA,
    // RDKit✔️✔️:                             bool cleanupSubstructures = true) {
    // RDKit✔️✔️:   v2::FileParsers::Mol2ParserParams ps;
    // RDKit✔️✔️:   ps.sanitize = sanitize;
    // RDKit✔️✔️:   ps.removeHs = removeHs;
    // RDKit✔️✔️:   ps.variant = variant;
    // RDKit✔️✔️:   ps.cleanupSubstructures = cleanupSubstructures;
    // RDKit✔️✔️:   return v2::FileParsers::MolFromMol2File(fName, ps).release();
    // RDKit✔️✔️: }
    // END RDKIT CPP INLINE third_party/rdkit/Code/GraphMol/FileParsers/FileParsers.h :: Mol2FileToMol
    mol_from_mol2_file_like_rdkit(path, params)
}

fn assign_chiral_types_from_3d_for_mol2_like_rdkit(molecule: &mut Molecule) {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2DataStream 3D chirality call
    // RDKit✔️✔️:   // set chirality prior to sanitization since it happens from 3D and it's not
    // RDKit✔️✔️:   // possible anymore once the hydrogens are removed
    // RDKit✔️✔️:   // FIX: for now this is only for the first conformer - need to be changed once
    // RDKit✔️✔️:   // we use multiconformer files
    // RDKit✔️✔️:   MolOps::assignChiralTypesFrom3D(*res);
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2DataStream 3D chirality call
    crate::notation::smiles::assign_chiral_types_from_3d(molecule, 0);
}

fn sanitize_mol2_molecule_like_rdkit(
    molecule: Molecule,
    params: Mol2ReadParams,
) -> Result<Molecule, Mol2ReadError> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2DataStream sanitized removeHs branch
    // RDKit✔️✔️:   if (res && params.sanitize) {
    // RDKit✔️✔️:     MolOps::cleanUp(*res);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     try {
    // RDKit✔️✔️:       // when we sanitize for mol2, we skip the cleanup organometallic step since it's
    // RDKit✔️✔️:       // not really compatible with the semantics of mol2 files
    // RDKit✔️✔️:       constexpr auto sanitizeFlags = MolOps::SanitizeFlags::SANITIZE_ALL ^
    // RDKit✔️✔️:                             MolOps::SanitizeFlags::SANITIZE_CLEANUP_ORGANOMETALLICS;
    // RDKit✔️✔️:       if (params.removeHs) {
    // RDKit✔️✔️:         // Bond stereo detection must happen before H removal, or
    // RDKit✔️✔️:         // else we might be removing stereogenic H atoms in double
    // RDKit✔️✔️:         // bonds (e.g. imines). But before we run stereo detection,
    // RDKit✔️✔️:         // we need to run mol cleanup so don't have trouble with
    // RDKit✔️✔️:         // e.g. nitro groups. Sadly, this a;; means we will find
    // RDKit✔️✔️:         // run both cleanup and ring finding twice (a fast find
    // RDKit✔️✔️:         // rings in bond stereo detection, and another in
    // RDKit✔️✔️:         // sanitization's SSSR symmetrization).
    // RDKit✔️✔️:         unsigned int failedOp = 0;
    // RDKit✔️✔️:         MolOps::sanitizeMol(*res, failedOp, MolOps::SanitizeFlags::SANITIZE_CLEANUP);
    // RDKit✔️✔️:         MolOps::detectBondStereochemistry(*res);
    // RDKit✔️✔️:         MolOps::RemoveHsParameters rhp;
    // RDKit✔️✔️:         bool sanitize = false;
    // RDKit✔️✔️:         MolOps::removeHs(*res, rhp, sanitize);
    // RDKit✔️✔️:         MolOps::sanitizeMol(*res, failedOp, sanitizeFlags);
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2DataStream sanitized removeHs branch
    if !params.sanitize {
        return Ok(molecule);
    }
    let mut molecule = molecule
        .sanitize_with_ops(crate::SanitizeOps::CLEANUP)
        .map_err(mol2_operation_error)?;
    if params.remove_hs {
        if let Some(conf_id) = molecule.conformers_3d().first().map(|conf| conf.id()) {
            crate::notation::smiles::set_double_bond_neighbor_directions(&mut molecule, conf_id)
                .map_err(|err| {
                    Mol2ReadError::Parse(format!("bond stereo detection failed: {err}"))
                })?;
        }
        molecule = molecule
            .without_hydrogens_with_sanitize(false)
            .map_err(mol2_operation_error)?;
        molecule = molecule
            .sanitize_with_ops(mol2_sanitize_flags_like_rdkit())
            .map_err(mol2_operation_error)?;
    } else {
        // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2DataStream sanitized keepHs branch
        // RDKit✔️✔️:       } else {
        // RDKit✔️✔️:         unsigned int failedOp;
        // RDKit✔️✔️:         MolOps::sanitizeMol(*res, failedOp, sanitizeFlags);
        // RDKit✔️✔️:         MolOps::detectBondStereochemistry(*res);
        // RDKit✔️✔️:       }
        // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2DataStream sanitized keepHs branch
        molecule = molecule
            .sanitize_with_ops(mol2_sanitize_flags_like_rdkit())
            .map_err(mol2_operation_error)?;
        if let Some(conf_id) = molecule.conformers_3d().first().map(|conf| conf.id()) {
            crate::notation::smiles::set_double_bond_neighbor_directions(&mut molecule, conf_id)
                .map_err(|err| {
                    Mol2ReadError::Parse(format!("bond stereo detection failed: {err}"))
                })?;
        }
    }
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2DataStream sanitized finalization
    // RDKit✔️✔️:     res->updatePropertyCache(false);
    // RDKit✔️✔️:     MolOps::assignStereochemistry(*res, true, true);
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp :: MolFromMol2DataStream sanitized finalization
    crate::notation::smiles::assign_stereochemistry_cleanup_subset(&mut molecule, true)
        .map_err(|err| Mol2ReadError::Parse(format!("stereochemistry assignment failed: {err}")))?;
    Ok(molecule)
}

fn mol2_sanitize_flags_like_rdkit() -> crate::SanitizeOps {
    crate::SanitizeOps::CLEANUP
        | crate::SanitizeOps::PROPERTIES
        | crate::SanitizeOps::SYMMRINGS
        | crate::SanitizeOps::KEKULIZE
        | crate::SanitizeOps::FIND_RADICALS
        | crate::SanitizeOps::SET_AROMATICITY
        | crate::SanitizeOps::SET_CONJUGATION
        | crate::SanitizeOps::SET_HYBRIDIZATION
        | crate::SanitizeOps::CLEANUP_CHIRALITY
        | crate::SanitizeOps::ADJUST_HYDROGENS
        | crate::SanitizeOps::CLEANUP_ATROPISOMERS
}

fn mol2_operation_error(error: crate::OperationError) -> Mol2ReadError {
    match error {
        crate::OperationError::UnsupportedFeature { source, .. } => {
            Mol2ReadError::UnsupportedFeature(source)
        }
        other => Mol2ReadError::Parse(other.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_header_from_block(input: &str) -> Result<Mol2MoleculeHeader, Mol2ReadError> {
        let offsets = scan_mol2_sections_like_rdkit(input)?;
        parse_mol2_molecule_header_like_rdkit(input, offsets.molecule_start)
    }

    fn read_unity_atom_attr_from_block(
        input: &str,
        n_atoms: u32,
        n_bonds: u32,
    ) -> Result<Mol2UnityAtomAttrBlock, Mol2ReadError> {
        let offsets = scan_mol2_sections_like_rdkit(input)?;
        let atom_block = parse_mol2_atom_block_like_rdkit(input, offsets.atom_start, n_atoms)?;
        let bond_block = if n_bonds == 0 {
            Mol2BondBlock {
                builder: atom_block.builder,
                n_bad_bonds: 0,
            }
        } else {
            parse_mol2_bond_block_like_rdkit(
                input,
                offsets
                    .bond_start
                    .ok_or_else(|| Mol2ReadError::Parse("No BOND block found".to_string()))?,
                n_bonds,
                atom_block,
            )?
        };
        read_formal_charges_from_attr_like_rdkit(
            input,
            offsets.charge_start.ok_or_else(|| {
                Mol2ReadError::Parse("No UNITY_ATOM_ATTR block found".to_string())
            })?,
            bond_block,
        )
    }

    fn parse_bond_block_from_block(
        input: &str,
        n_atoms: u32,
        n_bonds: u32,
    ) -> Result<Mol2BondBlock, Mol2ReadError> {
        let offsets = scan_mol2_sections_like_rdkit(input)?;
        let atom_block = parse_mol2_atom_block_like_rdkit(input, offsets.atom_start, n_atoms)?;
        parse_mol2_bond_block_like_rdkit(
            input,
            offsets
                .bond_start
                .ok_or_else(|| Mol2ReadError::Parse("No BOND block found".to_string()))?,
            n_bonds,
            atom_block,
        )
    }

    fn read_unsanitized_mol2_from_block(input: &str) -> Result<Option<Mol2Record>, Mol2ReadError> {
        mol_from_mol2_data_stream_like_rdkit(
            input,
            Mol2ReadParams {
                sanitize: false,
                remove_hs: false,
                variant: Mol2Type::Corina,
                cleanup_substructures: true,
            },
        )
    }

    #[test]
    fn mol2_header_trims_name_right_only_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\n  named mol \t  \n1 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n";
        let header = parse_header_from_block(input).unwrap();

        assert_eq!(header.name, "  named mol");
    }

    #[test]
    fn mol2_header_rejects_empty_counts_line_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n \t \nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n";
        let err = parse_header_from_block(input).unwrap_err();

        assert_eq!(err.to_string(), "MOL2 parse failed: Empty counts line");
    }

    #[test]
    fn mol2_header_rejects_invalid_atom_count_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\nabc 1\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n";
        let err = parse_header_from_block(input).unwrap_err();

        assert_eq!(
            err.to_string(),
            "MOL2 parse failed: Cannot convert abc to unsigned int"
        );
    }

    #[test]
    fn mol2_header_defaults_missing_bond_count_to_zero_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n3\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n";
        let header = parse_header_from_block(input).unwrap();

        assert_eq!(header.n_atoms, 3);
        assert_eq!(header.n_bonds, 0);
    }

    #[test]
    fn mol2_header_rejects_zero_atoms_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n";
        let err = parse_header_from_block(input).unwrap_err();

        assert_eq!(err.to_string(), "MOL2 parse failed: molecule has no atoms");
    }

    #[test]
    fn mol2_header_ignores_molecule_type_line_like_rdkit() {
        let input =
            "@<TRIPOS>MOLECULE\nmol\n2 1\nUNRECOGNIZED_MOL_TYPE\nGASTEIGER\n@<TRIPOS>ATOM\n";
        let header = parse_header_from_block(input).unwrap();

        assert_eq!(header.n_atoms, 2);
        assert_eq!(header.n_bonds, 1);
        assert_eq!(header.charge_type, "GASTEIGER");
    }

    #[test]
    fn mol2_header_trims_and_retains_charge_type_property_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n1 0\nSMALL\n \t USER_CHARGES \t \n@<TRIPOS>ATOM\n";
        let header = parse_header_from_block(input).unwrap();

        assert_eq!(header.charge_type, "USER_CHARGES");
    }

    #[test]
    fn mol2_atom_line_reads_required_fields_like_rdkit() {
        let atom = parse_mol2_atom_line_like_rdkit("17 C7 -1.25 2.5 0 C.ar")
            .unwrap()
            .unwrap();

        assert_eq!(atom.atom_name, "C7");
        assert_eq!(atom.position, [-1.25, 2.5, 0.0]);
        assert_eq!(atom.sybyl_atom_type, "C.ar");
        assert_eq!(atom.sybyl_symbol, "C");
        assert!(atom.no_implicit);
        assert_eq!(atom.partial_charge, None);
    }

    #[test]
    fn mol2_atom_line_rejects_coordinate_errors_like_rdkit() {
        let err = parse_mol2_atom_line_like_rdkit("1 C1 not-a-number 0 0 C.3").unwrap_err();

        assert_eq!(
            err.to_string(),
            "MOL2 parse failed: Cannot process mol2 coordinates."
        );
    }

    #[test]
    fn mol2_atom_line_retains_tripos_atom_name_like_rdkit() {
        let atom = parse_mol2_atom_line_like_rdkit("1 Zn42 0 0 0 Zn")
            .unwrap()
            .unwrap();

        assert_eq!(atom.atom_name, "Zn42");
        assert_eq!(atom.atom_spec.prop("_TriposAtomName"), Some("Zn42"));
    }

    #[test]
    fn mol2_atom_line_retains_tripos_atom_type_like_rdkit() {
        let atom = parse_mol2_atom_line_like_rdkit("1 N1 0 0 0 N.am")
            .unwrap()
            .unwrap();

        assert_eq!(atom.sybyl_atom_type, "N.am");
        assert_eq!(atom.sybyl_symbol, "N");
        assert_eq!(atom.atom_spec.prop("_TriposAtomType"), Some("N.am"));
    }

    #[test]
    fn mol2_atom_line_retains_tripos_partial_charge_like_rdkit() {
        let atom = parse_mol2_atom_line_like_rdkit("1 O1 0 0 0 O.2 12 RES -0.55 BACKBONE")
            .unwrap()
            .unwrap();

        assert_eq!(atom.partial_charge.as_deref(), Some("-0.55"));
        assert_eq!(atom.atom_spec.prop("_TriposPartialCharge"), Some("-0.55"));
    }

    #[test]
    fn mol2_atom_line_sets_no_implicit_hydrogen_state_like_rdkit() {
        let atom = parse_mol2_atom_line_like_rdkit("1 H1 0 0 0 H")
            .unwrap()
            .unwrap();

        assert!(atom.no_implicit);
        assert!(atom.atom_spec.no_implicit());
    }

    #[test]
    fn mol2_atom_type_removes_lp_atoms_like_rdkit() {
        let atom = parse_mol2_atom_line_like_rdkit("1 LP1 0 0 0 LP").unwrap();

        assert_eq!(atom, None);
    }

    #[test]
    fn mol2_atom_type_any_query_matches_anything_like_rdkit() {
        let atom = parse_mol2_atom_line_like_rdkit("1 A1 0 0 0 ANY")
            .unwrap()
            .unwrap();

        assert_eq!(atom.atom_spec.element(), Element::DUMMY);
        assert_eq!(
            atom.atom_spec.query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::Any))
        );
    }

    #[test]
    fn mol2_atom_type_du_query_matches_anything_like_rdkit() {
        let atom = parse_mol2_atom_line_like_rdkit("1 DU1 0 0 0 Du")
            .unwrap()
            .unwrap();

        assert_eq!(atom.atom_spec.element(), Element::DUMMY);
        assert_eq!(
            atom.atom_spec.query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::Any))
        );
    }

    #[test]
    fn mol2_atom_type_hev_query_matches_non_hydrogen_like_rdkit() {
        let atom = parse_mol2_atom_line_like_rdkit("1 Q1 0 0 0 HEV")
            .unwrap()
            .unwrap();

        assert_eq!(atom.atom_spec.element(), Element::H);
        assert_eq!(
            atom.atom_spec.query(),
            Some(&QueryNode::not(QueryNode::predicate(
                AtomQueryPredicate::AtomicNumber(1)
            )))
        );
    }

    #[test]
    fn mol2_atom_type_het_query_matches_n_o_p_s_like_rdkit() {
        let atom = parse_mol2_atom_line_like_rdkit("1 HET1 0 0 0 HET")
            .unwrap()
            .unwrap();

        assert_eq!(atom.atom_spec.element(), Element::N);
        assert_eq!(
            atom.atom_spec.query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::AtomicNumberIn(
                vec![7, 8, 15, 16]
            )))
        );
    }

    #[test]
    fn mol2_atom_type_hal_query_matches_f_cl_br_i_like_rdkit() {
        let atom = parse_mol2_atom_line_like_rdkit("1 HAL1 0 0 0 HAL")
            .unwrap()
            .unwrap();

        assert_eq!(atom.atom_spec.element(), Element::F);
        assert_eq!(
            atom.atom_spec.query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::AtomicNumberIn(
                vec![9, 17, 35, 53]
            )))
        );
    }

    #[test]
    fn mol2_atom_type_ordinary_symbols_use_periodic_table_like_rdkit() {
        let atom = parse_mol2_atom_line_like_rdkit("1 CL1 0 0 0 Cl")
            .unwrap()
            .unwrap();

        assert_eq!(atom.atom_spec.element(), Element::CL);
        assert_eq!(atom.atom_spec.query(), None);
    }

    #[test]
    fn mol2_atom_type_rejects_unknown_symbol_like_rdkit() {
        let err = parse_mol2_atom_line_like_rdkit("1 XX1 0 0 0 Xx").unwrap_err();

        assert_eq!(err.to_string(), "MOL2 parse failed: Element 'Xx' not found");
    }

    #[test]
    fn mol2_atom_block_preserves_lp_index_correspondence_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n3 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 1 2 3 C.3\n2 LP1 4 5 6 LP\n3 H1 7 8 9 H\n@<TRIPOS>BOND\n";
        let offsets = scan_mol2_sections_like_rdkit(input).unwrap();
        let block = parse_mol2_atom_block_like_rdkit(input, offsets.atom_start, 3).unwrap();
        let molecule = block.builder.clone().build().unwrap();

        assert_eq!(molecule.num_atoms(), 2);
        assert_eq!(block.idx_corresp[0], Some(AtomId::new(0)));
        assert_eq!(block.idx_corresp[1], None);
        assert_eq!(block.idx_corresp[2], Some(AtomId::new(1)));
        assert_eq!(
            molecule.conformers_3d()[0].coordinates(),
            &[[1.0, 2.0, 3.0], [7.0, 8.0, 9.0]]
        );
    }

    #[test]
    fn mol2_atom_block_creates_3d_conformer_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n2 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 -1.25 2.5 0 C.3\n2 O1 3 4.5 -6 O.2\n";
        let offsets = scan_mol2_sections_like_rdkit(input).unwrap();
        let block = parse_mol2_atom_block_like_rdkit(input, offsets.atom_start, 2).unwrap();
        let molecule = block.builder.build().unwrap();

        assert_eq!(molecule.conformers_3d().len(), 1);
        assert_eq!(molecule.conformers_3d()[0].id(), 0);
        assert!(molecule.conformers_3d()[0].is_3d());
        assert_eq!(
            molecule.conformers_3d()[0].coordinates(),
            &[[-1.25, 2.5, 0.0], [3.0, 4.5, -6.0]]
        );
    }

    #[test]
    fn mol2_atom_block_detects_explicit_hydrogen_atoms_like_rdkit() {
        let without_h =
            "@<TRIPOS>MOLECULE\nmol\n1 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n";
        let without_h_offsets = scan_mol2_sections_like_rdkit(without_h).unwrap();
        let without_h_block =
            parse_mol2_atom_block_like_rdkit(without_h, without_h_offsets.atom_start, 1).unwrap();

        let with_h =
            "@<TRIPOS>MOLECULE\nmol\n1 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 H1 0 0 0 H\n";
        let with_h_offsets = scan_mol2_sections_like_rdkit(with_h).unwrap();
        let with_h_block =
            parse_mol2_atom_block_like_rdkit(with_h, with_h_offsets.atom_start, 1).unwrap();

        assert!(!without_h_block.has_h_atoms);
        assert!(with_h_block.has_h_atoms);
    }

    #[test]
    fn mol2_atom_block_rejects_premature_eof_like_rdkit() {
        let input =
            "@<TRIPOS>MOLECULE\nmol\n2 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n";
        let offsets = scan_mol2_sections_like_rdkit(input).unwrap();
        let err = parse_mol2_atom_block_like_rdkit(input, offsets.atom_start, 2).unwrap_err();

        assert_eq!(err.to_string(), "MOL2 parse failed: premature EOF");
    }

    #[test]
    fn mol2_atom_block_postcondition_counts_non_lp_atoms_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n4 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 LP1 0 0 0 LP\n2 C1 1 0 0 C.3\n3 LP2 2 0 0 LP\n4 O1 3 0 0 O.2\n";
        let offsets = scan_mol2_sections_like_rdkit(input).unwrap();
        let block = parse_mol2_atom_block_like_rdkit(input, offsets.atom_start, 4).unwrap();
        let molecule = block.builder.build().unwrap();

        assert_eq!(molecule.num_atoms(), 2);
        assert_eq!(
            block
                .idx_corresp
                .iter()
                .filter(|atom| atom.is_some())
                .count(),
            2
        );
        assert_eq!(
            molecule.conformers_3d()[0].coordinates(),
            &[[1.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
        );
    }

    #[test]
    fn mol2_bond_line_converts_one_based_indices_like_rdkit() {
        let idx_corresp = vec![Some(AtomId::new(10)), Some(AtomId::new(11))];
        let bond = parse_mol2_bond_line_like_rdkit("7 1 2 1", &idx_corresp)
            .unwrap()
            .unwrap();

        assert_eq!(bond.bond_spec.begin(), AtomId::new(10));
        assert_eq!(bond.bond_spec.end(), AtomId::new(11));
        assert_eq!(bond.zero_based_begin_idx, 0);
        assert_eq!(bond.zero_based_end_idx, 1);
    }

    #[test]
    fn mol2_bond_line_skips_lp_endpoints_like_rdkit() {
        let idx_corresp = vec![Some(AtomId::new(0)), None, Some(AtomId::new(1))];

        assert_eq!(
            parse_mol2_bond_line_like_rdkit("1 1 2 1", &idx_corresp).unwrap(),
            None
        );
        assert_eq!(
            parse_mol2_bond_line_like_rdkit("2 2 3 1", &idx_corresp).unwrap(),
            None
        );
    }

    #[test]
    fn mol2_bond_line_rejects_index_mismatch_like_rdkit() {
        let idx_corresp = vec![Some(AtomId::new(0)), Some(AtomId::new(1))];
        let err = parse_mol2_bond_line_like_rdkit("1 9 10 1", &idx_corresp).unwrap_err();

        assert_eq!(err.to_string(), "MOL2 parse failed: index mismatch");
    }

    #[test]
    fn mol2_bond_line_maps_supported_bond_types_like_rdkit() {
        let idx_corresp = vec![Some(AtomId::new(0)), Some(AtomId::new(1))];
        let cases = [
            ("1", BondOrder::Single),
            ("am", BondOrder::Single),
            ("2", BondOrder::Double),
            ("3", BondOrder::Triple),
            ("ar", BondOrder::Aromatic),
            ("du", BondOrder::Unspecified),
            ("un", BondOrder::Unspecified),
        ];

        for (sybyl_type, expected_order) in cases {
            let line = format!("1 1 2 {sybyl_type}");
            let bond = parse_mol2_bond_line_like_rdkit(&line, &idx_corresp)
                .unwrap()
                .unwrap();

            assert_eq!(bond.bond_spec.order(), expected_order);
            assert_eq!(bond.sybyl_bond_type, sybyl_type);
        }
    }

    #[test]
    fn mol2_bond_line_skips_nc_and_unsupported_types_like_rdkit() {
        let idx_corresp = vec![Some(AtomId::new(0)), Some(AtomId::new(1))];

        assert_eq!(
            parse_mol2_bond_line_like_rdkit("1 1 2 nc", &idx_corresp).unwrap(),
            None
        );
        assert_eq!(
            parse_mol2_bond_line_like_rdkit("2 1 2 weird", &idx_corresp).unwrap(),
            None
        );
    }

    #[test]
    fn mol2_bond_block_sets_aromatic_bond_and_atom_flags_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n2 1\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.ar\n2 C2 1 0 0 C.ar\n@<TRIPOS>BOND\n1 1 2 ar\n";
        let offsets = scan_mol2_sections_like_rdkit(input).unwrap();
        let atom_block = parse_mol2_atom_block_like_rdkit(input, offsets.atom_start, 2).unwrap();
        let bond_block =
            parse_mol2_bond_block_like_rdkit(input, offsets.bond_start.unwrap(), 1, atom_block)
                .unwrap();
        let molecule = bond_block.builder.build().unwrap();

        assert_eq!(molecule.num_bonds(), 1);
        assert_eq!(molecule.bonds()[0].order(), BondOrder::Aromatic);
        assert!(molecule.bonds()[0].is_aromatic());
        assert!(molecule.atoms()[0].is_aromatic());
        assert!(molecule.atoms()[1].is_aromatic());
    }

    #[test]
    fn mol2_bond_block_skips_bad_bonds_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n3 2\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n2 LP1 1 0 0 LP\n3 O1 2 0 0 O.2\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 nc\n";
        let offsets = scan_mol2_sections_like_rdkit(input).unwrap();
        let atom_block = parse_mol2_atom_block_like_rdkit(input, offsets.atom_start, 3).unwrap();
        let bond_block =
            parse_mol2_bond_block_like_rdkit(input, offsets.bond_start.unwrap(), 2, atom_block)
                .unwrap();
        let molecule = bond_block.builder.build().unwrap();

        assert_eq!(bond_block.n_bad_bonds, 2);
        assert_eq!(molecule.num_bonds(), 0);
    }

    #[test]
    fn mol2_bond_block_rejects_premature_eof_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n2 2\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n2 O1 1 0 0 O.2\n@<TRIPOS>BOND\n1 1 2 1\n";
        let offsets = scan_mol2_sections_like_rdkit(input).unwrap();
        let atom_block = parse_mol2_atom_block_like_rdkit(input, offsets.atom_start, 2).unwrap();
        let err =
            parse_mol2_bond_block_like_rdkit(input, offsets.bond_start.unwrap(), 2, atom_block)
                .unwrap_err();

        assert_eq!(err.to_string(), "MOL2 parse failed: premature EOF");
    }

    #[test]
    fn mol2_bond_block_postcondition_counts_non_bad_bonds_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n3 3\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n2 C2 1 0 0 C.3\n3 O1 2 0 0 O.2\n@<TRIPOS>BOND\n1 1 2 1\n2 2 3 weird\n3 1 3 2\n";
        let offsets = scan_mol2_sections_like_rdkit(input).unwrap();
        let atom_block = parse_mol2_atom_block_like_rdkit(input, offsets.atom_start, 3).unwrap();
        let bond_block =
            parse_mol2_bond_block_like_rdkit(input, offsets.bond_start.unwrap(), 3, atom_block)
                .unwrap();
        let molecule = bond_block.builder.build().unwrap();

        assert_eq!(bond_block.n_bad_bonds, 1);
        assert_eq!(molecule.num_bonds(), 2);
        assert_eq!(molecule.bonds()[0].order(), BondOrder::Single);
        assert_eq!(molecule.bonds()[1].order(), BondOrder::Double);
    }

    #[test]
    fn mol2_unity_atom_attr_assigns_atomexpr_formal_charges_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n2 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.3\n2 O1 1 0 0 O.2\n@<TRIPOS>UNITY_ATOM_ATTR\n1 1\nAtomExpr 1\n2 1\nAtomExpr -1\n@<TRIPOS>SUBSTRUCTURE\n";
        let attr_block = read_unity_atom_attr_from_block(input, 2, 0).unwrap();
        let molecule = attr_block.builder.build().unwrap();

        assert_eq!(molecule.atoms()[0].formal_charge(), 1);
        assert_eq!(molecule.atoms()[1].formal_charge(), -1);
    }

    #[test]
    fn mol2_unity_atom_attr_rejects_malformed_header_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n1 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.3\n@<TRIPOS>UNITY_ATOM_ATTR\nnot-an-index 1\nAtomExpr 1\n";
        let err = read_unity_atom_attr_from_block(input, 1, 0).unwrap_err();

        assert_eq!(
            err.to_string(),
            "MOL2 parse failed: Cannot process mol2 UnityAtomAttr."
        );
    }

    #[test]
    fn mol2_unity_atom_attr_rejects_malformed_formal_charge_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n1 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.3\n@<TRIPOS>UNITY_ATOM_ATTR\n1 1\nAtomExpr not-a-charge\n@<TRIPOS>SUBSTRUCTURE\n";
        let err = read_unity_atom_attr_from_block(input, 1, 0).unwrap_err();

        assert_eq!(
            err.to_string(),
            "MOL2 parse failed: Cannot process mol2 formal charge."
        );
    }

    #[test]
    fn mol2_unity_atom_attr_stops_at_next_section_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n1 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.3\n@<TRIPOS>UNITY_ATOM_ATTR\n1 1\nAtomExpr 1\n@<TRIPOS>SUBSTRUCTURE\n1 <0> 1 RESIDUE 1 A **** 0 ROOT\n";
        let attr_block = read_unity_atom_attr_from_block(input, 1, 0).unwrap();
        let molecule = attr_block.builder.build().unwrap();

        assert_eq!(molecule.atoms()[0].formal_charge(), 1);
        assert!(input[attr_block.next_offset..].starts_with("1 <0> 1 RESIDUE"));
    }

    #[test]
    fn mol2_unity_atom_attr_stops_at_blank_line_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n2 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.3\n2 O1 1 0 0 O.2\n@<TRIPOS>UNITY_ATOM_ATTR\n1 1\nAtomExpr 1\n\n2 1\nAtomExpr -1\n";
        let attr_block = read_unity_atom_attr_from_block(input, 2, 0).unwrap();
        let molecule = attr_block.builder.build().unwrap();

        assert_eq!(molecule.atoms()[0].formal_charge(), 1);
        assert_eq!(molecule.atoms()[1].formal_charge(), 0);
        assert!(input[attr_block.next_offset..].starts_with("2 1\n"));
    }

    #[test]
    fn mol2_unity_atom_attr_stops_at_comment_line_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n2 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.3\n2 O1 1 0 0 O.2\n@<TRIPOS>UNITY_ATOM_ATTR\n1 1\nAtomExpr 1\n# done with attributes\n2 1\nAtomExpr -1\n";
        let attr_block = read_unity_atom_attr_from_block(input, 2, 0).unwrap();
        let molecule = attr_block.builder.build().unwrap();

        assert_eq!(molecule.atoms()[0].formal_charge(), 1);
        assert_eq!(molecule.atoms()[1].formal_charge(), 0);
        assert!(input[attr_block.next_offset..].starts_with("2 1\n"));
    }

    #[test]
    fn mol2_unity_atom_attr_rejects_premature_eof_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n1 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.3\n@<TRIPOS>UNITY_ATOM_ATTR\n";
        let err = read_unity_atom_attr_from_block(input, 1, 0).unwrap_err();

        assert_eq!(err.to_string(), "MOL2 parse failed: premature EOF");
    }

    #[test]
    fn mol2_guess_charges_assigns_four_valent_non_carbon_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n5 4\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.3\n2 C1 1 0 0 C.3\n3 C2 -1 0 0 C.3\n4 C3 0 1 0 C.3\n5 C4 0 -1 0 C.3\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n3 1 4 1\n4 1 5 1\n";
        let mut bond_block = parse_bond_block_from_block(input, 5, 4).unwrap();

        guess_formal_charges_like_rdkit(&mut bond_block.builder).unwrap();

        assert_eq!(bond_block.builder.atoms()[0].formal_charge(), 1);
        assert_eq!(bond_block.builder.atoms()[1].formal_charge(), 0);
        assert_eq!(bond_block.builder.atoms()[2].formal_charge(), 0);
        assert_eq!(bond_block.builder.atoms()[3].formal_charge(), 0);
        assert_eq!(bond_block.builder.atoms()[4].formal_charge(), 0);
    }

    #[test]
    fn mol2_guess_charges_skips_query_atoms_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n1 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 Q1 0 0 0 Du\n@<TRIPOS>BOND\n";
        let mut atom_block = parse_mol2_atom_block_like_rdkit(
            input,
            scan_mol2_sections_like_rdkit(input).unwrap().atom_start,
            1,
        )
        .unwrap();

        guess_formal_charges_like_rdkit(&mut atom_block.builder).unwrap();

        assert!(atom_block.builder.atoms()[0].query().is_some());
        assert_eq!(atom_block.builder.atoms()[0].formal_charge(), 0);
    }

    #[test]
    fn mol2_guess_charges_skips_aromatic_five_member_non_ar_type_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n5 5\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 O1 0 1 0 O.2\n2 C1 1 0 0 C.2\n3 C2 0 -1 0 C.2\n4 C3 -1 0 0 C.2\n5 C4 0 0 1 C.2\n@<TRIPOS>BOND\n1 1 2 ar\n2 2 3 ar\n3 3 4 ar\n4 4 5 ar\n5 5 1 ar\n";
        let mut bond_block = parse_bond_block_from_block(input, 5, 5).unwrap();

        guess_formal_charges_like_rdkit(&mut bond_block.builder).unwrap();

        assert!(bond_block.builder.atoms()[0].is_aromatic());
        assert_eq!(
            bond_block.builder.atoms()[0].prop("_TriposAtomType"),
            Some("O.2")
        );
        assert_eq!(bond_block.builder.atoms()[0].formal_charge(), 0);
    }

    #[test]
    fn mol2_guess_charges_skips_n_ar_with_three_aromatic_bonds_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n4 3\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.ar\n2 C1 1 0 0 C.ar\n3 C2 -1 0 0 C.ar\n4 C3 0 1 0 C.ar\n@<TRIPOS>BOND\n1 1 2 ar\n2 1 3 ar\n3 1 4 ar\n";
        let mut bond_block = parse_bond_block_from_block(input, 4, 3).unwrap();

        guess_formal_charges_like_rdkit(&mut bond_block.builder).unwrap();

        assert_eq!(bond_block.builder.atoms()[0].formal_charge(), 0);
    }

    #[test]
    fn mol2_guess_charges_uses_multi_valence_positive_preference_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n4 3\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 S1 0 0 0 S.3\n2 C1 1 0 0 C.3\n3 O1 -1 0 0 O.3\n4 O2 0 1 0 O.3\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n3 1 4 1\n";
        let mut bond_block = parse_bond_block_from_block(input, 4, 3).unwrap();

        guess_formal_charges_like_rdkit(&mut bond_block.builder).unwrap();

        assert_eq!(bond_block.builder.atoms()[0].formal_charge(), 1);
        assert_eq!(bond_block.builder.atoms()[2].formal_charge(), -1);
        assert_eq!(bond_block.builder.atoms()[3].formal_charge(), -1);
    }

    #[test]
    fn mol2_guess_charges_clamps_aromatic_charge_magnitude_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n4 3\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 O1 0 0 0 O.ar\n2 C1 1 0 0 C.ar\n3 C2 -1 0 0 C.ar\n4 C3 0 1 0 C.ar\n@<TRIPOS>BOND\n1 1 2 ar\n2 1 3 ar\n3 1 4 ar\n";
        let mut bond_block = parse_bond_block_from_block(input, 4, 3).unwrap();

        guess_formal_charges_like_rdkit(&mut bond_block.builder).unwrap();

        assert!(bond_block.builder.atoms()[0].is_aromatic());
        assert_eq!(bond_block.builder.atoms()[0].formal_charge(), 1);
    }

    #[test]
    fn mol2_guess_charges_hands_nitro_assignment_to_cleanup_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n4 3\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.3\n2 O1 1 0 0 O.2\n3 O2 -1 0 0 O.2\n4 C1 0 1 0 C.3\n@<TRIPOS>BOND\n1 1 2 2\n2 1 3 2\n3 1 4 1\n";
        let mut bond_block = parse_bond_block_from_block(input, 4, 3).unwrap();

        guess_formal_charges_like_rdkit(&mut bond_block.builder).unwrap();

        assert_eq!(bond_block.builder.atoms()[0].formal_charge(), 1);
        assert_eq!(bond_block.builder.atoms()[1].formal_charge(), 0);
        assert_eq!(bond_block.builder.atoms()[2].formal_charge(), -1);
        assert_eq!(bond_block.builder.bonds()[0].order(), BondOrder::Double);
        assert_eq!(bond_block.builder.bonds()[1].order(), BondOrder::Single);
        assert_eq!(bond_block.builder.bonds()[2].order(), BondOrder::Single);
    }

    #[test]
    fn mol2_data_stream_unsanitized_builds_atom_block_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n2 1\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n2 N1 1 0 0 N.3\n@<TRIPOS>BOND\n1 1 2 1\n";

        let record = read_unsanitized_mol2_from_block(input).unwrap().unwrap();

        assert_eq!(record.molecule.num_atoms(), 2);
        assert_eq!(record.molecule.properties().name(), Some("mol"));
        assert_eq!(
            record.molecule.prop("_TriposChargeType"),
            Some("NO_CHARGES")
        );
        assert_eq!(
            record.molecule.atoms()[0].prop("_TriposAtomName"),
            Some("C1")
        );
        assert_eq!(
            record.molecule.atoms()[1].prop("_TriposAtomType"),
            Some("N.3")
        );
        assert_eq!(
            record.molecule.conformers_3d()[0].coordinates()[1],
            [1.0, 0.0, 0.0]
        );
    }

    #[test]
    fn mol2_data_stream_unsanitized_builds_bond_block_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n2 1\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.2\n2 C2 1 0 0 C.2\n@<TRIPOS>BOND\n1 1 2 2\n";

        let record = read_unsanitized_mol2_from_block(input).unwrap().unwrap();

        assert_eq!(record.molecule.num_bonds(), 1);
        assert_eq!(record.molecule.bonds()[0].begin(), AtomId::new(0));
        assert_eq!(record.molecule.bonds()[0].end(), AtomId::new(1));
        assert_eq!(record.molecule.bonds()[0].order(), BondOrder::Double);
    }

    #[test]
    fn mol2_data_stream_unsanitized_accepts_no_bonds_like_rdkit() {
        let input =
            "@<TRIPOS>MOLECULE\nmol\n1 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 O1 0 0 0 O.3\n";

        let record = read_unsanitized_mol2_from_block(input).unwrap().unwrap();

        assert_eq!(record.molecule.num_atoms(), 1);
        assert_eq!(record.molecule.num_bonds(), 0);
        assert_eq!(record.molecule.atoms()[0].atomic_number(), 8);
        assert_eq!(record.molecule.atoms()[0].formal_charge(), -2);
    }

    #[test]
    fn mol2_data_stream_unsanitized_runs_cleanup_when_enabled_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n3 2\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.2\n2 O1 1 0 0 O.co2\n3 O2 -1 0 0 O.co2\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 2\n";

        let record = read_unsanitized_mol2_from_block(input).unwrap().unwrap();

        assert_eq!(record.molecule.atoms()[1].formal_charge(), -1);
        assert_eq!(record.molecule.atoms()[2].formal_charge(), 0);
        assert_eq!(record.molecule.bonds()[0].order(), BondOrder::Single);
        assert_eq!(record.molecule.bonds()[1].order(), BondOrder::Double);
    }

    #[test]
    fn mol2_data_stream_unsanitized_skips_cleanup_when_disabled_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n3 2\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.2\n2 O1 1 0 0 O.co2\n3 O2 -1 0 0 O.co2\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 2\n";
        let params = Mol2ReadParams {
            sanitize: false,
            remove_hs: false,
            variant: Mol2Type::Corina,
            cleanup_substructures: false,
        };

        let record = mol_from_mol2_data_stream_like_rdkit(input, params)
            .unwrap()
            .unwrap();

        assert_eq!(record.molecule.atoms()[1].formal_charge(), -1);
        assert_eq!(record.molecule.atoms()[2].formal_charge(), 0);
        assert_eq!(record.molecule.bonds()[0].order(), BondOrder::Single);
        assert_eq!(record.molecule.bonds()[1].order(), BondOrder::Double);
    }

    #[test]
    fn mol2_data_stream_unsanitized_unity_atom_attr_overrides_charge_guess_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n1 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.3\n@<TRIPOS>UNITY_ATOM_ATTR\n1 1\nAtomExpr -1\n";

        let record = read_unsanitized_mol2_from_block(input).unwrap().unwrap();

        assert_eq!(record.molecule.atoms()[0].formal_charge(), -1);
    }

    #[test]
    fn mol2_data_stream_unsanitized_returns_null_for_failed_cleanup_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n3 2\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 O1 0 0 0 O.co2\n2 C1 1 0 0 C.3\n3 C2 -1 0 0 C.3\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n";

        let record = read_unsanitized_mol2_from_block(input).unwrap();

        assert!(record.is_none());
    }

    #[test]
    fn mol2_chirality_3d_assigns_before_hydrogen_removal_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n5 4\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n2 F1 1 0 0 F\n3 Cl1 0 1 0 Cl\n4 Br1 0 0 1 Br\n5 H1 -1 -1 -1 H\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n3 1 4 1\n4 1 5 1\n";
        let mut record = read_unsanitized_mol2_from_block(input).unwrap().unwrap();

        assert_eq!(record.molecule.num_atoms(), 5);
        assert_eq!(
            record.molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::Unspecified
        );
        assert_eq!(record.molecule.atoms()[4].atomic_number(), 1);

        assign_chiral_types_from_3d_for_mol2_like_rdkit(&mut record.molecule);

        assert_eq!(record.molecule.num_atoms(), 5);
        assert_eq!(record.molecule.atoms()[4].atomic_number(), 1);
        assert_eq!(
            record.molecule.atoms()[0].chiral_tag(),
            crate::ChiralTag::TetrahedralCcw
        );
        assert_eq!(
            record.molecule.atoms()[0].prop("_NonExplicit3DChirality"),
            Some("1")
        );
    }

    #[test]
    fn mol2_sanitize_remove_hs_runs_cleanup_before_stereo_detection_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n4 3\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.3\n2 O1 1 0 0 O.2\n3 O2 -1 0 0 O.2\n4 H1 0 1 0 H\n@<TRIPOS>BOND\n1 1 2 2\n2 1 3 2\n3 1 4 1\n";
        let record = read_unsanitized_mol2_from_block(input).unwrap().unwrap();

        let sanitized =
            sanitize_mol2_molecule_like_rdkit(record.molecule, Mol2ReadParams::default()).unwrap();

        assert_eq!(sanitized.num_atoms(), 3);
        assert_eq!(sanitized.atoms()[0].formal_charge(), 1);
        assert_eq!(sanitized.atoms()[1].formal_charge(), 0);
        assert_eq!(sanitized.atoms()[2].formal_charge(), -1);
        assert_eq!(sanitized.bonds()[0].order(), BondOrder::Double);
        assert_eq!(sanitized.bonds()[1].order(), BondOrder::Single);
    }

    #[test]
    fn mol2_sanitize_remove_hs_detects_bond_stereo_before_hydrogen_removal_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n6 5\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.2\n2 C2 1 0 0 C.2\n3 F1 0 1 1 F\n4 Cl1 1 -1 -1 Cl\n5 H1 0 -1 -1 H\n6 H2 1 1 1 H\n@<TRIPOS>BOND\n1 1 2 2\n2 1 3 1\n3 2 4 1\n4 1 5 1\n5 2 6 1\n";
        let record = read_unsanitized_mol2_from_block(input).unwrap().unwrap();

        let sanitized =
            sanitize_mol2_molecule_like_rdkit(record.molecule, Mol2ReadParams::default()).unwrap();

        assert_eq!(sanitized.num_atoms(), 4);
        assert!(matches!(
            sanitized.bonds()[1].direction(),
            crate::BondDirection::EndUpRight | crate::BondDirection::EndDownRight
        ));
        assert!(matches!(
            sanitized.bonds()[2].direction(),
            crate::BondDirection::EndUpRight | crate::BondDirection::EndDownRight
        ));
    }

    #[test]
    fn mol2_sanitize_remove_hs_uses_remove_hs_without_resanitize_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n2 1\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n2 H1 1 0 0 H\n@<TRIPOS>BOND\n1 1 2 1\n";
        let record = read_unsanitized_mol2_from_block(input).unwrap().unwrap();
        let cleanup_then_remove_hs = record
            .molecule
            .clone()
            .sanitize_with_ops(crate::SanitizeOps::CLEANUP)
            .unwrap()
            .without_hydrogens_with_sanitize(false)
            .unwrap();

        assert!(cleanup_then_remove_hs.derived_cache().valence.is_none());
        assert!(cleanup_then_remove_hs.derived_cache().rings.is_none());

        let sanitized =
            sanitize_mol2_molecule_like_rdkit(record.molecule, Mol2ReadParams::default()).unwrap();

        assert_eq!(sanitized.num_atoms(), 1);
        assert!(sanitized.derived_cache().valence.is_some());
        assert!(sanitized.derived_cache().rings.is_some());
    }

    #[test]
    fn mol2_sanitize_remove_hs_skips_organometallic_cleanup_like_rdkit() {
        let flags = mol2_sanitize_flags_like_rdkit();

        assert!(!flags.contains(crate::SanitizeOps::CLEANUP_ORGANOMETALLICS));
        assert!(flags.contains(crate::SanitizeOps::CLEANUP));
        assert!(flags.contains(crate::SanitizeOps::PROPERTIES));
        assert!(flags.contains(crate::SanitizeOps::SET_HYBRIDIZATION));
    }

    #[test]
    fn mol2_sanitize_remove_hs_runs_final_sanitize_property_cache_and_stereo_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n5 4\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n2 F1 1 0 0 F\n3 Cl1 0 1 0 Cl\n4 Br1 0 0 1 Br\n5 H1 -1 -1 -1 H\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n3 1 4 1\n4 1 5 1\n";
        let mut record = read_unsanitized_mol2_from_block(input).unwrap().unwrap();
        assign_chiral_types_from_3d_for_mol2_like_rdkit(&mut record.molecule);

        let sanitized =
            sanitize_mol2_molecule_like_rdkit(record.molecule, Mol2ReadParams::default()).unwrap();

        assert_eq!(sanitized.num_atoms(), 4);
        assert!(sanitized.derived_cache().valence.is_some());
        assert!(sanitized.derived_cache().rings.is_some());
        assert_eq!(sanitized.prop("_StereochemDone"), Some("1"));
        assert_eq!(
            sanitized.atoms()[0].chiral_tag(),
            crate::ChiralTag::TetrahedralCcw
        );
    }

    #[test]
    fn mol2_sanitize_keep_hs_preserves_explicit_hydrogens_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n2 1\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n2 H1 1 0 0 H\n@<TRIPOS>BOND\n1 1 2 1\n";
        let record = read_unsanitized_mol2_from_block(input).unwrap().unwrap();
        let params = Mol2ReadParams {
            sanitize: true,
            remove_hs: false,
            variant: Mol2Type::Corina,
            cleanup_substructures: true,
        };

        let sanitized = sanitize_mol2_molecule_like_rdkit(record.molecule, params).unwrap();

        assert_eq!(sanitized.num_atoms(), 2);
        assert_eq!(sanitized.atoms()[1].atomic_number(), 1);
    }

    #[test]
    fn mol2_sanitize_keep_hs_skips_organometallic_cleanup_like_rdkit() {
        let flags = mol2_sanitize_flags_like_rdkit();

        assert!(!flags.contains(crate::SanitizeOps::CLEANUP_ORGANOMETALLICS));
        assert!(flags.contains(crate::SanitizeOps::CLEANUP));
        assert!(flags.contains(crate::SanitizeOps::CLEANUP_ATROPISOMERS));
    }

    #[test]
    fn mol2_sanitize_keep_hs_detects_bond_stereo_after_sanitize_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n6 5\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.2\n2 C2 1 0 0 C.2\n3 F1 0 1 1 F\n4 Cl1 1 -1 -1 Cl\n5 H1 0 -1 -1 H\n6 H2 1 1 1 H\n@<TRIPOS>BOND\n1 1 2 2\n2 1 3 1\n3 2 4 1\n4 1 5 1\n5 2 6 1\n";
        let record = read_unsanitized_mol2_from_block(input).unwrap().unwrap();
        let params = Mol2ReadParams {
            sanitize: true,
            remove_hs: false,
            variant: Mol2Type::Corina,
            cleanup_substructures: true,
        };

        let sanitized = sanitize_mol2_molecule_like_rdkit(record.molecule, params).unwrap();

        assert_eq!(sanitized.num_atoms(), 6);
        assert!(matches!(
            sanitized.bonds()[1].direction(),
            crate::BondDirection::EndUpRight | crate::BondDirection::EndDownRight
        ));
        assert!(matches!(
            sanitized.bonds()[2].direction(),
            crate::BondDirection::EndUpRight | crate::BondDirection::EndDownRight
        ));
    }

    #[test]
    fn mol2_sanitize_keep_hs_runs_property_cache_and_stereo_assignment_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n5 4\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n2 F1 1 0 0 F\n3 Cl1 0 1 0 Cl\n4 Br1 0 0 1 Br\n5 H1 -1 -1 -1 H\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n3 1 4 1\n4 1 5 1\n";
        let mut record = read_unsanitized_mol2_from_block(input).unwrap().unwrap();
        assign_chiral_types_from_3d_for_mol2_like_rdkit(&mut record.molecule);
        let params = Mol2ReadParams {
            sanitize: true,
            remove_hs: false,
            variant: Mol2Type::Corina,
            cleanup_substructures: true,
        };

        let sanitized = sanitize_mol2_molecule_like_rdkit(record.molecule, params).unwrap();

        assert_eq!(sanitized.num_atoms(), 5);
        assert!(sanitized.derived_cache().valence.is_some());
        assert!(sanitized.derived_cache().rings.is_some());
        assert_eq!(sanitized.prop("_StereochemDone"), Some("1"));
        assert_eq!(
            sanitized.atoms()[0].chiral_tag(),
            crate::ChiralTag::TetrahedralCcw
        );
    }

    #[test]
    fn mol2_nitro_cleanup_fixes_two_double_bond_oxygen_neighbors_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n4 3\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.3\n2 O1 1 0 0 O.2\n3 O2 -1 0 0 O.2\n4 C1 0 1 0 C.3\n@<TRIPOS>BOND\n1 1 2 2\n2 1 3 2\n3 1 4 1\n";
        let mut bond_block = parse_bond_block_from_block(input, 4, 3).unwrap();

        fix_nitro_substructure_and_charge_like_rdkit(&mut bond_block.builder, AtomId::new(0))
            .unwrap();

        assert_eq!(bond_block.builder.atoms()[0].formal_charge(), 1);
        assert_eq!(bond_block.builder.atoms()[1].formal_charge(), 0);
        assert_eq!(bond_block.builder.atoms()[2].formal_charge(), -1);
        assert_eq!(bond_block.builder.bonds()[0].order(), BondOrder::Double);
        assert_eq!(bond_block.builder.bonds()[1].order(), BondOrder::Single);
        assert_eq!(bond_block.builder.bonds()[2].order(), BondOrder::Single);
    }

    #[test]
    fn mol2_nitro_cleanup_leaves_non_matching_nitrogen_environments_like_rdkit() {
        let one_double = "@<TRIPOS>MOLECULE\nmol\n3 2\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.3\n2 O1 1 0 0 O.2\n3 O2 -1 0 0 O.2\n@<TRIPOS>BOND\n1 1 2 2\n2 1 3 1\n";
        let mut one_double_block = parse_bond_block_from_block(one_double, 3, 2).unwrap();

        fix_nitro_substructure_and_charge_like_rdkit(&mut one_double_block.builder, AtomId::new(0))
            .unwrap();

        assert_eq!(one_double_block.builder.atoms()[0].formal_charge(), 0);
        assert_eq!(one_double_block.builder.atoms()[1].formal_charge(), 0);
        assert_eq!(one_double_block.builder.atoms()[2].formal_charge(), 0);
        assert_eq!(
            one_double_block.builder.bonds()[0].order(),
            BondOrder::Double
        );
        assert_eq!(
            one_double_block.builder.bonds()[1].order(),
            BondOrder::Single
        );

        let three_double = "@<TRIPOS>MOLECULE\nmol\n4 3\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.3\n2 O1 1 0 0 O.2\n3 O2 -1 0 0 O.2\n4 O3 0 1 0 O.2\n@<TRIPOS>BOND\n1 1 2 2\n2 1 3 2\n3 1 4 2\n";
        let mut three_double_block = parse_bond_block_from_block(three_double, 4, 3).unwrap();

        fix_nitro_substructure_and_charge_like_rdkit(
            &mut three_double_block.builder,
            AtomId::new(0),
        )
        .unwrap();

        assert_eq!(three_double_block.builder.atoms()[0].formal_charge(), 0);
        assert_eq!(three_double_block.builder.atoms()[1].formal_charge(), 0);
        assert_eq!(three_double_block.builder.atoms()[2].formal_charge(), 0);
        assert_eq!(three_double_block.builder.atoms()[3].formal_charge(), 0);
        assert_eq!(
            three_double_block.builder.bonds()[0].order(),
            BondOrder::Double
        );
        assert_eq!(
            three_double_block.builder.bonds()[1].order(),
            BondOrder::Double
        );
        assert_eq!(
            three_double_block.builder.bonds()[2].order(),
            BondOrder::Double
        );
    }

    #[test]
    fn mol2_n_oxide_helper_counts_hydrogen_neighbors_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n4 3\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.3\n2 H1 1 0 0 H\n3 H2 -1 0 0 H\n4 C1 0 1 0 C.3\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n3 1 4 1\n";
        let bond_block = parse_bond_block_from_block(input, 4, 3).unwrap();
        let mut to_mod_idx = None;

        let no_h_neighbors = check_no_h_neighbors_n_oxide_like_rdkit(
            &bond_block.builder,
            AtomId::new(0),
            &mut to_mod_idx,
        )
        .unwrap();

        assert_eq!(no_h_neighbors, 2);
        assert_eq!(to_mod_idx, None);
    }

    #[test]
    fn mol2_n_oxide_helper_detects_terminal_oxygen_neighbor_like_rdkit() {
        let terminal_input = "@<TRIPOS>MOLECULE\nmol\n2 1\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.3\n2 O1 1 0 0 O.2\n@<TRIPOS>BOND\n1 1 2 1\n";
        let terminal_block = parse_bond_block_from_block(terminal_input, 2, 1).unwrap();
        let mut terminal_to_mod_idx = None;

        let terminal_h_neighbors = check_no_h_neighbors_n_oxide_like_rdkit(
            &terminal_block.builder,
            AtomId::new(0),
            &mut terminal_to_mod_idx,
        )
        .unwrap();

        assert_eq!(terminal_h_neighbors, 0);
        assert_eq!(terminal_to_mod_idx, Some(AtomId::new(0)));

        let nonterminal_input = "@<TRIPOS>MOLECULE\nmol\n3 2\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.3\n2 O1 1 0 0 O.2\n3 C1 2 0 0 C.3\n@<TRIPOS>BOND\n1 1 2 1\n2 2 3 1\n";
        let nonterminal_block = parse_bond_block_from_block(nonterminal_input, 3, 2).unwrap();
        let mut nonterminal_to_mod_idx = None;

        let nonterminal_h_neighbors = check_no_h_neighbors_n_oxide_like_rdkit(
            &nonterminal_block.builder,
            AtomId::new(0),
            &mut nonterminal_to_mod_idx,
        )
        .unwrap();

        assert_eq!(nonterminal_h_neighbors, 0);
        assert_eq!(nonterminal_to_mod_idx, None);
    }

    #[test]
    fn mol2_cleanup_o_co2_assigns_n4_formal_charge_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n1 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.4\n@<TRIPOS>BOND\n";
        let mut atom_block = parse_mol2_atom_block_like_rdkit(
            input,
            scan_mol2_sections_like_rdkit(input).unwrap().atom_start,
            1,
        )
        .unwrap();

        assert!(clean_up_mol2_substructures_like_rdkit(&mut atom_block.builder).unwrap());
        assert_eq!(atom_block.builder.atoms()[0].formal_charge(), 1);
    }

    #[test]
    fn mol2_cleanup_o_co2_fixes_carboxylate_first_single_second_double_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n3 2\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.2\n2 O1 1 0 0 O.co2\n3 O2 -1 0 0 O.co2\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n";
        let mut bond_block = parse_bond_block_from_block(input, 3, 2).unwrap();

        assert!(clean_up_mol2_substructures_like_rdkit(&mut bond_block.builder).unwrap());
        assert_eq!(bond_block.builder.atoms()[0].formal_charge(), 0);
        assert_eq!(bond_block.builder.atoms()[1].formal_charge(), -1);
        assert_eq!(bond_block.builder.atoms()[2].formal_charge(), 0);
        assert_eq!(bond_block.builder.bonds()[0].order(), BondOrder::Single);
        assert_eq!(bond_block.builder.bonds()[1].order(), BondOrder::Double);
        assert!(!bond_block.builder.atoms()[0].is_aromatic());
        assert!(!bond_block.builder.atoms()[1].is_aromatic());
        assert!(!bond_block.builder.atoms()[2].is_aromatic());
    }

    #[test]
    fn mol2_cleanup_o_co2_fixes_sulfonate_first_single_second_double_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n3 2\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 S1 0 0 0 S.o2\n2 O1 1 0 0 O.co2\n3 O2 -1 0 0 O.co2\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n";
        let mut bond_block = parse_bond_block_from_block(input, 3, 2).unwrap();

        assert!(clean_up_mol2_substructures_like_rdkit(&mut bond_block.builder).unwrap());
        assert_eq!(bond_block.builder.atoms()[0].formal_charge(), 0);
        assert_eq!(bond_block.builder.atoms()[1].formal_charge(), -1);
        assert_eq!(bond_block.builder.atoms()[2].formal_charge(), 0);
        assert_eq!(bond_block.builder.bonds()[0].order(), BondOrder::Single);
        assert_eq!(bond_block.builder.bonds()[1].order(), BondOrder::Double);
    }

    #[test]
    fn mol2_cleanup_o_co2_fixes_phosphate_first_double_rest_single_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n4 3\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 P1 0 0 0 P.3\n2 O1 1 0 0 O.co2\n3 O2 -1 0 0 O.co2\n4 O3 0 1 0 O.co2\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n3 1 4 1\n";
        let mut bond_block = parse_bond_block_from_block(input, 4, 3).unwrap();

        assert!(clean_up_mol2_substructures_like_rdkit(&mut bond_block.builder).unwrap());
        assert_eq!(bond_block.builder.atoms()[0].formal_charge(), 0);
        assert_eq!(bond_block.builder.atoms()[1].formal_charge(), 0);
        assert_eq!(bond_block.builder.atoms()[2].formal_charge(), -1);
        assert_eq!(bond_block.builder.atoms()[3].formal_charge(), -1);
        assert_eq!(bond_block.builder.bonds()[0].order(), BondOrder::Double);
        assert_eq!(bond_block.builder.bonds()[1].order(), BondOrder::Single);
        assert_eq!(bond_block.builder.bonds()[2].order(), BondOrder::Single);
        assert!(!bond_block.builder.atoms()[0].is_aromatic());
        assert!(!bond_block.builder.atoms()[1].is_aromatic());
        assert!(!bond_block.builder.atoms()[2].is_aromatic());
        assert!(!bond_block.builder.atoms()[3].is_aromatic());
    }

    #[test]
    fn mol2_cleanup_o_co2_returns_false_for_degree_greater_than_one_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n3 2\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 O1 0 0 0 O.co2\n2 C1 1 0 0 C.2\n3 C2 -1 0 0 C.2\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n";
        let mut bond_block = parse_bond_block_from_block(input, 3, 2).unwrap();

        assert!(!clean_up_mol2_substructures_like_rdkit(&mut bond_block.builder).unwrap());
    }

    #[test]
    fn mol2_cleanup_o_co2_returns_false_for_unsupported_neighbor_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n2 1\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 O1 0 0 0 O.co2\n2 N1 1 0 0 N.3\n@<TRIPOS>BOND\n1 1 2 1\n";
        let mut bond_block = parse_bond_block_from_block(input, 2, 1).unwrap();

        assert!(!clean_up_mol2_substructures_like_rdkit(&mut bond_block.builder).unwrap());
    }

    #[test]
    fn mol2_cleanup_c_cat_two_n_uses_n_oxide_precedence_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n5 4\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.cat\n2 N1 1 0 0 N.3\n3 N2 -1 0 0 N.3\n4 O1 2 0 0 O.2\n5 H1 -2 0 0 H\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n3 2 4 1\n4 3 5 1\n";
        let mut bond_block = parse_bond_block_from_block(input, 5, 4).unwrap();

        assert!(clean_up_mol2_substructures_like_rdkit(&mut bond_block.builder).unwrap());
        assert_eq!(bond_block.builder.atoms()[1].formal_charge(), 1);
        assert_eq!(bond_block.builder.atoms()[2].formal_charge(), 0);
        assert_eq!(bond_block.builder.bonds()[0].order(), BondOrder::Double);
        assert_eq!(bond_block.builder.bonds()[1].order(), BondOrder::Single);
        assert!(!bond_block.builder.atoms()[0].is_aromatic());
        assert!(!bond_block.builder.atoms()[1].is_aromatic());
        assert!(!bond_block.builder.atoms()[2].is_aromatic());
    }

    #[test]
    fn mol2_cleanup_c_cat_two_n_uses_hydrogen_count_precedence_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n6 5\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.cat\n2 N1 1 0 0 N.3\n3 N2 -1 0 0 N.3\n4 H1 2 0 0 H\n5 H2 1 1 0 H\n6 H3 -2 0 0 H\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n3 2 4 1\n4 2 5 1\n5 3 6 1\n";
        let mut bond_block = parse_bond_block_from_block(input, 6, 5).unwrap();

        assert!(clean_up_mol2_substructures_like_rdkit(&mut bond_block.builder).unwrap());
        assert_eq!(bond_block.builder.atoms()[1].formal_charge(), 1);
        assert_eq!(bond_block.builder.atoms()[2].formal_charge(), 0);
        assert_eq!(bond_block.builder.bonds()[0].order(), BondOrder::Double);
        assert_eq!(bond_block.builder.bonds()[1].order(), BondOrder::Single);
    }

    #[test]
    fn mol2_cleanup_c_cat_two_n_uses_ring_precedence_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n5 5\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.cat\n2 N1 1 0 0 N.3\n3 N2 -1 0 0 N.3\n4 C2 2 0 0 C.3\n5 C3 1 1 0 C.3\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n3 2 4 1\n4 4 5 1\n5 5 2 1\n";
        let mut bond_block = parse_bond_block_from_block(input, 5, 5).unwrap();

        assert!(clean_up_mol2_substructures_like_rdkit(&mut bond_block.builder).unwrap());
        assert_eq!(bond_block.builder.atoms()[1].formal_charge(), 1);
        assert_eq!(bond_block.builder.atoms()[2].formal_charge(), 0);
        assert_eq!(bond_block.builder.bonds()[0].order(), BondOrder::Double);
        assert_eq!(bond_block.builder.bonds()[1].order(), BondOrder::Single);
    }

    #[test]
    fn mol2_cleanup_c_cat_two_n_uses_second_n_fallback_precedence_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n3 2\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.cat\n2 N1 1 0 0 N.3\n3 N2 -1 0 0 N.3\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n";
        let mut bond_block = parse_bond_block_from_block(input, 3, 2).unwrap();

        assert!(clean_up_mol2_substructures_like_rdkit(&mut bond_block.builder).unwrap());
        assert_eq!(bond_block.builder.atoms()[1].formal_charge(), 0);
        assert_eq!(bond_block.builder.atoms()[2].formal_charge(), 1);
        assert_eq!(bond_block.builder.bonds()[0].order(), BondOrder::Single);
        assert_eq!(bond_block.builder.bonds()[1].order(), BondOrder::Double);
    }

    #[test]
    fn mol2_cleanup_c_cat_three_n_uses_lowest_heavy_atom_degree_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n6 5\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.cat\n2 N1 1 0 0 N.3\n3 N2 -1 0 0 N.3\n4 N3 0 1 0 N.3\n5 C2 2 0 0 C.3\n6 C3 0 2 0 C.3\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n3 1 4 1\n4 2 5 1\n5 4 6 1\n";
        let mut bond_block = parse_bond_block_from_block(input, 6, 5).unwrap();

        assert!(clean_up_mol2_substructures_like_rdkit(&mut bond_block.builder).unwrap());
        assert_eq!(bond_block.builder.atoms()[1].formal_charge(), 0);
        assert_eq!(bond_block.builder.atoms()[2].formal_charge(), 1);
        assert_eq!(bond_block.builder.atoms()[3].formal_charge(), 0);
        assert_eq!(bond_block.builder.bonds()[0].order(), BondOrder::Single);
        assert_eq!(bond_block.builder.bonds()[1].order(), BondOrder::Double);
        assert_eq!(bond_block.builder.bonds()[2].order(), BondOrder::Single);
    }

    #[test]
    fn mol2_cleanup_c_cat_three_n_returns_false_for_bad_nitrogen_counts_like_rdkit() {
        let one_n_input = "@<TRIPOS>MOLECULE\nmol\n3 2\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.cat\n2 N1 1 0 0 N.3\n3 C2 -1 0 0 C.3\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n";
        let mut one_n_block = parse_bond_block_from_block(one_n_input, 3, 2).unwrap();

        assert!(!clean_up_mol2_substructures_like_rdkit(&mut one_n_block.builder).unwrap());

        let four_n_input = "@<TRIPOS>MOLECULE\nmol\n5 4\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.cat\n2 N1 1 0 0 N.3\n3 N2 -1 0 0 N.3\n4 N3 0 1 0 N.3\n5 N4 0 -1 0 N.3\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n3 1 4 1\n4 1 5 1\n";
        let mut four_n_block = parse_bond_block_from_block(four_n_input, 5, 4).unwrap();

        assert!(!clean_up_mol2_substructures_like_rdkit(&mut four_n_block.builder).unwrap());
    }

    #[test]
    fn mol2_cleanup_c_cat_three_n_sets_already_fixed_neighbor_bond_single_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n7 6\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.cat\n2 N1 1 0 0 N.3\n3 N2 -1 0 0 N.3\n4 N3 0 1 0 N.3\n5 C2 2 0 0 C.cat\n6 N4 3 0 0 N.3\n7 N5 2 1 0 N.3\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n3 1 4 1\n4 2 5 1\n5 5 6 1\n6 5 7 1\n";
        let mut bond_block = parse_bond_block_from_block(input, 7, 6).unwrap();

        assert!(clean_up_mol2_substructures_like_rdkit(&mut bond_block.builder).unwrap());
        assert_eq!(bond_block.builder.atoms()[1].formal_charge(), 0);
        assert_eq!(bond_block.builder.atoms()[2].formal_charge(), 1);
        assert_eq!(bond_block.builder.atoms()[5].formal_charge(), 1);
        assert_eq!(bond_block.builder.bonds()[3].order(), BondOrder::Single);
        assert_eq!(bond_block.builder.bonds()[4].order(), BondOrder::Double);
        assert_eq!(bond_block.builder.bonds()[5].order(), BondOrder::Single);
    }

    #[test]
    fn mol2_cleanup_c_cat_three_n_keeps_first_lowest_heavy_atom_degree_tie_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n4 3\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.cat\n2 N1 1 0 0 N.3\n3 N2 -1 0 0 N.3\n4 N3 0 1 0 N.3\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n3 1 4 1\n";
        let mut bond_block = parse_bond_block_from_block(input, 4, 3).unwrap();

        assert!(clean_up_mol2_substructures_like_rdkit(&mut bond_block.builder).unwrap());
        assert_eq!(bond_block.builder.atoms()[1].formal_charge(), 1);
        assert_eq!(bond_block.builder.atoms()[2].formal_charge(), 0);
        assert_eq!(bond_block.builder.atoms()[3].formal_charge(), 0);
        assert_eq!(bond_block.builder.bonds()[0].order(), BondOrder::Double);
        assert_eq!(bond_block.builder.bonds()[1].order(), BondOrder::Single);
        assert_eq!(bond_block.builder.bonds()[2].order(), BondOrder::Single);
    }

    #[test]
    fn mol2_section_rejects_missing_molecule_block_like_rdkit() {
        let err = scan_mol2_sections_like_rdkit("@<TRIPOS>ATOM\n").unwrap_err();
        assert_eq!(
            err.to_string(),
            "MOL2 parse failed: No MOLECULE block found in Mol2 data"
        );
    }

    #[test]
    fn mol2_section_rejects_missing_atom_block_like_rdkit() {
        let err = scan_mol2_sections_like_rdkit("@<TRIPOS>MOLECULE\nname\n").unwrap_err();
        assert_eq!(
            err.to_string(),
            "MOL2 parse failed: No ATOM block found in Mol2 data"
        );
    }

    #[test]
    fn mol2_section_records_required_and_optional_block_offsets_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nmol\n1 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n@<TRIPOS>BOND\n@<TRIPOS>UNITY_ATOM_ATTR\n";
        let offsets = scan_mol2_sections_like_rdkit(input).unwrap();

        assert!(input[offsets.molecule_start..].starts_with("mol\n"));
        assert!(input[offsets.atom_start..].starts_with("1 C1 0 0 0 C.3\n"));
        assert_eq!(
            offsets.bond_start.map(|offset| &input[offset..]),
            Some("@<TRIPOS>UNITY_ATOM_ATTR\n")
        );
        assert_eq!(offsets.charge_start, Some(input.len()));
    }

    #[test]
    fn mol2_section_stops_at_repeated_molecule_block_like_rdkit() {
        let input = "@<TRIPOS>MOLECULE\nfirst\n1 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n@<TRIPOS>MOLECULE\nsecond\n";
        let offsets = scan_mol2_sections_like_rdkit(input).unwrap();

        assert!(input[offsets.molecule_start..].starts_with("first\n"));
        assert!(input[offsets.atom_start..].starts_with("1 C1 0 0 0 C.3\n"));
        assert_eq!(offsets.bond_start, None);
        assert_eq!(offsets.charge_start, None);
    }
}
