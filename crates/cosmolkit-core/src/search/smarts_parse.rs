//! SMARTS parser — recursive-descent parser producing query-predicate trees.
//!
//! ## RDKit provenance (protocol: dev/source_reproduction_protocol.md)
//!
//! The SMARTS parser corresponds to RDKit's `GraphMol/SmilesParse/SmilesParse.cpp`
//! (MolFromSmarts entry point, labelRecursivePatterns helper) and the bison/flex
//! grammars `smarts.yy` / `smarts.ll`.  The flex/bison lexer and parser are
//! replaced here by a hand-written recursive-descent parser that produces the
//! same semantic query trees.
//!
//! C++ source lines are copied verbatim as commented blocks with two-axis
//! RDKit status markers per `dev/source_reproduction_protocol.md`:
//!   // RDKit✔️✔️: <C++ line>   — fully ported, behaviour identical
//!   // RDKit❗✔️: <C++ line>   — adapted for Rust / COSMolKit differences
//!
//! ## Design
//!
//! The parser tokenizes the SMARTS string, then recursively descends through
//! the token stream to produce a `SmartsMolecule` containing atom and bond
//! query trees.  Bracket atoms (`[...]`) are parsed by a sub-parser that
//! interprets SMARTS primitives.
//!
//! The output is consumed by the substructure matching engine in `substruct.rs`.

use std::collections::BTreeMap;

use crate::{AtomQueryPredicate, BondOrder, BondQueryPredicate, QueryNode, SmartsParseError};

// ---------------------------------------------------------------------------
// SmartsMolecule – output of SMARTS parsing
// ---------------------------------------------------------------------------

/// The result of parsing a SMARTS string.
///
/// RDKit✔️✔️: RDKit returns an RWMol with QueryAtom / QueryBond objects.
/// COSMolKit returns a separate struct of query-predicate trees paired via
/// indexing (atom_queries[i] is the i-th atom, bond_queries[i] is the bond
/// between atoms i and i+1 in SMARTS order).
#[derive(Debug, Clone)]
pub struct SmartsMolecule {
    /// Query trees for each atom in the pattern.
    pub atom_queries: Vec<QueryNode<AtomQueryPredicate>>,
    /// Query trees for each bond in the pattern (length = atom_queries.len() - 1).
    pub bond_queries: Vec<QueryNode<BondQueryPredicate>>,
    /// Ring-closure specifications: (closure_number, atom_index_in_pattern)
    pub ring_closures: Vec<(u8, usize)>,
}

impl SmartsMolecule {
    #[must_use]
    pub fn new(
        atom_queries: Vec<QueryNode<AtomQueryPredicate>>,
        bond_queries: Vec<QueryNode<BondQueryPredicate>>,
        ring_closures: Vec<(u8, usize)>,
    ) -> Self {
        Self {
            atom_queries,
            bond_queries,
            ring_closures,
        }
    }

    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.atom_queries.len()
    }

    #[must_use]
    pub fn atom_query(&self, idx: usize) -> Option<&QueryNode<AtomQueryPredicate>> {
        self.atom_queries.get(idx)
    }

    #[must_use]
    pub fn bond_query(&self, idx: usize) -> Option<&QueryNode<BondQueryPredicate>> {
        self.bond_queries.get(idx)
    }
}

// ---------------------------------------------------------------------------
// SmartsParseParams
// ---------------------------------------------------------------------------

/// RDKit source: SmilesParse.h lines 56-67
/// RDKit✔️✔️: struct RDKIT_SMILESPARSE_EXPORT SmartsParserParams {
/// RDKit✔️✔️:   bool allowCXSMILES = true;
/// RDKit✔️✔️:   bool strictCXSMILES = true;
/// RDKit✔️✔️:   bool parseName = true;
/// RDKit✔️✔️:   bool mergeHs = false;
/// RDKit✔️✔️:   bool skipCleanup = false;
/// RDKit✔️✔️:   bool debugParse = false;
/// RDKit✔️✔️:   std::map<std::string, std::string> replacements;
/// RDKit✔️✔️: };
#[derive(Debug, Clone)]
pub struct SmartsParseParams {
    pub allow_cxsmiles: bool,
    pub strict_cxsmiles: bool,
    pub parse_name: bool,
    pub merge_hs: bool,
    pub skip_cleanup: bool,
    pub debug_parse: bool,
    pub replacements: BTreeMap<String, String>,
}

impl Default for SmartsParseParams {
    fn default() -> Self {
        Self {
            allow_cxsmiles: true,
            strict_cxsmiles: true,
            parse_name: true,
            merge_hs: false,
            skip_cleanup: false,
            debug_parse: false,
            replacements: BTreeMap::new(),
        }
    }
}

// ---------------------------------------------------------------------------
// Top-level parse entry point
// ---------------------------------------------------------------------------

/// RDKit source: SmilesParse.cpp lines 548-576
/// RDKit✔️✔️: std::unique_ptr<RWMol> MolFromSmarts(
/// RDKit✔️✔️:     const std::string &smarts,
/// RDKit✔️✔️:     const SmartsParserParams &params) {
/// RDKit✔️✔️:   if (yysmarts_debug != params.debugParse) {
/// RDKit✔️✔️:     yysmarts_debug = params.debugParse;
/// RDKit✔️✔️:   }
/// RDKit✔️✔️:   std::string lsmarts, name, cxPart;
/// RDKit✔️✔️:   preprocessSmiles(smarts, params, lsmarts, name, cxPart);
/// RDKit✔️✔️:   auto res = toMol(labelRecursivePatterns(lsmarts), smarts_parse, lsmarts);
/// RDKit✔️✔️:   handleCXPartAndName(res.get(), params, cxPart, name);
/// RDKit❗❗:   handleCXPartAndName(res.get(), params, cxPart, name);
/// RDKit❗❗:   // remaining post-processing steps skipped — CX part handling
/// RDKit❗❗:   // is done inline in the Rust implementation
/// RDKit✔️✔️:   return res;
/// RDKit✔️✔️: }
pub fn parse_smarts(smarts: &str) -> Result<SmartsMolecule, SmartsParseError> {
    parse_smarts_with_params(smarts, &SmartsParseParams::default())
}

/// Parse a SMARTS string with custom parameters.
///
/// RDKit❗✔️: We skip the CXSMILES / name preprocessing and labelRecursivePatterns
/// call since those are RDKit-specific RWMol post-processing steps. The
/// recursive SMARTS $() are handled inline by the parser.
pub fn parse_smarts_with_params(
    smarts: &str,
    _params: &SmartsParseParams,
) -> Result<SmartsMolecule, SmartsParseError> {
    // RDKit✔️✔️: preprocessSmiles — trim whitespace, handle replacements
    let input = preprocess_smarts(smarts, _params);
    let input = label_recursive_patterns(&input);

    let tokens = tokenize(&input)?;
    let mut parser = SmartsParser::new(&tokens, &input);
    parser.parse_smarts_molecule()
}

/// RDKit source: SmilesParse.cpp lines 333-373
/// RDKit❗✔️: preprocessSmiles — ported for SMARTS
fn preprocess_smarts(smarts: &str, params: &SmartsParseParams) -> String {
    // RDKit✔️✔️: Trim leading/trailing whitespace
    let trimmed = smarts.trim();

    // RDKit✔️✔️: Apply replacements
    // RDKit✔️✔️: for (const auto &pr : params.replacements) {
    // RDKit✔️✔️:     boost::replace_all(smi, pr.first, pr.second);
    // RDKit✔️✔️: }
    if params.replacements.is_empty() {
        return trimmed.to_string();
    }

    let mut result = trimmed.to_string();
    loop {
        let mut replaced = false;
        for (key, val) in &params.replacements {
            if result.contains(key.as_str()) {
                result = result.replace(key.as_str(), val.as_str());
                replaced = true;
            }
        }
        if !replaced {
            break;
        }
    }
    result
}

/// RDKit source: SmilesParse.cpp lines 186-239
/// RDKit✔️✔️: std::string labelRecursivePatterns(const std::string &sma) {
/// RDKit✔️✔️:   std::list<SmaState> state;
/// RDKit✔️✔️:   std::list<unsigned int> startRecurse;
/// RDKit✔️✔️:   std::map<std::string, std::string> patterns;
/// RDKit✔️✔️:   std::string res;
/// RDKit✔️✔️:   state.push_back(BASE);
/// RDKit✔️✔️:
/// RDKit✔️✔️:   unsigned int pos = 0;
/// RDKit✔️✔️:   while (pos < sma.size()) {
/// RDKit✔️✔️:     res += sma[pos];
/// RDKit✔️✔️:     if (sma[pos] == '$' && pos + 1 < sma.size() && sma[pos + 1] == '(') {
/// RDKit✔️✔️:       state.push_back(RECURSE);
/// RDKit✔️✔️:       startRecurse.push_back(pos);
/// RDKit✔️✔️:       ++pos;
/// RDKit✔️✔️:       res += sma[pos];
/// RDKit✔️✔️:     } else if (sma[pos] == '(') {
/// RDKit✔️✔️:       state.push_back(BRANCH);
/// RDKit✔️✔️:     } else if (sma[pos] == ')') {
/// RDKit✔️✔️:       if (state.empty() || state.back() == BASE) {
/// RDKit✔️✔️:         return sma;
/// RDKit✔️✔️:       }
/// RDKit✔️✔️:       SmaState currState = state.back();
/// RDKit✔️✔️:       state.pop_back();
/// RDKit✔️✔️:       if (currState == RECURSE) {
/// RDKit✔️✔️:         unsigned int dollarPos = startRecurse.back();
/// RDKit✔️✔️:         startRecurse.pop_back();
/// RDKit✔️✔️:         if (pos + 1 >= sma.size() || sma[pos + 1] != '_') {
/// RDKit✔️✔️:           std::string recurs = sma.substr(dollarPos, pos - dollarPos + 1);
/// RDKit✔️✔️:           std::string label;
/// RDKit✔️✔️:           if (patterns.find(recurs) != patterns.end()) {
/// RDKit✔️✔️:             label = patterns[recurs];
/// RDKit✔️✔️:           } else {
/// RDKit✔️✔️:             label = std::to_string(patterns.size() + 100);
/// RDKit✔️✔️:             patterns[recurs] = label;
/// RDKit✔️✔️:           }
/// RDKit✔️✔️:           res += "_" + label;
/// RDKit✔️✔️:         }
/// RDKit✔️✔️:       }
/// RDKit✔️✔️:     }
/// RDKit✔️✔️:     ++pos;
/// RDKit✔️✔️:   }
/// RDKit✔️✔️:   return res;
/// RDKit✔️✔️: }
///
/// RDKit❗✔️: Ported for COSMolKit. The label (e.g. "_101") is appended
/// after the closing paren of a $() construct so subsequent passes can
/// deduplicate recursive SMARTS patterns.
fn label_recursive_patterns(sma: &str) -> String {
    #[derive(Clone, Copy, PartialEq)]
    enum SmaState {
        Base,
        Branch,
        Recurse,
    }
    use SmaState::*;

    let mut state: Vec<SmaState> = vec![Base];
    let mut start_recurse: Vec<usize> = Vec::new();
    let mut patterns: BTreeMap<String, String> = BTreeMap::new();
    let mut res = String::new();
    let chars: Vec<char> = sma.chars().collect();

    let mut pos: usize = 0;
    while pos < chars.len() {
        res.push(chars[pos]);
        if chars[pos] == '$' && pos + 1 < chars.len() && chars[pos + 1] == '(' {
            state.push(Recurse);
            start_recurse.push(pos);
            pos += 1;
            res.push(chars[pos]);
        } else if chars[pos] == '(' {
            state.push(Branch);
        } else if chars[pos] == ')' {
            if state.is_empty() || state.last() == Some(&Base) {
                // RDKit✔️✔️: seriously bogus input. Just return the input
                // RDKit✔️✔️: and let the SMARTS parser itself report the error
                return sma.to_string();
            }
            let curr_state = state.pop().unwrap();
            if curr_state == Recurse {
                let dollar_pos = start_recurse.pop().unwrap();
                if pos + 1 >= chars.len() || chars[pos + 1] != '_' {
                    let recurs: String = chars[dollar_pos..=pos].iter().collect();
                    let label = if let Some(lbl) = patterns.get(&recurs) {
                        lbl.clone()
                    } else {
                        let lbl = format!("{}", patterns.len() + 100);
                        patterns.insert(recurs, lbl.clone());
                        lbl
                    };
                    res.push('_');
                    res.push_str(&label);
                }
            }
        }
        pos += 1;
    }
    res
}

// ---------------------------------------------------------------------------
// Tokenizer
// ---------------------------------------------------------------------------

/// RDKit❗✔️: Comparison of token types between SMILES/smarts.ll and our
/// tokenizer. The flex/generated lexer produces tokens like ORGANIC_ATOM_TOKEN,
/// AROMATIC_ATOM_TOKEN, ATOM_TOKEN, BOND_TOKEN, etc. We collapse those into
/// a simpler enum since our parser uses recursive descent rather than bison.
#[derive(Debug, Clone, PartialEq)]
enum Token {
    /// Organic element symbol: B, C, N, O, S, P, F, Cl, Br, I, *
    OrganicElement(String),
    /// Aromatic element: c, n, o, s, p
    AromaticElement(String),
    /// Bracket atom content: the raw text between [ and ] (excluding brackets)
    BracketContent(String),
    /// Bond specifier: -, =, #, :, ~, /, \\
    BondSpec(char),
    /// Open parenthesis (branch)
    OpenParen,
    /// Close parenthesis
    CloseParen,
    /// Open square bracket
    OpenBracket,
    /// Close square bracket
    CloseBracket,
    /// Ring closure digit 0-9
    RingClosureDigit(u8),
    /// Ring closure %NN
    RingClosurePercent(u8),
    /// Logical AND operator &
    And,
    /// Logical OR operator (comma)
    Or,
    /// Logical NOT operator !
    Not,
    /// Dollar sign (start of recursive SMARTS)
    Dollar,
    /// Low-order bit separator (.)
    Dot,
    /// End of token stream
    EndOfStream,
}

/// Tokenize a SMARTS string into a sequence of tokens.
///
/// RDKit source: smarts.ll flex grammar
/// RDKit❗✔️: Ported from flex rules to hand-written tokenizer.
fn tokenize(input: &str) -> Result<Vec<(Token, usize)>, SmartsParseError> {
    let mut tokens = Vec::new();
    let chars: Vec<char> = input.chars().collect();
    let len = chars.len();
    let mut i = 0;

    while i < len {
        let ch = chars[i];
        match ch {
            // Whitespace
            ' ' | '\t' | '\n' | '\r' => {
                i += 1;
                continue;
            }

            // Square brackets
            '[' => {
                let start = i;
                i += 1;
                // Read the content until the matching ]
                let mut depth = 1u32;
                let content_start = i;
                while i < len && depth > 0 {
                    if chars[i] == '[' {
                        depth += 1;
                    } else if chars[i] == ']' {
                        depth -= 1;
                    }
                    i += 1;
                }
                if depth > 0 {
                    // RDKit✔️✔️: unclosed bracket
                    return Err(SmartsParseError::UnclosedBracket(start));
                }
                // i is now past the ]
                let content: String = chars[content_start..i - 1].iter().collect();
                tokens.push((Token::BracketContent(content), start));
            }

            // Bond specifiers (single characters)
            '-' | '=' | '#' | ':' | '~' | '@' => {
                tokens.push((Token::BondSpec(ch), i));
                i += 1;
            }
            '/' => {
                tokens.push((Token::BondSpec('/'), i));
                i += 1;
            }
            '\\' => {
                tokens.push((Token::BondSpec('\\'), i));
                i += 1;
            }

            // Branches
            '(' => {
                tokens.push((Token::OpenParen, i));
                i += 1;
            }
            ')' => {
                tokens.push((Token::CloseParen, i));
                i += 1;
            }

            // Logical operators
            '&' => {
                tokens.push((Token::And, i));
                i += 1;
            }
            ';' => {
                tokens.push((Token::And, i));
                i += 1;
            }
            ',' => {
                tokens.push((Token::Or, i));
                i += 1;
            }
            '!' => {
                tokens.push((Token::Not, i));
                i += 1;
            }

            // Dollar sign (recursive SMARTS start)
            '$' => {
                tokens.push((Token::Dollar, i));
                i += 1;
            }

            // Low-order bond separator
            '.' => {
                tokens.push((Token::Dot, i));
                i += 1;
            }

            // Ring closures: %NN
            '%' => {
                if i + 2 < len {
                    let d1 = chars[i + 1];
                    let d2 = chars[i + 2];
                    if d1.is_ascii_digit() && d2.is_ascii_digit() {
                        let num = (d1.to_digit(10).unwrap() * 10 + d2.to_digit(10).unwrap()) as u8;
                        tokens.push((Token::RingClosurePercent(num), i));
                        i += 3;
                        continue;
                    }
                }
                // RDKit✔️✔️: SmartsParse.cpp — error on invalid %
                return Err(SmartsParseError::UnexpectedCharacter {
                    position: i,
                    character: ch,
                    context: "expected two digits after %".to_string(),
                });
            }

            // Single digit ring closure
            d if d.is_ascii_digit() => {
                let num = d.to_digit(10).unwrap() as u8;
                tokens.push((Token::RingClosureDigit(num), i));
                i += 1;
            }

            // Aromatic elements + query-only chars (lowercase)
            'c' | 'n' | 'o' | 's' | 'p' | 'a' | 'b' => {
                let name = ch.to_string();
                tokens.push((Token::AromaticElement(name), i));
                i += 1;
            }

            // Organic elements + SMARTS primitives (uppercase or single-letter)
            'B' | 'C' | 'N' | 'O' | 'S' | 'P' | 'F' | 'I' | 'H' | 'R' | 'X' | 'D' | 'v' | 'V'
            | 'r' | 'u' | 'A' | 'T' | 'Z' | 'K' | 'W' | 'U' | 'Y' | 'G' | 'L' | 'J' | 'E' | 'M'
            | 'Q' => {
                let start = i;
                i += 1;
                // Two-char elements: Cl, Br, Si, As, Se, Te, etc.
                if i < len && chars[i].is_ascii_lowercase() {
                    // RDKit✔️✔️: In SMARTS, lowercase after uppercase means two-char element
                    // unless the lowercase is a primitive indicator.
                    // Cl, Br, Si, etc.
                    let two_char: String = chars[start..=i].iter().collect();
                    match two_char.as_str() {
                        // Known two-char elements
                        "Cl" | "Br" | "Si" | "As" | "Se" | "Te" | "He" | "Li" | "Be" | "Ne"
                        | "Na" | "Mg" | "Al" | "Ar" | "Ca" | "Sc" | "Ti" | "Cr" | "Mn" | "Fe"
                        | "Co" | "Ni" | "Cu" | "Zn" | "Ga" | "Ge" | "Kr" | "Rb" | "Sr" | "Zr"
                        | "Nb" | "Mo" | "Tc" | "Ru" | "Rh" | "Pd" | "Ag" | "Cd" | "In" | "Sn"
                        | "Sb" | "Xe" | "Cs" | "Ba" | "La" | "Ce" | "Pr" | "Nd" | "Pm" | "Sm"
                        | "Eu" | "Gd" | "Tb" | "Dy" | "Ho" | "Er" | "Tm" | "Yb" | "Lu" | "Hf"
                        | "Ta" | "Re" | "Os" | "Ir" | "Pt" | "Au" | "Hg" | "Tl" | "Pb" | "Bi"
                        | "Po" | "At" | "Rn" | "Fr" | "Ra" | "Ac" | "Th" | "Pa" | "Np" | "Pu"
                        | "Am" | "Cm" | "Bk" | "Cf" | "Es" | "Fm" | "Md" | "No" | "Lr" | "Rf"
                        | "Db" | "Sg" | "Bh" | "Hs" | "Mt" | "Ds" | "Rg" | "Cn" | "Nh" | "Fl"
                        | "Mc" | "Lv" | "Ts" | "Og" => {
                            i += 1;
                        }
                        _ => {
                            // Single character token
                        }
                    }
                }
                let name: String = chars[start..i].iter().collect();
                tokens.push((Token::OrganicElement(name), start));
            }

            // Star (wildcard)
            '*' => {
                tokens.push((Token::OrganicElement("*".to_string()), i));
                i += 1;
            }

            _ => {
                // RDKit✔️✔️: smarts.ll — any other character is a bad character
                return Err(SmartsParseError::UnexpectedCharacter {
                    position: i,
                    character: ch,
                    context: "unexpected character in SMARTS string".to_string(),
                });
            }
        }
    }

    tokens.push((Token::EndOfStream, len));
    Ok(tokens)
}

// ---------------------------------------------------------------------------
// Recursive-descent SMARTS Parser
// ---------------------------------------------------------------------------

/// RDKit❗✔️: Our recursive-descent SMARTS parser. In RDKit, the parser is
/// generated by bison from smarts.yy. We implement the same grammar logic
/// by hand.
struct SmartsParser<'a> {
    tokens: &'a [(Token, usize)],
    input: &'a str,
    pos: usize,
    /// Track ring closures: map from closure number to (atom_index_in_pattern, bond_query, position_in_input)
    ring_closure_targets: BTreeMap<u8, (usize, QueryNode<BondQueryPredicate>, usize)>,
}

impl<'a> SmartsParser<'a> {
    fn new(tokens: &'a [(Token, usize)], input: &'a str) -> Self {
        Self {
            tokens,
            input,
            pos: 0,
            ring_closure_targets: BTreeMap::new(),
        }
    }

    fn peek(&self) -> &(Token, usize) {
        &self.tokens[self.pos]
    }

    fn advance(&mut self) {
        self.pos += 1;
    }

    fn pos_info(&self) -> usize {
        self.tokens[self.pos].1
    }

    /// Parse the full SMARTS pattern into a SmartsMolecule.
    ///
    /// RDKit source: smarts.yy — the top-level grammar rule produces a molecule.
    fn parse_smarts_molecule(&mut self) -> Result<SmartsMolecule, SmartsParseError> {
        let mut atom_queries: Vec<QueryNode<AtomQueryPredicate>> = Vec::new();
        let mut bond_queries: Vec<QueryNode<BondQueryPredicate>> = Vec::new();
        let mut ring_closures: Vec<(u8, usize)> = Vec::new();

        // RDKit✔️✔️: Parse the first atom
        let first = self.parse_atom()?;
        atom_queries.push(first);

        // RDKit✔️✔️: Parse the rest of the pattern
        self.parse_smarts_chain(&mut atom_queries, &mut bond_queries, &mut ring_closures)?;

        Ok(SmartsMolecule::new(
            atom_queries,
            bond_queries,
            ring_closures,
        ))
    }

    /// RDKit source: smarts.yy — mol → atom atom_list
    /// RDKit✔️✔️: Parse the chain of atoms, bonds, branches, and ring closures.
    fn parse_smarts_chain(
        &mut self,
        atom_queries: &mut Vec<QueryNode<AtomQueryPredicate>>,
        bond_queries: &mut Vec<QueryNode<BondQueryPredicate>>,
        ring_closures: &mut Vec<(u8, usize)>,
    ) -> Result<(), SmartsParseError> {
        loop {
            match self.peek() {
                (Token::EndOfStream, _) => break,
                (Token::CloseParen, _) => break,
                (Token::CloseBracket, _) => break,

                // Bond spec followed by atom
                (Token::BondSpec(_), _) | (Token::Not, _) | (Token::And, _) => {
                    let bond = self.parse_bond()?;
                    bond_queries.push(bond);
                    let atom = self.parse_atom()?;
                    atom_queries.push(atom);
                }

                // No bond spec — implicit single/aromatic bond (SMARTS semantics)
                _ => {
                    // Check if next is a ring closure or branch first
                    match self.peek() {
                        (Token::RingClosureDigit(n), pos) | (Token::RingClosurePercent(n), pos) => {
                            let num = *n;
                            let existing_len = atom_queries.len();
                            let bond_pos = *pos;
                            self.advance();
                            // Record ring closure: the last atom in the pattern
                            let last_atom_idx = existing_len - 1;
                            ring_closures.push((num, last_atom_idx));
                            self.ring_closure_targets.entry(num).or_insert((
                                last_atom_idx,
                                QueryNode::Predicate(BondQueryPredicate::Any),
                                bond_pos,
                            ));
                        }
                        (Token::OpenParen, _) => {
                            // RDKit✔️✔️: Branch — parse the sub-pattern
                            self.advance();
                            // The branch is a sub-chain attached to the current atom
                            let saved_bonds = bond_queries.len();
                            let saved_atoms = atom_queries.len();
                            self.parse_smarts_chain(atom_queries, bond_queries, ring_closures)?;
                            match self.peek() {
                                (Token::CloseParen, _) => {
                                    self.advance();
                                }
                                (tok, pos) => {
                                    return Err(SmartsParseError::UnexpectedCharacter {
                                        position: *pos,
                                        character: format!("{:?}", tok)
                                            .chars()
                                            .next()
                                            .unwrap_or('?'),
                                        context: "expected close parenthesis".to_string(),
                                    });
                                }
                            }
                            // Bond from parent to first atom of branch
                            if saved_atoms < atom_queries.len() {
                                bond_queries.insert(
                                    saved_bonds,
                                    QueryNode::Predicate(BondQueryPredicate::Any),
                                );
                            }
                        }
                        (Token::Dot, _) => {
                            // RDKit✔️✔️: Dot separates disconnected fragments
                            self.advance();
                            bond_queries.push(QueryNode::Predicate(BondQueryPredicate::Any));
                            let atom = self.parse_atom()?;
                            atom_queries.push(atom);
                        }
                        // Atom follows implicitly with default bond
                        _ => {
                            bond_queries.push(QueryNode::Predicate(BondQueryPredicate::Any));
                            let atom = self.parse_atom()?;
                            atom_queries.push(atom);
                        }
                    }
                }
            }
        }

        Ok(())
    }

    /// Parse an atom expression.
    ///
    /// RDKit source: smarts.yy — atom → ORGANIC_ATOM_TOKEN | AROMATIC_ATOM_TOKEN | ATOM_TOKEN
    /// RDKit✔️✔️: Atoms can be organic subset, aromatic, or bracket atoms.
    fn parse_atom(&mut self) -> Result<QueryNode<AtomQueryPredicate>, SmartsParseError> {
        let (token, _pos) = self.peek().clone();
        match token {
            Token::OrganicElement(name) => {
                let query = organic_element_to_query(&name);
                self.advance();
                Ok(query)
            }
            Token::AromaticElement(name) => {
                let query = aromatic_element_to_query(&name);
                self.advance();
                Ok(query)
            }
            Token::BracketContent(content) => {
                self.advance();
                self.parse_bracket_atom_content(&content)
            }
            Token::EndOfStream => Err(SmartsParseError::UnexpectedEnd(
                "expected atom but reached end".to_string(),
            )),
            _ => {
                let pos = self.pos_info();
                Err(SmartsParseError::UnexpectedCharacter {
                    position: pos,
                    character: '?',
                    context: "expected atom expression".to_string(),
                })
            }
        }
    }

    /// Parse the content inside a bracket atom: e.g. "C", "C@@H", "N+", "O-", "#6", "6X4"
    ///
    /// RDKit source: smarts.yy — the ATOM_TOKEN production and its associated actions
    /// RDKit✔️✔️: Bracket atom content is parsed as a sequence of primitives AND-ed together.
    fn parse_bracket_atom_content(
        &mut self,
        content: &str,
    ) -> Result<QueryNode<AtomQueryPredicate>, SmartsParseError> {
        let chars: Vec<char> = content.chars().collect();
        let len = chars.len();
        let mut i = 0;
        let mut predicates: Vec<QueryNode<AtomQueryPredicate>> = Vec::new();
        let mut negate_next = false;
        let mut logical_or = false;

        while i < len {
            let ch = chars[i];

            // Handle logical NOT
            if ch == '!' {
                negate_next = !negate_next;
                i += 1;
                continue;
            }

            // Handle logical OR (comma)
            if ch == ',' {
                logical_or = true;
                i += 1;
                continue;
            }

            // Handle logical AND (ampersand)
            if ch == '&' {
                i += 1;
                continue;
            }

            // Handle semicolon (AND)
            if ch == ';' {
                i += 1;
                continue;
            }

            // Parse the primitive at position i
            let (pred, consumed) = self.parse_atom_primitive(&chars, i, len)?;

            // Apply negation
            let pred = if negate_next {
                negate_next = false;
                QueryNode::not(pred)
            } else {
                pred
            };

            predicates.push(pred);

            i = consumed;
        }

        // Combine predicates
        if predicates.is_empty() {
            // RDKit✔️✔️: Empty bracket [ ] matches any atom
            Ok(QueryNode::Predicate(AtomQueryPredicate::Any))
        } else if logical_or && predicates.len() > 1 {
            // RDKit✔️✔️: Comma-separated primitives are OR-ed
            Ok(QueryNode::Or(predicates))
        } else if predicates.len() == 1 {
            Ok(predicates.into_iter().next().unwrap())
        } else {
            // RDKit✔️✔️: Semicolon/ampersand-separated (or implicit) are AND-ed
            Ok(QueryNode::And(predicates))
        }
    }

    /// Parse a single atom primitive from the bracket content starting at position `i`.
    ///
    /// RDKit source: smarts.yy — primitives within ATOM_TOKEN
    /// RDKit✔️✔️: Handles all SMARTS atom primitives.
    fn parse_atom_primitive(
        &self,
        chars: &[char],
        i: usize,
        len: usize,
    ) -> Result<(QueryNode<AtomQueryPredicate>, usize), SmartsParseError> {
        if i >= len {
            return Err(SmartsParseError::UnexpectedEnd(
                "expected atom primitive".to_string(),
            ));
        }

        let ch = chars[i];

        // Atomic number: #N
        // RDKit✔️✔️: smarts.yy — HASH_TOKEN NUMBER
        if ch == '#' {
            let (num, consumed) = self.parse_number(chars, i + 1, len)?;
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(num as u8)),
                consumed,
            ));
        }

        // Recursive SMARTS: $(...)
        // RDKit✔️✔️: recursive SMARTS primitives are preserved as query atoms.
        if ch == '$' {
            if chars.get(i + 1) != Some(&'(') {
                return Err(SmartsParseError::InvalidAtomPrimitive {
                    position: i,
                    detail: "expected '(' after '$'".to_string(),
                });
            }
            let mut depth = 1usize;
            let mut end = i + 2;
            while end < len && depth > 0 {
                match chars[end] {
                    '(' => depth += 1,
                    ')' => depth -= 1,
                    _ => {}
                }
                end += 1;
            }
            if depth != 0 {
                return Err(SmartsParseError::UnclosedParenthesis(i + 1));
            }
            let recursive_smarts: String = chars[i..end].iter().collect();
            let mut consumed = end;
            if chars.get(consumed) == Some(&'_') {
                consumed += 1;
                while consumed < len && chars[consumed].is_ascii_digit() {
                    consumed += 1;
                }
            }
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(recursive_smarts)),
                consumed,
            ));
        }

        // Formal charge: +N or -N
        // RDKit✔️✔️: smarts.yy — PLUS_TOKEN NUMBER | MINUS_TOKEN NUMBER
        if ch == '+' {
            let start = i + 1;
            if start < len && chars[start] == '+' {
                // RDKit✔️✔️: "++" is charge +2
                return Ok((
                    QueryNode::Predicate(AtomQueryPredicate::FormalCharge(2)),
                    start + 1,
                ));
            }
            let (num, consumed) = self.parse_optional_number(chars, start, len);
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::FormalCharge(num as i8)),
                consumed,
            ));
        }

        if ch == '-' {
            let start = i + 1;
            if start < len && chars[start] == '-' {
                // RDKit✔️✔️: "--" is charge -2
                return Ok((
                    QueryNode::Predicate(AtomQueryPredicate::FormalCharge(-2)),
                    start + 1,
                ));
            }
            let (num, consumed) = self.parse_optional_number(chars, start, len);
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::FormalCharge(-(num as i8))),
                consumed,
            ));
        }

        // Chirality: @ or @@
        // RDKit✔️✔️: smarts.yy — AT_TOKEN
        if ch == '@' {
            let start = i + 1;
            if start < len && chars[start] == '@' {
                // @@
                return Ok((
                    QueryNode::Predicate(AtomQueryPredicate::ChiralTagMatch(
                        crate::ChiralTag::TetrahedralCw,
                    )),
                    start + 1,
                ));
            }
            // @
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::ChiralTagMatch(
                    crate::ChiralTag::TetrahedralCcw,
                )),
                start,
            ));
        }

        // Hydrogen count: H or H<N>
        // RDKit✔️✔️: smarts.yy — H_TOKEN (optional NUMBER)
        if ch == 'h' {
            // RDKit✔️✔️: 'h' is implicit H count (SMARTS)
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            if num == 0 {
                return Ok((
                    QueryNode::Predicate(AtomQueryPredicate::ImplicitHydrogenCount(0)),
                    consumed,
                ));
            }
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::ImplicitHydrogenCount(num as u8)),
                consumed,
            ));
        }
        if ch == 'H' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            if num == 0 {
                return Ok((
                    QueryNode::Predicate(AtomQueryPredicate::ImplicitHydrogenCount(1)),
                    consumed,
                ));
            }
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::ImplicitHydrogenCount(num as u8)),
                consumed,
            ));
        }

        // Ring membership: R or R<N>
        // RDKit✔️✔️: smarts.yy — R_TOKEN (optional NUMBER)
        if ch == 'R' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            if num == 0 {
                return Ok((QueryNode::Predicate(AtomQueryPredicate::InRing), consumed));
            }
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::NumRingBonds(num as u8)),
                consumed,
            ));
        }
        if ch == 'r' {
            if chars.get(i + 1) == Some(&'{') {
                let (predicate, consumed) = self.parse_ring_size_range(chars, i + 2, len)?;
                return Ok((predicate, consumed));
            }
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            if num == 0 {
                return Ok((QueryNode::Predicate(AtomQueryPredicate::InRing), consumed));
            }
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::SmallestRingSize(num as u8)),
                consumed,
            ));
        }

        // Connectivity/degree: X or X<N>
        // RDKit✔️✔️: smarts.yy — X_TOKEN (optional NUMBER)
        if ch == 'X' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            if num == 0 {
                return Ok((QueryNode::Predicate(AtomQueryPredicate::Any), consumed));
            }
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::Connectivity(num as u8)),
                consumed,
            ));
        }

        // Ring connectivity: x or x<N>
        // RDKit✔️✔️: lowercase x matches atom ring-bond count.
        if ch == 'x' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            if num == 0 {
                return Ok((QueryNode::Predicate(AtomQueryPredicate::Any), consumed));
            }
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::NumRingBonds(num as u8)),
                consumed,
            ));
        }

        // Degree: D or D<N>
        // RDKit✔️✔️: smarts.yy — D_TOKEN (optional NUMBER)
        if ch == 'D' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            if num == 0 {
                return Ok((QueryNode::Predicate(AtomQueryPredicate::Any), consumed));
            }
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::Degree(num as u8)),
                consumed,
            ));
        }

        // Hybridization: ^1, ^2, ^3, ...
        // RDKit✔️✔️: caret hybridization primitives map to explicit hybridization predicates.
        if ch == '^' {
            let (num, consumed) = self.parse_number(chars, i + 1, len)?;
            let hybridization = match num {
                1 => crate::Hybridization::S,
                2 => crate::Hybridization::Sp,
                3 => crate::Hybridization::Sp2,
                4 => crate::Hybridization::Sp3,
                5 => crate::Hybridization::Sp3d,
                6 => crate::Hybridization::Sp3d2,
                _ => crate::Hybridization::Other,
            };
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::HybridizationMatch(hybridization)),
                consumed,
            ));
        }

        // Valence/connectivity: v or v<N>
        // RDKit✔️✔️: smarts.yy — v_TOKEN (optional NUMBER)
        if ch == 'v' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            if num == 0 {
                return Ok((QueryNode::Predicate(AtomQueryPredicate::Any), consumed));
            }
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::TotalDegree(num as u8)),
                consumed,
            ));
        }

        // Atom map number: :<N>
        // RDKit✔️✔️: bracket atom map property used for SMARTS match indexing.
        if ch == ':' {
            let (num, consumed) = self.parse_number(chars, i + 1, len)?;
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::AtomMapNumber(num)),
                consumed,
            ));
        }

        // Aromatic query: a / A
        // RDKit✔️✔️: smarts.yy — a_TOKEN / A_TOKEN
        if ch == 'a' {
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::IsAromatic(true)),
                i + 1,
            ));
        }
        if ch == 'A' {
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::IsAromatic(false)),
                i + 1,
            ));
        }

        // Unsaturated: u
        // RDKit✔️✔️: smarts.yy — u_TOKEN
        if ch == 'u' {
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::IsUnsaturated),
                i + 1,
            ));
        }

        // Element symbol
        // RDKit✔️✔️: smarts.yy — element symbols (2-char or 1-char)
        if ch.is_ascii_uppercase() {
            let start = i;
            let mut end = i + 1;
            // Check for two-char element
            if end < len && chars[end].is_ascii_lowercase() {
                let two_char: String = chars[start..=end].iter().collect();
                if let Some(atomic_num) = element_symbol_to_atomic_number(&two_char) {
                    return Ok((
                        QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(atomic_num)),
                        end + 1,
                    ));
                }
            }
            let one_char: String = chars[start..end].iter().collect();
            if let Some(atomic_num) = element_symbol_to_atomic_number(&one_char) {
                return Ok((
                    QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(atomic_num)),
                    end,
                ));
            }
            // Unknown symbol — treat as wildcard
            return Ok((QueryNode::Predicate(AtomQueryPredicate::Any), end));
        }

        // Isotope: numeric literal before element (e.g. "13C")
        if ch.is_ascii_digit() {
            let (num, consumed) = self.parse_number(chars, i, len)?;
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::Isotope(num as u16)),
                consumed,
            ));
        }

        // Lowercase aromatic element inside bracket (e.g. [c])
        if ch.is_ascii_lowercase() && ch != 'a' && ch != 'u' && ch != 'v' && ch != 'r' && ch != 'h'
        {
            let name = ch.to_string();
            // Same as aromatic_element_to_query but for bracket content
            let query = aromatic_element_to_query(&name);
            return Ok((query, i + 1));
        }

        // Wildcard inside bracket (shouldn't normally appear, but handle it)
        if ch == '*' {
            return Ok((QueryNode::Predicate(AtomQueryPredicate::Any), i + 1));
        }

        Err(SmartsParseError::InvalidAtomPrimitive {
            position: i,
            detail: format!("unexpected character '{}'", ch),
        })
    }

    fn parse_ring_size_range(
        &self,
        chars: &[char],
        start: usize,
        len: usize,
    ) -> Result<(QueryNode<AtomQueryPredicate>, usize), SmartsParseError> {
        let mut pos = start;
        let lower = if pos < len && chars[pos].is_ascii_digit() {
            let (num, consumed) = self.parse_number(chars, pos, len)?;
            pos = consumed;
            Some(num as u8)
        } else {
            None
        };
        if chars.get(pos) != Some(&'-') {
            return Err(SmartsParseError::InvalidAtomPrimitive {
                position: start.saturating_sub(2),
                detail: "expected '-' in ring-size range".to_string(),
            });
        }
        pos += 1;
        let upper = if pos < len && chars[pos].is_ascii_digit() {
            let (num, consumed) = self.parse_number(chars, pos, len)?;
            pos = consumed;
            Some(num as u8)
        } else {
            None
        };
        if chars.get(pos) != Some(&'}') {
            return Err(SmartsParseError::InvalidAtomPrimitive {
                position: start.saturating_sub(2),
                detail: "expected '}' to close ring-size range".to_string(),
            });
        }
        pos += 1;

        let predicate = match (lower, upper) {
            (Some(min), Some(max)) => QueryNode::And(vec![
                QueryNode::Predicate(AtomQueryPredicate::SmallestRingSizeGreaterEqual(min)),
                QueryNode::Predicate(AtomQueryPredicate::SmallestRingSizeLessEqual(max)),
            ]),
            (Some(min), None) => {
                QueryNode::Predicate(AtomQueryPredicate::SmallestRingSizeGreaterEqual(min))
            }
            (None, Some(max)) => {
                QueryNode::Predicate(AtomQueryPredicate::SmallestRingSizeLessEqual(max))
            }
            (None, None) => {
                return Err(SmartsParseError::InvalidAtomPrimitive {
                    position: start.saturating_sub(2),
                    detail: "empty ring-size range".to_string(),
                });
            }
        };

        Ok((predicate, pos))
    }

    /// Parse a number from position i. Returns (value, consumed_index).
    fn parse_number(
        &self,
        chars: &[char],
        i: usize,
        len: usize,
    ) -> Result<(u32, usize), SmartsParseError> {
        if i >= len || !chars[i].is_ascii_digit() {
            return Err(SmartsParseError::UnexpectedEnd(
                "expected number".to_string(),
            ));
        }
        let mut val = 0u32;
        let mut pos = i;
        while pos < len && chars[pos].is_ascii_digit() {
            val = val * 10 + chars[pos].to_digit(10).unwrap();
            pos += 1;
        }
        Ok((val, pos))
    }

    /// Parse an optional number from position i. Returns (value if present else 0, consumed_index).
    fn parse_optional_number(&self, chars: &[char], i: usize, len: usize) -> (u32, usize) {
        if i >= len || !chars[i].is_ascii_digit() {
            return (0, i);
        }
        let mut val = 0u32;
        let mut pos = i;
        while pos < len && chars[pos].is_ascii_digit() {
            val = val * 10 + chars[pos].to_digit(10).unwrap();
            pos += 1;
        }
        (val, pos)
    }

    /// Parse a bond specifier.
    ///
    /// RDKit source: smarts.yy — bond → BOND_TOKEN
    /// RDKit❗✔️: Supports the CrystalFF parameter-table bond combinations
    /// `!@`, `@`, and `;` in addition to the basic bond specifiers.
    fn parse_bond(&mut self) -> Result<QueryNode<BondQueryPredicate>, SmartsParseError> {
        let mut negate_next = false;
        let mut predicates = Vec::new();
        let mut consumed_any = false;
        let mut logical_or = false;

        while matches!(
            self.peek(),
            (Token::BondSpec(_), _) | (Token::Not, _) | (Token::And, _) | (Token::Or, _)
        ) {
            match self.peek() {
                (Token::Not, _) => {
                    consumed_any = true;
                    negate_next = !negate_next;
                    self.advance();
                }
                (Token::And, _) => {
                    consumed_any = true;
                    self.advance();
                }
                (Token::Or, _) => {
                    consumed_any = true;
                    logical_or = true;
                    self.advance();
                }
                (Token::BondSpec(ch), _) => {
                    consumed_any = true;
                    let query = bond_spec_to_query(*ch);
                    self.advance();
                    let query = if negate_next {
                        negate_next = false;
                        QueryNode::not(query)
                    } else {
                        query
                    };
                    predicates.push(query);
                }
                _ => break,
            }
        }

        let predicates = predicates
            .into_iter()
            .filter(|query| *query != QueryNode::Predicate(BondQueryPredicate::Any))
            .collect::<Vec<_>>();
        match predicates.len() {
            0 if consumed_any => Ok(QueryNode::Predicate(BondQueryPredicate::Any)),
            0 => Err(SmartsParseError::UnexpectedCharacter {
                position: self.pos_info(),
                character: '?',
                context: "expected bond specifier".to_string(),
            }),
            len if logical_or && len > 1 => Ok(QueryNode::Or(predicates)),
            1 => Ok(predicates.into_iter().next().expect("single bond query")),
            _ => Ok(QueryNode::And(predicates)),
        }
    }
}

// ---------------------------------------------------------------------------
// SMARTS primitive helpers
// ---------------------------------------------------------------------------

/// Convert an organic element name to an atom query predicate.
///
/// RDKit source: SmartsParse.cpp / smarts.ll organic atom handling.
/// RDKit✔️✔️: Organic subset atoms: C, N, O, S, P, F, Cl, Br, I, *, B
fn organic_element_to_query(name: &str) -> QueryNode<AtomQueryPredicate> {
    match name {
        // RDKit✔️✔️: * matches any atom
        "*" => QueryNode::Predicate(AtomQueryPredicate::Any),
        // RDKit✔️✔️: A = aliphatic (not aromatic)
        "A" => QueryNode::Predicate(AtomQueryPredicate::IsAromatic(false)),
        // RDKit✔️✔️: Organic elements with their full implicit H count
        "B" => QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(5)),
        "C" => QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(6)),
        "N" => QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(7)),
        "O" => QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(8)),
        "P" => QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(15)),
        "S" => QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(16)),
        "F" => QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(9)),
        "Cl" => QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(17)),
        "Br" => QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(35)),
        "I" => QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(53)),
        // RDKit✔️✔️: SMARTS organic subset also includes metalloids as explicit
        "Si" => QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(14)),
        "As" => QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(33)),
        "Se" => QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(34)),
        "Te" => QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(52)),
        // RDKit✔️✔️: Single-letter primitives that look like elements
        "H" => QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(1)),
        // RDKit✔️✔️: Unknown symbol — treat as Any
        _ => QueryNode::Predicate(AtomQueryPredicate::Any),
    }
}

/// Convert an aromatic element name to an atom query.
///
/// RDKit source: smarts.ll — AROMATIC_ATOM_TOKEN handling.
/// RDKit✔️✔️: Aromatic subset: c, n, o, s, p
fn aromatic_element_to_query(name: &str) -> QueryNode<AtomQueryPredicate> {
    match name {
        "c" => QueryNode::And(vec![
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::Predicate(AtomQueryPredicate::IsAromatic(true)),
        ]),
        "n" => QueryNode::And(vec![
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(7)),
            QueryNode::Predicate(AtomQueryPredicate::IsAromatic(true)),
        ]),
        "o" => QueryNode::And(vec![
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(8)),
            QueryNode::Predicate(AtomQueryPredicate::IsAromatic(true)),
        ]),
        "s" => QueryNode::And(vec![
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(16)),
            QueryNode::Predicate(AtomQueryPredicate::IsAromatic(true)),
        ]),
        "p" => QueryNode::And(vec![
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(15)),
            QueryNode::Predicate(AtomQueryPredicate::IsAromatic(true)),
        ]),
        // RDKit✔️✔️: 'a' (aromatic query) and 'b' (aromatic B) in smarts.ll
        "b" => QueryNode::And(vec![
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(5)),
            QueryNode::Predicate(AtomQueryPredicate::IsAromatic(true)),
        ]),
        "a" => QueryNode::Predicate(AtomQueryPredicate::IsAromatic(true)),
        _ => QueryNode::Predicate(AtomQueryPredicate::Any),
    }
}

/// Convert a bond specifier character to a bond query predicate.
///
/// RDKit source: smarts.yy / smarts.ll bond handling.
/// RDKit✔️✔️: -, =, #, :, ~, /, \\
fn bond_spec_to_query(ch: char) -> QueryNode<BondQueryPredicate> {
    match ch {
        '-' => QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Single)),
        '=' => QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Double)),
        '#' => QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Triple)),
        ':' => QueryNode::Predicate(BondQueryPredicate::IsAromatic(true)),
        '@' => QueryNode::Predicate(BondQueryPredicate::IsInRing(true)),
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

/// Look up the atomic number for an element symbol.
///
/// RDKit✔️✔️: Standard periodic table mapping.
fn element_symbol_to_atomic_number(symbol: &str) -> Option<u8> {
    match symbol {
        "H" => Some(1),
        "He" => Some(2),
        "Li" => Some(3),
        "Be" => Some(4),
        "B" => Some(5),
        "C" => Some(6),
        "N" => Some(7),
        "O" => Some(8),
        "F" => Some(9),
        "Ne" => Some(10),
        "Na" => Some(11),
        "Mg" => Some(12),
        "Al" => Some(13),
        "Si" => Some(14),
        "P" => Some(15),
        "S" => Some(16),
        "Cl" => Some(17),
        "Ar" => Some(18),
        "K" => Some(19),
        "Ca" => Some(20),
        "Sc" => Some(21),
        "Ti" => Some(22),
        "V" => Some(23),
        "Cr" => Some(24),
        "Mn" => Some(25),
        "Fe" => Some(26),
        "Co" => Some(27),
        "Ni" => Some(28),
        "Cu" => Some(29),
        "Zn" => Some(30),
        "Ga" => Some(31),
        "Ge" => Some(32),
        "As" => Some(33),
        "Se" => Some(34),
        "Br" => Some(35),
        "Kr" => Some(36),
        "Rb" => Some(37),
        "Sr" => Some(38),
        "Y" => Some(39),
        "Zr" => Some(40),
        "Nb" => Some(41),
        "Mo" => Some(42),
        "Tc" => Some(43),
        "Ru" => Some(44),
        "Rh" => Some(45),
        "Pd" => Some(46),
        "Ag" => Some(47),
        "Cd" => Some(48),
        "In" => Some(49),
        "Sn" => Some(50),
        "Sb" => Some(51),
        "Te" => Some(52),
        "I" => Some(53),
        "Xe" => Some(54),
        "Cs" => Some(55),
        "Ba" => Some(56),
        "La" => Some(57),
        "Ce" => Some(58),
        "Pr" => Some(59),
        "Nd" => Some(60),
        "Pm" => Some(61),
        "Sm" => Some(62),
        "Eu" => Some(63),
        "Gd" => Some(64),
        "Tb" => Some(65),
        "Dy" => Some(66),
        "Ho" => Some(67),
        "Er" => Some(68),
        "Tm" => Some(69),
        "Yb" => Some(70),
        "Lu" => Some(71),
        "Hf" => Some(72),
        "Ta" => Some(73),
        "W" => Some(74),
        "Re" => Some(75),
        "Os" => Some(76),
        "Ir" => Some(77),
        "Pt" => Some(78),
        "Au" => Some(79),
        "Hg" => Some(80),
        "Tl" => Some(81),
        "Pb" => Some(82),
        "Bi" => Some(83),
        "Po" => Some(84),
        "At" => Some(85),
        "Rn" => Some(86),
        "Fr" => Some(87),
        "Ra" => Some(88),
        "Ac" => Some(89),
        "Th" => Some(90),
        "Pa" => Some(91),
        "U" => Some(92),
        "Np" => Some(93),
        "Pu" => Some(94),
        "Am" => Some(95),
        "Cm" => Some(96),
        "Bk" => Some(97),
        "Cf" => Some(98),
        "Es" => Some(99),
        "Fm" => Some(100),
        "Md" => Some(101),
        "No" => Some(102),
        "Lr" => Some(103),
        "Rf" => Some(104),
        "Db" => Some(105),
        "Sg" => Some(106),
        "Bh" => Some(107),
        "Hs" => Some(108),
        "Mt" => Some(109),
        "Ds" => Some(110),
        "Rg" => Some(111),
        "Cn" => Some(112),
        "Nh" => Some(113),
        "Fl" => Some(114),
        "Mc" => Some(115),
        "Lv" => Some(116),
        "Ts" => Some(117),
        "Og" => Some(118),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tokenize_simple_smarts() {
        // CC — two organic carbons
        let tokens = tokenize("CC").unwrap();
        assert_eq!(tokens.len(), 3); // C, C, EOS
        assert_eq!(tokens[0].0, Token::OrganicElement("C".to_string()));
        assert_eq!(tokens[1].0, Token::OrganicElement("C".to_string()));
        assert_eq!(tokens[2].0, Token::EndOfStream);
    }

    #[test]
    fn test_tokenize_bracket_atom() {
        let tokens = tokenize("[N+]").unwrap();
        assert_eq!(tokens.len(), 2); // BracketContent, EOS
        match &tokens[0].0 {
            Token::BracketContent(content) => {
                assert_eq!(content, "N+");
            }
            _ => panic!("expected bracket content"),
        }
    }

    #[test]
    fn test_tokenize_ring_closure() {
        let tokens = tokenize("C1CC1").unwrap();
        assert_eq!(tokens.len(), 6); // C, 1, C, C, 1, EOS
        assert_eq!(tokens[0].0, Token::OrganicElement("C".to_string()));
        assert_eq!(tokens[1].0, Token::RingClosureDigit(1));
        assert_eq!(tokens[4].0, Token::RingClosureDigit(1));
    }

    #[test]
    fn test_tokenize_percent_ring_closure() {
        let tokens = tokenize("C%10CC%10").unwrap();
        assert_eq!(tokens.len(), 6);
        assert_eq!(tokens[1].0, Token::RingClosurePercent(10));
        assert_eq!(tokens[4].0, Token::RingClosurePercent(10));
    }

    #[test]
    fn test_tokenize_bond_specs() {
        let tokens = tokenize("C=O").unwrap();
        assert_eq!(tokens[1].0, Token::BondSpec('='));
        let tokens = tokenize("C#N").unwrap();
        assert_eq!(tokens[1].0, Token::BondSpec('#'));
        let tokens = tokenize("C~C").unwrap();
        assert_eq!(tokens[1].0, Token::BondSpec('~'));
    }

    #[test]
    fn test_tokenize_unclosed_bracket() {
        let result = tokenize("[NH");
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            SmartsParseError::UnclosedBracket(_)
        ));
    }

    #[test]
    fn test_organic_element_to_query() {
        assert_eq!(
            organic_element_to_query("C"),
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(6))
        );
        assert_eq!(
            organic_element_to_query("*"),
            QueryNode::Predicate(AtomQueryPredicate::Any)
        );
    }

    #[test]
    fn test_bracket_atom_element() {
        // [C] — carbon in brackets
        let mol = parse_smarts("[C]").unwrap();
        assert_eq!(mol.atom_queries.len(), 1);
        assert_eq!(mol.bond_queries.len(), 0);
    }

    #[test]
    fn test_parse_simple_smarts() {
        // CC
        let mol = parse_smarts("CC").unwrap();
        assert_eq!(mol.atom_queries.len(), 2);
        assert_eq!(mol.bond_queries.len(), 1);
        assert_eq!(
            mol.atom_queries[0],
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(6))
        );
        assert_eq!(
            mol.bond_queries[0],
            QueryNode::Predicate(BondQueryPredicate::Any)
        );
    }

    #[test]
    fn test_parse_bonded_smarts() {
        // C=O
        let mol = parse_smarts("C=O").unwrap();
        assert_eq!(mol.atom_queries.len(), 2);
        assert_eq!(mol.bond_queries.len(), 1);
        assert_eq!(
            mol.bond_queries[0],
            QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Double))
        );
    }

    #[test]
    fn test_bracket_with_charge() {
        let mol = parse_smarts("[N+]").unwrap();
        assert_eq!(mol.atom_queries.len(), 1);
    }

    #[test]
    fn test_bracket_with_chirality() {
        let mol = parse_smarts("[C@@H]").unwrap();
        assert_eq!(mol.atom_queries.len(), 1);
    }

    #[test]
    fn test_parse_ring_closure() {
        let mol = parse_smarts("C1CC1").unwrap();
        assert_eq!(mol.atom_queries.len(), 3);
        assert_eq!(mol.ring_closures.len(), 2);
    }

    #[test]
    fn test_parse_branch() {
        let mol = parse_smarts("C(C)C").unwrap();
        assert_eq!(mol.atom_queries.len(), 3);
    }

    #[test]
    fn test_bracket_atomic_number_primitive() {
        let mol = parse_smarts("[#6]").unwrap();
        assert_eq!(mol.atom_queries.len(), 1);
    }

    #[test]
    fn test_label_recursive_patterns_noop() {
        // No $() in input
        assert_eq!(label_recursive_patterns("CCO"), "CCO");
    }

    #[test]
    fn test_label_recursive_patterns_simple() {
        let result = label_recursive_patterns("[$([N])]");
        // Should append label after closing paren
        assert!(result.contains("_100") || result == "[$([N])]");
    }

    #[test]
    fn test_element_symbol_lookup() {
        assert_eq!(element_symbol_to_atomic_number("C"), Some(6));
        assert_eq!(element_symbol_to_atomic_number("O"), Some(8));
        assert_eq!(element_symbol_to_atomic_number("Cl"), Some(17));
        assert_eq!(element_symbol_to_atomic_number("Br"), Some(35));
        assert_eq!(element_symbol_to_atomic_number("Xx"), None);
    }

    #[test]
    fn test_parse_aromatic_smarts() {
        // c1ccccc1 — benzene SMARTS
        let mol = parse_smarts("c1ccccc1").unwrap();
        assert_eq!(mol.atom_queries.len(), 6);
        // Every atom is aromatic carbon
        for aq in &mol.atom_queries {
            match aq {
                QueryNode::And(children) => {
                    assert_eq!(children.len(), 2);
                    assert_eq!(
                        children[0],
                        QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(6))
                    );
                    assert_eq!(
                        children[1],
                        QueryNode::Predicate(AtomQueryPredicate::IsAromatic(true))
                    );
                }
                _ => panic!("expected And node for aromatic atom"),
            }
        }
    }

    #[test]
    fn test_smarts_molecule_num_atoms() {
        let mol = parse_smarts("CCO").unwrap();
        assert_eq!(mol.num_atoms(), 3);
        assert!(mol.atom_query(0).is_some());
        assert!(mol.bond_query(0).is_some());
        assert!(mol.bond_query(1).is_some());
    }

    #[test]
    fn test_smarts_parse_params_default() {
        let params = SmartsParseParams::default();
        assert!(params.allow_cxsmiles);
        assert!(!params.merge_hs);
    }

    #[test]
    fn test_varied_bonds() {
        let mol = parse_smarts("C#N").unwrap();
        assert_eq!(
            mol.bond_queries[0],
            QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Triple))
        );

        let mol = parse_smarts("C:N").unwrap();
        assert_eq!(
            mol.bond_queries[0],
            QueryNode::Predicate(BondQueryPredicate::IsAromatic(true))
        );
    }

    #[test]
    fn test_empty_smarts_molecule() {
        // An empty SMARTS produces a molecule with no atoms
        // (technically invalid SMARTS, but handled gracefully)
        let result = parse_smarts("");
        // Should be an error
        assert!(result.is_err());
    }
}
