use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;
use std::sync::atomic::{AtomicBool, Ordering};

use crate::{
    AdjacencyList, Atom, AtomId, AtomQueryPredicate, AtomSpec, BondDirection, BondId, BondOrder,
    BondQueryPredicate, BondSpec, BondStereo, ChiralTag, Conformer3D, Element, Molecule,
    MoleculeBuilder, QueryNode, RemoveHsParams, SmilesParseError, StereoError, StereoGroup,
    StereoGroupKind, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind,
};

mod cx;
mod stereo;

pub(crate) use self::stereo::{
    assign_chiral_types_from_3d_for_testing, assign_chiral_types_from_bond_dirs,
    assign_stereochemistry_cleanup_subset, clear_all_bond_dir_flags,
    set_double_bond_neighbor_directions, set_double_bond_neighbor_directions_from_stereo,
};
use self::{cx::*, stereo::*};

// RDKit source-reproduction note:
//
// This module is the source-level alignment frame for the future SMILES parser
// port from `third_party/rdkit/Code/GraphMol/SmilesParse`. The functions below
// follow `dev/source_reproduction_protocol.md`. They intentionally contain
// RDKit C++ comments with unresolved status markers before any behavior is
// implemented. Do not upgrade markers or add semantic behavior without checking
// the corresponding RDKit source line and the local Rust code.

const SMILES_START_PROP: &str = "_SmilesStart";
const CXSMILES_BOND_IDX_PROP: &str = "_cxsmilesBondIdx";
const UNSPECIFIED_ORDER_PROP: &str = "_unspecifiedOrder";
static YYSMILES_DEBUG: AtomicBool = AtomicBool::new(false);

fn atom_spec_from_atom(atom: &Atom) -> AtomSpec {
    let mut spec = AtomSpec::new(atom.element())
        .with_formal_charge(atom.formal_charge())
        .with_explicit_hydrogens(atom.explicit_hydrogens())
        .with_chiral_tag(atom.chiral_tag())
        .with_aromatic(atom.is_aromatic())
        .with_no_implicit(atom.no_implicit())
        .with_radical_electrons(atom.radical_electrons())
        .with_hybridization(atom.hybridization());
    if let Some(chiral_permutation) = atom.chiral_permutation() {
        spec = spec.with_chiral_permutation(chiral_permutation);
    }
    if atom.unknown_stereo() {
        spec = spec.with_unknown_stereo(true);
    }
    if let Some(mol_parity) = atom.mol_parity() {
        spec = spec.with_mol_parity(mol_parity);
    }
    if let Some(mol_inversion_flag) = atom.mol_inversion_flag() {
        spec = spec.with_mol_inversion_flag(mol_inversion_flag);
    }
    if !atom.implicit_hydrogen() {
        spec = spec.with_implicit_hydrogen(false);
    }
    if !atom.tracked_isotopic_hydrogens().is_empty() {
        spec = spec.with_tracked_isotopic_hydrogens(atom.tracked_isotopic_hydrogens().to_vec());
    }
    if let Some(isotope) = atom.isotope() {
        spec = spec.with_isotope(isotope);
    }
    if let Some(atom_map) = atom.atom_map() {
        spec = spec.with_atom_map(atom_map);
    }
    if let Some(query) = atom.query().cloned() {
        spec = spec.with_query(query);
    }
    for (key, value) in atom.props() {
        spec = spec.with_prop(key.clone(), value.clone());
    }
    if let Some(info) = atom.pdb_residue_info().cloned() {
        spec = spec.with_pdb_residue_info(info);
    }
    spec
}

fn clear_atom_chemical_props(atom: &mut Atom) {
    // BEGIN RDKIT CPP FUNCTION ClearAtomChemicalProps
    // RDKit✔️✔️: void ClearAtomChemicalProps(RDKit::Atom *atom) {
    // RDKit✔️✔️:   TEST_ASSERT(atom);
    // RDKit✔️✔️:   atom->setIsotope(0);
    // RDKit✔️✔️:   atom->setFormalCharge(0);
    // RDKit✔️✔️:   atom->setNumExplicitHs(0);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ClearAtomChemicalProps
    atom.set_isotope(None);
    atom.set_formal_charge(0);
    atom.set_explicit_hydrogens(0);
}

fn report_parse_error(message: &str, throw_it: bool) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION ReportParseError
    // RDKit✔️✔️: void ReportParseError(const char *message, bool throwIt) {
    // RDKit✔️✔️:   PRECONDITION(message, "bad message");
    // RDKit✔️✔️:   if (!throwIt) {
    // RDKit✔️✔️:     BOOST_LOG(rdErrorLog) << "SMILES Parse Error: " << message << std::endl;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     throw SmilesParseException(message);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ReportParseError
    if throw_it {
        Err(SmilesParseError::ParseError(message.to_string()))
    } else {
        eprintln!("SMILES Parse Error: {message}");
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SmilesParseParams {
    sanitize: bool,
    allow_cxsmiles: bool,
    strict_cxsmiles: bool,
    parse_name: bool,
    remove_hs: bool,
    skip_cleanup: bool,
    debug_parse: bool,
    replacements: BTreeMap<String, String>,
}

impl Default for SmilesParseParams {
    fn default() -> Self {
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

impl SmilesParseParams {
    pub(crate) fn with_sanitize(sanitize: bool) -> Self {
        Self {
            sanitize,
            // RDKit's Python-facing `MolFromSmiles(..., sanitize=False)`
            // keeps explicit hydrogens and does not leave `_StereochemDone`
            // on the returned molecule. The explicit `removeHs=true` path is
            // still exercised separately through `mol_from_smiles()` params.
            remove_hs: sanitize,
            ..Default::default()
        }
    }

    #[cfg(test)]
    fn without_cxsmiles_for_test() -> Self {
        Self {
            allow_cxsmiles: false,
            ..Default::default()
        }
    }

    #[cfg(test)]
    fn without_parse_name_for_test() -> Self {
        Self {
            parse_name: false,
            ..Default::default()
        }
    }

    #[cfg(test)]
    fn non_strict_cxsmiles_for_test() -> Self {
        Self {
            strict_cxsmiles: false,
            ..Default::default()
        }
    }

    #[cfg(test)]
    fn with_replacements_for_test(replacements: BTreeMap<String, String>) -> Self {
        Self {
            replacements,
            ..Default::default()
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
struct PreprocessedSmiles {
    smiles: String,
    name: String,
    cx_part: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SmilesLexer<'a> {
    input: &'a str,
    scan_start: usize,
    scan_end: usize,
    scan_cursor: usize,
    current_token_position: usize,
    start_token: Option<SmilesToken>,
    in_atom_state: bool,
}

impl<'a> SmilesLexer<'a> {
    fn new(input: &'a str) -> Self {
        // BEGIN RDKIT CPP FUNCTION setup_smiles_string
        // RDKit✔️✔️: size_t setup_smiles_string(const std::string &text,yyscan_t yyscanner){
        // RDKit✔️✔️: //  YY_BUFFER_STATE buff=yysmiles__scan_string(text.c_str()+pos,yyscanner);
        // RDKit✔️✔️:   // Faster implementation of yysmiles__scan_string that handles trimming
        // RDKit✔️✔️:   YY_BUFFER_STATE b;
        // RDKit✔️✔️:   char *buf;
        // RDKit✔️✔️:   yyconst char * yybytes = text.c_str();
        // RDKit✔️✔️:   yy_size_t _yybytes_len=text.size(), n, start, end;
        // RDKit✔️✔️:   /* Get memory for full buffer, including space for trailing EOB's. */
        // RDKit✔️✔️:   n = _yybytes_len + 2;
        // RDKit✔️✔️:   buf = (char *) yysmiles_alloc(n ,yyscanner );
        // RDKit✔️✔️:   if ( ! buf ) {
        // RDKit✔️✔️:     smiles_lexer_error( "out of dynamic memory in yysmiles__scan_bytes()" );
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // ltrim
        // RDKit✔️✔️:
        // RDKit✔️✔️:   for(start = 0 ; start < _yybytes_len; ++start) {
        // RDKit✔️✔️:     if (yybytes[start] > 32) { break; }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   for(end = _yybytes_len ; end > start; --end) {
        // RDKit✔️✔️:     if (yybytes[end] > 32) { break; }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   _yybytes_len = end-start+1;
        // RDKit✔️✔️:   n = _yybytes_len + 2;
        // RDKit✔️✔️:   memcpy(buf, yybytes+start, _yybytes_len);
        // RDKit✔️✔️:
        // RDKit✔️✔️:
        // RDKit✔️✔️:   buf[_yybytes_len] = buf[_yybytes_len+1] = YY_END_OF_BUFFER_CHAR;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   b = yysmiles__scan_buffer(buf,n ,yyscanner);
        // RDKit✔️✔️:   if ( ! b ) {
        // RDKit✔️✔️:     smiles_lexer_error( "bad buffer in yysmiles__scan_bytes()" );
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   /* It's okay to grow etc. this buffer, and we should throw it
        // RDKit✔️✔️:    * away when we're done.
        // RDKit✔️✔️:    */
        // RDKit✔️✔️:   b->yy_is_our_buffer = 1;
        // RDKit✔️✔️:
        // RDKit✔️✔️:
        // RDKit✔️✔️:   POSTCONDITION(b,"invalid buffer");
        // RDKit✔️✔️:   return start;
        // RDKit✔️✔️:
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION setup_smiles_string
        let (scan_start, scan_end) = setup_smiles_scan_range(input);
        Self {
            input,
            scan_start,
            scan_end,
            scan_cursor: scan_start,
            current_token_position: 0,
            start_token: None,
            in_atom_state: false,
        }
    }

    #[cfg(test)]
    fn with_start_token_for_test(input: &'a str, start_token: SmilesToken) -> Self {
        let mut lexer = Self::new(input);
        lexer.start_token = Some(start_token);
        lexer
    }

    #[allow(dead_code)]
    fn scan_input(&self) -> &'a str {
        &self.input[self.scan_start..self.scan_end]
    }

    #[allow(dead_code)]
    fn next_token(&mut self) -> Result<SmilesToken, SmilesParseError> {
        // BEGIN RDKIT CPP LEXER RULES smiles.ll token rules
        // RDKit✔️✔️: @{
        // RDKit✔️✔️:   if (start_token)
        // RDKit✔️✔️:     {
        // RDKit✔️✔️:       int t = start_token;
        // RDKit✔️✔️:       start_token = 0;
        // RDKit✔️✔️:       return t;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️: @}
        // RDKit✔️✔️:
        // RDKit✔️✔️: @[' ']*TH { yylval->chiraltype = Atom::ChiralType::CHI_TETRAHEDRAL; return CHI_CLASS_TOKEN; }
        // RDKit✔️✔️: @[' ']*AL { yylval->chiraltype = Atom::ChiralType::CHI_ALLENE; return CHI_CLASS_TOKEN; }
        // RDKit✔️✔️: @[' ']*SP { yylval->chiraltype = Atom::ChiralType::CHI_SQUAREPLANAR; return CHI_CLASS_TOKEN; }
        // RDKit✔️✔️: @[' ']*TB { yylval->chiraltype = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL; return CHI_CLASS_TOKEN; }
        // RDKit✔️✔️: @[' ']*OH { yylval->chiraltype = Atom::ChiralType::CHI_OCTAHEDRAL; return CHI_CLASS_TOKEN; }
        // RDKit✔️✔️:
        // RDKit✔️✔️: @		{ return AT_TOKEN; }
        // RDKit✔️✔️:
        // RDKit✔️✔️: <IN_ATOM_STATE>He	{ yylval->atom = new Atom(2); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Li	{ yylval->atom = new Atom(3); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Be	{ yylval->atom = new Atom(4); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Ne	{ yylval->atom = new Atom(10); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Na	{ yylval->atom = new Atom(11); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Mg	{ yylval->atom = new Atom(12); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Al	{ yylval->atom = new Atom(13); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Si	{ yylval->atom = new Atom(14); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Ar	{ yylval->atom = new Atom(18); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>K	{ yylval->atom = new Atom(19); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Ca	{ yylval->atom = new Atom(20); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Sc	{ yylval->atom = new Atom(21); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Ti	{ yylval->atom = new Atom(22); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>V	{ yylval->atom = new Atom(23); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Cr	{ yylval->atom = new Atom(24); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Mn	{ yylval->atom = new Atom(25); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Fe	{ yylval->atom = new Atom(26); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Co	{ yylval->atom = new Atom(27); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Ni	{ yylval->atom = new Atom(28); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Cu	{ yylval->atom = new Atom(29); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Zn	{ yylval->atom = new Atom(30); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Ga	{ yylval->atom = new Atom(31); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Ge	{ yylval->atom = new Atom(32); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>As	{ yylval->atom = new Atom(33); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Se	{ yylval->atom = new Atom(34); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Kr	{ yylval->atom = new Atom(36); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Rb	{ yylval->atom = new Atom(37); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Sr	{ yylval->atom = new Atom(38); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Y	{ yylval->atom = new Atom(39); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Zr	{ yylval->atom = new Atom(40); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Nb	{ yylval->atom = new Atom(41); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Mo	{ yylval->atom = new Atom(42); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Tc	{ yylval->atom = new Atom(43); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Ru	{ yylval->atom = new Atom(44); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Rh	{ yylval->atom = new Atom(45); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Pd	{ yylval->atom = new Atom(46); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Ag	{ yylval->atom = new Atom(47); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Cd	{ yylval->atom = new Atom(48); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>In	{ yylval->atom = new Atom(49); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Sn	{ yylval->atom = new Atom(50); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Sb	{ yylval->atom = new Atom(51); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Te	{ yylval->atom = new Atom(52); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Xe	{ yylval->atom = new Atom(54); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Cs	{ yylval->atom = new Atom(55); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Ba	{ yylval->atom = new Atom(56); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>La	{ yylval->atom = new Atom(57); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Ce	{ yylval->atom = new Atom(58); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Pr	{ yylval->atom = new Atom(59); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Nd	{ yylval->atom = new Atom(60); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Pm	{ yylval->atom = new Atom(61); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Sm	{ yylval->atom = new Atom(62); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Eu	{ yylval->atom = new Atom(63); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Gd	{ yylval->atom = new Atom(64); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Tb	{ yylval->atom = new Atom(65); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Dy	{ yylval->atom = new Atom(66); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Ho	{ yylval->atom = new Atom(67); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Er	{ yylval->atom = new Atom(68); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Tm	{ yylval->atom = new Atom(69); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Yb	{ yylval->atom = new Atom(70); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Lu	{ yylval->atom = new Atom(71); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Hf	{ yylval->atom = new Atom(72); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Ta	{ yylval->atom = new Atom(73); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>W	{ yylval->atom = new Atom(74); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Re	{ yylval->atom = new Atom(75); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Os	{ yylval->atom = new Atom(76); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Ir	{ yylval->atom = new Atom(77); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Pt	{ yylval->atom = new Atom(78); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Au	{ yylval->atom = new Atom(79); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Hg	{ yylval->atom = new Atom(80); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Tl	{ yylval->atom = new Atom(81); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Pb	{ yylval->atom = new Atom(82); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Bi	{ yylval->atom = new Atom(83); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Po	{ yylval->atom = new Atom(84); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>At	{ yylval->atom = new Atom(85); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Rn	{ yylval->atom = new Atom(86); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Fr	{ yylval->atom = new Atom(87); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Ra	{ yylval->atom = new Atom(88); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Ac	{ yylval->atom = new Atom(89); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Th	{ yylval->atom = new Atom(90); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Pa	{ yylval->atom = new Atom(91); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>U	{ yylval->atom = new Atom(92); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Np	{ yylval->atom = new Atom(93); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Pu	{ yylval->atom = new Atom(94); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Am	{ yylval->atom = new Atom(95); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Cm	{ yylval->atom = new Atom(96); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Bk	{ yylval->atom = new Atom(97); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Cf	{ yylval->atom = new Atom(98); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Es	{ yylval->atom = new Atom(99); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Fm	{ yylval->atom = new Atom(100); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Md	{ yylval->atom = new Atom(101); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>No	{ yylval->atom = new Atom(102); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Lr	{ yylval->atom = new Atom(103); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Rf	{ yylval->atom = new Atom(104); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Db	{ yylval->atom = new Atom(105); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Sg	{ yylval->atom = new Atom(106); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Bh	{ yylval->atom = new Atom(107); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Hs	{ yylval->atom = new Atom(108); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Mt	{ yylval->atom = new Atom(109); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Ds	{ yylval->atom = new Atom(110); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Rg	{ yylval->atom = new Atom(111); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Cn	{ yylval->atom = new Atom(112); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Nh	{ yylval->atom = new Atom(113); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Fl	{ yylval->atom = new Atom(114); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Mc	{ yylval->atom = new Atom(115); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Lv	{ yylval->atom = new Atom(116); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Ts	{ yylval->atom = new Atom(117); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Og	{ yylval->atom = new Atom(118); return ATOM_TOKEN; }
        // RDKit✔️✔️:
        // RDKit✔️✔️: <IN_ATOM_STATE>Uun	{ yylval->atom = new Atom(110); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Uuu	{ yylval->atom = new Atom(111); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Uub	{ yylval->atom = new Atom(112); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Uut	{ yylval->atom = new Atom(113); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Uuq	{ yylval->atom = new Atom(114); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Uup	{ yylval->atom = new Atom(115); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Uuh	{ yylval->atom = new Atom(116); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Uus	{ yylval->atom = new Atom(117); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>Uuo	{ yylval->atom = new Atom(118); return ATOM_TOKEN; }
        // RDKit✔️✔️:
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
        // RDKit✔️✔️:
        // RDKit✔️✔️: H			{
        // RDKit✔️✔️: 				return H_TOKEN;
        // RDKit✔️✔️: 			}
        // RDKit✔️✔️:
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
        // RDKit✔️✔️:
        // RDKit✔️✔️: <IN_ATOM_STATE>si   {	yylval->atom = new Atom( 14 );
        // RDKit✔️✔️: 			yylval->atom->setIsAromatic(true);
        // RDKit✔️✔️: 				return AROMATIC_ATOM_TOKEN;
        // RDKit✔️✔️: 			}
        // RDKit✔️✔️: <IN_ATOM_STATE>as   {	yylval->atom = new Atom( 33 );
        // RDKit✔️✔️: 			yylval->atom->setIsAromatic(true);
        // RDKit✔️✔️: 				return AROMATIC_ATOM_TOKEN;
        // RDKit✔️✔️: 			}
        // RDKit✔️✔️: <IN_ATOM_STATE>se   {	yylval->atom = new Atom( 34 );
        // RDKit✔️✔️: 			yylval->atom->setIsAromatic(true);
        // RDKit✔️✔️: 				return AROMATIC_ATOM_TOKEN;
        // RDKit✔️✔️: 			}
        // RDKit✔️✔️: <IN_ATOM_STATE>te   {	yylval->atom = new Atom( 52 );
        // RDKit✔️✔️: 			yylval->atom->setIsAromatic(true);
        // RDKit✔️✔️: 				return AROMATIC_ATOM_TOKEN;
        // RDKit✔️✔️: 			}
        // RDKit✔️✔️:
        // RDKit✔️✔️: \* 	            {   yylval->atom = new Atom( 0 );
        // RDKit✔️✔️: 		            yylval->atom->setProp(common_properties::dummyLabel,
        // RDKit✔️✔️:                                                         std::string("*"));
        // RDKit✔️✔️:                                 // must be ORGANIC_ATOM_TOKEN because
        // RDKit✔️✔️:                                 // we aren't in square brackets:
        // RDKit✔️✔️: 				return ORGANIC_ATOM_TOKEN;
        // RDKit✔️✔️: 			}
        // RDKit✔️✔️:
        // RDKit✔️✔️: <IN_ATOM_STATE>\: 	{ return COLON_TOKEN; }
        // RDKit✔️✔️:
        // RDKit✔️✔️: <IN_ATOM_STATE>\# 	{ return HASH_TOKEN; }
        // RDKit✔️✔️:
        // RDKit✔️✔️: %{ // Biovia quoted heavy atom workaround %}
        // RDKit✔️✔️: <IN_ATOM_STATE>\'Rf\'	{ yylval->atom = new Atom(104); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>\'Db\'	{ yylval->atom = new Atom(105); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>\'Sg\'	{ yylval->atom = new Atom(106); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>\'Bh\'	{ yylval->atom = new Atom(107); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>\'Hs\'	{ yylval->atom = new Atom(108); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>\'Mt\'	{ yylval->atom = new Atom(109); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>\'Ds\'	{ yylval->atom = new Atom(110); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>\'Rg\'	{ yylval->atom = new Atom(111); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>\'Cn\'	{ yylval->atom = new Atom(112); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>\'Nh\'	{ yylval->atom = new Atom(113); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>\'Fl\'	{ yylval->atom = new Atom(114); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>\'Mc\'	{ yylval->atom = new Atom(115); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>\'Lv\'	{ yylval->atom = new Atom(116); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>\'Ts\'	{ yylval->atom = new Atom(117); return ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>\'Og\'	{ yylval->atom = new Atom(118); return ATOM_TOKEN; }
        // RDKit✔️✔️:
        // RDKit✔️✔️: \=	{ yylval->bond = new Bond(Bond::DOUBLE); return BOND_TOKEN; }
        // RDKit✔️✔️: \#	{ yylval->bond = new Bond(Bond::TRIPLE); return BOND_TOKEN; }
        // RDKit✔️✔️: \:	{ yylval->bond = new Bond(Bond::AROMATIC);
        // RDKit✔️✔️: 	  yylval->bond->setIsAromatic(true); return BOND_TOKEN; }
        // RDKit✔️✔️: \$	{ yylval->bond = new Bond(Bond::QUADRUPLE); return BOND_TOKEN; }
        // RDKit✔️✔️: \-\>	{ yylval->bond = new Bond(Bond::DATIVER); return BOND_TOKEN; }
        // RDKit✔️✔️: \<\-	{ yylval->bond = new Bond(Bond::DATIVEL); return BOND_TOKEN; }
        // RDKit✔️✔️: \~	{ yylval->bond = new QueryBond();
        // RDKit✔️✔️: 	  yylval->bond->setQuery(makeBondNullQuery()); return BOND_TOKEN;  }
        // RDKit✔️✔️:
        // RDKit✔️✔️: [\\]{1,2}    { yylval->bond = new Bond(Bond::UNSPECIFIED);
        // RDKit✔️✔️: 	yylval->bond->setProp(RDKit::common_properties::_unspecifiedOrder,1);
        // RDKit✔️✔️: 	yylval->bond->setBondDir(Bond::ENDDOWNRIGHT); return BOND_TOKEN;  }
        // RDKit✔️✔️:
        // RDKit✔️✔️: [\/]    { yylval->bond = new Bond(Bond::UNSPECIFIED);
        // RDKit✔️✔️: 	yylval->bond->setProp(RDKit::common_properties::_unspecifiedOrder,1);
        // RDKit✔️✔️: 	yylval->bond->setBondDir(Bond::ENDUPRIGHT); return BOND_TOKEN;  }
        // RDKit✔️✔️:
        // RDKit✔️✔️: \-			{ return MINUS_TOKEN; }
        // RDKit✔️✔️: \+			{ return PLUS_TOKEN; }
        // RDKit✔️✔️: \(       	{ return GROUP_OPEN_TOKEN; }
        // RDKit✔️✔️: \)       	{ return GROUP_CLOSE_TOKEN; }
        // RDKit✔️✔️: \[			{ BEGIN IN_ATOM_STATE; return ATOM_OPEN_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>\]	{ BEGIN INITIAL; return ATOM_CLOSE_TOKEN; }
        // RDKit✔️✔️: \.       	{ return SEPARATOR_TOKEN; }
        // RDKit✔️✔️: \%              { return PERCENT_TOKEN; }
        // RDKit✔️✔️: [0]		{ yylval->ival = 0; return ZERO_TOKEN; }
        // RDKit✔️✔️: [1-9]		{ yylval->ival = yytext[0] - '0'; return NONZERO_DIGIT_TOKEN; }
        // RDKit✔️✔️: \n		return 0;
        // RDKit✔️✔️: <<EOF>>		{ return EOS_TOKEN; }
        // RDKit✔️✔️: .		return BAD_CHARACTER;
        // END RDKIT CPP LEXER RULES smiles.ll token rules
        if let Some(token) = self.start_token.take() {
            return Ok(token);
        }

        if self.scan_cursor >= self.scan_end {
            return Ok(SmilesToken::Eos);
        }

        let remaining = &self.input[self.scan_cursor..self.scan_end];
        let (len, token) = self.match_next_token(remaining);
        self.scan_cursor += len;
        self.current_token_position += len;
        if token == SmilesToken::Eos {
            self.scan_cursor = self.scan_end;
        }
        Ok(token)
    }

    fn match_next_token(&mut self, remaining: &str) -> (usize, SmilesToken) {
        if self.in_atom_state {
            if let Some((symbol, atomic_number)) = match_bracket_atom_symbol(remaining) {
                return (
                    symbol.len(),
                    SmilesToken::Atom(SmilesAtomToken::new(atomic_number)),
                );
            }
            if let Some((symbol, atomic_number)) = match_quoted_biovia_atom_symbol(remaining) {
                return (
                    symbol.len(),
                    SmilesToken::Atom(SmilesAtomToken::new(atomic_number)),
                );
            }
            if let Some((symbol, atomic_number)) = match_bracket_aromatic_atom_symbol(remaining) {
                return (
                    symbol.len(),
                    SmilesToken::AromaticAtom(SmilesAtomToken::aromatic(atomic_number)),
                );
            }
            if remaining.starts_with(']') {
                self.in_atom_state = false;
                return (1, SmilesToken::AtomClose);
            }
        }

        if let Some((len, chiral_tag)) = match_chiral_class(remaining) {
            return (len, SmilesToken::ChiralClass(chiral_tag));
        }
        if remaining.starts_with('@') {
            return (1, SmilesToken::At);
        }
        if let Some((symbol, atomic_number)) = match_organic_atom_symbol(remaining) {
            return (
                symbol.len(),
                SmilesToken::OrganicAtom(SmilesAtomToken::new(atomic_number)),
            );
        }
        if remaining.starts_with('H') {
            return (1, SmilesToken::H);
        }
        if let Some((symbol, atomic_number)) = match_aromatic_atom_symbol(remaining) {
            return (
                symbol.len(),
                SmilesToken::AromaticAtom(SmilesAtomToken::aromatic(atomic_number)),
            );
        }
        if remaining.starts_with('*') {
            return (1, SmilesToken::OrganicAtom(SmilesAtomToken::dummy()));
        }
        if remaining.starts_with("->") {
            return (
                2,
                SmilesToken::Bond(SmilesBondToken::new(BondOrder::DativeRight)),
            );
        }
        if remaining.starts_with("<-") {
            return (
                2,
                SmilesToken::Bond(SmilesBondToken::new(BondOrder::DativeLeft)),
            );
        }
        if remaining.starts_with('\\') {
            let len = if remaining.starts_with("\\\\") { 2 } else { 1 };
            return (
                len,
                SmilesToken::Bond(SmilesBondToken::directional(BondDirection::EndDownRight)),
            );
        }
        if remaining.starts_with('/') {
            return (
                1,
                SmilesToken::Bond(SmilesBondToken::directional(BondDirection::EndUpRight)),
            );
        }

        let ch = remaining
            .chars()
            .next()
            .expect("remaining input is non-empty");
        match ch {
            '=' => (
                1,
                SmilesToken::Bond(SmilesBondToken::new(BondOrder::Double)),
            ),
            '#' if self.in_atom_state => (1, SmilesToken::Hash),
            '#' => (
                1,
                SmilesToken::Bond(SmilesBondToken::new(BondOrder::Triple)),
            ),
            ':' if self.in_atom_state => (1, SmilesToken::Colon),
            ':' => (1, SmilesToken::Bond(SmilesBondToken::aromatic())),
            '$' => (
                1,
                SmilesToken::Bond(SmilesBondToken::new(BondOrder::Quadruple)),
            ),
            '~' => (1, SmilesToken::Bond(SmilesBondToken::null_query())),
            '-' => (1, SmilesToken::Minus),
            '+' => (1, SmilesToken::Plus),
            '(' => (1, SmilesToken::GroupOpen),
            ')' => (1, SmilesToken::GroupClose),
            '[' => {
                self.in_atom_state = true;
                (1, SmilesToken::AtomOpen)
            }
            '.' => (1, SmilesToken::Separator),
            '%' => (1, SmilesToken::Percent),
            '0' => (1, SmilesToken::Zero),
            '1'..='9' => (1, SmilesToken::NonzeroDigit(i32::from(ch as u8 - b'0'))),
            '\n' => (1, SmilesToken::Eos),
            _ => (ch.len_utf8(), SmilesToken::BadCharacter(ch)),
        }
    }
}

fn match_chiral_class(input: &str) -> Option<(usize, ChiralTag)> {
    let rest = input.strip_prefix('@')?;
    let spaces = rest.bytes().take_while(|byte| *byte == b' ').count();
    let label = &rest[spaces..];
    let tag = if label.starts_with("TH") {
        ChiralTag::Tetrahedral
    } else if label.starts_with("AL") {
        ChiralTag::Allene
    } else if label.starts_with("SP") {
        ChiralTag::SquarePlanar
    } else if label.starts_with("TB") {
        ChiralTag::TrigonalBipyramidal
    } else if label.starts_with("OH") {
        ChiralTag::Octahedral
    } else {
        return None;
    };
    Some((1 + spaces + 2, tag))
}

fn match_organic_atom_symbol(input: &str) -> Option<(&'static str, u8)> {
    for (symbol, atomic_number) in [("Cl", 17), ("Br", 35)] {
        if input.starts_with(symbol) {
            return Some((symbol, atomic_number));
        }
    }
    match input.as_bytes().first().copied()? {
        b'B' => Some(("B", 5)),
        b'C' => Some(("C", 6)),
        b'N' => Some(("N", 7)),
        b'O' => Some(("O", 8)),
        b'P' => Some(("P", 15)),
        b'S' => Some(("S", 16)),
        b'F' => Some(("F", 9)),
        b'I' => Some(("I", 53)),
        _ => None,
    }
}

fn match_aromatic_atom_symbol(input: &str) -> Option<(&'static str, u8)> {
    match input.as_bytes().first().copied()? {
        b'b' => Some(("b", 5)),
        b'c' => Some(("c", 6)),
        b'n' => Some(("n", 7)),
        b'o' => Some(("o", 8)),
        b'p' => Some(("p", 15)),
        b's' => Some(("s", 16)),
        _ => None,
    }
}

fn match_bracket_aromatic_atom_symbol(input: &str) -> Option<(&'static str, u8)> {
    for (symbol, atomic_number) in [("si", 14), ("as", 33), ("se", 34), ("te", 52)] {
        if input.starts_with(symbol) {
            return Some((symbol, atomic_number));
        }
    }
    match_aromatic_atom_symbol(input)
}

fn allow_nontetrahedral_chirality() -> bool {
    true
}

fn match_bracket_atom_symbol(input: &str) -> Option<(&'static str, u8)> {
    match input.as_bytes() {
        [b'U', b'u', b'n', ..] => Some(("Uun", 110)),
        [b'U', b'u', b'u', ..] => Some(("Uuu", 111)),
        [b'U', b'u', b'b', ..] => Some(("Uub", 112)),
        [b'U', b'u', b't', ..] => Some(("Uut", 113)),
        [b'U', b'u', b'q', ..] => Some(("Uuq", 114)),
        [b'U', b'u', b'p', ..] => Some(("Uup", 115)),
        [b'U', b'u', b'h', ..] => Some(("Uuh", 116)),
        [b'U', b'u', b's', ..] => Some(("Uus", 117)),
        [b'U', b'u', b'o', ..] => Some(("Uuo", 118)),
        [b'H', b'e', ..] => Some(("He", 2)),
        [b'L', b'i', ..] => Some(("Li", 3)),
        [b'B', b'e', ..] => Some(("Be", 4)),
        [b'N', b'e', ..] => Some(("Ne", 10)),
        [b'N', b'a', ..] => Some(("Na", 11)),
        [b'M', b'g', ..] => Some(("Mg", 12)),
        [b'A', b'l', ..] => Some(("Al", 13)),
        [b'S', b'i', ..] => Some(("Si", 14)),
        [b'A', b'r', ..] => Some(("Ar", 18)),
        [b'C', b'a', ..] => Some(("Ca", 20)),
        [b'S', b'c', ..] => Some(("Sc", 21)),
        [b'T', b'i', ..] => Some(("Ti", 22)),
        [b'C', b'r', ..] => Some(("Cr", 24)),
        [b'M', b'n', ..] => Some(("Mn", 25)),
        [b'F', b'e', ..] => Some(("Fe", 26)),
        [b'C', b'o', ..] => Some(("Co", 27)),
        [b'N', b'i', ..] => Some(("Ni", 28)),
        [b'C', b'u', ..] => Some(("Cu", 29)),
        [b'Z', b'n', ..] => Some(("Zn", 30)),
        [b'G', b'a', ..] => Some(("Ga", 31)),
        [b'G', b'e', ..] => Some(("Ge", 32)),
        [b'A', b's', ..] => Some(("As", 33)),
        [b'S', b'e', ..] => Some(("Se", 34)),
        [b'K', b'r', ..] => Some(("Kr", 36)),
        [b'R', b'b', ..] => Some(("Rb", 37)),
        [b'S', b'r', ..] => Some(("Sr", 38)),
        [b'Z', b'r', ..] => Some(("Zr", 40)),
        [b'N', b'b', ..] => Some(("Nb", 41)),
        [b'M', b'o', ..] => Some(("Mo", 42)),
        [b'T', b'c', ..] => Some(("Tc", 43)),
        [b'R', b'u', ..] => Some(("Ru", 44)),
        [b'R', b'h', ..] => Some(("Rh", 45)),
        [b'P', b'd', ..] => Some(("Pd", 46)),
        [b'A', b'g', ..] => Some(("Ag", 47)),
        [b'C', b'd', ..] => Some(("Cd", 48)),
        [b'I', b'n', ..] => Some(("In", 49)),
        [b'S', b'n', ..] => Some(("Sn", 50)),
        [b'S', b'b', ..] => Some(("Sb", 51)),
        [b'T', b'e', ..] => Some(("Te", 52)),
        [b'X', b'e', ..] => Some(("Xe", 54)),
        [b'C', b's', ..] => Some(("Cs", 55)),
        [b'B', b'a', ..] => Some(("Ba", 56)),
        [b'L', b'a', ..] => Some(("La", 57)),
        [b'C', b'e', ..] => Some(("Ce", 58)),
        [b'P', b'r', ..] => Some(("Pr", 59)),
        [b'N', b'd', ..] => Some(("Nd", 60)),
        [b'P', b'm', ..] => Some(("Pm", 61)),
        [b'S', b'm', ..] => Some(("Sm", 62)),
        [b'E', b'u', ..] => Some(("Eu", 63)),
        [b'G', b'd', ..] => Some(("Gd", 64)),
        [b'T', b'b', ..] => Some(("Tb", 65)),
        [b'D', b'y', ..] => Some(("Dy", 66)),
        [b'H', b'o', ..] => Some(("Ho", 67)),
        [b'E', b'r', ..] => Some(("Er", 68)),
        [b'T', b'm', ..] => Some(("Tm", 69)),
        [b'Y', b'b', ..] => Some(("Yb", 70)),
        [b'L', b'u', ..] => Some(("Lu", 71)),
        [b'H', b'f', ..] => Some(("Hf", 72)),
        [b'T', b'a', ..] => Some(("Ta", 73)),
        [b'R', b'e', ..] => Some(("Re", 75)),
        [b'O', b's', ..] => Some(("Os", 76)),
        [b'I', b'r', ..] => Some(("Ir", 77)),
        [b'P', b't', ..] => Some(("Pt", 78)),
        [b'A', b'u', ..] => Some(("Au", 79)),
        [b'H', b'g', ..] => Some(("Hg", 80)),
        [b'T', b'l', ..] => Some(("Tl", 81)),
        [b'P', b'b', ..] => Some(("Pb", 82)),
        [b'B', b'i', ..] => Some(("Bi", 83)),
        [b'P', b'o', ..] => Some(("Po", 84)),
        [b'A', b't', ..] => Some(("At", 85)),
        [b'R', b'n', ..] => Some(("Rn", 86)),
        [b'F', b'r', ..] => Some(("Fr", 87)),
        [b'R', b'a', ..] => Some(("Ra", 88)),
        [b'A', b'c', ..] => Some(("Ac", 89)),
        [b'T', b'h', ..] => Some(("Th", 90)),
        [b'P', b'a', ..] => Some(("Pa", 91)),
        [b'N', b'p', ..] => Some(("Np", 93)),
        [b'P', b'u', ..] => Some(("Pu", 94)),
        [b'A', b'm', ..] => Some(("Am", 95)),
        [b'C', b'm', ..] => Some(("Cm", 96)),
        [b'B', b'k', ..] => Some(("Bk", 97)),
        [b'C', b'f', ..] => Some(("Cf", 98)),
        [b'E', b's', ..] => Some(("Es", 99)),
        [b'F', b'm', ..] => Some(("Fm", 100)),
        [b'M', b'd', ..] => Some(("Md", 101)),
        [b'N', b'o', ..] => Some(("No", 102)),
        [b'L', b'r', ..] => Some(("Lr", 103)),
        [b'R', b'f', ..] => Some(("Rf", 104)),
        [b'D', b'b', ..] => Some(("Db", 105)),
        [b'S', b'g', ..] => Some(("Sg", 106)),
        [b'B', b'h', ..] => Some(("Bh", 107)),
        [b'H', b's', ..] => Some(("Hs", 108)),
        [b'M', b't', ..] => Some(("Mt", 109)),
        [b'D', b's', ..] => Some(("Ds", 110)),
        [b'R', b'g', ..] => Some(("Rg", 111)),
        [b'C', b'n', ..] => Some(("Cn", 112)),
        [b'N', b'h', ..] => Some(("Nh", 113)),
        [b'F', b'l', ..] => Some(("Fl", 114)),
        [b'M', b'c', ..] => Some(("Mc", 115)),
        [b'L', b'v', ..] => Some(("Lv", 116)),
        [b'T', b's', ..] => Some(("Ts", 117)),
        [b'O', b'g', ..] => Some(("Og", 118)),
        [b'K', ..] => Some(("K", 19)),
        [b'V', ..] => Some(("V", 23)),
        [b'Y', ..] => Some(("Y", 39)),
        [b'W', ..] => Some(("W", 74)),
        [b'U', ..] => Some(("U", 92)),
        _ => None,
    }
}

fn match_quoted_biovia_atom_symbol(input: &str) -> Option<(&'static str, u8)> {
    match input.as_bytes() {
        [b'\'', b'R', b'f', b'\'', ..] => Some(("'Rf'", 104)),
        [b'\'', b'D', b'b', b'\'', ..] => Some(("'Db'", 105)),
        [b'\'', b'S', b'g', b'\'', ..] => Some(("'Sg'", 106)),
        [b'\'', b'B', b'h', b'\'', ..] => Some(("'Bh'", 107)),
        [b'\'', b'H', b's', b'\'', ..] => Some(("'Hs'", 108)),
        [b'\'', b'M', b't', b'\'', ..] => Some(("'Mt'", 109)),
        [b'\'', b'D', b's', b'\'', ..] => Some(("'Ds'", 110)),
        [b'\'', b'R', b'g', b'\'', ..] => Some(("'Rg'", 111)),
        [b'\'', b'C', b'n', b'\'', ..] => Some(("'Cn'", 112)),
        [b'\'', b'N', b'h', b'\'', ..] => Some(("'Nh'", 113)),
        [b'\'', b'F', b'l', b'\'', ..] => Some(("'Fl'", 114)),
        [b'\'', b'M', b'c', b'\'', ..] => Some(("'Mc'", 115)),
        [b'\'', b'L', b'v', b'\'', ..] => Some(("'Lv'", 116)),
        [b'\'', b'T', b's', b'\'', ..] => Some(("'Ts'", 117)),
        [b'\'', b'O', b'g', b'\'', ..] => Some(("'Og'", 118)),
        _ => None,
    }
}

#[cfg(test)]
const BRACKET_ATOM_SYMBOLS: &[(&str, u8)] = &[
    ("Uun", 110),
    ("Uuu", 111),
    ("Uub", 112),
    ("Uut", 113),
    ("Uuq", 114),
    ("Uup", 115),
    ("Uuh", 116),
    ("Uus", 117),
    ("Uuo", 118),
    ("He", 2),
    ("Li", 3),
    ("Be", 4),
    ("Ne", 10),
    ("Na", 11),
    ("Mg", 12),
    ("Al", 13),
    ("Si", 14),
    ("Ar", 18),
    ("Ca", 20),
    ("Sc", 21),
    ("Ti", 22),
    ("Cr", 24),
    ("Mn", 25),
    ("Fe", 26),
    ("Co", 27),
    ("Ni", 28),
    ("Cu", 29),
    ("Zn", 30),
    ("Ga", 31),
    ("Ge", 32),
    ("As", 33),
    ("Se", 34),
    ("Kr", 36),
    ("Rb", 37),
    ("Sr", 38),
    ("Zr", 40),
    ("Nb", 41),
    ("Mo", 42),
    ("Tc", 43),
    ("Ru", 44),
    ("Rh", 45),
    ("Pd", 46),
    ("Ag", 47),
    ("Cd", 48),
    ("In", 49),
    ("Sn", 50),
    ("Sb", 51),
    ("Te", 52),
    ("Xe", 54),
    ("Cs", 55),
    ("Ba", 56),
    ("La", 57),
    ("Ce", 58),
    ("Pr", 59),
    ("Nd", 60),
    ("Pm", 61),
    ("Sm", 62),
    ("Eu", 63),
    ("Gd", 64),
    ("Tb", 65),
    ("Dy", 66),
    ("Ho", 67),
    ("Er", 68),
    ("Tm", 69),
    ("Yb", 70),
    ("Lu", 71),
    ("Hf", 72),
    ("Ta", 73),
    ("Re", 75),
    ("Os", 76),
    ("Ir", 77),
    ("Pt", 78),
    ("Au", 79),
    ("Hg", 80),
    ("Tl", 81),
    ("Pb", 82),
    ("Bi", 83),
    ("Po", 84),
    ("At", 85),
    ("Rn", 86),
    ("Fr", 87),
    ("Ra", 88),
    ("Ac", 89),
    ("Th", 90),
    ("Pa", 91),
    ("Np", 93),
    ("Pu", 94),
    ("Am", 95),
    ("Cm", 96),
    ("Bk", 97),
    ("Cf", 98),
    ("Es", 99),
    ("Fm", 100),
    ("Md", 101),
    ("No", 102),
    ("Lr", 103),
    ("Rf", 104),
    ("Db", 105),
    ("Sg", 106),
    ("Bh", 107),
    ("Hs", 108),
    ("Mt", 109),
    ("Ds", 110),
    ("Rg", 111),
    ("Cn", 112),
    ("Nh", 113),
    ("Fl", 114),
    ("Mc", 115),
    ("Lv", 116),
    ("Ts", 117),
    ("Og", 118),
    ("K", 19),
    ("V", 23),
    ("Y", 39),
    ("W", 74),
    ("U", 92),
];

#[cfg(test)]
const QUOTED_BIOVIA_ATOM_SYMBOLS: &[(&str, u8)] = &[
    ("'Rf'", 104),
    ("'Db'", 105),
    ("'Sg'", 106),
    ("'Bh'", 107),
    ("'Hs'", 108),
    ("'Mt'", 109),
    ("'Ds'", 110),
    ("'Rg'", 111),
    ("'Cn'", 112),
    ("'Nh'", 113),
    ("'Fl'", 114),
    ("'Mc'", 115),
    ("'Lv'", 116),
    ("'Ts'", 117),
    ("'Og'", 118),
];

fn setup_smiles_scan_range(text: &str) -> (usize, usize) {
    let bytes = text.as_bytes();
    let mut start = 0;
    while start < bytes.len() && bytes[start] <= 32 {
        start += 1;
    }

    let mut end = bytes.len();
    while end > start && bytes[end - 1] <= 32 {
        end -= 1;
    }

    (start, end)
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(dead_code)]
enum SmilesToken {
    StartMol,
    StartAtom,
    StartBond,
    AromaticAtom(SmilesAtomToken),
    Atom(SmilesAtomToken),
    OrganicAtom(SmilesAtomToken),
    NonzeroDigit(i32),
    Zero,
    GroupOpen,
    GroupClose,
    Separator,
    LoopConnector,
    Minus,
    Plus,
    H,
    At,
    Percent,
    Colon,
    Hash,
    Bond(SmilesBondToken),
    ChiralClass(ChiralTag),
    AtomOpen,
    AtomClose,
    BadCharacter(char),
    Eos,
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[allow(dead_code)]
struct SmilesAtomToken {
    spec: AtomSpec,
}

#[allow(dead_code)]
impl SmilesAtomToken {
    fn new(atomic_number: u8) -> Self {
        Self {
            spec: AtomSpec::new(Element::from_atomic_number(atomic_number).unwrap()),
        }
    }

    fn aromatic(atomic_number: u8) -> Self {
        Self {
            spec: AtomSpec::new(Element::from_atomic_number(atomic_number).unwrap())
                .with_aromatic(true),
        }
    }

    fn dummy() -> Self {
        Self {
            spec: AtomSpec::new(Element::DUMMY).with_prop("dummyLabel", "*"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
struct SmilesBondToken {
    order: BondOrder,
    is_aromatic: bool,
    direction: BondDirection,
    explicit_unspecified_order: bool,
    is_null_query: bool,
}

#[allow(dead_code)]
impl SmilesBondToken {
    fn new(order: BondOrder) -> Self {
        Self {
            order,
            is_aromatic: false,
            direction: BondDirection::None,
            explicit_unspecified_order: false,
            is_null_query: false,
        }
    }

    fn aromatic() -> Self {
        Self {
            order: BondOrder::Aromatic,
            is_aromatic: true,
            direction: BondDirection::None,
            explicit_unspecified_order: false,
            is_null_query: false,
        }
    }

    fn directional(direction: BondDirection) -> Self {
        Self {
            order: BondOrder::Unspecified,
            is_aromatic: false,
            direction,
            explicit_unspecified_order: true,
            is_null_query: false,
        }
    }

    fn null_query() -> Self {
        Self {
            order: BondOrder::Unspecified,
            is_aromatic: false,
            direction: BondDirection::None,
            explicit_unspecified_order: false,
            is_null_query: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SmilesParser<'a> {
    lexer: SmilesLexer<'a>,
    lookahead: Option<SmilesToken>,
}

impl<'a> SmilesParser<'a> {
    fn new(lexer: SmilesLexer<'a>) -> Self {
        Self {
            lexer,
            lookahead: None,
        }
    }

    fn parse_mol(&mut self, state: &mut SmilesBuildState) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP FUNCTION smiles_parse
        // RDKit✔️✔️: int smiles_parse(const std::string &inp, std::vector<RDKit::RWMol *> &molVect) {
        // RDKit✔️✔️:   auto start_tok = static_cast<int>(START_MOL);
        // RDKit✔️✔️:   Atom *atom = nullptr;
        // RDKit✔️✔️:   Bond *bond = nullptr;
        // RDKit✔️✔️:   return smiles_parse_helper(inp, molVect, atom, bond, start_tok);
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION smiles_parse

        // BEGIN RDKIT CPP FUNCTION smiles_parse_helper
        // RDKit✔️✔️: int smiles_parse_helper(const std::string &inp,
        // RDKit✔️✔️:                         std::vector<RDKit::RWMol *> &molVect, Atom *&atom,
        // RDKit✔️✔️:                         Bond *&bond, int start_tok) {
        // RDKit✔️✔️:     return generic_parse_helper<yysmiles_lex_init,
        // RDKit✔️✔️:                                 setup_smiles_string,
        // RDKit✔️✔️:                                 yysmiles_lex_destroy>(yysmiles_parse,
        // RDKit✔️✔️:                                                       inp,
        // RDKit✔️✔️:                                                       molVect,
        // RDKit✔️✔️:                                                       atom,
        // RDKit✔️✔️:                                                       bond,
        // RDKit✔️✔️:                                                       start_tok,
        // RDKit✔️✔️:                                                       "SMILES");
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION smiles_parse_helper

        // BEGIN RDKIT CPP FUNCTION generic_parse_helper
        // RDKit✔️✔️: template<int(*lex_init)(void**),
        // RDKit✔️✔️:          size_t(*string_setup)(const std::string &, void *),
        // RDKit✔️✔️:          int(*lex_destroy)(void*),
        // RDKit✔️✔️:          typename T>
        // RDKit✔️✔️: int generic_parse_helper(T parser,
        // RDKit✔️✔️:                          const std::string &inp,
        // RDKit✔️✔️:                          std::vector<RDKit::RWMol *> &molVect,
        // RDKit✔️✔️:                          Atom *&atom,
        // RDKit✔️✔️:                          Bond *&bond,
        // RDKit✔️✔️:                          int start_tok,
        // RDKit✔️✔️:                          const std::string& input_type) {
        // RDKit✔️✔️:   std::vector<std::pair<unsigned int, unsigned int>> branchPoints;
        // RDKit✔️✔️:   void *scanner;
        // RDKit✔️✔️:   int res = 1;  // initialize with fail code
        // RDKit✔️✔️:
        // RDKit✔️✔️:   TEST_ASSERT(!lex_init(&scanner));
        // RDKit✔️✔️:   size_t ltrim = 0;
        // RDKit✔️✔️:   try {
        // RDKit✔️✔️:     ltrim = string_setup(inp, scanner);
        // RDKit✔️✔️:     unsigned numAtomsParsed = 0;
        // RDKit✔️✔️:     unsigned numBondsParsed = 0;
        // RDKit✔️✔️:     // NOTE: This variable will be used to point to the location of the
        // RDKit✔️✔️:     // offending token if we encounter a syntax error
        // RDKit✔️✔️:     unsigned int current_token_position = 0;
        // RDKit✔️✔️:     res = parser(inp.c_str() + ltrim, &molVect, atom, bond,
        // RDKit✔️✔️:                          numAtomsParsed, numBondsParsed, branchPoints, scanner,
        // RDKit✔️✔️:                          start_tok, current_token_position);
        // RDKit✔️✔️:   } catch (...) {
        // RDKit✔️✔️:     lex_destroy(scanner);
        // RDKit✔️✔️:     throw;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   lex_destroy(scanner);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (!branchPoints.empty()) {
        // RDKit✔️✔️:     auto input = inp.c_str() + ltrim;
        // RDKit✔️✔️:     // If there are multiple unclosed brackets, we want to report them all at
        // RDKit✔️✔️:     // once.
        // RDKit✔️✔️:     for (auto [_, open_bracket_position] : branchPoints) {
        // RDKit✔️✔️:       SmilesParseOps::detail::printSyntaxErrorMessage(
        // RDKit✔️✔️:           input, "extra open parentheses", open_bracket_position, input_type);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (res == 1 || !branchPoints.empty()) {
        // RDKit✔️✔️:     throw SmilesParseException("Failed parsing " + input_type + " '" + inp + "'");
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION generic_parse_helper
        let first_atom = self.parse_simple_atomd()?;
        state.add_first_atom(first_atom)?;

        loop {
            match self.next_token()? {
                SmilesToken::Eos => {
                    if state.branch_stack.is_empty() {
                        return Ok(());
                    }
                    return Err(SmilesParseError::ParseError(
                        "extra open parentheses".to_string(),
                    ));
                }
                SmilesToken::OrganicAtom(atom)
                | SmilesToken::AromaticAtom(atom)
                | SmilesToken::Atom(atom) => {
                    state.add_atom_connected_to_active(atom)?;
                }
                SmilesToken::AtomOpen => {
                    let atom = self.parse_bracket_atomd()?;
                    state.add_atom_connected_to_active(atom)?;
                }
                SmilesToken::Bond(bond) => {
                    if self.next_token_starts_ring_number()? {
                        let ring_number = self.parse_ring_number()?;
                        state.add_explicit_bond_ring_marker(bond, ring_number)?;
                    } else {
                        let atom = self.parse_simple_atomd()?;
                        state.add_explicit_bond_to_atom(bond, atom)?;
                    }
                }
                SmilesToken::Minus => {
                    if self.next_token_starts_ring_number()? {
                        let ring_number = self.parse_ring_number()?;
                        state.add_single_bond_ring_marker(ring_number)?;
                    } else {
                        let atom = self.parse_simple_atomd()?;
                        state.add_single_bond_to_atom(atom)?;
                    }
                }
                SmilesToken::Separator => {
                    let atom = self.parse_simple_atomd()?;
                    state.add_disconnected_atom(atom)?;
                }
                SmilesToken::GroupOpen => {
                    let open_position =
                        SmilesBuildState::branch_open_token(self.lexer.current_token_position);
                    match self.next_token()? {
                        SmilesToken::Bond(bond) => {
                            let atom = self.parse_simple_atomd()?;
                            state.add_branch_explicit_bond(open_position, bond, atom)?;
                        }
                        SmilesToken::Minus => {
                            let atom = self.parse_simple_atomd()?;
                            state.add_branch_single_bond(open_position, atom)?;
                        }
                        SmilesToken::OrganicAtom(atom)
                        | SmilesToken::AromaticAtom(atom)
                        | SmilesToken::Atom(atom) => {
                            state.add_branch_atom_connected_to_active(open_position, atom)?;
                        }
                        SmilesToken::AtomOpen => {
                            let atom = self.parse_bracket_atomd()?;
                            state.add_branch_atom_connected_to_active(open_position, atom)?;
                        }
                        other => {
                            return Err(SmilesParseError::ParseError(format!(
                                "expected branch atom or bond, got {other:?}"
                            )));
                        }
                    }
                }
                SmilesToken::GroupClose => {
                    state.close_branch()?;
                }
                SmilesToken::Zero => {
                    state.add_ring_marker(0)?;
                }
                SmilesToken::NonzeroDigit(digit) => {
                    state.add_ring_marker(u32::try_from(digit).expect("digit is non-negative"))?;
                }
                SmilesToken::Percent => {
                    let ring_number = self.parse_percent_ring_number()?;
                    state.add_ring_marker(ring_number)?;
                }
                SmilesToken::BadCharacter(_) => {
                    return Err(SmilesParseError::ParseError("syntax error".to_string()));
                }
                _ => {
                    return Err(SmilesParseError::ParseError(
                        "unrecognized SMILES token".to_string(),
                    ));
                }
            }
        }
    }

    fn parse_simple_atomd(&mut self) -> Result<SmilesAtomToken, SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy atomd/simple_atom
        // RDKit✔️✔️: atomd:	simple_atom
        // RDKit✔️✔️: simple_atom:       ORGANIC_ATOM_TOKEN
        // RDKit✔️✔️:                 | AROMATIC_ATOM_TOKEN
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
        // END RDKIT CPP GRAMMAR ACTION smiles.yy atomd/simple_atom
        match self.next_token()? {
            SmilesToken::OrganicAtom(atom)
            | SmilesToken::AromaticAtom(atom)
            | SmilesToken::Atom(atom) => Ok(atom),
            SmilesToken::AtomOpen => self.parse_bracket_atomd(),
            SmilesToken::BadCharacter(_) => {
                Err(SmilesParseError::ParseError("syntax error".to_string()))
            }
            other => Err(SmilesParseError::ParseError(format!(
                "expected atom token, got {other:?}"
            ))),
        }
    }

    fn next_token(&mut self) -> Result<SmilesToken, SmilesParseError> {
        if let Some(token) = self.lookahead.take() {
            return Ok(token);
        }
        self.lexer.next_token()
    }

    fn peek_token(&mut self) -> Result<SmilesToken, SmilesParseError> {
        if self.lookahead.is_none() {
            self.lookahead = Some(self.lexer.next_token()?);
        }
        Ok(self
            .lookahead
            .clone()
            .expect("lookahead was just initialized"))
    }

    fn parse_bracket_atomd(&mut self) -> Result<SmilesAtomToken, SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy bracket atomd
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
        // END RDKIT CPP GRAMMAR ACTION smiles.yy bracket atomd
        let mut atom = self.parse_charge_element()?;
        atom.spec = atom.spec.with_no_implicit(true);
        match self.next_token()? {
            SmilesToken::Colon => {
                let map_number = self.parse_number()?;
                let map_number = u32::try_from(map_number).map_err(|_| {
                    SmilesParseError::ParseError(format!(
                        "atom map number out of range: {map_number}"
                    ))
                })?;
                atom.spec = atom.spec.with_atom_map(map_number);
                self.expect_atom_close()?;
                Ok(atom)
            }
            SmilesToken::AtomClose => Ok(atom),
            other => Err(SmilesParseError::ParseError(format!(
                "expected bracket atom close or atom map, got {other:?}"
            ))),
        }
    }

    fn parse_charge_element(&mut self) -> Result<SmilesAtomToken, SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy charge_element
        // RDKit✔️✔️: charge_element:	h_element
        // RDKit✔️✔️: | h_element PLUS_TOKEN { $1->setFormalCharge(1); }
        // RDKit✔️✔️: | h_element PLUS_TOKEN PLUS_TOKEN { $1->setFormalCharge(2); }
        // RDKit✔️✔️: | h_element PLUS_TOKEN number { $1->setFormalCharge($3); }
        // RDKit✔️✔️: | h_element MINUS_TOKEN { $1->setFormalCharge(-1); }
        // RDKit✔️✔️: | h_element MINUS_TOKEN MINUS_TOKEN { $1->setFormalCharge(-2); }
        // RDKit✔️✔️: | h_element MINUS_TOKEN number { $1->setFormalCharge(-$3); }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy charge_element
        let mut atom = self.parse_h_element()?;
        match self.peek_token()? {
            SmilesToken::Plus => {
                self.next_token()?;
                let charge = match self.peek_token()? {
                    SmilesToken::Plus => {
                        self.next_token()?;
                        2
                    }
                    SmilesToken::Zero | SmilesToken::NonzeroDigit(_) => self.parse_number()?,
                    _ => 1,
                };
                atom.spec = atom.spec.with_formal_charge(i8_from_rdkit_int(charge)?);
            }
            SmilesToken::Minus => {
                self.next_token()?;
                let charge = match self.peek_token()? {
                    SmilesToken::Minus => {
                        self.next_token()?;
                        -2
                    }
                    SmilesToken::Zero | SmilesToken::NonzeroDigit(_) => -self.parse_number()?,
                    _ => -1,
                };
                atom.spec = atom.spec.with_formal_charge(i8_from_rdkit_int(charge)?);
            }
            _ => {}
        }
        Ok(atom)
    }

    fn parse_h_element(&mut self) -> Result<SmilesAtomToken, SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy h_element
        // RDKit✔️✔️: h_element:      H_TOKEN { $$ = new Atom(1); }
        // RDKit✔️✔️:                 | number H_TOKEN { $$ = new Atom(1); $$->setIsotope($1); }
        // RDKit✔️✔️:                 | H_TOKEN H_TOKEN { $$ = new Atom(1); $$->setNumExplicitHs(1); }
        // RDKit✔️✔️:                 | number H_TOKEN H_TOKEN { $$ = new Atom(1); $$->setIsotope($1); $$->setNumExplicitHs(1);}
        // RDKit✔️✔️:                 | H_TOKEN H_TOKEN number { $$ = new Atom(1); $$->setNumExplicitHs($3); }
        // RDKit✔️✔️:                 | number H_TOKEN H_TOKEN number { $$ = new Atom(1); $$->setIsotope($1); $$->setNumExplicitHs($4);}
        // RDKit✔️✔️:                 | chiral_element
        // RDKit✔️✔️: 		| chiral_element H_TOKEN		{ $$ = $1; $1->setNumExplicitHs(1);}
        // RDKit✔️✔️: 		| chiral_element H_TOKEN number	{ $$ = $1; $1->setNumExplicitHs($3);}
        // END RDKIT CPP GRAMMAR ACTION smiles.yy h_element
        match self.peek_token()? {
            SmilesToken::H => {
                self.next_token()?;
                let mut atom = SmilesAtomToken::new(1);
                if matches!(self.peek_token()?, SmilesToken::H) {
                    self.next_token()?;
                    let explicit_hydrogens = if matches!(
                        self.peek_token()?,
                        SmilesToken::Zero | SmilesToken::NonzeroDigit(_)
                    ) {
                        self.parse_number()?
                    } else {
                        1
                    };
                    atom.spec = atom
                        .spec
                        .with_explicit_hydrogens(u8_from_rdkit_int(explicit_hydrogens)?);
                }
                Ok(atom)
            }
            SmilesToken::Zero | SmilesToken::NonzeroDigit(_) => {
                let number = self.parse_number()?;
                if matches!(self.peek_token()?, SmilesToken::H) {
                    self.next_token()?;
                    let mut atom = SmilesAtomToken::new(1);
                    atom.spec = atom.spec.with_isotope(u16_from_rdkit_int(number)?);
                    if matches!(self.peek_token()?, SmilesToken::H) {
                        self.next_token()?;
                        let explicit_hydrogens = if matches!(
                            self.peek_token()?,
                            SmilesToken::Zero | SmilesToken::NonzeroDigit(_)
                        ) {
                            self.parse_number()?
                        } else {
                            1
                        };
                        atom.spec = atom
                            .spec
                            .with_explicit_hydrogens(u8_from_rdkit_int(explicit_hydrogens)?);
                    }
                    Ok(atom)
                } else {
                    let atom = self.parse_chiral_element_after_leading_number(number)?;
                    self.parse_optional_explicit_hydrogen_suffix(atom)
                }
            }
            _ => {
                let atom = self.parse_chiral_element()?;
                self.parse_optional_explicit_hydrogen_suffix(atom)
            }
        }
    }

    fn parse_optional_explicit_hydrogen_suffix(
        &mut self,
        mut atom: SmilesAtomToken,
    ) -> Result<SmilesAtomToken, SmilesParseError> {
        if matches!(self.peek_token()?, SmilesToken::H) {
            self.next_token()?;
            let explicit_hydrogens = if matches!(
                self.peek_token()?,
                SmilesToken::Zero | SmilesToken::NonzeroDigit(_)
            ) {
                self.parse_number()?
            } else {
                1
            };
            atom.spec = atom
                .spec
                .with_explicit_hydrogens(u8_from_rdkit_int(explicit_hydrogens)?);
        }
        Ok(atom)
    }

    fn parse_chiral_element(&mut self) -> Result<SmilesAtomToken, SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy chiral_element
        // RDKit✔️✔️: chiral_element:	 element
        // RDKit✔️✔️: | element AT_TOKEN { $1->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW); }
        // RDKit✔️✔️: | element AT_TOKEN AT_TOKEN { $1->setChiralTag(Atom::CHI_TETRAHEDRAL_CW); }
        // RDKit✔️✔️: | element CHI_CLASS_TOKEN { $1->setChiralTag($2); $1->setProp(common_properties::_chiralPermutation,0); }
        // RDKit✔️✔️: | element CHI_CLASS_TOKEN number {
        // RDKit✔️✔️:     if($3==0){
        // RDKit✔️✔️:       yyerror(input,molList,branchPoints,scanner,start_token, current_token_position,
        // RDKit✔️✔️:             "chiral permutation cannot be zero");
        // RDKit✔️✔️:       yyErrorCleanup(molList);
        // RDKit✔️✔️:       delete $1;
        // RDKit✔️✔️:       YYABORT;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     $1->setChiralTag($2); $1->setProp(common_properties::_chiralPermutation,$3);
        // RDKit✔️✔️: };
        // END RDKIT CPP GRAMMAR ACTION smiles.yy chiral_element
        let atom = self.parse_element()?;
        self.parse_chiral_suffix(atom)
    }

    fn parse_chiral_element_after_leading_number(
        &mut self,
        number: i32,
    ) -> Result<SmilesAtomToken, SmilesParseError> {
        let atom = self.parse_element_after_leading_number(number)?;
        self.parse_chiral_suffix(atom)
    }

    fn parse_chiral_suffix(
        &mut self,
        mut atom: SmilesAtomToken,
    ) -> Result<SmilesAtomToken, SmilesParseError> {
        match self.peek_token()? {
            SmilesToken::At => {
                self.next_token()?;
                if matches!(self.peek_token()?, SmilesToken::At) {
                    self.next_token()?;
                    atom.spec = atom.spec.with_chiral_tag(ChiralTag::TetrahedralCw);
                } else {
                    atom.spec = atom.spec.with_chiral_tag(ChiralTag::TetrahedralCcw);
                }
            }
            SmilesToken::ChiralClass(chiral_tag) => {
                self.next_token()?;
                atom.spec = atom
                    .spec
                    .with_chiral_tag(chiral_tag)
                    .with_chiral_permutation(0);
                if matches!(
                    self.peek_token()?,
                    SmilesToken::Zero | SmilesToken::NonzeroDigit(_)
                ) {
                    let permutation = self.parse_number()?;
                    if permutation == 0 {
                        return Err(SmilesParseError::ParseError(
                            "chiral permutation cannot be zero".to_string(),
                        ));
                    }
                    atom.spec = atom.spec.with_chiral_permutation(
                        u32::try_from(permutation)
                            .expect("SMILES number parser returns non-negative values"),
                    );
                }
            }
            _ => {}
        }
        Ok(atom)
    }

    fn parse_element(&mut self) -> Result<SmilesAtomToken, SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy element
        // RDKit✔️✔️: element:	simple_atom
        // RDKit✔️✔️: 		|	number simple_atom { $2->setIsotope( $1 ); $$ = $2; }
        // RDKit✔️✔️: 		|	ATOM_TOKEN
        // RDKit✔️✔️: 		|	number ATOM_TOKEN	   { $2->setIsotope( $1 ); $$ = $2; }
        // RDKit✔️✔️: 		|	HASH_TOKEN	number   { $$ = new Atom($2); }
        // RDKit✔️✔️: 		|	number HASH_TOKEN	number   { $$ = new Atom($3); $$->setIsotope($1); }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy element
        match self.next_token()? {
            SmilesToken::OrganicAtom(atom)
            | SmilesToken::AromaticAtom(atom)
            | SmilesToken::Atom(atom) => Ok(atom),
            SmilesToken::Hash => {
                let atomic_number = self.parse_number()?;
                Ok(SmilesAtomToken::new(u8_from_rdkit_int(atomic_number)?))
            }
            SmilesToken::Zero => self.parse_element_after_leading_number(0),
            SmilesToken::NonzeroDigit(digit) => self.parse_number_tail_then_element(digit),
            other => Err(SmilesParseError::ParseError(format!(
                "expected element token, got {other:?}"
            ))),
        }
    }

    fn parse_number_tail_then_element(
        &mut self,
        first_digit: i32,
    ) -> Result<SmilesAtomToken, SmilesParseError> {
        let number = self.parse_number_tail(first_digit)?;
        self.parse_element_after_leading_number(number)
    }

    fn parse_element_after_leading_number(
        &mut self,
        number: i32,
    ) -> Result<SmilesAtomToken, SmilesParseError> {
        match self.next_token()? {
            SmilesToken::OrganicAtom(mut atom)
            | SmilesToken::AromaticAtom(mut atom)
            | SmilesToken::Atom(mut atom) => {
                atom.spec = atom.spec.with_isotope(u16_from_rdkit_int(number)?);
                Ok(atom)
            }
            SmilesToken::Hash => {
                let atomic_number = self.parse_number()?;
                let mut atom = SmilesAtomToken::new(u8_from_rdkit_int(atomic_number)?);
                atom.spec = atom.spec.with_isotope(u16_from_rdkit_int(number)?);
                Ok(atom)
            }
            other => Err(SmilesParseError::ParseError(format!(
                "expected element after isotope, got {other:?}"
            ))),
        }
    }

    fn parse_number(&mut self) -> Result<i32, SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy number/nonzero_number
        // RDKit✔️✔️: number:  ZERO_TOKEN
        // RDKit✔️✔️: | nonzero_number
        // RDKit✔️✔️: ;
        // RDKit✔️✔️: nonzero_number:  NONZERO_DIGIT_TOKEN
        // RDKit✔️✔️: | nonzero_number digit {
        // RDKit✔️✔️:   if($1 >= std::numeric_limits<std::int32_t>::max()/10 ||
        // RDKit✔️✔️:      ($1 == std::numeric_limits<std::int32_t>::max()/10 && $2 > std::numeric_limits<std::int32_t>::max()%10)) {
        // RDKit✔️✔️:     yyerror(input,molList,branchPoints,scanner,start_token,current_token_position,"number too large");
        // RDKit✔️✔️:     YYABORT;
        // RDKit✔️✔️:   }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy number/nonzero_number
        match self.next_token()? {
            SmilesToken::Zero => Ok(0),
            SmilesToken::NonzeroDigit(digit) => self.parse_number_tail(digit),
            other => Err(SmilesParseError::ParseError(format!(
                "expected number, got {other:?}"
            ))),
        }
    }

    fn parse_number_tail(&mut self, first_digit: i32) -> Result<i32, SmilesParseError> {
        let mut value = first_digit;
        while let SmilesToken::Zero | SmilesToken::NonzeroDigit(_) = self.peek_token()? {
            let digit = match self.next_token()? {
                SmilesToken::Zero => 0,
                SmilesToken::NonzeroDigit(digit) => digit,
                _ => unreachable!("peeked digit token changed"),
            };
            value = value
                .checked_mul(10)
                .and_then(|value| value.checked_add(digit))
                .ok_or_else(|| SmilesParseError::ParseError("number too large".to_string()))?;
        }
        Ok(value)
    }

    fn expect_atom_close(&mut self) -> Result<(), SmilesParseError> {
        match self.next_token()? {
            SmilesToken::AtomClose => Ok(()),
            other => Err(SmilesParseError::ParseError(format!(
                "expected atom close, got {other:?}"
            ))),
        }
    }

    fn next_token_starts_ring_number(&mut self) -> Result<bool, SmilesParseError> {
        Ok(matches!(
            self.peek_token()?,
            SmilesToken::Zero | SmilesToken::NonzeroDigit(_) | SmilesToken::Percent
        ))
    }

    fn parse_ring_number(&mut self) -> Result<u32, SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy ring_number
        // RDKit✔️✔️: ring_number:  digit
        // RDKit✔️✔️: | PERCENT_TOKEN NONZERO_DIGIT_TOKEN digit { $$ = $2*10+$3; }
        // RDKit✔️✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit GROUP_CLOSE_TOKEN { $$ = $3; }
        // RDKit✔️✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit GROUP_CLOSE_TOKEN { $$ = $3*10+$4; }
        // RDKit✔️✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit GROUP_CLOSE_TOKEN { $$ = $3*100+$4*10+$5; }
        // RDKit✔️✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit GROUP_CLOSE_TOKEN { $$ = $3*1000+$4*100+$5*10+$6; }
        // RDKit✔️✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit digit GROUP_CLOSE_TOKEN { $$ = $3*10000+$4*1000+$5*100+$6*10+$7; }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy ring_number
        match self.next_token()? {
            SmilesToken::Zero => Ok(0),
            SmilesToken::NonzeroDigit(digit) => {
                Ok(u32::try_from(digit).expect("digit is non-negative"))
            }
            SmilesToken::Percent => self.parse_percent_ring_number(),
            other => Err(SmilesParseError::ParseError(format!(
                "expected ring number, got {other:?}"
            ))),
        }
    }

    fn parse_percent_ring_number(&mut self) -> Result<u32, SmilesParseError> {
        match self.next_token()? {
            SmilesToken::NonzeroDigit(tens) => {
                let ones = self.parse_digit()?;
                Ok(u32::try_from(tens * 10 + ones).expect("ring number is non-negative"))
            }
            SmilesToken::GroupOpen => {
                let mut digits = Vec::new();
                while !matches!(self.peek_token()?, SmilesToken::GroupClose) {
                    if digits.len() == 5 {
                        return Err(SmilesParseError::ParseError(
                            "ring number too large".to_string(),
                        ));
                    }
                    digits.push(self.parse_digit()?);
                }
                self.next_token()?;
                if digits.is_empty() {
                    return Err(SmilesParseError::ParseError(
                        "empty ring number".to_string(),
                    ));
                }
                let mut value = 0_u32;
                for digit in digits {
                    value = value
                        .checked_mul(10)
                        .and_then(|value| value.checked_add(u32::try_from(digit).ok()?))
                        .ok_or_else(|| {
                            SmilesParseError::ParseError("ring number too large".to_string())
                        })?;
                }
                Ok(value)
            }
            other => Err(SmilesParseError::ParseError(format!(
                "expected percent ring number, got {other:?}"
            ))),
        }
    }

    fn parse_digit(&mut self) -> Result<i32, SmilesParseError> {
        match self.next_token()? {
            SmilesToken::Zero => Ok(0),
            SmilesToken::NonzeroDigit(digit) => Ok(digit),
            other => Err(SmilesParseError::ParseError(format!(
                "expected digit, got {other:?}"
            ))),
        }
    }
}

fn i8_from_rdkit_int(value: i32) -> Result<i8, SmilesParseError> {
    i8::try_from(value)
        .map_err(|_| SmilesParseError::ParseError(format!("value out of i8 range: {value}")))
}

fn u8_from_rdkit_int(value: i32) -> Result<u8, SmilesParseError> {
    u8::try_from(value)
        .map_err(|_| SmilesParseError::ParseError(format!("value out of u8 range: {value}")))
}

fn u16_from_rdkit_int(value: i32) -> Result<u16, SmilesParseError> {
    u16::try_from(value)
        .map_err(|_| SmilesParseError::ParseError(format!("value out of u16 range: {value}")))
}

fn canonical_bond_pair(begin: AtomId, end: AtomId) -> (usize, usize) {
    let begin = begin.index();
    let end = end.index();
    if begin <= end {
        (begin, end)
    } else {
        (end, begin)
    }
}

fn chiral_permutation_limit(chiral_tag: ChiralTag) -> Option<i32> {
    match chiral_tag {
        ChiralTag::Tetrahedral | ChiralTag::Allene => Some(2),
        ChiralTag::SquarePlanar => Some(3),
        ChiralTag::TrigonalBipyramidal => Some(20),
        ChiralTag::Octahedral => Some(30),
        _ => None,
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct BranchPoint {
    atom: AtomId,
    open_position: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PendingBond {
    token: SmilesBondToken,
    cx_smiles_bond_idx: Option<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RingOpening {
    atom: AtomId,
    pending_bond: Option<PendingBond>,
    input_position: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RingClosureRecord {
    ring_number: u32,
    bond: Option<BondId>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PendingRingClosure {
    ring_number: u32,
    opening_atom: AtomId,
    closing_atom: AtomId,
    opening_pending_bond: PendingBond,
    closing_pending_bond: PendingBond,
}

#[derive(Debug)]
#[allow(dead_code)]
struct SmilesBuildState {
    builder: MoleculeBuilder,
    active_atom: Option<AtomId>,
    atom_aromatic: Vec<bool>,
    bond_pairs: BTreeSet<(usize, usize)>,
    explicit_unspecified_bonds: Vec<BondId>,
    branch_stack: Vec<BranchPoint>,
    ring_openings: BTreeMap<u32, RingOpening>,
    ring_closures_by_atom: BTreeMap<AtomId, Vec<RingClosureRecord>>,
    pending_ring_closures: Vec<PendingRingClosure>,
    smiles_start_atoms: BTreeSet<AtomId>,
    cx_bond_order: Vec<BondId>,
    next_cx_smiles_bond_idx: usize,
    temporary_chiral_permutations: BTreeMap<AtomId, u32>,
    cx_stereo_group_tracker: BTreeMap<(u8, u32), usize>,
}

#[allow(dead_code)]
impl SmilesBuildState {
    fn new() -> Self {
        Self {
            builder: MoleculeBuilder::new(),
            active_atom: None,
            atom_aromatic: Vec::new(),
            bond_pairs: BTreeSet::new(),
            explicit_unspecified_bonds: Vec::new(),
            branch_stack: Vec::new(),
            ring_openings: BTreeMap::new(),
            ring_closures_by_atom: BTreeMap::new(),
            pending_ring_closures: Vec::new(),
            smiles_start_atoms: BTreeSet::new(),
            cx_bond_order: Vec::new(),
            next_cx_smiles_bond_idx: 0,
            temporary_chiral_permutations: BTreeMap::new(),
            cx_stereo_group_tracker: BTreeMap::new(),
        }
    }

    fn cleanup_after_parse_error(&mut self) {
        // BEGIN RDKIT CPP FUNCTION CleanupAfterParseError
        // RDKit✔️✔️: void CleanupAfterParseError(RWMol *mol) {
        // RDKit✔️✔️:   PRECONDITION(mol, "no molecule");
        // RDKit✔️✔️:   // blow out any partial bonds:
        // RDKit✔️✔️:   for (auto markI : *mol->getBondBookmarks()) {
        // RDKit✔️✔️:     RWMol::BOND_PTR_LIST &bonds = markI.second;
        // RDKit✔️✔️:     for (auto &bond : bonds) {
        // RDKit✔️✔️:       delete bond;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION CleanupAfterParseError
        self.ring_openings.clear();
        self.ring_closures_by_atom.clear();
        self.pending_ring_closures.clear();
    }

    fn add_frag_to_mol(
        &mut self,
        frag: SmilesBuildState,
        bond_order: BondOrder,
        bond_dir: BondDirection,
    ) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP FUNCTION AddFragToMol
        // RDKit✔️✔️: void AddFragToMol(RWMol *mol, RWMol *frag, Bond::BondType bondOrder,
        // RDKit✔️✔️:                   Bond::BondDir bondDir) {
        // RDKit✔️✔️:   PRECONDITION(mol, "no molecule");
        // RDKit✔️✔️:   PRECONDITION(frag, "no fragment");
        // RDKit✔️✔️:   PRECONDITION(mol->getActiveAtom(), "no active atom");
        // RDKit✔️✔️:   Atom *lastAt = mol->getActiveAtom();
        // RDKit✔️✔️:   int nOrigAtoms = mol->getNumAtoms();
        // RDKit✔️✔️:   int nOrigBonds = mol->getNumBonds();
        // RDKit✔️✔️:   mol->insertMol(*frag);
        // RDKit✔️✔️:   // update ring-closure order information on the added atoms
        // RDKit✔️✔️:   // and copy fragment bookmarks/partial bonds into the destination.
        // RDKit✔️✔️:   // When bondOrder != IONIC, connect the former active atom to the first
        // RDKit✔️✔️:   // fragment atom using the requested order and direction.
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION AddFragToMol
        let last_active = self
            .active_atom
            .ok_or_else(|| SmilesParseError::ParseError("no active atom".to_string()))?;
        let first_old_atom = frag
            .builder
            .atoms()
            .first()
            .map(Atom::id)
            .ok_or_else(|| SmilesParseError::ParseError("fragment has no atoms".to_string()))?;

        let mut atom_mapping = Vec::with_capacity(frag.builder.atoms().len());
        for atom in frag.builder.atoms() {
            let new_atom = self.builder.add_atom(atom_spec_from_atom(atom));
            atom_mapping.push(new_atom);
            self.atom_aromatic.push(atom.is_aromatic());
            if frag.smiles_start_atoms.contains(&atom.id()) {
                self.smiles_start_atoms.insert(new_atom);
            }
        }

        let mut bond_mapping = Vec::with_capacity(frag.builder.bonds().len());
        for bond in frag.builder.bonds() {
            let begin = atom_mapping[bond.begin().index()];
            let end = atom_mapping[bond.end().index()];
            let mut spec = BondSpec::new(begin, end, bond.order())
                .with_aromatic(bond.is_aromatic())
                .with_conjugated(bond.is_conjugated())
                .with_direction(bond.direction())
                .with_stereo(bond.stereo());
            if let Some([stereo_begin, stereo_end]) = bond.stereo_atoms() {
                spec = spec.with_stereo_atoms(
                    atom_mapping[stereo_begin.index()],
                    atom_mapping[stereo_end.index()],
                );
            }
            if bond.unknown_stereo() {
                spec = spec.with_unknown_stereo(true);
            }
            if let Some(query) = bond.query().cloned() {
                spec = spec.with_query(query);
            }
            for (key, value) in bond.props() {
                if key == CXSMILES_BOND_IDX_PROP {
                    continue;
                }
                spec = spec.with_prop(key.clone(), value.clone());
            }
            spec = spec.with_prop(CXSMILES_BOND_IDX_PROP, self.cx_bond_order.len().to_string());
            let bond_id = self
                .builder
                .add_bond(spec)
                .map_err(|error| SmilesParseError::ParseError(error.to_string()))?;
            self.bond_pairs.insert(canonical_bond_pair(begin, end));
            if bond.prop(UNSPECIFIED_ORDER_PROP).is_some() {
                self.explicit_unspecified_bonds.push(bond_id);
            }
            self.cx_bond_order.push(bond_id);
            bond_mapping.push(bond_id);
        }

        for (old_atom, records) in frag.ring_closures_by_atom {
            let new_atom = atom_mapping[old_atom.index()];
            let mut new_records = Vec::with_capacity(records.len());
            for record in records {
                new_records.push(RingClosureRecord {
                    ring_number: record.ring_number,
                    bond: record.bond.map(|bond| bond_mapping[bond.index()]),
                });
            }
            self.ring_closures_by_atom.insert(new_atom, new_records);
        }
        for (ring_number, opening) in frag.ring_openings {
            self.ring_openings.insert(
                ring_number,
                RingOpening {
                    atom: atom_mapping[opening.atom.index()],
                    pending_bond: opening.pending_bond,
                    input_position: opening.input_position,
                },
            );
        }
        for closure in frag.pending_ring_closures {
            self.pending_ring_closures.push(PendingRingClosure {
                ring_number: closure.ring_number,
                opening_atom: atom_mapping[closure.opening_atom.index()],
                closing_atom: atom_mapping[closure.closing_atom.index()],
                opening_pending_bond: closure.opening_pending_bond,
                closing_pending_bond: closure.closing_pending_bond,
            });
        }
        for (atom, permutation) in frag.temporary_chiral_permutations {
            self.temporary_chiral_permutations
                .insert(atom_mapping[atom.index()], permutation);
        }
        for (key, value) in frag.cx_stereo_group_tracker {
            self.cx_stereo_group_tracker.insert(key, value);
        }

        if bond_order != BondOrder::Ionic {
            let mut token = SmilesBondToken::new(bond_order);
            token.direction = bond_dir;
            if bond_order == BondOrder::Unspecified {
                let begin_aromatic = self
                    .atom_aromatic
                    .get(last_active.index())
                    .copied()
                    .unwrap_or(false);
                let end_aromatic = self
                    .atom_aromatic
                    .get(atom_mapping[first_old_atom.index()].index())
                    .copied()
                    .unwrap_or(false);
                token.order = get_unspecified_bond_type_for_atoms(begin_aromatic, end_aromatic);
                token.explicit_unspecified_order = true;
            }
            let bond_id =
                self.add_bond_from_token(last_active, atom_mapping[first_old_atom.index()], token)?;
            self.cx_bond_order.push(bond_id);
            self.active_atom = Some(atom_mapping[first_old_atom.index()]);
        } else {
            self.active_atom = Some(atom_mapping[first_old_atom.index()]);
        }
        Ok(())
    }

    fn add_first_atom(&mut self, mut atom: SmilesAtomToken) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy mol: atomd
        // RDKit✔️✔️: mol: atomd {
        // RDKit✔️✔️:   int sz     = molList->size();
        // RDKit✔️✔️:   molList->resize( sz + 1);
        // RDKit✔️✔️:   (*molList)[ sz ] = new RWMol();
        // RDKit✔️✔️:   RDKit::RWMol *curMol = (*molList)[ sz ];
        // RDKit✔️✔️:   $1->setProp(RDKit::common_properties::_SmilesStart,1);
        // RDKit✔️✔️:   curMol->addAtom($1, true, true);
        // RDKit✔️✔️:   //delete $1;
        // RDKit✔️✔️:   $$ = sz;
        // RDKit✔️✔️: }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy mol: atomd
        atom.spec = atom.spec.with_prop(SMILES_START_PROP, "1");
        let atom_id = self.append_atom(atom);
        self.active_atom = Some(atom_id);
        self.smiles_start_atoms.insert(atom_id);
        Ok(())
    }

    fn add_atom_connected_to_active(
        &mut self,
        atom: SmilesAtomToken,
    ) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy mol atomd
        // RDKit✔️✔️: | mol atomd       {
        // RDKit✔️✔️:   RWMol *mp = (*molList)[$$];
        // RDKit✔️✔️:   Atom *a1 = mp->getActiveAtom();
        // RDKit✔️✔️:   int atomIdx1=a1->getIdx();
        // RDKit✔️✔️:   int atomIdx2=mp->addAtom($2,true,true);
        // RDKit✔️✔️:   mp->addBond(atomIdx1,atomIdx2,
        // RDKit✔️✔️: 	      SmilesParseOps::GetUnspecifiedBondType(mp,a1,mp->getAtomWithIdx(atomIdx2)));
        // RDKit✔️✔️:   mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit✔️✔️:   //delete $2;
        // RDKit✔️✔️: }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy mol atomd
        let begin = self
            .active_atom
            .ok_or_else(|| SmilesParseError::ParseError("no active atom".to_string()))?;
        let begin_aromatic = self
            .atom_aromatic
            .get(begin.index())
            .copied()
            .unwrap_or(false);
        let end_aromatic = atom.spec.is_aromatic();
        let end = self.append_atom(atom);
        let order = get_unspecified_bond_type_for_atoms(begin_aromatic, end_aromatic);
        let bond_id = self.add_bond_with_cx_index(begin, end, order)?;
        self.cx_bond_order.push(bond_id);
        self.active_atom = Some(end);
        Ok(())
    }

    fn add_explicit_bond_to_atom(
        &mut self,
        bond: SmilesBondToken,
        atom: SmilesAtomToken,
    ) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy mol BOND_TOKEN atomd
        // RDKit✔️✔️: | mol BOND_TOKEN atomd  {
        // RDKit✔️✔️:   RWMol *mp = (*molList)[$$];
        // RDKit✔️✔️:   int atomIdx1 = mp->getActiveAtom()->getIdx();
        // RDKit✔️✔️:   int atomIdx2 = mp->addAtom($3,true,true);
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
        // RDKit✔️✔️:   $2->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit✔️✔️:   mp->addBond($2,true);
        // RDKit✔️✔️:   //delete $3;
        // RDKit✔️✔️: }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy mol BOND_TOKEN atomd
        let atom_idx1 = self
            .active_atom
            .ok_or_else(|| SmilesParseError::ParseError("no active atom".to_string()))?;
        let atom_idx2 = self.append_atom(atom);
        let bond_id = self.add_bond_from_token(atom_idx1, atom_idx2, bond)?;
        self.cx_bond_order.push(bond_id);
        self.active_atom = Some(atom_idx2);
        Ok(())
    }

    fn add_single_bond_to_atom(&mut self, atom: SmilesAtomToken) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy mol MINUS_TOKEN atomd
        // RDKit✔️✔️: | mol MINUS_TOKEN atomd {
        // RDKit✔️✔️:   RWMol *mp = (*molList)[$$];
        // RDKit✔️✔️:   int atomIdx1 = mp->getActiveAtom()->getIdx();
        // RDKit✔️✔️:   int atomIdx2 = mp->addAtom($3,true,true);
        // RDKit✔️✔️:   mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
        // RDKit✔️✔️:   mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit✔️✔️:   //delete $3;
        // RDKit✔️✔️: }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy mol MINUS_TOKEN atomd
        let begin = self
            .active_atom
            .ok_or_else(|| SmilesParseError::ParseError("no active atom".to_string()))?;
        let end = self.append_atom(atom);
        let bond_id = self.add_bond_with_cx_index(begin, end, BondOrder::Single)?;
        self.cx_bond_order.push(bond_id);
        self.active_atom = Some(end);
        Ok(())
    }

    fn add_disconnected_atom(&mut self, mut atom: SmilesAtomToken) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy mol SEPARATOR_TOKEN atomd
        // RDKit✔️✔️: | mol SEPARATOR_TOKEN atomd {
        // RDKit✔️✔️:   RWMol *mp = (*molList)[$$];
        // RDKit✔️✔️:   $3->setProp(RDKit::common_properties::_SmilesStart,1,true);
        // RDKit✔️✔️:   mp->addAtom($3,true,true);
        // RDKit✔️✔️: }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy mol SEPARATOR_TOKEN atomd
        atom.spec = atom.spec.with_prop(SMILES_START_PROP, "1");
        let atom_id = self.append_atom(atom);
        self.active_atom = Some(atom_id);
        self.smiles_start_atoms.insert(atom_id);
        Ok(())
    }

    fn branch_open_token(current_token_position: usize) -> usize {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy branch_open_token
        // RDKit✔️✔️: branch_open_token: GROUP_OPEN_TOKEN { $$ = current_token_position; };
        // END RDKIT CPP GRAMMAR ACTION smiles.yy branch_open_token
        current_token_position
    }

    fn add_branch_atom_connected_to_active(
        &mut self,
        open_position: usize,
        atom: SmilesAtomToken,
    ) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy mol branch_open_token atomd
        // RDKit✔️✔️: | mol branch_open_token atomd {
        // RDKit✔️✔️:   RWMol *mp = (*molList)[$$];
        // RDKit✔️✔️:   Atom *a1 = mp->getActiveAtom();
        // RDKit✔️✔️:   int atomIdx1=a1->getIdx();
        // RDKit✔️✔️:   int atomIdx2=mp->addAtom($3,true,true);
        // RDKit✔️✔️:   mp->addBond(atomIdx1,atomIdx2,
        // RDKit✔️✔️: 	      SmilesParseOps::GetUnspecifiedBondType(mp,a1,mp->getAtomWithIdx(atomIdx2)));
        // RDKit✔️✔️:   mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit✔️✔️:   branchPoints.push_back({atomIdx1, $2});
        // RDKit✔️✔️: }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy mol branch_open_token atomd
        let branch_root = self
            .active_atom
            .ok_or_else(|| SmilesParseError::ParseError("no active atom".to_string()))?;
        self.add_atom_connected_to_active(atom)?;
        self.branch_stack.push(BranchPoint {
            atom: branch_root,
            open_position,
        });
        Ok(())
    }

    fn add_branch_explicit_bond(
        &mut self,
        open_position: usize,
        bond: SmilesBondToken,
        atom: SmilesAtomToken,
    ) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy mol branch_open_token BOND_TOKEN atomd
        // RDKit✔️✔️: | mol branch_open_token BOND_TOKEN atomd  {
        // RDKit✔️✔️:   RWMol *mp = (*molList)[$$];
        // RDKit✔️✔️:   int atomIdx1 = mp->getActiveAtom()->getIdx();
        // RDKit✔️✔️:   int atomIdx2 = mp->addAtom($4,true,true);
        // RDKit✔️✔️:   if( $3->getBondType() == Bond::DATIVER ){
        // RDKit✔️✔️:     $3->setBeginAtomIdx(atomIdx1);
        // RDKit✔️✔️:     $3->setEndAtomIdx(atomIdx2);
        // RDKit✔️✔️:     $3->setBondType(Bond::DATIVE);
        // RDKit✔️✔️:   }else if ( $3->getBondType() == Bond::DATIVEL ){
        // RDKit✔️✔️:     $3->setBeginAtomIdx(atomIdx2);
        // RDKit✔️✔️:     $3->setEndAtomIdx(atomIdx1);
        // RDKit✔️✔️:     $3->setBondType(Bond::DATIVE);
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     $3->setBeginAtomIdx(atomIdx1);
        // RDKit✔️✔️:     $3->setEndAtomIdx(atomIdx2);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   $3->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit✔️✔️:   mp->addBond($3,true);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   branchPoints.push_back({atomIdx1, $2});
        // RDKit✔️✔️: }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy mol branch_open_token BOND_TOKEN atomd
        let branch_root = self
            .active_atom
            .ok_or_else(|| SmilesParseError::ParseError("no active atom".to_string()))?;
        self.add_explicit_bond_to_atom(bond, atom)?;
        self.branch_stack.push(BranchPoint {
            atom: branch_root,
            open_position,
        });
        Ok(())
    }

    fn add_branch_single_bond(
        &mut self,
        open_position: usize,
        atom: SmilesAtomToken,
    ) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy mol branch_open_token MINUS_TOKEN atomd
        // RDKit✔️✔️: | mol branch_open_token MINUS_TOKEN atomd {
        // RDKit✔️✔️:   RWMol *mp = (*molList)[$$];
        // RDKit✔️✔️:   int atomIdx1 = mp->getActiveAtom()->getIdx();
        // RDKit✔️✔️:   int atomIdx2=mp->addAtom($4,true,true);
        // RDKit✔️✔️:   mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
        // RDKit✔️✔️:   mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit✔️✔️:   branchPoints.push_back({atomIdx1, $2});
        // RDKit✔️✔️: }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy mol branch_open_token MINUS_TOKEN atomd
        let branch_root = self
            .active_atom
            .ok_or_else(|| SmilesParseError::ParseError("no active atom".to_string()))?;
        self.add_single_bond_to_atom(atom)?;
        self.branch_stack.push(BranchPoint {
            atom: branch_root,
            open_position,
        });
        Ok(())
    }

    fn close_branch(&mut self) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy mol GROUP_CLOSE_TOKEN
        // RDKit✔️✔️: | mol GROUP_CLOSE_TOKEN {
        // RDKit✔️✔️:   if(branchPoints.empty()){
        // RDKit✔️✔️:      yyerror(input,molList,branchPoints,scanner,start_token,current_token_position,"extra close parentheses");
        // RDKit✔️✔️:      yyErrorCleanup(molList);
        // RDKit✔️✔️:      YYABORT;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   RWMol *mp = (*molList)[$$];
        // RDKit✔️✔️:
        // RDKit✔️✔️:   mp->setActiveAtom(branchPoints.back().first);
        // RDKit✔️✔️:   branchPoints.pop_back();
        // RDKit✔️✔️: }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy mol GROUP_CLOSE_TOKEN
        let branch_point = self
            .branch_stack
            .pop()
            .ok_or_else(|| SmilesParseError::ParseError("extra close parentheses".to_string()))?;
        self.active_atom = Some(branch_point.atom);
        Ok(())
    }

    fn add_ring_marker(&mut self, ring_number: u32) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy mol ring_number
        // RDKit✔️✔️: | mol ring_number {
        // RDKit✔️✔️:   RWMol * mp = (*molList)[$$];
        // RDKit✔️✔️:   Atom *atom=mp->getActiveAtom();
        // RDKit✔️✔️:   mp->setAtomBookmark(atom,$2);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   Bond *newB = mp->createPartialBond(atom->getIdx(),
        // RDKit✔️✔️: 				     Bond::UNSPECIFIED);
        // RDKit✔️✔️:   mp->setBondBookmark(newB,$2);
        // RDKit✔️✔️:   newB->setProp(RDKit::common_properties::_unspecifiedOrder,1);
        // RDKit✔️✔️:   if(!(mp->getAllBondsWithBookmark($2).size()%2)){
        // RDKit✔️✔️:     newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   INT_VECT tmp;
        // RDKit✔️✔️:   atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
        // RDKit✔️✔️:   tmp.push_back(-($2+1));
        // RDKit✔️✔️:   atom->setProp(RDKit::common_properties::_RingClosures,tmp);
        // RDKit✔️✔️: }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy mol ring_number
        let pending_bond = PendingBond {
            token: SmilesBondToken {
                order: BondOrder::Unspecified,
                is_aromatic: false,
                direction: BondDirection::None,
                explicit_unspecified_order: true,
                is_null_query: false,
            },
            cx_smiles_bond_idx: None,
        };
        self.add_ring_marker_with_pending_bond(ring_number, pending_bond)
    }

    fn add_explicit_bond_ring_marker(
        &mut self,
        bond: SmilesBondToken,
        ring_number: u32,
    ) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy mol BOND_TOKEN ring_number
        // RDKit✔️✔️: | mol BOND_TOKEN ring_number {
        // RDKit✔️✔️:   RWMol * mp = (*molList)[$$];
        // RDKit✔️✔️:   Atom *atom=mp->getActiveAtom();
        // RDKit✔️✔️:   Bond *newB = mp->createPartialBond(atom->getIdx(),
        // RDKit✔️✔️: 				     $2->getBondType());
        // RDKit✔️✔️:   if($2->hasProp(RDKit::common_properties::_unspecifiedOrder)){
        // RDKit✔️✔️:     newB->setProp(RDKit::common_properties::_unspecifiedOrder,1);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   newB->setBondDir($2->getBondDir());
        // RDKit✔️✔️:   mp->setAtomBookmark(atom,$3);
        // RDKit✔️✔️:   mp->setBondBookmark(newB,$3);
        // RDKit✔️✔️:   if(!(mp->getAllBondsWithBookmark($3).size()%2)){
        // RDKit✔️✔️:     newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);
        // RDKit✔️✔️:   INT_VECT tmp;
        // RDKit✔️✔️:   atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
        // RDKit✔️✔️:   tmp.push_back(-($3+1));
        // RDKit✔️✔️:   atom->setProp(RDKit::common_properties::_RingClosures,tmp);
        // RDKit✔️✔️:   delete $2;
        // RDKit✔️✔️: }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy mol BOND_TOKEN ring_number
        self.add_ring_marker_with_pending_bond(
            ring_number,
            PendingBond {
                token: bond,
                cx_smiles_bond_idx: None,
            },
        )
    }

    fn add_single_bond_ring_marker(&mut self, ring_number: u32) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy mol MINUS_TOKEN ring_number
        // RDKit✔️✔️: | mol MINUS_TOKEN ring_number {
        // RDKit✔️✔️:   RWMol * mp = (*molList)[$$];
        // RDKit✔️✔️:   Atom *atom=mp->getActiveAtom();
        // RDKit✔️✔️:   Bond *newB = mp->createPartialBond(atom->getIdx(),
        // RDKit✔️✔️: 				     Bond::SINGLE);
        // RDKit✔️✔️:   mp->setAtomBookmark(atom,$3);
        // RDKit✔️✔️:   mp->setBondBookmark(newB,$3);
        // RDKit✔️✔️:   if(!(mp->getAllBondsWithBookmark($3).size()%2)){
        // RDKit✔️✔️:     newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);
        // RDKit✔️✔️:   INT_VECT tmp;
        // RDKit✔️✔️:   atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
        // RDKit✔️✔️:   tmp.push_back(-($3+1));
        // RDKit✔️✔️:   atom->setProp(RDKit::common_properties::_RingClosures,tmp);
        // RDKit✔️✔️: }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy mol MINUS_TOKEN ring_number
        self.add_ring_marker_with_pending_bond(
            ring_number,
            PendingBond {
                token: SmilesBondToken::new(BondOrder::Single),
                cx_smiles_bond_idx: None,
            },
        )
    }

    fn add_ring_marker_with_pending_bond(
        &mut self,
        ring_number: u32,
        mut pending_bond: PendingBond,
    ) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy ring closure _cxsmilesBondIdx assignment
        // RDKit✔️✔️:   if(!(mp->getAllBondsWithBookmark($2).size()%2)){
        // RDKit✔️✔️:     newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if(!(mp->getAllBondsWithBookmark($3).size()%2)){
        // RDKit✔️✔️:     newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit✔️✔️:   }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy ring closure _cxsmilesBondIdx assignment
        let atom = self
            .active_atom
            .ok_or_else(|| SmilesParseError::ParseError("no active atom".to_string()))?;
        self.check_ring_closure_branch_status(atom)?;
        if self.ring_openings.contains_key(&ring_number) {
            pending_bond.cx_smiles_bond_idx = Some(self.next_cx_smiles_bond_idx);
            self.next_cx_smiles_bond_idx += 1;
        }
        if let Some(opening) = self.ring_openings.remove(&ring_number) {
            let opening_pending_bond = opening.pending_bond.ok_or_else(|| {
                SmilesParseError::ParseError(format!("missing ring bond {ring_number}"))
            })?;
            self.pending_ring_closures.push(PendingRingClosure {
                ring_number,
                opening_atom: opening.atom,
                closing_atom: atom,
                opening_pending_bond,
                closing_pending_bond: pending_bond,
            });
            self.ring_closures_by_atom
                .entry(atom)
                .or_default()
                .push(RingClosureRecord {
                    ring_number,
                    bond: None,
                });
        } else {
            self.ring_openings.insert(
                ring_number,
                RingOpening {
                    atom,
                    pending_bond: Some(pending_bond),
                    input_position: 0,
                },
            );
            self.ring_closures_by_atom
                .entry(atom)
                .or_default()
                .push(RingClosureRecord {
                    ring_number,
                    bond: None,
                });
        }
        Ok(())
    }

    fn close_ring_opening(
        &mut self,
        opening_atom: AtomId,
        closing_atom: AtomId,
        opening_pending_bond: PendingBond,
        current_pending_bond: PendingBond,
    ) -> Result<BondId, SmilesParseError> {
        // BEGIN RDKIT CPP FUNCTION CloseMolRings bond selection + bond orientation section
        // RDKit✔️✔️:           // figure out which (if either) bond has a specified type, we'll
        // RDKit✔️✔️:           // keep that one.  We also need to update the end atom index to
        // RDKit✔️✔️:           // match FIX: daylight barfs when you give it multiple specs for the
        // RDKit✔️✔️:           // closure
        // RDKit✔️✔️:           //   bond, we'll just take the first one and ignore others
        // RDKit✔️✔️:           //   NOTE: we used to do this the other way (take the last
        // RDKit✔️✔️:           //   specification),
        // RDKit✔️✔️:           //   but that turned out to be troublesome in odd cases like
        // RDKit✔️✔️:           //   C1CC11CC1.
        // RDKit✔️✔️:           if (!bond1->hasProp(common_properties::_unspecifiedOrder)) {
        // RDKit✔️✔️:             matchedBond = bond1;
        // RDKit✔️✔️:             matchedBond->setEndAtomIdx(atom2->getIdx());
        // RDKit✔️✔️:           } else {
        // RDKit✔️✔️:             matchedBond = bond2;
        // RDKit✔️✔️:             matchedBond->setEndAtomIdx(atom1->getIdx());
        // RDKit✔️✔️:           }
        // END RDKIT CPP FUNCTION CloseMolRings bond selection + bond orientation section
        //
        // Rust stores pending bond specs as a fixed-size value token
        // (`BondOrder`/`BondDirection`/bool flags only), so choosing the
        // retained token here remains O(1) with no heap allocation.
        if opening_atom == closing_atom {
            return Err(SmilesParseError::ParseError(format!(
                "duplicated ring closure bonds atom {} to itself",
                opening_atom.index()
            )));
        }
        if self
            .bond_pairs
            .contains(&canonical_bond_pair(opening_atom, closing_atom))
        {
            return Err(SmilesParseError::ParseError(format!(
                "ring closure duplicates bond between atom {} and atom {}",
                opening_atom.index(),
                closing_atom.index()
            )));
        }
        let opening_aromatic = self
            .atom_aromatic
            .get(opening_atom.index())
            .copied()
            .unwrap_or(false);
        let closing_aromatic = self
            .atom_aromatic
            .get(closing_atom.index())
            .copied()
            .unwrap_or(false);
        let use_closing_token = opening_pending_bond.token.explicit_unspecified_order;
        let opening_token = opening_pending_bond.token;
        let closing_token = current_pending_bond.token;
        let mut token = if use_closing_token {
            closing_token
        } else {
            opening_token
        };
        let (begin_atom, end_atom) = if use_closing_token {
            (closing_atom, opening_atom)
        } else {
            (opening_atom, closing_atom)
        };
        token.direction = if use_closing_token {
            swap_bond_direction_if_needed(
                token.direction,
                opening_token.direction,
                closing_atom,
                opening_atom,
            )
        } else {
            swap_bond_direction_if_needed(
                token.direction,
                closing_token.direction,
                opening_atom,
                closing_atom,
            )
        };
        if token.order == BondOrder::Unspecified {
            token.order = get_unspecified_bond_type_for_atoms(opening_aromatic, closing_aromatic);
            token.explicit_unspecified_order = false;
        }
        // BEGIN RDKIT CPP FUNCTION CloseMolRings _cxsmilesBondIdx transfer
        // RDKit✔️✔️:          // we use the _cxsmilesBondIdx value from the second one, if it's
        // RDKit✔️✔️:          // there
        // RDKit✔️✔️:          if (bond2->hasProp("_cxsmilesBondIdx")) {
        // RDKit✔️✔️:            bond1->setProp("_cxsmilesBondIdx",
        // RDKit✔️✔️:                           bond2->getProp<unsigned int>("_cxsmilesBondIdx"));
        // RDKit✔️✔️:          }
        // END RDKIT CPP FUNCTION CloseMolRings _cxsmilesBondIdx transfer
        let cx_smiles_bond_idx = current_pending_bond
            .cx_smiles_bond_idx
            .or(opening_pending_bond.cx_smiles_bond_idx);
        let bond_id =
            self.add_bond_from_token_with_cx_idx(begin_atom, end_atom, token, cx_smiles_bond_idx)?;
        self.cx_bond_order.push(bond_id);
        Ok(bond_id)
    }

    fn check_ring_closure_branch_status(&mut self, atom: AtomId) -> Result<(), SmilesParseError> {
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
        let num_atoms = self.builder.atoms().len();
        let degree = self.atom_degree(atom);
        let should_invert = atom.index() != num_atoms.saturating_sub(1)
            && (degree == 1
                || (degree == 2 && atom.index() != 0)
                || (degree == 3 && atom.index() == 0));
        if !should_invert {
            return Ok(());
        }
        let chiral_tag = self
            .builder
            .atoms()
            .get(atom.index())
            .map(|atom| atom.chiral_tag())
            .ok_or_else(|| SmilesParseError::ParseError(format!("missing atom {atom}")))?;
        let inverted = match chiral_tag {
            ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
            ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
            _ => return Ok(()),
        };
        self.builder
            .atom_mut(atom)
            .ok_or_else(|| SmilesParseError::ParseError(format!("missing atom {atom}")))?
            .set_chiral_tag(inverted);
        Ok(())
    }

    fn atom_degree(&self, atom: AtomId) -> usize {
        self.builder.degree(atom)
    }

    fn finish_parse(&mut self) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP FUNCTION toMol post-parse section
        // RDKit✔️✔️:     func(inp, molVect);
        // RDKit✔️✔️:     if (!molVect.empty()) {
        // RDKit✔️✔️:       res.reset(molVect[0]);
        // RDKit✔️✔️:       SmilesParseOps::CloseMolRings(res.get(), false);
        // RDKit✔️✔️:       SmilesParseOps::CheckChiralitySpecifications(res.get(), true);
        // RDKit✔️✔️:       SmilesParseOps::SetUnspecifiedBondTypes(res.get());
        // RDKit✔️✔️:       SmilesParseOps::AdjustAtomChiralityFlags(res.get());
        // RDKit✔️✔️:       // No sense leaving this bookmark intact:
        // RDKit✔️✔️:       if (res->hasAtomBookmark(ci_RIGHTMOST_ATOM)) {
        // RDKit✔️✔️:         res->clearAtomBookmark(ci_RIGHTMOST_ATOM);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       molVect[0] = nullptr;  // NOTE: to avoid leaks on failures, this should
        // RDKit✔️✔️:                              // occur last in this if.
        // RDKit✔️✔️:     }
        // END RDKIT CPP FUNCTION toMol post-parse section
        self.close_mol_rings()?;
        self.check_chirality_specifications()?;
        self.set_unspecified_bond_types()?;
        self.adjust_atom_chirality_flags()?;
        Ok(())
    }

    fn close_mol_rings(&mut self) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP FUNCTION CloseMolRings
        // RDKit✔️✔️: void CloseMolRings(RWMol *mol, bool toleratePartials) {
        // RDKit✔️✔️:   //  Here's what we want to do here:
        // RDKit✔️✔️:   //    loop through the molecule's atom bookmarks
        // RDKit✔️✔️:   //    for each bookmark:
        // RDKit✔️✔️:   //       connect pairs of atoms sharing that bookmark
        // RDKit✔️✔️:   //          left to right (in the order in which they were
        // RDKit✔️✔️:   //          inserted into the molecule).
        // RDKit✔️✔️:   //       whilst doing this, we have to be cognizant of the fact that
        // RDKit✔️✔️:   //          there may well be partial bonds in the molecule which need
        // RDKit✔️✔️:   //          to be tied in as well.  WOO HOO! IT'S A BIG MESS!
        // RDKit✔️✔️:   PRECONDITION(mol, "no molecule");
        // RDKit✔️✔️: };
        // END RDKIT CPP FUNCTION CloseMolRings
        // BEGIN RDKIT CPP FUNCTION CloseMolRings matched-closure realization
        // RDKit❗✔️: while (bookmarkIt != mol->getAtomBookmarks()->end()) {
        // RDKit❗✔️:   ...
        // RDKit❗✔️:   // connect pairs of atoms sharing that bookmark
        // RDKit❗✔️:   // left to right (in the order in which they were inserted)
        // RDKit❗✔️:   RWMol::BOND_PTR_LIST bonds = mol->getAllBondsWithBookmark(bookmark.first);
        // RDKit❗✔️:   ...
        // RDKit❗✔️:   bondIdx = mol->addBond(matchedBond, true);
        // RDKit❗✔️:   ...
        // RDKit❗✔️:   *closurePos = bondIdx - 1;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION CloseMolRings matched-closure realization
        let mut pending_ring_closures = std::mem::take(&mut self.pending_ring_closures);
        pending_ring_closures.sort_by_key(|closure| closure.ring_number);
        for closure in pending_ring_closures {
            let bond_id = self.close_ring_opening(
                closure.opening_atom,
                closure.closing_atom,
                closure.opening_pending_bond,
                closure.closing_pending_bond,
            )?;
            if let Some(records) = self.ring_closures_by_atom.get_mut(&closure.opening_atom)
                && let Some(record) = records.iter_mut().find(|record| {
                    record.ring_number == closure.ring_number && record.bond.is_none()
                })
            {
                record.bond = Some(bond_id);
            }
            if let Some(records) = self.ring_closures_by_atom.get_mut(&closure.closing_atom)
                && let Some(record) = records.iter_mut().find(|record| {
                    record.ring_number == closure.ring_number && record.bond.is_none()
                })
            {
                record.bond = Some(bond_id);
            }
        }
        if self.ring_openings.is_empty() {
            return Ok(());
        }
        Err(SmilesParseError::ParseError("unclosed ring".to_string()))
    }

    fn check_chirality_specifications(&mut self) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP FUNCTION CheckChiralitySpecifications
        // RDKit✔️✔️: static const std::map<int, int> permutationLimits = {
        // RDKit✔️✔️:     {RDKit::Atom::ChiralType::CHI_TETRAHEDRAL, 2},
        // RDKit✔️✔️:     {RDKit::Atom::ChiralType::CHI_ALLENE, 2},
        // RDKit✔️✔️:     {RDKit::Atom::ChiralType::CHI_SQUAREPLANAR, 3},
        // RDKit✔️✔️:     {RDKit::Atom::ChiralType::CHI_OCTAHEDRAL, 30},
        // RDKit✔️✔️:     {RDKit::Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL, 20}};
        // RDKit✔️✔️: bool checkChiralPermutation(int chiralTag, int permutation) {
        // RDKit✔️✔️:   if (chiralTag > RDKit::Atom::ChiralType::CHI_OTHER &&
        // RDKit✔️✔️:       permutationLimits.find(chiralTag) != permutationLimits.end() &&
        // RDKit✔️✔️:       (permutation < 0 || permutation > permutationLimits.at(chiralTag))) {
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return true;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: void CheckChiralitySpecifications(RDKit::RWMol *mol, bool strict) {
        // RDKit✔️✔️:   PRECONDITION(mol, "no molecule");
        // RDKit✔️✔️:   for (const auto atom : mol->atoms()) {
        // RDKit✔️✔️:     int permutation;
        // RDKit✔️✔️:     if (atom->getChiralTag() > RDKit::Atom::ChiralType::CHI_OTHER &&
        // RDKit✔️✔️:         permutationLimits.find(atom->getChiralTag()) !=
        // RDKit✔️✔️:             permutationLimits.end() &&
        // RDKit✔️✔️:         atom->getPropIfPresent(common_properties::_chiralPermutation,
        // RDKit✔️✔️:                                permutation)) {
        // RDKit✔️✔️:       if (!checkChiralPermutation(atom->getChiralTag(), permutation)) {
        // RDKit✔️✔️:         throw SmilesParseException(error);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       if (atom->getChiralTag() == RDKit::Atom::ChiralType::CHI_TETRAHEDRAL) {
        // RDKit✔️✔️:         if (permutation == 0 || permutation == 1) {
        // RDKit✔️✔️:           atom->setChiralTag(RDKit::Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
        // RDKit✔️✔️:           atom->clearProp(common_properties::_chiralPermutation);
        // RDKit✔️✔️:         } else if (permutation == 2) {
        // RDKit✔️✔️:           atom->setChiralTag(RDKit::Atom::ChiralType::CHI_TETRAHEDRAL_CW);
        // RDKit✔️✔️:           atom->clearProp(common_properties::_chiralPermutation);
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION CheckChiralitySpecifications
        for atom in self.builder.atoms_mut() {
            let Some(permutation) = atom.chiral_permutation() else {
                continue;
            };
            // RDKit stores the parsed permutation as a signed int. COSMolKit
            // keeps the modeled state as `u32`, so the negative branch of the
            // source check is unrepresentable here and the upper-bound check is
            // the only reachable validation path.
            if let Some(limit) = chiral_permutation_limit(atom.chiral_tag())
                && permutation > u32::try_from(limit).expect("chiral permutation limit is positive")
            {
                return Err(SmilesParseError::ParseError(format!(
                    "Invalid chiral specification on atom {}",
                    atom.id().index()
                )));
            }
            if atom.chiral_tag() == ChiralTag::Tetrahedral {
                if permutation == 0 || permutation == 1 {
                    atom.set_chiral_tag(ChiralTag::TetrahedralCcw);
                    atom.set_chiral_permutation(None);
                } else if permutation == 2 {
                    atom.set_chiral_tag(ChiralTag::TetrahedralCw);
                    atom.set_chiral_permutation(None);
                }
            }
        }
        Ok(())
    }

    fn set_unspecified_bond_types(&mut self) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP FUNCTION SetUnspecifiedBondTypes
        // RDKit✔️✔️: void SetUnspecifiedBondTypes(RWMol *mol) {
        // RDKit✔️✔️:   PRECONDITION(mol, "no molecule");
        // RDKit✔️✔️:   for (auto bond : mol->bonds()) {
        // RDKit✔️✔️:     if (bond->hasProp(RDKit::common_properties::_unspecifiedOrder)) {
        // RDKit✔️✔️:       bond->setBondType(GetUnspecifiedBondType(mol, bond->getBeginAtom(),
        // RDKit✔️✔️:                                                bond->getEndAtom()));
        // RDKit✔️✔️:       if (bond->getBondType() == Bond::AROMATIC) {
        // RDKit✔️✔️:         bond->setIsAromatic(true);
        // RDKit✔️✔️:       } else {
        // RDKit✔️✔️:         bond->setIsAromatic(false);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION SetUnspecifiedBondTypes
        for bond_id in self.explicit_unspecified_bonds.iter().copied() {
            let (begin, end) = {
                let bond = self.builder.bond(bond_id).ok_or_else(|| {
                    SmilesParseError::ParseError(format!("missing bond {bond_id}"))
                })?;
                (bond.begin(), bond.end())
            };
            let begin_aromatic = self
                .atom_aromatic
                .get(begin.index())
                .copied()
                .unwrap_or(false);
            let end_aromatic = self
                .atom_aromatic
                .get(end.index())
                .copied()
                .unwrap_or(false);
            let order = get_unspecified_bond_type_for_atoms(begin_aromatic, end_aromatic);
            let bond = self
                .builder
                .bond_mut(bond_id)
                .ok_or_else(|| SmilesParseError::ParseError(format!("missing bond {bond_id}")))?;
            bond.set_order(order);
            bond.set_aromatic(order == BondOrder::Aromatic);
        }
        Ok(())
    }

    fn adjust_atom_chirality_flags(&mut self) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP FUNCTION AdjustAtomChiralityFlags
        // RDKit✔️✔️: void AdjustAtomChiralityFlags(RWMol *mol) {
        // RDKit✔️✔️:   PRECONDITION(mol, "no molecule");
        // RDKit✔️✔️:   for (auto atom : mol->atoms()) {
        // RDKit✔️✔️:     Atom::ChiralType chiralType = atom->getChiralTag();
        // RDKit✔️✔️:     if (chiralType == Atom::CHI_TETRAHEDRAL_CW ||
        // RDKit✔️✔️:         chiralType == Atom::CHI_TETRAHEDRAL_CCW) {
        // RDKit✔️✔️:       INT_LIST bondOrdering;
        // RDKit✔️✔️:       unsigned int numClosures = GetBondOrdering(bondOrdering, mol, atom);
        // RDKit✔️✔️:
        // RDKit✔️✔️:       // ok, we now have the SMILES ordering of the bonds, figure out the
        // RDKit✔️✔️:       // permutation order.
        // RDKit✔️✔️:       //
        // RDKit✔️✔️:       //  This whole thing is necessary because the ring-closure bonds
        // RDKit✔️✔️:       //  in the SMILES come before the bonds to the other neighbors, but
        // RDKit✔️✔️:       //  they come after the neighbors in the molecule we build.
        // RDKit✔️✔️:       //  A crude example:
        // RDKit✔️✔️:       //   in F[C@](Cl)(Br)I the C-Cl bond is index 1 in both SMILES
        // RDKit✔️✔️:       //         and as built
        // RDKit✔️✔️:       //   in F[C@]1(Br)I.Cl1 the C-Cl bond is index 1 in the SMILES
        // RDKit✔️✔️:       //         and index 3 as built.
        // RDKit✔️✔️:       //
        // RDKit✔️✔️:       int nSwaps = atom->getPerturbationOrder(bondOrdering);
        // RDKit✔️✔️:       // FIX: explain this one:
        // RDKit✔️✔️:       // At least part of what's going on here for degree 3 atoms:
        // RDKit✔️✔️:       //   - The first part: if we're at the beginning of the SMILES and have
        // RDKit✔️✔️:       //      an explicit H, we need to add a swap.
        // RDKit✔️✔️:       //      This is to reflect that [C@](Cl)(F)C is equivalent to Cl[C@@](F)C
        // RDKit✔️✔️:       //      but [C@H](Cl)(F)C is fine as-is (The H-C bond is the one you look
        // RDKit✔️✔️:       //      down).
        // RDKit✔️✔️:       //   - The second part is more complicated and deals with situations like
        // RDKit✔️✔️:       //      F[C@]1CCO1. In this case we otherwise end up looking like we need
        // RDKit✔️✔️:       //      to invert the chirality, which is bogus. The chirality here needs
        // RDKit✔️✔️:       //      to remain @ just as it does in F[C@](Cl)CCO1
        // RDKit✔️✔️:       //   - We have to be careful with the second part to not sweep things like
        // RDKit✔️✔️:       //      C[S@]2(=O).Cl2 into the same bin (was github #760). We detect
        // RDKit✔️✔️:       //      those cases by looking for unsaturated atoms
        // RDKit✔️✔️:       //
        // RDKit✔️✔️:       if (Canon::chiralAtomNeedsTagInversion(
        // RDKit✔️✔️:               *mol, atom, atom->hasProp(common_properties::_SmilesStart),
        // RDKit✔️✔️:               numClosures)) {
        // RDKit✔️✔️:         ++nSwaps;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       // std::cerr << "nswaps " << atom->getIdx() << " " << nSwaps
        // RDKit✔️✔️:       //           << std::endl;
        // RDKit✔️✔️:       // std::copy(bondOrdering.begin(), bondOrdering.end(),
        // RDKit✔️✔️:       //           std::ostream_iterator<int>(std::cerr, ", "));
        // RDKit✔️✔️:       // std::cerr << std::endl;
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
        let atom_ids: Vec<AtomId> = self.builder.atoms().iter().map(|atom| atom.id()).collect();
        let mut atoms_to_invert = Vec::new();
        let mut nontetra_permutations = Vec::new();
        for atom_id in atom_ids {
            let atom =
                self.builder.atoms().get(atom_id.index()).ok_or_else(|| {
                    SmilesParseError::ParseError(format!("missing atom {atom_id}"))
                })?;
            match atom.chiral_tag() {
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw => {
                    let (bond_ordering, num_closures) = self.get_bond_ordering(atom_id)?;
                    let mut swaps = self.perturbation_order(atom_id, &bond_ordering)?;
                    if self.chiral_atom_needs_tag_inversion(atom_id, num_closures)? {
                        swaps = swaps.saturating_add(1);
                    }
                    if swaps % 2 == 1 {
                        atoms_to_invert.push(atom_id);
                    }
                }
                ChiralTag::SquarePlanar
                | ChiralTag::TrigonalBipyramidal
                | ChiralTag::Octahedral => {
                    let (mut bond_ordering, _) = self.get_bond_ordering(atom_id)?;
                    let Some(ref_max) = nontetrahedral_max_neighbors(atom.chiral_tag()) else {
                        continue;
                    };
                    if bond_ordering.len() < ref_max {
                        let missing_count = ref_max - bond_ordering.len();
                        let insert_at = if self.smiles_start_atoms.contains(&atom_id) {
                            0
                        } else {
                            1.min(bond_ordering.len())
                        };
                        bond_ordering.splice(insert_at..insert_at, vec![-1; missing_count]);
                    }
                    let permutation =
                        self.nontetrahedral_chiral_permutation(atom_id, &bond_ordering, true)?;
                    nontetra_permutations.push((atom_id, permutation));
                }
                _ => {}
            }
        }
        for atom_id in atoms_to_invert {
            let atom = self
                .builder
                .atom_mut(atom_id)
                .ok_or_else(|| SmilesParseError::ParseError(format!("missing atom {atom_id}")))?;
            match atom.chiral_tag() {
                ChiralTag::TetrahedralCw => atom.set_chiral_tag(ChiralTag::TetrahedralCcw),
                ChiralTag::TetrahedralCcw => atom.set_chiral_tag(ChiralTag::TetrahedralCw),
                _ => {}
            }
        }
        for (atom_id, permutation) in nontetra_permutations {
            self.builder
                .atom_mut(atom_id)
                .ok_or_else(|| SmilesParseError::ParseError(format!("missing atom {atom_id}")))?
                .set_chiral_permutation(Some(permutation));
        }
        Ok(())
    }

    fn get_bond_ordering(&self, atom_id: AtomId) -> Result<(Vec<i32>, usize), SmilesParseError> {
        // BEGIN RDKIT CPP FUNCTION GetBondOrdering
        // RDKit✔️✔️: unsigned int GetBondOrdering(INT_LIST &bondOrdering, const RDKit::RWMol *mol,
        // RDKit✔️✔️:                              const RDKit::Atom *atom) {
        // RDKit✔️✔️:   PRECONDITION(mol, "no mol");
        // RDKit✔️✔️:   PRECONDITION(atom, "no atom");
        // RDKit✔️✔️:   INT_VECT ringClosures;
        // RDKit✔️✔️:   atom->getPropIfPresent(common_properties::_RingClosures, ringClosures);
        // RDKit✔️✔️:   std::list<SIZET_PAIR> neighbors;
        // RDKit✔️✔️:   neighbors.emplace_back(atom->getIdx(), -1);
        // RDKit✔️✔️:   std::list<size_t> bondOrder;
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
        let ring_closures: Vec<i32> = self
            .ring_closures_by_atom
            .get(&atom_id)
            .into_iter()
            .flatten()
            .filter_map(|record| record.bond.map(|bond| bond.index() as i32))
            .collect();
        let mut neighbors = Vec::<(usize, i32)>::new();
        neighbors.push((atom_id.index(), -1));
        for &bond_id in self.builder.neighbor_bonds(atom_id) {
            let bond = self
                .builder
                .bond(bond_id)
                .ok_or_else(|| SmilesParseError::ParseError(format!("missing bond {bond_id}")))?;
            let neighbor = if bond.begin() == atom_id {
                bond.end()
            } else {
                bond.begin()
            };
            let bond_idx = bond_id.index() as i32;
            if !ring_closures.contains(&bond_idx) {
                neighbors.push((neighbor.index(), bond_idx));
            }
        }
        neighbors.sort_by_key(|(neighbor, _)| *neighbor);
        let self_pos = neighbors
            .iter()
            .position(|(neighbor, _)| *neighbor == atom_id.index())
            .ok_or_else(|| SmilesParseError::ParseError("weird atom ordering".to_string()))?;
        let mut bond_ordering = Vec::new();
        for (idx, (_, bond_idx)) in neighbors.into_iter().enumerate() {
            if idx == self_pos {
                bond_ordering.extend(ring_closures.iter().copied());
            } else {
                bond_ordering.push(bond_idx);
            }
        }
        Ok((bond_ordering, ring_closures.len()))
    }

    fn perturbation_order(
        &self,
        atom_id: AtomId,
        probe: &[i32],
    ) -> Result<usize, SmilesParseError> {
        // BEGIN RDKIT CPP FUNCTION Atom::getPerturbationOrder / countSwapsToInterconvert
        // RDKit✔️✔️: int Atom::getPerturbationOrder(const INT_LIST &probe) const {
        // RDKit✔️✔️:   INT_LIST ref;
        // RDKit✔️✔️:   for (const auto bnd : getOwningMol().atomBonds(this)) {
        // RDKit✔️✔️:     ref.push_back(bnd->getIdx());
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return static_cast<int>(countSwapsToInterconvert(probe, ref));
        // RDKit✔️✔️: }
        // RDKit✔️✔️: unsigned int countSwapsToInterconvert(const T &ref, T probe) {
        // RDKit✔️✔️:   PRECONDITION(ref.size() == probe.size(), "size mismatch");
        // RDKit✔️✔️:   typename T::const_iterator refIt = ref.begin();
        // RDKit✔️✔️:   typename T::iterator probeIt = probe.begin();
        // RDKit✔️✔️:   typename T::iterator probeIt2;
        // RDKit✔️✔️:   unsigned int nSwaps = 0;
        // RDKit✔️✔️:   while (refIt != ref.end()) {
        // RDKit✔️✔️:     if ((*probeIt) != (*refIt)) {
        // RDKit✔️✔️:       bool foundIt = false;
        // RDKit✔️✔️:       probeIt2 = probeIt;
        // RDKit✔️✔️:       while ((*probeIt2) != (*refIt) && probeIt2 != probe.end()) {
        // RDKit✔️✔️:         ++probeIt2;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       if (probeIt2 != probe.end()) {
        // RDKit✔️✔️:         foundIt = true;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       CHECK_INVARIANT(foundIt, "could not find probe element");
        // RDKit✔️✔️:       std::swap(*probeIt, *probeIt2);
        // RDKit✔️✔️:       nSwaps++;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     ++probeIt;
        // RDKit✔️✔️:     ++refIt;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return nSwaps;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION Atom::getPerturbationOrder / countSwapsToInterconvert
        let mut storage_order = self
            .builder
            .neighbor_bonds(atom_id)
            .iter()
            .copied()
            .map(|bond_id| bond_id.index() as i32)
            .collect::<Vec<_>>();
        if probe.len() != storage_order.len() {
            return Err(SmilesParseError::ParseError("size mismatch".to_string()));
        }
        let mut swaps = 0usize;
        for (idx, expected) in probe.iter().copied().enumerate() {
            if storage_order[idx] == expected {
                continue;
            }
            let Some(found_idx) = storage_order[idx..]
                .iter()
                .position(|value| *value == expected)
                .map(|offset| idx + offset)
            else {
                return Err(SmilesParseError::ParseError(
                    "could not find probe element".to_string(),
                ));
            };
            storage_order.swap(idx, found_idx);
            swaps = swaps.saturating_add(1);
        }
        Ok(swaps)
    }

    fn chiral_atom_needs_tag_inversion(
        &self,
        atom_id: AtomId,
        num_closures: usize,
    ) -> Result<bool, SmilesParseError> {
        // BEGIN RDKIT CPP FUNCTION Canon::chiralAtomNeedsTagInversion
        // RDKit✔️✔️: bool chiralAtomNeedsTagInversion(const RDKit::ROMol &mol,
        // RDKit✔️✔️:                                  const RDKit::Atom *atom, bool isAtomFirst,
        // RDKit✔️✔️:                                  size_t numClosures) {
        // RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
        // RDKit✔️✔️:   return atom->getDegree() == 3 &&
        // RDKit✔️✔️:          ((isAtomFirst && atom->getNumExplicitHs() == 1) ||
        // RDKit✔️✔️:           (!details::atomHasFourthValence(atom) && numClosures == 1 &&
        // RDKit✔️✔️:            !details::isUnsaturated(atom, mol)));
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION Canon::chiralAtomNeedsTagInversion
        let atom = self
            .builder
            .atoms()
            .get(atom_id.index())
            .ok_or_else(|| SmilesParseError::ParseError(format!("missing atom {atom_id}")))?;
        let degree = self.atom_degree(atom_id);
        let is_atom_first = self.smiles_start_atoms.contains(&atom_id);
        Ok(degree == 3
            && ((is_atom_first && atom.explicit_hydrogens() == 1)
                || (!self.atom_has_fourth_valence(atom_id)?
                    && num_closures == 1
                    && !self.is_unsaturated(atom_id))))
    }

    fn atom_has_fourth_valence(&self, atom_id: AtomId) -> Result<bool, SmilesParseError> {
        // BEGIN RDKIT CPP FUNCTION Canon::details::atomHasFourthValence
        // RDKit✔️✔️: bool atomHasFourthValence(const Atom *atom) {
        // RDKit✔️✔️:   if (atom->getNumExplicitHs() == 1 ||
        // RDKit✔️✔️:       (!atom->needsUpdatePropertyCache() &&
        // RDKit✔️✔️:        atom->getValence(Atom::ValenceType::IMPLICIT) == 1)) {
        // RDKit✔️✔️:     return true;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (atom->hasQuery()) {
        // RDKit✔️✔️:     return hasSingleHQuery(atom->getQuery());
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return false;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION Canon::details::atomHasFourthValence
        //
        // NOTE: COSMolKit does not implement per-atom property cache
        // (needsUpdatePropertyCache / getValence(IMPLICIT)). During SMILES
        // canonicalization the explicit H count is already set by the
        // property-cache equivalent, so the implicit-valence branch collapses
        // to the explicit-H state already stored on the atom.
        let atom = self
            .builder
            .atoms()
            .get(atom_id.index())
            .ok_or_else(|| SmilesParseError::ParseError(format!("missing atom {atom_id}")))?;
        Ok(atom.explicit_hydrogens() == 1
            || atom.query().is_some_and(atom_query_has_single_h_count))
    }

    fn is_unsaturated(&self, atom_id: AtomId) -> bool {
        // BEGIN RDKIT CPP FUNCTION Canon::details::isUnsaturated
        // RDKit✔️✔️: bool isUnsaturated(const Atom *atom, const ROMol &mol) {
        // RDKit✔️✔️:   for (auto bond : mol.atomBonds(atom)) {
        // RDKit✔️✔️:     // can't just check for single bonds, because dative bonds also have an
        // RDKit✔️✔️:     // order of 1
        // RDKit✔️✔️:     if (bond->getBondTypeAsDouble() > 1) {
        // RDKit✔️✔️:       return true;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return false;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION Canon::details::isUnsaturated
        self.builder
            .neighbor_bonds(atom_id)
            .iter()
            .copied()
            .any(|bond_id| {
                self.builder
                    .bond(bond_id)
                    .is_some_and(|bond| bond_order_as_double(bond.order()) > 1.0)
            })
    }

    fn nontetrahedral_chiral_permutation(
        &self,
        atom_id: AtomId,
        probe: &[i32],
        inverse: bool,
    ) -> Result<u32, SmilesParseError> {
        // BEGIN RDKIT CPP FUNCTION Chirality::getChiralPermutation
        // RDKit✔️✔️: unsigned int getChiralPermutation(const Atom *cen, const INT_LIST &probe,
        // RDKit✔️✔️:                                   bool inverse) {
        // RDKit✔️✔️:   int perm;
        // RDKit✔️✔️:   if (!cen->getPropIfPresent(common_properties::_chiralPermutation, perm) ||
        // RDKit✔️✔️:       perm <= 0) {
        // RDKit✔️✔️:     return 0;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   decltype(&swap_octahedral) swap_func = nullptr;
        // RDKit✔️✔️:   switch (cen->getChiralTag()) {
        // RDKit✔️✔️:     case Atom::ChiralType::CHI_OCTAHEDRAL:
        // RDKit✔️✔️:       if (probe.size() > 6) {
        // RDKit✔️✔️:         return 0;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       swap_func = swap_octahedral;
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     case Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
        // RDKit✔️✔️:       if (probe.size() > 5) {
        // RDKit✔️✔️:         return 0;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       swap_func = swap_trigonalbipyramidal;
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     case Atom::ChiralType::CHI_SQUAREPLANAR:
        // RDKit✔️✔️:       if (probe.size() > 4) {
        // RDKit✔️✔️:         return 0;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       swap_func = swap_squareplanar;
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     default:
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   unsigned int nbrIdx = 0;
        // RDKit✔️✔️:   std::vector<int> order(cen->getOwningMol().getNumBonds(), -1);
        // RDKit✔️✔️:   for (const auto bnd : cen->getOwningMol().atomBonds(cen)) {
        // RDKit✔️✔️:     order[bnd->getIdx()] = nbrIdx++;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   std::vector<unsigned int> nbrPerm(nbrIdx);
        // RDKit✔️✔️:   std::iota(nbrPerm.begin(), nbrPerm.end(), 0);
        // RDKit✔️✔️:   std::vector<unsigned int> probePerm(probe.size());
        // RDKit✔️✔️:   nbrIdx = 0;
        // RDKit✔️✔️:   for (auto v : probe) {
        // RDKit✔️✔️:     probePerm[nbrIdx++] = v < 0 ? -1 : order[v];
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (nbrPerm.size() < nbrIdx) {
        // RDKit✔️✔️:     nbrPerm.insert(nbrPerm.end(), nbrIdx - nbrPerm.size(), -1);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (inverse) {
        // RDKit✔️✔️:     std::swap(nbrPerm, probePerm);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   for (unsigned int i = 0; i < probePerm.size() - 1; ++i) {
        // RDKit✔️✔️:     auto pval = probePerm[i];
        // RDKit✔️✔️:     if (nbrPerm[i] == pval) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     auto tgt = std::find(nbrPerm.begin() + i, nbrPerm.end(), pval);
        // RDKit✔️✔️:     perm = swap_func(perm, i, tgt - nbrPerm.begin());
        // RDKit✔️✔️:     std::swap(*tgt, nbrPerm[i]);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return perm;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION Chirality::getChiralPermutation
        let atom = self
            .builder
            .atoms()
            .get(atom_id.index())
            .ok_or_else(|| SmilesParseError::ParseError(format!("missing atom {atom_id}")))?;
        let mut perm = atom.chiral_permutation().unwrap_or(0);
        if perm == 0 {
            return Ok(0);
        }
        let swap_func: fn(u32, usize, usize) -> u32 = match atom.chiral_tag() {
            ChiralTag::SquarePlanar => {
                if probe.len() > 4 {
                    return Ok(0);
                }
                swap_squareplanar
            }
            ChiralTag::TrigonalBipyramidal => {
                if probe.len() > 5 {
                    return Ok(0);
                }
                swap_trigonalbipyramidal
            }
            ChiralTag::Octahedral => {
                if probe.len() > 6 {
                    return Ok(0);
                }
                swap_octahedral
            }
            _ => return Ok(0),
        };
        let mut order = vec![-1isize; self.builder.bonds().len()];
        let mut nbr_idx = 0isize;
        for &bond_id in self.builder.neighbor_bonds(atom_id) {
            order[bond_id.index()] = nbr_idx;
            nbr_idx += 1;
        }
        let mut nbr_perm = (0..nbr_idx).collect::<Vec<_>>();
        let mut probe_perm = Vec::with_capacity(probe.len());
        for value in probe {
            if *value < 0 {
                probe_perm.push(-1);
            } else {
                probe_perm.push(order[*value as usize]);
            }
        }
        if nbr_perm.len() < probe_perm.len() {
            nbr_perm.extend(std::iter::repeat_n(-1, probe_perm.len() - nbr_perm.len()));
        }
        if nbr_perm.len() != probe_perm.len() {
            return Err(SmilesParseError::ParseError(
                "probe vector size does not match".to_string(),
            ));
        }
        if inverse {
            std::mem::swap(&mut nbr_perm, &mut probe_perm);
        }
        for i in 0..probe_perm.len().saturating_sub(1) {
            let pval = probe_perm[i];
            if nbr_perm[i] == pval {
                continue;
            }
            let Some(target) = nbr_perm[i..]
                .iter()
                .position(|value| *value == pval)
                .map(|offset| i + offset)
            else {
                return Err(SmilesParseError::ParseError(
                    "could not find probe element".to_string(),
                ));
            };
            perm = swap_func(perm, i, target);
            nbr_perm.swap(i, target);
        }
        Ok(perm)
    }

    fn is_empty(&self) -> bool {
        self.active_atom.is_none()
            && self.atom_aromatic.is_empty()
            && self.bond_pairs.is_empty()
            && self.explicit_unspecified_bonds.is_empty()
            && self.branch_stack.is_empty()
            && self.ring_openings.is_empty()
            && self.ring_closures_by_atom.is_empty()
            && self.pending_ring_closures.is_empty()
            && self.smiles_start_atoms.is_empty()
            && self.cx_bond_order.is_empty()
            && self.next_cx_smiles_bond_idx == 0
            && self.temporary_chiral_permutations.is_empty()
    }

    fn append_atom(&mut self, atom: SmilesAtomToken) -> AtomId {
        let is_aromatic = atom.spec.is_aromatic();
        let atom_id = self.builder.add_atom(atom.spec);
        self.atom_aromatic.push(is_aromatic);
        atom_id
    }

    fn add_bond_with_cx_index(
        &mut self,
        begin: AtomId,
        end: AtomId,
        order: BondOrder,
    ) -> Result<BondId, SmilesParseError> {
        let bond_idx = self.next_cx_smiles_bond_idx.to_string();
        self.next_cx_smiles_bond_idx += 1;
        let spec = BondSpec::new(begin, end, order)
            .with_aromatic(order == BondOrder::Aromatic)
            .with_prop(CXSMILES_BOND_IDX_PROP, bond_idx);
        let bond_id = self
            .builder
            .add_bond(spec)
            .map_err(|error| SmilesParseError::ParseError(error.to_string()))?;
        self.bond_pairs.insert(canonical_bond_pair(begin, end));
        Ok(bond_id)
    }

    fn add_bond_from_token(
        &mut self,
        atom_idx1: AtomId,
        atom_idx2: AtomId,
        bond: SmilesBondToken,
    ) -> Result<BondId, SmilesParseError> {
        self.add_bond_from_token_with_cx_idx(atom_idx1, atom_idx2, bond, None)
    }

    fn add_bond_from_token_with_cx_idx(
        &mut self,
        atom_idx1: AtomId,
        atom_idx2: AtomId,
        bond: SmilesBondToken,
        cx_smiles_bond_idx: Option<usize>,
    ) -> Result<BondId, SmilesParseError> {
        let (begin, end, order) = match bond.order {
            BondOrder::DativeRight => (atom_idx1, atom_idx2, BondOrder::Dative),
            BondOrder::DativeLeft => (atom_idx2, atom_idx1, BondOrder::Dative),
            other => (atom_idx1, atom_idx2, other),
        };
        let cx_smiles_bond_idx = cx_smiles_bond_idx.unwrap_or_else(|| {
            let idx = self.next_cx_smiles_bond_idx;
            self.next_cx_smiles_bond_idx += 1;
            idx
        });
        let mut spec = BondSpec::new(begin, end, order)
            .with_aromatic(bond.is_aromatic || order == BondOrder::Aromatic)
            .with_direction(bond.direction)
            .with_prop(CXSMILES_BOND_IDX_PROP, cx_smiles_bond_idx.to_string());
        if bond.explicit_unspecified_order {
            spec = spec.with_prop(UNSPECIFIED_ORDER_PROP, "1");
        }
        if bond.is_null_query {
            spec = spec.with_query(QueryNode::predicate(BondQueryPredicate::Any));
        }
        let bond_id = self
            .builder
            .add_bond(spec)
            .map_err(|error| SmilesParseError::ParseError(error.to_string()))?;
        self.bond_pairs.insert(canonical_bond_pair(begin, end));
        if bond.explicit_unspecified_order {
            self.explicit_unspecified_bonds.push(bond_id);
        }
        Ok(bond_id)
    }

    fn into_molecule(self) -> Result<Molecule, SmilesParseError> {
        self.builder
            .build()
            .map_err(|error| SmilesParseError::ParseError(error.to_string()))
    }

    fn cleanup_after_parsing(&mut self) {
        // BEGIN RDKIT CPP FUNCTION CleanupAfterParsing
        // RDKit✔️✔️: void CleanupAfterParsing(RWMol *mol) {
        // RDKit✔️✔️:   PRECONDITION(mol, "no molecule");
        // RDKit✔️✔️:   for (auto atom : mol->atoms()) {
        // RDKit✔️✔️:     atom->clearProp(common_properties::_RingClosures);
        // RDKit✔️✔️:     atom->clearProp(common_properties::_SmilesStart);
        // RDKit✔️✔️:     std::string label;
        // RDKit✔️✔️:     if (atom->getAtomicNum() == 0 &&
        // RDKit✔️✔️:         atom->getPropIfPresent(common_properties::atomLabel, label)) {
        // RDKit✔️✔️:       // marvinsketch can output higher labels than _AP1 and _AP2, but they
        // RDKit✔️✔️:       // aren't part of the MOL file spec so we don't treat them as attachment
        // RDKit✔️✔️:       // points
        // RDKit✔️✔️:       if (label == "_AP1") {
        // RDKit✔️✔️:         atom->setProp(common_properties::_fromAttachPoint, 1);
        // RDKit✔️✔️:       } else if (label == "_AP2") {
        // RDKit✔️✔️:         atom->setProp(common_properties::_fromAttachPoint, 2);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   for (auto bond : mol->bonds()) {
        // RDKit✔️✔️:     bond->clearProp(common_properties::_unspecifiedOrder);
        // RDKit✔️✔️:     bond->clearProp("_cxsmilesBondIdx");
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   for (auto sg : RDKit::getSubstanceGroups(*mol)) {
        // RDKit✔️✔️:     sg.clearProp("_cxsmilesindex");
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (!Chirality::getAllowNontetrahedralChirality()) {
        // END RDKIT CPP FUNCTION CleanupAfterParsing
        for atom in self.builder.atoms_mut() {
            atom.clear_prop("_RingClosures");
            atom.clear_prop(SMILES_START_PROP);
            if atom.atomic_number() == 0 {
                match atom.prop("atomLabel") {
                    Some("_AP1") => atom.set_prop("_fromAttachPoint", "1"),
                    Some("_AP2") => atom.set_prop("_fromAttachPoint", "2"),
                    _ => {}
                }
            }
        }
        for bond in self.builder.bonds_mut() {
            bond.clear_prop(UNSPECIFIED_ORDER_PROP);
            bond.clear_prop(CXSMILES_BOND_IDX_PROP);
        }
        // substance group cleanup
        for i in 0..self.builder.substance_groups_len() {
            if let Some(sg) = self.builder.substance_group_mut(i) {
                sg.clear_prop("_cxsmilesindex")
            }
        }
        // RDKit only strips non-tetrahedral chirality when the global
        // allow-nontetrahedral switch is disabled. COSMolKit currently models
        // the default RDKit behavior with that switch enabled.
        if !allow_nontetrahedral_chirality() {
            for atom in self.builder.atoms_mut() {
                if atom.chiral_permutation().is_some() {
                    atom.set_chiral_permutation(None);
                }
                match atom.chiral_tag() {
                    ChiralTag::Allene
                    | ChiralTag::SquarePlanar
                    | ChiralTag::TrigonalBipyramidal
                    | ChiralTag::Octahedral
                    | ChiralTag::Other => {
                        atom.set_chiral_tag(ChiralTag::Unspecified);
                    }
                    _ => {}
                }
            }
        }
    }

    fn set_name(&mut self, name: &str) {
        self.builder = std::mem::take(&mut self.builder).with_name(name);
    }

    fn set_property(&mut self, key: &str, value: &str) {
        self.builder = std::mem::take(&mut self.builder).with_property(key, value);
    }
}

pub(crate) fn mol_from_smiles(
    smiles: &str,
    params: &SmilesParseParams,
) -> Result<Molecule, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION MolFromSmiles
    // RDKit✔️✔️: std::unique_ptr<RWMol> MolFromSmiles(const std::string &smiles,
    // RDKit✔️✔️:                                      const SmilesParserParams &params) {
    // RDKit✔️✔️:   // Calling MolFromSmiles in a multithreaded context is generally safe *unless*
    // RDKit✔️✔️:   // the value of debugParse is different for different threads. The if
    // RDKit✔️✔️:   // statement below avoids a TSAN warning in the case where multiple threads
    // RDKit✔️✔️:   // all use the same value for debugParse.
    // RDKit✔️✔️:   if (yysmiles_debug != params.debugParse) {
    // RDKit✔️✔️:     yysmiles_debug = params.debugParse;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::string lsmiles, name, cxPart;
    // RDKit✔️✔️:   preprocessSmiles(smiles, params, lsmiles, name, cxPart);
    // RDKit✔️✔️:   // strip any leading/trailing whitespace:
    // RDKit✔️✔️:   // boost::trim_if(smi,boost::is_any_of(" \t\r\n"));
    // RDKit✔️✔️:   auto res = toMol(lsmiles, smiles_parse, lsmiles);
    // RDKit✔️✔️:   if (!res) {
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   handleCXPartAndName(res.get(), params, cxPart, name);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION MolFromSmiles
    if YYSMILES_DEBUG.load(Ordering::Relaxed) != params.debug_parse {
        YYSMILES_DEBUG.store(params.debug_parse, Ordering::Relaxed);
    }

    let preprocessed = preprocess_smiles(smiles, params)?;
    let mut state = to_mol(&preprocessed.smiles)?;
    handle_cx_part_and_name(
        &mut state,
        params,
        &preprocessed.cx_part,
        &preprocessed.name,
    )?;
    apply_smiles_postprocessing(&mut state, params)?;
    let _requires_full_sanitize_postprocess = state.should_run_full_sanitize_postprocess(params);
    let mut mol = state.into_molecule()?;
    // RDKit MolOps::sanitizeMol runs through the registered operations.
    // COSMolKit applies equivalent operations on the built Molecule.
    if params.sanitize {
        mol = mol.sanitized_with_ops(crate::SanitizeOps::ALL).map_err(
            |e: crate::OperationError| {
                SmilesParseError::ParseError(format!("sanitize during smiles parse failed: {e}"))
            },
        )?;
    }

    let (first_2d_conf_id, first_3d_conf_id) = first_2d_and_3d_conformer_ids(&mol);

    // RDKit✔️✔️: Post-sanitize: handle CX wedged-bond atom stereo detection
    // RDKit✔️✔️: BEGIN RDKIT CPP FUNCTION MolOps::assignChiralTypesFromBondDirs
    // RDKit✔️✔️: void assignChiralTypesFromBondDirs(ROMol &mol, const int confId,
    // RDKit✔️✔️:                                    const bool replaceExistingTags) {
    // RDKit✔️✔️:   if (!mol.getNumConformers()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto conf = mol.getConformer(confId);
    // RDKit✔️✔️:   boost::dynamic_bitset<> atomsSet(mol.getNumAtoms(), 0);
    // RDKit✔️✔️:   for (auto &bond : mol.bonds()) {
    // RDKit✔️✔️:     const Bond::BondDir dir = bond->getBondDir();
    // RDKit✔️✔️:     Atom *atom = bond->getBeginAtom();
    // RDKit✔️✔️:     if (dir == Bond::UNKNOWN) {
    // RDKit✔️✔️:       if (atomsSet[atom->getIdx()] || replaceExistingTags) {
    // RDKit✔️✔️:         atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit✔️✔️:         atomsSet.set(atom->getIdx());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       if (dir == Bond::BEGINWEDGE || dir == Bond::BEGINDASH) {
    // RDKit✔️✔️:         if (atomsSet[atom->getIdx()] ||
    // RDKit✔️✔️:             (!replaceExistingTags &&
    // RDKit✔️✔️:              atom->getChiralTag() != Atom::CHI_UNSPECIFIED)) {
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (atom->needsUpdatePropertyCache()) {
    // RDKit✔️✔️:           atom->updatePropertyCache(false);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         Atom::ChiralType code =
    // RDKit✔️✔️:             Chirality::atomChiralTypeFromBondDirPseudo3D(mol, bond, &conf)
    // RDKit✔️✔️:                 .value_or(Atom::ChiralType::CHI_UNSPECIFIED);
    // RDKit✔️✔️:         if (code != Atom::ChiralType::CHI_UNSPECIFIED) {
    // RDKit✔️✔️:           atomsSet.set(atom->getIdx());
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         atom->setChiralTag(code);
    // RDKit✔️✔️:
    // RDKit✔️✔️:         // within the RD representation, if a three-coordinate atom
    // RDKit✔️✔️:         // is chiral and has an implicit H, that H needs to be made explicit:
    // RDKit✔️✔️:         if (atom->getDegree() == 3 && !atom->getNumExplicitHs() &&
    // RDKit✔️✔️:             atom->getNumImplicitHs() == 1) {
    // RDKit✔️✔️:           atom->setNumExplicitHs(1);
    // RDKit✔️✔️:           atom->updatePropertyCache();
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: END RDKIT CPP FUNCTION MolOps::assignChiralTypesFromBondDirs
    // RDKit✔️✔️: SmilesParse::MolFromSmiles caller:
    // RDKit✔️✔️: if (res->hasProp(detail::_needsDetectAtomStereo)) {
    // RDKit✔️✔️:   res->clearProp(detail::_needsDetectAtomStereo);
    // RDKit✔️✔️:   if (conf) {
    // RDKit✔️✔️:     MolOps::assignChiralTypesFromBondDirs(*res, conf->getId());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    if mol.prop("_needsDetectAtomStereo").is_some() {
        mol.properties_mut().clear_prop("_needsDetectAtomStereo");
        if let Some(conf_id) = first_2d_conf_id {
            assign_chiral_types_from_bond_dirs(&mut mol, conf_id);
        }
    }

    // RDKit✔️✔️: 3D atom stereo assignment
    // RDKit✔️✔️: if (!conf && conf3d) {
    // RDKit✔️✔️:   res->updatePropertyCache(false);
    // RDKit✔️✔️:   MolOps::assignChiralTypesFrom3D(*res, conf3d->getId(), true);
    // RDKit✔️✔️: }
    if first_2d_conf_id.is_none()
        && let Some(conf_id) = first_3d_conf_id
    {
        let _ = assign_stereochemistry_from_3d(&mut mol, conf_id);
    }

    // RDKit✔️✔️: Atropisomer detection
    // RDKit✔️✔️: if (conf) {
    // RDKit✔️✔️:   Atropisomers::detectAtropisomerChirality(*res, conf);
    // RDKit✔️✔️: } else if (conf3d) {
    // RDKit✔️✔️:   Atropisomers::detectAtropisomerChirality(*res, conf3d);
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   Atropisomers::detectAtropisomerChirality(*res, nullptr);
    // RDKit✔️✔️: }
    if let Some(conf_id) = first_2d_conf_id {
        let atropo = atropisomer_stereo_from_conformer(&mol, conf_id);
        if !atropo.is_empty() {
            apply_atropisomer_stereo_assignments(&mut mol, atropo);
        }
    } else if let Some(conf_id) = first_3d_conf_id {
        let atropo = atropisomer_stereo_from_conformer(&mol, conf_id);
        if !atropo.is_empty() {
            apply_atropisomer_stereo_assignments(&mut mol, atropo);
        }
    } else {
        let atropo = atropisomer_stereo_without_conformer(&mol);
        if !atropo.is_empty() {
            apply_coordinate_free_atropisomer_assignments(&mut mol, atropo);
        }
    }

    // RDKit✔️✔️: CX wiggly-bond double-bond stereo detection and final
    // RDKit✔️✔️: stereochemistry assignment only run inside the
    // RDKit✔️✔️: (params.sanitize || params.removeHs) branch.
    // RDKit✔️✔️: if (res && (params.sanitize || params.removeHs)) {
    // RDKit✔️✔️:   if (res->hasProp(detail::_needsDetectBondStereo)) {
    // RDKit✔️✔️:     if (conf || conf3d) {
    // RDKit✔️✔️:       MolOps::clearSingleBondDirFlags(*res);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     MolOps::setDoubleBondNeighborDirections(*res, conf ? conf : conf3d);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res->clearProp(detail::_needsDetectBondStereo);
    // RDKit✔️✔️:   MolOps::assignStereochemistry(*res, cleanIt, force, flagPossible);
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   MolOps::clearSingleBondDirFlags(*res, true);
    // RDKit✔️✔️: }
    if params.sanitize || params.remove_hs {
        if mol.prop("_needsDetectBondStereo").is_some() {
            if first_2d_conf_id.is_some() || first_3d_conf_id.is_some() {
                clear_single_bond_dir_flags(&mut mol, false);
            }
            let _ = set_double_bond_neighbor_directions_impl(
                &mut mol,
                first_2d_conf_id.or(first_3d_conf_id),
            );
        }
        mol.properties_mut().clear_prop("_needsDetectBondStereo");
        let _ = assign_stereochemistry_cleanup_subset(&mut mol, true);
    } else {
        // RDKit github #337 path: after atom stereo perception, wedge-style
        // single-bond directions are no longer needed, but the CX double-bond
        // stereo marker property is preserved when sanitize/removeHs are off.
        clear_single_bond_dir_flags(&mut mol, true);
    }

    // RDKit✔️✔️: _NeedsQueryScan query completion
    // RDKit✔️✔️: if (res->hasProp(common_properties::_NeedsQueryScan)) {
    // RDKit✔️✔️:   if (!params.sanitize) {
    // RDKit✔️✔️:     MolOps::fastFindRings(*res);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   QueryOps::completeMolQueries(*res, 0xDEADBEEF);
    // RDKit✔️✔️: }
    // COSMolKit: for the currently modeled SMILES CX query-scan inputs, this
    // completes both `rb:*` and `s:*` sentinel placeholders and clears
    // `_NeedsQueryScan` once no unresolved scan work remains.
    if mol.prop("_NeedsQueryScan").is_some() {
        if !params.sanitize {
            // RDKit runs fastFindRings() here before query completion. This
            // subset computes ring-bond counts directly from the graph, so no
            // extra ring-cache mutation is required for the modeled `rb:*`
            // SMILES queries.
        }
        complete_smiles_query_scan_subset(&mut mol);
    }

    Ok(mol)
}

fn to_mol(inp: &str) -> Result<SmilesBuildState, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION toMol
    // RDKit✔️✔️: std::unique_ptr<RWMol> toMol(const std::string &inp,
    // RDKit✔️✔️:                              int func(const std::string &,
    // RDKit✔️✔️:                                       std::vector<RDKit::RWMol *> &),
    // RDKit✔️✔️:                              const std::string &origInp) {
    // RDKit✔️✔️:   // empty strings produce empty molecules:
    // RDKit✔️✔️:   if (inp.empty()) {
    // RDKit✔️✔️:     return std::make_unique<RWMol>();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::unique_ptr<RWMol> res;
    // RDKit✔️✔️:   std::vector<RDKit::RWMol *> molVect;
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     func(inp, molVect);
    // RDKit✔️✔️:     if (!molVect.empty()) {
    // RDKit✔️✔️:       res.reset(molVect[0]);
    // RDKit✔️✔️:       SmilesParseOps::CloseMolRings(res.get(), false);
    // RDKit✔️✔️:       SmilesParseOps::CheckChiralitySpecifications(res.get(), true);
    // RDKit✔️✔️:       SmilesParseOps::SetUnspecifiedBondTypes(res.get());
    // RDKit✔️✔️:       SmilesParseOps::AdjustAtomChiralityFlags(res.get());
    // RDKit✔️✔️:       // No sense leaving this bookmark intact:
    // RDKit✔️✔️:       if (res->hasAtomBookmark(ci_RIGHTMOST_ATOM)) {
    // RDKit✔️✔️:         res->clearAtomBookmark(ci_RIGHTMOST_ATOM);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       molVect[0] = nullptr;  // NOTE: to avoid leaks on failures, this should
    // RDKit✔️✔️:                              // occur last in this if.
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } catch (SmilesParseException &e) {
    // RDKit✔️✔️:     std::string nm = "SMILES";
    // RDKit✔️✔️:     if (func == smarts_parse) {
    // RDKit✔️✔️:       nm = "SMARTS";
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     BOOST_LOG(rdErrorLog) << nm << " Parse Error: " << e.what()
    // RDKit✔️✔️:                           << " for input: '" << origInp << "'" << std::endl;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // reset res so that we return a nullptr. We don't want to reset(),
    // RDKit✔️✔️:     // because that would delete the mol and leak any unmatched
    // RDKit✔️✔️:     // ring closure bonds. These will be cleaned up in the loop below.
    // RDKit✔️✔️:     res.release();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto *molPtr : molVect) {
    // RDKit✔️✔️:     if (molPtr) {
    // RDKit✔️✔️:       // Clean-up the bond bookmarks when not calling CloseMolRings
    // RDKit✔️✔️:       SmilesParseOps::CleanupAfterParseError(molPtr);
    // RDKit✔️✔️:       delete molPtr;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION toMol
    if inp.is_empty() {
        return Ok(SmilesBuildState::new());
    }

    let mut state = SmilesBuildState::new();
    let result = {
        let lexer = SmilesLexer::new(inp);
        let mut parser = SmilesParser::new(lexer);
        parser.parse_mol(&mut state)?;
        state.finish_parse()?;
        Ok(())
    };
    match result {
        Ok(()) => Ok(state),
        Err(err) => {
            state.cleanup_after_parse_error();
            let _ = report_parse_error(&format!("{err} for input: '{inp}'"), false);
            Err(err)
        }
    }
}

fn preprocess_smiles(
    smiles: &str,
    params: &SmilesParseParams,
) -> Result<PreprocessedSmiles, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION preprocessSmiles
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️: void preprocessSmiles(const std::string &smiles, const T &params,
    // RDKit✔️✔️:                       std::string &lsmiles, std::string &name,
    // RDKit✔️✔️:                       std::string &cxPart) {
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
    // RDKit✔️✔️:
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
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION preprocessSmiles
    let mut lsmiles = String::new();
    let mut name = String::new();
    let mut cx_part = String::new();

    if params.parse_name && !params.allow_cxsmiles {
        if let Some(sidx) = smiles.find([' ', '\t'])
            && sidx != 0
        {
            lsmiles = smiles[..sidx].to_string();
            name = smiles[sidx..].trim().to_string();
        }
    } else if params.allow_cxsmiles
        && let Some(sidx) = smiles.find([' ', '\t'])
        && sidx != 0
    {
        lsmiles = smiles[..sidx].to_string();
        cx_part = smiles[sidx..].trim().to_string();
    }

    if lsmiles.is_empty() {
        lsmiles = smiles.to_string();
    }

    if !params.replacements.is_empty() {
        let mut smi = lsmiles;
        let mut loop_again = true;
        while loop_again {
            loop_again = false;
            for (from, to) in &params.replacements {
                if smi.contains(from) {
                    loop_again = true;
                    smi = smi.replace(from, to);
                }
            }
        }
        lsmiles = smi;
    }

    Ok(PreprocessedSmiles {
        smiles: lsmiles,
        name,
        cx_part,
    })
}

fn apply_smiles_postprocessing(
    state: &mut SmilesBuildState,
    params: &SmilesParseParams,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION MolFromSmiles post-CX processing section
    // RDKit✔️✔️:   // get a conformer
    // RDKit✔️✔️:   const Conformer *conf = nullptr, *conf3d = nullptr;
    // RDKit✔️✔️:   if (res && res->getNumConformers() > 0) {
    // RDKit✔️✔️:     for (unsigned int confId = 0; confId < res->getNumConformers(); ++confId) {
    // RDKit✔️✔️:       auto *testConf = &res->getConformer(confId);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (res && (params.sanitize || params.removeHs)) {
    // RDKit✔️✔️:     if (params.removeHs) {
    // RDKit✔️✔️:       MolOps::RemoveHsParameters rhp;
    // RDKit✔️✔️:       rhp.updateExplicitCount = true;
    // RDKit✔️✔️:       MolOps::removeHs(*res, rhp, params.sanitize);
    // RDKit✔️✔️:     } else if (params.sanitize) {
    // RDKit✔️✔️:       MolOps::sanitizeMol(*res);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION MolFromSmiles post-CX processing section
    let _ = params;
    if state.is_empty() {
        return Ok(());
    }
    if params.remove_hs {
        state.remove_hs_update_explicit_count(&RemoveHsParams {
            update_explicit_count: true,
            ..RemoveHsParams::default()
        })?;
    }
    if !params.sanitize {
        if !params.skip_cleanup {
            state.cleanup_after_parsing();
        }
        return Ok(());
    }
    // Sanitization is split across the active reader path:
    // this state-level hook performs the pre-build removeHs() branch, while
    // mol_from_smiles() performs the molecule-level sanitize operations.
    if !params.skip_cleanup {
        state.cleanup_after_parsing();
    }
    Ok(())
}

impl SmilesBuildState {
    fn remove_hs_update_explicit_count(
        &mut self,
        params: &RemoveHsParams,
    ) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP FUNCTION removeHs / molRemoveH / shouldRemoveH
        // RDKit✔️✔️: void removeHs(RWMol &mol, const RemoveHsParameters &ps, bool sanitize) {
        // RDKit✔️✔️:   if (ps.removeAndTrackIsotopes) {
        // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
        // RDKit✔️✔️:     if (shouldRemoveH(mol, atom, ps)) {
        // RDKit✔️✔️:       atomsToRemove.set(atom->getIdx());
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   for (int idx = mol.getNumAtoms() - 1; idx >= 0; --idx) {
        // RDKit✔️✔️:     if (atomsToRemove[idx]) {
        // RDKit✔️✔️:       molRemoveH(mol, idx, ps.updateExplicitCount);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // RDKit✔️✔️: void molRemoveH(RWMol &mol, unsigned int idx, bool updateExplicitCount) {
        // RDKit✔️✔️:   auto atom = mol.getAtomWithIdx(idx);
        // RDKit✔️✔️:   PRECONDITION(atom->getAtomicNum() == 1, "idx corresponds to a non-Hydrogen");
        // RDKit✔️✔️:   for (const auto bond : mol.atomBonds(atom)) {
        // RDKit✔️✔️:     Atom *heavyAtom = bond->getOtherAtom(atom);
        // RDKit✔️✔️:     if (updateExplicitCount || heavyAtom->getNoImplicit() ||
        // RDKit✔️✔️:         heavyAtom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
        // RDKit✔️✔️:       heavyAtom->setNumExplicitHs(heavyAtom->getNumExplicitHs() + 1);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       // this is a special case related to Issue 228 and the
        // RDKit✔️✔️:       // "disappearing Hydrogen" problem discussed in MolOps::adjustHs
        // RDKit✔️✔️:       //
        // RDKit✔️✔️:       // If we remove a hydrogen from an aromatic N or P, or if
        // RDKit✔️✔️:       // the heavy atom it is connected to is not in its default
        // RDKit✔️✔️:       // valence state, we need to be *sure* to increment the
        // RDKit✔️✔️:       // explicit count, even if the H itself isn't marked as explicit
        // RDKit✔️✔️:       const INT_VECT &defaultVs =
        // RDKit✔️✔️:           PeriodicTable::getTable()->getValenceList(heavyAtomNum);
        // RDKit✔️✔️:       if (((heavyAtomNum == 7 || heavyAtomNum == 15 ||
        // RDKit✔️✔️:             may_need_extra_H(mol, heavyAtom)) &&
        // RDKit✔️✔️:            isAromaticAtom(*heavyAtom)) ||
        // RDKit✔️✔️:           (std::find(defaultVs.begin() + 1, defaultVs.end(),
        // RDKit✔️✔️:                      heavyAtom->getTotalValence()) != defaultVs.end())) {
        // RDKit✔️✔️:         heavyAtom->setNumExplicitHs(heavyAtom->getNumExplicitHs() + 1);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     // One other consequence of removing the H from the graph is
        // RDKit✔️✔️:     // that we may change the ordering of the bonds about a
        // RDKit✔️✔️:     // chiral center.  This may change the chiral label at that
        // RDKit✔️✔️:     // atom.  We deal with that by explicitly checking here:
        // RDKit✔️✔️:     if (heavyAtom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
        // RDKit✔️✔️:       INT_LIST neighborIndices;
        // RDKit✔️✔️:       for (const auto &nbnd : mol.atomBonds(heavyAtom)) {
        // RDKit✔️✔️:         if (nbnd->getIdx() != bond->getIdx()) {
        // RDKit✔️✔️:           neighborIndices.push_back(nbnd->getIdx());
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       neighborIndices.push_back(bond->getIdx());
        // RDKit✔️✔️:
        // RDKit✔️✔️:       int nSwaps = heavyAtom->getPerturbationOrder(neighborIndices);
        // RDKit✔️✔️:       if (nSwaps % 2) {
        // RDKit✔️✔️:         heavyAtom->invertChirality();
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     // If we are removing a H atom that defines bond stereo (e.g. imines),
        // RDKit✔️✔️:     // Then also remove the bond stereo information, as it is no longer valid.
        // RDKit✔️✔️:     if (heavyAtom->getDegree() == 2) {
        // RDKit✔️✔️:       for (auto &nbnd : mol.atomBonds(heavyAtom)) {
        // RDKit✔️✔️:         if (nbnd != bond) {
        // RDKit✔️✔️:           if (nbnd->getStereo() > Bond::STEREOANY) {
        // RDKit✔️✔️:             nbnd->setStereo(Bond::STEREONONE);
        // RDKit✔️✔️:             nbnd->getStereoAtoms().clear();
        // RDKit✔️✔️:           }
        // RDKit✔️✔️:           break;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   mol.removeAtom(atom, clearProps);
        // RDKit✔️✔️: }
        // RDKit✔️✔️: bool shouldRemoveH(const RWMol &mol, const Atom *atom,
        // RDKit✔️✔️:                    const RemoveHsParameters &ps) {
        // RDKit✔️✔️:   if (atom->getAtomicNum() != 1) {
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (!ps.removeWithQuery && atom->hasQuery()) {
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (!ps.removeDegreeZero && !atom->getDegree()) {
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (!ps.removeHigherDegrees && atom->getDegree() > 1) {
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (!ps.removeIsotopes && !ps.removeAndTrackIsotopes && atom->getIsotope()) {
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (!ps.removeNonimplicit && !atom->hasProp(common_properties::isImplicit)) {
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (!ps.removeMapped && atom->getAtomMapNum()) {
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (!ps.removeHydrides && atom->getFormalCharge() == -1) {
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (atom->getDegree() &&
        // RDKit✔️✔️:       (!ps.removeDummyNeighbors || !ps.removeOnlyHNeighbors ||
        // RDKit✔️✔️:        !ps.removeWithWedgedBond)) {
        // RDKit✔️✔️:     bool onlyHNeighbors = true;
        // RDKit✔️✔️:     for (const auto nbr : mol.atomNeighbors(atom)) {
        // RDKit✔️✔️:       if (!ps.removeDummyNeighbors && nbr->getAtomicNum() < 1) {
        // RDKit✔️✔️:         return false;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       if (!ps.removeOnlyHNeighbors && nbr->getAtomicNum() != 1) {
        // RDKit✔️✔️:         onlyHNeighbors = false;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       if (!ps.removeWithWedgedBond) {
        // RDKit✔️✔️:         const auto bnd = mol.getBondBetweenAtoms(atom->getIdx(), nbr->getIdx());
        // RDKit✔️✔️:         if (bnd->getBondDir() == Bond::BEGINDASH ||
        // RDKit✔️✔️:             bnd->getBondDir() == Bond::BEGINWEDGE) {
        // RDKit✔️✔️:           if (ps.showWarnings) {
        // RDKit✔️✔️:             BOOST_LOG(rdWarningLog) << "WARNING: not removing hydrogen atom "
        // RDKit✔️✔️:                                        "with wedged bond"
        // RDKit✔️✔️:                                     << std::endl;
        // RDKit✔️✔️:           }
        // RDKit✔️✔️:           return false;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       // Check to see if the neighbor has a double bond and we're the only
        // RDKit✔️✔️:       // neighbor at this end.  This was part of github #1810
        // RDKit✔️✔️:       if (!ps.removeDefiningBondStereo && nbr->getDegree() == 2) {
        // RDKit✔️✔️:         for (const auto bnd : mol.atomBonds(nbr)) {
        // RDKit✔️✔️:           if (bnd->getBondType() == Bond::DOUBLE &&
        // RDKit✔️✔️:               (bnd->getStereo() > Bond::STEREOANY ||
        // RDKit✔️✔️:                mol.getBondBetweenAtoms(atom->getIdx(), nbr->getIdx())
        // RDKit✔️✔️:                        ->getBondDir() > Bond::NONE)) {
        // RDKit✔️✔️:             return false;
        // RDKit✔️✔️:           }
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     if (removeIt && (!ps.removeOnlyHNeighbors && onlyHNeighbors)) {
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION removeHs / molRemoveH / shouldRemoveH
        if params.remove_and_track_isotopes {
            let need_remove_hs = self
                .builder
                .atoms()
                .iter()
                .any(|atom| atom.atomic_number() == 1 && atom.isotope().is_none());
            if need_remove_hs {
                let first_pass_params = RemoveHsParams {
                    remove_and_track_isotopes: false,
                    remove_isotopes: false,
                    ..params.clone()
                };
                self.remove_hs_update_explicit_count(&first_pass_params)?;
            }

            let isotopic_hydrogens = tracked_smiles_isotopic_hydrogens(&self.builder)?;
            for atom in self.builder.atoms_mut() {
                if !atom.tracked_isotopic_hydrogens().is_empty() {
                    atom.set_tracked_isotopic_hydrogens(Vec::new());
                }
            }
            for (atom_id, isotopes) in isotopic_hydrogens {
                let atom = self.builder.atom_mut(atom_id).ok_or_else(|| {
                    SmilesParseError::ParseError(format!("missing atom {atom_id}"))
                })?;
                atom.set_tracked_isotopic_hydrogens(isotopes);
            }
        }

        let atom_count = self.builder.atoms().len();
        let mut neighbors = vec![Vec::<(AtomId, BondId)>::new(); atom_count];
        for bond in self.builder.bonds() {
            neighbors[bond.begin().index()].push((bond.end(), bond.id()));
            neighbors[bond.end().index()].push((bond.begin(), bond.id()));
        }

        let mut removal_candidates = Vec::<AtomId>::new();
        for atom in self.builder.atoms() {
            if smiles_should_remove_hydrogen(&self.builder, &neighbors, atom.id(), params)? {
                removal_candidates.push(atom.id());
            }
        }

        let mut atoms_to_remove = removal_candidates;
        if params.remove_in_sgroups {
            filter_smiles_sgroup_emptying_hydrogens(&self.builder, &mut atoms_to_remove);
        }
        let cached_total_valence = smiles_cached_total_valence(&self.builder)?;

        let mut atoms_to_invert = Vec::new();
        let mut atoms_to_set_unknown_stereo = Vec::new();
        let mut bonds_to_clear_stereo = Vec::new();
        let mut bond_direction_updates = Vec::new();
        let mut bond_stereo_updates = Vec::new();
        let mut atoms_to_increment_explicit_h = Vec::new();
        for atom_id in &atoms_to_remove {
            for (neighbor, removed_bond) in &neighbors[atom_id.index()] {
                let hydrogen_bond = self.builder.bond(*removed_bond).ok_or_else(|| {
                    SmilesParseError::ParseError(format!("missing bond {removed_bond}"))
                })?;
                let chiral_tag = self
                    .builder
                    .atoms()
                    .get(neighbor.index())
                    .ok_or_else(|| {
                        SmilesParseError::ParseError(format!("missing atom {neighbor}"))
                    })?
                    .chiral_tag();
                if chiral_tag != ChiralTag::Unspecified {
                    let mut neighbor_indices = neighbors[neighbor.index()]
                        .iter()
                        .filter_map(|(_, bond_id)| {
                            (*bond_id != *removed_bond).then_some(bond_id.index() as i32)
                        })
                        .collect::<Vec<_>>();
                    neighbor_indices.push(removed_bond.index() as i32);
                    if self.perturbation_order(*neighbor, &neighbor_indices)? % 2 == 1 {
                        atoms_to_invert.push(*neighbor);
                    }
                }
                if neighbors[neighbor.index()].len() == 2 {
                    if let Some((_, remaining_bond)) = neighbors[neighbor.index()]
                        .iter()
                        .find(|(_, bond_id)| *bond_id != *removed_bond)
                        && self.builder.bond(*remaining_bond).is_some_and(|bond| {
                            !matches!(bond.stereo(), BondStereo::None | BondStereo::Any)
                        })
                    {
                        bonds_to_clear_stereo.push(*remaining_bond);
                    }
                }
                match hydrogen_bond.direction() {
                    BondDirection::Unknown => {
                        if hydrogen_bond.begin() == *neighbor {
                            atoms_to_set_unknown_stereo.push(*neighbor);
                        }
                    }
                    BondDirection::EndDownRight | BondDirection::EndUpRight => {
                        if let Some(update) = smiles_copy_hydrogen_bond_direction_to_neighbor(
                            &self.builder,
                            &neighbors,
                            *neighbor,
                            *removed_bond,
                        )? {
                            bond_direction_updates.push(update);
                        }
                    }
                    BondDirection::None
                    | BondDirection::BeginWedge
                    | BondDirection::BeginDash
                    | BondDirection::EitherDouble => {}
                }
                if let Some(update) = smiles_adjust_stereo_atoms_if_required(
                    &self.builder,
                    &neighbors,
                    *atom_id,
                    *neighbor,
                )? {
                    bond_stereo_updates.push(update);
                }
                if smiles_should_increment_explicit_hydrogens(
                    &self.builder,
                    *neighbor,
                    params.update_explicit_count,
                    cached_total_valence.as_deref(),
                    &neighbors,
                )? {
                    atoms_to_increment_explicit_h.push(*neighbor);
                }
            }
        }
        for atom_id in atoms_to_increment_explicit_h {
            let atom = self
                .builder
                .atom_mut(atom_id)
                .ok_or_else(|| SmilesParseError::ParseError(format!("missing atom {atom_id}")))?;
            atom.set_explicit_hydrogens(atom.explicit_hydrogens().saturating_add(1));
        }
        for atom_id in atoms_to_invert {
            let atom = self
                .builder
                .atom_mut(atom_id)
                .ok_or_else(|| SmilesParseError::ParseError(format!("missing atom {atom_id}")))?;
            match atom.chiral_tag() {
                ChiralTag::TetrahedralCw => atom.set_chiral_tag(ChiralTag::TetrahedralCcw),
                ChiralTag::TetrahedralCcw => atom.set_chiral_tag(ChiralTag::TetrahedralCw),
                _ => {}
            }
        }
        for atom_id in atoms_to_set_unknown_stereo {
            let atom = self
                .builder
                .atom_mut(atom_id)
                .ok_or_else(|| SmilesParseError::ParseError(format!("missing atom {atom_id}")))?;
            atom.set_unknown_stereo(true);
        }
        for bond_id in bonds_to_clear_stereo {
            let bond = self
                .builder
                .bond_mut(bond_id)
                .ok_or_else(|| SmilesParseError::ParseError(format!("missing bond {bond_id}")))?;
            bond.set_stereo(BondStereo::None);
            bond.set_stereo_atoms(None);
        }
        for (bond_id, direction) in bond_direction_updates {
            let bond = self
                .builder
                .bond_mut(bond_id)
                .ok_or_else(|| SmilesParseError::ParseError(format!("missing bond {bond_id}")))?;
            bond.set_direction(direction);
        }
        for (bond_id, stereo, stereo_atoms) in bond_stereo_updates {
            let bond = self
                .builder
                .bond_mut(bond_id)
                .ok_or_else(|| SmilesParseError::ParseError(format!("missing bond {bond_id}")))?;
            bond.set_stereo(stereo);
            bond.set_stereo_atoms(Some(stereo_atoms));
        }

        if atoms_to_remove.is_empty() {
            return Ok(());
        }
        let mapping = self.builder.remove_atoms_for_construction(&atoms_to_remove);
        self.remap_after_atom_removal(&mapping);
        Ok(())
    }

    fn requires_full_sanitize_mol(&self) -> bool {
        self.builder.atoms().iter().any(|atom| {
            atom.is_aromatic()
                || atom.formal_charge() != 0
                || atom.radical_electrons() != 0
                || atom.query().is_some()
                || !matches!(
                    atom.chiral_tag(),
                    ChiralTag::Unspecified | ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
                )
                || atom.chiral_permutation().is_some()
        }) || self.builder.bonds().iter().any(|bond| {
            bond.is_aromatic()
                || bond.query().is_some()
                || bond.stereo() != BondStereo::None
                || bond.direction() != BondDirection::None
        }) || self.builder.substance_groups_len() != 0
            || self.builder.stereo_groups_len() != 0
            || self.builder.conformers_3d_len() != 0
    }

    fn should_run_full_sanitize_postprocess(&self, params: &SmilesParseParams) -> bool {
        params.sanitize && self.requires_full_sanitize_mol()
    }

    fn remap_after_atom_removal(&mut self, mapping: &[Option<AtomId>]) {
        self.atom_aromatic = self
            .atom_aromatic
            .iter()
            .enumerate()
            .filter_map(|(old_idx, aromatic)| mapping[old_idx].map(|_| *aromatic))
            .collect();
        self.active_atom = self.active_atom.and_then(|atom| mapping[atom.index()]);
        self.smiles_start_atoms = self
            .smiles_start_atoms
            .iter()
            .filter_map(|atom| mapping[atom.index()])
            .collect();
        self.ring_openings.clear();
        self.ring_closures_by_atom.clear();
        self.bond_pairs = self
            .builder
            .bonds()
            .iter()
            .map(|bond| canonical_bond_pair(bond.begin(), bond.end()))
            .collect();
    }
}

fn smiles_should_remove_hydrogen(
    builder: &MoleculeBuilder,
    neighbors: &[Vec<(AtomId, BondId)>],
    atom_id: AtomId,
    params: &RemoveHsParams,
) -> Result<bool, SmilesParseError> {
    let atom = builder
        .atoms()
        .get(atom_id.index())
        .ok_or_else(|| SmilesParseError::ParseError(format!("missing atom {atom_id}")))?;
    if atom.atomic_number() != 1 {
        return Ok(false);
    }
    if !params.remove_with_query && atom.query().is_some() {
        return Ok(false);
    }
    let degree = neighbors[atom_id.index()].len();
    if !params.remove_degree_zero && degree == 0 {
        return Ok(false);
    }
    if !params.remove_higher_degrees && degree > 1 {
        return Ok(false);
    }
    if !params.remove_isotopes && !params.remove_and_track_isotopes && atom.isotope().is_some() {
        return Ok(false);
    }
    if !params.remove_nonimplicit && !atom.implicit_hydrogen() {
        return Ok(false);
    }
    if !params.remove_mapped && atom.atom_map().is_some() {
        return Ok(false);
    }
    if params.remove_in_sgroups {
        if smiles_hydrogen_has_special_sgroup_role(builder, atom_id)? {
            return Ok(false);
        }
    } else if builder
        .substance_groups()
        .iter()
        .any(|substance_group| smiles_sgroup_includes_atom(substance_group, atom_id))
    {
        return Ok(false);
    }
    if !params.remove_hydrides && atom.formal_charge() == -1 {
        return Ok(false);
    }
    if degree != 0
        && (!params.remove_dummy_neighbors
            || !params.remove_defining_bond_stereo
            || !params.remove_only_h_neighbors
            || !params.remove_nontetrahedral_neighbors
            || !params.remove_with_wedged_bond)
    {
        let mut only_h_neighbors = true;
        for (neighbor, bond_id) in &neighbors[atom_id.index()] {
            let neighbor_atom = builder
                .atoms()
                .get(neighbor.index())
                .ok_or_else(|| SmilesParseError::ParseError(format!("missing atom {neighbor}")))?;
            if !params.remove_dummy_neighbors && neighbor_atom.atomic_number() < 1 {
                return Ok(false);
            }
            if !params.remove_nontetrahedral_neighbors
                && crate::chemistry::stereo::has_non_tetrahedral_stereo(neighbor_atom)
            {
                return Ok(false);
            }
            if !params.remove_only_h_neighbors && neighbor_atom.atomic_number() != 1 {
                only_h_neighbors = false;
            }
            if !params.remove_with_wedged_bond {
                let hydrogen_bond = builder.bond(*bond_id).ok_or_else(|| {
                    SmilesParseError::ParseError(format!("missing bond {bond_id}"))
                })?;
                if matches!(
                    hydrogen_bond.direction(),
                    BondDirection::BeginDash | BondDirection::BeginWedge
                ) {
                    return Ok(false);
                }
            }
            if !params.remove_defining_bond_stereo && neighbors[neighbor.index()].len() == 2 {
                let hydrogen_bond = builder.bond(*bond_id).ok_or_else(|| {
                    SmilesParseError::ParseError(format!("missing bond {bond_id}"))
                })?;
                for (_, neighbor_bond_id) in &neighbors[neighbor.index()] {
                    let neighbor_bond = builder.bond(*neighbor_bond_id).ok_or_else(|| {
                        SmilesParseError::ParseError(format!("missing bond {neighbor_bond_id}"))
                    })?;
                    if neighbor_bond.order() == BondOrder::Double
                        && (!matches!(neighbor_bond.stereo(), BondStereo::None | BondStereo::Any)
                            || hydrogen_bond.direction() != BondDirection::None)
                    {
                        return Ok(false);
                    }
                }
            }
        }
        if !params.remove_only_h_neighbors && only_h_neighbors {
            return Ok(false);
        }
    }
    Ok(true)
}

fn smiles_cached_total_valence(
    builder: &MoleculeBuilder,
) -> Result<Option<Vec<i32>>, SmilesParseError> {
    if builder.atoms().is_empty() {
        return Ok(None);
    }
    let molecule = builder
        .clone()
        .build()
        .map_err(|error| SmilesParseError::ParseError(error.to_string()))?;
    let assignment =
        crate::assign_valence_with_options(&molecule, crate::ValenceModel::RdkitLike, false)
            .map_err(|error| SmilesParseError::ParseError(error.to_string()))?;
    Ok(Some(
        assignment
            .explicit_valence
            .iter()
            .zip(assignment.implicit_hydrogens.iter())
            .map(|(explicit, implicit)| *explicit + *implicit)
            .collect(),
    ))
}

fn smiles_should_increment_explicit_hydrogens(
    builder: &MoleculeBuilder,
    heavy_atom: AtomId,
    update_explicit_count: bool,
    cached_total_valence: Option<&[i32]>,
    neighbors: &[Vec<(AtomId, BondId)>],
) -> Result<bool, SmilesParseError> {
    let atom = builder
        .atoms()
        .get(heavy_atom.index())
        .ok_or_else(|| SmilesParseError::ParseError(format!("missing atom {heavy_atom}")))?;
    if update_explicit_count || atom.no_implicit() || atom.chiral_tag() != ChiralTag::Unspecified {
        return Ok(true);
    }
    let Some(total_valence) =
        cached_total_valence.and_then(|values| values.get(heavy_atom.index()))
    else {
        return Ok(false);
    };
    let default_valences = crate::rdkit_valence_list(atom.atomic_number())
        .map_err(|error| SmilesParseError::ParseError(error.to_string()))?;
    let non_default_valence = default_valences
        .is_some_and(|values| values.iter().skip(1).any(|value| *value == *total_valence));
    Ok(((atom.atomic_number() == 7
        || atom.atomic_number() == 15
        || smiles_may_need_extra_h(builder, heavy_atom, cached_total_valence, neighbors)?)
        && atom.is_aromatic())
        || non_default_valence)
}

fn smiles_copy_hydrogen_bond_direction_to_neighbor(
    builder: &MoleculeBuilder,
    neighbors: &[Vec<(AtomId, BondId)>],
    heavy_atom: AtomId,
    hydrogen_bond: BondId,
) -> Result<Option<(BondId, BondDirection)>, SmilesParseError> {
    let hydrogen_bond_data = builder
        .bond(hydrogen_bond)
        .ok_or_else(|| SmilesParseError::ParseError(format!("missing bond {hydrogen_bond}")))?;
    let mut found_direction = false;
    let mut other_single_bond = None;
    for (_, bond_id) in &neighbors[heavy_atom.index()] {
        if *bond_id == hydrogen_bond {
            continue;
        }
        let bond = builder
            .bond(*bond_id)
            .ok_or_else(|| SmilesParseError::ParseError(format!("missing bond {bond_id}")))?;
        if bond.order() != BondOrder::Single {
            continue;
        }
        if bond.direction() == BondDirection::None {
            other_single_bond = Some((bond.id(), bond.begin()));
        } else {
            found_direction = true;
        }
    }
    if found_direction {
        return Ok(None);
    }
    let Some((other_bond_id, other_bond_begin)) = other_single_bond else {
        return Ok(None);
    };
    let mut direction = hydrogen_bond_data.direction();
    if other_bond_begin == heavy_atom && hydrogen_bond_data.begin() == heavy_atom {
        direction = opposite_dir(direction);
    }
    Ok(Some((other_bond_id, direction)))
}

fn smiles_adjust_stereo_atoms_if_required(
    builder: &MoleculeBuilder,
    neighbors: &[Vec<(AtomId, BondId)>],
    atom: AtomId,
    heavy_atom: AtomId,
) -> Result<Option<(BondId, BondStereo, [AtomId; 2])>, SmilesParseError> {
    if neighbors[heavy_atom.index()].len() == 2 {
        return Ok(None);
    }
    for (_, bond_id) in &neighbors[heavy_atom.index()] {
        let bond = builder
            .bond(*bond_id)
            .ok_or_else(|| SmilesParseError::ParseError(format!("missing bond {bond_id}")))?;
        if bond.order() != BondOrder::Double
            || matches!(bond.stereo(), BondStereo::None | BondStereo::Any)
        {
            continue;
        }
        let Some(mut stereo_atoms) = bond.stereo_atoms() else {
            continue;
        };
        if stereo_atoms[0] != atom && stereo_atoms[1] != atom {
            continue;
        }
        let double_neighbor = if bond.begin() == heavy_atom {
            bond.end()
        } else {
            bond.begin()
        };
        for (replacement, _) in &neighbors[heavy_atom.index()] {
            if *replacement == double_neighbor || *replacement == atom {
                continue;
            }
            if stereo_atoms[0] == atom {
                stereo_atoms[0] = *replacement;
            } else {
                stereo_atoms[1] = *replacement;
            }
            let stereo = match bond.stereo() {
                BondStereo::Cis => BondStereo::Trans,
                BondStereo::Trans => BondStereo::Cis,
                other => other,
            };
            return Ok(Some((bond.id(), stereo, stereo_atoms)));
        }
    }
    Ok(None)
}

fn smiles_may_need_extra_h(
    builder: &MoleculeBuilder,
    atom_id: AtomId,
    cached_total_valence: Option<&[i32]>,
    neighbors: &[Vec<(AtomId, BondId)>],
) -> Result<bool, SmilesParseError> {
    let Some(total_valence) = cached_total_valence.and_then(|values| values.get(atom_id.index()))
    else {
        return Ok(false);
    };
    let mut single_bonds = 0usize;
    let mut aromatic_bonds = 0usize;
    for (_, bond_id) in &neighbors[atom_id.index()] {
        match builder
            .bond(*bond_id)
            .ok_or_else(|| SmilesParseError::ParseError(format!("missing bond {bond_id}")))?
            .order()
        {
            BondOrder::Single => single_bonds += 1,
            BondOrder::Aromatic => aromatic_bonds += 1,
            _ => return Ok(false),
        }
    }
    Ok(single_bonds == 1 && aromatic_bonds == 2 && *total_valence == 3)
}

fn tracked_smiles_isotopic_hydrogens(
    builder: &MoleculeBuilder,
) -> Result<Vec<(AtomId, Vec<u16>)>, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION getIsoMap
    // RDKit✔️✔️: std::map<unsigned int, std::vector<unsigned int>> getIsoMap(const ROMol &mol) {
    // RDKit✔️✔️:   std::map<unsigned int, std::vector<unsigned int>> isoMap;
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     if (atom->hasProp(common_properties::_isotopicHs)) {
    // RDKit✔️✔️:       atom->clearProp(common_properties::_isotopicHs);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     auto ba = bond->getBeginAtom();
    // RDKit✔️✔️:     auto ea = bond->getEndAtom();
    // RDKit✔️✔️:     int ha = -1;
    // RDKit✔️✔️:     unsigned int iso;
    // RDKit✔️✔️:     if (ba->getAtomicNum() == 1 && ba->getIsotope() &&
    // RDKit✔️✔️:         ea->getAtomicNum() != 1) {
    // RDKit✔️✔️:       ha = ea->getIdx();
    // RDKit✔️✔️:       iso = ba->getIsotope();
    // RDKit✔️✔️:     } else if (ea->getAtomicNum() == 1 && ea->getIsotope() &&
    // RDKit✔️✔️:                ba->getAtomicNum() != 1) {
    // RDKit✔️✔️:       ha = ba->getIdx();
    // RDKit✔️✔️:       iso = ea->getIsotope();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (ha == -1) { continue; }
    // RDKit✔️✔️:     auto &v = isoMap[ha];
    // RDKit✔️✔️:     v.push_back(iso);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return isoMap;
    // RDKit✔️✔️: }
    let mut map = BTreeMap::<AtomId, Vec<u16>>::new();
    for bond in builder.bonds() {
        let begin = builder.atoms().get(bond.begin().index()).ok_or_else(|| {
            SmilesParseError::ParseError(format!("missing atom {}", bond.begin()))
        })?;
        let end = builder
            .atoms()
            .get(bond.end().index())
            .ok_or_else(|| SmilesParseError::ParseError(format!("missing atom {}", bond.end())))?;
        let tracked = if begin.atomic_number() == 1
            && begin.isotope().is_some()
            && end.atomic_number() != 1
        {
            Some((end.id(), begin.isotope().expect("checked above")))
        } else if end.atomic_number() == 1 && end.isotope().is_some() && begin.atomic_number() != 1
        {
            Some((begin.id(), end.isotope().expect("checked above")))
        } else {
            None
        };
        if let Some((heavy_atom, isotope)) = tracked {
            map.entry(heavy_atom).or_default().push(isotope);
        }
    }
    Ok(map.into_iter().collect())
    // END RDKIT CPP FUNCTION getIsoMap
}

fn smiles_hydrogen_has_special_sgroup_role(
    builder: &MoleculeBuilder,
    hydrogen: AtomId,
) -> Result<bool, SmilesParseError> {
    for substance_group in builder.substance_groups() {
        for bond_id in substance_group.bonds() {
            let bond = builder
                .bond(*bond_id)
                .ok_or_else(|| SmilesParseError::ParseError(format!("missing bond {bond_id}")))?;
            let begin_in_group = substance_group.atoms().contains(&bond.begin());
            let end_in_group = substance_group.atoms().contains(&bond.end());
            let is_xbond = match (begin_in_group, end_in_group) {
                (true, true) => false,
                (true, false) | (false, true) => true,
                (false, false) => {
                    return Err(SmilesParseError::ParseError(format!(
                        "substance-group bond {bond_id} is neither crossing nor contained"
                    )));
                }
            };
            if is_xbond && (bond.begin() == hydrogen || bond.end() == hydrogen) {
                return Ok(true);
            }
        }
        if substance_group
            .attach_points()
            .iter()
            .any(|attach_point| attach_point.atom == hydrogen)
        {
            return Ok(true);
        }
        for cstate in substance_group.cstates() {
            let bond = builder.bond(cstate.bond).ok_or_else(|| {
                SmilesParseError::ParseError(format!("missing bond {}", cstate.bond))
            })?;
            if bond.begin() == hydrogen || bond.end() == hydrogen {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

fn smiles_sgroup_includes_atom(substance_group: &SubstanceGroup, atom: AtomId) -> bool {
    substance_group.atoms().contains(&atom)
        || substance_group.parent_atoms().contains(&atom)
        || substance_group.attach_points().iter().any(|attach_point| {
            attach_point.atom == atom || attach_point.leaving_atom == Some(atom)
        })
}

fn filter_smiles_sgroup_emptying_hydrogens(
    builder: &MoleculeBuilder,
    atoms_to_remove: &mut Vec<AtomId>,
) {
    let original_atoms_to_remove = atoms_to_remove.clone();
    atoms_to_remove.retain(|atom| {
        builder.substance_groups().iter().all(|substance_group| {
            let atoms = substance_group.atoms();
            atoms.is_empty()
                || atoms
                    .iter()
                    .any(|member| member != atom && !original_atoms_to_remove.contains(member))
        })
    });
}

fn get_unspecified_bond_type_for_atoms(atom1_aromatic: bool, atom2_aromatic: bool) -> BondOrder {
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

fn bond_order_as_double(order: BondOrder) -> f64 {
    match order {
        BondOrder::Single
        | BondOrder::Dative
        | BondOrder::DativeOne
        | BondOrder::DativeLeft
        | BondOrder::DativeRight => 1.0,
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Quadruple => 4.0,
        BondOrder::Quintuple => 5.0,
        BondOrder::Hextuple => 6.0,
        BondOrder::OneAndHalf | BondOrder::Aromatic => 1.5,
        BondOrder::TwoAndHalf => 2.5,
        BondOrder::ThreeAndHalf => 3.5,
        BondOrder::FourAndHalf => 4.5,
        BondOrder::FiveAndHalf => 5.5,
        BondOrder::Null
        | BondOrder::Ionic
        | BondOrder::Hydrogen
        | BondOrder::ThreeCenter
        | BondOrder::Other
        | BondOrder::Zero
        | BondOrder::Unspecified => 0.0,
    }
}

fn atom_query_has_single_h_count(query: &QueryNode<AtomQueryPredicate>) -> bool {
    // BEGIN RDKIT CPP FUNCTION Canon::details::hasSingleHQuery
    // RDKit✔️✔️: bool hasSingleHQuery(const Atom::QUERYATOM_QUERY *q) {
    // RDKit✔️✔️:   // list queries are series of nested ors of AtomAtomicNum queries
    // RDKit✔️✔️:   PRECONDITION(q, "bad query");
    // RDKit✔️✔️:   bool res = false;
    // RDKit✔️✔️:   const auto &descr = q->getDescription();
    // RDKit✔️✔️:   if (descr == "AtomAnd") {
    // RDKit✔️✔️:     for (auto cIt = q->beginChildren(); cIt != q->endChildren(); ++cIt) {
    // RDKit✔️✔️:       const auto &cDescr = (*cIt)->getDescription();
    // RDKit✔️✔️:       if (cDescr == "AtomHCount") {
    // RDKit✔️✔️:         return !(*cIt)->getNegation() &&
    // RDKit✔️✔️:                ((ATOM_EQUALS_QUERY *)(*cIt).get())->getVal() == 1;
    // RDKit✔️✔️:       } else if (cDescr == "AtomAnd") {
    // RDKit✔️✔️:         res = hasSingleHQuery((*cIt).get());
    // RDKit✔️✔️:         if (res) {
    // RDKit✔️✔️:           return true;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Canon::details::hasSingleHQuery
    match query {
        QueryNode::And(children) => {
            for child in children {
                match child {
                    QueryNode::Predicate(AtomQueryPredicate::ImplicitHydrogenCount(1)) => {
                        return true;
                    }
                    QueryNode::And(_) if atom_query_has_single_h_count(child) => {
                        return true;
                    }
                    _ => {}
                }
            }
            false
        }
        _ => false,
    }
}

const SWAP_SQUAREPLANAR_TABLE: [[u8; 6]; 4] = [
    [0, 0, 0, 0, 0, 0],
    [3, 1, 2, 2, 1, 3],
    [2, 3, 1, 1, 3, 2],
    [1, 2, 3, 3, 2, 1],
];

const SWAP_TRIGONALBIPYRAMIDAL_TABLE: [[u8; 10]; 21] = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [9, 20, 17, 2, 2, 2, 7, 2, 6, 3],
    [11, 15, 18, 1, 1, 1, 8, 1, 5, 4],
    [10, 19, 4, 18, 4, 8, 4, 5, 4, 1],
    [12, 16, 3, 17, 3, 7, 3, 6, 3, 2],
    [13, 6, 16, 20, 7, 6, 6, 3, 2, 6],
    [14, 5, 19, 15, 8, 5, 5, 4, 1, 5],
    [8, 14, 10, 11, 5, 4, 1, 8, 8, 8],
    [7, 13, 12, 9, 6, 3, 2, 7, 7, 7],
    [1, 11, 11, 8, 15, 18, 11, 11, 14, 10],
    [3, 12, 7, 12, 16, 12, 17, 13, 12, 9],
    [2, 9, 9, 7, 20, 17, 9, 9, 13, 12],
    [4, 10, 8, 10, 19, 10, 18, 14, 10, 11],
    [5, 8, 14, 14, 14, 19, 15, 10, 11, 14],
    [6, 7, 13, 13, 13, 16, 20, 12, 9, 13],
    [20, 2, 20, 6, 9, 20, 13, 17, 20, 16],
    [19, 4, 5, 19, 10, 14, 19, 19, 18, 15],
    [18, 18, 1, 4, 18, 11, 10, 15, 19, 18],
    [17, 17, 2, 3, 17, 9, 12, 20, 16, 17],
    [16, 3, 6, 16, 12, 13, 16, 16, 17, 20],
    [15, 1, 15, 5, 11, 15, 14, 18, 15, 19],
];

const SWAP_OCTAHEDRAL_TABLE: [[u8; 15]; 31] = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [17, 16, 30, 21, 2, 14, 2, 10, 25, 8, 2, 22, 4, 7, 3],
    [7, 3, 25, 22, 1, 4, 1, 8, 30, 10, 1, 21, 14, 17, 16],
    [18, 2, 29, 16, 22, 15, 16, 26, 11, 9, 21, 16, 6, 5, 1],
    [15, 18, 19, 28, 14, 2, 8, 14, 27, 14, 10, 24, 1, 6, 5],
    [14, 17, 20, 15, 27, 16, 9, 28, 15, 15, 23, 11, 7, 3, 4],
    [16, 14, 18, 26, 24, 17, 29, 18, 13, 19, 12, 18, 3, 4, 7],
    [2, 15, 17, 23, 25, 18, 30, 12, 17, 20, 17, 13, 5, 1, 6],
    [23, 26, 11, 12, 10, 10, 4, 2, 29, 1, 14, 20, 10, 13, 9],
    [24, 25, 10, 11, 13, 11, 5, 30, 16, 3, 19, 15, 12, 11, 8],
    [20, 29, 9, 13, 8, 8, 14, 1, 26, 2, 4, 23, 8, 12, 11],
    [19, 30, 8, 9, 12, 9, 15, 25, 3, 16, 24, 5, 13, 9, 10],
    [22, 27, 13, 8, 11, 13, 28, 7, 18, 21, 6, 17, 9, 10, 13],
    [21, 28, 12, 10, 9, 12, 27, 17, 6, 22, 18, 7, 11, 8, 12],
    [5, 6, 24, 27, 4, 1, 10, 4, 28, 4, 8, 19, 2, 18, 15],
    [4, 7, 23, 5, 28, 3, 11, 27, 5, 5, 20, 9, 17, 16, 14],
    [6, 1, 26, 3, 21, 5, 3, 29, 9, 11, 22, 3, 18, 15, 2],
    [1, 5, 7, 20, 30, 6, 25, 13, 7, 23, 7, 12, 15, 2, 18],
    [3, 4, 6, 29, 19, 7, 26, 6, 12, 24, 13, 6, 16, 14, 17],
    [11, 24, 4, 30, 18, 25, 23, 24, 22, 6, 9, 14, 21, 24, 20],
    [10, 23, 5, 17, 29, 26, 24, 21, 23, 7, 15, 8, 23, 22, 19],
    [13, 22, 28, 1, 16, 27, 22, 20, 24, 12, 3, 2, 19, 23, 22],
    [12, 21, 27, 2, 3, 28, 21, 23, 19, 13, 16, 1, 24, 20, 21],
    [8, 20, 15, 7, 26, 29, 19, 22, 20, 17, 5, 10, 20, 21, 24],
    [9, 19, 14, 25, 6, 30, 20, 19, 21, 18, 11, 4, 22, 19, 23],
    [30, 9, 2, 24, 7, 19, 17, 11, 1, 29, 30, 28, 27, 30, 26],
    [29, 8, 16, 6, 23, 20, 18, 3, 10, 30, 27, 29, 29, 28, 25],
    [28, 12, 22, 14, 5, 21, 13, 15, 4, 28, 26, 30, 25, 29, 28],
    [27, 13, 21, 4, 15, 22, 12, 5, 14, 27, 29, 25, 30, 26, 27],
    [26, 10, 3, 18, 20, 23, 6, 16, 8, 25, 28, 26, 26, 27, 30],
    [25, 11, 1, 19, 17, 24, 7, 9, 2, 26, 25, 27, 28, 25, 29],
];

fn swap_squareplanar(perm: u32, x: usize, y: usize) -> u32 {
    if perm as usize >= SWAP_SQUAREPLANAR_TABLE.len() {
        return 0;
    }
    swap_nontetrahedral(perm, x, y, &[0, 2, 3], 3, |perm, swapidx| {
        SWAP_SQUAREPLANAR_TABLE[perm][swapidx]
    })
}

fn swap_trigonalbipyramidal(perm: u32, x: usize, y: usize) -> u32 {
    if perm as usize >= SWAP_TRIGONALBIPYRAMIDAL_TABLE.len() {
        return 0;
    }
    swap_nontetrahedral(perm, x, y, &[0, 3, 5, 6], 4, |perm, swapidx| {
        SWAP_TRIGONALBIPYRAMIDAL_TABLE[perm][swapidx]
    })
}

fn swap_octahedral(perm: u32, x: usize, y: usize) -> u32 {
    if perm as usize >= SWAP_OCTAHEDRAL_TABLE.len() {
        return 0;
    }
    swap_nontetrahedral(perm, x, y, &[0, 4, 7, 9, 10], 5, |perm, swapidx| {
        SWAP_OCTAHEDRAL_TABLE[perm][swapidx]
    })
}

fn swap_nontetrahedral(
    perm: u32,
    x: usize,
    y: usize,
    offset: &[usize],
    max_index: usize,
    table_lookup: impl Fn(usize, usize) -> u8,
) -> u32 {
    if x == y {
        return perm;
    }
    let swapidx = if x < y {
        if y > max_index {
            return 0;
        }
        offset[x] + (y - 1)
    } else {
        if x > max_index {
            return 0;
        }
        offset[y] + (x - 1)
    };
    table_lookup(perm as usize, swapidx) as u32
}

pub(crate) fn nontetrahedral_max_neighbors(chiral_tag: ChiralTag) -> Option<usize> {
    match chiral_tag {
        ChiralTag::SquarePlanar => Some(4),
        ChiralTag::TrigonalBipyramidal => Some(5),
        ChiralTag::Octahedral => Some(6),
        _ => None,
    }
}

pub(crate) fn nontetrahedral_chiral_permutation_for_probe(
    molecule: &Molecule,
    atom_id: AtomId,
    probe: &[Option<BondId>],
    inverse: bool,
) -> Result<u32, String> {
    // BEGIN RDKIT CPP FUNCTION Chirality::getChiralPermutation
    // RDKit✔️✔️: unsigned int getChiralPermutation(const Atom *cen, const INT_LIST &probe,
    // RDKit✔️✔️:                                   bool inverse) {
    // RDKit✔️✔️:   int perm;
    // RDKit✔️✔️:   if (!cen->getPropIfPresent(common_properties::_chiralPermutation, perm) ||
    // RDKit✔️✔️:       perm <= 0) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   decltype(&swap_octahedral) swap_func = nullptr;
    // RDKit✔️✔️:   switch (cen->getChiralTag()) { ... }
    // RDKit✔️✔️:   std::vector<int> order(cen->getOwningMol().getNumBonds(), -1);
    // RDKit✔️✔️:   for (const auto bnd : cen->getOwningMol().atomBonds(cen)) {
    // RDKit✔️✔️:     order[bnd->getIdx()] = nbrIdx++;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::vector<unsigned int> nbrPerm(nbrIdx);
    // RDKit✔️✔️:   std::iota(nbrPerm.begin(), nbrPerm.end(), 0);
    // RDKit✔️✔️:   std::vector<unsigned int> probePerm(probe.size());
    // RDKit✔️✔️:   for (auto v : probe) {
    // RDKit✔️✔️:     probePerm[nbrIdx++] = v < 0 ? -1 : order[v];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (nbrPerm.size() < nbrIdx) {
    // RDKit✔️✔️:     nbrPerm.insert(nbrPerm.end(), nbrIdx - nbrPerm.size(), -1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (inverse) {
    // RDKit✔️✔️:     std::swap(nbrPerm, probePerm);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (unsigned int i = 0; i < probePerm.size() - 1; ++i) { ... }
    // RDKit✔️✔️:   return perm;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Chirality::getChiralPermutation
    let atom = molecule
        .atoms()
        .get(atom_id.index())
        .ok_or_else(|| format!("missing atom {atom_id}"))?;
    let mut perm = atom.chiral_permutation().unwrap_or(0);
    if perm == 0 {
        return Ok(0);
    }
    let swap_func: fn(u32, usize, usize) -> u32 = match atom.chiral_tag() {
        ChiralTag::SquarePlanar => {
            if probe.len() > 4 {
                return Ok(0);
            }
            swap_squareplanar
        }
        ChiralTag::TrigonalBipyramidal => {
            if probe.len() > 5 {
                return Ok(0);
            }
            swap_trigonalbipyramidal
        }
        ChiralTag::Octahedral => {
            if probe.len() > 6 {
                return Ok(0);
            }
            swap_octahedral
        }
        _ => return Ok(0),
    };

    let mut order = vec![-1isize; molecule.num_bonds()];
    let mut nbr_idx = 0isize;
    for bond in molecule.bonds() {
        if bond.begin() == atom_id || bond.end() == atom_id {
            order[bond.id().index()] = nbr_idx;
            nbr_idx += 1;
        }
    }
    let mut nbr_perm = (0..nbr_idx).collect::<Vec<_>>();
    let mut probe_perm = Vec::with_capacity(probe.len());
    for bond in probe {
        match bond {
            Some(bond_id) => probe_perm.push(order[bond_id.index()]),
            None => probe_perm.push(-1),
        }
    }
    if nbr_perm.len() < probe_perm.len() {
        nbr_perm.extend(std::iter::repeat_n(-1, probe_perm.len() - nbr_perm.len()));
    }
    if inverse {
        std::mem::swap(&mut nbr_perm, &mut probe_perm);
    }
    for i in 0..probe_perm.len().saturating_sub(1) {
        let pval = probe_perm[i];
        if nbr_perm[i] == pval {
            continue;
        }
        let Some(target) = nbr_perm[i..]
            .iter()
            .position(|value| *value == pval)
            .map(|offset| i + offset)
        else {
            return Err("could not find probe element".to_string());
        };
        perm = swap_func(perm, i, target);
        nbr_perm.swap(i, target);
    }
    Ok(perm)
}

pub fn assign_double_bond_stereo_from_directions(
    molecule: &mut Molecule,
) -> Result<(), StereoError> {
    // BEGIN RDKIT CPP FUNCTION setBondStereoFromDirections caller path
    // RDKit✔️✔️: MolOps::setBondStereoFromDirections(*res);
    // RDKit✔️✔️: // clear _needsDetectBondStereo and assign E/Z plus stereo atoms from
    // RDKit✔️✔️: // neighboring directed single bonds.
    // END RDKIT CPP FUNCTION setBondStereoFromDirections caller path
    set_bond_stereo_from_directions(molecule)
}

#[cfg(test)]
mod tests;
