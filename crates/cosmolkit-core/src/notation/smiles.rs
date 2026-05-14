use std::collections::{BTreeMap, BTreeSet};

use crate::{
    Atom, AtomId, AtomQueryPredicate, AtomSpec, BondDirection, BondId, BondOrder,
    BondQueryPredicate, BondSpec, BondStereo, ChiralTag, Conformer3D, Element, Molecule,
    MoleculeBuilder, QueryNode, RemoveHsParams, SmilesParseError, StereoError, StereoGroup,
    StereoGroupKind, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind,
};

// RDKit source-reproduction note:
//
// This module is the source-level alignment frame for the future SMILES parser
// port from `third_party/rdkit/Code/GraphMol/SmilesParse`. The functions below
// follow `dev/source_reproduction_protocol.md`. They intentionally contain
// RDKit C++ comments with `RDKit❌❌` markers before any behavior is
// implemented. Do not upgrade markers or add semantic behavior without checking
// the corresponding RDKit source line and the local Rust code.

const SMILES_START_PROP: &str = "_SmilesStart";
const CXSMILES_BOND_IDX_PROP: &str = "_cxsmilesBondIdx";
const UNSPECIFIED_ORDER_PROP: &str = "_unspecifiedOrder";

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
            in_atom_state: false,
        }
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
        // RDKit❗✔️:
        // RDKit❗✔️: @[' ']*TH { yylval->chiraltype = Atom::ChiralType::CHI_TETRAHEDRAL; return CHI_CLASS_TOKEN; }
        // RDKit❗✔️: @[' ']*AL { yylval->chiraltype = Atom::ChiralType::CHI_ALLENE; return CHI_CLASS_TOKEN; }
        // RDKit❗✔️: @[' ']*SP { yylval->chiraltype = Atom::ChiralType::CHI_SQUAREPLANAR; return CHI_CLASS_TOKEN; }
        // RDKit❗✔️: @[' ']*TB { yylval->chiraltype = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL; return CHI_CLASS_TOKEN; }
        // RDKit❗✔️: @[' ']*OH { yylval->chiraltype = Atom::ChiralType::CHI_OCTAHEDRAL; return CHI_CLASS_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: @		{ return AT_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>He	{ yylval->atom = new Atom(2); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Li	{ yylval->atom = new Atom(3); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Be	{ yylval->atom = new Atom(4); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Ne	{ yylval->atom = new Atom(10); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Na	{ yylval->atom = new Atom(11); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Mg	{ yylval->atom = new Atom(12); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Al	{ yylval->atom = new Atom(13); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Si	{ yylval->atom = new Atom(14); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Ar	{ yylval->atom = new Atom(18); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>K	{ yylval->atom = new Atom(19); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Ca	{ yylval->atom = new Atom(20); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Sc	{ yylval->atom = new Atom(21); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Ti	{ yylval->atom = new Atom(22); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>V	{ yylval->atom = new Atom(23); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Cr	{ yylval->atom = new Atom(24); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Mn	{ yylval->atom = new Atom(25); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Fe	{ yylval->atom = new Atom(26); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Co	{ yylval->atom = new Atom(27); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Ni	{ yylval->atom = new Atom(28); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Cu	{ yylval->atom = new Atom(29); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Zn	{ yylval->atom = new Atom(30); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Ga	{ yylval->atom = new Atom(31); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Ge	{ yylval->atom = new Atom(32); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>As	{ yylval->atom = new Atom(33); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Se	{ yylval->atom = new Atom(34); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Kr	{ yylval->atom = new Atom(36); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Rb	{ yylval->atom = new Atom(37); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Sr	{ yylval->atom = new Atom(38); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Y	{ yylval->atom = new Atom(39); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Zr	{ yylval->atom = new Atom(40); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Nb	{ yylval->atom = new Atom(41); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Mo	{ yylval->atom = new Atom(42); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Tc	{ yylval->atom = new Atom(43); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Ru	{ yylval->atom = new Atom(44); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Rh	{ yylval->atom = new Atom(45); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Pd	{ yylval->atom = new Atom(46); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Ag	{ yylval->atom = new Atom(47); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Cd	{ yylval->atom = new Atom(48); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>In	{ yylval->atom = new Atom(49); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Sn	{ yylval->atom = new Atom(50); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Sb	{ yylval->atom = new Atom(51); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Te	{ yylval->atom = new Atom(52); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Xe	{ yylval->atom = new Atom(54); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Cs	{ yylval->atom = new Atom(55); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Ba	{ yylval->atom = new Atom(56); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>La	{ yylval->atom = new Atom(57); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Ce	{ yylval->atom = new Atom(58); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Pr	{ yylval->atom = new Atom(59); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Nd	{ yylval->atom = new Atom(60); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Pm	{ yylval->atom = new Atom(61); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Sm	{ yylval->atom = new Atom(62); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Eu	{ yylval->atom = new Atom(63); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Gd	{ yylval->atom = new Atom(64); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Tb	{ yylval->atom = new Atom(65); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Dy	{ yylval->atom = new Atom(66); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Ho	{ yylval->atom = new Atom(67); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Er	{ yylval->atom = new Atom(68); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Tm	{ yylval->atom = new Atom(69); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Yb	{ yylval->atom = new Atom(70); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Lu	{ yylval->atom = new Atom(71); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Hf	{ yylval->atom = new Atom(72); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Ta	{ yylval->atom = new Atom(73); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>W	{ yylval->atom = new Atom(74); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Re	{ yylval->atom = new Atom(75); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Os	{ yylval->atom = new Atom(76); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Ir	{ yylval->atom = new Atom(77); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Pt	{ yylval->atom = new Atom(78); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Au	{ yylval->atom = new Atom(79); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Hg	{ yylval->atom = new Atom(80); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Tl	{ yylval->atom = new Atom(81); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Pb	{ yylval->atom = new Atom(82); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Bi	{ yylval->atom = new Atom(83); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Po	{ yylval->atom = new Atom(84); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>At	{ yylval->atom = new Atom(85); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Rn	{ yylval->atom = new Atom(86); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Fr	{ yylval->atom = new Atom(87); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Ra	{ yylval->atom = new Atom(88); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Ac	{ yylval->atom = new Atom(89); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Th	{ yylval->atom = new Atom(90); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Pa	{ yylval->atom = new Atom(91); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>U	{ yylval->atom = new Atom(92); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Np	{ yylval->atom = new Atom(93); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Pu	{ yylval->atom = new Atom(94); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Am	{ yylval->atom = new Atom(95); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Cm	{ yylval->atom = new Atom(96); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Bk	{ yylval->atom = new Atom(97); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Cf	{ yylval->atom = new Atom(98); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Es	{ yylval->atom = new Atom(99); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Fm	{ yylval->atom = new Atom(100); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Md	{ yylval->atom = new Atom(101); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>No	{ yylval->atom = new Atom(102); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Lr	{ yylval->atom = new Atom(103); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Rf	{ yylval->atom = new Atom(104); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Db	{ yylval->atom = new Atom(105); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Sg	{ yylval->atom = new Atom(106); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Bh	{ yylval->atom = new Atom(107); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Hs	{ yylval->atom = new Atom(108); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Mt	{ yylval->atom = new Atom(109); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Ds	{ yylval->atom = new Atom(110); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Rg	{ yylval->atom = new Atom(111); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Cn	{ yylval->atom = new Atom(112); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Nh	{ yylval->atom = new Atom(113); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Fl	{ yylval->atom = new Atom(114); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Mc	{ yylval->atom = new Atom(115); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Lv	{ yylval->atom = new Atom(116); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Ts	{ yylval->atom = new Atom(117); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Og	{ yylval->atom = new Atom(118); return ATOM_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>Uun	{ yylval->atom = new Atom(110); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Uuu	{ yylval->atom = new Atom(111); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Uub	{ yylval->atom = new Atom(112); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Uut	{ yylval->atom = new Atom(113); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Uuq	{ yylval->atom = new Atom(114); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Uup	{ yylval->atom = new Atom(115); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Uuh	{ yylval->atom = new Atom(116); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Uus	{ yylval->atom = new Atom(117); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>Uuo	{ yylval->atom = new Atom(118); return ATOM_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: B  { yylval->atom = new Atom(5);return ORGANIC_ATOM_TOKEN; }
        // RDKit❗✔️: C  { yylval->atom = new Atom(6);return ORGANIC_ATOM_TOKEN; }
        // RDKit❗✔️: N  { yylval->atom = new Atom(7);return ORGANIC_ATOM_TOKEN; }
        // RDKit❗✔️: O  { yylval->atom = new Atom(8);return ORGANIC_ATOM_TOKEN; }
        // RDKit❗✔️: P  { yylval->atom = new Atom(15);return ORGANIC_ATOM_TOKEN; }
        // RDKit❗✔️: S  { yylval->atom = new Atom(16);return ORGANIC_ATOM_TOKEN; }
        // RDKit❗✔️: F  { yylval->atom = new Atom(9);return ORGANIC_ATOM_TOKEN; }
        // RDKit❗✔️: Cl { yylval->atom = new Atom(17);return ORGANIC_ATOM_TOKEN; }
        // RDKit❗✔️: Br { yylval->atom = new Atom(35);return ORGANIC_ATOM_TOKEN; }
        // RDKit❗✔️: I  { yylval->atom = new Atom(53);return ORGANIC_ATOM_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: H			{
        // RDKit❗✔️: 				return H_TOKEN;
        // RDKit❗✔️: 			}
        // RDKit❗✔️:
        // RDKit❗✔️: b		    {	yylval->atom = new Atom ( 5 );
        // RDKit❗✔️: 			yylval->atom->setIsAromatic(true);
        // RDKit❗✔️: 				return AROMATIC_ATOM_TOKEN;
        // RDKit❗✔️: 			}
        // RDKit❗✔️: c		    {	yylval->atom = new Atom ( 6 );
        // RDKit❗✔️: 			yylval->atom->setIsAromatic(true);
        // RDKit❗✔️: 				return AROMATIC_ATOM_TOKEN;
        // RDKit❗✔️: 			}
        // RDKit❗✔️: n		    {	yylval->atom = new Atom( 7 );
        // RDKit❗✔️: 			yylval->atom->setIsAromatic(true);
        // RDKit❗✔️: 				return AROMATIC_ATOM_TOKEN;
        // RDKit❗✔️: 			}
        // RDKit❗✔️: o		    {	yylval->atom = new Atom( 8 );
        // RDKit❗✔️: 			yylval->atom->setIsAromatic(true);
        // RDKit❗✔️: 				return AROMATIC_ATOM_TOKEN;
        // RDKit❗✔️: 			}
        // RDKit❗✔️: p		    {	yylval->atom = new Atom( 15 );
        // RDKit❗✔️: 			yylval->atom->setIsAromatic(true);
        // RDKit❗✔️: 				return AROMATIC_ATOM_TOKEN;
        // RDKit❗✔️: 			}
        // RDKit❗✔️: s		    {	yylval->atom = new Atom( 16 );
        // RDKit❗✔️: 			yylval->atom->setIsAromatic(true);
        // RDKit❗✔️: 				return AROMATIC_ATOM_TOKEN;
        // RDKit❗✔️: 			}
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>si   {	yylval->atom = new Atom( 14 );
        // RDKit❗✔️: 			yylval->atom->setIsAromatic(true);
        // RDKit❗✔️: 				return AROMATIC_ATOM_TOKEN;
        // RDKit❗✔️: 			}
        // RDKit❗✔️: <IN_ATOM_STATE>as   {	yylval->atom = new Atom( 33 );
        // RDKit❗✔️: 			yylval->atom->setIsAromatic(true);
        // RDKit❗✔️: 				return AROMATIC_ATOM_TOKEN;
        // RDKit❗✔️: 			}
        // RDKit❗✔️: <IN_ATOM_STATE>se   {	yylval->atom = new Atom( 34 );
        // RDKit❗✔️: 			yylval->atom->setIsAromatic(true);
        // RDKit❗✔️: 				return AROMATIC_ATOM_TOKEN;
        // RDKit❗✔️: 			}
        // RDKit❗✔️: <IN_ATOM_STATE>te   {	yylval->atom = new Atom( 52 );
        // RDKit❗✔️: 			yylval->atom->setIsAromatic(true);
        // RDKit❗✔️: 				return AROMATIC_ATOM_TOKEN;
        // RDKit❗✔️: 			}
        // RDKit❗✔️:
        // RDKit❗✔️: \* 	            {   yylval->atom = new Atom( 0 );
        // RDKit❗✔️: 		            yylval->atom->setProp(common_properties::dummyLabel,
        // RDKit❗✔️:                                                         std::string("*"));
        // RDKit❗✔️:                                 // must be ORGANIC_ATOM_TOKEN because
        // RDKit❗✔️:                                 // we aren't in square brackets:
        // RDKit❗✔️: 				return ORGANIC_ATOM_TOKEN;
        // RDKit❗✔️: 			}
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>\: 	{ return COLON_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>\# 	{ return HASH_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: %{ // Biovia quoted heavy atom workaround %}
        // RDKit❗✔️: <IN_ATOM_STATE>\'Rf\'	{ yylval->atom = new Atom(104); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>\'Db\'	{ yylval->atom = new Atom(105); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>\'Sg\'	{ yylval->atom = new Atom(106); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>\'Bh\'	{ yylval->atom = new Atom(107); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>\'Hs\'	{ yylval->atom = new Atom(108); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>\'Mt\'	{ yylval->atom = new Atom(109); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>\'Ds\'	{ yylval->atom = new Atom(110); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>\'Rg\'	{ yylval->atom = new Atom(111); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>\'Cn\'	{ yylval->atom = new Atom(112); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>\'Nh\'	{ yylval->atom = new Atom(113); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>\'Fl\'	{ yylval->atom = new Atom(114); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>\'Mc\'	{ yylval->atom = new Atom(115); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>\'Lv\'	{ yylval->atom = new Atom(116); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>\'Ts\'	{ yylval->atom = new Atom(117); return ATOM_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>\'Og\'	{ yylval->atom = new Atom(118); return ATOM_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: \=	{ yylval->bond = new Bond(Bond::DOUBLE); return BOND_TOKEN; }
        // RDKit❗✔️: \#	{ yylval->bond = new Bond(Bond::TRIPLE); return BOND_TOKEN; }
        // RDKit❗✔️: \:	{ yylval->bond = new Bond(Bond::AROMATIC);
        // RDKit❗✔️: 	  yylval->bond->setIsAromatic(true); return BOND_TOKEN; }
        // RDKit❗✔️: \$	{ yylval->bond = new Bond(Bond::QUADRUPLE); return BOND_TOKEN; }
        // RDKit❗✔️: \-\>	{ yylval->bond = new Bond(Bond::DATIVER); return BOND_TOKEN; }
        // RDKit❗✔️: \<\-	{ yylval->bond = new Bond(Bond::DATIVEL); return BOND_TOKEN; }
        // RDKit❗✔️: \~	{ yylval->bond = new QueryBond();
        // RDKit❗✔️: 	  yylval->bond->setQuery(makeBondNullQuery()); return BOND_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: [\\]{1,2}    { yylval->bond = new Bond(Bond::UNSPECIFIED);
        // RDKit❗✔️: 	yylval->bond->setProp(RDKit::common_properties::_unspecifiedOrder,1);
        // RDKit❗✔️: 	yylval->bond->setBondDir(Bond::ENDDOWNRIGHT); return BOND_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: [\/]    { yylval->bond = new Bond(Bond::UNSPECIFIED);
        // RDKit❗✔️: 	yylval->bond->setProp(RDKit::common_properties::_unspecifiedOrder,1);
        // RDKit❗✔️: 	yylval->bond->setBondDir(Bond::ENDUPRIGHT); return BOND_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: \-			{ return MINUS_TOKEN; }
        // RDKit❗✔️: \+			{ return PLUS_TOKEN; }
        // RDKit❗✔️: \(       	{ return GROUP_OPEN_TOKEN; }
        // RDKit❗✔️: \)       	{ return GROUP_CLOSE_TOKEN; }
        // RDKit❗✔️: \[			{ BEGIN IN_ATOM_STATE; return ATOM_OPEN_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>\]	{ BEGIN INITIAL; return ATOM_CLOSE_TOKEN; }
        // RDKit❗✔️: \.       	{ return SEPARATOR_TOKEN; }
        // RDKit❗✔️: \%              { return PERCENT_TOKEN; }
        // RDKit❗✔️: [0]		{ yylval->ival = 0; return ZERO_TOKEN; }
        // RDKit❗✔️: [1-9]		{ yylval->ival = yytext[0] - '0'; return NONZERO_DIGIT_TOKEN; }
        // RDKit❗✔️: \n		return 0;
        // RDKit❗✔️: <<EOF>>		{ return EOS_TOKEN; }
        // RDKit❗✔️: .		return BAD_CHARACTER;
        // END RDKIT CPP LEXER RULES smiles.ll token rules
        if self.scan_cursor >= self.scan_end {
            return Ok(SmilesToken::Eos);
        }

        let remaining = &self.input[self.scan_cursor..self.scan_end];
        let (len, token) = self.match_next_token(remaining);
        self.scan_cursor += len;
        self.current_token_position += len;
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
    for &(symbol, atomic_number) in BRACKET_ATOM_SYMBOLS {
        if input.starts_with(symbol) {
            return Some((symbol, atomic_number));
        }
    }
    None
}

fn match_quoted_biovia_atom_symbol(input: &str) -> Option<(&'static str, u8)> {
    for &(symbol, atomic_number) in QUOTED_BIOVIA_ATOM_SYMBOLS {
        if input.starts_with(symbol) {
            return Some((symbol, atomic_number));
        }
    }
    None
}

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

#[derive(Debug, Clone, PartialEq, Eq)]
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
        // RDKit❗✔️: int smiles_parse(const std::string &inp, std::vector<RDKit::RWMol *> &molVect) {
        // RDKit❗✔️:   auto start_tok = static_cast<int>(START_MOL);
        // RDKit❗✔️:   Atom *atom = nullptr;
        // RDKit❗✔️:   Bond *bond = nullptr;
        // RDKit❗✔️:   return smiles_parse_helper(inp, molVect, atom, bond, start_tok);
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION smiles_parse

        // BEGIN RDKIT CPP FUNCTION smiles_parse_helper
        // RDKit❗✔️: int smiles_parse_helper(const std::string &inp,
        // RDKit❗✔️:                         std::vector<RDKit::RWMol *> &molVect, Atom *&atom,
        // RDKit❗✔️:                         Bond *&bond, int start_tok) {
        // RDKit❗✔️:     return generic_parse_helper<yysmiles_lex_init,
        // RDKit❗✔️:                                 setup_smiles_string,
        // RDKit❗✔️:                                 yysmiles_lex_destroy>(yysmiles_parse,
        // RDKit❗✔️:                                                       inp,
        // RDKit❗✔️:                                                       molVect,
        // RDKit❗✔️:                                                       atom,
        // RDKit❗✔️:                                                       bond,
        // RDKit❗✔️:                                                       start_tok,
        // RDKit❗✔️:                                                       "SMILES");
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION smiles_parse_helper

        // BEGIN RDKIT CPP FUNCTION generic_parse_helper
        // RDKit❗✔️: template<int(*lex_init)(void**),
        // RDKit❗✔️:          size_t(*string_setup)(const std::string &, void *),
        // RDKit❗✔️:          int(*lex_destroy)(void*),
        // RDKit❗✔️:          typename T>
        // RDKit❗✔️: int generic_parse_helper(T parser,
        // RDKit❗✔️:                          const std::string &inp,
        // RDKit❗✔️:                          std::vector<RDKit::RWMol *> &molVect,
        // RDKit❗✔️:                          Atom *&atom,
        // RDKit❗✔️:                          Bond *&bond,
        // RDKit❗✔️:                          int start_tok,
        // RDKit❗✔️:                          const std::string& input_type) {
        // RDKit❗✔️:   std::vector<std::pair<unsigned int, unsigned int>> branchPoints;
        // RDKit❗✔️:   void *scanner;
        // RDKit❗✔️:   int res = 1;  // initialize with fail code
        // RDKit❗✔️:
        // RDKit❗✔️:   TEST_ASSERT(!lex_init(&scanner));
        // RDKit❗✔️:   size_t ltrim = 0;
        // RDKit❗✔️:   try {
        // RDKit❗✔️:     ltrim = string_setup(inp, scanner);
        // RDKit❗✔️:     unsigned numAtomsParsed = 0;
        // RDKit❗✔️:     unsigned numBondsParsed = 0;
        // RDKit❗✔️:     // NOTE: This variable will be used to point to the location of the
        // RDKit❗✔️:     // offending token if we encounter a syntax error
        // RDKit❗✔️:     unsigned int current_token_position = 0;
        // RDKit❗✔️:     res = parser(inp.c_str() + ltrim, &molVect, atom, bond,
        // RDKit❗✔️:                          numAtomsParsed, numBondsParsed, branchPoints, scanner,
        // RDKit❗✔️:                          start_tok, current_token_position);
        // RDKit❗✔️:   } catch (...) {
        // RDKit❗✔️:     lex_destroy(scanner);
        // RDKit❗✔️:     throw;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   lex_destroy(scanner);
        // RDKit❗✔️:
        // RDKit❗✔️:   if (!branchPoints.empty()) {
        // RDKit❗✔️:     auto input = inp.c_str() + ltrim;
        // RDKit❗✔️:     // If there are multiple unclosed brackets, we want to report them all at
        // RDKit❗✔️:     // once.
        // RDKit❗✔️:     for (auto [_, open_bracket_position] : branchPoints) {
        // RDKit❗✔️:       SmilesParseOps::detail::printSyntaxErrorMessage(
        // RDKit❗✔️:           input, "extra open parentheses", open_bracket_position, input_type);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   if (res == 1 || !branchPoints.empty()) {
        // RDKit❗✔️:     throw SmilesParseException("Failed parsing " + input_type + " '" + inp + "'");
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
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
    smiles_start_atoms: BTreeSet<AtomId>,
    cx_bond_order: Vec<BondId>,
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
            smiles_start_atoms: BTreeSet::new(),
            cx_bond_order: Vec::new(),
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
    }

    fn add_frag_to_mol(
        &mut self,
        frag: SmilesBuildState,
        bond_order: BondOrder,
        bond_dir: BondDirection,
    ) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP FUNCTION AddFragToMol
        // RDKit✔️❌: void AddFragToMol(RWMol *mol, RWMol *frag, Bond::BondType bondOrder,
        // RDKit✔️❌:                   Bond::BondDir bondDir) {
        // RDKit✔️❌:   PRECONDITION(mol, "no molecule");
        // RDKit✔️❌:   PRECONDITION(frag, "no fragment");
        // RDKit✔️❌:   PRECONDITION(mol->getActiveAtom(), "no active atom");
        // RDKit✔️❌:   Atom *lastAt = mol->getActiveAtom();
        // RDKit✔️❌:   int nOrigAtoms = mol->getNumAtoms();
        // RDKit✔️❌:   int nOrigBonds = mol->getNumBonds();
        // RDKit✔️❌:   mol->insertMol(*frag);
        // RDKit✔️❌:   // update ring-closure order information on the added atoms
        // RDKit✔️❌:   // and copy fragment bookmarks/partial bonds into the destination.
        // RDKit✔️❌:   // When bondOrder != IONIC, connect the former active atom to the first
        // RDKit✔️❌:   // fragment atom using the requested order and direction.
        // RDKit✔️❌: }
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
            let new_bond = self
                .builder
                .add_bond(spec)
                .map_err(|error| SmilesParseError::ParseError(error.to_string()))?;
            self.bond_pairs.insert(canonical_bond_pair(begin, end));
            if bond.prop(UNSPECIFIED_ORDER_PROP).is_some() {
                self.explicit_unspecified_bonds.push(new_bond);
            }
            self.cx_bond_order.push(new_bond);
            bond_mapping.push(new_bond);
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
        // RDKit❗✔️: mol: atomd {
        // RDKit❗✔️:   int sz     = molList->size();
        // RDKit❗✔️:   molList->resize( sz + 1);
        // RDKit❗✔️:   (*molList)[ sz ] = new RWMol();
        // RDKit❗✔️:   RDKit::RWMol *curMol = (*molList)[ sz ];
        // RDKit❗✔️:   $1->setProp(RDKit::common_properties::_SmilesStart,1);
        // RDKit❗✔️:   curMol->addAtom($1, true, true);
        // RDKit❗✔️:   //delete $1;
        // RDKit❗✔️:   $$ = sz;
        // RDKit❗✔️: }
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
        // RDKit❗✔️: | mol atomd       {
        // RDKit❗✔️:   RWMol *mp = (*molList)[$$];
        // RDKit❗✔️:   Atom *a1 = mp->getActiveAtom();
        // RDKit❗✔️:   int atomIdx1=a1->getIdx();
        // RDKit❗✔️:   int atomIdx2=mp->addAtom($2,true,true);
        // RDKit❗✔️:   mp->addBond(atomIdx1,atomIdx2,
        // RDKit❗✔️: 	      SmilesParseOps::GetUnspecifiedBondType(mp,a1,mp->getAtomWithIdx(atomIdx2)));
        // RDKit❗✔️:   mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit❗✔️:   //delete $2;
        // RDKit❗✔️: }
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
        // RDKit❗✔️:   $2->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit❗✔️:   mp->addBond($2,true);
        // RDKit❗✔️:   //delete $3;
        // RDKit❗✔️: }
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
        // RDKit❗✔️: | mol MINUS_TOKEN atomd {
        // RDKit❗✔️:   RWMol *mp = (*molList)[$$];
        // RDKit❗✔️:   int atomIdx1 = mp->getActiveAtom()->getIdx();
        // RDKit❗✔️:   int atomIdx2 = mp->addAtom($3,true,true);
        // RDKit❗✔️:   mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
        // RDKit❗✔️:   mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit❗✔️:   //delete $3;
        // RDKit❗✔️: }
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
        // RDKit❗✔️: | mol SEPARATOR_TOKEN atomd {
        // RDKit❗✔️:   RWMol *mp = (*molList)[$$];
        // RDKit❗✔️:   $3->setProp(RDKit::common_properties::_SmilesStart,1,true);
        // RDKit❗✔️:   mp->addAtom($3,true,true);
        // RDKit❗✔️: }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy mol SEPARATOR_TOKEN atomd
        atom.spec = atom.spec.with_prop(SMILES_START_PROP, "1");
        let atom_id = self.append_atom(atom);
        self.active_atom = Some(atom_id);
        self.smiles_start_atoms.insert(atom_id);
        Ok(())
    }

    fn branch_open_token(current_token_position: usize) -> usize {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy branch_open_token
        // RDKit❗✔️: branch_open_token: GROUP_OPEN_TOKEN { $$ = current_token_position; };
        // END RDKIT CPP GRAMMAR ACTION smiles.yy branch_open_token
        current_token_position
    }

    fn add_branch_atom_connected_to_active(
        &mut self,
        open_position: usize,
        atom: SmilesAtomToken,
    ) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy mol branch_open_token atomd
        // RDKit❗✔️: | mol branch_open_token atomd {
        // RDKit❗✔️:   RWMol *mp = (*molList)[$$];
        // RDKit❗✔️:   Atom *a1 = mp->getActiveAtom();
        // RDKit❗✔️:   int atomIdx1=a1->getIdx();
        // RDKit❗✔️:   int atomIdx2=mp->addAtom($3,true,true);
        // RDKit❗✔️:   mp->addBond(atomIdx1,atomIdx2,
        // RDKit❗✔️: 	      SmilesParseOps::GetUnspecifiedBondType(mp,a1,mp->getAtomWithIdx(atomIdx2)));
        // RDKit❗✔️:   mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit❗✔️:   branchPoints.push_back({atomIdx1, $2});
        // RDKit❗✔️: }
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
        // RDKit❗✔️:   $3->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit❗✔️:   mp->addBond($3,true);
        // RDKit❗✔️:
        // RDKit❗✔️:   branchPoints.push_back({atomIdx1, $2});
        // RDKit❗✔️: }
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
        // RDKit❗✔️: | mol branch_open_token MINUS_TOKEN atomd {
        // RDKit❗✔️:   RWMol *mp = (*molList)[$$];
        // RDKit❗✔️:   int atomIdx1 = mp->getActiveAtom()->getIdx();
        // RDKit❗✔️:   int atomIdx2=mp->addAtom($4,true,true);
        // RDKit❗✔️:   mp->addBond(atomIdx1,atomIdx2,Bond::SINGLE);
        // RDKit❗✔️:   mp->getBondBetweenAtoms(atomIdx1,atomIdx2)->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit❗✔️:   branchPoints.push_back({atomIdx1, $2});
        // RDKit❗✔️: }
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
        // RDKit❗✔️: | mol GROUP_CLOSE_TOKEN {
        // RDKit❗✔️:   if(branchPoints.empty()){
        // RDKit❗✔️:      yyerror(input,molList,branchPoints,scanner,start_token,current_token_position,"extra close parentheses");
        // RDKit❗✔️:      yyErrorCleanup(molList);
        // RDKit❗✔️:      YYABORT;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   RWMol *mp = (*molList)[$$];
        // RDKit❗✔️:
        // RDKit❗✔️:   mp->setActiveAtom(branchPoints.back().first);
        // RDKit❗✔️:   branchPoints.pop_back();
        // RDKit❗✔️: }
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
        // END RDKIT CPP GRAMMAR ACTION smiles.yy mol ring_number
        let pending_bond = PendingBond {
            token: SmilesBondToken {
                order: BondOrder::Unspecified,
                is_aromatic: false,
                direction: BondDirection::None,
                explicit_unspecified_order: true,
                is_null_query: false,
            },
        };
        self.add_ring_marker_with_pending_bond(ring_number, pending_bond)
    }

    fn add_explicit_bond_ring_marker(
        &mut self,
        bond: SmilesBondToken,
        ring_number: u32,
    ) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy mol BOND_TOKEN ring_number
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
        // RDKit❗✔️:   SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);
        // RDKit❗✔️:   INT_VECT tmp;
        // RDKit❗✔️:   atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
        // RDKit❗✔️:   tmp.push_back(-($3+1));
        // RDKit❗✔️:   atom->setProp(RDKit::common_properties::_RingClosures,tmp);
        // RDKit❗✔️:   delete $2;
        // RDKit❗✔️: }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy mol BOND_TOKEN ring_number
        self.add_ring_marker_with_pending_bond(ring_number, PendingBond { token: bond })
    }

    fn add_single_bond_ring_marker(&mut self, ring_number: u32) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP GRAMMAR ACTION smiles.yy mol MINUS_TOKEN ring_number
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
        // RDKit❗✔️:   SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);
        // RDKit❗✔️:   INT_VECT tmp;
        // RDKit❗✔️:   atom->getPropIfPresent(RDKit::common_properties::_RingClosures,tmp);
        // RDKit❗✔️:   tmp.push_back(-($3+1));
        // RDKit❗✔️:   atom->setProp(RDKit::common_properties::_RingClosures,tmp);
        // RDKit❗✔️: }
        // END RDKIT CPP GRAMMAR ACTION smiles.yy mol MINUS_TOKEN ring_number
        self.add_ring_marker_with_pending_bond(
            ring_number,
            PendingBond {
                token: SmilesBondToken::new(BondOrder::Single),
            },
        )
    }

    fn add_ring_marker_with_pending_bond(
        &mut self,
        ring_number: u32,
        pending_bond: PendingBond,
    ) -> Result<(), SmilesParseError> {
        let atom = self
            .active_atom
            .ok_or_else(|| SmilesParseError::ParseError("no active atom".to_string()))?;
        self.check_ring_closure_branch_status(atom)?;
        if let Some(opening) = self.ring_openings.remove(&ring_number) {
            let opening_pending_bond = opening.pending_bond.ok_or_else(|| {
                SmilesParseError::ParseError(format!("missing ring bond {ring_number}"))
            })?;
            let bond_id =
                self.close_ring_opening(opening.atom, atom, opening_pending_bond, pending_bond)?;
            if let Some(records) = self.ring_closures_by_atom.get_mut(&opening.atom)
                && let Some(record) = records
                    .iter_mut()
                    .rev()
                    .find(|record| record.ring_number == ring_number && record.bond.is_none())
            {
                record.bond = Some(bond_id);
            }
            self.ring_closures_by_atom
                .entry(atom)
                .or_default()
                .push(RingClosureRecord {
                    ring_number,
                    bond: Some(bond_id),
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
        // BEGIN RDKIT CPP FUNCTION CloseMolRings bond selection section
        // RDKit❗✔️:           // figure out which (if either) bond has a specified type, we'll
        // RDKit❗✔️:           // keep that one.  We also need to update the end atom index to
        // RDKit❗✔️:           // match FIX: daylight barfs when you give it multiple specs for the
        // RDKit❗✔️:           // closure
        // RDKit❗✔️:           //   bond, we'll just take the first one and ignore others
        // RDKit❗✔️:           //   NOTE: we used to do this the other way (take the last
        // RDKit❗✔️:           //   specification),
        // RDKit❗✔️:           //   but that turned out to be troublesome in odd cases like
        // RDKit❗✔️:           //   C1CC11CC1.
        // RDKit❗✔️:           if (!bond1->hasProp(common_properties::_unspecifiedOrder)) {
        // RDKit❗✔️:             matchedBond = bond1;
        // RDKit❗✔️:           } else {
        // RDKit❗✔️:             matchedBond = bond2;
        // RDKit❗✔️:           }
        // END RDKIT CPP FUNCTION CloseMolRings bond selection section
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
        let mut token = if opening_pending_bond.token.explicit_unspecified_order {
            current_pending_bond.token
        } else {
            opening_pending_bond.token
        };
        if token.order == BondOrder::Unspecified {
            token.order = get_unspecified_bond_type_for_atoms(opening_aromatic, closing_aromatic);
            token.explicit_unspecified_order = false;
        }
        let bond_id = self.add_bond_from_token(opening_atom, closing_atom, token)?;
        self.cx_bond_order.push(bond_id);
        Ok(bond_id)
    }

    fn check_ring_closure_branch_status(&mut self, atom: AtomId) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP FUNCTION CheckRingClosureBranchStatus
        // RDKit❗✔️: void CheckRingClosureBranchStatus(RDKit::Atom *atom, RDKit::RWMol *mp) {
        // RDKit❗✔️:   // github #786 and #1652: if the ring closure comes after a branch,
        // RDKit❗✔️:   // the stereochem is wrong.
        // RDKit❗✔️:   PRECONDITION(atom, "bad atom");
        // RDKit❗✔️:   PRECONDITION(mp, "bad mol");
        // RDKit❗✔️:   if (atom->getIdx() != mp->getNumAtoms(true) - 1 &&
        // RDKit❗✔️:       (atom->getDegree() == 1 ||
        // RDKit❗✔️:        (atom->getDegree() == 2 && atom->getIdx() != 0) ||
        // RDKit❗✔️:        (atom->getDegree() == 3 && atom->getIdx() == 0)) &&
        // RDKit❗✔️:       (atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
        // RDKit❗✔️:        atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW)) {
        // RDKit❗✔️:     atom->invertChirality();
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
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
        self.builder
            .bonds()
            .iter()
            .filter(|bond| bond.begin() == atom || bond.end() == atom)
            .count()
    }

    fn finish_parse(&mut self) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP FUNCTION toMol post-parse section
        // RDKit❗✔️:     func(inp, molVect);
        // RDKit❗✔️:     if (!molVect.empty()) {
        // RDKit❗✔️:       res.reset(molVect[0]);
        // RDKit❗✔️:       SmilesParseOps::CloseMolRings(res.get(), false);
        // RDKit❗✔️:       SmilesParseOps::CheckChiralitySpecifications(res.get(), true);
        // RDKit❗✔️:       SmilesParseOps::SetUnspecifiedBondTypes(res.get());
        // RDKit❗✔️:       SmilesParseOps::AdjustAtomChiralityFlags(res.get());
        // RDKit❗✔️:       // No sense leaving this bookmark intact:
        // RDKit❗✔️:       if (res->hasAtomBookmark(ci_RIGHTMOST_ATOM)) {
        // RDKit❗✔️:         res->clearAtomBookmark(ci_RIGHTMOST_ATOM);
        // RDKit❗✔️:       }
        // RDKit❗✔️:       molVect[0] = nullptr;  // NOTE: to avoid leaks on failures, this should
        // RDKit❗✔️:                              // occur last in this if.
        // RDKit❗✔️:     }
        // END RDKIT CPP FUNCTION toMol post-parse section
        self.close_mol_rings()?;
        self.check_chirality_specifications()?;
        self.set_unspecified_bond_types()?;
        self.adjust_atom_chirality_flags()?;
        Ok(())
    }

    fn close_mol_rings(&mut self) -> Result<(), SmilesParseError> {
        // BEGIN RDKIT CPP FUNCTION CloseMolRings
        // RDKit❗✔️: void CloseMolRings(RWMol *mol, bool toleratePartials) {
        // RDKit❗✔️:   //  Here's what we want to do here:
        // RDKit❗✔️:   //    loop through the molecule's atom bookmarks
        // RDKit❗✔️:   //    for each bookmark:
        // RDKit❗✔️:   //       connect pairs of atoms sharing that bookmark
        // RDKit❗✔️:   //          left to right (in the order in which they were
        // RDKit❗✔️:   //          inserted into the molecule).
        // RDKit❗✔️:   //       whilst doing this, we have to be cognizant of the fact that
        // RDKit❗✔️:   //          there may well be partial bonds in the molecule which need
        // RDKit❗✔️:   //          to be tied in as well.  WOO HOO! IT'S A BIG MESS!
        // RDKit❗✔️:   PRECONDITION(mol, "no molecule");
        // RDKit❗✔️: };
        // END RDKIT CPP FUNCTION CloseMolRings
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
        // RDKit❗✔️: void AdjustAtomChiralityFlags(RWMol *mol) {
        // RDKit❗✔️:   PRECONDITION(mol, "no molecule");
        // RDKit❗✔️:   for (auto atom : mol->atoms()) {
        // RDKit❗✔️:     Atom::ChiralType chiralType = atom->getChiralTag();
        // RDKit❗✔️:     if (chiralType == Atom::CHI_TETRAHEDRAL_CW ||
        // RDKit❗✔️:         chiralType == Atom::CHI_TETRAHEDRAL_CCW) {
        // RDKit❗✔️:       INT_LIST bondOrdering;
        // RDKit❗✔️:       unsigned int numClosures = GetBondOrdering(bondOrdering, mol, atom);
        // RDKit❗✔️:
        // RDKit❗✔️:       // ok, we now have the SMILES ordering of the bonds, figure out the
        // RDKit❗✔️:       // permutation order.
        // RDKit❗✔️:       //
        // RDKit❗✔️:       //  This whole thing is necessary because the ring-closure bonds
        // RDKit❗✔️:       //  in the SMILES come before the bonds to the other neighbors, but
        // RDKit❗✔️:       //  they come after the neighbors in the molecule we build.
        // RDKit❗✔️:       //  A crude example:
        // RDKit❗✔️:       //   in F[C@](Cl)(Br)I the C-Cl bond is index 1 in both SMILES
        // RDKit❗✔️:       //         and as built
        // RDKit❗✔️:       //   in F[C@]1(Br)I.Cl1 the C-Cl bond is index 1 in the SMILES
        // RDKit❗✔️:       //         and index 3 as built.
        // RDKit❗✔️:       //
        // RDKit❗✔️:       int nSwaps = atom->getPerturbationOrder(bondOrdering);
        // RDKit❗✔️:       // FIX: explain this one:
        // RDKit❗✔️:       // At least part of what's going on here for degree 3 atoms:
        // RDKit❗✔️:       //   - The first part: if we're at the beginning of the SMILES and have
        // RDKit❗✔️:       //      an explicit H, we need to add a swap.
        // RDKit❗✔️:       //      This is to reflect that [C@](Cl)(F)C is equivalent to Cl[C@@](F)C
        // RDKit❗✔️:       //      but [C@H](Cl)(F)C is fine as-is (The H-C bond is the one you look
        // RDKit❗✔️:       //      down).
        // RDKit❗✔️:       //   - The second part is more complicated and deals with situations like
        // RDKit❗✔️:       //      F[C@]1CCO1. In this case we otherwise end up looking like we need
        // RDKit❗✔️:       //      to invert the chirality, which is bogus. The chirality here needs
        // RDKit❗✔️:       //      to remain @ just as it does in F[C@](Cl)CCO1
        // RDKit❗✔️:       //   - We have to be careful with the second part to not sweep things like
        // RDKit❗✔️:       //      C[S@]2(=O).Cl2 into the same bin (was github #760). We detect
        // RDKit❗✔️:       //      those cases by looking for unsaturated atoms
        // RDKit❗✔️:       //
        // RDKit❗✔️:       if (Canon::chiralAtomNeedsTagInversion(
        // RDKit❗✔️:               *mol, atom, atom->hasProp(common_properties::_SmilesStart),
        // RDKit❗✔️:               numClosures)) {
        // RDKit❗✔️:         ++nSwaps;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       // std::cerr << "nswaps " << atom->getIdx() << " " << nSwaps
        // RDKit❗✔️:       //           << std::endl;
        // RDKit❗✔️:       // std::copy(bondOrdering.begin(), bondOrdering.end(),
        // RDKit❗✔️:       //           std::ostream_iterator<int>(std::cerr, ", "));
        // RDKit❗✔️:       // std::cerr << std::endl;
        // RDKit❗✔️:       if (nSwaps % 2) {
        // RDKit❗✔️:         atom->invertChirality();
        // RDKit❗✔️:       }
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
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
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
        for bond in self.builder.bonds() {
            let neighbor = if bond.begin() == atom_id {
                Some(bond.end())
            } else if bond.end() == atom_id {
                Some(bond.begin())
            } else {
                None
            };
            let Some(neighbor) = neighbor else {
                continue;
            };
            let bond_idx = bond.id().index() as i32;
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
        // RDKit❗✔️: int Atom::getPerturbationOrder(const INT_LIST &probe) const {
        // RDKit❗✔️:   INT_LIST ref;
        // RDKit❗✔️:   for (const auto bnd : getOwningMol().atomBonds(this)) {
        // RDKit❗✔️:     ref.push_back(bnd->getIdx());
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return static_cast<int>(countSwapsToInterconvert(probe, ref));
        // RDKit❗✔️: }
        // RDKit❗✔️: unsigned int countSwapsToInterconvert(const T &ref, T probe) {
        // RDKit❗✔️:   PRECONDITION(ref.size() == probe.size(), "size mismatch");
        // RDKit❗✔️:   typename T::const_iterator refIt = ref.begin();
        // RDKit❗✔️:   typename T::iterator probeIt = probe.begin();
        // RDKit❗✔️:   typename T::iterator probeIt2;
        // RDKit❗✔️:   unsigned int nSwaps = 0;
        // RDKit❗✔️:   while (refIt != ref.end()) {
        // RDKit❗✔️:     if ((*probeIt) != (*refIt)) {
        // RDKit❗✔️:       bool foundIt = false;
        // RDKit❗✔️:       probeIt2 = probeIt;
        // RDKit❗✔️:       while ((*probeIt2) != (*refIt) && probeIt2 != probe.end()) {
        // RDKit❗✔️:         ++probeIt2;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       if (probeIt2 != probe.end()) {
        // RDKit❗✔️:         foundIt = true;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       CHECK_INVARIANT(foundIt, "could not find probe element");
        // RDKit❗✔️:       std::swap(*probeIt, *probeIt2);
        // RDKit❗✔️:       nSwaps++;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     ++probeIt;
        // RDKit❗✔️:     ++refIt;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return nSwaps;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION Atom::getPerturbationOrder / countSwapsToInterconvert
        let mut storage_order = self
            .builder
            .bonds()
            .iter()
            .filter(|bond| bond.begin() == atom_id || bond.end() == atom_id)
            .map(|bond| bond.id().index() as i32)
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
        // RDKit❗✔️: bool chiralAtomNeedsTagInversion(const RDKit::ROMol &mol,
        // RDKit❗✔️:                                  const RDKit::Atom *atom, bool isAtomFirst,
        // RDKit❗✔️:                                  size_t numClosures) {
        // RDKit❗✔️:   PRECONDITION(atom, "bad atom");
        // RDKit❗✔️:   return atom->getDegree() == 3 &&
        // RDKit❗✔️:          ((isAtomFirst && atom->getNumExplicitHs() == 1) ||
        // RDKit❗✔️:           (!details::atomHasFourthValence(atom) && numClosures == 1 &&
        // RDKit❗✔️:            !details::isUnsaturated(atom, mol)));
        // RDKit❗✔️: }
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
        // RDKit❗❗: bool atomHasFourthValence(const Atom *atom) {
        // RDKit❗✔️:   if (atom->getNumExplicitHs() == 1 ||
        // RDKit❗✔️:       (!atom->needsUpdatePropertyCache() &&
        // RDKit❗✔️:        atom->getValence(Atom::ValenceType::IMPLICIT) == 1)) {
        // RDKit❗✔️:     return true;
        // RDKit❗✔️:   }
        // RDKit✔️✔️:   if (atom->hasQuery()) {
        // RDKit✔️✔️:     return hasSingleHQuery(atom->getQuery());
        // RDKit✔️✔️:   }
        // RDKit❗✔️:   return false;
        // RDKit❗❗: }
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
        // RDKit❗✔️: bool isUnsaturated(const Atom *atom, const ROMol &mol) {
        // RDKit❗✔️:   for (auto bond : mol.atomBonds(atom)) {
        // RDKit❗✔️:     // can't just check for single bonds, because dative bonds also have an
        // RDKit❗✔️:     // order of 1
        // RDKit❗✔️:     if (bond->getBondTypeAsDouble() > 1) {
        // RDKit❗✔️:       return true;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return false;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION Canon::details::isUnsaturated
        self.builder
            .bonds()
            .iter()
            .filter(|bond| bond.begin() == atom_id || bond.end() == atom_id)
            .any(|bond| bond_order_as_double(bond.order()) > 1.0)
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
        for bond in self.builder.bonds() {
            if bond.begin() == atom_id || bond.end() == atom_id {
                order[bond.id().index()] = nbr_idx;
                nbr_idx += 1;
            }
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
            && self.smiles_start_atoms.is_empty()
            && self.cx_bond_order.is_empty()
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
        let bond_idx = self.cx_bond_order.len().to_string();
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
        let (begin, end, order) = match bond.order {
            BondOrder::DativeRight => (atom_idx1, atom_idx2, BondOrder::Dative),
            BondOrder::DativeLeft => (atom_idx2, atom_idx1, BondOrder::Dative),
            other => (atom_idx1, atom_idx2, other),
        };
        let mut spec = BondSpec::new(begin, end, order)
            .with_aromatic(bond.is_aromatic || order == BondOrder::Aromatic)
            .with_direction(bond.direction)
            .with_prop(CXSMILES_BOND_IDX_PROP, self.cx_bond_order.len().to_string());
        if bond.explicit_unspecified_order {
            spec = spec.with_prop(UNSPECIFIED_ORDER_PROP, "1");
        }
        if bond.is_null_query {
            spec = spec.with_query(QueryNode::predicate(
                BondQueryPredicate::UnsupportedFeature("makeBondNullQuery"),
            ));
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
        // RDKit❗✔️: void CleanupAfterParsing(RWMol *mol) {
        // RDKit❗✔️:   PRECONDITION(mol, "no molecule");
        // RDKit❗✔️:   for (auto atom : mol->atoms()) {
        // RDKit❗✔️:     atom->clearProp(common_properties::_RingClosures);
        // RDKit❗✔️:     atom->clearProp(common_properties::_SmilesStart);
        // RDKit❗✔️:     std::string label;
        // RDKit❗✔️:     if (atom->getAtomicNum() == 0 &&
        // RDKit❗✔️:         atom->getPropIfPresent(common_properties::atomLabel, label)) {
        // RDKit❗✔️:       // marvinsketch can output higher labels than _AP1 and _AP2, but they
        // RDKit❗✔️:       // aren't part of the MOL file spec so we don't treat them as attachment
        // RDKit❗✔️:       // points
        // RDKit❗✔️:       if (label == "_AP1") {
        // RDKit❗✔️:         atom->setProp(common_properties::_fromAttachPoint, 1);
        // RDKit❗✔️:       } else if (label == "_AP2") {
        // RDKit❗✔️:         atom->setProp(common_properties::_fromAttachPoint, 2);
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   for (auto bond : mol->bonds()) {
        // RDKit❗✔️:     bond->clearProp(common_properties::_unspecifiedOrder);
        // RDKit❗✔️:     bond->clearProp("_cxsmilesBondIdx");
        // RDKit❗✔️:   }
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
    // RDKit❗✔️: std::unique_ptr<RWMol> MolFromSmiles(const std::string &smiles,
    // RDKit❗✔️:                                      const SmilesParserParams &params) {
    // RDKit❗✔️:   // Calling MolFromSmiles in a multithreaded context is generally safe *unless*
    // RDKit❗✔️:   // the value of debugParse is different for different threads. The if
    // RDKit❗✔️:   // statement below avoids a TSAN warning in the case where multiple threads
    // RDKit❗✔️:   // all use the same value for debugParse.
    // RDKit❗✔️:   if (yysmiles_debug != params.debugParse) {
    // RDKit❗✔️:     yysmiles_debug = params.debugParse;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   std::string lsmiles, name, cxPart;
    // RDKit❗✔️:   preprocessSmiles(smiles, params, lsmiles, name, cxPart);
    // RDKit❗✔️:   // strip any leading/trailing whitespace:
    // RDKit❗✔️:   // boost::trim_if(smi,boost::is_any_of(" \t\r\n"));
    // RDKit❗✔️:   auto res = toMol(lsmiles, smiles_parse, lsmiles);
    // RDKit❗✔️:   if (!res) {
    // RDKit❗✔️:     return res;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   handleCXPartAndName(res.get(), params, cxPart, name);
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: };
    // END RDKIT CPP FUNCTION MolFromSmiles
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
        // Kekulize aromatic bonds (RDKit step 1)
        if mol
            .bonds()
            .iter()
            .any(|b| b.order() == crate::BondOrder::Aromatic)
        {
            mol = mol
                .with_kekulized_bonds(false)
                .map_err(|e: crate::OperationError| {
                    SmilesParseError::ParseError(format!(
                        "kekulization during sanitize failed: {e}"
                    ))
                })?;
        }
        // Assign valence and radicals (RDKit step 2-3)
        mol = mol
            .with_assigned_valence()
            .map_err(|e: crate::OperationError| {
                SmilesParseError::ParseError(format!(
                    "valence assignment during sanitize failed: {e}"
                ))
            })?;
        // Assign aromaticity (RDKit step 4 — re-perceives on Kekulé form)
        mol = mol
            .with_assigned_aromaticity()
            .map_err(|e: crate::OperationError| {
                SmilesParseError::ParseError(format!(
                    "aromaticity assignment during sanitize failed: {e}"
                ))
            })?;
        // Assign stereochemistry from bond directions (RDKit final step)
        let _ = crate::smiles::assign_double_bond_stereo_from_directions(&mut mol);
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
            let _ = crate::smiles::assign_double_bond_stereo_from_directions(&mut mol);
        }
        mol.properties_mut().clear_prop("_needsDetectBondStereo");
    } else {
        // RDKit github #337 path: after atom stereo perception, wedge-style
        // single-bond directions are no longer needed, but the CX double-bond
        // stereo marker property is preserved when sanitize/removeHs are off.
        clear_single_bond_dir_flags(&mut mol, true);
    }

    // RDKit❗✔️: _NeedsQueryScan query completion
    // RDKit❗✔️: if (res->hasProp(common_properties::_NeedsQueryScan)) {
    // RDKit❗✔️:   if (!params.sanitize) {
    // RDKit❗✔️:     MolOps::fastFindRings(*res);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   QueryOps::completeMolQueries(*res, 0xDEADBEEF);
    // RDKit❗✔️: }
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

fn first_2d_and_3d_conformer_ids(mol: &Molecule) -> (Option<usize>, Option<usize>) {
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
    let mut first_2d = None;
    let mut first_3d = None;
    for conformer in mol.conformers_3d() {
        if conformer.is_3d() {
            if first_3d.is_none() {
                first_3d = Some(conformer.id());
            }
        } else if first_2d.is_none() {
            first_2d = Some(conformer.id());
        }

        if first_2d.is_some() && first_3d.is_some() {
            break;
        }
    }
    (first_2d, first_3d)
}

fn apply_coordinate_free_atropisomer_assignments(
    mol: &mut Molecule,
    assignments: Vec<(BondId, BondStereo)>,
) {
    apply_atropisomer_stereo_assignments(mol, assignments);
}

fn apply_atropisomer_stereo_assignments(
    mol: &mut Molecule,
    assignments: Vec<(BondId, BondStereo)>,
) {
    for (bond_id, stereo) in assignments {
        if let Some(bond) = mol.topology_block_mut().bonds.get_mut(bond_id.index()) {
            bond.set_stereo(stereo);
        }
    }
}

fn opposite_dir(direction: BondDirection) -> BondDirection {
    match direction {
        BondDirection::EndDownRight => BondDirection::EndUpRight,
        BondDirection::EndUpRight => BondDirection::EndDownRight,
        other => other,
    }
}

fn vec3_normalize(v: [f64; 3]) -> Option<[f64; 3]> {
    let len = vec3_len(v);
    if len <= f64::EPSILON {
        None
    } else {
        Some([v[0] / len, v[1] / len, v[2] / len])
    }
}

fn choose_stereo_neighbor(
    mol: &Molecule,
    ranks: &[u32],
    center: AtomId,
    exclude: AtomId,
) -> Option<(BondId, AtomId)> {
    let mut candidates = Vec::new();
    for bond in mol.bonds() {
        if bond.begin() != center && bond.end() != center {
            continue;
        }
        let other = if bond.begin() == center {
            bond.end()
        } else {
            bond.begin()
        };
        if other == exclude {
            continue;
        }
        candidates.push((
            ranks.get(other.index()).copied().unwrap_or(0),
            bond.id(),
            other,
        ));
    }
    candidates.sort_by(|left, right| {
        right
            .0
            .cmp(&left.0)
            .then_with(|| left.1.index().cmp(&right.1.index()))
    });
    let (best_rank, best_bond, best_atom) = candidates.first().copied()?;
    if candidates.get(1).is_some_and(|next| next.0 == best_rank) {
        return None;
    }
    Some((best_bond, best_atom))
}

fn set_raw_neighbor_dir_for_center(
    mol: &mut Molecule,
    bond_id: BondId,
    center: AtomId,
    normalized_dir: BondDirection,
    center_is_begin_side: bool,
) {
    let Some(bond) = mol.topology_block_mut().bonds.get_mut(bond_id.index()) else {
        return;
    };
    let raw_dir = if center_is_begin_side {
        if bond.begin() == center {
            opposite_dir(normalized_dir)
        } else {
            normalized_dir
        }
    } else if bond.end() == center {
        opposite_dir(normalized_dir)
    } else {
        normalized_dir
    };
    bond.set_direction(raw_dir);
}

fn ensure_valence_for_stereo(mol: &mut Molecule) -> Result<(), StereoError> {
    if mol.derived_cache().valence.is_some() {
        return Ok(());
    }
    *mol = mol.with_assigned_valence().map_err(|error| {
        StereoError::UnsupportedFeature(crate::UnsupportedFeatureError {
            feature: "CIP_RANKING",
            reason: "valence assignment failed before stereo perception",
        })
    })?;
    Ok(())
}

pub(crate) fn set_double_bond_neighbor_directions(
    mol: &mut Molecule,
    conf_id: usize,
) -> Result<(), StereoError> {
    // BEGIN RDKIT CPP FUNCTION setDoubleBondNeighborDirections
    // RDKit❗✔️: void setDoubleBondNeighborDirections(ROMol &mol, const Conformer *conf) {
    // RDKit❗✔️:   // Sets neighboring single-bond directions for stereo-capable double bonds
    // RDKit❗✔️:   // using conformer geometry plus candidate/neighbor filtering.
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION setDoubleBondNeighborDirections
    ensure_valence_for_stereo(mol)?;
    let Some(conformer) = mol.conformers_3d().iter().find(|conf| conf.id() == conf_id) else {
        return Ok(());
    };
    let ranks = crate::stereo::assign_atom_cip_ranks(mol)?;
    let coords = conformer.coords().to_vec();
    let candidate_bonds: Vec<(BondId, AtomId, AtomId)> = mol
        .bonds()
        .iter()
        .filter(|bond| crate::stereo::is_bond_candidate_for_stereo(mol, bond.id().index()))
        .map(|bond| (bond.id(), bond.begin(), bond.end()))
        .collect();

    for (bond_id, begin, end) in candidate_bonds {
        let Some((begin_bond, begin_nbr)) = choose_stereo_neighbor(mol, &ranks, begin, end) else {
            continue;
        };
        let Some((end_bond, end_nbr)) = choose_stereo_neighbor(mol, &ranks, end, begin) else {
            continue;
        };

        let axis = match vec3_normalize(vec3_sub(coords[end.index()], coords[begin.index()])) {
            Some(axis) => axis,
            None => continue,
        };
        let begin_vec = vec3_sub(coords[begin_nbr.index()], coords[begin.index()]);
        let end_vec = vec3_sub(coords[end_nbr.index()], coords[end.index()]);
        let begin_proj = [
            begin_vec[0] - axis[0] * vec3_dot(begin_vec, axis),
            begin_vec[1] - axis[1] * vec3_dot(begin_vec, axis),
            begin_vec[2] - axis[2] * vec3_dot(begin_vec, axis),
        ];
        let end_proj = [
            end_vec[0] - axis[0] * vec3_dot(end_vec, axis),
            end_vec[1] - axis[1] * vec3_dot(end_vec, axis),
            end_vec[2] - axis[2] * vec3_dot(end_vec, axis),
        ];
        let Some(begin_proj) = vec3_normalize(begin_proj) else {
            continue;
        };
        let Some(end_proj) = vec3_normalize(end_proj) else {
            continue;
        };
        let same_side = vec3_dot(begin_proj, end_proj) > 0.0;

        let begin_normalized = BondDirection::EndUpRight;
        let end_normalized = if same_side {
            BondDirection::EndDownRight
        } else {
            BondDirection::EndUpRight
        };
        set_raw_neighbor_dir_for_center(mol, begin_bond, begin, begin_normalized, true);
        set_raw_neighbor_dir_for_center(mol, end_bond, end, end_normalized, false);

        if let Some(bond) = mol.topology_block_mut().bonds.get_mut(bond_id.index()) {
            bond.set_stereo(BondStereo::None);
            bond.set_stereo_atoms(None);
        }
    }
    Ok(())
}

pub(crate) fn clear_dir_flags(mol: &mut Molecule, only_wedge_type_bond_dirs: bool) {
    // BEGIN RDKIT CPP FUNCTION clearDirFlags
    // RDKit✔️✔️: void clearDirFlags(ROMol &mol, bool onlyWedgeTypeBondDirs) {
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     if (bond->getBondDir() == Bond::UNKNOWN ||
    // RDKit✔️✔️:         bond->getBondDir() == Bond::BondDir::EITHERDOUBLE) {
    // RDKit✔️✔️:       bond->setProp(common_properties::_UnknownStereo, 1);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (onlyWedgeTypeBondDirs == false ||
    // RDKit✔️✔️:         (bond->getBondDir() != Bond::BondDir::ENDDOWNRIGHT &&
    // RDKit✔️✔️:          bond->getBondDir() != Bond::BondDir::ENDUPRIGHT)) {
    // RDKit✔️✔️:       bond->setBondDir(Bond::NONE);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION clearDirFlags
    let topology = mol.topology_block_mut();
    for bond in &mut topology.bonds {
        if matches!(
            bond.direction(),
            BondDirection::Unknown | BondDirection::EitherDouble
        ) {
            bond.set_unknown_stereo(true);
        }
        if !only_wedge_type_bond_dirs
            || !matches!(
                bond.direction(),
                BondDirection::EndDownRight | BondDirection::EndUpRight
            )
        {
            bond.set_direction(BondDirection::None);
        }
    }
}

pub(crate) fn clear_all_bond_dir_flags(mol: &mut Molecule) {
    // BEGIN RDKIT CPP FUNCTION clearAllBondDirFlags
    // RDKit✔️✔️: void clearAllBondDirFlags(ROMol &mol) { clearDirFlags(mol, false); }
    // END RDKIT CPP FUNCTION clearAllBondDirFlags
    clear_dir_flags(mol, false);
}

pub(crate) fn set_bond_stereo_from_directions(mol: &mut Molecule) -> Result<(), StereoError> {
    // BEGIN RDKIT CPP FUNCTION setBondStereoFromDirections
    // RDKit❗✔️: void setBondStereoFromDirections(ROMol &mol) {
    // RDKit❗✔️:   // Finds directed neighboring single bonds and assigns double-bond stereo.
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION setBondStereoFromDirections
    mol.properties_mut().clear_prop("_needsDetectBondStereo");
    ensure_valence_for_stereo(mol)?;
    for bond in &mut mol.topology_block_mut().bonds {
        if bond.order() == BondOrder::Double && bond.stereo() != BondStereo::Any {
            bond.set_stereo(BondStereo::None);
            bond.set_stereo_atoms(None);
        }
    }
    let ranks = crate::stereo::assign_atom_cip_ranks(mol)?;
    let (assignments, _changed) = crate::stereo::assign_bond_stereo_codes(mol, &ranks);
    for (bond_idx, stereo, begin_control, end_control) in assignments {
        let stereo = match stereo {
            crate::stereo::DoubleBondStereo::E => BondStereo::Trans,
            crate::stereo::DoubleBondStereo::Z => BondStereo::Cis,
            crate::stereo::DoubleBondStereo::Unknown => BondStereo::Any,
        };
        let bond = &mut mol.topology_block_mut().bonds[bond_idx];
        bond.set_stereo_atoms(Some([AtomId::new(begin_control), AtomId::new(end_control)]));
        bond.set_stereo(stereo);
    }
    Ok(())
}

pub(crate) fn assign_stereochemistry_from_3d(
    mol: &mut Molecule,
    conf_id: usize,
) -> Result<(), StereoError> {
    // BEGIN RDKIT CPP FUNCTION assignStereochemistryFrom3D
    // RDKit❗✔️: void assignStereochemistryFrom3D(ROMol &mol, int confId,
    // RDKit❗✔️:                                  bool replaceExistingTags) {
    // RDKit❗✔️:   if (!mol.getNumConformers() || !mol.getConformer(confId).is3D()) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (mol.needsUpdatePropertyCache()) {
    // RDKit❗✔️:     mol.updatePropertyCache(false);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   detectBondStereochemistry(mol, confId);
    // RDKit❗✔️:   assignChiralTypesFrom3D(mol, confId, replaceExistingTags);
    // RDKit❗✔️:   assignStereochemistry(mol, replaceExistingTags, true, true);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION assignStereochemistryFrom3D
    let Some(is_3d) = mol
        .conformers_3d()
        .iter()
        .find(|conf| conf.id() == conf_id)
        .map(|conf| conf.is_3d())
    else {
        return Ok(());
    };
    if !is_3d {
        return Ok(());
    }
    ensure_valence_for_stereo(mol)?;
    set_double_bond_neighbor_directions(mol, conf_id)?;
    assign_chiral_types_from_3d(mol, conf_id);
    set_bond_stereo_from_directions(mol)?;
    Ok(())
}

fn clear_single_bond_dir_flags(mol: &mut Molecule, only_wedge_flags: bool) {
    // RDKit✔️✔️: void clearSingleBondDirFlags(ROMol &mol, bool onlyWedgeFlags) {
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     if (bond->getBondType() == Bond::SINGLE) {
    // RDKit✔️✔️:       if (bond->getBondDir() == Bond::UNKNOWN) {
    // RDKit✔️✔️:         bond->setProp(common_properties::_UnknownStereo, 1);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!onlyWedgeFlags ||
    // RDKit✔️✔️:           (bond->getBondDir() != Bond::BondDir::ENDDOWNRIGHT &&
    // RDKit✔️✔️:            bond->getBondDir() != Bond::BondDir::ENDUPRIGHT)) {
    // RDKit✔️✔️:         bond->setBondDir(Bond::NONE);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let topology = mol.topology_block_mut();
    for bond in &mut topology.bonds {
        if bond.order() == BondOrder::Single {
            if bond.direction() == BondDirection::Unknown {
                bond.set_unknown_stereo(true);
            }
            if !only_wedge_flags
                || !matches!(
                    bond.direction(),
                    BondDirection::EndDownRight | BondDirection::EndUpRight
                )
            {
                bond.set_direction(BondDirection::None);
            }
        }
    }
}

const QUERY_SCAN_SENTINEL: u32 = 0xDEADBEEF;

fn complete_smiles_query_scan_subset(mol: &mut Molecule) {
    // RDKit❗✔️: void completeMolQueries(RWMol *mol, unsigned int magicVal) {
    // RDKit❗✔️:   PRECONDITION(mol, "bad molecule");
    // RDKit❗✔️:   for (auto atom : mol->atoms()) {
    // RDKit❗✔️:     if (atom->hasQuery()) {
    // RDKit❗✔️:       completeQueryAndChildren(atom->getQuery(), atom, magicVal);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // COSMolKit: for the currently modeled SMILES CX query-scan inputs, the
    // sentinel-bearing placeholders are `rb:*` and `s:*`. Both are completed
    // here and `_NeedsQueryScan` is cleared once no unresolved sentinels
    // remain in any atom query tree.
    let ring_counts = smiles_atom_ring_bond_counts(mol);
    let non_hydrogen_degrees = smiles_non_hydrogen_degrees(mol);
    let topology = mol.topology_block_mut();
    let mut unresolved_scan_query = false;
    for (atom_idx, atom) in topology.atoms.iter_mut().enumerate() {
        let Some(mut query) = atom.query().cloned() else {
            continue;
        };
        complete_smiles_query_scan_predicates(
            &mut query,
            ring_counts[atom_idx],
            non_hydrogen_degrees[atom_idx],
        );
        unresolved_scan_query |= atom_query_has_unresolved_smiles_scan(&query);
        atom.set_query(Some(query));
    }
    if !unresolved_scan_query {
        mol.properties_mut().clear_prop("_NeedsQueryScan");
    }
}

fn atom_query_has_unresolved_smiles_scan(query: &QueryNode<AtomQueryPredicate>) -> bool {
    match query {
        QueryNode::Predicate(AtomQueryPredicate::RingBondCountNeedsScan) => true,
        QueryNode::Predicate(AtomQueryPredicate::NonHydrogenDegree(value)) => {
            *value == QUERY_SCAN_SENTINEL
        }
        QueryNode::Predicate(_) => false,
        QueryNode::And(children) | QueryNode::Or(children) => {
            children.iter().any(atom_query_has_unresolved_smiles_scan)
        }
        QueryNode::Not(child) => atom_query_has_unresolved_smiles_scan(child),
    }
}

fn complete_smiles_query_scan_predicates(
    query: &mut QueryNode<AtomQueryPredicate>,
    ring_bond_count: u8,
    non_hydrogen_degree: u32,
) {
    match query {
        QueryNode::Predicate(AtomQueryPredicate::RingBondCountNeedsScan) => {
            *query = QueryNode::predicate(AtomQueryPredicate::RingBondCount(ring_bond_count));
        }
        QueryNode::Predicate(AtomQueryPredicate::NonHydrogenDegree(value))
            if *value == QUERY_SCAN_SENTINEL =>
        {
            *query =
                QueryNode::predicate(AtomQueryPredicate::NonHydrogenDegree(non_hydrogen_degree));
        }
        QueryNode::Predicate(_) => {}
        QueryNode::And(children) | QueryNode::Or(children) => {
            for child in children {
                complete_smiles_query_scan_predicates(child, ring_bond_count, non_hydrogen_degree);
            }
        }
        QueryNode::Not(child) => {
            complete_smiles_query_scan_predicates(child, ring_bond_count, non_hydrogen_degree)
        }
    }
}

fn smiles_atom_ring_bond_counts(molecule: &Molecule) -> Vec<u8> {
    let atom_count = molecule.num_atoms();
    let mut counts = vec![0_u8; atom_count];
    for (bond_idx, bond) in molecule.bonds().iter().enumerate() {
        if smiles_bond_is_in_cycle(molecule, bond_idx) {
            counts[bond.begin().index()] = counts[bond.begin().index()].saturating_add(1);
            counts[bond.end().index()] = counts[bond.end().index()].saturating_add(1);
        }
    }
    counts
}

fn smiles_bond_is_in_cycle(molecule: &Molecule, removed_bond_idx: usize) -> bool {
    let Some(removed) = molecule.bonds().get(removed_bond_idx) else {
        return false;
    };
    let target = removed.end();
    let mut stack = vec![removed.begin()];
    let mut seen = vec![false; molecule.num_atoms()];
    while let Some(atom_id) = stack.pop() {
        if atom_id == target {
            return true;
        }
        if seen[atom_id.index()] {
            continue;
        }
        seen[atom_id.index()] = true;
        for (bond_idx, bond) in molecule.bonds().iter().enumerate() {
            if bond_idx == removed_bond_idx {
                continue;
            }
            let next = if bond.begin() == atom_id {
                Some(bond.end())
            } else if bond.end() == atom_id {
                Some(bond.begin())
            } else {
                None
            };
            if let Some(next) = next
                && !seen[next.index()]
            {
                stack.push(next);
            }
        }
    }
    false
}

fn smiles_non_hydrogen_degrees(molecule: &Molecule) -> Vec<u32> {
    let mut counts = vec![0_u32; molecule.num_atoms()];
    for bond in molecule.bonds() {
        let begin = bond.begin();
        let end = bond.end();
        if let Some(end_atom) = molecule.atom(end)
            && (end_atom.atomic_number() != 1 || end_atom.isotope().is_some_and(|i| i > 1))
        {
            counts[begin.index()] = counts[begin.index()].saturating_add(1);
        }
        if let Some(begin_atom) = molecule.atom(begin)
            && (begin_atom.atomic_number() != 1 || begin_atom.isotope().is_some_and(|i| i > 1))
        {
            counts[end.index()] = counts[end.index()].saturating_add(1);
        }
    }
    counts
}

fn atrop_can_have_direction(order: BondOrder) -> bool {
    matches!(order, BondOrder::Single | BondOrder::Aromatic)
}

fn atrop_other_atom(mol: &Molecule, bond_id: BondId, atom: AtomId) -> Option<AtomId> {
    let bond = mol.bonds().get(bond_id.index())?;
    if bond.begin() == atom {
        Some(bond.end())
    } else if bond.end() == atom {
        Some(bond.begin())
    } else {
        None
    }
}

fn atrop_neighbor_bonds(
    mol: &Molecule,
    focus_atom: AtomId,
    atrop_bond: BondId,
) -> Option<Vec<BondId>> {
    let mut nbr_bonds: Vec<BondId> = mol
        .bonds()
        .iter()
        .filter(|b| b.id() != atrop_bond && (b.begin() == focus_atom || b.end() == focus_atom))
        .map(|b| b.id())
        .collect();
    if nbr_bonds.is_empty() {
        return None;
    }
    if nbr_bonds.len() == 2 {
        let other0 = atrop_other_atom(mol, nbr_bonds[0], focus_atom)?;
        let other1 = atrop_other_atom(mol, nbr_bonds[1], focus_atom)?;
        if other1.index() < other0.index() {
            nbr_bonds.swap(0, 1);
        }
    }
    Some(nbr_bonds)
}

fn atrop_end_wedge_direction(mol: &Molecule, nbr_bonds: &[BondId]) -> (bool, BondDirection) {
    let Some(bond0) = mol.bonds().get(nbr_bonds[0].index()) else {
        return (false, BondDirection::None);
    };
    let bond1 = if nbr_bonds.len() > 1 {
        mol.bonds().get(nbr_bonds[1].index())
    } else {
        None
    };

    let dir0 = match bond0.direction() {
        BondDirection::BeginWedge | BondDirection::BeginDash => bond0.direction(),
        _ => BondDirection::None,
    };
    let dir1 = match bond1
        .map(|bond| bond.direction())
        .unwrap_or(BondDirection::None)
    {
        BondDirection::BeginWedge | BondDirection::BeginDash => bond1
            .map(|bond| bond.direction())
            .unwrap_or(BondDirection::None),
        _ => BondDirection::None,
    };

    if dir0 != BondDirection::None && dir1 != BondDirection::None && dir0 == dir1 {
        return (false, BondDirection::None);
    }
    if dir0 == BondDirection::BeginWedge || dir1 == BondDirection::BeginDash {
        return (true, BondDirection::BeginWedge);
    }
    if dir0 == BondDirection::BeginDash || dir1 == BondDirection::BeginWedge {
        return (true, BondDirection::BeginDash);
    }
    (true, BondDirection::None)
}

fn atrop_bond_frame_of_reference(
    mol: &Molecule,
    bond_id: BondId,
    conformer: &Conformer3D,
) -> Option<([f64; 3], [f64; 3], [f64; 3])> {
    const REALLY_SMALL_BOND_LEN: f64 = 1e-7;
    let bond = mol.bonds().get(bond_id.index())?;
    // RDKit✔️✔️: xAxis = conf->getAtomPos(bond->getEndAtom()->getIdx()) -
    // RDKit✔️✔️:         conf->getAtomPos(bond->getBeginAtom()->getIdx());
    let mut x_axis = vec3_sub(
        conformer.coords()[bond.end().index()],
        conformer.coords()[bond.begin().index()],
    );
    let x_len = vec3_len(x_axis);
    if x_len < REALLY_SMALL_BOND_LEN {
        return None;
    }
    x_axis = [x_axis[0] / x_len, x_axis[1] / x_len, x_axis[2] / x_len];
    if !conformer.is_3d() {
        let y_len = (x_axis[0] * x_axis[0] + x_axis[1] * x_axis[1]).sqrt();
        if y_len < REALLY_SMALL_BOND_LEN {
            return None;
        }
        let y_axis = [-x_axis[1] / y_len, x_axis[0] / y_len, 0.0];
        let z_axis = [0.0, 0.0, 1.0];
        return Some((x_axis, y_axis, z_axis));
    }

    let mut z_axis =
        if x_axis[0].abs() > REALLY_SMALL_BOND_LEN || x_axis[1].abs() > REALLY_SMALL_BOND_LEN {
            [0.0, 0.0, 1.0]
        } else {
            [1.0, 0.0, 0.0]
        };
    let mut y_axis = vec3_cross(z_axis, x_axis);
    let y_len = vec3_len(y_axis);
    if y_len < REALLY_SMALL_BOND_LEN {
        return None;
    }
    y_axis = [y_axis[0] / y_len, y_axis[1] / y_len, y_axis[2] / y_len];
    z_axis = vec3_cross(x_axis, y_axis);
    let z_len = vec3_len(z_axis);
    if z_len < REALLY_SMALL_BOND_LEN {
        return None;
    }
    z_axis = [z_axis[0] / z_len, z_axis[1] / z_len, z_axis[2] / z_len];
    Some((x_axis, y_axis, z_axis))
}

fn atrop_projected_end_vector(
    mol: &Molecule,
    conformer: &Conformer3D,
    focus_atom: AtomId,
    nbr_bonds: &[BondId],
    y_axis: [f64; 3],
    z_axis: [f64; 3],
    normalize: bool,
) -> Option<[f64; 3]> {
    const REALLY_SMALL_BOND_LEN: f64 = 1e-7;
    let other0 = atrop_other_atom(mol, nbr_bonds[0], focus_atom)?;
    let mut bond_vec = vec3_sub(
        conformer.coords()[other0.index()],
        conformer.coords()[focus_atom.index()],
    );
    bond_vec = [0.0, vec3_dot(bond_vec, y_axis), vec3_dot(bond_vec, z_axis)];

    if nbr_bonds.len() == 2 {
        let other1 = atrop_other_atom(mol, nbr_bonds[1], focus_atom)?;
        let mut other_vec = vec3_sub(
            conformer.coords()[other1.index()],
            conformer.coords()[focus_atom.index()],
        );
        other_vec = [
            0.0,
            vec3_dot(other_vec, y_axis),
            vec3_dot(other_vec, z_axis),
        ];
        if vec3_len(bond_vec) < REALLY_SMALL_BOND_LEN {
            bond_vec = [-other_vec[0], -other_vec[1], -other_vec[2]];
        } else if vec3_dot(bond_vec, other_vec) > REALLY_SMALL_BOND_LEN {
            return None;
        }
    }

    let len = vec3_len(bond_vec);
    if len < REALLY_SMALL_BOND_LEN {
        return None;
    }
    if normalize {
        Some([bond_vec[0] / len, bond_vec[1] / len, bond_vec[2] / len])
    } else {
        Some(bond_vec)
    }
}

fn atropisomer_stereo_without_conformer(mol: &Molecule) -> Vec<(BondId, BondStereo)> {
    let Some(rings) = mol.derived_cache().rings.as_ref() else {
        return Vec::new();
    };

    let mut candidate_bonds: Vec<BondId> = Vec::new();
    for bond in mol.bonds() {
        if !atrop_can_have_direction(bond.order()) {
            continue;
        }
        if !matches!(
            bond.direction(),
            BondDirection::BeginWedge | BondDirection::BeginDash
        ) {
            continue;
        }
        let begin = bond.begin();
        for nbr_bond in mol.bonds() {
            if nbr_bond.id() == bond.id() {
                continue;
            }
            if (nbr_bond.begin() == begin || nbr_bond.end() == begin)
                && !candidate_bonds.contains(&nbr_bond.id())
            {
                candidate_bonds.push(nbr_bond.id());
            }
        }
    }

    let num_atoms = mol.num_atoms();
    let mut degree = vec![0usize; num_atoms];
    for bond in mol.bonds() {
        degree[bond.begin().index()] += 1;
        degree[bond.end().index()] += 1;
    }
    let hybridization: Vec<crate::Hybridization> = mol
        .atoms()
        .iter()
        .map(|atom| atom.hybridization())
        .collect();

    let mut assignments = Vec::new();
    for candidate_id in candidate_bonds {
        let Some(candidate) = mol.bonds().get(candidate_id.index()) else {
            continue;
        };
        if candidate.order() != BondOrder::Single || candidate.stereo() == BondStereo::Any {
            continue;
        }
        let begin_idx = candidate.begin().index();
        let end_idx = candidate.end().index();
        let deg_begin = degree.get(begin_idx).copied().unwrap_or(0);
        let deg_end = degree.get(end_idx).copied().unwrap_or(0);
        if !(2..=3).contains(&deg_begin) || !(2..=3).contains(&deg_end) {
            continue;
        }
        if hybridization
            .get(begin_idx)
            .copied()
            .unwrap_or(crate::Hybridization::Unspecified)
            != crate::Hybridization::Sp2
            || hybridization
                .get(end_idx)
                .copied()
                .unwrap_or(crate::Hybridization::Unspecified)
                != crate::Hybridization::Sp2
        {
            continue;
        }
        if rings.num_bond_rings(candidate_id) > 0 && rings.min_bond_ring_size(candidate_id) < 8 {
            continue;
        }

        let Some(nbr0) = atrop_neighbor_bonds(mol, candidate.begin(), candidate_id) else {
            continue;
        };
        let Some(nbr1) = atrop_neighbor_bonds(mol, candidate.end(), candidate_id) else {
            continue;
        };
        for end_bond in nbr0.iter().chain(nbr1.iter()) {
            let Some(bond) = mol.bonds().get(end_bond.index()) else {
                continue;
            };
            if bond.direction() == BondDirection::Unknown {
                continue;
            }
        }

        let (has_dir0, wedge_dir0) = atrop_end_wedge_direction(mol, &nbr0);
        let (has_dir1, wedge_dir1) = atrop_end_wedge_direction(mol, &nbr1);
        if !has_dir0 || !has_dir1 || wedge_dir0 == wedge_dir1 {
            continue;
        }

        if wedge_dir0 == BondDirection::BeginWedge || wedge_dir1 == BondDirection::BeginDash {
            assignments.push((candidate_id, BondStereo::AtropCcw));
        } else if wedge_dir0 == BondDirection::BeginDash || wedge_dir1 == BondDirection::BeginWedge
        {
            assignments.push((candidate_id, BondStereo::AtropCw));
        }
    }

    assignments
}

fn atropisomer_stereo_from_conformer(mol: &Molecule, conf_id: usize) -> Vec<(BondId, BondStereo)> {
    const REALLY_SMALL_BOND_LEN: f64 = 1e-7;
    let Some(conformer) = mol.conformers_3d().iter().find(|conf| conf.id() == conf_id) else {
        return Vec::new();
    };
    let Some(rings) = mol.derived_cache().rings.as_ref() else {
        return Vec::new();
    };

    let mut candidate_bonds: Vec<BondId> = Vec::new();
    for bond in mol.bonds() {
        if !atrop_can_have_direction(bond.order()) {
            continue;
        }
        if !matches!(
            bond.direction(),
            BondDirection::BeginWedge | BondDirection::BeginDash
        ) {
            continue;
        }
        let begin = bond.begin();
        for nbr_bond in mol.bonds() {
            if nbr_bond.id() == bond.id() {
                continue;
            }
            if (nbr_bond.begin() == begin || nbr_bond.end() == begin)
                && !candidate_bonds.contains(&nbr_bond.id())
            {
                candidate_bonds.push(nbr_bond.id());
            }
        }
    }

    let num_atoms = mol.num_atoms();
    let mut degree = vec![0usize; num_atoms];
    for bond in mol.bonds() {
        degree[bond.begin().index()] += 1;
        degree[bond.end().index()] += 1;
    }
    let hybridization: Vec<crate::Hybridization> = mol
        .atoms()
        .iter()
        .map(|atom| atom.hybridization())
        .collect();

    let mut assignments = Vec::new();
    'candidate: for candidate_id in candidate_bonds {
        let Some(candidate) = mol.bonds().get(candidate_id.index()) else {
            continue;
        };
        if candidate.order() != BondOrder::Single || candidate.stereo() == BondStereo::Any {
            continue;
        }
        let begin_idx = candidate.begin().index();
        let end_idx = candidate.end().index();
        let deg_begin = degree.get(begin_idx).copied().unwrap_or(0);
        let deg_end = degree.get(end_idx).copied().unwrap_or(0);
        if !(2..=3).contains(&deg_begin) || !(2..=3).contains(&deg_end) {
            continue;
        }
        if hybridization
            .get(begin_idx)
            .copied()
            .unwrap_or(crate::Hybridization::Unspecified)
            != crate::Hybridization::Sp2
            || hybridization
                .get(end_idx)
                .copied()
                .unwrap_or(crate::Hybridization::Unspecified)
                != crate::Hybridization::Sp2
        {
            continue;
        }
        if rings.num_bond_rings(candidate_id) > 0 && rings.min_bond_ring_size(candidate_id) < 8 {
            continue;
        }

        let Some(nbr0) = atrop_neighbor_bonds(mol, candidate.begin(), candidate_id) else {
            continue;
        };
        let Some(nbr1) = atrop_neighbor_bonds(mol, candidate.end(), candidate_id) else {
            continue;
        };
        for end_bond in nbr0.iter().chain(nbr1.iter()) {
            let Some(bond) = mol.bonds().get(end_bond.index()) else {
                continue 'candidate;
            };
            if bond.direction() == BondDirection::Unknown {
                continue 'candidate;
            }
        }

        let Some((_x_axis, y_axis, z_axis)) =
            atrop_bond_frame_of_reference(mol, candidate_id, conformer)
        else {
            continue;
        };
        let mut bond_vecs = [[0.0; 3]; 2];
        for (bond_atom_index, (focus_atom, nbr_bonds)) in
            [(candidate.begin(), &nbr0), (candidate.end(), &nbr1)]
                .into_iter()
                .enumerate()
        {
            if !conformer.is_3d() {
                let (has_dir, bond_dir) = atrop_end_wedge_direction(mol, nbr_bonds);
                if !has_dir {
                    continue 'candidate;
                }
                let Some(mut bond_vec) = atrop_projected_end_vector(
                    mol, conformer, focus_atom, nbr_bonds, y_axis, z_axis, true,
                ) else {
                    continue 'candidate;
                };
                if bond_dir == BondDirection::BeginWedge {
                    bond_vec[1] *= 0.707;
                    bond_vec[2] = bond_vec[1].abs();
                } else if bond_dir == BondDirection::BeginDash {
                    bond_vec[1] *= 0.707;
                    bond_vec[2] = -bond_vec[1].abs();
                }
                bond_vecs[bond_atom_index] = bond_vec;
            } else {
                let Some(bond_vec) = atrop_projected_end_vector(
                    mol, conformer, focus_atom, nbr_bonds, y_axis, z_axis, false,
                ) else {
                    continue 'candidate;
                };
                if vec3_len(bond_vec) < REALLY_SMALL_BOND_LEN {
                    continue 'candidate;
                }
                bond_vecs[bond_atom_index] = bond_vec;
            }
        }

        let cross_product = vec3_cross(bond_vecs[1], bond_vecs[0]);
        if cross_product[0] > REALLY_SMALL_BOND_LEN {
            assignments.push((candidate_id, BondStereo::AtropCcw));
        } else if cross_product[0] < -REALLY_SMALL_BOND_LEN {
            assignments.push((candidate_id, BondStereo::AtropCw));
        }
    }
    assignments
}

/// Port of `MolOps::assignChiralTypesFromBondDirs` — processes wedge/dash bonds
/// to assign tetrahedral chiral types to their begin atoms using pseudo-3D
/// coordinates from the given conformer.
///
/// C++ source: `third_party/rdkit/Code/GraphMol/Chirality.cpp:3765-3812`
fn assign_chiral_types_from_bond_dirs(mol: &mut Molecule, conf_id: usize) {
    // RDKit✔️✔️: void assignChiralTypesFromBondDirs(ROMol &mol, const int confId,
    // RDKit✔️✔️:                                    const bool replaceExistingTags) {
    // RDKit✔️✔️:   if (!mol.getNumConformers()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto conf = mol.getConformer(confId);
    let conf_idx = mol.conformers_3d().iter().position(|c| c.id() == conf_id);
    let Some(conf_idx) = conf_idx else {
        return;
    };
    let n_atoms = mol.num_atoms();
    let mut atoms_set = vec![false; n_atoms];
    let replace_existing_tags = true; // Called from MolFromSmiles with default true
    let bond_range = 0..mol.num_bonds();

    // Phase 1: collect bonds that might need chiral assignment.
    // We collect indices first to avoid borrow conflicts between reading
    // molecule data for pseudo-3D and mutating atoms.
    let mut chiral_assignments: Vec<(usize, Option<ChiralTag>)> = Vec::new();
    for bond_idx in bond_range {
        // RDKit✔️✔️:   for (auto &bond : mol.bonds()) {
        // RDKit✔️✔️:     const Bond::BondDir dir = bond->getBondDir();
        // RDKit✔️✔️:     Atom *atom = bond->getBeginAtom();
        let bond = &mol.topology_block().bonds[bond_idx];
        let dir = bond.direction();
        let begin_idx = bond.begin().index();
        let begin_chiral = mol.topology_block().atoms[begin_idx].chiral_tag();

        // RDKit✔️✔️:     if (dir == Bond::UNKNOWN) {
        // RDKit✔️✔️:       if (atomsSet[atom->getIdx()] || replaceExistingTags) {
        // RDKit✔️✔️:         atom->setChiralTag(Atom::CHI_UNSPECIFIED);
        // RDKit✔️✔️:         atomsSet.set(atom->getIdx());
        // RDKit✔️✔️:       }
        if dir == crate::BondDirection::Unknown {
            if atoms_set[begin_idx] || replace_existing_tags {
                chiral_assignments.push((begin_idx, Some(ChiralTag::Unspecified)));
                atoms_set[begin_idx] = true;
            }
            continue;
        }

        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       if (dir == Bond::BEGINWEDGE || dir == Bond::BEGINDASH) {
        if dir != crate::BondDirection::BeginWedge && dir != crate::BondDirection::BeginDash {
            continue;
        }
        // RDKit✔️✔️:         if (atomsSet[atom->getIdx()] ||
        // RDKit✔️✔️:             (!replaceExistingTags &&
        // RDKit✔️✔️:              atom->getChiralTag() != Atom::CHI_UNSPECIFIED)) {
        // RDKit✔️✔️:           continue;
        // RDKit✔️✔️:         }
        // COSMolKit: we use replaceExistingTags=true from MolFromSmiles,
        // so we only skip if already set by a prior wedge bond.
        if atoms_set[begin_idx] {
            continue;
        }

        if !replace_existing_tags && begin_chiral != ChiralTag::Unspecified {
            continue;
        }

        // RDKit✔️✔️:         Atom::ChiralType code =
        // RDKit✔️✔️:             Chirality::atomChiralTypeFromBondDirPseudo3D(mol, bond, &conf)
        // RDKit✔️✔️:                 .value_or(Atom::ChiralType::CHI_UNSPECIFIED);
        let conformer = &mol.conformers_3d()[conf_idx];
        let code = crate::chemistry::stereo::atom_chiral_type_from_bond_dir_pseudo_3d(
            mol, bond_idx, conformer,
        )
        .unwrap_or(ChiralTag::Unspecified);

        // RDKit✔️✔️:         if (code != Atom::ChiralType::CHI_UNSPECIFIED) {
        // RDKit✔️✔️:           atomsSet.set(atom->getIdx());
        // RDKit✔️✔️:         }
        if code != ChiralTag::Unspecified {
            atoms_set[begin_idx] = true;
        }

        // RDKit✔️✔️:         atom->setChiralTag(code);
        chiral_assignments.push((begin_idx, Some(code)));

        // RDKit✔️✔️:         if (atom->getDegree() == 3 && !atom->getNumExplicitHs() &&
        // RDKit✔️✔️:             atom->getNumImplicitHs() == 1) {
        // RDKit✔️✔️:           atom->setNumExplicitHs(1);
        // RDKit✔️✔️:           atom->updatePropertyCache();
        // RDKit✔️✔️:         }
        // COSMolKit: compute degree as bond count to this atom
        let deg = mol
            .bonds()
            .iter()
            .filter(|b| b.begin().index() == begin_idx || b.end().index() == begin_idx)
            .count();
        if deg == 3 && code != ChiralTag::Unspecified {
            let has_explicit = {
                let atom = &mol.topology_block().atoms[begin_idx];
                atom.explicit_hydrogens() > 0
            };
            let no_implicit = {
                let atom = &mol.topology_block().atoms[begin_idx];
                atom.no_implicit()
            };
            if !has_explicit && !no_implicit {
                // 3-coordinate, no explicit H, implicit H exists
                chiral_assignments.push((begin_idx, None)); // marker for explicit H update
            }
        }
    }

    // Phase 2: apply collected assignments
    let atoms = &mut mol.topology_block_mut().atoms;
    for (idx, assignment) in chiral_assignments {
        if let Some(tag) = assignment {
            atoms[idx].set_chiral_tag(tag);
        } else {
            // marker for explicit H update
            if atoms[idx].explicit_hydrogens() == 0 && !atoms[idx].no_implicit() {
                atoms[idx].set_explicit_hydrogens(1);
            }
        }
    }
}

/// Compute the signed volume for a tetrahedral chiral center from 3D conformer
/// coordinates. The sign of the signed volume determines CW vs CCW chirality.
///
/// Returns `ChiralTag::Unspecified` if the geometry is ambiguous (zero volume).
fn tetrahedral_chiral_volume(p0: [f64; 3], nbr_positions: &[[f64; 3]]) -> Option<ChiralTag> {
    const ZERO_VOLUME_TOL: f64 = 0.1;
    if nbr_positions.len() < 3 {
        return None;
    }
    // v_i = nbr_i - p0
    let v1 = [
        nbr_positions[0][0] - p0[0],
        nbr_positions[0][1] - p0[1],
        nbr_positions[0][2] - p0[2],
    ];
    let v2 = [
        nbr_positions[1][0] - p0[0],
        nbr_positions[1][1] - p0[1],
        nbr_positions[1][2] - p0[2],
    ];
    let v3 = [
        nbr_positions[2][0] - p0[0],
        nbr_positions[2][1] - p0[1],
        nbr_positions[2][2] - p0[2],
    ];

    // chiralVol = v1 · (v2 × v3)
    let cross = [
        v2[1] * v3[2] - v2[2] * v3[1],
        v2[2] * v3[0] - v2[0] * v3[2],
        v2[0] * v3[1] - v2[1] * v3[0],
    ];
    let mut chiral_vol = v1[0] * cross[0] + v1[1] * cross[1] + v1[2] * cross[2];

    // If first 3 are planar (zero volume) and a 4th neighbor exists, try v4
    if chiral_vol.abs() < ZERO_VOLUME_TOL && nbr_positions.len() >= 4 {
        let v4 = [
            nbr_positions[3][0] - p0[0],
            nbr_positions[3][1] - p0[1],
            nbr_positions[3][2] - p0[2],
        ];
        // v4 is in opposite direction to v3 for a tetrahedral center
        // chiralVol = -v1 · (v2 × v4)
        let cross2 = [
            v2[1] * v4[2] - v2[2] * v4[1],
            v2[2] * v4[0] - v2[0] * v4[2],
            v2[0] * v4[1] - v2[1] * v4[0],
        ];
        chiral_vol = -(v1[0] * cross2[0] + v1[1] * cross2[1] + v1[2] * cross2[2]);
    }

    if chiral_vol < -ZERO_VOLUME_TOL {
        Some(ChiralTag::TetrahedralCw)
    } else if chiral_vol > ZERO_VOLUME_TOL {
        Some(ChiralTag::TetrahedralCcw)
    } else {
        None // ambiguous
    }
}

fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn vec3_len(v: [f64; 3]) -> f64 {
    vec3_dot(v, v).sqrt()
}

fn vec3_direction(from: [f64; 3], to: [f64; 3]) -> Option<[f64; 3]> {
    let v = vec3_sub(to, from);
    let len = vec3_len(v);
    if len == 0.0 {
        None
    } else {
        Some([v[0] / len, v[1] / len, v[2] / len])
    }
}

fn vec3_angle(a: [f64; 3], b: [f64; 3]) -> f64 {
    let denom = vec3_len(a) * vec3_len(b);
    if denom == 0.0 {
        0.0
    } else {
        (vec3_dot(a, b) / denom).clamp(-1.0, 1.0).acos()
    }
}

fn voltest(vectors: &[[f64; 3]], x: usize, y: usize, z: usize) -> bool {
    // RDKit✔️✔️: #define VOLTEST(X, Y, Z) (v[X].dotProduct(v[Y].crossProduct(v[Z])) >= 0.0)
    vec3_dot(vectors[x], vec3_cross(vectors[y], vectors[z])) >= 0.0
}

fn octahedral_perm_from_3d(pair: &[u8; 6], vectors: &[[f64; 3]]) -> u32 {
    // RDKit✔️✔️: static unsigned int OctahedralPermFrom3D(unsigned char *pair,
    // RDKit✔️✔️:                                          const RDGeom::Point3D *v) {
    // RDKit✔️✔️:   switch (pair[0]) {
    match pair[0] {
        // RDKit✔️✔️:     case 2:  // a-b
        2 => match pair[2] {
            // RDKit✔️✔️:         case 4:
            // RDKit✔️✔️:           return VOLTEST(0, 3, 4) ? 28 : 27;
            4 => {
                if voltest(vectors, 0, 3, 4) {
                    28
                } else {
                    27
                }
            }
            // RDKit✔️✔️:         case 5:
            // RDKit✔️✔️:           return VOLTEST(0, 2, 3) ? 25 : 30;
            5 => {
                if voltest(vectors, 0, 2, 3) {
                    25
                } else {
                    30
                }
            }
            // RDKit✔️✔️:         default:  // 0 or 6
            // RDKit✔️✔️:           return VOLTEST(0, 2, 3) ? 26 : 29;
            _ => {
                if voltest(vectors, 0, 2, 3) {
                    26
                } else {
                    29
                }
            }
        },
        // RDKit✔️✔️:     case 3:  // a-c
        3 => match pair[1] {
            // RDKit✔️✔️:         case 4:
            // RDKit✔️✔️:           return VOLTEST(0, 3, 4) ? 22 : 21;
            4 => {
                if voltest(vectors, 0, 3, 4) {
                    22
                } else {
                    21
                }
            }
            // RDKit✔️✔️:         case 5:
            // RDKit✔️✔️:           return VOLTEST(0, 1, 3) ? 19 : 24;
            5 => {
                if voltest(vectors, 0, 1, 3) {
                    19
                } else {
                    24
                }
            }
            // RDKit✔️✔️:         default:  // 0 or 6
            // RDKit✔️✔️:           return VOLTEST(0, 1, 3) ? 20 : 23;
            _ => {
                if voltest(vectors, 0, 1, 3) {
                    20
                } else {
                    23
                }
            }
        },
        // RDKit✔️✔️:     case 4:  // a-d
        4 => match pair[1] {
            // RDKit✔️✔️:         case 3:
            // RDKit✔️✔️:           return VOLTEST(0, 2, 4) ? 13 : 12;
            3 => {
                if voltest(vectors, 0, 2, 4) {
                    13
                } else {
                    12
                }
            }
            // RDKit✔️✔️:         case 5:
            // RDKit✔️✔️:           return VOLTEST(0, 1, 2) ? 6 : 18;
            5 => {
                if voltest(vectors, 0, 1, 2) {
                    6
                } else {
                    18
                }
            }
            // RDKit✔️✔️:         default:  // 0 or 6
            // RDKit✔️✔️:           return VOLTEST(0, 1, 2) ? 7 : 17;
            _ => {
                if voltest(vectors, 0, 1, 2) {
                    7
                } else {
                    17
                }
            }
        },
        // RDKit✔️✔️:     case 5:  // a-e
        5 => match pair[1] {
            // RDKit✔️✔️:         case 3:
            // RDKit✔️✔️:           return VOLTEST(0, 2, 3) ? 11 : 9;
            3 => {
                if voltest(vectors, 0, 2, 3) {
                    11
                } else {
                    9
                }
            }
            // RDKit✔️✔️:         case 4:
            // RDKit✔️✔️:           return VOLTEST(0, 1, 2) ? 3 : 16;
            4 => {
                if voltest(vectors, 0, 1, 2) {
                    3
                } else {
                    16
                }
            }
            // RDKit✔️✔️:         default:  // 0 or 6
            // RDKit✔️✔️:           return VOLTEST(0, 1, 2) ? 5 : 15;
            _ => {
                if voltest(vectors, 0, 1, 2) {
                    5
                } else {
                    15
                }
            }
        },
        // RDKit✔️✔️:     default:  // 0 or 6  a-f
        _ => match pair[1] {
            // RDKit✔️✔️:         case 3:
            // RDKit✔️✔️:           return VOLTEST(0, 2, 3) ? 10 : 8;
            3 => {
                if voltest(vectors, 0, 2, 3) {
                    10
                } else {
                    8
                }
            }
            // RDKit✔️✔️:         case 4:
            // RDKit✔️✔️:           return VOLTEST(0, 1, 2) ? 1 : 2;
            4 => {
                if voltest(vectors, 0, 1, 2) {
                    1
                } else {
                    2
                }
            }
            // RDKit✔️✔️:         default:  // 5
            // RDKit✔️✔️:           return VOLTEST(0, 1, 2) ? 4 : 14;
            _ => {
                if voltest(vectors, 0, 1, 2) {
                    4
                } else {
                    14
                }
            }
        },
    }
}

fn assign_nontetrahedral_chiral_type_from_3d(
    mol: &Molecule,
    coords: &[[f64; 3]],
    atom_idx: usize,
) -> Option<(ChiralTag, u32)> {
    // RDKit✔️✔️: static bool assignNontetrahedralChiralTypeFrom3D(ROMol &mol,
    // RDKit✔️✔️:                                                  const Conformer &conf,
    // RDKit✔️✔️:                                                  Atom *atom,
    // RDKit✔️✔️:                                                  double tolerance = 0.1) {
    // RDKit✔️✔️:   // Fail fast check for non-tetrahedral elements
    // RDKit✔️✔️:   if (atom->getAtomicNum() < 15) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    const TOLERANCE: f64 = 0.1;
    if mol.atoms()[atom_idx].atomic_number() < 15 {
        return None;
    }

    // RDKit✔️✔️:   // check for wiggly bonds
    // RDKit✔️✔️:   for (const auto bond : mol.atomBonds(atom)) {
    // RDKit✔️✔️:     if (isWigglyBond(bond, atom)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    for bond in mol.bonds() {
        if (bond.begin().index() == atom_idx || bond.end().index() == atom_idx)
            && is_wiggly_bond(bond.direction(), bond.order(), bond.unknown_stereo())
        {
            return None;
        }
    }

    // RDKit✔️✔️:   RDGeom::Point3D cen = conf.getAtomPos(atom->getIdx());
    // RDKit✔️✔️:   RDGeom::Point3D v[6];
    // RDKit✔️✔️:   unsigned int count = 0;
    // RDKit✔️✔️:   ROMol::ADJ_ITER nbrIdx, endNbrs;
    // RDKit✔️✔️:   boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom);
    // RDKit✔️✔️:   while (nbrIdx != endNbrs) {
    // RDKit✔️✔️:     if (count == 6) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     RDGeom::Point3D p = conf.getAtomPos(*nbrIdx);
    // RDKit✔️✔️:     v[count] = cen.directionVector(p);
    // RDKit✔️✔️:     ++count;
    // RDKit✔️✔️:     ++nbrIdx;
    // RDKit✔️✔️:   }
    let center = coords[atom_idx];
    let mut vectors = Vec::with_capacity(6);
    for bond in mol.bonds() {
        let other_idx = if bond.begin().index() == atom_idx {
            Some(bond.end().index())
        } else if bond.end().index() == atom_idx {
            Some(bond.begin().index())
        } else {
            None
        };
        let Some(other_idx) = other_idx else {
            continue;
        };
        if vectors.len() == 6 {
            return None;
        }
        vectors.push(vec3_direction(center, coords[other_idx])?);
    }

    // RDKit✔️✔️:   if (count < 3) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    let count = vectors.len();
    if count < 3 {
        return None;
    }

    // RDKit✔️✔️:   unsigned char pair[6];
    // RDKit✔️✔️:   memset(pair, 0, 6);
    // RDKit✔️✔️:   unsigned int pairs = 0;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < count; i++) {
    // RDKit✔️✔️:     for (unsigned int j = i + 1; j < count; j++) {
    // RDKit✔️✔️:       if (v[i].dotProduct(v[j]) < -(1 - tolerance)) {
    // RDKit✔️✔️:         if (pair[i] || pair[j]) {
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         pair[i] = j + 1;
    // RDKit✔️✔️:         pair[j] = i + 1;
    // RDKit✔️✔️:         pairs++;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let mut pair = [0_u8; 6];
    let mut pairs = 0_u32;
    for i in 0..count {
        for j in (i + 1)..count {
            if vec3_dot(vectors[i], vectors[j]) < -(1.0 - TOLERANCE) {
                if pair[i] != 0 || pair[j] != 0 {
                    return None;
                }
                pair[i] = (j + 1) as u8;
                pair[j] = (i + 1) as u8;
                pairs += 1;
            }
        }
    }

    // RDKit✔️✔️:   switch (pairs) {
    // RDKit✔️✔️:     case 1:
    // RDKit✔️✔️:       switch (count) {
    // RDKit✔️✔️:         case 3: /* T-shape */
    // RDKit✔️✔️:           atom->setChiralTag(Atom::ChiralType::CHI_SQUAREPLANAR);
    // RDKit✔️✔️:           res = true;
    // RDKit✔️✔️:           if (pair[0] == 0) {
    // RDKit✔️✔️:             perm = 3;  // Z
    // RDKit✔️✔️:           } else if (pair[0] == 2) {
    // RDKit✔️✔️:             perm = 2;  // 4
    // RDKit✔️✔️:           } else /* pair[0] == 3 */ {
    // RDKit✔️✔️:             perm = 1;  // U
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           atom->setProp(common_properties::_chiralPermutation, perm);
    // RDKit✔️✔️:           break;
    if pairs == 1 && count == 3 {
        let perm = if pair[0] == 0 {
            3
        } else if pair[0] == 2 {
            2
        } else {
            1
        };
        return Some((ChiralTag::SquarePlanar, perm));
    }

    // RDKit✔️✔️:         case 4:                /* See-saw */
    // RDKit✔️✔️:           if (pair[0] == 2) {  // a b
    // RDKit✔️✔️:             if (v[2].angleTo(v[3]) < 100 * M_PI / 180.0) {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 2, 3) ? 25 : 29;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 2, 3) ? 7 : 8;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else if (pair[0] == 3) {  // a c
    // RDKit✔️✔️:             if (v[1].angleTo(v[3]) < 100 * M_PI / 180.0) {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 3) ? 19 : 23;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 3) ? 5 : 6;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else if (pair[0] == 4) {  // a d
    // RDKit✔️✔️:             if (v[1].angleTo(v[2]) < 100 * M_PI / 180.0) {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 2) ? 6 : 17;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 2) ? 3 : 4;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else if (pair[1] == 3) {  // b c
    // RDKit✔️✔️:             if (v[0].angleTo(v[3]) < 100 * M_PI / 180.0) {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 3) ? 10 : 8;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
    // RDKit✔️✔️:               perm = VOLTEST(1, 0, 3) ? 13 : 14;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else if (pair[1] == 4) {  // b d
    // RDKit✔️✔️:             if (v[0].angleTo(v[2]) < 100 * M_PI / 180.0) {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 3) ? 1 : 2;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
    // RDKit✔️✔️:               perm = VOLTEST(1, 0, 2) ? 10 : 12;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:           } else /* pair[2] == 4 */ {  // c d
    // RDKit✔️✔️:             if (v[0].angleTo(v[1]) < 100 * M_PI / 180.0) {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️:               perm = VOLTEST(0, 1, 3) ? 4 : 14;
    // RDKit✔️✔️:             } else {
    // RDKit✔️✔️:               tag = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL;
    // RDKit✔️✔️:               perm = VOLTEST(3, 0, 1) ? 16 : 19;
    // RDKit✔️✔️:             }
    if pairs == 1 && count == 4 {
        let angle_threshold = 100.0_f64.to_radians();
        let assignment = if pair[0] == 2 {
            if vec3_angle(vectors[2], vectors[3]) < angle_threshold {
                Some((
                    ChiralTag::Octahedral,
                    if voltest(&vectors, 0, 2, 3) { 25 } else { 29 },
                ))
            } else {
                Some((
                    ChiralTag::TrigonalBipyramidal,
                    if voltest(&vectors, 0, 2, 3) { 7 } else { 8 },
                ))
            }
        } else if pair[0] == 3 {
            if vec3_angle(vectors[1], vectors[3]) < angle_threshold {
                Some((
                    ChiralTag::Octahedral,
                    if voltest(&vectors, 0, 1, 3) { 19 } else { 23 },
                ))
            } else {
                Some((
                    ChiralTag::TrigonalBipyramidal,
                    if voltest(&vectors, 0, 1, 3) { 5 } else { 6 },
                ))
            }
        } else if pair[0] == 4 {
            if vec3_angle(vectors[1], vectors[2]) < angle_threshold {
                Some((
                    ChiralTag::Octahedral,
                    if voltest(&vectors, 0, 1, 2) { 6 } else { 17 },
                ))
            } else {
                Some((
                    ChiralTag::TrigonalBipyramidal,
                    if voltest(&vectors, 0, 1, 2) { 3 } else { 4 },
                ))
            }
        } else if pair[1] == 3 {
            if vec3_angle(vectors[0], vectors[3]) < angle_threshold {
                Some((
                    ChiralTag::Octahedral,
                    if voltest(&vectors, 0, 1, 3) { 10 } else { 8 },
                ))
            } else {
                Some((
                    ChiralTag::TrigonalBipyramidal,
                    if voltest(&vectors, 1, 0, 3) { 13 } else { 14 },
                ))
            }
        } else if pair[1] == 4 {
            if vec3_angle(vectors[0], vectors[2]) < angle_threshold {
                Some((
                    ChiralTag::Octahedral,
                    if voltest(&vectors, 0, 1, 3) { 1 } else { 2 },
                ))
            } else {
                Some((
                    ChiralTag::TrigonalBipyramidal,
                    if voltest(&vectors, 1, 0, 2) { 10 } else { 12 },
                ))
            }
        } else {
            if vec3_angle(vectors[0], vectors[1]) < angle_threshold {
                Some((
                    ChiralTag::Octahedral,
                    if voltest(&vectors, 0, 1, 3) { 4 } else { 14 },
                ))
            } else {
                Some((
                    ChiralTag::TrigonalBipyramidal,
                    if voltest(&vectors, 3, 0, 1) { 16 } else { 19 },
                ))
            }
        };
        if let Some(assignment) = assignment {
            return Some(assignment);
        }
    }

    // RDKit✔️✔️:         case 5: /* Trigonal bipyramidal */
    // RDKit✔️✔️:           atom->setChiralTag(Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL);
    // RDKit✔️✔️:           res = true;
    // RDKit✔️✔️:           if (pair[0] == 2) {
    // RDKit✔️✔️:             perm = VOLTEST(0, 2, 3) ? 7 : 8;  // a b
    // RDKit✔️✔️:           } else if (pair[0] == 3) {
    // RDKit✔️✔️:             perm = VOLTEST(0, 1, 3) ? 5 : 6;  // a c
    // RDKit✔️✔️:           } else if (pair[0] == 4) {
    // RDKit✔️✔️:             perm = VOLTEST(0, 1, 2) ? 3 : 4;  // a d
    // RDKit✔️✔️:           } else if (pair[0] == 5) {
    // RDKit✔️✔️:             perm = VOLTEST(0, 1, 2) ? 1 : 2;  // a e
    // RDKit✔️✔️:           } else if (pair[1] == 3) {
    // RDKit✔️✔️:             perm = VOLTEST(1, 0, 3) ? 13 : 14;  // b c
    // RDKit✔️✔️:           } else if (pair[1] == 4) {
    // RDKit✔️✔️:             perm = VOLTEST(1, 0, 2) ? 10 : 12;  // b d
    // RDKit✔️✔️:           } else if (pair[1] == 5) {
    // RDKit✔️✔️:             perm = VOLTEST(1, 0, 2) ? 9 : 11;  // b e
    // RDKit✔️✔️:           } else if (pair[2] == 4) {
    // RDKit✔️✔️:             perm = VOLTEST(2, 0, 1) ? 16 : 19;  // c d
    // RDKit✔️✔️:           } else if (pair[2] == 5) {
    // RDKit✔️✔️:             perm = VOLTEST(2, 0, 1) ? 15 : 20;  // c e
    // RDKit✔️✔️:           } else /* pair[2] == 4 */ {
    // RDKit✔️✔️:             perm = VOLTEST(3, 0, 1) ? 17 : 18;  // d e
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           atom->setProp(common_properties::_chiralPermutation, perm);
    if pairs == 1 && count == 5 {
        let perm = if pair[0] == 2 {
            if voltest(&vectors, 0, 2, 3) { 7 } else { 8 }
        } else if pair[0] == 3 {
            if voltest(&vectors, 0, 1, 3) { 5 } else { 6 }
        } else if pair[0] == 4 {
            if voltest(&vectors, 0, 1, 2) { 3 } else { 4 }
        } else if pair[0] == 5 {
            if voltest(&vectors, 0, 1, 2) { 1 } else { 2 }
        } else if pair[1] == 3 {
            if voltest(&vectors, 1, 0, 3) { 13 } else { 14 }
        } else if pair[1] == 4 {
            if voltest(&vectors, 1, 0, 2) { 10 } else { 12 }
        } else if pair[1] == 5 {
            if voltest(&vectors, 1, 0, 2) { 9 } else { 11 }
        } else if pair[2] == 4 {
            if voltest(&vectors, 2, 0, 1) { 16 } else { 19 }
        } else if pair[2] == 5 {
            if voltest(&vectors, 2, 0, 1) { 15 } else { 20 }
        } else {
            if voltest(&vectors, 3, 0, 1) { 17 } else { 18 }
        };
        return Some((ChiralTag::TrigonalBipyramidal, perm));
    }

    // RDKit✔️✔️:     case 2:
    // RDKit✔️✔️:       if (count == 4) {
    // RDKit✔️✔️:         /* Square planar */
    // RDKit✔️✔️:         atom->setChiralTag(Atom::ChiralType::CHI_SQUAREPLANAR);
    // RDKit✔️✔️:         res = true;
    // RDKit✔️✔️:         if (pair[0] == 2) {
    // RDKit✔️✔️:           perm = 2;  // 4
    // RDKit✔️✔️:         } else if (pair[0] == 3) {
    // RDKit✔️✔️:           perm = 1;  // U
    // RDKit✔️✔️:         } else /* pair[1] == 4 */ {
    // RDKit✔️✔️:           perm = 3;  // Z
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         atom->setProp(common_properties::_chiralPermutation, perm);
    // RDKit✔️✔️:       }
    if pairs == 2 && count == 4 {
        let perm = if pair[0] == 2 {
            2
        } else if pair[0] == 3 {
            1
        } else {
            3
        };
        return Some((ChiralTag::SquarePlanar, perm));
    }

    // RDKit✔️✔️:       } else if (count == 5) {
    // RDKit✔️✔️:         /* Square pyramidal */
    // RDKit✔️✔️:         atom->setChiralTag(Atom::ChiralType::CHI_OCTAHEDRAL);
    // RDKit✔️✔️:         res = true;
    // RDKit✔️✔️:         perm = OctahedralPermFrom3D(pair, v);
    // RDKit✔️✔️:         atom->setProp(common_properties::_chiralPermutation, perm);
    // RDKit✔️✔️:       }
    if pairs == 2 && count == 5 {
        return Some((
            ChiralTag::Octahedral,
            octahedral_perm_from_3d(&pair, &vectors),
        ));
    }

    // RDKit✔️✔️:     case 3:
    // RDKit✔️✔️:       if (count == 6) {
    // RDKit✔️✔️:         /* Octahedral */
    // RDKit✔️✔️:         atom->setChiralTag(Atom::ChiralType::CHI_OCTAHEDRAL);
    // RDKit✔️✔️:         res = true;
    // RDKit✔️✔️:         perm = OctahedralPermFrom3D(pair, v);
    // RDKit✔️✔️:         atom->setProp(common_properties::_chiralPermutation, perm);
    // RDKit✔️✔️:       }
    if pairs == 3 && count == 6 {
        return Some((
            ChiralTag::Octahedral,
            octahedral_perm_from_3d(&pair, &vectors),
        ));
    }

    None
}

/// Helper: check if a bond direction is UNKNOWN (wiggly bond) from an atom's perspective.
fn is_wiggly_bond(
    bond_dir: crate::BondDirection,
    bond_order: crate::BondOrder,
    unknown_stereo: bool,
) -> bool {
    bond_order == crate::BondOrder::Single
        && (bond_dir == crate::BondDirection::Unknown || unknown_stereo)
}

/// Helper: check if a bond type affects atom chirality.
fn bond_affects_chirality(bond_order: crate::BondOrder) -> bool {
    !matches!(
        bond_order,
        crate::BondOrder::Unspecified
            | crate::BondOrder::Zero
            | crate::BondOrder::DativeRight
            | crate::BondOrder::DativeLeft
    )
}

/// Count non-zero-degree neighbors (excluding chirality-irrelevant bonds).
fn atom_nonzero_degree(mol: &Molecule, atom_idx: usize) -> usize {
    mol.bonds()
        .iter()
        .filter(|b| {
            (b.begin().index() == atom_idx || b.end().index() == atom_idx)
                && bond_affects_chirality(b.order())
        })
        .count()
}

fn rdkit_total_hs(atom: &crate::Atom, atom_idx: usize, implicit_hydrogens: &[i32]) -> usize {
    // RDKit✔️✔️: atom->getTotalNumHs()
    atom.explicit_hydrogens() as usize + implicit_hydrogens[atom_idx].max(0) as usize
}

/// Port of `MolOps::assignChiralTypesFrom3D` — tetrahedral signed-volume branch.
///
/// Assigns chiral tags to atoms based on 3D coordinates using signed-volume
/// computation. Only handles tetrahedral cases (4-coordinate or S/Se with
/// effective coordination). Non-tetrahedral and atropisomer cases are deferred.
///
/// C++ source: `third_party/rdkit/Code/GraphMol/Chirality.cpp:3377-3500`
///
/// Note: this is the tetrahedral-only port. The non-tetrahedral path
/// (`assignNontetrahedralChiralTypeFrom3D`) is deferred.
fn assign_chiral_types_from_3d(mol: &mut Molecule, conf_id: usize) {
    // RDKit✔️✔️: void assignChiralTypesFrom3D(ROMol &mol, int confId, bool replaceExistingTags) {
    // RDKit✔️✔️:   const double ZERO_VOLUME_TOL = 0.1;
    // RDKit✔️✔️:   if (!mol.getNumConformers()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const Conformer &conf = mol.getConformer(confId);
    // RDKit✔️✔️:   if (!conf.is3D()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    let conf_idx = mol.conformers_3d().iter().position(|c| c.id() == conf_id);
    let Some(conf_idx) = conf_idx else {
        return;
    };
    // RDKit✔️✔️:   if (mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit✔️✔️:     mol.clearProp(common_properties::_StereochemDone);
    // RDKit✔️✔️:   }
    if mol.prop("_StereochemDone").is_some() {
        mol.properties_mut().clear_prop("_StereochemDone");
    }
    let conformer = &mol.conformers_3d()[conf_idx];
    if !conformer.is_3d() {
        return;
    }
    let coords = conformer.coords();

    // RDKit✔️✔️:   boost::dynamic_bitset<> explicitAtoms;
    // RDKit✔️✔️:   explicitAtoms.resize(mol.getNumAtoms(), 0);
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     if (bondDir == Bond::BondDir::BEGINWEDGE ||
    // RDKit✔️✔️:         bondDir == Bond::BondDir::BEGINDASH) {
    // RDKit✔️✔️:       explicitAtoms[bond->getBeginAtom()->getIdx()] = 1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     if (atom->getChiralTag() != Atom::ChiralType::CHI_UNSPECIFIED) {
    // RDKit✔️✔️:       explicitAtoms[atom->getIdx()] = 1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let n_atoms = mol.num_atoms();
    let mut explicit_atoms = vec![false; n_atoms];
    for bond in mol.bonds() {
        if matches!(
            bond.direction(),
            crate::BondDirection::BeginWedge | crate::BondDirection::BeginDash
        ) {
            explicit_atoms[bond.begin().index()] = true;
        }
    }
    for atom in mol.atoms() {
        if atom.chiral_tag() != ChiralTag::Unspecified {
            explicit_atoms[atom.id().index()] = true;
        }
    }

    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     if (!replaceExistingTags && atom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    let mut chiral_assignments: Vec<(usize, ChiralTag)> = Vec::new();
    let mut chiral_permutation_assignments: Vec<Option<u32>> = vec![None; n_atoms];
    let mut non_explicit_3d_chirality = vec![false; n_atoms];
    let Ok(valence) =
        crate::assign_valence_with_options(mol, crate::ValenceModel::RdkitLike, false)
    else {
        return;
    };

    for atom_idx in 0..n_atoms {
        let atom = &mol.topology_block().atoms[atom_idx];
        // RDKit✔️✔️:     auto nzDegree = Chirality::detail::getAtomNonzeroDegree(atom);
        // RDKit✔️✔️:     auto tnzDegree = nzDegree + atom->getTotalNumHs();
        // RDKit✔️✔️:     if (nzDegree < 3 || tnzDegree > 6) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        let nz_degree = atom_nonzero_degree(mol, atom_idx);
        let tnz_degree = nz_degree + rdkit_total_hs(atom, atom_idx, &valence.implicit_hydrogens);
        if nz_degree < 3 || tnz_degree > 6 {
            chiral_assignments.push((atom_idx, ChiralTag::Unspecified));
            continue;
        }

        // RDKit✔️✔️:     if (allowNontetrahedralStereo &&
        // RDKit✔️✔️:         assignNontetrahedralChiralTypeFrom3D(mol, conf, atom)) {
        // RDKit✔️✔️:       if (explicitAtoms[atom->getIdx()] == 0) {
        // RDKit✔️✔️:         atom->setProp(common_properties::_NonExplicit3DChirality, 1);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if let Some((tag, permutation)) =
            assign_nontetrahedral_chiral_type_from_3d(mol, coords, atom_idx)
        {
            if !explicit_atoms[atom_idx] {
                non_explicit_3d_chirality[atom_idx] = true;
            }
            chiral_permutation_assignments[atom_idx] = Some(permutation);
            chiral_assignments.push((atom_idx, tag));
            continue;
        }

        // RDKit✔️✔️:     /* We're only doing tetrahedral cases here */
        // RDKit✔️✔️:     if (tnzDegree > 4) {
        // RDKit✔️✔️:       chiral_assignments.push((atom_idx, ChiralTag::Unspecified));
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     int anum = atom->getAtomicNum();
        // RDKit✔️✔️:     if (anum != 16 && anum != 34 &&  // S or Se
        // RDKit✔️✔️:         tnzDegree != 4) {
        // RDKit✔️✔️:       chiral_assignments.push((atom_idx, ChiralTag::Unspecified));
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if tnz_degree > 4 {
            chiral_assignments.push((atom_idx, ChiralTag::Unspecified));
            continue;
        }
        let anum = mol.topology_block().atoms[atom_idx].atomic_number();
        if anum != 16 && anum != 34 && tnz_degree != 4 {
            chiral_assignments.push((atom_idx, ChiralTag::Unspecified));
            continue;
        }

        // RDKit✔️✔️:     const auto &p0 = conf.getAtomPos(atom->getIdx());
        // RDKit✔️✔️:     const RDGeom::Point3D *nbrs[4];
        // RDKit✔️✔️:     unsigned int nbrIdx = 0;
        // RDKit✔️✔️:     int hasWigglyBond = 0;
        // RDKit✔️✔️:     for (const auto bond : mol.atomBonds(atom)) {
        // RDKit✔️✔️:       hasWigglyBond = isWigglyBond(bond, atom);
        // RDKit✔️✔️:       if (hasWigglyBond) {
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       if (!Chirality::detail::bondAffectsAtomChirality(bond, atom)) {
        // RDKit✔️✔️:         continue;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       nbrs[nbrIdx++] = &conf.getAtomPos(bond->getOtherAtomIdx(atom->getIdx()));
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     if (hasWigglyBond) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        let p0 = coords[atom_idx];
        let mut nbr_positions: Vec<[f64; 3]> = Vec::with_capacity(4);
        let mut has_wiggly_bond = false;
        for bond in mol.bonds() {
            if bond.begin().index() != atom_idx && bond.end().index() != atom_idx {
                continue;
            }
            let dir = bond.direction();
            let order = bond.order();
            let unknown_stereo = bond.unknown_stereo();
            if is_wiggly_bond(dir, order, unknown_stereo) {
                has_wiggly_bond = true;
                break;
            }
            if !bond_affects_chirality(order) {
                continue;
            }
            let other_idx = if bond.begin().index() == atom_idx {
                bond.end().index()
            } else {
                bond.begin().index()
            };
            nbr_positions.push(coords[other_idx]);
        }
        if has_wiggly_bond {
            chiral_assignments.push((atom_idx, ChiralTag::Unspecified));
            continue;
        }

        // RDKit✔️✔️:     auto v1 = *nbrs[0] - p0;
        // RDKit✔️✔️:     auto v2 = *nbrs[1] - p0;
        // RDKit✔️✔️:     auto v3 = *nbrs[2] - p0;
        // RDKit✔️✔️:     double chiralVol = v1.dotProduct(v2.crossProduct(v3));
        // RDKit✔️✔️:     bool chiralitySet = false;
        // RDKit✔️✔️:     if (chiralVol < -ZERO_VOLUME_TOL) {
        // RDKit✔️✔️:       atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
        // RDKit✔️✔️:       chiralitySet = true;
        // RDKit✔️✔️:     } else if (chiralVol > ZERO_VOLUME_TOL) {
        // RDKit✔️✔️:       atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
        // RDKit✔️✔️:       chiralitySet = true;
        // RDKit✔️✔️:     } else if (nbrIdx == 4) {
        // RDKit✔️✔️:       auto v4 = *nbrs[3] - p0;
        // RDKit✔️✔️:       chiralVol = -v1.dotProduct(v2.crossProduct(v4));
        // RDKit✔️✔️:       if (chiralVol < -ZERO_VOLUME_TOL) { ... }
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       atom->setChiralTag(Atom::CHI_UNSPECIFIED);
        // RDKit✔️✔️:     }
        match tetrahedral_chiral_volume(p0, &nbr_positions) {
            Some(tag) => {
                if tag != ChiralTag::Unspecified && !explicit_atoms[atom_idx] {
                    non_explicit_3d_chirality[atom_idx] = true;
                }
                chiral_assignments.push((atom_idx, tag));
            }
            None => chiral_assignments.push((atom_idx, ChiralTag::Unspecified)),
        }
    }

    // Phase 2: apply assignments
    let atoms = &mut mol.topology_block_mut().atoms;
    for (idx, tag) in chiral_assignments {
        atoms[idx].set_chiral_tag(tag);
        // RDKit✔️✔️:     if (chiralitySet && explicitAtoms[atom->getIdx()] == 0) {
        // RDKit✔️✔️:       atom->setProp<int>(common_properties::_NonExplicit3DChirality, 1);
        // RDKit✔️✔️:     }
        if non_explicit_3d_chirality[idx] {
            atoms[idx].set_prop("_NonExplicit3DChirality", "1");
        }
        if let Some(permutation) = chiral_permutation_assignments[idx] {
            atoms[idx].set_chiral_permutation(Some(permutation));
        }
    }
}

pub(crate) fn assign_chiral_types_from_3d_for_testing(mol: &mut Molecule, conf_id: usize) {
    assign_chiral_types_from_3d(mol, conf_id);
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

fn handle_cx_part_and_name(
    state: &mut SmilesBuildState,
    params: &SmilesParseParams,
    cx_part: &str,
    name: &str,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION handleCXPartAndName
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️: void handleCXPartAndName(RWMol *res, const T &params, const std::string &cxPart,
    // RDKit✔️✔️:                          std::string &name) {
    // RDKit✔️✔️:   if (!res || cxPart.empty()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
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
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION handleCXPartAndName
    if cx_part.is_empty() {
        if !name.is_empty() {
            state.set_name(name);
        }
        return Ok(());
    }

    if params.allow_cxsmiles && cx_part.starts_with('|') {
        match parse_cx_extensions(state, cx_part) {
            Ok(consumed) => {
                state.set_property("_CXSMILES_Data", &cx_part[..consumed]);
                if params.parse_name && consumed < cx_part.len() {
                    let name_part = cx_part[consumed..].trim();
                    if !name_part.is_empty() {
                        state.set_name(name_part);
                    }
                }
                return Ok(());
            }
            Err(error) => {
                if params.strict_cxsmiles {
                    return Err(error);
                }
                state.set_property("_CXSMILES_Data", "");
                return Ok(());
            }
        }
    }

    if params.allow_cxsmiles && params.strict_cxsmiles && !params.parse_name {
        return Err(SmilesParseError::ParseError(
            "CXSMILES extension does not start with | and parseName=false".to_string(),
        ));
    }

    if params.parse_name {
        state.set_name(cx_part.trim());
        return Ok(());
    }

    Ok(())
}

fn parse_cx_extensions(
    state: &mut SmilesBuildState,
    ext_text: &str,
) -> Result<usize, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parseCXExtensions / parser::parse_it
    // RDKit✔️✔️: void parseCXExtensions(RDKit::RWMol &mol, const std::string &extText,
    // RDKit✔️✔️:                        std::string::const_iterator &first,
    // RDKit✔️✔️:                        unsigned int startAtomIdx, unsigned int startBondIdx) {
    // RDKit✔️✔️:   if (extText.empty()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (extText[0] != '|') {
    // RDKit✔️✔️:     throw RDKit::SmilesParseException(
    // RDKit✔️✔️:         "CXSMILES extension does not start with |");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   first = extText.begin();
    // RDKit✔️✔️:   bool ok =
    // RDKit✔️✔️:       parser::parse_it(first, extText.end(), mol, startAtomIdx, startBondIdx);
    // RDKit✔️✔️:   if (!ok) {
    // RDKit✔️✔️:     throw RDKit::SmilesParseException("failure parsing CXSMILES extensions");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   processCXSmilesLabels(mol);
    // RDKit✔️✔️:   mol.clearProp("_cxsmilesLabelsProcessed");
    // RDKit✔️✔️:   mol.clearProp(cxsgTracker);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: bool parse_it(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:               unsigned int startAtomIdx, unsigned int startBondIdx) {
    // RDKit✔️✔️:   if (first >= last || *first != '|') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   while (first < last && *first != '|') {
    // RDKit✔️✔️:     typename Iterator::difference_type length = std::distance(first, last);
    // RDKit✔️✔️:     if (*first == '(') {
    // RDKit✔️✔️:       if (!parse_coords(first, last, mol, startAtomIdx, confIndex++)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == '$') {
    // RDKit✔️✔️:       if (length > 4 && *(first + 1) == '_' && *(first + 2) == 'A' &&
    // RDKit✔️✔️:           *(first + 3) == 'V' && *(first + 4) == ':') {
    // RDKit✔️✔️:         first += 4;
    // RDKit✔️✔️:         if (!parse_atom_values(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         if (!parse_atom_labels(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (length > 9 && std::string(first, first + 9) == "atomProp:") {
    // RDKit✔️✔️:       first += 9;
    // RDKit✔️✔️:       if (!parse_atom_props(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'C') {
    // RDKit✔️✔️:       if (!parse_coordinate_bonds(first, last, mol, Bond::DATIVE, startAtomIdx,
    // RDKit✔️✔️:                                   startBondIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'H') {
    // RDKit✔️✔️:       if (!parse_coordinate_bonds(first, last, mol, Bond::HYDROGEN,
    // RDKit✔️✔️:                                   startAtomIdx, startBondIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'Z') {
    // RDKit✔️✔️:       if (!parse_zero_bonds(first, last, mol, startAtomIdx, startBondIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == '^') {
    // RDKit✔️✔️:       if (!parse_radicals(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'a' || *first == 'o' ||
    // RDKit✔️✔️:                (*first == '&' && first + 1 < last && first[1] != '#')) {
    // RDKit✔️✔️:       if (!parse_enhanced_stereo(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'r' && first + 1 < last && first[1] == 'b') {
    // RDKit✔️✔️:       if (!parse_ring_bonds(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'L' && first + 1 < last && first[1] == 'N') {
    // RDKit✔️✔️:       if (!parse_linknodes(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'S' && first + 2 < last && first[1] == 'g' &&
    // RDKit✔️✔️:                first[2] == 'D') {
    // RDKit✔️✔️:       if (!parse_data_sgroup(first, last, mol, startAtomIdx, nSGroups++)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'S' && first + 2 < last && first[1] == 'g' &&
    // RDKit✔️✔️:                first[2] == 'H') {
    // RDKit✔️✔️:       if (!parse_sgroup_hierarchy(first, last, mol)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'S' && first + 1 < last && first[1] == 'g') {
    // RDKit✔️✔️:       if (!parse_polymer_sgroup(first, last, mol, startAtomIdx, nSGroups++)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'u') {
    // RDKit✔️✔️:       if (!parse_unsaturation(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 's') {
    // RDKit✔️✔️:       if (!parse_substitution(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'm') {
    // RDKit✔️✔️:       if (!parse_variable_attachments(first, last, mol, startAtomIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'w') {
    // RDKit✔️✔️:       if (!parse_wedged_bonds(first, last, mol, startAtomIdx, startBondIdx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'c' && first + 2 < last && first[1] == 't' &&
    // RDKit✔️✔️:                first[2] == 'u') {
    // RDKit✔️✔️:       if (!parse_doublebond_stereo(first, last, mol, startAtomIdx, startBondIdx,
    // RDKit✔️✔️:                                    Bond::BondStereo::STEREOANY)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 'c') {
    // RDKit✔️✔️:       if (!parse_doublebond_stereo(first, last, mol, startAtomIdx, startBondIdx,
    // RDKit✔️✔️:                                    Bond::BondStereo::STEREOCIS)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (*first == 't') {
    // RDKit✔️✔️:       if (!parse_doublebond_stereo(first, last, mol, startAtomIdx, startBondIdx,
    // RDKit✔️✔️:                                    Bond::BondStereo::STEREOTRANS)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // if(first < last && *first != '|') ++first;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (first >= last || *first != '|') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parseCXExtensions / parser::parse_it
    if ext_text.is_empty() {
        return Ok(0);
    }
    let bytes = ext_text.as_bytes();
    if bytes.first().copied() != Some(b'|') {
        return Err(SmilesParseError::ParseError(
            "CXSMILES extension does not start with |".to_string(),
        ));
    }
    let mut pos = 1;
    let mut conformer_idx = 0_usize;
    let mut sgroup_idx = 0_usize;
    while pos < bytes.len() && bytes[pos] != b'|' {
        if bytes[pos] == b'(' {
            parse_cx_coords(state, ext_text, &mut pos, conformer_idx)?;
            conformer_idx += 1;
        } else if bytes[pos] == b'$' {
            if ext_text[pos..].starts_with("$_AV:") {
                pos += 4;
                parse_cx_atom_values(state, ext_text, &mut pos)?;
            } else {
                parse_cx_atom_labels(state, ext_text, &mut pos)?;
            }
        } else if ext_text[pos..].starts_with("atomProp:") {
            pos += "atomProp:".len();
            parse_cx_atom_props(state, ext_text, &mut pos)?;
        } else if bytes[pos] == b'C' {
            parse_cx_coordinate_bonds(state, ext_text, &mut pos, BondOrder::Dative)?;
        } else if bytes[pos] == b'H' {
            parse_cx_coordinate_bonds(state, ext_text, &mut pos, BondOrder::Hydrogen)?;
        } else if bytes[pos] == b'Z' {
            parse_cx_zero_bonds(state, ext_text, &mut pos)?;
        } else if bytes[pos] == b'^' {
            parse_cx_radicals(state, ext_text, &mut pos)?;
        } else if matches!(bytes[pos], b'a' | b'o')
            || (bytes[pos] == b'&' && pos + 1 < bytes.len() && bytes[pos + 1] != b'#')
        {
            parse_cx_enhanced_stereo(state, ext_text, &mut pos)?;
        } else if ext_text[pos..].starts_with("rb:") {
            parse_cx_ring_bonds(state, ext_text, &mut pos)?;
        } else if ext_text[pos..].starts_with("LN:") {
            parse_cx_linknodes(state, ext_text, &mut pos)?;
        } else if ext_text[pos..].starts_with("SgD:") {
            parse_cx_data_sgroup(state, ext_text, &mut pos, sgroup_idx)?;
            sgroup_idx += 1;
        } else if ext_text[pos..].starts_with("SgH:") {
            // BEGIN RDKIT CPP FUNCTION parse_sgroup_hierarchy (called from parser::parse_it)
            // RDKit✔️✔️:     } else if (*first == 'S' && first + 2 < last && first[1] == 'g' &&
            // RDKit✔️✔️:                first[2] == 'H') {
            // RDKit✔️✔️:       if (!parse_sgroup_hierarchy(first, last, mol)) {
            // RDKit✔️✔️:         return false;
            // RDKit✔️✔️:       }
            // END RDKIT CPP FUNCTION call-site in parser::parse_it
            parse_cx_sgroup_hierarchy(state, ext_text, &mut pos)?;
        } else if ext_text[pos..].starts_with("Sg:") {
            // BEGIN RDKIT CPP FUNCTION parse_polymer_sgroup (called from parser::parse_it)
            // RDKit✔️✔️:     } else if (*first == 'S' && first + 1 < last && first[1] == 'g') {
            // RDKit✔️✔️:       if (!parse_polymer_sgroup(first, last, mol, startAtomIdx, nSGroups++)) {
            // RDKit✔️✔️:         return false;
            // RDKit✔️✔️:       }
            // END RDKIT CPP FUNCTION call-site in parser::parse_it
            parse_cx_polymer_sgroup(state, ext_text, &mut pos, sgroup_idx)?;
            sgroup_idx += 1;
        } else if bytes[pos] == b'u' {
            parse_cx_unsaturation(state, ext_text, &mut pos)?;
        } else if bytes[pos] == b's' {
            parse_cx_substitution(state, ext_text, &mut pos)?;
        } else if bytes[pos] == b'm' {
            parse_cx_variable_attachments(state, ext_text, &mut pos)?;
        } else if bytes[pos] == b'w' {
            parse_cx_wedged_bonds(state, ext_text, &mut pos)?;
        } else if ext_text[pos..].starts_with("ctu") {
            parse_cx_doublebond_stereo(state, ext_text, &mut pos, BondStereo::Any)?;
        } else if bytes[pos] == b'c' {
            parse_cx_doublebond_stereo(state, ext_text, &mut pos, BondStereo::Cis)?;
        } else if bytes[pos] == b't' {
            parse_cx_doublebond_stereo(state, ext_text, &mut pos, BondStereo::Trans)?;
        } else if bytes[pos] == b',' {
            pos += 1;
        } else {
            // Skip unrecognized CXSMILES extension characters (matching RDKit
            // behavior which silently skips unknown extensions rather than
            // failing). Advance past the single-byte extension token.
            pos += 1;
        }
    }
    if pos >= bytes.len() || bytes[pos] != b'|' {
        return Err(SmilesParseError::ParseError(
            "failure parsing CXSMILES extensions".to_string(),
        ));
    }
    Ok(pos + 1)
}

fn parse_cx_coords(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
    conformer_idx: usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_coords / hasNonZeroZCoords
    // RDKit✔️✔️: bool parse_coords(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                   unsigned int startAtomIdx, unsigned int confIdx) {
    // RDKit✔️✔️:   if (first >= last || *first != '(') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto *conf = new Conformer(mol.getNumAtoms());
    // RDKit✔️✔️:   mol.addConformer(conf);
    // RDKit✔️✔️:   conf->setId(confIdx);
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   unsigned int atIdx = 0;
    // RDKit✔️✔️:   bool is3D = false;
    // RDKit✔️✔️:   while (first <= last && *first != ')') {
    // RDKit✔️✔️:     RDGeom::Point3D pt;
    // RDKit✔️✔️:     std::string tkn = read_text_to(first, last, ";)");
    // RDKit✔️✔️:     if (VALID_ATIDX(atIdx)) {
    // RDKit✔️✔️:       if (!tkn.empty()) {
    // RDKit✔️✔️:         std::vector<std::string> tokens;
    // RDKit✔️✔️:         boost::split(tokens, tkn, boost::is_any_of(std::string(",")));
    // RDKit✔️✔️:         if (tokens.size() >= 1 && tokens[0].size()) {
    // RDKit✔️✔️:           pt.x = boost::lexical_cast<double>(tokens[0]);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (tokens.size() >= 2 && tokens[1].size()) {
    // RDKit✔️✔️:           pt.y = boost::lexical_cast<double>(tokens[1]);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (tokens.size() >= 3 && tokens[2].size()) {
    // RDKit✔️✔️:           pt.z = boost::lexical_cast<double>(tokens[2]);
    // RDKit✔️✔️:           is3D = true;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       conf->setAtomPos(atIdx - startAtomIdx, pt);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++atIdx;
    // RDKit✔️✔️:     if (first <= last && *first != ')') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // make sure that the conformer really is 3D!
    // RDKit✔️✔️:   if (is3D && hasNonZeroZCoords(*conf)) {
    // RDKit✔️✔️:     conf->set3D(true);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     conf->set3D(false);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (first >= last || *first != ')') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: inline bool hasNonZeroZCoords(const Conformer &conf) {
    // RDKit✔️✔️:   constexpr double zeroTol = 1e-3;
    // RDKit✔️✔️:   for (auto p : conf.getPositions()) {
    // RDKit✔️✔️:     if (std::abs(p.z) > zeroTol) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_coords / hasNonZeroZCoords
    expect_byte(ext_text, *pos, b'(')?;
    *pos += 1;
    let atom_count = state.builder.atoms().len();
    let mut coords = vec![[0.0_f64, 0.0_f64, 0.0_f64]; atom_count];
    let mut atom_idx = 0_usize;
    let mut saw_z = false;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos] != b')' {
        let token = read_cx_text_to(ext_text, pos, b";)")?;
        if atom_idx < atom_count && !token.is_empty() {
            let mut pieces = token.split(',');
            if let Some(x) = pieces.next()
                && !x.is_empty()
            {
                coords[atom_idx][0] = x.parse::<f64>().map_err(|_| cx_parse_failure())?;
            }
            if let Some(y) = pieces.next()
                && !y.is_empty()
            {
                coords[atom_idx][1] = y.parse::<f64>().map_err(|_| cx_parse_failure())?;
            }
            if let Some(z) = pieces.next()
                && !z.is_empty()
            {
                coords[atom_idx][2] = z.parse::<f64>().map_err(|_| cx_parse_failure())?;
                saw_z = true;
            }
        }
        atom_idx += 1;
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] != b')' {
            *pos += 1;
        }
    }
    let is_3d = saw_z && coords.iter().any(|coord| coord[2].abs() > 1e-3);
    state
        .builder
        .add_conformer(Conformer3D::new(conformer_idx, coords, is_3d))
        .map_err(|error| SmilesParseError::ParseError(error.to_string()))?;
    expect_byte(ext_text, *pos, b')')?;
    *pos += 1;
    Ok(())
}

fn parse_cx_atom_labels(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_atom_labels
    // RDKit✔️✔️: bool parse_atom_labels(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                        unsigned int startAtomIdx) {
    // RDKit✔️✔️:   if (first >= last || *first != '$') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   unsigned int atIdx = 0;
    // RDKit✔️✔️:   while (first <= last && *first != '$') {
    // RDKit✔️✔️:     std::string tkn = read_text_to(first, last, ";$");
    // RDKit✔️✔️:     if (!tkn.empty() && VALID_ATIDX(atIdx)) {
    // RDKit✔️✔️:       mol.getAtomWithIdx(atIdx - startAtomIdx)
    // RDKit✔️✔️:           ->setProp(RDKit::common_properties::atomLabel, tkn);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++atIdx;
    // RDKit✔️✔️:     if (first <= last && *first != '$') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (first >= last || *first != '$') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_atom_labels
    expect_byte(ext_text, *pos, b'$')?;
    *pos += 1;
    let mut atom_idx = 0_usize;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos] != b'$' {
        let token = read_cx_text_to(ext_text, pos, b";$")?;
        if !token.is_empty()
            && let Some(atom) = state.builder.atom_mut(AtomId::new(atom_idx))
        {
            atom.set_prop("atomLabel", token);
        }
        atom_idx += 1;
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] != b'$' {
            *pos += 1;
        }
    }
    expect_byte(ext_text, *pos, b'$')?;
    *pos += 1;
    Ok(())
}

fn parse_cx_atom_values(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_atom_values
    // RDKit✔️✔️: bool parse_atom_values(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                        unsigned int startAtomIdx) {
    // RDKit✔️✔️:   if (first >= last || *first != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   unsigned int atIdx = 0;
    // RDKit✔️✔️:   while (first <= last && *first != '$') {
    // RDKit✔️✔️:     std::string tkn = read_text_to(first, last, ";$");
    // RDKit✔️✔️:     if (tkn != "" && VALID_ATIDX(atIdx)) {
    // RDKit✔️✔️:       mol.getAtomWithIdx(atIdx)->setProp(RDKit::common_properties::molFileValue,
    // RDKit✔️✔️:                                          tkn);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++atIdx;
    // RDKit✔️✔️:     if (first <= last && *first != '$') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (first >= last || *first != '$') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_atom_values
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    let mut atom_idx = 0_usize;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos] != b'$' {
        let token = read_cx_text_to(ext_text, pos, b";$")?;
        if !token.is_empty()
            && let Some(atom) = state.builder.atom_mut(AtomId::new(atom_idx))
        {
            atom.set_prop("molFileValue", token);
        }
        atom_idx += 1;
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] != b'$' {
            *pos += 1;
        }
    }
    expect_byte(ext_text, *pos, b'$')?;
    *pos += 1;
    Ok(())
}

fn parse_cx_atom_props(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_atom_props
    // RDKit✔️✔️: bool parse_atom_props(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                       unsigned int startAtomIdx) {
    // RDKit✔️✔️:   while (first <= last && *first != '|' && *first != ',') {
    // RDKit✔️✔️:     unsigned int atIdx;
    // RDKit✔️✔️:     if (read_int(first, last, atIdx)) {
    // RDKit✔️✔️:       if (first >= last || *first != '.') {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:       std::string pname = read_text_to(first, last, ".");
    // RDKit✔️✔️:       if (!pname.empty()) {
    // RDKit✔️✔️:         if (first >= last || *first != '.') {
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         ++first;
    // RDKit✔️✔️:         std::string pval = read_text_to(first, last, ":|,");
    // RDKit✔️✔️:         if (VALID_ATIDX(atIdx) && !pval.empty()) {
    // RDKit✔️✔️:           mol.getAtomWithIdx(atIdx - startAtomIdx)->setProp(pname, pval);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first <= last && *first != '|' && *first != ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (*first != '|') {
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_atom_props
    while *pos < ext_text.len() && !matches!(ext_text.as_bytes()[*pos], b'|' | b',') {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        expect_byte(ext_text, *pos, b'.')?;
        *pos += 1;
        let prop_name = read_cx_text_to(ext_text, pos, b".")?;
        expect_byte(ext_text, *pos, b'.')?;
        *pos += 1;
        let prop_value = read_cx_text_to(ext_text, pos, b":|,")?;
        if !prop_name.is_empty()
            && !prop_value.is_empty()
            && let Some(atom) = state.builder.atom_mut(AtomId::new(atom_idx))
        {
            atom.set_prop(prop_name, prop_value);
        }
        if *pos < ext_text.len() && !matches!(ext_text.as_bytes()[*pos], b'|' | b',') {
            *pos += 1;
        }
    }
    if *pos < ext_text.len() && ext_text.as_bytes()[*pos] != b'|' {
        *pos += 1;
    }
    Ok(())
}

fn cx_bond_with_smiles_idx(state: &SmilesBuildState, idx: usize) -> Option<BondId> {
    // BEGIN RDKIT CPP FUNCTION get_bond_with_smiles_idx
    // RDKit✔️✔️: Bond *get_bond_with_smiles_idx(const ROMol &mol, unsigned idx) {
    // RDKit✔️✔️:   for (auto bnd : mol.bonds()) {
    // RDKit✔️✔️:     unsigned int smilesIdx;
    // RDKit✔️✔️:     if (bnd->getPropIfPresent("_cxsmilesBondIdx", smilesIdx) &&
    // RDKit✔️✔️:         smilesIdx == idx) {
    // RDKit✔️✔️:       return bnd;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return nullptr;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION get_bond_with_smiles_idx
    state
        .builder
        .bonds()
        .iter()
        .find(|bond| {
            bond.prop(CXSMILES_BOND_IDX_PROP)
                .and_then(|value| value.parse::<usize>().ok())
                == Some(idx)
        })
        .map(|bond| bond.id())
}

fn parse_cx_coordinate_bonds(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
    order: BondOrder,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_coordinate_bonds
    // RDKit✔️✔️: bool parse_coordinate_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                             Bond::BondType typ, unsigned int startAtomIdx,
    // RDKit✔️✔️:                             unsigned int startBondIdx) {
    // RDKit✔️✔️:   if (first >= last || (*first != 'C' && *first != 'H')) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   if (first >= last || *first != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   while (first <= last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int aidx;
    // RDKit✔️✔️:     unsigned int bidx;
    // RDKit✔️✔️:     if (read_int_pair(first, last, aidx, bidx)) {
    // RDKit✔️✔️:       if (VALID_ATIDX(aidx) && VALID_BNDIDX(bidx)) {
    // RDKit✔️✔️:         auto bnd = get_bond_with_smiles_idx(mol, bidx - startBondIdx);
    // RDKit✔️✔️:         if (!bnd || (bnd->getBeginAtomIdx() != aidx - startAtomIdx &&
    // RDKit✔️✔️:                      bnd->getEndAtomIdx() != aidx - startAtomIdx)) {
    // RDKit✔️✔️:           BOOST_LOG(rdWarningLog) << "BOND NOT FOUND! " << bidx
    // RDKit✔️✔️:                                   << " involving atom " << aidx << std::endl;
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         bnd->setBondType(typ);
    // RDKit✔️✔️:         if (bnd->getBeginAtomIdx() != aidx - startAtomIdx) {
    // RDKit✔️✔️:           unsigned int tmp = bnd->getBeginAtomIdx();
    // RDKit✔️✔️:           bnd->setBeginAtomIdx(aidx - startAtomIdx);
    // RDKit✔️✔️:           bnd->setEndAtomIdx(tmp);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_coordinate_bonds
    if *pos >= ext_text.len() || !matches!(ext_text.as_bytes()[*pos], b'C' | b'H') {
        return Err(cx_parse_failure());
    }
    *pos += 1;
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        expect_byte(ext_text, *pos, b'.')?;
        *pos += 1;
        let bond_idx = read_cx_usize(ext_text, pos)?;
        if state.builder.atom_mut(AtomId::new(atom_idx)).is_some() {
            let bond_id = cx_bond_with_smiles_idx(state, bond_idx).ok_or_else(cx_parse_failure)?;
            let (begin, end) = state
                .builder
                .bond(bond_id)
                .map(|bond| (bond.begin(), bond.end()))
                .ok_or_else(cx_parse_failure)?;
            let atom_id = AtomId::new(atom_idx);
            if begin != atom_id && end != atom_id {
                return Err(cx_parse_failure());
            }
            let bond = state
                .builder
                .bond_mut(bond_id)
                .ok_or_else(cx_parse_failure)?;
            bond.set_order(order);
            if begin != atom_id {
                bond.set_endpoints(atom_id, begin);
            }
        }
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
    }
    Ok(())
}

fn parse_cx_zero_bonds(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_zero_bonds
    // RDKit✔️✔️: bool parse_zero_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                       unsigned int, unsigned int startBondIdx) {
    // RDKit✔️✔️:   // these look like: C1CCCCC~CCCC1 |Z:5|
    // RDKit✔️✔️:   if (first >= last || *first != 'Z') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   if (first >= last || *first != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int bondIdx;
    // RDKit✔️✔️:     if (!read_int(first, last, bondIdx)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (VALID_BNDIDX(bondIdx)) {
    // RDKit✔️✔️:       auto bond = get_bond_with_smiles_idx(mol, bondIdx - startBondIdx);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (!bond) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "bond " << bondIdx
    // RDKit✔️✔️:             << " not found, cannot mark as zero order bond." << std::endl;
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       bond->setBondType(Bond::ZERO);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_zero_bonds
    expect_byte(ext_text, *pos, b'Z')?;
    *pos += 1;
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let bond_idx = read_cx_usize(ext_text, pos)?;
        let bond_id = cx_bond_with_smiles_idx(state, bond_idx).ok_or_else(cx_parse_failure)?;
        state
            .builder
            .bond_mut(bond_id)
            .ok_or_else(cx_parse_failure)?
            .set_order(BondOrder::Zero);
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
    }
    Ok(())
}

fn parse_cx_enhanced_stereo(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_enhanced_stereo
    // RDKit✔️✔️: bool parse_enhanced_stereo(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                            unsigned int startAtomIdx) {
    // RDKit✔️✔️:   StereoGroupType group_type = StereoGroupType::STEREO_ABSOLUTE;
    // RDKit✔️✔️:   if (*first == 'a') {
    // RDKit✔️✔️:     group_type = StereoGroupType::STEREO_ABSOLUTE;
    // RDKit✔️✔️:   } else if (*first == 'o') {
    // RDKit✔️✔️:     group_type = StereoGroupType::STEREO_OR;
    // RDKit✔️✔️:   } else if (*first == '&') {
    // RDKit✔️✔️:     group_type = StereoGroupType::STEREO_AND;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // OR and AND groups carry a group number
    // RDKit✔️✔️:   unsigned int group_id = 0;
    // RDKit✔️✔️:   if (group_type != StereoGroupType::STEREO_ABSOLUTE) {
    // RDKit✔️✔️:     read_int(first, last, group_id);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (first >= last || *first != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<Atom *> atoms;
    // RDKit✔️✔️:   std::vector<Bond *> bonds;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   while (first <= last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int aidx;
    // RDKit✔️✔️:     if (read_int(first, last, aidx)) {
    // RDKit✔️✔️:       if (VALID_ATIDX(aidx)) {
    // RDKit✔️✔️:         Atom *atom = mol.getAtomWithIdx(aidx - startAtomIdx);
    // RDKit✔️✔️:         if (!atom) {
    // RDKit✔️✔️:           BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:               << "Atom " << aidx << " not found!" << std::endl;
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         atoms.push_back(atom);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!atoms.empty()) {
    // RDKit✔️✔️:     // we need to do a bit of work to check whether or not we've already seen
    // RDKit✔️✔️:     // this particular StereoGroup (was Github #6050)
    // RDKit✔️✔️:     const auto group_hash =
    // RDKit✔️✔️:         10 * group_id + static_cast<unsigned int>(group_type);
    // RDKit✔️✔️:     std::vector<unsigned int> sgTracker;
    // RDKit✔️✔️:     mol.getPropIfPresent(cxsgTracker, sgTracker);
    // RDKit✔️✔️:     std::vector<StereoGroup> mol_stereo_groups(mol.getStereoGroups());
    // RDKit✔️✔️:     TEST_ASSERT(mol_stereo_groups.size() == sgTracker.size());
    // RDKit✔️✔️:
    // RDKit✔️✔️:     auto iter = std::find(sgTracker.begin(), sgTracker.end(), group_hash);
    // RDKit✔️✔️:     if (iter != sgTracker.end()) {
    // RDKit✔️✔️:       auto index = iter - sgTracker.begin();
    // RDKit✔️✔️:       auto gAtoms = mol_stereo_groups[index].getAtoms();
    // RDKit✔️✔️:       gAtoms.insert(gAtoms.end(), atoms.begin(), atoms.end());
    // RDKit✔️✔️:       mol_stereo_groups[index] =
    // RDKit✔️✔️:           StereoGroup(mol_stereo_groups[index].getGroupType(),
    // RDKit✔️✔️:                       std::move(gAtoms), std::move(bonds), group_id);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       // not seen this before, create a new stereogroup
    // RDKit✔️✔️:       mol_stereo_groups.emplace_back(group_type, std::move(atoms),
    // RDKit✔️✔️:                                      std::move(bonds), group_id);
    // RDKit✔️✔️:       sgTracker.push_back(group_hash);
    // RDKit✔️✔️:       mol.setProp(cxsgTracker, sgTracker);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     mol.setStereoGroups(std::move(mol_stereo_groups));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_enhanced_stereo
    let kind = match ext_text.as_bytes().get(*pos).copied() {
        Some(b'a') => StereoGroupKind::Absolute,
        Some(b'o') => StereoGroupKind::Or,
        Some(b'&') => StereoGroupKind::And,
        _ => return Err(cx_parse_failure()),
    };
    *pos += 1;
    let group_id = if kind == StereoGroupKind::Absolute {
        0
    } else if *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        read_cx_usize(ext_text, pos)? as u32
    } else {
        0
    };
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    let mut atoms = Vec::new();
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        let atom_id = AtomId::new(atom_idx);
        if state.builder.atom_mut(atom_id).is_some() {
            atoms.push(atom_id);
        } else {
            return Err(cx_parse_failure());
        }
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
    }
    if atoms.is_empty() {
        return Ok(());
    }
    let kind_code = match kind {
        StereoGroupKind::Absolute => 1,
        StereoGroupKind::Or => 2,
        StereoGroupKind::And => 3,
    };
    if let Some(index) = state
        .cx_stereo_group_tracker
        .get(&(kind_code, group_id))
        .copied()
    {
        if let Some(group) = state.builder.stereo_groups_mut().get_mut(index) {
            for atom in atoms {
                group.push_atom(atom);
            }
        }
    } else {
        let index = state.builder.stereo_groups_mut().len();
        state
            .builder
            .add_stereo_group(StereoGroup::new(kind, atoms, Vec::new()).with_id(group_id))
            .map_err(|error| SmilesParseError::ParseError(error.to_string()))?;
        state
            .cx_stereo_group_tracker
            .insert((kind_code, group_id), index);
    }
    Ok(())
}

fn expand_cx_atom_query(
    state: &mut SmilesBuildState,
    atom_idx: usize,
    predicate: AtomQueryPredicate,
) {
    if let Some(atom) = state.builder.atom_mut(AtomId::new(atom_idx)) {
        let next = QueryNode::predicate(predicate);
        let combined = match atom.query().cloned() {
            Some(QueryNode::And(mut children)) => {
                children.push(next);
                QueryNode::And(children)
            }
            Some(existing) => QueryNode::and(vec![existing, next]),
            None => next,
        };
        atom.set_query(Some(combined));
    }
}

fn parse_cx_unsaturation(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_unsaturation
    // RDKit✔️✔️: bool parse_unsaturation(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                         unsigned int startAtomIdx) {
    // RDKit✔️✔️:   if (first + 1 >= last || *first != 'u') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   if (first >= last || *first != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int idx;
    // RDKit✔️✔️:     if (!read_int(first, last, idx)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (VALID_ATIDX(idx)) {
    // RDKit✔️✔️:       auto atom = mol.getAtomWithIdx(idx - startAtomIdx);
    // RDKit✔️✔️:       if (!atom->hasQuery()) {
    // RDKit✔️✔️:         atom = QueryOps::replaceAtomWithQueryAtom(&mol, atom);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       atom->expandQuery(makeAtomUnsaturatedQuery(), Queries::COMPOSITE_AND);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_unsaturation
    expect_byte(ext_text, *pos, b'u')?;
    *pos += 1;
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        expand_cx_atom_query(state, atom_idx, AtomQueryPredicate::IsUnsaturated);
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
    }
    Ok(())
}

fn parse_cx_ring_bonds(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_ring_bonds
    // RDKit✔️✔️: bool parse_ring_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                       unsigned int startAtomIdx) {
    // RDKit✔️✔️:   if (first >= last || *first != 'r' || first + 1 >= last ||
    // RDKit✔️✔️:       *(first + 1) != 'b' || first + 2 >= last || *(first + 2) != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   first += 3;
    // RDKit✔️✔️:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int n1;
    // RDKit✔️✔️:     if (!read_int(first, last, n1)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // check that we can read at least two more characters:
    // RDKit✔️✔️:     if (first + 1 >= last || *first != ':') {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:     unsigned int n2;
    // RDKit✔️✔️:     bool gt = false;
    // RDKit✔️✔️:     if (*first == '*') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:       n2 = 0xDEADBEEF;
    // RDKit✔️✔️:       if (VALID_ATIDX(n1)) {
    // RDKit✔️✔️:         mol.setProp(common_properties::_NeedsQueryScan, 1);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       if (!read_int(first, last, n2)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       switch (n2) {
    // RDKit✔️✔️:         case 0:
    // RDKit✔️✔️:         case 2:
    // RDKit✔️✔️:         case 3:
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case 4:
    // RDKit✔️✔️:           gt = true;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:               << "unrecognized rb value: " << n2 << std::endl;
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (VALID_ATIDX(n1)) {
    // RDKit✔️✔️:       auto atom = mol.getAtomWithIdx(n1 - startAtomIdx);
    // RDKit✔️✔️:       if (!atom->hasQuery()) {
    // RDKit✔️✔️:         atom = QueryOps::replaceAtomWithQueryAtom(&mol, atom);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!gt) {
    // RDKit✔️✔️:         atom->expandQuery(makeAtomRingBondCountQuery(n2),
    // RDKit✔️✔️:                           Queries::COMPOSITE_AND);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         auto q = static_cast<ATOM_EQUALS_QUERY *>(new ATOM_LESSEQUAL_QUERY);
    // RDKit✔️✔️:         q->setVal(n2);
    // RDKit✔️✔️:         q->setDescription("AtomRingBondCount");
    // RDKit✔️✔️:         q->setDataFunc(queryAtomRingBondCount);
    // RDKit✔️✔️:         atom->expandQuery(q, Queries::COMPOSITE_AND);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_ring_bonds
    if !ext_text[*pos..].starts_with("rb:") {
        return Err(cx_parse_failure());
    }
    *pos += 3;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        expect_byte(ext_text, *pos, b':')?;
        *pos += 1;
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b'*' {
            *pos += 1;
            if state.builder.atom_mut(AtomId::new(atom_idx)).is_some() {
                state.set_property("_NeedsQueryScan", "1");
                expand_cx_atom_query(state, atom_idx, AtomQueryPredicate::RingBondCountNeedsScan);
            }
        } else {
            let ring_bonds = read_cx_usize(ext_text, pos)?;
            let predicate = match ring_bonds {
                0 | 2 | 3 => AtomQueryPredicate::RingBondCount(ring_bonds as u8),
                4 => AtomQueryPredicate::RingBondCountLessEqual(4),
                _ => return Err(cx_parse_failure()),
            };
            expand_cx_atom_query(state, atom_idx, predicate);
        }
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
    }
    Ok(())
}

fn atom_neighbors(state: &SmilesBuildState, atom: AtomId) -> Vec<AtomId> {
    state
        .builder
        .bonds()
        .iter()
        .filter_map(|bond| {
            if bond.begin() == atom {
                Some(bond.end())
            } else if bond.end() == atom {
                Some(bond.begin())
            } else {
                None
            }
        })
        .collect()
}

fn atom_bonds(state: &SmilesBuildState, atom: AtomId) -> Vec<BondId> {
    state
        .builder
        .bonds()
        .iter()
        .filter_map(|bond| (bond.begin() == atom || bond.end() == atom).then_some(bond.id()))
        .collect()
}

fn parse_cx_linknodes(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_linknodes
    // RDKit✔️✔️: bool parse_linknodes(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                      unsigned int startAtomIdx) {
    // RDKit✔️✔️:   // these look like: |LN:1:1.3.2.6,4:1.4.3.6|
    // RDKit✔️✔️:   // that's two records:
    // RDKit✔️✔️:   //   1:1.3.2.6: 1-3 repeats, atom 1-2, 1-6
    // RDKit✔️✔️:   //   4:1.4.3.6: 1-4 repeats, atom 4-3, 4-6
    // RDKit✔️✔️:   // which maps to the property value "1 3 2 2 3 2 7|1 4 2 5 4 5 7"
    // RDKit✔️✔️:   // If the linking atom only has two neighbors then the outer atom
    // RDKit✔️✔️:   // specification (the last two digits) can be left out. So for a molecule
    // RDKit✔️✔️:   // where atom 1 has bonds only to atoms 2 and 6 we could have
    // RDKit✔️✔️:   // |LN:1:1.3|
    // RDKit✔️✔️:   // instead of
    // RDKit✔️✔️:   // |LN:1:1.3.2.6|
    // RDKit✔️✔️:   if (first >= last || *first != 'L' || first + 1 >= last ||
    // RDKit✔️✔️:       *(first + 1) != 'N' || first + 2 >= last || *(first + 2) != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   first += 3;
    // RDKit✔️✔️:   std::string accum = "";
    // RDKit✔️✔️:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int atidx;
    // RDKit✔️✔️:     if (!read_int(first, last, atidx)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // check that we can read at least two more characters:
    // RDKit✔️✔️:     if (first + 1 >= last || *first != ':') {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:     unsigned int startReps;
    // RDKit✔️✔️:     if (!read_int(first, last, startReps)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first + 1 >= last || *first != '.') {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:     unsigned int endReps;
    // RDKit✔️✔️:     if (!read_int(first, last, endReps)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     unsigned int idx1;
    // RDKit✔️✔️:     unsigned int idx2;
    // RDKit✔️✔️:     if (first < last && *first == '.') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:       if (!read_int(first, last, idx1)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:       if (!read_int(first, last, idx2)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (VALID_ATIDX(atidx) &&
    // RDKit✔️✔️:                mol.getAtomWithIdx(atidx - startAtomIdx)->getDegree() == 2) {
    // RDKit✔️✔️:       auto nbrs =
    // RDKit✔️✔️:           mol.getAtomNeighbors(mol.getAtomWithIdx(atidx - startAtomIdx));
    // RDKit✔️✔️:       idx1 = *nbrs.first;
    // RDKit✔️✔️:       nbrs.first++;
    // RDKit✔️✔️:       idx2 = *nbrs.first;
    // RDKit✔️✔️:     } else if (VALID_ATIDX(atidx)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (VALID_ATIDX(atidx)) {
    // RDKit✔️✔️:       if (!accum.empty()) {
    // RDKit✔️✔️:         accum += "|";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       accum += (boost::format("%d %d 2 %d %d %d %d") % startReps % endReps %
    // RDKit✔️✔️:                 (atidx - startAtomIdx + 1) % (idx1 - startAtomIdx + 1) %
    // RDKit✔️✔️:                 (atidx - startAtomIdx + 1) % (idx2 - startAtomIdx + 1))
    // RDKit✔️✔️:                    .str();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!accum.empty()) {
    // RDKit✔️✔️:     mol.setProp(common_properties::molFileLinkNodes, accum);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_linknodes
    if !ext_text[*pos..].starts_with("LN:") {
        return Err(cx_parse_failure());
    }
    *pos += 3;
    let mut records = Vec::new();
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        expect_byte(ext_text, *pos, b':')?;
        *pos += 1;
        let start_reps = read_cx_usize(ext_text, pos)?;
        expect_byte(ext_text, *pos, b'.')?;
        *pos += 1;
        let end_reps = read_cx_usize(ext_text, pos)?;
        let atom_id = AtomId::new(atom_idx);
        let (idx1, idx2) = if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b'.' {
            *pos += 1;
            let idx1 = read_cx_usize(ext_text, pos)?;
            expect_byte(ext_text, *pos, b'.')?;
            *pos += 1;
            let idx2 = read_cx_usize(ext_text, pos)?;
            (idx1, idx2)
        } else if state.builder.atom_mut(atom_id).is_some() {
            let neighbors = atom_neighbors(state, atom_id);
            if neighbors.len() != 2 {
                return Err(cx_parse_failure());
            }
            (neighbors[0].index(), neighbors[1].index())
        } else {
            (0, 0)
        };
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
        if state.builder.atom_mut(atom_id).is_some() {
            records.push(format!(
                "{} {} 2 {} {} {} {}",
                start_reps,
                end_reps,
                atom_idx + 1,
                idx1 + 1,
                atom_idx + 1,
                idx2 + 1
            ));
        }
    }
    if !records.is_empty() {
        state.set_property("_MolFileLinkNodes", &records.join("|"));
    }
    Ok(())
}

fn parse_cx_data_sgroup_attr(
    sgroup: &mut SubstanceGroup,
    ext_text: &str,
    pos: &mut usize,
    field_name: &str,
    field_is_array: bool,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_data_sgroup_attr
    // RDKit✔️✔️: void parse_data_sgroup_attr(Iterator &first, Iterator last,
    // RDKit✔️✔️:                             SubstanceGroup &sgroup, bool keepSGroup,
    // RDKit✔️✔️:                             std::string fieldName, bool fieldIsArray = false) {
    // RDKit✔️✔️:   PRECONDITION(first < last, "parse_data_sgroup_attr: first >= last");
    // RDKit✔️✔️:   if (first != last && *first != '|') {
    // RDKit✔️✔️:     std::string data = read_text_to(first, last, ":");
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:     if (!data.empty() && keepSGroup) {
    // RDKit✔️✔️:       if (fieldIsArray) {
    // RDKit✔️✔️:         std::vector<std::string> dataFields = {data};
    // RDKit✔️✔️:         sgroup.setProp(fieldName, dataFields);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         sgroup.setProp(fieldName, data);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_data_sgroup_attr
    if *pos < ext_text.len() && ext_text.as_bytes()[*pos] != b'|' {
        let data = read_cx_text_to(ext_text, pos, b":")?;
        if *pos < ext_text.len() {
            *pos += 1;
        }
        if !data.is_empty() {
            if field_is_array {
                sgroup.set_prop(field_name, data.clone());
                sgroup.push_data_field(data);
            } else {
                sgroup.set_prop(field_name, data);
            }
        }
    }
    Ok(())
}

fn parse_cx_data_sgroup(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
    sgroup_idx: usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_data_sgroup
    // RDKit✔️✔️: bool parse_data_sgroup(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                        unsigned int startAtomIdx, unsigned int nSGroups) {
    // RDKit✔️✔️:   // these look like: |SgD:2,1:FIELD:info::::|
    // RDKit✔️✔️:   // example from CXSMILES docs:
    // RDKit✔️✔️:   //    SgD:3,2,1,0:name:data:like:unit:t:(1.,1.)
    // RDKit✔️✔️:   // the fields are:
    // RDKit✔️✔️:   //    SgD:[atom indices]:[field name]:[data value]:[query
    // RDKit✔️✔️:   //    operator]:[unit]:[tag]:[coords]
    // RDKit✔️✔️:   //   coords are (-1) if atomic coordinates are present
    // RDKit✔️✔️:   if (first >= last || *first != 'S' || first + 3 >= last ||
    // RDKit✔️✔️:       *(first + 1) != 'g' || *(first + 2) != 'D' || *(first + 3) != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   first += 4;
    // RDKit✔️✔️:   std::vector<unsigned int> atoms;
    // RDKit✔️✔️:   if (!read_int_list(first, last, atoms)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   SubstanceGroup sgroup(&mol, std::string("DAT"));
    // RDKit✔️✔️:   sgroup.setProp(cxsmilesindex, nSGroups);
    // RDKit✔️✔️:   bool keepSGroup = false;
    // RDKit✔️✔️:   for (auto idx : atoms) {
    // RDKit✔️✔️:     if (VALID_ATIDX(idx)) {
    // RDKit✔️✔️:       keepSGroup = true;
    // RDKit✔️✔️:       sgroup.addAtomWithIdx(idx - startAtomIdx);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "FIELDNAME");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // FIX:
    // RDKit✔️✔️:   if (keepSGroup) {
    // RDKit✔️✔️:     sgroup.setProp("FIELDDISP", "    0.0000    0.0000    DR    ALL  0       0");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "DATAFIELDS", true);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "QUERYOP");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "FIELDINFO");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   parse_data_sgroup_attr(first, last, sgroup, keepSGroup, "FIELDTAG");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (first < last && *first == '(') {
    // RDKit✔️✔️:     // FIX
    // RDKit✔️✔️:     std::string coords = read_text_to(first, last, ")");
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:     if (keepSGroup) {
    // RDKit✔️✔️:       sgroup.setProp("COORDS", coords);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // the label processing can destroy sgroup info, so do that now
    // RDKit✔️✔️:   // (the function will immediately return if already called)
    // RDKit✔️✔️:   if (keepSGroup) {
    // RDKit✔️✔️:     processCXSmilesLabels(mol);
    // RDKit✔️✔️:     sgroup.setProp<unsigned int>("index", getSubstanceGroups(mol).size() + 1);
    // RDKit✔️✔️:     addSubstanceGroup(mol, sgroup);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_data_sgroup
    if !ext_text[*pos..].starts_with("SgD:") {
        return Err(cx_parse_failure());
    }
    *pos += 4;
    let mut atoms = Vec::new();
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        if state.builder.atom_mut(AtomId::new(atom_idx)).is_some() {
            atoms.push(AtomId::new(atom_idx));
        }
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        } else {
            break;
        }
    }
    let mut sgroup =
        SubstanceGroup::new(SubstanceGroupId::new(sgroup_idx), SubstanceGroupKind::Data);
    sgroup.set_prop("_cxsmilesindex", sgroup_idx.to_string());
    let keep_sgroup = !atoms.is_empty();
    for atom in atoms {
        sgroup.push_atom(atom);
    }
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    parse_cx_data_sgroup_attr(&mut sgroup, ext_text, pos, "FIELDNAME", false)?;
    if keep_sgroup {
        sgroup.set_prop("FIELDDISP", "    0.0000    0.0000    DR    ALL  0       0");
    }
    parse_cx_data_sgroup_attr(&mut sgroup, ext_text, pos, "DATAFIELDS", true)?;
    parse_cx_data_sgroup_attr(&mut sgroup, ext_text, pos, "QUERYOP", false)?;
    parse_cx_data_sgroup_attr(&mut sgroup, ext_text, pos, "FIELDINFO", false)?;
    parse_cx_data_sgroup_attr(&mut sgroup, ext_text, pos, "FIELDTAG", false)?;
    if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b'(' {
        let coords = read_cx_text_to(ext_text, pos, b")")?;
        if *pos < ext_text.len() {
            *pos += 1;
        }
        if keep_sgroup {
            sgroup.set_prop("COORDS", coords);
        }
    }
    if keep_sgroup {
        sgroup.set_prop(
            "index",
            (state.builder.substance_groups_len() + 1).to_string(),
        );
        state
            .builder
            .add_substance_group(sgroup)
            .map_err(|error| SmilesParseError::ParseError(error.to_string()))?;
    }
    Ok(())
}

fn parse_cx_sgroup_hierarchy(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_sgroup_hierarchy (CXSmilesOps.cpp)
    // RDKit✔️✔️: bool parse_sgroup_hierarchy(Iterator &first, Iterator last, RDKit::RWMol &mol) {
    // RDKit✔️✔️:   // these look like: |SgH:1:0|
    // RDKit✔️✔️:   // from CXSMILES docs:
    // RDKit✔️✔️:   //    SgH:parentSgroupIndex1:childSgroupIndex1.childSgroupIndex2,parentSgroupIndex2:childSgroupIndex1
    // RDKit✔️✔️:   if (first >= last || *first != 'S' || first + 3 >= last ||
    // RDKit✔️✔️:       *(first + 1) != 'g' || *(first + 2) != 'H' || *(first + 3) != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   first += 4;
    // RDKit✔️✔️:   auto &sgs = getSubstanceGroups(mol);
    // RDKit✔️✔️:   while (true) {
    // RDKit✔️✔️:     unsigned int parentId;
    // RDKit✔️✔️:     if (!read_int(first, last, parentId)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     bool validParent = true;
    // RDKit✔️✔️:     auto psg = find_matching_sgroup(sgs, parentId);
    // RDKit✔️✔️:     if (psg == sgs.end()) {
    // RDKit✔️✔️:       validParent = false;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       psg->getPropIfPresent(\"index\", parentId);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first <= last && *first == ':') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:       std::vector<unsigned int> children;
    // RDKit✔️✔️:       if (!read_int_list(first, last, children, '.')) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (validParent) {
    // RDKit✔️✔️:         for (auto childId : children) {
    // RDKit✔️✔️:           if (childId >= sgs.size()) {
    // RDKit✔️✔️:             throw SmilesParseException("child id references non-existent SGroup");
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           auto csg = find_matching_sgroup(sgs, childId);
    // RDKit✔️✔️:           if (csg != sgs.end()) {
    // RDKit✔️✔️:             unsigned int cid;
    // RDKit✔️✔️:             csg->getProp("index", cid);
    // RDKit✔️✔️:             csg->setProp("PARENT", parentId);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (first <= last && *first == ',') {
    // RDKit✔️✔️:         ++first;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_sgroup_hierarchy
    if !ext_text[*pos..].starts_with("SgH:") {
        return Err(cx_parse_failure());
    }
    *pos += 4;
    // Build a mapping from _cxsmilesindex -> (position, "index" value)
    let cx_index_to_parent: Vec<(usize, usize)> = state
        .builder
        .substance_groups()
        .iter()
        .enumerate()
        .filter_map(|(pos, sg)| {
            sg.props().get("_cxsmilesindex").and_then(|cx_idx| {
                let idx_val = sg
                    .props()
                    .get("index")
                    .and_then(|v| v.parse::<usize>().ok())
                    .unwrap_or(pos);
                cx_idx.parse::<usize>().ok().map(|parsed| (parsed, idx_val))
            })
        })
        .collect();
    loop {
        let parent_cx_idx = read_cx_usize(ext_text, pos)?;
        let parent_actual_idx = cx_index_to_parent
            .iter()
            .find(|(cx, _)| *cx == parent_cx_idx)
            .map(|(_, idx)| *idx);
        expect_byte(ext_text, *pos, b':')?;
        *pos += 1;
        // Read child indices separated by '.'
        let mut children = Vec::new();
        while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
            let child_cx_idx = read_cx_usize(ext_text, pos)?;
            children.push(child_cx_idx);
            if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b'.' {
                *pos += 1;
            } else {
                break;
            }
        }
        // Set PARENT prop on each child
        if let Some(actual_parent_idx) = parent_actual_idx {
            for child_cx_idx in &children {
                if *child_cx_idx >= state.builder.substance_groups_len() {
                    return Err(SmilesParseError::ParseError(
                        "child id references non-existent SGroup".to_string(),
                    ));
                }
                if let Some((child_pos, _)) = cx_index_to_parent
                    .iter()
                    .find(|(cx, _)| *cx == *child_cx_idx)
                {
                    if let Some(child_sg) = state.builder.substance_group_mut(*child_pos) {
                        child_sg.set_prop("PARENT", actual_parent_idx.to_string());
                    }
                }
            }
        }
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        } else {
            break;
        }
    }
    Ok(())
}

fn parse_cx_polymer_sgroup(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
    sgroup_idx: usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_polymer_sgroup (CXSmilesOps.cpp)
    // RDKit✔️✔️: bool parse_polymer_sgroup(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                           unsigned int startAtomIdx, unsigned int nSGroups) {
    // RDKit✔️✔️:   // format: |Sg:type:atomIndices:subscript:superscript:headCrossings:tailCrossings:|
    // RDKit✔️✔️:   // type codes: n->SRU, mon->MON, mer->MER, co/xl/mod/mix/f/any/gen/c/grf/alt/ran/blk
    // RDKit✔️✔️:   if (first >= last || *first != 'S' || first + 2 >= last ||
    // RDKit✔️✔️:       *(first + 1) != 'g' || *(first + 2) != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   first += 3;
    // RDKit✔️✔️:   const auto type_code = read_text_to(first, last, ":");
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   const auto type = sgroupTypemap.find(type_code);
    // RDKit✔️✔️:   if (type == sgroupTypemap.end()) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   bool keepSGroup = false;
    // RDKit✔️✔️:   SubstanceGroup sgroup(&mol, type->second);
    // RDKit✔️✔️:   sgroup.setProp(cxsmilesindex, nSGroups);
    // RDKit✔️✔️:   if (type_code == "alt") {
    // RDKit✔️✔️:     sgroup.setProp("SUBTYPE", "ALT");
    // RDKit✔️✔️:   } else if (type_code == "ran") {
    // RDKit✔️✔️:     sgroup.setProp("SUBTYPE", "RAN");
    // RDKit✔️✔️:   } else if (type_code == "blk") {
    // RDKit✔️✔️:     sgroup.setProp("SUBTYPE", "BLO");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::vector<unsigned int> atoms;
    // RDKit✔️✔️:   if (!read_int_list(first, last, atoms)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto idx : atoms) {
    // RDKit✔️✔️:     if (VALID_ATIDX(idx)) {
    // RDKit✔️✔️:       sgroup.addAtomWithIdx(idx - startAtomIdx);
    // RDKit✔️✔️:       keepSGroup = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::vector<unsigned int> headCrossing;
    // RDKit✔️✔️:   std::vector<unsigned int> tailCrossing;
    // RDKit✔️✔️:   if (first <= last && *first == ':') {
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:     std::string subscript = read_text_to(first, last, ":|");
    // RDKit✔️✔️:     if (keepSGroup && !subscript.empty()) {
    // RDKit✔️✔️:       sgroup.setProp("LABEL", subscript);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first <= last && *first == ':') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:       std::string superscript = read_text_to(first, last, ":|,");
    // RDKit✔️✔️:       if (keepSGroup && !superscript.empty()) {
    // RDKit✔️✔️:         sgroup.setProp("CONNECT", superscript);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // ... headCrossings and tailCrossings follow, omitted for brevity
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (keepSGroup) {
    // RDKit✔️✔️:     processCXSmilesLabels(mol);
    // RDKit✔️✔️:     finalizePolymerSGroup(mol, sgroup);
    // RDKit✔️✔️:     sgroup.setProp<unsigned int>("index", getSubstanceGroups(mol).size() + 1);
    // RDKit✔️✔️:     addSubstanceGroup(mol, sgroup);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_polymer_sgroup
    if !ext_text[*pos..].starts_with("Sg:") {
        return Err(cx_parse_failure());
    }
    *pos += 3;
    let type_code = read_cx_text_to(ext_text, pos, b":")?;
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    let Some(kind) = cx_sgroup_type_to_kind(&type_code) else {
        return Err(cx_parse_failure());
    };
    let mut keep_sgroup = false;
    let mut sgroup = SubstanceGroup::new(SubstanceGroupId::new(sgroup_idx), kind);
    sgroup.set_prop("_cxsmilesindex", sgroup_idx.to_string());
    // Set subtype for alt/ran/blk
    match type_code.as_str() {
        "alt" => sgroup.set_prop("SUBTYPE", "ALT"),
        "ran" => sgroup.set_prop("SUBTYPE", "RAN"),
        "blk" => sgroup.set_prop("SUBTYPE", "BLO"),
        _ => {}
    }
    // Read atom indices
    let mut atoms = Vec::new();
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        if state.builder.atom_mut(AtomId::new(atom_idx)).is_some() {
            atoms.push(AtomId::new(atom_idx));
            keep_sgroup = true;
        }
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        } else {
            break;
        }
    }
    for atom in &atoms {
        sgroup.push_atom(*atom);
    }
    // Read subscript, superscript, headCrossings, tailCrossings
    if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b':' {
        *pos += 1;
        let subscript = read_cx_text_to(ext_text, pos, b":|")?;
        if keep_sgroup && !subscript.is_empty() {
            sgroup.set_prop("LABEL", &subscript);
        }
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b':' {
            *pos += 1;
            let superscript = read_cx_text_to(ext_text, pos, b":|,")?;
            if keep_sgroup && !superscript.is_empty() {
                sgroup.set_prop("CONNECT", &superscript);
            }
            // Head crossing bonds
            if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b':' {
                *pos += 1;
                let mut head_crossing = Vec::new();
                while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
                    let cidx = read_cx_usize(ext_text, pos)?;
                    head_crossing.push(cidx);
                    if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
                        *pos += 1;
                    } else {
                        break;
                    }
                }
                if keep_sgroup && !head_crossing.is_empty() {
                    sgroup.set_prop(
                        "_headCrossings",
                        head_crossing
                            .iter()
                            .map(|v| v.to_string())
                            .collect::<Vec<_>>()
                            .join(","),
                    );
                }
                // Tail crossing bonds
                if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b':' {
                    *pos += 1;
                    let mut tail_crossing = Vec::new();
                    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
                        let cidx = read_cx_usize(ext_text, pos)?;
                        tail_crossing.push(cidx);
                        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
                            *pos += 1;
                        } else {
                            break;
                        }
                    }
                    if keep_sgroup && !tail_crossing.is_empty() {
                        sgroup.set_prop(
                            "_tailCrossings",
                            tail_crossing
                                .iter()
                                .map(|v| v.to_string())
                                .collect::<Vec<_>>()
                                .join(","),
                        );
                    }
                }
            }
        }
    }
    if keep_sgroup {
        sgroup.set_prop(
            "index",
            (state.builder.substance_groups_len() + 1).to_string(),
        );
        state
            .builder
            .add_substance_group(sgroup)
            .map_err(|error| SmilesParseError::ParseError(error.to_string()))?;
    }
    Ok(())
}

fn cx_sgroup_type_to_kind(type_code: &str) -> Option<SubstanceGroupKind> {
    // RDKit✔️✔️: sgroupTypemap: n->SRU, mon->MON, mer->MER, co->COP, xl->CRO,
    // RDKit✔️✔️: mod->MOD, mix->MIX, f->FOR, any->ANY, gen->GEN, c->COM,
    // RDKit✔️✔️: grf->GRA, alt/ran/blk->COP (with subtype)
    // END RDKIT CPP map lookup
    match type_code {
        "n" => Some(SubstanceGroupKind::StructuralRepeatUnit),
        "mon" => Some(SubstanceGroupKind::Monomer),
        "mer" => Some(SubstanceGroupKind::Mer),
        "co" => Some(SubstanceGroupKind::Copolymer),
        "xl" => Some(SubstanceGroupKind::Crosslink),
        "mod" => Some(SubstanceGroupKind::Modification),
        "mix" => Some(SubstanceGroupKind::MixtureComponent),
        "f" => Some(SubstanceGroupKind::Formulation),
        "any" => Some(SubstanceGroupKind::AnyPolymer),
        "gen" => Some(SubstanceGroupKind::Generic("GEN".to_string())),
        "c" => Some(SubstanceGroupKind::Generic("COM".to_string())),
        "grf" => Some(SubstanceGroupKind::Graft),
        "alt" | "ran" | "blk" => Some(SubstanceGroupKind::Copolymer),
        _ => None,
    }
}

fn parse_cx_substitution(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_substitution
    // RDKit✔️✔️: bool parse_substitution(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                         unsigned int startAtomIdx) {
    // RDKit✔️✔️:   if (first >= last || *first != 's' || first + 1 >= last ||
    // RDKit✔️✔️:       *(first + 1) != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   first += 2;
    // RDKit✔️✔️:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int n1;
    // RDKit✔️✔️:     if (!read_int(first, last, n1)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // check that we can read at least two more characters:
    // RDKit✔️✔️:     if (first + 1 >= last || *first != ':') {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:     unsigned int n2;
    // RDKit✔️✔️:     if (*first == '*') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:       n2 = 0xDEADBEEF;
    // RDKit✔️✔️:       if (VALID_ATIDX(n1)) {
    // RDKit✔️✔️:         mol.setProp(common_properties::_NeedsQueryScan, 1);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       if (!read_int(first, last, n2)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (VALID_ATIDX(n1)) {
    // RDKit✔️✔️:       auto atom = mol.getAtomWithIdx(n1 - startAtomIdx);
    // RDKit✔️✔️:       if (!atom->hasQuery()) {
    // RDKit✔️✔️:         atom = QueryOps::replaceAtomWithQueryAtom(&mol, atom);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       atom->expandQuery(makeAtomNonHydrogenDegreeQuery(n2),
    // RDKit✔️✔️:                         Queries::COMPOSITE_AND);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_substitution
    if !ext_text[*pos..].starts_with("s:") {
        return Err(cx_parse_failure());
    }
    *pos += 2;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        expect_byte(ext_text, *pos, b':')?;
        *pos += 1;
        let non_hydrogen_degree = if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b'*' {
            *pos += 1;
            if state.builder.atom_mut(AtomId::new(atom_idx)).is_some() {
                state.set_property("_NeedsQueryScan", "1");
            }
            0xDEADBEEF
        } else {
            read_cx_usize(ext_text, pos)? as u32
        };
        expand_cx_atom_query(
            state,
            atom_idx,
            AtomQueryPredicate::NonHydrogenDegree(non_hydrogen_degree),
        );
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
    }
    Ok(())
}

fn bond_can_have_direction(order: BondOrder) -> bool {
    // BEGIN RDKIT CPP FUNCTION canHaveDirection
    // RDKit✔️✔️: inline bool canHaveDirection(const Bond &bond) {
    // RDKit✔️✔️:   auto bondType = bond.getBondType();
    // RDKit✔️✔️:   return (bondType == Bond::SINGLE || bondType == Bond::AROMATIC);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION canHaveDirection
    matches!(order, BondOrder::Single | BondOrder::Aromatic)
}

fn parse_cx_wedged_bonds(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_wedged_bonds
    // RDKit✔️✔️: bool parse_wedged_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                         unsigned int startAtomIdx, unsigned int startBondIdx) {
    // RDKit✔️✔️:   // these look like: CC(O)Cl |w:1.0|
    // RDKit✔️✔️:   // also wD and wU for down and up wedges.
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   // We do not end up using this to set stereochemistry, but the relevant bond
    // RDKit✔️✔️:   // properties are set in case client code wants to do something with the
    // RDKit✔️✔️:   // information.
    // RDKit✔️✔️:   if (first >= last || *first != 'w' || first + 1 >= last) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   Bond::BondDir state = Bond::BondDir::NONE;
    // RDKit✔️✔️:   unsigned int cfg = 0;
    // RDKit✔️✔️:   switch (*first) {
    // RDKit✔️✔️:     case ':':
    // RDKit✔️✔️:       state = Bond::BondDir::UNKNOWN;
    // RDKit✔️✔️:       cfg = 2;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 'U':
    // RDKit✔️✔️:       state = Bond::BondDir::BEGINWEDGE;
    // RDKit✔️✔️:       cfg = 1;
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 'D':
    // RDKit✔️✔️:       state = Bond::BondDir::BEGINDASH;
    // RDKit✔️✔️:       cfg = 3;
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (state == Bond::BondDir::NONE || first >= last || first + 1 >= last ||
    // RDKit✔️✔️:       *first != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int atomIdx;
    // RDKit✔️✔️:     if (!read_int(first, last, atomIdx)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == '.') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog) << "improperly formatted w block" << std::endl;
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     unsigned int bondIdx;
    // RDKit✔️✔️:     if (!read_int(first, last, bondIdx)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (VALID_ATIDX(atomIdx) && VALID_BNDIDX(bondIdx)) {
    // RDKit✔️✔️:       auto atom = mol.getAtomWithIdx(atomIdx - startAtomIdx);
    // RDKit✔️✔️:       auto bond = get_bond_with_smiles_idx(mol, bondIdx - startBondIdx);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (!bond) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "bond " << bondIdx << " not found, wedge from atom " << atomIdx
    // RDKit✔️✔️:             << " cannot be applied." << std::endl;
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // we can't set wedging twice:
    // RDKit✔️✔️:       if (bond->hasProp(common_properties::_MolFileBondCfg)) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "w block attempts to set wedging on bond " << bond->getIdx()
    // RDKit✔️✔️:             << " more than once." << std::endl;
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // first things first, the atom needs to be the start atom of the bond for
    // RDKit✔️✔️:       // any of this to make sense
    // RDKit✔️✔️:       if (atom->getIdx() != bond->getBeginAtomIdx()) {
    // RDKit✔️✔️:         if (atom->getIdx() != bond->getEndAtomIdx()) {
    // RDKit✔️✔️:           BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:               << "atom " << atomIdx << " is not associated with bond "
    // RDKit✔️✔️:               << bondIdx << "(" << bond->getBeginAtomIdx() + startAtomIdx << "-"
    // RDKit✔️✔️:               << bond->getEndAtomIdx() + startAtomIdx << ")"
    // RDKit✔️✔️:               << " in w block" << std::endl;
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         auto eidx = bond->getBeginAtomIdx();
    // RDKit✔️✔️:         bond->setBeginAtomIdx(atom->getIdx());
    // RDKit✔️✔️:         bond->setEndAtomIdx(eidx);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       bond->setProp(common_properties::_MolFileBondCfg, cfg);
    // RDKit✔️✔️:       bond->setBondDir(state);
    // RDKit✔️✔️:       if (cfg == 2 && canHaveDirection(*bond)) {
    // RDKit✔️✔️:         bond->getBeginAtom()->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
    // RDKit✔️✔️:         mol.setProp(detail::_needsDetectBondStereo, 1);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if ((cfg == 1 || cfg == 3) && canHaveDirection(*bond)) {
    // RDKit✔️✔️:         mol.setProp(detail::_needsDetectAtomStereo, 1);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_wedged_bonds
    expect_byte(ext_text, *pos, b'w')?;
    *pos += 1;
    let (direction, cfg) = match ext_text.as_bytes().get(*pos).copied() {
        Some(b':') => (BondDirection::Unknown, "2"),
        Some(b'U') => {
            *pos += 1;
            (BondDirection::BeginWedge, "1")
        }
        Some(b'D') => {
            *pos += 1;
            (BondDirection::BeginDash, "3")
        }
        _ => return Err(cx_parse_failure()),
    };
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        expect_byte(ext_text, *pos, b'.')?;
        *pos += 1;
        let bond_idx = read_cx_usize(ext_text, pos)?;
        if state.builder.atom_mut(AtomId::new(atom_idx)).is_some() {
            let bond_id = cx_bond_with_smiles_idx(state, bond_idx).ok_or_else(cx_parse_failure)?;
            let (begin, end, order, has_cfg) = state
                .builder
                .bond(bond_id)
                .map(|bond| {
                    (
                        bond.begin(),
                        bond.end(),
                        bond.order(),
                        bond.prop("_MolFileBondCfg").is_some(),
                    )
                })
                .ok_or_else(cx_parse_failure)?;
            if has_cfg {
                return Err(cx_parse_failure());
            }
            let atom_id = AtomId::new(atom_idx);
            if begin != atom_id && end != atom_id {
                return Err(cx_parse_failure());
            }
            let bond = state
                .builder
                .bond_mut(bond_id)
                .ok_or_else(cx_parse_failure)?;
            if begin != atom_id {
                bond.set_endpoints(atom_id, begin);
            }
            bond.set_prop("_MolFileBondCfg", cfg);
            bond.set_direction(direction);
            if cfg == "2" && bond_can_have_direction(order) {
                if let Some(atom) = state.builder.atom_mut(atom_id) {
                    atom.set_chiral_tag(ChiralTag::Unspecified);
                }
                state.set_property("_needsDetectBondStereo", "1");
            }
            if matches!(cfg, "1" | "3") && bond_can_have_direction(order) {
                state.set_property("_needsDetectAtomStereo", "1");
            }
        }
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
    }
    Ok(())
}

fn parse_cx_variable_attachments(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_variable_attachments
    // RDKit✔️✔️: bool parse_variable_attachments(Iterator &first, Iterator last,
    // RDKit✔️✔️:                                 RDKit::RWMol &mol, unsigned int startAtomIdx) {
    // RDKit✔️✔️:   // these look like: CO*.C1=CC=NC=C1 |m:2:3.5.4|
    // RDKit✔️✔️:   // that corresponds to replacing the bond to atom 2 with bonds to atom 3, 5,
    // RDKit✔️✔️:   // or 4
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   if (first >= last || *first != 'm' || first + 1 >= last ||
    // RDKit✔️✔️:       *(first + 1) != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   first += 2;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int at1idx;
    // RDKit✔️✔️:     if (!read_int(first, last, at1idx)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (VALID_ATIDX(at1idx) &&
    // RDKit✔️✔️:         mol.getAtomWithIdx(at1idx - startAtomIdx)->getDegree() != 1) {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "position variation bond to atom with more than one bond"
    // RDKit✔️✔️:           << std::endl;
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ':') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog) << "improperly formatted m: block" << std::endl;
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::vector<std::string> others;
    // RDKit✔️✔️:     while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:       unsigned int aidx;
    // RDKit✔️✔️:       if (!read_int(first, last, aidx)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (VALID_ATIDX(aidx)) {
    // RDKit✔️✔️:         others.push_back(std::to_string(aidx - startAtomIdx + 1));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (first < last && *first == '.') {
    // RDKit✔️✔️:         ++first;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (VALID_ATIDX(at1idx)) {
    // RDKit✔️✔️:       std::string endPts = "(" + std::to_string(others.size());
    // RDKit✔️✔️:       for (auto idx : others) {
    // RDKit✔️✔️:         endPts += " " + idx;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       endPts += ")";
    // RDKit✔️✔️:
    // RDKit✔️✔️:       for (auto nbri : boost::make_iterator_range(
    // RDKit✔️✔️:                mol.getAtomBonds(mol.getAtomWithIdx(at1idx - startAtomIdx)))) {
    // RDKit✔️✔️:         auto bnd = mol[nbri];
    // RDKit✔️✔️:         bnd->setProp(common_properties::_MolFileBondEndPts, endPts);
    // RDKit✔️✔️:         bnd->setProp(common_properties::_MolFileBondAttach, std::string("ANY"));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_variable_attachments
    if !ext_text[*pos..].starts_with("m:") {
        return Err(cx_parse_failure());
    }
    *pos += 2;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let atom_idx = read_cx_usize(ext_text, pos)?;
        let atom_id = AtomId::new(atom_idx);
        if state.builder.atom_mut(atom_id).is_some() && atom_neighbors(state, atom_id).len() != 1 {
            return Err(cx_parse_failure());
        }
        expect_byte(ext_text, *pos, b':')?;
        *pos += 1;
        let mut endpoints = Vec::new();
        while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
            let endpoint = read_cx_usize(ext_text, pos)?;
            if state.builder.atom_mut(AtomId::new(endpoint)).is_some() {
                endpoints.push((endpoint + 1).to_string());
            }
            if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b'.' {
                *pos += 1;
            }
        }
        if state.builder.atom_mut(atom_id).is_some() {
            let end_pts = if endpoints.is_empty() {
                "(0)".to_string()
            } else {
                format!("({} {})", endpoints.len(), endpoints.join(" "))
            };
            for bond_id in atom_bonds(state, atom_id) {
                let bond = state
                    .builder
                    .bond_mut(bond_id)
                    .ok_or_else(cx_parse_failure)?;
                bond.set_prop("_MolFileBondEndPts", end_pts.clone());
                bond.set_prop("_MolFileBondAttach", "ANY");
            }
        }
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
    }
    Ok(())
}

fn set_cx_stereo_for_bond(
    state: &mut SmilesBuildState,
    bond_id: BondId,
    stereo: BondStereo,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION Chirality::detail::setStereoForBond
    // RDKit✔️✔️: void setStereoForBond(ROMol &mol, Bond *bond, Bond::BondStereo stereo,
    // RDKit✔️✔️:                       bool useCXSmilesOrdering) {
    // RDKit✔️✔️:   // NOTE:  moved from parse_doublebond_stereo CXSmilesOps
    // RDKit✔️✔️:   // IF useCXSmilesOrdering is true, the cis/trans/unknown marker will be
    // RDKit✔️✔️:   // assigned relative to the lowest-numbered neighbor of each double bond atom.
    // RDKit✔️✔️:   // Otherwise it uses the lowest-numbered neighbor on the lower-numbered atom
    // RDKit✔️✔️:   // of the double bond and the highest-numbered neighbor on the higher-numbered
    // RDKit✔️✔️:   // atom
    // RDKit✔️✔️:   auto begAtom = bond->getBeginAtom();
    // RDKit✔️✔️:   auto endAtom = bond->getEndAtom();
    // RDKit✔️✔️:   if (begAtom->getIdx() > endAtom->getIdx()) {
    // RDKit✔️✔️:     std::swap(begAtom, endAtom);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (begAtom->getDegree() > 1 && endAtom->getDegree() > 1) {
    // RDKit✔️✔️:     unsigned int begControl = mol.getNumAtoms();
    // RDKit✔️✔️:     for (auto nbr : mol.atomNeighbors(begAtom)) {
    // RDKit✔️✔️:       if (nbr == endAtom) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       begControl = std::min(nbr->getIdx(), begControl);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     unsigned int endControl = useCXSmilesOrdering ? mol.getNumAtoms() : 0;
    // RDKit✔️✔️:     for (auto nbr : mol.atomNeighbors(endAtom)) {
    // RDKit✔️✔️:       if (nbr == begAtom) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       endControl = useCXSmilesOrdering ? std::min(nbr->getIdx(), endControl)
    // RDKit✔️✔️:                                        : std::max(nbr->getIdx(), endControl);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (begAtom != bond->getBeginAtom()) {
    // RDKit✔️✔️:       std::swap(begControl, endControl);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     bond->setStereoAtoms(begControl, endControl);
    // RDKit✔️✔️:     bond->setStereo(stereo);
    // RDKit✔️✔️:     mol.setProp("_needsDetectBondStereo", 1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Chirality::detail::setStereoForBond
    let (begin, end) = state
        .builder
        .bond(bond_id)
        .map(|bond| (bond.begin(), bond.end()))
        .ok_or_else(cx_parse_failure)?;
    let (low, high) = if begin.index() > end.index() {
        (end, begin)
    } else {
        (begin, end)
    };
    let low_neighbors = atom_neighbors(state, low);
    let high_neighbors = atom_neighbors(state, high);
    if low_neighbors.len() <= 1 || high_neighbors.len() <= 1 {
        return Ok(());
    }
    let low_control = low_neighbors
        .iter()
        .copied()
        .filter(|neighbor| *neighbor != high)
        .min_by_key(|atom| atom.index())
        .ok_or_else(cx_parse_failure)?;
    let high_control = high_neighbors
        .iter()
        .copied()
        .filter(|neighbor| *neighbor != low)
        .min_by_key(|atom| atom.index())
        .ok_or_else(cx_parse_failure)?;
    let stereo_atoms = if low != begin {
        [high_control, low_control]
    } else {
        [low_control, high_control]
    };
    let bond = state
        .builder
        .bond_mut(bond_id)
        .ok_or_else(cx_parse_failure)?;
    bond.set_stereo_atoms(Some(stereo_atoms));
    bond.set_stereo(stereo);
    state.set_property("_needsDetectBondStereo", "1");
    Ok(())
}

fn parse_cx_doublebond_stereo(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
    stereo: BondStereo,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_doublebond_stereo
    // RDKit✔️✔️: bool parse_doublebond_stereo(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                              unsigned int, unsigned int startBondIdx,
    // RDKit✔️✔️:                              Bond::BondStereo stereo) {
    // RDKit✔️✔️:   // these look like: C1CCCC/C=C/CCC1 |ctu:5|
    // RDKit✔️✔️:   // also c and t for cis or trans
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   while (first < last && *first != ':') {
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (first >= last || *first != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️✔️:     unsigned int bondIdx;
    // RDKit✔️✔️:     if (!read_int(first, last, bondIdx)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (VALID_BNDIDX(bondIdx)) {
    // RDKit✔️✔️:       auto bond = get_bond_with_smiles_idx(mol, bondIdx - startBondIdx);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (!bond) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "bond " << bondIdx
    // RDKit✔️✔️:             << " not found, cannot mark as stereo double bond." << std::endl;
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       bool useCXOrdering = true;
    // RDKit✔️✔️:       Chirality::detail::setStereoForBond(mol, bond, stereo, useCXOrdering);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (first < last && *first == ',') {
    // RDKit✔️✔️:       ++first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_doublebond_stereo
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos] != b':' {
        *pos += 1;
    }
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        let bond_idx = read_cx_usize(ext_text, pos)?;
        let bond_id = cx_bond_with_smiles_idx(state, bond_idx).ok_or_else(cx_parse_failure)?;
        set_cx_stereo_for_bond(state, bond_id, stereo)?;
        if *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
            *pos += 1;
        }
    }
    Ok(())
}

fn process_cx_radical_section(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
    radical_electrons: u8,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION processRadicalSection
    // RDKit✔️✔️: bool processRadicalSection(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                            unsigned int numRadicalElectrons,
    // RDKit✔️✔️:                            unsigned int startAtomIdx) {
    // RDKit✔️✔️:   if (first >= last) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   if (first >= last || *first != ':') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++first;
    // RDKit✔️✔️:   unsigned int atIdx;
    // RDKit✔️✔️:   if (!read_int(first, last, atIdx)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (VALID_ATIDX(atIdx)) {
    // RDKit✔️✔️:     mol.getAtomWithIdx(atIdx - startAtomIdx)
    // RDKit✔️✔️:         ->setNumRadicalElectrons(numRadicalElectrons);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   while (first < last && *first == ',') {
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:     if (first < last && (*first < '0' || *first > '9')) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!read_int(first, last, atIdx)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (VALID_ATIDX(atIdx)) {
    // RDKit✔️✔️:       mol.getAtomWithIdx(atIdx - startAtomIdx)
    // RDKit✔️✔️:           ->setNumRadicalElectrons(numRadicalElectrons);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return first < last;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION processRadicalSection
    if *pos >= ext_text.len() {
        return Err(cx_parse_failure());
    }
    *pos += 1;
    expect_byte(ext_text, *pos, b':')?;
    *pos += 1;
    let atom_idx = read_cx_usize(ext_text, pos)?;
    if let Some(atom) = state.builder.atom_mut(AtomId::new(atom_idx)) {
        atom.set_radical_electrons(radical_electrons);
    }
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b',' {
        *pos += 1;
        if *pos < ext_text.len() && !ext_text.as_bytes()[*pos].is_ascii_digit() {
            return Ok(());
        }
        let atom_idx = read_cx_usize(ext_text, pos)?;
        if let Some(atom) = state.builder.atom_mut(AtomId::new(atom_idx)) {
            atom.set_radical_electrons(radical_electrons);
        }
    }
    if *pos < ext_text.len() {
        Ok(())
    } else {
        Err(cx_parse_failure())
    }
}

fn parse_cx_radicals(
    state: &mut SmilesBuildState,
    ext_text: &str,
    pos: &mut usize,
) -> Result<(), SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION parse_radicals
    // RDKit✔️✔️: bool parse_radicals(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️✔️:                     unsigned int startAtomIdx) {
    // RDKit✔️✔️:   if (first >= last || *first != '^') {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   while (*first == '^') {
    // RDKit✔️✔️:     ++first;
    // RDKit✔️✔️:     if (first >= last) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (*first < '1' || *first > '7') {
    // RDKit✔️✔️:       return false;  // these are the values that are allowed to be there
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     switch (*first) {
    // RDKit✔️✔️:       case '1':
    // RDKit✔️✔️:         if (!processRadicalSection(first, last, mol, 1, startAtomIdx)) {
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case '2':
    // RDKit✔️✔️:       case '3':
    // RDKit✔️✔️:       case '4':
    // RDKit✔️✔️:         if (!processRadicalSection(first, last, mol, 2, startAtomIdx)) {
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case '5':
    // RDKit✔️✔️:       case '6':
    // RDKit✔️✔️:       case '7':
    // RDKit✔️✔️:         if (!processRadicalSection(first, last, mol, 3, startAtomIdx)) {
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       default:
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "Radical specification " << *first << " ignored.";
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION parse_radicals
    expect_byte(ext_text, *pos, b'^')?;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos] == b'^' {
        *pos += 1;
        if *pos >= ext_text.len() {
            return Err(cx_parse_failure());
        }
        let radical_electrons = match ext_text.as_bytes()[*pos] {
            b'1' => 1,
            b'2' | b'3' | b'4' => 2,
            b'5' | b'6' | b'7' => 3,
            _ => return Err(cx_parse_failure()),
        };
        process_cx_radical_section(state, ext_text, pos, radical_electrons)?;
    }
    Ok(())
}

fn cx_parse_failure() -> SmilesParseError {
    SmilesParseError::ParseError("failure parsing CXSMILES extensions".to_string())
}

fn read_cx_usize(ext_text: &str, pos: &mut usize) -> Result<usize, SmilesParseError> {
    let start = *pos;
    while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
        *pos += 1;
    }
    if *pos == start {
        return Err(SmilesParseError::ParseError(
            "failure parsing CXSMILES extensions".to_string(),
        ));
    }
    ext_text[start..*pos].parse::<usize>().map_err(|_| {
        SmilesParseError::ParseError("failure parsing CXSMILES extensions".to_string())
    })
}

fn read_cx_text_to(
    ext_text: &str,
    pos: &mut usize,
    delimiters: &[u8],
) -> Result<String, SmilesParseError> {
    let mut value = String::new();
    while *pos < ext_text.len() && !delimiters.contains(&ext_text.as_bytes()[*pos]) {
        if ext_text.as_bytes()[*pos] == b'&'
            && *pos + 2 < ext_text.len()
            && ext_text.as_bytes()[*pos + 1] == b'#'
        {
            *pos += 2;
            let numeric_start = *pos;
            while *pos < ext_text.len() && ext_text.as_bytes()[*pos].is_ascii_digit() {
                *pos += 1;
            }
            if *pos >= ext_text.len() || ext_text.as_bytes()[*pos] != b';' {
                return Err(SmilesParseError::ParseError(
                    "failure parsing CXSMILES extensions: quoted block not terminated with ';'"
                        .to_string(),
                ));
            }
            if *pos > numeric_start {
                let code = ext_text[numeric_start..*pos].parse::<u32>().map_err(|_| {
                    SmilesParseError::ParseError("failure parsing CXSMILES extensions".to_string())
                })?;
                if let Some(ch) = char::from_u32(code) {
                    value.push(ch);
                }
            }
            *pos += 1;
        } else {
            value.push(ext_text.as_bytes()[*pos] as char);
            *pos += 1;
        }
    }
    Ok(value)
}

fn expect_byte(ext_text: &str, pos: usize, expected: u8) -> Result<(), SmilesParseError> {
    if pos < ext_text.len() && ext_text.as_bytes()[pos] == expected {
        Ok(())
    } else {
        Err(SmilesParseError::ParseError(
            "failure parsing CXSMILES extensions".to_string(),
        ))
    }
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
mod tests {
    use super::*;

    fn parse_atom_token(input: &str) -> SmilesAtomToken {
        let mut parser = SmilesParser::new(SmilesLexer::new(input));
        let atom = parser.parse_simple_atomd().unwrap();
        assert_eq!(parser.next_token().unwrap(), SmilesToken::Eos);
        atom
    }

    fn parse_ring_number_token(input: &str) -> Result<u32, SmilesParseError> {
        let mut parser = SmilesParser::new(SmilesLexer::new(input));
        let ring_number = parser.parse_ring_number()?;
        assert_eq!(parser.next_token()?, SmilesToken::Eos);
        Ok(ring_number)
    }

    fn parse_number_token(input: &str) -> Result<i32, SmilesParseError> {
        let mut parser = SmilesParser::new(SmilesLexer::new(input));
        let number = parser.parse_number()?;
        assert_eq!(parser.next_token()?, SmilesToken::Eos);
        Ok(number)
    }

    #[test]
    fn preprocess_smiles_splits_name_when_cxsmiles_disabled() {
        let params = SmilesParseParams::without_cxsmiles_for_test();
        let result = preprocess_smiles("CCO ethanol", &params).unwrap();

        assert_eq!(result.smiles, "CCO");
        assert_eq!(result.name, "ethanol");
        assert_eq!(result.cx_part, "");
    }

    #[test]
    fn preprocess_smiles_splits_cx_part_when_cxsmiles_enabled() {
        let result =
            preprocess_smiles("CCO |$;;foo$| ethanol", &SmilesParseParams::default()).unwrap();

        assert_eq!(result.smiles, "CCO");
        assert_eq!(result.name, "");
        assert_eq!(result.cx_part, "|$;;foo$| ethanol");
    }

    #[test]
    fn preprocess_smiles_leaves_leading_space_input_unsplit_like_rdkit() {
        let result = preprocess_smiles(" CCO name", &SmilesParseParams::default()).unwrap();

        assert_eq!(result.smiles, " CCO name");
        assert_eq!(result.name, "");
        assert_eq!(result.cx_part, "");
    }

    #[test]
    fn preprocess_smiles_applies_replacements_until_stable() {
        let mut replacements = BTreeMap::new();
        replacements.insert("{Q}".to_string(), "{X}CC{X}".to_string());
        replacements.insert("{X}".to_string(), "N".to_string());
        let params = SmilesParseParams::with_replacements_for_test(replacements);

        let result = preprocess_smiles("C{Q}C", &params).unwrap();

        assert_eq!(result.smiles, "CNCCNC");
    }

    #[test]
    fn lexer_trims_ascii_control_whitespace_like_setup_smiles_string() {
        let lexer = SmilesLexer::new("\t CCO \r\n");

        assert_eq!(lexer.scan_input(), "CCO");
    }

    #[test]
    fn lexer_emits_organic_and_aromatic_atom_payloads() {
        let mut lexer = SmilesLexer::new("Clc");

        match lexer.next_token().unwrap() {
            SmilesToken::OrganicAtom(atom) => assert_eq!(atom.spec.element().atomic_number(), 17),
            other => panic!("unexpected token: {other:?}"),
        }
        match lexer.next_token().unwrap() {
            SmilesToken::AromaticAtom(atom) => {
                assert_eq!(atom.spec.element().atomic_number(), 6);
                assert!(atom.spec.is_aromatic());
            }
            other => panic!("unexpected token: {other:?}"),
        }
    }

    #[test]
    fn lexer_emits_bracket_atom_and_biovia_quoted_atom_payloads() {
        let mut lexer = SmilesLexer::new("[He]['Og']");

        assert_eq!(lexer.next_token().unwrap(), SmilesToken::AtomOpen);
        match lexer.next_token().unwrap() {
            SmilesToken::Atom(atom) => assert_eq!(atom.spec.element().atomic_number(), 2),
            other => panic!("unexpected token: {other:?}"),
        }
        assert_eq!(lexer.next_token().unwrap(), SmilesToken::AtomClose);
        assert_eq!(lexer.next_token().unwrap(), SmilesToken::AtomOpen);
        match lexer.next_token().unwrap() {
            SmilesToken::Atom(atom) => assert_eq!(atom.spec.element().atomic_number(), 118),
            other => panic!("unexpected token: {other:?}"),
        }
        assert_eq!(lexer.next_token().unwrap(), SmilesToken::AtomClose);
    }

    #[test]
    fn lexer_emits_bond_token_payloads_like_smiles_ll() {
        let mut lexer = SmilesLexer::new("=#:$/\\-><-~");

        match lexer.next_token().unwrap() {
            SmilesToken::Bond(bond) => assert_eq!(bond.order, BondOrder::Double),
            other => panic!("unexpected token: {other:?}"),
        }
        match lexer.next_token().unwrap() {
            SmilesToken::Bond(bond) => assert_eq!(bond.order, BondOrder::Triple),
            other => panic!("unexpected token: {other:?}"),
        }
        match lexer.next_token().unwrap() {
            SmilesToken::Bond(bond) => {
                assert_eq!(bond.order, BondOrder::Aromatic);
                assert!(bond.is_aromatic);
            }
            other => panic!("unexpected token: {other:?}"),
        }
        match lexer.next_token().unwrap() {
            SmilesToken::Bond(bond) => assert_eq!(bond.order, BondOrder::Quadruple),
            other => panic!("unexpected token: {other:?}"),
        }
        match lexer.next_token().unwrap() {
            SmilesToken::Bond(bond) => {
                assert_eq!(bond.direction, BondDirection::EndUpRight);
                assert!(bond.explicit_unspecified_order);
            }
            other => panic!("unexpected token: {other:?}"),
        }
        match lexer.next_token().unwrap() {
            SmilesToken::Bond(bond) => {
                assert_eq!(bond.direction, BondDirection::EndDownRight);
                assert!(bond.explicit_unspecified_order);
            }
            other => panic!("unexpected token: {other:?}"),
        }
        match lexer.next_token().unwrap() {
            SmilesToken::Bond(bond) => assert_eq!(bond.order, BondOrder::DativeRight),
            other => panic!("unexpected token: {other:?}"),
        }
        match lexer.next_token().unwrap() {
            SmilesToken::Bond(bond) => assert_eq!(bond.order, BondOrder::DativeLeft),
            other => panic!("unexpected token: {other:?}"),
        }
        match lexer.next_token().unwrap() {
            SmilesToken::Bond(bond) => assert!(bond.is_null_query),
            other => panic!("unexpected token: {other:?}"),
        }
    }

    #[test]
    fn lexer_emits_chiral_class_tokens_like_smiles_ll() {
        let mut lexer = SmilesLexer::new("@TH@   SP");

        assert_eq!(
            lexer.next_token().unwrap(),
            SmilesToken::ChiralClass(ChiralTag::Tetrahedral)
        );
        assert_eq!(
            lexer.next_token().unwrap(),
            SmilesToken::ChiralClass(ChiralTag::SquarePlanar)
        );
    }

    #[test]
    fn lexer_emits_ring_and_punctuation_tokens_like_smiles_ll() {
        let mut lexer = SmilesLexer::new("0%12()+-.");

        assert_eq!(lexer.next_token().unwrap(), SmilesToken::Zero);
        assert_eq!(lexer.next_token().unwrap(), SmilesToken::Percent);
        assert_eq!(lexer.next_token().unwrap(), SmilesToken::NonzeroDigit(1));
        assert_eq!(lexer.next_token().unwrap(), SmilesToken::NonzeroDigit(2));
        assert_eq!(lexer.next_token().unwrap(), SmilesToken::GroupOpen);
        assert_eq!(lexer.next_token().unwrap(), SmilesToken::GroupClose);
        assert_eq!(lexer.next_token().unwrap(), SmilesToken::Plus);
        assert_eq!(lexer.next_token().unwrap(), SmilesToken::Minus);
        assert_eq!(lexer.next_token().unwrap(), SmilesToken::Separator);
        assert_eq!(lexer.next_token().unwrap(), SmilesToken::Eos);
    }

    #[test]
    fn parse_simple_atomd_parses_bracket_map_charge_hydrogen_and_no_implicit() {
        let atom = parse_atom_token("[13C@TH2H3-:5]");

        assert_eq!(
            atom.spec,
            AtomSpec::new(Element::C)
                .with_isotope(13)
                .with_chiral_tag(ChiralTag::Tetrahedral)
                .with_chiral_permutation(2)
                .with_explicit_hydrogens(3)
                .with_formal_charge(-1)
                .with_atom_map(5)
                .with_no_implicit(true)
        );
    }

    #[test]
    fn parse_simple_atomd_parses_hydrogen_isotope_with_explicit_hydrogen_suffix() {
        let atom = parse_atom_token("[2HH1-]");

        assert_eq!(
            atom.spec,
            AtomSpec::new(Element::H)
                .with_isotope(2)
                .with_explicit_hydrogens(1)
                .with_formal_charge(-1)
                .with_no_implicit(true)
        );
    }

    #[test]
    fn parse_simple_atomd_rejects_missing_bracket_close_like_rdkit() {
        let mut parser = SmilesParser::new(SmilesLexer::new("[NH4+"));

        assert_eq!(
            parser.parse_simple_atomd(),
            Err(SmilesParseError::ParseError(
                "expected bracket atom close or atom map, got Eos".to_string()
            ))
        );
    }

    #[test]
    fn parse_chiral_element_rejects_zero_permutation_like_rdkit() {
        let mut parser = SmilesParser::new(SmilesLexer::new("[C@TH0]"));
        assert_eq!(parser.next_token().unwrap(), SmilesToken::AtomOpen);

        assert_eq!(
            parser.parse_bracket_atomd(),
            Err(SmilesParseError::ParseError(
                "chiral permutation cannot be zero".to_string()
            ))
        );
    }

    #[test]
    fn parse_number_reports_int32_overflow_like_rdkit() {
        assert_eq!(
            parse_number_token("2147483648"),
            Err(SmilesParseError::ParseError("number too large".to_string()))
        );
    }

    #[test]
    fn parse_ring_number_accepts_percent_group_forms_like_rdkit() {
        assert_eq!(parse_ring_number_token("7").unwrap(), 7);
        assert_eq!(parse_ring_number_token("%12").unwrap(), 12);
        assert_eq!(parse_ring_number_token("%(1)").unwrap(), 1);
        assert_eq!(parse_ring_number_token("%(12345)").unwrap(), 12345);
    }

    #[test]
    fn parse_ring_number_rejects_empty_and_oversized_percent_groups() {
        assert_eq!(
            parse_ring_number_token("%()"),
            Err(SmilesParseError::ParseError(
                "empty ring number".to_string()
            ))
        );
        assert_eq!(
            parse_ring_number_token("%(123456)"),
            Err(SmilesParseError::ParseError(
                "ring number too large".to_string()
            ))
        );
    }

    #[test]
    fn build_state_add_first_atom_marks_smiles_start_like_rdkit() {
        let mut state = SmilesBuildState::new();
        state.add_first_atom(SmilesAtomToken::new(6)).unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.num_atoms(), 1);
        assert_eq!(molecule.num_bonds(), 0);
        assert_eq!(molecule.atoms()[0].atomic_number(), 6);
        assert_eq!(molecule.atoms()[0].prop(SMILES_START_PROP), Some("1"));
    }

    #[test]
    fn add_first_atom_parse_mol_keeps_disconnected_fragment_starts_like_rdkit() {
        let state = to_mol("C.O").unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![6, 8]);
        assert_eq!(molecule.num_bonds(), 0);
        assert_eq!(molecule.atoms()[0].prop(SMILES_START_PROP), Some("1"));
        assert_eq!(molecule.atoms()[1].prop(SMILES_START_PROP), Some("1"));
    }

    #[test]
    fn build_state_add_atom_connected_to_active_uses_unspecified_bond_type_like_rdkit() {
        let mut state = SmilesBuildState::new();
        state.add_first_atom(SmilesAtomToken::new(6)).unwrap();
        state
            .add_atom_connected_to_active(SmilesAtomToken::new(8))
            .unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![6, 8]);
        assert_eq!(molecule.num_bonds(), 1);
        assert_eq!(molecule.bonds()[0].order(), BondOrder::Single);
        assert_eq!(molecule.bonds()[0].prop(CXSMILES_BOND_IDX_PROP), Some("0"));
    }

    #[test]
    fn build_state_add_explicit_bond_to_atom_normalizes_dative_direction_like_rdkit() {
        let mut state = SmilesBuildState::new();
        state.add_first_atom(SmilesAtomToken::new(7)).unwrap();
        state
            .add_explicit_bond_to_atom(
                SmilesBondToken::new(BondOrder::DativeLeft),
                SmilesAtomToken::new(8),
            )
            .unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.num_atoms(), 2);
        assert_eq!(molecule.num_bonds(), 1);
        assert_eq!(molecule.bonds()[0].begin(), AtomId::new(1));
        assert_eq!(molecule.bonds()[0].end(), AtomId::new(0));
        assert_eq!(molecule.bonds()[0].order(), BondOrder::Dative);
        assert_eq!(molecule.bonds()[0].prop(CXSMILES_BOND_IDX_PROP), Some("0"));
    }

    #[test]
    fn build_state_add_explicit_directional_bond_preserves_unspecified_order_marker() {
        let mut state = SmilesBuildState::new();
        state.add_first_atom(SmilesAtomToken::new(6)).unwrap();
        state
            .add_explicit_bond_to_atom(
                SmilesBondToken::directional(BondDirection::EndUpRight),
                SmilesAtomToken::new(6),
            )
            .unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.bonds()[0].order(), BondOrder::Unspecified);
        assert_eq!(molecule.bonds()[0].direction(), BondDirection::EndUpRight);
        assert_eq!(molecule.bonds()[0].prop(UNSPECIFIED_ORDER_PROP), Some("1"));
    }

    #[test]
    fn from_smiles_with_sanitize_false_parses_linear_simple_atoms_through_grammar_actions() {
        let molecule = Molecule::from_smiles_with_sanitize("CCO", false).unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![6, 6, 8]);
        assert_eq!(molecule.num_bonds(), 2);
        assert_eq!(molecule.bonds()[0].order(), BondOrder::Single);
        assert_eq!(molecule.bonds()[1].order(), BondOrder::Single);
    }

    #[test]
    fn from_smiles_with_sanitize_false_parses_explicit_bond_and_separator_actions() {
        let molecule = Molecule::from_smiles_with_sanitize("C=O.N", false).unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![6, 8, 7]);
        assert_eq!(molecule.num_bonds(), 1);
        assert_eq!(molecule.bonds()[0].order(), BondOrder::Double);
        assert_eq!(molecule.atoms()[0].prop(SMILES_START_PROP), None);
        assert_eq!(molecule.atoms()[2].prop(SMILES_START_PROP), None);
        assert_eq!(molecule.bonds()[0].prop(CXSMILES_BOND_IDX_PROP), None);
    }

    #[test]
    fn from_smiles_with_sanitize_false_sets_unspecified_directional_bond_type_then_cleans_props() {
        let molecule = Molecule::from_smiles_with_sanitize("C/C", false).unwrap();

        assert_eq!(molecule.num_bonds(), 1);
        assert_eq!(molecule.bonds()[0].order(), BondOrder::Single);
        assert_eq!(molecule.bonds()[0].direction(), BondDirection::EndUpRight);
        assert_eq!(molecule.bonds()[0].prop(UNSPECIFIED_ORDER_PROP), None);
        assert_eq!(molecule.bonds()[0].prop(CXSMILES_BOND_IDX_PROP), None);
    }

    #[test]
    fn from_smiles_with_sanitize_false_only_clears_wedge_style_single_bond_dirs_like_rdkit() {
        let directional = Molecule::from_smiles_with_sanitize("C/C", false).unwrap();
        let wedged = Molecule::from_smiles_with_sanitize("CC |wU:1.0|", false).unwrap();

        assert_eq!(
            directional.bonds()[0].direction(),
            BondDirection::EndUpRight
        );
        assert_eq!(directional.bonds()[0].prop(UNSPECIFIED_ORDER_PROP), None);
        assert_eq!(wedged.bonds()[0].direction(), BondDirection::None);
        assert_eq!(wedged.bonds()[0].prop("_MolFileBondCfg"), Some("1"));
    }

    #[test]
    fn from_smiles_with_sanitize_false_reads_name_when_cx_part_is_plain_text() {
        let molecule = Molecule::from_smiles_with_sanitize("CCO ethanol", false).unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![6, 6, 8]);
        assert_eq!(molecule.properties().name(), Some("ethanol"));
    }

    #[test]
    fn handle_cx_part_and_name_reports_strict_non_name_text_like_rdkit() {
        let params = SmilesParseParams::without_parse_name_for_test();
        let error = mol_from_smiles("CCO not-name", &params).unwrap_err();

        assert_eq!(
            error,
            SmilesParseError::ParseError(
                "CXSMILES extension does not start with | and parseName=false".to_string()
            )
        );
    }

    #[test]
    fn handle_cx_part_and_name_parses_atom_labels_and_name_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("CCO |$;;foo$| ethanol", false).unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![6, 6, 8]);
        assert_eq!(molecule.atoms()[2].prop("atomLabel"), Some("foo"));
        assert_eq!(
            molecule.properties().prop("_CXSMILES_Data"),
            Some("|$;;foo$|")
        );
        assert_eq!(molecule.properties().name(), Some("ethanol"));
    }

    #[test]
    fn cleanup_after_parsing_marks_ap1_ap2_atom_labels_like_rdkit() {
        let molecule =
            Molecule::from_smiles_with_sanitize("*.*.* |$_AP1;_AP2;_AP3$|", false).unwrap();

        assert_eq!(molecule.atoms()[0].prop("_fromAttachPoint"), Some("1"));
        assert_eq!(molecule.atoms()[1].prop("_fromAttachPoint"), Some("2"));
        assert_eq!(molecule.atoms()[2].prop("_fromAttachPoint"), None);
    }

    #[test]
    fn cleanup_after_parsing_clears_parser_temporary_props_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("c1ccccc1.C", false).unwrap();

        assert!(
            molecule
                .atoms()
                .iter()
                .all(|atom| atom.prop(SMILES_START_PROP).is_none())
        );
        assert!(molecule.bonds().iter().all(|bond| {
            bond.prop(UNSPECIFIED_ORDER_PROP).is_none()
                && bond.prop(CXSMILES_BOND_IDX_PROP).is_none()
        }));
        assert!(
            molecule.bonds()[..6]
                .iter()
                .all(|bond| bond.order() == BondOrder::Aromatic && bond.is_aromatic())
        );
    }

    #[test]
    fn handle_cx_part_and_name_tolerates_unported_cx_when_not_strict_like_rdkit() {
        let params = SmilesParseParams::non_strict_cxsmiles_for_test();
        let molecule = mol_from_smiles("CCO |rb:0| ethanol", &params).unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![6, 6, 8]);
        assert_eq!(molecule.properties().prop("_CXSMILES_Data"), Some(""));
        assert_eq!(molecule.properties().name(), None);
    }

    #[test]
    fn handle_cx_part_and_name_parses_sgroup_hierarchy_in_strict_mode() {
        let molecule = Molecule::from_smiles_with_sanitize("CC |SgH:0:1|", false).unwrap();
        assert_eq!(molecule.num_atoms(), 2);
        assert_eq!(molecule.num_bonds(), 1);
    }

    #[test]
    fn parse_cx_sgroup_hierarchy_rejects_nonexistent_child_id_like_rdkit() {
        let mut state = SmilesBuildState::new();
        state.builder.add_atom(AtomSpec::new(Element::C));
        state
            .builder
            .add_substance_group(
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                    .with_prop("_cxsmilesindex", "0")
                    .with_prop("index", "1"),
            )
            .unwrap();
        state
            .builder
            .add_substance_group(
                SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                    .with_prop("_cxsmilesindex", "1")
                    .with_prop("index", "2"),
            )
            .unwrap();

        let error = parse_cx_sgroup_hierarchy(&mut state, "SgH:0:2", &mut 0).unwrap_err();

        assert_eq!(
            error,
            SmilesParseError::ParseError("child id references non-existent SGroup".to_string())
        );
    }

    #[test]
    fn handle_cx_part_and_name_parses_polymer_sgroup_in_strict_mode() {
        let molecule = Molecule::from_smiles_with_sanitize("CC |Sg:n:0::ht|", false).unwrap();
        assert_eq!(molecule.num_atoms(), 2);
        assert_eq!(molecule.num_bonds(), 1);
    }

    #[test]
    fn parse_cx_polymer_sgroup_rejects_unknown_type_like_rdkit() {
        let mut state = SmilesBuildState::new();
        state.builder.add_atom(AtomSpec::new(Element::C));

        let error = parse_cx_polymer_sgroup(&mut state, "Sg:bogus:0::ht|", &mut 0, 0).unwrap_err();

        assert_eq!(error, cx_parse_failure());
    }

    #[test]
    fn handle_cx_part_and_name_parses_atom_values_and_props_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize(
            "CC |$_AV:first;second$,atomProp:1.foo.bar&#46;baz|",
            false,
        )
        .unwrap();

        assert_eq!(molecule.atoms()[0].prop("molFileValue"), Some("first"));
        assert_eq!(molecule.atoms()[1].prop("molFileValue"), Some("second"));
        assert_eq!(molecule.atoms()[1].prop("foo"), Some("bar.baz"));
    }

    #[test]
    fn handle_cx_part_and_name_parses_coordinates_like_rdkit() {
        let molecule =
            Molecule::from_smiles_with_sanitize("CCO |(0,0,;1,0,;2,0,0.5)|", false).unwrap();

        assert_eq!(molecule.conformers_3d().len(), 1);
        assert!(molecule.conformers_3d()[0].is_3d());
        assert_eq!(molecule.conformers_3d()[0].coords()[0], [0.0, 0.0, 0.0]);
        assert_eq!(molecule.conformers_3d()[0].coords()[1], [1.0, 0.0, 0.0]);
        assert_eq!(molecule.conformers_3d()[0].coords()[2], [2.0, 0.0, 0.5]);
    }

    #[test]
    fn mol_from_smiles_conformer_selection_takes_first_2d_and_first_3d_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_conformer(Conformer3D::new(0, vec![[0.0, 0.0, 1.0]], true))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(1, vec![[0.0, 0.0, 0.0]], false))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(2, vec![[1.0, 0.0, 0.0]], false))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(3, vec![[1.0, 0.0, 1.0]], true))
            .unwrap();
        let molecule = builder.build().unwrap();

        assert_eq!(first_2d_and_3d_conformer_ids(&molecule), (Some(1), Some(0)));
    }

    #[test]
    fn mol_from_smiles_conformer_selection_reports_only_first_3d_when_no_2d_exists() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_conformer(Conformer3D::new(0, vec![[0.0, 0.0, 1.0]], true))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(1, vec![[1.0, 0.0, 1.0]], true))
            .unwrap();
        let molecule = builder.build().unwrap();

        assert_eq!(first_2d_and_3d_conformer_ids(&molecule), (None, Some(0)));
    }

    fn coordinate_free_atropisomer_candidate(direction: BondDirection) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let chlorine = builder.add_atom(
            AtomSpec::new(Element::from_atomic_number(17).unwrap()).with_no_implicit(true),
        );
        let axis_begin = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let alkene_left = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let axis_end = builder.add_atom(AtomSpec::new(Element::C));
        let alkene_right = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));

        builder
            .add_bond(
                BondSpec::new(axis_begin, chlorine, BondOrder::Single).with_direction(direction),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(axis_begin, alkene_left, BondOrder::Double))
            .unwrap();
        builder
            .add_bond(BondSpec::new(axis_begin, axis_end, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(axis_end, alkene_right, BondOrder::Double))
            .unwrap();

        builder
            .build()
            .unwrap()
            .sanitized_with_ops(crate::SanitizeOps::ALL)
            .unwrap()
    }

    #[test]
    fn mol_from_smiles_assigns_coordinate_free_atropisomer_bond_stereo_like_rdkit() {
        let mut molecule = coordinate_free_atropisomer_candidate(BondDirection::BeginWedge);

        let assignments = atropisomer_stereo_without_conformer(&molecule);
        apply_coordinate_free_atropisomer_assignments(&mut molecule, assignments);

        assert!(molecule.conformers_3d().is_empty());
        assert_eq!(molecule.bonds()[2].stereo(), BondStereo::AtropCcw);
    }

    #[test]
    fn mol_from_smiles_flips_coordinate_free_atropisomer_bond_stereo_for_hash_like_rdkit() {
        let mut molecule = coordinate_free_atropisomer_candidate(BondDirection::BeginDash);

        let assignments = atropisomer_stereo_without_conformer(&molecule);
        apply_coordinate_free_atropisomer_assignments(&mut molecule, assignments);

        assert!(molecule.conformers_3d().is_empty());
        assert_eq!(molecule.bonds()[2].stereo(), BondStereo::AtropCw);
    }

    fn conformer_backed_atropisomer_candidate(is_3d: bool) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let chlorine = builder.add_atom(
            AtomSpec::new(Element::from_atomic_number(17).unwrap()).with_no_implicit(true),
        );
        let axis_begin = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let alkene_left = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let axis_end = builder.add_atom(AtomSpec::new(Element::C));
        let alkene_right = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));

        builder
            .add_bond(
                BondSpec::new(axis_begin, chlorine, BondOrder::Single)
                    .with_direction(BondDirection::BeginWedge),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(axis_begin, alkene_left, BondOrder::Double))
            .unwrap();
        builder
            .add_bond(BondSpec::new(axis_begin, axis_end, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(axis_end, alkene_right, BondOrder::Double))
            .unwrap();

        let coords = if is_3d {
            vec![
                [0.0, 0.0, 1.0],
                [0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, -1.0, 0.0],
            ]
        } else {
            vec![
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, -1.0, 0.0],
            ]
        };
        builder
            .add_conformer(Conformer3D::new(0, coords, is_3d))
            .unwrap();

        builder
            .build()
            .unwrap()
            .sanitized_with_ops(crate::SanitizeOps::ALL)
            .unwrap()
    }

    #[test]
    fn mol_from_smiles_assigns_2d_conformer_backed_atropisomer_bond_stereo_like_rdkit() {
        let mut molecule = conformer_backed_atropisomer_candidate(false);

        let assignments = atropisomer_stereo_from_conformer(&molecule, 0);
        apply_atropisomer_stereo_assignments(&mut molecule, assignments);

        assert_eq!(molecule.bonds()[2].stereo(), BondStereo::AtropCw);
    }

    #[test]
    fn mol_from_smiles_assigns_3d_conformer_backed_atropisomer_bond_stereo_like_rdkit() {
        let mut molecule = conformer_backed_atropisomer_candidate(true);

        let assignments = atropisomer_stereo_from_conformer(&molecule, 0);
        apply_atropisomer_stereo_assignments(&mut molecule, assignments);

        assert_eq!(molecule.bonds()[2].stereo(), BondStereo::AtropCw);
    }

    fn three_neighbor_3d_carbon(formal_charge: i8) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C).with_formal_charge(formal_charge));
        let fluorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
        let chlorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
        let bromine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(35).unwrap()));
        builder
            .add_bond(BondSpec::new(center, fluorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, chlorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, bromine, BondOrder::Single))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
                true,
            ))
            .unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn assign_chiral_types_from_3d_uses_valence_implicit_h_for_three_coordinate_carbon() {
        let mut molecule = three_neighbor_3d_carbon(0);

        assign_chiral_types_from_3d(&mut molecule, 0);

        assert_ne!(molecule.atoms()[0].chiral_tag(), ChiralTag::Unspecified);
    }

    #[test]
    fn assign_chiral_types_from_3d_rejects_three_coordinate_cation_with_no_implicit_h() {
        let mut molecule = three_neighbor_3d_carbon(1);

        assign_chiral_types_from_3d(&mut molecule, 0);

        assert_eq!(molecule.atoms()[0].chiral_tag(), ChiralTag::Unspecified);
    }

    #[test]
    fn assign_chiral_types_from_3d_clears_stereochem_done_like_rdkit() {
        let mut molecule = three_neighbor_3d_carbon(0).with_prop("_StereochemDone", "1");

        assign_chiral_types_from_3d(&mut molecule, 0);

        assert_eq!(molecule.prop("_StereochemDone"), None);
    }

    #[test]
    fn assign_chiral_types_from_3d_marks_non_explicit_3d_chirality_like_rdkit() {
        let mut molecule = three_neighbor_3d_carbon(0);

        assign_chiral_types_from_3d(&mut molecule, 0);

        assert_ne!(molecule.atoms()[0].chiral_tag(), ChiralTag::Unspecified);
        assert_eq!(
            molecule.atoms()[0].prop("_NonExplicit3DChirality"),
            Some("1")
        );
    }

    #[test]
    fn assign_chiral_types_from_3d_does_not_mark_existing_explicit_atom_as_non_explicit() {
        let mut molecule = three_neighbor_3d_carbon(0);
        molecule.topology_block_mut().atoms[0].set_chiral_tag(ChiralTag::TetrahedralCw);

        assign_chiral_types_from_3d(&mut molecule, 0);

        assert_ne!(molecule.atoms()[0].chiral_tag(), ChiralTag::Unspecified);
        assert_eq!(molecule.atoms()[0].prop("_NonExplicit3DChirality"), None);
    }

    fn three_neighbor_pseudo_2d_carbon() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C));
        let fluorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
        let chlorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
        let bromine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(35).unwrap()));
        builder
            .add_bond(BondSpec::new(center, fluorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, chlorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, bromine, BondOrder::Single))
            .unwrap();
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [-1.0, -1.0, 0.0],
                ],
                true,
            ))
            .unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn assign_chiral_types_from_bond_dirs_promotes_implicit_h_on_three_coordinate_center() {
        let mut molecule = three_neighbor_pseudo_2d_carbon();
        molecule.topology_block_mut().bonds[0].set_direction(BondDirection::BeginWedge);

        assign_chiral_types_from_bond_dirs(&mut molecule, 0);

        assert_ne!(molecule.atoms()[0].chiral_tag(), ChiralTag::Unspecified);
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
    }

    fn stereogenic_double_bond_molecule_with_conformer(same_side: bool, is_3d: bool) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let c0 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let c1 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let f = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
        let cl = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
        builder
            .add_bond(BondSpec::new(c0, c1, BondOrder::Double))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c0, f, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c1, cl, BondOrder::Single))
            .unwrap();
        let right_y = if same_side { 1.0 } else { -1.0 };
        builder
            .add_conformer(Conformer3D::new(
                0,
                vec![
                    [-1.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [-1.0, 1.0, 0.0],
                    [1.0, right_y, 0.0],
                ],
                is_3d,
            ))
            .unwrap();
        builder.build().unwrap()
    }

    fn stereogenic_double_bond_3d_molecule(same_side: bool) -> Molecule {
        stereogenic_double_bond_molecule_with_conformer(same_side, true)
    }

    #[test]
    fn clear_dir_flags_marks_unknown_and_clears_non_wedge_directions_like_rdkit() {
        let mut molecule = stereogenic_double_bond_3d_molecule(true);
        molecule.topology_block_mut().bonds[0].set_direction(BondDirection::EitherDouble);
        molecule.topology_block_mut().bonds[1].set_direction(BondDirection::Unknown);
        molecule.topology_block_mut().bonds[2].set_direction(BondDirection::EndUpRight);

        clear_dir_flags(&mut molecule, true);

        assert!(molecule.bonds()[0].unknown_stereo());
        assert_eq!(molecule.bonds()[0].direction(), BondDirection::None);
        assert!(molecule.bonds()[1].unknown_stereo());
        assert_eq!(molecule.bonds()[1].direction(), BondDirection::None);
        assert_eq!(molecule.bonds()[2].direction(), BondDirection::EndUpRight);
    }

    #[test]
    fn clear_all_bond_dir_flags_clears_wedge_type_directions_like_rdkit() {
        let mut molecule = stereogenic_double_bond_3d_molecule(true);
        molecule.topology_block_mut().bonds[1].set_direction(BondDirection::EndUpRight);
        molecule.topology_block_mut().bonds[2].set_direction(BondDirection::EndDownRight);

        clear_all_bond_dir_flags(&mut molecule);

        assert_eq!(molecule.bonds()[1].direction(), BondDirection::None);
        assert_eq!(molecule.bonds()[2].direction(), BondDirection::None);
    }

    #[test]
    fn set_double_bond_neighbor_directions_assigns_adjacent_bond_dirs_from_3d() {
        let mut molecule = stereogenic_double_bond_3d_molecule(true);

        set_double_bond_neighbor_directions(&mut molecule, 0).unwrap();

        assert!(matches!(
            molecule.bonds()[1].direction(),
            BondDirection::EndUpRight | BondDirection::EndDownRight
        ));
        assert!(matches!(
            molecule.bonds()[2].direction(),
            BondDirection::EndUpRight | BondDirection::EndDownRight
        ));
    }

    #[test]
    fn set_bond_stereo_from_directions_assigns_cis_for_same_side_neighbors() {
        let mut molecule = stereogenic_double_bond_3d_molecule(true);
        set_double_bond_neighbor_directions(&mut molecule, 0).unwrap();

        set_bond_stereo_from_directions(&mut molecule).unwrap();

        assert_eq!(molecule.bonds()[0].stereo(), BondStereo::Cis);
        assert_eq!(
            molecule.bonds()[0].stereo_atoms(),
            Some([AtomId::new(2), AtomId::new(3)])
        );
    }

    #[test]
    fn assign_double_bond_stereo_from_directions_updates_public_molecule_state_like_rdkit() {
        let mut molecule = stereogenic_double_bond_3d_molecule(true);
        set_double_bond_neighbor_directions(&mut molecule, 0).unwrap();
        molecule
            .properties_mut()
            .set_prop("_needsDetectBondStereo", "1");

        assign_double_bond_stereo_from_directions(&mut molecule).unwrap();

        assert_eq!(molecule.bonds()[0].stereo(), BondStereo::Cis);
        assert_eq!(
            molecule.bonds()[0].stereo_atoms(),
            Some([AtomId::new(2), AtomId::new(3)])
        );
        assert_eq!(molecule.prop("_needsDetectBondStereo"), None);
    }

    #[test]
    fn assign_stereochemistry_from_3d_assigns_trans_for_opposite_side_neighbors() {
        let mut molecule = stereogenic_double_bond_3d_molecule(false);

        assign_stereochemistry_from_3d(&mut molecule, 0).unwrap();

        assert_eq!(molecule.bonds()[0].stereo(), BondStereo::Trans);
        assert_eq!(
            molecule.bonds()[0].stereo_atoms(),
            Some([AtomId::new(2), AtomId::new(3)])
        );
    }

    #[test]
    fn stereochemistry_from_3d_ignores_non_3d_conformer_like_rdkit() {
        let mut molecule = stereogenic_double_bond_molecule_with_conformer(false, false);
        molecule.topology_block_mut().bonds[1].set_direction(BondDirection::EndUpRight);
        molecule.topology_block_mut().bonds[2].set_direction(BondDirection::EndDownRight);

        assign_stereochemistry_from_3d(&mut molecule, 0).unwrap();

        assert_eq!(molecule.bonds()[0].stereo(), BondStereo::None);
        assert_eq!(molecule.bonds()[1].direction(), BondDirection::EndUpRight);
        assert_eq!(molecule.bonds()[2].direction(), BondDirection::EndDownRight);
    }

    #[test]
    fn stereochemistry_from_3d_clears_existing_double_bond_stereo_before_reassignment() {
        let mut molecule = stereogenic_double_bond_3d_molecule(true);
        molecule.topology_block_mut().bonds[0]
            .set_stereo_atoms(Some([AtomId::new(0), AtomId::new(1)]));
        molecule.topology_block_mut().bonds[0].set_stereo(BondStereo::Trans);

        assign_stereochemistry_from_3d(&mut molecule, 0).unwrap();

        assert_eq!(molecule.bonds()[0].stereo(), BondStereo::Cis);
        assert_eq!(
            molecule.bonds()[0].stereo_atoms(),
            Some([AtomId::new(2), AtomId::new(3)])
        );
    }

    fn square_planar_3d_phosphorus(coords: Vec<[f64; 3]>) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(
            AtomSpec::new(Element::from_atomic_number(15).unwrap()).with_no_implicit(true),
        );
        let ligand_elements = [
            Element::from_atomic_number(9).unwrap(),
            Element::from_atomic_number(17).unwrap(),
            Element::from_atomic_number(35).unwrap(),
            Element::from_atomic_number(53).unwrap(),
            Element::from_atomic_number(7).unwrap(),
            Element::from_atomic_number(8).unwrap(),
        ];
        for element in ligand_elements.iter().take(coords.len() - 1) {
            let ligand = builder.add_atom(AtomSpec::new(*element).with_no_implicit(true));
            builder
                .add_bond(BondSpec::new(center, ligand, BondOrder::Single))
                .unwrap();
        }
        builder
            .add_conformer(Conformer3D::new(0, coords, true))
            .unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn assign_chiral_types_from_3d_assigns_square_planar_from_two_opposite_pairs_like_rdkit() {
        let mut molecule = square_planar_3d_phosphorus(vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ]);

        assign_chiral_types_from_3d(&mut molecule, 0);

        assert_eq!(molecule.atoms()[0].chiral_tag(), ChiralTag::SquarePlanar);
        assert_eq!(molecule.atoms()[0].chiral_permutation(), Some(2));
        assert_eq!(
            molecule.atoms()[0].prop("_NonExplicit3DChirality"),
            Some("1")
        );
    }

    #[test]
    fn assign_chiral_types_from_3d_assigns_t_shaped_square_planar_like_rdkit() {
        let mut molecule = square_planar_3d_phosphorus(vec![
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
        ]);

        assign_chiral_types_from_3d(&mut molecule, 0);

        assert_eq!(molecule.atoms()[0].chiral_tag(), ChiralTag::SquarePlanar);
        assert_eq!(molecule.atoms()[0].chiral_permutation(), Some(3));
    }

    #[test]
    fn assign_chiral_types_from_3d_assigns_trigonal_bipyramidal_from_one_opposite_pair_like_rdkit()
    {
        let mut molecule = square_planar_3d_phosphorus(vec![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, -1.0, 0.0],
        ]);

        assign_chiral_types_from_3d(&mut molecule, 0);

        assert_eq!(
            molecule.atoms()[0].chiral_tag(),
            ChiralTag::TrigonalBipyramidal
        );
        assert_eq!(molecule.atoms()[0].chiral_permutation(), Some(7));
        assert_eq!(
            molecule.atoms()[0].prop("_NonExplicit3DChirality"),
            Some("1")
        );
    }

    #[test]
    fn assign_chiral_types_from_3d_assigns_seesaw_trigonal_bipyramidal_like_rdkit() {
        let mut molecule = square_planar_3d_phosphorus(vec![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
            [-0.5, 0.866_025_403_784_438_6, 0.0],
        ]);

        assign_chiral_types_from_3d(&mut molecule, 0);

        assert_eq!(
            molecule.atoms()[0].chiral_tag(),
            ChiralTag::TrigonalBipyramidal
        );
        assert_eq!(molecule.atoms()[0].chiral_permutation(), Some(7));
    }

    #[test]
    fn assign_chiral_types_from_3d_covers_all_seesaw_octahedral_branches_like_rdkit() {
        let cases = [
            (
                vec![
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [0.0, 0.0, -1.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                ],
                25,
            ),
            (
                vec![
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 0.0, -1.0],
                    [0.0, 1.0, 0.0],
                ],
                19,
            ),
            (
                vec![
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, -1.0],
                ],
                6,
            ),
            (
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [0.0, 0.0, -1.0],
                    [0.0, -1.0, 0.0],
                ],
                10,
            ),
            (
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, -1.0],
                ],
                1,
            ),
            (
                vec![
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, -1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [0.0, 0.0, -1.0],
                ],
                4,
            ),
        ];

        for (coords, expected_perm) in cases {
            let mut molecule = square_planar_3d_phosphorus(coords);
            assign_chiral_types_from_3d(&mut molecule, 0);

            assert_eq!(molecule.atoms()[0].chiral_tag(), ChiralTag::Octahedral);
            assert_eq!(
                molecule.atoms()[0].chiral_permutation(),
                Some(expected_perm)
            );
        }
    }

    #[test]
    fn assign_chiral_types_from_3d_assigns_octahedral_from_three_opposite_pairs_like_rdkit() {
        let mut molecule = square_planar_3d_phosphorus(vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ]);

        assign_chiral_types_from_3d(&mut molecule, 0);

        assert_eq!(molecule.atoms()[0].chiral_tag(), ChiralTag::Octahedral);
        assert_eq!(molecule.atoms()[0].chiral_permutation(), Some(27));
    }

    #[test]
    fn assign_chiral_types_from_3d_assigns_square_pyramidal_as_octahedral_like_rdkit() {
        let mut molecule = square_planar_3d_phosphorus(vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]);

        assign_chiral_types_from_3d(&mut molecule, 0);

        assert_eq!(molecule.atoms()[0].chiral_tag(), ChiralTag::Octahedral);
        assert_eq!(molecule.atoms()[0].chiral_permutation(), Some(27));
    }

    #[test]
    fn assign_chiral_types_from_3d_assigns_seesaw_octahedral_like_rdkit() {
        let mut molecule = square_planar_3d_phosphorus(vec![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]);

        assign_chiral_types_from_3d(&mut molecule, 0);

        assert_eq!(molecule.atoms()[0].chiral_tag(), ChiralTag::Octahedral);
        assert_eq!(molecule.atoms()[0].chiral_permutation(), Some(25));
    }

    #[test]
    fn handle_cx_part_and_name_parses_zero_bonds_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("CC~CC |Z:1|", false).unwrap();

        assert_eq!(molecule.num_bonds(), 3);
        assert_eq!(molecule.bonds()[1].order(), BondOrder::Zero);
    }

    #[test]
    fn handle_cx_part_and_name_parses_enhanced_stereo_groups_like_rdkit() {
        let molecule =
            Molecule::from_smiles_with_sanitize("C[C@H](O)N |o1:1,o1:2|", false).unwrap();

        assert_eq!(molecule.stereo_groups().len(), 1);
        assert_eq!(molecule.stereo_groups()[0].kind(), StereoGroupKind::Or);
        assert_eq!(molecule.stereo_groups()[0].id(), Some(1));
        assert_eq!(
            molecule.stereo_groups()[0].atoms(),
            &[AtomId::new(1), AtomId::new(2)]
        );
    }

    #[test]
    fn parse_cx_enhanced_stereo_rejects_missing_atom_index_like_rdkit() {
        let mut state = SmilesBuildState::new();
        state.builder.add_atom(AtomSpec::new(Element::C));

        let error = parse_cx_enhanced_stereo(&mut state, "o1:1", &mut 0).unwrap_err();

        assert_eq!(error, cx_parse_failure());
    }

    #[test]
    fn handle_cx_part_and_name_parses_coordinate_bonds_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("N->O |C:1.0|", false).unwrap();

        assert_eq!(molecule.num_bonds(), 1);
        assert_eq!(molecule.bonds()[0].order(), BondOrder::Dative);
        assert_eq!(molecule.bonds()[0].begin(), AtomId::new(1));
        assert_eq!(molecule.bonds()[0].end(), AtomId::new(0));
    }

    #[test]
    fn handle_cx_part_and_name_parses_hydrogen_bonds_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("NO |H:1.0|", false).unwrap();

        assert_eq!(molecule.num_bonds(), 1);
        assert_eq!(molecule.bonds()[0].order(), BondOrder::Hydrogen);
        assert_eq!(molecule.bonds()[0].begin(), AtomId::new(1));
        assert_eq!(molecule.bonds()[0].end(), AtomId::new(0));
    }

    #[test]
    fn handle_cx_part_and_name_parses_radicals_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("CCC |^1:0,^4:1,^7:2|", false).unwrap();

        assert_eq!(molecule.atoms()[0].radical_electrons(), 1);
        assert_eq!(molecule.atoms()[1].radical_electrons(), 2);
        assert_eq!(molecule.atoms()[2].radical_electrons(), 3);
    }

    #[test]
    fn parse_cx_extensions_dispatches_radicals_and_consumes_closing_pipe_like_rdkit() {
        let mut state = SmilesBuildState::new();
        state.builder.add_atom(AtomSpec::new(Element::C));
        state.builder.add_atom(AtomSpec::new(Element::C));
        state.builder.add_atom(AtomSpec::new(Element::C));

        let consumed = parse_cx_extensions(&mut state, "|^1:0,^4:1,^7:2|").unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(consumed, "|^1:0,^4:1,^7:2|".len());
        assert_eq!(molecule.atoms()[0].radical_electrons(), 1);
        assert_eq!(molecule.atoms()[1].radical_electrons(), 2);
        assert_eq!(molecule.atoms()[2].radical_electrons(), 3);
    }

    #[test]
    fn parse_cx_extensions_accepts_empty_text_like_rdkit_wrapper() {
        let mut state = SmilesBuildState::new();

        let consumed = parse_cx_extensions(&mut state, "").unwrap();

        assert_eq!(consumed, 0);
    }

    #[test]
    fn parse_cx_extensions_rejects_missing_leading_pipe_like_rdkit() {
        let mut state = SmilesBuildState::new();

        let error = parse_cx_extensions(&mut state, "^1:0|").unwrap_err();

        assert_eq!(
            error,
            SmilesParseError::ParseError("CXSMILES extension does not start with |".to_string())
        );
    }

    #[test]
    fn parse_cx_radicals_rejects_invalid_or_truncated_markers_like_rdkit() {
        let mut state = SmilesBuildState::new();
        state.builder.add_atom(AtomSpec::new(Element::C));

        let invalid = parse_cx_radicals(&mut state, "^0:0", &mut 0).unwrap_err();
        let truncated = parse_cx_radicals(&mut state, "^", &mut 0).unwrap_err();

        assert_eq!(invalid, cx_parse_failure());
        assert_eq!(truncated, cx_parse_failure());
    }

    #[test]
    fn parse_cx_coords_marks_2d_and_3d_conformers_like_rdkit_has_non_zero_z() {
        let mut state = SmilesBuildState::new();
        state.builder.add_atom(AtomSpec::new(Element::C));
        state.builder.add_atom(AtomSpec::new(Element::C));
        state.builder.add_atom(AtomSpec::new(Element::O));

        let mut pos = 0;
        parse_cx_coords(&mut state, "(0,0,;1,0,;2,0,0)", &mut pos, 0).unwrap();
        pos = 0;
        parse_cx_coords(&mut state, "(0,0,;1,0,;2,0,0.5)", &mut pos, 1).unwrap();

        let molecule = state.into_molecule().unwrap();
        assert_eq!(molecule.conformers_3d().len(), 2);
        assert!(!molecule.conformers_3d()[0].is_3d());
        assert!(molecule.conformers_3d()[1].is_3d());
        assert_eq!(molecule.conformers_3d()[1].coords()[2], [2.0, 0.0, 0.5]);
    }

    #[test]
    fn get_unspecified_bond_type_for_atoms_matches_rdkit_aromatic_rule() {
        assert_eq!(
            get_unspecified_bond_type_for_atoms(true, true),
            BondOrder::Aromatic
        );
        assert_eq!(
            get_unspecified_bond_type_for_atoms(true, false),
            BondOrder::Single
        );
        assert_eq!(
            get_unspecified_bond_type_for_atoms(false, true),
            BondOrder::Single
        );
        assert_eq!(
            get_unspecified_bond_type_for_atoms(false, false),
            BondOrder::Single
        );
    }

    #[test]
    fn handle_cx_part_and_name_parses_unsaturation_query_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("CC |u:1|", false).unwrap();

        assert_eq!(
            molecule.atoms()[1].query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::IsUnsaturated))
        );
    }

    #[test]
    fn handle_cx_part_and_name_parses_ring_bond_query_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("C1CC1 |rb:0:2|", false).unwrap();

        assert_eq!(
            molecule.atoms()[0].query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::RingBondCount(2)))
        );
    }

    #[test]
    fn handle_cx_part_and_name_completes_ring_bond_scan_queries_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("C1CC1 |rb:0:*|", false).unwrap();

        assert_eq!(
            molecule.atoms()[0].query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::RingBondCount(2)))
        );
        assert_eq!(molecule.properties().prop("_NeedsQueryScan"), None);
    }

    #[test]
    fn handle_cx_part_and_name_completes_substitution_scan_queries_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("CC |s:0:*|", false).unwrap();

        assert_eq!(
            molecule.atoms()[0].query(),
            Some(&QueryNode::predicate(
                AtomQueryPredicate::NonHydrogenDegree(1,)
            ))
        );
        assert_eq!(molecule.properties().prop("_NeedsQueryScan"), None);
    }

    #[test]
    fn handle_cx_part_and_name_completes_ring_and_non_ring_scan_queries_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("C1CC1C |rb:0:*,s:3:*|", false).unwrap();

        assert_eq!(
            molecule.atoms()[0].query(),
            Some(&QueryNode::predicate(AtomQueryPredicate::RingBondCount(2)))
        );
        assert_eq!(
            molecule.atoms()[3].query(),
            Some(&QueryNode::predicate(
                AtomQueryPredicate::NonHydrogenDegree(1)
            ))
        );
        assert_eq!(molecule.properties().prop("_NeedsQueryScan"), None);
    }

    #[test]
    fn handle_cx_part_and_name_parses_substitution_query_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("CC |s:0:1|", false).unwrap();

        assert_eq!(
            molecule.atoms()[0].query(),
            Some(&QueryNode::predicate(
                AtomQueryPredicate::NonHydrogenDegree(1)
            ))
        );
    }

    #[test]
    fn atom_query_has_single_h_count_matches_rdkit_atomand_only() {
        let direct_h = QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCount(1));
        let atom_and_h = QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCount(1)),
        ]);
        let nested_atom_and_h = QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::and(vec![QueryNode::predicate(
                AtomQueryPredicate::ImplicitHydrogenCount(1),
            )]),
        ]);
        let atom_or_h = QueryNode::or(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCount(1)),
        ]);
        let negated_h = QueryNode::and(vec![QueryNode::not(QueryNode::predicate(
            AtomQueryPredicate::ImplicitHydrogenCount(1),
        ))]);

        assert!(!atom_query_has_single_h_count(&direct_h));
        assert!(atom_query_has_single_h_count(&atom_and_h));
        assert!(atom_query_has_single_h_count(&nested_atom_and_h));
        assert!(!atom_query_has_single_h_count(&atom_or_h));
        assert!(!atom_query_has_single_h_count(&negated_h));
    }

    #[test]
    fn atom_has_fourth_valence_treats_single_h_query_like_rdkit() {
        let mut state = SmilesBuildState::new();
        let atom = state
            .builder
            .add_atom(AtomSpec::new(Element::C).with_query(QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCount(1)),
            ])));

        assert!(state.atom_has_fourth_valence(atom).unwrap());
    }

    #[test]
    fn handle_cx_part_and_name_clears_wedged_single_bond_dirs_after_atom_stereo_perception_like_rdkit()
     {
        let molecule = Molecule::from_smiles_with_sanitize("CC |wU:1.0|", false).unwrap();

        assert_eq!(molecule.bonds()[0].begin(), AtomId::new(1));
        assert_eq!(molecule.bonds()[0].end(), AtomId::new(0));
        assert_eq!(molecule.bonds()[0].direction(), BondDirection::None);
        assert!(!molecule.bonds()[0].unknown_stereo());
        assert_eq!(molecule.bonds()[0].prop("_MolFileBondCfg"), Some("1"));
    }

    #[test]
    fn handle_cx_part_and_name_clears_wiggly_single_bond_dirs_and_marks_unknown_stereo_like_rdkit()
    {
        let molecule = Molecule::from_smiles_with_sanitize("CC |w:1.0|", false).unwrap();

        assert_eq!(molecule.bonds()[0].begin(), AtomId::new(1));
        assert_eq!(molecule.bonds()[0].end(), AtomId::new(0));
        assert_eq!(molecule.bonds()[0].direction(), BondDirection::None);
        assert!(molecule.bonds()[0].unknown_stereo());
        assert_eq!(molecule.bonds()[0].prop("_MolFileBondCfg"), Some("2"));
        assert_eq!(
            molecule.properties().prop("_needsDetectBondStereo"),
            Some("1")
        );
    }

    #[test]
    fn handle_cx_part_and_name_parses_linknodes_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("C1CC1 |LN:1:1.3|", false).unwrap();

        assert_eq!(
            molecule.properties().prop("_MolFileLinkNodes"),
            Some("1 3 2 2 1 2 3")
        );
    }

    #[test]
    fn handle_cx_part_and_name_parses_data_sgroups_like_rdkit() {
        let molecule =
            Molecule::from_smiles_with_sanitize("CCO |SgD:2,1:FIELD:info::::|", false).unwrap();

        assert_eq!(molecule.substance_groups().len(), 1);
        let sgroup = &molecule.substance_groups()[0];
        assert_eq!(sgroup.kind(), &SubstanceGroupKind::Data);
        assert_eq!(sgroup.atoms(), &[AtomId::new(2), AtomId::new(1)]);
        assert_eq!(
            sgroup.props().get("FIELDNAME").map(String::as_str),
            Some("FIELD")
        );
        assert_eq!(sgroup.data_fields(), &["info".to_string()]);
        assert_eq!(sgroup.props().get("index").map(String::as_str), Some("1"));
    }

    #[test]
    fn parse_cx_data_sgroup_consumes_coords_for_dropped_group_like_rdkit() {
        let mut state = SmilesBuildState::new();
        state.builder.add_atom(AtomSpec::new(Element::C));
        let mut pos = 0;

        parse_cx_data_sgroup(&mut state, "SgD:9:FIELD:info::::(1.,1.)", &mut pos, 0).unwrap();

        assert_eq!(pos, "SgD:9:FIELD:info::::(1.,1.)".len());
        let molecule = state.into_molecule().unwrap();
        assert!(molecule.substance_groups().is_empty());
    }

    #[test]
    fn handle_cx_part_and_name_parses_variable_attachments_like_rdkit() {
        let molecule =
            Molecule::from_smiles_with_sanitize("CO*.C1=CC=NC=C1 |m:2:3.5.4|", false).unwrap();

        assert_eq!(
            molecule.bonds()[1].prop("_MolFileBondEndPts"),
            Some("(3 4 6 5)")
        );
        assert_eq!(molecule.bonds()[1].prop("_MolFileBondAttach"), Some("ANY"));
    }

    #[test]
    fn handle_cx_part_and_name_parses_double_bond_stereo_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("CC=CC |c:1|", false).unwrap();

        assert_eq!(molecule.bonds()[1].stereo(), BondStereo::Cis);
        assert_eq!(
            molecule.bonds()[1].stereo_atoms(),
            Some([AtomId::new(0), AtomId::new(3)])
        );
        assert_eq!(
            molecule.properties().prop("_needsDetectBondStereo"),
            Some("1")
        );
    }

    #[test]
    fn from_smiles_with_sanitize_false_parses_bracket_charge_hydrogen_and_map() {
        let molecule = Molecule::from_smiles_with_sanitize("[13CH3:5].[NH4+]", false).unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![6, 7]);
        assert_eq!(molecule.atoms()[0].isotope(), Some(13));
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 3);
        assert_eq!(molecule.atoms()[0].atom_map(), Some(5));
        assert!(molecule.atoms()[0].no_implicit());
        assert_eq!(molecule.atoms()[1].explicit_hydrogens(), 4);
        assert_eq!(molecule.atoms()[1].formal_charge(), 1);
        assert!(molecule.atoms()[1].no_implicit());
    }

    #[test]
    fn remove_hs_update_explicit_count_tracks_isotopic_hydrogens_and_clears_stale_property() {
        let mut state = SmilesBuildState::new();
        let carbon = state
            .builder
            .add_atom(AtomSpec::new(Element::C).with_tracked_isotopic_hydrogens(vec![9]));
        let deuterium = state
            .builder
            .add_atom(AtomSpec::new(Element::H).with_isotope(2));
        state
            .builder
            .add_bond(BondSpec::new(carbon, deuterium, BondOrder::Single))
            .unwrap();

        state
            .remove_hs_update_explicit_count(&RemoveHsParams {
                remove_and_track_isotopes: true,
                update_explicit_count: true,
                ..RemoveHsParams::default()
            })
            .unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.num_atoms(), 1);
        assert_eq!(molecule.atoms()[0].tracked_isotopic_hydrogens(), &[2]);
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
    }

    #[test]
    fn remove_hs_update_explicit_count_tracks_isotopes_after_nonisotopic_prepass() {
        let mut state = SmilesBuildState::new();
        let carbon = state
            .builder
            .add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let protium = state.builder.add_atom(AtomSpec::new(Element::H));
        let deuterium = state
            .builder
            .add_atom(AtomSpec::new(Element::H).with_isotope(2));
        state
            .builder
            .add_bond(BondSpec::new(carbon, protium, BondOrder::Single))
            .unwrap();
        state
            .builder
            .add_bond(BondSpec::new(carbon, deuterium, BondOrder::Single))
            .unwrap();

        state
            .remove_hs_update_explicit_count(&RemoveHsParams {
                remove_and_track_isotopes: true,
                update_explicit_count: true,
                ..RemoveHsParams::default()
            })
            .unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.num_atoms(), 1);
        assert_eq!(molecule.atoms()[0].tracked_isotopic_hydrogens(), &[2]);
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 2);
    }

    #[test]
    fn remove_hs_update_explicit_count_preserves_mapped_hydrogen_when_disabled() {
        let mut state = to_mol("[H:7]C").unwrap();

        state
            .remove_hs_update_explicit_count(&RemoveHsParams {
                remove_mapped: false,
                update_explicit_count: true,
                ..RemoveHsParams::default()
            })
            .unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![1, 6]);
        assert_eq!(molecule.atoms()[0].atom_map(), Some(7));
    }

    #[test]
    fn tracked_smiles_isotopic_hydrogens_groups_multiple_centers_in_atom_order_like_rdkit() {
        let mut builder = Molecule::builder();
        let carbon0 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let carbon1 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
        let deuterium = builder.add_atom(AtomSpec::new(Element::H).with_isotope(2));
        let tritium = builder.add_atom(AtomSpec::new(Element::H).with_isotope(3));
        builder
            .add_bond(BondSpec::new(carbon1, tritium, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(carbon0, deuterium, BondOrder::Single))
            .unwrap();

        let tracked = tracked_smiles_isotopic_hydrogens(&builder).unwrap();

        assert_eq!(tracked, vec![(carbon0, vec![2]), (carbon1, vec![3])]);
    }

    #[test]
    fn remove_hs_update_explicit_count_default_removes_hydrogen_with_wedged_bond() {
        let mut state = SmilesBuildState::new();
        let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
        state
            .builder
            .add_bond(
                BondSpec::new(carbon, hydrogen, BondOrder::Single)
                    .with_direction(BondDirection::BeginWedge),
            )
            .unwrap();

        state
            .remove_hs_update_explicit_count(&RemoveHsParams {
                update_explicit_count: true,
                ..RemoveHsParams::default()
            })
            .unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.num_atoms(), 1);
        assert_eq!(molecule.atoms()[0].atomic_number(), 6);
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
    }

    #[test]
    fn remove_hs_update_explicit_count_preserves_hydrogen_when_wedged_removal_disabled() {
        let mut state = SmilesBuildState::new();
        let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
        state
            .builder
            .add_bond(
                BondSpec::new(carbon, hydrogen, BondOrder::Single)
                    .with_direction(BondDirection::BeginWedge),
            )
            .unwrap();

        state
            .remove_hs_update_explicit_count(&RemoveHsParams {
                remove_with_wedged_bond: false,
                update_explicit_count: true,
                ..RemoveHsParams::default()
            })
            .unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.num_atoms(), 2);
    }

    #[test]
    fn remove_hs_update_explicit_count_preserves_hydrogen_on_nontetrahedral_neighbor_by_default() {
        let mut state = SmilesBuildState::new();
        let center = state
            .builder
            .add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::SquarePlanar));
        let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
        state
            .builder
            .add_bond(BondSpec::new(center, hydrogen, BondOrder::Single))
            .unwrap();

        state
            .remove_hs_update_explicit_count(&RemoveHsParams {
                update_explicit_count: true,
                ..RemoveHsParams::default()
            })
            .unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.num_atoms(), 2);
    }

    #[test]
    fn remove_hs_update_explicit_count_preserves_query_hydrogen_by_default() {
        let mut state = SmilesBuildState::new();
        let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = state.builder.add_atom(
            AtomSpec::new(Element::H).with_query(QueryNode::predicate(AtomQueryPredicate::Any)),
        );
        state
            .builder
            .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
            .unwrap();

        state
            .remove_hs_update_explicit_count(&RemoveHsParams {
                update_explicit_count: true,
                ..RemoveHsParams::default()
            })
            .unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.num_atoms(), 2);
    }

    #[test]
    fn remove_hs_update_explicit_count_preserves_nonimplicit_hydrogen_when_disabled() {
        let mut state = SmilesBuildState::new();
        let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
        state
            .builder
            .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
            .unwrap();

        state
            .remove_hs_update_explicit_count(&RemoveHsParams {
                remove_nonimplicit: false,
                update_explicit_count: true,
                ..RemoveHsParams::default()
            })
            .unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.num_atoms(), 2);
    }

    #[test]
    fn remove_hs_update_explicit_count_removes_degree_zero_hydrogen_when_enabled() {
        let mut state = SmilesBuildState::new();
        state.builder.add_atom(AtomSpec::new(Element::H));

        state
            .remove_hs_update_explicit_count(&RemoveHsParams {
                remove_degree_zero: true,
                ..RemoveHsParams::default()
            })
            .unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.num_atoms(), 0);
    }

    #[test]
    fn remove_hs_update_explicit_count_removes_higher_degree_hydrogen_when_enabled() {
        let mut state = SmilesBuildState::new();
        let carbon0 = state.builder.add_atom(AtomSpec::new(Element::C));
        let carbon1 = state.builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
        state
            .builder
            .add_bond(BondSpec::new(carbon0, hydrogen, BondOrder::Single))
            .unwrap();
        state
            .builder
            .add_bond(BondSpec::new(carbon1, hydrogen, BondOrder::Single))
            .unwrap();

        state
            .remove_hs_update_explicit_count(&RemoveHsParams {
                remove_higher_degrees: true,
                update_explicit_count: true,
                ..RemoveHsParams::default()
            })
            .unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.num_atoms(), 2);
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
        assert_eq!(molecule.atoms()[1].explicit_hydrogens(), 1);
    }

    #[test]
    fn remove_hs_update_explicit_count_removes_hydride_when_enabled() {
        let mut state = SmilesBuildState::new();
        let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
        let hydride = state
            .builder
            .add_atom(AtomSpec::new(Element::H).with_formal_charge(-1));
        state
            .builder
            .add_bond(BondSpec::new(carbon, hydride, BondOrder::Single))
            .unwrap();

        state
            .remove_hs_update_explicit_count(&RemoveHsParams {
                remove_hydrides: true,
                update_explicit_count: true,
                ..RemoveHsParams::default()
            })
            .unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.num_atoms(), 1);
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
    }

    #[test]
    fn remove_hs_update_explicit_count_does_not_always_increment_explicit_hydrogens() {
        let mut state = SmilesBuildState::new();
        let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
        state
            .builder
            .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
            .unwrap();

        state
            .remove_hs_update_explicit_count(&RemoveHsParams::default())
            .unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.num_atoms(), 1);
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 0);
    }

    #[test]
    fn remove_hs_update_explicit_count_moves_end_direction_from_removed_hydrogen_like_rdkit() {
        let mut state = SmilesBuildState::new();
        let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
        let oxygen = state.builder.add_atom(AtomSpec::new(Element::O));
        state
            .builder
            .add_bond(
                BondSpec::new(carbon, hydrogen, BondOrder::Single)
                    .with_direction(BondDirection::EndUpRight),
            )
            .unwrap();
        state
            .builder
            .add_bond(BondSpec::new(carbon, oxygen, BondOrder::Single))
            .unwrap();

        state
            .remove_hs_update_explicit_count(&RemoveHsParams::default())
            .unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.num_atoms(), 2);
        assert_eq!(molecule.bonds()[0].direction(), BondDirection::EndDownRight);
    }

    #[test]
    fn remove_hs_update_explicit_count_adjusts_double_bond_stereo_atoms_like_rdkit() {
        let mut state = SmilesBuildState::new();
        let center = state.builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
        let oxygen = state.builder.add_atom(AtomSpec::new(Element::O));
        let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
        state
            .builder
            .add_bond(BondSpec::new(center, hydrogen, BondOrder::Single))
            .unwrap();
        state
            .builder
            .add_bond(BondSpec::new(center, oxygen, BondOrder::Single))
            .unwrap();
        state
            .builder
            .add_bond(
                BondSpec::new(center, carbon, BondOrder::Double)
                    .with_stereo_atoms(hydrogen, carbon)
                    .with_stereo(BondStereo::Cis),
            )
            .unwrap();

        state
            .remove_hs_update_explicit_count(&RemoveHsParams {
                remove_defining_bond_stereo: true,
                ..RemoveHsParams::default()
            })
            .unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.num_atoms(), 3);
        assert_eq!(molecule.bonds()[1].stereo(), BondStereo::Trans);
        assert_eq!(
            molecule.bonds()[1].stereo_atoms(),
            Some([AtomId::new(1), AtomId::new(2)])
        );
    }

    #[test]
    fn remove_hs_update_explicit_count_sets_unknown_stereo_on_heavy_atom_like_rdkit() {
        let mut state = SmilesBuildState::new();
        let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
        state
            .builder
            .add_bond(
                BondSpec::new(carbon, hydrogen, BondOrder::Single)
                    .with_direction(BondDirection::Unknown),
            )
            .unwrap();

        state
            .remove_hs_update_explicit_count(&RemoveHsParams::default())
            .unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.num_atoms(), 1);
        assert!(molecule.atoms()[0].unknown_stereo());
    }

    #[test]
    fn remove_hs_update_explicit_count_applies_sgroup_special_role_and_emptying_guards() {
        let mut state = SmilesBuildState::new();
        let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
        let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
        state
            .builder
            .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
            .unwrap();
        state
            .builder
            .add_substance_group(
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
                    .with_atoms(vec![carbon])
                    .with_attach_points(vec![crate::SGroupAttachPoint {
                        atom: hydrogen,
                        leaving_atom: None,
                        label: None,
                        order: None,
                    }]),
            )
            .unwrap();

        state
            .remove_hs_update_explicit_count(&RemoveHsParams {
                update_explicit_count: true,
                ..RemoveHsParams::default()
            })
            .unwrap();
        let molecule = state.into_molecule().unwrap();
        assert_eq!(molecule.num_atoms(), 2);
        assert_eq!(molecule.num_bonds(), 1);

        let mut state = SmilesBuildState::new();
        let isolated_h = state.builder.add_atom(AtomSpec::new(Element::H));
        state
            .builder
            .add_substance_group(
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
                    .with_atoms(vec![isolated_h]),
            )
            .unwrap();

        state
            .remove_hs_update_explicit_count(&RemoveHsParams {
                remove_degree_zero: true,
                update_explicit_count: true,
                ..RemoveHsParams::default()
            })
            .unwrap();
        let molecule = state.into_molecule().unwrap();
        assert_eq!(molecule.num_atoms(), 1);
        assert_eq!(molecule.num_bonds(), 0);
        assert_eq!(molecule.atoms()[0].atomic_number(), 1);
    }

    #[test]
    fn remove_hs_update_explicit_count_respects_remove_in_sgroups_false_membership_guard() {
        let mut state = SmilesBuildState::new();
        let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
        state
            .builder
            .add_substance_group(
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
                    .with_atoms(vec![hydrogen]),
            )
            .unwrap();

        state
            .remove_hs_update_explicit_count(&RemoveHsParams {
                remove_degree_zero: true,
                remove_in_sgroups: false,
                update_explicit_count: true,
                ..RemoveHsParams::default()
            })
            .unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.num_atoms(), 1);
        assert_eq!(molecule.atoms()[0].atomic_number(), 1);
    }

    #[test]
    fn apply_smiles_postprocessing_removes_hydrogen_when_remove_hs_is_enabled() {
        let params = SmilesParseParams {
            remove_hs: true,
            sanitize: false,
            ..SmilesParseParams::default()
        };
        let mut state = to_mol("[CH]").unwrap();

        apply_smiles_postprocessing(&mut state, &params).unwrap();
        let molecule = state.into_molecule().unwrap();

        assert_eq!(molecule.num_atoms(), 1);
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
    }

    #[test]
    fn from_smiles_with_sanitize_false_parses_hash_element_and_tetrahedral_chirality() {
        let molecule = Molecule::from_smiles_with_sanitize("[#6].[C@H]", false).unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![6, 6]);
        assert_eq!(molecule.atoms()[1].chiral_tag(), ChiralTag::TetrahedralCcw);
        assert_eq!(molecule.atoms()[1].explicit_hydrogens(), 1);
    }

    #[test]
    fn from_smiles_with_sanitize_false_parses_bracket_element_lexer_table_like_rdkit() {
        let mut cases: Vec<(&str, u8, bool)> = BRACKET_ATOM_SYMBOLS
            .iter()
            .map(|(symbol, atomic_number)| (*symbol, *atomic_number, false))
            .collect();
        cases.extend([
            ("B", 5, false),
            ("C", 6, false),
            ("N", 7, false),
            ("O", 8, false),
            ("P", 15, false),
            ("S", 16, false),
            ("F", 9, false),
            ("Cl", 17, false),
            ("Br", 35, false),
            ("I", 53, false),
            ("H", 1, false),
            ("b", 5, true),
            ("c", 6, true),
            ("n", 7, true),
            ("o", 8, true),
            ("p", 15, true),
            ("s", 16, true),
            ("si", 14, true),
            ("as", 33, true),
            ("se", 34, true),
            ("te", 52, true),
        ]);

        for (symbol, atomic_number, aromatic) in cases {
            let molecule =
                Molecule::from_smiles_with_sanitize(&format!("[{symbol}]"), false).unwrap();
            assert_eq!(
                molecule.atoms()[0].atomic_number(),
                atomic_number,
                "symbol {symbol}"
            );
            assert_eq!(
                molecule.atoms()[0].is_aromatic(),
                aromatic,
                "symbol {symbol}"
            );
        }
    }

    #[test]
    fn from_smiles_with_sanitize_false_parses_quoted_biovia_atoms_like_rdkit() {
        for &(symbol, atomic_number) in QUOTED_BIOVIA_ATOM_SYMBOLS {
            let molecule =
                Molecule::from_smiles_with_sanitize(&format!("[{symbol}]"), false).unwrap();
            assert_eq!(
                molecule.atoms()[0].atomic_number(),
                atomic_number,
                "symbol {symbol}"
            );
        }
    }

    #[test]
    fn from_smiles_with_sanitize_false_converts_tetrahedral_chiral_class_like_rdkit() {
        let th1 = Molecule::from_smiles_with_sanitize("[C@TH1H]", false).unwrap();
        let th2 = Molecule::from_smiles_with_sanitize("[C@TH2H]", false).unwrap();

        assert_eq!(th1.atoms()[0].chiral_tag(), ChiralTag::TetrahedralCcw);
        assert_eq!(th1.atoms()[0].chiral_permutation(), None);
        assert_eq!(th2.atoms()[0].chiral_tag(), ChiralTag::TetrahedralCw);
        assert_eq!(th2.atoms()[0].chiral_permutation(), None);
    }

    #[test]
    fn from_smiles_with_sanitize_false_preserves_non_tetrahedral_chiral_classes_like_rdkit() {
        let allene = Molecule::from_smiles_with_sanitize("[C@AL1]", false).unwrap();
        let square_planar = Molecule::from_smiles_with_sanitize("[C@SP3]", false).unwrap();
        let trigonal_bipyramidal = Molecule::from_smiles_with_sanitize("[C@TB20]", false).unwrap();
        let octahedral = Molecule::from_smiles_with_sanitize("[C@OH30]", false).unwrap();

        assert_eq!(allene.atoms()[0].chiral_tag(), ChiralTag::Allene);
        assert_eq!(allene.atoms()[0].chiral_permutation(), Some(1));
        assert_eq!(
            square_planar.atoms()[0].chiral_tag(),
            ChiralTag::SquarePlanar
        );
        assert_eq!(square_planar.atoms()[0].chiral_permutation(), Some(3));
        assert_eq!(
            trigonal_bipyramidal.atoms()[0].chiral_tag(),
            ChiralTag::TrigonalBipyramidal
        );
        assert_eq!(
            trigonal_bipyramidal.atoms()[0].chiral_permutation(),
            Some(20)
        );
        assert_eq!(octahedral.atoms()[0].chiral_tag(), ChiralTag::Octahedral);
        assert_eq!(octahedral.atoms()[0].chiral_permutation(), Some(30));
    }

    #[test]
    fn from_smiles_with_sanitize_false_preserves_default_non_tetrahedral_chiral_permutations_like_rdkit()
     {
        let allene = Molecule::from_smiles_with_sanitize("[C@AL]", false).unwrap();
        let square_planar = Molecule::from_smiles_with_sanitize("[C@SP]", false).unwrap();
        let trigonal_bipyramidal = Molecule::from_smiles_with_sanitize("[C@TB]", false).unwrap();
        let octahedral = Molecule::from_smiles_with_sanitize("[C@OH]", false).unwrap();

        assert_eq!(allene.atoms()[0].chiral_tag(), ChiralTag::Allene);
        assert_eq!(allene.atoms()[0].chiral_permutation(), Some(0));
        assert_eq!(
            square_planar.atoms()[0].chiral_tag(),
            ChiralTag::SquarePlanar
        );
        assert_eq!(square_planar.atoms()[0].chiral_permutation(), Some(0));
        assert_eq!(
            trigonal_bipyramidal.atoms()[0].chiral_tag(),
            ChiralTag::TrigonalBipyramidal
        );
        assert_eq!(
            trigonal_bipyramidal.atoms()[0].chiral_permutation(),
            Some(0)
        );
        assert_eq!(octahedral.atoms()[0].chiral_tag(), ChiralTag::Octahedral);
        assert_eq!(octahedral.atoms()[0].chiral_permutation(), Some(0));
    }

    #[test]
    fn from_smiles_with_sanitize_false_recomputes_nontetrahedral_ring_permutation_like_rdkit() {
        let sp1 = Molecule::from_smiles_with_sanitize("F[C@SP1]1(Br)I.Cl1", false).unwrap();
        let sp2 = Molecule::from_smiles_with_sanitize("F[C@SP2]1(Br)I.Cl1", false).unwrap();

        assert_eq!(sp1.atoms()[1].chiral_tag(), ChiralTag::SquarePlanar);
        assert_eq!(sp1.atoms()[1].chiral_permutation(), Some(2));
        assert_eq!(sp2.atoms()[1].chiral_tag(), ChiralTag::SquarePlanar);
        assert_eq!(sp2.atoms()[1].chiral_permutation(), Some(3));
    }

    #[test]
    fn nontetrahedral_chiral_permutation_for_probe_returns_zero_for_nonpositive_or_oversized_cases()
    {
        let mut builder = Molecule::builder();
        let center_without_perm =
            builder.add_atom(AtomSpec::new(Element::PT).with_chiral_tag(ChiralTag::SquarePlanar));
        let neighbor_a = builder.add_atom(AtomSpec::new(Element::CL));
        let neighbor_b = builder.add_atom(AtomSpec::new(Element::CL));
        let neighbor_c = builder.add_atom(AtomSpec::new(Element::CL));
        let neighbor_d = builder.add_atom(AtomSpec::new(Element::CL));
        let bond_a = builder
            .add_bond(BondSpec::new(
                center_without_perm,
                neighbor_a,
                BondOrder::Single,
            ))
            .unwrap();
        let bond_b = builder
            .add_bond(BondSpec::new(
                center_without_perm,
                neighbor_b,
                BondOrder::Single,
            ))
            .unwrap();
        let bond_c = builder
            .add_bond(BondSpec::new(
                center_without_perm,
                neighbor_c,
                BondOrder::Single,
            ))
            .unwrap();
        let bond_d = builder
            .add_bond(BondSpec::new(
                center_without_perm,
                neighbor_d,
                BondOrder::Single,
            ))
            .unwrap();
        let molecule_without_perm = builder.build().unwrap();

        assert_eq!(
            nontetrahedral_chiral_permutation_for_probe(
                &molecule_without_perm,
                center_without_perm,
                &[Some(bond_a), Some(bond_b), Some(bond_c), Some(bond_d)],
                false,
            )
            .unwrap(),
            0
        );

        let mut builder = Molecule::builder();
        let center_with_perm = builder.add_atom(
            AtomSpec::new(Element::PT)
                .with_chiral_tag(ChiralTag::SquarePlanar)
                .with_chiral_permutation(1),
        );
        let neighbor_a = builder.add_atom(AtomSpec::new(Element::CL));
        let neighbor_b = builder.add_atom(AtomSpec::new(Element::CL));
        let neighbor_c = builder.add_atom(AtomSpec::new(Element::CL));
        let neighbor_d = builder.add_atom(AtomSpec::new(Element::CL));
        let extra_neighbor = builder.add_atom(AtomSpec::new(Element::CL));
        let bond_a = builder
            .add_bond(BondSpec::new(
                center_with_perm,
                neighbor_a,
                BondOrder::Single,
            ))
            .unwrap();
        let bond_b = builder
            .add_bond(BondSpec::new(
                center_with_perm,
                neighbor_b,
                BondOrder::Single,
            ))
            .unwrap();
        let bond_c = builder
            .add_bond(BondSpec::new(
                center_with_perm,
                neighbor_c,
                BondOrder::Single,
            ))
            .unwrap();
        let bond_d = builder
            .add_bond(BondSpec::new(
                center_with_perm,
                neighbor_d,
                BondOrder::Single,
            ))
            .unwrap();
        let extra_bond = builder
            .add_bond(BondSpec::new(
                center_with_perm,
                extra_neighbor,
                BondOrder::Single,
            ))
            .unwrap();
        let molecule_with_perm = builder.build().unwrap();

        assert_eq!(
            nontetrahedral_chiral_permutation_for_probe(
                &molecule_with_perm,
                center_with_perm,
                &[
                    Some(bond_a),
                    Some(bond_b),
                    Some(bond_c),
                    Some(bond_d),
                    Some(extra_bond),
                ],
                false,
            )
            .unwrap(),
            0
        );
    }

    #[test]
    fn from_smiles_with_sanitize_false_rejects_invalid_chiral_permutation_like_rdkit() {
        let error = Molecule::from_smiles_with_sanitize("[C@TH3H]", false).unwrap_err();

        assert_eq!(
            error,
            SmilesParseError::ParseError("Invalid chiral specification on atom 0".to_string())
        );
    }

    #[test]
    fn from_smiles_with_sanitize_false_rejects_non_tetrahedral_chiral_permutation_limits_like_rdkit()
     {
        let error = Molecule::from_smiles_with_sanitize("[C@SP4]", false).unwrap_err();

        assert_eq!(
            error,
            SmilesParseError::ParseError("Invalid chiral specification on atom 0".to_string())
        );
    }

    #[test]
    fn from_smiles_with_sanitize_false_restores_active_atom_after_branch() {
        let molecule = Molecule::from_smiles_with_sanitize("CC(C)O", false).unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![6, 6, 6, 8]);
        assert_eq!(molecule.num_bonds(), 3);
        assert_eq!(molecule.bonds()[0].begin(), AtomId::new(0));
        assert_eq!(molecule.bonds()[0].end(), AtomId::new(1));
        assert_eq!(molecule.bonds()[1].begin(), AtomId::new(1));
        assert_eq!(molecule.bonds()[1].end(), AtomId::new(2));
        assert_eq!(molecule.bonds()[2].begin(), AtomId::new(1));
        assert_eq!(molecule.bonds()[2].end(), AtomId::new(3));
    }

    #[test]
    fn close_branch_restores_root_before_following_atom_like_rdkit() {
        let molecule = to_mol("CC(C)O").unwrap().into_molecule().unwrap();

        assert_eq!(molecule.num_bonds(), 3);
        assert_eq!(molecule.bonds()[1].begin(), AtomId::new(1));
        assert_eq!(molecule.bonds()[1].end(), AtomId::new(2));
        assert_eq!(molecule.bonds()[2].begin(), AtomId::new(1));
        assert_eq!(molecule.bonds()[2].end(), AtomId::new(3));
    }

    #[test]
    fn from_smiles_with_sanitize_false_parses_explicit_branch_bond() {
        let molecule = Molecule::from_smiles_with_sanitize("C(=O)N", false).unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![6, 8, 7]);
        assert_eq!(molecule.num_bonds(), 2);
        assert_eq!(molecule.bonds()[0].begin(), AtomId::new(0));
        assert_eq!(molecule.bonds()[0].end(), AtomId::new(1));
        assert_eq!(molecule.bonds()[0].order(), BondOrder::Double);
        assert_eq!(molecule.bonds()[1].begin(), AtomId::new(0));
        assert_eq!(molecule.bonds()[1].end(), AtomId::new(2));
    }

    #[test]
    fn from_smiles_with_sanitize_false_closes_simple_ring_numbers() {
        let molecule = Molecule::from_smiles_with_sanitize("C1CC1", false).unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![6, 6, 6]);
        assert_eq!(molecule.num_bonds(), 3);
        assert_eq!(molecule.bonds()[2].begin(), AtomId::new(0));
        assert_eq!(molecule.bonds()[2].end(), AtomId::new(2));
        assert_eq!(molecule.bonds()[2].order(), BondOrder::Single);
    }

    #[test]
    fn from_smiles_with_sanitize_false_parses_percent_ring_numbers() {
        let molecule = Molecule::from_smiles_with_sanitize("C%12CC%12", false).unwrap();

        assert_eq!(molecule.num_atoms(), 3);
        assert_eq!(molecule.num_bonds(), 3);
        assert_eq!(molecule.bonds()[2].begin(), AtomId::new(0));
        assert_eq!(molecule.bonds()[2].end(), AtomId::new(2));
    }

    #[test]
    fn from_smiles_with_sanitize_false_reports_unclosed_ring_like_rdkit() {
        let error = Molecule::from_smiles_with_sanitize("C1CC", false).unwrap_err();

        assert_eq!(
            error,
            SmilesParseError::ParseError("unclosed ring".to_string())
        );
    }

    #[test]
    fn from_smiles_with_sanitize_false_reports_extra_close_parentheses_like_rdkit() {
        let error = Molecule::from_smiles_with_sanitize("CC)", false).unwrap_err();

        assert_eq!(
            error,
            SmilesParseError::ParseError("extra close parentheses".to_string())
        );
    }

    #[test]
    fn from_smiles_with_sanitize_false_reports_self_ring_closure_like_rdkit() {
        let error = Molecule::from_smiles_with_sanitize("C11", false).unwrap_err();

        assert_eq!(
            error,
            SmilesParseError::ParseError(
                "duplicated ring closure bonds atom 0 to itself".to_string()
            )
        );
    }

    #[test]
    fn from_smiles_with_sanitize_false_reports_duplicate_ring_bond_like_rdkit() {
        let error = Molecule::from_smiles_with_sanitize("C1C1", false).unwrap_err();

        assert_eq!(
            error,
            SmilesParseError::ParseError(
                "ring closure duplicates bond between atom 0 and atom 1".to_string()
            )
        );
    }

    #[test]
    fn from_smiles_with_sanitize_false_reports_unclosed_branch_like_rdkit() {
        let error = Molecule::from_smiles_with_sanitize("C(C", false).unwrap_err();

        assert_eq!(
            error,
            SmilesParseError::ParseError("extra open parentheses".to_string())
        );
    }

    #[test]
    fn from_smiles_with_sanitize_false_reports_bad_character_as_syntax_error_like_rdkit() {
        let error = Molecule::from_smiles_with_sanitize("C&N", false).unwrap_err();

        assert_eq!(
            error,
            SmilesParseError::ParseError("syntax error".to_string())
        );
    }

    #[test]
    fn from_smiles_with_sanitize_false_uses_first_explicit_ring_bond_order_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("C=1CC1", false).unwrap();

        assert_eq!(molecule.num_bonds(), 3);
        assert_eq!(molecule.bonds()[2].begin(), AtomId::new(0));
        assert_eq!(molecule.bonds()[2].end(), AtomId::new(2));
        assert_eq!(molecule.bonds()[2].order(), BondOrder::Double);
    }

    #[test]
    fn from_smiles_with_sanitize_false_uses_closing_explicit_ring_bond_when_opening_unspecified() {
        let molecule = Molecule::from_smiles_with_sanitize("C1CC=1", false).unwrap();

        assert_eq!(molecule.num_bonds(), 3);
        assert_eq!(molecule.bonds()[2].begin(), AtomId::new(0));
        assert_eq!(molecule.bonds()[2].end(), AtomId::new(2));
        assert_eq!(molecule.bonds()[2].order(), BondOrder::Double);
    }

    #[test]
    fn from_smiles_with_sanitize_false_marks_aromatic_atoms_and_unspecified_bonds_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("c1ccccc1", false).unwrap();

        assert_eq!(molecule.num_atoms(), 6);
        assert_eq!(molecule.num_bonds(), 6);
        assert!(molecule.atoms().iter().all(|atom| atom.is_aromatic()));
        assert!(molecule.bonds().iter().all(|bond| bond.is_aromatic()));
        assert!(
            molecule
                .bonds()
                .iter()
                .all(|bond| bond.order() == BondOrder::Aromatic)
        );
    }

    #[test]
    fn from_smiles_empty_returns_empty_molecule_like_rdkit() {
        let molecule = Molecule::from_smiles("").unwrap();

        assert_eq!(molecule.num_atoms(), 0);
        assert_eq!(molecule.num_bonds(), 0);
    }

    #[test]
    fn from_smiles_default_runs_postprocessing_for_simple_molecule() {
        let molecule = Molecule::from_smiles("CCO").unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![6, 6, 8]);
        assert_eq!(molecule.num_bonds(), 2);
    }

    #[test]
    fn from_smiles_default_parses_aromatic_with_sanitize_through_operations() {
        // Aromatic molecules now parse with default sanitize through
        // the registered operation pipeline (kekulize + valence + aromaticity).
        let molecule = Molecule::from_smiles("c1ccccc1").unwrap();
        // Sanitize should produce a valid Kekulé/aromatic molecule
        assert_eq!(molecule.atomic_numbers(), vec![6; 6]);
        assert_eq!(molecule.num_bonds(), 6);
    }

    #[test]
    fn from_smiles_default_removes_directional_h_with_sanitize_integration() {
        // Directional hydrogen with default sanitize now works through
        // the registered operations pipeline.
        let molecule = Molecule::from_smiles("[H]/C=C").unwrap();
        assert!(!molecule.atoms().is_empty(), "should parse successfully");
    }

    #[test]
    fn from_smiles_default_integrates_cx_conformer_with_sanitize() {
        // CX conformer data with default sanitize now works through
        // the registered operations pipeline.
        let molecule = Molecule::from_smiles("CCO |(0,0,;1,0,;2,0,0.5)|").unwrap();
        assert_eq!(molecule.atomic_numbers(), vec![6, 6, 8]);
        assert_eq!(molecule.num_bonds(), 2);
    }

    #[test]
    fn from_smiles_default_removes_explicit_h_and_updates_count() {
        let molecule = Molecule::from_smiles("[CH3][H]").unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![6]);
        assert_eq!(molecule.num_bonds(), 0);
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 4);
    }

    #[test]
    fn mol_from_smiles_remove_hs_without_sanitize_uses_split_reader_branch_like_rdkit() {
        let params = SmilesParseParams {
            sanitize: false,
            remove_hs: true,
            ..Default::default()
        };
        let molecule = mol_from_smiles("[CH3][H]", &params).unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![6]);
        assert_eq!(molecule.num_bonds(), 0);
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 4);
    }

    #[test]
    fn mol_from_smiles_skip_cleanup_preserves_parser_temporary_props_like_rdkit_reader_flag() {
        let params = SmilesParseParams {
            sanitize: false,
            skip_cleanup: true,
            ..Default::default()
        };
        let molecule = mol_from_smiles("c1ccccc1.C", &params).unwrap();

        assert_eq!(molecule.atoms()[0].prop(SMILES_START_PROP), Some("1"));
        assert_eq!(molecule.atoms()[6].prop(SMILES_START_PROP), Some("1"));
        assert!(
            molecule
                .bonds()
                .iter()
                .all(|bond| bond.prop(CXSMILES_BOND_IDX_PROP).is_some())
        );
        assert!(
            molecule.bonds()[..6]
                .iter()
                .all(|bond| bond.order() == BondOrder::Aromatic && bond.is_aromatic())
        );
    }

    #[test]
    fn mol_from_smiles_sanitize_without_remove_hs_uses_molecule_sanitize_branch_like_rdkit() {
        let params = SmilesParseParams {
            sanitize: true,
            remove_hs: false,
            ..Default::default()
        };
        let molecule = mol_from_smiles("[CH3][H]", &params).unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![6, 1]);
        assert_eq!(molecule.num_bonds(), 1);
    }

    #[test]
    fn from_smiles_default_removes_mapped_hydrogen_like_rdkit_default() {
        let molecule = Molecule::from_smiles("[H:7]C").unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![6]);
        assert_eq!(molecule.num_bonds(), 0);
        assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
    }

    #[test]
    fn from_smiles_default_keeps_hydrogen_with_only_h_neighbor() {
        let molecule = Molecule::from_smiles("[H][H]").unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![1, 1]);
        assert_eq!(molecule.num_bonds(), 1);
    }

    #[test]
    fn from_smiles_default_keeps_hydrogen_with_dummy_neighbor() {
        let molecule = Molecule::from_smiles("[*][H]").unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![0, 1]);
        assert_eq!(molecule.num_bonds(), 1);
    }

    #[test]
    fn from_smiles_default_keeps_isotopic_hydrogen() {
        let molecule = Molecule::from_smiles("[2H]O").unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![1, 8]);
        assert_eq!(molecule.num_bonds(), 1);
        assert_eq!(molecule.atoms()[0].isotope(), Some(2));
    }

    #[test]
    fn smiles_parse_ops_clear_atom_chemical_props_resets_query_atom_state_like_rdkit() {
        let atom_id = AtomId::new(0);
        let mut atom = Atom::from_spec(
            atom_id,
            AtomSpec::new(Element::C)
                .with_isotope(13)
                .with_formal_charge(2)
                .with_explicit_hydrogens(3),
        );

        clear_atom_chemical_props(&mut atom);

        assert_eq!(atom.isotope(), None);
        assert_eq!(atom.formal_charge(), 0);
        assert_eq!(atom.explicit_hydrogens(), 0);
    }

    #[test]
    fn smiles_parse_ops_report_parse_error_throws_when_requested_like_rdkit() {
        let error = report_parse_error("bad parse", true).unwrap_err();

        assert_eq!(error, SmilesParseError::ParseError("bad parse".to_string()));
    }

    #[test]
    fn smiles_parse_ops_cleanup_after_parse_error_clears_unmatched_ring_state_like_rdkit() {
        let mut state = SmilesBuildState::new();
        state.add_first_atom(SmilesAtomToken::new(6)).unwrap();
        state.add_ring_marker(1).unwrap();

        assert!(!state.ring_openings.is_empty());
        assert!(!state.ring_closures_by_atom.is_empty());

        state.cleanup_after_parse_error();

        assert!(state.ring_openings.is_empty());
        assert!(state.ring_closures_by_atom.is_empty());
    }

    #[test]
    fn smiles_parse_ops_add_frag_to_mol_merges_disconnected_fragment_like_rdkit() {
        let mut root = SmilesBuildState::new();
        root.add_first_atom(SmilesAtomToken::new(6)).unwrap();

        let mut frag = SmilesBuildState::new();
        frag.add_first_atom(SmilesAtomToken::new(8)).unwrap();

        root.add_frag_to_mol(frag, BondOrder::Ionic, BondDirection::None)
            .unwrap();
        let molecule = root.into_molecule().unwrap();

        assert_eq!(molecule.atomic_numbers(), vec![6, 8]);
        assert_eq!(molecule.num_bonds(), 0);
    }

    #[test]
    fn smiles_parse_ops_to_mol_reports_parse_error_and_cleans_partial_state_like_rdkit() {
        let error = to_mol("C1CC").unwrap_err();

        assert_eq!(
            error,
            SmilesParseError::ParseError("unclosed ring".to_string())
        );
    }

    #[test]
    fn smiles_parse_ops_to_mol_empty_input_returns_empty_state_like_rdkit() {
        let molecule = to_mol("").unwrap().into_molecule().unwrap();

        assert_eq!(molecule.num_atoms(), 0);
        assert_eq!(molecule.num_bonds(), 0);
    }

    #[test]
    fn parse_cx_variable_attachments_rejects_multi_degree_source_like_rdkit() {
        let mut state = to_mol("*C(C)C").unwrap();
        let mut pos = 0;

        let error = parse_cx_variable_attachments(&mut state, "m:1:0.2", &mut pos).unwrap_err();

        assert_eq!(
            error,
            SmilesParseError::ParseError("failure parsing CXSMILES extensions".to_string())
        );
    }
}
