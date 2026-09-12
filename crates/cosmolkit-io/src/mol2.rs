//! RDKit source-backed MOL2 parsing over detached model values.
//!
//! This module owns MOL2 syntax, Tripos-specific substructure cleanup, and
//! formal-charge inference. Generic sanitization, hydrogen removal, and stereo
//! finalization remain separate chemistry/runtime stages.

use cosmolkit_core::{
    bond_valence_contrib, find_sssr_from_parts, periodic_table_outer_electrons, rdkit_valence_list,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, Conformer3D, CoordinateBlock,
    CoordinateDimension, MoleculeProperties, TopologyBlock,
};
use cosmolkit_types::{BondOrder, Element};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum Mol2ReadError {
    #[error("MOL2 parse failed: {0}")]
    Parse(String),
    #[error("detached MOL2 reader does not support {feature}")]
    Unsupported { feature: &'static str },
    #[error("invalid detached MOL2 topology: {0}")]
    Topology(#[from] cosmolkit_model::TopologyValidationError),
    #[error("invalid detached MOL2 coordinates: {0}")]
    Coordinates(#[from] cosmolkit_model::CoordinateValidationError),
    #[error("invalid atom property: {0}")]
    AtomProperty(#[from] cosmolkit_model::AtomPropertyError),
    #[error("invalid molecule property: {0}")]
    MoleculeProperty(#[from] cosmolkit_model::MoleculePropertyError),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Mol2Type {
    // RDKit source: FileParsers.h `Mol2Type`
    // RDKit✔️✔️: typedef enum {
    // RDKit✔️✔️:   CORINA = 0  //!< supports output from Corina and some dbtranslate output
    // RDKit✔️✔️: } Mol2Type;
    #[default]
    Corina,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Mol2ReadParams {
    pub variant: Mol2Type,
    pub cleanup_substructures: bool,
}

impl Default for Mol2ReadParams {
    fn default() -> Self {
        // RDKit source: FileParsers.h `Mol2ParserParams`
        // RDKit✔️✔️:   Mol2Type variant = Mol2Type::CORINA; /**< the atom type definitions to use */
        // RDKit✔️✔️:   bool cleanupSubstructures =
        // RDKit✔️✔️:       true; /**< toggles recognition and cleanup of common substructures */
        Self {
            variant: Mol2Type::Corina,
            cleanup_substructures: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Mol2Record {
    pub topology: TopologyBlock,
    pub coordinates: CoordinateBlock,
    pub properties: MoleculeProperties,
}

#[derive(Debug, Clone, Copy)]
struct SectionOffsets {
    molecule_start: usize,
    atom_start: usize,
    bond_start: Option<usize>,
    charge_start: Option<usize>,
}

#[derive(Debug)]
struct MoleculeHeader {
    name: String,
    atom_count: u32,
    bond_count: u32,
    charge_type: String,
}

#[derive(Debug)]
struct ParsedAtom {
    spec: AtomSpec,
    position: [f64; 3],
}

#[derive(Debug, Default)]
struct DetachedBuilder {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    adjacency: Vec<Vec<BondId>>,
    coordinates: Vec<[f64; 3]>,
    properties: MoleculeProperties,
}

impl DetachedBuilder {
    fn add_atom(&mut self, spec: AtomSpec, position: [f64; 3]) -> AtomId {
        let id = AtomId::new(self.atoms.len());
        self.atoms.push(Atom::from_spec(id, spec));
        self.adjacency.push(Vec::new());
        self.coordinates.push(position);
        id
    }

    fn add_bond(&mut self, spec: BondSpec) -> Result<BondId, Mol2ReadError> {
        if spec.begin().index() >= self.atoms.len()
            || spec.end().index() >= self.atoms.len()
            || spec.begin() == spec.end()
        {
            return Err(Mol2ReadError::Parse("index mismatch".to_owned()));
        }
        let id = BondId::new(self.bonds.len());
        let begin = spec.begin();
        let end = spec.end();
        self.bonds.push(Bond::from_spec(id, spec));
        self.adjacency[begin.index()].push(id);
        self.adjacency[end.index()].push(id);
        Ok(id)
    }

    fn atoms(&self) -> &[Atom] {
        &self.atoms
    }

    fn atom_mut(&mut self, atom: AtomId) -> Option<&mut Atom> {
        self.atoms.get_mut(atom.index())
    }

    fn bonds(&self) -> &[Bond] {
        &self.bonds
    }

    fn bond(&self, bond: BondId) -> Option<&Bond> {
        self.bonds.get(bond.index())
    }

    fn bond_mut(&mut self, bond: BondId) -> Option<&mut Bond> {
        self.bonds.get_mut(bond.index())
    }

    fn neighbor_bonds(&self, atom: AtomId) -> &[BondId] {
        self.adjacency.get(atom.index()).map_or(&[], Vec::as_slice)
    }

    fn degree(&self, atom: AtomId) -> usize {
        self.neighbor_bonds(atom).len()
    }

    fn set_bond_order(&mut self, bond: BondId, order: BondOrder) -> Result<(), Mol2ReadError> {
        let bond = self
            .bond_mut(bond)
            .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_owned()))?;
        bond.set_order(order);
        Ok(())
    }

    fn finish(self) -> Result<Mol2Record, Mol2ReadError> {
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(self.atoms.len(), &self.bonds),
            atoms: self.atoms,
            bonds: self.bonds,
            ..TopologyBlock::default()
        };
        topology.validate()?;
        let coordinates = CoordinateBlock {
            conformers_2d: Vec::new(),
            conformers_3d: vec![Conformer3D::new(0, self.coordinates, true)],
            source_coordinate_dim: Some(CoordinateDimension::ThreeD),
        };
        coordinates.validate_for_atom_count(topology.atoms.len())?;
        Ok(Mol2Record {
            topology,
            coordinates,
            properties: self.properties,
        })
    }
}

fn get_line_at_with_eof<'a>(input: &'a str, offset: &mut usize) -> (&'a str, bool) {
    // RDKit source: RDGeneral/StreamOps.h `getLine`
    // RDKit✔️✔️: inline std::string getLine(std::istream *inStream) {
    // RDKit✔️✔️:   std::string res;
    // RDKit✔️✔️:   std::getline(*inStream, res);
    // RDKit✔️✔️:   if (!res.empty() && (res.back() == '\r')) {
    // RDKit✔️✔️:     res.resize(res.length() - 1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    if *offset >= input.len() {
        return ("", true);
    }
    let remaining = &input[*offset..];
    let (line, next, eof) = match remaining.find('\n') {
        Some(index) => (&remaining[..index], *offset + index + 1, false),
        None => (remaining, input.len(), true),
    };
    *offset = next;
    (line.strip_suffix('\r').unwrap_or(line), eof)
}

fn get_line_at<'a>(input: &'a str, offset: &mut usize) -> &'a str {
    get_line_at_with_eof(input, offset).0
}

fn get_checked_line_at<'a>(input: &'a str, offset: &mut usize) -> Result<&'a str, Mol2ReadError> {
    let (line, eof) = get_line_at_with_eof(input, offset);
    if eof {
        Err(Mol2ReadError::Parse("premature EOF".to_owned()))
    } else {
        Ok(line)
    }
}

fn parse_unsigned(token: &str) -> Option<u32> {
    let (negative, digits) = match token.as_bytes().first() {
        Some(b'+') => (false, &token[1..]),
        Some(b'-') => (true, &token[1..]),
        _ => (false, token),
    };
    if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
        return None;
    }
    let value = digits.parse::<u64>().ok()?;
    let value = u32::try_from(value).ok()?;
    Some(if negative {
        0_u32.wrapping_sub(value)
    } else {
        value
    })
}

fn scan_sections(input: &str) -> Result<SectionOffsets, Mol2ReadError> {
    // RDKit source: Mol2FileParser.cpp lines 824-857
    // RDKit✔️✔️:   std::streampos molStart = 0, atomStart = 0, bondStart = 0, chargeStart = 0;
    // RDKit✔️✔️:   while (!inStream.eof() && !inStream.fail()) {
    // RDKit✔️✔️:     tempStr = getLine(inStream);
    // RDKit✔️✔️:     if (inStream.eof()) {
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (tempStr != "" && tempStr[0] == '@') {
    // RDKit✔️✔️:       tokenizer tokens(tempStr, sep);
    // RDKit✔️✔️:       std::string firstToken = *tokens.begin();
    // RDKit✔️✔️:       if (firstToken == "@<TRIPOS>MOLECULE") {
    // RDKit✔️✔️:         if (!molStart) {
    // RDKit✔️✔️:           molStart = inStream.tellg();
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else if (firstToken == "@<TRIPOS>ATOM") {
    // RDKit✔️✔️:         atomStart = inStream.tellg();
    // RDKit✔️✔️:       } else if (firstToken == "@<TRIPOS>BOND") {
    // RDKit✔️✔️:         bondStart = inStream.tellg();
    // RDKit✔️✔️:       } else if (firstToken == "@<TRIPOS>UNITY_ATOM_ATTR") {
    // RDKit✔️✔️:         chargeStart = inStream.tellg();
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let mut molecule_start = None;
    let mut atom_start = None;
    let mut bond_start = None;
    let mut charge_start = None;
    let mut offset = 0;
    while offset < input.len() {
        let (line, eof) = get_line_at_with_eof(input, &mut offset);
        if eof {
            break;
        }
        if line.starts_with('@') {
            match line.split_whitespace().next().unwrap_or("") {
                "@<TRIPOS>MOLECULE" if molecule_start.is_none() => molecule_start = Some(offset),
                "@<TRIPOS>MOLECULE" => break,
                "@<TRIPOS>ATOM" => atom_start = Some(offset),
                "@<TRIPOS>BOND" => bond_start = Some(offset),
                "@<TRIPOS>UNITY_ATOM_ATTR" => charge_start = Some(offset),
                _ => {}
            }
        }
    }
    // RDKit✔️✔️:   if (!molStart) {
    // RDKit✔️✔️:     throw FileParseException("No MOLECULE block found in Mol2 data");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!atomStart) {
    // RDKit✔️✔️:     throw FileParseException("No ATOM block found in Mol2 data");
    // RDKit✔️✔️:   }
    Ok(SectionOffsets {
        molecule_start: molecule_start.ok_or_else(|| {
            Mol2ReadError::Parse("No MOLECULE block found in Mol2 data".to_owned())
        })?,
        atom_start: atom_start
            .ok_or_else(|| Mol2ReadError::Parse("No ATOM block found in Mol2 data".to_owned()))?,
        bond_start,
        charge_start,
    })
}

fn parse_header(input: &str, start: usize) -> Result<MoleculeHeader, Mol2ReadError> {
    // RDKit source: Mol2FileParser.cpp lines 863-902
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
    let mut offset = start;
    let name = get_line_at(input, &mut offset)
        .trim_end_matches(|character: char| character.is_ascii_whitespace())
        .to_owned();
    let counts = get_line_at(input, &mut offset);
    let mut tokens = counts.split_whitespace();
    let atom_token = tokens
        .next()
        .ok_or_else(|| Mol2ReadError::Parse("Empty counts line".to_owned()))?;
    let atom_count = parse_unsigned(atom_token).ok_or_else(|| {
        Mol2ReadError::Parse(format!("Cannot convert {atom_token} to unsigned int"))
    })?;
    let bond_count = if let Some(token) = tokens.next() {
        parse_unsigned(token).ok_or_else(|| {
            Mol2ReadError::Parse(format!("Cannot convert {token} to unsigned int"))
        })?
    } else {
        0
    };
    // RDKit✔️✔️:   if (nAtoms == 0) {
    // RDKit✔️✔️:     throw FileParseException("molecule has no atoms");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   tempStr = getLine(inStream);  // mol_type - ignore
    // RDKit✔️✔️:   tempStr = getLine(inStream);
    // RDKit✔️✔️:   boost::trim(tempStr);
    // RDKit✔️✔️:   res->setProp("_TriposChargeType", tempStr);
    if atom_count == 0 {
        return Err(Mol2ReadError::Parse("molecule has no atoms".to_owned()));
    }
    let _molecule_type = get_line_at(input, &mut offset);
    let charge_type = get_line_at(input, &mut offset).trim().to_owned();
    Ok(MoleculeHeader {
        name,
        atom_count,
        bond_count,
        charge_type,
    })
}

fn atom_spec_from_sybyl_symbol(symbol: &str) -> Result<Option<AtomSpec>, Mol2ReadError> {
    // RDKit source: Mol2FileParser.cpp lines 588-625
    // RDKit✔️✔️:   // bad symbols:
    // RDKit✔️✔️:   // LP is not an atom so remove it ...
    // RDKit✔️✔️:   if (symb == "LP") {
    // RDKit✔️✔️:     delete res;
    // RDKit✔️✔️:     return nullptr;
    if symbol == "LP" {
        return Ok(None);
    }
    // Query atoms cannot be installed in a concrete TopologyBlock. Preserve
    // the source capability boundary instead of replacing their predicates.
    // RDKit❌❌:   } else if (symb == "ANY" || symb == "Du") {
    // RDKit❌❌:     // queryAtoms
    // RDKit❌❌:     // according to the SYBYL spec, these match anything
    // RDKit❌❌:     auto *query = new QueryAtom(0);
    // RDKit❌❌:     query->setQuery(makeAtomNullQuery());
    // RDKit❌❌:   } else if (symb == "HEV") {
    // RDKit❌❌:     auto *query = new QueryAtom(1);
    // RDKit❌❌:     query->getQuery()->setNegation(true);
    // RDKit❌❌:   } else if (symb == "HET") {
    // RDKit❌❌:     // Tripos: N,O,P,S
    // RDKit❌❌:   } else if (symb == "HAL") {
    // RDKit❌❌:     // Tripos: F,Cl,Br,I
    if matches!(symbol, "ANY" | "Du" | "HEV" | "HET" | "HAL") {
        return Err(Mol2ReadError::Unsupported {
            feature: "MOL2 query atom types ANY/Du/HEV/HET/HAL in concrete TopologyBlock",
        });
    }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res->setAtomicNum(PeriodicTable::getTable()->getAtomicNumber(symb));
    // RDKit✔️✔️:   }
    let element = Element::from_symbol(symbol)
        .ok_or_else(|| Mol2ReadError::Parse(format!("Element '{symbol}' not found")))?;
    Ok(Some(AtomSpec::new(element)))
}

fn parse_atom_line(line: &str) -> Result<Option<ParsedAtom>, Mol2ReadError> {
    // RDKit source: Mol2FileParser.cpp lines 528-587, 626-643
    // RDKit✔️✔️: Atom *ParseMol2FileAtomLine(const std::string atomLine, RDGeom::Point3D &pos) {
    // RDKit✔️✔️:   typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
    // RDKit✔️✔️:   boost::char_separator<char> sep(" \t\n");
    // RDKit✔️✔️:   std::string tAN, tAT;
    // RDKit✔️✔️:   tokenizer tokens(atomLine, sep);
    // RDKit✔️✔️:   tokenizer::const_iterator itemIt = tokens.begin();
    // RDKit✔️✔️:   if (itemIt == tokens.end()) {
    // RDKit✔️✔️:     throw FileParseException("no info in mol2 atom line");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto *res = new Atom();
    // RDKit✔️✔️:   // skip TriposAtomId
    // RDKit✔️✔️:   ++itemIt;
    let mut tokens = line.split_whitespace();
    let _atom_id = tokens
        .next()
        .ok_or_else(|| Mol2ReadError::Parse("no info in mol2 atom line".to_owned()))?;
    let atom_name = tokens
        .next()
        .ok_or_else(|| Mol2ReadError::Parse("premature end of mol2 atom line".to_owned()))?;
    let coordinate = |token: Option<&str>| -> Result<f64, Mol2ReadError> {
        token
            .ok_or_else(|| Mol2ReadError::Parse("premature end of mol2 atom line".to_owned()))?
            .parse()
            .map_err(|_| Mol2ReadError::Parse("Cannot process mol2 coordinates.".to_owned()))
    };
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     pos.x = boost::lexical_cast<double>(*itemIt);
    // RDKit✔️✔️:     pos.y = boost::lexical_cast<double>(*itemIt);
    // RDKit✔️✔️:     pos.z = boost::lexical_cast<double>(*itemIt);
    // RDKit✔️✔️:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️✔️:     throw FileParseException("Cannot process mol2 coordinates.");
    // RDKit✔️✔️:   }
    let position = [
        coordinate(tokens.next())?,
        coordinate(tokens.next())?,
        coordinate(tokens.next())?,
    ];
    let atom_type = tokens
        .next()
        .ok_or_else(|| Mol2ReadError::Parse("premature end of mol2 atom line".to_owned()))?;
    // RDKit✔️✔️:   tAT = *itemIt;
    // RDKit✔️✔️:   std::string symb = (*itemIt).substr(0, (*itemIt).find('.'));
    // RDKit✔️✔️:   res->setProp("_TriposAtomName", tAN);
    // RDKit✔️✔️:   res->setProp(common_properties::_TriposAtomType, tAT);
    // RDKit✔️✔️:   // no implicit hydrogens for mol2 files
    // RDKit✔️✔️:   res->setNoImplicit(true);
    let symbol = atom_type
        .split_once('.')
        .map_or(atom_type, |(symbol, _)| symbol);
    let Some(spec) = atom_spec_from_sybyl_symbol(symbol)? else {
        return Ok(None);
    };
    let mut spec = spec
        .with_prop("_TriposAtomName", atom_name)?
        .with_prop("_TriposAtomType", atom_type)?
        .with_no_implicit(true);
    let _substructure_id = tokens.next();
    let _substructure_name = tokens.next();
    // RDKit✔️✔️:   // the Partial charge in the file
    // RDKit✔️✔️:   if (itemIt != tokens.end()) {
    // RDKit✔️✔️:     res->setProp("_TriposPartialCharge", *itemIt);
    // RDKit✔️✔️:   }
    if let Some(partial_charge) = tokens.next() {
        spec = spec.with_prop("_TriposPartialCharge", partial_charge)?;
    }
    Ok(Some(ParsedAtom { spec, position }))
}

fn parse_atom_block(
    input: &str,
    start: usize,
    count: u32,
) -> Result<(DetachedBuilder, Vec<Option<AtomId>>), Mol2ReadError> {
    // RDKit source: Mol2FileParser.cpp lines 723-776
    // RDKit✔️✔️:   std::vector<RDGeom::Point3D> threeDPs;
    // RDKit✔️✔️:   threeDPs.reserve(nAtoms);
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit✔️✔️:     std::string tempStr = getLine(inStream);
    // RDKit✔️✔️:     if (inStream->eof()) {
    // RDKit✔️✔️:       throw FileParseException("premature EOF");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     RDGeom::Point3D pos;
    // RDKit✔️✔️:     Atom *atom = ParseMol2FileAtomLine(tempStr, pos);
    // RDKit✔️✔️:     // if atom is NULL then we hit LP
    // RDKit✔️✔️:     if (atom) {
    // RDKit✔️✔️:       int aid = res->addAtom(atom, false, true);
    // RDKit✔️✔️:       idxCorresp[i] = aid;
    // RDKit✔️✔️:       threeDPs.push_back(pos);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       ++nLP;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let count = usize::try_from(count)
        .map_err(|_| Mol2ReadError::Parse("atom count does not fit platform usize".to_owned()))?;
    let mut builder = DetachedBuilder::default();
    let mut correspondence = vec![None; count];
    let mut offset = start;
    for slot in &mut correspondence {
        let line = get_checked_line_at(input, &mut offset)?;
        if let Some(parsed) = parse_atom_line(line)? {
            *slot = Some(builder.add_atom(parsed.spec, parsed.position));
        }
    }
    Ok((builder, correspondence))
}

fn parse_bond_line(
    line: &str,
    correspondence: &[Option<AtomId>],
) -> Result<Option<BondSpec>, Mol2ReadError> {
    // RDKit source: Mol2FileParser.cpp lines 646-720
    // RDKit✔️✔️:   tokenizer tokens(bondLine, sep);
    // RDKit✔️✔️:   tokenizer::const_iterator itemIt = tokens.begin();
    // RDKit✔️✔️:   if (itemIt == tokens.end()) {
    // RDKit✔️✔️:     throw FileParseException("no info in mol2 bond line");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     // tripos bond id skip
    // RDKit✔️✔️:     ++itemIt;
    // RDKit✔️✔️:     idx1 = boost::lexical_cast<unsigned int>(*itemIt);
    // RDKit✔️✔️:     idx2 = boost::lexical_cast<unsigned int>(*itemIt);
    // RDKit✔️✔️:   } catch (boost::bad_lexical_cast &) {
    // RDKit✔️✔️:     throw FileParseException("Cannot process mol2 bonds.");
    // RDKit✔️✔️:   }
    let mut tokens = line.split_whitespace();
    let _bond_id = tokens
        .next()
        .ok_or_else(|| Mol2ReadError::Parse("no info in mol2 bond line".to_owned()))?;
    let index = |token: Option<&str>| -> Result<usize, Mol2ReadError> {
        let token =
            token.ok_or_else(|| Mol2ReadError::Parse("no info in mol2 bond line".to_owned()))?;
        let index = parse_unsigned(token)
            .ok_or_else(|| Mol2ReadError::Parse("Cannot process mol2 bonds.".to_owned()))?
            .wrapping_sub(1);
        usize::try_from(index).map_err(|_| Mol2ReadError::Parse("index mismatch".to_owned()))
    };
    let begin_index = index(tokens.next())?;
    let end_index = index(tokens.next())?;
    let bond_type = tokens
        .next()
        .ok_or_else(|| Mol2ReadError::Parse("no info in mol2 bond line".to_owned()))?;
    // RDKit❗✔️:   if (!(idx1 < idxCorresp.size() || idx2 < idxCorresp.size())) {
    // RDKit❗✔️:     throw FileParseException("index mismatch");
    // RDKit❗✔️:   }
    let begin = correspondence
        .get(begin_index)
        .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_owned()))?;
    let end = correspondence
        .get(end_index)
        .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_owned()))?;
    let (Some(begin), Some(end)) = (*begin, *end) else {
        // RDKit✔️✔️:   if (idxCorresp[idx1] < 0 || idxCorresp[idx2] < 0) {
        // RDKit✔️✔️:     return nullptr;
        // RDKit✔️✔️:   }
        return Ok(None);
    };
    // RDKit✔️✔️:   if (tBType == "1" || tBType == "am") {
    // RDKit✔️✔️:     type = Bond::SINGLE;
    // RDKit✔️✔️:   } else if (tBType == "2") {
    // RDKit✔️✔️:     type = Bond::DOUBLE;
    // RDKit✔️✔️:   } else if (tBType == "3") {
    // RDKit✔️✔️:     type = Bond::TRIPLE;
    // RDKit✔️✔️:   } else if (tBType == "ar") {
    // RDKit✔️✔️:     type = Bond::AROMATIC;
    // RDKit✔️✔️:   } else if (tBType == "du" || tBType == "un") {
    // RDKit✔️✔️:     type = Bond::UNSPECIFIED;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return nullptr;
    // RDKit✔️✔️:   }
    let order = match bond_type {
        "1" | "am" => BondOrder::Single,
        "2" => BondOrder::Double,
        "3" => BondOrder::Triple,
        "ar" => BondOrder::Aromatic,
        "du" | "un" => BondOrder::Unspecified,
        _ => return Ok(None),
    };
    let mut spec = BondSpec::new(begin, end, order);
    if order == BondOrder::Aromatic {
        spec = spec.with_aromatic(true);
    }
    Ok(Some(spec))
}

fn parse_bond_block(
    input: &str,
    start: usize,
    count: u32,
    correspondence: &[Option<AtomId>],
    builder: &mut DetachedBuilder,
) -> Result<(), Mol2ReadError> {
    // RDKit source: Mol2FileParser.cpp lines 779-810
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
    let mut offset = start;
    for _ in 0..count {
        let line = get_checked_line_at(input, &mut offset)?;
        let Some(spec) = parse_bond_line(line, correspondence)? else {
            continue;
        };
        if spec.order() == BondOrder::Aromatic {
            builder
                .atom_mut(spec.begin())
                .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_owned()))?
                .set_aromatic(true);
            builder
                .atom_mut(spec.end())
                .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_owned()))?
                .set_aromatic(true);
        }
        builder.add_bond(spec)?;
    }
    Ok(())
}

fn read_unity_atom_attributes(
    input: &str,
    start: usize,
    builder: &mut DetachedBuilder,
) -> Result<(), Mol2ReadError> {
    // RDKit source: Mol2FileParser.cpp lines 77-137
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
    let mut offset = start;
    let mut line = get_checked_line_at(input, &mut offset)
        .map_err(|_| Mol2ReadError::Parse("premature EOF in readFormalCharges".to_owned()))?;
    loop {
        let mut tokens = line.split_whitespace();
        let atom_index = tokens
            .next()
            .and_then(parse_unsigned)
            .ok_or_else(|| Mol2ReadError::Parse("Cannot process mol2 UnityAtomAttr.".to_owned()))?;
        let attribute_count = tokens
            .next()
            .and_then(parse_unsigned)
            .ok_or_else(|| Mol2ReadError::Parse("Cannot process mol2 UnityAtomAttr.".to_owned()))?;
        for _ in 0..attribute_count {
            let attribute = get_checked_line_at(input, &mut offset).map_err(|_| {
                Mol2ReadError::Parse("premature EOF in readFormalCharges".to_owned())
            })?;
            let mut tokens = attribute.split_whitespace();
            if tokens.next() == Some("AtomExpr") {
                let charge = tokens.next().ok_or_else(|| {
                    Mol2ReadError::Parse("Cannot process mol2 formal charge.".to_owned())
                })?;
                // RDKit✔️✔️:         if ((*itemIt).find("=") == std::string::npos) {
                // RDKit✔️✔️:           formCharge = boost::lexical_cast<int>(*itemIt);
                // RDKit✔️✔️:           // assign the charge
                // RDKit✔️✔️:           res->getAtomWithIdx(atomIdx - 1)->setFormalCharge(formCharge);
                if !charge.contains('=') {
                    let charge = charge.parse::<i32>().map_err(|_| {
                        Mol2ReadError::Parse("Cannot process mol2 formal charge.".to_owned())
                    })?;
                    let charge = i8::try_from(charge).map_err(|_| {
                        Mol2ReadError::Parse("Cannot process mol2 formal charge.".to_owned())
                    })?;
                    let index = usize::try_from(atom_index.wrapping_sub(1))
                        .map_err(|_| Mol2ReadError::Parse("index mismatch".to_owned()))?;
                    builder
                        .atom_mut(AtomId::new(index))
                        .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_owned()))?
                        .set_formal_charge(charge);
                }
            }
        }
        if offset >= input.len() {
            break;
        }
        line = get_line_at(input, &mut offset);
        // RDKit✔️✔️:       if (tempStr == "" || tempStr[0] == '@' || tempStr[0] == '#') {
        // RDKit✔️✔️:         readNextAtomAttribs = false;
        // RDKit✔️✔️:       }
        if line.is_empty() || line.starts_with(['@', '#']) {
            break;
        }
    }
    Ok(())
}

fn fix_nitro_substructure(
    builder: &mut DetachedBuilder,
    atom: AtomId,
) -> Result<(), Mol2ReadError> {
    // RDKit source: Mol2FileParser.cpp lines 55-75
    // RDKit✔️✔️: void fixNitroSubstructureAndCharge(RWMol *res, unsigned int atIdx) {
    // RDKit✔️✔️:   unsigned int noODblNeighbors = 0;
    // RDKit✔️✔️:   unsigned int toModIdx = 0;
    // RDKit✔️✔️:   while (nbrIdxIt != nbrEndIdxIt) {
    // RDKit✔️✔️:     Bond *curBond = res->getBondBetweenAtoms(atIdx, *nbrIdxIt);
    // RDKit✔️✔️:     if (res->getAtomWithIdx(*nbrIdxIt)->getAtomicNum() == 8 &&
    // RDKit✔️✔️:         curBond->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:       ++noODblNeighbors;
    // RDKit✔️✔️:       toModIdx = *nbrIdxIt;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (noODblNeighbors == 2) {
    // RDKit✔️✔️:     res->getBondBetweenAtoms(atIdx, toModIdx)->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:     res->getAtomWithIdx(atIdx)->setFormalCharge(1);
    // RDKit✔️✔️:     res->getAtomWithIdx(toModIdx)->setFormalCharge(-1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let mut matches = Vec::new();
    for bond_id in builder.neighbor_bonds(atom) {
        let bond = builder
            .bond(*bond_id)
            .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_owned()))?;
        let neighbor = if bond.begin() == atom {
            bond.end()
        } else {
            bond.begin()
        };
        if builder.atoms()[neighbor.index()].atomic_number() == 8
            && bond.order() == BondOrder::Double
        {
            matches.push((neighbor, *bond_id));
        }
    }
    if matches.len() == 2 {
        let (oxygen, bond) = matches[1];
        builder.set_bond_order(bond, BondOrder::Single)?;
        builder.atom_mut(atom).unwrap().set_formal_charge(1);
        builder.atom_mut(oxygen).unwrap().set_formal_charge(-1);
    }
    Ok(())
}

fn guess_formal_charges(builder: &mut DetachedBuilder) -> Result<(), Mol2ReadError> {
    // RDKit source: Mol2FileParser.cpp lines 139-264
    // RDKit✔️✔️: void guessFormalCharges(RWMol *res) {
    // RDKit✔️✔️:   // FIX: this whole thing has problems with positively charged pyridines et al.
    // RDKit✔️✔️:   for (RWMol::AtomIterator atomIt = res->beginAtoms();
    // RDKit✔️✔️:        atomIt != res->endAtoms(); ++atomIt) {
    // RDKit✔️✔️:     Atom *at = (*atomIt);
    // RDKit✔️✔️:     if (at->getFormalCharge() == 0 && at->getSymbol() != "C" &&
    // RDKit✔️✔️:         !(at->hasQuery())) {
    let mut ring_info = None;
    for index in 0..builder.atoms.len() {
        let atom_id = AtomId::new(index);
        if builder.atoms[index].formal_charge() != 0 || builder.atoms[index].atomic_number() == 6 {
            continue;
        }
        // RDKit✔️✔️:       int noAromBonds = 0;
        // RDKit✔️✔️:       double accum = 0;
        // RDKit✔️✔️:       for (const auto bnd : res->atomBonds(at)) {
        // RDKit✔️✔️:         accum += bnd->getValenceContrib(at);
        // RDKit✔️✔️:         if (bnd->getBondType() == Bond::AROMATIC) {
        // RDKit✔️✔️:           ++noAromBonds;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        let mut aromatic_bonds = 0;
        let mut accumulated = 0.0;
        for bond_id in builder.neighbor_bonds(atom_id) {
            let bond = builder.bond(*bond_id).unwrap();
            accumulated += bond_valence_contrib(bond, atom_id)
                .map_err(|error| Mol2ReadError::Parse(error.to_string()))?;
            aromatic_bonds += i32::from(bond.order() == BondOrder::Aromatic);
        }
        if ring_info.is_none() {
            let adjacency = AdjacencyList::from_topology(builder.atoms.len(), &builder.bonds);
            ring_info = Some(
                find_sssr_from_parts(builder.atoms.len(), &builder.bonds, &adjacency)
                    .map_err(|error| Mol2ReadError::Parse(error.to_string()))?,
            );
        }
        let tripos_type = builder.atoms[index]
            .prop("_TriposAtomType")
            .ok_or_else(|| Mol2ReadError::Parse("Missing _TriposAtomType".to_owned()))?;
        // RDKit✔️✔️:       if (at->getIsAromatic() && tATT.find("ar") == std::string::npos &&
        // RDKit✔️✔️:           res->getRingInfo()->isAtomInRingOfSize(at->getIdx(), 5)) {
        // RDKit✔️✔️:         continue;
        // RDKit✔️✔️:       }
        if builder.atoms[index].is_aromatic()
            && !tripos_type.contains("ar")
            && ring_info
                .as_ref()
                .is_some_and(|rings| rings.is_atom_in_ring_of_size(atom_id, 5))
        {
            continue;
        }
        // RDKit✔️✔️:       if (noAromBonds == 3 && tATT == "N.ar") {
        // RDKit✔️✔️:         continue;
        // RDKit✔️✔️:       }
        if aromatic_bonds == 3 && tripos_type == "N.ar" {
            continue;
        }
        // RDKit✔️✔️:       auto expVal = static_cast<int>(std::round(accum + 0.1));
        // RDKit✔️✔️:       const auto &valens =
        // RDKit✔️✔️:           PeriodicTable::getTable()->getValenceList(at->getAtomicNum());
        // RDKit✔️✔️:       int nElectrons =
        // RDKit✔️✔️:           PeriodicTable::getTable()->getNouterElecs(at->getAtomicNum());
        let explicit_valence = (accumulated + 0.1).round() as i32;
        let valences = rdkit_valence_list(builder.atoms[index].atomic_number())
            .map_err(|error| Mol2ReadError::Parse(error.to_string()))?
            .ok_or_else(|| Mol2ReadError::Parse("Missing RDKit valence list".to_owned()))?;
        let outer_electrons = periodic_table_outer_electrons(builder.atoms[index].atomic_number())
            .map_err(|error| Mol2ReadError::Parse(error.to_string()))?;
        let mut charge = if outer_electrons >= 4 {
            explicit_valence - valences[0]
        } else {
            valences[0] - explicit_valence
        };
        // RDKit✔️✔️:       if (assignChg > 0 && nElectrons >= 4) {
        // RDKit✔️✔️:         for (auto vi : valens) {
        // RDKit✔️✔️:           assignChg = expVal - vi;
        // RDKit✔️✔️:           if (vi <= expVal && assignChg < 2) {
        // RDKit✔️✔️:             break;
        // RDKit✔️✔️:           }
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        if charge > 0 && outer_electrons >= 4 {
            for &valence in valences {
                charge = explicit_valence - valence;
                if valence <= explicit_valence && charge < 2 {
                    break;
                }
            }
        }
        if charge != 0 {
            // RDKit✔️✔️:         if (at->getIsAromatic() && abs(assignChg) > 1) {
            // RDKit✔️✔️:           at->setFormalCharge((assignChg > 0) -
            // RDKit✔️✔️:                               (assignChg < 0));
            // RDKit✔️✔️:         } else {
            // RDKit✔️✔️:           at->setFormalCharge(assignChg);
            // RDKit✔️✔️:         }
            let assigned = if builder.atoms[index].is_aromatic() && charge.unsigned_abs() > 1 {
                i8::from(charge > 0) - i8::from(charge < 0)
            } else {
                i8::try_from(charge)
                    .map_err(|_| Mol2ReadError::Parse("formal charge out of range".to_owned()))?
            };
            builder.atoms[index].set_formal_charge(assigned);
            // RDKit✔️✔️:         if (assignChg == 2 && expVal == 5 && at->getSymbol() == "N") {
            // RDKit✔️✔️:           fixNitroSubstructureAndCharge(res, at->getIdx());
            // RDKit✔️✔️:         }
            if charge == 2 && explicit_valence == 5 && builder.atoms[index].atomic_number() == 7 {
                fix_nitro_substructure(builder, atom_id)?;
            }
        }
    }
    Ok(())
}

fn check_no_h_neighbors_n_oxide(
    builder: &DetachedBuilder,
    atom: AtomId,
    atom_to_modify: &mut Option<AtomId>,
) -> Result<u32, Mol2ReadError> {
    // RDKit source: Mol2FileParser.cpp lines 267-287
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
    if atom.index() >= builder.atoms().len() {
        return Err(Mol2ReadError::Parse("index mismatch".to_owned()));
    }
    let mut hydrogen_neighbors = 0_u32;
    for bond_id in builder.neighbor_bonds(atom).iter().copied() {
        let bond = builder
            .bond(bond_id)
            .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_owned()))?;
        let neighbor = if bond.begin() == atom {
            bond.end()
        } else {
            bond.begin()
        };
        let neighbor_atom = builder
            .atoms()
            .get(neighbor.index())
            .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_owned()))?;
        if neighbor_atom.atomic_number() == 1 {
            hydrogen_neighbors = hydrogen_neighbors.wrapping_add(1);
        } else if neighbor_atom.atomic_number() == 8 && builder.degree(neighbor) == 1 {
            *atom_to_modify = Some(atom);
        }
    }
    Ok(hydrogen_neighbors)
}

fn cleanup_substructures(builder: &mut DetachedBuilder) -> Result<bool, Mol2ReadError> {
    // RDKit source: Mol2FileParser.cpp lines 289-523
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
    // RDKit✔️✔️:       // for the structural fix - the last nitrogen with the least number of
    // RDKit✔️✔️:       // heavy atoms will get the double bond and the positive charge.
    // RDKit✔️✔️:       // remember : this is not canonical!
    // RDKit✔️✔️:       // first - set the C.cat as fixed
    // RDKit✔️✔️:       isFixed[idx] = 1;
    // RDKit✔️✔️:       ROMol::ADJ_ITER nbrIdxIt, endNbrsIdxIt, tmpIdxIt;
    // RDKit✔️✔️:       unsigned int lowestDeg = 100;
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
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     idx++;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    let mut fixed = vec![false; builder.atoms().len()];
    for index in 0..builder.atoms().len() {
        if fixed[index] {
            continue;
        }
        let atom = AtomId::new(index);
        let tripos_type = builder.atoms()[index]
            .prop("_TriposAtomType")
            .ok_or_else(|| Mol2ReadError::Parse("Missing _TriposAtomType".to_owned()))?
            .to_owned();
        if tripos_type == "N.4" {
            builder
                .atom_mut(atom)
                .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_owned()))?
                .set_formal_charge(1);
        } else if tripos_type == "O.co2" {
            if builder.degree(atom) != 1 {
                return Ok(false);
            }
            let bond_id = builder.neighbor_bonds(atom)[0];
            let bond = builder
                .bond(bond_id)
                .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_owned()))?;
            let neighbor = if bond.begin() == atom {
                bond.end()
            } else {
                bond.begin()
            };
            let neighbor_type = builder.atoms()[neighbor.index()]
                .prop("_TriposAtomType")
                .ok_or_else(|| Mol2ReadError::Parse("Missing _TriposAtomType".to_owned()))?
                .to_owned();
            if neighbor_type == "P.3" {
                let bond = builder
                    .bond_mut(bond_id)
                    .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_owned()))?;
                bond.set_order(BondOrder::Double);
                bond.set_aromatic(false);
                builder.atom_mut(atom).unwrap().set_aromatic(false);
                fixed[index] = true;
                let oxygen_bonds = builder.neighbor_bonds(neighbor).to_vec();
                for oxygen_bond_id in oxygen_bonds {
                    let oxygen_bond = builder.bond(oxygen_bond_id).unwrap();
                    let oxygen = if oxygen_bond.begin() == neighbor {
                        oxygen_bond.end()
                    } else {
                        oxygen_bond.begin()
                    };
                    if builder.atoms()[oxygen.index()].atomic_number() == 8
                        && !fixed[oxygen.index()]
                        && builder.atoms()[oxygen.index()].prop("_TriposAtomType") == Some("O.co2")
                    {
                        let oxygen_bond = builder.bond_mut(oxygen_bond_id).unwrap();
                        oxygen_bond.set_order(BondOrder::Single);
                        oxygen_bond.set_aromatic(false);
                        let oxygen_atom = builder.atom_mut(oxygen).unwrap();
                        oxygen_atom.set_formal_charge(-1);
                        oxygen_atom.set_aromatic(false);
                        fixed[oxygen.index()] = true;
                    }
                }
                builder.atom_mut(neighbor).unwrap().set_aromatic(false);
                fixed[neighbor.index()] = true;
            } else if neighbor_type == "C.2" || neighbor_type == "S.o2" {
                if !fixed[neighbor.index()] {
                    let bond = builder.bond_mut(bond_id).unwrap();
                    bond.set_order(BondOrder::Single);
                    bond.set_aromatic(false);
                    let oxygen = builder.atom_mut(atom).unwrap();
                    oxygen.set_formal_charge(-1);
                    oxygen.set_aromatic(false);
                    builder.atom_mut(neighbor).unwrap().set_aromatic(false);
                    fixed[index] = true;
                    fixed[neighbor.index()] = true;
                } else {
                    let bond = builder.bond_mut(bond_id).unwrap();
                    bond.set_order(BondOrder::Double);
                    bond.set_aromatic(false);
                    builder.atom_mut(atom).unwrap().set_aromatic(false);
                    fixed[index] = true;
                }
            } else {
                return Ok(false);
            }
        } else if tripos_type == "C.cat" {
            fixed[index] = true;
            let neighbor_pairs = builder
                .neighbor_bonds(atom)
                .iter()
                .copied()
                .map(|bond_id| {
                    let bond = builder.bond(bond_id).unwrap();
                    let neighbor = if bond.begin() == atom {
                        bond.end()
                    } else {
                        bond.begin()
                    };
                    (neighbor, bond_id)
                })
                .collect::<Vec<_>>();
            let nitrogen_neighbors = neighbor_pairs
                .iter()
                .copied()
                .filter(|(neighbor, _)| builder.atoms()[neighbor.index()].atomic_number() == 7)
                .collect::<Vec<_>>();
            if !(2..=3).contains(&nitrogen_neighbors.len()) {
                return Ok(false);
            }
            if nitrogen_neighbors.len() == 2 {
                let (first_nitrogen, first_bond) = nitrogen_neighbors[0];
                let (second_nitrogen, second_bond) = nitrogen_neighbors[1];
                for (nitrogen, bond_id) in nitrogen_neighbors.iter().copied() {
                    let bond = builder.bond_mut(bond_id).unwrap();
                    bond.set_order(BondOrder::Single);
                    bond.set_aromatic(false);
                    builder.atom_mut(nitrogen).unwrap().set_aromatic(false);
                    fixed[nitrogen.index()] = true;
                }
                let mut atom_to_modify = None;
                let first_hydrogens =
                    check_no_h_neighbors_n_oxide(builder, first_nitrogen, &mut atom_to_modify)?;
                let second_hydrogens =
                    check_no_h_neighbors_n_oxide(builder, second_nitrogen, &mut atom_to_modify)?;
                let atom_to_modify = if let Some(atom_to_modify) = atom_to_modify {
                    atom_to_modify
                } else if first_hydrogens != second_hydrogens {
                    if first_hydrogens > second_hydrogens {
                        first_nitrogen
                    } else {
                        second_nitrogen
                    }
                } else {
                    let adjacency =
                        AdjacencyList::try_from_topology(builder.atoms().len(), builder.bonds())
                            .map_err(|error| Mol2ReadError::Parse(error.to_string()))?;
                    let rings =
                        find_sssr_from_parts(builder.atoms().len(), builder.bonds(), &adjacency)
                            .map_err(|error| Mol2ReadError::Parse(error.to_string()))?;
                    if rings.num_atom_rings(first_nitrogen) > rings.num_atom_rings(second_nitrogen)
                    {
                        first_nitrogen
                    } else {
                        second_nitrogen
                    }
                };
                let bond_to_modify = if atom_to_modify == first_nitrogen {
                    first_bond
                } else {
                    second_bond
                };
                let bond = builder.bond_mut(bond_to_modify).unwrap();
                bond.set_order(BondOrder::Double);
                bond.set_aromatic(false);
                builder
                    .atom_mut(atom_to_modify)
                    .unwrap()
                    .set_formal_charge(1);
                builder.atom_mut(atom).unwrap().set_aromatic(false);
            } else {
                let mut lowest_degree = 100_u32;
                let mut atom_to_modify = None;
                let mut bond_to_modify = None;
                for (neighbor, bond_id) in neighbor_pairs {
                    if !fixed[neighbor.index()] {
                        let mut heavy_atom_degree = 0_u32;
                        for neighbor_bond_id in builder.neighbor_bonds(neighbor).iter().copied() {
                            let neighbor_bond = builder.bond(neighbor_bond_id).unwrap();
                            let next_neighbor = if neighbor_bond.begin() == neighbor {
                                neighbor_bond.end()
                            } else {
                                neighbor_bond.begin()
                            };
                            if builder.atoms()[next_neighbor.index()].atomic_number() > 1 {
                                if builder.atoms()[next_neighbor.index()].prop("_TriposAtomType")
                                    == Some("C.cat")
                                {
                                    heavy_atom_degree = heavy_atom_degree.wrapping_add(2);
                                } else {
                                    heavy_atom_degree = heavy_atom_degree.wrapping_add(1);
                                }
                            }
                        }
                        if heavy_atom_degree < lowest_degree {
                            atom_to_modify = Some(neighbor);
                            bond_to_modify = Some(bond_id);
                            lowest_degree = heavy_atom_degree;
                        }
                        let bond = builder.bond_mut(bond_id).unwrap();
                        bond.set_order(BondOrder::Single);
                        bond.set_aromatic(false);
                        builder.atom_mut(neighbor).unwrap().set_aromatic(false);
                        fixed[neighbor.index()] = true;
                    } else {
                        let bond = builder.bond_mut(bond_id).unwrap();
                        bond.set_order(BondOrder::Single);
                        bond.set_aromatic(false);
                    }
                }
                let atom_to_modify = atom_to_modify
                    .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_owned()))?;
                let bond_to_modify = bond_to_modify
                    .ok_or_else(|| Mol2ReadError::Parse("index mismatch".to_owned()))?;
                let bond = builder.bond_mut(bond_to_modify).unwrap();
                bond.set_order(BondOrder::Double);
                bond.set_aromatic(false);
                builder
                    .atom_mut(atom_to_modify)
                    .unwrap()
                    .set_formal_charge(1);
                builder.atom_mut(atom).unwrap().set_aromatic(false);
            }
        }
    }
    Ok(true)
}

/// Parse a MOL2 block into detached values and apply MOL2-specific charge rules.
pub fn read_mol2_detached_with_params(
    input: &str,
    params: Mol2ReadParams,
) -> Result<Option<Mol2Record>, Mol2ReadError> {
    let offsets = scan_sections(input)?;
    let header = parse_header(input, offsets.molecule_start)?;
    let (mut builder, correspondence) =
        parse_atom_block(input, offsets.atom_start, header.atom_count)?;
    if header.bond_count > 0 {
        let start = offsets
            .bond_start
            .ok_or_else(|| Mol2ReadError::Parse("No BOND block found".to_owned()))?;
        parse_bond_block(
            input,
            start,
            header.bond_count,
            &correspondence,
            &mut builder,
        )?;
    }
    builder.properties = MoleculeProperties::default()
        .with_name(header.name)
        .with_prop("_TriposChargeType", header.charge_type)?;
    if let Some(start) = offsets.charge_start {
        read_unity_atom_attributes(input, start, &mut builder)?;
    } else {
        if params.cleanup_substructures && !cleanup_substructures(&mut builder)? {
            return Ok(None);
        }
        guess_formal_charges(&mut builder)?;
    }
    builder.finish().map(Some)
}

/// Parse a MOL2 block with the source defaults.
pub fn read_mol2_detached(input: &str) -> Result<Option<Mol2Record>, Mol2ReadError> {
    read_mol2_detached_with_params(input, Mol2ReadParams::default())
}

#[cfg(test)]
mod tests {
    use super::*;

    const BASIC: &str = "@<TRIPOS>MOLECULE\nexample   \n3 2\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0.0 0.0 0.0 C.2 1 MOL 0.25\n2 O1 1.2 0.0 0.0 O.2 1 MOL -0.25\n3 H1 -0.5 0.0 0.0 H 1 MOL 0.0\n@<TRIPOS>BOND\n1 1 2 2\n2 1 3 1\n";

    fn rdkit_fixture(name: &str) -> &'static str {
        match name {
            "3505.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/3505.mol2"
            ),
            "Canion.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/Canion.mol2"
            ),
            "EZ_mol2_issue114.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/EZ_mol2_issue114.mol2"
            ),
            "Issue3399798.2.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/Issue3399798.2.mol2"
            ),
            "Issue3399798.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/Issue3399798.mol2"
            ),
            "Noxide.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/Noxide.mol2"
            ),
            "Sulfonate.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/Sulfonate.mol2"
            ),
            "badSubstPyridine.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/badSubstPyridine.mol2"
            ),
            "benzene.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/benzene.mol2"
            ),
            "chargedAmidine.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/chargedAmidine.mol2"
            ),
            "chargedAmidineEC.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/chargedAmidineEC.mol2"
            ),
            "chargedAmidineRWH.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/chargedAmidineRWH.mol2"
            ),
            "dbtranslateCharged.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/dbtranslateCharged.mol2"
            ),
            "dbtranslateUncharged.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/dbtranslateUncharged.mol2"
            ),
            "dbtranslateUnchargedRing.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/dbtranslateUnchargedRing.mol2"
            ),
            "fusedRing.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/fusedRing.mol2"
            ),
            "github438_1.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/github438_1.mol2"
            ),
            "github438_2.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/github438_2.mol2"
            ),
            "highlySymmetricGuanidine.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/highlySymmetricGuanidine.mol2"
            ),
            "lonePairMol.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/lonePairMol.mol2"
            ),
            "pyrazole_pyridine.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/pyrazole_pyridine.mol2"
            ),
            "pyridiniumPhenyl.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/pyridiniumPhenyl.mol2"
            ),
            "sulfonAmide.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/sulfonAmide.mol2"
            ),
            "symmetricGuanidine.mol2" => include_str!(
                "../../../third_party/rdkit/Code/GraphMol/FileParsers/test_data/symmetricGuanidine.mol2"
            ),
            _ => panic!("unknown fixture {name}"),
        }
    }

    #[test]
    fn parses_concrete_atom_bond_coordinate_and_property_state() {
        let record = read_mol2_detached_with_params(
            BASIC,
            Mol2ReadParams {
                cleanup_substructures: false,
                ..Mol2ReadParams::default()
            },
        )
        .expect("MOL2")
        .expect("record");
        assert_eq!(record.topology.atoms.len(), 3);
        assert_eq!(record.topology.bonds.len(), 2);
        assert_eq!(record.topology.bonds[0].order(), BondOrder::Double);
        assert_eq!(record.properties.name(), Some("example"));
        assert_eq!(
            record.properties.prop("_TriposChargeType"),
            Some("NO_CHARGES")
        );
        assert_eq!(record.topology.atoms[0].prop("_TriposAtomName"), Some("C1"));
        assert_eq!(
            record.topology.atoms[0].prop("_TriposPartialCharge"),
            Some("0.25")
        );
        assert_eq!(
            record.coordinates.conformers_3d[0].coordinates()[1],
            [1.2, 0.0, 0.0]
        );
        assert!(record.topology.atoms.iter().all(Atom::no_implicit));
    }

    #[test]
    fn skips_lone_pairs_and_bonds_to_them() {
        let input = "@<TRIPOS>MOLECULE\nlp\n3 2\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n2 LP1 1 0 0 LP\n3 O1 2 0 0 O.2\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 2\n@<TRIPOS>UNITY_ATOM_ATTR\n1 1\nAtomExpr 0\n";
        let record = read_mol2_detached(input).expect("MOL2").expect("record");
        assert_eq!(record.topology.atoms.len(), 2);
        assert_eq!(record.topology.bonds.len(), 1);
        assert_eq!(record.topology.bonds[0].end(), AtomId::new(1));
    }

    #[test]
    fn unity_atom_attributes_override_formal_charge_guessing() {
        let input = "@<TRIPOS>MOLECULE\ncharged\n1 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 N1 0 0 0 N.4\n@<TRIPOS>UNITY_ATOM_ATTR\n1 1\nAtomExpr -1\n";
        let record = read_mol2_detached(input).expect("MOL2").expect("record");
        assert_eq!(record.topology.atoms[0].formal_charge(), -1);
    }

    #[test]
    fn aromatic_bonds_mark_bond_and_atoms() {
        let input = "@<TRIPOS>MOLECULE\naromatic\n2 1\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.ar\n2 N1 1 0 0 N.ar\n@<TRIPOS>BOND\n1 1 2 ar\n@<TRIPOS>UNITY_ATOM_ATTR\n1 1\nAtomExpr 0\n";
        let record = read_mol2_detached(input).expect("MOL2").expect("record");
        assert!(record.topology.bonds[0].is_aromatic());
        assert!(record.topology.atoms[0].is_aromatic());
        assert!(record.topology.atoms[1].is_aromatic());
    }

    #[test]
    fn query_atom_types_fail_closed_at_concrete_boundary() {
        let input =
            "@<TRIPOS>MOLECULE\nquery\n1 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 A1 0 0 0 ANY\n";
        assert!(matches!(
            read_mol2_detached_with_params(
                input,
                Mol2ReadParams {
                    cleanup_substructures: false,
                    ..Mol2ReadParams::default()
                }
            ),
            Err(Mol2ReadError::Unsupported { .. })
        ));
    }

    #[test]
    fn source_errors_cover_missing_sections_zero_atoms_and_premature_eof() {
        assert!(
            read_mol2_detached("")
                .unwrap_err()
                .to_string()
                .contains("No MOLECULE")
        );
        assert!(
            read_mol2_detached("@<TRIPOS>MOLECULE\nnone\n0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n")
                .unwrap_err()
                .to_string()
                .contains("no atoms")
        );
        assert!(
            read_mol2_detached_with_params(
                "@<TRIPOS>MOLECULE\nshort\n1 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n",
                Mol2ReadParams {
                    cleanup_substructures: false,
                    ..Mol2ReadParams::default()
                }
            )
            .unwrap_err()
            .to_string()
            .contains("premature EOF")
        );
        assert!(
            scan_sections("@<TRIPOS>MOLECULE\nmol\n1 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM")
                .unwrap_err()
                .to_string()
                .contains("No ATOM"),
            "RDKit stops before interpreting a final unterminated section line"
        );
    }

    #[test]
    fn rdkit_mol2_corpus_matches_unsanitized_graph_charge_and_aromatic_state() {
        // RDKit source: testMol2ToMol.cpp `testGeneral`, `testGithub438`
        // RDKit✔️✔️: TEST_ASSERT(m->getAtomWithIdx(8)->getFormalCharge() == 1);
        // RDKit✔️✔️: TEST_ASSERT(m->getAtomWithIdx(9)->getFormalCharge() == -1);
        // RDKit✔️✔️: TEST_ASSERT(mol->getAtomWithIdx(0)->getFormalCharge() == 1);
        // RDKit✔️✔️: TEST_ASSERT(mol->getAtomWithIdx(0)->getFormalCharge() == 2);
        // Expected values below were also checked against Python RDKit with
        // sanitize=false, removeHs=false, and cleanupSubstructures=true so
        // indices cover the detached parser's pre-finalization boundary.
        let cases: &[(&str, usize, usize, usize, usize, &[(usize, i8)])] = &[
            ("3505.mol2", 34, 33, 0, 0, &[(9, 1)]),
            ("Canion.mol2", 16, 15, 0, 0, &[]),
            ("EZ_mol2_issue114.mol2", 9, 8, 0, 0, &[]),
            ("Issue3399798.2.mol2", 50, 52, 12, 12, &[]),
            ("Issue3399798.mol2", 50, 52, 12, 12, &[]),
            ("Noxide.mol2", 27, 27, 6, 6, &[(8, 1), (9, -1)]),
            ("Sulfonate.mol2", 8, 7, 0, 0, &[(2, -1)]),
            ("badSubstPyridine.mol2", 15, 15, 6, 6, &[(5, 1)]),
            ("benzene.mol2", 12, 12, 6, 6, &[]),
            ("chargedAmidine.mol2", 22, 23, 6, 6, &[(9, 1)]),
            ("chargedAmidineEC.mol2", 13, 13, 0, 0, &[(3, 1)]),
            ("chargedAmidineRWH.mol2", 19, 20, 6, 6, &[(6, 1)]),
            ("dbtranslateCharged.mol2", 18, 18, 6, 6, &[(17, 1)]),
            ("dbtranslateUncharged.mol2", 20, 20, 6, 6, &[]),
            ("dbtranslateUnchargedRing.mol2", 8, 8, 5, 5, &[]),
            ("fusedRing.mol2", 22, 24, 14, 16, &[]),
            ("github438_1.mol2", 6, 4, 0, 0, &[(0, 1)]),
            ("github438_2.mol2", 6, 4, 0, 0, &[(0, 2)]),
            (
                "highlySymmetricGuanidine.mol2",
                65,
                64,
                0,
                0,
                &[(4, 1), (13, 1), (22, 1)],
            ),
            ("lonePairMol.mol2", 11, 10, 0, 0, &[]),
            ("pyrazole_pyridine.mol2", 9, 9, 0, 0, &[]),
            ("pyridiniumPhenyl.mol2", 22, 23, 12, 12, &[(5, 1)]),
            ("sulfonAmide.mol2", 10, 9, 0, 0, &[]),
            ("symmetricGuanidine.mol2", 28, 27, 0, 0, &[(1, 1), (8, 1)]),
        ];
        for &(name, atom_count, bond_count, aromatic_atoms, aromatic_bonds, charges) in cases {
            let record = read_mol2_detached(rdkit_fixture(name))
                .unwrap_or_else(|error| panic!("{name}: {error}"))
                .unwrap_or_else(|| panic!("{name}: rejected cleanup"));
            assert_eq!(record.topology.atoms.len(), atom_count, "{name}: atoms");
            assert_eq!(record.topology.bonds.len(), bond_count, "{name}: bonds");
            assert_eq!(
                record
                    .topology
                    .atoms
                    .iter()
                    .filter(|atom| atom.is_aromatic())
                    .count(),
                aromatic_atoms,
                "{name}: aromatic atoms"
            );
            assert_eq!(
                record
                    .topology
                    .bonds
                    .iter()
                    .filter(|bond| bond.is_aromatic())
                    .count(),
                aromatic_bonds,
                "{name}: aromatic bonds"
            );
            let actual_charges = record
                .topology
                .atoms
                .iter()
                .enumerate()
                .filter_map(|(index, atom)| {
                    (atom.formal_charge() != 0).then_some((index, atom.formal_charge()))
                })
                .collect::<Vec<_>>();
            assert_eq!(actual_charges, charges, "{name}: formal charges");
        }
    }

    #[test]
    fn cleanup_parameter_matches_rdkit_guanidinium_branch() {
        // RDKit source: testMol2ToMol.cpp `testDisableCleanup`
        // RDKit✔️✔️: TEST_ASSERT(mol->getBondBetweenAtoms(3, 12)->getBondType() ==
        // RDKit✔️✔️:             Bond::SINGLE);
        // RDKit✔️✔️: TEST_ASSERT(mol->getAtomWithIdx(12)->getFormalCharge() == 0);
        // RDKit✔️✔️: TEST_ASSERT(mol->getBondBetweenAtoms(3, 12)->getBondType() ==
        // RDKit✔️✔️:             Bond::DOUBLE);
        // RDKit✔️✔️: TEST_ASSERT(mol->getAtomWithIdx(12)->getFormalCharge() == 1);
        let input = rdkit_fixture("3505.mol2");
        let cleaned = read_mol2_detached(input).unwrap().unwrap();
        let raw = read_mol2_detached_with_params(
            input,
            Mol2ReadParams {
                cleanup_substructures: false,
                ..Mol2ReadParams::default()
            },
        )
        .unwrap()
        .unwrap();
        let matching_bond = |record: &Mol2Record| {
            record
                .topology
                .bonds
                .iter()
                .find(|bond| {
                    [bond.begin().index(), bond.end().index()].contains(&3)
                        && [bond.begin().index(), bond.end().index()].contains(&12)
                })
                .unwrap()
                .order()
        };
        assert_eq!(matching_bond(&cleaned), BondOrder::Single);
        assert_eq!(cleaned.topology.atoms[12].formal_charge(), 0);
        assert_eq!(matching_bond(&raw), BondOrder::Double);
        assert_eq!(raw.topology.atoms[12].formal_charge(), 1);
    }
}
