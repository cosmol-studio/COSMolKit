// BEGIN RDKIT CPP FUNCTION: MolPickler::pickleMol (MolPickler.cpp)
// RDKit✔️✔️: void MolPickler::pickleMol(const ROMol &mol, std::ostream &ss) {
// RDKit✔️✔️:   // Serializes molecule state: atoms, bonds, coordinates,
// RDKit✔️✔️:   // substance groups, stereo groups, and properties.
// RDKit✔️✔️:   MolPickler::pickleMol(mol, ss, PicklerOps::defaultOps());
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION: MolPickler::pickleMol
//
// COSMolKit MolPickler — compact binary serialization for Molecule.
// Uses a forward-compatible versioned format with little-endian encoding.

use std::collections::BTreeMap;

use crate::{
    Atom, AtomId, Bond, BondDirection, BondId, BondOrder, BondStereo, ChiralTag, Conformer3D,
    CoordinateDimension, Hybridization, Molecule, SGroupAttachPoint, SGroupBondRole, SGroupBracket,
    SGroupBracketStyle, SGroupCState, SGroupConnection, SGroupData, SGroupDisplay, SdfPropertyList,
    SdfPropertyListTarget, StereoGroup, StereoGroupKind, SubstanceGroup, SubstanceGroupId,
    SubstanceGroupKind,
};

// ──────────────────────────────────────────────
// Format version
// ──────────────────────────────────────────────
const PICKLE_VERSION: u8 = 1;

// ──────────────────────────────────────────────
// Error type
// ──────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PickleError {
    #[error("unexpected end of data while reading pickle")]
    UnexpectedEof,
    #[error("unsupported pickle version: {0}")]
    UnsupportedVersion(u8),
    #[error("data length mismatch: expected {expected}, got {actual}")]
    DataLengthMismatch { expected: usize, actual: usize },
    #[error("invalid enum value: {value} for {type_name}")]
    InvalidEnumValue { value: u8, type_name: &'static str },
    #[error("invalid molecule state after unpickling: {0}")]
    InvalidMolecule(String),
    #[error("too many atoms: {0}")]
    TooManyAtoms(usize),
    #[error("too many bonds: {0}")]
    TooManyBonds(usize),
    #[error("string too long: {0}")]
    StringTooLong(usize),
}

// ──────────────────────────────────────────────
// Binary writer + reader helpers
// ──────────────────────────────────────────────

struct PickleWriter {
    buf: Vec<u8>,
}

impl PickleWriter {
    fn new() -> Self {
        Self { buf: Vec::new() }
    }

    fn into_inner(self) -> Vec<u8> {
        self.buf
    }

    fn write_u8(&mut self, v: u8) {
        self.buf.push(v);
    }

    fn write_u32(&mut self, v: u32) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    fn write_i32(&mut self, v: i32) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    fn write_i8(&mut self, v: i8) {
        self.buf.push(v as u8);
    }

    fn write_f64(&mut self, v: f64) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    fn write_bool(&mut self, v: bool) {
        self.buf.push(if v { 1 } else { 0 });
    }

    fn write_string(&mut self, s: &str) {
        let bytes = s.as_bytes();
        if bytes.len() > u32::MAX as usize {
            // Truncate silently to avoid panic
            let truncated = &bytes[..u32::MAX as usize];
            self.write_u32(truncated.len() as u32);
            self.buf.extend_from_slice(truncated);
        } else {
            self.write_u32(bytes.len() as u32);
            self.buf.extend_from_slice(bytes);
        }
    }

    fn write_option_string(&mut self, s: Option<&str>) {
        match s {
            Some(val) => {
                self.write_bool(true);
                self.write_string(val);
            }
            None => {
                self.write_bool(false);
            }
        }
    }

    fn write_props(&mut self, props: &BTreeMap<String, String>) {
        self.write_u32(props.len() as u32);
        for (key, value) in props {
            self.write_string(key);
            self.write_string(value);
        }
    }
}

struct PickleReader<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> PickleReader<'a> {
    fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0 }
    }

    fn ensure(&self, n: usize) -> Result<(), PickleError> {
        if self.pos + n > self.data.len() {
            Err(PickleError::UnexpectedEof)
        } else {
            Ok(())
        }
    }

    fn read_u8(&mut self) -> Result<u8, PickleError> {
        self.ensure(1)?;
        let v = self.data[self.pos];
        self.pos += 1;
        Ok(v)
    }

    fn read_u32(&mut self) -> Result<u32, PickleError> {
        self.ensure(4)?;
        let bytes: [u8; 4] = self.data[self.pos..self.pos + 4].try_into().unwrap();
        self.pos += 4;
        Ok(u32::from_le_bytes(bytes))
    }

    fn read_i32(&mut self) -> Result<i32, PickleError> {
        self.ensure(4)?;
        let bytes: [u8; 4] = self.data[self.pos..self.pos + 4].try_into().unwrap();
        self.pos += 4;
        Ok(i32::from_le_bytes(bytes))
    }

    fn read_i8(&mut self) -> Result<i8, PickleError> {
        Ok(self.read_u8()? as i8)
    }

    fn read_f64(&mut self) -> Result<f64, PickleError> {
        self.ensure(8)?;
        let bytes: [u8; 8] = self.data[self.pos..self.pos + 8].try_into().unwrap();
        self.pos += 8;
        Ok(f64::from_le_bytes(bytes))
    }

    fn read_bool(&mut self) -> Result<bool, PickleError> {
        Ok(self.read_u8()? != 0)
    }

    fn read_string(&mut self) -> Result<String, PickleError> {
        let len = self.read_u32()? as usize;
        if len > 10_000_000 {
            return Err(PickleError::StringTooLong(len));
        }
        self.ensure(len)?;
        let s = std::str::from_utf8(&self.data[self.pos..self.pos + len])
            .map_err(|_| PickleError::InvalidMolecule("invalid UTF-8 in pickle string".into()))?;
        self.pos += len;
        Ok(s.to_string())
    }

    fn read_option_string(&mut self) -> Result<Option<String>, PickleError> {
        if self.read_bool()? {
            Ok(Some(self.read_string()?))
        } else {
            Ok(None)
        }
    }

    fn read_props(&mut self) -> Result<BTreeMap<String, String>, PickleError> {
        let count = self.read_u32()? as usize;
        if count > 1_000_000 {
            return Err(PickleError::StringTooLong(count));
        }
        let mut props = BTreeMap::new();
        for _ in 0..count {
            let key = self.read_string()?;
            let value = self.read_string()?;
            props.insert(key, value);
        }
        Ok(props)
    }
}

// ──────────────────────────────────────────────
// Enum serialization helpers
// ──────────────────────────────────────────────

fn write_chiral_tag(w: &mut PickleWriter, tag: ChiralTag) {
    let code: u8 = match tag {
        ChiralTag::Unspecified => 0,
        ChiralTag::TetrahedralCw => 1,
        ChiralTag::TetrahedralCcw => 2,
        ChiralTag::Other => 3,
        ChiralTag::Tetrahedral => 4,
        ChiralTag::Allene => 5,
        ChiralTag::SquarePlanar => 6,
        ChiralTag::TrigonalBipyramidal => 7,
        ChiralTag::Octahedral => 8,
    };
    w.write_u8(code);
}

fn read_chiral_tag(r: &mut PickleReader) -> Result<ChiralTag, PickleError> {
    match r.read_u8()? {
        0 => Ok(ChiralTag::Unspecified),
        1 => Ok(ChiralTag::TetrahedralCw),
        2 => Ok(ChiralTag::TetrahedralCcw),
        3 => Ok(ChiralTag::Other),
        4 => Ok(ChiralTag::Tetrahedral),
        5 => Ok(ChiralTag::Allene),
        6 => Ok(ChiralTag::SquarePlanar),
        7 => Ok(ChiralTag::TrigonalBipyramidal),
        8 => Ok(ChiralTag::Octahedral),
        v => Err(PickleError::InvalidEnumValue {
            value: v,
            type_name: "ChiralTag",
        }),
    }
}

fn write_hybridization(w: &mut PickleWriter, h: Hybridization) {
    let code: u8 = match h {
        Hybridization::Unspecified => 0,
        Hybridization::S => 1,
        Hybridization::Sp => 2,
        Hybridization::Sp2 => 3,
        Hybridization::Sp3 => 4,
        Hybridization::Sp2d => 5,
        Hybridization::Sp3d => 6,
        Hybridization::Sp3d2 => 7,
        Hybridization::Other => 8,
    };
    w.write_u8(code);
}

fn read_hybridization(r: &mut PickleReader) -> Result<Hybridization, PickleError> {
    match r.read_u8()? {
        0 => Ok(Hybridization::Unspecified),
        1 => Ok(Hybridization::S),
        2 => Ok(Hybridization::Sp),
        3 => Ok(Hybridization::Sp2),
        4 => Ok(Hybridization::Sp3),
        5 => Ok(Hybridization::Sp2d),
        6 => Ok(Hybridization::Sp3d),
        7 => Ok(Hybridization::Sp3d2),
        8 => Ok(Hybridization::Other),
        v => Err(PickleError::InvalidEnumValue {
            value: v,
            type_name: "Hybridization",
        }),
    }
}

fn write_bond_order(w: &mut PickleWriter, order: BondOrder) {
    let code: u8 = match order {
        BondOrder::Null | BondOrder::Unspecified => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Quintuple => 5,
        BondOrder::Hextuple => 6,
        BondOrder::OneAndHalf => 7,
        BondOrder::TwoAndHalf => 8,
        BondOrder::ThreeAndHalf => 9,
        BondOrder::FourAndHalf => 10,
        BondOrder::FiveAndHalf => 11,
        BondOrder::Aromatic => 12,
        BondOrder::Ionic => 13,
        BondOrder::Dative => 14,
        BondOrder::DativeOne => 15,
        BondOrder::DativeLeft => 16,
        BondOrder::DativeRight => 17,
        BondOrder::Hydrogen => 18,
        BondOrder::ThreeCenter => 19,
        BondOrder::Other => 20,
        BondOrder::Zero => 21,
    };
    w.write_u8(code);
}

fn read_bond_order(r: &mut PickleReader) -> Result<BondOrder, PickleError> {
    match r.read_u8()? {
        0 => Ok(BondOrder::Unspecified),
        1 => Ok(BondOrder::Single),
        2 => Ok(BondOrder::Double),
        3 => Ok(BondOrder::Triple),
        4 => Ok(BondOrder::Quadruple),
        5 => Ok(BondOrder::Quintuple),
        6 => Ok(BondOrder::Hextuple),
        7 => Ok(BondOrder::OneAndHalf),
        8 => Ok(BondOrder::TwoAndHalf),
        9 => Ok(BondOrder::ThreeAndHalf),
        10 => Ok(BondOrder::FourAndHalf),
        11 => Ok(BondOrder::FiveAndHalf),
        12 => Ok(BondOrder::Aromatic),
        13 => Ok(BondOrder::Ionic),
        14 => Ok(BondOrder::Dative),
        15 => Ok(BondOrder::DativeOne),
        16 => Ok(BondOrder::DativeLeft),
        17 => Ok(BondOrder::DativeRight),
        18 => Ok(BondOrder::Hydrogen),
        19 => Ok(BondOrder::ThreeCenter),
        20 => Ok(BondOrder::Other),
        21 => Ok(BondOrder::Zero),
        v => Err(PickleError::InvalidEnumValue {
            value: v,
            type_name: "BondOrder",
        }),
    }
}

fn write_bond_direction(w: &mut PickleWriter, dir: BondDirection) {
    let code: u8 = match dir {
        BondDirection::None => 0,
        BondDirection::BeginWedge => 1,
        BondDirection::BeginDash => 2,
        BondDirection::EndUpRight => 3,
        BondDirection::EndDownRight => 4,
        BondDirection::EitherDouble => 5,
        BondDirection::Unknown => 6,
    };
    w.write_u8(code);
}

fn read_bond_direction(r: &mut PickleReader) -> Result<BondDirection, PickleError> {
    match r.read_u8()? {
        0 => Ok(BondDirection::None),
        1 => Ok(BondDirection::BeginWedge),
        2 => Ok(BondDirection::BeginDash),
        3 => Ok(BondDirection::EndUpRight),
        4 => Ok(BondDirection::EndDownRight),
        5 => Ok(BondDirection::EitherDouble),
        6 => Ok(BondDirection::Unknown),
        v => Err(PickleError::InvalidEnumValue {
            value: v,
            type_name: "BondDirection",
        }),
    }
}

fn write_bond_stereo(w: &mut PickleWriter, stereo: BondStereo) {
    let code: u8 = match stereo {
        BondStereo::None => 0,
        BondStereo::Any => 1,
        BondStereo::Z => 2,
        BondStereo::E => 3,
        BondStereo::Cis => 4,
        BondStereo::Trans => 5,
        BondStereo::AtropCw => 6,
        BondStereo::AtropCcw => 7,
    };
    w.write_u8(code);
}

fn read_bond_stereo(r: &mut PickleReader) -> Result<BondStereo, PickleError> {
    match r.read_u8()? {
        0 => Ok(BondStereo::None),
        1 => Ok(BondStereo::Any),
        2 => Ok(BondStereo::Z),
        3 => Ok(BondStereo::E),
        4 => Ok(BondStereo::Cis),
        5 => Ok(BondStereo::Trans),
        6 => Ok(BondStereo::AtropCw),
        7 => Ok(BondStereo::AtropCcw),
        v => Err(PickleError::InvalidEnumValue {
            value: v,
            type_name: "BondStereo",
        }),
    }
}

// BEGIN RDKIT CPP FUNCTION: MolPickler::pickleAtomProperties (MolPickler.cpp)
// RDKit✔️✔️: void MolPickler::pickleAtomProperties(const ROMol &mol, std::ostream &ss) {
// RDKit✔️✔️:   for (const auto atom : mol.atoms()) {
// RDKit✔️✔️:     pickleAtom(atom, ss);
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION: MolPickler::pickleAtomProperties

fn write_atom(w: &mut PickleWriter, atom: &Atom) {
    // Atomic number
    w.write_u8(atom.atomic_number());

    // Formal charge
    w.write_i8(atom.formal_charge());

    // Isotope
    if let Some(isotope) = atom.isotope() {
        w.write_bool(true);
        w.write_u32(isotope as u32);
    } else {
        w.write_bool(false);
    }

    // Chiral tag
    write_chiral_tag(w, atom.chiral_tag());

    // Chiral permutation
    if let Some(perm) = atom.chiral_permutation() {
        w.write_bool(true);
        w.write_u32(perm);
    } else {
        w.write_bool(false);
    }

    // Unknown stereo
    w.write_bool(atom.unknown_stereo());

    // Mol parity
    if let Some(parity) = atom.mol_parity() {
        w.write_bool(true);
        w.write_i32(parity);
    } else {
        w.write_bool(false);
    }

    // Mol inversion flag
    if let Some(inv) = atom.mol_inversion_flag() {
        w.write_bool(true);
        w.write_i32(inv);
    } else {
        w.write_bool(false);
    }

    // Radical electrons
    w.write_u8(atom.radical_electrons());

    // Is aromatic
    w.write_bool(atom.is_aromatic());

    // Hybridization
    write_hybridization(w, atom.hybridization());

    // Atom map
    if let Some(map) = atom.atom_map() {
        w.write_bool(true);
        w.write_u32(map);
    } else {
        w.write_bool(false);
    }

    // No implicit
    w.write_bool(atom.no_implicit());

    // Implicit hydrogen flag
    w.write_bool(atom.implicit_hydrogen());

    // Explicit hydrogens count
    w.write_u8(atom.explicit_hydrogens());

    // Tracked isotopic hydrogens
    let tracked_isotopes = atom.tracked_isotopic_hydrogens();
    w.write_u32(tracked_isotopes.len() as u32);
    for &iso in tracked_isotopes {
        w.write_u32(iso as u32);
    }

    // Has query flag — we serialize presence but not the full query tree
    w.write_bool(atom.query().is_some());

    // Properties
    w.write_props(atom.props());

    // PDB residue info presence flag
    w.write_bool(atom.pdb_residue_info().is_some());
}

fn read_bond(r: &mut PickleReader) -> Result<Bond, PickleError> {
    let begin_idx = r.read_u32()? as usize;
    let end_idx = r.read_u32()? as usize;
    let order = read_bond_order(r)?;
    let stereo = read_bond_stereo(r)?;
    let direction = read_bond_direction(r)?;
    let is_aromatic = r.read_bool()?;
    let is_conjugated = r.read_bool()?;

    // Stereo atoms
    let stereo_atoms = if r.read_bool()? {
        let sa_begin = AtomId::new(r.read_u32()? as usize);
        let sa_end = AtomId::new(r.read_u32()? as usize);
        Some([sa_begin, sa_end])
    } else {
        None
    };

    let unknown_stereo = r.read_bool()?;

    // Query presence flag (skip full query tree deserialization)
    let _has_query = r.read_bool()?;

    let props = r.read_props()?;

    // Reconstruct via builder-like approach using spec
    // Note: BondSpec is the construction payload; we need id, begin, end
    let spec = crate::BondSpec::new(AtomId::new(begin_idx), AtomId::new(end_idx), order)
        .with_stereo(stereo)
        .with_direction(direction)
        .with_aromatic(is_aromatic)
        .with_conjugated(is_conjugated)
        .with_unknown_stereo(unknown_stereo);

    let spec = if let Some([sa_begin, sa_end]) = stereo_atoms {
        spec.with_stereo_atoms(sa_begin, sa_end)
    } else {
        spec
    };

    // Properties
    let mut spec = spec;
    for (key, value) in &props {
        spec = spec.with_prop(key.clone(), value.clone());
    }

    // Creating from_spec requires BondId
    Ok(Bond::from_spec(BondId::new(begin_idx), spec))
}

// BEGIN RDKIT CPP FUNCTION: MolPickler::pickleBond (MolPickler.cpp)
// RDKit✔️✔️: void MolPickler::pickleBond(const Bond *bond, std::ostream &ss) {
// RDKit✔️✔️:   ss << bond->getBeginAtomIdx();
// RDKit✔️✔️:   ss << bond->getEndAtomIdx();
// RDKit✔️✔️:   ss << bond->getBondType();
// RDKit✔️✔️:   ss << bond->getStereo();
// RDKit✔️✔️:   // ... additional bond properties
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION: MolPickler::pickleBond

fn write_bond(w: &mut PickleWriter, bond: &Bond) {
    w.write_u32(bond.begin().index() as u32);
    w.write_u32(bond.end().index() as u32);
    write_bond_order(w, bond.order());
    write_bond_stereo(w, bond.stereo());
    write_bond_direction(w, bond.direction());
    w.write_bool(bond.is_aromatic());
    w.write_bool(bond.is_conjugated());

    // Stereo atoms
    if let Some([sa_begin, sa_end]) = bond.stereo_atoms() {
        w.write_bool(true);
        w.write_u32(sa_begin.index() as u32);
        w.write_u32(sa_end.index() as u32);
    } else {
        w.write_bool(false);
    }

    w.write_bool(bond.unknown_stereo());

    // Query presence flag
    w.write_bool(bond.query().is_some());

    // Properties
    w.write_props(bond.props());
}

// ──────────────────────────────────────────────
// Substance group serialization
// ──────────────────────────────────────────────

fn write_substance_group_kind(w: &mut PickleWriter, kind: &SubstanceGroupKind) {
    let (code, generic_name): (u8, Option<&str>) = match kind {
        SubstanceGroupKind::Data => (0, None),
        SubstanceGroupKind::Superatom => (1, None),
        SubstanceGroupKind::MultipleGroup => (2, None),
        SubstanceGroupKind::StructuralRepeatUnit => (3, None),
        SubstanceGroupKind::Monomer => (4, None),
        SubstanceGroupKind::Copolymer => (5, None),
        SubstanceGroupKind::Crosslink => (6, None),
        SubstanceGroupKind::Graft => (7, None),
        SubstanceGroupKind::Modification => (8, None),
        SubstanceGroupKind::Mer => (9, None),
        SubstanceGroupKind::AnyPolymer => (10, None),
        SubstanceGroupKind::MixtureComponent => (11, None),
        SubstanceGroupKind::Mixture => (12, None),
        SubstanceGroupKind::Formulation => (13, None),
        SubstanceGroupKind::Generic(name) => (14, Some(name.as_str())),
    };
    w.write_u8(code);
    if let Some(name) = generic_name {
        w.write_string(name);
    }
}

fn read_substance_group_kind(r: &mut PickleReader) -> Result<SubstanceGroupKind, PickleError> {
    match r.read_u8()? {
        0 => Ok(SubstanceGroupKind::Data),
        1 => Ok(SubstanceGroupKind::Superatom),
        2 => Ok(SubstanceGroupKind::MultipleGroup),
        3 => Ok(SubstanceGroupKind::StructuralRepeatUnit),
        4 => Ok(SubstanceGroupKind::Monomer),
        5 => Ok(SubstanceGroupKind::Copolymer),
        6 => Ok(SubstanceGroupKind::Crosslink),
        7 => Ok(SubstanceGroupKind::Graft),
        8 => Ok(SubstanceGroupKind::Modification),
        9 => Ok(SubstanceGroupKind::Mer),
        10 => Ok(SubstanceGroupKind::AnyPolymer),
        11 => Ok(SubstanceGroupKind::MixtureComponent),
        12 => Ok(SubstanceGroupKind::Mixture),
        13 => Ok(SubstanceGroupKind::Formulation),
        14 => {
            let name = r.read_string()?;
            Ok(SubstanceGroupKind::Generic(name))
        }
        v => Err(PickleError::InvalidEnumValue {
            value: v,
            type_name: "SubstanceGroupKind",
        }),
    }
}

fn write_sgroup_connection(w: &mut PickleWriter, conn: Option<&SGroupConnection>) {
    match conn {
        None => w.write_u8(0),
        Some(SGroupConnection::HeadToHead) => w.write_u8(1),
        Some(SGroupConnection::HeadToTail) => w.write_u8(2),
        Some(SGroupConnection::Either) => w.write_u8(3),
        Some(SGroupConnection::Unknown(s)) => {
            w.write_u8(4);
            w.write_string(s);
        }
    }
}

fn read_sgroup_connection(r: &mut PickleReader) -> Result<Option<SGroupConnection>, PickleError> {
    match r.read_u8()? {
        0 => Ok(None),
        1 => Ok(Some(SGroupConnection::HeadToHead)),
        2 => Ok(Some(SGroupConnection::HeadToTail)),
        3 => Ok(Some(SGroupConnection::Either)),
        4 => Ok(Some(SGroupConnection::Unknown(r.read_string()?))),
        v => Err(PickleError::InvalidEnumValue {
            value: v,
            type_name: "SGroupConnection",
        }),
    }
}

fn write_sgroup_bracket_style(w: &mut PickleWriter, style: Option<&SGroupBracketStyle>) {
    match style {
        None => w.write_u8(0),
        Some(SGroupBracketStyle::Bracket) => w.write_u8(1),
        Some(SGroupBracketStyle::Parenthesis) => w.write_u8(2),
        Some(SGroupBracketStyle::None) => w.write_u8(3),
        Some(SGroupBracketStyle::Unknown(s)) => {
            w.write_u8(4);
            w.write_string(s);
        }
    }
}

fn read_sgroup_bracket_style(
    r: &mut PickleReader,
) -> Result<Option<SGroupBracketStyle>, PickleError> {
    match r.read_u8()? {
        0 => Ok(None),
        1 => Ok(Some(SGroupBracketStyle::Bracket)),
        2 => Ok(Some(SGroupBracketStyle::Parenthesis)),
        3 => Ok(Some(SGroupBracketStyle::None)),
        4 => Ok(Some(SGroupBracketStyle::Unknown(r.read_string()?))),
        v => Err(PickleError::InvalidEnumValue {
            value: v,
            type_name: "SGroupBracketStyle",
        }),
    }
}

fn write_sgroup_display(w: &mut PickleWriter, display: Option<&SGroupDisplay>) {
    match display {
        None => w.write_bool(false),
        Some(d) => {
            w.write_bool(true);
            // Brackets
            w.write_u32(d.brackets.len() as u32);
            for bracket in &d.brackets {
                w.write_f64(bracket.p1[0]);
                w.write_f64(bracket.p1[1]);
                w.write_f64(bracket.p2[0]);
                w.write_f64(bracket.p2[1]);
            }
            // Field position
            match d.field_position {
                Some(pos) => {
                    w.write_bool(true);
                    w.write_f64(pos[0]);
                    w.write_f64(pos[1]);
                }
                None => w.write_bool(false),
            }
            // Display tag
            w.write_option_string(d.display_tag.as_deref());
        }
    }
}

fn read_sgroup_display(r: &mut PickleReader) -> Result<Option<SGroupDisplay>, PickleError> {
    if !r.read_bool()? {
        return Ok(None);
    }
    let mut display = SGroupDisplay::default();
    let bracket_count = r.read_u32()? as usize;
    for _ in 0..bracket_count {
        let p1x = r.read_f64()?;
        let p1y = r.read_f64()?;
        let p2x = r.read_f64()?;
        let p2y = r.read_f64()?;
        display.brackets.push(SGroupBracket {
            p1: [p1x, p1y],
            p2: [p2x, p2y],
        });
    }
    if r.read_bool()? {
        let fx = r.read_f64()?;
        let fy = r.read_f64()?;
        display.field_position = Some([fx, fy]);
    }
    display.display_tag = r.read_option_string()?;
    Ok(Some(display))
}

fn write_sgroup_data(w: &mut PickleWriter, data: Option<&SGroupData>) {
    match data {
        None => w.write_bool(false),
        Some(d) => {
            w.write_bool(true);
            w.write_option_string(d.field_name.as_deref());
            w.write_option_string(d.field_type.as_deref());
            w.write_option_string(d.field_info.as_deref());
            w.write_option_string(d.field_display.as_deref());
            w.write_option_string(d.units.as_deref());
            w.write_option_string(d.query_type.as_deref());
            w.write_option_string(d.query_op.as_deref());
            w.write_u32(d.values.len() as u32);
            for v in &d.values {
                w.write_string(v);
            }
        }
    }
}

fn read_sgroup_data(r: &mut PickleReader) -> Result<Option<SGroupData>, PickleError> {
    if !r.read_bool()? {
        return Ok(None);
    }
    let mut data = SGroupData::default();
    data.field_name = r.read_option_string()?;
    data.field_type = r.read_option_string()?;
    data.field_info = r.read_option_string()?;
    data.field_display = r.read_option_string()?;
    data.units = r.read_option_string()?;
    data.query_type = r.read_option_string()?;
    data.query_op = r.read_option_string()?;
    let val_count = r.read_u32()? as usize;
    for _ in 0..val_count {
        data.values.push(r.read_string()?);
    }
    Ok(Some(data))
}

fn write_substance_group(w: &mut PickleWriter, sg: &SubstanceGroup) {
    // id index
    w.write_u32(sg.id().index() as u32);

    // rdkit sequence id
    if let Some(seq_id) = sg.rdkit_sequence_id() {
        w.write_bool(true);
        w.write_u32(seq_id);
    } else {
        w.write_bool(false);
    }

    // external id
    if let Some(ext_id) = sg.external_id() {
        w.write_bool(true);
        w.write_u32(ext_id);
    } else {
        w.write_bool(false);
    }

    // kind
    write_substance_group_kind(w, sg.kind());

    // atoms
    let atoms = sg.atoms();
    w.write_u32(atoms.len() as u32);
    for &a in atoms {
        w.write_u32(a.index() as u32);
    }

    // bonds
    let bonds = sg.bonds();
    w.write_u32(bonds.len() as u32);
    for &b in bonds {
        w.write_u32(b.index() as u32);
    }

    // bond roles — only write non-default roles
    let mut roles_written = 0u32;
    for &b in bonds {
        if sg.bond_role(b) == SGroupBondRole::Contained {
            roles_written += 1;
        }
    }
    w.write_u32(roles_written);
    for &b in bonds {
        if sg.bond_role(b) == SGroupBondRole::Contained {
            w.write_u32(b.index() as u32);
        }
    }

    // parent atoms
    let parent_atoms = sg.parent_atoms();
    w.write_u32(parent_atoms.len() as u32);
    for &a in parent_atoms {
        w.write_u32(a.index() as u32);
    }

    // parent
    if let Some(parent) = sg.parent() {
        w.write_bool(true);
        w.write_u32(parent.index() as u32);
    } else {
        w.write_bool(false);
    }

    // label
    w.write_option_string(sg.label());

    // connection
    write_sgroup_connection(w, sg.connection());

    // subtype
    w.write_option_string(sg.subtype());

    // bracket style
    write_sgroup_bracket_style(w, sg.bracket_style());

    // expansion state
    w.write_option_string(sg.expansion_state());

    // class
    w.write_option_string(sg.class());

    // component number
    if let Some(cn) = sg.component_number() {
        w.write_bool(true);
        w.write_u32(cn);
    } else {
        w.write_bool(false);
    }

    // display
    write_sgroup_display(w, sg.display());

    // data
    write_sgroup_data(w, sg.data());

    // attach points
    let attach_pts = sg.attach_points();
    w.write_u32(attach_pts.len() as u32);
    for ap in attach_pts {
        w.write_u32(ap.atom.index() as u32);
        if let Some(la) = ap.leaving_atom {
            w.write_bool(true);
            w.write_u32(la.index() as u32);
        } else {
            w.write_bool(false);
        }
        w.write_option_string(ap.label.as_deref());
        if let Some(order) = ap.order {
            w.write_bool(true);
            w.write_u32(order);
        } else {
            w.write_bool(false);
        }
    }

    // cstates
    let cstates = sg.cstates();
    w.write_u32(cstates.len() as u32);
    for cs in cstates {
        w.write_u32(cs.bond.index() as u32);
        w.write_f64(cs.vector[0]);
        w.write_f64(cs.vector[1]);
    }

    // props
    w.write_props(sg.props());

    // data fields
    let data_fields = sg.data_fields();
    w.write_u32(data_fields.len() as u32);
    for df in data_fields {
        w.write_string(df);
    }
}

// ──────────────────────────────────────────────
// Stereo group serialization
// ──────────────────────────────────────────────

fn write_stereo_group_kind(w: &mut PickleWriter, kind: StereoGroupKind) {
    let code: u8 = match kind {
        StereoGroupKind::Absolute => 0,
        StereoGroupKind::Or => 1,
        StereoGroupKind::And => 2,
    };
    w.write_u8(code);
}

fn read_stereo_group_kind(r: &mut PickleReader) -> Result<StereoGroupKind, PickleError> {
    match r.read_u8()? {
        0 => Ok(StereoGroupKind::Absolute),
        1 => Ok(StereoGroupKind::Or),
        2 => Ok(StereoGroupKind::And),
        v => Err(PickleError::InvalidEnumValue {
            value: v,
            type_name: "StereoGroupKind",
        }),
    }
}

fn write_stereo_group(w: &mut PickleWriter, sg: &StereoGroup) {
    // id
    if let Some(id) = sg.id() {
        w.write_bool(true);
        w.write_u32(id);
    } else {
        w.write_bool(false);
    }

    // kind
    write_stereo_group_kind(w, sg.kind());

    // atoms
    let atoms = sg.atoms();
    w.write_u32(atoms.len() as u32);
    for &a in atoms {
        w.write_u32(a.index() as u32);
    }

    // bonds
    let bonds = sg.bonds();
    w.write_u32(bonds.len() as u32);
    for &b in bonds {
        w.write_u32(b.index() as u32);
    }
}

// ──────────────────────────────────────────────
// Main serialization / deserialization
// ──────────────────────────────────────────────

/// Serialize a `Molecule` to a compact binary format.
///
/// Format structure:
/// - Version byte
/// - Atom count + atom data
/// - Bond count + bond data
/// - 2D coordinates (optional)
/// - 3D conformers
/// - Substance groups
/// - Stereo groups
/// - Molecule properties
///
/// # Errors
///
/// Returns `PickleError` if serialization encounters an internal issue.
pub fn mol_to_binary(mol: &Molecule) -> Result<Vec<u8>, PickleError> {
    let mut w = PickleWriter::new();

    // Version
    w.write_u8(PICKLE_VERSION);

    // ── Atoms ──
    let atoms = mol.atoms();
    if atoms.len() > u32::MAX as usize {
        return Err(PickleError::TooManyAtoms(atoms.len()));
    }
    w.write_u32(atoms.len() as u32);
    for atom in atoms {
        write_atom(&mut w, atom);
    }

    // ── Bonds ──
    let bonds = mol.bonds();
    if bonds.len() > u32::MAX as usize {
        return Err(PickleError::TooManyBonds(bonds.len()));
    }
    w.write_u32(bonds.len() as u32);
    for bond in bonds {
        write_bond(&mut w, bond);
    }

    // ── 2D Coordinates ──
    if let Some(coords_2d) = mol.coords_2d() {
        w.write_bool(true);
        if coords_2d.len() > u32::MAX as usize {
            return Err(PickleError::TooManyAtoms(coords_2d.len()));
        }
        w.write_u32(coords_2d.len() as u32);
        for &[x, y] in coords_2d {
            w.write_f64(x);
            w.write_f64(y);
        }
    } else {
        w.write_bool(false);
    }

    // ── 3D Conformers ──
    let conformers = mol.conformers_3d();
    w.write_u32(conformers.len() as u32);
    for conf in conformers {
        w.write_u32(conf.id() as u32);
        let coords = conf.coords();
        w.write_u32(coords.len() as u32);
        for &[x, y, z] in coords {
            w.write_f64(x);
            w.write_f64(y);
            w.write_f64(z);
        }
        w.write_bool(conf.is_3d());
        w.write_props(conf.props());
    }

    // ── Source coordinate dimension ──
    match mol.source_coordinate_dim() {
        None => w.write_u8(0),
        Some(CoordinateDimension::TwoD) => w.write_u8(1),
        Some(CoordinateDimension::ThreeD) => w.write_u8(2),
    }

    // ── Substance Groups ──
    let sgroups = mol.substance_groups();
    w.write_u32(sgroups.len() as u32);
    for sg in sgroups {
        write_substance_group(&mut w, sg);
    }

    // ── Stereo Groups ──
    let stereo_groups = mol.stereo_groups();
    w.write_u32(stereo_groups.len() as u32);
    for sg in stereo_groups {
        write_stereo_group(&mut w, sg);
    }

    // ── Molecule Properties ──
    let props = mol.properties();
    w.write_option_string(props.name());
    w.write_props(props.props());

    // SDF data fields
    let sdf_fields = props.sdf_data_fields();
    w.write_u32(sdf_fields.len() as u32);
    for (key, value) in sdf_fields {
        w.write_string(key);
        w.write_string(value);
    }

    // SDF property lists
    let sdf_prop_lists = props.sdf_property_lists();
    w.write_u32(sdf_prop_lists.len() as u32);
    for plist in sdf_prop_lists {
        match plist.target() {
            SdfPropertyListTarget::Atom => w.write_u8(0),
            SdfPropertyListTarget::Bond => w.write_u8(1),
        }
        w.write_string(plist.name());
        let values = plist.values();
        w.write_u32(values.len() as u32);
        for v in values {
            w.write_option_string(v.as_deref());
        }
    }

    // ── Finish ──
    Ok(w.into_inner())
}

// BEGIN RDKIT CPP FUNCTION: MolPickler::unpickleMol (MolPickler.cpp)
// RDKit✔️✔️: void MolPickler::unpickleMol(std::istream &ss, RWMol &mol) {
// RDKit✔️✔️:   // Deserializes molecule state from binary stream.
// RDKit✔️✔️:   // First byte is the pickle version.
// RDKit✔️✔️:   int version;
// RDKit✔️✔️:   ss >> version;
// RDKit✔️✔️:   // Version-dependent dispatch...
// RDKit✔️✔️:   MolPickler::unpickleMol(ss, mol, version);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION: MolPickler::unpickleMol

/// Deserialize a `Molecule` from binary data produced by `mol_to_binary`.
///
/// # Errors
///
/// Returns `PickleError` if the data is corrupt, has an unsupported version,
/// or produces an invalid molecule state.
pub fn mol_from_binary(data: &[u8]) -> Result<Molecule, PickleError> {
    let mut r = PickleReader::new(data);

    // Version check
    let version = r.read_u8()?;
    if version > PICKLE_VERSION {
        return Err(PickleError::UnsupportedVersion(version));
    }

    // ── Atoms ──
    let atom_count = r.read_u32()? as usize;
    if atom_count > 1_000_000 {
        return Err(PickleError::TooManyAtoms(atom_count));
    }

    let mut atom_specs = Vec::with_capacity(atom_count);
    for i in 0..atom_count {
        let atomic_number = r.read_u8()?;
        let formal_charge = r.read_i8()?;
        let isotope = if r.read_bool()? {
            Some(r.read_u32()? as u16)
        } else {
            None
        };
        let chiral_tag = read_chiral_tag(&mut r)?;
        let chiral_perm = if r.read_bool()? {
            Some(r.read_u32()?)
        } else {
            None
        };
        let unknown_stereo = r.read_bool()?;
        let mol_parity = if r.read_bool()? {
            Some(r.read_i32()?)
        } else {
            None
        };
        let mol_inv_flag = if r.read_bool()? {
            Some(r.read_i32()?)
        } else {
            None
        };
        let radical_electrons = r.read_u8()?;
        let is_aromatic = r.read_bool()?;
        let hybridization = read_hybridization(&mut r)?;
        let atom_map = if r.read_bool()? {
            Some(r.read_u32()?)
        } else {
            None
        };
        let no_implicit = r.read_bool()?;
        let implicit_hydrogen = r.read_bool()?;
        let explicit_hydrogens = r.read_u8()?;
        let tracked_isotope_count = r.read_u32()? as usize;
        let mut tracked_isotopic_hydrogens = Vec::with_capacity(tracked_isotope_count);
        for _ in 0..tracked_isotope_count {
            tracked_isotopic_hydrogens.push(r.read_u32()? as u16);
        }
        let _has_query = r.read_bool()?;
        let props = r.read_props()?;
        let _has_pdb_info = r.read_bool()?;

        let element =
            crate::Element::from_atomic_number(atomic_number).unwrap_or(crate::Element::DUMMY);

        let mut spec = crate::AtomSpec::new(element)
            .with_formal_charge(formal_charge)
            .with_chiral_tag(chiral_tag)
            .with_radical_electrons(radical_electrons)
            .with_aromatic(is_aromatic)
            .with_hybridization(hybridization)
            .with_no_implicit(no_implicit)
            .with_implicit_hydrogen(implicit_hydrogen)
            .with_explicit_hydrogens(explicit_hydrogens)
            .with_unknown_stereo(unknown_stereo);

        if let Some(iso) = isotope {
            spec = spec.with_isotope(iso);
        }
        if let Some(perm) = chiral_perm {
            spec = spec.with_chiral_permutation(perm);
        }
        if let Some(map) = atom_map {
            spec = spec.with_atom_map(map);
        }
        if let Some(parity) = mol_parity {
            spec = spec.with_mol_parity(parity);
        }
        if let Some(inv) = mol_inv_flag {
            spec = spec.with_mol_inversion_flag(inv);
        }
        if !tracked_isotopic_hydrogens.is_empty() {
            spec = spec.with_tracked_isotopic_hydrogens(tracked_isotopic_hydrogens);
        }

        // Properties
        for (key, value) in &props {
            spec = spec.with_prop(key.clone(), value.clone());
        }

        atom_specs.push((i, spec));
    }

    // ── Bonds ──
    let bond_count = r.read_u32()? as usize;
    if bond_count > 1_000_000 {
        return Err(PickleError::TooManyBonds(bond_count));
    }

    let mut bond_specs = Vec::with_capacity(bond_count);
    for _ in 0..bond_count {
        let bond = read_bond(&mut r)?;
        // Store just the spec components we need
        bond_specs.push((
            bond.begin(),
            bond.end(),
            bond.order(),
            bond.stereo(),
            bond.direction(),
            bond.is_aromatic(),
            bond.is_conjugated(),
            bond.stereo_atoms(),
            bond.unknown_stereo(),
            bond.props().clone(),
        ));
    }

    // ── 2D Coordinates ──
    let mut coords_2d: Option<Vec<[f64; 2]>> = None;
    if r.read_bool()? {
        let coord_count = r.read_u32()? as usize;
        let mut coords = Vec::with_capacity(coord_count);
        for _ in 0..coord_count {
            let x = r.read_f64()?;
            let y = r.read_f64()?;
            coords.push([x, y]);
        }
        coords_2d = Some(coords);
    }

    // ── 3D Conformers ──
    let conformer_count = r.read_u32()? as usize;
    let mut conformers_3d = Vec::with_capacity(conformer_count);
    for _ in 0..conformer_count {
        let conf_id = r.read_u32()? as usize;
        let coord_count = r.read_u32()? as usize;
        let mut coords = Vec::with_capacity(coord_count);
        for _ in 0..coord_count {
            let x = r.read_f64()?;
            let y = r.read_f64()?;
            let z = r.read_f64()?;
            coords.push([x, y, z]);
        }
        let is_3d = r.read_bool()?;
        let props = r.read_props()?;

        let mut conformer = Conformer3D::new(conf_id, coords, is_3d);
        // Replay props using the builder pattern — Conformer3D doesn't expose set_prop directly
        // but has with_prop (consumes self). We must reconstruct.
        // Actually Conformer3D::new creates with empty props, so we use the with_prop pattern
        // but since with_prop takes self, we need to re-architect for the loop.
        // Use crate::Conformer3D constructor then... actually let's just make a new one each time.
        for (k, v) in &props {
            conformer = conformer.with_prop(k.clone(), v.clone());
        }
        conformers_3d.push(conformer);
    }

    // ── Source coordinate dimension ──
    let _source_coordinate_dim = match r.read_u8()? {
        0 => None,
        1 => Some(CoordinateDimension::TwoD),
        2 => Some(CoordinateDimension::ThreeD),
        _ => None,
    };

    // ── Substance Groups ──
    let sgroup_count = r.read_u32()? as usize;
    let mut sgroups = Vec::with_capacity(sgroup_count);
    for _ in 0..sgroup_count {
        let id = SubstanceGroupId::new(r.read_u32()? as usize);
        let has_rdkit_seq = r.read_bool()?;
        let rdkit_seq = if has_rdkit_seq {
            Some(r.read_u32()?)
        } else {
            None
        };
        let has_ext_id = r.read_bool()?;
        let ext_id = if has_ext_id {
            Some(r.read_u32()?)
        } else {
            None
        };
        let kind = read_substance_group_kind(&mut r)?;
        let atom_count_sg = r.read_u32()? as usize;
        let mut atoms_sg = Vec::with_capacity(atom_count_sg);
        for _ in 0..atom_count_sg {
            atoms_sg.push(AtomId::new(r.read_u32()? as usize));
        }
        let bond_count_sg = r.read_u32()? as usize;
        let mut bonds_sg = Vec::with_capacity(bond_count_sg);
        for _ in 0..bond_count_sg {
            bonds_sg.push(BondId::new(r.read_u32()? as usize));
        }
        let role_count = r.read_u32()? as usize;
        let mut bond_roles = BTreeMap::new();
        for _ in 0..role_count {
            let b = BondId::new(r.read_u32()? as usize);
            bond_roles.insert(b, SGroupBondRole::Contained);
        }
        let parent_atom_count = r.read_u32()? as usize;
        let mut parent_atoms = Vec::with_capacity(parent_atom_count);
        for _ in 0..parent_atom_count {
            parent_atoms.push(AtomId::new(r.read_u32()? as usize));
        }
        let has_parent = r.read_bool()?;
        let parent = if has_parent {
            Some(SubstanceGroupId::new(r.read_u32()? as usize))
        } else {
            None
        };
        let label = r.read_option_string()?;
        let connection = read_sgroup_connection(&mut r)?;
        let subtype = r.read_option_string()?;
        let bracket_style = read_sgroup_bracket_style(&mut r)?;
        let expansion_state = r.read_option_string()?;
        let class = r.read_option_string()?;
        let has_component_number = r.read_bool()?;
        let component_number = if has_component_number {
            Some(r.read_u32()?)
        } else {
            None
        };
        let display = read_sgroup_display(&mut r)?;
        let data = read_sgroup_data(&mut r)?;

        // Attach points
        let ap_count = r.read_u32()? as usize;
        let mut attach_points = Vec::with_capacity(ap_count);
        for _ in 0..ap_count {
            let ap_atom = AtomId::new(r.read_u32()? as usize);
            let leaving = if r.read_bool()? {
                Some(AtomId::new(r.read_u32()? as usize))
            } else {
                None
            };
            let ap_label = r.read_option_string()?;
            let has_ap_order = r.read_bool()?;
            let ap_order = if has_ap_order {
                Some(r.read_u32()?)
            } else {
                None
            };
            attach_points.push(SGroupAttachPoint {
                atom: ap_atom,
                leaving_atom: leaving,
                label: ap_label,
                order: ap_order,
            });
        }

        // CStates
        let cs_count = r.read_u32()? as usize;
        let mut cstates = Vec::with_capacity(cs_count);
        for _ in 0..cs_count {
            let cs_bond = BondId::new(r.read_u32()? as usize);
            let cs_x = r.read_f64()?;
            let cs_y = r.read_f64()?;
            cstates.push(SGroupCState {
                bond: cs_bond,
                vector: [cs_x, cs_y],
            });
        }

        let props = r.read_props()?;
        let data_field_count = r.read_u32()? as usize;
        let mut data_fields = Vec::with_capacity(data_field_count);
        for _ in 0..data_field_count {
            data_fields.push(r.read_string()?);
        }

        let mut sg = SubstanceGroup::new(id, kind)
            .with_atoms(atoms_sg)
            .with_bonds(bonds_sg)
            .with_parent_atoms(parent_atoms)
            .with_attach_points(attach_points)
            .with_cstates(cstates);

        if let Some(seq) = rdkit_seq {
            sg = sg.with_rdkit_sequence_id(seq);
        }
        if let Some(eid) = ext_id {
            sg = sg.with_external_id(eid);
        }
        if let Some(p) = parent {
            sg = sg.with_parent(p);
        }
        if let Some(l) = label {
            sg = sg.with_label(l);
        }
        if let Some(conn) = connection {
            sg = sg.with_connection(conn);
        }
        if let Some(st) = subtype {
            sg = sg.with_subtype(st);
        }
        if let Some(bs) = bracket_style {
            sg = sg.with_bracket_style(bs);
        }
        if let Some(disp) = display {
            sg = sg.with_display(disp);
        }
        if let Some(es) = expansion_state {
            sg = sg.with_expansion_state(es);
        }
        if let Some(c) = class {
            sg = sg.with_class(c);
        }
        if let Some(cn) = component_number {
            sg = sg.with_component_number(cn);
        }
        if let Some(d) = data {
            sg = sg.with_data(d);
        }
        for (key, value) in &props {
            sg = sg.with_prop(key.clone(), value.clone());
        }
        for df in &data_fields {
            sg = sg.with_data_field(df.clone());
        }
        // Apply bond roles (need to call after bonds are set)
        for (bond, role) in &bond_roles {
            if *role == SGroupBondRole::Contained {
                sg = sg.with_bond_role(*bond, SGroupBondRole::Contained);
            }
        }

        sgroups.push(sg);
    }

    // ── Stereo Groups ──
    let stereo_group_count = r.read_u32()? as usize;
    let mut stereo_groups = Vec::with_capacity(stereo_group_count);
    for _ in 0..stereo_group_count {
        let has_id = r.read_bool()?;
        let sg_id = if has_id { Some(r.read_u32()?) } else { None };
        let kind = read_stereo_group_kind(&mut r)?;
        let atom_count_sg = r.read_u32()? as usize;
        let mut atoms_sg = Vec::with_capacity(atom_count_sg);
        for _ in 0..atom_count_sg {
            atoms_sg.push(AtomId::new(r.read_u32()? as usize));
        }
        let bond_count_sg = r.read_u32()? as usize;
        let mut bonds_sg = Vec::with_capacity(bond_count_sg);
        for _ in 0..bond_count_sg {
            bonds_sg.push(BondId::new(r.read_u32()? as usize));
        }
        let mut sg = StereoGroup::new(kind, atoms_sg, bonds_sg);
        if let Some(id) = sg_id {
            sg = sg.with_id(id);
        }
        stereo_groups.push(sg);
    }

    // ── Molecule Properties ──
    let prop_name = r.read_option_string()?;
    let props = r.read_props()?;
    let sdf_field_count = r.read_u32()? as usize;
    let mut sdf_data_fields = Vec::with_capacity(sdf_field_count);
    for _ in 0..sdf_field_count {
        let key = r.read_string()?;
        let value = r.read_string()?;
        sdf_data_fields.push((key, value));
    }
    let sdf_plist_count = r.read_u32()? as usize;
    let mut sdf_property_lists = Vec::with_capacity(sdf_plist_count);
    for _ in 0..sdf_plist_count {
        let target = match r.read_u8()? {
            0 => SdfPropertyListTarget::Atom,
            1 => SdfPropertyListTarget::Bond,
            _ => {
                return Err(PickleError::InvalidEnumValue {
                    value: r.read_u8()?,
                    type_name: "SdfPropertyListTarget",
                });
            }
        };
        let name = r.read_string()?;
        let val_count = r.read_u32()? as usize;
        let mut values = Vec::with_capacity(val_count);
        for _ in 0..val_count {
            values.push(r.read_option_string()?);
        }
        sdf_property_lists.push(SdfPropertyList::new(target, name, values));
    }

    // ── Build Molecule ──
    // Use MoleculeBuilder to construct the molecule
    let mut builder = crate::MoleculeBuilder::new();

    for (_, spec) in &atom_specs {
        let _ = builder.add_atom(spec.clone());
    }

    // Add bonds
    for (
        begin,
        end,
        order,
        stereo,
        direction,
        is_aromatic,
        is_conjugated,
        stereo_atoms,
        unknown_stereo,
        bond_props,
    ) in &bond_specs
    {
        let mut bspec = crate::BondSpec::new(*begin, *end, *order)
            .with_stereo(*stereo)
            .with_direction(*direction)
            .with_aromatic(*is_aromatic)
            .with_conjugated(*is_conjugated)
            .with_unknown_stereo(*unknown_stereo);
        if let Some([sa_begin, sa_end]) = stereo_atoms {
            bspec = bspec.with_stereo_atoms(*sa_begin, *sa_end);
        }
        for (key, value) in bond_props {
            bspec = bspec.with_prop(key.clone(), value.clone());
        }
        builder
            .add_bond(bspec)
            .map_err(|e| PickleError::InvalidMolecule(e.to_string()))?;
    }

    // 2D coordinates
    if let Some(coords) = &coords_2d {
        builder
            .set_2d_coordinates(coords.clone())
            .map_err(|e| PickleError::InvalidMolecule(e.to_string()))?;
    }

    // 3D conformers
    for conformer in &conformers_3d {
        if conformer.coords().len() != atom_count {
            return Err(PickleError::InvalidMolecule(format!(
                "3D conformer row count mismatch: rows={}, atom_count={}",
                conformer.coords().len(),
                atom_count
            )));
        }
        builder
            .add_conformer(conformer.clone())
            .map_err(|e| PickleError::InvalidMolecule(e.to_string()))?;
    }

    // Substance groups
    for sg in &sgroups {
        builder
            .add_substance_group(sg.clone())
            .map_err(|e| PickleError::InvalidMolecule(e.to_string()))?;
    }

    // Stereo groups
    for sg in &stereo_groups {
        builder
            .add_stereo_group(sg.clone())
            .map_err(|e| PickleError::InvalidMolecule(e.to_string()))?;
    }

    // Molecule properties — construct complete MoleculeProperties and set it
    let mut mol_props = crate::MoleculeProperties::default();
    if let Some(name) = &prop_name {
        mol_props = mol_props.with_name(name.clone());
    }
    for (key, value) in &props {
        mol_props = mol_props.with_prop(key.clone(), value.clone());
    }
    for (key, value) in &sdf_data_fields {
        mol_props = mol_props.with_sdf_data_field(key.clone(), value.clone());
    }
    for plist in &sdf_property_lists {
        mol_props = mol_props.with_sdf_property_list(plist.clone());
    }
    builder = builder.with_properties(mol_props);

    // Build
    builder
        .build()
        .map_err(|e| PickleError::InvalidMolecule(e.to_string()))
}

// ──────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        AtomSpec, BondOrder, BondSpec, BondStereo, ChiralTag, Element, Hybridization,
        MoleculeBuilder, SdfPropertyList, SdfPropertyListTarget, StereoGroup, StereoGroupKind,
    };

    fn build_simple_methane() -> Molecule {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(c));
        for _ in 0..4 {
            builder.add_atom(AtomSpec::new(h));
        }
        for i in 0..4 {
            builder
                .add_bond(BondSpec::new(
                    AtomId::new(0),
                    AtomId::new(i + 1),
                    BondOrder::Single,
                ))
                .unwrap();
        }
        builder.build().expect("build methane")
    }

    fn build_simple_ethanol() -> Molecule {
        let c = Element::from_atomic_number(6).unwrap();
        let o = Element::from_atomic_number(8).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = MoleculeBuilder::new();
        // C-C-O-H backbone: indices 0,1,2,3
        builder.add_atom(AtomSpec::new(c));
        builder.add_atom(
            AtomSpec::new(c)
                .with_chiral_tag(ChiralTag::TetrahedralCw)
                .with_chiral_permutation(0),
        );
        builder.add_atom(AtomSpec::new(o));
        builder.add_atom(AtomSpec::new(h));
        // 3 H on C0
        for _ in 0..3 {
            builder.add_atom(AtomSpec::new(h));
        }
        // 2 H on C1
        for _ in 0..2 {
            builder.add_atom(AtomSpec::new(h));
        }
        // 1 H on O
        builder.add_atom(AtomSpec::new(h));

        // C-C single
        builder
            .add_bond(BondSpec::new(
                AtomId::new(0),
                AtomId::new(1),
                BondOrder::Single,
            ))
            .unwrap();
        // C-O single
        builder
            .add_bond(BondSpec::new(
                AtomId::new(1),
                AtomId::new(2),
                BondOrder::Single,
            ))
            .unwrap();
        // O-H single
        builder
            .add_bond(BondSpec::new(
                AtomId::new(2),
                AtomId::new(3),
                BondOrder::Single,
            ))
            .unwrap();
        // 3 C-H on C0
        for i in 0..3 {
            builder
                .add_bond(BondSpec::new(
                    AtomId::new(0),
                    AtomId::new(4 + i),
                    BondOrder::Single,
                ))
                .unwrap();
        }
        // 2 C-H on C1
        for i in 0..2 {
            builder
                .add_bond(BondSpec::new(
                    AtomId::new(1),
                    AtomId::new(7 + i),
                    BondOrder::Single,
                ))
                .unwrap();
        }
        // O-H already on index 3
        builder
            .add_bond(BondSpec::new(
                AtomId::new(2),
                AtomId::new(9),
                BondOrder::Single,
            ))
            .unwrap();

        builder.build().expect("build ethanol")
    }

    #[test]
    fn test_empty_molecule_roundtrip() {
        let mol = Molecule::new();
        let data = mol_to_binary(&mol).unwrap();
        let mol2 = mol_from_binary(&data).unwrap();
        assert_eq!(mol, mol2, "empty molecule roundtrip failed");
    }

    #[test]
    fn test_methane_roundtrip() {
        let mol = build_simple_methane();
        assert_eq!(mol.num_atoms(), 5);
        assert_eq!(mol.num_bonds(), 4);

        let data = mol_to_binary(&mol).unwrap();
        let mol2 = mol_from_binary(&data).unwrap();
        assert_eq!(mol, mol2, "methane roundtrip failed");
    }

    #[test]
    fn test_methane_with_props_roundtrip() {
        let mut mol = build_simple_methane();
        mol = mol.with_name("methane_test");
        mol = mol.with_prop("key1", "value1");
        mol = mol.with_prop("key2", "value2");

        let data = mol_to_binary(&mol).unwrap();
        let mol2 = mol_from_binary(&data).unwrap();

        assert_eq!(mol.properties().name(), Some("methane_test"));
        assert_eq!(mol2.properties().name(), Some("methane_test"));
        assert_eq!(mol2.prop("key1"), Some("value1"));
        assert_eq!(mol2.prop("key2"), Some("value2"));
        assert_eq!(mol, mol2, "methane with properties roundtrip failed");
    }

    #[test]
    fn test_methane_with_2d_coords() {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(c));
        for _ in 0..4 {
            builder.add_atom(AtomSpec::new(h));
        }
        for i in 0..4 {
            builder
                .add_bond(BondSpec::new(
                    AtomId::new(0),
                    AtomId::new(i + 1),
                    BondOrder::Single,
                ))
                .unwrap();
        }
        builder
            .set_2d_coordinates(vec![
                [0.0, 0.0],
                [1.0, 0.0],
                [-0.5, 0.866],
                [-0.5, -0.866],
                [0.0, 1.0],
            ])
            .unwrap();
        let mol = builder.build().expect("build methane with coords");
        assert!(mol.coords_2d().is_some());
        assert_eq!(mol.coords_2d().unwrap().len(), 5);

        let data = mol_to_binary(&mol).unwrap();
        let mol2 = mol_from_binary(&data).unwrap();
        assert_eq!(mol, mol2, "methane with 2D coords roundtrip failed");
        assert!(mol2.coords_2d().is_some());
    }

    #[test]
    fn test_methane_with_3d_conformer() {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(c));
        for _ in 0..4 {
            builder.add_atom(AtomSpec::new(h));
        }
        for i in 0..4 {
            builder
                .add_bond(BondSpec::new(
                    AtomId::new(0),
                    AtomId::new(i + 1),
                    BondOrder::Single,
                ))
                .unwrap();
        }
        builder
            .add_3d_conformer(vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [-0.5, 0.866, 0.0],
                [-0.5, -0.866, 0.0],
                [0.0, 1.0, 0.0],
            ])
            .unwrap();
        let mol = builder.build().expect("build methane with 3D");
        assert_eq!(mol.conformers_3d().len(), 1);

        let data = mol_to_binary(&mol).unwrap();
        let mol2 = mol_from_binary(&data).unwrap();
        assert_eq!(mol, mol2, "methane with 3D conformer roundtrip failed");
    }

    #[test]
    fn test_ethanol_roundtrip() {
        let mol = build_simple_ethanol();
        assert_eq!(mol.num_atoms(), 10);

        let data = mol_to_binary(&mol).unwrap();
        let mol2 = mol_from_binary(&data).unwrap();
        assert_eq!(mol, mol2, "ethanol roundtrip failed");
    }

    #[test]
    fn test_roundtrip_with_bond_props() {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(c));
        builder.add_atom(AtomSpec::new(h));
        builder
            .add_bond(
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                    .with_prop("wiberg", "0.85")
                    .with_stereo(BondStereo::Z)
                    .with_aromatic(false)
                    .with_conjugated(true),
            )
            .unwrap();
        let mol = builder.build().expect("build molecule");

        let data = mol_to_binary(&mol).unwrap();
        let mol2 = mol_from_binary(&data).unwrap();

        let bond = &mol2.bonds()[0];
        assert_eq!(bond.prop("wiberg"), Some("0.85"));
        assert_eq!(bond.stereo(), BondStereo::Z);
        assert!(!bond.is_aromatic());
        assert!(bond.is_conjugated());
        assert_eq!(mol, mol2, "bond props roundtrip failed");
    }

    #[test]
    fn test_roundtrip_with_atom_props() {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(
            AtomSpec::new(c)
                .with_chiral_tag(ChiralTag::TetrahedralCw)
                .with_chiral_permutation(42)
                .with_isotope(13)
                .with_formal_charge(1)
                .with_radical_electrons(0)
                .with_hybridization(Hybridization::Sp3)
                .with_atom_map(5)
                .with_aromatic(false)
                .with_prop("test_key", "test_val"),
        );
        builder.add_atom(AtomSpec::new(h));
        builder
            .add_bond(BondSpec::new(
                AtomId::new(0),
                AtomId::new(1),
                BondOrder::Single,
            ))
            .unwrap();
        let mol = builder.build().expect("build molecule");

        let data = mol_to_binary(&mol).unwrap();
        let mol2 = mol_from_binary(&data).unwrap();

        let atom = &mol2.atoms()[0];
        assert_eq!(atom.atomic_number(), 6);
        assert_eq!(atom.chiral_tag(), ChiralTag::TetrahedralCw);
        assert_eq!(atom.chiral_permutation(), Some(42));
        assert_eq!(atom.isotope(), Some(13));
        assert_eq!(atom.formal_charge(), 1);
        assert_eq!(atom.radical_electrons(), 0);
        assert_eq!(atom.hybridization(), Hybridization::Sp3);
        assert_eq!(atom.atom_map(), Some(5));
        assert!(!atom.is_aromatic());
        assert_eq!(atom.prop("test_key"), Some("test_val"));
        assert_eq!(mol, mol2, "atom props roundtrip failed");
    }

    #[test]
    fn test_roundtrip_with_sdf_property_lists() {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(c));
        builder.add_atom(AtomSpec::new(h));
        builder
            .add_bond(BondSpec::new(
                AtomId::new(0),
                AtomId::new(1),
                BondOrder::Single,
            ))
            .unwrap();

        let plist = SdfPropertyList::new(
            SdfPropertyListTarget::Atom,
            "test_list",
            vec![Some("val1".into()), Some("val2".into())],
        );

        // Build MoleculeProperties separately since builder doesn't expose with_sdf_property_list
        let base = builder.build().expect("build molecule");

        // Construct the molecule with properties directly via internal API
        let mut mol_props = crate::MoleculeProperties::default();
        mol_props = mol_props.with_sdf_property_list(plist);
        let topology = crate::molecule::TopologyBlock {
            atoms: base.atoms().to_vec(),
            bonds: base.bonds().to_vec(),
            substance_groups: vec![],
            stereo_groups: vec![],
        };
        let coord_block = crate::molecule::CoordinateBlock {
            coords_2d: None,
            conformers_3d: vec![],
            source_coordinate_dim: None,
        };
        let mol = crate::Molecule::from_blocks(topology, coord_block, mol_props)
            .expect("build molecule with property lists");
        assert_eq!(mol.properties().sdf_property_lists().len(), 1);

        let data = mol_to_binary(&mol).unwrap();
        let mol2 = mol_from_binary(&data).unwrap();
        assert_eq!(mol, mol2, "SDF property list roundtrip failed");
    }

    #[test]
    fn test_roundtrip_with_stereo_groups() {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(c));
        builder.add_atom(AtomSpec::new(h));
        builder
            .add_bond(BondSpec::new(
                AtomId::new(0),
                AtomId::new(1),
                BondOrder::Single,
            ))
            .unwrap();

        let sg =
            StereoGroup::new(StereoGroupKind::Absolute, vec![AtomId::new(0)], vec![]).with_id(1);
        builder.add_stereo_group(sg).unwrap();

        let mol = builder.build().expect("build molecule");
        assert_eq!(mol.stereo_groups().len(), 1);

        let data = mol_to_binary(&mol).unwrap();
        let mol2 = mol_from_binary(&data).unwrap();
        assert_eq!(mol, mol2, "stereo group roundtrip failed");
    }

    #[test]
    fn test_invalid_version() {
        let data = vec![0xFF, 0x00, 0x00, 0x00, 0x00];
        let result = mol_from_binary(&data);
        assert!(result.is_err(), "expected error for unsupported version");
        match result {
            Err(PickleError::UnsupportedVersion(v)) => assert_eq!(v, 0xFF),
            _ => panic!("expected UnsupportedVersion error"),
        }
    }

    #[test]
    fn test_truncated_data() {
        let data = vec![0x01];
        let result = mol_from_binary(&data);
        assert!(result.is_err(), "expected error for truncated data");
        match result {
            Err(PickleError::UnexpectedEof) => {}
            _ => panic!("expected UnexpectedEof error"),
        }
    }

    #[test]
    fn test_methane_with_sdf_data_fields() {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(c));
        for _ in 0..4 {
            builder.add_atom(AtomSpec::new(h));
        }
        for i in 0..4 {
            builder
                .add_bond(BondSpec::new(
                    AtomId::new(0),
                    AtomId::new(i + 1),
                    BondOrder::Single,
                ))
                .unwrap();
        }
        builder = builder.with_sdf_data_field("PUBCHEM_IUPAC_NAME", "methane");
        builder = builder.with_sdf_data_field("PUBCHEM_MOLECULAR_FORMULA", "CH4");
        let mol = builder.build().expect("build methane");

        let data = mol_to_binary(&mol).unwrap();
        let mol2 = mol_from_binary(&data).unwrap();
        assert_eq!(mol, mol2, "methane with SDF data fields roundtrip failed");
    }

    #[test]
    fn test_ethanol_with_stereo_atoms() {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = MoleculeBuilder::new();

        // Build ethylene-like with stereo atoms on double bond
        // C=C with stereo
        builder.add_atom(AtomSpec::new(c)); // 0
        builder.add_atom(AtomSpec::new(c)); // 1
        builder.add_atom(AtomSpec::new(h)); // 2
        builder.add_atom(AtomSpec::new(h)); // 3
        builder.add_atom(AtomSpec::new(h)); // 4
        builder.add_atom(AtomSpec::new(h)); // 5

        // Bonds for substituents
        builder
            .add_bond(
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double)
                    .with_stereo(BondStereo::Z)
                    .with_stereo_atoms(AtomId::new(2), AtomId::new(4)),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(
                AtomId::new(0),
                AtomId::new(2),
                BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(BondSpec::new(
                AtomId::new(0),
                AtomId::new(3),
                BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(BondSpec::new(
                AtomId::new(1),
                AtomId::new(4),
                BondOrder::Single,
            ))
            .unwrap();
        builder
            .add_bond(BondSpec::new(
                AtomId::new(1),
                AtomId::new(5),
                BondOrder::Single,
            ))
            .unwrap();

        let mol = builder.build().expect("build ethylene");
        assert_eq!(
            mol.bonds()[0].stereo_atoms(),
            Some([AtomId::new(2), AtomId::new(4)])
        );

        let data = mol_to_binary(&mol).unwrap();
        let mol2 = mol_from_binary(&data).unwrap();
        assert_eq!(mol, mol2, "ethylene roundtrip failed");
        assert_eq!(
            mol2.bonds()[0].stereo_atoms(),
            Some([AtomId::new(2), AtomId::new(4)])
        );
    }

    #[test]
    fn test_multiple_conformers() {
        let c = Element::from_atomic_number(6).unwrap();
        let h = Element::from_atomic_number(1).unwrap();
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(c));
        builder.add_atom(AtomSpec::new(h));
        builder
            .add_bond(BondSpec::new(
                AtomId::new(0),
                AtomId::new(1),
                BondOrder::Single,
            ))
            .unwrap();

        builder
            .add_3d_conformer(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
            .unwrap();
        builder
            .add_3d_conformer(vec![[0.5, 0.5, 0.5], [1.5, 0.5, 0.5]])
            .unwrap();

        let mol = builder.build().expect("build with conformers");
        assert_eq!(mol.conformers_3d().len(), 2);

        let data = mol_to_binary(&mol).unwrap();
        let mol2 = mol_from_binary(&data).unwrap();
        assert_eq!(mol, mol2, "multiple conformers roundtrip failed");
    }
}
