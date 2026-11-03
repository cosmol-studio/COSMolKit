//! Same-version viewer transfer. Wire records are private IO data, not a second BIO model.
//! Preserve hierarchy, atom display values, coordinate bits and structural annotations.
//! Original input text, refinement/publication metadata and crystal/assembly state are omitted.
use cosmolkit_bio::*;
use cosmolkit_types::Element;
use musli::{Decode, Encode};

const MAGIC: &[u8; 8] = b"CKPROTX1";
const VERSION: &str = env!("CARGO_PKG_VERSION");
const MAX_BYTES: usize = 256 * 1024 * 1024;
const MAX_ROWS: usize = 1_000_000;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum ProteinTransferError {
    #[error("protein transfer requires COSMolKit {expected}, but the sender uses {actual}")]
    VersionMismatch { expected: String, actual: String },
    #[error("invalid protein transfer: {message}")]
    InvalidData { message: String },
    #[error("invalid transferred protein: {0}")]
    Projection(#[from] ProteinProjectionError),
}

impl From<BioStructureError> for ProteinTransferError {
    fn from(error: BioStructureError) -> Self {
        ProteinProjectionError::Structure(error).into()
    }
}

fn bad(message: impl Into<String>) -> ProteinTransferError {
    ProteinTransferError::InvalidData {
        message: message.into(),
    }
}

// Closed vocabulary conversion is exhaustive on the native side; a new native
// variant cannot silently turn into a default or an unchecked discriminant.
macro_rules! vocabulary {
    ($wire:ident, $native:ident, $($variant:ident),+ $(,)?) => {
        #[derive(Encode, Decode)]
        enum $wire { $($variant),+ }
        impl From<$native> for $wire {
            fn from(value: $native) -> Self {
                match value { $($native::$variant => Self::$variant),+ }
            }
        }
        impl From<$wire> for $native {
            fn from(value: $wire) -> Self {
                match value { $($wire::$variant => Self::$variant),+ }
            }
        }
    };
}
vocabulary!(
    Format,
    BioCoordinateFormat,
    Unknown,
    Detect,
    Pdb,
    Mmcif,
    Mmjson,
    ChemComp
);
vocabulary!(
    ChainType,
    ChainKind,
    Protein,
    Dna,
    Rna,
    ProteinDnaComplex,
    ProteinRnaComplex,
    LigandOnly,
    WaterOnly,
    Mixed,
    Unknown
);
vocabulary!(
    ResidueType,
    ResidueInfoKind,
    Unknown,
    Aa,
    Aad,
    Paa,
    Maa,
    Rna,
    Dna,
    Buf,
    Hoh,
    Pyr,
    Ket,
    Els
);
vocabulary!(
    EntityType, EntityKind, Unknown, Polymer, NonPolymer, Branched, Water
);
vocabulary!(
    PolymerType,
    PolymerKind,
    Unknown,
    PeptideL,
    PeptideD,
    Dna,
    Rna,
    DnaRnaHybrid,
    SaccharideD,
    SaccharideL,
    Pna,
    CyclicPseudoPeptide,
    Other
);
vocabulary!(
    CalcFlag,
    BioCalcFlag,
    NotSet,
    NoHydrogen,
    Determined,
    Calculated,
    Dummy
);
vocabulary!(
    ConnectionType,
    BioConnectionKind,
    Covale,
    Disulf,
    Hydrog,
    MetalC,
    Unknown
);
vocabulary!(Asu, BioAsu, Same, Different, Any);
vocabulary!(
    HelixClass,
    BioHelixClass,
    UnknownHelix,
    RAlpha,
    ROmega,
    RPi,
    RGamma,
    R310,
    LAlpha,
    LOmega,
    LGamma,
    Helix27,
    HelixPolyProlineNone
);

/// Fixed-capacity source names avoid an allocation for every atom/residue name.
#[derive(Encode, Decode)]
struct ShortText {
    bytes: [u8; 4],
    len: u8,
}
impl ShortText {
    fn new(bytes: &[u8]) -> Self {
        // Callers supply only existing validated four-byte BIO names.
        let mut stored = [0; 4];
        stored[..bytes.len()].copy_from_slice(bytes);
        Self {
            bytes: stored,
            len: bytes.len() as u8,
        }
    }
    fn bytes(&self) -> Result<&[u8], ProteinTransferError> {
        self.bytes
            .get(..usize::from(self.len))
            .filter(|bytes| bytes.is_ascii())
            .ok_or_else(|| bad("invalid bounded source name"))
    }
}

#[derive(Encode, Decode)]
struct Model {
    start: u32,
    len: u32,
    number: Option<i32>,
}
#[derive(Encode, Decode)]
struct Chain {
    model: u32,
    entity: Option<u32>,
    start: u32,
    len: u32,
    kind: ChainType,
    auth_id: Option<ShortText>,
    label_id: Option<String>,
}
#[derive(Encode, Decode)]
struct Residue {
    chain: u32,
    start: u32,
    len: u32,
    name: ShortText,
    kind: ResidueType,
    entity_kind: EntityType,
    entity: Option<u32>,
    het_flag: Option<u8>,
    sequence: Option<(i32, Option<u8>)>,
    label_sequence: Option<i32>,
    segment: Option<[u8; 4]>,
    subchain: Option<String>,
    label_entity: Option<String>,
}
#[derive(Encode, Decode)]
struct Atom {
    residue: u32,
    name: ShortText,
    element: u8,
    isotope: Option<u16>,
    altloc: Option<u8>,
    charge: i8,
    calc_flag: CalcFlag,
    occupancy: f64,
    b_iso: f64,
    anisou: [f64; 6],
    tls_group: i16,
    fraction: f64,
    serial: Option<i32>,
}
#[derive(Encode, Decode)]
struct Entity {
    kind: EntityType,
    polymer: PolymerType,
    microhetero: bool,
    sequence: Vec<String>,
    subchains: Vec<String>,
    source: String,
}
#[derive(Encode, Decode)]
struct Address {
    chain: ShortText,
    sequence: Option<i32>,
    insertion: Option<u8>,
    segment: ShortText,
    residue: ShortText,
    atom: String,
    altloc: u8,
}
impl From<&AtomAddress> for Address {
    fn from(value: &AtomAddress) -> Self {
        let residue = value.residue();
        Self {
            chain: ShortText::new(value.chain_name().as_bytes()),
            sequence: residue.sequence_number(),
            insertion: residue.insertion_code(),
            segment: ShortText::new(residue.segment().as_bytes()),
            residue: ShortText::new(residue.name().as_bytes()),
            atom: value.logical_atom_name().into(),
            altloc: value.altloc(),
        }
    }
}
impl Address {
    fn restore(self) -> Result<AtomAddress, ProteinTransferError> {
        let chain =
            PdbChainId::from_ascii(self.chain.bytes()?).ok_or_else(|| bad("invalid chain name"))?;
        let name = ResidueName::from_ascii(self.residue.bytes()?)
            .ok_or_else(|| bad("invalid residue name"))?;
        let residue =
            ResidueAddress::new(self.sequence, self.insertion, self.segment.bytes()?, name)
                .ok_or_else(|| bad("invalid residue address"))?;
        Ok(AtomAddress::new(
            chain,
            residue,
            self.atom,
            Some(self.altloc),
        ))
    }
}
#[derive(Encode, Decode)]
struct Connection {
    name: String,
    link: String,
    kind: ConnectionType,
    asu: Asu,
    first: Address,
    second: Address,
    distance: f64,
    symmetry: [i16; 4],
}
#[derive(Encode, Decode)]
struct CisPeptide {
    first: Address,
    second: Address,
    model: i32,
    altloc: u8,
    angle: f64,
}
#[derive(Encode, Decode)]
struct Modification {
    address: Address,
    parent: String,
    id: String,
    details: String,
}
#[derive(Encode, Decode)]
struct Helix {
    start: Address,
    end: Address,
    class: HelixClass,
    length: i32,
}
#[derive(Encode, Decode)]
struct Strand {
    start: Address,
    end: Address,
    hbond2: Address,
    hbond1: Address,
    sense: i32,
    name: String,
}
#[derive(Encode, Decode)]
struct Sheet {
    name: String,
    strands: Vec<Strand>,
}
#[derive(Encode, Decode)]
struct ProteinState {
    format: Format,
    name: String,
    models: Vec<Model>,
    chains: Vec<Chain>,
    residues: Vec<Residue>,
    atoms: Vec<Atom>,
    coordinates: Vec<[f64; 3]>,
    entities: Vec<Entity>,
    connections: Vec<Connection>,
    cispeps: Vec<CisPeptide>,
    modifications: Vec<Modification>,
    helices: Vec<Helix>,
    sheets: Vec<Sheet>,
    serial_connections: std::collections::BTreeMap<i32, Vec<i32>>,
}

impl ProteinState {
    fn capture(data: &BioStructureData) -> Self {
        Self {
            format: data.input_format().into(),
            name: data.source_state().name.clone(),
            serial_connections: data.source_state().conect_map.clone(),
            models: data
                .models()
                .iter()
                .map(|row| Model {
                    start: row.chain_span().start(),
                    len: row.chain_span().len(),
                    number: row.source_model_number(),
                })
                .collect(),
            chains: data
                .chains()
                .iter()
                .map(|row| Chain {
                    model: row.model_id().value(),
                    entity: row.entity_id().map(BioEntityId::value),
                    start: row.residue_span().start(),
                    len: row.residue_span().len(),
                    kind: row.kind().into(),
                    auth_id: row
                        .source()
                        .auth_chain_id()
                        .map(|id| ShortText::new(id.as_bytes())),
                    label_id: row.source().label_asym_id().map(str::to_owned),
                })
                .collect(),
            residues: data
                .residues()
                .iter()
                .map(|row| Residue {
                    chain: row.chain_id().value(),
                    start: row.atom_span().start(),
                    len: row.atom_span().len(),
                    name: ShortText::new(row.name().as_bytes()),
                    kind: row.residue_info_kind().into(),
                    entity_kind: row.entity_kind().into(),
                    entity: row.entity_id().map(BioEntityId::value),
                    het_flag: row.het_flag(),
                    sequence: row
                        .source()
                        .seq_id()
                        .map(|id| (id.seq_num(), id.ins_code())),
                    label_sequence: row.source().label_seq_id(),
                    segment: row.source().segment_id().copied(),
                    subchain: row.source().subchain_id().map(str::to_owned),
                    label_entity: row.source().label_entity_id().map(str::to_owned),
                })
                .collect(),
            atoms: data
                .atoms()
                .iter()
                .map(|row| Atom {
                    residue: row.residue_id().value(),
                    name: ShortText::new(row.name().as_bytes()),
                    element: row.element().atomic_number(),
                    isotope: row.isotope_mass_number(),
                    altloc: row.altloc().map(AltLocLabel::value),
                    charge: row.formal_charge(),
                    calc_flag: row.calc_flag().into(),
                    occupancy: row.occupancy(),
                    b_iso: row.b_iso(),
                    anisou: *row.anisou(),
                    tls_group: row.tls_group_id(),
                    fraction: row.fraction(),
                    serial: row.source().serial().map(PdbAtomSerial::value),
                })
                .collect(),
            coordinates: data.coordinates().positions().to_vec(),
            entities: data
                .entities()
                .iter()
                .map(|row| Entity {
                    kind: row.kind().into(),
                    polymer: row.polymer_kind().into(),
                    microhetero: row.reflects_microhetero(),
                    sequence: row.full_sequence().to_vec(),
                    subchains: row.subchains().to_vec(),
                    source: row.source().source_entity_id().into(),
                })
                .collect(),
            connections: data
                .connections()
                .iter()
                .map(|row| Connection {
                    name: row.name.clone(),
                    link: row.link_id.clone(),
                    kind: row.kind.into(),
                    asu: row.asu.into(),
                    first: (&row.partner1).into(),
                    second: (&row.partner2).into(),
                    distance: row.reported_distance,
                    symmetry: row.reported_sym,
                })
                .collect(),
            cispeps: data
                .cispeps()
                .iter()
                .map(|row| CisPeptide {
                    first: (&row.partner_c).into(),
                    second: (&row.partner_n).into(),
                    model: row.model_num,
                    altloc: row.only_altloc,
                    angle: row.reported_angle,
                })
                .collect(),
            modifications: data
                .mod_residues()
                .iter()
                .map(|row| Modification {
                    address: (&AtomAddress::new(row.chain_name, row.res_id, "", None)).into(),
                    parent: row.parent_comp_id.clone(),
                    id: row.mod_id.clone(),
                    details: row.details.clone(),
                })
                .collect(),
            helices: data
                .helices()
                .iter()
                .map(|row| Helix {
                    start: (&row.start).into(),
                    end: (&row.end).into(),
                    class: row.pdb_helix_class.into(),
                    length: row.length,
                })
                .collect(),
            sheets: data
                .sheets()
                .iter()
                .map(|row| Sheet {
                    name: row.name.clone(),
                    strands: row
                        .strands
                        .iter()
                        .map(|row| Strand {
                            start: (&row.start).into(),
                            end: (&row.end).into(),
                            hbond2: (&row.hbond_atom2).into(),
                            hbond1: (&row.hbond_atom1).into(),
                            sense: row.sense,
                            name: row.name.clone(),
                        })
                        .collect(),
                })
                .collect(),
        }
    }

    fn restore(self) -> Result<BioStructureData, ProteinTransferError> {
        check_counts([
            self.models.len(),
            self.chains.len(),
            self.residues.len(),
            self.atoms.len(),
            self.coordinates.len(),
            self.entities.len(),
            self.connections.len(),
            self.cispeps.len(),
            self.modifications.len(),
            self.helices.len(),
            self.sheets.len(),
        ])?;
        let parts = BioStructureParts {
            input_format: self.format.into(),
            models: self
                .models
                .into_iter()
                .map(|row| {
                    Ok(BioModelRow::new(
                        BioRowSpan::new(row.start, row.len)?,
                        row.number,
                    ))
                })
                .collect::<Result<_, ProteinTransferError>>()?,
            chains: self
                .chains
                .into_iter()
                .map(|row| {
                    let auth = row
                        .auth_id
                        .map(|name| {
                            PdbChainId::from_ascii(name.bytes()?)
                                .ok_or_else(|| bad("invalid chain name"))
                        })
                        .transpose()?;
                    Ok(BioChainRow::new(
                        BioModelId::new(row.model),
                        row.entity.map(BioEntityId::new),
                        BioRowSpan::new(row.start, row.len)?,
                        row.kind.into(),
                        ChainSourceIds::new(auth, row.label_id),
                    ))
                })
                .collect::<Result<_, ProteinTransferError>>()?,
            residues: self
                .residues
                .into_iter()
                .map(|row| {
                    let name = ResidueName::from_ascii(row.name.bytes()?)
                        .ok_or_else(|| bad("invalid residue name"))?;
                    let source = ResidueSourceIds::new(
                        row.sequence
                            .map(|(number, code)| PdbSeqId::new(number, code)),
                        row.label_sequence,
                        row.segment,
                        row.subchain,
                        row.label_entity,
                    )
                    .ok_or_else(|| bad("invalid residue source identity"))?;
                    Ok(BioResidueRow::new(
                        BioChainId::new(row.chain),
                        BioRowSpan::new(row.start, row.len)?,
                        name,
                        row.kind.into(),
                        row.entity_kind.into(),
                        row.entity.map(BioEntityId::new),
                        row.het_flag,
                        source,
                        BioSiftsUnpResidue::default(),
                    ))
                })
                .collect::<Result<_, ProteinTransferError>>()?,
            atoms: self
                .atoms
                .into_iter()
                .map(|row| {
                    let name = AtomName::from_ascii(row.name.bytes()?)
                        .ok_or_else(|| bad("invalid atom name"))?;
                    let element = Element::from_atomic_number(row.element)
                        .ok_or_else(|| bad("invalid element"))?;
                    let altloc = row.altloc.map(AltLocLabel::new);
                    Ok(BioAtomRow::new(
                        BioResidueId::new(row.residue),
                        name,
                        element,
                        row.isotope,
                        altloc,
                        row.charge,
                        row.calc_flag.into(),
                        row.occupancy,
                        row.b_iso,
                        row.anisou,
                        row.tls_group,
                        row.fraction,
                        AtomSourceIds::new(row.serial.map(PdbAtomSerial::new)),
                    ))
                })
                .collect::<Result<_, ProteinTransferError>>()?,
            coordinates: BioCoordinateBlock::new(self.coordinates),
            entities: self
                .entities
                .into_iter()
                .map(|row| {
                    BioEntityRow::new(
                        row.kind.into(),
                        row.polymer.into(),
                        row.microhetero,
                        row.sequence,
                        Vec::new(),
                        Vec::new(),
                        row.subchains,
                        EntitySourceIds::new(row.source),
                    )
                })
                .collect(),
            connections: self
                .connections
                .into_iter()
                .map(|row| {
                    Ok(BioConnection {
                        name: row.name,
                        link_id: row.link,
                        kind: row.kind.into(),
                        asu: row.asu.into(),
                        partner1: row.first.restore()?,
                        partner2: row.second.restore()?,
                        reported_distance: row.distance,
                        reported_sym: row.symmetry,
                    })
                })
                .collect::<Result<_, ProteinTransferError>>()?,
            cispeps: self
                .cispeps
                .into_iter()
                .map(|row| {
                    Ok(BioCisPep {
                        partner_c: row.first.restore()?,
                        partner_n: row.second.restore()?,
                        model_num: row.model,
                        only_altloc: row.altloc,
                        reported_angle: row.angle,
                    })
                })
                .collect::<Result<_, ProteinTransferError>>()?,
            mod_residues: self
                .modifications
                .into_iter()
                .map(|row| {
                    let address = row.address.restore()?;
                    Ok(BioModRes {
                        chain_name: address.chain_name(),
                        res_id: address.residue(),
                        parent_comp_id: row.parent,
                        mod_id: row.id,
                        details: row.details,
                    })
                })
                .collect::<Result<_, ProteinTransferError>>()?,
            helices: self
                .helices
                .into_iter()
                .map(|row| {
                    Ok(BioHelix {
                        start: row.start.restore()?,
                        end: row.end.restore()?,
                        pdb_helix_class: row.class.into(),
                        length: row.length,
                    })
                })
                .collect::<Result<_, ProteinTransferError>>()?,
            sheets: self
                .sheets
                .into_iter()
                .map(|row| {
                    Ok(BioSheet {
                        name: row.name,
                        strands: row
                            .strands
                            .into_iter()
                            .map(|row| {
                                Ok(BioStrand::new(
                                    row.start.restore()?,
                                    row.end.restore()?,
                                    row.hbond2.restore()?,
                                    row.hbond1.restore()?,
                                    row.sense,
                                    row.name,
                                ))
                            })
                            .collect::<Result<_, ProteinTransferError>>()?,
                    })
                })
                .collect::<Result<_, ProteinTransferError>>()?,
            metadata: BioMetadata::default(),
            source_state: BioStructureSourceState {
                name: self.name,
                conect_map: self.serial_connections,
                ..Default::default()
            },
            crystal: None,
            ncs_operators: Vec::new(),
            assemblies: Vec::new(),
        };
        let data = BioStructureData::from_parts(parts)?;
        validate_protein_structure(&data)?;
        Ok(data)
    }
}

fn check_counts(counts: impl IntoIterator<Item = usize>) -> Result<(), ProteinTransferError> {
    if counts.into_iter().any(|count| count > MAX_ROWS) {
        return Err(bad("protein transfer table exceeds 1000000 rows"));
    }
    Ok(())
}

/// Internal detached value encoder; no pointers, runtime handles or text-format conversion.
#[doc(hidden)]
pub fn encode_protein_transfer(data: &BioStructureData) -> Result<Vec<u8>, ProteinTransferError> {
    validate_protein_structure(data)?;
    check_counts([
        data.models().len(),
        data.chains().len(),
        data.residues().len(),
        data.atoms().len(),
        data.entities().len(),
        data.connections().len(),
        data.cispeps().len(),
        data.mod_residues().len(),
        data.helices().len(),
        data.sheets().len(),
    ])?;
    let state = ProteinState::capture(data);
    let mut output = Vec::new();
    output.extend_from_slice(MAGIC);
    let length =
        u16::try_from(VERSION.len()).map_err(|_| bad("version exceeds header capacity"))?;
    output.extend_from_slice(&length.to_le_bytes());
    output.extend_from_slice(VERSION.as_bytes());
    let length_offset = output.len();
    output.extend_from_slice(&[0; 4]);
    let start = output.len();
    musli::wire::to_writer(&mut output, &state)
        .map_err(|error| bad(format!("wire encode: {error}")))?;
    if output.len() > MAX_BYTES {
        return Err(bad("protein transfer exceeds 256 MiB"));
    }
    let length = (output.len() - start) as u32;
    output[length_offset..start].copy_from_slice(&length.to_le_bytes());
    Ok(output)
}

/// Reject the sender version before decoding any structural data.
#[doc(hidden)]
pub fn decode_protein_transfer(data: &[u8]) -> Result<BioStructureData, ProteinTransferError> {
    if data.len() > MAX_BYTES {
        return Err(bad("protein transfer exceeds 256 MiB"));
    }
    let mut input = data;
    fn take<'a>(input: &mut &'a [u8], count: usize) -> Result<&'a [u8], ProteinTransferError> {
        let (value, rest) = input
            .split_at_checked(count)
            .ok_or_else(|| bad("truncated protein transfer"))?;
        *input = rest;
        Ok(value)
    }
    if take(&mut input, 8)? != MAGIC {
        return Err(bad("not a COSMolKit protein transfer"));
    }
    let length = u16::from_le_bytes(take(&mut input, 2)?.try_into().unwrap()) as usize;
    let version =
        std::str::from_utf8(take(&mut input, length)?).map_err(|_| bad("version is not UTF-8"))?;
    if version != VERSION {
        return Err(ProteinTransferError::VersionMismatch {
            expected: VERSION.into(),
            actual: version.into(),
        });
    }
    let length = u32::from_le_bytes(take(&mut input, 4)?.try_into().unwrap()) as usize;
    let payload = take(&mut input, length)?;
    if !input.is_empty() {
        return Err(bad("trailing data after protein transfer"));
    }
    let state: ProteinState =
        musli::wire::from_slice(payload).map_err(|error| bad(format!("wire decode: {error}")))?;
    state.restore()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn restored_viewer_data_must_satisfy_structure_and_protein_invariants() {
        let data = crate::read_pdb_bio_structure(
            "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nEND\n",
            "<transfer test>", &crate::BioPdbReadParams::default(),
        ).unwrap();
        let mut state = ProteinState::capture(&data);
        state.residues[0].kind = ResidueType::Hoh;
        assert!(matches!(
            state.restore(),
            Err(ProteinTransferError::Projection(
                ProteinProjectionError::NonAminoAcidResidue { index: 0 }
            ))
        ));
        let mut state = ProteinState::capture(&data);
        state.models[0].start = 1;
        assert!(matches!(
            state.restore(),
            Err(ProteinTransferError::Projection(
                ProteinProjectionError::Structure(_)
            ))
        ));
    }
}
