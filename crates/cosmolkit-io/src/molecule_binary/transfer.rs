//! Same-version value transfer, sharing the archive's molecular wire schema.
//! No runtime layout, pointer, derived cache, or historical decoder is on this path.
use super::*;

const TRANSFER_MAGIC: &[u8; 8] = b"CKMOLTX1";
const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Failure to transfer a molecule between independently compiled CK modules.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum MoleculeTransferError {
    /// The sender and receiver must have exactly the same complete CK version.
    #[error("molecule transfer requires COSMolKit {expected}, but the sender uses {actual}")]
    VersionMismatch { expected: String, actual: String },
    /// Retains the precise wire/structural validation error and its cause.
    #[error("invalid molecule transfer: {0}")]
    Codec(#[from] PickleError),
}

fn bad(message: &str) -> MoleculeTransferError {
    PickleError::InvalidArchive(message.into()).into()
}

/// Encode molecular state and all coordinate sets, without derived caches.
pub fn encode_molecule_transfer(input: &BinaryInput<'_>) -> Result<Vec<u8>, MoleculeTransferError> {
    input.validate()?;
    rows(input.num_atoms())?;
    rows(input.bonds().len())?;
    let molecule = MoleculeState::from_input(input)?;
    let version_length = u16::try_from(VERSION.len())
        .map_err(|_| bad("transfer version exceeds header capacity"))?;
    let mut output = Vec::new();
    output.extend_from_slice(TRANSFER_MAGIC);
    write_u16_le(&mut output, version_length);
    output.extend_from_slice(VERSION.as_bytes());
    let length_offset = output.len();
    write_u32_le(&mut output, 0);
    let payload_start = output.len();
    // Encode into the final buffer, not a second complete payload allocation.
    musli::wire::to_writer(&mut output, &molecule)
        .map_err(|e| PickleError::InvalidArchive(format!("transfer encode: {e}")))?;
    if output.len() > MAX_ARCHIVE_BYTES {
        return Err(bad("molecule transfer exceeds 256 MiB"));
    }
    let payload_length = (output.len() - payload_start) as u32;
    output[length_offset..payload_start].copy_from_slice(&payload_length.to_le_bytes());
    Ok(output)
}

/// Check the complete sender version before decoding any molecular state.
pub fn decode_molecule_transfer(data: &[u8]) -> Result<BinaryRecord, MoleculeTransferError> {
    if data.len() > MAX_ARCHIVE_BYTES {
        return Err(bad("molecule transfer exceeds 256 MiB"));
    }
    let mut reader = PickleReader::new(data);
    if reader.read_exact_slice(TRANSFER_MAGIC.len())? != TRANSFER_MAGIC {
        return Err(bad("not a COSMolKit molecule transfer"));
    }
    let version_length = usize::from(read_u16_le(&mut reader)?);
    let version = std::str::from_utf8(reader.read_exact_slice(version_length)?)
        .map_err(|_| bad("transfer version is not UTF-8"))?;
    if version != VERSION {
        return Err(MoleculeTransferError::VersionMismatch {
            expected: VERSION.into(),
            actual: version.into(),
        });
    }
    let payload_length = read_u32_le(&mut reader)? as usize;
    let payload = reader.read_exact_slice(payload_length)?;
    if reader.remaining() != 0 {
        return Err(bad("trailing data after molecule transfer"));
    }
    let molecule: MoleculeState = musli::wire::from_slice(payload)
        .map_err(|e| PickleError::InvalidArchive(format!("transfer decode: {e}")))?;
    Ok(molecule.into_record(DerivedState {
        valid_bits: 0,
        rings: None,
        ring_families: None,
        valence: None,
    })?)
}
