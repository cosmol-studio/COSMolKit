//! Small manifest preflight and a single streaming checksum/decompression pass.
use crate::Result;
use flate2::{Compression, read::ZlibDecoder, write::ZlibEncoder};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    fs::File,
    io::{BufReader, Cursor, Read, Write},
    path::Path,
};

pub const MAGIC: &[u8; 8] = b"CKREF001";
pub const MAX_FRAME_BYTES: usize = 128 * 1024 * 1024;
pub const MAX_RECORD_BYTES: usize = MAX_FRAME_BYTES - 4;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Identity {
    pub rdkit: String,
    pub task: String,
    pub profiles: serde_json::Value,
    pub input_sha256: String,
    pub records: u64,
}

impl Identity {
    pub fn key(&self) -> Result<String> {
        Ok(hex(&Sha256::digest(serde_json::to_vec(self)?)))
    }
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Manifest {
    pub format: u32,
    pub identity: Identity,
    pub frames: u64,
    pub bytes: u64,
    pub sha256: String,
}

pub fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

pub fn frame(records: &[Vec<u8>]) -> Result<Vec<u8>> {
    if records.is_empty() {
        return Err("empty reference frame".into());
    }
    let mut raw = Vec::new();
    for record in records {
        if record.len() > MAX_RECORD_BYTES || raw.len() + 4 + record.len() > MAX_FRAME_BYTES {
            return Err("reference exceeds the bounded frame/record budget".into());
        }
        raw.extend((record.len() as u32).to_le_bytes());
        raw.extend(record);
    }
    let mut encoder = ZlibEncoder::new(Vec::new(), Compression::new(1));
    encoder.write_all(&raw)?;
    let compressed = encoder.finish()?;
    if compressed.len() > MAX_FRAME_BYTES {
        return Err("compressed reference frame exceeds budget".into());
    }
    let mut output = Vec::with_capacity(44 + compressed.len());
    output.extend((compressed.len() as u32).to_le_bytes());
    output.extend((raw.len() as u32).to_le_bytes());
    output.extend(u32::try_from(records.len())?.to_le_bytes());
    output.extend(Sha256::digest(&compressed));
    output.extend(compressed);
    Ok(output)
}

pub struct References {
    reader: BufReader<File>,
    manifest: Manifest,
    digest: Sha256,
    frames: u64,
    records: u64,
}

impl References {
    /// No full-reference copy, decode, hash or Python process during preflight.
    pub fn open(folder: &Path, identity: &Identity) -> Result<Self> {
        let metadata = File::open(folder.join("manifest.json"))?;
        let manifest: Manifest = serde_json::from_reader(metadata.take(4 * 1024 * 1024))?;
        if manifest.format != 1 || manifest.identity != *identity {
            return Err("reference identity changed; run prepare for this selection".into());
        }
        let file = File::open(folder.join("reference.bin"))?;
        if file.metadata()?.len() != manifest.bytes || manifest.frames == 0 || identity.records == 0
        {
            return Err("incomplete reference; run prepare for this selection".into());
        }
        let mut reader = BufReader::new(file);
        let mut magic = [0; 8];
        reader.read_exact(&mut magic)?;
        if &magic != MAGIC {
            return Err("unknown reference format".into());
        }
        let mut digest = Sha256::new();
        digest.update(magic);
        Ok(Self {
            reader,
            manifest,
            digest,
            frames: 0,
            records: 0,
        })
    }

    pub fn next_frame(&mut self) -> Result<Option<Vec<Vec<u8>>>> {
        if self.frames == self.manifest.frames {
            if self.records != self.manifest.identity.records
                || hex(&self.digest.clone().finalize()) != self.manifest.sha256
            {
                return Err("reference checksum/record count mismatch".into());
            }
            let mut extra = [0];
            if self.reader.read(&mut extra)? != 0 {
                return Err("trailing reference data".into());
            }
            return Ok(None);
        }
        let mut header = [0; 44];
        self.reader.read_exact(&mut header)?;
        let stored = u32::from_le_bytes(header[0..4].try_into()?) as usize;
        let expanded = u32::from_le_bytes(header[4..8].try_into()?) as usize;
        let count = u32::from_le_bytes(header[8..12].try_into()?) as usize;
        if stored > MAX_FRAME_BYTES
            || expanded > MAX_FRAME_BYTES
            || count == 0
            || count > expanded / 4
        {
            return Err("invalid reference frame allocation/count".into());
        }
        let mut compressed = vec![0; stored];
        self.reader.read_exact(&mut compressed)?;
        if Sha256::digest(&compressed).as_slice() != &header[12..44] {
            return Err("reference frame checksum mismatch".into());
        }
        self.digest.update(header);
        self.digest.update(&compressed);
        let mut raw = Vec::with_capacity(expanded);
        let mut decoder = ZlibDecoder::new(compressed.as_slice());
        (&mut decoder)
            .take(expanded as u64 + 1)
            .read_to_end(&mut raw)?;
        if raw.len() != expanded || decoder.total_in() != stored as u64 {
            return Err("reference decompression length mismatch".into());
        }
        let mut cursor = Cursor::new(raw.as_slice());
        let mut records = Vec::with_capacity(count);
        for _ in 0..count {
            let mut length = [0; 4];
            cursor.read_exact(&mut length)?;
            let length = u32::from_le_bytes(length) as usize;
            if length > MAX_RECORD_BYTES || length > raw.len() - cursor.position() as usize {
                return Err("invalid reference record length".into());
            }
            let mut record = vec![0; length];
            cursor.read_exact(&mut record)?;
            records.push(record);
        }
        if cursor.position() != expanded as u64 {
            return Err("trailing reference frame data".into());
        }
        self.frames += 1;
        self.records += count as u64;
        Ok(Some(records))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn identity() -> Identity {
        Identity {
            rdkit: "2026.03.6".into(),
            task: "example".into(),
            profiles: serde_json::json!([{}]),
            input_sha256: "input".into(),
            records: 2,
        }
    }
    fn fixture() -> tempfile::TempDir {
        let dir = tempfile::tempdir().unwrap();
        let mut data = MAGIC.to_vec();
        data.extend(frame(&[vec![0, 1], vec![2, 3]]).unwrap());
        let manifest = Manifest {
            format: 1,
            identity: identity(),
            frames: 1,
            bytes: data.len() as u64,
            sha256: hex(&Sha256::digest(&data)),
        };
        std::fs::write(dir.path().join("reference.bin"), data).unwrap();
        std::fs::write(
            dir.path().join("manifest.json"),
            serde_json::to_vec(&manifest).unwrap(),
        )
        .unwrap();
        dir
    }

    #[test]
    fn references_stream_and_validate_once() {
        let dir = fixture();
        let mut reader = References::open(dir.path(), &identity()).unwrap();
        assert_eq!(
            reader.next_frame().unwrap().unwrap(),
            [vec![0, 1], vec![2, 3]]
        );
        assert!(reader.next_frame().unwrap().is_none());
    }

    #[test]
    fn cache_identity_changes_only_with_reference_options_or_inputs() {
        let original = identity();
        let mut changed = original.clone();
        changed.rdkit.push('x');
        assert_ne!(original.key().unwrap(), changed.key().unwrap());
        changed = original.clone();
        changed.profiles = serde_json::json!([{"radius": 3}]);
        assert_ne!(original.key().unwrap(), changed.key().unwrap());
        changed = original.clone();
        changed.input_sha256.push('x');
        assert_ne!(original.key().unwrap(), changed.key().unwrap());
        assert_eq!(original.key().unwrap(), original.clone().key().unwrap());
    }

    #[test]
    fn stale_missing_and_corrupt_references_fail() {
        let dir = fixture();
        let mut stale = identity();
        stale.records = 3;
        assert!(References::open(dir.path(), &stale).is_err());
        let path = dir.path().join("reference.bin");
        let mut data = std::fs::read(&path).unwrap();
        let last = data.len() - 1;
        data[last] ^= 1;
        std::fs::write(path, data).unwrap();
        let mut reader = References::open(dir.path(), &identity()).unwrap();
        assert!(reader.next_frame().is_err());
        assert!(References::open(&dir.path().join("missing"), &identity()).is_err());
    }

    #[test]
    fn hostile_lengths_fail_before_allocation() {
        let dir = fixture();
        let path = dir.path().join("reference.bin");
        let mut data = std::fs::read(&path).unwrap();
        data[8..12].copy_from_slice(&u32::MAX.to_le_bytes());
        std::fs::write(path, data).unwrap();
        let mut reader = References::open(dir.path(), &identity()).unwrap();
        assert!(reader.next_frame().is_err());
        assert!(frame(&[]).is_err());
    }
}
