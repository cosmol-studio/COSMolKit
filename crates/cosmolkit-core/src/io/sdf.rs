use std::io::BufRead;

use crate::Molecule;

/// One record extracted from an SDF stream.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfRecord {
    pub molecule: Molecule,
    pub data_fields: Vec<(String, String)>,
}

/// Errors returned by SDF reading APIs.
#[derive(Debug)]
pub enum SdfReadError {
    Io(std::io::Error),
    Parse(String),
    NotImplemented,
}

impl From<std::io::Error> for SdfReadError {
    fn from(value: std::io::Error) -> Self {
        Self::Io(value)
    }
}

/// Streaming-capable SDF reader.
pub struct SdfReader<R> {
    reader: R,
}

impl<R: BufRead> SdfReader<R> {
    /// Create a reader from any buffered input stream.
    #[must_use]
    pub fn new(reader: R) -> Self {
        Self { reader }
    }

    /// Read next SDF record from stream.
    pub fn next_record(&mut self) -> Result<Option<SdfRecord>, SdfReadError> {
        let _ = &mut self.reader;
        unimplemented!("SDF streaming reader is not implemented yet")
    }
}
