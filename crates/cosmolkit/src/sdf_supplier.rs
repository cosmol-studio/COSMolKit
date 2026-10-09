//! Thin live-graph projection of the IO owner's framed and indexed records.
use crate::{SdfError, SdfReadParams, SdfRecord};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub use cosmolkit_io::SdfRecordMetadata;

pub(crate) fn data_params(params: &SdfReadParams) -> cosmolkit_io::SdfDataReadParams {
    cosmolkit_io::SdfDataReadParams {
        strict_parsing: params.strict_parsing,
        process_property_lists: params.process_property_lists,
        coordinate_mode: params.coordinate_mode,
        mol_post: Some(cosmolkit_io::MolPostParams {
            sanitize: params.sanitize,
            remove_hs: params.remove_hs,
            expand_attachment_points: params.expand_attachment_points,
        }),
    }
}

pub struct SdfRecordStream<R> {
    inner: cosmolkit_io::SdfGraphReader<R>,
    params: SdfReadParams,
}
impl<R: BufRead> SdfRecordStream<R> {
    pub fn new(reader: R) -> Self {
        Self::with_params(reader, SdfReadParams::default())
    }
    pub fn with_params(reader: R, params: SdfReadParams) -> Self {
        Self {
            inner: cosmolkit_io::SdfGraphReader::with_params(reader, data_params(&params)),
            params,
        }
    }
    pub fn next_record(&mut self) -> Result<Option<SdfRecord>, SdfError> {
        let index = self.inner.records_consumed();
        self.inner
            .next_record()?
            .map(|parsed| SdfRecord::from_parsed(parsed, index))
            .transpose()
    }
    pub fn is_end(&self) -> bool {
        self.inner.is_end()
    }
    pub fn records_consumed(&self) -> usize {
        self.inner.records_consumed()
    }
    pub fn bytes_consumed(&self) -> u64 {
        self.inner.bytes_consumed()
    }
    pub fn lines_consumed(&self) -> usize {
        self.inner.lines_consumed()
    }
}
impl<R: BufRead> Iterator for SdfRecordStream<R> {
    type Item = Result<SdfRecord, SdfError>;
    fn next(&mut self) -> Option<Self::Item> {
        self.next_record().transpose()
    }
}
impl SdfRecordStream<BufReader<File>> {
    pub fn open(path: &str) -> Result<Self, crate::MolecularIoError> {
        Self::open_with_params(path, &SdfReadParams::default())
    }
    pub fn open_with_params(
        path: &str,
        params: &SdfReadParams,
    ) -> Result<Self, crate::MolecularIoError> {
        let file = crate::molecular_io::open_file(path)?;
        Ok(Self::with_params(BufReader::new(file), *params))
    }
}

/// An index contains only detached offsets; accessing a graph reads that record
/// and performs the existing checked construction lazily.
#[derive(Debug, Clone)]
pub struct SdfDataset {
    inner: std::sync::Arc<cosmolkit_io::SdfGraphDataset>,
    params: SdfReadParams,
}
impl SdfDataset {
    pub fn open(path: &str) -> Result<Self, crate::MolecularIoError> {
        Self::open_with_params(path, &SdfReadParams::default())
    }
    pub fn open_with_params(
        path: &str,
        params: &SdfReadParams,
    ) -> Result<Self, crate::MolecularIoError> {
        let path = crate::molecular_io::expanded_path(path)?;
        // Preserve real IO errors before the detached index is constructed.
        File::open(&path).map_err(|source| crate::MolecularIoError::Io {
            path: path.clone(),
            source,
        })?;
        let inner = cosmolkit_io::SdfGraphDataset::open_with_params(path, data_params(params))
            .map_err(|error| crate::MolecularIoError::Sdf(SdfError::from(error)))?;
        Ok(Self {
            inner: std::sync::Arc::new(inner),
            params: *params,
        })
    }
    #[cfg(feature = "cap-batch")]
    pub(crate) fn detached_dataset(&self) -> &cosmolkit_io::SdfGraphDataset {
        &self.inner
    }
    pub fn len(&self) -> usize {
        self.inner.len()
    }
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    pub fn path(&self) -> &Path {
        self.inner.path()
    }
    pub fn metadata(&self, index: usize) -> Option<&SdfRecordMetadata> {
        self.inner.metadata(index)
    }
    pub fn record(&self, index: usize) -> Result<SdfRecord, SdfError> {
        self.record_with_params(index, &self.params)
    }
    pub fn record_with_params(
        &self,
        index: usize,
        params: &SdfReadParams,
    ) -> Result<SdfRecord, SdfError> {
        let parsed = self.inner.record_with_params(index, data_params(params))?;
        SdfRecord::from_parsed(parsed, index)
    }
    pub fn record_text(&self, index: usize) -> Result<String, SdfError> {
        Ok(self.inner.record_text(index)?)
    }
    pub fn iter(&self) -> SdfDatasetIterator {
        SdfDatasetIterator {
            dataset: self.clone(),
            position: 0,
        }
    }
}
#[derive(Debug, Clone)]
pub struct SdfDatasetIterator {
    dataset: SdfDataset,
    position: usize,
}
impl Iterator for SdfDatasetIterator {
    type Item = Result<SdfRecord, SdfError>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.position >= self.dataset.len() {
            return None;
        }
        let index = self.position;
        self.position += 1;
        Some(self.dataset.record(index))
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        let n = self.dataset.len().saturating_sub(self.position);
        (n, Some(n))
    }
}
impl ExactSizeIterator for SdfDatasetIterator {}
impl std::iter::FusedIterator for SdfDatasetIterator {}

/// Reusable forward-file source. Opening a source resolves its path without
/// touching the file; batches() opens a fresh stream for each iteration.
#[derive(Debug, Clone)]
pub struct SdfReader {
    path: std::path::PathBuf,
    params: SdfReadParams,
}
impl SdfReader {
    pub fn open(path: &str) -> Result<Self, crate::MolecularIoError> {
        Self::open_with_params(path, &SdfReadParams::default())
    }
    pub fn open_with_params(
        path: &str,
        params: &SdfReadParams,
    ) -> Result<Self, crate::MolecularIoError> {
        Ok(Self {
            path: crate::molecular_io::expanded_path(path)?,
            params: *params,
        })
    }
    pub fn path(&self) -> &Path {
        &self.path
    }
    pub fn params(&self) -> &SdfReadParams {
        &self.params
    }
}
