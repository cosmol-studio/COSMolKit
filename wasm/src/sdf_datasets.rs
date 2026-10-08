//! Fixed file-backed supplier carriers; all indexing and framing remain canonical.
use crate::SdfRecord;
use cosmolkit as ck;
use std::fs::File;
use std::io::BufReader;
#[derive(Clone, Debug)]
pub struct SdfDataset {
    pub(crate) inner: ck::SdfDataset,
}
#[derive(Clone, Debug)]
pub struct SdfDatasetIterator {
    pub(crate) inner: ck::SdfDatasetIterator,
}
pub struct SdfRecordStream {
    pub(crate) inner: ck::SdfRecordStream<BufReader<File>>,
}
#[derive(Clone, Debug)]
pub struct SdfReader {
    pub(crate) inner: ck::SdfReader,
}
impl SdfDataset {
    pub fn open(path: &str) -> Result<Self, ck::MolecularIoError> {
        // COSMolKit❗✔️: ck::SdfDataset::open(path)
        ck::SdfDataset::open(path).map(|inner| Self { inner })
    }
    pub fn open_with_params(
        path: &str,
        params: &ck::SdfReadParams,
    ) -> Result<Self, ck::MolecularIoError> {
        // COSMolKit❗✔️: ck::SdfDataset::open_with_params(path,params)
        ck::SdfDataset::open_with_params(path, params).map(|inner| Self { inner })
    }
}
impl SdfRecordStream {
    pub fn open(path: &str) -> Result<Self, ck::MolecularIoError> {
        // COSMolKit❗✔️: ck::SdfRecordStream::open(path)
        ck::SdfRecordStream::open(path).map(|inner| Self { inner })
    }
    pub fn open_with_params(
        path: &str,
        params: &ck::SdfReadParams,
    ) -> Result<Self, ck::MolecularIoError> {
        // COSMolKit❗✔️: ck::SdfRecordStream::open_with_params(path,params)
        ck::SdfRecordStream::open_with_params(path, params).map(|inner| Self { inner })
    }
}
impl SdfReader {
    pub fn open(path: &str) -> Result<Self, ck::MolecularIoError> {
        // COSMolKit❗✔️: ck::SdfReader::open(path)
        ck::SdfReader::open(path).map(|inner| Self { inner })
    }
    pub fn open_with_params(
        path: &str,
        params: &ck::SdfReadParams,
    ) -> Result<Self, ck::MolecularIoError> {
        // COSMolKit❗✔️: ck::SdfReader::open_with_params(path,params)
        ck::SdfReader::open_with_params(path, params).map(|inner| Self { inner })
    }
}
impl SdfDataset {
    pub fn len(&self) -> usize {
        self.inner.len()
    }
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    pub fn path(&self) -> &std::path::Path {
        self.inner.path()
    }
    pub fn metadata(&self, index: usize) -> Option<ck::SdfRecordMetadata> {
        self.inner.metadata(index).cloned()
    }
    pub fn record(&self, index: usize) -> Result<SdfRecord, ck::SdfError> {
        self.inner.record(index).map(|inner| SdfRecord { inner })
    }
    pub fn record_with_params(
        &self,
        index: usize,
        params: &ck::SdfReadParams,
    ) -> Result<SdfRecord, ck::SdfError> {
        self.inner
            .record_with_params(index, params)
            .map(|inner| SdfRecord { inner })
    }
    pub fn record_text(&self, index: usize) -> Result<String, ck::SdfError> {
        self.inner.record_text(index)
    }
    pub fn iter(&self) -> SdfDatasetIterator {
        SdfDatasetIterator {
            inner: self.inner.iter(),
        }
    }
}
impl Iterator for SdfDatasetIterator {
    type Item = Result<SdfRecord, ck::SdfError>;
    fn next(&mut self) -> Option<Self::Item> {
        self.inner
            .next()
            .map(|r| r.map(|inner| SdfRecord { inner }))
    }
}
impl SdfRecordStream {
    pub fn next_record(&mut self) -> Result<Option<SdfRecord>, ck::SdfError> {
        self.inner
            .next_record()
            .map(|r| r.map(|inner| SdfRecord { inner }))
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
impl SdfReader {
    pub fn path(&self) -> &std::path::Path {
        self.inner.path()
    }
    pub fn params(&self) -> ck::SdfReadParams {
        *self.inner.params()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const FIRST: &str = "first\n  COSMolKit         2D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n>  <NOTE>\none\n\n$$$$\n";
    #[test]
    fn indexed_records_metadata_iteration_and_stream_counters_match_owner() {
        let second = FIRST.replacen("first", "second", 1).replace("\n", "\r\n");
        let text = format!("{FIRST}{second}");
        let path =
            std::env::temp_dir().join(format!("cosmolkit-suppliers-{}.sdf", std::process::id()));
        std::fs::write(&path, &text).unwrap();
        let path = path.to_str().unwrap();
        let dataset = SdfDataset::open(path).unwrap();
        let owner = ck::SdfDataset::open(path).unwrap();
        assert_eq!(dataset.len(), 2);
        assert!(!dataset.is_empty());
        assert_eq!(dataset.path(), owner.path());
        assert_eq!(dataset.metadata(2), None);
        assert_eq!(dataset.metadata(usize::MAX), None);
        for index in [0, 1] {
            let metadata = dataset.metadata(index).unwrap();
            assert_eq!(&metadata, owner.metadata(index).unwrap());
            assert_eq!(metadata.index(), index);
            assert_eq!(
                metadata.title(),
                Some(if index == 0 { "first" } else { "second" })
            );
            let begin = if index == 0 { 0 } else { FIRST.len() };
            let len = if index == 0 {
                FIRST.len()
            } else {
                second.len()
            };
            assert_eq!(metadata.byte_offset(), begin as u64);
            assert_eq!(metadata.byte_len(), len as u64);
            assert_eq!(metadata.byte_range(), (begin as u64, (begin + len) as u64));
            let line_offset = if index == 0 { 0 } else { FIRST.lines().count() };
            assert_eq!(
                metadata.line_range(),
                (line_offset, line_offset + FIRST.lines().count())
            );
            assert_eq!(
                dataset.record_text(index).unwrap(),
                owner.record_text(index).unwrap()
            );
            assert_eq!(dataset.record(index).unwrap().index(), index);
            let params = ck::SdfReadParams {
                sanitize: false,
                remove_hydrogens: false,
                ..Default::default()
            };
            assert_eq!(
                dataset
                    .record_with_params(index, &params)
                    .unwrap()
                    .molecule()
                    .unwrap()
                    .num_atoms(),
                1
            );
        }
        assert_eq!(
            dataset.record(2).unwrap_err().to_string(),
            owner.record(2).unwrap_err().to_string()
        );
        assert_eq!(
            dataset.record_text(2).unwrap_err().to_string(),
            owner.record_text(2).unwrap_err().to_string()
        );
        let mut iter = dataset.iter();
        drop(dataset);
        assert_eq!(iter.next().unwrap().unwrap().title(), Some("first".into()));
        assert_eq!(iter.next().unwrap().unwrap().title(), Some("second".into()));
        assert!(iter.next().is_none());
        assert!(iter.next().is_none());
        let params = ck::SdfReadParams::default();
        let mut stream = SdfRecordStream::open_with_params(path, &params).unwrap();
        let mut canonical = ck::SdfRecordStream::open_with_params(path, &params).unwrap();
        assert_eq!(stream.records_consumed(), 0);
        assert_eq!(stream.bytes_consumed(), 0);
        assert_eq!(stream.lines_consumed(), 0);
        assert!(!stream.is_end());
        for index in 0..2 {
            assert_eq!(stream.next_record().unwrap().unwrap().index(), index);
            assert_eq!(canonical.next_record().unwrap().unwrap().index(), index);
            assert_eq!(stream.bytes_consumed(), canonical.bytes_consumed());
            assert_eq!(stream.lines_consumed(), canonical.lines_consumed());
            assert_eq!(stream.records_consumed(), canonical.records_consumed());
            assert_eq!(stream.is_end(), canonical.is_end());
        }
        assert!(stream.next_record().unwrap().is_none());
        assert!(stream.is_end());
        assert_eq!(stream.records_consumed(), 2);
        assert_eq!(stream.bytes_consumed(), text.len() as u64);
        assert_eq!(stream.lines_consumed(), text.lines().count());
        assert!(stream.next_record().unwrap().is_none());
        assert!(
            SdfRecordStream::open(path)
                .unwrap()
                .next_record()
                .unwrap()
                .is_some()
        );
        let reader = SdfReader::open_with_params(path, &params).unwrap();
        assert_eq!(reader.path(), std::path::Path::new(path));
        assert_eq!(reader.params(), params);
        assert_eq!(SdfReader::open(path).unwrap().params(), params);
        std::fs::remove_file(path).unwrap();
    }
    #[test]
    fn lazy_reader_missing_files_and_malformed_records_preserve_errors() {
        let path = std::env::temp_dir().join(format!(
            "cosmolkit-supplier-errors-{}.sdf",
            std::process::id()
        ));
        let path = path.to_str().unwrap();
        assert!(SdfReader::open(path).is_ok());
        assert!(
            matches!(SdfDataset::open(path),Err(ck::MolecularIoError::Io{source,..})if source.kind()==std::io::ErrorKind::NotFound)
        );
        assert!(
            matches!(SdfRecordStream::open(path),Err(ck::MolecularIoError::Io{source,..})if source.kind()==std::io::ErrorKind::NotFound)
        );
        let text = format!("{FIRST}bad\nproducer\n\ninvalid counts\njunk\n$$$$\n");
        std::fs::write(path, text).unwrap();
        let ds = SdfDataset::open(path).unwrap();
        assert_eq!(ds.len(), 2);
        assert!(ds.record(0).is_ok());
        assert!(ds.record(1).is_err());
        let mut it = ds.iter();
        assert!(it.next().unwrap().is_ok());
        assert!(it.next().unwrap().is_err());
        assert!(it.next().is_none());
        let mut stream = SdfRecordStream::open(path).unwrap();
        let mut owner = ck::SdfRecordStream::open(path).unwrap();
        assert!(stream.next_record().unwrap().is_some());
        assert!(owner.next_record().unwrap().is_some());
        assert_eq!(
            stream.next_record().unwrap_err().to_string(),
            owner.next_record().unwrap_err().to_string()
        );
        assert_eq!(stream.records_consumed(), owner.records_consumed());
        assert_eq!(stream.bytes_consumed(), owner.bytes_consumed());
        assert_eq!(stream.lines_consumed(), owner.lines_consumed());
        assert_eq!(stream.is_end(), owner.is_end());
        assert!(stream.next_record().unwrap().is_none());
        std::fs::remove_file(path).unwrap();
    }
}
