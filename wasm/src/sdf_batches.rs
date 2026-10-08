//! Full SDF batch transports, reusing the canonical batch and supplier owners.
use crate::{BatchRecord, MoleculeBatch, SdfDataset, SdfReader, SdfRecordStream};
use cosmolkit as ck;
use std::{fs::File, io::BufReader};
pub struct SdfBatchIterator {
    inner: ck::SdfBatchIterator,
}
pub struct SdfReaderBatchIterator {
    inner: ck::SdfReaderBatchIterator<BufReader<File>>,
}
impl MoleculeBatch {
    pub fn from_sdf_records(text: &str) -> Result<Self, ck::BatchValidationError> {
        // COSMolKit❗✔️: ck::MoleculeBatch::from_sdf_records(text)
        ck::MoleculeBatch::from_sdf_records(text).map(|inner| Self { inner })
    }
    pub fn from_sdf_records_with_params(
        text: &str,
        read: &ck::SdfReadParams,
        mode: ck::BatchErrorMode,
        n_jobs: Option<usize>,
    ) -> Result<Self, ck::BatchValidationError> {
        // COSMolKit❗✔️: ck::MoleculeBatch::from_sdf_records_with_params(text,read,mode,n_jobs)
        ck::MoleculeBatch::from_sdf_records_with_params(text, read, mode, n_jobs)
            .map(|inner| Self { inner })
    }
    pub fn read_sdf(path: &str) -> Result<Self, ck::BatchValidationError> {
        // COSMolKit❗✔️: ck::MoleculeBatch::read_sdf(path)
        ck::MoleculeBatch::read_sdf(path).map(|inner| Self { inner })
    }
    pub fn read_sdf_with_params(
        path: &str,
        read: &ck::SdfReadParams,
        mode: ck::BatchErrorMode,
        n_jobs: Option<usize>,
        progress_bar: bool,
    ) -> Result<Self, ck::BatchValidationError> {
        // COSMolKit❗✔️: ck::MoleculeBatch::read_sdf_with_params(path,read,mode,n_jobs,progress_bar)
        ck::MoleculeBatch::read_sdf_with_params(path, read, mode, n_jobs, progress_bar)
            .map(|inner| Self { inner })
    }
    pub fn from_dataset_indices(
        dataset: &SdfDataset,
        indices: &[usize],
        mode: ck::BatchErrorMode,
    ) -> Result<Self, ck::BatchValidationError> {
        // COSMolKit❗✔️: ck::MoleculeBatch::from_dataset_indices(&dataset.inner,indices,mode)
        ck::MoleculeBatch::from_dataset_indices(&dataset.inner, indices, mode)
            .map(|inner| Self { inner })
    }
    pub fn get(&self, index: usize) -> Option<BatchRecord> {
        self.inner
            .get(index)
            .cloned()
            .map(|inner| BatchRecord { inner })
    }
    pub fn records(&self) -> Vec<BatchRecord> {
        self.inner
            .records()
            .iter()
            .cloned()
            .map(|inner| BatchRecord { inner })
            .collect()
    }
    pub fn to_sdf(&self, path: &str) -> Result<ck::BatchExportReport, ck::BatchValidationError> {
        // COSMolKit❗✔️: self.inner.to_sdf(path)
        self.inner.to_sdf(path)
    }
    pub fn to_sdf_files(
        &self,
        path: &str,
    ) -> Result<ck::BatchExportReport, ck::BatchValidationError> {
        // COSMolKit❗✔️: self.inner.to_sdf_files(path)
        self.inner.to_sdf_files(path)
    }
    pub fn to_sdf_with_params(
        &self,
        path: &str,
        params: &ck::BatchExportParams,
        report_path: Option<&str>,
    ) -> Result<ck::BatchExportReport, ck::BatchValidationError> {
        // COSMolKit❗✔️: self.inner.to_sdf_with_params(path,params,report_path)
        self.inner.to_sdf_with_params(path, params, report_path)
    }
    pub fn to_sdf_files_with_params(
        &self,
        directory: &str,
        params: &ck::BatchExportParams,
        filenames: Option<&[Option<String>]>,
        report_path: Option<&str>,
    ) -> Result<ck::BatchExportReport, ck::BatchValidationError> {
        // COSMolKit❗✔️: self.inner.to_sdf_files_with_params(directory,params,filenames,report_path)
        self.inner
            .to_sdf_files_with_params(directory, params, filenames, report_path)
    }
}
impl SdfDataset {
    pub fn batches(
        &self,
        size: usize,
        mode: ck::BatchErrorMode,
        n_jobs: Option<usize>,
        progress_bar: bool,
    ) -> Result<SdfBatchIterator, ck::BatchValidationError> {
        // COSMolKit❗✔️: self.inner.batches(size,mode,n_jobs,progress_bar)
        self.inner
            .batches(size, mode, n_jobs, progress_bar)
            .map(|inner| SdfBatchIterator { inner })
    }
}
impl SdfReader {
    pub fn batches(
        &self,
        size: usize,
        mode: ck::BatchErrorMode,
        n_jobs: Option<usize>,
    ) -> Result<SdfReaderBatchIterator, ck::BatchValidationError> {
        // COSMolKit❗✔️: self.inner.batches(size,mode,n_jobs)
        self.inner
            .batches(size, mode, n_jobs)
            .map(|inner| SdfReaderBatchIterator { inner })
    }
}
impl SdfRecordStream {
    pub fn batches(
        self,
        size: usize,
        mode: ck::BatchErrorMode,
        n_jobs: Option<usize>,
    ) -> Result<SdfReaderBatchIterator, ck::BatchValidationError> {
        // COSMolKit❗✔️: self.inner.batches(size,mode,n_jobs)
        self.inner
            .batches(size, mode, n_jobs)
            .map(|inner| SdfReaderBatchIterator { inner })
    }
}
impl SdfBatchIterator {
    pub fn next_batch(&mut self) -> Result<Option<MoleculeBatch>, ck::BatchValidationError> {
        self.inner
            .next_batch()
            .map(|b| b.map(|inner| MoleculeBatch { inner }))
    }
}
impl SdfReaderBatchIterator {
    pub fn next_batch(&mut self) -> Result<Option<MoleculeBatch>, ck::BatchValidationError> {
        self.inner
            .next_batch()
            .map(|b| b.map(|inner| MoleculeBatch { inner }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const RECORD: &str = "batch\n  COSMolKit         2D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n";
    fn masks(a: &MoleculeBatch, b: &ck::MoleculeBatch) {
        assert_eq!(a.valid_mask(), b.valid_mask());
        assert_eq!(
            a.errors()
                .iter()
                .map(|e| (e.index, e.operation, e.message.clone()))
                .collect::<Vec<_>>(),
            b.errors()
                .iter()
                .map(|e| (e.index, e.operation, e.message.clone()))
                .collect::<Vec<_>>()
        );
    }
    #[test]
    fn sdf_text_records_index_selection_and_all_iterators_preserve_order() {
        let text = RECORD.repeat(3);
        let batch = MoleculeBatch::from_sdf_records(&text).unwrap();
        assert_eq!(batch.len(), 3);
        assert_eq!(batch.records().len(), 3);
        assert!(batch.get(3).is_none());
        assert!(batch.get(0).unwrap().molecule_value().is_some());
        let root =
            std::env::temp_dir().join(format!("cosmolkit-sdf-batch-read-{}", std::process::id()));
        std::fs::create_dir(&root).unwrap();
        let path = root.join("input.sdf");
        std::fs::write(&path, &text).unwrap();
        let path = path.to_str().unwrap();
        let read = ck::SdfReadParams::default();
        for mode in [ck::BatchErrorMode::Strict, ck::BatchErrorMode::KeepErrors] {
            for jobs in [None, Some(1), Some(2)] {
                let b =
                    MoleculeBatch::from_sdf_records_with_params(&text, &read, mode, jobs).unwrap();
                let owner =
                    ck::MoleculeBatch::from_sdf_records_with_params(&text, &read, mode, jobs)
                        .unwrap();
                masks(&b, &owner);
                for progress in [false, true] {
                    let b = MoleculeBatch::read_sdf_with_params(path, &read, mode, jobs, progress)
                        .unwrap();
                    masks(
                        &b,
                        &ck::MoleculeBatch::read_sdf_with_params(path, &read, mode, jobs, progress)
                            .unwrap(),
                    );
                    let ds = SdfDataset::open_with_params(path, &read).unwrap();
                    let mut it = ds.batches(2, mode, jobs, progress).unwrap();
                    assert_eq!(it.next_batch().unwrap().unwrap().len(), 2);
                    assert_eq!(it.next_batch().unwrap().unwrap().len(), 1);
                    assert!(it.next_batch().unwrap().is_none());
                    assert!(it.next_batch().unwrap().is_none());
                }
                let reader = SdfReader::open_with_params(path, &read).unwrap();
                let mut iter = reader.batches(2, mode, jobs).unwrap();
                assert_eq!(iter.next_batch().unwrap().unwrap().len(), 2);
                assert_eq!(iter.next_batch().unwrap().unwrap().len(), 1);
                assert!(iter.next_batch().unwrap().is_none());
                let mut stream = SdfRecordStream::open(path)
                    .unwrap()
                    .batches(2, mode, jobs)
                    .unwrap();
                assert_eq!(stream.next_batch().unwrap().unwrap().len(), 2);
                assert_eq!(stream.next_batch().unwrap().unwrap().len(), 1);
                assert!(stream.next_batch().unwrap().is_none());
            }
        }
        assert_eq!(MoleculeBatch::read_sdf(path).unwrap().len(), 3);
        let ds = SdfDataset::open(path).unwrap();
        let selected = MoleculeBatch::from_dataset_indices(
            &ds,
            &[2, 0, 2, 99],
            ck::BatchErrorMode::KeepErrors,
        )
        .unwrap();
        assert_eq!(selected.valid_mask(), [true, true, true, false]);
        assert_eq!(selected.errors()[0].index, 99);
        assert!(
            MoleculeBatch::from_dataset_indices(&ds, &[99], ck::BatchErrorMode::Strict).is_err()
        );
        assert!(
            ds.batches(0, ck::BatchErrorMode::Strict, None, false)
                .is_err()
        );
        assert!(
            SdfReader::open(path)
                .unwrap()
                .batches(1, ck::BatchErrorMode::Strict, Some(0))
                .is_err()
        );
        assert!(
            MoleculeBatch::from_sdf_records_with_params(
                &text,
                &read,
                ck::BatchErrorMode::Strict,
                Some(0)
            )
            .is_err()
        );
        std::fs::remove_file(path).unwrap();
        std::fs::remove_dir(root).unwrap();
    }
    #[test]
    fn exports_reports_formats_and_strict_validation_match_owner() {
        let batch = MoleculeBatch::from_sdf_records(&RECORD.repeat(2)).unwrap();
        let root =
            std::env::temp_dir().join(format!("cosmolkit-sdf-batch-export-{}", std::process::id()));
        std::fs::create_dir(&root).unwrap();
        let path = root.join("default.sdf");
        let r = batch.to_sdf(path.to_str().unwrap()).unwrap();
        assert_eq!(
            (r.total(), r.success(), r.failed(), r.errors().len()),
            (2, 2, 0, 0)
        );
        let dir = root.join("default-files");
        assert_eq!(
            batch.to_sdf_files(dir.to_str().unwrap()).unwrap().success(),
            2
        );
        for (format, index) in [(ck::SdfFormat::V2000, 0), (ck::SdfFormat::V3000, 1)] {
            let p = ck::BatchExportParams {
                format,
                errors: ck::BatchErrorMode::Strict,
                n_jobs: Some(2),
                progress_bar: Some(false),
            };
            let output = root.join(format!("out-{index}.sdf"));
            let report = root.join(format!("report-{index}.json"));
            let r = batch
                .to_sdf_with_params(output.to_str().unwrap(), &p, Some(report.to_str().unwrap()))
                .unwrap();
            assert_eq!(r.success(), 2);
            assert!(
                std::fs::read_to_string(output)
                    .unwrap()
                    .contains(if index == 0 { "V2000" } else { "V3000" })
            );
            assert!(std::fs::read_to_string(report).unwrap().contains("written"));
            let directory = root.join(format!("files-{index}"));
            let names = [Some("a.sdf".into()), Some("b.sdf".into())];
            let r = batch
                .to_sdf_files_with_params(directory.to_str().unwrap(), &p, Some(&names), None)
                .unwrap();
            assert_eq!(r.success(), 2);
            assert!(directory.join("a.sdf").exists());
            assert!(directory.join("b.sdf").exists());
        }
        let bad = crate::BatchRecord::error(&ck::BatchError::new(7, "test", "invalid"));
        let mixed = MoleculeBatch::from_records(
            vec![batch.get(0).unwrap(), bad],
            ck::BatchErrorMode::KeepErrors,
        )
        .unwrap();
        let fail = root.join("strict-must-not-open.sdf");
        assert!(mixed.to_sdf(fail.to_str().unwrap()).is_err());
        assert!(!fail.exists());
        let keep = ck::BatchExportParams {
            errors: ck::BatchErrorMode::KeepErrors,
            ..Default::default()
        };
        let r = mixed
            .to_sdf_with_params(fail.to_str().unwrap(), &keep, None)
            .unwrap();
        assert_eq!(
            (r.total(), r.success(), r.failed(), r.skipped),
            (2, 1, 0, 1)
        );
        assert!(r.errors().is_empty());
        for entry in std::fs::read_dir(&root).unwrap() {
            let p = entry.unwrap().path();
            if p.is_dir() {
                for f in std::fs::read_dir(&p).unwrap() {
                    std::fs::remove_file(f.unwrap().path()).unwrap();
                }
                std::fs::remove_dir(p).unwrap();
            } else {
                std::fs::remove_file(p).unwrap();
            }
        }
        std::fs::remove_dir(root).unwrap();
    }
}
