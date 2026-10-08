//! Language values wrapping the sole public live batch and molecule owners.
use crate::Molecule;
use cosmolkit as ck;

#[derive(Clone)]
pub struct BatchRecord {
    pub(crate) inner: ck::BatchRecord,
}

impl BatchRecord {
    pub fn molecule(value: &Molecule) -> Self {
        Self {
            inner: ck::BatchRecord::Molecule(value.inner.borrow().clone()),
        }
    }
    pub fn error(value: &ck::BatchError) -> Self {
        Self {
            inner: ck::BatchRecord::Error(value.clone()),
        }
    }
    pub fn molecule_value(&self) -> Option<Molecule> {
        match &self.inner {
            ck::BatchRecord::Molecule(value) => Some(Molecule {
                inner: value.clone().into(),
            }),
            ck::BatchRecord::Error(_) => None,
        }
    }
    pub fn error_value(&self) -> Option<ck::BatchError> {
        match &self.inner {
            ck::BatchRecord::Molecule(_) => None,
            ck::BatchRecord::Error(error) => Some(error.clone()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn batch_foundation_preserves_indices_settings_and_record_value_independence() {
        let input = ["CCO", "[", "O", "C("].map(str::to_owned);
        let params = ck::BatchParams {
            errors: ck::BatchErrorMode::KeepErrors,
            ..Default::default()
        };
        let batch = MoleculeBatch::from_smiles_list_with_params(
            &input,
            &ck::SmilesParseParams::default(),
            &params,
        )
        .unwrap();
        assert_eq!(batch.len(), 4);
        assert!(!batch.is_empty());
        assert_eq!(batch.valid_mask(), [true, false, true, false]);
        assert_eq!(batch.invalid_mask(), [false, true, false, true]);
        assert_eq!(batch.valid_count(), 2);
        assert_eq!(batch.invalid_count(), 2);
        assert_eq!(
            batch.errors().iter().map(|e| e.index).collect::<Vec<_>>(),
            [1, 3]
        );
        assert_eq!(
            batch
                .to_list()
                .iter()
                .map(Option::is_some)
                .collect::<Vec<_>>(),
            batch.valid_mask()
        );
        assert_eq!(batch.parallel_jobs(), None);
        assert_eq!(batch.progress_bar(), None);
        let configured = batch
            .with_parallel_jobs(Some(2))
            .unwrap()
            .with_progress_bar(Some(false));
        assert_eq!(configured.parallel_jobs(), Some(2));
        assert_eq!(configured.progress_bar(), Some(false));
        assert_eq!(batch.parallel_jobs(), None);
        assert!(batch.with_parallel_jobs(Some(0)).is_err());
        assert_eq!(batch.with_valid_records().len(), 2);
        let strict = MoleculeBatch::from_smiles_list(&input).err().unwrap();
        assert_eq!(strict.errors, 2);
        assert_eq!(
            strict
                .record_errors
                .iter()
                .map(|e| e.index)
                .collect::<Vec<_>>(),
            [1, 3]
        );
        let source = Molecule::from_smiles_with_sanitize("C1=CC=CC=C1", false).unwrap();
        let record = BatchRecord::molecule(&source);
        let error = BatchRecord::error(&batch.errors()[0]);
        assert!(record.error_value().is_none());
        assert!(error.molecule_value().is_none());
        let copied = MoleculeBatch::from_records(
            vec![record.clone(), error.clone()],
            ck::BatchErrorMode::KeepErrors,
        )
        .unwrap();
        source.assign_aromaticity_().unwrap();
        assert_eq!(source.to_smiles().unwrap(), "c1ccccc1".into());
        assert_eq!(
            record.molecule_value().unwrap().to_smiles().unwrap(),
            "C1=CC=CC=C1".into()
        );
        assert_eq!(
            copied.to_list()[0].as_ref().unwrap().to_smiles().unwrap(),
            "C1=CC=CC=C1".into()
        );
        assert_eq!(error.error_value().unwrap().index, 1);
        assert!(MoleculeBatch::from_records(vec![error], ck::BatchErrorMode::Strict).is_err());
        assert!(
            MoleculeBatch::from_records(Vec::new(), ck::BatchErrorMode::Strict)
                .unwrap()
                .is_empty()
        );
    }
}

#[derive(Clone)]
pub struct MoleculeBatch {
    pub(crate) inner: ck::MoleculeBatch,
}

impl MoleculeBatch {
    pub fn from_records(
        records: Vec<BatchRecord>,
        mode: ck::BatchErrorMode,
    ) -> Result<Self, ck::BatchValidationError> {
        ck::MoleculeBatch::from_records(
            records.into_iter().map(|value| value.inner).collect(),
            mode,
        )
        .map(|inner| Self { inner })
    }
    #[cfg(feature = "smiles")]
    pub fn from_smiles_list(smiles: &[String]) -> Result<Self, ck::BatchValidationError> {
        ck::MoleculeBatch::from_smiles_list(smiles).map(|inner| Self { inner })
    }
    #[cfg(feature = "smiles")]
    pub fn from_smiles_list_with_params(
        smiles: &[String],
        parse: &ck::SmilesParseParams,
        params: &ck::BatchParams,
    ) -> Result<Self, ck::BatchValidationError> {
        ck::MoleculeBatch::from_smiles_list_with_params(smiles, parse, params)
            .map(|inner| Self { inner })
    }
    pub fn len(&self) -> usize {
        self.inner.len()
    }
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    pub fn valid_mask(&self) -> Vec<bool> {
        self.inner.valid_mask()
    }
    pub fn invalid_mask(&self) -> Vec<bool> {
        self.inner.invalid_mask()
    }
    pub fn valid_count(&self) -> usize {
        self.inner.valid_count()
    }
    pub fn invalid_count(&self) -> usize {
        self.inner.invalid_count()
    }
    pub fn errors(&self) -> Vec<ck::BatchError> {
        self.inner.errors()
    }
    pub fn parallel_jobs(&self) -> Option<usize> {
        self.inner.parallel_jobs()
    }
    pub fn progress_bar(&self) -> Option<bool> {
        self.inner.progress_bar()
    }
    pub fn to_list(&self) -> Vec<Option<Molecule>> {
        self.inner
            .to_list()
            .into_iter()
            .map(|value| {
                value.map(|inner| Molecule {
                    inner: inner.into(),
                })
            })
            .collect()
    }
    pub fn with_valid_records(&self) -> Self {
        Self {
            inner: self.inner.with_valid_records(),
        }
    }
    pub fn with_parallel_jobs(
        &self,
        n_jobs: Option<usize>,
    ) -> Result<Self, ck::BatchValidationError> {
        self.inner
            .clone()
            .with_parallel_jobs(n_jobs)
            .map(|inner| Self { inner })
    }
    pub fn with_progress_bar(&self, progress_bar: Option<bool>) -> Self {
        Self {
            inner: self.inner.clone().with_progress_bar(progress_bar),
        }
    }
}

impl MoleculeBatch {
    #[cfg(feature = "sanitize")]
    pub fn sanitize(&self) -> Result<Self, ck::BatchValidationError> {
        self.inner.sanitize().map(|inner| Self { inner })
    }
    #[cfg(feature = "sanitize")]
    pub fn sanitize_with_params(
        &self,
        options: &ck::SanitizeParams,
        params: &ck::BatchParams,
    ) -> Result<Self, ck::BatchValidationError> {
        self.inner
            .sanitize_with_params(options, params)
            .map(|inner| Self { inner })
    }
    #[cfg(feature = "hydrogens")]
    pub fn with_hydrogens(&self) -> Result<Self, ck::BatchValidationError> {
        self.inner.with_hydrogens().map(|inner| Self { inner })
    }
    #[cfg(feature = "hydrogens")]
    pub fn with_hydrogens_with_params(
        &self,
        options: &ck::AddHsParams,
        params: &ck::BatchParams,
    ) -> Result<Self, ck::BatchValidationError> {
        self.inner
            .with_hydrogens_with_params(options, params)
            .map(|inner| Self { inner })
    }
    #[cfg(feature = "hydrogens")]
    pub fn without_hydrogens(&self) -> Result<Self, ck::BatchValidationError> {
        self.inner.without_hydrogens().map(|inner| Self { inner })
    }
    #[cfg(feature = "hydrogens")]
    pub fn without_hydrogens_with_params(
        &self,
        options: &ck::RemoveHsParams,
        params: &ck::BatchParams,
    ) -> Result<Self, ck::BatchValidationError> {
        self.inner
            .without_hydrogens_with_params(options, params)
            .map(|inner| Self { inner })
    }
    #[cfg(feature = "kekulize")]
    pub fn with_kekulized_bonds(&self) -> Result<Self, ck::BatchValidationError> {
        self.inner
            .with_kekulized_bonds()
            .map(|inner| Self { inner })
    }
    #[cfg(feature = "kekulize")]
    pub fn with_kekulized_bonds_with_params(
        &self,
        options: &ck::KekulizeParams,
        params: &ck::BatchParams,
    ) -> Result<Self, ck::BatchValidationError> {
        self.inner
            .with_kekulized_bonds_with_params(options, params)
            .map(|inner| Self { inner })
    }
    #[cfg(feature = "depict")]
    pub fn with_2d_coordinates(&self) -> Result<Self, ck::BatchValidationError> {
        self.inner.with_2d_coordinates().map(|inner| Self { inner })
    }
    #[cfg(feature = "depict")]
    pub fn with_2d_coordinates_with_params(
        &self,
        options: &ck::Coordinate2DParams,
        params: &ck::BatchParams,
    ) -> Result<Self, ck::BatchValidationError> {
        self.inner
            .with_2d_coordinates_with_params(options, params)
            .map(|inner| Self { inner })
    }
}

#[cfg(test)]
mod transform_tests {
    use super::*;
    fn snapshot(batch: &MoleculeBatch) -> Vec<Option<(usize, ck::PropertyText, Vec<f64>)>> {
        batch
            .to_list()
            .iter()
            .map(|v| {
                v.as_ref().map(|m| {
                    (
                        m.num_atoms() as usize,
                        m.to_smiles().unwrap(),
                        m.coordinates_2d(),
                    )
                })
            })
            .collect()
    }
    fn facade_snapshot(
        batch: ck::MoleculeBatch,
    ) -> Vec<Option<(usize, ck::PropertyText, Vec<f64>)>> {
        batch
            .to_list()
            .iter()
            .map(|v| {
                v.as_ref().map(|m| {
                    (
                        m.num_atoms(),
                        m.to_smiles().unwrap(),
                        m.coordinates_2d()
                            .map_or(Vec::new(), |v| v.iter().flatten().copied().collect()),
                    )
                })
            })
            .collect()
    }
    #[test]
    fn batch_transforms_match_all_ten_facade_calls_and_preserve_input() {
        let input = ["CCO".to_owned(), "c1ccccc1".to_owned()];
        let batch = MoleculeBatch::from_smiles_list(&input).unwrap();
        let source = ck::MoleculeBatch::from_smiles_list(&input).unwrap();
        let before = snapshot(&batch);
        let execution = ck::BatchParams::default();
        macro_rules! compare {
            ($default:ident,$with:ident,$params:expr) => {{
                assert_eq!(
                    snapshot(&batch.$default().unwrap()),
                    facade_snapshot(source.$default().unwrap())
                );
                assert_eq!(
                    snapshot(&batch.$with(&$params, &execution).unwrap()),
                    facade_snapshot(source.$with(&$params, &execution).unwrap())
                );
                assert_eq!(snapshot(&batch), before);
            }};
        }
        compare!(
            sanitize,
            sanitize_with_params,
            ck::SanitizeParams::default()
        );
        compare!(
            with_hydrogens,
            with_hydrogens_with_params,
            ck::AddHsParams::default()
        );
        compare!(
            without_hydrogens,
            without_hydrogens_with_params,
            ck::RemoveHsParams::default()
        );
        compare!(
            with_kekulized_bonds,
            with_kekulized_bonds_with_params,
            ck::KekulizeParams::default()
        );
        compare!(
            with_2d_coordinates,
            with_2d_coordinates_with_params,
            ck::Coordinate2DParams::default()
        );
        let invalid = ck::AddHsParams {
            only_on_atoms: Some(vec![ck::AtomId::new(99)]),
            ..Default::default()
        };
        let failure = batch
            .with_hydrogens_with_params(&invalid, &execution)
            .err()
            .unwrap();
        assert_eq!(failure.errors, 2);
        assert_eq!(failure.record_errors[0].index, 0);
        let kept = batch
            .with_hydrogens_with_params(
                &invalid,
                &ck::BatchParams {
                    errors: ck::BatchErrorMode::KeepErrors,
                    ..Default::default()
                },
            )
            .unwrap();
        assert_eq!(kept.valid_mask(), [false, false]);
        assert_eq!(snapshot(&batch), before);
        let partial = MoleculeBatch::from_smiles_list_with_params(
            &["C".into(), "[".into()],
            &ck::SmilesParseParams::default(),
            &ck::BatchParams {
                errors: ck::BatchErrorMode::KeepErrors,
                ..Default::default()
            },
        )
        .unwrap();
        let kept = partial
            .with_hydrogens_with_params(
                &ck::AddHsParams::default(),
                &ck::BatchParams {
                    errors: ck::BatchErrorMode::KeepErrors,
                    ..Default::default()
                },
            )
            .unwrap();
        assert_eq!(kept.valid_mask(), [true, false]);
        assert_eq!(kept.errors()[0].operation, "batch.from_smiles_list");
        let empty = MoleculeBatch::from_smiles_list(&[]).unwrap();
        assert!(empty.sanitize().unwrap().is_empty());
    }
}

impl MoleculeBatch {
    #[cfg(feature = "smiles")]
    pub fn to_smiles_list(
        &self,
    ) -> Result<Vec<Option<ck::PropertyText>>, ck::BatchValidationError> {
        self.inner.to_smiles_list()
    }
    #[cfg(feature = "smiles")]
    pub fn to_smiles_list_with_params(
        &self,
        options: &ck::SmilesWriteParams,
        params: &ck::BatchQueryParams,
    ) -> Result<Vec<Option<ck::PropertyText>>, ck::BatchValidationError> {
        self.inner.to_smiles_list_with_params(options, params)
    }
    #[cfg(feature = "conformer")]
    pub fn dg_bounds_matrix_list(
        &self,
    ) -> Result<Vec<Option<Vec<Vec<f64>>>>, ck::BatchValidationError> {
        self.inner.dg_bounds_matrix_list()
    }
    #[cfg(feature = "conformer")]
    pub fn dg_bounds_matrix_list_with_params(
        &self,
        params: &ck::BatchQueryParams,
    ) -> Result<Vec<Option<Vec<Vec<f64>>>>, ck::BatchValidationError> {
        self.inner.dg_bounds_matrix_list_with_params(params)
    }
    #[cfg(feature = "depict")]
    pub fn to_svg_list(
        &self,
        width: u32,
        height: u32,
    ) -> Result<Vec<Option<String>>, ck::BatchValidationError> {
        self.inner.to_svg_list(width, height)
    }
    #[cfg(feature = "depict")]
    pub fn to_svg_list_with_params(
        &self,
        width: u32,
        height: u32,
        params: &ck::BatchQueryParams,
    ) -> Result<Vec<Option<String>>, ck::BatchValidationError> {
        self.inner.to_svg_list_with_params(width, height, params)
    }
}

#[cfg(test)]
mod query_tests {
    use super::*;
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };
    #[test]
    fn batch_queries_match_six_facade_calls_and_notify_every_row() {
        let input = ["CCO".into(), "[".into(), "O".into()];
        let parse = ck::SmilesParseParams::default();
        let keep = ck::BatchParams {
            errors: ck::BatchErrorMode::KeepErrors,
            ..Default::default()
        };
        let batch = MoleculeBatch::from_smiles_list_with_params(&input, &parse, &keep).unwrap();
        let source =
            ck::MoleculeBatch::from_smiles_list_with_params(&input, &parse, &keep).unwrap();
        assert_eq!(
            batch.to_smiles_list().unwrap(),
            source.to_smiles_list().unwrap()
        );
        assert_eq!(
            batch.dg_bounds_matrix_list().unwrap(),
            source.dg_bounds_matrix_list().unwrap()
        );
        assert_eq!(
            batch.to_svg_list(120, 100).unwrap(),
            source.to_svg_list(120, 100).unwrap()
        );
        let ticks = Arc::new(AtomicUsize::new(0));
        let target = Arc::clone(&ticks);
        let params = ck::BatchQueryParams {
            progress_callback: Some(Arc::new(move || {
                target.fetch_add(1, Ordering::SeqCst);
            })),
            ..Default::default()
        };
        let write = ck::SmilesWriteParams {
            all_bonds_explicit: true,
            ..Default::default()
        };
        assert_eq!(
            batch.to_smiles_list_with_params(&write, &params).unwrap(),
            source
                .to_smiles_list_with_params(&write, &ck::BatchQueryParams::default())
                .unwrap()
        );
        assert_eq!(ticks.swap(0, Ordering::SeqCst), 3);
        assert_eq!(
            batch.dg_bounds_matrix_list_with_params(&params).unwrap(),
            source
                .dg_bounds_matrix_list_with_params(&ck::BatchQueryParams::default())
                .unwrap()
        );
        assert_eq!(ticks.swap(0, Ordering::SeqCst), 3);
        assert_eq!(
            batch.to_svg_list_with_params(120, 100, &params).unwrap(),
            source
                .to_svg_list_with_params(120, 100, &ck::BatchQueryParams::default())
                .unwrap()
        );
        assert_eq!(ticks.swap(0, Ordering::SeqCst), 3);
        let failure = batch.to_svg_list_with_params(0, 100, &params).unwrap_err();
        assert_eq!(failure.errors, 2);
        assert_eq!(
            failure
                .record_errors
                .iter()
                .map(|e| e.index)
                .collect::<Vec<_>>(),
            [0, 2]
        );
        assert_eq!(ticks.load(Ordering::SeqCst), 3);
        assert_eq!(batch.valid_mask(), [true, false, true]);
        assert!(batch.to_smiles_list().unwrap()[1].is_none());
        assert!(
            MoleculeBatch::from_smiles_list(&[])
                .unwrap()
                .dg_bounds_matrix_list()
                .unwrap()
                .is_empty()
        );
    }
}

#[cfg(feature = "depict")]
impl MoleculeBatch {
    pub fn to_images(
        &self,
        directory: &str,
    ) -> Result<ck::BatchExportReport, ck::BatchValidationError> {
        self.inner.to_images(std::path::Path::new(directory))
    }
    pub fn to_images_with_params(
        &self,
        directory: &str,
        options: &ck::BatchImageParams,
    ) -> Result<ck::BatchExportReport, ck::BatchValidationError> {
        self.inner
            .to_images_with_params(std::path::Path::new(directory), options)
    }
}

#[cfg(test)]
mod image_tests {
    use super::*;
    #[test]
    fn batch_images_write_real_png_svg_and_source_reports_with_typed_failures() {
        let base = std::path::PathBuf::from(
            "/home/datahouse-raid5/wjt/COSMolKit_p3/target/agent-handoff/wasm-registry-binding-v1",
        )
        .join(format!(
            "batch-images-native-fixtures-{}",
            std::process::id()
        ));
        assert!(!base.exists());
        let valid = MoleculeBatch::from_smiles_list(&["CCO".into()]).unwrap();
        let report = valid
            .to_images(base.join("default").to_str().unwrap())
            .unwrap();
        assert_eq!(report.total(), 1);
        assert_eq!(report.success(), 1);
        assert_eq!(report.failed(), 0);
        assert!(report.errors().is_empty());
        assert!(
            std::fs::read(base.join("default/mol_0.png"))
                .unwrap()
                .starts_with(b"\x89PNG\r\n\x1a\n")
        );
        let keep = ck::BatchParams {
            errors: ck::BatchErrorMode::KeepErrors,
            ..Default::default()
        };
        let partial = MoleculeBatch::from_smiles_list_with_params(
            &["CCO".into(), "[".into(), "O".into()],
            &ck::SmilesParseParams::default(),
            &keep,
        )
        .unwrap();
        let options = ck::BatchImageParams {
            format: "svg".into(),
            width: 120,
            height: 100,
            execution: keep,
            filenames: Some(vec![Some("ethanol".into()), None, Some("water.svg".into())]),
            report_path: Some(base.join("counts.json")),
        };
        let report = partial
            .to_images_with_params(base.join("custom").to_str().unwrap(), &options)
            .unwrap();
        assert_eq!(
            (
                report.total(),
                report.success(),
                report.failed(),
                report.skipped
            ),
            (3, 2, 0, 1)
        );
        for path in ["custom/ethanol.svg", "custom/water.svg"] {
            assert!(
                std::fs::read_to_string(base.join(path))
                    .unwrap()
                    .contains("<svg")
            );
        }
        assert!(!base.join("custom/mol_1.svg").exists());
        assert!(
            std::fs::read_to_string(base.join("counts.json"))
                .unwrap()
                .contains("\"written\": 2")
        );
        report.write_report(&base.join("counts.CSV")).unwrap();
        assert_eq!(
            std::fs::read_to_string(base.join("counts.CSV")).unwrap(),
            "written,skipped,failed\n2,1,0\n"
        );
        assert!(
            report
                .write_report(&base.join("missing/counts.json"))
                .is_err()
        );
        let invalid = ck::BatchImageParams {
            format: "gif".into(),
            execution: keep,
            ..Default::default()
        };
        let failure = partial
            .to_images_with_params(base.join("bad-format").to_str().unwrap(), &invalid)
            .unwrap();
        assert_eq!(
            (failure.written, failure.skipped, failure.failed),
            (0, 1, 2)
        );
        assert!(
            std::error::Error::source(&failure.errors[0])
                .unwrap()
                .downcast_ref::<ck::BatchImageError>()
                .is_some()
        );
        assert_eq!(
            partial
                .to_images(base.join("strict").to_str().unwrap())
                .err()
                .unwrap()
                .errors,
            1
        );
        let invalid = ck::BatchImageParams {
            execution: ck::BatchParams {
                n_jobs: Some(0),
                ..Default::default()
            },
            ..Default::default()
        };
        assert!(
            valid
                .to_images_with_params(base.join("zero").to_str().unwrap(), &invalid)
                .is_err()
        );
        assert!(!base.join("zero").exists());
        let invalid = ck::BatchImageParams {
            filenames: Some(vec![]),
            ..Default::default()
        };
        assert!(
            valid
                .to_images_with_params(base.join("wrong-names").to_str().unwrap(), &invalid)
                .is_err()
        );
    }
}

#[cfg(feature = "fingerprints")]
impl MoleculeBatch {
    pub fn fingerprint_atom_pair_list(
        &self,
    ) -> Result<Vec<Option<ck::Fingerprint>>, ck::BatchValidationError> {
        self.inner.fingerprint_atom_pair_list()
    }
    pub fn fingerprint_atom_pair_list_with_params(
        &self,
        options: &ck::AtomPairFingerprintParams,
        params: &ck::BatchQueryParams,
    ) -> Result<Vec<Option<ck::Fingerprint>>, ck::BatchValidationError> {
        self.inner
            .fingerprint_atom_pair_list_with_params(options, params)
    }
    pub fn fingerprint_atom_pair_sparse_count_list(
        &self,
    ) -> Result<Vec<Option<ck::SparseCountFingerprint>>, ck::BatchValidationError> {
        self.inner.fingerprint_atom_pair_sparse_count_list()
    }
    pub fn fingerprint_atom_pair_sparse_count_list_with_params(
        &self,
        options: &ck::AtomPairFingerprintParams,
        params: &ck::BatchQueryParams,
    ) -> Result<Vec<Option<ck::SparseCountFingerprint>>, ck::BatchValidationError> {
        self.inner
            .fingerprint_atom_pair_sparse_count_list_with_params(options, params)
    }
    pub fn fingerprint_atom_pair_count_list(
        &self,
    ) -> Result<Vec<Option<ck::SparseCountFingerprint32>>, ck::BatchValidationError> {
        self.inner.fingerprint_atom_pair_count_list()
    }
    pub fn fingerprint_atom_pair_count_list_with_params(
        &self,
        options: &ck::AtomPairFingerprintParams,
        params: &ck::BatchQueryParams,
    ) -> Result<Vec<Option<ck::SparseCountFingerprint32>>, ck::BatchValidationError> {
        self.inner
            .fingerprint_atom_pair_count_list_with_params(options, params)
    }
    pub fn fingerprint_atom_pair_sparse_bits_list(
        &self,
    ) -> Result<Vec<Option<ck::SparseBitFingerprint>>, ck::BatchValidationError> {
        self.inner.fingerprint_atom_pair_sparse_bits_list()
    }
    pub fn fingerprint_atom_pair_sparse_bits_list_with_params(
        &self,
        options: &ck::AtomPairFingerprintParams,
        params: &ck::BatchQueryParams,
    ) -> Result<Vec<Option<ck::SparseBitFingerprint>>, ck::BatchValidationError> {
        self.inner
            .fingerprint_atom_pair_sparse_bits_list_with_params(options, params)
    }
    pub fn fingerprint_atom_pair_with_output_list(
        &self,
    ) -> Result<Vec<Option<ck::BatchFingerprintOutput>>, ck::BatchValidationError> {
        self.inner.fingerprint_atom_pair_with_output_list()
    }
    pub fn fingerprint_atom_pair_with_output_list_with_params(
        &self,
        options: &ck::AtomPairFingerprintParams,
        collect_additional_output: bool,
        params: &ck::BatchQueryParams,
    ) -> Result<Vec<Option<ck::BatchFingerprintOutput>>, ck::BatchValidationError> {
        self.inner
            .fingerprint_atom_pair_with_output_list_with_params(
                options,
                collect_additional_output,
                params,
            )
    }
}

#[cfg(all(test, feature = "fingerprints"))]
mod atom_pair_tests {
    use super::*;
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };
    #[test]
    fn every_atom_pair_batch_projection_preserves_scalar_values_null_rows_and_output_policy() {
        let batch = MoleculeBatch::from_smiles_list_with_params(
            &["CCO".into(), "[".into(), "O".into()],
            &ck::SmilesParseParams::default(),
            &ck::BatchParams {
                errors: ck::BatchErrorMode::KeepErrors,
                ..Default::default()
            },
        )
        .unwrap();
        let options = ck::AtomPairFingerprintParams::default();
        let execution = ck::BatchQueryParams::default();
        let scalar = ck::Molecule::from_smiles("CCO").unwrap();
        let bits = batch.fingerprint_atom_pair_list().unwrap();
        assert_eq!(
            bits[0].as_ref().unwrap(),
            &scalar.atom_pair_fingerprint().unwrap()
        );
        assert!(bits[1].is_none());
        assert_eq!(
            bits,
            batch
                .fingerprint_atom_pair_list_with_params(&options, &execution)
                .unwrap()
        );
        let sparse = batch.fingerprint_atom_pair_sparse_count_list().unwrap();
        assert_eq!(
            sparse[0].as_ref().unwrap(),
            &scalar.atom_pair_sparse_count_fingerprint().unwrap()
        );
        assert_eq!(
            sparse,
            batch
                .fingerprint_atom_pair_sparse_count_list_with_params(&options, &execution)
                .unwrap()
        );
        let count = batch.fingerprint_atom_pair_count_list().unwrap();
        assert_eq!(
            count,
            batch
                .fingerprint_atom_pair_count_list_with_params(&options, &execution)
                .unwrap()
        );
        let sparse_bits = batch.fingerprint_atom_pair_sparse_bits_list().unwrap();
        assert_eq!(
            sparse_bits,
            batch
                .fingerprint_atom_pair_sparse_bits_list_with_params(&options, &execution)
                .unwrap()
        );
        let outputs = batch.fingerprint_atom_pair_with_output_list().unwrap();
        assert_eq!(
            outputs[0].as_ref().unwrap().fingerprint(),
            bits[0].as_ref().unwrap()
        );
        assert_eq!(
            outputs,
            batch
                .fingerprint_atom_pair_with_output_list_with_params(&options, true, &execution)
                .unwrap()
        );
        let additional = outputs[0].as_ref().unwrap().additional_output().unwrap();
        assert_eq!(additional.atom_counts().unwrap().len(), 3);
        assert!(additional.bit_paths().is_none());
        let no_output = batch
            .fingerprint_atom_pair_with_output_list_with_params(&options, false, &execution)
            .unwrap();
        assert!(matches!(
            no_output[0].as_ref().unwrap().additional_output(),
            Err(ck::BatchFingerprintOutputError::MissingAdditionalOutput {
                fingerprint_kind: "AtomPair"
            })
        ));
        let ticks = Arc::new(AtomicUsize::new(0));
        let observed = ticks.clone();
        let callback = ck::BatchQueryParams {
            progress_callback: Some(Arc::new(move || {
                observed.fetch_add(1, Ordering::SeqCst);
            })),
            ..Default::default()
        };
        batch
            .fingerprint_atom_pair_list_with_params(&options, &callback)
            .unwrap();
        assert_eq!(ticks.load(Ordering::SeqCst), 3);
        assert_eq!(batch.valid_mask(), [true, false, true]);
    }
}

#[cfg(feature = "fingerprints")]
impl MoleculeBatch {
    pub fn fingerprint_layered_list(
        &self,
    ) -> Result<Vec<Option<ck::Fingerprint>>, ck::BatchValidationError> {
        self.inner.fingerprint_layered_list()
    }
    pub fn fingerprint_layered_list_with_params(
        &self,
        options: &ck::LayeredFingerprintParams,
        params: &ck::BatchQueryParams,
    ) -> Result<Vec<Option<ck::Fingerprint>>, ck::BatchValidationError> {
        self.inner
            .fingerprint_layered_list_with_params(options, params)
    }
    pub fn fingerprint_layered_with_output_list(
        &self,
    ) -> Result<Vec<Option<ck::LayeredFingerprintResult>>, ck::BatchValidationError> {
        self.inner.fingerprint_layered_with_output_list()
    }
    pub fn fingerprint_layered_with_output_list_with_params(
        &self,
        options: &ck::LayeredFingerprintParams,
        params: &ck::BatchQueryParams,
    ) -> Result<Vec<Option<ck::LayeredFingerprintResult>>, ck::BatchValidationError> {
        self.inner
            .fingerprint_layered_with_output_list_with_params(options, params)
    }
    pub fn pattern_fingerprint_list(
        &self,
    ) -> Result<Vec<Option<ck::Fingerprint>>, ck::BatchValidationError> {
        self.inner.pattern_fingerprint_list()
    }
    pub fn pattern_fingerprint_list_with_params(
        &self,
        options: &ck::PatternFingerprintParams,
        params: &ck::BatchQueryParams,
    ) -> Result<Vec<Option<ck::Fingerprint>>, ck::BatchValidationError> {
        self.inner
            .pattern_fingerprint_list_with_params(options, params)
    }
}

#[cfg(all(test, feature = "fingerprints"))]
mod layered_pattern_tests {
    use super::*;
    #[test]
    fn all_six_layered_pattern_batch_calls_preserve_scalar_results_and_seeded_counts() {
        let batch = MoleculeBatch::from_smiles_list_with_params(
            &["CCO".into(), "[".into(), "O".into()],
            &ck::SmilesParseParams::default(),
            &ck::BatchParams {
                errors: ck::BatchErrorMode::KeepErrors,
                ..Default::default()
            },
        )
        .unwrap();
        let scalar = ck::Molecule::from_smiles("CCO").unwrap();
        let options = ck::LayeredFingerprintParams::default();
        let query = ck::BatchQueryParams::default();
        let dense = batch.fingerprint_layered_list().unwrap();
        assert_eq!(
            dense[0].as_ref().unwrap(),
            &scalar.layered_fingerprint().unwrap()
        );
        assert!(dense[1].is_none());
        assert_eq!(
            dense,
            batch
                .fingerprint_layered_list_with_params(&options, &query)
                .unwrap()
        );
        let output = batch.fingerprint_layered_with_output_list().unwrap();
        assert_eq!(
            output[0].as_ref().unwrap(),
            &scalar.layered_fingerprint_with_output().unwrap()
        );
        assert_eq!(
            output,
            batch
                .fingerprint_layered_with_output_list_with_params(&options, &query)
                .unwrap()
        );
        assert!(output[0].as_ref().unwrap().atom_counts().is_none());
        let pattern = batch.pattern_fingerprint_list().unwrap();
        assert_eq!(
            pattern[0].as_ref().unwrap(),
            &scalar.pattern_fingerprint().unwrap()
        );
        assert_eq!(
            pattern,
            batch
                .pattern_fingerprint_list_with_params(
                    &ck::PatternFingerprintParams::default(),
                    &query
                )
                .unwrap()
        );
        let single = MoleculeBatch::from_smiles_list(&["CCO".into()]).unwrap();
        let seeded = ck::LayeredFingerprintParams {
            atom_counts: Some(vec![10, 20, 30]),
            set_only_bits: Some(ck::Fingerprint::from_on_bits(2048, [674]).unwrap()),
            ..Default::default()
        };
        let result = single
            .fingerprint_layered_with_output_list_with_params(&seeded, &query)
            .unwrap();
        assert_eq!(
            result[0].as_ref().unwrap(),
            &scalar
                .layered_fingerprint_with_output_with_params(&seeded)
                .unwrap()
        );
        assert_eq!(
            result[0].as_ref().unwrap().atom_counts(),
            Some([11, 22, 31].as_slice())
        );
        assert_eq!(seeded.atom_counts, Some(vec![10, 20, 30]));
        assert_eq!(batch.valid_mask(), [true, false, true]);
    }
}

#[cfg(feature = "fingerprints")]
impl MoleculeBatch {
    pub fn fingerprint_morgan_list(
        &self,
    ) -> Result<Vec<Option<ck::Fingerprint>>, ck::BatchValidationError> {
        self.inner.fingerprint_morgan_list()
    }
    pub fn fingerprint_morgan_list_with_params(
        &self,
        options: &ck::MorganFingerprintParams,
        params: &ck::BatchQueryParams,
    ) -> Result<Vec<Option<ck::Fingerprint>>, ck::BatchValidationError> {
        self.inner
            .fingerprint_morgan_list_with_params(options, params)
    }
    pub fn fingerprint_morgan_with_output_list(
        &self,
    ) -> Result<Vec<Option<ck::BatchFingerprintOutput>>, ck::BatchValidationError> {
        self.inner.fingerprint_morgan_with_output_list()
    }
    pub fn fingerprint_morgan_with_output_list_with_params(
        &self,
        options: &ck::MorganFingerprintParams,
        collect_additional_output: bool,
        params: &ck::BatchQueryParams,
    ) -> Result<Vec<Option<ck::BatchFingerprintOutput>>, ck::BatchValidationError> {
        self.inner.fingerprint_morgan_with_output_list_with_params(
            options,
            collect_additional_output,
            params,
        )
    }
    pub fn fingerprint_morgan_list_with_generator_params(
        &self,
        options: &ck::MorganParams,
        atom_invariants: Option<&ck::MorganAtomInvariantsGenerator>,
        bond_invariants: Option<&ck::MorganBondInvariantsGenerator>,
        call: &ck::MorganCallParams,
        params: &ck::BatchQueryParams,
    ) -> Result<Vec<Option<ck::Fingerprint>>, ck::BatchValidationError> {
        self.inner.fingerprint_morgan_list_with_generator_params(
            options,
            atom_invariants,
            bond_invariants,
            call,
            params,
        )
    }
    pub fn fingerprint_morgan_with_output_list_with_generator_params(
        &self,
        options: &ck::MorganParams,
        atom_invariants: Option<&ck::MorganAtomInvariantsGenerator>,
        bond_invariants: Option<&ck::MorganBondInvariantsGenerator>,
        call: &ck::MorganCallParams,
        collect_additional_output: bool,
        params: &ck::BatchQueryParams,
    ) -> Result<Vec<Option<ck::BatchFingerprintOutput>>, ck::BatchValidationError> {
        self.inner
            .fingerprint_morgan_with_output_list_with_generator_params(
                options,
                atom_invariants,
                bond_invariants,
                call,
                collect_additional_output,
                params,
            )
    }
}

#[cfg(all(test, feature = "fingerprints"))]
mod morgan_tests {
    use super::*;
    #[test]
    fn all_six_morgan_batch_calls_preserve_scalar_values_providers_and_optional_output() {
        let batch = MoleculeBatch::from_smiles_list_with_params(
            &["CCO".into(), "[".into(), "O".into()],
            &ck::SmilesParseParams::default(),
            &ck::BatchParams {
                errors: ck::BatchErrorMode::KeepErrors,
                ..Default::default()
            },
        )
        .unwrap();
        let scalar = ck::Molecule::from_smiles("CCO").unwrap();
        let generator = ck::MorganParams {
            radius: 2,
            ..Default::default()
        };
        let options = ck::MorganFingerprintParams {
            generator: generator.clone(),
            ..Default::default()
        };
        let call = ck::MorganCallParams::default();
        let query = ck::BatchQueryParams::default();
        let dense = batch.fingerprint_morgan_list().unwrap();
        assert_eq!(
            dense[0].as_ref().unwrap(),
            &scalar
                .morgan_fingerprint_with_params(&options, None)
                .unwrap()
        );
        assert!(dense[1].is_none());
        assert_eq!(
            dense,
            batch
                .fingerprint_morgan_list_with_params(&options, &query)
                .unwrap()
        );
        assert_eq!(
            dense,
            batch
                .fingerprint_morgan_list_with_generator_params(
                    &generator, None, None, &call, &query
                )
                .unwrap()
        );
        let output = batch.fingerprint_morgan_with_output_list().unwrap();
        assert_eq!(
            output,
            batch
                .fingerprint_morgan_with_output_list_with_params(&options, true, &query)
                .unwrap()
        );
        assert_eq!(
            output,
            batch
                .fingerprint_morgan_with_output_list_with_generator_params(
                    &generator, None, None, &call, true, &query
                )
                .unwrap()
        );
        assert_eq!(
            output[0].as_ref().unwrap().fingerprint(),
            dense[0].as_ref().unwrap()
        );
        assert_eq!(
            output[0]
                .as_ref()
                .unwrap()
                .additional_output()
                .unwrap()
                .atom_counts()
                .unwrap()
                .len(),
            3
        );
        let missing = batch
            .fingerprint_morgan_with_output_list_with_params(&options, false, &query)
            .unwrap();
        assert!(matches!(
            missing[0].as_ref().unwrap().additional_output(),
            Err(ck::BatchFingerprintOutputError::MissingAdditionalOutput {
                fingerprint_kind: "Morgan"
            })
        ));
        let atom = ck::MorganAtomInvariantsGenerator::features(Some(vec![
            ck::search::from_smarts("[#8]").unwrap(),
        ]));
        let bond = ck::MorganBondInvariantsGenerator::new(false, true);
        let configured =
            ck::MorganFingerprintGenerator::new(Some(&generator), Some(&atom), Some(&bond))
                .unwrap();
        let values = batch
            .fingerprint_morgan_list_with_generator_params(
                &generator,
                Some(&atom),
                Some(&bond),
                &call,
                &query,
            )
            .unwrap();
        assert_eq!(
            values[0].as_ref().unwrap(),
            &scalar
                .morgan_fingerprint_with_generator(&configured, Some(&call), None)
                .unwrap()
        );
        assert_eq!(batch.valid_mask(), [true, false, true]);
    }
}
