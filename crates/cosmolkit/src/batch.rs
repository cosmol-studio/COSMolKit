//! The sole live batch value; domain chemistry is invoked only through public scalar APIs.
#[cfg(feature = "cap-fingerprints")]
use crate::BatchFingerprintOutput;
use crate::Molecule;
pub use cosmolkit_batch::{BatchErrorMode, BatchRecordError as BatchError, BatchValidationError};

/// Per-call scheduling and error policy; overrides never mutate stored batch configuration.
#[derive(Debug, Clone, Copy, Default)]
pub struct BatchParams {
    /// None inherits the batch policy; constructors without a batch use Strict.
    pub errors: Option<BatchErrorMode>,
    /// Inherit a stored batch setting when absent, otherwise default to one worker.
    pub n_jobs: Option<usize>,
    pub progress_bar: Option<bool>,
}
#[derive(Debug, Clone)]
pub enum BatchRecord {
    Molecule(Molecule),
    Error(BatchError),
}
#[derive(Debug, Clone, Default)]
pub struct MoleculeBatch {
    pub(super) records: Vec<BatchRecord>,
    pub(super) error_mode: BatchErrorMode,
    pub(super) n_jobs: Option<usize>,
    pub(super) progress_bar: Option<bool>,
}

impl MoleculeBatch {
    pub fn from_records(
        records: Vec<BatchRecord>,
        mode: BatchErrorMode,
    ) -> Result<Self, BatchValidationError> {
        let errors = records
            .iter()
            .filter_map(|r| match r {
                BatchRecord::Error(e) => Some(e.clone()),
                _ => None,
            })
            .collect();
        cosmolkit_batch::validate_record_errors(errors, mode)?;
        Ok(Self {
            records,
            error_mode: mode,
            n_jobs: None,
            progress_bar: None,
        })
    }
    pub fn len(&self) -> usize {
        self.records.len()
    }
    /// The error policy inherited by subsequent transforms and exports.
    pub fn error_mode(&self) -> BatchErrorMode {
        self.error_mode
    }
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
    pub fn get(&self, index: usize) -> Option<&BatchRecord> {
        self.records.get(index)
    }
    pub fn records(&self) -> &[BatchRecord] {
        &self.records
    }
    pub fn valid_mask(&self) -> Vec<bool> {
        self.records
            .iter()
            .map(|r| matches!(r, BatchRecord::Molecule(_)))
            .collect()
    }
    pub fn errors(&self) -> Vec<BatchError> {
        self.records
            .iter()
            .filter_map(|r| match r {
                BatchRecord::Error(e) => Some(e.clone()),
                _ => None,
            })
            .collect()
    }
    pub fn valid_count(&self) -> usize {
        self.records
            .iter()
            .filter(|r| matches!(r, BatchRecord::Molecule(_)))
            .count()
    }
    pub fn invalid_count(&self) -> usize {
        self.len() - self.valid_count()
    }
    pub fn with_parallel_jobs(
        mut self,
        n_jobs: Option<usize>,
    ) -> Result<Self, BatchValidationError> {
        if n_jobs == Some(0) {
            return Err(BatchValidationError::parameter(
                "n_jobs",
                "n_jobs must be >= 1",
            ));
        }
        self.n_jobs = n_jobs;
        Ok(self)
    }
    pub fn parallel_jobs(&self) -> Option<usize> {
        self.n_jobs
    }
    pub fn with_progress_bar(mut self, progress_bar: Option<bool>) -> Self {
        self.progress_bar = progress_bar;
        self
    }
    pub fn progress_bar(&self) -> Option<bool> {
        self.progress_bar
    }
    pub fn invalid_mask(&self) -> Vec<bool> {
        self.valid_mask().into_iter().map(|valid| !valid).collect()
    }
    pub fn iter(&self) -> std::slice::Iter<'_, BatchRecord> {
        self.records.iter()
    }
    pub fn to_list(&self) -> Vec<Option<Molecule>> {
        self.records
            .iter()
            .map(|record| match record {
                BatchRecord::Molecule(molecule) => Some(molecule.clone()),
                BatchRecord::Error(_) => None,
            })
            .collect()
    }
    pub fn with_valid_records(&self) -> Self {
        // pub fn filter_valid(&self) -> Self {
        //         Self {
        //             records: self
        //                 .records
        //                 .iter()
        //                 .filter(|record| matches!(record, BatchRecord::Molecule(_)))
        //                 .cloned()
        //                 .collect(),
        //             n_jobs: self.n_jobs,
        //             progress_bar: self.progress_bar,
        //         }
        //     }

        Self {
            records: self
                .records
                .iter()
                .filter(|r| matches!(r, BatchRecord::Molecule(_)))
                .cloned()
                .collect(),
            n_jobs: self.n_jobs,
            progress_bar: self.progress_bar,
            error_mode: self.error_mode,
        }
    }
    #[cfg(feature = "cap-smiles")]
    pub fn from_smiles_list(smiles: &[String]) -> Result<Self, BatchValidationError> {
        Self::from_smiles_list_with_params(
            smiles,
            &crate::SmilesParseParams::default(),
            &BatchParams::default(),
        )
    }
    #[cfg(feature = "cap-smiles")]
    pub fn from_smiles_list_with_params(
        smiles: &[String],
        parse: &crate::SmilesParseParams,
        params: &BatchParams,
    ) -> Result<Self, BatchValidationError> {
        // pub fn from_smiles_list_with_sanitize_and_options(
        //         smiles: &[String],
        //         sanitize: bool,
        //         errors: BatchErrorMode,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> Result<Self, BatchValidationError> {
        //         with_progress_bar_for_option(progress_bar, smiles.len(), "Parsing SMILES", |progress| {
        //             let records = run_with_parallel_jobs_option(n_jobs, || {
        //                 smiles
        //                     .par_iter()
        //                     .enumerate()
        //                     .map(|(index, smiles)| {
        //                         let out = match Molecule::from_smiles_with_sanitize(smiles, sanitize) {
        //                             Ok(molecule) => BatchRecord::Molecule(molecule),
        //                             Err(error) => BatchRecord::Error(BatchRecordError::new(
        //                                 index,
        //                                 "batch.from_smiles_list",
        //                                 error.to_string(),
        //                             )),
        //                         };
        //                         tick_progress(progress);
        //                         out
        //                     })
        //                     .collect()
        //             });
        //             Self::from_records_with_mode(records, errors)
        //         })
        //     }

        let records = cosmolkit_batch::run_indexed(
            smiles.len(),
            params.n_jobs,
            params.progress_bar,
            "Parsing SMILES",
            |index| match Molecule::from_smiles_with_params(&smiles[index], parse) {
                Ok(molecule) => BatchRecord::Molecule(molecule),
                Err(error) => BatchRecord::Error(BatchError::with_source(
                    index,
                    "batch.from_smiles_list",
                    error,
                )),
            },
        )?;
        Self::from_records(records, params.errors.unwrap_or_default())
    }
    fn transform<E: std::error::Error + Send + Sync + 'static>(
        &self,
        operation: &'static str,
        message: &'static str,
        params: &BatchParams,
        apply: impl Fn(&Molecule) -> Result<Molecule, E> + Sync + Send,
    ) -> Result<Self, BatchValidationError> {
        // fn transform_with_options<F>(
        //         &self,
        //         operation: &'static str,
        //         errors: BatchErrorMode,
        //         transform: F,
        //         progress: BatchProgress<'_>,
        //         n_jobs: Option<usize>,
        //     ) -> Result<Self, BatchValidationError>
        //     where
        //         F: Fn(&Molecule) -> Result<Molecule, crate::OperationError> + Sync + Send,
        //     {
        //         let records: Vec<BatchRecord> = self.run_with_parallel_jobs(n_jobs, || {
        //             self.records
        //                 .par_iter()
        //                 .enumerate()
        //                 .map(|(index, record)| {
        //                     let out = match record {
        //                         BatchRecord::Molecule(molecule) => match transform(molecule) {
        //                             Ok(molecule) => BatchRecord::Molecule(molecule),
        //                             Err(error) => BatchRecord::Error(BatchRecordError::new(
        //                                 index,
        //                                 operation,
        //                                 error.to_string(),
        //                             )),
        //                         },
        //                         BatchRecord::Error(error) => BatchRecord::Error(error.clone()),
        //                     };
        //                     tick_progress(progress);
        //                     out
        //                 })
        //                 .collect()
        //         });
        //         if errors.raise_on_errors() {
        //             let record_errors = record_errors_from_records(&records);
        //             if !record_errors.is_empty() {
        //                 return Err(BatchValidationError::from_record_errors(record_errors));
        //             }
        //         }
        //         Ok(Self {
        //             records,
        //             n_jobs: self.n_jobs,
        //             progress_bar: self.progress_bar,
        //         })
        //     }

        let records = cosmolkit_batch::run_indexed(
            self.len(),
            params.n_jobs.or(self.n_jobs),
            params.progress_bar.or(self.progress_bar),
            message,
            |index| match &self.records[index] {
                BatchRecord::Molecule(molecule) => match apply(molecule) {
                    Ok(molecule) => BatchRecord::Molecule(molecule),
                    Err(error) => {
                        BatchRecord::Error(BatchError::with_source(index, operation, error))
                    }
                },
                BatchRecord::Error(error) => BatchRecord::Error(error.clone()),
            },
        )?;
        // Explicit overrides apply to the returned chain, never to the source batch.
        let mut result = Self::from_records(records, params.errors.unwrap_or(self.error_mode))?;
        result.n_jobs = self.n_jobs;
        result.progress_bar = self.progress_bar;
        Ok(result)
    }
    #[cfg(feature = "cap-sanitize")]
    pub fn sanitize(&self) -> Result<Self, BatchValidationError> {
        self.sanitize_with_params(&crate::SanitizeParams::default(), &BatchParams::default())
    }
    #[cfg(feature = "cap-sanitize")]
    pub fn sanitize_with_params(
        &self,
        options: &crate::SanitizeParams,
        params: &BatchParams,
    ) -> Result<Self, BatchValidationError> {
        // pub fn sanitize_with_options(
        //         &self,
        //         errors: BatchErrorMode,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> Result<Self, BatchValidationError> {
        //         self.with_progress_bar_for(
        //             progress_bar,
        //             self.records.len(),
        //             "Sanitizing molecules",
        //             |progress| {
        //                 self.transform_with_options(
        //                     "batch.sanitize",
        //                     errors,
        //                     Molecule::sanitize,
        //                     progress,
        //                     n_jobs,
        //                 )
        //             },
        //         )
        //     }
        self.transform(
            "batch.sanitize",
            "Sanitizing molecules",
            params,
            |molecule| molecule.sanitize_with_params(options),
        )
    }
    #[cfg(feature = "cap-hydrogens")]
    pub fn with_hydrogens(&self) -> Result<Self, BatchValidationError> {
        self.with_hydrogens_with_params(&crate::AddHsParams::default(), &BatchParams::default())
    }
    #[cfg(feature = "cap-hydrogens")]
    pub fn with_hydrogens_with_params(
        &self,
        options: &crate::AddHsParams,
        params: &BatchParams,
    ) -> Result<Self, BatchValidationError> {
        // pub fn with_hydrogens_with_options(
        //         &self,
        //         errors: BatchErrorMode,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> Result<Self, BatchValidationError> {
        //         self.with_progress_bar_for(
        //             progress_bar,
        //             self.records.len(),
        //             "Adding hydrogens",
        //             |progress| {
        //                 self.transform_with_options(
        //                     "batch.with_hydrogens",
        //                     errors,
        //                     Molecule::with_hydrogens,
        //                     progress,
        //                     n_jobs,
        //                 )
        //             },
        //         )
        //     }
        self.transform(
            "batch.with_hydrogens",
            "Adding hydrogens",
            params,
            |molecule| molecule.with_hydrogens_with_params(options),
        )
    }
    #[cfg(feature = "cap-hydrogens")]
    pub fn without_hydrogens(&self) -> Result<Self, BatchValidationError> {
        self.without_hydrogens_with_params(
            &crate::RemoveHsParams::default(),
            &BatchParams::default(),
        )
    }
    #[cfg(feature = "cap-hydrogens")]
    pub fn without_hydrogens_with_params(
        &self,
        options: &crate::RemoveHsParams,
        params: &BatchParams,
    ) -> Result<Self, BatchValidationError> {
        // pub fn without_hydrogens_with_options(
        //         &self,
        //         errors: BatchErrorMode,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> Result<Self, BatchValidationError> {
        //         self.with_progress_bar_for(
        //             progress_bar,
        //             self.records.len(),
        //             "Removing hydrogens",
        //             |progress| {
        //                 self.transform_with_options(
        //                     "batch.without_hydrogens",
        //                     errors,
        //                     Molecule::without_hydrogens,
        //                     progress,
        //                     n_jobs,
        //                 )
        //             },
        //         )
        //     }
        self.transform(
            "batch.without_hydrogens",
            "Removing hydrogens",
            params,
            |molecule| molecule.without_hydrogens_with_params(options),
        )
    }
    #[cfg(feature = "cap-kekulize")]
    pub fn with_kekulized_bonds(&self) -> Result<Self, BatchValidationError> {
        self.with_kekulized_bonds_with_params(
            &crate::KekulizeParams::default(),
            &BatchParams::default(),
        )
    }
    #[cfg(feature = "cap-kekulize")]
    pub fn with_kekulized_bonds_with_params(
        &self,
        options: &crate::KekulizeParams,
        params: &BatchParams,
    ) -> Result<Self, BatchValidationError> {
        // pub fn with_kekulized_bonds_with_options(
        //         &self,
        //         clear_aromatic_flags: bool,
        //         errors: BatchErrorMode,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> Result<Self, BatchValidationError> {
        //         self.with_progress_bar_for(
        //             progress_bar,
        //             self.records.len(),
        //             "Kekulizing molecules",
        //             |progress| {
        //                 self.transform_with_options(
        //                     "batch.with_kekulized_bonds",
        //                     errors,
        //                     move |molecule| molecule.with_kekulized_bonds(clear_aromatic_flags),
        //                     progress,
        //                     n_jobs,
        //                 )
        //             },
        //         )
        //     }
        self.transform(
            "batch.with_kekulized_bonds",
            "Kekulizing molecules",
            params,
            |molecule| molecule.with_kekulized_bonds_with_params(options),
        )
    }
    #[cfg(feature = "cap-depict")]
    pub fn with_2d_coordinates(&self) -> Result<Self, BatchValidationError> {
        self.with_2d_coordinates_with_params(
            &crate::Coordinate2DParams::default(),
            &BatchParams::default(),
        )
    }
    #[cfg(feature = "cap-depict")]
    pub fn with_2d_coordinates_with_params(
        &self,
        options: &crate::Coordinate2DParams,
        params: &BatchParams,
    ) -> Result<Self, BatchValidationError> {
        // pub fn with_2d_coordinates_with_options(
        //         &self,
        //         errors: BatchErrorMode,
        //         n_jobs: Option<usize>,
        //         progress_bar: Option<bool>,
        //     ) -> Result<Self, BatchValidationError> {
        //         self.with_progress_bar_for(
        //             progress_bar,
        //             self.records.len(),
        //             "Computing 2D coordinates",
        //             |progress| {
        //                 self.transform_with_options(
        //                     "batch.with_2d_coordinates",
        //                     errors,
        //                     Molecule::with_2d_coordinates,
        //                     progress,
        //                     n_jobs,
        //                 )
        //             },
        //         )
        //     }
        self.transform(
            "batch.with_2d_coordinates",
            "Computing 2D coordinates",
            params,
            |molecule| molecule.with_2d_coordinates_with_params(options),
        )
    }
}

/// Runtime overrides for ordered read results; newly failed queries always raise.
#[derive(Clone, Default)]
pub struct BatchQueryParams {
    /// Inherit a stored batch setting when absent, otherwise default to one worker.
    pub n_jobs: Option<usize>,
    pub progress_bar: Option<bool>,
    /// One notification after every completed input row, including invalid rows.
    pub progress_callback: Option<std::sync::Arc<dyn Fn() + Send + Sync>>,
}
impl std::fmt::Debug for BatchQueryParams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BatchQueryParams")
            .field("n_jobs", &self.n_jobs)
            .field("progress_bar", &self.progress_bar)
            .field("has_progress_callback", &self.progress_callback.is_some())
            .finish()
    }
}
impl MoleculeBatch {
    fn collect<T: Send, E: std::error::Error + Send + Sync + 'static>(
        &self,
        operation: &'static str,
        message: &'static str,
        params: &BatchQueryParams,
        apply: impl Fn(&Molecule) -> Result<T, E> + Sync + Send,
    ) -> Result<Vec<Option<T>>, BatchValidationError> {
        // COSMolKit❗✔️: pinned d892ec3 properties/batch.rs::collect_optional_values_with_options; same indexed order/all-result collection and scalar owner calls.
        // fn collect_optional_values_with_options<T, F>(
        //         &self,
        //         operation: &'static str,
        //         collect: F,
        //         progress: BatchProgress<'_>,
        //         n_jobs: Option<usize>,
        //     ) -> Result<Vec<Option<T>>, BatchValidationError>
        //     where
        //         T: Send,
        //         F: Fn(&Molecule) -> Result<T, String> + Sync + Send,
        //     {
        //         let pairs: Vec<(Option<T>, Option<BatchRecordError>)> =
        //             self.run_with_parallel_jobs(n_jobs, || {
        //                 self.records
        //                     .par_iter()
        //                     .enumerate()
        //                     .map(|(index, record)| {
        //                         let out = match record {
        //                             BatchRecord::Molecule(molecule) => match collect(molecule) {
        //                                 Ok(value) => (Some(value), None),
        //                                 Err(error) => {
        //                                     (None, Some(BatchRecordError::new(index, operation, error)))
        //                                 }
        //                             },
        //                             BatchRecord::Error(_) => (None, None),
        //                         };
        //                         tick_progress(progress);
        //                         out
        //                     })
        //                     .collect()
        //             });
        //         let mut values = Vec::with_capacity(pairs.len());
        //         let mut record_errors = Vec::new();
        //         for (value, error) in pairs {
        //             values.push(value);
        //             if let Some(error) = error {
        //                 record_errors.push(error);
        //             }
        //         }
        //         if !record_errors.is_empty() {
        //             return Err(BatchValidationError::from_record_errors(record_errors));
        //         }
        //         Ok(values)
        //     }

        let outcomes = cosmolkit_batch::run_indexed(
            self.len(),
            params.n_jobs.or(self.n_jobs),
            params.progress_bar.or(self.progress_bar),
            message,
            |index| {
                let result = match &self.records[index] {
                    BatchRecord::Molecule(molecule) => apply(molecule)
                        .map(Some)
                        .map_err(|e| BatchError::with_source(index, operation, e)),
                    BatchRecord::Error(_) => Ok(None),
                };
                if let Some(progress) = &params.progress_callback {
                    progress();
                }
                result
            },
        )?;
        let mut values = Vec::with_capacity(outcomes.len());
        let mut errors = Vec::new();
        for outcome in outcomes {
            match outcome {
                Ok(value) => values.push(value),
                Err(error) => {
                    values.push(None);
                    errors.push(error);
                }
            }
        }
        cosmolkit_batch::validate_record_errors(errors, BatchErrorMode::Strict)?;
        Ok(values)
    }
    #[cfg(feature = "cap-fingerprints")]
    /// Ordered source torsion fingerprints using the batch's stored runtime policy.
    pub fn fingerprint_topological_torsion_list(
        &self,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        self.fingerprint_topological_torsion_list_with_params(
            &crate::TopologicalTorsionFingerprintParams::default(),
            &BatchQueryParams::default(),
        )
    }
    #[cfg(feature = "cap-fingerprints")]
    /// Per-call overrides share the canonical ordered scheduler and typed errors.
    pub fn fingerprint_topological_torsion_list_with_params(
        &self,
        options: &crate::TopologicalTorsionFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        // ROOT CK-d41: project-defined thin batch query. Chemistry remains in
        // Molecule's canonical torsion query and the sole fingerprints owner.
        // One immutable borrow per input; collection/error/tick costs are the
        // accepted shared batch mechanism below, with no chemistry duplication.
        self.collect(
            "batch.fingerprint_topological_torsion",
            "Computing Topological Torsion fingerprints",
            params,
            |m| m.fingerprint_topological_torsion_with_params(options, None),
        )
    }

    #[cfg(feature = "cap-smiles")]
    pub fn to_smiles_list(&self) -> Result<Vec<Option<crate::PropertyText>>, BatchValidationError> {
        self.to_smiles_list_with_params(
            &crate::SmilesWriteParams::default(),
            &BatchQueryParams::default(),
        )
    }
    #[cfg(feature = "cap-smiles")]
    pub fn to_smiles_list_with_params(
        &self,
        options: &crate::SmilesWriteParams,
        params: &BatchQueryParams,
    ) -> Result<Vec<Option<crate::PropertyText>>, BatchValidationError> {
        // COSMolKit❗✔️: pinned d892ec3 properties/batch.rs::to_smiles_optional_list_with_params_and_runtime; same indexed order/all-result collection and scalar owner calls.
        // fn to_smiles_optional_list_with_params_and_runtime(
        //         &self,
        //         params: &SmilesWriteParams,
        //         progress: BatchProgress<'_>,
        //         n_jobs: Option<usize>,
        //     ) -> Result<Vec<Option<String>>, BatchValidationError> {
        //         self.collect_optional_values_with_options(
        //             "to_smiles",
        //             |molecule| {
        //                 molecule
        //                     .to_smiles_with_params(params)
        //                     .map_err(|error| error.to_string())
        //             },
        //             progress,
        //             n_jobs,
        //         )
        //     }

        self.collect("to_smiles", "Writing SMILES", params, |m| {
            m.to_smiles_with_params(options)
        })
    }
    #[cfg(feature = "cap-conformer")]
    pub fn dg_bounds_matrix_list(
        &self,
    ) -> Result<Vec<Option<Vec<Vec<f64>>>>, BatchValidationError> {
        self.dg_bounds_matrix_list_with_params(&BatchQueryParams::default())
    }
    #[cfg(feature = "cap-conformer")]
    pub fn dg_bounds_matrix_list_with_params(
        &self,
        params: &BatchQueryParams,
    ) -> Result<Vec<Option<Vec<Vec<f64>>>>, BatchValidationError> {
        // COSMolKit❗✔️: pinned d892ec3 properties/batch.rs::dg_bounds_matrix_list_with_runtime; same indexed order/all-result collection and scalar owner calls.
        // fn dg_bounds_matrix_list_with_runtime(
        //         &self,
        //         progress: BatchProgress<'_>,
        //         n_jobs: Option<usize>,
        //     ) -> Result<Vec<Option<Vec<Vec<f64>>>>, BatchValidationError> {
        //         self.collect_optional_values_with_options(
        //             "batch.dg_bounds_matrix",
        //             |molecule| {
        //                 molecule
        //                     .dg_bounds_matrix()
        //                     .map_err(|error| error.to_string())
        //             },
        //             progress,
        //             n_jobs,
        //         )
        //     }

        self.collect(
            "batch.dg_bounds_matrix",
            "Computing DG bounds matrices",
            params,
            Molecule::dg_bounds_matrix,
        )
    }
    #[cfg(feature = "cap-depict")]
    pub fn to_svg_list(
        &self,
        width: u32,
        height: u32,
    ) -> Result<Vec<Option<String>>, BatchValidationError> {
        self.to_svg_list_with_params(width, height, &BatchQueryParams::default())
    }
    #[cfg(feature = "cap-depict")]
    pub fn to_svg_list_with_params(
        &self,
        width: u32,
        height: u32,
        params: &BatchQueryParams,
    ) -> Result<Vec<Option<String>>, BatchValidationError> {
        // COSMolKit❗✔️: pinned d892ec3 properties/batch.rs::to_svg_list_with_runtime; same indexed order/all-result collection and scalar owner calls.
        // fn to_svg_list_with_runtime(
        //         &self,
        //         width: u32,
        //         height: u32,
        //         progress: BatchProgress<'_>,
        //         n_jobs: Option<usize>,
        //     ) -> Result<Vec<Option<String>>, BatchValidationError> {
        //         self.collect_optional_values_with_options(
        //             "batch.to_svg",
        //             |molecule| {
        //                 molecule
        //                     .to_svg(width, height)
        //                     .map_err(|error| error.to_string())
        //             },
        //             progress,
        //             n_jobs,
        //         )
        //     }

        self.collect("batch.to_svg", "Drawing SVG molecules", params, |m| {
            m.to_svg(width, height)
        })
    }
}

#[cfg(feature = "cap-depict")]
#[derive(Debug, Clone)]
pub struct BatchImageParams {
    pub format: String,
    pub width: u32,
    pub height: u32,
    pub execution: BatchParams,
    pub filenames: Option<Vec<Option<String>>>,
    pub report_path: Option<std::path::PathBuf>,
}
#[cfg(feature = "cap-depict")]
impl Default for BatchImageParams {
    fn default() -> Self {
        Self {
            format: "png".into(),
            width: 300,
            height: 300,
            execution: BatchParams::default(),
            filenames: None,
            report_path: None,
        }
    }
}
#[cfg(feature = "cap-depict")]
#[derive(Debug)]
pub enum BatchImageError {
    InvalidFormat(String),
    Write(crate::DrawingWriteError),
}
#[cfg(feature = "cap-depict")]
impl std::fmt::Display for BatchImageError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidFormat(format) => write!(
                f,
                "unsupported image format '{format}', expected 'png' or 'svg'"
            ),
            Self::Write(error) => error.fmt(f),
        }
    }
}
#[cfg(feature = "cap-depict")]
impl std::error::Error for BatchImageError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Write(e) => Some(e),
            _ => None,
        }
    }
}
#[cfg(feature = "cap-depict")]
impl MoleculeBatch {
    pub fn write_images(
        &self,
        directory: &std::path::Path,
    ) -> Result<crate::BatchExportReport, BatchValidationError> {
        self.write_images_with_params(directory, &BatchImageParams::default())
    }
    pub fn write_images_with_params(
        &self,
        directory: &std::path::Path,
        options: &BatchImageParams,
    ) -> Result<crate::BatchExportReport, BatchValidationError> {
        // COSMolKit❗✔️: pinned d892ec3 properties/batch.rs::write_images_with_runtime; same indexed order/all-result collection and scalar owner calls.
        // fn write_images_with_runtime(
        //         &self,
        //         output_dir: &Path,
        //         format: &str,
        //         width: u32,
        //         height: u32,
        //         errors: BatchErrorMode,
        //         filenames: Option<&[String]>,
        //         progress: BatchProgress<'_>,
        //         n_jobs: Option<usize>,
        //     ) -> Result<BatchExportReport, BatchValidationError> {
        //         fs::create_dir_all(output_dir).map_err(|_| {
        //             BatchValidationError::unsupported("batch.write_images", "directory creation failed")
        //         })?;
        //         let paths = output_paths(output_dir, self.records.len(), format, filenames)?;
        //         let outcomes: Vec<Result<bool, BatchRecordError>> =
        //             self.run_with_parallel_jobs(n_jobs, || {
        //                 self.records
        //                     .par_iter()
        //                     .enumerate()
        //                     .map(|(index, record)| {
        //                         let out = match record {
        //                             BatchRecord::Molecule(molecule) => {
        //                                 write_one_image(molecule, &paths[index], format, width, height)
        //                                     .map(|()| true)
        //                                     .map_err(|error| {
        //                                         BatchRecordError::new(index, "batch.write_images", error)
        //                                     })
        //                             }
        //                             BatchRecord::Error(_) => Ok(false),
        //                         };
        //                         tick_progress(progress);
        //                         out
        //                     })
        //                     .collect()
        //             });
        //         let mut written = 0usize;
        //         let mut skipped = 0usize;
        //         let mut record_errors = Vec::new();
        //         for outcome in outcomes {
        //             match outcome {
        //                 Ok(true) => written += 1,
        //                 Ok(false) => skipped += 1,
        //                 Err(error) => record_errors.push(error),
        //             }
        //         }
        //         if errors.raise_on_errors() && (skipped != 0 || !record_errors.is_empty()) {
        //             let mut all_errors = record_errors.clone();
        //             all_errors.extend(record_errors_from_records(&self.records));
        //             return Err(BatchValidationError::from_record_errors(all_errors));
        //         }
        //         Ok(BatchExportReport {
        //             written,
        //             skipped,
        //             failed: record_errors.len(),
        //             errors: record_errors,
        //         })
        //     }

        let params = options.execution;
        if params.n_jobs == Some(0) {
            return Err(BatchValidationError::parameter(
                "n_jobs",
                "n_jobs must be >= 1",
            ));
        }
        if let Some(names) = &options.filenames
            && names.len() != self.len()
        {
            return Err(BatchValidationError::parameter(
                "filenames",
                "filenames length must match batch length",
            ));
        }
        std::fs::create_dir_all(directory).map_err(|e| {
            BatchValidationError::from_record_errors(vec![BatchError::with_source(
                0,
                "batch.write_images",
                e,
            )])
        })?;
        let paths = cosmolkit_batch::output_paths_with_extension(
            directory,
            self.len(),
            &options.format,
            options.filenames.as_deref(),
        )
        .map_err(|e| {
            let duplicate = e
                .record_errors
                .first()
                .is_some_and(|e| e.message == "duplicate output filename");
            let mut error = BatchError::with_source(0, "batch.output_paths", e);
            error.message = if duplicate {
                "duplicate output filename"
            } else {
                "invalid filename"
            }
            .into();
            BatchValidationError::from_record_errors(vec![error])
        })?;
        let outcomes = cosmolkit_batch::run_indexed(
            self.len(),
            params.n_jobs.or(self.n_jobs),
            params.progress_bar.or(self.progress_bar),
            "Writing molecule images",
            |index| match &self.records[index] {
                BatchRecord::Error(error) => Err(error.clone()),
                BatchRecord::Molecule(molecule) => {
                    let write = match options.format.as_str() {
                        "png" => molecule
                            .write_png(&paths[index], options.width, options.height)
                            .map_err(BatchImageError::Write),
                        "svg" => molecule
                            .write_svg(&paths[index], options.width, options.height)
                            .map_err(BatchImageError::Write),
                        other => Err(BatchImageError::InvalidFormat(other.to_owned())),
                    };
                    write.map_err(|e| BatchError::with_source(index, "batch.write_images", e))
                }
            },
        )?;
        let mut written = 0;
        let mut errors = Vec::new();
        for outcome in outcomes {
            match outcome {
                Ok(()) => written += 1,
                Err(error) => errors.push(error),
            }
        }
        // User-approved report contract: original and new errors are retained
        // once per unsuccessful record, in input order, with their own causes.
        if params.errors.unwrap_or(self.error_mode) == BatchErrorMode::Strict && !errors.is_empty()
        {
            return Err(BatchValidationError::from_record_errors(errors));
        }
        let report = crate::BatchExportReport {
            written,
            failed: errors.len(),
            errors,
        };
        if let Some(path) = &options.report_path {
            cosmolkit_batch::write_export_report(path, &report)?;
        }
        Ok(report)
    }
}

#[cfg(feature = "cap-fingerprints")]
impl MoleculeBatch {
    fn atom_pair_preflight<'a>(
        options: &'a crate::AtomPairFingerprintParams,
        operation: &'static str,
    ) -> Result<cosmolkit_fingerprints::AtomPairBatchArguments<'a>, BatchValidationError> {
        // COSMolKit❗✔️: pinned d892ec3 properties/batch.rs::atom_pair_generator_for_batch; same indexed order/all-result collection and scalar owner calls.
        // fn atom_pair_generator_for_batch(
        //     params: &crate::AtomPairFingerprintParams,
        //     operation: &'static str,
        // ) -> Result<crate::AtomPairFingerprintGenerator, BatchValidationError> {
        //     crate::AtomPairFingerprintGenerator::new(params).map_err(|error| {
        //         BatchValidationError::simple(1, operation, format!("invalid configuration: {error}"))
        //     })
        // }

        cosmolkit_fingerprints::AtomPairBatchArguments::new(&options.generator).map_err(|cause| {
            let mut error = BatchError::with_source(0, operation, cause);
            error.message = format!("invalid configuration: {}", error.message);
            BatchValidationError::from_record_errors(vec![error])
        })
    }
    pub fn fingerprint_atom_pair_list(
        &self,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        self.fingerprint_atom_pair_list_with_params(
            &crate::AtomPairFingerprintParams::default(),
            &BatchQueryParams::default(),
        )
    }
    pub fn fingerprint_atom_pair_list_with_params(
        &self,
        options: &crate::AtomPairFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        // COSMolKit❗✔️: pinned d892ec3 properties/batch.rs::atom_pair_fingerprint_list_with_runtime.
        // Cost: one immutable owner argument construction before scheduling;
        // count_bounds is copied once, O(bounds_length), and borrowed by every
        // row. Per-record chemistry and separately marked AO snapshot costs
        // remain in their existing detached owners.
        // fn atom_pair_fingerprint_list_with_runtime(
        //         &self,
        //         params: &crate::AtomPairFingerprintParams,
        //         progress: BatchProgress<'_>,
        //         n_jobs: Option<usize>,
        //     ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        //         let generator = atom_pair_generator_for_batch(params, "batch.atom_pair_fingerprint")?;
        //         self.collect_optional_values_with_options(
        //             "batch.atom_pair_fingerprint",
        //             |molecule| {
        //                 generator
        //                     .fingerprint(
        //                         molecule,
        //                         &mut crate::fingerprint::atom_pair_function_arguments(params),
        //                     )
        //                     .map_err(|error| error.to_string())
        //             },
        //             progress,
        //             n_jobs,
        //         )
        //     }
        let generator = Self::atom_pair_preflight(options, "batch.fingerprint_atom_pair")?;
        self.collect(
            "batch.fingerprint_atom_pair",
            "Computing AtomPair fingerprints",
            params,
            |m| {
                m.with_atom_pair_batch_input(options, |input, call| {
                    generator.bits(input, call, None)
                })
            },
        )
    }
    pub fn fingerprint_atom_pair_sparse_count_list(
        &self,
    ) -> Result<Vec<Option<crate::SparseCountFingerprint>>, BatchValidationError> {
        self.fingerprint_atom_pair_sparse_count_list_with_params(
            &crate::AtomPairFingerprintParams::default(),
            &BatchQueryParams::default(),
        )
    }
    pub fn fingerprint_atom_pair_sparse_count_list_with_params(
        &self,
        options: &crate::AtomPairFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Vec<Option<crate::SparseCountFingerprint>>, BatchValidationError> {
        // COSMolKit❗✔️: pinned d892ec3 properties/batch.rs::atom_pair_sparse_count_fingerprint_list_with_runtime.
        // Cost: one immutable owner argument construction before scheduling;
        // count_bounds is copied once, O(bounds_length), and borrowed by every
        // row. Per-record chemistry and separately marked AO snapshot costs
        // remain in their existing detached owners.
        // fn atom_pair_sparse_count_fingerprint_list_with_runtime(
        //         &self,
        //         params: &crate::AtomPairFingerprintParams,
        //         progress: BatchProgress<'_>,
        //         n_jobs: Option<usize>,
        //     ) -> Result<Vec<Option<crate::SparseCountFingerprint>>, BatchValidationError> {
        //         let generator =
        //             atom_pair_generator_for_batch(params, "batch.atom_pair_sparse_count_fingerprint")?;
        //         self.collect_optional_values_with_options(
        //             "batch.atom_pair_sparse_count_fingerprint",
        //             |molecule| {
        //                 generator
        //                     .sparse_count_fingerprint(
        //                         molecule,
        //                         &mut crate::fingerprint::atom_pair_function_arguments(params),
        //                     )
        //                     .map_err(|error| error.to_string())
        //             },
        //             progress,
        //             n_jobs,
        //         )
        //     }
        let generator =
            Self::atom_pair_preflight(options, "batch.fingerprint_atom_pair_sparse_count")?;
        self.collect(
            "batch.fingerprint_atom_pair_sparse_count",
            "Computing AtomPair fingerprints",
            params,
            |m| {
                m.with_atom_pair_batch_input(options, |input, call| {
                    generator.sparse_count(input, call, None)
                })
            },
        )
    }
    pub fn fingerprint_atom_pair_count_list(
        &self,
    ) -> Result<Vec<Option<crate::SparseCountFingerprint32>>, BatchValidationError> {
        self.fingerprint_atom_pair_count_list_with_params(
            &crate::AtomPairFingerprintParams::default(),
            &BatchQueryParams::default(),
        )
    }
    pub fn fingerprint_atom_pair_count_list_with_params(
        &self,
        options: &crate::AtomPairFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Vec<Option<crate::SparseCountFingerprint32>>, BatchValidationError> {
        // COSMolKit❗✔️: pinned d892ec3 properties/batch.rs::atom_pair_count_fingerprint_list_with_runtime.
        // Cost: one immutable owner argument construction before scheduling;
        // count_bounds is copied once, O(bounds_length), and borrowed by every
        // row. Per-record chemistry and separately marked AO snapshot costs
        // remain in their existing detached owners.
        // fn atom_pair_count_fingerprint_list_with_runtime(
        //         &self,
        //         params: &crate::AtomPairFingerprintParams,
        //         progress: BatchProgress<'_>,
        //         n_jobs: Option<usize>,
        //     ) -> Result<Vec<Option<crate::SparseCountFingerprint>>, BatchValidationError> {
        //         let generator = atom_pair_generator_for_batch(params, "batch.atom_pair_count_fingerprint")?;
        //         self.collect_optional_values_with_options(
        //             "batch.atom_pair_count_fingerprint",
        //             |molecule| {
        //                 generator
        //                     .count_fingerprint(
        //                         molecule,
        //                         &mut crate::fingerprint::atom_pair_function_arguments(params),
        //                     )
        //                     .map_err(|error| error.to_string())
        //             },
        //             progress,
        //             n_jobs,
        //         )
        //     }
        let generator = Self::atom_pair_preflight(options, "batch.fingerprint_atom_pair_count")?;
        self.collect(
            "batch.fingerprint_atom_pair_count",
            "Computing AtomPair fingerprints",
            params,
            |m| {
                m.with_atom_pair_batch_input(options, |input, call| {
                    generator.count(input, call, None)
                })
            },
        )
    }
    pub fn fingerprint_atom_pair_sparse_bits_list(
        &self,
    ) -> Result<Vec<Option<crate::SparseBitFingerprint>>, BatchValidationError> {
        self.fingerprint_atom_pair_sparse_bits_list_with_params(
            &crate::AtomPairFingerprintParams::default(),
            &BatchQueryParams::default(),
        )
    }
    pub fn fingerprint_atom_pair_sparse_bits_list_with_params(
        &self,
        options: &crate::AtomPairFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Vec<Option<crate::SparseBitFingerprint>>, BatchValidationError> {
        // COSMolKit❗✔️: pinned d892ec3 properties/batch.rs::atom_pair_sparse_bit_fingerprint_list_with_runtime.
        // Cost: one immutable owner argument construction before scheduling;
        // count_bounds is copied once, O(bounds_length), and borrowed by every
        // row. Per-record chemistry and separately marked AO snapshot costs
        // remain in their existing detached owners.
        // fn atom_pair_sparse_bit_fingerprint_list_with_runtime(
        //         &self,
        //         params: &crate::AtomPairFingerprintParams,
        //         progress: BatchProgress<'_>,
        //         n_jobs: Option<usize>,
        //     ) -> Result<Vec<Option<crate::SparseBitFingerprint>>, BatchValidationError> {
        //         let generator =
        //             atom_pair_generator_for_batch(params, "batch.atom_pair_sparse_bit_fingerprint")?;
        //         self.collect_optional_values_with_options(
        //             "batch.atom_pair_sparse_bit_fingerprint",
        //             |molecule| {
        //                 generator
        //                     .sparse_bit_fingerprint(
        //                         molecule,
        //                         &mut crate::fingerprint::atom_pair_function_arguments(params),
        //                     )
        //                     .map_err(|error| error.to_string())
        //             },
        //             progress,
        //             n_jobs,
        //         )
        //     }
        let generator =
            Self::atom_pair_preflight(options, "batch.atom_pair_sparse_bit_fingerprint")?;
        self.collect(
            "batch.atom_pair_sparse_bit_fingerprint",
            "Computing AtomPair fingerprints",
            params,
            |m| {
                m.with_atom_pair_batch_input(options, |input, call| {
                    generator.sparse_bits(input, call, None)
                })
            },
        )
    }
    pub fn fingerprint_layered_list(
        &self,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        self.fingerprint_layered_list_with_params(
            &crate::LayeredFingerprintParams::default(),
            &BatchQueryParams::default(),
        )
    }
    pub fn fingerprint_layered_list_with_params(
        &self,
        options: &crate::LayeredFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        // COSMolKit❗✔️: pinned d892ec3 properties/batch.rs::layered_fingerprint_list_with_runtime; same indexed order/all-result collection and scalar owner calls.
        // fn layered_fingerprint_list_with_runtime(
        //         &self,
        //         params: &crate::LayeredFingerprintParams,
        //         progress: BatchProgress<'_>,
        //         n_jobs: Option<usize>,
        //     ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        //         self.collect_optional_values_with_options(
        //             "batch.layered_fingerprint",
        //             |molecule| {
        //                 molecule
        //                     .layered_fingerprint(params)
        //                     .map_err(|error| error.to_string())
        //             },
        //             progress,
        //             n_jobs,
        //         )
        //     }
        self.collect(
            "batch.fingerprint_layered",
            "Computing Layered fingerprints",
            params,
            |m| m.fingerprint_layered_with_params(options),
        )
    }
    pub fn fingerprint_layered_with_output_list(
        &self,
    ) -> Result<Vec<Option<crate::LayeredFingerprintResult>>, BatchValidationError> {
        self.fingerprint_layered_with_output_list_with_params(
            &crate::LayeredFingerprintParams::default(),
            &BatchQueryParams::default(),
        )
    }
    pub fn fingerprint_layered_with_output_list_with_params(
        &self,
        options: &crate::LayeredFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Vec<Option<crate::LayeredFingerprintResult>>, BatchValidationError> {
        // COSMolKit❗✔️: pinned d892ec3 properties/batch.rs::layered_fingerprint_with_output_list_with_runtime; same indexed order/all-result collection and scalar owner calls.
        // fn layered_fingerprint_with_output_list_with_runtime(
        //         &self,
        //         params: &crate::LayeredFingerprintParams,
        //         progress: BatchProgress<'_>,
        //         n_jobs: Option<usize>,
        //     ) -> Result<Vec<Option<crate::LayeredFingerprintResult>>, BatchValidationError> {
        //         self.collect_optional_values_with_options(
        //             "batch.layered_fingerprint_with_output",
        //             |molecule| {
        //                 molecule
        //                     .layered_fingerprint_with_output(params)
        //                     .map_err(|error| error.to_string())
        //             },
        //             progress,
        //             n_jobs,
        //         )
        //     }
        self.collect(
            "batch.fingerprint_layered_with_output",
            "Computing Layered fingerprints",
            params,
            |m| m.fingerprint_layered_with_output_with_params(options),
        )
    }
    pub fn fingerprint_pattern_list(
        &self,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        self.fingerprint_pattern_list_with_params(
            &crate::PatternFingerprintParams::default(),
            &BatchQueryParams::default(),
        )
    }
    pub fn fingerprint_pattern_list_with_params(
        &self,
        options: &crate::PatternFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        // COSMolKit❗✔️: pinned d892ec3 properties/batch.rs::pattern_fingerprint_list_with_runtime; same indexed order/all-result collection and scalar owner calls.
        // fn pattern_fingerprint_list_with_runtime(
        //         &self,
        //         params: &crate::PatternFingerprintParams,
        //         progress: BatchProgress<'_>,
        //         n_jobs: Option<usize>,
        //     ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        //         self.collect_optional_values_with_options(
        //             "batch.pattern_fingerprint",
        //             |molecule| {
        //                 molecule
        //                     .pattern_fingerprint(params)
        //                     .map_err(|error| error.to_string())
        //             },
        //             progress,
        //             n_jobs,
        //         )
        //     }
        self.collect(
            "batch.fingerprint_pattern",
            "Computing Pattern fingerprints",
            params,
            |m| m.fingerprint_pattern_with_params(options),
        )
    }
    pub fn fingerprint_morgan_list(
        &self,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        let options = crate::MorganFingerprintParams::default();
        self.fingerprint_morgan_list_with_params(&options, &BatchQueryParams::default())
    }
    pub fn fingerprint_morgan_list_with_params(
        &self,
        options: &crate::MorganFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        // COSMolKit❗✔️: pinned d892ec3 properties/batch.rs::morgan_fingerprint_list_with_runtime; same indexed order/all-result collection and scalar owner calls.
        // fn morgan_fingerprint_list_with_runtime(
        //         &self,
        //         params: &crate::MorganFingerprintParams,
        //         progress: BatchProgress<'_>,
        //         n_jobs: Option<usize>,
        //     ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        //         self.collect_optional_values_with_options(
        //             "batch.morgan_fingerprint",
        //             |molecule| {
        //                 molecule
        //                     .morgan_fingerprint(params)
        //                     .map_err(|error| error.to_string())
        //             },
        //             progress,
        //             n_jobs,
        //         )
        //     }
        let (atom, call) = Self::morgan_batch_arguments(options);
        self.fingerprint_morgan_list_with_generator_params(
            &options.generator,
            Some(&atom),
            None,
            &call,
            params,
        )
    }
}

#[cfg(feature = "cap-fingerprints")]
impl MoleculeBatch {
    fn additional_output(collect: bool) -> Option<crate::FingerprintAdditionalOutput> {
        // COSMolKit❗✔️: pinned d892ec3 fingerprint/atom_pair.rs::atom_pair_function_arguments; reuse canonical owner, preserve allocation mask and per-record preflight.
        // pub(crate) fn atom_pair_function_arguments(
        //     params: &AtomPairFingerprintParams,
        // ) -> FingerprintFuncArguments {
        //     let mut arguments = FingerprintFuncArguments {
        //         from_atoms: params.from_atoms.clone(),
        //         ignore_atoms: params.ignore_atoms.clone(),
        //         conf_id: params.conformer_id,
        //         custom_atom_invariants: params.custom_atom_invariants.clone(),
        //         ..Default::default()
        //     };
        //     if params.collect_additional_output {
        //         let mut additional_output = AdditionalOutput::new();
        //         additional_output.allocate_atom_counts();
        //         additional_output.allocate_atom_to_bits();
        //         additional_output.allocate_bit_info_map();
        //         additional_output.allocate_atoms_per_bit();
        //         arguments.additional_output = Some(additional_output);
        //     }
        //     arguments
        // }

        collect.then(|| {
            let mut output = crate::FingerprintAdditionalOutput::new();
            output.allocate_atom_counts();
            output.allocate_atom_to_bits();
            output.allocate_bit_info_map();
            output.allocate_atoms_per_bit();
            output
        })
    }
    pub fn fingerprint_atom_pair_with_output_list(
        &self,
    ) -> Result<Vec<Option<BatchFingerprintOutput>>, BatchValidationError> {
        self.fingerprint_atom_pair_with_output_list_with_params(
            &crate::AtomPairFingerprintParams::default(),
            true,
            &BatchQueryParams::default(),
        )
    }
    pub fn fingerprint_atom_pair_with_output_list_with_params(
        &self,
        options: &crate::AtomPairFingerprintParams,
        collect_additional_output: bool,
        params: &BatchQueryParams,
    ) -> Result<Vec<Option<BatchFingerprintOutput>>, BatchValidationError> {
        // COSMolKit❗✔️: pinned d892ec3 properties/batch.rs::atom_pair_fingerprint_with_output_list_with_runtime.
        // Cost: one immutable owner argument construction before scheduling;
        // count_bounds is copied once, O(bounds_length), and borrowed by every
        // row. Per-record chemistry and separately marked AO snapshot costs
        // remain in their existing detached owners.
        // fn atom_pair_fingerprint_with_output_list_with_runtime(
        //         &self,
        //         params: &crate::AtomPairFingerprintParams,
        //         progress: BatchProgress<'_>,
        //         n_jobs: Option<usize>,
        //     ) -> Result<Vec<Option<crate::AtomPairFingerprintOutput>>, BatchValidationError> {
        //         let generator =
        //             atom_pair_generator_for_batch(params, "batch.atom_pair_fingerprint_with_output")?;
        //         self.collect_optional_values_with_options(
        //             "batch.atom_pair_fingerprint_with_output",
        //             |molecule| {
        //                 let mut arguments = crate::fingerprint::atom_pair_function_arguments(params);
        //                 let fingerprint = generator
        //                     .fingerprint(molecule, &mut arguments)
        //                     .map_err(|error| error.to_string())?;
        //                 Ok(crate::AtomPairFingerprintOutput {
        //                     fingerprint,
        //                     additional_output: arguments.additional_output,
        //                 })
        //             },
        //             progress,
        //             n_jobs,
        //         )
        //     }

        let generator =
            Self::atom_pair_preflight(options, "batch.fingerprint_atom_pair_with_output")?;
        self.collect(
            "batch.fingerprint_atom_pair_with_output",
            "Computing AtomPair fingerprints",
            params,
            |m| {
                let mut output = Self::additional_output(collect_additional_output);
                let fingerprint = m.with_atom_pair_batch_input(options, |input, call| {
                    generator.bits(input, call, output.as_mut())
                })?;
                Ok::<_, crate::AtomPairReadError>(cosmolkit_fingerprints::batch_fingerprint_output(
                    fingerprint,
                    output.as_ref(),
                    "AtomPair",
                ))
            },
        )
    }
    pub fn fingerprint_morgan_with_output_list(
        &self,
    ) -> Result<Vec<Option<BatchFingerprintOutput>>, BatchValidationError> {
        let options = crate::MorganFingerprintParams::default();
        self.fingerprint_morgan_with_output_list_with_params(
            &options,
            true,
            &BatchQueryParams::default(),
        )
    }
    pub fn fingerprint_morgan_with_output_list_with_params(
        &self,
        options: &crate::MorganFingerprintParams,
        collect_additional_output: bool,
        params: &BatchQueryParams,
    ) -> Result<Vec<Option<BatchFingerprintOutput>>, BatchValidationError> {
        // COSMolKit❗✔️: pinned d892ec3 properties/batch.rs::morgan_fingerprint_with_output_list_with_runtime; same indexed order/all-result collection and scalar owner calls.
        // fn morgan_fingerprint_with_output_list_with_runtime(
        //         &self,
        //         params: &crate::MorganFingerprintParams,
        //         progress: BatchProgress<'_>,
        //         n_jobs: Option<usize>,
        //     ) -> Result<Vec<Option<crate::MorganFingerprintOutput>>, BatchValidationError> {
        //         self.collect_optional_values_with_options(
        //             "batch.morgan_fingerprint_with_output",
        //             |molecule| {
        //                 molecule
        //                     .morgan_fingerprint_with_output(params)
        //                     .map_err(|error| error.to_string())
        //             },
        //             progress,
        //             n_jobs,
        //         )
        //     }

        let (atom, call) = Self::morgan_batch_arguments(options);
        self.fingerprint_morgan_with_output_list_with_generator_params(
            &options.generator,
            Some(&atom),
            None,
            &call,
            collect_additional_output,
            params,
        )
    }
    fn morgan_batch_arguments(
        options: &crate::MorganFingerprintParams,
    ) -> (
        crate::MorganAtomInvariantsGenerator,
        crate::MorganCallParams,
    ) {
        // Transport only: source provider selection and owned per-call arguments
        // from the verbatim morgan_fingerprint_with_output anchor below.
        let atom = match &options.invariants {
            crate::MorganInvariants::Connectivity => {
                crate::MorganAtomInvariantsGenerator::connectivity(
                    options.generator.include_ring_membership,
                )
            }
            crate::MorganInvariants::Features => {
                crate::MorganAtomInvariantsGenerator::features(None)
            }
            crate::MorganInvariants::FeaturePatterns(patterns) => {
                crate::MorganAtomInvariantsGenerator::features(Some(patterns.clone()))
            }
        };
        let call = crate::MorganCallParams::new(
            options.from_atoms.clone(),
            options.ignore_atoms.clone(),
            options.custom_atom_invariants.clone(),
            options.custom_bond_invariants.clone(),
            options.conformer_id,
        );
        (atom, call)
    }
    fn morgan_batch_generator(
        options: &crate::MorganParams,
        atom_invariants: Option<&crate::MorganAtomInvariantsGenerator>,
        bond_invariants: Option<&crate::MorganBondInvariantsGenerator>,
    ) -> Result<crate::MorganFingerprintGenerator, crate::MorganReadError> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855
        // properties/fingerprint.rs::morgan_fingerprint_with_output, full defining wrapper.
        // pub fn morgan_fingerprint_with_output(
        //     molecule: &Molecule,
        //     params: &MorganFingerprintParams,
        // ) -> Result<MorganFingerprintOutput, FingerprintError> {
        //     validate_morgan_params(params)?;
        //     if params.n_bits == 0 {
        //         return Err(FingerprintError::EmptyFingerprint);
        //     }
        //     let atom_invariants_generator = match params.atom_invariants_generator.clone() {
        //         MorganAtomInvariantsGenerator::Connectivity {
        //             include_ring_membership,
        //         } => Some(MorganAtomInvariantsGenerator::Connectivity {
        //             include_ring_membership,
        //         }),
        //         MorganAtomInvariantsGenerator::Feature => Some(MorganAtomInvariantsGenerator::Feature),
        //     };
        //     let bond_invariants_generator =
        //         Some(params.bond_invariants_generator.clone().unwrap_or_else(|| {
        //             MorganBondInvariantsGenerator {
        //                 use_bond_types: params.use_bond_types,
        //                 use_chirality: params.use_chirality,
        //             }
        //         }));
        //     let mut generator = getMorganGeneratorWithParams(
        //         params.radius,
        //         params.count_simulation,
        //         params.use_chirality,
        //         params.use_bond_types,
        //         params.only_nonzero_invariants,
        //         params.include_redundant_environments,
        //         atom_invariants_generator,
        //         bond_invariants_generator,
        //         params.n_bits as u32,
        //         params.count_bounds.clone(),
        //         true,
        //         true,
        //     )?;
        //     // RDKit✔️✔️: fpgen->getOptions()->d_numBitsPerFeature = nBitsPerHash;
        //     generator
        //         .fingerprint_arguments
        //         .fingerprint_arguments
        //         .d_num_bits_per_feature = params.num_bits_per_feature;
        //
        //     let mut args = FingerprintFuncArguments {
        //         from_atoms: params.from_atoms.clone(),
        //         ignore_atoms: params.ignore_atoms.clone(),
        //         custom_atom_invariants: params.custom_atom_invariants.clone(),
        //         custom_bond_invariants: params.custom_bond_invariants.clone(),
        //         ..Default::default()
        //     };
        //     if params.collect_additional_output {
        //         let mut additional_output = AdditionalOutput::new();
        //         additional_output.allocate_atom_counts();
        //         additional_output.allocate_atom_to_bits();
        //         additional_output.allocate_bit_info_map();
        //         additional_output.allocate_atoms_per_bit();
        //         args.additional_output = Some(additional_output);
        //     }
        //
        //     let fingerprint = generator.getFingerprint(molecule, &mut args)?;
        //     let additional_output = args
        //         .additional_output
        //         .map(morgan_additional_output_from_rdkit_output);
        //
        //     Ok(MorganFingerprintOutput {
        //         fingerprint,
        //         additional_output,
        //     })
        // }
        // Source constructs generic arguments with one bit, then directly assigns
        // the requested value (including zero). Reuse sole canonical construction
        // and its existing live setting; no fingerprint chemistry is implemented here.
        // Complexity: same per-valid-record construction and parameter clone as source.
        cosmolkit_fingerprints::validate_morgan_batch_params(options)
            .map_err(crate::MorganReadError::Generator)?;
        let mut construction = options.clone();
        construction.bits_per_feature = 1;
        let generator = crate::MorganFingerprintGenerator::new(
            Some(&construction),
            atom_invariants,
            bond_invariants,
        )?;
        generator
            .settings()
            .set_bits_per_feature(options.bits_per_feature)?;
        Ok(generator)
    }
    pub fn fingerprint_morgan_list_with_generator_params(
        &self,
        options: &crate::MorganParams,
        atom_invariants: Option<&crate::MorganAtomInvariantsGenerator>,
        bond_invariants: Option<&crate::MorganBondInvariantsGenerator>,
        call: &crate::MorganCallParams,
        params: &BatchQueryParams,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        // COSMolKit❗✔️: pinned d892ec3 fingerprint.rs::morgan_fingerprint_with_output; reuse canonical owner, preserve allocation mask and per-record preflight.
        // pub fn morgan_fingerprint_with_output(
        //     molecule: &Molecule,
        //     params: &MorganFingerprintParams,
        // ) -> Result<MorganFingerprintOutput, FingerprintError> {
        //     validate_morgan_params(params)?;
        //     if params.n_bits == 0 {
        //         return Err(FingerprintError::EmptyFingerprint);
        //     }
        //     let atom_invariants_generator = match params.atom_invariants_generator.clone() {
        //         MorganAtomInvariantsGenerator::Connectivity {
        //             include_ring_membership,
        //         } => Some(MorganAtomInvariantsGenerator::Connectivity {
        //             include_ring_membership,
        //         }),
        //         MorganAtomInvariantsGenerator::Feature => Some(MorganAtomInvariantsGenerator::Feature),
        //     };
        //     let bond_invariants_generator =
        //         Some(params.bond_invariants_generator.clone().unwrap_or_else(|| {
        //             MorganBondInvariantsGenerator {
        //                 use_bond_types: params.use_bond_types,
        //                 use_chirality: params.use_chirality,
        //             }
        //         }));
        //     let mut generator = getMorganGeneratorWithParams(
        //         params.radius,
        //         params.count_simulation,
        //         params.use_chirality,
        //         params.use_bond_types,
        //         params.only_nonzero_invariants,
        //         params.include_redundant_environments,
        //         atom_invariants_generator,
        //         bond_invariants_generator,
        //         params.n_bits as u32,
        //         params.count_bounds.clone(),
        //         true,
        //         true,
        //     )?;
        //     // RDKit✔️✔️: fpgen->getOptions()->d_numBitsPerFeature = nBitsPerHash;
        //     generator
        //         .fingerprint_arguments
        //         .fingerprint_arguments
        //         .d_num_bits_per_feature = params.num_bits_per_feature;
        //
        //     let mut args = FingerprintFuncArguments {
        //         from_atoms: params.from_atoms.clone(),
        //         ignore_atoms: params.ignore_atoms.clone(),
        //         custom_atom_invariants: params.custom_atom_invariants.clone(),
        //         custom_bond_invariants: params.custom_bond_invariants.clone(),
        //         ..Default::default()
        //     };
        //     if params.collect_additional_output {
        //         let mut additional_output = AdditionalOutput::new();
        //         additional_output.allocate_atom_counts();
        //         additional_output.allocate_atom_to_bits();
        //         additional_output.allocate_bit_info_map();
        //         additional_output.allocate_atoms_per_bit();
        //         args.additional_output = Some(additional_output);
        //     }
        //
        //     let fingerprint = generator.getFingerprint(molecule, &mut args)?;
        //     let additional_output = args
        //         .additional_output
        //         .map(morgan_additional_output_from_rdkit_output);
        //
        //     Ok(MorganFingerprintOutput {
        //         fingerprint,
        //         additional_output,
        //     })
        // }
        self.collect(
            "batch.fingerprint_morgan",
            "Computing Morgan fingerprints",
            params,
            |m| {
                let generator =
                    Self::morgan_batch_generator(options, atom_invariants, bond_invariants)?;
                m.fingerprint_morgan_with_generator(&generator, Some(call), None)
            },
        )
    }
    pub fn fingerprint_morgan_with_output_list_with_generator_params(
        &self,
        options: &crate::MorganParams,
        atom_invariants: Option<&crate::MorganAtomInvariantsGenerator>,
        bond_invariants: Option<&crate::MorganBondInvariantsGenerator>,
        call: &crate::MorganCallParams,
        collect_additional_output: bool,
        params: &BatchQueryParams,
    ) -> Result<Vec<Option<BatchFingerprintOutput>>, BatchValidationError> {
        // COSMolKit❗✔️: pinned d892ec3 fingerprint.rs::morgan_fingerprint_with_output; reuse canonical owner, preserve allocation mask and per-record preflight.
        // pub fn morgan_fingerprint_with_output(
        //     molecule: &Molecule,
        //     params: &MorganFingerprintParams,
        // ) -> Result<MorganFingerprintOutput, FingerprintError> {
        //     validate_morgan_params(params)?;
        //     if params.n_bits == 0 {
        //         return Err(FingerprintError::EmptyFingerprint);
        //     }
        //     let atom_invariants_generator = match params.atom_invariants_generator.clone() {
        //         MorganAtomInvariantsGenerator::Connectivity {
        //             include_ring_membership,
        //         } => Some(MorganAtomInvariantsGenerator::Connectivity {
        //             include_ring_membership,
        //         }),
        //         MorganAtomInvariantsGenerator::Feature => Some(MorganAtomInvariantsGenerator::Feature),
        //     };
        //     let bond_invariants_generator =
        //         Some(params.bond_invariants_generator.clone().unwrap_or_else(|| {
        //             MorganBondInvariantsGenerator {
        //                 use_bond_types: params.use_bond_types,
        //                 use_chirality: params.use_chirality,
        //             }
        //         }));
        //     let mut generator = getMorganGeneratorWithParams(
        //         params.radius,
        //         params.count_simulation,
        //         params.use_chirality,
        //         params.use_bond_types,
        //         params.only_nonzero_invariants,
        //         params.include_redundant_environments,
        //         atom_invariants_generator,
        //         bond_invariants_generator,
        //         params.n_bits as u32,
        //         params.count_bounds.clone(),
        //         true,
        //         true,
        //     )?;
        //     // RDKit✔️✔️: fpgen->getOptions()->d_numBitsPerFeature = nBitsPerHash;
        //     generator
        //         .fingerprint_arguments
        //         .fingerprint_arguments
        //         .d_num_bits_per_feature = params.num_bits_per_feature;
        //
        //     let mut args = FingerprintFuncArguments {
        //         from_atoms: params.from_atoms.clone(),
        //         ignore_atoms: params.ignore_atoms.clone(),
        //         custom_atom_invariants: params.custom_atom_invariants.clone(),
        //         custom_bond_invariants: params.custom_bond_invariants.clone(),
        //         ..Default::default()
        //     };
        //     if params.collect_additional_output {
        //         let mut additional_output = AdditionalOutput::new();
        //         additional_output.allocate_atom_counts();
        //         additional_output.allocate_atom_to_bits();
        //         additional_output.allocate_bit_info_map();
        //         additional_output.allocate_atoms_per_bit();
        //         args.additional_output = Some(additional_output);
        //     }
        //
        //     let fingerprint = generator.getFingerprint(molecule, &mut args)?;
        //     let additional_output = args
        //         .additional_output
        //         .map(morgan_additional_output_from_rdkit_output);
        //
        //     Ok(MorganFingerprintOutput {
        //         fingerprint,
        //         additional_output,
        //     })
        // }
        self.collect(
            "batch.fingerprint_morgan_with_output",
            "Computing Morgan fingerprints",
            params,
            |m| {
                let generator =
                    Self::morgan_batch_generator(options, atom_invariants, bond_invariants)?;
                let mut output = Self::additional_output(collect_additional_output);
                let fingerprint =
                    m.fingerprint_morgan_with_generator(&generator, Some(call), output.as_mut())?;
                Ok::<_, crate::MorganReadError>(cosmolkit_fingerprints::batch_fingerprint_output(
                    fingerprint,
                    output.as_ref(),
                    "Morgan",
                ))
            },
        )
    }
}
