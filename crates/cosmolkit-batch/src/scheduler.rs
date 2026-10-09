//! Ordered scheduling over detached indexed work; no molecule or runtime capability.
use crate::{BatchProgressBar, BatchRecordError, BatchValidationError};
use rayon::prelude::*;
pub(crate) fn parallel<R: Send>(
    n_jobs: Option<usize>,
    f: impl FnOnce() -> R + Send,
) -> Result<R, BatchValidationError> {
    // COSMolKit❗✔️: d892ec3 properties/batch.rs source; indexed Rayon scheduling and one selected pool.
    // fn run_with_parallel_jobs_option<R: Send>(
    //     n_jobs: Option<usize>,
    //     f: impl FnOnce() -> R + Send,
    // ) -> R {
    //     match n_jobs.map(|value| value.max(1)) {
    //         Some(n_jobs) => rayon::ThreadPoolBuilder::new()
    //             .num_threads(n_jobs)
    //             .build()
    //             .expect("batch rayon thread pool build must succeed")
    //             .install(f),
    //         None => f(),
    //     }
    // }
    // Python source rejects zero; pool failures retain their real structural cause.

    if n_jobs == Some(0) {
        return Err(BatchValidationError::parameter(
            "n_jobs",
            "n_jobs must be >= 1",
        ));
    }
    // wasm32 has no native worker threads: retain Rayon's serial fallback for
    // the default/one-worker case instead of trying to spawn a thread pool.
    #[cfg(target_arch = "wasm32")]
    if n_jobs.is_none() || n_jobs == Some(1) {
        return Ok(f());
    }
    // Public scheduling policy: an omitted worker count uses one thread. The
    // facade resolves stored batch overrides before calling this boundary.
    rayon::ThreadPoolBuilder::new()
        .num_threads(n_jobs.unwrap_or(1))
        .build()
        .map_err(|e| {
            BatchValidationError::from_record_errors(vec![BatchRecordError::with_source(
                0, "n_jobs", e,
            )])
        })
        .map(|pool| pool.install(f))
}
/// Complete every indexed item in the selected pool, then return input order.
/// The callback belongs to the caller; this owner has no chemistry authority.
pub fn run_indexed<T: Send>(
    total: usize,
    n_jobs: Option<usize>,
    progress_bar: Option<bool>,
    message: &'static str,
    work: impl Fn(usize) -> T + Send + Sync,
) -> Result<Vec<T>, BatchValidationError> {
    // COSMolKit❗✔️: d892ec3 properties/batch.rs::transform_with_options:
    //             self.records
    //                 .par_iter()
    //                 .enumerate()
    //                 .map(|(index, record)| {
    //                     tick_progress(progress);
    //                     out
    //                 })
    //                 .collect()
    // An indexed range preserves the same ordered collection and one work/tick
    // per row, including failures. One result allocation; no molecule cloning.
    let progress = progress_bar
        .unwrap_or(false)
        .then(|| BatchProgressBar::new(total, message));
    let result = parallel(n_jobs, || {
        (0..total)
            .into_par_iter()
            .map(|index| {
                let value = work(index);
                if let Some(progress) = &progress {
                    progress.inc(1);
                }
                value
            })
            .collect()
    });
    if let Some(progress) = progress {
        progress.finish();
    }
    result
}
