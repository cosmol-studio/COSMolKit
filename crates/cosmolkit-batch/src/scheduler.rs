//! Ordered scheduling over detached indexed work; no molecule or runtime capability.
use crate::{BatchProgressBar, BatchRecordError, BatchValidationError};
use rayon::prelude::*;
use std::sync::{Arc, Mutex};

// Keep one idle pool, not one pool for every worker count ever requested.
// In-flight calls hold their own Arc when a later call replaces the cache.
static POOL: Mutex<Option<Arc<rayon::ThreadPool>>> = Mutex::new(None);

fn pool_failure(source: impl std::error::Error + Send + Sync + 'static) -> BatchValidationError {
    BatchValidationError::from_record_errors(vec![BatchRecordError::with_source(
        0, "n_jobs", source,
    )])
}

fn selected_pool(
    cache: &Mutex<Option<Arc<rayon::ThreadPool>>>,
    workers: usize,
) -> Result<Arc<rayon::ThreadPool>, BatchValidationError> {
    let mut cached = cache
        .lock()
        .map_err(|_| pool_failure(std::io::Error::other("batch thread-pool cache poisoned")))?;
    if let Some(pool) = cached
        .as_ref()
        .filter(|pool| pool.current_num_threads() == workers)
    {
        return Ok(Arc::clone(pool));
    }
    let pool = Arc::new(
        rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .map_err(pool_failure)?,
    );
    let previous = cached.replace(Arc::clone(&pool));
    drop(cached);
    // Pool teardown must not hold the cache lock: existing work can reenter.
    drop(previous);
    Ok(pool)
}

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
    let workers = n_jobs.unwrap_or(1);
    // A caller already inside an equally sized Rayon pool can reuse it, including
    // nested batches. Do not inherit a different width or change the default of one.
    if rayon::current_thread_index().is_some() && rayon::current_num_threads() == workers {
        return Ok(f());
    }
    selected_pool(&POOL, workers).map(|pool| pool.install(f))
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_worker_count_reuses_pool_and_replacement_is_bounded() {
        let cache = Mutex::new(None);
        let first = selected_pool(&cache, 2).unwrap();
        assert_eq!(first.install(rayon::current_num_threads), 2);
        let again = selected_pool(&cache, 2).unwrap();
        assert!(Arc::ptr_eq(&first, &again));
        let old = Arc::downgrade(&first);
        let next = selected_pool(&cache, 3).unwrap();
        assert_eq!(next.install(rayon::current_num_threads), 3);
        assert!(!Arc::ptr_eq(&first, &next));
        assert_eq!(first.install(|| 7), 7);
        drop(first);
        drop(again);
        assert!(old.upgrade().is_none());
    }

    #[test]
    fn pool_install_can_reenter_after_replacement() {
        let cache = Mutex::new(None);
        let pool = selected_pool(&cache, 2).unwrap();
        let value = pool.install(|| {
            selected_pool(&cache, 3)
                .unwrap()
                .install(|| selected_pool(&cache, 2).unwrap().install(|| 42))
        });
        assert_eq!(value, 42);
    }

    #[test]
    fn public_scheduling_preserves_defaults_order_and_item_errors() {
        assert_eq!(parallel(None, rayon::current_num_threads).unwrap(), 1);
        assert_eq!(parallel(Some(3), rayon::current_num_threads).unwrap(), 3);
        assert!(parallel(Some(0), || ()).is_err());
        let values = run_indexed(24, Some(3), None, "test", |index| {
            if index == 7 {
                Err(index)
            } else {
                Ok(index * 2)
            }
        })
        .unwrap();
        assert_eq!(values.len(), 24);
        for (index, value) in values.into_iter().enumerate() {
            assert_eq!(value, if index == 7 { Err(7) } else { Ok(index * 2) });
        }
    }

    #[test]
    fn nested_batch_reuses_equal_width_without_inheriting_other_widths() {
        let outer = rayon::ThreadPoolBuilder::new()
            .num_threads(3)
            .build()
            .unwrap();
        outer.install(|| {
            let current = std::thread::current().id();
            assert_eq!(
                parallel(Some(3), || std::thread::current().id()).unwrap(),
                current
            );
            assert_eq!(parallel(None, rayon::current_num_threads).unwrap(), 1);
            assert_eq!(parallel(Some(4), rayon::current_num_threads).unwrap(), 4);
            assert!(parallel(Some(0), || ()).is_err());
        });
    }
}
