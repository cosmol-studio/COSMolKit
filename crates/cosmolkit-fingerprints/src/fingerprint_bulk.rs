//! Shared source strided fingerprint worker scheduling, with ordered None slots.
use std::{fmt, num::NonZeroUsize};
#[derive(Debug)]
pub enum FingerprintWorkerError {
    ThreadCount(cosmolkit_core::ThreadCountError),
    ThreadSpawn(std::io::Error),
    Panic,
    Protocol,
}
impl fmt::Display for FingerprintWorkerError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ThreadCount(e) => e.fmt(f),
            Self::ThreadSpawn(e) => write!(f, "fingerprint worker creation failed: {e}"),
            Self::Panic => f.write_str("fingerprint worker panicked"),
            Self::Protocol => {
                f.write_str("fingerprint worker returned an incomplete result sequence")
            }
        }
    }
}
impl std::error::Error for FingerprintWorkerError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::ThreadCount(e) => Some(e),
            Self::ThreadSpawn(e) => Some(e),
            Self::Panic | Self::Protocol => None,
        }
    }
}
pub(crate) fn workers(request: i32) -> Result<NonZeroUsize, FingerprintWorkerError> {
    cosmolkit_core::rdkit_thread_count(request)
        .map(|n| NonZeroUsize::new(n.get() as usize).expect("nonzero u32 workers"))
        .map_err(FingerprintWorkerError::ThreadCount)
}
pub(crate) fn ordered_bulk<I, T, E, F>(
    inputs: &[Option<I>],
    workers: NonZeroUsize,
    action: F,
) -> Result<Vec<Option<T>>, E>
where
    I: Sync,
    T: Send,
    E: Send + From<FingerprintWorkerError>,
    F: Fn(&I) -> Result<T, E> + Sync,
{
    // RDKit❗✔️: template <typename ReturnType, typename FuncType>
    // RDKit❗✔️: std::vector<std::unique_ptr<ReturnType>> mtgetFingerprints(
    // RDKit❗✔️:     FuncType func, const std::vector<const ROMol *> &mols, int numThreads) {
    // RDKit❗✔️:   std::vector<std::uint32_t> *fromAtoms = nullptr;
    // RDKit❗✔️:   std::vector<std::uint32_t> *ignoreAtoms = nullptr;
    // RDKit❗✔️:   std::vector<std::uint32_t> *customAtomInvariants = nullptr;
    // RDKit❗✔️:   std::vector<std::uint32_t> *customBondInvariants = nullptr;
    // RDKit❗✔️:   int confId = -1;
    // RDKit❗✔️:   AdditionalOutput *additionalOutput = nullptr;
    // RDKit❗✔️:   FingerprintFuncArguments args(fromAtoms, ignoreAtoms, confId,
    // RDKit❗✔️:                                 additionalOutput, customAtomInvariants,
    // RDKit❗✔️:                                 customBondInvariants);
    // RDKit❗✔️:
    // RDKit❗✔️:   std::vector<std::unique_ptr<ReturnType>> result;
    // RDKit❗✔️:   auto numThreadsToUse = getNumThreadsToUse(numThreads);
    // RDKit❗✔️:   unsigned int nmols = mols.size();
    // RDKit❗✔️:   result.reserve(nmols);
    // RDKit❗✔️:   if (numThreadsToUse == 1) {
    // RDKit❗✔️:     for (auto i = 0u; i < nmols; ++i) {
    // RDKit❗✔️:       if (!mols[i]) {
    // RDKit❗✔️:         result.emplace_back(std::unique_ptr<ReturnType>());
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         result.emplace_back(std::move(func(*mols[i], args)));
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗✔️:   else {
    // RDKit❗✔️:     std::vector<std::vector<std::unique_ptr<ReturnType>>> accum(
    // RDKit❗✔️:         numThreadsToUse);
    // RDKit❗✔️:     std::vector<std::thread> tg;
    // RDKit❗✔️:     for (auto ti = 0u; ti < numThreadsToUse; ++ti) {
    // RDKit❗✔️:       auto lfunc = [&](unsigned int tidx) {
    // RDKit❗✔️:         for (auto midx = tidx; midx < mols.size(); midx += numThreadsToUse) {
    // RDKit❗✔️:           if (!mols[midx]) {
    // RDKit❗✔️:             accum[tidx].emplace_back(std::unique_ptr<ReturnType>());
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             accum[tidx].emplace_back(std::move(func(*mols[midx], args)));
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       };
    // RDKit❗✔️:       tg.emplace_back(std::thread(lfunc, ti));
    // RDKit❗✔️:     }
    // RDKit❗✔️:     for (auto &thread : tg) {
    // RDKit❗✔️:       if (thread.joinable()) {
    // RDKit❗✔️:         thread.join();
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     for (auto midx = 0u; midx < mols.size(); ++midx) {
    // RDKit❗✔️:       auto tidx = midx % numThreadsToUse;
    // RDKit❗✔️:       auto jidx = midx / numThreadsToUse;
    // RDKit❗✔️:       result.emplace_back(std::move(accum[tidx][jidx]));
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: #endif
    // RDKit❗✔️:   return result;
    // RDKit❗✔️: }
    // Behavior: source worker strides and every missing slot are retained.
    // All launched workers join before structured launch/panic/calculation errors.
    // Complexity: O(N) ordered moves and O(W) worker bookkeeping; each callback
    // borrows the same configured state, without molecule/config reconstruction.
    let n = workers.get();
    if n == 1 {
        return inputs
            .iter()
            .map(|i| i.as_ref().map(&action).transpose())
            .collect();
    }
    std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(n);
        let mut launch_error = None;
        for ti in 0..n {
            let action = &action;
            match std::thread::Builder::new().spawn_scoped(scope, move || {
                (ti..inputs.len())
                    .step_by(n)
                    .map(|midx| inputs[midx].as_ref().map(action).transpose())
                    .collect::<Vec<_>>()
            }) {
                Ok(h) => handles.push(h),
                Err(e) => {
                    launch_error = Some(e);
                    break;
                }
            }
        }
        let mut accum = Vec::with_capacity(handles.len());
        let mut panic = false;
        for h in handles {
            match h.join() {
                Ok(rows) => accum.push(rows.into_iter()),
                Err(_) => panic = true,
            }
        }
        if let Some(e) = launch_error {
            return Err(FingerprintWorkerError::ThreadSpawn(e).into());
        }
        if panic {
            return Err(FingerprintWorkerError::Panic.into());
        }
        let mut result = Vec::with_capacity(inputs.len());
        for midx in 0..inputs.len() {
            result.push(
                accum[midx % n]
                    .next()
                    .ok_or_else(|| E::from(FingerprintWorkerError::Protocol))??,
            );
        }
        Ok(result)
    })
}
