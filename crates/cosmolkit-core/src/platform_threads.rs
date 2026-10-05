//! Source RDThreads.h normalization and narrow pure-Rust CPU observation.
use std::{fmt, num::NonZeroU32};
#[derive(Debug)]
pub enum ThreadCountError {
    UndefinedSignedNegation,
    ObservationIo(std::io::Error),
    InvalidCpuList {
        segment: usize,
        reason: &'static str,
    },
    AutomaticObservationUnavailable {
        platform: &'static str,
    },
}
impl fmt::Display for ThreadCountError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UndefinedSignedNegation => {
                f.write_str("RDThreads::getNumThreadsToUse negation of INT_MIN is source-undefined")
            }
            Self::ObservationIo(e) => write!(f, "CPU online observation failed: {e}"),
            Self::InvalidCpuList { segment, reason } => {
                write!(f, "invalid CPU online list at segment {segment}: {reason}")
            }
            Self::AutomaticObservationUnavailable { platform } => write!(
                f,
                "source-backed automatic hardware count unavailable on {platform}; use an explicit positive thread count"
            ),
        }
    }
}
impl std::error::Error for ThreadCountError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::ObservationIo(e) => Some(e),
            _ => None,
        }
    }
}
/// Normalize source-defined thread selection against explicit observed hardware.
/// Observation zero is the C++ hardware_concurrency unknown sentinel.
pub fn rdkit_threads_with_observed_hardware(
    target: i32,
    observed_hardware: u32,
) -> Result<NonZeroU32, ThreadCountError> {
    // RDKit❗✔️: inline unsigned int getNumThreadsToUse(int target) {
    // RDKit❗✔️:   if (target >= 1) {
    // RDKit❗✔️:     return static_cast<unsigned int>(target);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   unsigned int res = std::thread::hardware_concurrency();
    // RDKit❗✔️:   if (res > rdcast<unsigned int>(-target)) {
    // RDKit❗✔️:     return res + target;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     return 1;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // Same unsigned comparison/subtraction with no allocation. Only source
    // actually executed signed INT_MIN negation has an explicit safety error.
    if target >= 1 {
        return Ok(NonZeroU32::new(target as u32).expect("positive source branch"));
    }
    let magnitude = target
        .checked_neg()
        .ok_or(ThreadCountError::UndefinedSignedNegation)? as u32;
    let count = if observed_hardware > magnitude {
        observed_hardware - magnitude
    } else {
        1
    };
    Ok(NonZeroU32::new(count).expect("strict greater branch or one"))
}
/// Count the kernel's sorted online CPU mask ranges, without expanding IDs.
/// Kernel ABI: https://docs.kernel.org/core-api/cpu_hotplug.html#using-cpu-hotplug
/// `online` represents cpu_online_mask, distinct from a process affinity mask.
pub fn cpu_online_count(list: &str) -> Result<u32, ThreadCountError> {
    let list = list.trim();
    let mut previous_end = None;
    let mut count = 0_u32;
    for (segment, item) in list.split(',').enumerate() {
        let invalid = |reason| ThreadCountError::InvalidCpuList { segment, reason };
        let mut range = item.split('-');
        let first = range.next().unwrap_or("");
        let parse = |text: &str| -> Result<u32, ThreadCountError> {
            if text.is_empty() || !text.bytes().all(|b| b.is_ascii_digit()) {
                return Err(invalid("CPU id must contain decimal digits"));
            }
            text.parse().map_err(|_| invalid("CPU id exceeds u32"))
        };
        let begin = parse(first)?;
        let end = match range.next() {
            None => begin,
            Some(last) => parse(last)?,
        };
        if range.next().is_some() {
            return Err(invalid("CPU interval has more than one separator"));
        }
        if begin > end {
            return Err(invalid("CPU interval is reversed"));
        }
        if previous_end.is_some_and(|last| begin <= last) {
            return Err(invalid("CPU intervals overlap or are unordered"));
        }
        let width = end
            .checked_sub(begin)
            .and_then(|n| n.checked_add(1))
            .ok_or_else(|| invalid("CPU count exceeds u32"))?;
        count = count
            .checked_add(width)
            .ok_or_else(|| invalid("CPU count exceeds u32"))?;
        previous_end = Some(end);
    }
    Ok(count)
}
/// Read system online CPUs, not process quota/affinity. No cached observation:
/// CPU hotplug affects the next automatic thread-selection call.
pub fn observe_hardware_threads() -> Result<u32, ThreadCountError> {
    #[cfg(target_os = "linux")]
    {
        let online = std::fs::read_to_string("/sys/devices/system/cpu/online")
            .map_err(ThreadCountError::ObservationIo)?;
        cpu_online_count(&online)
    }
    #[cfg(not(target_os = "linux"))]
    {
        Err(ThreadCountError::AutomaticObservationUnavailable {
            platform: std::env::consts::OS,
        })
    }
}
/// Positive source worker counts never require system observation.
pub fn rdkit_thread_count(target: i32) -> Result<NonZeroU32, ThreadCountError> {
    if target >= 1 {
        return rdkit_threads_with_observed_hardware(target, 0);
    }
    if target == i32::MIN {
        return Err(ThreadCountError::UndefinedSignedNegation);
    }
    rdkit_threads_with_observed_hardware(target, observe_hardware_threads()?)
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn source_normalization_positive_explicit_and_unknown_observation() {
        for hardware in [0, 1, 224, u32::MAX] {
            for target in [1, 2, 7, i32::MAX] {
                assert_eq!(
                    rdkit_threads_with_observed_hardware(target, hardware)
                        .unwrap()
                        .get(),
                    target as u32
                );
            }
        }
        for target in [0, -1, -223, i32::MIN + 1] {
            assert_eq!(
                rdkit_threads_with_observed_hardware(target, 0)
                    .unwrap()
                    .get(),
                1
            );
        }
        assert!(matches!(
            rdkit_threads_with_observed_hardware(i32::MIN, 224),
            Err(ThreadCountError::UndefinedSignedNegation)
        ));
    }
    #[test]
    fn source_normalization_hardware_comparison_boundaries() {
        for (target, hardware, expected) in [
            (0, 224, 224),
            (-1, 224, 223),
            (-223, 224, 1),
            (-224, 224, 1),
            (-225, 224, 1),
            (0, 1, 1),
            (-1, 1, 1),
            (-1, u32::MAX, u32::MAX - 1),
            (-2147483647, u32::MAX, 2147483648),
        ] {
            assert_eq!(
                rdkit_threads_with_observed_hardware(target, hardware)
                    .unwrap()
                    .get(),
                expected
            );
        }
    }
    #[test]
    fn kernel_cpu_online_ranges_count_holes_without_expansion() {
        for (list, count) in [
            ("0", 1),
            ("0-223\n", 224),
            ("0-3,8-11,17", 9),
            ("2,5,7-10", 6),
            ("0-4294967294", u32::MAX),
            ("4294967295", 1),
        ] {
            assert_eq!(cpu_online_count(list).unwrap(), count, "{list}");
        }
    }
    #[test]
    fn kernel_cpu_online_malformed_or_overflow_masks_fail_structurally() {
        for list in [
            "",
            " ",
            "-1",
            "+1",
            "1-",
            "3-2",
            "1--2",
            "0,0",
            "0-2,2-4",
            "2,1",
            "1,,3",
            "0-4294967295",
            "4294967296",
            "0,4294967295,4294967296",
            "1 2",
        ] {
            assert!(
                matches!(
                    cpu_online_count(list),
                    Err(ThreadCountError::InvalidCpuList { .. })
                ),
                "{list}"
            );
        }
    }
}
