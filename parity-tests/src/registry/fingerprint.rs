//! Sparse fingerprint parameters and validation; not cache policy.
use super::{FingerprintInput, FingerprintValue, Input, Operation, Pair, Width};
use std::collections::BTreeSet;

pub const WIDTHS: &[Width] = &[Width::U32, Width::U64];

pub fn validate(cases: &[Pair]) -> Result<(), String> {
    let mut ids = BTreeSet::new();
    for case in cases {
        if case.id.is_empty() || !ids.insert(&case.id) {
            return Err("empty/duplicate case ID".into());
        }
        for width in WIDTHS {
            let max = match width {
                Width::U32 => u32::MAX as u64,
                Width::U64 => u64::MAX,
            };
            if case.length > max {
                return Err(format!("{}: length does not fit {width:?}", case.id));
            }
        }
        for values in [&case.left, &case.right] {
            let mut keys = BTreeSet::new();
            for &(key, _) in values {
                if !keys.insert(key) || key >= case.length {
                    return Err(format!("{}: invalid/duplicate index", case.id));
                }
            }
        }
    }
    Ok(())
}

pub fn expand(cases: &[Pair], operation: Operation) -> Vec<Input> {
    cases
        .iter()
        .flat_map(|case| {
            WIDTHS.iter().map(move |&width| {
                Input::Fingerprint(FingerprintInput {
                    case: case.clone(),
                    operation,
                    width,
                })
            })
        })
        .collect()
}

pub fn validate_output(input: &FingerprintInput, output: &FingerprintValue) -> Result<(), String> {
    if output.length != input.case.length
        || output.entries.windows(2).any(|w| w[0].0 >= w[1].0)
        || output
            .entries
            .iter()
            .any(|&(key, _)| key >= input.case.length)
    {
        return Err(format!("{}: malformed reference output", input.case.id));
    }
    Ok(())
}
