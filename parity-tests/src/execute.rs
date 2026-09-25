use super::registry::{Input, Operation, Record, Value, Width};
use cosmolkit::{SparseCountFingerprint, SparseCountFingerprint32};

pub fn run(input: &Input) -> Result<Record, String> {
    macro_rules! execute {
        ($ty:ty, $key:ty) => {{
            let build = |entries: &[(u64, i32)]| -> Result<$ty, String> {
                let mut v = <$ty>::new(input.case.length as $key);
                for &(key, count) in entries {
                    if count == 0 {
                        v.set_value(key as $key, 1).map_err(|e| e.to_string())?;
                    }
                }
                // Materialize stored zeros first; never offset nonzero counts,
                // so i32 extrema need no arbitrary corpus restriction.
                v = v.with_added_scalar(-1).map_err(|e| e.to_string())?;
                for &(key, count) in entries {
                    if count != 0 {
                        v.set_value(key as $key, count).map_err(|e| e.to_string())?;
                    }
                }
                Ok(v)
            };
            let left = build(&input.case.left)?;
            let right = build(&input.case.right)?;
            let before = (left.clone(), right.clone());
            let result = match input.operation {
                Operation::FuzzyAnd => left.fuzzy_and(&right),
                Operation::FuzzyOr => left.fuzzy_or(&right),
            }
            .map_err(|e| e.to_string())?;
            if (left, right) != before {
                return Err("operands changed".into());
            }
            Value {
                length: result.length() as u64,
                entries: result
                    .nonzero_elements()
                    .iter()
                    .map(|(&k, &v)| (k as u64, v))
                    .collect(),
            }
        }};
    }
    let output = match input.width {
        Width::U32 => execute!(SparseCountFingerprint32, u32),
        Width::U64 => execute!(SparseCountFingerprint, u64),
    };
    Ok(Record {
        input: input.clone(),
        output,
    })
}
