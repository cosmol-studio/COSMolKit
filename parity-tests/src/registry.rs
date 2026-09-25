//! The pilot's only operation/width/corpus declaration. No timing matrix.
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Operation {
    FuzzyAnd,
    FuzzyOr,
}

impl Operation {
    pub fn name(self) -> &'static str {
        match self {
            Self::FuzzyAnd => "fuzzy_and",
            Self::FuzzyOr => "fuzzy_or",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Width {
    U32,
    U64,
}

pub struct Task {
    pub operation: Operation,
    pub widths: &'static [Width],
}

pub const TASKS: &[Task] = &[
    Task {
        operation: Operation::FuzzyAnd,
        widths: &[Width::U32, Width::U64],
    },
    Task {
        operation: Operation::FuzzyOr,
        widths: &[Width::U32, Width::U64],
    },
];

pub const RDKIT_VERSION: &str = "2026.03.1";

// Only typed integer indices/counts cross the oracle boundary.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Pair {
    pub id: String,
    pub length: u64,
    pub left: Vec<(u64, i32)>,
    pub right: Vec<(u64, i32)>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Input {
    pub case: Pair,
    pub operation: Operation,
    pub width: Width,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Value {
    pub length: u64,
    pub entries: Vec<(u64, i32)>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Record {
    pub input: Input,
    pub output: Value,
}

pub fn select(name: Option<&str>) -> Result<Vec<&'static Task>, String> {
    let selected: Vec<_> = TASKS
        .iter()
        .filter(|t| name.is_none_or(|n| n == t.operation.name()))
        .collect();
    if selected.is_empty() {
        return Err(format!("unknown task: {name:?}"));
    }
    Ok(selected)
}

pub fn validate(cases: &[Pair], tasks: &[&Task]) -> Result<(), String> {
    use std::collections::BTreeSet;
    if cases.is_empty() {
        return Err("corpus is empty".into());
    }
    let mut ids = BTreeSet::new();
    for c in cases {
        if c.id.is_empty() || !ids.insert(&c.id) {
            return Err("empty/duplicate case ID".into());
        }
        for width in tasks.iter().flat_map(|t| t.widths) {
            let max = match width {
                Width::U32 => u32::MAX as u64,
                Width::U64 => u64::MAX,
            };
            if c.length > max {
                return Err(format!("{}: length does not fit {width:?}", c.id));
            }
            for values in [&c.left, &c.right] {
                let mut keys = BTreeSet::new();
                for &(key, _) in values {
                    if !keys.insert(key) || key >= c.length {
                        return Err(format!("{}: invalid/duplicate index", c.id));
                    }
                }
            }
        }
    }
    Ok(())
}

pub fn expand(cases: &[Pair], task: &Task) -> Vec<Input> {
    cases
        .iter()
        .flat_map(|case| {
            task.widths.iter().map(move |&width| Input {
                case: case.clone(),
                operation: task.operation,
                width,
            })
        })
        .collect()
}

pub fn builtin() -> Vec<Pair> {
    // Each row names a source branch, not an arbitrary sampling size.
    [
        ("empty_both", vec![], vec![]),
        ("left_empty_right_tail", vec![], vec![(2, 3)]),
        ("right_empty_remove_left", vec![(2, -3)], vec![]),
        (
            "shared_signed_min_max",
            vec![(1, 5), (3, -2), (8, 4)],
            vec![(1, 3), (3, -4), (9, 7)],
        ),
        (
            "disjoint_interleaved",
            vec![(1, 2), (3, -1)],
            vec![(0, 4), (2, -3), (4, 5)],
        ),
        (
            "explicit_zero_shared_and_exclusive",
            vec![(1, 0), (3, 0)],
            vec![(1, -2), (4, 0)],
        ),
    ]
    .into_iter()
    .map(|(id, left, right)| Pair {
        id: id.into(),
        length: 16,
        left,
        right,
    })
    .collect()
}
