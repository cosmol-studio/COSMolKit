//! Ordered dense fingerprint results, including failed input positions.

use crate::Fingerprint;

/// An ordered, immutable fingerprint collection. Missing rows are preserved.
/// The first present row defines the matrix width; empty selections retain it.
#[derive(Clone, Debug, PartialEq)]
pub struct FingerprintBatch {
    values: Vec<Option<Fingerprint>>,
    n_bits: u32,
}

impl FingerprintBatch {
    pub fn new(values: Vec<Option<Fingerprint>>) -> Self {
        let n_bits = values
            .iter()
            .flatten()
            .next()
            .map_or(0, Fingerprint::n_bits);
        Self { values, n_bits }
    }

    pub fn len(&self) -> usize {
        self.values.len()
    }
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }
    pub fn n_bits(&self) -> u32 {
        self.n_bits
    }
    pub fn values(&self) -> &[Option<Fingerprint>] {
        &self.values
    }
    pub fn into_values(self) -> Vec<Option<Fingerprint>> {
        self.values
    }

    /// Select rows in the supplied order without dropping missing inputs.
    pub fn select(&self, indices: &[usize]) -> Option<Self> {
        let values = indices
            .iter()
            .map(|&i| self.values.get(i).cloned())
            .collect::<Option<Vec<_>>>()?;
        Some(Self {
            values,
            n_bits: self.n_bits,
        })
    }

    /// Validate a rectangular, complete matrix before allocating an export.
    pub fn validate_matrix(&self) -> Result<(), FingerprintBatchError> {
        for (index, value) in self.values.iter().enumerate() {
            let value = value
                .as_ref()
                .ok_or(FingerprintBatchError::MissingRow { index })?;
            if value.n_bits() != self.n_bits {
                return Err(FingerprintBatchError::WidthMismatch {
                    index,
                    actual: value.n_bits(),
                    expected: self.n_bits,
                });
            }
        }
        Ok(())
    }
}

impl From<Vec<Option<Fingerprint>>> for FingerprintBatch {
    fn from(values: Vec<Option<Fingerprint>>) -> Self {
        Self::new(values)
    }
}
impl FromIterator<Option<Fingerprint>> for FingerprintBatch {
    fn from_iter<T: IntoIterator<Item = Option<Fingerprint>>>(values: T) -> Self {
        Self::new(values.into_iter().collect())
    }
}
impl AsRef<[Option<Fingerprint>]> for FingerprintBatch {
    fn as_ref(&self) -> &[Option<Fingerprint>] {
        self.values()
    }
}
impl PartialEq<Vec<Option<Fingerprint>>> for FingerprintBatch {
    fn eq(&self, other: &Vec<Option<Fingerprint>>) -> bool {
        self.values() == other
    }
}
impl PartialEq<FingerprintBatch> for Vec<Option<Fingerprint>> {
    fn eq(&self, other: &FingerprintBatch) -> bool {
        self == other.values()
    }
}
impl std::ops::Deref for FingerprintBatch {
    type Target = [Option<Fingerprint>];
    fn deref(&self) -> &Self::Target {
        self.values()
    }
}
impl IntoIterator for FingerprintBatch {
    type Item = Option<Fingerprint>;
    type IntoIter = std::vec::IntoIter<Self::Item>;
    fn into_iter(self) -> Self::IntoIter {
        self.values.into_iter()
    }
}
impl<'a> IntoIterator for &'a FingerprintBatch {
    type Item = &'a Option<Fingerprint>;
    type IntoIter = std::slice::Iter<'a, Option<Fingerprint>>;
    fn into_iter(self) -> Self::IntoIter {
        self.values.iter()
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FingerprintBatchError {
    MissingRow {
        index: usize,
    },
    WidthMismatch {
        index: usize,
        actual: u32,
        expected: u32,
    },
}
impl std::fmt::Display for FingerprintBatchError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MissingRow { index } => write!(f, "fingerprint row {index} is None"),
            Self::WidthMismatch {
                index,
                actual,
                expected,
            } => write!(
                f,
                "fingerprint row {index} has {actual} bits; expected {expected}"
            ),
        }
    }
}
impl std::error::Error for FingerprintBatchError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn selection_preserves_positions_width_and_matrix_errors() {
        let fingerprint = Fingerprint::from_on_bits(65, [0, 64]).unwrap();
        let batch = FingerprintBatch::new(vec![None, Some(fingerprint.clone())]);
        assert_eq!(batch.n_bits(), 65);
        assert_eq!(
            batch.validate_matrix(),
            Err(FingerprintBatchError::MissingRow { index: 0 })
        );
        let selected = batch.select(&[1, 1]).unwrap();
        assert_eq!(
            selected.values(),
            &[Some(fingerprint.clone()), Some(fingerprint)]
        );
        assert!(selected.validate_matrix().is_ok());
        let empty = batch.select(&[]).unwrap();
        assert!(empty.is_empty());
        assert_eq!(empty.n_bits(), 65);
        assert!(empty.validate_matrix().is_ok());
        assert!(batch.select(&[2]).is_none());
        let mixed =
            FingerprintBatch::new(vec![Some(Fingerprint::new(8)), Some(Fingerprint::new(9))]);
        assert_eq!(
            mixed.validate_matrix(),
            Err(FingerprintBatchError::WidthMismatch {
                index: 1,
                actual: 9,
                expected: 8
            })
        );
    }
}
