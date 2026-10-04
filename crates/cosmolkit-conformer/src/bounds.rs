//! Dense ordered distance bounds used by the private distance-geometry code.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum BoundsMatrixError {
    DimensionOverflow {
        dimension: usize,
    },
    AllocationFailed {
        elements: usize,
    },
    DataLengthTooShort {
        expected: usize,
        actual: usize,
    },
    IndexOutOfBounds {
        axis: MatrixAxis,
        index: usize,
        dimension: usize,
    },
    InvalidBound {
        kind: BoundKind,
        row: usize,
        column: usize,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum MatrixAxis {
    Row,
    Column,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum BoundKind {
    Upper,
    Lower,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct BoundsMatrix {
    data: Vec<f64>,
    dimension: usize,
}

impl BoundsMatrix {
    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::BoundsMatrix (BoundsMatrix.h:31-32)
    // RDKit❗✔️: explicit BoundsMatrix(unsigned int N)
    // RDKit❗✔️:     : RDNumeric::SquareMatrix<double>(N, 0.0) {}
    // BEGIN RDKIT CPP FUNCTION RDNumeric::SquareMatrix::SquareMatrix (SquareMatrix.h:25)
    // RDKit❗✔️: SquareMatrix(unsigned int N, TYPE val) : Matrix<TYPE>(N, N, val) {}
    // BEGIN RDKIT CPP FUNCTION RDNumeric::Matrix::Matrix (Matrix.h:41-49)
    // RDKit❗✔️: Matrix(unsigned int nRows, unsigned int nCols, TYPE val)
    // RDKit❗✔️:     : d_nRows(nRows), d_nCols(nCols), d_dataSize(nRows * nCols) {
    // RDKit❗✔️:   TYPE *data = new TYPE[d_dataSize];
    // RDKit❗✔️:   unsigned int i;
    // RDKit❗✔️:   for (i = 0; i < d_dataSize; i++) {
    // RDKit❗✔️:     data[i] = val;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   d_data.reset(data);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDNumeric::Matrix::Matrix
    // END RDKIT CPP FUNCTION RDNumeric::SquareMatrix::SquareMatrix
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::BoundsMatrix
    pub(crate) fn new(dimension: usize) -> Result<Self, BoundsMatrixError> {
        let elements = dimension
            .checked_mul(dimension)
            .ok_or(BoundsMatrixError::DimensionOverflow { dimension })?;
        let mut data = Vec::new();
        data.try_reserve_exact(elements)
            .map_err(|_| BoundsMatrixError::AllocationFailed { elements })?;
        data.resize(elements, 0.0);
        Ok(Self { data, dimension })
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::BoundsMatrix (BoundsMatrix.h:33-34)
    // RDKit❗✔️: BoundsMatrix(unsigned int N, DATA_SPTR data)
    // RDKit❗✔️:     : RDNumeric::SquareMatrix<double>(N, data) {}
    // BEGIN RDKIT CPP FUNCTION RDNumeric::SquareMatrix::SquareMatrix (SquareMatrix.h:27-28)
    // RDKit❗✔️: SquareMatrix(unsigned int N, typename Matrix<TYPE>::DATA_SPTR data)
    // RDKit❗✔️:     : Matrix<TYPE>(N, N, data) {}
    // BEGIN RDKIT CPP FUNCTION RDNumeric::Matrix::Matrix (Matrix.h:56-59)
    // RDKit❗✔️: Matrix(unsigned int nRows, unsigned int nCols, DATA_SPTR data)
    // RDKit❗✔️:     : d_nRows(nRows), d_nCols(nCols), d_dataSize(nRows * nCols) {
    // RDKit❗✔️:   d_data = data;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDNumeric::Matrix::Matrix
    // END RDKIT CPP FUNCTION RDNumeric::SquareMatrix::SquareMatrix
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::BoundsMatrix
    // The only pinned wrapper caller copies the caller's array into a new
    // buffer and transfers its sole shared owner here (rdDistGeom.cpp:88-98).
    // Taking that buffer by Vec value preserves this call path without raw
    // pointer aliasing or a per-access lock.
    pub(crate) fn from_data(
        dimension: usize,
        mut data: Vec<f64>,
    ) -> Result<Self, BoundsMatrixError> {
        let elements = dimension
            .checked_mul(dimension)
            .ok_or(BoundsMatrixError::DimensionOverflow { dimension })?;
        if data.len() < elements {
            return Err(BoundsMatrixError::DataLengthTooShort {
                expected: elements,
                actual: data.len(),
            });
        }
        data.truncate(elements);
        Ok(Self { data, dimension })
    }

    fn index(&self, row: usize, column: usize) -> Result<usize, BoundsMatrixError> {
        if row >= self.dimension {
            return Err(BoundsMatrixError::IndexOutOfBounds {
                axis: MatrixAxis::Row,
                index: row,
                dimension: self.dimension,
            });
        }
        if column >= self.dimension {
            return Err(BoundsMatrixError::IndexOutOfBounds {
                axis: MatrixAxis::Column,
                index: column,
                dimension: self.dimension,
            });
        }
        Ok(row * self.dimension + column)
    }

    /// Reads one dense matrix entry without repeating the caller's index checks.
    ///
    /// # Safety
    ///
    /// `row` and `column` must both be less than `dimension()`. The backing
    /// vector must retain the exact square length established by `new` or
    /// `from_data`.
    #[inline]
    pub(super) unsafe fn get_val_unchecked(&self, row: usize, column: usize) -> f64 {
        let index = row * self.dimension + column;
        // SAFETY: the caller proves both coordinates are in range; constructors
        // check dimension multiplication and keep at least that many elements.
        unsafe { *self.data.get_unchecked(index) }
    }

    /// Writes one dense matrix entry without applying bound validation.
    ///
    /// # Safety
    ///
    /// `row` and `column` must both be less than `dimension()`. The backing
    /// vector must retain the exact square length established by `new` or
    /// `from_data`.
    #[inline]
    pub(super) unsafe fn set_val_unchecked(&mut self, row: usize, column: usize, value: f64) {
        let index = row * self.dimension + column;
        // SAFETY: the caller proves both coordinates are in range; constructors
        // check dimension multiplication and keep at least that many elements.
        unsafe { *self.data.get_unchecked_mut(index) = value };
    }

    // BEGIN RDKIT CPP FUNCTION RDNumeric::Matrix::getVal (Matrix.h:94-99)
    // RDKit❗✔️: inline virtual TYPE getVal(unsigned int i, unsigned int j) const {
    // RDKit❗✔️:   PRECONDITION(i < d_nRows, "bad index");
    // RDKit❗✔️:   PRECONDITION(j < d_nCols, "bad index");
    // RDKit❗✔️:   unsigned int id = i * d_nCols + j;
    // RDKit❗✔️:   return d_data[id];
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDNumeric::Matrix::getVal
    pub(crate) fn get_val(&self, row: usize, column: usize) -> Result<f64, BoundsMatrixError> {
        Ok(self.data[self.index(row, column)?])
    }

    // BEGIN RDKIT CPP FUNCTION RDNumeric::Matrix::setVal (Matrix.h:102-108)
    // RDKit❗✔️: inline virtual void setVal(unsigned int i, unsigned int j, TYPE val) {
    // RDKit❗✔️:   PRECONDITION(i < d_nRows, "bad index");
    // RDKit❗✔️:   PRECONDITION(j < d_nCols, "bad index");
    // RDKit❗✔️:   unsigned int id = i * d_nCols + j;
    // RDKit❗✔️:   d_data[id] = val;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDNumeric::Matrix::setVal
    pub(crate) fn set_val(
        &mut self,
        row: usize,
        column: usize,
        value: f64,
    ) -> Result<(), BoundsMatrixError> {
        let index = self.index(row, column)?;
        self.data[index] = value;
        Ok(())
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::getUpperBound (BoundsMatrix.h:37-43)
    // RDKit❗✔️: inline double getUpperBound(unsigned int i, unsigned int j) const {
    // RDKit❗✔️:   if (i < j) {
    // RDKit❗✔️:     return getVal(i, j);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     return getVal(j, i);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::getUpperBound
    pub(crate) fn get_upper(&self, row: usize, column: usize) -> Result<f64, BoundsMatrixError> {
        if row < column {
            self.get_val(row, column)
        } else {
            self.get_val(column, row)
        }
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setUpperBound (BoundsMatrix.h:46-53)
    // RDKit❗✔️: inline void setUpperBound(unsigned int i, unsigned int j, double val) {
    // RDKit❗✔️:   CHECK_INVARIANT(val >= 0.0, "Negative upper bound");
    // RDKit❗✔️:   if (i < j) {
    // RDKit❗✔️:     setVal(i, j, val);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     setVal(j, i, val);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setUpperBound
    pub(crate) fn set_upper(
        &mut self,
        row: usize,
        column: usize,
        value: f64,
    ) -> Result<(), BoundsMatrixError> {
        if !(value >= 0.0) {
            return Err(BoundsMatrixError::InvalidBound {
                kind: BoundKind::Upper,
                row,
                column,
            });
        }
        if row < column {
            self.set_val(row, column, value)
        } else {
            self.set_val(column, row, value)
        }
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setUpperBoundIfBetter (BoundsMatrix.h:57-62)
    // RDKit❗✔️: inline void setUpperBoundIfBetter(unsigned int i, unsigned int j,
    // RDKit❗✔️:                                   double val) {
    // RDKit❗✔️:   if ((val < getUpperBound(i, j)) && (val > getLowerBound(i, j))) {
    // RDKit❗✔️:     setUpperBound(i, j, val);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setUpperBoundIfBetter
    pub(crate) fn set_upper_if_better(
        &mut self,
        row: usize,
        column: usize,
        value: f64,
    ) -> Result<(), BoundsMatrixError> {
        if value < self.get_upper(row, column)? && value > self.get_lower(row, column)? {
            self.set_upper(row, column, value)?;
        }
        Ok(())
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::getLowerBound (BoundsMatrix.h:84-90)
    // RDKit❗✔️: inline double getLowerBound(unsigned int i, unsigned int j) const {
    // RDKit❗✔️:   if (i < j) {
    // RDKit❗✔️:     return getVal(j, i);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     return getVal(i, j);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::getLowerBound
    pub(crate) fn get_lower(&self, row: usize, column: usize) -> Result<f64, BoundsMatrixError> {
        if row < column {
            self.get_val(column, row)
        } else {
            self.get_val(row, column)
        }
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setLowerBound (BoundsMatrix.h:65-72)
    // RDKit❗✔️: inline void setLowerBound(unsigned int i, unsigned int j, double val) {
    // RDKit❗✔️:   CHECK_INVARIANT(val >= 0.0, "Negative lower bound");
    // RDKit❗✔️:   if (i < j) {
    // RDKit❗✔️:     setVal(j, i, val);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     setVal(i, j, val);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setLowerBound
    pub(crate) fn set_lower(
        &mut self,
        row: usize,
        column: usize,
        value: f64,
    ) -> Result<(), BoundsMatrixError> {
        if !(value >= 0.0) {
            return Err(BoundsMatrixError::InvalidBound {
                kind: BoundKind::Lower,
                row,
                column,
            });
        }
        if row < column {
            self.set_val(column, row, value)
        } else {
            self.set_val(row, column, value)
        }
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setLowerBoundIfBetter (BoundsMatrix.h:76-81)
    // RDKit❗✔️: inline void setLowerBoundIfBetter(unsigned int i, unsigned int j,
    // RDKit❗✔️:                                   double val) {
    // RDKit❗✔️:   if ((val > getLowerBound(i, j)) && (val < getUpperBound(i, j))) {
    // RDKit❗✔️:     setLowerBound(i, j, val);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setLowerBoundIfBetter
    pub(crate) fn set_lower_if_better(
        &mut self,
        row: usize,
        column: usize,
        value: f64,
    ) -> Result<(), BoundsMatrixError> {
        if value > self.get_lower(row, column)? && value < self.get_upper(row, column)? {
            self.set_lower(row, column, value)?;
        }
        Ok(())
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::checkValid (BoundsMatrix.h:94-103)
    // RDKit❗✔️: inline bool checkValid() const {
    // RDKit❗✔️:   unsigned int i, j;
    // RDKit❗✔️:   for (i = 1; i < d_nRows; i++) {
    // RDKit❗✔️:     for (j = 0; j < i; j++) {
    // RDKit❗✔️:       if (getUpperBound(i, j) < getLowerBound(i, j)) {
    // RDKit❗✔️:         return false;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return true;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::checkValid
    pub(crate) fn check_valid(&self) -> bool {
        for row in 1..self.dimension {
            for column in 0..row {
                let upper = self
                    .get_upper(row, column)
                    .expect("check_valid visits only matrix indices");
                let lower = self
                    .get_lower(row, column)
                    .expect("check_valid visits only matrix indices");
                if upper < lower {
                    return false;
                }
            }
        }
        true
    }

    pub(crate) fn dimension(&self) -> usize {
        self.dimension
    }
}

#[cfg(test)]
mod tests {
    use super::{BoundKind, BoundsMatrix, BoundsMatrixError, MatrixAxis};

    #[test]
    fn cf3d_c01_zero_one_two_points_keep_dense_ordered_triangles() {
        let empty = BoundsMatrix::new(0).unwrap();
        assert_eq!(empty.dimension(), 0);
        assert!(empty.check_valid());
        assert!(matches!(
            empty.get_val(0, 0),
            Err(BoundsMatrixError::IndexOutOfBounds {
                axis: MatrixAxis::Row,
                ..
            })
        ));

        let mut one = BoundsMatrix::new(1).unwrap();
        assert_eq!(one.get_val(0, 0).unwrap(), 0.0);
        one.set_upper(0, 0, 4.0).unwrap();
        assert_eq!(one.get_upper(0, 0).unwrap(), 4.0);
        assert_eq!(one.get_lower(0, 0).unwrap(), 4.0);
        one.set_lower(0, 0, 7.0).unwrap();
        assert_eq!(one.get_upper(0, 0).unwrap(), 7.0);
        assert_eq!(one.get_lower(0, 0).unwrap(), 7.0);
        assert!(one.check_valid());

        let mut two = BoundsMatrix::new(2).unwrap();
        assert_eq!(two.get_val(0, 1).unwrap(), 0.0);
        assert_eq!(two.get_val(1, 0).unwrap(), 0.0);
        two.set_upper(0, 1, 8.0).unwrap();
        two.set_lower(0, 1, 2.0).unwrap();
        assert_eq!(two.get_val(0, 1).unwrap(), 8.0);
        assert_eq!(two.get_val(1, 0).unwrap(), 2.0);
        assert_eq!(two.get_upper(0, 1).unwrap(), 8.0);
        assert_eq!(two.get_upper(1, 0).unwrap(), 8.0);
        assert_eq!(two.get_lower(0, 1).unwrap(), 2.0);
        assert_eq!(two.get_lower(1, 0).unwrap(), 2.0);

        two.set_upper(1, 0, 9.0).unwrap();
        assert_eq!(two.get_val(0, 1).unwrap(), 9.0);
        assert_eq!(two.get_val(1, 0).unwrap(), 2.0);
        assert!(two.check_valid());

        // RDKit BoundsMatrix.h:57-62 and 76-81 use strict inequalities.
        two.set_upper_if_better(0, 1, 5.0).unwrap();
        two.set_lower_if_better(0, 1, 3.0).unwrap();
        assert_eq!(two.get_upper(0, 1).unwrap(), 5.0);
        assert_eq!(two.get_lower(0, 1).unwrap(), 3.0);
        assert!(two.check_valid());
        two.set_lower(0, 1, 6.0).unwrap();
        assert!(!two.check_valid());
    }

    #[test]
    fn cf3d_c01_if_better_uses_strict_source_bounds() {
        // RDKit BoundsMatrix.h:57-62: val < upper && val > lower.
        let mut bounds = BoundsMatrix::new(2).unwrap();
        bounds.set_lower(0, 1, 2.0).unwrap();
        bounds.set_upper(0, 1, 8.0).unwrap();
        for value in [2.0, 8.0, 1.0, 9.0, f64::NAN] {
            bounds.set_upper_if_better(0, 1, value).unwrap();
            assert_eq!(bounds.get_upper(0, 1).unwrap(), 8.0);
        }
        bounds.set_upper_if_better(0, 1, 5.0).unwrap();
        assert_eq!(bounds.get_upper(0, 1).unwrap(), 5.0);

        // RDKit BoundsMatrix.h:76-81: val > lower && val < upper.
        for value in [2.0, 5.0, 1.0, 6.0, f64::NAN] {
            bounds.set_lower_if_better(0, 1, value).unwrap();
            assert_eq!(bounds.get_lower(0, 1).unwrap(), 2.0);
        }
        bounds.set_lower_if_better(0, 1, 3.0).unwrap();
        assert_eq!(bounds.get_lower(0, 1).unwrap(), 3.0);
    }

    #[test]
    fn cf3d_c01_checked_indexing_reports_first_invalid_axis() {
        let mut bounds = BoundsMatrix::new(2).unwrap();
        assert!(matches!(
            bounds.get_val(2, 0),
            Err(BoundsMatrixError::IndexOutOfBounds {
                axis: MatrixAxis::Row,
                index: 2,
                dimension: 2,
            })
        ));
        assert!(matches!(
            bounds.get_val(0, 2),
            Err(BoundsMatrixError::IndexOutOfBounds {
                axis: MatrixAxis::Column,
                index: 2,
                dimension: 2,
            })
        ));
        assert!(matches!(
            bounds.set_val(0, 2, 1.0),
            Err(BoundsMatrixError::IndexOutOfBounds {
                axis: MatrixAxis::Column,
                ..
            })
        ));
    }

    #[test]
    fn cf3d_c01_data_constructor_and_clone_preserve_independent_storage() {
        // The pinned wrapper supplies a row-major copied array to the
        // DATA_SPTR constructor; the Rust private owner consumes that buffer.
        let original = BoundsMatrix::from_data(2, vec![1.0, 2.0, 3.0, 4.0]).unwrap();
        assert_eq!(original.get_val(0, 0).unwrap(), 1.0);
        assert_eq!(original.get_val(0, 1).unwrap(), 2.0);
        assert_eq!(original.get_val(1, 0).unwrap(), 3.0);
        assert_eq!(original.get_val(1, 1).unwrap(), 4.0);

        let mut cloned = original.clone();
        cloned.set_val(0, 1, 20.0).unwrap();
        assert_eq!(cloned.get_val(0, 1).unwrap(), 20.0);
        assert_eq!(original.get_val(0, 1).unwrap(), 2.0);

        assert_eq!(
            BoundsMatrix::from_data(2, vec![1.0, 2.0, 3.0]),
            Err(BoundsMatrixError::DataLengthTooShort {
                expected: 4,
                actual: 3,
            })
        );
    }

    #[test]
    fn cf3d_c01_bound_checks_match_nonnegative_source_invariant() {
        let mut bounds = BoundsMatrix::new(2).unwrap();
        for value in [-1.0, f64::NAN] {
            assert!(matches!(
                bounds.set_upper(0, 1, value),
                Err(BoundsMatrixError::InvalidBound {
                    kind: BoundKind::Upper,
                    row: 0,
                    column: 1,
                })
            ));
            assert!(matches!(
                bounds.set_lower(0, 1, value),
                Err(BoundsMatrixError::InvalidBound {
                    kind: BoundKind::Lower,
                    row: 0,
                    column: 1,
                })
            ));
        }

        bounds.set_upper(0, 1, -0.0).unwrap();
        assert!(bounds.get_upper(0, 1).unwrap().is_sign_negative());
        bounds.set_lower(0, 1, -0.0).unwrap();
        assert!(bounds.get_lower(0, 1).unwrap().is_sign_negative());
    }

    #[test]
    fn cf3d_c01_check_valid_keeps_source_scan_and_nan_comparison() {
        let mut bounds = BoundsMatrix::new(2).unwrap();
        bounds.set_upper(0, 1, 4.0).unwrap();
        bounds.set_lower(0, 1, 4.0).unwrap();
        assert!(bounds.check_valid());

        bounds.set_val(1, 0, f64::NAN).unwrap();
        assert!(bounds.check_valid());
    }
}
