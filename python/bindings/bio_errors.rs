use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

crate::enums::binding_enum!(BioStructureErrorKind, {
    RowIndexTooLarge,
    RowSpanOverflow,
    RowSpanOutOfBounds,
    TableTooLarge,
    NonContiguousSpan,
    IncompleteCoverage,
    ParentMismatch,
    RowReferenceOutOfBounds,
    CoordinateCountMismatch,
    EntitySubchainMismatch,
    EmptyResidueSpan,
    ImpossibleCrystalAngle,
    AtomNotFound,
});

crate::error_framework::error_projection! {record
pub struct BioStructureError {
    inner: ck::BioStructureError,
}
}

crate::error_framework::error_projection! {methods impl BioStructureError {
    #[getter]
    fn domain(&self) -> &'static str {
        "bio"
    }
    #[getter]
    fn variant(&self) -> BioStructureErrorKind {
        use BioStructureErrorKind as K;
        use ck::BioStructureError as E;
        match &self.inner {
            E::RowIndexTooLarge { value: _ } => K::RowIndexTooLarge,
            E::RowSpanOverflow { start: _, len: _ } => K::RowSpanOverflow,
            E::RowSpanOutOfBounds {
                start: _,
                len: _,
                table_len: _,
            } => K::RowSpanOutOfBounds,
            E::TableTooLarge { table: _, len: _ } => K::TableTooLarge,
            E::NonContiguousSpan {
                table: _,
                expected_start: _,
                actual_start: _,
            } => K::NonContiguousSpan,
            E::IncompleteCoverage {
                table: _,
                covered: _,
                table_len: _,
            } => K::IncompleteCoverage,
            E::ParentMismatch { table: _, index: _ } => K::ParentMismatch,
            E::RowReferenceOutOfBounds {
                table: _,
                index: _,
                table_len: _,
            } => K::RowReferenceOutOfBounds,
            E::CoordinateCountMismatch {
                atom_count: _,
                coordinate_count: _,
            } => K::CoordinateCountMismatch,
            E::EntitySubchainMismatch {
                entity_id: _,
                subchain: _,
            } => K::EntitySubchainMismatch,
            E::EmptyResidueSpan { operation: _ } => K::EmptyResidueSpan,
            E::ImpossibleCrystalAngle => K::ImpossibleCrystalAngle,
            E::AtomNotFound => K::AtomNotFound,
        }
    }
    #[getter]
    fn source_message(&self) -> Option<String> {
        self.inner.source().map(ToString::to_string)
    }
    #[getter]
    fn value(&self) -> Option<usize> {
        match self.inner {
            ck::BioStructureError::RowIndexTooLarge { value } => Some(value),
            _ => None,
        }
    }
    #[getter]
    fn start(&self) -> Option<u32> {
        match self.inner {
            ck::BioStructureError::RowSpanOverflow { start, .. }
            | ck::BioStructureError::RowSpanOutOfBounds { start, .. } => Some(start),
            _ => None,
        }
    }
    #[getter]
    fn len(&self) -> Option<usize> {
        match self.inner {
            ck::BioStructureError::RowSpanOverflow { len, .. }
            | ck::BioStructureError::RowSpanOutOfBounds { len, .. } => Some(len as usize),
            ck::BioStructureError::TableTooLarge { len, .. } => Some(len),
            _ => None,
        }
    }
    #[getter]
    fn table_len(&self) -> Option<usize> {
        match self.inner {
            ck::BioStructureError::RowSpanOutOfBounds { table_len, .. }
            | ck::BioStructureError::IncompleteCoverage { table_len, .. }
            | ck::BioStructureError::RowReferenceOutOfBounds { table_len, .. } => Some(table_len),
            _ => None,
        }
    }
    #[getter]
    fn table(&self) -> Option<&'static str> {
        match self.inner {
            ck::BioStructureError::TableTooLarge { table, .. }
            | ck::BioStructureError::NonContiguousSpan { table, .. }
            | ck::BioStructureError::IncompleteCoverage { table, .. }
            | ck::BioStructureError::ParentMismatch { table, .. }
            | ck::BioStructureError::RowReferenceOutOfBounds { table, .. } => Some(table),
            _ => None,
        }
    }
    #[getter]
    fn expected_start(&self) -> Option<u32> {
        match self.inner {
            ck::BioStructureError::NonContiguousSpan { expected_start, .. } => Some(expected_start),
            _ => None,
        }
    }
    #[getter]
    fn actual_start(&self) -> Option<u32> {
        match self.inner {
            ck::BioStructureError::NonContiguousSpan { actual_start, .. } => Some(actual_start),
            _ => None,
        }
    }
    #[getter]
    fn covered(&self) -> Option<u32> {
        match self.inner {
            ck::BioStructureError::IncompleteCoverage { covered, .. } => Some(covered),
            _ => None,
        }
    }
    #[getter]
    fn index(&self) -> Option<u32> {
        match self.inner {
            ck::BioStructureError::ParentMismatch { index, .. }
            | ck::BioStructureError::RowReferenceOutOfBounds { index, .. } => Some(index),
            _ => None,
        }
    }
    #[getter]
    fn atom_count(&self) -> Option<usize> {
        match self.inner {
            ck::BioStructureError::CoordinateCountMismatch { atom_count, .. } => Some(atom_count),
            _ => None,
        }
    }
    #[getter]
    fn coordinate_count(&self) -> Option<usize> {
        match self.inner {
            ck::BioStructureError::CoordinateCountMismatch {
                coordinate_count, ..
            } => Some(coordinate_count),
            _ => None,
        }
    }
    #[getter]
    fn entity_id(&self) -> Option<u32> {
        match self.inner {
            ck::BioStructureError::EntitySubchainMismatch { entity_id, .. } => {
                Some(entity_id.value())
            }
            _ => None,
        }
    }
    #[getter]
    fn subchain(&self) -> Option<&str> {
        match &self.inner {
            ck::BioStructureError::EntitySubchainMismatch { subchain, .. } => Some(subchain),
            _ => None,
        }
    }
    #[getter]
    fn operation(&self) -> Option<&'static str> {
        match self.inner {
            ck::BioStructureError::EmptyResidueSpan { operation } => Some(operation),
            _ => None,
        }
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}

impl BioStructureError {
    pub(crate) fn from_rust(inner: ck::BioStructureError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => crate::error_framework::from_payload(value),
            Err(error) => error,
        })
    }
}
