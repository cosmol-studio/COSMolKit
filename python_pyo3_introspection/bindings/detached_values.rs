use crate::chemical_values::EnumValue;
use crate::coordinates::{Conformer2D, Conformer3D};
use crate::embedding_values::CoordinateDimension;
use cosmolkit as ck;
use pyo3::prelude::*;

#[derive(FromPyObject, IntoPyObject)]
pub enum DimensionInput {
    Value(CoordinateDimension),
    Name(String),
}
impl DimensionInput {
    fn into_rust(self) -> PyResult<ck::CoordinateDimension> {
        match self {
            Self::Value(value) => Ok(value.into()),
            Self::Name(value) => CoordinateDimension::from_name(&value).map(Into::into),
        }
    }
}

#[pyclass(module = "fixed_projection._native", frozen, from_py_object)]
#[derive(Clone)]
pub struct CoordinateBlock {
    pub(crate) inner: ck::CoordinateBlock,
}
#[pymethods]
impl CoordinateBlock {
    #[new]
    #[pyo3(signature=(conformers_2d, conformers_3d, source_coordinate_dim, source_conformer_order))]
    fn new(
        conformers_2d: Vec<Conformer2D>,
        conformers_3d: Vec<Conformer3D>,
        source_coordinate_dim: Option<DimensionInput>,
        source_conformer_order: Option<Vec<DimensionInput>>,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: ck::CoordinateBlock {
                conformers_2d: conformers_2d.into_iter().map(|value| value.inner).collect(),
                conformers_3d: conformers_3d.into_iter().map(|value| value.inner).collect(),
                source_coordinate_dim: source_coordinate_dim
                    .map(DimensionInput::into_rust)
                    .transpose()?,
                source_conformer_order: source_conformer_order
                    .map(|values| values.into_iter().map(DimensionInput::into_rust).collect())
                    .transpose()?,
            },
        })
    }
    #[getter]
    fn conformers_2d(&self) -> Vec<Conformer2D> {
        self.inner
            .conformers_2d
            .iter()
            .cloned()
            .map(|inner| Conformer2D { inner })
            .collect()
    }
    #[getter]
    fn conformers_3d(&self) -> Vec<Conformer3D> {
        self.inner
            .conformers_3d
            .iter()
            .cloned()
            .map(|inner| Conformer3D { inner })
            .collect()
    }
    #[getter]
    fn source_coordinate_dim(&self) -> Option<CoordinateDimension> {
        self.inner.source_coordinate_dim.map(Into::into)
    }
    #[getter]
    fn source_conformer_order(&self) -> Option<Vec<CoordinateDimension>> {
        self.inner
            .source_conformer_order
            .as_ref()
            .map(|values| values.iter().copied().map(Into::into).collect())
    }
}

#[pyclass(module = "fixed_projection._native", frozen, from_py_object)]
#[derive(Clone)]
pub struct SdfRecordMetadata {
    pub(crate) inner: ck::SdfRecordMetadata,
}
#[pymethods]
impl SdfRecordMetadata {
    #[new]
    #[pyo3(signature=(index, byte_offset, byte_len, line_offset, line_len, title))]
    fn new(
        index: usize,
        byte_offset: u64,
        byte_len: u64,
        line_offset: usize,
        line_len: usize,
        title: Option<String>,
    ) -> Self {
        Self {
            inner: ck::SdfRecordMetadata {
                index,
                byte_offset,
                byte_len,
                line_offset,
                line_len,
                title,
            },
        }
    }
    fn index(&self) -> usize {
        self.inner.index()
    }
    fn byte_offset(&self) -> u64 {
        self.inner.byte_offset()
    }
    fn byte_len(&self) -> u64 {
        self.inner.byte_len()
    }
    fn byte_range(&self) -> (u64, u64) {
        self.inner.byte_range()
    }
    fn line_range(&self) -> (usize, usize) {
        self.inner.line_range()
    }
    fn title(&self) -> Option<&str> {
        self.inner.title()
    }
}
