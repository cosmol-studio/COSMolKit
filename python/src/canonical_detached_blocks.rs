//! Detached snapshots used by the facade's checked construction boundary.
//!
//! These wrappers retain the complete Rust values. Accessors copy readonly
//! fields into Python values; they never expose mutable live molecule storage.

use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct TopologyBlock {
    pub(crate) inner: ck::TopologyBlock,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TopologyBlock {
    #[getter]
    fn atoms(&self) -> Vec<crate::canonical_atom_bond::Atom> {
        self.inner
            .atoms
            .iter()
            .enumerate()
            .map(|(index, atom)| {
                crate::canonical_atom_bond::Atom::from_detached(
                    atom.clone(),
                    self.inner.adjacency.neighbors_of(index).len(),
                )
            })
            .collect()
    }
    #[getter]
    fn bonds(&self) -> Vec<crate::canonical_atom_bond::Bond> {
        self.inner
            .bonds
            .iter()
            .cloned()
            .map(|inner| crate::canonical_atom_bond::Bond { inner })
            .collect()
    }
    #[getter]
    fn substance_groups(&self) -> Vec<crate::canonical_group_values::SubstanceGroup> {
        self.inner
            .substance_groups
            .iter()
            .cloned()
            .map(|inner| crate::canonical_group_values::SubstanceGroup { inner })
            .collect()
    }
    #[getter]
    fn stereo_groups(&self) -> Vec<crate::canonical_group_values::StereoGroup> {
        self.inner
            .stereo_groups
            .iter()
            .cloned()
            .map(|inner| crate::canonical_group_values::StereoGroup { inner })
            .collect()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct Conformer2D {
    inner: ck::Conformer2D,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl Conformer2D {
    fn id(&self) -> usize {
        self.inner.id()
    }
    fn coordinates(&self) -> Vec<[f64; 2]> {
        self.inner.coordinates().to_vec()
    }
    fn props(&self, py: Python<'_>) -> PyResult<std::collections::BTreeMap<String, String>> {
        self.inner
            .props()
            .iter()
            .map(|(key, value)| {
                Ok((
                    crate::canonical_sdf::decode_source_text(py, key)?,
                    crate::canonical_sdf::decode_source_text(py, value)?,
                ))
            })
            .collect()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct CoordinateBlock {
    pub(crate) inner: ck::CoordinateBlock,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl CoordinateBlock {
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
    fn conformers_3d(&self) -> Vec<crate::mmff_binding::Conformer3D> {
        self.inner
            .conformers_3d
            .iter()
            .cloned()
            .map(|inner| crate::mmff_binding::Conformer3D { inner })
            .collect()
    }
    #[getter]
    fn source_coordinate_dim(&self) -> Option<crate::canonical_sdf::CoordinateDimension> {
        self.inner.source_coordinate_dim.map(Into::into)
    }
    #[getter]
    fn source_conformer_order(&self) -> Option<Vec<crate::canonical_sdf::CoordinateDimension>> {
        self.inner
            .source_conformer_order
            .as_ref()
            .map(|order| order.iter().copied().map(Into::into).collect())
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<TopologyBlock>()?;
    module.add_class::<CoordinateBlock>()?;
    module.add_class::<Conformer2D>()?;
    Ok(())
}
