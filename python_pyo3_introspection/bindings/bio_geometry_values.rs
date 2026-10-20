use crate::bio_values::{BioCrystalCell, BioTransform};
use crate::chemical_values::{EnumValue, enum_value};
use crate::coordinates;
use cosmolkit as ck;
use pyo3::prelude::*;

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct BioCoordinateBlock {
        pub(crate) inner: ck::BioCoordinateBlock,
    }
}
#[pymethods]
impl BioCoordinateBlock {
    #[new]
    fn new(positions: Vec<[f64; 3]>) -> Self {
        Self {
            inner: ck::BioCoordinateBlock::new(positions),
        }
    }
    #[doc = r#"{"role":"return_type","imports":["numpy"],"type":"numpy.ndarray[tuple[int, int], numpy.dtype[numpy.float64]]"}"#]
    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f64>> {
        coordinates::array(py, self.inner.positions())
    }
    #[getter]
    fn len(&self) -> usize {
        self.inner.len()
    }
    #[getter]
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct BioCrystalInfo {
        pub(crate) inner: ck::BioCrystalInfo,
    }
}
#[pymethods]
impl BioCrystalInfo {
    #[new]
    #[pyo3(signature=(cell, space_group_hm, z_pdb, orthogonal, fractional, explicit_matrices, cs_count, symmetry_images))]
    fn new(
        py: Python<'_>,
        cell: &BioCrystalCell,
        space_group_hm: Option<String>,
        z_pdb: Option<String>,
        orthogonal: BioTransform,
        fractional: BioTransform,
        explicit_matrices: bool,
        cs_count: i16,
        symmetry_images: Vec<BioTransform>,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: ck::BioCrystalInfo::new(
                cell.to_rust(py)?,
                space_group_hm,
                z_pdb,
                orthogonal.inner,
                fractional.inner,
                explicit_matrices,
                cs_count,
                symmetry_images
                    .into_iter()
                    .map(|value| value.inner)
                    .collect(),
            ),
        })
    }
    #[getter]
    fn space_group_number(&self) -> Option<i32> {
        self.inner.space_group_number()
    }
    #[getter]
    fn is_crystal(&self) -> bool {
        self.inner.is_crystal()
    }
    #[getter]
    fn cell(&self, py: Python<'_>) -> PyResult<BioCrystalCell> {
        BioCrystalCell::from_rust(py, self.inner.cell())
    }
    #[getter]
    fn orthogonal(&self) -> BioTransform {
        BioTransform {
            inner: *self.inner.orthogonal(),
        }
    }
    #[getter]
    fn fractional(&self) -> BioTransform {
        BioTransform {
            inner: *self.inner.fractional(),
        }
    }
    #[getter]
    fn volume(&self) -> f64 {
        self.inner.volume()
    }
    #[getter]
    fn reciprocal_lengths(&self) -> Vec<f64> {
        self.inner.reciprocal_lengths().to_vec()
    }
    #[getter]
    fn reciprocal_cosines(&self) -> Vec<f64> {
        self.inner.reciprocal_cosines().to_vec()
    }
    #[getter]
    fn explicit_matrices(&self) -> bool {
        self.inner.explicit_matrices()
    }
    #[getter]
    fn cs_count(&self) -> i16 {
        self.inner.cs_count()
    }
    #[getter]
    fn symmetry_images(&self) -> Vec<BioTransform> {
        self.inner
            .symmetry_images()
            .iter()
            .copied()
            .map(|inner| BioTransform { inner })
            .collect()
    }
    #[getter]
    fn space_group_hm(&self) -> Option<&str> {
        self.inner.space_group_hm()
    }
    #[getter]
    fn z_pdb(&self) -> Option<&str> {
        self.inner.z_pdb()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct BioNcsOperator {
        pub(crate) inner: ck::BioNcsOperator,
    }
}
#[pymethods]
impl BioNcsOperator {
    #[new]
    fn new(id: String, given: bool, transform: BioTransform) -> Self {
        Self {
            inner: ck::BioNcsOperator::new(id, given, transform.inner),
        }
    }
    #[getter]
    fn id(&self) -> &str {
        &self.inner.id
    }
    #[getter]
    fn given(&self) -> bool {
        self.inner.given
    }
    #[getter]
    fn transform(&self) -> BioTransform {
        BioTransform {
            inner: self.inner.transform,
        }
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct BioAssemblyOperator {
        pub(crate) inner: ck::BioAssemblyOperator,
    }
}
#[pymethods]
impl BioAssemblyOperator {
    #[new]
    #[pyo3(signature=(name, operator_type, transform))]
    fn new(name: Option<String>, operator_type: Option<String>, transform: BioTransform) -> Self {
        Self {
            inner: ck::BioAssemblyOperator::new(name, operator_type, transform.inner),
        }
    }
    #[getter]
    fn name(&self) -> Option<&str> {
        self.inner.name.as_deref()
    }
    #[getter]
    fn operator_type(&self) -> Option<&str> {
        self.inner.operator_type.as_deref()
    }
    #[getter]
    fn transform(&self) -> BioTransform {
        BioTransform {
            inner: self.inner.transform,
        }
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct BioAssemblyGenerator {
        pub(crate) inner: ck::BioAssemblyGenerator,
    }
}
#[pymethods]
impl BioAssemblyGenerator {
    #[new]
    #[pyo3(signature=(*, chains=ck::BioAssemblyGenerator::default().chains, subchains=ck::BioAssemblyGenerator::default().subchains, operators=None))]
    fn new(
        chains: Vec<String>,
        subchains: Vec<String>,
        operators: Option<Vec<BioAssemblyOperator>>,
    ) -> Self {
        Self {
            inner: ck::BioAssemblyGenerator::new(
                chains,
                subchains,
                operators
                    .map(|values| values.into_iter().map(|value| value.inner).collect())
                    .unwrap_or_else(|| ck::BioAssemblyGenerator::default().operators),
            ),
        }
    }
    #[getter]
    fn chains(&self) -> Vec<String> {
        self.inner.chains.clone()
    }
    #[getter]
    fn subchains(&self) -> Vec<String> {
        self.inner.subchains.clone()
    }
    #[getter]
    fn operators(&self) -> Vec<BioAssemblyOperator> {
        self.inner
            .operators
            .iter()
            .cloned()
            .map(|inner| BioAssemblyOperator { inner })
            .collect()
    }
}

enum_value!(BioAssemblySpecialKind, { NotApplicable, CompleteIcosahedral, RepresentativeHelical, CompletePoint });
#[derive(FromPyObject, IntoPyObject)]
pub enum BioAssemblySpecialKindInput {
    Enum(BioAssemblySpecialKind),
    Name(String),
}
impl BioAssemblySpecialKindInput {
    fn into_rust(self) -> PyResult<ck::BioAssemblySpecialKind> {
        match self {
            Self::Enum(value) => Ok(value.into()),
            Self::Name(name) => BioAssemblySpecialKind::from_name(&name).map(Into::into),
        }
    }
}
crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct BioAssembly {
        pub(crate) inner: ck::BioAssembly,
    }
}
#[pymethods]
impl BioAssembly {
    #[new]
    fn new(
        name: String,
        author_determined: bool,
        software_determined: bool,
        special_kind: BioAssemblySpecialKindInput,
        oligomeric_count: i32,
        oligomeric_details: String,
        software_name: String,
        buried_surface_area: f64,
        surface_area: f64,
        solvent_free_energy_change: f64,
        generators: Vec<BioAssemblyGenerator>,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: ck::BioAssembly::new(
                name,
                author_determined,
                software_determined,
                special_kind.into_rust()?,
                oligomeric_count,
                oligomeric_details,
                software_name,
                buried_surface_area,
                surface_area,
                solvent_free_energy_change,
                generators.into_iter().map(|value| value.inner).collect(),
            ),
        })
    }
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }
    #[getter]
    fn author_determined(&self) -> bool {
        self.inner.author_determined
    }
    #[getter]
    fn software_determined(&self) -> bool {
        self.inner.software_determined
    }
    #[getter]
    fn special_kind(&self) -> BioAssemblySpecialKind {
        self.inner.special_kind.into()
    }
    #[getter]
    fn oligomeric_count(&self) -> i32 {
        self.inner.oligomeric_count
    }
    #[getter]
    fn oligomeric_details(&self) -> &str {
        &self.inner.oligomeric_details
    }
    #[getter]
    fn software_name(&self) -> &str {
        &self.inner.software_name
    }
    #[getter]
    fn buried_surface_area(&self) -> f64 {
        self.inner.buried_surface_area
    }
    #[getter]
    fn surface_area(&self) -> f64 {
        self.inner.surface_area
    }
    #[getter]
    fn solvent_free_energy_change(&self) -> f64 {
        self.inner.solvent_free_energy_change
    }
    #[getter]
    fn generators(&self) -> Vec<BioAssemblyGenerator> {
        self.inner
            .generators
            .iter()
            .cloned()
            .map(|inner| BioAssemblyGenerator { inner })
            .collect()
    }
}
