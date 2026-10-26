use crate::bio_values::{BioCrystalCell, BioTransform};
use crate::chemical_values::{EnumValue, enum_value};
use crate::coordinates;
use cosmolkit as ck;
use pyo3::prelude::*;

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Cartesian atom coordinates stored as three-component rows; positions returns an independent NumPy array.
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
    /// Independent float64 NumPy array with shape (num_atoms, 3), in stored atom order.
    fn positions<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f64>> {
        coordinates::array(py, self.inner.positions())
    }
    #[getter]
    /// Number of stored Cartesian coordinate rows.
    fn len(&self) -> usize {
        self.inner.len()
    }
    #[getter]
    /// Whether there are no coordinate rows.
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Unit-cell geometry, coordinate transforms and retained crystallographic symmetry metadata.
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
    /// International space-group number inferred from the stored designation, or None if unresolved.
    fn space_group_number(&self) -> Option<i32> {
        self.inner.space_group_number()
    }
    #[getter]
    /// Whether the native cell information represents a crystallographic unit cell.
    fn is_crystal(&self) -> bool {
        self.inner.is_crystal()
    }
    #[getter]
    /// Independent copy of the stored unit-cell parameters.
    fn cell(&self, py: Python<'_>) -> PyResult<BioCrystalCell> {
        BioCrystalCell::from_rust(py, self.inner.cell())
    }
    #[getter]
    /// Transform from fractional to Cartesian coordinates.
    fn orthogonal(&self) -> BioTransform {
        BioTransform {
            inner: *self.inner.orthogonal(),
        }
    }
    #[getter]
    /// Transform from Cartesian to fractional coordinates.
    fn fractional(&self) -> BioTransform {
        BioTransform {
            inner: *self.inner.fractional(),
        }
    }
    #[getter]
    /// Unit-cell volume computed from the stored cell geometry.
    fn volume(&self) -> f64 {
        self.inner.volume()
    }
    #[getter]
    /// Three reciprocal-cell axis lengths, in a*, b*, c* order.
    fn reciprocal_lengths(&self) -> Vec<f64> {
        self.inner.reciprocal_lengths().to_vec()
    }
    #[getter]
    /// Cosines of the three reciprocal-cell angles.
    fn reciprocal_cosines(&self) -> Vec<f64> {
        self.inner.reciprocal_cosines().to_vec()
    }
    #[getter]
    /// Whether the coordinate conversion matrices were explicitly supplied.
    fn explicit_matrices(&self) -> bool {
        self.inner.explicit_matrices()
    }
    #[getter]
    /// Stored number of crystallographic symmetry images.
    fn cs_count(&self) -> i16 {
        self.inner.cs_count()
    }
    #[getter]
    /// Independent list of retained symmetry-image transforms.
    fn symmetry_images(&self) -> Vec<BioTransform> {
        self.inner
            .symmetry_images()
            .iter()
            .copied()
            .map(|inner| BioTransform { inner })
            .collect()
    }
    #[getter]
    /// Hermann-Mauguin space-group designation, or None if absent.
    fn space_group_hm(&self) -> Option<&str> {
        self.inner.space_group_hm()
    }
    #[getter]
    /// Source PDB Z value as text, or None if absent.
    fn z_pdb(&self) -> Option<&str> {
        self.inner.z_pdb()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Non-crystallographic symmetry operator with its source identifier and transform.
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
    /// Identifier of this non-crystallographic symmetry operator in the source structure.
    fn id(&self) -> &str {
        &self.inner.id
    }
    #[getter]
    /// Whether this NCS operator is marked as already represented by the input coordinates.
    fn given(&self) -> bool {
        self.inner.given
    }
    #[getter]
    /// Affine transform associated with this symmetry or assembly operator.
    fn transform(&self) -> BioTransform {
        BioTransform {
            inner: self.inner.transform,
        }
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Named affine operator used to generate a biological assembly.
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
    /// Source identifier or name of this assembly or operator.
    fn name(&self) -> Option<&str> {
        self.inner.name.as_deref()
    }
    #[getter]
    /// Source classification of the assembly operator, or None if absent.
    fn operator_type(&self) -> Option<&str> {
        self.inner.operator_type.as_deref()
    }
    #[getter]
    /// Affine transform associated with this symmetry or assembly operator.
    fn transform(&self) -> BioTransform {
        BioTransform {
            inner: self.inner.transform,
        }
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Chain and subchain selections together with the operators that generate an assembly.
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
    /// Independent list of chain identifiers to which the generator applies.
    fn chains(&self) -> Vec<String> {
        self.inner.chains.clone()
    }
    #[getter]
    /// Independent list of subchain identifiers to which the generator applies.
    fn subchains(&self) -> Vec<String> {
        self.inner.subchains.clone()
    }
    #[getter]
    /// Independent list of assembly operators applied by this generator.
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
    /// Biological assembly definition with source provenance, reported metrics and generation rules.
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
    /// Source identifier or name of this assembly or operator.
    fn name(&self) -> &str {
        &self.inner.name
    }
    #[getter]
    /// Whether the source marks this assembly as author-determined.
    fn author_determined(&self) -> bool {
        self.inner.author_determined
    }
    #[getter]
    /// Whether the source marks this assembly as software-determined.
    fn software_determined(&self) -> bool {
        self.inner.software_determined
    }
    #[getter]
    /// Stored special assembly classification.
    fn special_kind(&self) -> BioAssemblySpecialKind {
        self.inner.special_kind.into()
    }
    #[getter]
    /// Reported number of components in the oligomeric assembly.
    fn oligomeric_count(&self) -> i32 {
        self.inner.oligomeric_count
    }
    #[getter]
    /// Source description of the oligomeric assembly.
    fn oligomeric_details(&self) -> &str {
        &self.inner.oligomeric_details
    }
    #[getter]
    /// Name of the software used to determine this assembly.
    fn software_name(&self) -> &str {
        &self.inner.software_name
    }
    #[getter]
    /// Buried surface-area value reported for this assembly.
    fn buried_surface_area(&self) -> f64 {
        self.inner.buried_surface_area
    }
    #[getter]
    /// Surface-area value reported for this assembly.
    fn surface_area(&self) -> f64 {
        self.inner.surface_area
    }
    #[getter]
    /// Solvent free-energy change reported for this assembly.
    fn solvent_free_energy_change(&self) -> f64 {
        self.inner.solvent_free_energy_change
    }
    #[getter]
    /// Independent list of generators defining which chains and transforms build the assembly.
    fn generators(&self) -> Vec<BioAssemblyGenerator> {
        self.inner
            .generators
            .iter()
            .cloned()
            .map(|inner| BioAssemblyGenerator { inner })
            .collect()
    }
}
