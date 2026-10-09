//! Read-only transport of canonical atom/bond values and calculation results.
use ::cosmolkit as ck;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

pyo3::create_exception!(cosmolkit, CipDescriptorError, PyValueError);
pyo3::create_exception!(cosmolkit, PropertyValueError, PyValueError);
pyo3::create_exception!(cosmolkit, ValenceError, PyValueError);

pub(crate) fn cip_pyerr(py: Python<'_>, source: ck::CipDescriptorError) -> PyErr {
    let kind = match &source {
        ck::CipDescriptorError::InvalidNeighborOrder { .. } => "InvalidNeighborOrder",
        ck::CipDescriptorError::InvalidStoredDescriptor { .. } => "InvalidStoredDescriptor",
        ck::CipDescriptorError::Property(..) => "Property",
    };
    crate::canonical_values::annotate(
        py,
        CipDescriptorError::new_err(source.to_string()),
        "cip",
        kind,
        &source,
    )
}
pub(crate) fn property_pyerr(py: Python<'_>, source: ck::PropertyValueError) -> PyErr {
    let e = crate::canonical_values::annotate(
        py,
        PropertyValueError::new_err(source.to_string()),
        "property",
        "KindMismatch",
        &source,
    );
    let fields = || -> PyResult<()> {
        e.value(py).setattr(
            "_expected",
            Py::new(
                py,
                crate::canonical_property_values::PropertyValueKind::from(source.expected()),
            )?,
        )?;
        e.value(py).setattr(
            "_actual",
            Py::new(
                py,
                crate::canonical_property_values::PropertyValueKind::from(source.actual()),
            )?,
        )
    };
    if let Err(error) = fields() {
        return error;
    }
    e
}
pub(crate) fn valence_pyerr(py: Python<'_>, source: ck::ValenceError) -> PyErr {
    use ck::ValenceError as E;
    let kind = match &source {
        E::PiElectronExplicitValenceCacheNotInitialized { .. } => {
            "PiElectronExplicitValenceCacheNotInitialized"
        }
        E::PiElectronInvariant { .. } => "PiElectronInvariant",
        E::InvalidValence { .. } => "InvalidValence",
        E::InvalidTopology { .. } => "InvalidTopology",
        E::AtomOutOfRange { .. } => "AtomOutOfRange",
        E::AdjacencyAtomOutOfRange { .. } => "AdjacencyAtomOutOfRange",
        E::AdjacencyBondOutOfRange { .. } => "AdjacencyBondOutOfRange",
        E::AdjacencyEndpointMismatch { .. } => "AdjacencyEndpointMismatch",
        E::InvalidExplicitValenceInput { .. } => "InvalidExplicitValenceInput",
        E::PeriodicTableLookup { .. } => "PeriodicTableLookup",
        E::ExplicitValenceCacheNotInitialized { .. } => "ExplicitValenceCacheNotInitialized",
        E::ImplicitValenceCacheNotInitialized { .. } => "ImplicitValenceCacheNotInitialized",
        E::HydrogenCountOverflow { .. } => "HydrogenCountOverflow",
        E::BadBondType { .. } => "BadBondType",
    };
    let error = crate::canonical_values::annotate(
        py,
        ValenceError::new_err(source.to_string()),
        "valence",
        kind,
        &source,
    );
    let attributes = || -> PyResult<()> {
        let value = error.value(py);
        match &source {
            E::PiElectronExplicitValenceCacheNotInitialized { atom }
            | E::ExplicitValenceCacheNotInitialized { atom }
            | E::ImplicitValenceCacheNotInitialized { atom } => {
                value.setattr("atom", atom.index())?
            }
            E::PiElectronInvariant {
                atom,
                explicit_valence,
                physical_bonds,
            } => {
                value.setattr("atom", atom.index())?;
                value.setattr("explicit_valence", *explicit_valence)?;
                value.setattr("physical_bonds", *physical_bonds)?;
            }
            E::InvalidValence {
                atom,
                atomic_number,
                formal_charge,
                phase,
                calculated,
                reason,
                message,
            } => {
                value.setattr("atom", atom.index())?;
                value.setattr("atomic_number", *atomic_number)?;
                value.setattr("formal_charge", *formal_charge)?;
                value.setattr("phase", format!("{phase:?}"))?;
                value.setattr("calculated", *calculated)?;
                value.setattr("reason", *reason)?;
                value.setattr("message", message.as_str())?;
            }
            E::AtomOutOfRange { atom, atom_count } => {
                value.setattr("atom", atom.index())?;
                value.setattr("atom_count", *atom_count)?;
            }
            E::AdjacencyAtomOutOfRange {
                atom,
                neighbor_atom,
                atom_count,
            } => {
                value.setattr("atom", atom.index())?;
                value.setattr("neighbor_atom", *neighbor_atom)?;
                value.setattr("atom_count", *atom_count)?;
            }
            E::AdjacencyBondOutOfRange {
                atom,
                bond,
                bond_count,
            } => {
                value.setattr("atom", atom.index())?;
                value.setattr("bond", bond.index())?;
                value.setattr("bond_count", *bond_count)?;
            }
            E::AdjacencyEndpointMismatch {
                atom,
                neighbor_atom,
                bond,
                begin,
                end,
            } => {
                value.setattr("atom", atom.index())?;
                value.setattr("neighbor_atom", *neighbor_atom)?;
                value.setattr("bond", bond.index())?;
                value.setattr("begin", begin.index())?;
                value.setattr("end", end.index())?;
            }
            E::InvalidExplicitValenceInput { atom, value: input } => {
                value.setattr("atom", atom.index())?;
                value.setattr("value", *input)?;
            }
            E::PeriodicTableLookup {
                atomic_number,
                field,
            } => {
                value.setattr("atomic_number", *atomic_number)?;
                value.setattr("field", *field)?;
            }
            E::HydrogenCountOverflow {
                atom,
                explicit,
                implicit,
                neighbor_hydrogens,
            } => {
                value.setattr("atom", atom.index())?;
                value.setattr("explicit", *explicit)?;
                value.setattr("implicit", *implicit)?;
                value.setattr("neighbor_hydrogens", *neighbor_hydrogens)?;
            }
            E::BadBondType { bond, order } => {
                value.setattr("bond", bond.map(|id| id.index()))?;
                value.setattr("order", enum_member(py, "BondOrder", *order as i64)?)?;
            }
            E::InvalidTopology { .. } => {}
        }
        Ok(())
    };
    match attributes() {
        Ok(()) => error,
        Err(error) => error,
    }
}
pub(crate) fn enum_member<'py>(
    py: Python<'py>,
    name: &str,
    code: i64,
) -> PyResult<Bound<'py, PyAny>> {
    py.import("cosmolkit")?.getattr(name)?.call1((code,))
}

/// Dynamic IntEnum inputs use the already registered source vocabulary.
/// Retain existing explicit integer-code calls; never guess an unknown name.
pub(crate) fn enum_code(value: &Bound<'_, PyAny>, name: &str) -> PyResult<i64> {
    use pyo3::exceptions::PyTypeError;
    use pyo3::types::{PyInt, PyString};
    let class = value.py().import("cosmolkit")?.getattr(name)?;
    if value.is_instance_of::<PyString>() {
        let text: String = value.extract()?;
        let members = class.getattr("__members__")?;
        for item in members.call_method0("items")?.try_iter()? {
            let (member_name, member): (String, Bound<'_, PyAny>) = item?.extract()?;
            if text == member_name.to_ascii_lowercase() {
                return member.extract();
            }
        }
        return Err(PyValueError::new_err(format!(
            "invalid {name} string: {text:?}"
        )));
    }
    if value.is_instance(&class)? || value.is_instance_of::<PyInt>() {
        return value.extract();
    }
    Err(PyTypeError::new_err(format!(
        "expected {name}, str or integer code"
    )))
}
fn descriptor_member<'py>(
    py: Python<'py>,
    descriptor: Option<ck::CipDescriptor>,
) -> PyResult<Option<Bound<'py, PyAny>>> {
    descriptor
        .map(|x| {
            py.import("cosmolkit")?
                .getattr("CipDescriptor")?
                .call1((x.as_str(),))
        })
        .transpose()
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AtomMetadata {
    pub(crate) inner: ck::AtomMetadata,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomMetadata {
    fn degree(&self) -> usize {
        self.inner.degree
    }
    fn explicit_valence(&self) -> i32 {
        self.inner.explicit_valence
    }
    fn implicit_hydrogens(&self) -> i32 {
        self.inner.implicit_hydrogens
    }
    fn total_hydrogens(&self) -> i32 {
        self.inner.total_hydrogens
    }
    fn total_valence(&self) -> i32 {
        self.inner.total_valence
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct Atom {
    pub(crate) inner: ck::Atom,
    pub(crate) degree: usize,
    // An actual owner error is retained; raw atom fields remain readable.
    pub(crate) metadata: Result<ck::AtomMetadata, ck::ValenceError>,
}

impl Atom {
    /// Detached Atom values have no runtime valence assignment. Preserve that
    /// absence explicitly; reading raw construction state must not run chemistry.
    pub(crate) fn from_detached(inner: ck::Atom, degree: usize) -> Self {
        let metadata =
            Err(ck::ValenceError::ExplicitValenceCacheNotInitialized { atom: inner.id() });
        Self {
            inner,
            degree,
            metadata,
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl Atom {
    fn template_attachment_order(
        &self,
    ) -> Option<crate::canonical_group_values::TemplateAttachmentOrder> {
        self.inner
            .template_attachment_order()
            .cloned()
            .map(|inner| crate::canonical_group_values::TemplateAttachmentOrder { inner })
    }
    fn id(&self) -> usize {
        self.inner.id().index()
    }
    fn degree(&self) -> usize {
        self.degree
    }
    fn element(&self) -> crate::canonical_element_metadata::Element {
        crate::canonical_element_metadata::Element {
            inner: self.inner.element(),
        }
    }
    #[gen_stub(override_return_type(type_repr = "ChiralTag"))]
    fn chiral_tag<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "ChiralTag", self.inner.chiral_tag_code())
    }
    #[gen_stub(override_return_type(type_repr="typing.Optional[CipDescriptor]",imports=("typing")))]
    fn cip_descriptor<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyAny>>> {
        descriptor_member(
            py,
            self.inner.cip_descriptor().map_err(|e| cip_pyerr(py, e))?,
        )
    }
    fn cip_neighbor_order(&self, py: Python<'_>) -> PyResult<Option<Vec<u32>>> {
        self.inner
            .cip_neighbor_order()
            .map_err(|e| cip_pyerr(py, e))
    }
    fn cip_rank(&self, py: Python<'_>) -> PyResult<Option<u32>> {
        self.inner.cip_rank().map_err(|e| property_pyerr(py, e))
    }
    fn atomic_number(&self) -> u8 {
        self.inner.atomic_number()
    }
    fn formal_charge(&self) -> i8 {
        self.inner.formal_charge()
    }
    fn explicit_hydrogens(&self) -> u8 {
        self.inner.explicit_hydrogens()
    }
    fn chiral_tag_code(&self) -> i64 {
        self.inner.chiral_tag_code()
    }
    fn chiral_tag_name(&self) -> &'static str {
        self.inner.chiral_tag_name()
    }
    fn isotope(&self) -> Option<u16> {
        self.inner.isotope()
    }
    fn atom_map(&self) -> Option<u32> {
        self.inner.atom_map()
    }
    fn is_aromatic(&self) -> bool {
        self.inner.is_aromatic()
    }
    fn no_implicit(&self) -> bool {
        self.inner.no_implicit()
    }
    #[gen_stub(override_return_type(type_repr = "Hybridization"))]
    fn hybridization<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "Hybridization", self.inner.hybridization().rdkit_code())
    }
    fn radical_electrons(&self) -> u8 {
        self.inner.radical_electrons()
    }
    fn explicit_valence(&self, py: Python<'_>) -> PyResult<i32> {
        self.metadata
            .as_ref()
            .map(|m| m.explicit_valence)
            .map_err(|e| valence_pyerr(py, e.clone()))
    }
    fn implicit_hydrogens(&self, py: Python<'_>) -> PyResult<i32> {
        self.metadata
            .as_ref()
            .map(|m| m.implicit_hydrogens)
            .map_err(|e| valence_pyerr(py, e.clone()))
    }
    fn total_hydrogens(&self, py: Python<'_>) -> PyResult<i32> {
        self.metadata
            .as_ref()
            .map(|m| m.total_hydrogens)
            .map_err(|e| valence_pyerr(py, e.clone()))
    }
    fn total_valence(&self, py: Python<'_>) -> PyResult<i32> {
        self.metadata
            .as_ref()
            .map(|m| m.total_valence)
            .map_err(|e| valence_pyerr(py, e.clone()))
    }
    fn __repr__(&self) -> String {
        format!(
            "Atom(id={}, atomic_number={}, formal_charge={}, chiral_tag='{}', isotope={}, is_aromatic={}, degree={})",
            self.id(),
            self.atomic_number(),
            self.formal_charge(),
            self.chiral_tag_name(),
            self.isotope()
                .map(|x| x.to_string())
                .unwrap_or_else(|| "None".to_owned()),
            self.is_aromatic(),
            self.degree()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct Bond {
    pub(crate) inner: ck::Bond,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl Bond {
    fn id(&self) -> usize {
        self.inner.id().index()
    }
    fn begin(&self) -> usize {
        self.inner.begin().index()
    }
    fn end(&self) -> usize {
        self.inner.end().index()
    }
    fn stereo_atoms(&self) -> Option<[usize; 2]> {
        self.inner.stereo_atoms().map(|x| x.map(ck::AtomId::index))
    }
    #[gen_stub(override_return_type(type_repr="typing.Optional[CipDescriptor]",imports=("typing")))]
    fn cip_descriptor<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyAny>>> {
        descriptor_member(
            py,
            self.inner.cip_descriptor().map_err(|e| cip_pyerr(py, e))?,
        )
    }
    fn cip_neighbor_order(&self, py: Python<'_>) -> PyResult<Option<Vec<u32>>> {
        self.inner
            .cip_neighbor_order()
            .map_err(|e| cip_pyerr(py, e))
    }
    #[gen_stub(override_return_type(type_repr = "BondOrder"))]
    fn order<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "BondOrder", self.inner.order_code())
    }
    #[gen_stub(override_return_type(type_repr = "BondDirection"))]
    fn direction<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "BondDirection", self.inner.direction_code())
    }
    #[gen_stub(override_return_type(type_repr = "BondStereo"))]
    fn stereo<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        enum_member(py, "BondStereo", self.inner.stereo_code())
    }
    fn order_code(&self) -> i64 {
        self.inner.order_code()
    }
    fn order_name(&self) -> &'static str {
        self.inner.order_name()
    }
    fn direction_code(&self) -> i64 {
        self.inner.direction_code()
    }
    fn direction_name(&self) -> &'static str {
        self.inner.direction_name()
    }
    fn stereo_code(&self) -> i64 {
        self.inner.stereo_code()
    }
    fn stereo_name(&self) -> &'static str {
        self.inner.stereo_name()
    }
    fn is_aromatic(&self) -> bool {
        self.inner.is_aromatic()
    }
    fn is_conjugated(&self) -> bool {
        self.inner.is_conjugated()
    }
    fn __repr__(&self) -> String {
        format!(
            "Bond(id={}, begin={}, end={}, order='{}', direction='{}', stereo='{}')",
            self.id(),
            self.begin(),
            self.end(),
            self.order_name(),
            self.direction_name(),
            self.stereo_name()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct CipLabelOptions {
    pub(crate) inner: ck::CipLabelOptions,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl CipLabelOptions {
    #[new]
    #[pyo3(signature=(*,atoms=None,bonds=None,max_recursive_iterations=0))]
    fn new(
        atoms: Option<Vec<usize>>,
        bonds: Option<Vec<usize>>,
        max_recursive_iterations: u32,
    ) -> Self {
        let mut inner =
            ck::CipLabelOptions::default().with_max_recursive_iterations(max_recursive_iterations);
        if let Some(atoms) = atoms {
            inner = inner.with_atoms(atoms.into_iter().map(ck::AtomId::new));
        }
        if let Some(bonds) = bonds {
            inner = inner.with_bonds(bonds.into_iter().map(ck::BondId::new));
        }
        Self { inner }
    }
    #[getter]
    fn atoms(&self) -> Option<Vec<usize>> {
        self.inner
            .atoms()
            .map(|x| x.iter().map(|x| x.index()).collect())
    }
    #[getter]
    fn bonds(&self) -> Option<Vec<usize>> {
        self.inner
            .bonds()
            .map(|x| x.iter().map(|x| x.index()).collect())
    }
    #[getter]
    fn max_recursive_iterations(&self) -> u32 {
        self.inner.max_recursive_iterations()
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    crate::canonical_error_accessors::attach(
        module.py().get_type::<PropertyValueError>().as_any(),
        &[("expected", "_expected"), ("actual", "_actual")],
    )?;
    module.add_class::<Atom>()?;
    module.add_class::<Bond>()?;
    module.add_class::<AtomMetadata>()?;
    module.add_class::<CipLabelOptions>()?;
    module.add(
        "CipDescriptorError",
        module.py().get_type::<CipDescriptorError>(),
    )?;
    module.add(
        "PropertyValueError",
        module.py().get_type::<PropertyValueError>(),
    )?;
    module.add("ValenceError", module.py().get_type::<ValenceError>())?;
    let enum_module = module.py().import("enum")?;
    let members = PyDict::new(module.py());
    for code in 0..=ck::Hybridization::Other.rdkit_code() {
        let value = ck::Hybridization::from_rdkit_code(code).ok_or_else(|| {
            PyRuntimeError::new_err(format!(
                "canonical Hybridization lacks declared code {code}"
            ))
        })?;
        members.set_item(value.rdkit_name(), value.rdkit_code())?;
    }
    let kwargs = PyDict::new(module.py());
    kwargs.set_item("module", "cosmolkit")?;
    module.add(
        "Hybridization",
        enum_module
            .getattr("IntEnum")?
            .call(("Hybridization", members), Some(&kwargs))?,
    )?;
    let members = PyDict::new(module.py());
    for code in 0..=ck::BondOrder::Zero.rdkit_code() {
        let value = ck::BondOrder::from_rdkit_code(code).ok_or_else(|| {
            PyRuntimeError::new_err(format!("canonical BondOrder lacks declared code {code}"))
        })?;
        members.set_item(value.rdkit_name(), value.rdkit_code())?;
    }
    let kwargs = PyDict::new(module.py());
    kwargs.set_item("module", "cosmolkit")?;
    module.add(
        "BondOrder",
        enum_module
            .getattr("IntEnum")?
            .call(("BondOrder", members), Some(&kwargs))?,
    )?;
    let members = PyDict::new(module.py());
    for code in 0..=ck::ChiralTag::Octahedral.rdkit_code() {
        let value = ck::ChiralTag::from_rdkit_code(code).ok_or_else(|| {
            PyRuntimeError::new_err(format!("canonical ChiralTag lacks declared code {code}"))
        })?;
        members.set_item(value.rdkit_name(), value.rdkit_code())?;
    }
    let kwargs = PyDict::new(module.py());
    kwargs.set_item("module", "cosmolkit")?;
    module.add(
        "ChiralTag",
        enum_module
            .getattr("IntEnum")?
            .call(("ChiralTag", members), Some(&kwargs))?,
    )?;
    let members = PyDict::new(module.py());
    for code in 0..=ck::BondDirection::Unknown.rdkit_code() {
        let value = ck::BondDirection::from_rdkit_code(code).ok_or_else(|| {
            PyRuntimeError::new_err(format!(
                "canonical BondDirection lacks declared code {code}"
            ))
        })?;
        members.set_item(value.rdkit_name(), value.rdkit_code())?;
    }
    let kwargs = PyDict::new(module.py());
    kwargs.set_item("module", "cosmolkit")?;
    module.add(
        "BondDirection",
        enum_module
            .getattr("IntEnum")?
            .call(("BondDirection", members), Some(&kwargs))?,
    )?;
    let members = PyDict::new(module.py());
    for code in 0..=ck::BondStereo::AtropCcw.rdkit_code() {
        let value = ck::BondStereo::from_rdkit_code(code).ok_or_else(|| {
            PyRuntimeError::new_err(format!("canonical BondStereo lacks declared code {code}"))
        })?;
        members.set_item(value.rdkit_name(), value.rdkit_code())?;
    }
    let kwargs = PyDict::new(module.py());
    kwargs.set_item("module", "cosmolkit")?;
    module.add(
        "BondStereo",
        enum_module
            .getattr("IntEnum")?
            .call(("BondStereo", members), Some(&kwargs))?,
    )?;
    let members = PyDict::new(module.py());
    for value in [
        ck::CipDescriptor::R,
        ck::CipDescriptor::S,
        ck::CipDescriptor::LowerR,
        ck::CipDescriptor::LowerS,
        ck::CipDescriptor::E,
        ck::CipDescriptor::Z,
        ck::CipDescriptor::LowerE,
        ck::CipDescriptor::LowerZ,
        ck::CipDescriptor::M,
        ck::CipDescriptor::P,
        ck::CipDescriptor::LowerM,
        ck::CipDescriptor::LowerP,
    ] {
        members.set_item(value.as_str(), value.as_str())?;
    }
    let kwargs = PyDict::new(module.py());
    kwargs.set_item("module", "cosmolkit")?;
    kwargs.set_item("type", module.py().get_type::<pyo3::types::PyString>())?;
    module.add(
        "CipDescriptor",
        enum_module
            .getattr("Enum")?
            .call(("CipDescriptor", members), Some(&kwargs))?,
    )?;
    Ok(())
}
