//! Direct projections of canonical potential-stereo parameters and result rows.
use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

pyo3::create_exception!(
    cosmolkit,
    PotentialStereoError,
    pyo3::exceptions::PyValueError,
    "Potential stereocenters could not be identified for the supplied graph."
);

/// Preserve the canonical error variant and its scalar carrier fields. This
/// projection neither retries perception nor turns a present malformed value
/// into an absent property.
pub(crate) fn error_pyerr(py: Python<'_>, source: &ck::PotentialStereoError) -> PyErr {
    use ck::PotentialStereoError as E;
    let kind = match source {
        E::InvalidPropertyKind { .. } => "InvalidPropertyKind",
        E::InvalidRingStereoReference { .. } => "InvalidRingStereoReference",
        E::EmptyRingStereoReferences { .. } => "EmptyRingStereoReferences",
        E::Ring(_) => "Ring",
        E::InvalidTopology(_) => "InvalidTopology",
        E::InvalidValence { .. } => "InvalidValence",
        E::InvalidValenceValue { .. } => "InvalidValenceValue",
        E::AtomOutOfRange { .. } => "AtomOutOfRange",
        E::InvalidRingInfo { .. } => "InvalidRingInfo",
        E::InvalidAtomDegree { .. } => "InvalidAtomDegree",
        E::InvalidBondDegree { .. } => "InvalidBondDegree",
        E::InvalidStereoReferences { .. } => "InvalidStereoReferences",
        E::UnsupportedBondOrder { .. } => "UnsupportedBondOrder",
        E::InvalidChiralPermutation { .. } => "InvalidChiralPermutation",
        E::AtropisomerDependencyUnavailable { .. } => "AtropisomerDependencyUnavailable",
        E::RefinementDidNotConverge { .. } => "RefinementDidNotConverge",
        E::Valence(_) => "Valence",
        E::StereoOrder(_) => "StereoOrder",
        E::DoubleStereo(_) => "DoubleStereo",
        E::BondValue(_) => "BondValue",
    };
    let error = crate::canonical_values::annotate(
        py,
        PotentialStereoError::new_err(source.to_string()),
        "stereo",
        kind,
        source,
    );
    let attributes = || -> PyResult<()> {
        let value = error.value(py);
        match source {
            E::InvalidPropertyKind {
                atom,
                property,
                kind,
            } => {
                value.setattr("atom", atom.index())?;
                value.setattr("property", *property)?;
                value.setattr("property_kind", format!("{kind:?}"))?;
            }
            E::InvalidRingStereoReference {
                atom,
                value: reference,
            } => {
                value.setattr("atom", atom.index())?;
                value.setattr("value", *reference)?;
            }
            E::EmptyRingStereoReferences { atom } => value.setattr("atom", atom.index())?,
            E::InvalidValence {
                field,
                actual,
                atom_count,
            } => {
                value.setattr("field", *field)?;
                value.setattr("actual", *actual)?;
                value.setattr("atom_count", *atom_count)?;
            }
            E::InvalidValenceValue {
                field,
                atom,
                value: invalid,
            } => {
                value.setattr("field", *field)?;
                value.setattr("atom", atom.index())?;
                value.setattr("value", *invalid)?;
            }
            E::AtomOutOfRange { atom, atom_count } => {
                value.setattr("atom", atom.index())?;
                value.setattr("atom_count", *atom_count)?;
            }
            E::InvalidRingInfo {
                reason,
                row,
                value: invalid,
                limit,
            } => {
                value.setattr("reason", *reason)?;
                value.setattr("row", *row)?;
                value.setattr("value", *invalid)?;
                value.setattr("limit", *limit)?;
            }
            E::InvalidAtomDegree { atom, degree } => {
                value.setattr("atom", atom.index())?;
                value.setattr("degree", *degree)?;
            }
            E::InvalidBondDegree {
                bond,
                endpoint,
                degree,
            } => {
                value.setattr("bond", bond.index())?;
                value.setattr("endpoint", *endpoint)?;
                value.setattr("degree", *degree)?;
            }
            E::InvalidStereoReferences { bond, reason } => {
                value.setattr("bond", bond.index())?;
                value.setattr("reason", *reason)?;
            }
            E::UnsupportedBondOrder { bond, order } => {
                value.setattr("bond", bond.index())?;
                value.setattr(
                    "order",
                    py.import("cosmolkit")?
                        .getattr("BondOrder")?
                        .call1((order.rdkit_code(),))?,
                )?;
            }
            E::InvalidChiralPermutation { atom, permutation } => {
                value.setattr("atom", atom.index())?;
                value.setattr("permutation", *permutation)?;
            }
            E::AtropisomerDependencyUnavailable { bond } => value.setattr("bond", bond.index())?,
            E::RefinementDidNotConverge { iterations } => {
                value.setattr("iterations", *iterations)?
            }
            E::Ring(_)
            | E::InvalidTopology(_)
            | E::Valence(_)
            | E::StereoOrder(_)
            | E::DoubleStereo(_)
            | E::BondValue(_) => {}
        }
        Ok(())
    };
    match attributes() {
        Ok(()) => error,
        Err(error) => error,
    }
}

fn stereo_type(value: ck::PotentialStereoType) -> &'static str {
    match value {
        ck::PotentialStereoType::AtomTetrahedral => "atom_tetrahedral",
        ck::PotentialStereoType::AtomSquarePlanar => "atom_squareplanar",
        ck::PotentialStereoType::AtomTrigonalBipyramidal => "atom_trigonalbipyramidal",
        ck::PotentialStereoType::AtomOctahedral => "atom_octahedral",
        ck::PotentialStereoType::BondDouble => "bond_double",
        ck::PotentialStereoType::BondCumuleneEven => "bond_cumulene_even",
        ck::PotentialStereoType::BondAtropisomer => "bond_atropisomer",
    }
}
fn specified(value: ck::PotentialStereoSpecified) -> &'static str {
    match value {
        ck::PotentialStereoSpecified::Unspecified => "unspecified",
        ck::PotentialStereoSpecified::Specified => "specified",
        ck::PotentialStereoSpecified::Unknown => "unknown",
    }
}
fn descriptor(value: ck::PotentialStereoDescriptor) -> &'static str {
    match value {
        ck::PotentialStereoDescriptor::None => "none",
        ck::PotentialStereoDescriptor::TetrahedralClockwise => "tetrahedral_clockwise",
        ck::PotentialStereoDescriptor::TetrahedralCounterclockwise => {
            "tetrahedral_counterclockwise"
        }
        ck::PotentialStereoDescriptor::BondCis => "bond_cis",
        ck::PotentialStereoDescriptor::BondTrans => "bond_trans",
        ck::PotentialStereoDescriptor::BondAtropCw => "bond_atrop_cw",
        ck::PotentialStereoDescriptor::BondAtropCcw => "bond_atrop_ccw",
    }
}
fn enum_value(py: Python<'_>, name: &str, value: &str) -> PyResult<Py<PyAny>> {
    Ok(py
        .import("cosmolkit")?
        .getattr(name)?
        .call1((value,))?
        .unbind())
}
/// Writable configuration for potential stereocenter perception.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct PotentialStereoParams {
    pub(crate) inner: ck::PotentialStereoParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl PotentialStereoParams {
    /// Configure potential stereocenter perception; omitted fields use the defaults shown in the signature.
    #[new]
    #[pyo3(signature=(*, clean=false, flag_possible=true, allow_nontetrahedral=true))]
    fn new(clean: bool, flag_possible: bool, allow_nontetrahedral: bool) -> Self {
        Self {
            inner: ck::PotentialStereoParams {
                clean,
                flag_possible,
                allow_nontetrahedral,
            },
        }
    }
    /// Whether invalid stereo annotations are cleaned during perception.
    #[getter]
    fn clean(&self) -> bool {
        self.inner.clean
    }
    /// Whether potential stereocenters are flagged even if unspecified.
    #[getter]
    fn flag_possible(&self) -> bool {
        self.inner.flag_possible
    }
    /// Whether non-tetrahedral stereochemistry is considered.
    #[getter]
    fn allow_nontetrahedral(&self) -> bool {
        self.inner.allow_nontetrahedral
    }
    fn __repr__(&self) -> String {
        format!(
            "PotentialStereoParams(clean={}, flag_possible={}, allow_nontetrahedral={})",
            self.inner.clean, self.inner.flag_possible, self.inner.allow_nontetrahedral
        )
    }
}
/// Index and category of a potential atom/bond stereo feature.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct PotentialStereoCenter {
    inner: ck::PotentialStereoCenter,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl PotentialStereoCenter {
    /// Classification/discriminant of this value as defined by its owning type.
    #[getter]
    fn kind(&self) -> &'static str {
        match self.inner {
            ck::PotentialStereoCenter::Atom(_) => "atom",
            ck::PotentialStereoCenter::Bond(_) => "bond",
        }
    }
    /// Original zero-based input record index; failed records retain their index.
    #[getter]
    fn index(&self) -> usize {
        match self.inner {
            ck::PotentialStereoCenter::Atom(id) => id.index(),
            ck::PotentialStereoCenter::Bond(id) => id.index(),
        }
    }
    fn __repr__(&self) -> String {
        format!(
            "PotentialStereoCenter(kind='{}', index={})",
            self.kind(),
            self.index()
        )
    }
}
/// Potential stereo feature with its center, specification state, descriptor and controlling atoms.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct PotentialStereoInfo {
    inner: ck::PotentialStereoInfo,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl PotentialStereoInfo {
    /// Type of potential stereogenic atom or bond.
    #[getter]
    fn stereo_type(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        enum_value(
            py,
            "PotentialStereoType",
            stereo_type(self.inner.stereo_type),
        )
    }
    /// Whether the stereochemistry is specified, unspecified or unknown.
    #[getter]
    fn specified(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        enum_value(
            py,
            "PotentialStereoSpecified",
            specified(self.inner.specified),
        )
    }
    /// Atom or bond index on which the stereo feature is centered.
    #[getter]
    fn centered_on(&self) -> PotentialStereoCenter {
        PotentialStereoCenter {
            inner: self.inner.centered_on,
        }
    }
    /// Local stereo descriptor retained by the perception result.
    #[getter]
    fn descriptor(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        enum_value(
            py,
            "PotentialStereoDescriptor",
            descriptor(self.inner.descriptor),
        )
    }
    /// Permutation code for non-tetrahedral stereochemistry, when applicable.
    #[getter]
    fn permutation(&self) -> u32 {
        self.inner.permutation
    }
    /// Atom indices defining the local stereo feature.
    #[getter]
    fn controlling_atoms(&self) -> Vec<Option<usize>> {
        self.inner
            .controlling_atoms
            .iter()
            .map(|v| v.map(ck::AtomId::index))
            .collect()
    }
    fn __repr__(&self) -> String {
        format!(
            "PotentialStereoInfo(stereo_type='{}', specified='{}', centered_on={}, descriptor='{}', permutation={})",
            stereo_type(self.inner.stereo_type),
            specified(self.inner.specified),
            self.centered_on().__repr__(),
            descriptor(self.inner.descriptor),
            self.inner.permutation
        )
    }
}
/// Relative stereochemical orientation between two ring stereocenters.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct RingStereoRelation {
    inner: ck::RingStereoRelation,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl RingStereoRelation {
    /// Zero-based atom index of the first stereocenter in this ring relationship.
    #[getter]
    fn atom(&self) -> usize {
        self.inner.atom.index()
    }
    /// Index of the other stereocenter in this ring relationship.
    #[getter]
    fn other(&self) -> usize {
        self.inner.other.index()
    }
    /// Whether the linked ring stereocenters have the same relative orientation.
    #[getter]
    fn same_orientation(&self) -> bool {
        self.inner.same_orientation
    }
}
/// Potential-stereochemistry analysis with stereocenter rows, atom ranks and ring relationships. The cleaned_molecule is present only when clean=True requested a cleaned result.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct PotentialStereoResult {
    pub(crate) inner: ck::PotentialStereoResult,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl PotentialStereoResult {
    /// Return the perceived PotentialStereoInfo rows for potential atom and bond stereocenters.
    #[getter]
    fn stereo(&self) -> Vec<PotentialStereoInfo> {
        self.inner
            .stereo
            .iter()
            .cloned()
            .map(|inner| PotentialStereoInfo { inner })
            .collect()
    }
    /// Atom-indexed ranks used by stereo perception.
    #[getter]
    fn atom_ranks(&self) -> Vec<u32> {
        self.inner.atom_ranks.clone()
    }
    /// Relative stereochemical relationships between ring centers.
    #[getter]
    fn ring_relations(&self) -> Vec<RingStereoRelation> {
        self.inner
            .ring_relations
            .iter()
            .copied()
            .map(|inner| RingStereoRelation { inner })
            .collect()
    }
    /// Molecule value with the requested stereo cleanup applied; source input is unchanged.
    #[getter]
    fn cleaned_molecule(&self) -> Option<crate::drawing_binding::Molecule> {
        self.inner
            .cleaned_molecule
            .clone()
            .map(|inner| crate::drawing_binding::Molecule { inner })
    }
    fn __len__(&self) -> usize {
        self.inner.stereo.len()
    }
    fn __repr__(&self) -> String {
        format!("PotentialStereoResult(records={})", self.inner.stereo.len())
    }
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "PotentialStereoError",
        module.py().get_type::<PotentialStereoError>(),
    )?;
    module.add_class::<PotentialStereoParams>()?;
    module.add_class::<PotentialStereoCenter>()?;
    module.add_class::<PotentialStereoInfo>()?;
    module.add_class::<RingStereoRelation>()?;
    module.add_class::<PotentialStereoResult>()?;
    let enum_ = module.py().import("enum")?.getattr("Enum")?;
    let string_type = module.py().get_type::<pyo3::types::PyString>();
    macro_rules! vocabulary {($name:literal,$function:ident,[$($variant:path),+])=>{{
        let members=pyo3::types::PyDict::new(module.py());
        $(let value=$function($variant); members.set_item(value.to_ascii_uppercase(),value)?;)+
        let kwargs=pyo3::types::PyDict::new(module.py());kwargs.set_item("module","cosmolkit")?;
        kwargs.set_item("type", &string_type)?;
        let vocabulary = enum_.call(($name,members),Some(&kwargs))?;
        vocabulary.setattr("__str__", string_type.getattr("__str__")?)?;
        vocabulary.setattr("__format__", string_type.getattr("__format__")?)?;
        module.add($name,vocabulary)?;
    }}}
    vocabulary!(
        "PotentialStereoType",
        stereo_type,
        [
            ck::PotentialStereoType::AtomTetrahedral,
            ck::PotentialStereoType::AtomSquarePlanar,
            ck::PotentialStereoType::AtomTrigonalBipyramidal,
            ck::PotentialStereoType::AtomOctahedral,
            ck::PotentialStereoType::BondDouble,
            ck::PotentialStereoType::BondCumuleneEven,
            ck::PotentialStereoType::BondAtropisomer
        ]
    );
    vocabulary!(
        "PotentialStereoSpecified",
        specified,
        [
            ck::PotentialStereoSpecified::Unspecified,
            ck::PotentialStereoSpecified::Specified,
            ck::PotentialStereoSpecified::Unknown
        ]
    );
    vocabulary!(
        "PotentialStereoDescriptor",
        descriptor,
        [
            ck::PotentialStereoDescriptor::None,
            ck::PotentialStereoDescriptor::TetrahedralClockwise,
            ck::PotentialStereoDescriptor::TetrahedralCounterclockwise,
            ck::PotentialStereoDescriptor::BondCis,
            ck::PotentialStereoDescriptor::BondTrans,
            ck::PotentialStereoDescriptor::BondAtropCw,
            ck::PotentialStereoDescriptor::BondAtropCcw
        ]
    );
    for (name, description) in [
        (
            "PotentialStereoType",
            "Type of a potential stereocenter, distinguishing atom geometries from double, cumulene and atropisomeric bonds.",
        ),
        (
            "PotentialStereoSpecified",
            "Whether potential stereochemistry is unspecified, specified or explicitly unknown.",
        ),
        (
            "PotentialStereoDescriptor",
            "Perceived local stereo orientation, including tetrahedral handedness, cis/trans and atropisomeric directions.",
        ),
    ] {
        module.getattr(name)?.setattr("__doc__", description)?;
    }
    Ok(())
}
