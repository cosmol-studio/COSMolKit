use crate::chemical_values::{EnumValue, enum_value};
use crate::native::Molecule;
use cosmolkit as ck;
use pyo3::prelude::*;

enum_value!(PotentialStereoType, {
    AtomTetrahedral [ATOM_TETRAHEDRAL], AtomSquarePlanar [ATOM_SQUAREPLANAR], AtomTrigonalBipyramidal [ATOM_TRIGONALBIPYRAMIDAL], AtomOctahedral [ATOM_OCTAHEDRAL],
    BondDouble [BOND_DOUBLE], BondCumuleneEven [BOND_CUMULENE_EVEN], BondAtropisomer [BOND_ATROPISOMER]
});
enum_value!(PotentialStereoSpecified, { Unspecified [UNSPECIFIED], Specified [SPECIFIED], Unknown [UNKNOWN] });
enum_value!(PotentialStereoDescriptor, {
    #[pyo3(name="None_")] None [NONE], TetrahedralClockwise [TETRAHEDRAL_CLOCKWISE], TetrahedralCounterclockwise [TETRAHEDRAL_COUNTERCLOCKWISE],
    BondCis [BOND_CIS], BondTrans [BOND_TRANS], BondAtropCw [BOND_ATROP_CW], BondAtropCcw [BOND_ATROP_CCW]
});

macro_rules! enum_input {
    ($input:ident, $value:ident) => {
        impl $input {
            fn into_rust(self) -> PyResult<ck::$value> {
                match self {
                    Self::Enum(value) => Ok(value.into()),
                    Self::Name(name) => $value::from_name(&name).map(Into::into),
                }
            }
        }
    };
}
#[derive(FromPyObject, IntoPyObject)]
pub enum StereoTypeInput {
    Enum(PotentialStereoType),
    Name(String),
}
#[derive(FromPyObject, IntoPyObject)]
pub enum StereoSpecifiedInput {
    Enum(PotentialStereoSpecified),
    Name(String),
}
#[derive(FromPyObject, IntoPyObject)]
pub enum StereoDescriptorInput {
    Enum(PotentialStereoDescriptor),
    Name(String),
}
enum_input!(StereoTypeInput, PotentialStereoType);
enum_input!(StereoSpecifiedInput, PotentialStereoSpecified);
enum_input!(StereoDescriptorInput, PotentialStereoDescriptor);

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone, Copy)]
    /// Atom-or-bond location of a potential stereochemical element.
    pub struct PotentialStereoCenter {
        inner: ck::PotentialStereoCenter,
    }
}
impl From<ck::PotentialStereoCenter> for PotentialStereoCenter {
    fn from(value: ck::PotentialStereoCenter) -> Self {
        Self { inner: value }
    }
}
impl From<PotentialStereoCenter> for ck::PotentialStereoCenter {
    fn from(value: PotentialStereoCenter) -> Self {
        value.inner
    }
}
#[pymethods]
impl PotentialStereoCenter {
    #[getter]
    /// Whether the stereochemical center is an atom or a bond, returned as 'atom' or 'bond'.
    fn kind(&self) -> &'static str {
        match self.inner {
            ck::PotentialStereoCenter::Atom(_) => "atom",
            ck::PotentialStereoCenter::Bond(_) => "bond",
        }
    }
    #[getter]
    /// Zero-based atom or bond index, interpreted according to kind.
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

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Geometric type, specification state and controlling atoms for a potential stereochemical element.
    pub struct PotentialStereoInfo {
        pub(crate) inner: ck::PotentialStereoInfo,
    }
}
#[pymethods]
impl PotentialStereoInfo {
    #[new]
    fn new(
        stereo_type: StereoTypeInput,
        specified: StereoSpecifiedInput,
        centered_on: PotentialStereoCenter,
        descriptor: StereoDescriptorInput,
        permutation: u32,
        controlling_atoms: Vec<Option<usize>>,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: ck::PotentialStereoInfo {
                stereo_type: stereo_type.into_rust()?,
                specified: specified.into_rust()?,
                centered_on: centered_on.into(),
                descriptor: descriptor.into_rust()?,
                permutation,
                controlling_atoms: controlling_atoms
                    .into_iter()
                    .map(|value| value.map(ck::AtomId::new))
                    .collect(),
            },
        })
    }
    #[getter]
    /// Geometric type of the potential stereochemical element.
    fn stereo_type(&self) -> PotentialStereoType {
        self.inner.stereo_type.into()
    }
    #[getter]
    /// Whether stereochemistry is specified, unspecified or explicitly unknown.
    fn specified(&self) -> PotentialStereoSpecified {
        self.inner.specified.into()
    }
    #[getter]
    /// Atom or bond identifying the stereochemical center.
    fn centered_on(&self) -> PotentialStereoCenter {
        self.inner.centered_on.into()
    }
    #[getter]
    /// Native stereochemical descriptor for this element.
    fn descriptor(&self) -> PotentialStereoDescriptor {
        self.inner.descriptor.into()
    }
    #[getter]
    /// Permutation index used for non-tetrahedral stereochemistry.
    fn permutation(&self) -> u32 {
        self.inner.permutation
    }
    #[getter]
    /// Ordered controlling atom indices, with None for source sentinel positions.
    fn controlling_atoms(&self) -> Vec<Option<usize>> {
        self.inner
            .controlling_atoms
            .iter()
            .map(|value| value.map(ck::AtomId::index))
            .collect()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Relative stereochemical orientation of a pair of ring atoms.
    pub struct RingStereoRelation {
        pub(crate) inner: ck::RingStereoRelation,
    }
}
#[pymethods]
impl RingStereoRelation {
    #[new]
    fn new(atom: usize, other: usize, same_orientation: bool) -> Self {
        Self {
            inner: ck::RingStereoRelation {
                atom: ck::AtomId::new(atom),
                other: ck::AtomId::new(other),
                same_orientation,
            },
        }
    }
    #[getter]
    /// First atom participating in the ring stereochemistry relationship.
    fn atom(&self) -> usize {
        self.inner.atom.index()
    }
    #[getter]
    /// Second atom participating in the ring stereochemistry relationship.
    fn other(&self) -> usize {
        self.inner.other.index()
    }
    #[getter]
    /// Whether the two related ring atoms have the same orientation.
    fn same_orientation(&self) -> bool {
        self.inner.same_orientation
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Potential stereochemical elements, atom ranks, ring relationships and optional cleaned molecule.
    pub struct PotentialStereoResult {
        pub(crate) inner: ck::PotentialStereoResult,
    }
}
#[pymethods]
impl PotentialStereoResult {
    #[new]
    fn new(
        stereo: Vec<PotentialStereoInfo>,
        atom_ranks: Vec<u32>,
        ring_relations: Vec<RingStereoRelation>,
        cleaned_molecule: Option<&Molecule>,
    ) -> Self {
        Self {
            inner: ck::PotentialStereoResult {
                stereo: stereo.into_iter().map(|value| value.inner).collect(),
                atom_ranks,
                ring_relations: ring_relations
                    .into_iter()
                    .map(|value| value.inner)
                    .collect(),
                cleaned_molecule: cleaned_molecule.map(|value| value.inner.clone()),
            },
        }
    }
    #[getter]
    /// Independent list of identified potential stereochemical elements.
    fn stereo(&self) -> Vec<PotentialStereoInfo> {
        self.inner
            .stereo
            .iter()
            .cloned()
            .map(|inner| PotentialStereoInfo { inner })
            .collect()
    }
    #[getter]
    /// Native atom ranks used during potential stereochemistry detection.
    fn atom_ranks(&self) -> Vec<u32> {
        self.inner.atom_ranks.clone()
    }
    #[getter]
    /// Independent list of ring stereochemistry relationships.
    fn ring_relations(&self) -> Vec<RingStereoRelation> {
        self.inner
            .ring_relations
            .iter()
            .cloned()
            .map(|inner| RingStereoRelation { inner })
            .collect()
    }
    #[getter]
    /// COW molecular value with cleaned stereochemistry, or None when no cleaned value was produced.
    fn cleaned_molecule(&self) -> Option<Molecule> {
        self.inner
            .cleaned_molecule
            .as_ref()
            .cloned()
            .map(|inner| Molecule { inner })
    }
}
