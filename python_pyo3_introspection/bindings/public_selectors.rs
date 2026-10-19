use crate::chemical_values::enum_value;
use cosmolkit as ck;
use pyo3::prelude::*;

enum_value!(RotatableBondsOptions, { Default, NonStrict, Strict, StrictLinkages });
enum_value!(BioPdbReadStage, { Stream, Record, Finalization });
enum_value!(BioMmcifReadStage, {
    CifDocument, CoordinateBlock, CrystalCell, Refinement, Tls, Experimental,
    Reflections, Software, Ncs, FractionalTransform, Origx, AnisotropicU, AtomSites,
    EntitySequence, Helices, Sheets, Connections, CisPeptides, ModifiedResidues,
    Assemblies, SiftsUnp, CcdRestoration, Materialization, StructureValidation
});
enum_value!(ParityPolicy, { NotApplicable, RequiredWhenSupported, RequiredNow });
enum_value!(MolecularForceFieldErrorKind, {
    MissingConformer, Parameterization, InvalidParameterization, InvalidAtomIndex,
    InvalidFixedAtom, CoordinateCount, CoordinateShape, NonFiniteCoordinate,
    InvalidTolerance, InvalidConstraintBounds, InvalidConstraintRange,
    Preparation, Construction, Initialization, Rings, Kernel
});

macro_rules! flag_value {
    ($name:ident, {$($flag:ident),+ $(,)?}) => {
        #[pyclass(module="cosmolkit", frozen, from_py_object, eq)]
        #[derive(Clone, Copy, PartialEq)]
        pub struct $name { pub(crate) inner: ck::$name }
        #[pymethods]
        impl $name {
            $(#[classattr] const $flag: Self = Self { inner: ck::$name::$flag };)+
            #[staticmethod] fn from_bits_retain(bits: u32) -> Self { Self { inner: ck::$name::from_bits_retain(bits) } }
            fn bits(&self) -> u32 { self.inner.bits() }
            fn contains(&self, other: &Self) -> bool { self.inner.contains(other.inner) }
            fn __index__(&self) -> u32 { self.inner.bits() }
            fn __or__(&self, other: &Self) -> Self { Self { inner: self.inner | other.inner } }
        }
    };
}
flag_value!(LayeredFingerprintLayers, {
    TOPOLOGY, BOND_ORDER, ATOM_TYPE, RING_PRESENCE, RING_SIZE, AROMATICITY,
    ACTIVE, SUBSTRUCTURE, ALL_SOURCE_BITS
});

// Kept as an ordinary Rust item: PyO3's derives require non-macro item hygiene.
#[pyclass(module = "cosmolkit", frozen, from_py_object, eq)]
#[derive(Clone, Copy, PartialEq)]
pub struct AvalonFingerprintFlags {
    pub(crate) inner: ck::AvalonFingerprintFlags,
}
#[pymethods]
impl AvalonFingerprintFlags {
    #[classattr]
    const RING_PATTERN: Self = Self {
        inner: ck::AvalonFingerprintFlags::RING_PATTERN,
    };
    #[classattr]
    const RING_PATH: Self = Self {
        inner: ck::AvalonFingerprintFlags::RING_PATH,
    };
    #[classattr]
    const ATOM_SYMBOL_PATH: Self = Self {
        inner: ck::AvalonFingerprintFlags::ATOM_SYMBOL_PATH,
    };
    #[classattr]
    const ATOM_CLASS_PATH: Self = Self {
        inner: ck::AvalonFingerprintFlags::ATOM_CLASS_PATH,
    };
    #[classattr]
    const ATOM_COUNT: Self = Self {
        inner: ck::AvalonFingerprintFlags::ATOM_COUNT,
    };
    #[classattr]
    const AUGMENTED_ATOM: Self = Self {
        inner: ck::AvalonFingerprintFlags::AUGMENTED_ATOM,
    };
    #[classattr]
    const HCOUNT_PATH: Self = Self {
        inner: ck::AvalonFingerprintFlags::HCOUNT_PATH,
    };
    #[classattr]
    const HCOUNT_CLASS_PATH: Self = Self {
        inner: ck::AvalonFingerprintFlags::HCOUNT_CLASS_PATH,
    };
    #[classattr]
    const HCOUNT_PAIR: Self = Self {
        inner: ck::AvalonFingerprintFlags::HCOUNT_PAIR,
    };
    #[classattr]
    const BOND_PATH: Self = Self {
        inner: ck::AvalonFingerprintFlags::BOND_PATH,
    };
    #[classattr]
    const AUGMENTED_BOND: Self = Self {
        inner: ck::AvalonFingerprintFlags::AUGMENTED_BOND,
    };
    #[classattr]
    const RING_SIZE_COUNTS: Self = Self {
        inner: ck::AvalonFingerprintFlags::RING_SIZE_COUNTS,
    };
    #[classattr]
    const DEGREE_PATH: Self = Self {
        inner: ck::AvalonFingerprintFlags::DEGREE_PATH,
    };
    #[classattr]
    const CLASS_SPIDERS: Self = Self {
        inner: ck::AvalonFingerprintFlags::CLASS_SPIDERS,
    };
    #[classattr]
    const FEATURE_PAIRS: Self = Self {
        inner: ck::AvalonFingerprintFlags::FEATURE_PAIRS,
    };
    #[classattr]
    const SCAFFOLD_IDS: Self = Self {
        inner: ck::AvalonFingerprintFlags::SCAFFOLD_IDS,
    };
    #[classattr]
    const SCAFFOLD_COLORS: Self = Self {
        inner: ck::AvalonFingerprintFlags::SCAFFOLD_COLORS,
    };
    #[classattr]
    const SCAFFOLD_LINKS: Self = Self {
        inner: ck::AvalonFingerprintFlags::SCAFFOLD_LINKS,
    };
    #[classattr]
    const SHORTCUT_LABELS: Self = Self {
        inner: ck::AvalonFingerprintFlags::SHORTCUT_LABELS,
    };
    #[classattr]
    const ALL_FEATURES: Self = Self {
        inner: ck::AvalonFingerprintFlags::ALL_FEATURES,
    };
    #[classattr]
    const NON_SSS_BITS: Self = Self {
        inner: ck::AvalonFingerprintFlags::NON_SSS_BITS,
    };
    #[classattr]
    const SIMILARITY: Self = Self {
        inner: ck::AvalonFingerprintFlags::SIMILARITY,
    };
    #[staticmethod]
    fn from_bits(bits: u32) -> Option<Self> {
        ck::AvalonFingerprintFlags::from_bits(bits).map(|inner| Self { inner })
    }
    #[staticmethod]
    fn from_bits_retain(bits: u32) -> Self {
        Self {
            inner: ck::AvalonFingerprintFlags::from_bits_retain(bits),
        }
    }
    fn bits(&self) -> u32 {
        self.inner.bits()
    }
    fn contains(&self, other: &Self) -> bool {
        self.inner.contains(other.inner)
    }
    fn __index__(&self) -> u32 {
        self.inner.bits()
    }
    fn __or__(&self, other: &Self) -> Self {
        Self {
            inner: self.inner | other.inner,
        }
    }
}
