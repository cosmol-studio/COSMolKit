use cosmolkit as ck;
use pyo3::prelude::*;

macro_rules! flags {
    ($name:ident { $($flag:ident),+ $(,)? }) => {
        crate::configuration::type_projection! {
            #[pyclass(module="cosmolkit", frozen, from_py_object, eq)]
            #[derive(Clone, Copy, PartialEq)]
            pub struct $name { pub(crate) inner: ck::$name }
        }
        #[pymethods]
        impl $name {
            $(#[classattr] #[allow(non_snake_case)] fn $flag() -> Self { Self { inner: ck::$name::$flag } })+
            #[staticmethod]
            fn from_bits(bits: u32) -> Self { Self { inner: ck::$name::from_bits(bits) } }
            #[getter]
            fn bits(&self) -> u32 { self.inner.bits() }
            fn contains(&self, flags: &Self) -> bool { self.inner.contains(flags.inner) }
            #[getter]
            fn is_empty(&self) -> bool { self.inner.is_empty() }
            fn __index__(&self) -> u32 { self.inner.bits() }
            fn __or__(&self, other: &Self) -> Self { Self { inner: self.inner | other.inner } }
            fn __and__(&self, other: &Self) -> Self { Self { inner: self.inner & other.inner } }
        }
    };
}
flags!(PdbReadFlags {
    NONE,
    KEEP_ALTERNATE_AND_PSEUDO_ATOMS,
    STANDARD_RESIDUE_BOND_ORDERS
});
flags!(PdbWriteFlags {
    NONE,
    OMIT_SINGLE_CONECT,
    BIDIRECTIONAL_CONECT,
    NO_BOND_ORDER,
    MASTER,
    TER
});
