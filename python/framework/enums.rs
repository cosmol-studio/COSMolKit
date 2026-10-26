//! Protocols for native, payload-free enums. No domain vocabulary lives here.
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyMappingProxy};
use std::marker::PhantomData;

crate::callbacks::typed_callback!(
    EnumMemberLookup,
    [
        pyo3::type_hint_identifier!("builtins", "type"),
        pyo3::type_hint_identifier!("builtins", "str")
    ],
    pyo3::type_hint_identifier!("builtins", "object")
);

pub(crate) struct Member<T> {
    pub name: &'static str,
    pub aliases: &'static [&'static str],
    pub value: T,
}

pub(crate) trait EnumMembers: pyo3::PyClass<BaseType = PyAny> + Copy + 'static {
    const MEMBERS: &'static [Member<Self>];
    fn number(self) -> isize;
}

pub(crate) struct Names<T>(pub PhantomData<T>);
pub(crate) struct Values<T>(pub PhantomData<T>);

impl<'py, T: EnumMembers> IntoPyObject<'py> for Names<T> {
    type Target = PyMappingProxy;
    type Output = Bound<'py, PyMappingProxy>;
    type Error = PyErr;
    const OUTPUT_TYPE: pyo3::inspect::PyStaticExpr = pyo3::type_hint_subscript!(
        pyo3::type_hint_identifier!("collections.abc", "Mapping"),
        pyo3::type_hint_identifier!("builtins", "str"),
        <Py<T> as IntoPyObject<'py>>::OUTPUT_TYPE
    );
    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        let mapping = PyDict::new(py);
        for member in T::MEMBERS {
            let value = Py::new(py, member.value)?;
            for name in std::iter::once(&member.name).chain(member.aliases) {
                if mapping.contains(name)? {
                    return Err(pyo3::exceptions::PyValueError::new_err(format!(
                        "duplicate enum name {name}"
                    )));
                }
                mapping.set_item(name, &value)?;
            }
        }
        Ok(PyMappingProxy::new(py, mapping.as_mapping()))
    }
}

impl<'py, T: EnumMembers> IntoPyObject<'py> for Values<T> {
    type Target = PyMappingProxy;
    type Output = Bound<'py, PyMappingProxy>;
    type Error = PyErr;
    const OUTPUT_TYPE: pyo3::inspect::PyStaticExpr = pyo3::type_hint_subscript!(
        pyo3::type_hint_identifier!("collections.abc", "Mapping"),
        pyo3::type_hint_identifier!("builtins", "int"),
        <Py<T> as IntoPyObject<'py>>::OUTPUT_TYPE
    );
    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        let mapping = PyDict::new(py);
        for member in T::MEMBERS {
            mapping.set_item(member.value.number(), Py::new(py, member.value)?)?;
        }
        Ok(PyMappingProxy::new(py, mapping.as_mapping()))
    }
}

macro_rules! enum_protocol {
    (@select [] $name:ident, $members:tt, methods = $methods:tt) => {
        $crate::enums::enum_protocol!($name, $members, methods = $methods);
    };
    (@select [custom] $name:ident, $members:tt, methods = $methods:tt) => {
        $crate::enums::enum_protocol!(custom $name, $members, methods = $methods);
    };
    ($name:ident, {$($variant:ident => ($canonical:expr, [$($alias:ident),*])),+ $(,)?}
        $(, methods = {$($method:item)*})?) => {
        $crate::enums::enum_protocol!(@emit $name,
            {$($variant => ($canonical, [$($alias),*])),+}, {
                #[new]
                fn new(value: isize) -> pyo3::PyResult<Self> {
                    match value {
                        $(value if value == Self::$variant as isize => Ok(Self::$variant),)+
                        _ => Err(pyo3::exceptions::PyValueError::new_err(format!(
                            "invalid {} value: {value}", stringify!($name)
                        ))),
                    }
                }
                fn __reduce__<'py>(
                    &self,
                    py: pyo3::Python<'py>,
                ) -> (pyo3::Bound<'py, pyo3::types::PyType>, (isize,)) {
                    (py.get_type::<Self>(), (*self as isize,))
                }
                $($($method)*)?
            });
    };
    (custom $name:ident, $members:tt, methods = {$($method:item)*}) => {
        $crate::enums::enum_protocol!(@emit $name, $members, {
            // Custom constructors need not accept numeric discriminants.
            // Canonical members already exist on the public class; no new
            // reconstruction entrypoint or parallel member table is needed.
            fn __reduce__<'py>(
                &self,
                py: pyo3::Python<'py>,
            ) -> pyo3::PyResult<(
                $crate::enums::EnumMemberLookup,
                (pyo3::Bound<'py, pyo3::types::PyType>, &'static str),
            )> {
                use pyo3::types::PyAnyMethods;
                let getattr = py.import("builtins")?.getattr("getattr")?
                    .extract::<$crate::enums::EnumMemberLookup>()?;
                Ok((getattr, (py.get_type::<Self>(), self.name())))
            }
            $($method)*
        });
    };
    (@emit $name:ident, {$($variant:ident => ($canonical:expr, [$($alias:ident),*])),+ $(,)?}, {$($method:item)*}) => {
        impl $crate::enums::EnumMembers for $name {
            const MEMBERS: &'static [$crate::enums::Member<Self>] = &[
                $($crate::enums::Member {
                    name: $canonical,
                    aliases: &[$(stringify!($alias)),*],
                    value: Self::$variant,
                }),+
            ];
            fn number(self) -> isize { self as isize }
        }
        #[pyo3::pymethods]
        impl $name {
            $($(#[classattr]
            #[allow(non_upper_case_globals)]
            const $alias: Self = Self::$variant;)*)+
            /// Canonical member name; aliases of the same value share this name.
            #[getter]
            fn name(&self) -> &'static str {
                match self { $(Self::$variant => $canonical),+ }
            }
            #[classattr]
            fn names() -> $crate::enums::Names<Self> {
                $crate::enums::Names(std::marker::PhantomData)
            }
            #[classattr]
            fn values() -> $crate::enums::Values<Self> {
                $crate::enums::Values(std::marker::PhantomData)
            }
            fn __hash__(&self, py: pyo3::Python<'_>) -> pyo3::PyResult<isize> {
                use pyo3::{IntoPyObject, types::PyAnyMethods};
                (*self as isize).into_pyobject(py)?.hash()
            }
            $($method)*
        }
    };
}
pub(crate) use enum_protocol;

macro_rules! public_name {
    ($variant:ident;) => {
        stringify!($variant)
    };
    ($variant:ident; #[pyo3(name = $name:literal)]) => {
        $name
    };
}
pub(crate) use public_name;

// Error discriminants without a directly corresponding native enum still have
// one declaration. Existing discriminant expressions are preserved verbatim.
macro_rules! binding_enum {
    ($name:ident, {$($variant:ident $(= $value:expr)?),+ $(,)?}) => {
        #[doc = concat!("Discriminants of ", stringify!($name), ".\n\n",
            "Members: ", $("``", stringify!($variant), "``, ",)+
            "\nCompare members by value. names and values are readonly lookup mappings; ",
            "pickle preserves the enum type and value.")]
        #[pyo3::pyclass(module = "cosmolkit", eq, eq_int, from_py_object)]
        #[derive(Clone, Copy, PartialEq)]
        pub enum $name { $($variant $(= $value)?),+ }
        $crate::enums::enum_protocol!($name, {
            $($variant => (stringify!($variant), [])),+
        });
    };
}
pub(crate) use binding_enum;
