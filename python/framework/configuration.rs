//! Field projection is selected by Rust types, not a handwritten parameter field list.
use crate::parameters::{Container, Dict, List, Scalar};
use std::collections::BTreeMap;

pub(crate) trait FieldProjection: Sized + Clone {
    type Adapter: Container<Native = Self>;
}

/// A binding declares its conversion and native field association together.
/// No business types or runtime lookup table belong in this macro.
macro_rules! type_projection {
    (value(eq) $($declaration:tt)*) => {
        $crate::configuration::type_projection! {
            #[pyo3::pyclass(module = "cosmolkit", frozen, from_py_object, eq)]
            #[derive(Clone, PartialEq)]
            $($declaration)*
        }
    };
    (value(eq, hash) $($declaration:tt)*) => {
        $crate::configuration::type_projection! {
            #[pyo3::pyclass(module = "cosmolkit", frozen, from_py_object, eq, hash)]
            #[derive(Clone, PartialEq, Eq, Hash)]
            $($declaration)*
        }
    };
    (@association $(#[$($attribute:tt)*])* $visibility:vis struct $wrapper:ident {
        $field_visibility:vis inner: $native:ty $(,)?
    }) => {
        $crate::configuration::type_projection! {
            impl Configuration for $wrapper {
                type Native = $native;
                fn from_native(_py: pyo3::Python<'_>, inner: Self::Native) -> pyo3::PyResult<Self> {
                    Ok(Self { inner })
                }
                fn to_native(&self, _py: pyo3::Python<'_>) -> pyo3::PyResult<Self::Native> {
                    Ok(self.inner.clone())
                }
            }
        }
    };
    (impl Configuration for $wrapper:ty {
        type Native = $native:ty;
        $($body:tt)*
    }) => {
        impl crate::parameters::Configuration for $wrapper {
            type Native = $native;
            $($body)*
        }
        impl crate::configuration::FieldProjection for $native {
            type Adapter = crate::parameters::ValueField<$wrapper>;
        }
    };
    (impl NativeScalar for $native:ty {
        $($body:tt)*
    }) => {
        impl crate::parameters::NativeScalar for $native { $($body)* }
        impl crate::configuration::FieldProjection for $native {
            type Adapter = crate::parameters::NativeScalarField<Self>;
        }
    };
    (impl Container for $adapter:ty {
        type Native = $native:ty;
        $($body:tt)*
    }) => {
        impl crate::parameters::Container for $adapter {
            type Native = $native;
            $($body)*
        }
        impl crate::configuration::FieldProjection for $native {
            type Adapter = $adapter;
        }
    };
    (@association $($invalid:tt)*) => {
        compile_error!("value projection requires a single native inner field; stateful owners need an explicit adapter");
    };
    ($($declaration:tt)*) => {
        // Preserve the original item tokens, including attribute spans. Rebuilding
        // #[pyclass(...)] changes the hygiene of PyO3-generated comparison code.
        $($declaration)*
        $crate::configuration::type_projection!(@association $($declaration)*);
    };
}
pub(crate) use type_projection;
macro_rules! scalars {
    ($($ty:ty),+ $(,)?) => {$(
        impl FieldProjection for $ty { type Adapter = Scalar<Self>; }
    )+};
}
scalars!(
    bool,
    u8,
    u16,
    u32,
    u64,
    usize,
    i8,
    i16,
    i32,
    i64,
    isize,
    f32,
    f64,
    String,
    std::path::PathBuf,
    (u32, u32),
    (i32, i32)
);
impl<T, const N: usize> FieldProjection for [T; N]
where
    T: Clone + 'static,
    for<'a, 'py> T: pyo3::FromPyObject<'a, 'py>,
    for<'py> T: pyo3::IntoPyObject<'py>,
{
    type Adapter = Self;
}
impl<T: FieldProjection> FieldProjection for Vec<T> {
    type Adapter = List<T::Adapter>;
}
impl<T: FieldProjection> FieldProjection for Option<T> {
    type Adapter = Option<T::Adapter>;
}
impl<K: FieldProjection + Ord, V: FieldProjection> FieldProjection for BTreeMap<K, V> {
    type Adapter = Dict<K::Adapter, V::Adapter>;
}

#[cfg(test)]
pub mod fixture {
    use std::collections::BTreeMap;

    #[derive(Clone, cosmolkit::Configuration)]
    pub struct ProjectionChild {
        pub radius: u32,
    }
    impl Default for ProjectionChild {
        fn default() -> Self {
            Self { radius: 3 }
        }
    }

    #[derive(Clone, Default, cosmolkit::Configuration)]
    pub struct ProjectionExample {
        pub generator: ProjectionChild,
        pub optional: Option<ProjectionChild>,
        pub from_atoms: Vec<u32>,
        pub replacements: BTreeMap<String, String>,
        pub limits: [u32; 2],
    }

    pub use __ck_configuration_ProjectionChild;
    pub use __ck_configuration_ProjectionExample;
}

#[cfg(test)]
mod tests {
    use super::fixture as native;
    use pyo3::prelude::*;
    cosmolkit::project_configuration!(native::ProjectionChild);
    cosmolkit::project_configuration!(native::ProjectionExample);

    #[derive(Clone, PartialEq)]
    struct PartialValue(f64);

    super::type_projection! {
        value(eq)
        struct PartialValueWrapper { inner: PartialValue }
    }

    #[test]
    fn value_projection_preserves_native_equality_and_hash_contracts() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::new(py, "value_projection_sample")?;
            module.add_class::<crate::elements::Element>()?;
            module.add(
                "a",
                Py::new(
                    py,
                    PartialValueWrapper {
                        inner: PartialValue(2.0),
                    },
                )?,
            )?;
            module.add(
                "b",
                Py::new(
                    py,
                    PartialValueWrapper {
                        inner: PartialValue(2.0),
                    },
                )?,
            )?;
            module.add(
                "nan",
                Py::new(
                    py,
                    PartialValueWrapper {
                        inner: PartialValue(f64::NAN),
                    },
                )?,
            )?;
            py.run(
                pyo3::ffi::c_str!(
                    r#"
chlorine = Element.from_symbol("Cl")
assert chlorine is not Element.CL
assert chlorine == Element.CL and not (chlorine != Element.CL)
assert chlorine != Element.C
assert chlorine != 17 and chlorine != "Cl"
assert chlorine.__eq__(17) is NotImplemented
assert hash(chlorine) == hash(Element.CL)
assert len({chlorine, Element.CL, Element.C}) == 2
assert {chlorine: "chlorine"}[Element.CL] == "chlorine"
assert a is not b and a == b and not (a != b)
assert nan != nan
assert type(a).__hash__ is None
try:
    hash(a)
except TypeError:
    pass
else:
    raise AssertionError("equality-only values must be unhashable")
"#
                ),
                Some(&module.dict()),
                None,
            )?;
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn native_declared_fields_use_atomic_assignment_and_immutable_reads() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::new(py, "configuration_sample")?;
            module.add_class::<ProjectionChild>()?;
            module.add_class::<ProjectionExample>()?;
            py.run(
                pyo3::ffi::c_str!(
                    r#"
source = ProjectionChild(radius=3)
indices = [0]
rules = {"{A}": "O"}
p = ProjectionExample(generator=source, from_atoms=indices, replacements=rules)
source.radius = 9
indices.append(9)
rules["{A}"] = "C"
assert p.generator.radius == 3 and p.from_atoms == (0,)
assert p.replacements == {"{A}": "O"}
child = p.generator
child.radius = 0
assert p.generator.radius == 0
p.generator = p.generator
p.generator = ProjectionChild(radius=4)
assert child.radius == 4
child.radius = 0
assert p.generator.radius == 0
p.optional = source
optional = p.optional
optional.radius = 6
assert p.optional.radius == 6 and source.radius == 9
p.optional = None
try:
    optional.radius = 1
except ReferenceError:
    pass
else:
    raise AssertionError("removed optional field remained writable")
p.from_atoms = [0, 1]
p.replacements = {"{B}": "N"}
p.limits = [2, 4]
assert p.limits == (2, 4)
for action, kind in [
    (lambda: p.from_atoms.append(2), AttributeError),
    (lambda: p.replacements.__setitem__("{B}", "C"), AttributeError),
    (lambda: setattr(p, "from_atoms", ["bad"]), TypeError),
    (lambda: setattr(p, "from_atoms", [-1]), OverflowError),
    (lambda: setattr(p, "replacements", {1: "C"}), TypeError),
    (lambda: setattr(p, "limits", [9]), ValueError),
]:
    try:
        action()
    except kind:
        pass
    else:
        raise AssertionError("invalid modification was accepted")
assert p.from_atoms == (0, 1) and p.replacements == {"{B}": "N"}
assert p.limits == (2, 4)
"#
                ),
                Some(&module.dict()),
                None,
            )?;
            let object = module.getattr("p")?.extract::<Py<ProjectionExample>>()?;
            let captured = object.borrow(py).to_rust(py)?;
            assert_eq!(captured.generator.radius, 0);
            assert_eq!(captured.from_atoms, [0, 1]);
            assert_eq!(captured.replacements["{B}"], "N");
            module.getattr("child")?.setattr("radius", 2)?;
            assert_eq!(captured.generator.radius, 0);
            assert_eq!(object.borrow(py).to_rust(py)?.generator.radius, 2);
            assert_eq!(object.borrow(py).to_rust(py)?.limits, [2, 4]);
            module.dict().del_item("p")?;
            drop(object);
            module.getattr("child")?.setattr("radius", 7)?;
            assert_eq!(
                module
                    .getattr("child")?
                    .getattr("radius")?
                    .extract::<u32>()?,
                7
            );
            Ok(())
        })
        .unwrap();
    }
}
