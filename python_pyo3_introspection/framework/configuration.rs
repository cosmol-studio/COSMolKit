//! Field projection is selected by Rust types, not a handwritten parameter field list.
use crate::parameters::{Container, Dict, List, Scalar};
use std::collections::BTreeMap;

pub(crate) trait FieldProjection: Sized + Clone {
    type Adapter: Container<Native = Self>;
}
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
    T: Clone,
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

macro_rules! value_fields {
    ($($native:ty => $projection:ty),+ $(,)?) => {$(
        impl FieldProjection for $native {
            type Adapter = crate::parameters::ValueField<$projection>;
        }
    )+};
}
value_fields! {
    cosmolkit::Fingerprint => crate::native::Fingerprint,
    cosmolkit::AlignmentAtomMap => crate::collection_params::AlignmentAtomMap,
    cosmolkit::ReactionCoordinateSelection => crate::construction_values::ReactionCoordinateSelection,
    cosmolkit::CxSmilesFields => crate::notation_params::CxSmilesFields,
    cosmolkit::CxCoordinateSelection => crate::notation_params::CxCoordinateSelection,
    cosmolkit::TautomerScoreTerm => crate::tautomer_values::TautomerScoreTerm,
    cosmolkit::MorganInvariants => crate::morgan_configuration::MorganInvariants,
    cosmolkit::MolCoordinateSelection => crate::io_batch_params::MolCoordinateSelection,
    cosmolkit::SanitizeOperations => crate::sanitize_values::SanitizeOperations,
}
impl FieldProjection for cosmolkit::AtomId {
    type Adapter = crate::parameters::NativeScalarField<Self>;
}
impl FieldProjection for cosmolkit::BondId {
    type Adapter = crate::parameters::NativeScalarField<Self>;
}
impl FieldProjection for cosmolkit::LayeredFingerprintLayers {
    type Adapter = crate::parameters::NativeScalarField<Self>;
}
impl crate::parameters::NativeScalar for cosmolkit::AvalonFingerprintFlags {
    type Repr = u32;
    fn from_scalar(value: u32) -> Self {
        Self::from_bits_retain(value)
    }
    fn scalar(self) -> u32 {
        self.bits()
    }
}
impl FieldProjection for cosmolkit::AvalonFingerprintFlags {
    type Adapter = crate::parameters::NativeScalarField<Self>;
}

/// Counted native byte text is exposed as Python text only by strict decoding.
/// This shares the existing explicit UTF-8 boundary; no lossy replacement.
pub(crate) struct PropertyTextField;
impl Container for PropertyTextField {
    type Native = cosmolkit::PropertyText;
    const TYPE: pyo3::inspect::PyStaticExpr = pyo3::type_hint_identifier!("builtins", "str");
    const INPUT: pyo3::inspect::PyStaticExpr = Self::TYPE;
    fn extract(value: pyo3::Borrowed<'_, '_, pyo3::PyAny>) -> pyo3::PyResult<Self::Native> {
        Ok(value.extract::<String>()?.into())
    }
    fn output<'py>(
        py: pyo3::Python<'py>,
        value: &Self::Native,
    ) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::PyAny>> {
        use pyo3::IntoPyObjectExt;
        crate::property_values::decode_source_text(py, value)?.into_bound_py_any(py)
    }
}
impl FieldProjection for cosmolkit::PropertyText {
    type Adapter = PropertyTextField;
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
assert p.generator.radius == 3
p.generator = child
assert p.generator.radius == 0
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
            assert_eq!(object.borrow(py).inner.generator.radius, 0);
            assert_eq!(object.borrow(py).inner.limits, [2, 4]);
            Ok(())
        })
        .unwrap();
    }
}
