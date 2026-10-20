//! A parameter declaration generates construction and writable fields.
//! It deliberately does not generate domain calls or a second default policy.
use pyo3::BoundObject;
use pyo3::prelude::*;
use pyo3::type_hint_union;
use pyo3::types::{PyDict, PyDictMethods, PyMapping, PyMappingProxy, PyTuple};
use std::collections::BTreeMap;

/// Typed value conversion at the language boundary. Configuration storage stays
/// native; container reads are immutable snapshots, never editable mirrors.
pub(crate) trait Container: Sized {
    const OPTIONAL: bool = false;
    type Native: Clone + 'static;
    const TYPE: pyo3::inspect::PyStaticExpr;
    const INPUT: pyo3::inspect::PyStaticExpr;
    fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<Self::Native>;
    fn output<'py>(py: Python<'py>, value: &Self::Native) -> PyResult<Bound<'py, PyAny>>;
    fn output_field<'py, V: NativeView<Self::Native> + 'static>(
        py: Python<'py>,
        field: V,
    ) -> PyResult<Bound<'py, PyAny>> {
        let mut result = None;
        field.read(py, &mut |value| result = Some(Self::output(py, value)))?;
        result.unwrap()
    }
}

/// A short, checked borrow on each access, never a retained reference into a
/// Python allocation. The parent handle keeps the owner alive.
pub(crate) trait NativeView<T>: Send + Sync {
    fn read(&self, py: Python<'_>, visit: &mut dyn FnMut(&T)) -> PyResult<()>;
    fn write(&self, py: Python<'_>, visit: &mut dyn FnMut(&mut T)) -> PyResult<()>;
    fn traverse(&self, visit: pyo3::gc::PyVisit<'_>) -> Result<(), pyo3::PyTraverseError>;
}

pub(crate) enum NativeField<T> {
    Owned(T),
    View(Box<dyn NativeView<T>>),
}
impl<T> NativeField<T> {
    pub(crate) fn traverse(
        &self,
        visit: pyo3::gc::PyVisit<'_>,
    ) -> Result<(), pyo3::PyTraverseError> {
        match self {
            Self::Owned(_) => Ok(()),
            Self::View(view) => view.traverse(visit),
        }
    }
    pub(crate) fn read<R>(&self, py: Python<'_>, f: impl FnOnce(&T) -> R) -> PyResult<R> {
        match self {
            Self::Owned(value) => Ok(f(value)),
            Self::View(view) => {
                let mut f = Some(f);
                let mut result = None;
                view.read(py, &mut |value| result = Some(f.take().unwrap()(value)))?;
                Ok(result.unwrap())
            }
        }
    }
    pub(crate) fn write<R>(&mut self, py: Python<'_>, f: impl FnOnce(&mut T) -> R) -> PyResult<R> {
        match self {
            Self::Owned(value) => Ok(f(value)),
            Self::View(view) => {
                let mut f = Some(f);
                let mut result = None;
                view.write(py, &mut |value| result = Some(f.take().unwrap()(value)))?;
                Ok(result.unwrap())
            }
        }
    }
}

pub(crate) trait NativeOwner:
    pyo3::PyClass<Frozen = pyo3::pyclass::boolean_struct::False>
{
    type Value;
    fn read_native(&self, py: Python<'_>, visit: &mut dyn FnMut(&Self::Value)) -> PyResult<()>;
    fn write_native(
        &mut self,
        py: Python<'_>,
        visit: &mut dyn FnMut(&mut Self::Value),
    ) -> PyResult<()>;
}
pub(crate) struct ChildField<P: NativeOwner, T> {
    parent: Py<P>,
    get: fn(&P::Value) -> &T,
    get_mut: fn(&mut P::Value) -> &mut T,
}
impl<P: NativeOwner, T> NativeView<T> for ChildField<P, T> {
    fn traverse(&self, visit: pyo3::gc::PyVisit<'_>) -> Result<(), pyo3::PyTraverseError> {
        visit.call(&self.parent)
    }
    fn read(&self, py: Python<'_>, visit: &mut dyn FnMut(&T)) -> PyResult<()> {
        self.parent
            .try_borrow(py)?
            .read_native(py, &mut |parent| visit((self.get)(parent)))
    }
    fn write(&self, py: Python<'_>, visit: &mut dyn FnMut(&mut T)) -> PyResult<()> {
        self.parent
            .try_borrow_mut(py)?
            .write_native(py, &mut |parent| visit((self.get_mut)(parent)))
    }
}

pub(crate) struct FieldOutput<T: Container, P: NativeOwner>(ChildField<P, T::Native>);
impl<T: Container, P: NativeOwner> FieldOutput<T, P> {
    pub(crate) fn new(
        parent: Py<P>,
        get: fn(&P::Value) -> &T::Native,
        get_mut: fn(&mut P::Value) -> &mut T::Native,
    ) -> Self {
        Self(ChildField {
            parent,
            get,
            get_mut,
        })
    }
}
impl<'py, T: Container, P: NativeOwner> IntoPyObject<'py> for FieldOutput<T, P> {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;
    const OUTPUT_TYPE: pyo3::inspect::PyStaticExpr = T::TYPE;
    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        T::output_field(py, self.0)
    }
}

pub(crate) trait ConfigurationView: Configuration {
    fn from_field(field: NativeField<Self::Native>) -> Self;
}

impl<T> Container for Vec<T>
where
    T: Clone + 'static,
    for<'a, 'py> T: FromPyObject<'a, 'py>,
    for<'py> T: IntoPyObject<'py>,
{
    type Native = Self;
    const TYPE: pyo3::inspect::PyStaticExpr = pyo3::type_hint_subscript!(
        pyo3::type_hint_identifier!("builtins", "tuple"),
        <T as IntoPyObject<'_>>::OUTPUT_TYPE,
        pyo3::inspect::PyStaticExpr::Constant {
            value: pyo3::inspect::PyStaticConstant::Ellipsis
        }
    );
    const INPUT: pyo3::inspect::PyStaticExpr = <Self as FromPyObject<'_, '_>>::INPUT_TYPE;
    fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<Self::Native> {
        value.extract()
    }

    fn output<'py>(py: Python<'py>, value: &Self::Native) -> PyResult<Bound<'py, PyAny>> {
        Ok(PyTuple::new(py, value.clone())?.into_any())
    }
}

// Fixed arrays use immutable tuples for reads and PyO3's native length/type
// conversion for assignment.
impl<T, const N: usize> Container for [T; N]
where
    T: Clone + 'static,
    for<'a, 'py> T: FromPyObject<'a, 'py>,
    for<'py> T: IntoPyObject<'py>,
{
    type Native = Self;
    const TYPE: pyo3::inspect::PyStaticExpr = <Vec<T> as Container>::TYPE;
    const INPUT: pyo3::inspect::PyStaticExpr = <Self as FromPyObject<'_, '_>>::INPUT_TYPE;
    fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<Self> {
        value.extract()
    }
    fn output<'py>(py: Python<'py>, value: &Self::Native) -> PyResult<Bound<'py, PyAny>> {
        Ok(PyTuple::new(py, value.clone())?.into_any())
    }
}

impl<T: Container> Container for Option<T> {
    const OPTIONAL: bool = true;
    type Native = Option<T::Native>;
    const TYPE: pyo3::inspect::PyStaticExpr = pyo3::type_hint_union!(
        T::TYPE,
        pyo3::inspect::PyStaticExpr::Constant {
            value: pyo3::inspect::PyStaticConstant::None
        }
    );
    const INPUT: pyo3::inspect::PyStaticExpr = pyo3::type_hint_union!(
        T::INPUT,
        pyo3::inspect::PyStaticExpr::Constant {
            value: pyo3::inspect::PyStaticConstant::None
        }
    );
    fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<Self::Native> {
        if value.is_none() {
            Ok(None)
        } else {
            T::extract(value).map(Some)
        }
    }

    fn output<'py>(py: Python<'py>, value: &Self::Native) -> PyResult<Bound<'py, PyAny>> {
        match value {
            Some(value) => T::output(py, value),
            None => Ok(py.None().into_bound(py)),
        }
    }
    fn output_field<'py, V: NativeView<Self::Native> + 'static>(
        py: Python<'py>,
        field: V,
    ) -> PyResult<Bound<'py, PyAny>> {
        let mut present = false;
        field.read(py, &mut |value| present = value.is_some())?;
        if present {
            T::output_field(
                py,
                OptionalField::<V, T::Native>(field, std::marker::PhantomData),
            )
        } else {
            Ok(py.None().into_bound(py))
        }
    }
}

struct OptionalField<V, T>(V, std::marker::PhantomData<fn() -> T>);
impl<V: NativeView<Option<T>>, T> NativeView<T> for OptionalField<V, T> {
    fn traverse(&self, visit: pyo3::gc::PyVisit<'_>) -> Result<(), pyo3::PyTraverseError> {
        self.0.traverse(visit)
    }
    fn read(&self, py: Python<'_>, visit: &mut dyn FnMut(&T)) -> PyResult<()> {
        let mut present = false;
        self.0.read(py, &mut |value| {
            if let Some(value) = value {
                present = true;
                visit(value);
            }
        })?;
        if present {
            Ok(())
        } else {
            Err(pyo3::exceptions::PyReferenceError::new_err(
                "configuration field is now None",
            ))
        }
    }
    fn write(&self, py: Python<'_>, visit: &mut dyn FnMut(&mut T)) -> PyResult<()> {
        let mut present = false;
        self.0.write(py, &mut |value| {
            if let Some(value) = value {
                present = true;
                visit(value);
            }
        })?;
        if present {
            Ok(())
        } else {
            Err(pyo3::exceptions::PyReferenceError::new_err(
                "configuration field is now None",
            ))
        }
    }
}

pub(crate) struct ContainerOutput<'a, T: Container>(pub(crate) &'a T::Native);
impl<'py, T: Container> IntoPyObject<'py> for ContainerOutput<'_, T> {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;
    const OUTPUT_TYPE: pyo3::inspect::PyStaticExpr = T::TYPE;
    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        T::output(py, self.0)
    }
}

/// Projection of a pure native input value, never an executor or diagnostics.
pub(crate) trait Configuration: pyo3::PyClass<BaseType = PyAny> {
    type Native: Clone + 'static;
    const INPUT: pyo3::inspect::PyStaticExpr = <Py<Self> as FromPyObject<'_, '_>>::INPUT_TYPE;
    fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<Self::Native> {
        let instance = value.extract::<Py<Self>>()?;
        let native = instance.try_borrow(value.py())?.to_native(value.py())?;
        Ok(native)
    }
    fn from_native(py: Python<'_>, value: Self::Native) -> PyResult<Self>;
    fn to_native(&self, py: Python<'_>) -> PyResult<Self::Native>;
}

pub(crate) struct ConfigurationField<T>(std::marker::PhantomData<T>);
pub(crate) struct LiveConfigurationField<T>(std::marker::PhantomData<T>);
impl<T: ConfigurationView> Container for LiveConfigurationField<T> {
    type Native = T::Native;
    const TYPE: pyo3::inspect::PyStaticExpr = <Py<T> as IntoPyObject<'_>>::OUTPUT_TYPE;
    const INPUT: pyo3::inspect::PyStaticExpr = T::INPUT;
    fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<Self::Native> {
        T::extract(value)
    }
    fn output<'py>(py: Python<'py>, value: &Self::Native) -> PyResult<Bound<'py, PyAny>> {
        ConfigurationField::<T>::output(py, value)
    }
    fn output_field<'py, V: NativeView<Self::Native> + 'static>(
        py: Python<'py>,
        field: V,
    ) -> PyResult<Bound<'py, PyAny>> {
        Ok(
            Py::new(py, T::from_field(NativeField::View(Box::new(field))))?
                .into_bound(py)
                .into_any(),
        )
    }
}
impl<T: Configuration> Container for ConfigurationField<T> {
    type Native = T::Native;
    const TYPE: pyo3::inspect::PyStaticExpr = <Py<T> as IntoPyObject<'_>>::OUTPUT_TYPE;
    const INPUT: pyo3::inspect::PyStaticExpr = T::INPUT;
    fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<Self::Native> {
        T::extract(value)
    }
    fn output<'py>(py: Python<'py>, value: &Self::Native) -> PyResult<Bound<'py, PyAny>> {
        Ok(Py::new(py, T::from_native(py, value.clone())?)?
            .into_bound(py)
            .into_any())
    }
}

/// A constructor/setter captures a value, not the caller's editable object.
pub(crate) struct ContainerInput<T: Container>(pub(crate) T::Native);
impl<'a, 'py, T: Container> FromPyObject<'a, 'py> for ContainerInput<T> {
    type Error = PyErr;
    const INPUT_TYPE: pyo3::inspect::PyStaticExpr = T::INPUT;
    fn extract(value: pyo3::Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        T::extract(value).map(Self)
    }
}
impl<'py, T: Container> IntoPyObject<'py> for ContainerInput<T> {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;
    const OUTPUT_TYPE: pyo3::inspect::PyStaticExpr = T::TYPE;
    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        T::output(py, &self.0)
    }
}
pub(crate) type ConfigurationInput<T> = ContainerInput<ConfigurationField<T>>;

/// Omission uses the field's native default. Explicit None is a value for a
/// genuine Option; over concrete storage it selects the documented default.
pub(crate) enum FieldInput<T: Container> {
    Default,
    Value(T::Native),
}
impl<T: Container> Default for FieldInput<T> {
    fn default() -> Self {
        Self::Default
    }
}
impl<T: Container> FieldInput<T> {
    pub(crate) fn or_default(self, default: impl FnOnce() -> T::Native) -> T::Native {
        match self {
            Self::Default => default(),
            Self::Value(value) => value,
        }
    }
}
impl<'a, 'py, T: Container> FromPyObject<'a, 'py> for FieldInput<T> {
    type Error = PyErr;
    const INPUT_TYPE: pyo3::inspect::PyStaticExpr = pyo3::type_hint_union!(
        T::INPUT,
        pyo3::inspect::PyStaticExpr::Constant {
            value: pyo3::inspect::PyStaticConstant::None
        }
    );
    fn extract(value: pyo3::Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        if value.is_none() && !T::OPTIONAL {
            Ok(Self::Default)
        } else {
            T::extract(value).map(Self::Value)
        }
    }
}
impl<'py, T: Container> IntoPyObject<'py> for FieldInput<T> {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;
    const OUTPUT_TYPE: pyo3::inspect::PyStaticExpr = T::TYPE;
    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        match self {
            Self::Default => Ok(py.None().into_bound(py)),
            Self::Value(value) => T::output(py, &value),
        }
    }
}

/// Readonly value leaves do not need a retained editable child object.
pub(crate) struct ValueField<T>(std::marker::PhantomData<T>);
impl<T> Container for ValueField<T>
where
    T: Configuration + pyo3::PyClass<Frozen = pyo3::pyclass::boolean_struct::True>,
    T::Native: Clone,
{
    type Native = T::Native;
    const TYPE: pyo3::inspect::PyStaticExpr = <Py<T> as IntoPyObject<'_>>::OUTPUT_TYPE;
    const INPUT: pyo3::inspect::PyStaticExpr = T::INPUT;
    fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<Self::Native> {
        T::extract(value)
    }
    fn output<'py>(py: Python<'py>, value: &Self::Native) -> PyResult<Bound<'py, PyAny>> {
        Ok(Py::new(py, T::from_native(py, value.clone())?)?
            .into_bound(py)
            .into_any())
    }
}

pub(crate) trait NativeScalar: Copy + 'static {
    type Repr: Copy;
    fn from_scalar(index: Self::Repr) -> Self;
    fn scalar(self) -> Self::Repr;
}
pub(crate) struct NativeScalarField<T>(std::marker::PhantomData<T>);
impl<T: NativeScalar> Container for NativeScalarField<T>
where
    for<'a, 'py> T::Repr: FromPyObject<'a, 'py>,
    for<'py> T::Repr: IntoPyObject<'py>,
{
    type Native = T;
    const TYPE: pyo3::inspect::PyStaticExpr = <T::Repr as IntoPyObject<'_>>::OUTPUT_TYPE;
    const INPUT: pyo3::inspect::PyStaticExpr = <T::Repr as FromPyObject<'_, '_>>::INPUT_TYPE;
    fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<T> {
        value
            .extract::<T::Repr>()
            .map(T::from_scalar)
            .map_err(Into::into)
    }
    fn output<'py>(py: Python<'py>, value: &T) -> PyResult<Bound<'py, PyAny>> {
        value
            .scalar()
            .into_pyobject(py)
            .map(|value| value.into_any().into_bound())
            .map_err(Into::into)
    }
}

/// Immutable scalar leaves inside a mutable collection (for example map keys
/// or coordinate tuples). This is ordinary PyO3 conversion, not validation.
pub(crate) struct Scalar<T>(std::marker::PhantomData<T>);
impl<T> Container for Scalar<T>
where
    T: Clone + 'static,
    for<'a, 'py> T: FromPyObject<'a, 'py>,
    for<'py> T: IntoPyObject<'py>,
{
    type Native = T;
    const TYPE: pyo3::inspect::PyStaticExpr = <T as IntoPyObject<'_>>::OUTPUT_TYPE;
    const INPUT: pyo3::inspect::PyStaticExpr = <T as FromPyObject<'_, '_>>::INPUT_TYPE;
    fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<T> {
        value.extract().map_err(Into::into)
    }
    fn output<'py>(py: Python<'py>, value: &T) -> PyResult<Bound<'py, PyAny>> {
        value
            .clone()
            .into_pyobject(py)
            .map(|value| value.into_any().into_bound())
            .map_err(Into::into)
    }
}

/// A native vector with recursively projected immutable reads.
pub(crate) struct List<T>(std::marker::PhantomData<T>);
impl<T: Container> Container for List<T> {
    type Native = Vec<T::Native>;
    const TYPE: pyo3::inspect::PyStaticExpr = pyo3::type_hint_subscript!(
        pyo3::type_hint_identifier!("builtins", "tuple"),
        T::TYPE,
        pyo3::inspect::PyStaticExpr::Constant {
            value: pyo3::inspect::PyStaticConstant::Ellipsis
        }
    );
    const INPUT: pyo3::inspect::PyStaticExpr = pyo3::type_hint_subscript!(
        pyo3::type_hint_identifier!("collections.abc", "Sequence"),
        T::INPUT
    );
    fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<Self::Native> {
        value.cast::<pyo3::types::PySequence>()?;
        value
            .try_iter()?
            .map(|value| T::extract(value?.as_borrowed()))
            .collect()
    }
    fn output<'py>(py: Python<'py>, value: &Self::Native) -> PyResult<Bound<'py, PyAny>> {
        let values = value
            .iter()
            .map(|value| T::output(py, value))
            .collect::<PyResult<Vec<_>>>()?;
        Ok(PyTuple::new(py, values)?.into_any())
    }
}

pub(crate) struct Dict<K, V>(std::marker::PhantomData<(K, V)>);
impl<K: Container, V: Container> Container for Dict<K, V>
where
    K::Native: Ord,
{
    type Native = BTreeMap<K::Native, V::Native>;
    const TYPE: pyo3::inspect::PyStaticExpr = pyo3::type_hint_subscript!(
        pyo3::type_hint_identifier!("collections.abc", "Mapping"),
        K::TYPE,
        V::TYPE
    );
    const INPUT: pyo3::inspect::PyStaticExpr = pyo3::type_hint_subscript!(
        pyo3::type_hint_identifier!("collections.abc", "Mapping"),
        K::INPUT,
        V::INPUT
    );
    fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<Self::Native> {
        let values = PyDict::new(value.py());
        let mapping = value.cast::<PyMapping>()?;
        values.update(&mapping)?;
        values
            .iter()
            .map(|(key, value)| {
                Ok((
                    K::extract(key.as_borrowed())?,
                    V::extract(value.as_borrowed())?,
                ))
            })
            .collect()
    }
    fn output<'py>(py: Python<'py>, value: &Self::Native) -> PyResult<Bound<'py, PyAny>> {
        let output = PyDict::new(py);
        for (key, value) in value {
            output.set_item(K::output(py, key)?, V::output(py, value)?)?;
        }
        Ok(PyMappingProxy::new(py, output.as_mapping()).into_any())
    }
}

pub(crate) trait NativeEnum:
    crate::chemical_values::EnumValue + pyo3::PyClass<BaseType = PyAny> + Copy
{
    type Native: Copy + 'static;
    const INPUT: pyo3::inspect::PyStaticExpr = pyo3::type_hint_union!(
        <Py<Self> as IntoPyObject<'_>>::OUTPUT_TYPE,
        pyo3::type_hint_identifier!("builtins", "str")
    );
    fn from_native(value: Self::Native) -> Self;
    fn into_native(self) -> Self::Native;
    fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<Self::Native> {
        if value.is_instance_of::<pyo3::types::PyString>() {
            Self::from_name(&value.extract::<String>()?).map(Self::into_native)
        } else {
            let object = value.extract::<Py<Self>>()?;
            let native = *object.try_borrow(value.py())?;
            Ok(native.into_native())
        }
    }
}
pub(crate) struct EnumField<T>(std::marker::PhantomData<T>);
impl<T: NativeEnum> Container for EnumField<T> {
    type Native = T::Native;
    const TYPE: pyo3::inspect::PyStaticExpr = <Py<T> as IntoPyObject<'_>>::OUTPUT_TYPE;
    const INPUT: pyo3::inspect::PyStaticExpr = T::INPUT;
    fn extract(value: pyo3::Borrowed<'_, '_, PyAny>) -> PyResult<Self::Native> {
        T::extract(value)
    }
    fn output<'py>(py: Python<'py>, value: &Self::Native) -> PyResult<Bound<'py, PyAny>> {
        Ok(Py::new(py, T::from_native(*value))?
            .into_bound(py)
            .into_any())
    }
}
