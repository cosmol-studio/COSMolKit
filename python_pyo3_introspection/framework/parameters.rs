//! A parameter declaration generates construction and writable fields.
//! It deliberately does not generate domain calls or a second default policy.
macro_rules! parameters {
    ($name:ident, $module:ident::$rust:ident, {$($field:ident: $ty:ty),+ $(,)?}) => {
        #[pyo3::pyclass(module = "fixed_projection._native")]
        pub struct $name { pub(crate) inner: $module::$rust }

        paste::paste! {
            #[pyo3::pymethods]
            impl $name {
                #[new]
                #[pyo3(signature = (*, $($field=$module::$rust::default().$field),+))]
                fn new($($field: $ty),+) -> Self {
                    // The exhaustive literal rejects new domain fields at compile time.
                    Self { inner: $module::$rust { $($field),+ } }
                }
                $(
                    #[getter]
                    fn $field(&self) -> $ty { self.inner.$field.clone() }
                    #[setter]
                    fn [<set_ $field>](&mut self, value: $ty) { self.inner.$field = value; }
                )+
            }
        }
    };
}
pub(crate) use parameters;
