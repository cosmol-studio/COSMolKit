//! Batch call forms share scalar keyword construction and one native query path.
//! Accumulate methods before applying PyO3, which owns one impl per class.
macro_rules! batch_fingerprint_methods {
    ([$($families:tt)*] impl $owner:ident { $($rest:tt)* }) => {
        batch_fingerprint_methods!(@families [$owner] [$($rest)*] [] ; $($families)*);
    };
    (@families [$owner:ident] [$($rest:tt)*] [$($methods:tt)*] ;) => {
        #[pyo3::pymethods]
        impl $owner { $($rest)* $($methods)* }
    };
    (@families [$owner:ident] [$($rest:tt)*] [$($methods:tt)*] ;
        $params:ty, $py:ident, $factory:expr,
        { $($fields:tt)* } { $($queries:tt)* } $($tail:tt)*) => {
        batch_fingerprint_methods!(@queries [$owner] [$($rest)*] [$($methods)*]
            [$($tail)*] [$params, $py, $factory, {$($fields)*}] ; $($queries)*);
    };
    (@queries [$owner:ident] [$($rest:tt)*] [$($methods:tt)*]
        [$($tail:tt)*] [$($family:tt)*] ;) => {
        batch_fingerprint_methods!(@families [$owner] [$($rest)*] [$($methods)*] ; $($tail)*);
    };
    (@queries [$owner:ident] [$($rest:tt)*] [$($methods:tt)*]
        [$($tail:tt)*]
        [$params:ty, $py:ident, $factory:expr,
         {$($field:ident: $ty:ty = $default:expr),* $(,)?}] ;
        $name:ident, $params_branch:ident, $keywords_branch:ident, $native:ident,
        $result:ty, $convert:expr, {$($extra:ident: $extra_ty:ty = $extra_default:expr),* $(,)?};
        $($queries:tt)*) => {
        batch_fingerprint_methods!(@queries [$owner] [$($rest)*] [
            $($methods)*
            #[pyo3(signature=(params=Argument::Missing, *, execution=None, n_jobs=None, progress_bar=None,
                $($extra=$extra_default,)* $($field=Argument::Missing),*))]
            fn $name(&self, $py: Python<'_>,
                params: Argument<Option<PyRef<'_, $params>>>,
                execution: Option<&crate::batch_query_params::BatchQueryParams>,
                n_jobs: Option<usize>, progress_bar: Option<bool>,
                $($extra: $extra_ty,)* $($field: Argument<$ty>),*
            ) -> PyResult<$result> {
                crate::arguments::configuration_arguments!(params; $($field),*);
                if let Some(params) = params.as_deref() {
                    self.$params_branch($py, Some(params), execution, n_jobs, progress_bar, $($extra),*)
                } else {
                    self.$keywords_branch($py, execution, n_jobs, progress_bar, $($extra,)* $($field),*)
                }
            }
            #[doc = concat!(r#"{"role":"overload","target":"MoleculeBatch."#, stringify!($name), r#""}"#)]
            #[pyo3(signature=(params=None, *, execution=None, n_jobs=None, progress_bar=None, $($extra=$extra_default),*))]
            fn $params_branch(&self, $py: Python<'_>, params: Option<&$params>,
                execution: Option<&crate::batch_query_params::BatchQueryParams>,
                n_jobs: Option<usize>, progress_bar: Option<bool>,
                $($extra: $extra_ty),*
            ) -> PyResult<$result> {
                let selected = params.map(|v| v.to_rust($py)).transpose()?.unwrap_or_default();
                crate::batch_calls::query($py, execution, n_jobs, progress_bar, |runtime| {
                    self.inner.$native(&selected, $($extra,)* runtime)
                }).map($convert)
            }
            #[doc = concat!(r#"{"role":"overload","target":"MoleculeBatch."#, stringify!($name), r#""}"#)]
            #[pyo3(signature=(*, execution=None, n_jobs=None, progress_bar=None, $($extra=$extra_default,)* $($field=$default),*))]
            fn $keywords_branch(&self, $py: Python<'_>,
                execution: Option<&crate::batch_query_params::BatchQueryParams>,
                n_jobs: Option<usize>, progress_bar: Option<bool>,
                $($extra: $extra_ty,)* $($field: $ty),*
            ) -> PyResult<$result> {
                let selected = $factory;
                crate::batch_calls::query($py, execution, n_jobs, progress_bar, |runtime| {
                    self.inner.$native(&selected, $($extra,)* runtime)
                }).map($convert)
            }
            fn $native(&self, $py: Python<'_>, options: &$params,
                $($extra: $extra_ty,)* execution: &crate::batch_query_params::BatchQueryParams
            ) -> PyResult<$result> {
                self.$params_branch($py, Some(options), Some(execution), None, None, $($extra),*)
            }
        ] [$($tail)*] [$params, $py, $factory, {$($field: $ty = $default),*}] ; $($queries)*);
    };
}

pub(crate) use batch_fingerprint_methods;
