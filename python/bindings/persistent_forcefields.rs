//! Owned evaluators forward the public facade; no force-field state is mirrored.
use crate::chemical_params::{ForceFieldMinimizeParams, MmffForceFieldParams, UffForceFieldParams};
use crate::error_framework::Source;
use crate::forcefield_execution_errors::MmffBuilderError;
use crate::io_source_errors::ErrorRef;
use crate::mmff_properties::MmffMolPropertiesError;
use crate::native::Molecule;
use crate::native_numerical_values::{ForceFieldEnergyGradient, ForceFieldMinimizeOutcome};
use crate::operation_errors::OperationError;
use crate::params::{
    ForceFieldAngleConstraintParams, ForceFieldDistanceConstraintParams,
    ForceFieldPositionConstraintParams, ForceFieldTorsionConstraintParams,
};
use crate::public_selectors::MolecularForceFieldErrorKind;
use crate::retained_error::RetainedError;
use crate::ring_errors::RingFindingError;
use crate::uff_parameter_errors::{ForceFieldKernelError, UffBuilderError, UffParameterError};
use cosmolkit as ck;
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

type FieldSource = Source<
    crate::error_framework::ErrorObject<MmffMolPropertiesError>,
    Source<
        crate::error_framework::ErrorObject<MmffBuilderError>,
        Source<
            crate::error_framework::ErrorObject<UffParameterError>,
            Source<
                crate::error_framework::ErrorObject<UffBuilderError>,
                crate::error_framework::ErrorObject<ForceFieldKernelError>,
            >,
        >,
    >,
>;

fn field_source(node: ErrorRef) -> PyErr {
    if let Some(error) = node.error().downcast_ref::<ck::MmffMolPropertiesError>() {
        return MmffMolPropertiesError::from_rust(error.clone());
    }
    if node.error().is::<ck::MmffBuilderError>() {
        return MmffBuilderError::from_ref(node);
    }
    if node.error().is::<ck::UffParameterError>() {
        return UffParameterError::from_ref(node);
    }
    if let Some(error) = node.error().downcast_ref::<ck::UffBuilderError>() {
        return UffBuilderError::from_rust(error.clone());
    }
    if let Some(error) = node.error().downcast_ref::<ck::ForceFieldKernelError>() {
        return ForceFieldKernelError::from_rust(*error);
    }
    unreachable!("complete native persistent force-field source graph")
}

// Both real record errors have the same readonly public accessors. Keep their
// actual native owner, not a duplicate enum/payload reconstructed from text.
macro_rules! field_error {
    ($name:ident, $source:ty, $convert:expr) => {
        crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
        pub struct $name {
            inner: RetainedError<ck::$name>,
        }
}
        impl $name {
            pub(crate) fn from_rust(inner: ck::$name) -> PyErr {
                Self::from_retained(RetainedError::new(inner))
            }
            pub(crate) fn from_ref(node: ErrorRef) -> PyErr {
                Self::from_retained(RetainedError::from_ref(node))
            }
            fn from_retained(inner: RetainedError<ck::$name>) -> PyErr {
                let cause = inner.node().source().map($convert);
                Python::attach(|py| match Py::new(py, Self { inner }) {
                    Ok(value) => {
                        let error = crate::error_framework::from_payload(value);
                        error.set_cause(py, cause);
                        error
                    }
                    Err(error) => error,
                })
            }
        }
        crate::error_framework::error_projection! {methods impl $name {
            #[getter]
            /// Native subsystem reporting this failure, preserved independently of the human-readable message.
            fn domain(&self) -> &'static str {
                "forcefield"
            }
            #[getter]
            /// Typed native persistent-force-field failure category.
            fn kind(&self) -> MolecularForceFieldErrorKind {
                self.inner.kind().into()
            }
            #[getter]
            /// Requested conformer identifier when relevant, otherwise None.
            fn requested(&self) -> Option<usize> {
                self.inner.requested()
            }
            #[getter]
            /// Atom index implicated in the failure, or None for failures without an atom.
            fn atom_index(&self) -> Option<usize> {
                self.inner.atom_index()
            }
            #[getter]
            /// Coordinate component implicated in the failure, or None when inapplicable.
            fn component(&self) -> Option<usize> {
                self.inner.component()
            }
            #[getter]
            /// Observed size or index reported by the native failure, or None when inapplicable.
            fn actual(&self) -> Option<usize> {
                self.inner.actual()
            }
            #[getter]
            /// Required size reported by the native failure, or None when inapplicable.
            fn expected(&self) -> Option<usize> {
                self.inner.expected()
            }
            #[staticmethod]
            #[doc = r#"{"role":"error_source"}"#]
            /// Actual underlying exception, identical to __cause__.
            fn source(slf: Bound<'_, pyo3::exceptions::PyBaseException>) -> PyResult<$source> {
                Ok(slf.getattr("__cause__")?.extract()?)
            }
            fn __str__(&self) -> String {
                self.inner.to_string()
            }
        }
}
    };
}
field_error!(ForceFieldError, Option<FieldSource>, field_source);
field_error!(
    MmffForceFieldError,
    crate::error_framework::ErrorObject<ForceFieldError>,
    ForceFieldError::from_ref
);

type UffFieldSource = Source<
    crate::error_framework::ErrorObject<ForceFieldError>,
    Source<
        crate::error_framework::ErrorObject<RingFindingError>,
        crate::error_framework::ErrorObject<OperationError>,
    >,
>;

crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
pub struct UffForceFieldError {
    inner: RetainedError<ck::UffForceFieldError>,
}
}
impl UffForceFieldError {
    pub(crate) fn from_rust(inner: ck::UffForceFieldError) -> PyErr {
        let inner = RetainedError::new(inner);
        let cause = inner.node().source().map(|node| {
            if node.error().is::<ck::ForceFieldError>() {
                ForceFieldError::from_ref(node)
            } else if let Some(error) = node.error().downcast_ref::<ck::RingFindingError>() {
                RingFindingError::from_rust(error.clone())
            } else if node.error().is::<ck::OperationError>() {
                OperationError::from_ref(node)
            } else {
                unreachable!("complete native persistent UFF factory source graph")
            }
        });
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => {
                let error = crate::error_framework::from_payload(value);
                error.set_cause(py, cause);
                error
            }
            Err(error) => error,
        })
    }
}
crate::error_framework::error_projection! {methods impl UffForceFieldError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "forcefield"
    }
    #[getter]
    /// Typed native persistent-force-field failure category.
    fn kind(&self) -> MolecularForceFieldErrorKind {
        self.inner.kind().into()
    }
    #[getter]
    /// Requested conformer identifier when relevant, otherwise None.
    fn requested(&self) -> Option<usize> {
        self.inner.requested()
    }
    #[staticmethod]
    #[doc = r#"{"role":"error_source"}"#]
    /// Actual underlying exception, identical to __cause__.
    fn source(slf: Bound<'_, pyo3::exceptions::PyBaseException>) -> PyResult<Option<UffFieldSource>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}

#[pyclass(module = "cosmolkit")]
/// Persistent native force-field evaluator owning positions, constraints and minimization state.
/// Position edits and minimization affect this evaluator, not the molecule used to construct it.
pub struct MolecularForceField {
    pub(crate) inner: ck::MolecularForceField,
}

// A single pymethods block is emitted; no PyO3 inventory/dependency feature is
// needed merely to group four repetitive argument-forwarding signatures.
macro_rules! evaluator_methods {
    (constraints [$(($method:ident, $params_branch:ident, $keywords_branch:ident, $target:literal,
       $params:ident, [$($atom:ident),+], {$($field:ident: $ty:ty),+})),+ $(,)?]
     methods {$($ordinary:item)*}) => {
        #[pymethods]
        impl MolecularForceField {
            $($ordinary)*
            $(
            #[doc = concat!("Add the constraint specified by ", stringify!($params),
                " to this force-field evaluator in place. Atom arguments are zero-based indices. ",
                "Supply the parameter object or the required keyword fields, not both; ",
                "the source molecule is unchanged.")]
            #[pyo3(signature=($($atom,)+ params=None, *, $($field=None),+))]
            fn $method(&mut self, $($atom: usize,)+ params: Option<&$params>, $($field: Option<$ty>),+) -> PyResult<()> {
                let configured = false $(|| $field.is_some())+;
                if params.is_some() && configured { return Err(mixed()); }
                if let Some(params) = params { self.$params_branch($($atom,)+ params) }
                else { self.$keywords_branch($($atom,)+ $(required($field, stringify!($field))?),+) }
            }
            #[doc = $target]
            fn $params_branch(&mut self, $($atom: usize,)+ params: &$params) -> PyResult<()> {
                let params = Python::attach(|py| params.to_rust(py))?;
                self.inner.$method($(ck::AtomId::new($atom),)+ &params).map_err(ForceFieldError::from_rust)
            }
            #[doc = $target]
            #[pyo3(signature=($($atom,)+ *, $($field),+))]
            fn $keywords_branch(&mut self, $($atom: usize,)+ $($field: $ty),+) -> PyResult<()> {
                let params = $params { inner: crate::parameters::NativeField::Owned(ck::$params::new($($field),+)) };
                self.$params_branch($($atom,)+ &params)
            }
            )+
        }
    };
}
evaluator_methods! {
    constraints [
        (add_distance_constraint_, _distance_params, _distance_keywords,
         r#"{"role":"overload","target":"MolecularForceField.add_distance_constraint_"}"#,
         ForceFieldDistanceConstraintParams, [atom1, atom2], {relative: bool, min_distance: f64, max_distance: f64, force_constant: f64}),
        (add_angle_constraint_, _angle_params, _angle_keywords,
         r#"{"role":"overload","target":"MolecularForceField.add_angle_constraint_"}"#,
         ForceFieldAngleConstraintParams, [atom1, atom2, atom3], {relative: bool, min_angle_deg: f64, max_angle_deg: f64, force_constant: f64}),
        (add_torsion_constraint_, _torsion_params, _torsion_keywords,
         r#"{"role":"overload","target":"MolecularForceField.add_torsion_constraint_"}"#,
         ForceFieldTorsionConstraintParams, [atom1, atom2, atom3, atom4], {relative: bool, min_dihedral_deg: f64, max_dihedral_deg: f64, force_constant: f64}),
        (add_position_constraint_, _position_params, _position_keywords,
         r#"{"role":"overload","target":"MolecularForceField.add_position_constraint_"}"#,
         ForceFieldPositionConstraintParams, [atom], {max_displacement: f64, force_constant: f64})
    ]
    methods {
    /// Return the force field's current Cartesian position for the specified atom.
    fn position(&self, atom_id: usize) -> PyResult<[f64; 3]> {
        self.inner
            .position(ck::AtomId::new(atom_id))
            .map_err(ForceFieldError::from_rust)
    }
    /// Update one atom's position in this force-field evaluator, not in the source molecule.
    fn set_position_(&mut self, atom_id: usize, position: [f64; 3]) -> PyResult<()> {
        self.inner
            .set_position_(ck::AtomId::new(atom_id), position)
            .map_err(ForceFieldError::from_rust)
    }
    #[getter]
    /// Independent list of the evaluator's current Cartesian position triples.
    fn positions(&self) -> Vec<[f64; 3]> {
        self.inner.positions()
    }
    /// Replace the evaluator's positions from one triple per atom; native shape and value errors propagate.
    fn set_positions_(&mut self, positions: Vec<[f64; 3]>) -> PyResult<()> {
        self.inner
            .set_positions_(&positions)
            .map_err(ForceFieldError::from_rust)
    }
    #[getter]
    /// Atom indices currently held fixed during minimization.
    fn fixed_atoms(&self) -> Vec<usize> {
        self.inner
            .fixed_atoms()
            .into_iter()
            .map(ck::AtomId::index)
            .collect()
    }
    /// Replace the set of atoms held fixed during minimization.
    fn set_fixed_atoms_(&mut self, atom_ids: Vec<usize>) -> PyResult<()> {
        self.inner
            .set_fixed_atoms_(
                &atom_ids
                    .into_iter()
                    .map(ck::AtomId::new)
                    .collect::<Vec<_>>(),
            )
            .map_err(ForceFieldError::from_rust)
    }
    /// Evaluate the native force-field energy at the evaluator's current positions.
    fn energy(&self, py: Python<'_>) -> PyResult<f64> {
        py.detach(|| self.inner.energy())
            .map_err(ForceFieldError::from_rust)
    }
    /// Evaluate the Cartesian gradient with the native fixed-atom treatment.
    fn gradient(&self, py: Python<'_>) -> PyResult<Vec<[f64; 3]>> {
        py.detach(|| self.inner.gradient())
            .map_err(ForceFieldError::from_rust)
    }
    /// Evaluate the Cartesian gradient without fixing gradient components for fixed atoms.
    fn gradient_unconstrained(&self, py: Python<'_>) -> PyResult<Vec<[f64; 3]>> {
        py.detach(|| self.inner.gradient_unconstrained())
            .map_err(ForceFieldError::from_rust)
    }
    /// Evaluate energy and per-atom Cartesian gradients together.
    fn energy_gradient(&self, py: Python<'_>) -> PyResult<ForceFieldEnergyGradient> {
        py.detach(|| self.inner.energy_gradient())
            .map(|inner| ForceFieldEnergyGradient { inner })
            .map_err(ForceFieldError::from_rust)
    }
#[pyo3(signature=(params=crate::arguments::Argument::Missing, *, max_iterations=crate::arguments::Argument::Missing, force_tolerance=crate::arguments::Argument::Missing, energy_tolerance=crate::arguments::Argument::Missing))]
    /// Minimize the evaluator's positions in place and return convergence, iteration and energy information.
    fn minimize_(
        &mut self,
        py: Python<'_>,
        params: crate::arguments::Argument<Option<PyRef<'_, ForceFieldMinimizeParams>>>,
        max_iterations: crate::arguments::Argument<Option<u32>>,
        force_tolerance: crate::arguments::Argument<Option<f64>>,
        energy_tolerance: crate::arguments::Argument<Option<f64>>,
    ) -> PyResult<ForceFieldMinimizeOutcome> {
        let keyword_form = max_iterations.is_supplied() || force_tolerance.is_supplied() || energy_tolerance.is_supplied();
        crate::arguments::configuration_arguments!(params; max_iterations, force_tolerance, energy_tolerance);
        let params = params.as_deref();
        if keyword_form { self._minimize_keywords(py, max_iterations, force_tolerance, energy_tolerance) } else { self._minimize_params(py, params) }
    }
    #[doc = r#"{"role":"overload","target":"MolecularForceField.minimize_"}"#]
    #[pyo3(signature=(params=None))]
    fn _minimize_params(
        &mut self,
        py: Python<'_>,
        params: Option<&ForceFieldMinimizeParams>,
    ) -> PyResult<ForceFieldMinimizeOutcome> {
let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

        py.detach(|| match params {
            Some(params) => self.inner.minimize_with_params_(&params),
            None => self.inner.minimize_(),
        })
        .map(|inner| ForceFieldMinimizeOutcome { inner })
        .map_err(ForceFieldError::from_rust)
    }
    #[doc = r#"{"role":"overload","target":"MolecularForceField.minimize_"}"#]
    #[pyo3(signature=(*, max_iterations=None, force_tolerance=None, energy_tolerance=None))]
    fn _minimize_keywords(
        &mut self,
        py: Python<'_>,
        max_iterations: Option<u32>,
        force_tolerance: Option<f64>,
        energy_tolerance: Option<f64>,
    ) -> PyResult<ForceFieldMinimizeOutcome> {
        let defaults = ck::ForceFieldMinimizeParams::default();
        let params = ForceFieldMinimizeParams {
            inner: crate::parameters::NativeField::Owned(ck::ForceFieldMinimizeParams::new(
                max_iterations.unwrap_or(defaults.max_iterations()),
                force_tolerance.unwrap_or(defaults.force_tolerance()),
                energy_tolerance.unwrap_or(defaults.energy_tolerance()),
            )),
        };
        self._minimize_params(py, Some(&params))
    }
    }
}

fn mixed() -> PyErr {
    PyTypeError::new_err("Params and keyword configuration are mutually exclusive")
}
fn required<T>(value: Option<T>, name: &'static str) -> PyResult<T> {
    value.ok_or_else(|| PyTypeError::new_err(format!("missing required keyword argument: {name}")))
}

impl Molecule {
    pub(crate) fn uff_field_params(
        &self,
        py: Python<'_>,
        params: Option<&UffForceFieldParams>,
    ) -> PyResult<MolecularForceField> {
        let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

        py.detach(|| match params {
            Some(params) => self.inner.uff_force_field_with_params(&params),
            None => self.inner.uff_force_field(),
        })
        .map(|inner| MolecularForceField { inner })
        .map_err(UffForceFieldError::from_rust)
    }
    pub(crate) fn uff_field_keywords(
        &self,
        py: Python<'_>,
        conformer_id: Option<usize>,
        vdw_threshold: Option<f64>,
        ignore_interfragment_interactions: Option<bool>,
    ) -> PyResult<MolecularForceField> {
        let defaults = ck::UffForceFieldParams::default();
        let params = UffForceFieldParams {
            inner: crate::parameters::NativeField::Owned(ck::UffForceFieldParams::new(
                conformer_id.or(defaults.conformer_id()),
                vdw_threshold.unwrap_or(defaults.vdw_threshold()),
                ignore_interfragment_interactions
                    .unwrap_or(defaults.ignore_interfragment_interactions()),
            )),
        };
        self.uff_field_params(py, Some(&params))
    }
    pub(crate) fn mmff_field_params(
        &self,
        py: Python<'_>,
        params: Option<&MmffForceFieldParams>,
    ) -> PyResult<MolecularForceField> {
        let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

        py.detach(|| match params {
            Some(params) => self.inner.mmff_force_field_with_params(&params),
            None => self.inner.mmff_force_field(),
        })
        .map(|inner| MolecularForceField { inner })
        .map_err(MmffForceFieldError::from_rust)
    }
    pub(crate) fn mmff_field_keywords(
        &self,
        py: Python<'_>,
        conformer_id: Option<usize>,
        mmff_variant: Option<crate::leaf_values::MmffVariantInput>,
        non_bonded_threshold: Option<f64>,
        ignore_interfragment_interactions: Option<bool>,
    ) -> PyResult<MolecularForceField> {
        let defaults = ck::MmffForceFieldParams::default();
        let params = MmffForceFieldParams {
            inner: crate::parameters::NativeField::Owned(ck::MmffForceFieldParams::new(
                conformer_id.or(defaults.conformer_id()),
                mmff_variant
                    .map(|value| value.0)
                    .unwrap_or(defaults.mmff_variant()),
                non_bonded_threshold.unwrap_or(defaults.non_bonded_threshold()),
                ignore_interfragment_interactions
                    .unwrap_or(defaults.ignore_interfragment_interactions()),
            )),
        };
        self.mmff_field_params(py, Some(&params))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::types::PyDict;
    use std::error::Error;

    fn bits(rows: &[[f64; 3]]) -> Vec<u64> {
        rows.iter().flatten().map(|value| value.to_bits()).collect()
    }

    #[test]
    fn persistent_uff_factory_preserves_native_options_preparation_causes_and_owned_lifetime() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let molecule = ck::Molecule::from_smiles("CCO").unwrap();
            let mut builder = molecule.to_builder();
            builder
                .add_3d_conformer(vec![[0., 0., 0.], [1.7, 0.1, 0.], [2.3, 1.1, 0.2]])
                .unwrap();
            let molecule = builder.build().unwrap().with_assigned_valence().unwrap();
            let before = molecule.to_binary().unwrap();
            let value = Py::new(
                py,
                Molecule {
                    inner: molecule.clone(),
                },
            )?;
            let mol = value.bind(py);
            let defaults = mol.call_method0("uff_force_field")?;
            same_field(py, &defaults, &molecule.uff_force_field().unwrap())?;
            let options = ck::UffForceFieldParams::new(Some(0), 5., false);
            let params = Py::new(
                py,
                UffForceFieldParams {
                    inner: crate::parameters::NativeField::Owned(options.clone()),
                },
            )?;
            let configured = mol.call_method1("uff_force_field", (&params,))?;
            same_field(
                py,
                &configured,
                &molecule.uff_force_field_with_params(&options).unwrap(),
            )?;
            let kw = PyDict::new(py);
            kw.set_item("conformer_id", 0)?;
            kw.set_item("vdw_threshold", 5.)?;
            kw.set_item("ignore_interfragment_interactions", false)?;
            let configured_kw = mol.call_method("uff_force_field", (), Some(&kw))?;
            same_field(
                py,
                &configured_kw,
                &molecule.uff_force_field_with_params(&options).unwrap(),
            )?;
            assert!(
                mol.call_method("uff_force_field", (&params,), Some(&kw))
                    .unwrap_err()
                    .is_instance_of::<PyTypeError>(py)
            );
            let missing_options = ck::UffForceFieldParams::new(Some(usize::MAX), 5., false);
            let expected = molecule
                .uff_force_field_with_params(&missing_options)
                .unwrap_err();
            let kw = PyDict::new(py);
            kw.set_item("conformer_id", usize::MAX)?;
            let missing = mol
                .call_method("uff_force_field", (), Some(&kw))
                .unwrap_err();
            assert!(missing.is_instance(py, <UffForceFieldError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert_eq!(missing.value(py).str()?.to_str()?, expected.to_string());
            assert_eq!(
                missing.value(py).getattr("requested")?.extract::<usize>()?,
                usize::MAX
            );
            assert!(missing.cause(py).is_none());
            assert!(missing.value(py).setattr("requested", 0).is_err());

            let parse = ck::SmilesParseParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            };
            let unprepared = ck::Molecule::from_smiles_with_params("CCO", &parse).unwrap();
            let mut builder = unprepared.to_builder();
            builder
                .add_3d_conformer(vec![[0., 0., 0.], [1.7, 0.1, 0.], [2.3, 1.1, 0.2]])
                .unwrap();
            let unprepared = builder.build().unwrap();
            let expected = unprepared.uff_force_field().unwrap_err();
            assert!(expected.source().unwrap().is::<ck::OperationError>());
            let text = expected.to_string();
            let unprepared = Py::new(py, Molecule { inner: unprepared })?;
            let failed = unprepared
                .bind(py)
                .call_method0("uff_force_field")
                .unwrap_err();
            assert_eq!(failed.value(py).str()?.to_str()?, text);
            let cause = failed.cause(py).unwrap();
            assert!(cause.is_instance(py, <OperationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert!(failed.value(py).getattr("source")?.is(cause.value(py)));
            drop(failed);
            drop(unprepared);
            assert_eq!(
                cause.value(py).getattr("state")?.extract::<String>()?,
                "valence"
            );
            assert_eq!(
                cause.value(py).getattr("field")?.extract::<String>()?,
                "assignment"
            );
            assert_eq!(cause.value(py).getattr("actual")?.extract::<usize>()?, 0);
            assert_eq!(cause.value(py).getattr("expected")?.extract::<usize>()?, 1);
            assert!(cause.cause(py).is_none());
            assert_eq!(molecule.to_binary().unwrap(), before);
            drop(value);
            same_field(py, &defaults, &molecule.uff_force_field().unwrap())?;
            Ok(())
        })
        .unwrap();
    }
    fn same_field(
        py: Python<'_>,
        projected: &Bound<'_, PyAny>,
        native: &ck::MolecularForceField,
    ) -> PyResult<()> {
        assert_eq!(
            bits(&projected.getattr("positions")?.extract::<Vec<[f64; 3]>>()?),
            bits(&native.positions())
        );
        assert_eq!(
            projected
                .call_method0("energy")?
                .extract::<f64>()?
                .to_bits(),
            native.energy().unwrap().to_bits()
        );
        assert_eq!(
            bits(
                &projected
                    .call_method0("gradient")?
                    .extract::<Vec<[f64; 3]>>()?
            ),
            bits(&native.gradient().unwrap())
        );
        assert_eq!(
            bits(
                &projected
                    .call_method0("gradient_unconstrained")?
                    .extract::<Vec<[f64; 3]>>()?
            ),
            bits(&native.gradient_unconstrained().unwrap())
        );
        let pair = projected.call_method0("energy_gradient")?;
        let expected = native.energy_gradient().unwrap();
        assert_eq!(
            pair.getattr("energy")?.extract::<f64>()?.to_bits(),
            expected.energy().to_bits()
        );
        assert_eq!(
            bits(&pair.getattr("gradient")?.extract::<Vec<[f64; 3]>>()?),
            bits(expected.gradient())
        );
        assert!(pair.setattr("energy", 0.).is_err());
        let _ = py;
        Ok(())
    }

    #[test]
    fn persistent_mmff_forwards_owned_state_constraints_minimization_and_concrete_causes() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let molecule = ck::Molecule::from_smiles("CCCC").unwrap();
            let mut builder = molecule.to_builder();
            builder
                .add_3d_conformer(vec![
                    [0., 0., 0.],
                    [1.7, 0.1, 0.],
                    [2.3, 1.1, 0.2],
                    [3.4, 1.5, 0.9],
                ])
                .unwrap();
            let molecule = builder.build().unwrap().with_assigned_valence().unwrap();
            let before = molecule.to_binary().unwrap();
            let projected = Py::new(
                py,
                Molecule {
                    inner: molecule.clone(),
                },
            )?;
            let mol = projected.bind(py);
            let mut native = molecule.mmff_force_field().unwrap();
            let field = mol.call_method0("mmff_force_field")?;
            same_field(py, &field, &native)?;
            let selection =
                ck::MmffForceFieldParams::new(Some(0), ck::MmffVariant::Mmff94s, 20., false);
            let params = Py::new(
                py,
                MmffForceFieldParams {
                    inner: crate::parameters::NativeField::Owned(selection.clone()),
                },
            )?;
            let configured = mol.call_method1("mmff_force_field", (&params,))?;
            same_field(
                py,
                &configured,
                &molecule.mmff_force_field_with_params(&selection).unwrap(),
            )?;
            let kwargs = PyDict::new(py);
            kwargs.set_item("conformer_id", 0)?;
            kwargs.set_item("mmff_variant", "MMFF94s")?;
            kwargs.set_item("non_bonded_threshold", 20.)?;
            kwargs.set_item("ignore_interfragment_interactions", false)?;
            let keywords = mol.call_method("mmff_force_field", (), Some(&kwargs))?;
            same_field(
                py,
                &keywords,
                &molecule.mmff_force_field_with_params(&selection).unwrap(),
            )?;
            assert!(
                mol.call_method("mmff_force_field", (&params,), Some(&kwargs))
                    .unwrap_err()
                    .is_instance_of::<PyTypeError>(py)
            );

            let snapshot = field
                .getattr("positions")?
                .extract::<Vec<[f64; 3]>>()?;
            field.call_method1("set_position_", (0, [0.1, 0.2, 0.3]))?;
            native
                .set_position_(ck::AtomId::new(0), [0.1, 0.2, 0.3])
                .unwrap();
            assert_eq!(snapshot[0], [0., 0., 0.]);
            assert_eq!(
                field
                    .call_method1("position", (0,))?
                    .extract::<[f64; 3]>()?,
                [0.1, 0.2, 0.3]
            );
            let rows = vec![
                [0.2, 0., 0.],
                [1.8, 0.2, 0.],
                [2.4, 1.2, 0.2],
                [3.5, 1.6, 0.9],
            ];
            field.call_method1("set_positions_", (rows.clone(),))?;
            native.set_positions_(&rows).unwrap();
            field.call_method1("set_fixed_atoms_", (vec![2, 0, 2],))?;
            native
                .set_fixed_atoms_(&[ck::AtomId::new(2), ck::AtomId::new(0), ck::AtomId::new(2)])
                .unwrap();
            assert_eq!(
                field.getattr("fixed_atoms")?.extract::<Vec<usize>>()?,
                vec![0, 2]
            );
            same_field(py, &field, &native)?;

            // Every constraint's Params and keyword branch executes a real
            // native append, in identical order, against the same owned rows.
            for direct in [false, true] {
                let params = Py::new(
                    py,
                    ForceFieldDistanceConstraintParams {
                        inner: crate::parameters::NativeField::Owned(
                            ck::ForceFieldDistanceConstraintParams::new(true, -0.1, 0.1, 3.),
                        ),
                    },
                )?;
                if direct {
                    let kw = PyDict::new(py);
                    kw.set_item("relative", true)?;
                    kw.set_item("min_distance", -0.1)?;
                    kw.set_item("max_distance", 0.1)?;
                    kw.set_item("force_constant", 3.)?;
                    field.call_method("add_distance_constraint_", (0, 3), Some(&kw))?;
                    assert!(
                        field
                            .call_method("add_distance_constraint_", (0, 3, &params), Some(&kw))
                            .unwrap_err()
                            .is_instance_of::<PyTypeError>(py)
                    );
                } else {
                    field.call_method1("add_distance_constraint_", (0, 3, &params))?;
                }
                native
                    .add_distance_constraint_(
                        ck::AtomId::new(0),
                        ck::AtomId::new(3),
                        &ck::ForceFieldDistanceConstraintParams::new(true, -0.1, 0.1, 3.),
                    )
                    .unwrap();
                let params = Py::new(
                    py,
                    ForceFieldAngleConstraintParams {
                        inner: crate::parameters::NativeField::Owned(
                            ck::ForceFieldAngleConstraintParams::new(false, 90., 120., 2.),
                        ),
                    },
                )?;
                if direct {
                    let kw = PyDict::new(py);
                    kw.set_item("relative", false)?;
                    kw.set_item("min_angle_deg", 90.)?;
                    kw.set_item("max_angle_deg", 120.)?;
                    kw.set_item("force_constant", 2.)?;
                    field.call_method("add_angle_constraint_", (0, 1, 2), Some(&kw))?;
                } else {
                    field.call_method1("add_angle_constraint_", (0, 1, 2, &params))?;
                }
                native
                    .add_angle_constraint_(
                        ck::AtomId::new(0),
                        ck::AtomId::new(1),
                        ck::AtomId::new(2),
                        &ck::ForceFieldAngleConstraintParams::new(false, 90., 120., 2.),
                    )
                    .unwrap();
                let params = Py::new(
                    py,
                    ForceFieldTorsionConstraintParams {
                        inner: crate::parameters::NativeField::Owned(
                            ck::ForceFieldTorsionConstraintParams::new(false, -30., 30., 1.),
                        ),
                    },
                )?;
                if direct {
                    let kw = PyDict::new(py);
                    kw.set_item("relative", false)?;
                    kw.set_item("min_dihedral_deg", -30.)?;
                    kw.set_item("max_dihedral_deg", 30.)?;
                    kw.set_item("force_constant", 1.)?;
                    field.call_method("add_torsion_constraint_", (0, 1, 2, 3), Some(&kw))?;
                } else {
                    field.call_method1("add_torsion_constraint_", (0, 1, 2, 3, &params))?;
                }
                native
                    .add_torsion_constraint_(
                        ck::AtomId::new(0),
                        ck::AtomId::new(1),
                        ck::AtomId::new(2),
                        ck::AtomId::new(3),
                        &ck::ForceFieldTorsionConstraintParams::new(false, -30., 30., 1.),
                    )
                    .unwrap();
                let params = Py::new(
                    py,
                    ForceFieldPositionConstraintParams {
                        inner: crate::parameters::NativeField::Owned(
                            ck::ForceFieldPositionConstraintParams::new(0.2, 4.),
                        ),
                    },
                )?;
                if direct {
                    let kw = PyDict::new(py);
                    kw.set_item("max_displacement", 0.2)?;
                    kw.set_item("force_constant", 4.)?;
                    field.call_method("add_position_constraint_", (1,), Some(&kw))?;
                } else {
                    field.call_method1("add_position_constraint_", (1, &params))?;
                }
                native
                    .add_position_constraint_(
                        ck::AtomId::new(1),
                        &ck::ForceFieldPositionConstraintParams::new(0.2, 4.),
                    )
                    .unwrap();
                same_field(py, &field, &native)?;
            }
            assert!(
                field
                    .call_method1("add_distance_constraint_", (0, 1))
                    .unwrap_err()
                    .is_instance_of::<PyTypeError>(py)
            );
            for direct in [false, true] {
                let effective = ck::ForceFieldMinimizeParams::new(2, 1e-4, 1e-6);
                let outcome = if direct {
                    let kw = PyDict::new(py);
                    kw.set_item("max_iterations", 2)?;
                    kw.set_item("force_tolerance", 1e-4)?;
                    kw.set_item("energy_tolerance", 1e-6)?;
                    field.call_method("minimize_", (), Some(&kw))?
                } else {
                    let params = Py::new(
                        py,
                        ForceFieldMinimizeParams {
                            inner: crate::parameters::NativeField::Owned(effective),
                        },
                    )?;
                    field.call_method1("minimize_", (params,))?
                };
                let expected = native
                    .minimize_with_params_(&ck::ForceFieldMinimizeParams::new(2, 1e-4, 1e-6))
                    .unwrap();
                assert_eq!(
                    outcome.getattr("converged")?.extract::<bool>()?,
                    expected.converged()
                );
                assert_eq!(
                    outcome.getattr("iterations")?.extract::<u32>()?,
                    expected.iterations()
                );
                assert_eq!(
                    outcome.getattr("energy")?.extract::<f64>()?.to_bits(),
                    expected.energy().to_bits()
                );
                same_field(py, &field, &native)?;
            }
            let expected = native.position(ck::AtomId::new(usize::MAX)).unwrap_err();
            let error = field.call_method1("position", (usize::MAX,)).unwrap_err();
            assert!(error.is_instance(py, <ForceFieldError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert_eq!(error.value(py).str()?.to_str()?, expected.to_string());
            assert_eq!(
                error.value(py).getattr("atom_index")?.extract::<usize>()?,
                usize::MAX
            );
            assert!(error.value(py).setattr("kind", "changed").is_err());
            assert!(error.cause(py).is_none());
            let saved = field
                .getattr("positions")?
                .extract::<Vec<[f64; 3]>>()?;
            let failure = field
                .call_method1("set_position_", (1, [f64::NAN, 0., 0.]))
                .unwrap_err();
            assert_eq!(
                failure.value(py).getattr("component")?.extract::<usize>()?,
                0
            );
            assert_eq!(
                bits(
                    &field
                        .getattr("positions")?
                        .extract::<Vec<[f64; 3]>>()?
                ),
                bits(&saved)
            );
            let kw = PyDict::new(py);
            kw.set_item("force_tolerance", 0.)?;
            let failure = field.call_method("minimize_", (), Some(&kw)).unwrap_err();
            assert!(
                failure
                    .cause(py)
                    .unwrap()
                    .is_instance(py, <ForceFieldKernelError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())
            );
            let bare = ck::Molecule::from_smiles("CC").unwrap();
            let expected = bare.mmff_force_field().unwrap_err().to_string();
            let bare = Py::new(py, Molecule { inner: bare })?;
            let missing = bare.bind(py).call_method0("mmff_force_field").unwrap_err();
            assert!(missing.is_instance(py, <MmffForceFieldError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert_eq!(missing.value(py).str()?.to_str()?, expected);
            let child = missing.cause(py).unwrap();
            drop(missing);
            drop(bare);
            assert!(child.is_instance(py, <ForceFieldError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert_eq!(child.value(py).str()?.to_str()?, expected);
            assert!(child.cause(py).is_none());
            assert_eq!(molecule.to_binary().unwrap(), before);
            drop(projected); // The evaluator independently owns native state.
            same_field(py, &field, &native)?;
            Ok(())
        })
        .unwrap();
    }
}
