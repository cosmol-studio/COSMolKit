use crate::chemical_values::{EnumValue, enum_value};
use crate::native::Molecule;
use cosmolkit as ck;
use pyo3::prelude::*;

enum_value!(ReactionRole, { Reactant, Product, Agent });
enum_value!(ReactionValidationSeverity, { Warning, Error });
enum_value!(ReactionValidationIssueKind, {
    MissingReactants, MissingProducts, DuplicateReactantMap, UnmappedReactant,
    DuplicateProductMap, MissingProductMapReactant, UnmappedProduct, UnmappedReactantMaps,
    MultipleCharge, MultipleHydrogenCount, MultipleMass, MultipleIsotope
});
#[derive(FromPyObject, IntoPyObject)]
pub enum RoleInput {
    Value(ReactionRole),
    Name(String),
}
#[derive(FromPyObject, IntoPyObject)]
pub enum SeverityInput {
    Value(ReactionValidationSeverity),
    Name(String),
}
#[derive(FromPyObject, IntoPyObject)]
pub enum IssueKindInput {
    Value(ReactionValidationIssueKind),
    Name(String),
}
macro_rules! input {
    ($input:ident, $value:ident) => {
        impl $input {
            fn into_rust(self) -> PyResult<ck::$value> {
                match self {
                    Self::Value(value) => Ok(value.into()),
                    Self::Name(name) => $value::from_name(&name).map(Into::into),
                }
            }
        }
    };
}
input!(RoleInput, ReactionRole);
input!(SeverityInput, ReactionValidationSeverity);
input!(IssueKindInput, ReactionValidationIssueKind);

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct ReactionValidationIssue {
    pub(crate) inner: ck::ReactionValidationIssue,
}
#[pymethods]
impl ReactionValidationIssue {
    #[new]
    #[pyo3(signature=(kind, severity, role, template, atom, map, maps, detail))]
    fn new(
        kind: IssueKindInput,
        severity: SeverityInput,
        role: Option<RoleInput>,
        template: Option<usize>,
        atom: Option<usize>,
        map: Option<i32>,
        maps: Vec<i32>,
        detail: String,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: ck::ReactionValidationIssue {
                kind: kind.into_rust()?,
                severity: severity.into_rust()?,
                role: role.map(RoleInput::into_rust).transpose()?,
                template,
                atom: atom.map(ck::AtomId::new),
                map,
                maps,
                detail,
            },
        })
    }
    #[getter]
    fn kind(&self) -> ReactionValidationIssueKind {
        self.inner.kind.into()
    }
    #[getter]
    fn severity(&self) -> ReactionValidationSeverity {
        self.inner.severity.into()
    }
    #[getter]
    fn role(&self) -> Option<ReactionRole> {
        self.inner.role.map(Into::into)
    }
    #[getter]
    fn template(&self) -> Option<usize> {
        self.inner.template
    }
    #[getter]
    fn atom(&self) -> Option<usize> {
        self.inner.atom.map(ck::AtomId::index)
    }
    #[getter]
    fn map(&self) -> Option<i32> {
        self.inner.map
    }
    #[getter]
    fn maps(&self) -> Vec<i32> {
        self.inner.maps.clone()
    }
    #[getter]
    fn detail(&self) -> &str {
        &self.inner.detail
    }
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct ReactionValidationReport {
    pub(crate) inner: ck::ReactionValidationReport,
}
#[pymethods]
impl ReactionValidationReport {
    #[new]
    fn new(warnings: Vec<ReactionValidationIssue>, errors: Vec<ReactionValidationIssue>) -> Self {
        Self {
            inner: ck::ReactionValidationReport {
                warnings: warnings.into_iter().map(|value| value.inner).collect(),
                errors: errors.into_iter().map(|value| value.inner).collect(),
            },
        }
    }
    #[getter]
    fn warnings(&self) -> Vec<ReactionValidationIssue> {
        self.inner
            .warnings
            .iter()
            .cloned()
            .map(|inner| ReactionValidationIssue { inner })
            .collect()
    }
    #[getter]
    fn errors(&self) -> Vec<ReactionValidationIssue> {
        self.inner
            .errors
            .iter()
            .cloned()
            .map(|inner| ReactionValidationIssue { inner })
            .collect()
    }
    #[getter]
    fn num_warnings(&self) -> usize {
        self.inner.num_warnings()
    }
    #[getter]
    fn num_errors(&self) -> usize {
        self.inner.num_errors()
    }
    #[getter]
    fn is_valid(&self) -> bool {
        self.inner.is_valid()
    }
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct ReactionApplyResult {
    pub(crate) inner: ck::ReactionApplyResult,
}
#[pymethods]
impl ReactionApplyResult {
    #[new]
    fn new(molecule: &Molecule, changed: bool) -> Self {
        Self {
            inner: ck::ReactionApplyResult::from((molecule.inner.clone(), changed)),
        }
    }
    #[getter]
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    #[getter]
    fn changed(&self) -> bool {
        self.inner.changed()
    }
}
