//! Public writer projection and its exact native fragment error closure.
use crate::canonical_errors::{CanonicalRankError, PropertyStringError};
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::{AtomPropertyError, BondValueError};
use crate::model_errors::StereoGroupError;
use crate::notation_params::{FragmentCxSmilesWriteParams, FragmentSmilesWriteParams};
use crate::ring_errors::RingFindingError;
use crate::smiles_parse_errors::SmilesParseError;
use crate::stereo_errors::LegacyStereoError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

#[derive(FromPyObject)]
pub enum FragmentSmilesInput {
    Atoms(Vec<usize>),
    Params(Py<FragmentSmilesWriteParams>),
}
#[derive(FromPyObject)]
pub enum FragmentCxSmilesInput {
    Atoms(Vec<usize>),
    Params(Py<FragmentCxSmilesWriteParams>),
}

pub(crate) fn text(
    py: Python<'_>,
    result: Result<ck::PropertyText, ck::SmilesWriteError>,
) -> PyResult<String> {
    result
        .map_err(SmilesWriteError::from_rust)
        .and_then(|value| crate::property_values::decode_source_text(py, &value))
}
pub(crate) fn texts(
    py: Python<'_>,
    result: Result<Vec<ck::PropertyText>, ck::SmilesWriteError>,
) -> PyResult<Vec<String>> {
    result
        .map_err(SmilesWriteError::from_rust)?
        .iter()
        .map(|value| crate::property_values::decode_source_text(py, value))
        .collect()
}

type FragmentSource = Source<
    Source<
        Py<StereoGroupError>,
        Source<Py<ValenceError>, Source<Py<RingFindingError>, Py<LegacyStereoError>>>,
    >,
    Source<
        Py<CanonicalRankError>,
        Source<
            Py<SmilesParseError>,
            Source<Py<PropertyStringError>, Source<Py<AtomPropertyError>, Py<BondValueError>>>,
        >,
    >,
>;
error_projection!(FragmentWriteInputError, FragmentWriteInputErrorKind, ck::FragmentWriteInputError, "smiles", {
    StereoGroup(_), NoAtomsProvided, RootAtomOutOfRange { atom_index, atom_count },
    RootAtomNotSelected { atom_index }, AtomSymbolsTooShort { len, expected },
    BondSymbolsTooShort { len, expected }, AtomOutOfRange { atom_index, atom_count },
    BondOutOfRange { bond_index: usize, bond_count: usize }, RootAtomRequiresSingleFragment { atom_index, fragment_count: usize },
    InvalidTopology(_), ValenceAssignmentLengthMismatch { atom_count, explicit_len: usize, implicit_len: usize },
    Valence(_), Ring(_), Stereo(_), Rank(_), Writer(_), CipProperty(_), AtomProperty(_), StereoBond(_),
}, cause = |inner: &ck::FragmentWriteInputError| match inner {
    ck::FragmentWriteInputError::StereoGroup(value) => Some(StereoGroupError::from_rust(*value)),
    ck::FragmentWriteInputError::Valence(value) => Some(ValenceError::from_rust(value.clone())),
    ck::FragmentWriteInputError::Ring(value) => Some(RingFindingError::from_rust(value.clone())),
    ck::FragmentWriteInputError::Stereo(value) => Some(LegacyStereoError::from_rust(value.clone())),
    ck::FragmentWriteInputError::Rank(value) => Some(CanonicalRankError::from_rust(value.clone())),
    ck::FragmentWriteInputError::Writer(value) => Some(SmilesParseError::from_rust(value.clone())),
    ck::FragmentWriteInputError::CipProperty(value) => Some(PropertyStringError::from_rust(*value)),
    ck::FragmentWriteInputError::AtomProperty(value) => Some(AtomPropertyError::from_rust(value.clone())),
    ck::FragmentWriteInputError::StereoBond(value) => Some(BondValueError::from_rust(value.clone())),
    ck::FragmentWriteInputError::NoAtomsProvided | ck::FragmentWriteInputError::RootAtomOutOfRange { .. }
    | ck::FragmentWriteInputError::RootAtomNotSelected { .. } | ck::FragmentWriteInputError::AtomSymbolsTooShort { .. }
    | ck::FragmentWriteInputError::BondSymbolsTooShort { .. } | ck::FragmentWriteInputError::AtomOutOfRange { .. }
    | ck::FragmentWriteInputError::BondOutOfRange { .. } | ck::FragmentWriteInputError::RootAtomRequiresSingleFragment { .. }
    | ck::FragmentWriteInputError::InvalidTopology(_) | ck::FragmentWriteInputError::ValenceAssignmentLengthMismatch { .. } => None,
}, source = Option<FragmentSource>, methods = {
    #[getter]
    fn atom_index(&self) -> Option<usize> { match self.inner {
        ck::FragmentWriteInputError::RootAtomOutOfRange { atom_index, .. }
        | ck::FragmentWriteInputError::RootAtomNotSelected { atom_index }
        | ck::FragmentWriteInputError::AtomOutOfRange { atom_index, .. }
        | ck::FragmentWriteInputError::RootAtomRequiresSingleFragment { atom_index, .. } => Some(atom_index), _ => None,
    } }
    #[getter]
    fn atom_count(&self) -> Option<usize> { match self.inner {
        ck::FragmentWriteInputError::RootAtomOutOfRange { atom_count, .. }
        | ck::FragmentWriteInputError::AtomOutOfRange { atom_count, .. }
        | ck::FragmentWriteInputError::ValenceAssignmentLengthMismatch { atom_count, .. } => Some(atom_count), _ => None,
    } }
    #[getter]
    fn len(&self) -> Option<usize> { match self.inner {
        ck::FragmentWriteInputError::AtomSymbolsTooShort { len, .. } | ck::FragmentWriteInputError::BondSymbolsTooShort { len, .. } => Some(len), _ => None,
    } }
    #[getter]
    fn expected(&self) -> Option<usize> { match self.inner {
        ck::FragmentWriteInputError::AtomSymbolsTooShort { expected, .. } | ck::FragmentWriteInputError::BondSymbolsTooShort { expected, .. } => Some(expected), _ => None,
    } }
    #[getter]
    fn detail(&self) -> Option<&str> { match &self.inner { ck::FragmentWriteInputError::InvalidTopology(value) => Some(value), _ => None } }
});

error_projection!(retained SmilesWriteError, SmilesWriteErrorKind, ck::SmilesWriteError, "smiles", {
    Write(_), Fragment(_),
}, cause = |inner: &crate::retained_error::RetainedError<ck::SmilesWriteError>| match inner.as_ref() {
    ck::SmilesWriteError::Write(value) => Some(SmilesParseError::from_rust(value.clone())),
    ck::SmilesWriteError::Fragment(value) => Some(FragmentWriteInputError::from_rust(value.clone())),
}, source = Option<Source<Py<SmilesParseError>, Py<FragmentWriteInputError>>>);

pub(crate) fn keywords(
    isomeric_smiles: Option<bool>,
    kekule: Option<bool>,
    canonical: Option<bool>,
    clean_stereo: Option<bool>,
    rooted_at_atom: Option<usize>,
    all_bonds_explicit: Option<bool>,
    all_hydrogens_explicit: Option<bool>,
    include_dative_bonds: Option<bool>,
    ignore_atom_map_numbers: Option<bool>,
) -> crate::chemical_params::SmilesWriteParams {
    let defaults = ck::SmilesWriteParams::default();
    crate::chemical_params::SmilesWriteParams {
        inner: ck::SmilesWriteParams {
            isomeric_smiles: isomeric_smiles.unwrap_or(defaults.isomeric_smiles),
            kekule: kekule.unwrap_or(defaults.kekule),
            canonical: canonical.unwrap_or(defaults.canonical),
            clean_stereo: clean_stereo.unwrap_or(defaults.clean_stereo),
            rooted_at_atom: rooted_at_atom.map(ck::AtomId::new),
            all_bonds_explicit: all_bonds_explicit.unwrap_or(defaults.all_bonds_explicit),
            all_hydrogens_explicit: all_hydrogens_explicit
                .unwrap_or(defaults.all_hydrogens_explicit),
            include_dative_bonds: include_dative_bonds.unwrap_or(defaults.include_dative_bonds),
            ignore_atom_map_numbers: ignore_atom_map_numbers
                .unwrap_or(defaults.ignore_atom_map_numbers),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemical_params::SmilesWriteParams;
    use crate::native::Molecule;
    use crate::notation_params::{CxSmilesFields, CxSmilesWriteParams};
    use crate::params::RandomSmilesWriteParams;
    use pyo3::types::PyDict;

    fn decoded(value: ck::PropertyText) -> String {
        String::from_utf8(value.into_bytes()).unwrap()
    }

    #[test]
    fn complete_writer_family_forwards_native_options_selections_random_order_and_errors() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            for input in ["CC(C)O", "C[C@H](F)Cl", "c1ccncc1"] {
                let molecule = Py::new(
                    py,
                    Molecule {
                        inner: ck::Molecule::from_smiles(input).unwrap(),
                    },
                )?;
                let before = molecule.borrow(py).inner.to_binary().unwrap();
                let count = molecule.borrow(py).inner.num_atoms();
                let default = molecule.borrow(py).inner.to_smiles().unwrap();
                assert_eq!(
                    molecule
                        .bind(py)
                        .call_method0("to_smiles")?
                        .extract::<String>()?,
                    decoded(default)
                );
                let native_params = ck::SmilesWriteParams {
                    isomeric_smiles: false,
                    kekule: true,
                    canonical: false,
                    clean_stereo: false,
                    rooted_at_atom: Some(ck::AtomId::new(0)),
                    all_bonds_explicit: true,
                    all_hydrogens_explicit: true,
                    include_dative_bonds: false,
                    ignore_atom_map_numbers: true,
                };
                let expected = decoded(
                    molecule
                        .borrow(py)
                        .inner
                        .to_smiles_with_params(&native_params)
                        .unwrap(),
                );
                let params = Py::new(
                    py,
                    SmilesWriteParams {
                        inner: native_params,
                    },
                )?;
                let kwargs = PyDict::new(py);
                kwargs.set_item("isomeric_smiles", native_params.isomeric_smiles)?;
                kwargs.set_item("kekule", native_params.kekule)?;
                kwargs.set_item("canonical", native_params.canonical)?;
                kwargs.set_item("clean_stereo", native_params.clean_stereo)?;
                kwargs.set_item("rooted_at_atom", 0)?;
                kwargs.set_item("all_bonds_explicit", native_params.all_bonds_explicit)?;
                kwargs.set_item(
                    "all_hydrogens_explicit",
                    native_params.all_hydrogens_explicit,
                )?;
                kwargs.set_item("include_dative_bonds", native_params.include_dative_bonds)?;
                kwargs.set_item(
                    "ignore_atom_map_numbers",
                    native_params.ignore_atom_map_numbers,
                )?;
                for actual in [
                    molecule
                        .bind(py)
                        .call_method1("to_smiles", (params.clone_ref(py),))?,
                    molecule
                        .bind(py)
                        .call_method1("to_smiles_with_params", (params.clone_ref(py),))?,
                    molecule
                        .bind(py)
                        .call_method("to_smiles", (), Some(&kwargs))?,
                ] {
                    assert_eq!(actual.extract::<String>()?, expected);
                }
                assert!(
                    molecule
                        .bind(py)
                        .call_method("to_smiles", (params.clone_ref(py),), Some(&kwargs))
                        .unwrap_err()
                        .is_instance_of::<pyo3::exceptions::PyTypeError>(py)
                );

                let cx_default = decoded(molecule.borrow(py).inner.to_cx_smiles().unwrap());
                assert_eq!(
                    molecule
                        .bind(py)
                        .call_method0("to_cx_smiles")?
                        .extract::<String>()?,
                    cx_default
                );
                let native_cx = ck::CxSmilesWriteParams {
                    smiles: native_params,
                    fields: ck::CxSmilesFields::NONE,
                    coordinate_selection: ck::CxCoordinateSelection::Auto,
                };
                let expected = decoded(
                    molecule
                        .borrow(py)
                        .inner
                        .to_cx_smiles_with_params(&native_cx)
                        .unwrap(),
                );
                let cx = Py::new(
                    py,
                    CxSmilesWriteParams {
                        inner: native_cx.clone(),
                    },
                )?;
                let cx_kwargs = PyDict::new(py);
                cx_kwargs.set_item("smiles", params.clone_ref(py))?;
                cx_kwargs.set_item(
                    "fields",
                    Py::new(
                        py,
                        CxSmilesFields {
                            inner: ck::CxSmilesFields::NONE,
                        },
                    )?,
                )?;
                for actual in [
                    molecule
                        .bind(py)
                        .call_method1("to_cx_smiles", (cx.clone_ref(py),))?,
                    molecule
                        .bind(py)
                        .call_method1("to_cx_smiles_with_params", (cx.clone_ref(py),))?,
                    molecule
                        .bind(py)
                        .call_method("to_cx_smiles", (), Some(&cx_kwargs))?,
                ] {
                    assert_eq!(actual.extract::<String>()?, expected);
                }
                assert!(
                    molecule
                        .bind(py)
                        .call_method("to_cx_smiles", (cx.clone_ref(py),), Some(&cx_kwargs))
                        .unwrap_err()
                        .is_instance_of::<pyo3::exceptions::PyTypeError>(py)
                );

                let atoms = vec![ck::AtomId::new(0), ck::AtomId::new(1)];
                let atom_indices = vec![0_usize, 1];
                let expected = decoded(
                    molecule
                        .borrow(py)
                        .inner
                        .to_fragment_smiles(&atoms)
                        .unwrap(),
                );
                assert_eq!(
                    molecule
                        .bind(py)
                        .call_method1("to_fragment_smiles", (atom_indices.clone(),))?
                        .extract::<String>()?,
                    expected
                );
                let expected = decoded(
                    molecule
                        .borrow(py)
                        .inner
                        .to_fragment_cx_smiles(&atoms)
                        .unwrap(),
                );
                assert_eq!(
                    molecule
                        .bind(py)
                        .call_method1("to_fragment_cx_smiles", (atom_indices.clone(),))?
                        .extract::<String>()?,
                    expected
                );
                // A present empty bond selection stays empty; no guessed default.
                let fragment = ck::FragmentSmilesWriteParams {
                    atoms: atoms.clone(),
                    smiles: native_params,
                    bonds: Some(vec![]),
                    atom_symbols: Some(vec!["C".into(); count]),
                    bond_symbols: None,
                };
                let expected = decoded(
                    molecule
                        .borrow(py)
                        .inner
                        .to_fragment_smiles_with_params(&fragment)
                        .unwrap(),
                );
                let fragment = Py::new(py, FragmentSmilesWriteParams { inner: fragment })?;
                let fragment_kwargs = PyDict::new(py);
                fragment_kwargs.set_item("smiles", params.clone_ref(py))?;
                fragment_kwargs.set_item("bonds", Vec::<usize>::new())?;
                fragment_kwargs.set_item("atom_symbols", vec!["C"; count])?;
                let named_params = PyDict::new(py);
                named_params.set_item("params", fragment.clone_ref(py))?;
                for actual in [
                    molecule
                        .bind(py)
                        .call_method1("to_fragment_smiles", (fragment.clone_ref(py),))?,
                    molecule
                        .bind(py)
                        .call_method("to_fragment_smiles", (), Some(&named_params))?,
                    molecule.bind(py).call_method1(
                        "to_fragment_smiles_with_params",
                        (fragment.clone_ref(py),),
                    )?,
                    molecule.bind(py).call_method(
                        "to_fragment_smiles",
                        (atom_indices.clone(),),
                        Some(&fragment_kwargs),
                    )?,
                ] {
                    assert_eq!(actual.extract::<String>()?, expected);
                }
                assert!(
                    molecule
                        .bind(py)
                        .call_method(
                            "to_fragment_smiles",
                            (fragment.clone_ref(py),),
                            Some(&fragment_kwargs)
                        )
                        .unwrap_err()
                        .is_instance_of::<pyo3::exceptions::PyTypeError>(py)
                );

                let fragment_cx = ck::FragmentCxSmilesWriteParams {
                    atoms: atoms.clone(),
                    cx: native_cx,
                    bonds: Some(vec![]),
                    atom_symbols: Some(vec!["C".into(); count]),
                    bond_symbols: None,
                };
                let expected = decoded(
                    molecule
                        .borrow(py)
                        .inner
                        .to_fragment_cx_smiles_with_params(&fragment_cx)
                        .unwrap(),
                );
                let fragment_cx = Py::new(py, FragmentCxSmilesWriteParams { inner: fragment_cx })?;
                let fragment_cx_kwargs = PyDict::new(py);
                fragment_cx_kwargs.set_item("cx", cx.clone_ref(py))?;
                fragment_cx_kwargs.set_item("bonds", Vec::<usize>::new())?;
                fragment_cx_kwargs.set_item("atom_symbols", vec!["C"; count])?;
                let named_params = PyDict::new(py);
                named_params.set_item("params", fragment_cx.clone_ref(py))?;
                for actual in [
                    molecule
                        .bind(py)
                        .call_method1("to_fragment_cx_smiles", (fragment_cx.clone_ref(py),))?,
                    molecule.bind(py).call_method(
                        "to_fragment_cx_smiles",
                        (),
                        Some(&named_params),
                    )?,
                    molecule.bind(py).call_method1(
                        "to_fragment_cx_smiles_with_params",
                        (fragment_cx.clone_ref(py),),
                    )?,
                    molecule.bind(py).call_method(
                        "to_fragment_cx_smiles",
                        (atom_indices.clone(),),
                        Some(&fragment_cx_kwargs),
                    )?,
                ] {
                    assert_eq!(actual.extract::<String>()?, expected);
                }
                assert!(
                    molecule
                        .bind(py)
                        .call_method(
                            "to_fragment_cx_smiles",
                            (fragment_cx.clone_ref(py),),
                            Some(&fragment_cx_kwargs)
                        )
                        .unwrap_err()
                        .is_instance_of::<pyo3::exceptions::PyTypeError>(py)
                );

                let expected = molecule
                    .borrow(py)
                    .inner
                    .to_random_smiles(16, 12345)
                    .unwrap()
                    .into_iter()
                    .map(decoded)
                    .collect::<Vec<_>>();
                assert_eq!(
                    molecule
                        .bind(py)
                        .call_method1("to_random_smiles", (16, 12345))?
                        .extract::<Vec<String>>()?,
                    expected
                );
                let random = ck::RandomSmilesWriteParams {
                    isomeric_smiles: false,
                    kekule: true,
                    all_bonds_explicit: true,
                    all_hydrogens_explicit: true,
                };
                let expected = molecule
                    .borrow(py)
                    .inner
                    .to_random_smiles_with_params(16, 12345, &random)
                    .unwrap()
                    .into_iter()
                    .map(decoded)
                    .collect::<Vec<_>>();
                let random = Py::new(py, RandomSmilesWriteParams { inner: random })?;
                let random_kwargs = PyDict::new(py);
                random_kwargs.set_item("isomeric_smiles", false)?;
                random_kwargs.set_item("kekule", true)?;
                random_kwargs.set_item("all_bonds_explicit", true)?;
                random_kwargs.set_item("all_hydrogens_explicit", true)?;
                for actual in [
                    molecule
                        .bind(py)
                        .call_method1("to_random_smiles", (16, 12345, random.clone_ref(py)))?,
                    molecule.bind(py).call_method1(
                        "to_random_smiles_with_params",
                        (16, 12345, random.clone_ref(py)),
                    )?,
                    molecule.bind(py).call_method(
                        "to_random_smiles",
                        (16, 12345),
                        Some(&random_kwargs),
                    )?,
                ] {
                    assert_eq!(actual.extract::<Vec<String>>()?, expected);
                }
                assert!(
                    molecule
                        .bind(py)
                        .call_method(
                            "to_random_smiles",
                            (16, 12345, random.clone_ref(py)),
                            Some(&random_kwargs)
                        )
                        .unwrap_err()
                        .is_instance_of::<pyo3::exceptions::PyTypeError>(py)
                );
                assert!(
                    molecule
                        .bind(py)
                        .call_method1("to_random_smiles", (0, u32::MAX))?
                        .extract::<Vec<String>>()?
                        .is_empty()
                );

                let root = PyDict::new(py);
                root.set_item("rooted_at_atom", usize::MAX)?;
                let error = molecule
                    .bind(py)
                    .call_method("to_smiles", (), Some(&root))
                    .unwrap_err();
                assert!(error.is_instance_of::<SmilesWriteError>(py));
                let source = error.value(py).getattr("source")?;
                assert!(source.is_instance_of::<SmilesParseError>());
                assert!(source.is(&error.value(py).getattr("__cause__")?));
                assert_eq!(
                    source.getattr("atom_index")?.extract::<usize>()?,
                    usize::MAX
                );
                let error = molecule
                    .bind(py)
                    .call_method1("to_fragment_smiles", (vec![usize::MAX],))
                    .unwrap_err();
                assert!(error.is_instance_of::<SmilesWriteError>(py));
                let source = error.value(py).getattr("source")?;
                assert!(source.is_instance_of::<FragmentWriteInputError>());
                assert_eq!(
                    source.getattr("atom_index")?.extract::<usize>()?,
                    usize::MAX
                );
                assert!(source.getattr("source")?.is_none());
                assert!(source.setattr("atom_index", 0).is_err());
                let empty = molecule
                    .bind(py)
                    .call_method1("to_fragment_smiles", (Vec::<usize>::new(),))
                    .unwrap_err();
                assert!(
                    empty
                        .value(py)
                        .getattr("source")?
                        .getattr("variant")?
                        .extract::<FragmentWriteInputErrorKind>()?
                        == FragmentWriteInputErrorKind::NoAtomsProvided
                );
                assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);
            }
            let bytes = ck::PropertyText::from_bytes(&[b'x', 0xff]);
            let error = text(py, Ok(bytes)).unwrap_err();
            assert!(error.is_instance_of::<pyo3::exceptions::PyUnicodeDecodeError>(py));
            Ok(())
        })
        .unwrap();
    }
}
