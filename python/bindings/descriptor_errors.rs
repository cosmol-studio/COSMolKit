//! Native descriptor read boundaries and algorithm/search error values.
use crate::error_framework::{Source, error_projection};
use crate::hydrogen_errors::HydrogenError;
use crate::match_errors::{QueryMatchContextError, SubstructMatchError};
use crate::model_errors::TopologyValidationError;
use crate::path_errors::PathError;
use crate::ring_errors::RingFindingError;
use crate::smarts_errors::SmartsParseError;
use crate::stereo_errors::LegacyStereoError;
use crate::topology_edit_errors::TopologyEditError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(DescriptorSearchCause, DescriptorSearchCauseKind, ck::DescriptorSearchCause, "descriptors", {
    Compile(_), Match(_), Context(_),
}, cause = |inner: &ck::DescriptorSearchCause| match inner {
    ck::DescriptorSearchCause::Compile(value) => Some(SmartsParseError::from_rust(value.clone())),
    ck::DescriptorSearchCause::Match(value) => Some(SubstructMatchError::from_rust(value.clone())),
    ck::DescriptorSearchCause::Context(value) => Some(QueryMatchContextError::from_rust(value.clone())),
}, source = Option<Source<crate::error_framework::ErrorObject<SmartsParseError>, Source<crate::error_framework::ErrorObject<SubstructMatchError>, crate::error_framework::ErrorObject<QueryMatchContextError>>>>);

pub(crate) type DescriptorSource = Source<
    crate::error_framework::ErrorObject<TopologyEditError>,
    Source<
        crate::error_framework::ErrorObject<TopologyValidationError>,
        Source<
            crate::error_framework::ErrorObject<PathError>,
            Source<
                crate::error_framework::ErrorObject<HydrogenError>,
                Source<
                    crate::error_framework::ErrorObject<ValenceError>,
                    Source<
                        crate::error_framework::ErrorObject<RingFindingError>,
                        Source<
                            crate::error_framework::ErrorObject<DescriptorSearchCause>,
                            crate::error_framework::ErrorObject<LegacyStereoError>,
                        >,
                    >,
                >,
            >,
        >,
    >,
>;

error_projection!(DescriptorError, DescriptorErrorKind, ck::DescriptorError, "descriptors", {
    TopologyEdit { function, source },
    MissingFinalHydrogenState { field },
    InvalidConnectivityPath { function, expected_rows: usize, actual_rows: Option<usize> },
    InvalidTopology { function, source },
    Path { function, source },
    InvalidHallKierContributionRows { actual, minimum: usize },
    MissingComputedScalar { function, include_sulfur_phosphorus: bool },
    MissingLabuteHydrogens { function },
    MissingLabuteAsa { function },
    MismatchedBinArrays { contribs_len: usize, bin_prop_len: usize, bins_len: usize },
    CrippenParamNumeric { field, cell: String },
    MissingCrippenMrContributions { function },
    MissingCrippenMr { function },
    InvalidCrippenOptionalRows { field, actual, expected },
    MissingCrippenDefaultPattern { row: u32 },
    Hydrogens { function, source },
    CountOverflow { function, field },
    Valence { function, source },
    InvalidValenceRows { function, field, actual, expected },
    Unsupported { function, detail: String },
    Ring { function, source },
    Search { function, source },
    Stereo { function, source },
}, cause = |inner: &ck::DescriptorError| match inner {
    ck::DescriptorError::TopologyEdit { source, .. } => Some(TopologyEditError::from_rust(source.clone())),
    ck::DescriptorError::InvalidTopology { source, .. } => Some(TopologyValidationError::from_rust(source.clone())),
    ck::DescriptorError::Path { source, .. } => Some(PathError::from_rust(source.clone())),
    ck::DescriptorError::Hydrogens { source, .. } => Some(HydrogenError::from_rust(source.clone())),
    ck::DescriptorError::Valence { source, .. } => Some(ValenceError::from_rust(source.clone())),
    ck::DescriptorError::Ring { source, .. } => Some(RingFindingError::from_rust(source.clone())),
    ck::DescriptorError::Search { source, .. } => Some(DescriptorSearchCause::from_rust(source.clone())),
    ck::DescriptorError::Stereo { source, .. } => Some(LegacyStereoError::from_rust(source.clone())),
    ck::DescriptorError::MissingFinalHydrogenState { .. }
    | ck::DescriptorError::InvalidConnectivityPath { .. }
    | ck::DescriptorError::InvalidHallKierContributionRows { .. }
    | ck::DescriptorError::MissingComputedScalar { .. }
    | ck::DescriptorError::MissingLabuteHydrogens { .. }
    | ck::DescriptorError::MissingLabuteAsa { .. }
    | ck::DescriptorError::MismatchedBinArrays { .. }
    | ck::DescriptorError::CrippenParamNumeric { .. }
    | ck::DescriptorError::MissingCrippenMrContributions { .. }
    | ck::DescriptorError::MissingCrippenMr { .. }
    | ck::DescriptorError::InvalidCrippenOptionalRows { .. }
    | ck::DescriptorError::MissingCrippenDefaultPattern { .. }
    | ck::DescriptorError::CountOverflow { .. }
    | ck::DescriptorError::InvalidValenceRows { .. }
    | ck::DescriptorError::Unsupported { .. } => None,
}, source = Option<DescriptorSource>, methods = {
    #[getter]
    #[doc = "Descriptor function reporting the failure. Present for `DescriptorError.TopologyEdit`, `DescriptorError.InvalidConnectivityPath`, `DescriptorError.InvalidTopology`, `DescriptorError.Path`, `DescriptorError.MissingComputedScalar`, `DescriptorError.MissingLabuteHydrogens`, `DescriptorError.MissingLabuteAsa`, `DescriptorError.MissingCrippenMrContributions`, `DescriptorError.MissingCrippenMr`, `DescriptorError.Hydrogens`, `DescriptorError.CountOverflow`, `DescriptorError.Valence`, `DescriptorError.InvalidValenceRows`, `DescriptorError.Unsupported`, `DescriptorError.Ring`, `DescriptorError.Search`, `DescriptorError.Stereo`; None when unavailable."] fn function(&self) -> Option<&'static str> { match self.inner {
        ck::DescriptorError::TopologyEdit { function, .. }
        | ck::DescriptorError::InvalidConnectivityPath { function, .. }
        | ck::DescriptorError::InvalidTopology { function, .. }
        | ck::DescriptorError::Path { function, .. }
        | ck::DescriptorError::MissingComputedScalar { function, .. }
        | ck::DescriptorError::MissingLabuteHydrogens { function, .. }
        | ck::DescriptorError::MissingLabuteAsa { function, .. }
        | ck::DescriptorError::MissingCrippenMrContributions { function, .. }
        | ck::DescriptorError::MissingCrippenMr { function, .. }
        | ck::DescriptorError::Hydrogens { function, .. }
        | ck::DescriptorError::CountOverflow { function, .. }
        | ck::DescriptorError::Valence { function, .. }
        | ck::DescriptorError::InvalidValenceRows { function, .. }
        | ck::DescriptorError::Unsupported { function, .. }
        | ck::DescriptorError::Ring { function, .. }
        | ck::DescriptorError::Search { function, .. }
        | ck::DescriptorError::Stereo { function, .. } => Some(function), _ => None,
    } }
    #[getter]
    #[doc = "Input field associated with the failure. Present for `DescriptorError.MissingFinalHydrogenState`, `DescriptorError.CrippenParamNumeric`, `DescriptorError.InvalidCrippenOptionalRows`, `DescriptorError.CountOverflow`, `DescriptorError.InvalidValenceRows`; None when unavailable."] fn field(&self) -> Option<&'static str> { match self.inner {
        ck::DescriptorError::MissingFinalHydrogenState { field, .. }
        | ck::DescriptorError::CrippenParamNumeric { field, .. }
        | ck::DescriptorError::InvalidCrippenOptionalRows { field, .. }
        | ck::DescriptorError::CountOverflow { field, .. }
        | ck::DescriptorError::InvalidValenceRows { field, .. } => Some(field), _ => None,
    } }
    #[getter]
    #[doc = "Observed count or value. Present for `DescriptorError.InvalidHallKierContributionRows`, `DescriptorError.InvalidCrippenOptionalRows`, `DescriptorError.InvalidValenceRows`; None when unavailable."] fn actual(&self) -> Option<usize> { match self.inner {
        ck::DescriptorError::InvalidHallKierContributionRows { actual, .. }
        | ck::DescriptorError::InvalidCrippenOptionalRows { actual, .. }
        | ck::DescriptorError::InvalidValenceRows { actual, .. } => Some(actual), _ => None,
    } }
    #[getter]
    #[doc = "Expected count or value. Present for `DescriptorError.InvalidCrippenOptionalRows`, `DescriptorError.InvalidValenceRows`; None when unavailable."] fn expected(&self) -> Option<usize> { match self.inner {
        ck::DescriptorError::InvalidCrippenOptionalRows { expected, .. }
        | ck::DescriptorError::InvalidValenceRows { expected, .. } => Some(expected), _ => None,
    } }
});

error_projection!(DescriptorReadError, DescriptorReadErrorKind, ck::DescriptorReadError, "descriptors", {
    MissingPreparedValence, CachePoisoned, MissingInitializedRings, Algorithm { source },
}, cause = |inner: &ck::DescriptorReadError| match inner {
    ck::DescriptorReadError::Algorithm { source } => Some(DescriptorError::from_rust(source.clone())),
    ck::DescriptorReadError::MissingPreparedValence | ck::DescriptorReadError::CachePoisoned
    | ck::DescriptorReadError::MissingInitializedRings => None,
}, source = Option<crate::error_framework::ErrorObject<DescriptorError>>);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn descriptor_producers_forward_all_native_options_bits_rows_and_errors() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            type Count = fn(&ck::Molecule) -> Result<u32, ck::DescriptorReadError>;
            type Scalar = fn(&ck::Molecule) -> Result<f64, ck::DescriptorReadError>;
            let counts: [(&str, Count); 13] = [
                ("num_heavy_atoms", ck::Molecule::num_heavy_atoms),
                ("total_atom_count", ck::Molecule::total_atom_count),
                ("lipinski_hba", ck::Molecule::lipinski_hba),
                ("lipinski_hbd", ck::Molecule::lipinski_hbd),
                ("num_heteroatoms", ck::Molecule::num_heteroatoms),
                ("num_hba", ck::Molecule::num_hba),
                ("num_hbd", ck::Molecule::num_hbd),
                ("num_amide_bonds", ck::Molecule::num_amide_bonds),
                ("num_atom_stereo_centers", ck::Molecule::num_atom_stereo_centers),
                ("num_unspecified_atom_stereo_centers", ck::Molecule::num_unspecified_atom_stereo_centers),
                ("num_spiro_atoms", ck::Molecule::num_spiro_atoms),
                ("num_bridgehead_atoms", ck::Molecule::num_bridgehead_atoms),
                ("num_rotatable_bonds", ck::Molecule::num_rotatable_bonds),
            ];
            let scalars: [(&str, Scalar); 26] = [
                ("fraction_csp3", ck::Molecule::fraction_csp3),
                ("qed", ck::Molecule::qed),
                ("slogp_vsa_1", ck::Molecule::slogp_vsa_1),
                ("slogp_vsa_2", ck::Molecule::slogp_vsa_2),
                ("slogp_vsa_3", ck::Molecule::slogp_vsa_3),
                ("slogp_vsa_4", ck::Molecule::slogp_vsa_4),
                ("slogp_vsa_5", ck::Molecule::slogp_vsa_5),
                ("slogp_vsa_6", ck::Molecule::slogp_vsa_6),
                ("slogp_vsa_7", ck::Molecule::slogp_vsa_7),
                ("slogp_vsa_8", ck::Molecule::slogp_vsa_8),
                ("slogp_vsa_9", ck::Molecule::slogp_vsa_9),
                ("slogp_vsa_10", ck::Molecule::slogp_vsa_10),
                ("slogp_vsa_11", ck::Molecule::slogp_vsa_11),
                ("slogp_vsa_12", ck::Molecule::slogp_vsa_12),
                ("smr_vsa_1", ck::Molecule::smr_vsa_1),
                ("smr_vsa_2", ck::Molecule::smr_vsa_2),
                ("smr_vsa_3", ck::Molecule::smr_vsa_3),
                ("smr_vsa_4", ck::Molecule::smr_vsa_4),
                ("smr_vsa_5", ck::Molecule::smr_vsa_5),
                ("smr_vsa_6", ck::Molecule::smr_vsa_6),
                ("smr_vsa_7", ck::Molecule::smr_vsa_7),
                ("smr_vsa_8", ck::Molecule::smr_vsa_8),
                ("smr_vsa_9", ck::Molecule::smr_vsa_9),
                ("smr_vsa_10", ck::Molecule::smr_vsa_10),
                ("molecular_weight", ck::Molecule::molecular_weight),
                ("exact_molecular_weight", ck::Molecule::exact_molecular_weight),
            ];
            let bits = |rows: Vec<f64>| rows.into_iter().map(f64::to_bits).collect::<Vec<_>>();
            for text in ["CC(=O)NCCO", "[2H]C([13CH3])O", "C1CCC2(CC1)CCCC2"] {
                for prepared in [true, false] {
                    let native = ck::Molecule::from_smiles_with_params(text, &ck::SmilesParseParams {
                        sanitize: prepared, remove_hs: prepared, ..Default::default()
                    }).unwrap();
                    let before = native.to_binary().unwrap();
                    let molecule = Py::new(py, crate::native::Molecule { inner: native })?;
                    for (name, run) in counts {
                        let expected = run(&molecule.borrow(py).inner);
                        let actual = molecule.bind(py).call_method0(name);
                        match expected {
                            Ok(value) => assert_eq!(actual?.extract::<u32>()?, value, "{name}"),
                            Err(expected) => {
                                let actual = actual.unwrap_err();
                                assert!(actual.is_instance(py, <DescriptorReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()), "{name}");
                                assert_eq!(actual.value(py).str()?.to_str()?, expected.to_string(), "{name}");
                            }
                        }
                    }
                    for (name, run) in scalars {
                        let expected = run(&molecule.borrow(py).inner);
                        let actual = molecule.bind(py).call_method0(name);
                        match expected {
                            Ok(value) => assert_eq!(actual?.extract::<f64>()?.to_bits(), value.to_bits(), "{name}"),
                            Err(expected) => {
                                let actual = actual.unwrap_err();
                                assert!(actual.is_instance(py, <DescriptorReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()), "{name}");
                                assert_eq!(actual.value(py).str()?.to_str()?, expected.to_string(), "{name}");
                            }
                        }
                    }
                    let formula = molecule.borrow(py).inner.molecular_formula();
                    match formula {
                        Ok(expected) => assert_eq!(molecule.bind(py).call_method0("molecular_formula")?.extract::<String>()?, expected),
                        Err(expected) => {
                            let actual = molecule.bind(py).call_method0("molecular_formula").unwrap_err();
                            assert!(actual.is_instance(py, <DescriptorReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
                            assert_eq!(actual.value(py).str()?.to_str()?, expected.to_string());
                        }
                    }
                    if prepared {
                        for option in [ck::RotatableBondsOptions::Default, ck::RotatableBondsOptions::NonStrict, ck::RotatableBondsOptions::Strict, ck::RotatableBondsOptions::StrictLinkages] {
                            let expected = molecule.borrow(py).inner.num_rotatable_bonds_with_params(&option).unwrap();
                            let option = crate::public_selectors::RotatableBondsOptions::from(option);
                            assert_eq!(molecule.bind(py).call_method1("num_rotatable_bonds", (option,))?.extract::<u32>()?, expected);
                            assert_eq!(molecule.bind(py).call_method1("num_rotatable_bonds_with_params", (option,))?.extract::<u32>()?, expected);
                        }
                        for (name, run) in [
                            ("molecular_weight", ck::Molecule::molecular_weight_with_params as fn(&ck::Molecule, bool) -> Result<f64, ck::DescriptorReadError>),
                            ("exact_molecular_weight", ck::Molecule::exact_molecular_weight_with_params),
                        ] {
                            for heavy in [false, true] {
                                let expected = run(&molecule.borrow(py).inner, heavy).unwrap();
                                let kwargs = pyo3::types::PyDict::new(py);
                                kwargs.set_item("only_heavy", heavy)?;
                                assert_eq!(molecule.bind(py).call_method(name, (), Some(&kwargs))?.extract::<f64>()?.to_bits(), expected.to_bits());
                                assert_eq!(molecule.bind(py).call_method1(format!("{name}_with_params"), (heavy,))?.extract::<f64>()?.to_bits(), expected.to_bits());
                            }
                        }
                        for option in [false, true] {
                            for force in [false, true] {
                                let expected = molecule.borrow(py).inner.molecular_formula_with_params(option, force).unwrap();
                                let kwargs = pyo3::types::PyDict::new(py);
                                kwargs.set_item("separate_isotopes", option)?;
                                kwargs.set_item("abbreviate_h_isotopes", force)?;
                                assert_eq!(molecule.bind(py).call_method("molecular_formula", (), Some(&kwargs))?.extract::<String>()?, expected);
                                assert_eq!(molecule.bind(py).call_method1("molecular_formula_with_params", (option, force))?.extract::<String>()?, expected);
                                for (name, key, run) in [
                                    ("labute_asa", "include_hydrogens", ck::Molecule::labute_asa_with_params as fn(&ck::Molecule, bool, bool) -> Result<f64, ck::DescriptorReadError>),
                                    ("tpsa", "include_sulfur_phosphorus", ck::Molecule::tpsa_with_params),
                                ] {
                                    let expected = run(&molecule.borrow(py).inner, option, force).unwrap();
                                    let kwargs = pyo3::types::PyDict::new(py);
                                    kwargs.set_item(key, option)?; kwargs.set_item("force", force)?;
                                    assert_eq!(molecule.bind(py).call_method(name, (), Some(&kwargs))?.extract::<f64>()?.to_bits(), expected.to_bits());
                                    assert_eq!(molecule.bind(py).call_method1(format!("{name}_with_params"), (option, force))?.extract::<f64>()?.to_bits(), expected.to_bits());
                                }
                                let kwargs = pyo3::types::PyDict::new(py);
                                kwargs.set_item("include_hydrogens", option)?; kwargs.set_item("force", force)?;
                                let expected = molecule.borrow(py).inner.crippen_descriptors_with_params(option, force).unwrap();
                                for actual in [
                                    molecule.bind(py).call_method("crippen_descriptors", (), Some(&kwargs))?,
                                    molecule.bind(py).call_method1("crippen_descriptors_with_params", (option, force))?,
                                ] {
                                    let actual = actual.extract::<PyRef<'_, crate::numerical_results::CrippenTotals>>()?;
                                    assert_eq!(actual.inner.logp.to_bits(), expected.logp.to_bits());
                                    assert_eq!(actual.inner.molar_refractivity.to_bits(), expected.molar_refractivity.to_bits());
                                }
                                let expected = molecule.borrow(py).inner.labute_asa_contributions_with_params(option, force).unwrap();
                                for actual in [
                                    molecule.bind(py).call_method("labute_asa_contributions", (), Some(&kwargs))?,
                                    molecule.bind(py).call_method1("labute_asa_contributions_with_params", (option, force))?,
                                ] {
                                    let actual = actual.extract::<PyRef<'_, crate::numerical_results::LabuteAsaContributions>>()?;
                                    assert_eq!(actual.inner.asa.to_bits(), expected.asa.to_bits());
                                    assert_eq!(actual.inner.hydrogen_contribution.to_bits(), expected.hydrogen_contribution.to_bits());
                                    assert_eq!(bits(actual.inner.atom_contributions.clone()), bits(expected.atom_contributions.clone()));
                                }
                            }
                        }
                        for (name, run) in [
                            ("slogp_vsa", ck::Molecule::slogp_vsa_with_params as fn(&ck::Molecule, Option<&[f64]>, bool) -> Result<Vec<f64>, ck::DescriptorReadError>),
                            ("smr_vsa", ck::Molecule::smr_vsa_with_params),
                        ] {
                            assert_eq!(bits(molecule.bind(py).call_method0(name)?.extract::<Vec<f64>>()?), bits(run(&molecule.borrow(py).inner, None, false).unwrap()));
                            for bins in [None, Some(vec![]), Some(vec![-1.0, 0.0, 0.0, 2.0])] {
                                for force in [false, true] {
                                    let expected = bits(run(&molecule.borrow(py).inner, bins.as_deref(), force).unwrap());
                                    let kwargs = pyo3::types::PyDict::new(py);
                                    kwargs.set_item("bins", bins.clone())?; kwargs.set_item("force", force)?;
                                    assert_eq!(bits(molecule.bind(py).call_method(name, (), Some(&kwargs))?.extract::<Vec<f64>>()?), expected);
                                    assert_eq!(bits(molecule.bind(py).call_method1(format!("{name}_with_params"), (bins.clone(), force))?.extract::<Vec<f64>>()?), expected);
                                }
                            }
                        }
                        // Default record/scalar paths remain native default calls, not independently recomputed values.
                        let expected = molecule.borrow(py).inner.crippen_descriptors().unwrap();
                        let actual = molecule.bind(py).call_method0("crippen_descriptors")?.extract::<PyRef<'_, crate::numerical_results::CrippenTotals>>()?;
                        assert_eq!(actual.inner.logp.to_bits(), expected.logp.to_bits());
                        assert_eq!(actual.inner.molar_refractivity.to_bits(), expected.molar_refractivity.to_bits());
                        let expected = molecule.borrow(py).inner.labute_asa_contributions().unwrap();
                        let actual = molecule.bind(py).call_method0("labute_asa_contributions")?.extract::<PyRef<'_, crate::numerical_results::LabuteAsaContributions>>()?;
                        assert_eq!(bits(actual.inner.atom_contributions.clone()), bits(expected.atom_contributions));
                        assert_eq!(molecule.bind(py).call_method0("labute_asa")?.extract::<f64>()?.to_bits(), molecule.borrow(py).inner.labute_asa().unwrap().to_bits());
                        assert_eq!(molecule.bind(py).call_method0("tpsa")?.extract::<f64>()?.to_bits(), molecule.borrow(py).inner.tpsa().unwrap().to_bits());
                    }
                    assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);
                }
            }
            Ok(())
        }).unwrap();
    }

    #[test]
    fn hall_kier_shape_and_mqns_keep_native_bits_order_and_input_requirements() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            type Query = fn(&ck::Molecule) -> Result<f64, ck::DescriptorReadError>;
            let scalar: [(&str, Query); 5] = [
                ("hall_kier_alpha", ck::Molecule::hall_kier_alpha),
                ("kappa_1", ck::Molecule::kappa_1),
                ("kappa_2", ck::Molecule::kappa_2),
                ("kappa_3", ck::Molecule::kappa_3),
                ("phi", ck::Molecule::phi),
            ];
            for text in ["CC(O)C", "c1ccncc1", "*C"] {
                for prepared in [true, false] {
                    let native = ck::Molecule::from_smiles_with_params(
                        text,
                        &ck::SmilesParseParams {
                            sanitize: prepared,
                            remove_hs: prepared,
                            ..Default::default()
                        },
                    )
                    .unwrap();
                    let before = native.to_binary().unwrap();
                    let molecule = Py::new(py, crate::native::Molecule { inner: native })?;
                    for (name, run) in scalar {
                        let expected = run(&molecule.borrow(py).inner).unwrap();
                        let actual = molecule.bind(py).call_method0(name)?.extract::<f64>()?;
                        assert_eq!(actual.to_bits(), expected.to_bits());
                    }
                    let (expected_scalar, expected_rows) = molecule
                        .borrow(py)
                        .inner
                        .hall_kier_alpha_with_contributions()
                        .unwrap();
                    let (actual_scalar, actual_rows) = molecule
                        .bind(py)
                        .call_method0("hall_kier_alpha_with_contributions")?
                        .extract::<(f64, Vec<f64>)>()?;
                    assert_eq!(actual_scalar.to_bits(), expected_scalar.to_bits());
                    assert_eq!(
                        actual_rows.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                        expected_rows
                            .iter()
                            .map(|x| x.to_bits())
                            .collect::<Vec<_>>()
                    );
                    assert_eq!(actual_rows.len(), molecule.borrow(py).inner.num_atoms());
                    for force in [false, true] {
                        let expected = molecule.borrow(py).inner.mqns(force);
                        let keywords = pyo3::types::PyDict::new(py);
                        keywords.set_item("force", force)?;
                        let actual = molecule.bind(py).call_method("mqns", (), Some(&keywords));
                        if prepared {
                            let expected = expected.unwrap();
                            assert_eq!(expected.len(), 42);
                            assert_eq!(actual?.extract::<Vec<u32>>()?, expected);
                            if !force {
                                assert_eq!(
                                    molecule
                                        .bind(py)
                                        .call_method0("mqns")?
                                        .extract::<Vec<u32>>()?,
                                    expected
                                );
                            }
                        } else {
                            assert!(matches!(
                                expected,
                                Err(ck::DescriptorReadError::MissingPreparedValence)
                            ));
                            let error = actual.unwrap_err();
                            assert!(error.is_instance(py, <DescriptorReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
                            assert!(
                                error
                                    .value(py)
                                    .getattr("variant")?
                                    .extract::<DescriptorReadErrorKind>()?
                                    == DescriptorReadErrorKind::MissingPreparedValence
                            );
                            assert!(error.value(py).getattr("source")?.is_none());
                        }
                    }
                    assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);
                }
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn connectivity_producers_keep_native_bits_force_and_prepared_valence_policy() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            type Query = fn(&ck::Molecule) -> Result<f64, ck::DescriptorReadError>;
            type Forced = fn(&ck::Molecule, bool) -> Result<f64, ck::DescriptorReadError>;
            let degree: [(&str, Query); 2] = [
                ("chi_0", ck::Molecule::chi_0),
                ("chi_1", ck::Molecule::chi_1),
            ];
            let valence: [(&str, Query, Forced); 10] = [
                (
                    "chi_0_v",
                    ck::Molecule::chi_0_v,
                    ck::Molecule::chi_0_v_with_params,
                ),
                (
                    "chi_1_v",
                    ck::Molecule::chi_1_v,
                    ck::Molecule::chi_1_v_with_params,
                ),
                (
                    "chi_2_v",
                    ck::Molecule::chi_2_v,
                    ck::Molecule::chi_2_v_with_params,
                ),
                (
                    "chi_3_v",
                    ck::Molecule::chi_3_v,
                    ck::Molecule::chi_3_v_with_params,
                ),
                (
                    "chi_4_v",
                    ck::Molecule::chi_4_v,
                    ck::Molecule::chi_4_v_with_params,
                ),
                (
                    "chi_0_n",
                    ck::Molecule::chi_0_n,
                    ck::Molecule::chi_0_n_with_params,
                ),
                (
                    "chi_1_n",
                    ck::Molecule::chi_1_n,
                    ck::Molecule::chi_1_n_with_params,
                ),
                (
                    "chi_2_n",
                    ck::Molecule::chi_2_n,
                    ck::Molecule::chi_2_n_with_params,
                ),
                (
                    "chi_3_n",
                    ck::Molecule::chi_3_n,
                    ck::Molecule::chi_3_n_with_params,
                ),
                (
                    "chi_4_n",
                    ck::Molecule::chi_4_n,
                    ck::Molecule::chi_4_n_with_params,
                ),
            ];
            for prepared in [true, false] {
                let native = ck::Molecule::from_smiles_with_params(
                    "CC(O)C",
                    &ck::SmilesParseParams {
                        sanitize: prepared,
                        remove_hs: prepared,
                        ..Default::default()
                    },
                )
                .unwrap();
                let before = native.to_binary().unwrap();
                let molecule = Py::new(py, crate::native::Molecule { inner: native })?;
                for (name, run) in degree {
                    let expected = run(&molecule.borrow(py).inner).unwrap();
                    let actual = molecule.bind(py).call_method0(name)?.extract::<f64>()?;
                    assert_eq!(actual.to_bits(), expected.to_bits());
                }
                for (name, run, force_run) in valence {
                    let expected = run(&molecule.borrow(py).inner);
                    if prepared {
                        let expected = expected.unwrap();
                        assert_eq!(
                            molecule
                                .bind(py)
                                .call_method0(name)?
                                .extract::<f64>()?
                                .to_bits(),
                            expected.to_bits()
                        );
                        for force in [false, true] {
                            let expected = force_run(&molecule.borrow(py).inner, force).unwrap();
                            let keywords = pyo3::types::PyDict::new(py);
                            keywords.set_item("force", force)?;
                            let actual = molecule
                                .bind(py)
                                .call_method(name, (), Some(&keywords))?
                                .extract::<f64>()?;
                            assert_eq!(actual.to_bits(), expected.to_bits());
                            let actual = molecule
                                .bind(py)
                                .call_method1(format!("{name}_with_params"), (force,))?
                                .extract::<f64>()?;
                            assert_eq!(actual.to_bits(), expected.to_bits());
                        }
                    } else {
                        assert!(matches!(
                            expected,
                            Err(ck::DescriptorReadError::MissingPreparedValence)
                        ));
                        let error = molecule.bind(py).call_method0(name).unwrap_err();
                        assert!(error.is_instance(py, <DescriptorReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
                        assert!(
                            error
                                .value(py)
                                .getattr("variant")?
                                .extract::<DescriptorReadErrorKind>()?
                                == DescriptorReadErrorKind::MissingPreparedValence
                        );
                    }
                }
                for (name, run) in [
                    (
                        "chi_n_v",
                        ck::Molecule::chi_n_v_with_params
                            as fn(&ck::Molecule, u32, bool) -> Result<f64, ck::DescriptorReadError>,
                    ),
                    ("chi_n_n", ck::Molecule::chi_n_n_with_params),
                ] {
                    if prepared {
                        let expected = run(&molecule.borrow(py).inner, 3, false).unwrap();
                        assert_eq!(
                            molecule
                                .bind(py)
                                .call_method1(name, (3,))?
                                .extract::<f64>()?
                                .to_bits(),
                            expected.to_bits()
                        );
                        for force in [false, true] {
                            let expected = run(&molecule.borrow(py).inner, 3, force).unwrap();
                            let keywords = pyo3::types::PyDict::new(py);
                            keywords.set_item("force", force)?;
                            assert_eq!(
                                molecule
                                    .bind(py)
                                    .call_method(name, (3,), Some(&keywords))?
                                    .extract::<f64>()?
                                    .to_bits(),
                                expected.to_bits()
                            );
                            assert_eq!(
                                molecule
                                    .bind(py)
                                    .call_method1(format!("{name}_with_params"), (3, force))?
                                    .extract::<f64>()?
                                    .to_bits(),
                                expected.to_bits()
                            );
                        }
                    } else {
                        let error = molecule.bind(py).call_method1(name, (3,)).unwrap_err();
                        assert!(error.is_instance(py, <DescriptorReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
                    }
                }
                // Native weight memoization is allowed; authoritative value blocks are not changed.
                assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);
            }
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn actual_ring_count_producers_keep_native_empty_state_policy_and_molecule_value() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            type Read = fn(&ck::Molecule) -> Result<u32, ck::DescriptorReadError>;
            let queries: [(&str, Read); 11] = [
                ("num_rings", ck::Molecule::num_rings),
                ("num_heterocycles", ck::Molecule::num_heterocycles),
                ("num_aromatic_rings", ck::Molecule::num_aromatic_rings),
                ("num_saturated_rings", ck::Molecule::num_saturated_rings),
                ("num_aliphatic_rings", ck::Molecule::num_aliphatic_rings),
                (
                    "num_aromatic_heterocycles",
                    ck::Molecule::num_aromatic_heterocycles,
                ),
                (
                    "num_aromatic_carbocycles",
                    ck::Molecule::num_aromatic_carbocycles,
                ),
                (
                    "num_aliphatic_heterocycles",
                    ck::Molecule::num_aliphatic_heterocycles,
                ),
                (
                    "num_aliphatic_carbocycles",
                    ck::Molecule::num_aliphatic_carbocycles,
                ),
                (
                    "num_saturated_heterocycles",
                    ck::Molecule::num_saturated_heterocycles,
                ),
                (
                    "num_saturated_carbocycles",
                    ck::Molecule::num_saturated_carbocycles,
                ),
            ];
            for (sanitize, remove_hs) in [(true, true), (false, true), (false, false)] {
                for text in ["c1ccncc1", "C1CCCCC1", "CC"] {
                    let native = ck::Molecule::from_smiles_with_params(
                        text,
                        &ck::SmilesParseParams {
                            sanitize,
                            remove_hs,
                            ..Default::default()
                        },
                    )
                    .unwrap();
                    let before = native.to_binary().unwrap();
                    let expected = queries
                        .iter()
                        .map(|(_, run)| run(&native))
                        .collect::<Vec<_>>();
                    if !sanitize && !remove_hs {
                        assert!(matches!(
                            expected[0],
                            Err(ck::DescriptorReadError::MissingInitializedRings)
                        ));
                        assert!(expected[1..].iter().all(|value| matches!(value, Ok(0))));
                    }
                    let molecule = Py::new(py, crate::native::Molecule { inner: native })?;
                    for ((name, _), expected) in queries.iter().zip(expected) {
                        let actual = molecule.bind(py).call_method0(*name);
                        match expected {
                            Ok(value) => {
                                assert_eq!(actual?.extract::<u32>()?, value, "{name} {text}")
                            }
                            Err(expected) => {
                                let error = actual.unwrap_err();
                                assert!(error.is_instance(py, <DescriptorReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
                                assert_eq!(error.value(py).str()?.to_str()?, expected.to_string());
                                assert!(error.value(py).getattr("source")?.is_none());
                                assert!(
                                    error
                                        .value(py)
                                        .getattr("variant")?
                                        .extract::<DescriptorReadErrorKind>()?
                                        == DescriptorReadErrorKind::MissingInitializedRings
                                );
                            }
                        }
                    }
                    assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);
                }
            }
            Ok(())
        })
        .unwrap();
    }

    use crate::property_numeric_errors::PropertyIntReadError;
    use crate::property_values::PropertyValueKind;

    #[test]
    fn descriptor_read_failures_keep_prepared_state_boundary_and_complete_search_cause() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let molecule = ck::Molecule::from_smiles_with_params(
                "c1ccccc1",
                &ck::SmilesParseParams {
                    sanitize: false,
                    remove_hs: false,
                    ..Default::default()
                },
            )
            .unwrap();
            let native = molecule.num_rings().unwrap_err();
            assert!(matches!(
                native,
                ck::DescriptorReadError::MissingInitializedRings
            ));
            assert!(native.source().is_none());
            let message = native.to_string();
            let error = DescriptorReadError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, message);
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<DescriptorReadErrorKind>()?
                    == DescriptorReadErrorKind::MissingInitializedRings
            );
            assert!(error.value(py).getattr("source")?.is_none());

            let native = ck::DescriptorReadError::Algorithm {
                source: ck::DescriptorError::Search {
                    function: "num_hba",
                    source: ck::DescriptorSearchCause::Match(
                        ck::SubstructMatchError::PropertyInteger {
                            property: "count",
                            source: ck::PropertyIntReadError::InvalidKind {
                                kind: ck::PropertyValueKind::Bool,
                            },
                        },
                    ),
                },
            };
            assert!(native.source().unwrap().is::<ck::DescriptorError>());
            let error = DescriptorReadError::from_rust(native);
            let algorithm = error.value(py).getattr("source")?;
            assert!(algorithm.is(&error.value(py).getattr("__cause__")?));
            assert!(algorithm.is_instance(<DescriptorError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert_eq!(
                algorithm.getattr("function")?.extract::<String>()?,
                "num_hba"
            );
            assert!(algorithm.setattr("function", "other").is_err());
            let search = algorithm.getattr("source")?;
            assert!(search.is_instance(<DescriptorSearchCause as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            let matching = search.getattr("source")?;
            assert!(matching.is_instance(<SubstructMatchError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            let numeric = matching.getattr("source")?;
            assert!(numeric.is_instance(<PropertyIntReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(
                numeric.getattr("kind")?.extract::<PropertyValueKind>()? == PropertyValueKind::Bool
            );

            let error = DescriptorError::from_rust(ck::DescriptorError::InvalidConnectivityPath {
                function: "chi",
                expected_rows: usize::MAX,
                actual_rows: None,
            });
            assert_eq!(
                error
                    .value(py)
                    .getattr("expected_rows")?
                    .extract::<usize>()?,
                usize::MAX
            );
            assert!(error.value(py).getattr("actual_rows")?.is_none());
            assert!(error.cause(py).is_none());
            Ok(())
        })
        .unwrap();
    }
}
