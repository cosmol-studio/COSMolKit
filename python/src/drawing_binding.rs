//! One canonical Molecule class, retaining the delivered drawing projections.

use crate::canonical_values::*;
use ::cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};

pyo3::create_exception!(cosmolkit, DrawingError, PyValueError);

pyo3::create_exception!(cosmolkit, OperationError, PyValueError);
pyo3::create_exception!(cosmolkit, DrawingWriteError, pyo3::exceptions::PyOSError);

// Transport actual source messages only; Python cannot retain Rust downcast identity.
fn source_pyerr(py: Python<'_>, source: &(dyn std::error::Error + 'static)) -> PyErr {
    let error = PyValueError::new_err(source.to_string());
    error.set_cause(py, source.source().map(|cause| source_pyerr(py, cause)));
    error
}

fn drawing_pyerr(py: Python<'_>, source: ck::DrawingError) -> PyErr {
    use ck::DrawingError as E;
    let error = DrawingError::new_err(source.to_string());
    let kind = match &source {
        E::Property(..) => "Property",
        E::Topology(..) => "Topology",
        E::Coordinates(..) => "Coordinates",
        E::Mapping(..) => "Mapping",
        E::Kekulize(..) => "Kekulize",
        E::Hydrogen(..) => "Hydrogen",
        E::Wedge(..) => "Wedge",
        E::Valence(..) => "Valence",
        E::CoordinateGeneration(..) => "CoordinateGeneration",
        E::SvgParse(..) => "SvgParse",
        E::PngEncode(..) => "PngEncode",
        E::StateRows { .. } => "StateRows",
        E::HydrogenAppend { .. } => "HydrogenAppend",
        E::InvalidDimensions { .. } => "InvalidDimensions",
        E::PixmapAllocation { .. } => "PixmapAllocation",
    };
    let attributes = || -> PyResult<()> {
        let value = error.value(py);
        value.setattr("domain", "drawing")?;
        value.setattr("kind", kind)?;
        match &source {
            E::InvalidDimensions { width, height } | E::PixmapAllocation { width, height } => {
                value.setattr("width", *width)?;
                value.setattr("height", *height)?;
            }
            E::StateRows {
                field,
                actual,
                expected,
            } => {
                value.setattr("field", *field)?;
                value.setattr("actual", *actual)?;
                value.setattr("expected", *expected)?;
            }
            E::HydrogenAppend { row, reason } => {
                value.setattr("row", *row)?;
                value.setattr("reason", *reason)?;
            }
            _ => {}
        }
        Ok(())
    };
    if let Err(attribute_error) = attributes() {
        return attribute_error;
    }
    error.set_cause(
        py,
        std::error::Error::source(&source).map(|cause| source_pyerr(py, cause)),
    );
    error
}

fn operation_pyerr(py: Python<'_>, source: ck::OperationError) -> PyErr {
    use ck::OperationError as E;
    let kind = match &source {
        E::UnsupportedFeature { .. } => "UnsupportedFeature",
        E::Unsupported { .. } => "Unsupported",
        E::OutputMismatch { .. } => "OutputMismatch",
        E::AccessDenied { .. } => "AccessDenied",
        E::BlockCheckedOut { .. } => "BlockCheckedOut",
        E::BlockNotCheckedOut { .. } => "BlockNotCheckedOut",
        E::IncompleteCommit { .. } => "IncompleteCommit",
        E::TopologyEditContract { .. } => "TopologyEditContract",
        E::MappingContract { .. } => "MappingContract",
        E::InvalidTopologyMapping { .. } => "InvalidTopologyMapping",
        E::AutoRemapContract { .. } => "AutoRemapContract",
        E::OperationContract { .. } => "OperationContract",
        E::SemanticPreconditionContract { .. } => "SemanticPreconditionContract",
        E::CoordinateAppendRequiresValues { .. } => "CoordinateAppendRequiresValues",
        E::DerivedEffectContract { .. } => "DerivedEffectContract",
        E::CipStateContract { .. } => "CipStateContract",
        E::InvalidPropertyList { .. } => "InvalidPropertyList",
        E::InvalidDerivedCache { .. } => "InvalidDerivedCache",
        E::InvalidAlgorithmResult { .. } => "InvalidAlgorithmResult",
        E::Algorithm { .. } => "Algorithm",
        E::InvalidTopology(..) => "InvalidTopology",
        E::InvalidTopologyEdit(..) => "InvalidTopologyEdit",
        E::InvalidCoordinates(..) => "InvalidCoordinates",
        E::InvalidProperty(..) => "InvalidProperty",
        E::Valence(..) => "Valence",
        E::Radical(..) => "Radical",
        E::Rings(..) => "Rings",
        E::PotentialStereo(..) => "PotentialStereo",
        E::Stereo(..) => "Stereo",
        E::CipLabeler(..) => "CipLabeler",
        E::Transform(..) => "Transform",
        E::Coordinate2D(..) => "Coordinate2D",
        E::Kekulize(..) => "Kekulize",
        E::Aromaticity(..) => "Aromaticity",
        E::Sanitize(..) => "Sanitize",
        E::Hydrogen(..) => "Hydrogen",
        E::UffOptimization(..) => "UffOptimization",
    };
    let error = OperationError::new_err(source.to_string());
    if let Err(attribute_error) = error
        .value(py)
        .setattr("domain", "operation")
        .and_then(|()| error.value(py).setattr("kind", kind))
    {
        return attribute_error;
    }
    error.set_cause(
        py,
        std::error::Error::source(&source).map(|cause| source_pyerr(py, cause)),
    );
    error
}

fn expand_user_path(path: &str) -> PyResult<std::path::PathBuf> {
    if path == "~" || path.starts_with("~/") {
        let home = std::env::var_os("HOME")
            .ok_or_else(|| PyValueError::new_err("cannot expand '~': HOME is not set"))?;
        let mut expanded = std::path::PathBuf::from(home);
        if let Some(rest) = path.strip_prefix("~/") {
            expanded.push(rest);
        }
        Ok(expanded)
    } else {
        Ok(std::path::PathBuf::from(path))
    }
}

fn drawing_write_pyerr(py: Python<'_>, source: ck::DrawingWriteError) -> PyErr {
    match source {
        ck::DrawingWriteError::Drawing(error) => drawing_pyerr(py, error),
        ck::DrawingWriteError::Io { path, source } => {
            let error = DrawingWriteError::new_err((
                source.raw_os_error(),
                source.to_string(),
                path.to_string_lossy().into_owned(),
            ));
            if let Err(attribute_error) = error
                .value(py)
                .setattr("domain", "drawing")
                .and_then(|()| error.value(py).setattr("kind", "Io"))
            {
                return attribute_error;
            }
            error
        }
    }
}

/// Immutable detached parameters projected from the public facade.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
struct Coordinate2DParams {
    inner: ck::Coordinate2DParams,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl Coordinate2DParams {
    #[new]
    #[pyo3(signature = (coordinate_map=None, *, canonical_orientation=false, clear_existing_2d=true, flips_per_sample=0, samples=0, sample_seed=0, permute_degree_four=false, force_rdkit=false, use_ring_templates=false))]
    fn new(
        coordinate_map: Option<std::collections::BTreeMap<usize, [f64; 2]>>,
        canonical_orientation: bool,
        clear_existing_2d: bool,
        flips_per_sample: u32,
        samples: u32,
        sample_seed: i32,
        permute_degree_four: bool,
        force_rdkit: bool,
        use_ring_templates: bool,
    ) -> Self {
        Self {
            inner: ck::Coordinate2DParams {
                coordinate_map: coordinate_map.unwrap_or_default(),
                canonical_orientation,
                clear_existing_2d,
                flips_per_sample,
                samples,
                sample_seed,
                permute_degree_four,
                force_rdkit,
                use_ring_templates,
            },
        }
    }

    #[getter]
    fn coordinate_map(&self) -> std::collections::BTreeMap<usize, [f64; 2]> {
        self.inner.coordinate_map.clone()
    }

    #[getter]
    fn canonical_orientation(&self) -> bool {
        self.inner.canonical_orientation
    }

    #[getter]
    fn clear_existing_2d(&self) -> bool {
        self.inner.clear_existing_2d
    }

    #[getter]
    fn flips_per_sample(&self) -> u32 {
        self.inner.flips_per_sample
    }

    #[getter]
    fn samples(&self) -> u32 {
        self.inner.samples
    }

    #[getter]
    fn sample_seed(&self) -> i32 {
        self.inner.sample_seed
    }

    #[getter]
    fn permute_degree_four(&self) -> bool {
        self.inner.permute_degree_four
    }

    #[getter]
    fn force_rdkit(&self) -> bool {
        self.inner.force_rdkit
    }

    #[getter]
    fn use_ring_templates(&self) -> bool {
        self.inner.use_ring_templates
    }
}

/// Immutable source definition selector, projected from the public facade.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", frozen, eq)]
#[derive(Clone, Copy, PartialEq)]
enum RotatableBondsOptions {
    Default,
    NonStrict,
    Strict,
    StrictLinkages,
}
impl RotatableBondsOptions {
    fn canonical(self) -> ck::RotatableBondsOptions {
        match self {
            Self::Default => ck::RotatableBondsOptions::Default,
            Self::NonStrict => ck::RotatableBondsOptions::NonStrict,
            Self::Strict => ck::RotatableBondsOptions::Strict,
            Self::StrictLinkages => ck::RotatableBondsOptions::StrictLinkages,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
struct CrippenTotals {
    inner: ck::CrippenTotals,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl CrippenTotals {
    #[getter]
    fn logp(&self) -> f64 {
        self.inner.logp
    }
    #[getter]
    fn molar_refractivity(&self) -> f64 {
        self.inner.molar_refractivity
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
struct LabuteAsaContributions {
    inner: ck::LabuteAsaContributions,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl LabuteAsaContributions {
    #[getter]
    fn asa(&self) -> f64 {
        self.inner.asa
    }
    #[getter]
    fn atom_contributions(&self) -> Vec<f64> {
        self.inner.atom_contributions.clone()
    }
    #[getter]
    fn hydrogen_contribution(&self) -> f64 {
        self.inner.hydrogen_contribution
    }
}
/// Python ownership wraps the ONE live runtime value, not detached chemistry.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
struct Molecule {
    inner: ck::Molecule,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl Molecule {
    /// Descriptor query through the canonical public Rust method.
    fn num_amide_bonds(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_amide_bonds()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Descriptor query through the canonical public Rust method.
    fn num_spiro_atoms(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_spiro_atoms()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Descriptor query through the canonical public Rust method.
    fn num_bridgehead_atoms(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_bridgehead_atoms()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Descriptor query through the canonical public Rust method.
    fn num_atom_stereo_centers(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_atom_stereo_centers()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Descriptor query through the canonical public Rust method.
    fn num_unspecified_atom_stereo_centers(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_unspecified_atom_stereo_centers()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Descriptor query through the canonical public Rust method.
    fn num_rotatable_bonds(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_rotatable_bonds()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Descriptor query through the canonical public Rust method.
    fn num_rotatable_bonds_with_params(
        &self,
        py: Python<'_>,
        params: RotatableBondsOptions,
    ) -> PyResult<u32> {
        self.inner
            .num_rotatable_bonds_with_params(&params.canonical())
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Descriptor query through the canonical public Rust method.
    fn molecular_weight_with_params(&self, py: Python<'_>, only_heavy: bool) -> PyResult<f64> {
        self.inner
            .molecular_weight_with_params(only_heavy)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Descriptor query through the canonical public Rust method.
    fn exact_molecular_weight_with_params(
        &self,
        py: Python<'_>,
        only_heavy: bool,
    ) -> PyResult<f64> {
        self.inner
            .exact_molecular_weight_with_params(only_heavy)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Descriptor query through the canonical public Rust method.
    fn molecular_formula_with_params(
        &self,
        py: Python<'_>,
        separate_isotopes: bool,
        abbreviate_h_isotopes: bool,
    ) -> PyResult<String> {
        self.inner
            .molecular_formula_with_params(separate_isotopes, abbreviate_h_isotopes)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Descriptor query through the canonical public Rust method.
    fn molecular_weight(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .molecular_weight()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Descriptor query through the canonical public Rust method.
    fn exact_molecular_weight(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .exact_molecular_weight()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Descriptor query through the canonical public Rust method.
    fn molecular_formula(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .molecular_formula()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Experimental read-only query using the canonical public facade.
    /// Read-only descriptor through the canonical public Rust method.
    fn num_heavy_atoms(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_heavy_atoms()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Read-only descriptor through the canonical public Rust method.
    fn total_atom_count(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .total_atom_count()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Read-only descriptor through the canonical public Rust method.
    fn lipinski_hba(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .lipinski_hba()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Read-only descriptor through the canonical public Rust method.
    fn lipinski_hbd(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .lipinski_hbd()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Read-only descriptor through the canonical public Rust method.
    fn fraction_csp3(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .fraction_csp3()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Read-only descriptor through the canonical public Rust method.
    fn num_heteroatoms(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_heteroatoms()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Read-only descriptor through the canonical public Rust method.
    fn num_hba(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_hba()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Read-only descriptor through the canonical public Rust method.
    fn num_hbd(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_hbd()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Read-only descriptor through the canonical public Rust method.
    fn num_rings(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_rings()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Read-only descriptor through the canonical public Rust method.
    fn num_heterocycles(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_heterocycles()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Read-only descriptor through the canonical public Rust method.
    fn num_aromatic_rings(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_aromatic_rings()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Read-only descriptor through the canonical public Rust method.
    fn num_saturated_rings(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_saturated_rings()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Read-only descriptor through the canonical public Rust method.
    fn num_aliphatic_rings(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_aliphatic_rings()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Read-only descriptor through the canonical public Rust method.
    fn num_aromatic_heterocycles(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_aromatic_heterocycles()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Read-only descriptor through the canonical public Rust method.
    fn num_aromatic_carbocycles(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_aromatic_carbocycles()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Read-only descriptor through the canonical public Rust method.
    fn num_aliphatic_heterocycles(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_aliphatic_heterocycles()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Read-only descriptor through the canonical public Rust method.
    fn num_aliphatic_carbocycles(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_aliphatic_carbocycles()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Read-only descriptor through the canonical public Rust method.
    fn num_saturated_heterocycles(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_saturated_heterocycles()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Read-only descriptor through the canonical public Rust method.
    fn num_saturated_carbocycles(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_saturated_carbocycles()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Degree-based Chi0 through the canonical public Rust method.
    fn chi_0(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_0()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Degree-based Chi1 through the canonical public Rust method.
    fn chi_1(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_1()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    fn hall_kier_alpha(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .hall_kier_alpha()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Experimental read-only query using the canonical public facade.
    fn hall_kier_alpha_with_contributions(&self, py: Python<'_>) -> PyResult<(f64, Vec<f64>)> {
        self.inner
            .hall_kier_alpha_with_contributions()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Experimental read-only query using the canonical public facade.
    fn kappa_1(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .kappa_1()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Experimental read-only query using the canonical public facade.
    fn kappa_2(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .kappa_2()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Experimental read-only query using the canonical public facade.
    fn kappa_3(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .kappa_3()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Experimental read-only query using the canonical public facade.
    fn phi(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .phi()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Experimental read-only query using the canonical public facade.
    /// All 42 components; the source ignores force and performs no cache write.
    #[pyo3(signature = (force=false))]
    fn mqns(&self, py: Python<'_>, force: bool) -> PyResult<Vec<u32>> {
        self.inner
            .mqns(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Experimental read-only query using the canonical public facade.
    /// Recomputes only; no vector-property caching or force parameter.
    /// Requires existing prepared valence only; missing/reset rings succeed.
    fn chi_0_v(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_0_v()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Experimental read-only query using the canonical public facade.
    /// Recomputes only; no vector-property caching or force parameter.
    /// Requires existing prepared valence only; missing/reset rings succeed.
    fn chi_1_v(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_1_v()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Experimental read-only query using the canonical public facade.
    /// Recomputes only; no vector-property caching or force parameter.
    /// Requires existing prepared valence only; missing/reset rings succeed.
    fn chi_2_v(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_2_v()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Experimental read-only query using the canonical public facade.
    /// Recomputes only; no vector-property caching or force parameter.
    /// Requires existing prepared valence only; missing/reset rings succeed.
    fn chi_3_v(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_3_v()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Experimental read-only query using the canonical public facade.
    /// Recomputes only; no vector-property caching or force parameter.
    /// Requires existing prepared valence only; missing/reset rings succeed.
    fn chi_4_v(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_4_v()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Experimental read-only query using the canonical public facade.
    /// Recomputes only; no vector-property caching or force parameter.
    /// Requires existing prepared valence only; missing/reset rings succeed.
    #[pyo3(signature = (order))]
    fn chi_n_v(&self, py: Python<'_>, order: u32) -> PyResult<f64> {
        self.inner
            .chi_n_v(order)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Experimental read-only query using the canonical public facade.
    /// Recomputes only; no vector-property caching or force parameter.
    /// Requires existing prepared valence only; missing/reset rings succeed.
    fn chi_0_n(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_0_n()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Experimental read-only query using the canonical public facade.
    /// Recomputes only; no vector-property caching or force parameter.
    /// Requires existing prepared valence only; missing/reset rings succeed.
    fn chi_1_n(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_1_n()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Experimental read-only query using the canonical public facade.
    /// Recomputes only; no vector-property caching or force parameter.
    /// Requires existing prepared valence only; missing/reset rings succeed.
    fn chi_2_n(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_2_n()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Experimental read-only query using the canonical public facade.
    /// Recomputes only; no vector-property caching or force parameter.
    /// Requires existing prepared valence only; missing/reset rings succeed.
    fn chi_3_n(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_3_n()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Experimental read-only query using the canonical public facade.
    /// Recomputes only; no vector-property caching or force parameter.
    /// Requires existing prepared valence only; missing/reset rings succeed.
    fn chi_4_n(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_4_n()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Experimental read-only query using the canonical public facade.
    /// Recomputes only; no vector-property caching or force parameter.
    /// Requires existing prepared valence only; missing/reset rings succeed.
    #[pyo3(signature = (order))]
    fn chi_n_n(&self, py: Python<'_>, order: u32) -> PyResult<f64> {
        self.inner
            .chi_n_n(order)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    #[staticmethod]
    fn new() -> Self {
        Self {
            inner: ck::Molecule::new(),
        }
    }

    #[staticmethod]
    fn from_smiles(py: Python<'_>, smiles: &str) -> PyResult<Self> {
        ck::Molecule::from_smiles(smiles)
            .map(|inner| Self { inner })
            .map_err(|error| smiles_pyerr(py, error))
    }

    #[staticmethod]
    fn from_smiles_with_params(
        py: Python<'_>,
        input: &str,
        params: &SmilesParseParams,
    ) -> PyResult<Self> {
        ck::Molecule::from_smiles_with_params(input, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| smiles_pyerr(py, e))
    }

    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }

    fn num_bonds(&self) -> usize {
        self.inner.num_bonds()
    }

    fn to_smiles(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .to_smiles()
            .map_err(|error| smiles_write_pyerr(py, error))
    }

    fn to_smiles_with_params(
        &self,
        py: Python<'_>,
        params: &SmilesWriteParams,
    ) -> PyResult<String> {
        self.inner
            .to_smiles_with_params(&params.inner)
            .map_err(|e| smiles_write_pyerr(py, e))
    }

    fn morgan_fingerprint(&self, py: Python<'_>) -> PyResult<Fingerprint> {
        self.inner
            .morgan_fingerprint()
            .map(|inner| Fingerprint { inner })
            .map_err(|e| morgan_pyerr(py, e))
    }

    fn morgan_sparse_fingerprint(&self, py: Python<'_>) -> PyResult<SparseBitFingerprint> {
        self.inner
            .morgan_sparse_fingerprint()
            .map(|inner| SparseBitFingerprint { inner })
            .map_err(|e| morgan_pyerr(py, e))
    }

    fn morgan_sparse_count_fingerprint(&self, py: Python<'_>) -> PyResult<SparseCountFingerprint> {
        self.inner
            .morgan_sparse_count_fingerprint()
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(|e| morgan_pyerr(py, e))
    }

    fn morgan_count_fingerprint(&self, py: Python<'_>) -> PyResult<SparseCountFingerprint32> {
        self.inner
            .morgan_count_fingerprint()
            .map(|inner| SparseCountFingerprint32 { inner })
            .map_err(|e| morgan_pyerr(py, e))
    }

    fn coordinates_2d(&self) -> Option<Vec<[f64; 2]>> {
        // Python receives an owned copy and cannot mutate the runtime block.
        self.inner.coordinates_2d().map(<[_]>::to_vec)
    }

    fn has_2d_coordinates(&self) -> bool {
        self.inner.has_2d_coordinates()
    }

    fn compute_2d_coordinates_(&mut self, py: Python<'_>) -> PyResult<()> {
        self.inner
            .compute_2d_coordinates_()
            .map_err(|error| operation_pyerr(py, error))
    }

    fn compute_2d_coordinates_with_params_(
        &mut self,
        py: Python<'_>,
        params: &Coordinate2DParams,
    ) -> PyResult<()> {
        self.inner
            .compute_2d_coordinates_with_params_(&params.inner)
            .map_err(|error| operation_pyerr(py, error))
    }

    fn with_2d_coordinates(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .with_2d_coordinates()
            .map(|inner| Self { inner })
            .map_err(|error| operation_pyerr(py, error))
    }

    /// Return the existing Rust sanitization operation's new molecule value.
    /// Its runtime owns descriptor-cache clearing and atomic commit.
    fn sanitize(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .sanitize()
            .map(|inner| Self { inner })
            .map_err(|error| operation_pyerr(py, error))
    }

    fn with_2d_coordinates_with_params(
        &self,
        py: Python<'_>,
        params: &Coordinate2DParams,
    ) -> PyResult<Self> {
        self.inner
            .with_2d_coordinates_with_params(&params.inner)
            .map(|inner| Self { inner })
            .map_err(|error| operation_pyerr(py, error))
    }

    /// Experimental SVG; required dimensions follow the current registry.
    #[pyo3(signature = (width, height))]
    fn to_svg(&self, py: Python<'_>, width: u32, height: u32) -> PyResult<String> {
        self.inner
            .to_svg(width, height)
            .map_err(|error| drawing_pyerr(py, error))
    }

    /// Experimental PNG bytes from the same public drawing owner.
    #[gen_stub(override_return_type(type_repr = "builtins.bytes", imports = ("builtins")))]
    #[pyo3(signature = (width, height))]
    fn to_png<'py>(
        &self,
        py: Python<'py>,
        width: u32,
        height: u32,
    ) -> PyResult<Bound<'py, PyBytes>> {
        let png = self
            .inner
            .to_png(width, height)
            .map_err(|error| drawing_pyerr(py, error))?;
        Ok(PyBytes::new(py, &png))
    }

    #[pyo3(signature = (path, width, height))]
    fn write_svg(&self, py: Python<'_>, path: &str, width: u32, height: u32) -> PyResult<()> {
        let path = expand_user_path(path)?;
        self.inner
            .write_svg(&path, width, height)
            .map_err(|error| drawing_write_pyerr(py, error))
    }

    #[pyo3(signature = (path, width, height))]
    fn write_png(&self, py: Python<'_>, path: &str, width: u32, height: u32) -> PyResult<()> {
        let path = expand_user_path(path)?;
        self.inner
            .write_png(&path, width, height)
            .map_err(|error| drawing_write_pyerr(py, error))
    }
    fn crippen_descriptors(&self, py: Python<'_>) -> PyResult<CrippenTotals> {
        let result = self
            .inner
            .crippen_descriptors()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(CrippenTotals { inner: result })
    }

    fn labute_asa(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .labute_asa()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn labute_asa_contributions(&self, py: Python<'_>) -> PyResult<LabuteAsaContributions> {
        let result = self
            .inner
            .labute_asa_contributions()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(LabuteAsaContributions { inner: result })
    }

    fn tpsa(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .tpsa()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn slogp_vsa(&self, py: Python<'_>) -> PyResult<Vec<f64>> {
        let result = self
            .inner
            .slogp_vsa()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn smr_vsa(&self, py: Python<'_>) -> PyResult<Vec<f64>> {
        let result = self
            .inner
            .smr_vsa()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn slogp_vsa_1(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_1()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn slogp_vsa_2(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_2()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn slogp_vsa_3(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_3()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn slogp_vsa_4(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_4()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn slogp_vsa_5(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_5()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn slogp_vsa_6(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_6()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn slogp_vsa_7(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_7()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn slogp_vsa_8(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_8()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn slogp_vsa_9(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_9()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn slogp_vsa_10(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_10()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn slogp_vsa_11(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_11()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn slogp_vsa_12(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_12()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn smr_vsa_1(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_1()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn smr_vsa_2(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_2()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn smr_vsa_3(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_3()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn smr_vsa_4(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_4()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn smr_vsa_5(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_5()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn smr_vsa_6(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_6()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn smr_vsa_7(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_7()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn smr_vsa_8(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_8()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn smr_vsa_9(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_9()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn smr_vsa_10(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_10()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    #[pyo3(signature = (include_hydrogens, force))]
    fn crippen_descriptors_with_params(
        &self,
        py: Python<'_>,
        include_hydrogens: bool,
        force: bool,
    ) -> PyResult<CrippenTotals> {
        let result = self
            .inner
            .crippen_descriptors_with_params(include_hydrogens, force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(CrippenTotals { inner: result })
    }

    #[pyo3(signature = (include_hydrogens, force))]
    fn labute_asa_with_params(
        &self,
        py: Python<'_>,
        include_hydrogens: bool,
        force: bool,
    ) -> PyResult<f64> {
        let result = self
            .inner
            .labute_asa_with_params(include_hydrogens, force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    #[pyo3(signature = (include_hydrogens, force))]
    fn labute_asa_contributions_with_params(
        &self,
        py: Python<'_>,
        include_hydrogens: bool,
        force: bool,
    ) -> PyResult<LabuteAsaContributions> {
        let result = self
            .inner
            .labute_asa_contributions_with_params(include_hydrogens, force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(LabuteAsaContributions { inner: result })
    }

    #[pyo3(signature = (include_sulfur_phosphorus, force))]
    fn tpsa_with_params(
        &self,
        py: Python<'_>,
        include_sulfur_phosphorus: bool,
        force: bool,
    ) -> PyResult<f64> {
        let result = self
            .inner
            .tpsa_with_params(include_sulfur_phosphorus, force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    #[pyo3(signature = (bins, force))]
    fn slogp_vsa_with_params(
        &self,
        py: Python<'_>,
        bins: Option<Vec<f64>>,
        force: bool,
    ) -> PyResult<Vec<f64>> {
        let result = self
            .inner
            .slogp_vsa_with_params(bins.as_deref(), force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    #[pyo3(signature = (bins, force))]
    fn smr_vsa_with_params(
        &self,
        py: Python<'_>,
        bins: Option<Vec<f64>>,
        force: bool,
    ) -> PyResult<Vec<f64>> {
        let result = self
            .inner
            .smr_vsa_with_params(bins.as_deref(), force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    fn qed(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .qed()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }
    #[pyo3(signature = (force))]
    fn chi_0_v_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_0_v_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    #[pyo3(signature = (force))]
    fn chi_1_v_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_1_v_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    #[pyo3(signature = (force))]
    fn chi_2_v_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_2_v_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    #[pyo3(signature = (force))]
    fn chi_3_v_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_3_v_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    #[pyo3(signature = (force))]
    fn chi_4_v_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_4_v_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    #[pyo3(signature = (order, force))]
    fn chi_n_v_with_params(&self, py: Python<'_>, order: u32, force: bool) -> PyResult<f64> {
        self.inner
            .chi_n_v_with_params(order, force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    #[pyo3(signature = (force))]
    fn chi_0_n_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_0_n_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    #[pyo3(signature = (force))]
    fn chi_1_n_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_1_n_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    #[pyo3(signature = (force))]
    fn chi_2_n_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_2_n_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    #[pyo3(signature = (force))]
    fn chi_3_n_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_3_n_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    #[pyo3(signature = (force))]
    fn chi_4_n_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_4_n_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    #[pyo3(signature = (order, force))]
    fn chi_n_n_with_params(&self, py: Python<'_>, order: u32, force: bool) -> PyResult<f64> {
        self.inner
            .chi_n_n_with_params(order, force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }
}

#[pymodule]
fn cosmolkit(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("__version__", ck::version())?;
    module.add(
        "_binding_profile",
        if cfg!(feature = "drawing-bindings") {
            "drawing-bindings"
        } else {
            "canonical-bootstrap"
        },
    )?;
    module.add("DrawingError", module.py().get_type::<DrawingError>())?;
    module.add(
        "DrawingWriteError",
        module.py().get_type::<DrawingWriteError>(),
    )?;
    module.add("OperationError", module.py().get_type::<OperationError>())?;
    module.add_class::<Molecule>()?;
    module.add_class::<RotatableBondsOptions>()?;
    module.add_class::<CrippenTotals>()?;
    module.add_class::<LabuteAsaContributions>()?;
    module.add_class::<Coordinate2DParams>()?;
    crate::canonical_descriptor_binding::register(module)?;
    crate::canonical_values::register(module)?;
    Ok(())
}
