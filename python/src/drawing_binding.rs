//! One canonical Molecule class, retaining the delivered drawing projections.

use crate::alignment_binding::*;
use crate::canonical_fingerprint_values::{
    AtomPairFingerprintParams, FingerprintAdditionalOutput, LegacyTopologicalTorsionParams,
    MorganFingerprintParams, TopologicalTorsionCallParams, TopologicalTorsionFingerprintGenerator,
    TopologicalTorsionFingerprintParams,
};
use crate::canonical_values::*;
use crate::conformer_binding::{EmbedMoleculeResult, EmbedMultipleConfsResult, EmbedParams};
use ::cosmolkit as ck;
use numpy::{IntoPyArray, ndarray::Array2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};

pyo3::create_exception!(
    cosmolkit,
    DrawingError,
    PyValueError,
    "Molecular depiction or image rendering failed for the graph or drawing options."
);

pyo3::create_exception!(
    cosmolkit,
    OperationError,
    PyValueError,
    "A molecular operation failed validation or execution without committing partial changes."
);
pyo3::create_exception!(
    cosmolkit,
    DrawingWriteError,
    pyo3::exceptions::PyOSError,
    "A rendered molecular image could not be written to the requested path."
);

// Preserve each recognized canonical domain cause through the shared projection.
pub(crate) fn source_pyerr(py: Python<'_>, source: &(dyn std::error::Error + 'static)) -> PyErr {
    crate::canonical_values::source_pyerr(py, source)
}

pub(crate) fn drawing_pyerr(
    py: Python<'_>,
    source: impl std::borrow::Borrow<ck::DrawingError>,
) -> PyErr {
    let source = source.borrow();
    use ck::DrawingError as E;
    let error = DrawingError::new_err(source.to_string());
    let kind = match source {
        E::Property(..) => "Property",
        E::Topology(..) => "Topology",
        E::Coordinates(..) => "Coordinates",
        E::Mapping(..) => "Mapping",
        E::Kekulize(..) => "Kekulize",
        E::Hydrogen(..) => "Hydrogen",
        E::Wedge(..) => "Wedge",
        E::Valence(..) => "Valence",
        E::CoordinateGeneration(..) => "CoordinateGeneration",
        E::SvgTextProjection(..) => "SvgTextProjection",
        E::PropertyString(..) => "PropertyString",
        E::VariationArray(..) => "VariationArray",
        E::VariationIndex { .. } => "VariationIndex",
        E::DataFieldDouble { .. } => "DataFieldDouble",
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
        match source {
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
        std::error::Error::source(source).map(|cause| source_pyerr(py, cause)),
    );
    error
}

pub(crate) fn operation_pyerr(
    py: Python<'_>,
    source: impl std::borrow::Borrow<ck::OperationError>,
) -> PyErr {
    let source = source.borrow();
    use ck::OperationError as E;
    let kind = match source {
        E::AtomPropertyIndex { .. } => "AtomPropertyIndex",
        E::ReservedAtomPropertyKey { .. } => "ReservedAtomPropertyKey",
        E::ReactionRun(..) => "ReactionRun",
        E::Scaffold(..) => "Scaffold",
        E::Fragments(..) => "Fragments",
        E::EmptyFragments => "EmptyFragments",
        E::ReactionApply(..) => "ReactionApply",
        E::Alignment(..) => "Alignment",
        E::Enumeration(..) => "Enumeration",
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
        E::AtomProperty(..) => "AtomProperty",
        E::BondProperty(..) => "BondProperty",
        E::InvalidReconstructionOrigin { .. } => "InvalidReconstructionOrigin",
        E::InvalidProperty(..) => "InvalidProperty",
        E::Valence(..) => "Valence",
        E::Radical(..) => "Radical",
        E::Rings(..) => "Rings",
        E::PotentialStereo(..) => "PotentialStereo",
        E::Stereo(..) => "Stereo",
        E::CipLabeler(..) => "CipLabeler",
        E::AtomCode(..) => "AtomCode",
        E::Transform(..) => "Transform",
        E::CoordinateInput(..) => "CoordinateInput",
        E::Coordinate2D(..) => "Coordinate2D",
        E::Kekulize(..) => "Kekulize",
        E::Aromaticity(..) => "Aromaticity",
        E::Sanitize(..) => "Sanitize",
        E::Hydrogen(..) => "Hydrogen",
        E::UffOptimization(..) => "UffOptimization",
        E::MmffOptimization(..) => "MmffOptimization",
        E::Tautomer(..) => "Tautomer",
        E::Conformer(..) => "Conformer",
    };
    let error = OperationError::new_err(source.to_string());
    if let Err(attribute_error) = error
        .value(py)
        .setattr("domain", "operation")
        .and_then(|()| error.value(py).setattr("kind", kind))
    {
        return attribute_error;
    }
    let context = match source {
        E::AtomPropertyIndex { atom, atom_count } => error
            .value(py)
            .setattr("atom_index", atom.index())
            .and_then(|()| error.value(py).setattr("atom_count", atom_count)),
        E::ReservedAtomPropertyKey { key } => error.value(py).setattr("key", key),
        _ => Ok(()),
    };
    if let Err(attribute_error) = context {
        return attribute_error;
    }
    error.set_cause(
        py,
        match source {
            E::ReactionRun(cause) => Some(crate::canonical_reaction::run_error(py, cause)),
            E::ReactionApply(cause) => Some(crate::canonical_reaction::apply_error(py, cause)),
            E::UffOptimization(cause) => Some(crate::uff_binding::optimization_pyerr(py, cause)),
            E::MmffOptimization(cause) => Some(crate::mmff_binding::optimization_pyerr(py, cause)),
            E::Tautomer(cause) => Some(crate::tautomer_binding::run_pyerr(py, cause)),
            E::Enumeration(cause) => Some(crate::canonical_stereoisomers::run_pyerr(py, cause)),
            E::Alignment(cause) => {
                Some(crate::alignment_binding::alignment_pyerr(py, cause.clone()))
            }
            E::PotentialStereo(cause) => {
                Some(crate::canonical_potential_stereo::error_pyerr(py, cause))
            }
            _ => std::error::Error::source(source).map(|cause| source_pyerr(py, cause)),
        },
    );
    error
}

fn expand_user_path(path: &str) -> PyResult<std::path::PathBuf> {
    crate::user_path::expand_user_path(path)
        .map_err(|source| PyValueError::new_err(source.to_string()))
}

pub(crate) fn drawing_write_pyerr(
    py: Python<'_>,
    source: impl std::borrow::Borrow<ck::DrawingWriteError>,
) -> PyErr {
    let source = source.borrow();
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

/// Writable configuration for 2D coordinate generation.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct Coordinate2DParams {
    pub(crate) inner: ck::Coordinate2DParams,
}

#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl Coordinate2DParams {
    /// Configure 2D coordinate generation; omitted fields use the defaults shown in the signature.
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

    /// Fixed coordinates indexed by atom index for coordinate generation.
    #[getter]
    fn coordinate_map(&self) -> std::collections::BTreeMap<usize, [f64; 2]> {
        self.inner.coordinate_map.clone()
    }

    /// Whether generated 2D coordinates are canonically oriented.
    #[getter]
    fn canonical_orientation(&self) -> bool {
        self.inner.canonical_orientation
    }

    /// Whether existing 2D conformers are cleared before storing generated coordinates.
    #[getter]
    fn clear_existing_2d(&self) -> bool {
        self.inner.clear_existing_2d
    }

    /// Number of rotatable-bond flips attempted per 2D layout sample.
    #[getter]
    fn flips_per_sample(&self) -> u32 {
        self.inner.flips_per_sample
    }

    /// Number of randomized 2D layout samples.
    #[getter]
    fn samples(&self) -> u32 {
        self.inner.samples
    }

    /// Random seed for sampled 2D coordinate generation.
    #[getter]
    fn sample_seed(&self) -> i32 {
        self.inner.sample_seed
    }

    /// Whether degree-four atom arrangements are permuted during 2D sampling.
    #[getter]
    fn permute_degree_four(&self) -> bool {
        self.inner.permute_degree_four
    }

    /// Whether to force the RDKit-style coordinate generation path.
    #[getter]
    fn force_rdkit(&self) -> bool {
        self.inner.force_rdkit
    }

    /// Whether ring templates are used for 2D layouts.
    #[getter]
    fn use_ring_templates(&self) -> bool {
        self.inner.use_ring_templates
    }
}

/// Selects which rotatable-bond definition the count functions use.
///
/// Mirrors RDKit ``NumRotatableBondsOptions`` (``Lipinski.h:39-43``). The C++
/// int discriminants (``Default = -1``, ``NonStrict = 0``, ``Strict = 1``,
/// ``StrictLinkages = 2``) are never observable through the source API (the
/// option is only compared for equality), so the Rust projection carries no
/// numeric discriminants.
///
/// The source also exposes a deprecated ``bool`` overload
/// (``Lipinski.cpp:185-187``) mapping ``strict == true`` to ``Strict`` and
/// ``false`` to ``NonStrict``. This binding exposes the enum options, not the
/// deprecated boolean overload.
///
/// Declared values: ``Default``, ``NonStrict``, ``Strict``, ``StrictLinkages``.
#[cosmolkit_macros::python_enum]
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

/// Wildman-Crippen logP and molar-refractivity scalar results.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
struct CrippenTotals {
    inner: ck::CrippenTotals,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl CrippenTotals {
    /// Wildman-Crippen octanol/water log partition coefficient.
    #[getter]
    fn logp(&self) -> f64 {
        self.inner.logp
    }
    /// Wildman-Crippen molar refractivity.
    #[getter]
    fn molar_refractivity(&self) -> f64 {
        self.inner.molar_refractivity
    }
}

/// Total Labute accessible surface area plus atom-indexed and hydrogen contributions.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
struct LabuteAsaContributions {
    inner: ck::LabuteAsaContributions,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl LabuteAsaContributions {
    /// Total Labute approximate accessible surface area.
    #[getter]
    fn asa(&self) -> f64 {
        self.inner.asa
    }
    /// Per-atom contributions in atom-index order.
    #[getter]
    fn atom_contributions(&self) -> Vec<f64> {
        self.inner.atom_contributions.clone()
    }
    /// Hydrogen contribution to the total Labute accessible surface area.
    #[getter]
    fn hydrogen_contribution(&self) -> f64 {
        self.inner.hydrogen_contribution
    }
}
/// Owned molecular graph with properties and separate 2D/3D conformers.
///
/// Queries do not change the graph. Value transformations return a new molecule;
/// methods ending in an underscore modify this object in place with copy-on-write
/// isolation from other molecules.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct Molecule {
    pub(crate) inner: ck::Molecule,
}
impl Molecule {
    pub(crate) fn from_inner(inner: ck::Molecule) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl Molecule {
    /// Parse InChI text into a new Molecule. sanitize and remove_hs default to True; params and keyword options are mutually exclusive. Invalid input raises InchiError.
    #[staticmethod]
    #[pyo3(signature=(text, params=None, *, sanitize=None, remove_hs=None))]
    fn from_inchi(
        py: Python<'_>,
        text: crate::text_input::TextInput<'_>,
        params: Option<&crate::canonical_inchi::InchiReadParams>,
        sanitize: Option<bool>,
        remove_hs: Option<bool>,
    ) -> PyResult<Self> {
        if params.is_some() && (sanitize.is_some() || remove_hs.is_some()) {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "params and keyword options are mutually exclusive",
            ));
        }
        let params = params.map(|p| p.inner).unwrap_or(ck::InchiReadParams {
            sanitize: sanitize.unwrap_or(true),
            remove_hs: remove_hs.unwrap_or(true),
        });
        ck::Molecule::from_inchi_with_params(&text.as_text()?, &params)
            .map(Self::from_inner)
            .map_err(|e| crate::canonical_inchi::error(py, e))
    }
    /// Parse InChI text into a new Molecule. sanitize and remove_hs default to True; params and keyword options are mutually exclusive. Invalid input raises InchiError. Uses the supplied configuration object.
    #[staticmethod]
    fn from_inchi_with_params(
        py: Python<'_>,
        text: crate::text_input::TextInput<'_>,
        params: &crate::canonical_inchi::InchiReadParams,
    ) -> PyResult<Self> {
        ck::Molecule::from_inchi_with_params(&text.as_text()?, &params.inner)
            .map(Self::from_inner)
            .map_err(|e| crate::canonical_inchi::error(py, e))
    }
    /// Return InChI text for this molecule without changing it. Supply InchiWriteParams or options, not both; failures raise InchiError.
    #[pyo3(signature=(params=None, *, options=None))]
    fn to_inchi(
        &self,
        py: Python<'_>,
        params: Option<&crate::canonical_inchi::InchiWriteParams>,
        options: Option<String>,
    ) -> PyResult<String> {
        self.inner
            .to_inchi_with_params(&crate::canonical_inchi::write_params(params, options)?)
            .map_err(|e| crate::canonical_inchi::error(py, e))
    }
    /// Return InChI text for this molecule without changing it. Supply InchiWriteParams or options, not both; failures raise InchiError. Uses the supplied configuration object.
    fn to_inchi_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_inchi::InchiWriteParams,
    ) -> PyResult<String> {
        self.inner
            .to_inchi_with_params(&params.inner)
            .map_err(|e| crate::canonical_inchi::error(py, e))
    }
    /// Return the InChIKey for this molecule using the selected InChI options; leave the molecule unchanged.
    #[pyo3(signature=(params=None, *, options=None))]
    fn to_inchi_key(
        &self,
        py: Python<'_>,
        params: Option<&crate::canonical_inchi::InchiWriteParams>,
        options: Option<String>,
    ) -> PyResult<String> {
        self.inner
            .to_inchi_key_with_params(&crate::canonical_inchi::write_params(params, options)?)
            .map_err(|e| crate::canonical_inchi::error(py, e))
    }
    /// Return the InChIKey for this molecule using the selected InChI options; leave the molecule unchanged. Uses the supplied configuration object.
    fn to_inchi_key_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_inchi::InchiWriteParams,
    ) -> PyResult<String> {
        self.inner
            .to_inchi_key_with_params(&params.inner)
            .map_err(|e| crate::canonical_inchi::error(py, e))
    }
    /// Run this molecule against the selected reactant template and return reaction product sets; the source molecule is unchanged.
    fn reaction_products(
        &self,
        py: Python<'_>,
        mut reaction: PyRefMut<'_, crate::canonical_reaction::Reaction>,
        reactant_template: usize,
    ) -> PyResult<Vec<Vec<Molecule>>> {
        self.inner
            .reaction_products(&mut reaction.inner, reactant_template)
            .map(crate::canonical_reaction::product_sets)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Run this molecule against the selected reactant template and return reaction product sets; the source molecule is unchanged. Uses the supplied configuration object.
    fn reaction_products_with_params(
        &self,
        py: Python<'_>,
        mut reaction: PyRefMut<'_, crate::canonical_reaction::Reaction>,
        reactant_template: usize,
        params: &crate::canonical_reaction::ReactionSingleRunParams,
    ) -> PyResult<Vec<Vec<Molecule>>> {
        self.inner
            .reaction_products_with_params(&mut reaction.inner, reactant_template, &params.inner)
            .map(crate::canonical_reaction::product_sets)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Run the reaction on the explicit reactant list and return product sets in enumeration order.
    fn reaction_products_from_inputs(
        &self,
        py: Python<'_>,
        mut reaction: PyRefMut<'_, crate::canonical_reaction::Reaction>,
        reactants: Vec<PyRef<'_, Molecule>>,
        params: &crate::canonical_reaction::ReactionRunParams,
    ) -> PyResult<Vec<Vec<Molecule>>> {
        let inputs: Vec<_> = reactants.iter().map(|m| &m.inner).collect();
        self.inner
            .reaction_products_from_inputs(&mut reaction.inner, &inputs, &params.inner)
            .map(crate::canonical_reaction::product_sets)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Return the product and changed flag from single-reactant in-place-compatible reaction application, without changing this source molecule.
    fn apply_reaction(
        &self,
        py: Python<'_>,
        mut reaction: PyRefMut<'_, crate::canonical_reaction::Reaction>,
    ) -> PyResult<crate::canonical_reaction::ReactionApplyResult> {
        self.inner
            .apply_reaction(&mut reaction.inner)
            .map(crate::canonical_reaction::apply_result)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Return the product and changed flag from single-reactant in-place-compatible reaction application, without changing this source molecule. Uses the supplied configuration object.
    fn apply_reaction_with_params(
        &self,
        py: Python<'_>,
        mut reaction: PyRefMut<'_, crate::canonical_reaction::Reaction>,
        params: &crate::canonical_reaction::ReactionApplyParams,
    ) -> PyResult<crate::canonical_reaction::ReactionApplyResult> {
        self.inner
            .apply_reaction_with_params(&mut reaction.inner, &params.inner)
            .map(crate::canonical_reaction::apply_result)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply the supported single-reactant reaction in place and return whether the molecule changed.
    fn apply_reaction_(
        &mut self,
        py: Python<'_>,
        mut reaction: PyRefMut<'_, crate::canonical_reaction::Reaction>,
    ) -> PyResult<bool> {
        self.inner
            .apply_reaction_(&mut reaction.inner)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Return the product and changed flag from single-reactant in-place-compatible reaction application, without changing this source molecule. Uses the supplied configuration object.
    fn apply_reaction_with_params_(
        &mut self,
        py: Python<'_>,
        mut reaction: PyRefMut<'_, crate::canonical_reaction::Reaction>,
        params: &crate::canonical_reaction::ReactionApplyParams,
    ) -> PyResult<bool> {
        self.inner
            .apply_reaction_with_params_(&mut reaction.inner, &params.inner)
            .map_err(|e| operation_pyerr(py, e))
    }

    /// Read typed user metadata; missing keys return None.
    #[gen_stub(override_return_type(
        type_repr = "bool | int | float | str | list[int] | list[str] | None"
    ))]
    fn atom_property(&self, py: Python<'_>, atom: usize, key: &str) -> PyResult<Option<Py<PyAny>>> {
        self.inner
            .atom_property(ck::AtomId::new(atom), key)
            .map_err(|e| operation_pyerr(py, e))?
            .map(|v| crate::native_property::NativeProperty::to_python(v, py))
            .transpose()
    }
    /// Return an independently editable molecule with one typed user property.
    fn with_atom_property(
        &self,
        py: Python<'_>,
        atom: usize,
        key: &str,
        value: crate::native_property::NativeProperty,
    ) -> PyResult<Self> {
        self.inner
            .with_atom_property(ck::AtomId::new(atom), key, &value.0)
            .map(Self::from_inner)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Set one typed user property in place; shared copies remain unchanged.
    fn set_atom_property_(
        &mut self,
        py: Python<'_>,
        atom: usize,
        key: &str,
        value: crate::native_property::NativeProperty,
    ) -> PyResult<()> {
        self.inner
            .set_atom_property_(ck::AtomId::new(atom), key, &value.0)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Return the string projection of an atom property; invalid indices/types raise an error.
    fn atom_property_string(
        &self,
        py: Python<'_>,
        id: usize,
        key: &str,
    ) -> PyResult<Option<String>> {
        self.inner
            .atom_property_string(ck::AtomId::new(id), key)
            .map_err(|source| crate::canonical_property_values::property_string_pyerr(py, source))?
            .map(|text| crate::canonical_sdf::decode_source_text(py, &text))
            .transpose()
    }

    /// Return the string projection of a bond property; invalid indices/types raise an error.
    fn bond_property_string(
        &self,
        py: Python<'_>,
        id: usize,
        key: &str,
    ) -> PyResult<Option<String>> {
        self.inner
            .bond_property_string(ck::BondId::new(id), key)
            .map_err(|source| crate::canonical_property_values::property_string_pyerr(py, source))?
            .map(|text| crate::canonical_sdf::decode_source_text(py, &text))
            .transpose()
    }
    /// Return molecules for the allowed stereoisomers using StereoisomerOptions; leave the source unchanged.
    fn enumerate_stereoisomers(
        &self,
        py: Python<'_>,
    ) -> PyResult<crate::canonical_stereoisomers::StereoisomerIterator> {
        crate::canonical_stereoisomers::enumerate_stereoisomers(self, py)
    }
    /// Return molecules for the allowed stereoisomers using StereoisomerOptions; leave the source unchanged. Uses the supplied configuration object.
    fn enumerate_stereoisomers_with_options(
        &self,
        py: Python<'_>,
        options: &crate::canonical_stereoisomers::StereoisomerOptions,
    ) -> PyResult<crate::canonical_stereoisomers::StereoisomerIterator> {
        crate::canonical_stereoisomers::enumerate_stereoisomers_with_options(self, py, options)
    }
    /// Enumerate stereoisomers using the supplied random-bit callback for sampling.
    fn enumerate_stereoisomers_with_random_bits(
        &self,
        py: Python<'_>,
        options: &crate::canonical_stereoisomers::StereoisomerOptions,
        callback: Py<PyAny>,
    ) -> PyResult<crate::canonical_stereoisomers::StereoisomerIterator> {
        crate::canonical_stereoisomers::enumerate_stereoisomers_with_random_bits(
            self, py, options, callback,
        )
    }
    /// Return the estimated number of stereoisomers for the configured enumeration scope.
    #[gen_stub(override_return_type(type_repr = "builtins.int", imports = ("builtins")))]
    fn stereoisomer_count(&self, py: Python<'_>) -> PyResult<num_bigint::BigUint> {
        crate::canonical_stereoisomers::stereoisomer_count(self, py)
    }
    /// Return the estimated number of stereoisomers for the configured enumeration scope. Uses the supplied configuration object.
    #[gen_stub(override_return_type(type_repr = "builtins.int", imports = ("builtins")))]
    fn stereoisomer_count_with_options(
        &self,
        py: Python<'_>,
        options: &crate::canonical_stereoisomers::StereoisomerOptions,
    ) -> PyResult<num_bigint::BigUint> {
        crate::canonical_stereoisomers::stereoisomer_count_with_options(self, py, options)
    }

    /// Parse XYZ without inferring bonds; preserve raw source coordinate values.
    #[staticmethod]
    fn from_xyz_block(py: Python<'_>, text: &str) -> PyResult<Self> {
        ck::Molecule::from_xyz_block(text)
            .map(Self::from_inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Read a XYZ file from a filesystem path and return a new Molecule; use from_* for in-memory text.
    #[staticmethod]
    fn read_xyz(py: Python<'_>, path: &str) -> PyResult<Self> {
        ck::Molecule::read_xyz(path)
            .map(Self::from_inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Return XYZ text for the selected stored 3D conformer; missing or ambiguous coordinates raise an error.
    fn to_xyz(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .to_xyz()
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Return XYZ text for the selected stored 3D conformer; missing or ambiguous coordinates raise an error. Uses the supplied configuration object.
    fn to_xyz_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_molecular_io::XyzWriteParams,
    ) -> PyResult<String> {
        self.inner
            .to_xyz_with_params(&params.inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Write XYZ output to the supplied filesystem path/directory using the selected options. This writes files rather than returning serialized text.
    fn write_xyz(&self, py: Python<'_>, path: &str) -> PyResult<()> {
        self.inner
            .write_xyz(path)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Write XYZ output to the supplied filesystem path/directory using the selected options. This writes files rather than returning serialized text.
    fn write_xyz_with_params(
        &self,
        py: Python<'_>,
        path: &str,
        params: &crate::canonical_molecular_io::XyzWriteParams,
    ) -> PyResult<()> {
        self.inner
            .write_xyz_with_params(path, &params.inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }

    /// Parse a MOL text block into a Molecule with the selected coordinate policy.
    #[staticmethod]
    fn from_mol(py: Python<'_>, text: &str) -> PyResult<Self> {
        ck::Molecule::from_mol(text)
            .map(Self::from_inner)
            .map_err(|e| crate::canonical_sdf::sdf_pyerr(py, e))
    }
    /// Parse a MOL text block into a Molecule with the selected coordinate policy. Uses the supplied configuration object.
    #[staticmethod]
    fn from_mol_with_params(
        py: Python<'_>,
        text: &str,
        params: &crate::canonical_sdf::SdfReadParams,
    ) -> PyResult<Self> {
        ck::Molecule::from_mol_with_params(text, &params.inner)
            .map(Self::from_inner)
            .map_err(|e| crate::canonical_sdf::sdf_pyerr(py, e))
    }
    /// Parse an SDF text record into a Molecule with the selected coordinate policy; use SdfRecord when the record may contain a query graph.
    #[staticmethod]
    fn from_sdf(py: Python<'_>, text: &str) -> PyResult<Self> {
        ck::Molecule::from_sdf(text)
            .map(Self::from_inner)
            .map_err(|e| crate::canonical_sdf::sdf_pyerr(py, e))
    }
    /// Parse an SDF text record into a Molecule with the selected coordinate policy; use SdfRecord when the record may contain a query graph. Uses the supplied configuration object.
    #[staticmethod]
    fn from_sdf_with_params(
        py: Python<'_>,
        text: &str,
        params: &crate::canonical_sdf::SdfReadParams,
    ) -> PyResult<Self> {
        ck::Molecule::from_sdf_with_params(text, &params.inner)
            .map(Self::from_inner)
            .map_err(|e| crate::canonical_sdf::sdf_pyerr(py, e))
    }
    /// Read a MOL file from a filesystem path and return a new Molecule; use from_* for in-memory text.
    #[staticmethod]
    fn read_mol(py: Python<'_>, path: &str) -> PyResult<Self> {
        ck::Molecule::read_mol(path)
            .map(Self::from_inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Read a MOL file from a filesystem path and return a new Molecule; use from_* for in-memory text.
    #[staticmethod]
    fn read_mol_with_params(
        py: Python<'_>,
        path: &str,
        params: &crate::canonical_sdf::SdfReadParams,
    ) -> PyResult<Self> {
        ck::Molecule::read_mol_with_params(path, &params.inner)
            .map(Self::from_inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Read a SDF file from a filesystem path and return a new Molecule; use from_* for in-memory text.
    #[staticmethod]
    fn read_sdf(py: Python<'_>, path: &str) -> PyResult<Self> {
        ck::Molecule::read_sdf(path)
            .map(Self::from_inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Read a SDF file from a filesystem path and return a new Molecule; use from_* for in-memory text.
    #[staticmethod]
    fn read_sdf_with_params(
        py: Python<'_>,
        path: &str,
        params: &crate::canonical_sdf::SdfReadParams,
    ) -> PyResult<Self> {
        ck::Molecule::read_sdf_with_params(path, &params.inner)
            .map(Self::from_inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Parse a Tripos MOL2 text block into a new Molecule using the read options.
    #[staticmethod]
    fn from_mol2(py: Python<'_>, text: &str) -> PyResult<Self> {
        ck::Molecule::from_mol2(text)
            .map(Self::from_inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Parse a Tripos MOL2 text block into a new Molecule using the read options. Uses the supplied configuration object.
    #[staticmethod]
    fn from_mol2_with_params(
        py: Python<'_>,
        text: &str,
        params: &crate::canonical_molecular_io::Mol2ReadParams,
    ) -> PyResult<Self> {
        ck::Molecule::from_mol2_with_params(text, &params.inner)
            .map(Self::from_inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Read a MOL2 file from a filesystem path and return a new Molecule; use from_* for in-memory text.
    #[staticmethod]
    fn read_mol2(py: Python<'_>, path: &str) -> PyResult<Self> {
        ck::Molecule::read_mol2(path)
            .map(Self::from_inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Read a MOL2 file from a filesystem path and return a new Molecule; use from_* for in-memory text.
    #[staticmethod]
    fn read_mol2_with_params(
        py: Python<'_>,
        path: &str,
        params: &crate::canonical_molecular_io::Mol2ReadParams,
    ) -> PyResult<Self> {
        ck::Molecule::read_mol2_with_params(path, &params.inner)
            .map(Self::from_inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Return a MOL text block using the selected writer options; does not write a file or install generated drawing coordinates on this object.
    fn to_mol(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .to_mol()
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Return a MOL text block using the selected writer options; does not write a file or install generated drawing coordinates on this object. Uses the supplied configuration object.
    fn to_mol_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_molecular_io::MolBlockWriteParams,
    ) -> PyResult<String> {
        self.inner
            .to_mol_with_params(&params.inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Return an SDF text record using the selected writer options; does not write a file.
    fn to_sdf(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .to_sdf()
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Return an SDF text record using the selected writer options; does not write a file. Uses the supplied configuration object.
    fn to_sdf_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_molecular_io::MolBlockWriteParams,
    ) -> PyResult<String> {
        self.inner
            .to_sdf_with_params(&params.inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Return SDF text with 2D coordinates. When needed and enabled, generate temporary 2D coordinates for export without installing them on the source molecule.
    fn to_sdf_2d(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .to_sdf_2d()
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Return SDF text with 2D coordinates. When needed and enabled, generate temporary 2D coordinates for export without installing them on the source molecule. Uses the supplied configuration object.
    fn to_sdf_2d_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_molecular_io::MolBlockWriteParams,
    ) -> PyResult<String> {
        self.inner
            .to_sdf_2d_with_params(&params.inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Return SDF text for the selected stored 3D conformer; does not fall back to 2D coordinates.
    fn to_sdf_3d(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .to_sdf_3d()
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Return SDF text for the selected stored 3D conformer; does not fall back to 2D coordinates. Uses the supplied configuration object.
    fn to_sdf_3d_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_molecular_io::MolBlockWriteParams,
    ) -> PyResult<String> {
        self.inner
            .to_sdf_3d_with_params(&params.inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Write MOL output to the supplied filesystem path/directory using the selected options. This writes files rather than returning serialized text.
    fn write_mol(&self, py: Python<'_>, path: &str) -> PyResult<()> {
        self.inner
            .write_mol(path)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Write MOL output to the supplied filesystem path/directory using the selected options. This writes files rather than returning serialized text.
    fn write_mol_with_params(
        &self,
        py: Python<'_>,
        path: &str,
        params: &crate::canonical_molecular_io::MolBlockWriteParams,
    ) -> PyResult<()> {
        self.inner
            .write_mol_with_params(path, &params.inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }

    /// Write SDF output to the supplied filesystem path/directory using the selected options. This writes files rather than returning serialized text.
    fn write_sdf(&self, py: Python<'_>, path: &str) -> PyResult<()> {
        self.inner
            .write_sdf(path)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Write SDF output to the supplied filesystem path/directory using the selected options. This writes files rather than returning serialized text.
    fn write_sdf_with_params(
        &self,
        py: Python<'_>,
        path: &str,
        params: &crate::canonical_molecular_io::MolBlockWriteParams,
    ) -> PyResult<()> {
        self.inner
            .write_sdf_with_params(path, &params.inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Write SDF FILES output to the supplied filesystem path/directory using the selected options. This writes files rather than returning serialized text.
    #[pyo3(signature=(directory, file_name=None))]
    fn write_sdf_files(
        &self,
        py: Python<'_>,
        directory: &str,
        file_name: Option<&str>,
    ) -> PyResult<std::path::PathBuf> {
        self.inner
            .write_sdf_files(directory, file_name)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Write SDF FILES output to the supplied filesystem path/directory using the selected options. This writes files rather than returning serialized text.
    fn write_sdf_files_with_params(
        &self,
        py: Python<'_>,
        directory: &str,
        file_name: Option<&str>,
        params: &crate::canonical_molecular_io::MolBlockWriteParams,
    ) -> PyResult<std::path::PathBuf> {
        self.inner
            .write_sdf_files_with_params(directory, file_name, &params.inner)
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Compute pattern fixed-width bit fingerprints for this molecule without changing its graph or coordinates.
    fn fingerprint_pattern(
        &self,
        py: Python<'_>,
    ) -> PyResult<crate::canonical_values::Fingerprint> {
        self.inner
            .fingerprint_pattern()
            .map(|inner| crate::canonical_values::Fingerprint { inner })
            .map_err(|e| crate::canonical_pattern::pattern_pyerr(py, e))
    }
    /// Compute pattern fixed-width bit fingerprints for this molecule without changing its graph or coordinates.
    fn fingerprint_pattern_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_pattern::PatternFingerprintParams,
    ) -> PyResult<crate::canonical_values::Fingerprint> {
        self.inner
            .fingerprint_pattern_with_params(&params.inner)
            .map(|inner| crate::canonical_values::Fingerprint { inner })
            .map_err(|e| crate::canonical_pattern::pattern_pyerr(py, e))
    }
    /// Compute topological fixed-width bit fingerprints for this molecule without changing its graph or coordinates.
    fn fingerprint_topological(
        &self,
        py: Python<'_>,
    ) -> PyResult<crate::canonical_values::Fingerprint> {
        self.inner
            .fingerprint_topological()
            .map(|inner| crate::canonical_values::Fingerprint { inner })
            .map_err(|e| crate::canonical_topological::topological_pyerr(py, e))
    }
    /// Compute topological fixed-width bit fingerprints for this molecule without changing its graph or coordinates.
    fn fingerprint_topological_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_topological::TopologicalFingerprintParams,
    ) -> PyResult<crate::canonical_values::Fingerprint> {
        self.inner
            .fingerprint_topological_with_params(&params.inner)
            .map(|inner| crate::canonical_values::Fingerprint { inner })
            .map_err(|e| crate::canonical_topological::topological_pyerr(py, e))
    }
    /// Compute topological fixed-width bit fingerprints for this molecule without changing its graph or coordinates; include the requested atom/bit environment metadata.
    fn fingerprint_topological_with_output(
        &self,
        py: Python<'_>,
    ) -> PyResult<crate::canonical_topological::TopologicalFingerprintResult> {
        self.inner
            .fingerprint_topological_with_output()
            .map(|inner| crate::canonical_topological::TopologicalFingerprintResult { inner })
            .map_err(|e| crate::canonical_topological::topological_pyerr(py, e))
    }
    /// Compute topological fixed-width bit fingerprints for this molecule without changing its graph or coordinates; include the requested atom/bit environment metadata.
    fn fingerprint_topological_with_output_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_topological::TopologicalFingerprintParams,
        request: &crate::canonical_topological::TopologicalFingerprintOutputRequest,
    ) -> PyResult<crate::canonical_topological::TopologicalFingerprintResult> {
        self.inner
            .fingerprint_topological_with_output_with_params(&params.inner, request.inner)
            .map(|inner| crate::canonical_topological::TopologicalFingerprintResult { inner })
            .map_err(|e| crate::canonical_topological::topological_pyerr(py, e))
    }
    /// Compute maccs fixed-width bit fingerprints for this molecule without changing its graph or coordinates.
    fn fingerprint_maccs(&self, py: Python<'_>) -> PyResult<crate::canonical_values::Fingerprint> {
        self.inner
            .fingerprint_maccs()
            .map(|inner| crate::canonical_values::Fingerprint { inner })
            .map_err(|e| crate::canonical_maccs::maccs_pyerr(py, e))
    }
    /// Compute maccs fixed-width bit fingerprints for this molecule without changing its graph or coordinates.
    fn fingerprint_maccs_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_maccs::MaccsFingerprintParams,
    ) -> PyResult<crate::canonical_values::Fingerprint> {
        self.inner
            .fingerprint_maccs_with_params(&params.inner)
            .map(|inner| crate::canonical_values::Fingerprint { inner })
            .map_err(|e| crate::canonical_maccs::maccs_pyerr(py, e))
    }
    /// Compute maccs fixed-width bit fingerprints for this molecule without changing its graph or coordinates.
    fn fingerprint_maccs_raw(
        &self,
        py: Python<'_>,
    ) -> PyResult<crate::canonical_values::Fingerprint> {
        self.inner
            .fingerprint_maccs_raw()
            .map(|inner| crate::canonical_values::Fingerprint { inner })
            .map_err(|e| crate::canonical_maccs::maccs_pyerr(py, e))
    }
    /// Compute layered fixed-width bit fingerprints for this molecule without changing its graph or coordinates.
    fn fingerprint_layered(
        &self,
        py: Python<'_>,
    ) -> PyResult<crate::canonical_values::Fingerprint> {
        self.inner
            .fingerprint_layered()
            .map(|inner| crate::canonical_values::Fingerprint { inner })
            .map_err(|e| crate::canonical_layered::layered_pyerr(py, e))
    }
    /// RDKit Python defaults; accepts configuration objects or keyword options.
    fn fingerprint_avalon(&self, py: Python<'_>) -> PyResult<crate::canonical_values::Fingerprint> {
        // RDKit External/AvalonTools/Wrap/pyAvalonTools.cpp:
        // python::arg("bitFlags") = AvalonTools::avalonSimilarityBits
        let params = ck::AvalonFingerprintParams {
            bit_flags: ck::AvalonFingerprintFlags::SIMILARITY,
            ..Default::default()
        };
        self.inner
            .fingerprint_avalon_with_params(&params)
            .map(|inner| crate::canonical_values::Fingerprint { inner })
            .map_err(|e| crate::canonical_avalon::error(py, e))
    }
    /// Compute avalon fixed-width bit fingerprints for this molecule without changing its graph or coordinates.
    fn fingerprint_avalon_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_avalon::AvalonFingerprintParams,
    ) -> PyResult<crate::canonical_values::Fingerprint> {
        self.inner
            .fingerprint_avalon_with_params(&params.inner)
            .map(|inner| crate::canonical_values::Fingerprint { inner })
            .map_err(|e| crate::canonical_avalon::error(py, e))
    }
    /// Compute layered fixed-width bit fingerprints for this molecule without changing its graph or coordinates.
    fn fingerprint_layered_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_layered::LayeredFingerprintParams,
    ) -> PyResult<crate::canonical_values::Fingerprint> {
        self.inner
            .fingerprint_layered_with_params(&params.inner)
            .map(|inner| crate::canonical_values::Fingerprint { inner })
            .map_err(|e| crate::canonical_layered::layered_pyerr(py, e))
    }
    /// Compute layered fixed-width bit fingerprints for this molecule without changing its graph or coordinates; include the requested atom/bit environment metadata.
    fn fingerprint_layered_with_output(
        &self,
        py: Python<'_>,
    ) -> PyResult<crate::canonical_layered::LayeredFingerprintResult> {
        self.inner
            .fingerprint_layered_with_output()
            .map(|inner| crate::canonical_layered::LayeredFingerprintResult { inner })
            .map_err(|e| crate::canonical_layered::layered_pyerr(py, e))
    }
    /// Compute layered fixed-width bit fingerprints for this molecule without changing its graph or coordinates; include the requested atom/bit environment metadata.
    fn fingerprint_layered_with_output_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_layered::LayeredFingerprintParams,
    ) -> PyResult<crate::canonical_layered::LayeredFingerprintResult> {
        self.inner
            .fingerprint_layered_with_output_with_params(&params.inner)
            .map(|inner| crate::canonical_layered::LayeredFingerprintResult { inner })
            .map_err(|e| crate::canonical_layered::layered_pyerr(py, e))
    }
    /// Return the packed topological-torsion identifier for the supplied path and atom codes.
    #[pyo3(signature=(path,size,atom_codes=None))]
    fn topological_torsion_path_score(
        &self,
        py: Python<'_>,
        path: Vec<usize>,
        size: usize,
        atom_codes: Option<Vec<u32>>,
    ) -> PyResult<u64> {
        self.inner
            .topological_torsion_path_score(&path, size, atom_codes.as_deref())
            .map_err(|e| crate::canonical_path_score::score_pyerr(py, e))
    }
    /// Return a detached snapshot of the stored properties.
    fn properties(&self) -> crate::canonical_property_values::MoleculeProperties {
        crate::canonical_property_values::MoleculeProperties {
            inner: self.inner.properties().clone(),
        }
    }
    /// Enumerate tautomers and return a TautomerEnumeration containing molecules, status and modified atom/bond indices.
    fn enumerate_tautomers(
        &self,
        py: Python<'_>,
    ) -> PyResult<crate::tautomer_binding::TautomerEnumeration> {
        crate::tautomer_binding::enumerate(py, self, None)
    }
    /// Enumerate tautomers and return a TautomerEnumeration containing molecules, status and modified atom/bond indices. Uses the supplied configuration object.
    fn enumerate_tautomers_with_params(
        &self,
        py: Python<'_>,
        params: &crate::tautomer_binding::TautomerParams,
    ) -> PyResult<crate::tautomer_binding::TautomerEnumeration> {
        crate::tautomer_binding::enumerate(py, self, Some(params))
    }
    /// Return the highest-ranked canonical tautomer using the configured scoring rules; leave the source unchanged.
    fn canonical_tautomer(&self, py: Python<'_>) -> PyResult<Self> {
        crate::tautomer_binding::canonical(py, self, None)
    }
    /// Return the highest-ranked canonical tautomer using the configured scoring rules; leave the source unchanged. Uses the supplied configuration object.
    fn canonical_tautomer_with_params(
        &self,
        py: Python<'_>,
        params: &crate::tautomer_binding::TautomerParams,
    ) -> PyResult<Self> {
        crate::tautomer_binding::canonical(py, self, Some(params))
    }
    /// Return ring, SMARTS-pattern and heteroatom-hydrogen score contributions for this tautomer.
    fn tautomer_score(&self, py: Python<'_>) -> PyResult<crate::tautomer_binding::TautomerScore> {
        self.inner
            .tautomer_score()
            .map(|inner| crate::tautomer_binding::TautomerScore { inner })
            .map_err(|error| operation_pyerr(py, error))
    }
    /// Return ring, SMARTS-pattern and heteroatom-hydrogen score contributions for this tautomer. Uses the supplied configuration object.
    fn tautomer_score_with_params(
        &self,
        py: Python<'_>,
        params: &crate::tautomer_binding::TautomerScoreParams,
    ) -> PyResult<crate::tautomer_binding::TautomerScore> {
        self.inner
            .tautomer_score_with_params(&params.inner)
            .map(|inner| crate::tautomer_binding::TautomerScore { inner })
            .map_err(|error| operation_pyerr(py, error))
    }
    /// Return the atom-pair atom code and associated molecule result for the requested atom.
    #[pyo3(signature=(atom_id, branch_subtract=0, include_chirality=false, use_legacy_stereo_perception=true))]
    fn with_atom_pair_atom_code(
        &self,
        py: Python<'_>,
        atom_id: usize,
        branch_subtract: u32,
        include_chirality: bool,
        use_legacy_stereo_perception: bool,
    ) -> PyResult<crate::canonical_fingerprint_values::AtomPairAtomCodeResult> {
        self.inner
            .with_atom_pair_atom_code(
                ck::AtomId::new(atom_id),
                branch_subtract,
                include_chirality,
                use_legacy_stereo_perception,
            )
            .map(|inner| crate::canonical_fingerprint_values::AtomPairAtomCodeResult { inner })
            .map_err(|error| operation_pyerr(py, error))
    }
    /// Compute morgan fixed-width bit fingerprints for this molecule without changing its graph or coordinates using the supplied generator and per-call options.
    #[pyo3(signature=(generator,*,params=None,output=None))]
    fn fingerprint_morgan_with_generator(
        &self,
        py: Python<'_>,
        generator: &crate::canonical_fingerprint_values::MorganFingerprintGenerator,
        params: Option<&crate::canonical_fingerprint_values::MorganCallParams>,
        mut output: Option<
            PyRefMut<'_, crate::canonical_fingerprint_values::FingerprintAdditionalOutput>,
        >,
    ) -> PyResult<crate::canonical_values::Fingerprint> {
        self.inner
            .fingerprint_morgan_with_generator(
                &generator.inner,
                params.map(|p| &p.inner),
                output.as_deref_mut().map(|o| &mut o.inner),
            )
            .map(|inner| crate::canonical_values::Fingerprint { inner })
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    /// Compute morgan folded feature counts for this molecule without changing its graph or coordinates using the supplied generator and per-call options.
    #[pyo3(signature=(generator,*,params=None,output=None))]
    fn fingerprint_morgan_count_with_generator(
        &self,
        py: Python<'_>,
        generator: &crate::canonical_fingerprint_values::MorganFingerprintGenerator,
        params: Option<&crate::canonical_fingerprint_values::MorganCallParams>,
        mut output: Option<
            PyRefMut<'_, crate::canonical_fingerprint_values::FingerprintAdditionalOutput>,
        >,
    ) -> PyResult<crate::canonical_values::SparseCountFingerprint32> {
        self.inner
            .fingerprint_morgan_count_with_generator(
                &generator.inner,
                params.map(|p| &p.inner),
                output.as_deref_mut().map(|o| &mut o.inner),
            )
            .map(|inner| crate::canonical_values::SparseCountFingerprint32 { inner })
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    /// Compute morgan sparse feature bits for this molecule without changing its graph or coordinates using the supplied generator and per-call options.
    #[pyo3(signature=(generator,*,params=None,output=None))]
    fn fingerprint_morgan_sparse_with_generator(
        &self,
        py: Python<'_>,
        generator: &crate::canonical_fingerprint_values::MorganFingerprintGenerator,
        params: Option<&crate::canonical_fingerprint_values::MorganCallParams>,
        mut output: Option<
            PyRefMut<'_, crate::canonical_fingerprint_values::FingerprintAdditionalOutput>,
        >,
    ) -> PyResult<crate::canonical_values::SparseBitFingerprint> {
        self.inner
            .fingerprint_morgan_sparse_with_generator(
                &generator.inner,
                params.map(|p| &p.inner),
                output.as_deref_mut().map(|o| &mut o.inner),
            )
            .map(|inner| crate::canonical_values::SparseBitFingerprint { inner })
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    /// Compute morgan sparse feature counts for this molecule without changing its graph or coordinates using the supplied generator and per-call options.
    #[pyo3(signature=(generator,*,params=None,output=None))]
    fn fingerprint_morgan_sparse_count_with_generator(
        &self,
        py: Python<'_>,
        generator: &crate::canonical_fingerprint_values::MorganFingerprintGenerator,
        params: Option<&crate::canonical_fingerprint_values::MorganCallParams>,
        mut output: Option<
            PyRefMut<'_, crate::canonical_fingerprint_values::FingerprintAdditionalOutput>,
        >,
    ) -> PyResult<crate::canonical_values::SparseCountFingerprint> {
        self.inner
            .fingerprint_morgan_sparse_count_with_generator(
                &generator.inner,
                params.map(|p| &p.inner),
                output.as_deref_mut().map(|o| &mut o.inner),
            )
            .map(|inner| crate::canonical_values::SparseCountFingerprint { inner })
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    /// Compute topological torsion fixed-width bit fingerprints for this molecule without changing its graph or coordinates using the supplied generator and per-call options.
    #[pyo3(signature=(generator,*,params=None,output=None))]
    fn fingerprint_topological_torsion_with_generator(
        &self,
        py: Python<'_>,
        generator: &TopologicalTorsionFingerprintGenerator,
        params: Option<&TopologicalTorsionCallParams>,
        mut output: Option<PyRefMut<'_, FingerprintAdditionalOutput>>,
    ) -> PyResult<Fingerprint> {
        self.inner
            .fingerprint_topological_torsion_with_generator(
                &generator.inner,
                params.map(|p| &p.inner),
                output.as_deref_mut().map(|o| &mut o.inner),
            )
            .map(|inner| Fingerprint { inner })
            .map_err(|e| topological_torsion_pyerr(py, e))
    }

    /// Compute topological torsion sparse feature bits for this molecule without changing its graph or coordinates using the supplied generator and per-call options.
    #[pyo3(signature=(generator,*,params=None,output=None))]
    fn fingerprint_topological_torsion_sparse_with_generator(
        &self,
        py: Python<'_>,
        generator: &TopologicalTorsionFingerprintGenerator,
        params: Option<&TopologicalTorsionCallParams>,
        mut output: Option<PyRefMut<'_, FingerprintAdditionalOutput>>,
    ) -> PyResult<SparseBitFingerprint> {
        self.inner
            .fingerprint_topological_torsion_sparse_with_generator(
                &generator.inner,
                params.map(|p| &p.inner),
                output.as_deref_mut().map(|o| &mut o.inner),
            )
            .map(|inner| SparseBitFingerprint { inner })
            .map_err(|e| topological_torsion_pyerr(py, e))
    }

    /// Compute topological torsion folded feature counts for this molecule without changing its graph or coordinates using the supplied generator and per-call options.
    #[pyo3(signature=(generator,*,params=None,output=None))]
    fn fingerprint_topological_torsion_count_with_generator(
        &self,
        py: Python<'_>,
        generator: &TopologicalTorsionFingerprintGenerator,
        params: Option<&TopologicalTorsionCallParams>,
        mut output: Option<PyRefMut<'_, FingerprintAdditionalOutput>>,
    ) -> PyResult<SparseCountFingerprint32> {
        self.inner
            .fingerprint_topological_torsion_count_with_generator(
                &generator.inner,
                params.map(|p| &p.inner),
                output.as_deref_mut().map(|o| &mut o.inner),
            )
            .map(|inner| SparseCountFingerprint32 { inner })
            .map_err(|e| topological_torsion_pyerr(py, e))
    }

    /// Compute topological torsion sparse feature counts for this molecule without changing its graph or coordinates using the supplied generator and per-call options.
    #[pyo3(signature=(generator,*,params=None,output=None))]
    fn fingerprint_topological_torsion_sparse_count_with_generator(
        &self,
        py: Python<'_>,
        generator: &TopologicalTorsionFingerprintGenerator,
        params: Option<&TopologicalTorsionCallParams>,
        mut output: Option<PyRefMut<'_, FingerprintAdditionalOutput>>,
    ) -> PyResult<SparseCountFingerprint> {
        self.inner
            .fingerprint_topological_torsion_sparse_count_with_generator(
                &generator.inner,
                params.map(|p| &p.inner),
                output.as_deref_mut().map(|o| &mut o.inner),
            )
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(|e| topological_torsion_pyerr(py, e))
    }

    /// Compute topological torsion sparse feature counts for this molecule without changing its graph or coordinates using the torsion-vector entry point.
    fn fingerprint_topological_torsion_sparse_count_legacy(
        &self,
        py: Python<'_>,
    ) -> PyResult<SparseCountFingerprint> {
        self.inner
            .fingerprint_topological_torsion_sparse_count_legacy()
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(|error| topological_torsion_pyerr(py, error))
    }

    /// Compute topological torsion sparse feature counts for this molecule without changing its graph or coordinates using the torsion-vector entry point.
    fn fingerprint_topological_torsion_sparse_count_legacy_with_params(
        &self,
        py: Python<'_>,
        params: &LegacyTopologicalTorsionParams,
    ) -> PyResult<SparseCountFingerprint> {
        self.inner
            .fingerprint_topological_torsion_sparse_count_legacy_with_params(&params.inner)
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(|error| topological_torsion_pyerr(py, error))
    }

    /// Compute topological torsion folded feature counts for this molecule without changing its graph or coordinates using the torsion-vector entry point.
    fn fingerprint_topological_torsion_count_legacy(
        &self,
        py: Python<'_>,
    ) -> PyResult<SparseCountFingerprint> {
        self.inner
            .fingerprint_topological_torsion_count_legacy()
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(|error| topological_torsion_pyerr(py, error))
    }

    /// Compute topological torsion folded feature counts for this molecule without changing its graph or coordinates using the torsion-vector entry point.
    fn fingerprint_topological_torsion_count_legacy_with_params(
        &self,
        py: Python<'_>,
        params: &LegacyTopologicalTorsionParams,
    ) -> PyResult<SparseCountFingerprint> {
        self.inner
            .fingerprint_topological_torsion_count_legacy_with_params(&params.inner)
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(|error| topological_torsion_pyerr(py, error))
    }

    /// Compute topological torsion fixed-width bit fingerprints for this molecule without changing its graph or coordinates using the torsion-vector entry point.
    fn fingerprint_topological_torsion_legacy(&self, py: Python<'_>) -> PyResult<Fingerprint> {
        self.inner
            .fingerprint_topological_torsion_legacy()
            .map(|inner| Fingerprint { inner })
            .map_err(|error| topological_torsion_pyerr(py, error))
    }

    /// Compute topological torsion fixed-width bit fingerprints for this molecule without changing its graph or coordinates using the torsion-vector entry point.
    fn fingerprint_topological_torsion_legacy_with_params(
        &self,
        py: Python<'_>,
        params: &LegacyTopologicalTorsionParams,
    ) -> PyResult<Fingerprint> {
        self.inner
            .fingerprint_topological_torsion_legacy_with_params(&params.inner)
            .map(|inner| Fingerprint { inner })
            .map_err(|error| topological_torsion_pyerr(py, error))
    }

    /// Compute the original native 64-bit hash in legacy CIP rank order.
    /// Requires already prepared valence; absent ordinary rings remain absent.
    fn molecular_hash(&self, py: Python<'_>) -> PyResult<u64> {
        self.inner
            .molecular_hash()
            .map_err(|error| crate::canonical_molecular_hash::error_pyerr(py, &error))
    }

    /// Compute the same native hash with exactly one supplied rank per atom.
    fn molecular_hash_with_ranks(&self, py: Python<'_>, ranks: Vec<u32>) -> PyResult<u64> {
        self.inner
            .molecular_hash_with_ranks(&ranks)
            .map_err(|error| crate::canonical_molecular_hash::error_pyerr(py, &error))
    }

    /// Return a native COSMolKit archive as bytes, including stored graph, coordinates, properties and retained derived state.
    fn to_binary(&self, py: Python<'_>) -> PyResult<Py<pyo3::types::PyBytes>> {
        let data = self
            .inner
            .to_binary()
            .map_err(|error| crate::canonical_binary::error_pyerr(py, &error))?;
        Ok(pyo3::types::PyBytes::new(py, &data).unbind())
    }

    /// Read a native COSMolKit archive into a new Molecule, dispatching supported archive versions automatically.
    #[staticmethod]
    fn from_binary(
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "builtins.bytes", imports = ("builtins")))]
        data: &[u8],
    ) -> PyResult<Self> {
        ck::Molecule::from_binary(data)
            .map(|inner| Self { inner })
            .map_err(|error| crate::canonical_binary::error_pyerr(py, &error))
    }

    fn __reduce__(&self, py: Python<'_>) -> PyResult<(Py<PyAny>, (Py<pyo3::types::PyBytes>,))> {
        let rebuild = py
            .import("cosmolkit")?
            .getattr("_molecule_from_binary")?
            .unbind();
        Ok((rebuild, (self.to_binary(py)?,)))
    }

    fn __reduce_ex__(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "typing.SupportsIndex", imports = ("typing")))]
        protocol: &Bound<'_, PyAny>,
    ) -> PyResult<(Py<PyAny>, (Py<pyo3::types::PyBytes>,))> {
        // Follow Python's pickle protocol input, including __index__ objects.
        py.import("operator")?
            .getattr("index")?
            .call1((protocol,))?;
        self.__reduce__(py)
    }

    /// Return the number of amide bonds.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn num_amide_bonds(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_amide_bonds()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the number of atoms shared by otherwise disjoint rings.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn num_spiro_atoms(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_spiro_atoms()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the number of bridgehead atoms in the ring system.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn num_bridgehead_atoms(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_bridgehead_atoms()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the number of specified and unspecified atom stereocenters.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn num_atom_stereo_centers(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_atom_stereo_centers()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the number of atom stereocenters without a specified configuration.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn num_unspecified_atom_stereo_centers(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_unspecified_atom_stereo_centers()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the rotatable bond count using the selected RotatableBondsOptions definition.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn num_rotatable_bonds(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_rotatable_bonds()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the rotatable bond count using the selected RotatableBondsOptions definition. Uses the explicit arguments shown in the signature.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn num_rotatable_bonds_with_params(
        &self,
        py: Python<'_>,
        params: RotatableBondsOptions,
    ) -> PyResult<u32> {
        self.inner
            .num_rotatable_bonds_with_params(&params.canonical())
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the average molecular weight in g/mol, including implicit hydrogens. Uses the explicit arguments shown in the signature.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn molecular_weight_with_params(&self, py: Python<'_>, only_heavy: bool) -> PyResult<f64> {
        self.inner
            .molecular_weight_with_params(only_heavy)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the monoisotopic molecular mass, respecting explicit isotope labels. Uses the explicit arguments shown in the signature.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn exact_molecular_weight_with_params(
        &self,
        py: Python<'_>,
        only_heavy: bool,
    ) -> PyResult<f64> {
        self.inner
            .exact_molecular_weight_with_params(only_heavy)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the molecular formula using the selected isotope and element-count formatting. Uses the explicit arguments shown in the signature.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
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

    /// Return the average molecular weight in g/mol, including implicit hydrogens.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn molecular_weight(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .molecular_weight()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the monoisotopic molecular mass, respecting explicit isotope labels.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn exact_molecular_weight(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .exact_molecular_weight()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the molecular formula using the selected isotope and element-count formatting.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn molecular_formula(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .molecular_formula()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the number of non-hydrogen graph atoms.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn num_heavy_atoms(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_heavy_atoms()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the total atom count including graph atoms and implicit hydrogens.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn total_atom_count(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .total_atom_count()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the Lipinski hydrogen-bond acceptor count.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn lipinski_hba(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .lipinski_hba()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the Lipinski hydrogen-bond donor count.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn lipinski_hbd(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .lipinski_hbd()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the fraction of carbon atoms with sp3 hybridization.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn fraction_csp3(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .fraction_csp3()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the number of graph atoms other than carbon and hydrogen.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn num_heteroatoms(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_heteroatoms()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the hydrogen-bond acceptor count using the descriptor SMARTS definitions.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn num_hba(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_hba()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the hydrogen-bond donor count using the descriptor SMARTS definitions.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn num_hbd(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_hbd()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the number of perceived rings; use the current perceived ring state without modifying the graph.
    fn num_rings(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_rings()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the number of heterocyclic rings; use the current perceived ring state without modifying the graph.
    fn num_heterocycles(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_heterocycles()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the number of aromatic rings; use the current perceived ring state without modifying the graph.
    fn num_aromatic_rings(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_aromatic_rings()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the number of fully saturated rings; use the current perceived ring state without modifying the graph.
    fn num_saturated_rings(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_saturated_rings()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the number of non-aromatic rings; use the current perceived ring state without modifying the graph.
    fn num_aliphatic_rings(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_aliphatic_rings()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the number of aromatic rings containing a heteroatom; use the current perceived ring state without modifying the graph.
    fn num_aromatic_heterocycles(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_aromatic_heterocycles()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the number of all-carbon aromatic rings; use the current perceived ring state without modifying the graph.
    fn num_aromatic_carbocycles(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_aromatic_carbocycles()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the number of non-aromatic rings containing a heteroatom; use the current perceived ring state without modifying the graph.
    fn num_aliphatic_heterocycles(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_aliphatic_heterocycles()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the number of all-carbon non-aromatic rings; use the current perceived ring state without modifying the graph.
    fn num_aliphatic_carbocycles(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_aliphatic_carbocycles()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the number of saturated rings containing a heteroatom; use the current perceived ring state without modifying the graph.
    fn num_saturated_heterocycles(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_saturated_heterocycles()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the number of all-carbon saturated rings; use the current perceived ring state without modifying the graph.
    fn num_saturated_carbocycles(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .num_saturated_carbocycles()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the degree-based zeroth-order Chi molecular connectivity index.
    fn chi_0(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_0()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the degree-based first-order Chi molecular connectivity index.
    fn chi_1(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_1()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the Hall-Kier alpha correction from the stored atomic hybridization state.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn hall_kier_alpha(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .hall_kier_alpha()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the Hall-Kier alpha value and atom-indexed contribution values.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn hall_kier_alpha_with_contributions(&self, py: Python<'_>) -> PyResult<(f64, Vec<f64>)> {
        self.inner
            .hall_kier_alpha_with_contributions()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the first-order Kier molecular shape index.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn kappa_1(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .kappa_1()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the second-order Kier molecular shape index.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn kappa_2(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .kappa_2()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the third-order Kier molecular shape index.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn kappa_3(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .kappa_3()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the Kier molecular flexibility index.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn phi(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .phi()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return all 42 molecular quantum numbers in their defined component order.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    #[pyo3(signature = (force=false))]
    fn mqns(&self, py: Python<'_>, force: bool) -> PyResult<Vec<u32>> {
        self.inner
            .mqns(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the 0th-order valence-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    fn chi_0_v(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_0_v()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the 1th-order valence-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    fn chi_1_v(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_1_v()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the 2th-order valence-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    fn chi_2_v(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_2_v()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the 3th-order valence-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    fn chi_3_v(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_3_v()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the 4th-order valence-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    fn chi_4_v(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_4_v()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the requested-order valence-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    #[pyo3(signature = (order))]
    fn chi_n_v(&self, py: Python<'_>, order: u32) -> PyResult<f64> {
        self.inner
            .chi_n_v(order)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the 0th-order n-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    fn chi_0_n(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_0_n()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the 1th-order n-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    fn chi_1_n(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_1_n()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the 2th-order n-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    fn chi_2_n(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_2_n()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the 3th-order n-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    fn chi_3_n(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_3_n()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the 4th-order n-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    fn chi_4_n(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .chi_4_n()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the requested-order n-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    #[pyo3(signature = (order))]
    fn chi_n_n(&self, py: Python<'_>, order: u32) -> PyResult<f64> {
        self.inner
            .chi_n_n(order)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }
    /// Return read-only atom values in graph order with recalculated valence metadata. Use atom_metadata(recalculate=False) to inspect an existing valid cache without recalculation.
    fn atoms(&self) -> Vec<crate::canonical_atom_bond::Atom> {
        let metadata = self.inner.atom_metadata(true);
        self.inner
            .atoms()
            .iter()
            .enumerate()
            .map(|(index, atom)| crate::canonical_atom_bond::Atom {
                inner: atom.clone(),
                degree: self.inner.topology().adjacency.neighbors_of(index).len(),
                metadata: metadata
                    .as_ref()
                    .map(|rows| rows[index].clone())
                    .map_err(Clone::clone),
            })
            .collect()
    }

    /// Return the read-only Atom at the zero-based graph index, or None if the index is out of range.
    fn atom(&self, atom_id: usize) -> Option<crate::canonical_atom_bond::Atom> {
        let inner = self.inner.atom(ck::AtomId::new(atom_id))?.clone();
        let metadata = self.inner.atom_metadata(true);
        Some(crate::canonical_atom_bond::Atom {
            inner,
            degree: self.inner.topology().adjacency.neighbors_of(atom_id).len(),
            metadata: metadata.map(|rows| rows[atom_id].clone()),
        })
    }

    /// Return the read-only Bond at the zero-based graph index, or None if the index is out of range.
    fn bond(&self, bond_id: usize) -> Option<crate::canonical_atom_bond::Bond> {
        self.inner
            .bond(ck::BondId::new(bond_id))
            .cloned()
            .map(|inner| crate::canonical_atom_bond::Bond { inner })
    }

    /// Return the stored value for the requested property key, or None if absent.
    fn property(&self, key: &str) -> Option<crate::canonical_property_values::PropertyValue> {
        self.inner
            .property(key)
            .cloned()
            .map(|inner| crate::canonical_property_values::PropertyValue { inner })
    }

    /// Return bond rows in graph order.
    fn bonds(&self) -> Vec<crate::canonical_atom_bond::Bond> {
        self.inner
            .bonds()
            .iter()
            .cloned()
            .map(|inner| crate::canonical_atom_bond::Bond { inner })
            .collect()
    }

    /// Return atom-indexed metadata. recalculate=True requests a fresh calculation; False requires an existing valid cache and raises ValenceError if it is absent or invalid. Does not modify graph topology.
    #[pyo3(signature = (recalculate=true))]
    fn atom_metadata(
        &self,
        py: Python<'_>,
        recalculate: bool,
    ) -> PyResult<Vec<crate::canonical_atom_bond::AtomMetadata>> {
        self.inner
            .atom_metadata(recalculate)
            .map(|rows| {
                rows.into_iter()
                    .map(|inner| crate::canonical_atom_bond::AtomMetadata { inner })
                    .collect()
            })
            .map_err(|e| crate::canonical_atom_bond::valence_pyerr(py, e))
    }

    fn __len__(&self) -> usize {
        self.inner.num_atoms()
    }

    fn __repr__(&self) -> String {
        format!(
            "Molecule(num_atoms={}, num_bonds={})",
            self.inner.num_atoms(),
            self.inner.num_bonds()
        )
    }

    /// In place, perform the selected chemical sanitization stages. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn sanitize_(&mut self, py: Python<'_>) -> PyResult<()> {
        self.inner.sanitize_().map_err(|e| operation_pyerr(py, e))
    }

    /// Apply to a new molecule and return the result: assign CIP stereochemical descriptors. The source molecule is unchanged.
    fn with_cip_labels(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .with_cip_labels()
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }

    /// Apply to a new molecule and return the result: assign CIP stereochemical descriptors. The source molecule is unchanged.
    fn with_cip_labels_with_options(
        &self,
        py: Python<'_>,
        options: &crate::canonical_atom_bond::CipLabelOptions,
    ) -> PyResult<Self> {
        self.inner
            .with_cip_labels_with_options(&options.inner)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }

    /// In place, assign CIP stereochemical descriptors. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn assign_cip_labels_(&mut self, py: Python<'_>) -> PyResult<()> {
        self.inner
            .assign_cip_labels_()
            .map_err(|e| operation_pyerr(py, e))
    }

    /// In place, assign CIP stereochemical descriptors. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn assign_cip_labels_with_options_(
        &mut self,
        py: Python<'_>,
        options: &crate::canonical_atom_bond::CipLabelOptions,
    ) -> PyResult<()> {
        self.inner
            .assign_cip_labels_with_options_(&options.inner)
            .map_err(|e| operation_pyerr(py, e))
    }

    /// Return candidate stereocenters/bonds and the associated stereo-perception result.
    fn potential_stereo(
        &self,
        py: Python<'_>,
    ) -> PyResult<crate::canonical_potential_stereo::PotentialStereoResult> {
        self.inner
            .potential_stereo()
            .map(|inner| crate::canonical_potential_stereo::PotentialStereoResult { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Return candidate stereocenters/bonds and the associated stereo-perception result. Uses the supplied configuration object.
    fn potential_stereo_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_potential_stereo::PotentialStereoParams,
    ) -> PyResult<crate::canonical_potential_stereo::PotentialStereoResult> {
        self.inner
            .potential_stereo_with_params(&params.inner)
            .map(|inner| crate::canonical_potential_stereo::PotentialStereoResult { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Return whether the molecule carries the CIP-computed marker.
    fn cip_computed(&self, py: Python<'_>) -> PyResult<bool> {
        self.inner
            .cip_computed()
            .map_err(|error| crate::canonical_atom_bond::property_pyerr(py, error))
    }

    /// Create an owned persistent MMFF evaluator from the selected stored conformer. Later evaluator edits do not change this molecule.
    #[pyo3(signature=(*,conformer_id=crate::persistent_forcefields::default_conformer(),mmff_variant=crate::persistent_forcefields::default_variant(),non_bonded_threshold=crate::persistent_forcefields::default_non_bonded(),ignore_interfragment_interactions=crate::persistent_forcefields::default_mmff_ignore()))]
    fn mmff_force_field(
        &self,
        py: Python<'_>,
        conformer_id: Option<usize>,
        mmff_variant: String,
        non_bonded_threshold: f64,
        ignore_interfragment_interactions: bool,
    ) -> PyResult<crate::persistent_forcefields::MolecularForceField> {
        let params = ::cosmolkit::MmffForceFieldParams::new(
            conformer_id,
            mmff_variant,
            non_bonded_threshold,
            ignore_interfragment_interactions,
        );
        self.inner
            .mmff_force_field_with_params(&params)
            .map(|inner| crate::persistent_forcefields::MolecularForceField { inner })
            .map_err(|e| crate::persistent_forcefields::mmff_pyerr(py, e))
    }
    /// Create an owned persistent MMFF evaluator from the selected stored conformer. Later evaluator edits do not change this molecule. Uses the supplied configuration object.
    fn mmff_force_field_with_params(
        &self,
        py: Python<'_>,
        params: &crate::persistent_forcefields::MmffForceFieldParams,
    ) -> PyResult<crate::persistent_forcefields::MolecularForceField> {
        self.inner
            .mmff_force_field_with_params(&params.inner)
            .map(|inner| crate::persistent_forcefields::MolecularForceField { inner })
            .map_err(|e| crate::persistent_forcefields::mmff_pyerr(py, e))
    }
    /// Create an owned persistent UFF evaluator from the selected stored conformer. Later evaluator edits do not change this molecule.
    #[pyo3(signature=(*,conformer_id=crate::persistent_forcefields::default_uff_conformer(),vdw_threshold=crate::persistent_forcefields::default_vdw(),ignore_interfragment_interactions=crate::persistent_forcefields::default_uff_ignore()))]
    fn uff_force_field(
        &self,
        py: Python<'_>,
        conformer_id: Option<usize>,
        vdw_threshold: f64,
        ignore_interfragment_interactions: bool,
    ) -> PyResult<crate::persistent_forcefields::MolecularForceField> {
        let params = ::cosmolkit::UffForceFieldParams::new(
            conformer_id,
            vdw_threshold,
            ignore_interfragment_interactions,
        );
        self.inner
            .uff_force_field_with_params(&params)
            .map(|inner| crate::persistent_forcefields::MolecularForceField { inner })
            .map_err(|e| crate::persistent_forcefields::uff_pyerr(py, e))
    }
    /// Create an owned persistent UFF evaluator from the selected stored conformer. Later evaluator edits do not change this molecule. Uses the supplied configuration object.
    fn uff_force_field_with_params(
        &self,
        py: Python<'_>,
        params: &crate::persistent_forcefields::UffForceFieldParams,
    ) -> PyResult<crate::persistent_forcefields::MolecularForceField> {
        self.inner
            .uff_force_field_with_params(&params.inner)
            .map(|inner| crate::persistent_forcefields::MolecularForceField { inner })
            .map_err(|e| crate::persistent_forcefields::uff_pyerr(py, e))
    }
    /// Evaluate UFF energy in kcal/mol and atom-ordered Cartesian energy derivatives at the selected stored coordinates.
    fn uff_energy_gradient(
        &self,
        py: Python<'_>,
    ) -> PyResult<crate::uff_binding::UffEnergyGradient> {
        self.inner
            .uff_energy_gradient()
            .map(|inner| crate::uff_binding::UffEnergyGradient { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Evaluate UFF energy in kcal/mol and atom-ordered Cartesian energy derivatives at the selected stored coordinates. Uses the supplied configuration object.
    fn uff_energy_gradient_with_params(
        &self,
        py: Python<'_>,
        params: &crate::uff_binding::UffEvaluationParams,
    ) -> PyResult<crate::uff_binding::UffEnergyGradient> {
        self.inner
            .uff_energy_gradient_with_params(&params.inner)
            .map(|inner| crate::uff_binding::UffEnergyGradient { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Return whether UFF parameters are available for all atoms in this molecule.
    fn uff_has_all_molecule_params(&self, py: Python<'_>) -> PyResult<bool> {
        self.inner
            .uff_has_all_molecule_params()
            .map_err(|e| crate::uff_binding::parameter_query_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: calculate and install the atom valence assignment. The source molecule is unchanged.
    fn with_assigned_valence(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .with_assigned_valence()
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: calculate and install the atom valence assignment. The source molecule is unchanged.
    fn with_assigned_valence_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_valence::ValenceParams,
    ) -> PyResult<Self> {
        self.inner
            .with_assigned_valence_with_params(&params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, calculate and install the atom valence assignment. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn assign_valence_(&mut self, py: Python<'_>) -> PyResult<()> {
        self.inner
            .assign_valence_()
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, calculate and install the atom valence assignment. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn assign_valence_with_params_(
        &mut self,
        py: Python<'_>,
        params: &crate::canonical_valence::ValenceParams,
    ) -> PyResult<()> {
        self.inner
            .assign_valence_with_params_(&params.inner)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Return whether the atom has a valence violation under the current chemical state.
    fn has_valence_violation(&self, py: Python<'_>, atom_id: usize) -> PyResult<bool> {
        self.inner
            .has_valence_violation(ck::AtomId::new(atom_id))
            .map_err(|e| crate::canonical_atom_bond::valence_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: assign explicit single/double bonds to aromatic systems. The source molecule is unchanged.
    fn with_kekulized_bonds(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .with_kekulized_bonds()
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, assign explicit single/double bonds to aromatic systems. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn kekulize_bonds_(&mut self, py: Python<'_>) -> PyResult<()> {
        self.inner
            .kekulize_bonds_()
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: assign explicit single/double bonds to aromatic systems. The source molecule is unchanged.
    fn with_kekulized_bonds_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_chemistry_values::KekulizeParams,
    ) -> PyResult<Self> {
        self.inner
            .with_kekulized_bonds_with_params(&params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, assign explicit single/double bonds to aromatic systems. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn kekulize_bonds_with_params_(
        &mut self,
        py: Python<'_>,
        params: &crate::canonical_chemistry_values::KekulizeParams,
    ) -> PyResult<()> {
        self.inner
            .kekulize_bonds_with_params_(&params.inner)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: assign aromatic atom and bond flags using the selected aromaticity model. The source molecule is unchanged.
    fn with_assigned_aromaticity(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .with_assigned_aromaticity()
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, assign aromatic atom and bond flags using the selected aromaticity model. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn assign_aromaticity_(&mut self, py: Python<'_>) -> PyResult<()> {
        self.inner
            .assign_aromaticity_()
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: assign aromatic atom and bond flags using the selected aromaticity model. The source molecule is unchanged.
    fn with_assigned_aromaticity_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_chemistry_values::AromaticityParams,
    ) -> PyResult<Self> {
        self.inner
            .with_assigned_aromaticity_with_params(&params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, assign aromatic atom and bond flags using the selected aromaticity model. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn assign_aromaticity_with_params_(
        &mut self,
        py: Python<'_>,
        params: &crate::canonical_chemistry_values::AromaticityParams,
    ) -> PyResult<()> {
        self.inner
            .assign_aromaticity_with_params_(&params.inner)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: assign atomic radical electron counts. The source molecule is unchanged.
    fn with_assigned_radicals(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .with_assigned_radicals()
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, assign atomic radical electron counts. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn assign_radicals_(&mut self, py: Python<'_>) -> PyResult<()> {
        self.inner
            .assign_radicals_()
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: calculate and install ring membership and ring information. The source molecule is unchanged.
    fn with_assigned_rings(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .with_assigned_rings()
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, calculate and install ring membership and ring information. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn assign_rings_(&mut self, py: Python<'_>) -> PyResult<()> {
        self.inner
            .assign_rings_()
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: calculate and install unique ring families. The source molecule is unchanged.
    fn with_assigned_ring_families(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .with_assigned_ring_families()
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, calculate and install unique ring families. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn assign_ring_families_(&mut self, py: Python<'_>) -> PyResult<()> {
        self.inner
            .assign_ring_families_()
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: calculate and install unique ring families. The source molecule is unchanged.
    fn with_assigned_ring_families_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_chemistry_values::RingSearchParams,
    ) -> PyResult<Self> {
        self.inner
            .with_assigned_ring_families_with_params(&params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, calculate and install unique ring families. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn assign_ring_families_with_params_(
        &mut self,
        py: Python<'_>,
        params: &crate::canonical_chemistry_values::RingSearchParams,
    ) -> PyResult<()> {
        self.inner
            .assign_ring_families_with_params_(&params.inner)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: assign atom chiral tags from the selected stored conformer. The source molecule is unchanged.
    fn with_chiral_tags_from_structure(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .with_chiral_tags_from_structure()
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, assign atom chiral tags from the selected stored conformer. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn assign_chiral_tags_from_structure_(&mut self, py: Python<'_>) -> PyResult<()> {
        self.inner
            .assign_chiral_tags_from_structure_()
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: assign atom chiral tags from the selected stored conformer. The source molecule is unchanged.
    fn with_chiral_tags_from_structure_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_chemistry_values::StructureTagParams,
    ) -> PyResult<Self> {
        self.inner
            .with_chiral_tags_from_structure_with_params(&params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, assign atom chiral tags from the selected stored conformer. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn assign_chiral_tags_from_structure_with_params_(
        &mut self,
        py: Python<'_>,
        params: &crate::canonical_chemistry_values::StructureTagParams,
    ) -> PyResult<()> {
        self.inner
            .assign_chiral_tags_from_structure_with_params_(&params.inner)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: remove eligible explicit hydrogen atoms using the selected removal policy. The source molecule is unchanged.
    fn without_hydrogens(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .without_hydrogens()
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, remove eligible explicit hydrogen atoms using the selected removal policy. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn remove_hydrogens_(&mut self, py: Python<'_>) -> PyResult<()> {
        self.inner
            .remove_hydrogens_()
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: remove eligible explicit hydrogen atoms using the selected removal policy. The source molecule is unchanged.
    fn without_hydrogens_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_chemistry_values::RemoveHsParams,
    ) -> PyResult<Self> {
        self.inner
            .without_hydrogens_with_params(&params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, remove eligible explicit hydrogen atoms using the selected removal policy. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn remove_hydrogens_with_params_(
        &mut self,
        py: Python<'_>,
        params: &crate::canonical_chemistry_values::RemoveHsParams,
    ) -> PyResult<()> {
        self.inner
            .remove_hydrogens_with_params_(&params.inner)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: add explicit hydrogen atoms with the selected hydrogen/coordinate options. The source molecule is unchanged.
    fn with_hydrogens_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_chemistry_values::AddHsParams,
    ) -> PyResult<Self> {
        self.inner
            .with_hydrogens_with_params(&params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, add explicit hydrogen atoms with the selected hydrogen/coordinate options. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn add_hydrogens_(&mut self, py: Python<'_>) -> PyResult<()> {
        self.inner
            .add_hydrogens_()
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, add explicit hydrogen atoms with the selected hydrogen/coordinate options. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn add_hydrogens_with_params_(
        &mut self,
        py: Python<'_>,
        params: &crate::canonical_chemistry_values::AddHsParams,
    ) -> PyResult<()> {
        self.inner
            .add_hydrogens_with_params_(&params.inner)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: embed one 3D conformer using the selected distance-geometry parameters. The source molecule is unchanged.
    fn with_3d_conformer_with_params(
        &self,
        py: Python<'_>,
        params: &EmbedParams,
    ) -> PyResult<Self> {
        self.inner
            .with_3d_conformer_with_params(&params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, embed one 3D conformer using the selected distance-geometry parameters. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn embed_3d_conformer_with_params_(
        &mut self,
        py: Python<'_>,
        params: &EmbedParams,
    ) -> PyResult<()> {
        self.inner
            .embed_3d_conformer_with_params_(&params.inner)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: embed one 3D conformer and report its ID and updated embedding parameters. The source molecule is unchanged.
    fn with_3d_conformer_result_with_params(
        &self,
        py: Python<'_>,
        params: &EmbedParams,
    ) -> PyResult<EmbedMoleculeResult> {
        self.inner
            .with_3d_conformer_result_with_params(&params.inner)
            .map(|inner| EmbedMoleculeResult { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, embed one 3D conformer and report its ID and updated embedding parameters. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn embed_3d_conformer_result_with_params_(
        &mut self,
        py: Python<'_>,
        params: &EmbedParams,
    ) -> PyResult<EmbedMoleculeResult> {
        self.inner
            .embed_3d_conformer_result_with_params_(&params.inner)
            .map(|inner| EmbedMoleculeResult { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: embed the requested number of 3D conformers using the selected distance-geometry parameters. The source molecule is unchanged.
    fn with_3d_conformers_with_params(
        &self,
        py: Python<'_>,
        num_confs: u32,
        params: &EmbedParams,
    ) -> PyResult<Self> {
        self.inner
            .with_3d_conformers_with_params(num_confs, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, embed the requested number of 3D conformers using the selected distance-geometry parameters. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn embed_3d_conformers_with_params_(
        &mut self,
        py: Python<'_>,
        num_confs: u32,
        params: &EmbedParams,
    ) -> PyResult<()> {
        self.inner
            .embed_3d_conformers_with_params_(num_confs, &params.inner)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: embed multiple 3D conformers and report generated IDs and embedding parameters. The source molecule is unchanged.
    fn with_3d_conformers_result_with_params(
        &self,
        py: Python<'_>,
        num_confs: u32,
        params: &EmbedParams,
    ) -> PyResult<EmbedMultipleConfsResult> {
        self.inner
            .with_3d_conformers_result_with_params(num_confs, &params.inner)
            .map(|inner| EmbedMultipleConfsResult { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, embed multiple 3D conformers and report generated IDs and embedding parameters. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn embed_3d_conformers_result_with_params_(
        &mut self,
        py: Python<'_>,
        num_confs: u32,
        params: &EmbedParams,
    ) -> PyResult<EmbedMultipleConfsResult> {
        self.inner
            .embed_3d_conformers_result_with_params_(num_confs, &params.inner)
            .map(|inner| EmbedMultipleConfsResult { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Return the atom-ordered graph-distance matrix using the selected bond/atom weight options.
    fn distance_matrix(
        &self,
        py: Python<'_>,
    ) -> PyResult<crate::canonical_chemistry_values::DenseMatrix> {
        self.inner
            .distance_matrix()
            .map(|inner| crate::canonical_chemistry_values::DenseMatrix { inner })
            .map_err(|e| crate::canonical_chemistry_values::matrix_pyerr(py, e))
    }
    /// Return the atom-ordered graph-distance matrix using the selected bond/atom weight options. Uses the supplied configuration object.
    fn distance_matrix_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_chemistry_values::DistanceMatrixParams,
    ) -> PyResult<crate::canonical_chemistry_values::DenseMatrix> {
        self.inner
            .distance_matrix_with_params(&params.inner)
            .map(|inner| crate::canonical_chemistry_values::DenseMatrix { inner })
            .map_err(|e| crate::canonical_chemistry_values::matrix_pyerr(py, e))
    }
    /// Return the atom-ordered Euclidean distance matrix for the selected stored 3D conformer, in angstroms.
    fn distance_matrix_3d(
        &self,
        py: Python<'_>,
    ) -> PyResult<crate::canonical_chemistry_values::DenseMatrix> {
        self.inner
            .distance_matrix_3d()
            .map(|inner| crate::canonical_chemistry_values::DenseMatrix { inner })
            .map_err(|e| crate::canonical_chemistry_values::matrix_pyerr(py, e))
    }
    /// Return the atom-ordered Euclidean distance matrix for the selected stored 3D conformer, in angstroms. Uses the supplied configuration object.
    fn distance_matrix_3d_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_chemistry_values::DistanceMatrix3dParams,
    ) -> PyResult<crate::canonical_chemistry_values::DenseMatrix> {
        self.inner
            .distance_matrix_3d_with_params(&params.inner)
            .map(|inner| crate::canonical_chemistry_values::DenseMatrix { inner })
            .map_err(|e| crate::canonical_chemistry_values::matrix_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: replace the position of one atom in the selected stored conformer. The source molecule is unchanged.
    fn with_atom_position(
        &self,
        py: Python<'_>,
        atom: usize,
        position: [f64; 3],
    ) -> PyResult<Self> {
        self.inner
            .with_atom_position(ck::AtomId::new(atom), position)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: replace the position of one atom in the selected stored conformer. The source molecule is unchanged.
    fn with_atom_position_with_params(
        &self,
        py: Python<'_>,
        atom: usize,
        position: [f64; 3],
        params: &crate::canonical_chemistry_values::AtomPositionParams,
    ) -> PyResult<Self> {
        self.inner
            .with_atom_position_with_params(ck::AtomId::new(atom), position, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Replace the requested atom position in place using the selected stored conformer; invalid input raises an error before commit.
    fn set_atom_position_(
        &mut self,
        py: Python<'_>,
        atom: usize,
        position: [f64; 3],
    ) -> PyResult<()> {
        self.inner
            .set_atom_position_(ck::AtomId::new(atom), position)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Replace the requested atom position in place using the selected stored conformer; invalid input raises an error before commit.
    fn set_atom_position_with_params_(
        &mut self,
        py: Python<'_>,
        atom: usize,
        position: [f64; 3],
        params: &crate::canonical_chemistry_values::AtomPositionParams,
    ) -> PyResult<()> {
        self.inner
            .set_atom_position_with_params_(ck::AtomId::new(atom), position, &params.inner)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: perform the selected chemical sanitization stages. The source molecule is unchanged.
    fn sanitize_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_chemistry_values::SanitizeParams,
    ) -> PyResult<Self> {
        self.inner
            .sanitize_with_params(&params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, perform the selected chemical sanitization stages. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn sanitize_with_params_(
        &mut self,
        py: Python<'_>,
        params: &crate::canonical_chemistry_values::SanitizeParams,
    ) -> PyResult<()> {
        self.inner
            .sanitize_with_params_(&params.inner)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Return structured chemistry problems without committing a sanitized molecule.
    fn detect_chemistry_problems(
        &self,
        py: Python<'_>,
    ) -> PyResult<crate::canonical_chemistry_values::ChemistryProblemReport> {
        self.inner
            .detect_chemistry_problems()
            .map(|inner| crate::canonical_chemistry_values::ChemistryProblemReport { inner })
            .map_err(|e| crate::canonical_chemistry_values::sanitize_pyerr(py, e))
    }
    /// Return structured chemistry problems without committing a sanitized molecule. Uses the supplied configuration object.
    fn detect_chemistry_problems_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_chemistry_values::SanitizeParams,
    ) -> PyResult<crate::canonical_chemistry_values::ChemistryProblemReport> {
        self.inner
            .detect_chemistry_problems_with_params(&params.inner)
            .map(|inner| crate::canonical_chemistry_values::ChemistryProblemReport { inner })
            .map_err(|e| crate::canonical_chemistry_values::sanitize_pyerr(py, e))
    }
    /// Optimize the selected 3D conformer with UFF and return UffOptimizationResult containing the new molecule, convergence status and final energy. The source molecule is unchanged.
    fn with_uff_optimized(
        &self,
        py: Python<'_>,
    ) -> PyResult<crate::uff_binding::UffOptimizationResult> {
        self.inner
            .with_uff_optimized()
            .map(|inner| crate::uff_binding::UffOptimizationResult { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Optimize the selected 3D conformer with UFF and return UffOptimizationResult containing the new molecule, convergence status and final energy. The source molecule is unchanged. Uses the supplied configuration object.
    fn with_uff_optimized_with_params(
        &self,
        py: Python<'_>,
        params: &crate::uff_binding::UffOptimizationParams,
    ) -> PyResult<crate::uff_binding::UffOptimizationResult> {
        self.inner
            .with_uff_optimized_with_params(&params.inner)
            .map(|inner| crate::uff_binding::UffOptimizationResult { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Optimize stored 3D conformers with UFF and return a new molecule with per-conformer convergence statuses and energies. The source molecule is unchanged.
    fn with_uff_optimized_conformers(
        &self,
        py: Python<'_>,
    ) -> PyResult<crate::uff_binding::UffConformerOptimizationResult> {
        self.inner
            .with_uff_optimized_conformers()
            .map(|inner| crate::uff_binding::UffConformerOptimizationResult { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Optimize stored 3D conformers with UFF and return a new molecule with per-conformer convergence statuses and energies. The source molecule is unchanged. Uses the supplied configuration object.
    fn with_uff_optimized_conformers_with_params(
        &self,
        py: Python<'_>,
        params: &crate::uff_binding::UffConformerOptimizationParams,
    ) -> PyResult<crate::uff_binding::UffConformerOptimizationResult> {
        self.inner
            .with_uff_optimized_conformers_with_params(&params.inner)
            .map(|inner| crate::uff_binding::UffConformerOptimizationResult { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Evaluate MMFF energy in kcal/mol and atom-ordered Cartesian energy derivatives at the selected stored coordinates.
    fn mmff_energy_gradient(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<crate::mmff_binding::MmffEnergyGradient>> {
        self.inner
            .mmff_energy_gradient()
            .map(|value| value.map(|inner| crate::mmff_binding::MmffEnergyGradient { inner }))
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Evaluate MMFF energy in kcal/mol and atom-ordered Cartesian energy derivatives at the selected stored coordinates. Uses the supplied configuration object.
    fn mmff_energy_gradient_with_params(
        &self,
        py: Python<'_>,
        params: &crate::mmff_binding::MmffEvaluationParams,
    ) -> PyResult<Option<crate::mmff_binding::MmffEnergyGradient>> {
        self.inner
            .mmff_energy_gradient_with_params(&params.inner)
            .map(|value| value.map(|inner| crate::mmff_binding::MmffEnergyGradient { inner }))
            .map_err(|e| operation_pyerr(py, e))
    }

    /// Apply to a new molecule and return the result: store caller-supplied atom-ordered 2D coordinates. The source molecule is unchanged.
    fn with_2d_coordinate_block(
        &self,
        py: Python<'_>,
        coordinates: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let rows = crate::canonical_coordinate_input::matrix(
            py,
            coordinates,
            self.inner.num_atoms(),
            true,
        )?;
        self.inner
            .with_2d_coordinate_block(rows)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: store caller-supplied atom-ordered 2D coordinates. The source molecule is unchanged.
    fn with_2d_coordinate_block_with_params(
        &self,
        py: Python<'_>,
        coordinates: &Bound<'_, PyAny>,
        params: &crate::canonical_coordinate_input::Coordinate2DInputParams,
    ) -> PyResult<Self> {
        let rows = crate::canonical_coordinate_input::matrix(
            py,
            coordinates,
            self.inner.num_atoms(),
            true,
        )?;
        self.inner
            .with_2d_coordinate_block_with_params(rows, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, store caller-supplied atom-ordered 2D coordinates. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn set_2d_coordinates_(
        &mut self,
        py: Python<'_>,
        coordinates: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        let rows = crate::canonical_coordinate_input::matrix(
            py,
            coordinates,
            self.inner.num_atoms(),
            true,
        )?;
        self.inner
            .set_2d_coordinates_(rows)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, store caller-supplied atom-ordered 2D coordinates. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn set_2d_coordinates_with_params_(
        &mut self,
        py: Python<'_>,
        coordinates: &Bound<'_, PyAny>,
        params: &crate::canonical_coordinate_input::Coordinate2DInputParams,
    ) -> PyResult<()> {
        let rows = crate::canonical_coordinate_input::matrix(
            py,
            coordinates,
            self.inner.num_atoms(),
            true,
        )?;
        self.inner
            .set_2d_coordinates_with_params_(rows, &params.inner)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: replace coordinates of the selected stored 3D conformer. The source molecule is unchanged.
    fn with_3d_coordinates(
        &self,
        py: Python<'_>,
        coordinates: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let rows = crate::canonical_coordinate_input::matrix(
            py,
            coordinates,
            self.inner.num_atoms(),
            false,
        )?;
        self.inner
            .with_3d_coordinates(rows)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: replace coordinates of the selected stored 3D conformer. The source molecule is unchanged.
    fn with_3d_coordinates_with_params(
        &self,
        py: Python<'_>,
        coordinates: &Bound<'_, PyAny>,
        params: &crate::canonical_coordinate_input::Replace3DCoordinatesParams,
    ) -> PyResult<Self> {
        let rows = crate::canonical_coordinate_input::matrix(
            py,
            coordinates,
            self.inner.num_atoms(),
            false,
        )?;
        self.inner
            .with_3d_coordinates_with_params(rows, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, replace coordinates of the selected stored 3D conformer. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn set_3d_coordinates_(
        &mut self,
        py: Python<'_>,
        coordinates: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        let rows = crate::canonical_coordinate_input::matrix(
            py,
            coordinates,
            self.inner.num_atoms(),
            false,
        )?;
        self.inner
            .set_3d_coordinates_(rows)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, replace coordinates of the selected stored 3D conformer. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn set_3d_coordinates_with_params_(
        &mut self,
        py: Python<'_>,
        coordinates: &Bound<'_, PyAny>,
        params: &crate::canonical_coordinate_input::Replace3DCoordinatesParams,
    ) -> PyResult<()> {
        let rows = crate::canonical_coordinate_input::matrix(
            py,
            coordinates,
            self.inner.num_atoms(),
            false,
        )?;
        self.inner
            .set_3d_coordinates_with_params_(rows, &params.inner)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: append a conformer using caller-supplied atom-ordered 3D coordinates. The source molecule is unchanged.
    fn with_added_3d_conformer(
        &self,
        py: Python<'_>,
        coordinates: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let rows = crate::canonical_coordinate_input::matrix(
            py,
            coordinates,
            self.inner.num_atoms(),
            false,
        )?;
        self.inner
            .with_added_3d_conformer(rows)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: append a conformer using caller-supplied atom-ordered 3D coordinates. The source molecule is unchanged.
    fn with_added_3d_conformer_with_params(
        &self,
        py: Python<'_>,
        coordinates: &Bound<'_, PyAny>,
        params: &crate::canonical_coordinate_input::Coordinate3DInputParams,
    ) -> PyResult<Self> {
        let rows = crate::canonical_coordinate_input::matrix(
            py,
            coordinates,
            self.inner.num_atoms(),
            false,
        )?;
        self.inner
            .with_added_3d_conformer_with_params(rows, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, append a conformer using caller-supplied atom-ordered 3D coordinates. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn add_3d_conformer_(
        &mut self,
        py: Python<'_>,
        coordinates: &Bound<'_, PyAny>,
    ) -> PyResult<usize> {
        let rows = crate::canonical_coordinate_input::matrix(
            py,
            coordinates,
            self.inner.num_atoms(),
            false,
        )?;
        self.inner
            .add_3d_conformer_(rows)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, append a conformer using caller-supplied atom-ordered 3D coordinates. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn add_3d_conformer_with_params_(
        &mut self,
        py: Python<'_>,
        coordinates: &Bound<'_, PyAny>,
        params: &crate::canonical_coordinate_input::Coordinate3DInputParams,
    ) -> PyResult<usize> {
        let rows = crate::canonical_coordinate_input::matrix(
            py,
            coordinates,
            self.inner.num_atoms(),
            false,
        )?;
        self.inner
            .add_3d_conformer_with_params_(rows, &params.inner)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: replace the 3D conformer collection with one caller-supplied conformer. The source molecule is unchanged.
    fn with_only_3d_conformer(
        &self,
        py: Python<'_>,
        coordinates: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let rows = crate::canonical_coordinate_input::matrix(
            py,
            coordinates,
            self.inner.num_atoms(),
            false,
        )?;
        self.inner
            .with_only_3d_conformer(rows)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: replace the 3D conformer collection with one caller-supplied conformer. The source molecule is unchanged.
    fn with_only_3d_conformer_with_params(
        &self,
        py: Python<'_>,
        coordinates: &Bound<'_, PyAny>,
        params: &crate::canonical_coordinate_input::Coordinate3DInputParams,
    ) -> PyResult<Self> {
        let rows = crate::canonical_coordinate_input::matrix(
            py,
            coordinates,
            self.inner.num_atoms(),
            false,
        )?;
        self.inner
            .with_only_3d_conformer_with_params(rows, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, replace the 3D conformer collection with one caller-supplied conformer. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn set_only_3d_conformer_(
        &mut self,
        py: Python<'_>,
        coordinates: &Bound<'_, PyAny>,
    ) -> PyResult<usize> {
        let rows = crate::canonical_coordinate_input::matrix(
            py,
            coordinates,
            self.inner.num_atoms(),
            false,
        )?;
        self.inner
            .set_only_3d_conformer_(rows)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, replace the 3D conformer collection with one caller-supplied conformer. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn set_only_3d_conformer_with_params_(
        &mut self,
        py: Python<'_>,
        coordinates: &Bound<'_, PyAny>,
        params: &crate::canonical_coordinate_input::Coordinate3DInputParams,
    ) -> PyResult<usize> {
        let rows = crate::canonical_coordinate_input::matrix(
            py,
            coordinates,
            self.inner.num_atoms(),
            false,
        )?;
        self.inner
            .set_only_3d_conformer_with_params_(rows, &params.inner)
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: remove all stored 3D conformers while retaining the separate 2D state. The source molecule is unchanged.
    fn with_cleared_3d_conformers(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .with_cleared_3d_conformers()
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// In place, remove all stored 3D conformers while retaining the separate 2D state. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn clear_3d_conformers_(&mut self, py: Python<'_>) -> PyResult<()> {
        self.inner
            .clear_3d_conformers_()
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: add explicit hydrogen atoms with the selected hydrogen/coordinate options. The source molecule is unchanged.
    fn with_hydrogens(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .with_hydrogens()
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Return sanitized connected components in source order, copying conformers.
    /// The source molecule is unchanged; empty input returns an empty list.
    fn fragments(&self, py: Python<'_>) -> PyResult<Vec<Self>> {
        self.inner
            .fragments()
            .map(|values| values.into_iter().map(|inner| Self { inner }).collect())
            .map_err(|error| operation_pyerr(py, error))
    }
    /// Most atoms, last component on ties; empty input raises OperationError.
    fn largest_fragment(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .largest_fragment()
            .map(|inner| Self { inner })
            .map_err(|error| operation_pyerr(py, error))
    }

    /// MolHash MurckoScaffold: prune terminal atoms, returning a new molecule.
    fn murcko_scaffold(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .murcko_scaffold()
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// MolHash ExtendedMurcko: retain adjacent substituents as dummy atoms.
    fn net_scaffold(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .net_scaffold()
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// ChemTransforms MurckoDecompose, including ring-exocyclic double bonds.
    fn murcko_decompose(&self, py: Python<'_>) -> PyResult<Self> {
        self.inner
            .murcko_decompose()
            .map(|inner| Self { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Return tetrahedral stereochemical information for the current molecule.
    fn tetrahedral_stereo(&self, py: Python<'_>) -> PyResult<Vec<Py<PyAny>>> {
        self.inner
            .tetrahedral_stereo()
            .map_err(|e| crate::canonical_stereo_queries::error_pyerr(py, e))?
            .into_iter()
            .map(|row| crate::canonical_stereo_queries::tetrahedral_row(py, row))
            .collect()
    }
    /// Return the stereochemistry-perception result using the selected cleanup/assignment options.
    fn perceive_stereochemistry(&self, py: Python<'_>) -> PyResult<()> {
        self.inner
            .perceive_stereochemistry()
            .map_err(|e| crate::canonical_stereo_queries::error_pyerr(py, e))
    }
    /// Return atom indices and assigned or unspecified stereocenter labels using the selected options.
    #[pyo3(signature = (include_unassigned=false))]
    fn find_chiral_centers(
        &self,
        py: Python<'_>,
        include_unassigned: bool,
    ) -> PyResult<Vec<(usize, String)>> {
        self.inner
            .find_chiral_centers(include_unassigned)
            .map_err(|e| crate::canonical_stereo_queries::error_pyerr(py, e))
    }
    /// Return an editable builder initialized from this molecule; building changes does not modify the source.
    fn to_builder(&self) -> crate::canonical_builder::MoleculeBuilder {
        crate::canonical_builder::MoleculeBuilder {
            inner: self.inner.to_builder(),
        }
    }
    /// Return stored 3D conformer values in storage order, preserving each conformer ID.
    fn conformers_3d(&self) -> Vec<crate::mmff_binding::Conformer3D> {
        self.inner
            .conformers_3d()
            .iter()
            .cloned()
            .map(|inner| crate::mmff_binding::Conformer3D { inner })
            .collect()
    }
    /// Return whether every atom can be parameterized by the selected MMFF variant.
    fn mmff_has_all_molecule_params(&self, py: Python<'_>) -> PyResult<bool> {
        self.inner
            .mmff_has_all_molecule_params()
            .map_err(|e| crate::mmff_binding::properties_pyerr(py, e))
    }
    /// Return MMFF atom types, formal charges and partial charges for this molecule.
    fn mmff_properties(&self, py: Python<'_>) -> PyResult<crate::mmff_binding::MmffProperties> {
        self.inner
            .mmff_properties()
            .map(|inner| crate::mmff_binding::MmffProperties { inner })
            .map_err(|e| crate::mmff_binding::properties_pyerr(py, e))
    }
    /// Return MMFF atom types, formal charges and partial charges for this molecule. Uses the supplied configuration object.
    fn mmff_properties_with_params(
        &self,
        py: Python<'_>,
        params: &crate::mmff_binding::MmffPropertiesParams,
    ) -> PyResult<crate::mmff_binding::MmffProperties> {
        self.inner
            .mmff_properties_with_params(&params.inner)
            .map(|inner| crate::mmff_binding::MmffProperties { inner })
            .map_err(|e| crate::mmff_binding::properties_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: optimize the selected stored 3D conformer with the selected MMFF variant. The source molecule is unchanged.
    fn with_mmff_optimized(
        &self,
        py: Python<'_>,
    ) -> PyResult<crate::mmff_binding::MmffOptimizeMoleculeResult> {
        self.inner
            .with_mmff_optimized()
            .map(|inner| crate::mmff_binding::MmffOptimizeMoleculeResult { inner })
            .map_err(|error| operation_pyerr(py, error))
    }
    /// Apply to a new molecule and return the result: optimize the selected stored 3D conformer with the selected MMFF variant. The source molecule is unchanged.
    fn with_mmff_optimized_with_params(
        &self,
        py: Python<'_>,
        params: &crate::mmff_binding::MmffOptimizationParams,
    ) -> PyResult<crate::mmff_binding::MmffOptimizeMoleculeResult> {
        self.inner
            .with_mmff_optimized_with_params(&params.inner)
            .map(|inner| crate::mmff_binding::MmffOptimizeMoleculeResult { inner })
            .map_err(|e| operation_pyerr(py, e))
    }
    /// Apply to a new molecule and return the result: optimize stored 3D conformers with MMFF and return per-conformer outcomes. The source molecule is unchanged.
    fn with_mmff_optimized_conformers(
        &self,
        py: Python<'_>,
    ) -> PyResult<crate::mmff_binding::MmffOptimizeMoleculeConfsResult> {
        self.inner
            .with_mmff_optimized_conformers()
            .map(|inner| crate::mmff_binding::MmffOptimizeMoleculeConfsResult { inner })
            .map_err(|error| operation_pyerr(py, error))
    }
    /// Apply to a new molecule and return the result: optimize stored 3D conformers with MMFF and return per-conformer outcomes. The source molecule is unchanged.
    fn with_mmff_optimized_conformers_with_params(
        &self,
        py: Python<'_>,
        params: &crate::mmff_binding::MmffConformerOptimizationParams,
    ) -> PyResult<crate::mmff_binding::MmffOptimizeMoleculeConfsResult> {
        self.inner
            .with_mmff_optimized_conformers_with_params(&params.inner)
            .map(|inner| crate::mmff_binding::MmffOptimizeMoleculeConfsResult { inner })
            .map_err(|e| operation_pyerr(py, e))
    }

    /// Return SMARTS text describing this graph using the selected query writer options.
    fn to_smarts(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .to_smarts()
            .map_err(|error| crate::canonical_search::write_pyerr(py, error))
            .and_then(|text| crate::canonical_sdf::decode_source_text(py, &text))
    }

    /// Return SMARTS text describing this graph using the selected query writer options. Uses the supplied configuration object.
    fn to_smarts_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_search::SmartsWriteParams,
    ) -> PyResult<String> {
        self.inner
            .to_smarts_with_params(&params.inner)
            .map_err(|error| crate::canonical_search::write_pyerr(py, error))
            .and_then(|text| crate::canonical_sdf::decode_source_text(py, &text))
    }

    /// Return SMARTS text with supported CX annotations.
    fn to_cx_smarts(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .to_cx_smarts()
            .map_err(|error| crate::canonical_search::write_pyerr(py, error))
            .and_then(|text| crate::canonical_sdf::decode_source_text(py, &text))
    }

    /// Return SMARTS text with supported CX annotations. Uses the supplied configuration object.
    fn to_cx_smarts_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_search::SmartsWriteParams,
    ) -> PyResult<String> {
        self.inner
            .to_cx_smarts_with_params(&params.inner)
            .map_err(|error| crate::canonical_search::write_pyerr(py, error))
            .and_then(|text| crate::canonical_sdf::decode_source_text(py, &text))
    }

    /// Return the first query-to-target MatchResult, or None if there is no match. params and keyword configuration are mutually exclusive.
    fn substruct_match(
        &self,
        py: Python<'_>,
        query: &crate::canonical_search::QueryGraph,
    ) -> PyResult<Option<crate::canonical_search::MatchResult>> {
        self.inner
            .substruct_match(&query.inner)
            .map(|result| result.map(|inner| crate::canonical_search::MatchResult { inner }))
            .map_err(|e| crate::canonical_search::substruct_pyerr(py, e))
    }
    /// Return query-to-target MatchResult values subject to max_matches and uniquify. An unmatched query returns an empty list; callbacks and chirality use the supplied configuration.
    fn substruct_matches(
        &self,
        py: Python<'_>,
        query: &crate::canonical_search::QueryGraph,
    ) -> PyResult<Vec<crate::canonical_search::MatchResult>> {
        self.inner
            .substruct_matches(&query.inner)
            .map(|results| {
                results
                    .into_iter()
                    .map(|inner| crate::canonical_search::MatchResult { inner })
                    .collect()
            })
            .map_err(|e| crate::canonical_search::substruct_pyerr(py, e))
    }
    /// Return query-to-target MatchResult values subject to max_matches and uniquify. An unmatched query returns an empty list; callbacks and chirality use the supplied configuration. Uses the supplied configuration object.
    fn substruct_matches_with_params(
        &self,
        py: Python<'_>,
        query: &crate::canonical_search::QueryGraph,
        params: &crate::canonical_search::SubstructMatchParams,
    ) -> PyResult<Vec<crate::canonical_search::MatchResult>> {
        params
            .with_callbacks(py, &self.inner, |params| {
                self.inner
                    .substruct_matches_with_params(&query.inner, params)
            })
            .map(|results| {
                results
                    .into_iter()
                    .map(|inner| crate::canonical_search::MatchResult { inner })
                    .collect()
            })
    }
    /// Return the first query-to-target MatchResult, or None if there is no match. params and keyword configuration are mutually exclusive. Uses the supplied configuration object.
    fn substruct_match_with_params(
        &self,
        py: Python<'_>,
        query: &crate::canonical_search::QueryGraph,
        params: &crate::canonical_search::SubstructMatchParams,
    ) -> PyResult<Option<crate::canonical_search::MatchResult>> {
        params
            .with_callbacks(py, &self.inner, |params| {
                self.inner.substruct_match_with_params(&query.inner, params)
            })
            .map(|result| result.map(|inner| crate::canonical_search::MatchResult { inner }))
    }
    /// Return whether the query matches this molecule using the supplied matching configuration. Uses the supplied configuration object.
    fn has_substruct_match_with_params(
        &self,
        py: Python<'_>,
        query: &crate::canonical_search::QueryGraph,
        params: &crate::canonical_search::SubstructMatchParams,
    ) -> PyResult<bool> {
        params.with_callbacks(py, &self.inner, |params| {
            self.inner
                .has_substruct_match_with_params(&query.inner, params)
        })
    }
    /// Return whether the query matches this molecule using the supplied matching configuration.
    fn has_substruct_match(
        &self,
        py: Python<'_>,
        query: &crate::canonical_search::QueryGraph,
    ) -> PyResult<bool> {
        self.inner
            .has_substruct_match(&query.inner)
            .map_err(|e| crate::canonical_search::substruct_pyerr(py, e))
    }
    /// Match a reusable CompiledQuery against this molecule and return MatchResult values.
    fn substruct_matches_compiled(
        &self,
        py: Python<'_>,
        query: &crate::canonical_search::CompiledQuery,
    ) -> PyResult<Vec<crate::canonical_search::MatchResult>> {
        self.inner
            .substruct_matches_compiled(&query.inner)
            .map(|results| {
                results
                    .into_iter()
                    .map(|inner| crate::canonical_search::MatchResult { inner })
                    .collect()
            })
            .map_err(|e| crate::canonical_search::match_pyerr(py, e))
    }

    /// Construct a Molecule value from the supplied inputs.
    #[staticmethod]
    fn new() -> Self {
        Self {
            inner: ck::Molecule::new(),
        }
    }

    /// Construct a value from explicit detached parts; required structural consistency is checked at the public boundary.
    #[staticmethod]
    fn from_parts(
        py: Python<'_>,
        topology: &crate::canonical_detached_blocks::TopologyBlock,
        coordinates: &crate::canonical_detached_blocks::CoordinateBlock,
        properties: &crate::canonical_property_values::MoleculeProperties,
    ) -> PyResult<Self> {
        ck::Molecule::from_parts(
            topology.inner.clone(),
            coordinates.inner.clone(),
            properties.inner.clone(),
        )
        .map(Self::from_inner)
        .map_err(|error| operation_pyerr(py, error))
    }

    /// Returns the immutable topology value.
    fn topology(&self) -> crate::canonical_detached_blocks::TopologyBlock {
        crate::canonical_detached_blocks::TopologyBlock {
            inner: self.inner.topology().clone(),
        }
    }

    /// Parse SMILES text into a new Molecule with the selected parsing/sanitization options; invalid input raises SmilesError.
    #[staticmethod]
    fn from_smiles(py: Python<'_>, smiles: crate::text_input::TextInput<'_>) -> PyResult<Self> {
        ck::Molecule::from_smiles(&smiles.as_text()?)
            .map(|inner| Self { inner })
            .map_err(|error| smiles_pyerr(py, error))
    }

    /// Parse SMILES text into a new Molecule with the selected parsing/sanitization options; invalid input raises SmilesError. Uses the supplied configuration object.
    #[staticmethod]
    fn from_smiles_with_params(
        py: Python<'_>,
        input: crate::text_input::TextInput<'_>,
        params: &SmilesParseParams,
    ) -> PyResult<Self> {
        ck::Molecule::from_smiles_with_params(&input.as_text()?, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| smiles_pyerr(py, e))
    }

    /// Number of atoms in the graph or selected structure.
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }

    /// Number of bonds in the graph.
    fn num_bonds(&self) -> usize {
        self.inner.num_bonds()
    }

    /// Return a SMILES string using the selected writer settings; this does not modify the molecule.
    fn to_smiles(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .to_smiles()
            .map_err(|error| smiles_write_pyerr(py, error))
            .and_then(|text| crate::canonical_sdf::decode_source_text(py, &text))
    }

    /// Return a SMILES string using the selected writer settings; this does not modify the molecule. Uses the supplied configuration object.
    fn to_smiles_with_params(
        &self,
        py: Python<'_>,
        params: &SmilesWriteParams,
    ) -> PyResult<String> {
        self.inner
            .to_smiles_with_params(&params.inner)
            .map_err(|e| smiles_write_pyerr(py, e))
            .and_then(|text| crate::canonical_sdf::decode_source_text(py, &text))
    }

    /// Return SMILES with the selected CXSMILES annotations; leave the molecule unchanged.
    fn to_cx_smiles(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .to_cx_smiles()
            .map_err(|e| smiles_write_pyerr(py, e))
            .and_then(|text| crate::canonical_sdf::decode_source_text(py, &text))
    }
    /// Return SMILES for the explicitly selected atom/bond fragment; indices refer to this molecule.
    fn to_fragment_smiles(&self, py: Python<'_>, atoms: Vec<usize>) -> PyResult<String> {
        let atoms = atoms.into_iter().map(ck::AtomId::new).collect::<Vec<_>>();
        self.inner
            .to_fragment_smiles(&atoms)
            .map_err(|e| smiles_write_pyerr(py, e))
            .and_then(|text| crate::canonical_sdf::decode_source_text(py, &text))
    }
    /// Return fragment SMILES with the selected CX annotations; indices refer to this molecule.
    fn to_fragment_cx_smiles(&self, py: Python<'_>, atoms: Vec<usize>) -> PyResult<String> {
        let atoms = atoms.into_iter().map(ck::AtomId::new).collect::<Vec<_>>();
        self.inner
            .to_fragment_cx_smiles(&atoms)
            .map_err(|e| smiles_write_pyerr(py, e))
            .and_then(|text| crate::canonical_sdf::decode_source_text(py, &text))
    }
    /// Return SMILES with the selected CXSMILES annotations; leave the molecule unchanged. Uses the supplied configuration object.
    fn to_cx_smiles_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_smiles_writer::CxSmilesWriteParams,
    ) -> PyResult<String> {
        self.inner
            .to_cx_smiles_with_params(&params.inner)
            .map_err(|e| smiles_write_pyerr(py, e))
            .and_then(|text| crate::canonical_sdf::decode_source_text(py, &text))
    }
    /// Return SMILES for the explicitly selected atom/bond fragment; indices refer to this molecule. Uses the supplied configuration object.
    fn to_fragment_smiles_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_smiles_writer::FragmentSmilesWriteParams,
    ) -> PyResult<String> {
        self.inner
            .to_fragment_smiles_with_params(&params.inner)
            .map_err(|e| smiles_write_pyerr(py, e))
            .and_then(|text| crate::canonical_sdf::decode_source_text(py, &text))
    }
    /// Return fragment SMILES with the selected CX annotations; indices refer to this molecule. Uses the supplied configuration object.
    fn to_fragment_cx_smiles_with_params(
        &self,
        py: Python<'_>,
        params: &crate::canonical_smiles_writer::FragmentCxSmilesWriteParams,
    ) -> PyResult<String> {
        self.inner
            .to_fragment_cx_smiles_with_params(&params.inner)
            .map_err(|e| smiles_write_pyerr(py, e))
            .and_then(|text| crate::canonical_sdf::decode_source_text(py, &text))
    }
    /// Return randomized SMILES strings using the supplied seed and output options; leave the graph unchanged.
    fn to_random_smiles(&self, py: Python<'_>, count: u32, seed: u32) -> PyResult<Vec<String>> {
        self.inner
            .to_random_smiles(count, seed)
            .map_err(|e| smiles_write_pyerr(py, e))?
            .iter()
            .map(|text| crate::canonical_sdf::decode_source_text(py, text))
            .collect()
    }
    /// Return randomized SMILES strings using the supplied seed and output options; leave the graph unchanged. Uses the supplied configuration object.
    fn to_random_smiles_with_params(
        &self,
        py: Python<'_>,
        count: u32,
        seed: u32,
        params: &crate::canonical_smiles_writer::RandomSmilesWriteParams,
    ) -> PyResult<Vec<String>> {
        self.inner
            .to_random_smiles_with_params(count, seed, &params.inner)
            .map_err(|e| smiles_write_pyerr(py, e))?
            .iter()
            .map(|text| crate::canonical_sdf::decode_source_text(py, text))
            .collect()
    }
    /// Compute atom pair fixed-width bit fingerprints for this molecule without changing its graph or coordinates.
    fn fingerprint_atom_pair(&self, py: Python<'_>) -> PyResult<Fingerprint> {
        self.inner
            .fingerprint_atom_pair()
            .map(|inner| Fingerprint { inner })
            .map_err(|error| atom_pair_pyerr(py, error))
    }
    /// Compute atom pair fixed-width bit fingerprints for this molecule without changing its graph or coordinates.
    fn fingerprint_atom_pair_with_params(
        &self,
        py: Python<'_>,
        params: &AtomPairFingerprintParams,
        mut additional_output: Option<PyRefMut<'_, FingerprintAdditionalOutput>>,
    ) -> PyResult<Fingerprint> {
        self.inner
            .fingerprint_atom_pair_with_params(
                &params.inner,
                additional_output.as_mut().map(|output| &mut output.inner),
            )
            .map(|inner| Fingerprint { inner })
            .map_err(|error| atom_pair_pyerr(py, error))
    }
    /// Compute atom pair sparse feature bits for this molecule without changing its graph or coordinates.
    fn fingerprint_atom_pair_sparse(&self, py: Python<'_>) -> PyResult<SparseBitFingerprint> {
        self.inner
            .fingerprint_atom_pair_sparse()
            .map(|inner| SparseBitFingerprint { inner })
            .map_err(|error| atom_pair_pyerr(py, error))
    }
    /// Compute atom pair sparse feature bits for this molecule without changing its graph or coordinates.
    fn fingerprint_atom_pair_sparse_with_params(
        &self,
        py: Python<'_>,
        params: &AtomPairFingerprintParams,
        mut additional_output: Option<PyRefMut<'_, FingerprintAdditionalOutput>>,
    ) -> PyResult<SparseBitFingerprint> {
        self.inner
            .fingerprint_atom_pair_sparse_with_params(
                &params.inner,
                additional_output.as_mut().map(|output| &mut output.inner),
            )
            .map(|inner| SparseBitFingerprint { inner })
            .map_err(|error| atom_pair_pyerr(py, error))
    }
    /// Compute atom pair folded feature counts for this molecule without changing its graph or coordinates.
    fn fingerprint_atom_pair_count(&self, py: Python<'_>) -> PyResult<SparseCountFingerprint32> {
        self.inner
            .fingerprint_atom_pair_count()
            .map(|inner| SparseCountFingerprint32 { inner })
            .map_err(|error| atom_pair_pyerr(py, error))
    }
    /// Compute atom pair folded feature counts for this molecule without changing its graph or coordinates.
    fn fingerprint_atom_pair_count_with_params(
        &self,
        py: Python<'_>,
        params: &AtomPairFingerprintParams,
        mut additional_output: Option<PyRefMut<'_, FingerprintAdditionalOutput>>,
    ) -> PyResult<SparseCountFingerprint32> {
        self.inner
            .fingerprint_atom_pair_count_with_params(
                &params.inner,
                additional_output.as_mut().map(|output| &mut output.inner),
            )
            .map(|inner| SparseCountFingerprint32 { inner })
            .map_err(|error| atom_pair_pyerr(py, error))
    }
    /// Compute atom pair sparse feature counts for this molecule without changing its graph or coordinates.
    fn fingerprint_atom_pair_sparse_count(
        &self,
        py: Python<'_>,
    ) -> PyResult<SparseCountFingerprint> {
        self.inner
            .fingerprint_atom_pair_sparse_count()
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(|error| atom_pair_pyerr(py, error))
    }
    /// Compute atom pair sparse feature counts for this molecule without changing its graph or coordinates.
    fn fingerprint_atom_pair_sparse_count_with_params(
        &self,
        py: Python<'_>,
        params: &AtomPairFingerprintParams,
        mut additional_output: Option<PyRefMut<'_, FingerprintAdditionalOutput>>,
    ) -> PyResult<SparseCountFingerprint> {
        self.inner
            .fingerprint_atom_pair_sparse_count_with_params(
                &params.inner,
                additional_output.as_mut().map(|output| &mut output.inner),
            )
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(|error| atom_pair_pyerr(py, error))
    }
    /// Return feature identifiers for the enumerated topological torsion paths using the supplied atom-count/options.
    fn topological_torsion_ids(&self, py: Python<'_>) -> PyResult<Vec<u64>> {
        self.inner
            .topological_torsion_ids()
            .map_err(|error| topological_torsion_pyerr(py, error))
    }
    /// Return feature identifiers for the enumerated topological torsion paths using the supplied atom-count/options. Uses the supplied configuration object.
    fn topological_torsion_ids_with_params(
        &self,
        py: Python<'_>,
        torsion_atom_count: u32,
    ) -> PyResult<Vec<u64>> {
        self.inner
            .topological_torsion_ids_with_params(torsion_atom_count)
            .map_err(|error| topological_torsion_pyerr(py, error))
    }
    /// Compute topological torsion fixed-width bit fingerprints for this molecule without changing its graph or coordinates.
    fn fingerprint_topological_torsion(&self, py: Python<'_>) -> PyResult<Fingerprint> {
        self.inner
            .fingerprint_topological_torsion()
            .map(|inner| Fingerprint { inner })
            .map_err(|error| topological_torsion_pyerr(py, error))
    }
    /// Compute topological torsion fixed-width bit fingerprints for this molecule without changing its graph or coordinates.
    fn fingerprint_topological_torsion_with_params(
        &self,
        py: Python<'_>,
        params: &TopologicalTorsionFingerprintParams,
        mut additional_output: Option<PyRefMut<'_, FingerprintAdditionalOutput>>,
    ) -> PyResult<Fingerprint> {
        self.inner
            .fingerprint_topological_torsion_with_params(
                &params.inner,
                additional_output.as_mut().map(|output| &mut output.inner),
            )
            .map(|inner| Fingerprint { inner })
            .map_err(|error| topological_torsion_pyerr(py, error))
    }
    /// Compute topological torsion sparse feature bits for this molecule without changing its graph or coordinates.
    fn fingerprint_topological_torsion_sparse(
        &self,
        py: Python<'_>,
    ) -> PyResult<SparseBitFingerprint> {
        self.inner
            .fingerprint_topological_torsion_sparse()
            .map(|inner| SparseBitFingerprint { inner })
            .map_err(|error| topological_torsion_pyerr(py, error))
    }
    /// Compute topological torsion sparse feature bits for this molecule without changing its graph or coordinates.
    fn fingerprint_topological_torsion_sparse_with_params(
        &self,
        py: Python<'_>,
        params: &TopologicalTorsionFingerprintParams,
        mut additional_output: Option<PyRefMut<'_, FingerprintAdditionalOutput>>,
    ) -> PyResult<SparseBitFingerprint> {
        self.inner
            .fingerprint_topological_torsion_sparse_with_params(
                &params.inner,
                additional_output.as_mut().map(|output| &mut output.inner),
            )
            .map(|inner| SparseBitFingerprint { inner })
            .map_err(|error| topological_torsion_pyerr(py, error))
    }
    /// Compute topological torsion folded feature counts for this molecule without changing its graph or coordinates.
    fn fingerprint_topological_torsion_count(
        &self,
        py: Python<'_>,
    ) -> PyResult<SparseCountFingerprint32> {
        self.inner
            .fingerprint_topological_torsion_count()
            .map(|inner| SparseCountFingerprint32 { inner })
            .map_err(|error| topological_torsion_pyerr(py, error))
    }
    /// Compute topological torsion folded feature counts for this molecule without changing its graph or coordinates.
    fn fingerprint_topological_torsion_count_with_params(
        &self,
        py: Python<'_>,
        params: &TopologicalTorsionFingerprintParams,
        mut additional_output: Option<PyRefMut<'_, FingerprintAdditionalOutput>>,
    ) -> PyResult<SparseCountFingerprint32> {
        self.inner
            .fingerprint_topological_torsion_count_with_params(
                &params.inner,
                additional_output.as_mut().map(|output| &mut output.inner),
            )
            .map(|inner| SparseCountFingerprint32 { inner })
            .map_err(|error| topological_torsion_pyerr(py, error))
    }
    /// Compute topological torsion sparse feature counts for this molecule without changing its graph or coordinates.
    fn fingerprint_topological_torsion_sparse_count(
        &self,
        py: Python<'_>,
    ) -> PyResult<SparseCountFingerprint> {
        self.inner
            .fingerprint_topological_torsion_sparse_count()
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(|error| topological_torsion_pyerr(py, error))
    }
    /// Compute topological torsion sparse feature counts for this molecule without changing its graph or coordinates.
    fn fingerprint_topological_torsion_sparse_count_with_params(
        &self,
        py: Python<'_>,
        params: &TopologicalTorsionFingerprintParams,
        mut additional_output: Option<PyRefMut<'_, FingerprintAdditionalOutput>>,
    ) -> PyResult<SparseCountFingerprint> {
        self.inner
            .fingerprint_topological_torsion_sparse_count_with_params(
                &params.inner,
                additional_output.as_mut().map(|output| &mut output.inner),
            )
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(|error| topological_torsion_pyerr(py, error))
    }
    /// Compute morgan fixed-width bit fingerprints for this molecule without changing its graph or coordinates.
    fn fingerprint_morgan(&self, py: Python<'_>) -> PyResult<Fingerprint> {
        self.inner
            .fingerprint_morgan()
            .map(|inner| Fingerprint { inner })
            .map_err(|e| morgan_pyerr(py, e))
    }

    /// Compute morgan sparse feature bits for this molecule without changing its graph or coordinates.
    fn fingerprint_morgan_sparse(&self, py: Python<'_>) -> PyResult<SparseBitFingerprint> {
        self.inner
            .fingerprint_morgan_sparse()
            .map(|inner| SparseBitFingerprint { inner })
            .map_err(|e| morgan_pyerr(py, e))
    }

    /// Compute morgan sparse feature counts for this molecule without changing its graph or coordinates.
    fn fingerprint_morgan_sparse_count(&self, py: Python<'_>) -> PyResult<SparseCountFingerprint> {
        self.inner
            .fingerprint_morgan_sparse_count()
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(|e| morgan_pyerr(py, e))
    }

    /// Compute morgan folded feature counts for this molecule without changing its graph or coordinates.
    fn fingerprint_morgan_count(&self, py: Python<'_>) -> PyResult<SparseCountFingerprint32> {
        self.inner
            .fingerprint_morgan_count()
            .map(|inner| SparseCountFingerprint32 { inner })
            .map_err(|e| morgan_pyerr(py, e))
    }

    /// Compute morgan fixed-width bit fingerprints for this molecule without changing its graph or coordinates.
    fn fingerprint_morgan_with_params(
        &self,
        py: Python<'_>,
        params: &MorganFingerprintParams,
        mut additional_output: Option<PyRefMut<'_, FingerprintAdditionalOutput>>,
    ) -> PyResult<Fingerprint> {
        self.inner
            .fingerprint_morgan_with_params(
                &params.inner,
                additional_output.as_mut().map(|output| &mut output.inner),
            )
            .map(|inner| Fingerprint { inner })
            .map_err(|error| morgan_pyerr(py, error))
    }

    /// Compute morgan sparse feature bits for this molecule without changing its graph or coordinates.
    fn fingerprint_morgan_sparse_with_params(
        &self,
        py: Python<'_>,
        params: &MorganFingerprintParams,
        mut additional_output: Option<PyRefMut<'_, FingerprintAdditionalOutput>>,
    ) -> PyResult<SparseBitFingerprint> {
        self.inner
            .fingerprint_morgan_sparse_with_params(
                &params.inner,
                additional_output.as_mut().map(|output| &mut output.inner),
            )
            .map(|inner| SparseBitFingerprint { inner })
            .map_err(|error| morgan_pyerr(py, error))
    }

    /// Compute morgan folded feature counts for this molecule without changing its graph or coordinates.
    fn fingerprint_morgan_count_with_params(
        &self,
        py: Python<'_>,
        params: &MorganFingerprintParams,
        mut additional_output: Option<PyRefMut<'_, FingerprintAdditionalOutput>>,
    ) -> PyResult<SparseCountFingerprint32> {
        self.inner
            .fingerprint_morgan_count_with_params(
                &params.inner,
                additional_output.as_mut().map(|output| &mut output.inner),
            )
            .map(|inner| SparseCountFingerprint32 { inner })
            .map_err(|error| morgan_pyerr(py, error))
    }

    /// Compute morgan sparse feature counts for this molecule without changing its graph or coordinates.
    fn fingerprint_morgan_sparse_count_with_params(
        &self,
        py: Python<'_>,
        params: &MorganFingerprintParams,
        mut additional_output: Option<PyRefMut<'_, FingerprintAdditionalOutput>>,
    ) -> PyResult<SparseCountFingerprint> {
        self.inner
            .fingerprint_morgan_sparse_count_with_params(
                &params.inner,
                additional_output.as_mut().map(|output| &mut output.inner),
            )
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(|error| morgan_pyerr(py, error))
    }

    /// Return an independent float64 NumPy array of shape (num_atoms, 3) for the exact stored conformer ID. Missing or ambiguous selection raises Coordinate3DReadError.
    #[gen_stub(override_return_type(type_repr = "numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]]", imports = ("numpy", "typing")))]
    #[pyo3(signature = (conformer_id=0))]
    fn coordinates_3d<'py>(
        &self,
        py: Python<'py>,
        conformer_id: usize,
    ) -> PyResult<Bound<'py, numpy::PyArray2<f64>>> {
        self.inner
            .coordinates_3d(conformer_id)
            .map(|rows| crate::canonical_coordinate_input::coordinate_array(py, rows))
            .map_err(|error| crate::canonical_coordinate_input::read_pyerr(py, &error))
    }

    /// Return an independent float64 NumPy array of shape (num_atoms, 2), or None when no 2D conformer is stored. Never generates coordinates or falls back to 3D.
    #[gen_stub(override_return_type(type_repr = "typing.Optional[numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]]]", imports = ("numpy", "typing")))]
    fn coordinates_2d<'py>(&self, py: Python<'py>) -> Option<Bound<'py, numpy::PyArray2<f64>>> {
        self.inner
            .coordinates_2d()
            .map(|rows| crate::canonical_coordinate_input::coordinate_array(py, rows))
    }

    /// Return whether a 2D conformer is stored; do not generate coordinates or use a 3D fallback.
    fn has_2d_coordinates(&self) -> bool {
        self.inner.has_2d_coordinates()
    }

    /// In place, generate and store 2D drawing coordinates. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn compute_2d_coordinates_(&mut self, py: Python<'_>) -> PyResult<()> {
        self.inner
            .compute_2d_coordinates_()
            .map_err(|error| operation_pyerr(py, error))
    }

    /// Copy RDKit graph fields and 3D conformers, as in COSMolKit 0.3.0.
    /// None prepares valence only; True sanitizes; False leaves caches unset.
    #[classmethod]
    #[pyo3(signature = (rdmol, sanitize=None))]
    fn from_rdkit(
        _cls: &Bound<'_, pyo3::types::PyType>,
        rdmol: &Bound<'_, PyAny>,
        sanitize: Option<bool>,
    ) -> PyResult<Self> {
        crate::rdkit_binding::from_rdkit(rdmol, sanitize)
    }

    /// In place, generate and store 2D drawing coordinates. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn compute_2d_coordinates_with_params_(
        &mut self,
        py: Python<'_>,
        params: &Coordinate2DParams,
    ) -> PyResult<()> {
        self.inner
            .compute_2d_coordinates_with_params_(&params.inner)
            .map_err(|error| operation_pyerr(py, error))
    }

    /// Apply to a new molecule and return the result: generate and store 2D drawing coordinates. The source molecule is unchanged.
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

    /// Apply to a new molecule and return the result: generate and store 2D drawing coordinates. The source molecule is unchanged.
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

    /// Write SVG output to the supplied filesystem path/directory using the selected options. This writes files rather than returning serialized text.
    #[pyo3(signature = (path, width, height))]
    fn write_svg(&self, py: Python<'_>, path: &str, width: u32, height: u32) -> PyResult<()> {
        let path = expand_user_path(path)?;
        self.inner
            .write_svg(&path, width, height)
            .map_err(|error| drawing_write_pyerr(py, error))
    }

    /// Write PNG output to the supplied filesystem path/directory using the selected options. This writes files rather than returning serialized text.
    #[pyo3(signature = (path, width, height))]
    fn write_png(&self, py: Python<'_>, path: &str, width: u32, height: u32) -> PyResult<()> {
        let path = expand_user_path(path)?;
        self.inner
            .write_png(&path, width, height)
            .map_err(|error| drawing_write_pyerr(py, error))
    }
    /// Return Wildman-Crippen logP and molar refractivity as CrippenTotals.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn crippen_descriptors(&self, py: Python<'_>) -> PyResult<CrippenTotals> {
        let result = self
            .inner
            .crippen_descriptors()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(CrippenTotals { inner: result })
    }

    /// Return the Labute approximate accessible surface area.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn labute_asa(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .labute_asa()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return total Labute surface area, per-atom contributions and the hydrogen contribution as LabuteAsaContributions.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn labute_asa_contributions(&self, py: Python<'_>) -> PyResult<LabuteAsaContributions> {
        let result = self
            .inner
            .labute_asa_contributions()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(LabuteAsaContributions { inner: result })
    }

    /// Return the topological polar surface area using the selected sulfur/phosphorus inclusion policy.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn tpsa(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .tpsa()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return the 12 surface-area bins grouped by atomic Wildman-Crippen logP contributions.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn slogp_vsa(&self, py: Python<'_>) -> PyResult<Vec<f64>> {
        let result = self
            .inner
            .slogp_vsa()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return the 10 surface-area bins grouped by atomic Wildman-Crippen molar-refractivity contributions.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn smr_vsa(&self, py: Python<'_>) -> PyResult<Vec<f64>> {
        let result = self
            .inner
            .smr_vsa()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 1 of the logP-weighted accessible-surface-area descriptor.
    fn slogp_vsa_1(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_1()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 2 of the logP-weighted accessible-surface-area descriptor.
    fn slogp_vsa_2(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_2()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 3 of the logP-weighted accessible-surface-area descriptor.
    fn slogp_vsa_3(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_3()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 4 of the logP-weighted accessible-surface-area descriptor.
    fn slogp_vsa_4(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_4()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 5 of the logP-weighted accessible-surface-area descriptor.
    fn slogp_vsa_5(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_5()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 6 of the logP-weighted accessible-surface-area descriptor.
    fn slogp_vsa_6(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_6()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 7 of the logP-weighted accessible-surface-area descriptor.
    fn slogp_vsa_7(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_7()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 8 of the logP-weighted accessible-surface-area descriptor.
    fn slogp_vsa_8(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_8()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 9 of the logP-weighted accessible-surface-area descriptor.
    fn slogp_vsa_9(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_9()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 10 of the logP-weighted accessible-surface-area descriptor.
    fn slogp_vsa_10(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_10()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 11 of the logP-weighted accessible-surface-area descriptor.
    fn slogp_vsa_11(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_11()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 12 of the logP-weighted accessible-surface-area descriptor.
    fn slogp_vsa_12(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .slogp_vsa_12()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 1 of the molar-refractivity-weighted accessible-surface-area descriptor.
    fn smr_vsa_1(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_1()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 2 of the molar-refractivity-weighted accessible-surface-area descriptor.
    fn smr_vsa_2(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_2()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 3 of the molar-refractivity-weighted accessible-surface-area descriptor.
    fn smr_vsa_3(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_3()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 4 of the molar-refractivity-weighted accessible-surface-area descriptor.
    fn smr_vsa_4(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_4()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 5 of the molar-refractivity-weighted accessible-surface-area descriptor.
    fn smr_vsa_5(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_5()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 6 of the molar-refractivity-weighted accessible-surface-area descriptor.
    fn smr_vsa_6(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_6()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 7 of the molar-refractivity-weighted accessible-surface-area descriptor.
    fn smr_vsa_7(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_7()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 8 of the molar-refractivity-weighted accessible-surface-area descriptor.
    fn smr_vsa_8(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_8()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 9 of the molar-refractivity-weighted accessible-surface-area descriptor.
    fn smr_vsa_9(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_9()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return bin 10 of the molar-refractivity-weighted accessible-surface-area descriptor.
    fn smr_vsa_10(&self, py: Python<'_>) -> PyResult<f64> {
        let result = self
            .inner
            .smr_vsa_10()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))?;
        Ok(result)
    }

    /// Return Wildman-Crippen logP and molar refractivity as CrippenTotals. Uses the explicit arguments shown in the signature.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
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

    /// Return the Labute approximate accessible surface area. Uses the explicit arguments shown in the signature.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
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

    /// Return total Labute surface area, per-atom contributions and the hydrogen contribution as LabuteAsaContributions. Uses the explicit arguments shown in the signature.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
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

    /// Return the topological polar surface area using the selected sulfur/phosphorus inclusion policy. Uses the explicit arguments shown in the signature.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
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

    /// Return the 12 surface-area bins grouped by atomic Wildman-Crippen logP contributions. Uses the explicit arguments shown in the signature.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
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

    /// Return the 10 surface-area bins grouped by atomic Wildman-Crippen molar-refractivity contributions. Uses the explicit arguments shown in the signature.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
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

    /// Return the quantitative estimate of drug-likeness using the configured descriptor weights.
    ///
    /// This is a read-only molecular query; the source graph and coordinates are unchanged.
    fn qed(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .qed()
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }
    /// Return the 0th-order valence-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    #[pyo3(signature = (force))]
    fn chi_0_v_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_0_v_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the 1th-order valence-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    #[pyo3(signature = (force))]
    fn chi_1_v_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_1_v_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the 2th-order valence-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    #[pyo3(signature = (force))]
    fn chi_2_v_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_2_v_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the 3th-order valence-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    #[pyo3(signature = (force))]
    fn chi_3_v_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_3_v_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the 4th-order valence-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    #[pyo3(signature = (force))]
    fn chi_4_v_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_4_v_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the requested-order valence-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    #[pyo3(signature = (order, force))]
    fn chi_n_v_with_params(&self, py: Python<'_>, order: u32, force: bool) -> PyResult<f64> {
        self.inner
            .chi_n_v_with_params(order, force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the 0th-order n-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    #[pyo3(signature = (force))]
    fn chi_0_n_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_0_n_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the 1th-order n-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    #[pyo3(signature = (force))]
    fn chi_1_n_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_1_n_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the 2th-order n-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    #[pyo3(signature = (force))]
    fn chi_2_n_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_2_n_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the 3th-order n-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    #[pyo3(signature = (force))]
    fn chi_3_n_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_3_n_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the 4th-order n-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    #[pyo3(signature = (force))]
    fn chi_4_n_with_params(&self, py: Python<'_>, force: bool) -> PyResult<f64> {
        self.inner
            .chi_4_n_with_params(force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }

    /// Return the requested-order n-weighted Chi molecular connectivity index. Requires prepared valence; does not change the graph.
    #[pyo3(signature = (order, force))]
    fn chi_n_n_with_params(&self, py: Python<'_>, order: u32, force: bool) -> PyResult<f64> {
        self.inner
            .chi_n_n_with_params(order, force)
            .map_err(|error| crate::canonical_descriptor_binding::descriptor_pyerr(py, error))
    }
    /// Return an alignment result containing RMSD and the transform for the supplied atom correspondence; do not move either molecule.
    #[pyo3(signature = (reference, params=None))]
    #[doc = r#"
Compute the transform aligning this molecule to ``reference`` without mutation.

The returned result contains the RMSD, 4x4 transform, and selected atom map.
Use ``with_alignment_to()`` or ``align_to_()`` to apply the transform.
"#]
    fn alignment_transform_to(
        &self,
        reference: &Molecule,
        params: Option<&PyAlignmentParameters>,
    ) -> PyResult<PyAlignmentResult> {
        let params = params
            .map(|params| params.wrapper_parameters(self.inner.num_atoms()))
            .transpose()
            .map_err(|err| Python::attach(|py| crate::alignment_binding::alignment_pyerr(py, err)))?
            .unwrap_or_default();
        self.inner
            .alignment_transform_to_with_params(&reference.inner, &params)
            .map(Into::into)
            .map_err(|err| Python::attach(|py| crate::alignment_binding::alignment_pyerr(py, err)))
    }

    /// Return the best fitted alignment RMSD and transform across the allowed atom correspondences; leave the molecules unchanged.
    #[pyo3(signature = (reference, params=None))]
    #[doc = "Return the best source-compatible alignment result without mutating either molecule."]
    fn best_alignment_to(
        &self,
        reference: &Molecule,
        params: Option<&PyBestAlignmentParameters>,
    ) -> PyResult<PyAlignmentResult> {
        let params = params
            .map(PyBestAlignmentParameters::wrapper_parameters)
            .transpose()
            .map_err(|err| Python::attach(|py| crate::alignment_binding::alignment_pyerr(py, err)))?
            .unwrap_or_default();
        self.inner
            .best_alignment_to_with_params(&reference.inner, &params)
            .map(Into::into)
            .map_err(|err| Python::attach(|py| crate::alignment_binding::alignment_pyerr(py, err)))
    }

    /// Return the smallest fitted RMSD to the reference over the allowed atom mappings, in angstroms; leave both molecules unchanged.
    #[pyo3(signature = (reference, params=None))]
    #[doc = "Return the best aligned RMSD without changing either molecule's coordinates."]
    fn best_rmsd_to(
        &self,
        reference: &Molecule,
        params: Option<&PyBestAlignmentParameters>,
    ) -> PyResult<f64> {
        let params = params
            .map(PyBestAlignmentParameters::wrapper_parameters)
            .transpose()
            .map_err(|err| Python::attach(|py| crate::alignment_binding::alignment_pyerr(py, err)))?
            .unwrap_or_default();
        self.inner
            .best_rmsd_to_with_params(&reference.inner, &params)
            .map_err(|err| Python::attach(|py| crate::alignment_binding::alignment_pyerr(py, err)))
    }

    /// Return coordinate RMSD to the reference without fitting or modifying either molecule.
    #[pyo3(signature = (reference, params=None))]
    #[doc = r#"
Measure RMSD in the existing coordinate frame without alignment or mutation.

This method corresponds to RDKit ``CalcRMS`` semantics, including map
enumeration and optional terminal-group symmetrization.
"#]
    fn coordinate_rmsd_to(
        &self,
        reference: &Molecule,
        params: Option<&PyCoordinateRmsdParameters>,
    ) -> PyResult<f64> {
        let params = params
            .map(PyCoordinateRmsdParameters::core_parameters)
            .unwrap_or_default();
        self.inner
            .coordinate_rmsd_to_with_params(&reference.inner, &params)
            .map_err(|err| Python::attach(|py| crate::alignment_binding::alignment_pyerr(py, err)))
    }

    /// Return best RMSDs for the stored conformer pairs using the configured mapping/symmetry policy.
    #[pyo3(signature = (params=None))]
    #[doc = "Return best RMSD values for every ordered triangular conformer pair without mutation."]
    fn all_conformer_best_rmsds(
        &self,
        params: Option<&PyAllConformerRmsdParameters>,
    ) -> PyResult<Vec<PyConformerRmsd>> {
        let params = params
            .map(PyAllConformerRmsdParameters::wrapper_parameters)
            .transpose()
            .map_err(|err| Python::attach(|py| crate::alignment_binding::alignment_pyerr(py, err)))?
            .unwrap_or_default();
        self.inner
            .all_conformer_best_rmsds_with_params(&params)
            .map(|values| values.into_iter().map(Into::into).collect())
            .map_err(|err| Python::attach(|py| crate::alignment_binding::alignment_pyerr(py, err)))
    }

    /// Apply to a new molecule and return the result: align the selected conformer to the reference molecule. The source molecule is unchanged.
    #[pyo3(signature = (reference, params=None))]
    #[doc = r#"
Return a new molecule aligned to ``reference`` together with its alignment result.

The source and reference molecules remain unchanged.
"#]
    fn with_alignment_to(
        &self,
        reference: &Molecule,
        params: Option<&PyAlignmentParameters>,
    ) -> PyResult<(Molecule, PyAlignmentResult)> {
        let params = params
            .map(|params| params.wrapper_parameters(self.inner.num_atoms()))
            .transpose()
            .map_err(|err| {
                Python::attach(|py| {
                    operation_pyerr(py, ::cosmolkit::OperationError::Alignment(err))
                })
            })?
            .unwrap_or_default();
        self.inner
            .with_alignment_to_with_params(&reference.inner, &params)
            .map(|(molecule, result)| (Molecule { inner: molecule }, result.into()))
            .map_err(|err| Python::attach(|py| operation_pyerr(py, err)))
    }

    /// In place, align the selected conformer to the reference molecule. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    #[pyo3(signature = (reference, params=None))]
    #[doc = "Align this molecule to ``reference`` in place and return the applied result."]
    fn align_to_<'py>(
        mut slf: PyRefMut<'py, Self>,
        #[gen_stub(override_type(type_repr = "Molecule"))] reference: &Bound<'py, PyAny>,
        params: Option<&PyAlignmentParameters>,
    ) -> PyResult<PyAlignmentResult> {
        let params = params
            .map(|params| params.wrapper_parameters(slf.inner.num_atoms()))
            .transpose()
            .map_err(|err| {
                Python::attach(|py| {
                    operation_pyerr(py, ::cosmolkit::OperationError::Alignment(err))
                })
            })?
            .unwrap_or_default();
        let reference_inner = if slf.as_ptr() == reference.as_ptr() {
            slf.inner.clone()
        } else {
            reference.extract::<PyRef<'_, Molecule>>()?.inner.clone()
        };
        slf.inner
            .align_to_with_params_(&reference_inner, &params)
            .map(Into::into)
            .map_err(|err| Python::attach(|py| operation_pyerr(py, err)))
    }

    /// Apply to a new molecule and return the result: align the selected stored conformers to a reference conformer. The source molecule is unchanged.
    #[pyo3(signature = (params=None))]
    #[doc = "Return a molecule with aligned conformers and the ordered source RMS report."]
    fn with_aligned_conformers(
        &self,
        params: Option<&PyConformerAlignmentParameters>,
    ) -> PyResult<(Molecule, PyConformerAlignmentReport)> {
        let params = params
            .map(PyConformerAlignmentParameters::core_parameters)
            .unwrap_or_default();
        self.inner
            .with_aligned_conformers_with_params(&params)
            .map(|(molecule, report)| (Molecule { inner: molecule }, report.into()))
            .map_err(|err| Python::attach(|py| operation_pyerr(py, err)))
    }

    /// In place, align the selected stored conformers to a reference conformer. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    #[pyo3(signature = (params=None))]
    #[doc = "Align selected or all conformers in place and return the ordered source RMS report."]
    fn align_conformers_(
        &mut self,
        params: Option<&PyConformerAlignmentParameters>,
    ) -> PyResult<PyConformerAlignmentReport> {
        let params = params
            .map(PyConformerAlignmentParameters::core_parameters)
            .unwrap_or_default();
        self.inner
            .align_conformers_with_params_(&params)
            .map(Into::into)
            .map_err(|err| Python::attach(|py| operation_pyerr(py, err)))
    }
    /// Return an alignment result containing RMSD and the transform for the supplied atom correspondence; do not move either molecule. Uses the supplied configuration object.
    fn alignment_transform_to_with_params(
        &self,
        reference: &Molecule,
        params: &PyAlignmentParameters,
    ) -> PyResult<PyAlignmentResult> {
        self.alignment_transform_to(reference, Some(params))
    }

    /// Return the best fitted alignment RMSD and transform across the allowed atom correspondences; leave the molecules unchanged. Uses the supplied configuration object.
    fn best_alignment_to_with_params(
        &self,
        reference: &Molecule,
        params: &PyBestAlignmentParameters,
    ) -> PyResult<PyAlignmentResult> {
        self.best_alignment_to(reference, Some(params))
    }

    /// Return the smallest fitted RMSD to the reference over the allowed atom mappings, in angstroms; leave both molecules unchanged. Uses the supplied configuration object.
    fn best_rmsd_to_with_params(
        &self,
        reference: &Molecule,
        params: &PyBestAlignmentParameters,
    ) -> PyResult<f64> {
        self.best_rmsd_to(reference, Some(params))
    }

    /// Return coordinate RMSD to the reference without fitting or modifying either molecule. Uses the supplied configuration object.
    fn coordinate_rmsd_to_with_params(
        &self,
        reference: &Molecule,
        params: &PyCoordinateRmsdParameters,
    ) -> PyResult<f64> {
        self.coordinate_rmsd_to(reference, Some(params))
    }

    /// Apply to a new molecule and return the result: align the selected conformer to the reference molecule. The source molecule is unchanged.
    fn with_alignment_to_with_params(
        &self,
        reference: &Molecule,
        params: &PyAlignmentParameters,
    ) -> PyResult<(Molecule, PyAlignmentResult)> {
        self.with_alignment_to(reference, Some(params))
    }

    /// Return best RMSDs for the stored conformer pairs using the configured mapping/symmetry policy. Uses the supplied configuration object.
    fn all_conformer_best_rmsds_with_params(
        &self,
        params: &PyAllConformerRmsdParameters,
    ) -> PyResult<Vec<PyConformerRmsd>> {
        self.all_conformer_best_rmsds(Some(params))
    }
    /// Apply to a new molecule and return the result: align the selected stored conformers to a reference conformer. The source molecule is unchanged.
    fn with_aligned_conformers_with_params(
        &self,
        params: &PyConformerAlignmentParameters,
    ) -> PyResult<(Molecule, PyConformerAlignmentReport)> {
        self.with_aligned_conformers(Some(params))
    }
    /// In place, align the selected stored conformers to a reference conformer. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn align_conformers_with_params_(
        &mut self,
        params: &PyConformerAlignmentParameters,
    ) -> PyResult<PyConformerAlignmentReport> {
        self.align_conformers_(Some(params))
    }
    /// In place, align the selected conformer to the reference molecule. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    fn align_to_with_params_<'py>(
        slf: PyRefMut<'py, Self>,
        #[gen_stub(override_type(type_repr = "Molecule"))] reference: &Bound<'py, PyAny>,
        params: &PyAlignmentParameters,
    ) -> PyResult<PyAlignmentResult> {
        Self::align_to_(slf, reference, Some(params))
    }

    /// Return the number of stored 3D conformers; 2D conformers are not included.
    fn num_3d_conformers(&self) -> usize {
        self.inner.num_3d_conformers()
    }
    /// Return the atom-pair lower/upper bounds matrix used by distance-geometry embedding.
    #[gen_stub(override_return_type(type_repr="numpy.ndarray[typing.Any, typing.Any]",imports=("numpy","typing")))]
    fn dg_bounds_matrix<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let rows = self
            .inner
            .dg_bounds_matrix()
            .map_err(|e| crate::canonical_values::source_pyerr(py, &e))?;
        let n = rows.len();
        Array2::from_shape_vec((n, n), rows.into_iter().flatten().collect())
            .map(|a| a.into_pyarray(py).into_any())
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }
    /// Apply to a new molecule and return the result: embed one 3D conformer using the selected distance-geometry parameters. The source molecule is unchanged.
    #[pyo3(signature=(params=None))]
    fn with_3d_conformer(&self, py: Python<'_>, params: Option<&EmbedParams>) -> PyResult<Self> {
        let defaults;
        let options = match params {
            Some(params) => &params.inner,
            None => {
                defaults = ck::EmbedParams::etkdg_v3();
                &defaults
            }
        };
        let outcome = self
            .inner
            .with_3d_conformer_result_with_params(options)
            .map_err(|e| operation_pyerr(py, e))?;
        Ok(Self {
            inner: outcome.molecule,
        })
    }
    /// In place, embed one 3D conformer using the selected distance-geometry parameters. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    #[pyo3(signature=(params=None))]
    fn embed_3d_conformer_(
        &mut self,
        py: Python<'_>,
        params: Option<&EmbedParams>,
    ) -> PyResult<()> {
        let defaults;
        let options = match params {
            Some(params) => &params.inner,
            None => {
                defaults = ck::EmbedParams::etkdg_v3();
                &defaults
            }
        };
        let outcome = self
            .inner
            .embed_3d_conformer_result_with_params_(options)
            .map_err(|e| operation_pyerr(py, e))?;
        Ok(())
    }
    /// Apply to a new molecule and return the result: embed one 3D conformer and report its ID and updated embedding parameters. The source molecule is unchanged.
    #[pyo3(signature=(params=None))]
    fn with_3d_conformer_result(
        &self,
        py: Python<'_>,
        params: Option<&EmbedParams>,
    ) -> PyResult<EmbedMoleculeResult> {
        let defaults;
        let options = match params {
            Some(params) => &params.inner,
            None => {
                defaults = ck::EmbedParams::etkdg_v3();
                &defaults
            }
        };
        let outcome = self
            .inner
            .with_3d_conformer_result_with_params(options)
            .map_err(|e| operation_pyerr(py, e))?;
        Ok(EmbedMoleculeResult { inner: outcome })
    }
    /// In place, embed one 3D conformer and report its ID and updated embedding parameters. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    #[pyo3(signature=(params=None))]
    fn embed_3d_conformer_result_(
        &mut self,
        py: Python<'_>,
        params: Option<&EmbedParams>,
    ) -> PyResult<EmbedMoleculeResult> {
        let defaults;
        let options = match params {
            Some(params) => &params.inner,
            None => {
                defaults = ck::EmbedParams::etkdg_v3();
                &defaults
            }
        };
        let outcome = self
            .inner
            .embed_3d_conformer_result_with_params_(options)
            .map_err(|e| operation_pyerr(py, e))?;
        Ok(EmbedMoleculeResult { inner: outcome })
    }
    /// Apply to a new molecule and return the result: embed the requested number of 3D conformers using the selected distance-geometry parameters. The source molecule is unchanged.
    #[pyo3(signature=(num_confs, params=None))]
    fn with_3d_conformers(
        &self,
        py: Python<'_>,
        num_confs: u32,
        params: Option<&EmbedParams>,
    ) -> PyResult<Self> {
        let defaults;
        let options = match params {
            Some(params) => &params.inner,
            None => {
                defaults = ck::EmbedParams::etkdg_v3();
                &defaults
            }
        };
        let outcome = self
            .inner
            .with_3d_conformers_result_with_params(num_confs, options)
            .map_err(|e| operation_pyerr(py, e))?;
        Ok(Self {
            inner: outcome.molecule,
        })
    }
    /// In place, embed the requested number of 3D conformers using the selected distance-geometry parameters. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    #[pyo3(signature=(num_confs, params=None))]
    fn embed_3d_conformers_(
        &mut self,
        py: Python<'_>,
        num_confs: u32,
        params: Option<&EmbedParams>,
    ) -> PyResult<()> {
        let defaults;
        let options = match params {
            Some(params) => &params.inner,
            None => {
                defaults = ck::EmbedParams::etkdg_v3();
                &defaults
            }
        };
        let outcome = self
            .inner
            .embed_3d_conformers_result_with_params_(num_confs, options)
            .map_err(|e| operation_pyerr(py, e))?;
        Ok(())
    }
    /// Apply to a new molecule and return the result: embed multiple 3D conformers and report generated IDs and embedding parameters. The source molecule is unchanged.
    #[pyo3(signature=(num_confs, params=None))]
    fn with_3d_conformers_result(
        &self,
        py: Python<'_>,
        num_confs: u32,
        params: Option<&EmbedParams>,
    ) -> PyResult<EmbedMultipleConfsResult> {
        let defaults;
        let options = match params {
            Some(params) => &params.inner,
            None => {
                defaults = ck::EmbedParams::etkdg_v3();
                &defaults
            }
        };
        let outcome = self
            .inner
            .with_3d_conformers_result_with_params(num_confs, options)
            .map_err(|e| operation_pyerr(py, e))?;
        Ok(EmbedMultipleConfsResult { inner: outcome })
    }
    /// In place, embed multiple 3D conformers and report generated IDs and embedding parameters. Copy-on-write isolates other molecule values; errors do not commit partial changes.
    #[pyo3(signature=(num_confs, params=None))]
    fn embed_3d_conformers_result_(
        &mut self,
        py: Python<'_>,
        num_confs: u32,
        params: Option<&EmbedParams>,
    ) -> PyResult<EmbedMultipleConfsResult> {
        let defaults;
        let options = match params {
            Some(params) => &params.inner,
            None => {
                defaults = ck::EmbedParams::etkdg_v3();
                &defaults
            }
        };
        let outcome = self
            .inner
            .embed_3d_conformers_result_with_params_(num_confs, options)
            .map_err(|e| operation_pyerr(py, e))?;
        Ok(EmbedMultipleConfsResult { inner: outcome })
    }
}

#[pymodule]
pub(crate) fn cosmolkit(module: &Bound<'_, PyModule>) -> PyResult<()> {
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
    crate::canonical_registered_errors::register(module)?;
    module.add_class::<Molecule>()?;
    module.add_class::<RotatableBondsOptions>()?;
    module.add_class::<CrippenTotals>()?;
    module.add_class::<LabuteAsaContributions>()?;
    module.add_class::<Coordinate2DParams>()?;
    crate::canonical_descriptor_binding::register(module)?;
    crate::canonical_values::register(module)?;
    crate::canonical_smiles_writer::register(module)?;
    crate::canonical_valence::register(module)?;
    crate::canonical_chemistry_values::register(module)?;
    crate::canonical_path_score::register(module)?;
    crate::canonical_maccs::register(module)?;
    crate::canonical_layered::register(module)?;
    crate::canonical_avalon::register(module)?;
    crate::canonical_pattern::register(module)?;
    crate::canonical_topological::register(module)?;
    crate::canonical_element_metadata::register(module)?;
    crate::canonical_operation_metadata::register(module)?;
    crate::mmff_binding::register(module)?;
    crate::uff_binding::register(module)?;
    crate::persistent_forcefields::register(module)?;
    crate::alignment_binding::register(module)?;
    crate::canonical_search::register(module)?;
    crate::canonical_mcs::register(module)?;
    crate::canonical_reaction::register(module)?;
    crate::canonical_sdf::register(module)?;
    crate::canonical_batch::register(module)?;
    crate::canonical_batch_params::register(module)?;
    crate::canonical_batch_fingerprint_values::register(module)?;
    crate::canonical_stereoisomers::register(module)?;
    crate::canonical_atom_bond::register(module)?;
    crate::canonical_potential_stereo::register(module)?;
    crate::canonical_binary::register(module)?;
    crate::canonical_molecular_io::register(module)?;
    crate::canonical_inchi::register(module)?;
    crate::canonical_group_values::register(module)?;
    crate::canonical_sdf::register(module)?;
    crate::canonical_batch::register(module)?;
    crate::canonical_sdf_supplier::register(module)?;
    crate::canonical_molecular_hash::register(module)?;
    crate::canonical_builder::register(module)?;
    crate::canonical_detached_blocks::register(module)?;
    crate::canonical_group_values::register(module)?;
    crate::canonical_coordinate_input::register(module)?;
    crate::canonical_stereo_queries::register(module)?;
    crate::tautomer_binding::register(module)?;
    crate::canonical_property_values::register(module)?;
    crate::canonical_bio_residue::register(module)?;
    crate::canonical_bio_binding::register(module)?;
    crate::canonical_bio_values::register(module)?;
    crate::conformer_binding::register(module)?;
    crate::configuration_projection::register(module)?;
    Ok(())
}
