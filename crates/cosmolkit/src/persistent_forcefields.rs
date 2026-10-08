//! Public persistent force fields: owned detached state, no implicit molecule writeback.
use crate::{AtomId, BINDING_CONTRACT, BindingDefault, Molecule};
pub use cosmolkit_forcefields::{ForceFieldError, MolecularForceFieldErrorKind};
use std::{error::Error, fmt, sync::Mutex};

// Defaults have one declaration in the compiler-checked binding contract.
fn contract_default(constructor: &str, name: &str) -> &'static str {
    let entry = BINDING_CONTRACT
        .iter()
        .find(|entry| entry.semantic_id == constructor)
        .expect("registered persistent parameter constructor");
    let parameter = entry
        .callable
        .as_ref()
        .expect("constructor is callable")
        .parameters
        .iter()
        .find(|parameter| parameter.name == name)
        .expect("registered parameter");
    match parameter.default {
        BindingDefault::Value(value) => value,
        BindingDefault::Required => panic!("persistent parameter must declare its default"),
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct MmffForceFieldParams {
    conformer_id: Option<usize>,
    mmff_variant: String,
    non_bonded_threshold: f64,
    ignore_interfragment_interactions: bool,
}
impl MmffForceFieldParams {
    pub fn new(
        conformer_id: Option<usize>,
        mmff_variant: String,
        non_bonded_threshold: f64,
        ignore_interfragment_interactions: bool,
    ) -> Self {
        Self {
            conformer_id,
            mmff_variant,
            non_bonded_threshold,
            ignore_interfragment_interactions,
        }
    }
    pub fn conformer_id(&self) -> Option<usize> {
        self.conformer_id
    }
    pub fn mmff_variant(&self) -> &str {
        &self.mmff_variant
    }
    pub fn non_bonded_threshold(&self) -> f64 {
        self.non_bonded_threshold
    }
    pub fn ignore_interfragment_interactions(&self) -> bool {
        self.ignore_interfragment_interactions
    }
}
impl Default for MmffForceFieldParams {
    fn default() -> Self {
        static DEFAULT: std::sync::OnceLock<MmffForceFieldParams> = std::sync::OnceLock::new();
        DEFAULT
            .get_or_init(|| {
                assert_eq!(
                    contract_default("MmffForceFieldParams.new", "conformer_id"),
                    "none"
                );
                Self::new(
                    None,
                    contract_default("MmffForceFieldParams.new", "mmff_variant")
                        .trim_matches('"')
                        .into(),
                    contract_default("MmffForceFieldParams.new", "non_bonded_threshold")
                        .parse()
                        .expect("registered numeric default"),
                    contract_default(
                        "MmffForceFieldParams.new",
                        "ignore_interfragment_interactions",
                    )
                    .parse()
                    .expect("registered boolean default"),
                )
            })
            .clone()
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UffForceFieldParams {
    conformer_id: Option<usize>,
    vdw_threshold: f64,
    ignore_interfragment_interactions: bool,
}
impl UffForceFieldParams {
    pub fn new(
        conformer_id: Option<usize>,
        vdw_threshold: f64,
        ignore_interfragment_interactions: bool,
    ) -> Self {
        Self {
            conformer_id,
            vdw_threshold,
            ignore_interfragment_interactions,
        }
    }
    pub fn conformer_id(&self) -> Option<usize> {
        self.conformer_id
    }
    pub fn vdw_threshold(&self) -> f64 {
        self.vdw_threshold
    }
    pub fn ignore_interfragment_interactions(&self) -> bool {
        self.ignore_interfragment_interactions
    }
}
impl Default for UffForceFieldParams {
    fn default() -> Self {
        static DEFAULT: std::sync::OnceLock<UffForceFieldParams> = std::sync::OnceLock::new();
        *DEFAULT.get_or_init(|| {
            assert_eq!(
                contract_default("UffForceFieldParams.new", "conformer_id"),
                "none"
            );
            Self::new(
                None,
                contract_default("UffForceFieldParams.new", "vdw_threshold")
                    .parse()
                    .expect("registered numeric default"),
                contract_default(
                    "UffForceFieldParams.new",
                    "ignore_interfragment_interactions",
                )
                .parse()
                .expect("registered boolean default"),
            )
        })
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ForceFieldMinimizeParams {
    max_iterations: u32,
    force_tolerance: f64,
    energy_tolerance: f64,
}
impl ForceFieldMinimizeParams {
    pub fn new(max_iterations: u32, force_tolerance: f64, energy_tolerance: f64) -> Self {
        Self {
            max_iterations,
            force_tolerance,
            energy_tolerance,
        }
    }
    pub fn max_iterations(&self) -> u32 {
        self.max_iterations
    }
    pub fn force_tolerance(&self) -> f64 {
        self.force_tolerance
    }
    /// Forwarded to RDKit's optimizer, which currently ignores this criterion.
    pub fn energy_tolerance(&self) -> f64 {
        self.energy_tolerance
    }
}
impl Default for ForceFieldMinimizeParams {
    fn default() -> Self {
        static DEFAULT: std::sync::OnceLock<ForceFieldMinimizeParams> = std::sync::OnceLock::new();
        *DEFAULT.get_or_init(|| {
            Self::new(
                contract_default("ForceFieldMinimizeParams.new", "max_iterations")
                    .parse()
                    .expect("registered count default"),
                contract_default("ForceFieldMinimizeParams.new", "force_tolerance")
                    .parse()
                    .expect("registered tolerance default"),
                contract_default("ForceFieldMinimizeParams.new", "energy_tolerance")
                    .parse()
                    .expect("registered energy tolerance default"),
            )
        })
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct ForceFieldEnergyGradient {
    energy: f64,
    gradient: Vec<[f64; 3]>,
}
impl ForceFieldEnergyGradient {
    pub fn energy(&self) -> f64 {
        self.energy
    }
    pub fn gradient(&self) -> &[[f64; 3]] {
        &self.gradient
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ForceFieldMinimizeOutcome {
    converged: bool,
    iterations: u32,
    energy: f64,
}
impl ForceFieldMinimizeOutcome {
    pub fn converged(&self) -> bool {
        self.converged
    }
    pub fn iterations(&self) -> u32 {
        self.iterations
    }
    pub fn energy(&self) -> f64 {
        self.energy
    }
}
#[derive(Debug)]
pub struct MmffForceFieldError(ForceFieldError);
impl MmffForceFieldError {
    pub fn kind(&self) -> MolecularForceFieldErrorKind {
        self.0.kind()
    }
    pub fn requested(&self) -> Option<usize> {
        self.0.requested()
    }
    pub fn atom_index(&self) -> Option<usize> {
        self.0.atom_index()
    }
    pub fn component(&self) -> Option<usize> {
        self.0.component()
    }
    pub fn actual(&self) -> Option<usize> {
        self.0.actual()
    }
    pub fn expected(&self) -> Option<usize> {
        self.0.expected()
    }
}
impl fmt::Display for MmffForceFieldError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}
impl Error for MmffForceFieldError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(&self.0)
    }
}
#[derive(Debug)]
pub struct UffForceFieldError {
    kind: MolecularForceFieldErrorKind,
    requested: Option<usize>,
    source: Option<Box<dyn Error + Send + Sync>>,
}
impl UffForceFieldError {
    fn cause(
        kind: MolecularForceFieldErrorKind,
        source: impl Error + Send + Sync + 'static,
    ) -> Self {
        Self {
            kind,
            requested: None,
            source: Some(Box::new(source)),
        }
    }
    pub fn kind(&self) -> MolecularForceFieldErrorKind {
        self.kind
    }
    pub fn requested(&self) -> Option<usize> {
        self.requested
    }
}
impl fmt::Display for UffForceFieldError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.source {
            Some(source) => source.fmt(f),
            None => write!(
                f,
                "no stored 3D conformer for selector {:?}",
                self.requested
            ),
        }
    }
}
impl Error for UffForceFieldError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        self.source
            .as_deref()
            .map(|source| source as &(dyn Error + 'static))
    }
}

/// An owned interactive evaluator, independent of the source molecule.
///
/// Terms (including the nonbonded pair selection) are built once. Coordinate
/// updates invalidate distances without rebuilding terms or changing the source
/// molecule. Returned coordinates and gradients are independent snapshots.
#[derive(Debug)]
pub struct MolecularForceField {
    inner: Mutex<cosmolkit_forcefields::PreparedForceField>,
}
impl MolecularForceField {
    fn new(inner: cosmolkit_forcefields::PreparedForceField) -> Self {
        Self {
            inner: Mutex::new(inner),
        }
    }
    pub fn position(&self, atom_id: AtomId) -> Result<[f64; 3], ForceFieldError> {
        self.inner
            .lock()
            .expect("force-field state lock")
            .position(atom_id.index())
    }
    pub fn set_position_(
        &mut self,
        atom_id: AtomId,
        position: [f64; 3],
    ) -> Result<(), ForceFieldError> {
        self.inner
            .get_mut()
            .expect("force-field state lock")
            .set_position(atom_id.index(), position)
    }
    pub fn positions(&self) -> Vec<[f64; 3]> {
        self.inner
            .lock()
            .expect("force-field state lock")
            .positions()
    }
    pub fn set_positions_(&mut self, positions: &[[f64; 3]]) -> Result<(), ForceFieldError> {
        self.inner
            .get_mut()
            .expect("force-field state lock")
            .set_positions(positions)
    }
    pub fn fixed_atoms(&self) -> Vec<AtomId> {
        self.inner
            .lock()
            .expect("force-field state lock")
            .fixed_atoms()
            .into_iter()
            .map(AtomId::new)
            .collect()
    }
    pub fn set_fixed_atoms_(&mut self, atom_ids: &[AtomId]) -> Result<(), ForceFieldError> {
        self.inner
            .get_mut()
            .expect("force-field state lock")
            .set_fixed_atoms(&atom_ids.iter().map(|id| id.index()).collect::<Vec<_>>())
    }
    pub fn energy(&self) -> Result<f64, ForceFieldError> {
        self.inner.lock().expect("force-field state lock").energy()
    }
    pub fn gradient(&self) -> Result<Vec<[f64; 3]>, ForceFieldError> {
        self.inner
            .lock()
            .expect("force-field state lock")
            .gradient()
    }
    pub fn energy_gradient(&self) -> Result<ForceFieldEnergyGradient, ForceFieldError> {
        let (energy, gradient) = self
            .inner
            .lock()
            .expect("force-field state lock")
            .energy_gradient()?;
        Ok(ForceFieldEnergyGradient { energy, gradient })
    }
    pub fn minimize_(&mut self) -> Result<ForceFieldMinimizeOutcome, ForceFieldError> {
        self.minimize_with_params_(&ForceFieldMinimizeParams::default())
    }
    /// Relax the current owned coordinates, starting a new optimizer history.
    pub fn minimize_with_params_(
        &mut self,
        params: &ForceFieldMinimizeParams,
    ) -> Result<ForceFieldMinimizeOutcome, ForceFieldError> {
        let (converged, iterations, energy) = self
            .inner
            .get_mut()
            .expect("force-field state lock")
            .minimize(
                params.max_iterations,
                params.force_tolerance,
                params.energy_tolerance,
            )?;
        Ok(ForceFieldMinimizeOutcome {
            converged,
            iterations,
            energy,
        })
    }
}
impl Molecule {
    /// Build an owned MMFF evaluator from the first stored 3D conformer.
    /// Does not generate coordinates, add hydrogens, or modify this molecule.
    pub fn mmff_force_field(&self) -> Result<MolecularForceField, MmffForceFieldError> {
        self.mmff_force_field_with_params(&MmffForceFieldParams::default())
    }
    /// Build MMFF once from the selected stored 3D conformer ID.
    /// A missing ID is an error; `None` selects the first stored 3D conformer.
    pub fn mmff_force_field_with_params(
        &self,
        params: &MmffForceFieldParams,
    ) -> Result<MolecularForceField, MmffForceFieldError> {
        let options = cosmolkit_forcefields::MmffEvaluationParams {
            mmff_variant: params.mmff_variant.clone(),
            conformer_id: params.conformer_id,
            non_bonded_threshold: params.non_bonded_threshold,
            ignore_interfragment_interactions: params.ignore_interfragment_interactions,
        };
        cosmolkit_forcefields::prepare_mmff_force_field(
            self.topology(),
            self.coordinate_block_runtime(),
            self.properties(),
            &options,
            self.derived_cache_runtime().valid_ring_info(),
        )
        .map(MolecularForceField::new)
        .map_err(MmffForceFieldError)
    }
    /// Build an owned UFF evaluator from the first stored 3D conformer.
    /// Requires prepared valence and does not generate or commit coordinates.
    pub fn uff_force_field(&self) -> Result<MolecularForceField, UffForceFieldError> {
        self.uff_force_field_with_params(&UffForceFieldParams::default())
    }
    /// Build UFF once from the selected stored 3D conformer ID.
    /// A missing ID is an error; `None` selects the first stored 3D conformer.
    pub fn uff_force_field_with_params(
        &self,
        params: &UffForceFieldParams,
    ) -> Result<MolecularForceField, UffForceFieldError> {
        let coordinates = self.coordinate_block_runtime();
        let row = match params.conformer_id {
            None => coordinates.conformers_3d.first(),
            Some(id) => coordinates.conformers_3d.iter().find(|row| row.id() == id),
        };
        if row.is_none() {
            return Err(UffForceFieldError {
                kind: MolecularForceFieldErrorKind::MissingConformer,
                requested: params.conformer_id,
                source: None,
            });
        }
        let topology = self.topology();
        let cache = self.derived_cache_runtime();
        cache
            .validate_for_atom_count(topology.atoms.len())
            .map_err(|source| {
                UffForceFieldError::cause(MolecularForceFieldErrorKind::Preparation, source)
            })?;
        let valence = cache.valence_assignment().ok_or_else(|| {
            UffForceFieldError::cause(
                MolecularForceFieldErrorKind::Preparation,
                crate::OperationError::InvalidDerivedCache {
                    state: "valence",
                    field: "assignment",
                    actual: 0,
                    expected: 1,
                },
            )
        })?;
        let computed_rings;
        let rings = match cache.valid_ring_info() {
            Some(rings) => rings,
            None => {
                computed_rings = cosmolkit_core::symmetrized_sssr(
                    topology,
                    &cosmolkit_core::RingSearchParams::default(),
                )
                .map_err(|source| {
                    UffForceFieldError::cause(MolecularForceFieldErrorKind::Rings, source)
                })?;
                &computed_rings
            }
        };
        let options = cosmolkit_forcefields::UffEvaluationParams {
            conformer_id: params.conformer_id,
            vdw_threshold: params.vdw_threshold,
            ignore_interfragment_interactions: params.ignore_interfragment_interactions,
        };
        cosmolkit_forcefields::prepare_uff_force_field(
            topology,
            coordinates,
            valence,
            rings,
            self.properties(),
            &options,
        )
        .map(MolecularForceField::new)
        .map_err(|source| UffForceFieldError::cause(source.kind(), source))
    }
}
