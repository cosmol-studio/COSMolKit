//! Shared detached fingerprint-generator arguments and execution state.
//!
//! Source attribution: pinned RDKit `GraphMol/Fingerprints/
//! FingerprintGenerator.{cpp,h}`, `MorganGenerator.{cpp,h}`, and
//! `MorganFingerprints.cpp`; the Boost.Random RNG closure is pinned and
//! noticed in `THIRD_PARTY_NOTICES.md`.

use crate::Fingerprint;
use crate::FingerprintError;
use crate::MorganError;
use crate::additional_output::AdditionalOutput;
use crate::morgan::{
    MorganAtomEnvironment, MorganDistanceMatrixCache, MorganGenerator, MorganOutput,
    generate_morgan_environments,
};
use crate::prepared::prepare_morgan_environment;
use crate::rng::{BoostFingerprintRng, BoostUniformIntDistribution};
use crate::sparse_bits::SparseBitFingerprint;
use crate::sparse_counts::{SparseCountFingerprint, SparseCountFingerprint32};
use cosmolkit_core::{RingInfo, ValenceAssignment};
use cosmolkit_model::{MoleculeProperties, TopologyBlock};

/// Common source arguments shared by fingerprint generators.
///
/// This is internal construction state; public detached Morgan parameters are
/// kept separate from the generic source arguments.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct FingerprintArguments {
    pub(crate) count_simulation: bool,
    pub(crate) include_chirality: bool,
    pub(crate) count_bounds: Vec<u32>,
    pub(crate) fp_size: u32,
    pub(crate) bits_per_feature: u32,
}

impl Default for FingerprintArguments {
    fn default() -> Self {
        // RDKit✔️✔️: bool df_countSimulation = false;
        // RDKit✔️✔️: bool df_includeChirality = false;
        // RDKit✔️✔️: std::vector<std::uint32_t> d_countBounds;
        // RDKit✔️✔️: std::uint32_t d_fpSize = 2048;
        // RDKit✔️✔️: std::uint32_t d_numBitsPerFeature = 1;
        // RDKit✔️✔️: FingerprintArguments() = default;
        Self {
            count_simulation: false,
            include_chirality: false,
            count_bounds: Vec::new(),
            fp_size: 2048,
            bits_per_feature: 1,
        }
    }
}

impl FingerprintArguments {
    /// Construct with the two trailing defaults from the source declaration.
    pub(crate) fn new_with_source_default_options(
        count_simulation: bool,
        count_bounds: Vec<u32>,
        fp_size: u32,
    ) -> Result<Self, FingerprintError> {
        // RDKit✔️✔️: std::uint32_t numBitsPerFeature = 1,
        // RDKit✔️✔️: bool includeChirality = false);
        Self::new(count_simulation, count_bounds, fp_size, 1, false)
    }

    /// Construct arguments using the source's only two precondition checks.
    pub(crate) fn new(
        count_simulation: bool,
        count_bounds: Vec<u32>,
        fp_size: u32,
        bits_per_feature: u32,
        include_chirality: bool,
    ) -> Result<Self, FingerprintError> {
        // RDKit✔️🔝: FingerprintArguments::FingerprintArguments(
        // RDKit✔️🔝:     const bool countSimulation, const std::vector<std::uint32_t> countBounds,
        // RDKit✔️🔝:     std::uint32_t fpSize, std::uint32_t numBitsPerFeature,
        // RDKit✔️🔝:     bool includeChirality)
        // RDKit✔️🔝:     : df_countSimulation(countSimulation),
        // RDKit✔️🔝:       df_includeChirality(includeChirality),
        // RDKit✔️🔝:       d_countBounds(countBounds),
        // RDKit✔️🔝:       d_fpSize(fpSize),
        // RDKit✔️🔝:       d_numBitsPerFeature(numBitsPerFeature) {
        // RDKit✔️🔝:   PRECONDITION(!countSimulation || !countBounds.empty(),
        // RDKit✔️🔝:                "bad count bounds provided");
        // RDKit✔️🔝:   PRECONDITION(d_numBitsPerFeature > 0, "numBitsPerFeature must be >0");
        // RDKit✔️🔝: }
        // Moving the owned bounds preserves order and duplicates while
        // avoiding the source constructor's by-value vector copies.
        if count_simulation && count_bounds.is_empty() {
            return Err(FingerprintError::PreconditionViolation {
                what: "bad count bounds provided",
            });
        }
        if bits_per_feature == 0 {
            return Err(FingerprintError::PreconditionViolation {
                what: "numBitsPerFeature must be >0",
            });
        }

        Ok(Self {
            count_simulation,
            include_chirality,
            count_bounds,
            fp_size,
            bits_per_feature,
        })
    }
}

/// Borrowed per-call inputs shared by the detached Morgan entrypoints.
///
/// `AdditionalOutput` is a separate mutable output argument in the frozen
/// detached interface, rather than another field on this input bundle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct FingerprintFuncArguments<'a> {
    pub(crate) from_atoms: Option<&'a [u32]>,
    pub(crate) ignore_atoms: Option<&'a [u32]>,
    pub(crate) custom_atom_invariants: Option<&'a [u32]>,
    pub(crate) custom_bond_invariants: Option<&'a [u32]>,
    pub(crate) conformer_id: i32,
}

impl<'a> Default for FingerprintFuncArguments<'a> {
    fn default() -> Self {
        // RDKit✔️✔️: const std::vector<std::uint32_t> *fromAtoms = nullptr;
        // RDKit✔️✔️: const std::vector<std::uint32_t> *ignoreAtoms = nullptr;
        // RDKit✔️✔️: int confId = -1;
        // RDKit✔️✔️: const std::vector<std::uint32_t> *customAtomInvariants = nullptr;
        // RDKit✔️✔️: const std::vector<std::uint32_t> *customBondInvariants = nullptr;
        // The source's AdditionalOutput pointer is represented by the
        // detached API's separate `output: Option<&mut AdditionalOutput>`.
        Self {
            from_atoms: None,
            ignore_atoms: None,
            custom_atom_invariants: None,
            custom_bond_invariants: None,
            conformer_id: -1,
        }
    }
}

impl<'a> FingerprintFuncArguments<'a> {
    /// Store borrowed inputs as supplied, retaining `None` versus empty.
    pub(crate) fn new(
        from_atoms: Option<&'a [u32]>,
        ignore_atoms: Option<&'a [u32]>,
        custom_atom_invariants: Option<&'a [u32]>,
        custom_bond_invariants: Option<&'a [u32]>,
        conformer_id: i32,
    ) -> Self {
        // RDKit✔️✔️: : fromAtoms(fromAtoms_arg),
        // RDKit✔️✔️:   ignoreAtoms(ignoreAtoms_arg),
        // RDKit✔️✔️:   confId(confId_arg),
        // RDKit✔️✔️:   customAtomInvariants(customAtomInvariants_arg),
        // RDKit✔️✔️:   customBondInvariants(customBondInvariants_arg) {};
        // This O(1) holder preserves each source pointer's null/present
        // branch without copying any vector at the argument boundary.
        Self {
            from_atoms,
            ignore_atoms,
            custom_atom_invariants,
            custom_bond_invariants,
            conformer_id,
        }
    }
}

/// Select the source environment topology and source invariant vectors, then
/// invoke the environment stage while a conditional prepared copy is alive.
///
/// The consumer cannot return values borrowing the temporary prepared record;
/// source Morgan environments are consumed before RDKit's temporary `tmol`
/// leaves `getFingerprintHelper` for the same reason.
pub(crate) fn with_morgan_environment_inputs<Output, Consumer>(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    generator: &MorganGenerator,
    arguments: &FingerprintFuncArguments<'_>,
    consume_environments: Consumer,
) -> Result<Output, MorganError>
where
    Consumer: for<'stage> FnOnce(
        &'stage TopologyBlock,
        &'stage MoleculeProperties,
        &'stage [u32],
        &'stage [u32],
    ) -> Result<Output, MorganError>,
{
    with_morgan_environment_inputs_and_output(
        topology,
        properties,
        valence,
        rings,
        generator,
        arguments,
        None,
        |prepared_topology, prepared_properties, atom_invariants, bond_invariants, _output| {
            consume_environments(
                prepared_topology,
                prepared_properties,
                atom_invariants,
                bond_invariants,
            )
        },
    )
}

/// Apply the source's prepared-state, AdditionalOutput reset, and original-
/// input invariant stages before one source-ordered Morgan environment pass.
pub(super) fn with_morgan_environment_inputs_and_output<'output, Output, Consumer>(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    generator: &MorganGenerator,
    arguments: &FingerprintFuncArguments<'_>,
    additional_output: Option<&'output mut AdditionalOutput>,
    consume_environments: Consumer,
) -> Result<Output, MorganError>
where
    Consumer: for<'stage> FnOnce(
        &'stage TopologyBlock,
        &'stage MoleculeProperties,
        &'stage [u32],
        &'stage [u32],
        Option<&'output mut AdditionalOutput>,
    ) -> Result<Output, MorganError>,
{
    with_morgan_environment_inputs_and_output_with_atom_invariants(
        topology,
        properties,
        valence,
        rings,
        generator,
        arguments,
        additional_output,
        source_default_atom_invariants,
        consume_environments,
    )
}

/// Source-ordered preparation with an owner-selected atom-invariant generator.
/// The source custom slice still takes precedence and avoids calling the
/// selected generator entirely.
pub(super) fn with_morgan_environment_inputs_and_output_with_atom_invariants<
    'output,
    Output,
    AtomInvariantProvider,
    Consumer,
>(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    generator: &MorganGenerator,
    arguments: &FingerprintFuncArguments<'_>,
    additional_output: Option<&'output mut AdditionalOutput>,
    atom_invariant_provider: AtomInvariantProvider,
    consume_environments: Consumer,
) -> Result<Output, MorganError>
where
    AtomInvariantProvider: FnOnce(
        &TopologyBlock,
        &MoleculeProperties,
        &ValenceAssignment,
        &RingInfo,
        &MorganGenerator,
    ) -> Result<Vec<u32>, MorganError>,
    Consumer: for<'stage> FnOnce(
        &'stage TopologyBlock,
        &'stage MoleculeProperties,
        &'stage [u32],
        &'stage [u32],
        Option<&'output mut AdditionalOutput>,
    ) -> Result<Output, MorganError>,
{
    with_fingerprint_environment_inputs_and_output(
        topology,
        properties,
        valence,
        rings,
        &generator.fingerprint_arguments,
        arguments,
        additional_output,
        || atom_invariant_provider(topology, properties, valence, rings, generator),
        || Ok(generator.bond_invariants.get_bond_invariants(topology)),
        consume_environments,
    )
}

/// One source preparation/reset/custom-invariant composition shared by families.
pub(super) fn with_fingerprint_environment_inputs_and_output<
    'output,
    Output,
    Error,
    AtomInvariantProvider,
    BondInvariantProvider,
    Consumer,
>(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    common: &FingerprintArguments,
    arguments: &FingerprintFuncArguments<'_>,
    mut additional_output: Option<&'output mut AdditionalOutput>,
    atom_invariant_provider: AtomInvariantProvider,
    bond_invariant_provider: BondInvariantProvider,
    consume_environments: Consumer,
) -> Result<Output, Error>
where
    Error: From<MorganError>,
    AtomInvariantProvider: FnOnce() -> Result<Vec<u32>, Error>,
    BondInvariantProvider: FnOnce() -> Result<Vec<u32>, Error>,
    Consumer: for<'stage> FnOnce(
        &'stage TopologyBlock,
        &'stage MoleculeProperties,
        &'stage [u32],
        &'stage [u32],
        Option<&'output mut AdditionalOutput>,
    ) -> Result<Output, Error>,
{
    // BEGIN RDKIT CPP FUNCTION FingerprintGenerator::getFingerprintHelper original/prepared composition
    // RDKit❗✔️:   const ROMol *lmol = &mol;
    // RDKit❗✔️:   std::unique_ptr<ROMol> tmol;
    // RDKit❗✔️:   if (dp_fingerprintArguments->df_includeChirality &&
    // RDKit❗✔️:       !mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit❗✔️:     tmol = std::unique_ptr<ROMol>(new ROMol(mol));
    // RDKit❗✔️:     MolOps::assignStereochemistry(*tmol);
    // RDKit❗✔️:     lmol = tmol.get();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (args.additionalOutput) {
    // RDKit❗✔️:     reinitAdditionalOutput(*args.additionalOutput, mol.getNumAtoms());
    // RDKit❗✔️:   }
    // RDKit❗✔️:   std::unique_ptr<std::vector<std::uint32_t>> atomInvariants = nullptr;
    // RDKit❗✔️:   if (args.customAtomInvariants) {
    // RDKit❗✔️:     atomInvariants.reset(
    // RDKit❗✔️:         new std::vector<std::uint32_t>(*args.customAtomInvariants));
    // RDKit❗✔️:   } else if (dp_atomInvariantsGenerator) {
    // RDKit❗✔️:     atomInvariants.reset(dp_atomInvariantsGenerator->getAtomInvariants(mol));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   std::unique_ptr<std::vector<std::uint32_t>> bondInvariants = nullptr;
    // RDKit❗✔️:   if (args.customBondInvariants) {
    // RDKit❗✔️:     bondInvariants.reset(
    // RDKit❗✔️:         new std::vector<std::uint32_t>(*args.customBondInvariants));
    // RDKit❗✔️:   } else if (dp_bondInvariantsGenerator) {
    // RDKit❗✔️:     bondInvariants.reset(dp_bondInvariantsGenerator->getBondInvariants(mol));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   auto atomEnvironments = dp_atomEnvironmentGenerator->getEnvironments(
    // RDKit❗✔️:       *lmol, dp_fingerprintArguments, args.fromAtoms, args.ignoreAtoms,
    // RDKit❗✔️:       args.confId, args.additionalOutput, atomInvariants.get(),
    // RDKit❗✔️:       bondInvariants.get(), hashResults);
    // END RDKIT CPP FUNCTION FingerprintGenerator::getFingerprintHelper original/prepared composition
    // Behavior review: source preparation is selected first; present output
    // is reinitialized against the ORIGINAL atom count before custom-vector
    // copies or default invariant generation. Each absent vector is generated
    // from the ORIGINAL topology and caller-owned valence/rings. Only the
    // prepared topology/properties and complete vectors reach environment
    // generation. Each typed family consumer retains atom-then-bond length
    // precondition order and its fixed legacy profile; no additional chemistry path,
    // cache recomputation, zero-vector fallback, or original-input mutation
    // is introduced. Higher-ranked topology borrows keep a result from
    // outliving the temporary prepared copy.
    // Complexity review: source-compatible custom vectors are copied once;
    // absent vectors use existing O(A) atom and O(B) bond owners. Conditional
    // preparation keeps its borrowed fast path or one topology/properties
    // clone. Reinitialization is one O(A+M) pass over existing output fields;
    // no duplicate scan or copy is introduced.
    let prepared = prepare_morgan_environment(
        topology,
        properties,
        valence,
        rings,
        common.include_chirality,
    )
    .map_err(Error::from)?;

    if let Some(output) = additional_output.as_deref_mut() {
        output.reinitialize(topology.atoms.len());
    }

    let atom_invariants = match arguments.custom_atom_invariants {
        Some(custom) => custom.to_vec(),
        None => atom_invariant_provider()?,
    };
    let bond_invariants = match arguments.custom_bond_invariants {
        Some(custom) => custom.to_vec(),
        None => bond_invariant_provider()?,
    };

    consume_environments(
        prepared.topology(),
        prepared.properties(),
        &atom_invariants,
        &bond_invariants,
        additional_output,
    )
}

fn source_default_atom_invariants(
    topology: &TopologyBlock,
    _properties: &MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    generator: &MorganGenerator,
) -> Result<Vec<u32>, MorganError> {
    generator
        .atom_invariants
        .get_atom_invariants(topology, valence, rings)
}

/// Source entry for `FingerprintGenerator::getSparseCountFingerprint`.
/// Its `getFingerprintHelper` call uses the default `fpSize=0`, so Morgan
/// sparse counts retain raw environment IDs and the u64 index domain even
/// when the configured dense `fp_size` is nonzero.
pub(crate) fn get_sparse_count_fingerprint(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    generator: &MorganGenerator,
    arguments: &FingerprintFuncArguments<'_>,
    output: Option<&mut AdditionalOutput>,
) -> Result<SparseCountFingerprint, MorganError> {
    get_sparse_count_fingerprint_with_atom_invariants(
        topology,
        properties,
        valence,
        rings,
        generator,
        arguments,
        output,
        source_default_atom_invariants,
    )
}

pub(super) fn get_sparse_count_fingerprint_with_atom_invariants<AtomInvariantProvider>(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    generator: &MorganGenerator,
    arguments: &FingerprintFuncArguments<'_>,
    output: Option<&mut AdditionalOutput>,
    atom_invariant_provider: AtomInvariantProvider,
) -> Result<SparseCountFingerprint, MorganError>
where
    AtomInvariantProvider: FnOnce(
        &TopologyBlock,
        &MoleculeProperties,
        &ValenceAssignment,
        &RingInfo,
        &MorganGenerator,
    ) -> Result<Vec<u32>, MorganError>,
{
    // BEGIN RDKIT CPP FUNCTION FingerprintGenerator::getSparseCountFingerprint
    // RDKit❗🔝: std::unique_ptr<SparseIntVect<OutputType>>
    // RDKit❗🔝: FingerprintGenerator<OutputType>::getSparseCountFingerprint(
    // RDKit❗🔝:     const ROMol &mol, FingerprintFuncArguments &args) const {
    // RDKit❗🔝:   return getFingerprintHelper(mol, args);
    // RDKit❗🔝: }
    // END RDKIT CPP FUNCTION FingerprintGenerator::getSparseCountFingerprint
    // Behavior review: preserve the wrapper's default helper size, conditional
    // prepared-state handoff, original-input invariants, source environment
    // order, raw u32-origin IDs, count/error behavior, and optional-output
    // mutation. Count simulation flags do not project this wrapper's result.
    // Complexity review: environment traversal is O(E), sparse updates are
    // O(log K), and the configured generator is reused. The source heap-owns
    // each environment and three RNG wrappers when needed; existing Rust Vec
    // and inline RNG/distribution owners avoid those allocations without
    // adding scans or changing stream state.
    let arguments_copy = *arguments;
    with_morgan_environment_inputs_and_output_with_atom_invariants(
        topology,
        properties,
        valence,
        rings,
        generator,
        &arguments_copy,
        output,
        atom_invariant_provider,
        |prepared_topology, _prepared_properties, atom_invariants, bond_invariants, output| {
            let environments = generate_morgan_environments::<u64>(
                prepared_topology,
                generator,
                &arguments_copy,
                atom_invariants,
                bond_invariants,
            )?;
            accumulate_morgan_sparse_counts(
                environments,
                &generator.fingerprint_arguments,
                atom_invariants,
                bond_invariants,
                0,
                output,
            )
        },
    )
}

/// Source entry for `FingerprintGenerator::getSparseFingerprint`.
/// The u64 Morgan environment domain is clamped to SparseBitVect's u32 length;
/// count simulation first folds environment IDs into its effective size.
pub(crate) fn get_sparse_fingerprint(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    generator: &MorganGenerator,
    arguments: &FingerprintFuncArguments<'_>,
    output: Option<&mut AdditionalOutput>,
) -> Result<SparseBitFingerprint, MorganError> {
    get_sparse_fingerprint_with_atom_invariants(
        topology,
        properties,
        valence,
        rings,
        generator,
        arguments,
        output,
        source_default_atom_invariants,
    )
}

pub(super) fn get_sparse_fingerprint_with_atom_invariants<AtomInvariantProvider>(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    generator: &MorganGenerator,
    arguments: &FingerprintFuncArguments<'_>,
    output: Option<&mut AdditionalOutput>,
    atom_invariant_provider: AtomInvariantProvider,
) -> Result<SparseBitFingerprint, MorganError>
where
    AtomInvariantProvider: FnOnce(
        &TopologyBlock,
        &MoleculeProperties,
        &ValenceAssignment,
        &RingInfo,
        &MorganGenerator,
    ) -> Result<Vec<u32>, MorganError>,
{
    let arguments_copy = *arguments;
    project_sparse_fingerprint(
        &generator.fingerprint_arguments,
        u64::MAX,
        topology.atoms.len(),
        output,
        |fp_size, helper_output| {
            with_morgan_environment_inputs_and_output_with_atom_invariants(
                topology,
                properties,
                valence,
                rings,
                generator,
                &arguments_copy,
                helper_output,
                atom_invariant_provider,
                |prepared_topology,
                 _prepared_properties,
                 atom_invariants,
                 bond_invariants,
                 output| {
                    let environments = generate_morgan_environments::<u64>(
                        prepared_topology,
                        generator,
                        &arguments_copy,
                        atom_invariants,
                        bond_invariants,
                    )?;
                    accumulate_morgan_sparse_counts(
                        environments,
                        &generator.fingerprint_arguments,
                        atom_invariants,
                        bond_invariants,
                        fp_size,
                        output,
                    )
                },
            )
        },
    )
}

/// Source entry for `FingerprintGenerator::getCountFingerprint`.
/// Unlike the sparse-count variant, this hashes at the configured `fp_size`
/// and returns the source's 32-bit sparse-count value.
pub(crate) fn get_count_fingerprint(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    generator: &MorganGenerator,
    arguments: &FingerprintFuncArguments<'_>,
    output: Option<&mut AdditionalOutput>,
) -> Result<SparseCountFingerprint32, MorganError> {
    get_count_fingerprint_with_atom_invariants(
        topology,
        properties,
        valence,
        rings,
        generator,
        arguments,
        output,
        source_default_atom_invariants,
    )
}

pub(super) fn get_count_fingerprint_with_atom_invariants<AtomInvariantProvider>(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    generator: &MorganGenerator,
    arguments: &FingerprintFuncArguments<'_>,
    output: Option<&mut AdditionalOutput>,
    atom_invariant_provider: AtomInvariantProvider,
) -> Result<SparseCountFingerprint32, MorganError>
where
    AtomInvariantProvider: FnOnce(
        &TopologyBlock,
        &MoleculeProperties,
        &ValenceAssignment,
        &RingInfo,
        &MorganGenerator,
    ) -> Result<Vec<u32>, MorganError>,
{
    let arguments_copy = *arguments;
    project_count_fingerprint(
        &generator.fingerprint_arguments,
        topology.atoms.len(),
        output,
        |fp_size, helper_output| {
            with_morgan_environment_inputs_and_output_with_atom_invariants(
                topology,
                properties,
                valence,
                rings,
                generator,
                &arguments_copy,
                helper_output,
                atom_invariant_provider,
                |prepared_topology,
                 _prepared_properties,
                 atom_invariants,
                 bond_invariants,
                 output| {
                    let environments = generate_morgan_environments::<u64>(
                        prepared_topology,
                        generator,
                        &arguments_copy,
                        atom_invariants,
                        bond_invariants,
                    )?;
                    accumulate_morgan_sparse_counts(
                        environments,
                        &generator.fingerprint_arguments,
                        atom_invariants,
                        bond_invariants,
                        fp_size,
                        output,
                    )
                },
            )
        },
    )
}

/// Source entry for `FingerprintGenerator::getFingerprint`.
/// The configured size is the dense result length and, when count simulation
/// is enabled, also determines the environment-folding size.
pub(crate) fn get_fingerprint(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    generator: &MorganGenerator,
    arguments: &FingerprintFuncArguments<'_>,
    output: Option<&mut AdditionalOutput>,
) -> Result<Fingerprint, MorganError> {
    get_fingerprint_with_atom_invariants(
        topology,
        properties,
        valence,
        rings,
        generator,
        arguments,
        output,
        source_default_atom_invariants,
    )
}

pub(super) fn get_fingerprint_with_atom_invariants<AtomInvariantProvider>(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    generator: &MorganGenerator,
    arguments: &FingerprintFuncArguments<'_>,
    output: Option<&mut AdditionalOutput>,
    atom_invariant_provider: AtomInvariantProvider,
) -> Result<Fingerprint, MorganError>
where
    AtomInvariantProvider: FnOnce(
        &TopologyBlock,
        &MoleculeProperties,
        &ValenceAssignment,
        &RingInfo,
        &MorganGenerator,
    ) -> Result<Vec<u32>, MorganError>,
{
    let arguments_copy = *arguments;
    project_fingerprint(
        &generator.fingerprint_arguments,
        topology.atoms.len(),
        output,
        |fp_size, helper_output| {
            with_morgan_environment_inputs_and_output_with_atom_invariants(
                topology,
                properties,
                valence,
                rings,
                generator,
                &arguments_copy,
                helper_output,
                atom_invariant_provider,
                |prepared_topology,
                 _prepared_properties,
                 atom_invariants,
                 bond_invariants,
                 output| {
                    let environments = generate_morgan_environments::<u64>(
                        prepared_topology,
                        generator,
                        &arguments_copy,
                        atom_invariants,
                        bond_invariants,
                    )?;
                    accumulate_morgan_sparse_counts(
                        environments,
                        &generator.fingerprint_arguments,
                        atom_invariants,
                        bond_invariants,
                        fp_size,
                        output,
                    )
                },
            )
        },
    )
}

/// Shared source projection; chemistry is supplied by the family owner.
pub(super) fn project_sparse_fingerprint<Error, Consumer>(
    fingerprint_arguments: &FingerprintArguments,
    result_size: u64,
    atom_count: usize,
    output: Option<&mut AdditionalOutput>,
    consumer: Consumer,
) -> Result<SparseBitFingerprint, Error>
where
    Error: From<FingerprintError>,
    Consumer: FnOnce(u64, Option<&mut AdditionalOutput>) -> Result<SparseCountFingerprint, Error>,
{
    // BEGIN RDKIT CPP FUNCTION FingerprintGenerator::getSparseFingerprint
    // RDKit❗✔️: template <typename OutputType>
    // RDKit❗✔️: std::unique_ptr<SparseBitVect>
    // RDKit❗✔️: FingerprintGenerator<OutputType>::getSparseFingerprint(
    // RDKit❗✔️:     const ROMol &mol, FingerprintFuncArguments &args) const {
    // RDKit❗✔️:   // make sure the result will fit into SparseBitVect
    // RDKit❗✔️:   std::uint32_t resultSize =
    // RDKit❗✔️:       std::min((std::uint64_t)std::numeric_limits<std::uint32_t>::max(),
    // RDKit❗✔️:                (std::uint64_t)dp_atomEnvironmentGenerator->getResultSize());
    // RDKit❗✔️:   std::uint32_t effectiveSize = resultSize;
    // RDKit❗✔️:   if (dp_fingerprintArguments->df_countSimulation) {
    // RDKit❗✔️:     effectiveSize /= dp_fingerprintArguments->d_countBounds.size();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   AdditionalOutput countSimulationOutput;
    // RDKit❗✔️:   AdditionalOutput *origAO = nullptr;
    // RDKit❗✔️:   if (dp_fingerprintArguments->df_countSimulation && args.additionalOutput) {
    // RDKit❗✔️:     setupTempAdditionalOutput(args, countSimulationOutput, mol.getNumAtoms());
    // RDKit❗✔️:     origAO = args.additionalOutput;
    // RDKit❗✔️:     args.additionalOutput = &countSimulationOutput;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   auto tempResult = getFingerprintHelper(mol, args, effectiveSize);
    // RDKit❗✔️:   auto result = std::make_unique<SparseBitVect>(resultSize);
    // RDKit❗✔️:   for (auto val : tempResult->getNonzeroElements()) {
    // RDKit❗✔️:     if (dp_fingerprintArguments->df_countSimulation) {
    // RDKit❗✔️:       for (unsigned int i = 0;
    // RDKit❗✔️:            i < dp_fingerprintArguments->d_countBounds.size(); ++i) {
    // RDKit❗✔️:         const auto &bounds_count = dp_fingerprintArguments->d_countBounds;
    // RDKit❗✔️:         if (val.second >= static_cast<int>(bounds_count[i])) {
    // RDKit❗✔️:           OutputType nBitId = val.first * bounds_count.size() + i;
    // RDKit❗✔️:           result->setBit(nBitId);
    // RDKit❗✔️:           if (args.additionalOutput) {
    // RDKit❗✔️:             duplicateAdditionalOutputBit(*args.additionalOutput, *origAO,
    // RDKit❗✔️:                                          static_cast<OutputType>(val.first), nBitId);
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       result->setBit(val.first);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (origAO) {
    // RDKit❗✔️:     if (origAO->atomCounts) {
    // RDKit❗✔️:       *origAO->atomCounts = *countSimulationOutput.atomCounts;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     args.additionalOutput = origAO;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return result;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION FingerprintGenerator::getSparseFingerprint
    // Behavior review: this fixed u64 Morgan instantiation has source result
    // size min(u32::MAX, u64::MAX) == u32::MAX, independent of configured
    // fpSize. Count simulation divides only when enabled; its validated
    // nonempty bounds keep source order and duplicates. The existing helper
    // folds environment IDs at effective_size before the sparse count map is
    // thresholded. The temporary output exists only for count simulation plus
    // caller output, mirrors all five allocations, and duplicates provenance
    // in source bound order; only atomCounts are copied back after the loop.
    // `bound as i32` preserves the pinned signed comparison, and wrapping u64
    // multiply/add followed by u32 narrowing preserves source unsigned
    // threshold-bit arithmetic and SparseBitVect's index conversion.
    // Complexity review: source and Rust visit K ordered nonzero counts and N
    // bounds, O(KN), with ordered sparse-set insertion. For every qualifying
    // bound, existing AO duplication scans A atom rows and performs ordered-map
    // lookup/copy only for present source keys. Temporary AO allocation and
    // reset remain conditional on count simulation and a present output;
    // existing tree-backed maps and set retain source lookup complexity.
    let count_simulation = fingerprint_arguments.count_simulation;
    let count_bounds = &fingerprint_arguments.count_bounds;
    let result_size = result_size.min(u64::from(u32::MAX)) as u32;
    let effective_size = if count_simulation {
        (result_size as usize / count_bounds.len()) as u32
    } else {
        result_size
    };

    let mut original_output = output;
    let mut count_simulation_output = if count_simulation {
        original_output
            .as_deref_mut()
            .map(|original| original.setup_count_simulation_output(atom_count))
    } else {
        None
    };
    let helper_output = match count_simulation_output.as_mut() {
        Some(temporary) => Some(&mut *temporary),
        None => original_output.as_deref_mut(),
    };

    let sparse_counts = consumer(u64::from(effective_size), helper_output)?;

    let mut result = SparseBitFingerprint::new(result_size);
    for (&base_bit_id, &count) in sparse_counts.nonzero_elements() {
        if count_simulation {
            for (bound_index, &bound) in count_bounds.iter().enumerate() {
                if count >= bound as i32 {
                    let new_bit_id = base_bit_id
                        .wrapping_mul(count_bounds.len() as u64)
                        .wrapping_add(bound_index as u64);
                    result.set_bit(new_bit_id as u32)?;

                    if let Some(temporary_output) = count_simulation_output.as_ref() {
                        let original = original_output
                            .as_deref_mut()
                            .expect("temporary count output requires an original output");
                        temporary_output.duplicate_bit_to(original, base_bit_id, new_bit_id)?;
                    }
                }
            }
        } else {
            result.set_bit(base_bit_id as u32)?;
        }
    }

    if let Some(temporary_output) = count_simulation_output.as_ref() {
        let original = original_output
            .as_deref_mut()
            .expect("temporary count output requires an original output");
        if let Some(original_atom_counts) = original.atom_counts.as_mut() {
            let temporary_atom_counts = temporary_output
                .atom_counts
                .as_ref()
                .expect("the temporary output mirrors atomCounts allocation");
            original_atom_counts.clone_from(temporary_atom_counts);
        }
    }

    Ok(result)
}

/// Shared source projection; chemistry is supplied by the family owner.
pub(super) fn project_count_fingerprint<Error, Consumer>(
    fingerprint_arguments: &FingerprintArguments,

    atom_count: usize,
    output: Option<&mut AdditionalOutput>,
    consumer: Consumer,
) -> Result<SparseCountFingerprint32, Error>
where
    Error: From<FingerprintError>,
    Consumer: FnOnce(u64, Option<&mut AdditionalOutput>) -> Result<SparseCountFingerprint, Error>,
{
    // BEGIN RDKIT CPP FUNCTION FingerprintGenerator::getCountFingerprint
    // RDKit❗🔝: template <typename OutputType>
    // RDKit❗🔝: std::unique_ptr<SparseIntVect<std::uint32_t>>
    // RDKit❗🔝: FingerprintGenerator<OutputType>::getCountFingerprint(
    // RDKit❗🔝:     const ROMol &mol, FingerprintFuncArguments &args) const {
    // RDKit❗🔝:   auto tempResult =
    // RDKit❗🔝:       getFingerprintHelper(mol, args, dp_fingerprintArguments->d_fpSize);
    // RDKit❗🔝:   auto result = std::make_unique<SparseIntVect<std::uint32_t>>(
    // RDKit❗🔝:       dp_fingerprintArguments->d_fpSize);
    // RDKit❗🔝:   for (auto val : tempResult->getNonzeroElements()) {
    // RDKit❗🔝:     result->setVal(val.first, val.second);
    // RDKit❗🔝:   }
    // RDKit❗🔝:   return result;
    // RDKit❗🔝: }
    // END RDKIT CPP FUNCTION FingerprintGenerator::getCountFingerprint
    // Behavior review: pass the configured size into the shared helper so
    // every environment and additional bit is folded before collisions are
    // counted. The source's final u64-to-u32 index narrowing is explicit in
    // `as u32`, and the existing u32 sparse value performs the source index
    // check. Thus fp_size=0 returns empty only when there are no environments;
    // any nonzero count reaches the source IndexError category. The source
    // passes AdditionalOutput directly through the generic helper, with no
    // count-simulation temporary or threshold projection.
    // Complexity review: reuse the one preparation, environment, RNG, and
    // BTreeMap accumulator; then retain the source's second ordered-map copy
    // into the configured u32 sparse value. Sparse work remains O(E log K +
    // K log K), with source-shaped O(K log K) projection and no dense buffer.
    // Existing contiguous environments and inline RNG state remove the
    // source's per-environment and conditional RNG wrapper allocations.
    let fp_size = u64::from(fingerprint_arguments.fp_size);
    let sparse_counts = consumer(fp_size, output)?;

    let mut result = SparseCountFingerprint32::new(fingerprint_arguments.fp_size);
    for (&source_index, &count) in sparse_counts.nonzero_elements() {
        result.set_value(source_index as u32, count)?;
    }

    Ok(result)
}

/// Shared source projection; chemistry is supplied by the family owner.
pub(super) fn project_fingerprint<Error, Consumer>(
    fingerprint_arguments: &FingerprintArguments,

    atom_count: usize,
    output: Option<&mut AdditionalOutput>,
    consumer: Consumer,
) -> Result<Fingerprint, Error>
where
    Error: From<FingerprintError>,
    Consumer: FnOnce(u64, Option<&mut AdditionalOutput>) -> Result<SparseCountFingerprint, Error>,
{
    // BEGIN RDKIT CPP FUNCTION FingerprintGenerator::getFingerprint
    // RDKit❗🔝: template <typename OutputType>
    // RDKit❗🔝: std::unique_ptr<ExplicitBitVect>
    // RDKit❗🔝: FingerprintGenerator<OutputType>::getFingerprint(
    // RDKit❗🔝:     const ROMol &mol, FingerprintFuncArguments &args) const {
    // RDKit❗🔝:   std::uint32_t effectiveSize = dp_fingerprintArguments->d_fpSize;
    // RDKit❗🔝:   if (dp_fingerprintArguments->df_countSimulation) {
    // RDKit❗🔝:     if (dp_fingerprintArguments->d_countBounds.empty()) {
    // RDKit❗🔝:       throw ValueErrorException("Count bounds are empty");
    // RDKit❗🔝:     }
    // RDKit❗🔝:     if (dp_fingerprintArguments->d_countBounds.size() >= effectiveSize) {
    // RDKit❗🔝:       throw ValueErrorException("Count bounds size is >= fingerprint size");
    // RDKit❗🔝:     }
    // RDKit❗🔝:     effectiveSize /= dp_fingerprintArguments->d_countBounds.size();
    // RDKit❗🔝:   }
    // RDKit❗🔝:   AdditionalOutput countSimulationOutput;
    // RDKit❗🔝:   AdditionalOutput *origAO = nullptr;
    // RDKit❗🔝:   if (dp_fingerprintArguments->df_countSimulation && args.additionalOutput) {
    // RDKit❗🔝:     setupTempAdditionalOutput(args, countSimulationOutput, mol.getNumAtoms());
    // RDKit❗🔝:     origAO = args.additionalOutput;
    // RDKit❗🔝:     args.additionalOutput = &countSimulationOutput;
    // RDKit❗🔝:   }
    // RDKit❗🔝:   auto tempResult = getFingerprintHelper(mol, args, effectiveSize);
    // RDKit❗🔝:   auto result = std::make_unique<ExplicitBitVect>(dp_fingerprintArguments->d_fpSize);
    // RDKit❗🔝:   for (auto val : tempResult->getNonzeroElements()) {
    // RDKit❗🔝:     if (dp_fingerprintArguments->df_countSimulation) {
    // RDKit❗🔝:       for (unsigned int i = 0; i < dp_fingerprintArguments->d_countBounds.size(); ++i) {
    // RDKit❗🔝:         const auto &bounds_count = dp_fingerprintArguments->d_countBounds;
    // RDKit❗🔝:         if (val.second >= static_cast<int>(bounds_count[i])) {
    // RDKit❗🔝:           OutputType nBitId = val.first * bounds_count.size() + i;
    // RDKit❗🔝:           result->setBit(nBitId);
    // RDKit❗🔝:           if (args.additionalOutput) {
    // RDKit❗🔝:             duplicateAdditionalOutputBit(*args.additionalOutput, *origAO,
    // RDKit❗🔝:                                          static_cast<OutputType>(val.first), nBitId);
    // RDKit❗🔝:           }
    // RDKit❗🔝:         }
    // RDKit❗🔝:       }
    // RDKit❗🔝:     } else {
    // RDKit❗🔝:       result->setBit(val.first);
    // RDKit❗🔝:     }
    // RDKit❗🔝:   }
    // RDKit❗🔝:   if (origAO) {
    // RDKit❗🔝:     if (origAO->atomCounts) {
    // RDKit❗🔝:       *origAO->atomCounts = *countSimulationOutput.atomCounts;
    // RDKit❗🔝:     }
    // RDKit❗🔝:     args.additionalOutput = origAO;
    // RDKit❗🔝:   }
    // RDKit❗🔝:   return result;
    // RDKit❗🔝: }
    // END RDKIT CPP FUNCTION FingerprintGenerator::getFingerprint
    // Behavior review: count bounds are validated only for count simulation,
    // with empty bounds rejected before `bounds.len() >= fp_size`; the check
    // guarantees `effective_size >= 1` before division. Without simulation,
    // zero fp_size reaches helper hashing-disabled mode and a zero-length
    // dense value; empty environments return it while the first source-ordered
    // nonzero ID fails the dense setter with its typed index error. Dense
    // output has exactly configured length. A temporary AdditionalOutput is
    // created only for simulation plus a present caller output, mirrors all
    // five allocation states, and resets the original's four source-reset
    // fields before the helper resets/updates the temporary. Thresholds retain
    // source order and signed `count >= bound as i32`; folded threshold IDs use
    // source-width wrapping multiplication/addition before dense u32 setter
    // narrowing. Provenance duplication uses the one existing owner; only
    // atomCounts is copied back after successful projection.
    // Complexity review: reuse the existing O(E log K) environment-count
    // accumulator, then visit K ordered nonzero counts and N bounds (O(KN));
    // each accepted threshold performs one checked dense set and the existing
    // AO duplication cost. The dense value is allocated once at configured
    // length and bits are set during that same projection pass, with no
    // intermediate list of on-bit IDs or second scan.
    let count_simulation = fingerprint_arguments.count_simulation;
    let count_bounds = &fingerprint_arguments.count_bounds;
    let mut effective_size = fingerprint_arguments.fp_size;

    if count_simulation {
        if count_bounds.is_empty() {
            return Err(FingerprintError::InvalidArguments {
                reason: "Count bounds are empty",
            }
            .into());
        }
        if count_bounds.len() >= effective_size as usize {
            return Err(FingerprintError::InvalidArguments {
                reason: "Count bounds size is >= fingerprint size",
            }
            .into());
        }
        effective_size /= count_bounds.len() as u32;
    }

    let mut original_output = output;
    let mut count_simulation_output = if count_simulation {
        original_output
            .as_deref_mut()
            .map(|original| original.setup_count_simulation_output(atom_count))
    } else {
        None
    };
    let helper_output = match count_simulation_output.as_mut() {
        Some(temporary) => Some(&mut *temporary),
        None => original_output.as_deref_mut(),
    };

    let sparse_counts = consumer(u64::from(effective_size), helper_output)?;

    let mut result = Fingerprint::new(fingerprint_arguments.fp_size);
    for (&base_bit_id, &count) in sparse_counts.nonzero_elements() {
        if count_simulation {
            for (bound_index, &bound) in count_bounds.iter().enumerate() {
                if count >= bound as i32 {
                    let new_bit_id = base_bit_id
                        .wrapping_mul(count_bounds.len() as u64)
                        .wrapping_add(bound_index as u64);
                    result.set_bit(new_bit_id as u32)?;

                    if let Some(temporary_output) = count_simulation_output.as_ref() {
                        let original = original_output
                            .as_deref_mut()
                            .expect("temporary count output requires an original output");
                        temporary_output.duplicate_bit_to(original, base_bit_id, new_bit_id)?;
                    }
                }
            }
        } else {
            result.set_bit(base_bit_id as u32)?;
        }
    }

    if let Some(temporary_output) = count_simulation_output.as_ref() {
        let original = original_output
            .as_deref_mut()
            .expect("temporary count output requires an original output");
        if let Some(original_atom_counts) = original.atom_counts.as_mut() {
            let temporary_atom_counts = temporary_output
                .atom_counts
                .as_ref()
                .expect("the temporary output mirrors atomCounts allocation");
            original_atom_counts.clone_from(temporary_atom_counts);
        }
    }

    Ok(result)
}

/// Accumulate the source-ordered Morgan environments into the generic sparse
/// count result used by the later value projections.
pub(super) fn accumulate_morgan_sparse_counts<'a>(
    environments: Vec<MorganAtomEnvironment<'a, u64>>,
    fingerprint_arguments: &FingerprintArguments,
    atom_invariants: &[u32],
    bond_invariants: &[u32],
    fp_size: u64,
    mut output: Option<&mut AdditionalOutput>,
) -> Result<SparseCountFingerprint, MorganError> {
    // BEGIN RDKIT CPP FUNCTION FingerprintGenerator::getFingerprintHelper sparse accumulation
    // RDKit❗🔝:   auto res = std::make_unique<SparseIntVect<OutputType>>(
    // RDKit❗🔝:       fpSize ? fpSize : dp_atomEnvironmentGenerator->getResultSize());
    // END RDKIT CPP FUNCTION FingerprintGenerator::getFingerprintHelper sparse accumulation
    // Behavior review: the source helper selects `fpSize` or the environment
    // generator's result width before entering the shared accumulation pass.
    // Complexity review: one empty ordered result is allocated; the width
    // specialization is supplied by the caller without a conversion map.
    let result_length = if fp_size != 0 { fp_size } else { u64::MAX };
    let mut result = SparseCountFingerprint::new(result_length);
    accumulate_morgan_sparse_counts_into(
        environments,
        fingerprint_arguments,
        atom_invariants,
        bond_invariants,
        fp_size,
        output,
        &mut result,
    )?;

    Ok(result)
}

/// One sparse-count accumulation loop shared by both source index widths.
/// The legacy wrapper's u32 specialization writes directly into its final
/// index domain instead of copying a complete u64 map after generation.
pub(super) trait FingerprintEnvironment<Error> {
    type Output: MorganOutput;
    type State: Default;
    fn bit_id(
        &self,
        arguments: &FingerprintArguments,
        atom_invariants: &[u32],
        bond_invariants: &[u32],
        output: Option<&mut AdditionalOutput>,
        hash_results: bool,
        fp_size: u64,
    ) -> Result<Self::Output, Error>;
    fn update_output(
        &self,
        output: &mut AdditionalOutput,
        bit_id: u64,
        state: &mut Self::State,
    ) -> Result<(), Error>;
}

impl<'a, Output: MorganOutput> FingerprintEnvironment<MorganError>
    for MorganAtomEnvironment<'a, Output>
{
    type Output = Output;
    type State = MorganDistanceMatrixCache<'a>;
    fn bit_id(
        &self,
        arguments: &FingerprintArguments,
        atom_invariants: &[u32],
        bond_invariants: &[u32],
        output: Option<&mut AdditionalOutput>,
        hash_results: bool,
        fp_size: u64,
    ) -> Result<Output, MorganError> {
        Ok(self.get_bit_id(
            Some(arguments),
            Some(atom_invariants),
            Some(bond_invariants),
            output,
            hash_results,
            fp_size,
        ))
    }
    fn update_output(
        &self,
        output: &mut AdditionalOutput,
        bit_id: u64,
        state: &mut Self::State,
    ) -> Result<(), MorganError> {
        Ok(self.update_additional_output(output, bit_id, state)?)
    }
}

pub(super) fn accumulate_morgan_sparse_counts_into<'a, OutputType, Store>(
    environments: Vec<MorganAtomEnvironment<'a, OutputType>>,
    fingerprint_arguments: &FingerprintArguments,
    atom_invariants: &[u32],
    bond_invariants: &[u32],
    fp_size: u64,
    mut output: Option<&mut AdditionalOutput>,
    result: &mut Store,
) -> Result<(), MorganError>
where
    OutputType: MorganOutput,
    Store: MorganSparseCountStore<MorganError>,
{
    accumulate_sparse_counts_into(
        environments,
        fingerprint_arguments,
        atom_invariants,
        bond_invariants,
        fp_size,
        output,
        result,
    )
}

/// One source accumulation loop for all families and both source widths.
pub(super) fn accumulate_sparse_counts_into<Environment, Store, Error>(
    environments: Vec<Environment>,
    fingerprint_arguments: &FingerprintArguments,
    atom_invariants: &[u32],
    bond_invariants: &[u32],
    fp_size: u64,
    mut output: Option<&mut AdditionalOutput>,
    result: &mut Store,
) -> Result<(), Error>
where
    Error: From<FingerprintError>,
    Environment: FingerprintEnvironment<Error>,
    Store: MorganSparseCountStore<Error>,
{
    // BEGIN RDKIT CPP FUNCTION FingerprintGenerator::getFingerprintHelper sparse accumulation
    // RDKit❗✔️:   // define a mersenne twister with customized parameters.
    // RDKit❗✔️:   // The standard parameters (used to create boost::mt19937)
    // RDKit❗✔️:   // result in an RNG that's much too computationally intensive
    // RDKit❗✔️:   // to seed.
    // RDKit❗✔️:   // These are the parameters that have been used for the RDKit fingerprint.
    // RDKit❗✔️:   typedef boost::random::mersenne_twister<std::uint32_t, 32, 4, 2, 31,
    // RDKit❗✔️:                                           0x9908b0df, 11, 7, 0x9d2c5680, 15,
    // RDKit❗✔️:                                           0xefc60000, 18, 3346425566U>
    // RDKit❗✔️:       rng_type;
    // RDKit❗✔️:   typedef boost::uniform_int<> distrib_type;
    // RDKit❗✔️:   typedef boost::variate_generator<rng_type &, distrib_type> source_type;
    // RDKit❗✔️:   std::unique_ptr<rng_type> generator;
    // RDKit❗✔️:   //
    // RDKit❗✔️:   // if we generate arbitrarily sized ints then mod them down to the
    // RDKit❗✔️:   // appropriate size, we can guarantee that a fingerprint of
    // RDKit❗✔️:   // size x has the same bits set as one of size 2x that's been folded
    // RDKit❗✔️:   // in half.  This is a nice guarantee to have.
    // RDKit❗✔️:   //
    // RDKit❗✔️:   std::unique_ptr<distrib_type> dist;
    // RDKit❗✔️:   std::unique_ptr<source_type> randomSource;
    // RDKit❗✔️:   if (dp_fingerprintArguments->d_numBitsPerFeature > 1) {
    // RDKit❗✔️:     // we will only create the RNG if we're going to need it
    // RDKit❗✔️:     generator.reset(new rng_type(42u));
    // RDKit❗✔️:     dist.reset(new distrib_type(0, INT_MAX));
    // RDKit❗✔️:     randomSource.reset(new source_type(*generator, *dist));
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // iterate over every atom environment and generate bit-ids that will make
    // RDKit❗✔️:   // up the fingerprint
    // RDKit❗✔️:   for (const auto env : atomEnvironments) {
    // RDKit❗✔️:     OutputType seed = env->getBitId(dp_fingerprintArguments,
    // RDKit❗✔️:                                     atomInvariants.get(), bondInvariants.get(),
    // RDKit❗✔️:                                     args.additionalOutput, hashResults, fpSize);
    // RDKit❗✔️:
    // RDKit❗✔️:     auto bitId = seed;
    // RDKit❗✔️:     if (fpSize != 0) {
    // RDKit❗✔️:       bitId %= fpSize;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     res->setVal(bitId, res->getVal(bitId) + 1);
    // RDKit❗✔️:     if (args.additionalOutput) {
    // RDKit❗✔️:       env->updateAdditionalOutput(args.additionalOutput, bitId);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // do the additional bits if required:
    // RDKit❗✔️:     if (dp_fingerprintArguments->d_numBitsPerFeature > 1) {
    // RDKit❗✔️:       generator->seed(static_cast<rng_type::result_type>(seed));
    // RDKit❗✔️:
    // RDKit❗✔️:       for (boost::uint32_t bitN = 1;
    // RDKit❗✔️:            bitN < dp_fingerprintArguments->d_numBitsPerFeature; ++bitN) {
    // RDKit❗✔️:         bitId = (*randomSource)();
    // RDKit❗✔️:         if (fpSize != 0) {
    // RDKit❗✔️:           bitId %= fpSize;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         res->setVal(bitId, res->getVal(bitId) + 1);
    // RDKit❗✔️:         if (args.additionalOutput) {
    // RDKit❗✔️:           env->updateAdditionalOutput(args.additionalOutput, bitId);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     delete env;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION FingerprintGenerator::getFingerprintHelper sparse accumulation

    // Behavior review: source getBitId returns the environment's source-width
    // code unchanged. The base ID is folded only when fpSize is nonzero; its
    // count is written before optional AdditionalOutput. For extra bits, the
    // source u32 seed conversion/reseed precedes each ordered draw, fold, count
    // and metadata update. The generator/distribution live across the outer
    // environment loop, and the owned Rust environment drops at loop end.
    // Ownership review: C++ allocates `res` at FingerprintGenerator.cpp:369–371
    // and returns it at :435 within one helper. The Rust outer
    // accumulate_morgan_sparse_counts wrapper allocates and returns its result;
    // this `_into` helper mutates the supplied `&mut Store` and returns unit.
    // This keeps allocation/return responsibility with the wrapper and the
    // source-ordered accumulation work with the inner helper.
    // Complexity review: source ordered-map updates and metadata costs remain
    // O(E log K) plus optional output work; no extra map or environment pass
    // is introduced for the u32 legacy specialization.
    let hash_results = fp_size != 0;
    let mut random = (fingerprint_arguments.bits_per_feature > 1).then(|| {
        (
            BoostFingerprintRng::new(42),
            BoostUniformIntDistribution::rdkit_additional_bit(),
        )
    });
    let mut environment_state = Environment::State::default();

    for environment in environments {
        let seed_output = environment.bit_id(
            fingerprint_arguments,
            atom_invariants,
            bond_invariants,
            output.as_deref_mut(),
            hash_results,
            fp_size,
        )?;
        let seed = seed_output.into_u64();

        let mut bit_id = if fp_size != 0 { seed % fp_size } else { seed };
        increment_morgan_sparse_count(result, bit_id)?;
        if let Some(additional_output) = output.as_deref_mut() {
            environment.update_output(additional_output, bit_id, &mut environment_state)?;
        }

        if fingerprint_arguments.bits_per_feature > 1 {
            let (generator, distribution) = random
                .as_mut()
                .expect("the source creates the RNG when bitsPerFeature is greater than one");
            generator.seed(seed_output.into_u32());

            for _bit_number in 1..fingerprint_arguments.bits_per_feature {
                bit_id = u64::from(distribution.sample(generator) as u32);
                if fp_size != 0 {
                    bit_id %= fp_size;
                }
                increment_morgan_sparse_count(result, bit_id)?;
                if let Some(additional_output) = output.as_deref_mut() {
                    environment.update_output(additional_output, bit_id, &mut environment_state)?;
                }
            }
        }
    }
    Ok(())
}

pub(super) trait MorganSparseCountStore<Error> {
    fn value_at(&self, bit_id: u64) -> Result<i32, Error>;
    fn set_at(&mut self, bit_id: u64, value: i32) -> Result<(), Error>;
}

impl<Error: From<FingerprintError>> MorganSparseCountStore<Error> for SparseCountFingerprint {
    fn value_at(&self, bit_id: u64) -> Result<i32, Error> {
        Ok(self.value(bit_id)?)
    }

    fn set_at(&mut self, bit_id: u64, value: i32) -> Result<(), Error> {
        Ok(self.set_value(bit_id, value)?)
    }
}

impl<Error: From<FingerprintError>> MorganSparseCountStore<Error> for SparseCountFingerprint32 {
    fn value_at(&self, bit_id: u64) -> Result<i32, Error> {
        Ok(self.value(bit_id as u32)?)
    }

    fn set_at(&mut self, bit_id: u64, value: i32) -> Result<(), Error> {
        Ok(self.set_value(bit_id as u32, value)?)
    }
}

fn increment_morgan_sparse_count<
    Store: MorganSparseCountStore<Error>,
    Error: From<FingerprintError>,
>(
    result: &mut Store,
    bit_id: u64,
) -> Result<(), Error> {
    // BEGIN RDKIT CPP FUNCTION FingerprintGenerator::getFingerprintHelper count increment
    // RDKit❗✔️: res->setVal(bitId, res->getVal(bitId) + 1);
    // END RDKIT CPP FUNCTION FingerprintGenerator::getFingerprintHelper count increment
    let current = result.value_at(bit_id)?;
    let incremented = current.checked_add(1).ok_or(
        FingerprintError::UndefinedArithmetic {
            site: "FingerprintGenerator::getFingerprintHelper count increment",
        }
        .into(),
    )?;
    result.set_at(bit_id, incremented)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use cosmolkit_core::{ValenceParams, assign_valence, fast_find_rings};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Element, TopologyBlock,
    };
    use cosmolkit_smiles::{SmilesParseParams, parse_smiles};

    use crate::additional_output::AdditionalOutput;
    use crate::morgan::{MorganAtomEnvironment, MorganParams, get_morgan_generator};
    use crate::sparse_counts::SparseCountFingerprint;

    use super::{
        FingerprintArguments, FingerprintFuncArguments, accumulate_morgan_sparse_counts,
        get_count_fingerprint, get_fingerprint, get_sparse_count_fingerprint,
        get_sparse_fingerprint, increment_morgan_sparse_count,
    };
    use crate::{FingerprintError, MorganError};

    #[test]
    fn fingerprint_morgan_g01_constructor_parameter_product() {
        const FP_SIZES: [u32; 2] = [0, 2048];
        const BITS_PER_FEATURE: [u32; 3] = [0, 1, 2];
        const COUNT_SIMULATION: [bool; 2] = [false, true];
        const COUNT_BOUNDS: [&[u32]; 6] = [
            &[],
            &[1, 2, 4, 8],
            &[4, 1, 2],
            &[1, 1],
            &[0, 1],
            &[2, 0, 1, 0],
        ];

        // Fixed source-derived outcomes by bits-per-feature, simulation, and
        // bounds case. fpSize does not change this constructor's validation.
        const EXPECTED_ERRORS: [[[Option<&str>; 6]; 2]; 3] = [
            [
                [Some("numBitsPerFeature must be >0"); 6],
                [
                    Some("bad count bounds provided"),
                    Some("numBitsPerFeature must be >0"),
                    Some("numBitsPerFeature must be >0"),
                    Some("numBitsPerFeature must be >0"),
                    Some("numBitsPerFeature must be >0"),
                    Some("numBitsPerFeature must be >0"),
                ],
            ],
            [
                [None; 6],
                [
                    Some("bad count bounds provided"),
                    None,
                    None,
                    None,
                    None,
                    None,
                ],
            ],
            [
                [None; 6],
                [
                    Some("bad count bounds provided"),
                    None,
                    None,
                    None,
                    None,
                    None,
                ],
            ],
        ];

        let mut calls = 0;
        for &fp_size in &FP_SIZES {
            for (bits_index, &bits_per_feature) in BITS_PER_FEATURE.iter().enumerate() {
                for (simulation_index, &count_simulation) in COUNT_SIMULATION.iter().enumerate() {
                    for (bounds_index, &count_bounds) in COUNT_BOUNDS.iter().enumerate() {
                        let actual = FingerprintArguments::new(
                            count_simulation,
                            count_bounds.to_vec(),
                            fp_size,
                            bits_per_feature,
                            false,
                        );
                        calls += 1;

                        match EXPECTED_ERRORS[bits_index][simulation_index][bounds_index] {
                            Some(what) => assert_eq!(
                                actual,
                                Err(FingerprintError::PreconditionViolation { what }),
                                "fpSize={fp_size}, bitsPerFeature={bits_per_feature}, \
                                 countSimulation={count_simulation}, countBounds={count_bounds:?}"
                            ),
                            None => {
                                let actual = actual.unwrap_or_else(|error| {
                                    panic!(
                                        "unexpected error for fpSize={fp_size}, \
                                         bitsPerFeature={bits_per_feature}, \
                                         countSimulation={count_simulation}, \
                                         countBounds={count_bounds:?}: {error}"
                                    )
                                });
                                assert_eq!(actual.count_simulation, count_simulation);
                                assert!(!actual.include_chirality);
                                assert_eq!(actual.count_bounds.as_slice(), count_bounds);
                                assert_eq!(actual.fp_size, fp_size);
                                assert_eq!(actual.bits_per_feature, bits_per_feature);
                            }
                        }
                    }
                }
            }
        }

        assert_eq!(calls, 72);
    }

    #[test]
    fn fingerprint_morgan_g01_defaults_and_value_preservation() {
        let defaults = FingerprintArguments::default();
        assert!(!defaults.count_simulation);
        assert!(!defaults.include_chirality);
        assert!(defaults.count_bounds.is_empty());
        assert_eq!(defaults.fp_size, 2048);
        assert_eq!(defaults.bits_per_feature, 1);

        let defaulted =
            FingerprintArguments::new_with_source_default_options(false, vec![3, 0, 3], 0)
                .expect("source defaults remain valid without count simulation");
        assert!(!defaulted.count_simulation);
        assert!(!defaulted.include_chirality);
        assert_eq!(defaulted.count_bounds, [3, 0, 3]);
        assert_eq!(defaulted.fp_size, 0);
        assert_eq!(defaulted.bits_per_feature, 1);

        let explicit = FingerprintArguments::new(false, vec![0, 0], 2048, 2, true)
            .expect("explicit chirality and repeated zero bounds are preserved");
        assert!(explicit.include_chirality);
        assert_eq!(explicit.count_bounds, [0, 0]);
        assert_eq!(explicit.fp_size, 2048);
        assert_eq!(explicit.bits_per_feature, 2);
        let constructor_calls = 2;
        assert_eq!(constructor_calls, 2);
    }

    #[test]
    fn fingerprint_morgan_g08_sparse_accumulation_and_all_output_channels() {
        struct ExpectedCase {
            fp_size: u64,
            bits_per_feature: u32,
            emitted_by_environment: [&'static [u64]; 4],
            sparse_counts: &'static [(u64, i32)],
        }

        // The environment seeds and additional-bit draws are the fixed R02
        // literals from the pinned Boost 1.85 source closure. Each expected
        // row below freezes raw IDs for fpSize=0 and source modulo results for
        // the two nonzero sizes; repeated seed 42 and fpSize=1 retain collisions.
        const CASES: [ExpectedCase; 9] = [
            ExpectedCase {
                fp_size: 0,
                bits_per_feature: 1,
                emitted_by_environment: [&[42], &[0x8000_0000], &[0xffff_ffff], &[42]],
                sparse_counts: &[(42, 2), (0x8000_0000, 1), (0xffff_ffff, 1)],
            },
            ExpectedCase {
                fp_size: 0,
                bits_per_feature: 2,
                emitted_by_environment: [
                    &[42, 0x3825_9ee2],
                    &[0x8000_0000, 0x68b3_35f0],
                    &[0xffff_ffff, 0x0d24_1918],
                    &[42, 0x3825_9ee2],
                ],
                sparse_counts: &[
                    (42, 2),
                    (0x0d24_1918, 1),
                    (0x3825_9ee2, 2),
                    (0x68b3_35f0, 1),
                    (0x8000_0000, 1),
                    (0xffff_ffff, 1),
                ],
            },
            ExpectedCase {
                fp_size: 0,
                bits_per_feature: 3,
                emitted_by_environment: [
                    &[42, 0x3825_9ee2, 0x290a_d8a5],
                    &[0x8000_0000, 0x68b3_35f0, 0x46e7_998d],
                    &[0xffff_ffff, 0x0d24_1918, 0x5979_48e5],
                    &[42, 0x3825_9ee2, 0x290a_d8a5],
                ],
                sparse_counts: &[
                    (42, 2),
                    (0x0d24_1918, 1),
                    (0x290a_d8a5, 2),
                    (0x3825_9ee2, 2),
                    (0x46e7_998d, 1),
                    (0x5979_48e5, 1),
                    (0x68b3_35f0, 1),
                    (0x8000_0000, 1),
                    (0xffff_ffff, 1),
                ],
            },
            ExpectedCase {
                fp_size: 1,
                bits_per_feature: 1,
                emitted_by_environment: [&[0], &[0], &[0], &[0]],
                sparse_counts: &[(0, 4)],
            },
            ExpectedCase {
                fp_size: 1,
                bits_per_feature: 2,
                emitted_by_environment: [&[0, 0], &[0, 0], &[0, 0], &[0, 0]],
                sparse_counts: &[(0, 8)],
            },
            ExpectedCase {
                fp_size: 1,
                bits_per_feature: 3,
                emitted_by_environment: [&[0, 0, 0], &[0, 0, 0], &[0, 0, 0], &[0, 0, 0]],
                sparse_counts: &[(0, 12)],
            },
            ExpectedCase {
                fp_size: 5,
                bits_per_feature: 1,
                emitted_by_environment: [&[2], &[3], &[0], &[2]],
                sparse_counts: &[(0, 1), (2, 2), (3, 1)],
            },
            ExpectedCase {
                fp_size: 5,
                bits_per_feature: 2,
                emitted_by_environment: [&[2, 2], &[3, 1], &[0, 3], &[2, 2]],
                sparse_counts: &[(0, 1), (1, 1), (2, 4), (3, 2)],
            },
            ExpectedCase {
                fp_size: 5,
                bits_per_feature: 3,
                emitted_by_environment: [&[2, 2, 2], &[3, 1, 0], &[0, 3, 1], &[2, 2, 2]],
                sparse_counts: &[(0, 2), (1, 2), (2, 6), (3, 2)],
            },
        ];
        const SEEDS: [u32; 4] = [42, 0x8000_0000, u32::MAX, 42];
        const CENTERS: [u32; 4] = [0, 1, 2, 0];
        const LAYERS: [u32; 4] = [0, 1, 2, 1];
        const ATOMS_PER_BIT_ROWS: [&[i32]; 4] = [&[0], &[1, 0, 2], &[2, 0, 1], &[0, 1]];
        const COUNT_SIMULATION: [bool; 2] = [false, true];

        let atoms = (0..3)
            .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect();
        let bonds = [(0, 1), (1, 2)]
            .into_iter()
            .enumerate()
            .map(|(index, (begin, end))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect();
        let topology = TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("the fixed three-atom chain is valid");

        let mut calls = 0;
        for case in &CASES {
            for count_simulation in COUNT_SIMULATION {
                for output_mask in 0_u8..32 {
                    let fingerprint_arguments = FingerprintArguments::new(
                        count_simulation,
                        vec![0, 1, 5],
                        case.fp_size as u32,
                        case.bits_per_feature,
                        false,
                    )
                    .expect("the source parameter row is valid");
                    let environments = (0..SEEDS.len())
                        .map(|index| {
                            MorganAtomEnvironment::<u64>::new(
                                SEEDS[index],
                                CENTERS[index],
                                LAYERS[index],
                                &topology,
                            )
                        })
                        .collect();
                    let mut actual_output = AdditionalOutput {
                        atom_counts: (output_mask & 0b00001 != 0).then(|| vec![0; 3]),
                        atom_to_bits: (output_mask & 0b00010 != 0)
                            .then(|| vec![Vec::new(), Vec::new(), Vec::new()]),
                        bit_info_map: (output_mask & 0b00100 != 0).then(BTreeMap::new),
                        bit_paths: (output_mask & 0b01000 != 0).then(BTreeMap::new),
                        atoms_per_bit: (output_mask & 0b10000 != 0)
                            .then(|| BTreeMap::from([(u64::MAX, vec![vec![77]])])),
                    };
                    let mut expected_output = AdditionalOutput {
                        atom_counts: (output_mask & 0b00001 != 0).then(|| vec![0; 3]),
                        atom_to_bits: (output_mask & 0b00010 != 0)
                            .then(|| vec![Vec::new(), Vec::new(), Vec::new()]),
                        bit_info_map: (output_mask & 0b00100 != 0).then(BTreeMap::new),
                        bit_paths: (output_mask & 0b01000 != 0).then(BTreeMap::new),
                        atoms_per_bit: (output_mask & 0b10000 != 0)
                            .then(|| BTreeMap::from([(u64::MAX, vec![vec![77]])])),
                    };

                    // This expected-output pass consumes only the fixed source
                    // rows above. It deliberately starts after the pinned outer
                    // reinit stage: four fields are already empty/zero, while
                    // atomsPerBit retains its preexisting source value.
                    for (environment_index, emitted_ids) in
                        case.emitted_by_environment.iter().enumerate()
                    {
                        for &bit_id in *emitted_ids {
                            if let Some(atom_counts) = &mut expected_output.atom_counts {
                                atom_counts[CENTERS[environment_index] as usize] += 1;
                            }
                            if let Some(atom_to_bits) = &mut expected_output.atom_to_bits {
                                atom_to_bits[CENTERS[environment_index] as usize].push(bit_id);
                            }
                            if let Some(bit_info_map) = &mut expected_output.bit_info_map {
                                bit_info_map
                                    .entry(bit_id)
                                    .or_default()
                                    .push((CENTERS[environment_index], LAYERS[environment_index]));
                            }
                            if let Some(atoms_per_bit) = &mut expected_output.atoms_per_bit {
                                atoms_per_bit
                                    .entry(bit_id)
                                    .or_default()
                                    .push(ATOMS_PER_BIT_ROWS[environment_index].to_vec());
                            }
                        }
                    }

                    let result = accumulate_morgan_sparse_counts(
                        environments,
                        &fingerprint_arguments,
                        &[11, 12, 13],
                        &[21, 22],
                        case.fp_size,
                        Some(&mut actual_output),
                    )
                    .expect("fixed source environment IDs are valid sparse indices");
                    let expected_counts: BTreeMap<_, _> =
                        case.sparse_counts.iter().copied().collect();

                    assert_eq!(
                        result.length(),
                        if case.fp_size == 0 {
                            u64::MAX
                        } else {
                            case.fp_size
                        },
                        "result length, fpSize={}, bitsPerFeature={}, countSimulation={count_simulation}",
                        case.fp_size,
                        case.bits_per_feature
                    );
                    assert_eq!(
                        result.nonzero_elements(),
                        &expected_counts,
                        "fixed counts, fpSize={}, bitsPerFeature={}, countSimulation={count_simulation}",
                        case.fp_size,
                        case.bits_per_feature
                    );
                    assert_eq!(
                        actual_output, expected_output,
                        "all optional channels, fpSize={}, bitsPerFeature={}, countSimulation={count_simulation}, mask={output_mask:#07b}",
                        case.fp_size, case.bits_per_feature
                    );

                    if case.fp_size == 5 && case.bits_per_feature == 3 && output_mask == 0b1_1111 {
                        assert_eq!(
                            result.nonzero_elements(),
                            &BTreeMap::from([(0, 2), (1, 2), (2, 6), (3, 2),])
                        );
                        assert_eq!(
                            actual_output,
                            AdditionalOutput {
                                atom_counts: Some(vec![6, 3, 3]),
                                atom_to_bits: Some(vec![
                                    vec![2, 2, 2, 2, 2, 2],
                                    vec![3, 1, 0],
                                    vec![0, 3, 1],
                                ]),
                                bit_info_map: Some(BTreeMap::from([
                                    (0, vec![(1, 1), (2, 2)]),
                                    (1, vec![(1, 1), (2, 2)]),
                                    (2, vec![(0, 0), (0, 0), (0, 0), (0, 1), (0, 1), (0, 1)]),
                                    (3, vec![(1, 1), (2, 2)]),
                                ])),
                                bit_paths: Some(BTreeMap::new()),
                                atoms_per_bit: Some(BTreeMap::from([
                                    (0, vec![vec![1, 0, 2], vec![2, 0, 1]]),
                                    (1, vec![vec![1, 0, 2], vec![2, 0, 1]]),
                                    (
                                        2,
                                        vec![
                                            vec![0],
                                            vec![0],
                                            vec![0],
                                            vec![0, 1],
                                            vec![0, 1],
                                            vec![0, 1],
                                        ],
                                    ),
                                    (3, vec![vec![1, 0, 2], vec![2, 0, 1]]),
                                    (u64::MAX, vec![vec![77]]),
                                ])),
                            },
                            "full fixed AdditionalOutput rows preserve collisions and source append order"
                        );
                    }
                    calls += 1;
                }
            }
        }

        assert_eq!(calls, 576);
    }

    #[test]
    fn fingerprint_morgan_g02_borrowed_arguments_parameter_product() {
        const FROM_ATOMS: [(&str, Option<&[u32]>); 5] = [
            ("absent", None),
            ("empty", Some(&[])),
            ("single", Some(&[1])),
            ("duplicate", Some(&[1, 1])),
            ("out_of_range_for_three_atoms", Some(&[3])),
        ];
        const IGNORE_ATOMS: [(&str, Option<&[u32]>); 5] = [
            ("absent", None),
            ("empty", Some(&[])),
            ("single", Some(&[2])),
            ("duplicate", Some(&[2, 2])),
            ("out_of_range_for_three_atoms", Some(&[3])),
        ];
        const ATOM_INVARIANTS: [(&str, Option<&[u32]>); 5] = [
            ("absent", None),
            ("empty", Some(&[])),
            ("exactly_three", Some(&[11, 12, 13])),
            ("short_two", Some(&[11, 12])),
            ("long_four", Some(&[11, 12, 13, 14])),
        ];
        const BOND_INVARIANTS: [(&str, Option<&[u32]>); 5] = [
            ("absent", None),
            ("empty", Some(&[])),
            ("exactly_two", Some(&[21, 22])),
            ("short_one", Some(&[21])),
            ("long_three", Some(&[21, 22, 23])),
        ];
        const CONFORMER_IDS: [(&str, i32); 3] = [
            ("source_default", -1),
            ("zero", 0),
            ("nonexistent_id_sentinel", i32::MAX),
        ];

        let defaults = FingerprintFuncArguments::default();
        assert_eq!(defaults.from_atoms, None);
        assert_eq!(defaults.ignore_atoms, None);
        assert_eq!(defaults.custom_atom_invariants, None);
        assert_eq!(defaults.custom_bond_invariants, None);
        assert_eq!(defaults.conformer_id, -1);

        // This is a holder-level test: invalid fromAtoms is preserved but is
        // not passed to Boost's assertion/unchecked-index consumer here.
        // Morgan ignores ignoreAtoms and conformer_id at its source consumer;
        // these cases verify the detached call keeps the supplied values.
        let mut calls = 0;
        for &(from_name, from_atoms) in &FROM_ATOMS {
            for &(ignore_name, ignore_atoms) in &IGNORE_ATOMS {
                for &(atom_name, atom_invariants) in &ATOM_INVARIANTS {
                    for &(bond_name, bond_invariants) in &BOND_INVARIANTS {
                        for &(conf_name, conformer_id) in &CONFORMER_IDS {
                            let actual = FingerprintFuncArguments::new(
                                from_atoms,
                                ignore_atoms,
                                atom_invariants,
                                bond_invariants,
                                conformer_id,
                            );
                            assert_eq!(actual.from_atoms, from_atoms, "fromAtoms={from_name}");
                            assert_eq!(
                                actual.ignore_atoms, ignore_atoms,
                                "ignoreAtoms={ignore_name}"
                            );
                            assert_eq!(
                                actual.custom_atom_invariants, atom_invariants,
                                "customAtomInvariants={atom_name}"
                            );
                            assert_eq!(
                                actual.custom_bond_invariants, bond_invariants,
                                "customBondInvariants={bond_name}"
                            );
                            assert_eq!(actual.conformer_id, conformer_id, "confId={conf_name}");
                            calls += 1;
                        }
                    }
                }
            }
        }

        assert_eq!(calls, 1_875);
    }

    #[test]
    fn fingerprint_morgan_o01_sparse_count_full_source_product() {
        struct ExpectedCase {
            bits_per_feature: u32,
            emitted_by_environment: [&'static [u64]; 3],
            sparse_counts: &'static [(u64, i32)],
        }

        // Fixed Boost 1.85 rows reused from G08's pinned source vectors.
        // The source wrapper does not fold these raw u32 environment IDs,
        // regardless of configured fpSize or count-simulation parameters.
        const CASES: [ExpectedCase; 3] = [
            ExpectedCase {
                bits_per_feature: 1,
                emitted_by_environment: [&[0xffff_ffff], &[0x8000_0000], &[0xffff_ffff]],
                sparse_counts: &[(0x8000_0000, 1), (0xffff_ffff, 2)],
            },
            ExpectedCase {
                bits_per_feature: 2,
                emitted_by_environment: [
                    &[0xffff_ffff, 0x0d24_1918],
                    &[0x8000_0000, 0x68b3_35f0],
                    &[0xffff_ffff, 0x0d24_1918],
                ],
                sparse_counts: &[
                    (0x0d24_1918, 2),
                    (0x68b3_35f0, 1),
                    (0x8000_0000, 1),
                    (0xffff_ffff, 2),
                ],
            },
            ExpectedCase {
                bits_per_feature: 3,
                emitted_by_environment: [
                    &[0xffff_ffff, 0x0d24_1918, 0x5979_48e5],
                    &[0x8000_0000, 0x68b3_35f0, 0x46e7_998d],
                    &[0xffff_ffff, 0x0d24_1918, 0x5979_48e5],
                ],
                sparse_counts: &[
                    (0x0d24_1918, 2),
                    (0x46e7_998d, 1),
                    (0x5979_48e5, 2),
                    (0x68b3_35f0, 1),
                    (0x8000_0000, 1),
                    (0xffff_ffff, 2),
                ],
            },
        ];
        const ATOM_INVARIANTS: [u32; 3] = [u32::MAX, 0x8000_0000, u32::MAX];
        const FP_SIZES: [u32; 3] = [0, 1, 2048];
        const COUNT_SIMULATION: [bool; 2] = [false, true];
        const EMPTY_FROM_ATOMS: &[u32] = &[];
        const FROM_ATOMS: [Option<&[u32]>; 2] = [None, Some(EMPTY_FROM_ATOMS)];

        let record = parse_smiles("CCO", &SmilesParseParams::default())
            .expect("the fixed CCO owner input parses with source defaults");
        assert_eq!(record.topology.atoms.len(), 3);
        let original_topology = record.topology.clone();
        let original_properties = record.properties.clone();
        let valence = assign_valence(&record.topology, &ValenceParams::default())
            .expect("fixed CCO has source-valid valence");
        let rings = fast_find_rings(&record.topology).expect("fixed CCO has ring information");
        let mut calls = 0;

        for case in &CASES {
            for fp_size in FP_SIZES {
                for count_simulation in COUNT_SIMULATION {
                    let params = MorganParams {
                        radius: 0,
                        fp_size,
                        count_simulation,
                        bits_per_feature: case.bits_per_feature,
                        ..MorganParams::default()
                    };
                    let generator = get_morgan_generator(&params)
                        .expect("every fixed O01 configuration satisfies source preconditions");

                    for &from_atoms in &FROM_ATOMS {
                        for output_mask in 0_u8..32 {
                            let arguments = FingerprintFuncArguments::new(
                                from_atoms,
                                None,
                                Some(&ATOM_INVARIANTS),
                                None,
                                -1,
                            );
                            let mut actual_output = AdditionalOutput {
                                atom_counts: (output_mask & 0b00001 != 0).then(|| vec![91]),
                                atom_to_bits: (output_mask & 0b00010 != 0)
                                    .then(|| vec![vec![81], vec![82]]),
                                bit_info_map: (output_mask & 0b00100 != 0)
                                    .then(|| BTreeMap::from([(17, vec![(9, 9)])])),
                                bit_paths: (output_mask & 0b01000 != 0)
                                    .then(|| BTreeMap::from([(17, vec![vec![9]])])),
                                atoms_per_bit: (output_mask & 0b10000 != 0)
                                    .then(|| BTreeMap::from([(u64::MAX, vec![vec![77]])])),
                            };
                            let mut expected_output = AdditionalOutput {
                                atom_counts: (output_mask & 0b00001 != 0).then(|| vec![0; 3]),
                                atom_to_bits: (output_mask & 0b00010 != 0)
                                    .then(|| vec![Vec::new(), Vec::new(), Vec::new()]),
                                bit_info_map: (output_mask & 0b00100 != 0).then(BTreeMap::new),
                                bit_paths: (output_mask & 0b01000 != 0).then(BTreeMap::new),
                                atoms_per_bit: (output_mask & 0b10000 != 0)
                                    .then(|| BTreeMap::from([(u64::MAX, vec![vec![77]])])),
                            };

                            if from_atoms.is_none() {
                                for (center, emitted_ids) in
                                    case.emitted_by_environment.iter().enumerate()
                                {
                                    for &bit_id in *emitted_ids {
                                        if let Some(atom_counts) = &mut expected_output.atom_counts
                                        {
                                            atom_counts[center] += 1;
                                        }
                                        if let Some(atom_to_bits) =
                                            &mut expected_output.atom_to_bits
                                        {
                                            atom_to_bits[center].push(bit_id);
                                        }
                                        if let Some(bit_info_map) =
                                            &mut expected_output.bit_info_map
                                        {
                                            bit_info_map
                                                .entry(bit_id)
                                                .or_default()
                                                .push((center as u32, 0));
                                        }
                                        if let Some(atoms_per_bit) =
                                            &mut expected_output.atoms_per_bit
                                        {
                                            atoms_per_bit
                                                .entry(bit_id)
                                                .or_default()
                                                .push(vec![center as i32]);
                                        }
                                    }
                                }
                            }

                            let result = get_sparse_count_fingerprint(
                                &record.topology,
                                &record.properties,
                                &valence,
                                &rings,
                                &generator,
                                &arguments,
                                Some(&mut actual_output),
                            )
                            .expect("fixed source raw IDs fit the sparse u64 index domain");
                            let expected_counts: BTreeMap<u64, i32> = if from_atoms.is_none() {
                                case.sparse_counts.iter().copied().collect()
                            } else {
                                BTreeMap::new()
                            };

                            assert_eq!(
                                result.length(),
                                u64::MAX,
                                "the source wrapper keeps the u64 output domain despite configured fpSize={fp_size}, countSimulation={count_simulation}, bitsPerFeature={}",
                                case.bits_per_feature
                            );
                            assert_eq!(
                                result.nonzero_elements(),
                                &expected_counts,
                                "fixed raw sparse counts for fpSize={fp_size}, countSimulation={count_simulation}, bitsPerFeature={}, fromAtoms={from_atoms:?}",
                                case.bits_per_feature
                            );
                            assert_eq!(
                                actual_output, expected_output,
                                "complete optional output for fpSize={fp_size}, countSimulation={count_simulation}, bitsPerFeature={}, fromAtoms={from_atoms:?}, mask={output_mask:#07b}",
                                case.bits_per_feature
                            );

                            if fp_size == 1
                                && count_simulation
                                && case.bits_per_feature == 3
                                && from_atoms.is_none()
                                && output_mask == 0b1_1111
                            {
                                assert_eq!(
                                    result.nonzero_elements(),
                                    &BTreeMap::from([
                                        (0x0d24_1918, 2),
                                        (0x46e7_998d, 1),
                                        (0x5979_48e5, 2),
                                        (0x68b3_35f0, 1),
                                        (0x8000_0000, 1),
                                        (0xffff_ffff, 2),
                                    ])
                                );
                                assert_eq!(
                                    actual_output,
                                    AdditionalOutput {
                                        atom_counts: Some(vec![3, 3, 3]),
                                        atom_to_bits: Some(vec![
                                            vec![0xffff_ffff, 0x0d24_1918, 0x5979_48e5],
                                            vec![0x8000_0000, 0x68b3_35f0, 0x46e7_998d],
                                            vec![0xffff_ffff, 0x0d24_1918, 0x5979_48e5],
                                        ]),
                                        bit_info_map: Some(BTreeMap::from([
                                            (0x0d24_1918, vec![(0, 0), (2, 0)]),
                                            (0x46e7_998d, vec![(1, 0)]),
                                            (0x5979_48e5, vec![(0, 0), (2, 0)]),
                                            (0x68b3_35f0, vec![(1, 0)]),
                                            (0x8000_0000, vec![(1, 0)]),
                                            (0xffff_ffff, vec![(0, 0), (2, 0)]),
                                        ])),
                                        bit_paths: Some(BTreeMap::new()),
                                        atoms_per_bit: Some(BTreeMap::from([
                                            (0x0d24_1918, vec![vec![0], vec![2]]),
                                            (0x46e7_998d, vec![vec![1]]),
                                            (0x5979_48e5, vec![vec![0], vec![2]]),
                                            (0x68b3_35f0, vec![vec![1]]),
                                            (0x8000_0000, vec![vec![1]]),
                                            (0xffff_ffff, vec![vec![0], vec![2]]),
                                            (u64::MAX, vec![vec![77]]),
                                        ])),
                                    }
                                );
                            }
                            calls += 1;
                        }
                    }
                }
            }
        }

        assert_eq!(calls, 1_152);
        assert_eq!(record.topology, original_topology);
        assert_eq!(record.properties, original_properties);

        let mut terminal_index = SparseCountFingerprint::new(u64::MAX);
        terminal_index
            .set_value(u64::MAX, 1)
            .expect("source checkIndex accepts idx==length at the u64 maximum");
        assert_eq!(terminal_index.value(u64::MAX), Ok(1));
        let mut shorter_domain = SparseCountFingerprint::new(u64::MAX - 1);
        assert_eq!(
            shorter_domain.set_value(u64::MAX, 1),
            Err(FingerprintError::SparseIndexOutOfRange {
                index: u64::MAX,
                size: u64::MAX - 1,
            })
        );
    }

    #[test]
    fn fingerprint_morgan_o01_preconditions_reset_output_and_overflow_are_typed() {
        const SHORT_ATOM_INVARIANTS: [u32; 2] = [11, 12];
        const FULL_ATOM_INVARIANTS: [u32; 3] = [11, 12, 13];
        const SHORT_BOND_INVARIANTS: [u32; 1] = [21];

        let record = parse_smiles("CCO", &SmilesParseParams::default())
            .expect("the fixed CCO owner input parses with source defaults");
        let valence = assign_valence(&record.topology, &ValenceParams::default())
            .expect("fixed CCO has source-valid valence");
        let rings = fast_find_rings(&record.topology).expect("fixed CCO has ring information");
        let generator = get_morgan_generator(&MorganParams {
            radius: 0,
            ..MorganParams::default()
        })
        .expect("the fixed O01 generator satisfies source preconditions");
        let expected_reset_output = || AdditionalOutput {
            atom_counts: Some(vec![0; 3]),
            atom_to_bits: Some(vec![Vec::new(), Vec::new(), Vec::new()]),
            bit_info_map: Some(BTreeMap::new()),
            bit_paths: Some(BTreeMap::new()),
            atoms_per_bit: Some(BTreeMap::from([(u64::MAX, vec![vec![77]])])),
        };
        let mut error_cases = 0;

        for (atom_invariants, bond_invariants, expected_error) in [
            (
                &SHORT_ATOM_INVARIANTS[..],
                &SHORT_BOND_INVARIANTS[..],
                "bad atom invariants size",
            ),
            (
                &FULL_ATOM_INVARIANTS[..],
                &SHORT_BOND_INVARIANTS[..],
                "bad bond invariants size",
            ),
        ] {
            let arguments = FingerprintFuncArguments::new(
                None,
                None,
                Some(atom_invariants),
                Some(bond_invariants),
                -1,
            );
            let mut actual_output = AdditionalOutput {
                atom_counts: Some(vec![71]),
                atom_to_bits: Some(vec![vec![72]]),
                bit_info_map: Some(BTreeMap::from([(73, vec![(74, 75)])])),
                bit_paths: Some(BTreeMap::from([(76, vec![vec![77]])])),
                atoms_per_bit: Some(BTreeMap::from([(u64::MAX, vec![vec![77]])])),
            };

            let error = get_sparse_count_fingerprint(
                &record.topology,
                &record.properties,
                &valence,
                &rings,
                &generator,
                &arguments,
                Some(&mut actual_output),
            )
            .expect_err("short source-used invariant vectors fail at environment entry");
            assert!(
                matches!(
                    error,
                    MorganError::Fingerprint(FingerprintError::PreconditionViolation { what })
                        if what == expected_error
                ),
                "source atom-before-bond precondition: {error}"
            );
            assert_eq!(
                actual_output,
                expected_reset_output(),
                "source output reinitialization occurs before invariant length rejection"
            );
            error_cases += 1;
        }

        let mut overflowing = SparseCountFingerprint::new(u64::MAX);
        overflowing
            .set_value(u64::MAX, i32::MAX)
            .expect("the full u64 source index domain accepts its terminal index");
        let overflow = increment_morgan_sparse_count(&mut overflowing, u64::MAX);
        assert!(matches!(
            overflow,
            Err(MorganError::Fingerprint(
                FingerprintError::UndefinedArithmetic {
                    site: "FingerprintGenerator::getFingerprintHelper count increment"
                }
            ))
        ));
        assert_eq!(overflowing.value(u64::MAX), Ok(i32::MAX));
        error_cases += 1;

        assert_eq!(error_cases, 3);
    }

    #[test]
    fn fingerprint_morgan_o02_sparse_bits_count_simulation_and_output_product() {
        const BOUND_CASES: [(&str, &[u32]); 4] = [
            ("default", &[1, 2, 4, 8]),
            ("unsorted", &[4, 1, 8, 2]),
            ("duplicates", &[1, 1, 2, 2]),
            ("empty", &[]),
        ];
        const ATOM_INVARIANTS: [u32; 3] = [1, 1_073_741_824, 0x8000_0000];
        const NO_SIMULATION_BITS: &[i32] = &[i32::MIN, 1, 1_073_741_824];
        const NO_SIMULATION_ATOM_BITS: [&[u32]; 3] = [&[1], &[1_073_741_824], &[0x8000_0000]];
        const SIMULATION_BITS: [&[i32]; 3] = [&[4, 5, 8], &[5, 7, 9], &[4, 5, 6, 7, 8, 9]];
        const SIMULATION_ATOM_BITS: [[&[u32]; 3]; 3] = [
            [&[4, 5], &[4, 5], &[8]],
            [&[5, 7], &[5, 7], &[9]],
            [&[4, 5, 6, 7], &[4, 5, 6, 7], &[8, 9]],
        ];

        let record = parse_smiles("CCO", &SmilesParseParams::default())
            .expect("the fixed CCO owner input parses with source defaults");
        assert_eq!(record.topology.atoms.len(), 3);
        let original_topology = record.topology.clone();
        let original_properties = record.properties.clone();
        let valence = assign_valence(&record.topology, &ValenceParams::default())
            .expect("fixed CCO has source-valid valence");
        let rings = fast_find_rings(&record.topology).expect("fixed CCO has ring information");
        let arguments = FingerprintFuncArguments::new(None, None, Some(&ATOM_INVARIANTS), None, -1);

        let make_expected_output = |mask: u8, per_atom_bits: [&[u32]; 3]| {
            let mut expected_bit_info = BTreeMap::<u64, Vec<(u32, u32)>>::new();
            let mut expected_atoms_per_bit = BTreeMap::from([(u64::MAX, vec![vec![77]])]);
            for (atom_index, bit_ids) in per_atom_bits.iter().enumerate() {
                for &bit_id in *bit_ids {
                    let key = u64::from(bit_id);
                    expected_bit_info
                        .entry(key)
                        .or_default()
                        .push((atom_index as u32, 0));
                    expected_atoms_per_bit
                        .entry(key)
                        .or_default()
                        .push(vec![atom_index as i32]);
                }
            }

            AdditionalOutput {
                atom_counts: (mask & 0b00001 != 0).then(|| vec![1, 1, 1]),
                atom_to_bits: (mask & 0b00010 != 0).then(|| {
                    per_atom_bits
                        .iter()
                        .map(|bit_ids| bit_ids.iter().copied().map(u64::from).collect())
                        .collect()
                }),
                bit_info_map: (mask & 0b00100 != 0).then_some(expected_bit_info),
                bit_paths: (mask & 0b01000 != 0).then(BTreeMap::new),
                atoms_per_bit: (mask & 0b10000 != 0).then_some(expected_atoms_per_bit),
            }
        };

        let mut parameter_combinations = 0;
        let mut rejected_empty_simulation_configs = 0;
        let mut output_mask_calls = 0;
        let mut absent_output_calls = 0;
        for count_simulation in [false, true] {
            for (bounds_index, &(bounds_name, bounds)) in BOUND_CASES.iter().enumerate() {
                parameter_combinations += 1;
                let generator_result = get_morgan_generator(&MorganParams {
                    radius: 0,
                    fp_size: 1,
                    count_simulation,
                    count_bounds: bounds.to_vec(),
                    ..MorganParams::default()
                });

                if count_simulation && bounds.is_empty() {
                    let error = generator_result
                        .expect_err("the source rejects empty count bounds when simulation is on");
                    assert!(matches!(
                        error,
                        MorganError::Fingerprint(
                            FingerprintError::PreconditionViolation { what }
                        ) if what == "bad count bounds provided"
                    ));
                    rejected_empty_simulation_configs += 1;
                    continue;
                }

                let generator = generator_result
                    .expect("all other fixed O02 configurations satisfy source preconditions");
                let (expected_bits, expected_per_atom_bits) = if count_simulation {
                    (
                        SIMULATION_BITS[bounds_index],
                        SIMULATION_ATOM_BITS[bounds_index],
                    )
                } else {
                    (NO_SIMULATION_BITS, NO_SIMULATION_ATOM_BITS)
                };

                for output_mask in 0_u8..32 {
                    let mut actual_output = AdditionalOutput {
                        atom_counts: (output_mask & 0b00001 != 0).then(|| vec![91]),
                        atom_to_bits: (output_mask & 0b00010 != 0).then(|| vec![vec![81]]),
                        bit_info_map: (output_mask & 0b00100 != 0)
                            .then(|| BTreeMap::from([(17, vec![(9, 9)])])),
                        bit_paths: (output_mask & 0b01000 != 0)
                            .then(|| BTreeMap::from([(17, vec![vec![9]])])),
                        atoms_per_bit: (output_mask & 0b10000 != 0)
                            .then(|| BTreeMap::from([(u64::MAX, vec![vec![77]])])),
                    };
                    let result = get_sparse_fingerprint(
                        &record.topology,
                        &record.properties,
                        &valence,
                        &rings,
                        &generator,
                        &arguments,
                        Some(&mut actual_output),
                    )
                    .expect("fixed source bit indices fit the u32 sparse vector domain");

                    assert_eq!(result.n_bits(), u32::MAX);
                    assert_eq!(
                        result.on_bits().as_slice(),
                        expected_bits,
                        "fixed bits for countSimulation={count_simulation}, bounds={bounds_name}, mask={output_mask:#07b}"
                    );
                    assert_eq!(
                        actual_output,
                        make_expected_output(output_mask, expected_per_atom_bits),
                        "complete output for countSimulation={count_simulation}, bounds={bounds_name}, mask={output_mask:#07b}"
                    );
                    output_mask_calls += 1;
                }

                let result_without_output = get_sparse_fingerprint(
                    &record.topology,
                    &record.properties,
                    &valence,
                    &rings,
                    &generator,
                    &arguments,
                    None,
                )
                .expect("a missing output pointer does not change sparse bits");
                assert_eq!(result_without_output.n_bits(), u32::MAX);
                assert_eq!(result_without_output.on_bits().as_slice(), expected_bits);
                absent_output_calls += 1;
            }
        }

        assert_eq!(parameter_combinations, 8);
        assert_eq!(rejected_empty_simulation_configs, 1);
        assert_eq!(output_mask_calls, 224);
        assert_eq!(absent_output_calls, 7);
        assert_eq!(output_mask_calls + absent_output_calls, 231);
        assert_eq!(record.topology, original_topology);
        assert_eq!(record.properties, original_properties);
    }

    #[test]
    fn fingerprint_morgan_o03_count_fingerprint_size_and_zero_size_product() {
        const FP_SIZES: [u32; 4] = [1, 2, 128, 2048];
        const ORDINARY_COUNTS: [&[(u32, i32)]; 4] = [
            &[(0, 3)],
            &[(0, 1), (1, 2)],
            &[(33, 1), (39, 1), (80, 1)],
            &[(80, 1), (807, 1), (1057, 1)],
        ];
        const COLLISION_COUNTS: [&[(u32, i32)]; 4] = [
            &[(0, 3)],
            &[(0, 1), (1, 2)],
            &[(1, 1), (2, 1), (3, 1)],
            &[(1, 1), (2, 1), (3, 1)],
        ];
        const COLLISION_ATOM_INVARIANTS: [u32; 3] = [1, 2, 3];

        let empty_record = parse_smiles("", &SmilesParseParams::default())
            .expect("the fixed empty SMILES input parses to the empty molecule");
        let empty_valence = assign_valence(&empty_record.topology, &ValenceParams::default())
            .expect("the empty molecule has an empty valid valence assignment");
        let empty_rings = fast_find_rings(&empty_record.topology)
            .expect("the empty molecule has an empty ring assignment");
        let ordinary_record = parse_smiles("CCO", &SmilesParseParams::default())
            .expect("the fixed CCO reference input parses with source defaults");
        let ordinary_valence = assign_valence(&ordinary_record.topology, &ValenceParams::default())
            .expect("fixed CCO has source-valid valence");
        let ordinary_rings =
            fast_find_rings(&ordinary_record.topology).expect("fixed CCO has ring information");
        let ordinary_arguments = FingerprintFuncArguments::default();
        let collision_arguments =
            FingerprintFuncArguments::new(None, None, Some(&COLLISION_ATOM_INVARIANTS), None, -1);
        let empty_arguments = FingerprintFuncArguments::default();

        let mut hashed_calls = 0;
        for (size_index, &fp_size) in FP_SIZES.iter().enumerate() {
            for count_simulation in [false, true] {
                let generator = get_morgan_generator(&MorganParams {
                    radius: 0,
                    fp_size,
                    count_simulation,
                    count_bounds: vec![1, 2, 4, 8],
                    ..MorganParams::default()
                })
                .expect("all fixed O03 generator configurations satisfy source preconditions");

                let empty = get_count_fingerprint(
                    &empty_record.topology,
                    &empty_record.properties,
                    &empty_valence,
                    &empty_rings,
                    &generator,
                    &empty_arguments,
                    None,
                )
                .expect("an empty environment map projects to every configured length");
                assert_eq!(empty.length(), fp_size);
                assert!(empty.nonzero_elements().is_empty());
                hashed_calls += 1;

                let ordinary = get_count_fingerprint(
                    &ordinary_record.topology,
                    &ordinary_record.properties,
                    &ordinary_valence,
                    &ordinary_rings,
                    &generator,
                    &ordinary_arguments,
                    None,
                )
                .expect("the fixed source CCO IDs fit each configured count vector");
                assert_eq!(ordinary.length(), fp_size);
                let ordinary_counts: Vec<(u32, i32)> = ordinary
                    .nonzero_elements()
                    .iter()
                    .map(|(&index, &count)| (index, count))
                    .collect();
                assert_eq!(ordinary_counts.as_slice(), ORDINARY_COUNTS[size_index]);
                hashed_calls += 1;

                let collision = get_count_fingerprint(
                    &ordinary_record.topology,
                    &ordinary_record.properties,
                    &ordinary_valence,
                    &ordinary_rings,
                    &generator,
                    &collision_arguments,
                    None,
                )
                .expect("the fixed custom IDs fit each configured count vector");
                assert_eq!(collision.length(), fp_size);
                let collision_counts: Vec<(u32, i32)> = collision
                    .nonzero_elements()
                    .iter()
                    .map(|(&index, &count)| (index, count))
                    .collect();
                assert_eq!(collision_counts.as_slice(), COLLISION_COUNTS[size_index]);
                hashed_calls += 1;
            }
        }

        let zero_generator = |count_simulation| {
            get_morgan_generator(&MorganParams {
                radius: 0,
                fp_size: 0,
                count_simulation,
                count_bounds: vec![1, 2, 4, 8],
                ..MorganParams::default()
            })
            .expect("zero-size O03 configuration passes constructor preconditions")
        };
        let mut zero_size_empty_calls = 0;
        let mut zero_size_typed_errors = 0;
        for count_simulation in [false, true] {
            let generator = zero_generator(count_simulation);
            let empty = get_count_fingerprint(
                &empty_record.topology,
                &empty_record.properties,
                &empty_valence,
                &empty_rings,
                &generator,
                &empty_arguments,
                None,
            )
            .expect("the source zero-length output succeeds when there are no environments");
            assert_eq!(empty.length(), 0);
            assert!(empty.nonzero_elements().is_empty());
            zero_size_empty_calls += 1;

            let ordinary_error = get_count_fingerprint(
                &ordinary_record.topology,
                &ordinary_record.properties,
                &ordinary_valence,
                &ordinary_rings,
                &generator,
                &ordinary_arguments,
                None,
            )
            .expect_err("a nonempty zero-length source count projection must fail");
            assert!(
                matches!(
                    ordinary_error,
                    MorganError::Fingerprint(FingerprintError::SparseIndexOutOfRange {
                        index: 864_662_311,
                        size: 0,
                    })
                ),
                "the first ordered fixed CCO source index fails projection at size zero"
            );
            zero_size_typed_errors += 1;

            let collision_error = get_count_fingerprint(
                &ordinary_record.topology,
                &ordinary_record.properties,
                &ordinary_valence,
                &ordinary_rings,
                &generator,
                &collision_arguments,
                None,
            )
            .expect_err("a nonempty custom zero-length projection must fail");
            assert!(
                matches!(
                    collision_error,
                    MorganError::Fingerprint(FingerprintError::SparseIndexOutOfRange {
                        index: 1,
                        size: 0,
                    })
                ),
                "the first ordered custom source index fails projection at size zero"
            );
            zero_size_typed_errors += 1;
        }

        assert_eq!(hashed_calls, 24);
        assert_eq!(zero_size_empty_calls, 2);
        assert_eq!(zero_size_typed_errors, 4);
    }

    #[test]
    fn fingerprint_morgan_o04_dense_bits_count_simulation_and_output_product() {
        const FP_SIZES: [u32; 4] = [1, 2, 128, 2048];
        const SIMULATION_BOUNDS: [&[u32]; 4] = [&[1], &[1], &[1, 2, 4, 8], &[1, 2, 4, 8]];
        const CUSTOM_ATOM_INVARIANTS: [u32; 3] = [1, 2, 1];
        const NO_SIMULATION_BITS: [&[u32]; 4] = [&[0], &[0, 1], &[1, 2], &[1, 2]];
        const NO_SIMULATION_ROWS: [[&[u32]; 3]; 4] = [
            [&[0], &[0], &[0]],
            [&[1], &[0], &[1]],
            [&[1], &[2], &[1]],
            [&[1], &[2], &[1]],
        ];
        const SIMULATION_BITS: [&[u32]; 4] = [&[], &[0, 1], &[4, 5, 8], &[4, 5, 8]];
        const SIMULATION_ROWS: [[&[u32]; 3]; 4] = [
            [&[], &[], &[]],
            [&[1], &[0], &[1]],
            [&[4, 5], &[8], &[4, 5]],
            [&[4, 5], &[8], &[4, 5]],
        ];
        const ZERO_SIZE_RAW_ROWS: [&[u32]; 3] = [&[1], &[2], &[1]];
        const CCO_NO_SIMULATION_BITS: [&[u32]; 4] =
            [&[0], &[0, 1], &[33, 39, 80], &[80, 807, 1057]];
        const CCO_SIMULATION_BITS: [&[u32]; 4] = [&[], &[0, 1], &[4, 28, 64], &[132, 320, 1180]];

        let record = parse_smiles("CCO", &SmilesParseParams::default())
            .expect("the fixed CCO dense-fingerprint input parses with source defaults");
        assert_eq!(record.topology.atoms.len(), 3);
        let original_topology = record.topology.clone();
        let original_properties = record.properties.clone();
        let valence = assign_valence(&record.topology, &ValenceParams::default())
            .expect("fixed CCO has source-valid valence");
        let rings = fast_find_rings(&record.topology).expect("fixed CCO has ring information");
        let collision_arguments =
            FingerprintFuncArguments::new(None, None, Some(&CUSTOM_ATOM_INVARIANTS), None, -1);
        let default_arguments = FingerprintFuncArguments::default();

        let make_initial_output = |mask: u8| AdditionalOutput {
            atom_counts: (mask & 0b00001 != 0).then(|| vec![91]),
            atom_to_bits: (mask & 0b00010 != 0).then(|| vec![vec![81]]),
            bit_info_map: (mask & 0b00100 != 0).then(|| BTreeMap::from([(17, vec![(9, 9)])])),
            bit_paths: (mask & 0b01000 != 0).then(|| BTreeMap::from([(17, vec![vec![9]])])),
            atoms_per_bit: (mask & 0b10000 != 0)
                .then(|| BTreeMap::from([(u64::MAX, vec![vec![77]])])),
        };
        let make_expected_output = |mask: u8, per_atom_bits: &[&[u32]]| {
            let mut bit_info_map = BTreeMap::<u64, Vec<(u32, u32)>>::new();
            let mut atoms_per_bit = BTreeMap::from([(u64::MAX, vec![vec![77]])]);
            let mut atom_to_bits = Vec::with_capacity(per_atom_bits.len());
            for (atom_index, bit_ids) in per_atom_bits.iter().enumerate() {
                let output_ids: Vec<u64> = bit_ids.iter().copied().map(u64::from).collect();
                atom_to_bits.push(output_ids.clone());
                for bit_id in output_ids {
                    bit_info_map
                        .entry(bit_id)
                        .or_default()
                        .push((atom_index as u32, 0));
                    atoms_per_bit
                        .entry(bit_id)
                        .or_default()
                        .push(vec![atom_index as i32]);
                }
            }

            AdditionalOutput {
                atom_counts: (mask & 0b00001 != 0).then(|| vec![1; per_atom_bits.len()]),
                atom_to_bits: (mask & 0b00010 != 0).then_some(atom_to_bits),
                bit_info_map: (mask & 0b00100 != 0).then_some(bit_info_map),
                bit_paths: (mask & 0b01000 != 0).then(BTreeMap::new),
                atoms_per_bit: (mask & 0b10000 != 0).then_some(atoms_per_bit),
            }
        };

        let mut matrix_calls = 0;
        let mut matrix_successes = 0;
        let mut matrix_errors = 0;
        let mut absent_output_calls = 0;
        let mut absent_output_successes = 0;
        for (size_index, &fp_size) in FP_SIZES.iter().enumerate() {
            for count_simulation in [false, true] {
                let count_bounds = if count_simulation {
                    SIMULATION_BOUNDS[size_index].to_vec()
                } else {
                    Vec::new()
                };
                let generator = get_morgan_generator(&MorganParams {
                    radius: 0,
                    fp_size,
                    count_simulation,
                    count_bounds,
                    ..MorganParams::default()
                })
                .expect("fixed O04 generator options satisfy constructor preconditions");

                for output_mask in 0_u8..32 {
                    let mut actual_output = make_initial_output(output_mask);
                    let result = get_fingerprint(
                        &record.topology,
                        &record.properties,
                        &valence,
                        &rings,
                        &generator,
                        &collision_arguments,
                        Some(&mut actual_output),
                    );
                    matrix_calls += 1;

                    if count_simulation && size_index == 0 {
                        let error = result.expect_err(
                            "one bound at fpSize 1 fails before division or output reset",
                        );
                        assert!(
                            matches!(
                                error,
                                MorganError::Fingerprint(FingerprintError::InvalidArguments {
                                    reason: "Count bounds size is >= fingerprint size"
                                })
                            ),
                            "pinned size-bound ValueError for mask={output_mask:#07b}: {error}"
                        );
                        assert_eq!(actual_output, make_initial_output(output_mask));
                        matrix_errors += 1;
                    } else {
                        let actual = result.expect("fixed O04 bit IDs fit the dense result");
                        let expected_bits = if count_simulation {
                            SIMULATION_BITS[size_index]
                        } else {
                            NO_SIMULATION_BITS[size_index]
                        };
                        let expected_rows = if count_simulation {
                            &SIMULATION_ROWS[size_index]
                        } else {
                            &NO_SIMULATION_ROWS[size_index]
                        };
                        assert_eq!(actual.n_bits(), fp_size);
                        assert_eq!(
                            actual.on_bits().as_slice(),
                            expected_bits,
                            "fixed dense bits for fpSize={fp_size}, countSimulation={count_simulation}, mask={output_mask:#07b}"
                        );
                        assert_eq!(
                            actual_output,
                            make_expected_output(output_mask, expected_rows),
                            "complete AdditionalOutput for fpSize={fp_size}, countSimulation={count_simulation}, mask={output_mask:#07b}"
                        );
                        matrix_successes += 1;
                    }
                }

                let result = get_fingerprint(
                    &record.topology,
                    &record.properties,
                    &valence,
                    &rings,
                    &generator,
                    &collision_arguments,
                    None,
                );
                absent_output_calls += 1;
                if count_simulation && size_index == 0 {
                    assert!(matches!(
                        result,
                        Err(MorganError::Fingerprint(
                            FingerprintError::InvalidArguments {
                                reason: "Count bounds size is >= fingerprint size"
                            }
                        ))
                    ));
                } else {
                    let actual = result.expect("valid options work with a null output pointer");
                    let expected_bits = if count_simulation {
                        SIMULATION_BITS[size_index]
                    } else {
                        NO_SIMULATION_BITS[size_index]
                    };
                    assert_eq!(actual.n_bits(), fp_size);
                    assert_eq!(actual.on_bits().as_slice(), expected_bits);
                    absent_output_successes += 1;
                }
            }
        }
        assert_eq!(matrix_calls, 256);
        assert_eq!(matrix_successes, 224);
        assert_eq!(matrix_errors, 32);
        assert_eq!(absent_output_calls, 8);
        assert_eq!(absent_output_successes, 7);

        let empty_record = parse_smiles("", &SmilesParseParams::default())
            .expect("the fixed empty SMILES input parses to the empty molecule");
        let empty_valence = assign_valence(&empty_record.topology, &ValenceParams::default())
            .expect("the empty molecule has an empty valid valence assignment");
        let empty_rings = fast_find_rings(&empty_record.topology)
            .expect("the empty molecule has an empty ring assignment");
        let zero_size_generator = get_morgan_generator(&MorganParams {
            radius: 0,
            fp_size: 0,
            count_simulation: false,
            count_bounds: Vec::new(),
            ..MorganParams::default()
        })
        .expect("zero-size non-simulation options satisfy source construction");
        let mut zero_size_empty_calls = 0;
        for output_mask in 0_u8..32 {
            let mut actual_output = make_initial_output(output_mask);
            let result = get_fingerprint(
                &empty_record.topology,
                &empty_record.properties,
                &empty_valence,
                &empty_rings,
                &zero_size_generator,
                &default_arguments,
                Some(&mut actual_output),
            )
            .expect("zero-size dense output succeeds for an empty environment map");
            assert_eq!(result.n_bits(), 0);
            assert!(result.on_bits().is_empty());
            assert_eq!(actual_output, make_expected_output(output_mask, &[]));
            zero_size_empty_calls += 1;
        }
        assert_eq!(zero_size_empty_calls, 32);

        let zero_size_arguments =
            FingerprintFuncArguments::new(None, None, Some(&CUSTOM_ATOM_INVARIANTS), None, -1);
        let mut zero_size_nonempty_errors = 0;
        for output_mask in 0_u8..32 {
            let mut actual_output = make_initial_output(output_mask);
            let error = get_fingerprint(
                &record.topology,
                &record.properties,
                &valence,
                &rings,
                &zero_size_generator,
                &zero_size_arguments,
                Some(&mut actual_output),
            )
            .expect_err("the first raw source ID fails the zero-length dense bit set");
            assert!(matches!(
                error,
                MorganError::Fingerprint(FingerprintError::SparseIndexOutOfRange {
                    index: 1,
                    size: 0,
                })
            ));
            assert_eq!(
                actual_output,
                make_expected_output(output_mask, &ZERO_SIZE_RAW_ROWS),
                "helper output is complete before dense projection errors, mask={output_mask:#07b}"
            );
            zero_size_nonempty_errors += 1;
        }
        assert_eq!(zero_size_nonempty_errors, 32);

        let mut bound_error_calls = 0;
        for (fp_size, count_bounds, clear_bounds, expected_reason) in [
            (2, &[1][..], true, "Count bounds are empty"),
            (
                2,
                &[1, 2][..],
                false,
                "Count bounds size is >= fingerprint size",
            ),
            (0, &[1][..], true, "Count bounds are empty"),
            (
                0,
                &[1][..],
                false,
                "Count bounds size is >= fingerprint size",
            ),
        ] {
            let mut generator = get_morgan_generator(&MorganParams {
                radius: 0,
                fp_size,
                count_simulation: true,
                count_bounds: count_bounds.to_vec(),
                ..MorganParams::default()
            })
            .expect("the source constructor allows these nonempty count bounds");
            if clear_bounds {
                // RDKit's mutable getOptions() exposes this exact per-call
                // defensive branch after valid generator construction.
                generator.fingerprint_arguments.count_bounds.clear();
            }

            for output_mask in 0_u8..32 {
                let mut actual_output = make_initial_output(output_mask);
                let error = get_fingerprint(
                    &record.topology,
                    &record.properties,
                    &valence,
                    &rings,
                    &generator,
                    &collision_arguments,
                    Some(&mut actual_output),
                )
                .expect_err("invalid count-bound state is rejected before output setup");
                assert!(
                    matches!(
                        error,
                        MorganError::Fingerprint(FingerprintError::InvalidArguments { reason })
                            if reason == expected_reason
                    ),
                    "source-ordered bound validation for size={fp_size}, mask={output_mask:#07b}: {error}"
                );
                assert_eq!(actual_output, make_initial_output(output_mask));
                bound_error_calls += 1;
            }
        }
        assert_eq!(bound_error_calls, 128);

        let mut cco_reference_calls = 0;
        for (size_index, &fp_size) in FP_SIZES.iter().enumerate() {
            for count_simulation in [false, true] {
                if count_simulation && size_index == 0 {
                    continue;
                }
                let generator = get_morgan_generator(&MorganParams {
                    radius: 0,
                    fp_size,
                    count_simulation,
                    count_bounds: if count_simulation {
                        SIMULATION_BOUNDS[size_index].to_vec()
                    } else {
                        Vec::new()
                    },
                    ..MorganParams::default()
                })
                .expect("fixed CCO reference options satisfy source construction");
                let actual = get_fingerprint(
                    &record.topology,
                    &record.properties,
                    &valence,
                    &rings,
                    &generator,
                    &default_arguments,
                    None,
                )
                .expect("fixed CCO radius-zero reference bits fit configured dense output");
                let expected_bits = if count_simulation {
                    CCO_SIMULATION_BITS[size_index]
                } else {
                    CCO_NO_SIMULATION_BITS[size_index]
                };
                assert_eq!(actual.n_bits(), fp_size);
                assert_eq!(actual.on_bits().as_slice(), expected_bits);
                cco_reference_calls += 1;
            }
        }
        assert_eq!(cco_reference_calls, 7);
        assert_eq!(record.topology, original_topology);
        assert_eq!(record.properties, original_properties);
    }
}
