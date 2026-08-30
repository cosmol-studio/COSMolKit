use cosmolkit_experiment_batch::{map_ordered, parse_smiles_batch};
use cosmolkit_experiment_chemistry::{
    add_hydrogens, enumerate_stereoisomers, enumerate_tautomers, generate_conformer,
    remove_hydrogens,
};
use cosmolkit_experiment_io::{parse_smarts, parse_smiles, read_pdb, read_sdf, read_xyz};
use cosmolkit_experiment_macros::{mol_multi_op_body, mol_op_body, molecule_ops};
use cosmolkit_experiment_model::{
    Conformer3D, CoordinateBlock, MoleculeProperties, QueryAst, TopologyBlock,
};
use indicatif::ProgressBar;
use std::sync::atomic::{AtomicU64, Ordering};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OperationOutput {
    Single,
    Multiple,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AccessMode {
    None,
    Read,
    Write,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BlockSet(u8);

impl BlockSet {
    pub const NONE: Self = Self(0);
    pub const TOPOLOGY: Self = Self(1 << 0);
    pub const COORDINATES: Self = Self(1 << 1);
    pub const PROPERTIES: Self = Self(1 << 2);

    #[must_use]
    pub const fn union(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }

    #[must_use]
    pub const fn contains(self, other: Self) -> bool {
        (self.0 & other.0) == other.0
    }

    #[must_use]
    pub const fn intersects(self, other: Self) -> bool {
        (self.0 & other.0) != 0
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BlockAccess {
    pub topology: AccessMode,
    pub coordinates: AccessMode,
    pub properties: AccessMode,
}

impl BlockAccess {
    #[must_use]
    pub const fn read(self) -> BlockSet {
        let mut result = BlockSet::NONE;
        if matches!(self.topology, AccessMode::Read) {
            result = result.union(BlockSet::TOPOLOGY);
        }
        if matches!(self.coordinates, AccessMode::Read) {
            result = result.union(BlockSet::COORDINATES);
        }
        if matches!(self.properties, AccessMode::Read) {
            result = result.union(BlockSet::PROPERTIES);
        }
        result
    }

    #[must_use]
    pub const fn write(self) -> BlockSet {
        let mut result = BlockSet::NONE;
        if matches!(self.topology, AccessMode::Write) {
            result = result.union(BlockSet::TOPOLOGY);
        }
        if matches!(self.coordinates, AccessMode::Write) {
            result = result.union(BlockSet::COORDINATES);
        }
        if matches!(self.properties, AccessMode::Write) {
            result = result.union(BlockSet::PROPERTIES);
        }
        result
    }

    #[must_use]
    pub const fn can_read(self, block: BlockSet) -> bool {
        self.read().union(self.write()).contains(block)
    }

    #[must_use]
    pub const fn can_write(self, block: BlockSet) -> bool {
        self.write().contains(block)
    }

    #[must_use]
    pub const fn has_overlapping_read_write(self) -> bool {
        self.read().intersects(self.write())
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TopologyEditKind {
    None,
    Local,
    Compacting,
    Appending,
    Renumbering,
    Merge,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OperationDomain {
    Topology,
    Coordinate,
    Construction,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ParityPolicy {
    NotApplicable,
    RequiredWhenSupported,
    RequiredNow,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CipStatePolicy {
    Preserve,
    ClearComputed,
    Assign,
    TautomerSourceTransition,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MappingRequirement {
    None,
    Identity,
    Required,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SemanticPreconditionSet(u8);

impl SemanticPreconditionSet {
    pub const NONE: Self = Self(0);
    pub const NON_EMPTY_INPUT: Self = Self(1 << 0);
    pub const TRUSTED_TOPOLOGY: Self = Self(1 << 1);

    #[must_use]
    pub const fn contains(self, other: Self) -> bool {
        (self.0 & other.0) == other.0
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct DerivedEffects {
    pub recompute: &'static [&'static str],
    pub preserve: &'static [&'static str],
    pub invalidate: &'static [&'static str],
    pub operation_defined: &'static [&'static str],
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MoleculeOpSpec {
    pub method: &'static str,
    pub impl_fn: &'static str,
    pub domain: &'static str,
    pub kind: &'static str,
    pub output: OperationOutput,
    pub result_type: Option<&'static str>,
    pub access: BlockAccess,
    pub may_mutate: BlockSet,
    pub topology_edit: Option<&'static str>,
    pub auto_remap: BlockSet,
    pub derived_effects: DerivedEffects,
    pub semantic_preconditions: &'static [&'static str],
    pub requires_mapping: bool,
    pub support: &'static str,
    pub parity: &'static str,
    pub parity_profile: Option<&'static str>,
    pub io_roundtrip: &'static str,
    pub invariant_profile: &'static str,
    pub feature: &'static str,
    pub mapping_requirement: MappingRequirement,
    pub cip_state: CipStatePolicy,
}

molecule_ops! {
    op add_hydrogens() {
        method: with_hydrogens,
        impl_fn: add_hydrogens_impl,
        output: single,
        domain: topology,
        kind: strong,
        topology_edit: appending,
        access: { read: [], write: [topology, coordinates, properties] },
        may_mutate: [topology, coordinates, properties],
        auto_remap: [coordinates, properties],
        derived_effects: {
            recompute: [],
            preserve: [],
            invalidate: [rings, valence, stereo],
        },
        cip_state: clear_computed,
        semantic_preconditions: [],
        requires_mapping: required,
        feature: HYDROGENS_FEATURE,
        parity: required_now,
        parity_profile: "experiment",
        io_roundtrip: false,
        invariant_profile: "strong_topology",
        inplace: true,
        inplace_method: add_hydrogens_,
    }

    op remove_hydrogens() {
        method: without_hydrogens,
        impl_fn: remove_hydrogens_impl,
        output: single,
        domain: topology,
        kind: strong,
        topology_edit: compacting,
        access: { read: [], write: [topology, coordinates, properties] },
        may_mutate: [topology, coordinates, properties],
        auto_remap: [coordinates, properties],
        derived_effects: {
            recompute: [],
            preserve: [],
            invalidate: [rings, stereo],
            operation_defined: [valence],
        },
        cip_state: clear_computed,
        semantic_preconditions: [],
        requires_mapping: required,
        feature: HYDROGENS_FEATURE,
        parity: required_now,
        parity_profile: "experiment",
        io_roundtrip: false,
        invariant_profile: "strong_topology",
        inplace: true,
        inplace_method: remove_hydrogens_,
    }

    op tautomers() {
        method: tautomers,
        impl_fn: tautomers_impl,
        output: multiple,
        domain: topology,
        kind: strong,
        topology_edit: local,
        access: { read: [], write: [topology, coordinates, properties] },
        may_mutate: [topology, coordinates, properties],
        auto_remap: [],
        derived_effects: {
            recompute: [],
            preserve: [],
            invalidate: [stereo],
        },
        cip_state: tautomer_source_transition,
        semantic_preconditions: [],
        requires_mapping: none,
        feature: TAUTOMER_FEATURE,
        parity: required_now,
        parity_profile: "experiment",
        io_roundtrip: false,
        invariant_profile: "multiple_output",
        inplace: false,
    }

    op stereoisomers() {
        method: stereoisomers,
        impl_fn: stereoisomers_impl,
        output: multiple,
        domain: topology,
        kind: strong,
        topology_edit: local,
        access: { read: [], write: [topology, coordinates, properties] },
        may_mutate: [topology, coordinates, properties],
        auto_remap: [],
        derived_effects: {
            recompute: [],
            preserve: [],
            invalidate: [stereo],
        },
        cip_state: assign,
        semantic_preconditions: [],
        requires_mapping: none,
        feature: STEREOISOMER_FEATURE,
        parity: required_now,
        parity_profile: "experiment",
        io_roundtrip: false,
        invariant_profile: "multiple_output",
        inplace: false,
    }

    op conformer() {
        method: with_conformer,
        impl_fn: conformer_impl,
        output: single,
        domain: coordinate,
        kind: coordinate,
        topology_edit: none,
        access: { read: [topology, properties], write: [coordinates] },
        may_mutate: [coordinates],
        auto_remap: [],
        derived_effects: {
            recompute: [],
            preserve: [topology],
            invalidate: [coordinate_dependent_stereo],
        },
        cip_state: preserve,
        semantic_preconditions: [],
        requires_mapping: none,
        feature: CONFORMER_FEATURE,
        parity: required_now,
        parity_profile: "experiment",
        io_roundtrip: false,
        invariant_profile: "coordinates",
        inplace: false,
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum OperationError {
    InvalidInput {
        operation: &'static str,
        message: &'static str,
    },
    Precondition {
        operation: &'static str,
        requirement: &'static str,
    },
    AccessDenied {
        operation: &'static str,
        block: &'static str,
    },
    BlockCheckedOut {
        block: &'static str,
    },
    IncompleteCommit {
        block: &'static str,
    },
    ContractViolation(&'static str),
    InvariantViolation(&'static str),
    InvalidTopology(String),
    Algorithm(String),
}

impl std::fmt::Display for OperationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidInput { operation, message } => {
                write!(f, "{operation}: invalid input: {message}")
            }
            Self::Precondition {
                operation,
                requirement,
            } => write!(f, "{operation}: precondition failed: {requirement}"),
            Self::AccessDenied { operation, block } => {
                write!(f, "{operation} cannot access {block}")
            }
            Self::BlockCheckedOut { block } => write!(f, "{block} is already checked out"),
            Self::IncompleteCommit { block } => write!(f, "{block} was not committed"),
            Self::ContractViolation(message) => f.write_str(message),
            Self::InvariantViolation(message) => f.write_str(message),
            Self::InvalidTopology(message) | Self::Algorithm(message) => f.write_str(message),
        }
    }
}

impl std::error::Error for OperationError {}

#[derive(Clone, Debug, PartialEq)]
pub struct Molecule {
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
}

impl Molecule {
    fn from_parts(
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
    ) -> Self {
        Self {
            topology,
            coordinates,
            properties,
        }
    }

    #[must_use]
    pub fn from_smiles(source: &str) -> Result<Self, OperationError> {
        validate_input("parse_smiles", source)?;
        parse_smiles(source)
            .map(|(topology, coordinates, properties)| {
                Self::from_parts(topology, coordinates, properties)
            })
            .map_err(OperationError::Algorithm)
    }

    #[must_use]
    pub fn from_sdf(source: &str) -> Result<Vec<Self>, OperationError> {
        validate_input("read_sdf", source)?;
        read_sdf(source)
            .map(|parts| {
                parts
                    .into_iter()
                    .map(|(topology, coordinates, properties)| {
                        Self::from_parts(topology, coordinates, properties)
                    })
                    .collect()
            })
            .map_err(OperationError::Algorithm)
    }

    #[must_use]
    pub fn from_pdb(source: &str) -> Result<Self, OperationError> {
        validate_input("read_pdb", source)?;
        read_pdb(source)
            .map(|(topology, coordinates, properties)| {
                Self::from_parts(topology, coordinates, properties)
            })
            .map_err(OperationError::Algorithm)
    }

    #[must_use]
    pub fn from_xyz(source: &str) -> Result<Self, OperationError> {
        validate_input("read_xyz", source)?;
        read_xyz(source)
            .map(|(topology, coordinates, properties)| {
                Self::from_parts(topology, coordinates, properties)
            })
            .map_err(OperationError::Algorithm)
    }

    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.topology.atom_count()
    }

    #[must_use]
    pub fn coordinates_3d(&self) -> &[Conformer3D] {
        self.coordinates.conformers_3d()
    }

    #[must_use]
    pub fn properties(&self) -> &MoleculeProperties {
        &self.properties
    }
}

pub fn molecule_from_smarts(source: &str) -> Result<QueryAst, OperationError> {
    validate_input("parse_smarts", source)?;
    parse_smarts(source).map_err(OperationError::Algorithm)
}

#[mol_op_body(add_hydrogens, parts)]
fn add_hydrogens_impl() -> Result<(), OperationError> {
    let (topology, coordinates, properties) = parts.extract_all_writable()?;
    let (topology, coordinates, properties) =
        add_hydrogens(topology, coordinates, properties).map_err(OperationError::Algorithm)?;
    parts.commit_topology(topology)?;
    parts.commit_coordinates(coordinates)?;
    parts.commit_properties(properties)?;
    parts.record_topology_mapping()?;
    Ok(())
}

#[mol_op_body(remove_hydrogens, parts)]
fn remove_hydrogens_impl() -> Result<(), OperationError> {
    let (topology, coordinates, properties) = parts.extract_all_writable()?;
    let (topology, coordinates, properties) =
        remove_hydrogens(topology, coordinates, properties).map_err(OperationError::Algorithm)?;
    parts.commit_topology(topology)?;
    parts.commit_coordinates(coordinates)?;
    parts.commit_properties(properties)?;
    parts.record_topology_mapping()?;
    Ok(())
}

#[mol_multi_op_body(tautomers, parts)]
fn tautomers_impl() -> Result<(), OperationError> {
    let candidates = parts.with_source_read_parts(|read| {
        enumerate_tautomers(read.topology()?, read.coordinates()?, read.properties()?)
            .map_err(OperationError::Algorithm)
    })?;
    for (topology, coordinates, properties) in candidates {
        let branch = parts.derive_from_source(move |branch| {
            branch.replace_all_writable(topology, coordinates, properties)
        })?;
        parts.emit(branch)?;
    }
    Ok(())
}

#[mol_multi_op_body(stereoisomers, parts)]
fn stereoisomers_impl() -> Result<(), OperationError> {
    let candidates = parts.with_source_read_parts(|read| {
        enumerate_stereoisomers(read.topology()?, read.coordinates()?, read.properties()?)
            .map_err(OperationError::Algorithm)
    })?;
    for (topology, coordinates, properties) in candidates {
        let branch = parts.derive_from_source(move |branch| {
            branch.replace_all_writable(topology, coordinates, properties)
        })?;
        parts.emit(branch)?;
    }
    Ok(())
}

#[mol_op_body(conformer, parts)]
fn conformer_impl() -> Result<(), OperationError> {
    let conformer =
        generate_conformer(parts.begin_topology_read()?).map_err(OperationError::Algorithm)?;
    let mut coordinates = parts.begin_coordinates_mut()?;
    coordinates.add_conformer_3d(conformer);
    parts.commit_coordinates(coordinates)
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct MoleculeBatch {
    records: Vec<Molecule>,
    n_jobs: Option<usize>,
    progress_bar: Option<bool>,
}

impl MoleculeBatch {
    pub fn from_smiles<S>(sources: impl IntoIterator<Item = S>) -> Result<Self, OperationError>
    where
        S: AsRef<str> + Send,
    {
        Self::from_smiles_with_options(sources, None, None)
    }

    pub fn from_smiles_with_options<S>(
        sources: impl IntoIterator<Item = S>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Self, OperationError>
    where
        S: AsRef<str> + Send,
    {
        let n_jobs = n_jobs.map(|value| value.max(1));
        let sources: Vec<S> = sources.into_iter().collect();
        let records = with_batch_progress(progress_bar, sources.len(), |progress| {
            parse_smiles_batch(sources, n_jobs, progress)
        })
        .map_err(OperationError::Algorithm)?;
        Ok(Self {
            records: records
                .into_iter()
                .map(|record| {
                    Molecule::from_parts(record.topology, record.coordinates, record.properties)
                })
                .collect(),
            n_jobs,
            progress_bar,
        })
    }

    #[must_use]
    pub fn with_parallel_jobs(mut self, n_jobs: Option<usize>) -> Self {
        self.n_jobs = n_jobs.map(|value| value.max(1));
        self
    }

    #[must_use]
    pub fn parallel_jobs(&self) -> Option<usize> {
        self.n_jobs
    }

    #[must_use]
    pub fn with_progress_bar(mut self, progress_bar: Option<bool>) -> Self {
        self.progress_bar = progress_bar;
        self
    }

    #[must_use]
    pub fn progress_bar(&self) -> Option<bool> {
        self.progress_bar
    }

    pub fn add_hydrogens(&self) -> Result<Self, OperationError> {
        self.add_hydrogens_with_options(None, None)
    }

    pub fn add_hydrogens_with_options(
        &self,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Self, OperationError> {
        let effective_n_jobs = n_jobs.or(self.n_jobs).map(|value| value.max(1));
        let effective_progress_bar = progress_bar.or(self.progress_bar);
        let records =
            with_batch_progress(effective_progress_bar, self.records.len(), |progress| {
                map_ordered(
                    self.records.iter(),
                    effective_n_jobs,
                    progress,
                    Molecule::with_hydrogens,
                )
            })?;
        Ok(Self {
            records,
            n_jobs: self.n_jobs,
            progress_bar: self.progress_bar,
        })
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.records.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
}

fn with_batch_progress<R>(
    enabled: Option<bool>,
    total: usize,
    operation: impl FnOnce(cosmolkit_experiment_batch::BatchProgress<'_>) -> R,
) -> R {
    if enabled.unwrap_or(false) {
        let progress_bar = ProgressBar::new(total as u64);
        let tick = || progress_bar.inc(1);
        let result = operation(Some(&tick));
        progress_bar.finish();
        result
    } else {
        operation(None)
    }
}

pub(crate) struct OpParts<'a> {
    spec: &'static MoleculeOpSpec,
    topology: Option<TopologyBlock>,
    coordinates: Option<CoordinateBlock>,
    properties: Option<MoleculeProperties>,
    topology_checked_out: bool,
    coordinates_checked_out: bool,
    properties_checked_out: bool,
    topology_mapping_recorded: bool,
    in_place_target: Option<&'a mut Molecule>,
}

pub(crate) struct MoleculeReadParts<'a> {
    spec: &'static MoleculeOpSpec,
    molecule: &'a Molecule,
}

impl<'a> MoleculeReadParts<'a> {
    fn topology(&self) -> Result<&'a TopologyBlock, OperationError> {
        if !self.spec.access.can_read(BlockSet::TOPOLOGY) {
            return Err(OperationError::AccessDenied {
                operation: self.spec.method,
                block: "topology",
            });
        }
        Ok(&self.molecule.topology)
    }

    fn coordinates(&self) -> Result<&'a CoordinateBlock, OperationError> {
        if !self.spec.access.can_read(BlockSet::COORDINATES) {
            return Err(OperationError::AccessDenied {
                operation: self.spec.method,
                block: "coordinates",
            });
        }
        Ok(&self.molecule.coordinates)
    }

    fn properties(&self) -> Result<&'a MoleculeProperties, OperationError> {
        if !self.spec.access.can_read(BlockSet::PROPERTIES) {
            return Err(OperationError::AccessDenied {
                operation: self.spec.method,
                block: "properties",
            });
        }
        Ok(&self.molecule.properties)
    }
}

impl<'a> OpParts<'a> {
    fn new(source: &'a Molecule, spec: &'static MoleculeOpSpec) -> Result<Self, OperationError> {
        if spec.output != OperationOutput::Single {
            return Err(OperationError::ContractViolation(
                "single-output OpParts used with multiple-output operation",
            ));
        }
        Ok(Self {
            spec,
            topology: Some(source.topology.clone()),
            coordinates: Some(source.coordinates.clone()),
            properties: Some(source.properties.clone()),
            topology_checked_out: false,
            coordinates_checked_out: false,
            properties_checked_out: false,
            topology_mapping_recorded: false,
            in_place_target: None,
        })
    }

    fn new_multi_branch(
        source: &'a Molecule,
        spec: &'static MoleculeOpSpec,
    ) -> Result<Self, OperationError> {
        if spec.output != OperationOutput::Multiple {
            return Err(OperationError::ContractViolation(
                "multiple-output branch used with single-output operation",
            ));
        }
        Ok(Self {
            spec,
            topology: Some(source.topology.clone()),
            coordinates: Some(source.coordinates.clone()),
            properties: Some(source.properties.clone()),
            topology_checked_out: false,
            coordinates_checked_out: false,
            properties_checked_out: false,
            topology_mapping_recorded: false,
            in_place_target: None,
        })
    }

    fn new_in_place(
        target: &'a mut Molecule,
        spec: &'static MoleculeOpSpec,
    ) -> Result<Self, OperationError> {
        if spec.output != OperationOutput::Single {
            return Err(OperationError::ContractViolation(
                "multiple-output operations cannot run in place",
            ));
        }
        Ok(Self {
            spec,
            topology: Some(target.topology.clone()),
            coordinates: Some(target.coordinates.clone()),
            properties: Some(target.properties.clone()),
            topology_checked_out: false,
            coordinates_checked_out: false,
            properties_checked_out: false,
            topology_mapping_recorded: false,
            in_place_target: Some(target),
        })
    }

    fn ensure_access(&self, mode: AccessMode, block: &'static str) -> Result<(), OperationError> {
        if mode == AccessMode::None {
            return Err(OperationError::AccessDenied {
                operation: self.spec.method,
                block,
            });
        }
        Ok(())
    }

    fn begin_topology_read(&self) -> Result<&TopologyBlock, OperationError> {
        self.ensure_access(self.spec.access.topology, "topology")?;
        self.topology
            .as_ref()
            .ok_or(OperationError::BlockCheckedOut { block: "topology" })
    }

    fn begin_topology_mut(&mut self) -> Result<TopologyBlock, OperationError> {
        self.ensure_access(self.spec.access.topology, "topology")?;
        let topology = self
            .topology
            .take()
            .ok_or(OperationError::BlockCheckedOut { block: "topology" })?;
        self.topology_checked_out = true;
        Ok(topology)
    }

    fn commit_topology(&mut self, topology: TopologyBlock) -> Result<(), OperationError> {
        topology
            .validate()
            .map_err(|error| OperationError::InvalidTopology(error.to_string()))?;
        self.topology = Some(topology);
        self.topology_checked_out = false;
        Ok(())
    }

    fn begin_coordinates_mut(&mut self) -> Result<CoordinateBlock, OperationError> {
        self.ensure_access(self.spec.access.coordinates, "coordinates")?;
        let coordinates = self
            .coordinates
            .take()
            .ok_or(OperationError::BlockCheckedOut {
                block: "coordinates",
            })?;
        self.coordinates_checked_out = true;
        Ok(coordinates)
    }

    fn commit_coordinates(&mut self, coordinates: CoordinateBlock) -> Result<(), OperationError> {
        self.coordinates = Some(coordinates);
        self.coordinates_checked_out = false;
        Ok(())
    }

    fn begin_properties_mut(&mut self) -> Result<MoleculeProperties, OperationError> {
        self.ensure_access(self.spec.access.properties, "properties")?;
        let properties = self
            .properties
            .take()
            .ok_or(OperationError::BlockCheckedOut {
                block: "properties",
            })?;
        self.properties_checked_out = true;
        Ok(properties)
    }

    fn commit_properties(&mut self, properties: MoleculeProperties) -> Result<(), OperationError> {
        self.properties = Some(properties);
        self.properties_checked_out = false;
        Ok(())
    }

    fn extract_all_writable(
        &mut self,
    ) -> Result<(TopologyBlock, CoordinateBlock, MoleculeProperties), OperationError> {
        let topology = self.begin_topology_mut()?;
        let coordinates = self.begin_coordinates_mut()?;
        let properties = self.begin_properties_mut()?;
        Ok((topology, coordinates, properties))
    }

    fn replace_all_writable(
        &mut self,
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
    ) -> Result<(), OperationError> {
        if !self.spec.access.can_write(BlockSet::TOPOLOGY)
            || !self.spec.access.can_write(BlockSet::COORDINATES)
            || !self.spec.access.can_write(BlockSet::PROPERTIES)
        {
            return Err(OperationError::AccessDenied {
                operation: self.spec.method,
                block: "all writable blocks",
            });
        }
        self.topology = Some(topology);
        self.coordinates = Some(coordinates);
        self.properties = Some(properties);
        self.topology_checked_out = false;
        self.coordinates_checked_out = false;
        self.properties_checked_out = false;
        Ok(())
    }

    fn record_topology_mapping(&mut self) -> Result<(), OperationError> {
        if self.spec.requires_mapping {
            self.topology_mapping_recorded = true;
            Ok(())
        } else {
            Err(OperationError::ContractViolation(
                "operation does not require topology mapping",
            ))
        }
    }

    fn finish(self) -> Result<Molecule, OperationError> {
        if self.topology_checked_out {
            return Err(OperationError::IncompleteCommit { block: "topology" });
        }
        if self.coordinates_checked_out {
            return Err(OperationError::IncompleteCommit {
                block: "coordinates",
            });
        }
        if self.properties_checked_out {
            return Err(OperationError::IncompleteCommit {
                block: "properties",
            });
        }
        if self.spec.may_mutate != self.spec.access.write() && self.spec.kind != "construction" {
            return Err(OperationError::ContractViolation(
                "may_mutate does not match declared topology write access",
            ));
        }
        if self.spec.requires_mapping && !self.topology_mapping_recorded {
            return Err(OperationError::ContractViolation(
                "strong topology operation is missing mapping",
            ));
        }
        let topology = self
            .topology
            .ok_or(OperationError::IncompleteCommit { block: "topology" })?;
        let coordinates = self.coordinates.ok_or(OperationError::IncompleteCommit {
            block: "coordinates",
        })?;
        let properties = self.properties.ok_or(OperationError::IncompleteCommit {
            block: "properties",
        })?;
        topology
            .validate()
            .map_err(|error| OperationError::InvalidTopology(error.to_string()))?;
        Ok(Molecule::from_parts(topology, coordinates, properties))
    }

    fn abort_in_place(&mut self) {}

    fn finish_in_place(self) -> Result<(), OperationError> {
        let mut parts = self;
        let target = parts
            .in_place_target
            .take()
            .ok_or(OperationError::ContractViolation(
                "non-in-place operation used in-place finish",
            ))?;
        let molecule = parts.finish()?;
        *target = molecule;
        Ok(())
    }
}

static NEXT_MULTI_OUTPUT_EXECUTOR_ID: AtomicU64 = AtomicU64::new(1);

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct MoleculeBranchId {
    executor_id: u64,
    index: usize,
}

pub(crate) struct MultiOutputOpParts<'a> {
    spec: &'static MoleculeOpSpec,
    source: &'a Molecule,
    executor_id: u64,
    candidates: Vec<Molecule>,
    emitted: Vec<MoleculeBranchId>,
}

impl<'a> MultiOutputOpParts<'a> {
    fn new(source: &'a Molecule, spec: &'static MoleculeOpSpec) -> Result<Self, OperationError> {
        if spec.output != OperationOutput::Multiple {
            return Err(OperationError::ContractViolation(
                "single-output operation used with multiple-output runtime",
            ));
        }
        let executor_id = NEXT_MULTI_OUTPUT_EXECUTOR_ID
            .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |id| id.checked_add(1))
            .map_err(|_| {
                OperationError::ContractViolation("multi-output executor identity exhausted")
            })?;
        Ok(Self {
            spec,
            source,
            executor_id,
            candidates: Vec::new(),
            emitted: Vec::new(),
        })
    }

    fn with_source_read_parts<R>(
        &self,
        read: impl FnOnce(MoleculeReadParts<'_>) -> Result<R, OperationError>,
    ) -> Result<R, OperationError> {
        read(MoleculeReadParts {
            spec: self.spec,
            molecule: self.source,
        })
    }

    fn derive_from_source(
        &mut self,
        derive: impl FnOnce(&mut OpParts<'_>) -> Result<(), OperationError>,
    ) -> Result<MoleculeBranchId, OperationError> {
        let mut parts = OpParts::new_multi_branch(self.source, self.spec)?;
        derive(&mut parts)?;
        let molecule = parts.finish()?;
        let branch = MoleculeBranchId {
            executor_id: self.executor_id,
            index: self.candidates.len(),
        };
        self.candidates.push(molecule);
        Ok(branch)
    }

    #[allow(dead_code)]
    fn derive_from_branch(
        &mut self,
        parent: MoleculeBranchId,
        derive: impl FnOnce(&mut OpParts<'_>) -> Result<(), OperationError>,
    ) -> Result<MoleculeBranchId, OperationError> {
        let source = self.branch(parent)?.clone();
        let mut parts = OpParts::new_multi_branch(&source, self.spec)?;
        derive(&mut parts)?;
        let molecule = parts.finish()?;
        let branch = MoleculeBranchId {
            executor_id: self.executor_id,
            index: self.candidates.len(),
        };
        self.candidates.push(molecule);
        Ok(branch)
    }

    fn branch(&self, branch: MoleculeBranchId) -> Result<&Molecule, OperationError> {
        if branch.executor_id != self.executor_id {
            return Err(OperationError::ContractViolation(
                "multi-output operation used a branch from another execution",
            ));
        }
        self.candidates
            .get(branch.index)
            .ok_or(OperationError::ContractViolation(
                "multi-output operation used an unknown branch",
            ))
    }

    #[allow(dead_code)]
    fn with_branch_read_parts<R>(
        &self,
        branch: MoleculeBranchId,
        read: impl FnOnce(MoleculeReadParts<'_>) -> Result<R, OperationError>,
    ) -> Result<R, OperationError> {
        read(MoleculeReadParts {
            spec: self.spec,
            molecule: self.branch(branch)?,
        })
    }

    fn emit(&mut self, branch: MoleculeBranchId) -> Result<(), OperationError> {
        self.branch(branch)?;
        self.emitted.push(branch);
        Ok(())
    }

    fn finish(self) -> Result<Vec<Molecule>, OperationError> {
        if self.emitted.is_empty() {
            return Err(OperationError::ContractViolation(
                "multiple-output operation produced no finalized branches",
            ));
        }
        let mut remaining = vec![0usize; self.candidates.len()];
        for branch in &self.emitted {
            if branch.executor_id != self.executor_id {
                return Err(OperationError::ContractViolation(
                    "multi-output operation emitted a branch from another execution",
                ));
            }
            let count =
                remaining
                    .get_mut(branch.index)
                    .ok_or(OperationError::ContractViolation(
                        "multi-output operation emitted unknown branch",
                    ))?;
            *count += 1;
        }

        let mut candidates = self.candidates.into_iter().map(Some).collect::<Vec<_>>();
        let mut outputs = Vec::with_capacity(self.emitted.len());
        for branch in self.emitted {
            remaining[branch.index] -= 1;
            if remaining[branch.index] == 0 {
                outputs.push(candidates[branch.index].take().ok_or(
                    OperationError::ContractViolation("multi-output branch was already consumed"),
                )?);
            } else {
                outputs.push(candidates[branch.index].as_ref().cloned().ok_or(
                    OperationError::ContractViolation("multi-output branch was already consumed"),
                )?);
            }
        }
        Ok(outputs)
    }
}

#[cfg(test)]
fn find_spec(method: &'static str) -> &'static MoleculeOpSpec {
    MOLECULE_OPS
        .iter()
        .find(|spec| spec.method == method)
        .expect("experiment registry must contain every wrapper operation")
}

fn validate_input(operation: &'static str, source: &str) -> Result<(), OperationError> {
    if source.is_empty() {
        return Err(OperationError::InvalidInput {
            operation,
            message: "input must not be empty",
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registry_contains_complete_contract_fields_for_every_experiment_surface() {
        assert_eq!(MOLECULE_OPS.len(), 5);
        assert_eq!(SUPPORT_MATRIX.len(), MOLECULE_OPS.len());
        assert_eq!(OPERATION_INVARIANT_MATRIX.len(), MOLECULE_OPS.len());
        assert_eq!(PARITY_MATRIX.len(), MOLECULE_OPS.len());
        for spec in MOLECULE_OPS {
            assert!(!spec.method.is_empty());
            assert!(!spec.impl_fn.is_empty());
            assert!(!spec.domain.is_empty());
            assert!(!spec.kind.is_empty());
            assert!(!spec.support.is_empty());
            assert!(!spec.parity.is_empty());
            assert!(!spec.invariant_profile.is_empty());
            assert!(!spec.access.has_overlapping_read_write());
            assert_eq!(spec.may_mutate, spec.access.write());
            assert_eq!(
                spec.requires_mapping,
                matches!(spec.mapping_requirement, MappingRequirement::Required)
            );
            if !spec.derived_effects.operation_defined.is_empty() {
                assert_eq!(spec.method, "without_hydrogens");
                assert_eq!(spec.derived_effects.operation_defined, &["valence"]);
            }
            if spec.parity != "not_applicable" {
                assert!(spec.parity_profile.is_some());
            }
        }
    }

    #[test]
    fn runtime_composes_all_experiment_surfaces_through_ops() {
        let molecule = Molecule::from_smiles("CC").expect("SMILES");
        assert!(molecule.with_hydrogens().unwrap().num_atoms() > molecule.num_atoms());
        assert_eq!(
            molecule.without_hydrogens().unwrap().num_atoms(),
            molecule.num_atoms()
        );
        assert_eq!(molecule.tautomers().unwrap().len(), 1);
        assert_eq!(molecule.stereoisomers().unwrap().len(), 1);
        assert_eq!(molecule.with_conformer().unwrap().coordinates_3d().len(), 1);
        assert_eq!(molecule_from_smarts("[#6]").unwrap().source(), "[#6]");
        assert_eq!(Molecule::from_sdf("sdf").unwrap().len(), 1);
        assert_eq!(Molecule::from_pdb("pdb").unwrap().num_atoms(), 1);
        assert_eq!(Molecule::from_xyz("xyz").unwrap().num_atoms(), 1);

        let batch = MoleculeBatch::from_smiles(["C", "CC"]).unwrap();
        assert_eq!(batch.len(), 2);
        assert_eq!(batch.add_hydrogens().unwrap().len(), 2);

        let mut in_place = molecule.clone();
        in_place.add_hydrogens_().unwrap();
        assert!(in_place.num_atoms() > molecule.num_atoms());
        in_place.remove_hydrogens_().unwrap();
    }

    #[test]
    fn batch_configuration_controls_ordered_parallel_execution() {
        let batch =
            MoleculeBatch::from_smiles_with_options(["C", "CC", "CCC"], Some(2), Some(false))
                .unwrap();
        assert_eq!(batch.parallel_jobs(), Some(2));
        assert_eq!(batch.progress_bar(), Some(false));
        assert_eq!(
            batch
                .records
                .iter()
                .map(Molecule::num_atoms)
                .collect::<Vec<_>>(),
            [1, 2, 3]
        );

        let transformed = batch
            .add_hydrogens_with_options(Some(1), Some(false))
            .unwrap();
        assert_eq!(transformed.parallel_jobs(), Some(2));
        assert_eq!(transformed.progress_bar(), Some(false));
        assert_eq!(
            transformed
                .records
                .iter()
                .map(Molecule::num_atoms)
                .collect::<Vec<_>>(),
            [2, 3, 4]
        );

        assert_eq!(
            MoleculeBatch::from_smiles_with_options(["C"], Some(0), Some(false))
                .unwrap()
                .parallel_jobs(),
            Some(1)
        );
    }

    #[test]
    fn op_parts_closes_each_block_begin_commit_pair_before_finish() {
        let molecule = Molecule::from_smiles("CC").expect("SMILES");
        let mut parts = OpParts::new(&molecule, find_spec("with_hydrogens")).unwrap();

        let topology = parts.begin_topology_mut().unwrap();
        parts.commit_topology(topology).unwrap();
        let coordinates = parts.begin_coordinates_mut().unwrap();
        parts.commit_coordinates(coordinates).unwrap();
        let properties = parts.begin_properties_mut().unwrap();
        parts.commit_properties(properties).unwrap();
        parts.record_topology_mapping().unwrap();

        let finished = parts.finish().unwrap();
        assert_eq!(finished.num_atoms(), molecule.num_atoms());
    }

    #[test]
    fn operation_declarations_and_bodies_use_the_proc_macro_surface() {
        let source = include_str!("lib.rs");
        assert!(source.contains("molecule_ops!"));
        assert!(source.contains("#[mol_op_body(add_hydrogens, parts)]"));
        assert!(source.contains("#[mol_multi_op_body(tautomers, parts)]"));
        assert_eq!(ADD_HYDROGENS_SPEC.impl_fn, "add_hydrogens_impl");
        assert_eq!(TAUTOMERS_SPEC.output, OperationOutput::Multiple);
    }

    #[test]
    fn multi_output_runtime_derives_reads_and_emits_validated_branches() {
        let molecule = Molecule::from_smiles("CC").unwrap();
        let mut outputs = MultiOutputOpParts::new(&molecule, &TAUTOMERS_SPEC).unwrap();
        let first = outputs.derive_from_source(|_| Ok(())).unwrap();
        assert_eq!(
            outputs
                .with_branch_read_parts(first, |read| Ok(read.topology()?.atom_count()))
                .unwrap(),
            2
        );
        let second = outputs.derive_from_branch(first, |_| Ok(())).unwrap();
        outputs.emit(first).unwrap();
        outputs.emit(second).unwrap();
        assert_eq!(outputs.finish().unwrap().len(), 2);

        let mut other = MultiOutputOpParts::new(&molecule, &TAUTOMERS_SPEC).unwrap();
        assert!(other.emit(first).is_err());
    }
}
