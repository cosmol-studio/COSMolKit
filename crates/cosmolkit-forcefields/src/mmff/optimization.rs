//! Canonical detached MMFF optimization; builder and shared kernel remain private.
//! RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8 (BSD license).
use super::{
    builder::{MmffBuilderError, MmffConformerContext, construct_force_field_with_props},
    mol_properties::{MmffMolProperties, MmffMolPropertiesError},
};
use crate::kernel::{ForceField, ForceFieldKernelError};
use crate::uff::convenience::{
    OptimizationOutcome, OptimizationStageError, SerialConformer, SerialConformerOptimizationError,
};
#[cfg(not(target_family = "wasm"))]
use crate::uff::convenience::{PreparedConformerDispatchOutcome, UffThreadCountError};
use cosmolkit_model::{CoordinateBlock, MoleculeProperties, TopologyBlock};

#[derive(Clone, Debug, PartialEq)]
pub struct MmffOptimizationParams {
    pub mmff_variant: String,
    pub max_iterations: i32,
    pub non_bonded_threshold: f64,
    pub conformer_id: Option<usize>,
    pub ignore_interfragment_interactions: bool,
}
impl Default for MmffOptimizationParams {
    fn default() -> Self {
        // RDKit❗✔️: int MMFFOptimizeMolecule(ROMol &mol, std::string mmffVariant = "MMFF94",
        // RDKit❗✔️:                          int maxIters = 200, double nonBondedThresh = 100.0,
        // RDKit❗✔️:                          int confId = -1,
        // RDKit❗✔️:                          bool ignoreInterfragInteractions = true) {
        // Behavior: exact Wrap/rdForceFields.cpp defaults, -1 represented by None.
        // Complexity: fixed parameter storage and one variant string allocation.
        Self {
            mmff_variant: "MMFF94".into(),
            max_iterations: 200,
            non_bonded_threshold: 100.,
            conformer_id: None,
            ignore_interfragment_interactions: true,
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct MmffConformerOptimizationParams {
    pub num_threads: i32,
    pub max_iterations: i32,
    pub mmff_variant: String,
    pub non_bonded_threshold: f64,
    pub ignore_interfragment_interactions: bool,
}
impl Default for MmffConformerOptimizationParams {
    fn default() -> Self {
        // RDKit❗✔️:                                       int numThreads = 1, int maxIters = 1000,
        // RDKit❗✔️:                                       std::string mmffVariant = "MMFF94",
        // RDKit❗✔️:                                       double nonBondedThresh = 10.0,
        // RDKit❗✔️:                                       bool ignoreInterfragInteractions = true) {
        // Behavior: MMFF.h multi-conformer defaults; no numerical work.
        // Complexity: constant storage plus one variant string allocation.
        Self {
            num_threads: 1,
            max_iterations: 1000,
            mmff_variant: "MMFF94".into(),
            non_bonded_threshold: 10.,
            ignore_interfragment_interactions: true,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MmffOptimizeMoleculeConfResult {
    pub needs_more: i32,
    pub energy: f64,
}
impl MmffOptimizeMoleculeConfResult {
    pub const fn needs_more(&self) -> bool {
        self.needs_more > 0
    }
    pub const fn status_code(&self) -> i32 {
        self.needs_more
    }
    pub const fn energy(&self) -> f64 {
        self.energy
    }
}
#[derive(Debug)]
pub struct MmffSingleOutcome {
    pub status: i32,
    pub aromatic_ring_count: usize,
    pub acquired_rings: Option<cosmolkit_core::RingInfo>,
}
#[derive(Debug)]
pub struct MmffConformerOutcomes {
    pub conformer_results: Vec<MmffOptimizeMoleculeConfResult>,
    pub aromatic_ring_count: usize,
    pub acquired_rings: Option<cosmolkit_core::RingInfo>,
}
#[derive(Debug, thiserror::Error)]
enum Failure {
    #[error(transparent)]
    Properties(#[from] MmffMolPropertiesError),
    #[error(transparent)]
    Builder(#[from] MmffBuilderError),
    #[error(transparent)]
    Kernel(#[from] ForceFieldKernelError),
    #[error(transparent)]
    Stage(#[from] OptimizationStageError),
    #[error(transparent)]
    Serial(#[from] SerialConformerOptimizationError),
    #[cfg(not(target_family = "wasm"))]
    #[error(transparent)]
    Threads(#[from] UffThreadCountError),
}
/// Keeps concrete source-stage causes without exporting private machinery.
#[derive(Debug)]
pub struct MmffOptimizationError(Failure);
impl std::fmt::Display for MmffOptimizationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.0, f)
    }
}
impl std::error::Error for MmffOptimizationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.0)
    }
}
fn construct_selected_field<'a>(
    props: &MmffMolProperties,
    coordinates: &'a mut CoordinateBlock,
    properties: &MoleculeProperties,
    id: Option<usize>,
    threshold: f64,
    ignore: bool,
) -> Result<ForceField<'a>, Failure> {
    // RDKit❗❌:   Conformer &conf = mol.getConformer(confId);
    // RDKit❗❌:   for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit❗❌:     res->positions().push_back(&(conf.getAtomPos(i)));
    // RDKit❗❌:   }
    // Behavior: select only actual stored 3D IDs, borrow rows directly. Split
    // metadata feeds the one existing sanitized fragment-copy implementation.
    // Complexity: selected ID scan O(C); one O(A) handle buffer. Selected
    // property-map clone is additional metadata work required by the detached
    // conformer's current borrow API; no topology/coordinate block is copied.
    let selected = match id {
        None if !coordinates.conformers_3d.is_empty() => 0,
        Some(requested) => coordinates
            .conformers_3d
            .iter()
            .position(|r| r.id() == requested)
            .ok_or(MmffBuilderError::Missing3dConformer {
                conf_id: requested as isize,
            })?,
        None => return Err(MmffBuilderError::Missing3dConformer { conf_id: -1 }.into()),
    };
    let (before, rest) = coordinates.conformers_3d.split_at_mut(selected);
    let (conf, after) = rest.split_first_mut().expect("selected stored conformer");
    let id = conf.id();
    let is_3d = conf.is_3d();
    let metadata = conf.props().clone();
    let positions = conf
        .coordinates_mut()
        .iter_mut()
        .map(|r| r.as_mut_slice())
        .collect();
    let context = MmffConformerContext {
        two_d: &coordinates.conformers_2d,
        before,
        selected_id: id,
        selected_is_3d: is_3d,
        selected_props: &metadata,
        after,
        source_dimension: coordinates.source_coordinate_dim,
    };
    Ok(construct_force_field_with_props(
        &props.topology,
        props,
        positions,
        &context,
        properties,
        threshold,
        ignore,
    )?)
}
/// Reproduce the Python-facing source single-conformer status wrapper.
pub fn optimize_mmff_single(
    topology: &mut TopologyBlock,
    coordinates: &mut CoordinateBlock,
    properties: &MoleculeProperties,
    params: &MmffOptimizationParams,
    supplied_rings: Option<&cosmolkit_core::RingInfo>,
) -> Result<MmffSingleOutcome, MmffOptimizationError> {
    // BEGIN RDKIT CPP FUNCTION Wrap/rdForceFields.cpp:114-129
    // RDKit❗❌: int MMFFOptimizeMolecule(ROMol &mol, std::string mmffVariant = "MMFF94",
    // RDKit❗❌:                          int maxIters = 200, double nonBondedThresh = 100.0,
    // RDKit❗❌:                          int confId = -1,
    // RDKit❗❌:                          bool ignoreInterfragInteractions = true) {
    // RDKit❗❌:   int res = -1;
    // RDKit❗❌:   MMFF::MMFFMolProperties mmffMolProperties(mol, mmffVariant);
    // RDKit❗❌:   if (mmffMolProperties.isValid()) {
    // RDKit❗❌:     NOGIL gil;
    // RDKit❗❌:     std::unique_ptr<ForceFields::ForceField> ff(
    // RDKit❗❌:         MMFF::constructForceField(mol, &mmffMolProperties, nonBondedThresh,
    // RDKit❗❌:                                   confId, ignoreInterfragInteractions));
    // RDKit❗❌:     ff->initialize();
    // RDKit❗❌:     res = ff->minimize(maxIters);
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // Behavior: source invalid typing returns -1. All actual preparation and
    // build/minimize failures stay typed. Kernel writes the selected borrowed
    // conformer rows; MMFF-prepared topology is moved back even on error.
    // Complexity: reuse one preparation, builder, kernel; no state rollback clone.
    let props = MmffMolProperties::new_prepared(
        topology,
        properties.prop("_MMFFSanitized").is_some(),
        &params.mmff_variant,
        0,
        supplied_rings,
    )
    .map_err(|e| MmffOptimizationError(Failure::Properties(e)))?;
    let result = (|| -> Result<i32, Failure> {
        if !props.is_valid() {
            return Ok(-1);
        }
        let mut field = construct_selected_field(
            &props,
            coordinates,
            properties,
            params.conformer_id,
            params.non_bonded_threshold,
            params.ignore_interfragment_interactions,
        )?;
        field.initialize()?;
        Ok(field.minimize(params.max_iterations as u32, 1e-4, 1e-6)?)
    })();
    let aromatic_ring_count = props.aromatic_ring_count;
    *topology = props.topology;
    result
        .map(|status| MmffSingleOutcome {
            status,
            aromatic_ring_count,
            acquired_rings: props.acquired_rings,
        })
        .map_err(MmffOptimizationError)
}
/// Optimize stored conformers in source order through the shared ST/MT driver.
pub fn optimize_mmff_conformers(
    topology: &mut TopologyBlock,
    coordinates: &mut CoordinateBlock,
    properties: &MoleculeProperties,
    params: &MmffConformerOptimizationParams,
    supplied_rings: Option<&cosmolkit_core::RingInfo>,
) -> Result<MmffConformerOutcomes, MmffOptimizationError> {
    // RDKit❗❌: inline void MMFFOptimizeMoleculeConfs(ROMol &mol,
    // RDKit❗❌:                                       std::vector<std::pair<int, double>> &res,
    // RDKit❗❌:                                       int numThreads = 1, int maxIters = 1000,
    // RDKit❗❌:                                       std::string mmffVariant = "MMFF94",
    // RDKit❗❌:                                       double nonBondedThresh = 10.0,
    // RDKit❗❌:                                       bool ignoreInterfragInteractions = true) {
    // RDKit❗❌:   MMFF::MMFFMolProperties mmffMolProperties(mol, mmffVariant);
    // RDKit❗❌:   if (mmffMolProperties.isValid()) {
    // RDKit❗❌:     std::unique_ptr<ForceFields::ForceField> ff(
    // RDKit❗❌:         MMFF::constructForceField(mol, &mmffMolProperties, nonBondedThresh, -1,
    // RDKit❗❌:                                   ignoreInterfragInteractions));
    // RDKit❗❌:     ForceFieldsHelper::OptimizeMoleculeConfs(mol, *ff, res, numThreads,
    // RDKit❗❌:                                              maxIters);
    // RDKit❗❌:   } else {
    // RDKit❗❌:     res.resize(mol.getNumConformers());
    // RDKit❗❌:     for (unsigned int i = 0; i < mol.getNumConformers(); ++i) {
    // RDKit❗❌:       res[i] = std::make_pair(static_cast<int>(-1), static_cast<double>(-1));
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // Behavior: preparation and construction happen before result resize.
    // Invalid types produce (-1,-1) per stored conformer. Existing shared ST/MT
    // owns field copies, rebind, initialize/minimize/energy and thread routing.
    // Complexity: one builder and existing source-shaped worker copies; O(C)
    // borrowed-row metadata and public outcome projection are additional buffers.
    let props = MmffMolProperties::new_prepared(
        topology,
        properties.prop("_MMFFSanitized").is_some(),
        &params.mmff_variant,
        0,
        supplied_rings,
    )
    .map_err(|e| MmffOptimizationError(Failure::Properties(e)))?;
    let result = (|| -> Result<Vec<MmffOptimizeMoleculeConfResult>, Failure> {
        if !props.is_valid() {
            return Ok(vec![
                MmffOptimizeMoleculeConfResult {
                    needs_more: -1,
                    energy: -1.
                };
                coordinates.conformers_3d.len()
            ]);
        }
        let field = construct_selected_field(
            &props,
            coordinates,
            properties,
            None,
            params.non_bonded_threshold,
            params.ignore_interfragment_interactions,
        )?
        .release_position_borrows();
        let mut rows = coordinates
            .conformers_3d
            .iter_mut()
            .map(|conf| SerialConformer {
                id: conf.id(),
                positions: conf.coordinates_mut(),
            })
            .collect::<Vec<_>>();
        let mut results = Vec::<OptimizationOutcome>::new();
        #[cfg(not(target_family = "wasm"))]
        {
            let hardware = std::thread::available_parallelism().map_or(1, |n| n.get()) as u32;
            let outcome = crate::uff::convenience::optimize_prepared_conformers_dispatch(
                field,
                &mut rows,
                &mut results,
                props.topology.atoms.len(),
                params.max_iterations,
                params.num_threads,
                hardware,
                true,
            )?;
            match outcome {
                PreparedConformerDispatchOutcome::Serial(result) => result?,
                PreparedConformerDispatchOutcome::Workers(result) => {
                    for joined in result? {
                        match joined {
                            Ok(result) => result?,
                            Err(payload) => std::panic::resume_unwind(payload),
                        }
                    }
                }
            }
        }
        #[cfg(target_family = "wasm")]
        {
            results.resize(
                rows.len(),
                OptimizationOutcome {
                    status: 0,
                    energy: 0.,
                },
            );
            crate::uff::convenience::optimize_serial_conformers(
                field,
                &mut rows,
                &mut results,
                props.topology.atoms.len(),
                params.max_iterations as u32,
            )?;
        }
        Ok(results
            .into_iter()
            .map(|r| MmffOptimizeMoleculeConfResult {
                needs_more: r.status,
                energy: r.energy,
            })
            .collect())
    })();
    let aromatic_ring_count = props.aromatic_ring_count;
    *topology = props.topology;
    result
        .map(|conformer_results| MmffConformerOutcomes {
            conformer_results,
            aromatic_ring_count,
            acquired_rings: props.acquired_rings,
        })
        .map_err(MmffOptimizationError)
}

#[cfg(test)]
mod original_optimization_tests {
    use super::*;
    use cosmolkit_model::{
        AtomId, AtomSpec, BondId, BondOrder, BondSpec, Conformer3D, CoordinateBlock, Element,
    };
    // Test-only detached fixture storage. No live TestInput or domain owner is
    // reimplemented. Every old literal atom/bond/conformer row is retained.
    #[derive(Clone, Debug, Default)]
    struct TestInput {
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
    }
    impl std::ops::Deref for TestInput {
        type Target = TopologyBlock;
        fn deref(&self) -> &TopologyBlock {
            &self.topology
        }
    }
    impl TestInput {
        fn new() -> Self {
            Self::default()
        }
        fn conformers_3d(&self) -> &[Conformer3D] {
            &self.coordinates.conformers_3d
        }
        fn num_atoms(&self) -> usize {
            self.topology.atoms.len()
        }
        fn num_bonds(&self) -> usize {
            self.topology.bonds.len()
        }
    }
    struct TestFixtureBuilder(TestInput);
    impl TestFixtureBuilder {
        fn new() -> Self {
            Self(TestInput::default())
        }
        fn add_atom(&mut self, spec: AtomSpec) -> AtomId {
            {
                let id = AtomId::new(self.0.topology.atoms.len());
                self.0
                    .topology
                    .atoms
                    .push(cosmolkit_model::Atom::from_spec(id, spec));
                self.0.topology.adjacency = cosmolkit_model::AdjacencyList::from_topology(
                    self.0.topology.atoms.len(),
                    &self.0.topology.bonds,
                );
                id
            }
        }
        fn add_bond(
            &mut self,
            spec: BondSpec,
        ) -> Result<BondId, cosmolkit_model::TopologyEditError> {
            {
                let id = BondId::new(self.0.topology.bonds.len());
                self.0
                    .topology
                    .bonds
                    .push(cosmolkit_model::Bond::from_spec(id, spec));
                self.0.topology.adjacency = cosmolkit_model::AdjacencyList::from_topology(
                    self.0.topology.atoms.len(),
                    &self.0.topology.bonds,
                );
                Ok(id)
            }
        }
        fn add_3d_conformer(&mut self, rows: Vec<[f64; 3]>) -> Result<(), String> {
            let id = self.0.coordinates.conformers_3d.len();
            self.add_conformer(Conformer3D::new(id, rows, true))
        }
        fn add_conformer(&mut self, row: Conformer3D) -> Result<(), String> {
            row.validate_for_atom_count(self.0.num_atoms())
                .map_err(|e| e.to_string())?;
            self.0.coordinates.conformers_3d.push(row);
            Ok(())
        }
        fn build(self) -> Result<TestInput, String> {
            self.0.topology.validate().map_err(|e| e.to_string())?;
            Ok(self.0)
        }
    }
    fn bonded_pair_with_3d_conformer(
        first: AtomSpec,
        second: AtomSpec,
        order: BondOrder,
        coords: Vec<[f64; 3]>,
    ) -> TestInput {
        let mut builder = TestFixtureBuilder::new();
        let a0 = builder.add_atom(first);
        let a1 = builder.add_atom(second);
        builder
            .add_bond(BondSpec::new(a0, a1, order))
            .expect("test bond");
        builder.add_3d_conformer(coords).expect("test 3d conformer");
        builder.build().expect("test bonded pair with conformer")
    }

    fn single_atom_with_3d_conformer(atom: AtomSpec, coord: [f64; 3]) -> TestInput {
        let mut builder = TestFixtureBuilder::new();
        builder.add_atom(atom);
        builder
            .add_3d_conformer(vec![coord])
            .expect("test 3d conformer");
        builder.build().expect("test single atom with conformer")
    }

    fn bonded_pair_with_named_3d_conformers(
        first: AtomSpec,
        second: AtomSpec,
        order: BondOrder,
        first_coords: Vec<[f64; 3]>,
        second_coords: Vec<[f64; 3]>,
    ) -> TestInput {
        let mut builder = TestFixtureBuilder::new();
        let a0 = builder.add_atom(first);
        let a1 = builder.add_atom(second);
        builder
            .add_bond(BondSpec::new(a0, a1, order))
            .expect("test bond");
        builder
            .add_conformer(Conformer3D::new(0, first_coords, true))
            .expect("first conformer");
        builder
            .add_conformer(Conformer3D::new(7, second_coords, true))
            .expect("second conformer");
        builder
            .build()
            .expect("test bonded pair with named conformers")
    }

    fn empty_molecule_with_named_3d_conformer(id: usize) -> TestInput {
        let mut builder = TestFixtureBuilder::new();
        builder
            .add_conformer(Conformer3D::new(id, vec![], true))
            .expect("empty conformer");
        builder.build().expect("empty molecule with conformer")
    }

    // Test-only value projection reuses the detached owner, preserving all old
    // fixture coordinates and assertions. Live facade tests independently cover
    // runtime COW and commit; this adapter grants no production authority.
    #[derive(Debug)]
    struct SingleTestResult {
        molecule: TestInput,
        needs_more: i32,
    }
    #[derive(Debug)]
    struct MultiTestResult {
        molecule: TestInput,
        conformer_results: Vec<MmffOptimizeMoleculeConfResult>,
    }
    fn mmff_optimize_molecule(
        input: &TestInput,
        variant: &str,
        iters: i32,
        threshold: f64,
        id: isize,
        ignore: bool,
    ) -> Result<SingleTestResult, MmffOptimizationError> {
        let mut molecule = input.clone();
        let params = MmffOptimizationParams {
            mmff_variant: variant.into(),
            max_iterations: iters,
            non_bonded_threshold: threshold,
            conformer_id: if id == -1 {
                None
            } else {
                Some(usize::try_from(id).expect("negative selectors tested by Python boundary"))
            },
            ignore_interfragment_interactions: ignore,
        };
        let outcome = optimize_mmff_single(
            &mut molecule.topology,
            &mut molecule.coordinates,
            &MoleculeProperties::default(),
            &params,
            None,
        )?;
        Ok(SingleTestResult {
            molecule,
            needs_more: outcome.status,
        })
    }
    fn mmff_optimize_molecule_confs(
        input: &TestInput,
        threads: i32,
        iters: i32,
        variant: &str,
        threshold: f64,
        ignore: bool,
    ) -> Result<MultiTestResult, MmffOptimizationError> {
        let mut molecule = input.clone();
        let params = MmffConformerOptimizationParams {
            num_threads: threads,
            max_iterations: iters,
            mmff_variant: variant.into(),
            non_bonded_threshold: threshold,
            ignore_interfragment_interactions: ignore,
        };
        let outcome = optimize_mmff_conformers(
            &mut molecule.topology,
            &mut molecule.coordinates,
            &MoleculeProperties::default(),
            &params,
            None,
        )?;
        Ok(MultiTestResult {
            molecule,
            conformer_results: outcome.conformer_results,
        })
    }
    #[test]
    fn mmff_public_api_mmff_optimize_molecule_returns_value_style_result_for_typed_molecule() {
        let molecule = bonded_pair_with_3d_conformer(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        );
        let original_coords = molecule.conformers_3d()[0].coordinates().to_vec();

        let result = mmff_optimize_molecule(&molecule, "MMFF94", 25, 100.0, -1, true)
            .expect("typed MMFF molecule should optimize");

        assert_eq!(result.needs_more, 0);
        assert_ne!(
            result.molecule.conformers_3d()[0].coordinates(),
            original_coords
        );
        assert_eq!(molecule.conformers_3d()[0].coordinates(), original_coords);
    }

    #[test]
    fn mmff_public_api_mmff_optimize_molecule_returns_minus_one_for_missing_atom_params() {
        let molecule = single_atom_with_3d_conformer(AtomSpec::new(Element::HE), [0.0, 0.0, 0.0]);
        let original_coords = molecule.conformers_3d()[0].coordinates().to_vec();

        let result = mmff_optimize_molecule(&molecule, "MMFF94", 25, 100.0, -1, true)
            .expect("missing MMFF atom typing should map to wrapper -1 result");

        assert_eq!(result.needs_more, -1);
        assert_eq!(
            result.molecule.conformers_3d()[0].coordinates(),
            original_coords
        );
        assert_eq!(molecule.conformers_3d()[0].coordinates(), original_coords);
    }

    #[test]
    fn mmff_public_api_mmff_optimize_molecule_uses_rdkit_variant_parser_for_invalid_variant() {
        let molecule = bonded_pair_with_3d_conformer(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        );

        let reference = mmff_optimize_molecule(&molecule, "MMFF94", 25, 100.0, -1, true)
            .expect("reference MMFF94 optimize should run");
        let result = mmff_optimize_molecule(&molecule, "MMFF94S", 25, 100.0, -1, true)
            .expect("invalid uppercase MMFF variant should fall back like RDKit parser");

        assert_eq!(result.needs_more, reference.needs_more);
        assert_eq!(
            result.molecule.conformers_3d()[0].coordinates(),
            reference.molecule.conformers_3d()[0].coordinates()
        );
    }

    #[test]
    fn mmff_public_api_mmff_optimize_molecule_reports_max_iteration_limit() {
        let molecule = bonded_pair_with_3d_conformer(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]],
        );
        let original_coords = molecule.conformers_3d()[0].coordinates().to_vec();

        let result = mmff_optimize_molecule(&molecule, "MMFF94", 0, 100.0, -1, true)
            .expect("empty-typed MMFF optimize should run");

        assert_eq!(result.needs_more, 1);
        assert_eq!(
            result.molecule.conformers_3d()[0].coordinates(),
            original_coords
        );
    }

    #[test]
    fn mmff_public_api_mmff_optimize_molecule_updates_selected_named_conformer_only() {
        let molecule = bonded_pair_with_named_3d_conformers(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]],
        );
        let first_coords = molecule.conformers_3d()[0].coordinates().to_vec();
        let selected_coords = molecule.conformers_3d()[1].coordinates().to_vec();

        let result = mmff_optimize_molecule(&molecule, "MMFF94", 25, 100.0, 7, true)
            .expect("MMFF optimize should preserve unselected conformers");

        assert_eq!(result.needs_more, 0);
        assert_eq!(
            result.molecule.conformers_3d()[0].coordinates(),
            first_coords
        );
        assert_ne!(
            result.molecule.conformers_3d()[1].coordinates(),
            selected_coords
        );
        assert_eq!(molecule.conformers_3d()[1].coordinates(), selected_coords);
    }

    #[test]
    fn shared_forcefield_conformer_driver_mmff_returns_value_style_results_for_all_conformers() {
        let molecule = bonded_pair_with_named_3d_conformers(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [2.3, 0.0, 0.0]],
        );
        let original_first = molecule.conformers_3d()[0].coordinates().to_vec();
        let original_second = molecule.conformers_3d()[1].coordinates().to_vec();

        let result = mmff_optimize_molecule_confs(&molecule, 1, 25, "MMFF94", 100.0, true)
            .expect("missing MMFF atom typing should return wrapper-style -1 results");

        assert_eq!(result.conformer_results.len(), 2);
        assert!(
            result
                .conformer_results
                .iter()
                .all(|entry| entry.needs_more == 0 && entry.energy.is_finite())
        );
        assert_ne!(
            result.molecule.conformers_3d()[0].coordinates(),
            original_first
        );
        assert_ne!(
            result.molecule.conformers_3d()[1].coordinates(),
            original_second
        );
        assert_eq!(molecule.conformers_3d()[0].coordinates(), original_first);
        assert_eq!(molecule.conformers_3d()[1].coordinates(), original_second);
    }

    #[test]
    fn shared_forcefield_conformer_driver_mmff_returns_minus_one_for_missing_atom_params() {
        let molecule = single_atom_with_3d_conformer(AtomSpec::new(Element::HE), [0.0, 0.0, 0.0]);
        let original_coords = molecule.conformers_3d()[0].coordinates().to_vec();

        let result = mmff_optimize_molecule_confs(&molecule, 1, 25, "MMFF94", 100.0, true)
            .expect("missing MMFF atom typing should map to wrapper -1 result");

        assert_eq!(result.conformer_results.len(), 1);
        assert_eq!(result.conformer_results[0].needs_more, -1);
        assert_eq!(result.conformer_results[0].energy, -1.0);
        assert_eq!(
            result.molecule.conformers_3d()[0].coordinates(),
            original_coords
        );
        assert_eq!(molecule.conformers_3d()[0].coordinates(), original_coords);
    }

    #[test]
    fn mmff_public_api_mmff_optimize_molecule_confs_uses_rdkit_variant_parser_for_invalid_variant()
    {
        let molecule = bonded_pair_with_named_3d_conformers(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [2.2, 0.0, 0.0]],
        );

        let result = mmff_optimize_molecule_confs(&molecule, 1, 25, "MMFF94S", 100.0, true)
            .expect("invalid uppercase MMFF variant should fall back like RDKit parser");

        let reference = mmff_optimize_molecule_confs(&molecule, 1, 25, "MMFF94", 100.0, true)
            .expect("reference MMFF94 conformer optimize should run");

        assert_eq!(result.conformer_results.len(), 2);
        assert_eq!(result.conformer_results, reference.conformer_results);
        assert_eq!(
            result.molecule.conformers_3d()[0].coordinates(),
            reference.molecule.conformers_3d()[0].coordinates()
        );
        assert_eq!(
            result.molecule.conformers_3d()[1].coordinates(),
            reference.molecule.conformers_3d()[1].coordinates()
        );
    }

    #[test]
    fn shared_forcefield_conformer_driver_mmff_handles_non_positive_thread_request_like_non_threaded_rdkit_build()
     {
        let molecule = bonded_pair_with_named_3d_conformers(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [2.2, 0.0, 0.0]],
        );

        let zero_threads = mmff_optimize_molecule_confs(&molecule, 0, 5, "MMFF94", 100.0, true)
            .expect("zero-thread request should use non-threaded RDKit path");
        let negative_threads =
            mmff_optimize_molecule_confs(&molecule, -1, 5, "MMFF94", 100.0, true)
                .expect("negative-thread request should use non-threaded RDKit path");

        assert_eq!(zero_threads.conformer_results.len(), 2);
        assert_eq!(negative_threads.conformer_results.len(), 2);
        for (left, right) in zero_threads
            .conformer_results
            .iter()
            .zip(negative_threads.conformer_results.iter())
        {
            assert_eq!(left.needs_more, right.needs_more);
            assert!(left.energy.is_finite());
            assert!(right.energy.is_finite());
        }
    }

    #[test]
    fn shared_forcefield_conformer_driver_mmff_reports_max_iteration_limit() {
        let molecule = bonded_pair_with_3d_conformer(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]],
        );
        let original_coords = molecule.conformers_3d()[0].coordinates().to_vec();

        let result = mmff_optimize_molecule_confs(&molecule, 1, 0, "MMFF94", 100.0, true)
            .expect("current modeled MMFF wrapper path should return -1 before minimization");

        assert_eq!(result.conformer_results.len(), 1);
        assert_eq!(result.conformer_results[0].needs_more, 1);
        assert!(result.conformer_results[0].energy.is_finite());
        assert_eq!(
            result.molecule.conformers_3d()[0].coordinates(),
            original_coords
        );
    }

    #[test]
    fn shared_forcefield_conformer_driver_mmff_preserves_all_named_conformers() {
        let molecule = bonded_pair_with_named_3d_conformers(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]],
        );
        let first_coords = molecule.conformers_3d()[0].coordinates().to_vec();
        let second_coords = molecule.conformers_3d()[1].coordinates().to_vec();

        let result = mmff_optimize_molecule_confs(&molecule, 1, 25, "MMFF94", 100.0, true)
            .expect("MMFF conformer optimization should preserve current modeled coordinates");

        assert_eq!(result.conformer_results.len(), 2);
        assert_ne!(
            result.molecule.conformers_3d()[0].coordinates(),
            first_coords
        );
        assert_ne!(
            result.molecule.conformers_3d()[1].coordinates(),
            second_coords
        );
        assert_eq!(molecule.conformers_3d()[0].coordinates(), first_coords);
        assert_eq!(molecule.conformers_3d()[1].coordinates(), second_coords);
    }

    #[test]
    fn mmff_optimizer_worker_dispatch_matches_source_serial_rows() {
        let molecule = bonded_pair_with_named_3d_conformers(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0., 0., 0.], [2., 0., 0.]],
            vec![[0., 0., 0.], [2.3, 0., 0.]],
        );
        let serial = mmff_optimize_molecule_confs(&molecule, 1, 25, "MMFF94", 100., true).unwrap();
        let threaded =
            mmff_optimize_molecule_confs(&molecule, 2, 25, "MMFF94", 100., true).unwrap();
        assert_eq!(serial.conformer_results, threaded.conformer_results);
        assert_eq!(serial.molecule.coordinates, threaded.molecule.coordinates);
    }
}

/// Non-mutating MMFF field evaluation over one stored 3D conformer.
#[derive(Clone, Debug, PartialEq)]
pub struct MmffEvaluationParams {
    pub mmff_variant: String,
    pub non_bonded_threshold: f64,
    pub conformer_id: Option<usize>,
    pub ignore_interfragment_interactions: bool,
}
impl Default for MmffEvaluationParams {
    fn default() -> Self {
        Self {
            mmff_variant: "MMFF94".into(),
            non_bonded_threshold: 100.,
            conformer_id: None,
            ignore_interfragment_interactions: true,
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct MmffEnergyGradient {
    pub energy: f64,
    pub gradient: Vec<f64>,
}
impl MmffEnergyGradient {
    pub fn energy(&self) -> f64 {
        self.energy
    }
    pub fn gradient(&self) -> &[f64] {
        &self.gradient
    }
}
/// Returns None for source-invalid MMFF atom parameters; actual preparation,
/// construction and numerical failures remain concrete errors.
pub fn evaluate_mmff(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
    params: &MmffEvaluationParams,
    supplied_rings: Option<&cosmolkit_core::RingInfo>,
) -> Result<Option<MmffEnergyGradient>, MmffOptimizationError> {
    // Source gradient parity helper delegates to the same MMFF constructor and
    // ForceField::calcGrad. Source energy evaluation uses ForceField::calcEnergy.
    // Both full CPP bodies remain inside the unique existing kernel owners.
    // Behavior: no live state is modified; source-invalid types are explicit None.
    // Complexity: clone selected coordinate rows only (O(A)), one existing MMFF
    // preparation and builder. Other conformers/properties remain borrowed.
    (|| -> Result<Option<MmffEnergyGradient>, Failure> {
        let props = MmffMolProperties::new_prepared(
            topology,
            properties.prop("_MMFFSanitized").is_some(),
            &params.mmff_variant,
            0,
            supplied_rings,
        )?;
        if !props.is_valid() {
            return Ok(None);
        }
        let stored = &coordinates.conformers_3d;
        let selected = match params.conformer_id {
            None if !stored.is_empty() => 0,
            Some(id) => stored.iter().position(|row| row.id() == id).ok_or(
                MmffBuilderError::Missing3dConformer {
                    conf_id: id as isize,
                },
            )?,
            None => return Err(MmffBuilderError::Missing3dConformer { conf_id: -1 }.into()),
        };
        let row = &stored[selected];
        let mut positions = row.coordinates().to_vec();
        let context = MmffConformerContext {
            two_d: &coordinates.conformers_2d,
            before: &stored[..selected],
            selected_id: row.id(),
            selected_is_3d: row.is_3d(),
            selected_props: row.props(),
            after: &stored[selected + 1..],
            source_dimension: coordinates.source_coordinate_dim,
        };
        let mut field = construct_force_field_with_props(
            &props.topology,
            &props,
            positions.iter_mut().map(|p| p.as_mut_slice()).collect(),
            &context,
            properties,
            params.non_bonded_threshold,
            params.ignore_interfragment_interactions,
        )?;
        field.initialize()?;
        let energy = field.calc_energy_current(None)?;
        let mut gradient = vec![0.; 3 * topology.atoms.len()];
        field.calc_grad_current(&mut gradient)?;
        Ok(Some(MmffEnergyGradient { energy, gradient }))
    })()
    .map_err(MmffOptimizationError)
}
