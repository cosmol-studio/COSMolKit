use crate::parameters::parameters;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;

parameters!(MmffPropertiesParams, ck::MmffPropertiesParams, { mmff_variant: String });

parameters!(SmartsWriteParams, ck::SmartsWriteParams, {
    include_atom_maps: bool,
    isomeric_smiles: bool,
    include_dative_bonds: bool,
    rooted_at_atom: Option<usize>,
});

#[pyclass(module = "fixed_projection._native")]
pub struct ForceFieldDistanceConstraintParams {
    pub(crate) inner: ck::ForceFieldDistanceConstraintParams,
}
#[pymethods]
impl ForceFieldDistanceConstraintParams {
    #[new]
    #[pyo3(signature = (*, relative, min_distance, max_distance, force_constant))]
    fn new(relative: bool, min_distance: f64, max_distance: f64, force_constant: f64) -> Self {
        Self {
            inner: ck::ForceFieldDistanceConstraintParams::new(
                relative,
                min_distance,
                max_distance,
                force_constant,
            ),
        }
    }
    #[getter]
    fn relative(&self) -> bool {
        self.inner.relative()
    }
    #[setter]
    fn set_relative(&mut self, value: bool) {
        self.inner = ck::ForceFieldDistanceConstraintParams::new(
            value,
            self.inner.min_distance(),
            self.inner.max_distance(),
            self.inner.force_constant(),
        );
    }
    #[getter]
    fn min_distance(&self) -> f64 {
        self.inner.min_distance()
    }
    #[setter]
    fn set_min_distance(&mut self, value: f64) {
        self.inner = ck::ForceFieldDistanceConstraintParams::new(
            self.inner.relative(),
            value,
            self.inner.max_distance(),
            self.inner.force_constant(),
        );
    }
    #[getter]
    fn max_distance(&self) -> f64 {
        self.inner.max_distance()
    }
    #[setter]
    fn set_max_distance(&mut self, value: f64) {
        self.inner = ck::ForceFieldDistanceConstraintParams::new(
            self.inner.relative(),
            self.inner.min_distance(),
            value,
            self.inner.force_constant(),
        );
    }
    #[getter]
    fn force_constant(&self) -> f64 {
        self.inner.force_constant()
    }
    #[setter]
    fn set_force_constant(&mut self, value: f64) {
        self.inner = ck::ForceFieldDistanceConstraintParams::new(
            self.inner.relative(),
            self.inner.min_distance(),
            self.inner.max_distance(),
            value,
        );
    }
}

#[pyclass(module = "fixed_projection._native")]
pub struct ForceFieldAngleConstraintParams {
    pub(crate) inner: ck::ForceFieldAngleConstraintParams,
}
#[pymethods]
impl ForceFieldAngleConstraintParams {
    #[new]
    #[pyo3(signature = (*, relative, min_angle_deg, max_angle_deg, force_constant))]
    fn new(relative: bool, min_angle_deg: f64, max_angle_deg: f64, force_constant: f64) -> Self {
        Self {
            inner: ck::ForceFieldAngleConstraintParams::new(
                relative,
                min_angle_deg,
                max_angle_deg,
                force_constant,
            ),
        }
    }
    #[getter]
    fn relative(&self) -> bool {
        self.inner.relative()
    }
    #[setter]
    fn set_relative(&mut self, value: bool) {
        self.inner = ck::ForceFieldAngleConstraintParams::new(
            value,
            self.inner.min_angle_deg(),
            self.inner.max_angle_deg(),
            self.inner.force_constant(),
        );
    }
    #[getter]
    fn min_angle_deg(&self) -> f64 {
        self.inner.min_angle_deg()
    }
    #[setter]
    fn set_min_angle_deg(&mut self, value: f64) {
        self.inner = ck::ForceFieldAngleConstraintParams::new(
            self.inner.relative(),
            value,
            self.inner.max_angle_deg(),
            self.inner.force_constant(),
        );
    }
    #[getter]
    fn max_angle_deg(&self) -> f64 {
        self.inner.max_angle_deg()
    }
    #[setter]
    fn set_max_angle_deg(&mut self, value: f64) {
        self.inner = ck::ForceFieldAngleConstraintParams::new(
            self.inner.relative(),
            self.inner.min_angle_deg(),
            value,
            self.inner.force_constant(),
        );
    }
    #[getter]
    fn force_constant(&self) -> f64 {
        self.inner.force_constant()
    }
    #[setter]
    fn set_force_constant(&mut self, value: f64) {
        self.inner = ck::ForceFieldAngleConstraintParams::new(
            self.inner.relative(),
            self.inner.min_angle_deg(),
            self.inner.max_angle_deg(),
            value,
        );
    }
}

#[pyclass(module = "fixed_projection._native")]
pub struct ForceFieldTorsionConstraintParams {
    pub(crate) inner: ck::ForceFieldTorsionConstraintParams,
}
#[pymethods]
impl ForceFieldTorsionConstraintParams {
    #[new]
    #[pyo3(signature = (*, relative, min_dihedral_deg, max_dihedral_deg, force_constant))]
    fn new(
        relative: bool,
        min_dihedral_deg: f64,
        max_dihedral_deg: f64,
        force_constant: f64,
    ) -> Self {
        Self {
            inner: ck::ForceFieldTorsionConstraintParams::new(
                relative,
                min_dihedral_deg,
                max_dihedral_deg,
                force_constant,
            ),
        }
    }
    #[getter]
    fn relative(&self) -> bool {
        self.inner.relative()
    }
    #[setter]
    fn set_relative(&mut self, value: bool) {
        self.inner = ck::ForceFieldTorsionConstraintParams::new(
            value,
            self.inner.min_dihedral_deg(),
            self.inner.max_dihedral_deg(),
            self.inner.force_constant(),
        );
    }
    #[getter]
    fn min_dihedral_deg(&self) -> f64 {
        self.inner.min_dihedral_deg()
    }
    #[setter]
    fn set_min_dihedral_deg(&mut self, value: f64) {
        self.inner = ck::ForceFieldTorsionConstraintParams::new(
            self.inner.relative(),
            value,
            self.inner.max_dihedral_deg(),
            self.inner.force_constant(),
        );
    }
    #[getter]
    fn max_dihedral_deg(&self) -> f64 {
        self.inner.max_dihedral_deg()
    }
    #[setter]
    fn set_max_dihedral_deg(&mut self, value: f64) {
        self.inner = ck::ForceFieldTorsionConstraintParams::new(
            self.inner.relative(),
            self.inner.min_dihedral_deg(),
            value,
            self.inner.force_constant(),
        );
    }
    #[getter]
    fn force_constant(&self) -> f64 {
        self.inner.force_constant()
    }
    #[setter]
    fn set_force_constant(&mut self, value: f64) {
        self.inner = ck::ForceFieldTorsionConstraintParams::new(
            self.inner.relative(),
            self.inner.min_dihedral_deg(),
            self.inner.max_dihedral_deg(),
            value,
        );
    }
}

parameters!(MaccsFingerprintParams, ck::MaccsFingerprintParams, { n_bits: usize });

parameters!(McsAtomCompareParameters, ck::McsAtomCompareParameters, {
    match_valences: bool,
    match_chiral_tag: bool,
    match_formal_charge: bool,
    ring_matches_ring_only: bool,
    complete_rings_only: bool,
    match_isotope: bool,
    max_distance: f64,
});
parameters!(McsBondCompareParameters, ck::McsBondCompareParameters, {
    ring_matches_ring_only: bool,
    complete_rings_only: bool,
    match_fused_rings: bool,
    match_fused_rings_strict: bool,
    match_stereo: bool,
});
parameters!(MmffEvaluationParams, ck::MmffEvaluationParams, {
    mmff_variant: String,
    non_bonded_threshold: f64,
    conformer_id: Option<usize>,
    ignore_interfragment_interactions: bool,
});
parameters!(MmffOptimizationParams, ck::MmffOptimizationParams, {
    mmff_variant: String,
    max_iterations: i32,
    non_bonded_threshold: f64,
    conformer_id: Option<usize>,
    ignore_interfragment_interactions: bool,
});
parameters!(MmffConformerOptimizationParams, ck::MmffConformerOptimizationParams, {
    num_threads: i32,
    max_iterations: i32,
    mmff_variant: String,
    non_bonded_threshold: f64,
    ignore_interfragment_interactions: bool,
});
parameters!(Coordinate3DInputParams, ck::Coordinate3DInputParams, { is_3d: bool });
parameters!(Replace3DCoordinatesParams, ck::Replace3DCoordinatesParams, { conformer_id: usize });
parameters!(MorganCallParams, ck::MorganCallParams, {
    from_atoms: Option<Vec<u32>>,
    ignore_atoms: Option<Vec<u32>>,
    custom_atom_invariants: Option<Vec<u32>>,
    custom_bond_invariants: Option<Vec<u32>>,
    conformer_id: i32,
}, methods = {
    #[staticmethod]
    #[pyo3(name = "new", signature = (from_atoms=None, ignore_atoms=None, custom_atom_invariants=None, custom_bond_invariants=None, conformer_id=ck::MorganCallParams::default().conformer_id))]
    fn registered_new(
        from_atoms: Option<Vec<u32>>,
        ignore_atoms: Option<Vec<u32>>,
        custom_atom_invariants: Option<Vec<u32>>,
        custom_bond_invariants: Option<Vec<u32>>,
        conformer_id: i32,
    ) -> Self {
        Self { inner: ck::MorganCallParams::new(from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants, conformer_id) }
    }
});
parameters!(TopologicalTorsionCallParams, ck::TopologicalTorsionCallParams, {
    from_atoms: Option<Vec<u32>>,
    ignore_atoms: Option<Vec<u32>>,
    custom_atom_invariants: Option<Vec<u32>>,
    custom_bond_invariants: Option<Vec<u32>>,
    conformer_id: i32,
    use_legacy_stereo_perception: bool,
});

parameters!(RemoveHsParams, ck::RemoveHsParams, {
    remove_degree_zero: bool,
    remove_higher_degrees: bool,
    remove_only_h_neighbors: bool,
    remove_isotopes: bool,
    remove_and_track_isotopes: bool,
    remove_dummy_neighbors: bool,
    remove_defining_bond_stereo: bool,
    remove_with_wedged_bond: bool,
    remove_with_query: bool,
    remove_mapped: bool,
    remove_in_sgroups: bool,
    show_warnings: bool,
    remove_nonimplicit: bool,
    update_explicit_count: bool,
    remove_hydrides: bool,
    remove_nontetrahedral_neighbors: bool,
    sanitize: bool,
});
parameters!(KekulizeParams, ck::KekulizeParams, {
    mark_atoms_bonds: bool,
    canonical: bool,
    max_backtracks: u32,
});
parameters!(RingSearchParams, ck::RingSearchParams, {
    include_dative_bonds: bool,
    include_hydrogen_bonds: bool,
});
parameters!(StructureTagParams, ck::StructureTagParams, {
    conformer_id: i32,
    replace_existing_tags: bool,
});
parameters!(PotentialStereoParams, ck::PotentialStereoParams, {
    clean: bool,
    flag_possible: bool,
    allow_nontetrahedral: bool,
});
parameters!(DistanceMatrixParams, ck::DistanceMatrixParams, {
    use_bond_order: bool,
    use_atom_weights: bool,
});
parameters!(DistanceMatrix3dParams, ck::DistanceMatrix3dParams, {
    conformer_id: Option<usize>,
    use_atom_weights: bool,
});
parameters!(AtomPositionParams, ck::AtomPositionParams, {
    conformer_id: Option<usize>,
});
parameters!(RandomSmilesWriteParams, ck::RandomSmilesWriteParams, {
    isomeric_smiles: bool,
    kekule: bool,
    all_bonds_explicit: bool,
    all_hydrogens_explicit: bool,
});
parameters!(PatternFingerprintParams, ck::PatternFingerprintParams, {
    n_bits: usize,
    tautomeric: bool,
});
parameters!(TopologicalFingerprintOutputRequest, ck::TopologicalFingerprintOutputRequest, {
    atom_bits: bool,
    bit_info: bool,
});
parameters!(UffEvaluationParams, ck::UffEvaluationParams, {
    vdw_threshold: f64,
    conformer_id: Option<usize>,
    ignore_interfragment_interactions: bool,
});
parameters!(UffOptimizationParams, ck::UffOptimizationParams, {
    max_iterations: i32,
    vdw_threshold: f64,
    conformer_id: Option<usize>,
    ignore_interfragment_interactions: bool,
});
parameters!(UffConformerOptimizationParams, ck::UffConformerOptimizationParams, {
    num_threads: i32,
    max_iterations: i32,
    vdw_threshold: f64,
    ignore_interfragment_interactions: bool,
});
parameters!(XyzWriteParams, ck::XyzWriteParams, {
    conformer_id: Option<usize>,
    precision: u32,
});
parameters!(InchiReadParams, ck::InchiReadParams, {
    sanitize: bool,
    remove_hs: bool,
});
parameters!(InchiWriteParams, ck::InchiWriteParams, {
    options: String,
});
parameters!(ReactionParseParams, ck::ReactionParseParams, {
    use_smiles: bool,
    sanitize: bool,
    allow_cxsmiles: bool,
    strict_cxsmiles: bool,
}, string_maps = [replacements]);
parameters!(ReactionValidationParams, ck::ReactionValidationParams, {
    silent: bool,
});
parameters!(ReactionApplyParams, ck::ReactionApplyParams, {
    remove_unmatched_atoms: bool,
});
parameters!(ReactionTemplateRemovalParams, ck::ReactionTemplateRemovalParams, {
    threshold_unmapped_atoms: f64,
    move_to_agent_templates: bool,
});
parameters!(BioPdbReadParams, ck::BioPdbReadParams, {
    max_line_length: i32,
    check_non_ascii: bool,
    ignore_ter: bool,
    split_chain_on_ter: bool,
    skip_remarks: bool,
});
parameters!(BioPdbWriteParams, ck::BioPdbWriteParams, {
    ter_records: bool,
    numbered_ter: bool,
    ter_ignores_type: bool,
    preserve_serial: bool,
    end_record: bool,
});
parameters!(BioMoleculeParams, ck::BioMoleculeParams, {
    sanitize: bool,
    remove_hs: bool,
    flavor: u32,
    proximity_bonding: bool,
});

/// Required inputs stay required; Rust's canonical constructor owns the value.
#[pyclass(module = "fixed_projection._native")]
pub struct ForceFieldPositionConstraintParams {
    pub(crate) inner: ck::ForceFieldPositionConstraintParams,
}

#[pymethods]
impl ForceFieldPositionConstraintParams {
    #[new]
    #[pyo3(signature = (*, max_displacement, force_constant))]
    fn new(max_displacement: f64, force_constant: f64) -> Self {
        Self {
            inner: ck::ForceFieldPositionConstraintParams::new(max_displacement, force_constant),
        }
    }
    #[getter]
    fn max_displacement(&self) -> f64 {
        self.inner.max_displacement()
    }
    #[setter]
    fn set_max_displacement(&mut self, value: f64) {
        self.inner =
            ck::ForceFieldPositionConstraintParams::new(value, self.inner.force_constant());
    }
    #[getter]
    fn force_constant(&self) -> f64 {
        self.inner.force_constant()
    }
    #[setter]
    fn set_force_constant(&mut self, value: f64) {
        self.inner =
            ck::ForceFieldPositionConstraintParams::new(self.inner.max_displacement(), value);
    }
}
