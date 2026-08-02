//! Source-backed RDKit force-field core primitives.

pub(crate) mod core;
pub(crate) mod crystalff;
pub mod mmff;
pub mod uff;

pub use core::{
    AngleConstraintContrib, AngleConstraintContribs, AngleConstraintContribsParams, DihedralOutput,
    DistanceConstraintContrib, DistanceConstraintContribs, DistanceConstraintContribsParams,
    ForceField, ForceFieldContrib, ForceFieldSnapshot, ForceFieldVec3, PositionConstraintContrib,
    TorsionConstraintContrib, compute_dihedral_from_flat, compute_dihedral_from_points,
    compute_dihedral_from_position_vec, normalize_angle_deg,
};
pub use crystalff::{
    TorsionAngleContribM6, TorsionAngleContribs, TorsionAngleContribsParams, calc_torsion_energy,
    calc_torsion_energy_m6,
};

pub(crate) fn optimize_molecule_confs_non_threaded<R>(
    molecule: &mut crate::Molecule,
    force_field: &mut ForceField,
    num_threads: i32,
    max_iters: usize,
    mut make_result: impl FnMut(i32, f64) -> R,
) -> Vec<R> {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/ForceFieldHelpers/FFConvenience.h :: ForceFieldsHelper::OptimizeMoleculeConfs
    // RDKit✔️❌: inline void OptimizeMoleculeConfs(ROMol &mol, ForceFields::ForceField &ff,
    // RDKit✔️❌:                                   std::vector<std::pair<int, double>> &res,
    // RDKit✔️❌:                                   int numThreads = 1, int maxIters = 1000) {
    // RDKit✔️❌:   res.resize(mol.getNumConformers());
    // RDKit✔️❌:   numThreads = getNumThreadsToUse(numThreads);
    // RDKit✔️❌:   if (numThreads == 1) {
    // RDKit✔️❌:     detail::OptimizeMoleculeConfsST(mol, ff, res, maxIters);
    // RDKit✔️❌:   }
    // RDKit✔️❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit✔️❌:   else {
    // RDKit✔️❌:     detail::OptimizeMoleculeConfsMT(mol, ff, res, numThreads, maxIters);
    // RDKit✔️❌:   }
    // RDKit✔️❌: #endif
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/ForceFieldHelpers/FFConvenience.h :: ForceFieldsHelper::OptimizeMoleculeConfs
    //
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/RDGeneral/RDThreads.h :: getNumThreadsToUse without RDK_BUILD_THREADSAFE_SSS
    // RDKit✔️❌: inline unsigned int getNumThreadsToUse(int target) {
    // RDKit✔️❌:   RDUNUSED_PARAM(target);
    // RDKit✔️❌:   return 1;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/RDGeneral/RDThreads.h :: getNumThreadsToUse without RDK_BUILD_THREADSAFE_SSS
    //
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/ForceFieldHelpers/FFConvenience.h :: ForceFieldsHelper::detail::OptimizeMoleculeConfsST
    // RDKit✔️❌: inline void OptimizeMoleculeConfsST(ROMol &mol, ForceFields::ForceField &ff,
    // RDKit✔️❌:                                     std::vector<std::pair<int, double>> &res,
    // RDKit✔️❌:                                     int maxIters) {
    // RDKit✔️❌:   PRECONDITION(res.size() >= mol.getNumConformers(),
    // RDKit✔️❌:                "res.size() must be >= mol.getNumConformers()");
    // RDKit✔️❌:   unsigned int i = 0;
    // RDKit✔️❌:   for (ROMol::ConformerIterator cit = mol.beginConformers();
    // RDKit✔️❌:        cit != mol.endConformers(); ++cit, ++i) {
    // RDKit✔️❌:     for (unsigned int aidx = 0; aidx < mol.getNumAtoms(); ++aidx) {
    // RDKit✔️❌:       ff.positions()[aidx] = &(*cit)->getAtomPos(aidx);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     ff.initialize();
    // RDKit✔️❌:     int needsMore = ff.minimize(maxIters);
    // RDKit✔️❌:     double e = ff.calcEnergy();
    // RDKit✔️❌:     res[i] = std::make_pair(needsMore, e);
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/GraphMol/ForceFieldHelpers/FFConvenience.h :: ForceFieldsHelper::detail::OptimizeMoleculeConfsST
    let _ = num_threads;
    let conformer_count = molecule.conformers_3d().len();
    let mut results = Vec::with_capacity(conformer_count);
    for conf_index in 0..conformer_count {
        let start_coords = molecule.conformers_3d()[conf_index].coordinates().to_vec();
        force_field.positions_mut().clear();
        force_field.positions_mut().extend(
            start_coords
                .iter()
                .map(|coord| ForceFieldVec3::new(coord[0], coord[1], coord[2])),
        );
        force_field.initialize();
        let needs_more = force_field.minimize(max_iters, 1.0e-4, 1.0e-6);
        let energy = force_field.calc_energy_current(None);
        results.push(make_result(needs_more, energy));
        let coords = molecule.coordinate_block_mut().conformers_3d[conf_index].coordinates_mut();
        for (coord, position) in coords.iter_mut().zip(force_field.positions()) {
            *coord = [position.x, position.y, position.z];
        }
    }
    results
}
