// Row-major assertions intentionally retain explicit `row * width + column`
// expressions, including row or column zero, to keep matrix orientation clear.
#![allow(clippy::erasing_op, clippy::identity_op)]

use super::*;
use crate::builder::MoleculeBuilder;
use crate::chemistry::forcefield::uff::atom_typer::add_atom_charge_flags_for_uff;
use crate::chemistry::forcefield::{
    ForceField, ForceFieldContrib, ForceFieldSnapshot, ForceFieldVec3,
};
use crate::{
    AtomSpec, BondSpec, Element, Molecule, ValenceModel, assign_valence, read_mol2_from_str,
};
use std::collections::BTreeMap;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;
use std::time::{Duration, Instant};

static EMBED_POINTS_CALLBACK_COUNT: AtomicUsize = AtomicUsize::new(0);

fn embed_points_test_callback(iter: u32) {
    EMBED_POINTS_CALLBACK_COUNT.fetch_add(iter as usize, Ordering::SeqCst);
}

#[test]
fn embed_failure_causes_match_rdkit_ordinals_and_names() {
    let expected = [
        (EmbedFailureCause::InitialCoords, 0, "INITIAL_COORDS"),
        (
            EmbedFailureCause::FirstMinimization,
            1,
            "FIRST_MINIMIZATION",
        ),
        (
            EmbedFailureCause::CheckTetrahedralCenters,
            2,
            "CHECK_TETRAHEDRAL_CENTERS",
        ),
        (
            EmbedFailureCause::CheckChiralCenters,
            3,
            "CHECK_CHIRAL_CENTERS",
        ),
        (
            EmbedFailureCause::MinimizeFourthDimension,
            4,
            "MINIMIZE_FOURTH_DIMENSION",
        ),
        (EmbedFailureCause::EtkMinimization, 5, "ETK_MINIMIZATION"),
        (
            EmbedFailureCause::FinalChiralBounds,
            6,
            "FINAL_CHIRAL_BOUNDS",
        ),
        (
            EmbedFailureCause::FinalCenterInVolume,
            7,
            "FINAL_CENTER_IN_VOLUME",
        ),
        (EmbedFailureCause::LinearDoubleBond, 8, "LINEAR_DOUBLE_BOND"),
        (
            EmbedFailureCause::BadDoubleBondStereo,
            9,
            "BAD_DOUBLE_BOND_STEREO",
        ),
        (
            EmbedFailureCause::CheckChiralCenters2,
            10,
            "CHECK_CHIRAL_CENTERS2",
        ),
        (EmbedFailureCause::ExceededTimeout, 11, "EXCEEDED_TIMEOUT"),
        (EmbedFailureCause::EndOfEnum, 12, "END_OF_ENUM"),
    ];

    assert_eq!(EmbedFailureCause::ALL.len(), expected.len());
    for (idx, (cause, ordinal, name)) in expected.iter().copied().enumerate() {
        assert_eq!(EmbedFailureCause::ALL[idx], cause);
        assert_eq!(cause.rdkit_ordinal(), ordinal);
        assert_eq!(EmbedFailureCause::from_rdkit_ordinal(ordinal), Some(cause));
        assert_eq!(cause.rdkit_name(), name);
    }
    assert_eq!(EmbedFailureCause::from_rdkit_ordinal(13), None);
}

#[test]
fn embed_parameters_defaults_match_rdkit_constructor_defaults() {
    let params = EmbedParameters::default();

    assert_eq!(params.max_iterations, 0);
    assert_eq!(params.num_threads, 1);
    assert_eq!(params.random_seed, -1);
    assert!(params.clear_confs);
    assert!(!params.use_random_coords);
    assert_eq!(params.box_size_mult, 2.0);
    assert!(params.rand_neg_eig);
    assert_eq!(params.num_zero_fail, 1);
    assert!(params.coord_map.is_none());
    assert_eq!(params.optimizer_force_tol, 1e-3);
    assert!(!params.ignore_smoothing_failures);
    assert!(params.enforce_chirality);
    assert!(!params.use_exp_torsion_angle_prefs);
    assert!(!params.use_basic_knowledge);
    assert!(!params.verbose);
    assert_eq!(params.basin_thresh, 5.0);
    assert_eq!(params.prune_rms_thresh, -1.0);
    assert!(params.only_heavy_atoms_for_rms);
    assert_eq!(params.et_version, 2);
    assert!(params.bounds_mat.is_none());
    assert!(params.embed_fragments_separately);
    assert!(!params.use_small_ring_torsions);
    assert!(!params.use_macrocycle_torsions);
    assert!(!params.use_macrocycle14config);
    assert_eq!(params.timeout, 0);
    assert!(params.cpci.is_none());
    assert!(params.callback.is_none());
    assert!(params.force_trans_amides);
    assert!(params.use_symmetry_for_pruning);
    assert_eq!(params.bounds_mat_force_scaling, 1.0);
    assert!(!params.track_failures);
    assert!(params.failures.is_empty());
    assert!(!params.enable_sequential_random_seeds);
    assert!(params.symmetrize_conjugated_terminal_groups_for_pruning);
}

#[test]
fn embed_parameters_new_matches_default_constructor() {
    let from_new = EmbedParameters::new();
    let from_default = EmbedParameters::default();

    assert_eq!(from_new.max_iterations, from_default.max_iterations);
    assert_eq!(from_new.num_threads, from_default.num_threads);
    assert_eq!(from_new.random_seed, from_default.random_seed);
    assert_eq!(from_new.clear_confs, from_default.clear_confs);
    assert_eq!(from_new.use_random_coords, from_default.use_random_coords);
    assert_eq!(from_new.box_size_mult, from_default.box_size_mult);
    assert_eq!(from_new.rand_neg_eig, from_default.rand_neg_eig);
    assert_eq!(from_new.num_zero_fail, from_default.num_zero_fail);
    assert_eq!(
        from_new.optimizer_force_tol,
        from_default.optimizer_force_tol
    );
    assert_eq!(
        from_new.ignore_smoothing_failures,
        from_default.ignore_smoothing_failures
    );
    assert_eq!(from_new.enforce_chirality, from_default.enforce_chirality);
    assert_eq!(
        from_new.use_exp_torsion_angle_prefs,
        from_default.use_exp_torsion_angle_prefs
    );
    assert_eq!(
        from_new.use_basic_knowledge,
        from_default.use_basic_knowledge
    );
    assert_eq!(from_new.verbose, from_default.verbose);
    assert_eq!(from_new.basin_thresh, from_default.basin_thresh);
    assert_eq!(from_new.prune_rms_thresh, from_default.prune_rms_thresh);
    assert_eq!(
        from_new.only_heavy_atoms_for_rms,
        from_default.only_heavy_atoms_for_rms
    );
    assert_eq!(from_new.et_version, from_default.et_version);
    assert_eq!(
        from_new.embed_fragments_separately,
        from_default.embed_fragments_separately
    );
    assert_eq!(
        from_new.use_small_ring_torsions,
        from_default.use_small_ring_torsions
    );
    assert_eq!(
        from_new.use_macrocycle_torsions,
        from_default.use_macrocycle_torsions
    );
    assert_eq!(
        from_new.use_macrocycle14config,
        from_default.use_macrocycle14config
    );
    assert_eq!(from_new.timeout, from_default.timeout);
    assert_eq!(from_new.force_trans_amides, from_default.force_trans_amides);
    assert_eq!(
        from_new.use_symmetry_for_pruning,
        from_default.use_symmetry_for_pruning
    );
    assert_eq!(
        from_new.bounds_mat_force_scaling,
        from_default.bounds_mat_force_scaling
    );
    assert_eq!(from_new.track_failures, from_default.track_failures);
    assert_eq!(
        from_new.enable_sequential_random_seeds,
        from_default.enable_sequential_random_seeds
    );
    assert_eq!(
        from_new.symmetrize_conjugated_terminal_groups_for_pruning,
        from_default.symmetrize_conjugated_terminal_groups_for_pruning
    );
}

#[allow(clippy::too_many_arguments)]
fn assert_embed_parameters_preset(
    params: &EmbedParameters,
    max_iterations: u32,
    num_threads: i32,
    random_seed: i32,
    clear_confs: bool,
    use_random_coords: bool,
    box_size_mult: f64,
    rand_neg_eig: bool,
    num_zero_fail: u32,
    optimizer_force_tol: f64,
    ignore_smoothing_failures: bool,
    enforce_chirality: bool,
    use_exp_torsion_angle_prefs: bool,
    use_basic_knowledge: bool,
    verbose: bool,
    basin_thresh: f64,
    prune_rms_thresh: f64,
    only_heavy_atoms_for_rms: bool,
    et_version: u32,
    embed_fragments_separately: bool,
    use_small_ring_torsions: bool,
    use_macrocycle_torsions: bool,
    use_macrocycle14config: bool,
    timeout: u32,
) {
    assert_eq!(params.max_iterations, max_iterations);
    assert_eq!(params.num_threads, num_threads);
    assert_eq!(params.random_seed, random_seed);
    assert_eq!(params.clear_confs, clear_confs);
    assert_eq!(params.use_random_coords, use_random_coords);
    assert_eq!(params.box_size_mult, box_size_mult);
    assert_eq!(params.rand_neg_eig, rand_neg_eig);
    assert_eq!(params.num_zero_fail, num_zero_fail);
    assert!(params.coord_map.is_none());
    assert_eq!(params.optimizer_force_tol, optimizer_force_tol);
    assert_eq!(params.ignore_smoothing_failures, ignore_smoothing_failures);
    assert_eq!(params.enforce_chirality, enforce_chirality);
    assert_eq!(
        params.use_exp_torsion_angle_prefs,
        use_exp_torsion_angle_prefs
    );
    assert_eq!(params.use_basic_knowledge, use_basic_knowledge);
    assert_eq!(params.verbose, verbose);
    assert_eq!(params.basin_thresh, basin_thresh);
    assert_eq!(params.prune_rms_thresh, prune_rms_thresh);
    assert_eq!(params.only_heavy_atoms_for_rms, only_heavy_atoms_for_rms);
    assert_eq!(params.et_version, et_version);
    assert!(params.bounds_mat.is_none());
    assert_eq!(
        params.embed_fragments_separately,
        embed_fragments_separately
    );
    assert_eq!(params.use_small_ring_torsions, use_small_ring_torsions);
    assert_eq!(params.use_macrocycle_torsions, use_macrocycle_torsions);
    assert_eq!(params.use_macrocycle14config, use_macrocycle14config);
    assert_eq!(params.timeout, timeout);
    assert!(params.cpci.is_none());
    assert!(params.callback.is_none());
    assert!(params.force_trans_amides);
    assert!(params.use_symmetry_for_pruning);
    assert_eq!(params.bounds_mat_force_scaling, 1.0);
    assert!(!params.track_failures);
    assert!(params.failures.is_empty());
    assert!(!params.enable_sequential_random_seeds);
    assert!(params.symmetrize_conjugated_terminal_groups_for_pruning);
}

fn finite_difference_gradient(
    mut energy: impl FnMut(&[f64]) -> f64,
    pos: &[f64],
    step: f64,
) -> Vec<f64> {
    let mut grad = vec![0.0; pos.len()];
    let mut plus = pos.to_vec();
    let mut minus = pos.to_vec();
    for i in 0..pos.len() {
        plus[i] += step;
        minus[i] -= step;
        let ep = energy(&plus);
        let em = energy(&minus);
        grad[i] = (ep - em) / (2.0 * step);
        plus[i] = pos[i];
        minus[i] = pos[i];
    }
    grad
}

fn print_debug_point_positions(stage: &str, positions: &[Vec<f64>]) {
    let coords: Vec<String> = positions
        .iter()
        .map(|point| format!("[{:.15},{:.15},{:.15}]", point[0], point[1], point[2]))
        .collect();
    println!("{stage} coords={}", coords.join(";"));
}

fn print_debug_flat_vector(stage: &str, values: &[f64]) {
    let payload: Vec<String> = values.iter().map(|value| format!("{value:.15}")).collect();
    println!("{stage} values=[{}]", payload.join(","));
}

fn print_row1_etkdg_details(details: &CrystalFFDetails) {
    println!(
        "row1_details exp_torsion_atoms={:?}",
        details.exp_torsion_atoms
    );
    println!(
        "row1_details exp_torsion_angles={:?}",
        details.exp_torsion_angles
    );
    println!("row1_details improper_atoms={:?}", details.improper_atoms);
    println!("row1_details bonds={:?}", details.bonds);
    println!("row1_details angles={:?}", details.angles);
    println!(
        "row1_details bounds_mat_force_scaling={:.15} constrained_atoms={:?}",
        details.bounds_mat_force_scaling, details.constrained_atoms
    );
}

fn print_row1_forcefield_contrib_layout(
    positions: &[ForceFieldVec3],
    mmat: &BoundsMatrix,
    details: &CrystalFFDetails,
) {
    let mut ff = ForceField::new(3);
    ff.positions_mut().extend_from_slice(positions);
    let n = positions.len();
    let mut atom_pairs = vec![false; n * n];
    let mut is_improper_constrained = vec![false; n];

    let before_exp = ff.contribs().len();
    add_experimental_torsion_terms(&mut ff, details, &mut atom_pairs, n);
    println!(
        "helper=add_experimental_torsion_terms added={} total={}",
        ff.contribs().len() - before_exp,
        ff.contribs().len()
    );

    let before_improper = ff.contribs().len();
    add_improper_torsion_terms(
        &mut ff,
        10.0,
        &details.improper_atoms,
        &mut is_improper_constrained,
    );
    println!(
        "helper=add_improper_torsion_terms added={} total={}",
        ff.contribs().len() - before_improper,
        ff.contribs().len()
    );

    let before_12 = ff.contribs().len();
    add_12_terms(
        &mut ff,
        details,
        &mut atom_pairs,
        positions,
        KNOWN_DIST_FORCE_CONSTANT,
        n,
    );
    println!(
        "helper=add_12_terms added={} total={}",
        ff.contribs().len() - before_12,
        ff.contribs().len()
    );

    let before_13 = ff.contribs().len();
    add_13_terms(
        &mut ff,
        details,
        &mut atom_pairs,
        positions,
        KNOWN_DIST_FORCE_CONSTANT,
        &is_improper_constrained,
        true,
        mmat,
        n,
    );
    println!(
        "contrib_idx={}..{} helper=add_13_terms added={}",
        before_13,
        ff.contribs().len().saturating_sub(1),
        ff.contribs().len() - before_13
    );

    let before_long = ff.contribs().len();
    add_long_range_distance_constraints(
        &mut ff,
        details,
        &atom_pairs,
        positions,
        KNOWN_DIST_FORCE_CONSTANT,
        mmat,
        n,
    );
    println!(
        "contrib_idx={} helper=add_long_range_distance_constraints added={}",
        before_long,
        ff.contribs().len() - before_long
    );
}

fn print_row1_constraint_terms(
    positions: &[ForceFieldVec3],
    mmat: &BoundsMatrix,
    details: &CrystalFFDetails,
) {
    let n = positions.len();
    let mut atom_pairs = vec![false; n * n];
    let mut is_improper_constrained = vec![false; n];

    for t in 0..details.exp_torsion_atoms.len() {
        let i = details.exp_torsion_atoms[t][0] as usize;
        let l = details.exp_torsion_atoms[t][3] as usize;
        let (a, b) = if i < l { (i, l) } else { (l, i) };
        atom_pairs[a * n + b] = true;
    }
    for improper_atom in &details.improper_atoms {
        is_improper_constrained[improper_atom[1] as usize] = true;
    }
    for &(first, second) in &details.bonds {
        let i = first as usize;
        let j = second as usize;
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        atom_pairs[a * n + b] = true;
        let d = (positions[i] - positions[j]).length();
        println!(
            "term add12 pair=({i},{j}) lower={:.15} upper={:.15} force={:.15}",
            d - KNOWN_DIST_TOL,
            d + KNOWN_DIST_TOL,
            KNOWN_DIST_FORCE_CONSTANT
        );
    }
    for angle in &details.angles {
        let i = angle[0] as usize;
        let j = angle[1] as usize;
        let k = angle[2] as usize;
        let (a, b) = if i < k { (i, k) } else { (k, i) };
        atom_pairs[a * n + b] = true;
        if is_improper_constrained[j] {
            println!(
                "term add13_dist pair=({i},{k}) lower={:.15} upper={:.15} force={:.15}",
                mmat.get_lower(i, k),
                mmat.get_upper(i, k),
                KNOWN_DIST_FORCE_CONSTANT
            );
        } else {
            let d = (positions[i] - positions[k]).length();
            println!(
                "term add13_dist pair=({i},{k}) lower={:.15} upper={:.15} force={:.15}",
                d - KNOWN_DIST_TOL,
                d + KNOWN_DIST_TOL,
                KNOWN_DIST_FORCE_CONSTANT
            );
        }
    }
    for i in 1..n {
        for j in 0..i {
            if !atom_pairs[j * n + i] {
                println!(
                    "term long_range pair=({i},{j}) lower={:.15} upper={:.15} force={:.15}",
                    mmat.get_lower(i, j),
                    mmat.get_upper(i, j),
                    details.bounds_mat_force_scaling * 10.0
                );
            }
        }
    }
}

fn print_debug_snapshot(stage: &str, iter: usize, snapshot: &ForceFieldSnapshot) {
    let coords: Vec<String> = snapshot
        .positions
        .chunks_exact(3)
        .map(|point| format!("[{:.15},{:.15},{:.15}]", point[0], point[1], point[2]))
        .collect();
    println!(
        "{stage} iter={iter} energy={:.15} coords={}",
        snapshot.energy,
        coords.join(";")
    );
}

#[test]
#[ignore = "debug helper for row-1 mixed ETKDG parity investigation"]
fn debug_ethene_row1_mixed_forcefield_gradient_breakdown() {
    let mol = Molecule::from_smiles("C=C")
        .expect("parse")
        .with_hydrogens()
        .expect("add hs");
    let mut params = EmbedParameters::etkdg_v3();
    params.random_seed = 61453;
    params.num_threads = 1;
    params.timeout = 10;
    let mut bounds = BoundsMatrix::new(mol.num_atoms());
    set_topol_bounds(
        &mol,
        &mut bounds,
        true,
        false,
        params.use_macrocycle14config,
        true,
        true,
        true,
    )
    .expect("bounds");
    bounds
        .triangle_smooth(0.0)
        .then_some(())
        .expect("triangle smoothing");
    let details = rd_distgeom_get_exp_tors_helper_with_params(&mol, &params).expect("details");

    let rdkit_coords = [
        [0.5694019518400151, 0.12481515704039713, -0.3351477553709494],
        [
            -0.5612971409218513,
            -0.13445204858356205,
            0.3138676674580523,
        ],
        [
            0.8197156052251863,
            -0.34949872921938324,
            -1.2505184200228587,
        ],
        [1.2587391802948145, 0.8464403614851113, 0.08713063971203365],
        [
            -1.2853139846722026,
            -0.8443735337722904,
            -0.05990757486282442,
        ],
        [-0.8012456117873941, 0.35706879304042133, 1.2445754430830582],
    ];
    let positions: Vec<ForceFieldVec3> = rdkit_coords
        .into_iter()
        .map(|c| ForceFieldVec3::new(c[0], c[1], c[2]))
        .collect();
    let flat: Vec<f64> = positions.iter().flat_map(|p| [p.x, p.y, p.z]).collect();
    let mut ff = construct_3d_forcefield(&bounds, &positions, &details);
    ff.initialize();

    let mut full_grad = vec![0.0; flat.len()];
    ff.calc_grad(&flat, &mut full_grad);
    let full_fd = finite_difference_gradient(|p| ff.calc_energy(p), &flat, 1.0e-7);
    println!("full_field_energy={:.16}", ff.calc_energy(&flat));
    let mut full_max = 0.0_f64;
    for i in 0..flat.len() {
        full_max = full_max.max((full_grad[i] - full_fd[i]).abs());
    }
    println!("full_field_grad_max_abs_diff={full_max:.16e}");

    for (idx, contrib) in ff.contribs().iter().enumerate() {
        let mut analytic = vec![0.0; flat.len()];
        contrib.get_grad(&flat, &mut analytic);
        let fd = finite_difference_gradient(|p| contrib.get_energy(p), &flat, 1.0e-7);
        let energy = contrib.get_energy(&flat);
        let mut max_abs = 0.0_f64;
        let mut argmax = 0usize;
        for i in 0..flat.len() {
            let diff = (analytic[i] - fd[i]).abs();
            if diff > max_abs {
                max_abs = diff;
                argmax = i;
            }
        }
        println!(
            "contrib_idx={idx} energy={energy:.16} grad_max_abs_diff={max_abs:.16e} argmax={argmax}"
        );
        if max_abs > 1.0e-5 {
            println!("analytic={analytic:?}");
            println!("fd={fd:?}");
        }
    }
}

#[test]
#[ignore = "debug helper for row-1 mixed ETKDG parity investigation"]
fn debug_ethene_row1_mixed_forcefield_minimization_progress() {
    let mol = Molecule::from_smiles("C=C")
        .expect("parse")
        .with_hydrogens()
        .expect("add hs");
    let mut params = EmbedParameters::etdg();
    params.random_seed = 61453;
    params.num_threads = 1;
    params.timeout = 10;
    let (embedded, _) = embed_molecule(&mol, &mut params).expect("etdg embed");
    let start: Vec<ForceFieldVec3> = embedded.conformers_3d()[0]
        .coordinates()
        .iter()
        .map(|c| ForceFieldVec3::new(c[0], c[1], c[2]))
        .collect();

    let mut etkdg = EmbedParameters::etkdg_v3();
    etkdg.random_seed = 61453;
    etkdg.num_threads = 1;
    etkdg.timeout = 10;
    let details = rd_distgeom_get_exp_tors_helper_with_params(&mol, &etkdg).expect("details");
    let mut bounds = BoundsMatrix::new(mol.num_atoms());
    set_topol_bounds(
        &mol,
        &mut bounds,
        true,
        false,
        etkdg.use_macrocycle14config,
        true,
        true,
        true,
    )
    .expect("bounds");
    bounds
        .triangle_smooth(0.0)
        .then_some(())
        .expect("triangle smoothing");

    let mut ff = construct_3d_forcefield(&bounds, &start, &details);
    ff.initialize();
    println!("initial_energy={:.16}", ff.calc_energy_current(None));
    for pass in 0..5 {
        let res = ff.minimize(300, etkdg.optimizer_force_tol, 1.0e-6);
        let coords: Vec<_> = ff.positions().iter().map(|p| [p.x, p.y, p.z]).collect();
        println!(
            "pass={pass} res={res} energy={:.16}",
            ff.calc_energy_current(None)
        );
        println!("coords={coords:?}");
    }
}

#[test]
#[ignore = "debug helper for row-1 ETKDG stage-by-stage parity investigation"]
fn debug_ethene_row1_etkdg_stage_trace() {
    let mol = Molecule::from_smiles("C=C")
        .expect("parse")
        .with_hydrogens()
        .expect("add hs");
    let mut params = EmbedParameters::etkdg_v3();
    params.random_seed = 61453;
    params.num_threads = 1;
    params.timeout = 10;

    let mut etkdg_details = CrystalFFDetails::default();
    embedder_init_etkdg(&mol, &params, &mut etkdg_details).expect("init etkdg");
    print_row1_etkdg_details(&etkdg_details);

    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let coord_map = params.coord_map.as_ref();
    let setup_ok = embedder_setup_initial_bounds_matrix(
        &mol,
        &mut mmat,
        coord_map,
        &params,
        &mut etkdg_details,
    )
    .expect("bounds setup");
    assert!(setup_ok);
    print_row1_etkdg_details(&etkdg_details);

    let mut chiral_centers = Vec::new();
    let mut tetrahedral_carbons = Vec::new();
    embedder_find_chiral_sets(
        &mol,
        &mut chiral_centers,
        &mut tetrahedral_carbons,
        coord_map,
    );

    let mut double_bond_ends = Vec::new();
    let mut stereo_double_bonds = Vec::new();
    embedder_find_double_bonds(
        &mol,
        &mut double_bond_ends,
        &mut stereo_double_bonds,
        coord_map,
    );

    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: Some(&etkdg_details),
        double_bond_ends: Some(&double_bond_ends),
        stereo_double_bonds: &stereo_double_bonds,
    };

    let dim = if params.use_random_coords || !chiral_centers.is_empty() {
        4
    } else {
        3
    };
    let mut positions = vec![vec![0.0; dim]; mol.num_atoms()];
    let mut dist_mat = SymmMatrix::with_value(mol.num_atoms(), 0.0);
    let mut rng = RdkitDistgeomMinStdRand::new_from_embed_points_seed(params.random_seed);

    let got_initial =
        embedder_generate_initial_coords(&mut positions, &eargs, &params, &mut dist_mat, &mut rng)
            .expect("generate initial coords");
    assert!(got_initial);
    print_debug_point_positions("after_generate_initial_coords", &positions);

    let got_first = embedder_first_minimization(&mut positions, &eargs, &params);
    assert!(got_first);
    let first_positions = point_vectors_to_forcefield_vec3(&positions);
    let mut first_ff = construct_distgeom_forcefield(
        &mmat,
        &positions,
        &chiral_centers,
        1.0,
        0.1,
        None,
        params.basin_thresh,
        None,
    );
    first_ff.initialize();
    let first_energy = first_ff.calc_energy_current(None);
    println!(
        "after_first_minimization energy={:.15} per_atom={:.15}",
        first_energy,
        first_energy / first_positions.len() as f64
    );
    print_debug_point_positions("after_first_minimization", &positions);

    let before_etk = positions.clone();
    print_debug_point_positions("before_etk_minimization", &before_etk);
    let before_etk_positions = point_vectors_to_forcefield_vec3(&before_etk);
    let mut mixed_ff = construct_3d_forcefield(&mmat, &before_etk_positions, &etkdg_details);
    mixed_ff.initialize();
    print_row1_forcefield_contrib_layout(&before_etk_positions, &mmat, &etkdg_details);
    print_row1_constraint_terms(&before_etk_positions, &mmat, &etkdg_details);
    println!(
        "before_etk_minimization_energy energy={:.15}",
        mixed_ff.calc_energy_current(None)
    );
    let mut before_grad = vec![0.0; before_etk_positions.len() * 3];
    mixed_ff.calc_grad_current(&mut before_grad);
    print_debug_flat_vector("before_etk_minimization_grad", &before_grad);
    let mut snapshots = Vec::new();
    let _ = mixed_ff.minimize_with_snapshots(
        1,
        Some(&mut snapshots),
        300,
        params.optimizer_force_tol,
        1.0e-6,
    );
    for (iter, snapshot) in snapshots.iter().enumerate().take(120) {
        print_debug_snapshot("etk_bfgs_snapshot", iter + 1, snapshot);
    }

    let got_etk = embedder_minimize_with_exp_torsions(&mut positions, &eargs, &params);
    println!("after_etk_minimization planar={}", i32::from(got_etk));
    print_debug_point_positions("after_etk_minimization", &positions);
    let after_etk_positions = point_vectors_to_forcefield_vec3(&positions);
    let mut after_ff = construct_3d_forcefield(&mmat, &after_etk_positions, &etkdg_details);
    after_ff.initialize();
    println!(
        "after_etk_minimization_energy energy={:.15}",
        after_ff.calc_energy_current(None)
    );
    let mut after_grad = vec![0.0; after_etk_positions.len() * 3];
    after_ff.calc_grad_current(&mut after_grad);
    print_debug_flat_vector("after_etk_minimization_grad", &after_grad);
}

#[test]
#[ignore = "debug helper for row-1 mixed-forcefield checkpoint comparison"]
fn debug_ethene_row1_mixed_forcefield_checkpoint_energy_breakdown() {
    let mol = Molecule::from_smiles("C=C")
        .expect("parse")
        .with_hydrogens()
        .expect("add hs");
    let params = {
        let mut params = EmbedParameters::etkdg_v3();
        params.random_seed = 61453;
        params.num_threads = 1;
        params.timeout = 10;
        params
    };
    let mut details = CrystalFFDetails::default();
    embedder_init_etkdg(&mol, &params, &mut details).expect("init etkdg");
    let mut bounds = BoundsMatrix::new(mol.num_atoms());
    let setup_ok = embedder_setup_initial_bounds_matrix(
        &mol,
        &mut bounds,
        params.coord_map.as_ref(),
        &params,
        &mut details,
    )
    .expect("bounds setup");
    assert!(setup_ok);
    let mut chiral_centers = Vec::new();
    let mut tetrahedral_carbons = Vec::new();
    embedder_find_chiral_sets(
        &mol,
        &mut chiral_centers,
        &mut tetrahedral_carbons,
        params.coord_map.as_ref(),
    );
    let mut double_bond_ends = Vec::new();
    let mut stereo_double_bonds = Vec::new();
    embedder_find_double_bonds(
        &mol,
        &mut double_bond_ends,
        &mut stereo_double_bonds,
        params.coord_map.as_ref(),
    );
    let eargs = EmbedArgs {
        mmat: &bounds,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: Some(&details),
        double_bond_ends: Some(&double_bond_ends),
        stereo_double_bonds: &stereo_double_bonds,
    };
    let dim = if params.use_random_coords || !chiral_centers.is_empty() {
        4
    } else {
        3
    };
    let mut positions = vec![vec![0.0; dim]; mol.num_atoms()];
    let mut dist_mat = SymmMatrix::with_value(mol.num_atoms(), 0.0);
    let mut rng = RdkitDistgeomMinStdRand::new_from_embed_points_seed(params.random_seed);
    let got_initial =
        embedder_generate_initial_coords(&mut positions, &eargs, &params, &mut dist_mat, &mut rng)
            .expect("generate initial coords");
    assert!(got_initial);
    let got_first = embedder_first_minimization(&mut positions, &eargs, &params);
    assert!(got_first);
    let before_etk_positions = point_vectors_to_forcefield_vec3(&positions);
    let mut ff = construct_3d_forcefield(&bounds, &before_etk_positions, &details);
    ff.initialize();

    let checkpoints = [
        (
            "rdkit_iter1",
            [
                [
                    0.6393890874626489,
                    -0.4823664506293169,
                    -0.11411397341442445,
                ],
                [-0.6663053606581905, 0.7685743107744221, 0.47361592884366016],
                [1.3947193657309926, -0.5444070942789608, -0.548973067307731],
                [1.02264410677849, 0.9201430550357918, 0.5502878117551371],
                [-1.025044014902383, -1.134436471272739, 0.4286481501239048],
                [-1.3654031844115575, 0.472492650370804, -0.7894648500005472],
            ],
        ),
        (
            "rdkit_iter2",
            [
                [
                    0.6198122882493197,
                    -0.48668726967691534,
                    -0.13634582605517215,
                ],
                [-0.6620921321332481, 0.7721028221849652, 0.48190001466984145],
                [1.403910143856784, -0.5395274140418864, -0.5304450152056254],
                [1.0248558251268691, 0.9202424679650751, 0.5452851217818754],
                [-1.022622585062977, -1.1338984898673572, 0.4280878672796878],
                [-1.3638635400367474, 0.4677678834361196, -0.7884821624706077],
            ],
        ),
        (
            "rdkit_iter2_inner1_trial",
            [
                [0.4436210953293569, -0.5255746411053016, -0.3364324998219015],
                [-0.624173075408767, 0.8038594248798533, 0.5564567871054729],
                [1.486627146988907, -0.49561029190821643, -0.363692546286675],
                [1.0447612902622825, 0.9211371843286256, 0.5002609120225199],
                [-1.0008297165083231, -1.1290566572189196, 0.423045321681735],
                [-1.3500067406634559, 0.4252449810239599, -0.7796379747011519],
            ],
        ),
        (
            "rdkit_iter3",
            [
                [
                    0.5877070454404664,
                    -0.4412761330622238,
                    -0.12417718642276727,
                ],
                [-0.6550505877365503, 0.7546373534586957, 0.5143742972512187],
                [1.382920260652038, -0.5392572677333918, -0.49030220995890705],
                [1.0489528162110517, 0.9102790984626599, 0.49474385079631134],
                [-1.0183778958958274, -1.109322517431635, 0.37636111703800157],
                [-1.346151638671178, 0.42493946630589574, -0.7709998687038578],
            ],
        ),
        (
            "rdkit_iter4",
            [
                [
                    0.5945782026504608,
                    -0.38286741075090656,
                    -0.1382834305334702,
                ],
                [-0.6563798935334102, 0.7211672800794769, 0.5495529970118774],
                [
                    1.3615860327727195,
                    -0.5737077825722357,
                    -0.43852640889900663,
                ],
                [1.0767298638008433, 0.9282388706685269, 0.44714847321875857],
                [
                    -1.0273123705982818,
                    -1.1090270704804573,
                    0.33991186443480775,
                ],
                [-1.3492018350923309, 0.4161961130555967, -0.7598034952329669],
            ],
        ),
        (
            "rdkit_iter5",
            [
                [
                    0.30617829081320283,
                    -0.2824487117997652,
                    0.18999450067944895,
                ],
                [-0.5762354466460521, 0.19721581687982089, 0.5602906852395194],
                [
                    1.3268680369810031,
                    -0.5293442493884913,
                    -0.09494799647312154,
                ],
                [1.4356588033017199, 1.1814346233845391, -0.1128639886771219],
                [
                    -1.1114876167252854,
                    -0.9339666490256927,
                    -0.048594636608771914,
                ],
                [-1.3809820677245899, 0.3671091699495898, -0.4938785641599527],
            ],
        ),
        (
            "rdkit_iter56",
            [
                [0.581041357508865, -0.09149956537660466, -0.2168038920850982],
                [-0.6690575493352357, 0.08411148246009337, 0.235507082539928],
                [1.4860858479004584, -0.2315212262936618, -0.5571487197296577],
                [1.111788139149555, 0.9393306597875761, 0.5016014673613825],
                [
                    -1.0311017762331927,
                    -0.9696977496965792,
                    -0.5574542501005204,
                ],
                [-1.4787560190020383, 0.2692763991140623, 0.5942983120111444],
            ],
        ),
        (
            "rdkit_iter67",
            [
                [0.5659519333655081, 0.09456786688077949, -0.3498907494411227],
                [
                    -0.5647174092513038,
                    -0.11788255545663441,
                    0.31749519602622234,
                ],
                [0.8485670583229148, -0.3723581483215544, -1.2526070818707484],
                [1.2549339738080016, 0.8318163037371477, 0.09824732820612223],
                [
                    -1.2854319751027474,
                    -0.8396265429949134,
                    -0.0742921358030683,
                ],
                [-0.819303581163789, 0.403483076145871, 1.2610474428789906],
            ],
        ),
        (
            "rdkit_iter72",
            [
                [
                    0.5649408856122422,
                    0.10541291418263296,
                    -0.34411644938593777,
                ],
                [
                    -0.5594806507014879,
                    -0.12264172784088852,
                    0.3155697283270235,
                ],
                [
                    0.8253454438177249,
                    -0.37258269733839006,
                    -1.2553941399997504,
                ],
                [1.2542687371718235, 0.8361498335645288, 0.09531551656722793],
                [
                    -1.2869925434601595,
                    -0.8344113660702405,
                    -0.06306170803606982,
                ],
                [-0.7980818724616817, 0.3880730434930072, 1.2516870525239427],
            ],
        ),
        (
            "rdkit_final",
            [
                [0.5694019518400151, 0.12481515704039713, -0.3351477553709494],
                [
                    -0.5612971409218513,
                    -0.13445204858356205,
                    0.3138676674580523,
                ],
                [
                    0.8197156052251863,
                    -0.34949872921938324,
                    -1.2505184200228587,
                ],
                [1.2587391802948145, 0.8464403614851113, 0.08713063971203365],
                [
                    -1.2853139846722026,
                    -0.8443735337722904,
                    -0.05990757486282442,
                ],
                [-0.8012456117873941, 0.35706879304042133, 1.2445754430830582],
            ],
        ),
    ];

    for (label, coords) in checkpoints {
        let flat: Vec<f64> = coords.into_iter().flatten().collect();
        println!("{label} total_energy={:.15}", ff.calc_energy(&flat));
        let mut grad = vec![0.0; flat.len()];
        ff.calc_grad(&flat, &mut grad);
        print_debug_flat_vector(&format!("{label}_grad"), &grad);
        if label == "rdkit_iter1" {
            let bits: Vec<String> = grad
                .iter()
                .map(|value| format!("{:#018x}", value.to_bits()))
                .collect();
            println!("{label}_grad_bits values=[{}]", bits.join(","));
        }
        for (idx, contrib) in ff.contribs().iter().enumerate() {
            println!(
                "{label} contrib_idx={idx} energy={:.15}",
                contrib.get_energy(&flat)
            );
            let mut contrib_grad = vec![0.0; flat.len()];
            contrib.get_grad(&flat, &mut contrib_grad);
            print_debug_flat_vector(&format!("{label}_contrib_{idx}_grad"), &contrib_grad);
            if label == "rdkit_iter1" {
                let bits: Vec<String> = contrib_grad
                    .iter()
                    .map(|value| format!("{:#018x}", value.to_bits()))
                    .collect();
                println!(
                    "{label}_contrib_{idx}_grad_bits values=[{}]",
                    bits.join(",")
                );
            }
        }
    }
}

#[test]
#[ignore = "debug helper for row-16 ETKDG stage-by-stage parity investigation"]
fn debug_azide_row16_etkdg_stage_trace() {
    let mol = Molecule::from_smiles("[N-]=[N+]=N")
        .expect("parse")
        .with_hydrogens()
        .expect("add hs");
    let mut params = EmbedParameters::etkdg_v3();
    params.random_seed = 61453;
    params.num_threads = 1;
    params.timeout = 10;

    let mut etkdg_details = CrystalFFDetails::default();
    embedder_init_etkdg(&mol, &params, &mut etkdg_details).expect("init etkdg");

    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let coord_map = params.coord_map.as_ref();
    let setup_ok = embedder_setup_initial_bounds_matrix(
        &mol,
        &mut mmat,
        coord_map,
        &params,
        &mut etkdg_details,
    )
    .expect("bounds setup");
    assert!(setup_ok);

    let mut chiral_centers = Vec::new();
    let mut tetrahedral_carbons = Vec::new();
    embedder_find_chiral_sets(
        &mol,
        &mut chiral_centers,
        &mut tetrahedral_carbons,
        coord_map,
    );

    let mut double_bond_ends = Vec::new();
    let mut stereo_double_bonds = Vec::new();
    embedder_find_double_bonds(
        &mol,
        &mut double_bond_ends,
        &mut stereo_double_bonds,
        coord_map,
    );

    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: Some(&etkdg_details),
        double_bond_ends: Some(&double_bond_ends),
        stereo_double_bonds: &stereo_double_bonds,
    };

    println!(
        "row16_details exp_torsion_atoms={:?} exp_torsion_angles={:?} improper_atoms={:?} bonds={:?} angles={:?} stereo_double_bonds={:?}",
        etkdg_details.exp_torsion_atoms,
        etkdg_details.exp_torsion_angles,
        etkdg_details.improper_atoms,
        etkdg_details.bonds,
        etkdg_details.angles,
        stereo_double_bonds
    );

    let dim = if params.use_random_coords || !chiral_centers.is_empty() {
        4
    } else {
        3
    };
    let mut positions = vec![vec![0.0; dim]; mol.num_atoms()];
    let mut dist_mat = SymmMatrix::with_value(mol.num_atoms(), 0.0);
    let mut rng = RdkitDistgeomMinStdRand::new_from_embed_points_seed(params.random_seed);

    let got_initial =
        embedder_generate_initial_coords(&mut positions, &eargs, &params, &mut dist_mat, &mut rng)
            .expect("generate initial coords");
    assert!(got_initial);
    print_debug_point_positions("row16_after_generate_initial_coords", &positions);

    let got_first = embedder_first_minimization(&mut positions, &eargs, &params);
    assert!(got_first);
    let first_positions = point_vectors_to_forcefield_vec3(&positions);
    let mut first_ff = construct_distgeom_forcefield(
        &mmat,
        &positions,
        &chiral_centers,
        1.0,
        0.1,
        None,
        params.basin_thresh,
        None,
    );
    first_ff.initialize();
    let first_energy = first_ff.calc_energy_current(None);
    println!(
        "row16_after_first_minimization energy={:.15} per_atom={:.15}",
        first_energy,
        first_energy / first_positions.len() as f64
    );
    print_debug_point_positions("row16_after_first_minimization", &positions);

    let before_etk = positions.clone();
    print_debug_point_positions("row16_before_etk_minimization", &before_etk);
    let before_etk_positions = point_vectors_to_forcefield_vec3(&before_etk);
    let mut mixed_ff = construct_3d_forcefield(&mmat, &before_etk_positions, &etkdg_details);
    mixed_ff.initialize();
    print_row1_forcefield_contrib_layout(&before_etk_positions, &mmat, &etkdg_details);
    print_row1_constraint_terms(&before_etk_positions, &mmat, &etkdg_details);
    println!(
        "row16_before_etk_minimization_energy energy={:.15}",
        mixed_ff.calc_energy_current(None)
    );
    let mut before_grad = vec![0.0; before_etk_positions.len() * 3];
    mixed_ff.calc_grad_current(&mut before_grad);
    print_debug_flat_vector("row16_before_etk_minimization_grad", &before_grad);
    let before_flat: Vec<f64> = before_etk_positions
        .iter()
        .flat_map(|p| [p.x, p.y, p.z])
        .collect();
    for (idx, contrib) in mixed_ff.contribs().iter().enumerate() {
        let energy = contrib.get_energy(&before_flat);
        let mut contrib_grad = vec![0.0; before_flat.len()];
        contrib.get_grad(&before_flat, &mut contrib_grad);
        println!("row16_contrib idx={idx} energy={energy:.15} grad={contrib_grad:?}");
    }

    let got_etk = embedder_minimize_with_exp_torsions(&mut positions, &eargs, &params);
    println!("row16_after_etk_minimization planar={}", i32::from(got_etk));
    print_debug_point_positions("row16_after_etk_minimization", &positions);
    let after_etk_positions = point_vectors_to_forcefield_vec3(&positions);
    let mut after_ff = construct_3d_forcefield(&mmat, &after_etk_positions, &etkdg_details);
    after_ff.initialize();
    println!(
        "row16_after_etk_minimization_energy energy={:.15}",
        after_ff.calc_energy_current(None)
    );
    let mut after_grad = vec![0.0; after_etk_positions.len() * 3];
    after_ff.calc_grad_current(&mut after_grad);
    print_debug_flat_vector("row16_after_etk_minimization_grad", &after_grad);
}

#[test]
#[ignore = "debug helper for row-34 ETKDG stage-by-stage parity investigation"]
fn debug_row34_etkdg_stage_trace() {
    let mol = Molecule::from_smiles("CCCC1CCC(c2ccc(OCC)cc2)CC1")
        .expect("parse")
        .with_hydrogens()
        .expect("add hs");
    let mut params = EmbedParameters::etkdg_v3();
    params.random_seed = 61453;
    params.num_threads = 1;
    params.timeout = 10;

    let mut etkdg_details = CrystalFFDetails::default();
    embedder_init_etkdg(&mol, &params, &mut etkdg_details).expect("init etkdg");

    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let coord_map = params.coord_map.as_ref();
    let setup_ok = embedder_setup_initial_bounds_matrix(
        &mol,
        &mut mmat,
        coord_map,
        &params,
        &mut etkdg_details,
    )
    .expect("bounds setup");
    assert!(setup_ok);

    let mut chiral_centers = Vec::new();
    let mut tetrahedral_carbons = Vec::new();
    embedder_find_chiral_sets(
        &mol,
        &mut chiral_centers,
        &mut tetrahedral_carbons,
        coord_map,
    );

    let mut double_bond_ends = Vec::new();
    let mut stereo_double_bonds = Vec::new();
    embedder_find_double_bonds(
        &mol,
        &mut double_bond_ends,
        &mut stereo_double_bonds,
        coord_map,
    );

    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: Some(&etkdg_details),
        double_bond_ends: Some(&double_bond_ends),
        stereo_double_bonds: &stereo_double_bonds,
    };

    println!(
        "row34_details exp_torsion_atoms={:?} exp_torsion_angles={:?} improper_atoms={:?} bonds={:?} angles={:?} stereo_double_bonds={:?}",
        etkdg_details.exp_torsion_atoms,
        etkdg_details.exp_torsion_angles,
        etkdg_details.improper_atoms,
        etkdg_details.bonds,
        etkdg_details.angles,
        stereo_double_bonds
    );

    let mut positions = vec![vec![0.0; 3]; mol.num_atoms()];
    let mut dist_mat = SymmMatrix::with_value(mol.num_atoms(), 0.0);
    let mut rng = RdkitDistgeomMinStdRand::new_from_embed_points_seed(params.random_seed);

    let got_initial =
        embedder_generate_initial_coords(&mut positions, &eargs, &params, &mut dist_mat, &mut rng)
            .expect("generate initial coords");
    assert!(got_initial);
    print_debug_point_positions("row34_after_generate_initial_coords", &positions);

    let got_first = embedder_first_minimization(&mut positions, &eargs, &params);
    assert!(got_first);
    let first_positions = point_vectors_to_forcefield_vec3(&positions);
    let mut first_ff = construct_distgeom_forcefield(
        &mmat,
        &positions,
        &chiral_centers,
        1.0,
        0.1,
        None,
        params.basin_thresh,
        None,
    );
    first_ff.initialize();
    let first_energy = first_ff.calc_energy_current(None);
    println!(
        "row34_after_first_minimization energy={:.15} per_atom={:.15}",
        first_energy,
        first_energy / first_positions.len() as f64
    );
    print_debug_point_positions("row34_after_first_minimization", &positions);

    let before_etk = positions.clone();
    print_debug_point_positions("row34_before_etk_minimization", &before_etk);
    let before_etk_positions = point_vectors_to_forcefield_vec3(&before_etk);
    let mut mixed_ff = construct_3d_forcefield(&mmat, &before_etk_positions, &etkdg_details);
    mixed_ff.initialize();
    print_row1_forcefield_contrib_layout(&before_etk_positions, &mmat, &etkdg_details);
    print_row1_constraint_terms(&before_etk_positions, &mmat, &etkdg_details);
    println!(
        "row34_before_etk_minimization_energy energy={:.15}",
        mixed_ff.calc_energy_current(None)
    );
    let mut before_grad = vec![0.0; before_etk_positions.len() * 3];
    mixed_ff.calc_grad_current(&mut before_grad);
    print_debug_flat_vector("row34_before_etk_minimization_grad", &before_grad);
    let before_flat: Vec<f64> = before_etk_positions
        .iter()
        .flat_map(|p| [p.x, p.y, p.z])
        .collect();
    for (idx, contrib) in mixed_ff.contribs().iter().enumerate() {
        let energy = contrib.get_energy(&before_flat);
        let mut contrib_grad = vec![0.0; before_flat.len()];
        contrib.get_grad(&before_flat, &mut contrib_grad);
        println!("row34_contrib idx={idx} energy={energy:.15} grad={contrib_grad:?}");
    }

    let got_etk = embedder_minimize_with_exp_torsions(&mut positions, &eargs, &params);
    println!("row34_after_etk_minimization planar={}", i32::from(got_etk));
    print_debug_point_positions("row34_after_etk_minimization", &positions);
    let after_etk_positions = point_vectors_to_forcefield_vec3(&positions);
    let mut after_ff = construct_3d_forcefield(&mmat, &after_etk_positions, &etkdg_details);
    after_ff.initialize();
    println!(
        "row34_after_etk_minimization_energy energy={:.15}",
        after_ff.calc_energy_current(None)
    );
    let mut after_grad = vec![0.0; after_etk_positions.len() * 3];
    after_ff.calc_grad_current(&mut after_grad);
    print_debug_flat_vector("row34_after_etk_minimization_grad", &after_grad);
}

#[test]
#[ignore = "debug helper for row-34 embedPoints failure-cause investigation"]
fn debug_row34_embed_points_failure_trace() {
    let mol = Molecule::from_smiles("CCCC1CCC(c2ccc(OCC)cc2)CC1")
        .expect("parse")
        .with_hydrogens()
        .expect("add hs");
    let mut params = EmbedParameters::etkdg_v3();
    params.random_seed = 61453;
    params.num_threads = 1;
    params.timeout = 10;
    params.track_failures = true;

    let mut etkdg_details = CrystalFFDetails::default();
    embedder_init_etkdg(&mol, &params, &mut etkdg_details).expect("init etkdg");

    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let coord_map = params.coord_map.as_ref();
    let setup_ok = embedder_setup_initial_bounds_matrix(
        &mol,
        &mut mmat,
        coord_map,
        &params,
        &mut etkdg_details,
    )
    .expect("bounds setup");
    assert!(setup_ok);

    let mut chiral_centers = Vec::new();
    let mut tetrahedral_carbons = Vec::new();
    embedder_find_chiral_sets(
        &mol,
        &mut chiral_centers,
        &mut tetrahedral_carbons,
        coord_map,
    );

    let mut double_bond_ends = Vec::new();
    let mut stereo_double_bonds = Vec::new();
    embedder_find_double_bonds(
        &mol,
        &mut double_bond_ends,
        &mut stereo_double_bonds,
        coord_map,
    );

    println!(
        "row34_embed_points_setup chiral_centers={} tetrahedral_carbons={} double_bond_ends={} stereo_double_bonds={}",
        chiral_centers.len(),
        tetrahedral_carbons.len(),
        double_bond_ends.len(),
        stereo_double_bonds.len()
    );

    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: Some(&etkdg_details),
        double_bond_ends: Some(&double_bond_ends),
        stereo_double_bonds: &stereo_double_bonds,
    };

    let dim = if params.use_random_coords || !chiral_centers.is_empty() {
        4
    } else {
        3
    };
    let mut positions = vec![vec![0.0; dim]; mol.num_atoms()];
    let got_coords = embedder_embed_points(&mut positions, eargs, &mut params, 61453, None)
        .expect("embed points");
    println!("row34_embed_points_result got_coords={got_coords}");
    println!("row34_embed_points_failures={:?}", params.failures);
    if got_coords {
        print_debug_point_positions("row34_embed_points_final", &positions);
    }
}

#[test]
#[ignore = "debug helper for row-34 timeout stage timing"]
fn debug_row34_timeout_timing_trace() {
    use std::time::Instant;

    let mol = Molecule::from_smiles("CCCC1CCC(c2ccc(OCC)cc2)CC1")
        .expect("parse")
        .with_hydrogens()
        .expect("add hs");
    let mut params = EmbedParameters::etkdg_v3();
    params.random_seed = 61453;
    params.num_threads = 1;
    params.timeout = 10;
    params.track_failures = true;

    let t0 = Instant::now();
    let mut etkdg_details = CrystalFFDetails::default();
    embedder_init_etkdg(&mol, &params, &mut etkdg_details).expect("init etkdg");
    let t1 = Instant::now();

    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let coord_map = params.coord_map.as_ref();
    let setup_ok = embedder_setup_initial_bounds_matrix(
        &mol,
        &mut mmat,
        coord_map,
        &params,
        &mut etkdg_details,
    )
    .expect("bounds setup");
    let t2 = Instant::now();
    assert!(setup_ok);

    let mut chiral_centers = Vec::new();
    let mut tetrahedral_carbons = Vec::new();
    embedder_find_chiral_sets(
        &mol,
        &mut chiral_centers,
        &mut tetrahedral_carbons,
        coord_map,
    );
    let t3 = Instant::now();

    let mut double_bond_ends = Vec::new();
    let mut stereo_double_bonds = Vec::new();
    embedder_find_double_bonds(
        &mol,
        &mut double_bond_ends,
        &mut stereo_double_bonds,
        coord_map,
    );
    let t4 = Instant::now();

    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: Some(&etkdg_details),
        double_bond_ends: Some(&double_bond_ends),
        stereo_double_bonds: &stereo_double_bonds,
    };
    let dim = if params.use_random_coords || !chiral_centers.is_empty() {
        4
    } else {
        3
    };
    let mut positions = vec![vec![0.0; dim]; mol.num_atoms()];
    let got_coords = embedder_embed_points(&mut positions, eargs, &mut params, 61453, None)
        .expect("embed points");
    let t5 = Instant::now();

    let mut full_params = EmbedParameters::etkdg_v3();
    full_params.random_seed = 61453;
    full_params.num_threads = 1;
    full_params.timeout = 10;
    full_params.track_failures = true;
    let (_embedded, status) = embed_molecule(&mol, &mut full_params).expect("full embed");
    let t6 = Instant::now();

    println!(
        "row34_timing init_etkdg={:.6} bounds={:.6} chiral={:.6} double_bonds={:.6} embed_points={:.6} full_embed={:.6}",
        (t1 - t0).as_secs_f64(),
        (t2 - t1).as_secs_f64(),
        (t3 - t2).as_secs_f64(),
        (t4 - t3).as_secs_f64(),
        (t5 - t4).as_secs_f64(),
        (t6 - t5).as_secs_f64(),
    );
    println!(
        "row34_timing embed_points_ok={got_coords} full_status={status} embed_point_failures={:?} full_failures={:?}",
        params.failures, full_params.failures
    );
}

#[test]
#[ignore = "debug helper for row-20 ETKDG stage-by-stage parity investigation"]
fn debug_row20_etkdg_stage_trace() {
    let mol = Molecule::from_smiles("N[C@@H](C)C(=O)O")
        .expect("parse")
        .with_hydrogens()
        .expect("add hs");
    let mut params = EmbedParameters::etkdg_v3();
    params.random_seed = 61453;
    params.num_threads = 1;
    params.timeout = 10;

    let mut etkdg_details = CrystalFFDetails::default();
    embedder_init_etkdg(&mol, &params, &mut etkdg_details).expect("init etkdg");

    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let coord_map = params.coord_map.as_ref();
    let setup_ok = embedder_setup_initial_bounds_matrix(
        &mol,
        &mut mmat,
        coord_map,
        &params,
        &mut etkdg_details,
    )
    .expect("bounds setup");
    assert!(setup_ok);

    let mut chiral_centers = Vec::new();
    let mut tetrahedral_carbons = Vec::new();
    embedder_find_chiral_sets(
        &mol,
        &mut chiral_centers,
        &mut tetrahedral_carbons,
        coord_map,
    );

    let mut double_bond_ends = Vec::new();
    let mut stereo_double_bonds = Vec::new();
    embedder_find_double_bonds(
        &mol,
        &mut double_bond_ends,
        &mut stereo_double_bonds,
        coord_map,
    );

    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: Some(&etkdg_details),
        double_bond_ends: Some(&double_bond_ends),
        stereo_double_bonds: &stereo_double_bonds,
    };

    println!(
        "row20_details chiral_centers={:?} tetrahedral_carbons={:?} exp_torsion_atoms={:?} exp_torsion_angles={:?} improper_atoms={:?} bonds={:?} angles={:?} stereo_double_bonds={:?}",
        chiral_centers,
        tetrahedral_carbons,
        etkdg_details.exp_torsion_atoms,
        etkdg_details.exp_torsion_angles,
        etkdg_details.improper_atoms,
        etkdg_details.bonds,
        etkdg_details.angles,
        stereo_double_bonds
    );

    let dim = if params.use_random_coords || !chiral_centers.is_empty() {
        4
    } else {
        3
    };
    let mut positions = vec![vec![0.0; dim]; mol.num_atoms()];
    let mut dist_mat = SymmMatrix::with_value(mol.num_atoms(), 0.0);
    let mut rng = RdkitDistgeomMinStdRand::new_from_embed_points_seed(params.random_seed);

    let got_initial =
        embedder_generate_initial_coords(&mut positions, &eargs, &params, &mut dist_mat, &mut rng)
            .expect("generate initial coords");
    assert!(got_initial);
    print_debug_point_positions("row20_after_generate_initial_coords", &positions);
    let mut initial_ff = construct_distgeom_forcefield(
        &mmat,
        &positions,
        &chiral_centers,
        1.0,
        0.1,
        None,
        params.basin_thresh,
        None,
    );
    initial_ff.initialize();
    println!(
        "row20_before_first_minimization_energy energy={:.15}",
        initial_ff.calc_energy_current(None)
    );
    let mut initial_grad = vec![0.0; positions.len() * dim];
    initial_ff.calc_grad_current(&mut initial_grad);
    print_debug_flat_vector("row20_before_first_minimization_grad", &initial_grad);
    let initial_flat: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();
    for (idx, contrib) in initial_ff.contribs().iter().enumerate() {
        let energy = contrib.get_energy(&initial_flat);
        let mut contrib_grad = vec![0.0; initial_flat.len()];
        contrib.get_grad(&initial_flat, &mut contrib_grad);
        println!("row20_first_contrib idx={idx} energy={energy:.15} grad={contrib_grad:?}");
    }

    let got_first = embedder_first_minimization(&mut positions, &eargs, &params);
    assert!(got_first);
    let mut first_ff = construct_distgeom_forcefield(
        &mmat,
        &positions,
        &chiral_centers,
        1.0,
        0.1,
        None,
        params.basin_thresh,
        None,
    );
    first_ff.initialize();
    let first_energy = first_ff.calc_energy_current(None);
    println!(
        "row20_after_first_minimization energy={:.15} per_atom={:.15}",
        first_energy,
        first_energy / positions.len() as f64
    );
    print_debug_point_positions("row20_after_first_minimization", &positions);

    assert!(embedder_check_tetrahedral_centers(
        &positions, &eargs, &params
    ));
    assert!(embedder_check_chiral_centers(&positions, &eargs, &params));

    let got_fourth = embedder_minimize_fourth_dimension(&mut positions, &eargs, &mut params, None);
    assert!(got_fourth);
    let fourth_positions = point_vectors_to_forcefield_vec3(&positions);
    let mut fourth_ff = construct_distgeom_forcefield(
        &mmat,
        &positions,
        &chiral_centers,
        0.2,
        1.0,
        None,
        params.basin_thresh,
        None,
    );
    fourth_ff.initialize();
    let fourth_energy = fourth_ff.calc_energy_current(None);
    println!(
        "row20_after_fourth_dimension energy={:.15} per_atom={:.15}",
        fourth_energy,
        fourth_energy / fourth_positions.len() as f64
    );
    print_debug_point_positions("row20_after_fourth_dimension", &positions);

    let before_etk = positions.clone();
    print_debug_point_positions("row20_before_etk_minimization", &before_etk);
    let before_etk_positions = point_vectors_to_forcefield_vec3(&before_etk);
    let mut mixed_ff = construct_3d_forcefield(&mmat, &before_etk_positions, &etkdg_details);
    mixed_ff.initialize();
    print_row1_forcefield_contrib_layout(&before_etk_positions, &mmat, &etkdg_details);
    print_row1_constraint_terms(&before_etk_positions, &mmat, &etkdg_details);
    println!(
        "row20_before_etk_minimization_energy energy={:.15}",
        mixed_ff.calc_energy_current(None)
    );
    let mut before_grad = vec![0.0; before_etk_positions.len() * 3];
    mixed_ff.calc_grad_current(&mut before_grad);
    print_debug_flat_vector("row20_before_etk_minimization_grad", &before_grad);
    let before_flat: Vec<f64> = before_etk_positions
        .iter()
        .flat_map(|p| [p.x, p.y, p.z])
        .collect();
    for (idx, contrib) in mixed_ff.contribs().iter().enumerate() {
        let energy = contrib.get_energy(&before_flat);
        let mut contrib_grad = vec![0.0; before_flat.len()];
        contrib.get_grad(&before_flat, &mut contrib_grad);
        println!("row20_contrib idx={idx} energy={energy:.15} grad={contrib_grad:?}");
    }

    let got_etk = embedder_minimize_with_exp_torsions(&mut positions, &eargs, &params);
    println!("row20_after_etk_minimization planar={}", i32::from(got_etk));
    print_debug_point_positions("row20_after_etk_minimization", &positions);
    let after_etk_positions = point_vectors_to_forcefield_vec3(&positions);
    let mut after_ff = construct_3d_forcefield(&mmat, &after_etk_positions, &etkdg_details);
    after_ff.initialize();
    println!(
        "row20_after_etk_minimization_energy energy={:.15}",
        after_ff.calc_energy_current(None)
    );
    let mut after_grad = vec![0.0; after_etk_positions.len() * 3];
    after_ff.calc_grad_current(&mut after_grad);
    print_debug_flat_vector("row20_after_etk_minimization_grad", &after_grad);
}

#[test]
#[ignore = "debug helper for row-57 ETKDG stage-by-stage parity investigation"]
fn debug_row57_etkdg_stage_trace() {
    let mol = Molecule::from_smiles("CCCCCC(=O)NNC(N)=S")
        .expect("parse")
        .with_hydrogens()
        .expect("add hs");
    let mut params = EmbedParameters::etkdg_v3();
    params.random_seed = 61453;
    params.num_threads = 1;
    params.timeout = 10;

    let mut etkdg_details = CrystalFFDetails::default();
    embedder_init_etkdg(&mol, &params, &mut etkdg_details).expect("init etkdg");

    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let coord_map = params.coord_map.as_ref();
    let setup_ok = embedder_setup_initial_bounds_matrix(
        &mol,
        &mut mmat,
        coord_map,
        &params,
        &mut etkdg_details,
    )
    .expect("bounds setup");
    assert!(setup_ok);

    let mut chiral_centers = Vec::new();
    let mut tetrahedral_carbons = Vec::new();
    embedder_find_chiral_sets(
        &mol,
        &mut chiral_centers,
        &mut tetrahedral_carbons,
        coord_map,
    );

    let mut double_bond_ends = Vec::new();
    let mut stereo_double_bonds = Vec::new();
    embedder_find_double_bonds(
        &mol,
        &mut double_bond_ends,
        &mut stereo_double_bonds,
        coord_map,
    );

    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: Some(&etkdg_details),
        double_bond_ends: Some(&double_bond_ends),
        stereo_double_bonds: &stereo_double_bonds,
    };

    println!(
        "row57_details chiral_centers={:?} tetrahedral_carbons={:?} exp_torsion_atoms={:?} exp_torsion_angles={:?} improper_atoms={:?} bonds={:?} angles={:?} stereo_double_bonds={:?}",
        chiral_centers,
        tetrahedral_carbons,
        etkdg_details.exp_torsion_atoms,
        etkdg_details.exp_torsion_angles,
        etkdg_details.improper_atoms,
        etkdg_details.bonds,
        etkdg_details.angles,
        stereo_double_bonds
    );

    let dim = if params.use_random_coords || !chiral_centers.is_empty() {
        4
    } else {
        3
    };
    let mut positions = vec![vec![0.0; dim]; mol.num_atoms()];
    let mut dist_mat = SymmMatrix::with_value(mol.num_atoms(), 0.0);
    let mut rng = RdkitDistgeomMinStdRand::new_from_embed_points_seed(params.random_seed);

    let got_initial =
        embedder_generate_initial_coords(&mut positions, &eargs, &params, &mut dist_mat, &mut rng)
            .expect("generate initial coords");
    assert!(got_initial);
    print_debug_point_positions("row57_after_generate_initial_coords", &positions);
    let mut initial_ff = construct_distgeom_forcefield(
        &mmat,
        &positions,
        &chiral_centers,
        1.0,
        0.1,
        None,
        params.basin_thresh,
        None,
    );
    initial_ff.initialize();
    println!(
        "row57_before_first_minimization_energy energy={:.15}",
        initial_ff.calc_energy_current(None)
    );
    let mut initial_grad = vec![0.0; positions.len() * dim];
    initial_ff.calc_grad_current(&mut initial_grad);
    print_debug_flat_vector("row57_before_first_minimization_grad", &initial_grad);
    let initial_flat: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();
    for (idx, contrib) in initial_ff.contribs().iter().enumerate() {
        let energy = contrib.get_energy(&initial_flat);
        let mut contrib_grad = vec![0.0; initial_flat.len()];
        contrib.get_grad(&initial_flat, &mut contrib_grad);
        println!("row57_first_contrib idx={idx} energy={energy:.15} grad={contrib_grad:?}");
    }

    let got_first = embedder_first_minimization(&mut positions, &eargs, &params);
    assert!(got_first);
    let mut first_ff = construct_distgeom_forcefield(
        &mmat,
        &positions,
        &chiral_centers,
        1.0,
        0.1,
        None,
        params.basin_thresh,
        None,
    );
    first_ff.initialize();
    let first_energy = first_ff.calc_energy_current(None);
    println!(
        "row57_after_first_minimization energy={:.15} per_atom={:.15}",
        first_energy,
        first_energy / positions.len() as f64
    );
    print_debug_point_positions("row57_after_first_minimization", &positions);

    assert!(embedder_check_tetrahedral_centers(
        &positions, &eargs, &params
    ));
    assert!(embedder_check_chiral_centers(&positions, &eargs, &params));

    if !chiral_centers.is_empty() || params.use_random_coords {
        let got_fourth =
            embedder_minimize_fourth_dimension(&mut positions, &eargs, &mut params, None);
        assert!(got_fourth);
        let fourth_positions = point_vectors_to_forcefield_vec3(&positions);
        let mut fourth_ff = construct_distgeom_forcefield(
            &mmat,
            &positions,
            &chiral_centers,
            0.2,
            1.0,
            None,
            params.basin_thresh,
            None,
        );
        fourth_ff.initialize();
        let fourth_energy = fourth_ff.calc_energy_current(None);
        println!(
            "row57_after_fourth_dimension energy={:.15} per_atom={:.15}",
            fourth_energy,
            fourth_energy / fourth_positions.len() as f64
        );
        print_debug_point_positions("row57_after_fourth_dimension", &positions);
    }

    let before_etk = positions.clone();
    print_debug_point_positions("row57_before_etk_minimization", &before_etk);
    let before_etk_positions = point_vectors_to_forcefield_vec3(&before_etk);
    let mut mixed_ff = construct_3d_forcefield(&mmat, &before_etk_positions, &etkdg_details);
    mixed_ff.initialize();
    println!(
        "row57_before_etk_minimization_energy energy={:.15}",
        mixed_ff.calc_energy_current(None)
    );
    let mut before_grad = vec![0.0; before_etk_positions.len() * 3];
    mixed_ff.calc_grad_current(&mut before_grad);
    print_debug_flat_vector("row57_before_etk_minimization_grad", &before_grad);
    let before_flat: Vec<f64> = before_etk_positions
        .iter()
        .flat_map(|p| [p.x, p.y, p.z])
        .collect();
    for (idx, contrib) in mixed_ff.contribs().iter().enumerate() {
        let energy = contrib.get_energy(&before_flat);
        let mut contrib_grad = vec![0.0; before_flat.len()];
        contrib.get_grad(&before_flat, &mut contrib_grad);
        println!("row57_contrib idx={idx} energy={energy:.15} grad={contrib_grad:?}");
    }

    let got_etk = embedder_minimize_with_exp_torsions(&mut positions, &eargs, &params);
    println!("row57_after_etk_minimization planar={}", i32::from(got_etk));
    print_debug_point_positions("row57_after_etk_minimization", &positions);
    let after_etk_positions = point_vectors_to_forcefield_vec3(&positions);
    let mut after_ff = construct_3d_forcefield(&mmat, &after_etk_positions, &etkdg_details);
    after_ff.initialize();
    println!(
        "row57_after_etk_minimization_energy energy={:.15}",
        after_ff.calc_energy_current(None)
    );
    let mut after_grad = vec![0.0; after_etk_positions.len() * 3];
    after_ff.calc_grad_current(&mut after_grad);
    print_debug_flat_vector("row57_after_etk_minimization_grad", &after_grad);
}

#[test]
#[ignore = "debug helper for row-60 ETKDG stage-by-stage parity investigation"]
fn debug_row60_etkdg_stage_trace() {
    let mol = Molecule::from_smiles("Cc1cc(C(=O)Nc2ccc(Cl)cn2)ccc1[N+](=O)[O-]")
        .expect("parse")
        .with_hydrogens()
        .expect("add hs");
    let mut params = EmbedParameters::etkdg_v3();
    params.random_seed = 61453;
    params.num_threads = 1;
    params.timeout = 10;

    let mut etkdg_details = CrystalFFDetails::default();
    embedder_init_etkdg(&mol, &params, &mut etkdg_details).expect("init etkdg");

    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let coord_map = params.coord_map.as_ref();
    let setup_ok = embedder_setup_initial_bounds_matrix(
        &mol,
        &mut mmat,
        coord_map,
        &params,
        &mut etkdg_details,
    )
    .expect("bounds setup");
    assert!(setup_ok);

    let mut chiral_centers = Vec::new();
    let mut tetrahedral_carbons = Vec::new();
    embedder_find_chiral_sets(
        &mol,
        &mut chiral_centers,
        &mut tetrahedral_carbons,
        coord_map,
    );

    let mut double_bond_ends = Vec::new();
    let mut stereo_double_bonds = Vec::new();
    embedder_find_double_bonds(
        &mol,
        &mut double_bond_ends,
        &mut stereo_double_bonds,
        coord_map,
    );

    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: Some(&etkdg_details),
        double_bond_ends: Some(&double_bond_ends),
        stereo_double_bonds: &stereo_double_bonds,
    };

    println!(
        "row60_details chiral_centers={:?} tetrahedral_carbons={:?} exp_torsion_atoms={:?} exp_torsion_angles={:?} improper_atoms={:?} bonds={:?} angles={:?} stereo_double_bonds={:?}",
        chiral_centers,
        tetrahedral_carbons,
        etkdg_details.exp_torsion_atoms,
        etkdg_details.exp_torsion_angles,
        etkdg_details.improper_atoms,
        etkdg_details.bonds,
        etkdg_details.angles,
        stereo_double_bonds
    );

    let dim = if params.use_random_coords || !chiral_centers.is_empty() {
        4
    } else {
        3
    };
    let mut positions = vec![vec![0.0; dim]; mol.num_atoms()];
    let mut dist_mat = SymmMatrix::with_value(mol.num_atoms(), 0.0);
    let mut rng = RdkitDistgeomMinStdRand::new_from_embed_points_seed(params.random_seed);

    let got_initial =
        embedder_generate_initial_coords(&mut positions, &eargs, &params, &mut dist_mat, &mut rng)
            .expect("generate initial coords");
    assert!(got_initial);
    print_debug_point_positions("row60_after_generate_initial_coords", &positions);

    let got_first = embedder_first_minimization(&mut positions, &eargs, &params);
    assert!(got_first);
    let mut first_ff = construct_distgeom_forcefield(
        &mmat,
        &positions,
        &chiral_centers,
        1.0,
        0.1,
        None,
        params.basin_thresh,
        None,
    );
    first_ff.initialize();
    let first_energy = first_ff.calc_energy_current(None);
    println!(
        "row60_after_first_minimization energy={:.15} per_atom={:.15}",
        first_energy,
        first_energy / positions.len() as f64
    );
    print_debug_point_positions("row60_after_first_minimization", &positions);

    assert!(embedder_check_tetrahedral_centers(
        &positions, &eargs, &params
    ));
    assert!(embedder_check_chiral_centers(&positions, &eargs, &params));

    if !chiral_centers.is_empty() || params.use_random_coords {
        let got_fourth =
            embedder_minimize_fourth_dimension(&mut positions, &eargs, &mut params, None);
        assert!(got_fourth);
        let fourth_positions = point_vectors_to_forcefield_vec3(&positions);
        let mut fourth_ff = construct_distgeom_forcefield(
            &mmat,
            &positions,
            &chiral_centers,
            0.2,
            1.0,
            None,
            params.basin_thresh,
            None,
        );
        fourth_ff.initialize();
        let fourth_energy = fourth_ff.calc_energy_current(None);
        println!(
            "row60_after_fourth_dimension energy={:.15} per_atom={:.15}",
            fourth_energy,
            fourth_energy / fourth_positions.len() as f64
        );
        print_debug_point_positions("row60_after_fourth_dimension", &positions);
    }

    let before_etk = positions.clone();
    print_debug_point_positions("row60_before_etk_minimization", &before_etk);
    let before_etk_positions = point_vectors_to_forcefield_vec3(&before_etk);
    let mut mixed_ff = construct_3d_forcefield(&mmat, &before_etk_positions, &etkdg_details);
    mixed_ff.initialize();
    print_row1_forcefield_contrib_layout(&before_etk_positions, &mmat, &etkdg_details);
    print_row1_constraint_terms(&before_etk_positions, &mmat, &etkdg_details);
    println!(
        "row60_before_etk_minimization_energy energy={:.15}",
        mixed_ff.calc_energy_current(None)
    );
    let mut before_grad = vec![0.0; before_etk_positions.len() * 3];
    mixed_ff.calc_grad_current(&mut before_grad);
    print_debug_flat_vector("row60_before_etk_minimization_grad", &before_grad);
    let before_flat: Vec<f64> = before_etk_positions
        .iter()
        .flat_map(|p| [p.x, p.y, p.z])
        .collect();
    for (idx, contrib) in mixed_ff.contribs().iter().enumerate() {
        let energy = contrib.get_energy(&before_flat);
        let mut contrib_grad = vec![0.0; before_flat.len()];
        contrib.get_grad(&before_flat, &mut contrib_grad);
        println!("row60_contrib idx={idx} energy={energy:.15} grad={contrib_grad:?}");
    }

    let got_etk = embedder_minimize_with_exp_torsions(&mut positions, &eargs, &params);
    println!("row60_after_etk_minimization planar={}", i32::from(got_etk));
    print_debug_point_positions("row60_after_etk_minimization", &positions);
    let after_etk_positions = point_vectors_to_forcefield_vec3(&positions);
    let mut after_ff = construct_3d_forcefield(&mmat, &after_etk_positions, &etkdg_details);
    after_ff.initialize();
    println!(
        "row60_after_etk_minimization_energy energy={:.15}",
        after_ff.calc_energy_current(None)
    );
    let mut after_grad = vec![0.0; after_etk_positions.len() * 3];
    after_ff.calc_grad_current(&mut after_grad);
    print_debug_flat_vector("row60_after_etk_minimization_grad", &after_grad);
}

#[test]
#[ignore = "debug helper for row-89 ETKDG stage-by-stage parity investigation"]
fn debug_row89_etkdg_stage_trace() {
    unsafe {
        std::env::set_var("RDKIT_ROW89_TRACE", "1");
    }
    let mol = Molecule::from_smiles("COC(=O)c4ccccc4(NC(=O)n3c1ccccc1sc2ccccc23)")
        .expect("parse")
        .with_hydrogens()
        .expect("add hs");
    let mut params = EmbedParameters::etkdg_v3();
    params.random_seed = 61453;
    params.num_threads = 1;
    params.timeout = 10;

    let mut etkdg_details = CrystalFFDetails::default();
    embedder_init_etkdg(&mol, &params, &mut etkdg_details).expect("init etkdg");

    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let coord_map = params.coord_map.as_ref();
    let setup_ok = embedder_setup_initial_bounds_matrix(
        &mol,
        &mut mmat,
        coord_map,
        &params,
        &mut etkdg_details,
    )
    .expect("bounds setup");
    assert!(setup_ok);

    let mut chiral_centers = Vec::new();
    let mut tetrahedral_carbons = Vec::new();
    embedder_find_chiral_sets(
        &mol,
        &mut chiral_centers,
        &mut tetrahedral_carbons,
        coord_map,
    );

    let mut double_bond_ends = Vec::new();
    let mut stereo_double_bonds = Vec::new();
    embedder_find_double_bonds(
        &mol,
        &mut double_bond_ends,
        &mut stereo_double_bonds,
        coord_map,
    );

    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: Some(&etkdg_details),
        double_bond_ends: Some(&double_bond_ends),
        stereo_double_bonds: &stereo_double_bonds,
    };

    println!(
        "row89_details num_atoms={} chiral_centers={:?} tetrahedral_carbons={:?} exp_torsion_atoms_len={} improper_atoms_len={} stereo_double_bonds={:?}",
        mol.num_atoms(),
        chiral_centers,
        tetrahedral_carbons,
        etkdg_details.exp_torsion_atoms.len(),
        etkdg_details.improper_atoms.len(),
        stereo_double_bonds
    );

    let dim = if params.use_random_coords || !chiral_centers.is_empty() {
        4
    } else {
        3
    };
    let mut positions = vec![vec![0.0; dim]; mol.num_atoms()];
    let mut dist_mat = SymmMatrix::with_value(mol.num_atoms(), 0.0);
    let mut rng = RdkitDistgeomMinStdRand::new_from_embed_points_seed(params.random_seed);

    let got_initial =
        embedder_generate_initial_coords(&mut positions, &eargs, &params, &mut dist_mat, &mut rng)
            .expect("generate initial coords");
    println!("row89_after_generate_initial_coords ok={got_initial}");
    if got_initial {
        print_debug_point_positions("row89_after_generate_initial_coords", &positions);
    }

    let got_first = if got_initial {
        embedder_first_minimization(&mut positions, &eargs, &params)
    } else {
        false
    };
    println!("row89_after_first_minimization ok={got_first}");
    if got_first {
        print_debug_point_positions("row89_after_first_minimization", &positions);
    }

    let got_tetra = if got_first {
        embedder_check_tetrahedral_centers(&positions, &eargs, &params)
    } else {
        false
    };
    println!("row89_after_tetrahedral_checks ok={got_tetra}");

    let got_chiral = if got_tetra && params.enforce_chirality && !chiral_centers.is_empty() {
        embedder_check_chiral_centers(&positions, &eargs, &params)
    } else {
        got_tetra
    };
    println!("row89_after_chiral_checks ok={got_chiral}");

    let got_fourth = if got_chiral && (!chiral_centers.is_empty() || params.use_random_coords) {
        embedder_minimize_fourth_dimension(&mut positions, &eargs, &mut params, None)
    } else {
        got_chiral
    };
    println!("row89_after_fourth_dimension ok={got_fourth}");
    if got_fourth && (!chiral_centers.is_empty() || params.use_random_coords) {
        print_debug_point_positions("row89_after_fourth_dimension", &positions);
    }

    let got_etk = if got_fourth
        && (params.use_exp_torsion_angle_prefs || params.use_basic_knowledge)
    {
        print_debug_point_positions("row89_before_etk_minimization", &positions);
        let before_etk_positions = point_vectors_to_forcefield_vec3(&positions);
        let mut mixed_ff = construct_3d_forcefield(&mmat, &before_etk_positions, &etkdg_details);
        mixed_ff.initialize();
        println!(
            "row89_before_etk_minimization_energy energy={:.15}",
            mixed_ff.calc_energy_current(None)
        );
        let mut before_grad = vec![0.0; before_etk_positions.len() * 3];
        mixed_ff.calc_grad_current(&mut before_grad);
        print_debug_flat_vector("row89_before_etk_minimization_grad", &before_grad);
        let before_flat: Vec<f64> = before_etk_positions
            .iter()
            .flat_map(|p| [p.x, p.y, p.z])
            .collect();
        for (idx, contrib) in mixed_ff.contribs().iter().enumerate() {
            let energy = contrib.get_energy(&before_flat);
            let mut contrib_grad = vec![0.0; before_flat.len()];
            contrib.get_grad(&before_flat, &mut contrib_grad);
            println!("row89_contrib idx={idx} energy={energy:.15} grad={contrib_grad:?}");
        }
        embedder_minimize_with_exp_torsions(&mut positions, &eargs, &params)
    } else {
        got_fourth
    };
    println!("row89_after_etk_minimization ok={got_etk}");
    if got_etk {
        print_debug_point_positions("row89_after_etk_minimization", &positions);
        let after_etk_positions = point_vectors_to_forcefield_vec3(&positions);
        let mut after_ff = construct_3d_forcefield(&mmat, &after_etk_positions, &etkdg_details);
        after_ff.initialize();
        println!(
            "row89_after_etk_minimization_energy energy={:.15}",
            after_ff.calc_energy_current(None)
        );
        let mut after_grad = vec![0.0; after_etk_positions.len() * 3];
        after_ff.calc_grad_current(&mut after_grad);
        print_debug_flat_vector("row89_after_etk_minimization_grad", &after_grad);
    }
    unsafe {
        std::env::remove_var("RDKIT_ROW89_TRACE");
    }
}

#[test]
#[ignore = "debug helper for row-61 ETKDG stage-by-stage parity investigation"]
fn debug_row61_etkdg_stage_trace() {
    let mol = Molecule::from_smiles("COC(=O)c1cccc([N+](=O)[O-])c1C(=O)OC")
        .expect("parse")
        .with_hydrogens()
        .expect("add hs");
    let mut params = EmbedParameters::etkdg_v3();
    params.random_seed = 61453;
    params.num_threads = 1;
    params.timeout = 10;

    let mut etkdg_details = CrystalFFDetails::default();
    embedder_init_etkdg(&mol, &params, &mut etkdg_details).expect("init etkdg");

    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let coord_map = params.coord_map.as_ref();
    let setup_ok = embedder_setup_initial_bounds_matrix(
        &mol,
        &mut mmat,
        coord_map,
        &params,
        &mut etkdg_details,
    )
    .expect("bounds setup");
    assert!(setup_ok);

    let mut chiral_centers = Vec::new();
    let mut tetrahedral_carbons = Vec::new();
    embedder_find_chiral_sets(
        &mol,
        &mut chiral_centers,
        &mut tetrahedral_carbons,
        coord_map,
    );

    let mut double_bond_ends = Vec::new();
    let mut stereo_double_bonds = Vec::new();
    embedder_find_double_bonds(
        &mol,
        &mut double_bond_ends,
        &mut stereo_double_bonds,
        coord_map,
    );

    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: Some(&etkdg_details),
        double_bond_ends: Some(&double_bond_ends),
        stereo_double_bonds: &stereo_double_bonds,
    };

    println!(
        "row61_details num_atoms={} chiral_centers={:?} tetrahedral_carbons={:?} exp_torsion_atoms_len={} improper_atoms_len={} stereo_double_bonds={:?}",
        mol.num_atoms(),
        chiral_centers,
        tetrahedral_carbons,
        etkdg_details.exp_torsion_atoms.len(),
        etkdg_details.improper_atoms.len(),
        stereo_double_bonds
    );

    let dim = if params.use_random_coords || !chiral_centers.is_empty() {
        4
    } else {
        3
    };
    let mut positions = vec![vec![0.0; dim]; mol.num_atoms()];
    let mut dist_mat = SymmMatrix::with_value(mol.num_atoms(), 0.0);
    let mut rng = RdkitDistgeomMinStdRand::new_from_embed_points_seed(params.random_seed);
    let mut got_initial = false;
    let mut iter = 0_u32;
    while !got_initial && iter < 2 {
        iter += 1;
        println!("row61_debug_iter_start iter={iter}");
        let largest_distance = pick_random_dist_mat_with_rng(&mmat, &mut dist_mat, &mut rng);
        let dist_data = dist_mat.get_data();
        let preview_len = dist_data.len().min(12);
        print_debug_flat_vector(
            &format!("row61_iter{iter}_after_pick_random_distmat"),
            &dist_data[..preview_len],
        );
        println!(
            "row61_iter{iter}_after_pick_random_distmat_largest_distance={largest_distance:.15}"
        );

        got_initial = compute_initial_coords_with_rng(
            &dist_mat,
            &mut positions,
            &mut rng,
            params.rand_neg_eig,
            params.num_zero_fail as usize,
        )
        .expect("compute initial coords");
        println!("row61_iter{iter}_after_generate_initial_coords ok={got_initial}");
        if got_initial {
            print_debug_point_positions(
                &format!("row61_iter{iter}_after_generate_initial_coords"),
                &positions,
            );
        }
    }
    assert!(got_initial, "row61 should get initial coords by iter 2");

    let got_first = embedder_first_minimization(&mut positions, &eargs, &params);
    println!("row61_after_first_minimization ok={got_first}");
    if got_first {
        print_debug_point_positions("row61_after_first_minimization", &positions);
    }
    assert!(got_first);

    let got_tetra = embedder_check_tetrahedral_centers(&positions, &eargs, &params);
    println!("row61_after_tetrahedral_checks ok={got_tetra}");
    assert!(got_tetra);

    let got_chiral = if params.enforce_chirality && !chiral_centers.is_empty() {
        embedder_check_chiral_centers(&positions, &eargs, &params)
    } else {
        got_tetra
    };
    println!("row61_after_chiral_checks ok={got_chiral}");
    assert!(got_chiral);

    let got_fourth = if !chiral_centers.is_empty() || params.use_random_coords {
        embedder_minimize_fourth_dimension(&mut positions, &eargs, &mut params, None)
    } else {
        got_chiral
    };
    println!("row61_after_fourth_dimension ok={got_fourth}");
    if !chiral_centers.is_empty() || params.use_random_coords {
        print_debug_point_positions("row61_after_fourth_dimension", &positions);
    }
    assert!(got_fourth);

    print_debug_point_positions("row61_before_etk_minimization", &positions);
    let before_etk_positions = point_vectors_to_forcefield_vec3(&positions);
    let mut mixed_ff = construct_3d_forcefield(&mmat, &before_etk_positions, &etkdg_details);
    mixed_ff.initialize();
    print_row1_forcefield_contrib_layout(&before_etk_positions, &mmat, &etkdg_details);
    print_row1_constraint_terms(&before_etk_positions, &mmat, &etkdg_details);
    println!(
        "row61_before_etk_minimization_energy energy={:.15}",
        mixed_ff.calc_energy_current(None)
    );
    let mut before_grad = vec![0.0; before_etk_positions.len() * 3];
    mixed_ff.calc_grad_current(&mut before_grad);
    print_debug_flat_vector("row61_before_etk_minimization_grad", &before_grad);
    let before_flat: Vec<f64> = before_etk_positions
        .iter()
        .flat_map(|p| [p.x, p.y, p.z])
        .collect();
    for (idx, contrib) in mixed_ff.contribs().iter().enumerate() {
        let energy = contrib.get_energy(&before_flat);
        let mut contrib_grad = vec![0.0; before_flat.len()];
        contrib.get_grad(&before_flat, &mut contrib_grad);
        println!("row61_contrib idx={idx} energy={energy:.15} grad={contrib_grad:?}");
    }

    let got_etk = embedder_minimize_with_exp_torsions(&mut positions, &eargs, &params);
    println!("row61_after_etk_minimization ok={got_etk}");
    assert!(got_etk);
    print_debug_point_positions("row61_after_etk_minimization", &positions);
    let after_etk_positions = point_vectors_to_forcefield_vec3(&positions);
    let mut after_ff = construct_3d_forcefield(&mmat, &after_etk_positions, &etkdg_details);
    after_ff.initialize();
    println!(
        "row61_after_etk_minimization_energy energy={:.15}",
        after_ff.calc_energy_current(None)
    );
    let mut after_grad = vec![0.0; after_etk_positions.len() * 3];
    after_ff.calc_grad_current(&mut after_grad);
    print_debug_flat_vector("row61_after_etk_minimization_grad", &after_grad);
}

#[test]
#[ignore = "debug helper for row-64 timeout stage timing"]
fn debug_row64_timeout_timing_trace() {
    use std::time::Instant;

    let mol = Molecule::from_smiles(
        "CCCCCCCCCc1ccc(C(=O)Nc2ccc(NC(=O)c3ccc(CCCCCCCCC)cc3)c3c2C(=O)c2ccccc2C3=O)cc1",
    )
    .expect("parse")
    .with_hydrogens()
    .expect("add hs");
    let mut params = EmbedParameters::etkdg_v3();
    params.random_seed = 61453;
    params.num_threads = 1;
    params.timeout = 10;
    params.track_failures = true;
    params
        .failures
        .resize(EmbedFailureCause::EndOfEnum as usize, 0);
    params.failures.fill(0);

    let t0 = Instant::now();
    let mut etkdg_details = CrystalFFDetails::default();
    embedder_init_etkdg(&mol, &params, &mut etkdg_details).expect("init etkdg");
    let t1 = Instant::now();

    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let coord_map = params.coord_map.as_ref();
    let setup_ok = embedder_setup_initial_bounds_matrix(
        &mol,
        &mut mmat,
        coord_map,
        &params,
        &mut etkdg_details,
    )
    .expect("bounds setup");
    let t2 = Instant::now();
    assert!(setup_ok);

    let mut chiral_centers = Vec::new();
    let mut tetrahedral_carbons = Vec::new();
    embedder_find_chiral_sets(
        &mol,
        &mut chiral_centers,
        &mut tetrahedral_carbons,
        coord_map,
    );
    let t3 = Instant::now();

    let mut double_bond_ends = Vec::new();
    let mut stereo_double_bonds = Vec::new();
    embedder_find_double_bonds(
        &mol,
        &mut double_bond_ends,
        &mut stereo_double_bonds,
        coord_map,
    );
    let t4 = Instant::now();

    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: Some(&etkdg_details),
        double_bond_ends: Some(&double_bond_ends),
        stereo_double_bonds: &stereo_double_bonds,
    };
    let dim = if params.use_random_coords || !chiral_centers.is_empty() {
        4
    } else {
        3
    };
    let mut positions = vec![vec![0.0; dim]; mol.num_atoms()];
    let got_coords = embedder_embed_points(&mut positions, eargs, &mut params, 61453, None)
        .expect("embed points");
    let t5 = Instant::now();

    let mut full_params = EmbedParameters::etkdg_v3();
    full_params.random_seed = 61453;
    full_params.num_threads = 1;
    full_params.timeout = 10;
    full_params.track_failures = true;
    full_params
        .failures
        .resize(EmbedFailureCause::EndOfEnum as usize, 0);
    full_params.failures.fill(0);
    let (_embedded, status) = embed_molecule(&mol, &mut full_params).expect("full embed");
    let t6 = Instant::now();

    println!(
        "row64_timing init_etkdg={:.6} bounds={:.6} chiral={:.6} double_bonds={:.6} embed_points={:.6} full_embed={:.6}",
        (t1 - t0).as_secs_f64(),
        (t2 - t1).as_secs_f64(),
        (t3 - t2).as_secs_f64(),
        (t4 - t3).as_secs_f64(),
        (t5 - t4).as_secs_f64(),
        (t6 - t5).as_secs_f64(),
    );
    println!(
        "row64_timing embed_points_ok={got_coords} full_status={status} embed_point_failures={:?} full_failures={:?}",
        params.failures, full_params.failures
    );
}

#[test]
#[ignore = "debug helper for row-64 ETKDG stage-by-stage parity investigation"]
fn debug_row64_etkdg_stage_trace() {
    let mol = Molecule::from_smiles(
        "CCCCCCCCCc1ccc(C(=O)Nc2ccc(NC(=O)c3ccc(CCCCCCCCC)cc3)c3c2C(=O)c2ccccc2C3=O)cc1",
    )
    .expect("parse")
    .with_hydrogens()
    .expect("add hs");
    let mut params = EmbedParameters::etkdg_v3();
    params.random_seed = 61453;
    params.num_threads = 1;
    params.timeout = 10;

    let mut etkdg_details = CrystalFFDetails::default();
    embedder_init_etkdg(&mol, &params, &mut etkdg_details).expect("init etkdg");

    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let coord_map = params.coord_map.as_ref();
    let setup_ok = embedder_setup_initial_bounds_matrix(
        &mol,
        &mut mmat,
        coord_map,
        &params,
        &mut etkdg_details,
    )
    .expect("bounds setup");
    assert!(setup_ok);

    let mut chiral_centers = Vec::new();
    let mut tetrahedral_carbons = Vec::new();
    embedder_find_chiral_sets(
        &mol,
        &mut chiral_centers,
        &mut tetrahedral_carbons,
        coord_map,
    );

    let mut double_bond_ends = Vec::new();
    let mut stereo_double_bonds = Vec::new();
    embedder_find_double_bonds(
        &mol,
        &mut double_bond_ends,
        &mut stereo_double_bonds,
        coord_map,
    );

    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: Some(&etkdg_details),
        double_bond_ends: Some(&double_bond_ends),
        stereo_double_bonds: &stereo_double_bonds,
    };

    println!(
        "row64_details num_atoms={} chiral_centers={:?} tetrahedral_carbons={:?} exp_torsion_atoms_len={} improper_atoms_len={} stereo_double_bonds={:?}",
        mol.num_atoms(),
        chiral_centers,
        tetrahedral_carbons,
        etkdg_details.exp_torsion_atoms.len(),
        etkdg_details.improper_atoms.len(),
        stereo_double_bonds
    );

    let dim = if params.use_random_coords || !chiral_centers.is_empty() {
        4
    } else {
        3
    };
    let mut positions = vec![vec![0.0; dim]; mol.num_atoms()];
    let mut dist_mat = SymmMatrix::with_value(mol.num_atoms(), 0.0);
    let mut rng = RdkitDistgeomMinStdRand::new_from_embed_points_seed(params.random_seed);

    let largest_distance = pick_random_dist_mat_with_rng(&mmat, &mut dist_mat, &mut rng);
    let dist_data = dist_mat.get_data();
    let preview_len = dist_data.len().min(12);
    print_debug_flat_vector("row64_after_pick_random_distmat", &dist_data[..preview_len]);
    println!("row64_after_pick_random_distmat_largest_distance={largest_distance:.15}");

    let got_initial = compute_initial_coords_with_rng(
        &dist_mat,
        &mut positions,
        &mut rng,
        params.rand_neg_eig,
        params.num_zero_fail as usize,
    )
    .expect("compute initial coords");
    println!("row64_after_compute_initial_coords ok={got_initial}");
    if got_initial {
        print_debug_point_positions("row64_after_generate_initial_coords", &positions);
    }

    let got_first = if got_initial {
        embedder_first_minimization(&mut positions, &eargs, &params)
    } else {
        false
    };
    println!("row64_after_first_minimization ok={got_first}");
    if got_first {
        print_debug_point_positions("row64_after_first_minimization", &positions);
    }

    let got_tetra = if got_first {
        embedder_check_tetrahedral_centers(&positions, &eargs, &params)
    } else {
        false
    };
    println!("row64_after_tetrahedral_checks ok={got_tetra}");

    let got_chiral = if got_tetra && params.enforce_chirality && !chiral_centers.is_empty() {
        embedder_check_chiral_centers(&positions, &eargs, &params)
    } else {
        got_tetra
    };
    println!("row64_after_chiral_checks ok={got_chiral}");

    let got_fourth = if got_chiral && (!chiral_centers.is_empty() || params.use_random_coords) {
        embedder_minimize_fourth_dimension(&mut positions, &eargs, &mut params, None)
    } else {
        got_chiral
    };
    println!("row64_after_fourth_dimension ok={got_fourth}");
    if got_fourth && (!chiral_centers.is_empty() || params.use_random_coords) {
        print_debug_point_positions("row64_after_fourth_dimension", &positions);
    }

    let got_etk =
        if got_fourth && (params.use_exp_torsion_angle_prefs || params.use_basic_knowledge) {
            print_debug_point_positions("row64_before_etk_minimization", &positions);
            embedder_minimize_with_exp_torsions(&mut positions, &eargs, &params)
        } else {
            got_fourth
        };
    println!("row64_after_etk_minimization ok={got_etk}");
    if got_etk {
        print_debug_point_positions("row64_after_etk_minimization", &positions);
    }
}

#[test]
#[ignore = "debug helper for row-103 ETKDG stage-by-stage parity investigation"]
fn debug_row103_etkdg_stage_trace() {
    unsafe {
        std::env::set_var("RDKIT_ROW103_TRACE", "1");
    }
    let mol =
        Molecule::from_smiles("O=C(C1=C(C2CC2)N(C3=C4C=CC=NC4=CC=C3)N=C1)NC(N)=N.[H]Cl.[H]O[H]")
            .expect("parse")
            .with_hydrogens()
            .expect("add hs");
    let (fragments, frag_mapping) =
        molecule_fragments_for_embed(&mol, true).expect("fragment split");
    println!(
        "row103_fragment_summary total_atoms={} fragment_count={} frag_mapping={:?}",
        mol.num_atoms(),
        fragments.len(),
        frag_mapping
    );
    for (frag_idx, frag) in fragments.iter().enumerate() {
        println!(
            "row103_fragment frag_idx={} num_atoms={} smiles={}",
            frag_idx,
            frag.num_atoms(),
            frag.to_smiles(false)
                .unwrap_or_else(|_| "<smiles-error>".to_string())
        );
    }
    let mol = fragments
        .iter()
        .max_by_key(|frag| frag.num_atoms())
        .expect("at least one fragment")
        .clone();
    let mut params = EmbedParameters::etkdg_v3();
    params.random_seed = 61453;
    params.num_threads = 1;
    params.timeout = 10;

    let mut etkdg_details = CrystalFFDetails::default();
    embedder_init_etkdg(&mol, &params, &mut etkdg_details).expect("init etkdg");

    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let coord_map = params.coord_map.as_ref();
    let setup_ok = embedder_setup_initial_bounds_matrix(
        &mol,
        &mut mmat,
        coord_map,
        &params,
        &mut etkdg_details,
    )
    .expect("bounds setup");
    assert!(setup_ok);

    let mut chiral_centers = Vec::new();
    let mut tetrahedral_carbons = Vec::new();
    embedder_find_chiral_sets(
        &mol,
        &mut chiral_centers,
        &mut tetrahedral_carbons,
        coord_map,
    );

    let mut double_bond_ends = Vec::new();
    let mut stereo_double_bonds = Vec::new();
    embedder_find_double_bonds(
        &mol,
        &mut double_bond_ends,
        &mut stereo_double_bonds,
        coord_map,
    );

    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: Some(&etkdg_details),
        double_bond_ends: Some(&double_bond_ends),
        stereo_double_bonds: &stereo_double_bonds,
    };

    println!(
        "row103_details num_atoms={} chiral_centers={:?} tetrahedral_carbons={:?} exp_torsion_atoms_len={} improper_atoms_len={} stereo_double_bonds={:?}",
        mol.num_atoms(),
        chiral_centers,
        tetrahedral_carbons,
        etkdg_details.exp_torsion_atoms.len(),
        etkdg_details.improper_atoms.len(),
        stereo_double_bonds
    );

    let dim = if params.use_random_coords || !chiral_centers.is_empty() {
        4
    } else {
        3
    };
    let mut positions = vec![vec![0.0; dim]; mol.num_atoms()];
    let mut dist_mat = SymmMatrix::with_value(mol.num_atoms(), 0.0);
    let mut rng = RdkitDistgeomMinStdRand::new_from_embed_points_seed(params.random_seed);

    let largest_distance = pick_random_dist_mat_with_rng(&mmat, &mut dist_mat, &mut rng);
    let dist_data = dist_mat.get_data();
    let preview_len = dist_data.len().min(12);
    print_debug_flat_vector(
        "row103_after_pick_random_distmat",
        &dist_data[..preview_len],
    );
    println!("row103_after_pick_random_distmat_largest_distance={largest_distance:.15}");

    let got_initial = compute_initial_coords_with_rng(
        &dist_mat,
        &mut positions,
        &mut rng,
        params.rand_neg_eig,
        params.num_zero_fail as usize,
    )
    .expect("compute initial coords");
    println!("row103_after_compute_initial_coords ok={got_initial}");
    if got_initial {
        print_debug_point_positions("row103_after_generate_initial_coords", &positions);
    }

    let got_first = if got_initial {
        embedder_first_minimization(&mut positions, &eargs, &params)
    } else {
        false
    };
    println!("row103_after_first_minimization ok={got_first}");
    if got_first {
        print_debug_point_positions("row103_after_first_minimization", &positions);
    }

    let got_tetra = if got_first {
        embedder_check_tetrahedral_centers(&positions, &eargs, &params)
    } else {
        false
    };
    println!("row103_after_tetrahedral_checks ok={got_tetra}");

    let got_chiral = if got_tetra && params.enforce_chirality && !chiral_centers.is_empty() {
        embedder_check_chiral_centers(&positions, &eargs, &params)
    } else {
        got_tetra
    };
    println!("row103_after_chiral_checks ok={got_chiral}");

    let got_fourth = if got_chiral && (!chiral_centers.is_empty() || params.use_random_coords) {
        embedder_minimize_fourth_dimension(&mut positions, &eargs, &mut params, None)
    } else {
        got_chiral
    };
    println!("row103_after_fourth_dimension ok={got_fourth}");
    if got_fourth && (!chiral_centers.is_empty() || params.use_random_coords) {
        print_debug_point_positions("row103_after_fourth_dimension", &positions);
    }

    let got_etk =
        if got_fourth && (params.use_exp_torsion_angle_prefs || params.use_basic_knowledge) {
            print_debug_point_positions("row103_before_etk_minimization", &positions);
            embedder_minimize_with_exp_torsions(&mut positions, &eargs, &params)
        } else {
            got_fourth
        };
    println!("row103_after_etk_minimization ok={got_etk}");
    if got_etk {
        print_debug_point_positions("row103_after_etk_minimization", &positions);
    }
    unsafe {
        std::env::remove_var("RDKIT_ROW103_TRACE");
    }
}

#[test]
#[ignore = "source-bisection probe for audited ETKDGv3 experimental torsion minimization"]
fn debug_audited_etkdgv3_forcefield_first_step() {
    for smiles in [
        "O=C1Nc2ccc(Cl)cc2C(c2ccccc2Cl)=NC1O",
        "CN1C2CCC1CC1(CN=C(c3cn(C)c4ccccc34)O1)C2",
    ] {
        let mol = Molecule::from_smiles(smiles)
            .expect("parse")
            .with_hydrogens()
            .expect("add hs");
        let mut params = EmbedParameters::etkdg_v3();
        params.random_seed = 61453;
        params.max_iterations = 3;
        params.num_threads = 1;
        params.timeout = 0;

        let staged = embedder_stage_coords_for_test(
            &mol,
            &mut params,
            EmbedderTestStage::FourthDimensionCleaned,
        )
        .expect("stage embed");
        let positions: Vec<ForceFieldVec3> = staged.conformers_3d()[0]
            .coordinates()
            .iter()
            .map(|point| ForceFieldVec3::new(point[0], point[1], point[2]))
            .collect();

        let mut details = CrystalFFDetails::default();
        embedder_init_etkdg(&mol, &params, &mut details).expect("init etkdg");
        let mut bounds = BoundsMatrix::new(mol.num_atoms());
        assert!(
            embedder_setup_initial_bounds_matrix(
                &mol,
                &mut bounds,
                params.coord_map.as_ref(),
                &params,
                &mut details,
            )
            .expect("bounds setup")
        );

        let mut field = construct_3d_forcefield(&bounds, &positions, &details);
        field.initialize();
        let mut contrib_energies = Vec::new();
        let total_energy = field.calc_energy_current(Some(&mut contrib_energies));
        let mut gradient = vec![0.0; 3 * positions.len()];
        field.calc_grad_current(&mut gradient);

        println!("probe molecule={smiles}");
        println!(
            "probe boundary=etk_forcefield_initial total_energy_bits={:x} contrib_count={}",
            total_energy.to_bits(),
            contrib_energies.len()
        );
        for (idx, energy) in contrib_energies.iter().enumerate() {
            println!(
                "probe contribution={idx} energy_bits={:x}",
                energy.to_bits()
            );
        }
        for (idx, (position, grad)) in positions
            .iter()
            .flat_map(|point| [point.x, point.y, point.z])
            .zip(&gradient)
            .enumerate()
        {
            println!(
                "probe initial={idx} pos_bits={:x} grad_bits={:x}",
                position.to_bits(),
                grad.to_bits()
            );
        }

        let status = field.minimize(1, params.optimizer_force_tol, 1.0e-6);
        println!("probe boundary=etk_bfgs_first_step status={status}");
        for (idx, position) in field
            .positions()
            .iter()
            .flat_map(|point| [point.x, point.y, point.z])
            .enumerate()
        {
            println!("probe first_step={idx} pos_bits={:x}", position.to_bits());
        }
    }
}

#[test]
fn embed_parameter_presets_match_rdkit_global_parameters() {
    assert_embed_parameters_preset(
        &EmbedParameters::kdg(),
        0,
        1,
        -1,
        true,
        false,
        2.0,
        true,
        1,
        1e-3,
        false,
        true,
        false,
        true,
        false,
        5.0,
        -1.0,
        true,
        1,
        true,
        false,
        false,
        false,
        0,
    );
    assert_embed_parameters_preset(
        &EmbedParameters::etdg(),
        0,
        1,
        -1,
        true,
        false,
        2.0,
        true,
        1,
        1e-3,
        false,
        false,
        true,
        false,
        false,
        5.0,
        -1.0,
        true,
        1,
        true,
        false,
        false,
        false,
        0,
    );
    assert_embed_parameters_preset(
        &EmbedParameters::etdg_v2(),
        0,
        1,
        -1,
        true,
        false,
        2.0,
        true,
        1,
        1e-3,
        false,
        false,
        true,
        false,
        false,
        5.0,
        -1.0,
        true,
        2,
        true,
        false,
        false,
        false,
        0,
    );
    assert_embed_parameters_preset(
        &EmbedParameters::etkdg(),
        0,
        1,
        -1,
        true,
        false,
        2.0,
        true,
        1,
        1e-3,
        false,
        true,
        true,
        true,
        false,
        5.0,
        -1.0,
        true,
        1,
        true,
        false,
        false,
        false,
        0,
    );
    assert_embed_parameters_preset(
        &EmbedParameters::etkdg_v2(),
        0,
        1,
        -1,
        true,
        false,
        2.0,
        true,
        1,
        1e-3,
        false,
        true,
        true,
        true,
        false,
        5.0,
        -1.0,
        true,
        2,
        true,
        false,
        false,
        false,
        0,
    );
    assert_embed_parameters_preset(
        &EmbedParameters::etkdg_v3(),
        0,
        1,
        -1,
        true,
        false,
        2.0,
        true,
        1,
        1e-3,
        false,
        true,
        true,
        true,
        false,
        5.0,
        -1.0,
        true,
        2,
        true,
        false,
        true,
        true,
        0,
    );
    assert_embed_parameters_preset(
        &EmbedParameters::sr_etkdg_v3(),
        0,
        1,
        -1,
        true,
        false,
        2.0,
        true,
        1,
        1e-3,
        false,
        true,
        true,
        true,
        false,
        5.0,
        -1.0,
        true,
        2,
        true,
        true,
        false,
        false,
        0,
    );
}

#[test]
fn update_embed_parameters_from_json_empty_string_is_noop() {
    let mut params = EmbedParameters::etkdg_v3();
    let before = params.clone();

    params.update_from_json("").unwrap();

    assert_eq!(params.max_iterations, before.max_iterations);
    assert_eq!(params.et_version, before.et_version);
    assert_eq!(
        params.use_macrocycle_torsions,
        before.use_macrocycle_torsions
    );
    assert_eq!(params.coord_map, before.coord_map);
}

#[test]
fn update_embed_parameters_from_json_updates_only_present_scalar_fields() {
    let mut params = EmbedParameters::default();

    params
        .update_from_json(
            r#"{
              "maxIterations": 23,
              "randomSeed": 17,
              "useRandomCoords": true,
              "optimizerForceTol": 0.25,
              "useMacrocycleTorsions": true,
              "unknownIgnored": false
            }"#,
        )
        .unwrap();

    assert_eq!(params.max_iterations, 23);
    assert_eq!(params.random_seed, 17);
    assert!(params.use_random_coords);
    assert_eq!(params.optimizer_force_tol, 0.25);
    assert!(params.use_macrocycle_torsions);
    assert_eq!(params.num_threads, 1);
    assert!(params.clear_confs);
    assert!(!params.use_macrocycle14config);
}

#[test]
fn update_embed_parameters_from_json_updates_all_rdkit_macro_fields() {
    let mut params = EmbedParameters::default();

    params
        .update_from_json(
            r#"{
              "basinThresh": 1.25,
              "boundsMatForceScaling": 2.5,
              "boxSizeMult": 3.5,
              "clearConfs": false,
              "embedFragmentsSeparately": false,
              "enableSequentialRandomSeeds": true,
              "enforceChirality": false,
              "ETversion": 3,
              "forceTransAmides": false,
              "ignoreSmoothingFailures": true,
              "maxIterations": 101,
              "numThreads": 4,
              "numZeroFail": 5,
              "onlyHeavyAtomsForRMS": false,
              "optimizerForceTol": 0.125,
              "pruneRmsThresh": 0.75,
              "randNegEig": false,
              "randomSeed": 99,
              "symmetrizeConjugatedTerminalGroupsForPruning": false,
              "timeout": 44,
              "trackFailures": true,
              "useBasicKnowledge": true,
              "useExpTorsionAnglePrefs": true,
              "useMacrocycle14config": true,
              "useMacrocycleTorsions": true,
              "useRandomCoords": true,
              "useSmallRingTorsions": true,
              "useSymmetryForPruning": false,
              "verbose": true
            }"#,
        )
        .unwrap();

    assert_eq!(params.basin_thresh, 1.25);
    assert_eq!(params.bounds_mat_force_scaling, 2.5);
    assert_eq!(params.box_size_mult, 3.5);
    assert!(!params.clear_confs);
    assert!(!params.embed_fragments_separately);
    assert!(params.enable_sequential_random_seeds);
    assert!(!params.enforce_chirality);
    assert_eq!(params.et_version, 3);
    assert!(!params.force_trans_amides);
    assert!(params.ignore_smoothing_failures);
    assert_eq!(params.max_iterations, 101);
    assert_eq!(params.num_threads, 4);
    assert_eq!(params.num_zero_fail, 5);
    assert!(!params.only_heavy_atoms_for_rms);
    assert_eq!(params.optimizer_force_tol, 0.125);
    assert_eq!(params.prune_rms_thresh, 0.75);
    assert!(!params.rand_neg_eig);
    assert_eq!(params.random_seed, 99);
    assert!(!params.symmetrize_conjugated_terminal_groups_for_pruning);
    assert_eq!(params.timeout, 44);
    assert!(params.track_failures);
    assert!(params.use_basic_knowledge);
    assert!(params.use_exp_torsion_angle_prefs);
    assert!(params.use_macrocycle14config);
    assert!(params.use_macrocycle_torsions);
    assert!(params.use_random_coords);
    assert!(params.use_small_ring_torsions);
    assert!(!params.use_symmetry_for_pruning);
    assert!(params.verbose);
}

#[test]
fn update_embed_parameters_from_json_updates_coord_map() {
    let mut params = EmbedParameters::default();

    params
        .update_from_json(
            r#"{
              "coordMap": {
                "2": [1.0, 2.0, 3.0],
                "5": ["4.5", "5.5", "6.5"]
              }
            }"#,
        )
        .unwrap();

    let coord_map = params.coord_map.as_ref().unwrap();
    assert_eq!(coord_map.len(), 2);
    assert_eq!(coord_map[&2], ForceFieldVec3::new(1.0, 2.0, 3.0));
    assert_eq!(coord_map[&5], ForceFieldVec3::new(4.5, 5.5, 6.5));
}

#[test]
fn update_embed_parameters_from_json_rejects_invalid_json_and_field_types() {
    let mut params = EmbedParameters::default();

    assert!(params.update_from_json("{").is_err());
    assert!(params.update_from_json(r#"{"maxIterations": -1}"#).is_err());
    assert!(
        params
            .update_from_json(r#"{"coordMap": {"x": [1.0, 2.0, 3.0]}}"#)
            .is_err()
    );
    assert!(
        params
            .update_from_json(r#"{"coordMap": {"1": [1.0, 2.0]}}"#)
            .is_err()
    );
}

#[test]
fn embed_parameters_to_json_matches_rdkit_without_maps() {
    let json = EmbedParameters::kdg().to_json();
    let expected = r#"{"basinThresh":"5","boundsMatForceScaling":"1","boxSizeMult":"2","clearConfs":"true","embedFragmentsSeparately":"true","enableSequentialRandomSeeds":"false","enforceChirality":"true","ETversion":"1","forceTransAmides":"true","ignoreSmoothingFailures":"false","maxIterations":"0","numThreads":"1","numZeroFail":"1","onlyHeavyAtomsForRMS":"true","optimizerForceTol":"0.001","pruneRmsThresh":"-1","randNegEig":"true","randomSeed":"-1","symmetrizeConjugatedTerminalGroupsForPruning":"true","timeout":"0","trackFailures":"false","useBasicKnowledge":"true","useExpTorsionAnglePrefs":"false","useMacrocycle14config":"false","useMacrocycleTorsions":"false","useRandomCoords":"false","useSmallRingTorsions":"false","useSymmetryForPruning":"true","verbose":"false"}"#;

    assert_eq!(json, expected);
}

#[test]
fn embed_parameters_to_json_matches_rdkit_coord_map_shape() {
    let mut params = EmbedParameters::kdg();
    params
        .coord_map
        .get_or_insert_with(BTreeMap::new)
        .insert(3, ForceFieldVec3::new(1.1, 2.2, 3.3));

    let json = params.to_json();
    let expected = r#"{"basinThresh":"5","boundsMatForceScaling":"1","boxSizeMult":"2","clearConfs":"true","embedFragmentsSeparately":"true","enableSequentialRandomSeeds":"false","enforceChirality":"true","ETversion":"1","forceTransAmides":"true","ignoreSmoothingFailures":"false","maxIterations":"0","numThreads":"1","numZeroFail":"1","onlyHeavyAtomsForRMS":"true","optimizerForceTol":"0.001","pruneRmsThresh":"-1","randNegEig":"true","randomSeed":"-1","symmetrizeConjugatedTerminalGroupsForPruning":"true","timeout":"0","trackFailures":"false","useBasicKnowledge":"true","useExpTorsionAnglePrefs":"false","useMacrocycle14config":"false","useMacrocycleTorsions":"false","useRandomCoords":"false","useSmallRingTorsions":"false","useSymmetryForPruning":"true","verbose":"false","coordMap":{"3":["1.100000","2.200000","3.300000"]}}"#;

    assert_eq!(json, expected);
}

#[test]
fn embed_parameters_to_json_includes_bounds_matrix_values() {
    let mut params = EmbedParameters::kdg();
    let mut bounds = BoundsMatrix::new(2);
    bounds.set_val(0, 0, 0.0);
    bounds.set_val(0, 1, 2.5);
    bounds.set_val(1, 0, 1.25);
    bounds.set_val(1, 1, 0.0);
    params.bounds_mat = Some(Arc::new(bounds));

    let json = params.to_json();

    assert!(json.ends_with(r#","boundsMatrix":[["0","2.5"],["1.25","0"]]}"#));
}

#[test]
fn embed_parameters_to_json_round_trips_through_update_from_json() {
    let params = EmbedParameters::etkdg_v3();
    let json = params.to_json();
    let mut round_trip = EmbedParameters::default();

    round_trip.update_from_json(&json).unwrap();

    assert_eq!(round_trip.to_json(), json);
}

#[test]
fn embedder_have_opposite_sign_matches_rdkit_signbit_xor() {
    assert!(have_opposite_sign(-1.0, 1.0));
    assert!(have_opposite_sign(1.0, -1.0));
    assert!(!have_opposite_sign(1.0, 2.0));
    assert!(!have_opposite_sign(-1.0, -2.0));
    assert!(have_opposite_sign(-0.0, 0.0));
    assert!(!have_opposite_sign(-0.0, -1.0));
    assert!(!have_opposite_sign(0.0, 1.0));

    let negative_nan = f64::from_bits(0xfff8_0000_0000_0001);
    let positive_nan = f64::from_bits(0x7ff8_0000_0000_0001);
    assert!(have_opposite_sign(negative_nan, positive_nan));
    assert!(!have_opposite_sign(negative_nan, -1.0));
}

#[test]
fn embedder_failure_mutex_returns_single_global_mutex() {
    failmutex_create();

    let first = get_fail_mutex() as *const _;
    let second = get_fail_mutex() as *const _;
    let direct = failmutex_get() as *const _;

    assert_eq!(first, second);
    assert_eq!(first, direct);
}

#[test]
fn embedder_failure_mutex_serializes_concurrent_access() {
    let inside = Arc::new(AtomicUsize::new(0));
    let max_inside = Arc::new(AtomicUsize::new(0));
    let mut handles = Vec::new();

    for _ in 0..8 {
        let inside = Arc::clone(&inside);
        let max_inside = Arc::clone(&max_inside);
        handles.push(thread::spawn(move || {
            let _guard = get_fail_mutex()
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            let now_inside = inside.fetch_add(1, Ordering::SeqCst) + 1;
            max_inside.fetch_max(now_inside, Ordering::SeqCst);
            thread::yield_now();
            inside.fetch_sub(1, Ordering::SeqCst);
        }));
    }

    for handle in handles {
        handle.join().unwrap();
    }

    assert_eq!(max_inside.load(Ordering::SeqCst), 1);
}

#[test]
fn embedder_volume_test_accepts_well_separated_tetrahedral_center() {
    let positions = vec![
        ForceFieldVec3::new(0.0, 0.0, 0.0),
        ForceFieldVec3::new(1.0, 0.0, 0.0),
        ForceFieldVec3::new(0.0, 1.0, 0.0),
        ForceFieldVec3::new(0.0, 0.0, 1.0),
        ForceFieldVec3::new(-1.0, -1.0, -1.0),
    ];
    let chiral_set = ChiralSet::with_default_structure_flags(0, 1, 2, 3, 4, -1.0, 1.0);

    assert!(embedder_volume_test(&chiral_set, &positions));
}

#[test]
fn embedder_volume_test_rejects_flat_or_low_volume_center() {
    let flat_positions = vec![
        ForceFieldVec3::new(0.0, 0.0, 0.0),
        ForceFieldVec3::new(1.0, 0.0, 0.0),
        ForceFieldVec3::new(0.0, 1.0, 0.0),
        ForceFieldVec3::new(0.0, 0.0, 1.0),
        ForceFieldVec3::new(-1.0, -1.0, 0.0),
    ];
    let chiral_set = ChiralSet::with_default_structure_flags(0, 1, 2, 3, 4, -1.0, 1.0);

    assert!(!embedder_volume_test(&chiral_set, &flat_positions));

    let low_volume_positions = vec![
        ForceFieldVec3::new(0.0, 0.0, 0.0),
        ForceFieldVec3::new(1.0, 0.0, 0.0),
        ForceFieldVec3::new(0.0, 1.0, 0.0),
        ForceFieldVec3::new(0.0, 0.0, 1.0),
        ForceFieldVec3::new(-1.0, -1.0, -0.3),
    ];
    assert!(!embedder_volume_test(&chiral_set, &low_volume_positions));
}

#[test]
fn embedder_volume_test_uses_fused_small_ring_relaxed_threshold() {
    let positions = vec![
        ForceFieldVec3::new(0.0, 0.0, 0.0),
        ForceFieldVec3::new(1.0, 0.0, 0.0),
        ForceFieldVec3::new(0.0, 1.0, 0.0),
        ForceFieldVec3::new(0.0, 0.0, 1.0),
        ForceFieldVec3::new(-1.0, -1.0, -0.3),
    ];
    let regular = ChiralSet::with_default_structure_flags(0, 1, 2, 3, 4, -1.0, 1.0);
    let fused = ChiralSet::new(
        0,
        1,
        2,
        3,
        4,
        -1.0,
        1.0,
        ChiralSetStructureFlags::InFusedSmallRings as u64,
    );

    assert!(!embedder_volume_test(&regular, &positions));
    assert!(embedder_volume_test(&fused, &positions));
}

#[test]
fn embedder_same_side_matches_plane_side_and_tolerance_rules() {
    let v1 = ForceFieldVec3::new(0.0, 0.0, 0.0);
    let v2 = ForceFieldVec3::new(1.0, 0.0, 0.0);
    let v3 = ForceFieldVec3::new(0.0, 1.0, 0.0);
    let v4 = ForceFieldVec3::new(0.0, 0.0, 1.0);

    assert!(embedder_same_side(
        v1,
        v2,
        v3,
        v4,
        ForceFieldVec3::new(0.25, 0.25, 0.5),
        0.1
    ));
    assert!(!embedder_same_side(
        v1,
        v2,
        v3,
        v4,
        ForceFieldVec3::new(0.25, 0.25, -0.5),
        0.1
    ));
    assert!(!embedder_same_side(
        v1,
        v2,
        v3,
        v4,
        ForceFieldVec3::new(0.25, 0.25, 0.05),
        0.1
    ));
    assert!(!embedder_same_side(
        v1,
        v2,
        v3,
        ForceFieldVec3::new(0.0, 0.0, 0.05),
        ForceFieldVec3::new(0.25, 0.25, 0.5),
        0.1
    ));
}

#[test]
fn embedder_center_in_volume_checks_all_four_faces() {
    let positions = vec![
        ForceFieldVec3::new(0.1, 0.1, 0.1),
        ForceFieldVec3::new(0.0, 0.0, 0.0),
        ForceFieldVec3::new(1.0, 0.0, 0.0),
        ForceFieldVec3::new(0.0, 1.0, 0.0),
        ForceFieldVec3::new(0.0, 0.0, 1.0),
        ForceFieldVec3::new(2.0, 2.0, 2.0),
    ];

    assert!(embedder_center_in_volume_indices(
        0, 1, 2, 3, 4, &positions, 0.01
    ));
    assert!(!embedder_center_in_volume_indices(
        5, 1, 2, 3, 4, &positions, 0.01
    ));
}

#[test]
fn embedder_center_in_volume_chiral_set_overload_handles_three_coordinate_centers() {
    let positions = vec![
        ForceFieldVec3::new(10.0, 10.0, 10.0),
        ForceFieldVec3::new(0.0, 0.0, 0.0),
        ForceFieldVec3::new(1.0, 0.0, 0.0),
        ForceFieldVec3::new(0.0, 1.0, 0.0),
    ];
    let three_coordinate = ChiralSet::with_default_structure_flags(0, 1, 2, 3, 0, -1.0, 1.0);

    assert!(embedder_center_in_volume(
        &three_coordinate,
        &positions,
        0.1
    ));
}

#[test]
fn embedder_center_in_volume_chiral_set_overload_delegates_indices() {
    let positions = vec![
        ForceFieldVec3::new(0.1, 0.1, 0.1),
        ForceFieldVec3::new(0.0, 0.0, 0.0),
        ForceFieldVec3::new(1.0, 0.0, 0.0),
        ForceFieldVec3::new(0.0, 1.0, 0.0),
        ForceFieldVec3::new(0.0, 0.0, 1.0),
    ];
    let chiral_set = ChiralSet::with_default_structure_flags(0, 1, 2, 3, 4, -1.0, 1.0);

    assert!(embedder_center_in_volume(&chiral_set, &positions, 0.01));
}

#[test]
fn embedder_bounds_fulfilled_matches_rdkit_tolerance_rule() {
    let mut mmat = BoundsMatrix::new(3);
    mmat.set_lower(0, 1, 0.9).expect("set lower");
    mmat.set_upper(0, 1, 1.1).expect("set upper");
    mmat.set_lower(0, 2, 1.0).expect("set lower");
    mmat.set_upper(0, 2, 2.0).expect("set upper");
    mmat.set_lower(1, 2, 1.0).expect("set lower");
    mmat.set_upper(1, 2, 2.0).expect("set upper");
    let atoms = [0, 1, 2];
    let ok_positions = vec![
        ForceFieldVec3::new(0.0, 0.0, 0.0),
        ForceFieldVec3::new(1.15, 0.0, 0.0),
        ForceFieldVec3::new(0.0, 1.2, 0.0),
    ];
    let bad_positions = vec![
        ForceFieldVec3::new(0.0, 0.0, 0.0),
        ForceFieldVec3::new(1.25, 0.0, 0.0),
        ForceFieldVec3::new(0.0, 1.2, 0.0),
    ];

    assert!(embedder_bounds_fulfilled(&atoms, &mmat, &ok_positions));
    assert!(!embedder_bounds_fulfilled(&atoms, &mmat, &bad_positions));
    assert!(embedder_bounds_fulfilled(&[], &mmat, &bad_positions));
    assert!(embedder_bounds_fulfilled(&[0], &mmat, &bad_positions));
}

#[test]
fn embedder_generate_initial_coords_random_branch_applies_coord_map_and_zeroes_higher_dimensions() {
    let mut mmat = BoundsMatrix::new(3);
    mmat.set_lower(0, 1, 1.0).expect("set lower");
    mmat.set_upper(0, 1, 1.0).expect("set upper");
    let chiral_centers: VectChiralSet = Vec::new();
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };
    let mut coord_map = BTreeMap::new();
    coord_map.insert(1, ForceFieldVec3::new(7.0, 8.0, 9.0));
    let mut params = EmbedParameters::default();
    params.use_random_coords = true;
    params.box_size_mult = 2.0;
    params.coord_map = Some(coord_map);
    let mut positions = vec![vec![0.0; 4], vec![0.0; 4], vec![0.0; 4]];
    let mut expected_positions = positions.clone();
    let mut expected_rng = RdkitDistgeomMinStdRand::new(123);
    assert!(compute_random_coords_with_rng(
        &mut expected_positions,
        10.0,
        &mut expected_rng
    ));
    expected_positions[1][0] = 7.0;
    expected_positions[1][1] = 8.0;
    expected_positions[1][2] = 9.0;
    expected_positions[1][3] = 0.0;

    let mut dist_mat = SymmMatrix::new(3);
    let mut rng = RdkitDistgeomMinStdRand::new(123);

    assert!(
        embedder_generate_initial_coords(&mut positions, &eargs, &params, &mut dist_mat, &mut rng)
            .expect("generate initial coords")
    );
    assert_eq!(positions, expected_positions);
}

#[test]
fn embedder_generate_initial_coords_distance_matrix_branch_uses_bounds_matrix() {
    let mut mmat = BoundsMatrix::new(2);
    mmat.set_lower(0, 1, 2.0).expect("set lower");
    mmat.set_upper(0, 1, 2.0).expect("set upper");
    let chiral_centers: VectChiralSet = Vec::new();
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };
    let params = EmbedParameters::default();
    let mut positions = vec![vec![0.0; 3], vec![0.0; 3]];
    let mut dist_mat = SymmMatrix::new(2);
    let mut rng = RdkitDistgeomMinStdRand::new(222);

    assert!(
        embedder_generate_initial_coords(&mut positions, &eargs, &params, &mut dist_mat, &mut rng)
            .expect("generate initial coords")
    );
    assert_eq!(dist_mat.get_val(0, 1), 2.0);
    let distance = ((positions[0][0] - positions[1][0]).powi(2)
        + (positions[0][1] - positions[1][1]).powi(2)
        + (positions[0][2] - positions[1][2]).powi(2))
    .sqrt();
    assert!((distance - 2.0).abs() < 1.0e-6);
}

#[test]
fn embedder_first_minimization_keeps_exact_satisfied_two_point_bounds() {
    let mut mmat = BoundsMatrix::new(2);
    mmat.set_lower(0, 1, 1.0).expect("set lower");
    mmat.set_upper(0, 1, 1.0).expect("set upper");
    let chiral_centers: VectChiralSet = Vec::new();
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };
    let params = EmbedParameters::default();
    let mut positions = vec![vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0]];

    assert!(embedder_first_minimization(&mut positions, &eargs, &params));
    assert_eq!(positions, vec![vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0]]);
}

#[test]
fn embedder_check_tetrahedral_centers_requires_volume_and_center_in_volume() {
    let mmat = BoundsMatrix::new(5);
    let tet_set: ChiralSetPtr = Arc::new(ChiralSet::with_default_structure_flags(
        0, 1, 2, 3, 4, -1.0, 1.0,
    ));
    let chiral_centers: VectChiralSet = Vec::new();
    let tetrahedral_carbons: VectChiralSet = vec![Arc::clone(&tet_set)];
    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };
    let params = EmbedParameters::default();
    let good_positions = vec![
        vec![0.0, 0.0, 0.0],
        vec![1.0, 1.0, 1.0],
        vec![1.0, -1.0, -1.0],
        vec![-1.0, 1.0, -1.0],
        vec![-1.0, -1.0, 1.0],
    ];
    let outside_positions = vec![
        vec![4.0, 4.0, 4.0],
        vec![1.0, 1.0, 1.0],
        vec![1.0, -1.0, -1.0],
        vec![-1.0, 1.0, -1.0],
        vec![-1.0, -1.0, 1.0],
    ];

    assert!(embedder_check_tetrahedral_centers(
        &good_positions,
        &eargs,
        &params
    ));
    assert!(!embedder_check_tetrahedral_centers(
        &outside_positions,
        &eargs,
        &params
    ));
}

#[test]
fn embedder_check_chiral_centers_matches_rdkit_volume_bound_rule() {
    let mmat = BoundsMatrix::new(5);
    let good_set: ChiralSetPtr = Arc::new(ChiralSet::with_default_structure_flags(
        0, 1, 2, 3, 4, 0.5, 2.0,
    ));
    let failing_set: ChiralSetPtr = Arc::new(ChiralSet::with_default_structure_flags(
        0, 1, 2, 3, 4, 2.0, 3.0,
    ));
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let positions = vec![
        vec![0.0, 0.0, 0.0],
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
        vec![0.0, 0.0, 0.0],
    ];
    let params = EmbedParameters::default();
    let good_chiral_centers: VectChiralSet = vec![good_set];
    let good_args = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &good_chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };
    let failing_chiral_centers: VectChiralSet = vec![failing_set];
    let failing_args = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &failing_chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };

    assert!(embedder_check_chiral_centers(
        &positions, &good_args, &params
    ));
    assert!(!embedder_check_chiral_centers(
        &positions,
        &failing_args,
        &params
    ));
}

#[test]
fn embedder_minimize_fourth_dimension_keeps_exact_satisfied_two_point_bounds() {
    let mut mmat = BoundsMatrix::new(2);
    mmat.set_lower(0, 1, 1.0).expect("set lower");
    mmat.set_upper(0, 1, 1.0).expect("set upper");
    let chiral_centers: VectChiralSet = Vec::new();
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };
    let mut params = EmbedParameters::default();
    let mut positions = vec![vec![0.0, 0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0, 0.0]];

    assert!(embedder_minimize_fourth_dimension(
        &mut positions,
        &eargs,
        &mut params,
        None
    ));
    assert_eq!(
        positions,
        vec![vec![0.0, 0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0, 0.0]]
    );
}

#[test]
fn embedder_minimize_fourth_dimension_preserves_random_coord_map_fixed_points() {
    let mut mmat = BoundsMatrix::new(2);
    mmat.set_lower(0, 1, 1.0).expect("set lower");
    mmat.set_upper(0, 1, 1.0).expect("set upper");
    let chiral_centers: VectChiralSet = Vec::new();
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };
    let mut coord_map = BTreeMap::new();
    coord_map.insert(0, ForceFieldVec3::new(7.0, 8.0, 9.0));
    let mut params = EmbedParameters::default();
    params.use_random_coords = true;
    params.coord_map = Some(coord_map);
    let mut positions = vec![vec![7.0, 8.0, 9.0, 0.0], vec![10.0, 8.0, 9.0, 3.0]];

    assert!(embedder_minimize_fourth_dimension(
        &mut positions,
        &eargs,
        &mut params,
        None
    ));
    assert_eq!(positions[0], vec![7.0, 8.0, 9.0, 0.0]);
}

#[test]
fn embedder_minimize_fourth_dimension_returns_false_after_timeout_deadline() {
    let mut mmat = BoundsMatrix::new(2);
    mmat.set_lower(0, 1, 1.0).expect("set lower");
    mmat.set_upper(0, 1, 1.0).expect("set upper");
    let chiral_centers: VectChiralSet = Vec::new();
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };
    let mut params = EmbedParameters::default();
    let mut positions = vec![vec![0.0, 0.0, 0.0, 10.0], vec![3.0, 0.0, 0.0, -10.0]];

    assert!(!embedder_minimize_fourth_dimension(
        &mut positions,
        &eargs,
        &mut params,
        Some(Instant::now() - Duration::from_secs(1))
    ));
}

#[test]
fn embedder_minimize_with_exp_torsions_plain_etdg_preserves_random_coord_map_fixed_points() {
    let mut mmat = BoundsMatrix::new(2);
    mmat.set_lower(0, 1, 1.0).expect("set lower");
    mmat.set_upper(0, 1, 1.0).expect("set upper");
    let chiral_centers: VectChiralSet = Vec::new();
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let details = CrystalFFDetails::default();
    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: Some(&details),
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };
    let mut coord_map = BTreeMap::new();
    coord_map.insert(0, ForceFieldVec3::new(3.0, 4.0, 5.0));
    let mut params = EmbedParameters::default();
    params.use_random_coords = true;
    params.coord_map = Some(coord_map);
    let mut positions = vec![vec![3.0, 4.0, 5.0, 9.0], vec![8.0, 4.0, 5.0, -2.0]];

    assert!(embedder_minimize_with_exp_torsions(
        &mut positions,
        &eargs,
        &params
    ));
    assert_eq!(&positions[0][..3], &[3.0, 4.0, 5.0]);
    assert_eq!(positions[0][3], 9.0);
    assert_eq!(positions[1][3], -2.0);
}

#[test]
fn embedder_minimize_with_exp_torsions_basic_knowledge_accepts_empty_cpci() {
    let mut mmat = BoundsMatrix::new(2);
    mmat.set_lower(0, 1, 1.0).expect("set lower");
    mmat.set_upper(0, 1, 1.0).expect("set upper");
    let chiral_centers: VectChiralSet = Vec::new();
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let details = CrystalFFDetails::default();
    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: Some(&details),
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };
    let mut params = EmbedParameters::default();
    params.use_basic_knowledge = true;
    params.cpci = Some(BTreeMap::new());
    let mut positions = vec![vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0]];

    assert!(embedder_minimize_with_exp_torsions(
        &mut positions,
        &eargs,
        &params
    ));
}

#[test]
#[should_panic(expected = "bogus etkdgDetails pointer")]
fn embedder_minimize_with_exp_torsions_requires_etkdg_details() {
    let mmat = BoundsMatrix::new(1);
    let chiral_centers: VectChiralSet = Vec::new();
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };
    let params = EmbedParameters::default();
    let mut positions = vec![vec![0.0, 0.0, 0.0]];

    let _ = embedder_minimize_with_exp_torsions(&mut positions, &eargs, &params);
}

#[test]
fn embedder_double_bond_geometry_checks_rejects_linear_arrangement() {
    let mmat = BoundsMatrix::new(3);
    let chiral_centers: VectChiralSet = Vec::new();
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let double_bond_ends = vec![(0_usize, 1_usize, 2_usize)];
    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: Some(&double_bond_ends),
        stereo_double_bonds: &[],
    };
    let mut params = EmbedParameters::default();
    let linear = vec![
        vec![0.0, 0.0, 0.0],
        vec![1.0, 0.0, 0.0],
        vec![2.0, 0.0, 0.0],
    ];
    let bent = vec![
        vec![0.0, 0.0, 0.0],
        vec![1.0, 0.0, 0.0],
        vec![1.0, 1.0, 0.0],
    ];

    assert!(!embedder_double_bond_geometry_checks(
        &linear,
        &eargs,
        &mut params,
        1.0e-3
    ));
    assert!(embedder_double_bond_geometry_checks(
        &bent,
        &eargs,
        &mut params,
        1.0e-3
    ));
}

#[test]
fn embedder_double_bond_geometry_checks_accepts_missing_double_bond_ends() {
    let mmat = BoundsMatrix::new(3);
    let chiral_centers: VectChiralSet = Vec::new();
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };
    let mut params = EmbedParameters::default();
    let positions = vec![
        vec![0.0, 0.0, 0.0],
        vec![1.0, 0.0, 0.0],
        vec![2.0, 0.0, 0.0],
    ];

    assert!(embedder_double_bond_geometry_checks(
        &positions,
        &eargs,
        &mut params,
        1.0e-3
    ));
}

#[test]
fn embedder_double_bond_stereo_checks_uses_dihedral_sign_rule() {
    let mmat = BoundsMatrix::new(4);
    let chiral_centers: VectChiralSet = Vec::new();
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let trans_bonds = vec![(vec![0_usize, 1_usize, 2_usize, 3_usize], 1_i32)];
    let cis_bonds = vec![(vec![0_usize, 1_usize, 2_usize, 3_usize], -1_i32)];
    let trans_args = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &trans_bonds,
    };
    let cis_args = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &cis_bonds,
    };
    let mut params = EmbedParameters::default();
    let cis_positions = vec![
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 0.0],
        vec![1.0, 0.0, 0.0],
        vec![1.0, 1.0, 0.0],
    ];
    let trans_positions = vec![
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 0.0],
        vec![1.0, 0.0, 0.0],
        vec![1.0, -1.0, 0.0],
    ];

    assert!(embedder_double_bond_stereo_checks(
        &cis_positions,
        &cis_args,
        &mut params
    ));
    assert!(!embedder_double_bond_stereo_checks(
        &cis_positions,
        &trans_args,
        &mut params
    ));
    assert!(embedder_double_bond_stereo_checks(
        &trans_positions,
        &trans_args,
        &mut params
    ));
}

fn broad_bounds_matrix(size: usize) -> BoundsMatrix {
    let mut mmat = BoundsMatrix::new(size);
    for i in 1..size {
        for j in 0..i {
            mmat.set_lower(i, j, 0.0).expect("set lower");
            mmat.set_upper(i, j, 100.0).expect("set upper");
        }
    }
    mmat
}

#[test]
fn embedder_final_chiral_checks_accepts_valid_final_chirality() {
    let mmat = broad_bounds_matrix(5);
    let chiral_set: ChiralSetPtr = Arc::new(ChiralSet::with_default_structure_flags(
        0, 1, 2, 3, 4, -100.0, 100.0,
    ));
    let chiral_centers: VectChiralSet = vec![chiral_set];
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };
    let mut params = EmbedParameters::default();
    let mut positions = vec![
        vec![0.25, 0.25, 0.25],
        vec![0.0, 0.0, 0.0],
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ];

    assert!(embedder_final_chiral_checks(
        &mut positions,
        &eargs,
        &mut params
    ));
}

#[test]
fn embedder_final_chiral_checks_tracks_volume_failure() {
    let mmat = broad_bounds_matrix(5);
    let chiral_set: ChiralSetPtr = Arc::new(ChiralSet::with_default_structure_flags(
        0, 1, 2, 3, 4, 2.0, 3.0,
    ));
    let chiral_centers: VectChiralSet = vec![chiral_set];
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };
    let mut params = EmbedParameters {
        track_failures: true,
        failures: vec![0; EmbedFailureCause::EndOfEnum as usize],
        ..EmbedParameters::default()
    };
    let mut positions = vec![
        vec![0.1, 0.1, 0.1],
        vec![0.0, 0.0, 0.0],
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ];

    assert!(!embedder_final_chiral_checks(
        &mut positions,
        &eargs,
        &mut params
    ));
    assert_eq!(
        params.failures[EmbedFailureCause::CheckChiralCenters2 as usize],
        1
    );
}

#[test]
fn embedder_final_chiral_checks_tracks_bounds_failure() {
    let mut mmat = broad_bounds_matrix(5);
    mmat.set_upper(0, 1, 0.05).expect("set upper");
    let chiral_set: ChiralSetPtr = Arc::new(ChiralSet::with_default_structure_flags(
        0, 1, 2, 3, 4, -100.0, 100.0,
    ));
    let chiral_centers: VectChiralSet = vec![chiral_set];
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };
    let mut params = EmbedParameters {
        track_failures: true,
        failures: vec![0; EmbedFailureCause::EndOfEnum as usize],
        ..EmbedParameters::default()
    };
    let mut positions = vec![
        vec![0.25, 0.25, 0.25],
        vec![0.0, 0.0, 0.0],
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ];

    assert!(!embedder_final_chiral_checks(
        &mut positions,
        &eargs,
        &mut params
    ));
    assert_eq!(
        params.failures[EmbedFailureCause::FinalChiralBounds as usize],
        1
    );
}

#[test]
fn embedder_final_chiral_checks_tracks_center_in_volume_failure() {
    let mmat = broad_bounds_matrix(5);
    let chiral_set: ChiralSetPtr = Arc::new(ChiralSet::with_default_structure_flags(
        0, 1, 2, 3, 4, -100.0, 100.0,
    ));
    let chiral_centers: VectChiralSet = vec![chiral_set];
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };
    let mut params = EmbedParameters {
        track_failures: true,
        failures: vec![0; EmbedFailureCause::EndOfEnum as usize],
        ..EmbedParameters::default()
    };
    let mut positions = vec![
        vec![4.0, 4.0, 4.0],
        vec![0.0, 0.0, 0.0],
        vec![1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ];

    assert!(!embedder_final_chiral_checks(
        &mut positions,
        &eargs,
        &mut params
    ));
    assert_eq!(
        params.failures[EmbedFailureCause::FinalCenterInVolume as usize],
        1
    );
}

#[test]
fn embedder_embed_points_sets_default_iterations_and_runs_callback() {
    EMBED_POINTS_CALLBACK_COUNT.store(0, Ordering::SeqCst);
    let mut mmat = BoundsMatrix::new(2);
    mmat.set_lower(0, 1, 1.0).expect("set lower");
    mmat.set_upper(0, 1, 1.0).expect("set upper");
    let chiral_centers: VectChiralSet = Vec::new();
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };
    let mut coord_map = BTreeMap::new();
    coord_map.insert(0, ForceFieldVec3::new(0.0, 0.0, 0.0));
    coord_map.insert(1, ForceFieldVec3::new(1.0, 0.0, 0.0));
    let mut params = EmbedParameters {
        callback: Some(embed_points_test_callback),
        use_random_coords: true,
        coord_map: Some(coord_map),
        ..EmbedParameters::default()
    };
    let mut positions = vec![vec![0.0; 3], vec![0.0; 3]];

    assert!(
        embedder_embed_points(&mut positions, eargs, &mut params, 0, None).expect("embed points")
    );
    assert_eq!(params.max_iterations, 20);
    assert_eq!(EMBED_POINTS_CALLBACK_COUNT.load(Ordering::SeqCst), 1);
}

#[test]
fn embedder_embed_points_seed_zero_is_local_and_reproducible() {
    let mut mmat = BoundsMatrix::new(3);
    for i in 1..3 {
        for j in 0..i {
            mmat.set_lower(i, j, 1.0).expect("set lower");
            mmat.set_upper(i, j, 2.0).expect("set upper");
        }
    }
    let chiral_centers: VectChiralSet = Vec::new();
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let mut params_a = EmbedParameters::default();
    let mut params_b = EmbedParameters::default();
    let mut positions_a = vec![vec![0.0; 3], vec![0.0; 3], vec![0.0; 3]];
    let mut positions_b = positions_a.clone();
    let eargs_a = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };
    let eargs_b = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };

    assert!(
        embedder_embed_points(&mut positions_a, eargs_a, &mut params_a, 0, None)
            .expect("embed points")
    );
    assert!(
        embedder_embed_points(&mut positions_b, eargs_b, &mut params_b, 0, None)
            .expect("embed points")
    );
    assert_eq!(positions_a, positions_b);
}

#[test]
fn embedder_embed_points_timeout_before_first_iteration_returns_false_without_callback() {
    let mmat = BoundsMatrix::new(2);
    let chiral_centers: VectChiralSet = Vec::new();
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: None,
        stereo_double_bonds: &[],
    };
    let mut params = EmbedParameters::default();
    let mut positions = vec![vec![0.0; 3], vec![0.0; 3]];

    assert!(
        !embedder_embed_points(
            &mut positions,
            eargs,
            &mut params,
            1,
            Some(Instant::now() - Duration::from_secs(1))
        )
        .expect("embed points")
    );
    assert_eq!(params.max_iterations, 20);
}

#[test]
fn embedder_embed_points_tracks_linear_double_bond_failure() {
    let mmat = broad_bounds_matrix(3);
    let chiral_centers: VectChiralSet = Vec::new();
    let tetrahedral_carbons: VectChiralSet = Vec::new();
    let double_bond_ends = vec![(0_usize, 1_usize, 2_usize)];
    let eargs = EmbedArgs {
        mmat: &mmat,
        chiral_centers: &chiral_centers,
        tetrahedral_carbons: &tetrahedral_carbons,
        etkdg_details: None,
        double_bond_ends: Some(&double_bond_ends),
        stereo_double_bonds: &[],
    };
    let mut coord_map = BTreeMap::new();
    coord_map.insert(0, ForceFieldVec3::new(0.0, 0.0, 0.0));
    coord_map.insert(1, ForceFieldVec3::new(1.0, 0.0, 0.0));
    coord_map.insert(2, ForceFieldVec3::new(2.0, 0.0, 0.0));
    let mut params = EmbedParameters {
        max_iterations: 1,
        use_random_coords: true,
        coord_map: Some(coord_map),
        track_failures: true,
        failures: vec![0; EmbedFailureCause::EndOfEnum as usize],
        ..EmbedParameters::default()
    };
    let mut positions = vec![vec![0.0; 3], vec![0.0; 3], vec![0.0; 3]];

    assert!(
        !embedder_embed_points(&mut positions, eargs, &mut params, 7, None).expect("embed points")
    );
    assert_eq!(
        params.failures[EmbedFailureCause::LinearDoubleBond as usize],
        1
    );
}

#[test]
fn embedder_find_double_bonds_collects_ends_and_clears_outputs() {
    let mut builder = MoleculeBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::C));
    let a3 = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(BondSpec::new(a0, a1, BondOrder::Double))
        .expect("double");
    builder
        .add_bond(BondSpec::new(a0, a2, BondOrder::Single))
        .expect("single left");
    builder
        .add_bond(BondSpec::new(a1, a3, BondOrder::Single))
        .expect("single right");
    let mol = builder.build().expect("mol");
    let mut double_bond_ends = vec![(99, 98, 97)];
    let mut stereo_double_bonds = vec![(vec![9, 8, 7, 6], -1)];

    embedder_find_double_bonds(&mol, &mut double_bond_ends, &mut stereo_double_bonds, None);

    assert_eq!(double_bond_ends, vec![(2, 0, 1), (3, 1, 0)]);
    assert!(stereo_double_bonds.is_empty());
}

#[test]
fn embedder_find_double_bonds_skips_non_single_neighbor_only_for_degree_two_atom() {
    let mut builder = MoleculeBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::C));
    let a3 = builder.add_atom(AtomSpec::new(Element::C));
    let a4 = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(BondSpec::new(a0, a1, BondOrder::Double))
        .expect("central double");
    builder
        .add_bond(BondSpec::new(a0, a2, BondOrder::Double))
        .expect("degree two non-single left");
    builder
        .add_bond(BondSpec::new(a1, a3, BondOrder::Double))
        .expect("degree three non-single right");
    builder
        .add_bond(BondSpec::new(a1, a4, BondOrder::Single))
        .expect("degree three single right");
    let mol = builder.build().expect("mol");
    let mut double_bond_ends = Vec::new();
    let mut stereo_double_bonds = Vec::new();

    embedder_find_double_bonds(&mol, &mut double_bond_ends, &mut stereo_double_bonds, None);

    assert!(double_bond_ends.contains(&(3, 1, 0)));
    assert!(double_bond_ends.contains(&(4, 1, 0)));
    assert!(!double_bond_ends.contains(&(2, 0, 1)));
}

#[test]
fn embedder_find_double_bonds_collects_stereo_sign_and_honors_coord_map_skip() {
    let mut builder = MoleculeBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::C));
    let a3 = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(
            BondSpec::new(a0, a1, BondOrder::Double)
                .with_stereo(BondStereo::Z)
                .with_stereo_atoms(a2, a3),
        )
        .expect("stereo double");
    builder
        .add_bond(BondSpec::new(a0, a2, BondOrder::Single))
        .expect("single left");
    builder
        .add_bond(BondSpec::new(a1, a3, BondOrder::Single))
        .expect("single right");
    let mol = builder.build().expect("mol");
    let mut double_bond_ends = Vec::new();
    let mut stereo_double_bonds = Vec::new();

    embedder_find_double_bonds(&mol, &mut double_bond_ends, &mut stereo_double_bonds, None);
    assert_eq!(stereo_double_bonds, vec![(vec![2, 0, 1, 3], -1)]);

    let mut coord_map = BTreeMap::new();
    coord_map.insert(2, ForceFieldVec3::new(0.0, 0.0, 0.0));
    embedder_find_double_bonds(
        &mol,
        &mut double_bond_ends,
        &mut stereo_double_bonds,
        Some(&coord_map),
    );
    assert_eq!(stereo_double_bonds, vec![(vec![2, 0, 1, 3], -1)]);

    coord_map.insert(3, ForceFieldVec3::new(1.0, 0.0, 0.0));
    embedder_find_double_bonds(
        &mol,
        &mut double_bond_ends,
        &mut stereo_double_bonds,
        Some(&coord_map),
    );
    assert!(stereo_double_bonds.is_empty());
}

#[test]
fn embedder_find_double_bonds_uses_positive_sign_for_trans_and_e() {
    for stereo in [BondStereo::Trans, BondStereo::E] {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        let a2 = builder.add_atom(AtomSpec::new(Element::C));
        let a3 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(
                BondSpec::new(a0, a1, BondOrder::Double)
                    .with_stereo(stereo)
                    .with_stereo_atoms(a2, a3),
            )
            .expect("stereo double");
        builder
            .add_bond(BondSpec::new(a0, a2, BondOrder::Single))
            .expect("single left");
        builder
            .add_bond(BondSpec::new(a1, a3, BondOrder::Single))
            .expect("single right");
        let mol = builder.build().expect("mol");
        let mut double_bond_ends = Vec::new();
        let mut stereo_double_bonds = Vec::new();

        embedder_find_double_bonds(&mol, &mut double_bond_ends, &mut stereo_double_bonds, None);

        assert_eq!(stereo_double_bonds, vec![(vec![2, 0, 1, 3], 1)]);
    }
}

#[test]
fn embedder_find_chiral_sets_collects_tagged_tetrahedral_center_bounds() {
    let mut builder = MoleculeBuilder::new();
    let center =
        builder.add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCcw));
    let n1 = builder.add_atom(AtomSpec::new(Element::C));
    let n2 = builder.add_atom(AtomSpec::new(Element::C));
    let n3 = builder.add_atom(AtomSpec::new(Element::C));
    let n4 = builder.add_atom(AtomSpec::new(Element::C));
    for nbr in [n1, n2, n3, n4] {
        builder
            .add_bond(BondSpec::new(center, nbr, BondOrder::Single))
            .expect("bond");
    }
    let mol = builder.build().expect("mol");
    let mut chiral_centers = Vec::new();
    let mut tetrahedral_centers = Vec::new();

    embedder_find_chiral_sets(&mol, &mut chiral_centers, &mut tetrahedral_centers, None);

    assert_eq!(chiral_centers.len(), 1);
    assert!(tetrahedral_centers.is_empty());
    let cset = &chiral_centers[0];
    assert_eq!(
        (cset.idx0, cset.idx1, cset.idx2, cset.idx3, cset.idx4),
        (0, 1, 2, 3, 4)
    );
    assert_eq!(cset.volume_lower_bound, 5.0);
    assert_eq!(cset.volume_upper_bound, 100.0);
}

#[test]
fn embedder_find_chiral_sets_uses_center_as_fourth_neighbor_for_three_coordinate_tagged_atom() {
    let mut builder = MoleculeBuilder::new();
    let center =
        builder.add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw));
    let n1 = builder.add_atom(AtomSpec::new(Element::C));
    let n2 = builder.add_atom(AtomSpec::new(Element::C));
    let n3 = builder.add_atom(AtomSpec::new(Element::C));
    for nbr in [n1, n2, n3] {
        builder
            .add_bond(BondSpec::new(center, nbr, BondOrder::Single))
            .expect("bond");
    }
    let mol = builder.build().expect("mol");
    let mut chiral_centers = Vec::new();
    let mut tetrahedral_centers = Vec::new();

    embedder_find_chiral_sets(&mol, &mut chiral_centers, &mut tetrahedral_centers, None);

    assert_eq!(chiral_centers.len(), 1);
    let cset = &chiral_centers[0];
    assert_eq!(
        (cset.idx0, cset.idx1, cset.idx2, cset.idx3, cset.idx4),
        (0, 1, 2, 3, 0)
    );
    assert_eq!(cset.volume_lower_bound, -100.0);
    assert_eq!(cset.volume_upper_bound, -2.0);
}

#[test]
fn embedder_find_chiral_sets_collects_unmarked_tetrahedral_c_or_n_unless_coord_mapped() {
    let mut builder = MoleculeBuilder::new();
    let center = builder.add_atom(AtomSpec::new(Element::N));
    let ligands = [
        builder.add_atom(AtomSpec::new(Element::C)),
        builder.add_atom(AtomSpec::new(Element::C)),
        builder.add_atom(AtomSpec::new(Element::C)),
        builder.add_atom(AtomSpec::new(Element::C)),
    ];
    for nbr in ligands {
        builder
            .add_bond(BondSpec::new(center, nbr, BondOrder::Single))
            .expect("bond");
    }
    let mol = builder.build().expect("mol");
    let mut chiral_centers = Vec::new();
    let mut tetrahedral_centers = Vec::new();

    embedder_find_chiral_sets(&mol, &mut chiral_centers, &mut tetrahedral_centers, None);
    assert!(chiral_centers.is_empty());
    assert_eq!(tetrahedral_centers.len(), 1);
    assert_eq!(tetrahedral_centers[0].volume_lower_bound, 0.0);
    assert_eq!(tetrahedral_centers[0].volume_upper_bound, 0.0);

    let mut coord_map = BTreeMap::new();
    coord_map.insert(center.index() as i32, ForceFieldVec3::new(0.0, 0.0, 0.0));
    embedder_find_chiral_sets(
        &mol,
        &mut chiral_centers,
        &mut tetrahedral_centers,
        Some(&coord_map),
    );
    assert_eq!(tetrahedral_centers.len(), 1);
}

#[test]
fn embedder_find_chiral_sets_skips_hydrogen_even_when_tagged() {
    let mut builder = MoleculeBuilder::new();
    let center =
        builder.add_atom(AtomSpec::new(Element::H).with_chiral_tag(ChiralTag::TetrahedralCcw));
    let ligands = [
        builder.add_atom(AtomSpec::new(Element::C)),
        builder.add_atom(AtomSpec::new(Element::C)),
        builder.add_atom(AtomSpec::new(Element::C)),
    ];
    for nbr in ligands {
        builder
            .add_bond(BondSpec::new(center, nbr, BondOrder::Single))
            .expect("bond");
    }
    let mol = builder.build().expect("mol");
    let mut chiral_centers = Vec::new();
    let mut tetrahedral_centers = Vec::new();

    embedder_find_chiral_sets(&mol, &mut chiral_centers, &mut tetrahedral_centers, None);

    assert!(chiral_centers.is_empty());
    assert!(tetrahedral_centers.is_empty());
}

#[test]
fn embedder_find_chiral_sets_collects_atropisomer_chiral_set_with_sorted_neighbors() {
    let mut builder = MoleculeBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let left_hi = builder.add_atom(AtomSpec::new(Element::C));
    let left_lo = builder.add_atom(AtomSpec::new(Element::C));
    let right = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(
            BondSpec::new(a0, a1, BondOrder::Single)
                .with_stereo(BondStereo::AtropCcw)
                .with_stereo_atoms(left_lo, right),
        )
        .expect("atrop");
    builder
        .add_bond(BondSpec::new(a0, left_hi, BondOrder::Single))
        .expect("left high");
    builder
        .add_bond(BondSpec::new(a0, left_lo, BondOrder::Single))
        .expect("left low");
    builder
        .add_bond(BondSpec::new(a1, right, BondOrder::Single))
        .expect("right");
    let mol = builder.build().expect("mol");
    let mut chiral_centers = Vec::new();
    let mut tetrahedral_centers = Vec::new();

    embedder_find_chiral_sets(&mol, &mut chiral_centers, &mut tetrahedral_centers, None);

    assert_eq!(chiral_centers.len(), 1);
    let cset = &chiral_centers[0];
    assert_eq!(
        (cset.idx0, cset.idx1, cset.idx2, cset.idx3, cset.idx4),
        (0, 2, 3, 4, 0)
    );
    assert_eq!(cset.volume_lower_bound, -100.0);
    assert_eq!(cset.volume_upper_bound, -1.0);
}

#[test]
fn embedder_adjust_bounds_mat_from_coord_map_sets_exact_pair_distances() {
    let mut mmat = BoundsMatrix::new(4);
    let original_03_upper = mmat.get_upper(0, 3);
    let original_03_lower = mmat.get_lower(0, 3);
    let mut coord_map = BTreeMap::new();
    coord_map.insert(2, ForceFieldVec3::new(0.0, 0.0, 0.0));
    coord_map.insert(0, ForceFieldVec3::new(3.0, 4.0, 0.0));
    coord_map.insert(1, ForceFieldVec3::new(3.0, 4.0, 12.0));

    embedder_adjust_bounds_mat_from_coord_map(&mut mmat, 4, &coord_map).expect("coord map bounds");

    assert_eq!(mmat.get_upper(0, 2), 5.0);
    assert_eq!(mmat.get_lower(0, 2), 5.0);
    assert_eq!(mmat.get_upper(0, 1), 12.0);
    assert_eq!(mmat.get_lower(0, 1), 12.0);
    assert_eq!(mmat.get_upper(1, 2), 13.0);
    assert_eq!(mmat.get_lower(1, 2), 13.0);
    assert_eq!(mmat.get_upper(0, 3), original_03_upper);
    assert_eq!(mmat.get_lower(0, 3), original_03_lower);
}

#[test]
fn embedder_adjust_bounds_mat_from_coord_map_accepts_empty_and_singleton_maps() {
    let mut mmat = BoundsMatrix::new(2);
    let before = (mmat.get_upper(0, 1), mmat.get_lower(0, 1));
    let empty = BTreeMap::new();
    embedder_adjust_bounds_mat_from_coord_map(&mut mmat, 2, &empty)
        .expect("empty coord map bounds");
    assert_eq!((mmat.get_upper(0, 1), mmat.get_lower(0, 1)), before);

    let mut singleton = BTreeMap::new();
    singleton.insert(0, ForceFieldVec3::new(1.0, 2.0, 3.0));
    embedder_adjust_bounds_mat_from_coord_map(&mut mmat, 2, &singleton)
        .expect("singleton coord map bounds");
    assert_eq!((mmat.get_upper(0, 1), mmat.get_lower(0, 1)), before);
}

#[test]
fn embedder_init_etkdg_without_knowledge_sets_only_force_scaling() {
    let mol = Molecule::from_smiles("CCO").expect("mol");
    let mut details = CrystalFFDetails {
        atom_nums: vec![99],
        bounds_mat_force_scaling: 0.25,
        ..CrystalFFDetails::default()
    };
    let params = EmbedParameters {
        bounds_mat_force_scaling: 2.5,
        ..EmbedParameters::default()
    };

    embedder_init_etkdg(&mol, &params, &mut details).expect("init");

    assert_eq!(details.atom_nums, vec![99]);
    assert_eq!(details.bounds_mat_force_scaling, 2.5);
}

#[test]
fn embedder_init_etkdg_with_basic_knowledge_populates_atom_numbers() {
    let mol = Molecule::from_smiles("CCO").expect("mol");
    let mut details = CrystalFFDetails::default();
    let params = EmbedParameters {
        use_basic_knowledge: true,
        bounds_mat_force_scaling: 3.0,
        ..EmbedParameters::default()
    };

    embedder_init_etkdg(&mol, &params, &mut details).expect("init");

    assert_eq!(details.atom_nums, vec![6, 6, 8]);
    assert_eq!(details.bounds_mat_force_scaling, 3.0);
}

#[test]
fn embedder_init_etkdg_propagates_empty_molecule_failure() {
    let mol = MoleculeBuilder::new().build().expect("empty");
    let mut details = CrystalFFDetails::default();
    let params = EmbedParameters {
        use_exp_torsion_angle_prefs: true,
        ..EmbedParameters::default()
    };

    let err = embedder_init_etkdg(&mol, &params, &mut details).expect_err("empty error");

    assert!(matches!(
        err,
        DgBoundsError::GenerationFailed(message) if message == "RDKit CrystalFF molecule has no atoms"
    ));
}

#[test]
fn embedder_setup_initial_bounds_matrix_sets_topological_bounds_and_smooths() {
    let mol = Molecule::from_smiles("CCO").expect("mol");
    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let mut details = CrystalFFDetails::default();
    let params = EmbedParameters::default();

    let ok = embedder_setup_initial_bounds_matrix(&mol, &mut mmat, None, &params, &mut details)
        .expect("setup");

    assert!(ok);
    assert!(mmat.check_valid());
    assert!(details.bonds.is_empty());
    assert!(details.angles.is_empty());
    assert!(mmat.get_upper(0, 1) < DEFAULT_UPPER);
}

#[test]
fn embedder_setup_initial_bounds_matrix_records_etkdg_bonds_and_angles() {
    let mol = Molecule::from_smiles("CCO").expect("mol");
    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let mut details = CrystalFFDetails::default();
    let params = EmbedParameters {
        use_basic_knowledge: true,
        ..EmbedParameters::default()
    };

    let ok = embedder_setup_initial_bounds_matrix(&mol, &mut mmat, None, &params, &mut details)
        .expect("setup");

    assert!(ok);
    assert_eq!(details.bonds.len(), mol.num_bonds());
    assert!(!details.angles.is_empty());
}

#[test]
fn embedder_setup_initial_bounds_matrix_applies_coord_map_exact_bounds() {
    let mol = Molecule::from_smiles("CCO").expect("mol");
    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let mut details = CrystalFFDetails::default();
    let mut coord_map = BTreeMap::new();
    coord_map.insert(0, ForceFieldVec3::new(0.0, 0.0, 0.0));
    coord_map.insert(1, ForceFieldVec3::new(1.0, 0.0, 0.0));
    let params = EmbedParameters::default();

    let ok = embedder_setup_initial_bounds_matrix(
        &mol,
        &mut mmat,
        Some(&coord_map),
        &params,
        &mut details,
    )
    .expect("setup");

    assert!(ok);
    assert_eq!(mmat.get_upper(0, 1), 1.0);
    assert_eq!(mmat.get_lower(0, 1), 1.0);
}

#[test]
fn embedder_setup_initial_bounds_matrix_propagates_topology_errors() {
    let mol = MoleculeBuilder::new().build().expect("empty");
    let mut mmat = BoundsMatrix::new(0);
    let mut details = CrystalFFDetails::default();
    let params = EmbedParameters::default();

    let err = embedder_setup_initial_bounds_matrix(&mol, &mut mmat, None, &params, &mut details)
        .expect_err("empty molecule");

    assert!(matches!(
        err,
        DgBoundsError::GenerationFailed(message) if message == "molecule has no atoms"
    ));
}

#[test]
fn embedder_fill_atom_positions_copies_positions_in_match_order() {
    let mol = Molecule::from_smiles("CCO").expect("mol");
    let conf = Conformer3D::new(
        0,
        vec![[0.0, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]],
        true,
    );
    let mut pts = vec![ForceFieldVec3::default(); 3];

    embedder_fill_atom_positions(&mut pts, &conf, &mol, &[2, 0, 1]);

    assert_eq!(pts[0], ForceFieldVec3::new(6.0, 7.0, 8.0));
    assert_eq!(pts[1], ForceFieldVec3::new(0.0, 1.0, 2.0));
    assert_eq!(pts[2], ForceFieldVec3::new(3.0, 4.0, 5.0));
}

#[test]
#[should_panic(expected = "bad pts size")]
fn embedder_fill_atom_positions_rejects_size_mismatch() {
    let mol = Molecule::from_smiles("CCO").expect("mol");
    let conf = Conformer3D::new(0, vec![[0.0, 0.0, 0.0]; 3], true);
    let mut pts = vec![ForceFieldVec3::default(); 2];

    embedder_fill_atom_positions(&mut pts, &conf, &mol, &[0]);
}

fn carbon_single_atom_molecule() -> Molecule {
    let mut builder = MoleculeBuilder::new();
    builder.add_atom(AtomSpec::new(Element::C));
    builder.build().expect("single atom")
}

fn two_atom_molecule() -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(BondSpec::new(a0, a1, crate::BondOrder::Single))
        .expect("bond");
    builder.build().expect("two atoms")
}

fn terminal_group_symmetry_molecule() -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let o_minus = builder.add_atom(AtomSpec::new(Element::O).with_formal_charge(-1));
    let carbon = builder.add_atom(AtomSpec::new(Element::C));
    let oxygen = builder.add_atom(AtomSpec::new(Element::O));
    builder
        .add_bond(BondSpec::new(o_minus, carbon, crate::BondOrder::Single))
        .expect("single bond");
    builder
        .add_bond(BondSpec::new(carbon, oxygen, crate::BondOrder::Double))
        .expect("double bond");
    builder.build().expect("terminal-group symmetry molecule")
}

fn random_embed_params(seed: i32, num_threads: i32) -> EmbedParameters {
    EmbedParameters {
        max_iterations: 1,
        random_seed: seed,
        use_random_coords: true,
        num_threads,
        box_size_mult: -2.0,
        use_symmetry_for_pruning: false,
        symmetrize_conjugated_terminal_groups_for_pruning: false,
        ..EmbedParameters::default()
    }
}

#[test]
fn embedder_is_conf_far_from_rest_prunes_close_existing_conformer_and_keeps_far_one() {
    let mut builder = MoleculeBuilder::new();
    builder.add_atom(AtomSpec::new(Element::C));
    builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_conformer(Conformer3D::new(
            0,
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            true,
        ))
        .expect("conformer");
    let mol = builder.build().expect("mol");
    let close = Conformer3D::new(1, vec![[0.01, 0.0, 0.0], [1.01, 0.0, 0.0]], true);
    let far = Conformer3D::new(2, vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]], true);
    let matches = vec![vec![0, 1]];

    assert!(!embedder_is_conf_far_from_rest(&mol, &close, 0.5, &matches));
    assert!(embedder_is_conf_far_from_rest(&mol, &far, 0.5, &matches));
}

#[test]
fn embedder_mol_self_matches_uses_heavy_atoms_or_all_atoms_like_rdkit() {
    let mut builder = MoleculeBuilder::new();
    let c = builder.add_atom(AtomSpec::new(Element::C));
    let h = builder.add_atom(AtomSpec::new(Element::H));
    builder
        .add_bond(BondSpec::new(c, h, crate::BondOrder::Single))
        .expect("bond");
    let mol = builder.build().expect("mol");

    let heavy = EmbedParameters {
        only_heavy_atoms_for_rms: true,
        use_symmetry_for_pruning: false,
        ..EmbedParameters::default()
    };
    assert_eq!(
        embedder_get_mol_self_matches(&mol, &heavy).unwrap(),
        vec![vec![0]]
    );

    let all = EmbedParameters {
        only_heavy_atoms_for_rms: false,
        use_symmetry_for_pruning: false,
        ..EmbedParameters::default()
    };
    assert_eq!(
        embedder_get_mol_self_matches(&mol, &all).unwrap(),
        vec![vec![0, 1]]
    );
}

#[test]
fn symmetrize_terminal_atoms_for_pruning_matches_rdkit_terminal_group_query() {
    let mol = terminal_group_symmetry_molecule();
    let symmetrized = symmetrize_terminal_atoms_for_pruning(&mol).expect("symmetrized");

    assert_eq!(symmetrized.atoms()[0].formal_charge(), 0);
    assert_eq!(symmetrized.atoms()[2].formal_charge(), 0);
    assert_eq!(
        symmetrized.bonds()[0].query(),
        Some(&QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
            BondOrder::Single,
            BondOrder::Double,
        ])))
    );
    assert_eq!(
        symmetrized.bonds()[1].query(),
        Some(&QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
            BondOrder::Single,
            BondOrder::Double,
        ])))
    );

    let matches = get_substruct_matches_with_params(
        &mol,
        &symmetrized,
        &SubstructMatchParams {
            max_matches: 1000,
            uniquify: false,
            use_chirality: false,
            specified_stereo_query_matches_unspecified: false,
        },
    );
    let atom_mappings: Vec<_> = matches.into_iter().map(|m| m.atom_mapping).collect();
    assert_eq!(atom_mappings.len(), 2);
    assert!(atom_mappings.contains(&vec![0, 1, 2]));
    assert!(atom_mappings.contains(&vec![2, 1, 0]));
}

#[test]
fn embedder_mol_self_matches_symmetrize_terminal_groups_for_pruning() {
    let mol = terminal_group_symmetry_molecule();
    let params = EmbedParameters {
        prune_rms_thresh: 0.5,
        use_symmetry_for_pruning: true,
        symmetrize_conjugated_terminal_groups_for_pruning: true,
        ..EmbedParameters::default()
    };

    let matches = embedder_get_mol_self_matches(&mol, &params).expect("self matches");
    assert_eq!(matches.len(), 2);
    assert!(matches.contains(&vec![0, 1, 2]));
    assert!(matches.contains(&vec![2, 1, 0]));
}

#[test]
fn embedder_multiplication_overflows_matches_rdkit_boundaries() {
    assert!(!rdkit_embedder_multiplication_overflows(0, i32::MAX));
    assert!(!rdkit_embedder_multiplication_overflows(1, i32::MAX));
    assert!(!rdkit_embedder_multiplication_overflows(46_340, 46_341));
    assert!(rdkit_embedder_multiplication_overflows(46_341, 46_341));
    assert!(rdkit_embedder_multiplication_overflows(46_342, 46_341));
}

#[test]
fn embedder_helper_uses_thread_scheduling_and_seed_policy() {
    let mut mmat = BoundsMatrix::new(1);
    let mut confs = vec![
        Conformer3D::new(0, vec![[0.0, 0.0, 0.0]], true),
        Conformer3D::new(1, vec![[0.0, 0.0, 0.0]], true),
    ];
    let mut confs_ok = vec![true, true];
    let details = CrystalFFDetails::default();
    let mut params = random_embed_params(11, 2);
    let mut args = EmbedHelperArgs {
        confs_ok: &mut confs_ok,
        four_d: true,
        frag_mapping: None,
        confs: &mut confs,
        frag_idx: 0,
        mmat: &mmat,
        chiral_centers: &[],
        tetrahedral_carbons: &[],
        double_bond_ends: &[],
        stereo_double_bonds: &[],
        etkdg_details: &details,
    };

    embedder_embed_helper(1, 2, &mut args, &mut params, None).expect("embed helper");

    assert_eq!(args.confs[0].coordinates()[0], [0.0, 0.0, 0.0]);
    assert_ne!(args.confs[1].coordinates()[0], [0.0, 0.0, 0.0]);
}

#[test]
fn embed_multiple_confs_generates_value_style_conformers() {
    let mol = carbon_single_atom_molecule();
    let mut params = random_embed_params(7, 1);

    let (embedded, ids) = embed_multiple_confs(&mol, 2, &mut params).expect("embed");

    assert!(mol.conformers_3d().is_empty());
    assert_eq!(ids, vec![0, 1]);
    assert_eq!(embedded.conformers_3d().len(), 2);
}

#[test]
fn embed_multiple_confs_return_vector_returns_ids_without_mutating_input() {
    let mol = carbon_single_atom_molecule();
    let mut params = random_embed_params(7, 1);

    let ids = embed_multiple_confs_return_vector(&mol, 2, &mut params).expect("ids");

    assert_eq!(ids, vec![0, 1]);
    assert!(mol.conformers_3d().is_empty());
}

#[test]
fn embed_molecule_params_returns_first_id_or_minus_one() {
    let mol = carbon_single_atom_molecule();
    let mut params = random_embed_params(3, 1);

    let (embedded, id) = embed_molecule(&mol, &mut params).expect("embed");

    assert_eq!(id, 0);
    assert_eq!(embedded.conformers_3d().len(), 1);
}

#[test]
fn conformer_generation_parameter_parity_custom_bounds_matrix_size_check_matches_rdkit() {
    let mol = Molecule::from_smiles("CC").expect("ethane");
    let mut wrong_size_params = EmbedParameters::etkdg();
    wrong_size_params.random_seed = 42;
    wrong_size_params.num_threads = 1;
    wrong_size_params.bounds_mat = Some(Arc::new(BoundsMatrix::new(mol.num_atoms() + 1)));
    let err = embed_molecule(&mol, &mut wrong_size_params)
        .expect_err("RDKit custom bounds matrix size mismatch must error");
    assert!(
        err.to_string().contains(
            "size of boundsMat provided does not match the number of atoms in the molecule"
        ),
        "unexpected custom bounds error: {err}"
    );
}

#[test]
fn embed_molecule_legacy_overloads_wrapper_constructs_rdkit_parameter_defaults() {
    let mol = carbon_single_atom_molecule();

    let (embedded, id) = rd_distgeom_embed_molecule_wrapper(
        &mol,
        1,
        5,
        true,
        true,
        -2.0,
        true,
        1,
        BTreeMap::new(),
        1e-3,
        false,
        true,
        false,
        false,
        false,
        false,
        false,
        2,
        false,
    )
    .expect("wrapper");

    assert_eq!(id, 0);
    assert_eq!(embedded.conformers_3d().len(), 1);
}

#[test]
fn embed_multiple_confs_legacy_overloads_preserve_seed_policy() {
    assert_eq!(rdkit_embedder_conformer_seed(5, 0, false), 5);
    assert_eq!(rdkit_embedder_conformer_seed(5, 1, false), 10);
    assert_eq!(rdkit_embedder_conformer_seed(5, 1, true), 7);
}

#[test]
fn rd_distgeom_embed_molecule_wrapper_with_params_forwards_parameters() {
    let mol = carbon_single_atom_molecule();
    let mut params = random_embed_params(9, 1);

    let (embedded, id) =
        rd_distgeom_embed_molecule_wrapper_with_params(&mol, &mut params).expect("wrapper");

    assert_eq!(id, 0);
    assert_eq!(embedded.conformers_3d().len(), 1);
}

#[test]
fn rd_distgeom_embed_multiple_confs_wrapper_matches_core_entry_point() {
    let mol = carbon_single_atom_molecule();
    let mut params = random_embed_params(13, 1);

    let (embedded, ids) = embed_multiple_confs(&mol, 3, &mut params).expect("embed");

    assert_eq!(ids, vec![0, 1, 2]);
    assert_eq!(embedded.conformers_3d().len(), 3);
}

#[test]
fn rd_distgeom_parameter_factories_match_embed_parameter_presets() {
    assert_eq!(
        rd_distgeom_get_kdg().to_json(),
        EmbedParameters::kdg().to_json()
    );
    assert_eq!(
        rd_distgeom_get_etdg().to_json(),
        EmbedParameters::etdg().to_json()
    );
    assert_eq!(
        rd_distgeom_get_etdg_v2().to_json(),
        EmbedParameters::etdg_v2().to_json()
    );
    assert_eq!(
        rd_distgeom_get_etkdg().to_json(),
        EmbedParameters::etkdg().to_json()
    );
    assert_eq!(
        rd_distgeom_get_etkdg_v2().to_json(),
        EmbedParameters::etkdg_v2().to_json()
    );
    assert_eq!(
        rd_distgeom_get_etkdg_v3().to_json(),
        EmbedParameters::etkdg_v3().to_json()
    );
    assert_eq!(
        rd_distgeom_get_sr_etkdg_v3().to_json(),
        EmbedParameters::sr_etkdg_v3().to_json()
    );
}

#[test]
fn rd_distgeom_exp_tors_helper_returns_source_backed_details() {
    let mol = two_atom_molecule();
    let details =
        rd_distgeom_get_exp_tors_helper(&mol, false, false, false, false, 2, false).unwrap();

    assert!(details.exp_torsion_angles.is_empty());
    assert!(details.improper_atoms.is_empty());
}

#[test]
fn rd_distgeom_exp_tors_helper_with_params_forwards_parameter_fields() {
    let mol = two_atom_molecule();
    let params = EmbedParameters {
        use_exp_torsion_angle_prefs: false,
        use_small_ring_torsions: true,
        use_macrocycle_torsions: true,
        use_basic_knowledge: false,
        et_version: 2,
        verbose: true,
        ..EmbedParameters::default()
    };

    let from_params = rd_distgeom_get_exp_tors_helper_with_params(&mol, &params).unwrap();
    let direct = rd_distgeom_get_exp_tors_helper(&mol, false, true, true, false, 2, true).unwrap();

    assert_eq!(from_params.exp_torsion_angles, direct.exp_torsion_angles);
    assert_eq!(from_params.improper_atoms, direct.improper_atoms);
}

#[test]
fn rd_distgeom_embed_parameters_json_helper_wraps_source_backed_json() {
    let mut params = EmbedParameters::etkdg_v3();
    params.random_seed = 123;
    params.enable_sequential_random_seeds = true;

    assert_eq!(
        rd_distgeom_embed_parameters_to_json_helper(&params),
        params.to_json()
    );
}

#[test]
fn distgeom_chiral_set_constructor_stores_indices_bounds_and_flags() {
    let chiral_set = ChiralSet::new(
        10,
        11,
        12,
        13,
        14,
        -2.5,
        3.5,
        ChiralSetStructureFlags::InFusedSmallRings as u64,
    );

    assert_eq!(chiral_set.idx0, 10);
    assert_eq!(chiral_set.idx1, 11);
    assert_eq!(chiral_set.idx2, 12);
    assert_eq!(chiral_set.idx3, 13);
    assert_eq!(chiral_set.idx4, 14);
    assert_eq!(chiral_set.volume_lower_bound, -2.5);
    assert_eq!(chiral_set.volume_upper_bound, 3.5);
    assert_eq!(
        chiral_set.structure_flags,
        ChiralSetStructureFlags::InFusedSmallRings as u64
    );
}

#[test]
fn distgeom_chiral_set_default_structure_flags_match_rdkit_default_argument() {
    let chiral_set = ChiralSet::with_default_structure_flags(0, 1, 2, 3, 4, -1.0, 1.0);

    assert_eq!(chiral_set.structure_flags, 0);
}

#[test]
fn distgeom_chiral_set_getters_return_volume_bounds() {
    let chiral_set = ChiralSet::with_default_structure_flags(0, 1, 2, 3, 4, -7.0, -3.0);

    assert_eq!(chiral_set.get_lower_volume_bound(), -7.0);
    assert_eq!(chiral_set.get_upper_volume_bound(), -3.0);
}

#[test]
fn distgeom_chiral_set_allows_equal_volume_bounds() {
    let chiral_set = ChiralSet::with_default_structure_flags(0, 1, 2, 3, 4, 2.0, 2.0);

    assert_eq!(chiral_set.get_lower_volume_bound(), 2.0);
    assert_eq!(chiral_set.get_upper_volume_bound(), 2.0);
}

#[test]
#[should_panic(expected = "Inconsistent bounds")]
fn distgeom_chiral_set_rejects_lower_bound_above_upper_bound() {
    let _ = ChiralSet::with_default_structure_flags(0, 1, 2, 3, 4, 2.0, 1.0);
}

#[test]
fn distgeom_chiral_set_aliases_model_shared_pointer_vector() {
    let chiral_set: ChiralSetPtr = Arc::new(ChiralSet::with_default_structure_flags(
        0, 1, 2, 3, 4, -1.0, 1.0,
    ));
    let chiral_sets: VectChiralSet = vec![Arc::clone(&chiral_set)];

    assert_eq!(chiral_sets.len(), 1);
    assert!(Arc::ptr_eq(&chiral_set, &chiral_sets[0]));
}

#[test]
fn chiral_volume_flat_returns_signed_triple_product() {
    let pos = [
        1.0, 0.0, 0.0, //
        0.0, 1.0, 0.0, //
        0.0, 0.0, 1.0, //
        0.0, 0.0, 0.0,
    ];

    assert_eq!(calc_chiral_volume_flat(0, 1, 2, 3, &pos, 3), 1.0);
    assert_eq!(calc_chiral_volume_flat(0, 2, 1, 3, &pos, 3), -1.0);
}

#[test]
fn chiral_volume_flat_uses_idx4_as_reference_point() {
    let pos = [
        2.0, 2.0, 2.0, //
        3.0, 2.0, 2.0, //
        2.0, 3.0, 2.0, //
        2.0, 2.0, 3.0,
    ];

    assert_eq!(calc_chiral_volume_flat(1, 2, 3, 0, &pos, 3), 1.0);
}

#[test]
fn chiral_volume_flat_ignores_dimensions_after_first_three() {
    let pos = [
        1.0, 0.0, 0.0, 100.0, //
        0.0, 1.0, 0.0, 200.0, //
        0.0, 0.0, 1.0, 300.0, //
        0.0, 0.0, 0.0, 400.0,
    ];

    assert_eq!(calc_chiral_volume_flat(0, 1, 2, 3, &pos, 4), 1.0);
}

#[test]
fn chiral_volume_points_returns_signed_triple_product() {
    let pts = [
        ForceFieldVec3::new(1.0, 0.0, 0.0),
        ForceFieldVec3::new(0.0, 1.0, 0.0),
        ForceFieldVec3::new(0.0, 0.0, 1.0),
        ForceFieldVec3::new(0.0, 0.0, 0.0),
    ];

    assert_eq!(calc_chiral_volume_points(0, 1, 2, 3, &pts), 1.0);
    assert_eq!(calc_chiral_volume_points(0, 2, 1, 3, &pts), -1.0);
}

fn chiral_violation_forcefield() -> ForceField {
    let mut ff = ForceField::new(3);
    ff.positions_mut().push(ForceFieldVec3::new(1.0, 0.0, 0.0));
    ff.positions_mut().push(ForceFieldVec3::new(0.0, 1.0, 0.0));
    ff.positions_mut().push(ForceFieldVec3::new(0.0, 0.0, 1.0));
    ff.positions_mut().push(ForceFieldVec3::new(0.0, 0.0, 0.0));
    ff
}

fn chiral_violation_pos() -> Vec<f64> {
    vec![
        1.0, 0.0, 0.0, //
        0.0, 1.0, 0.0, //
        0.0, 0.0, 1.0, //
        0.0, 0.0, 0.0,
    ]
}

#[test]
fn chiral_violation_contribs_constructor_starts_empty() {
    let ff = chiral_violation_forcefield();
    let contribs = ChiralViolationContribs::new(&ff);

    assert!(contribs.empty());
    assert_eq!(contribs.size(), 0);
}

#[test]
fn chiral_violation_contribs_add_contrib_copies_chiral_set_bounds_and_weight() {
    let ff = chiral_violation_forcefield();
    let mut contribs = ChiralViolationContribs::new(&ff);
    let cset = ChiralSet::with_default_structure_flags(99, 0, 1, 2, 3, -0.5, 0.5);

    contribs.add_contrib(&cset, 2.5);

    assert!(!contribs.empty());
    assert_eq!(contribs.size(), 1);
    assert_eq!(
        contribs.contribs()[0],
        ChiralViolationContribsParams::new(0, 1, 2, 3, 0.5, -0.5, 2.5)
    );
}

#[test]
#[should_panic]
fn chiral_violation_contribs_add_contrib_rejects_out_of_range_indices() {
    let ff = chiral_violation_forcefield();
    let mut contribs = ChiralViolationContribs::new(&ff);
    let cset = ChiralSet::with_default_structure_flags(0, 0, 1, 2, 4, -0.5, 0.5);

    contribs.add_contrib(&cset, 1.0);
}

#[test]
fn chiral_violation_contribs_get_energy_returns_zero_inside_bounds() {
    let ff = chiral_violation_forcefield();
    let mut contribs = ChiralViolationContribs::new(&ff);
    let cset = ChiralSet::with_default_structure_flags(0, 0, 1, 2, 3, 0.5, 1.5);
    contribs.add_contrib(&cset, 2.0);

    assert_eq!(contribs.get_energy(&chiral_violation_pos()), 0.0);
}

#[test]
fn chiral_violation_contribs_get_energy_accumulates_lower_and_upper_violations() {
    let ff = chiral_violation_forcefield();
    let mut contribs = ChiralViolationContribs::new(&ff);
    let upper = ChiralSet::with_default_structure_flags(0, 0, 1, 2, 3, -0.5, 0.5);
    let lower = ChiralSet::with_default_structure_flags(0, 0, 2, 1, 3, -0.5, 0.5);
    contribs.add_contrib(&upper, 2.0);
    contribs.add_contrib(&lower, 3.0);

    assert_eq!(
        contribs.get_energy(&chiral_violation_pos()),
        2.0 * 0.5 * 0.5 + 3.0 * 0.5 * 0.5
    );
}

#[test]
fn chiral_violation_contribs_get_grad_returns_early_inside_bounds() {
    let ff = chiral_violation_forcefield();
    let mut contribs = ChiralViolationContribs::new(&ff);
    let cset = ChiralSet::with_default_structure_flags(0, 0, 1, 2, 3, 0.5, 1.5);
    contribs.add_contrib(&cset, 2.0);
    let mut grad = vec![10.0; 12];

    contribs.get_grad(&chiral_violation_pos(), &mut grad);

    assert_eq!(grad, vec![10.0; 12]);
}

#[test]
fn chiral_violation_contribs_get_grad_matches_source_formula_for_upper_violation() {
    let ff = chiral_violation_forcefield();
    let mut contribs = ChiralViolationContribs::new(&ff);
    let cset = ChiralSet::with_default_structure_flags(0, 0, 1, 2, 3, -1.0, 0.0);
    contribs.add_contrib(&cset, 2.0);
    let mut grad = vec![0.0; 12];

    contribs.get_grad(&chiral_violation_pos(), &mut grad);

    assert_eq!(
        grad,
        vec![
            2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0, -2.0, -2.0, -2.0
        ]
    );
}

fn dist_violation_forcefield() -> ForceField {
    let mut ff = ForceField::new(3);
    ff.positions_mut().push(ForceFieldVec3::new(0.0, 0.0, 0.0));
    ff.positions_mut().push(ForceFieldVec3::new(2.0, 0.0, 0.0));
    ff
}

fn dist_violation_pos(distance: f64) -> Vec<f64> {
    vec![0.0, 0.0, 0.0, distance, 0.0, 0.0]
}

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1.0e-12,
        "actual={actual} expected={expected}"
    );
}

#[test]
fn dist_violation_contribs_constructor_and_add_contrib_store_squared_bounds() {
    let ff = dist_violation_forcefield();
    let mut contribs = DistViolationContribs::new(&ff);
    assert!(contribs.empty());

    contribs.add_contrib(0, 1, 3.0, 1.5, 2.0);

    assert_eq!(contribs.size(), 1);
    assert_eq!(
        contribs.contribs()[0],
        DistViolationContribsParams::new(0, 1, 3.0, 1.5, 2.0)
    );
    assert_eq!(contribs.contribs()[0].ub2, 9.0);
    assert_eq!(contribs.contribs()[0].lb2, 2.25);
}

#[test]
fn dist_violation_contribs_distance_helpers_follow_source_dim_loop() {
    let pos = [0.0, 0.0, 0.0, 5.0, 1.0, 2.0, 2.0, 9.0];

    assert_eq!(dist_violation_distance2(0, 1, &pos, 4), 25.0);
    assert_eq!(dist_violation_distance(0, 1, &pos, 4), 25.0_f64.sqrt());
}

#[test]
fn dist_violation_contribs_get_energy_returns_zero_inside_bounds() {
    let ff = dist_violation_forcefield();
    let mut contribs = DistViolationContribs::new(&ff);
    contribs.add_contrib(0, 1, 3.0, 1.0, 2.0);

    assert_eq!(contribs.get_energy(&dist_violation_pos(2.0)), 0.0);
}

#[test]
fn dist_violation_contribs_get_energy_matches_upper_violation_formula() {
    let ff = dist_violation_forcefield();
    let mut contribs = DistViolationContribs::new(&ff);
    contribs.add_contrib(0, 1, 1.0, 0.0, 2.0);

    assert_eq!(contribs.get_energy(&dist_violation_pos(2.0)), 18.0);
}

#[test]
fn dist_violation_contribs_get_energy_matches_lower_violation_formula() {
    let ff = dist_violation_forcefield();
    let mut contribs = DistViolationContribs::new(&ff);
    contribs.add_contrib(0, 1, 10.0, 2.0, 3.0);

    assert_close(contribs.get_energy(&dist_violation_pos(1.0)), 1.08);
}

#[test]
fn dist_violation_contribs_get_grad_returns_early_inside_bounds() {
    let ff = dist_violation_forcefield();
    let mut contribs = DistViolationContribs::new(&ff);
    contribs.add_contrib(0, 1, 3.0, 1.0, 2.0);
    let mut grad = vec![7.0; 6];

    contribs.get_grad(&dist_violation_pos(2.0), &mut grad);

    assert_eq!(grad, vec![7.0; 6]);
}

#[test]
fn dist_violation_contribs_get_grad_matches_upper_violation_formula() {
    let ff = dist_violation_forcefield();
    let mut contribs = DistViolationContribs::new(&ff);
    contribs.add_contrib(0, 1, 1.0, 0.0, 2.0);
    let mut grad = vec![0.0; 6];

    contribs.get_grad(&dist_violation_pos(2.0), &mut grad);

    assert_eq!(grad, vec![-48.0, 0.0, 0.0, 48.0, 0.0, 0.0]);
}

#[test]
fn dist_violation_contribs_get_grad_matches_lower_violation_formula() {
    let ff = dist_violation_forcefield();
    let mut contribs = DistViolationContribs::new(&ff);
    contribs.add_contrib(0, 1, 10.0, 2.0, 3.0);
    let mut grad = vec![0.0; 6];

    contribs.get_grad(&dist_violation_pos(1.0), &mut grad);

    assert_close(grad[0], 2.304);
    assert_eq!(grad[1], 0.0);
    assert_eq!(grad[2], 0.0);
    assert_close(grad[3], -2.304);
    assert_eq!(grad[4], 0.0);
    assert_eq!(grad[5], 0.0);
}

fn fourth_dim_forcefield() -> ForceField {
    let mut ff = ForceField::new(4);
    ff.positions_mut().push(ForceFieldVec3::new(0.0, 0.0, 0.0));
    ff.positions_mut().push(ForceFieldVec3::new(1.0, 0.0, 0.0));
    ff
}

fn fourth_dim_pos() -> Vec<f64> {
    vec![
        0.0, 0.0, 0.0, 2.0, //
        1.0, 2.0, 3.0, -3.0,
    ]
}

#[test]
fn fourth_dim_contribs_default_starts_without_owner_and_empty() {
    let contribs = FourthDimContribs::default();

    assert!(contribs.empty());
    assert_eq!(contribs.size(), 0);
}

#[test]
fn fourth_dim_contribs_constructor_requires_four_dimensional_forcefield() {
    let ff = fourth_dim_forcefield();
    let contribs = FourthDimContribs::new(&ff);

    assert!(contribs.empty());
    assert_eq!(contribs.size(), 0);
}

#[test]
#[should_panic(expected = "force field has wrong dimension")]
fn fourth_dim_contribs_constructor_rejects_non_four_dimensional_forcefield() {
    let ff = ForceField::new(3);

    let _ = FourthDimContribs::new(&ff);
}

#[test]
fn fourth_dim_contribs_add_contrib_appends_index_and_weight() {
    let ff = fourth_dim_forcefield();
    let mut contribs = FourthDimContribs::new(&ff);

    contribs.add_contrib(1, 2.5);

    assert!(!contribs.empty());
    assert_eq!(contribs.size(), 1);
    assert_eq!(contribs.contribs()[0], FourthDimContribsParams::new(1, 2.5));
}

#[test]
fn fourth_dim_contribs_get_energy_accumulates_weighted_fourth_coordinate_squares() {
    let ff = fourth_dim_forcefield();
    let mut contribs = FourthDimContribs::new(&ff);
    contribs.add_contrib(0, 2.0);
    contribs.add_contrib(1, 3.0);

    assert_eq!(contribs.get_energy(&fourth_dim_pos()), 35.0);
}

#[test]
fn fourth_dim_contribs_get_grad_adds_source_weighted_fourth_coordinate_terms() {
    let ff = fourth_dim_forcefield();
    let mut contribs = FourthDimContribs::new(&ff);
    contribs.add_contrib(0, 2.0);
    contribs.add_contrib(1, 3.0);
    let mut grad = vec![10.0; 8];

    contribs.get_grad(&fourth_dim_pos(), &mut grad);

    assert_eq!(grad, vec![10.0, 10.0, 10.0, 14.0, 10.0, 10.0, 10.0, 1.0]);
}

#[test]
fn fourth_dim_contribs_copy_preserves_contribs_and_behavior() {
    let ff = fourth_dim_forcefield();
    let mut contribs = FourthDimContribs::new(&ff);
    contribs.add_contrib(1, 3.0);
    let copied = contribs.copy();

    assert_eq!(copied.get_energy(&fourth_dim_pos()), 27.0);
}

#[derive(Debug)]
struct FixedDoubleRng {
    values: Vec<f64>,
    idx: usize,
}

impl FixedDoubleRng {
    fn new(values: Vec<f64>) -> Self {
        Self { values, idx: 0 }
    }
}

impl RdkitDoubleRng for FixedDoubleRng {
    fn next_unit_f64(&mut self) -> f64 {
        let value = self.values[self.idx];
        self.idx += 1;
        value
    }
}

fn rdkit_minstd_next_raw(state: &mut u64) -> u32 {
    *state = (*state * RDKIT_RANDOM_MULTIPLIER) % RDKIT_RANDOM_MODULUS;
    *state as u32
}

fn rdkit_minstd_next_unit(state: &mut u64) -> f64 {
    let raw = rdkit_minstd_next_raw(state) as f64;
    (raw - 1.0) / (RDKIT_RANDOM_MODULUS as f64 - 1.0)
}

fn rdkit_minstd_seed_state(seed: i32) -> u64 {
    let mut state = seed as u64 % RDKIT_RANDOM_MODULUS;
    if state == 0 {
        state = 1;
    }
    state
}

fn pick_random_dist_mat_bounds() -> BoundsMatrix {
    let mut mmat = BoundsMatrix::new(3);
    mmat.set_lower(1, 0, 1.0).expect("set lower");
    mmat.set_upper(1, 0, 3.0).expect("set upper");
    mmat.set_lower(2, 0, 2.0).expect("set lower");
    mmat.set_upper(2, 0, 6.0).expect("set upper");
    mmat.set_lower(2, 1, 4.0).expect("set lower");
    mmat.set_upper(2, 1, 8.0).expect("set upper");
    mmat
}

#[test]
fn pick_random_dist_mat_rng_overload_fills_lower_triangle_in_source_order() {
    let mmat = pick_random_dist_mat_bounds();
    let mut dist_mat = SymmMatrix::new(3);
    let mut rng = FixedDoubleRng::new(vec![0.25, 0.5, 0.75]);

    let largest = pick_random_dist_mat_with_rng(&mmat, &mut dist_mat, &mut rng);

    assert_eq!(dist_mat.get_data(), &[0.0, 1.5, 0.0, 4.0, 7.0, 0.0]);
    assert_eq!(dist_mat.get_val(1, 0), 1.5);
    assert_eq!(dist_mat.get_val(0, 1), 1.5);
    assert_eq!(largest, 7.0);
}

#[test]
fn pick_random_dist_mat_seed_overload_reseeds_rdkit_minstd_rng() {
    let mmat = pick_random_dist_mat_bounds();
    let mut first = SymmMatrix::new(3);
    let mut second = SymmMatrix::new(3);

    let first_largest = pick_random_dist_mat(&mmat, &mut first, 1);
    let second_largest = pick_random_dist_mat(&mmat, &mut second, 1);

    assert_eq!(first.get_data(), second.get_data());
    assert_eq!(first_largest, second_largest);

    let raw = 48_271.0;
    let rval = (raw - 1.0) / (2_147_483_647.0 - 1.0);
    assert_eq!(first.get_val(1, 0), 1.0 + rval * 2.0);
}

#[test]
fn rdkit_minstd_rand_matches_boost_seed_modulus_boundary() {
    let mut direct = RdkitDistgeomMinStdRand::new(1);
    let mut wrapped = RdkitDistgeomMinStdRand::new(2_147_483_647);

    assert_eq!(direct.next_unit_f64(), wrapped.next_unit_f64());
}

#[test]
fn embedder_conformer_seed_policy_matches_rdkit_non_sequential_and_sequential_modes() {
    assert_eq!(rdkit_embedder_conformer_seed(-1, 0, false), -1);
    assert_eq!(rdkit_embedder_conformer_seed(0, 0, false), 0);
    assert_eq!(rdkit_embedder_conformer_seed(7, 0, false), 7);
    assert_eq!(rdkit_embedder_conformer_seed(7, 1, false), 14);
    assert_eq!(rdkit_embedder_conformer_seed(7, 2, false), 21);

    assert_eq!(rdkit_embedder_conformer_seed(0, 0, true), 1);
    assert_eq!(rdkit_embedder_conformer_seed(7, 0, true), 8);
    assert_eq!(rdkit_embedder_conformer_seed(7, 2, true), 10);
}

#[test]
fn embedder_conformer_seed_policy_matches_rdkit_overflow_hash_branch() {
    assert!(rdkit_embedder_multiplication_overflows(46_342, 46_341));
    assert_eq!(
        rdkit_embedder_conformer_seed(46_341, 46_341, false),
        143_656
    );
    assert_eq!(
        rdkit_embedder_conformer_seed(100_000, 100_000, false),
        1_410_365_413
    );
}

#[test]
fn vector_set_to_random_clock_seeded_path_returns_normalized_vector() {
    let vec = rdkit_vector_set_to_random(3, 0).expect("clock-seeded vector");
    let norm_sq = vec.iter().map(|value| value * value).sum::<f64>();

    assert_eq!(vec.len(), 3);
    assert!(vec.iter().all(|value| value.is_finite()));
    assert!((norm_sq - 1.0).abs() < 1.0e-9);
}

#[test]
fn power_eigen_solver_accepts_clock_seeded_path() {
    let mut mat = symm_matrix_from_distances(1, &[(0, 0, 1.0)]);
    let mut eigen_values = vec![0.0; 1];

    assert!(power_eigen_solver(1, &mut mat, &mut eigen_values, None, 0).expect("power solver"));
    assert!(eigen_values[0].is_finite());
}

#[test]
fn pick_random_dist_mat_seed_overload_preserves_unseeded_global_stream() {
    let mmat = pick_random_dist_mat_bounds();
    let mut seeded = SymmMatrix::new(3);
    let mut continued = SymmMatrix::new(3);

    pick_random_dist_mat(&mmat, &mut seeded, 7);
    pick_random_dist_mat(&mmat, &mut continued, -1);

    assert_ne!(seeded.get_data(), continued.get_data());
}

#[test]
#[should_panic(expected = "Size mismatch")]
fn pick_random_dist_mat_rejects_size_mismatch() {
    let mmat = pick_random_dist_mat_bounds();
    let mut dist_mat = SymmMatrix::new(2);
    let mut rng = FixedDoubleRng::new(vec![0.0]);

    pick_random_dist_mat_with_rng(&mmat, &mut dist_mat, &mut rng);
}

#[test]
#[should_panic]
fn pick_random_dist_mat_rejects_upper_bound_below_lower_bound() {
    let mut mmat = BoundsMatrix::new(2);
    mmat.set_lower(1, 0, 3.0).expect("set lower");
    mmat.set_upper(1, 0, 2.0).expect("set upper");
    let mut dist_mat = SymmMatrix::new(2);
    let mut rng = FixedDoubleRng::new(vec![0.0]);

    pick_random_dist_mat_with_rng(&mmat, &mut dist_mat, &mut rng);
}

#[test]
fn symm_matrix_set_val_uses_same_lower_triangular_storage_for_either_index_order() {
    let mut mat = SymmMatrix::new(3);

    mat.set_val(0, 2, 4.5);

    assert_eq!(mat.get_val(2, 0), 4.5);
    assert_eq!(mat.get_data()[3], 4.5);
}

fn symm_matrix_from_distances(n: usize, distances: &[(usize, usize, f64)]) -> SymmMatrix {
    let mut mat = SymmMatrix::new(n);
    for &(i, j, value) in distances {
        mat.set_val(i, j, value);
    }
    mat
}

fn point_distance(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(ai, bi)| {
            let d = ai - bi;
            d * d
        })
        .sum::<f64>()
        .sqrt()
}

#[test]
fn compute_initial_coords_embeds_two_point_distance() {
    let dist_mat = symm_matrix_from_distances(2, &[(1, 0, 2.0)]);
    let mut positions = vec![vec![0.0; 3], vec![0.0; 3]];

    assert!(
        compute_initial_coords(&dist_mat, &mut positions, false, 2, 11)
            .expect("compute initial coords")
    );

    assert_close(point_distance(&positions[0], &positions[1]), 2.0);
}

#[test]
fn compute_initial_coords_rng_overload_rejects_size_mismatch() {
    let dist_mat = symm_matrix_from_distances(2, &[(1, 0, 2.0)]);
    let mut positions = vec![vec![0.0; 3]];
    let mut rng = FixedDoubleRng::new(vec![0.25]);

    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        compute_initial_coords_with_rng(&dist_mat, &mut positions, &mut rng, false, 2)
            .expect("compute initial coords");
    }));

    assert!(result.is_err());
}

#[test]
fn compute_initial_coords_fails_when_zero_eigen_threshold_is_reached_for_more_than_three_points() {
    let dist_mat = SymmMatrix::new(4);
    let mut positions = vec![vec![0.0; 3]; 4];

    assert!(
        !compute_initial_coords(&dist_mat, &mut positions, false, 2, 17)
            .expect("compute initial coords")
    );
}

#[test]
fn compute_initial_coords_fails_negative_eigenvalue_without_randomization() {
    let dist_mat = symm_matrix_from_distances(3, &[(1, 0, 1.0), (2, 0, 1.0), (2, 1, 3.0)]);
    let mut positions = vec![vec![0.0; 3]; 3];

    assert!(
        !compute_initial_coords(&dist_mat, &mut positions, false, 2, 19)
            .expect("compute initial coords")
    );
}

#[test]
fn compute_initial_coords_randomizes_negative_eigenvalue_when_requested() {
    let dist_mat = symm_matrix_from_distances(3, &[(1, 0, 1.0), (2, 0, 1.0), (2, 1, 3.0)]);
    let mut positions = vec![vec![0.0; 3]; 3];
    let mut rng = FixedDoubleRng::new(vec![0.25, 0.75, 0.5, 0.125, 0.875, 0.625]);

    assert!(
        compute_initial_coords_with_rng(&dist_mat, &mut positions, &mut rng, true, 2)
            .expect("compute initial coords")
    );

    assert!(positions.iter().flatten().all(|coord| coord.is_finite()));
}

#[test]
fn compute_random_coords_rng_overload_fills_points_in_source_order_exactly() {
    let mut positions = vec![vec![0.0; 3], vec![0.0; 2], vec![0.0; 1]];
    let mut rng = FixedDoubleRng::new(vec![0.0, 0.25, 0.5, 0.75, 1.0, 0.125]);

    assert!(compute_random_coords_with_rng(
        &mut positions,
        4.0,
        &mut rng
    ));

    assert_eq!(
        positions,
        vec![vec![-2.0, -1.0, 0.0], vec![1.0, 2.0], vec![-1.5]]
    );
    assert_eq!(rng.idx, 6);
}

#[test]
fn compute_random_coords_seed_overload_reseeds_boost_minstd_rng_exactly() {
    let mut positions = vec![vec![0.0; 3], vec![0.0; 3]];
    let mut expected_state = rdkit_minstd_seed_state(42);
    let expected = (0..2)
        .map(|_| {
            (0..3)
                .map(|_| 6.0 * (rdkit_minstd_next_unit(&mut expected_state) - 0.5))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    assert!(compute_random_coords(&mut positions, 6.0, 42));

    assert_eq!(positions, expected);
}

#[test]
fn compute_random_coords_seed_boundary_matches_boost_seed_zero_adjustment() {
    let mut wrapped_positions = vec![vec![0.0; 2]];
    let mut direct_positions = vec![vec![0.0; 2]];

    assert!(compute_random_coords(
        &mut wrapped_positions,
        2.0,
        2_147_483_647
    ));
    assert!(compute_random_coords(&mut direct_positions, 2.0, 1));

    assert_eq!(wrapped_positions, direct_positions);
}

#[test]
fn compute_random_coords_empty_positions_returns_true_without_consuming_rng() {
    let mut positions: Vec<Vec<f64>> = Vec::new();
    let mut rng = FixedDoubleRng::new(vec![0.25]);

    assert!(compute_random_coords_with_rng(
        &mut positions,
        1.0,
        &mut rng
    ));

    assert!(positions.is_empty());
    assert_eq!(rng.idx, 0);
}

#[test]
#[should_panic(expected = "bad boxSize")]
fn compute_random_coords_rejects_non_positive_box_size() {
    let mut positions = vec![vec![0.0; 3]];
    let mut rng = FixedDoubleRng::new(vec![0.25]);

    compute_random_coords_with_rng(&mut positions, 0.0, &mut rng);
}

#[test]
fn construct_distgeom_forcefield_adds_distance_terms_for_basin_and_extra_weights() {
    let mut mmat = BoundsMatrix::new(3);
    mmat.set_lower(1, 0, 1.0).expect("set lower");
    mmat.set_upper(1, 0, 2.0).expect("set upper");
    mmat.set_lower(2, 0, 1.0).expect("set lower");
    mmat.set_upper(2, 0, 2.0).expect("set upper");
    mmat.set_lower(2, 1, 1.0).expect("set lower");
    mmat.set_upper(2, 1, 4.0).expect("set upper");
    let positions = vec![
        vec![0.0, 0.0, 0.0],
        vec![3.0, 0.0, 0.0],
        vec![0.0, 4.0, 0.0],
    ];
    let mut extra_weights = std::collections::BTreeMap::new();
    extra_weights.insert((2, 1), 2.0);
    let fixed_pts = vec![true, false, true];

    let mut field = construct_distgeom_forcefield(
        &mmat,
        &positions,
        &[],
        0.0,
        0.0,
        Some(&extra_weights),
        1.0,
        Some(&fixed_pts),
    );
    field.initialize();
    let mut contrib_energies = Vec::new();
    let energy = field.calc_energy_current(Some(&mut contrib_energies));

    assert_eq!(field.dimension(), 3);
    assert_eq!(
        field.positions(),
        &positions
            .iter()
            .map(|p| ForceFieldVec3::new(p[0], p[1], p[2]))
            .collect::<Vec<_>>()
    );
    assert_eq!(contrib_energies, vec![1.5625 + 0.6328125]);
    assert_eq!(energy, 2.1953125);
}

#[test]
fn construct_distgeom_forcefield_adds_chiral_and_fourth_dimension_terms() {
    let mmat = BoundsMatrix::new(4);
    let positions = vec![
        vec![1.0, 0.0, 0.0, 1.0],
        vec![0.0, 1.0, 0.0, 2.0],
        vec![0.0, 0.0, 1.0, 3.0],
        vec![0.0, 0.0, 0.0, 4.0],
    ];
    let cset = Arc::new(ChiralSet::with_default_structure_flags(
        99, 0, 1, 2, 3, -0.1, 0.1,
    ));

    let mut field =
        construct_distgeom_forcefield(&mmat, &positions, &[cset], 2.0, 0.5, None, 0.0, None);
    field.initialize();
    let mut contrib_energies = Vec::new();
    let energy = field.calc_energy_current(Some(&mut contrib_energies));
    let expected_chiral = 2.0 * (1.0_f64 - 0.1) * (1.0_f64 - 0.1);
    let expected_fourth = 0.5 * (1.0_f64 + 4.0 + 9.0 + 16.0);

    assert_eq!(field.dimension(), 4);
    assert_eq!(
        field.positions(),
        &[
            ForceFieldVec3::new4(1.0, 0.0, 0.0, 1.0),
            ForceFieldVec3::new4(0.0, 1.0, 0.0, 2.0),
            ForceFieldVec3::new4(0.0, 0.0, 1.0, 3.0),
            ForceFieldVec3::new4(0.0, 0.0, 0.0, 4.0),
        ]
    );
    assert_eq!(contrib_energies, vec![expected_chiral, expected_fourth]);
    assert_eq!(energy, expected_chiral + expected_fourth);
}

#[test]
fn construct_distgeom_forcefield_omits_empty_contrib_groups() {
    let mmat = BoundsMatrix::new(2);
    let positions = vec![vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0]];

    let mut field =
        construct_distgeom_forcefield(&mmat, &positions, &[], 0.0, 0.0, None, 0.0, None);
    field.initialize();
    let mut contrib_energies = Vec::new();

    assert_eq!(field.calc_energy_current(Some(&mut contrib_energies)), 0.0);
    assert!(contrib_energies.is_empty());
}

#[test]
#[should_panic]
fn construct_distgeom_forcefield_rejects_size_mismatch() {
    let mmat = BoundsMatrix::new(2);
    let positions = vec![vec![0.0, 0.0, 0.0]];

    let _ = construct_distgeom_forcefield(&mmat, &positions, &[], 0.0, 0.0, None, 0.0, None);
}

fn distgeom_improper_force_field() -> ForceField {
    let mut ff = ForceField::new(3);
    ff.positions_mut().push(ForceFieldVec3::new(0.1, -0.2, 0.3));
    ff.positions_mut().push(ForceFieldVec3::new(1.0, 0.0, 0.1));
    ff.positions_mut().push(ForceFieldVec3::new(-0.2, 1.1, 0.2));
    ff.positions_mut().push(ForceFieldVec3::new(0.2, -0.1, 1.3));
    ff
}

fn flatten_force_field_positions(ff: &ForceField) -> Vec<f64> {
    ff.positions()
        .iter()
        .flat_map(|point| [point.x, point.y, point.z])
        .collect()
}

#[test]
fn distgeom_improper_torsion_terms_add_three_source_permutations_and_mark_center() {
    let mut ff = distgeom_improper_force_field();
    let improper_atoms = vec![vec![0, 1, 2, 3, 6, 1]];
    let mut is_improper_constrained = vec![false; 4];

    add_improper_torsion_terms(&mut ff, 2.0, &improper_atoms, &mut is_improper_constrained);
    ff.initialize();
    let mut contrib_energies = Vec::new();
    let energy = ff.calc_energy_current(Some(&mut contrib_energies));

    let expected_ff = distgeom_improper_force_field();
    let mut expected =
        crate::chemistry::forcefield::uff::inversion::InversionContribs::new(&expected_ff);
    expected.add_contrib(0, 1, 2, 3, 6, true, 2.0);
    expected.add_contrib(0, 1, 3, 2, 6, true, 2.0);
    expected.add_contrib(2, 1, 3, 0, 6, true, 2.0);
    let expected_energy = expected.get_energy(&flatten_force_field_positions(&expected_ff));

    assert_eq!(is_improper_constrained, vec![false, true, false, false]);
    assert_eq!(contrib_energies.len(), 1);
    assert_close(contrib_energies[0], expected_energy);
    assert_close(energy, expected_energy);
}

#[test]
fn distgeom_improper_torsion_terms_accumulate_multiple_centers_in_one_contrib_group() {
    let mut ff = distgeom_improper_force_field();
    let improper_atoms = vec![vec![0, 1, 2, 3, 6, 0], vec![3, 2, 1, 0, 15, 0]];
    let mut is_improper_constrained = vec![false; 4];

    add_improper_torsion_terms(&mut ff, 0.5, &improper_atoms, &mut is_improper_constrained);
    ff.initialize();
    let mut contrib_energies = Vec::new();
    let energy = ff.calc_energy_current(Some(&mut contrib_energies));

    let expected_ff = distgeom_improper_force_field();
    let mut expected =
        crate::chemistry::forcefield::uff::inversion::InversionContribs::new(&expected_ff);
    expected.add_contrib(0, 1, 2, 3, 6, false, 0.5);
    expected.add_contrib(0, 1, 3, 2, 6, false, 0.5);
    expected.add_contrib(2, 1, 3, 0, 6, false, 0.5);
    expected.add_contrib(3, 2, 1, 0, 15, false, 0.5);
    expected.add_contrib(3, 2, 0, 1, 15, false, 0.5);
    expected.add_contrib(1, 2, 0, 3, 15, false, 0.5);
    let expected_energy = expected.get_energy(&flatten_force_field_positions(&expected_ff));

    assert_eq!(is_improper_constrained, vec![false, true, true, false]);
    assert_eq!(contrib_energies.len(), 1);
    assert_close(contrib_energies[0], expected_energy);
    assert_close(energy, expected_energy);
}

#[test]
fn distgeom_improper_torsion_terms_leave_force_field_unchanged_for_empty_input() {
    let mut ff = distgeom_improper_force_field();
    let mut is_improper_constrained = vec![false; 4];

    add_improper_torsion_terms(&mut ff, 1.0, &[], &mut is_improper_constrained);
    ff.initialize();
    let mut contrib_energies = Vec::new();

    assert_eq!(ff.calc_energy_current(Some(&mut contrib_energies)), 0.0);
    assert!(contrib_energies.is_empty());
    assert_eq!(is_improper_constrained, vec![false; 4]);
}

#[test]
#[should_panic]
fn distgeom_improper_torsion_terms_reject_short_improper_atom_record() {
    let mut ff = distgeom_improper_force_field();
    let improper_atoms = vec![vec![0, 1, 2, 3, 6]];
    let mut is_improper_constrained = vec![false; 4];

    add_improper_torsion_terms(&mut ff, 1.0, &improper_atoms, &mut is_improper_constrained);
}

fn distgeom_experimental_torsion_force_field() -> ForceField {
    let mut ff = ForceField::new(3);
    ff.positions_mut().push(ForceFieldVec3::new(0.0, 0.0, 0.0));
    ff.positions_mut().push(ForceFieldVec3::new(1.0, 0.0, 0.0));
    ff.positions_mut().push(ForceFieldVec3::new(1.0, 1.0, 0.0));
    ff.positions_mut().push(ForceFieldVec3::new(1.0, 1.0, 1.0));
    ff.positions_mut().push(ForceFieldVec3::new(2.0, 1.0, 1.0));
    ff
}

#[test]
fn distgeom_experimental_torsion_terms_mark_ordered_endpoint_pairs_and_add_contribs() {
    let mut ff = distgeom_experimental_torsion_force_field();
    let details = CrystalFFDetails {
        exp_torsion_atoms: vec![vec![0, 1, 2, 3], vec![4, 3, 2, 1]],
        exp_torsion_angles: vec![
            (
                vec![1, -1, 1, -1, 1, -1],
                vec![0.10, 0.20, 0.30, 0.40, 0.50, 0.60],
            ),
            (
                vec![-1, 1, -1, 1, -1, 1],
                vec![0.05, 0.15, 0.25, 0.35, 0.45, 0.55],
            ),
        ],
        ..CrystalFFDetails::default()
    };
    let mut atom_pairs = vec![false; 25];

    add_experimental_torsion_terms(&mut ff, &details, &mut atom_pairs, 5);
    ff.initialize();
    let mut contrib_energies = Vec::new();
    let energy = ff.calc_energy_current(Some(&mut contrib_energies));

    let expected_ff = distgeom_experimental_torsion_force_field();
    let mut expected =
        crate::chemistry::forcefield::crystalff::TorsionAngleContribs::new(&expected_ff);
    expected.add_contrib(
        0,
        1,
        2,
        3,
        vec![0.10, 0.20, 0.30, 0.40, 0.50, 0.60],
        vec![1, -1, 1, -1, 1, -1],
    );
    expected.add_contrib(
        4,
        3,
        2,
        1,
        vec![0.05, 0.15, 0.25, 0.35, 0.45, 0.55],
        vec![-1, 1, -1, 1, -1, 1],
    );
    let expected_energy = expected.get_energy(&flatten_force_field_positions(&expected_ff));

    assert!(atom_pairs[3]);
    assert!(atom_pairs[9]);
    assert_eq!(atom_pairs.iter().filter(|&&set| set).count(), 2);
    assert_eq!(contrib_energies.len(), 1);
    assert_close(contrib_energies[0], expected_energy);
    assert_close(energy, expected_energy);
}

#[test]
fn distgeom_experimental_torsion_terms_leave_force_field_unchanged_for_empty_input() {
    let mut ff = distgeom_experimental_torsion_force_field();
    let details = CrystalFFDetails::default();
    let mut atom_pairs = vec![false; 25];

    add_experimental_torsion_terms(&mut ff, &details, &mut atom_pairs, 5);
    ff.initialize();
    let mut contrib_energies = Vec::new();

    assert_eq!(ff.calc_energy_current(Some(&mut contrib_energies)), 0.0);
    assert!(contrib_energies.is_empty());
    assert_eq!(atom_pairs, vec![false; 25]);
}

#[test]
#[should_panic]
fn distgeom_experimental_torsion_terms_reject_short_torsion_atom_record() {
    let mut ff = distgeom_experimental_torsion_force_field();
    let details = CrystalFFDetails {
        exp_torsion_atoms: vec![vec![0, 1, 2]],
        exp_torsion_angles: vec![(vec![1; 6], vec![0.1; 6])],
        ..CrystalFFDetails::default()
    };
    let mut atom_pairs = vec![false; 25];

    add_experimental_torsion_terms(&mut ff, &details, &mut atom_pairs, 5);
}

#[test]
fn distgeom_add12_terms_mark_ordered_bond_pairs_and_add_distance_contribs() {
    let mut ff = distgeom_experimental_torsion_force_field();
    let positions = vec![
        ForceFieldVec3::new(0.0, 0.0, 0.0),
        ForceFieldVec3::new(2.0, 0.0, 0.0),
        ForceFieldVec3::new(1.0, 1.0, 0.0),
        ForceFieldVec3::new(1.0, 1.0, 1.0),
        ForceFieldVec3::new(2.0, 3.0, 0.0),
    ];
    let details = CrystalFFDetails {
        bonds: vec![(0, 1), (4, 1)],
        ..CrystalFFDetails::default()
    };
    let mut atom_pairs = vec![false; 25];

    add_12_terms(&mut ff, &details, &mut atom_pairs, &positions, 7.0, 5);
    ff.initialize();
    let mut contrib_energies = Vec::new();
    let energy = ff.calc_energy_current(Some(&mut contrib_energies));

    let mut expected_ff = distgeom_experimental_torsion_force_field();
    expected_ff.initialize();
    let mut expected = crate::chemistry::forcefield::DistanceConstraintContribs::new(&expected_ff);
    expected.add_contrib(0, 1, 1.99, 2.01, 7.0);
    let d41 = (positions[4] - positions[1]).length();
    expected.add_contrib(4, 1, d41 - KNOWN_DIST_TOL, d41 + KNOWN_DIST_TOL, 7.0);
    let expected_energy = expected.get_energy(&flatten_force_field_positions(&expected_ff));

    assert!(atom_pairs[1]);
    assert!(atom_pairs[9]);
    assert_eq!(atom_pairs.iter().filter(|&&set| set).count(), 2);
    assert_eq!(contrib_energies.len(), 1);
    assert_close(contrib_energies[0], expected_energy);
    assert_close(energy, expected_energy);
}

#[test]
fn distgeom_add12_terms_leave_force_field_unchanged_for_empty_input() {
    let mut ff = distgeom_experimental_torsion_force_field();
    let positions = ff.positions().to_vec();
    let details = CrystalFFDetails::default();
    let mut atom_pairs = vec![false; 25];

    add_12_terms(&mut ff, &details, &mut atom_pairs, &positions, 3.0, 5);
    ff.initialize();
    let mut contrib_energies = Vec::new();

    assert_eq!(ff.calc_energy_current(Some(&mut contrib_energies)), 0.0);
    assert!(contrib_energies.is_empty());
    assert_eq!(atom_pairs, vec![false; 25]);
}

#[test]
#[should_panic]
fn distgeom_add12_terms_reject_out_of_range_position_index() {
    let mut ff = distgeom_experimental_torsion_force_field();
    let positions = ff.positions().to_vec();
    let details = CrystalFFDetails {
        bonds: vec![(0, 5)],
        ..CrystalFFDetails::default()
    };
    let mut atom_pairs = vec![false; 30];

    add_12_terms(&mut ff, &details, &mut atom_pairs, &positions, 3.0, 6);
}

#[test]
fn distgeom_add13_terms_add_angle_improper_bounds_and_current_distance_contribs() {
    let mut ff = distgeom_experimental_torsion_force_field();
    let positions = vec![
        ForceFieldVec3::new(0.0, 0.0, 0.0),
        ForceFieldVec3::new(2.0, 0.0, 0.0),
        ForceFieldVec3::new(1.0, 1.0, 0.0),
        ForceFieldVec3::new(1.0, 1.0, 1.0),
        ForceFieldVec3::new(2.0, 3.0, 0.0),
    ];
    let details = CrystalFFDetails {
        angles: vec![vec![0, 1, 2, 1], vec![0, 1, 3, 0], vec![4, 3, 1, 0]],
        ..CrystalFFDetails::default()
    };
    let mut atom_pairs = vec![false; 25];
    let mut mmat = BoundsMatrix::new(5);
    mmat.set_lower(0, 3, 1.5).expect("set lower");
    mmat.set_upper(0, 3, 1.6).expect("set upper");
    let is_improper_constrained = vec![false, true, false, false, false];

    add_13_terms(
        &mut ff,
        &details,
        &mut atom_pairs,
        &positions,
        4.0,
        &is_improper_constrained,
        true,
        &mmat,
        5,
    );
    ff.initialize();
    let mut contrib_energies = Vec::new();
    let energy = ff.calc_energy_current(Some(&mut contrib_energies));

    let mut expected_ff = distgeom_experimental_torsion_force_field();
    expected_ff.initialize();
    let mut expected_angle =
        crate::chemistry::forcefield::AngleConstraintContribs::new(&expected_ff);
    expected_angle.add_contrib(0, 1, 2, 179.0, 180.0, 1.0);
    let expected_angle_energy =
        expected_angle.get_energy(&flatten_force_field_positions(&expected_ff));

    let mut expected_dist =
        crate::chemistry::forcefield::DistanceConstraintContribs::new(&expected_ff);
    expected_dist.add_contrib(0, 3, 1.5, 1.6, 4.0);
    let d41 = (positions[4] - positions[1]).length();
    expected_dist.add_contrib(4, 1, d41 - KNOWN_DIST_TOL, d41 + KNOWN_DIST_TOL, 4.0);
    let expected_dist_energy =
        expected_dist.get_energy(&flatten_force_field_positions(&expected_ff));

    assert!(atom_pairs[2]);
    assert!(atom_pairs[3]);
    assert!(atom_pairs[9]);
    assert_eq!(atom_pairs.iter().filter(|&&set| set).count(), 3);
    assert_eq!(contrib_energies.len(), 2);
    assert_close(contrib_energies[0], expected_angle_energy);
    assert_close(contrib_energies[1], expected_dist_energy);
    assert_close(energy, expected_angle_energy + expected_dist_energy);
}

#[test]
fn distgeom_add13_terms_leave_force_field_unchanged_for_empty_input() {
    let mut ff = distgeom_experimental_torsion_force_field();
    let positions = ff.positions().to_vec();
    let details = CrystalFFDetails::default();
    let mut atom_pairs = vec![false; 25];
    let mmat = BoundsMatrix::new(5);
    let is_improper_constrained = vec![false; 5];

    add_13_terms(
        &mut ff,
        &details,
        &mut atom_pairs,
        &positions,
        4.0,
        &is_improper_constrained,
        true,
        &mmat,
        5,
    );
    ff.initialize();
    let mut contrib_energies = Vec::new();

    assert_eq!(ff.calc_energy_current(Some(&mut contrib_energies)), 0.0);
    assert!(contrib_energies.is_empty());
    assert_eq!(atom_pairs, vec![false; 25]);
}

#[test]
#[should_panic]
fn distgeom_add13_terms_reject_short_angle_record() {
    let mut ff = distgeom_experimental_torsion_force_field();
    let positions = ff.positions().to_vec();
    let details = CrystalFFDetails {
        angles: vec![vec![0, 1, 2]],
        ..CrystalFFDetails::default()
    };
    let mut atom_pairs = vec![false; 25];
    let mmat = BoundsMatrix::new(5);
    let is_improper_constrained = vec![false; 5];

    add_13_terms(
        &mut ff,
        &details,
        &mut atom_pairs,
        &positions,
        4.0,
        &is_improper_constrained,
        true,
        &mmat,
        5,
    );
}

#[test]
fn distgeom_long_range_distance_constraints_skip_atom_pairs_and_use_bounds_or_tight_constraints() {
    let mut ff = distgeom_experimental_torsion_force_field();
    let positions = vec![
        ForceFieldVec3::new(0.0, 0.0, 0.0),
        ForceFieldVec3::new(2.0, 0.0, 0.0),
        ForceFieldVec3::new(1.0, 1.0, 0.0),
        ForceFieldVec3::new(1.0, 1.0, 1.0),
    ];
    let mut mmat = BoundsMatrix::new(4);
    for i in 1..4 {
        for j in 0..i {
            mmat.set_lower(i, j, 0.5 + i as f64 + j as f64 * 0.1)
                .expect("set lower");
            mmat.set_upper(i, j, 2.5 + i as f64 + j as f64 * 0.1)
                .expect("set upper");
        }
    }
    let details = CrystalFFDetails {
        bounds_mat_force_scaling: 2.5,
        constrained_atoms: vec![true, true, false, true],
        ..CrystalFFDetails::default()
    };
    let mut atom_pairs = vec![false; 16];
    atom_pairs[2] = true;

    add_long_range_distance_constraints(&mut ff, &details, &atom_pairs, &positions, 8.0, &mmat, 4);
    ff.initialize();
    let mut contrib_energies = Vec::new();
    let energy = ff.calc_energy_current(Some(&mut contrib_energies));

    let mut expected_ff = distgeom_experimental_torsion_force_field();
    expected_ff.initialize();
    let mut expected = crate::chemistry::forcefield::DistanceConstraintContribs::new(&expected_ff);
    let d10 = (positions[1] - positions[0]).length();
    expected.add_contrib(1, 0, d10 - KNOWN_DIST_TOL, d10 + KNOWN_DIST_TOL, 8.0);
    expected.add_contrib(2, 1, mmat.get_lower(2, 1), mmat.get_upper(2, 1), 25.0);
    let d30 = (positions[3] - positions[0]).length();
    expected.add_contrib(3, 0, d30 - KNOWN_DIST_TOL, d30 + KNOWN_DIST_TOL, 8.0);
    let d31 = (positions[3] - positions[1]).length();
    expected.add_contrib(3, 1, d31 - KNOWN_DIST_TOL, d31 + KNOWN_DIST_TOL, 8.0);
    expected.add_contrib(3, 2, mmat.get_lower(3, 2), mmat.get_upper(3, 2), 25.0);
    let expected_energy = expected.get_energy(&flatten_force_field_positions(&expected_ff));

    assert_eq!(contrib_energies.len(), 1);
    assert_close(contrib_energies[0], expected_energy);
    assert_close(energy, expected_energy);
}

#[test]
fn distgeom_long_range_distance_constraints_leave_force_field_unchanged_when_all_pairs_present() {
    let mut ff = distgeom_experimental_torsion_force_field();
    let positions = ff.positions().to_vec();
    let details = CrystalFFDetails {
        bounds_mat_force_scaling: 2.5,
        ..CrystalFFDetails::default()
    };
    let atom_pairs = vec![true; 25];
    let mmat = BoundsMatrix::new(5);

    add_long_range_distance_constraints(&mut ff, &details, &atom_pairs, &positions, 8.0, &mmat, 5);
    ff.initialize();
    let mut contrib_energies = Vec::new();

    assert_eq!(ff.calc_energy_current(Some(&mut contrib_energies)), 0.0);
    assert!(contrib_energies.is_empty());
}

#[test]
#[should_panic]
fn distgeom_long_range_distance_constraints_reject_short_atom_pair_bitset() {
    let mut ff = distgeom_experimental_torsion_force_field();
    let positions = ff.positions().to_vec();
    let details = CrystalFFDetails::default();
    let atom_pairs = vec![false; 3];
    let mmat = BoundsMatrix::new(5);

    add_long_range_distance_constraints(&mut ff, &details, &atom_pairs, &positions, 8.0, &mmat, 5);
}

fn construct_3d_forcefield_positions() -> Vec<ForceFieldVec3> {
    vec![
        ForceFieldVec3::new(0.0, 0.0, 0.0),
        ForceFieldVec3::new(1.2, 0.0, 0.0),
        ForceFieldVec3::new(1.2, 1.1, 0.0),
        ForceFieldVec3::new(1.2, 1.1, 1.3),
        ForceFieldVec3::new(2.4, 1.1, 1.3),
    ]
}

fn construct_3d_forcefield_bounds_matrix(num_atoms: usize) -> BoundsMatrix {
    let mut mmat = BoundsMatrix::new(num_atoms);
    for i in 1..num_atoms {
        for j in 0..i {
            mmat.set_lower(i, j, 0.75 + i as f64 * 0.1 + j as f64 * 0.01)
                .expect("set lower");
            mmat.set_upper(i, j, 3.25 + i as f64 * 0.1 + j as f64 * 0.01)
                .expect("set upper");
        }
    }
    mmat
}

fn construct_3d_forcefield_details() -> CrystalFFDetails {
    CrystalFFDetails {
        exp_torsion_atoms: vec![vec![0, 1, 2, 3]],
        exp_torsion_angles: vec![(
            vec![1, -1, 1, -1, 1, -1],
            vec![0.10, 0.20, 0.30, 0.40, 0.50, 0.60],
        )],
        improper_atoms: vec![vec![1, 2, 3, 4, 6, 0]],
        bonds: vec![(0, 1), (3, 4)],
        angles: vec![vec![0, 1, 2, 1], vec![1, 2, 4, 0]],
        bounds_mat_force_scaling: 1.7,
        constrained_atoms: vec![true, true, false, true, true],
        ..CrystalFFDetails::default()
    }
}

#[test]
fn construct_3d_forcefield_reproduces_source_helper_sequence_and_positions() {
    let positions = construct_3d_forcefield_positions();
    let mmat = construct_3d_forcefield_bounds_matrix(positions.len());
    let details = construct_3d_forcefield_details();

    let mut field = construct_3d_forcefield(&mmat, &positions, &details);
    field.initialize();
    let mut contrib_energies = Vec::new();
    let energy = field.calc_energy_current(Some(&mut contrib_energies));

    let mut expected = ForceField::new(3);
    expected.positions_mut().extend_from_slice(&positions);
    let mut atom_pairs = vec![false; positions.len() * positions.len()];
    let mut is_improper_constrained = vec![false; positions.len()];
    add_experimental_torsion_terms(&mut expected, &details, &mut atom_pairs, positions.len());
    add_improper_torsion_terms(
        &mut expected,
        10.0,
        &details.improper_atoms,
        &mut is_improper_constrained,
    );
    add_12_terms(
        &mut expected,
        &details,
        &mut atom_pairs,
        &positions,
        KNOWN_DIST_FORCE_CONSTANT,
        positions.len(),
    );
    add_13_terms(
        &mut expected,
        &details,
        &mut atom_pairs,
        &positions,
        KNOWN_DIST_FORCE_CONSTANT,
        &is_improper_constrained,
        true,
        &mmat,
        positions.len(),
    );
    add_long_range_distance_constraints(
        &mut expected,
        &details,
        &atom_pairs,
        &positions,
        KNOWN_DIST_FORCE_CONSTANT,
        &mmat,
        positions.len(),
    );
    expected.initialize();
    let mut expected_contrib_energies = Vec::new();
    let expected_energy = expected.calc_energy_current(Some(&mut expected_contrib_energies));

    assert_eq!(field.dimension(), 3);
    assert_eq!(field.positions(), positions);
    assert_eq!(contrib_energies.len(), expected_contrib_energies.len());
    assert_eq!(contrib_energies.len(), 6);
    for (observed, expected) in contrib_energies
        .iter()
        .zip(expected_contrib_energies.iter())
    {
        assert_close(*observed, *expected);
    }
    assert_close(energy, expected_energy);
}

#[test]
#[should_panic]
fn construct_3d_forcefield_rejects_bounds_position_size_mismatch() {
    let positions = construct_3d_forcefield_positions();
    let mmat = BoundsMatrix::new(positions.len() + 1);
    let details = CrystalFFDetails::default();

    let _ = construct_3d_forcefield(&mmat, &positions, &details);
}

#[test]
#[should_panic]
fn construct_3d_forcefield_rejects_torsion_atom_angle_size_mismatch() {
    let positions = construct_3d_forcefield_positions();
    let mmat = BoundsMatrix::new(positions.len());
    let details = CrystalFFDetails {
        exp_torsion_atoms: vec![vec![0, 1, 2, 3]],
        exp_torsion_angles: Vec::new(),
        ..CrystalFFDetails::default()
    };

    let _ = construct_3d_forcefield(&mmat, &positions, &details);
}

#[test]
fn construct_3d_forcefield_with_cpci_appends_electrostatic_terms_after_base_field() {
    let positions = vec![
        ForceFieldVec3::new(0.0, 0.0, 0.0),
        ForceFieldVec3::new(1.0, 0.0, 0.0),
        ForceFieldVec3::new(0.0, 2.0, 0.0),
    ];
    let mmat = construct_3d_forcefield_bounds_matrix(positions.len());
    let details = CrystalFFDetails::default();
    let base = construct_3d_forcefield(&mmat, &positions, &details);
    let cpci = std::collections::BTreeMap::from([((0, 1), 0.5), ((1, 2), -0.25)]);

    let mut field = construct_3d_forcefield_with_cpci(&mmat, &positions, &details, &cpci);
    field.initialize();
    let mut contrib_energies = Vec::new();
    let energy = field.calc_energy_current(Some(&mut contrib_energies));

    let mut expected_cpci_field = ForceField::new(3);
    expected_cpci_field
        .positions_mut()
        .extend_from_slice(&positions);
    expected_cpci_field.initialize();
    let mut expected_cpci =
        crate::chemistry::forcefield::mmff::nonbonded::NonbondedContrib::new(&expected_cpci_field);
    expected_cpci.add_term(0, 1, None, true, 0.5, 1, false);
    expected_cpci.add_term(1, 2, None, true, -0.25, 1, false);
    let expected_cpci_energy =
        expected_cpci.get_energy(&flatten_force_field_positions(&expected_cpci_field));

    let mut base = base;
    base.initialize();
    let mut base_contrib_energies = Vec::new();
    let base_energy = base.calc_energy_current(Some(&mut base_contrib_energies));

    assert_eq!(contrib_energies.len(), base_contrib_energies.len() + 1);
    for (observed, expected) in contrib_energies
        .iter()
        .take(base_contrib_energies.len())
        .zip(base_contrib_energies.iter())
    {
        assert_close(*observed, *expected);
    }
    assert_close(
        *contrib_energies.last().expect("CPCI contribution"),
        expected_cpci_energy,
    );
    assert_close(energy, base_energy + expected_cpci_energy);
}

#[test]
fn construct_plain_3d_forcefield_reproduces_source_helper_sequence_without_improper_knowledge() {
    let positions = construct_3d_forcefield_positions();
    let mmat = construct_3d_forcefield_bounds_matrix(positions.len());
    let details = construct_3d_forcefield_details();

    let mut field = construct_plain_3d_forcefield(&mmat, &positions, &details);
    field.initialize();
    let mut contrib_energies = Vec::new();
    let energy = field.calc_energy_current(Some(&mut contrib_energies));

    let mut expected = ForceField::new(3);
    expected.positions_mut().extend_from_slice(&positions);
    let mut atom_pairs = vec![false; positions.len() * positions.len()];
    let is_improper_constrained = vec![false; positions.len()];
    add_experimental_torsion_terms(&mut expected, &details, &mut atom_pairs, positions.len());
    add_12_terms(
        &mut expected,
        &details,
        &mut atom_pairs,
        &positions,
        KNOWN_DIST_FORCE_CONSTANT,
        positions.len(),
    );
    add_13_terms(
        &mut expected,
        &details,
        &mut atom_pairs,
        &positions,
        KNOWN_DIST_FORCE_CONSTANT,
        &is_improper_constrained,
        false,
        &mmat,
        positions.len(),
    );
    add_long_range_distance_constraints(
        &mut expected,
        &details,
        &atom_pairs,
        &positions,
        KNOWN_DIST_FORCE_CONSTANT,
        &mmat,
        positions.len(),
    );
    expected.initialize();
    let mut expected_contrib_energies = Vec::new();
    let expected_energy = expected.calc_energy_current(Some(&mut expected_contrib_energies));

    assert_eq!(field.dimension(), 3);
    assert_eq!(field.positions(), positions);
    assert_eq!(contrib_energies.len(), expected_contrib_energies.len());
    assert_eq!(contrib_energies.len(), 4);
    for (observed, expected) in contrib_energies
        .iter()
        .zip(expected_contrib_energies.iter())
    {
        assert_close(*observed, *expected);
    }
    assert_close(energy, expected_energy);
}

#[test]
#[should_panic]
fn construct_plain_3d_forcefield_rejects_bounds_position_size_mismatch() {
    let positions = construct_3d_forcefield_positions();
    let mmat = BoundsMatrix::new(positions.len() + 1);
    let details = CrystalFFDetails::default();

    let _ = construct_plain_3d_forcefield(&mmat, &positions, &details);
}

#[test]
#[should_panic]
fn construct_plain_3d_forcefield_rejects_torsion_atom_angle_size_mismatch() {
    let positions = construct_3d_forcefield_positions();
    let mmat = BoundsMatrix::new(positions.len());
    let details = CrystalFFDetails {
        exp_torsion_atoms: vec![vec![0, 1, 2, 3]],
        exp_torsion_angles: Vec::new(),
        ..CrystalFFDetails::default()
    };

    let _ = construct_plain_3d_forcefield(&mmat, &positions, &details);
}

#[test]
fn construct_3d_improper_forcefield_reproduces_parts_overload_terms() {
    let positions = construct_3d_forcefield_positions();
    let mmat = construct_3d_forcefield_bounds_matrix(positions.len());
    let improper_atoms = vec![vec![1, 2, 3, 4, 6, 0]];
    let angles = vec![vec![0, 1, 2, 1], vec![1, 2, 4, 0], vec![4, 3, 2, 1]];
    let atom_nums = vec![6, 7, 8, 9, 16];

    let mut field = construct_3d_improper_forcefield_from_parts(
        &mmat,
        &positions,
        &improper_atoms,
        &angles,
        &atom_nums,
    );
    field.initialize();
    let mut contrib_energies = Vec::new();
    let energy = field.calc_energy_current(Some(&mut contrib_energies));

    let mut expected = ForceField::new(3);
    expected.positions_mut().extend_from_slice(&positions);
    let mut is_improper_constrained = vec![false; positions.len()];
    add_improper_torsion_terms(
        &mut expected,
        10.0,
        &improper_atoms,
        &mut is_improper_constrained,
    );
    let mut expected_angles = crate::chemistry::forcefield::AngleConstraintContribs::new(&expected);
    expected_angles.add_contrib(0, 1, 2, 179.0, 180.0, 10.0);
    expected_angles.add_contrib(4, 3, 2, 179.0, 180.0, 10.0);
    expected.add_contrib(Box::new(expected_angles));
    expected.initialize();
    let mut expected_contrib_energies = Vec::new();
    let expected_energy = expected.calc_energy_current(Some(&mut expected_contrib_energies));

    assert_eq!(field.dimension(), 3);
    assert_eq!(field.positions(), positions);
    assert_eq!(contrib_energies.len(), expected_contrib_energies.len());
    assert_eq!(contrib_energies.len(), 2);
    for (observed, expected) in contrib_energies
        .iter()
        .zip(expected_contrib_energies.iter())
    {
        assert_close(*observed, *expected);
    }
    assert_close(energy, expected_energy);
}

#[test]
fn construct_3d_improper_forcefield_details_overload_delegates_to_parts_and_ignores_atom_nums() {
    let positions = construct_3d_forcefield_positions();
    let mmat = construct_3d_forcefield_bounds_matrix(positions.len());
    let details = CrystalFFDetails {
        improper_atoms: vec![vec![1, 2, 3, 4, 6, 0]],
        angles: vec![vec![0, 1, 2, 1], vec![1, 2, 4, 0], vec![4, 3, 2, 1]],
        atom_nums: vec![6, 7, 8, 9, 16],
        ..CrystalFFDetails::default()
    };

    let mut field = construct_3d_improper_forcefield(&mmat, &positions, &details);
    let mut expected = construct_3d_improper_forcefield_from_parts(
        &mmat,
        &positions,
        &details.improper_atoms,
        &details.angles,
        &[999, 998],
    );
    field.initialize();
    expected.initialize();
    let mut contrib_energies = Vec::new();
    let mut expected_contrib_energies = Vec::new();

    let energy = field.calc_energy_current(Some(&mut contrib_energies));
    let expected_energy = expected.calc_energy_current(Some(&mut expected_contrib_energies));

    assert_eq!(contrib_energies.len(), expected_contrib_energies.len());
    for (observed, expected) in contrib_energies
        .iter()
        .zip(expected_contrib_energies.iter())
    {
        assert_close(*observed, *expected);
    }
    assert_close(energy, expected_energy);
}

#[test]
#[should_panic]
fn construct_3d_improper_forcefield_rejects_bounds_position_size_mismatch() {
    let positions = construct_3d_forcefield_positions();
    let mmat = BoundsMatrix::new(positions.len() + 1);
    let details = CrystalFFDetails::default();

    let _ = construct_3d_improper_forcefield(&mmat, &positions, &details);
}

#[test]
#[should_panic]
fn construct_3d_improper_forcefield_rejects_short_angle_record() {
    let positions = construct_3d_forcefield_positions();
    let mmat = BoundsMatrix::new(positions.len());
    let details = CrystalFFDetails {
        angles: vec![vec![0, 1, 2]],
        ..CrystalFFDetails::default()
    };

    let _ = construct_3d_improper_forcefield(&mmat, &positions, &details);
}

#[test]
fn compute_initial_coords_power_solver_follows_source_seeded_diagonal_iteration() {
    let mut mat = symm_matrix_from_distances(3, &[(0, 0, 5.0), (1, 1, 3.0), (2, 2, 1.0)]);
    let mut eigen_values = vec![0.0; 2];
    let mut eigen_vectors = DoubleMatrix::new(2, 3);

    assert!(
        power_eigen_solver(2, &mut mat, &mut eigen_values, Some(&mut eigen_vectors), 23)
            .expect("power solver")
    );

    assert!((eigen_values[0] - 3.0).abs() < 1.0e-3);
    assert!((eigen_values[1] - 1.006_175_852_283_403_2).abs() < 1.0e-3);
}

fn run_set12_bounds(mol: &Molecule) -> (BoundsMatrix, ComputedData) {
    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let mut accum_data = ComputedData::new(mol.num_atoms(), mol.num_bonds());
    set_12_bounds(mol, &mut mmat, &mut accum_data).expect("set12Bounds");
    (mmat, accum_data)
}

#[test]
fn set12_bounds_uff_total_valence_includes_bracket_explicit_hydrogen() {
    let mol = Molecule::from_smiles("C[SH+][O-]").expect("explicit sulfur hydrogen");
    let assignment = assign_valence(&mol, ValenceModel::RdkitLike).expect("valence");

    assert_eq!(assignment.explicit_valence[1], 3);
    assert_eq!(assignment.implicit_hydrogens[1], 0);
    assert_eq!(atom_total_valence_for_uff(&assignment, 1), 3);

    let (bounds, accum_data) = run_set12_bounds(&mol);
    assert!((accum_data.bond_lengths[0] - 1.776_838_813_449_356_2).abs() < 1.0e-14);
    assert!((accum_data.bond_lengths[1] - 1.679_472_654_030_108).abs() < 1.0e-14);
    assert!((bounds.get_upper(0, 1) - 1.786_838_813_449_356_2).abs() < 1.0e-14);
    assert!((bounds.get_upper(1, 2) - 1.689_472_654_030_108).abs() < 1.0e-14);
}

fn run_set13_bounds(mol: &Molecule) -> (BoundsMatrix, ComputedData) {
    let (mut mmat, mut accum_data) = run_set12_bounds(mol);
    let rinfo = ring_info_for_distgeom(mol).expect("ring info");
    set_13_bounds(mol, &mut mmat, &mut accum_data, &rinfo).expect("set13Bounds");
    (mmat, accum_data)
}

fn run_set14_bounds(
    mol: &Molecule,
    use_macrocycle_14config: bool,
    force_trans_amides: bool,
) -> (BoundsMatrix, ComputedData, Vec<f64>) {
    let (mut mmat, mut accum_data) = run_set13_bounds(mol);
    let dmat = flatten_topological_distances(mol);
    set_14_bounds(
        mol,
        &mut mmat,
        &mut accum_data,
        &dmat,
        use_macrocycle_14config,
        force_trans_amides,
        &ring_info_for_distgeom(mol).expect("ring info"),
    )
    .expect("set14Bounds");
    (mmat, accum_data, dmat)
}

fn run_set15_bounds(
    mol: &Molecule,
    use_macrocycle_14config: bool,
    force_trans_amides: bool,
) -> (BoundsMatrix, ComputedData, Vec<f64>) {
    let (mut mmat, mut accum_data, dmat) =
        run_set14_bounds(mol, use_macrocycle_14config, force_trans_amides);
    set_15_bounds(mol, &mut mmat, &mut accum_data, &dmat).expect("set15Bounds");
    (mmat, accum_data, dmat)
}

fn run_set_topol_bounds(
    mol: &Molecule,
    set15bounds: bool,
    scale_vdw: bool,
    use_macrocycle_14config: bool,
    force_trans_amides: bool,
    set14bounds: bool,
    set13bounds: bool,
) -> BoundsMatrix {
    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    set_topol_bounds(
        mol,
        &mut mmat,
        set15bounds,
        scale_vdw,
        use_macrocycle_14config,
        force_trans_amides,
        set14bounds,
        set13bounds,
    )
    .expect("setTopolBounds");
    mmat
}

fn run_set_topol_bounds_with_outputs(
    mol: &Molecule,
    set15bounds: bool,
    scale_vdw: bool,
    use_macrocycle_14config: bool,
    force_trans_amides: bool,
    set14bounds: bool,
    set13bounds: bool,
) -> (BoundsMatrix, Vec<(i32, i32)>, Vec<Vec<i32>>) {
    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let mut bonds = Vec::new();
    let mut angles = Vec::new();
    set_topol_bounds_with_outputs(
        mol,
        &mut mmat,
        &mut bonds,
        &mut angles,
        set15bounds,
        scale_vdw,
        use_macrocycle_14config,
        force_trans_amides,
        set14bounds,
        set13bounds,
    )
    .expect("setTopolBounds with outputs");
    (mmat, bonds, angles)
}

fn run_set14_same_ring_pass_only(
    mol: &Molecule,
    use_macrocycle_14config: bool,
) -> (BoundsMatrix, ComputedData, Vec<f64>) {
    let (mut mmat, mut accum_data) = run_set13_bounds(mol);
    let dmat = flatten_topological_distances(mol);
    let rinfo = ring_info_for_distgeom(mol).expect("ring info");
    for bring in rinfo.bond_rings() {
        let r_size = bring.len();
        if r_size < 3 {
            continue;
        }
        let mut bid1 = bring[r_size - 1].index();
        for i in 0..r_size {
            let bid2 = bring[i].index();
            let bid3 = bring[(i + 1) % r_size].index();
            if r_size > 5 {
                if use_macrocycle_14config && r_size >= MIN_MACROCYCLE_RING_SIZE {
                    set_macrocycle_all_in_same_ring_14_bounds(
                        mol,
                        bid1,
                        bid2,
                        bid3,
                        &mut accum_data,
                        &mut mmat,
                    )
                    .expect("macrocycle same-ring 1-4 bounds");
                } else {
                    set_in_ring_14_bounds(
                        mol,
                        bid1,
                        bid2,
                        bid3,
                        &mut accum_data,
                        &mut mmat,
                        &dmat,
                        r_size,
                        &rinfo,
                    )
                    .expect("in-ring 1-4 bounds");
                }
            } else {
                record_14_path(mol, bid1, bid2, bid3, &mut accum_data).expect("record14Path");
            }
            bid1 = bid2;
        }
    }
    (mmat, accum_data, dmat)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Set14DispatchCase {
    TwoSameRing,
    TwoDiffRing,
    ShareRingBond,
    Chain,
}

fn find_same_ring_dispatch_triple(
    mol: &Molecule,
    use_macrocycle_14config: bool,
) -> Option<(usize, usize, usize, usize)> {
    let rinfo = ring_info_for_distgeom(mol).expect("ring info");
    for bring in rinfo.bond_rings() {
        let r_size = bring.len();
        if r_size < 3 {
            continue;
        }
        if r_size > 5 && (!use_macrocycle_14config || r_size >= MIN_MACROCYCLE_RING_SIZE) {
            let bid1 = bring[r_size - 1].index();
            let bid2 = bring[0].index();
            let bid3 = bring[1].index();
            return Some((bid1, bid2, bid3, r_size));
        }
    }
    None
}

fn find_dispatch_triple(
    mol: &Molecule,
    target: Set14DispatchCase,
    use_macrocycle_14config: bool,
) -> Option<(usize, usize, usize)> {
    let rinfo = ring_info_for_distgeom(mol).expect("ring info");
    let mut bid_is_macrocycle: HashSet<usize> = HashSet::new();
    let mut ring_bond_pairs: HashSet<u64> = HashSet::new();
    let mut done_paths: HashSet<u64> = HashSet::new();
    let nb = mol.num_bonds() as u64;

    for bring in rinfo.bond_rings() {
        let r_size = bring.len();
        if r_size < 3 {
            continue;
        }
        let mut bid1 = bring[r_size - 1].index();
        for i in 0..r_size {
            let bid2 = bring[i].index();
            let bid3 = bring[(i + 1) % r_size].index();
            ring_bond_pairs.insert(bid1 as u64 * nb + bid2 as u64);
            ring_bond_pairs.insert(bid2 as u64 * nb + bid1 as u64);
            done_paths.insert(bid1 as u64 * nb * nb + bid2 as u64 * nb + bid3 as u64);
            done_paths.insert(bid3 as u64 * nb * nb + bid2 as u64 * nb + bid1 as u64);
            if use_macrocycle_14config && r_size >= MIN_MACROCYCLE_RING_SIZE {
                bid_is_macrocycle.insert(bid2);
            }
            bid1 = bid2;
        }
    }

    for bond in mol.bonds() {
        let bid2 = bond.id().index();
        let aid2 = bond.begin().index();
        let aid3 = bond.end().index();
        for nbr1 in neighbors_for_atom(mol, aid2) {
            let Some(bid1) = bond_between_idx_simple(mol, aid2, nbr1) else {
                continue;
            };
            if bid1 == bid2 {
                continue;
            }
            for nbr3 in neighbors_for_atom(mol, aid3) {
                let Some(bid3) = bond_between_idx_simple(mol, aid3, nbr3) else {
                    continue;
                };
                if bid3 == bid2 {
                    continue;
                }
                let id1 = bid1 as u64 * nb * nb + bid2 as u64 * nb + bid3 as u64;
                let id2 = bid3 as u64 * nb * nb + bid2 as u64 * nb + bid1 as u64;
                if done_paths.contains(&id1) || done_paths.contains(&id2) {
                    continue;
                }
                let pid1 = bid1 as u64 * nb + bid2 as u64;
                let pid2 = bid2 as u64 * nb + bid1 as u64;
                let pid3 = bid2 as u64 * nb + bid3 as u64;
                let pid4 = bid3 as u64 * nb + bid2 as u64;
                let case = if ring_bond_pairs.contains(&pid1)
                    || ring_bond_pairs.contains(&pid2)
                    || ring_bond_pairs.contains(&pid3)
                    || ring_bond_pairs.contains(&pid4)
                {
                    Set14DispatchCase::TwoSameRing
                } else if (rinfo.num_bond_rings(BondId::new(bid1)) > 0
                    && rinfo.num_bond_rings(BondId::new(bid2)) > 0)
                    || (rinfo.num_bond_rings(BondId::new(bid2)) > 0
                        && rinfo.num_bond_rings(BondId::new(bid3)) > 0)
                {
                    Set14DispatchCase::TwoDiffRing
                } else if rinfo.num_bond_rings(BondId::new(bid2)) > 0 {
                    Set14DispatchCase::ShareRingBond
                } else {
                    Set14DispatchCase::Chain
                };
                if case == target {
                    return Some((bid1, bid2, bid3));
                }
            }
        }
    }
    let _ = bid_is_macrocycle;
    None
}

fn flatten_topological_distances(mol: &Molecule) -> Vec<f64> {
    compute_topological_distances(mol).to_vec()
}

fn single_atom_molecule(spec: AtomSpec) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    builder.add_atom(spec);
    builder.build().expect("single atom molecule")
}

fn uff_atom_label(mol: &Molecule, atom_index: usize) -> String {
    let assignment = assign_valence(mol, ValenceModel::RdkitLike).expect("valence");
    let mut atom_degree = vec![0usize; mol.num_atoms()];
    for bond in mol.bonds() {
        atom_degree[bond.begin().index()] += 1;
        atom_degree[bond.end().index()] += 1;
    }
    let conjugated = compute_conjugated_bonds_for_uff(mol, &assignment, &atom_degree);
    let mut atom_has_conjugated_bond = vec![false; mol.num_atoms()];
    for (bond_index, bond) in mol.bonds().iter().enumerate() {
        if conjugated[bond_index] {
            atom_has_conjugated_bond[bond.begin().index()] = true;
            atom_has_conjugated_bond[bond.end().index()] = true;
        }
    }
    let hybridizations =
        compute_hybridizations_for_uff(mol, &assignment, &atom_degree, &atom_has_conjugated_bond);
    get_atom_label_for_uff(
        mol,
        &assignment,
        &hybridizations,
        &atom_has_conjugated_bond,
        atom_index,
    )
    .expect("UFF atom label")
}

#[test]
fn topological_distance_matches_rdkit_shortest_path_lengths_for_chain() {
    let mol = Molecule::from_smiles("CCCC").expect("butane");
    let dist = compute_topological_distances(&mol);
    let n = mol.num_atoms();

    assert_eq!(dist[0], 0.0);
    assert_eq!(dist[1], 1.0);
    assert_eq!(dist[2], 2.0);
    assert_eq!(dist[3], 3.0);
    assert_eq!(dist[3 * n], 3.0);
    assert_eq!(dist[1 * n + 3], 2.0);
}

#[test]
fn topological_distance_keeps_local_inf_for_disconnected_pairs() {
    let mut builder = MoleculeBuilder::new();
    builder.add_atom(AtomSpec::new(Element::C));
    builder.add_atom(AtomSpec::new(Element::O));
    let mol = builder.build().expect("disconnected atoms");
    let dist = compute_topological_distances(&mol);

    assert_eq!(dist[0], 0.0);
    assert_eq!(dist[3], 0.0);
    assert_eq!(dist[1], LOCAL_INF_DIST);
    assert_eq!(dist[2], LOCAL_INF_DIST);
}

#[test]
fn collect_bonds_and_angles_flags_triple_bond_paths_like_rdkit() {
    let mol = Molecule::from_smiles("CC#N").expect("acetonitrile");
    let mut bonds = Vec::new();
    let mut angles = Vec::new();

    collect_bonds_and_angles(&mol, &mut bonds, &mut angles);

    assert_eq!(bonds, vec![(0, 1), (1, 2)]);
    assert_eq!(angles, vec![vec![0, 1, 2, 1]]);
}

#[test]
fn collect_bonds_and_angles_flags_consecutive_double_bonds_like_rdkit() {
    let mol = Molecule::from_smiles("C=C=C").expect("allene");
    let mut bonds = Vec::new();
    let mut angles = Vec::new();

    collect_bonds_and_angles(&mol, &mut bonds, &mut angles);

    assert_eq!(bonds, vec![(0, 1), (1, 2)]);
    assert_eq!(angles, vec![vec![0, 1, 2, 1]]);
}

#[test]
fn init_bounds_ptr_sets_default_bounds_and_keeps_diagonal_zero() {
    let mut mmat = BoundsMatrix {
        data: vec![vec![0.0; 3]; 3],
        n: 3,
    };

    init_bounds_mat_ptr(&mut mmat, 0.25, 42.0);

    for i in 0..3 {
        assert_eq!(mmat.get_lower(i, i), 0.0);
        assert_eq!(mmat.get_upper(i, i), 0.0);
    }
    for i in 0..3 {
        for j in 0..3 {
            if i == j {
                continue;
            }
            assert_eq!(mmat.get_lower(i, j), 0.25);
            assert_eq!(mmat.get_upper(i, j), 42.0);
        }
    }
}

#[test]
fn init_bounds_shared_matches_pointer_overload() {
    let mut from_ptr = BoundsMatrix {
        data: vec![vec![0.0; 4]; 4],
        n: 4,
    };
    let mut from_shared = BoundsMatrix {
        data: vec![vec![0.0; 4]; 4],
        n: 4,
    };

    init_bounds_mat_ptr(&mut from_ptr, DEFAULT_LOWER, DEFAULT_UPPER);
    init_bounds_mat_shared(&mut from_shared, DEFAULT_LOWER, DEFAULT_UPPER);

    assert_eq!(from_ptr.data, from_shared.data);
}

#[test]
fn init_bounds_new_uses_rdkit_default_min_and_max() {
    let mmat = BoundsMatrix::new(2);

    assert_eq!(mmat.get_lower(0, 1), DEFAULT_LOWER);
    assert_eq!(mmat.get_upper(0, 1), DEFAULT_UPPER);
    assert_eq!(mmat.get_lower(0, 0), 0.0);
    assert_eq!(mmat.get_upper(1, 1), 0.0);
}

#[test]
fn bounds_matrix_stores_upper_in_upper_triangle_and_lower_in_lower_triangle() {
    let mut mmat = BoundsMatrix {
        data: vec![vec![0.0; 3]; 3],
        n: 3,
    };

    mmat.set_upper(2, 0, 4.2).expect("set upper");
    mmat.set_lower(2, 0, 1.8).expect("set lower");

    assert_eq!(mmat.data[0][2], 4.2);
    assert_eq!(mmat.data[2][0], 1.8);
    assert_eq!(mmat.get_upper(0, 2), 4.2);
    assert_eq!(mmat.get_upper(2, 0), 4.2);
    assert_eq!(mmat.get_lower(0, 2), 1.8);
    assert_eq!(mmat.get_lower(2, 0), 1.8);
}

#[test]
fn bounds_matrix_export_preserves_rdkit_raw_triangle_layout() {
    let mut mmat = BoundsMatrix::new(2);
    mmat.set_upper(0, 1, 2.5).expect("set upper");
    mmat.set_lower(0, 1, 1.25).expect("set lower");

    let raw = mmat.to_vec_vec();
    assert_eq!(raw, vec![vec![0.0, 2.5], vec![1.25, 0.0]]);
}

#[test]
fn bounds_matrix_set_upper_if_better_only_tightens_within_current_interval() {
    let mut mmat = BoundsMatrix::new(2);
    mmat.set_lower(0, 1, 1.2).expect("set lower");
    mmat.set_upper(0, 1, 3.4).expect("set upper");

    mmat.set_upper_if_better(0, 1, 2.6)
        .expect("set upper if better");
    assert_eq!(mmat.get_upper(0, 1), 2.6);

    mmat.set_upper_if_better(0, 1, 2.9)
        .expect("set upper if better");
    assert_eq!(mmat.get_upper(0, 1), 2.6);

    mmat.set_upper_if_better(0, 1, 1.0)
        .expect("set upper if better");
    assert_eq!(mmat.get_upper(0, 1), 2.6);
}

#[test]
fn bounds_matrix_set_lower_if_better_only_raises_within_current_interval() {
    let mut mmat = BoundsMatrix::new(2);
    mmat.set_lower(0, 1, 1.2).expect("set lower");
    mmat.set_upper(0, 1, 3.4).expect("set upper");

    mmat.set_lower_if_better(0, 1, 2.1)
        .expect("set lower if better");
    assert_eq!(mmat.get_lower(0, 1), 2.1);

    mmat.set_lower_if_better(0, 1, 1.8)
        .expect("set lower if better");
    assert_eq!(mmat.get_lower(0, 1), 2.1);

    mmat.set_lower_if_better(0, 1, 3.6)
        .expect("set lower if better");
    assert_eq!(mmat.get_lower(0, 1), 2.1);
}

#[test]
fn bounds_matrix_check_valid_detects_crossed_bounds() {
    let mut valid = BoundsMatrix::new(3);
    valid.set_lower(0, 2, 1.3).expect("set lower");
    valid.set_upper(0, 2, 2.8).expect("set upper");
    assert!(valid.check_valid());

    let mut invalid = BoundsMatrix::new(3);
    invalid.set_lower(0, 2, 2.9).expect("set lower");
    invalid.set_upper(0, 2, 2.1).expect("set upper");
    assert!(!invalid.check_valid());
}

#[test]
fn check_and_set_bounds_sets_uninitialized_pair_conservatively() {
    let mut mmat = BoundsMatrix::new(3);

    mmat.check_and_set_bounds(0, 2, 1.2, 2.8)
        .expect("check and set bounds");

    assert_eq!(mmat.get_lower(0, 2), 1.2);
    assert_eq!(mmat.get_upper(0, 2), 2.8);
    assert_eq!(mmat.data[2][0], 1.2);
    assert_eq!(mmat.data[0][2], 2.8);
}

#[test]
fn check_and_set_bounds_only_tightens_lower_and_relaxes_upper_when_allowed() {
    let mut mmat = BoundsMatrix::new(3);
    mmat.set_lower(0, 2, 1.8).expect("set lower");
    mmat.set_upper(0, 2, 2.2).expect("set upper");

    mmat.check_and_set_bounds(0, 2, 1.4, 2.6)
        .expect("check and set bounds");
    assert_eq!(mmat.get_lower(0, 2), 1.4);
    assert_eq!(mmat.get_upper(0, 2), 2.6);

    mmat.check_and_set_bounds(0, 2, 1.6, 2.0)
        .expect("check and set bounds");
    assert_eq!(mmat.get_lower(0, 2), 1.4);
    assert_eq!(mmat.get_upper(0, 2), 2.6);
}

#[test]
fn check_and_set_bounds_with_mode_uses_overlap_when_ranges_are_consistent() {
    let mut mmat = BoundsMatrix::new(2);
    mmat.set_lower(0, 1, 1.3).expect("set lower");
    mmat.set_upper(0, 1, 3.1).expect("set upper");

    mmat.check_and_set_bounds_with_mode(0, 1, 1.8, 2.4, true)
        .expect("check and set bounds");

    assert_eq!(mmat.get_lower(0, 1), 1.8);
    assert_eq!(mmat.get_upper(0, 1), 2.4);
    assert_eq!(mmat.data[1][0], 1.8);
    assert_eq!(mmat.data[0][1], 2.4);
}

#[test]
fn check_and_set_bounds_with_mode_falls_back_to_conservative_union_when_disjoint() {
    let mut mmat = BoundsMatrix::new(2);
    mmat.set_lower(0, 1, 1.3).expect("set lower");
    mmat.set_upper(0, 1, 2.0).expect("set upper");

    mmat.check_and_set_bounds_with_mode(0, 1, 2.4, 3.0, true)
        .expect("check and set bounds");

    assert_eq!(mmat.get_lower(0, 1), 1.3);
    assert_eq!(mmat.get_upper(0, 1), 3.0);
    assert_eq!(mmat.data[1][0], 1.3);
    assert_eq!(mmat.data[0][1], 3.0);
}

#[test]
fn set12_bounds_uses_uff_rest_length_for_supported_atoms() {
    let mol = Molecule::from_smiles_with_sanitize("CC", false).expect("ethane skeleton");

    let (mmat, accum_data) = run_set12_bounds(&mol);

    let lower = mmat.get_lower(0, 1);
    let upper = mmat.get_upper(0, 1);
    let width = upper - lower;
    assert!((width - (2.0 * DIST12_DELTA)).abs() < 1e-9);
    assert!(accum_data.bond_lengths[0] > 1.4);
    assert!(accum_data.bond_lengths[0] < 1.6);
    assert!((lower - (accum_data.bond_lengths[0] - DIST12_DELTA)).abs() < 1e-9);
    assert!((upper - (accum_data.bond_lengths[0] + DIST12_DELTA)).abs() < 1e-9);
    assert!(accum_data.visited12_bounds[1]);
}

#[test]
fn set12_bounds_falls_back_to_vdw_when_uff_params_are_missing() {
    let mut builder = MoleculeBuilder::new();
    let carbon = builder.add_atom(AtomSpec::new(Element::C));
    let dummy = builder.add_atom(AtomSpec::new(Element::DUMMY));
    builder
        .add_bond(BondSpec::new(carbon, dummy, BondOrder::Single))
        .unwrap();
    let mol = builder.build().expect("dummy-carbon molecule");

    let (mmat, accum_data) = run_set12_bounds(&mol);
    let expected = (vdw_radius(6) + vdw_radius(0)) / 2.0;

    assert!((accum_data.bond_lengths[0] - expected).abs() < 1e-9);
    assert!((mmat.get_lower(0, 1) - (0.5 * expected)).abs() < 1e-9);
    assert!((mmat.get_upper(0, 1) - (1.5 * expected)).abs() < 1e-9);
}

#[test]
fn set12_bounds_marks_visited_pid_using_sorted_atom_indices() {
    let mol = Molecule::from_smiles_with_sanitize("CCO", false).expect("ethanol skeleton");

    let (_mmat, accum_data) = run_set12_bounds(&mol);
    let n = mol.num_atoms();

    assert!(accum_data.visited12_bounds[1]);
    assert!(accum_data.visited12_bounds[n + 2]);
    assert!(!accum_data.visited12_bounds[2]);
}

#[test]
fn set12_bounds_adds_extra_squish_for_conjugated_hetero_five_ring_bonds() {
    let mol = Molecule::from_smiles("s1cccc1").expect("thiophene");

    let (mmat, _accum_data) = run_set12_bounds(&mol);

    let sulfur_atom = mol
        .atoms()
        .iter()
        .position(|atom| atom.atomic_number() == 16)
        .expect("sulfur atom");
    let squished_width = mol
        .bonds()
        .iter()
        .find(|bond| bond.begin().index() == sulfur_atom || bond.end().index() == sulfur_atom)
        .map(|bond| {
            mmat.get_upper(bond.begin().index(), bond.end().index())
                - mmat.get_lower(bond.begin().index(), bond.end().index())
        })
        .expect("sulfur bond");
    let mut sulfur_adjacent = vec![false; mol.num_atoms()];
    sulfur_adjacent[sulfur_atom] = true;
    for bond in mol.bonds() {
        if bond.begin().index() == sulfur_atom {
            sulfur_adjacent[bond.end().index()] = true;
        } else if bond.end().index() == sulfur_atom {
            sulfur_adjacent[bond.begin().index()] = true;
        }
    }
    let unsquished_width = mol
        .bonds()
        .iter()
        .find(|bond| !sulfur_adjacent[bond.begin().index()] && !sulfur_adjacent[bond.end().index()])
        .map(|bond| {
            mmat.get_upper(bond.begin().index(), bond.end().index())
                - mmat.get_lower(bond.begin().index(), bond.end().index())
        })
        .expect("carbon-carbon bond");

    assert!((squished_width - (2.0 * (0.2 + DIST12_DELTA))).abs() < 1e-9);
    assert!((unsquished_width - (2.0 * DIST12_DELTA)).abs() < 1e-9);
}

#[test]
fn hbond_helpers_follow_rdkit_acceptor_donor_and_bound_hydrogen_rules() {
    let ammonia = Molecule::from_smiles_with_sanitize("N", false).expect("ammonia-like N");
    let ammonia_valence = assign_valence(&ammonia, ValenceModel::RdkitLike).expect("valence");
    let ammonia_total_hs =
        total_num_hydrogens_for_distgeom(&ammonia, &ammonia.atoms()[0], &ammonia_valence, 0);
    assert!(is_hbond_acceptor(ammonia.atoms()[0].atomic_number()));
    assert!(is_hbond_donor(&ammonia.atoms()[0], ammonia_total_hs));

    let methane = Molecule::from_smiles_with_sanitize("C", false).expect("methane-like C");
    let methane_valence = assign_valence(&methane, ValenceModel::RdkitLike).expect("valence");
    let methane_total_hs =
        total_num_hydrogens_for_distgeom(&methane, &methane.atoms()[0], &methane_valence, 0);
    assert!(!is_hbond_acceptor(methane.atoms()[0].atomic_number()));
    assert!(!is_hbond_donor(&methane.atoms()[0], methane_total_hs));

    let mut builder = MoleculeBuilder::new();
    let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
    let oxygen = builder.add_atom(AtomSpec::new(Element::O));
    builder
        .add_bond(BondSpec::new(hydrogen, oxygen, BondOrder::Single))
        .unwrap();
    let hbond_mol = builder.build().expect("explicit O-H");
    assert!(is_h_in_hbond_donor(&hbond_mol, hydrogen.index()));
    assert!(!is_h_in_hbond_donor(&hbond_mol, oxygen.index()));
}

#[test]
fn total_num_hydrogens_for_distgeom_counts_neighbor_hydrogens_like_rdkit() {
    let mut builder = MoleculeBuilder::new();
    let hydrogen = builder.add_atom(AtomSpec::new(Element::H).with_no_implicit(true));
    let nitrogen = builder.add_atom(AtomSpec::new(Element::N).with_no_implicit(true));
    let methyl = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
    builder
        .add_bond(BondSpec::new(hydrogen, nitrogen, BondOrder::Single))
        .expect("h-n");
    builder
        .add_bond(BondSpec::new(nitrogen, methyl, BondOrder::Single))
        .expect("n-c");
    let mol = builder.build().expect("explicit secondary amine fragment");
    let assignment = assign_valence(&mol, ValenceModel::RdkitLike).expect("valence");
    let total_hs = total_num_hydrogens_for_distgeom(
        &mol,
        &mol.atoms()[nitrogen.index()],
        &assignment,
        nitrogen.index(),
    );
    assert_eq!(total_hs, 1);
}

#[test]
fn atom_charge_flags_adds_default_copper_charge_from_formal_charge() {
    let mol = single_atom_molecule(
        AtomSpec::new(Element::CU)
            .with_formal_charge(1)
            .with_no_implicit(true),
    );
    let assignment = assign_valence(&mol, ValenceModel::RdkitLike).expect("valence");
    let hybridizations = vec![Hybridization::Unspecified];
    let total_valence = atom_total_valence_for_uff(&assignment, 0);
    let mut atom_key = "Cu".to_string();

    add_atom_charge_flags_for_uff(
        &mol.atoms()[0],
        0,
        total_valence,
        &mut atom_key,
        hybridizations[0],
        false,
    );

    assert_eq!(atom_key, "Cu+1");
}

#[test]
fn atom_charge_flags_adds_lanthanide_plus_three_when_tolerating_mismatch() {
    let mol = single_atom_molecule(
        AtomSpec::new(Element::CE)
            .with_hybridization(Hybridization::Sp3d2)
            .with_no_implicit(true),
    );
    let assignment = assign_valence(&mol, ValenceModel::RdkitLike).expect("valence");
    let hybridizations = vec![Hybridization::Sp3d2];
    let total_valence = atom_total_valence_for_uff(&assignment, 0);
    let mut atom_key = "Ce".to_string();

    add_atom_charge_flags_for_uff(
        &mol.atoms()[0],
        0,
        total_valence,
        &mut atom_key,
        hybridizations[0],
        true,
    );

    assert_eq!(atom_key, "Ce+3");
}

#[test]
fn atom_charge_flags_rewrites_rhenium_special_labels_when_tolerating_mismatch() {
    let mol = single_atom_molecule(
        AtomSpec::new(Element::RE)
            .with_hybridization(Hybridization::Sp3d)
            .with_no_implicit(true),
    );
    let assignment = assign_valence(&mol, ValenceModel::RdkitLike).expect("valence");
    let hybridizations = vec![Hybridization::Sp3d];
    let total_valence = atom_total_valence_for_uff(&assignment, 0);
    let mut atom_key = "Re6".to_string();

    add_atom_charge_flags_for_uff(
        &mol.atoms()[0],
        0,
        total_valence,
        &mut atom_key,
        hybridizations[0],
        true,
    );

    assert_eq!(atom_key, "Re6+5");
}

#[test]
fn atom_label_uses_sp2_suffix_for_nonaromatic_alkene_carbon() {
    let mol = Molecule::from_smiles("C=C").expect("ethene");

    let label = uff_atom_label(&mol, 0);

    assert_eq!(label, "C_2");
}

#[test]
fn atom_label_uses_aromatic_r_suffix_for_benzene_carbon() {
    let mol = Molecule::from_smiles("c1ccccc1").expect("benzene");

    let label = uff_atom_label(&mol, 0);

    assert_eq!(label, "C_R");
}

#[test]
fn atom_label_composes_charge_flags_after_base_label() {
    let mol = single_atom_molecule(
        AtomSpec::new(Element::CU)
            .with_formal_charge(1)
            .with_no_implicit(true),
    );

    let label = uff_atom_label(&mol, 0);

    assert!(label.starts_with("Cu"));
    assert!(label.ends_with("+1"));
}

#[test]
fn set_lower_bound_vdw_scales_15_16_and_longer_paths_like_rdkit() {
    let mut builder = MoleculeBuilder::new();
    for _ in 0..7 {
        builder.add_atom(AtomSpec::new(Element::C));
    }
    let mol = builder.build().expect("seven carbons");
    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let mut dmat = vec![0.0; mol.num_atoms() * mol.num_atoms()];
    dmat[4 * mol.num_atoms()] = 4.0;
    dmat[5 * mol.num_atoms()] = 5.0;
    dmat[6 * mol.num_atoms()] = 6.0;

    set_lower_bound_vdw(&mol, &mut mmat, true, &dmat).expect("setLowerBoundVDW");

    let vdw_sum = vdw_radius(6) + vdw_radius(6);
    assert!((mmat.get_lower(4, 0) - (VDW_SCALE_15 * vdw_sum)).abs() < 1e-9);
    assert!(
        (mmat.get_lower(5, 0) - ((VDW_SCALE_15 + 0.5 * (1.0 - VDW_SCALE_15)) * vdw_sum)).abs()
            < 1e-9
    );
    assert!((mmat.get_lower(6, 0) - vdw_sum).abs() < 1e-9);
}

#[test]
fn set_lower_bound_vdw_uses_hbond_floor_before_vdw_scaling() {
    let mut builder = MoleculeBuilder::new();
    let nitrogen = builder.add_atom(AtomSpec::new(Element::N));
    let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
    let _oxygen = builder.add_atom(AtomSpec::new(Element::O));
    builder
        .add_bond(BondSpec::new(nitrogen, hydrogen, BondOrder::Single))
        .unwrap();
    let mol = builder.build().expect("H-N ... O");
    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let mut dmat = vec![0.0; mol.num_atoms() * mol.num_atoms()];
    dmat[2 * mol.num_atoms() + 1] = 6.0;
    dmat[1 * mol.num_atoms() + 2] = 6.0;

    set_lower_bound_vdw(&mol, &mut mmat, true, &dmat).expect("setLowerBoundVDW");

    assert!((mmat.get_lower(2, 1) - H_BOND_LENGTH).abs() < 1e-9);
}

#[test]
fn set_ring_angle_matches_rdkit_ring_hybridization_special_cases() {
    let mut three_builder = MoleculeBuilder::new();
    three_builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2));
    let cyclopropane_like = three_builder.build().expect("3-ring test atom");

    let mut five_builder = MoleculeBuilder::new();
    five_builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3));
    let cyclopentane_like = five_builder.build().expect("5-ring test atom");

    let tri = set_ring_angle(&cyclopropane_like, 0, 3);
    let five = set_ring_angle(&cyclopentane_like, 0, 5);

    assert!((tri - std::f64::consts::PI / 3.0).abs() < 1e-9);
    assert!((five - (104.0_f64.to_radians())).abs() < 1e-9);
}

#[test]
fn compute13_dist_returns_bond_sum_for_linear_angle() {
    let dist = compute_13_dist(1.4, 1.5, std::f64::consts::PI);

    assert!((dist - 2.9).abs() < 1e-9);
}

#[test]
fn compute13_dist_returns_bond_difference_for_zero_angle() {
    let dist = compute_13_dist(1.5, 1.4, 0.0);

    assert!((dist - 0.1).abs() < 1e-9);
}

#[test]
fn set_13_bounds_helper_doubles_tolerance_for_larger_sp2_ring_atoms() {
    let thiophene = Molecule::from_smiles("s1cccc1").expect("thiophene");
    let benzene = Molecule::from_smiles("c1ccccc1").expect("benzene");

    let (_thiophene_mmat_12, thiophene_accum) = run_set12_bounds(&thiophene);
    let (_benzene_mmat_12, benzene_accum) = run_set12_bounds(&benzene);
    let thiophene_rings = thiophene.derived_cache().rings.as_ref().expect("rings");
    let benzene_rings = benzene.derived_cache().rings.as_ref().expect("rings");

    let sulfur = thiophene
        .atoms()
        .iter()
        .position(|atom| atom.atomic_number() == 16)
        .expect("sulfur");
    let sulfur_neighbors: Vec<usize> = thiophene
        .bonds()
        .iter()
        .filter_map(|bond| {
            if bond.begin().index() == sulfur {
                Some(bond.end().index())
            } else if bond.end().index() == sulfur {
                Some(bond.begin().index())
            } else {
                None
            }
        })
        .collect();
    let mut thiophene_mmat = BoundsMatrix::new(thiophene.num_atoms());
    let thiophene_angle = set_ring_angle(&thiophene, sulfur, 5);
    set_13_bounds_helper(
        sulfur_neighbors[0],
        sulfur,
        sulfur_neighbors[1],
        thiophene_angle,
        &thiophene_accum.bond_lengths,
        &mut thiophene_mmat,
        &thiophene,
        thiophene_rings,
    )
    .expect("thiophene 1-3 bounds");

    let mut benzene_mmat = BoundsMatrix::new(benzene.num_atoms());
    let benzene_angle = set_ring_angle(&benzene, 1, 6);
    set_13_bounds_helper(
        0,
        1,
        2,
        benzene_angle,
        &benzene_accum.bond_lengths,
        &mut benzene_mmat,
        &benzene,
        benzene_rings,
    )
    .expect("benzene 1-3 bounds");

    let thiophene_width = thiophene_mmat.get_upper(sulfur_neighbors[0], sulfur_neighbors[1])
        - thiophene_mmat.get_lower(sulfur_neighbors[0], sulfur_neighbors[1]);
    let benzene_width = benzene_mmat.get_upper(0, 2) - benzene_mmat.get_lower(0, 2);

    assert!((thiophene_width - (4.0 * DIST13_TOL)).abs() < 1e-9);
    assert!((benzene_width - (2.0 * DIST13_TOL)).abs() < 1e-9);
    assert!(is_larger_sp2_atom_idx(&thiophene, thiophene_rings, sulfur));
    assert!(!is_larger_sp2_atom_idx(&benzene, benzene_rings, 1));
}

#[test]
fn visited_bound_obeys_rdkit_dist_type_thresholds() {
    let mut accum = ComputedData::new(4, 2);
    accum.visited13_bounds[3] = true;
    accum.visited14_bounds[5] = true;

    assert!(!accum.visited_bound(3, DistType::Dist12));
    assert!(accum.visited_bound(3, DistType::Dist13));
    assert!(accum.visited_bound(3, DistType::Dist14));
    assert!(!accum.visited_bound(5, DistType::Dist13));
    assert!(accum.visited_bound(5, DistType::Dist14));
}

#[test]
fn set_13_bounds_sets_non_ring_sp3_bounds_for_propane_path() {
    let mol = Molecule::from_smiles("CCC").expect("propane");
    let (mmat, accum_data) = run_set13_bounds(&mol);

    let bid01 = bond_between_idx_simple(&mol, 0, 1).expect("0-1 bond");
    let bid12 = bond_between_idx_simple(&mol, 1, 2).expect("1-2 bond");
    let expected = compute_13_dist(
        accum_data.bond_lengths[bid01],
        accum_data.bond_lengths[bid12],
        109.5_f64.to_radians(),
    );

    assert!((mmat.get_lower(0, 2) - (expected - DIST13_TOL)).abs() < 1e-9);
    assert!((mmat.get_upper(0, 2) - (expected + DIST13_TOL)).abs() < 1e-9);
    assert!(
        (accum_data.get_bond_angle(mol.num_bonds(), bid01, bid12) - 109.5_f64.to_radians()).abs()
            < 1e-9
    );
    assert!(accum_data.visited13_bounds[2]);
}

#[test]
fn set_13_bounds_distributes_remaining_fused_ring_angle_like_rdkit() {
    let mol = Molecule::from_smiles("c1cccc2ccccc12").expect("naphthalene");
    let (_mmat, accum_data) = run_set13_bounds(&mol);
    let rings = mol.derived_cache().rings.as_ref().expect("rings");

    let fusion_atom = mol
        .atoms()
        .iter()
        .enumerate()
        .find_map(|(idx, atom)| {
            (atom.hybridization() == Hybridization::Sp2
                && rings.num_atom_rings(atom.id()) > 1
                && neighbors_for_atom(&mol, idx).len() == 3)
                .then_some(idx)
        })
        .expect("fusion atom");

    let neighbors = neighbors_for_atom(&mol, fusion_atom);
    let mut pair_angles = Vec::new();
    for left in 0..neighbors.len() {
        let bid1 = bond_between_idx_simple(&mol, fusion_atom, neighbors[left]).expect("bond 1");
        for right in 0..left {
            let bid2 =
                bond_between_idx_simple(&mol, fusion_atom, neighbors[right]).expect("bond 2");
            pair_angles.push(accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2));
        }
    }

    assert_eq!(pair_angles.len(), 3);
    for angle in pair_angles {
        assert!((angle - 120.0_f64.to_radians()).abs() < 1e-9);
    }
}

#[test]
fn set_13_bounds_uses_wide_bounds_for_non_ring_degree_five_center() {
    let mol = Molecule::from_smiles("FP(F)(F)(F)F").expect("PF5-like");
    let center = mol
        .atoms()
        .iter()
        .position(|atom| atom.atomic_number() == 15)
        .expect("phosphorus center");
    let ligands: Vec<usize> = neighbors_for_atom(&mol, center);
    assert_eq!(ligands.len(), 5);

    let (mmat, accum_data) = run_set13_bounds(&mol);
    let bid1 = bond_between_idx_simple(&mol, center, ligands[0]).expect("P-F1");
    let bid2 = bond_between_idx_simple(&mol, center, ligands[1]).expect("P-F2");
    let dmax = accum_data.bond_lengths[bid1] + accum_data.bond_lengths[bid2];

    assert!((mmat.get_lower(ligands[0], ligands[1]) - 1.0).abs() < 1e-9);
    assert!((mmat.get_upper(ligands[0], ligands[1]) - (dmax * 1.2)).abs() < 1e-9);
    assert!(accum_data.visited13_bounds[ligands[0] * mol.num_atoms() + ligands[1]]);
}

#[test]
fn chain_and_carbonyl_classification_helpers_follow_rdkit_patterns() {
    let propane = Molecule::from_smiles("CCC").expect("propane");
    assert!(check_h2_nx3_h1_ox2(&propane, 1));

    let ether = Molecule::from_smiles("COC").expect("dimethyl ether");
    let oxygen = ether
        .atoms()
        .iter()
        .position(|atom| atom.atomic_number() == 8)
        .expect("ether oxygen");
    assert!(check_h2_nx3_h1_ox2(&ether, oxygen));

    let butane = Molecule::from_smiles("CCCC").expect("butane");
    assert!(check_nh_ch_ch_nh(&butane, 0, 1, 2, 3));

    let acetate = Molecule::from_smiles("CC(=O)O").expect("acetate");
    let carbonyl = acetate
        .atoms()
        .iter()
        .enumerate()
        .find_map(|(idx, _)| is_carbonyl(&acetate, idx).then_some(idx))
        .expect("carbonyl carbon");
    assert!(is_carbonyl(&acetate, carbonyl));
    assert!(!is_carbonyl(&acetate, 0));
}

#[test]
fn amide_ester_classification_helpers_match_ester_patterns() {
    let ester = Molecule::from_smiles("COC(=O)C").expect("methyl acetate");
    let carbonyl = ester
        .atoms()
        .iter()
        .enumerate()
        .find_map(|(idx, _)| is_carbonyl(&ester, idx).then_some(idx))
        .expect("carbonyl carbon");

    let double_hetero = neighbors_for_atom(&ester, carbonyl)
        .into_iter()
        .find(|&nbr| {
            let bond_idx = bond_between_idx_simple(&ester, carbonyl, nbr).expect("bond");
            ester.bonds()[bond_idx].order() == BondOrder::Double
                && (ester.atoms()[nbr].atomic_number() == 8
                    || ester.atoms()[nbr].atomic_number() == 7)
        })
        .expect("double bonded hetero");
    let single_hetero = neighbors_for_atom(&ester, carbonyl)
        .into_iter()
        .find(|&nbr| {
            let bond_idx = bond_between_idx_simple(&ester, carbonyl, nbr).expect("bond");
            ester.bonds()[bond_idx].order() == BondOrder::Single
                && (ester.atoms()[nbr].atomic_number() == 8
                    || ester.atoms()[nbr].atomic_number() == 7)
        })
        .expect("single bonded hetero");
    let atom1 = neighbors_for_atom(&ester, single_hetero)
        .into_iter()
        .find(|&nbr| nbr != carbonyl)
        .expect("atom1");
    let bnd1 = bond_between_idx_simple(&ester, atom1, single_hetero).expect("bond1");
    let bnd3_double =
        bond_between_idx_simple(&ester, carbonyl, double_hetero).expect("bond3 double");
    let carbonyl_substituent = neighbors_for_atom(&ester, carbonyl)
        .into_iter()
        .find(|&nbr| nbr != single_hetero && nbr != double_hetero)
        .expect("carbonyl substituent");
    let bnd3_single =
        bond_between_idx_simple(&ester, carbonyl, carbonyl_substituent).expect("bond3 single");

    assert!(check_amide_ester_14(
        &ester,
        bnd1,
        bnd3_double,
        single_hetero,
        carbonyl,
        double_hetero,
    ));
    assert!(check_amide_ester_15(
        &ester,
        bnd1,
        bnd3_single,
        single_hetero,
        carbonyl,
    ));

    let tertiary_amide = Molecule::from_smiles("CN(C)C(=O)C").expect("tertiary amide");
    let tertiary_carbonyl = tertiary_amide
        .atoms()
        .iter()
        .enumerate()
        .find_map(|(idx, _)| is_carbonyl(&tertiary_amide, idx).then_some(idx))
        .expect("tertiary amide carbonyl");
    let tertiary_nitrogen = neighbors_for_atom(&tertiary_amide, tertiary_carbonyl)
        .into_iter()
        .find(|&nbr| {
            let bond_idx =
                bond_between_idx_simple(&tertiary_amide, tertiary_carbonyl, nbr).expect("bond");
            tertiary_amide.bonds()[bond_idx].order() == BondOrder::Single
                && tertiary_amide.atoms()[nbr].atomic_number() == 7
        })
        .expect("amide nitrogen");
    let tertiary_atom1 = neighbors_for_atom(&tertiary_amide, tertiary_nitrogen)
        .into_iter()
        .find(|&nbr| nbr != tertiary_carbonyl)
        .expect("substituent carbon");
    let tertiary_bnd1 = bond_between_idx_simple(&tertiary_amide, tertiary_atom1, tertiary_nitrogen)
        .expect("tertiary bond1");
    let tertiary_side = neighbors_for_atom(&tertiary_amide, tertiary_carbonyl)
        .into_iter()
        .find(|&nbr| {
            let bond_idx =
                bond_between_idx_simple(&tertiary_amide, tertiary_carbonyl, nbr).expect("bond");
            tertiary_amide.bonds()[bond_idx].order() == BondOrder::Single
                && nbr != tertiary_nitrogen
        })
        .expect("carbonyl side");
    let tertiary_bnd3 = bond_between_idx_simple(&tertiary_amide, tertiary_carbonyl, tertiary_side)
        .expect("tertiary bond3");

    assert!(!check_amide_ester_15(
        &tertiary_amide,
        tertiary_bnd1,
        tertiary_bnd3,
        tertiary_nitrogen,
        tertiary_carbonyl,
    ));
}

#[test]
fn macrocycle_all_in_same_ring_amide_ester_helper_matches_lactone_pattern() {
    let mol = Molecule::from_smiles("O=C1N(C)CCCCCCCC1").expect("macrocyclic tertiary lactam");
    let carbonyl = mol
        .atoms()
        .iter()
        .enumerate()
        .find_map(|(idx, _)| is_carbonyl(&mol, idx).then_some(idx))
        .expect("carbonyl carbon");
    let atm2 = neighbors_for_atom(&mol, carbonyl)
        .into_iter()
        .find(|&nbr| {
            let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
            mol.bonds()[bond_idx].order() == BondOrder::Single
                && mol.atoms()[nbr].atomic_number() == 7
        })
        .expect("ring amide nitrogen");
    let atm4 = neighbors_for_atom(&mol, carbonyl)
        .into_iter()
        .find(|&nbr| {
            let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
            mol.bonds()[bond_idx].order() == BondOrder::Single && nbr != atm2
        })
        .expect("ring carbon neighbor");
    let atm1 = neighbors_for_atom(&mol, atm2)
        .into_iter()
        .find(|&nbr| nbr != carbonyl && mol.atoms()[nbr].atomic_number() == 6)
        .expect("preceding ring atom");

    assert!(check_macrocycle_all_in_same_ring_amide_ester_14(
        &mol, atm1, atm2, carbonyl, atm4,
    ));
}

#[test]
fn macrocycle_two_in_same_ring_amide_ester_helper_matches_tertiary_lactam_pattern() {
    let mol = Molecule::from_smiles("O=C1N(C)CCCCCCCC1").expect("macrocyclic tertiary lactam");
    let carbonyl = mol
        .atoms()
        .iter()
        .enumerate()
        .find_map(|(idx, _)| is_carbonyl(&mol, idx).then_some(idx))
        .expect("carbonyl carbon");
    let oxygen = neighbors_for_atom(&mol, carbonyl)
        .into_iter()
        .find(|&nbr| {
            let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
            mol.bonds()[bond_idx].order() == BondOrder::Double
        })
        .expect("carbonyl oxygen");
    let nitrogen = neighbors_for_atom(&mol, carbonyl)
        .into_iter()
        .find(|&nbr| {
            let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
            mol.bonds()[bond_idx].order() == BondOrder::Single
                && mol.atoms()[nbr].atomic_number() == 7
        })
        .expect("amide nitrogen");
    let atom1 = neighbors_for_atom(&mol, nitrogen)
        .into_iter()
        .find(|&nbr| nbr != carbonyl && mol.atoms()[nbr].atomic_number() == 6)
        .expect("preceding ring carbon");
    let bnd1 = bond_between_idx_simple(&mol, atom1, nitrogen).expect("bond1");
    let bnd3 = bond_between_idx_simple(&mol, carbonyl, oxygen).expect("bond3");

    assert!(check_macrocycle_two_in_same_ring_amide_ester_14(
        &mol, bnd1, bnd3, atom1, nitrogen, carbonyl, oxygen,
    ));
}

#[test]
fn set_macrocycle_two_in_same_ring_14_bounds_uses_cis_for_macrocycle_amide_path() {
    let mol = Molecule::from_smiles("O=C1N(C)CCCCCCCC1").expect("macrocyclic tertiary lactam");
    let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
    let dmat = flatten_topological_distances(&mol);

    let carbonyl = mol
        .atoms()
        .iter()
        .enumerate()
        .find_map(|(idx, _)| is_carbonyl(&mol, idx).then_some(idx))
        .expect("carbonyl carbon");
    let oxygen = neighbors_for_atom(&mol, carbonyl)
        .into_iter()
        .find(|&nbr| {
            let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
            mol.bonds()[bond_idx].order() == BondOrder::Double
        })
        .expect("carbonyl oxygen");
    let nitrogen = neighbors_for_atom(&mol, carbonyl)
        .into_iter()
        .find(|&nbr| {
            let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
            mol.bonds()[bond_idx].order() == BondOrder::Single
                && mol.atoms()[nbr].atomic_number() == 7
        })
        .expect("amide nitrogen");
    let atom1 = neighbors_for_atom(&mol, nitrogen)
        .into_iter()
        .find(|&nbr| nbr != carbonyl && mol.atoms()[nbr].atomic_number() == 6)
        .expect("preceding ring carbon");

    let bid1 = bond_between_idx_simple(&mol, atom1, nitrogen).expect("b1");
    let bid2 = bond_between_idx_simple(&mol, nitrogen, carbonyl).expect("b2");
    let bid3 = bond_between_idx_simple(&mol, carbonyl, oxygen).expect("b3");
    set_macrocycle_two_in_same_ring_14_bounds(
        &mol,
        bid1,
        bid2,
        bid3,
        &mut accum_data,
        &mut mmat,
        &dmat,
    )
    .expect("macrocycle two-in-same-ring bounds");

    let path = accum_data.paths14.last().expect("path");
    assert_eq!(path.kind, Path14Kind::Cis);
    let expected = compute_14_dist_cis(
        accum_data.bond_lengths[bid1],
        accum_data.bond_lengths[bid2],
        accum_data.bond_lengths[bid3],
        accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2),
        accum_data.get_bond_angle(mol.num_bonds(), bid2, bid3),
    );
    assert!((mmat.get_lower(atom1, oxygen) - (expected - GEN_DIST_TOL)).abs() < 1e-9);
    assert!((mmat.get_upper(atom1, oxygen) - (expected + GEN_DIST_TOL)).abs() < 1e-9);
}

#[test]
fn set_macrocycle_all_in_same_ring_14_bounds_uses_trans_plus_point_one_for_macrocycle_amide() {
    let mol = Molecule::from_smiles("O=C1N(C)CCCCCCCC1").expect("macrocyclic tertiary lactam");
    let (mut mmat, mut accum_data) = run_set13_bounds(&mol);

    let carbonyl = mol
        .atoms()
        .iter()
        .enumerate()
        .find_map(|(idx, _)| is_carbonyl(&mol, idx).then_some(idx))
        .expect("carbonyl carbon");
    let nitrogen = neighbors_for_atom(&mol, carbonyl)
        .into_iter()
        .find(|&nbr| {
            let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
            mol.bonds()[bond_idx].order() == BondOrder::Single
                && mol.atoms()[nbr].atomic_number() == 7
        })
        .expect("amide nitrogen");
    let atm4 = neighbors_for_atom(&mol, carbonyl)
        .into_iter()
        .find(|&nbr| {
            let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
            mol.bonds()[bond_idx].order() == BondOrder::Single && nbr != nitrogen
        })
        .expect("ring carbon neighbor");
    let atom1 = neighbors_for_atom(&mol, nitrogen)
        .into_iter()
        .find(|&nbr| nbr != carbonyl && mol.atoms()[nbr].atomic_number() == 6)
        .expect("preceding ring carbon");

    let bid1 = bond_between_idx_simple(&mol, atom1, nitrogen).expect("b1");
    let bid2 = bond_between_idx_simple(&mol, nitrogen, carbonyl).expect("b2");
    let bid3 = bond_between_idx_simple(&mol, carbonyl, atm4).expect("b3");
    set_macrocycle_all_in_same_ring_14_bounds(&mol, bid1, bid2, bid3, &mut accum_data, &mut mmat)
        .expect("macrocycle all-in-same-ring bounds");

    let path = accum_data.paths14.last().expect("path");
    assert_eq!(path.kind, Path14Kind::Trans);
    let expected = compute_14_dist_trans(
        accum_data.bond_lengths[bid1],
        accum_data.bond_lengths[bid2],
        accum_data.bond_lengths[bid3],
        accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2),
        accum_data.get_bond_angle(mol.num_bonds(), bid2, bid3),
    ) + 0.1;
    assert!((mmat.get_lower(atom1, atm4) - (expected - GEN_DIST_TOL)).abs() < 1e-9);
    assert!((mmat.get_upper(atom1, atm4) - (expected + GEN_DIST_TOL)).abs() < 1e-9);
}

#[test]
fn set_macrocycle_all_in_same_ring_14_bounds_uses_other_for_plain_macrocycle_chain() {
    let mol = Molecule::from_smiles("C1CCCCCCCCC1").expect("cyclodecane");
    let (mut mmat, mut accum_data, _) = run_set14_same_ring_pass_only(&mol, true);
    let (bid1, bid2, bid3, _) =
        find_same_ring_dispatch_triple(&mol, true).expect("macrocycle path");

    let before_paths = accum_data.paths14.len();
    set_macrocycle_all_in_same_ring_14_bounds(&mol, bid1, bid2, bid3, &mut accum_data, &mut mmat)
        .expect("macrocycle all-in-same-ring bounds");

    let path = accum_data.paths14.get(before_paths).expect("new path");
    assert_eq!(path.kind, Path14Kind::Other);
}

#[test]
fn set_chain_14_bounds_uses_defined_double_bond_stereo_for_alkenes() {
    let trans = Molecule::from_smiles("C/C=C/C").expect("trans alkene");
    let (mut trans_mmat, mut trans_accum) = run_set13_bounds(&trans);
    let t_bid1 = bond_between_idx_simple(&trans, 0, 1).expect("0-1");
    let t_bid2 = bond_between_idx_simple(&trans, 1, 2).expect("1-2");
    let t_bid3 = bond_between_idx_simple(&trans, 2, 3).expect("2-3");
    set_chain_14_bounds(
        &trans,
        t_bid1,
        t_bid2,
        t_bid3,
        &mut trans_accum,
        &mut trans_mmat,
        false,
    )
    .expect("chain trans bounds");
    let trans_path = trans_accum.paths14.last().expect("trans path");
    assert_eq!(trans_path.kind, Path14Kind::Trans);
    let trans_expected = compute_14_dist_trans(
        trans_accum.bond_lengths[t_bid1],
        trans_accum.bond_lengths[t_bid2],
        trans_accum.bond_lengths[t_bid3],
        trans_accum.get_bond_angle(trans.num_bonds(), t_bid1, t_bid2),
        trans_accum.get_bond_angle(trans.num_bonds(), t_bid2, t_bid3),
    );
    assert!((trans_mmat.get_lower(0, 3) - (trans_expected - GEN_DIST_TOL)).abs() < 1e-9);
    assert!((trans_mmat.get_upper(0, 3) - (trans_expected + GEN_DIST_TOL)).abs() < 1e-9);

    let cis = Molecule::from_smiles("C/C=C\\C").expect("cis alkene");
    let (mut cis_mmat, mut cis_accum) = run_set13_bounds(&cis);
    let c_bid1 = bond_between_idx_simple(&cis, 0, 1).expect("0-1");
    let c_bid2 = bond_between_idx_simple(&cis, 1, 2).expect("1-2");
    let c_bid3 = bond_between_idx_simple(&cis, 2, 3).expect("2-3");
    set_chain_14_bounds(
        &cis,
        c_bid1,
        c_bid2,
        c_bid3,
        &mut cis_accum,
        &mut cis_mmat,
        false,
    )
    .expect("chain cis bounds");
    let cis_path = cis_accum.paths14.last().expect("cis path");
    assert_eq!(cis_path.kind, Path14Kind::Cis);
    let cis_expected = compute_14_dist_cis(
        cis_accum.bond_lengths[c_bid1],
        cis_accum.bond_lengths[c_bid2],
        cis_accum.bond_lengths[c_bid3],
        cis_accum.get_bond_angle(cis.num_bonds(), c_bid1, c_bid2),
        cis_accum.get_bond_angle(cis.num_bonds(), c_bid2, c_bid3),
    );
    assert!((cis_mmat.get_lower(0, 3) - (cis_expected - GEN_DIST_TOL)).abs() < 1e-9);
    assert!((cis_mmat.get_upper(0, 3) - (cis_expected + GEN_DIST_TOL)).abs() < 1e-9);
}

#[test]
fn set_chain_14_bounds_uses_ss_special_case() {
    let mol = Molecule::from_smiles("CSSC").expect("disulfide");
    let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
    let bid1 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
    let bid2 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
    let bid3 = bond_between_idx_simple(&mol, 2, 3).expect("2-3");
    set_chain_14_bounds(&mol, bid1, bid2, bid3, &mut accum_data, &mut mmat, false)
        .expect("chain S-S bounds");

    let path = accum_data.paths14.last().expect("path");
    assert_eq!(path.kind, Path14Kind::Other);
    let expected = compute_14_dist_3d(
        accum_data.bond_lengths[bid1],
        accum_data.bond_lengths[bid2],
        accum_data.bond_lengths[bid3],
        accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2),
        accum_data.get_bond_angle(mol.num_bonds(), bid2, bid3),
        std::f64::consts::PI / 2.0,
    );
    assert!((mmat.get_lower(0, 3) - (expected - GEN_DIST_TOL)).abs() < 1e-9);
    assert!((mmat.get_upper(0, 3) - (expected + GEN_DIST_TOL)).abs() < 1e-9);
}

#[test]
fn set_chain_14_bounds_honors_force_trans_amides_for_secondary_amide_h_paths() {
    let mut builder = Molecule::builder();
    let hydrogen = builder.add_atom(AtomSpec::new(Element::H)).index();
    let nitrogen =
        builder.add_atom(AtomSpec::new(Element::N).with_hybridization(Hybridization::Sp2));
    let nitrogen = nitrogen.index();
    let carbonyl =
        builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2));
    let carbonyl = carbonyl.index();
    let oxygen = builder.add_atom(AtomSpec::new(Element::O).with_hybridization(Hybridization::Sp2));
    let oxygen = oxygen.index();
    let n_methyl = builder.add_atom(AtomSpec::new(Element::C)).index();
    let carbonyl_methyl = builder.add_atom(AtomSpec::new(Element::C)).index();
    builder
        .add_bond(BondSpec::new(
            AtomId::new(hydrogen),
            AtomId::new(nitrogen),
            BondOrder::Single,
        ))
        .expect("h-n");
    builder
        .add_bond(BondSpec::new(
            AtomId::new(nitrogen),
            AtomId::new(carbonyl),
            BondOrder::Single,
        ))
        .expect("n-c");
    builder
        .add_bond(BondSpec::new(
            AtomId::new(nitrogen),
            AtomId::new(n_methyl),
            BondOrder::Single,
        ))
        .expect("n-c methyl");
    builder
        .add_bond(BondSpec::new(
            AtomId::new(carbonyl),
            AtomId::new(oxygen),
            BondOrder::Double,
        ))
        .expect("c-o");
    builder
        .add_bond(BondSpec::new(
            AtomId::new(carbonyl),
            AtomId::new(carbonyl_methyl),
            BondOrder::Single,
        ))
        .expect("c-c");
    let mol = builder.build().expect("secondary amide with explicit N-H");
    let (mut mmat, mut accum_data) = run_set12_bounds(&mol);

    let bid_hn = bond_between_idx_simple(&mol, hydrogen, nitrogen).expect("h-n");
    let bid_nc = bond_between_idx_simple(&mol, nitrogen, carbonyl).expect("n-c");
    let bid_co = bond_between_idx_simple(&mol, carbonyl, oxygen).expect("c-o");
    let bid_cm = bond_between_idx_simple(&mol, carbonyl, carbonyl_methyl).expect("c-m");
    let nb = mol.num_bonds();
    accum_data.set_bond_adj(nb, bid_hn, bid_nc, nitrogen as i32);
    accum_data.set_bond_adj(nb, bid_nc, bid_co, carbonyl as i32);
    accum_data.set_bond_adj(nb, bid_nc, bid_cm, carbonyl as i32);
    accum_data.set_bond_angle(
        nb,
        bid_hn,
        bid_nc,
        ideal_bond_angle(&mol.atoms()[nitrogen].hybridization(), None),
    );
    accum_data.set_bond_angle(
        nb,
        bid_nc,
        bid_co,
        ideal_bond_angle(&mol.atoms()[carbonyl].hybridization(), None),
    );
    accum_data.set_bond_angle(
        nb,
        bid_nc,
        bid_cm,
        ideal_bond_angle(&mol.atoms()[carbonyl].hybridization(), None),
    );

    set_chain_14_bounds(
        &mol,
        bid_hn,
        bid_nc,
        bid_co,
        &mut accum_data,
        &mut mmat,
        true,
    )
    .expect("chain amide14 bounds");
    let amide14_path = accum_data.paths14.last().expect("amide14 path");
    assert_eq!(amide14_path.kind, Path14Kind::Trans);
    let expected_14 = compute_14_dist_trans(
        accum_data.bond_lengths[bid_hn],
        accum_data.bond_lengths[bid_nc],
        accum_data.bond_lengths[bid_co],
        accum_data.get_bond_angle(mol.num_bonds(), bid_hn, bid_nc),
        accum_data.get_bond_angle(mol.num_bonds(), bid_nc, bid_co),
    );
    assert!((mmat.get_lower(hydrogen, oxygen) - (expected_14 - GEN_DIST_TOL)).abs() < 1e-9);
    assert!((mmat.get_upper(hydrogen, oxygen) - (expected_14 + GEN_DIST_TOL)).abs() < 1e-9);

    set_chain_14_bounds(
        &mol,
        bid_hn,
        bid_nc,
        bid_cm,
        &mut accum_data,
        &mut mmat,
        true,
    )
    .expect("chain amide15 bounds");
    let amide15_path = accum_data.paths14.last().expect("amide15 path");
    assert_eq!(amide15_path.kind, Path14Kind::Cis);
    let expected_15 = compute_14_dist_cis(
        accum_data.bond_lengths[bid_hn],
        accum_data.bond_lengths[bid_nc],
        accum_data.bond_lengths[bid_cm],
        accum_data.get_bond_angle(mol.num_bonds(), bid_hn, bid_nc),
        accum_data.get_bond_angle(mol.num_bonds(), bid_nc, bid_cm),
    );
    assert!(
        (mmat.get_lower(hydrogen, carbonyl_methyl) - (expected_15 - GEN_DIST_TOL)).abs() < 1e-9
    );
    assert!(
        (mmat.get_upper(hydrogen, carbonyl_methyl) - (expected_15 + GEN_DIST_TOL)).abs() < 1e-9
    );
}

#[test]
fn record_14_path_marks_sp2_sp2_ring_paths_as_cis() {
    let mol = Molecule::from_smiles("c1ccccc1").expect("benzene");
    let (_mmat, mut accum_data) = run_set13_bounds(&mol);
    let bid1 = bond_between_idx_simple(&mol, 5, 0).expect("5-0");
    let bid2 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
    let bid3 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
    record_14_path(&mol, bid1, bid2, bid3, &mut accum_data).expect("record14Path");

    let path = accum_data.paths14.last().expect("path");
    assert_eq!(path.bid1, bid1);
    assert_eq!(path.bid2, bid2);
    assert_eq!(path.bid3, bid3);
    assert_eq!(path.kind, Path14Kind::Cis);
    assert!(has_path_flag(
        &accum_data.cis_paths,
        path14_id(mol.num_bonds(), bid1, bid2, bid3)
    ));
    assert!(has_path_flag(
        &accum_data.cis_paths,
        path14_id(mol.num_bonds(), bid3, bid2, bid1)
    ));
}

#[test]
fn record_14_path_uses_other_for_non_sp2_path_without_cis_flags() {
    let mol = Molecule::from_smiles("CCCC").expect("butane");
    let (_mmat, mut accum_data) = run_set13_bounds(&mol);
    let bid1 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
    let bid2 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
    let bid3 = bond_between_idx_simple(&mol, 2, 3).expect("2-3");

    record_14_path(&mol, bid1, bid2, bid3, &mut accum_data).expect("record14Path");

    let path = accum_data.paths14.last().expect("path");
    assert_eq!(path.kind, Path14Kind::Other);
    assert!(!has_path_flag(
        &accum_data.cis_paths,
        path14_id(mol.num_bonds(), bid1, bid2, bid3)
    ));
    assert!(!has_path_flag(
        &accum_data.cis_paths,
        path14_id(mol.num_bonds(), bid3, bid2, bid1)
    ));
}

#[test]
fn set_in_ring_14_bounds_prefers_cis_for_small_sp2_ring_paths() {
    let mol = Molecule::from_smiles("c1ccccc1").expect("benzene");
    let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
    let dmat = flatten_topological_distances(&mol);
    let bid1 = bond_between_idx_simple(&mol, 5, 0).expect("5-0");
    let bid2 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
    let bid3 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
    let rinfo = ring_info_for_distgeom(&mol).expect("ring info");
    set_in_ring_14_bounds(
        &mol,
        bid1,
        bid2,
        bid3,
        &mut accum_data,
        &mut mmat,
        &dmat,
        6,
        &rinfo,
    )
    .expect("in-ring cis bounds");

    let path = accum_data.paths14.last().expect("path");
    assert_eq!(path.kind, Path14Kind::Cis);
    assert!(has_path_flag(
        &accum_data.cis_paths,
        path14_id(mol.num_bonds(), bid1, bid2, bid3)
    ));

    let aid1 = 5usize;
    let aid4 = 2usize;
    let pid = aid1.min(aid4) * mol.num_atoms() + aid1.max(aid4);
    let expected = compute_14_dist_cis(
        accum_data.bond_lengths[bid1],
        accum_data.bond_lengths[bid2],
        accum_data.bond_lengths[bid3],
        accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2),
        accum_data.get_bond_angle(mol.num_bonds(), bid2, bid3),
    );
    assert!(accum_data.visited14_bounds[pid]);
    assert!((mmat.get_lower(aid1, aid4) - (expected - GEN_DIST_TOL)).abs() < 1e-9);
    assert!((mmat.get_upper(aid1, aid4) - (expected + GEN_DIST_TOL)).abs() < 1e-9);
}

#[test]
fn set_two_in_same_ring_14_bounds_uses_trans_for_sp2_external_substituent_path() {
    let mol = Molecule::from_smiles("Cc1ccccc1").expect("toluene");
    let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
    let dmat = flatten_topological_distances(&mol);

    let bid_exo = bond_between_idx_simple(&mol, 0, 1).expect("exo bond");
    let bid_ring_12 = bond_between_idx_simple(&mol, 1, 2).expect("ring bond 1-2");
    let bid_ring_23 = bond_between_idx_simple(&mol, 2, 3).expect("ring bond 2-3");
    set_two_in_same_ring_14_bounds(
        &mol,
        bid_exo,
        bid_ring_12,
        bid_ring_23,
        &mut accum_data,
        &mut mmat,
        &dmat,
    )
    .expect("two-in-same-ring bounds");

    let path = accum_data.paths14.last().expect("path");
    assert_eq!(path.kind, Path14Kind::Trans);
    assert!(has_path_flag(
        &accum_data.trans_paths,
        path14_id(mol.num_bonds(), bid_exo, bid_ring_12, bid_ring_23)
    ));

    let aid1 = 0usize;
    let aid4 = 3usize;
    let expected = compute_14_dist_trans(
        accum_data.bond_lengths[bid_exo],
        accum_data.bond_lengths[bid_ring_12],
        accum_data.bond_lengths[bid_ring_23],
        accum_data.get_bond_angle(mol.num_bonds(), bid_exo, bid_ring_12),
        accum_data.get_bond_angle(mol.num_bonds(), bid_ring_12, bid_ring_23),
    );
    assert!((mmat.get_lower(aid1, aid4) - (expected - GEN_DIST_TOL)).abs() < 1e-9);
    assert!((mmat.get_upper(aid1, aid4) - (expected + GEN_DIST_TOL)).abs() < 1e-9);
}

#[test]
fn diff_ring14_and_share_ring14_delegate_to_in_ring_helper() {
    let mol = Molecule::from_smiles("c1ccccc1").expect("benzene");
    let dmat = flatten_topological_distances(&mol);
    let bid1 = bond_between_idx_simple(&mol, 5, 0).expect("5-0");
    let bid2 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
    let bid3 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
    let rinfo = ring_info_for_distgeom(&mol).expect("ring info");

    let (mut base_mmat, mut base_accum) = run_set13_bounds(&mol);
    set_in_ring_14_bounds(
        &mol,
        bid1,
        bid2,
        bid3,
        &mut base_accum,
        &mut base_mmat,
        &dmat,
        0,
        &rinfo,
    )
    .expect("base in-ring bounds");

    let (mut diff_mmat, mut diff_accum) = run_set13_bounds(&mol);
    set_two_in_diff_ring_14_bounds(
        &mol,
        bid1,
        bid2,
        bid3,
        &mut diff_accum,
        &mut diff_mmat,
        &dmat,
        &rinfo,
    )
    .expect("diff-ring bounds");

    let (mut share_mmat, mut share_accum) = run_set13_bounds(&mol);
    set_share_ring_bond_14_bounds(
        &mol,
        bid1,
        bid2,
        bid3,
        &mut share_accum,
        &mut share_mmat,
        &dmat,
        &rinfo,
    )
    .expect("share-ring-bond bounds");

    assert_eq!(
        diff_accum.paths14.last().map(|p| p.kind),
        base_accum.paths14.last().map(|p| p.kind)
    );
    assert_eq!(
        share_accum.paths14.last().map(|p| p.kind),
        base_accum.paths14.last().map(|p| p.kind)
    );
    assert!((diff_mmat.get_lower(5, 2) - base_mmat.get_lower(5, 2)).abs() < 1e-9);
    assert!((diff_mmat.get_upper(5, 2) - base_mmat.get_upper(5, 2)).abs() < 1e-9);
    assert!((share_mmat.get_lower(5, 2) - base_mmat.get_lower(5, 2)).abs() < 1e-9);
    assert!((share_mmat.get_upper(5, 2) - base_mmat.get_upper(5, 2)).abs() < 1e-9);
}

#[test]
fn compute_15_dist_helpers_match_rdkit_cis_and_trans_formulas() {
    let d1: f64 = 1.41;
    let d2: f64 = 1.52;
    let d3: f64 = 1.38;
    let d4: f64 = 1.47;
    let ang12: f64 = 1.91;
    let ang23: f64 = 2.04;
    let ang34: f64 = 1.88;

    let cis_dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let cis_dy14 = d3 * ang23.sin() - d1 * ang12.sin();
    let cis_d14 = (cis_dx14 * cis_dx14 + cis_dy14 * cis_dy14).sqrt();
    let cis_cval =
        ((d3 - d2 * ang23.cos() + d1 * (ang12 + ang23).cos()) / cis_d14).clamp(-1.0, 1.0);
    let cis_ang143 = cis_cval.acos();
    let expected_cis_cis = compute_13_dist(cis_d14, d4, ang34 - cis_ang143);
    let expected_cis_trans = compute_13_dist(cis_d14, d4, ang34 + cis_ang143);

    let trans_dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let trans_dy14 = d3 * ang23.sin() + d1 * ang12.sin();
    let trans_d14 = (trans_dx14 * trans_dx14 + trans_dy14 * trans_dy14).sqrt();
    let trans_cval =
        ((d3 - d2 * ang23.cos() + d1 * (ang12 - ang23).cos()) / trans_d14).clamp(-1.0, 1.0);
    let trans_ang143 = trans_cval.acos();
    let expected_trans_trans = compute_13_dist(trans_d14, d4, ang34 + trans_ang143);
    let expected_trans_cis = compute_13_dist(trans_d14, d4, ang34 - trans_ang143);

    assert!(
        (compute_15_dist_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34) - expected_cis_cis).abs()
            < 1e-12
    );
    assert!(
        (compute_15_dist_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34) - expected_cis_trans).abs()
            < 1e-12
    );
    assert!(
        (compute_15_dist_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34) - expected_trans_trans)
            .abs()
            < 1e-12
    );
    assert!(
        (compute_15_dist_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34) - expected_trans_cis).abs()
            < 1e-12
    );
}

#[test]
fn compute_15_dist_helpers_support_rdkit_reverse_argument_order() {
    let d1: f64 = 1.41;
    let d2: f64 = 1.52;
    let d3: f64 = 1.38;
    let d4: f64 = 1.47;
    let ang12: f64 = 1.91;
    let ang23: f64 = 2.04;
    let ang34: f64 = 1.88;

    let reversed_cis_cis = compute_15_dist_cis_cis(d4, d3, d2, d1, ang34, ang23, ang12);
    let reversed_cis_trans = compute_15_dist_cis_trans(d4, d3, d2, d1, ang34, ang23, ang12);
    let reversed_trans_cis = compute_15_dist_trans_cis(d4, d3, d2, d1, ang34, ang23, ang12);
    let reversed_trans_trans = compute_15_dist_trans_trans(d4, d3, d2, d1, ang34, ang23, ang12);

    assert_ne!(
        compute_15_dist_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34),
        reversed_cis_cis
    );
    assert_ne!(
        compute_15_dist_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34),
        reversed_cis_trans
    );
    assert_ne!(
        compute_15_dist_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34),
        reversed_trans_cis
    );
    assert_ne!(
        compute_15_dist_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34),
        reversed_trans_trans
    );
}

#[test]
fn set_15_bounds_helper_returns_immediately_for_visited_14_pair() {
    let mol = Molecule::from_smiles("CCCCC").expect("pentane");
    let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
    let dmat = flatten_topological_distances(&mol);
    let nb = mol.num_bonds();
    let na = mol.num_atoms();
    let bid1 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
    let bid2 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
    let bid3 = bond_between_idx_simple(&mol, 2, 3).expect("2-3");
    let pid = 0usize * na + 4usize;
    accum_data.visited14_bounds[pid] = true;
    let before_lower = mmat.get_lower(0, 4);
    let before_upper = mmat.get_upper(0, 4);

    set_15_bounds_helper(
        &mol,
        &mut mmat,
        &mut accum_data,
        &dmat,
        nb,
        na,
        bid1,
        bid2,
        bid3,
        Path14Kind::Other,
    )
    .expect("set15BoundsHelper visited skip");

    assert_eq!(mmat.get_lower(0, 4), before_lower);
    assert_eq!(mmat.get_upper(0, 4), before_upper);
    assert!(!accum_data.set15_atoms[0 * na + 4]);
    assert!(!accum_data.set15_atoms[4 * na + 0]);
}

#[test]
fn set_15_bounds_helper_uses_vdw_fallback_and_marks_set15_atoms_for_unknown_path() {
    let mol = Molecule::from_smiles("CCCCC").expect("pentane");
    let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
    let dmat = flatten_topological_distances(&mol);
    let nb = mol.num_bonds();
    let na = mol.num_atoms();
    let bid1 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
    let bid2 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
    let bid3 = bond_between_idx_simple(&mol, 2, 3).expect("2-3");

    set_15_bounds_helper(
        &mol,
        &mut mmat,
        &mut accum_data,
        &dmat,
        nb,
        na,
        bid1,
        bid2,
        bid3,
        Path14Kind::Other,
    )
    .expect("set15BoundsHelper vdw fallback");

    let expected_lower = VDW_SCALE_15 * (vdw_radius(6) + vdw_radius(6));
    assert!((mmat.get_lower(0, 4) - expected_lower).abs() < 1e-12);
    assert_eq!(mmat.get_upper(0, 4), MAX_UPPER);
    assert!(accum_data.set15_atoms[0 * na + 4]);
    assert!(accum_data.set15_atoms[4 * na + 0]);
}

#[test]
fn set_15_bounds_helper_uses_reversed_other_branch_formula_for_cis_path() {
    let mol = Molecule::from_smiles("CCCCC").expect("pentane");
    let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
    let dmat = flatten_topological_distances(&mol);
    let nb = mol.num_bonds();
    let na = mol.num_atoms();
    let bid1 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
    let bid2 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
    let bid3 = bond_between_idx_simple(&mol, 2, 3).expect("2-3");
    let bid4 = bond_between_idx_simple(&mol, 3, 4).expect("3-4");
    let path_id = bid2 as u64 * nb as u64 * nb as u64 + bid3 as u64 * nb as u64 + bid4 as u64;
    record_path_flag(&mut accum_data.cis_paths, path_id);

    set_15_bounds_helper(
        &mol,
        &mut mmat,
        &mut accum_data,
        &dmat,
        nb,
        na,
        bid1,
        bid2,
        bid3,
        Path14Kind::Other,
    )
    .expect("set15BoundsHelper cis path");

    let d1 = accum_data.bond_lengths[bid1];
    let d2 = accum_data.bond_lengths[bid2];
    let d3 = accum_data.bond_lengths[bid3];
    let d4 = accum_data.bond_lengths[bid4];
    let ang12 = accum_data.get_bond_angle(nb, bid1, bid2);
    let ang23 = accum_data.get_bond_angle(nb, bid2, bid3);
    let ang34 = accum_data.get_bond_angle(nb, bid3, bid4);
    let expected_lower = compute_15_dist_cis_cis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL;
    let expected_upper =
        compute_15_dist_cis_trans(d4, d3, d2, d1, ang34, ang23, ang12) + DIST15_TOL;

    assert!((mmat.get_lower(0, 4) - expected_lower).abs() < 1e-12);
    assert!((mmat.get_upper(0, 4) - expected_upper).abs() < 1e-12);
    assert!(accum_data.set15_atoms[0 * na + 4]);
    assert!(accum_data.set15_atoms[4 * na + 0]);
}

#[test]
fn set_15_bounds_helper_uses_reversed_other_branch_formula_for_trans_path() {
    let mol = Molecule::from_smiles("CCCCC").expect("pentane");
    let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
    let dmat = flatten_topological_distances(&mol);
    let nb = mol.num_bonds();
    let na = mol.num_atoms();
    let bid1 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
    let bid2 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
    let bid3 = bond_between_idx_simple(&mol, 2, 3).expect("2-3");
    let bid4 = bond_between_idx_simple(&mol, 3, 4).expect("3-4");
    let path_id = bid2 as u64 * nb as u64 * nb as u64 + bid3 as u64 * nb as u64 + bid4 as u64;
    record_path_flag(&mut accum_data.trans_paths, path_id);

    set_15_bounds_helper(
        &mol,
        &mut mmat,
        &mut accum_data,
        &dmat,
        nb,
        na,
        bid1,
        bid2,
        bid3,
        Path14Kind::Other,
    )
    .expect("set15BoundsHelper trans path");

    let d1 = accum_data.bond_lengths[bid1];
    let d2 = accum_data.bond_lengths[bid2];
    let d3 = accum_data.bond_lengths[bid3];
    let d4 = accum_data.bond_lengths[bid4];
    let ang12 = accum_data.get_bond_angle(nb, bid1, bid2);
    let ang23 = accum_data.get_bond_angle(nb, bid2, bid3);
    let ang34 = accum_data.get_bond_angle(nb, bid3, bid4);
    let expected_lower =
        compute_15_dist_trans_cis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL;
    let expected_upper =
        compute_15_dist_trans_trans(d4, d3, d2, d1, ang34, ang23, ang12) + DIST15_TOL;

    assert!((mmat.get_lower(0, 4) - expected_lower).abs() < 1e-12);
    assert!((mmat.get_upper(0, 4) - expected_upper).abs() < 1e-12);
    assert!(!has_path_flag(&accum_data.cis_paths, path_id));
    assert!(accum_data.set15_atoms[0 * na + 4]);
    assert!(accum_data.set15_atoms[4 * na + 0]);
}

#[test]
fn set_15_bounds_entrypoint_matches_two_helper_calls_for_single_path() {
    let mol = Molecule::from_smiles("CCCCC").expect("pentane");
    let dmat = flatten_topological_distances(&mol);
    let nb = mol.num_bonds();
    let na = mol.num_atoms();
    let bid1 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
    let bid2 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
    let bid3 = bond_between_idx_simple(&mol, 2, 3).expect("2-3");
    let bid4 = bond_between_idx_simple(&mol, 3, 4).expect("3-4");
    let path_id = bid2 as u64 * nb as u64 * nb as u64 + bid3 as u64 * nb as u64 + bid4 as u64;

    let (mut entry_mmat, mut entry_accum) = run_set13_bounds(&mol);
    entry_accum.paths14.push(Path14Configuration {
        bid1,
        bid2,
        bid3,
        kind: Path14Kind::Other,
    });
    record_path_flag(&mut entry_accum.cis_paths, path_id);
    set_15_bounds(&mol, &mut entry_mmat, &mut entry_accum, &dmat).expect("set15Bounds");

    let (mut helper_mmat, mut helper_accum) = run_set13_bounds(&mol);
    helper_accum.paths14.push(Path14Configuration {
        bid1,
        bid2,
        bid3,
        kind: Path14Kind::Other,
    });
    record_path_flag(&mut helper_accum.cis_paths, path_id);
    set_15_bounds_helper(
        &mol,
        &mut helper_mmat,
        &mut helper_accum,
        &dmat,
        nb,
        na,
        bid1,
        bid2,
        bid3,
        Path14Kind::Other,
    )
    .expect("set15BoundsHelper forward");
    set_15_bounds_helper(
        &mol,
        &mut helper_mmat,
        &mut helper_accum,
        &dmat,
        nb,
        na,
        bid3,
        bid2,
        bid1,
        Path14Kind::Other,
    )
    .expect("set15BoundsHelper reverse");

    assert_eq!(entry_mmat.get_lower(0, 4), helper_mmat.get_lower(0, 4));
    assert_eq!(entry_mmat.get_upper(0, 4), helper_mmat.get_upper(0, 4));
    assert_eq!(entry_accum.set15_atoms, helper_accum.set15_atoms);
}

#[test]
fn set_15_bounds_uses_paths14_produced_by_set_14_bounds() {
    let mol = Molecule::from_smiles("CCCCC").expect("pentane");
    let (mmat, accum_data, _dmat) = run_set15_bounds(&mol, false, false);

    assert!(!accum_data.paths14.is_empty());
    assert!(accum_data.set15_atoms[0 * mol.num_atoms() + 4]);
    assert!(accum_data.set15_atoms[4 * mol.num_atoms() + 0]);
    assert!(mmat.get_lower(0, 4) > 0.0);
    assert!(mmat.get_upper(0, 4) >= mmat.get_lower(0, 4));
}

#[test]
fn set_topol_bounds_can_disable_13_and_14_stages_like_rdkit() {
    let mol = Molecule::from_smiles("CCCC").expect("butane");
    let disabled = run_set_topol_bounds(&mol, false, true, false, false, false, false);
    let enabled = run_set_topol_bounds(&mol, false, true, false, false, true, true);

    assert_eq!(disabled.get_upper(0, 2), MAX_UPPER);
    assert_eq!(disabled.get_upper(0, 3), MAX_UPPER);
    assert!(enabled.get_upper(0, 2) < MAX_UPPER);
    assert!(enabled.get_upper(0, 3) < MAX_UPPER);
}

#[test]
fn set_topol_bounds_can_enable_15_stage_independently_like_rdkit() {
    let mol = Molecule::from_smiles("C/C=C/CC").expect("stereo pentene");
    let without_15 = run_set_topol_bounds(&mol, false, true, false, false, true, true);
    let with_15 = run_set_topol_bounds(&mol, true, true, false, false, true, true);

    assert_eq!(without_15.get_upper(0, 4), MAX_UPPER);
    assert!(with_15.get_upper(0, 4) < MAX_UPPER);
    assert!(with_15.get_lower(0, 4) > without_15.get_lower(0, 4));
}

#[test]
fn set_topol_bounds_ignores_scale_vdw_flag_like_current_rdkit_source() {
    let mol = Molecule::from_smiles("CCCCCC").expect("hexane");
    let scaled = run_set_topol_bounds(&mol, false, true, false, false, false, false);
    let unscaled = run_set_topol_bounds(&mol, false, false, false, false, false, false);

    assert_eq!(scaled.get_lower(0, 5), unscaled.get_lower(0, 5));
    assert_eq!(scaled.get_upper(0, 5), unscaled.get_upper(0, 5));
}

#[test]
fn set_topol_bounds_with_outputs_matches_first_overload_matrix() {
    let mol = Molecule::from_smiles("C/C=C/CC").expect("stereo pentene");
    let plain = run_set_topol_bounds(&mol, true, true, false, false, true, true);
    let (with_outputs, bonds, angles) =
        run_set_topol_bounds_with_outputs(&mol, true, true, false, false, true, true);

    assert_eq!(plain.data, with_outputs.data);
    assert_eq!(bonds.len(), mol.num_bonds());
    assert!(!angles.is_empty());
}

#[test]
fn set_topol_bounds_with_outputs_emits_exact_rdkit_bonds_and_angles_for_triple_bond_case() {
    let mol = Molecule::from_smiles("CC#N").expect("acetonitrile");
    let (_mmat, bonds, angles) =
        run_set_topol_bounds_with_outputs(&mol, false, true, false, false, false, false);

    assert_eq!(bonds, vec![(0, 1), (1, 2)]);
    assert_eq!(angles, vec![vec![0, 1, 2, 1]]);
}

#[test]
fn set_14_bounds_entrypoint_same_ring_matches_direct_helper() {
    let mol = Molecule::from_smiles("c1ccccc1").expect("benzene");
    let (mmat, accum_data, dmat) = run_set14_bounds(&mol, false, false);
    let (bid1, bid2, bid3, ring_size) =
        find_same_ring_dispatch_triple(&mol, false).expect("same-ring triple");
    let rinfo = ring_info_for_distgeom(&mol).expect("ring info");

    let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
    set_in_ring_14_bounds(
        &mol,
        bid1,
        bid2,
        bid3,
        &mut direct_accum,
        &mut direct_mmat,
        &dmat,
        ring_size,
        &rinfo,
    )
    .expect("direct same-ring bounds");

    let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2).expect("shared atom");
    let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3).expect("shared atom");
    let aid1 = if mol.bonds()[bid1].begin().index() == atm2 {
        mol.bonds()[bid1].end().index()
    } else {
        mol.bonds()[bid1].begin().index()
    };
    let aid4 = if mol.bonds()[bid3].begin().index() == atm3 {
        mol.bonds()[bid3].end().index()
    } else {
        mol.bonds()[bid3].begin().index()
    };
    assert!((mmat.get_lower(aid1, aid4) - direct_mmat.get_lower(aid1, aid4)).abs() < 1e-9);
    assert!((mmat.get_upper(aid1, aid4) - direct_mmat.get_upper(aid1, aid4)).abs() < 1e-9);
    assert!(accum_data.visited14_bounds[aid1.min(aid4) * mol.num_atoms() + aid1.max(aid4)]);
}

#[test]
fn set_14_bounds_entrypoint_two_same_ring_matches_direct_helper() {
    let mol = Molecule::from_smiles("Cc1ccccc1").expect("toluene");
    let (mmat, _accum_data, dmat) = run_set14_bounds(&mol, false, false);
    let (bid1, bid2, bid3) = find_dispatch_triple(&mol, Set14DispatchCase::TwoSameRing, false)
        .expect("two-same-ring triple");

    let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
    set_two_in_same_ring_14_bounds(
        &mol,
        bid1,
        bid2,
        bid3,
        &mut direct_accum,
        &mut direct_mmat,
        &dmat,
    )
    .expect("direct two-same-ring bounds");

    let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2).expect("shared atom");
    let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3).expect("shared atom");
    let aid1 = if mol.bonds()[bid1].begin().index() == atm2 {
        mol.bonds()[bid1].end().index()
    } else {
        mol.bonds()[bid1].begin().index()
    };
    let aid4 = if mol.bonds()[bid3].begin().index() == atm3 {
        mol.bonds()[bid3].end().index()
    } else {
        mol.bonds()[bid3].begin().index()
    };
    assert!((mmat.get_lower(aid1, aid4) - direct_mmat.get_lower(aid1, aid4)).abs() < 1e-9);
    assert!((mmat.get_upper(aid1, aid4) - direct_mmat.get_upper(aid1, aid4)).abs() < 1e-9);
}

#[test]
fn set_14_bounds_entrypoint_two_diff_ring_matches_direct_helper() {
    let mol = Molecule::from_smiles("C1CCC2(CC1)CCC3CCCCC23").expect("two-diff-ring polycycle");
    let (mmat, _accum_data, dmat) = run_set14_bounds(&mol, false, false);
    let (bid1, bid2, bid3) = find_dispatch_triple(&mol, Set14DispatchCase::TwoDiffRing, false)
        .expect("two-diff-ring triple");
    let rinfo = ring_info_for_distgeom(&mol).expect("ring info");

    let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
    set_two_in_diff_ring_14_bounds(
        &mol,
        bid1,
        bid2,
        bid3,
        &mut direct_accum,
        &mut direct_mmat,
        &dmat,
        &rinfo,
    )
    .expect("direct two-diff-ring bounds");

    let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2).expect("shared atom");
    let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3).expect("shared atom");
    let aid1 = if mol.bonds()[bid1].begin().index() == atm2 {
        mol.bonds()[bid1].end().index()
    } else {
        mol.bonds()[bid1].begin().index()
    };
    let aid4 = if mol.bonds()[bid3].begin().index() == atm3 {
        mol.bonds()[bid3].end().index()
    } else {
        mol.bonds()[bid3].begin().index()
    };
    assert!((mmat.get_lower(aid1, aid4) - direct_mmat.get_lower(aid1, aid4)).abs() < 1e-9);
    assert!((mmat.get_upper(aid1, aid4) - direct_mmat.get_upper(aid1, aid4)).abs() < 1e-9);
}

#[test]
fn set_14_bounds_entrypoint_share_ring_bond_matches_direct_helper() {
    let mol = Molecule::from_smiles("Cc1ccccc1C").expect("xylene");
    let (mmat, _accum_data, dmat) = run_set14_bounds(&mol, false, false);
    let (bid1, bid2, bid3) = find_dispatch_triple(&mol, Set14DispatchCase::ShareRingBond, false)
        .expect("share-ring-bond triple");
    let rinfo = ring_info_for_distgeom(&mol).expect("ring info");

    let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
    set_share_ring_bond_14_bounds(
        &mol,
        bid1,
        bid2,
        bid3,
        &mut direct_accum,
        &mut direct_mmat,
        &dmat,
        &rinfo,
    )
    .expect("direct share-ring-bond bounds");

    let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2).expect("shared atom");
    let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3).expect("shared atom");
    let aid1 = if mol.bonds()[bid1].begin().index() == atm2 {
        mol.bonds()[bid1].end().index()
    } else {
        mol.bonds()[bid1].begin().index()
    };
    let aid4 = if mol.bonds()[bid3].begin().index() == atm3 {
        mol.bonds()[bid3].end().index()
    } else {
        mol.bonds()[bid3].begin().index()
    };
    assert!((mmat.get_lower(aid1, aid4) - direct_mmat.get_lower(aid1, aid4)).abs() < 1e-9);
    assert!((mmat.get_upper(aid1, aid4) - direct_mmat.get_upper(aid1, aid4)).abs() < 1e-9);
}

#[test]
fn set_14_bounds_entrypoint_chain_matches_direct_helper() {
    let mol = Molecule::from_smiles("CSSC").expect("disulfide");
    let (mmat, _accum_data, _dmat) = run_set14_bounds(&mol, false, false);
    let (bid1, bid2, bid3) =
        find_dispatch_triple(&mol, Set14DispatchCase::Chain, false).expect("chain triple");

    let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
    set_chain_14_bounds(
        &mol,
        bid1,
        bid2,
        bid3,
        &mut direct_accum,
        &mut direct_mmat,
        false,
    )
    .expect("direct chain bounds");

    let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2).expect("shared atom");
    let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3).expect("shared atom");
    let aid1 = if mol.bonds()[bid1].begin().index() == atm2 {
        mol.bonds()[bid1].end().index()
    } else {
        mol.bonds()[bid1].begin().index()
    };
    let aid4 = if mol.bonds()[bid3].begin().index() == atm3 {
        mol.bonds()[bid3].end().index()
    } else {
        mol.bonds()[bid3].begin().index()
    };
    assert!((mmat.get_lower(aid1, aid4) - direct_mmat.get_lower(aid1, aid4)).abs() < 1e-9);
    assert!((mmat.get_upper(aid1, aid4) - direct_mmat.get_upper(aid1, aid4)).abs() < 1e-9);
}

#[test]
fn set_14_bounds_entrypoint_macrocycle_matches_direct_helper() {
    let mol = Molecule::from_smiles("O=C1N(C)CCCCCCCC1").expect("macrocyclic tertiary lactam");
    let (mmat, _accum_data, _dmat) = run_set14_same_ring_pass_only(&mol, true);
    let (bid1, bid2, bid3, _ring_size) =
        find_same_ring_dispatch_triple(&mol, true).expect("macrocycle same-ring triple");

    let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
    set_macrocycle_all_in_same_ring_14_bounds(
        &mol,
        bid1,
        bid2,
        bid3,
        &mut direct_accum,
        &mut direct_mmat,
    )
    .expect("direct macrocycle bounds");
    let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2).expect("shared atom");
    let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3).expect("shared atom");
    let atom1 = if mol.bonds()[bid1].begin().index() == atm2 {
        mol.bonds()[bid1].end().index()
    } else {
        mol.bonds()[bid1].begin().index()
    };
    let atm4 = if mol.bonds()[bid3].begin().index() == atm3 {
        mol.bonds()[bid3].end().index()
    } else {
        mol.bonds()[bid3].begin().index()
    };
    assert!((mmat.get_lower(atom1, atm4) - direct_mmat.get_lower(atom1, atm4)).abs() < 1e-9);
    assert!((mmat.get_upper(atom1, atm4) - direct_mmat.get_upper(atom1, atm4)).abs() < 1e-9);
}

#[test]
fn test_single_atom() {
    let mol = Molecule::from_smiles("C").expect("methane skeleton");
    let result = dg_bounds_matrix(&mol).expect("dg_bounds");
    assert_eq!(result.len(), 1);
    assert_eq!(result[0][0], 0.0);
}

#[test]
fn test_diatomic() {
    let mol = Molecule::from_smiles("CC").expect("ethane skeleton");
    let result = dg_bounds_matrix(&mol).expect("dg_bounds");
    assert_eq!(result.len(), 2);
    assert!(result[0][1] > 0.0);
    assert!(result[0][1] < 5.0);
}

#[test]
fn test_ethane() {
    let mol = Molecule::from_smiles("CC")
        .expect("ethane skeleton")
        .with_hydrogens()
        .expect("explicit hydrogens");
    let result = dg_bounds_matrix(&mol).expect("dg_bounds");
    assert_eq!(result.len(), 8);
    let upper = |i: usize, j: usize| {
        if i < j { result[i][j] } else { result[j][i] }
    };
    let carbon_pair = mol
        .bonds()
        .iter()
        .find_map(|bond| {
            let begin = bond.begin().index();
            let end = bond.end().index();
            (mol.atoms()[begin].atomic_number() == 6 && mol.atoms()[end].atomic_number() == 6)
                .then_some((begin, end))
        })
        .expect("ethane must contain a carbon-carbon bond");
    let carbon_hydrogen_pair = mol
        .bonds()
        .iter()
        .find_map(|bond| {
            let begin = bond.begin().index();
            let end = bond.end().index();
            if mol.atoms()[begin].atomic_number() == 6 && mol.atoms()[end].atomic_number() == 1 {
                Some((begin, end))
            } else if mol.atoms()[begin].atomic_number() == 1
                && mol.atoms()[end].atomic_number() == 6
            {
                Some((end, begin))
            } else {
                None
            }
        })
        .expect("ethane must contain a carbon-hydrogen bond");

    assert!(upper(carbon_pair.0, carbon_pair.1) > 1.0 && upper(carbon_pair.0, carbon_pair.1) < 3.0);
    assert!(
        upper(carbon_hydrogen_pair.0, carbon_hydrogen_pair.1) > 1.0
            && upper(carbon_hydrogen_pair.0, carbon_hydrogen_pair.1) < 3.0
    );
}

#[test]
fn test_empty() {
    let builder = MoleculeBuilder::new();
    let mol = builder.build().expect("build");
    let err = dg_bounds_matrix(&mol).expect_err("empty molecule must now follow setTopolBounds");
    assert!(matches!(
        err,
        DgBoundsError::GenerationFailed(message) if message == "molecule has no atoms"
    ));
}

#[test]
fn dg_bounds_matrix_returns_error_instead_of_panicking_for_3rj7_re_complex() {
    let mol2 = include_str!("../../../../../testdata/mol2/fixtures/3rj7_ligand.mol2");
    let record = read_mol2_from_str(mol2)
        .expect("3rj7 mol2 should parse")
        .expect("3rj7 mol2 should contain a molecule");

    let err = dg_bounds_matrix(&record.molecule).expect_err("invalid bounds must be reported");

    assert!(matches!(
        err,
        DgBoundsError::InvalidBounds(message)
            if message.contains("bad lower bound") && message.contains("atom pair")
    ));
}

#[test]
fn dg_bounds_matrix_matches_source_backed_set_topol_bounds_path() {
    let mol = Molecule::from_smiles("C/C=C/CC").expect("stereo pentene");
    let result = dg_bounds_matrix(&mol).expect("dg_bounds");
    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    set_topol_bounds(&mol, &mut mmat, true, false, false, false, true, true)
        .expect("setTopolBounds");
    assert!(mmat.triangle_smooth(0.0));

    assert_eq!(result, mmat.to_vec_vec());
}

#[test]
fn dg_bounds_matrix_uses_rdkit_wrapper_defaults() {
    let mol = Molecule::from_smiles("CCCCCC").expect("hexane");
    let from_default = dg_bounds_matrix(&mol).expect("default dg_bounds");
    let explicit_default =
        dg_bounds_matrix_with_options(&mol, true, false, true, false).expect("explicit");
    let scaled = dg_bounds_matrix_with_options(&mol, true, true, true, false).expect("scaled vdw");

    assert_eq!(from_default, explicit_default);
    assert_eq!(from_default, scaled);
}

#[test]
fn dg_bounds_matrix_with_options_can_skip_triangle_smoothing() {
    let mol = Molecule::from_smiles("C/C=C/CC").expect("stereo pentene");
    let unsmoothed =
        dg_bounds_matrix_with_options(&mol, true, false, false, false).expect("unsmoothed");
    let smoothed = dg_bounds_matrix_with_options(&mol, true, false, true, false).expect("smoothed");
    let mut manual_unsmoothed = BoundsMatrix::new(mol.num_atoms());
    set_topol_bounds(
        &mol,
        &mut manual_unsmoothed,
        true,
        false,
        false,
        false,
        true,
        true,
    )
    .expect("setTopolBounds");
    let mut manual_smoothed = BoundsMatrix {
        data: manual_unsmoothed.data.clone(),
        n: manual_unsmoothed.n,
    };
    assert!(manual_smoothed.triangle_smooth(0.0));

    assert_eq!(unsmoothed, manual_unsmoothed.data);
    assert_eq!(smoothed, manual_smoothed.data);
}

#[test]
fn dg_bounds_matrix_with_options_forwards_macrocycle14config_like_rdkit_wrapper() {
    let mol = Molecule::from_smiles("C1CCCCCCCCC1").expect("cyclodecane");
    let wrapper_without_macrocycle =
        dg_bounds_matrix_with_options(&mol, true, false, false, false).expect("wrapper plain");
    let wrapper_with_macrocycle =
        dg_bounds_matrix_with_options(&mol, true, false, false, true).expect("wrapper macrocycle");

    let mut manual_without_macrocycle = BoundsMatrix::new(mol.num_atoms());
    set_topol_bounds(
        &mol,
        &mut manual_without_macrocycle,
        true,
        false,
        false,
        false,
        true,
        true,
    )
    .expect("manual plain");

    let mut manual_with_macrocycle = BoundsMatrix::new(mol.num_atoms());
    set_topol_bounds(
        &mol,
        &mut manual_with_macrocycle,
        true,
        false,
        true,
        false,
        true,
        true,
    )
    .expect("manual macrocycle");

    assert_eq!(wrapper_without_macrocycle, manual_without_macrocycle.data);
    assert_eq!(wrapper_with_macrocycle, manual_with_macrocycle.data);
}

#[test]
fn get_atom_stereo_preserves_stereo_when_stereo_atoms_match_query_order() {
    let mut builder = MoleculeBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::C));
    let a3 = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(
            BondSpec::new(a1, a2, BondOrder::Double)
                .with_stereo(BondStereo::Cis)
                .with_stereo_atoms(a0, a3),
        )
        .expect("double bond");
    let mol = builder.build().expect("build");

    assert_eq!(
        get_atom_stereo(&mol.bonds()[0], a0.index(), a3.index()),
        BondStereo::Cis
    );
}

#[test]
fn get_atom_stereo_flips_stereo_when_stereo_atoms_reverse_query_order() {
    let mut builder = MoleculeBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::C));
    let a3 = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(
            BondSpec::new(a1, a2, BondOrder::Double)
                .with_stereo(BondStereo::Cis)
                .with_stereo_atoms(a0, a3),
        )
        .expect("double bond");
    let mol = builder.build().expect("build");

    assert_eq!(
        get_atom_stereo(&mol.bonds()[0], a3.index(), a0.index()),
        BondStereo::Cis
    );
}

#[test]
fn get_atom_stereo_flips_stereo_when_single_end_mismatches() {
    let mut builder = MoleculeBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    let a2 = builder.add_atom(AtomSpec::new(Element::C));
    let a3 = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(
            BondSpec::new(a1, a2, BondOrder::Double)
                .with_stereo(BondStereo::Cis)
                .with_stereo_atoms(a0, a3),
        )
        .expect("double bond");
    let mol = builder.build().expect("build");

    assert_eq!(
        get_atom_stereo(&mol.bonds()[0], a3.index(), a3.index()),
        BondStereo::Trans
    );
}

#[test]
fn triangle_smooth_shared_forwarding_matches_pointer_entrypoint() {
    let mut via_ptr = BoundsMatrix {
        data: vec![vec![0.0; 3]; 3],
        n: 3,
    };
    via_ptr.set_lower(0, 1, 0.5).expect("set lower");
    via_ptr.set_upper(0, 1, 1.0).expect("set upper");
    via_ptr.set_lower(1, 2, 4.0).expect("set lower");
    via_ptr.set_upper(1, 2, 5.0).expect("set upper");
    via_ptr.set_lower(0, 2, 1.0).expect("set lower");
    via_ptr.set_upper(0, 2, 10.0).expect("set upper");

    let mut via_shared = BoundsMatrix {
        data: via_ptr.data.clone(),
        n: via_ptr.n,
    };

    let ptr_result = triangle_smooth_bounds_ptr(&mut via_ptr, 0.0);
    let shared_result = triangle_smooth_bounds_shared(&mut via_shared, 0.0);

    assert_eq!(ptr_result, shared_result);
    assert_eq!(via_ptr.data, via_shared.data);
}

#[test]
fn triangle_smooth_uses_rdkit_difference_formula_to_tighten_lower_bound() {
    let mut mmat = BoundsMatrix {
        data: vec![vec![0.0; 3]; 3],
        n: 3,
    };
    mmat.set_lower(0, 1, 0.5).expect("set lower");
    mmat.set_upper(0, 1, 1.0).expect("set upper");
    mmat.set_lower(1, 2, 4.0).expect("set lower");
    mmat.set_upper(1, 2, 5.0).expect("set upper");
    mmat.set_lower(0, 2, 1.0).expect("set lower");
    mmat.set_upper(0, 2, 10.0).expect("set upper");

    assert!(mmat.triangle_smooth(0.0));

    assert_eq!(mmat.get_upper(0, 2), 6.0);
    assert_eq!(mmat.get_lower(0, 2), 3.0);
    assert_eq!(mmat.data[0][2], 6.0);
    assert_eq!(mmat.data[2][0], 3.0);
}

#[test]
fn triangle_smooth_reconciles_small_lower_upper_inversion_with_tolerance() {
    let mut mmat = BoundsMatrix {
        data: vec![vec![0.0; 3]; 3],
        n: 3,
    };
    mmat.set_lower(0, 1, 2.01).expect("set lower");
    mmat.set_upper(0, 1, 2.0).expect("set upper");
    mmat.set_lower(0, 2, 0.1).expect("set lower");
    mmat.set_upper(0, 2, 100.0).expect("set upper");
    mmat.set_lower(1, 2, 0.1).expect("set lower");
    mmat.set_upper(1, 2, 100.0).expect("set upper");

    assert!(mmat.triangle_smooth(0.01));
    assert_eq!(mmat.get_lower(0, 1), 2.01);
    assert_eq!(mmat.get_upper(0, 1), 2.01);
}

#[test]
fn triangle_smooth_fails_when_lower_exceeds_upper_beyond_tolerance() {
    let mut mmat = BoundsMatrix {
        data: vec![vec![0.0; 3]; 3],
        n: 3,
    };
    mmat.set_lower(0, 1, 2.2).expect("set lower");
    mmat.set_upper(0, 1, 2.0).expect("set upper");
    mmat.set_lower(0, 2, 0.1).expect("set lower");
    mmat.set_upper(0, 2, 100.0).expect("set upper");
    mmat.set_lower(1, 2, 0.1).expect("set lower");
    mmat.set_upper(1, 2, 100.0).expect("set upper");

    assert!(!mmat.triangle_smooth(0.01));
}

#[test]
fn test_ethane_bounds_consistency() {
    // Ethane (C2H6) — verify 1-2, 1-3, 1-4, and VDW bounds are reasonable
    let mol = Molecule::from_smiles("CC")
        .expect("ethane skeleton")
        .with_hydrogens()
        .expect("explicit hydrogens");
    let result = dg_bounds_matrix(&mol).expect("dg_bounds");
    let upper = |i: usize, j: usize| {
        if i < j { result[i][j] } else { result[j][i] }
    };
    assert_eq!(result.len(), 8);
    let carbons: Vec<_> = mol
        .atoms()
        .iter()
        .enumerate()
        .filter_map(|(idx, atom)| (atom.atomic_number() == 6).then_some(idx))
        .collect();
    assert_eq!(carbons.len(), 2);
    let (c1, c2) = mol
        .bonds()
        .iter()
        .find_map(|bond| {
            let begin = bond.begin().index();
            let end = bond.end().index();
            (mol.atoms()[begin].atomic_number() == 6 && mol.atoms()[end].atomic_number() == 6)
                .then_some((begin, end))
        })
        .expect("ethane must contain a carbon-carbon bond");
    let hydrogens_on = |carbon_idx: usize| -> Vec<usize> {
        mol.bonds()
            .iter()
            .filter_map(|bond| {
                let begin = bond.begin().index();
                let end = bond.end().index();
                if begin == carbon_idx && mol.atoms()[end].atomic_number() == 1 {
                    Some(end)
                } else if end == carbon_idx && mol.atoms()[begin].atomic_number() == 1 {
                    Some(begin)
                } else {
                    None
                }
            })
            .collect()
    };
    let h_on_c1 = hydrogens_on(c1);
    let h_on_c2 = hydrogens_on(c2);
    assert!(!h_on_c1.is_empty());
    assert!(!h_on_c2.is_empty());
    let h1 = h_on_c1[0];
    let h4 = h_on_c2[0];
    let lower = |i: usize, j: usize| {
        if i < j { result[j][i] } else { result[i][j] }
    };
    let cc_upper = upper(c1, c2);
    let cc_lower = lower(c1, c2);
    let ch_upper = upper(c1, h1);
    let ch_lower = lower(c1, h1);
    let cch_upper = upper(c1, h4);
    let hcch_upper = upper(h1, h4);

    // Directly bonded pairs should keep valid finite bounds.
    assert!(
        cc_lower > 0.0 && cc_upper > cc_lower,
        "C-C bounds should be finite and ordered, got [{cc_lower}, {cc_upper}]"
    );
    assert!(
        ch_lower > 0.0 && ch_upper > ch_lower,
        "C-H bounds should be finite and ordered, got [{ch_lower}, {ch_upper}]"
    );

    // Longer topological separations should not be tighter than direct bonds.
    assert!(
        cch_upper > ch_upper,
        "C-C-H 1-3 upper bound {cch_upper} should exceed direct C-H upper bound {ch_upper}"
    );

    assert!(
        hcch_upper > cch_upper,
        "H-C-C-H 1-4 upper bound {hcch_upper} should exceed C-C-H 1-3 upper bound {cch_upper}"
    );

    // Matrix symmetry check
    for i in 0..8 {
        for j in 0..8 {
            assert!(!result[i][j].is_nan(), "bounds[{i}][{j}] is NaN");
            assert!(
                result[i][j] >= 0.0 || i == j,
                "bounds[{i}][{j}] = {} should be >= 0",
                result[i][j]
            );
        }
    }
}

#[test]
fn test_ring_bounds_consistency() {
    // Cyclohexane ring skeleton — this exercises ring-aware angle handling
    // This tests ring-aware angle computation and triangle smoothing
    let mol = Molecule::from_smiles("C1CCCCC1").expect("cyclohexane");
    let result = dg_bounds_matrix(&mol).expect("dg_bounds");
    let upper = |i: usize, j: usize| {
        if i < j { result[i][j] } else { result[j][i] }
    };
    assert_eq!(result.len(), 6);

    // 1-2 (C-C bonds): upper bound ~1.51
    for i in 0..6 {
        let j = (i + 1) % 6;
        assert!(
            upper(i, j) > 1.40 && upper(i, j) < 1.60,
            "C-C 1-2 upper bound [{i}][{j}] = {} should be ~1.50",
            upper(i, j)
        );
    }

    // 1-3 (C-C-C in ring): should be larger than 1-2
    for i in 0..6 {
        let j = (i + 2) % 6;
        assert!(
            upper(i, j) > upper(i, (i + 1) % 6),
            "1-3 upper bound [{i}][{j}] = {} should exceed 1-2 bound {}",
            upper(i, j),
            upper(i, (i + 1) % 6)
        );
    }

    // 1-4 (C-C-C-C across hexagon): should be the largest intra-ring distance
    for i in 0..6 {
        let j = (i + 3) % 6;
        assert!(
            upper(i, j) > upper(i, (i + 1) % 6) && upper(i, j) < 5.0,
            "1-4 upper bound [{i}][{j}] = {} should be > 1-2 and < 5.0",
            upper(i, j)
        );
    }

    assert!(!upper(0, 1).is_nan());
    assert!(!upper(0, 1).is_infinite());
}
