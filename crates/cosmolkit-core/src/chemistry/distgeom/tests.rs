use super::*;
use crate::builder::MoleculeBuilder;
use crate::{AtomSpec, BondSpec, Element, Molecule, ValenceModel, assign_valence};

fn run_set12_bounds(mol: &Molecule) -> (BoundsMatrix, ComputedData) {
    let mut mmat = BoundsMatrix::new(mol.num_atoms());
    let mut accum_data = ComputedData::new(mol.num_atoms(), mol.num_bonds());
    set_12_bounds(mol, &mut mmat, &mut accum_data).expect("set12Bounds");
    (mmat, accum_data)
}

fn run_set13_bounds(mol: &Molecule) -> (BoundsMatrix, ComputedData) {
    let (mut mmat, mut accum_data) = run_set12_bounds(mol);
    set_13_bounds(mol, &mut mmat, &mut accum_data);
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
    );
    (mmat, accum_data, dmat)
}

fn run_set15_bounds(
    mol: &Molecule,
    use_macrocycle_14config: bool,
    force_trans_amides: bool,
) -> (BoundsMatrix, ComputedData, Vec<f64>) {
    let (mut mmat, mut accum_data, dmat) =
        run_set14_bounds(mol, use_macrocycle_14config, force_trans_amides);
    set_15_bounds(mol, &mut mmat, &mut accum_data, &dmat);
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
    let rinfo = ring_info_for_distgeom(mol);
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
                    );
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
                    );
                }
            } else {
                record_14_path(mol, bid1, bid2, bid3, &mut accum_data);
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
    let rinfo = ring_info_for_distgeom(mol);
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
    let rinfo = ring_info_for_distgeom(mol);
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
    let topo = compute_topological_distances(mol);
    let n = mol.num_atoms();
    let mut flat = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            flat[i * n + j] = topo[i][j] as f64;
        }
    }
    flat
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

    mmat.set_upper(2, 0, 4.2);
    mmat.set_lower(2, 0, 1.8);

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
    mmat.set_upper(0, 1, 2.5);
    mmat.set_lower(0, 1, 1.25);

    let raw = mmat.to_vec_vec();
    assert_eq!(raw, vec![vec![0.0, 2.5], vec![1.25, 0.0]]);
}

#[test]
fn check_and_set_bounds_sets_uninitialized_pair_conservatively() {
    let mut mmat = BoundsMatrix::new(3);

    mmat.check_and_set_bounds(0, 2, 1.2, 2.8);

    assert_eq!(mmat.get_lower(0, 2), 1.2);
    assert_eq!(mmat.get_upper(0, 2), 2.8);
    assert_eq!(mmat.data[2][0], 1.2);
    assert_eq!(mmat.data[0][2], 2.8);
}

#[test]
fn check_and_set_bounds_only_tightens_lower_and_relaxes_upper_when_allowed() {
    let mut mmat = BoundsMatrix::new(3);
    mmat.set_lower(0, 2, 1.8);
    mmat.set_upper(0, 2, 2.2);

    mmat.check_and_set_bounds(0, 2, 1.4, 2.6);
    assert_eq!(mmat.get_lower(0, 2), 1.4);
    assert_eq!(mmat.get_upper(0, 2), 2.6);

    mmat.check_and_set_bounds(0, 2, 1.6, 2.0);
    assert_eq!(mmat.get_lower(0, 2), 1.4);
    assert_eq!(mmat.get_upper(0, 2), 2.6);
}

#[test]
fn check_and_set_bounds_with_mode_uses_overlap_when_ranges_are_consistent() {
    let mut mmat = BoundsMatrix::new(2);
    mmat.set_lower(0, 1, 1.3);
    mmat.set_upper(0, 1, 3.1);

    mmat.check_and_set_bounds_with_mode(0, 1, 1.8, 2.4, true);

    assert_eq!(mmat.get_lower(0, 1), 1.8);
    assert_eq!(mmat.get_upper(0, 1), 2.4);
    assert_eq!(mmat.data[1][0], 1.8);
    assert_eq!(mmat.data[0][1], 2.4);
}

#[test]
fn check_and_set_bounds_with_mode_falls_back_to_conservative_union_when_disjoint() {
    let mut mmat = BoundsMatrix::new(2);
    mmat.set_lower(0, 1, 1.3);
    mmat.set_upper(0, 1, 2.0);

    mmat.check_and_set_bounds_with_mode(0, 1, 2.4, 3.0, true);

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

    set_lower_bound_vdw(&mol, &mut mmat, true, &dmat);

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

    set_lower_bound_vdw(&mol, &mut mmat, true, &dmat);

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
fn set_13_bounds_helper_doubles_tolerance_for_larger_sp2_ring_atoms() {
    let thiophene = Molecule::from_smiles("s1cccc1").expect("thiophene");
    let benzene = Molecule::from_smiles("c1ccccc1").expect("benzene");

    let (_thiophene_mmat_12, thiophene_accum) = run_set12_bounds(&thiophene);
    let (_benzene_mmat_12, benzene_accum) = run_set12_bounds(&benzene);

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
    );

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
    );

    let thiophene_width = thiophene_mmat.get_upper(sulfur_neighbors[0], sulfur_neighbors[1])
        - thiophene_mmat.get_lower(sulfur_neighbors[0], sulfur_neighbors[1]);
    let benzene_width = benzene_mmat.get_upper(0, 2) - benzene_mmat.get_lower(0, 2);

    assert!((thiophene_width - (4.0 * DIST13_TOL)).abs() < 1e-9);
    assert!((benzene_width - (2.0 * DIST13_TOL)).abs() < 1e-9);
    assert!(is_larger_sp2_atom_idx(&thiophene, sulfur));
    assert!(!is_larger_sp2_atom_idx(&benzene, 1));
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
    );

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
    set_macrocycle_all_in_same_ring_14_bounds(&mol, bid1, bid2, bid3, &mut accum_data, &mut mmat);

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
    );
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
    );
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
    set_chain_14_bounds(&mol, bid1, bid2, bid3, &mut accum_data, &mut mmat, false);

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
    );
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
    );
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
    record_14_path(&mol, bid1, bid2, bid3, &mut accum_data);

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
fn set_in_ring_14_bounds_prefers_cis_for_small_sp2_ring_paths() {
    let mol = Molecule::from_smiles("c1ccccc1").expect("benzene");
    let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
    let dmat = flatten_topological_distances(&mol);
    let bid1 = bond_between_idx_simple(&mol, 5, 0).expect("5-0");
    let bid2 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
    let bid3 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
    set_in_ring_14_bounds(&mol, bid1, bid2, bid3, &mut accum_data, &mut mmat, &dmat, 6);

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
    );

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
    );

    let (mut diff_mmat, mut diff_accum) = run_set13_bounds(&mol);
    set_two_in_diff_ring_14_bounds(
        &mol,
        bid1,
        bid2,
        bid3,
        &mut diff_accum,
        &mut diff_mmat,
        &dmat,
    );

    let (mut share_mmat, mut share_accum) = run_set13_bounds(&mol);
    set_share_ring_bond_14_bounds(
        &mol,
        bid1,
        bid2,
        bid3,
        &mut share_accum,
        &mut share_mmat,
        &dmat,
    );

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
    );

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
    );

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
    );

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
    set_15_bounds(&mol, &mut entry_mmat, &mut entry_accum, &dmat);

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
    );
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
    );

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
    );

    let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2);
    let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3);
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
    );

    let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2);
    let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3);
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

    let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
    set_two_in_diff_ring_14_bounds(
        &mol,
        bid1,
        bid2,
        bid3,
        &mut direct_accum,
        &mut direct_mmat,
        &dmat,
    );

    let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2);
    let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3);
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

    let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
    set_share_ring_bond_14_bounds(
        &mol,
        bid1,
        bid2,
        bid3,
        &mut direct_accum,
        &mut direct_mmat,
        &dmat,
    );

    let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2);
    let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3);
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
    );

    let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2);
    let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3);
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
    );
    let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2);
    let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3);
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
fn triangle_smooth_uses_rdkit_difference_formula_to_tighten_lower_bound() {
    let mut mmat = BoundsMatrix {
        data: vec![vec![0.0; 3]; 3],
        n: 3,
    };
    mmat.set_lower(0, 1, 0.5);
    mmat.set_upper(0, 1, 1.0);
    mmat.set_lower(1, 2, 4.0);
    mmat.set_upper(1, 2, 5.0);
    mmat.set_lower(0, 2, 1.0);
    mmat.set_upper(0, 2, 10.0);

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
    mmat.set_lower(0, 1, 2.01);
    mmat.set_upper(0, 1, 2.0);
    mmat.set_lower(0, 2, 0.1);
    mmat.set_upper(0, 2, 100.0);
    mmat.set_lower(1, 2, 0.1);
    mmat.set_upper(1, 2, 100.0);

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
    mmat.set_lower(0, 1, 2.2);
    mmat.set_upper(0, 1, 2.0);
    mmat.set_lower(0, 2, 0.1);
    mmat.set_upper(0, 2, 100.0);
    mmat.set_lower(1, 2, 0.1);
    mmat.set_upper(1, 2, 100.0);

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
