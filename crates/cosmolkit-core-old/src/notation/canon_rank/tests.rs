use super::*;
use crate::{AtomSpec, BondSpec, Element, MoleculeBuilder, StereoGroup};

#[test]
#[should_panic(expected = "size mismatch")]
fn shared_count_swaps_canon_preserves_size_precondition_panic() {
    let _ = count_swaps_to_interconvert(&[1_u8, 2], vec![1]);
}

#[test]
#[should_panic(expected = "could not find probe element")]
fn shared_count_swaps_canon_preserves_missing_element_invariant_panic() {
    let _ = count_swaps_to_interconvert(&[1_u8, 2], vec![1, 3]);
}

fn init_fragment_canon_atoms_for_kekulize<'a>(
    molecule: &Molecule,
    atoms_to_use: &[bool],
    bonds_to_use: &[bool],
    include_chirality: bool,
    include_stereo_groups: bool,
    atom_symbols: Option<&'a [String]>,
    bond_symbols: Option<&'a [String]>,
) -> Result<Vec<CanonAtom<'a>>, KekulizeError> {
    let view = CanonRankReadView::from_molecule(molecule)?;
    super::init_fragment_canon_atoms_for_kekulize(
        &view,
        atoms_to_use,
        bonds_to_use,
        include_chirality,
        include_stereo_groups,
        atom_symbols,
        bond_symbols,
    )
}

fn compare_ring_atoms_concerning_num_neighbors_for_kekulize(
    molecule: &Molecule,
    atoms: &mut [CanonAtom<'_>],
) -> Result<(), KekulizeError> {
    let view = CanonRankReadView::from_molecule(molecule)?;
    super::compare_ring_atoms_concerning_num_neighbors_for_kekulize(&view, atoms);
    Ok(())
}

fn make_canon_bond_holder<'a>(
    molecule: &Molecule,
    bond_idx: usize,
    other_idx: usize,
    include_chirality: bool,
) -> Result<CanonBondHolder<'a>, KekulizeError> {
    let view = CanonRankReadView::from_molecule(molecule)?;
    super::make_canon_bond_holder(&view, bond_idx, other_idx, include_chirality)
}

fn atropisomer_atoms_and_bonds_for_canonical_rank(
    molecule: &Molecule,
    bond_id: BondId,
) -> Option<[(AtomId, Vec<BondId>); 2]> {
    let view = CanonRankReadView::from_molecule(molecule).ok()?;
    super::atropisomer_atoms_and_bonds_for_canonical_rank(&view, bond_id)
}

#[allow(clippy::too_many_arguments)]
fn break_ties_for_kekulize(
    molecule: &Molecule,
    atoms: &mut [CanonAtom<'_>],
    mode: bool,
    compare_mode: CanonCompareMode,
    flags: CanonRankFlags,
    order: &mut [usize],
    count: &mut [usize],
    active_set: &mut isize,
    next: &mut [isize],
    changed: &mut [bool],
    touched_partitions: &mut [bool],
) {
    let view = CanonRankReadView::from_molecule(molecule).unwrap();
    super::break_ties_for_kekulize(
        &view,
        atoms,
        mode,
        compare_mode,
        flags,
        order,
        count,
        active_set,
        next,
        changed,
        touched_partitions,
    );
}

#[test]
fn rank_mol_atoms_returns_empty_for_empty_molecule() {
    let molecule = Molecule::default();

    let ranks = rank_mol_atoms(&molecule).unwrap();

    assert!(ranks.is_empty());
}

#[test]
fn rank_mol_atoms_break_ties_produces_unique_ranks() {
    let molecule = Molecule::from_smiles_with_sanitize("CC", false).unwrap();

    let ranks = rank_mol_atoms(&molecule).unwrap();
    let mut sorted = ranks.clone();
    sorted.sort_unstable();

    assert_eq!(sorted, vec![0, 1]);
}

#[test]
fn rank_mol_atoms_without_break_ties_preserves_symmetry_classes() {
    let molecule = Molecule::from_smiles_with_sanitize("CC", false).unwrap();
    let options = CanonicalRankOptions {
        break_ties: false,
        ..CanonicalRankOptions::rank_mol_default()
    };

    let ranks = rank_mol_atoms_with_options(&molecule, options).unwrap();

    assert_eq!(ranks, vec![0, 0]);
}

#[test]
fn rank_mol_options_gate_chirality_ring_stereo_like_rdkit_rank_mol_atoms() {
    let mut rank_mol_options = CanonicalRankOptions::rank_mol_default();
    rank_mol_options.include_ring_stereo = false;
    let rank_mol_flags = CanonRankFlags::from_fragment_options(rank_mol_options);

    let mut fragment_options = CanonicalRankOptions::kekulize_default();
    fragment_options.include_ring_stereo = false;
    let fragment_flags = CanonRankFlags::from_fragment_options(fragment_options);

    assert!(!rank_mol_flags.use_chirality_rings);
    assert!(fragment_flags.use_chirality_rings);
}

#[test]
fn rank_fragment_atoms_accepts_symbol_overrides() {
    let molecule = Molecule::from_smiles_with_sanitize("CC", false).unwrap();
    let atoms_to_use = vec![true; molecule.num_atoms()];
    let bonds_to_use = vec![true; molecule.num_bonds()];
    let atom_symbols = vec!["A".to_string(), "B".to_string()];
    let scope =
        FragmentRankScope::new(&atoms_to_use, &bonds_to_use).with_atom_symbols(&atom_symbols);

    let ranks =
        rank_fragment_atoms(&molecule, scope, CanonicalRankOptions::kekulize_default()).unwrap();

    assert_eq!(ranks.len(), molecule.num_atoms());
}

#[test]
fn rank_fragment_atoms_accepts_non_default_option_variants() {
    let molecule = Molecule::from_smiles_with_sanitize("CC", false).unwrap();
    let atoms_to_use = vec![true; molecule.num_atoms()];
    let bonds_to_use = vec![true; molecule.num_bonds()];
    let scope = FragmentRankScope::new(&atoms_to_use, &bonds_to_use);
    let mut options = CanonicalRankOptions::kekulize_default();
    options.include_atom_maps = false;
    options.include_isotopes = false;
    options.include_chiral_presence = true;
    options.include_chirality = false;
    options.include_ring_stereo = false;

    let ranks = rank_fragment_atoms(&molecule, scope, options).unwrap();

    assert_eq!(ranks.len(), molecule.num_atoms());
}

#[test]
fn atom_symbol_override_short_circuits_base_atom_comparison() {
    let mut builder = MoleculeBuilder::new();
    builder.add_atom(AtomSpec::new(Element::C));
    builder.add_atom(AtomSpec::new(Element::O));
    let molecule = builder.build().unwrap();
    let atoms_to_use = vec![true; molecule.num_atoms()];
    let bonds_to_use = vec![true; molecule.num_bonds()];
    let atom_symbols = vec!["X".to_string(), "X".to_string()];
    let mut atoms = init_fragment_canon_atoms_for_kekulize(
        &molecule,
        &atoms_to_use,
        &bonds_to_use,
        true,
        true,
        Some(&atom_symbols),
        None,
    )
    .unwrap();
    atoms[0].index = 0;
    atoms[1].index = 0;

    let ordering = compare_canon_atom_base_for_kekulize(
        &atoms,
        0,
        1,
        CanonRankFlags::from_fragment_options(CanonicalRankOptions::kekulize_default()),
    );

    assert_eq!(ordering, Ordering::Equal);
}

#[test]
fn compare_ring_atoms_concerning_num_neighbors_out_of_play_atoms_are_equal_in_atom_compare() {
    let mut builder = MoleculeBuilder::new();
    builder.add_atom(AtomSpec::new(Element::O));
    builder.add_atom(AtomSpec::new(Element::C));
    let molecule = builder.build().unwrap();
    let atoms_to_use = vec![false; molecule.num_atoms()];
    let bonds_to_use = vec![true; molecule.num_bonds()];
    let mut atoms = init_fragment_canon_atoms_for_kekulize(
        &molecule,
        &atoms_to_use,
        &bonds_to_use,
        true,
        true,
        None,
        None,
    )
    .unwrap();

    let ordering = compare_canon_atoms_for_kekulize(
        &mut atoms,
        0,
        1,
        CanonCompareMode::Atom,
        CanonRankFlags::from_fragment_options(CanonicalRankOptions::kekulize_default()),
    );

    assert_eq!(ordering, Ordering::Equal);
}

#[test]
fn compare_ring_atoms_concerning_num_neighbors_out_of_play_atoms_are_equal_in_special_chirality_compare()
 {
    let mut builder = MoleculeBuilder::new();
    builder.add_atom(AtomSpec::new(Element::C));
    builder.add_atom(AtomSpec::new(Element::C));
    let molecule = builder.build().unwrap();
    let atoms_to_use = vec![false; molecule.num_atoms()];
    let bonds_to_use = vec![true; molecule.num_bonds()];
    let mut atoms = init_fragment_canon_atoms_for_kekulize(
        &molecule,
        &atoms_to_use,
        &bonds_to_use,
        true,
        true,
        None,
        None,
    )
    .unwrap();
    atoms[0].bonds.push(CanonBondHolder {
        bond_type: BondOrder::Single,
        bond_stereo: 0,
        stype: BondStereo::None,
        controlling_atoms: [None; 4],
        nbr_sym_class: 1,
        nbr_idx: 1,
        p_symbol: None,
        bond_idx: 0,
    });
    atoms[1].bonds.push(CanonBondHolder {
        bond_type: BondOrder::Single,
        bond_stereo: 0,
        stype: BondStereo::None,
        controlling_atoms: [None; 4],
        nbr_sym_class: 2,
        nbr_idx: 0,
        p_symbol: None,
        bond_idx: 0,
    });

    let ordering = compare_special_chirality_atoms_for_kekulize(&mut atoms, 0, 1);

    assert_eq!(ordering, Ordering::Equal);
}

#[test]
fn compare_ring_atoms_concerning_num_neighbors_out_of_play_atoms_are_equal_in_special_symmetry_compare()
 {
    let mut builder = MoleculeBuilder::new();
    builder.add_atom(AtomSpec::new(Element::C));
    builder.add_atom(AtomSpec::new(Element::C));
    let molecule = builder.build().unwrap();
    let atoms_to_use = vec![false; molecule.num_atoms()];
    let bonds_to_use = vec![true; molecule.num_bonds()];
    let mut atoms = init_fragment_canon_atoms_for_kekulize(
        &molecule,
        &atoms_to_use,
        &bonds_to_use,
        true,
        true,
        None,
        None,
    )
    .unwrap();
    atoms[0].neighbor_num = vec![1, -1];
    atoms[1].neighbor_num = vec![2, -1];
    atoms[0].revisted_neighbors = vec![3];
    atoms[1].revisted_neighbors = vec![4];

    let ordering = compare_special_symmetry_atoms_for_kekulize(&mut atoms, 0, 1);

    assert_eq!(ordering, Ordering::Equal);
}

#[test]
fn rank_options_control_isotope_atom_map_and_chiral_presence_comparison() {
    let mut builder = MoleculeBuilder::new();
    builder.add_atom(AtomSpec::new(Element::C).with_isotope(13).with_atom_map(7));
    builder.add_atom(
        AtomSpec::new(Element::C)
            .with_isotope(12)
            .with_atom_map(1)
            .with_chiral_tag(ChiralTag::TetrahedralCw),
    );
    let molecule = builder.build().unwrap();
    let atoms_to_use = vec![true; molecule.num_atoms()];
    let bonds_to_use = vec![true; molecule.num_bonds()];
    let mut atoms = init_fragment_canon_atoms_for_kekulize(
        &molecule,
        &atoms_to_use,
        &bonds_to_use,
        true,
        true,
        None,
        None,
    )
    .unwrap();
    atoms[0].index = 0;
    atoms[1].index = 0;

    let mut options = CanonicalRankOptions::kekulize_default();
    assert_eq!(
        compare_canon_atom_base_for_kekulize(
            &atoms,
            0,
            1,
            CanonRankFlags::from_fragment_options(options)
        ),
        Ordering::Greater
    );

    options.include_atom_maps = false;
    options.include_isotopes = false;
    options.include_chirality = false;
    assert_eq!(
        compare_canon_atom_base_for_kekulize(
            &atoms,
            0,
            1,
            CanonRankFlags::from_fragment_options(options)
        ),
        Ordering::Equal
    );

    options.include_chiral_presence = true;
    assert_eq!(
        compare_canon_atom_base_for_kekulize(
            &atoms,
            0,
            1,
            CanonRankFlags::from_fragment_options(options)
        ),
        Ordering::Less
    );
}

#[test]
fn special_symmetry_ring_neighbor_state_is_populated_for_ring_atoms() {
    let molecule = Molecule::from_smiles_with_sanitize("c1ccc2ccccc2c1", false).unwrap();
    let atoms_to_use = vec![true; molecule.num_atoms()];
    let bonds_to_use = vec![true; molecule.num_bonds()];
    let mut atoms = init_fragment_canon_atoms_for_kekulize(
        &molecule,
        &atoms_to_use,
        &bonds_to_use,
        true,
        true,
        None,
        None,
    )
    .unwrap();

    compare_ring_atoms_concerning_num_neighbors_for_kekulize(&molecule, &mut atoms).unwrap();

    assert!(atoms.iter().any(|atom| !atom.neighbor_num.is_empty()));
    assert!(atoms.iter().any(|atom| !atom.revisted_neighbors.is_empty()));
}

#[test]
fn special_symmetry_ring_neighbor_state_skips_acyclic_atoms_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("Cc1ccccc1", false).unwrap();
    let atoms_to_use = vec![true; molecule.num_atoms()];
    let bonds_to_use = vec![true; molecule.num_bonds()];
    let mut atoms = init_fragment_canon_atoms_for_kekulize(
        &molecule,
        &atoms_to_use,
        &bonds_to_use,
        true,
        true,
        None,
        None,
    )
    .unwrap();

    compare_ring_atoms_concerning_num_neighbors_for_kekulize(&molecule, &mut atoms).unwrap();

    assert!(atoms[0].neighbor_num.is_empty());
    assert!(atoms[0].revisted_neighbors.is_empty());
    assert!(atoms[1..].iter().all(|atom| !atom.neighbor_num.is_empty()));
}

#[test]
fn rank_fragment_atoms_handles_fused_ring_special_symmetry_path() {
    let molecule = Molecule::from_smiles_with_sanitize("c1ccc2ccccc2c1", false).unwrap();
    let atoms_to_use = vec![true; molecule.num_atoms()];
    let bonds_to_use = vec![true; molecule.num_bonds()];
    let scope = FragmentRankScope::new(&atoms_to_use, &bonds_to_use);

    let ranks =
        rank_fragment_atoms(&molecule, scope, CanonicalRankOptions::kekulize_default()).unwrap();

    assert_eq!(ranks.len(), molecule.num_atoms());
}

#[test]
fn hanoi_merges_equal_pair_into_single_partition() {
    let molecule = Molecule::from_smiles_with_sanitize("CC", false).unwrap();
    let atoms_to_use = vec![true; molecule.num_atoms()];
    let bonds_to_use = vec![true; molecule.num_bonds()];
    let mut atoms = init_fragment_canon_atoms_for_kekulize(
        &molecule,
        &atoms_to_use,
        &bonds_to_use,
        true,
        true,
        None,
        None,
    )
    .unwrap();
    atoms[0].index = 0;
    atoms[1].index = 0;

    let mut order = vec![0usize, 1usize];
    let mut temp = vec![0usize; 2];
    let mut count = vec![0usize; 2];
    let changed = vec![true; 2];

    let in_temp = hanoi_order_for_kekulize(
        &mut order,
        &mut temp,
        &mut count,
        &changed,
        &mut atoms,
        CanonCompareMode::Atom,
        CanonRankFlags::from_fragment_options(CanonicalRankOptions::kekulize_default()),
    );

    assert!(!in_temp);
    assert_eq!(order, vec![0, 1]);
    assert_eq!(count, vec![2, 0]);
}

#[test]
fn hanoi_merges_odd_length_partition_order_like_rdkit() {
    let mut builder = MoleculeBuilder::new();
    for _ in 0..3 {
        builder.add_atom(AtomSpec::new(Element::C));
    }
    let molecule = builder.build().unwrap();
    let atoms_to_use = vec![true; molecule.num_atoms()];
    let bonds_to_use = vec![true; molecule.num_bonds()];
    let mut atoms = init_fragment_canon_atoms_for_kekulize(
        &molecule,
        &atoms_to_use,
        &bonds_to_use,
        true,
        true,
        None,
        None,
    )
    .unwrap();
    atoms[0].index = 2;
    atoms[1].index = 0;
    atoms[2].index = 1;

    let mut order = vec![0usize, 1usize, 2usize];
    let mut temp = vec![0usize; 3];
    let mut count = vec![0usize; 3];
    let changed = vec![true; 3];

    let in_temp = hanoi_order_for_kekulize(
        &mut order,
        &mut temp,
        &mut count,
        &changed,
        &mut atoms,
        CanonCompareMode::Atom,
        CanonRankFlags::from_fragment_options(CanonicalRankOptions::kekulize_default()),
    );

    assert!(in_temp);
    assert_eq!(temp, vec![1, 2, 0]);
    assert_eq!(count, vec![1, 1, 1]);
}

#[test]
fn hanoi_break_ties_revisits_first_partition_after_refine_reorders_head() {
    let mut builder = MoleculeBuilder::new();
    for _ in 0..4 {
        builder.add_atom(AtomSpec::new(Element::C));
    }
    let molecule = builder.build().unwrap();
    let atoms_to_use = vec![true; molecule.num_atoms()];
    let bonds_to_use = vec![true; molecule.num_bonds()];
    let mut atoms = init_fragment_canon_atoms_for_kekulize(
        &molecule,
        &atoms_to_use,
        &bonds_to_use,
        true,
        true,
        None,
        None,
    )
    .unwrap();

    atoms[0].index = 0;
    atoms[1].index = 0;
    atoms[2].index = 0;
    atoms[3].index = 0;
    atoms[0].atomic_number = 8;
    atoms[1].atomic_number = 6;
    atoms[2].atomic_number = 6;
    atoms[3].atomic_number = 6;
    atoms[0].degree = 0;
    atoms[1].degree = 0;
    atoms[2].degree = 0;
    atoms[3].degree = 1;
    atoms[0].nbr_ids.clear();
    atoms[1].nbr_ids.clear();
    atoms[2].nbr_ids.clear();
    atoms[3].nbr_ids = vec![0];
    atoms[0].bonds.clear();
    atoms[1].bonds.clear();
    atoms[2].bonds.clear();
    atoms[3].bonds.clear();

    let mut order = vec![0usize, 1usize, 2usize, 3usize];
    let mut count = vec![4usize, 0usize, 0usize, 0usize];
    let mut active_set = -1isize;
    let mut next = vec![-2isize; 4];
    let mut changed = vec![false; 4];
    let mut touched = vec![false; 4];

    break_ties_for_kekulize(
        &molecule,
        &mut atoms,
        true,
        CanonCompareMode::Atom,
        CanonRankFlags::from_fragment_options(CanonicalRankOptions::kekulize_default()),
        &mut order,
        &mut count,
        &mut active_set,
        &mut next,
        &mut changed,
        &mut touched,
    );

    assert_eq!(order, vec![1, 2, 0, 3]);
    assert_eq!(count, vec![1, 1, 1, 1]);
    assert_eq!(atoms[0].index, 2);
    assert_eq!(atoms[1].index, 0);
    assert_eq!(atoms[2].index, 1);
    assert_eq!(atoms[3].index, 3);
}

#[test]
fn canonical_rank_state_records_ring_stereo_neighbor_and_stereo_group_order() {
    let mut builder = MoleculeBuilder::new();
    let left = builder.add_atom(AtomSpec::new(Element::C));
    let center = builder.add_atom(
        AtomSpec::new(Element::C)
            .with_chiral_tag(ChiralTag::TetrahedralCw)
            .with_prop("_ringStereoAtoms", "1"),
    );
    let right = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(BondSpec::new(left, center, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(center, right, BondOrder::Single))
        .unwrap();
    builder
        .add_stereo_group(StereoGroup::new(
            StereoGroupKind::Or,
            vec![center],
            Vec::new(),
        ))
        .unwrap();
    let molecule = builder.build().unwrap();
    let atoms_to_use = vec![true; molecule.num_atoms()];
    let bonds_to_use = vec![true; molecule.num_bonds()];

    let atoms = init_fragment_canon_atoms_for_kekulize(
        &molecule,
        &atoms_to_use,
        &bonds_to_use,
        true,
        true,
        None,
        None,
    )
    .unwrap();

    assert!(atoms[center.index()].is_ring_stereo_atom);
    assert!(atoms[left.index()].has_ring_nbr);
    assert_eq!(atoms[center.index()].which_stereo_group, 1);
    assert_eq!(
        atoms[center.index()].type_of_stereo_group,
        CanonStereoGroupType::Or
    );
}

#[test]
fn init_canon_atoms_stereo_group_assignment_marks_group_atoms_outside_fragment_scope() {
    let mut builder = MoleculeBuilder::new();
    let a0 = builder.add_atom(AtomSpec::new(Element::C));
    let a1 = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_stereo_group(StereoGroup::new(
            StereoGroupKind::Or,
            vec![a0, a1],
            Vec::new(),
        ))
        .unwrap();
    let molecule = builder.build().unwrap();
    let atoms_to_use = vec![true, false];
    let bonds_to_use = vec![true; molecule.num_bonds()];

    let atoms = init_fragment_canon_atoms_for_kekulize(
        &molecule,
        &atoms_to_use,
        &bonds_to_use,
        true,
        true,
        None,
        None,
    )
    .unwrap();

    assert_eq!(atoms[a0.index()].which_stereo_group, 1);
    assert_eq!(atoms[a1.index()].which_stereo_group, 1);
    assert_eq!(
        atoms[a1.index()].type_of_stereo_group,
        CanonStereoGroupType::Or
    );
}

#[test]
fn make_canon_bond_holder_populates_cis_trans_controlling_atoms_like_rdkit() {
    let mut builder = MoleculeBuilder::new();
    let left = builder.add_atom(AtomSpec::new(Element::C));
    let right = builder.add_atom(AtomSpec::new(Element::C));
    let fluorine = builder.add_atom(AtomSpec::new(Element::F));
    let chlorine = builder.add_atom(AtomSpec::new(Element::CL));
    let bromine = builder.add_atom(AtomSpec::new(Element::BR));
    let iodine = builder.add_atom(AtomSpec::new(Element::I));
    builder
        .add_bond(BondSpec::new(left, fluorine, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(left, chlorine, BondOrder::Single))
        .unwrap();
    let double_bond = builder
        .add_bond(
            BondSpec::new(left, right, BondOrder::Double)
                .with_stereo(BondStereo::Cis)
                .with_stereo_atoms(fluorine, bromine),
        )
        .unwrap();
    builder
        .add_bond(BondSpec::new(right, bromine, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(right, iodine, BondOrder::Single))
        .unwrap();
    let molecule = builder.build().unwrap();

    let holder =
        make_canon_bond_holder(&molecule, double_bond.index(), right.index(), true).unwrap();

    assert_eq!(holder.stype, BondStereo::Cis);
    assert_eq!(
        holder.controlling_atoms,
        [
            Some(fluorine.index()),
            Some(chlorine.index()),
            Some(bromine.index()),
            Some(iodine.index()),
        ]
    );
}

#[test]
fn atropisomer_helper_orders_neighbor_bonds_by_other_atom_index() {
    let mut builder = MoleculeBuilder::new();
    let begin_low = builder.add_atom(AtomSpec::new(Element::C));
    let begin = builder.add_atom(AtomSpec::new(Element::C));
    let end = builder.add_atom(AtomSpec::new(Element::C));
    let end_low = builder.add_atom(AtomSpec::new(Element::C));
    let end_high = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(BondSpec::new(begin_low, begin, BondOrder::Single))
        .unwrap();
    let atrop_bond = builder
        .add_bond(BondSpec::new(begin, end, BondOrder::Single).with_stereo(BondStereo::AtropCw))
        .unwrap();
    let high_bond = builder
        .add_bond(BondSpec::new(end, end_high, BondOrder::Single))
        .unwrap();
    let low_bond = builder
        .add_bond(BondSpec::new(end, end_low, BondOrder::Single))
        .unwrap();
    let molecule = builder.build().unwrap();

    let atoms_and_bonds =
        atropisomer_atoms_and_bonds_for_canonical_rank(&molecule, atrop_bond).unwrap();

    assert_eq!(atoms_and_bonds[0].0, begin);
    assert_eq!(atoms_and_bonds[1].0, end);
    assert_eq!(atoms_and_bonds[1].1, vec![low_bond, high_bond]);
}

#[test]
fn bondholder_compare_stereo_flips_cis_trans_by_controlling_atom_rank() {
    let left = CanonBondHolder {
        bond_type: BondOrder::Single,
        bond_stereo: rdkit_bond_stereo_rank(BondStereo::Cis),
        stype: BondStereo::Cis,
        controlling_atoms: [Some(0), Some(1), Some(2), None],
        nbr_sym_class: 4,
        nbr_idx: 4,
        p_symbol: None,
        bond_idx: 0,
    };
    let right = CanonBondHolder {
        bond_type: BondOrder::Single,
        bond_stereo: rdkit_bond_stereo_rank(BondStereo::Cis),
        stype: BondStereo::Cis,
        controlling_atoms: [Some(0), None, Some(2), None],
        nbr_sym_class: 4,
        nbr_idx: 4,
        p_symbol: None,
        bond_idx: 1,
    };
    let atom_ranks = vec![0, 3, 1];

    let ordering = compare_canon_bond_holder(&left, &right, &atom_ranks);

    assert_eq!(ordering, Ordering::Greater);
}

#[test]
fn bond_symbol_override_precedes_bond_order_in_bondholder_comparison() {
    let left = CanonBondHolder {
        bond_type: BondOrder::Triple,
        bond_stereo: rdkit_bond_stereo_rank(BondStereo::None),
        stype: BondStereo::None,
        controlling_atoms: [None; 4],
        nbr_sym_class: 0,
        nbr_idx: 0,
        p_symbol: Some("A"),
        bond_idx: 0,
    };
    let right = CanonBondHolder {
        bond_type: BondOrder::Single,
        bond_stereo: rdkit_bond_stereo_rank(BondStereo::None),
        stype: BondStereo::None,
        controlling_atoms: [None; 4],
        nbr_sym_class: 0,
        nbr_idx: 0,
        p_symbol: Some("B"),
        bond_idx: 1,
    };

    let ordering = compare_canon_bond_holder(&left, &right, &[]);

    assert_eq!(ordering, Ordering::Less);
}

#[test]
fn cis_trans_bond_without_stereo_atoms_is_build_error() {
    let mut builder = MoleculeBuilder::new();
    let left = builder.add_atom(AtomSpec::new(Element::C));
    let right = builder.add_atom(AtomSpec::new(Element::C));
    let error = builder
        .add_bond(BondSpec::new(left, right, BondOrder::Double).with_stereo(BondStereo::Cis))
        .unwrap_err();

    assert!(matches!(
        error,
        crate::MoleculeBuildError::BondStereoAtomsRequired {
            stereo: BondStereo::Cis
        }
    ));
}

#[test]
#[ignore = "debug helper for RDKit row 142 canon_atom state parity"]
fn debug_row_142_canon_atom_state() {
    let input = "[C:12]12([CH:62]([CH3:65])[c:61]3[cH:64][cH:67][cH:68][cH:66][cH:63]3)[CH:20]4[c:30]5[c:40]6[c:49]7[c:57]8[c:60]([c:59]9[c:55]([c:47]([c:44]([c:52]9[c:51]([c:43]%10[c:35]%11[c:25]%12[c:19]%13%14)[c:53]8[c:45]%11[c:39]6[c:29]4%13)[c:34]([c:24]%15[c:15]%16[c:7]%17[c:3]%18%19)[c:33]%10[c:23]%16[c:16]%12[c:8]%18[c:11]%14[c:5]1%20)[c:37]([c:36]%21[c:26]%22[c:18]%23[c:10]%24[c:13]%25[c:6]%26%27)[c:27]%15[c:17]%22[c:9]%17[c:4]%24[c:1]%19[c:2]%20%26)[c:54]([c:46]%21[c:38]%28[c:28]%23[c:21]%25%29)[c:56]%30[c:48]%28[c:41]%31[c:31]%29[c:22]%32[c:14]2%27)[c:58]%30[c:50]7[c:42]%31[c:32]5%32";
    let mut molecule = Molecule::from_smiles(input).unwrap();
    let params = crate::notation::smiles_write::SmilesWriteParams {
        do_isomeric_smiles: false,
        do_kekule: true,
        canonical: true,
        clean_stereo: false,
        include_dative_bonds: false,
        ignore_atom_map_numbers: true,
        rooted_at_atom: Some(molecule.num_atoms() - 1),
        ..Default::default()
    };
    if params.ignore_atom_map_numbers {
        for atom in &mut molecule.topology_block_mut().atoms {
            atom.set_atom_map(None);
        }
    }
    let view = CanonRankReadView::from_molecule(&molecule).unwrap();
    let atoms_to_use = vec![true; molecule.num_atoms()];
    let bonds_to_use = vec![true; molecule.num_bonds()];
    let atoms = super::init_fragment_canon_atoms_for_kekulize(
        &view,
        &atoms_to_use,
        &bonds_to_use,
        false,
        false,
        None,
        None,
    )
    .unwrap();
    for (idx, atom) in atoms.iter().enumerate() {
        eprintln!(
            "atom {idx} degree={} totalHs={} isRingAtom={} hasRingNbr={} isRingStereoAtom={} nbrs={:?}",
            atom.degree,
            atom.total_num_hs,
            atom.is_ring_atom,
            atom.has_ring_nbr,
            atom.is_ring_stereo_atom,
            atom.nbr_ids
        );
    }
}

#[test]
#[ignore = "debug helper for RDKit row 142 fragment rank parity"]
fn debug_row_142_fragment_ranks() {
    let input = "[C:12]12([CH:62]([CH3:65])[c:61]3[cH:64][cH:67][cH:68][cH:66][cH:63]3)[CH:20]4[c:30]5[c:40]6[c:49]7[c:57]8[c:60]([c:59]9[c:55]([c:47]([c:44]([c:52]9[c:51]([c:43]%10[c:35]%11[c:25]%12[c:19]%13%14)[c:53]8[c:45]%11[c:39]6[c:29]4%13)[c:34]([c:24]%15[c:15]%16[c:7]%17[c:3]%18%19)[c:33]%10[c:23]%16[c:16]%12[c:8]%18[c:11]%14[c:5]1%20)[c:37]([c:36]%21[c:26]%22[c:18]%23[c:10]%24[c:13]%25[c:6]%26%27)[c:27]%15[c:17]%22[c:9]%17[c:4]%24[c:1]%19[c:2]%20%26)[c:54]([c:46]%21[c:38]%28[c:28]%23[c:21]%25%29)[c:56]%30[c:48]%28[c:41]%31[c:31]%29[c:22]%32[c:14]2%27)[c:58]%30[c:50]7[c:42]%31[c:32]5%32";
    let mut molecule = Molecule::from_smiles(input).unwrap();
    for atom in &mut molecule.topology_block_mut().atoms {
        atom.set_atom_map(None);
    }
    let atoms_to_use = vec![true; molecule.num_atoms()];
    let bonds_to_use = vec![true; molecule.num_bonds()];
    let ranks = super::rank_fragment_atoms(
        &molecule,
        FragmentRankScope::new(&atoms_to_use, &bonds_to_use),
        CanonicalRankOptions::kekulize_default(),
    )
    .unwrap();
    eprintln!("fragment-ranks={ranks:?}");
}
