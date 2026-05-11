use std::collections::HashSet;

use crate::{BondOrder, BondStereo, ChiralTag, Molecule, ValenceModel, assign_valence};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MorganFingerprintParams {
    pub radius: u32,
    pub n_bits: usize,
    pub use_chirality: bool,
    pub use_bond_types: bool,
    pub count_simulation: bool,
    pub count_bounds: Vec<u32>,
    pub only_nonzero_invariants: bool,
    pub include_ring_membership: bool,
    pub include_redundant_environments: bool,
    pub from_atoms: Option<Vec<usize>>,
    pub ignore_atoms: Option<Vec<usize>>,
    pub custom_atom_invariants: Option<Vec<u32>>,
    pub custom_bond_invariants: Option<Vec<u32>>,
    pub atom_invariants_generator: MorganAtomInvariantsGenerator,
    pub bond_invariants_generator: Option<MorganBondInvariantsGenerator>,
    pub num_bits_per_feature: u32,
    pub collect_additional_output: bool,
}

impl Default for MorganFingerprintParams {
    fn default() -> Self {
        Self {
            radius: 2,
            n_bits: 2048,
            use_chirality: false,
            use_bond_types: true,
            count_simulation: false,
            count_bounds: vec![1, 2, 4, 8],
            only_nonzero_invariants: false,
            include_ring_membership: true,
            include_redundant_environments: false,
            from_atoms: None,
            ignore_atoms: None,
            custom_atom_invariants: None,
            custom_bond_invariants: None,
            atom_invariants_generator: MorganAtomInvariantsGenerator::Connectivity {
                include_ring_membership: true,
            },
            bond_invariants_generator: None,
            num_bits_per_feature: 1,
            collect_additional_output: false,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MorganAtomInvariantsGenerator {
    Connectivity { include_ring_membership: bool },
    Feature,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MorganBondInvariantsGenerator {
    pub use_bond_types: bool,
    pub use_chirality: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct MorganAdditionalOutput {
    pub atom_counts: Vec<u32>,
    pub atom_to_bits: Vec<Vec<usize>>,
    pub bit_info_map: std::collections::BTreeMap<usize, Vec<(usize, u32)>>,
    pub atoms_per_bit: std::collections::BTreeMap<usize, Vec<Vec<usize>>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MorganFingerprintOutput {
    pub fingerprint: Fingerprint,
    pub additional_output: Option<MorganAdditionalOutput>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Fingerprint {
    bits: Vec<u64>,
    n_bits: usize,
}

impl Fingerprint {
    #[must_use]
    pub fn from_on_bits(n_bits: usize, on_bits: impl IntoIterator<Item = usize>) -> Self {
        let mut bits = vec![0; n_bits.div_ceil(64)];
        for bit in on_bits {
            assert!(
                bit < n_bits,
                "fingerprint bit {bit} is outside n_bits={n_bits}"
            );
            bits[bit / 64] |= 1u64 << (bit % 64);
        }
        Self { bits, n_bits }
    }

    #[must_use]
    pub fn n_bits(&self) -> usize {
        self.n_bits
    }

    #[must_use]
    pub fn on_bits(&self) -> Vec<usize> {
        let mut out = Vec::new();
        for (word_idx, word) in self.bits.iter().copied().enumerate() {
            let mut remaining = word;
            while remaining != 0 {
                let offset = remaining.trailing_zeros() as usize;
                let bit = word_idx * 64 + offset;
                if bit < self.n_bits {
                    out.push(bit);
                }
                remaining &= remaining - 1;
            }
        }
        out
    }

    pub fn tanimoto(&self, other: &Self) -> Result<f64, FingerprintError> {
        if self.n_bits != other.n_bits {
            return Err(FingerprintError::BitLengthMismatch {
                left: self.n_bits,
                right: other.n_bits,
            });
        }
        let mut intersection = 0u32;
        let mut union = 0u32;
        for (left, right) in self.bits.iter().zip(&other.bits) {
            intersection += (left & right).count_ones();
            union += (left | right).count_ones();
        }
        if union == 0 {
            Ok(0.0)
        } else {
            Ok(f64::from(intersection) / f64::from(union))
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum FingerprintError {
    #[error("Morgan fingerprint requires n_bits > 0")]
    EmptyFingerprint,
    #[error("Morgan num_bits_per_feature must be > 0")]
    EmptyNumBitsPerFeature,
    #[error("Morgan count bounds cannot be empty when count simulation is enabled")]
    EmptyCountBounds,
    #[error("Morgan count bounds size must be smaller than n_bits: {bounds_len} >= {n_bits}")]
    CountBoundsTooLarge { bounds_len: usize, n_bits: usize },
    #[error("custom atom invariants length mismatch: {actual} != {expected}")]
    CustomAtomInvariantLength { actual: usize, expected: usize },
    #[error("custom bond invariants length mismatch: {actual} != {expected}")]
    CustomBondInvariantLength { actual: usize, expected: usize },
    #[error("atom index {atom_index} is outside atom count {atom_count}")]
    AtomIndexOutOfRange {
        atom_index: usize,
        atom_count: usize,
    },
    #[error("Morgan fingerprint valence assignment failed: {0}")]
    Valence(String),
    #[error("fingerprint bit length mismatch: {left} != {right}")]
    BitLengthMismatch { left: usize, right: usize },
}

pub fn morgan_fingerprint(
    mol: &Molecule,
    params: &MorganFingerprintParams,
) -> Result<Fingerprint, FingerprintError> {
    Ok(morgan_fingerprint_with_output(mol, params)?.fingerprint)
}

pub fn morgan_fingerprint_with_output(
    mol: &Molecule,
    params: &MorganFingerprintParams,
) -> Result<MorganFingerprintOutput, FingerprintError> {
    if params.n_bits == 0 {
        return Err(FingerprintError::EmptyFingerprint);
    }
    if params.num_bits_per_feature == 0 {
        return Err(FingerprintError::EmptyNumBitsPerFeature);
    }
    if params.count_simulation {
        if params.count_bounds.is_empty() {
            return Err(FingerprintError::EmptyCountBounds);
        }
        if params.count_bounds.len() >= params.n_bits {
            return Err(FingerprintError::CountBoundsTooLarge {
                bounds_len: params.count_bounds.len(),
                n_bits: params.n_bits,
            });
        }
    }
    validate_atom_indices(mol, params.from_atoms.as_deref())?;
    validate_atom_indices(mol, params.ignore_atoms.as_deref())?;
    let atom_invariants = if let Some(custom) = params.custom_atom_invariants.as_ref() {
        if custom.len() != mol.atoms().len() {
            return Err(FingerprintError::CustomAtomInvariantLength {
                actual: custom.len(),
                expected: mol.atoms().len(),
            });
        }
        custom.clone()
    } else {
        atom_invariants_from_generator(mol, params)?
    };
    let bond_invariants = if let Some(custom) = params.custom_bond_invariants.as_ref() {
        if custom.len() != mol.bonds().len() {
            return Err(FingerprintError::CustomBondInvariantLength {
                actual: custom.len(),
                expected: mol.bonds().len(),
            });
        }
        custom.clone()
    } else {
        bond_invariants(mol, params)
    };
    let environments = morgan_environments(mol, params, &atom_invariants, &bond_invariants)?;
    Ok(fingerprint_from_environments(
        environments,
        params,
        mol.atoms().len(),
    ))
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct MorganAtomEnv {
    code: u32,
    atom_idx: usize,
    layer: u32,
    atoms_involved: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct BondNeighborhood(Vec<bool>);

impl PartialOrd for BondNeighborhood {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for BondNeighborhood {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.iter().rev().cmp(other.0.iter().rev())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct AccumTuple {
    neighborhood: BondNeighborhood,
    code: u32,
    atom_idx: usize,
}

fn validate_atom_indices(
    mol: &Molecule,
    atom_indices: Option<&[usize]>,
) -> Result<(), FingerprintError> {
    if let Some(atom_indices) = atom_indices {
        for &atom_index in atom_indices {
            if atom_index >= mol.atoms().len() {
                return Err(FingerprintError::AtomIndexOutOfRange {
                    atom_index,
                    atom_count: mol.atoms().len(),
                });
            }
        }
    }
    Ok(())
}

fn fingerprint_from_environments(
    environments: Vec<MorganAtomEnv>,
    params: &MorganFingerprintParams,
    atom_count: usize,
) -> MorganFingerprintOutput {
    let effective_size = if params.count_simulation {
        params.n_bits / params.count_bounds.len()
    } else {
        params.n_bits
    };
    let mut counts = std::collections::BTreeMap::<usize, u32>::new();
    let mut additional_output = params
        .collect_additional_output
        .then(|| MorganAdditionalOutput {
            atom_counts: vec![0; atom_count],
            atom_to_bits: vec![Vec::new(); atom_count],
            bit_info_map: std::collections::BTreeMap::new(),
            atoms_per_bit: std::collections::BTreeMap::new(),
        });
    for env in environments {
        let mut bit_id = env.code as usize % effective_size;
        *counts.entry(bit_id).or_insert(0) += 1;
        update_additional_output(&mut additional_output, &env, bit_id);
        if params.num_bits_per_feature > 1 {
            let mut generator = RdkitFingerprintRng::new(env.code);
            for _ in 1..params.num_bits_per_feature {
                bit_id = generator.next_int_max() as usize % effective_size;
                *counts.entry(bit_id).or_insert(0) += 1;
                update_additional_output(&mut additional_output, &env, bit_id);
            }
        }
    }
    let fingerprint = if params.count_simulation {
        let mut on_bits = Vec::new();
        for (bit_id, count) in counts {
            for (bound_idx, &bound) in params.count_bounds.iter().enumerate() {
                if count >= bound {
                    on_bits.push(bit_id * params.count_bounds.len() + bound_idx);
                }
            }
        }
        Fingerprint::from_on_bits(params.n_bits, on_bits)
    } else {
        Fingerprint::from_on_bits(params.n_bits, counts.into_keys())
    };
    MorganFingerprintOutput {
        fingerprint,
        additional_output,
    }
}

fn update_additional_output(
    additional_output: &mut Option<MorganAdditionalOutput>,
    env: &MorganAtomEnv,
    bit_id: usize,
) {
    if let Some(output) = additional_output {
        if output.atom_counts.len() <= env.atom_idx {
            output.atom_counts.resize(env.atom_idx + 1, 0);
        }
        if output.atom_to_bits.len() <= env.atom_idx {
            output.atom_to_bits.resize(env.atom_idx + 1, Vec::new());
        }
        output.atom_counts[env.atom_idx] += 1;
        output.atom_to_bits[env.atom_idx].push(bit_id);
        output
            .bit_info_map
            .entry(bit_id)
            .or_default()
            .push((env.atom_idx, env.layer));
        output
            .atoms_per_bit
            .entry(bit_id)
            .or_default()
            .push(env.atoms_involved.clone());
    }
}

fn atom_invariants_from_generator(
    mol: &Molecule,
    params: &MorganFingerprintParams,
) -> Result<Vec<u32>, FingerprintError> {
    match &params.atom_invariants_generator {
        MorganAtomInvariantsGenerator::Connectivity {
            include_ring_membership,
        } => connectivity_invariants(mol, *include_ring_membership),
        MorganAtomInvariantsGenerator::Feature => feature_invariants(mol),
    }
}

fn connectivity_invariants(
    mol: &Molecule,
    include_ring_membership: bool,
) -> Result<Vec<u32>, FingerprintError> {
    let valence = assign_valence(mol, ValenceModel::RdkitLike)
        .map_err(|err| FingerprintError::Valence(err.to_string()))?;
    let degree = explicit_degree(mol);
    let h_neighbors = explicit_hydrogen_neighbors(mol);
    let ring_atoms = include_ring_membership.then(|| ring_atom_flags(mol));
    let mut invariants = Vec::with_capacity(mol.atoms().len());

    for atom in mol.atoms() {
        let total_hs = atom.explicit_hydrogens as u32
            + valence.implicit_hydrogens[atom.index] as u32
            + h_neighbors[atom.index] as u32;
        let total_degree = degree[atom.index] as u32
            + atom.explicit_hydrogens as u32
            + valence.implicit_hydrogens[atom.index] as u32;
        let mut components = Vec::with_capacity(6);
        components.push(atom.atomic_num as u32);
        components.push(total_degree);
        components.push(total_hs);
        components.push(atom.formal_charge as i32 as u32);
        components.push(isotope_delta_mass(atom.atomic_num, atom.isotope) as u32);
        if ring_atoms
            .as_ref()
            .is_some_and(|ring_atoms| ring_atoms[atom.index])
        {
            components.push(1);
        }
        invariants.push(hash_vec_u32(&components));
    }
    Ok(invariants)
}

fn bond_invariants(mol: &Molecule, params: &MorganFingerprintParams) -> Vec<u32> {
    let (use_bond_types, use_chirality) =
        if let Some(generator) = params.bond_invariants_generator.as_ref() {
            (generator.use_bond_types, generator.use_chirality)
        } else {
            (params.use_bond_types, params.use_chirality)
        };
    let mut invariants = vec![0; mol.bonds().len()];
    for bond in mol.bonds() {
        let mut invariant = 1i32;
        if use_bond_types {
            if !use_chirality
                || !matches!(bond.order, BondOrder::Double)
                || matches!(bond.stereo, BondStereo::None)
            {
                invariant = rdkit_bond_type_value(bond.order);
            } else {
                let stereo_offset = 100;
                let bond_type_offset = 10;
                invariant = stereo_offset
                    + bond_type_offset * rdkit_bond_type_value(bond.order)
                    + rdkit_bond_stereo_value(bond.stereo);
            }
        }
        invariants[bond.index] = invariant as u32;
    }
    invariants
}

fn morgan_environments(
    mol: &Molecule,
    params: &MorganFingerprintParams,
    atom_invariants: &[u32],
    bond_invariants: &[u32],
) -> Result<Vec<MorganAtomEnv>, FingerprintError> {
    let n_atoms = mol.atoms().len();
    let adjacency = crate::AdjacencyList::from_topology(n_atoms, mol.bonds());
    let max_num_results = (params.radius as usize + 1) * n_atoms;
    let mut result = Vec::with_capacity(max_num_results);
    let stereo_props = if params.use_chirality {
        Some(mol.rdkit_legacy_stereo_atom_props(false))
    } else {
        None
    };

    let mut current_invariants = atom_invariants.to_vec();
    let mut next_layer_invariants = vec![0u32; n_atoms];
    let atom_order = if params.only_nonzero_invariants {
        let mut ordering = Vec::<(i32, usize)>::with_capacity(n_atoms);
        for (atom_idx, &invariant) in current_invariants.iter().enumerate() {
            if invariant == 0 {
                ordering.push((1, atom_idx));
            } else {
                ordering.push((0, atom_idx));
            }
        }
        ordering.sort_unstable();
        ordering
            .into_iter()
            .map(|(_, atom_idx)| atom_idx)
            .collect::<Vec<_>>()
    } else {
        (0..n_atoms).collect::<Vec<_>>()
    };
    let include_atoms = include_atom_flags(n_atoms, params.from_atoms.as_deref());

    for (atom_idx, &invariant) in current_invariants.iter().enumerate() {
        if include_atoms[atom_idx] && (!params.only_nonzero_invariants || invariant != 0) {
            result.push(MorganAtomEnv {
                code: invariant,
                atom_idx,
                layer: 0,
                atoms_involved: vec![atom_idx],
            });
        }
    }

    let mut neighborhoods = HashSet::<BondNeighborhood>::with_capacity(max_num_results);
    let mut atom_neighborhoods = vec![BondNeighborhood(vec![false; mol.bonds().len()]); n_atoms];
    let mut round_atom_neighborhoods = atom_neighborhoods.clone();
    let mut dead_atoms = vec![false; n_atoms];
    let mut chiral_atoms = vec![false; n_atoms];
    let mut neighborhood_invariants = Vec::<(i32, u32)>::with_capacity(8);

    for layer in 0..params.radius {
        let mut all_neighborhoods_this_round = Vec::<AccumTuple>::new();
        for &atom_idx in &atom_order {
            if dead_atoms[atom_idx] {
                continue;
            }
            let neighbors = adjacency.neighbors_of(atom_idx);
            if neighbors.is_empty() {
                dead_atoms[atom_idx] = true;
                continue;
            }

            neighborhood_invariants.clear();
            for neighbor in neighbors {
                let bond_idx = neighbor.bond_index;
                round_atom_neighborhoods[atom_idx].0[bond_idx] = true;
                let other_idx = neighbor.atom_index;
                for (idx, in_neighborhood) in atom_neighborhoods[other_idx].0.iter().enumerate() {
                    if *in_neighborhood {
                        round_atom_neighborhoods[atom_idx].0[idx] = true;
                    }
                }
                neighborhood_invariants.push((
                    bond_invariants[bond_idx] as i32,
                    current_invariants[other_idx],
                ));
            }

            neighborhood_invariants.sort_unstable();
            let mut invariant = layer;
            hash_combine_u32(&mut invariant, current_invariants[atom_idx]);
            let mut looks_chiral = mol.atoms()[atom_idx].chiral_tag != ChiralTag::Unspecified;
            for (neighbor_idx, neighbor_invariant) in neighborhood_invariants.iter().enumerate() {
                hash_combine_pair_i32_u32(&mut invariant, *neighbor_invariant);
                if params.use_chirality && looks_chiral && !chiral_atoms[atom_idx] {
                    let duplicated_neighbor_invariant = neighbor_idx > 0
                        && neighbor_invariant.1 == neighborhood_invariants[neighbor_idx - 1].1;
                    if neighbor_invariant.0 != rdkit_bond_type_value(BondOrder::Single)
                        || duplicated_neighbor_invariant
                    {
                        looks_chiral = false;
                    }
                }
            }
            if params.use_chirality && looks_chiral {
                chiral_atoms[atom_idx] = true;
                let cip = stereo_props
                    .as_ref()
                    .and_then(|props| props.get(atom_idx))
                    .and_then(|props| props.cip_code.as_deref())
                    .unwrap_or("");
                match cip {
                    "R" => hash_combine_u32(&mut invariant, 3),
                    "S" => hash_combine_u32(&mut invariant, 2),
                    _ => hash_combine_u32(&mut invariant, 1),
                }
            }
            next_layer_invariants[atom_idx] = invariant;
            all_neighborhoods_this_round.push(AccumTuple {
                neighborhood: round_atom_neighborhoods[atom_idx].clone(),
                code: invariant,
                atom_idx,
            });
        }

        all_neighborhoods_this_round.sort();
        for item in all_neighborhoods_this_round {
            if params.include_redundant_environments || !neighborhoods.contains(&item.neighborhood)
            {
                if (!params.only_nonzero_invariants || atom_invariants[item.atom_idx] != 0)
                    && include_atoms[item.atom_idx]
                {
                    result.push(MorganAtomEnv {
                        code: item.code,
                        atom_idx: item.atom_idx,
                        layer: layer + 1,
                        atoms_involved: atoms_within_radius(&adjacency, item.atom_idx, layer + 1),
                    });
                    neighborhoods.insert(item.neighborhood);
                }
            } else {
                dead_atoms[item.atom_idx] = true;
            }
        }

        std::mem::swap(&mut current_invariants, &mut next_layer_invariants);
        next_layer_invariants.fill(0);
        atom_neighborhoods = round_atom_neighborhoods.clone();
    }

    Ok(result)
}

fn include_atom_flags(n_atoms: usize, from_atoms: Option<&[usize]>) -> Vec<bool> {
    let mut include_atoms = vec![from_atoms.is_none(); n_atoms];
    if let Some(from_atoms) = from_atoms {
        for &atom_idx in from_atoms {
            include_atoms[atom_idx] = true;
        }
    }
    include_atoms
}

fn atoms_within_radius(adjacency: &crate::AdjacencyList, center: usize, radius: u32) -> Vec<usize> {
    let mut out = vec![center];
    if radius == 0 {
        return out;
    }
    let mut seen = vec![false; adjacency.neighbors.len()];
    let mut distances = vec![u32::MAX; adjacency.neighbors.len()];
    let mut queue = std::collections::VecDeque::<(usize, u32)>::new();
    seen[center] = true;
    distances[center] = 0;
    queue.push_back((center, 0));
    while let Some((atom_idx, dist)) = queue.pop_front() {
        if dist >= radius {
            continue;
        }
        for neighbor in adjacency.neighbors_of(atom_idx) {
            if seen[neighbor.atom_index] {
                continue;
            }
            seen[neighbor.atom_index] = true;
            let next_dist = dist + 1;
            distances[neighbor.atom_index] = next_dist;
            if next_dist <= radius {
                queue.push_back((neighbor.atom_index, next_dist));
            }
        }
    }
    for (atom_idx, &distance) in distances.iter().enumerate() {
        if atom_idx != center && distance <= radius {
            out.push(atom_idx);
        }
    }
    out
}

fn explicit_degree(mol: &Molecule) -> Vec<usize> {
    let mut degree = vec![0usize; mol.atoms().len()];
    for bond in mol.bonds() {
        degree[bond.begin_atom] += 1;
        degree[bond.end_atom] += 1;
    }
    degree
}

fn explicit_hydrogen_neighbors(mol: &Molecule) -> Vec<u8> {
    let mut count = vec![0u8; mol.atoms().len()];
    for bond in mol.bonds() {
        if mol.atoms()[bond.begin_atom].atomic_num == 1 {
            count[bond.end_atom] = count[bond.end_atom].saturating_add(1);
        }
        if mol.atoms()[bond.end_atom].atomic_num == 1 {
            count[bond.begin_atom] = count[bond.begin_atom].saturating_add(1);
        }
    }
    count
}

fn feature_invariants(mol: &Molecule) -> Result<Vec<u32>, FingerprintError> {
    let valence = assign_valence(mol, ValenceModel::RdkitLike)
        .map_err(|err| FingerprintError::Valence(err.to_string()))?;
    let adjacency = crate::AdjacencyList::from_topology(mol.atoms().len(), mol.bonds());
    let degree = explicit_degree(mol);
    let h_neighbors = explicit_hydrogen_neighbors(mol);
    let total_hs_by_atom = mol
        .atoms()
        .iter()
        .map(|atom| {
            atom.explicit_hydrogens as usize
                + h_neighbors[atom.index] as usize
                + valence.implicit_hydrogens[atom.index] as usize
        })
        .collect::<Vec<_>>();
    Ok(mol
        .atoms()
        .iter()
        .map(|atom| {
            let total_hs = total_hs_by_atom[atom.index];
            let explicit_valence = explicit_valence_rough(mol, atom.index)
                + valence.implicit_hydrogens[atom.index] as u32;
            let mut invariant = 0u32;
            if is_feature_donor(
                atom.atomic_num,
                total_hs,
                explicit_valence,
                atom.formal_charge,
            ) {
                invariant |= 1;
            }
            if is_feature_acceptor(mol, &adjacency, atom.index, total_hs, explicit_valence) {
                invariant |= 2;
            }
            if atom.is_aromatic {
                invariant |= 4;
            }
            if matches!(atom.atomic_num, 9 | 17 | 35 | 53) {
                invariant |= 8;
            }
            if is_feature_basic(mol, &adjacency, atom.index, total_hs, degree[atom.index]) {
                invariant |= 16;
            }
            if is_feature_acidic(mol, atom.index, &total_hs_by_atom) {
                invariant |= 32;
            }
            invariant
        })
        .collect())
}

fn explicit_valence_rough(mol: &Molecule, atom_index: usize) -> u32 {
    let atom = &mol.atoms()[atom_index];
    let mut valence = atom.explicit_hydrogens as u32;
    for bond in mol.bonds() {
        if bond.begin_atom == atom_index || bond.end_atom == atom_index {
            if matches!(bond.order, BondOrder::Dative) && bond.begin_atom == atom_index {
                continue;
            }
            valence += match bond.order {
                BondOrder::Single | BondOrder::Aromatic | BondOrder::Dative => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                BondOrder::Quadruple => 4,
                BondOrder::Hydrogen => 0,
                BondOrder::Null => 0,
            };
        }
    }
    valence
}

fn is_feature_donor(
    atomic_num: u8,
    total_hs: usize,
    explicit_valence: u32,
    formal_charge: i8,
) -> bool {
    (atomic_num == 7
        && total_hs > 0
        && matches!(explicit_valence, 3 | 4)
        && matches!(formal_charge, 0 | 1))
        || ((atomic_num == 8 || atomic_num == 16) && total_hs == 1 && formal_charge == 0)
}

fn is_feature_acceptor(
    mol: &Molecule,
    adjacency: &crate::AdjacencyList,
    atom_index: usize,
    total_hs: usize,
    explicit_valence: u32,
) -> bool {
    let atom = &mol.atoms()[atom_index];
    if matches!(atom.atomic_num, 8 | 16) && !atom.is_aromatic {
        if atom.formal_charge < 0 {
            return true;
        }
        if total_hs == 0 && explicit_valence == 2 {
            return true;
        }
        if total_hs == 1
            && explicit_valence == 2
            && !neighbor_has_double_bond_to_hetero(mol, adjacency, atom_index)
        {
            return true;
        }
    }
    if atom.atomic_num == 7
        && !atom.is_aromatic
        && explicit_valence == 3
        && !neighbor_has_double_bond_to_hetero(mol, adjacency, atom_index)
    {
        return true;
    }
    if atom.is_aromatic && atom.atomic_num == 7 && total_hs == 0 && atom.formal_charge == 0 {
        return true;
    }
    if atom.is_aromatic
        && matches!(atom.atomic_num, 8 | 16)
        && atom.formal_charge == 0
        && !aromatic_o_or_s_excluded_from_acceptor(mol, adjacency, atom_index)
    {
        return true;
    }
    false
}

fn aromatic_o_or_s_excluded_from_acceptor(
    mol: &Molecule,
    adjacency: &crate::AdjacencyList,
    atom_index: usize,
) -> bool {
    adjacency.neighbors_of(atom_index).iter().any(|neighbor| {
        let neighbor_atom = &mol.atoms()[neighbor.atom_index];
        if neighbor_atom.is_aromatic && neighbor_atom.atomic_num == 7 {
            return true;
        }
        neighbor_atom.is_aromatic
            && neighbor_atom.atomic_num == 6
            && adjacency
                .neighbors_of(neighbor.atom_index)
                .iter()
                .any(|second| {
                    second.atom_index != atom_index
                        && mol.atoms()[second.atom_index].is_aromatic
                        && mol.atoms()[second.atom_index].atomic_num == 7
                })
    })
}

fn is_feature_basic(
    mol: &Molecule,
    adjacency: &crate::AdjacencyList,
    atom_index: usize,
    total_hs: usize,
    degree: usize,
) -> bool {
    let atom = &mol.atoms()[atom_index];
    if atom.atomic_num != 7 {
        return false;
    }
    if atom.formal_charge > 0 {
        return true;
    }
    if atom.formal_charge != 0 || total_hs + degree == 0 {
        return false;
    }
    if mol.bonds().iter().any(|bond| {
        (bond.begin_atom == atom_index || bond.end_atom == atom_index)
            && !matches!(bond.order, BondOrder::Single)
    }) {
        return false;
    }
    if adjacency.neighbors_of(atom_index).iter().any(|neighbor| {
        let neighbor_atom = &mol.atoms()[neighbor.atom_index];
        if total_hs == 0 {
            neighbor_atom.atomic_num != 6 || neighbor_atom.is_aromatic
        } else {
            neighbor_atom.atomic_num != 6 && !neighbor_atom.is_aromatic
        }
    }) {
        return false;
    }
    !neighbor_has_carbonyl_carbon(mol, adjacency, atom_index)
}

fn is_feature_acidic(mol: &Molecule, atom_index: usize, total_hs_by_atom: &[usize]) -> bool {
    let atom = &mol.atoms()[atom_index];
    matches!(atom.atomic_num, 6 | 16)
        && mol.bonds().iter().any(|bond| {
            (bond.begin_atom == atom_index || bond.end_atom == atom_index)
                && matches!(bond.order, BondOrder::Double)
                && matches!(
                    mol.atoms()[other_atom_idx(bond, atom_index)].atomic_num,
                    8 | 16 | 15
                )
        })
        && mol.bonds().iter().any(|bond| {
            (bond.begin_atom == atom_index || bond.end_atom == atom_index)
                && matches!(bond.order, BondOrder::Single)
                && {
                    let other = &mol.atoms()[other_atom_idx(bond, atom_index)];
                    other.atomic_num == 8
                        && (total_hs_by_atom[other.index] > 0 || other.formal_charge < 0)
                }
        })
}

fn neighbor_has_double_bond_to_hetero(
    mol: &Molecule,
    adjacency: &crate::AdjacencyList,
    atom_index: usize,
) -> bool {
    adjacency.neighbors_of(atom_index).iter().any(|neighbor| {
        if !matches!(mol.bonds()[neighbor.bond_index].order, BondOrder::Single) {
            return false;
        }
        mol.bonds().iter().any(|bond| {
            (bond.begin_atom == neighbor.atom_index || bond.end_atom == neighbor.atom_index)
                && matches!(bond.order, BondOrder::Double)
                && matches!(
                    mol.atoms()[other_atom_idx(bond, neighbor.atom_index)].atomic_num,
                    7 | 8 | 15 | 16
                )
        })
    })
}

fn neighbor_has_carbonyl_carbon(
    mol: &Molecule,
    adjacency: &crate::AdjacencyList,
    atom_index: usize,
) -> bool {
    adjacency.neighbors_of(atom_index).iter().any(|neighbor| {
        let atom = &mol.atoms()[neighbor.atom_index];
        atom.atomic_num == 6
            && mol.bonds().iter().any(|bond| {
                (bond.begin_atom == neighbor.atom_index || bond.end_atom == neighbor.atom_index)
                    && matches!(bond.order, BondOrder::Double)
                    && mol.atoms()[other_atom_idx(bond, neighbor.atom_index)].atomic_num == 8
            })
    })
}

fn other_atom_idx(bond: &crate::Bond, atom_index: usize) -> usize {
    if bond.begin_atom == atom_index {
        bond.end_atom
    } else {
        bond.begin_atom
    }
}

fn ring_atom_flags(mol: &Molecule) -> Vec<bool> {
    let mut flags = vec![false; mol.atoms().len()];
    for ring in crate::distgeom::rdkit_atom_rings(mol) {
        for atom_idx in ring {
            flags[atom_idx] = true;
        }
    }
    flags
}

fn rdkit_bond_type_value(order: BondOrder) -> i32 {
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => 12,
        BondOrder::Dative => 17,
        BondOrder::Hydrogen => 0,
        BondOrder::Null => 0,
    }
}

fn rdkit_bond_stereo_value(stereo: BondStereo) -> i32 {
    match stereo {
        BondStereo::None => 0,
        BondStereo::Any => 1,
        BondStereo::Cis => 2,
        BondStereo::Trans => 3,
    }
}

fn isotope_delta_mass(atomic_num: u8, isotope: Option<u16>) -> i32 {
    let Some(isotope) = isotope else {
        return 0;
    };
    let Some(standard_weight) = crate::periodic_table::average_atomic_weight(atomic_num) else {
        return 0;
    };
    let isotope_mass = crate::periodic_table::exact_isotope_mass(atomic_num, isotope)
        .unwrap_or(f64::from(isotope));
    (isotope_mass - standard_weight) as i32
}

fn hash_vec_u32(values: &[u32]) -> u32 {
    let mut seed = 0u32;
    for &value in values {
        hash_combine_u32(&mut seed, value);
    }
    seed
}

fn hash_pair_i32_u32(value: (i32, u32)) -> u32 {
    let mut seed = 0u32;
    hash_combine_i32(&mut seed, value.0);
    hash_combine_u32(&mut seed, value.1);
    seed
}

fn hash_combine_pair_i32_u32(seed: &mut u32, value: (i32, u32)) {
    hash_combine_raw(seed, hash_pair_i32_u32(value));
}

fn hash_combine_i32(seed: &mut u32, value: i32) {
    hash_combine_raw(seed, value as u32);
}

fn hash_combine_u32(seed: &mut u32, value: u32) {
    hash_combine_raw(seed, value);
}

fn hash_combine_raw(seed: &mut u32, hash_value: u32) {
    *seed ^= hash_value
        .wrapping_add(0x9e37_79b9)
        .wrapping_add(seed.wrapping_shl(6))
        .wrapping_add(seed.wrapping_shr(2));
}

struct RdkitFingerprintRng {
    state: [u32; 4],
    index: usize,
}

impl RdkitFingerprintRng {
    fn new(seed: u32) -> Self {
        let mut state = [0u32; 4];
        state[0] = seed;
        for i in 1..4 {
            let prev = state[i - 1];
            state[i] = 1_812_433_253u32
                .wrapping_mul(prev ^ (prev >> 30))
                .wrapping_add(i as u32);
        }
        Self { state, index: 4 }
    }

    fn next_raw(&mut self) -> u32 {
        if self.index >= 4 {
            self.twist();
        }
        let mut y = self.state[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }

    fn next_int_max(&mut self) -> u32 {
        self.next_raw() >> 1
    }

    fn twist(&mut self) {
        const UPPER_MASK: u32 = 0x8000_0000;
        const LOWER_MASK: u32 = 0x7fff_ffff;
        for i in 0..4 {
            let y = (self.state[i] & UPPER_MASK) | (self.state[(i + 1) % 4] & LOWER_MASK);
            let mut next = self.state[(i + 2) % 4] ^ (y >> 1);
            if y & 1 != 0 {
                next ^= 0x9908_b0df;
            }
            self.state[i] = next;
        }
        self.index = 0;
    }
}

#[allow(dead_code)]
const RDKIT_MORGAN_SOURCE_CHECKLIST: &str = r#"
RDKit source reference copied from third_party/rdkit Release_2026_03_1.
Legend: [x] line or block reproduced in Rust above; [ ] not yet reproduced.

[x] MorganGenerator.cpp::MorganAtomInvGenerator::getAtomInvariants
    getConnectivityInvariants(mol, *atomInvariants, df_includeRingMembership);

[x] MorganWrapper.cpp::getMorganGenerator
    bool,  // includeRingMembership

[x] FingerprintUtil.cpp::getConnectivityInvariants
    components.push_back(atom->getAtomicNum());
    components.push_back(atom->getTotalDegree());
    components.push_back(atom->getTotalNumHs(true));
    components.push_back(atom->getFormalCharge());
    int deltaMass = static_cast<int>(
        atom->getMass() -
        PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum()));
    components.push_back(deltaMass);
    if (includeRingMembership &&
        atom->getOwningMol().getRingInfo()->numAtomRings(atom->getIdx())) {
      components.push_back(1);
    }
    invars[i] = vectHasher(components);

[x] MorganGenerator.cpp::MorganBondInvGenerator::getBondInvariants
    int32_t bondInvariant = 1;
    if (df_useBondTypes) {
      if (!df_useChirality || bond->getBondType() != Bond::DOUBLE ||
          bond->getStereo() == Bond::STEREONONE) {
        bondInvariant = static_cast<int32_t>(bond->getBondType());
      } else {
        auto bondStereo = static_cast<int32_t>(bond->getStereo());
[ ]     if (!Chirality::getUseLegacyStereoPerception()) { ... CIP E/Z remap ... }
        const int32_t stereoOffset = 100;
        const int32_t bondTypeOffset = 10;
        bondInvariant =
            stereoOffset +
            bondTypeOffset * static_cast<int32_t>(bond->getBondType()) +
            bondStereo;
      }
    }
    (*result)[bond->getIdx()] = static_cast<int32_t>(bondInvariant);

[x] FingerprintGenerator.cpp::getFingerprint
    std::uint32_t effectiveSize = dp_fingerprintArguments->d_fpSize;
[x] if (dp_fingerprintArguments->df_countSimulation) { ... }
    auto tempResult = getFingerprintHelper(mol, args, effectiveSize);
    auto result =
        std::make_unique<ExplicitBitVect>(dp_fingerprintArguments->d_fpSize);
    for (auto val : tempResult->getNonzeroElements()) {
[x]   if (dp_fingerprintArguments->df_countSimulation) { ... }
[x]   else {
        result->setBit(val.first);
      }
    }

[x] FingerprintGenerator.cpp::getFingerprintHelper
    std::unique_ptr<std::vector<std::uint32_t>> atomInvariants = nullptr;
[x] if (args.customAtomInvariants) { ... }
[x] else if (dp_atomInvariantsGenerator) {
      atomInvariants.reset(dp_atomInvariantsGenerator->getAtomInvariants(mol));
    }
    std::unique_ptr<std::vector<std::uint32_t>> bondInvariants = nullptr;
[x] if (args.customBondInvariants) { ... }
[x] else if (dp_bondInvariantsGenerator) {
      bondInvariants.reset(dp_bondInvariantsGenerator->getBondInvariants(mol));
    }
    auto atomEnvironments = dp_atomEnvironmentGenerator->getEnvironments(...);
    auto res = std::make_unique<SparseIntVect<OutputType>>(
        fpSize ? fpSize : dp_atomEnvironmentGenerator->getResultSize());
    for (const auto env : atomEnvironments) {
      OutputType seed = env->getBitId(...);
      auto bitId = seed;
      if (fpSize != 0) {
        bitId %= fpSize;
      }
      res->setVal(bitId, res->getVal(bitId) + 1);
[x]   if (args.additionalOutput) { env->updateAdditionalOutput(...); }
[x]   if (dp_fingerprintArguments->d_numBitsPerFeature > 1) { ... }
    }

[x] MorganGenerator.cpp::MorganEnvGenerator::getEnvironments non-chiral/default path
    currentInvariants = atomInvariants
    fromAtoms includeAtoms filter
    onlyNonzeroInvariants atom ordering and zero-invariant skips
    add round 0 invariants
    for layer in radius:
      collect neighbor (bondInvariant, currentInvariant)
      sort neighbor list
      invar = layer
      hash_combine(invar, currentInvariants[atomIdx])
      hash_combine(invar, each sorted neighbor pair)
      nextLayerInvariants[atomIdx] = invar
      sort neighborhoods this round
      add unseen exact bond neighborhoods; mark duplicates dead
      includeRedundantEnvironments bypasses neighborhood deduplication
      swap current and next invariants

[x] MorganGenerator.cpp::MorganEnvGenerator::getEnvironments chiral atom path
    CIPLabeler::assignCIPLabels(...)
    looksChiral / chiralAtoms checks
    hash_combine(invar, 3 for R / 2 for S / 1 fallback)

[x] MorganGenerator.cpp::MorganEnvGenerator::getEnvironments ignored parameter
    const std::vector<std::uint32_t> *,  // ignoreAtoms

[ ] Unsupported after this pass:
    custom SMARTS pattern list for MorganFeatureAtomInvGenerator,
    AdditionalOutput duplication through countSimulation.
"#;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn boost_hash_matches_rdkit_connectivity_invariant_for_methane() {
        assert_eq!(hash_vec_u32(&[6, 4, 4, 0, 0]), 2_246_733_040);
    }

    #[test]
    fn tanimoto_of_empty_fingerprints_matches_rdkit_explicit_bitvect() {
        let left = Fingerprint::from_on_bits(128, []);
        let right = Fingerprint::from_on_bits(128, []);
        assert_eq!(left.tanimoto(&right).expect("same length"), 0.0);
    }

    #[test]
    fn rdkit_fingerprint_rng_matches_boost_custom_mt19937() {
        let mut rng = RdkitFingerprintRng::new(2_246_733_040);
        assert_eq!(
            (0..5).map(|_| rng.next_int_max()).collect::<Vec<_>>(),
            vec![
                1_322_820_225,
                1_619_309_095,
                1_630_852_618,
                239_107_533,
                2_079_982_906
            ]
        );
    }

    #[test]
    fn source_checklist_tracks_unfinished_rdkit_branches() {
        assert!(RDKIT_MORGAN_SOURCE_CHECKLIST.contains("[ ]"));
        assert!(RDKIT_MORGAN_SOURCE_CHECKLIST.contains("custom SMARTS pattern"));
    }
}
