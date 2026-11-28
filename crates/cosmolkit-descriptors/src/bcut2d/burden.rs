//! Standard BCUT2D projections over parent-prepared atom properties.
//!
//! The matrix and eigensolver are private implementation details. The caller
//! supplies the four source-ordered property rows after preparing the graph.

use super::eigen::{Workspace as SymmetricEigenWorkspace, extrema as symmetric_eigen_extrema};
use crate::DescriptorError;
use cosmolkit_model::{BondOrder, TopologyBlock};

const NONBONDED_WEIGHT: f64 = 0.001;

/// Calculate the four standard BCUT2D high/low eigenvalue pairs.
///
/// `masses`, `charges`, `crippen_logp`, and `crippen_mr` must each contain one
/// value per atom, in the topology's atom order. Chemistry preparation and
/// contribution calculation remain owned by the caller.
pub(super) fn bcut2d(
    topology: &TopologyBlock,
    masses: &[f64],
    charges: &[f64],
    crippen_logp: &[f64],
    crippen_mr: &[f64],
) -> Result<[f64; 8], DescriptorError> {
    let atom_count = topology.atoms.len();

    let rows = [
        ("mass", masses),
        ("charge", charges),
        ("Crippen logP", crippen_logp),
        ("Crippen molar refractivity", crippen_mr),
    ];
    for (weighting, row) in rows {
        if row.len() != atom_count {
            return Err(DescriptorError::BcutDiagonalRows {
                field: weighting,
                actual: row.len(),
                expected: atom_count,
            });
        }
    }

    if atom_count == 0 {
        return Ok([0.0; 8]);
    }

    let mut burden = BurdenMatrix::new(topology)?;
    let mut workspace = SymmetricEigenWorkspace::new(atom_count)?;
    let mut result = [0.0; 8];

    for (row_index, (_, diagonal)) in rows.into_iter().enumerate() {
        burden.replace_diagonal(diagonal);
        let (highest, lowest) = symmetric_eigen_extrema(&burden.values, &mut workspace)?;
        result[row_index * 2] = highest;
        result[row_index * 2 + 1] = lowest;
    }

    Ok(result)
}

struct BurdenMatrix {
    dimension: usize,
    values: Vec<f64>,
}

impl BurdenMatrix {
    fn new(topology: &TopologyBlock) -> Result<Self, DescriptorError> {
        // RDKit 2026.09.1:
        // Code/GraphMol/Descriptors/BCUT.cpp
        // RDKit❗❗: std::unique_ptr<Eigen::MatrixXd> make_burden(const ROMol &m) {
        // RDKit❗❗:   auto num_atoms = m.getNumAtoms();
        // RDKit❗❗:   std::unique_ptr<Eigen::MatrixXd> burden(
        // RDKit❗❗:       new Eigen::MatrixXd(num_atoms, num_atoms));
        // RDKit❗❗:   for (unsigned int i = 0; i < num_atoms; ++i) {
        // RDKit❗❗:     for (unsigned int j = 0; j < num_atoms; ++j) {
        // RDKit❗❗:       (*burden)(i, j) = (*burden)(j, i) = 0.001;
        // RDKit❗❗:     }
        // RDKit❗❗:   }
        let dimension = topology.atoms.len();
        let mut values = allocate_dense_values(dimension, NONBONDED_WEIGHT)?;

        // RDKit❗❗:   for (auto &bond : m.bonds()) {
        // RDKit❗❗:     unsigned int i = bond->getBeginAtomIdx();
        // RDKit❗❗:     unsigned int j = bond->getEndAtomIdx();
        // RDKit❗❗:     double score = 0.0;
        for bond in &topology.bonds {
            let begin = bond.begin().index();
            let end = bond.end().index();
            let weight = bond_weight(bond.order()).ok_or(DescriptorError::UnsupportedBcutBond {
                bond: bond.id(),
                order: bond.order(),
            })?;

            // RDKit❗❗:     (*burden)(i, j) = (*burden)(j, i) = score;
            // RDKit❗❗:   }
            values[begin * dimension + end] = weight;
            values[end * dimension + begin] = weight;
        }

        // RDKit❗❗:   return burden;
        // RDKit❗❗: }
        Ok(Self { dimension, values })
    }

    fn replace_diagonal(&mut self, diagonal: &[f64]) {
        for (index, &value) in diagonal.iter().enumerate() {
            self.values[index * self.dimension + index] = value;
        }
    }
}

fn bond_weight(order: BondOrder) -> Option<f64> {
    // RDKit❗❗:     switch (bond->getBondType()) {
    // RDKit❗❗:       case Bond::AROMATIC:
    // RDKit❗❗:         // score = 0.15; orig burden
    // RDKit❗❗:         score = 0.8164965809277261;  // 1/sqrt(1.5)
    // RDKit❗❗:         break;
    // RDKit❗❗:       case Bond::SINGLE:
    // RDKit❗❗:         // score = 0.1;
    // RDKit❗❗:         score = 1.0;  // 1/sqrt(1.0)
    // RDKit❗❗:         break;
    // RDKit❗❗:       case Bond::DOUBLE:
    // RDKit❗❗:         // score = 0.2;
    // RDKit❗❗:         score = 0.7071067811865475;  // 1/sqrt(2.0)
    // RDKit❗❗:         break;
    // RDKit❗❗:       case Bond::TRIPLE:
    // RDKit❗❗:         // score = 0.3;
    // RDKit❗❗:         score = 0.5773502691896258;  // 1/sqrt(3);
    // RDKit❗❗:         break;
    // RDKit❗❗:       default:
    // RDKit❗❗:         CHECK_INVARIANT(
    // RDKit❗❗:             0, "Bond order must be Single, Double, Triple or Aromatic");
    // RDKit❗❗:     }
    match order {
        BondOrder::Single => Some(1.0),
        BondOrder::Double => Some(0.7071067811865475),
        BondOrder::Triple => Some(0.5773502691896258),
        BondOrder::Aromatic => Some(0.8164965809277261),
        _ => None,
    }
}

pub(super) fn allocate_dense_values(
    dimension: usize,
    fill: f64,
) -> Result<Vec<f64>, DescriptorError> {
    let len = dimension
        .checked_mul(dimension)
        .ok_or(DescriptorError::CountOverflow {
            function: "bcut2d",
            field: "matrix size",
        })?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(len)
        .map_err(|_| DescriptorError::BcutMatrixAllocation { dimension })?;
    values.resize(len, fill);
    Ok(values)
}

pub(super) fn allocate_vector(length: usize) -> Result<Vec<f64>, DescriptorError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(length)
        .map_err(|_| DescriptorError::BcutMatrixAllocation { dimension: length })?;
    values.resize(length, 0.0);
    Ok(values)
}

#[cfg(test)]
mod tests {
    use cosmolkit_model::{
        AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Element,
        TopologyBlock,
    };

    use super::DescriptorError;
    use super::{
        BurdenMatrix, NONBONDED_WEIGHT, SymmetricEigenWorkspace, bcut2d, bond_weight,
        symmetric_eigen_extrema,
    };

    fn graph(
        atom_count: usize,
        bonds: &[(usize, usize, cosmolkit_model::BondOrder)],
    ) -> TopologyBlock {
        let atoms = (0..atom_count)
            .map(|index| {
                Atom::from_spec(
                    AtomId::new(index),
                    AtomSpec::new(Element::C).with_no_implicit(true),
                )
            })
            .collect::<Vec<_>>();
        let bonds = bonds
            .iter()
            .enumerate()
            .map(|(index, &(begin, end, order))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
                )
            })
            .collect::<Vec<_>>();
        let adjacency = AdjacencyList::try_from_topology(atom_count, &bonds).unwrap();
        let topology = TopologyBlock {
            atoms,
            bonds,
            adjacency,
            ..Default::default()
        };
        topology.validate().unwrap();
        topology
    }

    fn parsed(smiles: &str) -> TopologyBlock {
        cosmolkit_smiles::parse_smiles(smiles, &Default::default())
            .expect("valid SMILES")
            .topology
    }

    fn properties(atom_count: usize, value: f64) -> Vec<f64> {
        vec![value; atom_count]
    }

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual:?} differs from {expected:?} by more than {tolerance:?}"
        );
    }

    #[test]
    fn eigen_scaling_householder_and_qr_follow_the_fixed_source_path() {
        // Eigen 3.3.4 MatrixXd / Lower / Packet2d, default solver options.
        // The same matrix is the mass-weighted Burden matrix of CCO.
        let matrix = [12.011, 1.0, 0.001, 1.0, 12.011, 1.0, 0.001, 1.0, 15.999];
        let mut workspace = SymmetricEigenWorkspace::new(3).unwrap();
        let (high, low) = symmetric_eigen_extrema(&matrix, &mut workspace).unwrap();
        assert_eq!(high.to_bits(), 0x4030_3fbd_e550_7372);
        assert_eq!(low.to_bits(), 0x4025_d0fb_6693_efe5);
        // A five-cycle distinguishes the Eigen 3.3.4 QR deflation from 3.4.
        let mut matrix = vec![NONBONDED_WEIGHT; 25];
        for i in 0..5 {
            matrix[i * 5 + i] = 12.011;
            matrix[i * 5 + (i + 1) % 5] = 1.0;
            matrix[((i + 1) % 5) * 5 + i] = 1.0;
        }
        let mut workspace = SymmetricEigenWorkspace::new(5).unwrap();
        let (high, low) = symmetric_eigen_extrema(&matrix, &mut workspace).unwrap();
        assert_eq!(high.to_bits(), 0x402c_06a7_ef9d_b22e);
        assert_eq!(low.to_bits(), 0x4024_c983_d915_0e79);
    }

    #[test]
    fn empty_graph_returns_eight_positive_zeros() {
        let topology = TopologyBlock::default();
        let values = bcut2d(&topology, &[], &[], &[], &[]).unwrap();

        assert_eq!(values, [0.0; 8]);
        assert!(
            values
                .iter()
                .all(|value| value.to_bits() == 0.0_f64.to_bits())
        );
    }

    #[test]
    fn one_atom_duplicates_each_supplied_diagonal() {
        let topology = graph(1, &[]);
        let values = bcut2d(&topology, &[12.011], &[-0.25], &[0.1441], &[2.503]).unwrap();

        assert_eq!(
            values,
            [12.011, 12.011, -0.25, -0.25, 0.1441, 0.1441, 2.503, 2.503]
        );
    }

    #[test]
    fn disconnected_pair_uses_nonbonded_weight_in_both_directions() {
        let topology = graph(2, &[]);
        let matrix = BurdenMatrix::new(&topology).unwrap();
        assert_eq!(matrix.values[1].to_bits(), NONBONDED_WEIGHT.to_bits());
        assert_eq!(matrix.values[2].to_bits(), NONBONDED_WEIGHT.to_bits());

        let values = bcut2d(&topology, &[12.011; 2], &[0.0; 2], &[0.0; 2], &[0.0; 2]).unwrap();
        close(values[0], 12.012, 1.0e-12);
        close(values[1], 12.010, 1.0e-12);
    }

    #[test]
    fn literal_bond_weights_produce_analytical_pair_spectra() {
        for (order, expected_weight) in [
            (BondOrder::Single, 1.0_f64),
            (BondOrder::Double, 0.7071067811865475_f64),
            (BondOrder::Triple, 0.5773502691896258_f64),
            (BondOrder::Aromatic, 0.8164965809277261_f64),
        ] {
            assert_eq!(
                bond_weight(order).unwrap().to_bits(),
                expected_weight.to_bits()
            );

            let topology = graph(2, &[(0, 1, order)]);
            let values = bcut2d(&topology, &[0.0; 2], &[0.0; 2], &[0.0; 2], &[0.0; 2]).unwrap();
            close(values[0], expected_weight, 1e-14);
            close(values[1], -expected_weight, 1e-14);
        }
    }

    #[test]
    fn aromatic_disconnected_matrix_is_symmetric_and_fully_coupled() {
        let topology = parsed("c1ccccc1.C");
        let burden = BurdenMatrix::new(&topology).unwrap();
        let dimension = burden.dimension;

        for row in 0..dimension {
            for column in 0..dimension {
                let actual = burden.values[row * dimension + column];
                let transpose = burden.values[column * dimension + row];
                assert_eq!(actual.to_bits(), transpose.to_bits());
            }
        }

        let bonded = topology
            .bonds
            .iter()
            .map(|bond| (bond.begin().index(), bond.end().index()))
            .collect::<Vec<_>>();
        for row in 0..dimension {
            for column in row + 1..dimension {
                let actual = burden.values[row * dimension + column];
                let is_bonded = bonded.iter().any(|&(left, right)| {
                    (left == row && right == column) || (left == column && right == row)
                });
                let expected = if is_bonded {
                    0.8164965809277261
                } else {
                    NONBONDED_WEIGHT
                };
                assert_eq!(
                    actual.to_bits(),
                    expected.to_bits(),
                    "pair ({row}, {column})"
                );
            }
        }
    }

    #[test]
    fn medium_aromatic_topology_runs_all_four_diagonal_rows() {
        let topology = parsed("N[C@@H](Cc1ccccc1)C(=O)O");
        assert!(topology.atoms.len() >= 12);
        let atom_count = topology.atoms.len();
        let values = bcut2d(
            &topology,
            &(0..atom_count)
                .map(|index| 10.0 + index as f64)
                .collect::<Vec<_>>(),
            &properties(atom_count, -0.1),
            &properties(atom_count, 0.2),
            &properties(atom_count, 1.5),
        )
        .unwrap();

        for pair in values.chunks_exact(2) {
            assert!(pair[0].is_finite() && pair[1].is_finite());
            assert!(pair[0] >= pair[1]);
        }
    }

    #[test]
    fn pair_row_length_errors_identify_the_first_weighting() {
        let topology = graph(2, &[]);
        let error = bcut2d(&topology, &[12.011], &[0.0; 2], &[0.0; 2], &[0.0; 2]).unwrap_err();

        assert_eq!(
            error,
            DescriptorError::BcutDiagonalRows {
                field: "mass",
                actual: 1,
                expected: 2,
            }
        );
    }

    #[test]
    fn unsupported_bond_order_is_a_structured_error() {
        let topology = graph(2, &[(0, 1, cosmolkit_model::BondOrder::Dative)]);
        let error = bcut2d(&topology, &[12.011; 2], &[0.0; 2], &[0.0; 2], &[0.0; 2]).unwrap_err();

        assert_eq!(
            error,
            DescriptorError::UnsupportedBcutBond {
                bond: BondId::new(0),
                order: cosmolkit_model::BondOrder::Dative,
            }
        );
    }

    #[test]
    fn private_symmetric_solver_matches_known_three_by_three_spectrum() {
        let matrix = [2.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 2.0];
        let mut workspace = SymmetricEigenWorkspace::new(3).unwrap();
        let (highest, lowest) = symmetric_eigen_extrema(&matrix, &mut workspace).unwrap();

        close(highest, 2.0 + 2.0_f64.sqrt(), 1.0e-14);
        close(lowest, 2.0 - 2.0_f64.sqrt(), 1.0e-14);
    }
}
