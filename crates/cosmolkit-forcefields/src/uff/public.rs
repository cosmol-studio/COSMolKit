//! Narrow detached prepared UFF entries; optimizer internals remain private.
use cosmolkit_core::{RingInfo, ValenceAssignment};
use cosmolkit_model::{CoordinateBlock, MoleculeProperties, TopologyBlock};
use std::{error::Error, fmt};

/// Options for optimization of one explicitly selected 3D conformer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UffSingleOptions {
    pub conformer_id: usize,
    pub max_iterations: i32,
    pub vdw_threshold: f64,
    pub ignore_interfragment_interactions: bool,
}

/// Scalar outcome; coordinates are written to the borrowed block.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UffSingleOutcome {
    pub status: i32,
    pub energy: f64,
}

/// Concrete source chain without exposing builder or optimizer machinery.
#[derive(Debug)]
pub struct UffSingleError(super::optimization::UffPreparedOptimizationError);

impl fmt::Display for UffSingleError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.0, formatter)
    }
}
impl Error for UffSingleError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(&self.0)
    }
}

/// Options for serial optimization of every stored 3D conformer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UffConformerOptions {
    pub max_iterations: i32,
    pub vdw_threshold: f64,
    pub ignore_interfragment_interactions: bool,
}

/// Scalar result for one stored 3D conformer; coordinates are written in place.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UffConformerOutcome {
    pub conformer_id: usize,
    pub status: i32,
    pub energy: f64,
}

/// Concrete source chain without exposing builder or optimizer machinery.
#[derive(Debug)]
pub struct UffConformerError(super::optimization::UffPreparedOptimizationError);

impl fmt::Display for UffConformerError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.0, formatter)
    }
}
impl Error for UffConformerError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(&self.0)
    }
}

/// Optimize one stored 3D conformer using trusted borrowed chemistry state.
pub fn optimize_uff_single_prepared(
    topology: &TopologyBlock,
    coordinates: &mut CoordinateBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    properties: &MoleculeProperties,
    options: UffSingleOptions,
) -> Result<UffSingleOutcome, UffSingleError> {
    // RDKit❗✔️:   std::unique_ptr<ForceFields::ForceField> ff(UFF::constructForceField(
    // RDKit❗✔️:       mol, vdwThresh, confId, ignoreInterfragInteractions));
    // RDKit❗✔️:   std::pair<int, double> res =
    // RDKit❗✔️:       ForceFieldsHelper::OptimizeMolecule(*ff, maxIters);
    // RDKit❗✔️:   return res;
    // Behavior: delegate once to the accepted UFF.h:43-47 owner, which holds
    // the full source bodies and typed stage/error semantics.
    // Complexity: no graph, valence, coordinate or contribution projection.
    // Diagnostics allocate only when the existing typer emits a diagnostic.
    let mut diagnostics = Vec::new();
    let outcome = super::optimization::optimize_prepared_uff_single(
        topology,
        coordinates,
        valence,
        rings,
        properties,
        &mut diagnostics,
        super::convenience::SingleConformerOptions {
            conformer_id: options.conformer_id,
            max_iterations: options.max_iterations,
            vdw_threshold: options.vdw_threshold,
            ignore_interfragment_interactions: options.ignore_interfragment_interactions,
        },
    )
    .map_err(UffSingleError)?;
    Ok(UffSingleOutcome {
        status: outcome.status,
        energy: outcome.energy,
    })
}

/// Optimize all stored 3D conformers serially using trusted borrowed chemistry state.
#[allow(clippy::too_many_arguments)]
pub fn optimize_uff_conformers_prepared(
    topology: &TopologyBlock,
    coordinates: &mut CoordinateBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    properties: &MoleculeProperties,
    options: UffConformerOptions,
) -> Result<Vec<UffConformerOutcome>, UffConformerError> {
    // RDKit❗❌:   std::unique_ptr<ForceFields::ForceField> ff(UFF::constructForceField(
    // RDKit❗❌:       mol, vdwThresh, -1, ignoreInterfragInteractions));
    // RDKit❗❌:   ForceFieldsHelper::OptimizeMoleculeConfs(mol, *ff, res, numThreads, maxIters);
    // RDKit❗❌: }
    // Behavior: use the existing serial source owner once, which preserves
    // construct-then-resize ordering and visits stored 3D conformers in order.
    // The returned scalar rows are paired with those same stored IDs.
    // Complexity: this facade adds only the required O(C) public result rows;
    // it creates no graph, coordinate, builder, or intermediate-row copy.
    let construction_conformer_id = coordinates
        .conformers_3d
        .first()
        .map_or(0, |conformer| conformer.id());
    let mut diagnostics = Vec::new();
    let mut results = Vec::new();
    super::optimization::optimize_prepared_uff_serial(
        topology,
        coordinates,
        &mut results,
        valence,
        rings,
        properties,
        &mut diagnostics,
        super::convenience::SingleConformerOptions {
            conformer_id: construction_conformer_id,
            max_iterations: options.max_iterations,
            vdw_threshold: options.vdw_threshold,
            ignore_interfragment_interactions: options.ignore_interfragment_interactions,
        },
    )
    .map_err(UffConformerError)?;

    debug_assert_eq!(results.len(), coordinates.conformers_3d.len());
    Ok(results
        .into_iter()
        .enumerate()
        .map(|(index, outcome)| UffConformerOutcome {
            conformer_id: coordinates.conformers_3d[index].id(),
            status: outcome.status,
            energy: outcome.energy,
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use std::error::Error as _;

    use cosmolkit_core::{RingSearchParams, ValenceAssignment, fast_find_rings, symmetrized_sssr};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Conformer2D, Conformer3D,
        CoordinateBlock, Element, Hybridization, MoleculeProperties, TopologyBlock,
    };

    use super::{UffConformerOptions, optimize_uff_conformers_prepared};

    const IDS: [usize; 3] = [7, 3, 11];
    const BASE_POINTS: [[f64; 3]; 4] = [
        [0.0, 0.0, 0.0],
        [1.9, 0.2, 0.0],
        [5.0, 1.0, 0.0],
        [6.7, 1.1, 0.3],
    ];
    const ONE_ITERATION_COORDINATE_BITS: [[[u64; 3]; 4]; 3] = [
        [
            [4598641110794173855, 4584026890252664668, 0],
            [4609993285745583870, 4595327571088798126, 0],
            [
                4617482136253981157,
                4607221623104106851,
                4583184447793597021,
            ],
            [
                4619062929510853352,
                4607593574458665015,
                4598605487821677197,
            ],
        ],
        [
            [4601035628630809199, 4586362091530870878, 0],
            [4610250340215625428, 4595544038010938386, 0],
            [
                4618178990165372210,
                4607703489880483840,
                4586949088337815106,
            ],
            [
                4619683378490468170,
                4608057463604035831,
                4598767838013282391,
            ],
        ],
        [
            [4603054482820045727, 4588378527603826956, 0],
            [4610507394685666986, 4595760504933078645, 0],
            [
                4618875844076763260,
                4608185356656860828,
                4589570881368892312,
            ],
            [
                4620303827470082988,
                4608521352749406646,
                4598930188204887588,
            ],
        ],
    ];

    fn cc_cc_topology() -> TopologyBlock {
        let atoms = (0..4)
            .map(|row| {
                Atom::from_spec(
                    AtomId::new(row),
                    AtomSpec::new(Element::from_atomic_number(6).expect("carbon is present"))
                        .with_hybridization(Hybridization::Sp3),
                )
            })
            .collect();
        let bonds = [(0, 1), (2, 3)]
            .into_iter()
            .enumerate()
            .map(|(row, (begin, end))| {
                Bond::from_spec(
                    BondId::new(row),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed CC.CC topology is valid")
    }

    fn initial_conformer_rows(count: usize) -> Vec<Vec<[f64; 3]>> {
        [1.0, 1.1, 1.2][..count]
            .iter()
            .map(|scale| {
                BASE_POINTS
                    .iter()
                    .map(|point| [point[0] * *scale, point[1] * *scale, point[2] * *scale])
                    .collect()
            })
            .collect()
    }

    fn test_valence() -> ValenceAssignment {
        ValenceAssignment {
            explicit_valence: vec![1; 4],
            implicit_hydrogens: vec![3; 4],
        }
    }

    fn options(max_iterations: i32) -> UffConformerOptions {
        UffConformerOptions {
            max_iterations,
            vdw_threshold: 10.0,
            ignore_interfragment_interactions: true,
        }
    }

    fn actual_coordinate_bits(coordinates: &CoordinateBlock) -> Vec<[[u64; 3]; 4]> {
        coordinates
            .conformers_3d
            .iter()
            .map(|conformer| {
                std::array::from_fn(|atom| {
                    std::array::from_fn(|axis| conformer.coordinates()[atom][axis].to_bits())
                })
            })
            .collect()
    }

    #[test]
    fn uff_all_facade_reference_rows_keep_stored_ids_and_2d() {
        const TWO_D: [[f64; 2]; 4] = [[17.0, -3.5], [-2.25, 18.125], [9.5, -0.25], [-8.0, 0.125]];
        let topology = cc_cc_topology();
        let valence = test_valence();
        let rings = fast_find_rings(&topology).expect("fixed CC.CC ring state");
        let properties = MoleculeProperties::default();

        for count in [1, 3] {
            for max_iterations in [0, 1] {
                let initial_rows = initial_conformer_rows(count);
                let mut coordinates = CoordinateBlock::default();
                for (row, points) in initial_rows.iter().enumerate() {
                    coordinates.conformers_3d.push(
                        Conformer3D::new(IDS[row], points.clone(), true)
                            .with_prop("source-row", format!("row-{row}")),
                    );
                }
                coordinates.conformers_2d.push(
                    Conformer2D::new(IDS[0], TWO_D.to_vec())
                        .with_prop("layout", "unchanged-two-dimensional-row"),
                );
                let initial_2d = coordinates.conformers_2d.clone();

                let outcomes = optimize_uff_conformers_prepared(
                    &topology,
                    &mut coordinates,
                    &valence,
                    &rings,
                    &properties,
                    options(max_iterations),
                )
                .expect("source-supported CC.CC rows optimize");

                assert_eq!(
                    outcomes
                        .iter()
                        .map(|outcome| outcome.conformer_id)
                        .collect::<Vec<_>>(),
                    IDS[..count]
                );
                assert!(outcomes.iter().all(|outcome| outcome.status == 1));
                let expected_energy_bits: &[u64] = match (count, max_iterations) {
                    (1, 0) => &[4634709622640018448],
                    (1, 1) => &[4622575255082429943],
                    (3, 0) => &[
                        4634709622640018448,
                        4640306196202131730,
                        4644374702382500982,
                    ],
                    (3, 1) => &[
                        4622575255082429943,
                        4628472156411023569,
                        4632421107459294556,
                    ],
                    _ => unreachable!(),
                };
                assert_eq!(
                    outcomes
                        .iter()
                        .map(|outcome| outcome.energy.to_bits())
                        .collect::<Vec<_>>(),
                    expected_energy_bits
                );

                let expected_coordinate_bits = if max_iterations == 0 {
                    initial_rows
                        .iter()
                        .map(|row| {
                            std::array::from_fn(|atom| {
                                std::array::from_fn(|axis| row[atom][axis].to_bits())
                            })
                        })
                        .collect::<Vec<_>>()
                } else {
                    ONE_ITERATION_COORDINATE_BITS[..count].to_vec()
                };
                assert_eq!(
                    actual_coordinate_bits(&coordinates),
                    expected_coordinate_bits
                );
                assert_eq!(coordinates.conformers_2d, initial_2d);
            }
        }
    }

    #[test]
    fn uff_all_cost_counts_construction_and_lazy_source_ordered_rows() {
        const TWO_D: [[f64; 2]; 4] = [[17.0, -3.5], [-2.25, 18.125], [9.5, -0.25], [-8.0, 0.125]];
        let topology = cc_cc_topology();
        let valence = test_valence();
        let warm_rings = fast_find_rings(&topology).expect("fixed CC.CC cached ring state");
        let properties = MoleculeProperties::default();

        for count in [0, 1, 3] {
            for max_iterations in [0, 1] {
                for use_warm_rings in [false, true] {
                    let cold_rings = (!use_warm_rings).then(|| {
                        symmetrized_sssr(&topology, &RingSearchParams::default())
                            .expect("fixed CC.CC temporary ring state")
                    });
                    let rings = cold_rings.as_ref().unwrap_or(&warm_rings);
                    let initial_rows = initial_conformer_rows(count);
                    let mut coordinates = CoordinateBlock::default();
                    for (row, points) in initial_rows.iter().enumerate() {
                        coordinates.conformers_3d.push(
                            Conformer3D::new(IDS[row], points.clone(), true)
                                .with_prop("source-row", format!("row-{row}")),
                        );
                    }
                    coordinates.conformers_2d.push(
                        Conformer2D::new(IDS[0], TWO_D.to_vec())
                            .with_prop("layout", "unchanged-two-dimensional-row"),
                    );
                    let initial_2d = coordinates.conformers_2d.clone();

                    super::super::convenience::uff_all_cost_probe_start();
                    let result = optimize_uff_conformers_prepared(
                        &topology,
                        &mut coordinates,
                        &valence,
                        rings,
                        &properties,
                        options(max_iterations),
                    );
                    let (construction_calls, rows_at_stage_entry, visited_ids) =
                        super::super::convenience::uff_all_cost_probe_finish()
                            .expect("the actual serial path records its cost probe");

                    assert_eq!(construction_calls, 1);
                    if count == 0 {
                        assert_eq!(rows_at_stage_entry, None);
                        assert!(visited_ids.is_empty());
                        let error = result.expect_err("source construction requires a 3D row");
                        assert!(matches!(
                            &error.0,
                            super::super::optimization::UffPreparedOptimizationError::Serial(
                                super::super::convenience::SerialUffOptimizationError::Construction(
                                    super::super::builder::AutomaticForceFieldConstructionError::Construction(
                                        super::super::builder::ForceFieldConstructionError::Builder(
                                            super::super::builder::UffBuilderError::SelectedThreeDimensionalConformerNotFound {
                                                conformer_id: 0
                                            }
                                        )
                                    )
                                )
                            )
                        ));
                        assert!(coordinates.conformers_3d.is_empty());
                        assert_eq!(coordinates.conformers_2d, initial_2d);
                        continue;
                    }

                    assert_eq!(rows_at_stage_entry, Some(0));
                    assert_eq!(visited_ids, IDS[..count]);
                    let outcomes = result.expect("source-supported CC.CC rows optimize");
                    assert_eq!(
                        outcomes
                            .iter()
                            .map(|outcome| outcome.conformer_id)
                            .collect::<Vec<_>>(),
                        IDS[..count]
                    );
                    assert!(outcomes.iter().all(|outcome| outcome.status == 1));
                    let expected_energy_bits: &[u64] = match (count, max_iterations) {
                        (1, 0) => &[4634709622640018448],
                        (1, 1) => &[4622575255082429943],
                        (3, 0) => &[
                            4634709622640018448,
                            4640306196202131730,
                            4644374702382500982,
                        ],
                        (3, 1) => &[
                            4622575255082429943,
                            4628472156411023569,
                            4632421107459294556,
                        ],
                        _ => unreachable!(),
                    };
                    assert_eq!(
                        outcomes
                            .iter()
                            .map(|outcome| outcome.energy.to_bits())
                            .collect::<Vec<_>>(),
                        expected_energy_bits
                    );
                    let expected_coordinate_bits = if max_iterations == 0 {
                        initial_rows
                            .iter()
                            .map(|row| {
                                std::array::from_fn(|atom| {
                                    std::array::from_fn(|axis| row[atom][axis].to_bits())
                                })
                            })
                            .collect::<Vec<_>>()
                    } else {
                        ONE_ITERATION_COORDINATE_BITS[..count].to_vec()
                    };
                    assert_eq!(
                        actual_coordinate_bits(&coordinates),
                        expected_coordinate_bits
                    );
                    assert_eq!(coordinates.conformers_2d, initial_2d);
                }
            }
        }
    }

    #[test]
    fn uff_all_facade_missing_3d_preserves_concrete_error_sources() {
        use super::super::{
            builder::{
                AutomaticForceFieldConstructionError, ForceFieldConstructionError, UffBuilderError,
            },
            convenience::SerialUffOptimizationError,
            optimization::UffPreparedOptimizationError,
        };

        fn assert_stored_child<Parent, Child>(parent: &Parent, stored: &Child)
        where
            Parent: std::error::Error,
            Child: std::error::Error + 'static,
        {
            let reported = std::error::Error::source(parent)
                .expect("the error exposes its concrete stored child");
            let stored_error: &(dyn std::error::Error + 'static) = stored;
            assert!(
                std::ptr::eq(reported, stored_error),
                "source is not the stored child: {} -> {}",
                std::any::type_name::<Parent>(),
                std::any::type_name::<Child>()
            );
        }

        let topology = cc_cc_topology();
        let valence = test_valence();
        let rings = fast_find_rings(&topology).expect("fixed CC.CC ring state");
        let properties = MoleculeProperties::default();
        let mut coordinates = CoordinateBlock::default();
        coordinates
            .conformers_2d
            .push(Conformer2D::new(7, vec![[17.0, -3.5]; 4]).with_prop("layout", "2d-only"));

        let error = optimize_uff_conformers_prepared(
            &topology,
            &mut coordinates,
            &valence,
            &rings,
            &properties,
            options(0),
        )
        .expect_err("the source default requires a stored 3D conformer");

        let prepared = &error.0;
        let serial = match prepared {
            UffPreparedOptimizationError::Serial(stored) => {
                assert_stored_child(&error, prepared);
                assert_stored_child(prepared, stored);
                stored
            }
            other => panic!("unexpected prepared error branch: {other:?}"),
        };
        let construction = match serial {
            SerialUffOptimizationError::Construction(stored) => {
                assert_stored_child(serial, stored);
                stored
            }
            other => panic!("unexpected serial error branch: {other:?}"),
        };
        let field_construction = match construction {
            AutomaticForceFieldConstructionError::Construction(stored) => {
                assert_stored_child(construction, stored);
                stored
            }
            other => panic!("unexpected automatic-construction branch: {other:?}"),
        };
        let builder_error = match field_construction {
            ForceFieldConstructionError::Builder(stored) => {
                assert_stored_child(field_construction, stored);
                stored
            }
            other => panic!("unexpected construction error branch: {other:?}"),
        };
        assert!(matches!(
            builder_error,
            UffBuilderError::SelectedThreeDimensionalConformerNotFound { .. }
        ));
        assert_eq!(error.to_string(), error.0.to_string());
        assert!(coordinates.conformers_3d.is_empty());
        assert_eq!(coordinates.conformers_2d[0].id(), 7);
    }
}
