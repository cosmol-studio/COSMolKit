use cosmolkit_core::{
    AddHsParams, HydrogenError, ValenceAssignment, ValenceError, add_hydrogens_with_source_valence,
};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, CoordinateBlock, Element,
    MoleculeProperties, TopologyBlock,
};

fn chain(no_implicit: bool) -> TopologyBlock {
    TopologyBlock::try_from_parts(
        [Element::C, Element::C, Element::O]
            .into_iter()
            .enumerate()
            .map(|(i, e)| {
                Atom::from_spec(
                    AtomId::new(i),
                    AtomSpec::new(e).with_no_implicit(no_implicit),
                )
            })
            .collect(),
        (0..2)
            .map(|i| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(i), AtomId::new(i + 1), BondOrder::Single),
                )
            })
            .collect(),
        vec![],
        vec![],
    )
    .unwrap()
}
fn stale() -> ValenceAssignment {
    ValenceAssignment {
        explicit_valence: vec![4, 4, 2],
        implicit_hydrogens: vec![0, 0, 0],
    }
}
fn run(
    source: TopologyBlock,
    params: AddHsParams,
    valence: Option<ValenceAssignment>,
) -> Result<cosmolkit_core::AddHydrogensResult, HydrogenError> {
    add_hydrogens_with_source_valence(
        source,
        CoordinateBlock::default(),
        MoleculeProperties::default(),
        &params,
        valence,
    )
}

#[test]
fn stale_source_counts_control_append_then_only_selected_scalars_refresh() {
    // AddHs.cpp snapshots cached implicit Hs before edits. These are the
    // RemoveHs(sanitize=false) ethanol scalars, not fresh CCO valence.
    let source = chain(false);
    let peer = source.clone();
    let input = stale();
    let original = input.clone();
    let out = run(source.clone(), AddHsParams::default(), Some(input.clone())).unwrap();
    assert_eq!(out.topology.atoms.len(), 3);
    assert_eq!(out.topology.bonds.len(), 2);
    assert_eq!(out.final_valence.explicit_valence, [1, 2, 1]);
    assert_eq!(out.final_valence.implicit_hydrogens, [3, 2, 1]);
    assert_eq!(source, peer);
    assert_eq!(input, original);
    let out = run(
        source,
        AddHsParams {
            only_on_atoms: Some(vec![AtomId::new(0)]),
            ..Default::default()
        },
        Some(input),
    )
    .unwrap();
    assert_eq!(out.topology.atoms.len(), 3);
    assert_eq!(out.final_valence.explicit_valence, [1, 4, 2]);
    assert_eq!(out.final_valence.implicit_hydrogens, [3, 0, 0]);
}

#[test]
fn selected_hydrogens_get_cache_rows_and_untouched_rows_keep_source_fields() {
    let out = run(
        chain(false),
        AddHsParams {
            only_on_atoms: Some(vec![AtomId::new(0)]),
            ..Default::default()
        },
        Some(ValenceAssignment {
            explicit_valence: vec![1, 66, 77],
            implicit_hydrogens: vec![3, 2, 1],
        }),
    )
    .unwrap();
    assert_eq!(out.topology.atoms.len(), 6);
    assert_eq!(out.topology.bonds.len(), 5);
    assert_eq!(out.final_valence.explicit_valence, [4, 66, 77, 1, 1, 1]);
    assert_eq!(out.final_valence.implicit_hydrogens, [0, 2, 1, 0, 0, 0]);
    assert!(
        out.topology.atoms[3..]
            .iter()
            .all(|a| a.element() == Element::H && a.implicit_hydrogen())
    );
    out.mapping.validate_for_counts(3, 6, 2, 5).unwrap();
}

#[test]
fn all_source_implicit_getters_run_before_explicit_only_and_empty_selection() {
    for only in [None, Some(vec![]), Some(vec![AtomId::new(0)])] {
        for explicit in [false, true] {
            let params = AddHsParams {
                explicit_only: explicit,
                only_on_atoms: only.clone(),
                ..Default::default()
            };
            assert_eq!(
                run(chain(false), params.clone(), None),
                Err(HydrogenError::Valence(
                    ValenceError::ImplicitValenceCacheNotInitialized {
                        atom: AtomId::new(0)
                    }
                ))
            );
            assert_eq!(
                run(
                    chain(false),
                    params,
                    Some(ValenceAssignment {
                        explicit_valence: vec![1, 2, 1],
                        implicit_hydrogens: vec![0, -1, 0]
                    })
                ),
                Err(HydrogenError::Valence(
                    ValenceError::ImplicitValenceCacheNotInitialized {
                        atom: AtomId::new(1)
                    }
                ))
            );
        }
    }
}

#[test]
fn no_implicit_short_circuit_keeps_uninitialized_unselected_source_scalars() {
    let out = run(
        chain(true),
        AddHsParams {
            only_on_atoms: Some(vec![AtomId::new(0)]),
            ..Default::default()
        },
        None,
    )
    .unwrap();
    assert_eq!(out.topology.atoms.len(), 3);
    assert_eq!(out.final_valence.explicit_valence, [1, -1, -1]);
    assert_eq!(out.final_valence.implicit_hydrogens, [0, -1, -1]);
    let out = run(
        chain(true),
        AddHsParams {
            only_on_atoms: Some(vec![]),
            ..Default::default()
        },
        None,
    )
    .unwrap();
    assert_eq!(out.final_valence.explicit_valence, [-1, -1, -1]);
    assert_eq!(out.final_valence.implicit_hydrogens, [-1, -1, -1]);
}

#[test]
fn source_scalar_getter_observes_signed_eight_bit_cache_before_checking_init() {
    let out = run(
        chain(false),
        AddHsParams {
            only_on_atoms: Some(vec![]),
            ..Default::default()
        },
        Some(ValenceAssignment {
            explicit_valence: vec![1, 2, 1],
            implicit_hydrogens: vec![256, 0, 0],
        }),
    )
    .unwrap();
    assert_eq!(out.topology.atoms.len(), 3);
    assert_eq!(out.final_valence.implicit_hydrogens, [256, 0, 0]);
    assert_eq!(
        run(
            chain(false),
            AddHsParams::default(),
            Some(ValenceAssignment {
                explicit_valence: vec![1, 2, 1],
                implicit_hydrogens: vec![255, 0, 0]
            })
        ),
        Err(HydrogenError::Valence(
            ValenceError::ImplicitValenceCacheNotInitialized {
                atom: AtomId::new(0)
            }
        ))
    );
}

#[test]
fn source_assignment_row_errors_keep_canonical_field_and_counts() {
    assert_eq!(
        run(
            chain(false),
            AddHsParams::default(),
            Some(ValenceAssignment {
                explicit_valence: vec![1],
                implicit_hydrogens: vec![3, 2, 1]
            })
        ),
        Err(HydrogenError::ValenceAssignmentLength {
            field: "explicit_valence",
            expected: 3,
            actual: 1
        })
    );
    assert_eq!(
        run(
            chain(false),
            AddHsParams::default(),
            Some(ValenceAssignment {
                explicit_valence: vec![1, 2, 1],
                implicit_hydrogens: vec![3]
            })
        ),
        Err(HydrogenError::ValenceAssignmentLength {
            field: "implicit_hydrogens",
            expected: 3,
            actual: 1
        })
    );
}

#[test]
fn selected_parent_refresh_uses_actual_signed_source_fields_before_implicit_calculation() {
    for (degree, expected_explicit, expected_implicit) in
        [(127, 127, 0), (128, -128, -124), (255, -1, 5), (256, 0, 4)]
    {
        let mut atoms = vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))];
        for i in 1..=degree {
            atoms.push(Atom::from_spec(
                AtomId::new(i),
                AtomSpec::new(Element::DUMMY).with_no_implicit(true),
            ));
        }
        let bonds = (0..degree)
            .map(|i| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(0), AtomId::new(i + 1), BondOrder::Single),
                )
            })
            .collect();
        let source = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
        let peer = source.clone();
        let input = ValenceAssignment {
            explicit_valence: vec![66; degree + 1],
            implicit_hydrogens: vec![0; degree + 1],
        };
        let out = run(
            source.clone(),
            AddHsParams {
                only_on_atoms: Some(vec![AtomId::new(0)]),
                ..Default::default()
            },
            Some(input.clone()),
        )
        .unwrap();
        // Source snapshots zero Hs before parent refresh: no H is appended.
        assert_eq!(out.topology.atoms.len(), degree + 1);
        assert_eq!(out.topology.bonds.len(), degree);
        assert_eq!(out.final_valence.explicit_valence[0], expected_explicit);
        assert_eq!(out.final_valence.implicit_hydrogens[0], expected_implicit);
        assert_eq!(
            &out.final_valence.explicit_valence[1..],
            &input.explicit_valence[1..]
        );
        assert_eq!(
            &out.final_valence.implicit_hydrogens[1..],
            &input.implicit_hydrogens[1..]
        );
        let assigned = cosmolkit_core::assign_valence_with_options_for_topology(
            &source,
            cosmolkit_core::ValenceModel::RdkitLike,
            false,
        )
        .unwrap();
        assert_eq!(assigned.explicit_valence[0], expected_explicit);
        assert_eq!(assigned.implicit_hydrogens[0], expected_implicit);
        assert_eq!(&assigned.explicit_valence[1..], vec![1; degree].as_slice());
        assert_eq!(
            &assigned.implicit_hydrogens[1..],
            vec![0; degree].as_slice()
        );
        assert_eq!(source, peer);
    }
}
