//! Exact original descriptor literals at canonical public query boundaries.
#![cfg(all(
    feature = "cap-descriptors",
    feature = "cap-smiles",
    feature = "cap-sanitize"
))]
use cosmolkit::{
    AtomSpec, BondOrder, BondSpec, DescriptorReadError, Element, LabuteAsaContributions, Molecule,
    MoleculeBuilder,
};
fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|x| x.to_bits()).collect()
}
fn hydrogen_first(value: LabuteAsaContributions) -> Vec<f64> {
    std::iter::once(value.hydrogen_contribution)
        .chain(value.atom_contributions)
        .collect()
}

fn explicit_methane() -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let carbon = builder.add_atom(AtomSpec::new(Element::C));
    for _ in 0..4 {
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
            .expect("explicit methane bond");
    }
    builder
        .build()
        .expect("explicit methane graph")
        .sanitize()
        .expect("sanitized explicit methane")
}

fn explicit_methylamine() -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let nitrogen = builder.add_atom(AtomSpec::new(Element::N));
    let carbon = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(BondSpec::new(nitrogen, carbon, BondOrder::Single))
        .expect("N-C bond");
    for _ in 0..2 {
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(BondSpec::new(nitrogen, hydrogen, BondOrder::Single))
            .expect("N-H bond");
    }
    builder
        .build()
        .expect("methylamine graph")
        .sanitize()
        .expect("sanitized methylamine")
}

const SLOGP_SCALARS: [fn(&Molecule) -> Result<f64, DescriptorReadError>; 12] = [
    Molecule::slogp_vsa_1,
    Molecule::slogp_vsa_2,
    Molecule::slogp_vsa_3,
    Molecule::slogp_vsa_4,
    Molecule::slogp_vsa_5,
    Molecule::slogp_vsa_6,
    Molecule::slogp_vsa_7,
    Molecule::slogp_vsa_8,
    Molecule::slogp_vsa_9,
    Molecule::slogp_vsa_10,
    Molecule::slogp_vsa_11,
    Molecule::slogp_vsa_12,
];

const SMR_SCALARS: [fn(&Molecule) -> Result<f64, DescriptorReadError>; 10] = [
    Molecule::smr_vsa_1,
    Molecule::smr_vsa_2,
    Molecule::smr_vsa_3,
    Molecule::smr_vsa_4,
    Molecule::smr_vsa_5,
    Molecule::smr_vsa_6,
    Molecule::smr_vsa_7,
    Molecule::smr_vsa_8,
    Molecule::smr_vsa_9,
    Molecule::smr_vsa_10,
];

#[test]
fn direct_and_standard_lipinski_counts_match_pinned_rdkit_distinct_definitions() {
    const CASES: [(&str, [u32; 8]); 7] = [
        ("", [0, 0, 0, 0, 0, 0, 0, 0]),
        ("CCO", [1, 1, 1, 1, 1, 0, 3, 9]),
        ("NC(=O)C", [2, 2, 1, 1, 2, 1, 4, 9]),
        ("NC(=O)N", [3, 4, 1, 2, 3, 2, 4, 8]),
        ("NCC(=O)O", [3, 3, 2, 2, 3, 0, 5, 10]),
        ("c1ncc[nH]1", [2, 1, 1, 1, 2, 0, 5, 9]),
        ("[Na+].[Cl-]", [0, 0, 0, 0, 2, 0, 2, 2]),
    ];

    for (smiles, expected) in CASES {
        let molecule = Molecule::from_smiles(smiles).expect("Lipinski fixture");
        let actual = [
            molecule.lipinski_hba().unwrap(),
            molecule.lipinski_hbd().unwrap(),
            molecule.num_hba().unwrap(),
            molecule.num_hbd().unwrap(),
            molecule.num_heteroatoms().unwrap(),
            molecule.num_amide_bonds().unwrap(),
            molecule.num_heavy_atoms().unwrap(),
            molecule.total_atom_count().unwrap(),
        ];
        assert_eq!(actual, expected, "{smiles:?} Lipinski count vector");
    }
}

#[test]
fn lipinski_hydrogen_counts_include_explicit_neighbors_without_double_counting_atoms() {
    let molecule = explicit_methylamine();
    assert_eq!(molecule.lipinski_hba().unwrap(), 1);
    assert_eq!(molecule.lipinski_hbd().unwrap(), 2);
    assert_eq!(molecule.num_hba().unwrap(), 1);
    assert_eq!(molecule.num_hbd().unwrap(), 1);
    assert_eq!(molecule.num_heteroatoms().unwrap(), 1);
    assert_eq!(molecule.num_amide_bonds().unwrap(), 0);
    assert_eq!(molecule.num_heavy_atoms().unwrap(), 2);
    assert_eq!(molecule.total_atom_count().unwrap(), 7);
}

#[test]
fn ring_classification_matches_pinned_rdkit_across_every_topology_family() {
    const CASES: [(&str, [u32; 11]); 11] = [
        ("CC", [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        ("c1ccccc1", [1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0]),
        ("c1ncccc1", [1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0]),
        ("C1CCCCC1", [1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1]),
        ("C1=CCCCC1", [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0]),
        ("O1CCCCC1", [1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0]),
        ("C1CCC2CCCCC2C1", [2, 0, 0, 2, 2, 0, 0, 0, 2, 0, 2]),
        ("C1CC2CCC1C2", [2, 0, 0, 2, 2, 0, 0, 0, 2, 0, 2]),
        ("C1CCC2(CC1)CCCCC2", [2, 0, 0, 2, 2, 0, 0, 0, 2, 0, 2]),
        ("C1CCCCCCCCCCC1", [1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1]),
        ("c1ccc2ncccc2c1", [2, 1, 2, 0, 0, 1, 1, 0, 0, 0, 0]),
    ];

    for (smiles, expected) in CASES {
        let molecule = Molecule::from_smiles(smiles).expect("ring fixture");
        let actual = [
            molecule.num_rings().unwrap(),
            molecule.num_heterocycles().unwrap(),
            molecule.num_aromatic_rings().unwrap(),
            molecule.num_saturated_rings().unwrap(),
            molecule.num_aliphatic_rings().unwrap(),
            molecule.num_aromatic_heterocycles().unwrap(),
            molecule.num_aromatic_carbocycles().unwrap(),
            molecule.num_aliphatic_heterocycles().unwrap(),
            molecule.num_aliphatic_carbocycles().unwrap(),
            molecule.num_saturated_heterocycles().unwrap(),
            molecule.num_saturated_carbocycles().unwrap(),
        ];
        assert_eq!(actual, expected, "{smiles:?} ring descriptor vector");
    }
}

#[test]
fn stereo_center_counts_match_pinned_rdkit_without_mutating_the_caller() {
    const CASES: [(&str, u32, u32); 7] = [
        ("", 0, 0),
        ("CC(F)Cl", 1, 1),
        ("C[C@H](F)Cl", 1, 0),
        ("CC(F)C", 0, 0),
        ("CC(F)Cl.C[C@@H](Br)I", 2, 1),
        ("FC(Cl)(Br)I", 1, 1),
        ("F[C@](Cl)(Br)I", 1, 0),
    ];

    for (smiles, expected_all, expected_unspecified) in CASES {
        let molecule = Molecule::from_smiles(smiles).expect("stereo-center fixture");
        let before = molecule.clone();
        assert_eq!(
            molecule.num_atom_stereo_centers().unwrap(),
            expected_all,
            "{smiles:?} all stereo centers"
        );
        assert_eq!(
            molecule.num_unspecified_atom_stereo_centers().unwrap(),
            expected_unspecified,
            "{smiles:?} unspecified stereo centers"
        );
        assert_eq!(
            molecule.num_atom_stereo_centers().unwrap(),
            expected_all,
            "{smiles:?} repeated all stereo centers"
        );
        assert_eq!(molecule, before, "{smiles:?} descriptor caller state");
    }

    let mut preassigned =
        Molecule::from_smiles("C[C@H](F)Cl").expect("preassigned stereo-center fixture");
    let mut properties = preassigned.properties().clone();
    properties
        .set_computed_prop("_StereochemDone", "1")
        .unwrap();
    preassigned = preassigned
        .to_builder()
        .with_properties(properties)
        .build()
        .unwrap()
        .with_assigned_valence()
        .unwrap()
        .with_assigned_rings()
        .unwrap();
    assert_eq!(preassigned.properties().prop("_StereochemDone"), Some("1"));
    let before = preassigned.clone();
    assert_eq!(preassigned.num_atom_stereo_centers().unwrap(), 1);
    assert_eq!(
        preassigned.num_unspecified_atom_stereo_centers().unwrap(),
        0
    );
    assert_eq!(preassigned, before);
}

#[test]
fn chi_descriptors_match_pinned_rdkit_across_path_and_element_branches() {
    #[derive(Clone, Copy)]
    enum Fixture {
        Smiles(&'static str),
        ExplicitMethane,
    }

    const CASES: [(&str, Fixture, [u64; 12]); 5] = [
        (
            "open path",
            Fixture::Smiles("CCCCCC"),
            [
                0x4013_504f_333f_9de6,
                0x4007_504f_333f_9de6,
                0x4013_504f_333f_9de6,
                0x4007_504f_333f_9de6,
                0x3ffb_504f_333f_9de4,
                0x3fee_a09e_667f_3bca,
                0x3fdf_ffff_ffff_fffd,
                0x4013_504f_333f_9de6,
                0x4007_504f_333f_9de6,
                0x3ffb_504f_333f_9de4,
                0x3fee_a09e_667f_3bca,
                0x3fdf_ffff_ffff_fffd,
            ],
        ),
        (
            "closed three-membered ring",
            Fixture::Smiles("C1CC1"),
            [
                0x4000_f876_ccdf_6cda,
                0x3ff8_0000_0000_0000,
                0x4000_f876_ccdf_6cd9,
                0x3ff7_ffff_ffff_fffe,
                0x3ff0_f876_ccdf_6cd8,
                0x3fd6_a09e_667f_3bcb,
                0x0000_0000_0000_0000,
                0x4000_f876_ccdf_6cd9,
                0x3ff7_ffff_ffff_fffe,
                0x3ff0_f876_ccdf_6cd8,
                0x3fd6_a09e_667f_3bcb,
                0x0000_0000_0000_0000,
            ],
        ),
        (
            "disconnected graph",
            Fixture::Smiles("CC.CCC"),
            [
                0x4012_d413_cccf_e77a,
                0x4003_504f_333f_9de6,
                0x4012_d413_cccf_e77a,
                0x4003_504f_333f_9de6,
                0x3fe6_a09e_667f_3bcc,
                0x0000_0000_0000_0000,
                0x0000_0000_0000_0000,
                0x4012_d413_cccf_e77a,
                0x4003_504f_333f_9de6,
                0x3fe6_a09e_667f_3bcc,
                0x0000_0000_0000_0000,
                0x0000_0000_0000_0000,
            ],
        ),
        (
            "explicit hydrogens",
            Fixture::ExplicitMethane,
            [
                0x4012_0000_0000_0000,
                0x4000_0000_0000_0000,
                0x3fe0_0000_0000_0000,
                0x0000_0000_0000_0000,
                0x0000_0000_0000_0000,
                0x0000_0000_0000_0000,
                0x0000_0000_0000_0000,
                0x4012_0000_0000_0000,
                0x4000_0000_0000_0000,
                0x0000_0000_0000_0000,
                0x0000_0000_0000_0000,
                0x0000_0000_0000_0000,
            ],
        ),
        (
            "heavy elements",
            Fixture::Smiles("[SiH3][GeH3]"),
            [
                0x4000_0000_0000_0000,
                0x3ff0_0000_0000_0000,
                0x4020_646e_1721_1cc0,
                0x402f_2d4a_4563_5640,
                0x0000_0000_0000_0000,
                0x0000_0000_0000_0000,
                0x0000_0000_0000_0000,
                0x4000_0000_0000_0000,
                0x3ff0_0000_0000_0000,
                0x0000_0000_0000_0000,
                0x0000_0000_0000_0000,
                0x0000_0000_0000_0000,
            ],
        ),
    ];

    for (case, fixture, expected) in CASES {
        let molecule = match fixture {
            Fixture::Smiles(smiles) => Molecule::from_smiles(smiles).unwrap_or_else(|error| {
                panic!("failed to parse {case} fixture {smiles:?}: {error}")
            }),
            Fixture::ExplicitMethane => explicit_methane(),
        };
        let actual = [
            molecule.chi_0().unwrap(),
            molecule.chi_1().unwrap(),
            molecule.chi_0_v_with_params(true).unwrap(),
            molecule.chi_1_v_with_params(true).unwrap(),
            molecule.chi_2_v_with_params(true).unwrap(),
            molecule.chi_3_v_with_params(true).unwrap(),
            molecule.chi_4_v_with_params(true).unwrap(),
            molecule.chi_0_n_with_params(true).unwrap(),
            molecule.chi_1_n_with_params(true).unwrap(),
            molecule.chi_2_n_with_params(true).unwrap(),
            molecule.chi_3_n_with_params(true).unwrap(),
            molecule.chi_4_n_with_params(true).unwrap(),
        ];
        for (index, (actual, expected)) in actual.into_iter().zip(expected).enumerate() {
            assert_eq!(
                actual.to_bits(),
                expected,
                "{case} fixed Chi field {index}: actual={actual:?}"
            );
        }
        for (order, expected_index) in [(2, 4), (3, 5), (4, 6)] {
            let actual = molecule.chi_n_v_with_params(order, true).unwrap();
            assert_eq!(
                actual.to_bits(),
                expected[expected_index],
                "{case} generic valence Chi order {order}: actual={actual:?}"
            );
        }
        for (order, expected_index) in [(2, 9), (3, 10), (4, 11)] {
            let actual = molecule.chi_n_n_with_params(order, true).unwrap();
            assert_eq!(
                actual.to_bits(),
                expected[expected_index],
                "{case} generic nVal Chi order {order}: actual={actual:?}"
            );
        }
    }
}

#[test]
fn kappa_and_phi_match_pinned_rdkit_across_boundary_and_ring_branches() {
    const CASES: [(&str, [u64; 5]); 6] = [
        ("", [0, 0, 0, 0, 0]),
        (
            "CCC",
            [
                0,
                0x4008_0000_0000_0000,
                0x4000_0000_0000_0000,
                0,
                0x4000_0000_0000_0000,
            ],
        ),
        (
            "CCCC",
            [
                0,
                0x4010_0000_0000_0000,
                0x4008_0000_0000_0000,
                0x4000_0000_0000_0000,
                0x4008_0000_0000_0000,
            ],
        ),
        (
            "C1CC1",
            [
                0,
                0x3ff5_5555_5555_5555,
                0x3fcc_71c7_1c71_c71c,
                0,
                0x3fb9_48b0_fcd6_e9e0,
            ],
        ),
        (
            "c1ccccc1",
            [
                0xbfe8_f5c2_8f5c_28f6,
                0x400b_4ae5_ac96_d07c,
                0x3ff9_b13b_4bc5_063b,
                0x3fe2_a303_c5f0_83cd,
                0x3fed_3790_5fe9_fff9,
            ],
        ),
        (
            "[Na+].[Cl-]",
            [
                0x3ff4_a3d7_0a3d_70a4,
                0x4024_bc52_e1e5_3581,
                0x4002_51eb_851e_b852,
                0x3fb0_b08a_703e_3bdb,
                0x4027_be07_dc40_0b58,
            ],
        ),
    ];

    for (smiles, expected) in CASES {
        let molecule = Molecule::from_smiles(smiles).expect("Kappa fixture");
        let actual = [
            molecule.hall_kier_alpha().unwrap(),
            molecule.kappa_1().unwrap(),
            molecule.kappa_2().unwrap(),
            molecule.kappa_3().unwrap(),
            molecule.phi().unwrap(),
        ];
        for (index, (actual, expected)) in actual.into_iter().zip(expected).enumerate() {
            assert_eq!(
                actual.to_bits(),
                expected,
                "{smiles:?} Hall-Kier/Kappa/Phi field {index}: actual={actual:?}"
            );
        }
    }
}

#[test]
fn mqn_matches_pinned_rdkit_complete_vectors_for_focused_branches() {
    const CASES: [(&str, [u32; 42]); 5] = [
        ("", [0; 42]),
        (
            "CCO",
            [
                2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 3, 2, 0, 0, 0, 0, 0, 0, 2, 1, 1, 1, 0, 0, 2, 1, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            ],
        ),
        (
            "c1ccccc1",
            [
                6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 6, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
            ],
        ),
        (
            "[Na+].[Cl-]",
            [
                0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            ],
        ),
        (
            "C1CC2CCC1C2",
            [
                7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 5, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 3, 2,
            ],
        ),
    ];

    for (smiles, expected) in CASES {
        let molecule = Molecule::from_smiles(smiles).expect("focused MQN fixture");
        assert_eq!(
            molecule.mqns(true).unwrap(),
            expected.to_vec(),
            "{smiles:?} MQN vector"
        );
    }
}

#[test]
fn mqn_complete_branch_fixture_exercises_and_matches_every_index() {
    const SMILES: &str = concat!(
        "CC(F)(Cl)Br.CPI.N1CCOCC1.[O-]C(=O)[NH3+].c1ncccc1.C#CC=C.",
        "C1#CC1.C1CCC2(CC1)CCCCC2.C1C2CC3CC1CC(C2)C3.C1CC1.C1CCC1.",
        "C1CCCC1.C1CCCCCC1.C1CCCCCCC1.C1CCCCCCCC1.C1CCCCCCCCCC1.CCCC.CS"
    );
    const EXPECTED: [u32; 42] = [
        93, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 105, 13, 2, 1, 82, 3, 1, 1, 10, 6, 4, 2, 1, 1, 15, 5, 1,
        1, 78, 4, 1, 2, 1, 1, 8, 1, 1, 1, 1, 11, 12,
    ];

    assert!(EXPECTED.iter().all(|&value| value != 0));
    let molecule = Molecule::from_smiles(SMILES).expect("complete MQN branch fixture");
    assert_eq!(molecule.mqns(true).unwrap(), EXPECTED.to_vec());
}

#[test]
fn labute_contributions_match_pinned_rdkit_bits_across_bond_and_hydrogen_branches() {
    const CASES: [(&str, &[u64]); 7] = [
        ("", &[0x0000_0000_0000_0000, 0x0000_0000_0000_0000]),
        (
            "CC",
            &[
                0x3ff4_1b85_1ccc_a042,
                0x401b_b1e8_2a1b_1454,
                0x401b_b1e8_2a1b_1454,
                0x402e_3558_cdb4_a85c,
            ],
        ),
        (
            "C=C",
            &[
                0x3ff4_1b85_1ccc_a042,
                0x401a_50d4_840e_2bfb,
                0x401a_50d4_840e_2bfb,
                0x402c_d445_27a7_c003,
            ],
        ),
        (
            "C#N",
            &[
                0x3ff4_c3cc_4ecd_cb56,
                0x401a_49de_40b5_ed92,
                0x4015_0c2d_4cba_cd34,
                0x402a_437f_5092_16ce,
            ],
        ),
        (
            "c1ccccc1",
            &[
                0x3ff0_87fd_7763_0132,
                0x4018_43f5_ba92_4bc9,
                0x4018_43f5_ba92_4bc9,
                0x4018_43f5_ba92_4bc9,
                0x4018_43f5_ba92_4bc9,
                0x4018_43f5_ba92_4bc9,
                0x4018_43f5_ba92_4bc9,
                0x4042_b738_37a8_d0e1,
            ],
        ),
        (
            "[H][H]",
            &[
                0x3ff7_c1bb_ef62_32e4,
                0x3ff7_c1bb_ef62_32e4,
                0x3ff7_c1bb_ef62_32e4,
                0x4011_d14c_f389_a62b,
            ],
        ),
        (
            "CCO",
            &[
                0x3ff4_2e34_49f8_93d3,
                0x401b_b1e8_2a1b_1454,
                0x401a_6d72_7738_75fb,
                0x4014_6d15_8473_e02a,
                0x4033_e5ff_4e11_63db,
            ],
        ),
    ];

    for (smiles, expected) in CASES {
        let molecule = Molecule::from_smiles(smiles).expect("Labute fixture");
        let contribution = molecule
            .labute_asa_contributions_with_params(true, true)
            .unwrap();
        let mut actual = Vec::with_capacity(contribution.atom_contributions.len() + 2);
        actual.push(contribution.hydrogen_contribution);
        actual.extend_from_slice(&contribution.atom_contributions);
        actual.push(contribution.asa);
        assert_eq!(bits(&actual), expected, "{smiles:?} Labute contributions");
        assert_eq!(
            bits(&hydrogen_first(
                molecule
                    .labute_asa_contributions_with_params(true, false)
                    .unwrap()
            )),
            expected[..expected.len() - 1],
            "{smiles:?} hydrogen-first helper"
        );
        assert_eq!(
            molecule
                .labute_asa_with_params(true, false)
                .unwrap()
                .to_bits(),
            *expected.last().unwrap(),
            "{smiles:?} Labute ASA"
        );
    }
}

#[test]
fn slogp_and_smr_vsa_match_pinned_rdkit_default_and_custom_vectors() {
    const DEFAULT_CASES: [(&str, &[u64], &[u64]); 3] = [
        ("", &[0; 12], &[0; 10]),
        (
            "CCO",
            &[
                0,
                0x4027_6d43_fdd6_2b12,
                0,
                0,
                0x401b_b1e8_2a1b_1454,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ],
            &[
                0x4014_6d15_8473_e02a,
                0,
                0,
                0,
                0x401b_b1e8_2a1b_1454,
                0x401a_6d72_7738_75fb,
                0,
                0,
                0,
                0,
            ],
        ),
        (
            "CC(=O)N",
            &[
                0x4016_ef46_86f2_339f,
                0x4017_a0f3_b914_a2a9,
                0x4013_2d9b_27d4_2d79,
                0,
                0x401b_b1e8_2a1b_1454,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ],
            &[
                0x4013_2d9b_27d4_2d79,
                0,
                0,
                0x4016_ef46_86f2_339f,
                0x401b_b1e8_2a1b_1454,
                0,
                0,
                0,
                0,
                0x4017_a0f3_b914_a2a9,
            ],
        ),
    ];

    for (smiles, expected_slogp, expected_smr) in DEFAULT_CASES {
        let molecule = Molecule::from_smiles(smiles).expect("VSA fixture");
        let slogp = molecule.slogp_vsa_with_params(None, true).unwrap();
        let smr = molecule.smr_vsa_with_params(None, true).unwrap();
        assert_eq!(bits(&slogp), expected_slogp, "{smiles:?} SlogP-VSA");
        assert_eq!(bits(&smr), expected_smr, "{smiles:?} SMR-VSA");
        for (index, expected) in expected_slogp.iter().copied().enumerate() {
            assert_eq!(
                (SLOGP_SCALARS[index])(&molecule).unwrap().to_bits(),
                expected,
                "{smiles:?} SlogP-VSA scalar {}",
                index + 1
            );
        }
        for (index, expected) in expected_smr.iter().copied().enumerate() {
            assert_eq!(
                (SMR_SCALARS[index])(&molecule).unwrap().to_bits(),
                expected,
                "{smiles:?} SMR-VSA scalar {}",
                index + 1
            );
        }
    }

    const CUSTOM_CASES: [(&[f64], &[u64], &[u64]); 3] = [
        (
            &[-0.2, 0.0, 0.2],
            &[0x4027_6d43_fdd6_2b12, 0, 0x401b_b1e8_2a1b_1454, 0],
            &[0, 0, 0, 0x4032_a31c_0971_da9e],
        ),
        (
            &[0.0, 0.0],
            &[0x4027_6d43_fdd6_2b12, 0, 0x401b_b1e8_2a1b_1454],
            &[0, 0, 0x4032_a31c_0971_da9e],
        ),
        (
            &[2.0, 0.0, 1.0],
            &[0x4027_6d43_fdd6_2b12, 0, 0x401b_b1e8_2a1b_1454, 0],
            &[0, 0, 0x4014_6d15_8473_e02a, 0x402b_0fad_50a9_c528],
        ),
    ];
    for (bins, expected_slogp, expected_smr) in CUSTOM_CASES {
        let molecule = Molecule::from_smiles("CCO").expect("custom VSA fixture");
        assert_eq!(
            bits(&molecule.slogp_vsa_with_params(Some(bins), true).unwrap()),
            expected_slogp
        );
        assert_eq!(
            bits(&molecule.smr_vsa_with_params(Some(bins), true).unwrap()),
            expected_smr
        );
    }
}

#[test]
fn qed_aromatic_tellurium_alert_matches_rdkit() {
    let mol = Molecule::from_smiles("c1cc[te]c1").expect("tellurophene must parse");
    assert_eq!(mol.qed().unwrap().to_bits(), 0x3fe0827cd08382b8);
}

#[test]
fn qed_preserves_an_isolated_proton_like_rdkit_remove_hs() {
    let mol = Molecule::from_smiles("C.[H+]").expect("disconnected proton must parse");

    assert_eq!(mol.qed().unwrap().to_bits(), 0x3fd72f3ad7393056);
}

#[test]
fn vsa_force_controls_shared_labute_cache_without_changing_vector_shape() {
    let molecule = Molecule::from_smiles("CC").expect("VSA force fixture");
    let _ = molecule
        .labute_asa_contributions_with_params(false, false)
        .unwrap();

    let cached_slogp = molecule.slogp_vsa_with_params(None, false).unwrap();
    assert_eq!(cached_slogp.len(), 12);
    assert_eq!(cached_slogp[4].to_bits(), 0x402b_ca6e_1564_c404);
    let forced_slogp = molecule.slogp_vsa_with_params(None, true).unwrap();
    assert_eq!(forced_slogp.len(), 12);
    assert_eq!(forced_slogp[4].to_bits(), 0x402b_b1e8_2a1b_1454);

    let molecule_smr = Molecule::from_smiles("CC").expect("SMR-VSA force fixture");
    let _ = molecule_smr
        .labute_asa_contributions_with_params(false, false)
        .unwrap();
    let cached_smr = molecule_smr.smr_vsa_with_params(None, false).unwrap();
    assert_eq!(cached_smr.len(), 10);
    assert_eq!(cached_smr[4].to_bits(), 0x402b_ca6e_1564_c404);
    let forced_smr = molecule_smr.smr_vsa_with_params(None, true).unwrap();
    assert_eq!(forced_smr.len(), 10);
    assert_eq!(forced_smr[4].to_bits(), 0x402b_b1e8_2a1b_1454);
}
