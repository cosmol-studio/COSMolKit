use std::collections::BTreeMap;

use cosmolkit::{
    BINDING_CONTRACT, BindingExposure, BindingItem, BindingKind, BindingOwner, BindingParity,
    BindingSupport, ChiralTag, Element, Molecule, SmilesError, SmilesParseParams,
};

#[test]
fn public_default_and_parameterized_constructors_are_directly_callable() {
    let short: fn(&str) -> Result<Molecule, SmilesError> = Molecule::from_smiles;
    let configured: for<'a, 'b> fn(
        &'a str,
        &'b SmilesParseParams,
    ) -> Result<Molecule, SmilesError> = Molecule::from_smiles_with_params;

    let molecule = short("CCO").expect("default SMILES construction");
    assert_eq!(molecule.num_atoms(), 3);
    assert_eq!(molecule.topology().atoms[0].element(), Element::C);

    let empty = configured("", &SmilesParseParams::default()).expect("empty molecule");
    assert_eq!(empty.num_atoms(), 0);
}

#[test]
fn public_parameters_preserve_every_pinned_default() {
    let params = SmilesParseParams::default();
    assert!(params.sanitize);
    assert!(params.allow_cxsmiles);
    assert!(params.strict_cxsmiles);
    assert!(params.parse_name);
    assert!(params.remove_hydrogens);
    assert!(!params.skip_cleanup);
    assert!(!params.debug_parse);
    assert!(params.replacements.is_empty());
}

#[test]
fn public_zero_isotope_projects_to_absence_without_losing_nonzero_isotopes() {
    for (input, isotope) in [("[C]", None), ("[0C]", None), ("[13C]", Some(13))] {
        for sanitize in [false, true] {
            let molecule = Molecule::from_smiles_with_params(
                input,
                &SmilesParseParams {
                    sanitize,
                    remove_hydrogens: false,
                    ..Default::default()
                },
            )
            .unwrap();
            assert_eq!(molecule.topology().atoms[0].isotope(), isotope, "{input}");
        }
    }
}

#[test]
fn public_cx_double_bond_stereo_finishes_after_sanitize_or_remove_hydrogens() {
    use cosmolkit::{AtomId, BondDirection as D, BondStereo as S};
    // Fixed RDKit 2026.03.1 MolFromSmiles: c/t/ctu through all four
    // sanitize/removeHs combinations, including both conformer dimensions.
    let cases = [
        (
            "CC=CC |c:1|",
            S::Cis,
            S::Z,
            [D::EndUpRight, D::None, D::EndDownRight],
        ),
        (
            "CC=CC |t:1|",
            S::Trans,
            S::E,
            [D::EndUpRight, D::None, D::EndUpRight],
        ),
        ("CC=CC |ctu:1|", S::Any, S::Any, [D::None; 3]),
        (
            "CC=CC |(0,1,;0,0,;1,0,;1,1,),c:1|",
            S::Cis,
            S::Z,
            [D::EndUpRight, D::None, D::EndDownRight],
        ),
        (
            "CC=CC |(0,1,1;0,0,1;1,0,1;1,1,1),c:1|",
            S::Cis,
            S::Z,
            [D::EndUpRight, D::None, D::EndDownRight],
        ),
    ];
    for (input, raw_stereo, final_stereo, final_directions) in cases {
        for sanitize in [false, true] {
            for remove_hydrogens in [false, true] {
                let params = SmilesParseParams {
                    sanitize,
                    remove_hydrogens,
                    ..Default::default()
                };
                let molecule = Molecule::from_smiles_with_params(input, &params).unwrap();
                let finalized = sanitize || remove_hydrogens;
                assert_eq!(
                    molecule.properties().prop("_needsDetectBondStereo"),
                    if finalized { None } else { Some("1") },
                    "{input} {params:?}"
                );
                let bonds = &molecule.topology().bonds;
                assert_eq!(
                    bonds[1].stereo(),
                    if finalized { final_stereo } else { raw_stereo },
                    "{input} {params:?}"
                );
                assert_eq!(
                    bonds[1].stereo_atoms(),
                    Some([AtomId::new(0), AtomId::new(3)])
                );
                assert_eq!(
                    bonds
                        .iter()
                        .map(|bond| bond.direction())
                        .collect::<Vec<_>>(),
                    if finalized {
                        final_directions
                    } else {
                        [D::None; 3]
                    },
                    "{input} {params:?}"
                );
            }
        }
    }
}

#[test]
fn public_cx_stereo_completion_uses_post_removal_atom_references() {
    use cosmolkit::{AtomId, BondDirection as D, BondStereo};
    // Source reference atom 0 is removed. Direction reconstruction must use
    // the remapped controls, never the pre-removal CX row indices.
    let molecule = Molecule::from_smiles("[H]C(C)=CC |c:2|").unwrap();
    assert_eq!(molecule.num_atoms(), 4);
    let bonds = &molecule.topology().bonds;
    assert_eq!(bonds[1].stereo(), BondStereo::E);
    assert_eq!(
        bonds[1].stereo_atoms(),
        Some([AtomId::new(1), AtomId::new(3)])
    );
    assert_eq!(
        bonds
            .iter()
            .map(|bond| bond.direction())
            .collect::<Vec<_>>(),
        [D::EndDownRight, D::None, D::EndUpRight]
    );
    assert_eq!(molecule.properties().prop("_needsDetectBondStereo"), None);
}

#[test]
fn parameterized_constructor_applies_replacements_and_post_parse_policy() {
    let retained = Molecule::from_smiles_with_params(
        "{W}",
        &SmilesParseParams {
            sanitize: false,
            remove_hydrogens: false,
            replacements: BTreeMap::from([("{W}".to_owned(), "[H]O[H]".to_owned())]),
            ..Default::default()
        },
    )
    .expect("unsanitized explicit-water graph");
    assert_eq!(retained.num_atoms(), 3);

    let removed = Molecule::from_smiles("[H]O[H]").expect("default RemoveHs hand-off");
    assert_eq!(removed.num_atoms(), 1);
    assert_eq!(removed.topology().atoms[0].element(), Element::O);
}

#[test]
fn public_constructor_preserves_cx_name_stereo_and_properties() {
    let molecule = Molecule::from_smiles_with_params(
        "F[C@](Cl)(Br)I |$fluoro;center;chloro;bromo;iodo$| sample name",
        &SmilesParseParams {
            sanitize: false,
            remove_hydrogens: false,
            ..Default::default()
        },
    )
    .expect("CX/name/stereo construction");

    assert_eq!(molecule.num_atoms(), 5);
    assert_eq!(
        molecule.topology().atoms[0].prop("atomLabel"),
        Some("fluoro")
    );
    assert_eq!(
        molecule.topology().atoms[1].chiral_tag(),
        ChiralTag::TetrahedralCcw
    );
    assert_eq!(molecule.properties().name(), Some("sample name"));
}

#[test]
fn public_error_retains_parse_category_and_offset() {
    let baseline = Molecule::from_smiles("CO").expect("baseline molecule");
    let error = Molecule::from_smiles("C?N").expect_err("invalid token");
    assert!(matches!(&error, SmilesError::Parse(_)));
    assert!(error.to_string().contains("'?' at byte 1"));
    assert_eq!(baseline.num_atoms(), 2, "failed construction is atomic");
}

#[test]
fn binding_contract_matches_public_smiles_surface() {
    let rows = BINDING_CONTRACT
        .iter()
        .filter(|row| {
            matches!(
                row.semantic_id,
                "types.SmilesParseParams"
                    | "types.SmilesError"
                    | "types.SmilesStereoError"
                    | "Molecule.from_smiles"
                    | "Molecule.from_smiles_with_params"
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(rows.len(), 5);
    assert_eq!(
        rows.iter().map(|row| row.semantic_id).collect::<Vec<_>>(),
        [
            "types.SmilesParseParams",
            "types.SmilesError",
            "types.SmilesStereoError",
            "Molecule.from_smiles",
            "Molecule.from_smiles_with_params",
        ]
    );

    for row in &rows {
        assert_eq!(row.feature, "smiles");
        assert_eq!(row.exposure, BindingExposure::Public);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
    }
    for row in &rows[..3] {
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.owner, BindingOwner::Type);
    }
    for row in &rows[3..] {
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Molecule);
        let callable = row.callable.expect("callable metadata");
        assert_eq!(callable.kind, BindingKind::Static);
        assert_eq!(callable.operation_semantic_id, None);
    }
    assert_eq!(rows[3].python_name, "from_smiles");
    assert_eq!(rows[3].javascript_name, "fromSmiles");
    assert_eq!(rows[4].python_name, "from_smiles_with_params");
    assert_eq!(rows[4].javascript_name, "fromSmilesWithParams");
}
