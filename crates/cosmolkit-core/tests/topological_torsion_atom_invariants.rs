use cosmolkit_core::fingerprint::AtomPairAtomInvGenerator;
use cosmolkit_core::{AtomSpec, BondOrder, BondSpec, ChiralTag, Element, Hybridization, Molecule};

fn invariants(molecule: &Molecule, chirality: bool, correction: bool) -> Vec<u32> {
    let molecule = molecule
        .with_assigned_valence_strict(false)
        .expect("RDKit atom-code fixtures require an initialized property cache");
    AtomPairAtomInvGenerator::new(chirality, correction)
        .getAtomInvariants(&molecule)
        .expect("atom-pair invariants")
}

fn one_atom(spec: AtomSpec) -> Molecule {
    let mut builder = Molecule::builder();
    builder.add_atom(spec);
    builder.build().expect("one-atom molecule")
}

fn carbon_star(degree: usize) -> Molecule {
    let mut builder = Molecule::builder();
    let center = builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3));
    for _ in 0..degree {
        let leaf =
            builder.add_atom(AtomSpec::new(Element::H).with_hybridization(Hybridization::Sp3));
        builder
            .add_bond(BondSpec::new(center, leaf, BondOrder::Single))
            .expect("star bond");
    }
    builder.build().expect("carbon star")
}

#[test]
fn topological_torsion_atom_invariants_cover_every_element_bucket() {
    let elements = [
        Element::B,
        Element::C,
        Element::N,
        Element::O,
        Element::F,
        Element::SI,
        Element::P,
        Element::S,
        Element::CL,
        Element::AS,
        Element::SE,
        Element::BR,
        Element::SB,
        Element::TE,
        Element::I,
        Element::H,
    ];
    let mut builder = Molecule::builder();
    for element in elements {
        builder.add_atom(AtomSpec::new(element).with_hybridization(Hybridization::Sp3));
    }
    let molecule = builder.build().expect("element buckets");

    assert_eq!(
        invariants(&molecule, false, false),
        (0..16).map(|bucket| bucket * 32).collect::<Vec<_>>()
    );
}

#[test]
fn topological_torsion_atom_invariants_cover_aromatic_hybridization_and_explicit_hydrogen() {
    let aromatic = one_atom(
        AtomSpec::new(Element::C)
            .with_aromatic(true)
            .with_hybridization(Hybridization::Sp3),
    );
    assert_eq!(invariants(&aromatic, false, false), vec![40]);

    let explicit_hydrogen = one_atom(
        AtomSpec::new(Element::C)
            .with_explicit_hydrogens(3)
            .with_hybridization(Hybridization::Sp2),
    );
    assert_eq!(invariants(&explicit_hydrogen, false, false), vec![32]);

    let mut builder = Molecule::builder();
    let carbon = builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2));
    let oxygen = builder.add_atom(AtomSpec::new(Element::O).with_hybridization(Hybridization::Sp2));
    builder
        .add_bond(BondSpec::new(carbon, oxygen, BondOrder::Double))
        .expect("double bond");
    let double_bond = builder.build().expect("C=O molecule");
    assert_eq!(invariants(&double_bond, false, false), vec![41, 105]);
}

#[test]
fn topological_torsion_atom_invariants_cover_dative_direction_and_modulo_boundaries() {
    let mut builder = Molecule::builder();
    let donor = builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2));
    let acceptor =
        builder.add_atom(AtomSpec::new(Element::N).with_hybridization(Hybridization::Sp2));
    builder
        .add_bond(BondSpec::new(donor, acceptor, BondOrder::Dative))
        .expect("dative bond");
    let dative = builder.build().expect("dative molecule");
    assert_eq!(invariants(&dative, false, false), vec![33, 65]);

    assert_eq!(invariants(&carbon_star(6), false, false)[0], 38);
    assert_eq!(invariants(&carbon_star(7), false, false)[0], 32);
    assert_eq!(invariants(&carbon_star(8), false, false)[0], 33);

    let mut builder = Molecule::builder();
    let carbon = builder.add_atom(AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp));
    let metal = builder.add_atom(AtomSpec::new(Element::CR).with_hybridization(Hybridization::Sp));
    builder
        .add_bond(BondSpec::new(carbon, metal, BondOrder::Quadruple))
        .expect("quadruple bond");
    let quadruple = builder.build().expect("quadruple-bond molecule");
    assert_eq!(invariants(&quadruple, false, false)[0], 33);
}

#[test]
fn topological_torsion_atom_invariants_cover_r_s_and_missing_cip_state() {
    let r = one_atom(
        AtomSpec::new(Element::C)
            .with_hybridization(Hybridization::Sp3)
            .with_chiral_tag(ChiralTag::TetrahedralCw)
            .with_prop("_CIPCode", "R"),
    );
    let s = one_atom(
        AtomSpec::new(Element::C)
            .with_hybridization(Hybridization::Sp3)
            .with_chiral_tag(ChiralTag::TetrahedralCcw)
            .with_prop("_CIPCode", "S"),
    );
    let missing = one_atom(
        AtomSpec::new(Element::C)
            .with_hybridization(Hybridization::Sp3)
            .with_chiral_tag(ChiralTag::TetrahedralCw),
    );

    assert_eq!(invariants(&r, true, false), vec![544]);
    assert_eq!(invariants(&s, true, false), vec![1056]);
    assert_eq!(invariants(&missing, true, false), vec![32]);
    assert_eq!(missing.prop("_CIPComputed"), None);
    assert_eq!(missing.atoms()[0].prop("_CIPCode"), None);
}

#[test]
fn topological_torsion_atom_invariants_cover_correction_clone_and_json() {
    let boron = one_atom(AtomSpec::new(Element::B).with_hybridization(Hybridization::Sp3));
    let generator = AtomPairAtomInvGenerator::new(false, true);
    assert_eq!(
        generator
            .getAtomInvariants(&boron)
            .expect("corrected boron"),
        vec![u32::MAX - 1]
    );
    assert_eq!(
        generator.infoString(),
        "AtomPairInvariantGenerator topologicalTorsionCorrection=1"
    );
    assert_eq!(generator.clone(), generator);
    assert_eq!(
        generator.toJSON(),
        r#"{"type":"AtomPairAtomInvGenerator","includeChirality":"false","topologicalTorsionCorrection":"true"}"#
    );

    let mut restored = AtomPairAtomInvGenerator::default();
    restored
        .fromJSON(&generator.toJSON())
        .expect("generator JSON round trip");
    assert_eq!(restored, generator);
    restored
        .fromJSON(r#"{"includeChirality":true}"#)
        .expect("partial generator JSON");
    assert!(restored.include_chirality);
    assert!(restored.topological_torsion_correction);
}
