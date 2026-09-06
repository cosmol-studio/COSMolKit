use cosmolkit_core::{
    EmbedParameters, Molecule,
    chemistry::distgeom::{embed_molecule, rd_distgeom_get_exp_tors_helper_with_params},
};

#[test]
#[ignore = "debug helper for row-1 conformer parity investigation"]
fn debug_ethene_with_h_bounds_matrix() {
    let mol = Molecule::from_smiles("C=C")
        .expect("parse")
        .with_hydrogens()
        .expect("add hs");
    let bounds = mol.dg_bounds_matrix().expect("bounds");
    for row in bounds {
        println!("{row:?}");
    }
}

#[test]
#[ignore = "debug helper for row-1 conformer parity investigation"]
fn debug_ethene_with_h_embed_coords() {
    let mol = Molecule::from_smiles("C=C")
        .expect("parse")
        .with_hydrogens()
        .expect("add hs");
    let mut exp_only = EmbedParameters::etkdg_v3();
    exp_only.use_basic_knowledge = false;
    let mut basic_only = EmbedParameters::etkdg_v3();
    basic_only.use_exp_torsion_angle_prefs = false;
    for (label, mut params) in [
        ("ETDG", EmbedParameters::etdg()),
        ("ExpOnly", exp_only),
        ("BasicOnly", basic_only),
        ("ETKDGv3", EmbedParameters::etkdg_v3()),
    ] {
        params.random_seed = 61453;
        params.num_threads = 1;
        params.timeout = 10;
        let (embedded, status) = embed_molecule(&mol, &mut params).expect("embed");
        println!("preset={label} status={status}");
        for coord in embedded.conformers_3d()[0].coordinates() {
            println!("{coord:?}");
        }
    }
}

#[test]
#[ignore = "debug helper for row-1 conformer parity investigation"]
fn debug_ethene_with_h_etkdg_details() {
    let mol = Molecule::from_smiles("C=C")
        .expect("parse")
        .with_hydrogens()
        .expect("add hs");
    let params = EmbedParameters::etkdg_v3();
    let details = rd_distgeom_get_exp_tors_helper_with_params(&mol, &params).expect("details");
    println!("exp_torsion_atoms={:?}", details.exp_torsion_atoms);
    println!("exp_torsion_angles={:?}", details.exp_torsion_angles);
    println!("improper_atoms={:?}", details.improper_atoms);
    println!("atom_nums={:?}", details.atom_nums);
    println!(
        "bounds_mat_force_scaling={:?} constrained_atoms={:?}",
        details.bounds_mat_force_scaling, details.constrained_atoms
    );
    let bonds: Vec<_> = mol
        .bonds()
        .iter()
        .map(|bond| (bond.begin().index(), bond.end().index(), bond.order()))
        .collect();
    println!("mol_bonds={bonds:?}");
    let mut angles = Vec::new();
    for bondi in mol.bonds() {
        for j in (bondi.id().index() + 1)..mol.num_bonds() {
            let bondj = &mol.bonds()[j];
            let aid11 = bondi.begin().index() as i32;
            let aid12 = bondi.end().index() as i32;
            let aid21 = bondj.begin().index() as i32;
            let aid22 = bondj.end().index() as i32;
            if aid11 != aid21 && aid11 != aid22 && aid12 != aid21 && aid12 != aid22 {
                continue;
            }
            let mut tmp = vec![0; 4];
            if aid12 == aid21 {
                tmp[0] = aid11;
                tmp[1] = aid12;
                tmp[2] = aid22;
            } else if aid12 == aid22 {
                tmp[0] = aid11;
                tmp[1] = aid12;
                tmp[2] = aid21;
            } else if aid11 == aid21 {
                tmp[0] = aid12;
                tmp[1] = aid11;
                tmp[2] = aid22;
            } else if aid11 == aid22 {
                tmp[0] = aid12;
                tmp[1] = aid11;
                tmp[2] = aid21;
            }
            angles.push(tmp);
        }
    }
    println!("collect_bonds_and_angles_angles={angles:?}");
}

#[test]
#[ignore = "source-bisection probe for the two audited ETKDGv3 coordinate divergences"]
fn debug_audited_etkdgv3_exp_torsion_details() {
    for smiles in [
        "O=C1Nc2ccc(Cl)cc2C(c2ccccc2Cl)=NC1O",
        "CN1C2CCC1CC1(CN=C(c3cn(C)c4ccccc34)O1)C2",
    ] {
        let molecule = Molecule::from_smiles(smiles)
            .expect("parse")
            .with_hydrogens()
            .expect("add hs");
        let params = EmbedParameters::etkdg_v3();
        let details = rd_distgeom_get_exp_tors_helper_with_params(&molecule, &params)
            .expect("experimental torsion details");

        println!("smiles={smiles}");
        println!("exp_torsion_atoms={:?}", details.exp_torsion_atoms);
        println!("exp_torsion_angles={:?}", details.exp_torsion_angles);
        println!("improper_atoms={:?}", details.improper_atoms);
        println!("atom_nums={:?}", details.atom_nums);
    }
}
