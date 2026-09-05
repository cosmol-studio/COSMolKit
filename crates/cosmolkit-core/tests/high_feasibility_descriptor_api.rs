use std::thread;

use cosmolkit_core::{
    Molecule, calc_chi_0, calc_chi_0n, calc_chi_0v, calc_chi_1, calc_chi_nv, calc_hall_kier_alpha,
    calc_hall_kier_alpha_with_contributions, calc_kappa_1, calc_kappa_2, calc_kappa_3, calc_labute_asa,
    calc_labute_asa_contributions, calc_lipinski_hba, calc_lipinski_hbd, calc_mqns, calc_num_atom_stereo_centers,
    calc_num_bridgehead_atoms, calc_num_heteroatoms, calc_num_rings, calc_num_spiro_atoms,
    calc_num_unspecified_atom_stereo_centers, calc_phi, calc_slogp_vsa, calc_slogp_vsa_1, calc_slogp_vsa_12,
    calc_slogp_vsa_with_bins, calc_smr_vsa, calc_smr_vsa_1, calc_smr_vsa_10, calc_smr_vsa_with_bins,
};

#[derive(Debug, Clone, PartialEq)]
struct DescriptorSnapshot {
    chi: [u64; 5],
    shape: [u64; 5],
    counts: [u32; 7],
    mqns: [u32; 42],
    labute: u64,
    slogp_vsa: [u64; 12],
    smr_vsa: [u64; 10],
}

fn bits<const N: usize>(values: [f64; N]) -> [u64; N] {
    values.map(f64::to_bits)
}

fn snapshot(molecule: &Molecule) -> DescriptorSnapshot {
    DescriptorSnapshot {
        chi: [
            calc_chi_0(molecule).to_bits(),
            calc_chi_1(molecule).to_bits(),
            calc_chi_0v(molecule, false).unwrap().to_bits(),
            calc_chi_0n(molecule, false).unwrap().to_bits(),
            calc_chi_nv(molecule, 3, false).unwrap().to_bits(),
        ],
        shape: [
            calc_hall_kier_alpha(molecule).to_bits(),
            calc_kappa_1(molecule).to_bits(),
            calc_kappa_2(molecule).unwrap().to_bits(),
            calc_kappa_3(molecule).unwrap().to_bits(),
            calc_phi(molecule).unwrap().to_bits(),
        ],
        counts: [
            calc_lipinski_hba(molecule).unwrap(),
            calc_lipinski_hbd(molecule).unwrap(),
            calc_num_heteroatoms(molecule).unwrap(),
            calc_num_rings(molecule).unwrap(),
            calc_num_spiro_atoms(molecule).unwrap(),
            calc_num_bridgehead_atoms(molecule).unwrap(),
            calc_num_atom_stereo_centers(molecule).unwrap()
                + calc_num_unspecified_atom_stereo_centers(molecule).unwrap(),
        ],
        mqns: calc_mqns(molecule).unwrap(),
        labute: calc_labute_asa(molecule, true, false).to_bits(),
        slogp_vsa: bits(calc_slogp_vsa(molecule, false).unwrap()),
        smr_vsa: bits(calc_smr_vsa(molecule, false).unwrap()),
    }
}

#[test]
fn rust_descriptor_facade_composes_vectors_contributions_clones_and_custom_bins() {
    let molecule = Molecule::from_smiles("CC(O)c1ccncc1").unwrap();
    let structural_snapshot = molecule.clone();

    let expected = snapshot(&molecule);
    assert_eq!(calc_slogp_vsa_1(&molecule).unwrap().to_bits(), expected.slogp_vsa[0]);
    assert_eq!(calc_slogp_vsa_12(&molecule).unwrap().to_bits(), expected.slogp_vsa[11]);
    assert_eq!(calc_smr_vsa_1(&molecule).unwrap().to_bits(), expected.smr_vsa[0]);
    assert_eq!(calc_smr_vsa_10(&molecule).unwrap().to_bits(), expected.smr_vsa[9]);

    let custom_bins = [-0.2, 0.0, 0.25, 0.25, 0.8];
    assert_eq!(
        calc_slogp_vsa_with_bins(&molecule, &custom_bins, true).unwrap().len(),
        custom_bins.len() + 1
    );
    assert_eq!(
        calc_smr_vsa_with_bins(&molecule, &custom_bins, true).unwrap().len(),
        custom_bins.len() + 1
    );

    let hall_kier = calc_hall_kier_alpha_with_contributions(&molecule);
    assert_eq!(hall_kier.atom_contributions.len(), molecule.num_atoms());
    assert_eq!(
        hall_kier
            .atom_contributions
            .iter()
            .fold(0.0, |sum, contribution| sum + contribution)
            .to_bits(),
        hall_kier.alpha.to_bits()
    );
    let labute = calc_labute_asa_contributions(&molecule, true, true);
    assert_eq!(labute.atom_contributions.len(), molecule.num_atoms());
    assert_eq!(
        (labute
            .atom_contributions
            .iter()
            .fold(0.0, |sum, contribution| sum + contribution)
            + labute.hydrogen_contribution)
            .to_bits(),
        labute.asa.to_bits()
    );

    let clone = molecule.clone();
    let _ = calc_smr_vsa(&clone, true).unwrap();
    let _ = calc_chi_nv(&clone, 5, true).unwrap();
    assert_eq!(snapshot(&clone), expected);
    assert_eq!(snapshot(&molecule), expected);
    assert_eq!(molecule, structural_snapshot);
}

#[test]
fn rust_descriptor_facade_is_deterministic_for_repeated_mixed_and_parallel_reads() {
    let molecule = Molecule::from_smiles("CC(O)c1ccncc1").unwrap();
    let expected = snapshot(&molecule);

    for _ in 0..4 {
        let _ = calc_mqns(&molecule).unwrap();
        let _ = calc_slogp_vsa(&molecule, false).unwrap();
        assert_eq!(snapshot(&molecule), expected);
    }

    thread::scope(|scope| {
        let handles = (0..8).map(|_| scope.spawn(|| snapshot(&molecule))).collect::<Vec<_>>();
        for handle in handles {
            assert_eq!(handle.join().unwrap(), expected);
        }
    });
}
