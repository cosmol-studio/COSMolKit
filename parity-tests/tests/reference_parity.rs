//! Independent libtest functions; Cargo parallelizes them by default.
use cosmolkit_parity_tests::testing::run_registered;

macro_rules! tests {
    ($($key:ident),+ $(,)?) => {
        const KEYS: &[&str] = &[$(stringify!($key)),+];
        $(#[test] fn $key() { run_registered(stringify!($key), KEYS).unwrap(); })+
    };
}

tests!(
    fuzzy_and_fingerprint_pairs,
    fuzzy_or_fingerprint_pairs,
    smiles_read_smiles,
    sanitize_smiles,
    kekulize_smiles,
    molecular_weight_smiles,
    exact_molecular_weight_smiles,
    molecular_formula_smiles,
    num_heavy_atoms_smiles,
    total_atom_count_smiles,
    lipinski_hba_smiles,
    lipinski_hbd_smiles,
    fraction_csp3_smiles,
    add_hydrogens_smiles,
    remove_hydrogens_smiles,
    coordinates_2d_smiles,
    distance_matrix_smiles,
);
