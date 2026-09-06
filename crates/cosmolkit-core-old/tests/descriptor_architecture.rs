use std::{
    fs,
    path::{Path, PathBuf},
};

use cosmolkit_core::{
    CrippenDescriptorValues, DescriptorResult, Molecule, NumRotatableBondsOptions,
};

fn core_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn rust_sources_below(root: &Path, sources: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(root)
        .unwrap_or_else(|error| panic!("read source directory {}: {error}", root.display()))
    {
        let path = entry.expect("read source directory entry").path();
        if path.is_dir() {
            rust_sources_below(&path, sources);
        } else if path.extension().is_some_and(|extension| extension == "rs") {
            sources.push(path);
        }
    }
}

fn definition_paths(definition: &str) -> Vec<String> {
    let root = core_root();
    let source_root = root.join("src");
    let mut sources = Vec::new();
    rust_sources_below(&source_root, &mut sources);
    sources.sort();

    let mut definitions = Vec::new();
    for path in sources {
        let source = fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("read source file {}: {error}", path.display()));
        for _ in source.match_indices(definition) {
            definitions.push(
                path.strip_prefix(&root)
                    .expect("source stays below crate root")
                    .to_string_lossy()
                    .replace('\\', "/"),
            );
        }
    }
    definitions
}

#[test]
fn descriptor_family_cores_have_single_owners() {
    let definitions = [
        (
            "fn find_all_paths_of_length_n(",
            &["src/chemistry/subgraph.rs"][..],
        ),
        ("fn rdkit_rb0(", &["src/chemistry/valence.rs"][..]),
        ("const RDKIT_RB0:", &["src/chemistry/valence.rs"][..]),
        ("fn add_hs_rdkit_rb0(", &[][..]),
        (
            "fn rdkit_crippen_atom_contribs(",
            &["src/properties/descriptors.rs"][..],
        ),
        (
            "fn hk_deltas(",
            &["src/properties/descriptors/connectivity.rs"][..],
        ),
        (
            "fn n_vals(",
            &["src/properties/descriptors/connectivity.rs"][..],
        ),
        (
            "fn calc_mqns_core(",
            &["src/properties/descriptors/mqn.rs"][..],
        ),
        (
            "fn get_labute_atom_contribs(",
            &["src/properties/descriptors/mol_surface.rs"][..],
        ),
        (
            "fn assign_contributions_to_bins(",
            &["src/properties/descriptors/mol_surface.rs"][..],
        ),
    ];

    for (definition, allowed_paths) in definitions {
        let paths = definition_paths(definition);
        let expected_count = usize::from(!allowed_paths.is_empty());
        assert_eq!(
            paths.len(),
            expected_count,
            "descriptor architecture requires exactly {expected_count} `{definition}` core, found {paths:?}"
        );
        for path in paths {
            assert!(
                allowed_paths.contains(&path.as_str()),
                "`{definition}` belongs to {allowed_paths:?}, found in {path}"
            );
        }
    }
}

#[test]
fn descriptor_facade_declares_every_family_boundary_once() {
    let facade = fs::read_to_string(core_root().join("src/properties/descriptors.rs"))
        .expect("read descriptor facade");
    for declaration in [
        "mod connectivity;",
        "mod lipinski;",
        "mod mol_surface;",
        "mod mqn;",
    ] {
        assert_eq!(
            facade.matches(declaration).count(),
            1,
            "descriptor facade must declare `{declaration}` exactly once"
        );
    }
}

#[test]
fn descriptor_facade_delegates_to_each_single_family_core() {
    let facade = fs::read_to_string(core_root().join("src/properties/descriptors.rs"))
        .expect("read descriptor facade");

    for delegation in [
        "connectivity::calc_chi_0(molecule)",
        "connectivity::calc_chi_nv(molecule, order, force)",
        "connectivity::$core(molecule, force)",
        "connectivity::calc_hall_kier_alpha(molecule, None)",
        "connectivity::calc_kappa_1(molecule)",
        "connectivity::calc_phi(molecule)",
        "lipinski::calc_lipinski_hba(molecule)",
        "lipinski::$core(molecule)",
        "lipinski::calc_num_spiro_atoms(molecule, None)",
        "lipinski::num_atom_stereo_centers(molecule)",
        "mqn::calc_mqns_core(molecule)",
        "mol_surface::calc_labute_asa(molecule, include_hydrogens, force)",
        "mol_surface::get_labute_atom_contribs(molecule, include_hydrogens, force)",
        "mol_surface::calc_slogp_vsa(molecule, None, force)",
        "mol_surface::calc_smr_vsa(molecule, None, force)",
    ] {
        assert!(
            facade.contains(delegation),
            "descriptor facade stopped delegating through `{delegation}`"
        );
    }

    assert_eq!(
        facade.matches("macro_rules! vsa_scalar_descriptor").count(),
        1,
        "VSA scalar projections require one facade macro"
    );
    assert!(
        facade.contains("Ok($vector(molecule, false)?[$index])"),
        "VSA scalar projections must index the shared vector result"
    );
}

#[test]
fn python_descriptor_facade_borrows_the_retained_molecule() {
    let python_source = fs::read_to_string(core_root().join("../../python/src/lib.rs"))
        .expect("read Python facade");
    let descriptor_start = python_source
        .find("fn calc_mol_wt(")
        .expect("Python descriptor facade start");
    let descriptor_end = python_source[descriptor_start..]
        .find("fn element_from_symbol")
        .map(|offset| descriptor_start + offset)
        .expect("Python descriptor facade end");
    let descriptor_facade = &python_source[descriptor_start..descriptor_end];

    assert!(
        descriptor_facade.contains("molecule: PyRef<'_, Molecule>"),
        "Python descriptor facade must borrow the retained Python Molecule"
    );
    assert!(
        !descriptor_facade.contains("molecule: &Molecule"),
        "`&Molecule` uses the pyclass FromPyObject clone path and discards observational cache writes"
    );
    assert!(
        !descriptor_facade.contains("molecule: Molecule"),
        "Python descriptor facade must not extract an owned Molecule clone"
    );
    assert_eq!(
        descriptor_facade
            .matches("fn $name(molecule: PyRef<'_, Molecule>)")
            .count(),
        3,
        "all scalar/count descriptor wrapper macros must borrow their Molecule"
    );
    assert_eq!(
        descriptor_facade
            .matches("fn $name(molecule: PyRef<'_, Molecule>, force: bool)")
            .count(),
        1,
        "the fixed-Chi wrapper macro must borrow its Molecule"
    );
}

#[test]
fn existing_descriptor_exports_keep_their_public_signatures() {
    use cosmolkit_core::properties::descriptors as module;

    let _: fn(&Molecule, bool) -> DescriptorResult<f64> = module::calc_mol_wt;
    let _: fn(&Molecule, bool) -> DescriptorResult<f64> = cosmolkit_core::calc_mol_wt;
    let _: fn(&Molecule, bool) -> DescriptorResult<f64> = module::calc_exact_mol_wt;
    let _: fn(&Molecule, bool) -> DescriptorResult<f64> = cosmolkit_core::calc_exact_mol_wt;
    let _: fn(&Molecule, bool, bool) -> DescriptorResult<String> = module::calc_mol_formula;
    let _: fn(&Molecule, bool, bool) -> DescriptorResult<String> = cosmolkit_core::calc_mol_formula;
    let _: fn(&Molecule) -> DescriptorResult<u32> = module::calc_num_hbd;
    let _: fn(&Molecule) -> DescriptorResult<u32> = cosmolkit_core::calc_num_hbd;
    let _: fn(&Molecule) -> DescriptorResult<u32> = module::calc_num_hba;
    let _: fn(&Molecule) -> DescriptorResult<u32> = cosmolkit_core::calc_num_hba;
    let _: fn(&Molecule) -> DescriptorResult<f64> = module::calc_fraction_csp3;
    let _: fn(&Molecule) -> DescriptorResult<f64> = cosmolkit_core::calc_fraction_csp3;
    let _: fn(&Molecule, bool, bool) -> DescriptorResult<CrippenDescriptorValues> =
        module::calc_crippen_descriptors;
    let _: fn(&Molecule, bool, bool) -> DescriptorResult<CrippenDescriptorValues> =
        cosmolkit_core::calc_crippen_descriptors;
    let _: fn(&Molecule, bool, bool) -> DescriptorResult<f64> = module::calc_tpsa;
    let _: fn(&Molecule, bool, bool) -> DescriptorResult<f64> = cosmolkit_core::calc_tpsa;
    let _: fn(&Molecule) -> DescriptorResult<u32> = module::calc_num_aromatic_rings;
    let _: fn(&Molecule) -> DescriptorResult<u32> = cosmolkit_core::calc_num_aromatic_rings;
    let _: fn(&Molecule, NumRotatableBondsOptions) -> DescriptorResult<u32> =
        module::calc_num_rotatable_bonds;
    let _: fn(&Molecule, NumRotatableBondsOptions) -> DescriptorResult<u32> =
        cosmolkit_core::calc_num_rotatable_bonds;
    let _: fn(&Molecule) -> DescriptorResult<f64> = module::calc_qed;
    let _: fn(&Molecule) -> DescriptorResult<f64> = cosmolkit_core::calc_qed;
}
