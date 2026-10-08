//! Native rustdoc checks of public feature gates. These snippets deliberately
//! compile as external users, not as internal operation implementations.
//! Cargo runs them for the selected feature set, without nested Cargo commands.

macro_rules! absent_methods {
    ($feature:literal; $($method:ident),+ $(,)?) => {$(
        #[cfg(not(feature = $feature))]
        #[doc = concat!(
            "Requires ", $feature, ".\n\n",
            "~~~compile_fail,E0599\n",
            "let _ = cosmolkit::Molecule::", stringify!($method), ";\n",
            "~~~\n"
        )]
        const $method: () = ();
    )+};
}

macro_rules! absent_types {
    ($feature:literal; $($ty:ident),+ $(,)?) => {$(
        #[cfg(not(feature = $feature))]
        #[doc = concat!(
            "Requires ", $feature, ".\n\n",
            "~~~compile_fail,E0432\n",
            "use cosmolkit::", stringify!($ty), ";\n",
            "~~~\n"
        )]
        const $ty: () = ();
    )+};
}

absent_methods!("cap-smiles"; from_smiles);
absent_methods!("cap-io"; from_sdf);
absent_methods!("cap-hydrogens"; with_hydrogens, add_hydrogens_);
absent_methods!("cap-kekulize"; with_kekulized_bonds);
absent_methods!("cap-sanitize"; sanitize);
absent_methods!("cap-aromaticity"; with_assigned_aromaticity);
absent_methods!("cap-stereo"; potential_stereo);
absent_methods!("cap-rings"; with_assigned_rings);
absent_methods!("cap-valence"; with_assigned_valence);
absent_methods!("cap-depict"; with_2d_coordinates, to_svg, to_png);
absent_types!("cap-depict"; DrawingError);
absent_methods!("cap-forcefields";
    with_uff_optimized, with_uff_optimized_with_params,
    with_uff_optimized_confs, with_uff_optimized_confs_with_params,
    uff_has_all_molecule_params,
);
absent_types!("cap-forcefields";
    UffParameterError, UffParameterErrorKind, UffParameterQueryError, ForceFieldError,
);
absent_methods!("cap-descriptors";
    molecular_weight, num_heavy_atoms, total_atom_count, lipinski_hba,
    lipinski_hbd, fraction_csp3, num_rings, num_heterocycles,
    num_aromatic_rings, num_saturated_rings, num_aliphatic_rings,
    num_aromatic_heterocycles, num_aromatic_carbocycles,
    num_aliphatic_heterocycles, num_aliphatic_carbocycles,
    num_saturated_heterocycles, num_saturated_carbocycles,
);
absent_types!("cap-descriptors"; DescriptorError, DescriptorReadError);

/// The implementation module is private even when the capability is enabled.
/// ~~~compile_fail,E0603
/// use cosmolkit::forcefields::UffParameterQueryError;
/// ~~~
#[cfg(feature = "cap-forcefields")]
const PRIVATE_FORCEFIELDS_MODULE: () = ();

// Check each removed symbol independently so one import cannot mask another.
macro_rules! removed_items {
    ($($item:ident),+ $(,)?) => {$(
        #[doc = concat!(
            "This legacy root symbol is not public.\n\n",
            "~~~compile_fail,E0432\nuse cosmolkit::", stringify!($item), ";\n~~~\n"
        )]
        const $item: () = ();
    )+};
}
mod removed {
    removed_items!(
        ForceFieldOptions,
        mmff_has_all_molecule_params,
        mmff_optimize,
        uff_has_all_molecule_params,
    );
}

absent_types!("cap-valence"; ValenceParams);
absent_types!("cap-rings"; RingSearchParams);
