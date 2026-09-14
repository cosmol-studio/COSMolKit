#[path = "../src/result.rs"]
mod result;

fn expand(source: &str) -> syn::Result<proc_macro2::TokenStream> {
    result::expand(syn::parse_str(source).unwrap())
}

#[test]
fn pending_derive_maps_only_the_marked_field_and_retains_metadata() {
    for field_type in ["M", "Option<M>"] {
        let output = expand(&format!("struct Report<M = crate::Molecule> {{ count: usize, #[pending_molecule] molecule: {field_type} }}")).unwrap();
        let code = output.to_string().replace(' ', "");
        assert!(code.contains("Report<crate::PendingMolecule<Access>>"));
        assert!(code.contains("typeFinished=Report<crate::Molecule>"));
        assert!(code.contains("finalizer.resolve("));
        assert!(code.contains("count"));
        assert!(!code.contains("from_blocks"));
        syn::parse2::<syn::ItemImpl>(output).unwrap();
    }
}

#[test]
fn pending_derive_rejects_ambiguous_or_unrestricted_result_shapes() {
    for source in [
        "struct Report { molecule: Molecule }",
        "struct Report<M> { #[pending_molecule] molecule: M }",
        "struct Report<M: Clone = Molecule> { #[pending_molecule] molecule: M }",
        "struct Report<M = Molecule, N = Molecule> { #[pending_molecule] molecule: M, other: N }",
        "struct Report<M = Other> { #[pending_molecule] molecule: M }",
        "enum Report<M = Molecule> { Value(M) }",
        "struct Report<M = Molecule>(M);",
        "struct Report<M = Molecule> { molecule: M }",
        "struct Report<M = Molecule> { #[pending_molecule] a: M, #[pending_molecule] b: M }",
        "struct Report<M = Molecule> { #[pending_molecule] molecule: Vec<M> }",
        "struct Report<M = Molecule> { #[pending_molecule] molecule: Molecule }",
        "struct Report<M = Molecule> { #[pending_molecule] molecule: Option<Molecule> }",
    ] {
        assert!(expand(source).is_err(), "unexpectedly accepted {source}");
    }
}
