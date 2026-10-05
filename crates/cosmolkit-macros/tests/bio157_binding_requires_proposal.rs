//! ROOT-authorized requirements grammar proposals; p1 acceptance remains pending.
#[path = "../src/binding.rs"]
mod binding;
#[path = "../src/status.rs"]
mod status;
fn entry(extra: &str) -> String {
    format!(
        r#"pub static API = [
        #[cfg(feature = "cap-io")]
        {{semantic_id:"Value.convert",item:callable,owner:type_,rust:crate::Value::convert,
          python:"convert",javascript:"convert",feature:"cap-io",{extra}
          kind:instance,parameters:[],output:usize,error:none,state:read_only,operation:none,
          signature:fn(&crate::Value)->usize}}
    ];"#
    )
}
fn generated(extra: &str) -> String {
    binding::expand_binding_contract(entry(extra).parse().unwrap())
        .unwrap()
        .to_string()
        .chars()
        .filter(|c| !c.is_whitespace())
        .collect()
}
fn rejected(source: String) -> String {
    binding::expand_binding_contract(source.parse().unwrap())
        .unwrap_err()
        .to_string()
}
#[test]
fn one_declaration_gates_the_row_and_signature_and_preserves_owner() {
    let out = generated("requires:[\"cap-bio\"],");
    let cfg = "#[cfg(all(feature=\"cap-io\",feature=\"cap-bio\"))]";
    assert_eq!(out.matches(cfg).count(), 2);
    assert!(out.contains(&format!("{cfg}crate::BindingContractEntry")));
    assert!(out.contains(&format!("{cfg}const__BINDING_ASSERT_API_0")));
    assert!(out.contains("feature:\"cap-io\",required_capabilities:&[\"cap-bio\"]"));
    assert!(!out.contains("structValue"));
}
#[test]
fn old_entries_keep_their_original_single_cfg_and_no_added_requirements() {
    let out = generated("");
    assert_eq!(out.matches("#[cfg(feature=\"cap-io\")]").count(), 2);
    assert!(out.contains("required_capabilities:&[]"));
    assert!(!out.contains("cfg(all("));
}
#[test]
fn requirements_reject_empty_duplicate_owner_and_non_cap_spellings() {
    for (extra, expected) in [
        ("requires:[],", "at least one capability"),
        ("requires:[\"\"],", "nonempty cap- capability spelling"),
        ("requires:[\"cap-io\"],", "repeat its owner"),
        (
            "requires:[\"cap-bio\",\"cap-bio\"],",
            "duplicate required capability",
        ),
        ("requires:[\"bio\"],", "cap- capability spelling"),
        ("requires:[\"cap-BIO\"],", "cap- capability spelling"),
        ("requires:[\"cap-\"],", "cap- capability spelling"),
        ("requires:[\"cap--bio\"],", "cap- capability spelling"),
        ("requires:[\"cap-bio-\"],", "cap- capability spelling"),
        ("requires:[\" cap-bio\"],", "cap- capability spelling"),
        (
            "requires:[\"cap-bio\"],requires:[\"cap-search\"],",
            "duplicate binding field",
        ),
    ] {
        assert!(rejected(entry(extra)).contains(expected), "{extra}");
    }
}
#[test]
fn requires_does_not_relax_original_compound_mismatch_or_extra_cfg_rejections() {
    let base = entry("requires:[\"cap-bio\"],");
    assert!(
        rejected(base.replace(
            "#[cfg(feature = \"cap-io\")]",
            "#[cfg(all(feature = \"cap-io\",feature = \"cap-bio\"))]"
        ))
        .contains("compound or non-feature")
    );
    assert!(
        rejected(base.replace(
            "#[cfg(feature = \"cap-io\")]",
            "#[cfg(feature = \"cap-bio\")]"
        ))
        .contains("cfg feature disagrees")
    );
    assert!(
        rejected(base.replace(
            "#[cfg(feature = \"cap-io\")]",
            "#[cfg(feature = \"cap-io\")] #[cfg(feature = \"cap-bio\")]"
        ))
        .contains("at most one cfg")
    );
}
