#[path = "../src/binding.rs"]
mod binding;
#[path = "../src/status.rs"]
mod status;

fn declaration(adapters: &str) -> String {
    format!(
        r#"pub static API = [
        {{semantic_id:"types.Value",item:type,owner:type_,rust:crate::Value,
          python:"Value",javascript:"Value",feature:"runtime",role:value,
          python_adapters:[{adapters}]}},
        {{semantic_id:"Value.new",item:callable,owner:type_,rust:crate::Value::new,
          python:"new",javascript:"new",feature:"runtime",kind:static_,parameters:[],
          output:crate::Value,error:none,state:value_returning,operation:none,
          signature:fn()->crate::Value}}
    ];"#
    )
}

#[test]
fn adapter_metadata_keeps_real_rust_assertions() {
    let output = binding::expand_binding_contract(
        declaration(r#"{name:from_object,targets:["Value.new"]}"#)
            .parse()
            .unwrap(),
    )
    .unwrap()
    .to_string();
    assert!(output.contains("API_PYTHON_ADAPTERS"));
    assert!(output.contains("BindingPythonAdapterContract"));
    assert!(output.contains("from_object"));
    assert!(output.contains("__BINDING_ASSERT_API_1"));
    assert!(!output.contains("Value :: from_object"));
}

#[test]
fn adapter_targets_cannot_be_fake_empty_repeated_or_canonical_aliases() {
    for (adapters, message) in [
        (
            r#"{name:from_object,targets:[]}"#,
            "requires registered Rust targets",
        ),
        (
            r#"{name:from_object,targets:["Value.absent"]}"#,
            "not a registered Rust callable",
        ),
        (
            r#"{name:from_object,targets:["Value.new","Value.new"]}"#,
            "duplicate Python adapter target",
        ),
        (
            r#"{name:new,targets:["Value.new"]}"#,
            "duplicates a canonical Rust callable",
        ),
        (
            r#"{name:from_object,targets:["Value.new"]},{name:from_object,targets:["Value.new"]}"#,
            "duplicate Python adapter",
        ),
    ] {
        let error =
            binding::expand_binding_contract(declaration(adapters).parse().unwrap()).unwrap_err();
        assert!(error.to_string().contains(message), "{error}");
    }
}

#[test]
fn adapter_availability_inherits_its_target_capability_gate() {
    let source = declaration(r#"{name:from_object,targets:["Value.new"]}"#)
        .replace(
            "{semantic_id:\"Value.new\"",
            "#[cfg(feature = \"cap-valence\")] {semantic_id:\"Value.new\"",
        )
        .replace(
            "javascript:\"new\",feature:\"runtime\"",
            "javascript:\"new\",feature:\"cap-valence\"",
        );
    let output = binding::expand_binding_contract(source.parse().unwrap())
        .unwrap()
        .to_string()
        .chars()
        .filter(|c| !c.is_whitespace())
        .collect::<String>();
    assert!(output.contains("#[cfg(feature=\"cap-valence\")]crate::BindingPythonAdapterContract"));
    assert_eq!(output.matches("#[cfg(feature=\"cap-valence\")]").count(), 3);
}
