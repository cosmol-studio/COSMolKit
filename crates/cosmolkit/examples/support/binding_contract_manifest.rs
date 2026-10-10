//! Export of the linked, cfg-filtered canonical registry.
//! Shared by Python argument projection, stub generation and WASM checks; never parse source
//! text or maintain another list of API names here.
use cosmolkit::{
    BINDING_CONTRACT, BINDING_CONTRACT_KEYWORDS, BINDING_CONTRACT_PROPERTIES,
    BINDING_CONTRACT_PYTHON_ADAPTERS, BINDING_CONTRACT_PYTHON_ALIASES,
    BINDING_CONTRACT_PYTHON_COLLECTIONS, BindingDefault, BindingItem, BindingOwner,
    BindingPropertyAccess, BindingTypeRole,
};
use serde_json::{Value, json};

fn compact(ty: &str) -> String {
    ty.split_whitespace().collect()
}

pub fn manifest() -> Value {
    let entries = BINDING_CONTRACT.iter().map(|entry| {
        let constructor = if entry.type_role == Some(BindingTypeRole::Parameter) {
            BINDING_CONTRACT.iter().find(|candidate| {
                compact(candidate.rust_path).ends_with("::new")
                    && candidate.callable.is_some_and(|callable| {
                        callable.receiver.is_none()
                            && compact(callable.output_type) == compact(entry.rust_path)
                    })
            })
        } else {
            None
        };
        let fields = constructor.and_then(|row| row.callable).map(|callable| {
            callable.parameters.iter().map(|parameter| json!({
                "name": parameter.name,
                "type": parameter.type_name,
                "default": match parameter.default {
                    BindingDefault::Required => Value::Null,
                    BindingDefault::Value(value) => json!(value),
                },
            })).collect::<Vec<_>>()
        });
        json!({
            "semantic_id": entry.semantic_id,
            "rust_path": entry.rust_path,
            "python_name": entry.python_name,
            "python_native": entry.python_native,
            "python_fields": entry.python_configuration.map(|fields| fields.iter().map(|field| json!({
                "name": field.name, "type": field.type_name,
                "rust_field": field.rust_field, "aliases": field.aliases,
                "callback": field.callback,
                "default": match field.default {
                    BindingDefault::Required => Value::Null,
                    BindingDefault::Value(value) => json!(value),
                },
            })).collect::<Vec<_>>()),
            "javascript_name": entry.javascript_name,
            "python_property": match entry.python_property {
                None => None,
                Some(BindingPropertyAccess::Getter) => Some("getter"),
                Some(BindingPropertyAccess::Setter) => Some("setter"),
            },
            "feature": entry.feature,
            "required_capabilities": entry.required_capabilities,
            // Native archive APIs are excluded by the public platform contract.
            "platform": if entry.feature == "cap-serialization" { "native" } else { "all" },
            "item": match entry.item {
                BindingItem::Callable => "callable",
                BindingItem::Type => "type",
            },
            "owner": match entry.owner {
                BindingOwner::Module => "module",
                BindingOwner::Molecule => "molecule",
                BindingOwner::Type => "type",
            },
            "role": entry.type_role.map(|role| match role {
                BindingTypeRole::Parameter => "parameter",
                BindingTypeRole::ParameterSelector => "parameter_selector",
                BindingTypeRole::Value => "value",
                BindingTypeRole::Result => "result",
                BindingTypeRole::Error => "error",
            }),
            "constructor": constructor.map(|row| row.semantic_id),
            "fields": fields,
            "parameters": entry.callable.map(|callable| {
                callable.parameters.iter().map(|parameter| json!({
                    "name": parameter.name,
                    "type": parameter.type_name,
                    "default": match parameter.default {
                        BindingDefault::Required => Value::Null,
                        BindingDefault::Value(value) => json!(value),
                    },
                })).collect::<Vec<_>>()
            }),
            "output": entry.callable.map(|callable| callable.output_type),
            "receiver": entry.callable.and_then(|callable| callable.receiver).map(|value| format!("{value:?}")),
            "properties": BINDING_CONTRACT_PROPERTIES.iter()
                .filter(|field| field.type_semantic_id == entry.semantic_id)
                .map(|field| json!({"name": field.name, "type": field.output_type}))
                .collect::<Vec<_>>(),
        })
    }).collect::<Vec<_>>();
    json!({
        "python_collections": BINDING_CONTRACT_PYTHON_COLLECTIONS.iter().map(|(name, output)| json!({
            "name": name, "rust_output": output,
        })).collect::<Vec<_>>(),
        "entries": entries,
        "keywords": BINDING_CONTRACT_KEYWORDS.iter().map(|row| json!({
            "semantic_id": row.semantic_id,
            "constructor": row.parameters_semantic_id,
            "target": row.target_semantic_id,
        })).collect::<Vec<_>>(),
        "python_adapters": BINDING_CONTRACT_PYTHON_ADAPTERS.iter().map(|row| json!({
            "type_semantic_id": row.type_semantic_id,
            "name": row.name,
            "targets": row.targets,
        })).collect::<Vec<_>>(),
        "python_aliases": BINDING_CONTRACT_PYTHON_ALIASES.iter().map(|(name, target)| json!({
            "name": name, "target": target,
        })).collect::<Vec<_>>(),
    })
}
