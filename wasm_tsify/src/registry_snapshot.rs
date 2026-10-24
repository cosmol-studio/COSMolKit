//! Inspection-only metadata from the registry compiled into this exact WASM.
use cosmolkit as ck;
use serde_json::json;
use wasm_bindgen::prelude::*;

fn base_path(path: &str) -> String {
    let mut depth = 0usize;
    path.split_whitespace()
        .collect::<String>()
        .replace("::<", "<")
        .chars()
        .filter(|&c| match c {
            '<' => {
                depth += 1;
                false
            }
            '>' => {
                depth -= 1;
                false
            }
            _ => depth == 0,
        })
        .collect()
}

fn owner(path: &str) -> Option<&'static ck::BindingContractEntry> {
    let mut candidates = ck::BINDING_CONTRACT.iter().filter(|entry| {
        entry.item == ck::BindingItem::Type && base_path(entry.rust_path) == base_path(path)
    });
    let first = candidates.next()?;
    candidates.next().is_none().then_some(first)
}

#[wasm_bindgen(js_name = __bindingRegistrySnapshot, skip_typescript)]
pub fn snapshot() -> String {
    let records: Vec<_> = ck::BINDING_CONTRACT.iter().map(|entry| {
        let is_type = entry.item == ck::BindingItem::Type;
        let module_function = match entry.javascript_projection {
            ck::BindingJavascriptProjection::ModuleFunction(name) => Some(name),
            _ => None,
        };
        let module = is_type || entry.owner == ck::BindingOwner::Module || module_function.is_some();
        let parent = entry.rust_path.rsplit_once("::").and_then(|(path, _)| owner(path))
            .or_else(|| entry.callable.as_ref().and_then(|call| {
                (call.kind == ck::BindingKind::Static && call.state_model == ck::StateModel::ValueReturning)
                    .then(|| owner(call.output_type)).flatten()
            }));
        let projection = if parent.is_some_and(|p| p.javascript_projection == ck::BindingJavascriptProjection::HostFilesystem) {
            ck::BindingJavascriptProjection::HostFilesystem
        } else { entry.javascript_projection };
        let (projection, native_type) = match projection {
            ck::BindingJavascriptProjection::Runtime => ("runtime", None),
            ck::BindingJavascriptProjection::ModuleFunction(_) => ("runtime", None),
            ck::BindingJavascriptProjection::TypeOnly => ("type_only", None),
            ck::BindingJavascriptProjection::NativeValue(value) => ("native", Some(value)),
            ck::BindingJavascriptProjection::ErrorDetails => ("error_details", None),
            ck::BindingJavascriptProjection::HostFilesystem => ("host_filesystem", None),
        };
        json!({
            "id": entry.semantic_id,
            "name": module_function.unwrap_or(entry.javascript_name),
            "owner": if module { None } else { parent.map(|p| p.javascript_name) },
            "unresolved_owner": !module && parent.is_none(),
            "kind": if is_type { "type" } else { "member" },
            "instance": !module && entry.callable.as_ref().is_some_and(|c| c.kind == ck::BindingKind::Instance),
            "projection": projection, "native_type": native_type,
        })
    }).collect();
    let conveniences: Vec<_> = ck::BINDING_CONTRACT_CONVENIENCES
        .iter()
        .filter(|entry| entry.javascript_target.is_some())
        .map(|entry| {
            json!({"id": entry.target, "name": entry.name,
            "owner": null, "kind": "member", "instance": false})
        })
        .collect();
    json!({"version": 2, "records": records, "conveniences": conveniences}).to_string()
}
