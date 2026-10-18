//! Read all canonical Python projections enabled in this compiled facade.
use cosmolkit as ck;
use serde_json::json;
use std::collections::BTreeMap;

fn normalized_path(path: &str) -> String {
    path.split_whitespace()
        .collect::<String>()
        .replace("::<", "<")
}

fn base_path(path: &str) -> String {
    let mut depth = 0usize;
    normalized_path(path)
        .chars()
        .filter(|&ch| match ch {
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

fn owner_type(path: &str) -> Option<&'static ck::BindingContractEntry> {
    let types = || {
        ck::BINDING_CONTRACT
            .iter()
            .filter(|candidate| candidate.item == ck::BindingItem::Type)
    };
    if let Some(exact) =
        types().find(|candidate| normalized_path(candidate.rust_path) == normalized_path(path))
    {
        return Some(exact);
    }
    // Rust member paths omit a type's lifetime/default generic arguments. Only
    // accept an unambiguous registered owner; never guess between instantiations.
    let mut candidates =
        types().filter(|candidate| base_path(candidate.rust_path) == base_path(path));
    let first = candidates.next()?;
    candidates.next().is_none().then_some(first)
}

fn path_for(entry: &ck::BindingContractEntry) -> Result<String, String> {
    if entry.item == ck::BindingItem::Type || entry.owner == ck::BindingOwner::Module {
        return Ok(entry.python_name.to_owned());
    }
    let rust_path = normalized_path(entry.rust_path);
    let owner_path = rust_path
        .rsplit_once("::")
        .ok_or_else(|| format!("No Rust owner for {}", entry.semantic_id))?
        .0;
    let owner = owner_type(owner_path)
        .or_else(|| {
            entry.callable.as_ref().and_then(|callable| {
                // A registered bound factory may call a facade module function;
                // its static, value-returning signature identifies the owner type.
                (callable.kind == ck::BindingKind::Static
                    && callable.state_model == ck::StateModel::ValueReturning)
                    .then(|| owner_type(callable.output_type))
                    .flatten()
            })
        })
        .ok_or_else(|| format!("No canonical owner type for {}", entry.semantic_id))?;
    Ok(format!("{}.{}", owner.python_name, entry.python_name))
}

fn main() -> Result<(), String> {
    if std::env::args().len() != 1 {
        return Err("Registry selection is not supported; all compiled records are checked".into());
    }
    let records: Vec<_> = ck::BINDING_CONTRACT
        .iter()
        .map(|entry| {
            // Only the canonical keyword contract combines explicit/short call forms.
            let public_entry = ck::BINDING_CONTRACT_KEYWORDS
                .iter()
                .find(|keyword| keyword.target_semantic_id == entry.semantic_id)
                .and_then(|keyword| {
                    ck::BINDING_CONTRACT
                        .iter()
                        .find(|candidate| candidate.semantic_id == keyword.semantic_id)
                })
                .unwrap_or(entry);
            let path = path_for(public_entry);
            json!({
                "id": entry.semantic_id,
                "path": path.as_ref().ok(), "error": path.as_ref().err(),
                "kind": if entry.item == ck::BindingItem::Type { "class" }
                    else if entry.python_property.is_some() { "property" } else { "callable" },
                "writable": entry.python_property == Some(ck::BindingPropertyAccess::Setter),
                "native": entry.python_native,
                "error_type": entry.type_role == Some(ck::BindingTypeRole::Error),
                "read_accessor": entry.callable.as_ref().is_some_and(|callable|
                    callable.kind == ck::BindingKind::Instance
                    && callable.receiver == Some(ck::BindingReceiver::Shared)
                    && callable.parameters.is_empty()
                    && callable.state_model == ck::StateModel::ReadOnly
                    && callable.error_type.is_none()),
                "configuration": entry.python_configuration.map(|fields|
                    fields.iter().map(|field| field.name).collect::<Vec<_>>()),
                "configuration_schema": entry.python_configuration.map(|fields|
                    fields.iter().map(|field| json!({
                        "name": field.name, "type": field.type_name,
                        "aliases": field.aliases,
                        "required": matches!(field.default, ck::BindingDefault::Required),
                        "default": match field.default {
                            ck::BindingDefault::Required => None,
                            ck::BindingDefault::Value(value) => Some(value),
                        },
                    })).collect::<Vec<_>>()),
            })
        })
        .collect();
    // A callable's declared error is already canonical type information. Derive
    // its direct projection, not an independently maintained error inventory.
    let mut declared_errors = BTreeMap::<String, Vec<&str>>::new();
    for entry in ck::BINDING_CONTRACT {
        if let Some(error) = entry
            .callable
            .as_ref()
            .and_then(|callable| callable.error_type)
        {
            declared_errors
                .entry(normalized_path(error))
                .or_default()
                .push(entry.semantic_id);
        }
    }
    let mut error_names = BTreeMap::<String, String>::new();
    let mut direct_errors = Vec::new();
    for (rust_path, declared_by) in declared_errors {
        let candidates: Vec<_> = ck::BINDING_CONTRACT
            .iter()
            .filter(|entry| {
                entry.item == ck::BindingItem::Type
                    && base_path(entry.rust_path) == base_path(&rust_path)
            })
            .collect();
        let registered = match candidates.as_slice() {
            [] => None,
            [entry] => Some(*entry),
            _ => return Err(format!("Ambiguous registered direct error: {rust_path}")),
        };
        let name = registered
            .map(|entry| entry.python_name)
            .unwrap_or_else(|| rust_path.rsplit("::").next().unwrap_or(&rust_path));
        if name.is_empty()
            || name.chars().enumerate().any(|(index, ch)| {
                !(ch == '_' || ch.is_ascii_alphabetic() || index > 0 && ch.is_ascii_digit())
            })
        {
            return Err(format!(
                "Direct error requires an explicit registered alias: {rust_path}"
            ));
        }
        if let Some(previous) = error_names.insert(name.to_owned(), rust_path.clone())
            && previous != rust_path
        {
            return Err(format!(
                "Ambiguous direct error projection {name}: {previous}, {rust_path}"
            ));
        }
        if ck::BINDING_CONTRACT.iter().any(|entry| {
            entry.item == ck::BindingItem::Type
                && entry.python_name == name
                && base_path(entry.rust_path) != base_path(&rust_path)
        }) {
            return Err(format!(
                "Direct error conflicts with a registered Python type: {name}"
            ));
        }
        direct_errors.push(json!({"path":name,"kind":"class","error_type":true,
            "rust_type":rust_path,"declared_by":declared_by,
            "native":registered.and_then(|entry| entry.python_native)}));
    }
    let conveniences: Vec<_> = ck::BINDING_CONTRACT_CONVENIENCES
        .iter()
        .map(|item| json!({"path": item.name, "kind": "callable"}))
        .collect();
    let owner_member = |owner_id: &str, name: &str, kind: &str| {
        let owner = ck::BINDING_CONTRACT
            .iter()
            .find(|entry| entry.semantic_id == owner_id);
        json!({
            "path": owner.map(|entry| format!("{}.{}", entry.python_name, name)),
            "kind": kind,
            "error": owner.is_none().then(|| format!("No canonical owner {owner_id} for {name}")),
        })
    };
    let adapters: Vec<_> = ck::BINDING_CONTRACT_PYTHON_ADAPTERS
        .iter()
        .map(|item| owner_member(item.type_semantic_id, item.name, "callable"))
        .collect();
    let properties: Vec<_> = ck::BINDING_CONTRACT_PROPERTIES
        .iter()
        .map(|item| owner_member(item.type_semantic_id, item.name, "property"))
        .collect();
    let collections: Vec<_> = ck::BINDING_CONTRACT_PYTHON_COLLECTIONS
        .iter()
        .map(|(name, _)| json!({"path": name, "kind": "class"}))
        .collect();
    println!(
        "{}",
        json!({"version": 1, "records": records,
        "conveniences": conveniences, "adapters": adapters,
        "properties": properties, "collections": collections,
        "direct_errors": direct_errors})
    );
    Ok(())
}
