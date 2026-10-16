//! Read all canonical Python projections enabled in this compiled facade.
use cosmolkit as ck;
use serde_json::json;

fn normalized_path(path: &str) -> String {
    path.split_whitespace()
        .collect::<String>()
        .replace("::<", "<")
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
    let owner = ck::BINDING_CONTRACT
        .iter()
        .find(|candidate| {
            candidate.item == ck::BindingItem::Type
                && normalized_path(candidate.rust_path) == owner_path
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
                "configuration": entry.python_configuration.map(|fields|
                    fields.iter().map(|field| field.name).collect::<Vec<_>>()),
            })
        })
        .collect();
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
        "properties": properties, "collections": collections})
    );
    Ok(())
}
