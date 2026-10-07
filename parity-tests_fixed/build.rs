use std::{env, fs, path::PathBuf};

fn main() {
    println!("cargo:rerun-if-env-changed=PARITY_CORPUS");
    println!("cargo:rerun-if-env-changed=PARITY_SPECIAL");
    for lane in ["corpus", "special"] {
        let selected = format!("expected/{lane}/selection.json");
        println!("cargo:rerun-if-changed={selected}");
        let environment = format!("PARITY_{}", lane.to_uppercase());
        let fallback = if lane == "corpus" {
            "smiles_5000"
        } else {
            "all"
        };
        let name = env::var(environment)
            .ok()
            .or_else(|| {
                fs::read(&selected)
                    .ok()
                    .and_then(|v| serde_json::from_slice::<String>(&v).ok())
            })
            .unwrap_or_else(|| fallback.into());
        println!(
            "cargo:rustc-env=PARITY_DEFAULT_{}={name}",
            lane.to_uppercase()
        );
        if lane == "corpus" {
            for format in ["smiles", "bio"] {
                println!("cargo:rustc-check-cfg=cfg(parity_corpus_{format})");
            }
            let manifest = format!("testdata/corpora/{name}.json");
            println!("cargo:rerun-if-changed={manifest}");
            let value: serde_json::Value =
                serde_json::from_slice(&fs::read(manifest).expect("unknown corpus manifest"))
                    .unwrap();
            let format = value["format"].as_str().unwrap();
            assert!(["smiles", "bio"].contains(&format));
            println!("cargo:rustc-cfg=parity_corpus_{format}");
        } else {
            assert!(
                [
                    "all",
                    "bio_mmcif_switches",
                    "structure_tags",
                    "tautomer_long_conjugated",
                    "tautomer_focused",
                    "molalign_focused"
                ]
                .contains(&name.as_str()),
                "unknown special regression: {name}"
            );
            for key in [
                "bio_mmcif_switches",
                "structure_tags",
                "tautomer_long_conjugated",
                "tautomer_focused",
                "molalign_focused",
            ] {
                println!("cargo:rustc-check-cfg=cfg(parity_special_{key})");
                if name == "all" || name == key {
                    println!("cargo:rustc-cfg=parity_special_{key}");
                }
            }
        }
    }
    let path = "testdata/special/structure_tags.json";
    println!("cargo:rerun-if-changed={path}");
    let fixture: serde_json::Value = serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
    let mut output = String::new();
    let mut names = std::collections::BTreeSet::new();
    for table in ["cases", "octahedral_switch_cases"] {
        for case in fixture[table].as_array().unwrap() {
            let id = case["case_id"].as_str().unwrap();
            let name: String = id
                .chars()
                .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
                .collect();
            assert!(
                names.insert(name.clone()),
                "duplicate generated test name: {name}"
            );
            output.push_str(&format!(
                "#[test]\nfn case_{name}() {{ compare_case({id:?}); }}\n"
            ));
        }
    }
    fs::write(
        PathBuf::from(env::var_os("OUT_DIR").unwrap()).join("structure_cases.rs"),
        output,
    )
    .unwrap();
}
