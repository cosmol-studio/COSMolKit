use std::path::PathBuf;

pub fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

pub fn profile_name() -> String {
    std::env::var("COSMOLKIT_PARITY_PROFILE").unwrap_or_else(|_| "smiles_small".to_string())
}

pub fn smiles_path() -> PathBuf {
    if let Ok(path) = std::env::var("COSMOLKIT_PARITY_SMILES") {
        return PathBuf::from(path);
    }
    match profile_name().as_str() {
        "small" | "smiles_small" => repo_root().join("tests/smiles.smi"),
        "strict" | "5000" | "smiles_5000" => repo_root().join("tests/smiles_5000.smi"),
        other => panic!(
            "unknown COSMOLKIT_PARITY_PROFILE '{other}'; expected smiles_small or smiles_5000"
        ),
    }
}

pub fn golden_dir() -> PathBuf {
    if let Ok(path) = std::env::var("COSMOLKIT_PARITY_GOLDEN_DIR") {
        return PathBuf::from(path);
    }
    match profile_name().as_str() {
        "small" | "smiles_small" => repo_root().join("tests/golden/smiles_small"),
        "strict" | "5000" | "smiles_5000" => repo_root().join("tests/golden/smiles_5000"),
        other => panic!(
            "unknown COSMOLKIT_PARITY_PROFILE '{other}'; expected smiles_small or smiles_5000"
        ),
    }
}

pub fn golden_path(file_name: &str) -> PathBuf {
    golden_dir().join(file_name)
}

pub fn count_smiles_rows() -> usize {
    let path = smiles_path();
    std::fs::read_to_string(&path)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", path.display()))
        .lines()
        .filter(|line| {
            let line = line.trim();
            !line.is_empty() && !line.starts_with('#')
        })
        .count()
}

pub fn regenerate_command() -> String {
    format!(
        ".venv/bin/python tests/scripts/gen_all_rdkit_goldens.py --python .venv/bin/python --profile {} --suite default --clean --jobs 4",
        profile_name()
    )
}
