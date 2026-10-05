//! Preparation-only adapter for the registered pinned-Gemmi coordinate tasks.
use crate::{Corpus, Record, Result, Task, digest, encode, read, registry, root};
use std::{
    io::Write,
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

const ADAPTER: &str = include_str!("../../tools/oracles/gemmi/bio_pdb_values.py");

pub(crate) fn handles(task: &Task) -> bool {
    task.operation == registry::Operation::BioPdbOutput
}

fn native_path() -> PathBuf {
    std::env::var_os("BIO_PDB_NATIVE_ORACLE")
        .map(PathBuf::from)
        .unwrap_or_else(|| root().join("target/gemmi-mmcif-writer-golden/pdb_coordinate_oracle_v2"))
}

pub(crate) fn provenance(task: &Task) -> Option<serde_json::Value> {
    if !handles(task) {
        return None;
    }
    let binary = read(&native_path()).ok()?;
    Some(serde_json::json!({
        "library": registry::BIO_PDB_REFERENCE.library,
        "version": registry::BIO_PDB_REFERENCE.version,
        "source_revision": registry::BIO_PDB_REFERENCE.commit,
        "native_binary_sha256": digest(&binary),
        "projection": "ATOM/HETATM/ANISOU/TER/MODEL/ENDMDL/END, original bytes retained"
    }))
}

pub(crate) fn generate(task: &Task, cases: &Corpus, python: &Path) -> Result<Vec<Record>> {
    let script = root().join("tools/oracles/gemmi/bio_pdb_values.py");
    if read(&script)? != ADAPTER.as_bytes() {
        return Err("Gemmi adapter source changed; rebuild the preparation binary".into());
    }
    let binary = native_path();
    // Read-only preparation identity check; comparison never runs this binary.
    read(&binary)?;
    let inputs = registry::expand(cases, task);
    let mut child = Command::new(python)
        .arg(script)
        .arg(&binary)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("Gemmi adapter: {e}"))?;
    let mut stdin = child
        .stdin
        .take()
        .ok_or("Gemmi adapter stdin unavailable")?;
    let payload = encode(&inputs)?;
    let writer = std::thread::spawn(move || stdin.write_all(&payload));
    let output = child.wait_with_output().map_err(|e| e.to_string())?;
    let sent = writer
        .join()
        .map_err(|_| "Gemmi adapter input writer panicked")?;
    if !output.status.success() {
        return Err(format!(
            "Gemmi adapter failed: {}",
            String::from_utf8_lossy(&output.stderr)
        ));
    }
    sent.map_err(|e| e.to_string())?;
    serde_json::from_slice(&output.stdout).map_err(|e| format!("Gemmi adapter output: {e}"))
}
