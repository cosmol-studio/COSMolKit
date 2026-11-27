//! Source-built references used only by the new descriptor suites.
use crate::Result;
use crate::storage::hex;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use sha2::{Digest, Sha256};
use std::{
    fs::File,
    io::Read,
    path::{Path, PathBuf},
};

pub(crate) const TASKS: &[(&str, &str)] = &[
    ("gasteiger_charges", "descriptor-foundations"),
    ("partial_charge_extrema", "descriptor-foundations"),
    ("peoe_vsa", "descriptor-foundations"),
    ("estate_indices", "descriptor-foundations"),
    ("estate_extrema", "descriptor-foundations"),
    ("estate_vsa", "descriptor-foundations"),
    ("vsa_estate", "descriptor-foundations"),
    ("bcut2d", "descriptor-foundations"),
    ("ipc", "descriptor-foundations-ipc"),
    ("avg_ipc", "descriptor-foundations-ipc"),
];

pub(crate) fn kind(task: &str) -> Option<&'static str> {
    TASKS
        .iter()
        .find(|(name, _)| *name == task)
        .map(|(_, suite)| {
            if *suite == "descriptor-foundations-ipc" {
                "exact-ipc"
            } else {
                "master"
            }
        })
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct SmilesInput {
    smiles: String,
}

pub(crate) fn smiles_inputs(rows: &[String]) -> Result<Vec<String>> {
    rows.iter()
        .map(|line| Ok(serde_json::from_str::<SmilesInput>(line)?.smiles))
        .collect()
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NativeBuild {
    binary_sha256: String,
    build_info: BuildInfo,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct BuildInfo {
    protocol: u32,
    kind: String,
    revision: String,
    adapter_sha256: String,
    cmake_sha256: String,
    compiler: String,
    build_type: String,
    flags: String,
    boost: String,
    eigen: String,
    simd: String,
    rounding: i64,
}

impl NativeBuild {
    pub(crate) fn load(root: &Path, kind: &str) -> Result<Self> {
        let source = root.join(kind).join("identity.json");
        let input = File::open(&source).map_err(|e| {
            format!(
                "{}: {e}; build the {kind} descriptor oracle first",
                source.display()
            )
        })?;
        let build: Self = serde_json::from_reader(input.take(16384))?;
        let pins: Value = serde_json::from_str(include_str!(
            "../../testdata/reference/descriptor_foundations.json"
        ))?;
        let info = &build.build_info;
        if pins[kind].as_str() != Some(info.revision.as_str())
            || info.kind != kind
            || info.eigen != "3.4.0"
            || info.protocol != 1
            || info.adapter_sha256
                != hex(&Sha256::digest(include_bytes!(
                    "../tools/descriptor_reference/reference.cpp"
                )))
            || info.cmake_sha256
                != hex(&Sha256::digest(include_bytes!(
                    "../tools/descriptor_reference/rdkit.cmake"
                )))
            || build.binary_sha256.len() != 64
            || !build.binary_sha256.bytes().all(|b| b.is_ascii_hexdigit())
        {
            return Err(format!(
                "stale or incompatible {kind} descriptor oracle identity; rebuild the oracle"
            )
            .into());
        }
        Ok(build)
    }

    pub(crate) fn revision(&self) -> &str {
        &self.build_info.revision
    }

    pub(crate) fn executable(&self, root: &Path) -> PathBuf {
        root.join(&self.build_info.kind).join(&self.binary_sha256)
    }
}
