use std::fs;
use std::path::{Path, PathBuf};

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn gather_rs_files(dir: &Path, out: &mut Vec<PathBuf>) {
    let entries = fs::read_dir(dir).expect("read_dir should succeed");
    for entry in entries {
        let entry = entry.expect("dir entry should be readable");
        let path = entry.path();
        if path.is_dir() {
            gather_rs_files(&path, out);
        } else if path.extension().is_some_and(|ext| ext == "rs") {
            out.push(path);
        }
    }
}

fn strip_after_cfg_test(content: &str) -> &str {
    // Keep runtime path strict and avoid false positives from test-only helpers.
    if let Some(pos) = content.find("#[cfg(test)]") {
        &content[..pos]
    } else {
        content
    }
}

fn assert_no_banned_patterns(path: &Path, content: &str) {
    let banned = ["py.import(\"rdkit", "py.import('rdkit", "Command::new("];
    for needle in banned {
        assert!(
            !content.contains(needle),
            "banned runtime dependency pattern `{}` found in {}",
            needle,
            path.display()
        );
    }
}

fn is_out_of_line_test_module(path: &Path) -> bool {
    path.file_name().is_some_and(|name| name == "tests.rs")
}

#[test]
fn rust_and_python_binding_runtime_code_must_not_call_rdkit_or_shell_out() {
    let root = repo_root();

    let mut rust_src_files = Vec::new();
    gather_rs_files(&root.join("crates/cosmolkit-core/src"), &mut rust_src_files);
    gather_rs_files(&root.join("crates/cosmolkit/src"), &mut rust_src_files);
    for path in rust_src_files {
        if is_out_of_line_test_module(&path) {
            continue;
        }
        let content = fs::read_to_string(&path).expect("source file should be readable");
        assert_no_banned_patterns(&path, strip_after_cfg_test(&content));
    }

    let mut py_binding_files = Vec::new();
    gather_rs_files(&root.join("python/src"), &mut py_binding_files);
    for path in py_binding_files {
        let content = fs::read_to_string(&path).expect("python binding file should be readable");
        assert_no_banned_patterns(&path, &content);
    }
}
