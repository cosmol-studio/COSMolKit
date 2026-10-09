"""Small projection-checker regressions; no WASM build or nested Cargo run."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[3]
_spec = importlib.util.spec_from_file_location("_wasm_contract", ROOT / "wasm/tools/wasm_binding/contract.py")
assert _spec is not None and _spec.loader is not None
_checker = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_checker)


def _document():
    return {"entries": [{
        "semantic_id": "types.Configuration", "javascript_name": "Configuration",
        "python_name": "Configuration", "rust_path": "crate::Configuration",
        "owner": "type", "item": "type", "role": "parameter",
        "fields": [{"name": "max_matches", "type": "usize", "default": "100"}],
        "properties": [],
    }]}


class ContractTests(unittest.TestCase):
    def test_typescript_probe_requires_actual_writable_field_and_export(self):
        probe = _checker.typescript_probe(_document())
        self.assertIn("void ck.Configuration;", probe)
        self.assertIn("p0.maxMatches = p0.maxMatches;", probe)
        invalid = _document()
        invalid["entries"][0]["fields"] = None
        with self.assertRaisesRegex(ValueError, "registered canonical constructor"):
            _checker.typescript_probe(invalid)


    def test_actual_javascript_export_and_descriptor_gate(self):
        executable = shutil.which("node") or shutil.which("bun")
        self.assertIsNotNone(executable, "Install Node or Bun to check JavaScript bindings")
        checker = (ROOT / "wasm/tools/wasm_binding/contract.mjs").as_uri()
        script = f"""
import {{checkContract}} from {json.dumps(checker)};
import assert from 'node:assert/strict';
const document = {json.dumps(_document())};
assert.ok(checkContract({{}}, document).some(error => error.includes('export missing')));
class Configuration {{ get maxMatches() {{return 100;}} }}
assert.ok(checkContract({{Configuration}}, document).some(error => error.includes('getter/setter missing')));
Object.defineProperty(Configuration.prototype, 'maxMatches', {{get(){{return 100;}},set(value){{}}}});
assert.deepEqual(checkContract({{Configuration}}, document), []);
"""
        subprocess.run([executable, "--input-type=module", "-e", script], check=True, capture_output=True, text=True)
