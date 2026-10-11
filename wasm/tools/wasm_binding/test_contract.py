"""Small projection-checker regressions; no WASM build or nested Cargo run."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

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
    def test_compiled_conveniences_keep_canonical_platform_gates(self):
        document = {"entries": [], "javascript_aliases": [
            {"name": "mol_from_binary", "target": "Molecule.fromBinary", "platform": "native"},
            {"name": "mol_from_smiles", "target": "Molecule.fromSmiles", "platform": "all"},
        ]}
        result = subprocess.CompletedProcess([], 0, json.dumps(document))
        with tempfile.TemporaryDirectory() as directory, patch.object(_checker.subprocess, "run", return_value=result) as run:
            projected = _checker.compiled_contract(ROOT, {"core"}, Path(directory) / "contract.json")
        self.assertEqual(projected["javascript_aliases"], document["javascript_aliases"][1:])
        self.assertIn("--no-default-features", run.call_args.args[0])
        self.assertEqual(run.call_args.args[0][-1], "core")

    def test_conveniences_are_direct_aliases_with_inherited_types(self):
        document = {"entries": [], "javascript_aliases": [
            {"name": "bio_from_mmcif", "target": "BioStructure.fromMmcif", "platform": "all"},
        ]}
        with tempfile.TemporaryDirectory() as directory:
            package = Path(directory)
            for suffix in ["js", "d.ts"]:
                (package / f"example.{suffix}").write_text("// generated\n")
            _checker.generate_conveniences(package, "example", document)
            self.assertIn("Alias of `BioStructure.fromMmcif`.", (package / "example.d.ts").read_text())
            self.assertIn("export const bioFromMmcif = BioStructure.fromMmcif;", (package / "example.js").read_text())
            self.assertIn("export const bioFromMmcif: typeof BioStructure.fromMmcif;", (package / "example.d.ts").read_text())
        probe = _checker.typescript_probe(document)
        self.assertIn("const alias0: typeof ck.BioStructure.fromMmcif = ck.bioFromMmcif;", probe)
        self.assertIn("const original0: typeof ck.bioFromMmcif = ck.BioStructure.fromMmcif;", probe)

    def test_local_contract_build_uses_unoptimized_profile(self):
        document = {"entries": [], "javascript_aliases": []}
        result = subprocess.CompletedProcess([], 0, json.dumps(document))
        with tempfile.TemporaryDirectory() as directory, patch.object(_checker.subprocess, "run", return_value=result) as run:
            _checker.compiled_contract(ROOT, {"core"}, Path(directory) / "contract.json", profile="dev")
        command = run.call_args.args[0]
        self.assertEqual(command[command.index("--profile") + 1], "dev")

    def test_convenience_gate_rejects_missing_exports_and_wrappers(self):
        executable = shutil.which("node") or shutil.which("bun")
        self.assertIsNotNone(executable, "Install Node or Bun to check JavaScript bindings")
        checker = (ROOT / "wasm/tools/wasm_binding/contract.mjs").as_uri()
        script = f"""
import {{checkConveniences}} from {json.dumps(checker)};
import assert from 'node:assert/strict';
class BioStructure {{static fromMmcif(text, options) {{return [text, options];}}}}
const document = {{javascript_aliases:[{{name:'bio_from_mmcif',target:'BioStructure.fromMmcif'}}]}};
assert.ok(checkConveniences({{BioStructure}}, document).some(error => error.includes('missing')));
assert.ok(checkConveniences({{BioStructure,bioFromMmcif:(...args)=>BioStructure.fromMmcif(...args)}}, document).some(error => error.includes('directly reference')));
assert.deepEqual(checkConveniences({{BioStructure,bioFromMmcif:BioStructure.fromMmcif}}, document), []);
"""
        subprocess.run([executable, "--input-type=module", "-e", script], check=True, capture_output=True, text=True)

    def test_static_configuration_call_uses_class_not_instance(self):
        document = _document()
        document["entries"] += [
            {"semantic_id": "Molecule.from_smiles", "javascript_name": "fromSmiles", "item": "callable",
             "owner": "molecule", "parameters": [{"name": "text", "type": "&str"}], "receiver": None},
            {"semantic_id": "Molecule.from_smiles_with_params", "javascript_name": "fromSmilesWithParams", "item": "callable",
             "owner": "molecule", "parameters": [{"name": "text", "type": "&str"},
                {"name": "params", "type": "&crate::Configuration"}], "receiver": None},
        ]
        probe = _checker.typescript_probe(document)
        self.assertIn("ck.Molecule.fromSmiles(null as unknown as Parameters<typeof ck.Molecule.fromSmiles>[0], configuration0_0)", probe)
        self.assertNotIn("declare const receiver0", probe)

    def test_numeric_word_casing_matches_generated_javascript(self):
        self.assertEqual(_checker.camel("use_2d"), "use2d")
        self.assertEqual(_checker.camel("clear_existing_2d"), "clearExisting2d")

    def test_configuration_probe_checks_every_configuration_object(self):
        document = _document()
        document["entries"] += [
            {"semantic_id": "Molecule.example", "javascript_name": "example", "item": "callable",
             "owner": "molecule", "parameters": [], "receiver": "Shared"},
            {"semantic_id": "Molecule.example_with_params", "javascript_name": "exampleWithParams", "item": "callable",
             "owner": "molecule", "parameters": [
                 {"name": "first", "type": "&crate::Configuration"},
                 {"name": "second", "type": "&crate::Configuration"}], "receiver": "Shared"},
        ]
        probe = _checker.typescript_probe(document)
        self.assertIn("receiver0.example(configuration0_0, configuration0_1)", probe)
        self.assertIn("maxMatches: configuration0_0.maxMatches", probe)
        self.assertIn("maxMatches: configuration0_1.maxMatches", probe)
    def test_registered_configuration_does_not_require_rust_new(self):
        document = _document()
        row = document["entries"][0]
        row["python_fields"] = row["fields"]
        row["fields"] = None
        _checker.configuration_fields(document)
        self.assertIn("p0.maxMatches = p0.maxMatches;", _checker.typescript_probe(document))
        row["fields"] = []  # A zero-argument Rust factory is not the language schema.
        _checker.configuration_fields(document)
        self.assertIn("p0.maxMatches = p0.maxMatches;", _checker.typescript_probe(document))

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
