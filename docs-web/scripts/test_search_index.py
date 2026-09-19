"""Search data must come from declared pages and preserve symbol destinations."""

import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from generate_search_index import generate_index


class SearchIndexTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.output = self.root / "records.json"
        self.index = {"docnames": ["api"], "titles": ["API Reference"],
                      "objects": {"cosmolkit": [[0, 0, 0, "", "Molecule"]]},
                      "objnames": {"0": ["py", "class", "Python class"]},
                      "alltitles": {"Molecules": [[0, "molecules"]]}}
        (self.root / "api.html").write_text('<head><meta name="description" content="Molecules &amp; atoms"></head><nav>Not article text</nav><article id="furo-main-content">Body text</article>', encoding="utf-8")
        self.write_index()

    def write_index(self):
        (self.root / "searchindex.js").write_text("Search.setIndex(" + json.dumps(self.index) + ")", encoding="utf-8")

    def test_page_body_and_symbol_anchor_are_preserved(self):
        self.assertEqual(generate_index(self.root, self.output), 3)
        records = json.loads(self.output.read_text(encoding="utf-8"))
        self.assertEqual(records[0]["url"], "/python/api")
        self.assertEqual(records[0]["text"], "Body text")
        self.assertEqual(records[0]["summary"], "Molecules & atoms")
        self.assertEqual(records[1]["url"], "/python/api#cosmolkit.Molecule")

    def test_undeclared_source_is_rejected(self):
        self.index["docnames"] = ["unknown"]
        self.write_index()
        with self.assertRaisesRegex(ValueError, "absent from routes.toml"):
            generate_index(self.root, self.output)

    def test_missing_source_and_invalid_index_fail(self):
        (self.root / "api.html").unlink()
        with self.assertRaises(FileNotFoundError):
            generate_index(self.root, self.output)
        (self.root / "searchindex.js").write_text("changed format", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "index format"):
            generate_index(self.root, self.output)

    def test_binding_tool_version_must_match_the_rust_dependency(self):
        from build_search_bundle import binding_tool
        (self.root / "Cargo.toml").write_text('[dependencies]\nwasm-bindgen = {version = "=0.2.128"}\n')
        with patch("build_search_bundle.shutil.which", return_value=None), patch.dict("os.environ", {"COSMOLKIT_WASM_BINDGEN": "missing"}):
            with self.assertRaisesRegex(RuntimeError, "0.2.128 is required"):
                binding_tool(self.root)

    def test_generated_assets_use_crate_relative_paths(self):
        from build_search_bundle import build_search_bundle
        from pathlib import PurePosixPath
        out = self.root / 'target/wasm32-unknown-unknown/wasm-release/build/docs/out'
        out.mkdir(parents=True)
        # Stub external compilation, but create the actual binding output files.
        def compile_assets(command, **kwargs):
            if command[0] == 'wasm-bindgen':
                (out / 'docs_search.js').write_text('// generated bindings')
                (out / 'docs_search_bg.wasm').write_bytes(b'\x00asm\x01\x00\x00\x00')
        script = str(self.root / 'scripts/build_search_bundle.py')
        if os.name == 'nt':
            script = '\\\\?\\' + script
        with patch('build_search_bundle.__file__', script), \
             patch('build_search_bundle.binding_tool', return_value='wasm-bindgen'), \
             patch('build_search_bundle.subprocess.run', side_effect=compile_assets):
            build_search_bundle(out)
        generated = (out / 'search_asset.rs').read_text(encoding='utf-8')
        for name in ('docs_search.js', 'docs_search_bg.wasm'):
            relative = (out / name).relative_to(self.root).as_posix()
            self.assertIn(f'asset!(r#"/{relative}"#)', generated)
            # Match Manganis's Unix resolution: trim '/' and join crate root.
            linux_root = PurePosixPath('/home/runner/work/COSMolKit/COSMolKit/docs-web')
            resolved = linux_root / ('/' + relative).lstrip('/')
            self.assertEqual(resolved.relative_to(linux_root).as_posix(), relative)
            self.assertTrue((self.root / relative).is_file())
