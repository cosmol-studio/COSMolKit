"""Check npm package shape without rebuilding or publishing the WASM API."""

import json
from pathlib import Path
import subprocess
import tempfile
import unittest

from run import prepare_package


class NpmPackageTests(unittest.TestCase):
    def test_metadata_and_packed_inline_js_dependencies(self):
        with tempfile.TemporaryDirectory() as temporary:
            package = Path(temporary)
            for name in ("wasm_wasm.js", "wasm_wasm.d.ts", "wasm_wasm_bg.wasm", "wasm_wasm_bg.wasm.d.ts"):
                (package / name).write_bytes(b"fixture")
            snippets = package / "snippets" / "fixture"
            snippets.mkdir(parents=True)
            (snippets / "inline0.js").write_text("export const value = 1;", encoding="utf-8")
            (package / "not-for-publication.log").write_text("fixture", encoding="utf-8")
            prepare_package(package, "wasm_wasm", {"version": "0.5.0-rc.13", "license": "MIT", "repository": "https://github.com/cosmol-studio/COSMolKit"})
            metadata = json.loads((package / "package.json").read_text())
            self.assertEqual(metadata["name"], "@cosmol-studio/cosmolkit")
            self.assertEqual(metadata["version"], "0.5.0-rc.13")
            self.assertEqual(metadata["type"], "module")
            self.assertEqual(metadata["exports"]["."]["default"], "./wasm_wasm.js")
            self.assertEqual(metadata["exports"]["."]["types"], "./wasm_wasm.d.ts")
            self.assertEqual(metadata["exports"]["./*"], "./*")
            packed = subprocess.run(["npm", "pack", "--dry-run", "--json"], cwd=package, text=True, capture_output=True, check=True)
            files = {item["path"] for item in json.loads(packed.stdout)[0]["files"]}
            self.assertEqual(files, {"package.json", "README.md", "LICENSE", "wasm_wasm.js", "wasm_wasm.d.ts", "wasm_wasm_bg.wasm", "wasm_wasm_bg.wasm.d.ts", "snippets/fixture/inline0.js"})


if __name__ == "__main__":
    unittest.main()
