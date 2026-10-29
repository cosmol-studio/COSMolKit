"""Generated native API documentation must retain the checked stub's types."""

import ast
import importlib.util
import io
import json
import pickle
from pathlib import Path
import sys
import subprocess
import tempfile
from types import SimpleNamespace
import unittest

from sphinx.application import Sphinx


ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"
sys.path.insert(0, str(DOCS))
from api_types import public_description, read_stub, return_signature, type_text


class DocumentationWorkflowTests(unittest.TestCase):
    def test_javascript_optional_parameters_survive_rendering(self):
        # One protocol fixture, not a business API inventory. Required union
        # types must stay required; defaults and rest parameters stay intact.
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "node_modules").symlink_to(DOCS / "node_modules", target_is_directory=True)
            (root / "sample.ts").write_text('''
export class Example {
    constructor(options?: { enabled: boolean } | null) {}
    choose(required: string | undefined, optional?: number): void {}
    defaults(count = 3, ...names: string[]): void {}
}
export function query(required: string | null, optional?: boolean): number { return 0; }
export const alias: typeof query = query;
export class FactoryOnly {
    private constructor() {}
    static create(): FactoryOnly { return new FactoryOnly(); }
}
''')
            hook = (DOCS / "sphinx-js.config.ts").read_text().replace(
                'new URL("../target/docs-js-package/enums.rst", import.meta.url)',
                f'new URL({json.dumps((root / "enums.rst").as_uri())})',
            )
            (root / "sphinx-js.config.ts").write_text(hook)
            (root / "tsconfig.json").write_text(json.dumps({
                "compilerOptions": {"target": "ES2022", "strict": True},
                "include": ["sample.ts"],
            }))
            (root / "conf.py").write_text(
                f"import sys\nsys.path.insert(0, {str(DOCS)!r})\n"
                "extensions = ['sphinx_js', 'api_layout']\n"
                "js_language = 'typescript'\njs_source_path = 'sample.ts'\n"
                "root_for_relative_js_paths = '.'\n"
                "jsdoc_tsconfig_path = 'tsconfig.json'\n"
                f"ts_sphinx_js_config = {str(root / 'sphinx-js.config.ts')!r}\n"
                "exclude_patterns = ['enums.rst', 'html', 'doctrees']\n"
            )
            (root / "index.rst").write_text(
                "API\n===\n\n.. js:autoclass:: Example\n   :members:\n\n"
                ".. js:autofunction:: query\n\n.. js:autofunction:: alias\n"
                "\n.. js:autoclass:: FactoryOnly\n"
            )
            result = subprocess.run([
                sys.executable, "-m", "sphinx", "-W", "--keep-going", "-E",
                "-b", "html", "-d", str(root / "doctrees"), str(root), str(root / "html"),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            from sphinx import addnodes
            tree = pickle.loads((root / "doctrees/index.doctree").read_bytes())
            signatures = [n.astext() for n in tree.findall(addnodes.desc_signature)]
            self.assertTrue(any("options?:" in s and "null" in s for s in signatures), signatures)
            for name in ("choose", "query", "alias"):
                signature = next(s for s in signatures if f"{name}(" in s)
                self.assertIn("optional?:", signature)
                self.assertIn("required:", signature)
                self.assertNotIn("required?", signature)
            default = next(s for s in signatures if "defaults(" in s)
            self.assertIn("count=3", default)
            self.assertNotIn("count?", default)
            self.assertIn("...names", default)
            factory = next(s for s in signatures if s.endswith("FactoryOnly"))
            self.assertNotIn("(", factory)

    def test_projection_metadata_preserves_help(self):
        spec = importlib.util.spec_from_file_location("binding_generate", ROOT / "python/generate.py")
        generate = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(generate)
        tree = ast.parse('''
class Sample:
    def values(self) -> object:
        \"\"\"{"role":"return_type","type":"list[int]"}

        Return an independent list of selected values.
        \"\"\"
    def select(self, value):
        \"\"\"Select one value without changing the source.\"\"\"
    def _select_int(self, value: int) -> int:
        \"\"\"{"role":"overload","target":"Sample.select"}\"\"\"
    def _select_text(self, value: str) -> str:
        \"\"\"{"role":"overload","target":"Sample.select"}

        Select by its exact text label.
        \"\"\"
''')
        generate.declared_return_types(tree)
        generate.declared_overloads(tree)
        methods = tree.body[0].body
        self.assertEqual(ast.unparse(methods[0].returns), "list[int]")
        self.assertEqual(ast.get_docstring(methods[0]), "Return an independent list of selected values.")
        self.assertEqual(ast.get_docstring(methods[1]), "Select one value without changing the source.")
        self.assertEqual(ast.get_docstring(methods[2]), "Select by its exact text label.")
        lines = ['{"role":"return_type","type":"list[int]"}', "", "Return selected values."]
        public_description(None, "method", "cosmolkit.Sample.values", None, None, lines)
        self.assertEqual(lines, ["Return selected values."])
        prose = ['{"example": 2}', "This is a documented data example."]
        public_description(None, "method", "cosmolkit.Sample.example", None, None, prose)
        self.assertEqual(prose, ['{"example": 2}', "This is a documented data example."])

    def test_generates_stub_before_documentation_consumers(self):
        workflow = (ROOT / ".github/workflows/docs-web.yml").read_text(encoding="utf-8")
        generate = workflow.index(
            "python python/build_stubs.py"
        )
        for consumer in (
            "unittest discover -s dev/tools -p 'test_python_docs.py'",
            "python -m sphinx -W",
        ):
            self.assertLess(generate, workflow.index(consumer))


class PythonApiTypesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.entries = read_stub(ROOT / "python/artifacts/smoke/package/cosmolkit/__init__.pyi")

    def test_factory_alias_reuses_method_types(self):
        self.assertEqual(
            {type_text(n.returns) for n in self.entries["mol_from_smiles"]},
            {type_text(n.returns) for n in self.entries["Molecule.from_smiles"]},
        )
        self.assertEqual({type_text(n.returns) for n in self.entries["mol_from_smiles"]}, {"Molecule"})

    def test_property_getter_is_not_replaced_by_setter(self):
        entries = self.entries["AddHsParams.explicit_only"]
        self.assertEqual(len(entries), 1)
        self.assertEqual(type_text(entries[0].returns), "bool")

    def test_existing_arguments_and_optional_results_are_preserved(self):
        app = SimpleNamespace(_python_api_types=self.entries)
        result = return_signature(app, "method", "cosmolkit.Molecule.atom", None, None, "(atom_id)", None)
        self.assertEqual(result, ("(atom_id)", "MoleculeAtom | None"))
        self.assertIsNone(return_signature(app, "method", "other.Molecule.atom", None, None, "()", None))

    def test_overload_arguments_remain_paired_with_results(self):
        app = SimpleNamespace(_python_api_types=self.entries)
        signature, result = return_signature(app, "method", "cosmolkit.MoleculeBatch.__getitem__", None, None, "(key)", None)
        self.assertIsNone(result)
        self.assertIn("index: typing.SupportsIndex", signature)
        self.assertIn("-> Molecule | None", signature)
        self.assertIn("index: slice | _typeshed.SupportsGetItem[int, typing.SupportsIndex] | _typeshed.SupportsGetItem[int, bool]", signature)
        self.assertIn("-> MoleculeBatch", signature)
        self.assertNotIn("self", signature)

    def test_sphinx_renders_native_function_method_and_property_types(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "conf.py").write_text(
                f"import sys\nsys.path.insert(0, {str(DOCS)!r})\n"
                "extensions = ['sphinx.ext.autodoc', 'api_types']\n"
                f"python_api_stub = {str(ROOT / 'python/artifacts/smoke/package/cosmolkit/__init__.pyi')!r}\n",
                encoding="utf-8",
            )
            (root / "index.rst").write_text(
                "API\n===\n\n.. autofunction:: cosmolkit.mol_from_smiles\n\n"
                ".. automethod:: cosmolkit.Molecule.to_smiles\n\n"
                ".. automethod:: cosmolkit.Molecule.atom\n\n"
                ".. autoattribute:: cosmolkit.AddHsParams.explicit_only\n\n"
                ".. automethod:: cosmolkit.MoleculeBatch.__getitem__\n\n"
                ".. autoproperty:: cosmolkit.PropertyValueError.expected\n",
                encoding="utf-8",
            )
            warnings = io.StringIO()
            app = Sphinx(str(root), str(root), str(root / "html"), str(root / "doctrees"),
                         "html", status=io.StringIO(), warning=warnings, freshenv=True,
                         warningiserror=True)
            app.build(force_all=True)
            self.assertEqual(app.statuscode, 0, warnings.getvalue())
            from sphinx import addnodes
            signatures = list(app.env.get_doctree("index").findall(addnodes.desc_signature))
            self.assertEqual(len(signatures), 7)
            self.assertIn("-> Molecule", signatures[0].astext())
            self.assertIn("-> str", signatures[1].astext())
            self.assertIn("Atom | None", signatures[2].astext())
            self.assertIn("bool", signatures[3].astext())
            self.assertIn("index: SupportsIndex", signatures[4].astext())
            self.assertIn("-> Molecule | None", signatures[4].astext())
            self.assertIn("index: slice", signatures[5].astext())
            self.assertIn("-> MoleculeBatch", signatures[5].astext())
            self.assertIn("PropertyValueKind", signatures[6].astext())


if __name__ == "__main__":
    unittest.main()
