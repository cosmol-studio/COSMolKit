"""API presentation must preserve generated declarations and member anchors."""

from pathlib import Path
from types import SimpleNamespace
import sys
import unittest

from docutils import nodes
from docutils.utils import new_document
from sphinx import addnodes

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "docs"))
from api_layout import api_toc, format_api, format_javascript_signature


def declaration(kind, name, *children):
    entry = addnodes.desc(objtype=kind)
    signature = addnodes.desc_signature(fullname=name, ids=[f"cosmolkit.{name}"])
    signature += addnodes.desc_name("", name)
    entry += signature
    entry += addnodes.desc_content("", *children)
    return entry


class ApiLayoutTests(unittest.TestCase):
    def setUp(self):
        self.app = SimpleNamespace(config=SimpleNamespace(
            api_class_order=["Molecule", "BioStructure", "Protein"]))
        self.tree = new_document("api")
        self.section = nodes.section(ids=["api-reference"])
        self.section += nodes.title("", "API Reference")
        self.tree += self.section
        self.method = declaration("method", "Molecule.to_smiles",
                                  nodes.paragraph("", "Return a SMILES string."))
        self.molecule = declaration("class", "Molecule",
                                    nodes.paragraph("", "A molecule."), self.method)
        self.section.extend([
            declaration("class", "Zebra"), declaration("class", "Protein"),
            self.molecule, declaration("class", "Atom"),
            declaration("function", "mol_from_smiles"),
            declaration("class", "BioStructure"),
            declaration("function", "bio_from_mmcif"),
        ])

    def test_functions_then_prioritized_classes(self):
        format_api(self.app, self.tree, "python/api")
        names = [node[0]["fullname"] for node in self.section
                 if isinstance(node, addnodes.desc)]
        self.assertEqual(names, ["bio_from_mmcif", "mol_from_smiles", "Molecule",
                                 "BioStructure", "Protein", "Atom", "Zebra"])

    def test_closed_members_preserve_content_and_anchors(self):
        format_api(self.app, self.tree, "python/api")
        raw = [node for node in self.tree.findall(nodes.raw) if not node.get("api_member_anchor")]
        self.assertEqual([node.astext() for node in raw], [
            '<details class="api-members"><summary>Members (1)</summary>', '</details>'])
        self.assertEqual(self.method[0]["ids"], ["cosmolkit.Molecule.to_smiles"])
        self.assertIn("Return a SMILES string.", self.method.astext())
        self.assertEqual(self.method.parent["classes"], ["api-member-list"])
        self.assertIn("A molecule.", self.molecule[1][0].astext())
        format_api(self.app, self.tree, "python/api")
        self.assertEqual(len(list(self.tree.findall(nodes.raw))), 3)
        self.assertEqual(sum(bool(node.get("api_member_anchor"))
                             for node in self.tree.findall(nodes.raw)), 1)

    def test_flat_toc_matches_body_without_members_or_page_title(self):
        format_api(self.app, self.tree, "python/api")
        captured = []
        self.app.builder = SimpleNamespace(render_partial=lambda toc:
            captured.append(toc) or {"fragment": "rendered"})
        context = {}
        api_toc(self.app, "python/api", "page.html", context, self.tree)
        toc = captured[0]
        self.assertEqual(context["toc"], "rendered")
        self.assertEqual(len(list(toc.findall(nodes.bullet_list))), 1)
        links = list(toc.findall(nodes.reference))
        self.assertEqual([node.astext() for node in links], [
            "bio_from_mmcif()", "mol_from_smiles()", "Molecule", "BioStructure",
            "Protein", "Atom", "Zebra"])
        self.assertEqual(links[2]["refuri"], "#cosmolkit.Molecule")

    def test_javascript_uses_same_layout(self):
        format_api(self.app, self.tree, "javascript/api")
        self.assertEqual(self.section[1][0]["fullname"], "bio_from_mmcif")
        self.assertEqual(len(list(self.tree.findall(nodes.raw))), 3)

    def test_other_pages_unchanged(self):
        before = self.tree.pformat()
        format_api(self.app, self.tree, "quickstart")
        self.assertEqual(self.tree.pformat(), before)

    def javascript_callable(self, with_prose=False):
        entry = declaration("function", "compileQuery")
        signature, content = entry
        parameters = addnodes.desc_parameterlist()
        parameters += addnodes.desc_parameter("", "query")
        signature += parameters
        fields = nodes.field_list()
        for name, children in [
            ("param query", [nodes.Text("Query to compile.")] if with_prose else []),
            ("type query", [nodes.inline("", "QueryGraph", classes=["sphinx_js-type"])]),
            ("returns", [nodes.inline("", "CompiledQuery", classes=["sphinx_js-type"])] +
             ([nodes.Text(" -- Reusable compiled query.")] if with_prose else [])),
            ("throws Error", [nodes.Text("If the query is invalid.")]),
        ]:
            fields += nodes.field("", nodes.field_name("", name),
                                  nodes.field_body("", nodes.paragraph("", "", *children)))
        content += fields
        return entry

    def test_javascript_types_move_to_signature_without_empty_fields(self):
        entry = self.javascript_callable()
        format_javascript_signature(self.app, "js", "function", entry[1])
        self.assertEqual(entry[0].astext(), "compileQuery(query: QueryGraph) -> CompiledQuery")
        self.assertEqual([field[0].astext() for field in entry[1][0]], ["throws Error"])
        self.assertIn("If the query is invalid.", entry[1].astext())

    def test_javascript_parameter_and_return_prose_are_retained(self):
        entry = self.javascript_callable(with_prose=True)
        format_javascript_signature(self.app, "js", "function", entry[1])
        self.assertEqual([field[0].astext() for field in entry[1][0]],
                         ["param query", "returns", "throws Error"])
        self.assertIn("Query to compile.", entry[1].astext())
        self.assertIn("Reusable compiled query.", entry[1].astext())
        self.assertNotIn(" -- ", entry[1].astext())
        # No duplicate type/return insertion if another extension repeats it.
        before = entry.pformat()
        format_javascript_signature(self.app, "js", "function", entry[1])
        self.assertEqual(entry.pformat(), before)

    def test_javascript_type_projection_does_not_change_python(self):
        entry = self.javascript_callable()
        before = entry.pformat()
        format_javascript_signature(self.app, "py", "function", entry[1])
        self.assertEqual(entry.pformat(), before)

    def test_javascript_public_namespace_preserves_link_targets(self):
        for kind, name in (("function", "bioFromMmcif"), ("class", "Molecule"),
                           ("attribute", "SdfCoordinateMode")):
            with self.subTest(kind=kind):
                entry = declaration(kind, name)
                signature = entry[0]
                signature["module"] = "wasm_wasm.d"
                signature["ids"] = [f"wasm_wasm.d.{name}"]
                signature.insert(0, addnodes.desc_addname("", "wasm_wasm.d."))
                format_javascript_signature(self.app, "js", kind, entry[1])
                self.assertEqual(signature.astext(), f"cosmolkit.{name}")
                self.assertEqual(signature["ids"], [f"wasm_wasm.d.{name}"])
                self.assertEqual(signature["module"], "wasm_wasm.d")

    def test_javascript_class_member_prefix_is_preserved(self):
        entry = declaration("function", "Molecule.toSmiles")
        entry[0]["module"] = "wasm_wasm.d"
        entry[0].insert(0, addnodes.desc_addname("", "Molecule."))
        format_javascript_signature(self.app, "js", "function", entry[1])
        self.assertEqual(entry[0][0].astext(), "Molecule.")

    def test_javascript_void_result_removes_empty_field_list(self):
        entry = declaration("function", "free")
        entry[0] += addnodes.desc_parameterlist()
        entry[1] += nodes.field_list("", nodes.field(
            "", nodes.field_name("", "returns"),
            nodes.field_body("", nodes.paragraph("", "", nodes.strong("", "void")))))
        format_javascript_signature(self.app, "js", "function", entry[1])
        self.assertEqual(entry[0].astext(), "free() -> void")
        self.assertEqual(len(entry[1]), 0)


if __name__ == "__main__":
    unittest.main()
