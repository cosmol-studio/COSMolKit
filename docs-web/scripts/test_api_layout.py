"""API presentation must preserve generated declarations and member anchors."""

from pathlib import Path
from types import SimpleNamespace
import sys
import unittest

from docutils import nodes
from docutils.utils import new_document
from sphinx import addnodes

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python/docs"))
from api_layout import api_toc, format_api


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
        format_api(self.app, self.tree, "api")
        names = [node[0]["fullname"] for node in self.section
                 if isinstance(node, addnodes.desc)]
        self.assertEqual(names, ["bio_from_mmcif", "mol_from_smiles", "Molecule",
                                 "BioStructure", "Protein", "Atom", "Zebra"])

    def test_closed_members_preserve_content_and_anchors(self):
        format_api(self.app, self.tree, "api")
        raw = [node for node in self.tree.findall(nodes.raw) if not node.get("api_member_anchor")]
        self.assertEqual([node.astext() for node in raw], [
            '<details class="api-members"><summary>Members (1)</summary>', '</details>'])
        self.assertEqual(self.method[0]["ids"], ["cosmolkit.Molecule.to_smiles"])
        self.assertIn("Return a SMILES string.", self.method.astext())
        self.assertEqual(self.method.parent["classes"], ["api-member-list"])
        self.assertIn("A molecule.", self.molecule[1][0].astext())
        format_api(self.app, self.tree, "api")
        self.assertEqual(len(list(self.tree.findall(nodes.raw))), 3)
        self.assertEqual(sum(bool(node.get("api_member_anchor"))
                             for node in self.tree.findall(nodes.raw)), 1)

    def test_flat_toc_matches_body_without_members_or_page_title(self):
        format_api(self.app, self.tree, "api")
        captured = []
        self.app.builder = SimpleNamespace(render_partial=lambda toc:
            captured.append(toc) or {"fragment": "rendered"})
        context = {}
        api_toc(self.app, "api", "page.html", context, self.tree)
        toc = captured[0]
        self.assertEqual(context["toc"], "rendered")
        self.assertEqual(len(list(toc.findall(nodes.bullet_list))), 1)
        links = list(toc.findall(nodes.reference))
        self.assertEqual([node.astext() for node in links], [
            "bio_from_mmcif()", "mol_from_smiles()", "Molecule", "BioStructure",
            "Protein", "Atom", "Zebra"])
        self.assertEqual(links[2]["refuri"], "#cosmolkit.Molecule")

    def test_javascript_uses_same_layout(self):
        format_api(self.app, self.tree, "javascript-api")
        self.assertEqual(self.section[1][0]["fullname"], "bio_from_mmcif")
        self.assertEqual(len(list(self.tree.findall(nodes.raw))), 3)

    def test_other_pages_unchanged(self):
        before = self.tree.pformat()
        format_api(self.app, self.tree, "quickstart")
        self.assertEqual(self.tree.pformat(), before)


if __name__ == "__main__":
    unittest.main()
