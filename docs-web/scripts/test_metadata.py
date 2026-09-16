from pathlib import Path
import tempfile
import unittest

from extract_sphinx_metadata import generate_metadata
from html_metadata import PageMetadata
from route_contract import PAGES


class MetadataTests(unittest.TestCase):
    def test_only_document_head_metadata_counts(self):
        page = PageMetadata()
        page.feed('<head><title>Document</title><meta name="description" content="Atoms &amp; bonds"></head><body><svg><title>Icon</title></svg><meta name="robots" content="noindex"></body>')
        self.assertEqual(page.titles, ["Document"])
        self.assertEqual(page.descriptions, ["Atoms & bonds"])
        self.assertEqual(page.robots, set())

    def test_sphinx_metadata_preserves_description_without_exporting_index_policy(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path / "guide.html").write_text('<head><meta name="description" content="Atoms &amp; bonds: &quot;#分子"><meta name="robots" content="noindex, follow"></head>', encoding="utf-8")
            generated = generate_metadata(path, ["guide"])
            self.assertIn('=> r##"Atoms & bonds: "#分子"##,', generated)
            self.assertNotIn('bool', generated)
            self.assertNotIn('true', generated)

    def test_missing_guide_metadata_is_not_replaced_with_generic_text(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path / "guide.html").write_text('<head></head>', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "description"):
                generate_metadata(path, ["guide"])

    def test_generated_utility_descriptions_and_contract_noindex(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            for page in ("search", "genindex", "py-modindex"):
                (path / f"{page}.html").write_text('<head></head>', encoding="utf-8")
                generated = generate_metadata(path, [page])
                self.assertIn("COSMolKit", generated)
                route = next(route for route in PAGES if route.get('docname') == page)
                self.assertFalse(route['indexable'])


if __name__ == "__main__":
    unittest.main()
