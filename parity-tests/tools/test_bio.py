import itertools
import unittest
from pathlib import Path
import bio


class BioTests(unittest.TestCase):
    def test_complete_inputs_and_pdb_switches(self):
        rows = (Path(__file__).parents[1] / "testdata/bio/small.jsonl").read_text().splitlines()
        names = ("ter_records", "numbered_ter", "ter_ignores_type", "preserve_serial", "end_record")
        count = 0
        for row in rows:
            structure = bio.parse(row)
            for flags in itertools.product((False, True), repeat=5):
                profile = dict(zip(names, flags), gemmi_version=bio.GEMMI)
                result = bio.pdb(structure, profile)
                self.assertIn(b"ATOM", result)
                self.assertEqual(result.endswith(b"END                                                                             \n"), profile["end_record"])
                count += 1
        self.assertEqual(count, 64)


if __name__ == "__main__":
    unittest.main()
