import json
import unittest
from pathlib import Path

import special

INPUT = Path(__file__).resolve().parents[1] / "testdata/special"


class SpecialTests(unittest.TestCase):
    def test_structure_keeps_all_77_actual_constructions(self):
        rows = [json.loads(line) for line in (INPUT / "structure_tags.jsonl").read_text().splitlines()]
        self.assertEqual(len(rows), 77)
        for row in rows:
            before = special.Chem.GetAllowNontetrahedralChirality()
            result = special.structure_tags(row)
            self.assertEqual(len(result), 3)
            self.assertEqual(len(result[1][0]), len(row["atoms"]))
            self.assertEqual(special.Chem.GetAllowNontetrahedralChirality(), before)

    def test_optimizer_keeps_native_molecule_owner_alive(self):
        row = json.loads((INPUT / "forcefield_optimizers.jsonl").read_text().splitlines()[0])
        for _ in range(3):
            result = special.optimizer(row)
            self.assertTrue(result[0])
            self.assertEqual(len(result[1][1]), 3 * len(result[1][2][0][2]))
            self.assertEqual(len(result[4][0]), 2)

    def test_alignment_retains_all_14_calls(self):
        rows = [json.loads(line) for line in (INPUT / "molalign_focused.jsonl").read_text().splitlines()]
        self.assertEqual(sum(len(row["case"]["calls"]) for row in rows), 14)
        for row in rows:
            self.assertEqual(len(special.alignment(row["case"])), len(row["case"]["calls"]))


if __name__ == "__main__":
    unittest.main()
