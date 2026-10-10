"""Native rejection discriminants retain their source exception payloads."""
import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools/oracles/rdkit"))
from fingerprint_values_pilot import molecular


class RejectionTests(unittest.TestCase):
    def test_parse_rejection_does_not_produce_a_descriptor(self):
        value = molecular({"case": {"smiles": "C==C"}, "profile": "Kappa1"})
        self.assertEqual(set(value), {"ParseRejected"})

    def test_sanitize_exception_payloads(self):
        cases = [
            ("C(C)(C)(C)(C)C", {"ExplicitValence": {"atom": 0, "message": "Explicit valence for atom # 0 C, 5, is greater than permitted"}}),
            ("c1cccc1", {"Kekulize": {"atoms": [0, 1, 2, 3, 4]}}),
            ("c1c(ccc2NC(CN=c(c21)(C)C)=O)O", {"ExplicitValence": {"atom": 9, "message": "Explicit valence for atom # 9 C, 5, is greater than permitted"}}),
        ]
        for smiles, expected in cases:
            value = molecular({"case": {"smiles": smiles}, "profile": "SanitizeAll"})
            self.assertEqual(value["SanitizeRejected"]["reason"], expected)


if __name__ == "__main__":
    unittest.main()
