"""Focused native-call contract checks; corpus comparisons remain Cargo tests."""
import copy
import json
import unittest
from unittest.mock import patch
from rdkit import Chem, rdBase
from reference import PACKAGE
from fingerprints import fingerprint_case
from pattern_profiles import pattern_case
from _generate_pattern_fingerprint_golden import make_atom_counts, make_set_only_bits


class PatternProfilesTests(unittest.TestCase):
    def test_historical_profiles_forward_every_native_option(self):
        pin = json.loads((PACKAGE / "testdata/reference/rdkit.json").read_text())
        self.assertEqual(rdBase.rdkitVersion, pin["version"])
        profiles = json.loads((PACKAGE / "testdata/profiles/fingerprints_historical_0_3.json").read_text())["profiles"]["Pattern"]
        self.assertEqual(len(profiles), 10)
        mol = Chem.MolFromSmiles("CCO")
        native = Chem.PatternFingerprint
        for profile in profiles:
            with self.subTest(profile=profile["name"]):
                params = profile["params"]["Pattern"]
                options = params["source_options"]
                row = {"Fingerprint": {"case": {"id": "one", "smiles": "CCO"},
                       "right": None, "seed": 1, "params": copy.deepcopy(profile["params"])}}
                original = copy.deepcopy(row)
                with patch("pattern_profiles.Chem.PatternFingerprint", wraps=native) as call:
                    actual = pattern_case(row)
                call.assert_called_once()
                kwargs = call.call_args.kwargs
                self.assertEqual(kwargs["fpSize"], params["fp_size"])
                self.assertEqual(kwargs["tautomerFingerprints"], params["tautomeric"])
                counts = make_atom_counts(mol, options["atom_counts"])
                mask = make_set_only_bits(params["fp_size"], options["set_only_bits"])
                self.assertEqual("atomCounts" in kwargs, counts is not None)
                self.assertEqual("setOnlyBits" in kwargs, mask is not None)
                if counts is not None:
                    self.assertEqual(kwargs["atomCounts"], counts)
                if mask is not None:
                    self.assertEqual(kwargs["setOnlyBits"].GetNumBits(), mask.GetNumBits())
                    self.assertEqual(list(kwargs["setOnlyBits"].GetOnBits()), list(mask.GetOnBits()))
                self.assertEqual(row, original)
                self.assertEqual(actual, fingerprint_case(row))

    def test_absent_options_keep_seeded_reference_behavior(self):
        row = {"Fingerprint": {"case": {"id": "one", "smiles": "CCO"},
               "right": None, "seed": 1,
               "params": {"Pattern": {"fp_size": 257, "tautomeric": True}}}}
        self.assertEqual(pattern_case(row), fingerprint_case(row))


if __name__ == "__main__":
    unittest.main()
