"""Pure option-transport regressions, without preparing chemistry references."""
import copy
import itertools
import unittest

from reference import writer_native_profile


class WriterProfileTests(unittest.TestCase):
    def test_all_768_profiles_preserve_values_order_and_input(self):
        fields = ["isomeric_smiles", "kekule", "canonical", "clean_stereo",
                  "all_bonds_explicit", "all_hs_explicit", "include_dative_bonds",
                  "ignore_atom_map_numbers"]
        observed = []
        expected = []
        for values in itertools.product([False, True], repeat=8):
            for root in ["none", "first", "last"]:
                profile = dict(zip(fields, values, strict=True), rooted_at_atom=root)
                before = copy.deepcopy(profile)
                native = writer_native_profile(profile)
                self.assertEqual(profile, before)
                self.assertEqual(len(native), 9)
                self.assertNotIn("isomeric_smiles", native)
                self.assertNotIn("kekule", native)
                native_fields = ["do_isomeric_smiles", "do_kekule", *fields[2:]]
                observed.append(tuple(native[field] for field in native_fields)
                                + (native["rooted_at_atom"],))
                expected.append(values + (None if root == "none" else root,))
        self.assertEqual(len(observed), 768)
        self.assertEqual(len(set(observed)), 768)
        self.assertEqual(observed, expected)

    def test_missing_current_options_fail_instead_of_defaulting(self):
        profile = {"isomeric_smiles": True, "kekule": False, "rooted_at_atom": "none"}
        for field in profile:
            incomplete = dict(profile)
            del incomplete[field]
            with self.subTest(field=field), self.assertRaises(KeyError):
                writer_native_profile(incomplete)


if __name__ == "__main__":
    unittest.main()
