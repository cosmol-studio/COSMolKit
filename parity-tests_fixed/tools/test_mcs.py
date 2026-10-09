"""Small preparation-adapter regressions, not another chemistry corpus."""
import unittest

import mcs
import reference


class McsPreparationTests(unittest.TestCase):
    def test_complete_parameter_projection(self):
        params = mcs.parameters({
            "timeout": 30, "store_all": True, "maximize_bonds": False,
            "threshold": 0.75, "initial_seed": "[#6]", "verbose": False,
            "atom_comparator": "isotopes", "bond_comparator": "order_exact",
            "atom_compare_parameters": {key: (3.0 if key == "max_distance" else True)
                                        for key in mcs.ATOM_FIELDS},
            "bond_compare_parameters": {key: True for key in mcs.BOND_FIELDS},
        })
        self.assertEqual(params.Timeout, 30)
        self.assertTrue(params.StoreAll)
        self.assertFalse(params.MaximizeBonds)
        self.assertEqual(params.Threshold, 0.75)
        self.assertEqual(params.InitialSeed, "[#6]")
        self.assertEqual(params.AtomTyper, mcs.ATOMS["isotopes"])
        self.assertEqual(params.BondTyper, mcs.BONDS["order_exact"])
        for name in mcs.ATOM_FIELDS.values():
            self.assertEqual(getattr(params.AtomCompareParameters, name),
                             3.0 if name == "MaxDistance" else True)
        for name in mcs.BOND_FIELDS.values():
            self.assertTrue(getattr(params.BondCompareParameters, name))
        with self.assertRaises(ValueError):
            mcs.parameters({"misspelled": True})

    def test_parallel_order_and_full_result_fields(self):
        recipes = [{"format": "smiles", "text": text, "sanitize": True, "remove_hs": True}
                   for text in ("CCO", "CCN")]
        rows = [({"case_id": f"case_{i}", "inputs": [0, 1],
                  "parameters": {"timeout": 30, "verbose": verbose}, "source": {}}, recipes)
                for i, verbose in enumerate((False, True))]
        expected = [mcs.mcs_case(row) for row in rows]
        for threads in (1, 2):
            actual = reference.parallel(mcs.mcs_case, rows, threads, lambda *_: None)
            self.assertEqual(actual, expected)
        self.assertEqual(expected[0]["result"]["smarts"], "[#6]-[#6]")
        self.assertEqual(expected[0]["result"]["query_matches"], [True, True])
        self.assertTrue(expected[0]["result"]["completed"])
        self.assertEqual(expected[0]["result"]["query"]["atom_count"], 2)

    def test_parse_errors_remain_explicit_rows(self):
        case = {"case_id": "invalid", "inputs": [0, 1], "parameters": {"timeout": 30}, "source": {}}
        row = mcs.mcs_case((case, [{"format": "smiles", "text": "C#", "sanitize": True, "remove_hs": True}]))
        self.assertEqual(row["status"], "error")
        self.assertEqual(row["case_id"], "invalid")
        self.assertEqual(row["error"]["type"], "ValueError")


if __name__ == "__main__":
    unittest.main()
