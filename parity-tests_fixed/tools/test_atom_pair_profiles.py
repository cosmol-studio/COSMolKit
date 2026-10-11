"""Native-call adapter checks, not historical output reuse."""
import json
import unittest
from rdkit import Chem
from reference import ROOT
from atom_pair_profiles import atom_pair_case
from _generate_atom_pair_fingerprint_golden import (
    make_generator, output_record, resolved_arguments,
)


class AtomPairCallContract(unittest.TestCase):
    def test_all_historical_arguments_and_four_output_forms(self):
        profile = json.loads((ROOT / "tools/testdata/rdkit/atom_pair_fingerprint_profile.json").read_text())
        calls = 0
        for smiles in ("CCO", "C[C@H](N)C(=O)O", "c1ccccc1", ""):
            mol = Chem.MolFromSmiles(smiles)
            for branch in profile["branches"]:
                for form, method, kind in (
                    ("SparseCount", "GetSparseCountFingerprint", "count"),
                    ("Count", "GetCountFingerprint", "count"),
                    ("SparseBits", "GetSparseFingerprint", "bits"),
                    ("Bits", "GetFingerprint", "bits"),
                ):
                    p = dict(form=form, min_distance=branch["minDistance"],
                             max_distance=branch["maxDistance"],
                             include_chirality=branch["includeChirality"], use_2d=branch["use2D"],
                             count_simulation=branch["countSimulation"], count_bounds=branch["countBounds"],
                             fp_size=branch["fpSize"], bits_per_feature=branch["numBitsPerFeature"],
                             from_first=branch.get("fromAtoms") == "first",
                             ignore_first=branch.get("ignoreAtoms") == "first",
                             index_plus_11=branch.get("customAtomInvariants") == "index_plus_11",
                             additional_output=branch.get("additionalOutput", False))
                    row = {"Fingerprint": {"case": {"id": "focused", "smiles": smiles},
                                           "right": None, "seed": 1, "params": {"AtomPair": p}}}
                    old = output_record(make_generator(branch), method, kind, mol,
                                        resolved_arguments(mol, branch), p["additional_output"])
                    observed = atom_pair_case(row)["output"]["Fingerprint"]["AtomPair"]
                    with self.subTest(smiles=smiles, branch=branch["name"], form=form):
                        if not old["ok"]:
                            self.assertTrue(old["error"].startswith("IndexError: "))
                            self.assertEqual(observed, {"IndexError": {
                                "index": int(old["error"].split(": ", 1)[1]),
                                "size": 1 << (23 + (4 if p["include_chirality"] else 0))}})
                        else:
                            expected = {"length": old["length"]}
                            if kind == "count":
                                expected["entries"] = sorted([int(k), v] for k, v in old["nonzero_elements"].items())
                            else:
                                expected["on_bits"] = old["on_bits"]
                            self.assertEqual(observed, {"Counts" if kind == "count" else "Bits": {
                                "fingerprint": expected,
                                "additional_output": old.get("additional_output")}})
                    calls += 1
        self.assertEqual(calls, 160)

    def test_invalid_parse_is_not_a_fingerprint_or_bounds_error(self):
        row = {"Fingerprint": {"case": {"id": "bad", "smiles": "CC("},
                               "right": None, "seed": 1, "params": {"AtomPair": {}}}}
        self.assertEqual(atom_pair_case(row)["output"],
                         {"Fingerprint": {"ParseRejected": {"right": False}}})


if __name__ == "__main__":
    unittest.main()
