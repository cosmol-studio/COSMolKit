"""Validate every historical Torsion argument against independent native calls."""
import json
import unittest
from rdkit import Chem
from reference import ROOT
from torsion_profiles import torsion_case
from _generate_topological_torsion_fingerprint_golden import (
    make_generator, call_with_optional_output,
)


class TorsionCallContract(unittest.TestCase):
    def test_all_historical_arguments_forms_and_provenance(self):
        recipes = json.loads((ROOT / "tools/testdata/rdkit/topological_torsion_fingerprint_profile.json").read_text())
        calls = 0
        for smiles in ("CCCC", "C[C@H](N)C(=O)O", "c1ccccc1", ""):
            mol = Chem.MolFromSmiles(smiles)
            for branch in recipes["corpus_branches"]:
                for form, method in (
                    ("SparseCount", "GetSparseCountFingerprint"), ("Count", "GetCountFingerprint"),
                    ("SparseBits", "GetSparseFingerprint"), ("Bits", "GetFingerprint"),
                ):
                    p = dict(form=form, torsion_atom_count=branch.get("torsionAtomCount", 4),
                             only_shortest_paths=branch.get("onlyShortestPaths", False),
                             include_chirality=branch.get("includeChirality", False),
                             count_simulation=branch.get("countSimulation", True),
                             count_bounds=branch.get("countBounds", [1, 2, 4, 8]),
                             fp_size=branch.get("fpSize", 2048), bits_per_feature=branch.get("numBitsPerFeature", 1),
                             from_first=branch.get("fromAtoms") == "first",
                             ignore_first=branch.get("ignoreAtoms") == "first",
                             index_plus_17=branch.get("customAtomInvariants") == "index_plus_17",
                             additional_output=branch.get("additionalOutput", False))
                    row = {"Fingerprint": {"case": {"id": "focused", "smiles": smiles},
                                           "right": None, "seed": 1, "params": {"TopologicalTorsion": p}}}
                    fp, extra = call_with_optional_output(make_generator(branch), method, mol, branch)
                    if form in ("SparseCount", "Count"):
                        expected = {"Counts": {"fingerprint": {
                            "length": fp.GetLength(),
                            "entries": sorted([k, v] for k, v in fp.GetNonzeroElements().items())},
                            "additional_output": extra}}
                    else:
                        expected = {"Bits": {"fingerprint": {"length": fp.GetNumBits(),
                            "on_bits": sorted(i & 0xffffffff for i in fp.GetOnBits())},
                            "additional_output": extra}}
                    with self.subTest(smiles=smiles, branch=branch["name"], form=form):
                        self.assertEqual(torsion_case(row)["output"]["Fingerprint"]["TopologicalTorsion"], expected)
                    calls += 1
        self.assertEqual(calls, 144)

    def test_invalid_parse_retains_the_original_operand(self):
        row = {"Fingerprint": {"case": {"id": "bad", "smiles": "CC("},
                               "right": None, "seed": 1, "params": {"TopologicalTorsion": {}}}}
        self.assertEqual(torsion_case(row)["output"],
                         {"Fingerprint": {"ParseRejected": {"right": False}}})


if __name__ == "__main__":
    unittest.main()
