"""Validate historical calls against the original ChEMBL recipe definitions."""
import ast
import unittest
from typing import Any
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from reference import ROOT
from morgan_profiles import morgan_case


def original_recipes():
    source = ROOT / "dev/tools/chembl_parity/audit_core.py"
    tree = ast.parse(source.read_text())
    selected = [node for node in tree.body if
                isinstance(node, ast.FunctionDef) and node.name in (
                    "morgan_generator_with_num_bits_per_feature", "dynamic_morgan_branches")
                or isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name)
                and node.target.id == "MORGAN_BRANCHES"]
    if len(selected) != 3:
        raise AssertionError("historical ChEMBL recipe source changed")
    namespace = {"Chem": Chem, "Any": Any,
                 "rdFingerprintGenerator": rdFingerprintGenerator}
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(source), "exec"), namespace)
    return namespace


class MorganHistoricalCalls(unittest.TestCase):
    def test_every_source_recipe_and_eligibility_against_independent_calls(self):
        source = original_recipes()
        count = 0
        for smiles, expected in (("", 10), ("C", 14), ("CC", 15), ("C.C", 14),
                                 ("C[C@H](N)C(=O)O", 15), ("F/C=C/F", 15), ("c1ccccc1", 15)):
            mol = Chem.MolFromSmiles(smiles)
            branches = [(name, generator, {}, ck) for name, (generator, ck)
                        in source["MORGAN_BRANCHES"].items()]
            branches.extend(source["dynamic_morgan_branches"](mol))
            self.assertEqual(len(branches), expected)
            for name, generator, kwargs, _ in branches:
                extra = rdFingerprintGenerator.AdditionalOutput()
                extra.AllocateAtomCounts()
                extra.AllocateAtomToBits()
                extra.AllocateBitInfoMap()
                extra.AllocateBitPaths()
                extra.AllocateAtomsPerBit()
                fp = generator.GetFingerprint(mol, additionalOutput=extra, **kwargs)
                def rows(mapping):
                    return [[k, [list(v) for v in values]] for k, values in sorted(mapping.items())]
                expected_value = {
                    "atom_count": mol.GetNumAtoms(), "bond_count": mol.GetNumBonds(),
                    "fingerprint": {"length": fp.GetNumBits(), "on_bits": list(fp.GetOnBits()),
                        "additional_output": {
                            "atom_counts": list(extra.GetAtomCounts()),
                            "atom_to_bits": [list(v) for v in extra.GetAtomToBits()],
                            "bit_info_map": rows(extra.GetBitInfoMap()),
                            "bit_paths": rows(extra.GetBitPaths()),
                            "atoms_per_bit": rows(extra.GetAtomsPerBit()),
                        }}}
                row = {"Molecular": {"case": {"id": "focused", "smiles": smiles},
                                     "profile": {"MorganHistorical": {"recipe": name}}}}
                with self.subTest(smiles=smiles, recipe=name):
                    self.assertEqual(morgan_case(row)["output"]["Molecular"]["MorganHistorical"], expected_value)
                count += 1
        self.assertEqual(count, 98)

    def test_invalid_parse_and_ordinary_delegation(self):
        row = {"Molecular": {"case": {"id": "bad", "smiles": "CC("},
                             "profile": {"MorganHistorical": {"recipe": "default"}}}}
        self.assertIn("ParseRejected", morgan_case(row)["output"]["Molecular"])
        row["Molecular"]["case"]["smiles"] = "CCO"
        row["Molecular"]["profile"] = {"Morgan": {"output": "DenseBits", "radius": 2,
            "include_chirality": False, "invariants": "Connectivity", "count_simulation": False}}
        self.assertIn("MorganDenseBits", morgan_case(row)["output"]["Molecular"])


if __name__ == "__main__":
    unittest.main()
