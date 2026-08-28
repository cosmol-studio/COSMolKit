from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from _generate_molecular_descriptors_golden import build_record, iter_smiles


class MolecularDescriptorGenerationTests(unittest.TestCase):
    def test_focused_input_represents_empty_smiles_without_ambiguous_blank_rows(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "focused.smi"
            path.write_text("# source\n\n<EMPTY>\nC\n", encoding="utf-8")
            self.assertEqual(list(iter_smiles(path)), ["", "C"])

    def test_high_feasibility_record_has_complete_fixed_shapes_and_bit_trees(self) -> None:
        record = build_record("CC(O)c1ccncc1")
        self.assertTrue(record["rdkit_ok"])

        values = record["high_feasibility_descriptors"]
        bits = record["high_feasibility_descriptor_bits"]
        assert isinstance(values, dict)
        assert isinstance(bits, dict)
        self.assertEqual(len(values["chi_nv_orders_0_6"]), 7)
        self.assertEqual(len(values["chi_nn_orders_0_6"]), 7)
        self.assertEqual(len(values["mqns"]), 42)
        self.assertEqual(len(values["slogp_vsa"]), 12)
        self.assertEqual(len(values["smr_vsa"]), 10)
        self.assertEqual(len(bits["slogp_vsa"]), 12)
        self.assertTrue(all(len(value) == 16 for value in bits["slogp_vsa"]))

        contributions = record["high_feasibility_contributions"]
        contribution_bits = record["high_feasibility_contribution_bits"]
        assert isinstance(contributions, dict)
        assert isinstance(contribution_bits, dict)
        self.assertEqual(
            len(contributions["hall_kier_alpha"]["atom_contributions"]), 9
        )
        self.assertEqual(
            len(
                contributions["labute_asa"]["include_hs_true"][
                    "atom_contributions"
                ]
            ),
            9,
        )
        self.assertEqual(
            len(contribution_bits["hall_kier_alpha"]["atom_contributions"]), 9
        )

        profiles = record["high_feasibility_cache_profiles"]
        profile_bits = record["high_feasibility_cache_profile_bits"]
        assert isinstance(profiles, dict)
        assert isinstance(profile_bits, dict)
        self.assertEqual(len(profiles["chi_nv"]["cold"]), 7)
        self.assertEqual(len(profiles["labute_asa_sequence"]), 4)
        self.assertEqual(len(profiles["slogp_vsa_custom_forced"]), 6)
        self.assertEqual(len(profile_bits["smr_vsa_default_warm"]), 10)

    def test_rejected_row_preserves_schema_without_descriptor_payloads(self) -> None:
        record = build_record("C1(")
        self.assertFalse(record["rdkit_ok"])
        for field in (
            "high_feasibility_descriptors",
            "high_feasibility_descriptor_bits",
            "high_feasibility_contributions",
            "high_feasibility_contribution_bits",
            "high_feasibility_cache_profiles",
            "high_feasibility_cache_profile_bits",
        ):
            self.assertIsNone(record[field])


if __name__ == "__main__":
    unittest.main()
