"""Reference-adapter boundaries, not a second force-field corpus."""
import struct
import unittest
from unittest.mock import patch

from reference import persistent_forcefield_case


def recipe(smiles, kind):
    bits = lambda value: struct.unpack("<Q", struct.pack("<d", value))[0]
    return {"PersistentForceField": {
        "case": {"id": "boundary", "smiles": smiles}, "kind": kind,
        "seed": 0x434B464620261008, "max_iterations": 2,
        "force_tolerance_bits": bits(1e-4),
        "energy_tolerance_bits": bits(1e-6), "preparation": None,
    }}


class PersistentForceFieldReferenceTests(unittest.TestCase):
    def test_syntax_rejection_has_no_geometry(self):
        for kind in ("Mmff", "Uff"):
            request = recipe("CC(", kind)
            record = persistent_forcefield_case(request)
            self.assertEqual(record["input"], request)
            self.assertEqual(record["output"]["PersistentForceField"], "ParseRejected")

    def test_exact_positions_do_not_depend_on_mol_writer(self):
        with patch("rdkit.Chem.MolToMolBlock", side_effect=AssertionError("unrelated writer")):
            for kind in ("Mmff", "Uff"):
                record = persistent_forcefield_case(recipe("CCO", kind))
                geometry = record["input"]["PersistentForceField"]["preparation"]
                initial = record["output"]["PersistentForceField"]["Evaluated"]["initial"]
                self.assertEqual(geometry["atom_count"], 9)
                self.assertEqual(initial["positions_bits"], geometry["coordinate_rows"][0]["xyz_bits"])

    def test_high_charge_is_evaluated_or_parameterized_not_transport_rejected(self):
        smiles = "[C+9]" + "(F)" * 12 + "F"
        mmff = persistent_forcefield_case(recipe(smiles, "Mmff"))
        uff = persistent_forcefield_case(recipe(smiles, "Uff"))
        self.assertEqual(mmff["output"]["PersistentForceField"], "Unavailable")
        self.assertIn("Evaluated", uff["output"]["PersistentForceField"])
        self.assertEqual(uff["input"]["PersistentForceField"]["preparation"]["atom_count"], 14)

    def test_source_reports_first_missing_center_when_several_are_missing(self):
        record = persistent_forcefield_case(recipe("C~C", "Uff"))
        self.assertEqual(record["output"]["PersistentForceField"], {
            "SourceTbpCenterParamsMissing": {"center_atom_index": 0}})


if __name__ == "__main__":
    unittest.main()
