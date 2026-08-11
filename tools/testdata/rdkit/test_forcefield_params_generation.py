from __future__ import annotations

import unittest

from rdkit import Chem
from rdkit.Chem import AllChem

from _generate_forcefield_params_golden import (
    FORCEFIELD_NONBONDED_THRESH,
    FORCEFIELD_OPT_MAX_ITERS,
    FORCEFIELD_PARITY_SEED,
    forcefield_initial_energy_result,
    forcefield_optimized_result,
    mmff_result,
)
from _generate_forcefield_coverage_golden import build_record as build_coverage_record


class ForceFieldParamsGenerationTests(unittest.TestCase):
    def test_coverage_record_excludes_embedding_and_optimizer_surfaces(self) -> None:
        record = build_coverage_record("CCO")

        self.assertTrue(record["rdkit_ok"])
        self.assertNotIn("embedded", record)
        self.assertEqual(
            set(record),
            {
                "smiles",
                "rdkit_ok",
                "uff",
                "mmff",
                "uff_explicit_h",
                "mmff_explicit_h",
                "error",
            },
        )
        self.assertTrue(record["uff"]["ok"])
        self.assertTrue(record["mmff"]["ok"])
        self.assertTrue(record["uff_explicit_h"]["ok"])
        self.assertTrue(record["mmff_explicit_h"]["ok"])

    def test_unsupported_mmff_surfaces_do_not_mutate_later_uff_input(self) -> None:
        mol = Chem.MolFromSmiles("c1cc[se]c1", sanitize=True)
        self.assertIsNotNone(mol)
        assert mol is not None
        self.assertEqual(
            AllChem.EmbedMolecule(
                mol,
                randomSeed=FORCEFIELD_PARITY_SEED,
                useRandomCoords=True,
            ),
            0,
        )
        parity_mol = Chem.MolFromSmiles(
            Chem.MolToCXSmiles(mol, isomericSmiles=True),
            sanitize=True,
        )
        self.assertIsNotNone(parity_mol)
        assert parity_mol is not None

        expected_smiles = "c1cc[se]c1"
        expected_aromatic = [True] * parity_mol.GetNumAtoms()
        self.assertEqual(Chem.MolToSmiles(parity_mol), expected_smiles)
        self.assertEqual(
            [atom.GetIsAromatic() for atom in parity_mol.GetAtoms()],
            expected_aromatic,
        )

        properties = mmff_result(parity_mol)
        self.assertTrue(properties["ok"])
        self.assertFalse(properties["has_all"])
        self.assertIsNone(properties["atom_types"])

        initial = forcefield_initial_energy_result(
            lambda candidate: AllChem.MMFFGetMoleculeForceField(
                candidate,
                AllChem.MMFFGetMoleculeProperties(candidate, mmffVariant="MMFF94"),
                nonBondedThresh=FORCEFIELD_NONBONDED_THRESH,
                confId=0,
                ignoreInterfragInteractions=True,
            ),
            parity_mol,
        )
        self.assertTrue(initial["ok"])
        self.assertEqual(initial["needs_more"], -1)
        self.assertEqual(initial["energy"], -1.0)
        self.assertEqual(Chem.MolToSmiles(parity_mol), expected_smiles)
        self.assertEqual(
            [atom.GetIsAromatic() for atom in parity_mol.GetAtoms()],
            expected_aromatic,
        )

        optimized = forcefield_optimized_result(
            lambda candidate: AllChem.UFFGetMoleculeForceField(
                candidate,
                vdwThresh=FORCEFIELD_NONBONDED_THRESH,
                confId=0,
                ignoreInterfragInteractions=True,
            ),
            parity_mol,
            FORCEFIELD_OPT_MAX_ITERS,
        )
        self.assertTrue(optimized["ok"])
        self.assertEqual(optimized["needs_more"], 0)
        self.assertAlmostEqual(optimized["energy"], 8.341918182756241e-13, places=18)
        self.assertAlmostEqual(optimized["coords"][0][0], 1.3370117772556993, places=12)
        self.assertNotAlmostEqual(optimized["coords"][0][0], 1.3606108238464545, places=6)


if __name__ == "__main__":
    unittest.main()
