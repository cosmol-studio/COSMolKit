"""Canonical Python descriptor parity using the unchanged original 5000 reference.

The original reference's 41-character revision spelling is retained provenance,
not an assertion that a new native RDKit oracle was built. No data is regenerated.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import struct

import cosmolkit as ck
import pytest

ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT / "testdata/descriptors/expected/rdkit/smiles_5000"
SCALARS = {
    "chi_0": "chi_0", "chi_1": "chi_1", "fraction_csp3": "fraction_csp3",
    "hall_kier_alpha": "hall_kier_alpha",
    "kappa_1": "kappa_1", "kappa_2": "kappa_2", "kappa_3": "kappa_3", "phi": "phi",
    **{f"chi_{order}_{kind}": f"chi_{order}{kind}" for kind in ("n", "v") for order in range(5)},
}
COUNTS = {
    **{name: name for name in ("num_amide_bonds", "num_spiro_atoms", "num_bridgehead_atoms", "num_atom_stereo_centers", "num_unspecified_atom_stereo_centers")},
    "num_heavy_atoms": "num_heavy_atoms",
    "total_atom_count": "num_atoms",
    "lipinski_hba": "lipinski_hba",
    "lipinski_hbd": "lipinski_hbd",
    "num_heteroatoms": "num_heteroatoms",
    "num_hba": "num_hba",
    "num_hbd": "num_hbd",
    "num_rings": "num_rings",
    "num_heterocycles": "num_heterocycles",
    "num_aromatic_rings": "num_aromatic_rings",
    "num_saturated_rings": "num_saturated_rings",
    "num_aliphatic_rings": "num_aliphatic_rings",
    "num_aromatic_heterocycles": "num_aromatic_heterocycles",
    "num_aromatic_carbocycles": "num_aromatic_carbocycles",
    "num_aliphatic_heterocycles": "num_aliphatic_heterocycles",
    "num_aliphatic_carbocycles": "num_aliphatic_carbocycles",
    "num_saturated_heterocycles": "num_saturated_heterocycles",
    "num_saturated_carbocycles": "num_saturated_carbocycles",
}
PROFILES = [(name, None) for name in SCALARS]
PROFILES += [(name, None) for name in COUNTS]
PROFILES += [("hall_kier_alpha_with_contributions", None), ("mqns", False), ("mqns", True)]
PROFILES += [(f"chi_n_{kind}", order) for kind in ("n", "v") for order in range(7)]

PROFILES += [("num_rotatable_bonds", None)]
PROFILES += [("num_rotatable_bonds_with_params", mode) for mode in ("Default", "NonStrict", "Strict", "StrictLinkages")]
PROFILES += [(name, None) for name in ("molecular_weight", "exact_molecular_weight", "molecular_formula")]
PROFILES += [(name, heavy) for name in ("molecular_weight_with_params", "exact_molecular_weight_with_params") for heavy in (False, True)]
PROFILES += [("molecular_formula_with_params", flags) for flags in ((False, False), (False, True), (True, False), (True, True))]



PROFILES += [(name, None) for name in ("crippen_descriptors", "labute_asa", "labute_asa_contributions", "tpsa", "slogp_vsa", "smr_vsa")]
PROFILES += [(f"{family}_{index}", None) for family, count in (("slogp_vsa", 12), ("smr_vsa", 10)) for index in range(1, count + 1)]
PROFILES += [(name + "_with_params", (include, force)) for name in ("crippen_descriptors", "labute_asa", "labute_asa_contributions", "tpsa") for include in (False, True) for force in (False, True)]
PROFILES += [(name + "_with_params", (boundaries, force)) for name in ("slogp_vsa", "smr_vsa") for boundaries, force in ((None, False), (None, True), ([-0.2, 0.0, 0.25, 0.25, 0.8], True))]

PROFILES += [("qed", None)]
PROFILES += [(f"chi_{order}_{kind}_with_params", force) for kind in ("n", "v") for order in range(5) for force in (False, True)]
PROFILES += [(f"chi_n_{kind}_with_params", (order, force)) for kind in ("n", "v") for order in range(7) for force in (False, True)]

def verified(path: Path, sha: str) -> bytes:
    data = path.read_bytes()
    assert hashlib.sha256(data).hexdigest() == sha, f"Changed original reference: {path}"
    return data


@pytest.fixture(scope="module")
def reference():
    # Complete input preflight happens before any CK chemistry call.
    raw = verified(ROOT / "testdata/smiles/corpus/smiles_5000.smi", "a4d579cd72621af27772256bb23ba796452276bb924fd20aac83625ffa67d849")
    manifest = json.loads(verified(DIRECTORY / "manifest.json", "7da152d6f7a454d2dddecb7356e847fd61367690cdab6c91cae69896212452aa"))
    golden = verified(DIRECTORY / "molecular_descriptors.jsonl", "9b3aa0b3c5ca537388503dc4d3d495a2d04bebe0d64bf09fda0225ce452274c8")
    rows = [json.loads(line) for line in golden.splitlines()]
    inputs = [line.decode().split()[0] for line in raw.splitlines() if line.strip()]
    assert len(rows) == len(inputs) == 5000
    assert manifest["input"]["sha256"] == hashlib.sha256(raw).hexdigest()
    assert manifest["outputs"][0]["sha256"] == hashlib.sha256(golden).hexdigest()
    for smiles, row in zip(inputs, rows, strict=True):
        assert row["smiles"] == smiles and row["rdkit_ok"] is True and row["error"] is None
        for field in SCALARS.values():
            assert len(row["high_feasibility_descriptor_bits"][field]) == 16
        for field in COUNTS.values():
            assert isinstance(row["high_feasibility_descriptors"][field], int)
        for field in ("mol_wt", "exact_mol_wt"):
            assert len(row["descriptor_bits"][field]) == 16
            assert len(row["descriptor_option_bits"][field + "_only_heavy"]) == 16
        for mode in ("default", "non_strict", "strict", "strict_linkages"):
            assert isinstance(row["descriptors"]["num_rotatable_bonds_" + mode], int)
        for field in ("formula", "formula_separate_isotopes", "formula_separate_isotopes_no_h_abbrev"):
            assert isinstance(row["descriptors"][field], str)
        for kind in ("n", "v"):
            assert len(row["high_feasibility_descriptor_bits"]["chi_nn_orders_0_6" if kind == "n" else "chi_nv_orders_0_6"]) == 7
    return rows


def bits(value):
    assert isinstance(value, float) and math.isfinite(value)
    return struct.pack(">d", value).hex()


@pytest.mark.parametrize("name,argument", PROFILES, ids=[f"{name}-{argument}" for name, argument in PROFILES])
def test_canonical_descriptor_original5000(reference, name, argument):
    for index, row in enumerate(reference, 1):
        molecule = ck.Molecule.from_smiles(row["smiles"])
        query = getattr(molecule, name)
        if name == "num_rotatable_bonds_with_params":
            actual = query(getattr(ck.RotatableBondsOptions, argument))
        elif name.endswith("_with_params") and isinstance(argument, tuple):
            actual = query(*argument)
        else:
            actual = query() if argument is None else query(argument)
        context = (name, argument, index, row["smiles"])
        if name == "qed":
            assert bits(actual) == row["descriptor_bits"]["qed"], context
        elif name.startswith("chi_") and name.endswith("_with_params"):
            base = name.removesuffix("_with_params")
            if base in ("chi_n_n", "chi_n_v"):
                field = "chi_nv_orders_0_6" if base.endswith("v") else "chi_nn_orders_0_6"
                expected = row["high_feasibility_descriptor_bits"][field][argument[0]]
            else:
                order, kind = base.split("_")[1:]
                expected = row["high_feasibility_descriptor_bits"]["chi_" + order + kind]
            assert bits(actual) == expected, context
        elif name.startswith("crippen_descriptors"):
            expected = {"logp": row["descriptor_bits"]["crippen_logp"], "molar_refractivity": row["descriptor_bits"]["crippen_mr"]} if argument is None else row["descriptor_option_bits"]["crippen"][f"include_hs_{str(argument[0]).lower()}_force_{str(argument[1]).lower()}"]
            assert isinstance(actual, ck.CrippenTotals)
            assert bits(actual.logp) == expected["logp"], context
            assert bits(actual.molar_refractivity) == expected["molar_refractivity"], context
        elif name.startswith("labute_asa"):
            include = True if argument is None else argument[0]
            expected = row["high_feasibility_contribution_bits"]["labute_asa"][f"include_hs_{str(include).lower()}"]
            if "contributions" in name:
                assert isinstance(actual, ck.LabuteAsaContributions)
                assert bits(actual.asa) == expected["asa"], context
                assert [bits(v) for v in actual.atom_contributions] == expected["atom_contributions"], context
                assert bits(actual.hydrogen_contribution) == expected["hydrogen_contribution"], context
            else:
                assert bits(actual) == expected["asa"], context
        elif name.startswith("tpsa"):
            expected = row["descriptor_bits"]["tpsa"] if argument is None else row["descriptor_option_bits"]["tpsa"][f"force_{str(argument[1]).lower()}_include_sandp_{str(argument[0]).lower()}"]
            assert bits(actual) == expected, context
        elif name.startswith(("slogp_vsa", "smr_vsa")):
            family = "slogp_vsa" if name.startswith("slogp_vsa") else "smr_vsa"
            if argument is not None and argument[0] is not None:
                expected = row["high_feasibility_cache_profile_bits"][family + "_custom_forced"]
            else:
                expected = row["high_feasibility_descriptor_bits"][family]
            suffix = name.removeprefix(family + "_")
            if suffix.isdigit():
                assert bits(actual) == expected[int(suffix) - 1], context
            else:
                assert isinstance(actual, list) and [bits(v) for v in actual] == expected, context
        elif name.startswith("num_rotatable_bonds"):
            suffix = {"Default": "default", "NonStrict": "non_strict", "Strict": "strict", "StrictLinkages": "strict_linkages"}[argument or "Default"]
            assert isinstance(actual, int) and not isinstance(actual, bool)
            assert actual == row["descriptors"]["num_rotatable_bonds_" + suffix], context
        elif name.startswith("molecular_formula"):
            field = "formula"
            if argument is not None and argument[0]:
                field = "formula_separate_isotopes" if argument[1] else "formula_separate_isotopes_no_h_abbrev"
            assert isinstance(actual, str) and actual == row["descriptors"][field], context
        elif name.startswith("molecular_weight") or name.startswith("exact_molecular_weight"):
            field = "exact_mol_wt" if name.startswith("exact_") else "mol_wt"
            expected = row["descriptor_option_bits"][field + "_only_heavy"] if argument is True else row["descriptor_bits"][field]
            assert bits(actual) == expected, context
        elif name in SCALARS:
            assert bits(actual) == row["high_feasibility_descriptor_bits"][SCALARS[name]], context
        elif name in COUNTS:
            assert isinstance(actual, int) and not isinstance(actual, bool)
            assert actual == row["high_feasibility_descriptors"][COUNTS[name]], context
        elif name == "mqns":
            assert isinstance(actual, list) and len(actual) == 42
            assert actual == row["high_feasibility_descriptors"]["mqns"], context
        elif name == "hall_kier_alpha_with_contributions":
            expected = row["high_feasibility_contribution_bits"]["hall_kier_alpha"]
            assert isinstance(actual, tuple) and len(actual) == 2
            assert bits(actual[0]) == expected["value"], context
            assert [bits(value) for value in actual[1]] == expected["atom_contributions"], context
        else:
            field = "chi_nv_orders_0_6" if name == "chi_n_v" else "chi_nn_orders_0_6"
            assert bits(actual) == row["high_feasibility_descriptor_bits"][field][argument], context


def test_labute_cache_source_sequence_original5000(reference):
    for index, row in enumerate(reference, 1):
        molecule = ck.Molecule.from_smiles(row["smiles"])
        for step in row["high_feasibility_cache_profile_bits"]["labute_asa_sequence"]:
            actual = molecule.labute_asa_with_params(step["include_hs"], step["force"])
            assert bits(actual) == step["value"], (index, row["smiles"], step)


def test_vsa_shared_cache_source_sequence_original5000(reference):
    for index, row in enumerate(reference, 1):
        molecule = ck.Molecule.from_smiles(row["smiles"])
        expected = row["high_feasibility_cache_profile_bits"]
        calls = (
            (molecule.slogp_vsa, (), "slogp_vsa_default_cold"),
            (molecule.slogp_vsa, (), "slogp_vsa_default_warm"),
            (molecule.slogp_vsa_with_params, (expected["custom_bins"], True), "slogp_vsa_custom_forced"),
            (molecule.smr_vsa, (), "smr_vsa_default_warm"),
            (molecule.smr_vsa_with_params, (expected["custom_bins"], True), "smr_vsa_custom_forced"),
        )
        # custom_bins in this bit tree is a bit-vector, so use the original
        # numeric boundary vector; repeated edges are retained exactly.
        boundaries = row["high_feasibility_cache_profiles"]["custom_bins"]
        for query, args, name in calls:
            if args:
                args = (boundaries, True)
            assert [bits(v) for v in query(*args)] == expected[name], (index, row["smiles"], name)


def test_chi_cache_source_sequence_original5000(reference):
    for index, row in enumerate(reference, 1):
        molecule = ck.Molecule.from_smiles(row["smiles"])
        for kind in ("n", "v"):
            expected = row["high_feasibility_cache_profile_bits"]["chi_n" + kind]
            query = getattr(molecule, f"chi_n_{kind}_with_params")
            for phase, force in (("cold", False), ("warm", False), ("forced", True)):
                assert [bits(query(order, force)) for order in range(7)] == expected[phase], (index, row["smiles"], kind, phase)
