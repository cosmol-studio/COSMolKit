#!/usr/bin/env python3
"""Generate RDKit molecular-descriptor golden data from a SMILES corpus."""

from __future__ import annotations

import argparse
import json
import struct
from pathlib import Path
from typing import Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, QED, rdMolDescriptors
from rdkit.Chem import GraphDescriptors

RDLogger.DisableLog("rdApp.*")


def iter_smiles(path: Path) -> Iterable[str]:
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        yield "" if line == "<EMPTY>" else line


def descriptor_set(mol: Chem.Mol) -> dict[str, object]:
    descriptors = {
        "mol_wt": Descriptors.MolWt(mol),
        "exact_mol_wt": Descriptors.ExactMolWt(mol),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "formula_separate_isotopes": rdMolDescriptors.CalcMolFormula(
            mol, separateIsotopes=True
        ),
        "formula_separate_isotopes_no_h_abbrev": rdMolDescriptors.CalcMolFormula(
            mol, separateIsotopes=True, abbreviateHIsotopes=False
        ),
        "num_hbd": rdMolDescriptors.CalcNumHBD(mol),
        "num_hba": rdMolDescriptors.CalcNumHBA(mol),
        "fraction_csp3": rdMolDescriptors.CalcFractionCSP3(mol),
        "crippen_logp": Crippen.MolLogP(mol),
        "crippen_mr": Crippen.MolMR(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "tpsa_include_sandp": rdMolDescriptors.CalcTPSA(
            mol, force=True, includeSandP=True
        ),
        "num_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "num_rotatable_bonds_default": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "num_rotatable_bonds_non_strict": rdMolDescriptors.CalcNumRotatableBonds(
            mol, rdMolDescriptors.NumRotatableBondsOptions.NonStrict
        ),
        "num_rotatable_bonds_strict": rdMolDescriptors.CalcNumRotatableBonds(
            mol, rdMolDescriptors.NumRotatableBondsOptions.Strict
        ),
        "num_rotatable_bonds_strict_linkages": rdMolDescriptors.CalcNumRotatableBonds(
            mol, rdMolDescriptors.NumRotatableBondsOptions.StrictLinkages
        ),
        "qed": QED.qed(mol),
    }
    return descriptors


def descriptor_bits(descriptors: dict[str, object]) -> dict[str, str]:
    float_fields = (
        "mol_wt",
        "exact_mol_wt",
        "fraction_csp3",
        "crippen_logp",
        "crippen_mr",
        "tpsa",
        "tpsa_include_sandp",
        "qed",
    )
    return {
        field: f"{struct.unpack('<Q', struct.pack('<d', descriptors[field]))[0]:016x}"
        for field in float_fields
    }


def float_bits(value: float) -> str:
    return f"{struct.unpack('<Q', struct.pack('<d', value))[0]:016x}"


def float_tree_bits(value: object) -> object:
    """Return a shape-preserving bit tree for every floating leaf."""
    if isinstance(value, float):
        return float_bits(value)
    if isinstance(value, list):
        return [float_tree_bits(item) for item in value]
    if isinstance(value, dict):
        return {key: float_tree_bits(item) for key, item in value.items()}
    return value


def high_feasibility_descriptor_set(mol: Chem.Mol) -> dict[str, object]:
    chi_nv = [rdMolDescriptors.CalcChiNv(mol, order, False) for order in range(7)]
    chi_nn = [rdMolDescriptors.CalcChiNn(mol, order, False) for order in range(7)]
    return {
        # Python's sum([]) returns int(0), while the descriptor surface is a
        # floating scalar. Normalize only the wrapper's dynamic empty type so
        # the golden records the C++/Rust f64 behavior uniformly.
        "chi_0": float(GraphDescriptors.Chi0(mol)),
        "chi_1": float(GraphDescriptors.Chi1(mol)),
        "chi_nv_orders_0_6": chi_nv,
        "chi_nn_orders_0_6": chi_nn,
        "chi_0v": rdMolDescriptors.CalcChi0v(mol, False),
        "chi_1v": rdMolDescriptors.CalcChi1v(mol, False),
        "chi_2v": rdMolDescriptors.CalcChi2v(mol, False),
        "chi_3v": rdMolDescriptors.CalcChi3v(mol, False),
        "chi_4v": rdMolDescriptors.CalcChi4v(mol, False),
        "chi_0n": rdMolDescriptors.CalcChi0n(mol, False),
        "chi_1n": rdMolDescriptors.CalcChi1n(mol, False),
        "chi_2n": rdMolDescriptors.CalcChi2n(mol, False),
        "chi_3n": rdMolDescriptors.CalcChi3n(mol, False),
        "chi_4n": rdMolDescriptors.CalcChi4n(mol, False),
        "hall_kier_alpha": rdMolDescriptors.CalcHallKierAlpha(mol),
        "kappa_1": rdMolDescriptors.CalcKappa1(mol),
        "kappa_2": rdMolDescriptors.CalcKappa2(mol),
        "kappa_3": rdMolDescriptors.CalcKappa3(mol),
        "phi": rdMolDescriptors.CalcPhi(mol),
        "lipinski_hba": rdMolDescriptors.CalcNumLipinskiHBA(mol),
        "lipinski_hbd": rdMolDescriptors.CalcNumLipinskiHBD(mol),
        "num_hba": rdMolDescriptors.CalcNumHBA(mol),
        "num_hbd": rdMolDescriptors.CalcNumHBD(mol),
        "num_heteroatoms": rdMolDescriptors.CalcNumHeteroatoms(mol),
        "num_amide_bonds": rdMolDescriptors.CalcNumAmideBonds(mol),
        "num_heavy_atoms": rdMolDescriptors.CalcNumHeavyAtoms(mol),
        "num_atoms": rdMolDescriptors.CalcNumAtoms(mol),
        "num_rings": rdMolDescriptors.CalcNumRings(mol),
        "num_heterocycles": rdMolDescriptors.CalcNumHeterocycles(mol),
        "num_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "num_saturated_rings": rdMolDescriptors.CalcNumSaturatedRings(mol),
        "num_aliphatic_rings": rdMolDescriptors.CalcNumAliphaticRings(mol),
        "num_aromatic_heterocycles": rdMolDescriptors.CalcNumAromaticHeterocycles(mol),
        "num_aromatic_carbocycles": rdMolDescriptors.CalcNumAromaticCarbocycles(mol),
        "num_aliphatic_heterocycles": rdMolDescriptors.CalcNumAliphaticHeterocycles(mol),
        "num_aliphatic_carbocycles": rdMolDescriptors.CalcNumAliphaticCarbocycles(mol),
        "num_saturated_heterocycles": rdMolDescriptors.CalcNumSaturatedHeterocycles(mol),
        "num_saturated_carbocycles": rdMolDescriptors.CalcNumSaturatedCarbocycles(mol),
        "num_spiro_atoms": rdMolDescriptors.CalcNumSpiroAtoms(mol),
        "num_bridgehead_atoms": rdMolDescriptors.CalcNumBridgeheadAtoms(mol),
        "num_atom_stereo_centers": rdMolDescriptors.CalcNumAtomStereoCenters(mol),
        "num_unspecified_atom_stereo_centers": (
            rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol)
        ),
        "fraction_csp3": rdMolDescriptors.CalcFractionCSP3(mol),
        "mqns": list(rdMolDescriptors.MQNs_(mol, False)),
        "labute_asa_include_hs_false": rdMolDescriptors.CalcLabuteASA(
            mol, False, False
        ),
        "labute_asa_include_hs_true": rdMolDescriptors.CalcLabuteASA(
            mol, True, True
        ),
        "slogp_vsa": list(rdMolDescriptors.SlogP_VSA_(mol, [], False)),
        "smr_vsa": list(rdMolDescriptors.SMR_VSA_(mol, [], False)),
    }


def high_feasibility_contributions(smiles: str) -> dict[str, object]:
    alpha_mol = fresh_mol(smiles)
    alpha_contributions = [0.0] * alpha_mol.GetNumAtoms()
    alpha = rdMolDescriptors.CalcHallKierAlpha(alpha_mol, alpha_contributions)

    labute: dict[str, object] = {}
    for include_hs in (False, True):
        mol = fresh_mol(smiles)
        atom_contributions, hydrogen_contribution = (
            rdMolDescriptors._CalcLabuteASAContribs(mol, include_hs)
        )
        key = f"include_hs_{str(include_hs).lower()}"
        labute[key] = {
            "asa": rdMolDescriptors.CalcLabuteASA(
                fresh_mol(smiles), include_hs, False
            ),
            "atom_contributions": list(atom_contributions),
            "hydrogen_contribution": hydrogen_contribution,
        }

    return {
        "hall_kier_alpha": {
            "value": alpha,
            "atom_contributions": alpha_contributions,
        },
        "labute_asa": labute,
    }


def high_feasibility_cache_profiles(smiles: str) -> dict[str, object]:
    def chi_profile(function: object) -> dict[str, object]:
        mol = fresh_mol(smiles)
        return {
            "cold": [function(mol, order, False) for order in range(7)],
            "warm": [function(mol, order, False) for order in range(7)],
            "forced": [function(mol, order, True) for order in range(7)],
        }

    labute_mol = fresh_mol(smiles)
    labute_sequence = [
        {
            "include_hs": False,
            "force": False,
            "value": rdMolDescriptors.CalcLabuteASA(labute_mol, False, False),
        },
        {
            "include_hs": True,
            "force": False,
            "value": rdMolDescriptors.CalcLabuteASA(labute_mol, True, False),
        },
        {
            "include_hs": True,
            "force": True,
            "value": rdMolDescriptors.CalcLabuteASA(labute_mol, True, True),
        },
        {
            "include_hs": False,
            "force": False,
            "value": rdMolDescriptors.CalcLabuteASA(labute_mol, False, False),
        },
    ]

    custom_bins = [-0.2, 0.0, 0.25, 0.25, 0.8]
    vsa_mol = fresh_mol(smiles)
    return {
        "chi_nv": chi_profile(rdMolDescriptors.CalcChiNv),
        "chi_nn": chi_profile(rdMolDescriptors.CalcChiNn),
        "labute_asa_sequence": labute_sequence,
        "slogp_vsa_default_cold": list(
            rdMolDescriptors.SlogP_VSA_(vsa_mol, [], False)
        ),
        "slogp_vsa_default_warm": list(
            rdMolDescriptors.SlogP_VSA_(vsa_mol, [], False)
        ),
        "slogp_vsa_custom_forced": list(
            rdMolDescriptors.SlogP_VSA_(vsa_mol, custom_bins, True)
        ),
        "smr_vsa_default_warm": list(rdMolDescriptors.SMR_VSA_(vsa_mol, [], False)),
        "smr_vsa_custom_forced": list(
            rdMolDescriptors.SMR_VSA_(vsa_mol, custom_bins, True)
        ),
        "custom_bins": custom_bins,
    }


def fresh_mol(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        raise RuntimeError("MolFromSmiles unexpectedly failed after initial success")
    return mol


def descriptor_option_matrix(smiles: str) -> tuple[dict[str, object], dict[str, object]]:
    mass_mol = fresh_mol(smiles)
    options: dict[str, object] = {
        "mol_wt_only_heavy": Descriptors.MolWt(mass_mol, onlyHeavy=True),
        "exact_mol_wt_only_heavy": Descriptors.ExactMolWt(
            mass_mol, onlyHeavy=True
        ),
        "crippen": {},
        "tpsa": {},
    }
    bits: dict[str, object] = {
        "mol_wt_only_heavy": float_bits(options["mol_wt_only_heavy"]),
        "exact_mol_wt_only_heavy": float_bits(
            options["exact_mol_wt_only_heavy"]
        ),
        "crippen": {},
        "tpsa": {},
    }

    crippen = options["crippen"]
    crippen_bits = bits["crippen"]
    assert isinstance(crippen, dict)
    assert isinstance(crippen_bits, dict)
    for include_hs in (False, True):
        for force in (False, True):
            key = f"include_hs_{str(include_hs).lower()}_force_{str(force).lower()}"
            # Use a fresh molecule for every branch. RDKit caches Crippen
            # properties on the molecule and a previous force=False call can
            # otherwise mask the includeHs argument of the next call.
            logp, mr = rdMolDescriptors.CalcCrippenDescriptors(
                fresh_mol(smiles), includeHs=include_hs, force=force
            )
            crippen[key] = {"logp": logp, "molar_refractivity": mr}
            crippen_bits[key] = {
                "logp": float_bits(logp),
                "molar_refractivity": float_bits(mr),
            }

    tpsa = options["tpsa"]
    tpsa_bits = bits["tpsa"]
    assert isinstance(tpsa, dict)
    assert isinstance(tpsa_bits, dict)
    for force in (False, True):
        for include_sandp in (False, True):
            key = (
                f"force_{str(force).lower()}_"
                f"include_sandp_{str(include_sandp).lower()}"
            )
            value = rdMolDescriptors.CalcTPSA(
                fresh_mol(smiles), force=force, includeSandP=include_sandp
            )
            tpsa[key] = value
            tpsa_bits[key] = float_bits(value)

    return options, bits


def build_record(smiles: str) -> dict[str, object]:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return {
            "smiles": smiles,
            "rdkit_ok": False,
            "descriptors": None,
            "descriptor_bits": None,
            "descriptor_options": None,
            "descriptor_option_bits": None,
            "high_feasibility_descriptors": None,
            "high_feasibility_descriptor_bits": None,
            "high_feasibility_contributions": None,
            "high_feasibility_contribution_bits": None,
            "high_feasibility_cache_profiles": None,
            "high_feasibility_cache_profile_bits": None,
            "error": "MolFromSmiles returned None",
        }
    descriptors = descriptor_set(mol)
    descriptor_options, descriptor_option_bits = descriptor_option_matrix(smiles)
    high_feasibility_descriptors = high_feasibility_descriptor_set(fresh_mol(smiles))
    high_feasibility_contribution_values = high_feasibility_contributions(smiles)
    high_feasibility_cache_profile_values = high_feasibility_cache_profiles(smiles)
    return {
        "smiles": smiles,
        "rdkit_ok": True,
        "descriptors": descriptors,
        "descriptor_bits": descriptor_bits(descriptors),
        "descriptor_options": descriptor_options,
        "descriptor_option_bits": descriptor_option_bits,
        "high_feasibility_descriptors": high_feasibility_descriptors,
        "high_feasibility_descriptor_bits": float_tree_bits(
            high_feasibility_descriptors
        ),
        "high_feasibility_contributions": high_feasibility_contribution_values,
        "high_feasibility_contribution_bits": float_tree_bits(
            high_feasibility_contribution_values
        ),
        "high_feasibility_cache_profiles": high_feasibility_cache_profile_values,
        "high_feasibility_cache_profile_bits": float_tree_bits(
            high_feasibility_cache_profile_values
        ),
        "error": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("testdata/smiles/corpus/smiles_small.smi"),
        help="input SMILES file (default: testdata/smiles/corpus/smiles_small.smi)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("testdata/descriptors/expected/rdkit/smiles_small/molecular_descriptors.jsonl"),
        help="output JSONL path (default: testdata/descriptors/expected/rdkit/smiles_small/molecular_descriptors.jsonl)",
    )
    args = parser.parse_args()

    smiles_rows = list(iter_smiles(args.input))
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w", encoding="utf-8") as f:
        for smiles in smiles_rows:
            f.write(json.dumps(build_record(smiles), ensure_ascii=True))
            f.write("\n")

    print(f"Wrote {len(smiles_rows)} records to {args.output}")


if __name__ == "__main__":
    main()
