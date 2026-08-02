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

RDLogger.DisableLog("rdApp.*")


def iter_smiles(path: Path) -> Iterable[str]:
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        yield line


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
            "error": "MolFromSmiles returned None",
        }
    descriptors = descriptor_set(mol)
    descriptor_options, descriptor_option_bits = descriptor_option_matrix(smiles)
    return {
        "smiles": smiles,
        "rdkit_ok": True,
        "descriptors": descriptors,
        "descriptor_bits": descriptor_bits(descriptors),
        "descriptor_options": descriptor_options,
        "descriptor_option_bits": descriptor_option_bits,
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
