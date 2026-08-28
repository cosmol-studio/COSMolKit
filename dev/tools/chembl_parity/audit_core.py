#!/usr/bin/env python3
# pyright: reportAttributeAccessIssue=false
"""Differential ChEMBL audit for COSMolKit against pinned RDKit."""

from __future__ import annotations

import argparse
import json
import os
import struct
import time
import warnings
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Callable, cast

import cosmolkit
from rdkit import Chem, DataStructs, RDLogger
from rdkit import rdBase
from rdkit.Chem import (
    Crippen,
    Descriptors,
    GraphDescriptors,
    MACCSkeys,
    QED,
    rdFingerprintGenerator,
)
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdCIPLabeler


RDLogger.DisableLog("rdApp.*")

EXPECTED_RDKIT_VERSION = "2026.03.1"


def error_value(error: BaseException) -> dict[str, str]:
    return {"type": type(error).__name__, "message": str(error)}


def validate_runtime_versions() -> None:
    if rdBase.rdkitVersion != EXPECTED_RDKIT_VERSION:
        raise RuntimeError(
            "RDKit version mismatch: expected "
            f"{EXPECTED_RDKIT_VERSION}, got {rdBase.rdkitVersion}"
        )
    expected_cosmolkit = os.environ.get("COSMOLKIT_EXPECTED_VERSION")
    if expected_cosmolkit and cosmolkit.__version__ != expected_cosmolkit:
        raise RuntimeError(
            "COSMolKit version mismatch: expected "
            f"{expected_cosmolkit}, got {cosmolkit.__version__}"
        )


def atomic_json(path: Path, value: Any) -> None:
    temporary = path.with_suffix(path.suffix + f".tmp-{os.getpid()}")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def f64_bits(value: float) -> str:
    return struct.pack(">d", float(value)).hex()


def rd_explicit_valence(atom: Chem.Atom) -> int:
    if hasattr(Chem, "ValenceType"):
        return int(atom.GetValence(Chem.ValenceType.EXPLICIT))
    return int(atom.GetExplicitValence())


def rd_implicit_valence(atom: Chem.Atom) -> int:
    if hasattr(Chem, "ValenceType"):
        return int(atom.GetValence(Chem.ValenceType.IMPLICIT))
    return int(atom.GetImplicitValence())


def rd_atom_state(atom: Chem.Atom) -> tuple[Any, ...]:
    explicit = rd_explicit_valence(atom)
    implicit = rd_implicit_valence(atom)
    return (
        atom.GetAtomicNum(),
        atom.GetFormalCharge(),
        str(atom.GetChiralTag()),
        atom.GetIsotope() or None,
        atom.GetIsAromatic(),
        atom.GetNumExplicitHs(),
        atom.GetDegree(),
        explicit,
        atom.GetNumImplicitHs(),
        atom.GetTotalNumHs(),
        explicit + implicit,
        atom.GetNumRadicalElectrons(),
    )


def ck_atom_state(atom: Any) -> tuple[Any, ...]:
    return (
        atom.atomic_num(),
        atom.formal_charge(),
        atom.chiral_tag_name(),
        atom.isotope(),
        atom.is_aromatic(),
        atom.explicit_hydrogens(),
        atom.degree(),
        atom.explicit_valence(),
        atom.implicit_hydrogens(),
        atom.total_num_hs(),
        atom.total_valence(),
        atom.num_radical_electrons(),
    )


def normalize_rd_stereo(value: str) -> str:
    return value.removeprefix("STEREO")


def rd_bond_state(bond: Chem.Bond) -> tuple[Any, ...]:
    return (
        bond.GetBeginAtomIdx(),
        bond.GetEndAtomIdx(),
        str(bond.GetBondType()),
        str(bond.GetBondDir()),
        normalize_rd_stereo(str(bond.GetStereo())),
        tuple(bond.GetStereoAtoms()),
        bond.GetIsAromatic(),
    )


def ck_bond_state(bond: Any) -> tuple[Any, ...]:
    return (
        bond.begin_atom_idx(),
        bond.end_atom_idx(),
        bond.bond_type_name(),
        bond.bond_dir_name(),
        bond.stereo_name(),
        tuple(bond.stereo_atoms()),
        bond.is_aromatic(),
    )


class Audit:
    def __init__(self, output: Path, max_examples: int) -> None:
        self.output = output
        self.max_examples = max_examples
        self.counts: Counter[str] = Counter()
        self.example_counts: Counter[str] = Counter()
        self.started = time.monotonic()
        output.parent.mkdir(parents=True, exist_ok=True)
        self.findings = output.with_suffix(".findings.jsonl").open("w", encoding="utf-8")

    def count(self, name: str, amount: int = 1) -> None:
        self.counts[name] += amount

    def mismatch(
        self,
        feature: str,
        record: dict[str, Any],
        rdkit: Any,
        cosmolkit_value: Any,
        detail: str | None = None,
    ) -> None:
        key = f"mismatch.{feature}"
        self.counts[key] += 1
        if self.example_counts[key] >= self.max_examples:
            return
        self.example_counts[key] += 1
        self.findings.write(
            json.dumps(
                {
                    "feature": feature,
                    "chembl_id": record["chembl_id"],
                    "row": record["row"],
                    "smiles": record["smiles"],
                    "rdkit": rdkit,
                    "cosmolkit": cosmolkit_value,
                    "detail": detail,
                },
                ensure_ascii=True,
                default=str,
            )
            + "\n"
        )
        self.findings.flush()

    def finish(self, processed: int, mode: str, shard: int, shards: int) -> None:
        self.findings.close()
        summary = {
            "mode": mode,
            "shard": shard,
            "shards": shards,
            "processed": processed,
            "elapsed_seconds": time.monotonic() - self.started,
            "counts": dict(sorted(self.counts.items())),
        }
        atomic_json(self.output, summary)
        print(json.dumps(summary, sort_keys=True))


def iter_records(path: Path, shard: int, shards: int, limit: int | None):
    selected = 0
    with path.open(encoding="utf-8") as source:
        for index, line in enumerate(source):
            if index % shards != shard:
                continue
            if limit is not None and selected >= limit:
                break
            selected += 1
            yield json.loads(line)


def compare_value(
    audit: Audit,
    feature: str,
    record: dict[str, Any],
    rd_value: Any,
    ck_value: Any,
) -> None:
    if rd_value != ck_value:
        audit.mismatch(feature, record, rd_value, ck_value)
    else:
        audit.count(f"match.{feature}")


def parse_pair(audit: Audit, record: dict[str, Any]):
    smiles = record["smiles"]
    rd_error: Any = None
    try:
        rd_mol = Chem.MolFromSmiles(smiles, sanitize=True)
    except Exception as error:
        rd_mol = None
        rd_error = {"type": type(error).__name__, "message": str(error)}
    ck_error: Any = None
    try:
        ck_mol = cosmolkit.Molecule.from_smiles(smiles)
    except Exception as error:
        ck_mol = None
        ck_error = {"type": type(error).__name__, "message": str(error)}
    if rd_mol is None and ck_mol is None:
        audit.count("parse.both_rejected")
        if rd_error is not None:
            audit.count("rdkit.parse_exception")
        else:
            audit.count("rdkit.parse_none")
        return None, None
    if rd_mol is None or ck_mol is None:
        audit.mismatch(
            "parse",
            record,
            {"status": "rejected", "error": rd_error}
            if rd_mol is None
            else {"status": "accepted"},
            {"status": "rejected", "error": ck_error}
            if ck_mol is None
            else {"status": "accepted"},
        )
        return rd_mol, ck_mol
    audit.count("rdkit.parse_ok")
    audit.count("cosmolkit.parse_ok")
    return rd_mol, ck_mol


def audit_core(audit: Audit, record: dict[str, Any]) -> None:
    rd_mol, ck_mol = parse_pair(audit, record)
    if rd_mol is None or ck_mol is None:
        return

    compare_value(audit, "graph.num_atoms", record, rd_mol.GetNumAtoms(), ck_mol.num_atoms())
    compare_value(audit, "graph.num_bonds", record, rd_mol.GetNumBonds(), ck_mol.num_bonds())

    rd_atoms = [rd_atom_state(atom) for atom in rd_mol.GetAtoms()]
    ck_atoms = [ck_atom_state(atom) for atom in ck_mol.atoms()]
    compare_value(audit, "graph.atom_state", record, rd_atoms, ck_atoms)
    rd_bonds = [rd_bond_state(bond) for bond in rd_mol.GetBonds()]
    ck_bonds = [ck_bond_state(bond) for bond in ck_mol.bonds()]
    compare_value(audit, "graph.bond_state", record, rd_bonds, ck_bonds)

    try:
        rd_smiles = Chem.MolToSmiles(rd_mol)
    except Exception as error:
        rd_smiles = f"ERROR: {type(error).__name__}: {error}"
    try:
        ck_smiles = ck_mol.to_smiles()
    except Exception as error:
        ck_smiles = f"ERROR: {type(error).__name__}: {error}"
    compare_value(audit, "smiles.default", record, rd_smiles, ck_smiles)

    archived_inchi = record.get("standard_inchi") or ""
    archived_key = record.get("standard_inchi_key") or ""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rd_inchi = cast(str, cast(object, Chem.MolToInchi(rd_mol)))
    except Exception as error:
        rd_inchi = f"ERROR: {type(error).__name__}: {error}"
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ck_inchi = ck_mol.to_inchi()
    except Exception as error:
        ck_inchi = f"ERROR: {type(error).__name__}: {error}"
    compare_value(audit, "inchi.mol_to_inchi", record, rd_inchi, ck_inchi)
    if archived_inchi:
        compare_value(audit, "external.chembl_inchi_vs_rdkit", record, archived_inchi, rd_inchi)

    if rd_inchi and not rd_inchi.startswith("ERROR:"):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                rd_key = Chem.InchiToInchiKey(rd_inchi)
        except Exception as error:
            rd_key = f"ERROR: {type(error).__name__}: {error}"
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
            ck_key = cosmolkit.inchi_to_key(rd_inchi)
        except Exception as error:
            ck_key = f"ERROR: {type(error).__name__}: {error}"
        compare_value(audit, "inchi.inchi_to_key", record, rd_key, ck_key)
        if archived_key:
            compare_value(audit, "external.chembl_key_vs_rdkit", record, archived_key, rd_key)

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                rd_from_inchi = Chem.MolFromInchi(rd_inchi, sanitize=True, removeHs=True)
            rd_roundtrip = (
                Chem.MolToSmiles(rd_from_inchi) if rd_from_inchi is not None else None
            )
        except Exception as error:
            rd_roundtrip = f"ERROR: {type(error).__name__}: {error}"
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ck_from_inchi = cosmolkit.Molecule.from_inchi(
                    rd_inchi, sanitize=True, remove_hs=True
                )
            ck_roundtrip = (
                ck_from_inchi.to_smiles() if ck_from_inchi is not None else None
            )
        except Exception as error:
            ck_roundtrip = f"ERROR: {type(error).__name__}: {error}"
        compare_value(audit, "inchi.mol_from_inchi", record, rd_roundtrip, ck_roundtrip)

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                rd_mol_key = Chem.MolToInchiKey(rd_mol)
        except Exception as error:
            rd_mol_key = f"ERROR: {type(error).__name__}: {error}"
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
            ck_mol_key = ck_mol.to_inchi_key()
        except Exception as error:
            ck_mol_key = f"ERROR: {type(error).__name__}: {error}"
        compare_value(audit, "inchi.mol_to_key", record, rd_mol_key, ck_mol_key)


def descriptor_values_rd(mol: Chem.Mol) -> dict[str, Any]:
    values = {
        "mol_wt": f64_bits(Descriptors.MolWt(mol)),
        "mol_wt_only_heavy": f64_bits(Descriptors.MolWt(mol, onlyHeavy=True)),
        "exact_mol_wt": f64_bits(Descriptors.ExactMolWt(mol)),
        "exact_mol_wt_only_heavy": f64_bits(
            Descriptors.ExactMolWt(mol, onlyHeavy=True)
        ),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "formula_isotopes": rdMolDescriptors.CalcMolFormula(mol, separateIsotopes=True),
        "formula_isotopes_abbreviated_h": rdMolDescriptors.CalcMolFormula(
            mol, separateIsotopes=True, abbreviateHIsotopes=True
        ),
        "formula_no_isotopes_abbreviated_h": rdMolDescriptors.CalcMolFormula(
            mol, separateIsotopes=False, abbreviateHIsotopes=True
        ),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "fraction_csp3": f64_bits(rdMolDescriptors.CalcFractionCSP3(mol)),
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "rot_default": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "rot_non_strict": rdMolDescriptors.CalcNumRotatableBonds(
            mol, rdMolDescriptors.NumRotatableBondsOptions.NonStrict
        ),
        "rot_strict": rdMolDescriptors.CalcNumRotatableBonds(
            mol, rdMolDescriptors.NumRotatableBondsOptions.Strict
        ),
        "rot_strict_linkages": rdMolDescriptors.CalcNumRotatableBonds(
            mol, rdMolDescriptors.NumRotatableBondsOptions.StrictLinkages
        ),
        "qed": f64_bits(QED.qed(mol)),
    }
    for include_hs in (False, True):
        for force in (False, True):
            suffix = f"include_hs_{int(include_hs)}.force_{int(force)}"
            branch_mol = Chem.Mol(mol)
            values[f"logp.{suffix}"] = f64_bits(
                Crippen.MolLogP(branch_mol, includeHs=include_hs, force=force)
            )
            values[f"mr.{suffix}"] = f64_bits(
                Crippen.MolMR(branch_mol, includeHs=include_hs, force=force)
            )
    for force in (False, True):
        for include_sandp in (False, True):
            suffix = f"force_{int(force)}.include_sandp_{int(include_sandp)}"
            branch_mol = Chem.Mol(mol)
            values[f"tpsa.{suffix}"] = f64_bits(
                rdMolDescriptors.CalcTPSA(
                    branch_mol, force=force, includeSandP=include_sandp
                )
            )
    return values


def descriptor_values_ck(mol: Any, branch_smiles: str | None = None) -> dict[str, Any]:
    values = {
        "mol_wt": f64_bits(cosmolkit.calc_mol_wt(mol)),
        "mol_wt_only_heavy": f64_bits(cosmolkit.calc_mol_wt(mol, only_heavy=True)),
        "exact_mol_wt": f64_bits(cosmolkit.calc_exact_mol_wt(mol)),
        "exact_mol_wt_only_heavy": f64_bits(
            cosmolkit.calc_exact_mol_wt(mol, only_heavy=True)
        ),
        "formula": cosmolkit.calc_mol_formula(mol),
        "formula_isotopes": cosmolkit.calc_mol_formula(mol, separate_isotopes=True),
        "formula_isotopes_abbreviated_h": cosmolkit.calc_mol_formula(
            mol, separate_isotopes=True, abbreviate_h_isotopes=True
        ),
        "formula_no_isotopes_abbreviated_h": cosmolkit.calc_mol_formula(
            mol, separate_isotopes=False, abbreviate_h_isotopes=True
        ),
        "hbd": cosmolkit.calc_num_hbd(mol),
        "hba": cosmolkit.calc_num_hba(mol),
        "fraction_csp3": f64_bits(cosmolkit.calc_fraction_csp3(mol)),
        "aromatic_rings": cosmolkit.calc_num_aromatic_rings(mol),
        "rot_default": cosmolkit.calc_num_rotatable_bonds(mol),
        "rot_non_strict": cosmolkit.calc_num_rotatable_bonds(mol, mode="non_strict"),
        "rot_strict": cosmolkit.calc_num_rotatable_bonds(mol, mode="strict"),
        "rot_strict_linkages": cosmolkit.calc_num_rotatable_bonds(mol, mode="strict_linkages"),
        "qed": f64_bits(cosmolkit.calc_qed(mol)),
    }
    for include_hs in (False, True):
        for force in (False, True):
            suffix = f"include_hs_{int(include_hs)}.force_{int(force)}"
            # The RDKit golden generator starts every Crippen option branch
            # from a fresh molecule. Recreate the same isolation on the CK
            # side; cache-order behavior is audited separately.
            branch_mol = (
                cosmolkit.Molecule.from_smiles(branch_smiles)
                if branch_smiles is not None
                else mol
            )
            logp, mr = cosmolkit.calc_crippen_descriptors(
                branch_mol, include_hs=include_hs, force=force
            )
            values[f"logp.{suffix}"] = f64_bits(logp)
            values[f"mr.{suffix}"] = f64_bits(mr)
    for force in (False, True):
        for include_sandp in (False, True):
            suffix = f"force_{int(force)}.include_sandp_{int(include_sandp)}"
            values[f"tpsa.{suffix}"] = f64_bits(
                cosmolkit.calc_tpsa(
                    mol, force=force, include_sandp=include_sandp
                )
            )
    return values


def fresh_rd_molecule(smiles: str) -> Chem.Mol:
    molecule = Chem.MolFromSmiles(smiles, sanitize=True)
    if molecule is None:
        raise RuntimeError("RDKit rejected a SMILES accepted by the paired parser")
    return molecule


def fresh_ck_molecule(smiles: str) -> Any:
    return cosmolkit.Molecule.from_smiles(smiles)


def high_feasibility_descriptor_values_rd(
    mol: Chem.Mol, smiles: str
) -> dict[str, Any]:
    alpha_contributions = [0.0] * mol.GetNumAtoms()
    alpha = rdMolDescriptors.CalcHallKierAlpha(mol, alpha_contributions)
    labute_contributions: dict[str, Any] = {}
    for include_hydrogens in (False, True):
        branch = fresh_rd_molecule(smiles)
        atom_contributions, hydrogen_contribution = (
            rdMolDescriptors._CalcLabuteASAContribs(branch, include_hydrogens)
        )
        labute_contributions[f"include_hydrogens_{include_hydrogens}"] = {
            "asa": rdMolDescriptors.CalcLabuteASA(
                fresh_rd_molecule(smiles), include_hydrogens, False
            ),
            "atom_contributions": list(atom_contributions),
            "hydrogen_contribution": hydrogen_contribution,
        }

    def chi_cache_profile(function: Callable[..., float]) -> dict[str, Any]:
        branch = fresh_rd_molecule(smiles)
        return {
            "cold": [function(branch, order, False) for order in range(7)],
            "warm": [function(branch, order, False) for order in range(7)],
            "forced": [function(branch, order, True) for order in range(7)],
        }

    labute_branch = fresh_rd_molecule(smiles)
    vsa_branch = fresh_rd_molecule(smiles)
    custom_bins = [-0.2, 0.0, 0.25, 0.25, 0.8]
    values: dict[str, Any] = {
        "descriptors": {
            "chi_0": float(GraphDescriptors.Chi0(mol)),
            "chi_1": float(GraphDescriptors.Chi1(mol)),
            "chi_nv_orders_0_6": [
                rdMolDescriptors.CalcChiNv(mol, order, False)
                for order in range(7)
            ],
            "chi_nn_orders_0_6": [
                rdMolDescriptors.CalcChiNn(mol, order, False)
                for order in range(7)
            ],
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
            "hall_kier_alpha": alpha,
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
            "num_aromatic_heterocycles": (
                rdMolDescriptors.CalcNumAromaticHeterocycles(mol)
            ),
            "num_aromatic_carbocycles": (
                rdMolDescriptors.CalcNumAromaticCarbocycles(mol)
            ),
            "num_aliphatic_heterocycles": (
                rdMolDescriptors.CalcNumAliphaticHeterocycles(mol)
            ),
            "num_aliphatic_carbocycles": (
                rdMolDescriptors.CalcNumAliphaticCarbocycles(mol)
            ),
            "num_saturated_heterocycles": (
                rdMolDescriptors.CalcNumSaturatedHeterocycles(mol)
            ),
            "num_saturated_carbocycles": (
                rdMolDescriptors.CalcNumSaturatedCarbocycles(mol)
            ),
            "num_spiro_atoms": rdMolDescriptors.CalcNumSpiroAtoms(mol),
            "num_bridgehead_atoms": rdMolDescriptors.CalcNumBridgeheadAtoms(mol),
            "num_atom_stereo_centers": (
                rdMolDescriptors.CalcNumAtomStereoCenters(mol)
            ),
            "num_unspecified_atom_stereo_centers": (
                rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol)
            ),
            "fraction_csp3": rdMolDescriptors.CalcFractionCSP3(mol),
            "mqns": list(rdMolDescriptors.MQNs_(mol, False)),
            "labute_asa_include_hydrogens_false": (
                rdMolDescriptors.CalcLabuteASA(mol, False, False)
            ),
            "labute_asa_include_hydrogens_true": (
                rdMolDescriptors.CalcLabuteASA(mol, True, True)
            ),
            "slogp_vsa": list(rdMolDescriptors.SlogP_VSA_(mol, [], False)),
            "smr_vsa": list(rdMolDescriptors.SMR_VSA_(mol, [], False)),
        },
        "contributions": {
            "hall_kier_alpha": {
                "value": alpha,
                "atom_contributions": alpha_contributions,
            },
            "labute_asa": labute_contributions,
        },
        "cache_profiles": {
            "chi_nv": chi_cache_profile(rdMolDescriptors.CalcChiNv),
            "chi_nn": chi_cache_profile(rdMolDescriptors.CalcChiNn),
            "labute_asa_sequence": [
                rdMolDescriptors.CalcLabuteASA(labute_branch, False, False),
                rdMolDescriptors.CalcLabuteASA(labute_branch, True, False),
                rdMolDescriptors.CalcLabuteASA(labute_branch, True, True),
                rdMolDescriptors.CalcLabuteASA(labute_branch, False, False),
            ],
            "slogp_vsa_default_cold": list(
                rdMolDescriptors.SlogP_VSA_(vsa_branch, [], False)
            ),
            "slogp_vsa_default_warm": list(
                rdMolDescriptors.SlogP_VSA_(vsa_branch, [], False)
            ),
            "slogp_vsa_custom_forced": list(
                rdMolDescriptors.SlogP_VSA_(vsa_branch, custom_bins, True)
            ),
            "smr_vsa_default_warm": list(
                rdMolDescriptors.SMR_VSA_(vsa_branch, [], False)
            ),
            "smr_vsa_custom_forced": list(
                rdMolDescriptors.SMR_VSA_(vsa_branch, custom_bins, True)
            ),
        },
    }
    return float_tree_bits(values)


def high_feasibility_descriptor_values_ck(mol: Any, smiles: str) -> dict[str, Any]:
    alpha, alpha_contributions = (
        cosmolkit.calc_hall_kier_alpha_with_contributions(mol)
    )
    labute_contributions: dict[str, Any] = {}
    for include_hydrogens in (False, True):
        asa, atom_contributions, hydrogen_contribution = (
            cosmolkit.calc_labute_asa_contributions(
                fresh_ck_molecule(smiles),
                include_hydrogens=include_hydrogens,
                force=False,
            )
        )
        labute_contributions[f"include_hydrogens_{include_hydrogens}"] = {
            "asa": asa,
            "atom_contributions": atom_contributions,
            "hydrogen_contribution": hydrogen_contribution,
        }

    def chi_cache_profile(function: Callable[..., float]) -> dict[str, Any]:
        branch = fresh_ck_molecule(smiles)
        return {
            "cold": [function(branch, order, force=False) for order in range(7)],
            "warm": [function(branch, order, force=False) for order in range(7)],
            "forced": [function(branch, order, force=True) for order in range(7)],
        }

    labute_branch = fresh_ck_molecule(smiles)
    vsa_branch = fresh_ck_molecule(smiles)
    custom_bins = [-0.2, 0.0, 0.25, 0.25, 0.8]
    values: dict[str, Any] = {
        "descriptors": {
            "chi_0": cosmolkit.calc_chi_0(mol),
            "chi_1": cosmolkit.calc_chi_1(mol),
            "chi_nv_orders_0_6": [
                cosmolkit.calc_chi_nv(mol, order, force=False)
                for order in range(7)
            ],
            "chi_nn_orders_0_6": [
                cosmolkit.calc_chi_nn(mol, order, force=False)
                for order in range(7)
            ],
            "chi_0v": cosmolkit.calc_chi_0v(mol, force=False),
            "chi_1v": cosmolkit.calc_chi_1v(mol, force=False),
            "chi_2v": cosmolkit.calc_chi_2v(mol, force=False),
            "chi_3v": cosmolkit.calc_chi_3v(mol, force=False),
            "chi_4v": cosmolkit.calc_chi_4v(mol, force=False),
            "chi_0n": cosmolkit.calc_chi_0n(mol, force=False),
            "chi_1n": cosmolkit.calc_chi_1n(mol, force=False),
            "chi_2n": cosmolkit.calc_chi_2n(mol, force=False),
            "chi_3n": cosmolkit.calc_chi_3n(mol, force=False),
            "chi_4n": cosmolkit.calc_chi_4n(mol, force=False),
            "hall_kier_alpha": alpha,
            "kappa_1": cosmolkit.calc_kappa_1(mol),
            "kappa_2": cosmolkit.calc_kappa_2(mol),
            "kappa_3": cosmolkit.calc_kappa_3(mol),
            "phi": cosmolkit.calc_phi(mol),
            "lipinski_hba": cosmolkit.calc_lipinski_hba(mol),
            "lipinski_hbd": cosmolkit.calc_lipinski_hbd(mol),
            "num_hba": cosmolkit.calc_num_hba(mol),
            "num_hbd": cosmolkit.calc_num_hbd(mol),
            "num_heteroatoms": cosmolkit.calc_num_heteroatoms(mol),
            "num_amide_bonds": cosmolkit.calc_num_amide_bonds(mol),
            "num_heavy_atoms": cosmolkit.calc_num_heavy_atoms(mol),
            "num_atoms": cosmolkit.calc_num_atoms(mol),
            "num_rings": cosmolkit.calc_num_rings(mol),
            "num_heterocycles": cosmolkit.calc_num_heterocycles(mol),
            "num_aromatic_rings": cosmolkit.calc_num_aromatic_rings(mol),
            "num_saturated_rings": cosmolkit.calc_num_saturated_rings(mol),
            "num_aliphatic_rings": cosmolkit.calc_num_aliphatic_rings(mol),
            "num_aromatic_heterocycles": (
                cosmolkit.calc_num_aromatic_heterocycles(mol)
            ),
            "num_aromatic_carbocycles": (
                cosmolkit.calc_num_aromatic_carbocycles(mol)
            ),
            "num_aliphatic_heterocycles": (
                cosmolkit.calc_num_aliphatic_heterocycles(mol)
            ),
            "num_aliphatic_carbocycles": (
                cosmolkit.calc_num_aliphatic_carbocycles(mol)
            ),
            "num_saturated_heterocycles": (
                cosmolkit.calc_num_saturated_heterocycles(mol)
            ),
            "num_saturated_carbocycles": (
                cosmolkit.calc_num_saturated_carbocycles(mol)
            ),
            "num_spiro_atoms": cosmolkit.calc_num_spiro_atoms(mol),
            "num_bridgehead_atoms": cosmolkit.calc_num_bridgehead_atoms(mol),
            "num_atom_stereo_centers": (
                cosmolkit.calc_num_atom_stereo_centers(mol)
            ),
            "num_unspecified_atom_stereo_centers": (
                cosmolkit.calc_num_unspecified_atom_stereo_centers(mol)
            ),
            "fraction_csp3": cosmolkit.calc_fraction_csp3(mol),
            "mqns": cosmolkit.calc_mqns(mol),
            "labute_asa_include_hydrogens_false": cosmolkit.calc_labute_asa(
                mol, include_hydrogens=False, force=False
            ),
            "labute_asa_include_hydrogens_true": cosmolkit.calc_labute_asa(
                mol, include_hydrogens=True, force=True
            ),
            "slogp_vsa": cosmolkit.calc_slogp_vsa(mol, force=False),
            "smr_vsa": cosmolkit.calc_smr_vsa(mol, force=False),
        },
        "contributions": {
            "hall_kier_alpha": {
                "value": alpha,
                "atom_contributions": alpha_contributions,
            },
            "labute_asa": labute_contributions,
        },
        "cache_profiles": {
            "chi_nv": chi_cache_profile(cosmolkit.calc_chi_nv),
            "chi_nn": chi_cache_profile(cosmolkit.calc_chi_nn),
            "labute_asa_sequence": [
                cosmolkit.calc_labute_asa(
                    labute_branch, include_hydrogens=False, force=False
                ),
                cosmolkit.calc_labute_asa(
                    labute_branch, include_hydrogens=True, force=False
                ),
                cosmolkit.calc_labute_asa(
                    labute_branch, include_hydrogens=True, force=True
                ),
                cosmolkit.calc_labute_asa(
                    labute_branch, include_hydrogens=False, force=False
                ),
            ],
            "slogp_vsa_default_cold": cosmolkit.calc_slogp_vsa(
                vsa_branch, force=False
            ),
            "slogp_vsa_default_warm": cosmolkit.calc_slogp_vsa(
                vsa_branch, force=False
            ),
            "slogp_vsa_custom_forced": cosmolkit.calc_slogp_vsa(
                vsa_branch, bins=custom_bins, force=True
            ),
            "smr_vsa_default_warm": cosmolkit.calc_smr_vsa(
                vsa_branch, force=False
            ),
            "smr_vsa_custom_forced": cosmolkit.calc_smr_vsa(
                vsa_branch, bins=custom_bins, force=True
            ),
        },
    }
    return float_tree_bits(values)


def float_tree_bits(value: Any) -> Any:
    if isinstance(value, float):
        return f64_bits(value)
    if isinstance(value, list):
        return [float_tree_bits(item) for item in value]
    if isinstance(value, dict):
        return {key: float_tree_bits(item) for key, item in value.items()}
    return value


def compare_observation_tree(
    audit: Audit,
    feature: str,
    record: dict[str, Any],
    rd_value: Any,
    ck_value: Any,
    *,
    location: str = "",
) -> None:
    if isinstance(rd_value, dict) and isinstance(ck_value, dict):
        if rd_value.keys() != ck_value.keys():
            audit.mismatch(
                f"{feature}.shape",
                record,
                sorted(rd_value),
                sorted(ck_value),
                location or None,
            )
            return
        for key in rd_value:
            compare_observation_tree(
                audit,
                f"{feature}.{key}",
                record,
                rd_value[key],
                ck_value[key],
                location=location,
            )
        return
    if isinstance(rd_value, list) and isinstance(ck_value, list):
        if len(rd_value) != len(ck_value):
            audit.mismatch(
                f"{feature}.shape",
                record,
                len(rd_value),
                len(ck_value),
                location or None,
            )
            return
        for index, (rd_item, ck_item) in enumerate(zip(rd_value, ck_value)):
            item_location = f"{location}[{index}]" if location else f"[{index}]"
            compare_observation_tree(
                audit,
                feature,
                record,
                rd_item,
                ck_item,
                location=item_location,
            )
        return
    if rd_value != ck_value:
        audit.mismatch(feature, record, rd_value, ck_value, location or None)
    else:
        audit.count(f"match.{feature}")


def audit_descriptors(audit: Audit, record: dict[str, Any]) -> None:
    rd_mol, ck_mol = parse_pair(audit, record)
    if rd_mol is None or ck_mol is None:
        return
    try:
        rd_values = descriptor_values_rd(rd_mol)
    except Exception as error:
        audit.mismatch("descriptors.rdkit_error", record, repr(error), None)
        return
    try:
        ck_values = descriptor_values_ck(ck_mol, record["smiles"])
    except Exception as error:
        audit.mismatch("descriptors.cosmolkit_error", record, rd_values, repr(error))
        return
    for field, expected in rd_values.items():
        compare_value(audit, f"descriptor.{field}", record, expected, ck_values[field])
    try:
        rd_high_feasibility = high_feasibility_descriptor_values_rd(
            rd_mol, record["smiles"]
        )
    except Exception as error:
        audit.mismatch(
            "descriptor.high_feasibility.rdkit_error", record, repr(error), None
        )
        return
    try:
        ck_high_feasibility = high_feasibility_descriptor_values_ck(
            ck_mol, record["smiles"]
        )
    except Exception as error:
        audit.mismatch(
            "descriptor.high_feasibility.cosmolkit_error",
            record,
            rd_high_feasibility,
            repr(error),
        )
        return
    compare_observation_tree(
        audit,
        "descriptor.high_feasibility",
        record,
        rd_high_feasibility,
        ck_high_feasibility,
    )


def morgan_generator_with_num_bits_per_feature(value: int) -> Any:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    generator.GetOptions().numBitsPerFeature = value
    return generator


MORGAN_BRANCHES: dict[str, tuple[Any, dict[str, Any]]] = {
    "default": (
        rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048),
        {},
    ),
    "chiral": (
        rdFingerprintGenerator.GetMorganGenerator(
            radius=2, fpSize=2048, includeChirality=True
        ),
        {"include_chirality": True},
    ),
    "radius3": (
        rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048),
        {"radius": 3},
    ),
    "no_bond_types": (
        rdFingerprintGenerator.GetMorganGenerator(
            radius=2, fpSize=1024, useBondTypes=False
        ),
        {"n_bits": 1024, "use_bond_types": False},
    ),
    "count_simulation": (
        rdFingerprintGenerator.GetMorganGenerator(
            radius=2,
            fpSize=2048,
            countSimulation=True,
            countBounds=[1, 2, 4, 8],
        ),
        {"count_simulation": True, "count_bounds": [1, 2, 4, 8]},
    ),
    "redundant_environments": (
        rdFingerprintGenerator.GetMorganGenerator(
            radius=2, fpSize=2048, includeRedundantEnvironments=True
        ),
        {"include_redundant_environments": True},
    ),
    "feature_invariants": (
        rdFingerprintGenerator.GetMorganGenerator(
            radius=2,
            fpSize=2048,
            atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen(),
        ),
        {"atom_invariants_generator": "feature"},
    ),
    "connectivity_without_ring_membership": (
        rdFingerprintGenerator.GetMorganGenerator(
            radius=2,
            fpSize=2048,
            atomInvariantsGenerator=rdFingerprintGenerator.GetMorganAtomInvGen(
                includeRingMembership=False
            ),
        ),
        {
            "atom_invariants_generator": "connectivity",
            "atom_invariants_include_ring_membership": False,
        },
    ),
    "chiral_bond_invariants": (
        rdFingerprintGenerator.GetMorganGenerator(
            radius=2,
            fpSize=2048,
            bondInvariantsGenerator=rdFingerprintGenerator.GetMorganBondInvGen(
                useBondTypes=True, useChirality=True
            ),
        ),
        {
            "bond_invariants_generator": "morgan",
            "bond_invariants_use_chirality": True,
        },
    ),
    "two_bits_per_feature": (
        morgan_generator_with_num_bits_per_feature(2),
        {"num_bits_per_feature": 2},
    ),
}


def rdkit_additional_output_record(output: Any) -> dict[str, Any]:
    return {
        "atom_counts": list(output.GetAtomCounts()),
        "atom_to_bits": [list(bits) for bits in output.GetAtomToBits()],
        "bit_info_map": {
            int(bit): [tuple(pair) for pair in pairs]
            for bit, pairs in output.GetBitInfoMap().items()
        },
        "atoms_per_bit": {
            int(bit): [list(atoms) for atoms in atoms_per_bit]
            for bit, atoms_per_bit in output.GetAtomsPerBit().items()
        },
    }


def cosmolkit_additional_output_record(output: Any) -> dict[str, Any]:
    return {
        "atom_counts": output.atom_counts(),
        "atom_to_bits": output.atom_to_bits(),
        "bit_info_map": output.bit_info_map(),
        "atoms_per_bit": output.atoms_per_bit(),
    }


def dynamic_morgan_branches(rd_mol: Chem.Mol) -> list[tuple[str, Any, dict[str, Any], dict[str, Any]]]:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    branches: list[tuple[str, Any, dict[str, Any], dict[str, Any]]] = []
    if rd_mol.GetNumAtoms():
        branches.extend(
            [
                ("from_first", generator, {"fromAtoms": [0]}, {"from_atoms": [0]}),
                ("ignore_first", generator, {"ignoreAtoms": [0]}, {"ignore_atoms": [0]}),
            ]
        )
        atom_invariants = [idx + 1 for idx in range(rd_mol.GetNumAtoms())]
        branches.append(
            (
                "custom_atom_invariants",
                generator,
                {"customAtomInvariants": atom_invariants},
                {"custom_atom_invariants": atom_invariants},
            )
        )
        sparse_atom_invariants = [
            0 if idx % 2 == 0 else idx + 1 for idx in range(rd_mol.GetNumAtoms())
        ]
        only_nonzero_generator = rdFingerprintGenerator.GetMorganGenerator(
            radius=2, fpSize=2048, onlyNonzeroInvariants=True
        )
        branches.append(
            (
                "only_nonzero_custom_invariants",
                only_nonzero_generator,
                {"customAtomInvariants": sparse_atom_invariants},
                {
                    "custom_atom_invariants": sparse_atom_invariants,
                    "only_nonzero_invariants": True,
                },
            )
        )
    if rd_mol.GetNumBonds():
        bond_invariants = [idx + 7 for idx in range(rd_mol.GetNumBonds())]
        branches.append(
            (
                "custom_bond_invariants",
                generator,
                {"customBondInvariants": bond_invariants},
                {"custom_bond_invariants": bond_invariants},
            )
        )
    return branches


def audit_fingerprints(audit: Audit, record: dict[str, Any]) -> None:
    rd_mol, ck_mol = parse_pair(audit, record)
    if rd_mol is None or ck_mol is None:
        return
    branches = [
        (name, generator, {}, ck_kwargs)
        for name, (generator, ck_kwargs) in MORGAN_BRANCHES.items()
    ]
    branches.extend(dynamic_morgan_branches(rd_mol))
    for name, generator, rd_kwargs, ck_kwargs in branches:
        try:
            rd_output = rdFingerprintGenerator.AdditionalOutput()
            rd_output.AllocateAtomCounts()
            rd_output.AllocateAtomToBits()
            rd_output.AllocateBitInfoMap()
            rd_output.AllocateAtomsPerBit()
            rd_bits = list(
                generator.GetFingerprint(
                    rd_mol, additionalOutput=rd_output, **rd_kwargs
                ).GetOnBits()
            )
            rd_additional = rdkit_additional_output_record(rd_output)
        except Exception as error:
            audit.mismatch(f"morgan.{name}.rdkit_error", record, repr(error), None)
            continue
        try:
            ck_result = ck_mol.fingerprint_morgan_with_output(**ck_kwargs)
            ck_bits = ck_result.fingerprint().on_bits()
            ck_additional = cosmolkit_additional_output_record(
                ck_result.additional_output()
            )
        except Exception as error:
            audit.mismatch(f"morgan.{name}.cosmolkit_error", record, rd_bits, repr(error))
            continue
        compare_value(audit, f"morgan.{name}", record, rd_bits, ck_bits)
        compare_value(
            audit,
            f"morgan.{name}.additional_output",
            record,
            rd_additional,
            ck_additional,
        )

    try:
        rd_maccs = [bit - 1 for bit in MACCSkeys.GenMACCSKeys(rd_mol).GetOnBits() if bit > 0]
        ck_maccs = ck_mol.maccs_fingerprint().on_bits()
        compare_value(audit, "maccs.default", record, rd_maccs, ck_maccs)
    except Exception as error:
        audit.mismatch("maccs.error", record, "success", repr(error))


def rd_cip_property(owner: Any, name: str, kind: str) -> Any:
    if not owner.HasProp(name):
        return None
    if kind == "string":
        return owner.GetProp(name)
    if kind == "unsigned":
        return int(owner.GetUnsignedProp(name))
    if kind == "vector":
        return [int(value) & ((1 << 32) - 1) for value in json.loads(owner.GetProp(name))]
    raise ValueError(f"unsupported CIP property kind: {kind}")


def rd_cip_state(mol: Chem.Mol) -> dict[str, Any]:
    return {
        # RDKit stores this property as an integer in failed-labeling paths
        # (0/1), while the normal accessor exposes it as a bool.  Normalize
        # the source-observable value before comparing it with COSMolKit's
        # boolean property view.
        "cip_computed": bool(
            mol.HasProp("_CIPComputed") and mol.GetBoolProp("_CIPComputed")
        ),
        "atoms": [
            {
                "code": rd_cip_property(atom, "_CIPCode", "string"),
                "neighbor_order": rd_cip_property(atom, "_CIPNeighborOrder", "vector"),
                "rank": rd_cip_property(atom, "_CIPRank", "unsigned"),
            }
            for atom in mol.GetAtoms()
        ],
        "bonds": [
            {
                "code": rd_cip_property(bond, "_CIPCode", "string"),
                "neighbor_order": rd_cip_property(bond, "_CIPNeighborOrder", "vector"),
            }
            for bond in mol.GetBonds()
        ],
    }


def ck_cip_state(mol: Any) -> dict[str, Any]:
    return {
        "cip_computed": mol.cip_computed(),
        "atoms": [
            {
                "code": atom.cip_descriptor(),
                "neighbor_order": atom.cip_neighbor_order(),
                "rank": atom.cip_rank(),
            }
            for atom in mol.atoms()
        ],
        "bonds": [
            {
                "code": bond.cip_descriptor(),
                "neighbor_order": bond.cip_neighbor_order(),
            }
            for bond in mol.bonds()
        ],
    }


def audit_ciplabeler(audit: Audit, record: dict[str, Any]) -> None:
    rd_mol, ck_mol = parse_pair(audit, record)
    if rd_mol is None or ck_mol is None:
        return
    if rd_mol.GetNumAtoms() > 80:
        audit.count("skip.ciplabeler.max_atoms")
        return
    atom_selection = next(
        ([atom.GetIdx()] for atom in rd_mol.GetAtoms() if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED),
        [],
    )
    bond_selection = next(
        ([bond.GetIdx()] for bond in rd_mol.GetBonds() if bond.GetStereo() != Chem.BondStereo.STEREONONE),
        [],
    )
    branches = (
        ("full", None, None, 0),
        ("selected_atom", atom_selection, None, 0),
        ("selected_bond", None, bond_selection, 0),
        ("empty_selection", [], [], 1),
    )
    for name, atoms, bonds, limit in branches:
        rd_call = Chem.Mol(rd_mol)
        try:
            rdCIPLabeler.AssignCIPLabels(rd_call, atoms, bonds, limit)
            rd_value: Any = {"status": "ok", "state": rd_cip_state(rd_call)}
        except Exception as error:
            rd_value = {"status": "error", "error_message": str(error), "state": rd_cip_state(rd_call)}
        ck_call = ck_mol
        try:
            ck_call = type(ck_mol).from_smiles(record["smiles"])
            ck_call.assign_cip_labels_(atoms=atoms, bonds=bonds, max_recursive_iterations=limit)
            ck_value: Any = {"status": "ok", "state": ck_cip_state(ck_call)}
        except Exception as error:
            ck_value = {
                "status": "error",
                "error_message": str(error).split(": modern CIP labeling error: ", 1)[-1],
                "state": ck_cip_state(ck_call),
            }
        compare_value(audit, f"ciplabeler.{name}.state", record, rd_value, ck_value)


MODES: dict[str, Callable[[Audit, dict[str, Any]], None]] = {
    "core": audit_core,
    "descriptors": audit_descriptors,
    "fingerprints": audit_fingerprints,
    "ciplabeler": audit_ciplabeler,
}


def main() -> None:
    validate_runtime_versions()

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mode", choices=sorted(MODES), required=True)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--max-examples", type=int, default=100)
    parser.add_argument("--max-atoms", type=int)
    args = parser.parse_args()
    if not 0 <= args.shard < args.shards:
        raise SystemExit("shard must satisfy 0 <= shard < shards")

    audit = Audit(args.output, args.max_examples)
    processed = 0
    try:
        for processed, record in enumerate(
            iter_records(args.input, args.shard, args.shards, args.limit), start=1
        ):
            if args.max_atoms is not None:
                try:
                    probe = Chem.MolFromSmiles(record["smiles"], sanitize=True)
                except Exception:
                    probe = None
                if probe is not None and probe.GetNumAtoms() > args.max_atoms:
                    audit.count("skip.max_atoms")
                    continue
            MODES[args.mode](audit, record)
            if processed % 1000 == 0:
                print(
                    json.dumps(
                        {
                            "mode": args.mode,
                            "shard": args.shard,
                            "processed": processed,
                            "mismatches": sum(
                                value
                                for key, value in audit.counts.items()
                                if key.startswith("mismatch.")
                            ),
                        }
                    ),
                    flush=True,
                )
    finally:
        audit.finish(processed, args.mode, args.shard, args.shards)


if __name__ == "__main__":
    main()
