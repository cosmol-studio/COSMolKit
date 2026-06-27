"""Generate RDKit golden data for ConfSeq's base ETKDG template step.

This intentionally covers the source-backed prefix of ConfSeq decoding:
TD chiral-token recovery -> get_p_chiral_mol_3d(..., is_op=False).
It does not depend on ConfSeq's Indigo-based bond-token mapping.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem


CASES = [
    {
        "case_id": "nonplanar_ring_aa2ar_l003582_candidate4",
        "in_smiles": "C c 1 c c c ( F ) c ( - C ^ | ( = O ) - N 2 - C - [ C @ @ H ] ( - C ^ | ( = O ) - N ^ | - [ C @ @ H ] 3 - C - C - C - [ C @ H ] ( C ) - [ C @ @ H ] - 3 C ) - ! C 3 ( - C - C - N - C - C - 3 ) - C - 2 ) c 1 F",
        "td_smiles": "C c 1 c c c ( F ) c ( <91> C <112> | ( = O ) <175> N 2 <-178> C <-100> [ C @ @ H ] ( <74> C <117> | ( = O ) <-177> N <119> | <-130> [ C @ @ H ] 3 <177> C <58> C <-54> C <56> [ C @ H ] ( C ) <-57> [ C @ @ H ] <179> 3 C ) <101> } C 3 ( <-171> C <54> C <-53> N <52> C <-53> C <178> 3 ) <25> C <160> 2 ) c 1 F",
    }
]


try:
    from indigo import Indigo
except ImportError:  # pragma: no cover - exercised only in incomplete dev envs.
    Indigo = None


def strip_confseq_topology_markers(smiles: str) -> str:
    stripped: list[str] = []
    in_bracket = False
    for token in smiles.split():
        if token == "[" and not in_bracket:
            in_bracket = True
            stripped.append(token)
        elif token == "]" and in_bracket:
            in_bracket = False
            stripped.append(token)
        elif token in {"^", "|", "!"} and not in_bracket:
            continue
        elif token == "-" and not in_bracket:
            continue
        else:
            stripped.append(token)
    return "".join(stripped)


def bond_token_atom_pairs(smiles: str) -> dict[int, tuple[int, int]]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"failed to parse stripped SMILES for bond mapping: {smiles}")
    smiles_be = Chem.MolToSmiles(mol, canonical=False, allBondsExplicit=True)
    if Indigo is None:
        raise RuntimeError(
            "gen_rdkit_confseq_embed_golden.py requires `epam.indigo`; "
            "ConfSeq maps tokens through Indigo molfile bond order"
        )
    mol_block = Indigo().loadMolecule(smiles_be).molfile()
    lines = mol_block.splitlines()
    num_atoms = int(lines[3][:3])
    num_bonds = int(lines[3][3:6])
    bond_lines = lines[4 + num_atoms : 4 + num_atoms + num_bonds]
    atom_pairs = [
        tuple(sorted((int(line[:3]) - 1, int(line[3:6]) - 1))) for line in bond_lines
    ]

    token_idx_to_atom_pair: dict[int, tuple[int, int]] = {}
    count = 0
    in_parentheses = 0
    for idx, ch in enumerate(smiles_be):
        if ch == "[":
            in_parentheses = 1
        elif ch == "]":
            in_parentheses = 0
        if ch in {"-", "=", "#", ":", "/", "\\"} and in_parentheses == 0:
            token_idx_to_atom_pair[idx] = atom_pairs[count]
            count += 1
    return token_idx_to_atom_pair


def complete_t_smiles(t_smiles: list[str], smiles_be: list[str]) -> list[str]:
    idx = 0
    while idx < len(smiles_be):
        if idx >= len(t_smiles):
            t_smiles.append("")
        elif smiles_be[idx] != t_smiles[idx] and ">" not in t_smiles[idx]:
            t_smiles.insert(idx, "-")
        idx += 1
    return t_smiles


def atom_pair_dihedral_literals(stripped_smiles: str, td_smiles: str) -> dict[tuple[int, int], str]:
    mapping = bond_token_atom_pairs(stripped_smiles)
    smiles_be = list(
        Chem.MolToSmiles(
            Chem.MolFromSmiles(stripped_smiles), canonical=False, allBondsExplicit=True
        )
    )
    completed = complete_t_smiles(td_smiles.split(), smiles_be)
    values: dict[tuple[int, int], str] = {}
    for token_idx, atom_pair in mapping.items():
        if token_idx < len(completed) and "<" in completed[token_idx]:
            values[atom_pair] = completed[token_idx]
    return values


def confseq_chiral_tokens(stripped_smiles: str, td_smiles: str) -> tuple[str, dict[int, str]]:
    td_smiles_clean = (
        td_smiles.replace("{ <", "{<")
        .replace("> {", ">{")
        .replace("} <", "}<")
        .replace("> }", ">}")
    )
    chiral: dict[int, str] = {}
    for atom_pair, raw in atom_pair_dihedral_literals(stripped_smiles, td_smiles_clean).items():
        if raw.startswith("{"):
            chiral[atom_pair[0]] = "CHI_TETRAHEDRAL_CW"
        elif raw.startswith("}"):
            chiral[atom_pair[0]] = "CHI_TETRAHEDRAL_CCW"
        if raw.endswith("{"):
            chiral[atom_pair[1]] = "CHI_TETRAHEDRAL_CW"
        elif raw.endswith("}"):
            chiral[atom_pair[1]] = "CHI_TETRAHEDRAL_CCW"
    return (
        td_smiles_clean.replace("{<", "<")
        .replace("}<", "<")
        .replace(">{", ">")
        .replace(">}", ">"),
        chiral,
    )


def confseq_embed_template(in_smiles: str, td_smiles: str) -> tuple[Chem.Mol, int]:
    for match in re.findall(r"<-?\d+>\s*\|", td_smiles):
        td_smiles = td_smiles.replace(" " + match, "")
    in_smiles = (
        in_smiles.replace("^ |", "")
        .replace(" !", "")
        .replace("/ -", "/")
        .replace("\\ -", "\\")
    )
    td_smiles = (
        td_smiles.replace("/ -", "/")
        .replace("\\ -", "\\")
        .replace("/ <", "/<")
        .replace("\\ <", "\\<")
        .replace("> /", ">/")
        .replace("> \\", ">\\")
    )
    smiles = strip_confseq_topology_markers(in_smiles)
    _, chiral = confseq_chiral_tokens(smiles, td_smiles)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("RDKit failed to parse ConfSeq stripped SMILES")

    original_charge: dict[int, int] = {}
    original_explicit_h: dict[int, int] = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        hyb = atom.GetHybridization()
        if symbol == "N" and hyb == Chem.rdchem.HybridizationType.SP3 and atom.GetFormalCharge() == 0:
            original_charge[atom.GetIdx()] = atom.GetFormalCharge()
            original_explicit_h[atom.GetIdx()] = atom.GetNumExplicitHs()
            atom.SetFormalCharge(1)
            atom.SetNumExplicitHs(1)
        if symbol == "S" and hyb == Chem.rdchem.HybridizationType.SP3 and atom.GetFormalCharge() == 1:
            original_charge[atom.GetIdx()] = atom.GetFormalCharge()
            original_explicit_h[atom.GetIdx()] = atom.GetNumExplicitHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(1)

    for idx, atom in enumerate(mol.GetAtoms(), start=1):
        atom.SetIsotope(idx)
    for atom_idx, tag_name in chiral.items():
        mol.GetAtomWithIdx(atom_idx).SetChiralTag(getattr(Chem.rdchem.ChiralType, tag_name))

    pre_add_hs_summary = molecule_summary(mol)
    mol_with_h = Chem.AddHs(mol)
    with_h_summary = molecule_summary(mol_with_h)
    params = AllChem.ETKDG()
    params.randomSeed = 0
    status = int(AllChem.EmbedMolecule(mol_with_h, params))
    if status < 0:
        mol_with_h.SetProp("_confseq_pre_add_hs_summary", json.dumps(pre_add_hs_summary))
        mol_with_h.SetProp("_confseq_with_h_summary", json.dumps(with_h_summary))
        return mol_with_h, status

    mol = Chem.RemoveHs(mol_with_h)
    for atom_idx, charge in original_charge.items():
        mol.GetAtomWithIdx(atom_idx).SetFormalCharge(charge)
    for atom_idx, explicit_h in original_explicit_h.items():
        mol.GetAtomWithIdx(atom_idx).SetNumExplicitHs(explicit_h)
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    mol.SetProp("_confseq_pre_add_hs_summary", json.dumps(pre_add_hs_summary))
    mol.SetProp("_confseq_with_h_summary", json.dumps(with_h_summary))
    return mol, status


def coords_for_mol(mol: Chem.Mol) -> list[list[float]]:
    conf = mol.GetConformer()
    return [[float(v) for v in conf.GetAtomPosition(i)] for i in range(mol.GetNumAtoms())]


def atom_summary(atom: Chem.Atom) -> dict[str, object]:
    return {
        "idx": int(atom.GetIdx()),
        "atomic_num": int(atom.GetAtomicNum()),
        "symbol": atom.GetSymbol(),
        "formal_charge": int(atom.GetFormalCharge()),
        "explicit_hs": int(atom.GetNumExplicitHs()),
        "isotope": int(atom.GetIsotope()),
        "chiral_tag": str(atom.GetChiralTag()),
        "hybridization": str(atom.GetHybridization()),
        "is_aromatic": bool(atom.GetIsAromatic()),
        "no_implicit": bool(atom.GetNoImplicit()),
    }


def bond_summary(bond: Chem.Bond) -> dict[str, object]:
    return {
        "idx": int(bond.GetIdx()),
        "begin": int(bond.GetBeginAtomIdx()),
        "end": int(bond.GetEndAtomIdx()),
        "order": str(bond.GetBondType()),
        "is_aromatic": bool(bond.GetIsAromatic()),
    }


def molecule_summary(mol: Chem.Mol) -> dict[str, object]:
    return {
        "num_atoms": int(mol.GetNumAtoms()),
        "num_bonds": int(mol.GetNumBonds()),
        "atoms": [atom_summary(atom) for atom in mol.GetAtoms()],
        "bonds": [bond_summary(bond) for bond in mol.GetBonds()],
    }


def build_record(case: dict[str, str]) -> dict[str, object]:
    record: dict[str, object] = {
        **case,
        "preset": "ETKDG",
        "random_seed": 0,
        "optimize_with_uff": False,
        "rdkit_ok": True,
        "status": None,
        "coords": None,
        "sdf": None,
        "pre_add_hs_summary": None,
        "with_h_summary": None,
        "error": None,
    }
    try:
        mol, status = confseq_embed_template(case["in_smiles"], case["td_smiles"])
        record["status"] = status
        record["pre_add_hs_summary"] = json.loads(mol.GetProp("_confseq_pre_add_hs_summary"))
        record["with_h_summary"] = json.loads(mol.GetProp("_confseq_with_h_summary"))
        if status >= 0:
            record["coords"] = coords_for_mol(mol)
            record["sdf"] = Chem.MolToMolBlock(mol) + "$$$$\n"
    except Exception as err:
        record["rdkit_ok"] = False
        record["error"] = str(err)
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/golden/confseq_embed_template.jsonl"),
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("tests/smiles.smi"),
        help="unused compatibility flag accepted for gen_all_rdkit_goldens.py",
    )
    parser.add_argument(
        "--sdf-output",
        type=Path,
        default=Path("tmp/confseq_nonplanar_compare/rdkit_confseq_embed_template.sdf"),
    )
    args = parser.parse_args()

    RDLogger.DisableLog("rdApp.*")
    records = [build_record(case) for case in CASES]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=True))
            handle.write("\n")

    args.sdf_output.parent.mkdir(parents=True, exist_ok=True)
    with args.sdf_output.open("w", encoding="utf-8") as handle:
        for record in records:
            if record["sdf"] is not None:
                handle.write(str(record["sdf"]))

    print(f"Wrote {len(records)} records to {args.output}")
    print(f"Wrote SDF to {args.sdf_output}")


if __name__ == "__main__":
    main()
