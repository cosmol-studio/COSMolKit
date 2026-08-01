#!/usr/bin/env python3
"""Emit RDKit atom/bond feature lines for one SMILES (direct + AddHs)."""

from __future__ import annotations

import argparse
from dataclasses import dataclass

from rdkit import Chem


@dataclass(frozen=True)
class AtomLine:
    idx: int
    atomic_num: int
    degree: int
    formal_charge: int
    num_hs: int
    num_radical_electrons: int
    is_aromatic: bool
    is_in_ring: bool
    explicit_valence: int
    implicit_hs: int
    total_valence: int

    def as_line(self) -> str:
        return (
            f"A\t{self.idx}\t{self.atomic_num}\t{self.degree}\t{self.formal_charge}\t"
            f"{self.num_hs}\t{self.num_radical_electrons}\t"
            f"{int(self.is_aromatic)}\t{int(self.is_in_ring)}\t"
            f"{self.explicit_valence}\t{self.implicit_hs}\t{self.total_valence}"
        )


@dataclass(frozen=True)
class BondLine:
    begin_atom: int
    end_atom: int
    bond_type: str

    def normalized(self) -> "BondLine":
        if self.bond_type != "DATIVE" and self.begin_atom > self.end_atom:
            return BondLine(self.end_atom, self.begin_atom, self.bond_type)
        return self

    def as_line(self) -> str:
        return f"B\t{self.begin_atom}\t{self.end_atom}\t{self.bond_type}"


def atom_lines(mol: Chem.Mol) -> list[AtomLine]:
    rows: list[AtomLine] = []
    for atom in mol.GetAtoms():
        rows.append(
            AtomLine(
                idx=atom.GetIdx(),
                atomic_num=atom.GetAtomicNum(),
                degree=atom.GetTotalDegree(),
                formal_charge=atom.GetFormalCharge(),
                num_hs=atom.GetTotalNumHs(),
                num_radical_electrons=atom.GetNumRadicalElectrons(),
                is_aromatic=atom.GetIsAromatic(),
                is_in_ring=atom.IsInRing(),
                explicit_valence=int(atom.GetValence(Chem.rdchem.ValenceType.EXPLICIT)),
                implicit_hs=atom.GetNumImplicitHs(),
                total_valence=int(atom.GetValence(Chem.rdchem.ValenceType.EXPLICIT))
                + int(atom.GetValence(Chem.rdchem.ValenceType.IMPLICIT)),
            )
        )
    return rows


def bond_lines(mol: Chem.Mol) -> list[BondLine]:
    rows: list[BondLine] = []
    for bond in mol.GetBonds():
        rows.append(
            BondLine(
                begin_atom=bond.GetBeginAtomIdx(),
                end_atom=bond.GetEndAtomIdx(),
                bond_type=str(bond.GetBondType()),
            ).normalized()
        )
    rows.sort(key=lambda b: (b.begin_atom, b.end_atom, b.bond_type))
    return rows


def emit_section(name: str, mol: Chem.Mol) -> None:
    print(f"## {name}")
    for row in atom_lines(mol):
        print(row.as_line())
    for row in bond_lines(mol):
        print(row.as_line())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smiles", required=True, help="input SMILES")
    args = parser.parse_args()

    mol = Chem.MolFromSmiles(args.smiles, sanitize=True)
    if mol is None:
        print("ERROR\tMolFromSmiles returned None")
        raise SystemExit(2)

    emit_section("direct", mol)
    emit_section("with_hs", Chem.AddHs(Chem.Mol(mol)))


if __name__ == "__main__":
    main()
