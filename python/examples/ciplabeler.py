from __future__ import annotations

from cosmolkit import Molecule


def main() -> None:
    molecule = Molecule.from_smiles("C[C@H](F)Cl")
    labeled = molecule.with_cip_labels()
    print("full descriptor:", labeled.atoms()[1].cip_descriptor())

    selected = Molecule.from_smiles("C[C@H](F)Cl")
    selected.assign_cip_labels_(atoms=[1])
    print("selected descriptor:", selected.atoms()[1].cip_descriptor())

    double_bond = Molecule.from_smiles("F/C=C/F")
    double_bond.assign_cip_labels_()
    print("double-bond descriptor:", double_bond.bonds()[1].cip_descriptor())


if __name__ == "__main__":
    main()
