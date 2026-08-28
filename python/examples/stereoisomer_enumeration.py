"""Typed potential-stereo analysis and lazy stereoisomer enumeration."""

from __future__ import annotations

from cosmolkit import Molecule, StereoisomerOptions


def main() -> None:
    source = Molecule.from_smiles("CC(F)C(Cl)Br")
    source_smiles = source.to_smiles()

    analysis = source.analyze_potential_stereo()
    print(
        "potential centers:",
        [(item.center_kind, item.center_index) for item in analysis.stereo_info],
    )

    options = StereoisomerOptions(max_isomers=4, rand=0xF00D)
    print("upper-bound count:", source.stereoisomer_count(options))
    print("outputs:", [isomer.to_smiles() for isomer in source.stereoisomers(options)])

    assert source.to_smiles() == source_smiles


if __name__ == "__main__":
    main()
