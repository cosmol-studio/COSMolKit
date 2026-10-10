"""Typed potential-stereo analysis and lazy stereoisomer enumeration."""

from __future__ import annotations

import cosmolkit as ck


def main() -> None:
    source = ck.mol_from_smiles("CC(F)C(Cl)Br")
    source_smiles = source.to_smiles()

    analysis = source.potential_stereo()
    print(
        "potential centers:",
        [(item.stereo_type, item.centered_on) for item in analysis.stereo],
    )

    options = ck.StereoisomerOptions(
        max_isomers=4,
        random_source=ck.StereoisomerRandomSource.from_integer_seed(0xF00D),
    )
    print("upper-bound count:", source.stereoisomer_count(options))
    print("outputs:", [isomer.to_smiles() for isomer in source.enumerate_stereoisomers(options)])

    assert source.to_smiles() == source_smiles


if __name__ == "__main__":
    main()
