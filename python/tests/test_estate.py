"""Public atom EState and surface-area projections."""

import cosmolkit as ck
import pytest


def test_estate_bins_and_cache_lifetime() -> None:
    molecule = ck.Molecule.from_smiles("CCO")
    indices = molecule.estate_indices()
    area = molecule.estate_vsa()
    assert indices == pytest.approx([121 / 72, 0.25, 545 / 72])
    assert area[5] == pytest.approx(6.923737199690624)
    assert molecule.vsa_estate()[2] == indices[2]

    molecule.labute_asa(include_hydrogens=False, force=True)
    copy = molecule.copy()
    assert molecule.estate_vsa_6() == copy.estate_vsa_6() == area[5]
    assert sum(molecule.labute_asa_contributions().atom_contributions) == pytest.approx(
        18.68276248
    )
    assert sum(copy.labute_asa_contributions().atom_contributions) == pytest.approx(
        18.63714656
    )

    molecule.remove_atom_(1)
    assert molecule.estate_vsa_6() == 0.0
    assert indices == pytest.approx([121 / 72, 0.25, 545 / 72])
