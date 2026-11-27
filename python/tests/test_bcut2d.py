"""BCUT2D values and structured failures on the packaged binding."""

import cosmolkit as ck
import pytest


def test_bcut_eigenvalues_and_unsupported_bond() -> None:
    ethanol = ck.Molecule.from_smiles("CCO")
    expected = [
        16.248991329327914, 10.908168035086371,
        1.3402935186707783, -1.52235869276019,
        1.2962426382141046, -1.5683635470355741,
        3.8358023355894861, 0.31840561317391147,
    ]
    assert ethanol.bcut2d() == pytest.approx(expected, rel=1e-12, abs=1e-12)

    with pytest.raises(ck.DescriptorReadError) as exc:
        ck.Molecule.from_smiles("[2H]C([2H])N->B").bcut2d()
    cause = exc.value.__cause__
    assert isinstance(cause, ck.DescriptorError)
    assert cause.variant == ck.DescriptorErrorKind.UnsupportedBcutBond
    assert cause.bond == 3 and cause.order == ck.BondOrder.Dative
