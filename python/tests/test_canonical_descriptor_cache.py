"""Canonical query lifecycle and custom boundary regressions from pinned native controls."""
import math
import pytest
import cosmolkit as ck

@pytest.mark.parametrize("name", ["slogp_vsa_with_params", "smr_vsa_with_params"])
def test_vsa_nan_boundary_follows_source_upper_bound(name):
    # Exact source-built RDKit release control on CCO, boundary [NaN].
    actual = getattr(ck.Molecule.from_smiles("CCO"), name)([math.nan], True)
    assert actual == [0.0, 18.637146559044247]


def test_query_cache_survives_coordinate_transform_and_clears_on_sanitize():
    # Source Crippen scalar memo is not keyed by includeHs. Native ROMol copy,
    # Compute2DCoords and SanitizeMol controls are retained in ROOT receipts.
    mol = ck.Molecule.from_smiles("CCO")
    no_h = mol.crippen_descriptors_with_params(False, False)
    assert (no_h.logp, no_h.molar_refractivity) == (-0.3487, 6.0798000000000005)
    depicted = mol.with_2d_coordinates()
    warm = depicted.crippen_descriptors_with_params(True, False)
    assert (warm.logp, warm.molar_refractivity) == (-0.3487, 6.0798000000000005)
    cleared = mol.sanitize()
    fresh = cleared.crippen_descriptors_with_params(True, False)
    assert (fresh.logp, fresh.molar_refractivity) == (-0.0014000000000000123, 12.759800000000002)
    original = mol.crippen_descriptors_with_params(True, False)
    assert (original.logp, original.molar_refractivity) == (-0.3487, 6.0798000000000005)
