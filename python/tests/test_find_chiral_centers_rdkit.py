"""Language shape, defaults, errors and isolation of modern center results."""
import pytest
import cosmolkit as ck


@pytest.mark.parametrize("smiles, expected", [
    ("C", []),
    ("C[C@H](O)Cl", [(1, "R")]),
    ("C[C@@H](O)Cl", [(1, "S")]),
    ("CC(O)Cl", [(1, "?")]),
    ("C1C[C@H](C)[C@H](C)[C@H](C)C1", [(2, "S"), (4, "s"), (6, "R")]),
    ("F[Pt@SP1](Cl)(Br)I", []),
])
def test_chiral_center_tuple_results_and_default_filter(
    smiles: str, expected: list[tuple[int, str]],
) -> None:
    # Pinned RDKit 2026.03.1, modern perception and includeCIP=True.
    source = ck.Molecule.from_smiles(smiles)
    before = source.to_binary()
    result = source.find_chiral_centers(include_unassigned=True)
    assert result == expected
    assert isinstance(result, list)
    assert all(isinstance(row, tuple) and isinstance(row[0], int) and isinstance(row[1], str) for row in result)
    assigned = [row for row in expected if row[1] != "?"]
    assert source.find_chiral_centers() == assigned
    assert source.find_chiral_centers(include_unassigned=False) == assigned
    result.append((999, "modified"))
    assert source.find_chiral_centers(True) == expected
    assert source.to_binary() == before


def test_chiral_center_failure_retains_typed_cause_and_input_state():
    source = ck.Molecule.from_smiles_with_params(
        "C[C@H]1CCCC[C@H]1C |atomProp:1._ringStereochemCand.malformed|",
        ck.SmilesParseParams(sanitize=False, remove_hs=False, skip_cleanup=True),
    )
    before = source.to_binary()
    with pytest.raises(ck.StereoReadError) as caught:
        _ = source.find_chiral_centers()
    error = caught.value
    assert error.domain == "stereo" and error.kind == "PotentialStereo"
    cause = error.__cause__
    assert isinstance(cause, ck.PotentialStereoError)
    assert cause.kind == "InvalidPropertyKind"
    assert cause.atom == 1
    assert cause.property == "_ringStereochemCand"
    assert source.to_binary() == before


def test_chiral_center_preserves_existing_cip_property_text():
    source = ck.Molecule.from_smiles_with_params(
        "CC(O)Cl |atomProp:1._CIPCode.foo|",
        ck.SmilesParseParams(sanitize=False, remove_hs=False, skip_cleanup=True),
    )
    before = source.to_binary()
    assert source.find_chiral_centers(True) == [(1, "foo")]
    assert source.find_chiral_centers() == []
    assert source.to_binary() == before
