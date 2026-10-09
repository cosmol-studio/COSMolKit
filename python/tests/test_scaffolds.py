"""Small binding regressions; bulk parity is owned by the parity test crate."""
import cosmolkit as ck
import pytest


@pytest.mark.parametrize(
    "text,expected",
    [
        ("", ("", "", "")),
        ("CCO", ("", "", "")),
        ("Cc1ccccc1", ("c1ccccc1", "*c1ccccc1", "c1ccccc1")),
        ("O=C1CCCCC1", ("C1CCCCC1", "*=C1CCCCC1", "O=C1CCCCC1")),
        ("Cn1cccc1", ("c1cc[nH]c1", "*n1cccc1", "c1cc[nH]c1")),
        ("C[C@H]1CCCCO1", ("C1CCOCC1", "*[C@H]1CCCCO1", "C1CCOCC1")),
    ],
)
def test_scaffold_value_bindings(text: str, expected: tuple[str, str, str]) -> None:
    mol = ck.Molecule.from_smiles(text)
    before = mol.to_binary()
    for transform, want in zip(
        (mol.murcko_scaffold, mol.net_scaffold, mol.murcko_decompose), expected, strict=True
    ):
        output = transform()
        assert isinstance(output, ck.Molecule)
        assert output.to_smiles() == want
        assert mol.to_binary() == before


def test_scaffold_raw_input_error_is_atomic() -> None:
    mol = ck.Molecule.from_smiles_with_params(
        "CCO", ck.SmilesParseParams(sanitize=False, remove_hs=False)
    )
    before = mol.to_binary()
    for transform in (mol.murcko_scaffold, mol.net_scaffold, mol.murcko_decompose):
        with pytest.raises(ck.OperationError) as caught:
            _ = transform()
        assert caught.value.kind == "Scaffold"
        assert caught.value.__cause__ is not None
        assert mol.to_binary() == before
