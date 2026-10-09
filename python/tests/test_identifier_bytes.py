"""Local argument-extraction regressions, not an oracle-dependent corpus."""
import ast
from pathlib import Path

import cosmolkit as ck
import pytest


@pytest.mark.parametrize("input", ["CCO", b"CCO"])
def test_smiles_identifier_text_forms(input):
    assert ck.Molecule.from_smiles(input).to_smiles() == "CCO"
    assert ck.Molecule.from_smiles_with_params(input, ck.SmilesParseParams()).to_smiles() == "CCO"


@pytest.mark.parametrize("input", ["InChI=1S/CH4/h1H4", b"InChI=1S/CH4/h1H4"])
def test_inchi_identifier_text_forms(input):
    assert ck.Molecule.from_inchi(input).to_smiles() == "C"
    assert ck.Molecule.from_inchi_with_params(input, ck.InchiReadParams()).to_smiles() == "C"
    assert ck.Molecule.from_inchi(input, sanitize=True, remove_hs=True).to_smiles() == "C"
    assert ck.inchi_to_key(input) == "VNWKTOKETHGBQD-UHFFFAOYSA-N"


@pytest.mark.parametrize("call", [
    ck.Molecule.from_smiles,
    lambda text: ck.Molecule.from_smiles_with_params(text, ck.SmilesParseParams()),
    ck.Molecule.from_inchi,
    lambda text: ck.Molecule.from_inchi_with_params(text, ck.InchiReadParams()),
    ck.inchi_to_key,
])
def test_identifier_inputs_are_not_coerced_or_lossily_decoded(call):
    class TextLike:
        def __str__(self):
            raise AssertionError("arbitrary str(object) coercion must not run")

    for invalid in (None, 42, TextLike(), bytearray(b"CCO")):
        with pytest.raises(TypeError):
            call(invalid)
    with pytest.raises(UnicodeDecodeError) as error:
        call(b"CCO\xff")
    assert error.value.object == b"CCO\xff"


def test_generated_identifier_stubs_include_bytes():
    stub = Path(__file__).parents[1] / "cosmolkit.pyi"
    tree = ast.parse(stub.read_text())
    molecule = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "Molecule")
    functions = [n for n in molecule.body if isinstance(n, ast.FunctionDef) and n.name in (
        "from_smiles", "from_smiles_with_params", "from_inchi", "from_inchi_with_params")]
    functions += [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "inchi_to_key"]
    assert {f.name for f in functions} == {
        "from_smiles", "from_smiles_with_params", "from_inchi", "from_inchi_with_params", "inchi_to_key"}
    for function in functions:
        args = [arg for arg in function.args.args if arg.arg not in ("cls", "self")]
        annotation = ast.unparse(args[0].annotation)
        assert "str" in annotation and "bytes" in annotation, (function.name, annotation)
