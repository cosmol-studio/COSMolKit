"""Proposed source fixed cases and complete error/stub projection checks."""
import ast
from pathlib import Path
import pytest
import cosmolkit as ck

def fields(value):
    result = [value.symbol, value.branch_count, value.pi_electrons]
    if value.chirality is not None:
        result.append(value.chirality)
    return result

@pytest.mark.parametrize("code,expected", [(41, ["C", 1, 1]), (42, ["C", 2, 1]), (43, ["C", 3, 1]), (105, ["O", 1, 1]), (97, ["O", 1, 0])])
def test_source_documented_codes(code, expected):
    assert fields(ck.AtomCodeExplanation.from_code(code)) == expected

def test_optional_chirality_unused_subtract_full_source_table_and_high_bits():
    for code, label in [(35, ""), (547, "R"), (1059, "S")]:
        assert fields(ck.AtomCodeExplanation.from_code(code, include_chirality=True)) == ["C", 3, 0, label]
        assert fields(ck.AtomCodeExplanation.from_code(code)) == ["C", 3, 0]
    for subtract in [-(1 << 63), -2, 0, 1, 2, (1 << 63) - 1]:
        assert fields(ck.AtomCodeExplanation.from_code(481, subtract, True)) == ["*", 1, 0, ""]
    assert fields(ck.AtomCodeExplanation.from_code(547 | (1 << 63), 0, True)) == ["C", 3, 0, "R"]
    assert fields(ck.AtomCodeExplanation.from_code((1 << 64) - 1)) == ["*", 7, 3]

def test_source_key_error_three_retains_class_arguments_and_fields():
    for code in [1536, (1 << 64) - 1]:
        with pytest.raises(ck.AtomCodeExplanationError) as raised:
            ck.AtomCodeExplanation.from_code(code, include_chirality=True)
        error = raised.value
        assert isinstance(error, KeyError)
        assert error.args == (3,) and str(error) == "3"
        assert (error.domain, error.kind, error.code) == ("Fingerprint", "UnknownChirality", 3)
        assert ck.AtomCodeExplanation.from_code(code).chirality is None

def test_owned_frozen_value_and_conversion_failures():
    value = ck.AtomCodeExplanation.from_code(41)
    for name in ["symbol", "branch_count", "pi_electrons", "chirality", "other"]:
        with pytest.raises(AttributeError): setattr(value, name, None)
    for code in [-1, 1 << 64]:
        with pytest.raises(OverflowError): ck.AtomCodeExplanation.from_code(code)
    with pytest.raises(TypeError): ck.AtomCodeExplanation.from_code("41")
    with pytest.raises(TypeError): ck.AtomCodeExplanation(41)
    with pytest.raises(OverflowError): ck.AtomCodeExplanation.from_code(41, 1 << 63)

def test_original_generated_stub_has_complete_error_and_value_protocols():
    tree = ast.parse((Path(__file__).resolve().parents[1] / "cosmolkit.pyi").read_text())
    classes = {node.name: node for node in tree.body if isinstance(node, ast.ClassDef)}
    for name, base in [("AtomPairReadError", "builtins.ValueError"), ("TopologicalTorsionReadError", "builtins.ValueError"), ("AtomCodeExplanationError", "builtins.KeyError")]:
        declaration = classes[name]
        assert [ast.unparse(b) for b in declaration.bases] == [base]
        attributes = {node.target.id: ast.unparse(node.annotation) for node in declaration.body if isinstance(node, ast.AnnAssign)}
        assert attributes["domain"] == attributes["kind"] == "builtins.str"
        if name == "AtomCodeExplanationError": assert attributes["code"] == "builtins.int"
        assert getattr(ck, name).__name__ == name
    value = classes["AtomCodeExplanation"]
    methods = {node.name: node for node in value.body if isinstance(node, ast.FunctionDef)}
    assert set(methods) == {"from_code", "symbol", "branch_count", "pi_electrons", "chirality"}
    constructor = methods["from_code"]
    assert [a.arg for a in constructor.args.args] == ["code", "branch_subtract", "include_chirality"]
    assert [ast.literal_eval(v) for v in constructor.args.defaults] == [0, False]
    for name in ["symbol", "branch_count", "pi_electrons", "chirality"]:
        assert [ast.unparse(d) for d in methods[name].decorator_list] == ["property"]
