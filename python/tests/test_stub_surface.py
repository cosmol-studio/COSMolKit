from __future__ import annotations

import ast
import inspect
from collections import Counter
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import cast

import cosmolkit
import pytest


def _stub_module() -> ast.Module:
    path = Path(__file__).resolve().parents[1] / "artifacts/smoke/package/cosmolkit/__init__.pyi"
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _stub_exports(module: ast.Module) -> set[str]:
    names: set[str] = set()
    declared = False
    for node in module.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "__all__"
            for target in node.targets
        ):
            assert not declared, "duplicate __all__ declaration"
            declared = True
        elif isinstance(node, ast.AugAssign) and isinstance(node.target, ast.Name) and node.target.id == "__all__":
            assert declared and isinstance(node.op, ast.Add)
        else:
            continue
        value = cast(object, ast.literal_eval(node.value))
        assert isinstance(value, list)
        for name in cast(list[object], value):
            assert isinstance(name, str)
            names.add(name)
    assert declared, "generated cosmolkit.pyi has no __all__ declaration"
    return names


def _stub_class(module: ast.Module, name: str) -> ast.ClassDef:
    matches = [
        node
        for node in module.body
        if isinstance(node, ast.ClassDef) and node.name == name
    ]
    assert len(matches) == 1
    return matches[0]


def _stub_function_aliases(module: ast.Module) -> set[str]:
    aliases: set[str] = set()
    for node in module.body:
        if not isinstance(node, ast.Assign) or not isinstance(node.value, ast.Attribute):
            continue
        assert len(node.targets) == 1 and isinstance(node.targets[0], ast.Name)
        name = node.targets[0].id
        assert name not in aliases, f"duplicate callable alias: {name}"
        assert isinstance(node.value.value, ast.Name)
        owner_name, method_name = node.value.value.id, node.value.attr
        owner = _stub_class(module, owner_name)
        assert any(isinstance(member, ast.FunctionDef) and member.name == method_name for member in owner.body)
        original = getattr(getattr(cosmolkit, owner_name), method_name)
        assert inspect.isroutine(original)
        assert inspect.signature(getattr(cosmolkit, name)) == inspect.signature(original)
        aliases.add(name)
    return aliases


def test_stub_export_collection_includes_additions() -> None:
    module = ast.parse('__all__ = ["Molecule"]\n__all__ += ["mol_from_smiles"]\n')
    assert _stub_exports(module) == {"Molecule", "mol_from_smiles"}


def test_stub_alias_collection_rejects_duplicate_declarations() -> None:
    declaration = "mol_from_smiles = Molecule.from_smiles\n"
    module = ast.parse("class Molecule:\n    def from_smiles(text: str): ...\n" + declaration)
    assert _stub_function_aliases(module) == {"mol_from_smiles"}
    module.body.extend(ast.parse(declaration).body)
    with pytest.raises(AssertionError, match="duplicate callable alias"):
        _stub_function_aliases(module)


def test_generated_stub_covers_every_public_runtime_function_once() -> None:
    module = _stub_module()
    declarations = [
        node.name
        for node in module.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and not node.name.startswith("_")
    ]
    for name, count in Counter(declarations).items():
        if count == 1:
            continue
        overloads = [node for node in module.body if isinstance(node, ast.FunctionDef) and node.name == name]
        assert count == 2  # optional parameter instance or keyword fields
        assert all("typing.overload" in [ast.unparse(value) for value in node.decorator_list] for node in overloads)

    runtime_functions = {
        name
        for name, value in cast(Mapping[str, object], vars(cosmolkit)).items()
        if not name.startswith("_") and inspect.isroutine(value)
    }
    aliases = _stub_function_aliases(module)
    assert aliases.isdisjoint(declarations), "alias duplicates a function declaration"
    stub_functions = set(declarations) | aliases
    exports = _stub_exports(module)

    assert runtime_functions == stub_functions
    assert runtime_functions <= exports
    assert "_rebuild_molecule_from_pickle" not in stub_functions


def _assert_registered_method_surface(
    module: ast.Module, class_name: str, name: str, arguments: list[str], result: str,
    configuration: type | None = None,
) -> None:
    declarations = [
        node for node in _stub_class(module, class_name).body
        if isinstance(node, ast.FunctionDef) and node.name == name
    ]
    owner = getattr(cosmolkit, class_name)
    runtime = getattr(owner, name)
    if len(declarations) > 1:
        assert len(declarations) == 2
        assert all("typing.overload" in [ast.unparse(value) for value in node.decorator_list] for node in declarations)
        assert configuration is not None
        fields = list(inspect.signature(configuration).parameters)
        shared = ["execution", "n_jobs", "progress_bar"] if class_name == "MoleculeBatch" else []
        forms = {}
        for method in declarations:
            assert method.args.vararg is None and method.args.kwarg is None
            assert [arg.arg for arg in method.args.posonlyargs] == ["self"]
            positional = tuple(arg.arg for arg in method.args.args)
            assert positional not in forms
            forms[positional] = method
        assert set(forms) == {(), ("params",)}
        keywords, params = forms[()], forms[("params",)]
        assert not keywords.args.defaults
        assert len(params.args.defaults) == 1
        assert ast.literal_eval(params.args.defaults[0]) is None
        assert ast.unparse(params.args.args[0].annotation) == configuration.__name__ + " | None"
        assert [arg.arg for arg in keywords.args.kwonlyargs] == shared + fields
        assert [arg.arg for arg in params.args.kwonlyargs] == shared
        for method in declarations:
            assert all(default is not None and ast.literal_eval(default) in (None, Ellipsis)
                       for default in method.args.kw_defaults)
        signature = inspect.signature(runtime)
        assert list(signature.parameters) == ["self", "params", *shared, *fields]
        assert signature.parameters["self"].kind == inspect.Parameter.POSITIONAL_ONLY
        assert signature.parameters["params"].kind == inspect.Parameter.POSITIONAL_OR_KEYWORD
        assert signature.parameters["params"].default is Ellipsis
        for field in shared + fields:
            assert signature.parameters[field].kind == inspect.Parameter.KEYWORD_ONLY
            assert signature.parameters[field].default in (None, Ellipsis)
    else:
        assert len(declarations) == 1
        method = declarations[0]
        assert [arg.arg for arg in method.args.posonlyargs + method.args.args] == arguments
        assert not method.args.defaults and not method.args.kwonlyargs
        assert method.args.vararg is None and method.args.kwarg is None
        if inspect.isdatadescriptor(runtime):
            assert [ast.unparse(d) for d in method.decorator_list] == ["builtins.property"]
        else:
            signature = inspect.signature(runtime)
            assert list(signature.parameters) == arguments
            assert all(p.default is inspect.Parameter.empty for p in signature.parameters.values())
    for method in declarations:
        assert method.returns is not None and ast.unparse(method.returns) == result


def test_assign_chiral_tags_methods_match_generated_stub_and_runtime_surface() -> None:
    module = _stub_module()
    # Canonical bindings use StructureTagParams instead of positional flags.
    for name, arguments, result in [
        ("with_chiral_tags_from_structure", ["self"], "Molecule"),
        ("assign_chiral_tags_from_structure_", ["self"], "None"),
        ("with_chiral_tags_from_structure_with_params", ["self", "params"], "Molecule"),
        ("assign_chiral_tags_from_structure_with_params_", ["self", "params"], "None"),
    ]:
        _assert_registered_method_surface(module, "Molecule", name, arguments, result, cosmolkit.StructureTagParams)


def test_layered_fingerprint_methods_match_generated_stub_and_runtime_surface() -> None:
    module = _stub_module()
    for class_name, name, arguments, result in [
        ("Molecule", "fingerprint_layered", ["self"], "Fingerprint"),
        ("Molecule", "fingerprint_layered_with_params", ["self", "params"], "Fingerprint"),
        ("Molecule", "fingerprint_layered_with_output", ["self"], "LayeredFingerprintResult"),
        ("Molecule", "fingerprint_layered_with_output_with_params", ["self", "params"], "LayeredFingerprintResult"),
        ("MoleculeBatch", "fingerprint_layered", ["self"], "FingerprintBatch"),
        ("MoleculeBatch", "fingerprint_layered_with_params", ["self", "options", "execution"], "FingerprintBatch"),
        ("MoleculeBatch", "fingerprint_layered_with_output", ["self"], "list[LayeredFingerprintResult | None]"),
        ("MoleculeBatch", "fingerprint_layered_with_output_with_params", ["self", "options", "execution"], "list[LayeredFingerprintResult | None]"),
    ]:
        _assert_registered_method_surface(module, class_name, name, arguments, result, cosmolkit.LayeredFingerprintParams)

    _assert_registered_method_surface(module, "LayeredFingerprintResult", "fingerprint", ["self"], "Fingerprint")
    _assert_registered_method_surface(module, "LayeredFingerprintResult", "atom_counts", ["self"], "list[int] | None")


def test_pattern_fingerprint_methods_match_generated_stub_and_runtime_surface() -> None:
    module = _stub_module()
    for class_name, name, arguments, result in [
        ("Molecule", "fingerprint_pattern", ["self"], "Fingerprint"),
        ("Molecule", "fingerprint_pattern_with_params", ["self", "params"], "Fingerprint"),
        ("MoleculeBatch", "fingerprint_pattern", ["self"], "FingerprintBatch"),
        ("MoleculeBatch", "fingerprint_pattern_with_params", ["self", "options", "execution"], "FingerprintBatch"),
    ]:
        _assert_registered_method_surface(module, class_name, name, arguments, result, cosmolkit.PatternFingerprintParams)


def test_typed_parameters_preserve_source_defaults_in_runtime_and_stub() -> None:
    module = _stub_module()
    cases = [
        ("StructureTagParams", {"conformer_id": -1, "replace_existing_tags": True}),
        ("LayeredFingerprintParams", {
            "layers": 4294967295, "min_path": 1, "max_path": 7, "fp_size": 2048,
            "atom_counts": None, "set_only_bits": None, "branched_paths": True, "from_atoms": None,
        }),
        ("PatternFingerprintParams", {"n_bits": 2048, "tautomeric": False}),
    ]
    for class_name, expected in cases:
        constructor = next(node for node in _stub_class(module, class_name).body if isinstance(node, ast.FunctionDef) and node.name == "__new__")
        defaults: list[object] = []
        for value in constructor.args.kw_defaults:
            assert value is not None
            defaults.append(cast(object, ast.literal_eval(value)))
        assert dict(zip(
            [arg.arg for arg in constructor.args.kwonlyargs],
            defaults, strict=True,
        )) == dict.fromkeys(expected, Ellipsis)
        constructor_fn = cast(Callable[[], object], getattr(cosmolkit, class_name))
        params = constructor_fn()
        assert {name: cast(object, getattr(params, name)) for name in expected} == expected
