from __future__ import annotations

import ast
import inspect
from collections import Counter
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import cast

import cosmolkit


def _stub_module() -> ast.Module:
    path = Path(__file__).resolve().parents[1] / "cosmolkit.pyi"
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def _stub_exports(module: ast.Module) -> set[str]:
    for node in module.body:
        if not isinstance(node, ast.Assign):
            continue
        if any(
            isinstance(target, ast.Name) and target.id == "__all__"
            for target in node.targets
        ):
            value = cast(object, ast.literal_eval(node.value))
            assert isinstance(value, list)
            names: set[str] = set()
            for name in cast(list[object], value):
                assert isinstance(name, str)
                names.add(name)
            return names
    raise AssertionError("generated cosmolkit.pyi has no __all__ declaration")


def _stub_class(module: ast.Module, name: str) -> ast.ClassDef:
    matches = [
        node
        for node in module.body
        if isinstance(node, ast.ClassDef) and node.name == name
    ]
    assert len(matches) == 1
    return matches[0]


def test_generated_stub_covers_every_public_runtime_function_once() -> None:
    module = _stub_module()
    declarations = [
        node.name
        for node in module.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    ]
    for name, count in Counter(declarations).items():
        if count == 1:
            continue
        overloads = [node for node in module.body if isinstance(node, ast.FunctionDef) and node.name == name]
        assert count == 3  # default, parameter instance, keyword fields
        assert all("typing.overload" in [ast.unparse(value) for value in node.decorator_list] for node in overloads)
        assert hasattr(getattr(cosmolkit, name), "_configuration_contract")

    runtime_functions = {
        name
        for name, value in cast(Mapping[str, object], vars(cosmolkit)).items()
        if not name.startswith("_") and inspect.isroutine(value)
    }
    stub_functions = set(declarations)
    exports = _stub_exports(module)

    assert runtime_functions == stub_functions
    assert runtime_functions <= exports
    assert "_rebuild_molecule_from_pickle" not in stub_functions


def _assert_registered_method_surface(
    module: ast.Module, class_name: str, name: str, arguments: list[str], result: str
) -> None:
    declarations = [
        node for node in _stub_class(module, class_name).body
        if isinstance(node, ast.FunctionDef) and node.name == name
    ]
    if len(declarations) > 1:
        assert len(declarations) == 3
        assert all("typing.overload" in [ast.unparse(value) for value in node.decorator_list] for node in declarations)
    defaults = [node for node in declarations
                if [arg.arg for arg in node.args.posonlyargs + node.args.args] == arguments
                and not node.args.kwonlyargs]
    assert len(defaults) == 1
    method = defaults[0]
    assert [arg.arg for arg in method.args.posonlyargs + method.args.args] == arguments
    assert not method.args.defaults and not method.args.kwonlyargs
    assert method.returns is not None and ast.unparse(method.returns) == result
    owner = cast(type[object], getattr(cosmolkit, class_name))
    runtime = cast(Callable[..., object], getattr(owner, name))
    signature = inspect.signature(runtime)
    assert list(signature.parameters) == arguments
    assert all(cast(object, parameter.default) is inspect.Parameter.empty for parameter in signature.parameters.values())


def test_assign_chiral_tags_methods_match_generated_stub_and_runtime_surface() -> None:
    module = _stub_module()
    # Canonical bindings use StructureTagParams instead of positional flags.
    for name, arguments, result in [
        ("with_chiral_tags_from_structure", ["self"], "Molecule"),
        ("assign_chiral_tags_from_structure_", ["self"], "None"),
        ("with_chiral_tags_from_structure_with_params", ["self", "params"], "Molecule"),
        ("assign_chiral_tags_from_structure_with_params_", ["self", "params"], "None"),
    ]:
        _assert_registered_method_surface(module, "Molecule", name, arguments, result)


def test_layered_fingerprint_methods_match_generated_stub_and_runtime_surface() -> None:
    module = _stub_module()
    for class_name, name, arguments, result in [
        ("Molecule", "fingerprint_layered", ["self"], "Fingerprint"),
        ("Molecule", "fingerprint_layered_with_params", ["self", "params"], "Fingerprint"),
        ("Molecule", "fingerprint_layered_with_output", ["self"], "LayeredFingerprintResult"),
        ("Molecule", "fingerprint_layered_with_output_with_params", ["self", "params"], "LayeredFingerprintResult"),
        ("MoleculeBatch", "fingerprint_layered_list", ["self"], "builtins.list[typing.Optional[Fingerprint]]"),
        ("MoleculeBatch", "fingerprint_layered_list_with_params", ["self", "options", "params"], "builtins.list[typing.Optional[Fingerprint]]"),
        ("MoleculeBatch", "fingerprint_layered_with_output_list", ["self"], "builtins.list[typing.Optional[LayeredFingerprintResult]]"),
        ("MoleculeBatch", "fingerprint_layered_with_output_list_with_params", ["self", "options", "params"], "builtins.list[typing.Optional[LayeredFingerprintResult]]"),
    ]:
        _assert_registered_method_surface(module, class_name, name, arguments, result)

    _assert_registered_method_surface(module, "LayeredFingerprintResult", "fingerprint", ["self"], "Fingerprint")
    _assert_registered_method_surface(module, "LayeredFingerprintResult", "atom_counts", ["self"], "typing.Optional[builtins.list[builtins.int]]")


def test_pattern_fingerprint_methods_match_generated_stub_and_runtime_surface() -> None:
    module = _stub_module()
    for class_name, name, arguments, result in [
        ("Molecule", "fingerprint_pattern", ["self"], "Fingerprint"),
        ("Molecule", "fingerprint_pattern_with_params", ["self", "params"], "Fingerprint"),
        ("MoleculeBatch", "fingerprint_pattern_list", ["self"], "builtins.list[typing.Optional[Fingerprint]]"),
        ("MoleculeBatch", "fingerprint_pattern_list_with_params", ["self", "options", "params"], "builtins.list[typing.Optional[Fingerprint]]"),
    ]:
        _assert_registered_method_surface(module, class_name, name, arguments, result)


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
        )) == expected
        constructor_fn = cast(Callable[[], object], getattr(cosmolkit, class_name))
        params = constructor_fn()
        assert {name: cast(object, getattr(params, name)) for name in expected} == expected
