from __future__ import annotations

import ast
import inspect
from collections import Counter
from pathlib import Path

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
            value = ast.literal_eval(node.value)
            assert isinstance(value, list)
            return set(value)
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
    duplicate_declarations = {
        name for name, count in Counter(declarations).items() if count != 1
    }
    assert not duplicate_declarations

    runtime_functions = {
        name
        for name, value in vars(cosmolkit).items()
        if not name.startswith("_") and inspect.isroutine(value)
    }
    stub_functions = set(declarations)
    exports = _stub_exports(module)

    assert runtime_functions == stub_functions
    assert runtime_functions <= exports
    assert "_rebuild_molecule_from_pickle" not in stub_functions


def test_assign_chiral_tags_methods_match_generated_stub_and_runtime_surface() -> None:
    molecule_class = _stub_class(_stub_module(), "Molecule")
    methods = [
        node
        for node in molecule_class.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name
        in {
            "with_chiral_tags_from_structure",
            "assign_chiral_tags_from_structure_",
        }
    ]

    assert Counter(method.name for method in methods) == {
        "with_chiral_tags_from_structure": 1,
        "assign_chiral_tags_from_structure_": 1,
    }
    for method in methods:
        assert [argument.arg for argument in method.args.args] == [
            "self",
            "conf_id",
            "replace_existing_tags",
        ]
        assert len(method.args.defaults) == 2
        assert all(
            isinstance(default, ast.Constant) and default.value is Ellipsis
            for default in method.args.defaults
        )

    assert str(inspect.signature(cosmolkit.Molecule.with_chiral_tags_from_structure)) == (
        "(self, /, conf_id=-1, replace_existing_tags=True)"
    )
    assert str(inspect.signature(cosmolkit.Molecule.assign_chiral_tags_from_structure_)) == (
        "(self, /, conf_id=-1, replace_existing_tags=True)"
    )
