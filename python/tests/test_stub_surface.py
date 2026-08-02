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
