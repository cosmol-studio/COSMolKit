"""Small required-configuration adapter regressions, independent of chemistry."""
import ast
import importlib.util
import inspect
from pathlib import Path
from types import SimpleNamespace

import pytest

_path = Path(__file__).resolve().parents[1] / "src/configuration_projection.py"
_spec = importlib.util.spec_from_file_location("_required_configuration_projection", _path)
assert _spec and _spec.loader
_projection = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_projection)


class Settings:
    def __init__(self, flag, scale):
        self.flag, self.scale = flag, scale


class Runner:
    def run(self, atom, params):
        return atom, params.flag, params.scale


def wrapped():
    fields = [{"name": "flag", "type": "builtins.bool", "default": None},
              {"name": "scale", "type": "builtins.float", "default": None}]
    row = {"semantic_id": "types.Settings", "python_name": "Settings", "rust_path": "crate::Settings", "role": "parameter", "python_fields": fields}
    target = {"semantic_id": "Runner.run", "python_name": "run", "parameters": [{"name": "atom", "type": "u32"}, {"name": "params", "type": "&crate::Settings"}]}
    runner = Runner()
    return _projection._wrapper(runner.run, runner.run, target, [(target["parameters"][1], row)], SimpleNamespace(Settings=Settings), [row, target])


def test_required_config_object_and_keywords_reach_one_actual_operation():
    run = wrapped()
    params = Settings(False, 3.0)
    assert run(7, params) == run(7, flag=False, scale=3.0) == (7, False, 3.0)
    params.scale = 4.0
    assert run(7, params) == (7, False, 4.0)


@pytest.mark.parametrize("keywords", [{}, {"flag": False}, {"scale": 3.0}, {"flag": False, "scale": 3.0, "unknown": 1}])
def test_required_config_missing_and_unknown_fields_are_rejected(keywords):
    with pytest.raises(TypeError):
        wrapped()(7, **keywords)


def test_required_config_mixed_fields_and_instance_are_rejected():
    with pytest.raises(TypeError, match="mutually exclusive"):
        wrapped()(7, Settings(False, 3.0), scale=2.0)


@pytest.mark.parametrize("defaults", [(None, None), ("False", None), ("False", "3.0")])
@pytest.mark.parametrize("native_short", [False, True])
def test_generated_call_forms_require_configuration_without_complete_defaults(defaults, native_short):
    fields = [
        {"name": "flag", "type": "builtins.bool", "default": defaults[0]},
        {"name": "scale", "type": "builtins.float", "default": defaults[1]},
    ]
    row = {"item": "type", "semantic_id": "types.Settings", "python_name": "Settings",
           "rust_path": "crate::Settings", "role": "parameter", "python_fields": fields}
    target = {"item": "callable", "owner": "module", "semantic_id": "run", "python_name": "run",
              "parameters": [{"name": "atom", "type": "u32"}, {"name": "params", "type": "&crate::Settings"}]}
    def run(atom, params):
        return atom, params

    def short(atom):
        return atom, Settings(False, float(atom))

    source = target
    if native_short:
        source = {**target, "parameters": target["parameters"][:1]}
        target = {**target, "semantic_id": "run_with_params", "python_name": "run_with_params"}
    document = {"entries": [row, source] + ([target] if native_short else []),
                "keywords": [{"semantic_id": "run", "target": target["semantic_id"]}]}
    module = SimpleNamespace(Settings=Settings)
    module.run = _projection._wrapper(short if native_short else run, run, target, [(target["parameters"][1], row)], module, document["entries"])
    if native_short:
        module.run_with_params = run
    stub = """
import builtins
import typing
class Settings:
    flag: builtins.bool
    scale: builtins.float
def run(atom: builtins.int, params: Settings) -> None: ...
"""
    if native_short:
        stub = stub.replace("def run(", "def run_with_params(") + "\ndef run(atom: builtins.int) -> None: ...\n"
    generated = _projection.declarations(module, stub, document)
    namespace = {}
    exec(generated, namespace)
    signatures = []
    for form in ast.parse(generated).body:
        if isinstance(form, ast.FunctionDef) and form.name == "run":
            form.decorator_list = []
            exec(compile(ast.Module(body=[form], type_ignores=[]), "<generated-stub>", "exec"), namespace)
            signatures.append(inspect.signature(namespace["run"]))

    def accepts(*args, **kwargs):
        for signature in signatures:
            try:
                signature.bind(*args, **kwargs)
            except TypeError:
                continue
            return True
        return False

    assert accepts(7) == (native_short or all(default is not None for default in defaults))
    assert accepts(7, flag=False) == (defaults[1] is not None)
    assert accepts(7, Settings(False, 3.0))
    assert accepts(7, flag=False, scale=3.0)
    assert not accepts(7, Settings(False, 3.0), scale=3.0)
