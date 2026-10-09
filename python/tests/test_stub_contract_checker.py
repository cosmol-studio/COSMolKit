"""Fixed regressions for the registry-driven stub-generation gate."""

from __future__ import annotations

import ast
import importlib.util
import json
from collections.abc import Callable
from pathlib import Path
from typing import cast

import pytest


_path = Path(__file__).resolve().parents[2] / "dev/tools/check_python_stub_contract.py"
_spec = importlib.util.spec_from_file_location("check_python_stub_contract", _path)
assert _spec is not None and _spec.loader is not None
_checker = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_checker)
check_contract = cast(Callable[[str, str], list[str]], _checker.check_contract)
check_runtime = _checker.check_runtime


@pytest.mark.parametrize("method", ["", "def __repr__(self): ...", "def __repr__(self) -> int: ...", "def __repr__(self, extra) -> str: ...", "def __repr__(self, *args) -> str: ...", "def __repr__(self, **kwargs) -> str: ..."])
def test_batch_configuration_gate_requires_typed_repr_declaration(method: str):
    row = dict(_entry("types.Settings", "Settings", "type", "type"), rust_path="crate::Settings", feature="cap-batch", role="parameter", fields=[])
    stub = "class Settings:\n    def __init__(self) -> None: ...\n"
    if method:
        stub += f"    {method}\n"
    errors = check_contract(stub, json.dumps({"entries": [row], "python_adapters": []}))
    assert len(errors) == 1
    assert "requires __repr__(self) -> str" in errors[0]


@pytest.mark.parametrize("feature,required", [("cap-batch", []), ("cap-io", ["cap-batch"])])
def test_batch_repr_requirement_follows_registered_capabilities_not_class_name(feature: str, required: list[str]):
    row = dict(_entry("types.Settings", "Settings", "type", "type"), rust_path="crate::Settings", feature=feature, required_capabilities=required, role="parameter", fields=[])
    stub = "class Settings:\n    def __init__(self) -> None: ...\n    def __repr__(self) -> builtins.str: ...\n"
    assert check_contract(stub, json.dumps({"entries": [row], "python_adapters": []})) == []


def test_batch_repr_runtime_gate_rejects_missing_incomplete_stale_and_mutating_repr():
    class Settings:
        def __init__(self, count=1):
            self.count = count

    row = {"python_name": "Settings", "fields": [{"name": "count"}]}

    def check(value):
        return _checker.check_configuration_repr(value, row, {Settings: row["fields"]})

    assert check(Settings())
    Settings.__repr__ = lambda self: "Settings()"
    assert check(Settings())
    Settings.__repr__ = lambda self: "Settings(count=1)"
    assert check(Settings()) == []
    assert check(Settings(2))

    def mutating_repr(self):
        self.count += 1
        return f"Settings(count={self.count})"

    Settings.__repr__ = mutating_repr
    assert any("changed configuration state" in error for error in check(Settings()))


@pytest.mark.parametrize("annotation", ["str", "str | bytes", "str | os.PathLike", "str | os.PathLike[bytes]", "str | os.PathLike[str] | object", "str | os.PathLike[str] | None"])
def test_sdf_path_gate_rejects_incomplete_or_overbroad_types(annotation: str):
    row = dict(_entry("MoleculeBatch.read_sdf", "read_sdf", "type"), parameters=[{"name": "path", "type": "&str", "default": None}])
    document = {"entries": [row], "python_adapters": []}
    stub = f"class MoleculeBatch:\n    def read_sdf(path: {annotation}) -> MoleculeBatch: ...\n"
    errors = check_contract(stub, json.dumps(document))
    assert len(errors) == 1
    assert "expected str | os.PathLike[str]" in errors[0]


@pytest.mark.parametrize("annotation", ["str | os.PathLike[str]", "typing.Union[builtins.str, os.PathLike[builtins.str]]", "str | os.PathLike[str] | pathlib.Path"])
def test_sdf_path_gate_accepts_only_text_filesystem_protocol(annotation: str):
    row = dict(_entry("MoleculeBatch.read_sdf", "read_sdf", "type"), parameters=[{"name": "path", "type": "&str", "default": None}])
    document = {"entries": [row], "python_adapters": []}
    stub = f"class MoleculeBatch:\n    def read_sdf(path: {annotation}) -> MoleculeBatch: ...\n"
    assert check_contract(stub, json.dumps(document)) == []


def test_sdf_path_gate_checks_actual_extraction_not_just_annotations():
    import os

    def valid(path):
        value = os.fspath(path)
        if not isinstance(value, str):
            raise TypeError("expected text path")

    assert _checker.check_path_input(valid, {}, "path") == []
    assert _checker.check_path_input(lambda path: str(path), {}, "path")
    assert _checker.check_path_input(lambda path: os.fspath(path), {}, "path")

    def rejects_valid_paths(path):
        _ = os.fspath(path)
        raise TypeError("broken conversion after fspath")

    assert _checker.check_path_input(rejects_valid_paths, {}, "path")


@pytest.mark.parametrize("annotation,valid", [("str | os.PathLike[str] | None", True), ("typing.Optional[typing.Union[str, os.PathLike[str]]]", True), ("str | os.PathLike[str]", False), ("str | os.PathLike[str] | bytes | None", False)])
def test_sdf_path_gate_preserves_optional_report_contract(annotation: str, valid: bool):
    row = dict(_entry("MoleculeBatch.write_sdf_with_params", "write_sdf_with_params", "type"), parameters=[
        {"name": "path", "type": "&str", "default": None},
        {"name": "report_path", "type": "Option<&str>", "default": None},
    ])
    stub = f"class MoleculeBatch:\n    def write_sdf_with_params(self, path: str | os.PathLike[str], report_path: {annotation}) -> None: ...\n"
    errors = check_contract(stub, json.dumps({"entries": [row], "python_adapters": []}))
    assert (not errors) == valid


def test_sdf_path_gate_follows_registered_path_to_python_directory_name():
    row = dict(_entry("MoleculeBatch.write_sdf_files", "write_sdf_files", "type"), parameters=[{"name": "path", "type": "&str", "default": None}])
    document = {"entries": [row], "python_adapters": []}
    stub = "class MoleculeBatch:\n    def write_sdf_files(self, out_dir: str | os.PathLike[str]) -> None: ...\n"
    assert check_contract(stub, json.dumps(document)) == []
    assert check_contract(stub.replace("out_dir", "unrelated"), json.dumps(document))


def test_native_scalar_projection_is_explicit_not_a_missing_class_exemption():
    import types
    row = _entry("types.BioAtomId", "BioAtomId", "type", "type")
    row["python_native"] = "builtins.int"
    document = {"entries": [row], "python_adapters": []}
    assert check_contract("", json.dumps(document)) == []
    assert check_runtime(types.SimpleNamespace(), document) == []
    del row["python_native"]
    assert len(check_contract("", json.dumps(document))) == 1
    assert len(check_runtime(types.SimpleNamespace(), document)) == 1
    row["python_native"] = "object"
    assert len(check_contract("", json.dumps(document))) == 1
    assert len(check_runtime(types.SimpleNamespace(), document)) == 1


def test_native_union_requires_declared_and_actual_component_classes():
    import types
    row = _entry("types.Record", "Record", "type", "type")
    row["python_native"] = "Molecule | BatchError"
    entries = [row] + [_entry("types." + name, name, "type", "type") for name in ("Molecule", "BatchError")]
    document = {"entries": entries, "python_adapters": []}
    stub = "class Molecule: ...\nclass BatchError: ...\n"
    assert check_contract(stub, json.dumps(document)) == []
    module = types.SimpleNamespace(Molecule=type("Molecule", (), {}), BatchError=type("BatchError", (), {}))
    assert check_runtime(module, document) == []
    assert check_contract("class Molecule: ...\n", json.dumps(document))
    del module.BatchError
    assert check_runtime(module, document)


def test_native_list_requires_its_declared_and_real_element_type():
    import types
    row = _entry("types.ProteinChainIter", "ProteinChainIter", "type", "type")
    row["python_native"] = "list[ProteinChainRef]"
    element = _entry("types.ProteinChainRef", "ProteinChainRef", "type", "type")
    document = {"entries": [row, element], "python_adapters": []}
    assert check_contract("class ProteinChainRef: ...\n", json.dumps(document)) == []
    module = types.SimpleNamespace(ProteinChainRef=type("ProteinChainRef", (), {}))
    assert check_runtime(module, document) == []
    assert check_contract("", json.dumps(document))
    assert check_runtime(types.SimpleNamespace(), document)
    row["python_native"] = "builtins.list"
    assert check_contract("class ProteinChainRef: ...\n", json.dumps(document))


def test_enum_name_is_checked_as_an_instance_descriptor_without_overriding_it():
    import enum
    import types
    class Kind(enum.IntEnum):
        AA = 1
    row = _entry("Kind.name", "name", "type")
    row["python_property"] = "getter"
    document = {"entries": [_entry("types.Kind", "Kind", "type", "type"), row], "python_adapters": []}
    assert check_runtime(types.SimpleNamespace(Kind=Kind), document) == []
    assert Kind.AA.name == "AA"


def test_optional_owned_getter_is_a_valid_optional_sequence_projection():
    node = lambda text: ast.parse(text, mode="eval").body
    assert _checker._read_type_matches(node("typing.Optional[list[int]]"), node("typing.Optional[typing.Sequence[int]]"))
    assert not _checker._read_type_matches(node("typing.Optional[list[int]]"), node("typing.Sequence[int]"))
    assert not _checker._read_type_matches(node("typing.Optional[list[str]]"), node("typing.Optional[typing.Sequence[int]]"))
    assert _checker._read_type_matches(node("dict[str, list[int]]"), node("typing.Optional[typing.Mapping[str, typing.Sequence[int]]]"))
    assert not _checker._read_type_matches(node("dict[int, list[int]]"), node("typing.Mapping[str, typing.Sequence[int]]"))
    assert not _checker._read_type_matches(node("dict[str, list[str]]"), node("typing.Mapping[str, typing.Sequence[int]]"))
    assert not _checker._read_type_matches(node("typing.Optional[dict[str, int]]"), node("typing.Mapping[str, int]"))


def test_dynamic_stubs_describe_only_real_exception_and_enum_exports():
    import enum
    import types
    class ParseError(ValueError):
        pass
    class Selection(str, enum.Enum):
        FIRST = "first"
    module = types.SimpleNamespace(ParseError=ParseError, Selection=Selection)
    entries = [_entry("types." + name, name, "type", "type") for name in ("ParseError", "Selection", "Missing")]
    text = _checker.dynamic_type_declarations(module, "", {"entries": entries})
    assert "class ParseError(builtins.ValueError)" in text
    assert "class Selection(builtins.str, enum.Enum)" in text
    assert "FIRST = 'first'" in text
    assert "Missing" not in text
    assert len(check_contract(text, json.dumps(entries))) == 1


def test_intenum_read_projection_requires_real_integer_inheritance():
    node = lambda text: ast.parse(text, mode="eval").body
    cls = ast.parse("class Format(enum.IntEnum): pass").body[0]
    assert _checker._read_type_matches(node("Format"), node("int"), {"Format": cls})
    assert not _checker._read_type_matches(node("Format"), node("str"), {"Format": cls})
    assert not _checker._read_type_matches(node("Format"), node("int"), {})
    cls = ast.parse("class Format(enum.Enum): pass").body[0]
    assert not _checker._read_type_matches(node("Format"), node("int"), {"Format": cls})


def test_enum_string_input_does_not_allow_string_output_or_arbitrary_unions():
    node = lambda text: ast.parse(text, mode="eval").body
    cls = ast.parse("class Mode:\n    First: typing.ClassVar[Mode]\n").body[0]
    classes = {"Mode": cls}
    assert _checker._read_type_matches(node("Mode"), node("Mode | builtins.str"), classes)
    assert not _checker._read_type_matches(node("str"), node("Mode | builtins.str"), classes)
    assert not _checker._read_type_matches(node("Mode"), node("Mode | int"), classes)
    assert not _checker._read_type_matches(node("Mode"), node("Mode | str"), {})


def test_configuration_snapshots_check_nested_values_and_float_bits_not_copy_identity():
    class Nested:
        def __init__(self, limit=100):
            self.limit = limit
    class Params:
        def __init__(self, limit=100):
            self.limit = limit
            self.zero = -0.0
        @property
        def nested(self):
            return Nested(self.limit)
    schemas = {
        Nested: [{"name": "limit"}],
        Params: [{"name": "nested"}, {"name": "zero"}],
    }
    value = Params()
    before = _checker._configuration_snapshot(value, schemas)
    assert value.nested != value.nested  # independent owned copies without __eq__
    assert _checker._configuration_snapshot(value, schemas) == before
    value.limit = 99
    assert _checker._configuration_snapshot(value, schemas) != before
    value.limit = 100
    value.zero = 0.0
    assert _checker._configuration_snapshot(value, schemas) != before


def _entry(semantic_id: str, name: str, owner: str = "molecule", item: str = "callable"):
    return {
        "semantic_id": semantic_id,
        "python_name": name,
        "owner": owner,
        "item": item,
        "feature": "cap-valence",
        "status": "experimental",
    }


def test_module_class_factory_and_inplace_declarations_are_required() -> None:
    entries = [
        _entry("search.parse_smarts", "parse_smarts", "module"),
        _entry("types.QueryGraph", "QueryGraph", "type", "type"),
        _entry("QueryGraph.from_smarts", "from_smarts", "type"),
        _entry("Molecule.assign_valence_", "assign_valence_"),
    ]
    stub = """
def parse_smarts(text): ...
class QueryGraph:
    @staticmethod
    def from_smarts(text): ...
class Molecule:
    def assign_valence_(self): ...
"""
    assert check_contract(stub, json.dumps(entries)) == []


def test_runtime_and_stub_both_omitting_inplace_methods_cannot_pass() -> None:
    entries = [
        _entry("Molecule.assign_valence_", "assign_valence_"),
        _entry("Molecule.assign_valence_with_params_", "assign_valence_with_params_"),
    ]
    assert check_contract(
        "class Molecule:\n    def with_assigned_valence(self): ...\n", json.dumps(entries)
    ) == [
        "Molecule.assign_valence_ -> cosmolkit.Molecule.assign_valence_ [feature=cap-valence]",
        "Molecule.assign_valence_with_params_ -> cosmolkit.Molecule.assign_valence_with_params_ [feature=cap-valence]",
    ]


def test_comments_strings_attributes_and_wrong_owners_do_not_count() -> None:
    entries = [_entry("Molecule.assign_valence_", "assign_valence_")]
    for stub in [
        "# def assign_valence_(self): ...\n",
        'text = "def assign_valence_(self): ..."\n',
        "class Molecule:\n    assign_valence_: object\n",
        "def assign_valence_(self): ...\n",
        "class Other:\n    def assign_valence_(self): ...\n",
        "class Molecule:\n    @property\n    def assign_valence_(self): ...\n",
        "class Molecule:\n    @property\n    def assign_valence_(self): ...\n    @assign_valence_.setter\n    def assign_valence_(self, value): ...\n",
    ]:
        assert len(check_contract(stub, json.dumps(entries))) == 1


def test_only_compiled_rows_are_required_without_status_exemptions() -> None:
    # The generator supplies only cfg-enabled rows. No second feature list here.
    entries = [_entry("Molecule.assign_valence_", "assign_valence_")]
    assert check_contract(
        "class Molecule:\n    def assign_valence_(self): ...\n", json.dumps(entries)
    ) == []
    assert len(check_contract("", json.dumps(entries))) == 1


def test_registered_type_projection_and_exact_constructor_name() -> None:
    entries = [
        _entry("types.RustParams", "PythonParams", "type", "type"),
        _entry("RustParams.new", "new", "type"),
    ]
    assert check_contract(
        "class PythonParams:\n    @staticmethod\n    def new(): ...\n", json.dumps(entries)
    ) == []
    # Do not invent an unregistered new -> __new__ compatibility mapping.
    assert len(check_contract(
        "class PythonParams:\n    def __new__(cls): ...\n", json.dumps(entries)
    )) == 1


def test_invalid_stub_and_registry_fail_explicitly() -> None:
    with pytest.raises(SyntaxError):
        _ = check_contract("class:", "[]")
    with pytest.raises(json.JSONDecodeError):
        _ = check_contract("", "not json")


def test_explicit_constructor_projection() -> None:
    entries = [_entry("Params.new", "__new__", "type")]
    assert check_contract("class Params:\n    def __new__(cls): ...\n", json.dumps(entries)) == []
    assert len(check_contract("class Params:\n    def new(cls): ...\n", json.dumps(entries))) == 1


def test_explicit_property_access_requires_descriptor_not_method() -> None:
    getter = {**_entry("Params.limit", "limit", "type"), "python_property": "getter"}
    setter = {**_entry("Params.set_limit", "limit", "type"), "python_property": "setter"}
    readonly = "class Params:\n    @property\n    def limit(self) -> int: ...\n"
    writable = readonly + "    @limit.setter\n    def limit(self, value: int) -> None: ...\n"
    assert check_contract(readonly, json.dumps([getter])) == []
    assert len(check_contract(readonly, json.dumps([setter]))) == 1
    assert check_contract(writable, json.dumps([getter, setter])) == []
    assert check_contract("class Params:\n    limit: int\n", json.dumps([getter, setter])) == []
    assert len(check_contract("class Params:\n    def limit(self): ...\n", json.dumps([getter, setter]))) == 2


def _configuration_document():
    parameter = {
        **_entry("types.SearchParams", "SearchParams", "type", "type"),
        "rust_path": "crate::SearchParams", "role": "parameter",
        "constructor": "SearchParams.new",
        "fields": [{"name": "limit", "type": "usize", "default": "100"}],
    }
    plain = {**_entry("Molecule.matches", "matches"), "parameters": []}
    configured = {
        **_entry("Molecule.matches_with_params", "matches_with_params"),
        "parameters": [{"name": "params", "type": "&crate::SearchParams", "default": None}],
    }
    return {"entries": [parameter, plain, configured], "python_adapters": []}


_CONFIGURATION_STUB = """
from typing import overload
class SearchParams:
    limit: int
    def __new__(cls, *, limit: int = 100) -> SearchParams: ...
class Molecule:
    @overload
    def matches(self, params: SearchParams, /) -> list[int]: ...
    @overload
    def matches(self, *, limit: int = 100) -> list[int]: ...
    def matches_with_params(self, params: SearchParams) -> list[int]: ...
"""


def test_parameter_contract_requires_writable_fields_not_a_params_suffix():
    document = _configuration_document()
    assert check_contract(_CONFIGURATION_STUB, json.dumps(document)) == []
    assert check_contract(_CONFIGURATION_STUB.replace("SearchParams", "Settings"), json.dumps(document).replace("SearchParams", "Settings")) == []
    for field in ("limit: Final[int]", "limit: ClassVar[int]", "@property\n    def limit(self) -> int: ..."):
        stub = _CONFIGURATION_STUB.replace("limit: int\n", field + "\n", 1)
        assert any("must be readable and writable" in error for error in check_contract(stub, json.dumps(document)))


def test_explicit_registry_configuration_schema_needs_no_rust_new_and_checks_types():
    document = _configuration_document()
    row = document["entries"][0]
    row["constructor"] = None
    row["fields"] = None
    row["python_fields"] = [{"name": "limit", "type": "builtins.int", "default": "100"}]
    assert check_contract(_CONFIGURATION_STUB, json.dumps(document)) == []
    row["python_fields"][0]["type"] = "builtins.str"
    assert any("constructor type differs from registry" in error for error in check_contract(_CONFIGURATION_STUB, json.dumps(document)))
    row["python_fields"][0] = {"name": "limit", "type": "builtins.int", "default": "DEFAULT_LIMIT"}
    assert any("default differs from registry" in error for error in check_contract(_CONFIGURATION_STUB, json.dumps(document)))
    row["python_fields"] = None
    assert any("registered configuration constructor schema" in error for error in check_contract(_CONFIGURATION_STUB, json.dumps(document)))


def test_configuration_signature_defaults_and_overloads_cannot_drift():
    document = json.dumps(_configuration_document())
    for old, new, expected in (
        ("limit: int = 100) -> SearchParams", "limit: int = 99) -> SearchParams", "default differs from registry"),
        ("limit: int = 100) -> list[int]", "limit: int = 99) -> list[int]", "type/default differs"),
        ("params: SearchParams, /", "params: object, /", "missing SearchParams instance"),
        ("params: SearchParams, /", "params: SearchParams, *, limit: int = 100", "mutually exclusive"),
        ("*, limit: int = 100) -> list[int]", "limit: int = 100) -> list[int]", "missing keyword-only"),
        ("@overload", "# not an overload", "explicit overload"),
    ):
        errors = check_contract(_CONFIGURATION_STUB.replace(old, new), document)
        assert any(expected in error for error in errors), errors


def test_parameter_schema_and_class_cannot_be_omitted_from_both_projections():
    document = _configuration_document()
    document["entries"] = document["entries"][:1]
    assert check_contract("", json.dumps(document))
    document["entries"][0]["fields"] = None
    assert any("registered configuration constructor schema" in error for error in check_contract(_CONFIGURATION_STUB, json.dumps(document)))


def test_actual_readonly_extension_style_descriptor_cannot_be_hidden_by_writable_stub():
    from types import SimpleNamespace

    class SearchParams:
        __slots__ = ("_limit",)

        def __init__(self, limit=100):
            self._limit = limit

        @property
        def limit(self):
            return self._limit

    document = _configuration_document()
    document["entries"] = document["entries"][:1]
    assert any("actual setter missing" in error for error in check_runtime(SimpleNamespace(SearchParams=SearchParams), document))
    def set_limit(self, value):
        if not isinstance(value, int):
            raise TypeError("limit must be an integer")
        self._limit = value
    SearchParams.limit = SearchParams.limit.setter(set_limit)
    assert check_runtime(SimpleNamespace(SearchParams=SearchParams), document) == []
    def broken_set_limit(self, value):
        self._limit = 0
        set_limit(self, value)
    SearchParams.limit = SearchParams.limit.setter(broken_set_limit)
    assert any("failed assignment changed" in error for error in check_runtime(SimpleNamespace(SearchParams=SearchParams), document))


def test_runtime_gate_rejects_noop_setter_even_when_old_value_and_invalid_input_checks_pass():
    from types import SimpleNamespace

    class SearchParams:
        __slots__ = ("_limit",)

        def __init__(self, limit=100):
            self._limit = limit

        @property
        def limit(self):
            return self._limit

        @limit.setter
        def limit(self, value):
            if not isinstance(value, int):
                raise TypeError("integer required")
            # Deliberately ignore every valid assignment.

    document = _configuration_document()
    document["entries"] = document["entries"][:1]
    errors = check_runtime(SimpleNamespace(SearchParams=SearchParams), document)
    assert any("effective assignment differs from constructor" in error for error in errors)


def test_runtime_gate_rejects_unregistered_dynamic_attributes():
    from types import SimpleNamespace

    class SearchParams:
        def __init__(self, limit=100):
            self.limit = limit

    document = _configuration_document()
    document["entries"] = document["entries"][:1]
    assert any("unknown configuration field was silently accepted" in error for error in check_runtime(SimpleNamespace(SearchParams=SearchParams), document))


def test_callable_configuration_requires_registered_behavior_not_just_a_typed_stub():
    document = _configuration_document()
    document["entries"] = document["entries"][:1]
    document["entries"][0]["python_fields"] = [{"name": "limit", "type": "typing.Optional[typing.Callable[..., bool]]", "default": "None"}]
    stub = "class SearchParams:\n    limit: typing.Optional[typing.Callable[..., bool]]\n    def __new__(cls, *, limit: typing.Optional[typing.Callable[..., bool]] = None) -> SearchParams: ...\n"
    assert any("requires an executable callback contract" in error for error in check_contract(stub, json.dumps(document)))


def test_callback_gate_detects_binding_that_accepts_and_discards_callbacks():
    from types import SimpleNamespace

    class SearchParams:
        def __init__(self, limit=None):
            self.limit = limit

    class Molecule:
        @classmethod
        def from_smiles(cls, text):
            return cls()

        def to_smiles(self):
            return "CC"

        def matches(self, query, **kwargs):
            return [SimpleNamespace(atom_mapping=lambda: [0, 1])]

        def matches_with_params(self, query, params):
            return [SimpleNamespace(atom_mapping=lambda: [0, 1])]  # Never invoke params.limit.

    document = _configuration_document()
    document["entries"][0]["fields"] = [{"name": "limit", "type": "typing.Optional[typing.Callable[..., bool]]", "default": "None", "callback": "final_match"}]
    module = SimpleNamespace(Molecule=Molecule, SearchParams=SearchParams, parse_smarts=lambda _: object())
    assert any("callback not invoked" in error for error in _checker.check_callback_runtime(module, document))
