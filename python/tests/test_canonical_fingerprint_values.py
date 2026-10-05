"""Fixed canonical value proposals; populated metadata awaits explicit calls."""

import ast
import inspect
from pathlib import Path
from typing import Callable, cast

import cosmolkit
import pytest

STUB = Path(__file__).resolve().parents[1] / "cosmolkit.pyi"


def test_fingerprint_additional_output_domain_explicit_public_name():
    output = cosmolkit.FingerprintAdditionalOutput()
    assert type(output).__name__ == "FingerprintAdditionalOutput"
    assert not hasattr(cosmolkit, "AdditionalOutput")
    assert repr(output).startswith("FingerprintAdditionalOutput(")
    assert [
        output.atom_counts(),
        output.atom_to_bits(),
        output.bit_info_map(),
        output.bit_paths(),
        output.atoms_per_bit(),
    ] == [None] * 5
    assert isinstance(cosmolkit.FingerprintAdditionalOutput.default(), type(output))
    names = {
        node.name
        for node in ast.parse(STUB.read_text()).body
        if isinstance(node, ast.ClassDef)
    }
    assert "FingerprintAdditionalOutput" in names
    assert "AdditionalOutput" not in names


DEFAULTS: dict[str, object] = {
    "radius": 3,
    "include_chirality": False,
    "use_bond_types": True,
    "include_ring_membership": True,
    "only_nonzero_invariants": False,
    "include_redundant_environments": False,
    "fp_size": 2048,
    "count_simulation": False,
    "count_bounds": [1, 2, 4, 8],
    "bits_per_feature": 1,
}
U32_FIELDS = ("radius", "fp_size", "bits_per_feature")
INVALID_U32 = (
    ("negative", -1, OverflowError),
    ("overflow", 4294967296, OverflowError),
    ("float", 1.5, TypeError),
    ("string", "3", TypeError),
)


def stub_class(name: str) -> ast.ClassDef:
    return next(node for node in ast.parse(STUB.read_text()).body
                if isinstance(node, ast.ClassDef) and node.name == name)


def stub_methods(name: str) -> dict[str, ast.FunctionDef]:
    return {node.name: node for node in stub_class(name).body if isinstance(node, ast.FunctionDef)}


def expression(node: ast.expr | None) -> str:
    assert node is not None
    return ast.unparse(node)


def construct_morgan(**fields: object) -> cosmolkit.MorganParams:
    constructor = cast(Callable[..., cosmolkit.MorganParams], cosmolkit.MorganParams)
    return constructor(**fields)


class TestMorganParams:
    def test_defaults(self):
        params = cosmolkit.MorganParams()
        assert cosmolkit._binding_profile == "canonical-bootstrap"
        assert {field: cast(object, getattr(params, field)) for field in DEFAULTS} == DEFAULTS
        for field, value in DEFAULTS.items():
            assert type(cast(object, getattr(params, field))) is type(value)
        assert cosmolkit.MorganParams(count_bounds=None).count_bounds == [1, 2, 4, 8]

    def test_all_fields(self):
        values = dict(radius=0, include_chirality=True, use_bond_types=False,
                      include_ring_membership=False, only_nonzero_invariants=True,
                      include_redundant_environments=True, fp_size=0,
                      count_simulation=True, count_bounds=[4294967295, 0, 0, 2],
                      bits_per_feature=0)
        params = construct_morgan(**values)
        assert {field: cast(object, getattr(params, field)) for field in values} == values

    @pytest.mark.parametrize("field,value", [(field, value) for field in U32_FIELDS for value in (0, 4294967295)],
                             ids=[f"{field}-{value}" for field in U32_FIELDS for value in ("zero", "max")])
    def test_u32_endpoints(self, field: str, value: int):
        params = construct_morgan(**{field: value})
        assert cast(object, getattr(params, field)) == value

    @pytest.mark.parametrize("field,value,error", [(field, value, error) for field in U32_FIELDS for _, value, error in INVALID_U32],
                             ids=[f"{field}-{name}" for field in U32_FIELDS for name, _, _ in INVALID_U32])
    def test_invalid_u32(self, field: str, value: object, error: type[Exception]):
        with pytest.raises(error):
            _ = construct_morgan(**{field: value})

    @pytest.mark.parametrize("value,error", [(value, error) for _, value, error in INVALID_U32],
                             ids=[name for name, _, _ in INVALID_U32])
    def test_invalid_bounds(self, value: object, error: type[Exception]):
        with pytest.raises(error):
            _ = construct_morgan(count_bounds=[value])

    def test_explicit_empty_bounds(self):
        assert cosmolkit.MorganParams(count_bounds=[]).count_bounds == []
        assert cosmolkit.MorganParams(count_bounds=[0, 4294967295]).count_bounds == [0, 4294967295]

    def test_owned_bounds(self):
        supplied = [8, 1, 1, 0]
        params = cosmolkit.MorganParams(count_bounds=supplied)
        supplied.append(99)
        first = params.count_bounds
        first[0] = 77
        first.append(123)
        assert params.count_bounds == [8, 1, 1, 0]
        assert first is not params.count_bounds
        default_result = cosmolkit.MorganParams().count_bounds
        default_result.clear()
        assert cosmolkit.MorganParams().count_bounds == [1, 2, 4, 8]

    @pytest.mark.parametrize("field", list(DEFAULTS), ids=list(DEFAULTS))
    def test_frozen(self, field: str):
        params = cosmolkit.MorganParams()
        with pytest.raises(AttributeError):
            setattr(params, field, DEFAULTS[field])
        with pytest.raises(AttributeError):
            delattr(params, field)
        assert cast(object, getattr(params, field)) == DEFAULTS[field]

    def test_signature(self):
        signature = inspect.signature(cosmolkit.MorganParams)
        assert list(signature.parameters) == list(DEFAULTS)
        expected = dict(DEFAULTS, count_bounds=None)
        assert {name: cast(object, param.default) for name, param in signature.parameters.items()} == expected
        assert all(param.kind is inspect.Parameter.KEYWORD_ONLY for param in signature.parameters.values())

    @pytest.mark.parametrize("case", ["positional", "unknown", "bounds_scalar"], ids=["positional", "unknown", "bounds_scalar"])
    def test_constructor_errors(self, case: str):
        constructor = cast(Callable[..., cosmolkit.MorganParams], cosmolkit.MorganParams)
        with pytest.raises(TypeError):
            if case == "positional":
                _ = constructor(3)
            elif case == "unknown":
                _ = constructor(unknown=True)
            else:
                _ = constructor(count_bounds=8)

    def test_stub(self):
        methods = stub_methods("MorganParams")
        assert set(methods) == {"__new__", *DEFAULTS}
        constructor = methods["__new__"]
        assert [arg.arg for arg in constructor.args.args] == ["cls"]
        assert [arg.arg for arg in constructor.args.kwonlyargs] == list(DEFAULTS)
        expected = dict(DEFAULTS, count_bounds=None)
        assert [ast.literal_eval(value) for value in constructor.args.kw_defaults if value is not None] == list(expected.values())
        assert expression(constructor.returns) == "MorganParams"
        for field, value in DEFAULTS.items():
            assert [ast.unparse(node) for node in methods[field].decorator_list] == ["property"]
            expected_type = "builtins.list[builtins.int]" if field == "count_bounds" else ("builtins.bool" if type(value) is bool else "builtins.int")
            assert expression(methods[field].returns) == expected_type
        assert [expression(arg.annotation) for arg in constructor.args.kwonlyargs] == [
            "typing.Optional[typing.Sequence[builtins.int]]" if field == "count_bounds" else ("builtins.bool" if type(value) is bool else "builtins.int")
            for field, value in DEFAULTS.items()
        ]


OUTPUT_FIELDS = ("atom_counts", "atom_to_bits", "bit_info_map", "bit_paths", "atoms_per_bit")
OUTPUT_METHODS = ("default", *("allocate_" + field for field in OUTPUT_FIELDS), *OUTPUT_FIELDS)
MASK_IDS = [f"mask{mask:02d}" for mask in range(32)]


def output_method(output: cosmolkit.FingerprintAdditionalOutput, name: str) -> Callable[..., object]:
    method = cast(object, getattr(output, name))
    assert callable(method)
    return method


def output_state(output: cosmolkit.FingerprintAdditionalOutput) -> list[object]:
    return [output_method(output, field)() for field in OUTPUT_FIELDS]


def allocated_output(mask: int) -> cosmolkit.FingerprintAdditionalOutput:
    output = cosmolkit.FingerprintAdditionalOutput.default()
    for index, field in enumerate(OUTPUT_FIELDS):
        if mask & (1 << index):
            assert output_method(output, "allocate_" + field)() is None
    return output


def expected_state(mask: int) -> list[object]:
    return [([] if index < 2 else {}) if mask & (1 << index) else None for index in range(5)]


class TestAdditionalOutput:
    def test_default(self):
        output = cosmolkit.FingerprintAdditionalOutput.default()
        assert isinstance(output, cosmolkit.FingerprintAdditionalOutput)
        assert output_state(output) == [None, None, None, None, None]

    @pytest.mark.parametrize("mask", range(32), ids=MASK_IDS)
    def test_allocation_masks(self, mask: int):
        output = allocated_output(mask)
        assert output_state(output) == expected_state(mask)
        for index, result in enumerate(output_state(output)):
            assert type(result) is (type(None) if not mask & (1 << index) else (list if index < 2 else dict))

    @pytest.mark.parametrize("mask,field", [(mask, field) for mask in range(32) for field in OUTPUT_FIELDS],
                             ids=[f"mask{mask:02d}-{field}" for mask in range(32) for field in OUTPUT_FIELDS])
    def test_reallocation(self, mask: int, field: str):
        output = allocated_output(mask)
        selected = OUTPUT_FIELDS.index(field)
        for _ in range(3):
            previous_results = output_state(output)
            assert output_method(output, "allocate_" + field)() is None
            assert output_state(output) == expected_state(mask | (1 << selected))
            # Allocations must not mutate previously detached read results.
            assert previous_results == expected_state(mask)
            mask |= 1 << selected

    @pytest.mark.parametrize("field", OUTPUT_FIELDS, ids=OUTPUT_FIELDS)
    def test_owned_results(self, field: str):
        output = allocated_output(31)
        first = output_method(output, field)()
        second = output_method(output, field)()
        assert first is not second
        if isinstance(first, list):
            owned_list = cast(list[object], first)
            owned_list.append(123 if field == "atom_counts" else [4294967295, 4294967296])
        else:
            assert isinstance(first, dict)
            owned_dict = cast(dict[int, object], first)
            owned_dict[18446744073709551615] = [(4294967295, 0)] if field == "bit_info_map" else [[-2147483648, 2147483647]]
        assert output_state(output) == [[], [], {}, {}, {}]
        assert second == ([] if field in OUTPUT_FIELDS[:2] else {})
        assert output_method(output, "allocate_" + field)() is None
        assert first != output_method(output, field)()
        assert output_state(output) == [[], [], {}, {}, {}]
        # Values inserted above are only Python copies; canonical populated
        # numeric/nested transport is explicitly not tested in this packet.

    def test_independent_defaults(self):
        first = cosmolkit.FingerprintAdditionalOutput.default()
        second = cosmolkit.FingerprintAdditionalOutput.default()
        assert first is not second
        first.allocate_atom_counts()
        first.allocate_atoms_per_bit()
        assert output_state(first) == [[], None, None, None, {}]
        assert output_state(second) == [None, None, None, None, None]

    def test_surface_signatures(self):
        surface = {name for name in vars(cosmolkit.FingerprintAdditionalOutput) if not name.startswith("_")}
        assert surface == set(OUTPUT_METHODS)
        for name in OUTPUT_METHODS:
            descriptor = cast(object, getattr(cosmolkit.FingerprintAdditionalOutput, name))
            assert callable(descriptor)
            signature = inspect.signature(descriptor)
            assert list(signature.parameters) == ([] if name == "default" else ["self"])
            assert all(cast(object, param.default) is inspect.Parameter.empty for param in signature.parameters.values())
            assert all(param.kind is inspect.Parameter.POSITIONAL_ONLY for param in signature.parameters.values())
        assert isinstance(vars(cosmolkit.FingerprintAdditionalOutput)["default"], staticmethod)

    @pytest.mark.parametrize("name,kind", [(name, kind) for name in OUTPUT_METHODS for kind in ("positional", "keyword")],
                             ids=[f"{name}-{kind}" for name in OUTPUT_METHODS for kind in ("positional", "keyword")])
    def test_argument_errors(self, name: str, kind: str):
        output = cosmolkit.FingerprintAdditionalOutput.default()
        method = output_method(output, name)
        with pytest.raises(TypeError):
            if kind == "positional":
                _ = method(0)
            else:
                _ = method(value=0)
        assert output_state(output) == [None, None, None, None, None]

    def test_constructor_and_no_aliases(self):
        constructor = cast(Callable[..., object], cosmolkit.FingerprintAdditionalOutput)
        # Source default field initializers and canonical FingerprintAdditionalOutput.new
        # project to a Python constructor; retain negative argument coverage.
        output = cast(cosmolkit.FingerprintAdditionalOutput, constructor())
        assert output_state(output) == [None, None, None, None, None]
        other = cast(cosmolkit.FingerprintAdditionalOutput, constructor())
        output.allocate_atom_counts()
        assert output_state(other) == [None, None, None, None, None]
        for args, kwargs in [((0,), {}), ((), {"value": 0})]:
            with pytest.raises(TypeError):
                _ = constructor(*args, **kwargs)
        assert not hasattr(cosmolkit.FingerprintAdditionalOutput, "new")
        for field in OUTPUT_FIELDS:
            assert not hasattr(cosmolkit.FingerprintAdditionalOutput, "get_" + field)
        output = cosmolkit.FingerprintAdditionalOutput.default()
        with pytest.raises(AttributeError):
            setattr(output, "extra", [])
        with pytest.raises(AttributeError):
            setattr(output, "atom_counts", [])
        assert output_state(output) == [None, None, None, None, None]
        assert repr(output) == "FingerprintAdditionalOutput(atom_to_bits=false, bit_info_map=false, bit_paths=false, atom_counts=false, atoms_per_bit=false)"
        for method in ("__repr__",):
            with pytest.raises(TypeError):
                _ = getattr(output, method)(0)

    def test_stub(self):
        methods = stub_methods("FingerprintAdditionalOutput")
        assert set(methods) == {"__new__", "__repr__", *OUTPUT_METHODS}
        expected_returns = {
            "atom_counts": "typing.Optional[builtins.list[builtins.int]]",
            "atom_to_bits": "typing.Optional[builtins.list[builtins.list[builtins.int]]]",
            "bit_info_map": "typing.Optional[builtins.dict[builtins.int, builtins.list[tuple[builtins.int, builtins.int]]]]",
            "bit_paths": "typing.Optional[builtins.dict[builtins.int, builtins.list[builtins.list[builtins.int]]]]",
            "atoms_per_bit": "typing.Optional[builtins.dict[builtins.int, builtins.list[builtins.list[builtins.int]]]]",
        }
        for name, method in methods.items():
            assert [arg.arg for arg in method.args.args] == (["cls"] if name == "__new__" else ([] if name == "default" else ["self"]))
            assert not method.args.defaults and not method.args.kwonlyargs
            assert not method.args.kw_defaults and method.args.vararg is None and method.args.kwarg is None
            expected = "FingerprintAdditionalOutput" if name in ("default", "__new__") else ("builtins.str" if name == "__repr__" else ("None" if name.startswith("allocate_") else expected_returns[name]))
            assert expression(method.returns) == expected
            assert [ast.unparse(node) for node in method.decorator_list] == (["staticmethod"] if name == "default" else [])
