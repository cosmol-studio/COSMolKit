"""Inspect actual generated selected declarations and native signatures."""

import ast
import inspect
from pathlib import Path
from typing import cast

import cosmolkit


STUB = Path(__file__).resolve().parents[1] / "cosmolkit.pyi"


def declarations() -> dict[str, ast.ClassDef]:
    return {node.name: node for node in ast.parse(STUB.read_text()).body if isinstance(node, ast.ClassDef)}


def required_expression(node: ast.expr | None) -> ast.expr:
    assert node is not None
    return node


def declared_fields(cls: ast.ClassDef) -> dict[str, str]:
    fields: dict[str, str] = {}
    for node in cls.body:
        if isinstance(node, ast.AnnAssign):
            assert isinstance(node.target, ast.Name)
            fields[node.target.id] = ast.unparse(node.annotation)
    return fields


def test_selected_stub_classes_and_methods():
    classes = declarations()
    assert set(classes) == {"Molecule", "Coordinate2DParams", "DrawingError", "DrawingWriteError", "OperationError", "SmilesParseParams", "SmilesWriteParams", "SmilesError", "SmilesWriteError", "MorganReadError", "FingerprintError", "Fingerprint", "SparseBitFingerprint", "SparseCountFingerprint", "SparseCountFingerprint32", "MorganParams", "FingerprintAdditionalOutput", "Element", "ElementInfo", "DescriptorReadError", "DescriptorError"}
    methods = {n.name: n for n in classes["Molecule"].body if isinstance(n, ast.FunctionDef)}
    expected = {"from_smiles": "Molecule", "num_atoms": "builtins.int", "num_bonds": "builtins.int",
                "to_smiles": "builtins.str", "coordinates_2d": "typing.Optional[builtins.list[builtins.list[builtins.float]]]",
                "has_2d_coordinates": "builtins.bool", "compute_2d_coordinates_": "None",
                "compute_2d_coordinates_with_params_": "None", "with_2d_coordinates": "Molecule", "with_2d_coordinates_with_params": "Molecule",
                "to_svg": "builtins.str", "to_png": "builtins.bytes", "write_svg": "None", "write_png": "None"}
    expected.update({"new": "Molecule", "from_smiles_with_params": "Molecule", "to_smiles_with_params": "builtins.str",
        "morgan_fingerprint": "Fingerprint", "morgan_sparse_fingerprint": "SparseBitFingerprint",
        "morgan_count_fingerprint": "SparseCountFingerprint32", "morgan_sparse_count_fingerprint": "SparseCountFingerprint"})
    expected.update({'hall_kier_alpha': 'builtins.float', 'hall_kier_alpha_with_contributions': 'tuple[builtins.float, builtins.list[builtins.float]]', 'kappa_1': 'builtins.float', 'kappa_2': 'builtins.float', 'kappa_3': 'builtins.float', 'phi': 'builtins.float', 'mqns': 'builtins.list[builtins.int]', 'chi_0_v': 'builtins.float', 'chi_1_v': 'builtins.float', 'chi_2_v': 'builtins.float', 'chi_3_v': 'builtins.float', 'chi_4_v': 'builtins.float', 'chi_n_v': 'builtins.float', 'chi_0_n': 'builtins.float', 'chi_1_n': 'builtins.float', 'chi_2_n': 'builtins.float', 'chi_3_n': 'builtins.float', 'chi_4_n': 'builtins.float', 'chi_n_n': 'builtins.float'})
    assert set(methods) == set(expected)
    for name, result in expected.items():
        assert ast.unparse(required_expression(methods[name].returns)) == result
    for name in ("to_svg", "to_png", "write_svg", "write_png"):
        args = methods[name].args
        assert not args.defaults and not args.kw_defaults
        names = [a.arg for a in args.args]
        assert names == (["self", "path", "width", "height"] if name.startswith("write_") else ["self", "width", "height"])
        assert [ast.unparse(required_expression(a.annotation)) for a in args.args[-2:]] == ["builtins.int", "builtins.int"]
        # Inspect the actual native descriptor, erasing only reflection's Any
        # result to object; assertions still compare the observed signature.
        descriptor = cast(object, getattr(cosmolkit.Molecule, name))
        assert callable(descriptor)
        runtime = inspect.signature(descriptor)
        assert list(runtime.parameters) == names
        assert all(cast(object, p.default) is inspect.Parameter.empty for p in runtime.parameters.values())
    for name in ("has_2d_coordinates", "compute_2d_coordinates_"):
        assert [a.arg for a in methods[name].args.args] == ["self"]
        assert list(inspect.signature(getattr(cosmolkit.Molecule, name)).parameters) == ["self"]
    inplace = methods["compute_2d_coordinates_with_params_"]
    assert [a.arg for a in inplace.args.args] == ["self", "params"]
    assert ast.unparse(required_expression(inplace.args.args[1].annotation)) == "Coordinate2DParams"
    configured = methods["with_2d_coordinates_with_params"]
    assert [a.arg for a in configured.args.args] == ["self", "params"]
    assert ast.unparse(required_expression(configured.args.args[1].annotation)) == "Coordinate2DParams"


def test_nine_parameter_properties_types_and_defaults():
    cls = declarations()["Coordinate2DParams"]
    methods = {n.name: n for n in cls.body if isinstance(n, ast.FunctionDef)}
    defaults = {"coordinate_map": None, "canonical_orientation": False, "clear_existing_2d": True,
                "flips_per_sample": 0, "samples": 0, "sample_seed": 0,
                "permute_degree_four": False, "force_rdkit": False, "use_ring_templates": False}
    assert set(methods) == {"__new__", *defaults}
    runtime = inspect.signature(cosmolkit.Coordinate2DParams)
    assert {name: cast(object, parameter.default) for name, parameter in runtime.parameters.items()} == defaults
    constructor = methods["__new__"].args
    assert [a.arg for a in constructor.args] == ["cls", "coordinate_map"]
    assert len(constructor.defaults) == 1 and ast.literal_eval(constructor.defaults[0]) is None
    assert [a.arg for a in constructor.kwonlyargs] == list(defaults)[1:]
    assert [ast.literal_eval(required_expression(v)) for v in constructor.kw_defaults] == list(defaults.values())[1:]
    for name in defaults:
        field = methods[name]
        assert [ast.unparse(d) for d in field.decorator_list] == ["property"]
        expected_type = "builtins.dict[builtins.int, builtins.list[builtins.float]]" if name == "coordinate_map" else ("builtins.bool" if type(defaults[name]) is bool else "builtins.int")
        assert ast.unparse(required_expression(field.returns)) == expected_type


def test_generated_exception_and_profile_declarations():
    classes = declarations()
    for name in ("DrawingError", "DrawingWriteError", "OperationError"):
        assert [ast.unparse(base) for base in classes[name].bases] == ["builtins.ValueError"]
        fields = declared_fields(classes[name])
        assert fields["domain"] == fields["kind"] == "builtins.str"
    fields = set(declared_fields(classes["DrawingError"]))
    assert fields == {"domain", "kind", "width", "height", "field", "actual", "expected", "row", "reason"}
    assert cosmolkit._binding_profile in {"drawing-bindings", "canonical-bootstrap"}
    assert issubclass(cosmolkit.DrawingError, ValueError)
    assert issubclass(cosmolkit.OperationError, ValueError)


def test_generated_drawing_write_error_projection():
    cls = declarations()["DrawingWriteError"]
    assert [ast.unparse(base) for base in cls.bases] == ["builtins.OSError"]
    assert declared_fields(cls) == {"domain": "builtins.str", "kind": "builtins.str"}
    assert issubclass(cosmolkit.DrawingWriteError, OSError)
