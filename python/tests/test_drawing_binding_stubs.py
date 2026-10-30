"""Inspect actual generated selected declarations and native signatures."""

import ast
import inspect
from pathlib import Path
from typing import cast

import cosmolkit


STUB = Path(__file__).resolve().parents[1] / "artifacts/smoke/package/cosmolkit/__init__.pyi"


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
        elif isinstance(node, ast.FunctionDef) and any(ast.unparse(d) == "builtins.property" for d in node.decorator_list):
            fields[node.name] = ast.unparse(required_expression(node.returns))
    return fields


def test_selected_stub_classes_and_methods():
    classes = declarations()
    expected_classes = {"Molecule", "Coordinate2DParams", "DrawingError", "DrawingWriteError", "OperationError", "SmilesParseParams", "SmilesWriteParams", "SmilesError", "SmilesWriteError", "MorganReadError", "FingerprintError", "Fingerprint", "SparseBitFingerprint", "SparseCountFingerprint", "SparseCountFingerprint32", "MorganGeneratorParams", "FingerprintAdditionalOutput", "Element", "ElementInfo", "DescriptorReadError", "DescriptorError"}
    # This focused regression does not forbid other registered domains.
    # The generator checks the complete enabled registry against the stub.
    assert expected_classes <= classes.keys()
    methods = {}
    for node in classes["Molecule"].body:
        if isinstance(node, ast.FunctionDef):
            methods.setdefault(node.name, node)
    expected = {"from_smiles": "Molecule", "num_atoms": "int", "num_bonds": "int",
                "to_smiles": "str", "coordinates_2d": "numpy.ndarray[tuple[int, int], numpy.dtype[numpy.float64]] | None",
                "has_2d_coordinates": "bool", "add_2d_coordinates_": "None",
                "add_2d_coordinates_with_params_": "None", "with_2d_coordinates": "Molecule", "with_2d_coordinates_with_params": "Molecule",
                "to_svg": "str", "to_png": "bytes", "write_svg": "None", "write_png": "None"}
    expected.update({"new": "Molecule", "from_smiles_with_params": "Molecule", "to_smiles_with_params": "str",
        "fingerprint_morgan": "Fingerprint", "fingerprint_morgan_sparse": "SparseBitFingerprint",
        "fingerprint_morgan_count": "SparseCountFingerprint32", "fingerprint_morgan_sparse_count": "SparseCountFingerprint"})
    expected.update({'hall_kier_alpha': 'float', 'hall_kier_alpha_with_contributions': 'tuple[float, list[float]]', 'kappa_1': 'float', 'kappa_2': 'float', 'kappa_3': 'float', 'phi': 'float', 'mqns': 'list[int]', 'chi_0_v': 'float', 'chi_1_v': 'float', 'chi_2_v': 'float', 'chi_3_v': 'float', 'chi_4_v': 'float', 'chi_n_v': 'float', 'chi_0_n': 'float', 'chi_1_n': 'float', 'chi_2_n': 'float', 'chi_3_n': 'float', 'chi_4_n': 'float', 'chi_n_n': 'float'})
    assert expected.keys() <= methods.keys()
    for name, result in expected.items():
        assert ast.unparse(required_expression(methods[name].returns)) == result
    for name in ("to_svg", "to_png", "write_svg", "write_png"):
        args = methods[name].args
        assert not args.defaults and not args.kw_defaults
        names = [a.arg for a in args.posonlyargs + args.args]
        assert names == (["self", "path", "width", "height"] if name.startswith("write_") else ["self", "width", "height"])
        assert [ast.unparse(required_expression(a.annotation)) for a in args.args[-2:]] == ["SupportsIndex", "SupportsIndex"]
        # Inspect the actual native descriptor, erasing only reflection's Any
        # result to object; assertions still compare the observed signature.
        descriptor = cast(object, getattr(cosmolkit.Molecule, name))
        assert callable(descriptor)
        runtime = inspect.signature(descriptor)
        assert list(runtime.parameters) == names
        assert all(cast(object, p.default) is inspect.Parameter.empty for p in runtime.parameters.values())
    assert "builtins.property" in [ast.unparse(d) for d in methods["has_2d_coordinates"].decorator_list]
    inspect.signature(cosmolkit.Molecule.add_2d_coordinates_).bind(object())
    inplace = methods["add_2d_coordinates_with_params_"]
    assert [a.arg for a in inplace.args.posonlyargs + inplace.args.args] == ["self", "params"]
    assert ast.unparse(required_expression(inplace.args.args[0].annotation)) == "Coordinate2DParams"
    configured = methods["with_2d_coordinates_with_params"]
    assert [a.arg for a in configured.args.posonlyargs + configured.args.args] == ["self", "params"]
    assert ast.unparse(required_expression(configured.args.args[0].annotation)) == "Coordinate2DParams"


def test_parameter_properties_types_and_defaults():
    cls = declarations()["Coordinate2DParams"]
    methods = {n.name: n for n in cls.body if isinstance(n, ast.FunctionDef)
               and not any(isinstance(d, ast.Attribute) and d.attr == "setter" for d in n.decorator_list)}
    defaults = {"coordinate_map": {}, "canonical_orientation": True, "clear_existing_2d": True,
                "flips_per_sample": 0, "samples": 0, "sample_seed": 0,
                "permute_degree_four": False, "force_rdkit": False, "use_ring_templates": False,
                "use_path_angle_expansion": True}
    assert set(methods) == {"__new__", *defaults}
    runtime = inspect.signature(cosmolkit.Coordinate2DParams)
    assert {name: cast(object, parameter.default) for name, parameter in runtime.parameters.items()} == dict.fromkeys(defaults, Ellipsis)
    params = cosmolkit.Coordinate2DParams()
    assert {name: getattr(params, name) for name in defaults} == defaults
    constructor = methods["__new__"].args
    assert [a.arg for a in constructor.posonlyargs] == ["cls"]
    assert [a.arg for a in constructor.kwonlyargs] == list(defaults)
    assert all(ast.literal_eval(required_expression(v)) is Ellipsis for v in constructor.kw_defaults)
    for name in defaults:
        setters = [n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == name
                   and [ast.unparse(d) for d in n.decorator_list] == [name + ".setter"]]
        assert len(setters) == 1
        field = methods[name]
        assert [ast.unparse(d) for d in field.decorator_list] == ["builtins.property"]
        expected_type = "Mapping[int, tuple[float, ...]]" if name == "coordinate_map" else ("bool" if type(defaults[name]) is bool else "int")
        assert ast.unparse(required_expression(field.returns)) == expected_type


def test_generated_exception_and_profile_declarations():
    classes = declarations()
    for name in ("DrawingError", "OperationError"):
        assert [ast.unparse(base) for base in classes[name].bases] == ["ValueError"]
        fields = declared_fields(classes[name])
        assert fields["domain"] == "str"
        assert fields["variant"] == name + "Kind"
    assert [ast.unparse(base) for base in classes["DrawingWriteError"].bases] == ["ValueError"]
    fields = set(declared_fields(classes["DrawingError"]))
    assert {"domain", "variant", "width", "height", "field", "actual", "expected", "row", "reason"} <= fields
    assert issubclass(cosmolkit.DrawingError, ValueError)
    assert issubclass(cosmolkit.OperationError, ValueError)


def test_generated_drawing_write_error_projection():
    cls = declarations()["DrawingWriteError"]
    assert [ast.unparse(base) for base in cls.bases] == ["ValueError"]
    fields = declared_fields(cls)
    assert fields["domain"] == "str"
    assert fields["variant"] == "DrawingWriteErrorKind"
    assert fields["path"] == "Path | None"
    assert fields["source"] == "DrawingError | IoError | None"
    assert issubclass(cosmolkit.DrawingWriteError, ValueError)
