"""Fixed canonical Python boundary cases; no oracle/reference preparation.

These are proposal tests. Native execution needs explicit ROOT authorization.
"""
from __future__ import annotations

import ast
import inspect
import json
import math
from pathlib import Path
from collections.abc import Callable
from typing import TypedDict, cast

import pytest
import cosmolkit as ck


class ScalarCase(TypedDict):
    canonical: str
    smiles: str
    expected: float
    tolerance_literal: str
    order: int | None


class ContributionCase(TypedDict):
    smiles: str
    alpha: float
    rows: list[float]


class MqnCase(TypedDict):
    smiles: str
    values: list[int]


class FixedCases(TypedDict):
    scalar_projection_cases: list[ScalarCase]
    hall_contributions: list[ContributionCase]
    mqn_cases: list[MqnCase]


DATA_PATH = Path(__file__).resolve().parents[2] / "testdata/descriptors/d02_canonical_fixed.json"
DATA = cast(FixedCases, json.loads(DATA_PATH.read_text()))
PARAMETERS: dict[str, str | None] = {'hall_kier_alpha': None, 'hall_kier_alpha_with_contributions': None, 'kappa_1': None, 'kappa_2': None, 'kappa_3': None, 'phi': None, 'mqns': 'force', 'chi_0_v': None, 'chi_1_v': None, 'chi_2_v': None, 'chi_3_v': None, 'chi_4_v': None, 'chi_n_v': 'order', 'chi_0_n': None, 'chi_1_n': None, 'chi_2_n': None, 'chi_3_n': None, 'chi_4_n': None, 'chi_n_n': 'order'}
PREPARED: list[str] = ['mqns', 'chi_0_v', 'chi_1_v', 'chi_2_v', 'chi_3_v', 'chi_4_v', 'chi_n_v', 'chi_0_n', 'chi_1_n', 'chi_2_n', 'chi_3_n', 'chi_4_n', 'chi_n_n']


def query(molecule: ck.Molecule, name: str) -> Callable[..., object]:
    value = cast(object, getattr(molecule, name))
    assert callable(value)
    return value


def close(actual: object, expected: float, tolerance: float) -> None:
    assert isinstance(actual, float)
    assert math.isfinite(actual) and abs(actual - expected) < tolerance


@pytest.mark.parametrize("case", DATA["scalar_projection_cases"], ids=[case["canonical"] for case in DATA["scalar_projection_cases"]])
def test_d02_scalar_queries_project_source_literals(case: ScalarCase):
    molecule = ck.Molecule.from_smiles(case["smiles"])
    before = (molecule.num_atoms(), molecule.num_bonds(), molecule.to_smiles(), molecule.coordinates_2d())
    method = query(molecule, case["canonical"])
    result = method() if case["order"] is None else method(case["order"])
    close(result, case["expected"], float(case["tolerance_literal"]))
    assert (molecule.num_atoms(), molecule.num_bonds(), molecule.to_smiles(), molecule.coordinates_2d()) == before


@pytest.mark.parametrize("case", DATA["hall_contributions"], ids=[case["smiles"] for case in DATA["hall_contributions"]])
def test_d02_contributions_tuple_atom_order_and_owned_rows(case: ContributionCase):
    molecule = ck.Molecule.from_smiles(case["smiles"])
    alpha, rows = molecule.hall_kier_alpha_with_contributions()
    close(alpha, case["alpha"], 1e-12)
    assert isinstance(rows, list) and len(rows) == molecule.num_atoms()
    for actual, expected in zip(rows, case["rows"], strict=True):
        close(actual, expected, 1e-12)
    rows[0] = 123.0
    _, fresh = molecule.hall_kier_alpha_with_contributions()
    assert fresh is not rows
    for actual, expected in zip(fresh, case["rows"], strict=True):
        close(actual, expected, 1e-12)


@pytest.mark.parametrize("case", DATA["mqn_cases"], ids=[case["smiles"] for case in DATA["mqn_cases"]])
def test_d02_mqns_42_order_force_default_and_owned_list(case: MqnCase):
    molecule = ck.Molecule.from_smiles(case["smiles"])
    expected = case["values"]
    assert len(expected) == 42
    before = molecule.to_smiles()
    result = molecule.mqns()
    assert isinstance(result, list) and all(type(x) is int for x in result)
    assert result == expected
    assert molecule.mqns(force=False) == expected
    assert molecule.mqns(force=True) == expected
    result[0] = 2**32 - 1
    assert molecule.mqns() == expected
    assert molecule.to_smiles() == before


@pytest.mark.parametrize("family", ["v", "n"])
def test_d02_generic_zero_differs_from_fixed_zero_and_wraps(family: str):
    molecule = ck.Molecule.from_smiles("C")
    assert query(molecule, f"chi_0_{family}")() == 0.0
    assert query(molecule, f"chi_n_{family}")(0) == 1.0
    assert query(molecule, f"chi_n_{family}")(2**32 - 1) == 0.0


@pytest.mark.parametrize("name", PREPARED)
def test_d02_raw_molecule_prepared_boundary_is_typed_and_preserving(name: str):
    molecule = ck.Molecule.from_smiles_with_params("CCC", ck.SmilesParseParams(sanitize=False, remove_hydrogens=False))
    before = (molecule.num_atoms(), molecule.num_bonds(), molecule.coordinates_2d())
    with pytest.raises(ck.DescriptorReadError) as caught:
        if PARAMETERS[name] == "force":
            _ = query(molecule, name)(False)
        elif PARAMETERS[name] == "order":
            _ = query(molecule, name)(2)
        else:
            _ = query(molecule, name)()
    error = caught.value
    assert isinstance(error, ValueError)
    assert error.domain == "descriptors" and error.kind == "MissingPreparedValence"
    assert error.__cause__ is None
    assert (molecule.num_atoms(), molecule.num_bonds(), molecule.coordinates_2d()) == before


@pytest.mark.parametrize("name,expected", [
    ("hall_kier_alpha", 0.0), ("hall_kier_alpha_with_contributions", (0.0, [0.0, 0.0, 0.0])),
    ("kappa_1", 3.0), ("kappa_2", 2.0), ("kappa_3", 0.0), ("phi", 2.0),
])
def test_d02_raw_topology_queries_do_not_prepare(name: str, expected: object):
    molecule = ck.Molecule.from_smiles_with_params("CCC", ck.SmilesParseParams(sanitize=False, remove_hydrogens=False))
    assert query(molecule, name)() == expected
    with pytest.raises(ck.DescriptorReadError) as caught:
        _ = molecule.chi_0_v()
    assert caught.value.kind == "MissingPreparedValence"


@pytest.mark.parametrize("name", list(PARAMETERS))
def test_d02_native_signatures_are_canonical_receiver_methods(name: str):
    descriptor = cast(object, getattr(ck.Molecule, name))
    assert callable(descriptor)
    signature = inspect.signature(descriptor)
    parameter = PARAMETERS[name]
    assert list(signature.parameters) == ["self"] + ([] if parameter is None else [parameter])
    for p in signature.parameters.values():
        observed_default = cast(object, p.default)
        if p.name == "force":
            assert observed_default is False
        else:
            assert observed_default is inspect.Parameter.empty
    assert not hasattr(ck, name)


@pytest.mark.parametrize("name", [n for n in PREPARED if n.startswith("chi_")])
def test_d02_chi_unmodeled_force_is_not_silently_accepted(name: str):
    molecule = ck.Molecule.from_smiles("C")
    with pytest.raises(TypeError):
        _ = query(molecule, name)(force=True)


@pytest.mark.parametrize("name", ["chi_n_v", "chi_n_n"])
@pytest.mark.parametrize("order", [-1, 2**32])
def test_d02_generic_order_rejects_values_outside_source_u32(name: str, order: int):
    molecule = ck.Molecule.from_smiles("C")
    with pytest.raises(OverflowError):
        _ = query(molecule, name)(order)


def test_d02_original_module_names_are_not_added_as_aliases():
    old_names = ["calc_hall_kier_alpha", "calc_hall_kier_alpha_with_contributions",
        "calc_kappa_1", "calc_kappa_2", "calc_kappa_3", "calc_phi", "calc_mqns", "calc_chi_0", "calc_chi_1",
        "calc_chi_nv", "calc_chi_nn"] + [f"calc_chi_{order}{family}" for family in ("v", "n") for order in range(5)]
    for name in old_names:
        assert not hasattr(ck, name)


def test_d02_generated_stub_query_parameters_and_exception_context():
    path = Path(__file__).resolve().parents[1] / "cosmolkit.pyi"
    classes = {n.name: n for n in ast.parse(path.read_text()).body if isinstance(n, ast.ClassDef)}
    methods = {n.name: n for n in classes["Molecule"].body if isinstance(n, ast.FunctionDef)}
    for name, parameter in PARAMETERS.items():
        args = methods[name].args
        assert [a.arg for a in args.args] == ["self"] + ([] if parameter is None else [parameter])
        assert not args.vararg and not args.kwarg and not args.kwonlyargs
        if parameter == "force":
            assert len(args.defaults) == 1 and ast.literal_eval(args.defaults[0]) is False
        else:
            assert not args.defaults
        if parameter is not None:
            assert args.args[-1].annotation is not None
            assert ast.unparse(args.args[-1].annotation) == ("builtins.bool" if parameter == "force" else "builtins.int")
    for name in ("DescriptorReadError", "DescriptorError"):
        cls = classes[name]
        assert [ast.unparse(x) for x in cls.bases] == ["builtins.ValueError"]
        fields = {cast(ast.Name, n.target).id: ast.unparse(n.annotation) for n in cls.body if isinstance(n, ast.AnnAssign)}
        assert fields["domain"] == fields["kind"] == "builtins.str"
        native = cast(object, getattr(ck, name))
        assert isinstance(native, type) and issubclass(native, ValueError)
    assert classes["DescriptorError"] is not classes["DescriptorReadError"]
