"""Fixed tests of installed Python conversions, options, errors and ownership.

Morgan literals come from the installed Rust facade's pinned source regressions,
not a live oracle. This is not a chemistry corpus or whole 0.3.0 acceptance.
"""

import ast
import inspect
import struct
import tomllib
from collections.abc import Callable
from pathlib import Path
from typing import cast

import cosmolkit as ck
import pytest


def call(value: object) -> Callable[..., object]:
    assert callable(value)
    return value


def state(mol: ck.Molecule) -> tuple[str, int, int, list[list[float]] | None]:
    return mol.to_smiles(), mol.num_atoms(), mol.num_bonds(), mol.coordinates_2d()


def test_real_default_module_and_one_class():
    assert ck._binding_profile == "canonical-bootstrap"
    manifest = tomllib.loads((Path(__file__).resolve().parents[1] / "Cargo.toml").read_text())
    assert ck.version() == ck.__version__ == manifest["package"]["version"]
    assert ck.Molecule.__module__ == "cosmolkit"
    assert state(ck.Molecule.new()) == ("", 0, 0, None)
    assert type(ck.Molecule.from_smiles("CCO")) is ck.Molecule
    assert type(ck.Molecule.from_smiles("CCO").with_2d_coordinates()) is ck.Molecule


@pytest.mark.parametrize("input,output,atoms,bonds", [
    ("", "", 0, 0), ("OCC", "CCO", 3, 2),
    ("c1ccccc1", "c1ccccc1", 6, 6), ("[13CH4]", "[13CH4]", 1, 0),
])
def test_smiles_constructor_and_default_writer(input: str, output: str, atoms: int, bonds: int):
    mol = ck.Molecule.from_smiles(input)
    explicit = ck.Molecule.from_smiles_with_params(input, ck.SmilesParseParams())
    expected = (output, atoms, bonds, None)
    assert state(mol) == state(explicit) == expected
    assert mol.to_smiles_with_params(ck.SmilesWriteParams()) == output
    assert state(mol) == expected


PARSE_DEFAULTS: dict[str, object] = dict(sanitize=True, allow_cxsmiles=True,
    strict_cxsmiles=True, parse_name=True, remove_hs=True,
    skip_cleanup=False, debug_parse=False, replacements={})
WRITE_DEFAULTS: dict[str, object] = dict(isomeric_smiles=True, kekule=False,
    canonical=True, clean_stereo=True, rooted_at_atom=None, all_bonds_explicit=False,
    all_hydrogens_explicit=False, include_dative_bonds=True, ignore_atom_map_numbers=False)


@pytest.mark.parametrize("factory,defaults", [
    (ck.SmilesParseParams, PARSE_DEFAULTS), (ck.SmilesWriteParams, WRITE_DEFAULTS),
])
def test_option_defaults_explicit_fields_and_writable_configuration(factory: object, defaults: dict[str, object]):
    params = call(factory)()
    assert {key: cast(object, getattr(params, key)) for key in defaults} == defaults
    explicit = {key: not val if type(val) is bool else val for key, val in defaults.items()}
    configured = call(factory)(**explicit)
    assert {key: cast(object, getattr(configured, key)) for key in explicit} == explicit
    for key in defaults:
        setattr(params, key, explicit[key])
        assert getattr(params, key) == explicit[key]
    assert {key: getattr(configured, key) for key in explicit} == explicit
    with pytest.raises(TypeError):
        _ = call(factory)(unknown_option=True)
    with pytest.raises(TypeError):
        _ = call(factory)(False)


def test_replacements_are_copied_and_forwarded():
    mapping = {"{X}": "O"}
    params = ck.SmilesParseParams(replacements=mapping)
    mapping["{X}"] = "N"
    returned = params.replacements
    returned["{X}"] = "F"
    assert params.replacements == {"{X}": "O"}
    mol = ck.Molecule.from_smiles_with_params("CC{X}", params)
    assert mol.to_smiles() == "CCO"


@pytest.mark.parametrize("remove_hs,expected_atoms", [(True, 1), (False, 5)])
def test_hydrogen_option_forwarding(remove_hs: bool, expected_atoms: int):
    mol = ck.Molecule.from_smiles_with_params("[H]C([H])([H])[H]",
        ck.SmilesParseParams(remove_hs=remove_hs))
    assert mol.num_atoms() == expected_atoms


@pytest.mark.parametrize("family,error_type,label", [
    ("morgan", ck.MorganReadError, "Morgan"),
    ("atom_pair", ck.AtomPairReadError, "AtomPair"),
    ("topological_torsion", ck.TopologicalTorsionReadError, "Topological Torsion"),
])
@pytest.mark.parametrize("suffix", [
    "", "_sparse", "_count", "_sparse_count",
])
def test_unsanitized_read_reports_missing_preparation(family: str, error_type: type[Exception], label: str, suffix: str):
    mol = ck.Molecule.from_smiles_with_params("CCO",
        ck.SmilesParseParams(sanitize=False, remove_hs=False))
    before = state(mol)
    with pytest.raises(error_type) as caught:
        _ = call(getattr(mol, f"fingerprint_{family}{suffix}"))()
    assert (caught.value.domain, caught.value.kind) == ("fingerprints", "Preparation")
    preparation = caught.value.__cause__
    assert isinstance(preparation, ck.FingerprintPreparationError)
    assert (preparation.domain, preparation.kind) == ("fingerprints", "MissingPreparedValence")
    assert str(preparation) == "Fingerprint preparation requires a valid prepared valence assignment"
    assert str(caught.value) == f"{label} preparation failed: {preparation}"
    assert preparation.__cause__ is None
    assert state(mol) == before


@pytest.mark.parametrize("input,params,fragment", [
    ("(", ck.SmilesParseParams(), "invalid SMILES syntax"),
    ("C1", ck.SmilesParseParams(), "unclosed ring index 1"),
    ("C", ck.SmilesParseParams(replacements={"": "C"}), "replacement key must not be empty"),
    ("{X}", ck.SmilesParseParams(replacements={"{X}": "{X}"}), "replacements do not converge"),
])
def test_constructor_errors_keep_category_and_real_cause(input: str, params: ck.SmilesParseParams, fragment: str):
    with pytest.raises(ck.SmilesError) as caught:
        _ = ck.Molecule.from_smiles_with_params(input, params)
    error = caught.value
    assert isinstance(error, ValueError)
    assert (error.domain, error.kind) == ("smiles", "Parse")
    assert fragment in str(error)
    assert isinstance(error.__cause__, ValueError)
    assert str(error) == "SMILES parsing failed: " + str(error.__cause__)


def test_writer_options_are_forwarded_and_root_errors_keep_cause():
    mol = ck.Molecule.from_smiles("CCO")
    before = state(mol)
    assert mol.to_smiles_with_params(ck.SmilesWriteParams(all_bonds_explicit=True)) == "C-C-O"
    assert mol.to_smiles_with_params(ck.SmilesWriteParams(rooted_at_atom=2)) == "OCC"
    chiral = ck.Molecule.from_smiles("F[C@H](Cl)Br")
    assert "@" in chiral.to_smiles()
    assert "@" not in chiral.to_smiles_with_params(ck.SmilesWriteParams(isomeric_smiles=False))
    with pytest.raises(ck.SmilesWriteError) as caught:
        _ = mol.to_smiles_with_params(ck.SmilesWriteParams(rooted_at_atom=99))
    error = caught.value
    assert (error.domain, error.kind) == ("smiles", "Write")
    assert str(error) == "SMILES serialization failed: root atom index 99 is out of range for 3 atoms"
    assert isinstance(error.__cause__, ValueError)
    assert str(error.__cause__) == "root atom index 99 is out of range for 3 atoms"
    assert state(mol) == before


@pytest.mark.parametrize("value", [-1, 2 ** (8 * struct.calcsize("P"))])
def test_root_usize_conversion(value: int):
    with pytest.raises(OverflowError):
        _ = ck.SmilesWriteParams(rooted_at_atom=value)


@pytest.mark.parametrize("method,result_type,width,values", [
    ("fingerprint_morgan", ck.Fingerprint, 2048, [80, 222, 294, 807, 1057, 1410]),
    ("fingerprint_morgan_sparse", ck.SparseBitFingerprint, 2**32-1,
        [-2049583024, -2048238559, -752510682, -276918910, 864662311, 1535166686]),
    ("fingerprint_morgan_count", ck.SparseCountFingerprint32, 2048,
        {80: 1, 222: 1, 294: 1, 807: 1, 1057: 1, 1410: 1}),
    ("fingerprint_morgan_sparse_count", ck.SparseCountFingerprint, 2**64-1,
        {864662311: 1, 1535166686: 1, 2245384272: 1, 2246728737: 1, 3542456614: 1, 4018048386: 1}),
])
def test_four_morgan_results_and_owned_containers(method: str, result_type: type, width: int, values: object):
    mol = ck.Molecule.from_smiles("CCO").with_2d_coordinates()
    before = state(mol)
    result = call(cast(object, getattr(mol, method)))()
    assert type(result) is result_type
    is_bits = method in ("fingerprint_morgan", "fingerprint_morgan_sparse")
    accessor = "on_bits" if is_bits else "nonzero_elements"
    width_accessor = "n_bits" if is_bits else "length"
    assert call(cast(object, getattr(result, width_accessor)))() == width
    returned = call(cast(object, getattr(result, accessor)))()
    assert returned == values
    assert str(result).startswith(result_type.__name__ + "(")
    if is_bits:
        assert len(cast(ck.Fingerprint, result)) == width
        cast(list[int], returned).clear()
    else:
        cast(dict[int, int], returned).clear()
        assert call(cast(object, getattr(result, "total_value")))() == 6
    assert call(cast(object, getattr(result, accessor)))() == values
    assert state(mol) == before
    empty = call(cast(object, getattr(ck.Molecule.from_smiles(""), method)))()
    assert call(cast(object, getattr(empty, accessor)))() == ([] if is_bits else {})


@pytest.mark.parametrize("factory,bits", [(ck.SparseCountFingerprint, 64), (ck.SparseCountFingerprint32, 32)])
def test_count_widths_values_errors_and_copies(factory: object, bits: int):
    result = cast(ck.SparseCountFingerprint, call(cast(object, getattr(factory, "new")))((1 << bits)-1))
    result.set_value((1 << bits)-1, -7)
    assert result.value((1 << bits)-1) == -7
    assert result.nonzero_elements() == {(1 << bits)-1: -7}
    result.set_value((1 << bits)-1, 0)
    assert result.nonzero_elements() == {}
    for invalid in [-1, (1 << bits)]:
        with pytest.raises(OverflowError):
            _ = call(result.set_value)(invalid, 1)
    with pytest.raises(OverflowError):
        _ = result.set_value(1, 2**31)
    with pytest.raises(OverflowError):
        _ = result.set_value(1, -2**31-1)
    small = cast(ck.SparseCountFingerprint, call(cast(object, getattr(factory, "new")))(3))
    with pytest.raises(ck.FingerprintError) as caught:
        _ = small.value(3)
    assert (caught.value.kind, caught.value.index, caught.value.size) == ("SparseIndexOutOfRange", 3, 3)
    assert caught.value.domain == "fingerprints" and caught.value.__cause__ is None


@pytest.mark.parametrize("factory", [ck.SparseCountFingerprint, ck.SparseCountFingerprint32])
def test_complete_count_arithmetic_delegates_and_preserves_inputs(factory: object):
    a = cast(ck.SparseCountFingerprint, call(cast(object, getattr(factory, "new")))(10))
    b = cast(ck.SparseCountFingerprint, call(cast(object, getattr(factory, "new")))(10))
    a.set_value(1, 5); a.set_value(3, -2); a.set_value(8, 4)
    b.set_value(1, 3); b.set_value(3, -4); b.set_value(9, 7)
    before = a.nonzero_elements(), b.nonzero_elements()
    assert a.fuzzy_and(b).nonzero_elements() == {1: 3, 3: -4}
    assert a.fuzzy_or(b).nonzero_elements() == {1: 5, 3: -2, 8: 4, 9: 7}
    assert a.with_added(b).nonzero_elements() == {1: 8, 3: -6, 8: 4, 9: 7}
    assert a.with_subtracted(b).nonzero_elements() == {1: 2, 3: 2, 8: 4, 9: -7}
    assert a.with_added_scalar(2).nonzero_elements() == {1: 7, 3: 0, 8: 6}
    assert a.with_subtracted_scalar(2).nonzero_elements() == {1: 3, 3: -4, 8: 2}
    assert a.with_multiplied_scalar(2).nonzero_elements() == {1: 10, 3: -4, 8: 8}
    assert a.with_divided_scalar(2).nonzero_elements() == {1: 2, 3: -1, 8: 2}
    assert a.total_value() == 7 and a.total_value(use_abs=True) == 11
    with pytest.raises(ck.FingerprintError) as caught:
        _ = a.with_divided_scalar(0)
    assert caught.value.kind == "UndefinedArithmetic"
    assert caught.value.site == "SparseIntVect::operator/=(int)"
    assert a.nonzero_elements() == before[0] and b.nonzero_elements() == before[1]
    empty = cast(ck.SparseCountFingerprint, call(cast(object, getattr(factory, "new")))(10))
    assert empty.with_divided_scalar(0).nonzero_elements() == {}
    mismatch = cast(ck.SparseCountFingerprint, call(cast(object, getattr(factory, "new")))(11))
    for method in [a.fuzzy_and, a.fuzzy_or, a.with_added, a.with_subtracted]:
        with pytest.raises(ck.FingerprintError) as caught:
            _ = method(mismatch)
        assert (caught.value.kind, caught.value.left, caught.value.right) == ("BitLengthMismatch", 10, 11)


def test_generated_stubs_match_installed_canonical_surface():
    stub = Path(__file__).resolve().parents[1] / "cosmolkit.pyi"
    tree = ast.parse(stub.read_text())
    classes = {n.name: n for n in tree.body if isinstance(n, ast.ClassDef)}
    assert {"Molecule", "SmilesParseParams", "SmilesWriteParams", "Fingerprint", "SparseBitFingerprint",
            "SparseCountFingerprint", "SparseCountFingerprint32", "SmilesError", "SmilesWriteError", "MorganReadError", "FingerprintPreparationError", "FingerprintError"} <= classes.keys()
    for name in ("SmilesError", "SmilesWriteError", "MorganReadError", "FingerprintPreparationError", "FingerprintError"):
        assert [ast.unparse(n) for n in classes[name].bases] == ["builtins.ValueError"]
    for name, defaults in [("SmilesParseParams", PARSE_DEFAULTS), ("SmilesWriteParams", WRITE_DEFAULTS)]:
        factory = cast(object, getattr(ck, name))
        signature = inspect.signature(call(factory))
        expected = dict(defaults)
        if name == "SmilesParseParams": expected["replacements"] = None
        assert {key: cast(object, p.default) for key, p in signature.parameters.items()} == expected
    signature = inspect.signature(ck.Molecule.from_smiles_with_params)
    assert list(signature.parameters) == ["input", "params"]
    assert all(cast(object, p.default) is inspect.Parameter.empty for p in signature.parameters.values())
