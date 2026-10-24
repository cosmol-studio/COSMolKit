"""Fixed examples of binding mechanisms, not a growing business-API suite.

Rust/PyO3 and the shared projection framework own their guarantees. These few
samples exercise the Python boundary; do not add cases for ordinary new APIs,
enumerate registry entries, or move chemistry tests here.
"""
from pathlib import Path
from operator import setitem
import sys
from tempfile import TemporaryDirectory

import numpy as np

sys.path.insert(0, str(Path(sys.argv[1]).resolve()))
import cosmolkit as ck
from cosmolkit import _native


def raises(kind, call):
    try:
        call()
    except kind as error:
        return error
    raise AssertionError(f"expected {kind.__name__}")


def native_enum_protocol():
    member = ck.BondOrder.SINGLE
    assert member.name == "Single" and repr(member) == "BondOrder.Single"
    assert ck.BondOrder(int(member)) == member and hash(member) == hash(int(member))
    assert {member: "bond"}[int(member)] == "bond"
    assert ck.BondOrder.names["SINGLE"] == ck.BondOrder.names["Single"] == member
    assert ck.BondOrder.values[int(member)] == member
    assert ck.BondOrder.names is ck.BondOrder.names
    assert len(ck.BondOrder.values) < len(ck.BondOrder.names)
    raises(TypeError, lambda: setitem(ck.BondOrder.names, "bad", member))
    raises(ValueError, lambda: ck.BondOrder(999))


def forwarding_and_configuration():
    assert ck.Molecule is _native.Molecule
    assert ck.mol_from_smiles is _native.mol_from_smiles
    molecule = ck.mol_from_smiles("CCO")
    params = ck.MaccsFingerprintParams()
    params.n_bits = 167
    assert params.n_bits == 167
    raises(TypeError, lambda: setattr(params, "n_bits", "wrong"))
    assert params.n_bits == 167
    raises(AttributeError, lambda: setattr(params, "unknown", 1))
    params.n_bits = 166
    expected = molecule.fingerprint_maccs(params).on_bits
    assert molecule.fingerprint_maccs(n_bits=166).on_bits == expected
    assert molecule.fingerprint_maccs(None).on_bits == expected
    assert molecule.fingerprint_maccs(n_bits=None).on_bits == expected
    raises(TypeError, lambda: molecule.fingerprint_maccs(params, n_bits=None))
    raises(TypeError, lambda: molecule.fingerprint_maccs(None, n_bits=166))

    source = {"{A}": "O"}
    options = ck.SmartsParseParams(replacements=source)
    source["{A}"] = "C"
    assert options.replacements == {"{A}": "O"}
    raises(AttributeError, lambda: options.replacements.update({"{A}": "N"}))
    options.replacements = {"{A}": "N"}
    raises(TypeError, lambda: setattr(options, "replacements", {1: "O"}))
    assert options.replacements == {"{A}": "N"}

    child = ck.McsAtomCompareParameters()
    options = ck.McsParameters(atom_compare_parameters=child)
    child.match_valences = not child.match_valences
    assert options.atom_compare_parameters.match_valences != child.match_valences
    options.atom_compare_parameters.match_valences = child.match_valences
    assert options.atom_compare_parameters.match_valences == child.match_valences
    view = options.atom_compare_parameters
    options.atom_compare_parameters = ck.McsAtomCompareParameters()
    assert view.match_valences is False
    view.match_valences = True
    assert options.atom_compare_parameters.match_valences is True
    options.atom_compare_parameters = options.atom_compare_parameters
    del options
    view.match_valences = False
    assert view.match_valences is False

    source = [0]
    options = ck.MorganCallParams(from_atoms=source)
    source.append(1)
    assert options.from_atoms == (0,)
    raises(AttributeError, lambda: options.from_atoms.append(2))
    raises(TypeError, lambda: setattr(options, "from_atoms", ["bad"]))
    assert options.from_atoms == (0,)
    options.from_atoms = [0, 2]
    assert options.from_atoms == (0, 2)

    selection = ck.SdfReadParams(coordinate_mode="require_3d")
    assert selection.coordinate_mode == ck.SdfCoordinateMode.Require3D
    selection.coordinate_mode = ck.SdfCoordinateMode.Preserve
    raises(ValueError, lambda: setattr(selection, "coordinate_mode", "unknown"))
    assert selection.coordinate_mode == ck.SdfCoordinateMode.Preserve


def indexing_and_numpy():
    fingerprint = ck.Fingerprint.from_on_bits(8, [1, 4])
    items = [fingerprint, None, fingerprint]
    batch = ck.FingerprintBatch(items)
    assert batch[0] is fingerprint and batch[-1] is fingerprint and batch[1] is None
    for selection in (slice(None, None, -1), slice(1, None, sys.maxsize),
                      slice(None, None, -(10**100)), slice(0, 0)):
        expected, actual = items[selection], batch[selection]
        assert isinstance(actual, ck.FingerprintBatch) and len(actual) == len(expected)
        assert all(actual[i] is item for i, item in enumerate(expected))
    raises(IndexError, lambda: batch[9])
    raises(TypeError, lambda: batch["wrong"])
    raises(ValueError, lambda: batch[::0])
    assert batch[[2, 0, 2]][0] is fingerprint and len(batch[[]]) == 0
    assert len(batch[[True, False, True]]) == 2
    raises(TypeError, lambda: batch[True])
    raises(TypeError, lambda: batch[[True, 0]])
    raises(ValueError, lambda: batch[[True]])
    raises(IndexError, lambda: batch[10**100])
    with TemporaryDirectory() as folder:
        path = Path(folder) / "index.sdf"
        path.write_text(ck.mol_from_smiles("CO").to_sdf())
        dataset = ck.SdfDataset.open(path)
        assert isinstance(dataset[-1], ck.SdfRecord)
        assert isinstance(dataset[:1], ck.MoleculeBatch)

    array = fingerprint.to_numpy()
    assert array.dtype == np.uint8 and array.shape == (8,)
    assert np.flatnonzero(array).tolist() == fingerprint.on_bits
    array[:] = 0
    assert fingerprint.on_bits == [1, 4]
    matrix = ck.FingerprintBatch([fingerprint, fingerprint]).to_numpy()
    assert matrix.dtype == np.uint8 and matrix.shape == (2, 8)
    assert matrix[0].tolist() == matrix[1].tolist() == fingerprint.to_numpy().tolist()


def structured_errors():
    molecule = ck.mol_from_smiles("CCO")
    leaf = raises(ck.MatrixError, molecule.distance_matrix_3d)
    assert leaf.variant == ck.MatrixErrorKind.No3dConformer and leaf.__cause__ is None
    fingerprint = ck.Fingerprint.from_on_bits(8, [1])
    payload = raises(ck.FingerprintError, lambda: fingerprint.get_bit(9))
    assert payload.variant == ck.FingerprintErrorKind.SparseIndexOutOfRange
    assert payload.index == 9 and payload.size == 8
    raises(AttributeError, lambda: setattr(payload, "index", 0))
    error = raises(ck.SmilesError, lambda: ck.mol_from_smiles("C(C"))
    cause = error.source
    assert isinstance(cause, ck.SmilesParseError) and cause is error.__cause__
    assert cause.variant == ck.SmilesParseErrorKind.Syntax
    assert cause.offset == 3 and cause.message == "unclosed branch"
    assert str(cause) == "invalid SMILES syntax at byte 3: unclosed branch"
    error.__cause__ = None
    assert error.source is None
    error.__cause__ = cause
    assert error.source is cause
    raises(TypeError, ck.SmilesError)
    record = raises(ck.PropertyValueError, ck.PropertyValue.Int(7).as_bool)
    assert record.expected == ck.PropertyValueKind.Bool
    assert record.actual == ck.PropertyValueKind.Int and not hasattr(record, "variant")


def value_ownership_and_callbacks():
    molecule = ck.mol_from_smiles("CCO")
    copied = molecule.copy()
    copied.set_atom_property_int_(0, "tracking_id", 42)
    assert copied.atom_property_int(0, "tracking_id") == 42
    assert molecule.atom_property_int(0, "tracking_id") is None

    retained = []
    def accept(target, indices):
        target.set_atom_property_bool_(indices[0], "callback", True)
        retained.append(target)
        return True
    query = ck.parse_smarts("O")
    assert molecule.has_substruct_match(query, ck.SubstructMatchParams(final_match=accept))
    assert retained and molecule.atom_property_bool(2, "callback") is None
    assert retained[0].atom_property_bool(2, "callback")
    original_error = RuntimeError("callback failure")
    def reject(target, indices):
        raise original_error
    assert raises(RuntimeError, lambda: molecule.has_substruct_match(
        query, ck.SubstructMatchParams(final_match=reject))) is original_error

    live, owned = [], []
    def snapshot(source, progress):
        live.append(source)
        owned.append(source.snapshot())
        return False
    source = ck.mol_from_smiles("CC=O")
    source.enumerate_tautomers(callback=snapshot)
    assert live and owned[0].num_atoms == source.num_atoms
    raises(ReferenceError, lambda: live[0].num_atoms)
    seed = (1 << 20000) + 17
    assert (ck.StereoisomerRandomSource.from_integer_seed(seed)
            == ck.StereoisomerRandomSource.from_integer_seed(seed))


def consuming_resource():
    with TemporaryDirectory() as folder:
        path = Path(folder) / "source.sdf"
        path.write_text(ck.mol_from_smiles("C").to_sdf(include_coordinates=False))
        stream = ck.SdfRecordStream.open(path)
        batches = stream.batches(1, n_jobs=1)
        raises(RuntimeError, stream.next_record)
        del stream
        assert len(batches.next_batch()) == 1
        assert batches.next_batch() is None


native_enum_protocol()
forwarding_and_configuration()
indexing_and_numpy()
structured_errors()
value_ownership_and_callbacks()
consuming_resource()
print("Fixed binding protocol samples: OK")
