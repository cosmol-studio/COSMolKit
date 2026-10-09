"""Regression tests for canonical reusable Morgan fingerprints.

A false parameter selects an explicit connectivity(false) provider; the native
wrapper flag is recorded separately because that wrapper ignores its value.
"""
import ast
import hashlib
import json
from pathlib import Path
import pytest
import cosmolkit as ck

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "testdata/fingerprints/morgan_reusable_provider_ring_native.jsonl"
RECORDS = [json.loads(line) for line in FIXTURE.read_text().splitlines()]
assert len(RECORDS) == 26
FORMS = [
    ("bits", "fingerprint_morgan_with_generator", "fingerprints"),
    ("count", "fingerprint_morgan_count_with_generator", "counts"),
    ("sparse_bits", "fingerprint_morgan_sparse_with_generator", "sparse_fingerprints"),
    ("sparse_count", "fingerprint_morgan_sparse_count_with_generator", "sparse_counts"),
]

def allocated_output():
    output = ck.FingerprintAdditionalOutput()
    output.allocate_atom_counts()
    output.allocate_atom_to_bits()
    output.allocate_bit_info_map()
    output.allocate_bit_paths()
    return output

def same_value(actual, expected):
    if "nonzero" in expected:
        assert actual.length() == expected["size"]
        assert actual.nonzero_elements() == dict(expected["nonzero"])
    else:
        assert actual.n_bits() == expected["size"]
        assert actual.on_bits() == expected["native_bits"]

def same_output(actual, expected):
    assert actual.atom_counts() == expected["atom_counts"]
    assert actual.atom_to_bits() == expected["atom_to_bits"]
    assert actual.bit_info_map() == {key: [tuple(pair) for pair in pairs] for key, pairs in expected["bit_info_map"]}
    assert actual.bit_paths() == dict(expected["bit_paths"])
    assert actual.atoms_per_bit() is None

def native_generator(record):
    if record["case"] == "feature-invalid-source-null-pattern-skip":
        return ck.MorganFingerprintGenerator.from_json(json.dumps(record["input_json"]))
    profile = record["profile"]
    # Exact C++ wrapper records: all bare profiles capture connectivity(true).
    atom = None if profile.startswith("bare_") else ck.MorganAtomInvariantsGenerator.connectivity(profile == "explicit_true")
    return ck.MorganFingerprintGenerator(atom_invariants=atom)

@pytest.mark.parametrize("record", RECORDS, ids=[f"{r['case']}/{r.get('profile','restore')}/{r['smiles']}" for r in RECORDS])
@pytest.mark.parametrize("form,single,bulk", FORMS)
def test_source_native_full_scalar_ao_metadata_restore_and_bulk(record, form, single, bulk):
    molecule = ck.Molecule.from_smiles(record["smiles"])
    before = molecule.to_smiles()
    generator = native_generator(record)
    assert generator.info_string() == record["info"]
    assert json.loads(generator.to_json()) == record["json"]
    output = allocated_output()
    same_value(getattr(molecule, single)(generator, output=output), record["outputs"][form])
    same_output(output, record["outputs"][form]["additional_output"])
    restored = ck.MorganFingerprintGenerator.from_json(generator.to_json())
    same_value(getattr(molecule, single)(restored), record["outputs"][form])
    for threads in (1, 3, 7):
        values = getattr(generator, bulk)([molecule, None, molecule, None], num_threads=threads)
        assert len(values) == 4 and values[1] is None and values[3] is None
        same_value(values[0], record["outputs"][form])
        same_value(values[2], record["outputs"][form])
        assert getattr(generator, bulk)([], num_threads=threads) == []
        assert getattr(generator, bulk)([None, None], num_threads=threads) == [None, None]
    assert molecule.to_smiles() == before

def test_live_settings_all_fields_alias_snapshot_copy_and_lifetime():
    generator = ck.MorganFingerprintGenerator()
    settings, alias = generator.settings(), generator.settings()
    snapshot = settings.params()
    fields = {"radius": 0, "only_nonzero_invariants": True, "include_redundant_environments": True,
              "include_chirality": True, "count_simulation": True, "fp_size": 1000,
              "bits_per_feature": 2, "count_bounds": [1, 3, 5]}
    for name, value in fields.items():
        setattr(settings, name, value)
        assert getattr(alias, name) == value
    fields["count_bounds"].append(99)
    copied = alias.count_bounds
    copied.append(77)
    assert settings.count_bounds == [1, 3, 5]
    assert snapshot.radius == 3 and snapshot.fp_size == 2048 and snapshot.count_bounds == [1, 2, 4, 8]
    assert json.loads(generator.to_json())["bondInvariantsGenerator"]["useChirality"] == "false"
    del generator
    alias.radius = 1
    assert settings.radius == 1
    with pytest.raises(TypeError):
        ck.MorganSettings()

@pytest.mark.parametrize("name,value,error", [("radius", -1, OverflowError), ("fp_size", 4294967296, OverflowError), ("bits_per_feature", 1.5, TypeError), ("count_bounds", [-1], OverflowError)])
def test_invalid_binding_settings_transport_preserves_previous_value(name, value, error):
    settings = ck.MorganFingerprintGenerator().settings()
    old = getattr(settings, name)
    with pytest.raises(error):
        setattr(settings, name, value)
    assert getattr(settings, name) == old

def test_frozen_call_params_own_every_list_and_preserve_none_empty():
    supplied = [0]
    params = ck.MorganCallParams(from_atoms=supplied, ignore_atoms=supplied,
                               custom_atom_invariants=[11, 12, 13], custom_bond_invariants=[1, 1], conformer_id=3)
    supplied.append(1)
    assert params.from_atoms == [0] and params.ignore_atoms == [0]
    assert params.custom_atom_invariants == [11, 12, 13] and params.custom_bond_invariants == [1, 1]
    assert params.conformer_id == 3
    for name in ("from_atoms", "ignore_atoms", "custom_atom_invariants", "custom_bond_invariants"):
        snapshot = getattr(params, name)
        snapshot.append(99)
        assert getattr(params, name) != snapshot
        original = getattr(params, name)
        setattr(params, name, [])
        assert getattr(params, name) == []
        setattr(params, name, original)
    params.conformer_id = 0
    assert params.conformer_id == 0
    assert ck.MorganCallParams().from_atoms is None
    assert ck.MorganCallParams(from_atoms=[]).from_atoms == []

@pytest.mark.parametrize("form,single,bulk", FORMS)
def test_present_empty_roots_reset_all_allocated_output(form, single, bulk):
    molecule = ck.Molecule.from_smiles("CCO")
    generator = ck.MorganFingerprintGenerator()
    output = allocated_output()
    getattr(molecule, single)(generator, output=output)
    empty = getattr(molecule, single)(generator, params=ck.MorganCallParams(from_atoms=[]), output=output)
    assert empty.nonzero_elements() == {} if hasattr(empty, "nonzero_elements") else empty.on_bits() == []
    assert output.atom_counts() == [0, 0, 0] and output.atom_to_bits() == [[], [], []]
    assert output.bit_info_map() == {} and output.bit_paths() == {}

def test_provider_query_copy_lifetime_none_empty_and_custom_precedence():
    patterns = [ck.QueryGraph.from_smarts("[C]"), ck.QueryGraph.from_smarts("[O]")]
    provider = ck.MorganAtomInvariantsGenerator.features(patterns)
    patterns.clear()
    generator = ck.MorganFingerprintGenerator(params=ck.MorganParams(radius=0), atom_invariants=provider)
    del provider
    molecule = ck.Molecule.from_smiles("CCO")
    assert molecule.fingerprint_morgan_sparse_count_with_generator(generator).nonzero_elements() == {1: 2, 2: 1}
    assert json.loads(generator.to_json())["atomInvariantsGenerator"]["patternSMARTS"] == ["C", "O"]
    empty = ck.MorganFingerprintGenerator(params=ck.MorganParams(radius=0), atom_invariants=ck.MorganAtomInvariantsGenerator.features([]))
    assert molecule.fingerprint_morgan_sparse_count_with_generator(empty).nonzero_elements() == {0: 3}
    params = ck.MorganCallParams(custom_atom_invariants=[11, 12, 13])
    assert molecule.fingerprint_morgan_sparse_count_with_generator(empty, params=params).nonzero_elements() == {11: 1, 12: 1, 13: 1}
    default = ck.MorganFingerprintGenerator(atom_invariants=ck.MorganAtomInvariantsGenerator.features())
    assert "patternSMARTS" not in json.loads(default.to_json())["atomInvariantsGenerator"]

def test_root_explicit_ring_decision_and_bond_provider_precedence():
    params = ck.MorganParams(include_ring_membership=False, use_bond_types=False)
    generator = ck.MorganFingerprintGenerator(params=params)
    explicit = ck.MorganFingerprintGenerator(params=params, atom_invariants=ck.MorganAtomInvariantsGenerator.connectivity(True), bond_invariants=ck.MorganBondInvariantsGenerator(use_bond_types=True, include_chirality=True))
    assert json.loads(generator.to_json())["atomInvariantsGenerator"]["includeRingMembership"] == "false"
    value = json.loads(explicit.to_json())
    assert value["atomInvariantsGenerator"]["includeRingMembership"] == "true"
    assert value["bondInvariantsGenerator"] == {"type": "MorganBondInvGenerator", "useBondTypes": "true", "useChirality": "true"}
    assert value["fingerprintArguments"]["includeChirality"] == "false"

def test_source_null_provider_preconditions_json_fallback_and_live_error_recovery():
    molecule = ck.Molecule.from_smiles("CCO")
    value = json.loads(ck.MorganFingerprintGenerator().to_json())
    del value["atomInvariantsGenerator"], value["bondInvariantsGenerator"]
    null = ck.MorganFingerprintGenerator.from_json(json.dumps(value))
    with pytest.raises(ck.MorganReadError, match="atom invariants") as caught:
        molecule.fingerprint_morgan_sparse_count_with_generator(null)
    assert caught.value.kind == "Generator" and caught.value.__cause__ is not None
    with pytest.raises(ck.MorganReadError, match="bond invariants"):
        molecule.fingerprint_morgan_sparse_count_with_generator(null, params=ck.MorganCallParams(custom_atom_invariants=[11, 12, 13]))
    assert molecule.fingerprint_morgan_sparse_count_with_generator(null, params=ck.MorganCallParams(custom_atom_invariants=[11, 12, 13], custom_bond_invariants=[1, 1])).nonzero_elements()
    generator = ck.MorganFingerprintGenerator()
    settings = generator.settings()
    settings.count_simulation = True
    settings.count_bounds = []
    with pytest.raises(ck.MorganReadError) as caught:
        generator.fingerprints([None, molecule], num_threads=3)
    assert caught.value.kind == "Generator" and caught.value.__cause__ is not None
    assert generator.sparse_counts([molecule], num_threads=3)[0].nonzero_elements()
    settings.count_bounds = [1, 2, 4, 8]
    assert generator.fingerprints([molecule], num_threads=3)[0].on_bits()
    with pytest.raises(ck.MorganReadError, match="INT_MIN"):
        generator.counts([molecule], num_threads=-2147483648)
    text = generator.to_json().replace('"fpSize":"2048"', '"fpSize":1,"fpSize":2')
    assert ck.MorganFingerprintGenerator.from_json(text).settings().fp_size == 1

def test_original_stub_projects_all_five_types_methods_and_settings():
    tree = ast.parse((ROOT / "python/cosmolkit.pyi").read_text())
    classes = {node.name: node for node in tree.body if isinstance(node, ast.ClassDef)}
    for name in ("MorganFingerprintGenerator", "MorganSettings", "MorganCallParams", "MorganAtomInvariantsGenerator", "MorganBondInvariantsGenerator"):
        assert name in classes
    methods = {node.name for node in classes["MorganFingerprintGenerator"].body if isinstance(node, ast.FunctionDef)}
    assert methods == {"__new__", "new", "from_json", "settings", "info_string", "to_json", "__repr__", "fingerprints", "counts", "sparse_fingerprints", "sparse_counts"}
    names = [node.name for node in classes["MorganSettings"].body if isinstance(node, ast.FunctionDef)]
    for field in ("radius", "only_nonzero_invariants", "include_redundant_environments", "include_chirality", "count_simulation", "fp_size", "bits_per_feature", "count_bounds"):
        assert names.count(field) == 2
    molecule = {node.name for node in classes["Molecule"].body if isinstance(node, ast.FunctionDef)}
    assert all(single in molecule for _, single, _ in FORMS)
