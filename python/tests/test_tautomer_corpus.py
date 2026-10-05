"""Public Python projection of the prepared Rust TAU corpus tasks.

Preparation is separate: tests consume authenticated native records and the
persisted Rust profiles. Missing, stale, failed or incomplete records fail.
"""
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import subprocess
import time
import traceback

import cosmolkit as ck
import pytest


ROOT = Path(__file__).resolve().parents[2]
TASKS = ("tautomer_enumeration_smiles", "tautomer_canonicalization_smiles")


def _digest(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@pytest.fixture(scope="module")
def prepared_tautomers():
    data = Path(os.environ.get("PARITY_DATA", ROOT / "target/parity-tests-tautomer5000"))
    binary = ROOT / "target/release/cosmolkit-parity-tests"
    # The existing registry owns source identity, matrix and global preflight.
    # This read-only command executes zero COSMolKit chemistry operations.
    result = subprocess.run([str(binary), "preflight", "--data", str(data)],
                            cwd=ROOT, text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Ready: 20000 cases; 0 Rust operation calls" in result.stdout
    suite = json.loads((data / "suite.json").read_text())
    assert suite["tasks"] == list(TASKS)
    prepared = {}
    for task in TASKS:
        matches = []
        for path in data.glob(f"{task}-*/manifest.json"):
            manifest = json.loads(path.read_text())
            # Suite may retain older immutable generations. Only a generation
            # with this build's registry identity can be selected.
            matches.append((path.parent, manifest))
        assert matches, f"missing prepared task {task}"
        # Current generation selection follows the existing registry's suite
        # selection: input identity is constant, stale references are rejected
        # by global preflight above. All retained candidates must agree on input.
        inputs = []
        for directory, manifest in matches:
            assert _digest(directory / "input.json") == manifest["input_sha256"]
            inputs.append(json.loads((directory / "input.json").read_text()))
        assert all(value == inputs[0] for value in inputs)
        directory, manifest = matches[-1]
        assert manifest["rows"] == len(inputs[0]) == 10000
        identity = manifest["imported_reference"]
        assert identity["records"] == 5000 and identity["branches"] == 10000
        native = ROOT / identity["path"]
        assert _digest(native) == identity["sha256"]
        assert _digest(ROOT / identity["source_identity_manifest"]) == identity["source_identity_manifest_sha256"]
        assert identity["source_revision"] == "351f8f378f8ad6bbd517980c38896e66bf907af8"
        prepared[task] = (inputs[0], native)
    output = data / "python-results" / f"run-{time.time_ns()}"
    output.mkdir(parents=True)
    return output, prepared


def _state(molecule):
    atoms = []
    for atom in molecule.atoms():
        descriptor = atom.cip_descriptor()
        atoms.append(dict(atomic_number=atom.atomic_number(), formal_charge=atom.formal_charge(),
                          explicit_hydrogens=atom.explicit_hydrogens(), no_implicit=atom.no_implicit(),
                          isotope=atom.isotope() if atom.isotope() is not None else 0,
                          radical_electrons=atom.radical_electrons(), aromatic=atom.is_aromatic(),
                          chiral_tag=atom.chiral_tag_name(), hybridization=atom.hybridization().name,
                          cip_code=None if descriptor is None else descriptor.value))
    bonds = [dict(begin=bond.begin(), end=bond.end(), bond_type=bond.order_name(),
                  aromatic=bond.is_aromatic(), conjugated=bond.is_conjugated(),
                  direction=bond.direction_name(), stereo=bond.stereo_name().removeprefix("STEREO"),
                  stereo_atoms=list(bond.stereo_atoms()) if bond.stereo_atoms() is not None else [])
             for bond in molecule.bonds()]
    return dict(isomeric_smiles=molecule.to_smiles(), atoms=atoms, bonds=bonds)


def _score(molecule):
    score = molecule.tautomer_score()
    return dict(ring=score.ring(), substructure=score.substructure(),
                hetero_hydrogen=score.hetero_hydrogen(), total=score.total())


def _status(value):
    # This PyO3 enum exposes its variants directly rather than enum.Enum.name.
    for name in ("Completed", "MaxTautomersReached", "MaxTransformsReached", "Canceled"):
        if value == getattr(ck.TautomerEnumerationStatus, name):
            return name
    raise AssertionError(f"unknown enumeration status {value!r}")


def _parameters(profile):
    # No Python matrix: these values come from serialized Rust task inputs.
    params = ck.TautomerParams() if profile["catalog"] == "current" else ck.TautomerParams.v1()
    for field, value in profile.items():
        if field != "catalog":
            getattr(params, "set_" + field)(value)
    return params


def _difference(expected, actual, path="result"):
    if type(expected) is not type(actual):
        return f"{path}: type {type(expected).__name__} != {type(actual).__name__}"
    if isinstance(expected, dict):
        if expected.keys() != actual.keys():
            return f"{path}: keys differ"
        for key in expected:
            difference = _difference(expected[key], actual[key], f"{path}.{key}")
            if difference:
                return difference
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            return f"{path}: length {len(expected)} != {len(actual)}"
        for index, (left, right) in enumerate(zip(expected, actual)):
            difference = _difference(left, right, f"{path}[{index}]")
            if difference:
                return difference
    elif expected != actual:
        return f"{path}: {expected!r} != {actual!r}"
    return None


def _evaluate(task, row, inputs, native, failure_dir):
    failures = []
    for item in inputs:
        observation = item["Molecular"]
        case = observation["case"]
        profile = next(iter(observation["profile"].values()))["parameters"]
        branch = "default" if profile["catalog"] == "current" else "v1"
        actual = None
        try:
            assert native["row"] == row and native["smiles"] == case["smiles"]
            assert case["id"] == f"line:{row + 1}"
            assert native["parse"]["ok"] and native["parse"]["error"] is None
            expected = native["branches"][branch]
            assert expected["ok"] and expected["error"] is None
            parameters = dict(expected["parameters"])
            assert parameters.pop("name") == branch
            assert parameters == profile
            # Same canonical parser defaults as the Rust molecular executor.
            source = ck.Molecule.from_smiles(case["smiles"])
            before, props_before = _state(source), source.properties().props()
            params = _parameters(profile)
            if task == "tautomer_enumeration_smiles":
                result = source.enumerate_tautomers_with_params(params)
                canonical = result.canonical_tautomer()
                actual = dict(ordered_smiles=result.canonical_smiles(), status=_status(result.status()),
                              modified_atoms=result.modified_atoms(), modified_bonds=result.modified_bonds(),
                              scores=[_score(value) for value in result],
                              molecule_states=[_state(value) for value in result],
                              canonical_smiles=canonical.to_smiles(), canonical_state=_state(canonical))
            else:
                canonical = source.canonical_tautomer_with_params(params)
                actual = dict(canonical_smiles=canonical.to_smiles(), canonical_state=_state(canonical),
                              canonical_score=_score(canonical))
            assert _state(source) == before and source.properties().props() == props_before
            chemical = {key: value for key, value in expected.items() if key not in ("parameters", "ok", "error")}
            difference = _difference(chemical, actual)
            assert difference is None, difference
        except Exception:
            failure = dict(row=row, case=case, profile=profile, branch=branch,
                           error=traceback.format_exc(), actual=actual)
            path = Path(failure_dir) / f"{task}-row-{row}-{branch}.json"
            path.write_text(json.dumps(failure, indent=2) + "\n")
            failures.append(dict(row=row, branch=branch, error=failure["error"], details=str(path)))
    return len(inputs), failures


@pytest.mark.parametrize("task", TASKS)
def test_tautomer_public_corpus_5000(prepared_tautomers, task):
    output, prepared = prepared_tautomers
    inputs, native = prepared[task]
    failures, branches, rows = [], 0, 0
    started = time.time()
    # Bound retained native rows; never load the 11 GB enumeration corpus at
    # once. Multiprocessing exercises only the compiled public Python binding.
    with ProcessPoolExecutor(max_workers=8, mp_context=multiprocessing.get_context("spawn")) as pool:
        pending = set()
        with native.open() as records:
            for row, line in enumerate(records):
                assert row < 5000
                pair = inputs[row * 2:row * 2 + 2]
                assert len(pair) == 2
                pending.add(pool.submit(_evaluate, task, row, pair, json.loads(line), str(output)))
                rows += 1
                if len(pending) >= 16:
                    done, pending = wait(pending, return_when=FIRST_COMPLETED)
                    for future in done:
                        count, failed = future.result()
                        branches += count
                        failures.extend(failed)
            while pending:
                done, pending = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    count, failed = future.result()
                    branches += count
                    failures.extend(failed)
    receipt = dict(task=task, rows=rows, branches=branches, failures=failures,
                   started=started, finished=time.time(), reference_sha256=_digest(native),
                   fields="all ordered keys, molecule atom/bond fields, status, modified sets, scores and canonical state")
    path = output / f"{task}.json"
    path.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"{task}: {rows} rows, {branches} branches, {len(failures)} failures; {path}")
    assert rows == 5000 and branches == 10000, receipt
    assert not failures, f"{len(failures)} failed branches; full evidence: {path}; first: {failures[:1]}"
