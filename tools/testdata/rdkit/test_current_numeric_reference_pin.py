"""Reference adapters must follow the shared pin, not a stale local version."""
import importlib
import json
from pathlib import Path

import pytest

ADAPTERS = [
    "_generate_conformer_generation_golden",
    "_generate_conformer_generation_library_golden",
    "_generate_forcefield_params_golden",
    "_generate_mmff_builtin_golden",
]


@pytest.mark.parametrize("name", ADAPTERS)
def test_numeric_adapter_uses_current_pin_and_rejects_wrong_runtime(name, monkeypatch):
    adapter = importlib.import_module(name)
    pin = json.loads(
        (Path(__file__).resolve().parents[3] / "testdata/reference/rdkit.json")
        .read_text(encoding="utf-8")
    )["python_distribution_version"]
    assert adapter.EXPECTED_RDKIT_VERSION == pin
    monkeypatch.setattr(adapter, "version", lambda _: pin)
    adapter.assert_rdkit_version()
    monkeypatch.setattr(adapter, "version", lambda _: "0.0.0")
    with pytest.raises(AssertionError, match="RDKit version mismatch"):
        adapter.assert_rdkit_version()


def test_conformer_recipes_keep_original_census_seed_and_attempt_budget():
    import _generate_conformer_generation_golden as fixed
    import _generate_conformer_generation_library_golden as library

    assert len(fixed.CASES) == 19
    assert len({case["case_id"] for case in fixed.CASES}) == 19
    assert library.CONFORMER_LIBRARY_SEED == 61453
    assert library.CONFORMER_LIBRARY_MAX_ITERATIONS == 3
    assert library.CONFORMER_LIBRARY_TIMEOUT == 0
    seeded_multi = next(case for case in fixed.CASES if case["case_id"] == "multi_etkdg_seeded")
    assert seeded_multi["num_confs"] == 10
    assert seeded_multi["attrs"] == {
        "randomSeed": 61453, "numThreads": 1, "trackFailures": True, "timeout": 1,
    }
