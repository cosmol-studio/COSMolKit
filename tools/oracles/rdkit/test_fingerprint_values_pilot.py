"""Generator scheduling proofs; not COSMolKit chemistry parity tests."""
import json
import pytest
import fingerprint_values_pilot as oracle


def profiles(name, **axes):
    import itertools
    return [{name: dict(zip(axes, values))}
            for values in itertools.product(*axes.values())] if axes else [name]


MATRICES = {
    "smiles_read": profiles("SmilesRead", sanitize=[False, True], remove_hs=[False, True]),
    "sanitize": profiles("SanitizeAll"),
    "kekulize": profiles("Kekulize", clear_aromatic_flags=[False, True]),
    "molecular_weight": profiles("MolecularWeight", only_heavy=[False, True]),
    "exact_molecular_weight": profiles("ExactMolecularWeight", only_heavy=[False, True]),
    "molecular_formula": profiles("MolecularFormula", separate_isotopes=[False, True], abbreviate_h_isotopes=[False, True]),
    "num_heavy_atoms": profiles("NumHeavyAtoms", remove_hs=[False, True]),
    "total_atom_count": profiles("TotalAtomCount", remove_hs=[False, True]),
    "lipinski_hba": profiles("LipinskiHBA", remove_hs=[False, True]),
    "lipinski_hbd": profiles("LipinskiHBD", remove_hs=[False, True]),
    "fraction_csp3": profiles("FractionCSP3", remove_hs=[False, True]),
    "add_hydrogens": profiles("AddHydrogens", explicit_only=[False, True]),
    "remove_hydrogens": profiles("RemoveHydrogens", sanitize=[False, True]),
    "coordinates_2d": profiles("Coordinates2dDefault"),
    "svg": profiles("SvgDefault"),
    "distance_matrix": profiles("DistanceMatrix", use_bond_order=[False, True], use_atom_weights=[False, True]),
    "fuzzy_and": [{"operation": "FuzzyAnd", "width": width} for width in ("U32", "U64")],
    "fuzzy_or": [{"operation": "FuzzyOr", "width": width} for width in ("U32", "U64")],
}


UFF_GENERATORS = {
    "uff_has_all_molecule_params": "generate_uff_has_all_molecule_params",
    "uff_optimize": "generate_uff_optimize",
    "uff_optimize_conformers": "generate_uff_optimize_conformers",
}

UFF_MATRICES = {
    "uff_has_all_molecule_params": [
        {"Coverage": {"add_hydrogens": False}},
        {"Coverage": {"add_hydrogens": True}},
    ],
    "uff_optimize": [
        {"Optimization": {
            "add_hydrogens": True,
            "max_iterations": 1,
            "vdw_threshold": 100,
            "ignore_interfragment_interactions": True,
            "conformer_id": None,
        }},
    ],
    "uff_optimize_conformers": [
        {"ConformerOptimization": {
            "add_hydrogens": True,
            "max_iterations": 1,
            "vdw_threshold": 100,
            "ignore_interfragment_interactions": True,
            "conformer_count": 2,
        }},
    ],
}

UFF_CASES = [
    {"id": "a", "smiles": "C"},
    {"id": "b", "smiles": "C"},
    {"id": "bad-first", "smiles": "invalid"},
    {"id": "bad-second", "smiles": "invalid"},
]


def inputs(name):
    if name.startswith("fuzzy"):
        return [{"id": str(i), "length": 16, "left": left, "right": right}
                for i, (left, right) in enumerate([
                    ([], []), ([], [[2, 3]]),
                    ([[1, -2147483648], [2, 2147483647], [3, 0]], [[1, 3], [3, 0]]),
                    ([[1, 5], [3, -2]], [[1, 3], [3, -4], [9, 7]])])]
    return [{"id": str(i), "smiles": smiles}
            for i, smiles in enumerate(["CCO", "[2H]O[2H]", "c1ccccc1", "[H]N([H])[H]", "invalid"])]


@pytest.mark.parametrize("name", MATRICES)
def test_generators_preserve_complete_rows_and_order(name):
    generate = oracle.GENERATORS["generate_" + name]
    corpus, parameters = inputs(name), MATRICES[name]
    original = json.dumps([corpus, parameters], sort_keys=True)
    serial = generate(corpus, parameters, 1)
    assert len(serial) == len(corpus) * len(parameters)
    for threads in (2, 4):
        assert generate(corpus, parameters, threads) == serial
    kind = "Fingerprint" if name.startswith("fuzzy") else "Molecular"
    assert [row["input"][kind]["case"]["id"] for row in serial] == [
        case["id"] for case in corpus for _ in parameters]
    assert json.dumps([corpus, parameters], sort_keys=True) == original


@pytest.mark.parametrize("name", ["fuzzy_and", "molecular_weight"])
@pytest.mark.parametrize("threads", [0, -1, 1.5, True, None])
def test_invalid_concurrency_rejected(name, threads):
    with pytest.raises(ValueError, match="positive integer"):
        oracle.GENERATORS["generate_" + name](inputs(name), MATRICES[name], threads)


@pytest.mark.parametrize("name", ["fuzzy_or", "smiles_read"])
@pytest.mark.parametrize("empty_cases,empty_parameters", [(True, False), (False, True), (True, True)])
def test_empty_work_rejected(name, empty_cases, empty_parameters):
    with pytest.raises(ValueError, match="empty"):
        oracle.GENERATORS["generate_" + name](
            [] if empty_cases else inputs(name), [] if empty_parameters else MATRICES[name], 2)


@pytest.mark.parametrize("parameters", [[], ["SvgDefault", "SvgDefault"], ["Coordinates2dDefault"], [{"SvgDefault": {"width": 301}}]])
def test_svg_rejects_every_nonfrozen_parameter_matrix(parameters):
    with pytest.raises(ValueError, match="frozen SvgDefault"):
        oracle.generate_svg(inputs("svg"), parameters, 1)
def _uff_record_bytes(rows):
    return (
        json.dumps([row["input"] for row in rows], sort_keys=True, separators=(",", ":")).encode(),
        json.dumps([row["output"] for row in rows], sort_keys=True, separators=(",", ":")).encode(),
    )


@pytest.mark.parametrize("name", UFF_GENERATORS)
def test_uff_generators_preserve_full_matrix_order_and_parallel_bytes(name):
    generate = oracle.GENERATORS[UFF_GENERATORS[name]]
    corpus = [dict(case) for case in UFF_CASES]
    parameters = [dict(profile) for profile in UFF_MATRICES[name]]
    assert len(parameters) == {
        "uff_has_all_molecule_params": 2,
        "uff_optimize": 1,
        "uff_optimize_conformers": 1,
    }[name]
    original_corpus = json.dumps(corpus, sort_keys=True)
    original_parameters = json.dumps(parameters, sort_keys=True)
    expected_order = [
        (case["id"], profile)
        for case in corpus
        for profile in parameters
    ]

    if name != "uff_has_all_molecule_params":
        stale = generate(
            [{"id": "stale-first", "smiles": "invalid"}],
            [parameters[0]],
            1,
        )
        assert "stale-first" in stale[0]["input"]["Uff"]["preparation"]["Rejected"]["detail"]

    serial = generate(corpus, parameters, 1)
    assert len(serial) == len(corpus) * len(parameters)
    serial_bytes = _uff_record_bytes(serial)
    assert [
        (row["input"]["Uff"]["case"]["id"], row["input"]["Uff"]["profile"])
        for row in serial
    ] == expected_order

    records_by_identity = {
        (row["input"]["Uff"]["case"]["id"], json.dumps(row["input"]["Uff"]["profile"], sort_keys=True)): row
        for row in serial
    }
    for profile in parameters:
        profile_key = json.dumps(profile, sort_keys=True)
        first = records_by_identity[("a", profile_key)]
        duplicate = records_by_identity[("b", profile_key)]
        assert first["input"]["Uff"]["preparation"] == duplicate["input"]["Uff"]["preparation"]
        assert first["output"] == duplicate["output"]
        if "Coverage" not in profile:
            bad_first = records_by_identity[("bad-first", profile_key)]
            bad_second = records_by_identity[("bad-second", profile_key)]
            first_rejection = bad_first["input"]["Uff"]["preparation"]["Rejected"]
            second_rejection = bad_second["input"]["Uff"]["preparation"]["Rejected"]
            assert first_rejection == second_rejection
            assert "bad-first" in first_rejection["detail"]

    assert json.dumps(corpus, sort_keys=True) == original_corpus
    assert json.dumps(parameters, sort_keys=True) == original_parameters
    for threads in (2, 4):
        parallel = generate(corpus, parameters, threads)
        assert len(parallel) == len(expected_order)
        assert [
            (row["input"]["Uff"]["case"]["id"], row["input"]["Uff"]["profile"])
            for row in parallel
        ] == expected_order
        assert _uff_record_bytes(parallel) == serial_bytes
        assert json.dumps(corpus, sort_keys=True) == original_corpus
        assert json.dumps(parameters, sort_keys=True) == original_parameters


@pytest.mark.parametrize("name", UFF_GENERATORS)
@pytest.mark.parametrize("threads", [0, -1, 1.5, True, None])
def test_uff_generators_reject_invalid_concurrency(name, threads):
    with pytest.raises(ValueError, match="positive integer"):
        oracle.GENERATORS[UFF_GENERATORS[name]](
            UFF_CASES, UFF_MATRICES[name], threads)


@pytest.mark.parametrize("name", UFF_GENERATORS)
@pytest.mark.parametrize("empty_cases,empty_parameters", [
    (True, False), (False, True), (True, True),
])
def test_uff_generators_reject_empty_work(name, empty_cases, empty_parameters):
    with pytest.raises(ValueError, match="empty"):
        oracle.GENERATORS[UFF_GENERATORS[name]](
            [] if empty_cases else UFF_CASES,
            [] if empty_parameters else UFF_MATRICES[name],
            2,
        )
