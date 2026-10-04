"""Generator scheduling proofs; not COSMolKit chemistry parity tests."""
import json
import pytest
import fingerprint_values_pilot as oracle


def profiles(name, **axes):
    import itertools
    return [{name: dict(zip(axes, values))}
            for values in itertools.product(*axes.values())] if axes else [name]


MATRICES = {
    "smiles_read": profiles("SmilesRead", sanitize=[False, True], remove_hydrogens=[False, True]),
    "sanitize": profiles("SanitizeAll"),
    "kekulize": profiles("Kekulize", clear_aromatic_flags=[False, True]),
    "molecular_weight": profiles("MolecularWeight", only_heavy=[False, True]),
    "exact_molecular_weight": profiles("ExactMolecularWeight", only_heavy=[False, True]),
    "molecular_formula": profiles("MolecularFormula", separate_isotopes=[False, True], abbreviate_h_isotopes=[False, True]),
    "num_heavy_atoms": profiles("NumHeavyAtoms", remove_hydrogens=[False, True]),
    "total_atom_count": profiles("TotalAtomCount", remove_hydrogens=[False, True]),
    "lipinski_hba": profiles("LipinskiHBA", remove_hydrogens=[False, True]),
    "lipinski_hbd": profiles("LipinskiHBD", remove_hydrogens=[False, True]),
    "fraction_csp3": profiles("FractionCSP3", remove_hydrogens=[False, True]),
    "add_hydrogens": profiles("AddHydrogens", explicit_only=[False, True]),
    "remove_hydrogens": profiles("RemoveHydrogens", sanitize=[False, True]),
    "coordinates_2d": profiles("Coordinates2dDefault"),
    "svg": profiles("SvgDefault"),
    "distance_matrix": profiles("DistanceMatrix", use_bond_order=[False, True], use_atom_weights=[False, True]),
    "fuzzy_and": [{"operation": "FuzzyAnd", "width": width} for width in ("U32", "U64")],
    "fuzzy_or": [{"operation": "FuzzyOr", "width": width} for width in ("U32", "U64")],
}


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
