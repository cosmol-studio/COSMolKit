"""Focused MCS language-boundary regressions; no external reference runtime."""
import pytest
import cosmolkit as ck


def molecules(*smiles):
    return [ck.Molecule.from_smiles(s) for s in smiles]


@pytest.mark.parametrize("smiles,atoms,bonds,smarts", [
    (("CCO", "CCN"), 2, 1, "[#6]-[#6]"),
    (("Cl", "Br"), 0, 0, ""),
    (("c1ccccc1O", "c1ccccc1N"), 6, 6, "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"),
])
def test_mcs_default_query_and_readonly_inputs(smiles, atoms, bonds, smarts):
    # Frozen with pinned RDKit 2026.03.6 FindMCS defaults.
    inputs = molecules(*smiles)
    before = [m.to_binary() for m in inputs]
    result = ck.maximum_common_substructure(inputs)
    assert isinstance(result, ck.McsResult)
    assert (result.atom_count, result.bond_count, result.completed, result.smarts) == (
        atoms, bonds, True, smarts,
    )
    assert [m.to_binary() for m in inputs] == before
    if atoms:
        assert isinstance(result.query, ck.QueryGraph)
        assert all(m.has_substruct_match(result.query) for m in inputs)
    else:
        assert result.query is None
    with pytest.raises(AttributeError):
        result.atom_count = 99


def test_mcs_configuration_forms_enum_strings_and_forwarding():
    inputs = molecules("CCO", "CCN")
    for comparator in (ck.McsAtomComparator.Any, "any"):
        params = ck.McsParameters(atom_comparator=comparator)
        assigned = ck.McsParameters()
        assigned.atom_comparator = comparator
        results = [
            ck.maximum_common_substructure(inputs, params),
            ck.maximum_common_substructure(inputs, assigned),
            ck.maximum_common_substructure(inputs, atom_comparator=comparator),
            ck.maximum_common_substructure_with_params(inputs, params),
        ]
        assert all((r.atom_count, r.bond_count, r.completed) == (3, 2, True) for r in results)
        assert all(r.smarts == "[#6]-[#6]-[#8,#7]" for r in results)
        assert all(m.has_substruct_match(results[0].query) for m in inputs)
    with pytest.raises(TypeError):
        ck.maximum_common_substructure(inputs, params, timeout=1)
    with pytest.raises(TypeError):
        ck.maximum_common_substructure(inputs, unknown_option=True)
    with pytest.raises(AttributeError):
        params.unknown_option = True
    previous = params.atom_comparator
    with pytest.raises(ValueError):
        params.atom_comparator = "not_a_comparator"
    assert params.atom_comparator == previous


@pytest.mark.parametrize("smiles,smarts", [
    (("CCO", "CCN"), "[#6]-[#6]-[#8,#7]"),
    (("c1ccccc1O", "c1ccccc1N"), "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#8,#7]"),
    (("CCOC", "CCNC"), "[#6]-[#6]-[#8,#7]-[#6]"),
    (("CCO", "CCN", "CCS"), "[#6]-[#6]-[#8,#7,#16]"),
    (("CCO", "C=CO"), "[#6]-,=[#6]-[#8]"),
])
def test_mcs_result_alternatives_preserve_exact_source_order(smiles, smarts):
    # RDKit 2026.03.6 MaximumCommonSubgraph.cpp uses expandQuery(..., OR)
    # with QueryAtom/QueryBond's default maintainOrder=true. Preserve the
    # exact text and query tree, rather than sorting chemically equivalent ORs.
    inputs = molecules(*smiles)
    result = ck.maximum_common_substructure(
        inputs, atom_comparator="any", bond_comparator="any",
    )
    assert result.completed is True
    assert result.smarts == smarts
    assert ck.write_smarts(result.query, ck.SmartsWriteParams()) == smarts
    assert all(m.has_substruct_match(result.query) for m in inputs)


def test_mcs_nested_configuration_mutation_is_effective_and_independent():
    inputs = molecules("[NH4+]", "N")
    assert ck.maximum_common_substructure(inputs).atom_count == 1
    original = ck.McsAtomCompareParameters()
    params = ck.McsParameters(atom_compare_parameters=original)
    params.atom_compare_parameters.match_formal_charge = True
    assert params.atom_compare_parameters.match_formal_charge is True
    assert original.match_formal_charge is False
    assert ck.maximum_common_substructure(inputs, params).atom_count == 0
    assert ck.maximum_common_substructure(inputs, match_formal_charge=True).atom_count == 0
    params.bond_compare_parameters.match_stereo = True
    params.threshold = 0.5
    assert params.atom_compare_parameters.match_formal_charge is True
    assert params.bond_compare_parameters.match_stereo is True


def test_mcs_result_lifetime_degenerate_results_and_typed_errors():
    inputs = molecules("CCO", "CCN")
    result = ck.maximum_common_substructure(inputs, store_all=True)
    assert (result.atom_count, result.bond_count, result.completed) == (2, 1, True)
    assert result.degenerate
    assert all(isinstance(q, ck.QueryGraph) for q in result.degenerate.values())
    result = ck.maximum_common_substructure(inputs)
    del inputs
    assert result.query.num_atoms() == 2
    for inputs, options in [([], {}), (molecules("CC"), {}), (molecules("CC", "CO"), {"threshold": 1.1})]:
        with pytest.raises(ck.McsError) as failure:
            ck.maximum_common_substructure(inputs, **options)
        assert failure.value.domain == "mcs"
        assert failure.value.kind == "State"
    with pytest.raises(OverflowError):
        ck.McsParameters(timeout=-1)
