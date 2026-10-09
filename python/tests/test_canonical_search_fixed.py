import pytest
import cosmolkit as ck
import ast
import importlib
import inspect
from pathlib import Path


def test_flat_search_functions_and_query_factory_replace_the_domain_submodule():
    from cosmolkit import parse_smarts
    query = parse_smarts('C(N)O')
    assert isinstance(query, ck.QueryGraph)
    assert query.num_atoms() == 3
    assert query.num_bonds() == 2
    assert len(query) == 3
    assert repr(query) == 'QueryGraph(atoms=3, bonds=2)'
    assert ck.QueryGraph.from_smarts('C(N)O').num_atoms() == 3
    assert not hasattr(ck, 'search')
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module('cosmolkit.search')
    assert not isinstance(query, ck.Molecule)


def test_parse_params_all_seven_defaults_are_original_source_defaults():
    p = ck.SmartsParseParams()
    assert (p.allow_cxsmiles,p.strict_cxsmiles,p.parse_name,p.merge_hs,p.skip_cleanup,p.debug_parse,p.replacements) == (True,True,True,False,False,False,{})
    p.merge_hs = True
    assert p.merge_hs is True


def test_parse_params_all_seven_explicit_fields_are_projected():
    p = ck.SmartsParseParams(allow_cxsmiles=False,strict_cxsmiles=False,parse_name=False,merge_hs=True,skip_cleanup=True,debug_parse=True,replacements={'{X}':'N'})
    assert (p.allow_cxsmiles,p.strict_cxsmiles,p.parse_name,p.merge_hs,p.skip_cleanup,p.debug_parse,p.replacements) == (False,False,False,True,True,True,{'{X}':'N'})


@pytest.mark.parametrize(('text','atoms','bonds'), [('C(N)O',3,2),('C1CC1',3,3),('N<-C',2,1),('[$(C),$(N)_7]',1,0)])
def test_original_smarts_branch_ring_dative_and_recursive_inputs(text,atoms,bonds):
    q = ck.parse_smarts(text)
    assert (q.num_atoms(),q.num_bonds()) == (atoms,bonds)
    out = ck.write_smarts(q,ck.SmartsWriteParams())
    q2 = ck.parse_smarts(out)
    assert (q2.num_atoms(),q2.num_bonds()) == (atoms,bonds)


def test_merge_hs_and_replacements_use_explicit_owner_params():
    assert ck.parse_smarts('C[H]').num_atoms() == 2
    assert ck.parse_smarts_with_params('C[H]',ck.SmartsParseParams(merge_hs=True)).num_atoms() == 1
    q = ck.parse_smarts_with_params('{X}{Y}',ck.SmartsParseParams(replacements={'{X}':'C','{Y}':'O'}))
    assert (q.num_atoms(),q.num_bonds()) == (2,1)


def test_name_and_cx_strict_error_preserve_priority_and_typed_category():
    assert ck.parse_smarts('CO compound').name() == 'compound'
    with pytest.raises(ck.SmartsParseError) as no_name:
        ck.parse_smarts_with_params('CO compound',ck.SmartsParseParams(parse_name=False))
    assert no_name.value.kind == 'CxSmiles'
    assert ck.parse_smarts_with_params('CO compound',ck.SmartsParseParams(parse_name=False,strict_cxsmiles=False)).name() is None
    with pytest.raises(ck.SmartsParseError) as caught:
        ck.parse_smarts('C |notCX| name')
    assert caught.value.domain == 'search'
    assert caught.value.kind == 'CxSmiles'


def test_parse_diagnostic_gap_is_typed_and_never_empty_query():
    with pytest.raises(ck.SmartsParseError) as caught:
        ck.parse_smarts_with_params('C',ck.SmartsParseParams(debug_parse=True))
    assert caught.value.kind == 'UnsupportedFeature'
    assert caught.value.domain == 'search'
    assert caught.value.feature
    with pytest.raises(ck.SmartsParseError) as syntax:
        ck.parse_smarts('C(')
    assert syntax.value.kind != 'UnsupportedFeature'


def test_write_params_cover_every_field_and_out_of_range_is_structured():
    p = ck.SmartsWriteParams()
    assert (p.include_atom_maps,p.isomeric_smiles,p.include_dative_bonds,p.rooted_at_atom) == (True,True,True,None)
    p = ck.SmartsWriteParams(include_atom_maps=False,isomeric_smiles=False,include_dative_bonds=False,rooted_at_atom=9)
    assert (p.include_atom_maps,p.isomeric_smiles,p.include_dative_bonds,p.rooted_at_atom) == (False,False,False,9)
    with pytest.raises(ck.SmartsWriteError) as caught:
        ck.write_smarts(ck.parse_smarts('C'),p)
    assert caught.value.domain == 'search'
    assert caught.value.kind == 'RootedAtomOutOfRange'


def test_match_result_transports_atom_pairs_bonds_and_first_match():
    mol = ck.Molecule.from_smiles('CCO')
    query = ck.parse_smarts('CO')
    matches = mol.substruct_matches(query)
    assert len(matches) == 1
    assert matches[0].atom_mapping() == [1,2]
    assert matches[0].bond_mapping() == [1]
    assert matches[0].atom_pairs() == [(0,1),(1,2)]
    assert mol.substruct_match(query).atom_mapping() == [1,2]
    assert mol.has_substruct_match(query) is True
    missing = ck.parse_smarts('N')
    assert mol.substruct_match(missing) is None
    assert mol.substruct_matches(missing) == []
    assert mol.has_substruct_match(missing) is False


def test_compiled_query_reuses_order_and_independent_result_values():
    q = ck.parse_smarts('CO')
    plan = ck.compile_query(q)
    assert (plan.num_atoms(),plan.num_bonds()) == (2,1)
    before = plan.atom_order()
    assert sorted(before) == [0,1]
    assert ck.Molecule.from_smiles('CCO').substruct_matches_compiled(plan)[0].atom_mapping() == [1,2]
    assert ck.Molecule.from_smiles('NN').substruct_matches_compiled(plan) == []
    assert plan.atom_order() == before
    assert plan.query().num_atoms() == q.num_atoms()


def test_substruct_params_original_scalar_defaults_and_overrides():
    p = ck.SubstructMatchParams()
    assert (p.max_matches,p.uniquify,p.use_chirality,p.use_enhanced_stereo,p.specified_stereo_query_matches_unspecified,p.use_query_query_matches,p.recursion_possible,p.max_recursive_matches,p.num_threads,p.aromatic_matches_conjugated,p.aromatic_matches_single_or_double,p.atom_properties,p.bond_properties,p.extra_atom_check_overrides_default_check,p.extra_bond_check_overrides_default_check,p.use_generic_matchers) == (1000,True,False,False,False,False,True,1000,1,False,False,[],[],False,False,False)
    p = ck.SubstructMatchParams(max_matches=3,uniquify=False,use_chirality=True,use_enhanced_stereo=True,specified_stereo_query_matches_unspecified=True,use_query_query_matches=True,recursion_possible=False,max_recursive_matches=7,num_threads=2,aromatic_matches_conjugated=True,aromatic_matches_single_or_double=True,atom_properties=['foo'],bond_properties=['bar'],extra_atom_check_overrides_default_check=True,extra_bond_check_overrides_default_check=True,use_generic_matchers=True)
    assert (p.max_matches,p.uniquify,p.use_chirality,p.use_enhanced_stereo,p.specified_stereo_query_matches_unspecified,p.use_query_query_matches,p.recursion_possible,p.max_recursive_matches,p.num_threads,p.aromatic_matches_conjugated,p.aromatic_matches_single_or_double,p.atom_properties,p.bond_properties,p.extra_atom_check_overrides_default_check,p.extra_bond_check_overrides_default_check,p.use_generic_matchers) == (3,False,True,True,True,True,False,7,2,True,True,['foo'],['bar'],True,True,True)


def test_original_match_limits_uniquify_and_orientation_branches():
    mol = ck.Molecule.from_smiles('CCCO')
    q = ck.parse_smarts('CC')
    assert len(mol.substruct_matches(q)) == 2
    assert len(mol.substruct_matches_with_params(q,ck.SubstructMatchParams(uniquify=False))) == 4
    limited = mol.substruct_matches_with_params(q,ck.SubstructMatchParams(max_matches=1))
    assert len(limited) == 1
    assert limited[0].atom_mapping() == mol.substruct_match(q).atom_mapping()


def test_original_chirality_and_recursive_matching_branches():
    mol = ck.Molecule.from_smiles('N[C@@H](C)C(=O)O')
    q = ck.parse_smarts('N[C@H](C)C(=O)O')
    assert len(mol.substruct_matches(q)) == 1
    assert mol.substruct_matches_with_params(q,ck.SubstructMatchParams(use_chirality=True)) == []
    recursive = ck.parse_smarts('[$(C-O)]')
    assert ck.Molecule.from_smiles('CCO').substruct_matches(recursive)[0].atom_mapping() == [1]
    assert ck.Molecule.from_smiles('CCO').substruct_matches_with_params(recursive,ck.SubstructMatchParams(recursion_possible=False)) == []


def test_original_generic_matcher_branch():
    q = ck.parse_smarts('C* |$;ALK_p$|')
    options = ck.SubstructMatchParams(use_generic_matchers=True)
    assert len(ck.Molecule.from_smiles('CC').substruct_matches_with_params(q,options)) == 1
    # An atomLabel alone has not undergone SetGenericQueriesFromProperties.
    assert len(ck.Molecule.from_smiles('CO').substruct_matches_with_params(q,options)) == 1
    active = ck.parse_smarts('C* |atomProp:1._QueryAtomGenericLabel.ALK|')
    assert len(ck.Molecule.from_smiles('CC').substruct_matches_with_params(active,options)) == 1
    assert ck.Molecule.from_smiles('CO').substruct_matches_with_params(active,options) == []


def test_original_atom_property_matching_branch():
    mol = ck.Molecule.from_smiles('CC |atomProp:0.foo.one:1.foo.two|')
    q = ck.parse_smarts('C |atomProp:0.foo.one|')
    assert len(mol.substruct_matches(q)) == 2
    matches = mol.substruct_matches_with_params(q,ck.SubstructMatchParams(atom_properties=['foo']))
    assert [m.atom_mapping() for m in matches] == [[0]]


@pytest.mark.parametrize('text', ['C(N)O', 'C1CC1', 'N<-C', '[$(C),$(N)_7]', 'CO compound'])
def test_both_smarts_factories_preserve_the_same_query_and_match_results(text):
    text_before = text
    queries = [ck.parse_smarts(text), ck.QueryGraph.from_smarts(text)]
    rows = [(q.num_atoms(), q.num_bonds(), q.name(), ck.write_smarts(q, ck.SmartsWriteParams())) for q in queries]
    assert rows[0] == rows[1]
    target = ck.Molecule.from_smiles('CCO')
    before = (target.to_smiles(), target.num_atoms(), target.num_bonds(), target.coordinates_2d())
    matches = [[m.atom_mapping() for m in target.substruct_matches(q)] for q in queries]
    assert matches[0] == matches[1]
    assert (target.to_smiles(), target.num_atoms(), target.num_bonds(), target.coordinates_2d()) == before
    assert text == text_before


def test_both_explicit_factories_forward_params_and_preserve_them():
    params = ck.SmartsParseParams(merge_hs=True, replacements={'{X}': 'C'})
    before = (params.merge_hs, params.replacements)
    for factory in (ck.parse_smarts_with_params, ck.QueryGraph.from_smarts_with_params):
        query = factory('{X}[H]', params)
        assert (query.num_atoms(), query.num_bonds()) == (1, 0)
        assert (params.merge_hs, params.replacements) == before


@pytest.mark.parametrize('text', ['C(', '[#6', 'C |notCX| name'])
def test_both_factories_preserve_typed_parse_error_details(text):
    errors = []
    for factory in (ck.parse_smarts, ck.QueryGraph.from_smarts):
        with pytest.raises(ck.SmartsParseError) as caught:
            factory(text)
        error = caught.value
        errors.append((type(error), str(error), error.domain, error.kind, error.__dict__, str(error.__cause__)))
    assert errors[0] == errors[1]


def test_both_explicit_factories_preserve_unsupported_error_details():
    params = ck.SmartsParseParams(debug_parse=True)
    errors = []
    for factory in (ck.parse_smarts_with_params, ck.QueryGraph.from_smarts_with_params):
        with pytest.raises(ck.SmartsParseError) as caught:
            factory('C', params)
        error = caught.value
        errors.append((type(error), str(error), error.domain, error.kind, error.feature))
        assert params.debug_parse is True
    assert errors[0] == errors[1]


def test_flat_functions_and_class_factories_have_real_generated_signatures():
    stub = ast.parse((Path(__file__).resolve().parents[1] / 'cosmolkit.pyi').read_text())
    functions = {n.name: n for n in stub.body if isinstance(n, ast.FunctionDef)}
    query = next(n for n in stub.body if isinstance(n, ast.ClassDef) and n.name == 'QueryGraph')
    methods = {}
    for node in query.body:
        if isinstance(node, ast.FunctionDef):
            methods.setdefault(node.name, node)
    for name in ('parse_smarts', 'parse_smarts_with_params', 'compile_query', 'write_smarts', 'write_cx_smarts'):
        assert name in functions
        assert callable(getattr(ck, name))
    for name, parameters in [('from_smarts', ['text']), ('from_smarts_with_params', ['text', 'params'])]:
        assert list(inspect.signature(getattr(ck.QueryGraph, name)).parameters) == parameters
        method = methods[name]
        assert [a.arg for a in method.args.args] == parameters
        assert [ast.unparse(d) for d in method.decorator_list] == (
            ['typing.overload', 'staticmethod'] if name == 'from_smarts' else ['staticmethod'])
        if name == 'from_smarts':
            assert len([node for node in query.body if isinstance(node, ast.FunctionDef) and node.name == name]) == 3
        assert ast.unparse(method.returns) == 'QueryGraph'
