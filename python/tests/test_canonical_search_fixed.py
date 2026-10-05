import pytest
import cosmolkit as ck
from cosmolkit import search


def test_search_module_has_canonical_constructor_and_no_legacy_top_level_alias():
    from cosmolkit.search import parse_smarts
    query = parse_smarts('C(N)O')
    assert isinstance(query, ck.QueryGraph)
    assert query.num_atoms() == 3
    assert query.num_bonds() == 2
    assert len(query) == 3
    assert repr(query) == 'QueryGraph(atoms=3, bonds=2)'
    assert not hasattr(ck, 'parse_smarts')
    assert not isinstance(query, ck.Molecule)


def test_parse_params_all_seven_defaults_are_original_source_defaults():
    p = ck.SmartsParseParams()
    assert (p.allow_cxsmiles,p.strict_cxsmiles,p.parse_name,p.merge_hs,p.skip_cleanup,p.debug_parse,p.replacements) == (True,True,True,False,False,False,{})
    with pytest.raises(AttributeError):
        p.merge_hs = True


def test_parse_params_all_seven_explicit_fields_are_projected():
    p = ck.SmartsParseParams(allow_cxsmiles=False,strict_cxsmiles=False,parse_name=False,merge_hs=True,skip_cleanup=True,debug_parse=True,replacements={'{X}':'N'})
    assert (p.allow_cxsmiles,p.strict_cxsmiles,p.parse_name,p.merge_hs,p.skip_cleanup,p.debug_parse,p.replacements) == (False,False,False,True,True,True,{'{X}':'N'})


@pytest.mark.parametrize(('text','atoms','bonds'), [('C(N)O',3,2),('C1CC1',3,3),('N<-C',2,1),('[$(C),$(N)_7]',1,0)])
def test_original_smarts_branch_ring_dative_and_recursive_inputs(text,atoms,bonds):
    q = search.parse_smarts(text)
    assert (q.num_atoms(),q.num_bonds()) == (atoms,bonds)
    out = search.write_smarts(q,ck.SmartsWriteParams())
    q2 = search.parse_smarts(out)
    assert (q2.num_atoms(),q2.num_bonds()) == (atoms,bonds)


def test_merge_hs_and_replacements_use_explicit_owner_params():
    assert search.parse_smarts('C[H]').num_atoms() == 2
    assert search.parse_smarts_with_params('C[H]',ck.SmartsParseParams(merge_hs=True)).num_atoms() == 1
    q = search.parse_smarts_with_params('{X}{Y}',ck.SmartsParseParams(replacements={'{X}':'C','{Y}':'O'}))
    assert (q.num_atoms(),q.num_bonds()) == (2,1)


def test_name_and_cx_strict_error_preserve_priority_and_typed_category():
    assert search.parse_smarts('CO compound').name() == 'compound'
    with pytest.raises(ck.SmartsParseError) as no_name:
        search.parse_smarts_with_params('CO compound',ck.SmartsParseParams(parse_name=False))
    assert no_name.value.kind == 'CxSmiles'
    assert search.parse_smarts_with_params('CO compound',ck.SmartsParseParams(parse_name=False,strict_cxsmiles=False)).name() is None
    with pytest.raises(ck.SmartsParseError) as caught:
        search.parse_smarts('C |notCX| name')
    assert caught.value.domain == 'search'
    assert caught.value.kind == 'CxSmiles'


def test_parse_diagnostic_gap_is_typed_and_never_empty_query():
    with pytest.raises(ck.SmartsParseError) as caught:
        search.parse_smarts_with_params('C',ck.SmartsParseParams(debug_parse=True))
    assert caught.value.kind == 'UnsupportedFeature'
    assert caught.value.domain == 'search'
    assert caught.value.feature
    with pytest.raises(ck.SmartsParseError) as syntax:
        search.parse_smarts('C(')
    assert syntax.value.kind != 'UnsupportedFeature'


def test_write_params_cover_every_field_and_out_of_range_is_structured():
    p = ck.SmartsWriteParams()
    assert (p.include_atom_maps,p.do_isomeric_smiles,p.include_dative_bonds,p.rooted_at_atom) == (True,True,True,None)
    p = ck.SmartsWriteParams(include_atom_maps=False,do_isomeric_smiles=False,include_dative_bonds=False,rooted_at_atom=9)
    assert (p.include_atom_maps,p.do_isomeric_smiles,p.include_dative_bonds,p.rooted_at_atom) == (False,False,False,9)
    with pytest.raises(ck.SmartsWriteError) as caught:
        search.write_smarts(search.parse_smarts('C'),p)
    assert caught.value.domain == 'search'
    assert caught.value.kind == 'RootedAtomOutOfRange'


def test_match_result_transports_atom_pairs_bonds_and_first_match():
    mol = ck.Molecule.from_smiles('CCO')
    query = search.parse_smarts('CO')
    matches = mol.substruct_matches(query)
    assert len(matches) == 1
    assert matches[0].atom_mapping() == [1,2]
    assert matches[0].bond_mapping() == [1]
    assert matches[0].atom_pairs() == [(0,1),(1,2)]
    assert mol.substruct_match(query).atom_mapping() == [1,2]
    assert mol.has_substruct_match(query) is True
    missing = search.parse_smarts('N')
    assert mol.substruct_match(missing) is None
    assert mol.substruct_matches(missing) == []
    assert mol.has_substruct_match(missing) is False


def test_compiled_query_reuses_order_and_independent_result_values():
    q = search.parse_smarts('CO')
    plan = search.compile_query(q)
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
    q = search.parse_smarts('CC')
    assert len(mol.substruct_matches(q)) == 2
    assert len(mol.substruct_matches_with_params(q,ck.SubstructMatchParams(uniquify=False))) == 4
    limited = mol.substruct_matches_with_params(q,ck.SubstructMatchParams(max_matches=1))
    assert len(limited) == 1
    assert limited[0].atom_mapping() == mol.substruct_match(q).atom_mapping()


def test_original_chirality_and_recursive_matching_branches():
    mol = ck.Molecule.from_smiles('N[C@@H](C)C(=O)O')
    q = search.parse_smarts('N[C@H](C)C(=O)O')
    assert len(mol.substruct_matches(q)) == 1
    assert mol.substruct_matches_with_params(q,ck.SubstructMatchParams(use_chirality=True)) == []
    recursive = search.parse_smarts('[$(C-O)]')
    assert ck.Molecule.from_smiles('CCO').substruct_matches(recursive)[0].atom_mapping() == [1]
    assert ck.Molecule.from_smiles('CCO').substruct_matches_with_params(recursive,ck.SubstructMatchParams(recursion_possible=False)) == []


def test_original_generic_matcher_branch():
    q = search.parse_smarts('C* |$;ALK_p$|')
    options = ck.SubstructMatchParams(use_generic_matchers=True)
    assert len(ck.Molecule.from_smiles('CC').substruct_matches_with_params(q,options)) == 1
    # An atomLabel alone has not undergone SetGenericQueriesFromProperties.
    assert len(ck.Molecule.from_smiles('CO').substruct_matches_with_params(q,options)) == 1
    active = search.parse_smarts('C* |atomProp:1._QueryAtomGenericLabel.ALK|')
    assert len(ck.Molecule.from_smiles('CC').substruct_matches_with_params(active,options)) == 1
    assert ck.Molecule.from_smiles('CO').substruct_matches_with_params(active,options) == []


def test_original_atom_property_matching_branch():
    mol = ck.Molecule.from_smiles('CC |atomProp:0.foo.one:1.foo.two|')
    q = search.parse_smarts('C |atomProp:0.foo.one|')
    assert len(mol.substruct_matches(q)) == 2
    matches = mol.substruct_matches_with_params(q,ck.SubstructMatchParams(atom_properties=['foo']))
    assert [m.atom_mapping() for m in matches] == [[0]]
