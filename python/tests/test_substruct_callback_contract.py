"""Small native callback regressions; no external corpus or reference files."""
import gc
import weakref

import pytest

import cosmolkit as ck


@pytest.mark.parametrize("method", ["substruct_matches", "substruct_match", "has_substruct_match"])
@pytest.mark.parametrize("field,alias", [("final_match", "extra_final_check"), ("atom_match", "extra_atom_check"), ("bond_match", "extra_bond_check")])
def test_callbacks_work_in_constructor_properties_keywords_and_explicit_methods(method, field, alias):
    molecule = ck.Molecule.from_smiles("CCC")
    query = ck.parse_smarts("CC")
    baseline = getattr(molecule, method)(query)
    assert baseline
    seen = []

    def reject(*args):
        seen.append(args)
        return False

    params = ck.SubstructMatchParams(**{field: reject})
    for call in (
        lambda: getattr(molecule, method)(query, params),
        lambda: getattr(molecule, method + "_with_params")(query, params),
        lambda: getattr(molecule, method)(query, **{field: reject}),
        lambda: getattr(molecule, method)(query, **{alias: reject}),
    ):
        seen.clear()
        assert not call()
        assert seen
    setattr(params, alias, None)
    assert getattr(params, field) is None
    assert getattr(molecule, method)(query, params)
    setattr(params, alias, reject)
    params.max_matches = 7
    assert getattr(params, field) is reject
    assert not getattr(molecule, method)(query, params)
    with pytest.raises(TypeError):
        getattr(molecule, method)(query, **{field: reject, alias: reject})


@pytest.mark.parametrize("field", ["final_match", "atom_match", "bond_match"])
@pytest.mark.parametrize("method", ["substruct_matches", "substruct_match", "has_substruct_match"])
def test_callback_exceptions_preserve_identity_stop_dispatch_and_do_not_poison_params(field, method):
    molecule = ck.Molecule.from_smiles("CCC")
    query = ck.parse_smarts("CC")
    original = LookupError("original callback failure")
    calls = []

    def fail(*args):
        calls.append(args)
        raise original

    params = ck.SubstructMatchParams(**{field: fail})
    with pytest.raises(LookupError) as captured:
        getattr(molecule, method)(query, params)
    assert captured.value is original
    assert len(calls) == 1
    setattr(params, field, lambda *_: True)
    assert getattr(molecule, method)(query, params)
    assert molecule.to_smiles() == "CCC"


def test_atom_and_bond_callbacks_receive_readonly_source_values_including_wildcards():
    molecule = ck.Molecule.from_smiles("CO")
    query = ck.parse_smarts("[#6,#8]")
    seen = []

    retained = []

    def oxygen_only(query_graph, query_atom, target, target_atom):
        assert isinstance(query_graph, ck.QueryGraph)
        assert query_graph.num_atoms == 1 and query_graph.num_bonds == 0
        assert isinstance(query_atom, ck.QueryAtom)
        assert isinstance(target_atom, ck.Atom)
        contextual = target.atom(target_atom.id)
        assert contextual.explicit_valence == 1
        assert contextual.total_hydrogens == (1 if target_atom.atomic_number == 8 else 3)
        assert target.num_atoms == 2
        target.add_hydrogens_()
        retained.append((query_graph, target, contextual))
        seen.append(query_atom)
        return target_atom.atomic_number == 8

    assert [list(match.atom_mapping) for match in molecule.substruct_matches(query, atom_match=oxygen_only)] == [[1]]
    assert seen[0].atomic_number == 0  # Query wildcard must not become carbon.
    assert seen[0].id == 0
    assert molecule.num_atoms == 2
    del molecule, query
    assert all(graph.num_atoms == 1 and target.num_atoms == 6 for graph, target, _ in retained)
    assert all(atom.degree == 1 for _, _, atom in retained)
    assert {tuple(match.atom_mapping) for match in ck.Molecule.from_smiles("CC=C").substruct_matches(
        ck.parse_smarts("C~C"), uniquify=False,
        bond_match=lambda query_bond, target_bond: target_bond.order_name == "DOUBLE",
    )} == {(1, 2), (2, 1)}


def test_final_callback_has_owned_cow_target_and_query_order_mapping():
    molecule = ck.Molecule.from_smiles("CCC")
    snapshots = []

    def accept_last(target, ids):
        assert isinstance(target, ck.Molecule)
        target.add_hydrogens_()
        snapshots.append(target)
        return tuple(ids) == (1, 2)

    assert [list(match.atom_mapping) for match in molecule.substruct_matches(
        ck.parse_smarts("CC"), uniquify=False, final_match=accept_last,
    )] == [[1, 2]]
    assert molecule.num_atoms == 3
    assert molecule.to_smiles() == "CCC"
    assert snapshots and all(target.num_atoms == 11 for target in snapshots)


def test_atom_callback_retains_recursive_query_and_inplace_reaction_input():
    seen = []

    def accept(query, query_atom, target, target_atom):
        seen.append((query, target))
        assert query_atom.id < query.num_atoms
        assert target_atom.id < target.num_atoms
        target.add_hydrogens_()
        return True

    molecule = ck.mol_from_smiles("CO")
    assert molecule.has_substruct_match(ck.parse_smarts("[$(CO)]"), atom_match=accept)
    assert {query.num_atoms for query, _ in seen} == {1, 2}
    assert molecule.num_atoms == 2
    assert all(target.num_atoms == 6 for _, target in seen)

    seen.clear()
    reaction = ck.Reaction.from_smirks("[C:1]>>[N:1]")
    reaction.set_match_params(ck.SubstructMatchParams(atom_match=accept))
    molecule = ck.mol_from_smiles("C")
    assert molecule.apply_reaction_(reaction)
    assert molecule.num_atoms == 1 and molecule.atom(0).atomic_number == 7
    del molecule, reaction
    assert seen and all(target.num_atoms == 5 for _, target in seen)
    assert all(target.atom(0).atomic_number == 6 for _, target in seen)


def test_callback_assignment_validation_and_cycle_collection():
    params = ck.SubstructMatchParams()
    original = lambda *_: True
    params.extra_final_check = original
    with pytest.raises(TypeError):
        params.extra_final_check = 42
    assert params.extra_final_check is original
    with pytest.raises(AttributeError):
        params.extra_final_chek = original
    callback = lambda *_args, keep=params: keep is not None
    owner = weakref.ref(callback)
    params.final_match = callback
    del params, callback
    gc.collect()
    assert owner() is None


def test_reentrant_matching_has_independent_error_slots():
    molecule = ck.Molecule.from_smiles("CCC")
    query = ck.parse_smarts("C")
    nested = ck.SubstructMatchParams(final_match=lambda *_: False)
    params = ck.SubstructMatchParams(final_match=lambda *_: not molecule.has_substruct_match(query, nested))
    assert len(molecule.substruct_matches(query, params)) == 3
