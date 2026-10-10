"""Fixed public reaction boundary regressions; no corpus or external oracle."""
import pytest
import cosmolkit as ck


def test_reaction_run_without_molecule_receiver():
    carbon = ck.Molecule.from_smiles("C")
    oxygen = ck.Molecule.from_smiles("O")
    before = [carbon.to_binary(), oxygen.to_binary()]
    rxn = ck.Reaction.from_smirks("[C:1].[O:2]>>[C:1][O:2]")
    groups = rxn.run([carbon, oxygen], ck.ReactionRunParams())
    assert [[mol.to_smiles() for mol in group] for group in groups] == [["CO"]]
    assert rxn.is_initialized()
    assert [carbon.to_binary(), oxygen.to_binary()] == before
    with pytest.raises(ValueError):
        rxn.run([carbon], ck.ReactionRunParams())
    assert [carbon.to_binary(), oxygen.to_binary()] == before
    duplicate = ck.Reaction.from_smirks("[C:1].[C:2]>>[C:1].[C:2]")
    assert [[mol.to_smiles() for mol in group] for group in
            duplicate.run([carbon, carbon], ck.ReactionRunParams())] == [["C", "C"]]


@pytest.mark.parametrize("form", ["default", "none", "positional", "named", "keywords"])
def test_reaction_run_configuration_forms_preserve_inputs(form):
    carbon = ck.Molecule.from_smiles("C").with_atom_property(0, "tracking_id", 42)
    oxygen = ck.Molecule.from_smiles("O")
    before = [carbon.to_binary(), oxygen.to_binary()]
    reaction = ck.Reaction.from_smirks("[C:1].[O:2]>>[C:1][O:2]")
    params = ck.ReactionRunParams(copy_atom_properties=True)
    args, kwargs = {
        "default": ((), {}),
        "none": ((), {"params": None}),
        "positional": ((params,), {}),
        "named": ((), {"params": params}),
        "keywords": ((), {"copy_atom_properties": True, "max_products": 1000,
                            "coordinate_selections": [ck.ReactionCoordinateSelection.auto()] * 2}),
    }[form]
    groups = reaction.run([carbon, oxygen], *args, **kwargs)
    assert [[mol.to_smiles() for mol in group] for group in groups] == [["CO"]]
    assert groups[0][0].atom_property(0, "tracking_id") == (
        42 if form in ("positional", "named", "keywords") else None
    )
    assert [carbon.to_binary(), oxygen.to_binary()] == before
    assert params.copy_atom_properties is True and params.max_products == 1000


def test_reaction_run_keywords_limit_products_and_reject_conflicts_before_execution():
    source = ck.Molecule.from_smiles("CC")
    before = source.to_binary()
    params = ck.ReactionRunParams(max_products=1)
    for args, kwargs in [((params,), {}), ((), {"max_products": 1}),
                         ((), {"params": None, "max_products": 1})]:
        reaction = ck.Reaction.from_smirks("[C:1]>>[N:1]")
        groups = reaction.run([source], *args, **kwargs)
        assert [[mol.to_smiles() for mol in group] for group in groups] == [["CN"]]
        assert reaction.is_initialized()
    reaction = ck.Reaction.from_smirks("[C:1]>>[N:1]")
    assert len(reaction.run([source])) == 2
    for args, kwargs, error in [
        ((params,), {"max_products": 1}, TypeError),
        ((), {"params": params, "copy_atom_properties": False}, TypeError),
        ((), {"unknown_option": True}, TypeError),
        ((), {"params": {}}, TypeError),
        ((), {"max_products": -1}, OverflowError),
    ]:
        fresh = ck.Reaction.from_smirks("[C:1]>>[N:1]")
        with pytest.raises(error):
            fresh.run([source], *args, **kwargs)
        assert not fresh.is_initialized()
    assert source.to_binary() == before
    assert params.max_products == 1


def test_sanitized_products_do_not_gain_nitrogen_stereo():
    # Fixed RDKit 2026.03.1 outputs; no external oracle or corpus dependency.
    source = ck.Molecule.from_smiles("C1C[C@H]2CC[C@H]2C1")
    original = source.to_smiles()
    reaction = ck.Reaction.from_smirks("[C;H1:1]>>[N:1]")
    groups = source.reaction_products_from_inputs(reaction, [source], ck.ReactionRunParams())
    assert len(groups) == 2
    for group, expected in zip(groups, ["C1C[C@@H]2CCN2C1", "C1C[C@H]2CCN2C1"], strict=True):
        assert len(group) == 1
        raw = group[0]
        assert raw.to_smiles() == expected
        sanitized = raw.sanitize()
        assert sanitized.to_smiles() == expected
        assert sanitized.without_hydrogens().to_smiles() == expected
        assert raw.to_smiles() == expected
    assert source.to_smiles() == original


def test_template_writer_clears_unpaired_directions():
    reaction = ck.Reaction.from_smirks("[C:1]>>[C:1]/C=C")
    params = ck.ReactionWriteParams(isomeric_smiles=False)
    for text in [reaction.to_smirks(), reaction.to_smirks_with_params(params),
                 reaction.to_cx_smirks(), reaction.to_cx_smirks_with_params(params)]:
        assert text == "[C:1]>>[C:1]C=C"


def test_reaction_values_writers_and_initialization():
    rxn = ck.Reaction.from_smirks("[C:1]>O>[N:1]")
    assert (rxn.num_reactant_templates(), rxn.num_product_templates(), rxn.num_agent_templates()) == (1, 1, 1)
    assert rxn.reactant_template(0).num_atoms() == 1
    assert len(rxn.reactant_templates()) == len(rxn.product_templates()) == len(rxn.agent_templates()) == 1
    assert not rxn.is_initialized()
    report = rxn.validate_with_params(ck.ReactionValidationParams(silent=True))
    assert report.is_valid and report.num_errors == 0 and report.errors == []
    assert rxn.with_initialized().is_initialized()
    assert not rxn.is_initialized()
    assert rxn.with_initialized_with_params(ck.ReactionValidationParams(silent=True)).is_initialized()
    for text in [rxn.to_smirks(), rxn.to_smirks_with_params(ck.ReactionWriteParams()),
                 rxn.to_cx_smirks(), rxn.to_cx_smirks_with_params(ck.ReactionWriteParams())]:
        assert isinstance(text, str)
        assert ck.parse_smirks(text).num_agent_templates() == 1


def test_reaction_templates_are_owned_values():
    query = ck.QueryGraph.from_smarts("[C:1]")
    rxn = ck.Reaction.from_templates([query], [query], [])
    extended = rxn.with_agent_template(query).with_reactant_template(query).with_product_template(query)
    assert (rxn.num_reactant_templates(), rxn.num_agent_templates(), rxn.num_product_templates()) == (1, 0, 1)
    assert (extended.num_reactant_templates(), extended.num_agent_templates(), extended.num_product_templates()) == (2, 1, 2)
    removed = extended.without_agents()
    assert len(removed.removed_templates) == 1
    assert removed.reaction.num_agent_templates() == 0
    assert extended.num_agent_templates() == 1
    for name in ["without_unmapped_reactants", "without_unmapped_products"]:
        assert isinstance(getattr(rxn, name)(), ck.ReactionTemplateRemoval)
        assert isinstance(getattr(rxn, name + "_with_params")(ck.ReactionTemplateRemovalParams()), ck.ReactionTemplateRemoval)
    # Reaction.h:376 defaults manual ChemicalReaction templates to false.
    changed = rxn.with_implicit_properties(True)
    assert changed.implicit_properties() and not rxn.implicit_properties()
    assert rxn.with_match_params(ck.SubstructMatchParams(max_matches=7)).match_params().max_matches == 7


def test_reaction_parameter_defaults_and_coordinate_payloads():
    params = ck.ReactionParseParams()
    assert (params.use_smiles, params.sanitize, params.replacements, params.allow_cxsmiles, params.strict_cxsmiles) == (False, False, {}, True, True)
    assert ck.ReactionRunParams().max_products == 1000
    assert ck.ReactionRunParams().coordinate_selections == []
    assert ck.ReactionApplyParams().remove_unmatched_atoms
    assert ck.ReactionTemplateRemovalParams().threshold_unmapped_atoms == 0.2
    auto = ck.ReactionCoordinateSelection.auto()
    two = ck.ReactionCoordinateSelection.two_d(7)
    three = ck.ReactionCoordinateSelection.three_d(42)
    assert auto.is_auto and auto.id is None
    assert two.is_2d and two.id == 7 and not two.is_3d
    assert three.is_3d and three.id == 42
    assert ck.ReactionSingleRunParams(coordinate_selection=three).coordinate_selection.id == 42
    selections = ck.ReactionRunParams(max_products=3, coordinate_selections=[two, three]).coordinate_selections
    assert [s.id for s in selections] == [7, 42]
    write = ck.ReactionWriteParams(cx_fields=ck.CxSmilesFields.NONE, rooted_at_atom=0, coordinate_selections=[two])
    assert write.cx_fields.bits() == 0 and write.rooted_at_atom == 0
    assert write.coordinate_selections[0].id == 7
    params.sanitize = True
    assert params.sanitize is True
    with pytest.raises(TypeError):
        _ = ck.ReactionSingleRunParams(coordinate_selection="Auto")  # pyright: ignore[reportArgumentType]
    assert ck.parse_smirks_with_params("{C}>>[N:1]", ck.ReactionParseParams(replacements={"{C}": "[C:1]"})).num_reactant_templates() == 1


@pytest.mark.parametrize("explicit", [False, True])
def test_product_groups_mutate_reaction_initialization_not_input(explicit: bool):
    source = ck.Molecule.from_smiles("C")
    original = source.to_smiles()
    rxn = ck.Reaction.from_smirks("[C:1]>>[N:1]")
    assert not rxn.is_initialized()
    sets = (source.reaction_products_with_params(rxn, 0, ck.ReactionSingleRunParams())
            if explicit else source.reaction_products(rxn, 0))
    assert rxn.is_initialized()
    assert len(sets) == len(sets[0]) == 1
    assert sets[0][0].to_smiles() == "N"
    assert source.to_smiles() == original
    assert ck.Molecule.from_smiles("O").reaction_products(rxn, 0) == []
    pair = ck.Reaction.from_smirks("[C:1].[C:2]>>[C:1].[C:2]")
    sets = source.reaction_products_from_inputs(pair, [source, source], ck.ReactionRunParams())
    assert len(sets) == 1 and len(sets[0]) == 2
    assert [m.to_smiles() for m in sets[0]] == ["C", "C"]
    assert source.to_smiles() == original


@pytest.mark.parametrize("explicit", [False, True])
def test_apply_value_inplace_and_failure_atomicity(explicit: bool):
    source = ck.Molecule.from_smiles("C")
    rxn = ck.Reaction.from_smirks("[C:1]>>[N:1]")
    result = (source.apply_reaction_with_params(rxn, ck.ReactionApplyParams())
              if explicit else source.apply_reaction(rxn))
    assert result.changed and result.molecule.to_smiles() == "N"
    assert source.to_smiles() == "C"
    changed = (source.apply_reaction_with_params_(rxn, ck.ReactionApplyParams())
               if explicit else source.apply_reaction_(rxn))
    assert changed is True and source.to_smiles() == "N"
    assert source.apply_reaction_(rxn) is False
    with pytest.raises(ck.OperationError) as caught:
        _ = source.apply_reaction_(ck.Reaction.from_smirks("[N:1]>>[N:1]C"))
    assert caught.value.kind == "ReactionApply"
    assert isinstance(caught.value.__cause__, ck.ReactionApplyError)
    assert caught.value.__cause__.kind == "AddsProductAtom"
    assert caught.value.__cause__.atom == 1
    assert source.to_smiles() == "N"


def test_typed_parse_validation_model_and_execution_errors():
    with pytest.raises(ck.ReactionParseError) as caught:
        _ = ck.Reaction.from_smirks("CC")
    assert (caught.value.kind, caught.value.count) == ("Separators", 0)
    rxn = ck.Reaction.from_smirks("[C:1]>>[N:1]")
    with pytest.raises(ck.ReactionModelError) as caught:
        _ = rxn.reactant_template(9)
    assert caught.value.kind == "TemplateIndex" and caught.value.index == 9 and caught.value.count == 1
    assert caught.value.role == ck.ReactionRole.Reactant
    with pytest.raises(ck.OperationError) as caught:
        _ = ck.Molecule.from_smiles("C").reaction_products(rxn, 99)
    assert caught.value.kind == "ReactionRun"
    assert isinstance(caught.value.__cause__, ck.ReactionRunError)
    assert caught.value.__cause__.kind == "ReactantTemplateIndex" and caught.value.__cause__.index == 99
    empty = ck.Reaction()
    report = empty.validate_with_params(ck.ReactionValidationParams(silent=True))
    assert not report.is_valid and report.num_errors == 2
    assert report.errors[0].kind == ck.ReactionValidationIssueKind.MissingReactants
    assert report.errors[0].severity == ck.ReactionValidationSeverity.Error
    assert report.errors[0].role == ck.ReactionRole.Reactant
    assert report.errors[0].atom is None
    assert isinstance(report.errors[0].detail, str)
    with pytest.raises(ck.ReactionInitializationError) as caught:
        _ = empty.with_initialized_with_params(ck.ReactionValidationParams(silent=True))
    assert caught.value.kind == "Invalid" and caught.value.report.num_errors == 2


def test_reaction_template_parse_errors_keep_concrete_search_cause():
    with pytest.raises(ck.ReactionParseError) as caught:
        _ = ck.Reaction.from_smirks("[C:1]>>[")
    assert caught.value.kind == "Smarts"
    assert caught.value.role == ck.ReactionRole.Product
    assert caught.value.template == 0 and caught.value.text == "["
    assert isinstance(caught.value.__cause__, ck.SmartsParseError)
    assert caught.value.__cause__.kind == "UnclosedBracket"
