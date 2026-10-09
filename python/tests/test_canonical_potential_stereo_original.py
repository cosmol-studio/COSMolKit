import enum
import pytest
import cosmolkit

def test_potential_stereo_analysis_returns_typed_isolated_state():
    molecule = cosmolkit.Molecule.from_smiles("FC(Cl)Br")
    before = molecule.to_smiles()

    analysis = molecule.potential_stereo_with_params(cosmolkit.PotentialStereoParams(clean=True))
    assert isinstance(analysis, cosmolkit.PotentialStereoResult)
    assert len(analysis) == 1
    assert analysis.cleaned_molecule is not molecule
    assert analysis.cleaned_molecule.to_smiles() == before
    assert molecule.to_smiles() == before

    record = analysis.stereo[0]
    assert isinstance(record, cosmolkit.PotentialStereoInfo)
    assert record.stereo_type == "atom_tetrahedral"
    assert record.specified == "unspecified"
    assert record.centered_on.kind == "atom"
    assert record.centered_on.index == 1
    assert record.descriptor == "none"
    assert record.permutation == 0
    assert record.controlling_atoms == [0, 2, 3]
    assert repr(record) == (
        "PotentialStereoInfo(stereo_type='atom_tetrahedral', "
        "specified='unspecified', centered_on=PotentialStereoCenter(kind='atom', index=1), "
        "descriptor='none', permutation=0)"
    )
    assert repr(analysis) == "PotentialStereoResult(records=1)"



def test_potential_stereo_default_canonical_cleanup_option_and_ordered_values():
    params=cosmolkit.PotentialStereoParams()
    assert params.clean is False
    assert params.flag_possible is True
    assert params.allow_nontetrahedral is True
    molecule=cosmolkit.Molecule.from_smiles("FC(Cl)Br")
    before=molecule.to_smiles()
    result=molecule.potential_stereo()
    explicit=molecule.potential_stereo_with_params(params)
    assert result.cleaned_molecule is None and explicit.cleaned_molecule is None
    assert len(result)==len(explicit)==1
    assert len(result.atom_ranks)==4
    assert result.ring_relations==[]
    assert repr(result.stereo[0])==repr(explicit.stereo[0])
    record=result.stereo[0]
    assert isinstance(record.stereo_type,cosmolkit.PotentialStereoType)
    assert isinstance(record.specified,cosmolkit.PotentialStereoSpecified)
    assert isinstance(record.descriptor,cosmolkit.PotentialStereoDescriptor)
    assert isinstance(record.centered_on,cosmolkit.PotentialStereoCenter)
    assert record.controlling_atoms==[0,2,3]
    copied=record.controlling_atoms
    copied.append(99)
    assert record.controlling_atoms==[0,2,3]
    assert molecule.to_smiles()==before
    params.clean = True
    assert params.clean is True
    for object_,field in [(record,"permutation"),(record.centered_on,"index"),(result,"stereo")]:
        with pytest.raises(AttributeError):setattr(object_,field,1)


def test_complete_potential_stereo_vocabulary_projection():
    for cls,values in [
        # Pinned RDKit Chirality.h StereoType includes Bond_Atropisomer.
        (cosmolkit.PotentialStereoType,["atom_tetrahedral","atom_squareplanar","atom_trigonalbipyramidal","atom_octahedral","bond_double","bond_cumulene_even","bond_atropisomer"]),
        (cosmolkit.PotentialStereoSpecified,["unspecified","specified","unknown"]),
        (cosmolkit.PotentialStereoDescriptor,["none","tetrahedral_clockwise","tetrahedral_counterclockwise","bond_cis","bond_trans","bond_atrop_cw","bond_atrop_ccw"]),
    ]:
        assert issubclass(cls, str)
        assert issubclass(cls, enum.Enum)
        assert [v.value for v in cls]==values
        for member in cls:
            assert cls(member.value) is member
            assert member == member.value
            assert str(member) == member.value
            assert format(member) == member.value
            assert hash(member) == hash(member.value)
            assert getattr(cls,member.value.upper()) is member
        with pytest.raises(ValueError):cls("undefined")


@pytest.mark.parametrize("ring", ["C[C@H]1CCCC[C@H]1C", "C[C@H]1CCC[C@H]1C"])
def test_potential_stereo_public_error_retains_property_kind_and_source_atomicity(ring):
    # The canonical CX string value reaches the source strict-bool ring cache
    # read. A malformed present property must remain a typed failure.
    molecule = cosmolkit.Molecule.from_smiles_with_params(
        ring + " |atomProp:1._ringStereochemCand.malformed|",
        cosmolkit.SmilesParseParams(sanitize=False, remove_hs=False, skip_cleanup=True),
    )
    params = cosmolkit.SmilesWriteParams(clean_stereo=False)
    before = molecule.to_smiles_with_params(params)
    tags = [atom.chiral_tag_code() for atom in molecule.atoms()]
    assert issubclass(cosmolkit.PotentialStereoError, ValueError)
    for clean in [False, True]:
        with pytest.raises(cosmolkit.OperationError) as caught:
            molecule.potential_stereo_with_params(cosmolkit.PotentialStereoParams(clean=clean))
        error = caught.value
        assert error.domain == "operation" and error.kind == "PotentialStereo"
        cause = error.__cause__
        assert isinstance(cause, cosmolkit.PotentialStereoError)
        assert cause.domain == "stereo" and cause.kind == "InvalidPropertyKind"
        assert cause.atom == 1
        assert cause.property == "_ringStereochemCand"
        assert cause.property_kind == "String"
        assert str(cause) == "atom 1 property _ringStereochemCand has invalid kind String"
        assert cause.__cause__ is None
        assert molecule.to_smiles_with_params(params) == before
        assert [atom.chiral_tag_code() for atom in molecule.atoms()] == tags
