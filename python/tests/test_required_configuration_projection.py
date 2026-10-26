"""Required inputs and mutually exclusive native configuration call forms."""
import cosmolkit as ck
import pytest


def test_required_config_object_and_keywords_reach_one_actual_operation():
    mol = ck.mol_from_smiles("CCO")
    params = ck.FragmentSmilesWriteParams(atoms=[0, 1])
    assert mol.to_fragment_smiles(params) == mol.to_fragment_smiles(atoms=[0, 1]) == "CC"
    params.atoms = [2]
    assert mol.to_fragment_smiles(params) == "O"


@pytest.mark.parametrize("keywords", [{}, {"unknown": 1}, {"atoms": [0], "unknown": 1}])
def test_required_config_missing_and_unknown_fields_are_rejected(keywords):
    with pytest.raises(TypeError):
        ck.FragmentSmilesWriteParams(**keywords)


def test_required_config_mixed_fields_and_instance_are_rejected():
    mol = ck.mol_from_smiles("CCO")
    params = ck.FragmentSmilesWriteParams(atoms=[0])
    with pytest.raises(TypeError):
        mol.to_fragment_smiles(params, bonds=None)
