import pytest
import cosmolkit as ck

IDENTIFIER = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
KEY = "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"


def test_inchi_roundtrip_parameters_and_errors():
    molecule = ck.Molecule.from_smiles("CCO")
    read = ck.InchiReadParams()
    write = ck.InchiWriteParams()
    assert molecule.to_inchi() == IDENTIFIER
    assert molecule.to_inchi_key() == KEY
    assert ck.inchi_to_key(IDENTIFIER) == KEY
    write.options = "-AuxNone"
    assert molecule.to_inchi(write) == IDENTIFIER
    assert molecule.to_inchi(options="-AuxNone") == IDENTIFIER
    assert molecule.to_inchi_with_params(write) == IDENTIFIER
    assert molecule.to_inchi_key_with_params(write) == KEY
    for sanitize in (False, True):
        for remove in (False, True):
            read.sanitize = sanitize
            read.remove_hs = remove
            assert ck.Molecule.from_inchi(IDENTIFIER, read).to_inchi() == IDENTIFIER
            assert ck.Molecule.from_inchi_with_params(IDENTIFIER, read).to_inchi() == IDENTIFIER
            assert ck.Molecule.from_inchi(IDENTIFIER, sanitize=sanitize, remove_hs=remove).to_inchi() == IDENTIFIER
    with pytest.raises(TypeError):
        molecule.to_inchi(write, options="")
    with pytest.raises(TypeError):
        ck.Molecule.from_inchi(IDENTIFIER, read, sanitize=False)
    with pytest.raises(ck.InchiError) as raised:
        ck.Molecule.from_inchi("not InChI")
    assert raised.value.domain == "inchi"
    assert isinstance(raised.value.kind, ck.InchiErrorKind)
    assert molecule.to_smiles() == "CCO"
