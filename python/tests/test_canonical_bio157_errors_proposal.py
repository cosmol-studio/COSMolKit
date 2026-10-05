"""Canonical typed cause proposals, pending independent acceptance."""
import pytest
import cosmolkit as ck
from test_canonical_bio157_bridge_proposal import atom
from test_canonical_bio157_protocols_proposal import PDB

@pytest.mark.parametrize("reader,error",[(ck.BioStructure.from_pdb,ck.BioPdbReadError),(ck.Protein.from_pdb,ck.ProteinReadError)])
def test_pdb_stages_and_protein_outer_error_are_preserved(reader,error):
    with pytest.raises(error) as observed: reader("ATOM\n")
    assert observed.value.domain=="bio"
    if error is ck.BioPdbReadError:
        assert observed.value.line_number==1
        assert observed.value.stage
    else:
        assert isinstance(observed.value.__cause__,ck.BioPdbReadError)
        assert observed.value.__cause__.line_number==1

def test_protein_general_loader_preserves_its_outer_type():
    with pytest.raises(ck.ProteinReadError) as observed:ck.Protein.from_text("data_x\n_atom_site.id\n")
    assert observed.value.kind=="Structure"
    assert isinstance(observed.value.__cause__,ck.BioReadError)

def test_file_oserror_retains_path_errno_and_domain(tmp_path):
    missing=tmp_path/"missing.pdb"
    with pytest.raises(OSError) as observed:ck.BioStructure.read(missing)
    assert observed.value.filename==str(missing)
    assert observed.value.errno==2
    assert observed.value.domain=="bio"
    assert observed.value.kind=="Io"
    with pytest.raises(OSError) as written:ck.BioStructure.from_pdb(PDB).write_mmcif(tmp_path/"missing"/"output.cif")
    assert written.value.domain=="bio"
    assert written.value.kind=="FileWrite"

def test_conversion_has_a_typed_sanitize_cause_and_retains_input():
    text="".join(atom(i," C  ","C",10.*i) for i in range(1,7))+"CONECT    1    2    3    4    5\nCONECT    1    6\nEND\n"
    source=ck.BioStructure.from_pdb(text)
    with pytest.raises(ck.BioMoleculeError) as observed:source.to_molecule(remove_hs=False,proximity_bonding=False)
    cause=observed.value.__cause__
    assert isinstance(cause,ck.BioMoleculeConversionError)
    assert cause.domain=="bio"
    assert cause.kind=="Sanitize"
    assert cause.__cause__ is not None
    assert source.num_atoms()==6
