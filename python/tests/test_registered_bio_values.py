"""Small live-binding regressions; no reference files or generated corpora."""

import struct

import pytest

import cosmolkit as ck
from test_bio_baseline_capabilities import PDB


def test_row_span_checks_native_overflow_and_preserves_bounds():
    span = ck.BioRowSpan(3, 4)
    assert (span.start(), span.len(), span.end(), span.is_empty()) == (3, 4, 7, False)
    assert ck.BioRowSpan(7, 0).is_empty()
    with pytest.raises(ck.BioStructureError) as caught:
        ck.BioRowSpan(2**32 - 1, 1)
    assert (caught.value.start, caught.value.len) == (2**32 - 1, 1)


def test_crystal_cell_retains_exact_detached_scalars():
    default = ck.BioCrystalCell()
    assert (default.a, default.b, default.c, default.alpha, default.beta, default.gamma) == (1, 1, 1, 90, 90, 90)
    cell = ck.BioCrystalCell(-0.0, 2.5, 3.25, 90, 100, 120)
    assert struct.pack("d", cell.a) == struct.pack("d", -0.0)
    assert (cell.b, cell.c, cell.alpha, cell.beta, cell.gamma) == (2.5, 3.25, 90, 100, 120)
    with pytest.raises(AttributeError):
        cell.a = 5


def test_sifts_and_entity_database_reference_retain_every_field():
    sifts = ck.BioSiftsUnpResidue(255, 254, 65535)
    assert (sifts.residue, sifts.accession_index, sifts.number) == (255, 254, 65535)
    assert ck.BioSiftsUnpResidue().residue is None
    ref = ck.BioEntityDbRef(
        db_name="UNP", accession_code="P123", id_code="ENTRY", isoform="2",
        seq_begin=ck.PdbSeqId(-3, ord("A")), seq_end=ck.PdbSeqId(9),
        db_begin=ck.PdbSeqId(1), db_end=ck.PdbSeqId(12, ord("B")),
        label_seq_begin=2, label_seq_end=8,
    )
    assert (ref.db_name, ref.accession_code, ref.id_code, ref.isoform) == ("UNP", "P123", "ENTRY", "2")
    assert [(x.seq_num(), x.ins_code()) for x in (ref.seq_begin, ref.seq_end, ref.db_begin, ref.db_end)] == [(-3, "A"), (9, None), (1, None), (12, "B")]
    assert (ref.label_seq_begin, ref.label_seq_end) == (2, 8)


def test_bio_kinds_are_real_typed_values_without_breaking_string_comparisons():
    structure = ck.BioStructure.from_pdb(PDB)
    # The protein view performs classification; detached input rows retain
    # the reader's original kind, rather than inventing classification here.
    assert isinstance(structure.chains()[0].kind(), ck.ChainKind)
    assert structure.protein().chains()[0].kind() is ck.ChainKind.Protein
    assert structure.protein().chains()[0].kind() == "Protein"
    assert str(structure.protein().chains()[0].kind()) == "Protein"
    residue = structure.residues()[0]
    assert residue.kind() is ck.ResidueKind.AminoAcid
    assert isinstance(residue.entity_kind(), ck.EntityKind)
    assert isinstance(structure.atoms()[0].calc_flag(), ck.BioCalcFlag)
    assert isinstance(residue.sifts_unp(), ck.BioSiftsUnpResidue)
    assert ck.PolymerKind.PeptideL == "PeptideL"
    for entity in structure.entities():
        assert isinstance(entity.kind(), ck.EntityKind)
        assert isinstance(entity.polymer_kind(), ck.PolymerKind)
        assert all(isinstance(ref, ck.BioEntityDbRef) for ref in entity.dbrefs())


def test_registered_read_error_is_retained_in_actual_exception_chain():
    with pytest.raises(Exception) as caught:
        ck.Molecule.from_xyz_block("2\ninvalid\nC 0 0 0\n")
    chain = []
    error = caught.value
    while error is not None:
        chain.append(error)
        error = error.__cause__
    assert any(isinstance(error, ck.XyzReadError) for error in chain), chain
