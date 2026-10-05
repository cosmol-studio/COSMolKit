"""Explicit fingerprint values, source similarity and platform length bounds."""
import sys
import pytest
import cosmolkit as ck

@pytest.mark.parametrize("width,bits,expected",[(0,[],[]),(1,[0,0],[0]),(65,[64,0,63,64,1],[0,1,63,64])])
def test_explicit_constructor_fixed_width_sort_dedup_and_owned_results(width,bits,expected):
    original=bits.copy()
    value=ck.Fingerprint.from_on_bits(width,bits)
    assert value.n_bits()==len(value)==width
    assert value.on_bits()==expected
    assert bits==original
    detached=value.on_bits();detached.clear()
    assert value.on_bits()==expected

@pytest.mark.parametrize("width,bits,index",[(0,[0],0),(7,[7,8],7),((1<<32)-1,[(1<<32)-1],(1<<32)-1)])
def test_first_out_of_range_code_keeps_source_value_error_and_context(width,bits,index):
    with pytest.raises(ck.FingerprintError) as raised:ck.Fingerprint.from_on_bits(width,bits)
    assert isinstance(raised.value,ValueError)
    assert raised.value.domain=="fingerprints"
    assert raised.value.kind=="SparseIndexOutOfRange"
    assert raised.value.index==index and raised.value.size==width

def test_exact_dense_tanimoto_and_original_empty_zero():
    a=ck.Fingerprint.from_on_bits(65,[0,1,64]);b=ck.Fingerprint.from_on_bits(65,[1,63,64])
    assert a.tanimoto(b)==b.tanimoto(a)==0.5
    assert a.tanimoto(a)==1.0
    empty=ck.Fingerprint.from_on_bits(65,[])
    assert empty.tanimoto(empty)==empty.tanimoto(a)==0.0
    with pytest.raises(ck.FingerprintError) as raised:empty.tanimoto(ck.Fingerprint.from_on_bits(64,[]))
    assert raised.value.kind=="BitLengthMismatch"
    assert (raised.value.left,raised.value.right)==(65,64)

@pytest.mark.parametrize("length",[0,16,(1<<32)-1,sys.maxsize])
def test_sparse_count_len_is_vector_length_without_cloning_entries(length):
    value=ck.SparseCountFingerprint.new(length)
    assert len(value)==value.length()==length
    if length: value.set_value(0,3)
    assert len(value)==length

def test_length_and_unsigned_constructor_transport_overflow_are_explicit():
    for length in [sys.maxsize+1,(1<<64)-1]:
        with pytest.raises(OverflowError) as raised:
            len(ck.SparseCountFingerprint.new(length))
        # Original usize conversion has its own error only above usize::MAX;
        # pinned PyO3 callback.rs130 uses PyOverflowError::new_err(()) for the
        # signed Py_ssize_t slot conversion, hence an empty message otherwise.
        expected = "fingerprint size exceeds Python platform size" if length > 2*sys.maxsize+1 else ""
        assert str(raised.value) == expected
    for width,bits in [(-1,[]),(1<<32,[]),(7,[-1]),(7,[1<<32])]:
        with pytest.raises(OverflowError):ck.Fingerprint.from_on_bits(width,bits)
    with pytest.raises(TypeError):ck.Fingerprint.from_on_bits(7,["0"])
