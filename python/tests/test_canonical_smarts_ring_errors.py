"""Pinned CloseMolRings errors retain exact diagnostics and structural context."""
import pytest
import cosmolkit

@pytest.mark.parametrize("smarts,kind,fields,diagnostic", [
    ("[C]11", "SelfRingClosure", {"ring": 1, "atom": 0}, "SMARTS parse error: duplicated ring closure 1 bonds atom 0 to itself"),
    ("[C]1[C]1", "DuplicateRingBond", {"ring": 1, "begin_atom": 0, "end_atom": 1}, "SMARTS parse error: ring closure 1 duplicates bond between atom 0 and atom 1"),
])
def test_source_ring_syntax_errors(smarts, kind, fields, diagnostic):
    for constructor in (cosmolkit.parse_smarts, cosmolkit.QueryGraph.from_smarts):
        with pytest.raises(cosmolkit.SmartsParseError) as caught:
            constructor(smarts)
        error = caught.value
        assert str(error) == diagnostic
        assert error.kind == kind
        assert error.domain == "search"
        for name, value in fields.items():
            assert getattr(error, name) == value


def test_source_unclosed_ring_preserves_separate_parse_error():
    with pytest.raises(cosmolkit.SmartsParseError) as caught:
        cosmolkit.parse_smarts("C1CC")
    assert caught.value.kind == "Parse"
    assert str(caught.value) == "SMARTS parse error: unclosed ring"
