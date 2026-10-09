"""Original scalar conditions through canonical default/Params projections."""
import cosmolkit
import pytest


def test_pattern_scalar_defaults_options_type_and_exact_bits():
    molecule = cosmolkit.Molecule.from_smiles("CC")
    default = molecule.fingerprint_pattern()
    explicit = molecule.fingerprint_pattern_with_params(
        cosmolkit.PatternFingerprintParams(n_bits=2048, tautomeric=False)
    )
    tautomeric = molecule.fingerprint_pattern_with_params(
        cosmolkit.PatternFingerprintParams(n_bits=2048, tautomeric=True)
    )
    assert isinstance(default, cosmolkit.Fingerprint)
    assert default.n_bits() == 2048
    assert default.on_bits() == [429, 778, 1022, 1061, 1236, 1295]
    assert explicit.on_bits() == default.on_bits()
    assert tautomeric.on_bits() == [429, 776, 778, 1022, 1061, 1236, 1295]


def test_pattern_scalar_preserves_input_and_reports_argument_errors():
    molecule = cosmolkit.Molecule.from_smiles("c1ccccc1O")
    before = molecule.to_smiles()
    first = molecule.fingerprint_pattern_with_params(
        cosmolkit.PatternFingerprintParams(n_bits=127, tautomeric=True)
    )
    second = molecule.fingerprint_pattern_with_params(
        cosmolkit.PatternFingerprintParams(n_bits=127, tautomeric=True)
    )
    assert first.on_bits() == second.on_bits()
    assert molecule.to_smiles() == before
    with pytest.raises(ValueError, match="fingerprint requires n_bits > 0"):
        molecule.fingerprint_pattern_with_params(cosmolkit.PatternFingerprintParams(n_bits=0))
    with pytest.raises(TypeError):
        cosmolkit.PatternFingerprintParams(n_bits="2048")
    with pytest.raises(TypeError):
        cosmolkit.PatternFingerprintParams(tautomeric="yes")
