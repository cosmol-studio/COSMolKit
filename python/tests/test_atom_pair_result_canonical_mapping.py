"""Atom-pair result semantics through canonical owned values.

Fingerprint results and optional additional output preserve independent values
and detached metadata snapshots. Historical wrapper identities differ.
"""
import cosmolkit as ck

BITS = [624, 1144, 1336, 1337, 1404, 1596]
PAIRS = {624: [(0, 3)], 1144: [(0, 2)], 1336: [(1, 2), (1, 3)],
         1337: [(1, 2), (1, 3)], 1404: [(0, 1)], 1596: [(2, 3)]}
TO_BITS = [[624, 1144, 1404], [1336, 1337, 1404],
           [1144, 1336, 1337, 1596], [624, 1336, 1337, 1596]]
ATOMS = {bit: [list(pair) for pair in pairs] for bit, pairs in PAIRS.items()}

def collect():
    output = ck.FingerprintAdditionalOutput()
    for method in ('allocate_atom_counts', 'allocate_atom_to_bits',
                   'allocate_bit_info_map', 'allocate_atoms_per_bit'):
        getattr(output, method)()
    result = ck.Molecule.from_smiles('CCCO').fingerprint_atom_pair_with_params(
        ck.AtomPairFingerprintParams(), output)
    return result, output

def fields(output):
    return (output.atom_counts(), output.atom_to_bits(), output.bit_info_map(),
            output.atoms_per_bit(), output.bit_paths())

def test_original_fingerprint_getter_independent_owned_value_mapping():
    result, output = collect()
    assert isinstance(result, ck.Fingerprint)
    assert result.n_bits() == 2048
    assert result.on_bits() == BITS
    detached = result.on_bits()
    detached.clear()
    assert result.on_bits() == BITS
    second = ck.Molecule.from_smiles('CCCO').fingerprint_atom_pair()
    assert second is not result
    assert second.on_bits() == result.on_bits() == BITS
    assert not hasattr(ck, 'AtomPairFingerprintResult')

def test_original_additional_output_getter_collected_and_absent_mapping():
    result, output = collect()
    assert fields(output) == ([3, 3, 3, 3], TO_BITS, PAIRS, ATOMS, None)
    snapshot = fields(output)
    snapshot[0][0] = 0
    snapshot[1][0].clear()
    snapshot[2][624].clear()
    snapshot[3][624][0].clear()
    assert fields(output) == ([3, 3, 3, 3], TO_BITS, PAIRS, ATOMS, None)
    absent = None
    plain = ck.Molecule.from_smiles('CCCO').fingerprint_atom_pair_with_params(
        ck.AtomPairFingerprintParams(), absent)
    assert plain.on_bits() == BITS
    assert absent is None
    assert fields(ck.FingerprintAdditionalOutput()) == (None, None, None, None, None)
    # Old wrapper ValueError accessor is replaced by explicit optional caller
    # output; literal old error identity is not a pass claimed by this mapping.

def test_original_repr_width_and_collected_flag_information_mapping():
    result, output = collect()
    assert repr(result) == 'Fingerprint(n_bits=2048)'
    assert repr(output) == ('FingerprintAdditionalOutput(atom_to_bits=true, '
        'bit_info_map=true, bit_paths=false, atom_counts=true, atoms_per_bit=true)')
    empty = ck.FingerprintAdditionalOutput()
    assert repr(empty) == ('FingerprintAdditionalOutput(atom_to_bits=false, '
        'bit_info_map=false, bit_paths=false, atom_counts=false, atoms_per_bit=false)')
    # The old single wrapper repr string does not exist. Its width and
    # collection information are observable through two canonical values.

def test_original_result_snapshot_remains_independent_after_caller_AO_reuse():
    first, output = collect()
    snapshot = fields(output)
    first_bits = first.on_bits()
    second = ck.Molecule.from_smiles('CCO').fingerprint_atom_pair_with_params(
        ck.AtomPairFingerprintParams(), output)
    assert output.atom_counts() == [2, 2, 2]
    assert len(output.atom_to_bits()) == 3
    assert snapshot == ([3, 3, 3, 3], TO_BITS, PAIRS, ATOMS, None)
    assert first.on_bits() == first_bits == BITS
    assert second is not first


def test_caller_output_presence_distinguishes_unallocated_collected_from_absent():
    # Old has_additional_output observed Option presence, independently of
    # which fields were allocated. Preserve both facts using existing values.
    output = ck.FingerprintAdditionalOutput()
    fingerprint = ck.Molecule.from_smiles('CCCO').fingerprint_atom_pair_with_params(
        ck.AtomPairFingerprintParams(), output)
    assert output is not None
    assert fields(output) == (None, None, None, None, None)
    assert fingerprint.on_bits() == BITS
    absent = None
    second = ck.Molecule.from_smiles('CCCO').fingerprint_atom_pair_with_params(
        ck.AtomPairFingerprintParams(), absent)
    assert absent is None
    assert second.on_bits() == fingerprint.on_bits()
