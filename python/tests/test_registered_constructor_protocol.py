"""Native constructor protocols declared by the canonical Rust registry."""

import cosmolkit as ck
import pytest


@pytest.mark.parametrize(
    "name,args,kwargs",
    [
        ("MorganBondInvariantsGenerator", (), {}),
        ("MorganCallParams", (), {}),
        ("MorganFingerprintGenerator", (), {}),
        ("TautomerScoreTerm", ("carbon", "[#6]", 3), {}),
        ("FingerprintAdditionalOutput", (), {}),
        ("AtomSpec", (ck.Element.C,), {}),
        ("BondSpec", (0, 1, ck.BondOrder.SINGLE), {}),
        ("TopologicalTorsionFingerprintGenerator", (), {}),
        ("LegacyTopologicalTorsionParams", (), {}),
        ("AlignmentAtomMap", (0, 1), {}),
        ("AlignmentParameters", (), {}),
        ("BestAlignmentParameters", (), {}),
        ("CoordinateRmsdParameters", (), {}),
        ("AllConformerRmsdParameters", (), {}),
        ("ConformerAlignmentParameters", (), {}),
        ("EmbedParams", (), {}),
    ],
)
def test_exact_registered_new_is_a_native_class_factory(name, args, kwargs):
    kind = getattr(ck, name)
    assert callable(kind.new)
    value = kind.new(*args, **kwargs)
    assert type(value) is kind
    assert type(kind(*args, **kwargs)) is kind


def test_builder_factories_retain_ids_bond_order_and_atom_metadata():
    builder = ck.MoleculeBuilder.new()
    carbon = builder.add_atom(ck.AtomSpec.new(ck.Element.C).with_atom_map(17))
    oxygen = builder.add_atom(ck.AtomSpec.new(ck.Element.O).with_isotope(18))
    builder.add_bond(ck.BondSpec.new(carbon, oxygen, ck.BondOrder.DOUBLE))
    molecule = builder.build()
    assert [atom.atomic_number() for atom in molecule.atoms()] == [6, 8]
    assert molecule.atoms()[0].atom_map() == 17
    assert molecule.atoms()[1].isotope() == 18
    assert molecule.bonds()[0].order() == ck.BondOrder.DOUBLE
    with pytest.raises(ValueError, match="invalid BondOrder code"):
        ck.BondSpec.new(0, 1, 999)


def test_fingerprint_output_factories_do_not_share_allocation_state():
    first = ck.FingerprintAdditionalOutput.new()
    second = ck.FingerprintAdditionalOutput.new()
    assert first.atom_counts() is None
    assert second.atom_counts() is None
    first.allocate_atom_counts()
    assert first.atom_counts() == []
    assert second.atom_counts() is None


def test_score_term_factory_preserves_supplied_value():
    value = ck.TautomerScoreTerm.new("label", "[#7]", -9)
    assert value.name() == "label"
    assert value.smarts() == "[#7]"
    assert value.score() == -9
    assert value == ck.TautomerScoreTerm("label", "[#7]", -9)


def test_generator_factories_pass_explicit_configuration_to_the_owner():
    morgan = ck.MorganFingerprintGenerator.new(params=ck.MorganParams(radius=3))
    assert morgan.settings().params().radius == 3
    torsion = ck.TopologicalTorsionFingerprintGenerator.new(
        params=ck.TopologicalTorsionParams(torsion_atom_count=5)
    )
    assert torsion.settings().params().torsion_atom_count == 5
    with pytest.raises(OverflowError):
        ck.LegacyTopologicalTorsionParams.new(torsion_atom_count=-1)


def test_alignment_factories_preserve_explicit_maps_and_options():
    atom_map = ck.AlignmentAtomMap.new(4, 2)
    assert (atom_map.probe_atom, atom_map.reference_atom) == (4, 2)
    value = ck.AlignmentParameters.new(
        probe_conformer_id=7,
        reference_conformer_id=9,
        atom_map=[atom_map],
        weights=[2.5],
        reflect=True,
        max_iterations=19,
    )
    assert value.probe_conformer_id == 7
    assert value.reference_conformer_id == 9
    assert value.atom_map[0].probe_atom == 4
    assert value.atom_map[0].reference_atom == 2
    assert value.weights == [2.5]
    assert value.reflect is True
    assert value.max_iterations == 19
    with pytest.raises(OverflowError):
        ck.AlignmentAtomMap.new(-1, 0)


def test_embed_new_projects_the_default_rust_factory():
    value = ck.EmbedParams.new()
    assert value.to_json() == ck.EmbedParams().to_json()
    updated = value.with_json('{"randomSeed":23,"numThreads":2}')
    assert updated.to_json() != value.to_json()
    assert value.to_json() == ck.EmbedParams.new().to_json()
