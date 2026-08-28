import random
from collections.abc import Iterable, Iterator
from typing import override

import cosmolkit
import pytest


def smiles(isomers: Iterable[cosmolkit.Molecule]) -> list[str]:
    return [molecule.to_smiles() for molecule in isomers]


def test_options_defaults_mutability_and_repr_use_project_native_names():
    options = cosmolkit.StereoisomerOptions()

    assert options.try_embedding is False
    assert options.only_unassigned is True
    assert options.only_stereo_groups is False
    assert options.max_isomers == 1024
    assert options.rand is None
    assert options.unique is True
    assert repr(options) == (
        "StereoisomerOptions(try_embedding=false, only_unassigned=true, "
        "max_isomers=1024, rand=None, unique=true, only_stereo_groups=false)"
    )

    options.try_embedding = True
    options.only_unassigned = False
    options.only_stereo_groups = True
    options.max_isomers = 7
    options.rand = 61453
    options.unique = False

    assert options.try_embedding is True
    assert options.only_unassigned is False
    assert options.only_stereo_groups is True
    assert options.max_isomers == 7
    assert options.rand == 61453
    assert options.unique is False
    assert "rand=..." in repr(options)


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        ("CCO", ["CCO"]),
        ("FC(Cl)Br", ["F[C@H](Cl)Br", "F[C@@H](Cl)Br"]),
        ("FC=CF", ["F/C=C/F", "F/C=C\\F"]),
        (
            "CC(F)C(Cl)Br",
            [
                "C[C@H](F)[C@H](Cl)Br",
                "C[C@@H](F)[C@H](Cl)Br",
                "C[C@H](F)[C@@H](Cl)Br",
                "C[C@@H](F)[C@@H](Cl)Br",
            ],
        ),
    ],
)
def test_focused_outputs_match_pinned_rdkit_exactly(source: str, expected: list[str]):
    molecule = cosmolkit.Molecule.from_smiles(source)
    before = molecule.to_smiles()

    iterator = molecule.stereoisomers()
    assert isinstance(iterator, cosmolkit.StereoisomerIterator)
    assert iterator.yielded_count == 0
    assert repr(iterator) == "StereoisomerIterator(yielded_count=0)"

    assert smiles(iterator) == expected
    assert iterator.yielded_count == len(expected)
    assert list(iterator) == []
    assert molecule.to_smiles() == before


def test_iterator_is_lazy_and_no_center_yields_one_isolated_value():
    molecule = cosmolkit.Molecule.from_smiles("CCO")
    iterator = molecule.stereoisomers()

    first = next(iterator)
    assert first is not molecule
    assert first.to_smiles() == "CCO"
    assert iterator.yielded_count == 1
    with pytest.raises(StopIteration):
        next(iterator)


def test_count_is_exact_for_arbitrary_center_counts_and_options():
    source = "Br" + "[CH](Cl)" * 70 + "F"
    molecule = cosmolkit.Molecule.from_smiles(source)

    assert molecule.stereoisomer_count() == 1 << 70
    assert isinstance(molecule.stereoisomer_count(), int)

    assigned = cosmolkit.Molecule.from_smiles("F[C@H](Cl)C(F)Br")
    assert assigned.stereoisomer_count() == 2
    assert assigned.stereoisomer_count(
        cosmolkit.StereoisomerOptions(only_unassigned=False)
    ) == 4


def test_only_unassigned_unique_and_only_stereo_groups_match_source_branches():
    assigned = cosmolkit.Molecule.from_smiles("F[C@H](Cl)C(F)Br")
    assert len(list(assigned.stereoisomers())) == 2
    assert len(
        list(
            assigned.stereoisomers(
                cosmolkit.StereoisomerOptions(only_unassigned=False)
            )
        )
    ) == 4

    symmetric = cosmolkit.Molecule.from_smiles("FC(Cl)C=CC=CC(F)Cl")
    assert len(list(symmetric.stereoisomers())) == 10
    assert len(
        list(symmetric.stereoisomers(cosmolkit.StereoisomerOptions(unique=False)))
    ) == 16

    ungrouped = cosmolkit.Molecule.from_smiles("FC(Cl)Br")
    grouped_only = cosmolkit.StereoisomerOptions(only_stereo_groups=True)
    assert smiles(ungrouped.stereoisomers(grouped_only)) == ["FC(Cl)Br"]


def test_seeded_default_and_custom_random_sources_match_pinned_rdkit_sequences():
    source = "Br" + "[CH](Cl)" * 6 + "F"
    molecule = cosmolkit.Molecule.from_smiles(source)
    seeded_expected = [
        "F[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)Br",
        "F[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)Br",
        "F[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)Br",
    ]
    default_expected = [
        "F[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)Br",
        "F[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)Br",
        "F[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)Br",
    ]

    assert smiles(
        molecule.stereoisomers(
            cosmolkit.StereoisomerOptions(max_isomers=3, rand=0xF00D)
        )
    ) == seeded_expected
    assert smiles(
        molecule.stereoisomers(
            cosmolkit.StereoisomerOptions(
                max_isomers=3,
                rand=random.Random(0xF00D),
            )
        )
    ) == seeded_expected
    assert smiles(
        molecule.stereoisomers(cosmolkit.StereoisomerOptions(max_isomers=3))
    ) == default_expected

    class ScriptedRandom(random.Random):
        def __init__(self) -> None:
            super().__init__(0)
            self.values: Iterator[int] = iter([0, 1, 2])
            self.requested_widths: list[int] = []

        @override
        def getrandbits(self, width: int) -> int:
            self.requested_widths.append(width)
            return next(self.values)

    scripted = ScriptedRandom()
    scripted_expected = [
        "F[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)Br",
        "F[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)Br",
        "F[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)Br",
    ]
    assert smiles(
        molecule.stereoisomers(
            cosmolkit.StereoisomerOptions(max_isomers=3, rand=scripted)
        )
    ) == scripted_expected
    assert scripted.requested_widths == [6, 6, 6]


def test_custom_random_errors_remain_visible_at_lazy_iteration_boundary():
    class FailingRandom(random.Random):
        @override
        def getrandbits(self, width: int) -> int:
            raise RuntimeError(f"fixture random failure for {width} bits")

    molecule = cosmolkit.Molecule.from_smiles("Br" + "[CH](Cl)" * 6 + "F")
    iterator = molecule.stereoisomers(
        cosmolkit.StereoisomerOptions(max_isomers=1, rand=FailingRandom())
    )

    with pytest.raises(ValueError, match="fixture random failure for 6 bits"):
        next(iterator)
    assert iterator.yielded_count == 0


def test_potential_stereo_analysis_returns_typed_isolated_state():
    molecule = cosmolkit.Molecule.from_smiles("FC(Cl)Br")
    before = molecule.to_smiles()

    analysis = molecule.analyze_potential_stereo()
    assert isinstance(analysis, cosmolkit.PotentialStereoAnalysis)
    assert len(analysis) == 1
    assert analysis.molecule is not molecule
    assert analysis.molecule.to_smiles() == before
    assert molecule.to_smiles() == before

    record = analysis.stereo_info[0]
    assert isinstance(record, cosmolkit.PotentialStereoInfo)
    assert record.stereo_type == "atom_tetrahedral"
    assert record.specified == "unspecified"
    assert record.center_kind == "atom"
    assert record.center_index == 1
    assert record.descriptor == "none"
    assert record.permutation == 0
    assert record.controlling_atoms == [0, 2, 3]
    assert repr(record) == (
        "PotentialStereoInfo(stereo_type='atom_tetrahedral', "
        "specified='unspecified', center_kind='atom', center_index=1, "
        "descriptor='none', permutation=0)"
    )
    assert repr(analysis) == "PotentialStereoAnalysis(records=1)"


def test_enumeration_composes_with_hydrogen_conformer_and_descriptor_apis():
    source = cosmolkit.Molecule.from_smiles("FC(Cl)Br").with_hydrogens()
    parameters = cosmolkit.EmbedParameters.etkdg_v3()
    parameters.random_seed = 42
    parameters.num_threads = 1
    source = source.with_3d_conformer(parameters)
    before_smiles = source.to_smiles()
    before_mass = cosmolkit.calc_mol_wt(source)
    before_coordinates = source.coordinates_3d().copy()

    outputs = list(source.stereoisomers())
    assert [molecule.to_smiles() for molecule in outputs] == [
        "[H][C@@](F)(Cl)Br",
        "[H][C@](F)(Cl)Br",
    ]
    assert all(len(molecule) == len(source) for molecule in outputs)
    assert all(molecule.num_conformers() == 1 for molecule in outputs)
    assert all(cosmolkit.calc_mol_wt(molecule) == before_mass for molecule in outputs)

    assert source.to_smiles() == before_smiles
    assert source.num_conformers() == 1
    assert (source.coordinates_3d() == before_coordinates).all()
