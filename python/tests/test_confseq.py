from __future__ import annotations

import importlib

import pytest

import cosmolkit

CORPUS_CONFSEQ = (
    "N C <112> | ( = O ) <84> C <111> | <-45> n 1 c c ( <21> N <123> | "
    "<173> C <116> | ( = O ) <120> c 2 c c c c c 2 <112> C <112> | <-66> "
    "N 2 <-172> C ( = O ) <-2> C <0> N <3> C <174> 2 = O ) c n 1"
)


def test_confseq_submodule_is_importable():
    confseq = importlib.import_module("cosmolkit.confseq")

    assert confseq.decode is cosmolkit.confseq.decode


def test_confseq_decode_returns_cosmolkit_molecule():
    mol = cosmolkit.confseq.decode("C C")

    assert isinstance(mol, cosmolkit.Molecule)
    assert mol.num_atoms() == 2
    assert mol.num_conformers() == 1


def test_confseq_decode_can_disable_uff_optimization():
    mol = cosmolkit.confseq.decode("C C", optimize_with_uff=False)

    assert isinstance(mol, cosmolkit.Molecule)
    assert mol.num_atoms() == 2
    assert mol.num_conformers() == 1


def test_confseq_decode_accepts_explicit_template_backend():
    dg = cosmolkit.confseq.decode(
        CORPUS_CONFSEQ,
        optimize_with_uff=False,
        template_backend="distance_geometry",
    )
    fast = cosmolkit.confseq.decode(
        CORPUS_CONFSEQ,
        optimize_with_uff=False,
        template_backend="fast_geometry",
    )

    assert isinstance(dg, cosmolkit.Molecule)
    assert isinstance(fast, cosmolkit.Molecule)
    assert dg.num_atoms() == 26
    assert fast.num_atoms() == 26
    assert dg.num_conformers() == 1
    assert fast.num_conformers() == 1


def test_confseq_decode_rejects_unknown_template_backend():
    with pytest.raises(ValueError, match="template_backend"):
        cosmolkit.confseq.decode(
            "C C",
            optimize_with_uff=False,
            template_backend="unknown",
        )


def test_confseq_decode_rejects_old_base_conformer_backend_name():
    with pytest.raises(ValueError, match="template_backend"):
        cosmolkit.confseq.decode(
            "C <90> | C C",
            optimize_with_uff=False,
            template_backend="base_conformer",
        )


def test_confseq_decode_batch_preserves_order_and_uses_local_cache():
    mols = cosmolkit.confseq.decode_batch(
        ["C C", "C C"],
        n_jobs=1,
        optimize_with_uff=False,
        template_backend="distance_geometry",
    )

    assert all(isinstance(mol, cosmolkit.Molecule) for mol in mols)
    concrete_mols = [mol for mol in mols if mol is not None]
    assert [mol.num_atoms() for mol in concrete_mols] == [2, 2]
    assert [mol.num_conformers() for mol in concrete_mols] == [1, 1]


def test_confseq_decode_batch_can_keep_errors():
    mols = cosmolkit.confseq.decode_batch(
        ["C C", "not smiles"],
        errors="keep",
    )

    assert isinstance(mols[0], cosmolkit.Molecule)
    assert mols[1] is None
