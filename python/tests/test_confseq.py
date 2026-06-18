from __future__ import annotations

import cosmolkit


def test_confseq_submodule_is_importable():
    import cosmolkit.confseq as confseq

    assert confseq.decode is cosmolkit.confseq.decode


def test_confseq_decode_returns_cosmolkit_molecule():
    mol = cosmolkit.confseq.decode("C C", "C C")

    assert isinstance(mol, cosmolkit.Molecule)
    assert mol.num_atoms() == 2
    assert mol.num_conformers() == 1


def test_confseq_decode_can_disable_uff_optimization():
    mol = cosmolkit.confseq.decode("C C", "C C", optimize_with_uff=False)

    assert isinstance(mol, cosmolkit.Molecule)
    assert mol.num_atoms() == 2
    assert mol.num_conformers() == 1


def test_confseq_decode_batch_preserves_order_and_uses_local_cache():
    mols = cosmolkit.confseq.decode_batch(
        ["C C", "C C"],
        ["C C", "C C"],
        n_jobs=1,
        optimize_with_uff=False,
    )

    assert [mol.num_atoms() for mol in mols] == [2, 2]
    assert [mol.num_conformers() for mol in mols] == [1, 1]


def test_confseq_decode_batch_can_keep_errors():
    mols = cosmolkit.confseq.decode_batch(
        ["C C", "not smiles"],
        ["C C", "not smiles"],
        errors="keep",
    )

    assert isinstance(mols[0], cosmolkit.Molecule)
    assert mols[1] is None
