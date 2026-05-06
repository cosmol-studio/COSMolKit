import pytest

import cosmolkit

try:
    from rdkit import Chem
    from rdkit import DataStructs
    from rdkit.Chem import rdFingerprintGenerator
except Exception:  # pragma: no cover - import failure is reported by skip.
    Chem = None
    DataStructs = None
    rdFingerprintGenerator = None


pytestmark = pytest.mark.skipif(
    Chem is None or DataStructs is None or rdFingerprintGenerator is None,
    reason="RDKit is required for Morgan fingerprint parity tests",
)


def rdkit_on_bits(smiles, **generator_kwargs):
    assert Chem is not None
    assert rdFingerprintGenerator is not None
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    generator = rdFingerprintGenerator.GetMorganGenerator(**generator_kwargs)
    return list(generator.GetFingerprint(mol).GetOnBits())


def assert_morgan_bits_equal(smiles, **kwargs):
    assert Chem is not None
    assert rdFingerprintGenerator is not None
    rdkit_kwargs = {
        "radius": kwargs.get("radius", 2),
        "fpSize": kwargs.get("n_bits", 2048),
        "includeChirality": kwargs.get("include_chirality", False),
        "useBondTypes": kwargs.get("use_bond_types", True),
        "countSimulation": kwargs.get("count_simulation", False),
        "onlyNonzeroInvariants": kwargs.get("only_nonzero_invariants", False),
        "includeRedundantEnvironments": kwargs.get(
            "include_redundant_environments", False
        ),
    }
    if kwargs.get("count_bounds") is not None:
        rdkit_kwargs["countBounds"] = kwargs["count_bounds"]

    atom_generator = kwargs.get("atom_invariants_generator")
    if atom_generator == "feature":
        rdkit_kwargs["atomInvariantsGenerator"] = (
            rdFingerprintGenerator.GetMorganFeatureAtomInvGen()
        )
    elif atom_generator == "connectivity":
        rdkit_kwargs["atomInvariantsGenerator"] = (
            rdFingerprintGenerator.GetMorganAtomInvGen(
                kwargs.get("atom_invariants_include_ring_membership", True)
            )
        )

    bond_generator = kwargs.get("bond_invariants_generator")
    if bond_generator is not None:
        rdkit_kwargs["bondInvariantsGenerator"] = rdFingerprintGenerator.GetMorganBondInvGen(
            kwargs.get("bond_invariants_use_bond_types", True),
            kwargs.get("bond_invariants_use_chirality", False),
        )

    generator = rdFingerprintGenerator.GetMorganGenerator(**rdkit_kwargs)
    if "num_bits_per_feature" in kwargs:
        generator.GetOptions().numBitsPerFeature = kwargs["num_bits_per_feature"]

    rdmol = Chem.MolFromSmiles(smiles, sanitize=True)
    fp_kwargs = {}
    if kwargs.get("from_atoms") is not None:
        fp_kwargs["fromAtoms"] = kwargs["from_atoms"]
    if kwargs.get("ignore_atoms") is not None:
        fp_kwargs["ignoreAtoms"] = kwargs["ignore_atoms"]
    if kwargs.get("custom_atom_invariants") is not None:
        fp_kwargs["customAtomInvariants"] = kwargs["custom_atom_invariants"]
    if kwargs.get("custom_bond_invariants") is not None:
        fp_kwargs["customBondInvariants"] = kwargs["custom_bond_invariants"]

    expected = list(generator.GetFingerprint(rdmol, **fp_kwargs).GetOnBits())
    observed = cosmolkit.Molecule.from_smiles(smiles).fingerprint_morgan(**kwargs)
    assert observed.n_bits() == kwargs.get("n_bits", 2048)
    assert observed.on_bits() == expected


def test_morgan_fingerprint_default_and_tanimoto_match_rdkit():
    assert Chem is not None
    assert DataStructs is not None
    assert rdFingerprintGenerator is not None
    assert_morgan_bits_equal("CCO", radius=2, n_bits=2048)

    ethanol = cosmolkit.Molecule.from_smiles("CCO").fingerprint_morgan()
    propanol = cosmolkit.Molecule.from_smiles("CCCO").fingerprint_morgan()

    rd_ethanol = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048).GetFingerprint(
        Chem.MolFromSmiles("CCO")
    )
    rd_propanol = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048).GetFingerprint(
        Chem.MolFromSmiles("CCCO")
    )
    expected = DataStructs.TanimotoSimilarity(rd_ethanol, rd_propanol)
    assert ethanol.tanimoto(propanol) == pytest.approx(expected)


def test_morgan_fingerprint_advanced_generators_and_counts_match_rdkit():
    assert_morgan_bits_equal(
        "N[C@@H](C)C(=O)O",
        radius=2,
        n_bits=512,
        include_chirality=True,
        atom_invariants_generator="feature",
        num_bits_per_feature=2,
    )
    assert_morgan_bits_equal(
        "c1ccccc1O",
        radius=3,
        n_bits=256,
        atom_invariants_generator="connectivity",
        atom_invariants_include_ring_membership=False,
        bond_invariants_generator="morgan",
        bond_invariants_use_bond_types=False,
        count_simulation=True,
        count_bounds=[1, 2, 4, 8],
    )


def test_morgan_fingerprint_custom_invariants_and_root_atoms_match_rdkit():
    assert Chem is not None
    smiles = "CC(C)O"
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    assert_morgan_bits_equal(
        smiles,
        radius=2,
        n_bits=256,
        from_atoms=[0],
        custom_atom_invariants=[idx + 1 for idx in range(mol.GetNumAtoms())],
        custom_bond_invariants=[idx + 7 for idx in range(mol.GetNumBonds())],
        only_nonzero_invariants=True,
        include_redundant_environments=True,
    )


def test_morgan_additional_output_matches_rdkit():
    assert Chem is not None
    assert rdFingerprintGenerator is not None
    smiles = "CC(C)O"
    rdmol = Chem.MolFromSmiles(smiles, sanitize=True)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=256)
    additional = rdFingerprintGenerator.AdditionalOutput()
    additional.AllocateAtomCounts()
    additional.AllocateAtomToBits()
    additional.AllocateBitInfoMap()
    additional.AllocateAtomsPerBit()
    expected_fp = generator.GetFingerprint(rdmol, additionalOutput=additional)

    result = cosmolkit.Molecule.from_smiles(smiles).fingerprint_morgan_with_output(
        radius=2, n_bits=256
    )
    output = result.additional_output()

    assert result.fingerprint().on_bits() == list(expected_fp.GetOnBits())
    assert output.atom_counts() == list(additional.GetAtomCounts())
    assert output.atom_to_bits() == [list(bits) for bits in additional.GetAtomToBits()]
    assert output.bit_info_map() == {
        bit: list(entries) for bit, entries in additional.GetBitInfoMap().items()
    }
    assert output.atoms_per_bit() == {
        bit: [list(atom_ids) for atom_ids in entries]
        for bit, entries in additional.GetAtomsPerBit().items()
    }


def test_morgan_batch_api_preserves_order_and_none_for_invalid_records():
    batch = cosmolkit.MoleculeBatch.from_smiles_list(
        ["CCO", "not-a-smiles", "CCCO"], errors="keep"
    )
    values = batch.fingerprint_morgan_list(n_bits=256, n_jobs=2)
    assert [value is not None for value in values] == [True, False, True]
    assert values[0] is not None
    assert values[0].on_bits() == cosmolkit.Molecule.from_smiles("CCO").fingerprint_morgan(
        n_bits=256
    ).on_bits()

    with_output = batch.fingerprint_morgan_with_output_list(n_bits=256, n_jobs=2)
    assert with_output[0] is not None
    assert with_output[0].additional_output().atom_counts()
    assert with_output[1] is None


def test_morgan_binding_rejects_unknown_generators():
    mol = cosmolkit.Molecule.from_smiles("CCO")
    with pytest.raises(ValueError, match="atom_invariants_generator"):
        mol.fingerprint_morgan(atom_invariants_generator="unknown")
    with pytest.raises(ValueError, match="bond_invariants_generator"):
        mol.fingerprint_morgan(bond_invariants_generator="unknown")
