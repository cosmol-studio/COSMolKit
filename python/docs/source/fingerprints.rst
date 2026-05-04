Fingerprints
============

COSMolKit exposes RDKit-style Morgan fingerprints as fixed-length bit vectors.
The Python ``Fingerprint`` object is a sparse view over that binary vector:
``on_bits()`` returns the bit indexes whose value is 1. It is not a dense
floating-point neural embedding.

Single Molecules
----------------

.. code-block:: python

   from cosmolkit import Molecule

   mol = Molecule.from_smiles("c1ccccc1O")
   fp = mol.fingerprint_morgan(radius=2, n_bits=2048)

   print(fp.n_bits())
   print(fp.on_bits())

Tanimoto similarity is computed directly on ``Fingerprint`` values:

.. code-block:: python

   phenol = Molecule.from_smiles("c1ccccc1O").fingerprint_morgan()
   benzene = Molecule.from_smiles("c1ccccc1").fingerprint_morgan()

   print(phenol.tanimoto(benzene))

Additional Output
-----------------

``fingerprint_morgan_with_output()`` returns a ``MorganFingerprintResult`` with
the fingerprint and RDKit-style provenance data:

.. code-block:: python

   result = mol.fingerprint_morgan_with_output(radius=2, n_bits=2048)
   output = result.additional_output()

   print(result.fingerprint().on_bits())
   print(output.atom_counts())
   print(output.atom_to_bits())
   print(output.bit_info_map())
   print(output.atoms_per_bit())

Supported Parameters
--------------------

The Python binding exposes the supported RDKit-style Morgan generator branches:

- ``radius`` and ``n_bits``
- ``include_chirality`` and ``use_bond_types``
- ``count_simulation`` and ``count_bounds``
- ``only_nonzero_invariants``
- ``include_redundant_environments``
- ``from_atoms`` and ``ignore_atoms``
- ``custom_atom_invariants`` and ``custom_bond_invariants``
- ``atom_invariants_generator="connectivity" | "morgan" | "feature" | "fcfp"``
- ``atom_invariants_include_ring_membership``
- ``bond_invariants_generator="morgan" | "default" | "bond"``
- ``bond_invariants_use_bond_types``
- ``bond_invariants_use_chirality``
- ``num_bits_per_feature``

Batch Fingerprints
------------------

``MoleculeBatch`` exposes matching batch APIs. Invalid records kept with
``errors="keep"`` produce ``None`` in the corresponding output position.

.. code-block:: python

   from cosmolkit import MoleculeBatch

   batch = MoleculeBatch.from_smiles_list(
       ["CCO", "not-smiles", "CCCO"],
       errors="keep",
   )
   fps = batch.fingerprint_morgan_list(n_bits=2048, n_jobs=8)

   print([fp.on_bits() if fp is not None else None for fp in fps])
