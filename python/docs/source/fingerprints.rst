Fingerprints
============

.. meta::
   :description: COSMolKit molecular fingerprint APIs for source-backed Morgan, MACCS, RDKit topological, Avalon, and AtomPair bit and count vectors with explicit parameters and RDKit parity boundaries.

COSMolKit exposes fixed-length bit vectors plus the source-defined sparse bit
and count forms used by AtomPair. The exposed Morgan and MACCS branches are
covered by strict RDKit bit-identical parity tests. The source-backed
topological, Avalon, and AtomPair implementations follow their pinned upstream
algorithms and exact maintained-corpus validation.
Similarity-shape correlation or structurally similar hashing is not a
compatibility claim. The Python ``Fingerprint`` object is a sparse view over
the binary vector: ``on_bits()`` returns the bit indexes whose value is 1. It
is not a dense floating-point neural embedding.

Source-backed topological and Avalon fingerprints
-------------------------------------------------

``Molecule.topological_fingerprint()`` and
``Molecule.avalon_fingerprint()`` execute source-backed Rust implementations.
The maintained topological matrix is exact across 5,000 rows and 14 profiles;
the Avalon matrix is exact across 5,000 rows and 23 profiles.

The two APIs are source-backed and return a fresh explicit bit vector:

- ``topological_fingerprint()`` requires the complete RDKit
  ``RDKFingerprintMol``/RDKitFP generator behavior, including branched-path
  enumeration, source random-bit generation, density folding, atom
  invariants, and the exposed path and atom-selection parameters.
- ``avalon_fingerprint()`` follows the complete Avalon/reaccs bit-vector path,
  including ``bitFlags``, ``isQuery``, hydrogen handling, aromaticity passes,
  and byte-rounded vector semantics. ``resetVect`` is an internal adapter
  detail and is not exposed on COSMolKit's value-returning API.

The focused option fixtures and maintained corpus are exact source
comparisons. Similarity correlation, partial bit agreement, or a heuristic
replacement is not an acceptance condition.

.. code-block:: python

   mol = Molecule.from_smiles("c1ccccc1O")

   topological = mol.topological_fingerprint(
       min_path=1,
       max_path=7,
       fp_size=2048,
       num_bits_per_feature=2,
   )
   avalon = mol.avalon_fingerprint(
       n_bits=512,
       is_query=False,
       bit_flags=0xF07FFF,
   )

   print(topological.on_bits())
   print(avalon.on_bits())

AtomPair fingerprints
---------------------

``Molecule`` exposes the four source generator result forms through one
project-native implementation:

- ``fingerprint_atom_pair_sparse_count()`` returns the source-width count map;
- ``fingerprint_atom_pair_count()`` returns the folded count map;
- ``fingerprint_atom_pair_sparse_bits()`` returns the source-width bit set;
- ``fingerprint_atom_pair()`` returns the fixed-width explicit bit vector.

The shared parameter surface covers 2D topological or 3D conformer distances,
minimum and maximum distance, chirality, count simulation and bounds,
``num_bits_per_feature``, ``from_atoms``, ``ignore_atoms``, conformer selection,
and custom atom invariants. A 3D call requires a conformer and sets
``use_2d=False`` explicitly.

.. code-block:: python

   mol = Molecule.from_smiles("CCCO")

   explicit = mol.fingerprint_atom_pair(n_bits=2048)
   sparse_count = mol.fingerprint_atom_pair_sparse_count()
   folded_count = mol.fingerprint_atom_pair_count(n_bits=2048)
   sparse_bits = mol.fingerprint_atom_pair_sparse_bits()

   print(explicit.on_bits())
   print(sparse_count.nonzero_elements())
   print(folded_count.nonzero_elements())
   print(sparse_bits.on_bits())

``fingerprint_atom_pair_with_output()`` additionally returns exact atom counts,
atom-to-bit rows, endpoint-pair bit information, and atoms-per-bit provenance:

.. code-block:: python

   result = mol.fingerprint_atom_pair_with_output()
   output = result.additional_output()

   print(result.fingerprint().on_bits())
   print(output.atom_counts())
   print(output.atom_to_bits())
   print(output.bit_info_map())
   print(output.atoms_per_bit())

The complete ChEMBL 37 validation compares 40 vectors and one full provenance
output for every one of the 2,897,804 mutually parseable records: 118,809,964
exact comparisons with zero mismatch. The maintained 5,000-row, ten-profile
matrix remains the committed continuous regression gate.

Topological provenance
-----------------------

``topological_fingerprint_with_output()`` returns a
``TopologicalFingerprintResult``. Request ``atom_bits`` and/or ``bit_info`` to
receive the matching source provenance outputs:

.. code-block:: python

   result = mol.topological_fingerprint_with_output(
       fp_size=2048,
       atom_bits=True,
       bit_info=True,
   )

   print(result.fingerprint().on_bits())
   print(result.atom_bits())
   print(result.bit_info())

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
the fingerprint and experimental provenance data:

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

The Python binding exposes the source-backed Morgan generator branches covered
by exact RDKit bit parity:

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

The list is the tested boundary, not a claim that every RDKit fingerprint
generator or every input-state preparation branch is implemented. Unsupported
branches propagate an error; they do not silently fall back to another
fingerprint algorithm.

Batch Fingerprints
------------------

``MoleculeBatch`` exposes matching batch APIs. Invalid records kept with
``errors="keep"`` produce ``None`` in the corresponding output position.

.. code-block:: python

   from cosmolkit import MoleculeBatch

   batch = MoleculeBatch.from_smiles_list(
       ["CCO", "not-smiles", "CCCO"],
       errors="keep",
   ).with_parallel_jobs(8)
   fps = batch.fingerprint_morgan_list(n_bits=2048)
   atom_pairs = batch.fingerprint_atom_pair_list(
       n_bits=2048,
       n_jobs=8,
       progress_bar=False,
   )

   print([fp.on_bits() if fp is not None else None for fp in fps])
   print([fp.on_bits() if fp is not None else None for fp in atom_pairs])
