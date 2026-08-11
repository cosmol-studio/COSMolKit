Fingerprints
============

COSMolKit exposes fingerprint APIs as fixed-length bit vectors. The exposed
Morgan and MACCS branches are covered by strict RDKit bit-identical parity
tests. Similarity-shape correlation or structurally similar hashing is not a
compatibility claim. ``RDKFingerprintMol``/topological and Avalon fingerprints
are not implemented: those calls raise an unsupported-feature error instead of
returning approximate path-hash bits. The Python ``Fingerprint`` object is a
sparse view over the binary vector: ``on_bits()``
returns the bit indexes whose value is 1. It is not a dense floating-point
neural embedding.

Unsupported APIs and Planned Ports
----------------------------------

.. warning::

   ``Molecule.topological_fingerprint()`` and
   ``Molecule.avalon_fingerprint()`` are unsupported in COSMolKit 0.2.11.
   Calling either method raises ``ValueError`` with an ``unsupported
   fingerprint option`` message. They do not return an empty vector, reuse a
   Morgan fingerprint, or run a locally invented path-hash approximation.

The two APIs remain present so callers receive an explicit, stable failure at
the intended public boundary while their source-exact implementations are
pending:

- ``topological_fingerprint()`` requires the complete RDKit
  ``RDKFingerprintMol``/RDKitFP generator behavior, including branched-path
  enumeration, source random-bit generation, density folding, atom
  invariants, and the exposed path and atom-selection parameters.
- ``avalon_fingerprint()`` requires the complete Avalon/reaccs behavior,
  including ``bitFlags``, ``isQuery``, ``resetVect``, hydrogen handling,
  tautomeric fingerprints, and byte-rounded vector semantics.

Neither API will be marked supported until its implementation is copied from
the pinned upstream source under the project's two-axis source-marker rules and
produces bit-identical results on focused option fixtures and the maintained
5000-molecule parity corpus. Similarity correlation, partial bit agreement, or
a heuristic replacement is not an acceptance condition.

The executable follow-up work is tracked in the `topological and Avalon
fingerprint source-port plan
<https://github.com/cosmol-studio/COSMolKit/blob/main/dev/plans/rdkit_topological_avalon_fingerprint_port_plan.md>`_.

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

   print([fp.on_bits() if fp is not None else None for fp in fps])
