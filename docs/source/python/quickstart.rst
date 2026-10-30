Quick Start
===========

.. meta::
   :description: Install COSMolKit, parse molecules, draw structures, calculate fingerprints and descriptors, and process batches.

Installation
------------

.. code-block:: bash

   pip install cosmolkit

Check the installed version and use its matching documentation:

.. code-block:: python

   import cosmolkit as ck

   print(ck.__version__)

Create and inspect a molecule
-----------------------------

.. code-block:: python

   mol = ck.mol_from_smiles("CCO")
   print(mol.to_smiles())       # CCO
   print(mol.num_atoms)         # 3
   for atom in mol.atoms:
       print(atom.id, atom.element, atom.degree)

Top-level constructors such as ``ck.mol_from_smiles`` accept the same options
as their class factories, such as ``ck.Molecule.from_smiles``.

Values and configuration
------------------------

Transforms return new molecules. Methods ending in ``_`` modify their receiver.
``copy()`` shares molecular storage until a modification needs a copy:

.. code-block:: python

   hydrogenated = mol.with_hydrogens()
   assert mol.num_atoms == 3
   assert hydrogenated.num_atoms == 9

   editable = mol.copy()
   editable.add_hydrogens_()
   assert mol.num_atoms == 3

Pass either configuration keywords or a reusable parameter object, not both:

.. code-block:: python

   fp = mol.fingerprint_morgan(radius=2, fp_size=2048)
   options = ck.MorganFingerprintParams(radius=2, fp_size=2048)
   options.from_atoms = [0, 1]
   selected_fp = mol.fingerprint_morgan(options)

Assign sequence fields as a whole; their getters are immutable tuples.
For example, use ``options.from_atoms = [0, 1]``, not
``options.from_atoms.append(1)``. See :doc:`fingerprints` for generators and
additional output.

Drawing, coordinates and files
------------------------------

.. code-block:: python

   drawn = mol.with_2d_coordinates()
   print(drawn.coordinates_2d.shape)  # (3, 2)
   svg = drawn.to_svg(width=400, height=300)
   drawn.write_svg("ethanol.svg", width=400, height=300)
   drawn.write_sdf("ethanol.sdf")
   restored = ck.Molecule.read_sdf("ethanol.sdf")
   print(restored.to_smiles())

Use :doc:`io` for streaming SDF, MOL2, coordinate input and file options.
Paths are literal: expand ``~`` yourself with ``pathlib.Path.expanduser()``.

Fingerprints, search and descriptors
------------------------------------

.. code-block:: python

   phenol = ck.mol_from_smiles("c1ccccc1O")
   benzene = ck.mol_from_smiles("c1ccccc1")
   left = phenol.fingerprint_morgan()
   right = benzene.fingerprint_morgan()

   print(left.tanimoto(right))
   print(left.to_numpy().shape)  # Dense uint8 bit vector
   print(left.on_bits)           # Set-bit indexes
   print(phenol.molecular_weight(), phenol.tpsa())

   query = ck.parse_smarts("[#6]-O")
   print(phenol.has_substruct_match(query))

See :doc:`molecule` for stereochemistry and search callbacks,
:doc:`fingerprints` for fingerprint families, and :doc:`descriptors` for
molecular properties.

Conformers and force fields
---------------------------

Add explicit hydrogens before generating all-atom geometry:

.. code-block:: python

   spatial = mol.with_hydrogens()
   params = ck.EmbedParams.etkdg_v3()
   params.random_seed = 42
   params.track_failures = True
   result = spatial.with_3d_conformer_result(params)
   print(result.conf_id, result.failure_counts)

   if result.ok:
       optimized = result.molecule.with_uff_optimized()
       print(optimized.energy, optimized.needs_more)

``needs_more`` means the optimizer did not converge within its iteration
limit. See :doc:`forcefields` for MMFF, multiple conformers and constraints.

Batch processing
----------------

.. code-block:: python

   batch = ck.mols_from_smiles_list(["CCO", "CCN", "invalid"], errors="keep")
   print(batch.valid_mask)
   print(batch.to_smiles())
   for error in batch.errors:
       print(error.index, error.message)

   valid = batch[batch.valid_mask]
   matrix = valid.fingerprint_morgan(radius=2, n_jobs=4).to_numpy()
   print(matrix.shape)  # (2, 2048)

Batches preserve input order. ``errors="keep"`` retains failed positions;
``errors="raise"`` raises instead. See :doc:`batch` for scheduling and exports.

Next steps
----------

- :doc:`reaction`: reaction templates and products.
- :doc:`protein`: PDB/mmCIF and model → chain → residue → atom traversal.
- :doc:`mcs`: maximum common substructure.
- :doc:`api`: complete signatures, parameter types and return values.

Reference versions, comparison rules and completed results are recorded in
`VALIDATION.md <https://github.com/cosmol-studio/COSMolKit/blob/main/VALIDATION.md>`__.
