Batch Workflows
===============

``MoleculeBatch`` is an ordered collection for processing many molecules with a
single API call.

.. code-block:: python

   from cosmolkit import MoleculeBatch

   batch = MoleculeBatch.from_smiles_list(
       ["CCO", "c1ccccc1", "not-smiles"],
       errors="keep",
   )

   prepared = batch.add_hydrogens(errors="keep").compute_2d_coords(errors="keep")

   print(prepared.valid_mask())
   print(prepared.errors())

Error Handling
--------------

Batch APIs accept ``errors``:

- ``"raise"`` raises an exception when a record fails.
- ``"keep"`` keeps failed records and exposes structured errors.
- ``"skip"`` omits failed records from the returned result or export.

Export Images
-------------

.. code-block:: python

   report = prepared.to_images(
       "molecule_images",
       format="png",
       size=(300, 300),
       errors="skip",
       report_path="image_errors.json",
   )

   print(report.total(), report.success(), report.failed())

Export SDF
----------

.. code-block:: python

   report = prepared.to_sdf(
       "prepared.sdf",
       format="v2000",
       errors="skip",
       report_path="sdf_errors.csv",
   )

Derived Outputs
---------------

.. code-block:: python

   smiles = prepared.to_smiles_list()
   rooted = prepared.to_smiles_list(rooted_at_atom=0)
   explicit = prepared.to_smiles_list(
       all_bonds_explicit=True,
       all_hs_explicit=True,
   )
   svgs = prepared.to_svg_list(width=300, height=300)
   bounds = prepared.dg_bounds_matrix_list()

SMILES Options
--------------

``to_smiles_list()`` accepts the same output-shaping options for every record:

- ``isomeric_smiles`` includes stereochemical and isotopic information.
- ``canonical`` returns canonical SMILES when enabled.
- ``kekule`` writes aromatic systems in Kekule form.
- ``clean_stereo`` normalizes stereo output where possible.
- ``all_bonds_explicit`` writes explicit bond symbols.
- ``all_hs_explicit`` writes explicit hydrogens.
- ``include_dative_bonds`` includes dative bond notation.
- ``ignore_atom_map_numbers`` omits atom map numbers from canonical decisions.
- ``rooted_at_atom`` starts traversal from a selected atom index.

Batch Chirality
---------------

Batch SMILES output preserves isomeric chirality by default:

.. code-block:: python

   chiral_batch = MoleculeBatch.from_smiles_list(
       ["F[C@H](Cl)Br", "F[C@@H](Cl)Br"],
       errors="raise",
   )

   print(chiral_batch.to_smiles_list(isomeric_smiles=True))
   print(chiral_batch.to_smiles_list(isomeric_smiles=False))

Use ``canonical=False`` when you want output to stay closer to each record's
input traversal while keeping the same CW/CCW chiral tag path:

.. code-block:: python

   print(chiral_batch.to_smiles_list(canonical=False))

Parallel Work
-------------

Batch methods that accept ``n_jobs`` can run across multiple worker threads:

.. code-block:: python

   prepared = batch.compute_2d_coords(errors="keep", n_jobs=8)
   smiles = prepared.to_smiles_list(n_jobs=8)
