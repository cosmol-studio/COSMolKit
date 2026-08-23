Batch Workflows
===============

.. meta::
   :description: Process ordered molecule collections with COSMolKit batch APIs, structured per-record errors, deterministic parallel jobs, transforms, descriptors, and exports.

``MoleculeBatch`` is an ordered collection for processing many molecules with a
single API call. Valid records keep their original input order through
transform, export, and filtering steps.

.. code-block:: python

   from cosmolkit import BatchErrorMode, BatchValidationError, MoleculeBatch

   batch = MoleculeBatch.from_smiles_list(
       ["CCO", "c1ccccc1", "not-smiles"],
       errors=BatchErrorMode.KEEP,
   ).with_parallel_jobs(8)

   prepared = batch.with_hydrogens(errors=BatchErrorMode.KEEP).with_2d_coordinates(
       errors=BatchErrorMode.KEEP,
   )

   print(prepared.valid_mask())
   print(prepared.errors())

Error Handling
--------------

Batch APIs accept ``errors``:

- ``"raise"`` raises an exception when a record fails.
- ``"keep"`` keeps failed records and exposes structured errors. Export methods
  write valid records and count invalid records as skipped in the returned
  report.

String modes remain supported, but Python callers can also pass
``BatchErrorMode`` enum members. Per-record ``BatchError`` values expose the
input index, operation name, and message:

.. code-block:: python

   for error in batch.errors():
       print(error.index(), error.operation(), error.message())

   try:
       MoleculeBatch.from_smiles_list(["C1CC"], errors=BatchErrorMode.RAISE)
   except BatchValidationError as exc:
       print(exc.error_count)

The read-only ``BATCH_ERROR_MODE_MAP`` converts external string names to enum
members when needed.

Batch Values
------------

``MoleculeBatch`` behaves like an ordered Python container. Valid records are
returned as ``Molecule`` objects and invalid kept records are returned as
``None``:

.. code-block:: python

   molecules = prepared.to_list()
   first = prepared[0]
   tail = prepared[5:]
   valid = prepared[prepared.valid_mask()]

   for molecule in prepared:
       if molecule is not None:
           print(molecule.to_smiles())

Integer indexing returns ``Molecule | None`` because kept invalid records are
represented as ``None``. Slices, integer index lists, and boolean masks return a
new ``MoleculeBatch`` and preserve both input order and the batch-level
parallel-job setting.

Export Images
-------------

.. code-block:: python

   report = prepared.to_images(
       "molecule_images",
       format="png",
       size=(300, 300),
       errors="keep",
       filenames=["ethanol.png", "benzene.png", "invalid.png"],
       report_path="image_errors.json",
   )

   print(report.total(), report.success(), report.failed())

``filenames`` is optional. Entries must match the batch length; ``None`` uses
the default zero-padded name for that record. Names are relative to the output
directory, and missing extensions are filled from ``format``.

Export SDF
----------

.. code-block:: python

   report = prepared.to_sdf(
       "prepared.sdf",
       format="v2000",
       errors="keep",
       report_path="sdf_errors.csv",
   )

Use ``to_sdf_files()`` when each valid record should be written to its own SDF
file:

.. code-block:: python

   report = prepared.to_sdf_files(
       "prepared_records",
       format="v2000",
       errors="keep",
       filenames=["ethanol", "benzene.sdf", "invalid.sdf"],
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
   fingerprints = prepared.fingerprint_morgan_list(n_bits=2048)
   atom_pairs = prepared.fingerprint_atom_pair_list(n_bits=2048)

Morgan and AtomPair batch APIs use the same source-backed generator cores as
the single-molecule APIs, preserve record order, and are covered by exact
pinned-RDKit parity tests. Topological Torsion exposes the same ordered bulk
execution through its generator methods. Morgan and AtomPair can also collect
provenance through the following ``MoleculeBatch`` conveniences:

.. code-block:: python

   results = prepared.fingerprint_morgan_with_output_list(
       radius=2,
       n_bits=2048,
   )

   for result in results:
       if result is not None:
           print(result.fingerprint().on_bits())
           print(result.additional_output().bit_info_map())

   atom_pair_results = prepared.fingerprint_atom_pair_with_output_list(
       n_bits=2048,
   )

   for result in atom_pair_results:
       if result is not None:
           print(result.fingerprint().on_bits())
           print(result.additional_output().bit_info_map())

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

``with_parallel_jobs()`` returns a new batch with a default worker count for
later parallel operations. Because molecule values use copy-on-write storage,
this configuration step does not duplicate the molecular data.

.. code-block:: python

   configured = batch.with_parallel_jobs(8)
   prepared = configured.with_2d_coordinates(errors="keep")
   smiles = prepared.to_smiles_list()

Method-level ``n_jobs`` still overrides the batch default for a single call:

.. code-block:: python

   svgs = prepared.to_svg_list(n_jobs=2)

``with_progress_bar()`` returns a new batch with a default Rust-side progress
bar setting. Progress is emitted by Rust to stderr, matching the usual terminal
stream for progress indicators, and method-level ``progress_bar`` overrides the
batch default for one call:

.. code-block:: python

   tracked = batch.with_parallel_jobs(8).with_progress_bar(True)
   prepared = tracked.with_2d_coordinates(errors="keep")
   smiles = prepared.to_smiles_list(progress_bar=False)

Pass ``None`` to clear the batch-level default and let rayon choose:

.. code-block:: python

   default_scheduled = prepared.with_parallel_jobs(None)
   quiet = prepared.with_progress_bar(None)
