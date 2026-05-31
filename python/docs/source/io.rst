File IO and Arrays
==================

SDF Files
---------

COSMolKit exposes multiple SDF reading styles because small files, large
seekable files, and stream inputs have different memory and access patterns.

Use ``Molecule.read_sdf()`` when you only need the first record:

.. code-block:: python

   from cosmolkit import Molecule

   mol = Molecule.read_sdf("input.sdf", coordinate_dim="auto")

Use ``MoleculeBatch.read_sdf()`` when you intentionally want the entire file as
one in-memory batch:

.. code-block:: python

   from cosmolkit import MoleculeBatch

   batch = MoleculeBatch.read_sdf(
       "input.sdf",
       coordinate_dim="auto",
       errors="keep",
       progress_bar=True,
   )

``progress_bar=True`` first builds a lightweight record index so the progress
bar has an accurate total. For very large files this is still an all-in-memory
batch; use ``SdfDataset.batches()`` when bounded memory matters.

Use ``SdfDataset`` for large seekable files when you want cheap ``len()``,
random access, metadata inspection, or chunked processing:

.. code-block:: python

   from cosmolkit import SdfDataset

   dataset = SdfDataset.open("large.sdf", coordinate_dim="auto")

   print(len(dataset))
   print(dataset.metadata(0).title())

   record = dataset[12345]
   print(record.title(), record.data_field("supplier_id"))

   head = dataset[:1024]
   print(head.valid_mask())

   for batch in dataset.batches(size=1024, errors="keep", progress_bar=True):
       fingerprints = batch.fingerprint_morgan_list()

``dataset[i]`` returns one ``SdfRecord`` with SDF metadata and data fields.
Slices, integer index lists, and boolean masks return ``MoleculeBatch`` objects.

Use ``SdfReader`` for one-pass stream-style processing where random access is
not needed:

.. code-block:: python

   from cosmolkit import SdfReader

   for batch in SdfReader.open("large.sdf").batches(size=1024, errors="keep"):
       smiles = batch.to_smiles_list(canonical=True)

``SdfReader`` does not know the final record count without pre-indexing the
file, so accurate record-count progress belongs to ``SdfDataset``.

Read a single-record MDL molfile with the same CTAB parser:

.. code-block:: python

   mol = Molecule.read_mol("input.mol", coordinate_dim="auto")
   mol = Molecule.read_mol_from_str(mol_text, coordinate_dim="2d")

``Molecule.read_mol()`` and ``Molecule.read_mol_from_str()`` follow RDKit
``MolFromMolBlock`` boundaries: they parse the molfile CTAB through the first
``M  END`` line and ignore unread trailing text, including SDF data fields and
``$$$$`` record separators. Use ``Molecule.read_sdf()``, ``SdfDataset``, or
``SdfReader`` when those SDF data fields are part of the requested input.

Write a molecule to SDF. SDF writing is explicit about the coordinate source:
``write_sdf()`` and ``to_2d_sdf_string()`` export 2D coordinates, generating
them when needed; ``to_3d_sdf_string()`` exports an existing 3D conformer and
raises if the molecule has no 3D coordinates.

.. code-block:: python

   mol = Molecule.from_smiles("CCO").with_2d_coordinates()
   mol.write_sdf(
       "python/examples/output/ethanol.sdf",
       format="v2000",
       include_stereo=True,
       kekulize=True,
   )

SDF Strings
-----------

.. code-block:: python

   text = mol.to_2d_sdf_string(format="v2000", include_stereo=True, kekulize=True)
   restored = Molecule.read_sdf_from_str(text, coordinate_dim="2d")

``Molecule.read_sdf_from_str()`` uses the SDF record parser and therefore
validates and parses data fields after ``M  END``. For molfile-only string input
where trailing SDF text should be ignored, use ``Molecule.read_mol_from_str()``.

Use ``to_3d_sdf_string()`` when you explicitly want to export an existing 3D
conformer. The molecule must already have 3D coordinates.

The ``format`` argument accepts ``"v2000"``, ``"v3000"``, or ``None`` for
automatic selection. ``include_stereo=False`` and ``kekulize=False`` are exposed
as RDKit parity parameters; branches that are not implemented yet fail with a
clear unsupported-path error instead of silently changing behavior.

For multi-record strings, use the batch API:

.. code-block:: python

   batch = MoleculeBatch.read_sdf_records_from_str(sdf_text, coordinate_dim="auto")

MOL2 Files
----------

Read Tripos MOL2 input with the source-ported RDKit ``Mol2FileToMol`` and
``Mol2BlockToMol`` profile:

.. code-block:: python

   mol = Molecule.read_mol2("input.mol2")
   mol = Molecule.read_mol2_from_str(
       mol2_text,
       sanitize=True,
       remove_hs=True,
       variant="corina",
       cleanup_substructures=True,
   )

The MOL2 reader exposes RDKit's ``Mol2ParserParams`` controls. ``variant`` is
currently limited to ``"corina"``, the only variant present in RDKit's public
MOL2 enum for this parser.

PDB and mmCIF Blocks
--------------------

Use ``Molecule.from_pdb_block()`` when you want a molecule state comparable to
RDKit ``Chem.MolFromPDBBlock`` for the modeled conversion profile:

.. code-block:: python

   pdb_mol = Molecule.from_pdb_block(
       pdb_text,
       sanitize=True,
       remove_hs=True,
       proximity_bonding=True,
   )

Use ``Molecule.from_mmcif_block()`` for mmCIF structural text. COSMolKit reads
the mmCIF structure and then applies the same molecule-conversion profile used
by ``from_pdb_block()``:

.. code-block:: python

   cif_mol = Molecule.from_mmcif_block(
       cif_text,
       sanitize=True,
       remove_hs=True,
       proximity_bonding=True,
   )

XYZ Blocks
----------

Use ``Molecule.from_xyz_block()`` when you want to read atom identities and
Cartesian coordinates from XYZ text:

.. code-block:: python

   xyz_mol = Molecule.from_xyz_block(xyz_text)

XYZ input contains coordinates but no bond table, so the returned molecule has
atoms and one 3D conformer without inferred bonds. This matches COSMolKit's
RDKit-aligned ``MolFromXYZBlock`` profile.

Coordinate Arrays
-----------------

``coords_2d()``, ``coords_3d()``, and ``dg_bounds_matrix()`` return NumPy
arrays. ``coords_2d()`` has shape ``(num_atoms, 3)`` with a zero-filled z
column; ``coords_3d()`` has shape ``(num_atoms, 3)`` for the selected 3D
conformer:

.. code-block:: python

   mol = Molecule.from_smiles("c1ccccc1O").with_2d_coordinates()

   coords = mol.coords_2d()
   bounds = mol.dg_bounds_matrix()

   print(coords.shape)
   print(bounds.shape)

Depiction Files
---------------

.. code-block:: python

   mol = Molecule.from_smiles("c1ccccc1O").with_2d_coordinates()

   mol.write_svg("python/examples/output/phenol.svg", width=400, height=300)
   mol.write_png("python/examples/output/phenol.png", width=400, height=300)
