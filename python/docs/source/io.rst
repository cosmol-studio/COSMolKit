File IO and Arrays
==================

SDF Files
---------

Read the first molecule from an SDF file:

.. code-block:: python

   from cosmolkit import Molecule, MoleculeBatch

   mol = Molecule.read_sdf("input.sdf", coordinate_dim="auto")

Read a single-record MDL molfile with the same CTAB parser:

.. code-block:: python

   mol = Molecule.read_mol("input.mol", coordinate_dim="auto")
   mol = Molecule.read_mol_from_str(mol_text, coordinate_dim="2d")

Write a molecule to SDF. SDF writing is explicit about the coordinate source:
``write_sdf()`` and ``to_2d_sdf_string()`` export 2D coordinates, generating
them when needed; ``to_3d_sdf_string()`` exports an existing 3D conformer and
raises if the molecule has no 3D coordinates.

.. code-block:: python

   mol = Molecule.from_smiles("CCO").with_2d_coords()
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

Use ``to_3d_sdf_string()`` when you explicitly want to export an existing 3D
conformer. The molecule must already have 3D coordinates.

The ``format`` argument accepts ``"v2000"``, ``"v3000"``, or ``None`` for
automatic selection. ``include_stereo=False`` and ``kekulize=False`` are exposed
as RDKit parity parameters; branches that are not implemented yet fail with a
clear unsupported-path error instead of silently changing behavior.

For multi-record strings, use the batch API:

.. code-block:: python

   batch = MoleculeBatch.read_sdf_records_from_str(sdf_text, coordinate_dim="auto")

Coordinate Arrays
-----------------

``coords_2d()``, ``coords_3d()``, and ``dg_bounds_matrix()`` return NumPy
arrays:

.. code-block:: python

   mol = Molecule.from_smiles("c1ccccc1O").with_2d_coords()

   coords = mol.coords_2d()
   bounds = mol.dg_bounds_matrix()

   print(coords.shape)
   print(bounds.shape)

Depiction Files
---------------

.. code-block:: python

   mol = Molecule.from_smiles("c1ccccc1O").with_2d_coords()

   mol.write_svg("python/examples/output/phenol.svg", width=400, height=300)
   mol.write_png("python/examples/output/phenol.png", width=400, height=300)
