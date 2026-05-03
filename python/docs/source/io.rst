File IO and Arrays
==================

SDF Files
---------

Read the first molecule from an SDF file:

.. code-block:: python

   from cosmolkit import Molecule

   mol = Molecule.read_sdf("input.sdf", coordinate_dim="auto")

Write a molecule to SDF:

.. code-block:: python

   mol = Molecule.from_smiles("CCO").with_2d_coords()
   mol.write_sdf("ethanol.sdf", format="v2000")

SDF Strings
-----------

.. code-block:: python

   text = mol.to_sdf_string(format="v2000")
   restored = Molecule.read_sdf_record_from_str(text, coordinate_dim="2d")

For multi-record strings, use ``read_sdf_records_from_str()``:

.. code-block:: python

   molecules = Molecule.read_sdf_records_from_str(sdf_text, coordinate_dim="auto")

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

   mol.write_svg("phenol.svg", width=400, height=300)
   mol.write_png("phenol.png", width=400, height=300)
