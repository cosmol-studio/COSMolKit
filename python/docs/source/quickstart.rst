Quick Start
===========

Copy-On-Write (COW) Molecule Values
-----------------------------------

COSMolKit molecules use copy-on-write (COW) value semantics. Transform methods
return a new ``Molecule`` and leave the original object unchanged:

.. code-block:: python

   from cosmolkit import Molecule

   mol = Molecule.from_smiles("CCO")
   mol_h = mol.with_hydrogens()

   assert mol is not mol_h
   print(mol.to_smiles())
   print(mol_h.to_smiles())

This is an intentional difference from common RDKit Python usage. Do not assume
that a transform mutates the existing object; always keep the returned
``Molecule``.

Create a molecule from SMILES and export a depiction:

.. code-block:: python

   from cosmolkit import Molecule

   mol = Molecule.from_smiles("c1ccccc1O")
   drawn = mol.with_2d_coords()

   print(mol.to_smiles())
   drawn.write_png("phenol.png", width=400, height=300)

Inspect atoms and bonds:

.. code-block:: python

   for atom in mol.atoms():
       print(atom.idx(), atom.atomic_num(), atom.is_aromatic())

   for bond in mol.bonds():
       print(bond.begin_atom_idx(), bond.end_atom_idx(), bond.bond_type())

Inspect chiral tags without converting to an ordered tetrahedral record:

.. code-block:: python

   chiral = Molecule.from_smiles("F[C@H](Cl)Br")

   print(chiral.to_smiles())
   print(chiral.to_smiles(isomeric_smiles=False))

   for atom in chiral.atoms():
       if atom.chiral_tag() != "CHI_UNSPECIFIED":
           print(atom.idx(), atom.chiral_tag())

Read and write SDF:

.. code-block:: python

   mol = Molecule.read_sdf("input.sdf", coordinate_dim="auto")
   mol.write_sdf("output.sdf", format="v2000")

Access coordinates as NumPy arrays:

.. code-block:: python

   mol2d = Molecule.from_smiles("CCO").with_2d_coords()
   coords = mol2d.coords_2d()

   print(coords.shape)

Process a list of molecules:

.. code-block:: python

   from cosmolkit import MoleculeBatch

   batch = MoleculeBatch.from_smiles_list(
       ["CCO", "c1ccccc1", "not-smiles"],
       errors="keep",
   )

   prepared = batch.compute_2d_coords(errors="keep")

   print(prepared.valid_mask())
   print(prepared.to_smiles_list(canonical=True))
