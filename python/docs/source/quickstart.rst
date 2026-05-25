Quick Start
===========

Copy-On-Write (COW) Molecule Values
-----------------------------------

COSMolKit molecules use copy-on-write (COW) value semantics. Transform methods
return a new ``Molecule`` and leave the original object unchanged:

.. code-block:: python

   from cosmolkit import BatchErrorMode, BondOrder, ChiralTag, Molecule

   mol = Molecule.from_smiles("CCO")
   mol_h = mol.with_hydrogens()

   assert mol is not mol_h
   print(mol.to_smiles())
   print(mol_h.to_smiles())

This is an intentional difference from common RDKit Python usage. Do not assume
that a transform mutates the existing object; always keep the returned
``Molecule``.

In-place molecule operations are explicit and always end with ``_``:

.. code-block:: python

   mol = Molecule.from_smiles("CCO")
   mol.add_hydrogens_()
   mol.compute_2d_coords_()

Create a molecule from SMILES and export a depiction:

.. code-block:: python

   from cosmolkit import Molecule

   mol = Molecule.from_smiles("c1ccccc1O")
   drawn = mol.with_2d_coords()

   print(mol.to_smiles())
   drawn.write_png("python/examples/output/phenol.png", width=400, height=300)

Inspect atoms and bonds:

.. code-block:: python

   from cosmolkit import BondOrder, Molecule

   mol = Molecule.from_smiles("c1ccccc1O")

   for atom in mol.atoms():
       print(atom.idx(), atom.atomic_num(), atom.is_aromatic())

   for bond in mol.bonds():
       if bond.bond_type() == BondOrder.SINGLE:
           print(bond.begin_atom_idx(), bond.end_atom_idx(), bond.bond_type().name)

Inspect chiral tags without converting to an ordered tetrahedral record:

.. code-block:: python

   chiral = Molecule.from_smiles("F[C@H](Cl)Br")

   print(chiral.to_smiles())
   print(chiral.to_smiles(isomeric_smiles=False))

   for atom in chiral.atoms():
       if atom.chiral_tag() != ChiralTag.CHI_UNSPECIFIED:
           print(atom.idx(), atom.chiral_tag().name)

Read and write SDF:

.. code-block:: python

   mol = Molecule.read_sdf("input.sdf", coordinate_dim="auto")
   mol.write_sdf("python/examples/output/output.sdf", format="v2000")

Access coordinates as NumPy arrays:

.. code-block:: python

   mol2d = Molecule.from_smiles("CCO").with_2d_coords()
   coords = mol2d.coords_2d()

   print(coords.shape)

Generate a Morgan fingerprint:

.. code-block:: python

   fp = Molecule.from_smiles("c1ccccc1O").fingerprint_morgan(
       radius=2,
       n_bits=2048,
   )

   print(fp.n_bits())
   print(fp.on_bits())

``on_bits()`` returns the sparse bit indexes set inside the fixed-length binary
fingerprint. It is not a dense neural embedding.

Process a list of molecules:

.. code-block:: python

   from cosmolkit import BatchErrorMode, MoleculeBatch

   batch = MoleculeBatch.from_smiles_list(
       ["CCO", "c1ccccc1", "not-smiles"],
       errors="keep",
   ).with_parallel_jobs(8)

   for error in batch.errors():
       print(error.index(), error.operation(), error.message())

   prepared = batch.compute_2d_coords(errors=BatchErrorMode.KEEP)
   fingerprints = prepared.fingerprint_morgan_list(n_bits=2048)

   print(prepared.valid_mask())
   print(prepared.to_smiles_list(canonical=True))
   print([mol.to_smiles() if mol is not None else None for mol in prepared])
   print([fp.on_bits() if fp is not None else None for fp in fingerprints])
