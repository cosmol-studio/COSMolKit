Quick Start
===========

Value-Style Molecule Values
---------------------------

COSMolKit molecules use value semantics. Transform methods return a new
``Molecule`` and leave the original object unchanged. Internally the library
uses copy-on-write (COW) storage to share unchanged data efficiently:

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
   mol.compute_2d_coordinates_()

Create a molecule from SMILES and export a depiction:

.. code-block:: python

   from cosmolkit import Molecule

   mol = Molecule.from_smiles("c1ccccc1O")
   drawn = mol.with_2d_coordinates()

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

Read and write the first SDF record:

.. code-block:: python

   mol = Molecule.read_sdf("input.sdf", coordinate_dim="auto")
   mol.write_sdf("python/examples/output/output.sdf", format="v2000")

Read MOL2 with the RDKit-style parser profile:

.. code-block:: python

   mol = Molecule.read_mol2("input.mol2")

Access coordinates as NumPy arrays:

.. code-block:: python

   mol2d = Molecule.from_smiles("CCO").with_2d_coordinates()
   coords = mol2d.coords_2d()

   print(coords.shape)

Generate a native 3D conformer with ETKDGv3:

.. code-block:: python

   from cosmolkit import EmbedParameters, Molecule

   mol = Molecule.from_smiles("CC(=O)NC").with_hydrogens()
   params = EmbedParameters.etkdg_v3()
   params.random_seed = 0xF00D
   params.num_threads = 1
   params.track_failures = True

   embedded = mol.with_3d_conformer(params)

   print(embedded.num_conformers())
   print(embedded.coords_3d().shape)
   print(params.failures)

Generate multiple conformers with RMS pruning:

.. code-block:: python

   params = EmbedParameters.etkdg()
   params.random_seed = 123
   params.num_threads = 1
   params.prune_rms_thresh = 0.5
   params.enable_sequential_random_seeds = True

   conformers = mol.with_3d_conformers(5, params)
   print(conformers.num_conformers())

Optimize an existing 3D conformer with UFF:

.. code-block:: python

   mol = Molecule.read_sdf("input_3d.sdf", coordinate_dim="3d")

   if mol.has_uff_params():
       result = mol.with_uff_optimized(max_iters=200)
       optimized = result.molecule()

       print(not result.needs_more())
       print(result.status_code())
       print(result.energy())
       print(optimized.coords_3d().shape)

   if mol.has_mmff_params():
       result = mol.with_mmff_optimized(mmff_variant="MMFF94", max_iters=200)
       print(not result.needs_more())
       print(result.status_code())

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

Parse SMARTS metadata:

.. code-block:: python

   import cosmolkit

   query = cosmolkit.parse_smarts("[#6]-O")

   print(query.num_atoms())
   print(query.num_bonds())

``parse_smarts()`` returns a ``SmartsMolecule`` parse tree value. Direct SMARTS
query matching is not yet a Python API.

Process a list of molecules:

.. code-block:: python

   from cosmolkit import BatchErrorMode, MoleculeBatch

   batch = MoleculeBatch.from_smiles_list(
       ["CCO", "c1ccccc1", "not-smiles"],
       errors=BatchErrorMode.KEEP,
   ).with_parallel_jobs(8)

   for error in batch.errors():
       print(error.index(), error.operation(), error.message())

   prepared = batch.with_2d_coordinates(errors=BatchErrorMode.KEEP)
   fingerprints = prepared.fingerprint_morgan_list(n_bits=2048)

   print(prepared.valid_mask())
   print(prepared.to_smiles_list(canonical=True))
   print([mol.to_smiles() if mol is not None else None for mol in prepared])
   print([fp.on_bits() if fp is not None else None for fp in fingerprints])

``MoleculeBatch`` preserves input order. When ``errors="keep"`` or
``BatchErrorMode.KEEP`` is used, invalid records stay aligned with their input
positions and appear as ``None`` in molecule-valued outputs.
