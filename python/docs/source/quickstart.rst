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

Assign atom chiral tags from a stored 3D conformer:

.. code-block:: python

   import numpy as np

   spatial = Molecule.from_smiles("C(F)(Cl)Br").with_only_3d_conformer(
       np.array(
           [
               [0.0, 0.0, 0.0],
               [1.0, 0.0, 0.0],
               [0.0, 1.0, 0.0],
               [0.0, 0.0, 1.0],
           ]
       )
   )
   assigned = spatial.with_chiral_tags_from_structure()

   assert spatial.atoms()[0].chiral_tag() == ChiralTag.CHI_UNSPECIFIED
   assert assigned.atoms()[0].chiral_tag() != ChiralTag.CHI_UNSPECIFIED

``with_chiral_tags_from_structure()`` is a stable, RDKit-parity operation. Its
in-place counterpart is ``assign_chiral_tags_from_structure_()``. The
value-style form leaves caller state unchanged on failure. Like other in-place
operations, the ``_`` form may retain valid partial changes after an error but
keeps the molecule's internal storage complete.

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
   coords = mol2d.coordinates_2d()

   print(coords.shape)

Round-trip a molecule through Python pickle:

.. code-block:: python

   import pickle

   mol = Molecule.from_smiles("F[C@H](Cl)[13CH3:7]").with_2d_coordinates()

   restored = pickle.loads(pickle.dumps(mol, protocol=pickle.HIGHEST_PROTOCOL))

   print(restored.to_smiles(canonical=False))
   print(restored.has_2d_coordinates())

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
   print(embedded.coordinates_3d().shape)
   print(params.failures)

``with_3d_conformer()`` follows RDKit's ETKDG behavior for trusted molecular
graphs. A molecule without explicit hydrogens is embedded as a heavy-atom-only
conformer; the operation does not fail and does not automatically add
hydrogens. Add hydrogens first when you need all-atom geometry, force-field
optimization, or hydrogen-bond-sensitive coordinates.

Generate multiple conformers with RMS pruning:

.. code-block:: python

   params = EmbedParameters.etkdg()
   params.random_seed = 123
   params.num_threads = 1
   params.prune_rms_thresh = 0.5
   params.enable_sequential_random_seeds = True

   conformers = mol.with_3d_conformers(5, params)
   print(conformers.num_conformers())

Decode a ConfSeq record with the reference distance-geometry backend and the
opt-in FastGeometry backend:

.. code-block:: python

   import cosmolkit

   # Tokenized corpus record: dude_aa2ar__L021087 candidate 0.
   # Spaces are ConfSeq token boundaries and must be preserved.
   confseq = (
       "N C <112> | ( = O ) <84> C <111> | <-45> n 1 c c ( <21> N <123> | "
       "<173> C <116> | ( = O ) <120> c 2 c c c c c 2 <112> C <112> | <-66> "
       "N 2 <-172> C ( = O ) <-2> C <0> N <3> C <174> 2 = O ) c n 1"
   )

   reference = cosmolkit.confseq.decode(
       confseq,
       optimize_with_uff=False,
       template_backend="distance_geometry",
   )

   fast = cosmolkit.confseq.decode(
       confseq,
       optimize_with_uff=False,
       template_backend="fast_geometry",
   )

   print(reference.num_conformers())
   print(fast.num_conformers())

Optimize an existing 3D conformer with UFF:

.. code-block:: python

   mol = Molecule.read_sdf("input_3d.sdf", coordinate_dim="3d")

   if mol.has_uff_params():
       result = mol.with_uff_optimized(max_iters=200)
       optimized = result.molecule()

       print(not result.needs_more())
       print(result.status_code())
       print(result.energy())
       print(optimized.coordinates_3d().shape)

   if mol.has_mmff_params():
       result = mol.with_mmff_optimized(mmff_variant="MMFF94", max_iters=200)
       print(not result.needs_more())
       print(result.status_code())

Generate source-backed fingerprints:

.. code-block:: python

   mol = Molecule.from_smiles("c1ccccc1O")
   morgan = mol.fingerprint_morgan(
       radius=2,
       n_bits=2048,
   )
   topological = mol.topological_fingerprint(fp_size=2048)
   avalon = mol.avalon_fingerprint(n_bits=512)

   print(morgan.on_bits())
   print(topological.on_bits())
   print(avalon.on_bits())

``on_bits()`` returns the sparse bit indexes set inside the fixed-length binary
fingerprint. The exposed Morgan, MACCS, RDKit topological, and Avalon branches
are covered by exact pinned-RDKit comparisons; structurally similar hashes or
similarity correlation are not compatibility claims. It is not a dense neural
embedding.

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
