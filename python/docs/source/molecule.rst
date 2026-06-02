Molecule Values
===============

``Molecule`` objects behave as value-style molecule values. Transformation
methods return new molecule objects and leave the original object unchanged.
Internally COSMolKit uses copy-on-write (COW) storage to share unchanged data
efficiently. This is intentionally different from common RDKit Python
workflows, where code often mutates an existing molecule or ``RWMol`` directly.

.. code-block:: python

   from cosmolkit import Molecule

   mol = Molecule.from_smiles("CCO")
   mol_h = mol.with_hydrogens()

   assert mol is not mol_h
   print(mol.to_smiles())
   print(mol_h.to_smiles())

Do not write code that assumes ``mol.with_hydrogens()`` changes ``mol``. Keep
the returned value and pass that value to later operations.

Common transformations include:

- ``with_hydrogens()``
- ``without_hydrogens()``
- ``with_kekulized_bonds()``
- ``with_2d_coordinates()``

In-Place Operations
-------------------

Performance-sensitive code can opt into explicit in-place mutation. Every
public ``Molecule`` in-place method ends with ``_``; the trailing underscore has
no other ``Molecule`` API meaning.

.. code-block:: python

   mol = Molecule.from_smiles("CCO")
   mol.add_hydrogens_()
   mol.compute_2d_coordinates_()

Common in-place operations include:

- ``add_hydrogens_()``
- ``remove_hydrogens_()``
- ``kekulize_()``
- ``sanitize_()``
- ``compute_2d_coordinates_()``

If an in-place method returns an error, the molecule is not guaranteed to equal
its pre-call value. Use the value-style method when failure-preserving behavior
is required.

SMILES Output
-------------

``to_smiles()`` returns a SMILES string:

.. code-block:: python

   mol = Molecule.from_smiles("F[C@H](Cl)Br")

   print(mol.to_smiles())
   print(mol.to_smiles(isomeric_smiles=False))

SMILES writer options are available on both single molecules and batches:

.. code-block:: python

   benzene = Molecule.from_smiles("c1ccccc1")
   ethanol = Molecule.from_smiles("CCO")

   print(benzene.to_smiles(kekule=True))
   print(ethanol.to_smiles(all_bonds_explicit=True))
   print(ethanol.to_smiles(canonical=False, rooted_at_atom=2))

Explicit Editing
----------------

Use ``Molecule.edit()`` when you want to stage changes and commit them as one
new molecule:

.. code-block:: python

   editor = mol.edit()
   cl = editor.add_atom("Cl")
   editor.add_bond(0, cl, order="single")

   edited = editor.commit()

Depictions
----------

Molecules with 2D coordinates can be exported as SVG or PNG:

.. code-block:: python

   mol = Molecule.from_smiles("c1ccccc1O").with_2d_coordinates()

   svg = mol.to_svg(width=400, height=300)
   mol.write_svg("python/examples/output/phenol.svg", width=400, height=300)
   mol.write_png("python/examples/output/phenol.png", width=400, height=300)

Stereo
------

COSMolKit keeps the atom-level CW/CCW chiral tag path available. This is the
closest representation to the explicit chiral information carried by SMILES or
RDKit atoms:

.. code-block:: python

   from cosmolkit import ChiralTag, Molecule

   mol = Molecule.from_smiles("F[C@H](Cl)Br")

   for atom in mol.atoms():
       if atom.chiral_tag() != ChiralTag.CHI_UNSPECIFIED:
           print(atom.idx(), atom.chiral_tag().name)

   print(mol.find_chiral_centers(include_unassigned=False))

Atom and bond enum-valued fields return Python ``IntEnum`` members, so callers
can compare or match against ``ChiralTag``, ``BondOrder``, ``BondDirection``,
and ``BondStereo`` instead of spelling chemistry states as strings. Read-only
maps such as ``BOND_ORDER_MAP`` and ``CHIRAL_TAG_MAP`` are available when a
string name from an external source needs to be converted to the enum member.

When code needs COSMolKit's ordered-ligand tetrahedral representation, use
``tetrahedral_stereo()`` as a separate view derived from those chiral tags:

.. code-block:: python

   mol = Molecule.from_smiles("F[C@H](Cl)Br")

   print(mol.tetrahedral_stereo())

Conformer Generation And Force-Field Optimization
-------------------------------------------------

Conformer generation APIs create native 3D conformers through the source-ported
distance-geometry path. The default value-style operation uses ETKDGv3 and
returns a new molecule value.

.. code-block:: python

   from cosmolkit import EmbedParameters, Molecule

   mol = Molecule.from_smiles("CC(=O)NC").with_hydrogens()

   params = EmbedParameters.etkdg_v3()
   params.random_seed = 0xF00D
   params.num_threads = 1
   params.track_failures = True

   embedded = mol.with_3d_conformer(params)

   print(embedded.num_conformers())
   print(embedded.coordinates_3d())
   print(params.failures)

For multi-conformer generation, explicit seeds are deterministic. RMS pruning,
sequential seed expansion, and terminal-group symmetrization for pruning follow
the source-ported RDKit path.

.. code-block:: python

   params = EmbedParameters.etkdg()
   params.random_seed = 123
   params.num_threads = 1
   params.prune_rms_thresh = 0.5
   params.enable_sequential_random_seeds = True

   pruned = mol.with_3d_conformers(5, params)
   print(pruned.num_conformers())

UFF and MMFF optimization APIs operate on existing or generated 3D conformers
and return new molecule values through result objects. They do not mutate the
source molecule.

.. code-block:: python

   from cosmolkit import Molecule

   mol = Molecule.from_smiles("CCO").with_hydrogens().with_3d_conformer()

   if mol.has_uff_params():
       result = mol.with_uff_optimized(max_iters=200)
       optimized = result.molecule()

       print(not result.needs_more())
       print(result.status_code())
       print(result.energy())
       print(optimized.coordinates_3d())

   if mol.has_mmff_params():
       result = mol.with_mmff_optimized(mmff_variant="MMFF94", max_iters=200)
       optimized = result.molecule()

       print(not result.needs_more())
       print(result.status_code())

Substructure And SMARTS
-----------------------

Substructure matching functions accept molecule queries:

.. code-block:: python

   import cosmolkit

   mol = Molecule.from_smiles("CCO")
   query = Molecule.from_smiles("CO")

   print(cosmolkit.has_substruct_match(mol, query))
   print(cosmolkit.get_substruct_match(mol, query).atom_mapping())

``parse_smarts()`` exposes the Rust SMARTS parser as parse metadata. It returns
a ``SmartsMolecule`` query-tree value. Direct SMARTS query matching is not yet
a Python API; Python substructure functions currently accept molecule queries.

.. code-block:: python

   smarts = cosmolkit.parse_smarts("[#6]-O")

   print(smarts.num_atoms())
   print(smarts.num_bonds())
