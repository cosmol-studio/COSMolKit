Molecule Values
===============

``Molecule`` objects behave as copy-on-write (COW) molecule values.
Transformation methods return new molecule objects and leave the original
object unchanged. This is intentionally different from common RDKit Python
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
- ``with_2d_coords()``

In-Place Operations
-------------------

Performance-sensitive code can opt into explicit in-place mutation. Every
public ``Molecule`` in-place method ends with ``_``; the trailing underscore has
no other ``Molecule`` API meaning.

.. code-block:: python

   mol = Molecule.from_smiles("CCO")
   mol.add_hydrogens_()
   mol.compute_2d_coords_()

Common in-place operations include:

- ``add_hydrogens_()``
- ``remove_hydrogens_()``
- ``kekulize_()``
- ``sanitize_()``
- ``compute_2d_coords_()``

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

   mol = Molecule.from_smiles("c1ccccc1O").with_2d_coords()

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
