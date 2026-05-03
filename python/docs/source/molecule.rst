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

SMILES Output
-------------

``to_smiles()`` returns a SMILES string:

.. code-block:: python

   mol = Molecule.from_smiles("F[C@H](Cl)Br")

   print(mol.to_smiles())
   print(mol.to_smiles(isomeric_smiles=False))

Explicit Editing
----------------

Use ``Molecule.edit()`` when you want to stage changes and commit them as one
new molecule:

.. code-block:: python

   editor = mol.edit()
   cl = editor.add_atom("Cl")
   editor.add_bond(0, cl, order="single")
   editor.set_atom_charge(cl, -1)

   edited = editor.commit()

Depictions
----------

Molecules with 2D coordinates can be exported as SVG or PNG:

.. code-block:: python

   mol = Molecule.from_smiles("c1ccccc1O").with_2d_coords()

   svg = mol.to_svg(width=400, height=300)
   mol.write_svg("phenol.svg", width=400, height=300)
   mol.write_png("phenol.png", width=400, height=300)

Stereo
------

COSMolKit keeps the atom-level CW/CCW chiral tag path available. This is the
closest representation to the explicit chiral information carried by SMILES or
RDKit atoms:

.. code-block:: python

   mol = Molecule.from_smiles("F[C@H](Cl)Br")

   for atom in mol.atoms():
       if atom.chiral_tag() != "CHI_UNSPECIFIED":
           print(atom.idx(), atom.chiral_tag())

   print(mol.find_chiral_centers(include_unassigned=False))

When code needs COSMolKit's ordered-ligand tetrahedral representation, use
``tetrahedral_stereo()`` as a separate view derived from those chiral tags:

.. code-block:: python

   mol = Molecule.from_smiles("F[C@H](Cl)Br")

   print(mol.tetrahedral_stereo())
