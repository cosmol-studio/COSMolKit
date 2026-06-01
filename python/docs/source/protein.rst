Protein Structures
==================

Use ``Protein`` when a workflow starts from PDB or mmCIF structural data and
needs protein-chain, residue, and atom traversal.

Read a PDB file directly:

.. code-block:: python

   from cosmolkit import Protein

   protein = Protein.from_pdb("1crn.pdb")

   print(protein.num_models())
   print(protein.num_chains())
   print(protein.num_residues())
   print(protein.num_atoms())

Read PDB text that is already in memory:

.. code-block:: python

   protein = Protein.from_pdb_str(pdb_text)

Read mmCIF input with the same high-level protein projection:

.. code-block:: python

   protein = Protein.from_mmcif("1crn.cif")
   protein = Protein.from_mmcif_str(cif_text, path="1crn.cif")

``Protein`` keeps amino-acid residues and excludes ligands, nucleic acids, and
waters by default. Use it for protein-focused traversal rather than low-level
mixed structural tables.

Chains, Residues, And Atoms
---------------------------

``Protein`` behaves like a chain collection. ``len(protein)`` returns the
number of protein chains, and ``protein[i]`` returns a ``ProteinChain``.

.. code-block:: python

   first_chain = protein[0]
   print(first_chain.index(), first_chain.kind(), len(first_chain))

   for chain in protein.chains():
       for residue in chain.residues():
           print(residue.index(), residue.name(), residue.kind(), len(residue))

           for atom in residue.atoms():
               print(atom.index(), atom.name(), atom.element(), atom.position())

``atom.position()`` returns ``None`` when the atom has no Cartesian coordinate
in the selected structure data; otherwise it returns ``(x, y, z)``.

Protein vs Molecule PDB APIs
----------------------------

Use ``Protein.from_pdb()`` or ``Protein.from_pdb_str()`` when the desired
object is a protein structural view:

.. code-block:: python

   protein = Protein.from_pdb("input.pdb")

Use ``Molecule.from_pdb_block()`` only when the desired object is a
RDKit-compatible molecule conversion from PDB text:

.. code-block:: python

   from cosmolkit import Molecule

   mol = Molecule.from_pdb_block(
       pdb_text,
       sanitize=True,
       remove_hs=True,
       proximity_bonding=True,
   )

The molecule conversion path is useful for cheminformatics-style molecule
operations. The ``Protein`` path is the ergonomic path for protein chain,
residue, and atom access.
