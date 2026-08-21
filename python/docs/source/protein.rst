Protein Structures
==================

.. meta::
   :description: Read PDB and mmCIF protein structures with COSMolKit and traverse models, chains, residues, atoms, residue metadata, and modified amino-acid identities.

Use ``Protein`` when a workflow starts from PDB or mmCIF structural data and
needs protein-chain, residue, and atom traversal.

Read a PDB file directly:

.. code-block:: python

   from cosmolkit import Protein, ResidueCode

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
           if residue.code() == ResidueCode.MET:
               print("methionine", residue.index(), residue.fasta_code())
           print(residue.index(), residue.name(), residue.code(), len(residue))

           for atom in residue.atoms():
               print(atom.index(), atom.name(), atom.element(), atom.position())

``atom.position()`` returns ``None`` when the atom has no Cartesian coordinate
in the selected structure data; otherwise it returns ``(x, y, z)``.

Residue Information
-------------------

``ProteinResidue.name()`` returns the raw residue name from the structure.
Use ``ProteinResidue.code()`` for enum matching against Gemmi's tabulated
residue vocabulary, and ``ProteinResidue.info()`` when you need the
source-derived classification fields. Sequence expansion follows Gemmi's
``expand_one_letter`` and ``expand_one_letter_sequence`` residue tables.

.. code-block:: python

   from cosmolkit import (
       ResidueCode,
       ResidueInfoKind,
       expand_one_letter_sequence,
       find_tabulated_residue,
   )

   info = find_tabulated_residue("MSE")
   assert info.code() == ResidueCode.MSE
   assert info.kind() == ResidueInfoKind.AA
   assert info.fasta_code() == "X"
   assert info.canonical_one_letter_code() == "M"
   assert info.parent_standard_code() == ResidueCode.MET
   assert info.is_modified_amino_acid()
   assert expand_one_letter_sequence("ACD(MSE)", ResidueInfoKind.AA) == [
       "ALA",
       "CYS",
       "ASP",
       "MSE",
   ]

``fasta_code()`` deliberately follows Gemmi and emits ``"X"`` for modified
residues. For rendering, secondary-structure logic, or residue-family tests,
use ``canonical_one_letter_code()`` or ``parent_standard_code()`` instead. For
example, ``HYP`` maps to ``PRO`` and ``SEP`` maps to ``SER`` without changing
the raw residue name returned by ``ProteinResidue.name()``.

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
