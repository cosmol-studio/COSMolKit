Protein Structures
==================

.. meta::
   :description: Read PDB and mmCIF protein structures with COSMolKit and traverse models, chains, residues, atoms, residue metadata, and modified amino-acid identities.

Use ``BioStructure`` when a workflow must retain the complete modeled PDB or
mmCIF hierarchy, including nucleic acids, ligands, waters, and entities. Use
``Protein`` when the workflow intentionally needs only amino-acid chains,
residues, and atoms.

Complete Structures
-------------------

Read the complete structural value before selecting a projection:

.. code-block:: python

   from cosmolkit import BioStructure

   structure = BioStructure.from_pdb("complex.pdb")
   print(structure.num_models(), structure.num_entities())

   for model in structure.models():
       for chain in model.chains():
           for residue in chain.residues():
               print(residue.name(), residue.kind())

``structure.protein()`` returns an amino-acid-only projection without changing
``structure``. Structural child objects share the parent storage; traversal
does not copy the complete structure for every model, chain, residue, or atom.

Serialize the complete structural value as Gemmi-aligned mmCIF without
mutating it. Writer options control category groups and CIF formatting:

.. code-block:: python

   mmcif_text = structure.to_mmcif()
   structure.write_mmcif("roundtrip.cif")

These writers belong only to ``BioStructure``. ``Protein`` and ``Molecule`` do
not expose structural mmCIF writer aliases because they intentionally preserve
different state boundaries. The writer emits the state represented by
``BioStructure``; arbitrary source categories not modeled by that value are
not claimed to round-trip.

Protein Projections
-------------------

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
waters. Use it for protein-focused traversal rather than mixed structural
data. When those removed rows matter, start from ``BioStructure`` instead.

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

BioStructure vs Protein vs Molecule
-----------------------------------

Use ``BioStructure.from_pdb()`` or ``BioStructure.from_mmcif()`` when the
complete modeled structural hierarchy must remain available.

Use ``Protein.from_pdb()`` or ``Protein.from_pdb_str()`` when the desired
object is intentionally a protein-only structural view:

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
