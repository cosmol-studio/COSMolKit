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

   import cosmolkit as ck

   structure = ck.BioStructure.read("complex.pdb")
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

   import cosmolkit as ck

   protein = ck.Protein.read("1crn.pdb")

   print(protein.num_models())
   print(protein.num_chains())
   print(protein.num_residues())
   print(protein.num_atoms())

Read PDB text that is already in memory:

.. code-block:: python

   protein = ck.Protein.from_pdb(pdb_text)

Read mmCIF input with the same high-level protein projection:

.. code-block:: python

   protein = ck.Protein.read("1crn.cif")
   protein = ck.Protein.from_mmcif(cif_text)

``Protein`` keeps amino-acid residues and excludes ligands, nucleic acids, and
waters. Use it for protein-focused traversal rather than mixed structural
data. When those removed rows matter, start from ``BioStructure`` instead.

Chains, Residues, And Atoms
---------------------------

``Protein`` behaves like a chain collection. ``len(protein)`` returns the
number of protein chains, and ``protein[i]`` returns a ``ProteinChainRef``.

.. code-block:: python

   first_chain = protein[0]
   print(first_chain.id(), first_chain.kind(), len(first_chain))

   for chain in protein.chains():
       for residue in chain.residues():
           if residue.code() == ck.ResidueCode.MET:
               print("methionine", residue.id(), residue.fasta_code())
           print(residue.id(), residue.name(), residue.code(), len(residue))

           for atom in residue.atoms():
               print(atom.id(), atom.name(), atom.element(), atom.position())

``atom.position()`` returns ``None`` when the atom has no Cartesian coordinate
in the selected structure data; otherwise it returns ``(x, y, z)``.

Residue Information
-------------------

``ProteinResidueRef.name()`` returns the raw residue name from the structure.
Use ``ProteinResidueRef.code()`` for enum matching against Gemmi's tabulated
residue vocabulary, and ``ProteinResidueRef.info()`` when you need the
source-derived classification fields. Sequence expansion follows Gemmi's
``expand_one_letter`` and ``expand_one_letter_sequence`` residue tables.

.. code-block:: python

   import cosmolkit as ck

   info = ck.find_residue_info("MSE")
   assert info.code() == ck.ResidueCode.MSE
   assert info.kind() == ck.ResidueInfoKind.AA
   assert info.fasta_code() == "X"
   assert info.canonical_one_letter_code() == "M"
   assert info.parent_standard_code() == ck.ResidueCode.MET
   assert info.is_modified_amino_acid()
   assert ck.expand_one_letter_sequence("ACD(MSE)", ck.ResidueInfoKind.AA) == [
       "ALA",
       "CYS",
       "ASP",
       "MSE",
   ]

``fasta_code()`` deliberately follows Gemmi and emits ``"X"`` for modified
residues. For rendering, secondary-structure logic, or residue-family tests,
use ``canonical_one_letter_code()`` or ``parent_standard_code()`` instead. For
example, ``HYP`` maps to ``PRO`` and ``SEP`` maps to ``SER`` without changing
the raw residue name returned by ``ProteinResidueRef.name()``.

BioStructure vs Protein vs Molecule
-----------------------------------

Use ``BioStructure.from_pdb(text)`` or ``BioStructure.from_mmcif(text)`` when the
complete modeled structural hierarchy must remain available.

Use ``Protein.read(path)`` or ``Protein.from_pdb(text)`` when the desired
object is intentionally a protein-only structural view:

.. code-block:: python

   protein = ck.Protein.read("input.pdb")

Read PDB text into a ``BioStructure`` and call ``to_molecule()`` when the
desired object is an RDKit-compatible molecular graph:

.. code-block:: python

   import cosmolkit as ck

   structure = ck.BioStructure.from_pdb(pdb_text)
   mol = structure.to_molecule(
       sanitize=True,
       remove_hs=True,
       proximity_bonding=True,
   )

The molecule conversion path is useful for cheminformatics-style molecule
operations. The ``Protein`` path is the ergonomic path for protein chain,
residue, and atom access.
