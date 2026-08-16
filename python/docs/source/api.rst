API Reference
=============

.. automodule:: cosmolkit
   :members: Molecule, Atom, Bond, BondOrder, BondDirection, BondStereo, ChiralTag, ResidueCode, ResidueInfoKind, ResidueInfo, find_tabulated_residue, find_tabulated_residue_idx, get_residue_info, residue_code_from_name, expand_one_letter, expand_protein_one_letter, expand_one_letter_sequence, expand_protein_one_letter_string, EmbedParameters, MoleculeBatch, BatchErrorMode, Fingerprint, MorganAdditionalOutput, MorganFingerprintResult, UffOptimizeMoleculeResult, UffOptimizeMoleculeConfResult, UffOptimizeMoleculeConfsResult, MmffOptimizeMoleculeResult, MmffOptimizeMoleculeConfResult, MmffOptimizeMoleculeConfsResult, BatchError, BatchExportReport, MoleculeEdit, SdfDataset, SdfReader, SdfRecord, SdfRecordMetadata, BatchValidationError, InchiError, InchiAllocationError, InchiUnsupportedStateError, InchiDiagnosticWarning, inchi_to_key, SmartsMolecule, parse_smarts, SubstructMatchResult, uff_has_all_molecule_params, uff_optimize_molecule, uff_optimize_molecule_confs, mmff_has_all_molecule_params, mmff_optimize_molecule, mmff_optimize_molecule_confs, has_substruct_match, get_substruct_match, get_substruct_matches, get_substruct_matches_with_params, calc_mol_wt, calc_exact_mol_wt, calc_mol_formula, calc_num_hbd, calc_num_hba, calc_fraction_csp3, calc_crippen_descriptors, calc_tpsa, calc_num_aromatic_rings, calc_num_rotatable_bonds, calc_qed, mol_to_binary, mol_from_binary, version
   :no-show-inheritance:

InChI API
---------

The public InChI surface is limited to four scalar calls. Nonfatal source
diagnostics are emitted as :class:`cosmolkit.InchiDiagnosticWarning` instances
with ``level`` and ``message`` fields. Failures expose ``operation``, ``kind``,
and ``detail`` through :class:`cosmolkit.InchiError`; allocation and unsupported
state failures use dedicated subclasses.

Exact parity applies to behavior defined by pinned official InChI v1.07.5 and
RDKit 2026.03.1. The official C ``NormalizeAndCompare`` initial-buffer
allocation-failure path is undefined; COSMolKit returns a deterministic
``InchiAllocationError`` instead of claiming an exact C result. MolBlock,
SDF/V3000, IXA, AuxInfo, INCHIGEN, version-query, and extended-polymer InChI
entry points are not part of this public surface.

Typical usage keeps molecule conversion on :class:`cosmolkit.Molecule`::

.. code-block:: python

   from cosmolkit import Molecule, inchi_to_key

   molecule = Molecule.from_smiles("CCO")
   inchi = molecule.to_inchi()
   key = molecule.to_inchi_key()
   assert inchi_to_key(inchi) == key

   restored = Molecule.from_inchi(inchi)
   assert restored is not None

.. py:method:: cosmolkit.Molecule.to_inchi(options="")

   Return the molecule's InChI string without mutating it.

.. py:method:: cosmolkit.Molecule.to_inchi_key(options="")

   Return the molecule's InChIKey without mutating it.

.. py:function:: cosmolkit.inchi_to_key(inchi)

   Return the InChIKey for an InChI string, or ``None`` when the source API
   rejects the input.

.. py:classmethod:: cosmolkit.Molecule.from_inchi(inchi, *, sanitize=True, remove_hs=True)

   Return a :class:`cosmolkit.Molecule`, or ``None`` when the source API returns
   no graph.

Protein API
-----------

.. autoclass:: cosmolkit.Protein
   :members:
   :undoc-members:
   :no-show-inheritance:

.. autoclass:: cosmolkit.ProteinChain
   :members:
   :undoc-members:
   :no-show-inheritance:

.. autoclass:: cosmolkit.ProteinResidue
   :members:
   :undoc-members:
   :no-show-inheritance:

.. autoclass:: cosmolkit.ProteinAtom
   :members:
   :undoc-members:
   :no-show-inheritance:
