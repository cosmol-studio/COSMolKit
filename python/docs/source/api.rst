API Reference
=============

.. meta::
   :description: Complete COSMolKit Python API reference for molecules, atoms, bonds, fingerprints, descriptors, conformers, force fields, InChI, structural biology, IO, and batch processing.

.. automodule:: cosmolkit
   :members: Molecule, Atom, Bond, Element, ElementInfo, element_from_symbol, get_element_info, BondOrder, BondDirection, BondStereo, ChiralTag, ResidueCode, ResidueInfoKind, ResidueInfo, find_tabulated_residue, find_tabulated_residue_idx, get_residue_info, residue_code_from_name, expand_one_letter, expand_protein_one_letter, expand_one_letter_sequence, expand_protein_one_letter_string, EmbedParameters, StereoisomerOptions, StereoisomerIterator, PotentialStereoAnalysis, PotentialStereoInfo, MoleculeBatch, BatchErrorMode, Fingerprint, MorganAdditionalOutput, MorganFingerprintResult, TopologicalFingerprintResult, UffOptimizeMoleculeResult, UffOptimizeMoleculeConfResult, UffOptimizeMoleculeConfsResult, MmffOptimizeMoleculeResult, MmffOptimizeMoleculeConfResult, MmffOptimizeMoleculeConfsResult, BatchError, BatchExportReport, MoleculeEdit, SdfDataset, SdfReader, SdfRecord, SdfRecordMetadata, BatchValidationError, InchiError, InchiAllocationError, InchiUnsupportedStateError, InchiDiagnosticWarning, inchi_to_key, SmartsParserParams, parse_smarts, parse_smarts_with_params, SubstructMatchResult, uff_has_all_molecule_params, uff_optimize_molecule, uff_optimize_molecule_confs, mmff_has_all_molecule_params, mmff_optimize_molecule, mmff_optimize_molecule_confs, has_substruct_match, get_substruct_match, get_substruct_matches, get_substruct_matches_with_params, calc_mol_wt, calc_exact_mol_wt, calc_mol_formula, calc_num_hbd, calc_num_hba, calc_fraction_csp3, calc_crippen_descriptors, calc_tpsa, calc_num_aromatic_rings, calc_num_rotatable_bonds, calc_qed, calc_chi_0, calc_chi_1, calc_chi_nv, calc_chi_nn, calc_chi_0v, calc_chi_1v, calc_chi_2v, calc_chi_3v, calc_chi_4v, calc_chi_0n, calc_chi_1n, calc_chi_2n, calc_chi_3n, calc_chi_4n, calc_hall_kier_alpha, calc_hall_kier_alpha_with_contributions, calc_kappa_1, calc_kappa_2, calc_kappa_3, calc_phi, calc_lipinski_hba, calc_lipinski_hbd, calc_num_heteroatoms, calc_num_amide_bonds, calc_num_heavy_atoms, calc_num_atoms, calc_num_rings, calc_num_heterocycles, calc_num_saturated_rings, calc_num_aliphatic_rings, calc_num_aromatic_heterocycles, calc_num_aromatic_carbocycles, calc_num_aliphatic_heterocycles, calc_num_aliphatic_carbocycles, calc_num_saturated_heterocycles, calc_num_saturated_carbocycles, calc_num_spiro_atoms, calc_num_bridgehead_atoms, calc_num_atom_stereo_centers, calc_num_unspecified_atom_stereo_centers, calc_mqns, calc_labute_asa, calc_labute_asa_contributions, calc_slogp_vsa, calc_slogp_vsa_1, calc_slogp_vsa_2, calc_slogp_vsa_3, calc_slogp_vsa_4, calc_slogp_vsa_5, calc_slogp_vsa_6, calc_slogp_vsa_7, calc_slogp_vsa_8, calc_slogp_vsa_9, calc_slogp_vsa_10, calc_slogp_vsa_11, calc_slogp_vsa_12, calc_smr_vsa, calc_smr_vsa_1, calc_smr_vsa_2, calc_smr_vsa_3, calc_smr_vsa_4, calc_smr_vsa_5, calc_smr_vsa_6, calc_smr_vsa_7, calc_smr_vsa_8, calc_smr_vsa_9, calc_smr_vsa_10, mol_to_binary, mol_from_binary, version
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

Structural API
--------------

.. autoclass:: cosmolkit.MmcifOutputGroups
   :members:
   :no-show-inheritance:

.. autoclass:: cosmolkit.MmcifWriteOptions
   :members:
   :no-show-inheritance:

.. autoclass:: cosmolkit.BioStructure
   :members:
   :undoc-members:
   :no-show-inheritance:

.. autoclass:: cosmolkit.StructureModel
   :members:
   :no-show-inheritance:

.. autoclass:: cosmolkit.StructureChain
   :members:
   :no-show-inheritance:

.. autoclass:: cosmolkit.StructureResidue
   :members:
   :no-show-inheritance:

.. autoclass:: cosmolkit.StructureAtom
   :members:
   :no-show-inheritance:

.. autoclass:: cosmolkit.StructureEntity
   :members:
   :no-show-inheritance:

Protein Projection API
----------------------

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
