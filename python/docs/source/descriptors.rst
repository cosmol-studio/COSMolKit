Molecular Descriptors
=====================

.. meta::
   :description: Calculate source-backed molecular properties, connectivity and shape indices, Lipinski counts, MQN, Labute ASA, and SlogP/SMR VSA descriptors with COSMolKit.

COSMolKit exposes the source-backed molecular descriptor functions from the
Rust core directly in Python. Descriptor calls are read-only and do not mutate
the input :class:`cosmolkit.Molecule`.

The documented descriptor surface is supported with pinned-RDKit parity.
Validation compares scalar values, complete vectors, atom contributions,
custom bins, and cold/warm/forced cache sequences, including exact
floating-point bit patterns. See ``VALIDATION.md`` in the source repository
for the complete evidence boundary.

Basic Descriptors
-----------------

.. code-block:: python

   import cosmolkit

   molecule = cosmolkit.Molecule.from_smiles("c1ccccc1O")

   formula = cosmolkit.calc_mol_formula(molecule)
   average_weight = cosmolkit.calc_mol_wt(molecule)
   exact_weight = cosmolkit.calc_exact_mol_wt(molecule)
   aromatic_rings = cosmolkit.calc_num_aromatic_rings(molecule)

   print(formula, average_weight, exact_weight, aromatic_rings)

The basic property functions are:

* :func:`cosmolkit.calc_mol_wt`
* :func:`cosmolkit.calc_exact_mol_wt`
* :func:`cosmolkit.calc_mol_formula`
* :func:`cosmolkit.calc_num_hbd`
* :func:`cosmolkit.calc_num_hba`
* :func:`cosmolkit.calc_fraction_csp3`
* :func:`cosmolkit.calc_crippen_descriptors`
* :func:`cosmolkit.calc_tpsa`
* :func:`cosmolkit.calc_num_aromatic_rings`
* :func:`cosmolkit.calc_num_rotatable_bonds`
* :func:`cosmolkit.calc_qed`

Connectivity And Shape
----------------------

Connectivity descriptors include graph-degree ``calc_chi_0()`` and
``calc_chi_1()``, generic order-N ``calc_chi_nv()`` and ``calc_chi_nn()``, and
the fixed ``calc_chi_0v()`` through ``calc_chi_4v()`` and ``calc_chi_0n()``
through ``calc_chi_4n()`` projections.

Hall-Kier and shape functions are:

* :func:`cosmolkit.calc_hall_kier_alpha`
* :func:`cosmolkit.calc_hall_kier_alpha_with_contributions`
* :func:`cosmolkit.calc_kappa_1`
* :func:`cosmolkit.calc_kappa_2`
* :func:`cosmolkit.calc_kappa_3`
* :func:`cosmolkit.calc_phi`

Lipinski And Ring Counts
------------------------

The extended count surface includes direct Lipinski N/O donor and acceptor
counts, heteroatoms, amide bonds, explicit heavy and total atom counts, SSSR
ring count, aromatic/aliphatic/saturated heterocycle and carbocycle counts,
spiro and bridgehead atoms, and possible or unspecified atom stereocenters.
All functions use the shared molecular graph, valence, ring, SMARTS, and stereo
implementations; there is no descriptor-local chemistry path.

MQN And Molecular Surface
-------------------------

``calc_mqns()`` returns the fixed source-order 42-entry molecular quantum
number vector. ``calc_labute_asa()`` returns the scalar surface area, while
``calc_labute_asa_contributions()`` returns the scalar, atom-index-aligned
contributions, and aggregate hydrogen contribution.

``calc_slogp_vsa()`` and ``calc_smr_vsa()`` return 12-bin and 10-bin vectors.
Pass ``bins=`` for custom boundaries; the result then has ``len(bins) + 1``
entries. Scalar projections ``calc_slogp_vsa_1()`` through
``calc_slogp_vsa_12()`` and ``calc_smr_vsa_1()`` through
``calc_smr_vsa_10()`` delegate to the same vector cores.

.. code-block:: python

   import cosmolkit

   molecule = cosmolkit.Molecule.from_smiles("CC(O)c1ccncc1")

   chi = [cosmolkit.calc_chi_nv(molecule, order) for order in range(5)]
   mqns = cosmolkit.calc_mqns(molecule)
   asa, atom_asa, hydrogen_asa = cosmolkit.calc_labute_asa_contributions(molecule)
   slogp_vsa = cosmolkit.calc_slogp_vsa(molecule)

   assert len(mqns) == 42
   assert len(atom_asa) == molecule.num_atoms()
   assert len(slogp_vsa) == 12

Formula Options
---------------

``calc_mol_formula()`` accepts ``separate_isotopes`` and
``abbreviate_h_isotopes``. When isotope separation and hydrogen abbreviation
are enabled, hydrogen-2 and hydrogen-3 are written as D and T.

Rotatable-Bond Modes
--------------------

``calc_num_rotatable_bonds()`` accepts ``mode="default"``,
``"non_strict"``, ``"strict"``, or ``"strict_linkages"``. Unknown modes raise
``ValueError``.
