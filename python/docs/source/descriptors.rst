Molecular Descriptors
=====================

.. meta::
   :description: Calculate source-backed molecular properties, connectivity and shape indices, Lipinski counts, MQN, Labute ASA, and SlogP/SMR VSA descriptors with COSMolKit.

COSMolKit exposes the source-backed molecular descriptor methods on
:class:`cosmolkit.Molecule`. Descriptor calls are read-only and do not mutate
the input :class:`cosmolkit.Molecule`.

The implementations follow the pinned RDKit source. See ``VALIDATION.md`` in
the source repository for validation coverage.

Basic Descriptors
-----------------

.. code-block:: python

   import cosmolkit as ck

   molecule = ck.mol_from_smiles("c1ccccc1O")

   formula = molecule.molecular_formula()
   average_weight = molecule.molecular_weight()
   exact_weight = molecule.exact_molecular_weight()
   aromatic_rings = molecule.num_aromatic_rings()

   print(formula, average_weight, exact_weight, aromatic_rings)

The basic property functions are:

* :meth:`cosmolkit.Molecule.molecular_weight`
* :meth:`cosmolkit.Molecule.exact_molecular_weight`
* :meth:`cosmolkit.Molecule.molecular_formula`
* :meth:`cosmolkit.Molecule.num_hbd`
* :meth:`cosmolkit.Molecule.num_hba`
* :meth:`cosmolkit.Molecule.fraction_csp3`
* :meth:`cosmolkit.Molecule.crippen_descriptors`
* :meth:`cosmolkit.Molecule.tpsa`
* :meth:`cosmolkit.Molecule.num_aromatic_rings`
* :meth:`cosmolkit.Molecule.num_rotatable_bonds`
* :meth:`cosmolkit.Molecule.qed`

Connectivity And Shape
----------------------

Connectivity descriptors include graph-degree ``chi_0()`` and
``chi_1()``, generic order-N ``chi_n_v()`` and ``chi_n_n()``, and
the fixed ``chi_0_v()`` through ``chi_4_v()`` and ``chi_0_n()``
through ``chi_4_n()`` projections.

Hall-Kier and shape functions are:

* :meth:`cosmolkit.Molecule.hall_kier_alpha`
* :meth:`cosmolkit.Molecule.hall_kier_alpha_with_contributions`
* :meth:`cosmolkit.Molecule.kappa_1`
* :meth:`cosmolkit.Molecule.kappa_2`
* :meth:`cosmolkit.Molecule.kappa_3`
* :meth:`cosmolkit.Molecule.phi`

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

``mqns()`` returns the fixed source-order 42-entry molecular quantum
number vector. ``labute_asa()`` returns the scalar surface area, while
``labute_asa_contributions()`` returns a ``LabuteAsaContributions`` value
with ``asa``, ``atom_contributions`` (in atom-index order), and
``hydrogen_contribution`` properties.

``slogp_vsa()`` and ``smr_vsa()`` return 12-bin and 10-bin vectors.
Use ``slogp_vsa_with_params(bins, force)`` or
``smr_vsa_with_params(bins, force)`` for custom boundaries; the result then
has ``len(bins) + 1`` entries. Scalar projections ``slogp_vsa_1()`` through
``slogp_vsa_12()`` and ``smr_vsa_1()`` through
``smr_vsa_10()`` delegate to the same vector cores.

.. code-block:: python

   import cosmolkit as ck

   molecule = ck.mol_from_smiles("CC(O)c1ccncc1")

   chi = [molecule.chi_n_v(order) for order in range(5)]
   mqns = molecule.mqns()
   surface = molecule.labute_asa_contributions()
   asa = surface.asa
   atom_asa = surface.atom_contributions
   hydrogen_asa = surface.hydrogen_contribution
   slogp_vsa = molecule.slogp_vsa()

   assert len(mqns) == 42
   assert len(atom_asa) == molecule.num_atoms()
   assert len(slogp_vsa) == 12

Formula Options
---------------

``molecular_formula()`` uses the default options. Use
``molecular_formula_with_params(separate_isotopes, abbreviate_h_isotopes)``
for explicit options. When isotope separation and hydrogen abbreviation
are enabled, hydrogen-2 and hydrogen-3 are written as D and T.

Rotatable-Bond Modes
--------------------

``num_rotatable_bonds()`` uses the default mode. For explicit selection, use
``num_rotatable_bonds_with_params(ck.RotatableBondsOptions.Strict)`` or pass
``"default"``, ``"non_strict"``, ``"strict"``, or ``"strict_linkages"`` to
that method. Unknown modes raise ``ValueError``.
