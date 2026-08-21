Molecular Descriptors
=====================

.. meta::
   :description: Calculate molecular weight, formula, H-bond counts, Crippen LogP and MR, TPSA, aromatic rings, rotatable bonds, fraction Csp3, and QED with COSMolKit.

COSMolKit exposes the source-backed molecular descriptor functions from the
Rust core directly in Python. Descriptor calls are read-only and do not mutate
the input :class:`cosmolkit.Molecule`.

The descriptor surface is experimental. Supported rows are compared
field-by-field, including exact floating-point bit patterns, against pinned
RDKit golden data. A source state outside the modeled boundary raises
``NotImplementedError`` instead of returning an approximate value.

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

The complete public descriptor set is:

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
