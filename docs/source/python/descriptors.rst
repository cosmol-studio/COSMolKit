Molecular Descriptors
=====================

.. meta::
   :description: Calculate molecular properties, Gasteiger charges, surface-area and EState descriptors, BCUT2D eigenvalues, and graph information content with COSMolKit.

Call descriptor methods on :class:`cosmolkit.Molecule`. They leave atoms,
bonds and coordinates unchanged. Some surface-area methods reuse previously
calculated atom areas, as described below. The repository's
`validation results <https://github.com/cosmol-studio/COSMolKit/blob/main/VALIDATION.md>`_
list the RDKit references and test corpora.

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

Gasteiger Partial Charges
-------------------------

``gasteiger_charges()`` returns one row per atom in input order. Each row
separates the charge on the atom from the total charge on its attached
hydrogens that are not graph vertices. Hydrogens already present as vertices
have their own rows. The default is 12 iterations with permissive parameter
lookup; strict lookup raises an error for an unsupported element and mode.

.. code-block:: python

   import cosmolkit as ck

   ethanol = ck.mol_from_smiles("CCO")
   strict = ck.GasteigerChargeParams(throw_on_param_failure=True)
   oxygen = ethanol.gasteiger_charges(strict)[2]
   assert round(oxygen.atom_charge, 3) == -0.397
   assert round(oxygen.attached_hydrogen_charge, 3) == 0.210
   assert len(ethanol.peoe_vsa()) == 14

Reuse ``strict`` across molecules, or pass settings as keywords for one call:
``ethanol.gasteiger_charges(throw_on_param_failure=True)``. Passing both a
parameter object and keywords raises ``TypeError``.

``peoe_vsa()`` groups Labute atom surface areas by atom charge. An earlier
call that explicitly recalculated the Labute areas without hydrogen
contributions changes those areas. The four charge extrema use only the
graph atoms' own charges; ``min_abs_partial_charge()`` is the smaller absolute
value of the signed minimum and maximum.

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
   assert len(atom_asa) == molecule.num_atoms
   assert len(slogp_vsa) == 12

BCUT2D
------

``bcut2d()`` returns four high/low eigenvalue pairs, in order: atomic mass,
Gasteiger charge, Crippen logP and Crippen molar refractivity. It removes
hydrogens and sanitizes a copy of the molecule. An unsupported bond order or
missing Gasteiger parameter raises ``DescriptorReadError`` for the whole query.
The unsupported-bond error reports the bond index in the input molecule.
With one atom left after hydrogen removal, both charge extrema equal that
atom's charge, even if it is infinite. With two or more atoms, a NaN or
infinite charge makes the two charge extrema NaN. The other pairs keep their values.
The sign bit of a NaN can differ from RDKit's for the same molecule.
The Rust eigensolver is shared with WASM; its final eigenvalues can differ
from RDKit's Eigen results in the last digits. The parity runner treats such
differences, including differing NaN bits, as failures.

.. code-block:: python

   from math import isclose
   from cosmolkit import Molecule

   values = Molecule.from_smiles("C.C").bcut2d()
   mass_high, mass_low = values[:2]
   assert isclose(mass_high, 12.012) and isclose(mass_low, 12.010)

Electrotopological state
------------------------

``estate_indices()`` returns one EState index per graph atom, in atom order.
``estate_vsa()`` groups Labute atom areas by EState (11 bins);
``vsa_estate()`` groups the EState indices by Labute atom area (10 bins).
The separate Labute hydrogen contribution is excluded in both directions.

.. code-block:: python

   from cosmolkit import Molecule

   ethanol = Molecule.from_smiles("CCO")
   indices = ethanol.estate_indices()
   assert indices[1] == 0.25
   assert ethanol.vsa_estate_3() == indices[2]
   assert ethanol.estate_vsa_6() == ethanol.estate_vsa()[5]

The default vectors refresh hydrogen-inclusive Labute rows. Numbered methods
reuse that direction's last vector; a molecule copy starts with no saved
EState vectors. The extrema (``min_estate_index()``, ``max_estate_index()``,
``min_abs_estate_index()``, ``max_abs_estate_index()``) raise a structured descriptor error on an empty graph.

Graph information content
-------------------------

``avg_ipc()`` returns the Shannon entropy in bits of the unweighted graph's
characteristic-polynomial coefficients, using their normalized magnitudes as
probabilities. ``ipc()`` multiplies that entropy by the sum of the coefficient
magnitudes. Bond order is ignored; a hydrogen represented as a graph atom
counts as a vertex. Neither query needs coordinates or sanitization.

.. code-block:: python

   from cosmolkit import Molecule

   molecule = Molecule.from_smiles("CC")
   assert molecule.ipc() == 2 and molecule.avg_ipc() == 1
   assert molecule.ipc([[0, 1], [0, 0]]) == 0
   other = Molecule.from_smiles("C.C").distance_matrix()
   assert molecule.avg_ipc(other) == 0

To replace the bonds for one call, supply a square matrix with one row and
column per atom, in graph order. Every entry exactly equal to one contributes
to the adjacency matrix; a diagonal one adds a self-loop. An asymmetric input
remains asymmetric. The default polynomial expands eigenvalues in increasing
magnitude order, following RDKit 2026.09.2 Python rounding and nonfinite behavior.
Both results are zero for empty or edgeless graphs.

Choose ``IpcMethod.ExactPolynomial`` explicitly to use the integer-coefficient
algorithm from RDKit PR #9884 (commit ``b9518c8a``):

.. code-block:: python

   import cosmolkit as ck

   molecule = ck.Molecule.from_smiles("C" * 120)
   value = molecule.avg_ipc_with_method(ck.IpcMethod.ExactPolynomial)
   print(value)  # 3.7592685780108055
   total = molecule.ipc_with_method(ck.IpcMethod.ExactPolynomial)

``IpcMethod.RdkitCompatible`` selects the same behavior as the default methods.
The exact method uses the molecule's unweighted graph, including explicit
hydrogen vertices. Its integer coefficients are exact; entropy is still
computed in floating point. Scaling keeps AvgIpc finite when the unscaled
coefficient sum overflows; total Ipc may be infinite. Exceeding the upstream
64-prime reconstruction table returns NaN. Existing matrix methods keep their
full-matrix contract and default algorithm.


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
