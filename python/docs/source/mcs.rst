Maximum common substructure
==========================

.. meta::
   :description: Find maximum common molecular substructures and inspect query graphs with COSMolKit.

This experimental search finds a common **query graph**, not a concrete
molecule. It requires at least two molecules and never changes its inputs.
Full RDKit FMCS parity and comparable performance have not been established.
Result queries preserve RDKit's source-first atom and bond alternative order;
the SMARTS comparison is exact, without sorting equivalent alternatives.

.. code-block:: python

   import cosmolkit as ck

   molecules = [ck.Molecule.from_smiles(s) for s in ("CCO", "CCN")]
   result = ck.maximum_common_substructure(molecules)
   print(result.atom_count, result.bond_count)  # 2 1
   print(result.smarts)                        # [#6]-[#6]
   assert all(m.has_substruct_match(result.query) for m in molecules)

Use a writable parameter object or configuration keywords:

.. code-block:: python

   params = ck.McsParameters(timeout=10)
   params.atom_compare_parameters.match_formal_charge = True
   result = ck.maximum_common_substructure(molecules, params)
   result = ck.maximum_common_substructure(molecules, timeout=10)

Do not combine a parameter instance with configuration keywords. Unknown
options raise an error. Enum inputs accept their member or lowercase spelling:
``atom_comparator=ck.McsAtomComparator.Any`` or ``atom_comparator="any"``.
Atom modes are ``any``, ``elements`` (default), ``isotopes`` and
``any_heavy_atom``; bond modes are ``any``, ``order`` (default) and
``order_exact``. Nested comparison objects expose the existing valence,
charge, isotope, chirality, distance, ring and bond-stereo options.

``threshold`` defaults to 1.0 (all inputs); ``timeout`` is in seconds, with
0 meaning no time limit. A timeout returns the best available partial result
with ``completed=False``. An empty completed result has zero counts, empty
SMARTS and ``query=None``. ``store_all=True`` retains equal-best alternatives
in ``degenerate``; the owner may then leave the singular query/SMARTS empty.
Result queries are independent snapshots. Options needing unavailable valence,
ring or 3D state raise ``McsError`` rather than silently preparing or modifying
the input molecule.

Rust uses ``cosmolkit::maximum_common_substructure(&[&a, &b])`` or
``maximum_common_substructure_with_params(&[&a, &b], &params)`` behind
``search``. JavaScript projects the same functions as
``maximumCommonSubstructure([a, b], params)`` and accepts a plain options
object such as ``{timeout: 10}``. Nested JavaScript parameter views are writable.

.. automodule:: cosmolkit
   :members: maximum_common_substructure, McsParameters, McsAtomCompareParameters, McsBondCompareParameters, McsAtomComparator, McsBondComparator, McsResult
   :no-index:
