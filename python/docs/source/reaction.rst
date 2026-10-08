Reactions
=========

The experimental reaction API projects the existing Rust implementation.
Python and JavaScript do not implement separate reaction chemistry.

Generate products
-----------------

.. code-block:: python

   import cosmolkit as ck

   reaction = ck.Reaction.from_smirks("[C:1]>>[N:1]")
   source = ck.Molecule.from_smiles("C")
   product_sets = source.reaction_products(reaction, 0)
   assert product_sets[0][0].to_smiles() == "N"
   assert source.to_smiles() == "C"
   assert reaction.is_initialized()

``ck.parse_smirks`` is an equivalent factory. Explicit parse/write options use
the ``_with_params`` methods and ``ReactionParseParams`` /
``ReactionWriteParams``. ``to_smirks()`` and ``to_cx_smirks()`` return strings.
As in the source overloads, the former fixes CX output off and the latter
fixes it on; ``include_cx`` does not override that choice. ``cx_fields``
selects the fields emitted by the CX overload.

Products remain a list of product sets: each inner list preserves product
template order for one matching combination. No match returns an empty outer
list. Execution may initialize the same reaction handle.

For multiple inputs, preserve reactant-template order. The receiver remains
the canonical reconstruction anchor:

.. code-block:: python

   pair = ck.Reaction.from_smirks("[C:1].[C:2]>>[C:1].[C:2]")
   products = source.reaction_products_from_inputs(
       pair, [source, source], ck.ReactionRunParams(max_products=1000)
   )
   assert [m.to_smiles() for m in products[0]] == ["C", "C"]

``ReactionSingleRunParams`` selects coordinates for one input;
``ReactionRunParams.coordinate_selections`` selects per-reactant coordinates.
Use ``ReactionCoordinateSelection.auto()``, ``two_d(id)`` or ``three_d(id)``.
Explicit IDs are conformer IDs, not vector positions.

Restricted application
----------------------

.. code-block:: python

   result = source.apply_reaction(reaction)
   assert result.changed
   assert result.molecule.to_smiles() == "N"
   assert source.to_smiles() == "C"
   changed = source.apply_reaction_(reaction)
   assert changed and source.to_smiles() == "N"

Value application returns ``ReactionApplyResult`` with read-only ``molecule``
and ``changed`` properties. The ``_`` variant changes the receiver and returns
the source-defined boolean, not a graph-difference approximation. Explicit
parameter forms accept ``ReactionApplyParams``. Application requires one
reactant and one product template and cannot add new product atoms; use
product generation for other supported reactions. Failure is atomic.

Templates, validation and errors
--------------------------------

``Reaction.from_templates(reactants, products, agents)`` accepts ``QueryGraph``
values. Ordered template queries and ``with_*_template`` return detached
values. ``without_agents`` and ``without_unmapped_*`` return
``ReactionTemplateRemoval`` with ``reaction`` and ``removed_templates``.

``validate()`` returns ``ReactionValidationReport`` with warning/error issues,
counts and ``is_valid``. ``with_initialized()`` returns an initialized copy
without changing the source reaction. Explicit forms accept
``ReactionValidationParams(silent=True)``.

Options/results are read-only. Errors retain ``domain``, ``kind`` and variant
context such as ``role``, ``template``, ``index``, ``count`` and ``atom``.
``OperationError.__cause__`` retains ``ReactionRunError`` or
``ReactionApplyError``. Invalid initialization includes a validation ``report``.
These bindings do not establish full Reaction/SMIRKS corpus parity; existing
source-backed limitations and order-sensitive discrepancies remain unchanged.

JavaScript
----------

Names use camelCase: ``Reaction.fromSmirks``, ``molecule.reactionProducts``,
``reaction.withInitialized`` and ``ReactionCoordinateSelection.threeD``.
Value application is ``molecule.applyReaction``; mutation is
``molecule.applyReaction_``. Typed error payloads are under ``error.detail``
and concrete causes under ``error.cause``.
