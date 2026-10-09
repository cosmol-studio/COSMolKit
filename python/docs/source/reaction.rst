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
   product_sets = reaction.run([source], ck.ReactionRunParams())
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

For multiple inputs, preserve reactant-template order. Execution belongs to
the reaction, not an unrelated molecule receiver:

.. code-block:: python

   pair = ck.Reaction.from_smirks("[C:1].[C:2]>>[C:1].[C:2]")
   products = pair.run([source, source], ck.ReactionRunParams(max_products=1000))
   assert [m.to_smiles() for m in products[0]] == ["C", "C"]

Rust uses ``use cosmolkit::{Reaction, ReactionRunParams};`` and the inherent
``reaction.run(&[&mol_a, &mol_b], &params)`` method; no extension trait or
domain-crate import is needed. The public reaction type, constructors,
execution method and registry entries require ``cap-reaction``.

``ReactionSingleRunParams`` selects coordinates for one input;
``ReactionRunParams.coordinate_selections`` selects per-reactant coordinates.
Use ``ReactionCoordinateSelection.auto()``, ``two_d(id)`` or ``three_d(id)``.
Explicit IDs are conformer IDs, not vector positions.

Typed atom metadata
-------------------

User metadata retains its type instead of being converted to text:

.. code-block:: python

   carbon = ck.Molecule.from_smiles("C")
   tagged = carbon.with_atom_property(0, "tracking_id", 42)
   assert carbon.atom_property(0, "tracking_id") is None
   assert tagged.atom_property(0, "tracking_id") == 42
   tagged.set_atom_property_(0, "tracking_id", 43)

Values support ``bool``, ``int`` (signed 32-bit or unsigned 32-bit), ``float``,
``str``, homogeneous ``list[int]`` (signed 32-bit), and ``list[str]``.
An empty list is stored as a string list. Missing keys return ``None``;
invalid atom indices, reserved keys and unsupported value types raise errors.
Setters accept nonempty user keys, not underscore-prefixed private keys,
computed properties, atom-map numbers or reaction bookkeeping keys.
Clones and value-returning operations retain COW isolation.

Reaction copying is an explicit CK extension, disabled by default:

.. code-block:: python

   oxygen = ck.Molecule.from_smiles("O").with_atom_property(0, "tracking_id", 84)
   reaction = ck.Reaction.from_smirks("[C:1].[O:2]>>[C:1][O:2]")
   params = ck.ReactionRunParams(copy_atom_properties=True)
   product = reaction.run([tagged, oxygen], params)[0][0]
   assert [product.atom_property(i, "tracking_id") for i in range(2)] == [43, 84]

Copying follows both the input-reactant index and source-atom index, not
destination positions. Each duplicate inherits its source's user metadata;
new atoms inherit none and deleted atoms contribute none. Existing product
values, including explicit template values, win over copied values. Private,
computed, CIP and reaction bookkeeping properties are not additionally copied.
``False`` means no extra copying; it does not clear properties retained by the
normal reaction algorithm. Chemical assignment follows the normal pipeline.

JavaScript uses ``withAtomProperty``, ``setAtomProperty`` and ``atomProperty``;
missing keys return ``null``. Numbers with integral values use the integer
property types (except negative zero); other numbers use doubles. Arrays have
the same homogeneous value restrictions as Python. Set
``params.copyAtomProperties = true`` before ``reaction.run(inputs, params)``.

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

Results are read-only; ``ReactionRunParams`` fields are writable configuration.
Errors retain ``domain``, ``kind`` and variant
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
