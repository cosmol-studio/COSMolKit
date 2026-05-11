# Redesigned cosmolkit-core

This crate is the new canonical COSMolKit core.

Rules for future changes:

- `Molecule` is a value object.
- Public molecule transformations return a new value.
- New molecules are built through `MoleculeBuilder`.
- Existing molecule transformations must go through registered operations.
- Any new derived state must be represented in `DerivedState` and checked by
  invariants.
- RDKit behavior belongs in a compatibility layer, not in canonical state.

Agent guardrail: if a future change would bypass these rules, the agent must
not proceed silently. It must ask the human author to confirm the design
exception before making the change.
