# Active Development Plans

This directory contains executable plans that still have unchecked `Step`
entries. A file belongs here because work remains, not merely because its name
contains `plan` or `checklist`.

The **only split-crate architecture migration plan and progress ledger** is
[`crate_architecture_completion_plan.md`](./crate_architecture_completion_plan.md),
with ownership defined by [`../crate_architecture.md`](../crate_architecture.md).
The domain plans below are source-port references for that migration; they must
not override its owners, naming/registration gates, or execution order.

- [`coordinate_2d_rdkit_full_port_checklist.md`](./coordinate_2d_rdkit_full_port_checklist.md)
- [`pdb_mmcif_gemmi_full_port_checklist.md`](./pdb_mmcif_gemmi_full_port_checklist.md)
- [`rdkit_assign_atom_chiral_tags_from_structure_full_port_plan.md`](./rdkit_assign_atom_chiral_tags_from_structure_full_port_plan.md)
- [`rdkit_atom_pair_fingerprint_full_port_plan.md`](./rdkit_atom_pair_fingerprint_full_port_plan.md)
- [`rdkit_mol2_full_port_plan.md`](./rdkit_mol2_full_port_plan.md)
- [`rdkit_tautomer_enumerator_full_port_plan.md`](./rdkit_tautomer_enumerator_full_port_plan.md)
- [`rdkit_topological_avalon_fingerprint_port_plan.md`](./rdkit_topological_avalon_fingerprint_port_plan.md)
- [`rdkit_topological_torsion_fingerprint_full_port_plan.md`](./rdkit_topological_torsion_fingerprint_full_port_plan.md)
- [`smiles_rdkit_full_port_checklist.md`](./smiles_rdkit_full_port_checklist.md)

Plans must follow [`../agent_plan_standard.md`](../agent_plan_standard.md).
When no executable unchecked step remains, move the plan to
[`../archive/plans/`](../archive/plans/) without rewriting its historical
completion record.
