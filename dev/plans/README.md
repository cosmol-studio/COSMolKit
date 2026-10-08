# Task Plans

COSMolKit 0.5.0 uses the split-crate architecture described in
[crate_architecture.md](../crate_architecture.md).
Task-specific inventories and dated receipts retain their stated scope;
they are not blanket validation of the current checkout.

Execute only the current user-authorized task and its explicitly assigned
plan. Each plan must state its task scope, dependencies, validation and
completion criteria. Its entries apply only to that authorized task.

Task plans follow [agent_plan_standard.md](../agent_plan_standard.md).
Preserve historical commands, failures and completion records in their original
scope. Completed or superseded plans belong in [archive/plans](../archive/plans/)
when an archival move is authorized. This index does not itself move files or
change any task's completion status.
