# InChI source audit record

This static audit was conducted on `audit/inchi-source`, based on
`da0cc6b7dcc4129ccdada015beee83c7f2160066` (0.5.0-rc.5).

- [Assignment](assignment.md): scope and static-only execution rules.
- [Plan](plan.md): function-level progress and next action.
- [Inventory](inventory.md): source mappings and per-function evidence.
- [Findings](findings.md): suspected differences and quarantine locations.

The assignment, plan and findings preserve the historical audit stage, which
inserted deliberate panic guards. Subsequent source repairs are recorded in
commit `3bac417`; integration retains the current public facade and bindings
instead of restoring the obsolete Python guards. These audit records are
static source evidence, not runtime test or benchmark results, and are not an
active implementation queue.
