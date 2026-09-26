# InChI source audit

This audit is isolated on `audit/inchi-source`, based on
`da0cc6b7dcc4129ccdada015beee83c7f2160066` (0.5.0-rc.5).

- [Assignment](assignment.md): scope and static-only execution rules.
- [Plan](plan.md): function-level progress and next action.
- [Inventory](inventory.md): source mappings and per-function evidence.
- [Findings](findings.md): suspected differences and quarantine locations.

The branch deliberately contains audit panics. It is not a release candidate
and must not be merged or published as production code. Findings are static
evidence, not runtime test results. Builds, tests and benchmarks are not run.
Reports may be committed on this audit branch; an agent may not commit or push
without explicit authorization. Moving the audit does not reset its progress
or certify existing findings.
