# cosmolkit-core-old

This crate is the frozen pre-redesign `cosmolkit-core`.

It exists only as reference material while the new `crates/cosmolkit-core`
is rebuilt around value-style molecule state, explicit topology edit reports,
and invariant-checked operations.

Do not add new features here. Do not make the new core depend on this crate.
If an agent believes an exception is necessary, it must stop and confirm with
the human author before bypassing that rule.
