# Boost-Compatible Double Formatting

## Approved implementation boundary

Approved on 2026-09-26 for S44-TYPED Step24. This authorizes an independent
pure-Rust implementation of the observable Double string-conversion contract;
it does not approve a behavioral difference, relaxed comparison, or an
alternative architecture. The owner remains the private formatting helper in
cosmolkit-core's property_string module, consumed by the existing detached
property conversion API.

Preserve the pinned RDKit rdvalue_tostring / LocaleSwitcher and Boost
lexical_cast precision and special-value behavior. At the system formatter
boundary, implement the public percent-g semantics using Rust's standard
numeric formatter instead of translating glibc internals. This narrow
authorization supersedes the requirement to reproduce and copy the glibc
printf_fp/MPN helper bodies for this operation only. Other source-reproduction
rules, including separate behavioral and performance review, remain binding.

Do not copy, translate, incorporate, or relicense glibc implementation code.
Do not add production libc FFI, gpoint, another formatting dependency, a new
crate, or a public formatting API. Retain legally reusable RDKit/Boost anchors
and notices at the actual implementation; describe the independent system
formatter boundary honestly. Do not claim a legally certified clean-room
process or copy glibc bodies into comments.

## Exact behavior

- Input is IEEE binary64; precision is 17 significant decimal digits.
- Match the pinned reference in C locale with round-to-nearest, ties-to-even.
  Production formatting is deterministic, not ambient-locale/fenv dependent.
  Non-default C rounding modes are outside this approved contract; do not
  claim equivalence for every floating-point environment.
- Preserve positive and negative zero and the exact Boost spellings/signs for
  infinities and NaNs. Inspect Boost get_inf_nan, not just snprintf.
- Cover every finite binary64 magnitude, including subnormal values, and
  decimal-exponent transitions. No finite-value exclusions or tolerances.
- Compare output bytes exactly, including exponent sign, minimum two exponent
  digits, decimal point, and removal of fractional trailing zeros.
- The format uses fixed notation iff the exponent X of the correctly rounded
  scientific representation satisfies -4 <= X < 17.
- Existing numeric parse acceptance and SDF field error behavior do not change.
- Boost's float path promotes to double with float-derived precision 9, but
  f32 support is not required or authorized by this Double-only repair.

## Required algorithm

For finite nonzero values, generate the magnitude's scientific representation
once with 16 fractional digits using Rust's precision-controlled formatter.
Read its rounded exponent and decimal digits; do not compute the notation
decision with log10, floating scaling, or a guessed threshold.

For fixed notation, reposition the decimal point in those already-rounded
digits and insert necessary zeros. For scientific notation, keep the rounded
mantissa and spell the exponent with an explicit sign and at least two digits.
Remove only fractional trailing zeros and an empty decimal point, never
significant integer zeros. Apply the original sign independently, including
negative zero. Handle Boost special values before finite conversion.

Do not parse the intermediate decimal string back to f64 or round it again.
Use a bounded stack formatting buffer where practical; calculate its bound
from binary64 exponent range and precision rather than swallowing overflow.
The returned String is the output allocation; avoid extra heap strings and
duplicate numeric conversions. Inspect the actual cost relative to pinned
Boost/glibc and keep performance conclusions separate from byte parity.

Rust's formatter supplies numeric rounding, not the entire Boost contract.
Record the tested Rust toolchain and source/documentation basis for that
dependency. Compiler upgrades do not automatically inherit validation claims.

## Evidence and test boundaries

Local fixed regressions belong in the core owner under property_string_double_.
They contain fixed bit patterns and expected strings; ordinary cargo tests do
not invoke an upstream oracle or generate reference data.

Retain all 35 failures from the prior gpoint 0.3.0 comparison, including signed
zero and signed power-of-ten neighbors. The diagnostic scripts and exact
deterministic selection are currently in /tmp/ck-float-compare-xYTtB3/.
Recover every case, not just the first examples printed by that script. These
are counterexamples to the candidate library, not exemptions for CK.

Cover signs, zero, ordinary decimals, ties/carry, subnormal endpoints, normal
extrema, powers-of-ten neighborhoods, notation boundaries, and NaN signs.
Reference outputs come from the pinned RDKit/Boost path, not CK.

For the current private-helper acceptance, reuse the existing out-of-repository
comparison as a one-off diagnostic against the actual owner implementation,
not a separately rewritten formatter. Do not commit a new corpus runner or
add an owner-crate corpus loop. No production/test-only public export is
authorized. A temporary external harness may include the actual owner source
and its existing detached dependencies without changing production visibility.
Keep harness/source hashes, seed, full input identities, counts and mismatches
in ignored evidence; this is not a registered public parity capability.

The reproducible diagnostic selection is the previous 20,814-case comparison
plus 1,000,000 random binary64 patterns using its recorded generator/seed.
Add signed 32-predecessor/32-successor neighborhoods around each finite
nonzero decimal power from 10^-323 through 10^308, retaining overlap and signed
zeros with explicit case IDs. Include signaling/quiet NaNs and payload/sign
variants. Full preparation and input/reference identity checks precede CK
comparison; no dropped failures or zero-case passes.

Use the installed RDKit 2026.03.6 / Boost 1.85 / verified glibc
2.43-2ubuntu2.4 reference, C locale and verified default rounding. Actual RDKit
Double-property string conversion is an acceptable oracle because it executes
the relevant Boost path. An external pinned C++ Boost oracle is also allowed
if built from already available sources without dependency/system changes.
A bare snprintf oracle does not replace Boost special-value verification.

Permanent corpus parity remains exclusively in top-level parity-tests through
public cosmolkit with full. Do not expose a helper merely to register a task.
Public workflow coverage is separate from this private implementation gate.
Cross-platform/WASM results must be reported separately; native success is not
a WASM execution claim. Random agreement is evidence, not an exhaustive proof.

## Acceptance and continuation

Implementation-method approval removes the need to request LGPL incorporation
for this route. It does not complete Step24 or grant permission to label
untested behavior parity. Require complete owner implementation, nonzero fixed
regression passes, zero byte mismatches in the diagnostic selection, and
source/algorithm/performance review before closing the refinement.

Resume S44-TYPED Steps25-50, S44-PROPS reclosure and the original SMILES queue
after this gate. Existing strict/core/runtime/frozen-runner requirements remain.
No Git, dependency changes, cross-workspace edits by the worker, broad policy
changes or support-state promotions are authorized.

## References

- Pinned Boost lexical_cast: 02e5821ab32c45fad719829e9644e5d681c9ba0b,
  converter_lexical_streams.hpp, lcast_precision.hpp and inf_nan.hpp.
- Pinned Boost core: 083b41c17e34f1fc9b43ab796b40d0d8bece685c,
  include/boost/core/snprintf.hpp.
- C11 draft N1570, section 7.21.6.1, g/G conversion rules.
- Rust std::fmt precision and localization:
  https://doc.rust-lang.org/std/fmt/index.html#precision
- Existing S44-BOOST receipt records reference loader/source provenance;
  it remains historical evidence, not an implementation license.

