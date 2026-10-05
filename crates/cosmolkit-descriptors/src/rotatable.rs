//! Rotatable-bond descriptor owner (RDKit `Lipinski.cpp`).

use crate::{DescriptorError, DescriptorInput, DescriptorResult};

/// Selects which rotatable-bond definition the count functions use.
///
/// Mirrors RDKit `NumRotatableBondsOptions` (`Lipinski.h:39-43`). The C++
/// int discriminants (`Default = -1`, `NonStrict = 0`, `Strict = 1`,
/// `StrictLinkages = 2`) are never observable through the source API (the
/// option is only compared for equality), so the Rust projection carries no
/// numeric discriminants.
///
/// The source also exposes a deprecated `bool` overload
/// (`Lipinski.cpp:185-187`) mapping `strict == true` to `Strict` and
/// `false` to `NonStrict`; per the packet freeze it is recorded here as
/// source documentation only and is deliberately not projected as a public
/// compatibility alias.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RotatableBondsOptions {
    /// The unspecified option; resolves to the pinned build default.
    ///
    /// RDKit resolves `Default` through `DefaultStrictDefinition`
    /// (`Lipinski.cpp:99-106`). The pinned RDKit build defines
    /// `RDK_USE_STRICT_ROTOR_DEFINITION` (`third_party/rdkit/CMakeLists.txt`
    /// line 55, option default `ON`), so `Default` resolves to
    /// [`RotatableBondsOptions::Strict`].
    Default,
    /// Original (loose) SMARTS definition.
    NonStrict,
    /// Stricter SMARTS definition excluding amides, esters, etc.
    Strict,
    /// Much stricter ring-linkage-aware arithmetic definition.
    StrictLinkages,
}

impl Default for RotatableBondsOptions {
    /// Source default argument is `useStrictDefinition = Default`
    /// (`Lipinski.h:57-58`); the Rust projection expresses the same
    /// unspecified option through the [`Default`] trait.
    fn default() -> Self {
        Self::Default
    }
}

/// Version of the rotatable-bond descriptor family (`Lipinski.cpp:111`).
pub const NUM_ROTATABLE_BONDS_VERSION: &str = "3.2.0";

/// Resolves the unspecified option to the pinned build default.
///
/// The source keeps `DefaultStrictDefinition` in an anonymous namespace
/// (`Lipinski.cpp:99-106`), so this helper is crate-private, exactly like
/// the source constant it projects.
pub(crate) fn resolve_options(option: RotatableBondsOptions) -> RotatableBondsOptions {
    // RDKit source (Lipinski.cpp:99-106, anonymous namespace):
    //   #ifdef RDK_USE_STRICT_ROTOR_DEFINITION
    //   const NumRotatableBondsOptions DefaultStrictDefinition = Strict;
    //   #else
    //   const NumRotatableBondsOptions DefaultStrictDefinition = NonStrict;
    //   #endif
    //
    // The pinned RDKit build turns RDK_USE_STRICT_ROTOR_DEFINITION ON
    // (third_party/rdkit/CMakeLists.txt:55 and 586-588 inject the define
    // when the option is ON), so DefaultStrictDefinition is Strict in this
    // projection.
    //
    // RDKit✔️✔️: if (strict == Default) {
    // RDKit✔️✔️:   strict = DefaultStrictDefinition;
    // RDKit✔️✔️: }
    // (Lipinski.cpp:113-115)
    //
    // Zero-allocation match, same single comparison-and-substitute shape as
    // the source; NonStrict/Strict/StrictLinkages pass through unchanged.
    match option {
        RotatableBondsOptions::Default => RotatableBondsOptions::Strict,
        resolved => resolved,
    }
}

/// Fixed SMARTS for the NonStrict rotatable-bond definition
/// (`Lipinski.cpp:116`). Two-atom query joined by the (single or
/// aromatic) non-ring bond primitive; each endpoint recursively excludes
/// any atom attached to a triple bond and any degree-1 atom. Public
/// because the pattern itself is frozen source data asserted by the
/// regressions.
pub const NON_STRICT_ROTATABLE_PATTERN: &str = "[!$(*#*)&!D1]-,:;!@[!$(*#*)&!D1]";

/// Fixed SMARTS for the Strict rotatable-bond definition
/// (`Lipinski.cpp:119-126`). The C++ source assembles this ONE string from
/// five adjacent string literals; both endpoints exclude triple-adjacent,
/// degree-1, CF3/CCl3/CBr3, tert-butyl, and methyl atoms, and the FIRST
/// endpoint additionally excludes amide/ester/amidinium carbonyl carbons
/// and their N/O/S partners. Public because the pattern itself is frozen
/// source data asserted by the regressions.
pub const STRICT_ROTATABLE_PATTERN: &str = "[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])&!$([CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-,:;!@[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])&!$([CH3])]";

/// Fixed SMARTS for the StrictLinkages base rotatable-bond query
/// (`Lipinski.cpp:144`). An endpoint is excluded only when it is
/// degree-1 AND not hydrogen — an explicit-H D1 endpoint does NOT block,
/// and triple-bond adjacency is NOT excluded here (unlike the NonStrict
/// pattern). Public because the pattern itself is frozen source data
/// asserted by the regressions.
pub const STRICT_LINKAGES_BASE_PATTERN: &str = "[!$([D1&!#1])]-,:;!@[!$([D1&!#1])]";

/// Fixed SMARTS for the StrictLinkages symmetric-ring correction
/// (`Lipinski.cpp:146-149`; the C++ source assembles ONE string from two
/// adjacent literals). Matches a chain single/aromatic bond between two
/// aromatic 6-ring ipso atoms whose recursive queries each require a
/// chain-bonded aromatic r6 partner plus two aromatic non-H ring
/// neighbors. Public because the pattern is frozen source data asserted
/// by the regressions.
pub const SYMMETRIC_RINGS_PATTERN: &str =
    "[a;r6;$(a(-,:;!@[a;r6])(a[!#1])a[!#1])]-,:;!@[a;r6;$(a(-,:;!@[a;r6])(a[!#1])a)]";

/// Fixed SMARTS for the StrictLinkages triple-bond correction
/// (`Lipinski.cpp:150`, used at lines 162-167). Aliphatic carbon,
/// triple bond, second atom carbon OR nitrogen. Despite the source
/// variable name there is NO degree constraint: any C#C or C#N triple
/// bond matches. Needed because the StrictLinkages base pattern —
/// unlike the NonStrict pattern — does not exclude triple-adjacent
/// atoms. Public because the pattern is frozen source data asserted by
/// the regressions.
pub const TERMINAL_TRIPLE_BONDS_PATTERN: &str = "C#[#6,#7]";

/// Fixed SMARTS for the StrictLinkages shared-atom amide correction
/// (`Lipinski.cpp:145`, consumed by the ordered-match loop at lines
/// 168-183). Acyclic carbonyl carbon with an `=O` branch, a
/// single-or-aromatic-bonded N, and that N single-or-aromatic-bonded to a
/// further C — i.e. N-substituted non-ring amide `C(=O)N-C` fragments.
/// Unlike the Q05 amide-count pattern `C(=[O;!R])N`, the `!R` sits on the
/// CARBON and the trailing N-C bond is REQUIRED. Public because the
/// pattern is frozen source data asserted by the regressions.
pub const NON_RING_AMIDES_PATTERN: &str = "[C&!R](=O)NC";

/// StrictLinkages staged arithmetic (complete chain R04+R05+R06: base
/// count, early return, symmetric-ring subtraction, triple-bond
/// subtraction, ordered shared-atom amide correction, final clamp).
///
/// Real caller stage for the StrictLinkages else-branch of
/// `calcNumRotatableBonds` (`Lipinski.cpp:129-189`); the public
/// [`RotatableBondsOptions::StrictLinkages`] arm delegates here.
///
/// Behavior review (R04): reproduces lines 144-163 — count the base
/// pattern; if zero, return 0 WITHOUT evaluating the symmetric-ring
/// pattern (source early return); otherwise subtract the symmetric-ring
/// count and clamp at zero.
///
/// Behavior review (R05): reproduces lines 162-167 — subtract the
/// terminal-triple-bond count (any C#C or C#N; the base pattern does
/// not exclude triple-adjacent atoms) and clamp at zero.
///
/// Behavior review (R06): reproduces lines 168-189 — evaluate the
/// non-ring-amide pattern ONCE into ordered uniquified matches, then walk
/// them in order with an atoms-seen bitset: a match is distinct only if
/// NONE of its target atoms was marked by an earlier match; EVERY atom of
/// EVERY match is marked even after an overlap or once `res` reaches zero
/// (no break, no reordering, no deduplication, no fragment counting, no
/// match-count subtraction); decrement only when distinct AND `res > 0`;
/// final defensive clamp; the unsigned cast happens at the public arm.
///
/// Complexity review: ONE prepared shared query context is built per
/// whole evaluation and reused by all four pattern evaluations through
/// the `_with_context` path (no per-pattern preparation, no valence/ring
/// recomputation) — matching the source's reuse of ROMol-internal state
/// across its matchers. Four matcher passes in source order with the
/// source's early exit; the atoms-seen set is a packed `Vec<u64>` of
/// `atom_count.div_ceil(64)` words with constant-time test-and-mark
/// (the `boost::dynamic_bitset` projection; NOT a byte-per-bool
/// `Vec<bool>`); per-match work is O(matched atoms). 🔝
/// retained-flyweight family deviation unchanged.
pub(crate) fn strict_linkages_stage(input: &DescriptorInput<'_>) -> DescriptorResult<i32> {
    // ONE prepared context for the whole StrictLinkages evaluation; the
    // four pattern evaluations below all reuse it (supervisory correction:
    // no per-pattern preparation or valence/ring recomputation).
    let context = crate::patterns::prepared_context(input, "num_rotatable_bonds")?;

    // RDKit source (Lipinski.cpp:144,146-149,154-158,160-163):
    //   pattern_flyweight rotBonds_matcher("[!$([D1&!#1])]-,:;!@[!$([D1&!#1])]");
    //   pattern_flyweight symRings_matcher(
    //       "[a;r6;$(a(-,:;!@[a;r6])(a[!#1])a[!#1])]-,:;!@[a;r6;$(a(-,:;!@[a;r6])("
    //       "a[!#1])a)]");
    //   int res = rotBonds_matcher.get().countMatches(mol);
    //   if (!res) {
    //     return 0;
    //   }
    //   // remove symmetrical rings:
    //   res -= symRings_matcher.get().countMatches(mol);
    //   if (res < 0) {
    //     res = 0;
    //   }
    // RDKit✔️✔️: int res = rotBonds_matcher.get().countMatches(mol);
    // RDKit✔️✔️: if (!res) {
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: // remove symmetrical rings:
    // RDKit✔️✔️: res -= symRings_matcher.get().countMatches(mol);
    // RDKit✔️✔️: if (res < 0) {
    // RDKit✔️✔️:   res = 0;
    // RDKit✔️✔️: }
    let base = crate::patterns::count_pattern_matches_with_context(
        input,
        "num_rotatable_bonds",
        STRICT_LINKAGES_BASE_PATTERN,
        &context,
    )? as i32;
    if base == 0 {
        return Ok(0);
    }
    let symmetric_rings = crate::patterns::count_pattern_matches_with_context(
        input,
        "num_rotatable_bonds",
        SYMMETRIC_RINGS_PATTERN,
        &context,
    )? as i32;
    let mut res = (base - symmetric_rings).max(0);

    // RDKit source (Lipinski.cpp:162-167):
    //   // remove triple bonds
    //   res -= terminalTripleBonds_matcher.get().countMatches(mol);
    //   if (res < 0) {
    //     res = 0;
    //   }
    // RDKit✔️✔️: // remove triple bonds
    // RDKit✔️✔️: res -= terminalTripleBonds_matcher.get().countMatches(mol);
    // RDKit✔️✔️: if (res < 0) {
    // RDKit✔️✔️:   res = 0;
    // RDKit✔️✔️: }
    let terminal_triples = crate::patterns::count_pattern_matches_with_context(
        input,
        "num_rotatable_bonds",
        TERMINAL_TRIPLE_BONDS_PATTERN,
        &context,
    )? as i32;
    res = (res - terminal_triples).max(0);

    // RDKit source (Lipinski.cpp:168-189):
    //   // removing amides is more complex
    //   boost::dynamic_bitset<> atomsSeen(mol.getNumAtoms());
    //   SubstructMatch(mol, *(nonRingAmides_matcher.get().getMatcher()), matches);
    //   for (const auto &iv : matches) {
    //     bool distinct = true;
    //     for (const auto &mIt : iv) {
    //       if (atomsSeen[mIt.second]) {
    //         distinct = false;
    //       }
    //       atomsSeen.set(mIt.second);
    //     }
    //     if (distinct && res > 0) {
    //       --res;
    //     }
    //   }
    //
    //   if (res < 0) {
    //     res = 0;
    //   }
    //   return static_cast<unsigned int>(res);
    // RDKit✔️✔️: // removing amides is more complex
    // RDKit✔️✔️: boost::dynamic_bitset<> atomsSeen(mol.getNumAtoms());
    // RDKit✔️✔️: SubstructMatch(mol, *(nonRingAmides_matcher.get().getMatcher()), matches);
    // RDKit✔️✔️: for (const auto &iv : matches) {
    // RDKit✔️✔️:   bool distinct = true;
    // RDKit✔️✔️:   for (const auto &mIt : iv) {
    // RDKit✔️✔️:     if (atomsSeen[mIt.second]) {
    // RDKit✔️✔️:       distinct = false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     atomsSeen.set(mIt.second);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (distinct && res > 0) {
    // RDKit✔️✔️:     --res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: if (res < 0) {
    // RDKit✔️✔️:   res = 0;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: return static_cast<unsigned int>(res);
    //
    // Rust projection of `boost::dynamic_bitset<>`: a packed Vec<u64> of
    // atom_count.div_ceil(64) words with constant-time test/set (a
    // byte-per-bool Vec<bool> is deliberately NOT used). `atom_mapping`
    // holds target atom indices indexed by query atom — exactly the
    // `mIt.second` projection. The loop keeps the source's exact shape:
    // test-then-mark EVERY atom of EVERY match in matcher order (no break
    // on overlap or res==0, no reorder/dedup), decrement only when
    // distinct AND res > 0. The unsigned cast is at the public arm.
    let amide_matches = crate::patterns::pattern_matches_with_context(
        input,
        "num_rotatable_bonds",
        NON_RING_AMIDES_PATTERN,
        &context,
    )?;
    let atom_count = input.topology().atoms.len();
    let mut atoms_seen = vec![0u64; atom_count.div_ceil(64)];
    for amide_match in &amide_matches {
        let mut distinct = true;
        for &atom in &amide_match.atom_mapping {
            let (word, bit) = (atom / 64, atom % 64);
            if atoms_seen[word] & (1u64 << bit) != 0 {
                distinct = false;
            }
            atoms_seen[word] |= 1u64 << bit;
        }
        if distinct && res > 0 {
            res -= 1;
        }
    }
    Ok(res.max(0))
}

/// Rotatable-bond count over prepared FINAL input.
///
/// Behavior review: reproduces the full `calcNumRotatableBonds` dispatch
/// (`Lipinski.cpp:113-189`) — resolve `Default` through the pinned
/// `DefaultStrictDefinition`; `NonStrict` and `Strict` count their fixed
/// SMARTS with default `SubstructMatchParameters` (uniquify);
/// `StrictLinkages` runs the complete staged arithmetic in
/// [`strict_linkages_stage`] (base, symmetric rings, triple bonds, ordered
/// shared-atom amide correction, clamps) and casts the final non-negative
/// value to `u32` (the source's `static_cast<unsigned int>`).
///
/// Complexity review: NonStrict/Strict take one retained-query acquisition
/// plus one matcher pass; StrictLinkages builds ONE prepared shared query
/// context per evaluation and reuses it for all four pattern passes (see
/// [`strict_linkages_stage`]). Same cost class as the source flyweight +
/// `SubstructMatch` with the per-call pattern construction removed by
/// construction (same recorded 🔝 deviation family as the `Q01` pattern
/// flyweight).
pub fn num_rotatable_bonds_prepared(
    input: &DescriptorInput<'_>,
    options: RotatableBondsOptions,
) -> DescriptorResult<u32> {
    // RDKit source (Lipinski.cpp:113-118):
    //   if (strict == Default) {
    //     strict = DefaultStrictDefinition;
    //   }
    //
    //   if (strict == NonStrict) {
    //     std::string pattern = "[!$(*#*)&!D1]-,:;!@[!$(*#*)&!D1]";
    //     pattern_flyweight m(pattern);
    //     return m.get().countMatches(mol);
    //   }
    // RDKit✔️✔️: if (strict == Default) {
    // RDKit✔️✔️:   strict = DefaultStrictDefinition;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (strict == NonStrict) {
    // RDKit✔️✔️:   std::string pattern = "[!$(*#*)&!D1]-,:;!@[!$(*#*)&!D1]";
    // RDKit✔️✔️:   pattern_flyweight m(pattern);
    // RDKit✔️✔️:   return m.get().countMatches(mol);
    // RDKit✔️✔️: }
    match resolve_options(options) {
        RotatableBondsOptions::NonStrict => crate::patterns::count_pattern_matches(
            input,
            "num_rotatable_bonds",
            NON_STRICT_ROTATABLE_PATTERN,
        ),
        // `resolve_options` already maps Default -> Strict (the pinned
        // DefaultStrictDefinition); the or-pattern keeps this match
        // exhaustive over the unresolved enum and evaluates Default exactly
        // like its resolved value.
        RotatableBondsOptions::Default | RotatableBondsOptions::Strict => {
            // RDKit source (Lipinski.cpp:118-128):
            //   } else if (strict == Strict) {
            //     std::string strict_pattern =
            //         "[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])("
            //         "[CH3])[CH3])&!$([CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]="
            //         "[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-,:;!@[!$"
            //         "(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])(["
            //         "CH3])[CH3])&!$([CH3])]";
            //     pattern_flyweight m(strict_pattern);
            //     return m.get().countMatches(mol);
            //   }
            // RDKit✔️✔️: } else if (strict == Strict) {
            // RDKit✔️✔️:   std::string strict_pattern =
            // RDKit✔️✔️:       "[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])("
            // RDKit✔️✔️:       "[CH3])[CH3])&!$([CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]="
            // RDKit✔️✔️:       "[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]-,:;!@[!$"
            // RDKit✔️✔️:       "(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])(["
            // RDKit✔️✔️:       "CH3])[CH3])&!$([CH3])]";
            // RDKit✔️✔️:   pattern_flyweight m(strict_pattern);
            // RDKit✔️✔️:   return m.get().countMatches(mol);
            // RDKit✔️✔️: }
            //
            // (C++ adjacent-string-literal concatenation resolves to the
            // single STRICT_ROTATABLE_PATTERN string above.)
            crate::patterns::count_pattern_matches(
                input,
                "num_rotatable_bonds",
                STRICT_ROTATABLE_PATTERN,
            )
        }
        RotatableBondsOptions::StrictLinkages => {
            // RDKit source (Lipinski.cpp:129,168-189 tail):
            //   } else {
            //     ...staged arithmetic anchored in strict_linkages_stage...
            //     if (res < 0) {
            //       res = 0;
            //     }
            //     return static_cast<unsigned int>(res);
            //   }
            // RDKit✔️✔️: } else {
            // RDKit✔️✔️:   if (res < 0) {
            // RDKit✔️✔️:     res = 0;
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   return static_cast<unsigned int>(res);
            // RDKit✔️✔️: }
            //
            // The staged function already clamps at every subtraction; the
            // final clamp + unsigned cast are reproduced here (the source
            // cast cannot produce a negative input to try_from after the
            // clamps, so the CountOverflow mapping only covers the u32
            // domain edge).
            let res = strict_linkages_stage(input)?;
            u32::try_from(res.max(0)).map_err(|_| DescriptorError::CountOverflow {
                function: "num_rotatable_bonds",
                field: "strict_linkages_result",
            })
        }
    }
}

#[cfg(test)]
mod original_condition_query_tests {
    use super::*;
    #[test]
    fn fixed_rotatable_bond_queries_are_singletons_across_repeated_and_parallel_reads() {
        for pattern in [
            NON_STRICT_ROTATABLE_PATTERN,
            STRICT_ROTATABLE_PATTERN,
            STRICT_LINKAGES_BASE_PATTERN,
            NON_RING_AMIDES_PATTERN,
            SYMMETRIC_RINGS_PATTERN,
            TERMINAL_TRIPLE_BONDS_PATTERN,
        ] {
            let first = crate::patterns::retained_pattern("num_rotatable_bonds", pattern).unwrap();
            let second = crate::patterns::retained_pattern("num_rotatable_bonds", pattern).unwrap();
            assert!(std::sync::Arc::ptr_eq(&first, &second));
            let address = std::sync::Arc::as_ptr(&first) as usize;
            std::thread::scope(|scope| {
                let handles = (0..8)
                    .map(|_| {
                        scope.spawn(move || {
                            std::sync::Arc::as_ptr(
                                &crate::patterns::retained_pattern("num_rotatable_bonds", pattern)
                                    .unwrap(),
                            ) as usize
                        })
                    })
                    .collect::<Vec<_>>();
                for handle in handles {
                    assert_eq!(handle.join().unwrap(), address);
                }
            });
        }
    }
}
