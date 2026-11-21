//! Ordered sites from the default RDKit Chem/BRICS.py FindBRICSBonds rules.
//
//  Copyright (c) 2009, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

use crate::matcher::{
    full_matches_with_compiled_query_and_context,
    try_get_substruct_atom_matches_with_compiled_query_and_context,
};
use crate::{
    CompiledQuery, SearchTarget, SubstructMatchError, SubstructMatchParams,
    build_query_match_context, parse_smarts,
};
use cosmolkit_core::BondCut;
use cosmolkit_model::{AtomId, BondId};
use std::{collections::HashSet, sync::LazyLock};

mod regular;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct OrientedSite {
    bond: BondId,
    ends: [(AtomId, u8); 2],
}

macro_rules! bond_match {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq)]
        pub struct $name(OrientedSite);
        impl $name {
            pub fn bond(&self) -> BondId {
                self.0.bond
            }
            pub fn ends(&self) -> [(AtomId, u8); 2] {
                self.0.ends
            }

            /// Keep each label on its matched atom's side. Rediscover sites after editing the molecule.
            pub fn to_bond_cut(&self) -> BondCut {
                BondCut::from_matched_ends(
                    self.0.bond,
                    self.0.ends.map(|(atom, label)| (atom, usize::from(label))),
                )
            }
        }
    };
}
bond_match!(
    /// A bond with oriented standard BRICS environment labels.
    BricsBondMatch
);
bond_match!(
    /// A bond with oriented regular r-BRICS environment labels.
    RBricsBondMatch
);

struct Rule {
    environments: [usize; 2],
    labels: [u8; 2],
    pattern: CompiledQuery,
}

struct Rules {
    environments: Vec<CompiledQuery>,
    bonds: Vec<Rule>,
}

const ENVIRONMENTS: &[(&str, &str)] = &[
    ("L1", "[C;D3]([#0,#6,#7,#8])(=O)"),
    ("L3", "[O;D2]-;!@[#0,#6,#1]"),
    ("L4", "[C;!D1;!$(C=*)]-;!@[#6]"),
    (
        "L5",
        "[N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]",
    ),
    ("L6", "[C;D3;!R](=O)-;!@[#0,#6,#7,#8]"),
    ("L7a", "[C;D2,D3]-[#6]"),
    ("L7b", "[C;D2,D3]-[#6]"),
    ("#L8", "[C;!R;!D1]-;!@[#6]"),
    ("L8", "[C;!R;!D1;!$(C!-*)]"),
    ("L9", "[n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]"),
    ("L10", "[N;R;$(N(@C(=O))@[C,N,O,S])]"),
    ("L11", "[S;D2](-;!@[#0,#6])"),
    ("L12", "[S;D4]([#6,#0])(=O)(=O)"),
    ("L13", "[C;$(C(-;@[C,N,O,S])-;@[N,O,S])]"),
    ("L14", "[c;$(c(:[c,n,o,s]):[n,o,s])]"),
    ("L14b", "[c;$(c(:[c,n,o,s]):[n,o,s])]"),
    ("L15", "[C;$(C(-;@C)-;@C)]"),
    ("L16", "[c;$(c(:c):c)]"),
    ("L16b", "[c;$(c(:c):c)]"),
];
struct RuleSet {
    environments: &'static [(&'static str, &'static str)],
    groups: &'static [&'static [(&'static str, &'static str, &'static str)]],
}

const STANDARD: RuleSet = RuleSet {
    environments: ENVIRONMENTS,
    groups: &[
        &[("L1", "L3", "-"), ("L1", "L5", "-"), ("L1", "L10", "-")],
        &[
            ("L3", "L4", "-"),
            ("L3", "L13", "-"),
            ("L3", "L14", "-"),
            ("L3", "L15", "-"),
            ("L3", "L16", "-"),
        ],
        &[("L4", "L5", "-"), ("L4", "L11", "-")],
        &[
            ("L5", "L12", "-"),
            ("L5", "L14", "-"),
            ("L5", "L16", "-"),
            ("L5", "L13", "-"),
            ("L5", "L15", "-"),
        ],
        &[
            ("L6", "L13", "-"),
            ("L6", "L14", "-"),
            ("L6", "L15", "-"),
            ("L6", "L16", "-"),
        ],
        &[("L7a", "L7b", "=")],
        &[
            ("L8", "L9", "-"),
            ("L8", "L10", "-"),
            ("L8", "L13", "-"),
            ("L8", "L14", "-"),
            ("L8", "L15", "-"),
            ("L8", "L16", "-"),
        ],
        &[
            ("L9", "L13", "-"),
            ("L9", "L14", "-"),
            ("L9", "L15", "-"),
            ("L9", "L16", "-"),
        ],
        &[
            ("L10", "L13", "-"),
            ("L10", "L14", "-"),
            ("L10", "L15", "-"),
            ("L10", "L16", "-"),
        ],
        &[
            ("L11", "L13", "-"),
            ("L11", "L14", "-"),
            ("L11", "L15", "-"),
            ("L11", "L16", "-"),
        ],
        &[
            ("L13", "L14", "-"),
            ("L13", "L15", "-"),
            ("L13", "L16", "-"),
        ],
        &[
            ("L14", "L14", "-"),
            ("L14", "L15", "-"),
            ("L14", "L16", "-"),
        ],
        &[("L15", "L16", "-")],
        &[("L16", "L16", "-")],
    ],
};

struct RuleDefinition {
    environments: [usize; 2],
    labels: [u8; 2],
    bond: &'static str,
}

fn definitions(
    set: &'static RuleSet,
) -> impl Iterator<Item = impl Iterator<Item = RuleDefinition>> {
    set.groups.iter().map(move |group| {
        group
            .iter()
            .map(move |&(left, right, bond)| RuleDefinition {
                environments: [left, right].map(|key| {
                    set.environments
                        .iter()
                        .position(|(name, _)| *name == key)
                        .expect("known BRICS environment")
                }),
                labels: [left, right].map(|key| {
                    key.trim_matches(|c: char| c.is_ascii_alphabetic())
                        .parse()
                        .expect("numeric BRICS label")
                }),
                // rBRICS_public.py, init_reactions: if "@" not in bnd: #Leili
                bond: match bond {
                    "-" => "-;!@",
                    "=" => "=;!@",
                    explicit => explicit,
                },
            })
    })
}

/// Standard BRICS reaction templates in native reaction-group order.
#[doc(hidden)]
pub fn brics_reaction_smarts() -> Vec<Vec<String>> {
    // RDKit source: rdkit/Chem/BRICS.py, smartsGps
    // RDKit❗✔️:     sma = '[$(%s):1]%s;!@[$(%s):2]>>[%s*]-[*:1].[%s*]-[*:2]' % (r1, bnd, r2, g1, g2)
    reaction_smarts(&STANDARD)
}

#[doc(hidden)]
pub fn r_brics_reaction_smarts() -> Vec<Vec<String>> {
    reaction_smarts(&regular::REGULAR)
}

fn reaction_smarts(set: &'static RuleSet) -> Vec<Vec<String>> {
    // r-BRICS❗✔️:         sma = '[$(%s):1]%s[$(%s):2]>>[%s*]-[*:1].[%s*]-[*:2]' % (r1, bnd, r2, g1, g2)
    definitions(set)
        .map(|group| {
            group
                .map(|rule| {
                    let [a, b] = rule.environments.map(|i| set.environments[i].1);
                    let [x, y] = rule.labels;
                    let bond = rule.bond;
                    format!("[$({a}):1]{bond}[$({b}):2]>>[{x}*]-[*:1].[{y}*]-[*:2]")
                })
                .collect()
        })
        .collect()
}

fn compile_rules(set: &'static RuleSet) -> Rules {
    // RDKit✔️✔️: for env, sma in environs.items():
    //   environMatchers[env] = Chem.MolFromSmarts(sma)
    // RDKit✔️✔️: for i, compats in enumerate(reactionDefs):
    //   tmp = []
    //   for i1, i2, bType in compats:
    //     e1 = environs['L%s' % i1]
    //     e2 = environs['L%s' % i2]
    //     patt = '[$(%s)]%s;!@[$(%s)]' % (e1, bType, e2)
    //     patt = Chem.MolFromSmarts(patt)
    //     tmp.append((i1, i2, bType, patt))
    //   bondMatchers.append(tmp)
    let compile = |smarts: &str| {
        CompiledQuery::compile(
            parse_smarts(smarts, &Default::default()).expect("valid built-in BRICS SMARTS"),
        )
        .expect("valid built-in BRICS query graph")
    };
    let bonds = definitions(set)
        .flatten()
        .map(|rule| {
            let [a, b] = rule.environments.map(|i| set.environments[i].1);
            let bond = rule.bond;
            Rule {
                environments: rule.environments,
                labels: rule.labels,
                pattern: compile(&format!("[$({a})]{bond}[$({b})]")),
            }
        })
        .collect();
    Rules {
        environments: set.environments.iter().map(|(_, s)| compile(s)).collect(),
        bonds,
    }
}

static RULES: LazyLock<Rules> = LazyLock::new(|| compile_rules(&STANDARD));
static R_BRICS_RULES: LazyLock<Rules> = LazyLock::new(|| compile_rules(&regular::REGULAR));

/// Return eligible bonds once each. No output molecule or fragment is constructed.
pub fn find_brics_bonds(
    target: &SearchTarget<'_>,
) -> Result<Vec<BricsBondMatch>, SubstructMatchError> {
    find_sites(target, &RULES, BricsBondMatch)
}

/// Find regular r-BRICS matches in upstream rule order.
pub fn find_r_brics_bonds(
    target: &SearchTarget<'_>,
) -> Result<Vec<RBricsBondMatch>, SubstructMatchError> {
    find_sites(target, &R_BRICS_RULES, RBricsBondMatch)
}

fn find_sites<T>(
    target: &SearchTarget<'_>,
    rules: &Rules,
    wrap: fn(OrientedSite) -> T,
) -> Result<Vec<T>, SubstructMatchError> {
    // r-BRICS❗✔️:   for env,patt in environMatchers.items(): #Leili
    //     envMatches[env]=mol.HasSubstructMatch(patt)
    let context = build_query_match_context(target);
    let first = SubstructMatchParams {
        max_matches: 1,
        ..Default::default()
    };
    let params = SubstructMatchParams::default();
    // RDKit✔️✔️: for env, patt in environMatchers.items():
    //     envMatches[env] = mol.HasSubstructMatch(patt)
    let mut present = Vec::with_capacity(rules.environments.len());
    for env in &rules.environments {
        present.push(
            !try_get_substruct_atom_matches_with_compiled_query_and_context(
                target,
                &env.query,
                &first,
                &env.compiled_graph,
                &context,
            )?
            .is_empty(),
        );
    }

    // RDKit✔️✔️: for i1, i2, bType, patt in compats:
    //       if not envMatches['L' + i1] or not envMatches['L' + i2]:
    //         continue
    //       matches = mol.GetSubstructMatches(patt)
    //       i1 = letter.sub('', i1)
    //       i2 = letter.sub('', i2)
    //       for match in matches:
    //         match = tuple(match)
    //         if match not in bondsDone and (match[1], match[0]) not in bondsDone:
    //           bondsDone.add(match)
    //           yield (((match[0], match[1]), (i1, i2)))
    let mut matches = Vec::new();
    let mut used = HashSet::new();
    for rule in &rules.bonds {
        if !rule.environments.iter().all(|&i| present[i]) {
            continue;
        }
        let pattern = &rule.pattern;
        for hit in full_matches_with_compiled_query_and_context(
            target,
            &pattern.query,
            &params,
            &pattern.compiled_graph,
            &context,
        )? {
            let bond = BondId::new(hit.bond_mapping[0]);
            if used.insert(bond) {
                matches.push(wrap(OrientedSite {
                    bond,
                    ends: std::array::from_fn(|i| {
                        (AtomId::new(hit.atom_mapping[i]), rule.labels[i])
                    }),
                }));
            }
        }
    }
    Ok(matches)
}
