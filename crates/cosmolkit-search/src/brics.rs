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

/// A source bond and its two (atom, environment) ends in BRICS rule-match order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BricsBondMatch {
    bond: BondId,
    ends: [(AtomId, u8); 2],
}

impl BricsBondMatch {
    pub fn bond(&self) -> BondId {
        self.bond
    }
    pub fn ends(&self) -> [(AtomId, u8); 2] {
        self.ends
    }

    /// Keep each label on its matched atom's side. Rediscover sites after editing the molecule.
    pub fn to_bond_cut(&self) -> BondCut {
        BondCut::from_matched_ends(
            self.bond,
            self.ends.map(|(atom, label)| (atom, usize::from(label))),
        )
    }
}

struct Rule {
    environments: [usize; 2],
    labels: [u8; 2],
    pattern: CompiledQuery,
}

struct Rules {
    environments: Vec<CompiledQuery>,
    bonds: Vec<Rule>,
}

fn compile_rules() -> Rules {
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
    let environments = [
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
    let groups: &[(&str, char, &[&str])] = &[
        ("L1", '-', &["L3", "L5", "L10"]),
        ("L3", '-', &["L4", "L13", "L14", "L15", "L16"]),
        ("L4", '-', &["L5", "L11"]),
        ("L5", '-', &["L12", "L14", "L16", "L13", "L15"]),
        ("L6", '-', &["L13", "L14", "L15", "L16"]),
        ("L7a", '=', &["L7b"]),
        ("L8", '-', &["L9", "L10", "L13", "L14", "L15", "L16"]),
        ("L9", '-', &["L13", "L14", "L15", "L16"]),
        ("L10", '-', &["L13", "L14", "L15", "L16"]),
        ("L11", '-', &["L13", "L14", "L15", "L16"]),
        ("L13", '-', &["L14", "L15", "L16"]),
        ("L14", '-', &["L14", "L15", "L16"]),
        ("L15", '-', &["L16"]),
        ("L16", '-', &["L16"]),
    ];
    let compile = |smarts: &str| {
        CompiledQuery::compile(
            parse_smarts(smarts, &Default::default()).expect("valid built-in BRICS SMARTS"),
        )
        .expect("valid built-in BRICS query graph")
    };
    let mut bonds = Vec::new();
    for &(left, bond, partners) in groups {
        for &right in partners {
            let indices = [left, right].map(|key| {
                environments
                    .iter()
                    .position(|(name, _)| *name == key)
                    .expect("known BRICS environment")
            });
            let [a, b] = indices.map(|i| environments[i].1);
            let labels = [left, right].map(|key| {
                key.trim_matches(|c: char| c.is_ascii_alphabetic())
                    .parse::<u8>()
                    .expect("numeric BRICS label")
            });
            bonds.push(Rule {
                environments: indices,
                labels,
                pattern: compile(&format!("[$({a})]{bond};!@[$({b})]")),
            });
        }
    }
    Rules {
        environments: environments.iter().map(|(_, s)| compile(s)).collect(),
        bonds,
    }
}

static RULES: LazyLock<Rules> = LazyLock::new(compile_rules);

/// Return eligible bonds once each. No output molecule or fragment is constructed.
pub fn find_brics_bonds(
    target: &SearchTarget<'_>,
) -> Result<Vec<BricsBondMatch>, SubstructMatchError> {
    let context = build_query_match_context(target);
    let first = SubstructMatchParams {
        max_matches: 1,
        ..Default::default()
    };
    let params = SubstructMatchParams::default();
    // RDKit✔️✔️: for env, patt in environMatchers.items():
    //     envMatches[env] = mol.HasSubstructMatch(patt)
    let mut present = Vec::with_capacity(RULES.environments.len());
    for env in &RULES.environments {
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
    for rule in &RULES.bonds {
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
                matches.push(BricsBondMatch {
                    bond,
                    ends: std::array::from_fn(|i| {
                        (AtomId::new(hit.atom_mapping[i]), rule.labels[i])
                    }),
                });
            }
        }
    }
    Ok(matches)
}
