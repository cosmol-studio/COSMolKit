//! Terminal products from RDKit BRICS and regular r-BRICS decomposition.
//
//  Copyright (c) 2009, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//  Copyright (c) 2022, IBM LLC, Leili Zhang, Vasu Rao, Wendy Cornell
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
//     * Neither the name of International Business Machine
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

use std::collections::{HashSet, VecDeque};
use std::rc::Rc;
use std::sync::LazyLock;

use crate::decomposition::{Product, Reactant, canonical_key};
use crate::{Reaction, ReactionInput, ReactionRunError, ReactionRunParams};
use cosmolkit_core::SanitizeError;
use cosmolkit_smiles::SmilesParseError;

#[derive(Debug, thiserror::Error)]
pub enum BricsDecompositionError {
    #[error("BRICS reaction failed: {0}")]
    Reaction(#[from] ReactionRunError),
    #[error("accepted BRICS product has no SMILES because sanitization failed: {0}")]
    UnkeyedProduct(#[source] SanitizeError),
    #[error("could not write BRICS product SMILES: {0}")]
    Smiles(#[from] SmilesParseError),
    #[error("BRICS product SMILES is not valid UTF-8: {0}")]
    Encoding(#[from] std::string::FromUtf8Error),
}

static REACTIONS: LazyLock<Vec<Vec<Reaction>>> =
    LazyLock::new(|| prepare_reactions(cosmolkit_search::brics_reaction_smarts()));
static R_BRICS_REACTIONS: LazyLock<Vec<Vec<Reaction>>> =
    LazyLock::new(|| prepare_reactions(cosmolkit_search::r_brics_reaction_smarts()));

fn prepare_reactions(groups: Vec<Vec<String>>) -> Vec<Vec<Reaction>> {
    // RDKit❗✔️: reactions = tuple([[Reactions.ReactionFromSmarts(y) for y in x] for x in smartsGps])
    groups
        .into_iter()
        .map(|group| {
            group
                .into_iter()
                .map(|smarts| {
                    let mut reaction =
                        crate::parse_smirks(&smarts).expect("valid built-in BRICS template");
                    crate::runner::initialize_for_run(&mut reaction)
                        .expect("valid built-in BRICS reaction");
                    reaction
                })
                .collect()
        })
        .collect()
}

/// Enumerate unique canonical isomeric SMILES of terminal regular r-BRICS products.
pub fn decompose_r_brics(
    source: ReactionInput<'_>,
) -> Result<Vec<String>, BricsDecompositionError> {
    terminal_products(
        source,
        &R_BRICS_REACTIONS,
        FailedSanitization::RejectSequence,
    )
}

/// Enumerate unique canonical isomeric SMILES of terminal standard BRICS products.
pub fn decompose_brics(source: ReactionInput<'_>) -> Result<Vec<String>, BricsDecompositionError> {
    terminal_products(source, &REACTIONS, FailedSanitization::RetainUnkeyedProduct)
}

enum FailedSanitization {
    RejectSequence,
    RetainUnkeyedProduct,
}

fn terminal_products(
    source: ReactionInput<'_>,
    reactions: &[Vec<Reaction>],
    failed: FailedSanitization,
) -> Result<Vec<String>, BricsDecompositionError> {
    // RDKit source: rdkit/Chem/BRICS.py, BRICSDecompose (default arguments).
    // RDKit❗❗:   mSmi = Chem.MolToSmiles(mol)
    //   activePool = {mSmi: mol}
    //   allNodes.add(mSmi)
    //   foundMols = {mSmi: mol}
    //   for gpIdx, reactionGp in enumerate(reactions):
    //     newPool = {}
    //     while activePool:
    //       matched = False
    //       nSmi = next(iter(activePool))
    //       mol = activePool.pop(nSmi)
    let root = Rc::<str>::from(canonical_key::<BricsDecompositionError>(source)?);
    let mut seen = HashSet::from([Rc::clone(&root)]);
    let mut active = VecDeque::from([(root, Reactant::Source(source))]);
    let params = ReactionRunParams::default();

    for group in reactions {
        let mut survivors = VecDeque::new();
        while let Some((key, molecule)) = active.pop_front() {
            let input = molecule.as_input();
            let mut matched = false;
            for reaction in group {
                // RDKit❗❗:         ps = reaction.RunReactants((mol, ))
                for mut products in
                    crate::runner::run_reactants_source(reaction, &[input], &params)?
                {
                    // RDKit❗❗:             tSeq = [(prod.GetNumAtoms(onlyExplicit=True), idx) for idx, prod in enumerate(prodSeq)]
                    //             tSeq.sort()
                    //             for nats, idx in tSeq:
                    //               prod = prodSeq[idx]
                    //               try:
                    //                 Chem.SanitizeMol(prod)
                    //               except Exception:
                    //                 continue
                    //               pSmi = Chem.MolToSmiles(prod)
                    //               if minFragmentSize > 0:
                    //                 nDummies = pSmi.count('*')
                    //                 if nats - nDummies < minFragmentSize:
                    //                   seqOk = False
                    //                   break
                    products.sort_by_key(|p| p.topology.atoms.len());
                    let mut candidates = Vec::with_capacity(products.len());
                    let mut accepted = true;
                    for raw in products {
                        let n = raw.topology.atoms.len();
                        match Product::sanitize(raw) {
                            Ok(product) => {
                                let smiles =
                                    canonical_key::<BricsDecompositionError>(product.as_input())?;
                                if n <= smiles.bytes().filter(|&b| b == b'*').count() {
                                    accepted = false;
                                    break;
                                }
                                candidates.push(Ok((Rc::<str>::from(smiles), product)));
                            }
                            // rBRICS_public.py::rBRICSDecompose
                            // r-BRICS❗✔️:               except:
                            //                 seqOk=False #Issue 1
                            //                 continue
                            Err(_) if matches!(failed, FailedSanitization::RejectSequence) => {
                                accepted = false
                            }
                            Err(err) => candidates.push(Err(err)),
                        }
                    }
                    // RDKit❗❗:             if seqOk:
                    //               matched = True
                    //               for nats, prod in prodSeq:
                    //                 pSmi = prod.pSmi
                    //                 # print('\t',nats,pSmi)
                    //                 if pSmi not in allNodes:
                    //                   if not singlePass:
                    //                     activePool[pSmi] = prod
                    //                   allNodes.add(pSmi)
                    if accepted {
                        matched = true;
                        for candidate in candidates {
                            let (smiles, product) =
                                candidate.map_err(BricsDecompositionError::UnkeyedProduct)?;
                            if seen.insert(Rc::clone(&smiles)) {
                                active.push_back((smiles, Reactant::Product(product)));
                            }
                        }
                    }
                }
            }
            // RDKit❗❗:       if singlePass or keepNonLeafNodes or not matched:
            //         newPool[nSmi] = mol
            if !matched {
                survivors.push_back((key, molecule));
            }
        }
        active = survivors;
    }
    // RDKit❗❗:       res = set(activePool.keys())
    Ok(active.into_iter().map(|(key, _)| key.to_string()).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_core::sanitize_owned_topology;

    #[test]
    fn terminal_membership_follows_the_staged_reactions() {
        for (smi, expected) in [
            ("C", vec!["C"]),
            ("CCCOCCC", vec!["[3*]O[3*]", "[4*]CCC"]),
            ("CCOCC.[Cl-]", vec!["[3*]O[3*]", "[4*]CC"]),
            (
                "CCOC(=O)NCC",
                vec!["[1*]C([1*])=O", "[3*]O[3*]", "[4*]CC", "[5*]N[5*]"],
            ),
            (
                "CCCOCc1cccc(-c2ccccn2)c1",
                vec![
                    "[14*]c1ccccn1",
                    "[16*]c1cccc([16*])c1",
                    "[3*]O[3*]",
                    "[4*]CCC",
                    "[4*]C[8*]",
                ],
            ),
        ] {
            let record = cosmolkit_smiles::parse_smiles(smi, &Default::default()).unwrap();
            let clean = sanitize_owned_topology(record.topology, &Default::default()).unwrap();
            let input = ReactionInput {
                topology: &clean.topology,
                coordinates: &record.coordinates,
                properties: &record.properties,
                rings: clean.final_rings.as_ref(),
                valence: clean.final_valence.as_ref(),
            };
            let mut terminal = decompose_brics(input).unwrap();
            terminal.sort();
            assert_eq!(terminal, expected, "{smi}");
        }
    }
}
