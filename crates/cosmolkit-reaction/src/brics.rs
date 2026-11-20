//! Default terminal products from RDKit Chem/BRICS.py.
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

use std::collections::{HashSet, VecDeque};
use std::rc::Rc;
use std::sync::LazyLock;

use cosmolkit_core::{RingInfo, SanitizeError, ValenceAssignment, sanitize_owned_topology};
use cosmolkit_smiles::{SmilesParseError, SmilesRecord, SmilesRecordView, write_smiles};

use crate::{Reaction, ReactionInput, ReactionProduct, ReactionRunError, ReactionRunParams};

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

static REACTIONS: LazyLock<Vec<Vec<Reaction>>> = LazyLock::new(|| {
    // RDKit❗✔️: reactions = tuple([[Reactions.ReactionFromSmarts(y) for y in x] for x in smartsGps])
    cosmolkit_search::brics_reaction_smarts()
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
});

struct Product {
    record: SmilesRecord,
    rings: Option<RingInfo>,
    valence: Option<ValenceAssignment>,
}

impl Product {
    fn sanitize(product: ReactionProduct) -> Result<Self, SanitizeError> {
        // RDKit source: Code/GraphMol/MolOps.cpp, sanitizeMol
        // RDKit❗✔️:   mol.clearComputedProps();
        let mut properties = product.properties;
        properties.clear_computed_props()?;
        let sanitized = sanitize_owned_topology(product.topology, &Default::default())?;
        Ok(Self {
            record: SmilesRecord {
                topology: sanitized.topology,
                coordinates: product.coordinates,
                properties,
            },
            rings: sanitized.final_rings,
            valence: sanitized.final_valence,
        })
    }

    fn as_input(&self) -> ReactionInput<'_> {
        ReactionInput {
            topology: &self.record.topology,
            properties: &self.record.properties,
            coordinates: &self.record.coordinates,
            rings: self.rings.as_ref(),
            valence: self.valence.as_ref(),
        }
    }
}

enum Reactant<'a> {
    Source(ReactionInput<'a>),
    Product(Product),
}

fn canonical_key(input: ReactionInput<'_>) -> Result<String, BricsDecompositionError> {
    let bytes = write_smiles(SmilesRecordView {
        topology: input.topology,
        coordinates: input.coordinates,
        properties: input.properties,
        rings: input.rings,
    })?;
    Ok(String::from_utf8(bytes.into_bytes())?)
}

/// Enumerate unique canonical isomeric SMILES of terminal standard BRICS products.
pub fn decompose_brics(source: ReactionInput<'_>) -> Result<Vec<String>, BricsDecompositionError> {
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
    let root = Rc::<str>::from(canonical_key(source)?);
    let mut seen = HashSet::from([Rc::clone(&root)]);
    let mut active = VecDeque::from([(root, Reactant::Source(source))]);
    let params = ReactionRunParams::default();

    for group in &*REACTIONS {
        let mut survivors = VecDeque::new();
        while let Some((key, molecule)) = active.pop_front() {
            let input = match &molecule {
                Reactant::Source(source) => *source,
                Reactant::Product(product) => product.as_input(),
            };
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
                                let smiles = canonical_key(product.as_input())?;
                                if n <= smiles.bytes().filter(|&b| b == b'*').count() {
                                    accepted = false;
                                    break;
                                }
                                candidates.push(Ok((Rc::<str>::from(smiles), product)));
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
