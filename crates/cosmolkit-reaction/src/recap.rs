//! RDKit RECAP decomposition with a per-molecule product hierarchy.
//
//  Copyright (c) 2007, Novartis Institutes for BioMedical Research Inc.
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

use crate::decomposition::{Product, Reactant, canonical_key};
use crate::{Reaction, ReactionInput, ReactionRunError, ReactionRunParams};
use cosmolkit_core::SanitizeError;
use cosmolkit_smiles::SmilesParseError;
use std::collections::{HashMap, HashSet, VecDeque};
use std::sync::{Arc, LazyLock};

#[derive(Debug, thiserror::Error)]
pub enum RecapDecompositionError {
    #[error("RECAP reaction failed: {0}")]
    Reaction(#[from] ReactionRunError),
    #[error("accepted RECAP product has no SMILES because sanitization failed: {0}")]
    UnkeyedProduct(#[source] SanitizeError),
    #[error("could not write RECAP product SMILES: {0}")]
    Smiles(#[from] SmilesParseError),
    #[error("RECAP product SMILES is not valid UTF-8: {0}")]
    Encoding(#[from] std::string::FromUtf8Error),
}

#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
pub enum RecapLookupError {
    #[error("unknown RECAP product: {product}")]
    UnknownProduct { product: String },
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct Neighbors {
    parents: HashSet<Arc<str>>,
    children: HashSet<Arc<str>>,
}

/// Direct parent–product relationships in the RECAP decomposition of one molecule.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RecapDecomposition {
    root: Arc<str>,
    nodes: HashMap<Arc<str>, Neighbors>,
}

impl RecapDecomposition {
    /// Canonical isomeric SMILES of the input molecule.
    pub fn root(&self) -> &str {
        &self.root
    }

    fn node(&self, product: &str) -> Result<&Neighbors, RecapLookupError> {
        self.nodes
            .get(product)
            .ok_or_else(|| RecapLookupError::UnknownProduct {
                product: product.into(),
            })
    }

    /// Unique immediate precursors of this canonical product, in unspecified order.
    pub fn parents(&self, product: &str) -> Result<Vec<String>, RecapLookupError> {
        Ok(self
            .node(product)?
            .parents
            .iter()
            .map(ToString::to_string)
            .collect())
    }

    /// Unique direct products of this canonical precursor, in unspecified order.
    pub fn children(&self, product: &str) -> Result<Vec<String>, RecapLookupError> {
        Ok(self
            .node(product)?
            .children
            .iter()
            .map(ToString::to_string)
            .collect())
    }

    /// Unique terminal descendants in unspecified order. A childless root yields no leaves.
    pub fn leaves(&self) -> Vec<String> {
        // RDKit source: rdkit/Chem/Recap.py, GetLeaves
        // RDKit✔️🔝:     for smi, child in self.children.items():
        //       if not len(child.children):
        //         res[smi] = child
        // All interned nodes are reachable. Scanning them once avoids revisiting shared paths.
        self.nodes
            .iter()
            .filter(|(key, node)| **key != self.root && node.children.is_empty())
            .map(|(key, _)| key.to_string())
            .collect()
    }

    fn connect(&mut self, parent: &Arc<str>, child: String) -> Option<Arc<str>> {
        // RDKit❗✔️:                 pNode.parents[nSmi] = weakref.proxy(node)
        //                 node.children[pSmi] = pNode
        let (child, new) = if let Some((key, _)) = self.nodes.get_key_value(child.as_str()) {
            (key.clone(), false)
        } else {
            let key = Arc::<str>::from(child);
            self.nodes.insert(key.clone(), Neighbors::default());
            (key, true)
        };
        self.nodes
            .get_mut(&child)
            .expect("interned product")
            .parents
            .insert(parent.clone());
        self.nodes
            .get_mut(parent)
            .expect("queued precursor")
            .children
            .insert(child.clone());
        new.then_some(child)
    }
}

static REACTIONS: LazyLock<Vec<Reaction>> = LazyLock::new(|| {
    // RDKit source: rdkit/Chem/Recap.py, reactionDefs (fece8caa).
    // RDKit✔️✔️: reactions = tuple([Reactions.ReactionFromSmarts(x) for x in reactionDefs])
    [
        "[#7;+0;D2,D3:1]!@C(!@=O)!@[#7;+0;D2,D3:2]>>*[#7:1].[#7:2]*",
        "[C;!$(C([#7])[#7]):1](=!@[O:2])!@[#7;+0;!D1:3]>>*[C:1]=[O:2].*[#7:3]",
        "[C:1](=!@[O:2])!@[O;+0:3]>>*[C:1]=[O:2].[O:3]*",
        "[N;!D1;+0;!$(N-C=[#7,#8,#15,#16])](-!@[*:1])-!@[*:2]>>*[*:1].[*:2]*",
        "[#7;R;D3;+0:1]-!@[*:2]>>*[#7:1].[*:2]*",
        "[#6:1]-!@[O;+0]-!@[#6:2]>>[#6:1]*.*[#6:2]",
        "[C:1]=!@[C:2]>>[C:1]*.*[C:2]",
        "[n;+0:1]-!@[C:2]>>[n:1]*.[C:2]*",
        "[O:3]=[C:4]-@[N;+0:1]-!@[C:2]>>[O:3]=[C:4]-[N:1]*.[C:2]*",
        "[c:1]-!@[c:2]>>[c:1]*.*[c:2]",
        "[n;+0:1]-!@[c:2]>>[n:1]*.*[c:2]",
        "[#7;+0;D2,D3:1]-!@[S:2](=[O:3])=[O:4]>>[#7:1]*.*[S:2](=[O:3])=[O:4]",
    ]
    .into_iter()
    .map(|smarts| {
        let mut reaction = crate::parse_smirks(smarts).expect("valid built-in RECAP template");
        crate::runner::initialize_for_run(&mut reaction).expect("valid built-in RECAP reaction");
        reaction
    })
    .collect()
});

/// Apply all default RECAP reactions to each newly discovered product.
pub fn decompose_recap(
    source: ReactionInput<'_>,
) -> Result<RecapDecomposition, RecapDecompositionError> {
    // RDKit source: rdkit/Chem/Recap.py, RecapDecompose
    // RDKit❗✔️:   mSmi = Chem.MolToSmiles(mol)
    // RDKit❗✔️:   activePool = {mSmi: res}
    //   allNodes[mSmi] = res
    //   while activePool:
    //     nSmi = next(iter(activePool))
    //     node = activePool.pop(nSmi)
    let root = Arc::<str>::from(canonical_key::<RecapDecompositionError>(source)?);
    let mut graph = RecapDecomposition {
        root: root.clone(),
        nodes: HashMap::from([(root.clone(), Neighbors::default())]),
    };
    let mut queue = VecDeque::from([(root, Reactant::Source(source))]);
    let params = ReactionRunParams::default();

    while let Some((parent, molecule)) = queue.pop_front() {
        for reaction in REACTIONS.iter() {
            // RDKit❗✔️:       ps = reaction.RunReactants((node.mol, ))
            'sequence: for mut products in
                crate::runner::run_reactants_source(reaction, &[molecule.as_input()], &params)?
            {
                // RDKit❗✔️:           tSeq = [(prod.GetNumAtoms(onlyExplicit=True), idx) for idx, prod in enumerate(prodSeq)]
                //           tSeq.sort()
                products.sort_by_key(|p| p.topology.atoms.len());
                let mut candidates = Vec::with_capacity(products.len());
                for raw in products {
                    // RDKit❗✔️:               Chem.SanitizeMol(prod)
                    //             except Exception:
                    //               continue
                    let product = match Product::sanitize(raw) {
                        Ok(value) => value,
                        Err(error) => {
                            candidates.push(Err(error));
                            continue;
                        }
                    };
                    let key = canonical_key::<RecapDecompositionError>(product.as_input())?;
                    // RDKit❗✔️:             elif pSmi.replace('*', '').replace('()', '') in ('', 'C', 'CC', 'CCC'):
                    //               seqOk = False
                    //               break
                    if matches!(
                        key.replace('*', "").replace("()", "").as_str(),
                        "" | "C" | "CC" | "CCC"
                    ) {
                        continue 'sequence;
                    }
                    candidates.push(Ok((key, product)));
                }

                // RDKit❗✔️:           if seqOk:
                //             for nats, prod in prodSeq:
                //               pSmi = prod.pSmi
                // RDKit❗✔️:                 activePool[pSmi] = pNode
                //                 allNodes[pSmi] = pNode
                for item in candidates {
                    let (key, product) = item.map_err(RecapDecompositionError::UnkeyedProduct)?;
                    if let Some(new) = graph.connect(&parent, key) {
                        queue.push_back((new, Reactant::Product(product)));
                    }
                }
            }
        }
    }
    Ok(graph)
}
