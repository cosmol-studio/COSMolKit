//! Registered molecule fragment operations retain native ordering and ownership.
use crate::{
    Molecule, fields::Field, fragment_errors::BondCutError, input,
    operation_errors::OperationError, storage::NativeField,
};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::prelude::*;
crate::fields::value_array!(ck::Molecule, Molecules, "ReadonlyArray<Molecule>");

#[wasm_bindgen]
pub struct BondCut(ck::BondCut);
crate::fields::value_field!(ck::BondCut, BondCut, BondCutInput, "BondCut");
crate::fields::value_array!(ck::BondCut, BondCuts, "ReadonlyArray<BondCut>");

#[cfg(feature = "cap-fragmentation")]
mod decomposition {
    use super::{BondCut, Molecule};
    use crate::{
        errors::{enum_error, exception},
        input,
        match_errors::SubstructMatchError,
        storage::ReadonlyValue,
    };
    use wasm_bindgen::JsCast;

    #[wasm_bindgen]
    extern "C" {
        #[wasm_bindgen(
            typescript_type = "ReactionRunError | SanitizeError | SmilesParseError | FromUtf8Error"
        )]
        #[derive(Clone)]
        pub type DecompositionCause;
    }

    fn brics_cause(
        owner: &ReadonlyValue<ck::BricsDecompositionError>,
    ) -> Option<DecompositionCause> {
        use ck::BricsDecompositionError as E;
        Some(match owner.as_ref() {
            E::Reaction(_) => crate::errors::child!(
                owner,
                E::Reaction,
                crate::reaction_execution_errors::ReactionRunError
            ),
            E::UnkeyedProduct(_) => crate::errors::child!(
                owner,
                E::UnkeyedProduct,
                crate::sanitize_errors::SanitizeError
            ),
            E::Smiles(_) => crate::errors::child!(
                owner,
                E::Smiles,
                crate::smiles_parse_errors::SmilesParseError
            ),
            E::Encoding(_) => {
                crate::errors::child!(owner, E::Encoding, crate::external_errors::FromUtf8Error)
            }
        })
    }
    enum_error!(BricsDecompositionError, BricsDecompositionErrorKind, ck::BricsDecompositionError, {
        Reaction => ck::BricsDecompositionError::Reaction(_),
        UnkeyedProduct => ck::BricsDecompositionError::UnkeyedProduct(_),
        Smiles => ck::BricsDecompositionError::Smiles(_),
        Encoding => ck::BricsDecompositionError::Encoding(_),
    }, {}, source = DecompositionCause, cause = brics_cause);

    fn recap_cause(
        owner: &ReadonlyValue<ck::RecapDecompositionError>,
    ) -> Option<DecompositionCause> {
        use ck::RecapDecompositionError as E;
        Some(match owner.as_ref() {
            E::Reaction(_) => crate::errors::child!(
                owner,
                E::Reaction,
                crate::reaction_execution_errors::ReactionRunError
            ),
            E::UnkeyedProduct(_) => crate::errors::child!(
                owner,
                E::UnkeyedProduct,
                crate::sanitize_errors::SanitizeError
            ),
            E::Smiles(_) => crate::errors::child!(
                owner,
                E::Smiles,
                crate::smiles_parse_errors::SmilesParseError
            ),
            E::Encoding(_) => {
                crate::errors::child!(owner, E::Encoding, crate::external_errors::FromUtf8Error)
            }
        })
    }
    enum_error!(RecapDecompositionError, RecapDecompositionErrorKind, ck::RecapDecompositionError, {
        Reaction => ck::RecapDecompositionError::Reaction(_),
        UnkeyedProduct => ck::RecapDecompositionError::UnkeyedProduct(_),
        Smiles => ck::RecapDecompositionError::Smiles(_),
        Encoding => ck::RecapDecompositionError::Encoding(_),
    }, {}, source = DecompositionCause, cause = recap_cause);

    enum_error!(RecapLookupError,RecapLookupErrorKind,ck::RecapLookupError,{
        UnknownProduct => ck::RecapLookupError::UnknownProduct{..},
    },{
        #[wasm_bindgen(getter)] pub fn product(&self)->String{
            let ck::RecapLookupError::UnknownProduct{product}=self.inner.as_ref();
            product.clone()
        }
    });

    #[wasm_bindgen]
    pub struct RecapDecomposition(ck::RecapDecomposition);

    #[wasm_bindgen]
    impl RecapDecomposition {
        #[wasm_bindgen(getter)]
        pub fn root(&self) -> String {
            self.0.root().into()
        }

        pub fn parents(&self, product: Ts<input::TextInput>) -> Result<Vec<String>, JsValue> {
            let product = input::read::<String>(product)?;
            self.0
                .parents(&product)
                .map_err(RecapLookupError::from_rust)
        }

        pub fn children(&self, product: Ts<input::TextInput>) -> Result<Vec<String>, JsValue> {
            let product = input::read::<String>(product)?;
            self.0
                .children(&product)
                .map_err(RecapLookupError::from_rust)
        }

        pub fn leaves(&self) -> Vec<String> {
            self.0.leaves()
        }
    }
    use cosmolkit as ck;
    use tsify::{Ts, Tsify};
    use wasm_bindgen::prelude::*;

    #[derive(serde::Serialize, Tsify)]
    #[serde(transparent)]
    pub struct BricsEnds([(usize, u8); 2]);

    macro_rules! project_brics_match {
        ($name:ident) => {
            #[wasm_bindgen]
            pub struct $name(ck::$name);

            #[wasm_bindgen]
            impl $name {
                #[wasm_bindgen(getter)]
                pub fn bond(&self) -> usize {
                    self.0.bond().index()
                }

                /// Two [atom index, environment label] pairs in rule-match order.
                #[wasm_bindgen(getter)]
                pub fn ends(&self) -> Result<tsify::Ts<BricsEnds>, JsValue> {
                    BricsEnds(self.0.ends().map(|(atom, label)| (atom.index(), label)))
                        .into_ts()
                        .map_err(JsError::from)
                        .map_err(Into::into)
                }

                #[wasm_bindgen(js_name=toBondCut)]
                pub fn to_bond_cut(&self) -> BondCut {
                    BondCut(self.0.to_bond_cut())
                }
            }
        };
    }
    project_brics_match!(BricsBondMatch);
    project_brics_match!(RBricsBondMatch);

    #[wasm_bindgen]
    impl Molecule {
        /// Return the RECAP product hierarchy of this molecule.
        #[wasm_bindgen(js_name=decomposeRecap)]
        pub fn decompose_recap(&self) -> Result<RecapDecomposition, JsValue> {
            self.0
                .decompose_recap()
                .map(RecapDecomposition)
                .map_err(RecapDecompositionError::from_rust)
        }

        /// Return unique terminal canonical SMILES in unspecified order.
        #[wasm_bindgen(js_name=decomposeBrics)]
        pub fn decompose_brics(&self) -> Result<Vec<String>, JsValue> {
            self.0
                .decompose_brics()
                .map_err(BricsDecompositionError::from_rust)
        }

        /// Find standard BRICS sites on this topology; rediscover them after editing.
        #[wasm_bindgen(js_name=findBricsBonds)]
        pub fn find_brics_bonds(&self) -> Result<Vec<BricsBondMatch>, JsValue> {
            self.0
                .find_brics_bonds()
                .map(|rows| rows.into_iter().map(BricsBondMatch).collect())
                .map_err(SubstructMatchError::from_rust)
        }
        /// Return unique terminal regular r-BRICS SMILES in unspecified order.
        #[wasm_bindgen(js_name=decomposeRBrics)]
        pub fn decompose_r_brics(&self) -> Result<Vec<String>, JsValue> {
            self.0
                .decompose_r_brics()
                .map_err(BricsDecompositionError::from_rust)
        }

        /// Find regular r-BRICS sites on this topology; rediscover them after editing.
        #[wasm_bindgen(js_name=findRBricsBonds)]
        pub fn find_r_brics_bonds(&self) -> Result<Vec<RBricsBondMatch>, JsValue> {
            self.0
                .find_r_brics_bonds()
                .map(|rows| rows.into_iter().map(RBricsBondMatch).collect())
                .map_err(SubstructMatchError::from_rust)
        }
    }
}

#[wasm_bindgen]
impl BondCut {
    /// Supply both attachment-side isotope labels, or omit both to label each dummy
    /// with the index of the atom it replaces. Effective labels must fit in 16 bits.
    #[wasm_bindgen(js_name=atBond)]
    pub fn at_bond(
        bond: Ts<input::IndexInput>,
        begin_label: Option<Ts<input::IndexInput>>,
        end_label: Option<Ts<input::IndexInput>>,
    ) -> Result<Self, JsValue> {
        ck::BondCut::at_bond(
            ck::BondId::new(input::read::<usize>(bond)?),
            begin_label.map(input::read::<usize>).transpose()?,
            end_label.map(input::read::<usize>).transpose()?,
        )
        .map(Self)
        .map_err(BondCutError::from_rust)
    }
}

#[wasm_bindgen]
impl Molecule {
    /// Return the cut molecule; use fragments() to extract connected components.
    #[wasm_bindgen(js_name=withBondsCut)]
    pub fn with_bonds_cut(&self, cuts: BondCuts) -> Result<Self, JsValue> {
        self.0
            .with_bonds_cut(&Vec::<ck::BondCut>::input(cuts)?)
            .map(Self)
            .map_err(OperationError::from_rust)
    }

    /// Apply the same selected cuts to this molecule in place.
    #[wasm_bindgen(js_name=cutBonds_)]
    pub fn cut_bonds_(&mut self, cuts: BondCuts) -> Result<(), JsValue> {
        self.0
            .cut_bonds_(&Vec::<ck::BondCut>::input(cuts)?)
            .map_err(OperationError::from_rust)
    }

    #[wasm_bindgen(js_name=fragments)]
    pub fn fragments(&self) -> Result<Molecules, JsValue> {
        Vec::<ck::Molecule>::get(NativeField::owned(
            self.0.fragments().map_err(OperationError::from_rust)?,
        ))
    }
    #[wasm_bindgen(js_name=largestFragment)]
    pub fn largest_fragment(&self) -> Result<Self, JsValue> {
        self.0
            .largest_fragment()
            .map(Self)
            .map_err(OperationError::from_rust)
    }
}
