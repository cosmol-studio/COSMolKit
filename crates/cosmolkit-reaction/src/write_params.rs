use crate::ReactionCoordinateSelection;
use cosmolkit_smiles::CxSmilesFields;
/// Source-consumed SMARTS/CX options; canonical sorts strings within each role.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReactionWriteParams {
    pub canonical: bool,
    pub isomeric_smiles: bool,
    pub rooted_at_atom: Option<usize>,
    pub include_dative_bonds: bool,
    pub include_cx: bool,
    pub cx_fields: CxSmilesFields,
    /// Per template in original reactant, agent, product insertion order.
    pub coordinate_selections: Vec<ReactionCoordinateSelection>,
}
impl Default for ReactionWriteParams {
    fn default() -> Self {
        // RDKit❗✔️: inline std::string ChemicalReactionToRxnSmarts(const ChemicalReaction &rxn) {
        // RDKit❗✔️:   SmilesWriteParams params;
        // RDKit❗✔️:   params.canonical = false;
        // RDKit❗✔️:   return ChemicalReactionToRxnSmarts(rxn, params);
        // RDKit❗✔️: }
        Self {
            canonical: false,
            isomeric_smiles: true,
            rooted_at_atom: None,
            include_dative_bonds: true,
            include_cx: false,
            cx_fields: CxSmilesFields::ALL,
            coordinate_selections: Vec::new(),
        }
    }
}
