use crate::ReactionCoordinateSelection;
use cosmolkit_smiles::CxSmilesFields;
/// Source-consumed SMARTS/CX options; canonical sorts strings within each role.
#[derive(Debug, Clone, PartialEq, Eq, cosmolkit_macros::Configuration)]
pub struct ReactionWriteParams {
    /// Canonicalize component output and sort strings within each reaction role.
    pub canonical: bool,
    /// Include isotope and stereochemical information.
    pub isomeric_smiles: bool,
    /// Optional starting atom index for template output.
    pub rooted_at_atom: Option<usize>,
    /// Write dative bonds using their explicit notation.
    pub include_dative_bonds: bool,
    /// Include the selected CXSMILES extension fields.
    pub include_cx: bool,
    #[configuration(none_default)]
    /// Bit mask of CXSMILES extensions to write.
    pub cx_fields: CxSmilesFields,
    /// Per template in original reactant, agent, product insertion order.
    #[configuration(none_default)]
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
