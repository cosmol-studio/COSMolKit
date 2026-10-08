//! Reaction operation projection; canonical runtime owns COW and atomic commits.
use crate::Molecule;
use cosmolkit as ck;
#[derive(Clone)]
pub struct ReactionApplyResult {
    molecule: Molecule,
    changed: bool,
}
impl ReactionApplyResult {
    pub fn molecule(&self) -> Molecule {
        self.molecule.clone()
    }
    pub fn changed(&self) -> bool {
        self.changed
    }
}
impl Molecule {
    pub fn reaction_products(
        &self,
        reaction: &mut ck::Reaction,
        reactant_template: usize,
    ) -> Result<Vec<Vec<Molecule>>, ck::OperationError> {
        self.inner
            .borrow()
            .reaction_products(reaction, reactant_template)
            .map(sets)
    }
    pub fn reaction_products_with_params(
        &self,
        reaction: &mut ck::Reaction,
        reactant_template: usize,
        params: &ck::ReactionSingleRunParams,
    ) -> Result<Vec<Vec<Molecule>>, ck::OperationError> {
        self.inner
            .borrow()
            .reaction_products_with_params(reaction, reactant_template, params)
            .map(sets)
    }
    pub fn reaction_products_from_inputs(
        &self,
        reaction: &mut ck::Reaction,
        reactants: &[&Molecule],
        params: &ck::ReactionRunParams,
    ) -> Result<Vec<Vec<Molecule>>, ck::OperationError> {
        let borrowed: Vec<_> = reactants.iter().map(|m| m.inner.borrow()).collect();
        let inputs: Vec<_> = borrowed.iter().map(|m| &**m).collect();
        self.inner
            .borrow()
            .reaction_products_from_inputs(reaction, &inputs, params)
            .map(sets)
    }
    pub fn apply_reaction(
        &self,
        reaction: &mut ck::Reaction,
    ) -> Result<ReactionApplyResult, ck::OperationError> {
        self.inner
            .borrow()
            .apply_reaction(reaction)
            .map(|result| ReactionApplyResult {
                molecule: Molecule {
                    inner: result.molecule().clone().into(),
                },
                changed: result.changed(),
            })
    }
    pub fn apply_reaction_with_params(
        &self,
        reaction: &mut ck::Reaction,
        params: &ck::ReactionApplyParams,
    ) -> Result<ReactionApplyResult, ck::OperationError> {
        self.inner
            .borrow()
            .apply_reaction_with_params(reaction, params)
            .map(|result| ReactionApplyResult {
                molecule: Molecule {
                    inner: result.molecule().clone().into(),
                },
                changed: result.changed(),
            })
    }
    pub fn apply_reaction_(&self, reaction: &mut ck::Reaction) -> Result<bool, ck::OperationError> {
        self.inner.borrow_mut().apply_reaction_(reaction)
    }
    pub fn apply_reaction_with_params_(
        &self,
        reaction: &mut ck::Reaction,
        params: &ck::ReactionApplyParams,
    ) -> Result<bool, ck::OperationError> {
        self.inner
            .borrow_mut()
            .apply_reaction_with_params_(reaction, params)
    }
}
fn sets(sets: Vec<Vec<ck::Molecule>>) -> Vec<Vec<Molecule>> {
    sets.into_iter()
        .map(|set| {
            set.into_iter()
                .map(|inner| Molecule {
                    inner: inner.into(),
                })
                .collect()
        })
        .collect()
}
