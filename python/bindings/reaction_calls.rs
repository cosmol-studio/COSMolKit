//! Per-invocation callback transport; reaction chemistry stays in Rust.
use crate::native::Molecule;
use crate::reaction::Reaction;
use cosmolkit as ck;
use pyo3::prelude::*;

struct MatchingScope<'a> {
    reaction: &'a mut ck::Reaction,
    original: ck::SubstructMatchParams,
}
impl<'a> MatchingScope<'a> {
    fn new(reaction: &'a mut ck::Reaction, params: &ck::SubstructMatchParams) -> Self {
        let original = reaction.match_params().clone();
        reaction.set_match_params(params);
        Self { reaction, original }
    }
}
impl Drop for MatchingScope<'_> {
    fn drop(&mut self) {
        self.reaction.set_match_params(&self.original);
    }
}

impl Reaction {
    pub(crate) fn needs_callback_target(&self, py: Python<'_>) -> bool {
        self.copy_match_owner()
            .is_some_and(|owner| owner.borrow(py).has_final_match())
    }
    pub(crate) fn invoke<T: Send>(
        &mut self,
        py: Python<'_>,
        targets: &[&ck::Molecule],
        run: impl FnOnce(&mut ck::Reaction) -> Result<T, ck::OperationError> + Send,
    ) -> PyResult<T> {
        let owner = self.copy_match_owner();
        let result = if let Some(owner) = owner {
            owner
                .borrow(py)
                .with_callback_values(py, targets, |params| {
                    let scope = MatchingScope::new(self.native_mut(), params);
                    py.detach(|| run(scope.reaction))
                })?
        } else {
            let native = self.native_mut();
            py.detach(|| run(native))
        };
        result.map_err(crate::operation_errors::OperationError::from_rust)
    }
}

pub(crate) fn products(values: Vec<Vec<ck::Molecule>>) -> Vec<Vec<Molecule>> {
    values
        .into_iter()
        .map(|set| set.into_iter().map(|inner| Molecule { inner }).collect())
        .collect()
}

pub(crate) fn run(
    reaction: &mut Reaction,
    py: Python<'_>,
    targets: &[Py<Molecule>],
    params: &ck::ReactionRunParams,
) -> PyResult<Vec<Vec<Molecule>>> {
    let owners = targets
        .iter()
        .map(|target| target.borrow(py))
        .collect::<Vec<_>>();
    let inputs: Vec<_> = owners.iter().map(|owner| &owner.inner).collect();
    reaction
        .invoke(py, &inputs, |native| native.run(&inputs, params))
        .map(products)
}

pub(crate) fn keywords(
    max_products: Option<u32>,
    coordinate_selections: Option<Vec<crate::construction_values::ReactionCoordinateSelection>>,
    copy_atom_properties: Option<bool>,
) -> ck::ReactionRunParams {
    let defaults = ck::ReactionRunParams::default();
    ck::ReactionRunParams {
        max_products: max_products.unwrap_or(defaults.max_products),
        coordinate_selections: coordinate_selections
            .map(|values| values.into_iter().map(|value| value.inner).collect())
            .unwrap_or(defaults.coordinate_selections),
        copy_atom_properties: copy_atom_properties.unwrap_or(defaults.copy_atom_properties),
    }
}

pub(crate) fn run_from_inputs(
    molecule: &Molecule,
    reaction: &mut Reaction,
    py: Python<'_>,
    targets: &[Py<Molecule>],
    params: &ck::ReactionRunParams,
) -> PyResult<Vec<Vec<Molecule>>> {
    let owners = targets
        .iter()
        .map(|value| value.try_borrow(py))
        .collect::<Result<Vec<_>, _>>()?;
    let inputs = owners.iter().map(|value| &value.inner).collect::<Vec<_>>();
    reaction
        .invoke(py, &inputs, |reaction| {
            molecule
                .inner
                .reaction_products_from_inputs(reaction, &inputs, params)
        })
        .map(products)
}

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::{
        ffi::c_str,
        types::{PyDict, PyModule},
    };

    #[test]
    fn reaction_execution_keeps_cow_targets_scoped_callbacks_results_and_sources() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module=PyModule::new(py,"reaction_invocation_probe")?;
            module.add_class::<Molecule>()?;
            module.add_class::<Reaction>()?;
            module.add_class::<crate::search_callbacks::SubstructMatchParams>()?;
            module.add_class::<crate::optional_configuration::ReactionRunParams>()?;
            module.add_class::<crate::search_reaction_params::ReactionSingleRunParams>()?;
            module.add_class::<crate::params::ReactionApplyParams>()?;
            module.add_class::<crate::construction_values::ReactionCoordinateSelection>()?;
            module.add("OperationError", <crate::operation_errors::OperationError as crate::error_framework::ExceptionPayload>::exception_type(py)?)?;
            let source=ck::Molecule::from_smiles("C").unwrap();
            // Builder publication intentionally starts with fresh derived
            // state. Prepare valence explicitly for the callback's AddHs;
            // neither copy nor the binding invents that chemical preparation.
            let original_a=source.to_builder().with_name("a".into()).build().unwrap().with_assigned_valence().unwrap();
            let original_b=source.to_builder().with_name("b".into()).build().unwrap().with_assigned_valence().unwrap();
            let a=Py::new(py,Molecule{inner:original_a.clone()})?;
            let b=Py::new(py,Molecule{inner:original_b.clone()})?;
            assert!(std::ptr::eq(a.borrow(py).inner.topology(), b.borrow(py).inner.topology()));
            let scope=PyDict::new(py);
            scope.set_item("ck",module)?;
            scope.set_item("a",a.clone_ref(py))?;
            scope.set_item("b",b.clone_ref(py))?;
            py.run(c_str!(r#"
seen = []
atoms = []
def final(target, indices):
    assert target is not a and target is not b
    assert target.num_atoms == 1
    assert tuple(indices) == (0,)
    seen.append(target)
    return True
def atom(query, target):
    atoms.append(target.atomic_number)
    return True
matching = ck.SubstructMatchParams(final_match=final, atom_match=atom)
rxn = ck.Reaction.from_smirks("[C:1].[C:2]>>[N:1].[O:2]")
rxn.set_match_params(matching)
before = rxn.match_params
before.max_matches = 3
assert rxn.match_params.max_matches != 3
params = ck.ReactionRunParams()
for operation in (
    lambda: rxn.run([a, b]),
    lambda: rxn.run([a, b], params),
    lambda: rxn.run([a, b], max_products=1000, copy_atom_properties=False),
    lambda: a.reaction_products_from_inputs(rxn, [a, b]),
    lambda: a.reaction_products_from_inputs(rxn, [a, b], params),
    lambda: a.reaction_products_from_inputs(rxn, [a, b], max_products=1000),
):
    seen.clear()
    products = operation()
    assert len(products) == 1 and len(products[0]) == 2
    assert [target.properties.name for target in seen] == ["a", "b"]
    assert atoms
    assert rxn.match_params.final_match is final
    assert rxn.match_params.atom_match is atom
    assert rxn.is_initialized
for operation in (
    lambda: rxn.run([a, b], params, max_products=2),
    lambda: a.reaction_products_from_inputs(rxn, [a, b], params, max_products=2),
):
    try: operation()
    except TypeError: pass
    else: raise AssertionError("mixed configuration accepted")
marker = LookupError("exact reaction callback error")
def fail(target, indices): raise marker
rxn.set_match_params(ck.SubstructMatchParams(final_match=fail))
try: rxn.run([a,b])
except LookupError as error: assert error is marker
else: raise AssertionError("callback error lost")
rxn.set_match_params(matching)
assert len(rxn.run([a,b])) == 1
try: rxn.run([a])
except ck.OperationError as error:
    assert error.cause is error.__cause__
    assert error.cause.expected == 2 and error.cause.actual == 1
else: raise AssertionError("native arity failure lost")

single = ck.Reaction.from_smirks("[C:1]>>[N:1]")
single.set_match_params(matching)
single_params = ck.ReactionSingleRunParams()
for operation in (
    lambda: a.reaction_products(single, 0),
    lambda: a.reaction_products(single, 0, single_params),
    lambda: a.reaction_products(single, 0, coordinate_selection=ck.ReactionCoordinateSelection.auto()),
    lambda: a.reaction_products_with_params(single, 0, single_params),
):
    seen.clear()
    assert len(operation()) == 1
    assert seen[0] is not a and seen[0].to_binary() == a.to_binary()
apply_params = ck.ReactionApplyParams()
for operation in (
    lambda: a.apply_reaction(single),
    lambda: a.apply_reaction(single, apply_params),
    lambda: a.apply_reaction(single, remove_unmatched_atoms=True),
    lambda: a.apply_reaction_with_params(single, apply_params),
):
    seen.clear()
    result = operation()
    assert result.changed is True
    assert result.molecule.atoms[0].atomic_number == 7
    assert seen[0] is not a and seen[0].to_binary() == a.to_binary()
    assert a.atoms[0].atomic_number == 6
view = single.reactant_template(0)
saved = single.match_params
assert saved.final_match is final

# The native in-place call is not replaced with value execution/writeback.
# Callback values may mutate and escape, while the receiver stays exclusive.
retained = []
def local_edit(target, indices):
    assert target is not active and target.num_atoms == 1
    assert tuple(indices) == (0,)
    try: active.num_atoms
    except RuntimeError: pass
    else: raise AssertionError("original receiver was reentered during mutation")
    target.add_hydrogens_()
    assert target.num_atoms == 5
    retained.append(target)
    return True
single.set_match_params(ck.SubstructMatchParams(final_match=local_edit))
for operation in (
    lambda: active.apply_reaction_(single),
    lambda: active.apply_reaction_(single, apply_params),
    lambda: active.apply_reaction_(single, remove_unmatched_atoms=True),
    lambda: active.apply_reaction_with_params_(single, apply_params),
):
    active = a.copy()
    assert operation() is True
    assert active.num_atoms == 1 and active.atoms[0].atomic_number == 7
    assert a.atoms[0].atomic_number == 6
    assert retained[-1].atoms[0].atomic_number == 6
    assert retained[-1].num_atoms == 5
active = a.copy()
try: active.apply_reaction_(single, apply_params, remove_unmatched_atoms=False)
except TypeError: pass
else: raise AssertionError("mixed in-place configuration accepted")
marker = LookupError("original in-place callback failure")
single.set_match_params(ck.SubstructMatchParams(final_match=fail))
try: active.apply_reaction_(single)
except LookupError as error: assert error is marker
else: raise AssertionError("original in-place exception lost")
single.set_match_params(ck.SubstructMatchParams(final_match=lambda target, indices: False))
before = active.to_binary()
assert active.apply_reaction_(single) is False
assert active.to_binary() == before
single.set_match_params(ck.SubstructMatchParams())
assert active.apply_reaction_(single) is True
assert retained[0].num_atoms == 5
# Restore the ordinary matching owner for the remaining lifetime checks.
single.set_match_params(matching)
"#),Some(&scope),None)?;

            let view=scope.get_item("view")?.unwrap().extract::<Py<crate::query_compiler::QueryGraph>>()?;
            let snapshot=view.borrow(py).inner().clone();
            scope.get_item("single")?.unwrap().call_method1("run",(vec![a.clone_ref(py)],))?;
            assert_eq!(view.borrow(py).inner(),&snapshot);
            py.run(c_str!("del single\nassert saved.final_match is final\nassert view.num_atoms == 1"),Some(&scope),None)?;

            let projected=scope.get_item("rxn")?.unwrap().extract::<Py<Reaction>>()?;
            let pointer=projected.borrow(py).inner() as *const ck::Reaction;
            let result=projected.bind(py).call_method1("run",(vec![a.clone_ref(py),b.clone_ref(py)],))?.extract::<Vec<Vec<Py<Molecule>>>>()?;
            assert_eq!(projected.borrow(py).inner() as *const ck::Reaction,pointer,"unique reactions must not be cloned per invocation");
            let expected=ck::Reaction::from_smirks("[C:1].[C:2]>>[N:1].[O:2]").unwrap().run(&[&original_a,&original_b],&ck::ReactionRunParams::default()).unwrap();
            assert_eq!(result.len(),expected.len());
            for (actual,expected) in result.iter().flatten().zip(expected.iter().flatten()) {
                assert_eq!(actual.borrow(py).inner.to_binary().unwrap(),expected.to_binary().unwrap());
            }
            assert_eq!(a.borrow(py).inner,original_a);
            assert_eq!(b.borrow(py).inner,original_b);
            Ok(())
        }).unwrap();
    }
}
