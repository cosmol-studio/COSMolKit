//! Typed projection configuration for native chemistry-state operations.
use cosmolkit as ck;
use pyo3::prelude::*;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::Molecule;
    use crate::operation_errors::OperationError;
    use pyo3::ffi::c_str;
    use pyo3::types::{PyDict, PyList, PyModule};
    #[test]
    fn chemistry_state_families_match_native_value_inplace_options_payloads_and_order() {
        Python::initialize();
        Python::attach(|py|->PyResult<()>{
 let native=ck::Molecule::from_smiles("Cc1ccccc1.C[C@H](F)Cl").unwrap();
 let mut embed=ck::EmbedParams::etkdg_v3();embed.random_seed=17;
 let native=native.with_3d_conformer_with_params(&embed).unwrap();
 let before=native.to_binary().unwrap();
 let cases=PyList::empty(py);
{
 let config=ck::SanitizeParams {operations:ck::SanitizeOperations::NONE};
 let kw=PyDict::new(py);kw.set_item("operations",Py::new(py,crate::sanitize_values::SanitizeOperations{inner:ck::SanitizeOperations::NONE})?)?;
 let expected=native.sanitize().unwrap().to_binary().unwrap();
 let configured=native.sanitize_with_params(&config).unwrap().to_binary().unwrap();
 cases.append(("sanitize","sanitize_","sanitize_with_params","sanitize_with_params_",Py::new(py,crate::sanitize_values::SanitizeParams{inner:config})?,kw,expected,configured))?;
}
{
 let config=ck::KekulizeParams::default();
 let kw=PyDict::new(py);kw.set_item("canonical",config.canonical)?;
 let expected=native.with_kekulized_bonds().unwrap().to_binary().unwrap();
 let configured=native.with_kekulized_bonds_with_params(&config).unwrap().to_binary().unwrap();
 cases.append(("with_kekulized_bonds","kekulize_bonds_","with_kekulized_bonds_with_params","kekulize_bonds_with_params_",Py::new(py,crate::params::KekulizeParams{inner:config})?,kw,expected,configured))?;
}
{
 let config=ck::AromaticityParams::default();
 let kw=PyDict::new(py);kw.set_item("model","rdkit")?;
 let expected=native.with_assigned_aromaticity().unwrap().to_binary().unwrap();
 let configured=native.with_assigned_aromaticity_with_params(&config).unwrap().to_binary().unwrap();
 cases.append(("with_assigned_aromaticity","assign_aromaticity_","with_assigned_aromaticity_with_params","assign_aromaticity_with_params_",Py::new(py,crate::selectors::AromaticityParams{inner:config})?,kw,expected,configured))?;
}
{
 let config=ck::ValenceParams::default();
 let kw=PyDict::new(py);kw.set_item("strict",config.strict)?;
 let expected=native.with_assigned_valence().unwrap().to_binary().unwrap();
 let configured=native.with_assigned_valence_with_params(&config).unwrap().to_binary().unwrap();
 cases.append(("with_assigned_valence","assign_valence_","with_assigned_valence_with_params","assign_valence_with_params_",Py::new(py,crate::selectors::ValenceParams{inner:config})?,kw,expected,configured))?;
}
{
 let config=ck::RingSearchParams::default();
 let kw=PyDict::new(py);kw.set_item("include_dative_bonds",config.include_dative_bonds)?;
 let expected=native.with_assigned_ring_families().unwrap().to_binary().unwrap();
 let configured=native.with_assigned_ring_families_with_params(&config).unwrap().to_binary().unwrap();
 cases.append(("with_assigned_ring_families","assign_ring_families_","with_assigned_ring_families_with_params","assign_ring_families_with_params_",Py::new(py,crate::params::RingSearchParams{inner:config})?,kw,expected,configured))?;
}
{
 let config=ck::StructureTagParams::default();
 let kw=PyDict::new(py);kw.set_item("replace_existing_tags",config.replace_existing_tags)?;
 let expected=native.with_chiral_tags_from_structure().unwrap().to_binary().unwrap();
 let configured=native.with_chiral_tags_from_structure_with_params(&config).unwrap().to_binary().unwrap();
 cases.append(("with_chiral_tags_from_structure","assign_chiral_tags_from_structure_","with_chiral_tags_from_structure_with_params","assign_chiral_tags_from_structure_with_params_",Py::new(py,crate::params::StructureTagParams{inner:config})?,kw,expected,configured))?;
}
{
 let config=ck::CipLabelOptions::default();
 let kw=PyDict::new(py);kw.set_item("max_recursive_iterations",config.max_recursive_iterations())?;
 let expected=native.with_cip_labels().unwrap().to_binary().unwrap();
 let configured=native.with_cip_labels_with_options(&config).unwrap().to_binary().unwrap();
 cases.append(("with_cip_labels","assign_cip_labels_","with_cip_labels_with_options","assign_cip_labels_with_options_",Py::new(py,crate::numerical_params::CipLabelOptions::from_rust(py, config)?)?,kw,expected,configured))?;
}
{
 let config=ck::AddHsParams {only_on_atoms:Some(vec![]),..Default::default()};
 let kw=PyDict::new(py);kw.set_item("only_on_atoms",Vec::<usize>::new())?;
 let expected=native.with_hydrogens().unwrap().to_binary().unwrap();
 let configured=native.with_hydrogens_with_params(&config).unwrap().to_binary().unwrap();
 cases.append(("with_hydrogens","add_hydrogens_","with_hydrogens_with_params","add_hydrogens_with_params_",Py::new(py,crate::chemical_params::AddHsParams::from_rust(py, config)?)?,kw,expected,configured))?;
}
{
 let config=ck::RemoveHsParams::default();
 let kw=PyDict::new(py);kw.set_item("sanitize",config.sanitize)?;
 let expected=native.without_hydrogens().unwrap().to_binary().unwrap();
 let configured=native.without_hydrogens_with_params(&config).unwrap().to_binary().unwrap();
 cases.append(("without_hydrogens","remove_hydrogens_","without_hydrogens_with_params","remove_hydrogens_with_params_",Py::new(py,crate::params::RemoveHsParams{inner:config})?,kw,expected,configured))?;
}

 let simple=PyList::empty(py);
simple.append(("with_assigned_radicals","assign_radicals_",native.with_assigned_radicals().unwrap().to_binary().unwrap()))?;
simple.append(("with_assigned_rings","assign_rings_",native.with_assigned_rings().unwrap().to_binary().unwrap()))?;

 let scaffolds=PyList::empty(py);
scaffolds.append(("murcko_scaffold",native.murcko_scaffold().unwrap().to_binary().unwrap()))?;
scaffolds.append(("net_scaffold",native.net_scaffold().unwrap().to_binary().unwrap()))?;
scaffolds.append(("murcko_decompose",native.murcko_decompose().unwrap().to_binary().unwrap()))?;
scaffolds.append(("largest_fragment",native.largest_fragment().unwrap().to_binary().unwrap()))?;

 let expected_fragments=native.fragments().unwrap().iter().map(|value|value.to_binary().unwrap()).collect::<Vec<_>>();
 let expected_potential=native.potential_stereo().unwrap();
 let expected_score=native.tautomer_score().unwrap();
 let empty_score=native.tautomer_score_with_params(&ck::TautomerScoreParams{terms:Some(vec![])}).unwrap();
 let sources=PyList::empty(py);
 for _ in 0..100{sources.append(Py::new(py,Molecule{inner:native.clone()})?)?;}
 let module=PyModule::new(py,"chemistry_state_test")?;
 module.add_class::<OperationError>()?;
 module.add_class::<crate::params::KekulizeParams>()?;
 module.add_class::<crate::params::PotentialStereoParams>()?;
 module.add_class::<crate::tautomer_values::TautomerScoreParams>()?;
 let locals=PyDict::new(py);
 locals.set_item("api",module)?;
 locals.set_item("sources",sources)?;
 locals.set_item("before",before)?;
 locals.set_item("cases",cases)?;
 locals.set_item("simple",simple)?;
 locals.set_item("scaffolds",scaffolds)?;
 locals.set_item("expected_fragments",expected_fragments)?;
 locals.set_item("expected_score",(expected_score.ring(),expected_score.substructure(),expected_score.hetero_hydrogen(),expected_score.total()))?;
 locals.set_item("empty_score",(empty_score.ring(),empty_score.substructure(),empty_score.hetero_hydrogen(),empty_score.total()))?;
 py.run(c_str!(r#"
def fresh(): return sources.pop()
for value, inplace, value_alias, inplace_alias, params, kw, expected, configured in cases:
    for name, mutates in [(value, False), (inplace, True)]:
        for form in ['default', 'params', 'keywords', 'alias']:
            source = fresh()
            if form == 'default': output = getattr(source, name)()
            elif form == 'params': output = getattr(source, name)(params)
            elif form == 'keywords': output = getattr(source, name)(**kw)
            else: output = getattr(source, inplace_alias if mutates else value_alias)(params)
            want = expected if form == 'default' else configured
            if mutates:
                assert output is None
                assert bytes(source.to_binary()) == bytes(want), (name, form)
            else:
                assert bytes(output.to_binary()) == bytes(want), (name, form)
                assert bytes(source.to_binary()) == bytes(before), name
        source = fresh()
        try: getattr(source, name)(params, **kw)
        except TypeError: pass
        else: raise AssertionError('mixed state configuration accepted')
        assert bytes(source.to_binary()) == bytes(before)
for value, inplace, expected in simple:
    source = fresh()
    output = getattr(source, value)()
    assert bytes(output.to_binary()) == bytes(expected)
    assert bytes(source.to_binary()) == bytes(before)
    assert getattr(source, inplace)() is None
    assert bytes(source.to_binary()) == bytes(expected)
for name, expected in scaffolds:
    source = fresh()
    output = getattr(source, name)()
    assert bytes(output.to_binary()) == bytes(expected)
    assert bytes(source.to_binary()) == bytes(before)
source = fresh()
assert [bytes(value.to_binary()) for value in source.fragments()] == [bytes(value) for value in expected_fragments]
assert bytes(source.to_binary()) == bytes(before)
potential = source.potential_stereo()
assert bytes(source.to_binary()) == bytes(before)
assert source.potential_stereo(api.PotentialStereoParams()).atom_ranks == potential.atom_ranks
assert source.potential_stereo(clean=True).atom_ranks == potential.atom_ranks
assert source.potential_stereo_with_params(api.PotentialStereoParams()).atom_ranks == potential.atom_ranks
def parts(score): return (score.ring(), score.substructure(), score.hetero_hydrogen(), score.total())
assert parts(source.tautomer_score()) == expected_score
assert parts(source.tautomer_score(api.TautomerScoreParams(terms=[]))) == empty_score
assert parts(source.tautomer_score(terms=[])) == empty_score
assert parts(source.tautomer_score_with_params(api.TautomerScoreParams(terms=[]))) == empty_score
assert bytes(source.to_binary()) == bytes(before)
try: source.add_hydrogens_(only_on_atoms=[2**63])
except api.OperationError as error:
    assert error.hydrogen is not None
    assert error.cause is error.__cause__
else: raise AssertionError('native out-of-range hydrogen input accepted')
assert bytes(source.to_binary()) == bytes(before)
"#),Some(&locals),None)?;
 let potential=locals.get_item("potential")?.unwrap().extract::<Py<crate::stereo_values::PotentialStereoResult>>()?;
 assert!(potential.borrow(py).inner==expected_potential);
 Ok(())
 }).unwrap();
    }
}

pub(crate) fn sanitize(
    operations: Option<crate::sanitize_values::SanitizeOperationsInput>,
) -> PyResult<crate::sanitize_values::SanitizeParams> {
    let defaults = ck::SanitizeParams::default();
    Ok(crate::sanitize_values::SanitizeParams {
        inner: ck::SanitizeParams {
            operations: operations
                .map(|v| v.into_rust())
                .transpose()?
                .unwrap_or(defaults.operations),
        },
    })
}
pub(crate) fn kekulize(
    mark_atoms_bonds: Option<bool>,
    canonical: Option<bool>,
    max_backtracks: Option<u32>,
) -> PyResult<crate::params::KekulizeParams> {
    let defaults = ck::KekulizeParams::default();
    Ok(crate::params::KekulizeParams {
        inner: ck::KekulizeParams {
            mark_atoms_bonds: mark_atoms_bonds.unwrap_or(defaults.mark_atoms_bonds),
            canonical: canonical.unwrap_or(defaults.canonical),
            max_backtracks: max_backtracks.unwrap_or(defaults.max_backtracks),
        },
    })
}
pub(crate) fn aromaticity(
    model: Option<crate::selectors::AromaticityModelInput>,
) -> PyResult<crate::selectors::AromaticityParams> {
    let defaults = ck::AromaticityParams::default();
    Ok(crate::selectors::AromaticityParams {
        inner: ck::AromaticityParams {
            model: model
                .map(|v| v.into_rust())
                .transpose()?
                .unwrap_or(defaults.model),
        },
    })
}
pub(crate) fn valence(
    model: Option<crate::selectors::ValenceModelInput>,
    strict: Option<bool>,
) -> PyResult<crate::selectors::ValenceParams> {
    let defaults = ck::ValenceParams::default();
    Ok(crate::selectors::ValenceParams {
        inner: ck::ValenceParams {
            model: model
                .map(|v| v.into_rust())
                .transpose()?
                .unwrap_or(defaults.model),
            strict: strict.unwrap_or(defaults.strict),
        },
    })
}
pub(crate) fn rings(
    include_dative_bonds: Option<bool>,
    include_hydrogen_bonds: Option<bool>,
) -> PyResult<crate::params::RingSearchParams> {
    let defaults = ck::RingSearchParams::default();
    Ok(crate::params::RingSearchParams {
        inner: ck::RingSearchParams {
            include_dative_bonds: include_dative_bonds.unwrap_or(defaults.include_dative_bonds),
            include_hydrogen_bonds: include_hydrogen_bonds
                .unwrap_or(defaults.include_hydrogen_bonds),
        },
    })
}
pub(crate) fn structure(
    conformer_id: Option<i32>,
    replace_existing_tags: Option<bool>,
) -> PyResult<crate::params::StructureTagParams> {
    let defaults = ck::StructureTagParams::default();
    Ok(crate::params::StructureTagParams {
        inner: ck::StructureTagParams {
            conformer_id: conformer_id.unwrap_or(defaults.conformer_id),
            replace_existing_tags: replace_existing_tags.unwrap_or(defaults.replace_existing_tags),
        },
    })
}
pub(crate) fn potential(
    clean: Option<bool>,
    flag_possible: Option<bool>,
    allow_nontetrahedral: Option<bool>,
) -> PyResult<crate::params::PotentialStereoParams> {
    let defaults = ck::PotentialStereoParams::default();
    Ok(crate::params::PotentialStereoParams {
        inner: ck::PotentialStereoParams {
            clean: clean.unwrap_or(defaults.clean),
            flag_possible: flag_possible.unwrap_or(defaults.flag_possible),
            allow_nontetrahedral: allow_nontetrahedral.unwrap_or(defaults.allow_nontetrahedral),
        },
    })
}
pub(crate) fn hydrogens(
    explicit_only: Option<bool>,
    add_coords: Option<bool>,
    add_residue_info: Option<bool>,
    skip_queries: Option<bool>,
    only_on_atoms: Option<Vec<usize>>,
) -> PyResult<ck::AddHsParams> {
    let defaults = ck::AddHsParams::default();
    Ok(ck::AddHsParams {
        explicit_only: explicit_only.unwrap_or(defaults.explicit_only),
        add_coords: add_coords.unwrap_or(defaults.add_coords),
        add_residue_info: add_residue_info.unwrap_or(defaults.add_residue_info),
        skip_queries: skip_queries.unwrap_or(defaults.skip_queries),
        only_on_atoms: only_on_atoms
            .map(|atoms| atoms.into_iter().map(ck::AtomId::new).collect())
            .or(defaults.only_on_atoms),
    })
}
pub(crate) fn remove_hydrogens(
    remove_degree_zero: Option<bool>,
    remove_higher_degrees: Option<bool>,
    remove_only_h_neighbors: Option<bool>,
    remove_isotopes: Option<bool>,
    remove_and_track_isotopes: Option<bool>,
    remove_dummy_neighbors: Option<bool>,
    remove_defining_bond_stereo: Option<bool>,
    remove_with_wedged_bond: Option<bool>,
    remove_with_query: Option<bool>,
    remove_mapped: Option<bool>,
    remove_in_sgroups: Option<bool>,
    show_warnings: Option<bool>,
    remove_nonimplicit: Option<bool>,
    update_explicit_count: Option<bool>,
    remove_hydrides: Option<bool>,
    remove_nontetrahedral_neighbors: Option<bool>,
    sanitize: Option<bool>,
) -> PyResult<crate::params::RemoveHsParams> {
    let defaults = ck::RemoveHsParams::default();
    Ok(crate::params::RemoveHsParams {
        inner: ck::RemoveHsParams {
            remove_degree_zero: remove_degree_zero.unwrap_or(defaults.remove_degree_zero),
            remove_higher_degrees: remove_higher_degrees.unwrap_or(defaults.remove_higher_degrees),
            remove_only_h_neighbors: remove_only_h_neighbors
                .unwrap_or(defaults.remove_only_h_neighbors),
            remove_isotopes: remove_isotopes.unwrap_or(defaults.remove_isotopes),
            remove_and_track_isotopes: remove_and_track_isotopes
                .unwrap_or(defaults.remove_and_track_isotopes),
            remove_dummy_neighbors: remove_dummy_neighbors
                .unwrap_or(defaults.remove_dummy_neighbors),
            remove_defining_bond_stereo: remove_defining_bond_stereo
                .unwrap_or(defaults.remove_defining_bond_stereo),
            remove_with_wedged_bond: remove_with_wedged_bond
                .unwrap_or(defaults.remove_with_wedged_bond),
            remove_with_query: remove_with_query.unwrap_or(defaults.remove_with_query),
            remove_mapped: remove_mapped.unwrap_or(defaults.remove_mapped),
            remove_in_sgroups: remove_in_sgroups.unwrap_or(defaults.remove_in_sgroups),
            show_warnings: show_warnings.unwrap_or(defaults.show_warnings),
            remove_nonimplicit: remove_nonimplicit.unwrap_or(defaults.remove_nonimplicit),
            update_explicit_count: update_explicit_count.unwrap_or(defaults.update_explicit_count),
            remove_hydrides: remove_hydrides.unwrap_or(defaults.remove_hydrides),
            remove_nontetrahedral_neighbors: remove_nontetrahedral_neighbors
                .unwrap_or(defaults.remove_nontetrahedral_neighbors),
            sanitize: sanitize.unwrap_or(defaults.sanitize),
        },
    })
}
pub(crate) fn cip(
    atoms: Option<Vec<usize>>,
    bonds: Option<Vec<usize>>,
    max_recursive_iterations: Option<u32>,
) -> PyResult<ck::CipLabelOptions> {
    let defaults = ck::CipLabelOptions::default();
    let mut inner = defaults.with_max_recursive_iterations(
        max_recursive_iterations
            .unwrap_or_else(|| ck::CipLabelOptions::default().max_recursive_iterations()),
    );
    if let Some(atoms) = atoms {
        inner = inner.with_atoms(atoms.into_iter().map(ck::AtomId::new));
    }
    if let Some(bonds) = bonds {
        inner = inner.with_bonds(bonds.into_iter().map(ck::BondId::new));
    }
    Ok(inner)
}
pub(crate) fn tautomer_score(
    py: Python<'_>,
    terms: Option<Vec<Py<crate::tautomer_values::TautomerScoreTerm>>>,
) -> PyResult<ck::TautomerScoreParams> {
    Ok(ck::TautomerScoreParams {
        terms: terms.map(|terms| {
            terms
                .iter()
                .map(|term| term.borrow(py).inner.clone())
                .collect()
        }),
    })
}
