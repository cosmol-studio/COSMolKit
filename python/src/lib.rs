//! Canonical Python projections of the public COSMolKit facade.
//!
//! Default and explicit drawing selections share one live Molecule class.
//! Historical adapters remain as source evidence, outside module compilation.

mod canonical_atom_bond;
mod canonical_batch;
mod canonical_batch_fingerprint_values;
mod canonical_batch_params;
mod canonical_binary;
mod canonical_bio_binding;
mod canonical_bio_metadata;
mod canonical_bio_residue;
mod canonical_bio_values;
mod canonical_builder;
mod canonical_chemistry_values;
mod canonical_coordinate_input;
mod canonical_descriptor_binding;
mod canonical_detached_blocks;
mod canonical_element_metadata;
mod canonical_error_accessors;
mod canonical_error_values;
mod canonical_fingerprint_values;
mod canonical_group_values;
mod canonical_inchi;
mod canonical_mcs;
mod canonical_molecular_hash;
mod canonical_operation_metadata;
mod canonical_potential_stereo;
mod canonical_reaction;
mod canonical_registered_errors;
mod canonical_sdf;
mod canonical_search;
mod canonical_smiles_writer;
mod canonical_stereo_queries;
mod canonical_valence;
mod canonical_values;
mod configuration_projection;
mod drawing_binding;
mod native_property;
mod text_input;
mod text_path;
mod user_path;

#[cfg(feature = "stubgen")]
pyo3_stub_gen::define_stub_info_gatherer!(stub_info);

mod mmff_binding;

mod tautomer_binding;

mod canonical_property_values;
mod uff_binding;

mod alignment_binding;
mod canonical_avalon;
mod canonical_layered;
mod canonical_maccs;
mod canonical_path_score;

mod canonical_pattern;
mod conformer_binding;

mod canonical_topological;

mod canonical_molecular_io;

mod canonical_sdf_supplier;

mod canonical_stereoisomers;

mod persistent_forcefields;

/// Construct the real extension module for the development contract gate.
#[cfg(feature = "stubgen")]
pub fn binding_contract_module(
    py: pyo3::Python<'_>,
) -> pyo3::PyResult<pyo3::Bound<'_, pyo3::types::PyModule>> {
    use pyo3::prelude::*;

    let module = pyo3::types::PyModule::new(py, "cosmolkit")?;
    drawing_binding::cosmolkit(&module)?;
    // Match normal extension loading: recursive imports in getters/setters must
    // resolve to this build, never to an installed package or a missing module.
    py.import("sys")?
        .getattr("modules")?
        .set_item("cosmolkit", &module)?;
    Ok(module)
}
mod rdkit_binding;

#[cfg(all(test, feature = "stubgen"))]
mod binding_contract_tests {
    use pyo3::{prelude::*, types::PyModule};

    #[test]
    fn contract_imports_use_current_module_without_an_install_or_with_a_stale_module() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let modules = py.import("sys")?.getattr("modules")?;
            let previous = modules.call_method1("get", ("cosmolkit",))?;
            let result = (|| -> PyResult<()> {
                for stale in [false, true] {
                    if stale {
                        modules.set_item("cosmolkit", PyModule::new(py, "cosmolkit")?)?;
                    } else {
                        modules.call_method1("pop", ("cosmolkit", py.None()))?;
                    }
                    let module = super::binding_contract_module(py)?;
                    assert!(py.import("cosmolkit")?.is(&module));
                    for (name, field) in [
                        ("BioReadParams", "format"),
                        ("BatchParams", "errors"),
                        ("BatchExportParams", "format"),
                        ("BatchImageParams", "format"),
                    ] {
                        let value = module.getattr(name)?.call0()?;
                        let before = value.getattr(field)?;
                        value.setattr(field, &before)?;
                        assert!(value.getattr(field)?.eq(&before)?);
                    }
                }
                Ok(())
            })();
            if previous.is_none() {
                modules.call_method1("pop", ("cosmolkit", py.None()))?;
            } else {
                modules.set_item("cosmolkit", previous)?;
            }
            result
        })
        .expect("contract checks must use the current in-memory extension");
    }
}
