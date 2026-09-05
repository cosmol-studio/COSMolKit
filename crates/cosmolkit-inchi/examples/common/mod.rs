use cosmolkit_inchi::{InchiMolecule, InchiToMolToolkit, InchiToolkitError, MolToInchiToolkit};

/// Exact adapter for the methane-only state used by the neutral-graph examples.
pub struct MethaneToolkit;

fn toolkit_error(call: &'static str, detail: &str) -> InchiToolkitError {
    InchiToolkitError {
        kind: "methane example boundary",
        message: format!("{call}: {detail}"),
    }
}

fn require_methane(molecule: &InchiMolecule, call: &'static str) -> Result<(), InchiToolkitError> {
    if molecule.atoms().len() != 1 || molecule.atoms()[0].atomic_number != 6 || !molecule.bonds().is_empty() {
        return Err(toolkit_error(
            call,
            "this example adapter accepts only one carbon atom with no bonds",
        ));
    }
    Ok(())
}

fn unsupported(call: &'static str) -> InchiToolkitError {
    toolkit_error(call, "operation is outside the methane example boundary")
}

impl MolToInchiToolkit for MethaneToolkit {
    fn needs_update_property_cache(&mut self, molecule: &InchiMolecule) -> Result<bool, InchiToolkitError> {
        require_methane(molecule, "needs_update_property_cache")?;
        Ok(false)
    }

    fn update_property_cache(&mut self, molecule: &mut InchiMolecule, strict: bool) -> Result<(), InchiToolkitError> {
        require_methane(molecule, "update_property_cache")?;
        if strict {
            return Err(unsupported("update_property_cache(strict=true)"));
        }
        Ok(())
    }

    fn kekulize(&mut self, molecule: &mut InchiMolecule, mark_atoms_bonds: bool) -> Result<(), InchiToolkitError> {
        require_methane(molecule, "kekulize")?;
        if mark_atoms_bonds {
            return Err(unsupported("kekulize(mark_atoms_bonds=true)"));
        }
        Ok(())
    }

    fn element_symbol(&mut self, atomic_number: i32) -> Result<Vec<u8>, InchiToolkitError> {
        match atomic_number {
            6 => Ok(b"C".to_vec()),
            _ => Err(unsupported("element_symbol")),
        }
    }

    fn atomic_weight(&mut self, atomic_number: i32) -> Result<f64, InchiToolkitError> {
        match atomic_number {
            6 => Ok(12.011),
            _ => Err(unsupported("atomic_weight")),
        }
    }

    fn total_num_hydrogens(&mut self, molecule: &InchiMolecule, atom_index: u32) -> Result<u32, InchiToolkitError> {
        require_methane(molecule, "total_num_hydrogens")?;
        molecule
            .atoms()
            .get(atom_index as usize)
            .map(|atom| atom.num_explicit_hydrogens)
            .ok_or_else(|| toolkit_error("total_num_hydrogens", "atom index is out of bounds"))
    }

    fn calc_implicit_valence(
        &mut self,
        _molecule: &mut InchiMolecule,
        _atom_index: u32,
    ) -> Result<i32, InchiToolkitError> {
        Err(unsupported("calc_implicit_valence"))
    }

    fn total_degree(&mut self, molecule: &InchiMolecule, atom_index: u32) -> Result<u32, InchiToolkitError> {
        require_methane(molecule, "total_degree")?;
        if atom_index == 0 {
            Ok(0)
        } else {
            Err(toolkit_error("total_degree", "atom index is out of bounds"))
        }
    }
}

impl InchiToMolToolkit for MethaneToolkit {
    fn atomic_number(&mut self, element: &[u8]) -> Result<i32, InchiToolkitError> {
        match element {
            b"C" => Ok(6),
            _ => Err(unsupported("atomic_number")),
        }
    }

    fn average_atomic_weight(&mut self, atomic_number: i32) -> Result<f64, InchiToolkitError> {
        match atomic_number {
            6 => Ok(12.011),
            _ => Err(unsupported("average_atomic_weight")),
        }
    }

    fn update_property_cache(&mut self, molecule: &mut InchiMolecule, strict: bool) -> Result<(), InchiToolkitError> {
        require_methane(molecule, "update_property_cache")?;
        if strict {
            return Err(unsupported("update_property_cache(strict=true)"));
        }
        Ok(())
    }

    fn assign_atom_cip_ranks(&mut self, _molecule: &mut InchiMolecule) -> Result<Vec<u32>, InchiToolkitError> {
        Err(unsupported("assign_atom_cip_ranks"))
    }

    fn remove_hydrogens(&mut self, _molecule: &mut InchiMolecule) -> Result<(), InchiToolkitError> {
        Err(unsupported("remove_hydrogens"))
    }

    fn sanitize_molecule(&mut self, _molecule: &mut InchiMolecule) -> Result<(), InchiToolkitError> {
        Err(unsupported("sanitize_molecule"))
    }

    fn assign_stereochemistry(
        &mut self,
        molecule: &mut InchiMolecule,
        clean_it: bool,
        force: bool,
    ) -> Result<(), InchiToolkitError> {
        require_methane(molecule, "assign_stereochemistry")?;
        if !clean_it || !force {
            return Err(unsupported("assign_stereochemistry"));
        }
        Ok(())
    }
}
