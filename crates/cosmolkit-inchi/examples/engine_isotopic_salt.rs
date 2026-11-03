//! Construct [13CH3]C(=O)[O-].[Na+] directly with official InChI records.
use std::ffi::CString;
use std::io;

use cosmolkit_inchi::engine::{self, Atom, Input, ReadInput};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut atoms = vec![Atom::default(); 5];
    for (atom, symbol) in atoms
        .iter_mut()
        .zip([b"C".as_slice(), b"C", b"O", b"O", b"Na"])
    {
        for (slot, byte) in atom.elname.iter_mut().zip(symbol) {
            *slot = *byte as i8;
        }
    }
    atoms[0].num_iso_H[0] = 3;
    atoms[0].isotopic_mass = 13; // Absolute mass, an official input encoding.
    atoms[3].charge = -1;
    atoms[4].charge = 1;

    // Only acetate has bonds. Na is a separate connected component.
    // Populate both adjacency lists; 1 = single, 2 = double.
    for (first, second, order) in [(0, 1, 1), (1, 2, 2), (1, 3, 1)] {
        for (source, target) in [(first, second), (second, first)] {
            let atom = &mut atoms[source];
            let slot = usize::try_from(atom.num_bonds)?;
            atom.neighbor[slot] = i16::try_from(target)?;
            atom.bond_type[slot] = order;
            atom.num_bonds += 1;
        }
    }
    let generated = engine::generate(Input {
        atoms: &atoms,
        stereo0d: &[],
        options: None,
    })?;
    if !matches!(generated.return_code, 0 | 1) {
        return Err(io::Error::other(format!("Generation failed: {generated:?}")).into());
    }
    let inchi = generated
        .inchi
        .ok_or_else(|| io::Error::other("GetINCHI returned no text"))?;
    println!("{}", std::str::from_utf8(&inchi)?);
    if let Some(aux_info) = &generated.aux_info {
        println!("{}", std::str::from_utf8(aux_info)?);
    }

    let text = CString::new(inchi.clone())?;
    let structure = engine::read(ReadInput {
        inchi: &text,
        options: None,
    })?;
    if !matches!(structure.return_code, 0 | 1) {
        return Err(io::Error::other(format!("Reading failed: {structure:?}")).into());
    }
    // Reading uses canonical atom order; do not assume the input indices.
    // Returned isotope masses may use ISOTOPIC_SHIFT_FLAG + mass difference,
    // which is also accepted directly by the generator.
    for (index, atom) in structure.atoms.iter().enumerate() {
        println!(
            "atom {index}: charge={}, raw isotope={}, H={:?}",
            atom.charge, atom.isotopic_mass, atom.num_iso_H
        );
    }
    let regenerated = engine::generate(Input {
        atoms: &structure.atoms,
        stereo0d: &structure.stereo0d,
        options: None,
    })?;
    if !matches!(regenerated.return_code, 0 | 1) {
        return Err(io::Error::other(format!("Regeneration failed: {regenerated:?}")).into());
    }
    assert_eq!(regenerated.inchi.as_deref(), Some(inchi.as_slice()));
    println!("Exact isotope/charge/component roundtrip.");
    Ok(())
}
