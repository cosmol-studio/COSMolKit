//! Generate both configurations of penta-2,3-diene without any coordinates.
use std::ffi::CString;
use std::io;

use cosmolkit_inchi::engine::{self, Atom, Input, ReadInput, Stereo0D};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Atom indices: CH3(0)-CH(1)=C(2)=CH(3)-CH3(4).
    let mut atoms = vec![Atom::default(); 5];
    for (atom, hydrogens) in atoms.iter_mut().zip([3, 1, 0, 1, 3]) {
        atom.elname[0] = b'C' as i8;
        atom.num_iso_H[0] = hydrogens;
        // All coordinates stay zero so the official engine uses Stereo0D.
    }
    for (first, second, order) in [(0, 1, 1), (1, 2, 2), (2, 3, 2), (3, 4, 1)] {
        for (source, target) in [(first, second), (second, first)] {
            let atom = &mut atoms[source];
            let slot = usize::try_from(atom.num_bonds)?;
            atom.neighbor[slot] = i16::try_from(target)?;
            atom.bond_type[slot] = order;
            atom.num_bonds += 1;
        }
    }

    // Official inchi_api.h allene convention: X-A=C=B-Y,
    // neighbor = [X, A, B, Y], central_atom = C, type = Allene (3).
    // Parities 1/2 mean ODD/EVEN for this ordering, not absolute R/S labels.
    let mut identifiers = Vec::new();
    for parity in [1, 2] {
        let stereo = Stereo0D {
            neighbor: [0, 1, 3, 4],
            central_atom: 2,
            type_: 3,
            parity,
        };
        let generated = engine::generate(Input {
            atoms: &atoms,
            stereo0d: &[stereo],
            options: None,
        })?;
        if !matches!(generated.return_code, 0 | 1) {
            return Err(io::Error::other(format!("Generation failed: {generated:?}")).into());
        }
        let inchi = generated
            .inchi
            .ok_or_else(|| io::Error::other("GetINCHI returned no text"))?;
        println!("parity {parity}: {}", std::str::from_utf8(&inchi)?);
        let text = CString::new(inchi.clone())?;
        let structure = engine::read(ReadInput {
            inchi: &text,
            options: None,
        })?;
        if !matches!(structure.return_code, 0 | 1) {
            return Err(io::Error::other(format!("Reading failed: {structure:?}")).into());
        }
        assert_eq!(structure.stereo0d.len(), 1);
        assert_eq!(structure.stereo0d[0].type_, 3);
        println!("Read raw allene: {:?}", structure.stereo0d[0]);
        let regenerated = engine::generate(Input {
            atoms: &structure.atoms,
            stereo0d: &structure.stereo0d,
            options: None,
        })?;
        if !matches!(regenerated.return_code, 0 | 1) {
            return Err(io::Error::other(format!("Regeneration failed: {regenerated:?}")).into());
        }
        assert_eq!(regenerated.inchi.as_deref(), Some(inchi.as_slice()));
        identifiers.push(inchi);
    }
    assert_ne!(identifiers[0], identifiers[1]);
    println!("Opposite allene configurations stay distinct through read/regenerate.");
    Ok(())
}
