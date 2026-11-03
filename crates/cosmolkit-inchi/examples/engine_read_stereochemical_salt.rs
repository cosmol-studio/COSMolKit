//! Read a polycyclic, seven-stereocenter alkaloid salt as official records.
use std::ffi::CStr;
use std::io;

use cosmolkit_inchi::engine::{self, Input, ReadInput};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Same input as the owning crate's guanidinium alkaloid hydrochloride
    // regression. No neutral graph, Toolkit or RDKit adapter is constructed.
    let inchi = c"InChI=1S/C22H33N3O4.ClH/c1-3-16-8-4-5-11-21(29-16)13-15-9-10-17-18(19(26)27)22(12-6-7-14(2)28-22)24-20(23-21)25(15)17;/h4,8,14-18H,3,5-7,9-13H2,1-2H3,(H2,23,24,26,27);1H/t14-,15+,16+,17-,18+,21+,22-;/m1./s1";
    let structure = engine::read(ReadInput {
        inchi,
        options: None,
    })?;
    println!("Read status: {}", structure.return_code);
    println!("Warning flags: {:?}", structure.warning_flags);
    for (name, bytes) in [("Message", &structure.message), ("Log", &structure.log)] {
        if let Some(bytes) = bytes {
            println!("{name}: {}", String::from_utf8_lossy(bytes));
        }
    }
    if !matches!(structure.return_code, 0 | 1) {
        return Err(io::Error::other("GetStructFromINCHI did not return a structure").into());
    }

    println!("\n{} raw atoms:", structure.atoms.len());
    for (index, atom) in structure.atoms.iter().enumerate() {
        let symbol = CStr::from_bytes_until_nul(bytemuck::cast_slice(&atom.elname))?;
        let count = usize::try_from(atom.num_bonds)?;
        println!(
            "{index}: {} charge={} isotope={} H={:?} neighbors={:?} bond_types={:?}",
            symbol.to_str()?,
            atom.charge,
            atom.isotopic_mass,
            atom.num_iso_H,
            &atom.neighbor[..count],
            &atom.bond_type[..count],
        );
    }
    println!("\n{} raw stereo records:", structure.stereo0d.len());
    for stereo in &structure.stereo0d {
        println!("{stereo:?}");
    }

    let generated = engine::generate(Input {
        atoms: &structure.atoms,
        stereo0d: &structure.stereo0d,
        options: None,
    })?;
    if !matches!(generated.return_code, 0 | 1) {
        return Err(io::Error::other(format!("Generation failed: {generated:?}")).into());
    }
    let regenerated = generated
        .inchi
        .ok_or_else(|| io::Error::other("GetINCHI returned no text"))?;
    assert_eq!(regenerated, inchi.to_bytes());
    println!(
        "\nExact InChI roundtrip: {}",
        String::from_utf8(regenerated)?
    );
    Ok(())
}
