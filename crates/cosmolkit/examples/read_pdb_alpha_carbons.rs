use std::{env, fs};

use cosmolkit::{BioStructure, ResidueKind};

const DEMO_PDB: &str = "\
ATOM      1  N   MET A   1      11.104  13.207   9.900  1.00 20.00           N
ATOM      2  CA  MET A   1      12.210  13.912  10.555  1.00 20.00           C
ATOM      3  C   MET A   1      13.470  13.079  10.413  1.00 20.00           C
ATOM      4  N   MSE A   2      14.530  13.650  10.980  1.00 20.00           N
ATOM      5  CA  MSE A   2      15.790  12.920  10.910  1.00 20.00           C
ATOM      6  C   MSE A   2      16.720  13.340   9.770  1.00 20.00           C
HETATM    7  O   HOH A   3      18.000  10.000   8.000  1.00 10.00           O
";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let pdb = match env::args().nth(1) {
        Some(path) => fs::read_to_string(path)?,
        None => DEMO_PDB.to_string(),
    };
    let structure = BioStructure::from_pdb_str(&pdb)?;
    let mmcif = structure.to_mmcif()?;
    let roundtrip = BioStructure::from_mmcif_str(&mmcif, "converted.cif")?;

    for alpha_carbon in alpha_carbons(&roundtrip) {
        println!(
            "{} {} {:.3} {:.3} {:.3}",
            alpha_carbon.residue_name,
            alpha_carbon.residue_index,
            alpha_carbon.position[0],
            alpha_carbon.position[1],
            alpha_carbon.position[2],
        );
    }

    Ok(())
}

struct AlphaCarbon<'a> {
    residue_index: usize,
    residue_name: &'a str,
    position: [f64; 3],
}

fn alpha_carbons(structure: &BioStructure) -> impl Iterator<Item = AlphaCarbon<'_>> {
    structure
        .residues()
        .iter()
        .enumerate()
        .filter(|(_, residue)| residue.kind == ResidueKind::AminoAcid)
        .filter_map(|(residue_index, residue)| {
            let start = residue.atom_span.start as usize;
            let end = residue.atom_span.end() as usize;
            (start..end)
                .find(|&atom_index| {
                    matches!(
                        structure.atoms()[atom_index].name.0,
                        [b' ', b'C', b'A', b' '] | [b'C', b'A', b' ', b' ']
                    )
                })
                .map(|atom_index| AlphaCarbon {
                    residue_index,
                    residue_name: residue.name.as_str(),
                    position: structure.coordinates().positions()[atom_index],
                })
        })
}
