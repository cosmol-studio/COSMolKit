//! File format I/O entry points for chemistry and biomolecular data.

pub mod molblock;
pub mod sdf;

/// Returns versions of core modules to ensure linkage works.
#[must_use]
pub fn dependency_versions() -> (&'static str, &'static str) {
    (crate::version(), crate::bio::version())
}

#[cfg(test)]
mod tests {
    use super::{
        dependency_versions,
        molblock::{self, SdfFormat},
        sdf::SdfReader,
    };
    use crate::{BondOrder, Molecule};
    use std::io::Cursor;

    #[test]
    fn dependencies_are_available() {
        let (core, bio) = dependency_versions();
        assert!(!core.is_empty());
        assert!(!bio.is_empty());
    }

    #[test]
    fn sdf_reader_returns_none_for_empty_stream() {
        let mut reader = SdfReader::new(Cursor::new(Vec::<u8>::new()));
        let record = reader
            .next_record()
            .expect("empty SDF stream should be readable");
        assert_eq!(record, None);
    }

    #[test]
    fn sdf_reader_reads_minimal_v2000_topology_record() {
        let mut mol = Molecule::from_smiles("CC").expect("SMILES parser should parse CC");
        mol.compute_2d_coords().expect("2D coords should compute");
        let sdf = molblock::mol_to_sdf_record_minimal(&mol).expect("writer should work");

        let mut reader = SdfReader::new(Cursor::new(sdf.into_bytes()));
        let record = reader
            .next_record()
            .expect("minimal V2000 SDF record should parse")
            .expect("record should exist");

        assert_eq!(record.molecule.atoms.len(), 2);
        assert_eq!(record.molecule.bonds.len(), 1);
        assert_eq!(record.molecule.atomic_numbers(), vec![6, 6]);
        assert_eq!(record.molecule.bonds[0].order, BondOrder::Single);
        assert_eq!(record.data_fields, Vec::<(String, String)>::new());
        assert!(
            reader
                .next_record()
                .expect("second read should reach EOF")
                .is_none()
        );
    }

    #[test]
    fn sdf_reader_reads_rdkit_style_data_fields() {
        let sdf = b"ethane
     COSMolKit      2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
> <NAME>
ethane

> <NOTE>
line1
 continued after leading space

$$$$
";
        let mut reader = SdfReader::new(Cursor::new(sdf.to_vec()));
        let record = reader
            .next_record()
            .expect("SDF data fields should parse")
            .expect("record should exist");

        assert_eq!(
            record.data_fields,
            vec![
                ("NAME".to_owned(), "ethane".to_owned()),
                (
                    "NOTE".to_owned(),
                    "line1\n continued after leading space".to_owned()
                ),
            ]
        );
    }

    #[test]
    fn sdf_reader_roundtrips_v2000_2d_records_written_by_cosmolkit() {
        for smiles in [
            "CC",
            "C=C",
            "C#N",
            "C~C",
            "[Na+].[Cl-]",
            "[NH4+]",
            "[O-][N+](=O)O",
            "[13CH3:7][C@H](F)Cl",
            "F[C@](Cl)(Br)I",
            "CN1CCCC1",
        ] {
            let mut mol = Molecule::from_smiles(smiles).expect("SMILES parser should parse");
            mol.compute_2d_coords().expect("2D coords should compute");
            let sdf = molblock::mol_to_2d_sdf_record(&mol, SdfFormat::V2000).expect("V2000 write");

            let mut reader = SdfReader::new(Cursor::new(sdf.into_bytes()));
            let record = reader
                .next_record()
                .expect("written V2000 SDF should parse")
                .expect("record should exist");
            assert!(
                record.molecule.coords_2d().is_some(),
                "read molecule should retain 2D coordinates for {smiles}"
            );

            let rewritten = molblock::mol_to_2d_sdf_record(&record.molecule, SdfFormat::V2000)
                .expect("rewritten V2000 should write");
            let mut rereader = SdfReader::new(Cursor::new(rewritten.into_bytes()));
            let rerecord = rereader
                .next_record()
                .expect("rewritten V2000 SDF should parse")
                .expect("rewritten record should exist");

            assert_molecule_graph_and_coords_equal(&record.molecule, &rerecord.molecule, smiles);
        }
    }

    #[test]
    fn sdf_reader_roundtrips_v3000_2d_records_written_by_cosmolkit() {
        for smiles in [
            "CC",
            "C=C",
            "C#N",
            "C~C",
            "[Na+].[Cl-]",
            "[NH4+]",
            "[O-][N+](=O)O",
            "[13CH3:7][C@H](F)Cl",
            "F[C@](Cl)(Br)I",
            "c1ccccc1",
            "[NH3]->[Cu+2]<-[NH3]",
        ] {
            let mut mol = Molecule::from_smiles(smiles).expect("SMILES parser should parse");
            mol.compute_2d_coords().expect("2D coords should compute");
            let sdf = molblock::mol_to_2d_sdf_record(&mol, SdfFormat::V3000).expect("V3000 write");

            let mut reader = SdfReader::new(Cursor::new(sdf.into_bytes()));
            let record = reader
                .next_record()
                .expect("written V3000 SDF should parse")
                .expect("record should exist");
            assert!(
                record.molecule.coords_2d().is_some(),
                "read molecule should retain 2D coordinates for {smiles}"
            );

            let rewritten = molblock::mol_to_2d_sdf_record(&record.molecule, SdfFormat::V3000)
                .expect("rewritten V3000 should write");
            let mut rereader = SdfReader::new(Cursor::new(rewritten.into_bytes()));
            let rerecord = rereader
                .next_record()
                .expect("rewritten V3000 SDF should parse")
                .expect("rewritten record should exist");

            assert_molecule_graph_and_coords_equal(&record.molecule, &rerecord.molecule, smiles);
        }
    }

    fn assert_molecule_graph_and_coords_equal(lhs: &Molecule, rhs: &Molecule, label: &str) {
        assert_eq!(lhs.atoms.len(), rhs.atoms.len(), "atom count for {label}");
        assert_eq!(lhs.bonds.len(), rhs.bonds.len(), "bond count for {label}");

        for (idx, (lhs_atom, rhs_atom)) in lhs.atoms.iter().zip(rhs.atoms.iter()).enumerate() {
            assert_eq!(
                lhs_atom.atomic_num, rhs_atom.atomic_num,
                "atomic_num at atom {idx} for {label}"
            );
            assert_eq!(
                lhs_atom.formal_charge, rhs_atom.formal_charge,
                "formal_charge at atom {idx} for {label}"
            );
            assert_eq!(
                lhs_atom.num_radical_electrons, rhs_atom.num_radical_electrons,
                "radicals at atom {idx} for {label}"
            );
            assert_eq!(
                lhs_atom.is_aromatic, rhs_atom.is_aromatic,
                "is_aromatic at atom {idx} for {label}"
            );
            assert_eq!(
                lhs_atom.isotope, rhs_atom.isotope,
                "isotope at atom {idx} for {label}"
            );
        }

        for (idx, (lhs_bond, rhs_bond)) in lhs.bonds.iter().zip(rhs.bonds.iter()).enumerate() {
            assert_eq!(
                (lhs_bond.begin_atom, lhs_bond.end_atom),
                (rhs_bond.begin_atom, rhs_bond.end_atom),
                "bond endpoints at bond {idx} for {label}"
            );
            assert_eq!(
                lhs_bond.order, rhs_bond.order,
                "bond order at bond {idx} for {label}"
            );
            assert_eq!(
                lhs_bond.direction, rhs_bond.direction,
                "bond direction at bond {idx} for {label}"
            );
            assert_eq!(
                lhs_bond.stereo, rhs_bond.stereo,
                "bond stereo at bond {idx} for {label}"
            );
        }

        let lhs_coords = lhs.coords_2d().expect("lhs coords should exist");
        let rhs_coords = rhs.coords_2d().expect("rhs coords should exist");
        assert_eq!(
            lhs_coords.len(),
            rhs_coords.len(),
            "coord count for {label}"
        );
        for (idx, (lhs_coord, rhs_coord)) in lhs_coords.iter().zip(rhs_coords.iter()).enumerate() {
            assert!(
                (lhs_coord.x - rhs_coord.x).abs() <= 1e-12
                    && (lhs_coord.y - rhs_coord.y).abs() <= 1e-12,
                "coord mismatch at atom {} for {}: lhs=({:.6},{:.6}) rhs=({:.6},{:.6})",
                idx,
                label,
                lhs_coord.x,
                lhs_coord.y,
                rhs_coord.x,
                rhs_coord.y
            );
        }
    }

    #[test]
    fn molblock_minimal_writer_emits_v2000_for_ethane() {
        let mut mol = Molecule::from_smiles("CC").expect("SMILES parser should parse CC");
        mol.compute_2d_coords().expect("2D coords should compute");
        let out = molblock::mol_to_v2000_block_minimal(&mol).expect("writer should work");
        assert!(out.contains("V2000"));
    }
}
