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
    use super::{dependency_versions, molblock, sdf::SdfReader};
    use crate::Molecule;
    use std::io::Cursor;

    #[test]
    fn dependencies_are_available() {
        let (core, bio) = dependency_versions();
        assert!(!core.is_empty());
        assert!(!bio.is_empty());
    }

    #[test]
    fn sdf_reader_type_exists() {
        let _reader = SdfReader::new(Cursor::new(Vec::<u8>::new()));
    }

    #[test]
    fn molblock_minimal_writer_exists() {
        let mol = Molecule::from_smiles("CC").expect("SMILES parser should parse CC");
        let out = molblock::mol_to_v2000_block_minimal(&mol).expect("writer should work");
        assert!(out.contains("V2000"));
    }
}
