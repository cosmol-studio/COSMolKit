//! Canonical MOL/SDF text and filesystem transports.
use crate::{Molecule, SdfRecord};
use cosmolkit as ck;
impl SdfRecord {
    pub fn from_query_graph(
        query: ck::QueryGraph,
        properties: ck::MoleculeProperties,
    ) -> Result<Self, ck::SdfError> {
        // COSMolKit❗✔️: ck::SdfRecord::from_query_graph(query,properties)
        ck::SdfRecord::from_query_graph(query, properties).map(|inner| Self { inner })
    }
    pub fn to_mol(&self) -> Result<String, ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.to_mol()
        self.inner.to_mol()
    }
    pub fn to_mol_with_params(
        &self,
        params: &ck::MolBlockWriteParams,
    ) -> Result<String, ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.to_mol_with_params(params)
        self.inner.to_mol_with_params(params)
    }
    pub fn to_sdf(&self) -> Result<String, ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.to_sdf()
        self.inner.to_sdf()
    }
    pub fn to_sdf_with_params(
        &self,
        params: &ck::MolBlockWriteParams,
    ) -> Result<String, ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.to_sdf_with_params(params)
        self.inner.to_sdf_with_params(params)
    }
}
impl Molecule {
    pub fn to_mol(&self) -> Result<String, ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.to_mol()
        self.inner.borrow().to_mol()
    }
    pub fn to_mol_with_params(
        &self,
        params: &ck::MolBlockWriteParams,
    ) -> Result<String, ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.to_mol_with_params(params)
        self.inner.borrow().to_mol_with_params(params)
    }
    pub fn to_sdf(&self) -> Result<String, ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.to_sdf()
        self.inner.borrow().to_sdf()
    }
    pub fn to_sdf_with_params(
        &self,
        params: &ck::MolBlockWriteParams,
    ) -> Result<String, ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.to_sdf_with_params(params)
        self.inner.borrow().to_sdf_with_params(params)
    }
    pub fn to_sdf_2d(&self) -> Result<String, ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.to_sdf_2d()
        self.inner.borrow().to_sdf_2d()
    }
    pub fn to_sdf_2d_with_params(
        &self,
        params: &ck::MolBlockWriteParams,
    ) -> Result<String, ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.to_sdf_2d_with_params(params)
        self.inner.borrow().to_sdf_2d_with_params(params)
    }
    pub fn to_sdf_3d(&self) -> Result<String, ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.to_sdf_3d()
        self.inner.borrow().to_sdf_3d()
    }
    pub fn to_sdf_3d_with_params(
        &self,
        params: &ck::MolBlockWriteParams,
    ) -> Result<String, ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.to_sdf_3d_with_params(params)
        self.inner.borrow().to_sdf_3d_with_params(params)
    }
    pub fn write_mol(&self, path: &str) -> Result<(), ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.write_mol(path)
        self.inner.borrow().write_mol(path)
    }
    pub fn write_mol_with_params(
        &self,
        path: &str,
        params: &ck::MolBlockWriteParams,
    ) -> Result<(), ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.write_mol_with_params(path,params)
        self.inner.borrow().write_mol_with_params(path, params)
    }
    pub fn write_sdf(&self, path: &str) -> Result<(), ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.write_sdf(path)
        self.inner.borrow().write_sdf(path)
    }
    pub fn write_sdf_with_params(
        &self,
        path: &str,
        params: &ck::MolBlockWriteParams,
    ) -> Result<(), ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.write_sdf_with_params(path,params)
        self.inner.borrow().write_sdf_with_params(path, params)
    }
    pub fn write_sdf_files(
        &self,
        directory: &str,
        file_name: Option<&str>,
    ) -> Result<std::path::PathBuf, ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.write_sdf_files(directory,file_name)
        self.inner.borrow().write_sdf_files(directory, file_name)
    }
    pub fn write_sdf_files_with_params(
        &self,
        directory: &str,
        file_name: Option<&str>,
        params: &ck::MolBlockWriteParams,
    ) -> Result<std::path::PathBuf, ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.write_sdf_files_with_params(directory,file_name,params)
        self.inner
            .borrow()
            .write_sdf_files_with_params(directory, file_name, params)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const SDF: &str = "writing\n  COSMolKit         2D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\nM  END\n>  <NOTE>\nfirst\n\n>  <NOTE>\nsecond\n\n$$$$\n";
    fn equal(a: Result<String, ck::MolecularIoError>, b: Result<String, ck::MolecularIoError>) {
        match (a, b) {
            (Ok(a), Ok(b)) => assert_eq!(a, b),
            (Err(a), Err(b)) => assert_eq!(a.to_string(), b.to_string()),
            _ => panic!("transport changed canonical result"),
        }
    }
    #[test]
    fn all_policies_and_query_records_delegate_exactly() {
        let m = Molecule::from_sdf(SDF).unwrap();
        let source = ck::Molecule::from_sdf(SDF).unwrap();
        let record = SdfRecord::from_sdf(SDF).unwrap();
        let owner = ck::SdfRecord::from_sdf(SDF).unwrap();
        let before = m.inner.borrow().clone();
        equal(m.to_mol(), source.to_mol());
        equal(m.to_sdf(), source.to_sdf());
        equal(m.to_sdf_2d(), source.to_sdf_2d());
        equal(m.to_sdf_3d(), source.to_sdf_3d());
        assert!(matches!(
            m.to_sdf_3d(),
            Err(ck::MolecularIoError::MolWrite(
                ck::MolWriteError::MissingCoordinate {
                    dimension: ck::CoordinateDimension::ThreeD,
                    id: 0
                }
            ))
        ));
        let spatial = Molecule::from_xyz_block("2\n\nC 0 0 0\nO 1.4 0 1\n").unwrap();
        assert!(spatial.to_sdf_3d().unwrap().contains("3D"));

        equal(record.to_mol(), owner.to_mol());
        equal(record.to_sdf(), owner.to_sdf());
        let q = ck::search::from_smarts("[#6,#7]-O").unwrap();
        let properties = ck::MoleculeProperties::default()
            .with_name("query")
            .with_sdf_data_field("NOTE", "one");
        let qr = SdfRecord::from_query_graph(q.clone(), properties.clone()).unwrap();
        let qo = ck::SdfRecord::from_query_graph(q, properties).unwrap();
        assert_eq!(qr.graph().kind(), "query_graph");
        equal(qr.to_mol(), qo.to_mol());
        equal(qr.to_sdf(), qo.to_sdf());
        for format in [ck::SdfFormat::V2000, ck::SdfFormat::V3000] {
            for force_2d in [false, true] {
                for include_stereo in [false, true] {
                    for kekulize in [false, true] {
                        for precision in [0, 2, 6] {
                            for include_coordinates in [false, true] {
                                let p = ck::MolBlockWriteParams {
                                    format,
                                    force_2d,
                                    include_stereo,
                                    kekulize,
                                    precision,
                                    coordinate_selection: ck::MolCoordinateSelection::Auto,
                                    include_coordinates,
                                };
                                equal(m.to_mol_with_params(&p), source.to_mol_with_params(&p));
                                equal(m.to_sdf_with_params(&p), source.to_sdf_with_params(&p));
                                equal(
                                    m.to_sdf_2d_with_params(&p),
                                    source.to_sdf_2d_with_params(&p),
                                );
                                equal(
                                    m.to_sdf_3d_with_params(&p),
                                    source.to_sdf_3d_with_params(&p),
                                );
                                equal(record.to_mol_with_params(&p), owner.to_mol_with_params(&p));
                                equal(record.to_sdf_with_params(&p), owner.to_sdf_with_params(&p));
                                equal(qr.to_mol_with_params(&p), qo.to_mol_with_params(&p));
                                equal(qr.to_sdf_with_params(&p), qo.to_sdf_with_params(&p));
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(*m.inner.borrow(), before);
        let missing = ck::MolBlockWriteParams {
            coordinate_selection: ck::MolCoordinateSelection::TwoD { id: 99 },
            ..Default::default()
        };
        assert!(matches!(
            m.to_mol_with_params(&missing),
            Err(ck::MolecularIoError::MolWrite(
                ck::MolWriteError::MissingCoordinate { id: 99, .. }
            ))
        ));
        let mismatch = ck::MolBlockWriteParams {
            coordinate_selection: ck::MolCoordinateSelection::ThreeD { id: 0 },
            ..Default::default()
        };
        assert!(matches!(
            m.to_sdf_2d_with_params(&mismatch),
            Err(ck::MolecularIoError::MolWrite(
                ck::MolWriteError::CoordinateDimensionMismatch
            ))
        ));
    }
    #[test]
    fn native_file_outputs_and_path_failures_are_source_ordered() {
        let m = Molecule::from_sdf(SDF).unwrap();
        let root = std::env::temp_dir().join(format!("cosmolkit-mol-sdf-{}", std::process::id()));
        std::fs::create_dir(&root).unwrap();
        let path = root.join("write.mol");
        let path = path.to_str().unwrap();
        let p = ck::MolBlockWriteParams {
            format: ck::SdfFormat::V3000,
            ..Default::default()
        };
        m.write_mol(path).unwrap();
        assert_eq!(std::fs::read_to_string(path).unwrap(), m.to_mol().unwrap());
        m.write_mol_with_params(path, &p).unwrap();
        assert_eq!(
            std::fs::read_to_string(path).unwrap(),
            m.to_mol_with_params(&p).unwrap()
        );
        m.write_sdf(path).unwrap();
        assert_eq!(std::fs::read_to_string(path).unwrap(), m.to_sdf().unwrap());
        m.write_sdf_with_params(path, &p).unwrap();
        assert_eq!(
            std::fs::read_to_string(path).unwrap(),
            m.to_sdf_with_params(&p).unwrap()
        );
        let named = m.write_sdf_files(root.to_str().unwrap(), None).unwrap();
        assert_eq!(named, root.join("molecule.sdf"));
        assert_eq!(
            std::fs::read_to_string(&named).unwrap(),
            m.to_sdf().unwrap()
        );
        let custom = m
            .write_sdf_files_with_params(root.to_str().unwrap(), Some("custom.sdf"), &p)
            .unwrap();
        assert_eq!(
            std::fs::read_to_string(&custom).unwrap(),
            m.to_sdf_with_params(&p).unwrap()
        );
        assert!(matches!(
            m.write_sdf_files(root.to_str().unwrap(), Some(" ")),
            Err(ck::MolecularIoError::Parameter {
                name: "file_name",
                ..
            })
        ));
        assert!(
            matches!(m.write_sdf_files(path,None),Err(ck::MolecularIoError::Io{source,..}) if source.kind()==std::io::ErrorKind::NotADirectory)
        );
        for p in [PathBuf::from(path), named, custom] {
            std::fs::remove_file(p).unwrap();
        }
        std::fs::remove_dir(root).unwrap();
    }
    use std::path::PathBuf;
}
