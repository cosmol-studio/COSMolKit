//! XYZ/MOL2 transport delegates all parsing, finalization and formatting.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn from_xyz_block(text: &str) -> Result<Self, ck::MolecularIoError> {
        // COSMolKit❗✔️: ck::Molecule::from_xyz_block(text)
        ck::Molecule::from_xyz_block(text).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn read_xyz(text: &str) -> Result<Self, ck::MolecularIoError> {
        // COSMolKit❗✔️: ck::Molecule::read_xyz(text)
        ck::Molecule::read_xyz(text).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn from_mol2(text: &str) -> Result<Self, ck::MolecularIoError> {
        // COSMolKit❗✔️: ck::Molecule::from_mol2(text)
        ck::Molecule::from_mol2(text).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn from_mol2_with_params(
        text: &str,
        params: &ck::Mol2ReadParams,
    ) -> Result<Self, ck::MolecularIoError> {
        // COSMolKit❗✔️: ck::Molecule::from_mol2_with_params(text,params)
        ck::Molecule::from_mol2_with_params(text, params).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn read_mol2(text: &str) -> Result<Self, ck::MolecularIoError> {
        // COSMolKit❗✔️: ck::Molecule::read_mol2(text)
        ck::Molecule::read_mol2(text).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn read_mol2_with_params(
        text: &str,
        params: &ck::Mol2ReadParams,
    ) -> Result<Self, ck::MolecularIoError> {
        // COSMolKit❗✔️: ck::Molecule::read_mol2_with_params(text,params)
        ck::Molecule::read_mol2_with_params(text, params).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn to_xyz(&self) -> Result<String, ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.to_xyz()
        self.inner.borrow().to_xyz()
    }
    pub fn write_xyz(&self, path: &str) -> Result<(), ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.write_xyz(path)
        self.inner.borrow().write_xyz(path)
    }
    pub fn to_xyz_with_params(
        &self,
        params: &ck::XyzWriteParams,
    ) -> Result<String, ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.to_xyz_with_params(params)
        self.inner.borrow().to_xyz_with_params(params)
    }
    pub fn write_xyz_with_params(
        &self,
        path: &str,
        params: &ck::XyzWriteParams,
    ) -> Result<(), ck::MolecularIoError> {
        // COSMolKit❗✔️: self.inner.write_xyz_with_params(path,params)
        self.inner.borrow().write_xyz_with_params(path, params)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::error::Error;
    const XYZ: &str = "2\nXYZ comment\nC 0.125 -1.5 2\nH 1 0 0\n";
    const MOL2: &str = "@<TRIPOS>MOLECULE\nbinding mol2\n3 2 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n1 C1 0.0 0.0 0.0 C.3 1 UNL 0.0\n2 O1 1.4 0.0 0.0 O.3 1 UNL 0.0\n3 H1 2.2 0.0 0.0 H 1 UNL 0.0\n@<TRIPOS>BOND\n1 1 2 1\n2 2 3 1\n";
    #[test]
    fn full_parameters_and_real_files_match_canonical_owner() {
        let m = Molecule::from_xyz_block(XYZ).unwrap();
        let source = ck::Molecule::from_xyz_block(XYZ).unwrap();
        assert_eq!(*m.inner.borrow(), source);
        assert_eq!(m.num_atoms(), 2);
        assert_eq!(m.num_bonds(), 0);
        let before = m.inner.borrow().clone();
        assert_eq!(m.to_xyz().unwrap(), source.to_xyz().unwrap());
        for precision in [0, 2, 6, 12] {
            let p = ck::XyzWriteParams {
                conformer_id: Some(0),
                precision,
            };
            assert_eq!(
                m.to_xyz_with_params(&p).unwrap(),
                source.to_xyz_with_params(&p).unwrap()
            );
        }
        let missing = ck::XyzWriteParams {
            conformer_id: Some(99),
            precision: 2,
        };
        let e = m.to_xyz_with_params(&missing).unwrap_err();
        assert!(matches!(
            e,
            ck::MolecularIoError::XyzWrite(ck::XyzWriteError::ConformerNotFound { id: 99 })
        ));
        assert_eq!(*m.inner.borrow(), before);
        let root = std::env::temp_dir().join(format!("cosmolkit-xyz-mol2-{}", std::process::id()));
        std::fs::create_dir(&root).unwrap();
        let xyz = root.join("input.xyz");
        let mol2 = root.join("input.mol2");
        let output = root.join("output.xyz");
        std::fs::write(&xyz, XYZ).unwrap();
        std::fs::write(&mol2, MOL2).unwrap();
        assert_eq!(
            *Molecule::read_xyz(xyz.to_str().unwrap())
                .unwrap()
                .inner
                .borrow(),
            source
        );
        m.write_xyz(output.to_str().unwrap()).unwrap();
        assert_eq!(
            std::fs::read_to_string(&output).unwrap(),
            source.to_xyz().unwrap()
        );
        let write = ck::XyzWriteParams {
            conformer_id: Some(0),
            precision: 2,
        };
        m.write_xyz_with_params(output.to_str().unwrap(), &write)
            .unwrap();
        assert_eq!(
            std::fs::read_to_string(&output).unwrap(),
            source.to_xyz_with_params(&write).unwrap()
        );
        assert_eq!(
            *Molecule::from_mol2(MOL2).unwrap().inner.borrow(),
            ck::Molecule::from_mol2(MOL2).unwrap()
        );
        assert_eq!(
            *Molecule::read_mol2(mol2.to_str().unwrap())
                .unwrap()
                .inner
                .borrow(),
            ck::Molecule::read_mol2(mol2.to_str().unwrap()).unwrap()
        );
        for sanitize in [false, true] {
            for remove_hydrogens in [false, true] {
                for cleanup_substructures in [false, true] {
                    let p = ck::Mol2ReadParams {
                        sanitize,
                        remove_hydrogens,
                        cleanup_substructures,
                        variant: ck::Mol2Type::Corina,
                    };
                    let expected = ck::Molecule::from_mol2_with_params(MOL2, &p).unwrap();
                    assert_eq!(
                        *Molecule::from_mol2_with_params(MOL2, &p)
                            .unwrap()
                            .inner
                            .borrow(),
                        expected
                    );
                    assert_eq!(
                        *Molecule::read_mol2_with_params(mol2.to_str().unwrap(), &p)
                            .unwrap()
                            .inner
                            .borrow(),
                        expected
                    );
                    assert_eq!(
                        expected.num_atoms(),
                        if sanitize && remove_hydrogens { 2 } else { 3 }
                    );
                }
            }
        }
        for path in [xyz, mol2, output] {
            std::fs::remove_file(path).unwrap();
        }
        std::fs::remove_dir(root).unwrap();
    }
    #[test]
    fn structured_failures_and_source_priority_are_retained() {
        for (text, kind) in [
            ("", "EmptyBlock"),
            ("bad\n", "AtomCount"),
            ("1\ncomment\n", "UnexpectedEof"),
            ("1\n\nC 0 0\n", "MissingCoordinates"),
            ("1\n\nZz bad 0 0\n", "Coordinate"),
            ("1\n\nZz 0 0 0\n", "AtomSymbol"),
        ] {
            let e = Molecule::from_xyz_block(text).err().unwrap();
            let owner = ck::Molecule::from_xyz_block(text).unwrap_err();
            assert_eq!(e.to_string(), owner.to_string());
            let detail = e.source().unwrap();
            assert!(format!("{detail:?}").starts_with(kind));
            if kind == "Coordinate" {
                let e = detail.downcast_ref::<ck::XyzReadError>().unwrap();
                assert!(matches!(e,ck::XyzReadError::Coordinate{value,line:2,..} if value=="bad"));
                assert!(e.source().is_some());
            }
        }
        assert!(
            matches!(Molecule::from_mol2("").err().unwrap(),ck::MolecularIoError::Mol2Read(ck::Mol2ReadError::Parse(message)) if message=="No MOLECULE block found in Mol2 data")
        );
        let unsupported = MOL2.replace("C.3", "ANY");
        let e = Molecule::from_mol2(&unsupported).err().unwrap();
        assert!(matches!(
            e,
            ck::MolecularIoError::Mol2Read(ck::Mol2ReadError::Unsupported { .. })
        ));
        let parse = MOL2.replace("C.3", "Zz");
        assert!(matches!(
            Molecule::from_mol2(&parse).err().unwrap(),
            ck::MolecularIoError::Mol2Read(ck::Mol2ReadError::Parse(_))
        ));
        let bad = format!(
            "@<TRIPOS>MOLECULE\nhypervalent\n6 5 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3\n{}@<TRIPOS>BOND\n{}",
            (0..5)
                .map(|i| format!("{} F{} {} 0 0 F\n", i + 2, i + 1, i + 1))
                .collect::<String>(),
            (0..5)
                .map(|i| format!("{} 1 {} 1\n", i + 1, i + 2))
                .collect::<String>()
        );
        let error = Molecule::from_mol2(&bad).err().unwrap();
        let owner = ck::Molecule::from_mol2(&bad).unwrap_err();
        assert_eq!(error.to_string(), owner.to_string());
        assert!(matches!(&error,ck::MolecularIoError::Mol2Post(e) if e.stage=="Sanitize"));
        let post = error
            .source()
            .unwrap()
            .downcast_ref::<ck::Mol2PostError>()
            .unwrap();
        assert!(post.source().is_some());
        let empty = Molecule::from_xyz_block("0\n\n").unwrap();
        assert_eq!(
            empty
                .to_xyz_with_params(&ck::XyzWriteParams {
                    conformer_id: Some(99),
                    precision: 6
                })
                .unwrap(),
            ""
        );
    }
}
