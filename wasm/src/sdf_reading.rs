//! Complete record and molecule reading transport through the canonical root.
use crate::Molecule;
use cosmolkit as ck;
#[derive(Clone, Debug)]
pub struct SdfRecord {
    pub(crate) inner: ck::SdfRecord,
}
#[derive(Clone, Debug)]
pub struct SdfGraph {
    pub(crate) inner: ck::SdfGraph,
}
impl SdfGraph {
    pub fn kind(&self) -> &str {
        match &self.inner {
            ck::SdfGraph::Molecule(_) => "molecule",
            ck::SdfGraph::Query(_) => "query_graph",
        }
    }
    pub fn molecule(&self) -> Option<Molecule> {
        match &self.inner {
            ck::SdfGraph::Molecule(inner) => Some(Molecule {
                inner: inner.clone().into(),
            }),
            ck::SdfGraph::Query(_) => None,
        }
    }
    pub fn query_graph(&self) -> Option<ck::QueryGraph> {
        match &self.inner {
            ck::SdfGraph::Query(inner) => Some(inner.clone()),
            ck::SdfGraph::Molecule(_) => None,
        }
    }
}
impl SdfRecord {
    pub fn from_sdf(text: &str) -> Result<Self, ck::SdfError> {
        // COSMolKit❗✔️: ck::SdfRecord::from_sdf(text)
        ck::SdfRecord::from_sdf(text).map(|inner| Self { inner })
    }
    pub fn from_sdf_with_params(
        text: &str,
        params: &ck::SdfReadParams,
    ) -> Result<Self, ck::SdfError> {
        // COSMolKit❗✔️: ck::SdfRecord::from_sdf_with_params(text,params)
        ck::SdfRecord::from_sdf_with_params(text, params).map(|inner| Self { inner })
    }
    pub fn graph(&self) -> SdfGraph {
        SdfGraph {
            inner: self.inner.graph().clone(),
        }
    }
    pub fn molecule(&self) -> Result<Molecule, ck::SdfError> {
        self.inner.molecule().map(|inner| Molecule {
            inner: inner.clone().into(),
        })
    }
    pub fn query_graph(&self) -> Result<ck::QueryGraph, ck::SdfError> {
        self.inner.query_graph().cloned()
    }
    pub fn properties(&self) -> ck::MoleculeProperties {
        self.inner.properties().clone()
    }
    pub fn substance_groups(&self) -> Vec<ck::SubstanceGroup> {
        self.inner.substance_groups().to_vec()
    }
    pub fn data_fields(&self) -> Vec<(ck::PropertyText, ck::PropertyText)> {
        self.inner.data_fields().to_vec()
    }
    pub fn source_coordinate_dim(&self) -> Option<ck::CoordinateDimension> {
        self.inner.source_coordinate_dim()
    }
    pub fn title(&self) -> Option<ck::PropertyText> {
        self.inner.title().cloned()
    }
    pub fn index(&self) -> usize {
        self.inner.index()
    }
    pub fn data_field(&self, name: &str) -> Option<ck::PropertyText> {
        self.inner.data_field(name).cloned()
    }
}
impl Molecule {
    pub fn from_sdf(text: &str) -> Result<Self, ck::SdfError> {
        // COSMolKit❗✔️: ck::Molecule::from_sdf(text)
        ck::Molecule::from_sdf(text).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn from_sdf_with_params(
        text: &str,
        params: &ck::SdfReadParams,
    ) -> Result<Self, ck::SdfError> {
        // COSMolKit❗✔️: ck::Molecule::from_sdf_with_params(text,params)
        ck::Molecule::from_sdf_with_params(text, params).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn from_mol(text: &str) -> Result<Self, ck::SdfError> {
        // COSMolKit❗✔️: ck::Molecule::from_mol(text)
        ck::Molecule::from_mol(text).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn from_mol_with_params(
        text: &str,
        params: &ck::SdfReadParams,
    ) -> Result<Self, ck::SdfError> {
        // COSMolKit❗✔️: ck::Molecule::from_mol_with_params(text,params)
        ck::Molecule::from_mol_with_params(text, params).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn read_mol(text: &str) -> Result<Self, ck::MolecularIoError> {
        // COSMolKit❗✔️: ck::Molecule::read_mol(text)
        ck::Molecule::read_mol(text).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn read_mol_with_params(
        text: &str,
        params: &ck::SdfReadParams,
    ) -> Result<Self, ck::MolecularIoError> {
        // COSMolKit❗✔️: ck::Molecule::read_mol_with_params(text,params)
        ck::Molecule::read_mol_with_params(text, params).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn read_sdf(text: &str) -> Result<Self, ck::MolecularIoError> {
        // COSMolKit❗✔️: ck::Molecule::read_sdf(text)
        ck::Molecule::read_sdf(text).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn read_sdf_with_params(
        text: &str,
        params: &ck::SdfReadParams,
    ) -> Result<Self, ck::MolecularIoError> {
        // COSMolKit❗✔️: ck::Molecule::read_sdf_with_params(text,params)
        ck::Molecule::read_sdf_with_params(text, params).map(|inner| Self {
            inner: inner.into(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const SDF: &str = "binding record\n  COSMolKit         2D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\nM  END\n>  <NOTE>\nfirst\nline\n\n>  <NOTE>\nsecond\n\n>  <atom.iprop.score>\n-2147483648 2147483647\n\n>  <atom.prop.label>\nalpha n/a\n\n>  <atom.dprop.weight>\n1.25 -2.5\n\n>  <atom.bprop.flag>\n1 0\n\n>  <bond.iprop.edge>\n7\n\n$$$$\n";
    const QUERY: &str = "binding record\n  COSMolKit         2D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\nM  END\n>  <NOTE>\nfirst\nline\n\n>  <NOTE>\nsecond\n\n>  <atom.iprop.score>\n-2147483648 2147483647\n\n>  <atom.prop.label>\nalpha n/a\n\n>  <atom.dprop.weight>\n1.25 -2.5\n\n>  <atom.bprop.flag>\n1 0\n\n>  <bond.iprop.edge>\n7\n\n$$$$\n";
    const GROUPS: &str = "binding record\n  COSMolKit         2D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\nM  STY  1   1 SUP\nM  SAL   1  1   1\nM  SBL   1  1   1\nM  SMT   1 Me\nM  SDI   1  4    0.0000    1.0000    2.0000    3.0000\nM  SBV   1   1    0.5000    0.2500\nM  END\n$$$$\n";
    #[test]
    fn sdf_reading_record_values_and_all_concrete_constructors_match_canonical() {
        let record = SdfRecord::from_sdf(SDF).unwrap();
        let core = ck::SdfRecord::from_sdf(SDF).unwrap();
        assert_eq!(record.title(), Some("binding record".into()));
        assert_eq!(record.index(), 0);
        assert_eq!(record.data_field("NOTE"), Some("first\nline".into()));
        assert_eq!(record.data_field("missing"), None);
        assert_eq!(record.data_fields(), core.data_fields());
        assert_eq!(record.properties(), *core.properties());
        assert_eq!(record.properties().sdf_property_lists().len(), 5);
        assert_eq!(
            record
                .properties()
                .computed_prop_names()
                .unwrap()
                .unwrap()
                .iter()
                .map(ck::PropertyText::as_bytes)
                .collect::<Vec<_>>(),
            // MolFileParser.cpp sanitizes (including aromaticity) before
            // assignStereochemistry; both set computed integer properties.
            vec![b"numArom".as_slice(), b"_StereochemDone".as_slice()]
        );
        assert_eq!(
            record.properties().prop("_StereochemDone"),
            Some(&ck::PropertyValue::Int(1))
        );
        assert_eq!(
            record.properties().prop("numArom"),
            Some(&ck::PropertyValue::Int(0))
        );
        assert_eq!(
            record.source_coordinate_dim(),
            Some(ck::CoordinateDimension::TwoD)
        );
        let graph = record.graph();
        assert_eq!(graph.kind(), "molecule");
        assert!(graph.query_graph().is_none());
        assert_eq!(
            *graph.molecule().unwrap().inner.borrow(),
            *core.molecule().unwrap()
        );
        assert_eq!(
            *Molecule::from_sdf(SDF).unwrap().inner.borrow(),
            *core.molecule().unwrap()
        );
        let mol = SDF.split('>').next().unwrap();
        assert_eq!(
            *Molecule::from_mol(mol).unwrap().inner.borrow(),
            ck::Molecule::from_mol(mol).unwrap()
        );
        for sanitize in [false, true] {
            for remove_hydrogens in [false, true] {
                for strict_parsing in [false, true] {
                    for expand_attachment_points in [false, true] {
                        for process_property_lists in [false, true] {
                            for coordinate_mode in [
                                ck::SdfCoordinateMode::Preserve,
                                ck::SdfCoordinateMode::Require2D,
                                ck::SdfCoordinateMode::Require3D,
                            ] {
                                let p = ck::SdfReadParams {
                                    sanitize,
                                    remove_hydrogens,
                                    strict_parsing,
                                    expand_attachment_points,
                                    process_property_lists,
                                    coordinate_mode,
                                };
                                let got = SdfRecord::from_sdf_with_params(SDF, &p).unwrap();
                                let expected =
                                    ck::SdfRecord::from_sdf_with_params(SDF, &p).unwrap();
                                assert_eq!(
                                    *got.molecule().unwrap().inner.borrow(),
                                    *expected.molecule().unwrap()
                                );
                                assert_eq!(got.properties(), *expected.properties());
                                assert_eq!(
                                    *Molecule::from_sdf_with_params(SDF, &p)
                                        .unwrap()
                                        .inner
                                        .borrow(),
                                    ck::Molecule::from_sdf_with_params(SDF, &p).unwrap()
                                );
                                assert_eq!(
                                    *Molecule::from_mol_with_params(mol, &p)
                                        .unwrap()
                                        .inner
                                        .borrow(),
                                    ck::Molecule::from_mol_with_params(mol, &p).unwrap()
                                );
                            }
                        }
                    }
                }
            }
        }
        let q = SdfRecord::from_sdf(QUERY).unwrap();
        assert_eq!(q.graph().kind(), "query_graph");
        assert!(q.graph().molecule().is_none());
        assert_eq!(q.graph().query_graph().unwrap().num_atoms(), 2);
        assert_eq!(q.query_graph().unwrap().num_bonds(), 1);
        assert!(matches!(
            q.molecule(),
            Err(ck::SdfError::WrongGraphKind {
                expected: "molecule",
                actual: "query_graph"
            })
        ));
        assert!(matches!(
            Molecule::from_sdf(QUERY),
            Err(ck::SdfError::QueryRecord)
        ));
        assert!(matches!(
            record.query_graph(),
            Err(ck::SdfError::WrongGraphKind {
                expected: "query_graph",
                actual: "molecule"
            })
        ));
    }
    #[test]
    fn sdf_reading_native_files_groups_and_exact_property_kinds() {
        let path = std::env::temp_dir().join(format!(
            "cosmolkit-binding-sdf-{}-{}.sdf",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&path, SDF).unwrap();
        let path_text = path.to_str().unwrap();
        let params = ck::SdfReadParams::default();
        assert_eq!(
            *Molecule::read_sdf(path_text).unwrap().inner.borrow(),
            ck::Molecule::read_sdf(path_text).unwrap()
        );
        assert_eq!(
            *Molecule::read_sdf_with_params(path_text, &params)
                .unwrap()
                .inner
                .borrow(),
            ck::Molecule::read_sdf_with_params(path_text, &params).unwrap()
        );
        assert_eq!(
            *Molecule::read_mol(path_text).unwrap().inner.borrow(),
            ck::Molecule::read_mol(path_text).unwrap()
        );
        assert_eq!(
            *Molecule::read_mol_with_params(path_text, &params)
                .unwrap()
                .inner
                .borrow(),
            ck::Molecule::read_mol_with_params(path_text, &params).unwrap()
        );
        std::fs::remove_file(&path).unwrap();
        assert!(matches!(
            Molecule::read_sdf(path_text),
            Err(ck::MolecularIoError::Io { .. })
        ));
        let g = SdfRecord::from_sdf(GROUPS).unwrap().substance_groups();
        assert_eq!(g.len(), 1);
        assert_eq!(g[0].kind(), &ck::SubstanceGroupKind::Superatom);
        assert_eq!(
            g[0].display().unwrap().brackets()[0].points(),
            &[[0., 1., 0.], [2., 3., 0.], [0., 0., 0.]]
        );
        assert_eq!(g[0].cstates()[0].vector(), &[0.5, 0.25, 0.]);
        let u = ck::PropertyValue::UInt(u32::MAX);
        assert_eq!(u.as_uint().unwrap(), u32::MAX);
        assert_eq!(
            u.as_int().unwrap_err().actual(),
            ck::PropertyValueKind::UInt
        );
        let v = ck::PropertyValue::IntVector(vec![i32::MIN, i32::MAX]);
        assert_eq!(v.as_int_vector().unwrap(), [i32::MIN, i32::MAX]);
    }
}
