//! Thin canonical depiction calls; runtime owns all coordinate installation.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn with_2d_coordinates(&self) -> Result<Self, ck::OperationError> {
        // COSMolKit❗✔️: .with_2d_coordinates()
        self.inner.borrow().with_2d_coordinates().map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn with_2d_coordinates_with_params(
        &self,
        params: &ck::Coordinate2DParams,
    ) -> Result<Self, ck::OperationError> {
        // COSMolKit❗✔️: .with_2d_coordinates_with_params(&params.inner)
        self.inner
            .borrow()
            .with_2d_coordinates_with_params(params)
            .map(|inner| Self {
                inner: inner.into(),
            })
    }
    pub fn to_svg(&self, width: u32, height: u32) -> Result<String, ck::DrawingError> {
        // COSMolKit❗✔️: .to_svg(width, height)
        self.inner.borrow().to_svg(width, height)
    }
    pub fn to_png(&self, width: u32, height: u32) -> Result<Vec<u8>, ck::DrawingError> {
        // COSMolKit❗✔️: .to_png(width, height)
        self.inner.borrow().to_png(width, height)
    }
    pub fn compute_2d_coordinates_(&self) -> Result<(), ck::OperationError> {
        // COSMolKit❗✔️: .compute_2d_coordinates_()
        self.inner.borrow_mut().compute_2d_coordinates_()
    }
    pub fn compute_2d_coordinates_with_params_(
        &self,
        params: &ck::Coordinate2DParams,
    ) -> Result<(), ck::OperationError> {
        // COSMolKit❗✔️: .compute_2d_coordinates_with_params_(&params.inner)
        self.inner
            .borrow_mut()
            .compute_2d_coordinates_with_params_(params)
    }
    pub fn write_svg(
        &self,
        path: &str,
        width: u32,
        height: u32,
    ) -> Result<(), ck::DrawingWriteError> {
        // COSMolKit❗✔️: .write_svg(std::path::Path::new(path), width, height)
        self.inner
            .borrow()
            .write_svg(std::path::Path::new(path), width, height)
    }
    pub fn write_png(
        &self,
        path: &str,
        width: u32,
        height: u32,
    ) -> Result<(), ck::DrawingWriteError> {
        // COSMolKit❗✔️: .write_png(std::path::Path::new(path), width, height)
        self.inner
            .borrow()
            .write_png(std::path::Path::new(path), width, height)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn depict_eight_calls_and_detached_parameters() {
        let m = Molecule::from_smiles("CCO").unwrap();
        let initial = m.to_smiles().unwrap();
        let p = ck::Coordinate2DParams {
            coordinate_map: [(0, [4., 5.])].into(),
            force_rdkit: true,
            ..Default::default()
        };
        let out = m.with_2d_coordinates().unwrap();
        assert!(m.coordinates_2d().is_empty());
        assert_eq!(out.coordinates_2d().len(), 6);
        let pinned = m.with_2d_coordinates_with_params(&p).unwrap();
        assert_eq!(&pinned.coordinates_2d()[..2], &[4., 5.]);
        assert!(m.coordinates_2d().is_empty());
        m.compute_2d_coordinates_().unwrap();
        assert_eq!(m.coordinates_2d().len(), 6);
        m.compute_2d_coordinates_with_params_(&p).unwrap();
        assert_eq!(&m.coordinates_2d()[..2], &[4., 5.]);
        let before = m.coordinates_2d();
        let p_before = p.clone();
        let svg = m.to_svg(120, 100).unwrap();
        let png = m.to_png(120, 100).unwrap();
        assert!(svg.contains("</svg>"));
        assert_eq!(&png[..8], &[137, 80, 78, 71, 13, 10, 26, 10]);
        let root =
            std::env::temp_dir().join(format!("cosmolkit-depict-binding-{}", std::process::id()));
        std::fs::create_dir_all(&root).unwrap();
        let s = root.join("m.svg");
        let pth = root.join("m.png");
        m.write_svg(s.to_str().unwrap(), 120, 100).unwrap();
        m.write_png(pth.to_str().unwrap(), 120, 100).unwrap();
        assert_eq!(std::fs::read(&s).unwrap(), svg.as_bytes());
        assert_eq!(std::fs::read(&pth).unwrap(), png);
        for path in [&s, &pth] {
            std::fs::write(path, b"original destination").unwrap();
        }
        assert!(matches!(
            m.write_svg(s.to_str().unwrap(), 0, 100),
            Err(ck::DrawingWriteError::Drawing(
                ck::DrawingError::InvalidDimensions {
                    width: 0,
                    height: 100
                }
            ))
        ));
        assert!(matches!(
            m.write_png(pth.to_str().unwrap(), 120, 0),
            Err(ck::DrawingWriteError::Drawing(
                ck::DrawingError::InvalidDimensions {
                    width: 120,
                    height: 0
                }
            ))
        ));
        for path in [&s, &pth] {
            assert_eq!(std::fs::read(path).unwrap(), b"original destination");
        }
        let bad = ck::Coordinate2DParams {
            coordinate_map: [(99, [1., 2.])].into(),
            ..Default::default()
        };
        assert!(matches!(
            m.compute_2d_coordinates_with_params_(&bad),
            Err(ck::OperationError::Coordinate2D(
                ck::Coordinate2DError::Fragment(ck::Coordinate2DLayoutError::AtomIndexOutOfRange {
                    atom: 99,
                    atom_count: 3
                })
            ))
        ));
        assert_eq!(m.coordinates_2d(), before);
        assert_eq!(m.to_smiles().unwrap(), initial);
        assert_eq!(p, p_before);
    }
}
