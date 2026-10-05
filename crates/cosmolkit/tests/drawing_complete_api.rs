#![cfg(all(feature = "cap-depict", feature = "cap-smiles"))]

use cosmolkit::{
    Conformer2D, Conformer3D, Coordinate2DParams, CoordinateBlock, DrawingWriteError, Molecule,
};

#[test]
fn draw_complete_value_inplace_and_presence_preserve_peer_and_3d() {
    let parsed = Molecule::from_smiles("CCO").unwrap();
    let coordinates = CoordinateBlock {
        conformers_3d: vec![Conformer3D::new(
            29,
            vec![[0.0, -0.0, 0.0], [1.0, 0.5, 0.25], [2.0, 1.0, 1.25]],
            true,
        )],
        ..CoordinateBlock::default()
    };
    let mut source = Molecule::from_parts(
        parsed.topology().clone(),
        coordinates,
        parsed.properties().clone(),
    )
    .unwrap();
    let peer = source.clone();
    let three_d = source.conformers_3d().to_vec();
    let peer_before = peer.to_builder().coordinates().clone();
    assert!(!source.has_2d_coordinates());
    let value = source.with_2d_coordinates().unwrap();
    assert!(value.has_2d_coordinates());
    assert!(!source.has_2d_coordinates());
    source.compute_2d_coordinates_().unwrap();
    assert_eq!(source.coordinates_2d(), value.coordinates_2d());
    assert_eq!(source.conformers_3d(), three_d);
    assert_eq!(peer.to_builder().coordinates(), &peer_before);
    let topology = source.topology().clone();
    let params = Coordinate2DParams {
        clear_existing_2d: false,
        ..Coordinate2DParams::default()
    };
    source.compute_2d_coordinates_with_params_(&params).unwrap();
    assert_eq!(source.to_builder().coordinates().conformers_2d.len(), 2);
    assert_eq!(source.topology(), &topology);
    assert_eq!(source.conformers_3d(), three_d);
    let absent = Molecule::new();
    assert!(!absent.has_2d_coordinates());
    let stored_empty = Molecule::from_parts(
        absent.topology().clone(),
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(17, vec![])],
            ..CoordinateBlock::default()
        },
        absent.properties().clone(),
    )
    .unwrap();
    assert!(stored_empty.has_2d_coordinates());
    assert_eq!(stored_empty.coordinates_2d(), Some([].as_slice()));
}

#[test]
fn draw_complete_files_equal_serialization_and_render_error_preserves_destination() {
    let mol = Molecule::from_smiles("CCO").unwrap();
    let dir = std::env::temp_dir().join(format!(
        "cosmolkit-p5-drawing-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir(&dir).unwrap();
    let svg = dir.join("drawing.svg");
    let png = dir.join("drawing.png");
    mol.write_svg(&svg, 120, 80).unwrap();
    mol.write_png(&png, 120, 80).unwrap();
    assert_eq!(
        std::fs::read(&svg).unwrap(),
        mol.to_svg(120, 80).unwrap().as_bytes()
    );
    assert_eq!(std::fs::read(&png).unwrap(), mol.to_png(120, 80).unwrap());
    for (path, png_output) in [(&svg, false), (&png, true)] {
        let before = std::fs::read(path).unwrap();
        let error = if png_output {
            mol.write_png(path, 0, 80)
        } else {
            mol.write_svg(path, 0, 80)
        }
        .unwrap_err();
        assert!(matches!(
            error,
            DrawingWriteError::Drawing(cosmolkit::DrawingError::InvalidDimensions {
                width: 0,
                height: 80
            })
        ));
        assert_eq!(std::fs::read(path).unwrap(), before);
    }
    let missing = dir.join("absent/drawing.svg");
    let error = mol.write_svg(&missing, 120, 80).unwrap_err();
    match error {
        DrawingWriteError::Io { path, source } => {
            assert_eq!(path, missing);
            assert_eq!(source.kind(), std::io::ErrorKind::NotFound);
        }
        other => panic!("expected filesystem failure, got {other}"),
    }
    assert!(!mol.has_2d_coordinates());
    std::fs::remove_file(svg).unwrap();
    std::fs::remove_file(png).unwrap();
    std::fs::remove_dir(dir).unwrap();
}
