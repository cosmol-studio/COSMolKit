//! Frozen legacy SVG rasterization with the identical embedded font/backends.
use crate::DrawingError;
use crate::draw::{EMBEDDED_DRAW_FONT_FAMILY, embedded_draw_font_data};

pub(crate) fn svg_to_png(svg: &[u8]) -> Result<Vec<u8>, DrawingError> {
    let mut opt = usvg::Options::default();
    opt.font_family = EMBEDDED_DRAW_FONT_FAMILY.to_owned();
    let fontdb = opt.fontdb_mut();
    fontdb.load_font_source(usvg::fontdb::Source::Binary(embedded_draw_font_data()));
    fontdb.set_sans_serif_family(EMBEDDED_DRAW_FONT_FAMILY);
    let tree = usvg::Tree::from_data(svg, &opt)?;
    let size = tree.size().to_int_size();
    let mut pixmap = tiny_skia::Pixmap::new(size.width(), size.height()).ok_or(
        DrawingError::PixmapAllocation {
            width: size.width(),
            height: size.height(),
        },
    )?;
    resvg::render(&tree, tiny_skia::Transform::default(), &mut pixmap.as_mut());
    pixmap
        .encode_png()
        .map_err(|error| DrawingError::PngEncode(Box::new(error)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn references(label: &str) -> PathBuf {
        let folder = if matches!(label, "benzene" | "pyridine") {
            "source_aromatic_351f8f3"
        } else {
            "legacy_0_3_0"
        };
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../testdata/depiction/fixtures")
            .join(folder)
    }

    fn check_frozen(label: &str, width: u32, height: u32) {
        let base = references(label);
        let svg = std::fs::read_to_string(base.join(format!("{label}.svg"))).unwrap();
        let expected_png = std::fs::read(base.join(format!("{label}.png"))).unwrap();
        let expected_rgba = std::fs::read(base.join(format!("{label}.rgba"))).unwrap();
        let actual = svg_to_png(svg.as_bytes()).unwrap();
        assert_eq!(&actual[..8], b"\x89PNG\r\n\x1a\n", "{label}");
        assert_eq!(actual, expected_png, "frozen encoded PNG: {label}");
        let decoded = tiny_skia::Pixmap::decode_png(&actual).unwrap();
        let reference = tiny_skia::Pixmap::decode_png(&expected_png).unwrap();
        assert_eq!(
            (decoded.width(), decoded.height()),
            (width, height),
            "{label}"
        );
        assert_eq!(
            decoded.data(),
            reference.data(),
            "decoded legacy PNG: {label}"
        );
        assert_eq!(decoded.data(), expected_rgba, "frozen RGBA: {label}");
    }

    #[test]
    fn drawing_raster_twelve_frozen_png_rgba_cases() {
        let labels = [
            "empty",
            "single_carbon",
            "explicit_h_methane",
            "ethane",
            "ethanol",
            "benzene",
            "pyridine",
            "carbonyl",
            "isotope",
            "formal_charge",
            "tetrahedral_wedge",
            "mapped_atoms",
        ];
        let mut calls = 0;
        for label in labels {
            check_frozen(label, 300, 300);
            calls += 1;
        }
        assert_eq!(calls, 12);
    }

    #[test]
    fn drawing_raster_svg_to_png_rasterizes_embedded_font_text() {
        // Original legacy font regression: same literal and every assertion.
        let svg = "<?xml version='1.0' encoding='iso-8859-1'?>\
                   <svg version='1.1' baseProfile='full' xmlns='http://www.w3.org/2000/svg' \
                   width='120px' height='80px' viewBox='0 0 120 80'>\
                   <rect width='120' height='80' fill='#FFFFFF'/>\
                   <text x='12' y='54' style='font-size:48px;font-style:normal;\
                   font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;\
                   text-anchor:start;fill:#000000'>O</text>\
                   </svg>";
        let png = svg_to_png(svg.as_bytes()).expect("text png rasterization");
        let pixmap = tiny_skia::Pixmap::decode_png(&png).expect("decode rendered png");
        let text_pixels = pixmap
            .pixels()
            .iter()
            .filter(|pixel| {
                pixel.alpha() == 255
                    && pixel.red() < 245
                    && pixel.green() < 245
                    && pixel.blue() < 245
            })
            .count();
        assert!(
            text_pixels > 0,
            "embedded-font text should produce non-background pixels"
        );
        assert_eq!(
            svg,
            std::fs::read_to_string(references("font_text").join("font_text.svg")).unwrap()
        );
        check_frozen("font_text", 120, 80);
    }

    #[test]
    fn drawing_raster_malformed_svg_and_invalid_size_typed_causes() {
        use std::error::Error;
        let malformed = svg_to_png(b"<svg").unwrap_err();
        assert!(matches!(
            malformed,
            DrawingError::SvgParse(usvg::Error::ParsingFailed(_))
        ));
        assert!(malformed.source().is_some());
        let invalid =
            svg_to_png(b"<svg xmlns='http://www.w3.org/2000/svg' width='0' height='80'/>")
                .unwrap_err();
        assert!(matches!(
            invalid,
            DrawingError::SvgParse(usvg::Error::InvalidSize)
        ));
        assert!(invalid.source().is_some());
    }

    #[test]
    fn drawing_raster_zero_public_canvas_dimensions_are_typed() {
        let topology = cosmolkit_model::TopologyBlock::default();
        let coordinates = cosmolkit_model::CoordinateBlock::default();
        let properties = cosmolkit_model::MoleculeProperties::default();
        for (width, height) in [(0, 80), (120, 0), (0, 0)] {
            let input = || crate::DrawingInput {
                topology: &topology,
                coordinates: &coordinates,
                properties: &properties,
                valence: None,
                rings: None,
            };
            let options = crate::DepictOptions { width, height };
            for error in [
                crate::render_svg(input(), &options).unwrap_err(),
                crate::render_png(input(), &options).unwrap_err(),
            ] {
                assert!(
                    matches!(error, DrawingError::InvalidDimensions { width:w, height:h } if w==width && h==height)
                );
            }
        }
    }
}
