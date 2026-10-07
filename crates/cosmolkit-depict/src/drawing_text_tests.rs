//! Fixed ASCII expectations from RDKit DrawText.cpp and DrawTextSVG.cpp at
//! revision 351f8f378f8ad6bbd517980c38896e66bf907af8c. No oracle or generator.
use super::*;
use cosmolkit_core::{RingFindType, RingInfo, ValenceAssignment};
use cosmolkit_model::{
    AtomId, AtomSpec, BondId, BondSpec, Conformer2D, Element, MoleculeProperties, PropertyValue,
    PropertyValueKind,
};
use std::error::Error;

#[test]
fn drawing_text_script_mode_transitions() {
    // Upstream has one mode, not a stack; either closing tag resets it.
    let cases = [
        ("", "", ""),
        ("C", "C", "N"),
        ("<sub>2</sub>C", "2C", "DN"),
        ("<sup>3</sup>C", "3C", "UN"),
        ("<sub>2", "2", "D"),
        ("<sup>3", "3", "U"),
        ("</sub>C</sup>", "C", "N"),
        ("<sub>1<sup>2</sup>3</sub>4", "1234", "DUNN"),
        ("<sup>1<sub>2</sup>3</sub>4", "1234", "UDNN"),
        ("<sup>1</sub>2", "12", "UN"),
        ("<sub>1</sup>2", "12", "DN"),
        ("<sup><sub></sup></sub>", "", ""),
        ("<sup>1<sup>2</sup>3", "123", "UUN"),
        ("<sub>1<sub>2</sub>3", "123", "DDN"),
        ("<sup><supX>1</sup>", "<supX>1", "UUUUUUU"),
        ("C<sub", "C<sub", "NNNNN"),
        ("<SUP>1</SUP>", "<SUP>1</SUP>", "NNNNNNNNNNNN"),
        ("<sub >2</sub >", "<sub >2</sub >", "NNNNNNNNNNNNNN"),
    ];
    for (input, glyphs, modes) in cases {
        let (actual_chars, actual_modes) = parse_draw_chars(input.as_bytes());
        assert_eq!(
            std::str::from_utf8(&actual_chars).expect("fixed ASCII glyph bytes"),
            glyphs,
            "{input:?}"
        );
        let actual_modes: String = actual_modes
            .iter()
            .map(|mode| match mode {
                TextDrawType::Normal => 'N',
                TextDrawType::Subscript => 'D',
                TextDrawType::Superscript => 'U',
            })
            .collect();
        assert_eq!(actual_modes, modes, "{input:?}");
    }
}

#[test]
fn drawing_text_entities_and_unknown_tags_are_literal() {
    for input in [
        "&amp;",
        "&lt;",
        "&gt;",
        "&quot;",
        "&apos;",
        "&#65;",
        "&#x41;",
        "&amp;lt;",
        "&unknown;",
        "&",
        "<",
        "<sup",
        "</sub",
        "<foo>C</foo>",
        "<lit>C</lit>",
        "A<lit>B</lit>",
        "</lit>",
        "<subscript>",
    ] {
        let (glyphs, modes) = parse_draw_chars(input.as_bytes());
        assert_eq!(
            std::str::from_utf8(&glyphs).expect("fixed ASCII glyph bytes"),
            input,
            "{input:?}"
        );
        assert!(modes.iter().all(|mode| *mode == TextDrawType::Normal));
        let rects = get_string_rects_unsplit(input.as_bytes(), 20.0);
        assert_eq!(
            rects
                .iter()
                .map(|rect| char::from(rect.ch))
                .collect::<String>(),
            input
        );
        assert!(
            rects
                .iter()
                .all(|rect| rect.width.is_finite() && rect.height.is_finite())
        );
    }
    let (glyphs, modes) = parse_draw_chars(b"<sup>&amp;</sup>&lt;");
    assert_eq!(
        std::str::from_utf8(&glyphs).expect("fixed ASCII glyph bytes"),
        "&amp;&lt;"
    );
    assert_eq!(
        modes,
        [
            TextDrawType::Superscript,
            TextDrawType::Superscript,
            TextDrawType::Superscript,
            TextDrawType::Superscript,
            TextDrawType::Superscript,
            TextDrawType::Normal,
            TextDrawType::Normal,
            TextDrawType::Normal,
            TextDrawType::Normal,
        ]
    );
}

#[test]
fn drawing_text_leading_literal_wrapper_belongs_to_label_splitter() {
    for orient in [
        OrientType::C,
        OrientType::E,
        OrientType::W,
        OrientType::N,
        OrientType::S,
    ] {
        for (input, pieces, glyphs) in [
            ("<lit>CH3</lit>discard", "CH3", "CH3"),
            ("<lit>CH3", "CH3", "CH3"),
            ("<lit>&amp;</lit>", "&amp;", "&amp;"),
            ("<lit>H<sub>2</sub></lit>", "H<sub>2</sub>", "H2"),
        ] {
            assert_eq!(
                atom_label_to_pieces(input.as_bytes(), orient),
                vec![pieces.as_bytes().to_vec()]
            );
            assert_eq!(
                get_string_rects(input.as_bytes(), orient, 20.0)
                    .iter()
                    .map(|rect| char::from(rect.ch))
                    .collect::<String>(),
                glyphs
            );
        }
    }
}

// Read every text element's serialized content; no escaping or markup parser
// from production is used to build expectations.
fn text_contents(svg: &str) -> (String, usize) {
    let mut contents = String::new();
    let mut count = 0;
    for line in svg.lines().filter(|line| line.starts_with("<text ")) {
        let start = line.find('>').unwrap() + 1;
        let end = line.rfind("</text>").unwrap();
        contents.push_str(&line[start..end]);
        count += 1;
    }
    (contents, count)
}

#[test]
fn drawing_text_svg_glyphs_escape_once_at_output() {
    for (input, expected, glyph_count) in [
        ("&<>\"'", "&amp;&lt;&gt;&quot;&apos;", 5),
        (
            "&amp;&lt;&gt;&quot;&apos;",
            "&amp;amp;&amp;lt;&amp;gt;&amp;quot;&amp;apos;",
            25,
        ),
        ("A<lit>B</lit>", "A&lt;lit&gt;B&lt;/lit&gt;", 13),
        ("<sup>&amp;</sup>", "&amp;amp;", 5),
        ("C<sub>2</sub><sup>+</sup>", "C2+", 3),
    ] {
        let label = AtomLabel::new(
            input,
            7,
            6,
            OrientType::E,
            DVec2::new(30.0, 40.0),
            DrawColour::new(0.0, 0.0, 0.0),
            20.0,
        );
        let mut svg = Vec::new();
        draw_atom_label_svg(&mut svg, &label, 20.0);
        assert_eq!(
            text_contents(std::str::from_utf8(&svg).expect("fixed source SVG UTF8")),
            (expected.to_string(), glyph_count),
            "label {input:?}"
        );
        assert!(
            std::str::from_utf8(&svg)
                .expect("fixed source SVG UTF8")
                .lines()
                .all(|line| line.contains("class='atom-7'"))
        );
        let annotation = DrawAnnotation::new(
            input.to_string(),
            TextAlignType::Middle,
            "note".to_string(),
            0.8,
            DVec2::new(30.0, 40.0),
            DrawColour::new(0.0, 0.0, 1.0),
            20.0,
            1.0,
        );
        let mut svg = Vec::new();
        draw_annotation_svg(&mut svg, &annotation, 20.0);
        assert_eq!(
            text_contents(std::str::from_utf8(&svg).expect("fixed source SVG UTF8")),
            (expected.to_string(), glyph_count),
            "note {input:?}"
        );
        if input == "C<sub>2</sub><sup>+</sup>" {
            assert_eq!(label.rects[1].draw_mode, TextDrawType::Subscript);
            assert_eq!(label.rects[2].draw_mode, TextDrawType::Superscript);
        }
    }
}

fn fixed_scene(location: usize, value: PropertyValue) -> (TopologyBlock, MoleculeProperties) {
    let mut atom = AtomSpec::new(Element::C);
    let mut bond = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single);
    let mut properties = MoleculeProperties::default();
    match location {
        0 => atom = atom.with_prop("atomNote", value).unwrap(),
        1 => bond = bond.with_prop("bondNote", value).unwrap(),
        2 => {
            properties = properties
                .with_prop("molNote", value.as_string().unwrap())
                .unwrap()
        }
        _ => panic!("invalid test location"),
    }
    (
        TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(AtomId::new(0), atom),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            ],
            vec![Bond::from_spec(BondId::new(0), bond)],
            vec![],
            vec![],
        )
        .unwrap(),
        properties,
    )
}

#[test]
fn drawing_text_fixed_prepared_note_product_preserves_input() {
    let mut calls = 0;
    for (note, expected, glyph_count) in [
        ("", "", 0),
        ("&amp;", "&amp;amp;", 5),
        ("&lt;", "&amp;lt;", 4),
        ("&gt;", "&amp;gt;", 4),
        ("&quot;", "&amp;quot;", 6),
        ("&apos;", "&amp;apos;", 6),
        ("&#65;", "&amp;#65;", 5),
        ("<sup>&amp;</sup>", "&amp;amp;", 5),
        ("C<sub>2</sub>", "C2", 2),
        ("<sup>1<sub>2</sup>3", "123", 3),
        ("A<lit>B</lit>", "A&lt;lit&gt;B&lt;/lit&gt;", 13),
        (
            "<lit>&lt;</lit>discard",
            "&lt;lit&gt;&amp;lt;&lt;/lit&gt;discard",
            22,
        ),
        ("<sup", "&lt;sup", 4),
        ("<foo>C</foo>", "&lt;foo&gt;C&lt;/foo&gt;", 12),
    ] {
        for location in 0..3 {
            let (topology, properties) = fixed_scene(location, PropertyValue::String(note.into()));
            let layout = Conformer2D::new(17, vec![[-0.0, 0.0], [1.5, 0.0]]);
            let valence = ValenceAssignment {
                explicit_valence: vec![1, 1],
                implicit_hydrogens: vec![3, 3],
            };
            let rings = RingInfo::new(RingFindType::SymmSssr, 2, 1);
            let baseline = (
                topology.clone(),
                layout.clone(),
                properties.clone(),
                valence.clone(),
                rings.clone(),
            );
            let coordinate_bits: Vec<_> = layout
                .coordinates()
                .iter()
                .flatten()
                .map(|x| x.to_bits())
                .collect();
            let input = PreparedDrawingInput {
                topology: &topology,
                layout: &layout,
                properties: &properties,
                valence: &valence,
                rings: &rings,
            };
            for (width, height) in [(120, 80), (300, 300)] {
                let mut previous = None;
                for _ in 0..2 {
                    let svg = render_prepared_svg(&input, width, height).unwrap();
                    calls += 1;
                    assert_eq!(
                        text_contents(std::str::from_utf8(&svg).expect("fixed source SVG UTF8")),
                        (expected.to_string(), glyph_count),
                        "note {note:?} location={location} canvas={width}x{height}"
                    );
                    let png = crate::raster::svg_to_png(&svg).unwrap();
                    let pixmap = tiny_skia::Pixmap::decode_png(&png).unwrap();
                    assert_eq!((pixmap.width(), pixmap.height()), (width, height));
                    if let Some((old_svg, old_png)) = &previous {
                        assert_eq!(&svg, old_svg);
                        assert_eq!(&png, old_png);
                    }
                    previous = Some((svg, png));
                    assert_eq!(
                        (&topology, &layout, &properties, &valence, &rings),
                        (
                            &baseline.0,
                            &baseline.1,
                            &baseline.2,
                            &baseline.3,
                            &baseline.4
                        )
                    );
                    assert_eq!(
                        layout
                            .coordinates()
                            .iter()
                            .flatten()
                            .map(|x| x.to_bits())
                            .collect::<Vec<_>>(),
                        coordinate_bits
                    );
                }
            }
        }
    }
    assert_eq!(calls, 168);
}

#[test]
fn drawing_text_property_and_svg_errors_keep_typed_causes() {
    for location in 0..2 {
        for value in [
            PropertyValue::Int(7),
            PropertyValue::Double(1.5),
            PropertyValue::Bool(true),
        ] {
            let kind = value.kind();
            let (topology, properties) = fixed_scene(location, value);
            let layout = Conformer2D::new(17, vec![[-0.0, 0.0], [1.5, 0.0]]);
            let valence = ValenceAssignment {
                explicit_valence: vec![1, 1],
                implicit_hydrogens: vec![3, 3],
            };
            let rings = RingInfo::new(RingFindType::SymmSssr, 2, 1);
            let baseline = (
                topology.clone(),
                layout.clone(),
                properties.clone(),
                valence.clone(),
                rings.clone(),
            );
            let input = PreparedDrawingInput {
                topology: &topology,
                layout: &layout,
                properties: &properties,
                valence: &valence,
                rings: &rings,
            };
            // DrawMol::extractAtomNotes/extractBondNotes read std::string;
            // Dict::getVal documents source lexical conversion of scalars.
            let svg = render_prepared_svg(&input, 300, 300).unwrap();
            let expected = match kind {
                PropertyValueKind::Int => "7",
                PropertyValueKind::Double => "1.5",
                PropertyValueKind::Bool => "1",
                _ => unreachable!("the original scalar fixture contains these three kinds"),
            };
            assert_eq!(
                text_contents(std::str::from_utf8(&svg).expect("fixed SVG UTF8")),
                (expected.to_owned(), expected.len())
            );
            assert_eq!((topology, layout, properties, valence, rings), baseline);
        }
    }
    let error = crate::raster::svg_to_png(b"<svg><text>&unknown;</text></svg>").unwrap_err();
    assert!(matches!(error, DrawingError::SvgParse(_)));
    assert!(error.source().is_some());
}

// DRAW03 proposal only; pinned source 351f8f378f8ad6bbd517980c38896e66bf907af8.
// Positive-width printable ASCII with source-defined script reference state.
#[test]
fn drawing_text_annotation_unsplit_declared_alignment_geometry() {
    // DrawAnnotation.cpp:71-81; DrawTextNotFT.cpp:27-85;
    // DrawTextSVG.cpp:74-142 and StringRect.h:62-79.
    // AB at effective size 16: both widths .6*16=9.6, height .8*16=12.8,
    // running x step 9.6*1.15=11.04. Normal bounds span 20.64.
    // Proposed tolerance is new-case-only and awaits p1/ROOT approval.
    for (align, expected_x) in [
        (TextAlignType::Start, [0.0, 11.04]),
        (TextAlignType::Middle, [-10.32, 0.72]),
        (TextAlignType::End, [-11.04, 0.0]),
    ] {
        let annotation = DrawAnnotation::new(
            "AB",
            align,
            "note".into(),
            0.8,
            DVec2::new(30.0, 40.0),
            DrawColour::new(0.0, 0.0, 1.0),
            20.0,
            1.0,
        );
        assert_eq!(annotation.rects.len(), 2);
        assert_eq!(
            annotation
                .rects
                .iter()
                .map(|r| char::from(r.ch))
                .collect::<String>(),
            "AB"
        );
        for (rect, x) in annotation.rects.iter().zip(expected_x) {
            for (actual, expected) in [
                (rect.trans.x, x),
                (rect.trans.y, 0.0),
                (rect.offset.x, 4.8),
                (rect.offset.y, 8.0),
                (rect.width, 9.6),
                (rect.height, 12.8),
                (rect.y_shift, 0.0),
                (rect.rect_corr, 0.0),
            ] {
                assert!(
                    (actual - expected).abs() <= 1.0e-12,
                    "alignment={align:?}: {actual:?} != {expected:?}"
                );
            }
        }
    }
}

#[test]
fn drawing_text_annotation_source_svg_all_alignments() {
    // DrawText.cpp:555-570 computes each baseline; DrawTextSVG.cpp:49-70
    // serializes source class/style/font truncation and formatDouble %.1f.
    for (align, expected) in [
        (
            TextAlignType::Start,
            concat!(
                "<text x='25.2' y='48.0' class='note' style='font-size:16px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#0000FF' >A</text>\n",
                "<text x='36.2' y='48.0' class='note' style='font-size:16px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#0000FF' >B</text>\n"
            ),
        ),
        (
            TextAlignType::Middle,
            concat!(
                "<text x='14.9' y='48.0' class='note' style='font-size:16px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#0000FF' >A</text>\n",
                "<text x='25.9' y='48.0' class='note' style='font-size:16px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#0000FF' >B</text>\n"
            ),
        ),
        (
            TextAlignType::End,
            concat!(
                "<text x='14.2' y='48.0' class='note' style='font-size:16px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#0000FF' >A</text>\n",
                "<text x='25.2' y='48.0' class='note' style='font-size:16px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#0000FF' >B</text>\n"
            ),
        ),
    ] {
        let annotation = DrawAnnotation::new(
            "AB",
            align,
            "note".into(),
            0.8,
            DVec2::new(30.0, 40.0),
            DrawColour::new(0.0, 0.0, 1.0),
            20.0,
            1.0,
        );
        let mut svg = Vec::new();
        draw_annotation_svg(&mut svg, &annotation, 20.0);
        assert_eq!(svg.as_slice(), expected.as_bytes());
    }
}

#[test]
fn drawing_text_annotation_source_draw_uses_final_alignment_and_script_baselines() {
    // DrawAnnotation::draw uses the FINAL align_ in drawString, not a stale
    // cached extraction alignment. Source extractBrackets mutates align_ after new.
    let mut shifted = DrawAnnotation::new(
        "AB",
        TextAlignType::End,
        "note".into(),
        0.8,
        DVec2::new(30.0, 40.0),
        DrawColour::new(0.0, 0.0, 1.0),
        20.0,
        1.0,
    );
    shifted.align = TextAlignType::Start;
    let mut shifted_svg = Vec::new();
    draw_annotation_svg(&mut shifted_svg, &shifted, 20.0);
    assert_eq!(
        std::str::from_utf8(&shifted_svg).expect("fixed source SVG UTF8"),
        concat!(
            "<text x='25.2' y='48.0' class='note' style='font-size:16px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#0000FF' >A</text>\n",
            "<text x='36.2' y='48.0' class='note' style='font-size:16px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#0000FF' >B</text>\n",
        )
    );
    // DrawText.cpp:392-436 shifts scripts using C height=12.8; :555-570
    // then emits baselines normal 48.0, sub 54.4, sup 41.6. Scripts share x.
    let mut annotation = DrawAnnotation::new(
        "C<sub>2</sub><sup>+</sup>",
        TextAlignType::End,
        "note".into(),
        0.8,
        DVec2::new(30.0, 40.0),
        DrawColour::new(0.0, 0.0, 1.0),
        20.0,
        1.0,
    );
    annotation.align = TextAlignType::Start;
    let mut svg = Vec::new();
    draw_annotation_svg(&mut svg, &annotation, 20.0);
    assert_eq!(
        std::str::from_utf8(&svg).expect("fixed source SVG UTF8"),
        concat!(
            "<text x='25.2' y='48.0' class='note' style='font-size:16px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#0000FF' >C</text>\n",
            "<text x='36.2' y='54.4' class='note' style='font-size:10px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#0000FF' >2</text>\n",
            "<text x='36.2' y='41.6' class='note' style='font-size:10px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#0000FF' >+</text>\n",
        )
    );
}

#[test]
fn drawing_text_annotation_source_literal_wrappers_are_not_split() {
    // DrawAnnotation.cpp:79-80 dontSplit=true. DrawText.cpp:501-505 bypasses
    // atomLabelToPieces, preserves wrappers. XML escaping only at SVG emission.
    for align in [
        TextAlignType::Start,
        TextAlignType::Middle,
        TextAlignType::End,
    ] {
        for (input, serialized, count) in [
            ("AB", "AB", 2),
            (
                "<lit>&lt;</lit>discard",
                "&lt;lit&gt;&amp;lt;&lt;/lit&gt;discard",
                22,
            ),
            ("<lit>CH3", "&lt;lit&gt;CH3", 8),
            ("A<lit>B</lit>", "A&lt;lit&gt;B&lt;/lit&gt;", 13),
            ("C<sub>2</sub><sup>+</sup>", "C2+", 3),
        ] {
            let annotation = DrawAnnotation::new(
                input,
                align,
                "note".into(),
                0.8,
                DVec2::new(30.0, 40.0),
                DrawColour::new(0.0, 0.0, 1.0),
                20.0,
                1.0,
            );
            assert_eq!(annotation.rects.len(), count);
            let mut svg = Vec::new();
            draw_annotation_svg(&mut svg, &annotation, 20.0);
            assert_eq!(
                text_contents(std::str::from_utf8(&svg).expect("fixed source SVG UTF8")),
                (serialized.into(), count)
            );
        }
    }
}
