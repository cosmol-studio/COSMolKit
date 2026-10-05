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
        let (actual_chars, actual_modes) = parse_draw_chars(input);
        assert_eq!(actual_chars.iter().collect::<String>(), glyphs, "{input:?}");
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
        let (glyphs, modes) = parse_draw_chars(input);
        assert_eq!(glyphs.iter().collect::<String>(), input, "{input:?}");
        assert!(modes.iter().all(|mode| *mode == TextDrawType::Normal));
        let rects = get_string_rects_unsplit(input, 20.0);
        assert_eq!(rects.iter().map(|rect| rect.ch).collect::<String>(), input);
        assert!(
            rects
                .iter()
                .all(|rect| rect.width.is_finite() && rect.height.is_finite())
        );
    }
    let (glyphs, modes) = parse_draw_chars("<sup>&amp;</sup>&lt;");
    assert_eq!(glyphs.iter().collect::<String>(), "&amp;&lt;");
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
            assert_eq!(atom_label_to_pieces(input, orient), vec![pieces]);
            assert_eq!(
                get_string_rects(input, orient, 20.0)
                    .iter()
                    .map(|rect| rect.ch)
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
        let mut svg = String::new();
        draw_atom_label_svg(&mut svg, &label, 20.0);
        assert_eq!(
            text_contents(&svg),
            (expected.to_string(), glyph_count),
            "label {input:?}"
        );
        assert!(svg.lines().all(|line| line.contains("class='atom-7'")));
        let annotation = DrawAnnotation::new(
            input.to_string(),
            TextAlignType::Middle,
            "note".to_string(),
            0.8,
            DVec2::new(30.0, 40.0),
            DrawColour::new(0.0, 0.0, 1.0),
            20.0,
        );
        let mut svg = String::new();
        draw_annotation_svg(&mut svg, &annotation, 20.0);
        assert_eq!(
            text_contents(&svg),
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
        ("<lit>&lt;</lit>discard", "&amp;lt;", 4),
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
                        text_contents(&svg),
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
            let error = render_prepared_svg(&input, 300, 300).unwrap_err();
            assert!(error.source().is_some());
            let DrawingError::Property(source) = error else {
                panic!("wrong typed error: {error}")
            };
            assert_eq!(source.expected(), PropertyValueKind::String);
            assert_eq!(source.actual(), kind);
            assert_eq!((topology, layout, properties, valence, rings), baseline);
        }
    }
    let error = crate::raster::svg_to_png("<svg><text>&unknown;</text></svg>").unwrap_err();
    assert!(matches!(error, DrawingError::SvgParse(_)));
    assert!(error.source().is_some());
}
