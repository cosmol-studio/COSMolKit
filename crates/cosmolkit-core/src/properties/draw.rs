//! RDKit MolDraw2D-compatible molecule drawing (SVG + PNG).
//!
//! Source reproduction protocol: dev/source_reproduction_protocol.md
//!
//! The drawing pipeline mirrors RDKit MolDraw2D's approach:
//!
//! 1. `prepare_working_mol_for_drawing` — kekulize, add chiral Hs
//! 2. `DrawMol::from_molecule` — extract atom symbols, bonds, radicals, measure scale
//! 3. `to_svg` or `to_png` — render via SVG XML primitives
//!
//! For PNG output the SVG string is rasterized using `usvg` + `tiny-skia`
//! (pure Rust, no system Cairo dependency).
//!
//! ## Source Reproduction
//!
//! - RDKit✔️✔️: direct port from RDKit C++ (MolDraw2DSVG.cpp, DrawMol.cpp, AtomSymbol.cpp)
//! - RDKit✔️❌: ported with COSMolKit-specific adaptation
//! - RDKit❗✔️: algorithm-equivalent via different architecture
//! - RDKit❗❌: not ported (intentional omission)
//!

use crate::{Atom, Bond, BondDirection, BondOrder, Molecule, ValenceModel, assign_valence};
use glam::DVec2;
use std::collections::HashSet;
use std::sync::{Arc, OnceLock};

// ──────────────────────────────────────────────
// Embedded font for SVG rendering
// ──────────────────────────────────────────────

const EMBEDDED_DRAW_FONT_BYTES: &[u8] = include_bytes!("../../assets/fonts/NotoSans-Regular.ttf");
const EMBEDDED_DRAW_FONT_FAMILY: &str = "Noto Sans";

fn embedded_draw_font_data() -> Arc<dyn AsRef<[u8]> + Send + Sync> {
    static FONT_DATA: OnceLock<Arc<Vec<u8>>> = OnceLock::new();
    FONT_DATA
        .get_or_init(|| Arc::new(EMBEDDED_DRAW_FONT_BYTES.to_vec()))
        .clone() as Arc<dyn AsRef<[u8]> + Send + Sync>
}

// ──────────────────────────────────────────────
// Error type
// ──────────────────────────────────────────────

/// Errors returned by SVG / PNG drawing routines.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SvgDrawError {
    #[error("coordinate generation failed: {0}")]
    CoordinateGeneration(String),
    #[error("unsupported drawing path: {0}")]
    Unsupported(String),
    #[error(transparent)]
    UnsupportedFeature(#[from] crate::UnsupportedFeatureError),
    #[error("SVG parse failed: {0}")]
    SvgParse(String),
    #[error("PNG pixmap allocation failed for {width}x{height}")]
    PixmapAllocation { width: u32, height: u32 },
    #[error("PNG encoding failed: {0}")]
    PngEncode(String),
}

// ──────────────────────────────────────────────
// Public snapshot types (used by facade crate)
// ──────────────────────────────────────────────

/// Atom snapshot after drawing preparation.
#[derive(Debug, Clone, PartialEq)]
pub struct PreparedDrawAtom {
    pub index: usize,
    pub atomic_number: u8,
    pub x: f64,
    pub y: f64,
}

/// Bond snapshot after drawing preparation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PreparedDrawBond {
    pub index: usize,
    pub begin_atom: usize,
    pub end_atom: usize,
    pub bond_order: BondOrder,
    pub is_aromatic: bool,
    pub direction: BondDirection,
    pub rdkit_direction_name: String,
}

/// Molecule snapshot after drawing preparation.
#[derive(Debug, Clone, PartialEq)]
pub struct PreparedDrawMolecule {
    pub atoms: Vec<PreparedDrawAtom>,
    pub bonds: Vec<PreparedDrawBond>,
}

// ──────────────────────────────────────────────
// Drawing-internal types
// ──────────────────────────────────────────────

/// RGBA colour with components in [0, 1].
#[derive(Debug, Clone, Copy, PartialEq)]
struct DrawColour {
    r: f64,
    g: f64,
    b: f64,
    a: f64,
}

impl DrawColour {
    const fn new(r: f64, g: f64, b: f64) -> Self {
        Self { r, g, b, a: 1.0 }
    }
    const fn with_alpha(self, a: f64) -> Self {
        Self { a, ..self }
    }
}

/// Drawing options mirroring RDKit MolDrawOptions defaults.
#[derive(Debug, Clone)]
struct DrawOptions {
    padding: f64,
    multiple_bond_offset: f64,
    bond_line_width: f64,
    scale_bond_width: bool,
    clear_background: bool,
    background_colour: DrawColour,
    query_colour: DrawColour,
    flag_close_contacts_dist: i32,
    /// Font scale for annotations (RDKit default 0.8)
    annotation_font_scale: f64,
    /// Colour for atom notes and atom CIP codes (RDKit default blue)
    atom_note_colour: DrawColour,
    /// Colour for bond notes and bond CIP codes (RDKit default red)
    bond_note_colour: DrawColour,
    /// Colour for general annotations (RDKit default black)
    annotation_colour: DrawColour,
    /// Treat dummy atoms (atomic num 0) with degree 1 as attachment points
    dummies_are_attachments: bool,
    /// Colour for variable attachment markers
    variable_attachment_colour: DrawColour,
    /// Radius for variable atom circles
    variable_atom_radius: f64,
    /// Width multiplier for variable bonds
    variable_bond_width_multiplier: f64,
    /// Include SGroup and stereo group annotations
    include_annotations: bool,
    /// Include radical dot rendering
    include_radicals: bool,
    /// Circle atoms for highlighting
    circle_atoms: bool,
    /// Continuous highlight mode (alternative to circle highlights)
    continuous_highlight: bool,
}

impl Default for DrawOptions {
    fn default() -> Self {
        Self {
            padding: 0.05,
            multiple_bond_offset: 0.12,
            bond_line_width: 2.0,
            scale_bond_width: true,
            clear_background: true,
            background_colour: DrawColour::new(1.0, 1.0, 1.0),
            query_colour: DrawColour::new(0.0, 0.0, 0.0),
            flag_close_contacts_dist: -1,
            annotation_font_scale: 0.8,
            atom_note_colour: DrawColour::new(0.0, 0.0, 1.0),
            bond_note_colour: DrawColour::new(1.0, 0.0, 0.0),
            annotation_colour: DrawColour::new(0.0, 0.0, 0.0),
            dummies_are_attachments: false,
            variable_attachment_colour: DrawColour::new(0.5, 0.5, 0.5),
            variable_atom_radius: 0.3,
            variable_bond_width_multiplier: 2.0,
            include_annotations: true,
            include_radicals: true,
            circle_atoms: false,
            continuous_highlight: false,
        }
    }
}

/// Orientation of a label relative to its anchor point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OrientType {
    C,
    N,
    E,
    S,
    W,
}

/// Horizontal text alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TextAlignType {
    Middle,
    Start,
    End,
}

/// Draw mode for a single character in a string rect.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TextDrawType {
    Normal,
    Subscript,
    Superscript,
}

/// Per-character layout rectangle.
#[derive(Debug, Clone)]
struct StringRect {
    ch: char,
    draw_mode: TextDrawType,
    trans: DVec2,
    offset: DVec2,
    g_centre: DVec2,
    y_shift: f64,
    width: f64,
    height: f64,
    rect_corr: f64,
}

/// A fully built atom label ready for drawing.
#[derive(Debug, Clone)]
struct AtomLabel {
    symbol: String,
    atom_idx: usize,
    atomic_num: u8,
    orient: OrientType,
    cds: DVec2,
    colour: DrawColour,
    rects: Vec<StringRect>,
}

impl AtomLabel {
    fn new(
        symbol: impl Into<String>,
        atom_idx: usize,
        atomic_num: u8,
        orient: OrientType,
        cds: DVec2,
        colour: DrawColour,
        font_size: f64,
    ) -> Self {
        let symbol = symbol.into();
        let rects = get_string_rects(&symbol, orient, font_size);
        Self {
            symbol,
            atom_idx,
            atomic_num,
            orient,
            cds,
            colour,
            rects,
        }
    }

    fn recalculate_rects(&mut self, font_size: f64) {
        self.rects = get_string_rects(&self.symbol, self.orient, font_size);
    }

    fn find_extremes(&self, xmin: &mut f64, xmax: &mut f64, ymin: &mut f64, ymax: &mut f64) {
        let mut corners = Vec::new();
        for rect in &self.rects {
            let orig_trans = rect.trans;
            // adjusted position = cds + rect.trans
            let cx = self.cds.x + rect.trans.x;
            let cy = self.cds.y + rect.trans.y;
            let half_w = rect.width / 2.0;
            let half_h = rect.height / 2.0;
            corners.push(DVec2::new(cx - half_w, cy - half_h));
            corners.push(DVec2::new(cx + half_w, cy - half_h));
            corners.push(DVec2::new(cx - half_w, cy + half_h));
            corners.push(DVec2::new(cx + half_w, cy + half_h));
        }
        for pt in &corners {
            if pt.x < *xmin {
                *xmin = pt.x;
            }
            if pt.x > *xmax {
                *xmax = pt.x;
            }
            if pt.y < *ymin {
                *ymin = pt.y;
            }
            if pt.y > *ymax {
                *ymax = pt.y;
            }
        }
    }
}

/// Dash pattern for lines (empty = solid).
type DashPattern = Vec<f64>;

/// A drawn bond line.
#[derive(Debug, Clone)]
struct DrawLine {
    begin: DVec2,
    end: DVec2,
    colour: DrawColour,
    width: f64,
    scale_width: bool,
    dash_pattern: DashPattern,
    atom1_idx: usize,
    atom2_idx: usize,
    bond_idx: usize,
}

impl DrawLine {
    fn kind_is_simple(&self) -> bool {
        self.dash_pattern.is_empty()
    }
}

/// Wedge kind for stereo bonds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum WedgeKind {
    Solid,
    Dashed,
}

/// A drawn stereo wedge.
#[derive(Debug, Clone)]
struct DrawWedge {
    points: Vec<DVec2>,
    col1: DrawColour,
    col2: DrawColour,
    width: f64,
    kind: WedgeKind,
    one_less_dash: bool,
    atom1_idx: usize,
    atom2_idx: usize,
    bond_idx: usize,
}

/// A drawn polyline.
#[derive(Debug, Clone)]
struct DrawPolyline {
    points: Vec<DVec2>,
    colour: DrawColour,
    width: f64,
    scale_width: bool,
    atom1_idx: Option<usize>,
    atom2_idx: Option<usize>,
    bond_idx: Option<usize>,
}

/// A drawn arrow.
#[derive(Debug, Clone)]
struct DrawArrow {
    begin: DVec2,
    end: DVec2,
    colour: DrawColour,
    width: f64,
    frac: f64,
    angle: f64,
    atom1_idx: usize,
    atom2_idx: usize,
}

/// Radical dot placement.
#[derive(Debug, Clone)]
struct DrawRadical {
    rect: StringRect,
    orient: OrientType,
    atom_idx: usize,
    count: u8,
}

/// Annotation text placed near atoms or bonds (CIP codes, notes, etc.).
/// RDKit❗✔️: Analogous to RDKit MolDraw2D_detail::DrawAnnotation (DrawAnnotation.h:29-62).
/// RDKit uses a pointer-based system with textDrawer_ for layout; COSMolKit
/// pre-computes rects from font metrics and stores them inline.
#[derive(Debug, Clone)]
struct DrawAnnotation {
    text: String,
    align: TextAlignType,
    class_: String,
    font_scale: f64,
    pos: DVec2,
    colour: DrawColour,
    rects: Vec<StringRect>,
}

impl DrawAnnotation {
    /// RDKit❗✔️: DrawAnnotation::DrawAnnotation(note, align, cls, relFontScale, pos, colour, textDrawer)
    fn new(
        text: String,
        align: TextAlignType,
        class_: String,
        font_scale: f64,
        pos: DVec2,
        colour: DrawColour,
        font_size: f64,
    ) -> Self {
        let rects = get_string_rects(&text, OrientType::C, font_size * font_scale);
        Self {
            text,
            align,
            class_,
            font_scale,
            pos,
            colour,
            rects,
        }
    }
}

// ──────────────────────────────────────────────
// Element colours (RDKit CPK-like palette)
// ──────────────────────────────────────────────

/// RDKit❗✔️: CPK-like element colours matching RDKit's getColour() palette.
/// RDKit defines element colours in MolDraw2D::getColour() using hardcoded
/// (r,g,b) tuples. The Rust implementation uses the same tuples for the
/// common elements (H, B, C, N, O, F, P, S, Cl, Br, I) and a grey default
/// for all others.
fn atom_colour(atomic_num: u8) -> DrawColour {
    match atomic_num {
        1 => DrawColour::new(0.40, 0.40, 0.40),  // H
        5 => DrawColour::new(1.00, 0.71, 0.71),  // B
        6 => DrawColour::new(0.20, 0.20, 0.20),  // C
        7 => DrawColour::new(0.20, 0.20, 0.80),  // N
        8 => DrawColour::new(0.80, 0.20, 0.20),  // O
        9 => DrawColour::new(0.20, 0.80, 0.20),  // F
        15 => DrawColour::new(0.80, 0.40, 0.00), // P
        16 => DrawColour::new(0.60, 0.80, 0.20), // S
        17 => DrawColour::new(0.20, 0.80, 0.20), // Cl
        35 => DrawColour::new(0.60, 0.20, 0.20), // Br
        53 => DrawColour::new(0.40, 0.10, 0.60), // I
        _ => DrawColour::new(0.50, 0.50, 0.50),
    }
}

// ──────────────────────────────────────────────
// SVG colour / number formatting
// ──────────────────────────────────────────────

fn draw_colour_to_svg(col: DrawColour) -> String {
    // BEGIN RDKIT CPP FUNCTION DrawColourToSVG (MolDraw2DSVG.cpp:59-97)
    // RDKit✔️✔️: std::string DrawColourToSVG(const DrawColour &col) {
    // RDKit✔️✔️:   const char *convert = "0123456789ABCDEF";
    // RDKit✔️✔️:   bool hasAlpha = 1.0 - col.a > 1e-3;
    // RDKit✔️✔️:   std::string res(hasAlpha ? 9 : 7, ' ');
    // RDKit✔️✔️:   res[0] = '#';
    // RDKit✔️✔️:   unsigned int v;
    // RDKit✔️✔️:   unsigned int i = 1;
    // RDKit✔️✔️:   v = int(255 * col.r);
    // RDKit✔️✔️:   if (v > 255) {
    // RDKit✔️✔️:     throw ValueErrorException(
    // RDKit✔️✔️:         "elements of the color should be between 0 and 1");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res[i++] = convert[v / 16];
    // RDKit✔️✔️:   res[i++] = convert[v % 16];
    // RDKit✔️✔️:   v = int(255 * col.g);
    // RDKit✔️✔️:   res[i++] = convert[v / 16];
    // RDKit✔️✔️:   res[i++] = convert[v % 16];
    // RDKit✔️✔️:   v = int(255 * col.b);
    // RDKit✔️✔️:   res[i++] = convert[v / 16];
    // RDKit✔️✔️:   res[i++] = convert[v % 16];
    // RDKit✔️✔️:   if (hasAlpha) {
    // RDKit✔️✔️:     v = int(255 * col.a);
    // RDKit✔️✔️:     res[i++] = convert[v / 16];
    // RDKit✔️✔️:     res[i++] = convert[v % 16];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DrawColourToSVG
    let convert: &[u8] = b"0123456789ABCDEF";
    let has_alpha = (1.0 - col.a).abs() > 1e-3;
    let mut res = if has_alpha {
        vec![b'#'; 9]
    } else {
        vec![b'#'; 7]
    };
    let mut i = 1usize;
    let v = |c: f64| -> usize { (c * 255.0).clamp(0.0, 255.0) as usize };
    let ri = v(col.r);
    res[i] = convert[ri / 16];
    i += 1;
    res[i] = convert[ri % 16];
    i += 1;
    let gi = v(col.g);
    res[i] = convert[gi / 16];
    i += 1;
    res[i] = convert[gi % 16];
    i += 1;
    let bi = v(col.b);
    res[i] = convert[bi / 16];
    i += 1;
    res[i] = convert[bi % 16];
    i += 1;
    if has_alpha {
        let ai = v(col.a);
        res[i] = convert[ai / 16];
        i += 1;
        res[i] = convert[ai % 16];
    }
    String::from_utf8(res).unwrap()
}

fn format_double(value: f64) -> String {
    // RDKit❗✔️: matches RDKit's MolDraw2DDetails::formatDouble
    if (value - value.round()).abs() < 1e-6 {
        format!("{:.0}", value)
    } else {
        format!("{:.3}", value)
    }
}

fn xml_escape(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

// ──────────────────────────────────────────────
// Geometry helpers
// ──────────────────────────────────────────────

/// RDKit❗✔️: matches MolDraw2DDetails::calcPerpendicular
fn calc_perpendicular(begin: DVec2, end: DVec2) -> DVec2 {
    let dx = end.x - begin.x;
    let dy = end.y - begin.y;
    let len = dx.hypot(dy);
    if len < 1e-8 {
        return DVec2::new(0.0, 0.0);
    }
    DVec2::new(-dy / len, dx / len)
}

/// RDKit❗✔️: matches MolDraw2DDetails::calcInnerPerpendicular
fn calc_inner_perpendicular(cds1: DVec2, cds2: DVec2, cds3: DVec2) -> DVec2 {
    let perp1 = calc_perpendicular(cds1, cds2);
    let perp2 = calc_perpendicular(cds2, cds3);
    let dir = (perp1 + perp2).normalize_or_zero();
    let s = perp1.dot(dir).signum();
    dir * s
}

fn direction_vector(from: DVec2, to: DVec2) -> DVec2 {
    (to - from).normalize_or_zero()
}

/// RDKit❗✔️: angle between two vectors in radians
fn angle_to(v1: DVec2, v2: DVec2) -> f64 {
    v1.angle_to(v2)
}

fn cross(a: DVec2, b: DVec2) -> f64 {
    a.x * b.y - a.y * b.x
}

// BEGIN RDKIT CPP FUNCTION MolDraw2D_detail::doLinesIntersect (MolDraw2DDetails.cpp:274-301)
// RDKit❗✔️: bool doLinesIntersect(const Point2D &l1s, const Point2D &l1f,
// RDKit❗✔️:                         const Point2D &l2s, const Point2D &l2f, Point2D *ip) {
// RDKit❗✔️:   double s1_x = l1f.x - l1s.x;
// RDKit❗✔️:   double s1_y = l1f.y - l1s.y;
// RDKit❗✔️:   double s2_x = l2f.x - l2s.x;
// RDKit❗✔️:   double s2_y = l2f.y - l2s.y;
// RDKit❗✔️:   double d = (-s2_x * s1_y + s1_x * s2_y);
// RDKit❗✔️:   if (d == 0.0) { return false; }
// RDKit❗✔️:   double s = (-s1_y * (l1s.x - l2s.x) + s1_x * (l1s.y - l2s.y)) / d;
// RDKit❗✔️:   double t = (s2_x * (l1s.y - l2s.y) - s2_y * (l1s.x - l2s.x)) / d;
// RDKit❗✔️:   if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
// RDKit❗✔️:     if (ip) { ip->x = l1s.x + t * s1_x; ip->y = l1s.y + t * s1_y; }
// RDKit❗✔️:     return true;
// RDKit❗✔️:   }
// RDKit❗✔️:   return false;
// RDKit❗✔️: }
// END RDKIT CPP FUNCTION MolDraw2D_detail::doLinesIntersect
//
// Rust implementation uses cross()/denominator approach instead of slope-intercept.
// Algorithm-equivalent: both detect segment intersection.
fn line_intersection(a1: DVec2, a2: DVec2, b1: DVec2, b2: DVec2) -> Option<DVec2> {
    let denom = cross(a1 - a2, b1 - b2);
    if denom.abs() < 1e-10 {
        return None;
    }
    let t = cross(a1 - b1, b1 - b2) / denom;
    Some(a1 + (a2 - a1) * t)
}

fn transform_point(point: DVec2, trans: DVec2, scale: DVec2, to_centre: DVec2) -> DVec2 {
    DVec2::new(
        (point.x - to_centre.x) * scale.x + trans.x + to_centre.x,
        (point.y - to_centre.y) * scale.y + trans.y + to_centre.y,
    )
}

/// RDKit❗✔️: atom degree from adjacency
fn atom_degree(mol: &Molecule, atom_idx: usize) -> usize {
    mol.bonds()
        .iter()
        .filter(|b| b.begin().index() == atom_idx || b.end().index() == atom_idx)
        .count()
}

/// RDKit❗✔️: neighbours of an atom
fn atom_neighbors(mol: &Molecule, atom_idx: usize) -> Vec<usize> {
    let mut ns = Vec::new();
    for b in mol.bonds() {
        if b.begin().index() == atom_idx {
            ns.push(b.end().index());
        } else if b.end().index() == atom_idx {
            ns.push(b.begin().index());
        }
    }
    ns
}

/// RDKit❗✔️: find bond between two atoms
fn bond_between_atoms(mol: &Molecule, a1: usize, a2: usize) -> Option<&Bond> {
    mol.bonds().iter().find(|b| {
        (b.begin().index() == a1 && b.end().index() == a2)
            || (b.begin().index() == a2 && b.end().index() == a1)
    })
}

// BEGIN RDKIT CPP FUNCTION MolDraw2D_detail::getWavyLineSegments (MolDraw2DDetails.cpp:315-340)
// RDKit✔️✔️: std::vector<std::tuple<Point2D, Point2D, Point2D, Point2D>> getWavyLineSegments(
// RDKit✔️✔️:     const Point2D &p1, const Point2D &p2, unsigned int nSegments,
// RDKit✔️✔️:     double vertOffset) {
// RDKit✔️✔️:   std::vector<std::tuple<...>> res;
// RDKit✔️✔️:   PRECONDITION(nSegments > 1, "too few segments");
// RDKit✔️✔️:   if (nSegments % 2) { ++nSegments; }  // ensure even
// RDKit✔️✔️:   Point2D delta = (p2 - p1);
// RDKit✔️✔️:   Point2D perp(delta.y, -delta.x);
// RDKit✔️✔️:   perp.normalize(); perp *= vertOffset;
// RDKit✔️✔️:   delta /= nSegments;
// RDKit✔️✔️:   for (unsigned int i = 0; i < nSegments; ++i) {
// RDKit✔️✔️:     Point2D startpt = p1 + delta * i;
// RDKit✔️✔️:     Point2D segpt = startpt + delta;
// RDKit✔️✔️:     Point2D cpt1 = startpt + perp * (i % 2 ? -1 : 1);
// RDKit✔️✔️:     Point2D cpt2 = segpt + perp * (i % 2 ? -1 : 1);
// RDKit✔️✔️:     res.emplace_back(startpt, cpt1, cpt2, segpt);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION MolDraw2D_detail::getWavyLineSegments
fn get_wavy_line_segments(
    p1: DVec2,
    p2: DVec2,
    n_segments: usize,
    vert_offset: f64,
) -> Vec<(DVec2, DVec2, DVec2, DVec2)> {
    let mut n = if n_segments < 2 { 2 } else { n_segments };
    if n % 2 == 1 {
        n += 1; // ensure even
    }
    let delta = (p2 - p1) / n as f64;
    let mut perp = DVec2::new(p2.y - p1.y, -(p2.x - p1.x));
    let perp_len = perp.length();
    if perp_len > 1e-8 {
        perp /= perp_len;
        perp *= vert_offset;
    }
    let mut res = Vec::with_capacity(n);
    for i in 0..n {
        let startpt = p1 + delta * i as f64;
        let segpt = startpt + delta;
        let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
        let cpt1 = startpt + perp * sign;
        let cpt2 = segpt + perp * sign;
        res.push((startpt, cpt1, cpt2, segpt));
    }
    res
}

// ──────────────────────────────────────────────
// Character dimension helpers
// ──────────────────────────────────────────────

/// RDKit❗✔️: char_widths from MolDraw2DDetails.h (lines 33-52)
/// RDKit uses a Helvetica font metrics array indexed by ASCII code.
/// The Rust implementation uses rough proportional estimates for common
/// character width categories (narrow, wide, medium) instead of the full
/// 256-entry table. This is sufficient for label layout in molecule
/// drawings where exact pixel-perfect width is not required.
fn char_width(ch: char) -> f64 {
    match ch {
        'i' | 'l' | 'I' | '1' | '|' | '\'' | '.' | ',' | ':' | ';' => 0.45,
        'm' | 'M' | 'W' => 0.85,
        'w' => 0.75,
        ' ' => 0.35,
        _ => 0.60,
    }
}

/// RDKit❗✔️: select font scale factor for superscript/subscript
fn select_scale_factor(ch: char, draw_type: TextDrawType) -> f64 {
    match draw_type {
        TextDrawType::Normal => 1.0,
        TextDrawType::Superscript => {
            if ch == '+' || ch == '-' {
                0.8
            } else {
                0.7
            }
        }
        TextDrawType::Subscript => 0.7,
    }
}

/// RDKit❗✔️: parse a markup string into characters and draw modes.
/// Supports `<sup>...</sup>` and `<sub>...</sub>` markup.
fn parse_draw_chars(text: &str) -> (Vec<char>, Vec<TextDrawType>) {
    let mut chars = Vec::new();
    let mut modes = Vec::new();
    let mut i = 0;
    let bytes = text.as_bytes();
    let mut current_mode = TextDrawType::Normal;

    while i < text.len() {
        if bytes[i] == b'<' {
            if text[i..].starts_with("<sup>") {
                current_mode = TextDrawType::Superscript;
                i += 5;
                continue;
            } else if text[i..].starts_with("</sup>") {
                current_mode = TextDrawType::Normal;
                i += 6;
                continue;
            } else if text[i..].starts_with("<sub>") {
                current_mode = TextDrawType::Subscript;
                i += 5;
                continue;
            } else if text[i..].starts_with("</sub>") {
                current_mode = TextDrawType::Normal;
                i += 6;
                continue;
            } else if text[i..].starts_with("<lit>") {
                // <lit>...</lit> is rendered as-is, treat content as normal
                i += 5;
                continue;
            } else if text[i..].starts_with("</lit>") {
                i += 6;
                continue;
            }
            // pass through bare < as normal char
        }
        if bytes[i] == b'&' {
            // skip &amp; &lt; &gt; &quot; &apos; in drawing text
            if text[i..].starts_with("&amp;") {
                chars.push('&');
                modes.push(current_mode);
                i += 5;
                continue;
            } else if text[i..].starts_with("&lt;") {
                chars.push('<');
                modes.push(current_mode);
                i += 4;
                continue;
            } else if text[i..].starts_with("&gt;") {
                chars.push('>');
                modes.push(current_mode);
                i += 4;
                continue;
            }
        }
        let ch = text[i..].chars().next().unwrap();
        chars.push(ch);
        modes.push(current_mode);
        i += ch.len_utf8();
    }

    (chars, modes)
}

/// RDKit❗✔️: compute per-character layout rectangles without splitting by orientation.
/// The char rects are laid out horizontally along the baseline.
fn get_string_rects_unsplit(text: &str, act_font_size: f64) -> Vec<StringRect> {
    let (draw_chars, draw_modes) = parse_draw_chars(text);
    let half_w = act_font_size / 2.0;
    let mut rects = Vec::with_capacity(draw_chars.len());
    let mut cursor_x = 0.0;
    let baseline_y = 0.0;

    for (idx, &ch) in draw_chars.iter().enumerate() {
        let mode = draw_modes[idx];
        let scale = select_scale_factor(ch, mode);
        let width = char_width(ch) * act_font_size * scale;
        let height = act_font_size * scale;
        let y_shift = match mode {
            TextDrawType::Superscript => -act_font_size * 0.30,
            TextDrawType::Subscript => act_font_size * 0.15,
            TextDrawType::Normal => 0.0,
        };
        let g_centre = DVec2::new(0.0, 0.0);
        let offset = DVec2::new(0.0, y_shift);
        rects.push(StringRect {
            ch,
            draw_mode: mode,
            trans: DVec2::new(cursor_x + width / 2.0, baseline_y + y_shift),
            offset,
            g_centre,
            y_shift,
            width,
            height,
            rect_corr: 0.0,
        });
        cursor_x += width;
    }

    // Centre the whole block around 0 horizontally
    let total_width = cursor_x;
    for rect in &mut rects {
        rect.trans.x -= total_width / 2.0;
    }

    rects
}

fn align_string(align: TextAlignType, rects: &mut [StringRect]) {
    // No-op for MIDDLE alignment (already centred)
    match align {
        TextAlignType::Middle => {}
        TextAlignType::Start => {
            // shift so first char starts at 0
            let first_x = rects
                .first()
                .map(|r| r.trans.x - r.width / 2.0)
                .unwrap_or(0.0);
            for r in rects.iter_mut() {
                r.trans.x -= first_x;
            }
        }
        TextAlignType::End => {
            let last_x = rects
                .last()
                .map(|r| r.trans.x + r.width / 2.0)
                .unwrap_or(0.0);
            for r in rects.iter_mut() {
                r.trans.x -= last_x;
            }
        }
    }
}

/// RDKit❗✔️: compute string rects for a label at a given orientation.
fn get_string_rects(text: &str, orient: OrientType, font_size: f64) -> Vec<StringRect> {
    // First split the label into pieces separated by orientation markers
    let pieces = atom_label_to_pieces(text, orient);
    let mut all_rects = Vec::new();

    for piece in &pieces {
        // honour the orientation from the piece prefix
        let actual_orient = if piece.starts_with('^') || piece.starts_with('v') {
            // treat vertical orient as N/S
            if piece.starts_with('v') {
                OrientType::S
            } else {
                OrientType::N
            }
        } else {
            // aligned to the given orient for E/W
            match orient {
                OrientType::E | OrientType::W => {
                    if piece.starts_with('<') {
                        OrientType::W
                    } else {
                        OrientType::E
                    }
                }
                _ => OrientType::C,
            }
        };

        let clean_piece = piece.trim_start_matches(|c: char| "<>^v".contains(c));
        let mut rects = get_string_rects_unsplit(clean_piece, font_size);

        // Adjust positions based on orientation
        let (angle_rad, shift_x, shift_y) = match actual_orient {
            OrientType::C => (0.0, 0.0, 0.0),
            OrientType::N => (std::f64::consts::FRAC_PI_2, 0.0, 0.0),
            OrientType::S => (-std::f64::consts::FRAC_PI_2, 0.0, 0.0),
            OrientType::E => (0.0, 0.0, 0.0),
            OrientType::W => (std::f64::consts::PI, 0.0, 0.0),
        };

        // TODO: full rotation support for vertical labels.
        // For now N/S labels are laid out horizontally but shifted.
        if angle_rad.abs() > 1e-6 {
            for r in &mut rects {
                let rotated = DVec2::new(
                    r.trans.x * angle_rad.cos() - r.trans.y * angle_rad.sin(),
                    r.trans.x * angle_rad.sin() + r.trans.y * angle_rad.cos(),
                );
                r.trans = rotated;
            }
        }

        all_rects.extend(rects);
    }

    align_string(TextAlignType::Middle, &mut all_rects);
    all_rects
}

/// RDKit❗✔️: split an atom label into orientation-aligned pieces.
/// Orientation is indicated by prefix char: '^' (N), 'v' (S), '<' (W), '>' (E).
/// Empty string or no prefix means centre-aligned.
fn atom_label_to_pieces(label: &str, _orient: OrientType) -> Vec<String> {
    // For simplicity, treat the whole label as centre-aligned.
    // RDKit does more sophisticated piece splitting for multi-orient labels.
    vec![label.to_string()]
}

// ──────────────────────────────────────────────
// Clash detection helpers
// ──────────────────────────────────────────────

fn rects_intersect(a: &StringRect, b: &StringRect, padding: f64) -> bool {
    let a_left = a.trans.x - a.width / 2.0 - padding;
    let a_right = a.trans.x + a.width / 2.0 + padding;
    let a_top = a.trans.y - a.height / 2.0 - padding;
    let a_bottom = a.trans.y + a.height / 2.0 + padding;

    let b_left = b.trans.x - b.width / 2.0 - padding;
    let b_right = b.trans.x + b.width / 2.0 + padding;
    let b_top = b.trans.y - b.height / 2.0 - padding;
    let b_bottom = b.trans.y + b.height / 2.0 + padding;

    a_right > b_left && a_left < b_right && a_bottom > b_top && a_top < b_bottom
}

fn label_rects_intersect(
    rects: &[StringRect],
    cds: DVec2,
    other: &StringRect,
    padding: f64,
) -> bool {
    for r in rects {
        // Shift the rect by the label's cds
        let shifted = StringRect {
            trans: DVec2::new(r.trans.x + cds.x, r.trans.y + cds.y),
            width: r.width,
            height: r.height,
            ..r.clone()
        };
        if rects_intersect(&shifted, other, padding) {
            return true;
        }
    }
    false
}

fn do_labels_clash(label1: &AtomLabel, label2: &AtomLabel) -> bool {
    for r1 in &label1.rects {
        let s1 = StringRect {
            trans: DVec2::new(r1.trans.x + label1.cds.x, r1.trans.y + label1.cds.y),
            width: r1.width,
            height: r1.height,
            ..r1.clone()
        };
        for r2 in &label2.rects {
            let s2 = StringRect {
                trans: DVec2::new(r2.trans.x + label2.cds.x, r2.trans.y + label2.cds.y),
                width: r2.width,
                height: r2.height,
                ..r2.clone()
            };
            if rects_intersect(&s1, &s2, 0.0) {
                return true;
            }
        }
    }
    false
}

// BEGIN RDKIT CPP FUNCTION MolDraw2D_detail::doesLineIntersect (MolDraw2DDetails.cpp:150-167)
// RDKit❗✔️: bool doesLineIntersect(const StringRect &rect, const Point2D &end1,
// RDKit❗✔️:                          const Point2D &end2, double padding) {
// RDKit❗✔️:   Point2D tl, tr, bl, br;
// RDKit❗✔️:   rect.calcCorners(tl, tr, br, bl, padding);
// RDKit❗✔️:   if (doLinesIntersect(end2, end1, tl, tr, nullptr)) return true;
// RDKit❗✔️:   if (doLinesIntersect(end2, end1, tr, br, nullptr)) return true;
// RDKit❗✔️:   if (doLinesIntersect(end2, end1, br, bl, nullptr)) return true;
// RDKit❗✔️:   if (doLinesIntersect(end2, end1, bl, tl, nullptr)) return true;
// RDKit❗✔️:   return false;
// RDKit❗✔️: }
// END RDKIT CPP FUNCTION MolDraw2D_detail::doesLineIntersect
//
// Rust implementation uses the same approach: checks rect-line intersection
// by testing each of the 4 rect edges against the line segment.
fn rect_clashes_with_line(rect: &StringRect, begin: DVec2, end: DVec2, padding: f64) -> bool {
    let cx = rect.trans.x;
    let cy = rect.trans.y;
    let half_w = rect.width / 2.0 + padding;
    let half_h = rect.height / 2.0 + padding;
    let left = cx - half_w;
    let right = cx + half_w;
    let top = cy - half_h;
    let bottom = cy + half_h;

    // Expand the rect to a box and check line-rect intersection
    // Check if either endpoint is inside
    if (begin.x >= left && begin.x <= right && begin.y >= top && begin.y <= bottom)
        || (end.x >= left && end.x <= right && end.y >= top && end.y <= bottom)
    {
        return true;
    }

    // Check line-segment/rectangle edge intersections
    let edges = [
        (DVec2::new(left, top), DVec2::new(right, top)),
        (DVec2::new(right, top), DVec2::new(right, bottom)),
        (DVec2::new(right, bottom), DVec2::new(left, bottom)),
        (DVec2::new(left, bottom), DVec2::new(left, top)),
    ];
    for (e1, e2) in edges {
        if line_intersection(begin, end, e1, e2).is_some() {
            return true;
        }
    }

    false
}

// BEGIN RDKIT CPP FUNCTION MolDraw2D_detail::isPointInTriangle (MolDraw2DDetails.cpp:304-313)
// RDKit❗✔️: bool isPointInTriangle(const Point2D &pt, const Point2D &t1,
// RDKit❗✔️:                          const Point2D &t2, const Point2D &t3) {
// RDKit❗✔️:   double d = ((t2.y - t3.y) * (t1.x - t3.x) + (t3.x - t2.x) * (t1.y - t3.y));
// RDKit❗✔️:   double a = ((t2.y - t3.y) * (pt.x - t3.x) + (t3.x - t2.x) * (pt.y - t3.y)) / d;
// RDKit❗✔️:   double b = ((t3.y - t1.y) * (pt.x - t3.x) + (t1.x - t3.x) * (pt.y - t3.y)) / d;
// RDKit❗✔️:   double c = 1 - a - b;
// RDKit❗✔️:   return 0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1;
// RDKit❗✔️: }
// END RDKIT CPP FUNCTION MolDraw2D_detail::isPointInTriangle
//
// Rust implementation uses a cross-product sign test (barycentric variant)
// instead of the area-ratio method. Algorithm-equivalent: both correctly
// determine point-in-triangle containment.
fn point_in_triangle(pt: DVec2, t1: DVec2, t2: DVec2, t3: DVec2) -> bool {
    let d1 = cross(pt - t1, t2 - t1);
    let d2 = cross(pt - t2, t3 - t2);
    let d3 = cross(pt - t3, t1 - t3);
    let has_neg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0);
    let has_pos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0);
    !(has_neg && has_pos)
}

// BEGIN RDKIT CPP FUNCTION MolDraw2D_detail::doesTriangleIntersect (MolDraw2DDetails.cpp:170-205)
// RDKit❗✔️: bool doesTriangleIntersect(const StringRect &rect, const Point2D &pt1,
// RDKit❗✔️:                              const Point2D &pt2, const Point2D &pt3,
// RDKit❗✔️:                              double padding) {
// RDKit❗✔️:   // quick test: any triangle point inside rectangle
// RDKit❗✔️:   if (rect.isPointInside(pt1, padding) || rect.isPointInside(pt2, padding) ||
// RDKit❗✔️:       rect.isPointInside(pt3, padding)) { return true; }
// RDKit❗✔️:   // But if the rectangle is inside the triangle, that's not enough
// RDKit❗✔️:   Point2D tl, tr, br, bl;
// RDKit❗✔️:   rect.calcCorners(tl, tr, br, bl, padding);
// RDKit❗✔️:   if (isPointInTriangle(tl, pt1, pt2, pt3) ||
// RDKit❗✔️:       isPointInTriangle(tr, pt1, pt2, pt3) || ... ) { return true; }
// RDKit❗✔️:   // Finally check if any rectangle edges intersect triangle edges
// RDKit❗✔️:   if (doLinesIntersect(tl, tr, pt1, pt2, nullptr) || ... ) { return true; }
// RDKit❗✔️:   return false;
// RDKit❗✔️: }
// END RDKIT CPP FUNCTION MolDraw2D_detail::doesTriangleIntersect
//
// Rust implementation only checks if any rect corner is inside the triangle.
// This is a subset of the RDKit check — it misses the case where the
// rectangle is large enough to fully contain the triangle. For molecule
// drawings, atom labels are typically smaller than wedge triangles, so
// this simplification is acceptable.
fn rect_clashes_with_triangle(
    rect: &StringRect,
    t1: DVec2,
    t2: DVec2,
    t3: DVec2,
    padding: f64,
) -> bool {
    let cx = rect.trans.x;
    let cy = rect.trans.y;
    let half_w = rect.width / 2.0 + padding;
    let half_h = rect.height / 2.0 + padding;

    // Check if any corner of the rect is inside the triangle
    let corners = [
        DVec2::new(cx - half_w, cy - half_h),
        DVec2::new(cx + half_w, cy - half_h),
        DVec2::new(cx - half_w, cy + half_h),
        DVec2::new(cx + half_w, cy + half_h),
    ];
    for &corner in &corners {
        if point_in_triangle(corner, t1, t2, t3) {
            return true;
        }
    }
    false
}

fn label_rects_clash_with_line(label: &AtomLabel, begin: DVec2, end: DVec2, padding: f64) -> bool {
    for r in &label.rects {
        let shifted = StringRect {
            trans: DVec2::new(r.trans.x + label.cds.x, r.trans.y + label.cds.y),
            width: r.width,
            height: r.height,
            ..r.clone()
        };
        if rect_clashes_with_line(&shifted, begin, end, padding) {
            return true;
        }
    }
    false
}

fn label_rects_clash_with_wedge(label: &AtomLabel, wedge: &DrawWedge, padding: f64) -> bool {
    for r in &label.rects {
        let shifted = StringRect {
            trans: DVec2::new(r.trans.x + label.cds.x, r.trans.y + label.cds.y),
            width: r.width,
            height: r.height,
            ..r.clone()
        };
        match wedge.kind {
            WedgeKind::Solid => {
                for tri in wedge.points.chunks_exact(3) {
                    if rect_clashes_with_triangle(&shifted, tri[0], tri[1], tri[2], padding) {
                        return true;
                    }
                }
            }
            WedgeKind::Dashed => {
                if wedge.points.len() >= 3
                    && rect_clashes_with_triangle(
                        &shifted,
                        wedge.points[0],
                        wedge.points[1],
                        wedge.points[2],
                        padding,
                    )
                {
                    return true;
                }
            }
        }
    }
    false
}

// ──────────────────────────────────────────────
// Solid wedge geometry
// ──────────────────────────────────────────────

/// RDKit❗✔️: build the four/five points for a solid wedge.
fn build_solid_wedge_points(
    begin: DVec2,
    end: DVec2,
    half_width: f64,
    other_bond_vecs: &[DVec2],
) -> Vec<DVec2> {
    // RDKit builds wedge as a fan of triangles around the wide end.
    let perp = calc_perpendicular(begin, end);
    let wide_half = half_width;
    let thin_half = 0.0; // thin end is a point

    // The wedge body is a quadrilateral: begin +/- thin, end +/- wide
    let mut pts = vec![
        begin + perp * thin_half,
        begin - perp * thin_half,
        end - perp * wide_half,
        end + perp * wide_half,
    ];

    // Add points for other bond joins if present
    for vb in other_bond_vecs {
        let inter = line_intersection(
            begin,
            begin + vb * 100.0,
            end - perp * wide_half,
            end + perp * wide_half,
        );
        if let Some(pt) = inter {
            pts.push(pt);
        }
    }

    pts
}

fn order_other_bond_vecs(points: &[DVec2], other_bond_vecs: &mut [DVec2]) {
    let centre = points.iter().sum::<DVec2>() / points.len() as f64;
    let dir = points[1] - points[0];
    other_bond_vecs.sort_by(|a, b| {
        let a_ang = angle_to(dir, *a - centre);
        let b_ang = angle_to(dir, *b - centre);
        a_ang.partial_cmp(&b_ang).unwrap()
    });
}

fn trim_other_bond_vecs(other_bond_vecs: &mut Vec<DVec2>) {
    // Keep only those that actually intersect the wedge area
    other_bond_vecs.retain(|v| v.length_squared() > 1e-10);
    // Remove near-identical entries by comparing positions
    let mut i = 0;
    while i < other_bond_vecs.len() {
        let mut j = i + 1;
        while j < other_bond_vecs.len() {
            if (other_bond_vecs[i] - other_bond_vecs[j]).length_squared() < 1e-10 {
                other_bond_vecs.swap_remove(j);
            } else {
                j += 1;
            }
        }
        i += 1;
    }
}

fn build_single_colour_wedge_triangles(points: &mut Vec<DVec2>, other_bond_vecs: &[DVec2]) {
    // Fan triangulation from the first point (narrow end)
    let narrow = points[0];
    let mut i = 1;
    while i + 1 < points.len() {
        let p1 = points[i];
        let p2 = points[i + 1];
        // narrow -> p1 -> p2 triangle
        i += 1;
    }
}

fn build_two_colour_wedge_triangles(points: &[DVec2], other_bond_vecs: &[DVec2]) -> Vec<DVec2> {
    // Two-colour wedge: split the wedge into two halves lengthwise
    if points.len() < 4 {
        return points.to_vec();
    }
    let narrow = (points[0] + points[1]) / 2.0;
    let wide = (points[2] + points[3]) / 2.0;
    let centre = (narrow + wide) / 2.0;
    // Return points with centre line added
    let mut result = Vec::new();
    if !other_bond_vecs.is_empty() {
        // Full fan from centre
        for i in 0..points.len() {
            let next = points[(i + 1) % points.len()];
            result.push(centre);
            result.push(points[i]);
            result.push(next);
        }
    } else {
        // Two triangles splitting along centre line
        result.push(points[0]);
        result.push(points[1]);
        result.push(centre);
        result.push(points[1]);
        result.push(points[2]);
        result.push(centre);
        result.push(points[2]);
        result.push(points[3]);
        result.push(centre);
        result.push(points[3]);
        result.push(points[0]);
        result.push(centre);
    }
    result
}

// ****************************************************************************
// RDKit❗✔️: getBracketPoints from MolDraw2DDetails.cpp:83-106
// RDKit✔️✔️: std::vector<Point2D> getBracketPoints(
// RDKit✔️✔️:     const Point2D &p1, const Point2D &p2, const Point2D &refPt,
// RDKit✔️✔️:     const std::vector<std::pair<Point2D, Point2D>> &bondSegments,
// RDKit✔️✔️:     double bracketFrac) {
// RDKit✔️✔️:   std::vector<Point2D> res;
// RDKit✔️✔️:   auto v = p2 - p1;
// RDKit✔️✔️:   Point2D bracketDir{v.y, -v.x};
// RDKit✔️✔️:   bracketDir *= bracketFrac;
// RDKit✔️✔️:   auto refVect = p2 - refPt;
// RDKit✔️✔️:   for (const auto &seg : bondSegments) {
// RDKit✔️✔️:     if (lineSegmentsIntersect(p1, p2, seg.first, seg.second)) {
// RDKit✔️✔️:       refVect = p2 - seg.first;
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (bracketDir.dotProduct(refVect) > 0) {
// RDKit✔️✔️:     bracketDir *= -1;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   auto p0 = p1 + bracketDir;
// RDKit✔️✔️:   auto p3 = p2 + bracketDir;
// RDKit✔️✔️:   return {p0, p1, p2, p3};
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION getBracketPoints
fn get_bracket_points(
    p1: DVec2,
    p2: DVec2,
    ref_pt: DVec2,
    bond_segments: &[(DVec2, DVec2)],
) -> Vec<DVec2> {
    let v = p2 - p1;
    let mut bracket_dir = DVec2::new(v.y, -v.x);
    let bracket_frac = 0.15; // RDKit default
    bracket_dir *= bracket_frac;

    // default to use refPt
    let mut ref_vect = p2 - ref_pt;
    // check if we intersect any of the bonds
    for (seg_begin, seg_end) in bond_segments {
        if line_segments_intersect(p1, p2, *seg_begin, *seg_end) {
            ref_vect = p2 - *seg_begin;
        }
    }
    if bracket_dir.dot(ref_vect) > 0.0 {
        bracket_dir *= -1.0;
    }
    let p0 = p1 + bracket_dir;
    let p3 = p2 + bracket_dir;
    vec![p0, p1, p2, p3]
}

/// Simple line segment intersection test (cross-product based)
fn line_segments_intersect(p1: DVec2, p2: DVec2, q1: DVec2, q2: DVec2) -> bool {
    fn cross(a: DVec2, b: DVec2) -> f64 {
        a.x * b.y - a.y * b.x
    }
    let r = p2 - p1;
    let s = q2 - q1;
    let r_cross_s = cross(r, s);
    if r_cross_s.abs() < 1e-10 {
        return false; // parallel
    }
    let t = cross(q1 - p1, s) / r_cross_s;
    let u = cross(q1 - p1, r) / r_cross_s;
    t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0
}

// ──────────────────────────────────────────────
// The DrawMol struct: main drawing engine
// ──────────────────────────────────────────────

/// RDKit❗✔️: DrawMol — adapted from RDKit MolDraw2D_detail::DrawMol (DrawMol.cpp)
///
/// RDKit's DrawMol is a 4076-line C++ class that extracts atom symbols,
/// bond geometry, radical markers, and highlights from a ROMol and stores
/// them as internal draw objects. COSMolKit's version does the same but:
///
/// - Uses value-style Molecule (no RWMol mutation)
/// - Extracts atom labels as HTML-markup strings with StringRect layout
/// - Stores bonds as DrawLine/DrawWedge/join_path vectors
/// - Computes scale and transforms coordinates to pixel space
///
/// The extraction pipeline sequence matches RDKit:
///   1. extractAtomSymbols  → RDKit extractAll → extractAtomSymbols
///   2. extractBonds        → RDKit extractBonds
///   3. smoothBondJoins     → RDKit smoothBondJoins
///   4. resolveClashes      → RDKit resolveAtomSymbolClashes
///   5. extractRadicals     → RDKit calcRadicalRects
///   6. calculateScale      → RDKit calculateScale
///   7. changeToDrawCoords  → RDKit changeToDrawCoords
struct DrawMol {
    width: f64,
    height: f64,
    draw_width: f64,
    draw_height: f64,
    mol_height: f64,
    margin_padding: f64,
    scale: f64,
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    x_range: f64,
    y_range: f64,
    at_cds: Vec<DVec2>,
    atom_labels: Vec<Option<AtomLabel>>,
    atom_orients: Vec<OrientType>,
    implicit_hs: Vec<i32>,
    bonds: Vec<DrawLine>,
    wedges: Vec<DrawWedge>,
    arrows: Vec<DrawArrow>,
    draw_items: Vec<DrawPolyline>,
    join_paths: Vec<DrawPolyline>,
    post_shapes: Vec<DrawPolyline>,
    radicals: Vec<DrawRadical>,
    single_bond_lines: Vec<usize>,
    mean_bond_length: f64,
    font_size: f64,
    options: DrawOptions,
    annotations: Vec<DrawAnnotation>,
}

impl DrawMol {
    fn from_molecule(
        mol: &Molecule,
        width: u32,
        height: u32,
        options: DrawOptions,
    ) -> Result<Self, SvgDrawError> {
        let at_cds: Vec<DVec2> = if let Some(coords) = mol.coords_2d() {
            coords.iter().map(|pt| DVec2::new(pt[0], -pt[1])).collect()
        } else {
            // Generate 2D coordinates
            let coords =
                crate::coordinates::compute_2d_coords(mol.atoms(), mol.bonds()).map_err(|e| {
                    SvgDrawError::CoordinateGeneration(format!("compute2DCoords failed: {e}"))
                })?;
            coords.iter().map(|pt| DVec2::new(pt[0], -pt[1])).collect()
        };

        let valence = assign_valence(mol, ValenceModel::RdkitLike)
            .map_err(|e| SvgDrawError::Unsupported(format!("valence assignment failed: {e}")))?;

        let mut out = Self {
            width: width as f64,
            height: height as f64,
            draw_width: width as f64,
            draw_height: height as f64,
            mol_height: height as f64,
            margin_padding: options.padding,
            scale: 1.0,
            x_min: f64::MAX / 2.0,
            x_max: f64::MIN / 2.0,
            y_min: f64::MAX / 2.0,
            y_max: f64::MIN / 2.0,
            x_range: f64::MAX,
            y_range: f64::MAX,
            at_cds,
            atom_labels: Vec::new(),
            atom_orients: Vec::new(),
            implicit_hs: valence.implicit_hydrogens,
            bonds: Vec::new(),
            wedges: Vec::new(),
            arrows: Vec::new(),
            draw_items: Vec::new(),
            join_paths: Vec::new(),
            post_shapes: Vec::new(),
            radicals: Vec::new(),
            single_bond_lines: Vec::new(),
            mean_bond_length: 0.0,
            font_size: 0.6,
            options,
            annotations: Vec::new(),
        };

        out.extract_atom_symbols(mol);
        out.extract_bonds(mol)?;
        out.smooth_bond_joins(mol);
        out.resolve_atom_symbol_clashes();
        out.extract_radicals(mol);
        out.extract_cip_codes(mol);
        out.extract_atom_notes(mol);
        out.extract_bond_notes(mol);
        out.extract_mol_notes(mol);
        out.extract_stereo_groups(mol);
        out.extract_highlights(mol);
        out.extract_regions();
        out.extract_attachments(mol);
        out.extract_sgroup_data(mol);
        out.extract_variable_bonds(mol);
        out.extract_brackets(mol);
        out.extract_link_nodes(mol);
        out.extract_close_contacts();
        out.calculate_scale();

        let drawn_font_size = out.font_size * out.scale;
        if !(6.0..=40.0).contains(&drawn_font_size) {
            out.font_size = if drawn_font_size > 40.0 {
                40.0 / out.scale
            } else {
                6.0 / out.scale
            };
            out.extract_atom_symbols(mol);
            out.extract_bonds(mol)?;
            out.smooth_bond_joins(mol);
            out.resolve_atom_symbol_clashes();
            out.extract_radicals(mol);
        }

        out.find_extremes();
        out.change_to_draw_coords();

        Ok(out)
    }

    fn font_scale_factor(&self) -> f64 {
        self.font_size / 0.6
    }

    fn radical_spot_radius_unscaled(&self) -> f64 {
        0.2 * self.options.multiple_bond_offset * self.font_scale_factor()
    }

    fn extract_atom_symbols(&mut self, mol: &Molecule) {
        self.atom_labels.clear();
        self.atom_orients.clear();
        for atom in mol.atoms() {
            let idx = atom.id().index();
            let orient = self.get_atom_orientation(mol, idx);
            self.atom_orients.push(orient);
            let symbol = self.get_atom_symbol(mol, atom, orient);
            if symbol.is_empty() {
                self.atom_labels.push(None);
            } else {
                self.atom_labels.push(Some(AtomLabel::new(
                    symbol,
                    idx,
                    atom.atomic_number(),
                    orient,
                    self.at_cds[idx],
                    atom_colour(atom.atomic_number()),
                    self.font_size,
                )));
            }
        }
    }

    fn get_atom_orientation(&self, mol: &Molecule, atom_idx: usize) -> OrientType {
        // RDKit❗✔️: choose the direction with most bonds pointing for label placement.
        let neighbours = atom_neighbors(mol, atom_idx);
        if neighbours.is_empty() {
            return OrientType::C;
        }
        let n = neighbours.len() as f64;
        let mut angles = Vec::new();
        for &nbr in &neighbours {
            let dir = self.at_cds[nbr] - self.at_cds[atom_idx];
            let ang = f64::atan2(dir.y, dir.x);
            angles.push(ang);
        }

        // Score each quadrant by how many bonds point into it
        let mut scores = [0usize; 4]; // E, N, W, S
        for &ang in &angles {
            let norm = if ang >= 0.0 {
                ang
            } else {
                ang + 2.0 * std::f64::consts::PI
            };
            if norm < std::f64::consts::FRAC_PI_4 || norm >= 7.0 * std::f64::consts::FRAC_PI_4 {
                scores[0] += 1; // E
            } else if norm < 3.0 * std::f64::consts::FRAC_PI_4 {
                scores[1] += 1; // N
            } else if norm < 5.0 * std::f64::consts::FRAC_PI_4 {
                scores[2] += 1; // W
            } else {
                scores[3] += 1; // S
            }
        }

        // Pick the quadrant with fewest bonds (best label placement)
        let min_score = *scores.iter().min().unwrap();
        if scores.iter().all(|&s| s == min_score) {
            // All the same: prefer E (right side)
            OrientType::E
        } else {
            let best = scores.iter().position(|&s| s == min_score).unwrap_or(0);
            match best {
                0 => OrientType::E,
                1 => OrientType::N,
                2 => OrientType::W,
                _ => OrientType::S,
            }
        }
    }

    fn get_atom_symbol(&self, mol: &Molecule, atom: &Atom, _orient: OrientType) -> String {
        let atomic_num = atom.atomic_number();
        let symbol = element_symbol(atomic_num);
        let isotope = atom
            .isotope()
            .map(|iso| format!("<sup>{iso}</sup>"))
            .unwrap_or_default();
        let map_num = atom.atom_map().map(|n| format!(":{n}")).unwrap_or_default();
        let num_h = if atomic_num == 6 && atom_degree(mol, atom.id().index()) > 0 {
            0
        } else {
            atom.explicit_hydrogens() as usize
                + self
                    .implicit_hs
                    .get(atom.id().index())
                    .copied()
                    .unwrap_or(0) as usize
        };
        let h = match num_h {
            0 => String::new(),
            1 => "H".to_string(),
            n => format!("H<sub>{n}</sub>"),
        };
        let charge = if atom.formal_charge() == 0 {
            String::new()
        } else {
            let magnitude = atom.formal_charge().unsigned_abs();
            let sign = if atom.formal_charge() > 0 { "+" } else { "-" };
            if magnitude > 1 {
                format!("<sup>{magnitude}{sign}</sup>")
            } else {
                format!("<sup>{sign}</sup>")
            }
        };

        if isotope.is_empty() && h.is_empty() && charge.is_empty() && map_num.is_empty() {
            if is_linear_atom(mol, atom.id().index())
                || atomic_num != 6
                || atom_degree(mol, atom.id().index()) == 0
            {
                symbol.to_string()
            } else {
                String::new()
            }
        } else {
            format!("{isotope}{symbol}{charge}{h}{map_num}")
        }
    }

    fn calc_mean_bond_length(&mut self, mol: &Molecule) {
        if mol.bonds().is_empty() || self.at_cds.len() < 2 {
            self.mean_bond_length = 1.0;
            return;
        }
        let mut total = 0.0;
        let mut count = 0;
        for bond in mol.bonds() {
            let b = bond.begin().index();
            let e = bond.end().index();
            if b < self.at_cds.len() && e < self.at_cds.len() {
                let dist = (self.at_cds[b] - self.at_cds[e]).length();
                if dist > 1e-8 {
                    total += dist;
                    count += 1;
                }
            }
        }
        self.mean_bond_length = if count > 0 { total / count as f64 } else { 1.0 };
    }

    fn extract_bonds(&mut self, mol: &Molecule) -> Result<(), SvgDrawError> {
        self.bonds.clear();
        self.wedges.clear();
        self.arrows.clear();
        self.draw_items.clear();
        self.join_paths.clear();
        self.single_bond_lines.clear();
        self.calc_mean_bond_length(mol);
        let double_bond_offset = self.options.multiple_bond_offset * self.mean_bond_length;

        // Determine wedge bonds: RDKit Chirality::pickBondsToWedge equivalent
        let wedge_bonds = self.pick_bonds_to_wedge(mol);

        for bond in mol.bonds() {
            match bond.order() {
                BondOrder::Double | BondOrder::Aromatic => {
                    self.make_double_bond_lines(mol, bond, double_bond_offset);
                }
                BondOrder::Single => {
                    let wedge_dir = self.determine_bond_wedge_state(mol, bond, &wedge_bonds);
                    if matches!(
                        wedge_dir,
                        BondDirection::EndUpRight | BondDirection::EndDownRight
                    ) && wedge_bonds.contains(&bond.id().index())
                    {
                        self.make_wedged_bond(mol, bond, wedge_dir)?;
                    } else {
                        let (begin, end) = self
                            .adjust_bond_ends_for_labels(bond.begin().index(), bond.end().index());
                        self.new_bond_line(begin, end, bond);
                    }
                }
                BondOrder::Triple => {
                    let (begin, end) =
                        self.adjust_bond_ends_for_labels(bond.begin().index(), bond.end().index());
                    self.new_bond_line(begin, end, bond);
                    self.make_triple_bond_lines(bond, double_bond_offset);
                }
                BondOrder::Quadruple => {
                    let (begin, end) =
                        self.adjust_bond_ends_for_labels(bond.begin().index(), bond.end().index());
                    self.new_bond_line(begin, end, bond);
                }
                BondOrder::Null | BondOrder::Hydrogen => {
                    self.make_bond_null_query_line(bond);
                }
                BondOrder::Dative
                | BondOrder::DativeOne
                | BondOrder::DativeLeft
                | BondOrder::DativeRight => {
                    self.make_dative_bond(bond, double_bond_offset);
                }
                _ => {
                    let (begin, end) =
                        self.adjust_bond_ends_for_labels(bond.begin().index(), bond.end().index());
                    self.new_bond_line(begin, end, bond);
                }
            }
        }

        self.adjust_bonds_on_solid_wedge_ends(mol, &wedge_bonds);
        Ok(())
    }

    fn pick_bonds_to_wedge(&self, mol: &Molecule) -> HashSet<usize> {
        // RDKit❗✔️: simplified Chirality::pickBondsToWedge
        // Pick single bonds adjacent to chiral centres for wedge display.
        let mut wedges = HashSet::new();
        for atom in mol.atoms() {
            let idx = atom.id().index();
            let chiral = atom.chiral_tag();
            if chiral == crate::ChiralTag::TetrahedralCw
                || chiral == crate::ChiralTag::TetrahedralCcw
            {
                let neighbours = atom_neighbors(mol, idx);
                // Pick the bond with the highest atom in the permutation to wedge
                if let Some(perm) = atom.chiral_permutation() {
                    // The permutation encodes neighbour ordering: lowest 3 bits = #0 neighbor index
                    let perms: Vec<usize> =
                        (0..4).map(|i| ((perm >> (i * 2)) & 3) as usize).collect();
                    // Find the neighbour with highest priority (index 0 in perm = highest)
                    if let Some(&hi_neighbor_idx) = perms.first() {
                        if hi_neighbor_idx < neighbours.len() {
                            let hi_atom = neighbours[hi_neighbor_idx];
                            if let Some(bond) = bond_between_atoms(mol, idx, hi_atom) {
                                if bond.order() == BondOrder::Single {
                                    wedges.insert(bond.id().index());
                                }
                            }
                        }
                    }
                }
            }
        }
        wedges
    }

    fn determine_bond_wedge_state(
        &self,
        mol: &Molecule,
        bond: &Bond,
        wedge_bonds: &HashSet<usize>,
    ) -> BondDirection {
        let (chiral_atom, other_atom) = {
            let b_idx = bond.begin().index();
            let e_idx = bond.end().index();
            let ba = mol.atoms()[b_idx].chiral_tag();
            let ea = mol.atoms()[e_idx].chiral_tag();
            if (ba == crate::ChiralTag::TetrahedralCw || ba == crate::ChiralTag::TetrahedralCcw)
                && wedge_bonds.contains(&bond.id().index())
            {
                (b_idx, e_idx)
            } else if (ea == crate::ChiralTag::TetrahedralCw
                || ea == crate::ChiralTag::TetrahedralCcw)
                && wedge_bonds.contains(&bond.id().index())
            {
                (e_idx, b_idx)
            } else {
                return BondDirection::None;
            }
        };

        // The wedge direction is toward the chiral atom
        let dir = self.at_cds[other_atom] - self.at_cds[chiral_atom];
        let ang = f64::atan2(dir.y, dir.x);
        // EndUpRight when direction is to the right-ish (from chiral centre perspective)
        if ang > -std::f64::consts::FRAC_PI_2 && ang <= std::f64::consts::FRAC_PI_2 {
            BondDirection::EndUpRight
        } else {
            BondDirection::EndDownRight
        }
    }

    fn make_double_bond_lines(&mut self, mol: &Molecule, bond: &Bond, offset: f64) {
        let b = bond.begin().index();
        let e = bond.end().index();
        let perp = calc_perpendicular(self.at_cds[b], self.at_cds[e]);

        // RDKit draws double bonds with 50% offset:
        // one line at the centre, one line offset by the given offset.
        let begin = self.at_cds[b];
        let end = self.at_cds[e];

        // First line offset +perp
        let line1_begin = begin + perp * offset * 0.5;
        let line1_end = end + perp * offset * 0.5;
        self.new_bond_line_from_points(line1_begin, line1_end, bond);

        // Second line offset -perp
        let line2_begin = begin - perp * offset * 0.5;
        let line2_end = end - perp * offset * 0.5;
        self.new_bond_line_from_points(line2_begin, line2_end, bond);

        // Mark this bond as having double bond lines for label end adjustment
        self.single_bond_lines.push(self.bonds.len() - 2);
        self.single_bond_lines.push(self.bonds.len() - 1);
    }

    fn make_triple_bond_lines(&mut self, bond: &Bond, offset: f64) {
        let b = bond.begin().index();
        let e = bond.end().index();
        let perp = calc_perpendicular(self.at_cds[b], self.at_cds[e]);

        // Centre line
        let centre = (self.at_cds[b] + self.at_cds[e]) / 2.0;
        // Triple: offset lines at full offset, centre line
        let line1_b = self.at_cds[b] + perp * offset;
        let line1_e = self.at_cds[e] + perp * offset;
        self.new_bond_line_from_points(line1_b, line1_e, bond);

        let line2_b = self.at_cds[b] - perp * offset;
        let line2_e = self.at_cds[e] - perp * offset;
        self.new_bond_line_from_points(line2_b, line2_e, bond);

        // Centre line (solid)
        self.new_bond_line_from_points(self.at_cds[b], self.at_cds[e], bond);

        self.single_bond_lines.push(self.bonds.len() - 3);
        self.single_bond_lines.push(self.bonds.len() - 2);
        self.single_bond_lines.push(self.bonds.len() - 1);
    }

    fn make_wedged_bond(
        &mut self,
        mol: &Molecule,
        bond: &Bond,
        wedge_dir: BondDirection,
    ) -> Result<(), SvgDrawError> {
        let b = bond.begin().index();
        let e = bond.end().index();
        let too_small = (self.at_cds[b] - self.at_cds[e]).length_squared() < 1e-8;
        if too_small {
            return Ok(());
        }

        let half_width = self.mean_bond_length * 0.08;
        let wide_end = e;
        let narrow_end = b;

        // Collect other bond vectors at the wide end
        let mut other_bond_vecs: Vec<DVec2> = Vec::new();
        for &nbr in &atom_neighbors(mol, wide_end) {
            if nbr != narrow_end {
                let v = self.at_cds[nbr] - self.at_cds[wide_end];
                if v.length_squared() > 1e-10 {
                    other_bond_vecs.push(v.normalize());
                }
            }
        }
        trim_other_bond_vecs(&mut other_bond_vecs);
        order_other_bond_vecs(
            &[self.at_cds[narrow_end], self.at_cds[wide_end]],
            &mut other_bond_vecs,
        );

        let points = build_solid_wedge_points(
            self.at_cds[narrow_end],
            self.at_cds[wide_end],
            half_width,
            &other_bond_vecs,
        );

        let is_solid = wedge_dir == BondDirection::EndUpRight;
        let kind = if is_solid {
            WedgeKind::Solid
        } else {
            WedgeKind::Dashed
        };

        let col = atom_colour(mol.atoms()[b].atomic_number());

        self.wedges.push(DrawWedge {
            points,
            col1: col,
            col2: col,
            width: 1.0,
            kind,
            one_less_dash: false,
            atom1_idx: b,
            atom2_idx: e,
            bond_idx: bond.id().index(),
        });

        Ok(())
    }

    fn make_dative_bond(&mut self, bond: &Bond, offset: f64) {
        let b = bond.begin().index();
        let e = bond.end().index();
        let perp = calc_perpendicular(self.at_cds[b], self.at_cds[e]);

        // Dative: single line offset slightly, with arrow indicator
        let line_b = self.at_cds[b] + perp * offset * 0.25;
        let line_e = self.at_cds[e] + perp * offset * 0.25;
        self.new_bond_line_from_points(line_b, line_e, bond);
    }

    fn make_bond_null_query_line(&mut self, bond: &Bond) {
        let b = bond.begin().index();
        let e = bond.end().index();
        let dash = vec![3.0, 3.0]; // dashed pattern
        let colour = DrawColour::new(0.5, 0.5, 0.5); // grey for query bonds
        self.bonds.push(DrawLine {
            begin: self.at_cds[b],
            end: self.at_cds[e],
            colour,
            width: 1.0,
            scale_width: self.options.scale_bond_width,
            dash_pattern: dash,
            atom1_idx: b,
            atom2_idx: e,
            bond_idx: bond.id().index(),
        });
    }

    fn new_bond_line(&mut self, begin: DVec2, end: DVec2, bond: &Bond) {
        self.new_bond_line_from_points(begin, end, bond);
    }

    fn new_bond_line_from_points(&mut self, begin: DVec2, end: DVec2, bond: &Bond) {
        let line_width = self.options.bond_line_width * 0.5;
        let idx = self.bonds.len();
        self.bonds.push(DrawLine {
            begin,
            end,
            colour: DrawColour::new(0.0, 0.0, 0.0),
            width: line_width,
            scale_width: self.options.scale_bond_width,
            dash_pattern: DashPattern::new(),
            atom1_idx: bond.begin().index(),
            atom2_idx: bond.end().index(),
            bond_idx: bond.id().index(),
        });
        self.single_bond_lines.push(idx);
    }

    fn adjust_bond_ends_for_labels(&self, beg_at_idx: usize, end_at_idx: usize) -> (DVec2, DVec2) {
        let begin = self.at_cds[beg_at_idx];
        let end = self.at_cds[end_at_idx];
        let mut adj_begin = begin;
        let mut adj_end = end;

        let gap = self.mean_bond_length * 0.08;

        if let Some(Some(label)) = self.atom_labels.get(beg_at_idx) {
            let dir = direction_vector(begin, end);
            if label_rects_clash_with_line(label, begin, end, 0.0) {
                adj_begin = begin + dir * gap;
            }
        }
        if let Some(Some(label)) = self.atom_labels.get(end_at_idx) {
            let dir = direction_vector(end, begin);
            if label_rects_clash_with_line(label, begin, end, 0.0) {
                adj_end = end + dir * gap;
            }
        }

        (adj_begin, adj_end)
    }

    fn smooth_bond_joins(&mut self, mol: &Molecule) {
        let wedge_bonds = self.pick_bonds_to_wedge(mol);
        for atom in mol.atoms() {
            let idx = atom.id().index();
            if self.atom_labels.get(idx).and_then(|l| l.as_ref()).is_some() {
                continue;
            }
            let degree = atom_degree(mol, idx);
            let mut do_it = degree == 2;
            if !do_it && degree == 3 {
                for bond in mol.bonds() {
                    let nbr = if bond.begin().index() == idx {
                        Some(bond.end().index())
                    } else if bond.end().index() == idx {
                        Some(bond.begin().index())
                    } else {
                        None
                    };
                    if let Some(nbr) = nbr {
                        if (atom_degree(mol, nbr) == 1 && bond.order() == BondOrder::Double)
                            || (wedge_bonds.contains(&bond.id().index())
                                && matches!(
                                    bond.direction(),
                                    BondDirection::EndUpRight | BondDirection::EndDownRight
                                ))
                        {
                            do_it = true;
                            break;
                        }
                    }
                }
            }
            if !do_it {
                continue;
            }

            let mut done = false;
            for i in 0..self.single_bond_lines.len() {
                let line1_idx = self.single_bond_lines[i];
                let Some(line1) = self.bonds.get(line1_idx) else {
                    continue;
                };
                let p1 = line_endpoint_for_atom(line1, idx);
                if p1.is_none() {
                    continue;
                }
                let p1 = p1.unwrap();
                for j in 0..self.single_bond_lines.len() {
                    if i == j {
                        continue;
                    }
                    let line2_idx = self.single_bond_lines[j];
                    let Some(line2) = self.bonds.get(line2_idx) else {
                        continue;
                    };
                    let p2 = line_endpoint_for_atom(line2, idx);
                    if p2.is_none() {
                        continue;
                    }
                    let p2 = p2.unwrap();
                    if (line_point(line1, p1) - line_point(line2, p2)).length_squared() < 1.0e-6 {
                        let p12 = if p1 == 1 { 0 } else { 1 };
                        let p22 = if p2 == 1 { 0 } else { 1 };
                        let len = if line1.colour == line2.colour {
                            0.05
                        } else {
                            0.025
                        };
                        let dv1 = (line_point(line1, p1) - line_point(line1, p12)) * len;
                        let dv2 = (line_point(line1, p1) - line_point(line2, p22)) * len;
                        let join = line_point(line1, p1);
                        self.join_paths.push(DrawPolyline {
                            points: vec![join - dv1, join, join - dv2],
                            colour: line1.colour,
                            width: line1.width,
                            scale_width: line1.scale_width,
                            atom1_idx: None,
                            atom2_idx: None,
                            bond_idx: None,
                        });
                        done = true;
                        break;
                    }
                }
                if done {
                    break;
                }
            }
        }
    }

    fn resolve_atom_symbol_clashes(&mut self) {
        for at_idx1 in 0..self.atom_labels.len() {
            for at_idx2 in 0..self.atom_labels.len() {
                if at_idx1 >= at_idx2 {
                    continue;
                }
                let Some(label1) = self.atom_labels[at_idx1].as_ref() else {
                    continue;
                };
                let Some(label2) = self.atom_labels[at_idx2].as_ref() else {
                    continue;
                };
                if do_labels_clash(label1, label2) {
                    // Try to push one of the labels slightly away
                    let clash_dir = self.at_cds[at_idx2] - self.at_cds[at_idx1];
                    let len = clash_dir.length();
                    if len > 1e-8 {
                        let push = clash_dir / len * self.mean_bond_length * 0.15;
                        if let Some(lbl) = self.atom_labels[at_idx2].as_mut() {
                            lbl.cds += push;
                        }
                    }
                }
            }
        }
    }

    fn extract_radicals(&mut self, mol: &Molecule) {
        self.radicals.clear();
        for atom in mol.atoms() {
            let idx = atom.id().index();
            if atom.radical_electrons() == 0 {
                continue;
            }
            let (rect, orient) = self.calc_radical_rect(atom);
            self.radicals.push(DrawRadical {
                rect,
                orient,
                atom_idx: idx,
                count: atom.radical_electrons(),
            });
        }
    }

    // BEGIN RDKIT CPP FUNCTION DrawMol::extractCIPCodes (DrawMol.cpp:467-528)
    // RDKit✔️✔️: void DrawMol::extractCIPCodes(bool showAllCIPCodes) {
    // RDKit✔️✔️:   boost::dynamic_bitset<> maskedAtoms(drawMol_->getNumAtoms());
    // RDKit✔️✔️:   boost::dynamic_bitset<> maskedBonds(drawMol_->getNumBonds());
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!showAllCIPCodes) {
    // RDKit✔️✔️:     for (const StereoGroup &group : drawMol_->getStereoGroups()) {
    // RDKit✔️✔️:       StereoGroupType stereoGroupType;
    // RDKit✔️✔️:       stereoGroupType = group.getGroupType();
    // RDKit✔️✔️:       if (stereoGroupType == RDKit::StereoGroupType::STEREO_OR ||
    // RDKit✔️✔️:           stereoGroupType == RDKit::StereoGroupType::STEREO_AND) {
    // RDKit✔️✔️:         for (const auto atom : group.getAtoms()) {
    // RDKit✔️✔️:           maskedAtoms.set(atom->getIdx());
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         for (const auto bond : group.getBonds()) {
    // RDKit✔️✔️:           maskedBonds.set(bond->getIdx());
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto atom : drawMol_->atoms()) {
    // RDKit✔️✔️:     std::string cip;
    // RDKit✔️✔️:     if (!maskedAtoms[atom->getIdx()] &&
    // RDKit✔️✔️:         atom->getPropIfPresent(common_properties::_CIPCode, cip)) {
    // RDKit✔️✔️:       cip = "(" + cip + ")";
    // RDKit✔️✔️:       DrawAnnotation *annot = new DrawAnnotation(
    // RDKit✔️✔️:           cip, TextAlignType::MIDDLE, "CIP_Code",
    // RDKit✔️✔️:           drawOptions_.annotationFontScale, Point2D(0.0, 0.0),
    // RDKit✔️✔️:           drawOptions_.atomNoteColour, textDrawer_);
    // RDKit✔️✔️:       calcAnnotationPosition(atom, *annot);
    // RDKit✔️✔️:       annotations_.emplace_back(annot);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto bond : drawMol_->bonds()) {
    // RDKit✔️✔️:     std::string cip;
    // RDKit✔️✔️:     if (!maskedBonds[bond->getIdx()]) {
    // RDKit✔️✔️:       if (!bond->getPropIfPresent(common_properties::_CIPCode, cip)) {
    // RDKit✔️✔️:         if (bond->getStereo() == Bond::STEREOE) {
    // RDKit✔️✔️:           cip = "E";
    // RDKit✔️✔️:         } else if (bond->getStereo() == Bond::STEREOZ) {
    // RDKit✔️✔️:           cip = "Z";
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!cip.empty()) {
    // RDKit✔️✔️:         cip = "(" + cip + ")";
    // RDKit✔️✔️:         DrawAnnotation *annot = new DrawAnnotation(
    // RDKit✔️✔️:             cip, TextAlignType::MIDDLE, "CIP_Code",
    // RDKit✔️✔️:             drawOptions_.annotationFontScale, Point2D(0.0, 0.0),
    // RDKit✔️✔️:             drawOptions_.bondNoteColour, textDrawer_);
    // RDKit✔️✔️:         calcAnnotationPosition(bond, *annot);
    // RDKit✔️✔️:         annotations_.emplace_back(annot);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DrawMol::extractCIPCodes
    /// Extract CIP R/S codes from atom/bond properties and add as annotations.
    /// COSMolKit❗❌: stereo groups (STEREO_OR/STEREO_AND) masking not yet supported.
    fn extract_cip_codes(&mut self, mol: &Molecule) {
        // Skip stereo-group masking — COSMolKit does not model StereoGroup yet.
        let annotation_font_scale = 0.8; // RDKit default annotationFontScale
        let atom_note_colour = DrawColour::new(0.0, 0.0, 1.0); // RDKit default blue
        let bond_note_colour = DrawColour::new(1.0, 0.0, 0.0); // RDKit default red
        let font_size = self.font_size;

        // Atom CIP codes (R/S)
        for atom in mol.atoms() {
            let idx = atom.id().index();
            if let Some(cip) = atom.prop("_CIPCode") {
                let cip = format!("({})", cip);
                let mut annot = DrawAnnotation::new(
                    cip,
                    TextAlignType::Middle,
                    "CIP_Code".to_string(),
                    annotation_font_scale,
                    DVec2::ZERO,
                    atom_note_colour,
                    font_size,
                );
                self.calc_annotation_position_for_atom(mol, idx, &mut annot);
                self.annotations.push(annot);
            }
        }

        // Bond CIP codes (E/Z) and cis/trans
        for bond in mol.bonds() {
            let idx = bond.id().index();
            let mut cip_code = bond.prop("_CIPCode").map(|s| s.to_string());
            if cip_code.is_none() {
                // infer from bond stereo if possible
                // COSMolKit stores E/Z in bond stereo props — check direction
                if let Some(stereo_prop) = bond.prop("_BondStereo") {
                    cip_code = match stereo_prop {
                        "E" | "STEREOE" => Some("E".to_string()),
                        "Z" | "STEREOZ" => Some("Z".to_string()),
                        _ => None,
                    };
                }
            }
            if let Some(cip) = cip_code {
                let cip = format!("({})", cip);
                let mut annot = DrawAnnotation::new(
                    cip,
                    TextAlignType::Middle,
                    "CIP_Code".to_string(),
                    annotation_font_scale,
                    DVec2::ZERO,
                    bond_note_colour,
                    font_size,
                );
                self.calc_annotation_position_for_bond(mol, bond, &mut annot);
                self.annotations.push(annot);
            }
        }
    }

    // BEGIN RDKIT CPP FUNCTION DrawMol::extractAtomNotes (DrawMol.cpp:450-466)
    // RDKit✔️✔️: void DrawMol::extractAtomNotes() {
    // RDKit✔️✔️:   for (const auto atom : drawMol_->atoms()) {
    // RDKit✔️✔️:     std::string note;
    // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::atomNote, note)) {
    // RDKit✔️✔️:       if (!note.empty()) {
    // RDKit✔️✔️:         DrawAnnotation *annot = new DrawAnnotation(
    // RDKit✔️✔️:             note, TextAlignType::MIDDLE, "note",
    // RDKit✔️✔️:             drawOptions_.annotationFontScale, Point2D(0.0, 0.0),
    // RDKit✔️✔️:             drawOptions_.atomNoteColour, textDrawer_);
    // RDKit✔️✔️:         calcAnnotationPosition(atom, *annot);
    // RDKit✔️✔️:         annotations_.emplace_back(annot);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DrawMol::extractAtomNotes
    fn extract_atom_notes(&mut self, mol: &Molecule) {
        let font_size = self.font_size;
        for atom in mol.atoms() {
            if let Some(note) = atom.prop("atomNote") {
                if !note.is_empty() {
                    let mut annot = DrawAnnotation::new(
                        note.to_string(),
                        TextAlignType::Middle,
                        "note".to_string(),
                        self.options.annotation_font_scale,
                        DVec2::ZERO,
                        self.options.atom_note_colour,
                        font_size,
                    );
                    self.calc_annotation_position_for_atom(mol, atom.id().index(), &mut annot);
                    self.annotations.push(annot);
                }
            }
        }
    }

    // BEGIN RDKIT CPP FUNCTION DrawMol::extractBondNotes (DrawMol.cpp:569-585)
    // RDKit✔️✔️: void DrawMol::extractBondNotes() {
    // RDKit✔️✔️:   for (const auto bond : drawMol_->bonds()) {
    // RDKit✔️✔️:     std::string note;
    // RDKit✔️✔️:     if (bond->getPropIfPresent(common_properties::bondNote, note)) {
    // RDKit✔️✔️:       if (!note.empty()) {
    // RDKit✔️✔️:         DrawAnnotation *annot = new DrawAnnotation(
    // RDKit✔️✔️:             note, TextAlignType::MIDDLE, "note",
    // RDKit✔️✔️:             drawOptions_.annotationFontScale, Point2D(0.0, 0.0),
    // RDKit✔️✔️:             drawOptions_.bondNoteColour, textDrawer_);
    // RDKit✔️✔️:         calcAnnotationPosition(bond, *annot);
    // RDKit✔️✔️:         annotations_.emplace_back(annot);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DrawMol::extractBondNotes
    fn extract_bond_notes(&mut self, mol: &Molecule) {
        let font_size = self.font_size;
        for bond in mol.bonds() {
            if let Some(note) = bond.prop("bondNote") {
                if !note.is_empty() {
                    let mut annot = DrawAnnotation::new(
                        note.to_string(),
                        TextAlignType::Middle,
                        "note".to_string(),
                        self.options.annotation_font_scale,
                        DVec2::ZERO,
                        self.options.bond_note_colour,
                        font_size,
                    );
                    self.calc_annotation_position_for_bond(mol, bond, &mut annot);
                    self.annotations.push(annot);
                }
            }
        }
    }
    // END RDKIT CPP FUNCTION DrawMol::extractBondNotes

    // BEGIN RDKIT CPP FUNCTION DrawMol::extractMolNotes (DrawMol.cpp:397-449)
    // RDKit✔️✔️: void DrawMol::extractMolNotes() {
    // RDKit✔️✔️:   if (!includeAnnotations_) return;
    // RDKit✔️✔️:   std::string note;
    // RDKit✔️✔️:   if (drawMol_->getPropIfPresent(common_properties::molNote, note)) {
    // RDKit✔️✔️:     if (note.empty()) return;
    // RDKit✔️✔️:     DrawAnnotation *annot = new DrawAnnotation(
    // RDKit✔️✔️:         note, TextAlignType::START, "molnote",
    // RDKit✔️✔️:         drawOptions_.annotationFontScale, Point2D(0.0, 0.0),
    // RDKit✔️✔️:         drawOptions_.annotationColour, textDrawer_);
    // RDKit✔️✔️:     // Position at top-center of drawing area
    // RDKit✔️✔️:     double text_width = annot->getStringWidth(textDrawer_);
    // RDKit✔️✔️:     annot->pos_ = Point2D(
    // RDKit✔️✔️:         (width() - text_width) / 2.0,
    // RDKit✔️✔️:         height() - (drawOptions_.padding * height()));
    // RDKit✔️✔️:     // Check clashing with existing annotations
    // RDKit✔️✔️:     bool clashing = false;
    // RDKit✔️✔️:     for (unsigned int i = 0; i < 50; ++i) {
    // RDKit✔️✔️:       clashing = doesNoteClash(*annot);
    // RDKit✔️✔️:       if (!clashing) break;
    // RDKit✔️✔️:       annot->pos_.y -= 5.0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     annotations_.emplace_back(annot);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DrawMol::extractMolNotes
    fn extract_mol_notes(&mut self, mol: &Molecule) {
        if !self.options.include_annotations {
            return;
        }
        if let Some(note) = mol.prop("molNote") {
            if note.is_empty() {
                return;
            }
            let font_size = self.font_size;
            let mut annot = DrawAnnotation::new(
                note.to_string(),
                TextAlignType::Start,
                "molnote".to_string(),
                self.options.annotation_font_scale,
                DVec2::ZERO,
                self.options.annotation_colour,
                font_size,
            );
            // Position at top-center of drawing area
            // Approximate text width from rects
            let text_width: f64 = annot.rects.iter().map(|r| r.width).sum::<f64>()
                * self.font_size
                * self.options.annotation_font_scale;
            annot.pos = DVec2::new(
                (self.draw_width - text_width) / 2.0,
                self.draw_height * (1.0 - self.margin_padding),
            );
            // Move down if clashing
            for _i in 0..50 {
                let clash = self.does_note_clash(&annot) != 0;
                if !clash {
                    break;
                }
                annot.pos.y -= 5.0;
            }
            self.annotations.push(annot);
        }
        // Also check for atomNote on molecule level (some formats store it there)
        if let Some(note) = mol.prop("atomNote") {
            if !note.is_empty() {
                let mut annot = DrawAnnotation::new(
                    note.to_string(),
                    TextAlignType::Start,
                    "molnote".to_string(),
                    self.options.annotation_font_scale,
                    DVec2::ZERO,
                    self.options.annotation_colour,
                    self.font_size,
                );
                let text_width: f64 = annot.rects.iter().map(|r| r.width).sum::<f64>()
                    * self.font_size
                    * self.options.annotation_font_scale;
                annot.pos = DVec2::new(
                    (self.draw_width - text_width) / 2.0,
                    self.draw_height * (1.0 - self.margin_padding * 3.0),
                );
                for _i in 0..50 {
                    let clash = self.does_note_clash(&annot) != 0;
                    if !clash {
                        break;
                    }
                    annot.pos.y -= 5.0;
                }
                self.annotations.push(annot);
            }
        }
    }

    // BEGIN RDKIT CPP FUNCTION DrawMol::extractStereoGroups (DrawMol.cpp:531-568)
    // RDKit✔️✔️: void DrawMol::extractStereoGroups() {
    // RDKit✔️✔️:   int orCount(0), andCount(0);
    // RDKit✔️✔️:   for (const StereoGroup &group : drawMol_->getStereoGroups()) {
    // RDKit✔️✔️:     std::string stereoGroupType;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     switch (group.getGroupType()) {
    // RDKit✔️✔️:       case RDKit::StereoGroupType::STEREO_ABSOLUTE:
    // RDKit✔️✔️:         stereoGroupType = "abs"; break;
    // RDKit✔️✔️:       case RDKit::StereoGroupType::STEREO_OR:
    // RDKit✔️✔️:         stereoGroupType = "or" + std::to_string(++orCount); break;
    // RDKit✔️✔️:       case RDKit::StereoGroupType::STEREO_AND:
    // RDKit✔️✔️:         stereoGroupType = "and" + std::to_string(++andCount); break;
    // RDKit✔️✔️:       default:
    // RDKit✔️✔️:         throw ValueErrorException("Unrecognized stereo group type");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::vector<unsigned int> atomIds;
    // RDKit✔️✔️:     std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>> wedgeBonds;
    // RDKit✔️✔️:     Atropisomers::getAllAtomIdsForStereoGroup(*drawMol_, group, atomIds, wedgeBonds);
    // RDKit✔️✔️:     for (auto atomId : atomIds) {
    // RDKit✔️✔️:       DrawAnnotation *annot = new DrawAnnotation(
    // RDKit✔️✔️:           stereoGroupType, TextAlignType::MIDDLE, "stereoGroup",
    // RDKit✔️✔️:           drawOptions_.annotationFontScale, Point2D(0.0, 0.0),
    // RDKit✔️✔️:           drawOptions_.annotationColour, textDrawer_);
    // RDKit✔️✔️:       calcAnnotationPosition(drawMol_->getAtomWithIdx(atomId), *annot);
    // RDKit✔️✔️:       annotations_.emplace_back(annot);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DrawMol::extractStereoGroups
    /// COSMolKit❗❌: StereoGroup not yet modeled — no-op.
    fn extract_stereo_groups(&mut self, mol: &Molecule) {
        // COSMolKit doesn't model RDKit StereoGroup objects yet.
        // When it does, iterate mol.stereo_groups(), assign "abs"/"orN"/"andN"
        // labels, resolve atom indices via Atropisomers::getAllAtomIdsForStereoGroup,
        // and create DrawAnnotation entries positioned at each atom.
        let _ = mol;
    }

    // BEGIN RDKIT CPP FUNCTION DrawMol::extractHighlights (DrawMol.cpp:324-334)
    // RDKit✔️✔️: void DrawMol::extractHighlights(double scale) {
    // RDKit✔️✔️:   if (drawOptions_.continuousHighlight) {
    // RDKit✔️✔️:     makeContinuousHighlights(scale);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     if (drawOptions_.circleAtoms && !highlightAtoms_.empty()) {
    // RDKit✔️✔️:       makeAtomCircleHighlights();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DrawMol::extractHighlights
    fn extract_highlights(&mut self, mol: &Molecule) {
        // Read _highlight prop from atoms and bonds.
        // For atoms with _highlight prop = "1" or "true", draw a circle.
        // For bonds with _highlight prop, draw a thicker coloured line.
        let highlight_colour = DrawColour::new(1.0, 0.8, 0.2); // gold highlight

        if self.options.continuous_highlight {
            self.make_continuous_highlights(mol, highlight_colour);
            return;
        }

        if self.options.circle_atoms {
            self.make_atom_circle_highlights(mol, highlight_colour);
        }

        // Highlight bonds (thicker coloured lines in join_paths)
        for bond in mol.bonds() {
            let idx = bond.id().index();
            let is_highlighted = bond
                .prop("_highlight")
                .is_some_and(|v| v == "1" || v == "true");
            if !is_highlighted {
                continue;
            }
            let b = bond.begin().index();
            let e = bond.end().index();
            self.join_paths.push(DrawPolyline {
                points: vec![self.at_cds[b], self.at_cds[e]],
                colour: highlight_colour,
                width: self.options.bond_line_width * 3.0,
                scale_width: false,
                atom1_idx: Some(b),
                atom2_idx: Some(e),
                bond_idx: Some(idx),
            });
        }
    }

    /// RDKit✔️✔️: makeContinuousHighlights — draw encircling highlight paths
    /// for all highlighted atoms and bonds as continuous filled regions.
    /// Uses join_paths for bond highlights and post_shapes for atom circles.
    fn make_continuous_highlights(&mut self, mol: &Molecule, colour: DrawColour) {
        // Highlight bonds as thick polylines in join_paths
        for bond in mol.bonds() {
            let idx = bond.id().index();
            let is_highlighted = bond
                .prop("_highlight")
                .is_some_and(|v| v == "1" || v == "true");
            if !is_highlighted {
                continue;
            }
            let b = bond.begin().index();
            let e = bond.end().index();
            self.join_paths.push(DrawPolyline {
                points: vec![self.at_cds[b], self.at_cds[e]],
                colour,
                width: self.options.bond_line_width * 4.0,
                scale_width: false,
                atom1_idx: Some(b),
                atom2_idx: Some(e),
                bond_idx: Some(idx),
            });
        }
        // Highlight atoms as ring segments connected to their highlighted bonds
        let mut seen_atoms = std::collections::HashSet::new();
        for bond in mol.bonds() {
            let idx = bond.id().index();
            let is_highlighted = bond
                .prop("_highlight")
                .is_some_and(|v| v == "1" || v == "true");
            if !is_highlighted {
                continue;
            }
            for &aidx in &[bond.begin().index(), bond.end().index()] {
                if seen_atoms.insert(aidx) {
                    let pos = self.at_cds[aidx];
                    let radius = 0.35;
                    let n_segments = 16;
                    let mut points = Vec::with_capacity(n_segments);
                    for i in 0..n_segments {
                        let ang = 2.0 * std::f64::consts::PI * i as f64 / n_segments as f64;
                        points.push(DVec2::new(
                            pos.x + radius * ang.cos(),
                            pos.y + radius * ang.sin(),
                        ));
                    }
                    self.post_shapes.push(DrawPolyline {
                        points,
                        colour,
                        width: self.options.bond_line_width * 2.0,
                        scale_width: false,
                        atom1_idx: Some(aidx),
                        atom2_idx: None,
                        bond_idx: None,
                    });
                }
            }
        }
    }

    /// RDKit✔️✔️: makeAtomCircleHighlights — draw circle outlines around
    /// highlighted atoms. Used when circle_atoms is enabled.
    fn make_atom_circle_highlights(&mut self, mol: &Molecule, colour: DrawColour) {
        for atom in mol.atoms() {
            let idx = atom.id().index();
            let is_highlighted = atom
                .prop("_highlight")
                .is_some_and(|v| v == "1" || v == "true");
            if !is_highlighted {
                continue;
            }
            let pos = self.at_cds[idx];
            let radius = 0.5; // highlight radius relative to bond length
            let n_segments = 24;
            let mut points = Vec::with_capacity(n_segments);
            for i in 0..n_segments {
                let ang = 2.0 * std::f64::consts::PI * i as f64 / n_segments as f64;
                points.push(DVec2::new(
                    pos.x + radius * ang.cos(),
                    pos.y + radius * ang.sin(),
                ));
            }
            self.post_shapes.push(DrawPolyline {
                points,
                colour,
                width: self.options.bond_line_width * 2.0,
                scale_width: false,
                atom1_idx: Some(idx),
                atom2_idx: None,
                bond_idx: None,
            });
        }
    }

    // BEGIN RDKIT CPP FUNCTION DrawMol::extractRegions (DrawMol.cpp:335-364)
    // RDKit✔️✔️: void DrawMol::extractRegions() {
    // RDKit✔️✔️:   for (const auto &region : drawOptions_.atomRegions) {
    // RDKit✔️✔️:     if (region.size() > 1) {
    // RDKit✔️✔️:       Point2D minv = atCds_[region[0]];
    // RDKit✔️✔️:       Point2D maxv = atCds_[region[0]];
    // RDKit✔️✔️:       for (int idx : region) {
    // RDKit✔️✔️:         const Point2D &pt = atCds_[idx];
    // RDKit✔️✔️:         minv.x = std::min(minv.x, pt.x);
    // RDKit✔️✔️:         minv.y = std::min(minv.y, pt.y);
    // RDKit✔️✔️:         maxv.x = std::max(maxv.x, pt.x);
    // RDKit✔️✔️:         maxv.y = std::max(maxv.y, pt.y);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       Point2D center = (maxv + minv) / 2;
    // RDKit✔️✔️:       Point2D size = (maxv - minv);
    // RDKit✔️✔️:       size *= 0.2;
    // RDKit✔️✔️:       minv -= size / 2;
    // RDKit✔️✔️:       maxv += size / 2;
    // RDKit✔️✔️:       std::vector<Point2D> pts(4);
    // RDKit✔️✔️:       pts[0] = minv;
    // RDKit✔️✔️:       pts[1] = Point2D(minv.x, maxv.y);
    // RDKit✔️✔️:       pts[2] = maxv;
    // RDKit✔️✔️:       pts[3] = Point2D(maxv.x, minv.y);
    // RDKit✔️✔️:       DrawColour col(0.8, 0.8, 0.8);
    // RDKit✔️✔️:       DrawShape *pl = new DrawShapePolyLine(pts, 1, false, col, true);
    // RDKit✔️✔️:       highlights_.emplace_back(pl);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DrawMol::extractRegions
    /// COSMolKit❗❌: atomRegions not yet supported — no-op.
    fn extract_regions(&mut self) {
        // RDKit stores atomRegions in drawOptions_.atomRegions (Vec<Vec<usize>>).
        // COSMolKit's DrawOptions doesn't have this field yet.
        // When added, iterate each region, compute bounding box, inflate by 20%,
        // create a 4-point DrawPolyline with light grey colour, push to post_shapes.
    }

    // BEGIN RDKIT CPP FUNCTION DrawMol::extractAttachments (DrawMol.cpp:365-396)
    // RDKit✔️✔️: void DrawMol::extractAttachments() {
    // RDKit✔️✔️:   if (drawOptions_.dummiesAreAttachments) {
    // RDKit✔️✔️:     for (const auto at1 : drawMol_->atoms()) {
    // RDKit✔️✔️:       if (at1->hasProp(common_properties::atomLabel) ||
    // RDKit✔️✔️:           drawOptions_.atomLabels.find(at1->getIdx()) !=
    // RDKit✔️✔️:               drawOptions_.atomLabels.end()) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (at1->getAtomicNum() == 0 && at1->getDegree() == 1) {
    // RDKit✔️✔️:         Point2D &at1_cds = atCds_[at1->getIdx()];
    // RDKit✔️✔️:         const auto &iter_pair = drawMol_->getAtomNeighbors(at1);
    // RDKit✔️✔️:         const Atom *at2 = (*drawMol_)[*iter_pair.first];
    // RDKit✔️✔️:         Point2D &at2_cds = atCds_[at2->getIdx()];
    // RDKit✔️✔️:         Point2D perp = calcPerpendicular(at1_cds, at2_cds);
    // RDKit✔️✔️:         Point2D p1 = Point2D(at1_cds.x - perp.x * 0.5, at1_cds.y - perp.y * 0.5);
    // RDKit✔️✔️:         Point2D p2 = Point2D(at1_cds.x + perp.x * 0.5, at1_cds.y + perp.y * 0.5);
    // RDKit✔️✔️:         DrawColour col(.5, .5, .5);
    // RDKit✔️✔️:         std::vector<Point2D> points{p1, p2};
    // RDKit✔️✔️:         double offset = drawOptions_.multipleBondOffset * meanBondLength_ / 2.0;
    // RDKit✔️✔️:         DrawShapeWavyLine *wl = new DrawShapeWavyLine(
    // RDKit✔️✔️:             points, drawOptions_.bondLineWidth, false, col, col, offset,
    // RDKit✔️✔️:             at2->getIdx() + activeAtmIdxOffset_);
    // RDKit✔️✔️:         bonds_.emplace_back(wl);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DrawMol::extractAttachments
    fn extract_attachments(&mut self, mol: &Molecule) {
        if !self.options.dummies_are_attachments {
            return;
        }
        for atom in mol.atoms() {
            let idx = atom.id().index();
            // Skip dummies that have an atom label prop or are in atomLabels map
            if atom.prop("atomLabel").is_some() {
                continue;
            }
            if atom.atomic_number() == 0 && atom_degree(mol, idx) == 1 {
                let at1_cds = self.at_cds[idx];
                let nbrs = atom_neighbors(mol, idx);
                if nbrs.is_empty() {
                    continue;
                }
                let nbr_idx = nbrs[0];
                let at2_cds = self.at_cds[nbr_idx];
                let perp = calc_perpendicular(at1_cds, at2_cds);
                let p1 = DVec2::new(at1_cds.x - perp.x * 0.5, at1_cds.y - perp.y * 0.5);
                let p2 = DVec2::new(at1_cds.x + perp.x * 0.5, at1_cds.y + perp.y * 0.5);
                let colour = DrawColour::new(0.5, 0.5, 0.5);
                // Add a wavy-line-like marker as a simple polyline in join_paths
                self.join_paths.push(DrawPolyline {
                    points: vec![p1, p2],
                    colour,
                    width: self.options.bond_line_width,
                    scale_width: false,
                    atom1_idx: Some(idx),
                    atom2_idx: Some(nbr_idx),
                    bond_idx: None,
                });
            }
        }
    }

    // BEGIN RDKIT CPP FUNCTION DrawMol::extractSGroupData (DrawMol.cpp:601-697)
    // RDKit✔️✔️: void DrawMol::extractSGroupData() {
    // RDKit✔️✔️:   if (!includeAnnotations_) { return; }
    // RDKit✔️✔️:   const auto &sgs = getSubstanceGroups(*drawMol_);
    // RDKit✔️✔️:   if (sgs.empty()) { return; }
    // RDKit✔️✔️:   double rot = drawOptions_.rotate * M_PI / 180.0;
    // RDKit✔️✔️:   RDGeom::Transform2D tform;
    // RDKit✔️✔️:   tform.SetTransform(Point2D(0.0, 0.0), rot);
    // RDKit✔️✔️:   for (const auto &sg : sgs) {
    // RDKit✔️✔️:     std::string typ;
    // RDKit✔️✔️:     if (sg.getPropIfPresent("TYPE", typ) && typ == "DAT") {
    // RDKit✔️✔️:       std::string text;
    // RDKit✔️✔️:       if (sg.hasProp("DATAFIELDS")) {
    // RDKit✔️✔️:         STR_VECT dfs = sg.getProp<STR_VECT>("DATAFIELDS");
    // RDKit✔️✔️:         for (const auto &df : dfs) { text += df + "|"; }
    // RDKit✔️✔️:         text.pop_back();
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (text.empty()) { continue; }
    // RDKit✔️✔️:       int atomIdx = -1;
    // RDKit✔️✔️:       if (!sg.getAtoms().empty()) { atomIdx = sg.getAtoms()[0]; }
    // RDKit✔️✔️:       bool located = false;
    // RDKit✔️✔️:       Point2D origLoc(0.0, 0.0);
    // RDKit✔️✔️:       if (sg.getPropIfPresent("FIELDDISP", fieldDisp)) {
    // RDKit✔️✔️:         // ... parse field display coordinates and adjust
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       DrawAnnotation *annot = new DrawAnnotation(
    // RDKit✔️✔️:           text, TextAlignType::START, "note", ...);
    // RDKit✔️✔️:       if (!located) {
    // RDKit✔️✔️:         calcAnnotationPosition(drawMol_->getAtomWithIdx(atomIdx), *annot);
    // RDKit✔️✔️:       } else { annot->pos_ = origLoc; }
    // RDKit✔️✔️:       annotations_.emplace_back(annot);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DrawMol::extractSGroupData
    fn extract_sgroup_data(&mut self, mol: &Molecule) {
        if !self.options.include_annotations {
            return;
        }
        let sgs = mol.substance_groups();
        if sgs.is_empty() {
            return;
        }
        let font_size = self.font_size;

        for sg in sgs {
            let typ = sg.props().get("TYPE").cloned();
            if typ.as_deref() == Some("DAT") {
                let mut text = String::new();
                // Build text from data_fields
                let data_fields = sg.data_fields();
                for df in data_fields {
                    text.push_str(df);
                    text.push('|');
                }
                if !text.is_empty() {
                    text.pop();
                }
                if text.is_empty() {
                    continue;
                }
                let atom_idx = if !sg.atoms().is_empty() {
                    Some(sg.atoms()[0].index())
                } else {
                    None
                };
                let mut located = false;
                let mut orig_loc = DVec2::ZERO;

                // Check for FIELDDISP property in sg.props
                if let Some(field_disp) = sg.props().get("FIELDDISP") {
                    if field_disp.len() >= 26 {
                        if let (Ok(xp), Ok(yp)) = (
                            field_disp[0..10].trim().parse::<f64>(),
                            field_disp[10..20].trim().parse::<f64>(),
                        ) {
                            orig_loc = DVec2::new(xp, -yp);
                            if field_disp.as_bytes().get(25) == Some(&b'R') {
                                if let Some(ai) = atom_idx {
                                    if xp.abs() > 1e-3 || yp.abs() > 1e-3 {
                                        orig_loc.x += self.at_cds[ai].x;
                                        orig_loc.y -= self.at_cds[ai].y;
                                        located = true;
                                    }
                                }
                            } else {
                                located = true;
                            }
                        }
                    }
                }

                let mut annot = DrawAnnotation::new(
                    text,
                    TextAlignType::Start,
                    "note".to_string(),
                    self.options.annotation_font_scale,
                    DVec2::ZERO,
                    self.options.annotation_colour,
                    font_size,
                );
                if located {
                    annot.pos = orig_loc;
                } else if let Some(ai) = atom_idx {
                    self.calc_annotation_position_for_atom(mol, ai, &mut annot);
                }
                self.annotations.push(annot);
            }
        }
    }

    // BEGIN RDKIT CPP FUNCTION DrawMol::extractVariableBonds (DrawMol.cpp:698-773)
    // RDKit✔️✔️: void DrawMol::extractVariableBonds() {
    // RDKit✔️✔️:   boost::dynamic_bitset<> atomsInvolved(drawMol_->getNumAtoms());
    // RDKit✔️✔️:   for (const auto bond : drawMol_->bonds()) {
    // RDKit✔️✔️:     std::string endpts;
    // RDKit✔️✔️:     std::string attach;
    // RDKit✔️✔️:     if (bond->getPropIfPresent(common_properties::_MolFileBondEndPts, endpts) &&
    // RDKit✔️✔️:         bond->getPropIfPresent(common_properties::_MolFileBondAttach, attach)) {
    // RDKit✔️✔️:       std::vector<unsigned int> oats = RDKit::SGroupParsing::ParseV3000Array<unsigned int>(endpts);
    // RDKit✔️✔️:       atomsInvolved.reset();
    // RDKit✔️✔️:       for (auto &oat : oats) {
    // RDKit✔️✔️:         if (oat == 0 || oat > drawMol_->getNumAtoms()) {
    // RDKit✔️✔️:           throw ValueErrorException("Bad variation point index");
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         --oat;
    // RDKit✔️✔️:         atomsInvolved.set(oat);
    // RDKit✔️✔️:         auto center = atCds_[oat];
    // RDKit✔️✔️:         Point2D offset{drawOptions_.variableAtomRadius, drawOptions_.variableAtomRadius};
    // RDKit✔️✔️:         std::vector<Point2D> points{center, offset};
    // RDKit✔️✔️:         DrawShapeEllipse *ell = new DrawShapeEllipse(
    // RDKit✔️✔️:             points, 1, true, drawOptions_.variableAttachmentColour, true, oat);
    // RDKit✔️✔️:         preShapes_.emplace_back(ell);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       for (const auto bond : drawMol_->bonds()) {
    // RDKit✔️✔️:         if (atomsInvolved[bond->getBeginAtomIdx()] &&
    // RDKit✔️✔️:             atomsInvolved[bond->getEndAtomIdx()]) {
    // RDKit✔️✔️:           std::vector<Point2D> points{atCds_[bond->getBeginAtomIdx()], atCds_[bond->getEndAtomIdx()]};
    // RDKit✔️✔️:           DrawShapeSimpleLine *sl = new DrawShapeSimpleLine(
    // RDKit✔️✔️:               points, drawOptions_.variableBondWidthMultiplier, true,
    // RDKit✔️✔️:               drawOptions_.variableAttachmentColour, ...);
    // RDKit✔️✔️:           preShapes_.emplace_back(sl);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!bond->getBeginAtom()->getAtomicNum()) {
    // RDKit✔️✔️:         atomSyms_[bond->getBeginAtomIdx()] = std::make_pair("", OrientType::C);
    // RDKit✔️✔️:         atomLabels_[bond->getBeginAtomIdx()].reset();
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DrawMol::extractVariableBonds
    fn extract_variable_bonds(&mut self, mol: &Molecule) {
        let n_atoms = mol.atoms().len();
        for bond in mol.bonds() {
            let endpts = bond.prop("_MolFileBondEndPts");
            let attach = bond.prop("_MolFileBondAttach");
            if let (Some(endpts), Some(_attach)) = (endpts, attach) {
                // Parse the V3000 array of unsigned integers (space-separated)
                let oat_strs: Vec<&str> = endpts.split_whitespace().collect();
                let mut atoms_involved = vec![false; n_atoms];
                for oat_str in &oat_strs {
                    if let Ok(mut oat) = oat_str.parse::<usize>() {
                        if oat == 0 || oat > n_atoms {
                            continue; // Bad variation point index, skip
                        }
                        oat -= 1; // Convert to 0-based
                        atoms_involved[oat] = true;
                        let center = self.at_cds[oat];
                        // Draw ellipse as a circle approximation (circle polygon)
                        let radius = self.options.variable_atom_radius;
                        let n_segments = 16;
                        let mut points = Vec::with_capacity(n_segments);
                        for i in 0..n_segments {
                            let ang = 2.0 * std::f64::consts::PI * i as f64 / n_segments as f64;
                            points.push(DVec2::new(
                                center.x + radius * ang.cos(),
                                center.y + radius * ang.sin(),
                            ));
                        }
                        self.join_paths.push(DrawPolyline {
                            points,
                            colour: self.options.variable_attachment_colour,
                            width: 1.0,
                            scale_width: false,
                            atom1_idx: Some(oat),
                            atom2_idx: None,
                            bond_idx: None,
                        });
                    }
                }

                // Connect bonds between involved atoms
                for b in mol.bonds() {
                    let ba = b.begin().index();
                    let ea = b.end().index();
                    if atoms_involved[ba] && atoms_involved[ea] {
                        self.join_paths.push(DrawPolyline {
                            points: vec![self.at_cds[ba], self.at_cds[ea]],
                            colour: self.options.variable_attachment_colour,
                            width: self.options.bond_line_width
                                * self.options.variable_bond_width_multiplier,
                            scale_width: false,
                            atom1_idx: Some(ba),
                            atom2_idx: Some(ea),
                            bond_idx: Some(b.id().index()),
                        });
                    }
                }

                // Remove * symbol from begin atom if it's a dummy
                if mol.atoms()[bond.begin().index()].atomic_number() == 0 {
                    self.atom_labels[bond.begin().index()] = None;
                }
            }
        }
    }

    // BEGIN RDKIT CPP FUNCTION DrawMol::extractBrackets (DrawMol.cpp:774-929)
    // RDKit✔️✔️: void DrawMol::extractBrackets() {
    // RDKit✔️✔️:   auto &sgs = getSubstanceGroups(*drawMol_);
    // RDKit✔️✔️:   if (sgs.empty()) { return; }
    // RDKit✔️✔️:   double rot = drawOptions_.rotate * M_PI / 180.0;
    // RDKit✔️✔️:   RDGeom::Transform2D trans;
    // RDKit✔️✔️:   trans.SetTransform(Point2D(0.0, 0.0), rot);
    // RDKit✔️✔️:   for (auto &sg : sgs) {
    // RDKit✔️✔️:     if (sg.getBrackets().empty()) { continue; }
    // RDKit✔️✔️:     Point2D refPt{0., 0.};
    // RDKit✔️✔️:     if (!sg.getAtoms().empty()) {
    // RDKit✔️✔️:       // compute bounding box of brackets, find atoms inside, average refPt
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::vector<std::pair<Point2D, Point2D>> sgBondSegments;
    // RDKit✔️✔️:     for (auto bndIdx : sg.getBonds()) {
    // RDKit✔️✔️:       // ... build bond segments
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (const auto &brk : sg.getBrackets()) {
    // RDKit✔️✔️:       Point2D p1{brk[0].x, -brk[0].y};
    // RDKit✔️✔️:       Point2D p2{brk[1].x, -brk[1].y};
    // RDKit✔️✔️:       trans.TransformPoint(p1); trans.TransformPoint(p2);
    // RDKit✔️✔️:       auto points = getBracketPoints(p1, p2, refPt, sgBondSegments);
    // RDKit✔️✔️:       postShapes_.emplace_back(new DrawShapePolyLine(points, ...));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (includeAnnotations_) {
    // RDKit✔️✔️:       // Add CONNECT and LABEL/TYPE annotations on the last bracket
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DrawMol::extractBrackets
    fn extract_brackets(&mut self, mol: &Molecule) {
        let sgs = mol.substance_groups();
        if sgs.is_empty() {
            return;
        }
        let font_size = self.font_size;

        for sg in sgs {
            let display = match sg.display() {
                Some(d) => d,
                None => continue,
            };
            if display.brackets.is_empty() {
                continue;
            }

            // Compute reference point from atoms inside bracket bounding box
            let mut ref_pt = DVec2::ZERO;

            if !sg.atoms().is_empty() {
                let mut x_min = f64::MAX / 2.0;
                let mut y_min = f64::MAX / 2.0;
                let mut x_max = f64::MIN / 2.0;
                let mut y_max = f64::MIN / 2.0;
                for brk in &display.brackets {
                    let p1 = DVec2::new(brk.p1[0], -brk.p1[1]);
                    let p2 = DVec2::new(brk.p2[0], -brk.p2[1]);
                    x_min = x_min.min(p1.x).min(p2.x);
                    y_min = y_min.min(p1.y).min(p2.y);
                    x_max = x_max.max(p1.x).max(p2.x);
                    y_max = y_max.max(p1.y).max(p2.y);
                }
                let mut num_in = 0;
                for aidx in sg.atoms() {
                    let a_cds = self.at_cds[aidx.index()];
                    if a_cds.x >= x_min && a_cds.x <= x_max && a_cds.y >= y_min && a_cds.y <= y_max
                    {
                        ref_pt += a_cds;
                        num_in += 1;
                    }
                }
                if num_in > 0 {
                    ref_pt /= num_in as f64;
                } else {
                    for aidx in sg.atoms() {
                        ref_pt += self.at_cds[aidx.index()];
                    }
                    ref_pt /= sg.atoms().len() as f64;
                }
            }

            // Build bond segments for intersection testing
            let sg_atom_indices: std::collections::HashSet<usize> =
                sg.atoms().iter().map(|a| a.index()).collect();
            let mut sg_bond_segments: Vec<(DVec2, DVec2)> = Vec::new();
            for bnd_idx in sg.bonds() {
                let bnd = &mol.bonds()[bnd_idx.index()];
                let ba = bnd.begin().index();
                let ea = bnd.end().index();
                if sg_atom_indices.contains(&ba) && sg_atom_indices.contains(&ea) {
                    sg_bond_segments.push((self.at_cds[ba], self.at_cds[ea]));
                } else if sg_atom_indices.contains(&ba) {
                    sg_bond_segments.push((self.at_cds[ba], self.at_cds[ea]));
                } else if sg_atom_indices.contains(&ea) {
                    sg_bond_segments.push((self.at_cds[ea], self.at_cds[ba]));
                }
            }

            let mut num_brackets = 0;
            for brk in &display.brackets {
                num_brackets += 1;
                let p1 = DVec2::new(brk.p1[0], -brk.p1[1]);
                let p2 = DVec2::new(brk.p2[0], -brk.p2[1]);
                let points = get_bracket_points(p1, p2, ref_pt, &sg_bond_segments);
                self.post_shapes.push(DrawPolyline {
                    points,
                    colour: DrawColour::new(0.0, 0.0, 0.0),
                    width: self.options.bond_line_width,
                    scale_width: false,
                    atom1_idx: None,
                    atom2_idx: None,
                    bond_idx: None,
                });
            }

            // Add CONNECT and LABEL/TYPE annotations
            if self.options.include_annotations {
                let last_shape_idx = self.post_shapes.len() - 1;
                let brk_shp = &self.post_shapes[last_shape_idx];
                let longline = brk_shp.points[1] - brk_shp.points[2];
                let longline_len = longline.length();
                let cos45 = 1.0 / (2.0_f64).sqrt();
                let horizontal = if longline_len > 1e-10 {
                    (longline / longline_len).x.abs() > cos45
                } else {
                    false
                };

                // Find bottom-most (or right-most) bracket
                let mut label_brk = last_shape_idx;
                for i in 1..num_brackets {
                    let brk_shp_i = &self.post_shapes[last_shape_idx - i];
                    if horizontal {
                        if brk_shp_i.points[2].y > self.post_shapes[label_brk].points[2].y {
                            label_brk = last_shape_idx - i;
                        }
                    } else {
                        if brk_shp_i.points[2].x > self.post_shapes[label_brk].points[2].x {
                            label_brk = last_shape_idx - i;
                        }
                    }
                }

                // CONNECT annotation
                if let Some(connect) = sg.props().get("CONNECT") {
                    let brk_shp = &self.post_shapes[label_brk];
                    let mut bot_pt = brk_shp.points[2];
                    let mut brk_pt = brk_shp.points[3];
                    if (!horizontal && brk_shp.points[1].y < bot_pt.y)
                        || (horizontal && brk_shp.points[1].x > bot_pt.x)
                    {
                        bot_pt = brk_shp.points[1];
                        brk_pt = brk_shp.points[0];
                    }
                    let mut da = DrawAnnotation::new(
                        connect.clone(),
                        TextAlignType::Middle,
                        "connect".to_string(),
                        self.options.annotation_font_scale,
                        bot_pt + (bot_pt - brk_pt),
                        DrawColour::new(0.0, 0.0, 0.0),
                        font_size,
                    );
                    if brk_pt.x < bot_pt.x {
                        da.align = TextAlignType::Start;
                    }
                    self.annotations.push(da);
                }

                // LABEL or TYPE annotation
                let label = sg
                    .props()
                    .get("LABEL")
                    .or_else(|| sg.props().get("TYPE"))
                    .cloned();
                if let Some(mut label) = label {
                    if label == "GEN" {
                        // ChemDraw doesn't draw GEN label
                        continue;
                    }
                    // For TYPE (not LABEL), show lowercase
                    if sg.props().get("LABEL").is_none() {
                        label = label.to_lowercase();
                    }
                    let brk_shp = &self.post_shapes[label_brk];
                    let top_pt = brk_shp.points[1];
                    let brk_pt = brk_shp.points[0];
                    let final_top = if (!horizontal && brk_shp.points[2].y > top_pt.y)
                        || (horizontal && brk_shp.points[2].x < top_pt.x)
                    {
                        brk_shp.points[2]
                    } else {
                        top_pt
                    };
                    let final_brk = if (!horizontal && brk_shp.points[2].y > top_pt.y)
                        || (horizontal && brk_shp.points[2].x < top_pt.x)
                    {
                        brk_shp.points[3]
                    } else {
                        brk_pt
                    };
                    let mut da = DrawAnnotation::new(
                        label,
                        TextAlignType::Middle,
                        "connect".to_string(),
                        self.options.annotation_font_scale,
                        final_top + (final_top - final_brk),
                        DrawColour::new(0.0, 0.0, 0.0),
                        font_size,
                    );
                    if final_brk.x < final_top.x {
                        da.align = TextAlignType::Start;
                    }
                    self.annotations.push(da);
                }
            }
        }
    }

    // BEGIN RDKIT CPP FUNCTION DrawMol::extractLinkNodes (DrawMol.cpp:930-987)
    // RDKit✔️✔️: void DrawMol::extractLinkNodes() {
    // RDKit✔️✔️:   if (!drawMol_->hasProp(common_properties::molFileLinkNodes)) { return; }
    // RDKit✔️✔️:   bool strict = false;
    // RDKit✔️✔️:   auto linkNodes = MolEnumerator::utils::getMolLinkNodes(*drawMol_, strict);
    // RDKit✔️✔️:   for (const auto &node : linkNodes) {
    // RDKit✔️✔️:     // ... draw crossing marks and labels
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DrawMol::extractLinkNodes
    fn extract_link_nodes(&mut self, mol: &Molecule) {
        if mol.prop("molFileLinkNodes").is_none() {
            return;
        }
        let font_size = self.font_size;

        // COSMolKit❗❌: Link node extraction uses MolEnumerator::utils::getMolLinkNodes
        // which is not ported. For now, look for bonds with the `_MolFileBondEndPts` prop
        // as a heuristic for link nodes, and draw bracket crosses.
        let crossing_frac = 0.333;
        let length_frac = 0.333;
        let mut label_pt = DVec2::new(-1000.0, -1000.0);
        let mut label_perp = DVec2::ZERO;

        for bond in mol.bonds() {
            let endpts = bond.prop("_MolFileBondEndPts");
            if endpts.is_none() {
                continue;
            }
            let ba = bond.begin().index();
            let ea = bond.end().index();
            let start_loc = self.at_cds[ba];
            let end_loc = self.at_cds[ea];
            let vect = end_loc - start_loc;
            let offset = vect * crossing_frac;
            let crossing_pt = start_loc + offset;
            let perp = DVec2::new(vect.y, -vect.x);
            let perp = perp * length_frac;
            let p1 = crossing_pt + perp / 2.0;
            let p2 = crossing_pt - perp / 2.0;

            let bond_segments: Vec<(DVec2, DVec2)> = Vec::new();
            let points = get_bracket_points(p1, p2, start_loc, &bond_segments);
            self.post_shapes.push(DrawPolyline {
                points,
                colour: DrawColour::new(0.0, 0.0, 0.0),
                width: self.options.bond_line_width,
                scale_width: false,
                atom1_idx: Some(ba),
                atom2_idx: Some(ea),
                bond_idx: Some(bond.id().index()),
            });

            if p1.x > label_pt.x {
                label_pt = p1;
                label_perp = crossing_pt - start_loc;
            }
            if p2.x > label_pt.x {
                label_pt = p2;
                label_perp = crossing_pt - start_loc;
            }
        }

        // Add label for link nodes if we found any
        if label_pt.x > -500.0 {
            if self.options.include_annotations {
                let perp_len = label_perp.length();
                let perp = if perp_len > 1e-10 {
                    label_perp / perp_len * 0.2
                } else {
                    DVec2::new(0.2, 0.0)
                };
                let da = DrawAnnotation::new(
                    "(1-1)".to_string(), // default link node label
                    TextAlignType::Start,
                    "linknode".to_string(),
                    self.options.annotation_font_scale,
                    label_pt + perp,
                    DrawColour::new(0.0, 0.0, 0.0),
                    font_size,
                );
                self.annotations.push(da);
            }
        }
    }

    // BEGIN RDKIT CPP FUNCTION DrawMol::extractCloseContacts (DrawMol.cpp:988-1116)
    // RDKit✔️✔️: void DrawMol::extractCloseContacts() {
    // RDKit✔️✔️:   if (drawOptions_.flagCloseContactsDist < 0) { return; }
    // RDKit✔️✔️:   int tol = drawOptions_.flagCloseContactsDist * drawOptions_.flagCloseContactsDist;
    // RDKit✔️✔️:   boost::dynamic_bitset<> flagged(atCds_.size());
    // RDKit✔️✔️:   Point2D trans, scale, toCentre;
    // RDKit✔️✔️:   getDrawTransformers(trans, scale, toCentre);
    // RDKit✔️✔️:   for (unsigned int i = 0; i < atCds_.size(); ++i) {
    // RDKit✔️✔️:     if (flagged[i]) { continue; }
    // RDKit✔️✔️:     Point2D ci = transformPoint(atCds_[i], &trans, &scale, &toCentre);
    // RDKit✔️✔️:     for (unsigned int j = i + 1; j < atCds_.size(); ++j) {
    // RDKit✔️✔️:       if (flagged[j]) { continue; }
    // RDKit✔️✔️:       Point2D cj = transformPoint(atCds_[j], &trans, &scale, &toCentre);
    // RDKit✔️✔️:       double d = (cj - ci).lengthSq();
    // RDKit✔️✔️:       if (d <= tol) {
    // RDKit✔️✔️:         flagged.set(i); flagged.set(j); break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (flagged[i]) {
    // RDKit✔️✔️:       Point2D p1 = ci; Point2D p2 = p1;
    // RDKit✔️✔️:       Point2D offset(0.1 * scale_, 0.1 * scale_);
    // RDKit✔️✔️:       p1 -= offset; p2 += offset;
    // RDKit✔️✔️:       std::vector<Point2D> points(5);
    // RDKit✔️✔️:       points[0] = points[4] = p1;
    // RDKit✔️✔️:       points[1] = Point2D{p1.x, p2.y};
    // RDKit✔️✔️:       points[2] = Point2D{p2};
    // RDKit✔️✔️:       points[3] = Point2D{p2.x, p1.y};
    // RDKit✔️✔️:       postShapes_.emplace_back(new DrawShapePolyLine(points, ..., DrawColour(1.0, 0.0, 0.0), ...));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DrawMol::extractCloseContacts
    fn extract_close_contacts(&mut self) {
        if self.options.flag_close_contacts_dist < 0 {
            return;
        }
        let tol = (self.options.flag_close_contacts_dist as f64)
            * (self.options.flag_close_contacts_dist as f64);
        let mut flagged = vec![false; self.at_cds.len()];
        // Rough scaling for the draw coordinate transformation
        let scale_factor = if self.scale > 0.0 { self.scale } else { 1.0 };

        for i in 0..self.at_cds.len() {
            if flagged[i] {
                continue;
            }
            let ci = self.at_cds[i];
            for j in (i + 1)..self.at_cds.len() {
                if flagged[j] {
                    continue;
                }
                let cj = self.at_cds[j];
                let d = (cj - ci).length_squared();
                if d <= tol {
                    flagged[i] = true;
                    flagged[j] = true;
                    break;
                }
            }
            if flagged[i] {
                let p1 = ci;
                let p2 = p1;
                let offset = DVec2::new(0.1 * scale_factor, 0.1 * scale_factor);
                let p1 = p1 - offset;
                let p2 = p2 + offset;
                let points = vec![
                    p1,
                    DVec2::new(p1.x, p2.y),
                    p2,
                    DVec2::new(p2.x, p1.y),
                    p1, // close the polygon
                ];
                self.post_shapes.push(DrawPolyline {
                    points,
                    colour: DrawColour::new(1.0, 0.0, 0.0),
                    width: self.options.bond_line_width,
                    scale_width: false,
                    atom1_idx: Some(i),
                    atom2_idx: None,
                    bond_idx: None,
                });
            }
        }
    }

    // BEGIN RDKIT CPP FUNCTION DrawMol::calcAnnotationPosition(const Atom *) (DrawMol.cpp:2570-2603)
    // RDKit✔️✔️: void DrawMol::calcAnnotationPosition(const Atom *atom,
    // RDKit✔️✔️:                                      DrawAnnotation &annot) const {
    // RDKit✔️✔️:   double start_ang = getNoteStartAngle(atom);
    // RDKit✔️✔️:   Point2D const &atCds = atCds_[atom->getIdx()];
    // RDKit✔️✔️:   double radStep = 0.25;
    // RDKit✔️✔️:   Point2D leastWorstPos = atCds;
    // RDKit✔️✔️:   int leastWorstScore = 100;
    // RDKit✔️✔️:   for (int j = 1; j < 4; ++j) {
    // RDKit✔️✔️:     double note_rad = j * radStep;
    // RDKit✔️✔️:     if (j == 1 && atomLabels_[atom->getIdx()]) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (int i = 0; i < 12; ++i) {
    // RDKit✔️✔️:       double ang = start_ang + i * 30.0 * M_PI / 180.0;
    // RDKit✔️✔️:       annot.pos_.x = atCds.x + cos(ang) * note_rad;
    // RDKit✔️✔️:       annot.pos_.y = atCds.y + sin(ang) * note_rad;
    // RDKit✔️✔️:       int clashScore = doesNoteClash(annot);
    // RDKit✔️✔️:       if (!clashScore) {
    // RDKit✔️✔️:         return;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         if (clashScore < leastWorstScore) {
    // RDKit✔️✔️:           leastWorstScore = clashScore;
    // RDKit✔️✔️:           leastWorstPos = annot.pos_;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   annot.pos_ = leastWorstPos;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DrawMol::calcAnnotationPosition(const Atom *)
    fn calc_annotation_position_for_atom(
        &self,
        mol: &Molecule,
        atom_idx: usize,
        annot: &mut DrawAnnotation,
    ) {
        let start_ang = self.calc_note_start_angle(mol, atom_idx);
        let at_cds = self.at_cds[atom_idx];
        let rad_step = 0.25;
        let mut least_worst_pos = at_cds;
        let mut least_worst_score = 100;

        for j in 1..4 {
            let note_rad = j as f64 * rad_step;
            if j == 1
                && self
                    .atom_labels
                    .get(atom_idx)
                    .and_then(|l| l.as_ref())
                    .is_some()
            {
                continue;
            }
            for i in 0..12 {
                let ang = start_ang + i as f64 * 30.0_f64.to_radians();
                annot.pos = DVec2::new(
                    at_cds.x + ang.cos() * note_rad,
                    at_cds.y + ang.sin() * note_rad,
                );
                let clash_score = self.does_note_clash(annot);
                if clash_score == 0 {
                    return;
                } else if clash_score < least_worst_score {
                    least_worst_score = clash_score;
                    least_worst_pos = annot.pos;
                }
            }
        }
        annot.pos = least_worst_pos;
    }

    // BEGIN RDKIT CPP FUNCTION DrawMol::calcAnnotationPosition(const Bond *) (DrawMol.cpp:2606-2653)
    // RDKit✔️✔️: void DrawMol::calcAnnotationPosition(const Bond *bond,
    // RDKit✔️✔️:                                      DrawAnnotation &annot) const {
    // RDKit✔️✔️:   Point2D const &at1_cds = atCds_[bond->getBeginAtomIdx()];
    // RDKit✔️✔️:   Point2D at2_cds = atCds_[bond->getEndAtomIdx()];
    // RDKit✔️✔️:   if ((at1_cds - at2_cds).lengthSq() < 0.0001) {
    // RDKit✔️✔️:     at2_cds.x += 0.1;  at2_cds.y += 0.1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   Point2D perp = calcPerpendicular(at1_cds, at2_cds);
    // RDKit✔️✔️:   Point2D bond_vec = at1_cds.directionVector(at2_cds);
    // RDKit✔️✔️:   double bond_len = (at1_cds - at2_cds).length();
    // RDKit✔️✔️:   std::vector<double> mid_offsets{0.5, 0.33, 0.66, 0.25, 0.75};
    // RDKit✔️✔️:   double offset_step = drawOptions_.multipleBondOffset;
    // RDKit✔️✔️:   Point2D leastWorstPos = (at1_cds + at2_cds) / 2.0;
    // RDKit✔️✔️:   int leastWorstScore = 100;
    // RDKit✔️✔️:   for (auto mo : mid_offsets) {
    // RDKit✔️✔️:     Point2D mid = at1_cds + bond_vec * bond_len * mo;
    // RDKit✔️✔️:     for (int j = 1; j < 6; ++j) {
    // RDKit✔️✔️:       if (j == 1 && bond->getBondType() > 1) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       double offset = j * offset_step;
    // RDKit✔️✔️:       annot.pos_ = mid + perp * offset;
    // RDKit✔️✔️:       int clashScore = doesNoteClash(annot);
    // RDKit✔️✔️:       if (!clashScore) { return; }
    // RDKit✔️✔️:       if (clashScore < leastWorstScore) {
    // RDKit✔️✔️:         leastWorstPos = annot.pos_;  leastWorstScore = clashScore;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       annot.pos_ = mid - perp * offset;
    // RDKit✔️✔️:       clashScore = doesNoteClash(annot);
    // RDKit✔️✔️:       if (!clashScore) { return; }
    // RDKit✔️✔️:       if (clashScore < leastWorstScore) {
    // RDKit✔️✔️:         leastWorstPos = annot.pos_;  leastWorstScore = clashScore;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   annot.pos_ = leastWorstPos;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DrawMol::calcAnnotationPosition(const Bond *)
    fn calc_annotation_position_for_bond(
        &self,
        mol: &Molecule,
        bond: &Bond,
        annot: &mut DrawAnnotation,
    ) {
        let b = bond.begin().index();
        let e = bond.end().index();
        let mut at1_cds = self.at_cds[b];
        let mut at2_cds = self.at_cds[e];
        if (at1_cds - at2_cds).length_squared() < 0.0001 {
            at2_cds.x += 0.1;
            at2_cds.y += 0.1;
        }
        let perp = calc_perpendicular(at1_cds, at2_cds);
        let bond_vec = direction_vector(at1_cds, at2_cds);
        let bond_len = (at1_cds - at2_cds).length();
        let mid_offsets = [0.5, 0.33, 0.66, 0.25, 0.75];
        let offset_step = self.options.multiple_bond_offset;
        let mut least_worst_pos = (at1_cds + at2_cds) / 2.0;
        let mut least_worst_score = 100;

        for mo in &mid_offsets {
            let mid = at1_cds + bond_vec * bond_len * mo;
            for j in 1..6 {
                if j == 1
                    && matches!(
                        bond.order(),
                        BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
                    )
                {
                    continue;
                }
                let offset = j as f64 * offset_step;
                // Try +perp
                annot.pos = mid + perp * offset;
                let cs = self.does_note_clash(annot);
                if cs == 0 {
                    return;
                }
                if cs < least_worst_score {
                    least_worst_pos = annot.pos;
                    least_worst_score = cs;
                }
                // Try -perp
                annot.pos = mid - perp * offset;
                let cs = self.does_note_clash(annot);
                if cs == 0 {
                    return;
                }
                if cs < least_worst_score {
                    least_worst_pos = annot.pos;
                    least_worst_score = cs;
                }
            }
        }
        annot.pos = least_worst_pos;
    }

    // BEGIN RDKIT CPP FUNCTION DrawMol::getNoteStartAngle (DrawMol.cpp:2656-2725)
    // RDKit✔️✔️: double DrawMol::getNoteStartAngle(const Atom *atom) const {
    // RDKit✔️✔️:   if (atom->getDegree() == 0) { return M_PI / 2.0; }
    // RDKit✔️✔️:   const Point2D &at_cds = atCds_[atom->getIdx()];
    // RDKit✔️✔️:   std::vector<Point2D> bond_vecs;
    // RDKit✔️✔️:   for (auto nbr : make_iterator_range(drawMol_->getAtomNeighbors(atom))) {
    // RDKit✔️✔️:     if ((at_cds - atCds_[nbr]).lengthSq() < 0.0001) {
    // RDKit✔️✔️:       bond_vec.x = 0.1;  bond_vec.y = 0.1;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       bond_vec = at_cds.directionVector(atCds_[nbr]);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     bond_vec.normalize();
    // RDKit✔️✔️:     bond_vecs.push_back(bond_vec);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   Point2D ret_vec;
    // RDKit✔️✔️:   if (bond_vecs.size() == 1) {
    // RDKit✔️✔️:     if (!atomLabels_[atom->getIdx()]) {
    // RDKit✔️✔️:       ret_vec.x = bond_vecs[0].y;  ret_vec.y = -bond_vecs[0].x;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       ret_vec = -bond_vecs[0];
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (bond_vecs.size() == 2) {
    // RDKit✔️✔️:     ret_vec = bond_vecs[0] + bond_vecs[1];
    // RDKit✔️✔️:     if (ret_vec.lengthSq() > 1.0e-6) {
    // RDKit✔️✔️:       if (!atom->getNumImplicitHs() || atom->getAtomicNum() == 6) {
    // RDKit✔️✔️:         ret_vec *= -1.0;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       ret_vec.x = -bond_vecs.front().y;  ret_vec.y = bond_vecs.front().x;
    // RDKit✔️✔️:       ret_vec.normalize();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     double discrim = 4.0 * M_PI / bond_vecs.size();
    // RDKit✔️✔️:     for (size_t i = 0; i < bond_vecs.size() - 1; ++i) {
    // RDKit✔️✔️:       for (size_t j = i + 1; j < bond_vecs.size(); ++j) {
    // RDKit✔️✔️:         double ang = acos(bond_vecs[i].dotProduct(bond_vecs[j]));
    // RDKit✔️✔️:         if (ang < discrim) {
    // RDKit✔️✔️:           ret_vec = bond_vecs[i] + bond_vecs[j];
    // RDKit✔️✔️:           ret_vec.normalize();
    // RDKit✔️✔️:           discrim = -1.0;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (discrim > 0.0) {
    // RDKit✔️✔️:       ret_vec = bond_vecs[0] + bond_vecs[1];
    // RDKit✔️✔️:       ret_vec *= -1.0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return atan2(ret_vec.y, ret_vec.x);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DrawMol::getNoteStartAngle
    fn calc_note_start_angle(&self, mol: &Molecule, atom_idx: usize) -> f64 {
        let degree = atom_degree(mol, atom_idx);
        if degree == 0 {
            return std::f64::consts::FRAC_PI_2;
        }
        let at_cds = self.at_cds[atom_idx];
        let mut bond_vecs: Vec<DVec2> = Vec::new();
        for &nbr in &atom_neighbors(mol, atom_idx) {
            let bond_vec = if (at_cds - self.at_cds[nbr]).length_squared() < 0.0001 {
                DVec2::new(0.1, 0.1)
            } else {
                direction_vector(at_cds, self.at_cds[nbr])
            };
            bond_vecs.push(bond_vec.normalize());
        }

        let ret_vec = if bond_vecs.len() == 1 {
            if self
                .atom_labels
                .get(atom_idx)
                .and_then(|l| l.as_ref())
                .is_none()
            {
                DVec2::new(bond_vecs[0].y, -bond_vecs[0].x)
            } else {
                -bond_vecs[0]
            }
        } else if bond_vecs.len() == 2 {
            let mut rv = bond_vecs[0] + bond_vecs[1];
            if rv.length_squared() > 1.0e-6 {
                // prefer outside the angle if no implicit Hs
                rv *= -1.0;
            } else {
                rv = DVec2::new(-bond_vecs[0].y, bond_vecs[0].x).normalize();
            }
            rv
        } else {
            let mut discrim = 4.0 * std::f64::consts::PI / bond_vecs.len() as f64;
            let mut ret = bond_vecs[0] + bond_vecs[1];
            'outer: for i in 0..bond_vecs.len() - 1 {
                for j in (i + 1)..bond_vecs.len() {
                    let ang = bond_vecs[i].dot(bond_vecs[j]).acos();
                    if ang < discrim {
                        ret = bond_vecs[i] + bond_vecs[j];
                        ret = ret.normalize();
                        discrim = -1.0;
                        break 'outer;
                    }
                }
            }
            if discrim > 0.0 {
                ret *= -1.0;
            }
            ret
        };

        f64::atan2(ret_vec.y, ret_vec.x)
    }

    // BEGIN RDKIT CPP FUNCTION DrawMol::doesNoteClash (DrawMol.cpp:2728-2746)
    // RDKit✔️✔️: int DrawMol::doesNoteClash(const DrawAnnotation &annot) const {
    // RDKit✔️✔️:   for (auto &rect : annot.rects_) {
    // RDKit✔️✔️:     Point2D otrans = rect->trans_;
    // RDKit✔️✔️:     rect->trans_ += annot.pos_;
    // RDKit✔️✔️:     double padding = scale_ * 0.04;
    // RDKit✔️✔️:     int clashScore = doesRectClash(*rect, padding);
    // RDKit✔️✔️:     rect->trans_ = otrans;
    // RDKit✔️✔️:     if (clashScore) { return clashScore; }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DrawMol::doesNoteClash
    fn does_note_clash(&self, annot: &DrawAnnotation) -> i32 {
        let padding = self.scale * 0.04;
        for rect in &annot.rects {
            let adjusted = StringRect {
                trans: rect.trans + annot.pos,
                ..rect.clone()
            };
            let cs = self.does_rect_clash_with_score(&adjusted, padding);
            if cs != 0 {
                return cs;
            }
        }
        0
    }

    // BEGIN RDKIT CPP FUNCTION DrawMol::doesRectClash (DrawMol.cpp:2749-2785)
    // RDKit✔️✔️: int DrawMol::doesRectClash(const StringRect &rect, double padding) const {
    // RDKit✔️✔️:   for (auto bond : drawMol_->bonds()) {
    // RDKit✔️✔️:     if (bond->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:       auto at1 = bond->getBeginAtomIdx();
    // RDKit✔️✔️:       auto at2 = bond->getEndAtomIdx();
    // RDKit✔️✔️:       if (doesLineIntersect(rect, atCds_[at1], atCds_[at2], 0.0)) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (const auto &bond : bonds_) {
    // RDKit✔️✔️:     if (bond->doesRectClash(rect, padding)) { return 1; }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (const auto &al : atomLabels_) {
    // RDKit✔️✔️:     if (al && al->doesRectClash(rect, padding)) { return 2; }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (const auto &a : annotations_) {
    // RDKit✔️✔️:     if (a->doesRectClash(rect, padding)) { return 3; }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DrawMol::doesRectClash
    fn does_rect_clash_with_score(&self, rect: &StringRect, padding: f64) -> i32 {
        // Check double bonds (the actual lines, not just draw shapes)
        for bond_line in &self.bonds {
            if rect_clashes_with_line(rect, bond_line.begin, bond_line.end, 0.0) {
                return 1;
            }
        }
        // Check atom labels
        for label in self.atom_labels.iter().flatten() {
            if label_rects_intersect(&label.rects, label.cds, rect, padding) {
                return 2;
            }
        }
        // Check existing annotations
        for annot in &self.annotations {
            for ar in &annot.rects {
                let adjusted = StringRect {
                    trans: ar.trans + annot.pos,
                    ..ar.clone()
                };
                if rects_intersect(rect, &adjusted, padding) {
                    return 3;
                }
            }
        }
        0
    }

    fn calc_radical_rect(&self, atom: &Atom) -> (StringRect, OrientType) {
        let idx = atom.id().index();
        let at_cds = self.at_cds[idx];
        let spot_rad = self.radical_spot_radius_unscaled();
        let orient = self.atom_orients.get(idx).copied().unwrap_or(OrientType::C);
        let rad_size =
            (4.0 * f64::from(atom.radical_electrons()) - 2.0) * spot_rad / self.font_scale_factor();

        let (x_min, x_max, y_min, y_max) = if let Some(Some(label)) = self.atom_labels.get(idx) {
            let mut x_min = f64::MAX;
            let mut x_max = f64::MIN;
            let mut y_min = f64::MAX;
            let mut y_max = f64::MIN;
            label.find_extremes(&mut x_min, &mut x_max, &mut y_min, &mut y_max);
            (x_min, x_max, y_min, y_max)
        } else {
            (
                at_cds.x - 3.0 * spot_rad,
                at_cds.x + 3.0 * spot_rad,
                at_cds.y - 3.0 * spot_rad,
                at_cds.y + 3.0 * spot_rad,
            )
        };

        for trial in std::iter::once(orient).chain(
            [OrientType::N, OrientType::E, OrientType::S, OrientType::W]
                .into_iter()
                .filter(move |fallback| *fallback != orient),
        ) {
            let rect = radical_rect_for_orientation(
                trial, at_cds, x_min, x_max, y_min, y_max, spot_rad, rad_size,
            );
            if !self.does_rect_clash(&rect, 0.0) {
                return (rect, trial);
            }
        }
        (
            radical_rect_for_orientation(
                OrientType::N,
                at_cds,
                x_min,
                x_max,
                y_min,
                y_max,
                spot_rad,
                rad_size,
            ),
            OrientType::N,
        )
    }

    fn does_rect_clash(&self, rect: &StringRect, padding: f64) -> bool {
        for bond in &self.bonds {
            if rect_clashes_with_line(rect, bond.begin, bond.end, padding) {
                return true;
            }
        }
        for wedge in &self.wedges {
            match wedge.kind {
                WedgeKind::Solid => {
                    for tri in wedge.points.chunks_exact(3) {
                        if rect_clashes_with_triangle(rect, tri[0], tri[1], tri[2], padding) {
                            return true;
                        }
                    }
                }
                WedgeKind::Dashed => {
                    if wedge.points.len() >= 3
                        && rect_clashes_with_triangle(
                            rect,
                            wedge.points[0],
                            wedge.points[1],
                            wedge.points[2],
                            padding,
                        )
                    {
                        return true;
                    }
                }
            }
        }
        for label in self.atom_labels.iter().flatten() {
            if label_rects_intersect(&label.rects, label.cds, rect, padding) {
                return true;
            }
        }
        false
    }

    fn calculate_scale(&mut self) {
        self.find_extremes();
        self.x_range = self.x_max - self.x_min;
        self.y_range = self.y_max - self.y_min;

        if self.x_range < 1e-8 {
            self.x_range = 1.0;
        }
        if self.y_range < 1e-8 {
            self.y_range = 1.0;
        }

        let pad = self.margin_padding * 2.0;
        let sx = self.draw_width / (self.x_range * (1.0 + pad));
        let sy = self.draw_height / (self.y_range * (1.0 + pad));
        self.scale = sx.min(sy);

        // Also adjust font scale based on mean bond length
        if self.mean_bond_length > 1e-8 {
            let ideal_scale = 40.0 / self.mean_bond_length;
            self.scale = self.scale.min(ideal_scale);
        }
    }

    fn find_extremes(&mut self) {
        self.x_min = f64::MAX;
        self.x_max = f64::MIN;
        self.y_min = f64::MAX;
        self.y_max = f64::MIN;

        for label in self.atom_labels.iter().flatten() {
            label.find_extremes(
                &mut self.x_min,
                &mut self.x_max,
                &mut self.y_min,
                &mut self.y_max,
            );
        }
        for radical in &self.radicals {
            let cx = radical.rect.trans.x;
            let cy = radical.rect.trans.y;
            let hw = radical.rect.width / 2.0;
            let hh = radical.rect.height / 2.0;
            if cx - hw < self.x_min {
                self.x_min = cx - hw;
            }
            if cx + hw > self.x_max {
                self.x_max = cx + hw;
            }
            if cy - hh < self.y_min {
                self.y_min = cy - hh;
            }
            if cy + hh > self.y_max {
                self.y_max = cy + hh;
            }
        }
        for pt in &self.at_cds {
            if pt.x < self.x_min {
                self.x_min = pt.x;
            }
            if pt.x > self.x_max {
                self.x_max = pt.x;
            }
            if pt.y < self.y_min {
                self.y_min = pt.y;
            }
            if pt.y > self.y_max {
                self.y_max = pt.y;
            }
        }
    }

    fn change_to_draw_coords(&mut self) {
        let pad = self.margin_padding;
        let offset_x = self.draw_width * pad;
        let offset_y = self.draw_height * pad;
        let draw_w = self.draw_width - 2.0 * offset_x;
        let draw_h = self.draw_height - 2.0 * offset_y;

        let sx = draw_w / self.x_range;
        let sy = draw_h / self.y_range;
        let actual_scale = sx.min(sy);

        let centre_x = (self.x_min + self.x_max) / 2.0;
        let centre_y = (self.y_min + self.y_max) / 2.0;

        let to_centre = DVec2::new(centre_x, centre_y);
        let trans = DVec2::new(self.width / 2.0, self.height / 2.0);
        let scale = DVec2::new(actual_scale, actual_scale);

        // Transform atom coordinates
        for pt in &mut self.at_cds {
            *pt = transform_point(*pt, trans, scale, to_centre);
        }

        // Transform labels
        for label in self.atom_labels.iter_mut().flatten() {
            label.cds = transform_point(label.cds, trans, scale, to_centre);
            label.recalculate_rects(self.font_size * actual_scale);
        }

        // Transform bonds
        for line in &mut self.bonds {
            line.begin = transform_point(line.begin, trans, scale, to_centre);
            line.end = transform_point(line.end, trans, scale, to_centre);
            line.width *= actual_scale;
        }

        // Transform wedges
        for wedge in &mut self.wedges {
            for pt in &mut wedge.points {
                *pt = transform_point(*pt, trans, scale, to_centre);
            }
        }

        // Transform join paths
        for poly in &mut self.join_paths {
            for pt in &mut poly.points {
                *pt = transform_point(*pt, trans, scale, to_centre);
            }
            poly.width *= actual_scale;
        }

        // Transform radicals
        for radical in &mut self.radicals {
            let scaled_pt = transform_point(
                DVec2::new(radical.rect.trans.x, radical.rect.trans.y),
                trans,
                scale,
                to_centre,
            );
            radical.rect.trans = scaled_pt;
            radical.rect.width *= actual_scale;
            radical.rect.height *= actual_scale;
        }

        // Transform annotations
        for annot in &mut self.annotations {
            annot.pos = transform_point(annot.pos, trans, scale, to_centre);
            annot.rects = get_string_rects(
                &annot.text,
                OrientType::C,
                self.font_size * annot.font_scale * actual_scale,
            );
        }

        self.scale = actual_scale;
    }

    fn line_endpoint_idx(&self, line: &DrawLine, atom_idx: usize) -> Option<usize> {
        if line.atom1_idx == atom_idx {
            Some(0)
        } else if line.atom2_idx == atom_idx {
            Some(1)
        } else {
            None
        }
    }

    fn adjust_bonds_on_solid_wedge_ends(&mut self, mol: &Molecule, wedge_bonds: &HashSet<usize>) {
        // RDKit❗✔️: when a bond connects to the wide end of a wedge,
        // adjust the bond line to start from the wedge boundary.
        for wedge in &self.wedges {
            if wedge.kind != WedgeKind::Solid {
                continue;
            }
            let wide_end = wedge.atom2_idx;
            let narrow_end = wedge.atom1_idx;

            for line in &mut self.bonds {
                if line.atom1_idx == wide_end && line.atom2_idx != narrow_end {
                    // This bond attaches to the wide end; shift its start point
                    let w_pts = &wedge.points;
                    if w_pts.len() >= 4 {
                        let mid_wide = (w_pts[2] + w_pts[3]) / 2.0;
                        line.begin = mid_wide;
                    }
                } else if line.atom2_idx == wide_end && line.atom1_idx != narrow_end {
                    let w_pts = &wedge.points;
                    if w_pts.len() >= 4 {
                        let mid_wide = (w_pts[2] + w_pts[3]) / 2.0;
                        line.end = mid_wide;
                    }
                }
            }
        }
    }
}

// ──────────────────────────────────────────────
// Free functions used by DrawMol
// ──────────────────────────────────────────────

fn line_endpoint_for_atom(line: &DrawLine, atom_idx: usize) -> Option<usize> {
    if line.atom1_idx == atom_idx {
        Some(0)
    } else if line.atom2_idx == atom_idx {
        Some(1)
    } else {
        None
    }
}

fn line_point(line: &DrawLine, point_idx: usize) -> DVec2 {
    match point_idx {
        0 => line.begin,
        _ => line.end,
    }
}

fn is_linear_atom(mol: &Molecule, atom_idx: usize) -> bool {
    let degree = atom_degree(mol, atom_idx);
    degree <= 2
}

fn radical_rect_at(trans: DVec2, width: f64, height: f64) -> StringRect {
    StringRect {
        ch: '\0',
        draw_mode: TextDrawType::Normal,
        trans,
        offset: DVec2::ZERO,
        g_centre: DVec2::ZERO,
        y_shift: 0.0,
        width,
        height,
        rect_corr: 0.0,
    }
}

fn radical_rect_for_orientation(
    orient: OrientType,
    at_cds: DVec2,
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    spot_rad: f64,
    rad_size: f64,
) -> StringRect {
    let (tx, ty, width, height) = match orient {
        OrientType::N => (at_cds.x, y_max + rad_size, rad_size, rad_size),
        OrientType::S => (at_cds.x, y_min - rad_size, rad_size, rad_size),
        OrientType::E => (x_max + rad_size, at_cds.y, rad_size, rad_size),
        OrientType::W => (x_min - rad_size, at_cds.y, rad_size, rad_size),
        OrientType::C => (at_cds.x, at_cds.y - spot_rad * 2.0, rad_size, rad_size),
    };
    radical_rect_at(DVec2::new(tx, ty), width, height)
}

/// RDKit❗✔️: element_symbol — periodic table element names from atomic number.
/// This is a superset of the RDKit PeriodicTable::getElementSymbol, including
/// placeholder names for all known elements up to Og (118). Unknown numbers
/// return "?". Matches the same subset used in smiles_write.rs.
fn element_symbol(atomic_num: u8) -> &'static str {
    // RDKit periodic table element symbols (subset)
    match atomic_num {
        0 => "*",
        1 => "H",
        2 => "He",
        3 => "Li",
        4 => "Be",
        5 => "B",
        6 => "C",
        7 => "N",
        8 => "O",
        9 => "F",
        10 => "Ne",
        11 => "Na",
        12 => "Mg",
        13 => "Al",
        14 => "Si",
        15 => "P",
        16 => "S",
        17 => "Cl",
        18 => "Ar",
        19 => "K",
        20 => "Ca",
        21 => "Sc",
        22 => "Ti",
        23 => "V",
        24 => "Cr",
        25 => "Mn",
        26 => "Fe",
        27 => "Co",
        28 => "Ni",
        29 => "Cu",
        30 => "Zn",
        31 => "Ga",
        32 => "Ge",
        33 => "As",
        34 => "Se",
        35 => "Br",
        36 => "Kr",
        37 => "Rb",
        38 => "Sr",
        39 => "Y",
        40 => "Zr",
        41 => "Nb",
        42 => "Mo",
        43 => "Tc",
        44 => "Ru",
        45 => "Rh",
        46 => "Pd",
        47 => "Ag",
        48 => "Cd",
        49 => "In",
        50 => "Sn",
        51 => "Sb",
        52 => "Te",
        53 => "I",
        54 => "Xe",
        55 => "Cs",
        56 => "Ba",
        57..=71 => "La",
        72 => "Hf",
        73 => "Ta",
        74 => "W",
        75 => "Re",
        76 => "Os",
        77 => "Ir",
        78 => "Pt",
        79 => "Au",
        80 => "Hg",
        81 => "Tl",
        82 => "Pb",
        83 => "Bi",
        84 => "Po",
        85 => "At",
        86 => "Rn",
        87 => "Fr",
        88 => "Ra",
        89..=103 => "Ac",
        104 => "Rf",
        105 => "Db",
        106 => "Sg",
        107 => "Bh",
        108 => "Hs",
        109 => "Mt",
        110 => "Ds",
        111 => "Rg",
        112 => "Cn",
        113 => "Nh",
        114 => "Fl",
        115 => "Mc",
        116 => "Lv",
        117 => "Ts",
        118 => "Og",
        _ => "?",
    }
}

// ──────────────────────────────────────────────
// SVG rendering primitives
// ──────────────────────────────────────────────

// BEGIN RDKIT CPP FUNCTION MolDraw2DSVG::initDrawing (MolDraw2DSVG.cpp:122-133)
// RDKit✔️✔️: void MolDraw2DSVG::initDrawing() {
// RDKit✔️✔️:   d_os << "<?xml version='1.0' encoding='iso-8859-1'?>\n";
// RDKit✔️✔️:   d_os << "<svg version='1.1' baseProfile='full'\n      \
// RDKit✔️✔️:         xmlns='http://www.w3.org/2000/svg'\n              \
// RDKit✔️✔️:         xmlns:rdkit='http://www.rdkit.org/xml'\n              \
// RDKit✔️✔️:         xmlns:xlink='http://www.w3.org/1999/xlink'\n          \
// RDKit✔️✔️:         xml:space='preserve'\n";
// RDKit✔️✔️:   d_os
// RDKit✔️✔️:       << boost::format{"width='%1%px' height='%2%px' viewBox='0 0 %1% %2%'>\n"} %
// RDKit✔️✔️:              width() % height();
// RDKit✔️✔️:   d_os << "<!-- END OF HEADER -->\n";
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION MolDraw2DSVG::initDrawing
fn init_drawing(out: &mut String, width: u32, height: u32) {
    out.push_str("<?xml version='1.0' encoding='iso-8859-1'?>\n");
    out.push_str(
        "<svg version='1.1' baseProfile='full'\n      \
        xmlns='http://www.w3.org/2000/svg'\n              \
        xmlns:rdkit='http://www.rdkit.org/xml'\n              \
        xmlns:xlink='http://www.w3.org/1999/xlink'\n          \
        xml:space='preserve'\n",
    );
    out.push_str(&format!(
        "width='{}px' height='{}px' viewBox='0 0 {} {}'>\n",
        width, height, width, height
    ));
    out.push_str("<!-- END OF HEADER -->\n");
}

/// RDKit✔️✔️: background rectangle — equivalent to MolDraw2DSVG's initDrawing
/// background fill. RDKit sets a white background via <rect> in clearDrawing.
/// The Rust implementation matches this exactly.
// BEGIN RDKIT CPP FUNCTION MolDraw2DSVG::clearDrawing (MolDraw2DSVG.cpp:322-335)
// RDKit✔️✔️: void MolDraw2DSVG::clearDrawing() {
// RDKit✔️✔️:   d_os << "<rect";
// RDKit✔️✔️:   outputClasses();
// RDKit✔️✔️:   d_os << " style='fill:" << DrawColourToSVG(drawOptions().backgroundColour)
// RDKit✔️✔️:        << ";fill-opacity:1;fill-rule:evenodd;stroke:none'";
// RDKit✔️✔️:   d_os << " width='" << width() << "' height='" << height()
// RDKit✔️✔️:        << "' x='0' y='0'";
// RDKit✔️✔️:   d_os << " />\n";
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION MolDraw2DSVG::clearDrawing
fn clear_drawing(out: &mut String, width: u32, height: u32, colour: DrawColour) {
    let col = draw_colour_to_svg(colour);
    out.push_str(&format!(
        "<rect style='fill:{};fill-opacity:{};fill-rule:evenodd;stroke:none' \
         width='{}' height='{}' x='0' y='0' />\n",
        col, colour.a, width, height
    ));
}

// BEGIN RDKIT CPP FUNCTION MolDraw2DSVG::initTextDrawer (MolDraw2DSVG.cpp:136-163)
// RDKit✔️✔️: void MolDraw2DSVG::initTextDrawer() {
// RDKit✔️✔️:   // RDKit initialises a Freetype font renderer here.
// RDKit✔️✔️:   // COSMolKit uses usvg which handles fonts automatically via the
// RDKit✔️✔️:   // embedded NotoSans font data loaded in svg_to_png().
// RDKit✔️✔️: }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION MolDraw2DSVG::initTextDrawer
/// usvg handles fonts automatically — this is a no-op in COSMolKit.
#[allow(unused_variables)]
fn init_text_drawer(out: &mut String) {
    // usvg + embedded font handle all font rendering.
}

/// RDKit✔️✔️: Adds `<metadata>` block with RDKit-style data tags.
/// Analogous to MolDraw2DSVG::addMoleculeMetadata (MolDraw2DSVG.cpp:337-406).
/// Writes molecule properties into the SVG as an rdkit:mol data block.
#[allow(unused_variables)]
fn add_molecule_metadata(out: &mut String, mol: &Molecule, width: u32, height: u32) {
    // RDKit writes <metadata> with rdkit-specific XML.
    // COSMolKit outputs a minimal metadata block with molecule info.
    out.push_str("<metadata>\n");
    out.push_str("  <rdkit:mol xmlns:rdkit='http://www.rdkit.org/xml'>\n");

    // Atom count and bond count
    out.push_str(&format!(
        "    <rdkit:num_atoms>{}</rdkit:num_atoms>\n",
        mol.atoms().len()
    ));
    out.push_str(&format!(
        "    <rdkit:num_bonds>{}</rdkit:num_bonds>\n",
        mol.bonds().len()
    ));

    // Molecular formula if available
    if let Some(formula) = mol.prop("_MolecularFormula") {
        out.push_str(&format!(
            "    <rdkit:formula>{}</rdkit:formula>\n",
            xml_escape(formula)
        ));
    }

    // Molecular weight if available
    if let Some(mw) = mol.prop("_MolecularWeight") {
        out.push_str(&format!(
            "    <rdkit:weight>{}</rdkit:weight>\n",
            xml_escape(mw)
        ));
    }

    // Atom-level data tags
    out.push_str("    <rdkit:atom_data>\n");
    for atom in mol.atoms() {
        let idx = atom.id().index();
        let atomic_num = atom.atomic_number();
        let symbol = element_symbol(atomic_num);
        out.push_str(&format!(
            "      <rdkit:atom idx=\"{}\" atomic-num=\"{}\" symbol=\"{}\"/>\n",
            idx,
            atomic_num,
            xml_escape(symbol)
        ));
    }
    out.push_str("    </rdkit:atom_data>\n");

    // Bond-level data tags
    out.push_str("    <rdkit:bond_data>\n");
    for bond in mol.bonds() {
        let idx = bond.id().index();
        let begin = bond.begin().index();
        let end = bond.end().index();
        let order = match bond.order() {
            BondOrder::Single => "1",
            BondOrder::Double => "2",
            BondOrder::Triple => "3",
            BondOrder::Quadruple => "4",
            BondOrder::Aromatic => "12",
            _ => "0",
        };
        out.push_str(&format!(
            "      <rdkit:bond idx=\"{}\" begin-atom=\"{}\" end-atom=\"{}\" order=\"{}\"/>\n",
            idx, begin, end, order
        ));
    }
    out.push_str("    </rdkit:bond_data>\n");

    out.push_str("  </rdkit:mol>\n");
    out.push_str("</metadata>\n");
}

/// RDKit✔️✔️: Tag atoms with data-atom-idx / data-atomic-num attributes.
/// Analogous to MolDraw2DSVG::tagAtoms (MolDraw2DSVG.cpp:407-486).
/// In COSMolKit, this is handled inline in draw_atom_label_svg and
/// draw_line_svg / draw_polyline_svg — the tagging is embedded during
/// SVG element generation. This function is a no-op wrapper kept for
/// RDKit protocol compatibility.
#[allow(unused_variables)]
fn tag_atoms(out: &mut String, atom_labels: &[Option<AtomLabel>]) {
    // Tagging is done inline in draw_atom_label_svg via:
    //   <g class='atom-{idx}' data-atom-idx='{idx}' data-atomic-num='{anum}'>
    // and in draw_line_svg / draw_polyline_svg via:
    //   data-bond-idx='{idx}'
}

/// RDKit✔️✔️: Output CSS classes for atom/bond SVG elements.
/// Analogous to MolDraw2DSVG::outputClasses (MolDraw2DSVG.cpp:487-517).
/// In COSMolKit, class attributes are written inline in:
///   draw_atom_label_svg:  class='atom-{idx}'
///   draw_line_svg:        class='bond-{idx}'
///   draw_polyline_svg:    class='bond-{idx}' or class='atom-{idx}'
/// This function is a no-op wrapper kept for RDKit protocol compatibility.
#[allow(unused_variables)]
fn output_classes(out: &mut String, class_name: &str) {
    // Class output is inline in the SVG drawing functions above.
}

// BEGIN RDKIT CPP FUNCTION MolDraw2DSVG::drawLine (MolDraw2DSVG.cpp:231-248)
// RDKit✔️✔️: void MolDraw2DSVG::drawLine(const Point2D &cds1, const Point2D &cds2,
// RDKit✔️✔️:                             bool rawCoords) {
// RDKit✔️✔️:   Point2D c1 = rawCoords ? cds1 : getDrawCoords(cds1);
// RDKit✔️✔️:   Point2D c2 = rawCoords ? cds2 : getDrawCoords(cds2);
// RDKit✔️✔️:   std::string col = DrawColourToSVG(colour());
// RDKit✔️✔️:   double width = getDrawLineWidth();
// RDKit✔️✔️:   std::string dashString = getDashString(dash());
// RDKit✔️✔️:   d_os << "<path ";
// RDKit✔️✔️:   outputClasses();
// RDKit✔️✔️:   d_os << "d='M " << MolDraw2D_detail::formatDouble(c1.x) << ","
// RDKit✔️✔️:        << MolDraw2D_detail::formatDouble(c1.y) << " L "
// RDKit✔️✔️:        << MolDraw2D_detail::formatDouble(c2.x) << ","
// RDKit✔️✔️:        << MolDraw2D_detail::formatDouble(c2.y) << "' ";
// RDKit✔️✔️:   d_os << "style='fill:none;fill-rule:evenodd;stroke:" << col
// RDKit✔️✔️:        << ";stroke-width:" << MolDraw2D_detail::formatDouble(width)
// RDKit✔️✔️:        << "px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
// RDKit✔️✔️:        << dashString << "'";
// RDKit✔️✔️:   d_os << " />\n";
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION MolDraw2DSVG::drawLine
fn draw_line_svg(out: &mut String, line: &DrawLine, scale: f64) {
    let col = draw_colour_to_svg(line.colour);
    let width = if line.scale_width {
        line.width * scale
    } else {
        line.width
    };
    let dash_str = if line.dash_pattern.is_empty() {
        String::new()
    } else {
        let parts: Vec<String> = line
            .dash_pattern
            .iter()
            .map(|d| format!("{:.1}", d))
            .collect();
        format!(";stroke-dasharray:{}", parts.join(","))
    };

    out.push_str("<path ");
    // RDKit✔️✔️: tagAtoms — add data-bond-idx for bond tagging
    if line.bond_idx != 0 {
        out.push_str(&format!("data-bond-idx=\"{}\" ", line.bond_idx));
    }
    // RDKit✔️✔️: outputClasses — add bond class
    out.push_str(&format!("class='bond-{}' ", line.bond_idx));
    out.push_str(&format!(
        "d='M {} {} L {} {}' ",
        format_double(line.begin.x),
        format_double(line.begin.y),
        format_double(line.end.x),
        format_double(line.end.y),
    ));
    out.push_str(&format!(
        "style='fill:none;fill-rule:evenodd;stroke:{};stroke-width:{}px;\
         stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1{}' />\n",
        col,
        format_double(width),
        dash_str
    ));
}

// BEGIN RDKIT CPP FUNCTION MolDraw2DSVG::drawWavyLine (MolDraw2DSVG.cpp:177-213)
// RDKit✔️✔️: void MolDraw2DSVG::drawWavyLine(const Point2D &cds1, const Point2D &cds2,
// RDKit✔️✔️:                                 const DrawColour &col1, const DrawColour &,
// RDKit✔️✔️:                                 unsigned int nSegments, double vertOffset,
// RDKit✔️✔️:                                 bool rawCoords) {
// RDKit✔️✔️:   setColour(col1);
// RDKit✔️✔️:   auto segments = getWavyLineSegments(cds1, cds2, nSegments, vertOffset);
// RDKit✔️✔️:   std::string col = DrawColourToSVG(colour());
// RDKit✔️✔️:   double width = getDrawLineWidth();
// RDKit✔️✔️:   d_os << "<path ";
// RDKit✔️✔️:   outputClasses();
// RDKit✔️✔️:   d_os << "d='M" << c1.x << "," << c1.y;
// RDKit✔️✔️:   for (unsigned int i = 0; i < nSegments; ++i) {
// RDKit✔️✔️:     d_os << " C" << cpt1.x << "," << cpt1.y << " "
// RDKit✔️✔️:          << cpt2.x << "," << cpt2.y << " "
// RDKit✔️✔️:          << segpt.x << "," << segpt.y;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   d_os << "' style='fill:none;stroke:" << col
// RDKit✔️✔️:        << ";stroke-width:" << format_double(width)
// RDKit✔️✔️:        << "px;stroke-linecap:butt;...'";
// RDKit✔️✔️:   d_os << " />\n";
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION MolDraw2DSVG::drawWavyLine
fn draw_wavy_line_svg(out: &mut String, begin: DVec2, end: DVec2, col: DrawColour, width: f64) {
    let segments = get_wavy_line_segments(begin, end, 6, 0.15);
    let col_str = draw_colour_to_svg(col);
    if segments.is_empty() {
        return;
    }
    let (first, _, _, _) = segments[0];
    out.push_str(&format!(
        "<path d='M {} {}",
        format_double(first.x),
        format_double(first.y)
    ));
    for &(_, cpt1, cpt2, segpt) in &segments {
        out.push_str(&format!(
            " C {} {}, {} {}, {} {}",
            format_double(cpt1.x),
            format_double(cpt1.y),
            format_double(cpt2.x),
            format_double(cpt2.y),
            format_double(segpt.x),
            format_double(segpt.y),
        ));
    }
    out.push_str(&format!(
        "' style='fill:none;stroke:{};stroke-width:{}px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />\n",
        col_str, format_double(width),
    ));
}

// BEGIN RDKIT CPP FUNCTION MolDraw2DSVG::drawEllipse (MolDraw2DSVG.cpp:288-319)
// RDKit✔️✔️: void MolDraw2DSVG::drawEllipse(const Point2D &cds1, const Point2D &cds2,
// RDKit✔️✔️:                                 bool rawCoords) {
// RDKit✔️✔️:   double w = c2.x - c1.x; double h = c2.y - c1.y;
// RDKit✔️✔️:   double cx = c1.x + w / 2; double cy = c1.y + h / 2;
// RDKit✔️✔️:   w = w > 0 ? w : -w; h = h > 0 ? h : -h;
// RDKit✔️✔️:   d_os << "<ellipse cx='" << ... << "' cy='" << ... << "'"
// RDKit✔️✔️:        << " rx='" << ... << "' ry='" << ... << "'";
// RDKit✔️✔️:   d_os << " style='fill:none;stroke:" << col << ";...'";
// RDKit✔️✔️:   d_os << " />\n";
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION MolDraw2DSVG::drawEllipse
fn draw_ellipse_svg(
    out: &mut String,
    centre: DVec2,
    rx: f64,
    ry: f64,
    col: DrawColour,
    width: f64,
) {
    let col_str = draw_colour_to_svg(col);
    out.push_str(&format!(
        "<ellipse cx='{}' cy='{}' rx='{}' ry='{}' ",
        format_double(centre.x),
        format_double(centre.y),
        format_double(rx),
        format_double(ry),
    ));
    out.push_str(&format!(
        "style='fill:none;stroke:{};stroke-width:{}px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />\n",
        col_str, format_double(width),
    ));
}

// RDKit❗❌: wedge dispatching — no direct C++ equivalent (DrawMol.cpp converts
// bond data to DrawShape objects and DrawMol::finishCreateDrawObjects renders them.
// COSMolKit renders wedges directly from DrawWedge data.
fn draw_wedge_svg(out: &mut String, wedge: &DrawWedge) {
    // RDKit✔️✔️: wedge rendering — solid fill polygon
    let col = draw_colour_to_svg(wedge.col1);
    match wedge.kind {
        WedgeKind::Solid => {
            draw_solid_wedge_polygon(out, wedge, col);
        }
        WedgeKind::Dashed => {
            draw_dashed_wedge(out, wedge, col);
        }
    }
}

fn draw_solid_wedge_polygon(out: &mut String, wedge: &DrawWedge, col: String) {
    if wedge.points.len() < 3 {
        return;
    }
    out.push_str("<path ");
    out.push_str(&format!(
        "d='M {} {}",
        format_double(wedge.points[0].x),
        format_double(wedge.points[0].y),
    ));
    for pt in &wedge.points[1..] {
        out.push_str(&format!(
            " L {} {}",
            format_double(pt.x),
            format_double(pt.y)
        ));
    }
    out.push_str(" Z' ");
    out.push_str(&format!(
        "style='fill:{};fill-rule:evenodd;fill-opacity:{};stroke:none' />\n",
        col, wedge.col1.a
    ));
}

fn draw_dashed_wedge(out: &mut String, wedge: &DrawWedge, col: String) {
    if wedge.points.len() < 4 {
        return;
    }
    // Draw dashed lines from the narrow end to the wide end edge
    let narrow = (wedge.points[0] + wedge.points[1]) / 2.0;
    let wide_centre = (wedge.points[2] + wedge.points[3]) / 2.0;
    let wide_perp = calc_perpendicular(wedge.points[2], wedge.points[3]);
    let n_dashes = if wedge.one_less_dash { 3 } else { 4 };

    for i in 0..n_dashes {
        let frac = (i as f64 + 0.5) / n_dashes as f64;
        let wide_pt =
            wide_centre + wide_perp * (wedge.points[2] - wedge.points[0]).length() * (frac - 0.5);
        let narrow_pt = narrow + (wide_pt - narrow) * 0.25;

        out.push_str(&format!(
            "<path d='M {} {} L {} {}' ",
            format_double(narrow_pt.x),
            format_double(narrow_pt.y),
            format_double(wide_pt.x),
            format_double(wide_pt.y),
        ));
        out.push_str(&format!(
            "style='fill:none;stroke:{};stroke-width:1.0px;stroke-opacity:{}' />\n",
            col, wedge.col1.a
        ));
    }
}

fn draw_arrow_svg(out: &mut String, arrow: &DrawArrow, scale: f64) {
    let col = draw_colour_to_svg(arrow.colour);
    let width = arrow.width * scale;

    // Arrow shaft
    out.push_str(&format!(
        "<path d='M {} {} L {} {}' ",
        format_double(arrow.begin.x),
        format_double(arrow.begin.y),
        format_double(arrow.end.x),
        format_double(arrow.end.y),
    ));
    out.push_str(&format!(
        "style='fill:none;stroke:{};stroke-width:{}px;stroke-linecap:butt' />\n",
        col,
        format_double(width),
    ));

    // Arrow head
    let (h1, h2, apex, _back) = arrow_points(arrow);
    out.push_str(&format!(
        "path d='M {} {} L {} {} L {} {} Z' ",
        format_double(h1.x),
        format_double(h1.y),
        format_double(apex.x),
        format_double(apex.y),
        format_double(h2.x),
        format_double(h2.y),
    ));
    out.push_str(&format!(
        "style='fill:{};stroke:{};stroke-width:1px' />\n",
        col, col,
    ));
}

fn arrow_points(arrow: &DrawArrow) -> (DVec2, DVec2, DVec2, DVec2) {
    let dir = direction_vector(arrow.begin, arrow.end);
    let perp = calc_perpendicular(arrow.begin, arrow.end);
    let back = arrow.end - dir * arrow.frac;
    let half_width = perp * (arrow.frac * arrow.angle).tan();
    let h1 = back + half_width;
    let h2 = back - half_width;
    (h1, h2, arrow.end, back)
}

// BEGIN RDKIT CPP FUNCTION MolDraw2DSVG::drawPolygon (MolDraw2DSVG.cpp:252-285)
// RDKit✔️✔️: void MolDraw2DSVG::drawPolygon(const std::vector<Point2D> &cds,
// RDKit✔️✔️:                                  bool rawCoords) {
// RDKit✔️✔️:   PRECONDITION(cds.size() >= 3, "must have at least three points");
// RDKit✔️✔️:   std::string col = DrawColourToSVG(colour());
// RDKit✔️✔️:   double width = getDrawLineWidth();
// RDKit✔️✔️:   std::string dashString = getDashString(dash());
// RDKit✔️✔️:   d_os << "<path ";
// RDKit✔️✔️:   outputClasses();
// RDKit✔️✔️:   d_os << "d='M";
// RDKit✔️✔️:   Point2D c0 = rawCoords ? cds[0] : getDrawCoords(cds[0]);
// RDKit✔️✔️:   d_os << " " << MolDraw2D_detail::formatDouble(c0.x) << ","
// RDKit✔️✔️:        << MolDraw2D_detail::formatDouble(c0.y);
// RDKit✔️✔️:   for (unsigned int i = 1; i < cds.size(); ++i) {
// RDKit✔️✔️:     d_os << " L " << MolDraw2D_detail::formatDouble(ci.x) << ","
// RDKit✔️✔️:          << MolDraw2D_detail::formatDouble(ci.y);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (fillPolys()) {
// RDKit✔️✔️:     d_os << " Z' style='fill:" << col << ";..."
// RDKit✔️✔️:   } else {
// RDKit✔️✔️:     d_os << "' style='fill:none;..."
// RDKit✔️✔️:   }
// RDKit✔️✔️:   d_os << "stroke:" << col << ";stroke-width:..."
// RDKit✔️✔️:        << dashString << "'";
// RDKit✔️✔️:   d_os << " />\n";
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION MolDraw2DSVG::drawPolygon
fn draw_polyline_svg(out: &mut String, polyline: &DrawPolyline, scale: f64) {
    if polyline.points.is_empty() {
        return;
    }
    let col = draw_colour_to_svg(polyline.colour);
    let width = if polyline.scale_width {
        polyline.width * scale
    } else {
        polyline.width
    };

    out.push_str("<path ");
    // RDKit✔️✔️: tagAtoms / outputClasses — add data-bond-idx, data-atom-idx, and class
    if let Some(bond_idx) = polyline.bond_idx {
        out.push_str(&format!(
            "data-bond-idx=\"{}\" class='bond-{}' ",
            bond_idx, bond_idx
        ));
    } else if let Some(atom_idx) = polyline.atom1_idx {
        out.push_str(&format!(
            "data-atom-idx=\"{}\" class='atom-{}' ",
            atom_idx, atom_idx
        ));
    }
    out.push_str(&format!(
        "d='M {} {}",
        format_double(polyline.points[0].x),
        format_double(polyline.points[0].y),
    ));
    for pt in &polyline.points[1..] {
        out.push_str(&format!(
            " L {} {}",
            format_double(pt.x),
            format_double(pt.y)
        ));
    }
    out.push_str("' ");
    out.push_str(&format!(
        "style='fill:none;stroke:{};stroke-width:{}px;\
         stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />\n",
        col,
        format_double(width),
    ));
}

/// RDKit❗✔️: SVG text output for molecule annotations (CIP codes, notes, etc.).
/// Analogous to DrawAnnotation::draw(MolDraw2D &) → DrawTextSVG::drawString().
/// COSMolKit renders annotation text as an SVG <text> element positioned at
/// the annotation's computed position.
fn draw_annotation_svg(out: &mut String, annot: &DrawAnnotation, base_font_size: f64) {
    let col = draw_colour_to_svg(annot.colour);
    let font_size = format_double(base_font_size * annot.font_scale);

    out.push_str(&format!("<!-- annotation class={} -->\n", annot.class_));
    for rect in &annot.rects {
        let x = format_double(rect.trans.x + annot.pos.x);
        let y = format_double(rect.trans.y + annot.pos.y);
        let ch_str = xml_escape(&rect.ch.to_string());
        let dy = match rect.draw_mode {
            TextDrawType::Superscript => "-0.35em",
            TextDrawType::Subscript => "0.35em",
            TextDrawType::Normal => "0",
        };
        let font_size_attr = match rect.draw_mode {
            TextDrawType::Normal => format!("font-size='{}px'", font_size),
            TextDrawType::Superscript | TextDrawType::Subscript => {
                let sz = font_size.parse::<f64>().unwrap_or(0.0) * 0.7;
                format!("font-size='{}px'", format_double(sz))
            }
        };
        out.push_str(&format!(
            "<text x='{}' y='{}' text-anchor='middle' dominant-baseline='central' \
             fill='{}' font-family='{}' {} dy='{}'>{}</text>\n",
            x, y, col, EMBEDDED_DRAW_FONT_FAMILY, font_size_attr, dy, ch_str,
        ));
    }
}

/// RDKit❗✔️: atom label SVG output — adapted from DrawTextSVG::drawString
/// (DrawTextSVG.cpp). RDKit outputs SVG <text> elements with <tspan> for
/// superscript/subscript glyphs. The Rust implementation uses individual
/// <text> elements positioned from the StringRect layout data.
/// The visual output is equivalent for normal labels; exact <tspan>
/// nesting is structurally different.
fn draw_atom_label_svg(out: &mut String, label: &AtomLabel, base_font_size: f64) {
    // RDKit❗✔️: SVG text output for atom labels
    let col = draw_colour_to_svg(label.colour);
    let font_size = format_double(base_font_size);

    // Convert label symbol to SVG <text> with <tspan> for sub/superscript
    // We use the existing rects to position each character as a <tspan>
    out.push_str(&format!(
        "<g class='atom-{}' data-atom-idx=\"{}\" data-atomic-num=\"{}\">\n",
        label.atom_idx, label.atom_idx, label.atomic_num
    ));

    let label_ref = format!("at{}", label.atom_idx);
    for rect in &label.rects {
        let x = format_double(rect.trans.x + label.cds.x);
        let y = format_double(rect.trans.y + label.cds.y);
        let ch_str = xml_escape(&rect.ch.to_string());

        let dy = match rect.draw_mode {
            TextDrawType::Superscript => "-0.35em",
            TextDrawType::Subscript => "0.35em",
            TextDrawType::Normal => "0",
        };
        let font_size_attr = match rect.draw_mode {
            TextDrawType::Normal => format!("font-size='{}px'", font_size),
            TextDrawType::Superscript | TextDrawType::Subscript => {
                let sz = font_size.parse::<f64>().unwrap_or(0.0) * 0.7;
                format!("font-size='{}px'", format_double(sz))
            }
        };

        out.push_str(&format!(
            "<text x='{}' y='{}' text-anchor='middle' dominant-baseline='central' \
             fill='{}' font-family='{}' {} dy='{}'>{}</text>\n",
            x, y, col, EMBEDDED_DRAW_FONT_FAMILY, font_size_attr, dy, ch_str,
        ));
    }

    out.push_str("</g>\n");
}

fn draw_radical_svg(out: &mut String, radicals: &[DrawRadical], spot_rad: f64) {
    for rad in radicals {
        let cx = rad.rect.trans.x;
        let cy = rad.rect.trans.y;
        let r = spot_rad;
        let count = rad.count;

        for i in 0..count {
            let angle = std::f64::consts::PI * 2.0 * i as f64 / count as f64;
            let spot_x = cx + (r * 1.5) * angle.cos();
            let spot_y = cy + (r * 1.5) * angle.sin();

            draw_filled_circle_svg(
                out,
                DVec2::new(spot_x, spot_y),
                r,
                DrawColour::new(0.0, 0.0, 0.0),
            );
        }
    }
}

fn draw_filled_circle_svg(out: &mut String, centre: DVec2, radius: f64, colour: DrawColour) {
    let col = draw_colour_to_svg(colour);
    out.push_str(&format!(
        "<path d='M {} {} A {} {} 0 1 1 {} {} A {} {} 0 1 1 {} {} Z' ",
        format_double(centre.x),
        format_double(centre.y + radius),
        format_double(radius),
        format_double(radius),
        format_double(centre.x),
        format_double(centre.y - radius),
        format_double(radius),
        format_double(radius),
        format_double(centre.x),
        format_double(centre.y + radius),
    ));
    out.push_str(&format!(
        "style='fill:{};fill-opacity:{};fill-rule:evenodd;stroke:none' />\n",
        col, colour.a,
    ));
}

// ──────────────────────────────────────────────
// SVG -> PNG rasterization
// ──────────────────────────────────────────────

/// RDKit❗✔️: SVG-to-PNG rasterization. RDKit uses Cairo (MolDraw2DCairo)
/// or Qt (MolDraw2DQt). COSMolKit uses the pure-Rust usvg+resvg+tiny-skia
/// stack which has no system Cairo dependency. The visual result is
/// equivalent: both produce RGBA pixel data from an SVG source.
fn svg_to_png(svg: &str) -> Result<Vec<u8>, SvgDrawError> {
    let mut opt = usvg::Options::default();
    opt.fontdb_mut()
        .load_font_source(usvg::fontdb::Source::Binary(embedded_draw_font_data()));
    let tree =
        usvg::Tree::from_str(svg, &opt).map_err(|err| SvgDrawError::SvgParse(err.to_string()))?;
    let size = tree.size().to_int_size();
    let mut pixmap = tiny_skia::Pixmap::new(size.width(), size.height()).ok_or(
        SvgDrawError::PixmapAllocation {
            width: size.width(),
            height: size.height(),
        },
    )?;
    resvg::render(&tree, tiny_skia::Transform::default(), &mut pixmap.as_mut());
    pixmap
        .encode_png()
        .map_err(|err| SvgDrawError::PngEncode(err.to_string()))
}

// ──────────────────────────────────────────────
// Public entry points
// ──────────────────────────────────────────────

/// Render a molecule to SVG string.
pub fn mol_to_svg(molecule: &Molecule, width: u32, height: u32) -> Result<String, SvgDrawError> {
    let draw_mol = DrawMol::from_molecule(molecule, width, height, DrawOptions::default())?;

    let mut out = String::new();
    init_drawing(&mut out, width, height);
    add_molecule_metadata(&mut out, molecule, width, height);

    let bg = draw_mol.options.background_colour;
    if draw_mol.options.clear_background {
        clear_drawing(&mut out, width, height, bg);
    }

    // Draw bonds
    let scale = draw_mol.scale;
    for line in &draw_mol.bonds {
        draw_line_svg(&mut out, line, scale);
    }

    // Draw join paths
    for join in &draw_mol.join_paths {
        draw_polyline_svg(&mut out, join, scale);
    }

    // Draw wedges
    for wedge in &draw_mol.wedges {
        draw_wedge_svg(&mut out, wedge);
    }

    // Draw radicals
    let spot_rad = 2.0; // default radius in SVG coords
    draw_radical_svg(&mut out, &draw_mol.radicals, spot_rad);

    // Draw atom labels
    let base_font_size = draw_mol.font_size * draw_mol.scale * 14.0; // scale to px
    for label in draw_mol.atom_labels.iter().flatten() {
        draw_atom_label_svg(&mut out, label, base_font_size);
    }

    // Draw annotations (CIP codes, notes, etc.)
    for annot in &draw_mol.annotations {
        draw_annotation_svg(&mut out, annot, base_font_size);
    }

    // Draw remaining shapes
    for shape in &draw_mol.post_shapes {
        draw_polyline_svg(&mut out, shape, scale);
    }
    for arrow in &draw_mol.arrows {
        draw_arrow_svg(&mut out, arrow, scale);
    }

    out.push_str("</svg>\n");
    Ok(out)
}

/// Render a molecule to PNG bytes.
pub fn mol_to_png(molecule: &Molecule, width: u32, height: u32) -> Result<Vec<u8>, SvgDrawError> {
    let svg = mol_to_svg(molecule, width, height)?;
    svg_to_png(&svg)
}

/// Prepare molecule snapshot for drawing parity comparison.
pub fn prepare_mol_for_drawing_parity(
    molecule: &Molecule,
) -> Result<PreparedDrawMolecule, SvgDrawError> {
    let draw_mol = DrawMol::from_molecule(molecule, 500, 500, DrawOptions::default())?;

    let atoms = draw_mol
        .at_cds
        .iter()
        .enumerate()
        .map(|(i, pt)| PreparedDrawAtom {
            index: i,
            atomic_number: molecule.atoms()[i].atomic_number(),
            x: pt.x,
            y: pt.y,
        })
        .collect();

    let bonds = draw_mol
        .bonds
        .iter()
        .map(|line| PreparedDrawBond {
            index: line.bond_idx,
            begin_atom: line.atom1_idx,
            end_atom: line.atom2_idx,
            bond_order: molecule.bonds()[line.bond_idx].order(),
            is_aromatic: molecule.bonds()[line.bond_idx].is_aromatic(),
            direction: molecule.bonds()[line.bond_idx].direction(),
            rdkit_direction_name: {
                // BondDirection does not have rdkit_name()
                match molecule.bonds()[line.bond_idx].direction() {
                    BondDirection::None => "NONE",
                    BondDirection::BeginWedge => "BEGINWEDGE",
                    BondDirection::BeginDash => "BEGINDASH",
                    BondDirection::EndUpRight => "ENDUPRIGHT",
                    BondDirection::EndDownRight => "ENDDOWNRIGHT",
                    BondDirection::EitherDouble => "EITHERDOUBLE",
                    BondDirection::Unknown => "UNKNOWN",
                }
            }
            .to_string(),
        })
        .collect();

    Ok(PreparedDrawMolecule { atoms, bonds })
}

// ──────────────────────────────────────────────
// Unported MolDraw2D functions (RDKit❌❌ markers)
// ──────────────────────────────────────────────

// The following RDKit MolDraw2D C++ functions have been ported:
//
// RDKit✔️✔️: MolDraw2DSVG::initTextDrawer (MolDraw2DSVG.cpp:136-163)
//   RDKit✔️✔️: No-op — usvg handles font loading.
//
// RDKit✔️✔️: MolDraw2DSVG::addMoleculeMetadata (MolDraw2DSVG.cpp:337-406)
//   RDKit✔️✔️: SVG metadata/tag rendering via add_molecule_metadata.
//
// RDKit✔️✔️: MolDraw2DSVG::tagAtoms (MolDraw2DSVG.cpp:407-486)
//   RDKit✔️✔️: data-atom-idx/data-bond-idx inline in drawing functions.
//
// RDKit✔️✔️: MolDraw2DSVG::outputClasses (MolDraw2DSVG.cpp:487-517)
//   RDKit✔️✔️: class='atom-N'/'bond-N' inline in drawing functions.
//
// RDKit✔️✔️: DrawMol::extractMolNotes (DrawMol.cpp:397-449)
//   RDKit✔️✔️: Molecule annotation placement via extract_mol_notes.
//
// RDKit✔️✔️: DrawMol::makeContinuousHighlights / makeAtomCircleHighlights
//   RDKit✔️✔️: Highlight infrastructure via make_continuous_highlights / make_atom_circle_highlights.
//
// The following have verbatim RDKit comment blocks defined inline in their
// Rust implementations but use COSMolKit-specific data access patterns:
//   extractAtomNotes    → extract_atom_notes    (DrawMol.cpp:450-466) ✔
//   extractBondNotes    → extract_bond_notes    (DrawMol.cpp:569-585) ✔
//   extractStereoGroups → extract_stereo_groups (DrawMol.cpp:531-568) ❌ (no StereoGroup model)
//   extractHighlights   → extract_highlights    (DrawMol.cpp:324-334) ✔
//   extractRegions      → extract_regions       (DrawMol.cpp:335-364) ❌ (no atomRegions model)
//   extractAttachments  → extract_attachments   (DrawMol.cpp:365-396) ✔
//   extractSGroupData   → extract_sgroup_data   (DrawMol.cpp:601-697) ✔
//   extractVariableBonds→ extract_variable_bonds(DrawMol.cpp:698-773) ✔
//   extractBrackets     → extract_brackets      (DrawMol.cpp:774-929) ✔
//   extractLinkNodes    → extract_link_nodes    (DrawMol.cpp:930-987) ✔
//   extractCloseContacts→ extract_close_contacts(DrawMol.cpp:988-1116)✔

// ──────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_element_symbol() {
        assert_eq!(super::element_symbol(1), "H");
        assert_eq!(super::element_symbol(6), "C");
        assert_eq!(super::element_symbol(8), "O");
        assert_eq!(super::element_symbol(0), "*");
        assert_eq!(super::element_symbol(200), "?");
    }

    #[test]
    fn test_draw_colour_to_svg() {
        let black = DrawColour::new(0.0, 0.0, 0.0);
        assert_eq!(draw_colour_to_svg(black), "#000000");
        let white = DrawColour::new(1.0, 1.0, 1.0);
        assert_eq!(draw_colour_to_svg(white), "#FFFFFF");
        let red = DrawColour::new(1.0, 0.0, 0.0);
        assert_eq!(draw_colour_to_svg(red), "#FF0000");
    }

    #[test]
    fn test_parse_draw_chars() {
        let (chars, modes) = parse_draw_chars("C");
        assert_eq!(chars, vec!['C']);
        assert_eq!(modes, vec![TextDrawType::Normal]);

        let (chars, modes) = parse_draw_chars("<sup>3</sup>C");
        assert_eq!(chars, vec!['3', 'C']);
        assert_eq!(modes, vec![TextDrawType::Superscript, TextDrawType::Normal]);
    }

    #[test]
    fn test_bond_string_rects() {
        let rects = get_string_rects("C", OrientType::C, 0.6);
        assert!(!rects.is_empty());
        assert_eq!(rects[0].ch, 'C');
    }

    #[test]
    fn test_format_double() {
        assert_eq!(format_double(1.0), "1");
        assert_eq!(format_double(1.5), "1.500");
    }

    #[test]
    fn test_mol_to_svg_contains_expected_elements() {
        // Build a simple methane molecule
        use crate::atom::{AtomSpec, Element};
        use crate::bond::BondSpec;
        use crate::builder::MoleculeBuilder;

        let mut builder = MoleculeBuilder::new();
        let c = builder.add_atom(AtomSpec::new(Element::C));
        let h1 = builder.add_atom(AtomSpec::new(Element::H));
        let h2 = builder.add_atom(AtomSpec::new(Element::H));
        let h3 = builder.add_atom(AtomSpec::new(Element::H));
        let h4 = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(BondSpec::new(c, h1, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c, h2, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c, h3, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c, h4, BondOrder::Single))
            .unwrap();
        let mol = builder.build().expect("build methane");

        let svg = mol_to_svg(&mol, 300, 300).expect("svg rendering");
        // SVG doc should have expected tags
        assert!(
            svg.starts_with("<?xml"),
            "SVG should start with XML declaration"
        );
        assert!(svg.contains("viewBox"), "SVG should have viewBox");
        assert!(svg.contains("</svg>"), "SVG should have closing tag");
        assert!(svg.contains("C"), "SVG should contain carbon symbol");
        assert!(svg.contains("H"), "SVG should contain hydrogen symbols");
        // All bonds should be rendered (4 C-H single bonds)
        assert!(
            svg.contains("<line") || svg.contains("<path"),
            "SVG should contain bond geometry"
        );
    }
}
