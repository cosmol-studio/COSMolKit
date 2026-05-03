//! RDKit MolDraw2D-compatible drawing entry points.

use crate::distgeom::{rdkit_atom_rings_for_chiral_hs, rdkit_bond_rings_preserve_order};
use crate::io::molblock::{
    determine_bond_wedge_state_rdkit_subset, pick_bonds_to_wedge_rdkit_subset,
};
use crate::kekulize::kekulize_in_place;
use crate::{
    Atom, Bond, BondDirection, BondOrder, BondStereo, ChiralTag, Molecule, ValenceModel,
    assign_valence,
};
use glam::DVec2;
use std::collections::HashMap;
use std::sync::{Arc, OnceLock};

const EMBEDDED_DRAW_FONT_BYTES: &[u8] = include_bytes!("../assets/fonts/NotoSans-Regular.ttf");
const EMBEDDED_DRAW_FONT_FAMILY: &str = "Noto Sans";

fn embedded_draw_font_data() -> Arc<dyn AsRef<[u8]> + Send + Sync> {
    static FONT_DATA: OnceLock<Arc<Vec<u8>>> = OnceLock::new();
    FONT_DATA
        .get_or_init(|| Arc::new(EMBEDDED_DRAW_FONT_BYTES.to_vec()))
        .clone() as Arc<dyn AsRef<[u8]> + Send + Sync>
}

/// Errors returned by SVG drawing routines.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SvgDrawError {
    /// 2D coordinate generation failed while preparing the drawing molecule.
    #[error("RDKit MolDraw2D prepareMolForDrawing coordinate generation failed: {0}")]
    CoordinateGeneration(String),
    /// The requested drawing path has not yet been ported from RDKit MolDraw2D.
    #[error("RDKit MolDraw2D SVG generation path is not implemented yet: {0}")]
    Unsupported(String),
    /// SVG parsing failed before rasterization.
    #[error("SVG rasterization parse failed: {0}")]
    SvgParse(String),
    /// PNG pixmap allocation failed.
    #[error("PNG rasterization pixmap allocation failed for {width}x{height}")]
    PixmapAllocation { width: u32, height: u32 },
    /// PNG encoding failed.
    #[error("PNG encoding failed: {0}")]
    PngEncode(String),
}

/// Atom snapshot after the RDKit `PrepareMolForDrawing()` preparation boundary.
#[derive(Debug, Clone, PartialEq)]
pub struct PreparedDrawAtom {
    pub index: usize,
    pub atomic_num: u8,
    pub x: f64,
    pub y: f64,
}

/// Bond snapshot after the RDKit `PrepareMolForDrawing()` preparation boundary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PreparedDrawBond {
    pub index: usize,
    pub begin_atom: usize,
    pub end_atom: usize,
    pub bond_type: BondOrder,
    pub is_aromatic: bool,
    pub direction: BondDirection,
    pub rdkit_direction_name: String,
}

/// Molecule snapshot after the RDKit `PrepareMolForDrawing()` preparation boundary.
#[derive(Debug, Clone, PartialEq)]
pub struct PreparedDrawMolecule {
    pub atoms: Vec<PreparedDrawAtom>,
    pub bonds: Vec<PreparedDrawBond>,
}

#[derive(Debug, Copy, Clone, PartialEq)]
struct DrawColour {
    r: f64,
    g: f64,
    b: f64,
    a: f64,
}

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
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum OrientType {
    C,
    N,
    E,
    S,
    W,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum TextAlignType {
    Middle,
    Start,
    End,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum TextDrawType {
    Normal,
    Subscript,
    Superscript,
}

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
        self.dash_pattern == DashPattern::NoDash
    }
}

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

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum WedgeKind {
    Solid,
    Dashed,
}

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
    bond_idx: usize,
}

#[derive(Debug, Clone)]
struct DrawRadical {
    rect: StringRect,
    orient: OrientType,
    atom_idx: usize,
    count: u8,
}

#[derive(Debug, Copy, Clone)]
enum DrawItem {
    Line(usize),
    Wedge(usize),
    Arrow(usize),
}

#[derive(Debug, Clone)]
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
    implicit_hs: Vec<u8>,
    bonds: Vec<DrawLine>,
    wedges: Vec<DrawWedge>,
    arrows: Vec<DrawArrow>,
    draw_items: Vec<DrawItem>,
    join_paths: Vec<DrawPolyline>,
    post_shapes: Vec<DrawPolyline>,
    radicals: Vec<DrawRadical>,
    single_bond_lines: Vec<usize>,
    mean_bond_length: f64,
    font_size: f64,
    options: DrawOptions,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum DashPattern {
    NoDash,
    ShortDashes,
}

impl Default for DrawOptions {
    fn default() -> Self {
        Self {
            padding: 0.05,
            multiple_bond_offset: 0.15,
            bond_line_width: 2.0,
            scale_bond_width: false,
            clear_background: true,
            background_colour: DrawColour {
                r: 1.0,
                g: 1.0,
                b: 1.0,
                a: 1.0,
            },
            query_colour: DrawColour {
                r: 0.5,
                g: 0.5,
                b: 0.5,
                a: 1.0,
            },
            flag_close_contacts_dist: 3,
        }
    }
}

/// Serialize a molecule to RDKit-style SVG.
pub fn mol_to_svg(mol: &Molecule, width: u32, height: u32) -> Result<String, SvgDrawError> {
    let options = DrawOptions::default();
    let draw_mol = DrawMol::from_molecule(mol, width, height, options.clone())?;
    let mut svg = String::new();
    init_drawing(&mut svg, width, height);
    if options.clear_background {
        clear_drawing(&mut svg, width, height, options.background_colour);
    }
    for item in &draw_mol.draw_items {
        match *item {
            DrawItem::Line(idx) => draw_line(&mut svg, &draw_mol.bonds[idx], draw_mol.scale),
            DrawItem::Wedge(idx) => draw_wedge(&mut svg, &draw_mol.wedges[idx]),
            DrawItem::Arrow(idx) => draw_arrow(&mut svg, &draw_mol.arrows[idx]),
        }
    }
    for polyline in &draw_mol.join_paths {
        draw_polyline(&mut svg, polyline, draw_mol.scale);
    }
    let final_font_size = draw_mol.font_size * draw_mol.scale;
    for label in draw_mol.atom_labels.iter().flatten() {
        draw_atom_label(&mut svg, label, final_font_size);
    }
    draw_radicals(
        &mut svg,
        &draw_mol.radicals,
        draw_mol.radical_spot_radius_unscaled() * draw_mol.scale,
    );
    for polyline in &draw_mol.post_shapes {
        draw_polyline(&mut svg, polyline, draw_mol.scale);
    }
    svg.push_str("</svg>\n");
    Ok(svg)
}

/// Rasterize the RDKit-style SVG drawing into PNG bytes.
pub fn mol_to_png(mol: &Molecule, width: u32, height: u32) -> Result<Vec<u8>, SvgDrawError> {
    let svg = mol_to_svg(mol, width, height)?;
    svg_to_png(&svg)
}

/// Return the prepared molecule state matching RDKit
/// `PrepareMolForDrawing(kekulize=true, addChiralHs=true, wedgeBonds=true, forceCoords=true)`.
pub fn prepare_mol_for_drawing_parity(
    mol: &Molecule,
) -> Result<PreparedDrawMolecule, SvgDrawError> {
    let mut draw_mol = prepare_working_mol_for_drawing(mol)?;
    if draw_mol.coords_2d().is_none() {
        draw_mol.compute_2d_coords().map_err(|err| {
            SvgDrawError::CoordinateGeneration(format!(
                "compute2DCoords(canonOrient=true) failed: {err}"
            ))
        })?;
    }
    let coords = draw_mol
        .coords_2d()
        .expect("coords were generated immediately above");
    let wedge_bonds = pick_bonds_to_wedge_rdkit_subset(&draw_mol);

    let atoms = draw_mol
        .atoms()
        .iter()
        .map(|atom| {
            let pt = coords[atom.index];
            PreparedDrawAtom {
                index: atom.index,
                atomic_num: atom.atomic_num,
                x: pt.x,
                y: pt.y,
            }
        })
        .collect();
    let bonds = draw_mol
        .bonds()
        .iter()
        .map(|bond| {
            let direction =
                determine_bond_wedge_state_rdkit_subset(bond, &wedge_bonds, coords, &draw_mol);
            let mut begin_atom = bond.begin_atom;
            let mut end_atom = bond.end_atom;
            let is_wedged_chiral_bond = matches!(
                direction,
                BondDirection::EndUpRight | BondDirection::EndDownRight
            ) && wedge_bonds.contains_key(&bond.index);
            if is_wedged_chiral_bond {
                if let Some(&from_atom) = wedge_bonds.get(&bond.index) {
                    if from_atom != begin_atom {
                        std::mem::swap(&mut begin_atom, &mut end_atom);
                    }
                }
            }
            let rdkit_direction_name = match (direction, is_wedged_chiral_bond) {
                (BondDirection::None, _) => "NONE",
                (BondDirection::EndUpRight, true) => "BEGINWEDGE",
                (BondDirection::EndDownRight, true) => "BEGINDASH",
                (BondDirection::EndUpRight, false) => "ENDUPRIGHT",
                (BondDirection::EndDownRight, false) => "ENDDOWNRIGHT",
            }
            .to_string();
            PreparedDrawBond {
                index: bond.index,
                begin_atom,
                end_atom,
                bond_type: bond.order,
                is_aromatic: bond.is_aromatic,
                direction,
                rdkit_direction_name,
            }
        })
        .collect();

    Ok(PreparedDrawMolecule { atoms, bonds })
}

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

fn prepare_working_mol_for_drawing(mol: &Molecule) -> Result<Molecule, SvgDrawError> {
    let mut draw_mol = mol.clone();
    // Mirrors MolDraw2DUtils::prepareMolForDrawing(kekulize=true):
    // KekulizeIfPossible(..., clearAromaticFlags=false) suppresses
    // kekulization failures and leaves the molecule unchanged.
    let _ = kekulize_in_place(&mut draw_mol, false);
    add_chiral_hs_for_drawing(&mut draw_mol);
    Ok(draw_mol)
}

impl DrawMol {
    fn from_molecule(
        mol: &Molecule,
        width: u32,
        height: u32,
        options: DrawOptions,
    ) -> Result<Self, SvgDrawError> {
        let mut draw_mol = prepare_working_mol_for_drawing(mol)?;
        if draw_mol.coords_2d().is_none() {
            draw_mol.compute_2d_coords().map_err(|err| {
                SvgDrawError::CoordinateGeneration(format!(
                    "compute2DCoords(canonOrient=true) failed: {err}"
                ))
            })?;
        }
        let coords = draw_mol
            .coords_2d()
            .expect("coords were generated immediately above");
        let at_cds = coords
            .iter()
            .map(|pt| DVec2::new(pt.x, -pt.y))
            .collect::<Vec<_>>();
        let implicit_hs = assign_valence(&draw_mol, ValenceModel::RdkitLike)
            .map(|assignment| assignment.implicit_hydrogens)
            .map_err(|err| {
                SvgDrawError::Unsupported(format!(
                    "DrawMol::extractAtomSymbols implicit-H valence assignment failed: {err}"
                ))
            })?;
        let mut out = Self {
            width: f64::from(width),
            height: f64::from(height),
            draw_width: f64::from(width),
            draw_height: f64::from(height),
            mol_height: f64::from(height),
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
            implicit_hs,
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
        };
        out.extract_atom_symbols(&draw_mol);
        out.extract_bonds(&draw_mol)?;
        out.smooth_bond_joins(&draw_mol);
        out.resolve_atom_symbol_clashes();
        out.extract_radicals(&draw_mol);
        out.calculate_scale();
        let drawn_font_size = out.font_size * out.scale;
        if !(6.0..=40.0).contains(&drawn_font_size) {
            out.font_size = if drawn_font_size > 40.0 {
                40.0 / out.scale
            } else {
                6.0 / out.scale
            };
            out.extract_atom_symbols(&draw_mol);
            out.extract_bonds(&draw_mol)?;
            out.smooth_bond_joins(&draw_mol);
            out.resolve_atom_symbol_clashes();
            out.extract_radicals(&draw_mol);
            out.find_extremes();
        }
        out.change_to_draw_coords();
        Ok(out)
    }

    fn extract_atom_symbols(&mut self, mol: &Molecule) {
        self.atom_labels.clear();
        self.atom_orients.clear();
        for atom in mol.atoms() {
            let orient = self.get_atom_orientation(mol, atom.index);
            self.atom_orients.push(orient);
            let symbol = self.get_atom_symbol(mol, atom, orient);
            if symbol.is_empty() {
                self.atom_labels.push(None);
            } else {
                self.atom_labels.push(Some(AtomLabel::new(
                    symbol,
                    atom.index,
                    atom.atomic_num,
                    orient,
                    self.at_cds[atom.index],
                    atom_colour(atom.atomic_num),
                    self.font_size,
                )));
            }
        }
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
        let wedge_bonds = pick_bonds_to_wedge_rdkit_subset(mol);
        let wedge_coords = self
            .at_cds
            .iter()
            .map(|pt| DVec2::new(pt.x, -pt.y))
            .collect::<Vec<_>>();

        for bond in mol.bonds() {
            match bond.order {
                BondOrder::Double | BondOrder::Aromatic => {
                    self.make_double_bond_lines(mol, bond, double_bond_offset);
                }
                BondOrder::Single => {
                    let wedge_dir = determine_bond_wedge_state_rdkit_subset(
                        bond,
                        &wedge_bonds,
                        &wedge_coords,
                        mol,
                    );
                    if matches!(
                        wedge_dir,
                        BondDirection::EndUpRight | BondDirection::EndDownRight
                    ) && wedge_bonds.contains_key(&bond.index)
                    {
                        self.make_wedged_bond(mol, bond, &wedge_bonds, wedge_dir)?;
                    } else {
                        let (begin, end) =
                            self.adjust_bond_ends_for_labels(bond.begin_atom, bond.end_atom);
                        self.new_bond_line(begin, end, bond);
                    }
                }
                BondOrder::Triple => {
                    let (begin, end) =
                        self.adjust_bond_ends_for_labels(bond.begin_atom, bond.end_atom);
                    self.new_bond_line(begin, end, bond);
                    self.make_triple_bond_lines(bond, double_bond_offset);
                }
                BondOrder::Quadruple => {
                    let (begin, end) =
                        self.adjust_bond_ends_for_labels(bond.begin_atom, bond.end_atom);
                    self.new_bond_line(begin, end, bond);
                }
                BondOrder::Null => {
                    self.make_bond_null_query_line(bond);
                }
                BondOrder::Dative => {
                    self.make_dative_bond(bond, double_bond_offset);
                }
            }
        }
        self.adjust_bonds_on_solid_wedge_ends(mol, &wedge_bonds, &wedge_coords);
        Ok(())
    }

    fn extract_radicals(&mut self, mol: &Molecule) {
        self.radicals.clear();
        for atom in mol.atoms() {
            if atom.num_radical_electrons == 0 {
                continue;
            }
            let (rect, orient) = self.calc_radical_rect(atom);
            self.radicals.push(DrawRadical {
                rect,
                orient,
                atom_idx: atom.index,
                count: atom.num_radical_electrons,
            });
        }
    }

    fn calc_radical_rect(&self, atom: &Atom) -> (StringRect, OrientType) {
        let spot_rad = self.radical_spot_radius_unscaled();
        let at_cds = self.at_cds[atom.index];
        let orient = self
            .atom_orients
            .get(atom.index)
            .copied()
            .unwrap_or(OrientType::C);
        let rad_size = (4.0 * f64::from(atom.num_radical_electrons) - 2.0) * spot_rad
            / self.font_scale_factor();
        let (x_min, x_max, y_min, y_max) =
            if let Some(label) = self.atom_labels.get(atom.index).and_then(Option::as_ref) {
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

    fn font_scale_factor(&self) -> f64 {
        self.font_size / 0.6
    }

    fn radical_spot_radius_unscaled(&self) -> f64 {
        0.2 * self.options.multiple_bond_offset * self.font_scale_factor()
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

    fn smooth_bond_joins(&mut self, mol: &Molecule) {
        let wedge_bonds = pick_bonds_to_wedge_rdkit_subset(mol);
        let wedge_coords = self
            .at_cds
            .iter()
            .map(|pt| DVec2::new(pt.x, -pt.y))
            .collect::<Vec<_>>();
        for atom in mol.atoms() {
            if self.atom_labels[atom.index].is_some() {
                continue;
            }
            let degree = atom_degree(mol, atom.index);
            let mut do_it = degree == 2;
            if !do_it && degree == 3 {
                for bond in mol.bonds() {
                    let nbr = if bond.begin_atom == atom.index {
                        Some(bond.end_atom)
                    } else if bond.end_atom == atom.index {
                        Some(bond.begin_atom)
                    } else {
                        None
                    };
                    if let Some(nbr) = nbr {
                        let wedge_dir = determine_bond_wedge_state_rdkit_subset(
                            bond,
                            &wedge_bonds,
                            &wedge_coords,
                            mol,
                        );
                        if (atom_degree(mol, nbr) == 1 && bond.order == BondOrder::Double)
                            || (wedge_bonds.contains_key(&bond.index)
                                && matches!(
                                    wedge_dir,
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
                let line1 = &self.bonds[line1_idx];
                let Some(p1) = line_endpoint_for_atom(line1, atom.index) else {
                    continue;
                };
                for j in 0..self.single_bond_lines.len() {
                    if i == j {
                        continue;
                    }
                    let line2_idx = self.single_bond_lines[j];
                    let line2 = &self.bonds[line2_idx];
                    let Some(p2) = line_endpoint_for_atom(line2, atom.index) else {
                        continue;
                    };
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
                if (label1.rects.len() > 1 || label2.rects.len() > 1)
                    && do_labels_clash(label1, label2)
                {
                    let mut idxs = [None, None];
                    if label1.rects.len() > 1 {
                        idxs[0] = Some(at_idx1);
                    }
                    if label2.rects.len() > 1 {
                        idxs[1] = Some(at_idx2);
                    }
                    if label1.rects.len() > label2.rects.len() {
                        idxs.swap(0, 1);
                    }
                    let first_ok = idxs[0]
                        .map(|idx| self.orient_atom_label(idx))
                        .unwrap_or(false);
                    if !first_ok && let Some(idx) = idxs[1] {
                        self.orient_atom_label(idx);
                    }
                }
            }
        }
    }

    fn orient_atom_label(&mut self, atom_idx: usize) -> bool {
        let Some(label) = self.atom_labels[atom_idx].as_ref() else {
            return false;
        };
        if label.rects.len() == 2 && label.rects[1].ch.is_ascii_lowercase() {
            return false;
        }
        let orig = label.orient;
        let pref = match orig {
            OrientType::S => 1,
            OrientType::W => 2,
            OrientType::E => 3,
            _ => 0,
        };
        const NEW_ORIENTS: [[OrientType; 3]; 4] = [
            [OrientType::S, OrientType::E, OrientType::W],
            [OrientType::N, OrientType::W, OrientType::E],
            [OrientType::E, OrientType::S, OrientType::N],
            [OrientType::W, OrientType::N, OrientType::S],
        ];
        let mut ok = false;
        for orient in NEW_ORIENTS[pref] {
            {
                let label = self.atom_labels[atom_idx]
                    .as_mut()
                    .expect("label existed above");
                label.orient = orient;
                label.recalculate_rects(self.font_size);
            }
            let label = self.atom_labels[atom_idx]
                .as_ref()
                .expect("label existed above")
                .clone();
            let clashes_labels = self.atom_labels.iter().enumerate().any(|(idx, other)| {
                idx != atom_idx
                    && other
                        .as_ref()
                        .is_some_and(|other| do_labels_clash(&label, other))
            });
            let clashes_bonds = self
                .bonds
                .iter()
                .any(|bond| label_rects_clash_with_line(&label, bond.begin, bond.end, 0.0))
                || self
                    .wedges
                    .iter()
                    .any(|wedge| label_rects_clash_with_wedge(&label, wedge, 0.0));
            if !clashes_labels && !clashes_bonds {
                ok = true;
                break;
            }
        }
        if !ok && let Some(label) = &mut self.atom_labels[atom_idx] {
            label.orient = orig;
            label.recalculate_rects(self.font_size);
        }
        ok
    }

    fn make_bond_null_query_line(&mut self, bond: &Bond) {
        // RDKit stores SMILES '~' bonds as query bonds with description
        // "BondNull"; MolDraw2D renders those using queryColour and
        // shortDashes without treating them as ordinary bond-order lines.
        self.push_bond_line(
            self.at_cds[bond.begin_atom],
            self.at_cds[bond.end_atom],
            self.options.query_colour,
            bond.begin_atom,
            bond.end_atom,
            bond.index,
            DashPattern::ShortDashes,
        );
    }

    fn make_dative_bond(&mut self, bond: &Bond, offset: f64) {
        let (end1, end2) = self.adjust_bond_ends_for_labels(bond.begin_atom, bond.end_atom);
        let mid = (end1 + end2) * 0.5;
        let col1 = atom_colour(self.atom_atomic_num(bond.begin_atom));
        let col2 = atom_colour(self.atom_atomic_num(bond.end_atom));
        self.push_bond_line(
            end1,
            mid,
            col1,
            bond.begin_atom,
            bond.end_atom,
            bond.index,
            DashPattern::NoDash,
        );
        let frac = 2.0 * offset / (end2 - end1).length();
        self.arrows.push(DrawArrow {
            begin: mid,
            end: end2,
            colour: col2,
            width: self.options.bond_line_width,
            frac,
            angle: std::f64::consts::PI / 12.0,
            atom1_idx: bond.begin_atom,
            atom2_idx: bond.end_atom,
            bond_idx: bond.index,
        });
        self.draw_items.push(DrawItem::Arrow(self.arrows.len() - 1));
    }

    fn make_wedged_bond(
        &mut self,
        mol: &Molecule,
        bond: &Bond,
        wedge_bonds: &HashMap<usize, usize>,
        wedge_dir: BondDirection,
    ) -> Result<(), SvgDrawError> {
        let Some(&from_atom) = wedge_bonds.get(&bond.index) else {
            return Ok(());
        };
        let to_atom = if bond.begin_atom == from_atom {
            bond.end_atom
        } else if bond.end_atom == from_atom {
            bond.begin_atom
        } else {
            return Err(SvgDrawError::Unsupported(format!(
                "DrawMol::makeWedgedBond wedge source atom {from_atom} is not on bond {}",
                bond.index
            )));
        };
        match wedge_dir {
            BondDirection::EndUpRight => self.make_solid_wedge(mol, from_atom, to_atom, bond),
            BondDirection::EndDownRight => self.make_dashed_wedge(mol, from_atom, to_atom, bond),
            BondDirection::None => {}
        }
        Ok(())
    }

    fn make_solid_wedge(&mut self, mol: &Molecule, from_atom: usize, to_atom: usize, bond: &Bond) {
        let old_mean = self.mean_bond_length;
        if self.atom_labels[from_atom].is_some() || self.atom_labels[to_atom].is_some() {
            self.mean_bond_length *= 2.0;
        }
        let (end1, end2) = self.adjust_bond_ends_for_labels(from_atom, to_atom);
        self.mean_bond_length = old_mean;

        let perp = calc_perpendicular(self.at_cds[from_atom], self.at_cds[to_atom]);
        let disp = perp * self.options.multiple_bond_offset * self.mean_bond_length / 2.0;
        let t1 = end2 + disp;
        let t2 = end2 - disp;
        let line_width = if self.options.bond_line_width < 1.0 {
            self.options.bond_line_width
        } else {
            self.options.bond_line_width / 2.0
        };
        let points = build_solid_wedge_points(
            vec![end1, t1, t2],
            atom_colour(self.atom_atomic_num(from_atom)),
            atom_colour(self.atom_atomic_num(to_atom)),
            false,
            self.find_other_bond_vecs(mol, to_atom, from_atom),
        );
        self.wedges.push(DrawWedge {
            points,
            col1: atom_colour(self.atom_atomic_num(from_atom)),
            col2: atom_colour(self.atom_atomic_num(to_atom)),
            width: line_width / 2.0,
            kind: WedgeKind::Solid,
            one_less_dash: false,
            atom1_idx: from_atom,
            atom2_idx: to_atom,
            bond_idx: bond.index,
        });
        self.draw_items.push(DrawItem::Wedge(self.wedges.len() - 1));
    }

    fn make_dashed_wedge(&mut self, mol: &Molecule, from_atom: usize, to_atom: usize, bond: &Bond) {
        let old_mean = self.mean_bond_length;
        if self.atom_labels[from_atom].is_some() || self.atom_labels[to_atom].is_some() {
            self.mean_bond_length *= 2.0;
        }
        let (end1, end2) = self.adjust_bond_ends_for_labels(from_atom, to_atom);
        self.mean_bond_length = old_mean;

        let perp = calc_perpendicular(self.at_cds[from_atom], self.at_cds[to_atom]);
        let disp = perp * self.options.multiple_bond_offset * self.mean_bond_length / 2.0;
        let t1 = end2 + disp;
        let t2 = end2 - disp;
        let line_width = if self.options.bond_line_width < 1.0 {
            self.options.bond_line_width
        } else {
            self.options.bond_line_width / 2.0
        };
        self.wedges.push(DrawWedge {
            points: vec![end1, t1, t2],
            col1: atom_colour(self.atom_atomic_num(from_atom)),
            col2: atom_colour(self.atom_atomic_num(to_atom)),
            width: line_width,
            kind: WedgeKind::Dashed,
            one_less_dash: atom_degree(mol, to_atom) > 1,
            atom1_idx: from_atom,
            atom2_idx: to_atom,
            bond_idx: bond.index,
        });
        self.draw_items.push(DrawItem::Wedge(self.wedges.len() - 1));
    }

    fn adjust_bonds_on_solid_wedge_ends(
        &mut self,
        mol: &Molecule,
        wedge_bonds: &HashMap<usize, usize>,
        wedge_coords: &[DVec2],
    ) {
        for bond in mol.bonds().iter().rev() {
            let wedge_dir =
                determine_bond_wedge_state_rdkit_subset(bond, wedge_bonds, wedge_coords, mol);
            if wedge_dir != BondDirection::EndUpRight || !wedge_bonds.contains_key(&bond.index) {
                continue;
            }
            let from_atom = wedge_bonds[&bond.index];
            let to_atom = if bond.begin_atom == from_atom {
                bond.end_atom
            } else {
                bond.begin_atom
            };
            let end_atom_degree = atom_degree(mol, to_atom);
            if end_atom_degree != 2 || self.atom_labels[to_atom].is_some() {
                continue;
            }
            let Some(third_atom) = self.other_neighbor(mol, to_atom, from_atom, 0) else {
                continue;
            };
            let Some(bond1) = bond_between_atoms(mol, to_atom, third_atom) else {
                continue;
            };
            if bond1.order == BondOrder::Triple {
                continue;
            }
            let b1 = direction_vector(self.at_cds[to_atom], self.at_cds[from_atom]);
            let b2 = direction_vector(self.at_cds[bond1.end_atom], self.at_cds[bond1.begin_atom]);
            if (1.0 - b1.dot(b2)).abs() < 0.001 {
                continue;
            }
            let Some(wedge_idx) = self.wedges.iter().position(|w| w.bond_idx == bond.index) else {
                continue;
            };
            let Some((bond_line_idx, _)) = self
                .bonds
                .iter()
                .enumerate()
                .filter(|(_, line)| line.bond_idx == bond1.index)
                .filter(|(_, line)| line.kind_is_simple())
                .min_by(|(_, a), (_, b)| {
                    let end_cds = self.at_cds[to_atom];
                    let da0 = (a.begin - end_cds).length_squared();
                    let da1 = (a.end - end_cds).length_squared();
                    let db0 = (b.begin - end_cds).length_squared();
                    let db1 = (b.end - end_cds).length_squared();
                    da0.min(da1).partial_cmp(&db0.min(db1)).unwrap()
                })
            else {
                continue;
            };
            let wedge = &self.wedges[wedge_idx];
            let p1 = if wedge.points.len() == 3 {
                1
            } else if wedge.points.len() == 9 {
                4
            } else {
                continue;
            };
            let p2 = if wedge.points.len() == 3 {
                2
            } else if wedge.points.len() == 9 {
                5
            } else {
                continue;
            };
            let p1_idx = if (self.at_cds[third_atom] - wedge.points[p1]).length_squared()
                < (self.at_cds[third_atom] - wedge.points[p2]).length_squared()
            {
                p2
            } else {
                p1
            };
            let p1_point = wedge.points[p1_idx];
            let bond_line = &mut self.bonds[bond_line_idx];
            if bond_line.atom1_idx == to_atom {
                bond_line.begin = p1_point;
            } else {
                bond_line.end = p1_point;
            }
        }
    }

    fn calc_mean_bond_length(&mut self, mol: &Molecule) {
        if mol.bonds().is_empty() {
            self.mean_bond_length = 0.0;
            return;
        }
        let total = mol
            .bonds()
            .iter()
            .map(|bond| (self.at_cds[bond.begin_atom] - self.at_cds[bond.end_atom]).length())
            .sum::<f64>();
        self.mean_bond_length = total / mol.bonds().len() as f64;
    }

    fn make_double_bond_lines(&mut self, mol: &Molecule, bond: &Bond, offset: f64) {
        let (end1, end2) = self.adjust_bond_ends_for_labels(bond.begin_atom, bond.end_atom);
        let sat1 = self.at_cds[bond.begin_atom];
        let sat2 = self.at_cds[bond.end_atom];
        self.at_cds[bond.begin_atom] = end1;
        self.at_cds[bond.end_atom] = end2;
        let (l1s, l1f, l2s, l2f) = self.calc_double_bond_lines(mol, bond, offset);
        self.at_cds[bond.begin_atom] = sat1;
        self.at_cds[bond.end_atom] = sat2;
        self.new_bond_line_from_points(l1s, l1f, bond);
        let col1 = atom_colour(self.atom_atomic_num(bond.begin_atom));
        let col2 = atom_colour(self.atom_atomic_num(bond.end_atom));
        let l1 = (l1s - l1f).length_squared();
        let l2 = (l2s - l2f).length_squared();
        if (atom_degree(mol, bond.begin_atom) == 1 || atom_degree(mol, bond.end_atom) == 1)
            && col1 != col2
            && (l1 - l2).abs() > 0.01
        {
            let mid_len = l1.sqrt() / 2.0;
            let not_mid = if atom_degree(mol, bond.begin_atom) == 1 {
                l2s + direction_vector(l2s, l2f) * mid_len
            } else {
                l2f + direction_vector(l2f, l2s) * mid_len
            };
            self.push_bond_line(
                l2s,
                not_mid,
                col1,
                bond.begin_atom,
                bond.end_atom,
                bond.index,
                DashPattern::NoDash,
            );
            self.push_bond_line(
                not_mid,
                l2f,
                col2,
                bond.begin_atom,
                bond.end_atom,
                bond.index,
                DashPattern::NoDash,
            );
        } else {
            self.new_bond_line_from_points(l2s, l2f, bond);
        }
    }

    fn make_triple_bond_lines(&mut self, bond: &Bond, offset: f64) {
        let (end1, end2) = self.adjust_bond_ends_for_labels(bond.begin_atom, bond.end_atom);
        let sat1 = self.at_cds[bond.begin_atom];
        let sat2 = self.at_cds[bond.end_atom];
        self.at_cds[bond.begin_atom] = end1;
        self.at_cds[bond.end_atom] = end2;
        let at1_cds = self.at_cds[bond.begin_atom];
        let at2_cds = self.at_cds[bond.end_atom];
        let perp = calc_perpendicular(at1_cds, at2_cds);
        self.new_bond_line_from_points(at1_cds + perp * offset, at2_cds + perp * offset, bond);
        self.new_bond_line_from_points(at1_cds - perp * offset, at2_cds - perp * offset, bond);
        self.at_cds[bond.begin_atom] = sat1;
        self.at_cds[bond.end_atom] = sat2;
    }

    fn calc_double_bond_lines(
        &mut self,
        mol: &Molecule,
        bond: &Bond,
        offset: f64,
    ) -> (DVec2, DVec2, DVec2, DVec2) {
        let at1_degree = atom_degree(mol, bond.begin_atom);
        let at2_degree = atom_degree(mol, bond.end_atom);
        if self.is_linear_atom(mol, bond.begin_atom)
            || self.is_linear_atom(mol, bond.end_atom)
            || (at1_degree == 1 && at2_degree == 1)
        {
            let at1_cds = self.at_cds[bond.begin_atom];
            let at2_cds = self.at_cds[bond.end_atom];
            let perp = calc_perpendicular(at1_cds, at2_cds) * offset * 0.5;
            (
                at1_cds + perp,
                at2_cds + perp,
                at1_cds - perp,
                at2_cds - perp,
            )
        } else if at1_degree == 1 || at2_degree == 1 {
            self.double_bond_terminal(mol, bond.begin_atom, bond.end_atom, offset)
        } else {
            let at1_cds = self.at_cds[bond.begin_atom];
            let at2_cds = self.at_cds[bond.end_atom];
            let mut l1s = at1_cds;
            let mut l1f = at2_cds;
            let (mut l2s, l2f) = if self.bond_ring_to_use(mol, bond).is_some() {
                self.bond_inside_ring(mol, bond, offset)
            } else if self.atom_labels[bond.begin_atom].is_some()
                && self.atom_labels[bond.end_atom].is_some()
            {
                let (dbl_l1s, dbl_l1f, dbl_l2s, dbl_l2f) =
                    self.double_bond_terminal(mol, bond.begin_atom, bond.end_atom, offset);
                l1s = dbl_l1s;
                l1f = dbl_l1f;
                (dbl_l2s, dbl_l2f)
            } else {
                self.bond_non_ring(mol, bond, offset)
            };
            let mut l1s_out = l1s;
            if !are_bonds_parallel(l1s, l1f, l2f, l2s, 1.0e-4) {
                std::mem::swap(&mut l1s_out, &mut l2s);
            }
            if bond.direction != BondDirection::None || bond.stereo == BondStereo::Any {
                std::mem::swap(&mut l1s_out, &mut l2s);
            }
            (l1s_out, l1f, l2s, l2f)
        }
    }

    fn bond_non_ring(&self, mol: &Molecule, bond: &Bond, offset: f64) -> (DVec2, DVec2) {
        let beg = bond.begin_atom;
        let end = bond.end_atom;
        let beg_trunc = self.atom_labels[beg].is_none();
        let end_trunc = self.atom_labels[end].is_none();
        let beg_degree = atom_degree(mol, beg);
        let end_degree = atom_degree(mol, end);
        if beg_degree == 2 && end_degree == 2 {
            let third = self
                .other_neighbor(mol, beg, end, 0)
                .expect("DrawMol::bondNonRing degree-2 branch requires begin-side neighbor");
            let fourth = self
                .other_neighbor(mol, end, beg, 0)
                .expect("DrawMol::bondNonRing degree-2 branch requires end-side neighbor");
            let l2s = self.double_bond_end(third, beg, end, offset, beg_trunc);
            let is_trans = are_bonds_trans(
                self.at_cds[third],
                self.at_cds[beg],
                self.at_cds[end],
                self.at_cds[fourth],
            );
            let l2f = if is_trans {
                let perp = calc_inner_perpendicular(
                    self.at_cds[end],
                    self.at_cds[beg],
                    self.at_cds[third],
                );
                self.at_cds[end] + perp * offset
            } else {
                self.double_bond_end(fourth, end, beg, offset, end_trunc)
            };
            return (l2s, l2f);
        }
        if beg_degree == 2 && end_degree > 2 {
            let third = self
                .other_neighbor(mol, beg, end, 0)
                .expect("DrawMol::bondNonRing degree-2/>2 branch requires begin-side neighbor");
            let mut fourth = self
                .other_neighbor(mol, end, beg, 0)
                .expect("DrawMol::bondNonRing degree-2/>2 branch requires end-side neighbor");
            let l2s = self.double_bond_end(third, beg, end, offset, beg_trunc);
            let is_trans = are_bonds_trans(
                self.at_cds[third],
                self.at_cds[beg],
                self.at_cds[end],
                self.at_cds[fourth],
            );
            if is_trans {
                fourth = self.non_colinear_neighbor(mol, end, beg);
            }
            let l2f = self.double_bond_end(fourth, end, beg, offset, end_trunc);
            return (l2s, l2f);
        }
        if beg_degree > 2 && end_degree == 2 {
            let mut third = self
                .other_neighbor(mol, beg, end, 0)
                .expect("DrawMol::bondNonRing >2/degree-2 branch requires begin-side neighbor");
            let fourth = self
                .other_neighbor(mol, end, beg, 0)
                .expect("DrawMol::bondNonRing >2/degree-2 branch requires end-side neighbor");
            let mut l2s = self.double_bond_end(third, beg, end, offset, beg_trunc);
            let is_trans = are_bonds_trans(
                self.at_cds[third],
                self.at_cds[beg],
                self.at_cds[end],
                self.at_cds[fourth],
            );
            if is_trans {
                third = self.non_colinear_neighbor(mol, beg, end);
                let trans_trunc = if bond.stereo == BondStereo::None {
                    beg_trunc
                } else {
                    end_trunc
                };
                l2s = self.double_bond_end(third, beg, end, offset, trans_trunc);
            }
            let l2f = self.double_bond_end(fourth, end, beg, offset, end_trunc);
            return (l2s, l2f);
        }
        if beg_degree > 2 && end_degree > 2 {
            let third = self
                .other_neighbor(mol, beg, end, 0)
                .expect("DrawMol::bondNonRing >2/>2 branch requires begin-side neighbor");
            let l2s = self.double_bond_end(third, beg, end, offset, beg_trunc);
            let mut fourth = self
                .other_neighbor(mol, end, beg, 0)
                .expect("DrawMol::bondNonRing >2/>2 branch requires end-side neighbor");
            let is_trans = are_bonds_trans(
                self.at_cds[third],
                self.at_cds[beg],
                self.at_cds[end],
                self.at_cds[fourth],
            );
            if is_trans {
                fourth = self.non_colinear_neighbor(mol, end, beg);
            }
            let l2f = self.double_bond_end(fourth, end, beg, offset, end_trunc);
            return (l2s, l2f);
        }
        panic!("DrawMol::bondNonRing reached unsupported degree branch");
    }

    fn non_colinear_neighbor(&self, mol: &Molecule, at1: usize, at2: usize) -> usize {
        let mut third = None;
        for i in 1..atom_degree(mol, at1) {
            third = self.other_neighbor(mol, at1, at2, i);
            if let Some(third_atom) = third
                && !are_bonds_parallel(
                    self.at_cds[at1],
                    self.at_cds[at2],
                    self.at_cds[at1],
                    self.at_cds[third_atom],
                    1.0e-4,
                )
            {
                return third_atom;
            }
        }
        third
            .or_else(|| self.other_neighbor(mol, at1, at2, 1))
            .expect("DrawMol::bondNonRing nonColinearNbor requires alternate neighbor")
    }

    fn bond_inside_ring(&self, mol: &Molecule, bond: &Bond, offset: f64) -> (DVec2, DVec2) {
        let Some(ring) = self.bond_ring_to_use(mol, bond) else {
            return self.bond_non_ring(mol, bond, offset);
        };
        let third = self
            .other_ring_atom(mol, bond.begin_atom, bond, &ring)
            .unwrap_or(bond.begin_atom);
        let fourth = self
            .other_ring_atom(mol, bond.end_atom, bond, &ring)
            .unwrap_or(bond.end_atom);
        let is_trans = are_bonds_trans(
            self.at_cds[third],
            self.at_cds[bond.begin_atom],
            self.at_cds[bond.end_atom],
            self.at_cds[fourth],
        );
        if is_trans {
            self.bond_non_ring(mol, bond, offset)
        } else {
            (
                self.double_bond_end(
                    third,
                    bond.begin_atom,
                    bond.end_atom,
                    offset,
                    self.atom_labels[bond.begin_atom].is_none(),
                ),
                self.double_bond_end(
                    fourth,
                    bond.end_atom,
                    bond.begin_atom,
                    offset,
                    self.atom_labels[bond.end_atom].is_none(),
                ),
            )
        }
    }

    fn bond_ring_to_use(&self, mol: &Molecule, bond: &Bond) -> Option<Vec<usize>> {
        let bond_rings = rdkit_bond_rings_preserve_order(mol);
        let mut rings = bond_rings
            .into_iter()
            .filter(|ring| ring.contains(&bond.index))
            .collect::<Vec<_>>();
        if rings.is_empty() {
            return None;
        }
        if rings.len() > 1 {
            for ring in &rings {
                if ring
                    .iter()
                    .all(|&bond_idx| mol.bonds()[bond_idx].is_aromatic == bond.is_aromatic)
                {
                    return Some(ring.clone());
                }
            }
        }
        Some(rings.remove(0))
    }

    fn other_ring_atom(
        &self,
        mol: &Molecule,
        bond_atom: usize,
        bond: &Bond,
        ring_bonds: &[usize],
    ) -> Option<usize> {
        for bond2 in mol.bonds() {
            if bond2.index == bond.index {
                continue;
            }
            if !ring_bonds.contains(&bond2.index) {
                continue;
            }
            if bond2.begin_atom == bond_atom {
                return Some(bond2.end_atom);
            }
            if bond2.end_atom == bond_atom {
                return Some(bond2.begin_atom);
            }
        }
        None
    }

    fn double_bond_terminal(
        &mut self,
        mol: &Molecule,
        mut at1: usize,
        mut at2: usize,
        mut offset: f64,
    ) -> (DVec2, DVec2, DVec2, DVec2) {
        let mut swapped = false;
        if atom_degree(mol, at1) > 1 && atom_degree(mol, at2) == 1 {
            std::mem::swap(&mut at1, &mut at2);
            swapped = true;
        }
        let at1_cds = self.at_cds[at1];
        let at2_cds = self.at_cds[at2];
        let (mut l1s, mut l1f, mut l2s, mut l2f) = if self.atom_labels[at2].is_some() {
            offset /= 2.0;
            let perp = calc_perpendicular(at1_cds, at2_cds) * offset;
            (
                at1_cds + perp,
                at2_cds + perp,
                at1_cds - perp,
                at2_cds - perp,
            )
        } else if atom_degree(mol, at2) > 2 {
            offset /= 2.0;
            let perp = calc_perpendicular(at1_cds, at2_cds) * offset;
            let l1s = at1_cds + perp;
            let mut l1f = at2_cds + perp;
            let l2s = at1_cds - perp;
            let mut l2f = at2_cds - perp;
            let bl = (l1s - l1f).length().max((l2s - l2f).length());
            l1f = l1s + direction_vector(l1s, l1f) * 2.0 * bl;
            l2f = l2s + direction_vector(l2s, l2f) * 2.0 * bl;
            for nbr in atom_neighbors(mol, at2) {
                let nbr_cds = self.at_cds[nbr];
                if let Some(ip) = line_intersection(l1s, l1f, at2_cds, nbr_cds) {
                    l1f = ip;
                }
                if let Some(ip) = line_intersection(l2s, l2f, at2_cds, nbr_cds) {
                    l2f = ip;
                }
            }
            (l1s, l1f, l2s, l2f)
        } else {
            let third_atom = self.other_neighbor(mol, at2, at1, 0).unwrap_or(at1);
            let perp = calc_inner_perpendicular(at1_cds, at2_cds, self.at_cds[third_atom]);
            let l1s = at1_cds;
            let l1f = at2_cds;
            let l2s = at1_cds + perp * offset;
            let mut l2f = self.double_bond_end(at1, at2, third_atom, offset, true);
            if direction_vector(l1s, l1f)
                .dot(direction_vector(l2s, l2f))
                .abs()
                < 0.9999
            {
                l2f = self.double_bond_end(at1, at2, third_atom, -offset, true);
            }
            if let Some(label) = &mut self.atom_labels[at1] {
                label.cds += perp * offset * 0.5;
            }
            (l1s, l1f, l2s, l2f)
        };
        if swapped {
            std::mem::swap(&mut l1s, &mut l1f);
            std::mem::swap(&mut l2s, &mut l2f);
        }
        (l1s, l1f, l2s, l2f)
    }

    fn double_bond_end(
        &self,
        at1: usize,
        at2: usize,
        at3: usize,
        offset: f64,
        trunc: bool,
    ) -> DVec2 {
        let v21 = direction_vector(self.at_cds[at2], self.at_cds[at1]);
        let v23 = direction_vector(self.at_cds[at2], self.at_cds[at3]);
        let mut v23_perp = DVec2::new(-v23.y, v23.x).normalize_or_zero();
        let mut bis = v21 + v23;
        if bis.length_squared() < 1.0e-6 {
            return self.at_cds[at2] - v23_perp * offset;
        }
        bis = bis.normalize();
        if v23_perp.dot(bis) < 0.0 {
            v23_perp *= -1.0;
        }
        if trunc
            && let Some(ip) = line_intersection(
                self.at_cds[at2],
                self.at_cds[at2] + bis,
                self.at_cds[at2] + v23_perp * offset,
                self.at_cds[at3] + v23_perp * offset,
            )
        {
            return ip;
        }
        self.at_cds[at2] + v23_perp * offset
    }

    fn is_linear_atom(&self, mol: &Molecule, atom_idx: usize) -> bool {
        if atom_degree(mol, atom_idx) != 2 {
            return false;
        }
        let mut nbrs = Vec::new();
        for bond in mol.bonds() {
            if bond.begin_atom == atom_idx {
                nbrs.push((bond.end_atom, bond.order));
            } else if bond.end_atom == atom_idx {
                nbrs.push((bond.begin_atom, bond.order));
            }
        }
        if nbrs.len() != 2 || nbrs[0].1 != nbrs[1].1 {
            return false;
        }
        direction_vector(self.at_cds[atom_idx], self.at_cds[nbrs[0].0]).dot(direction_vector(
            self.at_cds[atom_idx],
            self.at_cds[nbrs[1].0],
        )) < -0.95
    }

    fn other_neighbor(
        &self,
        mol: &Molecule,
        first_atom: usize,
        second_atom: usize,
        nbor_num: usize,
    ) -> Option<usize> {
        let mut count = 0;
        for bond in mol.bonds() {
            let nbr = if bond.begin_atom == first_atom {
                Some(bond.end_atom)
            } else if bond.end_atom == first_atom {
                Some(bond.begin_atom)
            } else {
                None
            };
            if let Some(nbr) = nbr
                && nbr != second_atom
            {
                if count == nbor_num {
                    return Some(nbr);
                }
                count += 1;
            }
        }
        None
    }

    fn find_other_bond_vecs(&self, mol: &Molecule, atom: usize, other_atom: usize) -> Vec<DVec2> {
        if atom_degree(mol, atom) == 1 || self.atom_labels[atom].is_some() {
            return Vec::new();
        }
        let mut other_bond_vecs = Vec::new();
        for i in 1..atom_degree(mol, atom) {
            let Some(third_atom) = self.other_neighbor(mol, atom, other_atom, i - 1) else {
                continue;
            };
            let Some(bond) = bond_between_atoms(mol, atom, third_atom) else {
                continue;
            };
            if bond.order == BondOrder::Triple {
                continue;
            }
            if bond.order == BondOrder::Double
                && atom_degree(mol, atom) > 2
                && atom_degree(mol, third_atom) == 1
            {
                continue;
            }
            other_bond_vecs.push(direction_vector(self.at_cds[atom], self.at_cds[third_atom]));
        }
        other_bond_vecs
    }

    fn new_bond_line(&mut self, begin: DVec2, end: DVec2, bond: &Bond) {
        self.new_bond_line_from_points(begin, end, bond);
    }

    fn new_bond_line_from_points(&mut self, begin: DVec2, end: DVec2, bond: &Bond) {
        let col1 = atom_colour(self.atom_atomic_num(bond.begin_atom));
        let col2 = atom_colour(self.atom_atomic_num(bond.end_atom));
        if col1 == col2 {
            self.push_bond_line(
                begin,
                end,
                col1,
                bond.begin_atom,
                bond.end_atom,
                bond.index,
                DashPattern::NoDash,
            );
        } else {
            let mid = (begin + end) / 2.0;
            self.push_bond_line(
                begin,
                mid,
                col1,
                bond.begin_atom,
                bond.end_atom,
                bond.index,
                DashPattern::NoDash,
            );
            self.push_bond_line(
                mid,
                end,
                col2,
                bond.begin_atom,
                bond.end_atom,
                bond.index,
                DashPattern::NoDash,
            );
        }
    }

    fn push_bond_line(
        &mut self,
        begin: DVec2,
        end: DVec2,
        colour: DrawColour,
        atom1_idx: usize,
        atom2_idx: usize,
        bond_idx: usize,
        dash_pattern: DashPattern,
    ) {
        self.bonds.push(DrawLine {
            begin,
            end,
            colour,
            width: self.options.bond_line_width,
            scale_width: self.options.scale_bond_width,
            dash_pattern,
            atom1_idx,
            atom2_idx,
            bond_idx,
        });
        self.draw_items.push(DrawItem::Line(self.bonds.len() - 1));
        if dash_pattern == DashPattern::NoDash {
            self.single_bond_lines.push(self.bonds.len() - 1);
        }
    }

    fn atom_atomic_num(&self, atom_idx: usize) -> u8 {
        self.atom_labels
            .get(atom_idx)
            .and_then(Option::as_ref)
            .map_or(6, |label| label.atomic_num)
    }

    fn adjust_bond_ends_for_labels(&self, beg_at_idx: usize, end_at_idx: usize) -> (DVec2, DVec2) {
        let padding = 0.033 * self.mean_bond_length;
        let mut beg_cds = self.at_cds[beg_at_idx];
        let mut end_cds = self.at_cds[end_at_idx];
        if let Some(label) = &self.atom_labels[beg_at_idx] {
            adjust_bond_end_for_string(end_cds, padding, &label.rects, label.cds, &mut beg_cds);
        }
        if let Some(label) = &self.atom_labels[end_at_idx] {
            adjust_bond_end_for_string(beg_cds, padding, &label.rects, label.cds, &mut end_cds);
        }
        (beg_cds, end_cds)
    }

    fn calculate_scale(&mut self) {
        self.find_extremes();
        self.draw_width = self.width * (1.0 - 2.0 * self.margin_padding);
        self.draw_height = self.height * (1.0 - 2.0 * self.margin_padding);
        self.mol_height = self.draw_height;
        let mut new_scale = 1.0;
        if self.x_range > 1e-4 || self.y_range > 1e-4 {
            new_scale = (self.draw_width / self.x_range).min(self.mol_height / self.y_range);
        }
        self.scale *= new_scale / self.scale;
    }

    fn find_extremes(&mut self) {
        self.x_min = f64::MAX / 2.0;
        self.x_max = f64::MIN / 2.0;
        self.y_min = f64::MAX / 2.0;
        self.y_max = f64::MIN / 2.0;
        for bond in &self.bonds {
            self.x_min = self.x_min.min(bond.begin.x).min(bond.end.x);
            self.x_max = self.x_max.max(bond.begin.x).max(bond.end.x);
            self.y_min = self.y_min.min(bond.begin.y).min(bond.end.y);
            self.y_max = self.y_max.max(bond.begin.y).max(bond.end.y);
        }
        for wedge in &self.wedges {
            for point in &wedge.points {
                self.x_min = self.x_min.min(point.x);
                self.x_max = self.x_max.max(point.x);
                self.y_min = self.y_min.min(point.y);
                self.y_max = self.y_max.max(point.y);
            }
        }
        for arrow in &self.arrows {
            for point in [arrow.begin, arrow.end] {
                self.x_min = self.x_min.min(point.x);
                self.x_max = self.x_max.max(point.x);
                self.y_min = self.y_min.min(point.y);
                self.y_max = self.y_max.max(point.y);
            }
        }
        for polyline in &self.join_paths {
            for point in &polyline.points {
                self.x_min = self.x_min.min(point.x);
                self.x_max = self.x_max.max(point.x);
                self.y_min = self.y_min.min(point.y);
                self.y_max = self.y_max.max(point.y);
            }
        }
        for polyline in &self.post_shapes {
            for point in &polyline.points {
                self.x_min = self.x_min.min(point.x);
                self.x_max = self.x_max.max(point.x);
                self.y_min = self.y_min.min(point.y);
                self.y_max = self.y_max.max(point.y);
            }
        }
        for label in self.atom_labels.iter().flatten() {
            label.find_extremes(
                &mut self.x_min,
                &mut self.x_max,
                &mut self.y_min,
                &mut self.y_max,
            );
        }
        for radical in &self.radicals {
            find_rect_extremes(
                &radical.rect,
                &mut self.x_min,
                &mut self.x_max,
                &mut self.y_min,
                &mut self.y_max,
            );
        }
        self.x_range = self.x_max - self.x_min;
        self.y_range = self.y_max - self.y_min;
        if self.x_range < 1e-4 {
            self.x_range = 2.0;
            self.x_min -= 1.0;
            self.x_max += 1.0;
        }
        if self.y_range < 1e-4 {
            self.y_range = 2.0;
            self.y_min -= 1.0;
            self.y_max += 1.0;
        }
    }

    fn change_to_draw_coords(&mut self) {
        let trans = DVec2::new(-self.x_min, -self.y_min);
        let scale = DVec2::splat(self.scale);
        let scaled_ranges = DVec2::new(self.scale * self.x_range, self.scale * self.y_range);
        let to_centre = DVec2::new(
            (self.draw_width - scaled_ranges.x) / 2.0 + self.width * self.margin_padding,
            (self.mol_height - scaled_ranges.y) / 2.0 + self.height * self.margin_padding,
        );
        for bond in &mut self.bonds {
            bond.begin = transform_point(bond.begin, trans, scale, to_centre);
            bond.end = transform_point(bond.end, trans, scale, to_centre);
        }
        for wedge in &mut self.wedges {
            for point in &mut wedge.points {
                *point = transform_point(*point, trans, scale, to_centre);
            }
        }
        for arrow in &mut self.arrows {
            arrow.begin = transform_point(arrow.begin, trans, scale, to_centre);
            arrow.end = transform_point(arrow.end, trans, scale, to_centre);
        }
        for polyline in &mut self.join_paths {
            for point in &mut polyline.points {
                *point = transform_point(*point, trans, scale, to_centre);
            }
        }
        for label in self.atom_labels.iter_mut().flatten() {
            label.cds = transform_point(label.cds, trans, scale, to_centre);
            label.recalculate_rects(self.font_size * self.scale);
        }
        for radical in &mut self.radicals {
            radical.rect.trans = transform_point(radical.rect.trans, trans, scale, to_centre);
            radical.rect.width *= self.scale;
            radical.rect.height *= self.scale;
        }
        self.extract_close_contacts(trans, scale, to_centre);
    }

    fn extract_close_contacts(&mut self, trans: DVec2, scale: DVec2, to_centre: DVec2) {
        self.post_shapes.clear();
        if self.options.flag_close_contacts_dist < 0 {
            return;
        }
        let tol = f64::from(
            self.options.flag_close_contacts_dist * self.options.flag_close_contacts_dist,
        );
        let mut flagged = vec![false; self.at_cds.len()];
        let offset = DVec2::splat(0.1 * self.scale);
        for i in 0..self.at_cds.len() {
            if flagged[i] {
                continue;
            }
            let ci = transform_point(self.at_cds[i], trans, scale, to_centre);
            for j in (i + 1)..self.at_cds.len() {
                if flagged[j] {
                    continue;
                }
                let cj = transform_point(self.at_cds[j], trans, scale, to_centre);
                if (cj - ci).length_squared() <= tol {
                    flagged[i] = true;
                    flagged[j] = true;
                    break;
                }
            }
            if flagged[i] {
                let p1 = ci - offset;
                let p2 = ci + offset;
                self.post_shapes.push(DrawPolyline {
                    points: vec![p1, DVec2::new(p1.x, p2.y), p2, DVec2::new(p2.x, p1.y), p1],
                    colour: DrawColour {
                        r: 1.0,
                        g: 0.0,
                        b: 0.0,
                        a: 1.0,
                    },
                    width: self.options.bond_line_width,
                    scale_width: false,
                    atom1_idx: Some(i),
                    atom2_idx: None,
                    bond_idx: None,
                });
            }
        }
    }

    fn get_atom_orientation(&self, mol: &Molecule, atom_idx: usize) -> OrientType {
        const VERT_SLOPE: f64 = 2.7474774194546216;
        let at1_cds = self.at_cds[atom_idx];
        let mut nbr_sum = DVec2::ZERO;
        for bond in mol.bonds() {
            let other = if bond.begin_atom == atom_idx {
                Some(bond.end_atom)
            } else if bond.end_atom == atom_idx {
                Some(bond.begin_atom)
            } else {
                None
            };
            if let Some(other) = other {
                nbr_sum += self.at_cds[other] - at1_cds;
            }
        }
        let degree = atom_degree(mol, atom_idx);
        if degree == 0 {
            let atomic_num = mol.atoms()[atom_idx].atomic_num;
            return if matches!(atomic_num, 8 | 9 | 16 | 17 | 34 | 35 | 52 | 53 | 84 | 85) {
                OrientType::W
            } else {
                OrientType::E
            };
        }
        let mut islope = 1000.0;
        if nbr_sum.x.abs() > 1e-4 {
            islope = nbr_sum.y / nbr_sum.x;
        }
        let mut orient = if islope.abs() <= VERT_SLOPE {
            if nbr_sum.x > 0.0 {
                OrientType::W
            } else {
                OrientType::E
            }
        } else if nbr_sum.y > 0.0 {
            OrientType::S
        } else {
            OrientType::N
        };
        if matches!(orient, OrientType::N | OrientType::S) && degree == 1 {
            orient = if islope.abs() > VERT_SLOPE {
                OrientType::E
            } else if nbr_sum.x > 0.0 {
                OrientType::W
            } else {
                OrientType::E
            };
        }
        orient
    }

    fn get_atom_symbol(&self, mol: &Molecule, atom: &Atom, orient: OrientType) -> String {
        let symbol = crate::periodic_table::element_symbol(atom.atomic_num).unwrap_or("?");
        let isotope = atom
            .isotope
            .map(|isotope| format!("<sup>{isotope}</sup>"))
            .unwrap_or_default();
        let map_num = atom
            .atom_map_num
            .map(|map_num| format!(":{map_num}"))
            .unwrap_or_default();
        let num_h = if atom.atomic_num == 6 && atom_degree(mol, atom.index) > 0 {
            0
        } else {
            atom.explicit_hydrogens as usize
                + self
                    .implicit_hs
                    .get(atom.index)
                    .copied()
                    .unwrap_or_default() as usize
        };
        let h = match num_h {
            0 => String::new(),
            1 => "H".to_string(),
            n => format!("H<sub>{n}</sub>"),
        };
        let charge = if atom.formal_charge == 0 {
            String::new()
        } else {
            let magnitude = atom.formal_charge.unsigned_abs();
            let sign = if atom.formal_charge > 0 { "+" } else { "-" };
            if magnitude > 1 {
                format!("<sup>{magnitude}{sign}</sup>")
            } else {
                format!("<sup>{sign}</sup>")
            }
        };
        let _ = orient;
        if isotope.is_empty() && h.is_empty() && charge.is_empty() && map_num.is_empty() {
            if self.is_linear_atom(mol, atom.index)
                || atom.atomic_num != 6
                || atom_degree(mol, atom.index) == 0
            {
                symbol.to_string()
            } else {
                String::new()
            }
        } else {
            format!("{isotope}{symbol}{charge}{h}{map_num}")
        }
    }
}

impl AtomLabel {
    fn new(
        symbol: String,
        atom_idx: usize,
        atomic_num: u8,
        orient: OrientType,
        cds: DVec2,
        colour: DrawColour,
        font_size: f64,
    ) -> Self {
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
        for rect in &self.rects {
            let mut shifted = rect.clone();
            shifted.trans += self.cds;
            let (tl, tr, br, bl) = shifted.calc_corners(0.0);
            *xmin = xmin.min(tr.x).min(bl.x);
            *ymin = ymin.min(tr.y).min(bl.y);
            *xmax = xmax.max(tr.x).max(bl.x);
            *ymax = ymax.max(tr.y).max(bl.y);
            let _ = (tl, br);
        }
    }
}

impl StringRect {
    fn new(
        ch: char,
        draw_mode: TextDrawType,
        offset: DVec2,
        g_centre: DVec2,
        width: f64,
        height: f64,
    ) -> Self {
        Self {
            ch,
            draw_mode,
            trans: DVec2::ZERO,
            offset,
            g_centre,
            y_shift: 0.0,
            width,
            height,
            rect_corr: 0.0,
        }
    }

    fn calc_centre(&self) -> DVec2 {
        let mut c = DVec2::new(
            (self.trans.x + self.g_centre.x) - self.offset.x,
            (self.trans.y + self.g_centre.y) - self.offset.y,
        );
        c.y -= self.y_shift;
        c
    }

    fn calc_corners(&self, padding: f64) -> (DVec2, DVec2, DVec2, DVec2) {
        let wb2 = padding + self.width / 2.0;
        let hb2 = padding + self.height / 2.0;
        let c = self.calc_centre();
        (
            DVec2::new(c.x - wb2, c.y - hb2),
            DVec2::new(c.x + wb2, c.y - hb2),
            DVec2::new(c.x + wb2, c.y + hb2),
            DVec2::new(c.x - wb2, c.y + hb2),
        )
    }
}

fn init_drawing(out: &mut String, width: u32, height: u32) {
    out.push_str("<?xml version='1.0' encoding='iso-8859-1'?>\n");
    out.push_str("<svg version='1.1' baseProfile='full'\n");
    out.push_str("              xmlns='http://www.w3.org/2000/svg'\n");
    out.push_str("                      xmlns:cosmolkit='https://www.cosmol.org'\n");
    out.push_str("                      xmlns:xlink='http://www.w3.org/1999/xlink'\n");
    out.push_str("                  xml:space='preserve'\n");
    out.push_str(&format!(
        "width='{width}px' height='{height}px' viewBox='0 0 {width} {height}'>\n"
    ));
    out.push_str("<!-- END OF HEADER -->\n");
}

fn clear_drawing(out: &mut String, width: u32, height: u32, colour: DrawColour) {
    out.push_str("<rect");
    out.push_str(&format!(
        " style='opacity:1.0;fill:{};stroke:none'",
        draw_colour_to_svg(colour)
    ));
    out.push_str(&format!(
        " width='{}' height='{}'",
        format_double(f64::from(width)),
        format_double(f64::from(height))
    ));
    out.push_str(" x='0.0' y='0.0'");
    out.push_str("> </rect>\n");
}

fn draw_line(out: &mut String, line: &DrawLine, scale: f64) {
    let width = if line.scale_width {
        (line.width * scale * 0.02).max(0.0)
    } else {
        line.width
    };
    out.push_str("<path ");
    out.push_str(&format!(
        "class='bond-{} atom-{} atom-{}' ",
        line.bond_idx, line.atom1_idx, line.atom2_idx
    ));
    out.push_str(&format!(
        "d='M {},{} L {},{}' ",
        format_double(line.begin.x),
        format_double(line.begin.y),
        format_double(line.end.x),
        format_double(line.end.y)
    ));
    out.push_str(&format!(
        "style='fill:none;fill-rule:evenodd;stroke:{};stroke-width:{}px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1",
        draw_colour_to_svg(line.colour),
        format_double(width)
    ));
    match line.dash_pattern {
        DashPattern::NoDash => {}
        DashPattern::ShortDashes => out.push_str(";stroke-dasharray:2,2"),
    }
    out.push('\'');
    out.push_str(" />\n");
}

fn draw_wedge(out: &mut String, wedge: &DrawWedge) {
    if wedge.points.len() < 3 {
        return;
    }
    if wedge.kind == WedgeKind::Dashed {
        draw_dashed_wedge(out, wedge);
        return;
    }
    if wedge.points.len() == 3 || wedge.points.len() == 9 {
        draw_solid_wedge_polygon(out, wedge, &wedge.points[0..3], wedge.col1);
    }
    if wedge.points.len() == 6 {
        draw_solid_wedge_polygon(
            out,
            wedge,
            &[
                wedge.points[0],
                wedge.points[1],
                wedge.points[2],
                wedge.points[5],
            ],
            wedge.col1,
        );
    } else if wedge.points.len() == 9 {
        draw_solid_wedge_polygon(
            out,
            wedge,
            &[
                wedge.points[4],
                wedge.points[5],
                wedge.points[6],
                wedge.points[7],
            ],
            wedge.col2,
        );
    }
}

fn draw_solid_wedge_polygon(
    out: &mut String,
    wedge: &DrawWedge,
    points: &[DVec2],
    colour: DrawColour,
) {
    if points.is_empty() {
        return;
    }
    out.push_str("<path ");
    out.push_str(&format!(
        "class='bond-{} atom-{} atom-{}' ",
        wedge.bond_idx, wedge.atom1_idx, wedge.atom2_idx
    ));
    out.push_str("d='M");
    out.push_str(&format!(
        " {},{}",
        format_double(points[0].x),
        format_double(points[0].y)
    ));
    for point in &points[1..] {
        out.push_str(&format!(
            " L {},{}",
            format_double(point.x),
            format_double(point.y)
        ));
    }
    out.push_str(" Z' ");
    out.push_str(&format!(
        "style='fill:{};fill-rule:evenodd;fill-opacity:{};stroke:{};stroke-width:{}px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:{};'",
        draw_colour_to_svg(colour),
        format_alpha(colour.a),
        draw_colour_to_svg(colour),
        format_double(wedge.width),
        format_alpha(colour.a)
    ));
    out.push_str(" />\n");
}

fn draw_dashed_wedge(out: &mut String, wedge: &DrawWedge) {
    let at1 = wedge.points[0];
    let end1 = wedge.points[1];
    let end2 = wedge.points[2];
    let midend = (end1 + end2) * 0.5;
    let mut e1 = direction_vector(at1, end1);
    let mut e2 = direction_vector(at1, end2);
    let mut dash_sep = 2.5 + wedge.width;
    let central_len = (at1 - midend).length();
    let mut n_dashes = (central_len / dash_sep).round() as usize;
    let num_dashes_needed = if wedge.one_less_dash { 4 } else { 3 };
    if n_dashes < num_dashes_needed {
        n_dashes = num_dashes_needed;
    }
    if n_dashes == 0 {
        draw_wedge_dash_line(out, wedge, end1, end2, wedge.col1);
        return;
    }
    dash_sep = central_len / n_dashes as f64;
    if wedge.one_less_dash {
        let end_len_by_2 = (end1 - end2).length() / 2.0;
        let central_line = direction_vector(at1, midend);
        let central_perp = DVec2::new(-central_line.y, central_line.x);
        let new_end1 = at1 + central_line * (central_len - dash_sep) + central_perp * end_len_by_2;
        let new_end2 = at1 + central_line * (central_len - dash_sep) - central_perp * end_len_by_2;
        e1 = direction_vector(at1, new_end1);
        e2 = direction_vector(at1, new_end2);
    }
    dash_sep *= (end1 - at1).length() / central_len;
    let extra = if wedge.one_less_dash { 0 } else { 1 };
    for i in 1..(n_dashes + extra) {
        let p1 = at1 + e1 * i as f64 * dash_sep;
        let p2 = at1 + e2 * i as f64 * dash_sep;
        let colour = if i > n_dashes / 2 {
            wedge.col2
        } else {
            wedge.col1
        };
        draw_wedge_dash_line(out, wedge, p1, p2, colour);
    }
}

fn draw_wedge_dash_line(
    out: &mut String,
    wedge: &DrawWedge,
    begin: DVec2,
    end: DVec2,
    colour: DrawColour,
) {
    out.push_str("<path ");
    out.push_str(&format!(
        "class='bond-{} atom-{} atom-{}' ",
        wedge.bond_idx, wedge.atom1_idx, wedge.atom2_idx
    ));
    out.push_str(&format!(
        "d='M {},{} L {},{}' ",
        format_double(begin.x),
        format_double(begin.y),
        format_double(end.x),
        format_double(end.y)
    ));
    out.push_str(&format!(
        "style='fill:none;fill-rule:evenodd;stroke:{};stroke-width:{}px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:{}'",
        draw_colour_to_svg(colour),
        format_double(wedge.width),
        format_alpha(colour.a)
    ));
    out.push_str(" />\n");
}

fn draw_arrow(out: &mut String, arrow: &DrawArrow) {
    let (line_end, tip, p1, p2) = arrow_points(arrow);
    let line = DrawLine {
        begin: arrow.begin,
        end: line_end,
        colour: arrow.colour,
        width: arrow.width,
        scale_width: false,
        dash_pattern: DashPattern::NoDash,
        atom1_idx: arrow.atom1_idx,
        atom2_idx: arrow.atom2_idx,
        bond_idx: arrow.bond_idx,
    };
    draw_line(out, &line, 1.0);
    let wedge = DrawWedge {
        points: vec![p1, tip, p2],
        col1: arrow.colour,
        col2: arrow.colour,
        width: arrow.width,
        kind: WedgeKind::Solid,
        one_less_dash: false,
        atom1_idx: arrow.atom1_idx,
        atom2_idx: arrow.atom2_idx,
        bond_idx: arrow.bond_idx,
    };
    draw_wedge(out, &wedge);
}

fn arrow_points(arrow: &DrawArrow) -> (DVec2, DVec2, DVec2, DVec2) {
    let mut arrow_end = arrow.end;
    let adjuster = 0.5 * arrow.width / arrow.angle.sin();
    let len = (arrow.begin - arrow.end).length();
    if len > 1.0e-6 {
        let adj_len = len - adjuster;
        arrow_end = arrow.begin + (arrow.end - arrow.begin) * adj_len / len;
    }
    let delta = arrow.begin - arrow_end;
    let cos_angle = arrow.angle.cos();
    let sin_angle = arrow.angle.sin();
    let frac = arrow.frac / cos_angle;
    let p1 = DVec2::new(
        arrow_end.x + frac * (delta.x * cos_angle + delta.y * sin_angle),
        arrow_end.y + frac * (delta.y * cos_angle - delta.x * sin_angle),
    );
    let p2 = DVec2::new(
        arrow_end.x + frac * (delta.x * cos_angle - delta.y * sin_angle),
        arrow_end.y + frac * (delta.y * cos_angle + delta.x * sin_angle),
    );
    (arrow_end, arrow_end, p1, p2)
}

fn draw_polyline(out: &mut String, polyline: &DrawPolyline, scale: f64) {
    let width = if polyline.scale_width {
        (polyline.width * scale * 0.02).max(0.0)
    } else {
        polyline.width
    };
    let mut points = polyline.points.iter();
    let Some(first) = points.next() else {
        return;
    };
    out.push_str("<path ");
    if polyline.atom1_idx.is_some() || polyline.atom2_idx.is_some() || polyline.bond_idx.is_some() {
        out.push_str("class='");
        let mut need_space = false;
        if let Some(bond_idx) = polyline.bond_idx {
            out.push_str(&format!("bond-{bond_idx}"));
            need_space = true;
        }
        if let Some(atom1_idx) = polyline.atom1_idx {
            if need_space {
                out.push(' ');
            }
            out.push_str(&format!("atom-{atom1_idx}"));
            need_space = true;
        }
        if let Some(atom2_idx) = polyline.atom2_idx
            && Some(atom2_idx) != polyline.atom1_idx
        {
            if need_space {
                out.push(' ');
            }
            out.push_str(&format!("atom-{atom2_idx}"));
        }
        out.push_str("' ");
    }
    out.push_str("d='M ");
    out.push_str(&format!(
        "{},{}",
        format_double(first.x),
        format_double(first.y)
    ));
    for point in points {
        out.push_str(&format!(
            " L {},{}",
            format_double(point.x),
            format_double(point.y)
        ));
    }
    out.push_str(&format!(
        "' style='fill:none;stroke:{};stroke-width:{}px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:{};' />\n",
        draw_colour_to_svg(polyline.colour),
        format_double(width),
        format_alpha(polyline.colour.a)
    ));
}

fn draw_atom_label(out: &mut String, label: &AtomLabel, base_font_size: f64) {
    for rect in &label.rects {
        let draw_cds = DVec2::new(
            label.cds.x + rect.trans.x - rect.offset.x,
            label.cds.y - rect.trans.y + rect.offset.y - rect.rect_corr - rect.y_shift,
        );
        let font_size =
            (base_font_size * select_scale_factor(rect.ch, rect.draw_mode) + 1e-9) as u32;
        out.push_str("<text");
        out.push_str(&format!(" x='{}'", format_double(draw_cds.x)));
        out.push_str(&format!(" y='{}'", format_double(draw_cds.y)));
        out.push_str(&format!(" class='atom-{}'", label.atom_idx));
        out.push_str(&format!(
            " style='font-size:{}px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:\"{}\",sans-serif;text-anchor:start;fill:{}'",
            font_size,
            EMBEDDED_DRAW_FONT_FAMILY,
            draw_colour_to_svg(label.colour)
        ));
        out.push_str(" >");
        out.push_str(&xml_escape(&rect.ch.to_string()));
        out.push_str("</text>\n");
    }
}

fn draw_radicals(out: &mut String, radicals: &[DrawRadical], spot_rad: f64) {
    for radical in radicals {
        let rect = &radical.rect;
        let dir = !matches!(
            radical.orient,
            OrientType::N | OrientType::S | OrientType::C
        );
        draw_radical_spots(
            out,
            radical.atom_idx,
            rect.trans,
            radical.count,
            if dir { rect.height } else { rect.width },
            dir,
            spot_rad,
        );
    }
}

fn draw_radical_spots(
    out: &mut String,
    atom_idx: usize,
    cds: DVec2,
    num_spots: u8,
    width: f64,
    vertical: bool,
    spot_rad: f64,
) {
    let mut ncds = cds;
    match num_spots {
        3 => {
            if vertical {
                ncds.y = cds.y - 0.6 * width + spot_rad;
            } else {
                ncds.x = cds.x - 0.6 * width + spot_rad;
            }
            draw_filled_circle_path(out, atom_idx, ncds, spot_rad);
            if vertical {
                ncds.y = cds.y + 0.6 * width - spot_rad;
            } else {
                ncds.x = cds.x + 0.6 * width - spot_rad;
            }
            draw_filled_circle_path(out, atom_idx, ncds, spot_rad);
            draw_filled_circle_path(out, atom_idx, cds, spot_rad);
        }
        1 => draw_filled_circle_path(out, atom_idx, cds, spot_rad),
        4 => {
            if vertical {
                ncds.y = cds.y + 6.0 * spot_rad;
            } else {
                ncds.x = cds.x + 6.0 * spot_rad;
            }
            draw_filled_circle_path(out, atom_idx, ncds, spot_rad);
            if vertical {
                ncds.y = cds.y - 6.0 * spot_rad;
            } else {
                ncds.x = cds.x - 6.0 * spot_rad;
            }
            draw_filled_circle_path(out, atom_idx, ncds, spot_rad);
            if vertical {
                ncds.y = cds.y + 2.0 * spot_rad;
            } else {
                ncds.x = cds.x + 2.0 * spot_rad;
            }
            draw_filled_circle_path(out, atom_idx, ncds, spot_rad);
            if vertical {
                ncds.y = cds.y - 2.0 * spot_rad;
            } else {
                ncds.x = cds.x - 2.0 * spot_rad;
            }
            draw_filled_circle_path(out, atom_idx, ncds, spot_rad);
        }
        2 => {
            if vertical {
                ncds.y = cds.y + 2.0 * spot_rad;
            } else {
                ncds.x = cds.x + 2.0 * spot_rad;
            }
            draw_filled_circle_path(out, atom_idx, ncds, spot_rad);
            if vertical {
                ncds.y = cds.y - 2.0 * spot_rad;
            } else {
                ncds.x = cds.x - 2.0 * spot_rad;
            }
            draw_filled_circle_path(out, atom_idx, ncds, spot_rad);
        }
        _ => {}
    }
}

fn draw_filled_circle_path(out: &mut String, atom_idx: usize, centre: DVec2, radius: f64) {
    let num_steps = 1 + (360.0_f64 / 5.0) as usize;
    let angle_incr = 360.0_f64 / num_steps as f64 * std::f64::consts::PI / 180.0;
    out.push_str("<path ");
    out.push_str(&format!("class='atom-{atom_idx}' "));
    out.push_str("d='M");
    for i in 0..=num_steps {
        let angle = i as f64 * angle_incr;
        let point = DVec2::new(
            centre.x + radius * angle.cos(),
            centre.y + radius * angle.sin(),
        );
        out.push_str(&format!(
            " {},{}",
            format_double(point.x),
            format_double(point.y)
        ));
        if i < num_steps {
            out.push_str(" L");
        }
    }
    out.push_str(&format!(
        " L {},{} Z' ",
        format_double(centre.x),
        format_double(centre.y)
    ));
    out.push_str(
        "style='fill:#000000;fill-rule:evenodd;fill-opacity:1;stroke:#000000;stroke-width:0.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' />\n",
    );
}

fn transform_point(point: DVec2, trans: DVec2, scale: DVec2, to_centre: DVec2) -> DVec2 {
    let mut out = point;
    out.x += trans.x;
    out.y += trans.y;
    out.x *= scale.x;
    out.y *= scale.y;
    out.x += to_centre.x;
    out.y += to_centre.y;
    out
}

fn calc_perpendicular(begin: DVec2, end: DVec2) -> DVec2 {
    let dir = begin - end;
    let len = dir.length();
    if len == 0.0 {
        return DVec2::ZERO;
    }
    DVec2::new(-dir.y / len, dir.x / len)
}

fn calc_inner_perpendicular(cds1: DVec2, cds2: DVec2, cds3: DVec2) -> DVec2 {
    let mut perp = calc_perpendicular(cds1, cds2);
    let v1 = cds1 - cds2;
    let v2 = cds2 - cds3;
    let obv = v1 - v2;
    if obv.dot(perp) < 0.0 {
        perp *= -1.0;
    }
    perp
}

fn direction_vector(from: DVec2, to: DVec2) -> DVec2 {
    (to - from).normalize_or_zero()
}

fn build_solid_wedge_points(
    mut points: Vec<DVec2>,
    col1: DrawColour,
    col2: DrawColour,
    split_bonds: bool,
    mut other_bond_vecs: Vec<DVec2>,
) -> Vec<DVec2> {
    if other_bond_vecs.len() > 2 {
        trim_other_bond_vecs(&mut other_bond_vecs);
    }
    if other_bond_vecs.len() == 2 {
        order_other_bond_vecs(&points, &mut other_bond_vecs);
    }
    if col1 != col2 || split_bonds {
        build_two_colour_wedge_triangles(&points, &other_bond_vecs)
    } else {
        build_single_colour_wedge_triangles(&mut points, &other_bond_vecs);
        points
    }
}

fn build_single_colour_wedge_triangles(points: &mut Vec<DVec2>, other_bond_vecs: &[DVec2]) {
    let point = points[0];
    let end1 = points[1];
    let end2 = points[2];
    let mid_end = (end1 + end2) / 2.0;
    let mut adjend1 = end1;
    let mut adjend2 = end2;
    points.clear();
    if other_bond_vecs.is_empty() {
        points.push(point);
        points.push(adjend1);
        points.push(adjend2);
    } else if other_bond_vecs.len() == 1 {
        let side1 = (end1 - point) * 2.0;
        if let Some(ip) = line_intersection_rdkit(
            point,
            point + side1,
            mid_end - other_bond_vecs[0],
            mid_end + other_bond_vecs[0],
        ) {
            adjend1 = ip;
        }
        let side2 = (end2 - point) * 2.0;
        if let Some(ip) = line_intersection_rdkit(
            point,
            point + side2,
            mid_end - other_bond_vecs[0],
            mid_end + other_bond_vecs[0],
        ) {
            adjend2 = ip;
        }
        points.push(point);
        points.push(adjend1);
        points.push(adjend2);
    } else if other_bond_vecs.len() == 2 {
        let side1 = (end1 - point) * 2.0;
        if let Some(ip) = line_intersection_rdkit(
            point,
            point + side1,
            mid_end - other_bond_vecs[0],
            mid_end + other_bond_vecs[0],
        ) {
            adjend1 = ip;
        }
        points.push(point);
        points.push(adjend1);
        points.push(mid_end);
        let side2 = (end2 - point) * 2.0;
        if let Some(ip) = line_intersection_rdkit(
            point,
            point + side2,
            mid_end - other_bond_vecs[1],
            mid_end + other_bond_vecs[1],
        ) {
            adjend2 = ip;
        }
        points.push(point);
        points.push(mid_end);
        points.push(adjend2);
    }
}

fn build_two_colour_wedge_triangles(points: &[DVec2], other_bond_vecs: &[DVec2]) -> Vec<DVec2> {
    let point = points[0];
    let end1 = points[1];
    let end2 = points[2];
    let mid_end = (end1 + end2) / 2.0;
    let mut adjend1 = end1;
    let mut adjend2 = end2;
    let e1 = end1 - point;
    let e2 = end2 - point;
    let mid1 = point + e1 * 0.5;
    let mid2 = point + e2 * 0.5;
    let mut out = vec![point, mid1, mid2];
    if other_bond_vecs.is_empty() {
        out.extend([mid1, adjend2, adjend1, mid1, mid2, adjend2]);
    } else if other_bond_vecs.len() == 1 {
        let side1 = (end1 - point) * 2.0;
        if let Some(ip) = line_intersection_rdkit(
            point,
            point + side1,
            mid_end - other_bond_vecs[0],
            mid_end + other_bond_vecs[0],
        ) {
            adjend1 = ip;
        }
        let side2 = (end2 - point) * 2.0;
        if let Some(ip) = line_intersection_rdkit(
            point,
            point + side2,
            mid_end - other_bond_vecs[0],
            mid_end + other_bond_vecs[0],
        ) {
            adjend2 = ip;
        }
        out.extend([mid1, adjend2, adjend1, mid1, mid2, adjend2]);
    } else if other_bond_vecs.len() == 2 {
        let side1 = (end1 - point) * 2.0;
        if let Some(ip) = line_intersection_rdkit(
            point,
            point + side1,
            mid_end - other_bond_vecs[0],
            mid_end + other_bond_vecs[0],
        ) {
            adjend1 = ip;
        }
        let side2 = (end2 - point) * 2.0;
        if let Some(ip) = line_intersection_rdkit(
            point,
            point + side2,
            mid_end - other_bond_vecs[1],
            mid_end + other_bond_vecs[1],
        ) {
            adjend2 = ip;
        }
        out.extend([
            mid1, adjend1, mid_end, mid_end, mid2, mid1, mid_end, adjend2, mid2,
        ]);
    }
    out
}

fn trim_other_bond_vecs(other_bond_vecs: &mut Vec<DVec2>) {
    if other_bond_vecs.len() < 3 {
        return;
    }
    let mut first_vec = 0;
    let mut second_vec = 1;
    let mut largest_ang = -361.0;
    for i in 0..other_bond_vecs.len() - 1 {
        for j in i + 1..other_bond_vecs.len() {
            let ang = angle_to(other_bond_vecs[i], other_bond_vecs[j]);
            if ang > largest_ang {
                first_vec = i;
                second_vec = j;
                largest_ang = ang;
            }
        }
    }
    *other_bond_vecs = vec![other_bond_vecs[first_vec], other_bond_vecs[second_vec]];
}

fn order_other_bond_vecs(points: &[DVec2], other_bond_vecs: &mut [DVec2]) {
    if other_bond_vecs.len() < 2 {
        return;
    }
    let mid = (points[1] + points[2]) / 2.0;
    let midp1 = direction_vector(mid, points[1]);
    let dot1 = midp1.dot(other_bond_vecs[0]);
    let dot2 = midp1.dot(other_bond_vecs[1]);
    if dot1 < dot2 {
        other_bond_vecs.swap(0, 1);
    }
}

fn angle_to(v1: DVec2, v2: DVec2) -> f64 {
    let t1 = v1.normalize_or_zero();
    let t2 = v2.normalize_or_zero();
    t1.dot(t2).clamp(-1.0, 1.0).acos()
}

fn atom_degree(mol: &Molecule, atom_idx: usize) -> usize {
    mol.bonds()
        .iter()
        .filter(|bond| bond.begin_atom == atom_idx || bond.end_atom == atom_idx)
        .count()
}

fn add_chiral_hs_for_drawing(mol: &mut Molecule) {
    let atom_rings = rdkit_atom_rings_for_chiral_hs(mol);
    let Ok(valence) = assign_valence(mol, ValenceModel::RdkitLike) else {
        return;
    };
    let mut chiral_atoms = Vec::new();
    for atom in mol.atoms() {
        let in_ring_count = atom_rings
            .iter()
            .filter(|ring| ring.contains(&atom.index))
            .count();
        let h_count = atom.explicit_hydrogens + valence.implicit_hydrogens[atom.index];
        if in_ring_count > 1
            && h_count > 0
            && matches!(
                atom.chiral_tag,
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
            )
        {
            chiral_atoms.push(atom.index);
        }
    }
    if chiral_atoms.is_empty() {
        return;
    }
    for atom in mol.atoms_mut() {
        atom.rdkit_cip_rank = None;
    }
    for atom_idx in chiral_atoms {
        if mol.atoms()[atom_idx].explicit_hydrogens > 0 {
            mol.atoms_mut()[atom_idx].explicit_hydrogens -= 1;
        }
        let h_idx = mol.add_atom(Atom {
            index: 0,
            atomic_num: 1,
            is_aromatic: false,
            formal_charge: 0,
            explicit_hydrogens: 0,
            no_implicit: false,
            num_radical_electrons: 0,
            chiral_tag: ChiralTag::Unspecified,
            isotope: None,
            atom_map_num: None,
            props: Default::default(),
            rdkit_cip_rank: None,
        });
        mol.add_bond(Bond {
            index: 0,
            begin_atom: atom_idx,
            end_atom: h_idx,
            order: BondOrder::Single,
            is_aromatic: false,
            direction: BondDirection::None,
            stereo: BondStereo::None,
            stereo_atoms: Vec::new(),
        });
    }
    mol.rebuild_adjacency();
}

fn atom_neighbors(mol: &Molecule, atom_idx: usize) -> Vec<usize> {
    mol.bonds()
        .iter()
        .filter_map(|bond| {
            if bond.begin_atom == atom_idx {
                Some(bond.end_atom)
            } else if bond.end_atom == atom_idx {
                Some(bond.begin_atom)
            } else {
                None
            }
        })
        .collect()
}

fn bond_between_atoms(mol: &Molecule, atom1: usize, atom2: usize) -> Option<&Bond> {
    mol.bonds().iter().find(|bond| {
        (bond.begin_atom == atom1 && bond.end_atom == atom2)
            || (bond.begin_atom == atom2 && bond.end_atom == atom1)
    })
}

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
    if point_idx == 0 { line.begin } else { line.end }
}

fn are_bonds_trans(at1: DVec2, at2: DVec2, at3: DVec2, at4: DVec2) -> bool {
    let v21 = at1 - at2;
    let v34 = at4 - at3;
    v21.dot(v34) < 0.0
}

fn are_bonds_parallel(at1: DVec2, at2: DVec2, at3: DVec2, at4: DVec2, tol: f64) -> bool {
    let v21 = direction_vector(at1, at2);
    let v34 = direction_vector(at4, at3);
    (1.0 - v21.dot(v34).abs()).abs() < tol
}

fn do_labels_clash(label1: &AtomLabel, label2: &AtomLabel) -> bool {
    for rect1 in &label1.rects {
        let mut shifted = rect1.clone();
        shifted.trans += label1.cds;
        if label_rects_intersect(&label2.rects, label2.cds, &shifted, 0.0) {
            return true;
        }
    }
    false
}

fn label_rects_intersect(
    rects: &[StringRect],
    cds: DVec2,
    other: &StringRect,
    padding: f64,
) -> bool {
    for rect in rects {
        let mut shifted = rect.clone();
        shifted.trans += cds;
        if rects_intersect(&shifted, other, padding) {
            return true;
        }
    }
    false
}

fn rects_intersect(a: &StringRect, b: &StringRect, padding: f64) -> bool {
    let (a_tl, a_tr, a_br, a_bl) = a.calc_corners(padding);
    let (b_tl, b_tr, b_br, b_bl) = b.calc_corners(padding);
    let a_min_x = a_tl.x.min(a_tr.x).min(a_br.x).min(a_bl.x);
    let a_max_x = a_tl.x.max(a_tr.x).max(a_br.x).max(a_bl.x);
    let a_min_y = a_tl.y.min(a_tr.y).min(a_br.y).min(a_bl.y);
    let a_max_y = a_tl.y.max(a_tr.y).max(a_br.y).max(a_bl.y);
    let b_min_x = b_tl.x.min(b_tr.x).min(b_br.x).min(b_bl.x);
    let b_max_x = b_tl.x.max(b_tr.x).max(b_br.x).max(b_bl.x);
    let b_min_y = b_tl.y.min(b_tr.y).min(b_br.y).min(b_bl.y);
    let b_max_y = b_tl.y.max(b_tr.y).max(b_br.y).max(b_bl.y);
    a_min_x <= b_max_x && a_max_x >= b_min_x && a_min_y <= b_max_y && a_max_y >= b_min_y
}

fn find_rect_extremes(
    rect: &StringRect,
    xmin: &mut f64,
    xmax: &mut f64,
    ymin: &mut f64,
    ymax: &mut f64,
) {
    *xmax = xmax.max(rect.trans.x + rect.width / 2.0);
    *xmin = xmin.min(rect.trans.x - rect.width / 2.0);
    *ymax = ymax.max(rect.trans.y + rect.height / 2.0);
    *ymin = ymin.min(rect.trans.y - rect.height / 2.0);
}

fn radical_rect_at(trans: DVec2, width: f64, height: f64) -> StringRect {
    StringRect {
        ch: ' ',
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

#[allow(clippy::too_many_arguments)]
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
    match orient {
        OrientType::N | OrientType::C => radical_rect_at(
            DVec2::new(at_cds.x, y_max + 0.5 * spot_rad * 3.0),
            rad_size,
            spot_rad * 3.0,
        ),
        OrientType::S => radical_rect_at(
            DVec2::new(at_cds.x, y_min - 0.5 * spot_rad * 3.0),
            rad_size,
            spot_rad * 3.0,
        ),
        OrientType::E => radical_rect_at(
            DVec2::new(x_max + 3.0 * spot_rad, at_cds.y),
            spot_rad * 1.5,
            rad_size,
        ),
        OrientType::W => radical_rect_at(
            DVec2::new(x_min - 3.0 * spot_rad, at_cds.y),
            spot_rad * 1.5,
            rad_size,
        ),
    }
}

fn rect_clashes_with_line(rect: &StringRect, begin: DVec2, end: DVec2, padding: f64) -> bool {
    let (tl, tr, br, bl) = rect.calc_corners(padding);
    line_intersection(begin, end, tl, tr).is_some()
        || line_intersection(begin, end, tr, br).is_some()
        || line_intersection(begin, end, br, bl).is_some()
        || line_intersection(begin, end, bl, tl).is_some()
}

fn label_rects_clash_with_line(label: &AtomLabel, begin: DVec2, end: DVec2, padding: f64) -> bool {
    for rect in &label.rects {
        let mut shifted = rect.clone();
        shifted.trans += label.cds;
        if rect_clashes_with_line(&shifted, begin, end, padding) {
            return true;
        }
    }
    false
}

fn label_rects_clash_with_wedge(label: &AtomLabel, wedge: &DrawWedge, padding: f64) -> bool {
    let padding = padding * wedge.width;
    for rect in &label.rects {
        let mut shifted = rect.clone();
        shifted.trans += label.cds;
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

fn rect_clashes_with_triangle(
    rect: &StringRect,
    pt1: DVec2,
    pt2: DVec2,
    pt3: DVec2,
    padding: f64,
) -> bool {
    if rect_contains_point(rect, pt1, padding)
        || rect_contains_point(rect, pt2, padding)
        || rect_contains_point(rect, pt3, padding)
    {
        return true;
    }
    let (tl, tr, br, bl) = rect.calc_corners(padding);
    if point_in_triangle(tl, pt1, pt2, pt3)
        || point_in_triangle(tr, pt1, pt2, pt3)
        || point_in_triangle(br, pt1, pt2, pt3)
        || point_in_triangle(bl, pt1, pt2, pt3)
    {
        return true;
    }
    line_intersection_rdkit(tl, tr, pt1, pt2).is_some()
        || line_intersection_rdkit(tl, tr, pt2, pt3).is_some()
        || line_intersection_rdkit(tl, tr, pt3, pt1).is_some()
        || line_intersection_rdkit(tr, br, pt1, pt2).is_some()
        || line_intersection_rdkit(tr, br, pt2, pt3).is_some()
        || line_intersection_rdkit(tr, br, pt3, pt1).is_some()
        || line_intersection_rdkit(br, bl, pt1, pt2).is_some()
        || line_intersection_rdkit(br, bl, pt2, pt3).is_some()
        || line_intersection_rdkit(br, bl, pt3, pt1).is_some()
        || line_intersection_rdkit(bl, tl, pt1, pt2).is_some()
        || line_intersection_rdkit(bl, tl, pt2, pt3).is_some()
        || line_intersection_rdkit(bl, tl, pt3, pt1).is_some()
}

fn rect_contains_point(rect: &StringRect, pt: DVec2, padding: f64) -> bool {
    let (mut tl, mut tr, mut br, mut bl) = rect.calc_corners(padding);
    if tl.y < bl.y {
        std::mem::swap(&mut tl, &mut bl);
        std::mem::swap(&mut tr, &mut br);
    }
    pt.x >= tl.x && pt.x <= br.x && pt.y >= br.y && pt.y <= tl.y
}

fn point_in_triangle(pt: DVec2, t1: DVec2, t2: DVec2, t3: DVec2) -> bool {
    let d = (t2.y - t3.y) * (t1.x - t3.x) + (t3.x - t2.x) * (t1.y - t3.y);
    let a = ((t2.y - t3.y) * (pt.x - t3.x) + (t3.x - t2.x) * (pt.y - t3.y)) / d;
    let b = ((t3.y - t1.y) * (pt.x - t3.x) + (t1.x - t3.x) * (pt.y - t3.y)) / d;
    let c = 1.0 - a - b;
    (0.0..=1.0).contains(&a) && (0.0..=1.0).contains(&b) && (0.0..=1.0).contains(&c)
}

fn get_string_rects(text: &str, orient: OrientType, font_size: f64) -> Vec<StringRect> {
    let text_bits = atom_label_to_pieces(text, orient);
    if orient == OrientType::W {
        let mut label = String::new();
        for bit in text_bits.iter().rev() {
            label.push_str(bit);
        }
        let mut rects = get_string_rects_unsplit(&label, font_size);
        align_string(TextAlignType::End, &mut rects);
        return rects;
    }
    if orient == OrientType::E {
        let mut label = String::new();
        for bit in &text_bits {
            label.push_str(bit);
        }
        let mut rects = get_string_rects_unsplit(&label, font_size);
        align_string(TextAlignType::Start, &mut rects);
        return rects;
    }
    let mut rects = Vec::new();
    let mut running_y = 0.0;
    let text_align = if orient == OrientType::C {
        TextAlignType::Middle
    } else {
        TextAlignType::Middle
    };
    for bit in &text_bits {
        let mut bit_rects = get_string_rects_unsplit(bit, font_size);
        align_string(text_align, &mut bit_rects);
        let max_height = bit_rects
            .iter()
            .map(|rect| rect.height)
            .fold(f64::MIN, f64::max);
        for rect in &mut bit_rects {
            rect.y_shift = running_y;
        }
        rects.extend(bit_rects);
        if orient == OrientType::N {
            running_y -= 1.1 * max_height;
        } else if orient == OrientType::S {
            running_y += 1.1 * max_height;
        }
    }
    rects
}

fn get_string_rects_unsplit(text: &str, act_font_size: f64) -> Vec<StringRect> {
    let (draw_chars, draw_modes) = parse_draw_chars(text);
    let max_width = draw_chars
        .iter()
        .map(|&c| char_width(c))
        .fold(0.0_f64, f64::max);
    let mut running_x = 0.0;
    let mut rects = Vec::new();
    for (&ch, &draw_mode) in draw_chars.iter().zip(draw_modes.iter()) {
        let mut char_width = 0.6 * act_font_size * char_width(ch) / max_width;
        let mut char_height = if ch == '+' {
            0.6 * act_font_size
        } else if ch == '-' {
            0.4 * act_font_size
        } else {
            0.8 * act_font_size
        };
        let cscale = select_scale_factor(ch, draw_mode);
        char_height *= cscale;
        char_width *= cscale;
        let mut offset = DVec2::new(char_width / 2.0, char_height / 2.0);
        if ch == '+' || ch == '-' {
            offset.y /= 2.0;
        }
        let g_centre = DVec2::new(char_width / 2.0, char_height / 2.0);
        let mut rect = StringRect::new(ch, draw_mode, offset, g_centre, char_width, char_height);
        rect.trans.x += running_x;
        if draw_mode != TextDrawType::Normal {
            running_x += char_width * 1.05;
        } else {
            running_x += char_width * 1.15;
        }
        rects.push(rect);
    }
    for rect in &mut rects {
        rect.g_centre.y = act_font_size - rect.g_centre.y;
        rect.offset.y = act_font_size / 2.0;
    }
    adjust_string_rects_for_super_sub_script(&mut rects);
    rects
}

fn align_string(align: TextAlignType, rects: &mut [StringRect]) {
    if rects.is_empty() {
        return;
    }
    // DrawTextSVG derives from DrawTextNotFT; its alignString() semantics
    // differ from the FreeType/base DrawText path and are required for SVG
    // coordinate parity.
    let num_norm = rects
        .iter()
        .filter(|rect| rect.draw_mode == TextDrawType::Normal)
        .count();
    let align = if align == TextAlignType::Middle && num_norm == 1 {
        TextAlignType::Start
    } else {
        align
    };
    let (align_trans, align_offset) = match align {
        TextAlignType::Start | TextAlignType::End => {
            let mut align_char = 0;
            for (idx, rect) in rects.iter().enumerate() {
                if rect.draw_mode == TextDrawType::Normal {
                    align_char = idx;
                    if align == TextAlignType::Start {
                        break;
                    }
                }
            }
            (rects[align_char].trans, rects[align_char].offset)
        }
        TextAlignType::Middle => {
            let mut x_min = f64::MAX;
            let mut x_max = f64::MIN;
            let mut align_offset = DVec2::ZERO;
            let mut num_norm = 0usize;
            for rect in rects.iter() {
                if rect.draw_mode == TextDrawType::Normal {
                    let (_tl, tr, _br, bl) = rect.calc_corners(0.0);
                    x_min = x_min.min(bl.x).min(tr.x);
                    x_max = x_max.max(bl.x).max(tr.x);
                    align_offset += rect.offset;
                    num_norm += 1;
                }
            }
            if num_norm > 0 {
                align_offset /= num_norm as f64;
            }
            (DVec2::new((x_max - x_min) / 2.0, 0.0), align_offset)
        }
    };
    for rect in rects {
        rect.trans -= align_trans;
        rect.offset = align_offset;
    }
}

fn parse_draw_chars(text: &str) -> (Vec<char>, Vec<TextDrawType>) {
    let chars = text.chars().collect::<Vec<_>>();
    let mut draw_chars = Vec::new();
    let mut draw_modes = Vec::new();
    let mut draw_mode = TextDrawType::Normal;
    let mut i = 0;
    while i < chars.len() {
        let rest = chars[i..].iter().collect::<String>();
        if rest.starts_with("<sub>") {
            draw_mode = TextDrawType::Subscript;
            i += 5;
            continue;
        }
        if rest.starts_with("<sup>") {
            draw_mode = TextDrawType::Superscript;
            i += 5;
            continue;
        }
        if rest.starts_with("</sub>") {
            draw_mode = TextDrawType::Normal;
            i += 6;
            continue;
        }
        if rest.starts_with("</sup>") {
            draw_mode = TextDrawType::Normal;
            i += 6;
            continue;
        }
        draw_chars.push(chars[i]);
        draw_modes.push(draw_mode);
        i += 1;
    }
    (draw_chars, draw_modes)
}

fn select_scale_factor(ch: char, draw_type: TextDrawType) -> f64 {
    match draw_type {
        TextDrawType::Normal => 1.0,
        TextDrawType::Subscript => 0.66,
        TextDrawType::Superscript => {
            if ch == '-' || ch == '+' {
                0.66
            } else {
                0.66
            }
        }
    }
}

fn adjust_string_rects_for_super_sub_script(rects: &mut [StringRect]) {
    let mut last_char: Option<usize> = None;
    for i in 0..rects.len() {
        match rects[i].draw_mode {
            TextDrawType::Superscript => {
                if last_char.is_none() {
                    last_char =
                        (i + 1..rects.len()).find(|&j| rects[j].draw_mode == TextDrawType::Normal);
                }
                if let Some(last) = last_char {
                    rects[i].rect_corr = rects[last].height;
                    rects[i].trans.y -= rects[i].rect_corr / 2.0;
                }
            }
            TextDrawType::Subscript => {
                if let Some(last) = last_char {
                    rects[i].rect_corr = -rects[last].height;
                    rects[i].trans.y -= rects[i].rect_corr / 2.0;
                }
            }
            TextDrawType::Normal => last_char = Some(i),
        }
    }

    for i in 1..rects.len() {
        let prev = rects[i - 1].draw_mode;
        let curr = rects[i].draw_mode;
        if (curr == TextDrawType::Subscript && prev == TextDrawType::Superscript)
            || (prev == TextDrawType::Subscript && curr == TextDrawType::Superscript)
        {
            let move_back = rects[i].trans.x - rects[i - 1].trans.x;
            for rect in &mut rects[i..] {
                rect.trans.x -= move_back;
            }
        }
    }
}

fn atom_label_to_pieces(label: &str, orient: OrientType) -> Vec<String> {
    if label.is_empty() {
        return Vec::new();
    }
    let chars = label.chars().collect::<Vec<_>>();
    let mut label_pieces = Vec::new();
    let mut next_piece = String::new();
    let mut i = 0;
    while i < chars.len() {
        let rest = chars[i..].iter().collect::<String>();
        if (rest.starts_with("<s") || chars[i] == ':' || chars[i].is_ascii_uppercase())
            && !next_piece.is_empty()
        {
            label_pieces.push(std::mem::take(&mut next_piece));
        }
        next_piece.push(chars[i]);
        i += 1;
    }
    if !next_piece.is_empty() {
        label_pieces.push(next_piece);
    }
    if label_pieces.len() < 2 {
        return label_pieces;
    }
    if orient == OrientType::E || orient == OrientType::S {
        for i in 0..label_pieces.len() {
            if label_pieces[i] == "<sup>+</sup>" || label_pieces[i] == "<sup>-</sup>" {
                let piece = std::mem::take(&mut label_pieces[i]);
                label_pieces.push(piece);
                break;
            }
        }
    }

    let mut final_pieces = Vec::new();
    let mut curr_piece = String::new();
    let mut had_symbol = false;
    for piece in label_pieces {
        if piece.is_empty() {
            continue;
        }
        if !piece.as_bytes()[0].is_ascii_uppercase() {
            curr_piece.push_str(&piece);
        } else if had_symbol {
            final_pieces.push(curr_piece);
            curr_piece = piece;
            had_symbol = true;
        } else {
            curr_piece.push_str(&piece);
            had_symbol = true;
        }
    }
    if !curr_piece.is_empty() {
        final_pieces.push(curr_piece);
    }
    final_pieces
}

fn adjust_bond_end_for_string(
    end2: DVec2,
    padding: f64,
    rects: &[StringRect],
    label_pos: DVec2,
    move_end: &mut DVec2,
) {
    for rect in rects {
        let mut shifted = rect.clone();
        shifted.trans += label_pos;
        let (tl, tr, br, bl) = shifted.calc_corners(padding);
        if let Some(ip) = line_intersection(*move_end, end2, tl, tr) {
            *move_end = ip;
        }
        if let Some(ip) = line_intersection(*move_end, end2, tr, br) {
            *move_end = ip;
        }
        if let Some(ip) = line_intersection(*move_end, end2, br, bl) {
            *move_end = ip;
        }
        if let Some(ip) = line_intersection(*move_end, end2, bl, tl) {
            *move_end = ip;
        }
    }
}

fn line_intersection(a1: DVec2, a2: DVec2, b1: DVec2, b2: DVec2) -> Option<DVec2> {
    let r = a2 - a1;
    let s = b2 - b1;
    let denom = cross(r, s);
    if denom.abs() < 1e-12 {
        return None;
    }
    let u = cross(b1 - a1, r) / denom;
    let t = cross(b1 - a1, s) / denom;
    if (-1e-12..=1.0 + 1e-12).contains(&t) && (-1e-12..=1.0 + 1e-12).contains(&u) {
        Some(a1 + r * t)
    } else {
        None
    }
}

fn line_intersection_rdkit(l1s: DVec2, l1f: DVec2, l2s: DVec2, l2f: DVec2) -> Option<DVec2> {
    let s1_x = l1f.x - l1s.x;
    let s1_y = l1f.y - l1s.y;
    let s2_x = l2f.x - l2s.x;
    let s2_y = l2f.y - l2s.y;
    let d = -s2_x * s1_y + s1_x * s2_y;
    if d == 0.0 {
        return None;
    }
    let s = (-s1_y * (l1s.x - l2s.x) + s1_x * (l1s.y - l2s.y)) / d;
    let t = (s2_x * (l1s.y - l2s.y) - s2_y * (l1s.x - l2s.x)) / d;
    if (0.0..=1.0).contains(&s) && (0.0..=1.0).contains(&t) {
        Some(DVec2::new(l1s.x + t * s1_x, l1s.y + t * s1_y))
    } else {
        None
    }
}

fn cross(a: DVec2, b: DVec2) -> f64 {
    a.x * b.y - a.y * b.x
}

fn atom_colour(atomic_num: u8) -> DrawColour {
    match atomic_num {
        0 => DrawColour {
            r: 0.1,
            g: 0.1,
            b: 0.1,
            a: 1.0,
        },
        7 => DrawColour {
            r: 0.0,
            g: 0.0,
            b: 1.0,
            a: 1.0,
        },
        8 => DrawColour {
            r: 1.0,
            g: 0.0,
            b: 0.0,
            a: 1.0,
        },
        9 => DrawColour {
            r: 0.2,
            g: 0.8,
            b: 0.8,
            a: 1.0,
        },
        15 => DrawColour {
            r: 1.0,
            g: 0.5,
            b: 0.0,
            a: 1.0,
        },
        16 => DrawColour {
            r: 0.8,
            g: 0.8,
            b: 0.0,
            a: 1.0,
        },
        17 => DrawColour {
            r: 0.0,
            g: 0.802,
            b: 0.0,
            a: 1.0,
        },
        35 => DrawColour {
            r: 0.5,
            g: 0.3,
            b: 0.1,
            a: 1.0,
        },
        53 => DrawColour {
            r: 0.63,
            g: 0.12,
            b: 0.94,
            a: 1.0,
        },
        _ => DrawColour {
            r: 0.0,
            g: 0.0,
            b: 0.0,
            a: 1.0,
        },
    }
}

fn char_width(ch: char) -> f64 {
    const CHAR_WIDTHS: &[u16] = &[
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 278, 278, 355, 556, 556, 889, 667, 222, 333, 333, 389, 584, 278, 333, 278, 278, 556,
        556, 556, 556, 556, 556, 556, 556, 556, 556, 278, 278, 584, 584, 584, 556, 1015, 667, 667,
        722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833, 722, 778, 667, 778, 722, 667, 611,
        722, 667, 944, 667, 667, 611, 278, 278, 278, 469, 556, 222, 556, 556, 500, 556, 556, 278,
        556, 556, 222, 222, 500, 222, 833, 556, 556, 556, 556, 333, 500, 278, 556, 500, 722, 500,
        500, 500, 334, 260, 334, 584, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 333, 556, 556, 167, 556, 556, 556, 556, 191, 333,
        556, 333, 333, 500, 500, 0, 556, 556, 556, 278, 0, 537, 350, 222, 333, 333, 556, 1000,
        1000, 0, 611, 0, 333, 333, 333, 333, 333, 333, 333, 333, 0, 333, 333, 0, 333, 333, 333,
        1000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1000, 0, 370, 0, 0, 0, 0, 556, 778,
        1000, 365, 0, 0, 0, 0, 0, 889, 0, 0, 0, 278, 0, 0, 222, 611, 944, 611, 0, 0, 834,
    ];
    let idx = ch as usize;
    if idx < CHAR_WIDTHS.len() && CHAR_WIDTHS[idx] != 0 {
        f64::from(CHAR_WIDTHS[idx])
    } else {
        556.0
    }
}

fn xml_escape(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

fn draw_colour_to_svg(colour: DrawColour) -> String {
    fn component(v: f64) -> u8 {
        (255.0 * v) as u8
    }
    if (1.0 - colour.a) > 1e-3 {
        format!(
            "#{:02X}{:02X}{:02X}{:02X}",
            component(colour.r),
            component(colour.g),
            component(colour.b),
            component(colour.a)
        )
    } else {
        format!(
            "#{:02X}{:02X}{:02X}",
            component(colour.r),
            component(colour.g),
            component(colour.b)
        )
    }
}

fn format_alpha(value: f64) -> String {
    if (value - value.round()).abs() < 1e-9 {
        format!("{:.0}", value)
    } else {
        format_double(value)
    }
}

fn format_double(value: f64) -> String {
    format!("{value:.1}")
}
