//! Source-shaped REACCS state and MolBlock parser used by Avalon.
//!
//! This module is intentionally private to the Avalon implementation. RDKit's
//! molecule adapter serializes through `MolToMolBlock(mol, true)` and the
//! pinned Avalon engine parses that text into these records before computing
//! any fingerprint feature.

use crate::chemistry::valence::rdkit_element_symbol;
use crate::io::molblock::{MolBlockWriteParams, SdfFormat, mol_to_mol_block_with_params};
use crate::{FingerprintError, Molecule};

const MAX_BUFFER: usize = 4001;
const MAX_ATOMS: usize = 999;
const MAX_BONDS: usize = 999;
const MDL_MAX_LINE: usize = 200;
const MAX_SYMBOL_LIST: usize = 80;
const V30_ATOM_SYMBOL_CAPACITY: usize = 100;
const V30_PROPERTY_CAPACITY: usize = 20;

#[derive(Debug, Clone, PartialEq)]
pub(super) struct Atom {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub atom_symbol: String,
    pub mass_difference: i32,
    pub charge: i32,
    pub radical: i32,
    pub stereo_parity: i32,
    pub query_h_count: i32,
    pub query_stereo_box: i32,
    pub dummy1: i32,
    pub dummy2: i32,
    pub dummy3: i32,
    pub dummy4: i32,
    pub mapping: i32,
    pub second_stereo_parity: i32,
    pub dummy5: i32,
    pub sub_desc: i32,
    pub value: f32,
    pub rsize_flags: u32,
    pub color: i32,
    pub atext: String,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct Bond {
    pub atoms: [i32; 2],
    pub bond_type: i32,
    pub stereo_symbol: i32,
    pub dummy: i32,
    pub rsize_flags: u32,
    pub topography: i32,
    pub reaction_mark: i32,
    pub value: f32,
    pub bond_type_flags: i32,
    pub color: i32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct SymbolList {
    pub atom: i32,
    pub inclusive: bool,
    pub symbols: String,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct STextLine {
    pub x: f32,
    pub y: f32,
    pub text: String,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct MoleculeState {
    pub name: String,
    pub user_initials: String,
    pub program_name: String,
    pub date: String,
    pub time: String,
    pub dimensionality: String,
    pub scale1: i32,
    pub scale2: f32,
    pub energy: f32,
    pub registry_number: i64,
    pub comment: String,
    pub dummy1: i32,
    pub chiral_flag: i32,
    pub version: String,
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub symbol_lists: Vec<SymbolList>,
    pub stext_lines: Vec<STextLine>,
    pub prop_lines: Vec<String>,
    pub color: i32,
}

impl Default for Atom {
    fn default() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            atom_symbol: String::new(),
            mass_difference: 0,
            charge: 0,
            radical: 0,
            stereo_parity: 0,
            query_h_count: 0,
            query_stereo_box: 0,
            dummy1: 0,
            dummy2: 0,
            dummy3: 0,
            dummy4: 0,
            mapping: 0,
            second_stereo_parity: 0,
            dummy5: 0,
            sub_desc: 0,
            value: 0.0,
            rsize_flags: 0,
            color: 0,
            atext: String::new(),
        }
    }
}

impl Default for Bond {
    fn default() -> Self {
        Self {
            atoms: [0; 2],
            bond_type: 0,
            stereo_symbol: 0,
            dummy: 0,
            rsize_flags: 0,
            topography: 0,
            reaction_mark: 0,
            value: 0.0,
            bond_type_flags: 0,
            color: 0,
        }
    }
}

impl Default for MoleculeState {
    fn default() -> Self {
        Self {
            name: String::new(),
            user_initials: String::new(),
            program_name: "DUMMY".to_string(),
            date: String::new(),
            time: String::new(),
            dimensionality: String::new(),
            scale1: 0,
            scale2: 0.0,
            energy: 0.0,
            registry_number: 0,
            comment: String::new(),
            dummy1: 0,
            chiral_flag: 0,
            version: String::new(),
            atoms: Vec::new(),
            bonds: Vec::new(),
            symbol_lists: Vec::new(),
            stext_lines: Vec::new(),
            prop_lines: Vec::new(),
            color: 0,
        }
    }
}

fn conversion_error(reason: impl Into<String>) -> FingerprintError {
    FingerprintError::AvalonConversion {
        reason: reason.into(),
    }
}

/// In-memory equivalent of Avalon's `Fortran_FILE` string path.
struct FortranString<'a> {
    source: &'a str,
    offset: usize,
    line_no: usize,
    current: Option<&'a str>,
}

impl<'a> FortranString<'a> {
    fn open(source: &'a str) -> Self {
        // Avalon❗🔝: Fortran_FILE *FortranStringOpen(char *string_buffer)
        // Avalon❗🔝: {
        // Avalon❗🔝:    int i;
        // Avalon❗🔝:    Fortran_FILE *fp;
        // Avalon❗🔝:    for (i=0, fp=fortran_files; i<NFILES; i++, fp++)
        // Avalon❗🔝:       if (!fortran_files[i].in_use) break;
        // Avalon❗🔝:    if (i == NFILES) return ((Fortran_FILE *)NULL);
        // Avalon❗🔝:    fp->filep = (FILE *)NULL;
        // Avalon❗🔝:    fp->source_string  = TypeAlloc(strlen(string_buffer)+1, char);
        // Avalon❗🔝:    fp->string_pointer = fp->source_string;
        // Avalon❗🔝:    if (!fp->source_string) return((Fortran_FILE *)NULL);
        // Avalon❗🔝:    strcpy(fp->source_string, string_buffer);
        // Avalon❗🔝:    fortran_files[i].status    = FORTRAN_NORMAL;
        // Avalon❗🔝:    fortran_files[i].in_use    = TRUE;
        // Avalon❗🔝:    fortran_files[i].line_no   = 0;
        // Avalon❗🔝:    fortran_files[i].buffer[0] = '\0';
        // Avalon❗🔝:    GetBuffer(&fortran_files[i]);
        // Avalon❗🔝:    return(&fortran_files[i]);
        // Avalon❗🔝: }
        // Borrowing the caller's immutable string removes the source `strcpy`
        // allocation while retaining the same cursor and line-buffer behavior.
        let mut stream = Self {
            source,
            offset: 0,
            line_no: 0,
            current: None,
        };
        stream.advance();
        stream
    }

    fn line(&self) -> Result<&'a str, FingerprintError> {
        self.current
            .ok_or_else(|| conversion_error("unexpected end of MolBlock"))
    }

    fn advance(&mut self) {
        // Avalon❗🔝: void GetBuffer(Fortran_FILE *fp)
        // Avalon❗🔝: {
        // Avalon❗🔝:    char *cp, *cph;
        // Avalon❗🔝:    if (fp->source_string)
        // Avalon❗🔝:    {
        // Avalon❗🔝:       if (fp->string_pointer[0] == '\0')
        // Avalon❗🔝:       { fp->buffer[0] = '\0'; fp->status = FORTRAN_EOF; }
        // Avalon❗🔝:       else
        // Avalon❗🔝:       {
        // Avalon❗🔝:          for (cp=fp->string_pointer, cph=fp->buffer;
        // Avalon❗🔝:               (*cp) && (*cp) != '\n' && (*cp) != '\r' &&
        // Avalon❗🔝:                  cp < fp->string_pointer + MAX_BUFFER;
        // Avalon❗🔝:               cp++)
        // Avalon❗🔝:          { (*cph) = (*cp); cph++; }
        // Avalon❗🔝:          if ((*cp) == '\r') cp++;
        // Avalon❗🔝:          if ((*cp) == '\n') cp++;
        // Avalon❗🔝:          fp->string_pointer = cp;
        // Avalon❗🔝:          (*cph) = '\0';
        // Avalon❗🔝:       }
        // Avalon❗🔝:       RemoveTrailingBlanks(fp->buffer);
        // Avalon❗🔝:       fp->line_no++;
        // Avalon❗🔝:       return;
        // Avalon❗🔝:    }
        if self.offset >= self.source.len() {
            self.current = None;
            return;
        }
        let bytes = self.source.as_bytes();
        let start = self.offset;
        let mut end = start;
        while end < bytes.len()
            && end - start < MAX_BUFFER
            && bytes[end] != b'\r'
            && bytes[end] != b'\n'
        {
            end += 1;
        }
        self.offset = end;
        if self.offset < bytes.len() && bytes[self.offset] == b'\r' {
            self.offset += 1;
        }
        if self.offset < bytes.len() && bytes[self.offset] == b'\n' {
            self.offset += 1;
        }
        self.line_no += 1;
        self.current = self.source.get(start..end).map(str::trim_end);
    }
}

fn field(line: &str, start: usize, end: usize) -> &str {
    if start >= line.len() {
        ""
    } else {
        line.get(start..end.min(line.len())).unwrap_or("")
    }
}

fn parse_i32(text: &str, label: &str) -> Result<i32, FingerprintError> {
    text.trim()
        .parse()
        .map_err(|_| conversion_error(format!("invalid {label}")))
}

fn parse_i32_or_zero(text: &str) -> i32 {
    text.trim().parse().unwrap_or(0)
}

fn parse_f32(text: &str, label: &str) -> Result<f32, FingerprintError> {
    text.trim()
        .parse()
        .map_err(|_| conversion_error(format!("invalid {label}")))
}

fn scan_decimal_field<'a>(
    text: &'a str,
    cursor: &mut usize,
    width: usize,
    fractional: bool,
) -> Option<&'a str> {
    let bytes = text.as_bytes();
    while *cursor < bytes.len() && bytes[*cursor].is_ascii_whitespace() {
        *cursor += 1;
    }
    let start = *cursor;
    let limit = start.saturating_add(width).min(bytes.len());
    let mut index = start;
    if index < limit && matches!(bytes[index], b'+' | b'-') {
        index += 1;
    }
    let integer_start = index;
    while index < limit && bytes[index].is_ascii_digit() {
        index += 1;
    }
    let mut has_digit = index > integer_start;
    if fractional && index < limit && bytes[index] == b'.' {
        index += 1;
        let fraction_start = index;
        while index < limit && bytes[index].is_ascii_digit() {
            index += 1;
        }
        has_digit |= index > fraction_start;
    }
    if !has_digit {
        return None;
    }
    if fractional && index < limit && matches!(bytes[index], b'e' | b'E') {
        let exponent = index;
        index += 1;
        if index < limit && matches!(bytes[index], b'+' | b'-') {
            index += 1;
        }
        let exponent_digits = index;
        while index < limit && bytes[index].is_ascii_digit() {
            index += 1;
        }
        if index == exponent_digits {
            index = exponent;
        }
    }
    *cursor = index;
    text.get(start..index)
}

fn split_charge_radical(charge_radical: i32) -> (i32, i32) {
    // Avalon✔️✔️: void SplitChargeRadical(int charge_radical, int *chargep, int *radicalp)
    // Avalon✔️✔️: {
    // Avalon✔️✔️:    if (charge_radical == NONE)
    // Avalon✔️✔️:    { *chargep = NONE; *radicalp = NONE; }
    // Avalon✔️✔️:    else if (charge_radical == RADICAL)
    // Avalon✔️✔️:    { *chargep = NONE; *radicalp = DOUBLET; }
    // Avalon✔️✔️:    else
    // Avalon✔️✔️:    { *chargep = RADICAL-charge_radical; *radicalp = NONE; }
    // Avalon✔️✔️: }
    match charge_radical {
        0 => (0, 0),
        4 => (0, 2),
        value => (4 - value, 0),
    }
}

fn read_v2000_atom(stream: &mut FortranString<'_>) -> Result<Atom, FingerprintError> {
    // Avalon❗✔️: int ReadREACCSAtom(Fortran_FILE *fp, struct reaccs_atom_t *ap)
    // Avalon❗✔️: {
    // Avalon❗✔️:    int nitems;
    // Avalon❗✔️:    char buffer[MAX_BUFFER+1];
    // Avalon❗✔️:    int charge_radical;
    // Avalon❗✔️:
    // Avalon❗✔️:    if (fp->status != FORTRAN_NORMAL) return(fp->status);
    // Avalon❗✔️:
    // Avalon❗✔️:    strncpy(buffer,fp->buffer,MAX_BUFFER);
    // Avalon❗✔️:
    // Avalon❗✔️:    if (buffer[31] == ' ') buffer[32] = '?'; /* to cope with empty atom types */
    // Avalon❗✔️:    nitems = sscanf(buffer,
    // Avalon❗✔️:                   "%10f%10f%10f %s",
    // Avalon❗✔️:                   &ap->x,  &ap->y,  &ap->z,
    // Avalon❗✔️:                    ap->atom_symbol);
    // Avalon❗✔️:    BlankToZero(buffer+34);
    // Avalon❗✔️:    ap->mass_difference = 0; charge_radical = 0;
    // Avalon❗✔️:    ap->stereo_parity = NONE;
    // Avalon❗✔️:    ap->query_H_count = 0;
    // Avalon❗✔️:    ap->query_stereo_box = 0;
    // Avalon❗✔️:    ap->dummy1 = 0; ap->dummy2 = 0; ap->dummy3 = 0; ap->dummy4 = 0;
    // Avalon❗✔️:    ap->sub_desc = NONE;
    // Avalon❗✔️:    ap->mapping = NONE;
    // Avalon❗✔️:    ap->second_stereo_parity = NONE;
    // Avalon❗✔️:    ap->dummy5 = 0;
    // Avalon❗✔️:    ap->atext[0]='\0';
    // Avalon❗✔️:    sscanf(buffer+34,
    // Avalon❗✔️:       "%2d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d",
    // Avalon❗✔️:       &ap->mass_difference,
    // Avalon❗✔️:       &charge_radical, &ap->stereo_parity,
    // Avalon❗✔️:       &ap->query_H_count,  &ap->query_stereo_box,
    // Avalon❗✔️:       &ap->dummy1, &ap->dummy2, &ap->dummy3, &ap->dummy4,
    // Avalon❗✔️:       &ap->mapping,  &ap->second_stereo_parity,
    // Avalon❗✔️:       &ap->dummy5);
    // Avalon❗✔️:
    // Avalon❗✔️:    if (nitems >= 4)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       SplitChargeRadical(charge_radical, &ap->charge, &ap->radical);
    // Avalon❗✔️:       GetBuffer(fp);
    // Avalon❗✔️:       return(FORTRAN_NORMAL);
    // Avalon❗✔️:    }
    // Avalon❗✔️:    else
    // Avalon❗✔️:    {
    // Avalon❗✔️:       ShowMessageI("incorrect # (%d) of arguments on atom line",
    // Avalon❗✔️:                    "ReadREACCSAtom:",
    // Avalon❗✔️:                    nitems);
    // Avalon❗✔️:       ShowMessageS("buffer ==\n%s\n","ReadREACCSAtom",fp->buffer);
    // Avalon❗✔️:       GetBuffer(fp);
    // Avalon❗✔️:       return(FORTRAN_ERROR);
    // Avalon❗✔️:    }
    // Avalon❗✔️: }
    let line = stream.line()?;
    let mut atom = Atom {
        x: parse_f32(field(line, 0, 10), "V2000 atom x")?,
        y: parse_f32(field(line, 10, 20), "V2000 atom y")?,
        z: parse_f32(field(line, 20, 30), "V2000 atom z")?,
        atom_symbol: field(line, 31, 34).trim().to_string(),
        ..Atom::default()
    };
    if atom.atom_symbol.is_empty() {
        atom.atom_symbol = "?".to_string();
    }
    atom.mass_difference = parse_i32_or_zero(field(line, 34, 36));
    let charge_radical = parse_i32_or_zero(field(line, 36, 39));
    atom.stereo_parity = parse_i32_or_zero(field(line, 39, 42));
    atom.query_h_count = parse_i32_or_zero(field(line, 42, 45));
    atom.query_stereo_box = parse_i32_or_zero(field(line, 45, 48));
    atom.dummy1 = parse_i32_or_zero(field(line, 48, 51));
    atom.dummy2 = parse_i32_or_zero(field(line, 51, 54));
    atom.dummy3 = parse_i32_or_zero(field(line, 54, 57));
    atom.dummy4 = parse_i32_or_zero(field(line, 57, 60));
    atom.mapping = parse_i32_or_zero(field(line, 60, 63));
    atom.second_stereo_parity = parse_i32_or_zero(field(line, 63, 66));
    atom.dummy5 = parse_i32_or_zero(field(line, 66, 69));
    (atom.charge, atom.radical) = split_charge_radical(charge_radical);
    stream.advance();
    Ok(atom)
}

fn read_v2000_bond(stream: &mut FortranString<'_>) -> Result<Bond, FingerprintError> {
    // Avalon❗✔️: int ReadREACCSBond(Fortran_FILE *fp, struct reaccs_bond_t *bp)
    // Avalon❗✔️: {
    // Avalon❗✔️:    int nitems, i, j, k;
    // Avalon❗✔️:    int bond_line_len, n_chars, pos;
    // Avalon❗✔️:    int *ptrarray[MAX_BONDLINE_FIELDS];
    // Avalon❗✔️:    char c;
    // Avalon❗✔️:    char buffer[BONDLINE_FIELD_LEN+1];
    // Avalon❗✔️:
    // Avalon❗✔️:    if (fp->status != FORTRAN_NORMAL) return(fp->status);
    // Avalon❗✔️:
    // Avalon❗✔️:    bp->stereo_symbol = 0;
    // Avalon❗✔️:    bp->dummy = 0;
    // Avalon❗✔️:    bp->topography = 0;
    // Avalon❗✔️:    bp->reaction_mark = NONE;
    // Avalon❗✔️:    ptrarray[0] = &bp->atoms[0];
    // Avalon❗✔️:    ptrarray[1] = &bp->atoms[1];
    // Avalon❗✔️:    ptrarray[2] = &bp->bond_type;
    // Avalon❗✔️:    ptrarray[3] = &bp->stereo_symbol;
    // Avalon❗✔️:    ptrarray[4] = &bp->dummy;
    // Avalon❗✔️:    ptrarray[5] = &bp->topography;
    // Avalon❗✔️:    ptrarray[6] = &bp->reaction_mark;
    // Avalon❗✔️:    bond_line_len = strlen(fp->buffer);
    // Avalon❗✔️:    nitems = bond_line_len ? (bond_line_len - 1) / BONDLINE_FIELD_LEN + 1 : 0;
    // Avalon❗✔️:    if (nitems > MAX_BONDLINE_FIELDS)
    // Avalon❗✔️:       nitems = MAX_BONDLINE_FIELDS;
    // Avalon❗✔️:    for (i = 0; i < nitems; ++i)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       pos = i * BONDLINE_FIELD_LEN;
    // Avalon❗✔️:       memset(buffer, 0, BONDLINE_FIELD_LEN + 1);
    // Avalon❗✔️:       n_chars = bond_line_len - pos;
    // Avalon❗✔️:       if (n_chars > BONDLINE_FIELD_LEN)
    // Avalon❗✔️:          n_chars = BONDLINE_FIELD_LEN;
    // Avalon❗✔️:       for (j = 0, k = 0; j < n_chars; ++j)
    // Avalon❗✔️:       {
    // Avalon❗✔️:          c = fp->buffer[pos + j];
    // Avalon❗✔️:          if (c != ' ')
    // Avalon❗✔️:             buffer[k++] = c;
    // Avalon❗✔️:       }
    // Avalon❗✔️:       sscanf(buffer, "%3d", ptrarray[i]);
    // Avalon❗✔️:    }
    // Avalon❗✔️:    if (nitems >= 3)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       GetBuffer(fp);
    // Avalon❗✔️:       return(FORTRAN_NORMAL);
    // Avalon❗✔️:    }
    // Avalon❗✔️:    else
    // Avalon❗✔️:    {
    // Avalon❗✔️:       ShowMessageI("incorrect # (%d) of arguments on bond line",
    // Avalon❗✔️:                    "ReadREACCSBond",
    // Avalon❗✔️:                    nitems);
    // Avalon❗✔️:       ShowMessageS("buffer ==\n%s\n","ReadREACCSBond",buffer);
    // Avalon❗✔️:       return(FORTRAN_ERROR);
    // Avalon❗✔️:    }
    // Avalon❗✔️: }
    let line = stream.line()?;
    let nitems = line.len().div_ceil(3).min(7);
    if nitems < 3 {
        return Err(conversion_error(
            "V2000 bond line has fewer than three fields",
        ));
    }
    let mut values = [0_i32; 7];
    for (index, value) in values.iter_mut().enumerate().take(nitems) {
        *value = parse_i32_or_zero(field(line, index * 3, index * 3 + 3));
    }
    stream.advance();
    Ok(Bond {
        atoms: [values[0], values[1]],
        bond_type: values[2],
        stereo_symbol: values[3],
        dummy: values[4],
        topography: values[5],
        reaction_mark: values[6],
        ..Bond::default()
    })
}

fn lookup_mass_difference(mass: f64, symbol: &str) -> i32 {
    // Avalon❗✔️: int LookupMassDifference(double mass, char *symbol)
    // Avalon❗✔️: {
    // Avalon❗✔️:    int i;
    // Avalon❗✔️:    for (i=0; i<sizeof(ptable); i++)
    // Avalon❗✔️:       if (0 == strcmp(ptable[i].symbol, symbol))
    // Avalon❗✔️:       {
    // Avalon❗✔️:          double diff = mass-ptable[i].mass;
    // Avalon❗✔️:          if (diff > 0) return (int)(diff+0.49);
    // Avalon❗✔️:          else return (int)(diff-0.49);
    // Avalon❗✔️:       }
    // Avalon❗✔️:    if (mass > 0) return (int)(mass+0.49);
    // Avalon❗✔️:    else return (int)(mass-0.49);
    // Avalon❗✔️: }
    // The source table uses historical Avalon atomic weights rather than the
    // pinned RDKit table. Keep it local so MASS conversion cannot drift when
    // another periodic-table consumer changes.
    // Avalon❗✔️: struct ptable_entry
    // Avalon❗✔️:    {
    // Avalon❗✔️:       char *symbol;             /* atomic symbol */
    // Avalon❗✔️:       int  valence;             /* common valence of element */
    // Avalon❗✔️:       float mass;               /* average mass of atoms */
    // Avalon❗✔️:       float electronegativity;  /* el.neg. according to Pauling */
    // Avalon❗✔️:       float heat_of_formation;  /* h.o.f. of element from atoms */
    // Avalon❗✔️:    };
    const WEIGHTS: &[(&str, f32)] = &[
        ("C", 12.01115),
        ("N", 14.00670),
        ("O", 15.99940),
        ("H", 1.00797),
        ("He", 4.00260),
        ("Li", 6.93900),
        ("Be", 9.01220),
        ("B", 10.81100),
        ("F", 18.99840),
        ("Ne", 20.18300),
        ("Na", 22.98980),
        ("Mg", 24.31200),
        ("Al", 26.98150),
        ("Si", 28.08600),
        ("P", 30.97380),
        ("S", 32.06400),
        ("Cl", 35.45300),
        ("Ar", 39.94800),
        ("K", 39.10200),
        ("Ca", 40.08000),
        ("Sc", 44.95600),
        ("Ti", 47.90000),
        ("V", 50.94200),
        ("Cr", 51.99600),
        ("Mn", 54.93800),
        ("Fe", 55.84700),
        ("Co", 58.93320),
        ("Ni", 58.71000),
        ("Cu", 63.54600),
        ("Zn", 65.37000),
        ("Ga", 69.72000),
        ("Ge", 72.59000),
        ("As", 74.92160),
        ("Se", 78.96000),
        ("Br", 79.90400),
        ("Kr", 83.80000),
        ("Rb", 85.47000),
        ("Sr", 87.62000),
        ("Y", 88.90500),
        ("Zr", 91.22000),
        ("Nb", 92.90600),
        ("Mo", 95.94000),
        ("Tc", 98.90620),
        ("Ru", 101.07000),
        ("Rh", 102.90500),
        ("Pd", 106.40000),
        ("Ag", 107.86800),
        ("Cd", 112.40000),
        ("In", 114.82000),
        ("Sn", 118.69000),
        ("Sb", 121.75000),
        ("Te", 127.60000),
        ("I", 126.90440),
        ("Xe", 131.30000),
        ("Cs", 132.90500),
        ("Ba", 137.33000),
        ("La", 138.91000),
        ("Ce", 140.12000),
        ("Pr", 140.90700),
        ("Nd", 144.24000),
        ("Pm", 145.00000),
        ("Sm", 150.35000),
        ("Eu", 151.96000),
        ("Gd", 157.25000),
        ("Tb", 158.92400),
        ("Dy", 162.50000),
        ("Ho", 164.93000),
        ("Er", 167.26000),
        ("Tm", 168.93400),
        ("Yb", 173.04000),
        ("Lu", 174.97000),
        ("Hf", 178.49000),
        ("Ta", 180.94800),
        ("W", 183.85000),
        ("Re", 186.20000),
        ("Os", 190.20000),
        ("Ir", 192.20000),
        ("Pt", 195.09000),
        ("Au", 196.96700),
        ("Hg", 200.59000),
        ("Tl", 204.37000),
        ("Pb", 207.19000),
        ("Bi", 208.98000),
        ("Po", 209.00000),
        ("At", 210.00000),
        ("Rn", 222.00000),
        ("Fr", 223.00000),
        ("Ra", 226.03000),
        ("Ac", 227.00000),
        ("Th", 232.03800),
        ("Pa", 231.04000),
        ("U", 238.03000),
        ("Np", 237.05000),
        ("Pu", 244.00000),
        ("Am", 243.00000),
        ("Cm", 247.00000),
        ("Bk", 247.00000),
        ("Cf", 251.00000),
        ("Es", 254.00000),
        ("Fm", 257.00000),
        ("Md", 258.00000),
        ("No", 259.00000),
        ("Lr", 260.00000),
        ("D", 2.01400),
        ("R", 0.00000),
    ];
    let reference = WEIGHTS
        .iter()
        .find_map(|(name, value)| (*name == symbol).then_some(f64::from(*value)));
    let diff = reference.map_or(mass, |value| mass - value);
    if diff > 0.0 {
        (diff + 0.49) as i32
    } else {
        (diff - 0.49) as i32
    }
}

fn parse_v30_symbol(
    symbol: &str,
    atom_index: usize,
    state: &mut MoleculeState,
) -> Result<String, FingerprintError> {
    // Avalon✔️✔️: struct symbol_list_t *ParseV30SymbolList(char *symbol,
    // Avalon✔️✔️:                                         int iatom,
    // Avalon✔️✔️:                                         struct reaccs_molecule_t *mp,
    // Avalon✔️✔️:                                         struct symbol_list_t *old_list)
    // Avalon✔️✔️: {
    // Avalon✔️✔️:    bra_from = strchr(symbol, '[');
    // Avalon✔️✔️:    bra_to = strchr(symbol, ']');
    // Avalon✔️✔️:    if (bra_from == NULL || bra_to == NULL || bra_to < bra_from)
    // Avalon✔️✔️:    { list = old_list; strcpy(mp->atom_array[iatom].atom_symbol, "Unk"); }
    // Avalon✔️✔️:    else
    // Avalon✔️✔️:    {
    // Avalon✔️✔️:       len = (bra_to-bra_from)-1;
    // Avalon✔️✔️:       list = TypeAlloc(1,struct symbol_list_t);
    // Avalon✔️✔️:       list->next = old_list; list->atom = iatom+1;
    // Avalon✔️✔️:       strncpy(list->string, bra_from+1, len);
    // Avalon✔️✔️:       list->string[len] = '\0';
    // Avalon✔️✔️:       if (0 == strncmp(symbol,"NOT",3)) list->logic = EXCLUSIVE;
    // Avalon✔️✔️:       else list->logic = INCLUSIVE;
    // Avalon✔️✔️:       strcpy(mp->atom_array[iatom].atom_symbol, "L");
    // Avalon✔️✔️:    }
    // Avalon✔️✔️:    return list;
    // Avalon✔️✔️: }
    let Some(open) = symbol.find('[') else {
        return Ok("Unk".to_string());
    };
    let Some(close) = symbol.find(']') else {
        return Ok("Unk".to_string());
    };
    if close < open {
        return Ok("Unk".to_string());
    }
    if close - open - 1 > MAX_SYMBOL_LIST {
        return Err(conversion_error(
            "V3000 atom list exceeds Avalon's symbol-list capacity",
        ));
    }
    state.symbol_lists.push(SymbolList {
        atom: i32::try_from(atom_index + 1).unwrap_or(i32::MAX),
        inclusive: !symbol.starts_with("NOT"),
        symbols: symbol[open + 1..close].to_string(),
    });
    Ok("L".to_string())
}

fn parse_v30_atom(
    line: &str,
    index: usize,
    state: &mut MoleculeState,
) -> Result<Atom, FingerprintError> {
    // Avalon❗✔️: int ReadV30Atom(Fortran_FILE *fp,
    // Avalon❗✔️:                 struct reaccs_atom_t *ap,
    // Avalon❗✔️:                 struct reaccs_molecule_t *mp)
    // Avalon❗✔️: {
    // Avalon❗✔️:    int nitems;
    // Avalon❗✔️:    char buffer[MAX_BUFFER+1];
    // Avalon❗✔️:    int charge_radical;
    // Avalon❗✔️:    char prop[10][20];
    // Avalon❗✔️:    int seq;
    // Avalon❗✔️:    int idummy;
    // Avalon❗✔️:    float mass;
    // Avalon❗✔️:    int i;
    // Avalon❗✔️:    char symbol[100];
    // Avalon❗✔️:    struct symbol_list_t *list;
    // Avalon❗✔️:
    // Avalon❗✔️:    if (fp->status != FORTRAN_NORMAL) return(fp->status);
    // Avalon❗✔️:
    // Avalon❗✔️:    strncpy(buffer,fp->buffer,MAX_BUFFER);
    // Avalon❗✔️:
    // Avalon❗✔️:    nitems = sscanf(buffer+strlen("M  V30"),
    // Avalon❗✔️:                   "%d %s %f %f %f %d %s %s %s %s %s %s %s %s %s %s",
    // Avalon❗✔️:                   &seq, symbol,
    // Avalon❗✔️:                   &ap->x,  &ap->y,  &ap->z,
    // Avalon❗✔️:                   &ap->mapping,
    // Avalon❗✔️:                   prop[0], prop[1], prop[2], prop[3], prop[4],
    // Avalon❗✔️:                   prop[5], prop[6], prop[7], prop[8], prop[9]);
    // Avalon❗✔️:
    // Avalon❗✔️:    if (nitems < 6)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       ShowMessageI("incorrect # (%d) of arguments on atom line",
    // Avalon❗✔️:                    "ReadV30Atom:",
    // Avalon❗✔️:                    nitems);
    // Avalon❗✔️:       ShowMessageS("buffer ==\n%s\n","ReadV30Atom",fp->buffer);
    // Avalon❗✔️:       GetBuffer(fp);
    // Avalon❗✔️:       return(FORTRAN_ERROR);
    // Avalon❗✔️:    }
    // Avalon❗✔️:
    // Avalon❗✔️:    if (strlen(symbol) <= 3)
    // Avalon❗✔️:       strcpy(ap->atom_symbol, symbol);
    // Avalon❗✔️:    else     // need to parse symbol
    // Avalon❗✔️:    {
    // Avalon❗✔️:       mp->symbol_lists =
    // Avalon❗✔️:           ParseV30SymbolList(symbol, ap-mp->atom_array, mp, mp->symbol_lists);
    // Avalon❗✔️:       mp->n_atom_lists = 0;
    // Avalon❗✔️:       for (list=mp->symbol_lists; !IsNULL(list); list = list->next) mp->n_atom_lists++;
    // Avalon❗✔️:    }
    // Avalon❗✔️:
    // Avalon❗✔️:    if (nitems >= 16)
    // Avalon❗✔️:       ShowMessageI("Warning %d properties for atom.\n",
    // Avalon❗✔️:                    "ReadV30Atom:", (nitems-6));
    // Avalon❗✔️:
    // Avalon❗✔️:    /* Processing properties */
    // Avalon❗✔️:    for (i=0; i<nitems-6; i++)
    // Avalon❗✔️: {
    // Avalon❗✔️:       if (1 == sscanf(prop[i], "CHG=%d", &ap->charge)) continue;
    // Avalon❗✔️:       if (1 == sscanf(prop[i], "RAD=%d", &ap->radical)) continue;
    // Avalon❗✔️:       if (1 == sscanf(prop[i], "HCOUNT=%d", &ap->query_H_count)) continue;
    // Avalon❗✔️:       if (1 == sscanf(prop[i], "CFG=%d", &ap->stereo_parity)) continue;
    // Avalon❗✔️:       if (1 == sscanf(prop[i], "VAL=%d", &ap->dummy1)) continue;
    // Avalon❗✔️:       if (1 == sscanf(prop[i], "SUBST=%d", &ap->sub_desc)) continue;
    // Avalon❗✔️:       if (1 == sscanf(prop[i], "MASS=%f", &mass))
    // Avalon❗✔️:       {
    // Avalon❗✔️:          ap->mass_difference = LookupMassDifference(mass, ap->atom_symbol);
    // Avalon❗✔️:          continue;
    // Avalon❗✔️:       }
    // Avalon❗✔️:       if (1 == sscanf(prop[i], "RBCNT=%d", &idummy)) continue;
    // Avalon❗✔️:       if (1 == sscanf(prop[i], "UNSAT=%d", &idummy)) continue;
    // Avalon❗✔️:       fprintf(stderr, "ignoring '%s' for atom %d\n", prop[i], seq);
    // Avalon❗✔️: }
    // Avalon❗✔️:    GetBuffer(fp);
    // Avalon❗✔️:    return(FORTRAN_NORMAL);
    // Avalon❗✔️: }
    let fields: Vec<&str> = line
        .strip_prefix("M  V30 ")
        .unwrap_or("")
        .split_whitespace()
        .collect();
    if fields.len() < 6 {
        return Err(conversion_error(
            "incorrect number of arguments on V3000 atom line",
        ));
    }
    let symbol = fields[1];
    if symbol.len() >= V30_ATOM_SYMBOL_CAPACITY {
        return Err(conversion_error(
            "V3000 atom symbol exceeds Avalon's source buffer capacity",
        ));
    }
    if fields
        .iter()
        .skip(6)
        .take(10)
        .any(|property| property.len() >= V30_PROPERTY_CAPACITY)
    {
        return Err(conversion_error(
            "V3000 atom property exceeds Avalon's source buffer capacity",
        ));
    }
    let atom_symbol = if symbol.len() <= 3 {
        symbol.to_string()
    } else {
        parse_v30_symbol(symbol, index, state)?
    };
    let mut atom = Atom {
        x: parse_f32(fields[2], "V3000 atom x")?,
        y: parse_f32(fields[3], "V3000 atom y")?,
        z: parse_f32(fields[4], "V3000 atom z")?,
        atom_symbol,
        mapping: parse_i32(fields[5], "V3000 atom mapping")?,
        ..Atom::default()
    };
    for property in fields.iter().skip(6).take(10) {
        if let Some(value) = property.strip_prefix("CHG=") {
            atom.charge = parse_i32(value, "V3000 atom charge")?;
        } else if let Some(value) = property.strip_prefix("RAD=") {
            atom.radical = parse_i32(value, "V3000 atom radical")?;
        } else if let Some(value) = property.strip_prefix("HCOUNT=") {
            atom.query_h_count = parse_i32(value, "V3000 atom HCOUNT")?;
        } else if let Some(value) = property.strip_prefix("CFG=") {
            atom.stereo_parity = parse_i32(value, "V3000 atom CFG")?;
        } else if let Some(value) = property.strip_prefix("VAL=") {
            atom.dummy1 = parse_i32(value, "V3000 atom VAL")?;
        } else if let Some(value) = property.strip_prefix("SUBST=") {
            atom.sub_desc = parse_i32(value, "V3000 atom SUBST")?;
        } else if let Some(value) = property.strip_prefix("MASS=") {
            // `%f` stores into Avalon's `float mass`; the function call then
            // performs C's exact float-to-double promotion.
            let mass = parse_f32(value, "V3000 atom MASS")?;
            atom.mass_difference = lookup_mass_difference(f64::from(mass), &atom.atom_symbol);
        }
    }
    Ok(atom)
}

fn parse_v30_bond(line: &str) -> Result<Bond, FingerprintError> {
    // Avalon❗✔️: int ReadV30Bond(Fortran_FILE *fp, struct reaccs_bond_t *bp)
    // Avalon❗✔️: {
    // Avalon❗✔️:    int nitems;
    // Avalon❗✔️:    char buffer[MAX_BUFFER+1];
    // Avalon❗✔️:    char prop[10][20];
    // Avalon❗✔️:    int seq;
    // Avalon❗✔️:    int i;
    // Avalon❗✔️:    char symbol[10];
    // Avalon❗✔️:
    // Avalon❗✔️:    if (fp->status != FORTRAN_NORMAL) return(fp->status);
    // Avalon❗✔️:
    // Avalon❗✔️:    strncpy(buffer,fp->buffer,MAX_BUFFER);
    // Avalon❗✔️:
    // Avalon❗✔️:    nitems = sscanf(buffer+strlen("M  V30"),
    // Avalon❗✔️:                    "%d %d %d %d %s %s %s %s %s %s %s %s %s %s",
    // Avalon❗✔️:                    &seq, &bp->bond_type,
    // Avalon❗✔️:                    &bp->atoms[0],  &bp->atoms[1],
    // Avalon❗✔️:                    prop[0], prop[1], prop[2], prop[3], prop[4],
    // Avalon❗✔️:                    prop[5], prop[6], prop[7], prop[8], prop[9]);
    // Avalon❗✔️:
    // Avalon❗✔️:    if (nitems < 4)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       ShowMessageI("incorrect # (%d) of arguments on bond line",
    // Avalon❗✔️:                    "ReadV30Atom:",
    // Avalon❗✔️:                    nitems);
    // Avalon❗✔️:       ShowMessageS("buffer ==\n%s\n","ReadV30Atom",fp->buffer);
    // Avalon❗✔️:       GetBuffer(fp);
    // Avalon❗✔️:       return(FORTRAN_ERROR);
    // Avalon❗✔️:    }
    // Avalon❗✔️:
    // Avalon❗✔️:    if (nitems >= 14)
    // Avalon❗✔️:       ShowMessageI("Warning %d properties for bond.\n",
    // Avalon❗✔️:                    "ReadV30Bond:", (nitems-6));
    // Avalon❗✔️:
    // Avalon❗✔️:    /* Processing properties */
    // Avalon❗✔️:    for (i=0; i<nitems-6; i++)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       fprintf(stderr, "ignoring '%s' for bond %d\n", prop[i], seq);
    // Avalon❗✔️:    }
    // Avalon❗✔️:    GetBuffer(fp);
    // Avalon❗✔️:    return(FORTRAN_NORMAL);
    // Avalon❗✔️: }
    let fields: Vec<&str> = line
        .strip_prefix("M  V30 ")
        .unwrap_or("")
        .split_whitespace()
        .collect();
    if fields.len() < 4 {
        return Err(conversion_error(
            "incorrect number of arguments on V3000 bond line",
        ));
    }
    if fields
        .iter()
        .skip(4)
        .take(10)
        .any(|property| property.len() >= V30_PROPERTY_CAPACITY)
    {
        return Err(conversion_error(
            "V3000 bond property exceeds Avalon's source buffer capacity",
        ));
    }
    Ok(Bond {
        bond_type: parse_i32(fields[1], "V3000 bond type")?,
        atoms: [
            parse_i32(fields[2], "V3000 bond begin")?,
            parse_i32(fields[3], "V3000 bond end")?,
        ],
        ..Bond::default()
    })
}

fn property_pairs(line: &str, prefix: &str) -> Result<Vec<(usize, i32)>, FingerprintError> {
    let fields: Vec<&str> = line[prefix.len()..].split_whitespace().collect();
    let count = fields
        .first()
        .ok_or_else(|| conversion_error("property line has no entry count"))?
        .parse::<usize>()
        .map_err(|_| conversion_error("invalid property entry count"))?;
    if count > 8 || fields.len() < 1 + count * 2 {
        return Err(conversion_error(
            "property line entry count exceeds its source fields",
        ));
    }
    (0..count)
        .map(|index| {
            let atom = fields[1 + index * 2]
                .parse::<usize>()
                .map_err(|_| conversion_error("invalid property atom index"))?;
            let value = fields[2 + index * 2]
                .parse::<i32>()
                .map_err(|_| conversion_error("invalid property value"))?;
            Ok((atom, value))
        })
        .collect()
}

fn property_atom<'a>(
    atoms: &'a mut [Atom],
    index: usize,
) -> Result<&'a mut Atom, FingerprintError> {
    atoms
        .get_mut(index.wrapping_sub(1))
        .ok_or_else(|| conversion_error("property atom index is outside the atom table"))
}

fn mdl_property_text(line: &str) -> String {
    line.chars().take(MDL_MAX_LINE).collect()
}

fn apply_isotope(atom: &mut Atom, atom_index: usize, mass: i32, prop_lines: &mut Vec<String>) {
    // Avalon❗✔️:         /* Just handle mass difference for substitution point labels and some well-known other atoms, other labelings are captured as property lines */
    // Avalon❗✔️:         for (n=0; n<nentries; n++)
    // Avalon❗✔️:            if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "R"))
    // Avalon❗✔️:               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n];
    // Avalon❗✔️:            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "H")  &&    1 <= tmp_vals[n]  &&  tmp_vals[n] <=   3)
    // Avalon❗✔️:               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-1;
    // Avalon❗✔️:            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol,"Li")  &&    6 <= tmp_vals[n]  &&  tmp_vals[n] <=   7)
    // Avalon❗✔️:               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-7;
    // Avalon❗✔️:            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "B")  &&   10 <= tmp_vals[n]  &&  tmp_vals[n] <=  11)
    // Avalon❗✔️:               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-11;
    // Avalon❗✔️:            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "C")  &&   12 <= tmp_vals[n]  &&  tmp_vals[n] <=  14)
    // Avalon❗✔️:               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-12;
    // Avalon❗✔️:            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "N")  &&   13 <= tmp_vals[n]  &&  tmp_vals[n] <=  15)
    // Avalon❗✔️:               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-14;
    // Avalon❗✔️:            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "O")  &&   16 <= tmp_vals[n]  &&  tmp_vals[n] <=  18)
    // Avalon❗✔️:               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-16;
    // Avalon❗✔️:            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "F")  &&   18 <= tmp_vals[n]  &&  tmp_vals[n] <=  19)
    // Avalon❗✔️:               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-19;
    // Avalon❗✔️:            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "P")  &&   31 <= tmp_vals[n]  &&  tmp_vals[n] <=  33)
    // Avalon❗✔️:               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-31;
    // Avalon❗✔️:            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "S")  &&   32 <= tmp_vals[n]  &&  tmp_vals[n] <=  36)
    // Avalon❗✔️:               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-32;
    // Avalon❗✔️:            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol,"Cl")  &&   35 <= tmp_vals[n]  &&  tmp_vals[n] <=  37)
    // Avalon❗✔️:               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-35;
    // Avalon❗✔️:            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol,"Co")  &&   56 <= tmp_vals[n]  &&  tmp_vals[n] <=  61)
    // Avalon❗✔️:               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-59;
    // Avalon❗✔️:            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol,"Br")  &&   79 <= tmp_vals[n]  &&  tmp_vals[n] <=  81)
    // Avalon❗✔️:               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-80;
    // Avalon❗✔️:            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "Y")  &&   88 <= tmp_vals[n]  &&  tmp_vals[n] <=  90)
    // Avalon❗✔️:               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-89;
    // Avalon❗✔️:            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "I")  &&  123 <= tmp_vals[n]  &&  tmp_vals[n] <= 131)
    // Avalon❗✔️:               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-127;
    // Avalon❗✔️:            else { // make sure unknown isotopic labelings get represented as properties
    // Avalon❗✔️:                hp = TypeAlloc(1,struct prop_line_t);
    // Avalon❗✔️:                sprintf(hp->text, "M  ISO  1 %3d %3d", tmp_ats[n], tmp_vals[n]);
    // Avalon❗✔️:                hp->text[MDL_MAXLINE] = '\0';
    // Avalon❗✔️:                hp->next = result; result = hp;
    // Avalon❗✔️:            }
    let difference = match atom.atom_symbol.as_str() {
        "R" => Some(mass),
        "H" if (1..=3).contains(&mass) => Some(mass - 1),
        "Li" if (6..=7).contains(&mass) => Some(mass - 7),
        "B" if (10..=11).contains(&mass) => Some(mass - 11),
        "C" if (12..=14).contains(&mass) => Some(mass - 12),
        "N" if (13..=15).contains(&mass) => Some(mass - 14),
        "O" if (16..=18).contains(&mass) => Some(mass - 16),
        "F" if (18..=19).contains(&mass) => Some(mass - 19),
        "P" if (31..=33).contains(&mass) => Some(mass - 31),
        "S" if (32..=36).contains(&mass) => Some(mass - 32),
        "Cl" if (35..=37).contains(&mass) => Some(mass - 35),
        "Co" if (56..=61).contains(&mass) => Some(mass - 59),
        "Br" if (79..=81).contains(&mass) => Some(mass - 80),
        "Y" if (88..=90).contains(&mass) => Some(mass - 89),
        "I" if (123..=131).contains(&mass) => Some(mass - 127),
        _ => None,
    };
    if let Some(difference) = difference {
        atom.mass_difference = difference;
    } else {
        prop_lines.push(format!("M  ISO  1 {atom_index:3} {mass:3}"));
    }
}

fn read_properties(
    stream: &mut FortranString<'_>,
    state: &mut MoleculeState,
    mut remaining: usize,
) -> Result<(), FingerprintError> {
    // Avalon❗✔️: struct prop_line_t * ReadProperties(Fortran_FILE *fp,
    // Avalon❗✔️:                                    struct reaccs_molecule_t *mp,
    // Avalon❗✔️:                                    int nprops)
    // Avalon❗✔️: /*
    // Avalon❗✔️:  * Reads the nprops property lines off the file *fp. Interpretable lines
    // Avalon❗✔️:  * are stored in *mp while a pointer to the other lines is returned.
    // Avalon❗✔️:  */
    // Avalon❗✔️: {
    // Avalon❗✔️:    int nentries, n;
    // Avalon❗✔️:    int tmp_ats[8];
    // Avalon❗✔️:    int tmp_vals[8];
    // Avalon❗✔️:    struct prop_line_t *hp, *hhp, *result;
    // Avalon❗✔️:    int atom;
    // Avalon❗✔️:    float value;
    // Avalon❗✔️:
    // Avalon❗✔️:    int old_nprops;
    // Avalon❗✔️:
    // Avalon❗✔️:    old_nprops = nprops;
    // Avalon❗✔️:    result = (struct prop_line_t *)NULL;
    // Avalon❗✔️:    while (nprops-- > 0)
    // Avalon❗✔️: {
    // Avalon❗✔️:       if (STRING_BEGINS(fp->buffer,"M  CHG"))    /* charge property */
    // Avalon❗✔️:       {
    // Avalon❗✔️:          sscanf(fp->buffer+strlen("M  CHG"), "%3d", &nentries);
    // Avalon❗✔️:          n = sscanf(fp->buffer+strlen("M  CHG")+3,
    // Avalon❗✔️:                     " %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
    // Avalon❗✔️:                     &tmp_ats[0],  &tmp_vals[0],  &tmp_ats[1],  &tmp_vals[1],
    // Avalon❗✔️:                     &tmp_ats[2],  &tmp_vals[2],  &tmp_ats[3],  &tmp_vals[3],
    // Avalon❗✔️:                     &tmp_ats[4],  &tmp_vals[4],  &tmp_ats[5],  &tmp_vals[5],
    // Avalon❗✔️:                     &tmp_ats[6],  &tmp_vals[6],  &tmp_ats[7],  &tmp_vals[7]);
    // Avalon❗✔️:          if (n != nentries*2)
    // Avalon❗✔️:          {
    // Avalon❗✔️:             ShowMessageI("n = %d","ReadProperties", n);
    // Avalon❗✔️:             ShowMessageI("nentries = %d","ReadProperties", nentries);
    // Avalon❗✔️:             ShowMessageS("buffer = '%s'","ReadProperties",fp->buffer);
    // Avalon❗✔️:          }
    // Avalon❗✔️:          for (n=0; n<nentries; n++)
    // Avalon❗✔️:             mp->atom_array[tmp_ats[n]-1].charge = tmp_vals[n];
    // Avalon❗✔️:       }
    // Avalon❗✔️:       else if (STRING_BEGINS(fp->buffer,"M  RAD")) /* radical property */
    // Avalon❗✔️:       {
    // Avalon❗✔️:          sscanf(fp->buffer+strlen("M  RAD"), "%3d", &nentries);
    // Avalon❗✔️:          n = sscanf(fp->buffer+strlen("M  RAD")+3,
    // Avalon❗✔️:                     " %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
    // Avalon❗✔️:                     &tmp_ats[0],  &tmp_vals[0],  &tmp_ats[1],  &tmp_vals[1],
    // Avalon❗✔️:                     &tmp_ats[2],  &tmp_vals[2],  &tmp_ats[3],  &tmp_vals[3],
    // Avalon❗✔️:                     &tmp_ats[4],  &tmp_vals[4],  &tmp_ats[5],  &tmp_vals[5],
    // Avalon❗✔️:                     &tmp_ats[6],  &tmp_vals[6],  &tmp_ats[7],  &tmp_vals[7]);
    // Avalon❗✔️:          if (n != nentries*2)
    // Avalon❗✔️:          {
    // Avalon❗✔️:             ShowMessageI("n = %d","ReadProperties", n);
    // Avalon❗✔️:             ShowMessageI("nentries = %d","ReadProperties", nentries);
    // Avalon❗✔️:             ShowMessageS("buffer = '%s'","ReadProperties",fp->buffer);
    // Avalon❗✔️:          }
    // Avalon❗✔️:          for (n=0; n<nentries; n++)
    // Avalon❗✔️:             mp->atom_array[tmp_ats[n]-1].radical = tmp_vals[n];
    // Avalon❗✔️:       }
    // Avalon❗✔️:       else if (STRING_BEGINS(fp->buffer,"M  SUB")) /* subst. degree property */
    // Avalon❗✔️:       {
    // Avalon❗✔️:          sscanf(fp->buffer+strlen("M  SUB"), "%3d", &nentries);
    // Avalon❗✔️:          n = sscanf(fp->buffer+strlen("M  SUB")+3,
    // Avalon❗✔️:                     " %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
    // Avalon❗✔️:                     &tmp_ats[0],  &tmp_vals[0],  &tmp_ats[1],  &tmp_vals[1],
    // Avalon❗✔️:                     &tmp_ats[2],  &tmp_vals[2],  &tmp_ats[3],  &tmp_vals[3],
    // Avalon❗✔️:                     &tmp_ats[4],  &tmp_vals[4],  &tmp_ats[5],  &tmp_vals[5],
    // Avalon❗✔️:                     &tmp_ats[6],  &tmp_vals[6],  &tmp_ats[7],  &tmp_vals[7]);
    // Avalon❗✔️:          if (n != nentries*2)
    // Avalon❗✔️:          {
    // Avalon❗✔️:             ShowMessageI("n = %d","ReadProperties", n);
    // Avalon❗✔️:             ShowMessageI("nentries = %d","ReadProperties", nentries);
    // Avalon❗✔️:             ShowMessageS("buffer = '%s'","ReadProperties",fp->buffer);
    // Avalon❗✔️:          }
    // Avalon❗✔️:          for (n=0; n<nentries; n++)
    // Avalon❗✔️:             mp->atom_array[tmp_ats[n]-1].sub_desc = tmp_vals[n];
    // Avalon❗✔️:       }
    // Avalon❗✔️:       else if (STRING_BEGINS(fp->buffer,"V  ")) /* value for PC files */
    // Avalon❗✔️:       {
    // Avalon❗✔️:          n = sscanf(fp->buffer+strlen("V  "), " %d %f", &atom,  &value);
    // Avalon❗✔️:          if (n == 2)
    // Avalon❗✔️:             mp->atom_array[atom-1].value = value;
    // Avalon❗✔️:          else
    // Avalon❗✔️:          {
    // Avalon❗✔️:             ShowMessage("value line ignored","ReadProperties");
    // Avalon❗✔️:             ShowMessageS("buffer = '%s'","ReadProperties",fp->buffer);
    // Avalon❗✔️:          }
    // Avalon❗✔️:       }
    // Avalon❗✔️:       else if (STRING_BEGINS(fp->buffer,"M  END"))
    // Avalon❗✔️:       {
    // Avalon❗✔️:          if (nprops != 0)
    // Avalon❗✔️:          {
    // Avalon❗✔️:             if (old_nprops != 999)
    // Avalon❗✔️:                ShowMessage("M  END line was not last property line",
    // Avalon❗✔️:                            "ReadProperties");
    // Avalon❗✔️:             nprops = 0;
    // Avalon❗✔️:          }
    // Avalon❗✔️:          // GetBuffer(fp);
    // Avalon❗✔️:          break;
    // Avalon❗✔️:       }
    // Avalon❗✔️:       else if (STRING_BEGINS(fp->buffer,"M  RGP"))      /* R-Group Lines */
    // Avalon❗✔️:       {
    // Avalon❗✔️:          hp = TypeAlloc(1,struct prop_line_t);
    // Avalon❗✔️:          strncpy(hp->text, fp->buffer, MDL_MAXLINE);
    // Avalon❗✔️:          hp->text[MDL_MAXLINE] = '\0';
    // Avalon❗✔️:          hp->next = result; result = hp;
    // Avalon❗✔️:       }
    // Avalon❗✔️:       else if (STRING_BEGINS(fp->buffer,"G  "))       /* Group Lines */
    // Avalon❗✔️:       {
    // Avalon❗✔️:          hp = TypeAlloc(1,struct prop_line_t);
    // Avalon❗✔️:          strncpy(hp->text, fp->buffer, MDL_MAXLINE);
    // Avalon❗✔️:          hp->text[MDL_MAXLINE] = '\0';
    // Avalon❗✔️:          hp->next = result; result = hp;
    // Avalon❗✔️:          GetBuffer(fp); nprops--;
    // Avalon❗✔️:          hp = TypeAlloc(1,struct prop_line_t);
    // Avalon❗✔️:          strncpy(hp->text, fp->buffer, MDL_MAXLINE);
    // Avalon❗✔️:          hp->text[MDL_MAXLINE] = '\0';
    // Avalon❗✔️:          hp->next = result; result = hp;
    // Avalon❗✔️:       }
    // Avalon❗✔️:       else if (STRING_BEGINS(fp->buffer,"A  "))       /* atom text Lines */
    // Avalon❗✔️:       {
    // Avalon❗✔️:          n = sscanf(fp->buffer+strlen("A  "), " %d", &atom);
    // Avalon❗✔️:          if (n==1  &&  atom <= mp->n_atoms  && 0 == strcmp(mp->atom_array[atom-1].atom_symbol, "R"))
    // Avalon❗✔️:          {                                      // interpret atom text as atom label
    // Avalon❗✔️:              GetBuffer(fp); nprops--;
    // Avalon❗✔️:              strncpy(mp->atom_array[atom-1].atext, fp->buffer, MDL_MAXLINE);
    // Avalon❗✔️:              mp->atom_array[atom-1].atext[MDL_MAXLINE] = '\0';
    // Avalon❗✔️:          }
    // Avalon❗✔️:          else                                   // just store as standard property
    // Avalon❗✔️:          {
    // Avalon❗✔️:              hp = TypeAlloc(1,struct prop_line_t);
    // Avalon❗✔️:              strncpy(hp->text, fp->buffer, MDL_MAXLINE);
    // Avalon❗✔️:              hp->text[MDL_MAXLINE] = '\0';
    // Avalon❗✔️:              hp->next = result; result = hp;
    // Avalon❗✔️:              GetBuffer(fp); nprops--;
    // Avalon❗✔️:              hp = TypeAlloc(1,struct prop_line_t);
    // Avalon❗✔️:              strncpy(hp->text, fp->buffer, MDL_MAXLINE);
    // Avalon❗✔️:              hp->text[MDL_MAXLINE] = '\0';
    // Avalon❗✔️:              hp->next = result; result = hp;
    // Avalon❗✔️:           }
    // Avalon❗✔️:       }
    // The `M  ISO` body is reproduced line-for-line in `apply_isotope`
    // immediately above; the source branch framing continues here.
    // Avalon❗✔️:       else if (STRING_BEGINS(fp->buffer,"M  ISO"))      /* isotopic labelling */
    // Avalon❗✔️:       {
    // Avalon❗✔️:          /* Just ignore these, because they are taken care of by      */
    // Avalon❗✔️:          /* the mass difference field. It would require a major       */
    // Avalon❗✔️:          /* change to modify the semantics now. Room for improvement! */
    // Avalon❗✔️:          sscanf(fp->buffer+strlen("M  ISO"), "%3d", &nentries);
    // Avalon❗✔️:          n = sscanf(fp->buffer+strlen("M  ISO")+3,
    // Avalon❗✔️:                     " %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
    // Avalon❗✔️:                     &tmp_ats[0],  &tmp_vals[0],  &tmp_ats[1],  &tmp_vals[1],
    // Avalon❗✔️:                     &tmp_ats[2],  &tmp_vals[2],  &tmp_ats[3],  &tmp_vals[3],
    // Avalon❗✔️:                     &tmp_ats[4],  &tmp_vals[4],  &tmp_ats[5],  &tmp_vals[5],
    // Avalon❗✔️:                     &tmp_ats[6],  &tmp_vals[6],  &tmp_ats[7],  &tmp_vals[7]);
    // Avalon❗✔️:          if (n != nentries*2)
    // Avalon❗✔️:          {
    // Avalon❗✔️:             ShowMessageI("n = %d","ReadProperties", n);
    // Avalon❗✔️:             ShowMessageI("nentries = %d","ReadProperties", nentries);
    // Avalon❗✔️:             ShowMessageS("buffer = '%s'","ReadProperties",fp->buffer);
    // Avalon❗✔️:          }
    // Avalon❗✔️:       }
    // Avalon❗✔️:       else if (STRING_BEGINS(fp->buffer,"M  ALS"))      /* atom type lists */
    // Avalon❗✔️:       {
    // Avalon❗✔️:          /* Just ignore these, because they are taken care of by      */
    // Avalon❗✔️:          /* the atom list lines read before.			      */
    // Avalon❗✔️:       }
    // Avalon❗✔️:       /* ignore short cut lines */
    // Avalon❗✔️:       else if (STRING_BEGINS(fp->buffer,"M  SLB") ||    // Sgroup Labels
    // Avalon❗✔️:                STRING_BEGINS(fp->buffer,"M  STY") ||    // Sgroup Type
    // Avalon❗✔️:                STRING_BEGINS(fp->buffer,"M  SAL") ||    // Sgroup Atom List
    // Avalon❗✔️:                STRING_BEGINS(fp->buffer,"M  SBL") ||    // Sgroup Bond List
    // Avalon❗✔️:                STRING_BEGINS(fp->buffer,"M  SDI") ||    // Sgroup Display Information, e.g. bracket positions
    // Avalon❗✔️:                STRING_BEGINS(fp->buffer,"M  SMT") ||    // Sgroup Subscript text
    // Avalon❗✔️:                STRING_BEGINS(fp->buffer,"M  SBV") ||    // Abbreviation SGroup bond and vector information
    // Avalon❗✔️:                STRING_BEGINS(fp->buffer,"M  SCL") ||    // Abbreviation Sgroup class
    // Avalon❗✔️:                STRING_BEGINS(fp->buffer,"M  SAP") ||    // Abbreviation Sgroup Attachment Point
    // Avalon❗✔️:                STRING_BEGINS(fp->buffer,"M  SDT") ||    // Data Sgroup Field Description
    // Avalon❗✔️:                STRING_BEGINS(fp->buffer,"M  SDD") ||    // Data Sgroup Display Information
    // Avalon❗✔️:                STRING_BEGINS(fp->buffer,"M  SED") ||    // Data Sgroup Data
    // Avalon❗✔️:                STRING_BEGINS(fp->buffer,"M  SDS"))      // Sgroup Expansion
    // Avalon❗✔️:       {
    // Avalon❗✔️:          /* Just ignore those line because they are display-only. */
    // Avalon❗✔️:       }
    // Avalon❗✔️:       else
    // Avalon❗✔️:       {
    // Avalon❗✔️:          ShowMessage("found unknown property","ReadProperties");
    // Avalon❗✔️:          ShowMessageS("buffer = '%s'","ReadProperties",fp->buffer);
    // Avalon❗✔️:          hp = TypeAlloc(1,struct prop_line_t);
    // Avalon❗✔️:          strncpy(hp->text, fp->buffer, MDL_MAXLINE);
    // Avalon❗✔️:          hp->text[MDL_MAXLINE] = '\0';
    // Avalon❗✔️:          hp->next = result; result = hp;
    // Avalon❗✔️:       }
    // Avalon❗✔️:       GetBuffer(fp);
    // Avalon❗✔️:    }
    // Avalon❗✔️:
    // Avalon❗✔️:    /* restore order of list */
    // Avalon❗✔️:    hhp = result; result = (struct prop_line_t *)NULL;
    // Avalon❗✔️:    while (hhp)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       hp = hhp; hhp = hp->next;
    // Avalon❗✔️:       hp->next = result; result = hp;
    // Avalon❗✔️:    }
    // Avalon❗✔️:
    // Avalon❗✔️:    if (fp->buffer[0] != 'M'  &&  fp->status == FORTRAN_NORMAL)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       ShowMessage("possible format error","ReadProperties");
    // Avalon❗✔️:       ShowMessageS("current line: '%s'","ReadProperties",fp->buffer);
    // Avalon❗✔️:    }
    // Avalon❗✔️:    else
    // Avalon❗✔️:        GetBuffer(fp);
    // Avalon❗✔️:    return (result);
    // Avalon❗✔️: }
    while remaining > 0 {
        let line = stream.line()?.to_string();
        remaining -= 1;
        if line.starts_with("M  END") {
            stream.advance();
            break;
        }
        if line.starts_with("M  CHG") || line.starts_with("M  RAD") || line.starts_with("M  SUB") {
            let prefix = &line[..6];
            for (index, value) in property_pairs(&line, prefix)? {
                let atom = property_atom(&mut state.atoms, index)?;
                match prefix {
                    "M  CHG" => atom.charge = value,
                    "M  RAD" => atom.radical = value,
                    "M  SUB" => atom.sub_desc = value,
                    _ => unreachable!(),
                }
            }
        } else if let Some(value_line) = line.strip_prefix("V  ") {
            let mut fields = value_line.split_whitespace();
            if let (Some(index), Some(value)) = (fields.next(), fields.next()) {
                let index = index
                    .parse::<usize>()
                    .map_err(|_| conversion_error("invalid value atom index"))?;
                property_atom(&mut state.atoms, index)?.value = parse_f32(value, "atom value")?;
            }
        } else if line.starts_with("M  RGP") {
            state.prop_lines.push(mdl_property_text(&line));
        } else if line.starts_with("G  ") {
            state.prop_lines.push(mdl_property_text(&line));
            stream.advance();
            remaining = remaining.saturating_sub(1);
            state.prop_lines.push(mdl_property_text(stream.line()?));
        } else if let Some(atom_text) = line.strip_prefix("A  ") {
            let index = atom_text
                .trim()
                .parse::<usize>()
                .map_err(|_| conversion_error("invalid atom-text index"))?;
            stream.advance();
            remaining = remaining.saturating_sub(1);
            let text: String = stream.line()?.chars().take(MDL_MAX_LINE).collect();
            let atom = property_atom(&mut state.atoms, index)?;
            if atom.atom_symbol == "R" {
                atom.atext = text;
            } else {
                state.prop_lines.push(mdl_property_text(&line));
                state.prop_lines.push(text);
            }
        } else if line.starts_with("M  ISO") {
            for (index, mass) in property_pairs(&line, "M  ISO")? {
                let atom = property_atom(&mut state.atoms, index)?;
                apply_isotope(atom, index, mass, &mut state.prop_lines);
            }
        } else if line.starts_with("M  ALS") || is_ignored_sgroup_property(&line) {
        } else {
            state.prop_lines.push(mdl_property_text(&line));
        }
        stream.advance();
    }
    Ok(())
}

fn is_ignored_sgroup_property(line: &str) -> bool {
    [
        "M  SLB", "M  STY", "M  SAL", "M  SBL", "M  SDI", "M  SMT", "M  SBV", "M  SCL", "M  SAP",
        "M  SDT", "M  SDD", "M  SED", "M  SDS",
    ]
    .iter()
    .any(|prefix| line.starts_with(prefix))
}

fn read_v2000_symbol_lists(
    stream: &mut FortranString<'_>,
    count: usize,
) -> Result<Vec<SymbolList>, FingerprintError> {
    // Avalon❗✔️: result = (struct symbol_list_t *)NULL;
    // Avalon❗✔️: for (i=0; i<nlists; i++)
    // Avalon❗✔️: {
    // Avalon❗✔️:    sscanf(fp->buffer,"%d%s%d",&index,type_string,&n_elements);
    // Avalon❗✔️:    for (j=0; j<n_elements; j++)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       if (j==0) strcpy(buffer,""); else strncat(buffer,",",MAX_BUFFER);
    // Avalon❗✔️:       sscanf(&fp->buffer[11+4*j],"%d",&type);
    // Avalon❗✔️:       strncat(buffer,IntToString(periodic_table,type),MAX_BUFFER);
    // Avalon❗✔️:    }
    // Avalon❗✔️:    slp = TypeAlloc(1,struct symbol_list_t);
    // Avalon❗✔️:    slp->next = result; result = slp;
    // Avalon❗✔️:    slp->atom = index; strcpy(slp->string, buffer);
    // Avalon❗✔️:    if (0 == strcmp(type_string,"F")) slp->logic = INCLUSIVE;
    // Avalon❗✔️:    else slp->logic = EXCLUSIVE;
    // Avalon❗✔️:    GetBuffer(fp);
    // Avalon❗✔️: }
    let mut lists = Vec::with_capacity(count);
    for _ in 0..count {
        let line = stream.line()?;
        let mut fields = line.split_whitespace();
        let atom = fields
            .next()
            .ok_or_else(|| conversion_error("atom list is missing atom index"))?
            .parse::<i32>()
            .map_err(|_| conversion_error("invalid atom-list index"))?;
        let inclusive = fields.next() == Some("F");
        let element_count = fields
            .next()
            .ok_or_else(|| conversion_error("atom list is missing element count"))?
            .parse::<usize>()
            .map_err(|_| conversion_error("invalid atom-list element count"))?;
        let mut symbols = Vec::with_capacity(element_count);
        for index in 0..element_count {
            let atomic_number = parse_i32_or_zero(field(line, 11 + 4 * index, 15 + 4 * index));
            let atomic_number = u8::try_from(atomic_number)
                .map_err(|_| conversion_error("invalid atom-list atomic number"))?;
            symbols.push(
                rdkit_element_symbol(atomic_number)
                    .map_err(|error| conversion_error(error.to_string()))?,
            );
        }
        lists.push(SymbolList {
            atom,
            inclusive,
            symbols: symbols.join(","),
        });
        stream.advance();
    }
    // The source prepends each record, so preserve its linked-list traversal
    // order instead of normalizing to file order.
    lists.reverse();
    Ok(lists)
}

fn parse_info_line(line: &str, state: &mut MoleculeState) {
    // Avalon❗✔️: strncpy(buffer,fp->buffer,MAX_BUFFER);
    // Avalon❗✔️: BlankToZero(buffer+22);
    // Avalon❗✔️: strcpy(mp->dimensionality,"2D");
    // Avalon❗✔️: mp->user_initials[0] = '\0'; mp->program_name[0] = '\0';
    // Avalon❗✔️: mp->date[0] = '\0'; mp->time[0] = '\0';
    // Avalon❗✔️: mp->scale1 = 0; mp->scale2 = 0; mp->energy = 0; mp->registry_number = 0;
    // Avalon❗✔️: nitems = sscanf(buffer,
    // Avalon❗✔️:                 "%2c%8c%6c%4c%2c%2d%10f%12f%6ld",
    // Avalon❗✔️:                 mp->user_initials,
    // Avalon❗✔️:                 mp->program_name,
    // Avalon❗✔️:                 mp->date, mp->time,
    // Avalon❗✔️:                 mp->dimensionality,
    // Avalon❗✔️:                &mp->scale1, &mp->scale2,
    // Avalon❗✔️:                &mp->energy,
    // Avalon❗✔️:                &mp->registry_number);
    // Avalon❗✔️: if (regno != 0) mp->registry_number = regno;
    // Avalon❗✔️:
    // Avalon❗✔️: mp->user_initials[2]  = '\0'; RemoveTrailingBlanks(mp->user_initials);
    // Avalon❗✔️: mp->program_name[8]   = '\0'; RemoveTrailingBlanks(mp->program_name);
    // Avalon❗✔️: mp->date[6]           = '\0'; RemoveTrailingBlanks(mp->date);
    // Avalon❗✔️: mp->time[4]           = '\0'; RemoveTrailingBlanks(mp->time);
    // Avalon❗✔️: mp->dimensionality[2] = '\0'; RemoveTrailingBlanks(mp->dimensionality);
    state.user_initials.clear();
    state.program_name.clear();
    state.date.clear();
    state.time.clear();
    state.dimensionality = "2D".to_string();
    state.scale1 = 0;
    state.scale2 = 0.0;
    state.energy = 0.0;
    state.registry_number = 0;
    if !line.is_empty() {
        state.user_initials = field(line, 0, 2).trim_end().to_string();
    }
    if line.len() > 2 {
        state.program_name = field(line, 2, 10).trim_end().to_string();
    }
    if line.len() > 10 {
        state.date = field(line, 10, 16).trim_end().to_string();
    }
    if line.len() > 16 {
        state.time = field(line, 16, 20).trim_end().to_string();
    }
    if line.len() <= 20 {
        return;
    }
    state.dimensionality = field(line, 20, 22).trim_end().to_string();
    if line.len() < 22 {
        return;
    }

    let mut transformed = line.as_bytes().to_vec();
    for byte in transformed.iter_mut().skip(22) {
        if *byte == b' ' {
            *byte = b'0';
        }
    }
    let Ok(transformed) = std::str::from_utf8(&transformed) else {
        return;
    };
    let mut cursor = 22;
    let Some(scale1) = scan_decimal_field(transformed, &mut cursor, 2, false) else {
        return;
    };
    state.scale1 = scale1.parse().unwrap_or(0);
    let Some(scale2) = scan_decimal_field(transformed, &mut cursor, 10, true) else {
        return;
    };
    state.scale2 = scale2.parse().unwrap_or(0.0);
    let Some(energy) = scan_decimal_field(transformed, &mut cursor, 12, true) else {
        return;
    };
    state.energy = energy.parse().unwrap_or(0.0);
    let Some(registry_number) = scan_decimal_field(transformed, &mut cursor, 6, false) else {
        return;
    };
    state.registry_number = registry_number.parse().unwrap_or(0);
}

fn parse_v2000_counts(
    line: &str,
) -> Result<(usize, usize, usize, i32, i32, usize, usize), FingerprintError> {
    // Avalon❗✔️: BlankToZero(buffer);
    // Avalon❗✔️: nitems = sscanf(buffer,"%3d%3d%3d%3d%3d%3d%*12c%3d",
    // Avalon❗✔️:                 &mp->n_atoms, &mp->n_bonds, &mp->n_atom_lists,
    // Avalon❗✔️:                 &mp->dummy1, &mp->chiral_flag, &n_stexts, &mp->n_props);
    // Avalon❗✔️: if (nitems >= 2  &&
    // Avalon❗✔️:     mp->n_atoms >=0  &&  mp->n_atoms <= MAXATOMS &&
    // Avalon❗✔️:     mp->n_bonds >=0  &&  mp->n_bonds <= MAXBONDS)
    // Avalon❗✔️: {
    // Avalon❗✔️:    GetBuffer(fp);
    // Avalon❗✔️:    if (fp->status != FORTRAN_NORMAL) return(fp->status);
    // Avalon❗✔️: }
    // Avalon❗✔️: else
    // Avalon❗✔️: {
    // Avalon❗✔️:    ShowMessage("incorrect syntax of arguments on atom/bond number line\n", "ReadREACCSMolecule");
    // Avalon❗✔️:    fprintf(stderr,"buffer = '%s'\n", fp->buffer);
    // Avalon❗✔️:    return(FORTRAN_ERROR);
    // Avalon❗✔️: }
    let n_atoms = usize::try_from(parse_i32(field(line, 0, 3), "V2000 atom count")?)
        .map_err(|_| conversion_error("negative V2000 atom count"))?;
    let n_bonds = usize::try_from(parse_i32(field(line, 3, 6), "V2000 bond count")?)
        .map_err(|_| conversion_error("negative V2000 bond count"))?;
    if n_atoms > MAX_ATOMS || n_bonds > MAX_BONDS {
        return Err(conversion_error(
            "V2000 atom or bond count exceeds Avalon limits",
        ));
    }
    Ok((
        n_atoms,
        n_bonds,
        usize::try_from(parse_i32_or_zero(field(line, 6, 9))).unwrap_or(0),
        parse_i32_or_zero(field(line, 9, 12)),
        parse_i32_or_zero(field(line, 12, 15)),
        usize::try_from(parse_i32_or_zero(field(line, 15, 18))).unwrap_or(0),
        usize::try_from(parse_i32_or_zero(field(line, 30, 33))).unwrap_or(0),
    ))
}

fn read_v3000(
    stream: &mut FortranString<'_>,
    state: &mut MoleculeState,
) -> Result<(), FingerprintError> {
    // Avalon❗✔️: if (0 == strcmp(mp->version, "V3000"))
    // Avalon❗✔️: {
    // Avalon❗✔️:     GetBuffer(fp);
    // Avalon❗✔️:     mp->n_atom_lists = 0;
    // Avalon❗✔️:     mp->dummy1       = 0;
    // Avalon❗✔️:     mp->chiral_flag  = 0;
    // Avalon❗✔️:     n_stexts         = 0;
    // Avalon❗✔️:     mp->n_props      = 0;
    // Avalon❗✔️:     mp->n_atoms      = 0;
    // Avalon❗✔️:     mp->n_atoms      = 0;
    // Avalon❗✔️:     while (fp->status == FORTRAN_NORMAL  &&
    // Avalon❗✔️:            0 != strcmp("M  END", fp->buffer))
    // Avalon❗✔️:    {
    // Avalon❗✔️:       if (0 == strncmp("M  V30 BEGIN CTAB", fp->buffer,
    // Avalon❗✔️:                              strlen("M  V30 BEGIN CTAB")))
    // Avalon❗✔️:          GetBuffer(fp);
    // Avalon❗✔️:       else if (0 == strncmp("M  V30 COUNTS ", fp->buffer,
    // Avalon❗✔️:                             strlen("M  V30 COUNTS ")))
    // Avalon❗✔️:       {
    // Avalon❗✔️:          regtext[0] = '\0';
    // Avalon❗✔️:          nitems = sscanf(fp->buffer + strlen("M  V30 COUNTS "),
    // Avalon❗✔️:                          "%d %d %d %d %d %s",
    // Avalon❗✔️:                                     &mp->n_atoms, &mp->n_bonds,
    // Avalon❗✔️:                                     &n_sgroups,
    // Avalon❗✔️:                                     &n_3d,
    // Avalon❗✔️:                                     &mp->chiral_flag,
    // Avalon❗✔️:                                     &regtext);
    // Avalon❗✔️:          mp->atom_array = TypeAlloc(mp->n_atoms,struct reaccs_atom_t);
    // Avalon❗✔️:          mp->bond_array = TypeAlloc(mp->n_bonds,struct reaccs_bond_t);
    // Avalon❗✔️:          GetBuffer(fp);
    // Avalon❗✔️:       }
    // Avalon❗✔️:       else if (0 == strncmp("M  V30 BEGIN ATOM", fp->buffer,
    // Avalon❗✔️:                             strlen("M  V30 BEGIN ATOM")))
    // Avalon❗✔️:       {
    // Avalon❗✔️:          GetBuffer(fp);
    // Avalon❗✔️:          i=0;
    // Avalon❗✔️:          while (fp->status == FORTRAN_NORMAL  &&
    // Avalon❗✔️:                 0 != strcmp("M  V30 END ATOM", fp->buffer))
    // Avalon❗✔️:          {
    // Avalon❗✔️:             status = ReadV30Atom(fp,&mp->atom_array[i++], mp);
    // Avalon❗✔️:             if (status != FORTRAN_NORMAL) return(status);
    // Avalon❗✔️:          }
    // Avalon❗✔️:          if (fp->status == FORTRAN_NORMAL)
    // Avalon❗✔️:              GetBuffer(fp); // skip 'M  V30 END ATOM'
    // Avalon❗✔️:       }
    // Avalon❗✔️:       else if (0 == strncmp("M  V30 BEGIN BOND", fp->buffer,
    // Avalon❗✔️:                             strlen("M  V30 BEGIN BOND")))
    // Avalon❗✔️:       {
    // Avalon❗✔️:          GetBuffer(fp);
    // Avalon❗✔️:          i=0;
    // Avalon❗✔️:          while (fp->status == FORTRAN_NORMAL  &&
    // Avalon❗✔️:                 0 != strcmp("M  V30 END BOND", fp->buffer))
    // Avalon❗✔️:          {
    // Avalon❗✔️:             status = ReadV30Bond(fp,&mp->bond_array[i++]);
    // Avalon❗✔️:             if (status != FORTRAN_NORMAL) return(status);
    // Avalon❗✔️:          }
    // Avalon❗✔️:          if (fp->status == FORTRAN_NORMAL)
    // Avalon❗✔️:              GetBuffer(fp); // skip 'M  V30 END BOND'
    // Avalon❗✔️:       }
    // Avalon❗✔️:       else if (0 == strncmp("M  V30 BEGIN COLLECTION", fp->buffer,
    // Avalon❗✔️:                             strlen("M  V30 BEGIN COLLECTION")))
    // Avalon❗✔️:       {
    // Avalon❗✔️:          while (fp->status == FORTRAN_NORMAL  &&
    // Avalon❗✔️:                 0 != strcmp("M  V30 END COLLECTION", fp->buffer))
    // Avalon❗✔️:          {
    // Avalon❗✔️:              GetBuffer(fp);
    // Avalon❗✔️:          }
    // Avalon❗✔️:          if (fp->status == FORTRAN_NORMAL)
    // Avalon❗✔️:              GetBuffer(fp); // skip 'M  V30 END COLLECTION'
    // Avalon❗✔️:       }
    // Avalon❗✔️:       else if (0 == strncmp("M  V30 END CTAB", fp->buffer,
    // Avalon❗✔️:                              strlen("M  V30 END CTAB")))
    // Avalon❗✔️:          GetBuffer(fp);
    // Avalon❗✔️:       else
    // Avalon❗✔️:       {
    // Avalon❗✔️:          GetBuffer(fp);
    // Avalon❗✔️:       }
    // Avalon❗✔️:    }
    // Avalon❗✔️:    if (fp->status != FORTRAN_NORMAL) return(fp->status);
    // Avalon❗✔️:    // parsing was in V3000 but still using old representation
    // Avalon❗✔️:    strcpy(mp->version, "V2000");
    // Avalon❗✔️: }
    stream.advance();
    let mut expected_atoms = None;
    let mut expected_bonds = None;
    while let Some(line) = stream.current {
        if line == "M  END" {
            stream.advance();
            break;
        }
        if line.starts_with("M  V30 COUNTS ") {
            let fields: Vec<&str> = line[13..].split_whitespace().collect();
            if fields.len() < 5 {
                return Err(conversion_error("V3000 COUNTS has fewer than five fields"));
            }
            expected_atoms = Some(
                fields[0]
                    .parse::<usize>()
                    .map_err(|_| conversion_error("invalid V3000 atom count"))?,
            );
            expected_bonds = Some(
                fields[1]
                    .parse::<usize>()
                    .map_err(|_| conversion_error("invalid V3000 bond count"))?,
            );
            state.chiral_flag = parse_i32(fields[4], "V3000 chiral flag")?;
            stream.advance();
        } else if line == "M  V30 BEGIN ATOM" {
            stream.advance();
            while stream.line()? != "M  V30 END ATOM" {
                let atom_line = stream.line()?.to_string();
                let index = state.atoms.len();
                let atom = parse_v30_atom(&atom_line, index, state)?;
                state.atoms.push(atom);
                stream.advance();
            }
            // `ParseV30SymbolList` prepends each list to the linked list.
            state.symbol_lists.reverse();
            stream.advance();
        } else if line == "M  V30 BEGIN BOND" {
            stream.advance();
            while stream.line()? != "M  V30 END BOND" {
                state.bonds.push(parse_v30_bond(stream.line()?)?);
                stream.advance();
            }
            stream.advance();
        } else if line == "M  V30 BEGIN COLLECTION" {
            while stream.line()? != "M  V30 END COLLECTION" {
                stream.advance();
            }
            stream.advance();
        } else {
            stream.advance();
        }
    }
    if expected_atoms != Some(state.atoms.len()) || expected_bonds != Some(state.bonds.len()) {
        return Err(conversion_error(
            "V3000 parsed atom or bond count differs from COUNTS",
        ));
    }
    state.version = "V2000".to_string();
    Ok(())
}

fn read_v2000(
    stream: &mut FortranString<'_>,
    counts: &str,
    state: &mut MoleculeState,
) -> Result<(), FingerprintError> {
    // Avalon❗✔️: else                                     // old MOL format
    // Avalon❗✔️: {
    // Avalon❗✔️:    BlankToZero(buffer);
    // Avalon❗✔️:     mp->n_atom_lists = 0;
    // Avalon❗✔️:     mp->dummy1       = 0;
    // Avalon❗✔️:     mp->chiral_flag  = 0;
    // Avalon❗✔️:     n_stexts         = 0;
    // Avalon❗✔️:     mp->n_props      = 0;
    // Avalon❗✔️:     mp->n_atoms      = 0;
    // Avalon❗✔️:     mp->n_atoms      = 0;
    // The fixed-column `sscanf` and count validation block is reproduced
    // line-for-line in `parse_v2000_counts` immediately above.
    // Avalon❗✔️:     if (mp->version[0] == '\0') mp->version[0] = ' ';
    // Avalon❗✔️:
    // Avalon❗✔️:     mp->atom_array = TypeAlloc(mp->n_atoms,struct reaccs_atom_t);
    // Avalon❗✔️:     for (i=0; i<mp->n_atoms; i++)
    // Avalon❗✔️:     {
    // Avalon❗✔️:        status = ReadREACCSAtom(fp,&mp->atom_array[i]);
    // Avalon❗✔️:        if (status != FORTRAN_NORMAL) return(status);
    // Avalon❗✔️:     }
    // Avalon❗✔️:
    // Avalon❗✔️:     mp->bond_array = TypeAlloc(mp->n_bonds,struct reaccs_bond_t);
    // Avalon❗✔️:     for (i=0; i<mp->n_bonds; i++)
    // Avalon❗✔️:     {
    // Avalon❗✔️:        status = ReadREACCSBond(fp,&mp->bond_array[i]);
    // Avalon❗✔️:        if (status != FORTRAN_NORMAL) return(status);
    // Avalon❗✔️:     }
    // Avalon❗✔️:
    // Avalon❗✔️:     mp->symbol_lists = ReadSymbolLists(fp,mp->n_atom_lists);
    // Avalon❗✔️:
    // Avalon❗✔️:     mp->stext_lines = (struct stext_line_t *)NULL;
    // Avalon❗✔️:     for (i=0; i<n_stexts; i++)
    // Avalon❗✔️:     {
    // Avalon❗✔️:        stp = TypeAlloc(1,struct stext_line_t);
    // Avalon❗✔️:        sscanf(fp->buffer,"%f %f",&stp->x, &stp->y);
    // Avalon❗✔️:        GetBuffer(fp);
    // Avalon❗✔️:        strncpy(stp->text,fp->buffer,MDL_MAXLINE);
    // Avalon❗✔️:        GetBuffer(fp);
    // Avalon❗✔️:        stp->next = mp->stext_lines;
    // Avalon❗✔️:        mp->stext_lines = stp;
    // Avalon❗✔️:     }
    // Avalon❗✔️:     stp = mp->stext_lines;
    // Avalon❗✔️:     mp->stext_lines = (struct stext_line_t *)NULL;
    // Avalon❗✔️:     while (stp)
    // Avalon❗✔️:     {
    // Avalon❗✔️:        stph = stp->next;
    // Avalon❗✔️:        stp->next = mp->stext_lines;
    // Avalon❗✔️:        mp->stext_lines = stp;
    // Avalon❗✔️:        stp = stph;
    // Avalon❗✔️:     }
    // Avalon❗✔️:
    // Avalon❗✔️:     mp->prop_lines = ReadProperties(fp,mp,mp->n_props);
    // Avalon❗✔️:     mp->n_props = CountPropLines(mp->prop_lines);
    // Avalon❗✔️: }
    let (n_atoms, n_bonds, n_atom_lists, dummy1, chiral_flag, n_stexts, n_props) =
        parse_v2000_counts(counts)?;
    state.dummy1 = dummy1;
    state.chiral_flag = chiral_flag;
    if state.version.is_empty() {
        state.version = " ".to_string();
    }
    stream.advance();
    state.atoms.reserve(n_atoms);
    for _ in 0..n_atoms {
        state.atoms.push(read_v2000_atom(stream)?);
    }
    state.bonds.reserve(n_bonds);
    for _ in 0..n_bonds {
        state.bonds.push(read_v2000_bond(stream)?);
    }
    state.symbol_lists = read_v2000_symbol_lists(stream, n_atom_lists)?;
    state.stext_lines.reserve(n_stexts);
    for _ in 0..n_stexts {
        let coordinates = stream.line()?;
        let mut fields = coordinates.split_whitespace();
        let x = fields
            .next()
            .map_or(Ok(0.0), |value| parse_f32(value, "S-text x"))?;
        let y = fields
            .next()
            .map_or(Ok(0.0), |value| parse_f32(value, "S-text y"))?;
        stream.advance();
        let text: String = stream.line()?.chars().take(MDL_MAX_LINE).collect();
        stream.advance();
        state.stext_lines.push(STextLine { x, y, text });
    }
    read_properties(stream, state, n_props)
}

pub(super) fn parse_mol_block(block: &str) -> Result<MoleculeState, FingerprintError> {
    // Avalon❗✔️: int ReadREACCSMolecule(Fortran_FILE *fp,
    // Avalon❗✔️:                        struct reaccs_molecule_t *mp,
    // Avalon❗✔️:                        const char *label)
    // Avalon❗✔️: {
    // Avalon❗✔️:    if (IsNULL(mp)) return(FORTRAN_ERROR);
    // Avalon❗✔️:    // set some safe defaults
    // Avalon❗✔️:    mp->name[0] = '\0';
    // Avalon❗✔️:    mp->scale1 = 0; mp->scale2 = 0; mp->energy = 0; mp->registry_number = 0;
    // Avalon❗✔️:    mp->user_initials[0] = '\0'; strcpy(mp->program_name,"DUMMY");
    // Avalon❗✔️:    mp->date[0] = '\0'; mp->time[0] = '\0'; mp->dimensionality[0]= '\0';
    // Avalon❗✔️:    mp->comment[0] = '\0'; mp->n_atom_lists = 0; mp->dummy1 = 0;
    // Avalon❗✔️:    mp->chiral_flag = 0; n_stexts = 0; mp->stext_lines = NULL;
    // Avalon❗✔️:    mp->n_props = 0; mp->version[0] = '\0';
    // Avalon❗✔️:    mp->atom_array = NULL; mp->bond_array = NULL;
    // Avalon❗✔️:    if (IsNULL(fp) || IsNULL(label)) return(FORTRAN_ERROR);
    // Avalon❗✔️:    if (fp->status != FORTRAN_NORMAL) return(fp->status);
    // Avalon❗✔️:    if (0 != strcmp(label,""))
    // Avalon❗✔️:    {
    // Avalon❗✔️:       if (!SearchString(fp,label,"$MFMT")) return (fp->status);
    // Avalon❗✔️:       if (STRING_BEGINS(fp->buffer,"$MFMT $MIREG"))
    // Avalon❗✔️:          sscanf(fp->buffer+strlen("$MFMT $MIREG"),"%ld",&regno);
    // Avalon❗✔️:       else
    // Avalon❗✔️:          regno = 0;
    // Avalon❗✔️:       GetBuffer(fp);                            /* skip header line */
    // Avalon❗✔️:       if (fp->status != FORTRAN_NORMAL)
    // Avalon❗✔️:        {
    // Avalon❗✔️: fprintf(stderr, "ReadREACCSMolecule: fp->status(2) = %d\n", fp->status);
    // Avalon❗✔️:           return(fp->status);
    // Avalon❗✔️:        }
    // Avalon❗✔️:    }
    // Avalon❗✔️:    if (fp->buffer[0] == '>' && strchr(fp->buffer,'<'))
    // Avalon❗✔️:    {
    // Avalon❗✔️:       mp->atom_array = TypeAlloc(0, struct reaccs_atom_t);
    // Avalon❗✔️:       mp->bond_array = TypeAlloc(0, struct reaccs_bond_t);
    // Avalon❗✔️:       return(fp->status);
    // Avalon❗✔️:    }
    // Avalon❗✔️:
    // Avalon❗✔️:    strncpy(mp->name,fp->buffer,MAXNAME);
    // Avalon❗✔️:    mp->name[MAXNAME] = '\0';
    // Avalon❗✔️:    GetBuffer(fp);
    // Avalon❗✔️:    if (fp->status != FORTRAN_NORMAL) return(fp->status);
    // The complete info-line source block is reproduced in `parse_info_line`.
    // Avalon❗✔️:    GetBuffer(fp);
    // Avalon❗✔️:    if (fp->status != FORTRAN_NORMAL) return(fp->status);
    // Avalon❗✔️:
    // Avalon❗✔️:    strncpy(mp->comment,fp->buffer,MDL_MAXLINE);
    // Avalon❗✔️:    GetBuffer(fp);
    // Avalon❗✔️:    if (fp->status != FORTRAN_NORMAL) return(fp->status);
    // Avalon❗✔️:    RemoveTrailingBlanks(mp->comment);
    // Avalon❗✔️:
    // Avalon❗✔️:    strncpy(buffer,fp->buffer,MAX_BUFFER);
    // Avalon❗✔️:    mp->version[0]   = '\0';
    // Avalon❗✔️:    if (strlen(buffer) > 34)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       strncpy(mp->version, buffer+34, 6);
    // Avalon❗✔️:       mp->version[5] = '\0';
    // Avalon❗✔️:       // fix a common glitch
    // Avalon❗✔️:       if (NULL != strstr(mp->version, "V200")) strcpy(mp->version, "V2000");
    // Avalon❗✔️:    }
    // The complete V3000 and V2000 branches are reproduced in `read_v3000`
    // and `read_v2000`, respectively.
    // Avalon❗✔️:    return(FORTRAN_NORMAL);
    // Avalon❗✔️: }
    let mut stream = FortranString::open(block);
    let mut state = MoleculeState::default();
    let first = stream.line()?;
    if first.starts_with('>') && first.contains('<') {
        return Ok(state);
    }
    state.name = first.chars().take(80).collect();
    stream.advance();
    parse_info_line(stream.line()?, &mut state);
    stream.advance();
    state.comment = stream.line()?.chars().take(MDL_MAX_LINE).collect();
    stream.advance();
    let counts = stream.line()?.to_string();
    if counts.len() > 34 {
        state.version = field(&counts, 34, 40).chars().take(5).collect();
        if state.version.contains("V200") {
            state.version = "V2000".to_string();
        }
    }
    if state.version == "V3000" {
        read_v3000(&mut stream, &mut state)?;
    } else {
        read_v2000(&mut stream, &counts, &mut state)?;
    }
    Ok(state)
}

pub(super) fn mol_to_reaccs(molecule: &Molecule) -> Result<MoleculeState, FingerprintError> {
    // Avalon❗✔️: struct reaccs_molecule_t *molToReaccs(const ROMol &mol) {
    // Avalon❗✔️:   std::string molB = MolToMolBlock(mol, true);
    // Avalon❗✔️:   Utils::LocaleSwitcher ls;
    // Avalon❗✔️:   struct reaccs_molecule_t *res = MolStr2Mol((char *)molB.c_str());
    // Avalon❗✔️:   POSTCONDITION(res, "could not build a molecule");
    // Avalon❗✔️:   return res;
    // Avalon❗✔️: }
    // Avalon❗✔️: struct reaccs_molecule_t * MolStr2Mol(char * MolStr)
    // Avalon❗✔️: {
    // Avalon❗✔️:    fp = FortranStringOpen(MolStr);
    // Avalon❗✔️:    mp = TypeAlloc(1, struct reaccs_molecule_t);
    // Avalon❗✔️:    if (ReadREACCSMolecule(fp, mp, "") != FORTRAN_NORMAL)
    // Avalon❗✔️:    { FreeMolecule(mp); mp = NULL; }
    // Avalon❗✔️:    FortranClose(fp);
    // Avalon❗✔️:    return mp;
    // Avalon❗✔️: }
    let params = MolBlockWriteParams {
        format: SdfFormat::V2000,
        force_2d: false,
        include_stereo: true,
        kekulize: true,
        precision: 6,
    };
    let block = mol_to_mol_block_with_params(molecule, &params)
        .map_err(|error| conversion_error(error.to_string()))?;
    parse_mol_block(&block)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_source_v2000_atom_defaults_and_properties() {
        let block = concat!(
            "fixture\n",
            "  RDKit          2D\n",
            "comment\n",
            "  1  0  0  0  0  0  0  0  0  0999 V2000\n",
            "    0.0000    0.0000    0.0000 C   0  0\n",
            "M  CHG  1   1  -1\n",
            "M  ISO  1   1  13\n",
            "M  END\n",
        );
        let state = parse_mol_block(block).unwrap();
        assert_eq!(state.name, "fixture");
        assert_eq!(state.atoms.len(), 1);
        assert_eq!(state.atoms[0].charge, -1);
        assert_eq!(state.atoms[0].mass_difference, 1);
        assert!(state.prop_lines.is_empty());
    }

    #[test]
    fn unsupported_isotope_property_is_preserved_in_source_form() {
        let block = concat!(
            "\n\n\n",
            "  1  0  0  0  0  0  0  0  0  0999 V2000\n",
            "    0.0000    0.0000    0.0000 Xe  0  0\n",
            "M  ISO  1   1 129\n",
            "M  END\n",
        );
        let state = parse_mol_block(block).unwrap();
        assert_eq!(state.atoms[0].mass_difference, 0);
        assert_eq!(state.prop_lines, ["M  ISO  1   1 129"]);
    }

    #[test]
    fn retained_property_lines_use_the_source_mdl_limit() {
        let long_property = format!("M  XYZ{}", "x".repeat(MDL_MAX_LINE));
        let block =
            format!("\n\n\n  0  0  0  0  0  0  0  0  0  0999 V2000\n{long_property}\nM  END\n");
        let state = parse_mol_block(&block).unwrap();
        assert_eq!(state.prop_lines.len(), 1);
        assert_eq!(state.prop_lines[0].len(), MDL_MAX_LINE);
        assert_eq!(state.prop_lines[0], long_property[..MDL_MAX_LINE]);
    }

    #[test]
    fn rdkit_adapter_path_converts_an_ordinary_molecule() {
        let molecule = Molecule::from_smiles("[13CH3][NH2+]C(=O)O").unwrap();
        let state = mol_to_reaccs(&molecule).unwrap();
        assert_eq!(state.atoms.len(), molecule.num_atoms());
        assert_eq!(state.bonds.len(), molecule.num_bonds());
        assert_eq!(state.atoms[0].mass_difference, 1);
        assert_eq!(state.atoms[1].charge, 1);
    }

    #[test]
    fn v3000_path_uses_source_atom_and_bond_fields() {
        let block = concat!(
            "fixture\n",
            "  RDKit          2D\n",
            "\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 2 1 0 0 1\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 C 0.0 0.0 0.0 7 MASS=13 CHG=-1 CFG=2\n",
            "M  V30 2 O 1.0 0.0 0.0 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 1 2 1 2 CFG=2\n",
            "M  V30 END BOND\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );
        let state = parse_mol_block(block).unwrap();
        assert_eq!(state.version, "V2000");
        assert_eq!(state.program_name, "RDKit");
        assert_eq!(state.time, " 2D");
        assert_eq!(state.dimensionality, "2D");
        assert_eq!(state.chiral_flag, 1);
        assert_eq!(state.atoms[0].mapping, 7);
        assert_eq!(state.atoms[0].mass_difference, 1);
        assert_eq!(state.atoms[0].charge, -1);
        assert_eq!(state.atoms[0].stereo_parity, 2);
        assert_eq!(state.bonds[0].bond_type, 2);
        assert_eq!(state.bonds[0].stereo_symbol, 0);
    }

    #[test]
    fn v3000_mass_uses_source_float_table_before_double_subtraction() {
        let block = concat!(
            "\n\n\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 1 0 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 N 0 0 0 0 MASS=14.5166994\n",
            "M  V30 END ATOM\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );
        let state = parse_mol_block(block).unwrap();
        assert_eq!(state.atoms[0].mass_difference, 1);
    }

    #[test]
    fn v3000_over_capacity_symbol_list_fails_before_source_overflow() {
        let symbol = format!("[{}]", "C".repeat(MAX_SYMBOL_LIST + 1));
        let block = format!(
            concat!(
                "\n\n\n",
                "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
                "M  V30 BEGIN CTAB\n",
                "M  V30 COUNTS 1 0 0 0 0\n",
                "M  V30 BEGIN ATOM\n",
                "M  V30 1 {} 0 0 0 0\n",
                "M  V30 END ATOM\n",
                "M  V30 END CTAB\n",
                "M  END\n",
            ),
            symbol
        );
        let error = parse_mol_block(&block).unwrap_err();
        assert!(matches!(
            error,
            FingerprintError::AvalonConversion { ref reason }
                if reason == "V3000 atom list exceeds Avalon's symbol-list capacity"
        ));
    }

    #[test]
    fn native_v2000_full_state_profile_matches_every_field() {
        // Golden state captured with `tmp/parity-audit/avalon_reaccs_state.c`
        // linked to the pinned wheel's `libRDKitavalon_clib`.
        let block = concat!(
            "full-state\n",
            "ABRDKit 010203042D 1  1.250000    2.500000    42\n",
            "native fixture\n",
            "  5  4  0  0  1  1  0  0  0  0999 V2000\n",
            "    0.0000    0.0000    0.0000 C   0  0\n",
            "    1.0000    0.0000    0.0000 N   0  0\n",
            "    2.0000    0.0000    0.0000 O   0  0\n",
            "    3.0000    0.0000    0.0000 H   0  0\n",
            "    4.0000    0.0000    0.0000 C   1  0\n",
            "  1  2  1  0\n",
            "  2  3  2  0\n",
            "  3  4  4  0\n",
            "  4  5  1  0\n",
            "1.5 2.5\n",
            "source text\n",
            "M  CHG  1   2   1\n",
            "M  RAD  1   3   2\n",
            "M  SUB  1   1   3\n",
            "M  XYZ native-property\n",
            "M  END\n",
        );
        let mut atoms = vec![Atom::default(); 5];
        for (index, atom) in atoms.iter_mut().enumerate() {
            atom.x = index as f32;
        }
        atoms[0].atom_symbol = "C".to_string();
        atoms[0].sub_desc = 3;
        atoms[1].atom_symbol = "N".to_string();
        atoms[1].charge = 1;
        atoms[2].atom_symbol = "O".to_string();
        atoms[2].radical = 2;
        atoms[3].atom_symbol = "H".to_string();
        atoms[4].atom_symbol = "C".to_string();
        atoms[4].mass_difference = 1;
        let bonds = [(1, 2, 1), (2, 3, 2), (3, 4, 4), (4, 5, 1)]
            .map(|(begin, end, bond_type)| Bond {
                atoms: [begin, end],
                bond_type,
                ..Bond::default()
            })
            .to_vec();
        let expected = MoleculeState {
            name: "full-state".to_string(),
            user_initials: "AB".to_string(),
            program_name: "RDKit 01".to_string(),
            date: "020304".to_string(),
            time: "2D 1".to_string(),
            dimensionality: String::new(),
            scale1: 1,
            scale2: 0.25,
            energy: 2.5,
            registry_number: 42,
            comment: "native fixture".to_string(),
            chiral_flag: 1,
            version: "V2000".to_string(),
            atoms,
            bonds,
            stext_lines: vec![STextLine {
                x: 1.5,
                y: 2.5,
                text: "source text".to_string(),
            }],
            prop_lines: vec!["M  XYZ native-property".to_string()],
            ..MoleculeState::default()
        };
        assert_eq!(parse_mol_block(block).unwrap(), expected);
    }

    #[test]
    fn native_v3000_query_and_explicit_h_full_state_match_every_field() {
        let block = concat!(
            "query\n",
            "  RDKit          2D\n",
            "\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 2 1 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 [C,N] 0 0 0 0 HCOUNT=2\n",
            "M  V30 2 H 1 0 0 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 1 1 1 2\n",
            "M  V30 END BOND\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        );
        let mut query_atom = Atom {
            atom_symbol: "L".to_string(),
            query_h_count: 2,
            ..Atom::default()
        };
        query_atom.x = 0.0;
        let expected = MoleculeState {
            name: "query".to_string(),
            user_initials: String::new(),
            program_name: "RDKit".to_string(),
            date: String::new(),
            time: " 2D".to_string(),
            dimensionality: "2D".to_string(),
            version: "V2000".to_string(),
            atoms: vec![
                query_atom,
                Atom {
                    x: 1.0,
                    atom_symbol: "H".to_string(),
                    ..Atom::default()
                },
            ],
            bonds: vec![Bond {
                atoms: [1, 2],
                bond_type: 1,
                ..Bond::default()
            }],
            symbol_lists: vec![SymbolList {
                atom: 1,
                inclusive: true,
                symbols: "C,N".to_string(),
            }],
            ..MoleculeState::default()
        };
        assert_eq!(parse_mol_block(block).unwrap(), expected);
    }

    #[test]
    fn native_v2000_maximum_atom_table_boundary_matches_full_state() {
        let mut block = String::from("\n\n\n999  0  0  0  0  0  0  0  0  0999 V2000\n");
        let atom_line = "    0.0000    0.0000    0.0000 C   0  0\n";
        for _ in 0..MAX_ATOMS {
            block.push_str(atom_line);
        }
        block.push_str("M  END\n");
        let expected = MoleculeState {
            dimensionality: "2D".to_string(),
            version: "V2000".to_string(),
            atoms: vec![
                Atom {
                    atom_symbol: "C".to_string(),
                    ..Atom::default()
                };
                MAX_ATOMS
            ],
            ..MoleculeState::default()
        };
        let mut expected = expected;
        expected.program_name.clear();
        assert_eq!(parse_mol_block(&block).unwrap(), expected);
    }
}
