use std::collections::HashMap;

use crate::{Atom, Bond, BondOrder, Molecule, SmilesParseError};

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum PendingBond {
    Unspecified,
    Single,
    Double,
    Triple,
    Quadruple,
    Aromatic,
    DirectionalSingle,
    DativeForward,
    DativeBackward,
    Null,
}

#[derive(Debug, Clone)]
struct RingOpen {
    atom_index: usize,
    bond: PendingBond,
}

pub(crate) fn parse_smiles(smiles: &str) -> Result<Molecule, SmilesParseError> {
    let mut parser = Parser::new(smiles);
    let mut mol = parser.parse_molecule()?;
    parser.skip_ascii_whitespace();
    if !parser.is_eof() {
        return Err(parser.error("unexpected trailing characters"));
    }
    mol.rebuild_adjacency();
    Ok(mol)
}

struct Parser<'a> {
    input: &'a str,
    pos: usize,
    ring_opens: HashMap<u32, RingOpen>,
}

impl<'a> Parser<'a> {
    fn new(input: &'a str) -> Self {
        Self {
            input,
            pos: 0,
            ring_opens: HashMap::new(),
        }
    }

    fn parse_molecule(&mut self) -> Result<Molecule, SmilesParseError> {
        self.skip_ascii_whitespace();
        let mut mol = Molecule::new();
        let first_atom = self.parse_atom()?;
        let current = mol.add_atom(first_atom);
        self.parse_chain(&mut mol, current)?;
        while self.consume_if('.') {
            let fragment_atom = self.parse_atom()?;
            let fragment_idx = mol.add_atom(fragment_atom);
            self.parse_chain(&mut mol, fragment_idx)?;
        }
        if !self.ring_opens.is_empty() {
            return Err(self.error("unclosed ring"));
        }
        Ok(mol)
    }

    fn parse_chain(
        &mut self,
        mol: &mut Molecule,
        mut current_atom: usize,
    ) -> Result<(), SmilesParseError> {
        loop {
            self.skip_ascii_whitespace();
            if self.is_eof() || self.peek_char() == Some(')') || self.peek_char() == Some('.') {
                return Ok(());
            }

            if self.peek_char() == Some('(') {
                self.consume_char();
                let branch_bond = self.parse_optional_bond();
                let branch_atom = self.parse_atom()?;
                let branch_atom_idx = mol.add_atom(branch_atom);
                self.add_resolved_bond(
                    mol,
                    current_atom,
                    branch_atom_idx,
                    branch_bond.unwrap_or(PendingBond::Unspecified),
                );
                self.parse_chain(mol, branch_atom_idx)?;
                self.expect_char(')')?;
                continue;
            }

            if let Some(ring_number) = self.parse_optional_ring_number() {
                self.handle_ring_closure(mol, current_atom, PendingBond::Unspecified, ring_number)?;
                continue;
            }

            if let Some(pending_bond) = self.parse_optional_bond() {
                if let Some(ring_number) = self.parse_optional_ring_number() {
                    self.handle_ring_closure(mol, current_atom, pending_bond, ring_number)?;
                    continue;
                }

                let next_atom = self.parse_atom()?;
                let next_idx = mol.add_atom(next_atom);
                self.add_resolved_bond(mol, current_atom, next_idx, pending_bond);
                current_atom = next_idx;
                continue;
            }

            let next_atom = self.parse_atom()?;
            let next_idx = mol.add_atom(next_atom);
            self.add_resolved_bond(mol, current_atom, next_idx, PendingBond::Unspecified);
            current_atom = next_idx;
        }
    }

    fn handle_ring_closure(
        &mut self,
        mol: &mut Molecule,
        current_atom: usize,
        bond: PendingBond,
        ring_number: u32,
    ) -> Result<(), SmilesParseError> {
        if let Some(open) = self.ring_opens.remove(&ring_number) {
            if open.atom_index == current_atom {
                return Err(self.error("duplicated ring closure bonds atom to itself"));
            }
            if mol.bonds.iter().any(|b| {
                (b.begin_atom == open.atom_index && b.end_atom == current_atom)
                    || (b.begin_atom == current_atom && b.end_atom == open.atom_index)
            }) {
                return Err(self.error("ring closure duplicates existing bond"));
            }
            let chosen = if open.bond != PendingBond::Unspecified {
                open.bond
            } else {
                bond
            };
            self.add_resolved_bond(mol, open.atom_index, current_atom, chosen);
            Ok(())
        } else {
            self.ring_opens.insert(
                ring_number,
                RingOpen {
                    atom_index: current_atom,
                    bond,
                },
            );
            Ok(())
        }
    }

    fn add_resolved_bond(
        &self,
        mol: &mut Molecule,
        begin_atom: usize,
        end_atom: usize,
        pending: PendingBond,
    ) {
        let order = match pending {
            PendingBond::Single => BondOrder::Single,
            PendingBond::Double => BondOrder::Double,
            PendingBond::Triple => BondOrder::Triple,
            PendingBond::Quadruple => BondOrder::Quadruple,
            PendingBond::Aromatic => BondOrder::Aromatic,
            PendingBond::DirectionalSingle => BondOrder::Single,
            PendingBond::DativeForward | PendingBond::DativeBackward => BondOrder::Dative,
            PendingBond::Null => BondOrder::Null,
            PendingBond::Unspecified => {
                let atom1 = &mol.atoms[begin_atom];
                let atom2 = &mol.atoms[end_atom];
                if atom1.is_aromatic && atom2.is_aromatic {
                    BondOrder::Aromatic
                } else {
                    BondOrder::Single
                }
            }
        };
        mol.add_bond(Bond {
            index: 0,
            begin_atom,
            end_atom,
            order,
        });
    }

    fn parse_atom(&mut self) -> Result<Atom, SmilesParseError> {
        self.skip_ascii_whitespace();
        match self.peek_char() {
            Some('[') => self.parse_bracket_atom(),
            Some(_) => self.parse_simple_atom(),
            None => Err(self.error("expected atom")),
        }
    }

    fn parse_simple_atom(&mut self) -> Result<Atom, SmilesParseError> {
        if let Some((atomic_num, aromatic, consumed)) = self.match_simple_atom() {
            self.pos += consumed;
            return Ok(Atom {
                index: 0,
                atomic_num,
                is_aromatic: aromatic,
                formal_charge: 0,
                explicit_hydrogens: 0,
                isotope: None,
            });
        }
        Err(self.error("unsupported atom token"))
    }

    fn parse_bracket_atom(&mut self) -> Result<Atom, SmilesParseError> {
        self.expect_char('[')?;
        let isotope = self.parse_optional_number().map(|v| v as u16);
        let mut atom = if self.consume_if('H') {
            Atom {
                index: 0,
                atomic_num: 1,
                is_aromatic: false,
                formal_charge: 0,
                explicit_hydrogens: 0,
                isotope,
            }
        } else if self.consume_if('*') {
            Atom {
                index: 0,
                atomic_num: 0,
                is_aromatic: false,
                formal_charge: 0,
                explicit_hydrogens: 0,
                isotope,
            }
        } else if let Some((atomic_num, aromatic, consumed)) = self.match_bracket_atom_symbol() {
            self.pos += consumed;
            Atom {
                index: 0,
                atomic_num,
                is_aromatic: aromatic,
                formal_charge: 0,
                explicit_hydrogens: 0,
                isotope,
            }
        } else {
            return Err(self.error("unsupported bracket atom"));
        };

        if self.consume_if('@') {
            self.consume_if('@');
        }
        if self.consume_if('H') {
            atom.explicit_hydrogens = self.parse_optional_number().unwrap_or(1) as u8;
        }
        if self.consume_if('+') {
            atom.formal_charge = if self.consume_if('+') {
                2
            } else {
                self.parse_optional_number().unwrap_or(1) as i8
            };
        } else if self.consume_if('-') {
            atom.formal_charge = if self.consume_if('-') {
                -2
            } else {
                -(self.parse_optional_number().unwrap_or(1) as i8)
            };
        }
        if self.consume_if(':') {
            let _ = self.parse_required_number()?;
        }
        self.expect_char(']')?;
        Ok(atom)
    }

    fn parse_optional_bond(&mut self) -> Option<PendingBond> {
        match self.peek_char()? {
            '<' if self.input[self.pos..].starts_with("<-") => {
                self.pos += 2;
                Some(PendingBond::DativeBackward)
            }
            '-' if self.input[self.pos..].starts_with("->") => {
                self.pos += 2;
                Some(PendingBond::DativeForward)
            }
            '-' => {
                self.consume_char();
                Some(PendingBond::Single)
            }
            '=' => {
                self.consume_char();
                Some(PendingBond::Double)
            }
            '#' => {
                self.consume_char();
                Some(PendingBond::Triple)
            }
            '$' => {
                self.consume_char();
                Some(PendingBond::Quadruple)
            }
            ':' => {
                self.consume_char();
                Some(PendingBond::Aromatic)
            }
            '~' => {
                self.consume_char();
                Some(PendingBond::Null)
            }
            '/' => {
                self.consume_char();
                Some(PendingBond::DirectionalSingle)
            }
            '\\' => {
                self.consume_char();
                Some(PendingBond::DirectionalSingle)
            }
            _ => None,
        }
    }

    fn parse_optional_ring_number(&mut self) -> Option<u32> {
        if self.consume_if('%') {
            if self.consume_if('(') {
                let number = self.parse_required_number().ok()?;
                if !self.consume_if(')') {
                    return None;
                }
                return Some(number as u32);
            }
            let d1 = self.parse_single_digit()? as u32;
            let d2 = self.parse_single_digit()? as u32;
            return Some(d1 * 10 + d2);
        }
        self.parse_single_digit().map(|v| v as u32)
    }

    fn parse_optional_number(&mut self) -> Option<u32> {
        let start = self.pos;
        let mut value = 0u32;
        let mut consumed = false;
        while let Some(ch) = self.peek_char() {
            if let Some(digit) = ch.to_digit(10) {
                consumed = true;
                value = value.saturating_mul(10).saturating_add(digit);
                self.consume_char();
            } else {
                break;
            }
        }
        if consumed {
            Some(value)
        } else {
            self.pos = start;
            None
        }
    }

    fn parse_required_number(&mut self) -> Result<u32, SmilesParseError> {
        self.parse_optional_number()
            .ok_or_else(|| self.error("expected number"))
    }

    fn parse_single_digit(&mut self) -> Option<u8> {
        let ch = self.peek_char()?;
        if ch.is_ascii_digit() {
            self.consume_char();
            return Some((ch as u8) - b'0');
        }
        None
    }

    fn match_simple_atom(&self) -> Option<(u8, bool, usize)> {
        let rest = &self.input[self.pos..];
        for (token, atomic_num, aromatic) in [
            ("Cl", 17, false),
            ("Br", 35, false),
            ("*", 0, false),
            ("B", 5, false),
            ("C", 6, false),
            ("N", 7, false),
            ("O", 8, false),
            ("P", 15, false),
            ("S", 16, false),
            ("F", 9, false),
            ("I", 53, false),
            ("b", 5, true),
            ("c", 6, true),
            ("n", 7, true),
            ("o", 8, true),
            ("p", 15, true),
            ("s", 16, true),
        ] {
            if rest.starts_with(token) {
                return Some((atomic_num, aromatic, token.len()));
            }
        }
        None
    }

    fn match_bracket_atom_symbol(&self) -> Option<(u8, bool, usize)> {
        let rest = &self.input[self.pos..];
        for (token, atomic_num, aromatic) in [
            ("se", 34, true),
            ("as", 33, true),
            ("te", 52, true),
            ("si", 14, true),
            ("Rh", 45, false),
            ("Cu", 29, false),
            ("Cl", 17, false),
            ("Br", 35, false),
            ("Si", 14, false),
            ("As", 33, false),
            ("Se", 34, false),
            ("Li", 3, false),
            ("Na", 11, false),
            ("Mg", 12, false),
            ("Al", 13, false),
            ("Ca", 20, false),
            ("Fe", 26, false),
            ("Zn", 30, false),
            ("*", 0, false),
            ("b", 5, true),
            ("c", 6, true),
            ("n", 7, true),
            ("o", 8, true),
            ("p", 15, true),
            ("s", 16, true),
            ("B", 5, false),
            ("C", 6, false),
            ("N", 7, false),
            ("O", 8, false),
            ("P", 15, false),
            ("S", 16, false),
            ("F", 9, false),
            ("I", 53, false),
        ] {
            if rest.starts_with(token) {
                return Some((atomic_num, aromatic, token.len()));
            }
        }
        None
    }

    fn expect_char(&mut self, expected: char) -> Result<(), SmilesParseError> {
        if self.consume_if(expected) {
            Ok(())
        } else {
            Err(self.error("unexpected token"))
        }
    }

    fn consume_if(&mut self, expected: char) -> bool {
        if self.peek_char() == Some(expected) {
            self.consume_char();
            true
        } else {
            false
        }
    }

    fn consume_char(&mut self) -> Option<char> {
        let ch = self.peek_char()?;
        self.pos += ch.len_utf8();
        Some(ch)
    }

    fn peek_char(&self) -> Option<char> {
        self.input[self.pos..].chars().next()
    }

    fn skip_ascii_whitespace(&mut self) {
        while matches!(self.peek_char(), Some(ch) if ch.is_ascii_whitespace()) {
            self.consume_char();
        }
    }

    fn is_eof(&self) -> bool {
        self.pos >= self.input.len()
    }

    fn error(&self, message: &str) -> SmilesParseError {
        let _ = self.pos;
        SmilesParseError::ParseError(message.to_owned())
    }
}
