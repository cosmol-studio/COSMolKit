//! Invariant checks for BioStructure.

use crate::bio::BioStructure;

/// Returns `Err(&'static str)` describing the first violated invariant.
pub fn enforce_bio_structure_invariants(s: &BioStructure) -> Result<(), &'static str> {
    // coordinates.len() == atoms.len()
    if s.coordinates.positions.len() != s.atoms.len() {
        return Err("coordinates.len() != atoms.len()");
    }

    // residue atom_span validity
    for residue in &s.residues {
        let span = residue.atom_span;
        if span.end() as usize > s.atoms.len() {
            return Err("residue atom_span exceeds atoms block length");
        }
    }

    // chain residue_span validity
    for chain in &s.chains {
        let span = chain.residue_span;
        if span.end() as usize > s.residues.len() {
            return Err("chain residue_span exceeds residues block length");
        }
        if let Some(entity_id) = chain.entity_id
            && entity_id.index() as usize >= s.entities.len()
        {
            return Err("chain.entity_id out of range");
        }
    }

    // model chain_span validity
    for model in &s.models {
        let span = model.chain_span;
        if span.end() as usize > s.chains.len() {
            return Err("model chain_span exceeds chains block length");
        }
    }

    // atom.residue_id agrees with the residue span that contains it
    for (atom_idx, atom) in s.atoms.iter().enumerate() {
        let rid = atom.residue_id.index() as usize;
        if rid >= s.residues.len() {
            return Err("atom.residue_id out of range");
        }
        let span = s.residues[rid].atom_span;
        if (atom_idx as u32) < span.start || (atom_idx as u32) >= span.end() {
            return Err("atom.residue_id does not match the residue span containing this atom");
        }
    }

    // residue.chain_id agrees with the chain span that contains it
    for (res_idx, residue) in s.residues.iter().enumerate() {
        let cid = residue.chain_id.index() as usize;
        if cid >= s.chains.len() {
            return Err("residue.chain_id out of range");
        }
        let span = s.chains[cid].residue_span;
        if (res_idx as u32) < span.start || (res_idx as u32) >= span.end() {
            return Err("residue.chain_id does not match the chain span containing this residue");
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bio::*;

    fn empty_atom(residue_id: ResidueId) -> AtomRow {
        AtomRow {
            residue_id,
            name: AtomName([b' ', b'C', b'A', b' ']),
            element: crate::Element::C,
            altloc: None,
            occupancy: None,
            b_iso: None,
            formal_charge: None,
            anisou: None,
            source: AtomSourceIds { serial: None },
        }
    }

    fn empty_residue(chain_id: ChainId, start: u32, len: u32) -> ResidueRow {
        ResidueRow {
            chain_id,
            atom_span: RowSpan::new(start, len),
            name: ResidueName([b'A', b'L', b'A', 0], 3),
            kind: ResidueKind::AminoAcid,
            source: ResidueSourceIds { seq_id: None },
        }
    }

    fn empty_chain(model_id: ModelId, start: u32, len: u32) -> ChainRow {
        ChainRow {
            model_id,
            entity_id: None,
            residue_span: RowSpan::new(start, len),
            kind: ChainKind::Protein,
            source: ChainSourceIds {
                auth_chain_id: None,
                label_asym_id: None,
            },
        }
    }

    #[test]
    fn empty_structure_is_valid() {
        assert!(enforce_bio_structure_invariants(&BioStructure::new()).is_ok());
    }

    #[test]
    fn coordinate_mismatch_is_rejected() {
        let mut s = BioStructure::new();
        s.coordinates.positions.push([0.0, 0.0, 0.0]);
        assert!(enforce_bio_structure_invariants(&s).is_err());
    }

    #[test]
    fn valid_single_atom_structure_passes() {
        let mut s = BioStructure::new();
        s.models.push(ModelRow {
            chain_span: RowSpan::new(0, 1),
            source_model_number: Some(1),
        });
        s.chains.push(empty_chain(ModelId::new(0), 0, 1));
        s.residues.push(empty_residue(ChainId::new(0), 0, 1));
        s.atoms.push(empty_atom(ResidueId::new(0)));
        s.coordinates.positions.push([1.0, 2.0, 3.0]);
        assert!(enforce_bio_structure_invariants(&s).is_ok());
    }
}
