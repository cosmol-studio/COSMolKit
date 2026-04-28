use crate::{ChiralTag, Molecule};

/// A ligand participating in a tetrahedral stereocenter.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum LigandRef {
    /// A graph atom by molecule atom index.
    Atom(usize),
    /// A hydrogen carried as an atom hydrogen count rather than a graph atom.
    ImplicitH,
}

/// Ordered tetrahedral stereochemistry for one atom center.
///
/// Specification: `tetrahedral_stereo_representation.md` at repository root.
///
/// For RDKit-style tetrahedral tags, the ordering is derived from the atom's
/// incoming bond order and chiral tag.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct TetrahedralStereo {
    pub center: usize,
    pub ligands: [LigandRef; 4],
}

impl Molecule {
    /// Return ordered tetrahedral stereocenters derived from RDKit-style atom
    /// chiral tags and molecule bond order.
    #[must_use]
    pub fn tetrahedral_stereo(&self) -> Vec<TetrahedralStereo> {
        let mut out = Vec::new();
        for atom in &self.atoms {
            let tag = atom.chiral_tag;
            if !matches!(tag, ChiralTag::TetrahedralCcw | ChiralTag::TetrahedralCw) {
                continue;
            }

            let mut ligands = Vec::with_capacity(4);
            for bond in &self.bonds {
                if bond.begin_atom == atom.index {
                    ligands.push(LigandRef::Atom(bond.end_atom));
                } else if bond.end_atom == atom.index {
                    ligands.push(LigandRef::Atom(bond.begin_atom));
                }
            }

            if ligands.len() == 3 {
                ligands.push(LigandRef::ImplicitH);
            }
            let Ok(mut ligands) = <[LigandRef; 4]>::try_from(ligands) else {
                continue;
            };

            // RDKit DistGeom treats CHI_TETRAHEDRAL_CCW as positive volume
            // for atom bond order, and CHI_TETRAHEDRAL_CW as negative. COSMolKit
            // stores ordered ligands with positive orientation, so CW is an odd
            // permutation of the RDKit atom-bond order.
            if matches!(tag, ChiralTag::TetrahedralCw) {
                ligands.swap(0, 1);
            }

            out.push(TetrahedralStereo {
                center: atom.index,
                ligands,
            });
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::{LigandRef, TetrahedralStereo};
    use crate::Molecule;

    #[test]
    fn tetrahedral_stereo_is_derived_from_chiral_tags() {
        let ccw = Molecule::from_smiles("F[C@](Cl)(Br)I").expect("parse chiral smiles");
        let cw = Molecule::from_smiles("F[C@@](Cl)(Br)I").expect("parse chiral smiles");

        let ccw_stereo = ccw.tetrahedral_stereo();
        let cw_stereo = cw.tetrahedral_stereo();

        assert_eq!(ccw_stereo.len(), 1);
        assert_eq!(cw_stereo.len(), 1);
        assert_eq!(ccw_stereo[0].center, 1);
        assert_eq!(cw_stereo[0].center, 1);
        assert_eq!(
            ccw_stereo[0].ligands,
            [
                LigandRef::Atom(0),
                LigandRef::Atom(2),
                LigandRef::Atom(3),
                LigandRef::Atom(4)
            ]
        );
        assert_eq!(
            cw_stereo[0].ligands,
            [
                LigandRef::Atom(2),
                LigandRef::Atom(0),
                LigandRef::Atom(3),
                LigandRef::Atom(4)
            ]
        );
    }

    #[test]
    fn tetrahedral_stereo_places_implicit_hydrogen_as_fourth_ligand() {
        let mol = Molecule::from_smiles("[13CH3:7][C@H](F)Cl").expect("parse chiral smiles");
        let stereo = mol.tetrahedral_stereo();

        assert_eq!(
            stereo,
            vec![TetrahedralStereo {
                center: 1,
                ligands: [
                    LigandRef::Atom(0),
                    LigandRef::Atom(2),
                    LigandRef::Atom(3),
                    LigandRef::ImplicitH
                ],
            }]
        );
    }
}
