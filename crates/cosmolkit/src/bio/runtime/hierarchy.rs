//! Owned hierarchy snapshots. Relationship resolution belongs to the facade;
//! detached rows and coordinate storage remain in the BIO value model.
use super::*;
use std::sync::Arc;

macro_rules! view {
    ($name:ident, $id:ident, $row:ident, $table:ident) => {
        #[doc = concat!("Read-only ", stringify!($name), " retaining its structural version. Later edits of the source use COW and do not retarget this object.")]
        #[derive(Debug, Clone)]
        pub struct $name {
            owner: Arc<BioStructure>,
            id: $id,
        }
        impl $name {
            pub fn id(&self) -> $id {
                self.id
            }
            /// Read the detached row within this retained structural version.
            pub fn value(&self) -> &$row {
                &self.owner.data.$table()[self.id.index()]
            }
        }
        impl std::ops::Deref for $name {
            type Target = $row;
            fn deref(&self) -> &Self::Target {
                self.value()
            }
        }
    };
}

view!(BioModel, BioModelId, BioModelRow, models);
view!(BioChain, BioChainId, BioChainRow, chains);
view!(BioResidue, BioResidueId, BioResidueRow, residues);
view!(BioAtom, BioAtomId, BioAtomRow, atoms);

impl BioModel {
    pub fn source_model_number(&self) -> Option<i32> {
        self.value().source_model_number()
    }
    pub fn chains(&self) -> Vec<BioChain> {
        let span = self.chain_span();
        (span.start()..span.end())
            .map(|index| BioChain {
                owner: Arc::clone(&self.owner),
                id: BioChainId::new(index),
            })
            .collect()
    }
}
impl BioChain {
    pub fn kind(&self) -> ChainKind {
        self.value().kind()
    }
    pub fn source(&self) -> &ChainSourceIds {
        self.value().source()
    }
    pub fn residues(&self) -> Vec<BioResidue> {
        let span = self.residue_span();
        (span.start()..span.end())
            .map(|index| BioResidue {
                owner: Arc::clone(&self.owner),
                id: BioResidueId::new(index),
            })
            .collect()
    }
}
impl BioResidue {
    pub fn name(&self) -> ResidueName {
        self.value().name()
    }
    pub fn kind(&self) -> ResidueKind {
        self.value().kind()
    }
    pub fn atoms(&self) -> Vec<BioAtom> {
        let span = self.atom_span();
        (span.start()..span.end())
            .map(|index| BioAtom {
                owner: Arc::clone(&self.owner),
                id: BioAtomId::new(index),
            })
            .collect()
    }
}
impl BioAtom {
    pub fn name(&self) -> AtomName {
        self.value().name()
    }
    pub fn element(&self) -> Element {
        self.value().element()
    }
    pub fn occupancy(&self) -> f64 {
        self.value().occupancy()
    }
    pub fn b_iso(&self) -> f64 {
        self.value().b_iso()
    }
    pub fn position(&self) -> [f64; 3] {
        // A published BioStructure has coordinates for every atom. Both the
        // row and coordinates are retained by this same validated snapshot.
        self.owner
            .data
            .atom_position(self.id)
            .expect("validated BIO atom")
    }
}

macro_rules! collection {
    ($method:ident, $single:ident, $view:ident, $id:ident) => {
        impl BioStructure {
            /// Retain one shared structural snapshot for all returned objects.
            pub fn $method(&self) -> Vec<$view> {
                let owner = Arc::new(self.clone());
                (0..self.data.$method().len())
                    .map(|index| $view {
                        owner: Arc::clone(&owner),
                        id: $id::new(index as u32),
                    })
                    .collect()
            }
            /// Retain the selected row's structural version without copying blocks.
            pub fn $single(&self, id: $id) -> Option<$view> {
                self.data.$method().get(id.index())?;
                Some($view {
                    owner: Arc::new(self.clone()),
                    id,
                })
            }
        }
    };
}
collection!(models, model, BioModel, BioModelId);
collection!(chains, chain, BioChain, BioChainId);
collection!(residues, residue, BioResidue, BioResidueId);
collection!(atoms, atom, BioAtom, BioAtomId);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn owned_hierarchy_retains_rows_coordinates_and_cow_snapshot() {
        let pdb = concat!(
            "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \n",
            "HETATM    2  O   HOH B   2       4.000   5.000   6.000  1.00 10.00           O  \n",
        );
        let mut structure = BioStructure::from_pdb(pdb).unwrap();
        let model = structure.models().remove(0);
        let chains = model.chains();
        assert_eq!(chains.len(), 2);
        let residue = chains[1].residues().remove(0);
        assert_eq!(residue.name().as_str(), "HOH");
        let atom = residue.atoms().remove(0);
        assert_eq!(atom.id(), BioAtomId::new(1));
        assert_eq!(atom.residue_id(), residue.id());
        assert_eq!(atom.position(), [4.0, 5.0, 6.0]);
        assert_eq!(structure.atoms()[1].value(), atom.value());
        assert!(Arc::ptr_eq(&model.owner, &atom.owner));
        assert!(Arc::ptr_eq(&structure.data.atoms, &atom.owner.data.atoms));
        assert!(Arc::ptr_eq(
            &structure.data.coordinates,
            &atom.owner.data.coordinates
        ));
        structure.translate_([10.0, 0.0, 0.0]).unwrap();
        assert_eq!(
            structure.atom(BioAtomId::new(1)).unwrap().position(),
            [14.0, 5.0, 6.0]
        );
        assert!(Arc::ptr_eq(&structure.data.atoms, &atom.owner.data.atoms));
        assert!(!Arc::ptr_eq(
            &structure.data.coordinates,
            &atom.owner.data.coordinates
        ));
        let selection = BioSelection::from_cid("//A").unwrap();
        structure.retain_selection_(&selection).unwrap();
        assert_eq!(structure.num_atoms(), 1);
        assert!(structure.atom(BioAtomId::new(1)).is_none());
        drop(structure);
        assert_eq!(atom.position(), [4.0, 5.0, 6.0]);
        assert_eq!(
            model.chains()[1].residues()[0].atoms()[0].position(),
            atom.position()
        );
    }
}
impl BioModel {
    pub fn num_chains(&self) -> usize {
        self.chain_span().len() as usize
    }
}
impl BioChain {
    pub fn num_residues(&self) -> usize {
        self.residue_span().len() as usize
    }
}
impl BioResidue {
    pub fn num_atoms(&self) -> usize {
        self.atom_span().len() as usize
    }
}
impl BioAtom {
    pub fn atomic_number(&self) -> u8 {
        self.element().atomic_number()
    }
}
impl BioModel {
    pub fn chain_span(&self) -> crate::BioRowSpan<crate::BioChainId> {
        self.value().chain_span()
    }
}
impl BioChain {
    pub fn model_id(&self) -> crate::BioModelId {
        self.value().model_id()
    }
}
impl BioChain {
    pub fn entity_id(&self) -> Option<crate::BioEntityId> {
        self.value().entity_id()
    }
}
impl BioChain {
    pub fn residue_span(&self) -> crate::BioRowSpan<crate::BioResidueId> {
        self.value().residue_span()
    }
}
impl BioResidue {
    pub fn chain_id(&self) -> crate::BioChainId {
        self.value().chain_id()
    }
}
impl BioResidue {
    pub fn atom_span(&self) -> crate::BioRowSpan<crate::BioAtomId> {
        self.value().atom_span()
    }
}
impl BioResidue {
    pub fn residue_info_kind(&self) -> crate::ResidueInfoKind {
        self.value().residue_info_kind()
    }
}
impl BioResidue {
    pub fn entity_kind(&self) -> crate::EntityKind {
        self.value().entity_kind()
    }
}
impl BioResidue {
    pub fn entity_id(&self) -> Option<crate::BioEntityId> {
        self.value().entity_id()
    }
}
impl BioResidue {
    pub fn het_flag(&self) -> Option<u8> {
        self.value().het_flag()
    }
}
impl BioResidue {
    pub fn source(&self) -> &crate::ResidueSourceIds {
        self.value().source()
    }
}
impl BioResidue {
    pub fn sifts_unp(&self) -> crate::BioSiftsUnpResidue {
        self.value().sifts_unp()
    }
}
impl BioAtom {
    pub fn residue_id(&self) -> crate::BioResidueId {
        self.value().residue_id()
    }
}
impl BioAtom {
    pub fn isotope_mass_number(&self) -> Option<u16> {
        self.value().isotope_mass_number()
    }
}
impl BioAtom {
    pub fn altloc(&self) -> Option<crate::AltLocLabel> {
        self.value().altloc()
    }
}
impl BioAtom {
    pub fn formal_charge(&self) -> i8 {
        self.value().formal_charge()
    }
}
impl BioAtom {
    pub fn calculation_flag(&self) -> crate::BioCalcFlag {
        self.value().calc_flag()
    }
}
impl BioAtom {
    pub fn anisou(&self) -> &[f64; 6] {
        self.value().anisou()
    }
}
impl BioAtom {
    pub fn tls_group_id(&self) -> i16 {
        self.value().tls_group_id()
    }
}
impl BioAtom {
    pub fn fraction(&self) -> f64 {
        self.value().fraction()
    }
}
impl BioAtom {
    pub fn source(&self) -> &crate::AtomSourceIds {
        self.value().source()
    }
}
