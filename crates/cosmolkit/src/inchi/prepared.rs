//! Native detached generation input. No neutral-graph roundtrip.
use super::toolkit::{bond_type, chirality, direction, stereo};
use crate::{AtomId, CoordinateBlock, CoordinateSourceConformer, PropertyValue, TopologyBlock};
use cosmolkit_core as core;
use cosmolkit_inchi as ic;
use ic::rdkit::{BondInput, CleanupGraph, Graph, PreparedGraph};

fn error(value: impl std::fmt::Display) -> ic::InchiToolkitError {
    ic::InchiToolkitError {
        kind: "ToolkitError",
        message: value.to_string(),
    }
}

pub(super) struct Prepared<'a> {
    topology: TopologyBlock,
    conformer: Option<CoordinateSourceConformer<'a>>,
    valence: core::ValenceAssignment,
}

impl<'a> Prepared<'a> {
    pub(super) fn new(
        topology: &TopologyBlock,
        coordinates: &'a CoordinateBlock,
        cached: Option<&core::ValenceAssignment>,
    ) -> Result<Self, crate::InchiError> {
        // Preserve coordinate admission before the source toolkit callbacks.
        let conformer = coordinates
            .first_source_conformer()
            .map_err(error)
            .map_err(super::toolkit::public_error)?;
        let prepare = || -> Result<_, ic::InchiToolkitError> {
            u32::try_from(topology.atoms.len()).map_err(error)?;
            u32::try_from(topology.bonds.len()).map_err(error)?;
            // RDKit✔️✔️: if (mol.needsUpdatePropertyCache()) {
            // RDKit✔️✔️:   m->updatePropertyCache(false);
            // RDKit✔️✔️: }
            let needs_update = cached.is_none_or(|v| {
                topology.atoms.iter().enumerate().any(|(i, a)| {
                    v.explicit_valence[i] < 0 || (!a.no_implicit() && v.implicit_hydrogens[i] < 0)
                })
            });
            let valence = if needs_update {
                core::assign_valence(
                    topology,
                    &core::ValenceParams {
                        strict: false,
                        ..Default::default()
                    },
                )
                .map_err(error)?
            } else {
                cached.expect("cache admission above").clone()
            };
            // RDKit✔️✔️: MolOps::Kekulize(*m, false, false);
            // The owner returns detached working topology; no intermediary graph.
            let kekule = core::kekulize(
                topology,
                &core::KekulizeParams {
                    mark_atoms_bonds: false,
                    canonical: false,
                    ..Default::default()
                },
            )
            .map_err(error)?;
            Ok((kekule.topology, valence))
        };
        let (topology, valence) = prepare().map_err(|value| crate::InchiError {
            operation: "mol_to_inchi",
            kind: crate::InchiErrorKind::Toolkit,
            detail: value.message,
        })?;
        Ok(Self {
            topology,
            conformer,
            valence,
        })
    }
}

impl Graph for Prepared<'_> {
    fn atom_count(&self) -> usize {
        self.topology.atoms.len()
    }
    fn bond_count(&self) -> usize {
        self.topology.bonds.len()
    }
    fn atom(&self, index: usize) -> ic::InchiAtom {
        let a = &self.topology.atoms[index];
        ic::InchiAtom {
            atomic_number: a.atomic_number().into(),
            formal_charge: a.formal_charge().into(),
            num_explicit_hydrogens: a.explicit_hydrogens().into(),
            is_aromatic: a.is_aromatic(),
            isotope: a.isotope().unwrap_or(0).into(),
            num_radical_electrons: a.radical_electrons().into(),
            no_implicit: a.no_implicit(),
            chiral_tag: chirality(a.chiral_tag()),
            cip_rank: match a.prop("_CIPRank") {
                Some(PropertyValue::UInt(v)) => Some(*v),
                _ => None,
            },
            cached_explicit_valence: Some(self.valence.explicit_valence[index]),
            cached_implicit_valence: Some(self.valence.implicit_hydrogens[index]),
        }
    }
    fn bond(&self, index: usize) -> BondInput {
        let b = &self.topology.bonds[index];
        BondInput {
            begin_atom_index: b.begin().index() as u32,
            end_atom_index: b.end().index() as u32,
            bond_type: bond_type(b.order()),
            direction: direction(b.direction()),
            stereo: stereo(b.stereo()),
            stereo_atoms: b
                .stereo_atom_references()
                .get(..2)
                .map(|s| [s[0].index() as u32, s[1].index() as u32]),
        }
    }
    fn degree(&self, index: usize) -> usize {
        self.topology.adjacency.neighbors_of(index).len()
    }
    fn neighbor(&self, index: usize, offset: usize) -> (u32, u32) {
        let n = &self.topology.adjacency.neighbors_of(index)[offset];
        (n.atom_index as u32, n.bond.index() as u32)
    }
}

impl CleanupGraph for Prepared<'_> {
    fn set_charge(&mut self, index: usize, charge: i8) {
        self.topology.atoms[index].set_formal_charge(charge);
    }
    fn set_bond_type(&mut self, index: usize, order: ic::InchiBondType) {
        self.topology.bonds[index].set_order(
            crate::BondOrder::from_rdkit_code(order as i64).expect("shared bond vocabulary"),
        );
    }
}

impl PreparedGraph for Prepared<'_> {
    fn conformer_len(&self) -> Option<usize> {
        self.conformer.as_ref().map(|c| match c {
            CoordinateSourceConformer::TwoD(c) => c.coordinates().len(),
            CoordinateSourceConformer::ThreeD(c) => c.coordinates().len(),
        })
    }
    fn position(&self, index: usize) -> [f64; 3] {
        match self
            .conformer
            .as_ref()
            .expect("encoder checks conformer presence")
        {
            CoordinateSourceConformer::TwoD(c) => {
                let [x, y] = c.coordinates()[index];
                [x, y, 0.0]
            }
            CoordinateSourceConformer::ThreeD(c) => c.coordinates()[index],
        }
    }
    fn element_symbol(&mut self, number: i32) -> Result<Vec<u8>, ic::InchiToolkitError> {
        ic::MolToInchiToolkit::element_symbol(&mut super::toolkit::Toolkit::default(), number)
    }
    fn atomic_weight(&mut self, number: i32) -> Result<f64, ic::InchiToolkitError> {
        ic::MolToInchiToolkit::atomic_weight(&mut super::toolkit::Toolkit::default(), number)
    }
    fn total_num_hydrogens(&mut self, index: u32) -> Result<u32, ic::InchiToolkitError> {
        super::toolkit::total_hydrogens(&self.atom(index as usize))
    }
    fn calc_implicit_valence(&mut self, index: u32) -> Result<i32, ic::InchiToolkitError> {
        let value = core::implicit_valence_for_atom(
            &self.topology,
            AtomId::new(index as usize),
            Some(self.valence.explicit_valence[index as usize]),
            true,
        )
        .map_err(error)?;
        self.valence.implicit_hydrogens[index as usize] = value;
        Ok(value)
    }
    fn total_degree(&mut self, index: u32) -> Result<u32, ic::InchiToolkitError> {
        Ok(u32::try_from(self.degree(index as usize)).map_err(error)?
            + self.total_num_hydrogens(index)?)
    }
}
