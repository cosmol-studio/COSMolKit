//! Model transport only. Chemistry is delegated to existing detached owners.
use crate::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, CoordinateBlock, Element, MoleculeProperties,
    PropertyValue, TopologyBlock,
};
use cosmolkit_core as core;
use cosmolkit_inchi as ic;
use ic::{InchiMolecule, InchiToolkitError as Error};

#[derive(Default)]
pub(super) struct Toolkit {
    pub(super) rings: Option<core::RingInfo>,
    pub(super) final_topology: Option<TopologyBlock>,
}

fn error(value: impl std::fmt::Display) -> Error {
    Error {
        kind: "ToolkitError",
        message: value.to_string(),
    }
}

pub(super) fn public_error(value: Error) -> crate::InchiError {
    crate::InchiError {
        operation: "inchi_model",
        kind: crate::InchiErrorKind::Toolkit,
        detail: value.message,
    }
}

fn element(number: i32) -> Result<Element, Error> {
    u8::try_from(number)
        .ok()
        .and_then(Element::from_atomic_number)
        .ok_or_else(|| error(format!("invalid atomic number: {number}")))
}

pub(super) fn model(graph: &InchiMolecule) -> Result<(TopologyBlock, CoordinateBlock), Error> {
    let atoms = graph
        .atoms()
        .iter()
        .enumerate()
        .map(|(i, atom)| {
            let mut spec = AtomSpec::new(element(atom.atomic_number)?)
                .with_formal_charge(i8::try_from(atom.formal_charge).map_err(error)?)
                .with_explicit_hydrogens(u8::try_from(atom.num_explicit_hydrogens).map_err(error)?)
                .with_aromatic(atom.is_aromatic)
                .with_isotope(u16::try_from(atom.isotope).map_err(error)?)
                .with_radical_electrons(u8::try_from(atom.num_radical_electrons).map_err(error)?)
                .with_no_implicit(atom.no_implicit)
                .with_chiral_tag(
                    crate::ChiralTag::from_rdkit_code(atom.chiral_tag as i64)
                        .ok_or_else(|| error("invalid chiral tag"))?,
                );
            if let Some(rank) = atom.cip_rank {
                spec = spec
                    .with_computed_prop("_CIPRank", PropertyValue::UInt(rank))
                    .map_err(error)?;
            }
            Ok(Atom::from_spec(AtomId::new(i), spec))
        })
        .collect::<Result<Vec<_>, Error>>()?;
    let bonds = graph
        .bonds()
        .iter()
        .enumerate()
        .map(|(i, bond)| {
            let mut spec = BondSpec::new(
                AtomId::new(bond.begin_atom_index() as usize),
                AtomId::new(bond.end_atom_index() as usize),
                crate::BondOrder::from_rdkit_code(bond.bond_type as i64)
                    .ok_or_else(|| error("invalid bond type"))?,
            )
            .with_aromatic(bond.is_aromatic)
            .with_direction(
                crate::BondDirection::from_rdkit_code(bond.direction as i64)
                    .ok_or_else(|| error("invalid bond direction"))?,
            )
            .with_stereo(
                crate::BondStereo::from_rdkit_code(bond.stereo as i64)
                    .ok_or_else(|| error("invalid bond stereo"))?,
            );
            if let [a, b] = bond.stereo_atoms.as_slice() {
                spec = spec.with_stereo_atoms(AtomId::new(*a as usize), AtomId::new(*b as usize));
            }
            let mut row = Bond::from_spec(BondId::new(i), spec);
            row.set_source_stereo_atom_references(
                bond.stereo_atoms
                    .iter()
                    .map(|&a| AtomId::new(a as usize))
                    .collect(),
            );
            Ok(row)
        })
        .collect::<Result<Vec<_>, Error>>()?;
    let topology =
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).map_err(error)?;
    let coordinates = CoordinateBlock {
        conformers_3d: graph
            .conformers()
            .iter()
            .enumerate()
            .map(|(id, rows)| crate::Conformer3D::new(id, rows.clone(), true))
            .collect(),
        ..Default::default()
    };
    Ok((topology, coordinates))
}

pub(super) fn cached_valence(graph: &InchiMolecule) -> Option<core::ValenceAssignment> {
    Some(core::ValenceAssignment {
        explicit_valence: graph
            .atoms()
            .iter()
            .map(|a| a.cached_explicit_valence)
            .collect::<Option<Vec<_>>>()?,
        implicit_hydrogens: graph
            .atoms()
            .iter()
            .map(|a| a.cached_implicit_valence)
            .collect::<Option<Vec<_>>>()?,
    })
}

pub(super) fn graph(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    valence: Option<&core::ValenceAssignment>,
) -> Result<InchiMolecule, Error> {
    // Only the source-selected FIRST conformer participates in MolToInchi.
    // Preserve its source ordering; no 3D preference or coordinate generation.
    let conformers = match coordinates.first_source_conformer().map_err(error)? {
        None => Vec::new(),
        Some(crate::CoordinateSourceConformer::TwoD(c)) => {
            vec![c.coordinates().iter().map(|p| [p[0], p[1], 0.0]).collect()]
        }
        Some(crate::CoordinateSourceConformer::ThreeD(c)) => vec![c.coordinates().to_vec()],
    };
    let atoms = topology
        .atoms
        .iter()
        .enumerate()
        .map(|(i, a)| ic::InchiAtom {
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
            cached_explicit_valence: valence.map(|v| v.explicit_valence[i]),
            cached_implicit_valence: valence.map(|v| v.implicit_hydrogens[i]),
        })
        .collect();
    let bonds = topology
        .bonds
        .iter()
        .map(|b| {
            let mut value = ic::InchiBond::new(
                u32::try_from(b.begin().index()).map_err(error)?,
                u32::try_from(b.end().index()).map_err(error)?,
                bond_type(b.order()),
            );
            value.direction = direction(b.direction());
            value.stereo = stereo(b.stereo());
            value.is_aromatic = b.is_aromatic();
            value.stereo_atoms = b
                .stereo_atom_references()
                .iter()
                .map(|a| u32::try_from(a.index()).map_err(error))
                .collect::<Result<_, _>>()?;
            Ok(value)
        })
        .collect::<Result<Vec<_>, Error>>()?;
    InchiMolecule::try_from_graph(atoms, bonds, conformers).map_err(|e| error(format!("{e:?}")))
}

fn install(
    graph_value: &mut InchiMolecule,
    topology: &TopologyBlock,
    valence: Option<&core::ValenceAssignment>,
) -> Result<(), Error> {
    // Topology-only callbacks retain all original coordinate rows.
    let updated = graph(topology, &CoordinateBlock::default(), valence)?;
    let (atoms, bonds, _) = updated.into_graph_parts();
    graph_value
        .replace_topology(atoms, bonds)
        .map_err(|e| error(format!("{e:?}")))
}

fn assignment(graph: &InchiMolecule) -> Result<core::ValenceAssignment, Error> {
    cached_valence(graph).ok_or_else(|| error("InChI callback requires the source property cache"))
}

impl ic::MolToInchiToolkit for Toolkit {
    fn needs_update_property_cache(&mut self, graph: &InchiMolecule) -> Result<bool, Error> {
        Ok(graph.atoms().iter().any(|a| {
            a.cached_explicit_valence.is_none_or(|v| v < 0)
                || (!a.no_implicit && a.cached_implicit_valence.is_none_or(|v| v < 0))
        }))
    }
    fn update_property_cache(
        &mut self,
        graph: &mut InchiMolecule,
        strict: bool,
    ) -> Result<(), Error> {
        let (topology, _) = model(graph)?;
        let value = core::assign_valence(
            &topology,
            &core::ValenceParams {
                strict,
                ..Default::default()
            },
        )
        .map_err(error)?;
        for (i, atom) in graph.atom_properties_mut().iter_mut().enumerate() {
            atom.cached_explicit_valence = Some(value.explicit_valence[i]);
            atom.cached_implicit_valence = Some(value.implicit_hydrogens[i]);
        }
        Ok(())
    }
    fn kekulize(&mut self, graph: &mut InchiMolecule, mark_atoms_bonds: bool) -> Result<(), Error> {
        let (topology, _) = model(graph)?;
        let value = core::kekulize(
            &topology,
            &core::KekulizeParams {
                mark_atoms_bonds,
                ..Default::default()
            },
        )
        .map_err(error)?;
        let valence = cached_valence(graph);
        install(graph, &value.topology, valence.as_ref())
    }
    fn element_symbol(&mut self, number: i32) -> Result<Vec<u8>, Error> {
        Ok(element(number)?.symbol().as_bytes().to_vec())
    }
    fn atomic_weight(&mut self, number: i32) -> Result<f64, Error> {
        Ok(core::element_info(element(number)?).atomic_weight)
    }
    fn total_num_hydrogens(&mut self, graph: &InchiMolecule, index: u32) -> Result<u32, Error> {
        let a = graph
            .atoms()
            .get(index as usize)
            .ok_or_else(|| error("atom index out of range"))?;
        // RDKit✔️✔️: int res = getNumExplicitHs() + getNumImplicitHs();
        // RDKit✔️✔️: if (df_noImplicit) {
        // RDKit✔️✔️:   return 0;
        // RDKit✔️✔️: }
        // Atom.cpp:286-307; direct scalar access preserves O(1) cost.
        let implicit = if a.no_implicit {
            0
        } else {
            u32::try_from(
                a.cached_implicit_valence
                    .ok_or_else(|| error("implicit valence is not initialized"))?,
            )
            .map_err(error)?
        };
        Ok(a.num_explicit_hydrogens + implicit)
    }
    fn calc_implicit_valence(
        &mut self,
        graph: &mut InchiMolecule,
        index: u32,
    ) -> Result<i32, Error> {
        let (topology, _) = model(graph)?;
        let explicit = graph
            .atoms()
            .get(index as usize)
            .ok_or_else(|| error("atom index out of range"))?
            .cached_explicit_valence;
        let value =
            core::implicit_valence_for_atom(&topology, AtomId::new(index as usize), explicit, true)
                .map_err(error)?;
        graph
            .atom_properties_mut()
            .get_mut(index as usize)
            .ok_or_else(|| error("atom index out of range"))?
            .cached_implicit_valence = Some(value);
        Ok(value)
    }
    fn total_degree(&mut self, graph: &InchiMolecule, index: u32) -> Result<u32, Error> {
        let degree = graph
            .atom_degree(index)
            .ok_or_else(|| error("atom index out of range"))?;
        Ok(u32::try_from(degree).map_err(error)? + self.total_num_hydrogens(graph, index)?)
    }
}

impl ic::InchiToMolToolkit for Toolkit {
    fn atomic_number(&mut self, symbol: &[u8]) -> Result<i32, Error> {
        let symbol = std::str::from_utf8(symbol).map_err(error)?;
        Element::from_symbol(symbol)
            .map(|e| i32::from(e.atomic_number()))
            .ok_or_else(|| error("unknown element"))
    }
    fn average_atomic_weight(&mut self, number: i32) -> Result<f64, Error> {
        ic::MolToInchiToolkit::atomic_weight(self, number)
    }
    fn update_property_cache(
        &mut self,
        graph: &mut InchiMolecule,
        strict: bool,
    ) -> Result<(), Error> {
        ic::MolToInchiToolkit::update_property_cache(self, graph, strict)
    }
    fn assign_atom_cip_ranks(&mut self, graph: &mut InchiMolecule) -> Result<Vec<u32>, Error> {
        let (topology, _) = model(graph)?;
        let ranks = core::assign_atom_cip_ranks(&topology, &assignment(graph)?).map_err(error)?;
        for (atom, &rank) in graph.atom_properties_mut().iter_mut().zip(&ranks) {
            atom.cip_rank = Some(rank);
        }
        Ok(ranks)
    }
    fn remove_hydrogens(&mut self, graph_value: &mut InchiMolecule) -> Result<(), Error> {
        let (topology, coordinates) = model(graph_value)?;
        let value =
            core::remove_hydrogens_impl(topology, coordinates, MoleculeProperties::default())
                .map_err(error)?;
        *graph_value = graph(
            &value.topology,
            &value.coordinates,
            value.final_valence.as_ref(),
        )?;
        self.rings = value.final_rings;
        Ok(())
    }
    fn sanitize_molecule(&mut self, graph: &mut InchiMolecule) -> Result<(), Error> {
        let (topology, _) = model(graph)?;
        let value =
            core::sanitize_topology(&topology, &core::SanitizeParams::default()).map_err(|e| {
                Error {
                    kind: "MolSanitizeException",
                    message: e.to_string(),
                }
            })?;
        install(graph, &value.topology, value.final_valence.as_ref())?;
        self.rings = value.final_rings;
        Ok(())
    }
    fn assign_stereochemistry(
        &mut self,
        graph: &mut InchiMolecule,
        clean_it: bool,
        force: bool,
    ) -> Result<(), Error> {
        // The InChI owner calls assignStereochemistry(mol, true, true).
        // There is no molecule-level done flag in its neutral graph.
        if !force {
            return Err(error("neutral InChI stereo callback requires force=true"));
        }
        let (topology, _) = model(graph)?;
        let valence = assignment(graph)?;
        let rings = match self.rings.take() {
            Some(rings) => rings,
            None => {
                core::find_sssr(&topology, &core::RingSearchParams::default()).map_err(error)?
            }
        };
        let value = core::assign_legacy_stereochemistry_with_assignments(
            topology, &valence, &rings, clean_it, false,
        )
        .map_err(error)?;
        install(graph, &value.topology, Some(&valence))?;
        self.rings = value.ring_update.or(Some(rings));
        self.final_topology = Some(value.topology);
        Ok(())
    }
}

fn bond_type(value: crate::BondOrder) -> ic::InchiBondType {
    match value {
        crate::BondOrder::Unspecified => ic::InchiBondType::Unspecified,
        crate::BondOrder::Single => ic::InchiBondType::Single,
        crate::BondOrder::Double => ic::InchiBondType::Double,
        crate::BondOrder::Triple => ic::InchiBondType::Triple,
        crate::BondOrder::Quadruple => ic::InchiBondType::Quadruple,
        crate::BondOrder::Quintuple => ic::InchiBondType::Quintuple,
        crate::BondOrder::Hextuple => ic::InchiBondType::Hextuple,
        crate::BondOrder::OneAndHalf => ic::InchiBondType::OneAndAHalf,
        crate::BondOrder::TwoAndHalf => ic::InchiBondType::TwoAndAHalf,
        crate::BondOrder::ThreeAndHalf => ic::InchiBondType::ThreeAndAHalf,
        crate::BondOrder::FourAndHalf => ic::InchiBondType::FourAndAHalf,
        crate::BondOrder::FiveAndHalf => ic::InchiBondType::FiveAndAHalf,
        crate::BondOrder::Aromatic => ic::InchiBondType::Aromatic,
        crate::BondOrder::Ionic => ic::InchiBondType::Ionic,
        crate::BondOrder::Hydrogen => ic::InchiBondType::Hydrogen,
        crate::BondOrder::ThreeCenter => ic::InchiBondType::ThreeCenter,
        crate::BondOrder::DativeOne => ic::InchiBondType::DativeOne,
        crate::BondOrder::Dative => ic::InchiBondType::Dative,
        crate::BondOrder::DativeLeft => ic::InchiBondType::DativeL,
        crate::BondOrder::DativeRight => ic::InchiBondType::DativeR,
        crate::BondOrder::Other => ic::InchiBondType::Other,
        crate::BondOrder::Zero => ic::InchiBondType::Zero,
    }
}

fn direction(value: crate::BondDirection) -> ic::InchiBondDirection {
    match value {
        crate::BondDirection::None => ic::InchiBondDirection::None,
        crate::BondDirection::BeginWedge => ic::InchiBondDirection::BeginWedge,
        crate::BondDirection::BeginDash => ic::InchiBondDirection::BeginDash,
        crate::BondDirection::EndDownRight => ic::InchiBondDirection::EndDownRight,
        crate::BondDirection::EndUpRight => ic::InchiBondDirection::EndUpRight,
        crate::BondDirection::EitherDouble => ic::InchiBondDirection::EitherDouble,
        crate::BondDirection::Unknown => ic::InchiBondDirection::Unknown,
    }
}

fn stereo(value: crate::BondStereo) -> ic::InchiBondStereo {
    match value {
        crate::BondStereo::None => ic::InchiBondStereo::None,
        crate::BondStereo::Any => ic::InchiBondStereo::Any,
        crate::BondStereo::Z => ic::InchiBondStereo::Z,
        crate::BondStereo::E => ic::InchiBondStereo::E,
        crate::BondStereo::Cis => ic::InchiBondStereo::Cis,
        crate::BondStereo::Trans => ic::InchiBondStereo::Trans,
        crate::BondStereo::AtropCw => ic::InchiBondStereo::AtropCw,
        crate::BondStereo::AtropCcw => ic::InchiBondStereo::AtropCcw,
    }
}

fn chirality(value: crate::ChiralTag) -> ic::InchiChiralTag {
    match value {
        crate::ChiralTag::Unspecified => ic::InchiChiralTag::Unspecified,
        crate::ChiralTag::TetrahedralCw => ic::InchiChiralTag::TetrahedralCw,
        crate::ChiralTag::TetrahedralCcw => ic::InchiChiralTag::TetrahedralCcw,
        crate::ChiralTag::Other => ic::InchiChiralTag::Other,
        crate::ChiralTag::Tetrahedral => ic::InchiChiralTag::Tetrahedral,
        crate::ChiralTag::Allene => ic::InchiChiralTag::Allene,
        crate::ChiralTag::SquarePlanar => ic::InchiChiralTag::SquarePlanar,
        crate::ChiralTag::TrigonalBipyramidal => ic::InchiChiralTag::TrigonalBipyramidal,
        crate::ChiralTag::Octahedral => ic::InchiChiralTag::Octahedral,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ic::MolToInchiToolkit;

    #[test]
    fn no_implicit_hydrogens_does_not_require_an_implicit_cache() {
        let mut graph = InchiMolecule::try_from_graph(
            vec![ic::InchiAtom {
                atomic_number: 6,
                no_implicit: true,
                num_explicit_hydrogens: 4,
                ..Default::default()
            }],
            Vec::new(),
            Vec::new(),
        )
        .unwrap();
        let mut toolkit = Toolkit::default();
        assert_eq!(toolkit.total_num_hydrogens(&graph, 0).unwrap(), 4);
        graph.atom_properties_mut()[0].no_implicit = false;
        assert!(toolkit.total_num_hydrogens(&graph, 0).is_err());
        graph.atom_properties_mut()[0].cached_implicit_valence = Some(2);
        assert_eq!(toolkit.total_num_hydrogens(&graph, 0).unwrap(), 6);
        assert_eq!(toolkit.total_degree(&graph, 0).unwrap(), 6);
    }
}
