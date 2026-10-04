//! Private actual-site snapshots; no chemistry or updated valence cache.
use super::*;
use std::cell::RefCell;

#[derive(Debug, Clone, PartialEq)]
pub(super) struct Observation {
    pub stage: &'static str,
    pub topology: TopologyBlock,
    pub coordinates: CoordinateBlock,
    pub properties: MoleculeProperties,
    pub rings: Option<RingInfo>,
    pub coordinate_identity: Vec<(usize, usize, bool, Vec<u64>)>,
}

pub(super) fn coordinate_identity(
    coordinates: &CoordinateBlock,
) -> Vec<(usize, usize, bool, Vec<u64>)> {
    coordinates
        .conformers_2d
        .iter()
        .map(|c| {
            (
                c.id(),
                2,
                false,
                c.coordinates()
                    .iter()
                    .flatten()
                    .map(|x| x.to_bits())
                    .collect(),
            )
        })
        .chain(coordinates.conformers_3d.iter().map(|c| {
            (
                c.id(),
                3,
                c.is_3d(),
                c.coordinates()
                    .iter()
                    .flatten()
                    .map(|x| x.to_bits())
                    .collect(),
            )
        }))
        .collect()
}

thread_local! {
    static ACTIVE: RefCell<Option<Vec<Observation>>> = const { RefCell::new(None) };
}

pub(super) struct Guard;

impl Guard {
    pub(super) fn arm() -> Result<Self, &'static str> {
        ACTIVE.with(|active| {
            let mut active = active.borrow_mut();
            if active.is_some() {
                return Err("nested preparation stage observer");
            }
            *active = Some(Vec::new());
            Ok(Self)
        })
    }

    pub(super) fn finish(self) -> Vec<Observation> {
        ACTIVE.with(|active| active.borrow_mut().take().unwrap_or_default())
    }
}

impl Drop for Guard {
    fn drop(&mut self) {
        ACTIVE.with(|active| {
            active.borrow_mut().take();
        });
    }
}

pub(super) fn observe(
    stage: &'static str,
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
    rings: Option<&RingInfo>,
) {
    ACTIVE.with(|active| {
        let mut active = active.borrow_mut();
        if let Some(observations) = active.as_mut() {
            observations.push(Observation {
                stage,
                topology: topology.clone(),
                coordinates: coordinates.clone(),
                properties: properties.clone(),
                rings: rings.cloned(),
                coordinate_identity: coordinate_identity(coordinates),
            });
        }
    });
}
