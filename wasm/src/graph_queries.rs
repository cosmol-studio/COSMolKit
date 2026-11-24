use crate::{
    Molecule,
    canonical_errors::CanonicalRankError,
    fields::{self, Field},
    input,
    matrix_calls::DenseMatrix,
    matrix_errors::MatrixError,
    path_errors::PathError,
    storage::NativeField,
};
use cosmolkit as ck;
use serde::Serialize;
use std::collections::BTreeMap;
use tsify::{Ts, Tsify};
use wasm_bindgen::prelude::*;
ck::__ck_configuration_AtomRankParams!(crate::project_parameters, ck::AtomRankParams, {}, {
    AtomRankParams
});
ck::__ck_configuration_PathEnumerationParams!(
    crate::project_parameters,
    ck::PathEnumerationParams,
    {},
    { PathEnumerationParams }
);
ck::__ck_configuration_SubgraphEnumerationParams!(
    crate::project_parameters,
    ck::SubgraphEnumerationParams,
    {},
    { SubgraphEnumerationParams }
);
ck::__ck_configuration_UniqueSubgraphSearchParams!(
    crate::project_parameters,
    ck::UniqueSubgraphSearchParams,
    {},
    { UniqueSubgraphSearchParams }
);
ck::__ck_configuration_AtomEnvironmentParams!(
    crate::project_parameters,
    ck::AtomEnvironmentParams,
    {},
    { AtomEnvironmentParams }
);
#[derive(Serialize, Tsify)]
#[serde(transparent)]
pub struct IndexRows(#[tsify(type = "ReadonlyArray<ReadonlyArray<number>>")] Vec<Vec<usize>>);
#[derive(Serialize, Tsify)]
#[serde(transparent)]
pub struct IndexGroups(
    #[tsify(type = "ReadonlyArray<ReadonlyArray<ReadonlyArray<number>>>")] Vec<Vec<Vec<usize>>>,
);
#[derive(Serialize, Tsify)]
#[serde(transparent)]
pub struct AtomDistances(#[tsify(type = "ReadonlyMap<number,number>")] BTreeMap<usize, usize>);
#[wasm_bindgen]
pub struct AtomEnvironment(ck::AtomEnvironment);
#[wasm_bindgen]
impl AtomEnvironment {
    #[wasm_bindgen(getter)]
    pub fn bonds(&self) -> Result<Ts<fields::IndexArray>, JsValue> {
        Vec::<usize>::get(NativeField::owned(
            self.0.bonds().iter().map(|i| i.index()).collect(),
        ))
    }
    #[wasm_bindgen(getter,js_name=atomDistances)]
    pub fn atom_distances(&self) -> Result<Ts<AtomDistances>, JsValue> {
        AtomDistances(
            self.0
                .atom_distances()
                .iter()
                .map(|(i, d)| (i.index(), *d))
                .collect(),
        )
        .into_ts()
        .map_err(|e| wasm_bindgen::JsError::from(e).into())
    }
}
fn path_indices(path: ck::GraphPath) -> Vec<usize> {
    match path {
        ck::GraphPath::Atoms(v) => v.into_iter().map(|i| i.index()).collect(),
        ck::GraphPath::Bonds(v) => bond_indices(v),
    }
}
fn bond_indices(v: Vec<ck::BondId>) -> Vec<usize> {
    v.into_iter().map(|i| i.index()).collect()
}
#[wasm_bindgen]
impl Molecule {
    /// Canonical ranks in source atom order. Sort indices by rank to obtain a permutation.
    #[wasm_bindgen(js_name=canonicalAtomRanks)]
    pub fn canonical_atom_ranks(
        &self,
        params: Option<AtomRankParamsInput>,
    ) -> Result<Ts<fields::IndexArray>, JsValue> {
        AtomRankParams::with_input(params.as_ref(), |p| {
            let v = self
                .0
                .canonical_atom_ranks_with_params(p)
                .map_err(CanonicalRankError::from_rust)?;
            Vec::<usize>::get(NativeField::owned(v))
        })
    }

    #[wasm_bindgen(js_name=canonicalAtomRanksWithParams)]
    pub fn canonical_atom_ranks_with_params(
        &self,
        params: AtomRankParamsInput,
    ) -> Result<Ts<fields::IndexArray>, JsValue> {
        self.canonical_atom_ranks(Some(params))
    }

    /// Rank selected atoms using only the specified bonds. Unselected atoms have rank -1.
    #[wasm_bindgen(js_name=canonicalAtomRanksInFragment)]
    pub fn canonical_atom_ranks_in_fragment(
        &self,
        atoms: Ts<fields::IndexArray>,
        bonds: Option<Ts<fields::IndexArray>>,
        atom_symbols: Option<Ts<fields::TextArray>>,
        params: Option<AtomRankParamsInput>,
    ) -> Result<Ts<fields::FloatArray>, JsValue> {
        AtomRankParams::with_input(params.as_ref(), |p| {
            let v = self
                .0
                .canonical_atom_ranks_in_fragment_with_params(
                    &Vec::<usize>::input(atoms)?
                        .into_iter()
                        .map(ck::AtomId::new)
                        .collect::<Vec<_>>(),
                    &bonds
                        .map(Vec::<usize>::input)
                        .transpose()?
                        .unwrap_or_default()
                        .into_iter()
                        .map(ck::BondId::new)
                        .collect::<Vec<_>>(),
                    atom_symbols
                        .map(Vec::<String>::input)
                        .transpose()?
                        .as_deref(),
                    p,
                )
                .map_err(CanonicalRankError::from_rust)?;
            Vec::<f64>::get(NativeField::owned(
                v.into_iter().map(|rank| rank as f64).collect(),
            ))
        })
    }

    #[wasm_bindgen(js_name=canonicalAtomRanksInFragmentWithParams)]
    pub fn canonical_atom_ranks_in_fragment_with_params(
        &self,
        atoms: Ts<fields::IndexArray>,
        bonds: Option<Ts<fields::IndexArray>>,
        atom_symbols: Option<Ts<fields::TextArray>>,
        params: AtomRankParamsInput,
    ) -> Result<Ts<fields::FloatArray>, JsValue> {
        self.canonical_atom_ranks_in_fragment(atoms, bonds, atom_symbols, Some(params))
    }

    /// Enumerate linear paths. Length counts bonds by default, or atoms in atom-index mode.
    #[wasm_bindgen(js_name=pathsOfLength)]
    pub fn paths_of_length(
        &self,
        length: Ts<input::IndexInput>,
        params: Option<PathEnumerationParamsInput>,
    ) -> Result<Ts<IndexRows>, JsValue> {
        PathEnumerationParams::with_input(params.as_ref(), |p| {
            let v = self
                .0
                .paths_of_length_with_params(input::read(length)?, p)
                .map_err(PathError::from_rust)?;
            IndexRows(v.into_iter().map(path_indices).collect())
                .into_ts()
                .map_err(|e| wasm_bindgen::JsError::from(e).into())
        })
    }

    #[wasm_bindgen(js_name=pathsOfLengthWithParams)]
    pub fn paths_of_length_with_params(
        &self,
        length: Ts<input::IndexInput>,
        params: PathEnumerationParamsInput,
    ) -> Result<Ts<IndexRows>, JsValue> {
        self.paths_of_length(length, Some(params))
    }

    /// Enumerate connected sets of the specified number of bonds, including branched sets.
    #[wasm_bindgen(js_name=subgraphsOfLength)]
    pub fn subgraphs_of_length(
        &self,
        length: Ts<input::IndexInput>,
        params: Option<SubgraphEnumerationParamsInput>,
    ) -> Result<Ts<IndexRows>, JsValue> {
        SubgraphEnumerationParams::with_input(params.as_ref(), |p| {
            let v = self
                .0
                .subgraphs_of_length_with_params(input::read(length)?, p)
                .map_err(PathError::from_rust)?;
            IndexRows(v.into_iter().map(bond_indices).collect())
                .into_ts()
                .map_err(|e| wasm_bindgen::JsError::from(e).into())
        })
    }

    #[wasm_bindgen(js_name=subgraphsOfLengthWithParams)]
    pub fn subgraphs_of_length_with_params(
        &self,
        length: Ts<input::IndexInput>,
        params: SubgraphEnumerationParamsInput,
    ) -> Result<Ts<IndexRows>, JsValue> {
        self.subgraphs_of_length(length, Some(params))
    }

    /// Enumerate connected bond sets for each length from min to max, inclusive.
    #[wasm_bindgen(js_name=subgraphsInRange)]
    pub fn subgraphs_in_range(
        &self,
        min: Ts<input::IndexInput>,
        max: Ts<input::IndexInput>,
        params: Option<SubgraphEnumerationParamsInput>,
    ) -> Result<Ts<IndexGroups>, JsValue> {
        SubgraphEnumerationParams::with_input(params.as_ref(), |p| {
            let v = self
                .0
                .subgraphs_in_range_with_params(input::read(min)?, input::read(max)?, p)
                .map_err(PathError::from_rust)?;
            IndexGroups(
                v.into_iter()
                    .map(|group| group.into_iter().map(bond_indices).collect())
                    .collect(),
            )
            .into_ts()
            .map_err(|e| wasm_bindgen::JsError::from(e).into())
        })
    }

    #[wasm_bindgen(js_name=subgraphsInRangeWithParams)]
    pub fn subgraphs_in_range_with_params(
        &self,
        min: Ts<input::IndexInput>,
        max: Ts<input::IndexInput>,
        params: SubgraphEnumerationParamsInput,
    ) -> Result<Ts<IndexGroups>, JsValue> {
        self.subgraphs_in_range(min, max, Some(params))
    }

    /// Return the first representative of each distinct connected bond subgraph.
    #[wasm_bindgen(js_name=uniqueSubgraphsOfLength)]
    pub fn unique_subgraphs_of_length(
        &self,
        length: Ts<input::IndexInput>,
        params: Option<UniqueSubgraphSearchParamsInput>,
    ) -> Result<Ts<IndexRows>, JsValue> {
        UniqueSubgraphSearchParams::with_input(params.as_ref(), |p| {
            let v = self
                .0
                .unique_subgraphs_of_length_with_params(input::read(length)?, p)
                .map_err(PathError::from_rust)?;
            IndexRows(v.into_iter().map(bond_indices).collect())
                .into_ts()
                .map_err(|e| wasm_bindgen::JsError::from(e).into())
        })
    }

    #[wasm_bindgen(js_name=uniqueSubgraphsOfLengthWithParams)]
    pub fn unique_subgraphs_of_length_with_params(
        &self,
        length: Ts<input::IndexInput>,
        params: UniqueSubgraphSearchParamsInput,
    ) -> Result<Ts<IndexRows>, JsValue> {
        self.unique_subgraphs_of_length(length, Some(params))
    }

    /// Return source bonds and distances to atoms in a rooted environment.
    #[wasm_bindgen(js_name=atomEnvironment)]
    pub fn atom_environment(
        &self,
        radius: Ts<input::IndexInput>,
        root: Ts<input::IndexInput>,
        params: Option<AtomEnvironmentParamsInput>,
    ) -> Result<AtomEnvironment, JsValue> {
        AtomEnvironmentParams::with_input(params.as_ref(), |p| {
            let v = self
                .0
                .atom_environment_with_params(
                    input::read(radius)?,
                    ck::AtomId::new(input::read(root)?),
                    p,
                )
                .map_err(PathError::from_rust)?;
            Ok(AtomEnvironment(v))
        })
    }

    #[wasm_bindgen(js_name=atomEnvironmentWithParams)]
    pub fn atom_environment_with_params(
        &self,
        radius: Ts<input::IndexInput>,
        root: Ts<input::IndexInput>,
        params: AtomEnvironmentParamsInput,
    ) -> Result<AtomEnvironment, JsValue> {
        self.atom_environment(radius, root, Some(params))
    }

    /// Return the bond index, or undefined for two valid atoms without a bond.
    #[wasm_bindgen(js_name=bondBetweenAtoms)]
    pub fn bond_between_atoms(
        &self,
        begin: Ts<input::IndexInput>,
        end: Ts<input::IndexInput>,
    ) -> Result<Option<usize>, JsValue> {
        let v = self
            .0
            .bond_between_atoms(
                ck::AtomId::new(input::read(begin)?),
                ck::AtomId::new(input::read(end)?),
            )
            .map_err(PathError::from_rust)?;
        Ok(v.map(|i| i.index()))
    }

    /// Return a shortest path as atom indices; disconnected atoms give an empty path.
    #[wasm_bindgen(js_name=shortestPath)]
    pub fn shortest_path(
        &self,
        begin: Ts<input::IndexInput>,
        end: Ts<input::IndexInput>,
    ) -> Result<Ts<fields::IndexArray>, JsValue> {
        let v = self
            .0
            .shortest_path(
                ck::AtomId::new(input::read(begin)?),
                ck::AtomId::new(input::read(end)?),
            )
            .map_err(PathError::from_rust)?;
        Vec::<usize>::get(NativeField::owned(
            v.into_iter().map(|i| i.index()).collect(),
        ))
    }

    /// Return adjacency in atom order, optionally weighted by each bond's valence contribution to the row atom.
    #[wasm_bindgen(js_name=adjacencyMatrix)]
    pub fn adjacency_matrix(
        &self,
        use_bond_order: Option<Ts<input::BoolInput>>,
    ) -> Result<DenseMatrix, JsValue> {
        let v = self
            .0
            .adjacency_matrix(
                use_bond_order
                    .map(input::read)
                    .transpose()?
                    .unwrap_or(false),
            )
            .map_err(MatrixError::from_rust)?;
        Ok(DenseMatrix(v))
    }
}
