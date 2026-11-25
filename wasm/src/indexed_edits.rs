use crate::{
    Molecule, MoleculeInput,
    construction_values::{AtomSpec, AtomSpecInput, BondSpec, BondSpecInput},
    fields::{Field, IndexArray},
    input,
    operation_errors::OperationError,
    storage::NativeField,
};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::prelude::*;
ck::__ck_configuration_SubgraphExtractionParams!(
    crate::project_parameters,
    ck::SubgraphExtractionParams,
    {},
    { SubgraphExtractionParams }
);
#[wasm_bindgen]
pub struct SubgraphSelection(ck::SubgraphSelection);
#[wasm_bindgen]
pub struct SubgraphMapping(ck::SubgraphMapping);
#[wasm_bindgen]
pub struct MappedSubgraph(ck::MappedSubgraph);
crate::fields::value_field!(
    ck::SubgraphSelection,
    SubgraphSelection,
    SubgraphSelectionInput,
    "SubgraphSelection"
);
#[wasm_bindgen]
impl SubgraphSelection {
    pub fn atoms(indices: Ts<IndexArray>) -> Result<Self, JsValue> {
        Ok(Self(ck::SubgraphSelection::Atoms(
            Vec::<usize>::input(indices)?
                .into_iter()
                .map(ck::AtomId::new)
                .collect(),
        )))
    }
    pub fn bonds(indices: Ts<IndexArray>) -> Result<Self, JsValue> {
        Ok(Self(ck::SubgraphSelection::Bonds(
            Vec::<usize>::input(indices)?
                .into_iter()
                .map(ck::BondId::new)
                .collect(),
        )))
    }
    #[wasm_bindgen(js_name=atomsAndBonds)]
    pub fn atoms_and_bonds(atoms: Ts<IndexArray>, bonds: Ts<IndexArray>) -> Result<Self, JsValue> {
        Ok(Self(ck::SubgraphSelection::AtomsAndBonds(
            Vec::<usize>::input(atoms)?
                .into_iter()
                .map(ck::AtomId::new)
                .collect(),
            Vec::<usize>::input(bonds)?
                .into_iter()
                .map(ck::BondId::new)
                .collect(),
        )))
    }
}
#[wasm_bindgen]
impl MappedSubgraph {
    /// Each read yields an independently editable and releasable handle.
    #[wasm_bindgen(getter)]
    pub fn molecule(&self) -> Molecule {
        Molecule(self.0.molecule().clone())
    }
    #[wasm_bindgen(getter)]
    pub fn mapping(&self) -> SubgraphMapping {
        SubgraphMapping(self.0.mapping().clone())
    }
}
#[wasm_bindgen]
impl SubgraphMapping {
    #[wasm_bindgen(getter,js_name=sourceAtomIndices)]
    pub fn source_atom_indices(&self) -> Result<Ts<IndexArray>, JsValue> {
        Vec::<usize>::get(NativeField::owned(
            self.0
                .source_atom_indices()
                .iter()
                .map(|i| i.index())
                .collect(),
        ))
    }
    #[wasm_bindgen(getter,js_name=sourceBondIndices)]
    pub fn source_bond_indices(&self) -> Result<Ts<IndexArray>, JsValue> {
        Vec::<usize>::get(NativeField::owned(
            self.0
                .source_bond_indices()
                .iter()
                .map(|i| i.index())
                .collect(),
        ))
    }
}
#[wasm_bindgen]
impl Molecule {
    /// Apply a complete new-to-old permutation to atom and coordinate rows. Bonds retain their row order. Returns a new molecule; this molecule is unchanged.
    #[wasm_bindgen(js_name=withAtomsRenumbered)]
    pub fn with_atoms_renumbered(&self, order: Ts<IndexArray>) -> Result<Molecule, JsValue> {
        self.0
            .with_atoms_renumbered(
                &Vec::<usize>::input(order)?
                    .into_iter()
                    .map(ck::AtomId::new)
                    .collect::<Vec<_>>(),
            )
            .map(Molecule)
            .map_err(OperationError::from_rust)
    }

    /// Apply a complete new-to-old permutation to atom and coordinate rows. Bonds retain their row order. Changes this molecule.
    #[wasm_bindgen(js_name=renumberAtoms_)]
    pub fn renumber_atoms_(&mut self, order: Ts<IndexArray>) -> Result<(), JsValue> {
        self.0
            .renumber_atoms_(
                &Vec::<usize>::input(order)?
                    .into_iter()
                    .map(ck::AtomId::new)
                    .collect::<Vec<_>>(),
            )
            .map_err(OperationError::from_rust)
    }

    /// Remove an atom and its incident bonds; surviving rows stay in source order. Returns a new molecule; this molecule is unchanged.
    #[wasm_bindgen(js_name=withoutAtom)]
    pub fn without_atom(&self, atom: Ts<input::IndexInput>) -> Result<Molecule, JsValue> {
        self.0
            .without_atom(ck::AtomId::new(input::read(atom)?))
            .map(Molecule)
            .map_err(OperationError::from_rust)
    }

    /// Remove an atom and its incident bonds; surviving rows stay in source order. Changes this molecule.
    #[wasm_bindgen(js_name=removeAtom_)]
    pub fn remove_atom_(&mut self, atom: Ts<input::IndexInput>) -> Result<(), JsValue> {
        self.0
            .remove_atom_(ck::AtomId::new(input::read(atom)?))
            .map_err(OperationError::from_rust)
    }

    /// Replace one atom slot without moving its coordinates. Set preserve_props to use the original property dictionary. Returns a new molecule; this molecule is unchanged.
    #[wasm_bindgen(js_name=withAtomReplaced)]
    pub fn with_atom_replaced(
        &self,
        atom: Ts<input::IndexInput>,
        replacement: AtomSpecInput,
        preserve_props: Option<Ts<input::BoolInput>>,
    ) -> Result<Molecule, JsValue> {
        AtomSpec::with_ref(&replacement, |replacement| {
            self.0
                .with_atom_replaced(
                    ck::AtomId::new(input::read(atom)?),
                    replacement,
                    preserve_props
                        .map(input::read)
                        .transpose()?
                        .unwrap_or(false),
                )
                .map(Molecule)
                .map_err(OperationError::from_rust)
        })
    }

    /// Replace one atom slot without moving its coordinates. Set preserve_props to use the original property dictionary. Changes this molecule.
    #[wasm_bindgen(js_name=replaceAtom_)]
    pub fn replace_atom_(
        &mut self,
        atom: Ts<input::IndexInput>,
        replacement: AtomSpecInput,
        preserve_props: Option<Ts<input::BoolInput>>,
    ) -> Result<(), JsValue> {
        AtomSpec::with_ref(&replacement, |replacement| {
            self.0
                .replace_atom_(
                    ck::AtomId::new(input::read(atom)?),
                    replacement,
                    preserve_props
                        .map(input::read)
                        .transpose()?
                        .unwrap_or(false),
                )
                .map_err(OperationError::from_rust)
        })
    }

    /// Replace a bond, keeping its endpoints; supplied stereo neighbors belong to the target graph. Returns a new molecule; this molecule is unchanged.
    #[wasm_bindgen(js_name=withBondReplaced)]
    pub fn with_bond_replaced(
        &self,
        bond: Ts<input::IndexInput>,
        replacement: BondSpecInput,
        preserve_props: Option<Ts<input::BoolInput>>,
        keep_sgroups: Option<Ts<input::BoolInput>>,
    ) -> Result<Molecule, JsValue> {
        BondSpec::with_ref(&replacement, |replacement| {
            self.0
                .with_bond_replaced(
                    ck::BondId::new(input::read(bond)?),
                    replacement,
                    preserve_props
                        .map(input::read)
                        .transpose()?
                        .unwrap_or(false),
                    keep_sgroups.map(input::read).transpose()?.unwrap_or(true),
                )
                .map(Molecule)
                .map_err(OperationError::from_rust)
        })
    }

    /// Replace a bond, keeping its endpoints; supplied stereo neighbors belong to the target graph. Changes this molecule.
    #[wasm_bindgen(js_name=replaceBond_)]
    pub fn replace_bond_(
        &mut self,
        bond: Ts<input::IndexInput>,
        replacement: BondSpecInput,
        preserve_props: Option<Ts<input::BoolInput>>,
        keep_sgroups: Option<Ts<input::BoolInput>>,
    ) -> Result<(), JsValue> {
        BondSpec::with_ref(&replacement, |replacement| {
            self.0
                .replace_bond_(
                    ck::BondId::new(input::read(bond)?),
                    replacement,
                    preserve_props
                        .map(input::read)
                        .transpose()?
                        .unwrap_or(false),
                    keep_sgroups.map(input::read).transpose()?.unwrap_or(true),
                )
                .map_err(OperationError::from_rust)
        })
    }

    /// Append the other graph without an inter-molecular bond. Match coordinate sets by dimension and ID. Returns a new molecule; this molecule is unchanged.
    #[wasm_bindgen(js_name=withMoleculeAppended)]
    pub fn with_molecule_appended(&self, right: MoleculeInput) -> Result<Molecule, JsValue> {
        Molecule::with_ref(&right, |right| {
            self.0
                .with_molecule_appended(right)
                .map(Molecule)
                .map_err(OperationError::from_rust)
        })
    }

    /// Append the other graph without an inter-molecular bond. Match coordinate sets by dimension and ID. Changes this molecule.
    #[wasm_bindgen(js_name=insertMolecule_)]
    pub fn insert_molecule_(&mut self, right: MoleculeInput) -> Result<(), JsValue> {
        Molecule::with_ref(&right, |right| {
            self.0
                .insert_molecule_(right)
                .map_err(OperationError::from_rust)
        })
    }

    /// Return a new molecule and local-to-source mappings for the selected atoms and bonds. SGroups are omitted unless preserveSgroups is enabled; unrepresentable references cause an error.
    #[wasm_bindgen(js_name=extractSubgraph)]
    pub fn extract_subgraph(
        &self,
        selection: SubgraphSelectionInput,
        params: Option<SubgraphExtractionParamsInput>,
    ) -> Result<MappedSubgraph, JsValue> {
        SubgraphSelection::with_ref(&selection, |s| {
            SubgraphExtractionParams::with_input(params.as_ref(), |p| {
                self.0
                    .extract_subgraph_with_params(s, p)
                    .map(MappedSubgraph)
                    .map_err(OperationError::from_rust)
            })
        })
    }
    #[wasm_bindgen(js_name=extractSubgraphWithParams)]
    pub fn extract_subgraph_with_params(
        &self,
        selection: SubgraphSelectionInput,
        params: SubgraphExtractionParamsInput,
    ) -> Result<MappedSubgraph, JsValue> {
        self.extract_subgraph(selection, Some(params))
    }

    /// Retain the selected graph in this molecule; return its local-to-source mapping. SGroups are omitted unless preserveSgroups is enabled; unrepresentable references cause an error.
    #[wasm_bindgen(js_name=retainSubgraph_)]
    pub fn retain_subgraph_(
        &mut self,
        selection: SubgraphSelectionInput,
        params: Option<SubgraphExtractionParamsInput>,
    ) -> Result<SubgraphMapping, JsValue> {
        SubgraphSelection::with_ref(&selection, |s| {
            SubgraphExtractionParams::with_input(params.as_ref(), |p| {
                self.0
                    .retain_subgraph_with_params_(s, p)
                    .map(SubgraphMapping)
                    .map_err(OperationError::from_rust)
            })
        })
    }
    #[wasm_bindgen(js_name=retainSubgraphWithParams_)]
    pub fn retain_subgraph_with_params_(
        &mut self,
        selection: SubgraphSelectionInput,
        params: SubgraphExtractionParamsInput,
    ) -> Result<SubgraphMapping, JsValue> {
        self.retain_subgraph_(selection, Some(params))
    }
}
