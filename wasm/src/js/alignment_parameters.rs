//! Mutable alignment parameter projections of the canonical public facade.
//! Constructor defaults and fields match python/src/alignment_binding.rs.
//! Only host-value conversion lives here; algorithms remain in cosmolkit.
use cosmolkit_wasm::rust as ck;
use js_sys::{Array, ArrayBuffer, Reflect};
use wasm_bindgen::prelude::*;

use crate::host_values::*;
fn weights_value(value: &JsValue) -> Result<Option<Vec<f64>>, JsValue> {
    if value.is_null() || value.is_undefined() {
        return Ok(None);
    }
    sequence(value, "weights")?
        .iter()
        .map(|v| v.as_f64().ok_or_else(|| type_error("weights element")))
        .collect::<Result<Vec<_>, _>>()
        .map(Some)
}
fn indices_value(value: &JsValue, name: &str) -> Result<Option<Vec<usize>>, JsValue> {
    if value.is_null() || value.is_undefined() {
        return Ok(None);
    }
    sequence(value, name)?
        .iter()
        .map(|v| usize_value(&v, name))
        .collect::<Result<Vec<_>, _>>()
        .map(Some)
}
fn atom_map_value(value: &JsValue) -> Result<Option<Vec<ck::AlignmentAtomMap>>, JsValue> {
    if value.is_null() || value.is_undefined() {
        return Ok(None);
    }
    sequence(value, "atomMap")?
        .iter()
        .map(|v| {
            let probe = usize_value(&Reflect::get(&v, &"probeAtom".into())?, "probeAtom")?;
            let reference =
                usize_value(&Reflect::get(&v, &"referenceAtom".into())?, "referenceAtom")?;
            Ok(ck::AlignmentAtomMap::new(probe, reference))
        })
        .collect::<Result<Vec<_>, JsValue>>()
        .map(Some)
}
fn atom_maps_value(value: &JsValue) -> Result<Vec<Vec<ck::AlignmentAtomMap>>, JsValue> {
    if value.is_null() || value.is_undefined() {
        return Ok(Vec::new());
    }
    sequence(value, "atomMaps")?
        .iter()
        .map(|v| atom_map_value(&v)?.ok_or_else(|| type_error("atomMaps row")))
        .collect()
}
fn weights_js(value: &Option<Vec<f64>>) -> JsValue {
    match value {
        None => JsValue::NULL,
        Some(v) => v
            .iter()
            .map(|n| JsValue::from_f64(*n))
            .collect::<Array>()
            .into(),
    }
}
fn indices_js(value: &Option<Vec<usize>>) -> JsValue {
    match value {
        None => JsValue::NULL,
        Some(v) => v
            .iter()
            .map(|n| JsValue::from_f64(*n as f64))
            .collect::<Array>()
            .into(),
    }
}
fn atom_map_js(value: &Option<Vec<ck::AlignmentAtomMap>>) -> JsValue {
    match value {
        None => JsValue::NULL,
        Some(v) => v
            .iter()
            .map(|inner| JsValue::from(AlignmentAtomMap { inner: *inner }))
            .collect::<Array>()
            .into(),
    }
}
fn atom_maps_js(value: &[Vec<ck::AlignmentAtomMap>]) -> JsValue {
    value
        .iter()
        .map(|v| atom_map_js(&Some(v.clone())))
        .collect::<Array>()
        .into()
}

#[wasm_bindgen]
pub struct AlignmentAtomMap {
    pub(crate) inner: ck::AlignmentAtomMap,
}

#[wasm_bindgen]
impl AlignmentAtomMap {
    #[wasm_bindgen(constructor)]
    pub fn constructor(
        #[wasm_bindgen(unchecked_param_type = "number")] probe_atom: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number")] reference_atom: JsValue,
    ) -> Result<AlignmentAtomMap, JsValue> {
        Ok(Self {
            inner: ck::AlignmentAtomMap::new(
                usize_value(&probe_atom, "probeAtom")?,
                usize_value(&reference_atom, "referenceAtom")?,
            ),
        })
    }
    #[wasm_bindgen(js_name = new)]
    pub fn new(
        #[wasm_bindgen(unchecked_param_type = "number")] probe_atom: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number")] reference_atom: JsValue,
    ) -> Result<AlignmentAtomMap, JsValue> {
        Ok(Self {
            inner: ck::AlignmentAtomMap::new(
                usize_value(&probe_atom, "probeAtom")?,
                usize_value(&reference_atom, "referenceAtom")?,
            ),
        })
    }
    #[wasm_bindgen(getter, js_name = probeAtom)]
    pub fn probe_atom(&self) -> usize {
        self.inner.probe_atom
    }
    #[wasm_bindgen(setter, js_name = probeAtom)]
    pub fn set_probe_atom(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] probe_atom: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.probe_atom = usize_value(&probe_atom, "probeAtom")?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = referenceAtom)]
    pub fn reference_atom(&self) -> usize {
        self.inner.reference_atom
    }
    #[wasm_bindgen(setter, js_name = referenceAtom)]
    pub fn set_reference_atom(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] reference_atom: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.reference_atom = usize_value(&reference_atom, "referenceAtom")?;
        Ok(())
    }
}

#[wasm_bindgen]
pub struct AlignmentParameters {
    pub(crate) inner: ck::AlignmentParameters,
}

#[wasm_bindgen]
impl AlignmentParameters {
    #[wasm_bindgen(constructor)]
    pub fn constructor(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] probe_conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] reference_conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "AlignmentAtomMap[] | null")]
        atom_map: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Float64Array | null")]
        weights: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] reflect: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_iterations: JsValue,
    ) -> Result<AlignmentParameters, JsValue> {
        let mut inner = ck::AlignmentParameters::default();
        if !probe_conformer_id.is_undefined() {
            inner.probe_conformer_id = i32_value(&probe_conformer_id, "probeConformerId")?;
        }
        if !reference_conformer_id.is_undefined() {
            inner.reference_conformer_id =
                i32_value(&reference_conformer_id, "referenceConformerId")?;
        }
        if !atom_map.is_undefined() {
            inner.atom_map = atom_map_value(&atom_map)?;
        }
        if !weights.is_undefined() {
            inner.weights = weights_value(&weights)?;
        }
        if !reflect.is_undefined() {
            inner.reflect = bool_value(&reflect, "reflect")?;
        }
        if !max_iterations.is_undefined() {
            inner.max_iterations = u32_value(&max_iterations, "maxIterations")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(js_name = new)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] probe_conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] reference_conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "AlignmentAtomMap[] | null")]
        atom_map: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Float64Array | null")]
        weights: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] reflect: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_iterations: JsValue,
    ) -> Result<AlignmentParameters, JsValue> {
        let mut inner = ck::AlignmentParameters::default();
        if !probe_conformer_id.is_undefined() {
            inner.probe_conformer_id = i32_value(&probe_conformer_id, "probeConformerId")?;
        }
        if !reference_conformer_id.is_undefined() {
            inner.reference_conformer_id =
                i32_value(&reference_conformer_id, "referenceConformerId")?;
        }
        if !atom_map.is_undefined() {
            inner.atom_map = atom_map_value(&atom_map)?;
        }
        if !weights.is_undefined() {
            inner.weights = weights_value(&weights)?;
        }
        if !reflect.is_undefined() {
            inner.reflect = bool_value(&reflect, "reflect")?;
        }
        if !max_iterations.is_undefined() {
            inner.max_iterations = u32_value(&max_iterations, "maxIterations")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter, js_name = probeConformerId)]
    pub fn probe_conformer_id(&self) -> i32 {
        self.inner.probe_conformer_id
    }
    #[wasm_bindgen(setter, js_name = probeConformerId)]
    pub fn set_probe_conformer_id(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] probe_conformer_id: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.probe_conformer_id = i32_value(&probe_conformer_id, "probeConformerId")?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = referenceConformerId)]
    pub fn reference_conformer_id(&self) -> i32 {
        self.inner.reference_conformer_id
    }
    #[wasm_bindgen(setter, js_name = referenceConformerId)]
    pub fn set_reference_conformer_id(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] reference_conformer_id: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.reference_conformer_id =
            i32_value(&reference_conformer_id, "referenceConformerId")?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = atomMap, unchecked_return_type = "AlignmentAtomMap[] | null")]
    pub fn atom_map(&self) -> JsValue {
        atom_map_js(&self.inner.atom_map)
    }
    #[wasm_bindgen(setter, js_name = atomMap)]
    pub fn set_atom_map(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "AlignmentAtomMap[] | null")] atom_map: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.atom_map = atom_map_value(&atom_map)?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = weights, unchecked_return_type = "number[] | null")]
    pub fn weights(&self) -> JsValue {
        weights_js(&self.inner.weights)
    }
    #[wasm_bindgen(setter, js_name = weights)]
    pub fn set_weights(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number[] | Float64Array | null")] weights: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.weights = weights_value(&weights)?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = reflect)]
    pub fn reflect(&self) -> bool {
        self.inner.reflect
    }
    #[wasm_bindgen(setter, js_name = reflect)]
    pub fn set_reflect(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] reflect: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.reflect = bool_value(&reflect, "reflect")?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = maxIterations)]
    pub fn max_iterations(&self) -> u32 {
        self.inner.max_iterations
    }
    #[wasm_bindgen(setter, js_name = maxIterations)]
    pub fn set_max_iterations(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] max_iterations: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.max_iterations = u32_value(&max_iterations, "maxIterations")?;
        Ok(())
    }
}

#[wasm_bindgen]
pub struct BestAlignmentParameters {
    pub(crate) inner: ck::BestAlignmentParameters,
}

#[wasm_bindgen]
impl BestAlignmentParameters {
    #[wasm_bindgen(constructor)]
    pub fn constructor(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] probe_conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] reference_conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "AlignmentAtomMap[][] | null")]
        atom_maps: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Float64Array | null")]
        weights: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] reflect: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_iterations: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_matches: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        symmetrize_conjugated_terminal_groups: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] ignore_hydrogens: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] num_threads: JsValue,
    ) -> Result<BestAlignmentParameters, JsValue> {
        let mut inner = ck::BestAlignmentParameters::default();
        if !probe_conformer_id.is_undefined() {
            inner.probe_conformer_id = i32_value(&probe_conformer_id, "probeConformerId")?;
        }
        if !reference_conformer_id.is_undefined() {
            inner.reference_conformer_id =
                i32_value(&reference_conformer_id, "referenceConformerId")?;
        }
        if !atom_maps.is_undefined() {
            inner.atom_maps = atom_maps_value(&atom_maps)?;
        }
        if !weights.is_undefined() {
            inner.weights = weights_value(&weights)?;
        }
        if !reflect.is_undefined() {
            inner.reflect = bool_value(&reflect, "reflect")?;
        }
        if !max_iterations.is_undefined() {
            inner.max_iterations = u32_value(&max_iterations, "maxIterations")?;
        }
        if !max_matches.is_undefined() {
            inner.max_matches = i32_value(&max_matches, "maxMatches")?;
        }
        if !symmetrize_conjugated_terminal_groups.is_undefined() {
            inner.symmetrize_conjugated_terminal_groups = bool_value(
                &symmetrize_conjugated_terminal_groups,
                "symmetrizeConjugatedTerminalGroups",
            )?;
        }
        if !ignore_hydrogens.is_undefined() {
            inner.ignore_hydrogens = bool_value(&ignore_hydrogens, "ignoreHydrogens")?;
        }
        if !num_threads.is_undefined() {
            inner.num_threads = i32_value(&num_threads, "numThreads")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(js_name = new)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] probe_conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] reference_conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "AlignmentAtomMap[][] | null")]
        atom_maps: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Float64Array | null")]
        weights: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] reflect: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_iterations: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_matches: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        symmetrize_conjugated_terminal_groups: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] ignore_hydrogens: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] num_threads: JsValue,
    ) -> Result<BestAlignmentParameters, JsValue> {
        let mut inner = ck::BestAlignmentParameters::default();
        if !probe_conformer_id.is_undefined() {
            inner.probe_conformer_id = i32_value(&probe_conformer_id, "probeConformerId")?;
        }
        if !reference_conformer_id.is_undefined() {
            inner.reference_conformer_id =
                i32_value(&reference_conformer_id, "referenceConformerId")?;
        }
        if !atom_maps.is_undefined() {
            inner.atom_maps = atom_maps_value(&atom_maps)?;
        }
        if !weights.is_undefined() {
            inner.weights = weights_value(&weights)?;
        }
        if !reflect.is_undefined() {
            inner.reflect = bool_value(&reflect, "reflect")?;
        }
        if !max_iterations.is_undefined() {
            inner.max_iterations = u32_value(&max_iterations, "maxIterations")?;
        }
        if !max_matches.is_undefined() {
            inner.max_matches = i32_value(&max_matches, "maxMatches")?;
        }
        if !symmetrize_conjugated_terminal_groups.is_undefined() {
            inner.symmetrize_conjugated_terminal_groups = bool_value(
                &symmetrize_conjugated_terminal_groups,
                "symmetrizeConjugatedTerminalGroups",
            )?;
        }
        if !ignore_hydrogens.is_undefined() {
            inner.ignore_hydrogens = bool_value(&ignore_hydrogens, "ignoreHydrogens")?;
        }
        if !num_threads.is_undefined() {
            inner.num_threads = i32_value(&num_threads, "numThreads")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter, js_name = probeConformerId)]
    pub fn probe_conformer_id(&self) -> i32 {
        self.inner.probe_conformer_id
    }
    #[wasm_bindgen(setter, js_name = probeConformerId)]
    pub fn set_probe_conformer_id(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] probe_conformer_id: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.probe_conformer_id = i32_value(&probe_conformer_id, "probeConformerId")?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = referenceConformerId)]
    pub fn reference_conformer_id(&self) -> i32 {
        self.inner.reference_conformer_id
    }
    #[wasm_bindgen(setter, js_name = referenceConformerId)]
    pub fn set_reference_conformer_id(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] reference_conformer_id: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.reference_conformer_id =
            i32_value(&reference_conformer_id, "referenceConformerId")?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = atomMaps, unchecked_return_type = "AlignmentAtomMap[][]")]
    pub fn atom_maps(&self) -> JsValue {
        atom_maps_js(&self.inner.atom_maps)
    }
    #[wasm_bindgen(setter, js_name = atomMaps)]
    pub fn set_atom_maps(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "AlignmentAtomMap[][] | null")] atom_maps: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.atom_maps = atom_maps_value(&atom_maps)?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = weights, unchecked_return_type = "number[] | null")]
    pub fn weights(&self) -> JsValue {
        weights_js(&self.inner.weights)
    }
    #[wasm_bindgen(setter, js_name = weights)]
    pub fn set_weights(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number[] | Float64Array | null")] weights: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.weights = weights_value(&weights)?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = reflect)]
    pub fn reflect(&self) -> bool {
        self.inner.reflect
    }
    #[wasm_bindgen(setter, js_name = reflect)]
    pub fn set_reflect(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] reflect: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.reflect = bool_value(&reflect, "reflect")?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = maxIterations)]
    pub fn max_iterations(&self) -> u32 {
        self.inner.max_iterations
    }
    #[wasm_bindgen(setter, js_name = maxIterations)]
    pub fn set_max_iterations(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] max_iterations: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.max_iterations = u32_value(&max_iterations, "maxIterations")?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = maxMatches)]
    pub fn max_matches(&self) -> i32 {
        self.inner.max_matches
    }
    #[wasm_bindgen(setter, js_name = maxMatches)]
    pub fn set_max_matches(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] max_matches: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.max_matches = i32_value(&max_matches, "maxMatches")?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = symmetrizeConjugatedTerminalGroups)]
    pub fn symmetrize_conjugated_terminal_groups(&self) -> bool {
        self.inner.symmetrize_conjugated_terminal_groups
    }
    #[wasm_bindgen(setter, js_name = symmetrizeConjugatedTerminalGroups)]
    pub fn set_symmetrize_conjugated_terminal_groups(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")]
        symmetrize_conjugated_terminal_groups: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.symmetrize_conjugated_terminal_groups = bool_value(
            &symmetrize_conjugated_terminal_groups,
            "symmetrizeConjugatedTerminalGroups",
        )?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = ignoreHydrogens)]
    pub fn ignore_hydrogens(&self) -> bool {
        self.inner.ignore_hydrogens
    }
    #[wasm_bindgen(setter, js_name = ignoreHydrogens)]
    pub fn set_ignore_hydrogens(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] ignore_hydrogens: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.ignore_hydrogens = bool_value(&ignore_hydrogens, "ignoreHydrogens")?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = numThreads)]
    pub fn num_threads(&self) -> i32 {
        self.inner.num_threads
    }
    #[wasm_bindgen(setter, js_name = numThreads)]
    pub fn set_num_threads(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] num_threads: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.num_threads = i32_value(&num_threads, "numThreads")?;
        Ok(())
    }
}

#[wasm_bindgen]
pub struct CoordinateRmsdParameters {
    pub(crate) inner: ck::CoordinateRmsdParameters,
}

#[wasm_bindgen]
impl CoordinateRmsdParameters {
    #[wasm_bindgen(constructor)]
    pub fn constructor(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] probe_conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] reference_conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "AlignmentAtomMap[][] | null")]
        atom_maps: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Float64Array | null")]
        weights: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_matches: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        symmetrize_conjugated_terminal_groups: JsValue,
    ) -> Result<CoordinateRmsdParameters, JsValue> {
        let mut inner = ck::CoordinateRmsdParameters::default();
        if !probe_conformer_id.is_undefined() {
            inner.probe_conformer_id = i32_value(&probe_conformer_id, "probeConformerId")?;
        }
        if !reference_conformer_id.is_undefined() {
            inner.reference_conformer_id =
                i32_value(&reference_conformer_id, "referenceConformerId")?;
        }
        if !atom_maps.is_undefined() {
            inner.atom_maps = atom_maps_value(&atom_maps)?;
        }
        if !weights.is_undefined() {
            inner.weights = weights_value(&weights)?;
        }
        if !max_matches.is_undefined() {
            inner.max_matches = i32_value(&max_matches, "maxMatches")?;
        }
        if !symmetrize_conjugated_terminal_groups.is_undefined() {
            inner.symmetrize_conjugated_terminal_groups = bool_value(
                &symmetrize_conjugated_terminal_groups,
                "symmetrizeConjugatedTerminalGroups",
            )?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(js_name = new)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] probe_conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] reference_conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "AlignmentAtomMap[][] | null")]
        atom_maps: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Float64Array | null")]
        weights: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_matches: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        symmetrize_conjugated_terminal_groups: JsValue,
    ) -> Result<CoordinateRmsdParameters, JsValue> {
        let mut inner = ck::CoordinateRmsdParameters::default();
        if !probe_conformer_id.is_undefined() {
            inner.probe_conformer_id = i32_value(&probe_conformer_id, "probeConformerId")?;
        }
        if !reference_conformer_id.is_undefined() {
            inner.reference_conformer_id =
                i32_value(&reference_conformer_id, "referenceConformerId")?;
        }
        if !atom_maps.is_undefined() {
            inner.atom_maps = atom_maps_value(&atom_maps)?;
        }
        if !weights.is_undefined() {
            inner.weights = weights_value(&weights)?;
        }
        if !max_matches.is_undefined() {
            inner.max_matches = i32_value(&max_matches, "maxMatches")?;
        }
        if !symmetrize_conjugated_terminal_groups.is_undefined() {
            inner.symmetrize_conjugated_terminal_groups = bool_value(
                &symmetrize_conjugated_terminal_groups,
                "symmetrizeConjugatedTerminalGroups",
            )?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter, js_name = probeConformerId)]
    pub fn probe_conformer_id(&self) -> i32 {
        self.inner.probe_conformer_id
    }
    #[wasm_bindgen(setter, js_name = probeConformerId)]
    pub fn set_probe_conformer_id(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] probe_conformer_id: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.probe_conformer_id = i32_value(&probe_conformer_id, "probeConformerId")?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = referenceConformerId)]
    pub fn reference_conformer_id(&self) -> i32 {
        self.inner.reference_conformer_id
    }
    #[wasm_bindgen(setter, js_name = referenceConformerId)]
    pub fn set_reference_conformer_id(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] reference_conformer_id: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.reference_conformer_id =
            i32_value(&reference_conformer_id, "referenceConformerId")?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = atomMaps, unchecked_return_type = "AlignmentAtomMap[][]")]
    pub fn atom_maps(&self) -> JsValue {
        atom_maps_js(&self.inner.atom_maps)
    }
    #[wasm_bindgen(setter, js_name = atomMaps)]
    pub fn set_atom_maps(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "AlignmentAtomMap[][] | null")] atom_maps: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.atom_maps = atom_maps_value(&atom_maps)?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = weights, unchecked_return_type = "number[] | null")]
    pub fn weights(&self) -> JsValue {
        weights_js(&self.inner.weights)
    }
    #[wasm_bindgen(setter, js_name = weights)]
    pub fn set_weights(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number[] | Float64Array | null")] weights: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.weights = weights_value(&weights)?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = maxMatches)]
    pub fn max_matches(&self) -> i32 {
        self.inner.max_matches
    }
    #[wasm_bindgen(setter, js_name = maxMatches)]
    pub fn set_max_matches(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] max_matches: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.max_matches = i32_value(&max_matches, "maxMatches")?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = symmetrizeConjugatedTerminalGroups)]
    pub fn symmetrize_conjugated_terminal_groups(&self) -> bool {
        self.inner.symmetrize_conjugated_terminal_groups
    }
    #[wasm_bindgen(setter, js_name = symmetrizeConjugatedTerminalGroups)]
    pub fn set_symmetrize_conjugated_terminal_groups(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")]
        symmetrize_conjugated_terminal_groups: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.symmetrize_conjugated_terminal_groups = bool_value(
            &symmetrize_conjugated_terminal_groups,
            "symmetrizeConjugatedTerminalGroups",
        )?;
        Ok(())
    }
}

#[wasm_bindgen]
pub struct AllConformerRmsdParameters {
    pub(crate) inner: ck::AllConformerRmsdParameters,
}

#[wasm_bindgen]
impl AllConformerRmsdParameters {
    #[wasm_bindgen(constructor)]
    pub fn constructor(
        #[wasm_bindgen(unchecked_optional_param_type = "AlignmentAtomMap[][] | null")]
        atom_maps: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Float64Array | null")]
        weights: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_matches: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        symmetrize_conjugated_terminal_groups: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] ignore_hydrogens: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] num_threads: JsValue,
    ) -> Result<AllConformerRmsdParameters, JsValue> {
        let mut inner = ck::AllConformerRmsdParameters::default();
        if !atom_maps.is_undefined() {
            inner.atom_maps = atom_maps_value(&atom_maps)?;
        }
        if !weights.is_undefined() {
            inner.weights = weights_value(&weights)?;
        }
        if !max_matches.is_undefined() {
            inner.max_matches = i32_value(&max_matches, "maxMatches")?;
        }
        if !symmetrize_conjugated_terminal_groups.is_undefined() {
            inner.symmetrize_conjugated_terminal_groups = bool_value(
                &symmetrize_conjugated_terminal_groups,
                "symmetrizeConjugatedTerminalGroups",
            )?;
        }
        if !ignore_hydrogens.is_undefined() {
            inner.ignore_hydrogens = bool_value(&ignore_hydrogens, "ignoreHydrogens")?;
        }
        if !num_threads.is_undefined() {
            inner.num_threads = i32_value(&num_threads, "numThreads")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(js_name = new)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "AlignmentAtomMap[][] | null")]
        atom_maps: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Float64Array | null")]
        weights: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_matches: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        symmetrize_conjugated_terminal_groups: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] ignore_hydrogens: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] num_threads: JsValue,
    ) -> Result<AllConformerRmsdParameters, JsValue> {
        let mut inner = ck::AllConformerRmsdParameters::default();
        if !atom_maps.is_undefined() {
            inner.atom_maps = atom_maps_value(&atom_maps)?;
        }
        if !weights.is_undefined() {
            inner.weights = weights_value(&weights)?;
        }
        if !max_matches.is_undefined() {
            inner.max_matches = i32_value(&max_matches, "maxMatches")?;
        }
        if !symmetrize_conjugated_terminal_groups.is_undefined() {
            inner.symmetrize_conjugated_terminal_groups = bool_value(
                &symmetrize_conjugated_terminal_groups,
                "symmetrizeConjugatedTerminalGroups",
            )?;
        }
        if !ignore_hydrogens.is_undefined() {
            inner.ignore_hydrogens = bool_value(&ignore_hydrogens, "ignoreHydrogens")?;
        }
        if !num_threads.is_undefined() {
            inner.num_threads = i32_value(&num_threads, "numThreads")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter, js_name = atomMaps, unchecked_return_type = "AlignmentAtomMap[][]")]
    pub fn atom_maps(&self) -> JsValue {
        atom_maps_js(&self.inner.atom_maps)
    }
    #[wasm_bindgen(setter, js_name = atomMaps)]
    pub fn set_atom_maps(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "AlignmentAtomMap[][] | null")] atom_maps: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.atom_maps = atom_maps_value(&atom_maps)?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = weights, unchecked_return_type = "number[] | null")]
    pub fn weights(&self) -> JsValue {
        weights_js(&self.inner.weights)
    }
    #[wasm_bindgen(setter, js_name = weights)]
    pub fn set_weights(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number[] | Float64Array | null")] weights: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.weights = weights_value(&weights)?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = maxMatches)]
    pub fn max_matches(&self) -> i32 {
        self.inner.max_matches
    }
    #[wasm_bindgen(setter, js_name = maxMatches)]
    pub fn set_max_matches(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] max_matches: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.max_matches = i32_value(&max_matches, "maxMatches")?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = symmetrizeConjugatedTerminalGroups)]
    pub fn symmetrize_conjugated_terminal_groups(&self) -> bool {
        self.inner.symmetrize_conjugated_terminal_groups
    }
    #[wasm_bindgen(setter, js_name = symmetrizeConjugatedTerminalGroups)]
    pub fn set_symmetrize_conjugated_terminal_groups(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")]
        symmetrize_conjugated_terminal_groups: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.symmetrize_conjugated_terminal_groups = bool_value(
            &symmetrize_conjugated_terminal_groups,
            "symmetrizeConjugatedTerminalGroups",
        )?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = ignoreHydrogens)]
    pub fn ignore_hydrogens(&self) -> bool {
        self.inner.ignore_hydrogens
    }
    #[wasm_bindgen(setter, js_name = ignoreHydrogens)]
    pub fn set_ignore_hydrogens(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] ignore_hydrogens: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.ignore_hydrogens = bool_value(&ignore_hydrogens, "ignoreHydrogens")?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = numThreads)]
    pub fn num_threads(&self) -> i32 {
        self.inner.num_threads
    }
    #[wasm_bindgen(setter, js_name = numThreads)]
    pub fn set_num_threads(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] num_threads: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.num_threads = i32_value(&num_threads, "numThreads")?;
        Ok(())
    }
}

#[wasm_bindgen]
pub struct ConformerAlignmentParameters {
    pub(crate) inner: ck::ConformerAlignmentParameters,
}

#[wasm_bindgen]
impl ConformerAlignmentParameters {
    #[wasm_bindgen(constructor)]
    pub fn constructor(
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        atom_indices: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        conformer_ids: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Float64Array | null")]
        weights: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] reflect: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_iterations: JsValue,
    ) -> Result<ConformerAlignmentParameters, JsValue> {
        let mut inner = ck::ConformerAlignmentParameters::default();
        if !atom_indices.is_undefined() {
            inner.atom_indices = indices_value(&atom_indices, "atomIndices")?;
        }
        if !conformer_ids.is_undefined() {
            inner.conformer_ids = indices_value(&conformer_ids, "conformerIds")?;
        }
        if !weights.is_undefined() {
            inner.weights = weights_value(&weights)?;
        }
        if !reflect.is_undefined() {
            inner.reflect = bool_value(&reflect, "reflect")?;
        }
        if !max_iterations.is_undefined() {
            inner.max_iterations = u32_value(&max_iterations, "maxIterations")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(js_name = new)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        atom_indices: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        conformer_ids: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Float64Array | null")]
        weights: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] reflect: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_iterations: JsValue,
    ) -> Result<ConformerAlignmentParameters, JsValue> {
        let mut inner = ck::ConformerAlignmentParameters::default();
        if !atom_indices.is_undefined() {
            inner.atom_indices = indices_value(&atom_indices, "atomIndices")?;
        }
        if !conformer_ids.is_undefined() {
            inner.conformer_ids = indices_value(&conformer_ids, "conformerIds")?;
        }
        if !weights.is_undefined() {
            inner.weights = weights_value(&weights)?;
        }
        if !reflect.is_undefined() {
            inner.reflect = bool_value(&reflect, "reflect")?;
        }
        if !max_iterations.is_undefined() {
            inner.max_iterations = u32_value(&max_iterations, "maxIterations")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter, js_name = atomIndices, unchecked_return_type = "number[] | null")]
    pub fn atom_indices(&self) -> JsValue {
        indices_js(&self.inner.atom_indices)
    }
    #[wasm_bindgen(setter, js_name = atomIndices)]
    pub fn set_atom_indices(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number[] | Uint32Array | null")]
        atom_indices: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.atom_indices = indices_value(&atom_indices, "atomIndices")?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = conformerIds, unchecked_return_type = "number[] | null")]
    pub fn conformer_ids(&self) -> JsValue {
        indices_js(&self.inner.conformer_ids)
    }
    #[wasm_bindgen(setter, js_name = conformerIds)]
    pub fn set_conformer_ids(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number[] | Uint32Array | null")]
        conformer_ids: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.conformer_ids = indices_value(&conformer_ids, "conformerIds")?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = weights, unchecked_return_type = "number[] | null")]
    pub fn weights(&self) -> JsValue {
        weights_js(&self.inner.weights)
    }
    #[wasm_bindgen(setter, js_name = weights)]
    pub fn set_weights(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number[] | Float64Array | null")] weights: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.weights = weights_value(&weights)?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = reflect)]
    pub fn reflect(&self) -> bool {
        self.inner.reflect
    }
    #[wasm_bindgen(setter, js_name = reflect)]
    pub fn set_reflect(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] reflect: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.reflect = bool_value(&reflect, "reflect")?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name = maxIterations)]
    pub fn max_iterations(&self) -> u32 {
        self.inner.max_iterations
    }
    #[wasm_bindgen(setter, js_name = maxIterations)]
    pub fn set_max_iterations(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] max_iterations: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.max_iterations = u32_value(&max_iterations, "maxIterations")?;
        Ok(())
    }
}
