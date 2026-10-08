//! Frozen canonical chemistry and depiction parameter transport.
use crate::host_values::*;
use cosmolkit_wasm::rust as ck;
use js_sys::{Array, Map};
use std::collections::BTreeMap;
use wasm_bindgen::prelude::*;
fn default_bool(value: &JsValue, default: bool, name: &str) -> Result<bool, JsValue> {
    if value.is_undefined() {
        Ok(default)
    } else {
        bool_value(value, name)
    }
}
fn default_u32(value: &JsValue, default: u32, name: &str) -> Result<u32, JsValue> {
    if value.is_undefined() {
        Ok(default)
    } else {
        u32_value(value, name)
    }
}
fn default_i32(value: &JsValue, default: i32, name: &str) -> Result<i32, JsValue> {
    if value.is_undefined() {
        Ok(default)
    } else {
        i32_value(value, name)
    }
}
#[wasm_bindgen]
pub struct RemoveHsParams {
    pub(crate) inner: ck::RemoveHsParams,
}
#[wasm_bindgen]
impl RemoveHsParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_degree_zero: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_higher_degrees: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_only_h_neighbors: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_isotopes: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        remove_and_track_isotopes: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_dummy_neighbors: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        remove_defining_bond_stereo: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_with_wedged_bond: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_with_query: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_mapped: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_in_sgroups: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] show_warnings: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_nonimplicit: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] update_explicit_count: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_hydrides: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_nontetrahedral_neighbors: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] sanitize: JsValue,
    ) -> Result<Self, JsValue> {
        let defaults = ck::RemoveHsParams::default();
        Ok(Self {
            inner: ck::RemoveHsParams {
                remove_degree_zero: default_bool(
                    &remove_degree_zero,
                    defaults.remove_degree_zero,
                    "removeDegreeZero",
                )?,
                remove_higher_degrees: default_bool(
                    &remove_higher_degrees,
                    defaults.remove_higher_degrees,
                    "removeHigherDegrees",
                )?,
                remove_only_h_neighbors: default_bool(
                    &remove_only_h_neighbors,
                    defaults.remove_only_h_neighbors,
                    "removeOnlyHNeighbors",
                )?,
                remove_isotopes: default_bool(
                    &remove_isotopes,
                    defaults.remove_isotopes,
                    "removeIsotopes",
                )?,
                remove_and_track_isotopes: default_bool(
                    &remove_and_track_isotopes,
                    defaults.remove_and_track_isotopes,
                    "removeAndTrackIsotopes",
                )?,
                remove_dummy_neighbors: default_bool(
                    &remove_dummy_neighbors,
                    defaults.remove_dummy_neighbors,
                    "removeDummyNeighbors",
                )?,
                remove_defining_bond_stereo: default_bool(
                    &remove_defining_bond_stereo,
                    defaults.remove_defining_bond_stereo,
                    "removeDefiningBondStereo",
                )?,
                remove_with_wedged_bond: default_bool(
                    &remove_with_wedged_bond,
                    defaults.remove_with_wedged_bond,
                    "removeWithWedgedBond",
                )?,
                remove_with_query: default_bool(
                    &remove_with_query,
                    defaults.remove_with_query,
                    "removeWithQuery",
                )?,
                remove_mapped: default_bool(
                    &remove_mapped,
                    defaults.remove_mapped,
                    "removeMapped",
                )?,
                remove_in_sgroups: default_bool(
                    &remove_in_sgroups,
                    defaults.remove_in_sgroups,
                    "removeInSgroups",
                )?,
                show_warnings: default_bool(
                    &show_warnings,
                    defaults.show_warnings,
                    "showWarnings",
                )?,
                remove_nonimplicit: default_bool(
                    &remove_nonimplicit,
                    defaults.remove_nonimplicit,
                    "removeNonimplicit",
                )?,
                update_explicit_count: default_bool(
                    &update_explicit_count,
                    defaults.update_explicit_count,
                    "updateExplicitCount",
                )?,
                remove_hydrides: default_bool(
                    &remove_hydrides,
                    defaults.remove_hydrides,
                    "removeHydrides",
                )?,
                remove_nontetrahedral_neighbors: default_bool(
                    &remove_nontetrahedral_neighbors,
                    defaults.remove_nontetrahedral_neighbors,
                    "removeNontetrahedralNeighbors",
                )?,
                sanitize: default_bool(&sanitize, defaults.sanitize, "sanitize")?,
            },
        })
    }
    #[wasm_bindgen(getter, js_name = removeDegreeZero)]
    pub fn remove_degree_zero(&self) -> bool {
        self.inner.remove_degree_zero
    }
    #[wasm_bindgen(getter, js_name = removeHigherDegrees)]
    pub fn remove_higher_degrees(&self) -> bool {
        self.inner.remove_higher_degrees
    }
    #[wasm_bindgen(getter, js_name = removeOnlyHNeighbors)]
    pub fn remove_only_h_neighbors(&self) -> bool {
        self.inner.remove_only_h_neighbors
    }
    #[wasm_bindgen(getter, js_name = removeIsotopes)]
    pub fn remove_isotopes(&self) -> bool {
        self.inner.remove_isotopes
    }
    #[wasm_bindgen(getter, js_name = removeAndTrackIsotopes)]
    pub fn remove_and_track_isotopes(&self) -> bool {
        self.inner.remove_and_track_isotopes
    }
    #[wasm_bindgen(getter, js_name = removeDummyNeighbors)]
    pub fn remove_dummy_neighbors(&self) -> bool {
        self.inner.remove_dummy_neighbors
    }
    #[wasm_bindgen(getter, js_name = removeDefiningBondStereo)]
    pub fn remove_defining_bond_stereo(&self) -> bool {
        self.inner.remove_defining_bond_stereo
    }
    #[wasm_bindgen(getter, js_name = removeWithWedgedBond)]
    pub fn remove_with_wedged_bond(&self) -> bool {
        self.inner.remove_with_wedged_bond
    }
    #[wasm_bindgen(getter, js_name = removeWithQuery)]
    pub fn remove_with_query(&self) -> bool {
        self.inner.remove_with_query
    }
    #[wasm_bindgen(getter, js_name = removeMapped)]
    pub fn remove_mapped(&self) -> bool {
        self.inner.remove_mapped
    }
    #[wasm_bindgen(getter, js_name = removeInSgroups)]
    pub fn remove_in_sgroups(&self) -> bool {
        self.inner.remove_in_sgroups
    }
    #[wasm_bindgen(getter, js_name = showWarnings)]
    pub fn show_warnings(&self) -> bool {
        self.inner.show_warnings
    }
    #[wasm_bindgen(getter, js_name = removeNonimplicit)]
    pub fn remove_nonimplicit(&self) -> bool {
        self.inner.remove_nonimplicit
    }
    #[wasm_bindgen(getter, js_name = updateExplicitCount)]
    pub fn update_explicit_count(&self) -> bool {
        self.inner.update_explicit_count
    }
    #[wasm_bindgen(getter, js_name = removeHydrides)]
    pub fn remove_hydrides(&self) -> bool {
        self.inner.remove_hydrides
    }
    #[wasm_bindgen(getter, js_name = removeNontetrahedralNeighbors)]
    pub fn remove_nontetrahedral_neighbors(&self) -> bool {
        self.inner.remove_nontetrahedral_neighbors
    }
    #[wasm_bindgen(getter, js_name = sanitize)]
    pub fn sanitize(&self) -> bool {
        self.inner.sanitize
    }
}
#[wasm_bindgen]
pub struct KekulizeParams {
    pub(crate) inner: ck::KekulizeParams,
}
#[wasm_bindgen]
impl KekulizeParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] mark_atoms_bonds: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] canonical: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_backtracks: JsValue,
    ) -> Result<Self, JsValue> {
        let defaults = ck::KekulizeParams::default();
        Ok(Self {
            inner: ck::KekulizeParams {
                mark_atoms_bonds: default_bool(
                    &mark_atoms_bonds,
                    defaults.mark_atoms_bonds,
                    "markAtomsBonds",
                )?,
                canonical: default_bool(&canonical, defaults.canonical, "canonical")?,
                max_backtracks: default_u32(
                    &max_backtracks,
                    defaults.max_backtracks,
                    "maxBacktracks",
                )?,
            },
        })
    }
    #[wasm_bindgen(getter, js_name = markAtomsBonds)]
    pub fn mark_atoms_bonds(&self) -> bool {
        self.inner.mark_atoms_bonds
    }
    #[wasm_bindgen(getter, js_name = canonical)]
    pub fn canonical(&self) -> bool {
        self.inner.canonical
    }
    #[wasm_bindgen(getter, js_name = maxBacktracks)]
    pub fn max_backtracks(&self) -> u32 {
        self.inner.max_backtracks
    }
}
#[wasm_bindgen]
pub struct AddHsParams {
    pub(crate) inner: ck::AddHsParams,
}
#[wasm_bindgen]
impl AddHsParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] explicit_only: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] add_coords: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] add_residue_info: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] skip_queries: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        only_on_atoms: JsValue,
    ) -> Result<Self, JsValue> {
        let defaults = ck::AddHsParams::default();
        Ok(Self {
            inner: ck::AddHsParams {
                explicit_only: default_bool(
                    &explicit_only,
                    defaults.explicit_only,
                    "explicitOnly",
                )?,
                add_coords: default_bool(&add_coords, defaults.add_coords, "addCoords")?,
                add_residue_info: default_bool(
                    &add_residue_info,
                    defaults.add_residue_info,
                    "addResidueInfo",
                )?,
                skip_queries: default_bool(&skip_queries, defaults.skip_queries, "skipQueries")?,
                only_on_atoms: if only_on_atoms.is_null() || only_on_atoms.is_undefined() {
                    None
                } else {
                    Some(
                        sequence(&only_on_atoms, "onlyOnAtoms")?
                            .iter()
                            .map(|v| usize_value(&v, "atom index").map(ck::AtomId::new))
                            .collect::<Result<_, _>>()?,
                    )
                },
            },
        })
    }
    #[wasm_bindgen(getter, js_name = explicitOnly)]
    pub fn explicit_only(&self) -> bool {
        self.inner.explicit_only
    }
    #[wasm_bindgen(getter, js_name = addCoords)]
    pub fn add_coords(&self) -> bool {
        self.inner.add_coords
    }
    #[wasm_bindgen(getter, js_name = addResidueInfo)]
    pub fn add_residue_info(&self) -> bool {
        self.inner.add_residue_info
    }
    #[wasm_bindgen(getter, js_name = skipQueries)]
    pub fn skip_queries(&self) -> bool {
        self.inner.skip_queries
    }
    #[wasm_bindgen(getter, js_name = onlyOnAtoms, unchecked_return_type = "number[] | null")]
    pub fn only_on_atoms(&self) -> JsValue {
        self.inner
            .only_on_atoms
            .as_ref()
            .map_or(JsValue::NULL, |v| {
                v.iter()
                    .map(|id| JsValue::from_f64(id.index() as f64))
                    .collect::<Array>()
                    .into()
            })
    }
}
#[wasm_bindgen]
pub struct Coordinate2DParams {
    pub(crate) inner: ck::Coordinate2DParams,
}
#[wasm_bindgen]
impl Coordinate2DParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(
            unchecked_optional_param_type = "Map<number, [number, number] | Float64Array> | null"
        )]
        coordinate_map: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] canonical_orientation: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] clear_existing_2d: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] flips_per_sample: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] samples: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] sample_seed: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] permute_degree_four: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] force_rdkit: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_ring_templates: JsValue,
    ) -> Result<Self, JsValue> {
        let defaults = ck::Coordinate2DParams::default();
        Ok(Self {
            inner: ck::Coordinate2DParams {
                canonical_orientation: default_bool(
                    &canonical_orientation,
                    defaults.canonical_orientation,
                    "canonicalOrientation",
                )?,
                clear_existing_2d: default_bool(
                    &clear_existing_2d,
                    defaults.clear_existing_2d,
                    "clearExisting2D",
                )?,
                flips_per_sample: default_u32(
                    &flips_per_sample,
                    defaults.flips_per_sample,
                    "flipsPerSample",
                )?,
                samples: default_u32(&samples, defaults.samples, "samples")?,
                sample_seed: default_i32(&sample_seed, defaults.sample_seed, "sampleSeed")?,
                permute_degree_four: default_bool(
                    &permute_degree_four,
                    defaults.permute_degree_four,
                    "permuteDegreeFour",
                )?,
                force_rdkit: default_bool(&force_rdkit, defaults.force_rdkit, "forceRdkit")?,
                use_ring_templates: default_bool(
                    &use_ring_templates,
                    defaults.use_ring_templates,
                    "useRingTemplates",
                )?,
                coordinate_map: coordinates(&coordinate_map)?,
            },
        })
    }
    #[wasm_bindgen(getter, js_name = canonicalOrientation)]
    pub fn canonical_orientation(&self) -> bool {
        self.inner.canonical_orientation
    }
    #[wasm_bindgen(getter, js_name = clearExisting2d)]
    pub fn clear_existing_2d(&self) -> bool {
        self.inner.clear_existing_2d
    }
    #[wasm_bindgen(getter, js_name = flipsPerSample)]
    pub fn flips_per_sample(&self) -> u32 {
        self.inner.flips_per_sample
    }
    #[wasm_bindgen(getter, js_name = samples)]
    pub fn samples(&self) -> u32 {
        self.inner.samples
    }
    #[wasm_bindgen(getter, js_name = sampleSeed)]
    pub fn sample_seed(&self) -> i32 {
        self.inner.sample_seed
    }
    #[wasm_bindgen(getter, js_name = permuteDegreeFour)]
    pub fn permute_degree_four(&self) -> bool {
        self.inner.permute_degree_four
    }
    #[wasm_bindgen(getter, js_name = forceRdkit)]
    pub fn force_rdkit(&self) -> bool {
        self.inner.force_rdkit
    }
    #[wasm_bindgen(getter, js_name = useRingTemplates)]
    pub fn use_ring_templates(&self) -> bool {
        self.inner.use_ring_templates
    }
    #[wasm_bindgen(getter, js_name = coordinateMap, unchecked_return_type = "Map<number, [number, number]>")]
    pub fn coordinate_map(&self) -> Map {
        let result = Map::new();
        for (key, value) in &self.inner.coordinate_map {
            let pair: Array = value.iter().map(|v| JsValue::from_f64(*v)).collect();
            result.set(&JsValue::from_f64(*key as f64), &pair);
        }
        result
    }
}
fn coordinates(value: &JsValue) -> Result<BTreeMap<usize, [f64; 2]>, JsValue> {
    if value.is_undefined() || value.is_null() {
        return Ok(BTreeMap::new());
    }
    let map = value
        .dyn_ref::<Map>()
        .ok_or_else(|| type_error("coordinateMap"))?;
    let mut output = BTreeMap::new();
    for entry in js_sys::try_iter(map)?.ok_or_else(|| type_error("coordinateMap"))? {
        let entry = Array::from(&entry?);
        let key = usize_value(&entry.get(0), "coordinateMap atom index")?;
        let pair = sequence(&entry.get(1), "coordinateMap coordinate")?;
        if pair.length() != 2 {
            return Err(type_error("2D coordinate pair"));
        }
        let x = pair
            .get(0)
            .as_f64()
            .ok_or_else(|| type_error("coordinate x"))?;
        let y = pair
            .get(1)
            .as_f64()
            .ok_or_else(|| type_error("coordinate y"))?;
        output.insert(key, [x, y]);
    }
    Ok(output)
}
#[wasm_bindgen]
pub struct SanitizeOperations {
    inner: ck::SanitizeOperations,
}
#[wasm_bindgen]
impl SanitizeOperations {
    #[wasm_bindgen(js_name = fromBits)]
    pub fn from_bits(
        #[wasm_bindgen(unchecked_param_type = "number")] bits: JsValue,
    ) -> Result<Self, JsValue> {
        ck::SanitizeOperations::from_bits(u32_value(&bits, "sanitize bits")?)
            .map(|inner| Self { inner })
            .map_err(|e| crate::transform_errors::sanitize_error(&e).unwrap_or_else(|e| e))
    }
    pub fn bits(&self) -> u32 {
        self.inner.bits()
    }
    pub fn contains(&self, operation: &Self) -> bool {
        self.inner.contains(operation.inner)
    }
    #[wasm_bindgen(js_name = isEmpty)]
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    pub fn or(&self, other: &Self) -> Self {
        Self {
            inner: self.inner | other.inner,
        }
    }
    pub fn and(&self, other: &Self) -> Self {
        Self {
            inner: self.inner & other.inner,
        }
    }
    #[wasm_bindgen(getter, js_name = NONE)]
    pub fn flag_none() -> Self {
        Self {
            inner: ck::SanitizeOperations::NONE,
        }
    }
    #[wasm_bindgen(getter, js_name = CLEANUP)]
    pub fn flag_cleanup() -> Self {
        Self {
            inner: ck::SanitizeOperations::CLEANUP,
        }
    }
    #[wasm_bindgen(getter, js_name = PROPERTIES)]
    pub fn flag_properties() -> Self {
        Self {
            inner: ck::SanitizeOperations::PROPERTIES,
        }
    }
    #[wasm_bindgen(getter, js_name = SYMM_RINGS)]
    pub fn flag_symm_rings() -> Self {
        Self {
            inner: ck::SanitizeOperations::SYMM_RINGS,
        }
    }
    #[wasm_bindgen(getter, js_name = KEKULIZE)]
    pub fn flag_kekulize() -> Self {
        Self {
            inner: ck::SanitizeOperations::KEKULIZE,
        }
    }
    #[wasm_bindgen(getter, js_name = FIND_RADICALS)]
    pub fn flag_find_radicals() -> Self {
        Self {
            inner: ck::SanitizeOperations::FIND_RADICALS,
        }
    }
    #[wasm_bindgen(getter, js_name = SET_AROMATICITY)]
    pub fn flag_set_aromaticity() -> Self {
        Self {
            inner: ck::SanitizeOperations::SET_AROMATICITY,
        }
    }
    #[wasm_bindgen(getter, js_name = SET_CONJUGATION)]
    pub fn flag_set_conjugation() -> Self {
        Self {
            inner: ck::SanitizeOperations::SET_CONJUGATION,
        }
    }
    #[wasm_bindgen(getter, js_name = SET_HYBRIDIZATION)]
    pub fn flag_set_hybridization() -> Self {
        Self {
            inner: ck::SanitizeOperations::SET_HYBRIDIZATION,
        }
    }
    #[wasm_bindgen(getter, js_name = CLEANUP_CHIRALITY)]
    pub fn flag_cleanup_chirality() -> Self {
        Self {
            inner: ck::SanitizeOperations::CLEANUP_CHIRALITY,
        }
    }
    #[wasm_bindgen(getter, js_name = ADJUST_HS)]
    pub fn flag_adjust_hs() -> Self {
        Self {
            inner: ck::SanitizeOperations::ADJUST_HS,
        }
    }
    #[wasm_bindgen(getter, js_name = CLEANUP_ORGANOMETALLICS)]
    pub fn flag_cleanup_organometallics() -> Self {
        Self {
            inner: ck::SanitizeOperations::CLEANUP_ORGANOMETALLICS,
        }
    }
    #[wasm_bindgen(getter, js_name = CLEANUP_ATROPISOMERS)]
    pub fn flag_cleanup_atropisomers() -> Self {
        Self {
            inner: ck::SanitizeOperations::CLEANUP_ATROPISOMERS,
        }
    }
    #[wasm_bindgen(getter, js_name = ALL)]
    pub fn flag_all() -> Self {
        Self {
            inner: ck::SanitizeOperations::ALL,
        }
    }
}
#[wasm_bindgen]
pub struct SanitizeParams {
    pub(crate) inner: ck::SanitizeParams,
}
#[wasm_bindgen]
impl SanitizeParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "SanitizeOperations | null")]
        operations: JsValue,
    ) -> Result<Self, JsValue> {
        let mut selected = ck::SanitizeOperations::ALL;
        if !operations.is_undefined() && !operations.is_null() {
            visit_sanitize_operations(&operations, &mut |v: &SanitizeOperations| {
                selected = v.inner
            })?;
        }
        Ok(Self {
            inner: ck::SanitizeParams {
                operations: selected,
            },
        })
    }
    #[wasm_bindgen(getter)]
    pub fn operations(&self) -> SanitizeOperations {
        SanitizeOperations {
            inner: self.inner.operations,
        }
    }
}
// This callback only copies a flag value. Any callback exception is a host
// class borrowing failure; preserve it as the cause of the argument TypeError.
#[wasm_bindgen(
    inline_js = "export function visit_sanitize_operations(value, visit) { try { visit(value); } catch (cause) { throw new TypeError('invalid SanitizeOperations', { cause }); } }"
)]
extern "C" {
    #[wasm_bindgen(catch)]
    fn visit_sanitize_operations(
        value: &JsValue,
        visit: &mut dyn FnMut(&SanitizeOperations),
    ) -> Result<(), JsValue>;
}
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum SanitizeStage {
    None = 0x000,
    Cleanup = 0x001,
    Properties = 0x002,
    SymmRings = 0x004,
    Kekulize = 0x008,
    FindRadicals = 0x010,
    SetAromaticity = 0x020,
    SetConjugation = 0x040,
    SetHybridization = 0x080,
    CleanupChirality = 0x100,
    AdjustHs = 0x200,
    CleanupOrganometallics = 0x400,
    CleanupAtropisomers = 0x800,
}
pub(crate) fn stage(value: ck::SanitizeStage) -> SanitizeStage {
    match value {
        ck::SanitizeStage::None => SanitizeStage::None,
        ck::SanitizeStage::Cleanup => SanitizeStage::Cleanup,
        ck::SanitizeStage::Properties => SanitizeStage::Properties,
        ck::SanitizeStage::SymmRings => SanitizeStage::SymmRings,
        ck::SanitizeStage::Kekulize => SanitizeStage::Kekulize,
        ck::SanitizeStage::FindRadicals => SanitizeStage::FindRadicals,
        ck::SanitizeStage::SetAromaticity => SanitizeStage::SetAromaticity,
        ck::SanitizeStage::SetConjugation => SanitizeStage::SetConjugation,
        ck::SanitizeStage::SetHybridization => SanitizeStage::SetHybridization,
        ck::SanitizeStage::CleanupChirality => SanitizeStage::CleanupChirality,
        ck::SanitizeStage::AdjustHs => SanitizeStage::AdjustHs,
        ck::SanitizeStage::CleanupOrganometallics => SanitizeStage::CleanupOrganometallics,
        ck::SanitizeStage::CleanupAtropisomers => SanitizeStage::CleanupAtropisomers,
    }
}
