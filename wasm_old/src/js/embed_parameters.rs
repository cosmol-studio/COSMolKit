//! Full frozen embedding configuration; Python-equivalent map and scalar conversion.
use crate::host_values::{bool_value, i32_value, point, sequence, type_error, u32_value};
use cosmolkit_wasm::rust as ck;
use std::collections::BTreeMap;
use wasm_bindgen::prelude::*;
fn number(v: &JsValue, name: &str) -> Result<f64, JsValue> {
    v.as_f64().ok_or_else(|| type_error(name))
}
fn coord_map(v: &JsValue) -> Result<Option<BTreeMap<i32, [f64; 3]>>, JsValue> {
    if v.is_null() || v.is_undefined() {
        return Ok(None);
    }
    if !v.is_instance_of::<js_sys::Map>() {
        return Err(type_error("coordMap"));
    }
    let mut out = BTreeMap::new();
    for entry in js_sys::Array::from(v).iter() {
        let pair = js_sys::Array::from(&entry);
        out.insert(
            i32_value(&pair.get(0), "coordMap key")?,
            point(&pair.get(1))?,
        );
    }
    Ok(Some(out))
}
fn cpci(v: &JsValue) -> Result<Option<BTreeMap<(u32, u32), f64>>, JsValue> {
    if v.is_null() || v.is_undefined() {
        return Ok(None);
    }
    if !v.is_instance_of::<js_sys::Map>() {
        return Err(type_error("cpci"));
    }
    let mut out = BTreeMap::new();
    for entry in js_sys::Array::from(v).iter() {
        let pair = js_sys::Array::from(&entry);
        let key = sequence(&pair.get(0), "cpci key")?;
        if key.length() != 2 {
            return Err(js_sys::RangeError::new("cpci key must contain two atom indices").into());
        }
        out.insert(
            (
                u32_value(&key.get(0), "cpci first atom")?,
                u32_value(&key.get(1), "cpci second atom")?,
            ),
            number(&pair.get(1), "cpci value")?,
        );
    }
    Ok(Some(out))
}
#[derive(Clone)]
#[wasm_bindgen]
pub struct EmbedParams {
    pub(crate) inner: ck::EmbedParams,
}
#[cosmolkit_wasm::javascript_options("EmbedParamsOptions", clone_snapshot)]
#[wasm_bindgen]
impl EmbedParams {
    #[wasm_bindgen(constructor)]
    #[allow(clippy::too_many_arguments)]
    pub fn construct(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_iterations: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] num_threads: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] random_seed: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] clear_conformers: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_random_coords: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] box_size_mult: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] rand_neg_eig: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] num_zero_fail: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "Map<number, number[] | Float64Array> | null"
        )]
        coord_map: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] optimizer_force_tol: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        ignore_smoothing_failures: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] enforce_chirality: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        use_exp_torsion_angle_prefs: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_basic_knowledge: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] verbose: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] basin_thresh: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] prune_rms_thresh: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        only_heavy_atoms_for_rms: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] et_version: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        embed_fragments_separately: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_small_ring_torsions: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_macrocycle_torsions: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_macrocycle14config: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] timeout: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "Map<readonly [number, number], number> | null"
        )]
        cpci: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] force_trans_amides: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        use_symmetry_for_pruning: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] bounds_mat_force_scaling: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] track_failures: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        enable_sequential_random_seeds: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        symmetrize_conjugated_terminal_groups_for_pruning: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::EmbedParams::new();
        // COSMolKit❗✔️: inner.max_iterations = max_iterations;
        if !max_iterations.is_undefined() {
            inner.max_iterations = u32_value(&max_iterations, "maxIterations")?;
        }
        // COSMolKit❗✔️: inner.num_threads = num_threads;
        if !num_threads.is_undefined() {
            inner.num_threads = i32_value(&num_threads, "numThreads")?;
        }
        // COSMolKit❗✔️: inner.random_seed = random_seed;
        if !random_seed.is_undefined() {
            inner.random_seed = i32_value(&random_seed, "randomSeed")?;
        }
        // COSMolKit❗✔️: inner.clear_confs = clear_confs;
        if !clear_conformers.is_undefined() {
            inner.clear_conformers = bool_value(&clear_conformers, "clearConformers")?;
        }
        // COSMolKit❗✔️: inner.use_random_coords = use_random_coords;
        if !use_random_coords.is_undefined() {
            inner.use_random_coords = bool_value(&use_random_coords, "useRandomCoords")?;
        }
        // COSMolKit❗✔️: inner.box_size_mult = box_size_mult;
        if !box_size_mult.is_undefined() {
            inner.box_size_mult = number(&box_size_mult, "boxSizeMult")?;
        }
        // COSMolKit❗✔️: inner.rand_neg_eig = rand_neg_eig;
        if !rand_neg_eig.is_undefined() {
            inner.rand_neg_eig = bool_value(&rand_neg_eig, "randNegEig")?;
        }
        // COSMolKit❗✔️: inner.num_zero_fail = num_zero_fail;
        if !num_zero_fail.is_undefined() {
            inner.num_zero_fail = u32_value(&num_zero_fail, "numZeroFail")?;
        }
        // COSMolKit❗✔️: inner.coord_map = coord_map;
        inner.coord_map = self::coord_map(&coord_map)?;
        // COSMolKit❗✔️: inner.optimizer_force_tol = optimizer_force_tol;
        if !optimizer_force_tol.is_undefined() {
            inner.optimizer_force_tol = number(&optimizer_force_tol, "optimizerForceTol")?;
        }
        // COSMolKit❗✔️: inner.ignore_smoothing_failures = ignore_smoothing_failures;
        if !ignore_smoothing_failures.is_undefined() {
            inner.ignore_smoothing_failures =
                bool_value(&ignore_smoothing_failures, "ignoreSmoothingFailures")?;
        }
        // COSMolKit❗✔️: inner.enforce_chirality = enforce_chirality;
        if !enforce_chirality.is_undefined() {
            inner.enforce_chirality = bool_value(&enforce_chirality, "enforceChirality")?;
        }
        // COSMolKit❗✔️: inner.use_exp_torsion_angle_prefs = use_exp_torsion_angle_prefs;
        if !use_exp_torsion_angle_prefs.is_undefined() {
            inner.use_exp_torsion_angle_prefs =
                bool_value(&use_exp_torsion_angle_prefs, "useExpTorsionAnglePrefs")?;
        }
        // COSMolKit❗✔️: inner.use_basic_knowledge = use_basic_knowledge;
        if !use_basic_knowledge.is_undefined() {
            inner.use_basic_knowledge = bool_value(&use_basic_knowledge, "useBasicKnowledge")?;
        }
        // COSMolKit❗✔️: inner.verbose = verbose;
        if !verbose.is_undefined() {
            inner.verbose = bool_value(&verbose, "verbose")?;
        }
        // COSMolKit❗✔️: inner.basin_thresh = basin_thresh;
        if !basin_thresh.is_undefined() {
            inner.basin_thresh = number(&basin_thresh, "basinThresh")?;
        }
        // COSMolKit❗✔️: inner.prune_rms_thresh = prune_rms_thresh;
        if !prune_rms_thresh.is_undefined() {
            inner.prune_rms_thresh = number(&prune_rms_thresh, "pruneRmsThresh")?;
        }
        // COSMolKit❗✔️: inner.only_heavy_atoms_for_rms = only_heavy_atoms_for_rms;
        if !only_heavy_atoms_for_rms.is_undefined() {
            inner.only_heavy_atoms_for_rms =
                bool_value(&only_heavy_atoms_for_rms, "onlyHeavyAtomsForRms")?;
        }
        // COSMolKit❗✔️: inner.et_version = et_version;
        if !et_version.is_undefined() {
            inner.et_version = u32_value(&et_version, "etVersion")?;
        }
        // COSMolKit❗✔️: inner.embed_fragments_separately = embed_fragments_separately;
        if !embed_fragments_separately.is_undefined() {
            inner.embed_fragments_separately =
                bool_value(&embed_fragments_separately, "embedFragmentsSeparately")?;
        }
        // COSMolKit❗✔️: inner.use_small_ring_torsions = use_small_ring_torsions;
        if !use_small_ring_torsions.is_undefined() {
            inner.use_small_ring_torsions =
                bool_value(&use_small_ring_torsions, "useSmallRingTorsions")?;
        }
        // COSMolKit❗✔️: inner.use_macrocycle_torsions = use_macrocycle_torsions;
        if !use_macrocycle_torsions.is_undefined() {
            inner.use_macrocycle_torsions =
                bool_value(&use_macrocycle_torsions, "useMacrocycleTorsions")?;
        }
        // COSMolKit❗✔️: inner.use_macrocycle14config = use_macrocycle14config;
        if !use_macrocycle14config.is_undefined() {
            inner.use_macrocycle14config =
                bool_value(&use_macrocycle14config, "useMacrocycle14Config")?;
        }
        // COSMolKit❗✔️: inner.timeout = timeout;
        if !timeout.is_undefined() {
            inner.timeout = u32_value(&timeout, "timeout")?;
        }
        // COSMolKit❗✔️: inner.cpci = cpci;
        inner.cpci = self::cpci(&cpci)?;
        // COSMolKit❗✔️: inner.force_trans_amides = force_trans_amides;
        if !force_trans_amides.is_undefined() {
            inner.force_trans_amides = bool_value(&force_trans_amides, "forceTransAmides")?;
        }
        // COSMolKit❗✔️: inner.use_symmetry_for_pruning = use_symmetry_for_pruning;
        if !use_symmetry_for_pruning.is_undefined() {
            inner.use_symmetry_for_pruning =
                bool_value(&use_symmetry_for_pruning, "useSymmetryForPruning")?;
        }
        // COSMolKit❗✔️: inner.bounds_mat_force_scaling = bounds_mat_force_scaling;
        if !bounds_mat_force_scaling.is_undefined() {
            inner.bounds_mat_force_scaling =
                number(&bounds_mat_force_scaling, "boundsMatForceScaling")?;
        }
        // COSMolKit❗✔️: inner.track_failures = track_failures;
        if !track_failures.is_undefined() {
            inner.track_failures = bool_value(&track_failures, "trackFailures")?;
        }
        // COSMolKit❗✔️: inner.enable_sequential_random_seeds = enable_sequential_random_seeds;
        if !enable_sequential_random_seeds.is_undefined() {
            inner.enable_sequential_random_seeds = bool_value(
                &enable_sequential_random_seeds,
                "enableSequentialRandomSeeds",
            )?;
        }
        // COSMolKit❗✔️: inner.symmetrize_conjugated_terminal_groups_for_pruning = symmetrize_conjugated_terminal_groups_for_pruning;
        if !symmetrize_conjugated_terminal_groups_for_pruning.is_undefined() {
            inner.symmetrize_conjugated_terminal_groups_for_pruning = bool_value(
                &symmetrize_conjugated_terminal_groups_for_pruning,
                "symmetrizeConjugatedTerminalGroupsForPruning",
            )?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(js_name=new)]
    pub fn factory_new() -> Self {
        // COSMolKit❗✔️: inner: ck::EmbedParams::new(),
        Self {
            inner: ck::EmbedParams::new(),
        }
    }
    #[wasm_bindgen(js_name=dg)]
    pub fn factory_dg() -> Self {
        // COSMolKit❗✔️: inner: ck::EmbedParams::dg(),
        Self {
            inner: ck::EmbedParams::dg(),
        }
    }
    #[wasm_bindgen(js_name=kdg)]
    pub fn factory_kdg() -> Self {
        // COSMolKit❗✔️: inner: ck::EmbedParams::kdg(),
        Self {
            inner: ck::EmbedParams::kdg(),
        }
    }
    #[wasm_bindgen(js_name=etdg)]
    pub fn factory_etdg() -> Self {
        // COSMolKit❗✔️: inner: ck::EmbedParams::etdg(),
        Self {
            inner: ck::EmbedParams::etdg(),
        }
    }
    #[wasm_bindgen(js_name=etdgV2)]
    pub fn factory_etdg_v2() -> Self {
        // COSMolKit❗✔️: inner: ck::EmbedParams::etdg_v2(),
        Self {
            inner: ck::EmbedParams::etdg_v2(),
        }
    }
    #[wasm_bindgen(js_name=etkdg)]
    pub fn factory_etkdg() -> Self {
        // COSMolKit❗✔️: inner: ck::EmbedParams::etkdg(),
        Self {
            inner: ck::EmbedParams::etkdg(),
        }
    }
    #[wasm_bindgen(js_name=etkdgV2)]
    pub fn factory_etkdg_v2() -> Self {
        // COSMolKit❗✔️: inner: ck::EmbedParams::etkdg_v2(),
        Self {
            inner: ck::EmbedParams::etkdg_v2(),
        }
    }
    #[wasm_bindgen(js_name=etkdgV3)]
    pub fn factory_etkdg_v3() -> Self {
        // COSMolKit❗✔️: inner: ck::EmbedParams::etkdg_v3(),
        Self {
            inner: ck::EmbedParams::etkdg_v3(),
        }
    }
    #[wasm_bindgen(js_name=srEtkdgV3)]
    pub fn factory_sr_etkdg_v3() -> Self {
        // COSMolKit❗✔️: inner: ck::EmbedParams::sr_etkdg_v3(),
        Self {
            inner: ck::EmbedParams::sr_etkdg_v3(),
        }
    }
    #[wasm_bindgen(getter,js_name=maxIterations)]
    pub fn max_iterations(&self) -> u32 {
        self.inner.max_iterations()
    }
    #[wasm_bindgen(getter,js_name=numThreads)]
    pub fn num_threads(&self) -> i32 {
        self.inner.num_threads()
    }
    #[wasm_bindgen(getter,js_name=randomSeed)]
    pub fn random_seed(&self) -> i32 {
        self.inner.random_seed()
    }
    #[wasm_bindgen(getter,js_name=clearConformers)]
    pub fn clear_conformers(&self) -> bool {
        self.inner.clear_conformers()
    }
    #[wasm_bindgen(getter,js_name=useRandomCoords)]
    pub fn use_random_coords(&self) -> bool {
        self.inner.use_random_coords()
    }
    #[wasm_bindgen(getter,js_name=boxSizeMult)]
    pub fn box_size_mult(&self) -> f64 {
        self.inner.box_size_mult()
    }
    #[wasm_bindgen(getter,js_name=randNegEig)]
    pub fn rand_neg_eig(&self) -> bool {
        self.inner.rand_neg_eig()
    }
    #[wasm_bindgen(getter,js_name=numZeroFail)]
    pub fn num_zero_fail(&self) -> u32 {
        self.inner.num_zero_fail()
    }
    #[wasm_bindgen(getter,js_name=coordMap,unchecked_return_type="Map<number, number[]> | null")]
    pub fn coord_map(&self) -> JsValue {
        let Some(value) = self.inner.coord_map() else {
            return JsValue::NULL;
        };
        let map = js_sys::Map::new();
        for (key, point) in value {
            let point: js_sys::Array = point.into_iter().map(|v| JsValue::from(*v)).collect();
            map.set(&(*key).into(), &point.into());
        }
        map.into()
    }
    #[wasm_bindgen(getter,js_name=optimizerForceTol)]
    pub fn optimizer_force_tol(&self) -> f64 {
        self.inner.optimizer_force_tol()
    }
    #[wasm_bindgen(getter,js_name=ignoreSmoothingFailures)]
    pub fn ignore_smoothing_failures(&self) -> bool {
        self.inner.ignore_smoothing_failures()
    }
    #[wasm_bindgen(getter,js_name=enforceChirality)]
    pub fn enforce_chirality(&self) -> bool {
        self.inner.enforce_chirality()
    }
    #[wasm_bindgen(getter,js_name=useExpTorsionAnglePrefs)]
    pub fn use_exp_torsion_angle_prefs(&self) -> bool {
        self.inner.use_exp_torsion_angle_prefs()
    }
    #[wasm_bindgen(getter,js_name=useBasicKnowledge)]
    pub fn use_basic_knowledge(&self) -> bool {
        self.inner.use_basic_knowledge()
    }
    #[wasm_bindgen(getter,js_name=verbose)]
    pub fn verbose(&self) -> bool {
        self.inner.verbose()
    }
    #[wasm_bindgen(getter,js_name=basinThresh)]
    pub fn basin_thresh(&self) -> f64 {
        self.inner.basin_thresh()
    }
    #[wasm_bindgen(getter,js_name=pruneRmsThresh)]
    pub fn prune_rms_thresh(&self) -> f64 {
        self.inner.prune_rms_thresh()
    }
    #[wasm_bindgen(getter,js_name=onlyHeavyAtomsForRms)]
    pub fn only_heavy_atoms_for_rms(&self) -> bool {
        self.inner.only_heavy_atoms_for_rms()
    }
    #[wasm_bindgen(getter,js_name=etVersion)]
    pub fn et_version(&self) -> u32 {
        self.inner.et_version()
    }
    #[wasm_bindgen(getter,js_name=embedFragmentsSeparately)]
    pub fn embed_fragments_separately(&self) -> bool {
        self.inner.embed_fragments_separately()
    }
    #[wasm_bindgen(getter,js_name=useSmallRingTorsions)]
    pub fn use_small_ring_torsions(&self) -> bool {
        self.inner.use_small_ring_torsions()
    }
    #[wasm_bindgen(getter,js_name=useMacrocycleTorsions)]
    pub fn use_macrocycle_torsions(&self) -> bool {
        self.inner.use_macrocycle_torsions()
    }
    #[wasm_bindgen(getter,js_name=useMacrocycle14config)]
    pub fn use_macrocycle14config(&self) -> bool {
        self.inner.use_macrocycle14config()
    }
    #[wasm_bindgen(getter,js_name=timeout)]
    pub fn timeout(&self) -> u32 {
        self.inner.timeout()
    }
    #[wasm_bindgen(getter,js_name=cpci,unchecked_return_type="Map<readonly [number, number], number> | null")]
    pub fn cpci(&self) -> JsValue {
        let Some(value) = self.inner.cpci() else {
            return JsValue::NULL;
        };
        let map = js_sys::Map::new();
        for ((a, b), value) in value {
            let key = js_sys::Array::of2(&(*a).into(), &(*b).into());
            map.set(&key.into(), &(*value).into());
        }
        map.into()
    }
    #[wasm_bindgen(getter,js_name=forceTransAmides)]
    pub fn force_trans_amides(&self) -> bool {
        self.inner.force_trans_amides()
    }
    #[wasm_bindgen(getter,js_name=useSymmetryForPruning)]
    pub fn use_symmetry_for_pruning(&self) -> bool {
        self.inner.use_symmetry_for_pruning()
    }
    #[wasm_bindgen(getter,js_name=boundsMatForceScaling)]
    pub fn bounds_mat_force_scaling(&self) -> f64 {
        self.inner.bounds_mat_force_scaling()
    }
    #[wasm_bindgen(getter,js_name=trackFailures)]
    pub fn track_failures(&self) -> bool {
        self.inner.track_failures()
    }
    #[wasm_bindgen(getter,js_name=enableSequentialRandomSeeds)]
    pub fn enable_sequential_random_seeds(&self) -> bool {
        self.inner.enable_sequential_random_seeds()
    }
    #[wasm_bindgen(getter,js_name=symmetrizeConjugatedTerminalGroupsForPruning)]
    pub fn symmetrize_conjugated_terminal_groups_for_pruning(&self) -> bool {
        self.inner
            .symmetrize_conjugated_terminal_groups_for_pruning()
    }
    #[wasm_bindgen(setter, js_name = maxIterations)]
    pub fn set_max_iterations(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.max_iterations = u32_value(&value, "maxIterations")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = numThreads)]
    pub fn set_num_threads(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.num_threads = i32_value(&value, "numThreads")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = randomSeed)]
    pub fn set_random_seed(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.random_seed = i32_value(&value, "randomSeed")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = clearConformers)]
    pub fn set_clear_conformers(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.clear_conformers = bool_value(&value, "clearConformers")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = useRandomCoords)]
    pub fn set_use_random_coords(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.use_random_coords = bool_value(&value, "useRandomCoords")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = boxSizeMult)]
    pub fn set_box_size_mult(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.box_size_mult = number(&value, "boxSizeMult")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = randNegEig)]
    pub fn set_rand_neg_eig(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.rand_neg_eig = bool_value(&value, "randNegEig")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = numZeroFail)]
    pub fn set_num_zero_fail(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.num_zero_fail = u32_value(&value, "numZeroFail")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = coordMap)]
    pub fn set_coord_map(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "Map<number, number[] | Float64Array> | null")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.coord_map = self::coord_map(&value)?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = optimizerForceTol)]
    pub fn set_optimizer_force_tol(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.optimizer_force_tol = number(&value, "optimizerForceTol")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = ignoreSmoothingFailures)]
    pub fn set_ignore_smoothing_failures(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.ignore_smoothing_failures = bool_value(&value, "ignoreSmoothingFailures")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = enforceChirality)]
    pub fn set_enforce_chirality(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.enforce_chirality = bool_value(&value, "enforceChirality")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = useExpTorsionAnglePrefs)]
    pub fn set_use_exp_torsion_angle_prefs(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.use_exp_torsion_angle_prefs = bool_value(&value, "useExpTorsionAnglePrefs")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = useBasicKnowledge)]
    pub fn set_use_basic_knowledge(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.use_basic_knowledge = bool_value(&value, "useBasicKnowledge")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = verbose)]
    pub fn set_verbose(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.verbose = bool_value(&value, "verbose")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = basinThresh)]
    pub fn set_basin_thresh(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.basin_thresh = number(&value, "basinThresh")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = pruneRmsThresh)]
    pub fn set_prune_rms_thresh(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.prune_rms_thresh = number(&value, "pruneRmsThresh")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = onlyHeavyAtomsForRms)]
    pub fn set_only_heavy_atoms_for_rms(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.only_heavy_atoms_for_rms = bool_value(&value, "onlyHeavyAtomsForRms")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = etVersion)]
    pub fn set_et_version(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.et_version = u32_value(&value, "etVersion")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = embedFragmentsSeparately)]
    pub fn set_embed_fragments_separately(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.embed_fragments_separately = bool_value(&value, "embedFragmentsSeparately")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = useSmallRingTorsions)]
    pub fn set_use_small_ring_torsions(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.use_small_ring_torsions = bool_value(&value, "useSmallRingTorsions")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = useMacrocycleTorsions)]
    pub fn set_use_macrocycle_torsions(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.use_macrocycle_torsions = bool_value(&value, "useMacrocycleTorsions")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = useMacrocycle14config)]
    pub fn set_use_macrocycle14config(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.use_macrocycle14config = bool_value(&value, "useMacrocycle14config")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = timeout)]
    pub fn set_timeout(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.timeout = u32_value(&value, "timeout")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = cpci)]
    pub fn set_cpci(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "Map<readonly [number, number], number> | null")]
        value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.cpci = self::cpci(&value)?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = forceTransAmides)]
    pub fn set_force_trans_amides(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.force_trans_amides = bool_value(&value, "forceTransAmides")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = useSymmetryForPruning)]
    pub fn set_use_symmetry_for_pruning(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.use_symmetry_for_pruning = bool_value(&value, "useSymmetryForPruning")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = boundsMatForceScaling)]
    pub fn set_bounds_mat_force_scaling(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.bounds_mat_force_scaling = number(&value, "boundsMatForceScaling")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = trackFailures)]
    pub fn set_track_failures(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.track_failures = bool_value(&value, "trackFailures")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = enableSequentialRandomSeeds)]
    pub fn set_enable_sequential_random_seeds(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.enable_sequential_random_seeds =
            bool_value(&value, "enableSequentialRandomSeeds")?;
        Ok(())
    }
    #[wasm_bindgen(setter, js_name = symmetrizeConjugatedTerminalGroupsForPruning)]
    pub fn set_symmetrize_conjugated_terminal_groups_for_pruning(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.symmetrize_conjugated_terminal_groups_for_pruning =
            bool_value(&value, "symmetrizeConjugatedTerminalGroupsForPruning")?;
        Ok(())
    }
    pub fn failures(&self) -> Vec<u32> {
        self.inner.failures().clone()
    }
    #[wasm_bindgen(js_name=toJson)]
    pub fn to_json(&self) -> String {
        self.inner.to_json()
    }
    #[wasm_bindgen(js_name=withJson)]
    pub fn with_json(&self, json: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: self.inner.with_json(json)
        self.inner
            .with_json(json)
            .map(|inner| Self { inner })
            .map_err(|e| crate::conformer_errors::params_error(&e).unwrap_or_else(|e| e))
    }
}
