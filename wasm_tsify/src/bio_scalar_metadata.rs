//! Native readonly metadata records; defaults and state remain native.
use crate::input;
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct BioReflectionsInfo(pub(crate) ck::BioReflectionsInfo);
crate::fields::value_field!(
    ck::BioReflectionsInfo,
    BioReflectionsInfo,
    BioReflectionsInfoInput,
    "BioReflectionsInfo"
);
crate::fields::value_array!(
    ck::BioReflectionsInfo,
    BioReflectionsInfoArray,
    "ReadonlyArray<BioReflectionsInfo>"
);
#[wasm_bindgen]
impl BioReflectionsInfo {
    #[wasm_bindgen(constructor)]
    pub fn new(
        resolution_high: Option<Ts<input::F64Input>>,
        resolution_low: Option<Ts<input::F64Input>>,
        completeness: Option<Ts<input::F64Input>>,
        redundancy: Option<Ts<input::F64Input>>,
        r_merge: Option<Ts<input::F64Input>>,
        r_sym: Option<Ts<input::F64Input>>,
        mean_i_over_sigma: Option<Ts<input::F64Input>>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioReflectionsInfo::default();
        Ok(Self(ck::BioReflectionsInfo {
            resolution_high: resolution_high
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.resolution_high),
            resolution_low: resolution_low
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.resolution_low),
            completeness: completeness
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.completeness),
            redundancy: redundancy
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.redundancy),
            r_merge: r_merge
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.r_merge),
            r_sym: r_sym
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.r_sym),
            mean_i_over_sigma: mean_i_over_sigma
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.mean_i_over_sigma),
        }))
    }
    #[wasm_bindgen(getter,js_name=resolutionHigh)]
    pub fn resolution_high(&self) -> f64 {
        self.0.resolution_high
    }
    #[wasm_bindgen(getter,js_name=resolutionLow)]
    pub fn resolution_low(&self) -> f64 {
        self.0.resolution_low
    }
    #[wasm_bindgen(getter,js_name=completeness)]
    pub fn completeness(&self) -> f64 {
        self.0.completeness
    }
    #[wasm_bindgen(getter,js_name=redundancy)]
    pub fn redundancy(&self) -> f64 {
        self.0.redundancy
    }
    #[wasm_bindgen(getter,js_name=rMerge)]
    pub fn r_merge(&self) -> f64 {
        self.0.r_merge
    }
    #[wasm_bindgen(getter,js_name=rSym)]
    pub fn r_sym(&self) -> f64 {
        self.0.r_sym
    }
    #[wasm_bindgen(getter,js_name=meanIOverSigma)]
    pub fn mean_i_over_sigma(&self) -> f64 {
        self.0.mean_i_over_sigma
    }
}
#[wasm_bindgen]
pub struct BioBasicRefinementInfo(pub(crate) ck::BioBasicRefinementInfo);
crate::fields::value_field!(
    ck::BioBasicRefinementInfo,
    BioBasicRefinementInfo,
    BioBasicRefinementInfoInput,
    "BioBasicRefinementInfo"
);
crate::fields::value_array!(
    ck::BioBasicRefinementInfo,
    BioBasicRefinementInfoArray,
    "ReadonlyArray<BioBasicRefinementInfo>"
);
#[wasm_bindgen]
impl BioBasicRefinementInfo {
    #[wasm_bindgen(constructor)]
    pub fn new(
        resolution_high: Option<Ts<input::F64Input>>,
        resolution_low: Option<Ts<input::F64Input>>,
        completeness: Option<Ts<input::F64Input>>,
        reflection_count: Option<Ts<input::I32Input>>,
        work_set_count: Option<Ts<input::I32Input>>,
        rfree_set_count: Option<Ts<input::I32Input>>,
        r_all: Option<Ts<input::F64Input>>,
        r_work: Option<Ts<input::F64Input>>,
        r_free: Option<Ts<input::F64Input>>,
        cc_fo_fc_work: Option<Ts<input::F64Input>>,
        cc_fo_fc_free: Option<Ts<input::F64Input>>,
        fsc_work: Option<Ts<input::F64Input>>,
        fsc_free: Option<Ts<input::F64Input>>,
        cc_intensity_work: Option<Ts<input::F64Input>>,
        cc_intensity_free: Option<Ts<input::F64Input>>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioBasicRefinementInfo::default();
        Ok(Self(ck::BioBasicRefinementInfo {
            resolution_high: resolution_high
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.resolution_high),
            resolution_low: resolution_low
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.resolution_low),
            completeness: completeness
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.completeness),
            reflection_count: reflection_count
                .map(input::read::<i32>)
                .transpose()?
                .unwrap_or(d.reflection_count),
            work_set_count: work_set_count
                .map(input::read::<i32>)
                .transpose()?
                .unwrap_or(d.work_set_count),
            rfree_set_count: rfree_set_count
                .map(input::read::<i32>)
                .transpose()?
                .unwrap_or(d.rfree_set_count),
            r_all: r_all
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.r_all),
            r_work: r_work
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.r_work),
            r_free: r_free
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.r_free),
            cc_fo_fc_work: cc_fo_fc_work
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.cc_fo_fc_work),
            cc_fo_fc_free: cc_fo_fc_free
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.cc_fo_fc_free),
            fsc_work: fsc_work
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.fsc_work),
            fsc_free: fsc_free
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.fsc_free),
            cc_intensity_work: cc_intensity_work
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.cc_intensity_work),
            cc_intensity_free: cc_intensity_free
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.cc_intensity_free),
        }))
    }
    #[wasm_bindgen(getter,js_name=resolutionHigh)]
    pub fn resolution_high(&self) -> f64 {
        self.0.resolution_high
    }
    #[wasm_bindgen(getter,js_name=resolutionLow)]
    pub fn resolution_low(&self) -> f64 {
        self.0.resolution_low
    }
    #[wasm_bindgen(getter,js_name=completeness)]
    pub fn completeness(&self) -> f64 {
        self.0.completeness
    }
    #[wasm_bindgen(getter,js_name=reflectionCount)]
    pub fn reflection_count(&self) -> i32 {
        self.0.reflection_count
    }
    #[wasm_bindgen(getter,js_name=workSetCount)]
    pub fn work_set_count(&self) -> i32 {
        self.0.work_set_count
    }
    #[wasm_bindgen(getter,js_name=rfreeSetCount)]
    pub fn rfree_set_count(&self) -> i32 {
        self.0.rfree_set_count
    }
    #[wasm_bindgen(getter,js_name=rAll)]
    pub fn r_all(&self) -> f64 {
        self.0.r_all
    }
    #[wasm_bindgen(getter,js_name=rWork)]
    pub fn r_work(&self) -> f64 {
        self.0.r_work
    }
    #[wasm_bindgen(getter,js_name=rFree)]
    pub fn r_free(&self) -> f64 {
        self.0.r_free
    }
    #[wasm_bindgen(getter,js_name=ccFoFcWork)]
    pub fn cc_fo_fc_work(&self) -> f64 {
        self.0.cc_fo_fc_work
    }
    #[wasm_bindgen(getter,js_name=ccFoFcFree)]
    pub fn cc_fo_fc_free(&self) -> f64 {
        self.0.cc_fo_fc_free
    }
    #[wasm_bindgen(getter,js_name=fscWork)]
    pub fn fsc_work(&self) -> f64 {
        self.0.fsc_work
    }
    #[wasm_bindgen(getter,js_name=fscFree)]
    pub fn fsc_free(&self) -> f64 {
        self.0.fsc_free
    }
    #[wasm_bindgen(getter,js_name=ccIntensityWork)]
    pub fn cc_intensity_work(&self) -> f64 {
        self.0.cc_intensity_work
    }
    #[wasm_bindgen(getter,js_name=ccIntensityFree)]
    pub fn cc_intensity_free(&self) -> f64 {
        self.0.cc_intensity_free
    }
}
#[wasm_bindgen]
pub struct BioRefinementRestraint(pub(crate) ck::BioRefinementRestraint);
crate::fields::value_field!(
    ck::BioRefinementRestraint,
    BioRefinementRestraint,
    BioRefinementRestraintInput,
    "BioRefinementRestraint"
);
crate::fields::value_array!(
    ck::BioRefinementRestraint,
    BioRefinementRestraintArray,
    "ReadonlyArray<BioRefinementRestraint>"
);
#[wasm_bindgen]
impl BioRefinementRestraint {
    #[wasm_bindgen(constructor)]
    pub fn new(
        name: Option<Ts<input::TextInput>>,
        count: Option<Ts<input::I32Input>>,
        weight: Option<Ts<input::F64Input>>,
        function: Option<Ts<input::TextInput>>,
        dev_ideal: Option<Ts<input::F64Input>>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioRefinementRestraint::default();
        Ok(Self(ck::BioRefinementRestraint {
            name: name
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.name),
            count: count
                .map(input::read::<i32>)
                .transpose()?
                .unwrap_or(d.count),
            weight: weight
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.weight),
            function: function
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.function),
            dev_ideal: dev_ideal
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.dev_ideal),
        }))
    }
    #[wasm_bindgen(getter,js_name=name)]
    pub fn name(&self) -> String {
        self.0.name.clone()
    }
    #[wasm_bindgen(getter,js_name=count)]
    pub fn count(&self) -> i32 {
        self.0.count
    }
    #[wasm_bindgen(getter,js_name=weight)]
    pub fn weight(&self) -> f64 {
        self.0.weight
    }
    #[wasm_bindgen(getter,js_name=function)]
    pub fn function(&self) -> String {
        self.0.function.clone()
    }
    #[wasm_bindgen(getter,js_name=devIdeal)]
    pub fn dev_ideal(&self) -> f64 {
        self.0.dev_ideal
    }
}
#[wasm_bindgen]
pub struct BioDiffractionInfo(pub(crate) ck::BioDiffractionInfo);
crate::fields::value_field!(
    ck::BioDiffractionInfo,
    BioDiffractionInfo,
    BioDiffractionInfoInput,
    "BioDiffractionInfo"
);
crate::fields::value_array!(
    ck::BioDiffractionInfo,
    BioDiffractionInfoArray,
    "ReadonlyArray<BioDiffractionInfo>"
);
#[wasm_bindgen]
impl BioDiffractionInfo {
    #[wasm_bindgen(constructor)]
    pub fn new(
        id: Option<Ts<input::TextInput>>,
        temperature: Option<Ts<input::F64Input>>,
        source: Option<Ts<input::TextInput>>,
        source_type: Option<Ts<input::TextInput>>,
        synchrotron: Option<Ts<input::TextInput>>,
        beamline: Option<Ts<input::TextInput>>,
        wavelengths: Option<Ts<input::TextInput>>,
        scattering_type: Option<Ts<input::TextInput>>,
        mono_or_laue: Option<Ts<input::U8Input>>,
        monochromator: Option<Ts<input::TextInput>>,
        collection_date: Option<Ts<input::TextInput>>,
        optics: Option<Ts<input::TextInput>>,
        detector: Option<Ts<input::TextInput>>,
        detector_make: Option<Ts<input::TextInput>>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioDiffractionInfo::default();
        Ok(Self(ck::BioDiffractionInfo {
            id: id.map(input::read::<String>).transpose()?.unwrap_or(d.id),
            temperature: temperature
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.temperature),
            source: source
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.source),
            source_type: source_type
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.source_type),
            synchrotron: synchrotron
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.synchrotron),
            beamline: beamline
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.beamline),
            wavelengths: wavelengths
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.wavelengths),
            scattering_type: scattering_type
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.scattering_type),
            mono_or_laue: mono_or_laue
                .map(input::read::<u8>)
                .transpose()?
                .unwrap_or(d.mono_or_laue),
            monochromator: monochromator
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.monochromator),
            collection_date: collection_date
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.collection_date),
            optics: optics
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.optics),
            detector: detector
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.detector),
            detector_make: detector_make
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.detector_make),
        }))
    }
    #[wasm_bindgen(getter,js_name=id)]
    pub fn id(&self) -> String {
        self.0.id.clone()
    }
    #[wasm_bindgen(getter,js_name=temperature)]
    pub fn temperature(&self) -> f64 {
        self.0.temperature
    }
    #[wasm_bindgen(getter,js_name=source)]
    pub fn source(&self) -> String {
        self.0.source.clone()
    }
    #[wasm_bindgen(getter,js_name=sourceType)]
    pub fn source_type(&self) -> String {
        self.0.source_type.clone()
    }
    #[wasm_bindgen(getter,js_name=synchrotron)]
    pub fn synchrotron(&self) -> String {
        self.0.synchrotron.clone()
    }
    #[wasm_bindgen(getter,js_name=beamline)]
    pub fn beamline(&self) -> String {
        self.0.beamline.clone()
    }
    #[wasm_bindgen(getter,js_name=wavelengths)]
    pub fn wavelengths(&self) -> String {
        self.0.wavelengths.clone()
    }
    #[wasm_bindgen(getter,js_name=scatteringType)]
    pub fn scattering_type(&self) -> String {
        self.0.scattering_type.clone()
    }
    #[wasm_bindgen(getter,js_name=monoOrLaue)]
    pub fn mono_or_laue(&self) -> u8 {
        self.0.mono_or_laue
    }
    #[wasm_bindgen(getter,js_name=monochromator)]
    pub fn monochromator(&self) -> String {
        self.0.monochromator.clone()
    }
    #[wasm_bindgen(getter,js_name=collectionDate)]
    pub fn collection_date(&self) -> String {
        self.0.collection_date.clone()
    }
    #[wasm_bindgen(getter,js_name=optics)]
    pub fn optics(&self) -> String {
        self.0.optics.clone()
    }
    #[wasm_bindgen(getter,js_name=detector)]
    pub fn detector(&self) -> String {
        self.0.detector.clone()
    }
    #[wasm_bindgen(getter,js_name=detectorMake)]
    pub fn detector_make(&self) -> String {
        self.0.detector_make.clone()
    }
}
