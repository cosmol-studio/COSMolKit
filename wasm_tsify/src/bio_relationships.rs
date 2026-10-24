//! Source-owned structural annotations and explicit addresses.
use crate::{bio_hierarchy::BioStructure, fields::Field, input, storage::NativeField};
use cosmolkit as ck;
use tsify::{Ts, Tsify};
use wasm_bindgen::prelude::*;
crate::enums::enum_value!(BioConnectionKind,"BioConnectionKind"{Covale="covale",Disulf="disulf",Hydrog="hydrog",MetalC="metal_c",Unknown="unknown"});
crate::enums::enum_value!(BioAsu,"BioAsu"{Same="same",Different="different",Any="any"});
crate::enums::enum_value!(BioHelixClass,"BioHelixClass"{UnknownHelix="unknown_helix",RAlpha="r_alpha",ROmega="r_omega",RPi="r_pi",RGamma="r_gamma",R310="r_310",LAlpha="l_alpha",LOmega="l_omega",LGamma="l_gamma",Helix27="helix_27",HelixPolyProlineNone="helix_poly_proline_none"});
#[derive(serde::Serialize, serde::Deserialize, Tsify)]
#[serde(transparent)]
pub struct ReportedSym(#[serde(deserialize_with = "input::exact_array")] pub [i16; 4]);
#[wasm_bindgen]
pub struct ResidueAddress(ck::ResidueAddress);
crate::fields::value_field!(
    ck::ResidueAddress,
    ResidueAddress,
    ResidueAddressInput,
    "ResidueAddress"
);
#[wasm_bindgen]
impl ResidueAddress {
    pub fn new(
        sequence_number: Option<Ts<input::I32Input>>,
        insertion_code: Option<Ts<input::U8Input>>,
        segment: Ts<input::TextInput>,
        name: Ts<input::TextInput>,
    ) -> Result<Option<Self>, JsValue> {
        let sequence_number = sequence_number.map(input::read::<i32>).transpose()?;
        let insertion_code = insertion_code.map(input::read::<u8>).transpose()?;
        let segment = input::read::<String>(segment)?;
        let name = input::read::<String>(name)?;
        Ok(ck::ResidueName::from_ascii(name.as_bytes())
            .and_then(|name| {
                ck::ResidueAddress::new(sequence_number, insertion_code, segment.as_bytes(), name)
            })
            .map(Self))
    }
    #[wasm_bindgen(getter,js_name=sequenceNumber)]
    pub fn sequence_number(&self) -> Option<i32> {
        self.0.sequence_number()
    }
    #[wasm_bindgen(getter,js_name=insertionCode)]
    pub fn insertion_code(&self) -> Option<u8> {
        self.0.insertion_code()
    }
    #[wasm_bindgen(getter)]
    pub fn segment(&self) -> String {
        self.0.segment().into()
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.0.name().as_str().into()
    }
    pub fn matches(&self, other: ResidueAddressInput) -> Result<bool, JsValue> {
        Self::with_ref(&other, |v| Ok(self.0.matches(v)))
    }
    #[wasm_bindgen(js_name=matchesWithoutSegment)]
    pub fn matches_without_segment(&self, other: ResidueAddressInput) -> Result<bool, JsValue> {
        Self::with_ref(&other, |v| Ok(self.0.matches_without_segment(v)))
    }
}
#[wasm_bindgen]
pub struct AtomAddress(ck::AtomAddress);
crate::fields::value_field!(
    ck::AtomAddress,
    AtomAddress,
    AtomAddressInput,
    "AtomAddress"
);
#[wasm_bindgen]
impl AtomAddress {
    #[wasm_bindgen(constructor)]
    pub fn new(
        chain_name: Ts<input::TextInput>,
        residue: ResidueAddressInput,
        logical_atom_name: Ts<input::TextInput>,
        altloc: Option<Ts<input::U8Input>>,
    ) -> Result<Self, JsValue> {
        let chain = ck::PdbChainId::from_ascii(input::read::<String>(chain_name)?.as_bytes())
            .ok_or_else(|| {
                JsValue::from(js_sys::RangeError::new(
                    "chain_name is not a valid PdbChainId",
                ))
            })?;
        Ok(Self(ck::AtomAddress::new(
            chain,
            <ck::ResidueAddress as Field>::input(residue)?,
            input::read::<String>(logical_atom_name)?,
            altloc.map(input::read::<u8>).transpose()?,
        )))
    }
    #[wasm_bindgen(getter,js_name=chainName)]
    pub fn chain_name(&self) -> String {
        self.0.chain_name().as_str().into()
    }
    #[wasm_bindgen(getter)]
    pub fn residue(&self) -> ResidueAddress {
        ResidueAddress(self.0.residue())
    }
    #[wasm_bindgen(getter,js_name=logicalAtomName)]
    pub fn logical_atom_name(&self) -> String {
        self.0.logical_atom_name().into()
    }
    #[wasm_bindgen(getter)]
    pub fn altloc(&self) -> u8 {
        self.0.altloc()
    }
}
#[wasm_bindgen]
pub struct BioConnection(ck::BioConnection);
crate::fields::value_field!(
    ck::BioConnection,
    BioConnection,
    BioConnectionInput,
    "BioConnection"
);
crate::fields::value_array!(
    ck::BioConnection,
    BioConnections,
    "ReadonlyArray<BioConnection>"
);
#[wasm_bindgen]
impl BioConnection {
    #[wasm_bindgen(constructor)]
    pub fn new(
        name: Option<Ts<input::TextInput>>,
        link_id: Option<Ts<input::TextInput>>,
        kind: Option<BioConnectionKindInput>,
        asu: Option<BioAsuInput>,
        partner1: Option<AtomAddressInput>,
        partner2: Option<AtomAddressInput>,
        reported_distance: Option<Ts<input::F64Input>>,
        reported_sym: Option<Ts<ReportedSym>>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioConnection::default();
        Ok(Self(ck::BioConnection {
            name: name
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.name),
            link_id: link_id
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.link_id),
            kind: kind
                .map(<ck::BioConnectionKind as Field>::input)
                .transpose()?
                .unwrap_or(d.kind),
            asu: asu
                .map(<ck::BioAsu as Field>::input)
                .transpose()?
                .unwrap_or(d.asu),
            partner1: partner1
                .map(<ck::AtomAddress as Field>::input)
                .transpose()?
                .unwrap_or(d.partner1),
            partner2: partner2
                .map(<ck::AtomAddress as Field>::input)
                .transpose()?
                .unwrap_or(d.partner2),
            reported_distance: reported_distance
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.reported_distance),
            reported_sym: reported_sym
                .map(input::decode)
                .transpose()?
                .map(|v| v.0)
                .unwrap_or(d.reported_sym),
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.0.name.clone()
    }
    #[wasm_bindgen(getter,js_name=linkId)]
    pub fn link_id(&self) -> String {
        self.0.link_id.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> BioConnectionKind {
        self.0.kind.into()
    }
    #[wasm_bindgen(getter)]
    pub fn asu(&self) -> BioAsu {
        self.0.asu.into()
    }
    #[wasm_bindgen(getter)]
    pub fn partner1(&self) -> AtomAddress {
        AtomAddress(self.0.partner1.clone())
    }
    #[wasm_bindgen(getter)]
    pub fn partner2(&self) -> AtomAddress {
        AtomAddress(self.0.partner2.clone())
    }
    #[wasm_bindgen(getter,js_name=reportedDistance)]
    pub fn reported_distance(&self) -> f64 {
        self.0.reported_distance
    }
    #[wasm_bindgen(getter,js_name=reportedSym)]
    pub fn reported_sym(&self) -> Result<Ts<ReportedSym>, JsValue> {
        ReportedSym(self.0.reported_sym)
            .into_ts()
            .map_err(JsError::from)
            .map_err(Into::into)
    }
}
#[wasm_bindgen]
pub struct BioCisPep(ck::BioCisPep);
crate::fields::value_field!(ck::BioCisPep, BioCisPep, BioCisPepInput, "BioCisPep");
crate::fields::value_array!(ck::BioCisPep, BioCisPeps, "ReadonlyArray<BioCisPep>");
#[wasm_bindgen]
impl BioCisPep {
    #[wasm_bindgen(constructor)]
    pub fn new(
        partner_c: Option<AtomAddressInput>,
        partner_n: Option<AtomAddressInput>,
        model_num: Option<Ts<input::I32Input>>,
        only_altloc: Option<Ts<input::U8Input>>,
        reported_angle: Option<Ts<input::F64Input>>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioCisPep::default();
        Ok(Self(ck::BioCisPep {
            partner_c: partner_c
                .map(<ck::AtomAddress as Field>::input)
                .transpose()?
                .unwrap_or(d.partner_c),
            partner_n: partner_n
                .map(<ck::AtomAddress as Field>::input)
                .transpose()?
                .unwrap_or(d.partner_n),
            model_num: model_num
                .map(input::read::<i32>)
                .transpose()?
                .unwrap_or(d.model_num),
            only_altloc: only_altloc
                .map(input::read::<u8>)
                .transpose()?
                .unwrap_or(d.only_altloc),
            reported_angle: reported_angle
                .map(input::read::<f64>)
                .transpose()?
                .unwrap_or(d.reported_angle),
        }))
    }
    #[wasm_bindgen(getter,js_name=partnerC)]
    pub fn partner_c(&self) -> AtomAddress {
        AtomAddress(self.0.partner_c.clone())
    }
    #[wasm_bindgen(getter,js_name=partnerN)]
    pub fn partner_n(&self) -> AtomAddress {
        AtomAddress(self.0.partner_n.clone())
    }
    #[wasm_bindgen(getter,js_name=modelNum)]
    pub fn model_num(&self) -> i32 {
        self.0.model_num
    }
    #[wasm_bindgen(getter,js_name=onlyAltloc)]
    pub fn only_altloc(&self) -> u8 {
        self.0.only_altloc
    }
    #[wasm_bindgen(getter,js_name=reportedAngle)]
    pub fn reported_angle(&self) -> f64 {
        self.0.reported_angle
    }
}
#[wasm_bindgen]
pub struct BioModRes(ck::BioModRes);
crate::fields::value_field!(ck::BioModRes, BioModRes, BioModResInput, "BioModRes");
crate::fields::value_array!(ck::BioModRes, BioModResidues, "ReadonlyArray<BioModRes>");
#[wasm_bindgen]
impl BioModRes {
    #[wasm_bindgen(constructor)]
    pub fn new(
        chain_name: Option<Ts<input::TextInput>>,
        res_id: Option<ResidueAddressInput>,
        parent_comp_id: Option<Ts<input::TextInput>>,
        mod_id: Option<Ts<input::TextInput>>,
        details: Option<Ts<input::TextInput>>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioModRes::default();
        let chain_name = chain_name
            .map(input::read::<String>)
            .transpose()?
            .map(|v| {
                ck::PdbChainId::from_ascii(v.as_bytes()).ok_or_else(|| {
                    JsValue::from(js_sys::RangeError::new(
                        "chain_name is not a valid PdbChainId",
                    ))
                })
            })
            .transpose()?
            .unwrap_or(d.chain_name);
        Ok(Self(ck::BioModRes {
            chain_name,
            res_id: res_id
                .map(<ck::ResidueAddress as Field>::input)
                .transpose()?
                .unwrap_or(d.res_id),
            parent_comp_id: parent_comp_id
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.parent_comp_id),
            mod_id: mod_id
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.mod_id),
            details: details
                .map(input::read::<String>)
                .transpose()?
                .unwrap_or(d.details),
        }))
    }
    #[wasm_bindgen(getter,js_name=chainName)]
    pub fn chain_name(&self) -> String {
        self.0.chain_name.as_str().into()
    }
    #[wasm_bindgen(getter,js_name=resId)]
    pub fn res_id(&self) -> ResidueAddress {
        ResidueAddress(self.0.res_id)
    }
    #[wasm_bindgen(getter,js_name=parentCompId)]
    pub fn parent_comp_id(&self) -> String {
        self.0.parent_comp_id.clone()
    }
    #[wasm_bindgen(getter,js_name=modId)]
    pub fn mod_id(&self) -> String {
        self.0.mod_id.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn details(&self) -> String {
        self.0.details.clone()
    }
}
#[wasm_bindgen]
pub struct BioHelix(ck::BioHelix);
crate::fields::value_field!(ck::BioHelix, BioHelix, BioHelixInput, "BioHelix");
crate::fields::value_array!(ck::BioHelix, BioHelices, "ReadonlyArray<BioHelix>");
#[wasm_bindgen]
impl BioHelix {
    #[wasm_bindgen(constructor)]
    pub fn new(
        start: Option<AtomAddressInput>,
        end: Option<AtomAddressInput>,
        pdb_helix_class: Option<BioHelixClassInput>,
        length: Option<Ts<input::I32Input>>,
    ) -> Result<Self, JsValue> {
        let d = ck::BioHelix::default();
        Ok(Self(ck::BioHelix {
            start: start
                .map(<ck::AtomAddress as Field>::input)
                .transpose()?
                .unwrap_or(d.start),
            end: end
                .map(<ck::AtomAddress as Field>::input)
                .transpose()?
                .unwrap_or(d.end),
            pdb_helix_class: pdb_helix_class
                .map(<ck::BioHelixClass as Field>::input)
                .transpose()?
                .unwrap_or(d.pdb_helix_class),
            length: length
                .map(input::read::<i32>)
                .transpose()?
                .unwrap_or(d.length),
        }))
    }
    #[wasm_bindgen(getter)]
    pub fn start(&self) -> AtomAddress {
        AtomAddress(self.0.start.clone())
    }
    #[wasm_bindgen(getter)]
    pub fn end(&self) -> AtomAddress {
        AtomAddress(self.0.end.clone())
    }
    #[wasm_bindgen(getter,js_name=pdbHelixClass)]
    pub fn pdb_helix_class(&self) -> BioHelixClass {
        self.0.pdb_helix_class.into()
    }
    #[wasm_bindgen(getter)]
    pub fn length(&self) -> i32 {
        self.0.length
    }
    #[wasm_bindgen(js_name=setHelixClassAsInt)]
    pub fn set_helix_class_as_int(&mut self, n: Ts<input::I32Input>) -> Result<(), JsValue> {
        self.0.set_helix_class_as_int(input::read::<i32>(n)?);
        Ok(())
    }
}
#[wasm_bindgen]
pub struct BioStrand(ck::BioStrand);
crate::fields::value_field!(ck::BioStrand, BioStrand, BioStrandInput, "BioStrand");
crate::fields::value_array!(ck::BioStrand, BioStrands, "ReadonlyArray<BioStrand>");
#[wasm_bindgen]
impl BioStrand {
    #[wasm_bindgen(constructor)]
    pub fn new(
        start: AtomAddressInput,
        end: AtomAddressInput,
        hbond_atom2: AtomAddressInput,
        hbond_atom1: AtomAddressInput,
        sense: Ts<input::I32Input>,
        name: Ts<input::TextInput>,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::BioStrand::new(
            <ck::AtomAddress as Field>::input(start)?,
            <ck::AtomAddress as Field>::input(end)?,
            <ck::AtomAddress as Field>::input(hbond_atom2)?,
            <ck::AtomAddress as Field>::input(hbond_atom1)?,
            input::read::<i32>(sense)?,
            input::read::<String>(name)?,
        )))
    }
    #[wasm_bindgen(getter)]
    pub fn start(&self) -> AtomAddress {
        AtomAddress(self.0.start.clone())
    }
    #[wasm_bindgen(getter)]
    pub fn end(&self) -> AtomAddress {
        AtomAddress(self.0.end.clone())
    }
    #[wasm_bindgen(getter,js_name=hbondAtom2)]
    pub fn hbond_atom2(&self) -> AtomAddress {
        AtomAddress(self.0.hbond_atom2.clone())
    }
    #[wasm_bindgen(getter,js_name=hbondAtom1)]
    pub fn hbond_atom1(&self) -> AtomAddress {
        AtomAddress(self.0.hbond_atom1.clone())
    }
    #[wasm_bindgen(getter)]
    pub fn sense(&self) -> i32 {
        self.0.sense
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.0.name.clone()
    }
}
#[wasm_bindgen]
pub struct BioSheet(ck::BioSheet);
crate::fields::value_field!(ck::BioSheet, BioSheet, BioSheetInput, "BioSheet");
crate::fields::value_array!(ck::BioSheet, BioSheets, "ReadonlyArray<BioSheet>");
#[wasm_bindgen]
impl BioSheet {
    #[wasm_bindgen(constructor)]
    pub fn new(
        sheet_id: Ts<input::TextInput>,
        strands: Option<BioStrands>,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::BioSheet::new(&input::read::<String>(sheet_id)?);
        if let Some(strands) = strands {
            inner.strands = <Vec<ck::BioStrand> as Field>::input(strands)?;
        }
        Ok(Self(inner))
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.0.name.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn strands(&self) -> Result<BioStrands, JsValue> {
        <Vec<ck::BioStrand> as Field>::get(NativeField::owned(self.0.strands.clone()))
    }
}
#[wasm_bindgen]
impl BioStructure {
    #[wasm_bindgen(getter)]
    pub fn connections(&self) -> Result<BioConnections, JsValue> {
        <Vec<ck::BioConnection> as Field>::get(NativeField::owned(self.0.connections().to_vec()))
    }
    #[wasm_bindgen(getter)]
    pub fn cispeps(&self) -> Result<BioCisPeps, JsValue> {
        <Vec<ck::BioCisPep> as Field>::get(NativeField::owned(self.0.cispeps().to_vec()))
    }
    #[wasm_bindgen(getter,js_name=modResidues)]
    pub fn mod_residues(&self) -> Result<BioModResidues, JsValue> {
        <Vec<ck::BioModRes> as Field>::get(NativeField::owned(self.0.mod_residues().to_vec()))
    }
    #[wasm_bindgen(getter)]
    pub fn helices(&self) -> Result<BioHelices, JsValue> {
        <Vec<ck::BioHelix> as Field>::get(NativeField::owned(self.0.helices().to_vec()))
    }
    #[wasm_bindgen(getter)]
    pub fn sheets(&self) -> Result<BioSheets, JsValue> {
        <Vec<ck::BioSheet> as Field>::get(NativeField::owned(self.0.sheets().to_vec()))
    }
}
