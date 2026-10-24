//! Native residue vocabulary and information, without local recognition tables.
use crate::input;
use crate::{
    errors::{error_projection, exception},
    fields::Field,
    storage::NativeField,
};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::prelude::*;

#[wasm_bindgen(js_name = residueInfoKindName)]
pub fn residue_info_kind_name(value: ResidueInfoKindInput) -> Result<String, JsValue> {
    Ok(ck::ResidueInfoKind::input(value)?.name().into())
}
crate::enums::enum_value!(ResidueCode,"ResidueCode"{ALA="ALA",ARG="ARG",ASN="ASN",ABA="ABA",ASP="ASP",ASX="ASX",CYS="CYS",CSH="CSH",GLN="GLN",GLU="GLU",GLX="GLX",GLY="GLY",HIS="HIS",ILE="ILE",LEU="LEU",LYS="LYS",MET="MET",MSE="MSE",ORN="ORN",PHE="PHE",PRO="PRO",SER="SER",THR="THR",TRP="TRP",TYR="TYR",UNK="UNK",VAL="VAL",SEC="SEC",PYL="PYL",SEP="SEP",TPO="TPO",PCA="PCA",CSO="CSO",PTR="PTR",KCX="KCX",CSD="CSD",LLP="LLP",CME="CME",MLY="MLY",DAL="DAL",TYS="TYS",OCS="OCS",M3L="M3L",FME="FME",ALY="ALY",HYP="HYP",CAS="CAS",CRO="CRO",CSX="CSX",DPR="DPR",DGL="DGL",DVA="DVA",CSS="CSS",DPN="DPN",DSN="DSN",DLE="DLE",HIC="HIC",NLE="NLE",MVA="MVA",MLZ="MLZ",CR2="CR2",SAR="SAR",DAR="DAR",DLY="DLY",YCM="YCM",NRQ="NRQ",CGU="CGU",R0TD="R0TD",MLE="MLE",DAS="DAS",DTR="DTR",CXM="CXM",TPQ="TPQ",DCY="DCY",DSG="DSG",DTY="DTY",DHI="DHI",MEN="MEN",DTH="DTH",SAC="SAC",DGN="DGN",AIB="AIB",SMC="SMC",IAS="IAS",CIR="CIR",BMT="BMT",DIL="DIL",FGA="FGA",PHI="PHI",CRQ="CRQ",SME="SME",GHP="GHP",MHO="MHO",NEP="NEP",TRQ="TRQ",TOX="TOX",ALC="ALC",R3FG="R3FG",SCH="SCH",MDO="MDO",MAA="MAA",GYS="GYS",MK8="MK8",CR8="CR8",KPI="KPI",SCY="SCY",DHA="DHA",OMY="OMY",CAF="CAF",R0AF="R0AF",SNN="SNN",MHS="MHS",MLU="MLU",SNC="SNC",PHD="PHD",B3E="B3E",MEA="MEA",MED="MED",OAS="OAS",GL3="GL3",FVA="FVA",PHL="PHL",CRF="CRF",OMZ="OMZ",BFD="BFD",MEQ="MEQ",DAB="DAB",AGM="AGM",PSU="PSU",R5MU="R5MU",R7MG="R7MG",OMG="OMG",UR3="UR3",OMC="OMC",R2MG="R2MG",H2U="H2U",R4SU="R4SU",OMU="OMU",R4OC="R4OC",MA6="MA6",M2G="M2G",R1MA="R1MA",R6MZ="R6MZ",CCC="CCC",R2MA="R2MA",R1MG="R1MG",R5BU="R5BU",MIA="MIA",DOC="DOC",R8OG="R8OG",R5CM="R5CM",R3DR="R3DR",BRU="BRU",CBR="CBR",HOH="HOH",DOD="DOD",HEM="HEM",SO4="SO4",GOL="GOL",EDO="EDO",NAG="NAG",PO4="PO4",ACT="ACT",PEG="PEG",MAN="MAN",FAD="FAD",BMA="BMA",ADP="ADP",DMS="DMS",ACE="ACE",NH2="NH2",MPD="MPD",MES="MES",NAD="NAD",NAP="NAP",TRS="TRS",ATP="ATP",PG4="PG4",GDP="GDP",FUC="FUC",FMT="FMT",GAL="GAL",PGE="PGE",FMN="FMN",PLP="PLP",EPE="EPE",SF4="SF4",BME="BME",CIT="CIT",BE7="BE7",MRD="MRD",MHA="MHA",BU3="BU3",PGO="PGO",BU2="BU2",PDO="PDO",BU1="BU1",PG6="PG6",R1BO="R1BO",PE7="PE7",PG5="PG5",TFP="TFP",DHD="DHD",PEU="PEU",TAU="TAU",SBT="SBT",SAL="SAL",IOH="IOH",IPA="IPA",PIG="PIG",B3P="B3P",BTB="BTB",NHE="NHE",C8E="C8E",OTE="OTE",PE4="PE4",XPE="XPE",PE8="PE8",P33="P33",N8E="N8E",R2OS="R2OS",R1PS="R1PS",CPS="CPS",DMX="DMX",MPO="MPO",GCD="GCD",DXG="DXG",CM5="CM5",ACA="ACA",ACN="ACN",CCN="CCN",GLC="GLC",DR6="DR6",NH4="NH4",AZI="AZI",BNG="BNG",BOG="BOG",BGC="BGC",BCN="BCN",BRO="BRO",CAC="CAC",CBX="CBX",ACY="ACY",CBM="CBM",CLO="CLO",R3CO="R3CO",NCO="NCO",CU1="CU1",CYN="CYN",MA4="MA4",TAR="TAR",GLO="GLO",MTL="MTL",SOR="SOR",DMU="DMU",DDQ="DDQ",DMF="DMF",DIO="DIO",DOX="DOX",R12P="R12P",SDS="SDS",LMT="LMT",EOH="EOH",EEE="EEE",EGL="EGL",FLO="FLO",TRT="TRT",FCY="FCY",FRU="FRU",GBL="GBL",GPX="GPX",HTO="HTO",HTG="HTG",B7G="B7G",C10="C10",R16D="R16D",HEZ="HEZ",IOD="IOD",IDO="IDO",ICI="ICI",ICT="ICT",TLA="TLA",LAT="LAT",LBT="LBT",LDA="LDA",MN3="MN3",MRY="MRY",MOH="MOH",BEQ="BEQ",C15="C15",MG8="MG8",POL="POL",NO3="NO3",JEF="JEF",P4C="P4C",CE1="CE1",DIA="DIA",CXE="CXE",IPH="IPH",PIN="PIN",R15P="R15P",CRY="CRY",PGR="PGR",PGQ="PGQ",SPD="SPD",SPK="SPK",SPM="SPM",SUC="SUC",TBU="TBU",TMA="TMA",TEP="TEP",SCN="SCN",TRE="TRE",ETF="ETF",R144="R144",UMQ="UMQ",URE="URE",YT3="YT3",ZN2="ZN2",FE2="FE2",R3NI="R3NI",SIA="SIA",XYP="XYP",A2G="A2G",GLA="GLA",NDG="NDG",NGA="NGA",A="A",C="C",G="G",I="I",U="U",N="N",F="F",K="K",DA="DA",DC="DC",DG="DG",DI="DI",DT="DT",DU="DU",DN="DN",AG="AG",AL="AL",BA="BA",BR="BR",CA="CA",CD="CD",CL="CL",CM="CM",CN="CN",CO="CO",CS="CS",CU="CU",FE="FE",HG="HG",LI="LI",MG="MG",MN="MN",NA="NA",NI="NI",NO="NO",PB="PB",RB="RB",SR="SR",Y1="Y1",ZN="ZN",UNKNOWN="UNKNOWN"});
crate::enums::enum_value!(ResidueInfoKind,"ResidueInfoKind"{Unknown="unknown",Aa="aa",Aad="aad",Paa="paa",Maa="maa",Rna="rna",Dna="dna",Buf="buf",Hoh="hoh",Pyr="pyr",Ket="ket",Els="els"});
#[wasm_bindgen]
pub struct ResidueInfo(pub(crate) ck::ResidueInfo);
crate::fields::value_field!(
    ck::ResidueInfo,
    ResidueInfo,
    ResidueInfoInput,
    "ResidueInfo"
);
#[wasm_bindgen]
impl ResidueInfo {
    #[wasm_bindgen(getter)]
    pub fn code(&self) -> ResidueCode {
        self.0.code.into()
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.0.name.into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> ResidueInfoKind {
        self.0.kind.into()
    }
    #[wasm_bindgen(getter,js_name=linkingType)]
    pub fn linking_type(&self) -> u8 {
        self.0.linking_type
    }
    #[wasm_bindgen(getter,js_name=oneLetterCode)]
    pub fn one_letter_code(&self) -> char {
        self.0.one_letter_code
    }
    #[wasm_bindgen(getter,js_name=hydrogenCount)]
    pub fn hydrogen_count(&self) -> u8 {
        self.0.hydrogen_count
    }
    #[wasm_bindgen(getter)]
    pub fn weight(&self) -> f32 {
        self.0.weight
    }
    #[wasm_bindgen(getter,js_name=canonicalOneLetterCode)]
    pub fn canonical_one_letter_code(&self) -> Option<char> {
        self.0.canonical_one_letter_code()
    }
    #[wasm_bindgen(getter,js_name=parentStandardCode)]
    pub fn parent_standard_code(&self) -> Option<ResidueCode> {
        self.0.parent_standard_code().map(Into::into)
    }
    #[wasm_bindgen(getter,js_name=isModifiedAminoAcid)]
    pub fn is_modified_amino_acid(&self) -> bool {
        self.0.is_modified_amino_acid()
    }
    #[wasm_bindgen(getter)]
    pub fn found(&self) -> bool {
        self.0.found()
    }
    #[wasm_bindgen(getter,js_name=isWater)]
    pub fn is_water(&self) -> bool {
        self.0.is_water()
    }
    #[wasm_bindgen(getter,js_name=isDna)]
    pub fn is_dna(&self) -> bool {
        self.0.is_dna()
    }
    #[wasm_bindgen(getter,js_name=isRna)]
    pub fn is_rna(&self) -> bool {
        self.0.is_rna()
    }
}
#[wasm_bindgen]
pub struct ResidueIdentity(ck::ResidueIdentity);
crate::fields::value_field!(
    ck::ResidueIdentity,
    ResidueIdentity,
    ResidueIdentityInput,
    "ResidueIdentity"
);
#[wasm_bindgen]
impl ResidueIdentity {
    #[wasm_bindgen(constructor)]
    pub fn new(name: Ts<input::TextInput>) -> Result<Self, JsValue> {
        Ok(Self(ck::ResidueIdentity::new(input::read::<String>(name)?)))
    }
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        self.0.name().into()
    }
    #[wasm_bindgen(getter)]
    pub fn code(&self) -> ResidueCode {
        self.0.code().into()
    }
    #[wasm_bindgen(getter)]
    pub fn info(&self) -> ResidueInfo {
        ResidueInfo(self.0.info())
    }
    #[wasm_bindgen(getter,js_name=isTabulated)]
    pub fn is_tabulated(&self) -> bool {
        self.0.is_tabulated()
    }
}
#[wasm_bindgen(js_name=residueInfoChecked)]
pub fn residue_info_checked(index: Ts<input::IndexInput>) -> Result<Option<ResidueInfo>, JsValue> {
    Ok(ck::residue_info_checked(input::read::<usize>(index)?).map(ResidueInfo))
}
#[wasm_bindgen(js_name=residueInfo)]
pub fn residue_info(index: Ts<input::IndexInput>) -> Result<ResidueInfo, JsValue> {
    residue_info_checked(index)?
        .ok_or_else(|| js_sys::RangeError::new("residue information index out of range").into())
}
#[wasm_bindgen(js_name=findResidueInfoIndex)]
pub fn find_residue_info_index(name: Ts<input::TextInput>) -> Result<usize, JsValue> {
    Ok(ck::find_residue_info_index(&input::read::<String>(name)?))
}
#[wasm_bindgen(js_name=findResidueInfo)]
pub fn find_residue_info(name: Ts<input::TextInput>) -> Result<ResidueInfo, JsValue> {
    Ok(ResidueInfo(ck::find_residue_info(&input::read::<String>(
        name,
    )?)))
}
#[wasm_bindgen(js_name=residueCode)]
pub fn residue_code(name: Ts<input::TextInput>) -> Result<ResidueCode, JsValue> {
    Ok(ck::residue_code(&input::read::<String>(name)?).into())
}
crate::errors::enum_error!(ResidueSequenceError,ResidueSequenceErrorKind,ck::ResidueSequenceError,{UnmatchedParenthesis=>ck::ResidueSequenceError::UnmatchedParenthesis,UnexpectedLetter=>ck::ResidueSequenceError::UnexpectedLetter{..}},{
#[wasm_bindgen(getter)]pub fn kind(&self)->Option<String>{match self.inner.as_ref(){ck::ResidueSequenceError::UnexpectedLetter{kind,..}=>Some((*kind).into()),_=>None}}
#[wasm_bindgen(getter)]pub fn letter(&self)->Option<char>{match self.inner.as_ref(){ck::ResidueSequenceError::UnexpectedLetter{letter,..}=>Some(*letter),_=>None}}
#[wasm_bindgen(getter,js_name=sourceCode)]pub fn source_code(&self)->Option<i16>{match self.inner.as_ref(){ck::ResidueSequenceError::UnexpectedLetter{source_code,..}=>Some(*source_code),_=>None}}
});
error_projection!(ResidueCodeParseError, ck::ResidueCodeParseError, {
    #[wasm_bindgen(getter)]
    pub fn input(&self) -> String {
        self.inner.input().into()
    }
});
#[wasm_bindgen(js_name=expandOneLetter)]
pub fn expand_one_letter(
    code: Ts<input::CharInput>,
    kind: ResidueInfoKindInput,
) -> Result<Option<String>, JsValue> {
    Ok(ck::expand_one_letter(
        input::read::<char>(code)?,
        <ck::ResidueInfoKind as Field>::input(kind)?,
    )
    .map(str::to_owned))
}
#[wasm_bindgen(js_name=expandOneLetterSequence)]
pub fn expand_one_letter_sequence(
    sequence: Ts<input::TextInput>,
    kind: ResidueInfoKindInput,
) -> Result<Ts<crate::fields::TextArray>, JsValue> {
    let v = ck::expand_one_letter_sequence(
        &input::read::<String>(sequence)?,
        <ck::ResidueInfoKind as Field>::input(kind)?,
    )
    .map_err(ResidueSequenceError::from_rust)?;
    <Vec<String> as Field>::get(NativeField::owned(v))
}
