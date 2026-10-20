//! Element values and periodic metadata from the public Rust facade.
use cosmolkit as ck;
use pyo3::prelude::*;

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct Element {
        pub(crate) inner: ck::Element,
    }
}

#[pymethods]
impl Element {
    #[classattr]
    const DUMMY: Self = Self {
        inner: ck::Element::DUMMY,
    };
    #[classattr]
    const H: Self = Self {
        inner: ck::Element::H,
    };
    #[classattr]
    const HE: Self = Self {
        inner: ck::Element::HE,
    };
    #[classattr]
    const LI: Self = Self {
        inner: ck::Element::LI,
    };
    #[classattr]
    const BE: Self = Self {
        inner: ck::Element::BE,
    };
    #[classattr]
    const B: Self = Self {
        inner: ck::Element::B,
    };
    #[classattr]
    const N: Self = Self {
        inner: ck::Element::N,
    };
    #[classattr]
    const F: Self = Self {
        inner: ck::Element::F,
    };
    #[classattr]
    const NE: Self = Self {
        inner: ck::Element::NE,
    };
    #[classattr]
    const NA: Self = Self {
        inner: ck::Element::NA,
    };
    #[classattr]
    const MG: Self = Self {
        inner: ck::Element::MG,
    };
    #[classattr]
    const AL: Self = Self {
        inner: ck::Element::AL,
    };
    #[classattr]
    const SI: Self = Self {
        inner: ck::Element::SI,
    };
    #[classattr]
    const P: Self = Self {
        inner: ck::Element::P,
    };
    #[classattr]
    const S: Self = Self {
        inner: ck::Element::S,
    };
    #[classattr]
    const AR: Self = Self {
        inner: ck::Element::AR,
    };
    #[classattr]
    const K: Self = Self {
        inner: ck::Element::K,
    };
    #[classattr]
    const CA: Self = Self {
        inner: ck::Element::CA,
    };
    #[classattr]
    const SC: Self = Self {
        inner: ck::Element::SC,
    };
    #[classattr]
    const TI: Self = Self {
        inner: ck::Element::TI,
    };
    #[classattr]
    const V: Self = Self {
        inner: ck::Element::V,
    };
    #[classattr]
    const CR: Self = Self {
        inner: ck::Element::CR,
    };
    #[classattr]
    const MN: Self = Self {
        inner: ck::Element::MN,
    };
    #[classattr]
    const FE: Self = Self {
        inner: ck::Element::FE,
    };
    #[classattr]
    const CO: Self = Self {
        inner: ck::Element::CO,
    };
    #[classattr]
    const NI: Self = Self {
        inner: ck::Element::NI,
    };
    #[classattr]
    const CU: Self = Self {
        inner: ck::Element::CU,
    };
    #[classattr]
    const ZN: Self = Self {
        inner: ck::Element::ZN,
    };
    #[classattr]
    const GA: Self = Self {
        inner: ck::Element::GA,
    };
    #[classattr]
    const GE: Self = Self {
        inner: ck::Element::GE,
    };
    #[classattr]
    const AS: Self = Self {
        inner: ck::Element::AS,
    };
    #[classattr]
    const SE: Self = Self {
        inner: ck::Element::SE,
    };
    #[classattr]
    const BR: Self = Self {
        inner: ck::Element::BR,
    };
    #[classattr]
    const KR: Self = Self {
        inner: ck::Element::KR,
    };
    #[classattr]
    const RB: Self = Self {
        inner: ck::Element::RB,
    };
    #[classattr]
    const SR: Self = Self {
        inner: ck::Element::SR,
    };
    #[classattr]
    const Y: Self = Self {
        inner: ck::Element::Y,
    };
    #[classattr]
    const ZR: Self = Self {
        inner: ck::Element::ZR,
    };
    #[classattr]
    const NB: Self = Self {
        inner: ck::Element::NB,
    };
    #[classattr]
    const MO: Self = Self {
        inner: ck::Element::MO,
    };
    #[classattr]
    const TC: Self = Self {
        inner: ck::Element::TC,
    };
    #[classattr]
    const RU: Self = Self {
        inner: ck::Element::RU,
    };
    #[classattr]
    const RH: Self = Self {
        inner: ck::Element::RH,
    };
    #[classattr]
    const PD: Self = Self {
        inner: ck::Element::PD,
    };
    #[classattr]
    const AG: Self = Self {
        inner: ck::Element::AG,
    };
    #[classattr]
    const CD: Self = Self {
        inner: ck::Element::CD,
    };
    #[classattr]
    const IN: Self = Self {
        inner: ck::Element::IN,
    };
    #[classattr]
    const SN: Self = Self {
        inner: ck::Element::SN,
    };
    #[classattr]
    const SB: Self = Self {
        inner: ck::Element::SB,
    };
    #[classattr]
    const TE: Self = Self {
        inner: ck::Element::TE,
    };
    #[classattr]
    const I: Self = Self {
        inner: ck::Element::I,
    };
    #[classattr]
    const XE: Self = Self {
        inner: ck::Element::XE,
    };
    #[classattr]
    const CS: Self = Self {
        inner: ck::Element::CS,
    };
    #[classattr]
    const BA: Self = Self {
        inner: ck::Element::BA,
    };
    #[classattr]
    const LA: Self = Self {
        inner: ck::Element::LA,
    };
    #[classattr]
    const CE: Self = Self {
        inner: ck::Element::CE,
    };
    #[classattr]
    const PR: Self = Self {
        inner: ck::Element::PR,
    };
    #[classattr]
    const ND: Self = Self {
        inner: ck::Element::ND,
    };
    #[classattr]
    const PM: Self = Self {
        inner: ck::Element::PM,
    };
    #[classattr]
    const SM: Self = Self {
        inner: ck::Element::SM,
    };
    #[classattr]
    const EU: Self = Self {
        inner: ck::Element::EU,
    };
    #[classattr]
    const GD: Self = Self {
        inner: ck::Element::GD,
    };
    #[classattr]
    const TB: Self = Self {
        inner: ck::Element::TB,
    };
    #[classattr]
    const DY: Self = Self {
        inner: ck::Element::DY,
    };
    #[classattr]
    const HO: Self = Self {
        inner: ck::Element::HO,
    };
    #[classattr]
    const ER: Self = Self {
        inner: ck::Element::ER,
    };
    #[classattr]
    const TM: Self = Self {
        inner: ck::Element::TM,
    };
    #[classattr]
    const YB: Self = Self {
        inner: ck::Element::YB,
    };
    #[classattr]
    const LU: Self = Self {
        inner: ck::Element::LU,
    };
    #[classattr]
    const HF: Self = Self {
        inner: ck::Element::HF,
    };
    #[classattr]
    const TA: Self = Self {
        inner: ck::Element::TA,
    };
    #[classattr]
    const W: Self = Self {
        inner: ck::Element::W,
    };
    #[classattr]
    const RE: Self = Self {
        inner: ck::Element::RE,
    };
    #[classattr]
    const OS: Self = Self {
        inner: ck::Element::OS,
    };
    #[classattr]
    const IR: Self = Self {
        inner: ck::Element::IR,
    };
    #[classattr]
    const PT: Self = Self {
        inner: ck::Element::PT,
    };
    #[classattr]
    const AU: Self = Self {
        inner: ck::Element::AU,
    };
    #[classattr]
    const HG: Self = Self {
        inner: ck::Element::HG,
    };
    #[classattr]
    const TL: Self = Self {
        inner: ck::Element::TL,
    };
    #[classattr]
    const PB: Self = Self {
        inner: ck::Element::PB,
    };
    #[classattr]
    const BI: Self = Self {
        inner: ck::Element::BI,
    };
    #[classattr]
    const PO: Self = Self {
        inner: ck::Element::PO,
    };
    #[classattr]
    const AT: Self = Self {
        inner: ck::Element::AT,
    };
    #[classattr]
    const RN: Self = Self {
        inner: ck::Element::RN,
    };
    #[classattr]
    const FR: Self = Self {
        inner: ck::Element::FR,
    };
    #[classattr]
    const RA: Self = Self {
        inner: ck::Element::RA,
    };
    #[classattr]
    const AC: Self = Self {
        inner: ck::Element::AC,
    };
    #[classattr]
    const TH: Self = Self {
        inner: ck::Element::TH,
    };
    #[classattr]
    const PA: Self = Self {
        inner: ck::Element::PA,
    };
    #[classattr]
    const U: Self = Self {
        inner: ck::Element::U,
    };
    #[classattr]
    const NP: Self = Self {
        inner: ck::Element::NP,
    };
    #[classattr]
    const PU: Self = Self {
        inner: ck::Element::PU,
    };
    #[classattr]
    const AM: Self = Self {
        inner: ck::Element::AM,
    };
    #[classattr]
    const CM: Self = Self {
        inner: ck::Element::CM,
    };
    #[classattr]
    const BK: Self = Self {
        inner: ck::Element::BK,
    };
    #[classattr]
    const CF: Self = Self {
        inner: ck::Element::CF,
    };
    #[classattr]
    const ES: Self = Self {
        inner: ck::Element::ES,
    };
    #[classattr]
    const FM: Self = Self {
        inner: ck::Element::FM,
    };
    #[classattr]
    const MD: Self = Self {
        inner: ck::Element::MD,
    };
    #[classattr]
    const NO: Self = Self {
        inner: ck::Element::NO,
    };
    #[classattr]
    const LR: Self = Self {
        inner: ck::Element::LR,
    };
    #[classattr]
    const RF: Self = Self {
        inner: ck::Element::RF,
    };
    #[classattr]
    const DB: Self = Self {
        inner: ck::Element::DB,
    };
    #[classattr]
    const SG: Self = Self {
        inner: ck::Element::SG,
    };
    #[classattr]
    const BH: Self = Self {
        inner: ck::Element::BH,
    };
    #[classattr]
    const HS: Self = Self {
        inner: ck::Element::HS,
    };
    #[classattr]
    const MT: Self = Self {
        inner: ck::Element::MT,
    };
    #[classattr]
    const DS: Self = Self {
        inner: ck::Element::DS,
    };
    #[classattr]
    const RG: Self = Self {
        inner: ck::Element::RG,
    };
    #[classattr]
    const CN: Self = Self {
        inner: ck::Element::CN,
    };
    #[classattr]
    const NH: Self = Self {
        inner: ck::Element::NH,
    };
    #[classattr]
    const FL: Self = Self {
        inner: ck::Element::FL,
    };
    #[classattr]
    const MC: Self = Self {
        inner: ck::Element::MC,
    };
    #[classattr]
    const LV: Self = Self {
        inner: ck::Element::LV,
    };
    #[classattr]
    const TS: Self = Self {
        inner: ck::Element::TS,
    };
    #[classattr]
    const OG: Self = Self {
        inner: ck::Element::OG,
    };
    #[classattr]
    const C: Self = Self {
        inner: ck::Element::C,
    };
    #[classattr]
    const O: Self = Self {
        inner: ck::Element::O,
    };
    #[classattr]
    const CL: Self = Self {
        inner: ck::Element::CL,
    };

    #[staticmethod]
    fn from_atomic_number(atomic_number: u8) -> Option<Self> {
        ck::Element::from_atomic_number(atomic_number).map(|inner| Self { inner })
    }
    #[staticmethod]
    fn from_symbol(symbol: &str) -> Option<Self> {
        ck::Element::from_symbol(symbol).map(|inner| Self { inner })
    }
    #[getter]
    fn atomic_number(&self) -> u8 {
        self.inner.atomic_number()
    }
    #[getter]
    fn symbol(&self) -> &'static str {
        self.inner.symbol()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen)]
    pub struct ElementInfo {
        inner: ck::ElementInfo,
    }
}

#[pymethods]
impl ElementInfo {
    #[getter]
    fn element(&self) -> Element {
        Element {
            inner: self.inner.element,
        }
    }
    #[getter]
    fn symbol(&self) -> &'static str {
        self.inner.symbol
    }
    #[getter]
    fn atomic_number(&self) -> u8 {
        self.inner.atomic_number
    }
    #[getter]
    fn period(&self) -> u8 {
        self.inner.period
    }
    #[getter]
    fn outer_electrons(&self) -> i32 {
        self.inner.outer_electrons
    }
    #[getter]
    fn valences(&self) -> Vec<i32> {
        self.inner.valences.to_vec()
    }
    #[getter]
    fn rb0(&self) -> f64 {
        self.inner.rb0
    }
    #[getter]
    fn atomic_weight(&self) -> f64 {
        self.inner.atomic_weight
    }
}

#[pyfunction]
pub fn element_info(element: &Element) -> ElementInfo {
    ElementInfo {
        inner: ck::element_info(element.inner),
    }
}
