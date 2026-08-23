//! Source-backed RDKit UFF parameter structures.

// Generated UFF tables preserve the exact source parameter literals.
#![allow(clippy::approx_constant)]

use std::collections::{BTreeMap, HashMap};
use std::fmt;
use std::num::ParseFloatError;
use std::sync::{Arc, OnceLock, RwLock};

include!(concat!(env!("OUT_DIR"), "/uff_defaults_generated.rs"));

// RDKit✔️✔️: constexpr double DEG2RAD = M_PI / 180.0;
pub const DEG2RAD: f64 = std::f64::consts::PI / 180.0;

// RDKit✔️✔️: constexpr double RAD2DEG = 180.0 / M_PI;
pub const RAD2DEG: f64 = 180.0 / std::f64::consts::PI;

// RDKit✔️✔️: const double lambda = 0.1332;  //!< scaling factor for rBO correction
pub const PARAMS_LAMBDA: f64 = 0.1332;

// RDKit✔️✔️: const double G = 332.06;       //!< bond force constant prefactor
pub const PARAMS_G: f64 = 332.06;

// RDKit✔️✔️: const double amideBondOrder =
// RDKit✔️✔️:     1.41;  //!< special case bond order for amide C-N bonds.
pub const PARAMS_AMIDE_BOND_ORDER: f64 = 1.41;

// RDKit✔️✔️: "#Atom	r1	theta0	x1	D1	zeta	Z1	Vi	Uj	"
// RDKit✔️✔️: "Xi	Hard	Radius\n"
pub const ATOMIC_PARAMS_TSV_COLUMNS: [&str; 12] = [
    "Atom", "r1", "theta0", "x1", "D1", "zeta", "Z1", "Vi", "Uj", "Xi", "Hard", "Radius",
];

pub fn is_double_zero(x: f64) -> bool {
    // RDKit✔️✔️: inline bool isDoubleZero(const double x) {
    // RDKit✔️✔️:   return ((x < 1.0e-10) && (x > -1.0e-10));
    // RDKit✔️✔️: }
    x < 1.0e-10 && x > -1.0e-10
}

pub fn clip_to_one(x: &mut f64) {
    // RDKit✔️✔️: inline void clipToOne(double &x) { x = std::clamp(x, -1.0, 1.0); }
    *x = x.clamp(-1.0, 1.0);
}

/// UFF parameters for bond stretching.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UffBond {
    // RDKit✔️✔️: double kb;
    pub kb: f64,
    // RDKit✔️✔️: double r0;
    pub r0: f64,
}

/// UFF parameters for angle bending.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UffAngle {
    // RDKit✔️✔️: double ka;
    pub ka: f64,
    // RDKit✔️✔️: double theta0;
    pub theta0: f64,
}

/// UFF parameters for torsions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UffTor {
    // RDKit✔️✔️: double V;
    pub v: f64,
}

/// UFF parameters for inversions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UffInv {
    // RDKit✔️✔️: double K;
    pub k: f64,
}

/// UFF parameters for van der Waals interactions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UffVdw {
    // RDKit✔️✔️: double x_ij;
    pub x_ij: f64,
    // RDKit✔️✔️: double D_ij;
    pub d_ij: f64,
}

/// Atomic parameters for the Universal Force Field.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AtomicParams {
    // RDKit✔️✔️: double r1;            //!<  valence bond radius
    pub r1: f64,
    // RDKit✔️✔️: double theta0;        //!< valence angle
    pub theta0: f64,
    // RDKit✔️✔️: double x1;            //!< vdW characteristic length
    pub x1: f64,
    // RDKit✔️✔️: double D1;            //!< vdW atomic energy
    pub d1: f64,
    // RDKit✔️✔️: double zeta;          //!< vdW scaling term
    pub zeta: f64,
    // RDKit✔️✔️: double Z1;            //!< effective charge
    pub z1: f64,
    // RDKit✔️✔️: double V1;            //!< sp3 torsional barrier parameter
    pub v1: f64,
    // RDKit✔️✔️: double U1;            //!< torsional contribution for sp2-sp3 bonds
    pub u1: f64,
    // RDKit✔️✔️: double GMP_Xi;        //!< GMP Electronegativity;
    pub gmp_xi: f64,
    // RDKit✔️✔️: double GMP_Hardness;  //!< GMP Hardness
    pub gmp_hardness: f64,
    // RDKit✔️✔️: double GMP_Radius;    //!< GMP Radius value
    pub gmp_radius: f64,
}

#[derive(Debug)]
pub struct ParamCollection {
    source_param_data: String,
    params: ParamStorage,
}

#[derive(Debug)]
enum ParamStorage {
    Static(&'static [(&'static str, AtomicParams)]),
    Owned(BTreeMap<String, AtomicParams>),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum UffParamError {
    MissingDefaultParamData,
    MalformedDefaultParamData {
        line: String,
    },
    MalformedLine {
        line_number: usize,
        column_count: usize,
    },
    ParseFloat {
        line_number: usize,
        column_name: &'static str,
        value: String,
    },
}

impl fmt::Display for UffParamError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingDefaultParamData => write!(f, "missing RDKit UFF default parameter data"),
            Self::MalformedDefaultParamData { line } => {
                write!(
                    f,
                    "malformed RDKit UFF default parameter string line: {line}"
                )
            }
            Self::MalformedLine {
                line_number,
                column_count,
            } => write!(
                f,
                "malformed UFF parameter line {line_number}: expected 12 tab-separated columns, found {column_count}"
            ),
            Self::ParseFloat {
                line_number,
                column_name,
                value,
            } => write!(
                f,
                "invalid UFF parameter float at line {line_number}, column {column_name}: {value}"
            ),
        }
    }
}

impl std::error::Error for UffParamError {}

impl ParamCollection {
    fn new_default_from_generated() -> Self {
        // RDKit✔️🔝: const std::string defaultParamData =
        // RDKit✔️🔝:     "#Atom	r1	theta0	x1	D1	zeta	Z1	Vi	Uj	"
        // RDKit✔️🔝:     "Xi	Hard	Radius\n"
        // RDKit✔️🔝:     "H_	0.354	180	2.886	0.044	12	0.712	0	0	"
        // Default RDKit table bytes are decoded at build time and emitted as
        // static Rust data. This preserves the source table and removes
        // runtime source scanning plus numeric parsing for the default path.
        Self {
            source_param_data: UFF_DEFAULT_PARAM_DATA.to_owned(),
            params: ParamStorage::Static(UFF_DEFAULT_ATOMIC_PARAMS),
        }
    }

    fn new_from_source(param_data: &str) -> Result<Self, UffParamError> {
        // RDKit✔️✔️: ParamCollection::ParamCollection(std::string paramData) {
        // RDKit✔️✔️:   if (paramData.empty()) {
        // RDKit✔️✔️:     paramData = defaultParamData;
        // RDKit✔️✔️:   }
        let source_param_data = if param_data.is_empty() {
            default_param_data()?.to_owned()
        } else {
            param_data.to_owned()
        };
        let params = parse_param_data(&source_param_data)?;
        Ok(Self {
            source_param_data,
            params: ParamStorage::Owned(params),
        })
    }

    pub fn get_params(param_data: &str) -> Result<Arc<Self>, UffParamError> {
        // RDKit✔️✔️: const ParamCollection *ParamCollection::getParams(
        // RDKit✔️✔️:     const std::string &paramData) {
        // RDKit✔️✔️:   const ParamCollection *res = &(param_flyweight(paramData).get());
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        if param_data.is_empty() {
            static DEFAULT: OnceLock<Arc<ParamCollection>> = OnceLock::new();
            return Ok(Arc::clone(
                DEFAULT.get_or_init(|| Arc::new(Self::new_default_from_generated())),
            ));
        }

        let registry = param_collection_registry();
        {
            let guard = registry
                .read()
                .expect("UFF ParamCollection registry lock poisoned");
            if let Some(existing) = guard.get(param_data) {
                return Ok(Arc::clone(existing));
            }
        }

        let mut guard = registry
            .write()
            .expect("UFF ParamCollection registry lock poisoned");
        if let Some(existing) = guard.get(param_data) {
            return Ok(Arc::clone(existing));
        }
        let collection = Arc::new(Self::new_from_source(param_data)?);
        guard.insert(param_data.to_owned(), Arc::clone(&collection));
        Ok(collection)
    }

    pub fn source_param_data(&self) -> &str {
        &self.source_param_data
    }

    pub fn len(&self) -> usize {
        match &self.params {
            ParamStorage::Static(params) => params.len(),
            ParamStorage::Owned(params) => params.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        match &self.params {
            ParamStorage::Static(params) => params.is_empty(),
            ParamStorage::Owned(params) => params.is_empty(),
        }
    }

    pub fn get(&self, symbol: &str) -> Option<&AtomicParams> {
        // RDKit✔️✔️: const AtomicParams *operator()(const std::string &symbol) const {
        // RDKit✔️✔️:   std::map<std::string, AtomicParams>::const_iterator res;
        // RDKit✔️✔️:   res = d_params.find(symbol);
        // RDKit✔️✔️:   if (res != d_params.end()) {
        // RDKit✔️✔️:     return &((*res).second);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return nullptr;
        // RDKit✔️✔️: }
        match &self.params {
            ParamStorage::Static(params) => params
                .binary_search_by(|(candidate, _)| candidate.cmp(&symbol))
                .ok()
                .map(|idx| &params[idx].1),
            ParamStorage::Owned(params) => params.get(symbol),
        }
    }
}

fn param_collection_registry() -> &'static RwLock<HashMap<String, Arc<ParamCollection>>> {
    static REGISTRY: OnceLock<RwLock<HashMap<String, Arc<ParamCollection>>>> = OnceLock::new();
    REGISTRY.get_or_init(|| RwLock::new(HashMap::new()))
}

pub fn default_param_data() -> Result<&'static str, UffParamError> {
    Ok(UFF_DEFAULT_PARAM_DATA)
}

fn parse_param_data(param_data: &str) -> Result<BTreeMap<String, AtomicParams>, UffParamError> {
    // RDKit✔️✔️:   std::istringstream inStream(paramData);
    // RDKit✔️✔️:   std::string inLine = RDKit::getLine(inStream);
    // RDKit✔️✔️:   while (!inStream.eof()) {
    // RDKit✔️✔️:     if (inLine[0] != '#') {
    let mut params = BTreeMap::new();
    for (line_idx, line) in param_data.lines().enumerate() {
        let line_number = line_idx + 1;
        if line.starts_with('#') {
            continue;
        }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() != ATOMIC_PARAMS_TSV_COLUMNS.len() {
            return Err(UffParamError::MalformedLine {
                line_number,
                column_count: columns.len(),
            });
        }
        // RDKit✔️✔️:       std::string label = *token;
        // RDKit✔️✔️:       ++token;
        let label = columns[0].to_owned();
        let param_obj = AtomicParams {
            // RDKit✔️✔️:       paramObj.r1 = boost::lexical_cast<double>(*token);
            r1: parse_param_float(line_number, "r1", columns[1])?,
            // RDKit✔️✔️:       paramObj.theta0 = boost::lexical_cast<double>(*token);
            // RDKit✔️✔️:       paramObj.theta0 = paramObj.theta0 * M_PI / 180.;
            theta0: parse_param_float(line_number, "theta0", columns[2])? * std::f64::consts::PI
                / 180.0,
            // RDKit✔️✔️:       paramObj.x1 = boost::lexical_cast<double>(*token);
            x1: parse_param_float(line_number, "x1", columns[3])?,
            // RDKit✔️✔️:       paramObj.D1 = boost::lexical_cast<double>(*token);
            d1: parse_param_float(line_number, "D1", columns[4])?,
            // RDKit✔️✔️:       paramObj.zeta = boost::lexical_cast<double>(*token);
            zeta: parse_param_float(line_number, "zeta", columns[5])?,
            // RDKit✔️✔️:       paramObj.Z1 = boost::lexical_cast<double>(*token);
            z1: parse_param_float(line_number, "Z1", columns[6])?,
            // RDKit✔️✔️:       paramObj.V1 = boost::lexical_cast<double>(*token);
            v1: parse_param_float(line_number, "Vi", columns[7])?,
            // RDKit✔️✔️:       paramObj.U1 = boost::lexical_cast<double>(*token);
            u1: parse_param_float(line_number, "Uj", columns[8])?,
            // RDKit✔️✔️:       paramObj.GMP_Xi = boost::lexical_cast<double>(*token);
            gmp_xi: parse_param_float(line_number, "Xi", columns[9])?,
            // RDKit✔️✔️:       paramObj.GMP_Hardness = boost::lexical_cast<double>(*token);
            gmp_hardness: parse_param_float(line_number, "Hard", columns[10])?,
            // RDKit✔️✔️:       paramObj.GMP_Radius = boost::lexical_cast<double>(*token);
            gmp_radius: parse_param_float(line_number, "Radius", columns[11])?,
        };
        // RDKit✔️✔️:       d_params[label] = paramObj;
        params.insert(label, param_obj);
    }
    Ok(params)
}

fn parse_param_float(
    line_number: usize,
    column_name: &'static str,
    value: &str,
) -> Result<f64, UffParamError> {
    value
        .parse::<f64>()
        .map_err(|_err: ParseFloatError| UffParamError::ParseFloat {
            line_number,
            column_name,
            value: value.to_owned(),
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1.0e-12;

    fn sample_atomic_params(marker: f64) -> AtomicParams {
        AtomicParams {
            r1: marker,
            theta0: marker + 1.0,
            x1: marker + 2.0,
            d1: marker + 3.0,
            zeta: marker + 4.0,
            z1: marker + 5.0,
            v1: marker + 6.0,
            u1: marker + 7.0,
            gmp_xi: marker + 8.0,
            gmp_hardness: marker + 9.0,
            gmp_radius: marker + 10.0,
        }
    }

    #[test]
    fn uff_params_structs_store_every_atomic_field_in_source_order() {
        let params = AtomicParams {
            r1: 1.0,
            theta0: 2.0,
            x1: 3.0,
            d1: 4.0,
            zeta: 5.0,
            z1: 6.0,
            v1: 7.0,
            u1: 8.0,
            gmp_xi: 9.0,
            gmp_hardness: 10.0,
            gmp_radius: 11.0,
        };

        assert_eq!(params.r1, 1.0);
        assert_eq!(params.theta0, 2.0);
        assert_eq!(params.x1, 3.0);
        assert_eq!(params.d1, 4.0);
        assert_eq!(params.zeta, 5.0);
        assert_eq!(params.z1, 6.0);
        assert_eq!(params.v1, 7.0);
        assert_eq!(params.u1, 8.0);
        assert_eq!(params.gmp_xi, 9.0);
        assert_eq!(params.gmp_hardness, 10.0);
        assert_eq!(params.gmp_radius, 11.0);
    }

    #[test]
    fn uff_params_structs_store_contrib_fields_in_source_order() {
        assert_eq!(UffBond { kb: 1.0, r0: 2.0 }, UffBond { kb: 1.0, r0: 2.0 });
        assert_eq!(
            UffAngle {
                ka: 3.0,
                theta0: 4.0,
            },
            UffAngle {
                ka: 3.0,
                theta0: 4.0,
            }
        );
        assert_eq!(UffTor { v: 5.0 }.v, 5.0);
        assert_eq!(UffInv { k: 6.0 }.k, 6.0);
        assert_eq!(
            UffVdw {
                x_ij: 7.0,
                d_ij: 8.0,
            },
            UffVdw {
                x_ij: 7.0,
                d_ij: 8.0,
            }
        );
    }

    #[test]
    fn uff_params_structs_preserve_source_constants_and_table_shape() {
        assert!((DEG2RAD - std::f64::consts::PI / 180.0).abs() < EPS);
        assert!((RAD2DEG - 180.0 / std::f64::consts::PI).abs() < EPS);
        assert_eq!(PARAMS_LAMBDA, 0.1332);
        assert_eq!(PARAMS_G, 332.06);
        assert_eq!(PARAMS_AMIDE_BOND_ORDER, 1.41);
        assert_eq!(
            ATOMIC_PARAMS_TSV_COLUMNS,
            [
                "Atom", "r1", "theta0", "x1", "D1", "zeta", "Z1", "Vi", "Uj", "Xi", "Hard",
                "Radius",
            ]
        );
    }

    #[test]
    fn uff_params_structs_cover_double_zero_boundaries() {
        assert!(is_double_zero(0.0));
        assert!(is_double_zero(9.999e-11));
        assert!(is_double_zero(-9.999e-11));
        assert!(!is_double_zero(1.0e-10));
        assert!(!is_double_zero(-1.0e-10));
    }

    #[test]
    fn uff_params_structs_cover_clip_to_one_boundaries() {
        let mut below = -1.5;
        clip_to_one(&mut below);
        assert_eq!(below, -1.0);

        let mut inside = 0.25;
        clip_to_one(&mut inside);
        assert_eq!(inside, 0.25);

        let mut above = 1.5;
        clip_to_one(&mut above);
        assert_eq!(above, 1.0);
    }

    #[test]
    fn uff_param_collection_get_params_reuses_same_param_data_instance() {
        let first = ParamCollection::get_params("A\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\n")
            .expect("valid UFF params");
        let second = ParamCollection::get_params("A\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\n")
            .expect("valid UFF params");

        assert!(Arc::ptr_eq(&first, &second));
        assert_eq!(
            first.source_param_data(),
            "A\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\n"
        );
        assert!(!first.is_empty());
        assert_eq!(first.len(), 1);
    }

    #[test]
    fn uff_param_collection_get_params_distinguishes_different_param_data() {
        let first = ParamCollection::get_params("A\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\n")
            .expect("valid UFF params");
        let second = ParamCollection::get_params("B\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\n")
            .expect("valid UFF params");

        assert!(!Arc::ptr_eq(&first, &second));
        assert_eq!(
            first.source_param_data(),
            "A\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\n"
        );
        assert_eq!(
            second.source_param_data(),
            "B\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\n"
        );
    }

    #[test]
    fn uff_param_collection_get_params_caches_empty_param_data_key() {
        let first = ParamCollection::get_params("").expect("valid default UFF params");
        let second = ParamCollection::get_params("").expect("valid default UFF params");

        assert!(Arc::ptr_eq(&first, &second));
        assert!(first.source_param_data().starts_with("#Atom\tr1\ttheta0"));
        assert!(first.get("C_3").is_some());
    }

    #[test]
    fn uff_param_collection_operator_returns_existing_atomic_params() {
        let expected = sample_atomic_params(10.0);
        let mut params = BTreeMap::new();
        params.insert("C_3".to_string(), expected);
        let collection = ParamCollection {
            source_param_data: String::new(),
            params: ParamStorage::Owned(params),
        };

        assert_eq!(collection.get("C_3"), Some(&expected));
    }

    #[test]
    fn uff_param_collection_operator_returns_none_for_missing_symbol() {
        let mut params = BTreeMap::new();
        params.insert("C_3".to_string(), sample_atomic_params(10.0));
        let collection = ParamCollection {
            source_param_data: String::new(),
            params: ParamStorage::Owned(params),
        };

        assert_eq!(collection.get("N_3"), None);
    }

    #[test]
    fn uff_param_collection_default_loader_parses_source_table() {
        let collection = ParamCollection::get_params("").expect("valid default UFF params");
        let c3 = collection.get("C_3").expect("C_3 default params");

        assert!(collection.len() > 100);
        assert_eq!(c3.r1, 0.757);
        assert!((c3.theta0 - 109.47 * DEG2RAD).abs() < EPS);
        assert_eq!(c3.x1, 3.851);
        assert_eq!(c3.d1, 0.105);
        assert_eq!(c3.zeta, 12.73);
        assert_eq!(c3.z1, 1.912);
        assert_eq!(c3.v1, 2.119);
        assert_eq!(c3.u1, 2.0);
        assert_eq!(c3.gmp_xi, 5.343);
        assert_eq!(c3.gmp_hardness, 5.063);
        assert_eq!(c3.gmp_radius, 0.759);
        assert!(collection.get("Lw6+3").is_some());
    }

    #[test]
    fn uff_param_collection_preserves_rdkit_theta0_conversion_evaluation_order() {
        // RDKit✔️✔️:       paramObj.theta0 = boost::lexical_cast<double>(*token);
        // RDKit✔️✔️:       paramObj.theta0 = paramObj.theta0 * M_PI / 180.;
        let theta0_degrees = std::hint::black_box(104.51);
        let rdkit_theta0 = theta0_degrees * std::f64::consts::PI / 180.0;
        let previously_reassociated_theta0 = theta0_degrees * std::hint::black_box(DEG2RAD);
        assert_ne!(
            rdkit_theta0.to_bits(),
            previously_reassociated_theta0.to_bits()
        );

        let defaults = ParamCollection::get_params("").expect("valid default UFF params");
        let default_o3 = defaults.get("O_3").expect("O_3 default params");
        assert_eq!(default_o3.theta0.to_bits(), rdkit_theta0.to_bits());

        let custom_data =
            "O_3\t0.658\t104.51\t3.5\t0.06\t14.085\t2.3\t0.018\t2\t8.741\t6.682\t0.669\n";
        let custom = ParamCollection::get_params(custom_data).expect("valid custom UFF params");
        let custom_o3 = custom.get("O_3").expect("O_3 custom params");
        assert_eq!(custom_o3.theta0.to_bits(), rdkit_theta0.to_bits());
    }

    #[test]
    fn uff_param_collection_loader_skips_comment_lines_and_overwrites_duplicate_labels() {
        let data = "#Atom\tr1\ttheta0\tx1\tD1\tzeta\tZ1\tVi\tUj\tXi\tHard\tRadius\n\
                   C_3\t1\t180\t3\t4\t5\t6\t7\t8\t9\t10\t11\n\
                   C_3\t2\t90\t4\t5\t6\t7\t8\t9\t10\t11\t12\n";
        let collection = ParamCollection::get_params(data).expect("valid duplicate UFF params");
        let c3 = collection.get("C_3").expect("C_3 params");

        assert_eq!(collection.len(), 1);
        assert_eq!(c3.r1, 2.0);
        assert!((c3.theta0 - 90.0 * DEG2RAD).abs() < EPS);
        assert_eq!(c3.gmp_radius, 12.0);
    }

    #[test]
    fn uff_param_collection_loader_reports_malformed_column_count() {
        let err = ParamCollection::get_params("C_3\t1\t2\n").expect_err("malformed column count");

        assert_eq!(
            err,
            UffParamError::MalformedLine {
                line_number: 1,
                column_count: 3,
            }
        );
    }

    #[test]
    fn uff_param_collection_loader_reports_float_parse_error() {
        let err = ParamCollection::get_params("C_3\tbad\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\n")
            .expect_err("invalid float");

        assert_eq!(
            err,
            UffParamError::ParseFloat {
                line_number: 1,
                column_name: "r1",
                value: "bad".to_string(),
            }
        );
    }
}
