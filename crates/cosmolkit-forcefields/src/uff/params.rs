// Copyright (C) 2004-2024 Greg Landrum and other RDKit contributors.
// Copyright (C) 2004-2021 Greg Landrum and other RDKit contributors.
// This file contains behavior ported from RDKit and remains subject to its
// BSD license, included at the root of the pinned RDKit source tree.

use std::collections::{BTreeMap, HashMap};
use std::fmt;
use std::sync::{Arc, OnceLock, RwLock};

const DEFAULT_PARAM_DATA: &str = include_str!("default_params.tsv");

// RDKit❗✔️: constexpr double DEG2RAD = M_PI / 180.0;
pub(crate) const DEG2RAD: f64 = std::f64::consts::PI / 180.0;

// RDKit❗✔️: constexpr double RAD2DEG = 180.0 / M_PI;
pub(crate) const RAD2DEG: f64 = 180.0 / std::f64::consts::PI;

// RDKit❗✔️: const double lambda = 0.1332;
// RDKit❗✔️: const double G = 332.06;
// RDKit❗✔️: const double amideBondOrder = 1.41;
pub(crate) const PARAMS_LAMBDA: f64 = 0.1332;
pub(crate) const PARAMS_G: f64 = 332.06;
pub(crate) const PARAMS_AMIDE_BOND_ORDER: f64 = 1.41;

pub(crate) fn is_double_zero(value: f64) -> bool {
    // RDKit❗✔️: inline bool isDoubleZero(const double x) {
    // RDKit❗✔️:   return ((x < 1.0e-10) && (x > -1.0e-10));
    // RDKit❗✔️: }
    value < 1.0e-10 && value > -1.0e-10
}

pub(crate) fn clip_to_one(value: &mut f64) {
    // RDKit❗✔️: inline void clipToOne(double &x) { x = std::clamp(x, -1.0, 1.0); }
    // The ordered comparisons reproduce std::clamp's source ordering, including
    // leaving NaN unchanged when both comparisons are false.
    if *value < -1.0 {
        *value = -1.0;
    } else if *value > 1.0 {
        *value = 1.0;
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct UffBond {
    // RDKit❗✔️: double kb;
    pub(crate) kb: f64,
    // RDKit❗✔️: double r0;
    pub(crate) r0: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct UffAngle {
    // RDKit❗✔️: double ka;
    pub(crate) ka: f64,
    // RDKit❗✔️: double theta0;
    pub(crate) theta0: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct UffTor {
    // RDKit❗✔️: double V;
    pub(crate) v: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct UffInv {
    // RDKit❗✔️: double K;
    pub(crate) k: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct UffVdw {
    // RDKit❗✔️: double x_ij;
    pub(crate) x_ij: f64,
    // RDKit❗✔️: double D_ij;
    pub(crate) d_ij: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct AtomicParams {
    // RDKit❗✔️: double r1;
    pub(crate) r1: f64,
    // RDKit❗✔️: double theta0;
    pub(crate) theta0: f64,
    // RDKit❗✔️: double x1;
    pub(crate) x1: f64,
    // RDKit❗✔️: double D1;
    pub(crate) d1: f64,
    // RDKit❗✔️: double zeta;
    pub(crate) zeta: f64,
    // RDKit❗✔️: double Z1;
    pub(crate) z1: f64,
    // RDKit❗✔️: double V1;
    pub(crate) v1: f64,
    // RDKit❗✔️: double U1;
    pub(crate) u1: f64,
    // RDKit❗✔️: double GMP_Xi;
    pub(crate) gmp_xi: f64,
    // RDKit❗✔️: double GMP_Hardness;
    pub(crate) gmp_hardness: f64,
    // RDKit❗✔️: double GMP_Radius;
    pub(crate) gmp_radius: f64,
}

impl AtomicParams {
    fn fields(self) -> [f64; 11] {
        [
            self.r1,
            self.theta0,
            self.x1,
            self.d1,
            self.zeta,
            self.z1,
            self.v1,
            self.u1,
            self.gmp_xi,
            self.gmp_hardness,
            self.gmp_radius,
        ]
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum UffParamError {
    EmptyLine {
        line_number: usize,
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
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyLine { line_number } => {
                write!(formatter, "empty UFF parameter line {line_number}")
            }
            Self::MalformedLine {
                line_number,
                column_count,
            } => write!(
                formatter,
                "malformed UFF parameter line {line_number}: found {column_count} tokens before the required field"
            ),
            Self::ParseFloat {
                line_number,
                column_name,
                value,
            } => write!(
                formatter,
                "invalid UFF parameter float at line {line_number}, column {column_name}: {value}"
            ),
        }
    }
}

impl std::error::Error for UffParamError {}

#[derive(Debug)]
pub(crate) struct ParamCollection {
    params: BTreeMap<cosmolkit_model::PropertyText, AtomicParams>,
}

impl ParamCollection {
    pub(crate) fn get_params(param_data: &str) -> Result<Arc<Self>, UffParamError> {
        // RDKit❗✔️: const ParamCollection *ParamCollection::getParams(
        // RDKit❗✔️:     const std::string &paramData) {
        // RDKit❗✔️:   const ParamCollection *res = &(param_flyweight(paramData).get());
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        let registry = param_collection_registry();
        {
            let guard = registry
                .read()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            if let Some(existing) = guard.get(param_data) {
                return Ok(Arc::clone(existing));
            }
        }
        let mut guard = registry
            .write()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        // Another thread may have populated the exact flyweight key after the
        // shared lookup above.
        if let Some(existing) = guard.get(param_data) {
            return Ok(Arc::clone(existing));
        }

        let effective_data = if param_data.is_empty() {
            // RDKit❗✔️:   if (paramData.empty()) {
            // RDKit❗✔️:     paramData = defaultParamData;
            // RDKit❗✔️:   }
            default_param_data()
        } else {
            param_data
        };
        let collection = Arc::new(Self {
            params: parse_param_data(effective_data)?,
        });
        guard.insert(param_data.to_owned(), Arc::clone(&collection));
        Ok(collection)
    }

    pub(crate) fn get(&self, symbol: impl AsRef<[u8]>) -> Option<&AtomicParams> {
        // RDKit❗✔️: const AtomicParams *operator()(const std::string &symbol) const {
        // RDKit❗✔️:   std::map<std::string, AtomicParams>::const_iterator res;
        // RDKit❗✔️:   res = d_params.find(symbol);
        // RDKit❗✔️:   if (res != d_params.end()) {
        // RDKit❗✔️:     return &((*res).second);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return nullptr;
        // RDKit❗✔️: }
        // Behavior: std::map string equality and ordering compare counted
        // bytes. The canonical map key and borrowed query retain every byte;
        // opaque keys miss only if the exact byte key is absent in source data.
        // Complexity: one O(log P) borrowed tree lookup, no decoding or clone.
        self.params.get(symbol.as_ref())
    }

    pub(crate) fn len(&self) -> usize {
        self.params.len()
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.params.is_empty()
    }
}

fn param_collection_registry() -> &'static RwLock<HashMap<String, Arc<ParamCollection>>> {
    static REGISTRY: OnceLock<RwLock<HashMap<String, Arc<ParamCollection>>>> = OnceLock::new();
    REGISTRY.get_or_init(|| RwLock::new(HashMap::new()))
}

fn default_param_data() -> &'static str {
    // RDKit❗✔️: extern const std::string defaultParamData;
    // RDKit❗✔️: const std::string defaultParamData =
    DEFAULT_PARAM_DATA
}

fn source_records(param_data: &str) -> impl Iterator<Item = (usize, &str)> {
    // RDKit❗✔️: inline std::string getLine(std::istream *inStream) {
    // RDKit❗✔️:   std::string res;
    // RDKit❗✔️:   std::getline(*inStream, res);
    // RDKit❗✔️:   if (!res.empty() && (res.back() == '\r')) {
    // RDKit❗✔️:     res.resize(res.length() - 1);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    //
    // RDKit❗✔️: std::string inLine = RDKit::getLine(inStream);
    // RDKit❗✔️: while (!inStream.eof()) {
    // A failed final getline sets eofbit before the loop body, so only complete
    // newline-terminated records are yielded here.
    param_data
        .split_inclusive('\n')
        .enumerate()
        .filter_map(|(index, record)| {
            let line = record.strip_suffix('\n')?;
            let line = line.strip_suffix('\r').unwrap_or(line);
            Some((index + 1, line))
        })
}

fn parse_param_data(
    param_data: &str,
) -> Result<BTreeMap<cosmolkit_model::PropertyText, AtomicParams>, UffParamError> {
    // RDKit❗✔️: std::istringstream inStream(paramData);
    // RDKit❗✔️: while (!inStream.eof()) {
    // RDKit❗✔️:   if (inLine[0] != '#') {
    let mut params = BTreeMap::new();
    for (line_number, line) in source_records(param_data) {
        if line.is_empty() {
            // The C++ source indexes inLine[0] here; the Rust port reports the
            // malformed record instead of reproducing undefined behavior.
            return Err(UffParamError::EmptyLine { line_number });
        }
        if line.as_bytes()[0] == b'#' {
            continue;
        }

        let mut tokens = line.split('\t').filter(|token| !token.is_empty());
        // RDKit❗✔️:       boost::char_separator<char> tabSep("\t");
        // RDKit❗✔️:       tokenizer tokens(inLine, tabSep);
        // RDKit❗✔️:       tokenizer::iterator token = tokens.begin();
        let label = tokens.next().ok_or(UffParamError::MalformedLine {
            line_number,
            column_count: 0,
        })?;
        let mut token_count = 1;
        let params_for_label = AtomicParams {
            // RDKit❗✔️:       paramObj.r1 = boost::lexical_cast<double>(*token);
            r1: next_param_float(&mut tokens, &mut token_count, line_number, "r1")?,
            // RDKit❗✔️:       paramObj.theta0 = boost::lexical_cast<double>(*token);
            // RDKit❗✔️:       paramObj.theta0 = paramObj.theta0 * M_PI / 180.;
            theta0: {
                let theta0 =
                    next_param_float(&mut tokens, &mut token_count, line_number, "theta0")?;
                theta0 * std::f64::consts::PI / 180.0
            },
            // RDKit❗✔️:       paramObj.x1 = boost::lexical_cast<double>(*token);
            x1: next_param_float(&mut tokens, &mut token_count, line_number, "x1")?,
            // RDKit❗✔️:       paramObj.D1 = boost::lexical_cast<double>(*token);
            d1: next_param_float(&mut tokens, &mut token_count, line_number, "D1")?,
            // RDKit❗✔️:       paramObj.zeta = boost::lexical_cast<double>(*token);
            zeta: next_param_float(&mut tokens, &mut token_count, line_number, "zeta")?,
            // RDKit❗✔️:       paramObj.Z1 = boost::lexical_cast<double>(*token);
            z1: next_param_float(&mut tokens, &mut token_count, line_number, "Z1")?,
            // RDKit❗✔️:       paramObj.V1 = boost::lexical_cast<double>(*token);
            v1: next_param_float(&mut tokens, &mut token_count, line_number, "V1")?,
            // RDKit❗✔️:       paramObj.U1 = boost::lexical_cast<double>(*token);
            u1: next_param_float(&mut tokens, &mut token_count, line_number, "U1")?,
            // RDKit❗✔️:       paramObj.GMP_Xi = boost::lexical_cast<double>(*token);
            gmp_xi: next_param_float(&mut tokens, &mut token_count, line_number, "GMP_Xi")?,
            // RDKit❗✔️:       paramObj.GMP_Hardness = boost::lexical_cast<double>(*token);
            gmp_hardness: next_param_float(
                &mut tokens,
                &mut token_count,
                line_number,
                "GMP_Hardness",
            )?,
            // RDKit❗✔️:       paramObj.GMP_Radius = boost::lexical_cast<double>(*token);
            gmp_radius: next_param_float(&mut tokens, &mut token_count, line_number, "GMP_Radius")?,
        };
        // RDKit❗✔️:       d_params[label] = paramObj;
        // Unconsumed tokens intentionally follow the source and are ignored.
        params.insert(cosmolkit_model::PropertyText::from(label), params_for_label);
    }
    Ok(params)
}

fn next_param_float<'a>(
    tokens: &mut impl Iterator<Item = &'a str>,
    token_count: &mut usize,
    line_number: usize,
    column_name: &'static str,
) -> Result<f64, UffParamError> {
    let value = tokens.next().ok_or(UffParamError::MalformedLine {
        line_number,
        column_count: *token_count,
    })?;
    *token_count += 1;
    value.parse::<f64>().map_err(|_| UffParamError::ParseFloat {
        line_number,
        column_name,
        value: value.to_owned(),
    })
}

#[cfg(test)]
mod tests {
    use super::{
        AtomicParams, DEFAULT_PARAM_DATA, DEG2RAD, PARAMS_AMIDE_BOND_ORDER, PARAMS_G,
        PARAMS_LAMBDA, ParamCollection, RAD2DEG, UffParamError, clip_to_one, default_param_data,
        is_double_zero,
    };
    use std::sync::Arc;

    const FNV1A_OFFSET_BASIS: u64 = 0xcbf29ce484222325;
    const FNV1A_PRIME: u64 = 0x100000001b3;

    fn fnv1a_update(hash: &mut u64, bytes: &[u8]) {
        for byte in bytes {
            *hash ^= u64::from(*byte);
            *hash = hash.wrapping_mul(FNV1A_PRIME);
        }
    }

    fn parsed_fields(params: &AtomicParams) -> [u64; 11] {
        params.fields().map(f64::to_bits)
    }

    #[test]
    fn cf3d_u01_default_table_matches_fixed_all_field_fingerprint() {
        let collection = ParamCollection::get_params("").expect("pinned default UFF table");

        assert_eq!(default_param_data(), DEFAULT_PARAM_DATA);
        assert_eq!(DEFAULT_PARAM_DATA.len(), 7_552);
        assert!(DEFAULT_PARAM_DATA.ends_with('\n'));
        assert_eq!(collection.len(), 127);
        assert_eq!(collection.is_empty(), false);

        let mut raw_table_hash = FNV1A_OFFSET_BASIS;
        fnv1a_update(&mut raw_table_hash, DEFAULT_PARAM_DATA.as_bytes());
        assert_eq!(raw_table_hash, 0x6e4eb44649d1d902);

        // BTreeMap iteration is the source std::map order. Hash every key and
        // all eleven parsed fields, including the source theta conversion.
        let mut parsed_table_hash = FNV1A_OFFSET_BASIS;
        for (label, params) in &collection.params {
            fnv1a_update(&mut parsed_table_hash, label.as_bytes());
            fnv1a_update(&mut parsed_table_hash, &[0xff]);
            for field in params.fields() {
                fnv1a_update(&mut parsed_table_hash, &field.to_bits().to_le_bytes());
            }
        }
        assert_eq!(parsed_table_hash, 0xcc8c145b973fd1f8);

        let c3 = collection.get("C_3").expect("exact C_3 key");
        assert_eq!(
            parsed_fields(c3),
            [
                0.757_f64,
                109.47_f64 * std::f64::consts::PI / 180.0,
                3.851,
                0.105,
                12.73,
                1.912,
                2.119,
                2.0,
                5.343,
                5.063,
                0.759,
            ]
            .map(f64::to_bits)
        );
    }

    #[test]
    fn cf3d_u01_custom_rows_drop_empty_tokens_ignore_extras_and_strip_crlf() {
        let data = concat!(
            "# header\r\n",
            "O_3\t0.658\t104.51\t\t3.5\t0.06\t14.085\t2.3\t0.018\t2\t8.741\t6.682\t0.669\tignored\r\n",
        );
        let collection = ParamCollection::get_params(data).expect("source-shaped custom row");
        let oxygen = collection.get("O_3").expect("exact O_3 key");

        assert_eq!(collection.len(), 1);
        assert_eq!(
            parsed_fields(oxygen),
            [
                0.658_f64,
                104.51_f64 * std::f64::consts::PI / 180.0,
                3.5,
                0.06,
                14.085,
                2.3,
                0.018,
                2.0,
                8.741,
                6.682,
                0.669,
            ]
            .map(f64::to_bits)
        );
    }

    #[test]
    fn cf3d_u01_duplicate_labels_replace_and_lookup_is_exact() {
        let data = concat!(
            "# comment\n",
            "Q_test\t1\t\t30\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\n",
            "Q_test\t20\t90\t22\t23\t24\t25\t26\t27\t28\t29\t30\t31\n",
        );
        let collection = ParamCollection::get_params(data).expect("duplicate source labels");
        let replacement = collection.get("Q_test").expect("exact Q_test key");

        assert_eq!(collection.len(), 1);
        assert_eq!(
            parsed_fields(replacement),
            [
                20.0_f64,
                90.0_f64 * std::f64::consts::PI / 180.0,
                22.0,
                23.0,
                24.0,
                25.0,
                26.0,
                27.0,
                28.0,
                29.0,
                30.0,
            ]
            .map(f64::to_bits)
        );
        assert_eq!(collection.get("q_test"), None);
        assert_eq!(collection.get("missing"), None);
    }

    #[test]
    fn cf3d_u01_cache_uses_exact_original_key_and_empty_selects_defaults() {
        const ROW: &str = "CACHE_test\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\n";
        const ROW_WITH_COMMENT: &str =
            "# same parsed map\nCACHE_test\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\n";

        let first = ParamCollection::get_params(ROW).expect("valid custom row");
        let same_key = ParamCollection::get_params(ROW).expect("same exact flyweight key");
        let different_key =
            ParamCollection::get_params(ROW_WITH_COMMENT).expect("distinct exact flyweight key");
        assert!(Arc::ptr_eq(&first, &same_key));
        assert!(!Arc::ptr_eq(&first, &different_key));
        assert_eq!(first.params, different_key.params);

        let default_first = ParamCollection::get_params("").expect("empty key selects defaults");
        let default_same = ParamCollection::get_params("").expect("same empty key");
        assert!(Arc::ptr_eq(&default_first, &default_same));
        assert_eq!(default_first.len(), 127);
        assert_eq!(
            default_first.get("C_3"),
            Some(default_first.get("C_3").unwrap())
        );

        let explicit_default = ParamCollection::get_params(default_param_data())
            .expect("explicit table is a separate original key");
        assert!(!Arc::ptr_eq(&default_first, &explicit_default));
        assert_eq!(default_first.params, explicit_default.params);
    }

    #[test]
    fn cf3d_u01_unterminated_final_record_follows_stream_eof() {
        const FIRST: &str = "FIRST\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\n";
        const LAST_UNTERMINATED: &str = "LAST\t12\t13\t14\t15\t16\t17\t18\t19\t20\t21\t22";

        let complete = ParamCollection::get_params(concat!(
            "FIRST\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\n",
            "LAST\t12\t13\t14\t15\t16\t17\t18\t19\t20\t21\t22\n",
        ))
        .expect("both newline-terminated records");
        assert_eq!(complete.len(), 2);

        let with_unterminated_tail =
            ParamCollection::get_params(&format!("{FIRST}{LAST_UNTERMINATED}"))
                .expect("valid terminated record followed by unterminated source tail");
        assert_eq!(with_unterminated_tail.len(), 1);
        assert!(with_unterminated_tail.get("FIRST").is_some());
        assert_eq!(with_unterminated_tail.get("LAST"), None);

        let only_unterminated = ParamCollection::get_params(LAST_UNTERMINATED)
            .expect("source loop skips an unterminated first/final line");
        assert!(only_unterminated.is_empty());
    }

    #[test]
    fn cf3d_u01_malformed_rows_return_typed_context() {
        assert_eq!(
            ParamCollection::get_params("\n").unwrap_err(),
            UffParamError::EmptyLine { line_number: 1 }
        );
        assert_eq!(
            ParamCollection::get_params("C_3\t1\t2\n").unwrap_err(),
            UffParamError::MalformedLine {
                line_number: 1,
                column_count: 3,
            }
        );
        assert_eq!(
            ParamCollection::get_params("C_3\tbad\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\n").unwrap_err(),
            UffParamError::ParseFloat {
                line_number: 1,
                column_name: "r1",
                value: "bad".to_owned(),
            }
        );
        assert!(matches!(
            ParamCollection::get_params(" # first byte is not #\n").unwrap_err(),
            UffParamError::MalformedLine {
                line_number: 1,
                column_count: 1,
            }
        ));
    }

    #[test]
    fn cf3d_u01_source_helpers_preserve_strict_boundaries() {
        assert_eq!(DEG2RAD, std::f64::consts::PI / 180.0);
        assert_eq!(RAD2DEG, 180.0 / std::f64::consts::PI);
        assert_eq!(PARAMS_LAMBDA, 0.1332);
        assert_eq!(PARAMS_G, 332.06);
        assert_eq!(PARAMS_AMIDE_BOND_ORDER, 1.41);

        for value in [0.0, -0.0, 0.999e-10, -0.999e-10] {
            assert!(is_double_zero(value));
        }
        for value in [
            1.0e-10,
            -1.0e-10,
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
        ] {
            assert!(!is_double_zero(value));
        }

        for (input, expected) in [
            (-f64::INFINITY, -1.0),
            (-1.5, -1.0),
            (-1.0, -1.0),
            (0.25, 0.25),
            (1.0, 1.0),
            (1.5, 1.0),
            (f64::INFINITY, 1.0),
        ] {
            let mut value = input;
            clip_to_one(&mut value);
            assert_eq!(value, expected);
        }

        let mut nan = f64::from_bits(0x7ff8000000000042);
        clip_to_one(&mut nan);
        assert_eq!(nan.to_bits(), 0x7ff8000000000042);
    }
}
