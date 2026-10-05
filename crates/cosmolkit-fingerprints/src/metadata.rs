//! Shared source property-tree fingerprint metadata transport.
use serde_json::Value;
use std::fmt;
#[derive(Debug)]
pub enum FingerprintJsonError {
    Parse(serde_json::Error),
    Invalid(String),
    UnsupportedComponent {
        component: &'static str,
        source_type: String,
    },
}
impl fmt::Display for FingerprintJsonError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(e) => e.fmt(f),
            Self::Invalid(s) => f.write_str(s),
            Self::UnsupportedComponent {
                component,
                source_type,
            } => write!(
                f,
                "unsupported source generator component {component}: {source_type}"
            ),
        }
    }
}
impl std::error::Error for FingerprintJsonError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Parse(e) => Some(e),
            Self::Invalid(_) | Self::UnsupportedComponent { .. } => None,
        }
    }
}
pub(crate) fn parse_object(json: &str) -> Result<Value, FingerprintJsonError> {
    let value: Value = serde_json::from_str(json).map_err(FingerprintJsonError::Parse)?;
    if !value.is_object() {
        return Err(FingerprintJsonError::Invalid("expected JSON object".into()));
    }
    Ok(value)
}
pub(crate) fn common_arguments_string(
    count_simulation: bool,
    fp_size: u32,
    bits_per_feature: u32,
    include_chirality: bool,
) -> String {
    // RDKit source: FingerprintGenerator.cpp lines 50-54
    // RDKit❗✔️: std::string FingerprintArguments::commonArgumentsString() const {
    // RDKit❗✔️:   return "Common arguments : countSimulation=" +
    // RDKit❗✔️:          std::to_string(df_countSimulation) +
    // RDKit❗✔️:          " fpSize=" + std::to_string(d_fpSize) +
    // RDKit❗✔️:          " bitsPerFeature=" + std::to_string(d_numBitsPerFeature) +
    // RDKit❗✔️:          " includeChirality=" + std::to_string(df_includeChirality);
    // RDKit❗✔️: }
    format!(
        "Common arguments : countSimulation={} fpSize={} bitsPerFeature={} includeChirality={}",
        count_simulation as u8, fp_size, bits_per_feature, include_chirality as u8
    )
}

pub(crate) fn common_arguments_json(
    count_simulation: bool,
    fp_size: u32,
    bits_per_feature: u32,
    include_chirality: bool,
    count_bounds: &[u32],
) -> String {
    // RDKit source: FingerprintGenerator.cpp lines 58-71
    // RDKit❗✔️: void FingerprintArguments::toJSON(boost::property_tree::ptree &pt) const {
    // RDKit❗✔️:   pt.put("countSimulation", df_countSimulation);
    // RDKit❗✔️:   pt.put("fpSize", d_fpSize);
    // RDKit❗✔️:   pt.put("numBitsPerFeature", d_numBitsPerFeature);
    // RDKit❗✔️:   pt.put("includeChirality", df_includeChirality);
    // RDKit❗✔️:
    // RDKit❗✔️:   boost::property_tree::ptree countBoundsNode;
    // RDKit❗✔️:   for (const auto &bound : d_countBounds) {
    // RDKit❗✔️:     boost::property_tree::ptree boundNode;
    // RDKit❗✔️:     boundNode.put("", bound);
    // RDKit❗✔️:     countBoundsNode.push_back(std::make_pair("", boundNode));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   pt.add_child("countBounds", countBoundsNode);
    // RDKit❗✔️: }
    let count_bounds = count_bounds
        .iter()
        .map(|bound| format!("\"{bound}\""))
        .collect::<Vec<_>>()
        .join(",");
    // Boost property_tree stores every leaf as text, so write_json quotes
    // numeric and boolean values even though fromJSON reads typed values.
    format!(
        "{{\"countSimulation\":\"{}\",\"fpSize\":\"{}\",\"numBitsPerFeature\":\"{}\",\"includeChirality\":\"{}\",\"countBounds\":{}}}",
        count_simulation,
        fp_size,
        bits_per_feature,
        include_chirality,
        if count_bounds.is_empty() {
            "\"\"".to_string()
        } else {
            format!("[{count_bounds}]")
        }
    )
}

pub(crate) fn common_arguments_from_json(
    value: &Value,
    count_simulation: &mut bool,
    fp_size: &mut u32,
    bits_per_feature: &mut u32,
    include_chirality: &mut bool,
    count_bounds: &mut Vec<u32>,
) -> Result<(), FingerprintJsonError> {
    // RDKit source: FingerprintGenerator.cpp lines 73-90
    // RDKit❗✔️: void FingerprintArguments::fromJSON(const boost::property_tree::ptree &pt) {
    // RDKit❗✔️:   df_countSimulation = pt.get<bool>("countSimulation", df_countSimulation);
    // RDKit❗✔️:   d_fpSize = pt.get<std::uint32_t>("fpSize", d_fpSize);
    // RDKit❗✔️:   d_numBitsPerFeature =
    // RDKit❗✔️:       pt.get<std::uint32_t>("numBitsPerFeature", d_numBitsPerFeature);
    // RDKit❗✔️:   df_includeChirality = pt.get<bool>("includeChirality", df_includeChirality);
    // RDKit❗✔️:
    // RDKit❗✔️:   d_countBounds.clear();
    // RDKit❗✔️:   auto countBoundsNode = pt.get_child_optional("countBounds");
    // RDKit❗✔️:   if (countBoundsNode) {
    // RDKit❗✔️:     for (const auto &boundNode : *countBoundsNode) {
    // RDKit❗✔️:       d_countBounds.push_back(boundNode.second.get_value<std::uint32_t>());
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    let object = value
        .as_object()
        .ok_or_else(|| FingerprintJsonError::Invalid("expected JSON object".to_string()))?;

    if let Some(field) = object.get("countSimulation") {
        *count_simulation = json_value_as_bool("countSimulation", field)?;
    }
    if let Some(field) = object.get("fpSize") {
        *fp_size = json_value_as_u32("fpSize", field)?;
    }
    if let Some(field) = object.get("numBitsPerFeature") {
        *bits_per_feature = json_value_as_u32("numBitsPerFeature", field)?;
    }
    if let Some(field) = object.get("includeChirality") {
        *include_chirality = json_value_as_bool("includeChirality", field)?;
    }

    count_bounds.clear();
    if let Some(field) = object.get("countBounds") {
        if field.as_str() == Some("") {
            return Ok(());
        }
        let bounds = field.as_array().ok_or_else(|| {
            FingerprintJsonError::Invalid("countBounds must be an array".to_string())
        })?;
        for bound in bounds {
            count_bounds.push(json_value_as_u32("countBounds entry", bound)?);
        }
    }
    Ok(())
}
pub(crate) fn json_value_as_bool(name: &str, value: &Value) -> Result<bool, FingerprintJsonError> {
    if let Some(flag) = value.as_bool() {
        return Ok(flag);
    }
    if let Some(number) = value.as_u64() {
        return match number {
            0 => Ok(false),
            1 => Ok(true),
            _ => Err(FingerprintJsonError::Invalid(format!(
                "{name} must be a boolean"
            ))),
        };
    }
    if let Some(text) = value.as_str() {
        if let Ok(flag) = text.parse::<bool>() {
            return Ok(flag);
        }
        if let Ok(number) = text.parse::<u64>() {
            return match number {
                0 => Ok(false),
                1 => Ok(true),
                _ => Err(FingerprintJsonError::Invalid(format!(
                    "{name} must be a boolean"
                ))),
            };
        }
    }
    Err(FingerprintJsonError::Invalid(format!(
        "{name} must be a boolean"
    )))
}

pub(crate) fn json_value_as_u32(name: &str, value: &Value) -> Result<u32, FingerprintJsonError> {
    if let Some(number) = value.as_u64() {
        return u32::try_from(number).map_err(|_| {
            FingerprintJsonError::Invalid(format!("{name} must be a 32-bit integer"))
        });
    }
    if let Some(text) = value.as_str() {
        return text.parse::<u32>().map_err(|_| {
            FingerprintJsonError::Invalid(format!("{name} must be a 32-bit integer"))
        });
    }
    Err(FingerprintJsonError::Invalid(format!(
        "{name} must be a 32-bit integer"
    )))
}
