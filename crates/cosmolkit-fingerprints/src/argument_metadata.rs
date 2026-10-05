//! Source argument metadata; immutable updates of canonical input values.
use crate::metadata::{
    bool_or, common_arguments_from_json, common_arguments_json, parse_object, u32_or,
};
use crate::{AtomPairParams, FingerprintJsonError, MorganParams, TopologicalTorsionParams};

impl MorganParams {
    #[must_use]
    pub fn info_string(&self) -> String {
        // RDKit❗✔️: std::string MorganArguments::infoString() const {
        // RDKit❗✔️:   return "MorganArguments onlyNonzeroInvariants=" +
        // RDKit❗✔️:          std::to_string(df_onlyNonzeroInvariants) +
        // RDKit❗✔️:          " radius=" + std::to_string(d_radius);
        // RDKit❗✔️: }
        // Behavior: argument fields only; provider settings are independent.
        // Complexity: one fixed-length format plus integer digits, no graph work.
        format!(
            "MorganArguments onlyNonzeroInvariants={} radius={}",
            self.only_nonzero_invariants as u8, self.radius
        )
    }
    #[must_use]
    pub fn to_json(&self) -> String {
        // RDKit❗✔️: void MorganArguments::toJSON(boost::property_tree::ptree &pt) const {
        // RDKit❗✔️:   pt.put("type", "MorganArguments");
        // RDKit❗✔️:   pt.put("onlyNonzeroInvariants", df_onlyNonzeroInvariants);
        // RDKit❗✔️:   pt.put("radius", d_radius);
        // RDKit❗✔️:   FingerprintArguments::toJSON(pt);
        // RDKit❗✔️: }
        // Behavior: source intentionally omits redundant environments and useBondTypes;
        // ring membership belongs to the independent invariant provider.
        // Complexity: O(countBounds) text serialization, sole shared common owner.
        let common = common_arguments_json(
            self.count_simulation,
            self.fp_size,
            self.bits_per_feature,
            self.include_chirality,
            &self.count_bounds,
        );
        morgan_arguments_json(self.radius, self.only_nonzero_invariants, &common)
    }
    /// Apply source JSON fields to a new immutable configuration value.
    pub fn with_json(&self, json: &str) -> Result<Self, FingerprintJsonError> {
        // RDKit❗✔️: void MorganArguments::fromJSON(const boost::property_tree::ptree &pt) {
        // RDKit❗✔️:   d_radius = pt.get<std::uint32_t>("radius", d_radius);
        // RDKit❗✔️:   df_onlyNonzeroInvariants =
        // RDKit❗✔️:       pt.get<bool>("onlyNonzeroInvariants", df_onlyNonzeroInvariants);
        // RDKit❗✔️:   FingerprintArguments::fromJSON(pt);
        // RDKit❗✔️: }
        // Modern d892 source fromJSON explicitly returns for whitespace-only JSON.
        // Behavior: source update order, absent bounds clears, no factory precondition;
        // immutable projection preserves the input on successful and failed updates.
        // Complexity: fixed scalar copy plus O(JSON+new countBounds), no old-vector or molecule clone.
        if json.trim().is_empty() {
            return Ok(self.clone());
        }
        self.with_json_value(&parse_object(json)?)
    }
    pub(crate) fn with_json_value(
        &self,
        value: &crate::metadata::SourceNode,
    ) -> Result<Self, FingerprintJsonError> {
        // RDKit❗✔️:   d_radius = pt.get<std::uint32_t>("radius", d_radius);
        // RDKit❗✔️:   df_onlyNonzeroInvariants =
        // RDKit❗✔️:       pt.get<bool>("onlyNonzeroInvariants", df_onlyNonzeroInvariants);
        // RDKit❗✔️:   FingerprintArguments::fromJSON(pt);
        // Behavior: same source node, no JSON normalization/reparse in reusable restoration.
        // Complexity: fixed fields plus O(new bounds); no old bounds clone.
        // Every nonempty source update clears countBounds. Copy only retained
        // scalar state; do not clone the old vector just to immediately clear it.
        let mut result = Self {
            radius: self.radius,
            include_chirality: self.include_chirality,
            use_bond_types: self.use_bond_types,
            include_ring_membership: self.include_ring_membership,
            only_nonzero_invariants: self.only_nonzero_invariants,
            include_redundant_environments: self.include_redundant_environments,
            fp_size: self.fp_size,
            count_simulation: self.count_simulation,
            bits_per_feature: self.bits_per_feature,
            count_bounds: Vec::new(),
        };
        result.radius = u32_or(&value, "radius", result.radius);
        result.only_nonzero_invariants = bool_or(
            &value,
            "onlyNonzeroInvariants",
            result.only_nonzero_invariants,
        );
        common_arguments_from_json(
            &value,
            &mut result.count_simulation,
            &mut result.fp_size,
            &mut result.bits_per_feature,
            &mut result.include_chirality,
            &mut result.count_bounds,
        )?;
        Ok(result)
    }
}
impl AtomPairParams {
    #[must_use]
    pub fn info_string(&self) -> String {
        // RDKit❗✔️: std::string AtomPairArguments::infoString() const {
        // RDKit❗✔️:   return "AtomPairArguments use2D=" + std::to_string(df_use2D) +
        // RDKit❗✔️:          " minDistance=" + std::to_string(d_minDistance) +
        // RDKit❗✔️:          " maxDistance=" + std::to_string(d_maxDistance);
        // RDKit❗✔️: }
        // Behavior: exact derived argument fields; complexity: fixed-field formatting.
        format!(
            "AtomPairArguments use2D={} minDistance={} maxDistance={}",
            self.use_2d as u8, self.min_distance, self.max_distance
        )
    }
    #[must_use]
    pub fn to_json(&self) -> String {
        // RDKit❗✔️: void AtomPairArguments::toJSON(boost::property_tree::ptree &pt) const {
        // RDKit❗✔️:   pt.put("type", "AtomPairArguments");
        // RDKit❗✔️:   pt.put("use2D", df_use2D);
        // RDKit❗✔️:   pt.put("minDistance", d_minDistance);
        // RDKit❗✔️:   pt.put("maxDistance", d_maxDistance);
        // RDKit❗✔️:   FingerprintArguments::toJSON(pt);
        // RDKit❗✔️: }
        // Behavior: quoted Boost leaf values; complexity: O(countBounds) output.
        let common = common_arguments_json(
            self.count_simulation,
            self.fp_size,
            self.bits_per_feature,
            self.include_chirality,
            &self.count_bounds,
        );
        format!(
            "{{\"type\":\"AtomPairArguments\",\"use2D\":\"{}\",\"minDistance\":\"{}\",\"maxDistance\":\"{}\",{}}}",
            self.use_2d,
            self.min_distance,
            self.max_distance,
            &common[1..common.len() - 1]
        )
    }
    pub fn with_json(&self, json: &str) -> Result<Self, FingerprintJsonError> {
        // RDKit❗✔️: void AtomPairArguments::fromJSON(const boost::property_tree::ptree &pt) {
        // RDKit❗✔️:   df_use2D = pt.get<bool>("use2D", df_use2D);
        // RDKit❗✔️:   d_minDistance = pt.get<unsigned int>("minDistance", d_minDistance);
        // RDKit❗✔️:   d_maxDistance = pt.get<unsigned int>("maxDistance", d_maxDistance);
        // RDKit❗✔️:   FingerprintArguments::fromJSON(pt);
        // RDKit❗✔️: }
        // Behavior: update order does not repeat constructor distance validation.
        // Complexity: fixed scalar copy and source-linear JSON/new-bounds traversal.
        if json.trim().is_empty() {
            return Ok(self.clone());
        }
        self.with_json_value(&parse_object(json)?)
    }
    pub(crate) fn with_json_value(
        &self,
        value: &crate::metadata::SourceNode,
    ) -> Result<Self, FingerprintJsonError> {
        // RDKit❗✔️:   df_use2D = pt.get<bool>("use2D", df_use2D);
        // RDKit❗✔️:   d_minDistance = pt.get<unsigned int>("minDistance", d_minDistance);
        // RDKit❗✔️:   d_maxDistance = pt.get<unsigned int>("maxDistance", d_maxDistance);
        // RDKit❗✔️:   FingerprintArguments::fromJSON(pt);
        // Behavior: same source node; complexity: fixed fields plus O(new bounds).
        // Every nonempty source update clears countBounds. Copy only retained
        // scalar state; do not clone the old vector just to immediately clear it.
        let mut result = Self {
            min_distance: self.min_distance,
            max_distance: self.max_distance,
            include_chirality: self.include_chirality,
            use_2d: self.use_2d,
            count_simulation: self.count_simulation,
            fp_size: self.fp_size,
            bits_per_feature: self.bits_per_feature,
            count_bounds: Vec::new(),
        };
        result.use_2d = bool_or(&value, "use2D", result.use_2d);
        result.min_distance = u32_or(&value, "minDistance", result.min_distance);
        result.max_distance = u32_or(&value, "maxDistance", result.max_distance);
        common_arguments_from_json(
            &value,
            &mut result.count_simulation,
            &mut result.fp_size,
            &mut result.bits_per_feature,
            &mut result.include_chirality,
            &mut result.count_bounds,
        )?;
        Ok(result)
    }
}
impl TopologicalTorsionParams {
    pub fn with_json(&self, json: &str) -> Result<Self, FingerprintJsonError> {
        // RDKit❗✔️: void TopologicalTorsionArguments::fromJSON(
        // RDKit❗✔️:     const boost::property_tree::ptree &pt) {
        // RDKit❗✔️:   d_torsionAtomCount = pt.get<uint32_t>("torsionAtomCount", d_torsionAtomCount);
        // RDKit❗✔️:   df_onlyShortestPaths =
        // RDKit❗✔️:       pt.get<bool>("onlyShortestPaths", df_onlyShortestPaths);
        // RDKit❗✔️:   FingerprintArguments::fromJSON(pt);
        // RDKit❗✔️: }
        // Behavior: immutable public projection delegates all source field updates.
        // Complexity: fixed scalar copy and existing source-linear parser; old bounds are not cloned.
        if json.trim().is_empty() {
            return Ok(self.clone());
        }
        let mut result = Self {
            torsion_atom_count: self.torsion_atom_count,
            only_shortest_paths: self.only_shortest_paths,
            include_chirality: self.include_chirality,
            count_simulation: self.count_simulation,
            fp_size: self.fp_size,
            bits_per_feature: self.bits_per_feature,
            count_bounds: Vec::new(),
        };
        result.from_json(json)?;
        Ok(result)
    }
}

/// Shared source Morgan argument fields, consuming no configuration snapshots.
pub(crate) fn morgan_arguments_json(
    radius: u32,
    only_nonzero_invariants: bool,
    common: &str,
) -> String {
    // RDKit❗✔️:   pt.put("type", "MorganArguments");
    // RDKit❗✔️:   pt.put("onlyNonzeroInvariants", df_onlyNonzeroInvariants);
    // RDKit❗✔️:   pt.put("radius", d_radius);
    // RDKit❗✔️:   FingerprintArguments::toJSON(pt);
    // Behavior: sole argument/header composition, source-quoted Boost leaves.
    // Complexity: source-linear serialized text output; no countBounds clone.
    format!(
        "{{\"type\":\"MorganArguments\",\"onlyNonzeroInvariants\":\"{}\",\"radius\":\"{}\",{}}}",
        only_nonzero_invariants,
        radius,
        &common[1..common.len() - 1]
    )
}
