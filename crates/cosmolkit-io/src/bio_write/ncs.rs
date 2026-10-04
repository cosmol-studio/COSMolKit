//! Private mmCIF `_struct_ncs_oper` emitter (BIO-NCS-WRITE1-28).
//!
//! Not enabled in the coordinate-only composer; the public writer still
//! does NOT emit NCS. Owns only the source-backed row construction.

use super::BioMmcifWriteError;
use crate::cif::{CifBlock, CifLoop, format_cif_f64};
use cosmolkit_bio::{BioStructureData, BioTransform};

pub(super) fn write_ncs_oper(
    data: &BioStructureData,
    block: &mut CifBlock,
) -> Result<(), BioMmcifWriteError> {
    // BEGIN GEMMI CPP FUNCTION gemmi::write_ncs_oper (to_mmcif.cpp:326-350)
    // Gemmi✔️✔️: void write_ncs_oper(const Structure& st, cif::Block& block) {
    // Gemmi✔️✔️:   // _struct_ncs_oper (MTRIX)
    // Gemmi✔️✔️:   if (st.ncs.empty())
    // Gemmi✔️✔️:     return;
    // Gemmi✔️✔️:   cif::Loop& ncs_oper = block.init_mmcif_loop("_struct_ncs_oper.",
    // Gemmi✔️✔️:       {"id", "code",
    // Gemmi✔️✔️:        "matrix[1][1]", "matrix[1][2]", "matrix[1][3]", "vector[1]",
    // Gemmi✔️✔️:        "matrix[2][1]", "matrix[2][2]", "matrix[2][3]", "vector[2]",
    // Gemmi✔️✔️:        "matrix[3][1]", "matrix[3][2]", "matrix[3][3]", "vector[3]"});
    // Gemmi✔️✔️:   auto add_op = [&ncs_oper](const NcsOp& op) {
    // Gemmi✔️❌:     ncs_oper.values.emplace_back(op.id);
    // Gemmi✔️❌:     ncs_oper.values.emplace_back(op.given ? "given" : "generate");
    // Gemmi✔️✔️:     for (int i = 0; i < 3; ++i) {
    // Gemmi✔️✔️:       for (int j = 0; j < 3; ++j)
    // Gemmi❗❌:         ncs_oper.values.emplace_back(to_str(op.tr.mat[i][j]));
    // Gemmi❗❌:       ncs_oper.values.emplace_back(to_str(op.tr.vec.at(i)));
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:   };
    // Gemmi✔️✔️:   auto identity = st.info.find("_struct_ncs_oper.id");
    // Gemmi✔️✔️:   if (identity != st.info.end() &&
    // Gemmi✔️✔️:       !in_vector_f([&](const NcsOp& op) { return op.id == identity->second; }, st.ncs))
    // Gemmi✔️❌:     add_op(NcsOp{identity->second, true, {}});
    // Gemmi✔️✔️:   for (const NcsOp& op : st.ncs)
    // Gemmi✔️✔️:     add_op(op);
    // Gemmi✔️✔️: }
    // END GEMMI CPP FUNCTION
    //
    // Helper bodies consulted (cross-file rule): in_vector_f (util.hpp:246-249,
    // one find_if scan); NcsOp (unitcell.hpp:121-126); Transform/Mat33/Vec3
    // defaults (math.hpp:59-62 Vec3_() zeroes, :127-131 identity Mat33) — the
    // identity insert reuses BioTransform::identity(); to_str(double)
    // (sprintf.hpp:36-40, "%.9g") delegated to crate::cif::format_cif_f64;
    // init_mmcif_loop/ensure/find/setup (cifdoc.hpp:506-509/:69-73/:1031-1050/
    // :941-966) delegated to CifBlock::init_mmcif_loop.
    //
    // Behavior review: empty operator table returns success BEFORE any block
    // setup or identity inspection; otherwise ONE canonical init_mmcif_loop
    // with the exact 14 suffixes replaces the category (existing loop cleared
    // in place at its source position / first-pair position retained).
    // The optional identity id is inserted FIRST iff present in source state
    // AND no operator carries that exact id (given=true, identity transform);
    // every original operator then follows in original order, duplicate ids
    // retained. Each row appends id RAW (never re-quoted or normalized),
    // code given/generate, then per row three matrix cells followed by the
    // vector component. Structural CIF errors (mixed-category loops)
    // propagate unchanged through BioMmcifWriteError::Cif BEFORE mutation.
    //
    // State review: the structure and operator list are only borrowed; no
    // collect/sort/deduplicate/clone of rows or of the whole data value;
    // source-state maps are read via the existing ncs_oper_identity_id view.
    //
    // Collection/storage cost review (axes kept separate): control flow and
    // the 14 amortized-O(1) pushes per row are complexity-equivalent (no
    // temporary per-row Vec). String creation is NOT claimed equivalent:
    // the Rust id clone, identity id.to_string() and given/generate
    // .to_string() are always heap-owning Strings, while the inspected
    // libstdc++ ABI stores std::string values ≤15 bytes inline (SSO) —
    // those lines carry second-axis ❌: for the inspected ABI and these
    // short nonempty values, the additional Rust allocation is known from
    // code inspection, not merely unresolved pending measurement. This is
    // ABI-scoped, not a wall-clock benchmark claim. Numeric cells delegate to
    // format_cif_f64, whose own in-body markers (behavior ❗ sampled-oracle
    // profile, Rust-only allocation measurement) remain the authority; the
    // frozen fixed literals verify only these finite inputs and upgrade
    // nothing about arbitrary floating output.
    //
    // Delegated numeric cost: format_cif_f64 carries the existing qualified
    // Gemmi to_str(double) profile debt (in-body markers in cif.rs); this
    // writer neither upgrades nor blanket-claims that formatter's
    // allocation/performance behavior.
    if data.ncs_operators().is_empty() {
        return Ok(());
    }
    let ncs_oper = block.init_mmcif_loop(
        "_struct_ncs_oper.",
        &[
            "id",
            "code",
            "matrix[1][1]",
            "matrix[1][2]",
            "matrix[1][3]",
            "vector[1]",
            "matrix[2][1]",
            "matrix[2][2]",
            "matrix[2][3]",
            "vector[2]",
            "matrix[3][1]",
            "matrix[3][2]",
            "matrix[3][3]",
            "vector[3]",
        ],
    )?;
    if let Some(identity_id) = data.ncs_oper_identity_id()
        // BEGIN GEMMI CPP FUNCTION gemmi::in_vector_f (util.hpp:246-249)
        // Gemmi✔️✔️: template <typename F, typename T>
        // Gemmi✔️✔️: bool in_vector_f(F f, const std::vector<T>& v) {
        // Gemmi✔️✔️:   return std::find_if(v.begin(), v.end(), f) != v.end();
        // Gemmi✔️✔️: }
        // END GEMMI CPP FUNCTION
        // Behavior review: the .any() below is this exact predicate
        // reproduced locally — one find_if-style linear scan with
        // first-match short-circuit over the borrowed operator slice,
        // called ONCE only when an identity id exists.
        // Complexity review: O(n) single pass over the operator table,
        // constant borrowed closure, no allocation.
        && !data
            .ncs_operators()
            .iter()
            .any(|operator| operator.id == identity_id)
    {
        append_op(
            ncs_oper,
            identity_id.to_string(),
            true,
            &BioTransform::identity(),
        );
    }
    for operator in data.ncs_operators() {
        append_op(
            ncs_oper,
            operator.id.clone(),
            operator.given,
            &operator.transform,
        );
    }
    Ok(())
}

/// One source-shaped `add_op` row: id, code, then per row the three matrix
/// cells followed by the vector component.
fn append_op(ncs_oper: &mut CifLoop, id: String, given: bool, transform: &BioTransform) {
    // BEGIN GEMMI CPP FUNCTION add_op lambda (to_mmcif.cpp:335-343)
    // Gemmi✔️✔️:   auto add_op = [&ncs_oper](const NcsOp& op) {
    // Gemmi✔️❌:     ncs_oper.values.emplace_back(op.id);
    // Gemmi✔️❌:     ncs_oper.values.emplace_back(op.given ? "given" : "generate");
    // Gemmi✔️✔️:     for (int i = 0; i < 3; ++i) {
    // Gemmi✔️✔️:       for (int j = 0; j < 3; ++j)
    // Gemmi❗❌:         ncs_oper.values.emplace_back(to_str(op.tr.mat[i][j]));
    // Gemmi❗❌:       ncs_oper.values.emplace_back(to_str(op.tr.vec.at(i)));
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:   };
    // END GEMMI CPP FUNCTION
    // Behavior review: this function is the complete lambda body — id
    // appended RAW first, code given/generate second, then for each of the
    // three rows the three matrix cells followed by the vector component,
    // in the source's exact order.
    // Cost review (separate axes): the control flow and 14 amortized-O(1)
    // pushes are complexity-equivalent; the id and code cells are newly
    // created Rust heap Strings where the inspected libstdc++ std::string
    // would use SSO inline storage for these short nonempty values (known
    // extra allocation, second-axis ❌ for the inspected ABI), and the numeric
    // cells delegate to format_cif_f64's existing ❗❌ profile — neither a
    // generic numeric behavior upgrade nor a blanket cost-equivalence claim.
    ncs_oper.append_raw_value(id);
    ncs_oper.append_raw_value(if given { "given" } else { "generate" }.to_string());
    for row in 0..3 {
        for column in 0..3 {
            ncs_oper.append_raw_value(format_cif_f64(transform.matrix()[row][column]));
        }
        ncs_oper.append_raw_value(format_cif_f64(transform.translation()[row]));
    }
}

#[cfg(test)]
mod bio_ncs_writer_tests {
    use super::super::BioMmcifWriteError;
    use super::write_ncs_oper;
    use crate::cif::{CifBlock, CifCheckLevel, CifItem, read_cif_document};
    use cosmolkit_bio::{
        BioCoordinateBlock, BioCoordinateFormat, BioNcsOperator, BioStructureData,
        BioStructureParts, BioStructureSourceState, BioTransform,
    };
    use std::collections::BTreeMap;

    const TAGS: [&str; 14] = [
        "_struct_ncs_oper.id",
        "_struct_ncs_oper.code",
        "_struct_ncs_oper.matrix[1][1]",
        "_struct_ncs_oper.matrix[1][2]",
        "_struct_ncs_oper.matrix[1][3]",
        "_struct_ncs_oper.vector[1]",
        "_struct_ncs_oper.matrix[2][1]",
        "_struct_ncs_oper.matrix[2][2]",
        "_struct_ncs_oper.matrix[2][3]",
        "_struct_ncs_oper.vector[2]",
        "_struct_ncs_oper.matrix[3][1]",
        "_struct_ncs_oper.matrix[3][2]",
        "_struct_ncs_oper.matrix[3][3]",
        "_struct_ncs_oper.vector[3]",
    ];

    const M: [[f64; 3]; 3] = [[1.25, -0.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]];
    const M_VEC: [f64; 3] = [0.125, -2.0, 3.5];
    const T: [[f64; 3]; 3] = [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]];
    const T_VEC: [f64; 3] = [9.0, 10.0, 11.0];

    const I_ROW: [&str; 14] = [
        "I", "given", "1", "0", "0", "0", "0", "1", "0", "0", "0", "0", "1", "0",
    ];
    const M_ROW_A: [&str; 14] = [
        "A", "generate", "1.25", "-0", "2", "0.125", "3", "4", "5", "-2", "6", "7", "8", "3.5",
    ];
    const T_ROW_B: [&str; 14] = [
        "B", "given", "2", "0", "0", "9", "0", "3", "0", "10", "0", "0", "4", "11",
    ];
    const T_ROW_A: [&str; 14] = [
        "A", "given", "2", "0", "0", "9", "0", "3", "0", "10", "0", "0", "4", "11",
    ];

    fn op(id: &str, given: bool, matrix: [[f64; 3]; 3], vector: [f64; 3]) -> BioNcsOperator {
        BioNcsOperator::new(id.to_string(), given, BioTransform::new(matrix, vector))
    }

    fn data_with(operators: Vec<BioNcsOperator>, identity: Option<&str>) -> BioStructureData {
        let mut info = BTreeMap::new();
        if let Some(id) = identity {
            info.insert("_struct_ncs_oper.id".to_string(), id.to_string());
        }
        BioStructureData::from_parts(BioStructureParts {
            input_format: BioCoordinateFormat::Mmcif,
            models: vec![],
            chains: vec![],
            residues: vec![],
            atoms: vec![],
            entities: vec![],
            connections: vec![],
            cispeps: vec![],
            mod_residues: vec![],
            helices: vec![],
            sheets: vec![],
            metadata: Default::default(),
            source_state: BioStructureSourceState {
                info,
                ..BioStructureSourceState::default()
            },
            coordinates: BioCoordinateBlock::default(),
            crystal: None,
            ncs_operators: operators,
            assemblies: vec![],
        })
        .expect("empty hierarchy with NCS operators is valid")
    }

    fn parse_block(text: &str, name: &str) -> CifBlock {
        read_cif_document(text, name, CifCheckLevel::Syntax)
            .expect("frozen block state parses under the Syntax owner")
            .blocks()[0]
            .clone()
    }

    const F0: &str = "data_f0\n_entry.id before\n_after.id tail\n";
    const F1: &str = "data_f1\n_entry.id before\n_struct_ncs_oper.id OLD\n_after.id tail\n";
    const F2: &str = "data_f2\n_entry.id before\nloop_\n_struct_ncs_oper.id\n_struct_ncs_oper.code\nOLD stale\n_after.id tail\n";
    const F3: &str = "data_f3\n_entry.id before\nloop_\n_STRUCT_NCS_OPER.id\n_STRUCT_NCS_OPER.code\nOLD stale\n_after.id tail\n";

    /// Frozen per-invocation input capture: whole data value, op-id bytes,
    /// given flags, all matrix/vector bits (signed zero included), the
    /// identity-id bytes, the whole source-state info map, whole block.
    struct Inputs {
        data: BioStructureData,
        ids: Vec<Vec<u8>>,
        givens: Vec<bool>,
        bits: Vec<[u64; 12]>,
        identity: Option<Vec<u8>>,
        info: Vec<(String, String)>,
        block: CifBlock,
    }

    fn capture(data: &BioStructureData, block: CifBlock) -> Inputs {
        let bits = data
            .ncs_operators()
            .iter()
            .map(|operator| {
                let mut row_bits = [0u64; 12];
                let mut slot = 0;
                for matrix_row in operator.transform.matrix() {
                    for value in matrix_row {
                        row_bits[slot] = value.to_bits();
                        slot += 1;
                    }
                }
                for value in operator.transform.translation() {
                    row_bits[slot] = value.to_bits();
                    slot += 1;
                }
                row_bits
            })
            .collect();
        Inputs {
            data: data.clone(),
            ids: data
                .ncs_operators()
                .iter()
                .map(|operator| operator.id.as_bytes().to_vec())
                .collect(),
            givens: data.ncs_operators().iter().map(|o| o.given).collect(),
            bits,
            identity: data.ncs_oper_identity_id().map(|id| id.as_bytes().to_vec()),
            info: data
                .source_state()
                .info
                .iter()
                .map(|(key, value)| (key.clone(), value.clone()))
                .collect(),
            block,
        }
    }

    /// Literal bit table for one operator: 9 matrix components then 3
    /// translation components, all explicit to_bits() of the frozen
    /// literal matrices/vectors (signed zero included).
    fn literal_bits(matrix: &[[f64; 3]; 3], vector: &[f64; 3]) -> [u64; 12] {
        let mut bits = [0u64; 12];
        let mut slot = 0;
        for matrix_row in matrix {
            for value in matrix_row {
                bits[slot] = value.to_bits();
                slot += 1;
            }
        }
        for value in vector {
            bits[slot] = value.to_bits();
            slot += 1;
        }
        bits
    }

    /// Collected (non-panicking) preservation discrepancies for one call.
    fn preserved_discrepancies(
        before: &Inputs,
        data: &BioStructureData,
        label: &str,
    ) -> Vec<String> {
        let mut found = Vec::new();
        if data != &before.data {
            found.push(format!("{label}: data value mutated"));
        }
        let ids: Vec<Vec<u8>> = data
            .ncs_operators()
            .iter()
            .map(|operator| operator.id.as_bytes().to_vec())
            .collect();
        if ids != before.ids {
            found.push(format!("{label}: op id bytes mutated"));
        }
        let givens: Vec<bool> = data.ncs_operators().iter().map(|o| o.given).collect();
        if givens != before.givens {
            found.push(format!("{label}: given flags mutated"));
        }
        let bits: Vec<[u64; 12]> = data
            .ncs_operators()
            .iter()
            .map(|operator| {
                let mut row_bits = [0u64; 12];
                let mut slot = 0;
                for matrix_row in operator.transform.matrix() {
                    for value in matrix_row {
                        row_bits[slot] = value.to_bits();
                        slot += 1;
                    }
                }
                for value in operator.transform.translation() {
                    row_bits[slot] = value.to_bits();
                    slot += 1;
                }
                row_bits
            })
            .collect();
        if bits != before.bits {
            found.push(format!("{label}: transform bits mutated"));
        }
        let identity = data.ncs_oper_identity_id().map(|id| id.as_bytes().to_vec());
        if identity != before.identity {
            found.push(format!("{label}: identity bytes mutated"));
        }
        let info: Vec<(String, String)> = data
            .source_state()
            .info
            .iter()
            .map(|(key, value)| (key.clone(), value.clone()))
            .collect();
        if info != before.info {
            found.push(format!("{label}: source-state map mutated"));
        }
        found
    }

    fn ncs_loop(block: &CifBlock) -> Option<(&[String], &[crate::cif::CifValue])> {
        for item in block.items() {
            if let CifItem::Loop(row) = item
                && row
                    .tags()
                    .first()
                    .is_some_and(|tag| tag.starts_with("_struct_ncs_oper."))
            {
                return Some((row.tags(), row.values()));
            }
        }
        None
    }

    /// Setup-time (pre-invocation) validation of the parsed block literal:
    /// actual tags, raw values, dimensions and positions per frozen state.
    /// A failure here is an invalid fixture — a real setup failure.
    fn verify_block_state(block: &CifBlock, name: &str) {
        let items = block.items();
        match name {
            "F0" => {
                assert_eq!(items.len(), 2, "F0 literal dimension");
                let CifItem::Pair(entry) = &items[0] else {
                    panic!("F0 entry pair position");
                };
                assert_eq!(
                    (entry.tag(), entry.value().map(|value| value.raw())),
                    ("_entry.id", Some("before"))
                );
                let CifItem::Pair(tail) = &items[1] else {
                    panic!("F0 tail pair position");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|value| value.raw())),
                    ("_after.id", Some("tail"))
                );
            }
            "F1" => {
                assert_eq!(items.len(), 3, "F1 literal dimension");
                let CifItem::Pair(entry) = &items[0] else {
                    panic!("F1 entry pair position");
                };
                assert_eq!(
                    (entry.tag(), entry.value().map(|value| value.raw())),
                    ("_entry.id", Some("before"))
                );
                let CifItem::Pair(old) = &items[1] else {
                    panic!("F1 old category pair position");
                };
                assert_eq!(
                    (old.tag(), old.value().map(|value| value.raw())),
                    ("_struct_ncs_oper.id", Some("OLD"))
                );
                let CifItem::Pair(tail) = &items[2] else {
                    panic!("F1 tail pair position");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|value| value.raw())),
                    ("_after.id", Some("tail"))
                );
            }
            "F2" | "F3" => {
                let prefix = if name == "F2" {
                    "_struct_ncs_oper"
                } else {
                    "_STRUCT_NCS_OPER"
                };
                assert_eq!(items.len(), 3, "{name} literal dimension");
                let CifItem::Pair(entry) = &items[0] else {
                    panic!("{name} entry pair position");
                };
                assert_eq!(
                    (entry.tag(), entry.value().map(|value| value.raw())),
                    ("_entry.id", Some("before"))
                );
                let CifItem::Loop(stale) = &items[1] else {
                    panic!("{name} stale loop position");
                };
                assert_eq!(
                    stale.tags(),
                    [format!("{prefix}.id"), format!("{prefix}.code")],
                    "{name} stale tags"
                );
                assert_eq!(
                    stale.values().iter().map(|v| v.raw()).collect::<Vec<_>>(),
                    ["OLD", "stale"],
                    "{name} stale raw values"
                );
                let CifItem::Pair(tail) = &items[2] else {
                    panic!("{name} tail pair position");
                };
                assert_eq!(
                    (tail.tag(), tail.value().map(|value| value.raw())),
                    ("_after.id", Some("tail"))
                );
            }
            other => panic!("unknown frozen block state {other}"),
        }
    }

    /// Setup-time (pre-invocation) per-call operator prerequisites checked
    /// against the literal id/order/given/bit tables and the identity/map.
    fn verify_operators(
        data: &BioStructureData,
        spec: &[(&str, bool, &[[f64; 3]; 3], &[f64; 3])],
        identity: Option<&str>,
        label: &str,
    ) {
        let operators = data.ncs_operators();
        assert_eq!(
            operators.len(),
            spec.len(),
            "{label}: operator count prerequisite"
        );
        for (index, (id, given, matrix, vector)) in spec.iter().enumerate() {
            let operator = &operators[index];
            assert_eq!(operator.id, *id, "{label}: op {index} id prerequisite");
            assert_eq!(
                operator.given, *given,
                "{label}: op {index} given prerequisite"
            );
            let mut actual = [0u64; 12];
            let mut slot = 0;
            for matrix_row in operator.transform.matrix() {
                for value in matrix_row {
                    actual[slot] = value.to_bits();
                    slot += 1;
                }
            }
            for value in operator.transform.translation() {
                actual[slot] = value.to_bits();
                slot += 1;
            }
            assert_eq!(
                actual,
                literal_bits(matrix, vector),
                "{label}: op {index} transform bits prerequisite"
            );
        }
        assert_eq!(
            data.ncs_oper_identity_id(),
            identity,
            "{label}: identity id prerequisite"
        );
        if let Some(id) = identity {
            let matching = data
                .source_state()
                .info
                .iter()
                .filter(|(key, value)| {
                    key.as_str() == "_struct_ncs_oper.id" && value.as_str() == id
                })
                .count();
            assert_eq!(
                data.source_state().info.len(),
                matching,
                "{label}: source-state map prerequisite"
            );
        } else {
            assert!(
                data.source_state().info.is_empty(),
                "{label}: source-state map prerequisite"
            );
        }
    }

    #[test]
    fn bio_ncs_writer_product64() {
        let mut calls = 0usize;
        let mut new_rows = 0usize;
        let mut inserted_identity = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        type Spec = (
            &'static str,
            bool,
            &'static [[f64; 3]; 3],
            &'static [f64; 3],
        );
        let shapes: [(&str, &[Spec], Option<&str>, &[&[&str; 14]]); 8] = [
            ("E0", &[], None, &[]),
            ("E1", &[], Some("I"), &[]),
            ("E2", &[("A", false, &M, &M_VEC)], None, &[&M_ROW_A]),
            (
                "E3",
                &[("A", false, &M, &M_VEC)],
                Some("I"),
                &[&I_ROW, &M_ROW_A],
            ),
            ("E4", &[("A", false, &M, &M_VEC)], Some("A"), &[&M_ROW_A]),
            (
                "E5",
                &[("A", false, &M, &M_VEC), ("B", true, &T, &T_VEC)],
                None,
                &[&M_ROW_A, &T_ROW_B],
            ),
            (
                "E6",
                &[("B", true, &T, &T_VEC), ("A", false, &M, &M_VEC)],
                Some("A"),
                &[&T_ROW_B, &M_ROW_A],
            ),
            (
                "E7",
                &[("A", false, &M, &M_VEC), ("A", true, &T, &T_VEC)],
                Some("I"),
                &[&I_ROW, &M_ROW_A, &T_ROW_A],
            ),
        ];
        let block_states: [(&str, &str, bool); 4] = [
            ("F0", F0, false),
            ("F1", F1, true),
            ("F2", F2, true),
            ("F3", F3, true),
        ];

        for (shape_name, operators, identity, expected_rows) in shapes {
            let built: Vec<BioNcsOperator> = operators
                .iter()
                .map(|(id, given, matrix, vector)| op(*id, *given, **matrix, **vector))
                .collect();
            let data = data_with(built, identity);

            for (block_name, block_text, category_between) in block_states {
                for repetition in 0..2 {
                    let label = format!("{shape_name}/{block_name}/{repetition}");
                    let mut block = parse_block(block_text, &label);
                    // Complete per-call prerequisites immediately BEFORE
                    // this invocation: literal operator tables (ids, order,
                    // given flags, 9+3 transform bits incl. the M signed
                    // zero), identity/map values, and the parsed block's
                    // actual tags/raw values/dimensions/positions.
                    verify_operators(&data, operators, identity, &label);
                    verify_block_state(&block, block_name);

                    let before = capture(&data, block.clone());
                    let result = write_ncs_oper(&data, &mut block);
                    // calls increments only AFTER the real invocation.
                    calls += 1;
                    // Preservation AFTER Result, BEFORE any branch/output
                    // handling — collected, never panicking.
                    discrepancies.extend(preserved_discrepancies(&before, &data, &label));

                    if expected_rows.is_empty() {
                        if result.is_err() {
                            discrepancies.push(format!(
                                "{label}: unexpected error {:?}",
                                result.err().map(|error| error.to_string())
                            ));
                        }
                        if block != before.block {
                            discrepancies.push(format!("{label}: block mutated"));
                        }
                        // E0/E1: no category may be invented; the untouched
                        // OLD/stale rows are NOT counted as newly emitted.
                    } else {
                        if let Err(error) = &result {
                            discrepancies.push(format!("{label}: unexpected error: {error}"));
                        }

                        if let Some((tags, values)) = ncs_loop(&block) {
                            if tags != &TAGS {
                                discrepancies.push(format!("{label}: tags {tags:?}"));
                            }
                            if values.len() % 14 != 0 {
                                discrepancies.push(format!(
                                    "{label}: incomplete width {} values",
                                    values.len()
                                ));
                            }
                            if values.len() != expected_rows.len() * 14 {
                                discrepancies.push(format!(
                                    "{label}: value count {} != {}",
                                    values.len(),
                                    expected_rows.len() * 14
                                ));
                            }
                            // OBSERVED row census from actual output values/width.
                            if !values.is_empty() && values.len() % 14 == 0 {
                                new_rows += values.len() / 14;
                            }
                            for (row_index, expected) in expected_rows.iter().enumerate() {
                                for column in 0..14 {
                                    let actual = values
                                        .get(row_index * 14 + column)
                                        .map(|value| value.raw());
                                    if actual != Some(expected[column]) {
                                        discrepancies.push(format!(
                                            "{label}: cell {row_index}/{column} {actual:?} != {}",
                                            expected[column]
                                        ));
                                    }
                                }
                            }
                            // OBSERVED inserted-identity census: count actual full
                            // 14-cell rows whose complete content is the identity
                            // row carrying this shape's identity id.
                            if let Some(identity_id) = identity {
                                for row_index in 0..values.len() / 14 {
                                    let actual: Vec<&str> = (0..14)
                                        .filter_map(|column| {
                                            values.get(row_index * 14 + column).map(|v| v.raw())
                                        })
                                        .collect();
                                    if actual.len() == 14
                                        && actual[0] == identity_id
                                        && actual[1] == "given"
                                        && actual[2..] == I_ROW[2..]
                                    {
                                        inserted_identity += 1;
                                    }
                                }
                            }
                        } else {
                            discrepancies.push(format!("{label}: no _struct_ncs_oper loop"));
                        }

                        // Check all unrelated items even on unexpected errors
                        // or a missing generated category; neither can skip
                        // independent post-call block checks for this call.
                        let after_items = block.items();
                        let loop_index = after_items.iter().position(|item| {
                            matches!(item, CifItem::Loop(row) if row
                            .tags()
                            .first()
                            .is_some_and(|tag| tag.starts_with("_struct_ncs_oper.")))
                        });
                        let expected_loop_index = if category_between { 1 } else { 2 };
                        if loop_index != Some(expected_loop_index) {
                            discrepancies.push(format!(
                                "{label}: loop position {loop_index:?} != {expected_loop_index}"
                            ));
                        }
                        if after_items.len() != 3 {
                            discrepancies
                                .push(format!("{label}: item count {} != 3", after_items.len()));
                        }
                        let entry_ok = after_items.first().is_some_and(|item| {
                            matches!(item, CifItem::Pair(pair) if pair.tag() == "_entry.id"
                            && pair.value().is_some_and(|value| value.raw() == "before"))
                        });
                        if !entry_ok {
                            discrepancies.push(format!("{label}: entry pair moved"));
                        }
                        // F0: the new loop is appended at the END (source
                        // setup_loop emplaces a fresh Loop item last), so the
                        // tail pair stays at index 1; category-carrying states
                        // replace the category at index 1 and keep the tail at 2.
                        let tail_index = if category_between { 2 } else { 1 };
                        let tail_ok = after_items.get(tail_index).is_some_and(|item| {
                            matches!(item, CifItem::Pair(pair) if pair.tag() == "_after.id"
                            && pair.value().is_some_and(|value| value.raw() == "tail"))
                        });
                        if !tail_ok {
                            discrepancies.push(format!("{label}: tail pair moved"));
                        }
                    }
                }
            }
        }

        assert_eq!(calls, 64, "exact 64 actual write_ncs_oper calls");
        assert_eq!(new_rows, 88, "exact 88 newly emitted rows");
        assert_eq!(inserted_identity, 16, "exact 16 inserted identity rows");
        assert!(
            discrepancies.is_empty(),
            "product discrepancies: {discrepancies:?}"
        );
    }

    #[test]
    fn bio_ncs_writer_mixed_error4() {
        let mixed_lower = "data_m0\nloop_\n_struct_ncs_oper.id\n_foreign.id\nOLD tail\n";
        let mixed_upper = "data_m1\nloop_\n_STRUCT_NCS_OPER.id\n_foreign.id\nOLD tail\n";
        let mut calls = 0usize;
        let mut discrepancies: Vec<String> = Vec::new();

        for (name, text) in [("lower", mixed_lower), ("upper", mixed_upper)] {
            for (shape, spec) in [
                ("E0", [].as_slice()),
                ("E2", [("A", false, &M, &M_VEC)].as_slice()),
            ] {
                let label = format!("{name}/{shape}");
                let expect_error = shape == "E2";
                let built: Vec<BioNcsOperator> = spec
                    .iter()
                    .map(|(id, given, matrix, vector)| op(*id, *given, **matrix, **vector))
                    .collect();
                let data = data_with(built, None);
                let mut block = parse_block(text, &label);
                // Per-call prerequisites: operator tables and the parsed
                // mixed loop's actual tags/raw values/line position.
                verify_operators(&data, spec, None, &label);
                assert_eq!(block.items().len(), 1, "{label}: literal dimension");
                let CifItem::Loop(mixed) = &block.items()[0] else {
                    panic!("{label}: mixed loop position");
                };
                let expected_tag = if name == "lower" {
                    "_struct_ncs_oper.id"
                } else {
                    "_STRUCT_NCS_OPER.id"
                };
                assert_eq!(
                    mixed.tags(),
                    [expected_tag, "_foreign.id"],
                    "{label}: mixed tags"
                );
                assert_eq!(
                    mixed.values().iter().map(|v| v.raw()).collect::<Vec<_>>(),
                    ["OLD", "tail"],
                    "{label}: mixed raw values"
                );
                assert_eq!(mixed.line(), 2, "{label}: mixed loop literal line");

                let before = capture(&data, block.clone());
                let result = write_ncs_oper(&data, &mut block);
                calls += 1;
                discrepancies.extend(preserved_discrepancies(&before, &data, &label));

                // WHOLE original block is checked BEFORE Result dispatch,
                // including unexpected Ok or a foreign error variant.
                if block != before.block {
                    discrepancies.push(format!("{label}: block mutated"));
                }

                if !expect_error {
                    if let Err(error) = &result {
                        discrepancies.push(format!("{label}: unexpected error: {error}"));
                    }
                } else {
                    match &result {
                        // Borrow the actual error without flattening/cloning.
                        Err(error @ BioMmcifWriteError::Cif(cif)) => {
                            if cif.kind() != crate::cif::CifReadErrorKind::InvalidLoop {
                                discrepancies.push(format!("{label}: kind {:?}", cif.kind()));
                            }
                            if cif.message() != "Tag _foreign.id in loop with _struct_ncs_oper." {
                                discrepancies.push(format!("{label}: message {}", cif.message()));
                            }
                            if cif.source() != "cif" {
                                discrepancies.push(format!("{label}: source {}", cif.source()));
                            }
                            if cif.line() != 2 {
                                discrepancies.push(format!("{label}: line {}", cif.line()));
                            }
                            if cif.column() != 1 {
                                discrepancies.push(format!("{label}: column {}", cif.column()));
                            }
                            // Borrowed std::error::Error::source downcast must be
                            // the SAME CifReadError instance (same address).
                            let borrowed_source =
                                std::error::Error::source(error).and_then(|source| {
                                    source.downcast_ref::<crate::cif::CifReadError>()
                                });
                            match borrowed_source {
                                Some(downcast) if std::ptr::eq(downcast, cif) => {}
                                other => discrepancies
                                    .push(format!("{label}: borrowed source identity {other:?}")),
                            }
                        }
                        Err(error) => {
                            discrepancies.push(format!("{label}: wrong error variant: {error}"))
                        }
                        Ok(()) => discrepancies.push(format!("{label}: expected error, got Ok")),
                    }
                }
            }
        }
        assert_eq!(calls, 4, "exact 4 mixed-error calls");
        assert!(
            discrepancies.is_empty(),
            "mixed-error discrepancies: {discrepancies:?}"
        );
    }
}
