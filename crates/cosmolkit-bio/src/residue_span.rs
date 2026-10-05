//! Source-defined sequence-address translations used by structural metadata writers.
use crate::{BioResidueRow, BioStructureError, PdbSeqId};
fn num(row: &BioResidueRow) -> Option<i32> {
    row.source()
        .seq_id()
        .map(|s| s.seq_num())
        .filter(|n| *n != i32::MIN)
}
fn less(a: Option<i32>, b: Option<i32>) -> bool {
    // Gemmi❗✔️:   bool operator<(const OptionalInt& o) const {
    // Gemmi❗✔️:     return has_value() && o.has_value() && value < o.value;
    // Gemmi❗✔️:   }
    matches!((a,b), (Some(a),Some(b)) if a < b)
}
fn plus(a: Option<i32>, b: Option<i32>) -> Option<i32> {
    // Gemmi❗✔️:   OptionalInt operator+(OptionalInt o) const {
    // Gemmi❗✔️:     return OptionalInt(has_value() && o.has_value() ? value + o.value : N);
    // Gemmi❗✔️:   }
    a.zip(b)
        .map(|(a, b)| a.wrapping_add(b))
        .filter(|n| *n != i32::MIN)
}
fn minus(a: Option<i32>, b: Option<i32>) -> Option<i32> {
    // Gemmi❗✔️:   OptionalInt operator-(OptionalInt o) const {
    // Gemmi❗✔️:     return OptionalInt(has_value() && o.has_value() ? value - o.value : N);
    // Gemmi❗✔️:   }
    a.zip(b)
        .map(|(a, b)| a.wrapping_sub(b))
        .filter(|n| *n != i32::MIN)
}
/// Translate label numbering using the pinned nearest endpoint and interpolation rules.
pub fn label_seq_id_to_auth(
    rows: &[&BioResidueRow],
    label: Option<i32>,
) -> Result<Option<PdbSeqId>, BioStructureError> {
    // Gemmi❗✔️:   SeqId label_seq_id_to_auth(SeqId::OptionalNum label_seq_id) const {
    // Gemmi❗✔️:     if (empty())
    // Gemmi❗✔️:       throw std::out_of_range("label_seq_id_to_auth(): empty span");
    // Gemmi❗✔️:     const auto* it = std::lower_bound(begin(), end(), label_seq_id,
    // Gemmi❗✔️:         [](const Residue& r, SeqId::OptionalNum v){ return r.label_seq < v; });
    // Gemmi❗✔️:     if (it == end())
    // Gemmi❗✔️:       --it;
    // Gemmi❗✔️:     else if (it->label_seq == label_seq_id)
    // Gemmi❗✔️:       return it->seqid;
    // Gemmi❗✔️:     else if (it != begin() &&
    // Gemmi❗✔️:              label_seq_id - (it-1)->label_seq < it->label_seq - label_seq_id)
    // Gemmi❗✔️:       --it;
    // Gemmi❗✔️:     return {it->seqid.num + (label_seq_id - it->label_seq), ' '};
    // Gemmi❗✔️:   }
    // Cost: the same binary lower_bound and at most one neighbor comparison; no allocation.
    if rows.is_empty() {
        return Err(BioStructureError::EmptyResidueSpan {
            operation: "label_seq_id_to_auth",
        });
    }
    let mut index = rows.partition_point(|row| less(row.source().label_seq_id(), label));
    if index == rows.len() {
        index -= 1;
    } else if rows[index].source().label_seq_id() == label {
        return Ok(rows[index].source().seq_id());
    } else if index != 0
        && less(
            minus(label, rows[index - 1].source().label_seq_id()),
            minus(rows[index].source().label_seq_id(), label),
        )
    {
        index -= 1;
    }
    Ok(plus(
        num(rows[index]),
        minus(label, rows[index].source().label_seq_id()),
    )
    .map(|n| PdbSeqId::new(n, None)))
}
/// Translate author numbering with the source exact-match pass followed by source interpolation.
pub fn auth_seq_id_to_label(
    rows: &[&BioResidueRow],
    auth: Option<PdbSeqId>,
) -> Result<Option<i32>, BioStructureError> {
    // Gemmi❗✔️:   SeqId::OptionalNum auth_seq_id_to_label(SeqId auth_seq_id) const {
    // Gemmi❗✔️:     if (empty())
    // Gemmi❗✔️:       throw std::out_of_range("auth_seq_id_to_label(): empty span");
    // Gemmi❗✔️:     for (const Residue& r : *this)
    // Gemmi❗✔️:       if (r.seqid == auth_seq_id)
    // Gemmi❗✔️:         return r.label_seq;
    // Gemmi❗✔️:     const_iterator it;
    // Gemmi❗✔️:     if (auth_seq_id.num < front().seqid.num) {
    // Gemmi❗✔️:       it = begin();
    // Gemmi❗✔️:     } else if (back().seqid.num < auth_seq_id.num) {
    // Gemmi❗✔️:       it = end() - 1;
    // Gemmi❗✔️:     } else {
    // Gemmi❗✔️:       it = std::lower_bound(begin(), end(), auth_seq_id.num,
    // Gemmi❗✔️:         [](const Residue& r, SeqId::OptionalNum v){ return r.seqid.num < v; });
    // Gemmi❗✔️:       while (it != end() && it->seqid.num == auth_seq_id.num &&
    // Gemmi❗✔️:              it->seqid.icode != auth_seq_id.icode)
    // Gemmi❗✔️:         ++it;
    // Gemmi❗✔️:       if (it == end())
    // Gemmi❗✔️:         --it;
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:     return it->label_seq + (auth_seq_id.num - it->seqid.num);
    // Gemmi❗✔️:   }
    // Gemmi❗✔️: };
    // Cost: one exact-match scan followed by the same lower_bound, with no temporary storage.
    if rows.is_empty() {
        return Err(BioStructureError::EmptyResidueSpan {
            operation: "auth_seq_id_to_label",
        });
    }
    let equal = |a: Option<PdbSeqId>, b: Option<PdbSeqId>| {
        a.map(|s| s.seq_num()).filter(|n| *n != i32::MIN)
            == b.map(|s| s.seq_num()).filter(|n| *n != i32::MIN)
            && ((a.and_then(|s| s.ins_code()).unwrap_or(b' ')
                ^ b.and_then(|s| s.ins_code()).unwrap_or(b' '))
                & !0x20)
                == 0
    };
    for row in rows {
        if equal(row.source().seq_id(), auth) {
            return Ok(row.source().label_seq_id());
        }
    }
    let auth_num = auth.map(|s| s.seq_num()).filter(|n| *n != i32::MIN);
    let index = if less(auth_num, num(rows[0])) {
        0
    } else if less(num(rows[rows.len() - 1]), auth_num) {
        rows.len() - 1
    } else {
        let mut index = rows.partition_point(|row| less(num(row), auth_num));
        while index != rows.len()
            && num(rows[index]) == auth_num
            && rows[index]
                .source()
                .seq_id()
                .and_then(|s| s.ins_code())
                .unwrap_or(b' ')
                != auth.and_then(|s| s.ins_code()).unwrap_or(b' ')
        {
            index += 1;
        }
        if index == rows.len() {
            index -= 1;
        }
        index
    };
    Ok(plus(
        rows[index].source().label_seq_id(),
        minus(auth_num, num(rows[index])),
    ))
}
