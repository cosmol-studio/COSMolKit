//! Pinned Gemmi crystallographic spacegroup data owned by the BIO domain.
//!
//! Main and alternate rows derive from Gemmi `symmetry.cpp` at commit
//! `5cc1c23c6007e0e6cbd69289c6f7c0bff50e943e`. The operation arrays reuse the
//! independently verified historical data only after exact ordered comparison
//! with each pinned `SpaceGroup::operations()` result (lane report Step22).
//! Gemmi source is MPL-2.0; see `third_party/gemmi/LICENSE.txt`.

pub(super) const GEMMI_OP_DEN: i32 = 24;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct BioSymOp {
    rot: [[i32; 3]; 3],
    tran: [i32; 3],
}

#[derive(Clone, Copy, Debug)]
struct BioSpaceGroup {
    number: i32,
    ccp4: i32,
    hm: &'static str,
    ext: u8,
    qualifier: &'static str,
    hall: &'static str,
    basisop_idx: i32,
    ops: &'static [BioSymOp],
}

#[derive(Clone, Copy, Debug)]
struct BioAltName {
    hm: &'static str,
    ext: u8,
    pos: usize,
}

#[derive(Debug, PartialEq, Eq)]
enum ParsedSpaceGroupName {
    Numeric { ccp4: i32 },
    Named(NormalizedSpaceGroupName),
}

#[derive(Debug, PartialEq, Eq)]
struct NormalizedSpaceGroupName {
    bytes: Vec<u8>,
    first: u8,
    start: usize,
    c_string_end: usize,
}

// Gemmi❗✔️: static const char* skip_space(const char* p) {
// Gemmi❗✔️:   if (p)
// Gemmi❗✔️:     while (*p == ' ' || *p == '\t' || *p == '_') // '_' can be used as space
// Gemmi❗✔️:       ++p;
// Gemmi❗✔️:   return p;
// Gemmi❗✔️: }
// Behavior review: the byte cursor stops at the first NUL, matching C-string
// traversal; only space, tab, and underscore are skipped.
// Complexity review: advances one cursor over the skipped prefix with no
// allocation or repeated scan.
fn skip_space(bytes: &[u8], mut index: usize, c_string_end: usize) -> usize {
    while index < c_string_end && matches!(bytes[index], b' ' | b'\t' | b'_') {
        index += 1;
    }
    index
}

// Gemmi❗✔️: inline const SpaceGroup* find_spacegroup_by_number(int ccp4) noexcept {
// Gemmi❗✔️:   if (ccp4 == 0)
// Gemmi❗✔️:     return &spacegroup_tables::main[0];
// Gemmi❗✔️:   for (const SpaceGroup& sg : spacegroup_tables::main)
// Gemmi❗✔️:     if (sg.ccp4 == ccp4)
// Gemmi❗✔️:       return &sg;
// Gemmi❗✔️:   return nullptr;
// Gemmi❗✔️: }
// Behavior review: input zero selects main row zero; otherwise the first
// ordered row with the exact CCP4 key wins. The international number and row
// position are not lookup keys.
// Complexity review: one ordered scan of the fixed 564-row table, matching
// Gemmi's O(n) lookup and returning without allocating.
fn find_spacegroup_by_number(ccp4: i32) -> Option<&'static BioSpaceGroup> {
    if ccp4 == 0 {
        return GEMMI_SPACEGROUPS.first();
    }
    GEMMI_SPACEGROUPS.iter().find(|group| group.ccp4 == ccp4)
}

// Gemmi❗✔️: const char* p = skip_space(name.c_str());
// Gemmi❗✔️: if (*p >= '0' && *p <= '9') { // handle numbers
// Gemmi❗✔️:   char *endptr;
// Gemmi❗✔️:   long n = std::strtol(p, &endptr, 10);
// Gemmi❗✔️:   return *endptr == '\0' ? find_spacegroup_by_number(n) : nullptr;
// Gemmi❗✔️: }
// Behavior review: `strtol` positive overflow saturates at `LONG_MAX` while
// still consuming the digit string; Gemmi ignores `errno` and narrows its long
// to the `int` lookup parameter. Rust retains that ABI-width saturation and
// narrowing; out-of-range long-to-int parity remains subject to the fixed
// reference cases because C++ leaves that conversion implementation-defined.
// Complexity review: digit recognition/parsing is one input scan; lookup is
// one fixed ordered table scan. No per-entry conversion or secondary index is
// introduced.
fn parse_spacegroup_name(name: &str) -> Option<ParsedSpaceGroupName> {
    let input = name.as_bytes();
    let c_string_end = input
        .iter()
        .position(|&byte| byte == 0)
        .unwrap_or(input.len());
    let p = skip_space(input, 0, c_string_end);
    if p == c_string_end {
        return None;
    }

    if input[p] >= b'0' && input[p] <= b'9' {
        let digits_end = input[p..c_string_end]
            .iter()
            .position(|byte| !byte.is_ascii_digit())
            .map_or(c_string_end, |offset| p + offset);
        if digits_end != c_string_end {
            return None;
        }
        let digits = std::str::from_utf8(&input[p..digits_end]).ok()?;
        let number = digits
            .parse::<std::ffi::c_long>()
            .unwrap_or(std::ffi::c_long::MAX);
        return Some(ParsedSpaceGroupName::Numeric {
            ccp4: number as i32,
        });
    }

    normalize_spacegroup_name(input, p, c_string_end).map(ParsedSpaceGroupName::Named)
}

// Gemmi❗✔️: char first = *p & ~0x20; // to uppercase
// Gemmi❗✔️: if (first == '\0')
// Gemmi❗✔️:   return nullptr;
// Gemmi❗✔️: if (first == 'H')
// Gemmi❗✔️:   first = 'R';
// Gemmi❗✔️: p = skip_space(p+1);
// Gemmi❗✔️: size_t start = p - name.c_str();
// Gemmi❗✔️: // change letters to lower case, except the letter after :
// Gemmi❗✔️: for (size_t i = start; i < name.size(); ++i) {
// Gemmi❗✔️:   if (name[i] >= 'A' && name[i] <= 'Z')
// Gemmi❗✔️:     name[i] |= 0x20;  // to lowercase
// Gemmi❗✔️:   else if (name[i] == ':')
// Gemmi❗✔️:     while (++i < name.size())
// Gemmi❗✔️:       if (name[i] >= 'a' && name[i] <= 'z')
// Gemmi❗✔️:         name[i] &= ~0x20;  // to uppercase
// Gemmi❗✔️: }
// Gemmi❗✔️: // allow names ending with R or H, such as R3R instead of R3:R
// Gemmi❗✔️: if (name.back() == 'h' || name.back() == 'r') {
// Gemmi❗✔️:   name.back() &= ~0x20;  // to uppercase
// Gemmi❗✔️:   name.insert(name.end() - 1, ':');
// Gemmi❗✔️: }
// Behavior review: normalization is byte-for-byte ASCII-only, retains the
// source's C-string boundary for lookup, and applies the trailing H/R rule to
// the full owned string as Gemmi does. No whitespace classes beyond
// `skip_space` are normalized.
// Complexity review: one owned byte copy and one linear normalization pass;
// this matches the source's by-value `std::string` mutation and adds no scan
// per table row.
fn normalize_spacegroup_name(
    input: &[u8],
    first_position: usize,
    mut c_string_end: usize,
) -> Option<NormalizedSpaceGroupName> {
    let mut bytes = input.to_vec();
    let mut first = bytes[first_position] & !0x20;
    if first == 0 {
        return None;
    }
    if first == b'H' {
        first = b'R';
    }
    let start = skip_space(&bytes, first_position + 1, c_string_end);

    let mut index = start;
    while index < bytes.len() {
        if bytes[index].is_ascii_uppercase() {
            bytes[index] |= 0x20;
        } else if bytes[index] == b':' {
            index += 1;
            while index < bytes.len() {
                if bytes[index].is_ascii_lowercase() {
                    bytes[index] &= !0x20;
                }
                index += 1;
            }
            break;
        }
        index += 1;
    }

    let old_len = bytes.len();
    if matches!(bytes.last(), Some(b'h' | b'r')) {
        let last = bytes.last_mut()?;
        *last &= !0x20;
        bytes.insert(old_len - 1, b':');
        if c_string_end == old_len {
            c_string_end += 1;
        }
    }

    Some(NormalizedSpaceGroupName {
        bytes,
        first,
        start,
        c_string_end,
    })
}

// Gemmi❗✔️:   for (const SpaceGroup& sg : spacegroup_tables::main)
// Gemmi❗✔️:     if (sg.hm[0] == first) {
// Gemmi❗✔️:       if (sg.hm[2] == *p) {
// Gemmi❗✔️:         const char* a = skip_space(p + 1);
// Gemmi❗✔️:         const char* b = skip_space(sg.hm + 3);
// Gemmi❗✔️:         while ((*a == *b && *b != '\0') ||
// Gemmi❗✔️:                (*a == '3' && *b == '-' && b == sg.hm + 4 && *++b == '3')) {
// Gemmi❗✔️:           a = skip_space(a+1);
// Gemmi❗✔️:           b = skip_space(b+1);
// Gemmi❗✔️:         }
// Gemmi❗✔️:         if (*b == '\0') {
// Gemmi❗✔️:           if (*a == '\0') {
// Gemmi❗✔️:             if (sg.ext == 'H' && (alpha == 0. ? prefer_R : gamma < 1.125 * alpha))
// Gemmi❗✔️:               return &sg + 1;
// Gemmi❗✔️:             if (sg.ext == '1' && prefer_2)
// Gemmi❗✔️:               return &sg + 1;
// Gemmi❗✔️:             return &sg;
// Gemmi❗✔️:           }
// Gemmi❗✔️:           if (*a == ':' && *skip_space(a+1) == sg.ext)
// Gemmi❗✔️:             return &sg;
// Gemmi❗✔️:         }
// Gemmi❗✔️:       } else if (sg.hm[2] == '1' && sg.hm[3] == ' ') {
// Gemmi❗✔️:         // check monoclinic short names, matching P2 to "P 1 2 1";
// Gemmi❗✔️:         // as an exception "B 2" == "B 1 1 2" (like in the PDB)
// Gemmi❗✔️:         const char* b = sg.hm + 4;
// Gemmi❗✔️:         if (*b != '1' || (first == 'B' && *++b == ' ' && *++b != '1')) {
// Gemmi❗✔️:           char end = (b == sg.hm + 4 ? ' ' : '\0');
// Gemmi❗✔️:           const char* a = skip_space(p);
// Gemmi❗✔️:           while (*a == *b && *b != end) {
// Gemmi❗✔️:             ++a;
// Gemmi❗✔️:             ++b;
// Gemmi❗✔️:           }
// Gemmi❗✔️:           if (*skip_space(a) == '\0' && *b == end)
// Gemmi❗✔️:             return &sg;
// Gemmi❗✔️:         }
// Gemmi❗✔️:       }
// Gemmi❗✔️:     }
// Behavior review: this helper implements the complete pinned main-table
// branch, including its short-name branch only after the full-name character
// test fails. The caller passes source preference flags; Structure uses null.
// Its B short-name conditional mutates the comparison cursor before selecting
// the endpoint, and ordered table traversal preserves the first winner.
// The source's byte cursor, first-match order, cubic spelling, adjacent H/R and
// origin settings, and explicit extension test are preserved.
// Alternate names and full structure dispatch belong to the other helpers.
// Steps38/40 and 44/46 pinned regressions now pass for their named main and
// monoclinic cases. The marker remains conservative for unenumerated spellings.
// Complexity review: one ordered scan of the same 564 fixed rows, comparing
// symbols in place without allocating or preprocessing each table row. Each
// candidate comparison advances its two cursors monotonically, matching the
// source's row-scan complexity and avoiding repeated prefix rescans.
fn find_main_spacegroup_by_name(
    name: &NormalizedSpaceGroupName,
    alpha: f64,
    gamma: f64,
    prefer_r: bool,
    prefer_2: bool,
) -> Option<&'static BioSpaceGroup> {
    let input = &name.bytes;
    for (group_index, group) in GEMMI_SPACEGROUPS.iter().enumerate() {
        let group_name = group.hm.as_bytes();
        let group_first = byte_or_nul(group_name, 0, group_name.len());
        let group_third = byte_or_nul(group_name, 2, group_name.len());
        if group_first != name.first {
            continue;
        }

        if group_third != byte_or_nul(input, name.start, name.c_string_end) {
            if group_third == b'1' && byte_or_nul(group_name, 3, group_name.len()) == b' ' {
                let mut b = 4;
                let mut short_name_matches = byte_or_nul(group_name, b, group_name.len()) != b'1';
                if !short_name_matches && name.first == b'B' {
                    b += 1;
                    if byte_or_nul(group_name, b, group_name.len()) == b' ' {
                        b += 1;
                        short_name_matches = byte_or_nul(group_name, b, group_name.len()) != b'1';
                    }
                }

                if short_name_matches {
                    let end = if b == 4 { b' ' } else { 0 };
                    let mut input_index = skip_space(input, name.start, name.c_string_end);
                    while byte_or_nul(input, input_index, name.c_string_end)
                        == byte_or_nul(group_name, b, group_name.len())
                        && byte_or_nul(group_name, b, group_name.len()) != end
                    {
                        input_index += 1;
                        b += 1;
                    }

                    let input_end = skip_space(input, input_index, name.c_string_end);
                    if byte_or_nul(input, input_end, name.c_string_end) == 0
                        && byte_or_nul(group_name, b, group_name.len()) == end
                    {
                        return Some(group);
                    }
                }
            }
            continue;
        }

        let mut input_index = skip_space(input, name.start + 1, name.c_string_end);
        let mut group_index_byte = skip_space(group_name, 3, group_name.len());
        loop {
            let input_byte = byte_or_nul(input, input_index, name.c_string_end);
            let group_byte = byte_or_nul(group_name, group_index_byte, group_name.len());
            let equal_symbol_byte = input_byte == group_byte && group_byte != 0;
            let legacy_cubic_three = input_byte == b'3'
                && group_byte == b'-'
                && group_index_byte == 4
                && byte_or_nul(group_name, group_index_byte + 1, group_name.len()) == b'3';
            if !equal_symbol_byte && !legacy_cubic_three {
                break;
            }

            // The source's `*++b` consumes the '-' test byte before the loop
            // body advances past the following `3`.
            if legacy_cubic_three {
                group_index_byte += 1;
            }
            input_index = skip_space(input, input_index + 1, name.c_string_end);
            group_index_byte = skip_space(group_name, group_index_byte + 1, group_name.len());
        }

        if byte_or_nul(group_name, group_index_byte, group_name.len()) != 0 {
            continue;
        }

        let input_tail = byte_or_nul(input, input_index, name.c_string_end);
        if input_tail == 0 {
            if group.ext == b'H'
                && if alpha == 0.0 {
                    prefer_r
                } else {
                    gamma < 1.125 * alpha
                }
            {
                return GEMMI_SPACEGROUPS.get(group_index + 1);
            }
            if group.ext == b'1' && prefer_2 {
                return GEMMI_SPACEGROUPS.get(group_index + 1);
            }
            return Some(group);
        }

        if input_tail == b':' {
            let extension_index = skip_space(input, input_index + 1, name.c_string_end);
            if byte_or_nul(input, extension_index, name.c_string_end) == group.ext {
                return Some(group);
            }
        }
    }

    None
}

// Gemmi❗✔️:   for (const SpaceGroupAltName& sg : spacegroup_tables::alt_names)
// Gemmi❗✔️:     if (sg.hm[0] == first && sg.hm[2] == *p) {
// Gemmi❗✔️:       const char* a = skip_space(p + 1);
// Gemmi❗✔️:       const char* b = skip_space(sg.hm + 3);
// Gemmi❗✔️:       while (*a == *b && *b != '\0') {
// Gemmi❗✔️:         a = skip_space(a+1);
// Gemmi❗✔️:         b = skip_space(b+1);
// Gemmi❗✔️:       }
// Gemmi❗✔️:       if (*b == '\0' &&
// Gemmi❗✔️:           (*a == '\0' || (*a == ':' && *skip_space(a+1) == sg.ext)))
// Gemmi❗✔️:         return &spacegroup_tables::main[sg.pos];
// Gemmi❗✔️:     }
// Gemmi❗✔️:   return nullptr;
// Behavior review: alternate entries run in source order after all main rows;
// the C-string terminator, extension check and target main-row index are
// retained. Step44/46 pinned cases pass for the tested aliases/extensions;
// the marker remains conservative for unenumerated alias/tail combinations.
// Complexity review: one fixed 28-row scan only on a main-table miss, with
// monotone cursors and no allocation per row.
fn find_alternate_spacegroup_by_name(
    name: &NormalizedSpaceGroupName,
) -> Option<&'static BioSpaceGroup> {
    let input = &name.bytes;
    for alternate in GEMMI_ALT_NAMES {
        let alternate_name = alternate.hm.as_bytes();
        if alternate_name[0] != name.first
            || alternate_name[2] != byte_or_nul(input, name.start, name.c_string_end)
        {
            continue;
        }

        let mut input_index = skip_space(input, name.start + 1, name.c_string_end);
        let mut alternate_index = skip_space(alternate_name, 3, alternate_name.len());
        while byte_or_nul(input, input_index, name.c_string_end)
            == byte_or_nul(alternate_name, alternate_index, alternate_name.len())
            && byte_or_nul(alternate_name, alternate_index, alternate_name.len()) != 0
        {
            input_index = skip_space(input, input_index + 1, name.c_string_end);
            alternate_index = skip_space(alternate_name, alternate_index + 1, alternate_name.len());
        }

        if byte_or_nul(alternate_name, alternate_index, alternate_name.len()) == 0 {
            let input_tail = byte_or_nul(input, input_index, name.c_string_end);
            if input_tail == 0 {
                return GEMMI_SPACEGROUPS.get(alternate.pos);
            }
            if input_tail == b':' {
                let extension_index = skip_space(input, input_index + 1, name.c_string_end);
                if byte_or_nul(input, extension_index, name.c_string_end) == alternate.ext {
                    return GEMMI_SPACEGROUPS.get(alternate.pos);
                }
            }
        }
    }

    None
}

// Gemmi❗✔️:   const char* p = skip_space(name.c_str());
// Gemmi❗✔️:   if (*p >= '0' && *p <= '9') { // handle numbers
// Gemmi❗✔️:     char *endptr;
// Gemmi❗✔️:     long n = std::strtol(p, &endptr, 10);
// Gemmi❗✔️:     return *endptr == '\0' ? find_spacegroup_by_number(n) : nullptr;
// Gemmi❗✔️:   }
// Behavior review: `parse_spacegroup_name` owns the pinned token/normalization
// rules. Numeric keys resolve by CCP4; named keys scan the complete ordered
// main table, then the alternate table, with the source caller's null
// preference. No public preference or facade API is introduced.
// Complexity review: the same linear CCP4 scan or ordered main scan plus an
// alternate scan only on a main miss; no tables are copied or sorted.
fn find_spacegroup_by_name(name: &str, alpha: f64, gamma: f64) -> Option<&'static BioSpaceGroup> {
    match parse_spacegroup_name(name)? {
        ParsedSpaceGroupName::Numeric { ccp4 } => find_spacegroup_by_number(ccp4),
        ParsedSpaceGroupName::Named(normalized) => {
            find_main_spacegroup_by_name(&normalized, alpha, gamma, false, false)
                .or_else(|| find_alternate_spacegroup_by_name(&normalized))
        }
    }
}

// Gemmi❗✔️:   const SpaceGroup* find_spacegroup() const {
// Gemmi❗✔️:     if (!cell.is_crystal())
// Gemmi❗✔️:       return nullptr;
// Gemmi❗✔️:     return find_spacegroup_by_name(spacegroup_hm, cell.alpha, cell.gamma);
// Gemmi❗✔️:   }
// Behavior review: `BioCrystalInfo::is_crystal` is the canonical BIO mapping
// of the source predicate. Absent name state maps to the source empty-string
// no-match result; this remains private to BIO and does not expose runtime.
// Complexity review: constant-time crystal/accessor checks followed by the
// same lookup; the crystal and static table are not cloned.
fn find_structure_spacegroup(crystal: &super::BioCrystalInfo) -> Option<&'static BioSpaceGroup> {
    if !crystal.is_crystal() {
        return None;
    }
    let cell = crystal.cell();
    find_spacegroup_by_name(
        crystal.space_group_hm().unwrap_or_default(),
        cell.alpha,
        cell.gamma,
    )
}

pub(super) fn structure_space_group_number(crystal: &super::BioCrystalInfo) -> Option<i32> {
    // Gemmi❗✔️: if (const SpaceGroup* sg = st.find_spacegroup())
    // Gemmi❗✔️:   span.set_pair("_symmetry.Int_Tables_number", std::to_string(sg->number));
    // Thin scalar projection of the unique source-backed lookup, with no copied table.
    find_structure_spacegroup(crystal).map(|group| group.number)
}

fn set_crystal_symmetry_images(crystal: &mut super::BioCrystalInfo) {
    // Gemmi❗✔️:   void set_cell_images_from_groupops(const GroupOps& group_ops) {
    // Gemmi❗✔️:     images.clear();
    // Gemmi❗✔️:     cs_count = (short) group_ops.order() - 1;
    // Gemmi❗✔️:     images.reserve(cs_count);
    // Gemmi❗✔️:     for (Op op : group_ops)
    // Gemmi❗✔️:       if (op != Op::identity())
    // Gemmi❗✔️:         images.push_back(Transform{rot_as_mat33(op), tran_as_vec3(op)});
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   void set_cell_images_from_spacegroup(const SpaceGroup* sg) {
    // Gemmi❗✔️:     if (sg) {
    // Gemmi❗✔️:       set_cell_images_from_groupops(sg->operations());
    // Gemmi❗✔️:     } else {
    // Gemmi❗✔️:       images.clear();
    // Gemmi❗✔️:       cs_count = 0;
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   GroupOps operations() const { return symops_from_hall(hall); }
    // Gemmi❗✔️:   int order() const { return static_cast<int>(sym_ops.size()*cen_ops.size()); }
    // Gemmi❗✔️:   static constexpr Op identity() {
    // Gemmi❗✔️:     return {{{{DEN,0,0}, {0,DEN,0}, {0,0,DEN}}}, {0,0,0}, ' '};
    // Gemmi❗✔️:   }
    // Gemmi❗✔️: inline bool operator==(const Op& a, const Op& b) {
    // Gemmi❗✔️:   return a.rot == b.rot && a.tran == b.tran;
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline Mat33 rot_as_mat33(const Op::Rot& rot) {
    // Gemmi❗✔️:   double mult = 1.0 / Op::DEN;
    // Gemmi❗✔️:   return Mat33(mult * rot[0][0], mult * rot[0][1], mult * rot[0][2],
    // Gemmi❗✔️:                mult * rot[1][0], mult * rot[1][1], mult * rot[1][2],
    // Gemmi❗✔️:                mult * rot[2][0], mult * rot[2][1], mult * rot[2][2]);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline Mat33 rot_as_mat33(const Op& op) { return rot_as_mat33(op.rot); }
    // Gemmi❗✔️: inline Vec3 tran_as_vec3(const Op& op) {
    // Gemmi❗✔️:   double mult = 1.0 / Op::DEN;
    // Gemmi❗✔️:   return Vec3(mult * op.tran[0], mult * op.tran[1], mult * op.tran[2]);
    // Gemmi❗✔️: }
    //
    // Behavior review: `find_structure_spacegroup` is the pinned Structure
    // crystal/name gate. `BioSpaceGroup::ops` is the complete, ordered result
    // of `SpaceGroup::operations()`; Step22 verified every row and all 7,740
    // operation tuples against that source. Preserve its source order, count
    // formula, explicit equality-based identity filter, and reciprocal-denominator
    // multiplication order. No match clears both stored fields.
    // Complexity review: one O(number of operations) pass, clearing/reusing
    // the existing vector allocation and producing one fixed-size transform
    // per nonidentity operation, matching Gemmi's operation loop.
    let Some(space_group) = find_structure_spacegroup(crystal) else {
        crystal.symmetry_images.clear();
        crystal.cs_count = 0;
        return;
    };

    crystal.symmetry_images.clear();
    let operation_count = space_group.ops.len();
    let cs_count = i16::try_from(operation_count - 1)
        .expect("verified pinned Gemmi operation counts fit BioCrystalInfo::cs_count");
    crystal.cs_count = cs_count;
    crystal
        .symmetry_images
        .reserve(usize::try_from(cs_count).expect("space-group operation count is nonnegative"));

    let identity = BioSymOp {
        rot: [
            [GEMMI_OP_DEN, 0, 0],
            [0, GEMMI_OP_DEN, 0],
            [0, 0, GEMMI_OP_DEN],
        ],
        tran: [0, 0, 0],
    };
    let multiplier = 1.0 / f64::from(GEMMI_OP_DEN);
    for operation in space_group.ops.iter().copied() {
        if operation == identity {
            continue;
        }
        let rotation = operation
            .rot
            .map(|row| row.map(|component| multiplier * f64::from(component)));
        let translation = operation
            .tran
            .map(|component| multiplier * f64::from(component));
        crystal
            .symmetry_images
            .push(super::BioTransform::new(rotation, translation));
    }
}

fn add_ncs_images_to_crystal(crystal: &mut super::BioCrystalInfo, ncs: &[super::BioNcsOperator]) {
    // Gemmi❗✔️:   void add_ncs_images_to_cs_images(const std::vector<NcsOp>& ncs) {
    // Gemmi❗✔️:     assert(cs_count == (short) images.size());
    // Gemmi❗✔️:     for (const NcsOp& ncs_op : ncs)
    // Gemmi❗✔️:       if (!ncs_op.given) {
    // Gemmi❗✔️:         // We need it to operates on fractional, not orthogonal coordinates.
    // Gemmi❗✔️:         FTransform f = frac.combine(ncs_op.tr.combine(orth));
    // Gemmi❗✔️:         images.push_back(f);
    // Gemmi❗✔️:         for (int i = 0; i < cs_count; ++i)
    // Gemmi❗✔️:           images.push_back(images[i].combine(f));
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:   }
    //
    // Behavior review: the source's signed `i < cs_count` loop ranges only
    // over the original crystallographic images. Every nongiven NCS transform
    // is appended first, followed by its ordered compositions with that same
    // original prefix; `cs_count` is not changed. The nested composition order
    // is exactly fractional ∘ (NCS ∘ orthogonal).
    // Complexity review: O(ncs * (1 + cs_count)) fixed-size transform
    // compositions and vector appends, matching Gemmi's loops and storage.
    debug_assert_eq!(crystal.cs_count, crystal.symmetry_images.len() as i16);
    for ncs_operator in ncs {
        if !ncs_operator.given {
            let fractional_ncs = crystal
                .fractional
                .combine(&ncs_operator.transform.combine(&crystal.orthogonal));
            crystal.symmetry_images.push(fractional_ncs);
            for index in 0..crystal.cs_count {
                let crystal_ncs_image =
                    crystal.symmetry_images[index as usize].combine(&fractional_ncs);
                crystal.symmetry_images.push(crystal_ncs_image);
            }
        }
    }
}

pub fn setup_cell_images(crystal: &mut super::BioCrystalInfo, ncs: &[super::BioNcsOperator]) {
    // Gemmi❗✔️: inline void Structure::setup_cell_images() {
    // Gemmi❗✔️:   const SpaceGroup* sg = find_spacegroup();
    // Gemmi❗✔️:   cell.set_cell_images_from_spacegroup(sg);
    // Gemmi❗✔️:   cell.add_ncs_images_to_cs_images(ncs);
    // Gemmi❗✔️: }
    //
    // Behavior review: the first helper performs the pinned crystal/name
    // lookup and replaces the ordered crystallographic image state; the next
    // helper expands the supplied NCS rows from that prepared state. This
    // preserves Gemmi's call order without moving parser or hierarchy state
    // into BIO.
    // Complexity review: this wrapper adds no traversal, allocation or copy;
    // the two called owners retain their source-shaped table and NCS loops.
    set_crystal_symmetry_images(crystal);
    add_ncs_images_to_crystal(crystal, ncs);
}

fn byte_or_nul(bytes: &[u8], index: usize, end: usize) -> u8 {
    // Gemmi❗✔️: if (*b == '\0') {
    // Gemmi❗✔️:   if (*a == '\0') {
    // Gemmi❗✔️:     return &sg;
    // Gemmi❗✔️:   }
    // Behavior review: the C++ cursors dereference their terminating NUL;
    // this bounded byte accessor supplies that same sentinel when an index
    // reaches the Rust slice's source-string end.
    // Complexity review: constant-time indexed access, with no allocation.
    if index < end { bytes[index] } else { 0 }
}

// Gemmi✔️✔️: GroupOps operations() const { return symops_from_hall(hall); }
// Exact source rows preserve source number/ccp4 identifiers and name metadata.
// The fixed operation vectors use denominator 24 and retain source order.
const OPS_0: &[BioSymOp] = &[BioSymOp {
    rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
    tran: [0, 0, 0],
}];

const OPS_1: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_2: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_3: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_4: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_5: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_6: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_7: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
];

const OPS_8: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_9: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_10: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_11: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_12: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_13: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_14: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_15: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_16: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_17: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_18: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_19: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_20: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_21: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_22: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_23: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
];

const OPS_24: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_25: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_26: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_27: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_28: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_29: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_30: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_31: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_32: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_33: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_34: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_35: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_36: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_37: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_38: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_39: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_40: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_41: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_42: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_43: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_44: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_45: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_46: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_47: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_48: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_49: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_50: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_51: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_52: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_53: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_54: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_55: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_56: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_57: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_58: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_59: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_60: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
];

const OPS_61: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_62: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_63: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_64: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_65: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_66: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_67: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_68: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_69: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_70: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_71: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_72: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_73: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_74: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
];

const OPS_75: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_76: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_77: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_78: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_79: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_80: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_81: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_82: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_83: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_84: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_85: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_86: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_87: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_88: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_89: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_90: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_91: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_92: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_93: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_94: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_95: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_96: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_97: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_98: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_99: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_100: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_101: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_102: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_103: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_104: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_105: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_106: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_107: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_108: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
];

const OPS_109: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_110: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_111: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_112: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_113: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_114: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_115: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_116: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_117: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_118: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_119: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_120: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_121: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_122: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_123: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
];

const OPS_124: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_125: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_126: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_127: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_128: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_129: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_130: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_131: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_132: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_133: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_134: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_135: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_136: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_137: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_138: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_139: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_140: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_141: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_142: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_143: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_144: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_145: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_146: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_147: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_148: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_149: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_150: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_151: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_152: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_153: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_154: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_155: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_156: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_157: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_158: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_159: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_160: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_161: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_162: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_163: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_164: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_165: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_166: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_167: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_168: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_169: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_170: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_171: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_172: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_173: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_174: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_175: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_176: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_177: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_178: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_179: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_180: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_181: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_182: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_183: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_184: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_185: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_186: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_187: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_188: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_189: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_190: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_191: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_192: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_193: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_194: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_195: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_196: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_197: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_198: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_199: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_200: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_201: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_202: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_203: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_204: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_205: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_206: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_207: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_208: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_209: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_210: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_211: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
];

const OPS_212: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
];

const OPS_213: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_214: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_215: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_216: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_217: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_218: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_219: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_220: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_221: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_222: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_223: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_224: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_225: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_226: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_227: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_228: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_229: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_230: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_231: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_232: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_233: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_234: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_235: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_236: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_237: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_238: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_239: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_240: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_241: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_242: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_243: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_244: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_245: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_246: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_247: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_248: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_249: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_250: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_251: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_252: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_253: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_254: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_255: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_256: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_257: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_258: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_259: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_260: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_261: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_262: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_263: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_264: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_265: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_266: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_267: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_268: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_269: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_270: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_271: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_272: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_273: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_274: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_275: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_276: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_277: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_278: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_279: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_280: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_281: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_282: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_283: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_284: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_285: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_286: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_287: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_288: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_289: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_290: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_291: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_292: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_293: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_294: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_295: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_296: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_297: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_298: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_299: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_300: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_301: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_302: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_303: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_304: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_305: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_306: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_307: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_308: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_309: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_310: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_311: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_312: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_313: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_314: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_315: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_316: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_317: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_318: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_319: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_320: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_321: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_322: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_323: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_324: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_325: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_326: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_327: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_328: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_329: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_330: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_331: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_332: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_333: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_334: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
];

const OPS_335: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 12, 18],
    },
];

const OPS_336: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_337: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_338: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_339: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_340: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
];

const OPS_341: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_342: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_343: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_344: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_345: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_346: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_347: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_348: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_349: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 18],
    },
];

const OPS_350: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_351: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 6],
    },
];

const OPS_352: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_353: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
];

const OPS_354: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_355: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_356: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_357: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
];

const OPS_358: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_359: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_360: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_361: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_362: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_363: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_364: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 18],
    },
];

const OPS_365: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_366: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_367: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 18],
    },
];

const OPS_368: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_369: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
];

const OPS_370: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_371: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 6],
    },
];

const OPS_372: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_373: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_374: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_375: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_376: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_377: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_378: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_379: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_380: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_381: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_382: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_383: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_384: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_385: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
];

const OPS_386: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 18],
    },
];

const OPS_387: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_388: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_389: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_390: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_391: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_392: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
];

const OPS_393: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_394: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_395: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_396: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_397: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_398: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
];

const OPS_399: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_400: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_401: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_402: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_403: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_404: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_405: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_406: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_407: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_408: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_409: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_410: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_411: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_412: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_413: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_414: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_415: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_416: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_417: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_418: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_419: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_420: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_421: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_422: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_423: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_424: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_425: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
];

const OPS_426: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 6],
    },
];

const OPS_427: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 6],
    },
];

const OPS_428: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
];

const OPS_429: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_430: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
];

const OPS_431: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
];

const OPS_432: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
];

const OPS_433: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_434: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_435: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
];

const OPS_436: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_437: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_438: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_439: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_440: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
];

const OPS_441: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_442: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
];

const OPS_443: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
];

const OPS_444: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_445: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_446: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_447: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_448: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_449: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
];

const OPS_450: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_451: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 20],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 20],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 20],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 4],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 4],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 4],
    },
];

const OPS_452: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_453: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_454: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_455: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_456: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_457: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
];

const OPS_458: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_459: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [16, 8, 20],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [16, 8, 20],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [16, 8, 20],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [16, 8, 8],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [16, 8, 20],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [16, 8, 20],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [16, 8, 20],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [8, 16, 4],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [8, 16, 4],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [8, 16, 4],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [8, 16, 16],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [8, 16, 4],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [8, 16, 4],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [8, 16, 4],
    },
];

const OPS_460: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_461: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_462: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 4],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 20],
    },
];

const OPS_463: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 20],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 4],
    },
];

const OPS_464: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
];

const OPS_465: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
];

const OPS_466: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_467: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_468: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_469: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
];

const OPS_470: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_471: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 4],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 20],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 20],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 4],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_472: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 20],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 4],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 4],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 20],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_473: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_474: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 8],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 16],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_475: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_476: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_477: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_478: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_479: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_480: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_481: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_482: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_483: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
];

const OPS_484: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_485: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_486: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
];

const OPS_487: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [24, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, -24, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [-24, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 24, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
];

const OPS_488: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_489: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
];

const OPS_490: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_491: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
];

const OPS_492: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
];

const OPS_493: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_494: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_495: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
];

const OPS_496: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
];

const OPS_497: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [18, 18, 6],
    },
];

const OPS_498: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 6, 18],
    },
];

const OPS_499: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_500: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
];

const OPS_501: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 12],
    },
];

const OPS_502: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_503: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_504: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
];

const OPS_505: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_506: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_507: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
];

const OPS_508: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
];

const OPS_509: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
];

const OPS_510: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_511: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
];

const OPS_512: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_513: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_514: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
];

const OPS_515: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
];

const OPS_516: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_517: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_518: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
];

const OPS_519: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
];

const OPS_520: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_521: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
];

const OPS_522: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
];

const OPS_523: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
];

const OPS_524: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
];

const OPS_525: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 18, 12],
    },
];

const OPS_526: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 6, 6],
    },
];

const OPS_527: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 12, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 12, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 18, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 18, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 0, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 12, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 12, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 6, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 6, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 0, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 0, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 0, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 6, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 6, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 12, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 0, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 0, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 18, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 18, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 12, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 12, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 12, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 18, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 18, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 0, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 12, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 12, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 6, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 18, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 6, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 0, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 0, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 0, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 6, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 6, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 12, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [6, 6, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [6, 18, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 0, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 0, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 18, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [18, 6, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 18, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [18, 12, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 12, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [18, 18, 0],
    },
];

const OPS_528: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
];

const OPS_529: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [18, 18, 18],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [6, 18, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, 24, 0], [-24, 0, 0]],
        tran: [6, 6, 18],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [0, -24, 0], [24, 0, 0]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, 24, 0], [24, 0, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [0, -24, 0], [-24, 0, 0]],
        tran: [18, 6, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, -24], [0, -24, 0]],
        tran: [6, 18, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, -24], [0, 24, 0]],
        tran: [18, 6, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 0, 24], [0, 24, 0]],
        tran: [6, 6, 6],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 0, 24], [0, -24, 0]],
        tran: [18, 18, 6],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
];

const OPS_530: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_531: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
];

const OPS_532: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_533: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_534: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_535: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_536: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
];

const OPS_537: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 12],
    },
];

const OPS_538: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, 0, 24], [-24, 0, 0], [0, -24, 0]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [24, 0, 0], [0, -24, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 0, -24], [-24, 0, 0], [0, 24, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, 24], [-24, 0, 0]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [0, 0, -24], [24, 0, 0]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [0, 0, -24], [-24, 0, 0]],
        tran: [12, 0, 0],
    },
];

const OPS_539: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
];

const OPS_540: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_541: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_542: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_543: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_544: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
];

const OPS_545: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_546: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_547: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_548: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_549: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
];

const OPS_550: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_551: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 12],
    },
];

const OPS_552: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 12],
    },
];

const OPS_553: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_554: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_555: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 0, 6],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [6, 12, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 0, 18],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [18, 12, 6],
    },
];

const OPS_556: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_557: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
];

const OPS_558: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_559: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
];

const OPS_560: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
];

const OPS_561: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

const OPS_562: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 0],
    },
];

const OPS_563: &[BioSymOp] = &[
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 0, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [0, 12, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 0, 12],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, -24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [-24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, 24, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [24, 0, 0], [0, 0, -24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[-24, 0, 0], [0, 24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, 24, 0], [24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[24, 0, 0], [0, -24, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
    BioSymOp {
        rot: [[0, -24, 0], [-24, 0, 0], [0, 0, 24]],
        tran: [12, 12, 0],
    },
];

// Gemmi✔️✔️: const SpaceGroup spacegroup_tables::main[564] = {
static GEMMI_SPACEGROUPS: &[BioSpaceGroup] = &[
    BioSpaceGroup {
        number: 1,
        ccp4: 1,
        hm: "P 1",
        ext: 0,
        qualifier: "",
        hall: "P 1",
        basisop_idx: 0,
        ops: OPS_0,
    },
    BioSpaceGroup {
        number: 2,
        ccp4: 2,
        hm: "P -1",
        ext: 0,
        qualifier: "",
        hall: "-P 1",
        basisop_idx: 0,
        ops: OPS_1,
    },
    BioSpaceGroup {
        number: 3,
        ccp4: 3,
        hm: "P 1 2 1",
        ext: 0,
        qualifier: "b",
        hall: "P 2y",
        basisop_idx: 0,
        ops: OPS_2,
    },
    BioSpaceGroup {
        number: 3,
        ccp4: 1003,
        hm: "P 1 1 2",
        ext: 0,
        qualifier: "c",
        hall: "P 2",
        basisop_idx: 1,
        ops: OPS_3,
    },
    BioSpaceGroup {
        number: 3,
        ccp4: 0,
        hm: "P 2 1 1",
        ext: 0,
        qualifier: "a",
        hall: "P 2x",
        basisop_idx: 2,
        ops: OPS_4,
    },
    BioSpaceGroup {
        number: 4,
        ccp4: 4,
        hm: "P 1 21 1",
        ext: 0,
        qualifier: "b",
        hall: "P 2yb",
        basisop_idx: 0,
        ops: OPS_5,
    },
    BioSpaceGroup {
        number: 4,
        ccp4: 1004,
        hm: "P 1 1 21",
        ext: 0,
        qualifier: "c",
        hall: "P 2c",
        basisop_idx: 1,
        ops: OPS_6,
    },
    BioSpaceGroup {
        number: 4,
        ccp4: 0,
        hm: "P 21 1 1",
        ext: 0,
        qualifier: "a",
        hall: "P 2xa",
        basisop_idx: 2,
        ops: OPS_7,
    },
    BioSpaceGroup {
        number: 5,
        ccp4: 5,
        hm: "C 1 2 1",
        ext: 0,
        qualifier: "b1",
        hall: "C 2y",
        basisop_idx: 0,
        ops: OPS_8,
    },
    BioSpaceGroup {
        number: 5,
        ccp4: 2005,
        hm: "A 1 2 1",
        ext: 0,
        qualifier: "b2",
        hall: "A 2y",
        basisop_idx: 3,
        ops: OPS_9,
    },
    BioSpaceGroup {
        number: 5,
        ccp4: 4005,
        hm: "I 1 2 1",
        ext: 0,
        qualifier: "b3",
        hall: "I 2y",
        basisop_idx: 4,
        ops: OPS_10,
    },
    BioSpaceGroup {
        number: 5,
        ccp4: 0,
        hm: "A 1 1 2",
        ext: 0,
        qualifier: "c1",
        hall: "A 2",
        basisop_idx: 1,
        ops: OPS_11,
    },
    BioSpaceGroup {
        number: 5,
        ccp4: 1005,
        hm: "B 1 1 2",
        ext: 0,
        qualifier: "c2",
        hall: "B 2",
        basisop_idx: 5,
        ops: OPS_12,
    },
    BioSpaceGroup {
        number: 5,
        ccp4: 0,
        hm: "I 1 1 2",
        ext: 0,
        qualifier: "c3",
        hall: "I 2",
        basisop_idx: 6,
        ops: OPS_13,
    },
    BioSpaceGroup {
        number: 5,
        ccp4: 0,
        hm: "B 2 1 1",
        ext: 0,
        qualifier: "a1",
        hall: "B 2x",
        basisop_idx: 2,
        ops: OPS_14,
    },
    BioSpaceGroup {
        number: 5,
        ccp4: 0,
        hm: "C 2 1 1",
        ext: 0,
        qualifier: "a2",
        hall: "C 2x",
        basisop_idx: 7,
        ops: OPS_15,
    },
    BioSpaceGroup {
        number: 5,
        ccp4: 0,
        hm: "I 2 1 1",
        ext: 0,
        qualifier: "a3",
        hall: "I 2x",
        basisop_idx: 8,
        ops: OPS_16,
    },
    BioSpaceGroup {
        number: 6,
        ccp4: 6,
        hm: "P 1 m 1",
        ext: 0,
        qualifier: "b",
        hall: "P -2y",
        basisop_idx: 0,
        ops: OPS_17,
    },
    BioSpaceGroup {
        number: 6,
        ccp4: 1006,
        hm: "P 1 1 m",
        ext: 0,
        qualifier: "c",
        hall: "P -2",
        basisop_idx: 1,
        ops: OPS_18,
    },
    BioSpaceGroup {
        number: 6,
        ccp4: 0,
        hm: "P m 1 1",
        ext: 0,
        qualifier: "a",
        hall: "P -2x",
        basisop_idx: 2,
        ops: OPS_19,
    },
    BioSpaceGroup {
        number: 7,
        ccp4: 7,
        hm: "P 1 c 1",
        ext: 0,
        qualifier: "b1",
        hall: "P -2yc",
        basisop_idx: 0,
        ops: OPS_20,
    },
    BioSpaceGroup {
        number: 7,
        ccp4: 0,
        hm: "P 1 n 1",
        ext: 0,
        qualifier: "b2",
        hall: "P -2yac",
        basisop_idx: 9,
        ops: OPS_21,
    },
    BioSpaceGroup {
        number: 7,
        ccp4: 0,
        hm: "P 1 a 1",
        ext: 0,
        qualifier: "b3",
        hall: "P -2ya",
        basisop_idx: 3,
        ops: OPS_22,
    },
    BioSpaceGroup {
        number: 7,
        ccp4: 0,
        hm: "P 1 1 a",
        ext: 0,
        qualifier: "c1",
        hall: "P -2a",
        basisop_idx: 1,
        ops: OPS_23,
    },
    BioSpaceGroup {
        number: 7,
        ccp4: 0,
        hm: "P 1 1 n",
        ext: 0,
        qualifier: "c2",
        hall: "P -2ab",
        basisop_idx: 10,
        ops: OPS_24,
    },
    BioSpaceGroup {
        number: 7,
        ccp4: 1007,
        hm: "P 1 1 b",
        ext: 0,
        qualifier: "c3",
        hall: "P -2b",
        basisop_idx: 5,
        ops: OPS_25,
    },
    BioSpaceGroup {
        number: 7,
        ccp4: 0,
        hm: "P b 1 1",
        ext: 0,
        qualifier: "a1",
        hall: "P -2xb",
        basisop_idx: 2,
        ops: OPS_26,
    },
    BioSpaceGroup {
        number: 7,
        ccp4: 0,
        hm: "P n 1 1",
        ext: 0,
        qualifier: "a2",
        hall: "P -2xbc",
        basisop_idx: 11,
        ops: OPS_27,
    },
    BioSpaceGroup {
        number: 7,
        ccp4: 0,
        hm: "P c 1 1",
        ext: 0,
        qualifier: "a3",
        hall: "P -2xc",
        basisop_idx: 7,
        ops: OPS_28,
    },
    BioSpaceGroup {
        number: 8,
        ccp4: 8,
        hm: "C 1 m 1",
        ext: 0,
        qualifier: "b1",
        hall: "C -2y",
        basisop_idx: 0,
        ops: OPS_29,
    },
    BioSpaceGroup {
        number: 8,
        ccp4: 0,
        hm: "A 1 m 1",
        ext: 0,
        qualifier: "b2",
        hall: "A -2y",
        basisop_idx: 3,
        ops: OPS_30,
    },
    BioSpaceGroup {
        number: 8,
        ccp4: 0,
        hm: "I 1 m 1",
        ext: 0,
        qualifier: "b3",
        hall: "I -2y",
        basisop_idx: 4,
        ops: OPS_31,
    },
    BioSpaceGroup {
        number: 8,
        ccp4: 0,
        hm: "A 1 1 m",
        ext: 0,
        qualifier: "c1",
        hall: "A -2",
        basisop_idx: 1,
        ops: OPS_32,
    },
    BioSpaceGroup {
        number: 8,
        ccp4: 1008,
        hm: "B 1 1 m",
        ext: 0,
        qualifier: "c2",
        hall: "B -2",
        basisop_idx: 5,
        ops: OPS_33,
    },
    BioSpaceGroup {
        number: 8,
        ccp4: 0,
        hm: "I 1 1 m",
        ext: 0,
        qualifier: "c3",
        hall: "I -2",
        basisop_idx: 6,
        ops: OPS_34,
    },
    BioSpaceGroup {
        number: 8,
        ccp4: 0,
        hm: "B m 1 1",
        ext: 0,
        qualifier: "a1",
        hall: "B -2x",
        basisop_idx: 2,
        ops: OPS_35,
    },
    BioSpaceGroup {
        number: 8,
        ccp4: 0,
        hm: "C m 1 1",
        ext: 0,
        qualifier: "a2",
        hall: "C -2x",
        basisop_idx: 7,
        ops: OPS_36,
    },
    BioSpaceGroup {
        number: 8,
        ccp4: 0,
        hm: "I m 1 1",
        ext: 0,
        qualifier: "a3",
        hall: "I -2x",
        basisop_idx: 8,
        ops: OPS_37,
    },
    BioSpaceGroup {
        number: 9,
        ccp4: 9,
        hm: "C 1 c 1",
        ext: 0,
        qualifier: "b1",
        hall: "C -2yc",
        basisop_idx: 0,
        ops: OPS_38,
    },
    BioSpaceGroup {
        number: 9,
        ccp4: 0,
        hm: "A 1 n 1",
        ext: 0,
        qualifier: "b2",
        hall: "A -2yab",
        basisop_idx: 12,
        ops: OPS_39,
    },
    BioSpaceGroup {
        number: 9,
        ccp4: 0,
        hm: "I 1 a 1",
        ext: 0,
        qualifier: "b3",
        hall: "I -2ya",
        basisop_idx: 13,
        ops: OPS_40,
    },
    BioSpaceGroup {
        number: 9,
        ccp4: 0,
        hm: "A 1 a 1",
        ext: 0,
        qualifier: "-b1",
        hall: "A -2ya",
        basisop_idx: 3,
        ops: OPS_41,
    },
    BioSpaceGroup {
        number: 9,
        ccp4: 0,
        hm: "C 1 n 1",
        ext: 0,
        qualifier: "-b2",
        hall: "C -2yac",
        basisop_idx: 14,
        ops: OPS_42,
    },
    BioSpaceGroup {
        number: 9,
        ccp4: 0,
        hm: "I 1 c 1",
        ext: 0,
        qualifier: "-b3",
        hall: "I -2yc",
        basisop_idx: 4,
        ops: OPS_43,
    },
    BioSpaceGroup {
        number: 9,
        ccp4: 0,
        hm: "A 1 1 a",
        ext: 0,
        qualifier: "c1",
        hall: "A -2a",
        basisop_idx: 1,
        ops: OPS_44,
    },
    BioSpaceGroup {
        number: 9,
        ccp4: 0,
        hm: "B 1 1 n",
        ext: 0,
        qualifier: "c2",
        hall: "B -2ab",
        basisop_idx: 15,
        ops: OPS_45,
    },
    BioSpaceGroup {
        number: 9,
        ccp4: 0,
        hm: "I 1 1 b",
        ext: 0,
        qualifier: "c3",
        hall: "I -2b",
        basisop_idx: 16,
        ops: OPS_46,
    },
    BioSpaceGroup {
        number: 9,
        ccp4: 1009,
        hm: "B 1 1 b",
        ext: 0,
        qualifier: "-c1",
        hall: "B -2b",
        basisop_idx: 5,
        ops: OPS_47,
    },
    BioSpaceGroup {
        number: 9,
        ccp4: 0,
        hm: "A 1 1 n",
        ext: 0,
        qualifier: "-c2",
        hall: "A -2ab",
        basisop_idx: 10,
        ops: OPS_48,
    },
    BioSpaceGroup {
        number: 9,
        ccp4: 0,
        hm: "I 1 1 a",
        ext: 0,
        qualifier: "-c3",
        hall: "I -2a",
        basisop_idx: 6,
        ops: OPS_49,
    },
    BioSpaceGroup {
        number: 9,
        ccp4: 0,
        hm: "B b 1 1",
        ext: 0,
        qualifier: "a1",
        hall: "B -2xb",
        basisop_idx: 2,
        ops: OPS_50,
    },
    BioSpaceGroup {
        number: 9,
        ccp4: 0,
        hm: "C n 1 1",
        ext: 0,
        qualifier: "a2",
        hall: "C -2xac",
        basisop_idx: 17,
        ops: OPS_51,
    },
    BioSpaceGroup {
        number: 9,
        ccp4: 0,
        hm: "I c 1 1",
        ext: 0,
        qualifier: "a3",
        hall: "I -2xc",
        basisop_idx: 18,
        ops: OPS_52,
    },
    BioSpaceGroup {
        number: 9,
        ccp4: 0,
        hm: "C c 1 1",
        ext: 0,
        qualifier: "-a1",
        hall: "C -2xc",
        basisop_idx: 7,
        ops: OPS_53,
    },
    BioSpaceGroup {
        number: 9,
        ccp4: 0,
        hm: "B n 1 1",
        ext: 0,
        qualifier: "-a2",
        hall: "B -2xab",
        basisop_idx: 11,
        ops: OPS_54,
    },
    BioSpaceGroup {
        number: 9,
        ccp4: 0,
        hm: "I b 1 1",
        ext: 0,
        qualifier: "-a3",
        hall: "I -2xb",
        basisop_idx: 8,
        ops: OPS_55,
    },
    BioSpaceGroup {
        number: 10,
        ccp4: 10,
        hm: "P 1 2/m 1",
        ext: 0,
        qualifier: "b",
        hall: "-P 2y",
        basisop_idx: 0,
        ops: OPS_56,
    },
    BioSpaceGroup {
        number: 10,
        ccp4: 1010,
        hm: "P 1 1 2/m",
        ext: 0,
        qualifier: "c",
        hall: "-P 2",
        basisop_idx: 1,
        ops: OPS_57,
    },
    BioSpaceGroup {
        number: 10,
        ccp4: 0,
        hm: "P 2/m 1 1",
        ext: 0,
        qualifier: "a",
        hall: "-P 2x",
        basisop_idx: 2,
        ops: OPS_58,
    },
    BioSpaceGroup {
        number: 11,
        ccp4: 11,
        hm: "P 1 21/m 1",
        ext: 0,
        qualifier: "b",
        hall: "-P 2yb",
        basisop_idx: 0,
        ops: OPS_59,
    },
    BioSpaceGroup {
        number: 11,
        ccp4: 1011,
        hm: "P 1 1 21/m",
        ext: 0,
        qualifier: "c",
        hall: "-P 2c",
        basisop_idx: 1,
        ops: OPS_60,
    },
    BioSpaceGroup {
        number: 11,
        ccp4: 0,
        hm: "P 21/m 1 1",
        ext: 0,
        qualifier: "a",
        hall: "-P 2xa",
        basisop_idx: 2,
        ops: OPS_61,
    },
    BioSpaceGroup {
        number: 12,
        ccp4: 12,
        hm: "C 1 2/m 1",
        ext: 0,
        qualifier: "b1",
        hall: "-C 2y",
        basisop_idx: 0,
        ops: OPS_62,
    },
    BioSpaceGroup {
        number: 12,
        ccp4: 0,
        hm: "A 1 2/m 1",
        ext: 0,
        qualifier: "b2",
        hall: "-A 2y",
        basisop_idx: 3,
        ops: OPS_63,
    },
    BioSpaceGroup {
        number: 12,
        ccp4: 0,
        hm: "I 1 2/m 1",
        ext: 0,
        qualifier: "b3",
        hall: "-I 2y",
        basisop_idx: 4,
        ops: OPS_64,
    },
    BioSpaceGroup {
        number: 12,
        ccp4: 0,
        hm: "A 1 1 2/m",
        ext: 0,
        qualifier: "c1",
        hall: "-A 2",
        basisop_idx: 1,
        ops: OPS_65,
    },
    BioSpaceGroup {
        number: 12,
        ccp4: 1012,
        hm: "B 1 1 2/m",
        ext: 0,
        qualifier: "c2",
        hall: "-B 2",
        basisop_idx: 5,
        ops: OPS_66,
    },
    BioSpaceGroup {
        number: 12,
        ccp4: 0,
        hm: "I 1 1 2/m",
        ext: 0,
        qualifier: "c3",
        hall: "-I 2",
        basisop_idx: 6,
        ops: OPS_67,
    },
    BioSpaceGroup {
        number: 12,
        ccp4: 0,
        hm: "B 2/m 1 1",
        ext: 0,
        qualifier: "a1",
        hall: "-B 2x",
        basisop_idx: 2,
        ops: OPS_68,
    },
    BioSpaceGroup {
        number: 12,
        ccp4: 0,
        hm: "C 2/m 1 1",
        ext: 0,
        qualifier: "a2",
        hall: "-C 2x",
        basisop_idx: 7,
        ops: OPS_69,
    },
    BioSpaceGroup {
        number: 12,
        ccp4: 0,
        hm: "I 2/m 1 1",
        ext: 0,
        qualifier: "a3",
        hall: "-I 2x",
        basisop_idx: 8,
        ops: OPS_70,
    },
    BioSpaceGroup {
        number: 13,
        ccp4: 13,
        hm: "P 1 2/c 1",
        ext: 0,
        qualifier: "b1",
        hall: "-P 2yc",
        basisop_idx: 0,
        ops: OPS_71,
    },
    BioSpaceGroup {
        number: 13,
        ccp4: 0,
        hm: "P 1 2/n 1",
        ext: 0,
        qualifier: "b2",
        hall: "-P 2yac",
        basisop_idx: 9,
        ops: OPS_72,
    },
    BioSpaceGroup {
        number: 13,
        ccp4: 0,
        hm: "P 1 2/a 1",
        ext: 0,
        qualifier: "b3",
        hall: "-P 2ya",
        basisop_idx: 3,
        ops: OPS_73,
    },
    BioSpaceGroup {
        number: 13,
        ccp4: 0,
        hm: "P 1 1 2/a",
        ext: 0,
        qualifier: "c1",
        hall: "-P 2a",
        basisop_idx: 1,
        ops: OPS_74,
    },
    BioSpaceGroup {
        number: 13,
        ccp4: 0,
        hm: "P 1 1 2/n",
        ext: 0,
        qualifier: "c2",
        hall: "-P 2ab",
        basisop_idx: 10,
        ops: OPS_75,
    },
    BioSpaceGroup {
        number: 13,
        ccp4: 1013,
        hm: "P 1 1 2/b",
        ext: 0,
        qualifier: "c3",
        hall: "-P 2b",
        basisop_idx: 5,
        ops: OPS_76,
    },
    BioSpaceGroup {
        number: 13,
        ccp4: 0,
        hm: "P 2/b 1 1",
        ext: 0,
        qualifier: "a1",
        hall: "-P 2xb",
        basisop_idx: 2,
        ops: OPS_77,
    },
    BioSpaceGroup {
        number: 13,
        ccp4: 0,
        hm: "P 2/n 1 1",
        ext: 0,
        qualifier: "a2",
        hall: "-P 2xbc",
        basisop_idx: 11,
        ops: OPS_78,
    },
    BioSpaceGroup {
        number: 13,
        ccp4: 0,
        hm: "P 2/c 1 1",
        ext: 0,
        qualifier: "a3",
        hall: "-P 2xc",
        basisop_idx: 7,
        ops: OPS_79,
    },
    BioSpaceGroup {
        number: 14,
        ccp4: 14,
        hm: "P 1 21/c 1",
        ext: 0,
        qualifier: "b1",
        hall: "-P 2ybc",
        basisop_idx: 0,
        ops: OPS_80,
    },
    BioSpaceGroup {
        number: 14,
        ccp4: 2014,
        hm: "P 1 21/n 1",
        ext: 0,
        qualifier: "b2",
        hall: "-P 2yn",
        basisop_idx: 9,
        ops: OPS_81,
    },
    BioSpaceGroup {
        number: 14,
        ccp4: 3014,
        hm: "P 1 21/a 1",
        ext: 0,
        qualifier: "b3",
        hall: "-P 2yab",
        basisop_idx: 3,
        ops: OPS_82,
    },
    BioSpaceGroup {
        number: 14,
        ccp4: 0,
        hm: "P 1 1 21/a",
        ext: 0,
        qualifier: "c1",
        hall: "-P 2ac",
        basisop_idx: 1,
        ops: OPS_83,
    },
    BioSpaceGroup {
        number: 14,
        ccp4: 0,
        hm: "P 1 1 21/n",
        ext: 0,
        qualifier: "c2",
        hall: "-P 2n",
        basisop_idx: 10,
        ops: OPS_84,
    },
    BioSpaceGroup {
        number: 14,
        ccp4: 1014,
        hm: "P 1 1 21/b",
        ext: 0,
        qualifier: "c3",
        hall: "-P 2bc",
        basisop_idx: 5,
        ops: OPS_85,
    },
    BioSpaceGroup {
        number: 14,
        ccp4: 0,
        hm: "P 21/b 1 1",
        ext: 0,
        qualifier: "a1",
        hall: "-P 2xab",
        basisop_idx: 2,
        ops: OPS_86,
    },
    BioSpaceGroup {
        number: 14,
        ccp4: 0,
        hm: "P 21/n 1 1",
        ext: 0,
        qualifier: "a2",
        hall: "-P 2xn",
        basisop_idx: 11,
        ops: OPS_87,
    },
    BioSpaceGroup {
        number: 14,
        ccp4: 0,
        hm: "P 21/c 1 1",
        ext: 0,
        qualifier: "a3",
        hall: "-P 2xac",
        basisop_idx: 7,
        ops: OPS_88,
    },
    BioSpaceGroup {
        number: 15,
        ccp4: 15,
        hm: "C 1 2/c 1",
        ext: 0,
        qualifier: "b1",
        hall: "-C 2yc",
        basisop_idx: 0,
        ops: OPS_89,
    },
    BioSpaceGroup {
        number: 15,
        ccp4: 0,
        hm: "A 1 2/n 1",
        ext: 0,
        qualifier: "b2",
        hall: "-A 2yab",
        basisop_idx: 12,
        ops: OPS_90,
    },
    BioSpaceGroup {
        number: 15,
        ccp4: 0,
        hm: "I 1 2/a 1",
        ext: 0,
        qualifier: "b3",
        hall: "-I 2ya",
        basisop_idx: 13,
        ops: OPS_91,
    },
    BioSpaceGroup {
        number: 15,
        ccp4: 0,
        hm: "A 1 2/a 1",
        ext: 0,
        qualifier: "-b1",
        hall: "-A 2ya",
        basisop_idx: 3,
        ops: OPS_92,
    },
    BioSpaceGroup {
        number: 15,
        ccp4: 0,
        hm: "C 1 2/n 1",
        ext: 0,
        qualifier: "-b2",
        hall: "-C 2yac",
        basisop_idx: 19,
        ops: OPS_93,
    },
    BioSpaceGroup {
        number: 15,
        ccp4: 0,
        hm: "I 1 2/c 1",
        ext: 0,
        qualifier: "-b3",
        hall: "-I 2yc",
        basisop_idx: 4,
        ops: OPS_94,
    },
    BioSpaceGroup {
        number: 15,
        ccp4: 0,
        hm: "A 1 1 2/a",
        ext: 0,
        qualifier: "c1",
        hall: "-A 2a",
        basisop_idx: 1,
        ops: OPS_95,
    },
    BioSpaceGroup {
        number: 15,
        ccp4: 0,
        hm: "B 1 1 2/n",
        ext: 0,
        qualifier: "c2",
        hall: "-B 2ab",
        basisop_idx: 15,
        ops: OPS_96,
    },
    BioSpaceGroup {
        number: 15,
        ccp4: 0,
        hm: "I 1 1 2/b",
        ext: 0,
        qualifier: "c3",
        hall: "-I 2b",
        basisop_idx: 16,
        ops: OPS_97,
    },
    BioSpaceGroup {
        number: 15,
        ccp4: 1015,
        hm: "B 1 1 2/b",
        ext: 0,
        qualifier: "-c1",
        hall: "-B 2b",
        basisop_idx: 5,
        ops: OPS_98,
    },
    BioSpaceGroup {
        number: 15,
        ccp4: 0,
        hm: "A 1 1 2/n",
        ext: 0,
        qualifier: "-c2",
        hall: "-A 2ab",
        basisop_idx: 10,
        ops: OPS_99,
    },
    BioSpaceGroup {
        number: 15,
        ccp4: 0,
        hm: "I 1 1 2/a",
        ext: 0,
        qualifier: "-c3",
        hall: "-I 2a",
        basisop_idx: 6,
        ops: OPS_100,
    },
    BioSpaceGroup {
        number: 15,
        ccp4: 0,
        hm: "B 2/b 1 1",
        ext: 0,
        qualifier: "a1",
        hall: "-B 2xb",
        basisop_idx: 2,
        ops: OPS_101,
    },
    BioSpaceGroup {
        number: 15,
        ccp4: 0,
        hm: "C 2/n 1 1",
        ext: 0,
        qualifier: "a2",
        hall: "-C 2xac",
        basisop_idx: 17,
        ops: OPS_102,
    },
    BioSpaceGroup {
        number: 15,
        ccp4: 0,
        hm: "I 2/c 1 1",
        ext: 0,
        qualifier: "a3",
        hall: "-I 2xc",
        basisop_idx: 18,
        ops: OPS_103,
    },
    BioSpaceGroup {
        number: 15,
        ccp4: 0,
        hm: "C 2/c 1 1",
        ext: 0,
        qualifier: "-a1",
        hall: "-C 2xc",
        basisop_idx: 7,
        ops: OPS_104,
    },
    BioSpaceGroup {
        number: 15,
        ccp4: 0,
        hm: "B 2/n 1 1",
        ext: 0,
        qualifier: "-a2",
        hall: "-B 2xab",
        basisop_idx: 11,
        ops: OPS_105,
    },
    BioSpaceGroup {
        number: 15,
        ccp4: 0,
        hm: "I 2/b 1 1",
        ext: 0,
        qualifier: "-a3",
        hall: "-I 2xb",
        basisop_idx: 8,
        ops: OPS_106,
    },
    BioSpaceGroup {
        number: 16,
        ccp4: 16,
        hm: "P 2 2 2",
        ext: 0,
        qualifier: "",
        hall: "P 2 2",
        basisop_idx: 0,
        ops: OPS_107,
    },
    BioSpaceGroup {
        number: 17,
        ccp4: 17,
        hm: "P 2 2 21",
        ext: 0,
        qualifier: "",
        hall: "P 2c 2",
        basisop_idx: 0,
        ops: OPS_108,
    },
    BioSpaceGroup {
        number: 17,
        ccp4: 1017,
        hm: "P 21 2 2",
        ext: 0,
        qualifier: "cab",
        hall: "P 2a 2a",
        basisop_idx: 1,
        ops: OPS_109,
    },
    BioSpaceGroup {
        number: 17,
        ccp4: 2017,
        hm: "P 2 21 2",
        ext: 0,
        qualifier: "bca",
        hall: "P 2 2b",
        basisop_idx: 2,
        ops: OPS_110,
    },
    BioSpaceGroup {
        number: 18,
        ccp4: 18,
        hm: "P 21 21 2",
        ext: 0,
        qualifier: "",
        hall: "P 2 2ab",
        basisop_idx: 0,
        ops: OPS_111,
    },
    BioSpaceGroup {
        number: 18,
        ccp4: 3018,
        hm: "P 2 21 21",
        ext: 0,
        qualifier: "cab",
        hall: "P 2bc 2",
        basisop_idx: 1,
        ops: OPS_112,
    },
    BioSpaceGroup {
        number: 18,
        ccp4: 2018,
        hm: "P 21 2 21",
        ext: 0,
        qualifier: "bca",
        hall: "P 2ac 2ac",
        basisop_idx: 2,
        ops: OPS_113,
    },
    BioSpaceGroup {
        number: 19,
        ccp4: 19,
        hm: "P 21 21 21",
        ext: 0,
        qualifier: "",
        hall: "P 2ac 2ab",
        basisop_idx: 0,
        ops: OPS_114,
    },
    BioSpaceGroup {
        number: 20,
        ccp4: 20,
        hm: "C 2 2 21",
        ext: 0,
        qualifier: "",
        hall: "C 2c 2",
        basisop_idx: 0,
        ops: OPS_115,
    },
    BioSpaceGroup {
        number: 20,
        ccp4: 0,
        hm: "A 21 2 2",
        ext: 0,
        qualifier: "cab",
        hall: "A 2a 2a",
        basisop_idx: 1,
        ops: OPS_116,
    },
    BioSpaceGroup {
        number: 20,
        ccp4: 0,
        hm: "B 2 21 2",
        ext: 0,
        qualifier: "bca",
        hall: "B 2 2b",
        basisop_idx: 2,
        ops: OPS_117,
    },
    BioSpaceGroup {
        number: 21,
        ccp4: 21,
        hm: "C 2 2 2",
        ext: 0,
        qualifier: "",
        hall: "C 2 2",
        basisop_idx: 0,
        ops: OPS_118,
    },
    BioSpaceGroup {
        number: 21,
        ccp4: 0,
        hm: "A 2 2 2",
        ext: 0,
        qualifier: "cab",
        hall: "A 2 2",
        basisop_idx: 1,
        ops: OPS_119,
    },
    BioSpaceGroup {
        number: 21,
        ccp4: 0,
        hm: "B 2 2 2",
        ext: 0,
        qualifier: "bca",
        hall: "B 2 2",
        basisop_idx: 2,
        ops: OPS_120,
    },
    BioSpaceGroup {
        number: 22,
        ccp4: 22,
        hm: "F 2 2 2",
        ext: 0,
        qualifier: "",
        hall: "F 2 2",
        basisop_idx: 0,
        ops: OPS_121,
    },
    BioSpaceGroup {
        number: 23,
        ccp4: 23,
        hm: "I 2 2 2",
        ext: 0,
        qualifier: "",
        hall: "I 2 2",
        basisop_idx: 0,
        ops: OPS_122,
    },
    BioSpaceGroup {
        number: 24,
        ccp4: 24,
        hm: "I 21 21 21",
        ext: 0,
        qualifier: "",
        hall: "I 2b 2c",
        basisop_idx: 0,
        ops: OPS_123,
    },
    BioSpaceGroup {
        number: 25,
        ccp4: 25,
        hm: "P m m 2",
        ext: 0,
        qualifier: "",
        hall: "P 2 -2",
        basisop_idx: 0,
        ops: OPS_124,
    },
    BioSpaceGroup {
        number: 25,
        ccp4: 0,
        hm: "P 2 m m",
        ext: 0,
        qualifier: "cab",
        hall: "P -2 2",
        basisop_idx: 1,
        ops: OPS_125,
    },
    BioSpaceGroup {
        number: 25,
        ccp4: 0,
        hm: "P m 2 m",
        ext: 0,
        qualifier: "bca",
        hall: "P -2 -2",
        basisop_idx: 2,
        ops: OPS_126,
    },
    BioSpaceGroup {
        number: 26,
        ccp4: 26,
        hm: "P m c 21",
        ext: 0,
        qualifier: "",
        hall: "P 2c -2",
        basisop_idx: 0,
        ops: OPS_127,
    },
    BioSpaceGroup {
        number: 26,
        ccp4: 0,
        hm: "P c m 21",
        ext: 0,
        qualifier: "ba-c",
        hall: "P 2c -2c",
        basisop_idx: 7,
        ops: OPS_128,
    },
    BioSpaceGroup {
        number: 26,
        ccp4: 0,
        hm: "P 21 m a",
        ext: 0,
        qualifier: "cab",
        hall: "P -2a 2a",
        basisop_idx: 1,
        ops: OPS_129,
    },
    BioSpaceGroup {
        number: 26,
        ccp4: 0,
        hm: "P 21 a m",
        ext: 0,
        qualifier: "-cba",
        hall: "P -2 2a",
        basisop_idx: 3,
        ops: OPS_130,
    },
    BioSpaceGroup {
        number: 26,
        ccp4: 0,
        hm: "P b 21 m",
        ext: 0,
        qualifier: "bca",
        hall: "P -2 -2b",
        basisop_idx: 2,
        ops: OPS_131,
    },
    BioSpaceGroup {
        number: 26,
        ccp4: 0,
        hm: "P m 21 b",
        ext: 0,
        qualifier: "a-cb",
        hall: "P -2b -2",
        basisop_idx: 5,
        ops: OPS_132,
    },
    BioSpaceGroup {
        number: 27,
        ccp4: 27,
        hm: "P c c 2",
        ext: 0,
        qualifier: "",
        hall: "P 2 -2c",
        basisop_idx: 0,
        ops: OPS_133,
    },
    BioSpaceGroup {
        number: 27,
        ccp4: 0,
        hm: "P 2 a a",
        ext: 0,
        qualifier: "cab",
        hall: "P -2a 2",
        basisop_idx: 1,
        ops: OPS_134,
    },
    BioSpaceGroup {
        number: 27,
        ccp4: 0,
        hm: "P b 2 b",
        ext: 0,
        qualifier: "bca",
        hall: "P -2b -2b",
        basisop_idx: 2,
        ops: OPS_135,
    },
    BioSpaceGroup {
        number: 28,
        ccp4: 28,
        hm: "P m a 2",
        ext: 0,
        qualifier: "",
        hall: "P 2 -2a",
        basisop_idx: 0,
        ops: OPS_136,
    },
    BioSpaceGroup {
        number: 28,
        ccp4: 0,
        hm: "P b m 2",
        ext: 0,
        qualifier: "ba-c",
        hall: "P 2 -2b",
        basisop_idx: 7,
        ops: OPS_137,
    },
    BioSpaceGroup {
        number: 28,
        ccp4: 0,
        hm: "P 2 m b",
        ext: 0,
        qualifier: "cab",
        hall: "P -2b 2",
        basisop_idx: 1,
        ops: OPS_138,
    },
    BioSpaceGroup {
        number: 28,
        ccp4: 0,
        hm: "P 2 c m",
        ext: 0,
        qualifier: "-cba",
        hall: "P -2c 2",
        basisop_idx: 3,
        ops: OPS_139,
    },
    BioSpaceGroup {
        number: 28,
        ccp4: 0,
        hm: "P c 2 m",
        ext: 0,
        qualifier: "bca",
        hall: "P -2c -2c",
        basisop_idx: 2,
        ops: OPS_140,
    },
    BioSpaceGroup {
        number: 28,
        ccp4: 0,
        hm: "P m 2 a",
        ext: 0,
        qualifier: "a-cb",
        hall: "P -2a -2a",
        basisop_idx: 5,
        ops: OPS_141,
    },
    BioSpaceGroup {
        number: 29,
        ccp4: 29,
        hm: "P c a 21",
        ext: 0,
        qualifier: "",
        hall: "P 2c -2ac",
        basisop_idx: 0,
        ops: OPS_142,
    },
    BioSpaceGroup {
        number: 29,
        ccp4: 0,
        hm: "P b c 21",
        ext: 0,
        qualifier: "ba-c",
        hall: "P 2c -2b",
        basisop_idx: 7,
        ops: OPS_143,
    },
    BioSpaceGroup {
        number: 29,
        ccp4: 0,
        hm: "P 21 a b",
        ext: 0,
        qualifier: "cab",
        hall: "P -2b 2a",
        basisop_idx: 1,
        ops: OPS_144,
    },
    BioSpaceGroup {
        number: 29,
        ccp4: 0,
        hm: "P 21 c a",
        ext: 0,
        qualifier: "-cba",
        hall: "P -2ac 2a",
        basisop_idx: 3,
        ops: OPS_145,
    },
    BioSpaceGroup {
        number: 29,
        ccp4: 0,
        hm: "P c 21 b",
        ext: 0,
        qualifier: "bca",
        hall: "P -2bc -2c",
        basisop_idx: 2,
        ops: OPS_146,
    },
    BioSpaceGroup {
        number: 29,
        ccp4: 0,
        hm: "P b 21 a",
        ext: 0,
        qualifier: "a-cb",
        hall: "P -2a -2ab",
        basisop_idx: 5,
        ops: OPS_147,
    },
    BioSpaceGroup {
        number: 30,
        ccp4: 30,
        hm: "P n c 2",
        ext: 0,
        qualifier: "",
        hall: "P 2 -2bc",
        basisop_idx: 0,
        ops: OPS_148,
    },
    BioSpaceGroup {
        number: 30,
        ccp4: 0,
        hm: "P c n 2",
        ext: 0,
        qualifier: "ba-c",
        hall: "P 2 -2ac",
        basisop_idx: 7,
        ops: OPS_149,
    },
    BioSpaceGroup {
        number: 30,
        ccp4: 0,
        hm: "P 2 n a",
        ext: 0,
        qualifier: "cab",
        hall: "P -2ac 2",
        basisop_idx: 1,
        ops: OPS_150,
    },
    BioSpaceGroup {
        number: 30,
        ccp4: 0,
        hm: "P 2 a n",
        ext: 0,
        qualifier: "-cba",
        hall: "P -2ab 2",
        basisop_idx: 3,
        ops: OPS_151,
    },
    BioSpaceGroup {
        number: 30,
        ccp4: 0,
        hm: "P b 2 n",
        ext: 0,
        qualifier: "bca",
        hall: "P -2ab -2ab",
        basisop_idx: 2,
        ops: OPS_152,
    },
    BioSpaceGroup {
        number: 30,
        ccp4: 0,
        hm: "P n 2 b",
        ext: 0,
        qualifier: "a-cb",
        hall: "P -2bc -2bc",
        basisop_idx: 5,
        ops: OPS_153,
    },
    BioSpaceGroup {
        number: 31,
        ccp4: 31,
        hm: "P m n 21",
        ext: 0,
        qualifier: "",
        hall: "P 2ac -2",
        basisop_idx: 0,
        ops: OPS_154,
    },
    BioSpaceGroup {
        number: 31,
        ccp4: 0,
        hm: "P n m 21",
        ext: 0,
        qualifier: "ba-c",
        hall: "P 2bc -2bc",
        basisop_idx: 7,
        ops: OPS_155,
    },
    BioSpaceGroup {
        number: 31,
        ccp4: 0,
        hm: "P 21 m n",
        ext: 0,
        qualifier: "cab",
        hall: "P -2ab 2ab",
        basisop_idx: 1,
        ops: OPS_156,
    },
    BioSpaceGroup {
        number: 31,
        ccp4: 0,
        hm: "P 21 n m",
        ext: 0,
        qualifier: "-cba",
        hall: "P -2 2ac",
        basisop_idx: 3,
        ops: OPS_157,
    },
    BioSpaceGroup {
        number: 31,
        ccp4: 0,
        hm: "P n 21 m",
        ext: 0,
        qualifier: "bca",
        hall: "P -2 -2bc",
        basisop_idx: 2,
        ops: OPS_158,
    },
    BioSpaceGroup {
        number: 31,
        ccp4: 0,
        hm: "P m 21 n",
        ext: 0,
        qualifier: "a-cb",
        hall: "P -2ab -2",
        basisop_idx: 5,
        ops: OPS_159,
    },
    BioSpaceGroup {
        number: 32,
        ccp4: 32,
        hm: "P b a 2",
        ext: 0,
        qualifier: "",
        hall: "P 2 -2ab",
        basisop_idx: 0,
        ops: OPS_160,
    },
    BioSpaceGroup {
        number: 32,
        ccp4: 0,
        hm: "P 2 c b",
        ext: 0,
        qualifier: "cab",
        hall: "P -2bc 2",
        basisop_idx: 1,
        ops: OPS_161,
    },
    BioSpaceGroup {
        number: 32,
        ccp4: 0,
        hm: "P c 2 a",
        ext: 0,
        qualifier: "bca",
        hall: "P -2ac -2ac",
        basisop_idx: 2,
        ops: OPS_162,
    },
    BioSpaceGroup {
        number: 33,
        ccp4: 33,
        hm: "P n a 21",
        ext: 0,
        qualifier: "",
        hall: "P 2c -2n",
        basisop_idx: 0,
        ops: OPS_163,
    },
    BioSpaceGroup {
        number: 33,
        ccp4: 0,
        hm: "P b n 21",
        ext: 0,
        qualifier: "ba-c",
        hall: "P 2c -2ab",
        basisop_idx: 7,
        ops: OPS_164,
    },
    BioSpaceGroup {
        number: 33,
        ccp4: 0,
        hm: "P 21 n b",
        ext: 0,
        qualifier: "cab",
        hall: "P -2bc 2a",
        basisop_idx: 1,
        ops: OPS_165,
    },
    BioSpaceGroup {
        number: 33,
        ccp4: 0,
        hm: "P 21 c n",
        ext: 0,
        qualifier: "-cba",
        hall: "P -2n 2a",
        basisop_idx: 3,
        ops: OPS_166,
    },
    BioSpaceGroup {
        number: 33,
        ccp4: 0,
        hm: "P c 21 n",
        ext: 0,
        qualifier: "bca",
        hall: "P -2n -2ac",
        basisop_idx: 2,
        ops: OPS_167,
    },
    BioSpaceGroup {
        number: 33,
        ccp4: 0,
        hm: "P n 21 a",
        ext: 0,
        qualifier: "a-cb",
        hall: "P -2ac -2n",
        basisop_idx: 5,
        ops: OPS_168,
    },
    BioSpaceGroup {
        number: 34,
        ccp4: 34,
        hm: "P n n 2",
        ext: 0,
        qualifier: "",
        hall: "P 2 -2n",
        basisop_idx: 0,
        ops: OPS_169,
    },
    BioSpaceGroup {
        number: 34,
        ccp4: 0,
        hm: "P 2 n n",
        ext: 0,
        qualifier: "cab",
        hall: "P -2n 2",
        basisop_idx: 1,
        ops: OPS_170,
    },
    BioSpaceGroup {
        number: 34,
        ccp4: 0,
        hm: "P n 2 n",
        ext: 0,
        qualifier: "bca",
        hall: "P -2n -2n",
        basisop_idx: 2,
        ops: OPS_171,
    },
    BioSpaceGroup {
        number: 35,
        ccp4: 35,
        hm: "C m m 2",
        ext: 0,
        qualifier: "",
        hall: "C 2 -2",
        basisop_idx: 0,
        ops: OPS_172,
    },
    BioSpaceGroup {
        number: 35,
        ccp4: 0,
        hm: "A 2 m m",
        ext: 0,
        qualifier: "cab",
        hall: "A -2 2",
        basisop_idx: 1,
        ops: OPS_173,
    },
    BioSpaceGroup {
        number: 35,
        ccp4: 0,
        hm: "B m 2 m",
        ext: 0,
        qualifier: "bca",
        hall: "B -2 -2",
        basisop_idx: 2,
        ops: OPS_174,
    },
    BioSpaceGroup {
        number: 36,
        ccp4: 36,
        hm: "C m c 21",
        ext: 0,
        qualifier: "",
        hall: "C 2c -2",
        basisop_idx: 0,
        ops: OPS_175,
    },
    BioSpaceGroup {
        number: 36,
        ccp4: 0,
        hm: "C c m 21",
        ext: 0,
        qualifier: "ba-c",
        hall: "C 2c -2c",
        basisop_idx: 7,
        ops: OPS_176,
    },
    BioSpaceGroup {
        number: 36,
        ccp4: 0,
        hm: "A 21 m a",
        ext: 0,
        qualifier: "cab",
        hall: "A -2a 2a",
        basisop_idx: 1,
        ops: OPS_177,
    },
    BioSpaceGroup {
        number: 36,
        ccp4: 0,
        hm: "A 21 a m",
        ext: 0,
        qualifier: "-cba",
        hall: "A -2 2a",
        basisop_idx: 3,
        ops: OPS_178,
    },
    BioSpaceGroup {
        number: 36,
        ccp4: 0,
        hm: "B b 21 m",
        ext: 0,
        qualifier: "bca",
        hall: "B -2 -2b",
        basisop_idx: 2,
        ops: OPS_179,
    },
    BioSpaceGroup {
        number: 36,
        ccp4: 0,
        hm: "B m 21 b",
        ext: 0,
        qualifier: "a-cb",
        hall: "B -2b -2",
        basisop_idx: 5,
        ops: OPS_180,
    },
    BioSpaceGroup {
        number: 37,
        ccp4: 37,
        hm: "C c c 2",
        ext: 0,
        qualifier: "",
        hall: "C 2 -2c",
        basisop_idx: 0,
        ops: OPS_181,
    },
    BioSpaceGroup {
        number: 37,
        ccp4: 0,
        hm: "A 2 a a",
        ext: 0,
        qualifier: "cab",
        hall: "A -2a 2",
        basisop_idx: 1,
        ops: OPS_182,
    },
    BioSpaceGroup {
        number: 37,
        ccp4: 0,
        hm: "B b 2 b",
        ext: 0,
        qualifier: "bca",
        hall: "B -2b -2b",
        basisop_idx: 2,
        ops: OPS_183,
    },
    BioSpaceGroup {
        number: 38,
        ccp4: 38,
        hm: "A m m 2",
        ext: 0,
        qualifier: "",
        hall: "A 2 -2",
        basisop_idx: 0,
        ops: OPS_184,
    },
    BioSpaceGroup {
        number: 38,
        ccp4: 0,
        hm: "B m m 2",
        ext: 0,
        qualifier: "ba-c",
        hall: "B 2 -2",
        basisop_idx: 7,
        ops: OPS_185,
    },
    BioSpaceGroup {
        number: 38,
        ccp4: 0,
        hm: "B 2 m m",
        ext: 0,
        qualifier: "cab",
        hall: "B -2 2",
        basisop_idx: 1,
        ops: OPS_186,
    },
    BioSpaceGroup {
        number: 38,
        ccp4: 0,
        hm: "C 2 m m",
        ext: 0,
        qualifier: "-cba",
        hall: "C -2 2",
        basisop_idx: 3,
        ops: OPS_187,
    },
    BioSpaceGroup {
        number: 38,
        ccp4: 0,
        hm: "C m 2 m",
        ext: 0,
        qualifier: "bca",
        hall: "C -2 -2",
        basisop_idx: 2,
        ops: OPS_188,
    },
    BioSpaceGroup {
        number: 38,
        ccp4: 0,
        hm: "A m 2 m",
        ext: 0,
        qualifier: "a-cb",
        hall: "A -2 -2",
        basisop_idx: 5,
        ops: OPS_189,
    },
    BioSpaceGroup {
        number: 39,
        ccp4: 39,
        hm: "A b m 2",
        ext: 0,
        qualifier: "",
        hall: "A 2 -2b",
        basisop_idx: 0,
        ops: OPS_190,
    },
    BioSpaceGroup {
        number: 39,
        ccp4: 0,
        hm: "B m a 2",
        ext: 0,
        qualifier: "ba-c",
        hall: "B 2 -2a",
        basisop_idx: 7,
        ops: OPS_191,
    },
    BioSpaceGroup {
        number: 39,
        ccp4: 0,
        hm: "B 2 c m",
        ext: 0,
        qualifier: "cab",
        hall: "B -2a 2",
        basisop_idx: 1,
        ops: OPS_192,
    },
    BioSpaceGroup {
        number: 39,
        ccp4: 0,
        hm: "C 2 m b",
        ext: 0,
        qualifier: "-cba",
        hall: "C -2a 2",
        basisop_idx: 3,
        ops: OPS_193,
    },
    BioSpaceGroup {
        number: 39,
        ccp4: 0,
        hm: "C m 2 a",
        ext: 0,
        qualifier: "bca",
        hall: "C -2a -2a",
        basisop_idx: 2,
        ops: OPS_194,
    },
    BioSpaceGroup {
        number: 39,
        ccp4: 0,
        hm: "A c 2 m",
        ext: 0,
        qualifier: "a-cb",
        hall: "A -2b -2b",
        basisop_idx: 5,
        ops: OPS_195,
    },
    BioSpaceGroup {
        number: 40,
        ccp4: 40,
        hm: "A m a 2",
        ext: 0,
        qualifier: "",
        hall: "A 2 -2a",
        basisop_idx: 0,
        ops: OPS_196,
    },
    BioSpaceGroup {
        number: 40,
        ccp4: 0,
        hm: "B b m 2",
        ext: 0,
        qualifier: "ba-c",
        hall: "B 2 -2b",
        basisop_idx: 7,
        ops: OPS_197,
    },
    BioSpaceGroup {
        number: 40,
        ccp4: 0,
        hm: "B 2 m b",
        ext: 0,
        qualifier: "cab",
        hall: "B -2b 2",
        basisop_idx: 1,
        ops: OPS_198,
    },
    BioSpaceGroup {
        number: 40,
        ccp4: 0,
        hm: "C 2 c m",
        ext: 0,
        qualifier: "-cba",
        hall: "C -2c 2",
        basisop_idx: 3,
        ops: OPS_199,
    },
    BioSpaceGroup {
        number: 40,
        ccp4: 0,
        hm: "C c 2 m",
        ext: 0,
        qualifier: "bca",
        hall: "C -2c -2c",
        basisop_idx: 2,
        ops: OPS_200,
    },
    BioSpaceGroup {
        number: 40,
        ccp4: 0,
        hm: "A m 2 a",
        ext: 0,
        qualifier: "a-cb",
        hall: "A -2a -2a",
        basisop_idx: 5,
        ops: OPS_201,
    },
    BioSpaceGroup {
        number: 41,
        ccp4: 41,
        hm: "A b a 2",
        ext: 0,
        qualifier: "",
        hall: "A 2 -2ab",
        basisop_idx: 0,
        ops: OPS_202,
    },
    BioSpaceGroup {
        number: 41,
        ccp4: 0,
        hm: "B b a 2",
        ext: 0,
        qualifier: "ba-c",
        hall: "B 2 -2ab",
        basisop_idx: 7,
        ops: OPS_203,
    },
    BioSpaceGroup {
        number: 41,
        ccp4: 0,
        hm: "B 2 c b",
        ext: 0,
        qualifier: "cab",
        hall: "B -2ab 2",
        basisop_idx: 1,
        ops: OPS_204,
    },
    BioSpaceGroup {
        number: 41,
        ccp4: 0,
        hm: "C 2 c b",
        ext: 0,
        qualifier: "-cba",
        hall: "C -2ac 2",
        basisop_idx: 3,
        ops: OPS_205,
    },
    BioSpaceGroup {
        number: 41,
        ccp4: 0,
        hm: "C c 2 a",
        ext: 0,
        qualifier: "bca",
        hall: "C -2ac -2ac",
        basisop_idx: 2,
        ops: OPS_206,
    },
    BioSpaceGroup {
        number: 41,
        ccp4: 0,
        hm: "A c 2 a",
        ext: 0,
        qualifier: "a-cb",
        hall: "A -2ab -2ab",
        basisop_idx: 5,
        ops: OPS_207,
    },
    BioSpaceGroup {
        number: 42,
        ccp4: 42,
        hm: "F m m 2",
        ext: 0,
        qualifier: "",
        hall: "F 2 -2",
        basisop_idx: 0,
        ops: OPS_208,
    },
    BioSpaceGroup {
        number: 42,
        ccp4: 0,
        hm: "F 2 m m",
        ext: 0,
        qualifier: "cab",
        hall: "F -2 2",
        basisop_idx: 1,
        ops: OPS_209,
    },
    BioSpaceGroup {
        number: 42,
        ccp4: 0,
        hm: "F m 2 m",
        ext: 0,
        qualifier: "bca",
        hall: "F -2 -2",
        basisop_idx: 2,
        ops: OPS_210,
    },
    BioSpaceGroup {
        number: 43,
        ccp4: 43,
        hm: "F d d 2",
        ext: 0,
        qualifier: "",
        hall: "F 2 -2d",
        basisop_idx: 0,
        ops: OPS_211,
    },
    BioSpaceGroup {
        number: 43,
        ccp4: 0,
        hm: "F 2 d d",
        ext: 0,
        qualifier: "cab",
        hall: "F -2d 2",
        basisop_idx: 1,
        ops: OPS_212,
    },
    BioSpaceGroup {
        number: 43,
        ccp4: 0,
        hm: "F d 2 d",
        ext: 0,
        qualifier: "bca",
        hall: "F -2d -2d",
        basisop_idx: 2,
        ops: OPS_213,
    },
    BioSpaceGroup {
        number: 44,
        ccp4: 44,
        hm: "I m m 2",
        ext: 0,
        qualifier: "",
        hall: "I 2 -2",
        basisop_idx: 0,
        ops: OPS_214,
    },
    BioSpaceGroup {
        number: 44,
        ccp4: 0,
        hm: "I 2 m m",
        ext: 0,
        qualifier: "cab",
        hall: "I -2 2",
        basisop_idx: 1,
        ops: OPS_215,
    },
    BioSpaceGroup {
        number: 44,
        ccp4: 0,
        hm: "I m 2 m",
        ext: 0,
        qualifier: "bca",
        hall: "I -2 -2",
        basisop_idx: 2,
        ops: OPS_216,
    },
    BioSpaceGroup {
        number: 45,
        ccp4: 45,
        hm: "I b a 2",
        ext: 0,
        qualifier: "",
        hall: "I 2 -2c",
        basisop_idx: 0,
        ops: OPS_217,
    },
    BioSpaceGroup {
        number: 45,
        ccp4: 0,
        hm: "I 2 c b",
        ext: 0,
        qualifier: "cab",
        hall: "I -2a 2",
        basisop_idx: 1,
        ops: OPS_218,
    },
    BioSpaceGroup {
        number: 45,
        ccp4: 0,
        hm: "I c 2 a",
        ext: 0,
        qualifier: "bca",
        hall: "I -2b -2b",
        basisop_idx: 2,
        ops: OPS_219,
    },
    BioSpaceGroup {
        number: 46,
        ccp4: 46,
        hm: "I m a 2",
        ext: 0,
        qualifier: "",
        hall: "I 2 -2a",
        basisop_idx: 0,
        ops: OPS_220,
    },
    BioSpaceGroup {
        number: 46,
        ccp4: 0,
        hm: "I b m 2",
        ext: 0,
        qualifier: "ba-c",
        hall: "I 2 -2b",
        basisop_idx: 7,
        ops: OPS_221,
    },
    BioSpaceGroup {
        number: 46,
        ccp4: 0,
        hm: "I 2 m b",
        ext: 0,
        qualifier: "cab",
        hall: "I -2b 2",
        basisop_idx: 1,
        ops: OPS_222,
    },
    BioSpaceGroup {
        number: 46,
        ccp4: 0,
        hm: "I 2 c m",
        ext: 0,
        qualifier: "-cba",
        hall: "I -2c 2",
        basisop_idx: 3,
        ops: OPS_223,
    },
    BioSpaceGroup {
        number: 46,
        ccp4: 0,
        hm: "I c 2 m",
        ext: 0,
        qualifier: "bca",
        hall: "I -2c -2c",
        basisop_idx: 2,
        ops: OPS_224,
    },
    BioSpaceGroup {
        number: 46,
        ccp4: 0,
        hm: "I m 2 a",
        ext: 0,
        qualifier: "a-cb",
        hall: "I -2a -2a",
        basisop_idx: 5,
        ops: OPS_225,
    },
    BioSpaceGroup {
        number: 47,
        ccp4: 47,
        hm: "P m m m",
        ext: 0,
        qualifier: "",
        hall: "-P 2 2",
        basisop_idx: 0,
        ops: OPS_226,
    },
    BioSpaceGroup {
        number: 48,
        ccp4: 48,
        hm: "P n n n",
        ext: 49,
        qualifier: "",
        hall: "P 2 2 -1n",
        basisop_idx: 20,
        ops: OPS_227,
    },
    BioSpaceGroup {
        number: 48,
        ccp4: 0,
        hm: "P n n n",
        ext: 50,
        qualifier: "",
        hall: "-P 2ab 2bc",
        basisop_idx: 0,
        ops: OPS_228,
    },
    BioSpaceGroup {
        number: 49,
        ccp4: 49,
        hm: "P c c m",
        ext: 0,
        qualifier: "",
        hall: "-P 2 2c",
        basisop_idx: 0,
        ops: OPS_229,
    },
    BioSpaceGroup {
        number: 49,
        ccp4: 0,
        hm: "P m a a",
        ext: 0,
        qualifier: "cab",
        hall: "-P 2a 2",
        basisop_idx: 1,
        ops: OPS_230,
    },
    BioSpaceGroup {
        number: 49,
        ccp4: 0,
        hm: "P b m b",
        ext: 0,
        qualifier: "bca",
        hall: "-P 2b 2b",
        basisop_idx: 2,
        ops: OPS_231,
    },
    BioSpaceGroup {
        number: 50,
        ccp4: 50,
        hm: "P b a n",
        ext: 49,
        qualifier: "",
        hall: "P 2 2 -1ab",
        basisop_idx: 21,
        ops: OPS_232,
    },
    BioSpaceGroup {
        number: 50,
        ccp4: 0,
        hm: "P b a n",
        ext: 50,
        qualifier: "",
        hall: "-P 2ab 2b",
        basisop_idx: 0,
        ops: OPS_233,
    },
    BioSpaceGroup {
        number: 50,
        ccp4: 0,
        hm: "P n c b",
        ext: 49,
        qualifier: "cab",
        hall: "P 2 2 -1bc",
        basisop_idx: 22,
        ops: OPS_234,
    },
    BioSpaceGroup {
        number: 50,
        ccp4: 0,
        hm: "P n c b",
        ext: 50,
        qualifier: "cab",
        hall: "-P 2b 2bc",
        basisop_idx: 1,
        ops: OPS_235,
    },
    BioSpaceGroup {
        number: 50,
        ccp4: 0,
        hm: "P c n a",
        ext: 49,
        qualifier: "bca",
        hall: "P 2 2 -1ac",
        basisop_idx: 23,
        ops: OPS_236,
    },
    BioSpaceGroup {
        number: 50,
        ccp4: 0,
        hm: "P c n a",
        ext: 50,
        qualifier: "bca",
        hall: "-P 2a 2c",
        basisop_idx: 2,
        ops: OPS_237,
    },
    BioSpaceGroup {
        number: 51,
        ccp4: 51,
        hm: "P m m a",
        ext: 0,
        qualifier: "",
        hall: "-P 2a 2a",
        basisop_idx: 0,
        ops: OPS_238,
    },
    BioSpaceGroup {
        number: 51,
        ccp4: 0,
        hm: "P m m b",
        ext: 0,
        qualifier: "ba-c",
        hall: "-P 2b 2",
        basisop_idx: 7,
        ops: OPS_239,
    },
    BioSpaceGroup {
        number: 51,
        ccp4: 0,
        hm: "P b m m",
        ext: 0,
        qualifier: "cab",
        hall: "-P 2 2b",
        basisop_idx: 1,
        ops: OPS_240,
    },
    BioSpaceGroup {
        number: 51,
        ccp4: 0,
        hm: "P c m m",
        ext: 0,
        qualifier: "-cba",
        hall: "-P 2c 2c",
        basisop_idx: 3,
        ops: OPS_241,
    },
    BioSpaceGroup {
        number: 51,
        ccp4: 0,
        hm: "P m c m",
        ext: 0,
        qualifier: "bca",
        hall: "-P 2c 2",
        basisop_idx: 2,
        ops: OPS_242,
    },
    BioSpaceGroup {
        number: 51,
        ccp4: 0,
        hm: "P m a m",
        ext: 0,
        qualifier: "a-cb",
        hall: "-P 2 2a",
        basisop_idx: 5,
        ops: OPS_243,
    },
    BioSpaceGroup {
        number: 52,
        ccp4: 52,
        hm: "P n n a",
        ext: 0,
        qualifier: "",
        hall: "-P 2a 2bc",
        basisop_idx: 0,
        ops: OPS_244,
    },
    BioSpaceGroup {
        number: 52,
        ccp4: 0,
        hm: "P n n b",
        ext: 0,
        qualifier: "ba-c",
        hall: "-P 2b 2n",
        basisop_idx: 7,
        ops: OPS_245,
    },
    BioSpaceGroup {
        number: 52,
        ccp4: 0,
        hm: "P b n n",
        ext: 0,
        qualifier: "cab",
        hall: "-P 2n 2b",
        basisop_idx: 1,
        ops: OPS_246,
    },
    BioSpaceGroup {
        number: 52,
        ccp4: 0,
        hm: "P c n n",
        ext: 0,
        qualifier: "-cba",
        hall: "-P 2ab 2c",
        basisop_idx: 3,
        ops: OPS_247,
    },
    BioSpaceGroup {
        number: 52,
        ccp4: 0,
        hm: "P n c n",
        ext: 0,
        qualifier: "bca",
        hall: "-P 2ab 2n",
        basisop_idx: 2,
        ops: OPS_248,
    },
    BioSpaceGroup {
        number: 52,
        ccp4: 0,
        hm: "P n a n",
        ext: 0,
        qualifier: "a-cb",
        hall: "-P 2n 2bc",
        basisop_idx: 5,
        ops: OPS_249,
    },
    BioSpaceGroup {
        number: 53,
        ccp4: 53,
        hm: "P m n a",
        ext: 0,
        qualifier: "",
        hall: "-P 2ac 2",
        basisop_idx: 0,
        ops: OPS_250,
    },
    BioSpaceGroup {
        number: 53,
        ccp4: 0,
        hm: "P n m b",
        ext: 0,
        qualifier: "ba-c",
        hall: "-P 2bc 2bc",
        basisop_idx: 7,
        ops: OPS_251,
    },
    BioSpaceGroup {
        number: 53,
        ccp4: 0,
        hm: "P b m n",
        ext: 0,
        qualifier: "cab",
        hall: "-P 2ab 2ab",
        basisop_idx: 1,
        ops: OPS_252,
    },
    BioSpaceGroup {
        number: 53,
        ccp4: 0,
        hm: "P c n m",
        ext: 0,
        qualifier: "-cba",
        hall: "-P 2 2ac",
        basisop_idx: 3,
        ops: OPS_253,
    },
    BioSpaceGroup {
        number: 53,
        ccp4: 0,
        hm: "P n c m",
        ext: 0,
        qualifier: "bca",
        hall: "-P 2 2bc",
        basisop_idx: 2,
        ops: OPS_254,
    },
    BioSpaceGroup {
        number: 53,
        ccp4: 0,
        hm: "P m a n",
        ext: 0,
        qualifier: "a-cb",
        hall: "-P 2ab 2",
        basisop_idx: 5,
        ops: OPS_255,
    },
    BioSpaceGroup {
        number: 54,
        ccp4: 54,
        hm: "P c c a",
        ext: 0,
        qualifier: "",
        hall: "-P 2a 2ac",
        basisop_idx: 0,
        ops: OPS_256,
    },
    BioSpaceGroup {
        number: 54,
        ccp4: 0,
        hm: "P c c b",
        ext: 0,
        qualifier: "ba-c",
        hall: "-P 2b 2c",
        basisop_idx: 7,
        ops: OPS_257,
    },
    BioSpaceGroup {
        number: 54,
        ccp4: 0,
        hm: "P b a a",
        ext: 0,
        qualifier: "cab",
        hall: "-P 2a 2b",
        basisop_idx: 1,
        ops: OPS_258,
    },
    BioSpaceGroup {
        number: 54,
        ccp4: 0,
        hm: "P c a a",
        ext: 0,
        qualifier: "-cba",
        hall: "-P 2ac 2c",
        basisop_idx: 3,
        ops: OPS_259,
    },
    BioSpaceGroup {
        number: 54,
        ccp4: 0,
        hm: "P b c b",
        ext: 0,
        qualifier: "bca",
        hall: "-P 2bc 2b",
        basisop_idx: 2,
        ops: OPS_260,
    },
    BioSpaceGroup {
        number: 54,
        ccp4: 0,
        hm: "P b a b",
        ext: 0,
        qualifier: "a-cb",
        hall: "-P 2b 2ab",
        basisop_idx: 5,
        ops: OPS_261,
    },
    BioSpaceGroup {
        number: 55,
        ccp4: 55,
        hm: "P b a m",
        ext: 0,
        qualifier: "",
        hall: "-P 2 2ab",
        basisop_idx: 0,
        ops: OPS_262,
    },
    BioSpaceGroup {
        number: 55,
        ccp4: 0,
        hm: "P m c b",
        ext: 0,
        qualifier: "cab",
        hall: "-P 2bc 2",
        basisop_idx: 1,
        ops: OPS_263,
    },
    BioSpaceGroup {
        number: 55,
        ccp4: 0,
        hm: "P c m a",
        ext: 0,
        qualifier: "bca",
        hall: "-P 2ac 2ac",
        basisop_idx: 2,
        ops: OPS_264,
    },
    BioSpaceGroup {
        number: 56,
        ccp4: 56,
        hm: "P c c n",
        ext: 0,
        qualifier: "",
        hall: "-P 2ab 2ac",
        basisop_idx: 0,
        ops: OPS_265,
    },
    BioSpaceGroup {
        number: 56,
        ccp4: 0,
        hm: "P n a a",
        ext: 0,
        qualifier: "cab",
        hall: "-P 2ac 2bc",
        basisop_idx: 1,
        ops: OPS_266,
    },
    BioSpaceGroup {
        number: 56,
        ccp4: 0,
        hm: "P b n b",
        ext: 0,
        qualifier: "bca",
        hall: "-P 2bc 2ab",
        basisop_idx: 2,
        ops: OPS_267,
    },
    BioSpaceGroup {
        number: 57,
        ccp4: 57,
        hm: "P b c m",
        ext: 0,
        qualifier: "",
        hall: "-P 2c 2b",
        basisop_idx: 0,
        ops: OPS_268,
    },
    BioSpaceGroup {
        number: 57,
        ccp4: 0,
        hm: "P c a m",
        ext: 0,
        qualifier: "ba-c",
        hall: "-P 2c 2ac",
        basisop_idx: 7,
        ops: OPS_269,
    },
    BioSpaceGroup {
        number: 57,
        ccp4: 0,
        hm: "P m c a",
        ext: 0,
        qualifier: "cab",
        hall: "-P 2ac 2a",
        basisop_idx: 1,
        ops: OPS_270,
    },
    BioSpaceGroup {
        number: 57,
        ccp4: 0,
        hm: "P m a b",
        ext: 0,
        qualifier: "-cba",
        hall: "-P 2b 2a",
        basisop_idx: 3,
        ops: OPS_271,
    },
    BioSpaceGroup {
        number: 57,
        ccp4: 0,
        hm: "P b m a",
        ext: 0,
        qualifier: "bca",
        hall: "-P 2a 2ab",
        basisop_idx: 2,
        ops: OPS_272,
    },
    BioSpaceGroup {
        number: 57,
        ccp4: 0,
        hm: "P c m b",
        ext: 0,
        qualifier: "a-cb",
        hall: "-P 2bc 2c",
        basisop_idx: 5,
        ops: OPS_273,
    },
    BioSpaceGroup {
        number: 58,
        ccp4: 58,
        hm: "P n n m",
        ext: 0,
        qualifier: "",
        hall: "-P 2 2n",
        basisop_idx: 0,
        ops: OPS_274,
    },
    BioSpaceGroup {
        number: 58,
        ccp4: 0,
        hm: "P m n n",
        ext: 0,
        qualifier: "cab",
        hall: "-P 2n 2",
        basisop_idx: 1,
        ops: OPS_275,
    },
    BioSpaceGroup {
        number: 58,
        ccp4: 0,
        hm: "P n m n",
        ext: 0,
        qualifier: "bca",
        hall: "-P 2n 2n",
        basisop_idx: 2,
        ops: OPS_276,
    },
    BioSpaceGroup {
        number: 59,
        ccp4: 59,
        hm: "P m m n",
        ext: 49,
        qualifier: "",
        hall: "P 2 2ab -1ab",
        basisop_idx: 21,
        ops: OPS_277,
    },
    BioSpaceGroup {
        number: 59,
        ccp4: 1059,
        hm: "P m m n",
        ext: 50,
        qualifier: "",
        hall: "-P 2ab 2a",
        basisop_idx: 0,
        ops: OPS_278,
    },
    BioSpaceGroup {
        number: 59,
        ccp4: 0,
        hm: "P n m m",
        ext: 49,
        qualifier: "cab",
        hall: "P 2bc 2 -1bc",
        basisop_idx: 22,
        ops: OPS_279,
    },
    BioSpaceGroup {
        number: 59,
        ccp4: 0,
        hm: "P n m m",
        ext: 50,
        qualifier: "cab",
        hall: "-P 2c 2bc",
        basisop_idx: 1,
        ops: OPS_280,
    },
    BioSpaceGroup {
        number: 59,
        ccp4: 0,
        hm: "P m n m",
        ext: 49,
        qualifier: "bca",
        hall: "P 2ac 2ac -1ac",
        basisop_idx: 23,
        ops: OPS_281,
    },
    BioSpaceGroup {
        number: 59,
        ccp4: 0,
        hm: "P m n m",
        ext: 50,
        qualifier: "bca",
        hall: "-P 2c 2a",
        basisop_idx: 2,
        ops: OPS_282,
    },
    BioSpaceGroup {
        number: 60,
        ccp4: 60,
        hm: "P b c n",
        ext: 0,
        qualifier: "",
        hall: "-P 2n 2ab",
        basisop_idx: 0,
        ops: OPS_283,
    },
    BioSpaceGroup {
        number: 60,
        ccp4: 0,
        hm: "P c a n",
        ext: 0,
        qualifier: "ba-c",
        hall: "-P 2n 2c",
        basisop_idx: 7,
        ops: OPS_284,
    },
    BioSpaceGroup {
        number: 60,
        ccp4: 0,
        hm: "P n c a",
        ext: 0,
        qualifier: "cab",
        hall: "-P 2a 2n",
        basisop_idx: 1,
        ops: OPS_285,
    },
    BioSpaceGroup {
        number: 60,
        ccp4: 0,
        hm: "P n a b",
        ext: 0,
        qualifier: "-cba",
        hall: "-P 2bc 2n",
        basisop_idx: 3,
        ops: OPS_286,
    },
    BioSpaceGroup {
        number: 60,
        ccp4: 0,
        hm: "P b n a",
        ext: 0,
        qualifier: "bca",
        hall: "-P 2ac 2b",
        basisop_idx: 2,
        ops: OPS_287,
    },
    BioSpaceGroup {
        number: 60,
        ccp4: 0,
        hm: "P c n b",
        ext: 0,
        qualifier: "a-cb",
        hall: "-P 2b 2ac",
        basisop_idx: 5,
        ops: OPS_288,
    },
    BioSpaceGroup {
        number: 61,
        ccp4: 61,
        hm: "P b c a",
        ext: 0,
        qualifier: "",
        hall: "-P 2ac 2ab",
        basisop_idx: 0,
        ops: OPS_289,
    },
    BioSpaceGroup {
        number: 61,
        ccp4: 0,
        hm: "P c a b",
        ext: 0,
        qualifier: "ba-c",
        hall: "-P 2bc 2ac",
        basisop_idx: 3,
        ops: OPS_290,
    },
    BioSpaceGroup {
        number: 62,
        ccp4: 62,
        hm: "P n m a",
        ext: 0,
        qualifier: "",
        hall: "-P 2ac 2n",
        basisop_idx: 0,
        ops: OPS_291,
    },
    BioSpaceGroup {
        number: 62,
        ccp4: 0,
        hm: "P m n b",
        ext: 0,
        qualifier: "ba-c",
        hall: "-P 2bc 2a",
        basisop_idx: 7,
        ops: OPS_292,
    },
    BioSpaceGroup {
        number: 62,
        ccp4: 0,
        hm: "P b n m",
        ext: 0,
        qualifier: "cab",
        hall: "-P 2c 2ab",
        basisop_idx: 1,
        ops: OPS_293,
    },
    BioSpaceGroup {
        number: 62,
        ccp4: 0,
        hm: "P c m n",
        ext: 0,
        qualifier: "-cba",
        hall: "-P 2n 2ac",
        basisop_idx: 3,
        ops: OPS_294,
    },
    BioSpaceGroup {
        number: 62,
        ccp4: 0,
        hm: "P m c n",
        ext: 0,
        qualifier: "bca",
        hall: "-P 2n 2a",
        basisop_idx: 2,
        ops: OPS_295,
    },
    BioSpaceGroup {
        number: 62,
        ccp4: 0,
        hm: "P n a m",
        ext: 0,
        qualifier: "a-cb",
        hall: "-P 2c 2n",
        basisop_idx: 5,
        ops: OPS_296,
    },
    BioSpaceGroup {
        number: 63,
        ccp4: 63,
        hm: "C m c m",
        ext: 0,
        qualifier: "",
        hall: "-C 2c 2",
        basisop_idx: 0,
        ops: OPS_297,
    },
    BioSpaceGroup {
        number: 63,
        ccp4: 0,
        hm: "C c m m",
        ext: 0,
        qualifier: "ba-c",
        hall: "-C 2c 2c",
        basisop_idx: 7,
        ops: OPS_298,
    },
    BioSpaceGroup {
        number: 63,
        ccp4: 0,
        hm: "A m m a",
        ext: 0,
        qualifier: "cab",
        hall: "-A 2a 2a",
        basisop_idx: 1,
        ops: OPS_299,
    },
    BioSpaceGroup {
        number: 63,
        ccp4: 0,
        hm: "A m a m",
        ext: 0,
        qualifier: "-cba",
        hall: "-A 2 2a",
        basisop_idx: 3,
        ops: OPS_300,
    },
    BioSpaceGroup {
        number: 63,
        ccp4: 0,
        hm: "B b m m",
        ext: 0,
        qualifier: "bca",
        hall: "-B 2 2b",
        basisop_idx: 2,
        ops: OPS_301,
    },
    BioSpaceGroup {
        number: 63,
        ccp4: 0,
        hm: "B m m b",
        ext: 0,
        qualifier: "a-cb",
        hall: "-B 2b 2",
        basisop_idx: 5,
        ops: OPS_302,
    },
    BioSpaceGroup {
        number: 64,
        ccp4: 64,
        hm: "C m c a",
        ext: 0,
        qualifier: "",
        hall: "-C 2ac 2",
        basisop_idx: 0,
        ops: OPS_303,
    },
    BioSpaceGroup {
        number: 64,
        ccp4: 0,
        hm: "C c m b",
        ext: 0,
        qualifier: "ba-c",
        hall: "-C 2ac 2ac",
        basisop_idx: 7,
        ops: OPS_304,
    },
    BioSpaceGroup {
        number: 64,
        ccp4: 0,
        hm: "A b m a",
        ext: 0,
        qualifier: "cab",
        hall: "-A 2ab 2ab",
        basisop_idx: 1,
        ops: OPS_305,
    },
    BioSpaceGroup {
        number: 64,
        ccp4: 0,
        hm: "A c a m",
        ext: 0,
        qualifier: "-cba",
        hall: "-A 2 2ab",
        basisop_idx: 3,
        ops: OPS_306,
    },
    BioSpaceGroup {
        number: 64,
        ccp4: 0,
        hm: "B b c m",
        ext: 0,
        qualifier: "bca",
        hall: "-B 2 2ab",
        basisop_idx: 2,
        ops: OPS_307,
    },
    BioSpaceGroup {
        number: 64,
        ccp4: 0,
        hm: "B m a b",
        ext: 0,
        qualifier: "a-cb",
        hall: "-B 2ab 2",
        basisop_idx: 5,
        ops: OPS_308,
    },
    BioSpaceGroup {
        number: 65,
        ccp4: 65,
        hm: "C m m m",
        ext: 0,
        qualifier: "",
        hall: "-C 2 2",
        basisop_idx: 0,
        ops: OPS_309,
    },
    BioSpaceGroup {
        number: 65,
        ccp4: 0,
        hm: "A m m m",
        ext: 0,
        qualifier: "cab",
        hall: "-A 2 2",
        basisop_idx: 1,
        ops: OPS_310,
    },
    BioSpaceGroup {
        number: 65,
        ccp4: 0,
        hm: "B m m m",
        ext: 0,
        qualifier: "bca",
        hall: "-B 2 2",
        basisop_idx: 2,
        ops: OPS_311,
    },
    BioSpaceGroup {
        number: 66,
        ccp4: 66,
        hm: "C c c m",
        ext: 0,
        qualifier: "",
        hall: "-C 2 2c",
        basisop_idx: 0,
        ops: OPS_312,
    },
    BioSpaceGroup {
        number: 66,
        ccp4: 0,
        hm: "A m a a",
        ext: 0,
        qualifier: "cab",
        hall: "-A 2a 2",
        basisop_idx: 1,
        ops: OPS_313,
    },
    BioSpaceGroup {
        number: 66,
        ccp4: 0,
        hm: "B b m b",
        ext: 0,
        qualifier: "bca",
        hall: "-B 2b 2b",
        basisop_idx: 2,
        ops: OPS_314,
    },
    BioSpaceGroup {
        number: 67,
        ccp4: 67,
        hm: "C m m a",
        ext: 0,
        qualifier: "",
        hall: "-C 2a 2",
        basisop_idx: 0,
        ops: OPS_315,
    },
    BioSpaceGroup {
        number: 67,
        ccp4: 0,
        hm: "C m m b",
        ext: 0,
        qualifier: "ba-c",
        hall: "-C 2a 2a",
        basisop_idx: 14,
        ops: OPS_316,
    },
    BioSpaceGroup {
        number: 67,
        ccp4: 0,
        hm: "A b m m",
        ext: 0,
        qualifier: "cab",
        hall: "-A 2b 2b",
        basisop_idx: 1,
        ops: OPS_317,
    },
    BioSpaceGroup {
        number: 67,
        ccp4: 0,
        hm: "A c m m",
        ext: 0,
        qualifier: "-cba",
        hall: "-A 2 2b",
        basisop_idx: 3,
        ops: OPS_318,
    },
    BioSpaceGroup {
        number: 67,
        ccp4: 0,
        hm: "B m c m",
        ext: 0,
        qualifier: "bca",
        hall: "-B 2 2a",
        basisop_idx: 2,
        ops: OPS_319,
    },
    BioSpaceGroup {
        number: 67,
        ccp4: 0,
        hm: "B m a m",
        ext: 0,
        qualifier: "a-cb",
        hall: "-B 2a 2",
        basisop_idx: 5,
        ops: OPS_320,
    },
    BioSpaceGroup {
        number: 68,
        ccp4: 68,
        hm: "C c c a",
        ext: 49,
        qualifier: "",
        hall: "C 2 2 -1ac",
        basisop_idx: 24,
        ops: OPS_321,
    },
    BioSpaceGroup {
        number: 68,
        ccp4: 0,
        hm: "C c c a",
        ext: 50,
        qualifier: "",
        hall: "-C 2a 2ac",
        basisop_idx: 0,
        ops: OPS_322,
    },
    BioSpaceGroup {
        number: 68,
        ccp4: 0,
        hm: "C c c b",
        ext: 49,
        qualifier: "ba-c",
        hall: "C 2 2 -1ac",
        basisop_idx: 24,
        ops: OPS_323,
    },
    BioSpaceGroup {
        number: 68,
        ccp4: 0,
        hm: "C c c b",
        ext: 50,
        qualifier: "ba-c",
        hall: "-C 2a 2c",
        basisop_idx: 21,
        ops: OPS_324,
    },
    BioSpaceGroup {
        number: 68,
        ccp4: 0,
        hm: "A b a a",
        ext: 49,
        qualifier: "cab",
        hall: "A 2 2 -1ab",
        basisop_idx: 25,
        ops: OPS_325,
    },
    BioSpaceGroup {
        number: 68,
        ccp4: 0,
        hm: "A b a a",
        ext: 50,
        qualifier: "cab",
        hall: "-A 2a 2b",
        basisop_idx: 1,
        ops: OPS_326,
    },
    BioSpaceGroup {
        number: 68,
        ccp4: 0,
        hm: "A c a a",
        ext: 49,
        qualifier: "-cba",
        hall: "A 2 2 -1ab",
        basisop_idx: 25,
        ops: OPS_327,
    },
    BioSpaceGroup {
        number: 68,
        ccp4: 0,
        hm: "A c a a",
        ext: 50,
        qualifier: "-cba",
        hall: "-A 2ab 2b",
        basisop_idx: 3,
        ops: OPS_328,
    },
    BioSpaceGroup {
        number: 68,
        ccp4: 0,
        hm: "B b c b",
        ext: 49,
        qualifier: "bca",
        hall: "B 2 2 -1ab",
        basisop_idx: 26,
        ops: OPS_329,
    },
    BioSpaceGroup {
        number: 68,
        ccp4: 0,
        hm: "B b c b",
        ext: 50,
        qualifier: "bca",
        hall: "-B 2ab 2b",
        basisop_idx: 2,
        ops: OPS_330,
    },
    BioSpaceGroup {
        number: 68,
        ccp4: 0,
        hm: "B b a b",
        ext: 49,
        qualifier: "a-cb",
        hall: "B 2 2 -1ab",
        basisop_idx: 26,
        ops: OPS_331,
    },
    BioSpaceGroup {
        number: 68,
        ccp4: 0,
        hm: "B b a b",
        ext: 50,
        qualifier: "a-cb",
        hall: "-B 2b 2ab",
        basisop_idx: 5,
        ops: OPS_332,
    },
    BioSpaceGroup {
        number: 69,
        ccp4: 69,
        hm: "F m m m",
        ext: 0,
        qualifier: "",
        hall: "-F 2 2",
        basisop_idx: 0,
        ops: OPS_333,
    },
    BioSpaceGroup {
        number: 70,
        ccp4: 70,
        hm: "F d d d",
        ext: 49,
        qualifier: "",
        hall: "F 2 2 -1d",
        basisop_idx: 27,
        ops: OPS_334,
    },
    BioSpaceGroup {
        number: 70,
        ccp4: 0,
        hm: "F d d d",
        ext: 50,
        qualifier: "",
        hall: "-F 2uv 2vw",
        basisop_idx: 0,
        ops: OPS_335,
    },
    BioSpaceGroup {
        number: 71,
        ccp4: 71,
        hm: "I m m m",
        ext: 0,
        qualifier: "",
        hall: "-I 2 2",
        basisop_idx: 0,
        ops: OPS_336,
    },
    BioSpaceGroup {
        number: 72,
        ccp4: 72,
        hm: "I b a m",
        ext: 0,
        qualifier: "",
        hall: "-I 2 2c",
        basisop_idx: 0,
        ops: OPS_337,
    },
    BioSpaceGroup {
        number: 72,
        ccp4: 0,
        hm: "I m c b",
        ext: 0,
        qualifier: "cab",
        hall: "-I 2a 2",
        basisop_idx: 1,
        ops: OPS_338,
    },
    BioSpaceGroup {
        number: 72,
        ccp4: 0,
        hm: "I c m a",
        ext: 0,
        qualifier: "bca",
        hall: "-I 2b 2b",
        basisop_idx: 2,
        ops: OPS_339,
    },
    BioSpaceGroup {
        number: 73,
        ccp4: 73,
        hm: "I b c a",
        ext: 0,
        qualifier: "",
        hall: "-I 2b 2c",
        basisop_idx: 0,
        ops: OPS_340,
    },
    BioSpaceGroup {
        number: 73,
        ccp4: 0,
        hm: "I c a b",
        ext: 0,
        qualifier: "ba-c",
        hall: "-I 2a 2b",
        basisop_idx: 28,
        ops: OPS_341,
    },
    BioSpaceGroup {
        number: 74,
        ccp4: 74,
        hm: "I m m a",
        ext: 0,
        qualifier: "",
        hall: "-I 2b 2",
        basisop_idx: 0,
        ops: OPS_342,
    },
    BioSpaceGroup {
        number: 74,
        ccp4: 0,
        hm: "I m m b",
        ext: 0,
        qualifier: "ba-c",
        hall: "-I 2a 2a",
        basisop_idx: 28,
        ops: OPS_343,
    },
    BioSpaceGroup {
        number: 74,
        ccp4: 0,
        hm: "I b m m",
        ext: 0,
        qualifier: "cab",
        hall: "-I 2c 2c",
        basisop_idx: 1,
        ops: OPS_344,
    },
    BioSpaceGroup {
        number: 74,
        ccp4: 0,
        hm: "I c m m",
        ext: 0,
        qualifier: "-cba",
        hall: "-I 2 2b",
        basisop_idx: 3,
        ops: OPS_345,
    },
    BioSpaceGroup {
        number: 74,
        ccp4: 0,
        hm: "I m c m",
        ext: 0,
        qualifier: "bca",
        hall: "-I 2 2a",
        basisop_idx: 2,
        ops: OPS_346,
    },
    BioSpaceGroup {
        number: 74,
        ccp4: 0,
        hm: "I m a m",
        ext: 0,
        qualifier: "a-cb",
        hall: "-I 2c 2",
        basisop_idx: 5,
        ops: OPS_347,
    },
    BioSpaceGroup {
        number: 75,
        ccp4: 75,
        hm: "P 4",
        ext: 0,
        qualifier: "",
        hall: "P 4",
        basisop_idx: 0,
        ops: OPS_348,
    },
    BioSpaceGroup {
        number: 76,
        ccp4: 76,
        hm: "P 41",
        ext: 0,
        qualifier: "",
        hall: "P 4w",
        basisop_idx: 0,
        ops: OPS_349,
    },
    BioSpaceGroup {
        number: 77,
        ccp4: 77,
        hm: "P 42",
        ext: 0,
        qualifier: "",
        hall: "P 4c",
        basisop_idx: 0,
        ops: OPS_350,
    },
    BioSpaceGroup {
        number: 78,
        ccp4: 78,
        hm: "P 43",
        ext: 0,
        qualifier: "",
        hall: "P 4cw",
        basisop_idx: 0,
        ops: OPS_351,
    },
    BioSpaceGroup {
        number: 79,
        ccp4: 79,
        hm: "I 4",
        ext: 0,
        qualifier: "",
        hall: "I 4",
        basisop_idx: 0,
        ops: OPS_352,
    },
    BioSpaceGroup {
        number: 80,
        ccp4: 80,
        hm: "I 41",
        ext: 0,
        qualifier: "",
        hall: "I 4bw",
        basisop_idx: 0,
        ops: OPS_353,
    },
    BioSpaceGroup {
        number: 81,
        ccp4: 81,
        hm: "P -4",
        ext: 0,
        qualifier: "",
        hall: "P -4",
        basisop_idx: 0,
        ops: OPS_354,
    },
    BioSpaceGroup {
        number: 82,
        ccp4: 82,
        hm: "I -4",
        ext: 0,
        qualifier: "",
        hall: "I -4",
        basisop_idx: 0,
        ops: OPS_355,
    },
    BioSpaceGroup {
        number: 83,
        ccp4: 83,
        hm: "P 4/m",
        ext: 0,
        qualifier: "",
        hall: "-P 4",
        basisop_idx: 0,
        ops: OPS_356,
    },
    BioSpaceGroup {
        number: 84,
        ccp4: 84,
        hm: "P 42/m",
        ext: 0,
        qualifier: "",
        hall: "-P 4c",
        basisop_idx: 0,
        ops: OPS_357,
    },
    BioSpaceGroup {
        number: 85,
        ccp4: 85,
        hm: "P 4/n",
        ext: 49,
        qualifier: "",
        hall: "P 4ab -1ab",
        basisop_idx: 29,
        ops: OPS_358,
    },
    BioSpaceGroup {
        number: 85,
        ccp4: 0,
        hm: "P 4/n",
        ext: 50,
        qualifier: "",
        hall: "-P 4a",
        basisop_idx: 0,
        ops: OPS_359,
    },
    BioSpaceGroup {
        number: 86,
        ccp4: 86,
        hm: "P 42/n",
        ext: 49,
        qualifier: "",
        hall: "P 4n -1n",
        basisop_idx: 30,
        ops: OPS_360,
    },
    BioSpaceGroup {
        number: 86,
        ccp4: 0,
        hm: "P 42/n",
        ext: 50,
        qualifier: "",
        hall: "-P 4bc",
        basisop_idx: 0,
        ops: OPS_361,
    },
    BioSpaceGroup {
        number: 87,
        ccp4: 87,
        hm: "I 4/m",
        ext: 0,
        qualifier: "",
        hall: "-I 4",
        basisop_idx: 0,
        ops: OPS_362,
    },
    BioSpaceGroup {
        number: 88,
        ccp4: 88,
        hm: "I 41/a",
        ext: 49,
        qualifier: "",
        hall: "I 4bw -1bw",
        basisop_idx: 31,
        ops: OPS_363,
    },
    BioSpaceGroup {
        number: 88,
        ccp4: 0,
        hm: "I 41/a",
        ext: 50,
        qualifier: "",
        hall: "-I 4ad",
        basisop_idx: 0,
        ops: OPS_364,
    },
    BioSpaceGroup {
        number: 89,
        ccp4: 89,
        hm: "P 4 2 2",
        ext: 0,
        qualifier: "",
        hall: "P 4 2",
        basisop_idx: 0,
        ops: OPS_365,
    },
    BioSpaceGroup {
        number: 90,
        ccp4: 90,
        hm: "P 4 21 2",
        ext: 0,
        qualifier: "",
        hall: "P 4ab 2ab",
        basisop_idx: 0,
        ops: OPS_366,
    },
    BioSpaceGroup {
        number: 91,
        ccp4: 91,
        hm: "P 41 2 2",
        ext: 0,
        qualifier: "",
        hall: "P 4w 2c",
        basisop_idx: 0,
        ops: OPS_367,
    },
    BioSpaceGroup {
        number: 92,
        ccp4: 92,
        hm: "P 41 21 2",
        ext: 0,
        qualifier: "",
        hall: "P 4abw 2nw",
        basisop_idx: 0,
        ops: OPS_368,
    },
    BioSpaceGroup {
        number: 93,
        ccp4: 93,
        hm: "P 42 2 2",
        ext: 0,
        qualifier: "",
        hall: "P 4c 2",
        basisop_idx: 0,
        ops: OPS_369,
    },
    BioSpaceGroup {
        number: 94,
        ccp4: 94,
        hm: "P 42 21 2",
        ext: 0,
        qualifier: "",
        hall: "P 4n 2n",
        basisop_idx: 0,
        ops: OPS_370,
    },
    BioSpaceGroup {
        number: 95,
        ccp4: 95,
        hm: "P 43 2 2",
        ext: 0,
        qualifier: "",
        hall: "P 4cw 2c",
        basisop_idx: 0,
        ops: OPS_371,
    },
    BioSpaceGroup {
        number: 96,
        ccp4: 96,
        hm: "P 43 21 2",
        ext: 0,
        qualifier: "",
        hall: "P 4nw 2abw",
        basisop_idx: 0,
        ops: OPS_372,
    },
    BioSpaceGroup {
        number: 97,
        ccp4: 97,
        hm: "I 4 2 2",
        ext: 0,
        qualifier: "",
        hall: "I 4 2",
        basisop_idx: 0,
        ops: OPS_373,
    },
    BioSpaceGroup {
        number: 98,
        ccp4: 98,
        hm: "I 41 2 2",
        ext: 0,
        qualifier: "",
        hall: "I 4bw 2bw",
        basisop_idx: 0,
        ops: OPS_374,
    },
    BioSpaceGroup {
        number: 99,
        ccp4: 99,
        hm: "P 4 m m",
        ext: 0,
        qualifier: "",
        hall: "P 4 -2",
        basisop_idx: 0,
        ops: OPS_375,
    },
    BioSpaceGroup {
        number: 100,
        ccp4: 100,
        hm: "P 4 b m",
        ext: 0,
        qualifier: "",
        hall: "P 4 -2ab",
        basisop_idx: 0,
        ops: OPS_376,
    },
    BioSpaceGroup {
        number: 101,
        ccp4: 101,
        hm: "P 42 c m",
        ext: 0,
        qualifier: "",
        hall: "P 4c -2c",
        basisop_idx: 0,
        ops: OPS_377,
    },
    BioSpaceGroup {
        number: 102,
        ccp4: 102,
        hm: "P 42 n m",
        ext: 0,
        qualifier: "",
        hall: "P 4n -2n",
        basisop_idx: 0,
        ops: OPS_378,
    },
    BioSpaceGroup {
        number: 103,
        ccp4: 103,
        hm: "P 4 c c",
        ext: 0,
        qualifier: "",
        hall: "P 4 -2c",
        basisop_idx: 0,
        ops: OPS_379,
    },
    BioSpaceGroup {
        number: 104,
        ccp4: 104,
        hm: "P 4 n c",
        ext: 0,
        qualifier: "",
        hall: "P 4 -2n",
        basisop_idx: 0,
        ops: OPS_380,
    },
    BioSpaceGroup {
        number: 105,
        ccp4: 105,
        hm: "P 42 m c",
        ext: 0,
        qualifier: "",
        hall: "P 4c -2",
        basisop_idx: 0,
        ops: OPS_381,
    },
    BioSpaceGroup {
        number: 106,
        ccp4: 106,
        hm: "P 42 b c",
        ext: 0,
        qualifier: "",
        hall: "P 4c -2ab",
        basisop_idx: 0,
        ops: OPS_382,
    },
    BioSpaceGroup {
        number: 107,
        ccp4: 107,
        hm: "I 4 m m",
        ext: 0,
        qualifier: "",
        hall: "I 4 -2",
        basisop_idx: 0,
        ops: OPS_383,
    },
    BioSpaceGroup {
        number: 108,
        ccp4: 108,
        hm: "I 4 c m",
        ext: 0,
        qualifier: "",
        hall: "I 4 -2c",
        basisop_idx: 0,
        ops: OPS_384,
    },
    BioSpaceGroup {
        number: 109,
        ccp4: 109,
        hm: "I 41 m d",
        ext: 0,
        qualifier: "",
        hall: "I 4bw -2",
        basisop_idx: 0,
        ops: OPS_385,
    },
    BioSpaceGroup {
        number: 110,
        ccp4: 110,
        hm: "I 41 c d",
        ext: 0,
        qualifier: "",
        hall: "I 4bw -2c",
        basisop_idx: 0,
        ops: OPS_386,
    },
    BioSpaceGroup {
        number: 111,
        ccp4: 111,
        hm: "P -4 2 m",
        ext: 0,
        qualifier: "",
        hall: "P -4 2",
        basisop_idx: 0,
        ops: OPS_387,
    },
    BioSpaceGroup {
        number: 112,
        ccp4: 112,
        hm: "P -4 2 c",
        ext: 0,
        qualifier: "",
        hall: "P -4 2c",
        basisop_idx: 0,
        ops: OPS_388,
    },
    BioSpaceGroup {
        number: 113,
        ccp4: 113,
        hm: "P -4 21 m",
        ext: 0,
        qualifier: "",
        hall: "P -4 2ab",
        basisop_idx: 0,
        ops: OPS_389,
    },
    BioSpaceGroup {
        number: 114,
        ccp4: 114,
        hm: "P -4 21 c",
        ext: 0,
        qualifier: "",
        hall: "P -4 2n",
        basisop_idx: 0,
        ops: OPS_390,
    },
    BioSpaceGroup {
        number: 115,
        ccp4: 115,
        hm: "P -4 m 2",
        ext: 0,
        qualifier: "",
        hall: "P -4 -2",
        basisop_idx: 0,
        ops: OPS_391,
    },
    BioSpaceGroup {
        number: 116,
        ccp4: 116,
        hm: "P -4 c 2",
        ext: 0,
        qualifier: "",
        hall: "P -4 -2c",
        basisop_idx: 0,
        ops: OPS_392,
    },
    BioSpaceGroup {
        number: 117,
        ccp4: 117,
        hm: "P -4 b 2",
        ext: 0,
        qualifier: "",
        hall: "P -4 -2ab",
        basisop_idx: 0,
        ops: OPS_393,
    },
    BioSpaceGroup {
        number: 118,
        ccp4: 118,
        hm: "P -4 n 2",
        ext: 0,
        qualifier: "",
        hall: "P -4 -2n",
        basisop_idx: 0,
        ops: OPS_394,
    },
    BioSpaceGroup {
        number: 119,
        ccp4: 119,
        hm: "I -4 m 2",
        ext: 0,
        qualifier: "",
        hall: "I -4 -2",
        basisop_idx: 0,
        ops: OPS_395,
    },
    BioSpaceGroup {
        number: 120,
        ccp4: 120,
        hm: "I -4 c 2",
        ext: 0,
        qualifier: "",
        hall: "I -4 -2c",
        basisop_idx: 0,
        ops: OPS_396,
    },
    BioSpaceGroup {
        number: 121,
        ccp4: 121,
        hm: "I -4 2 m",
        ext: 0,
        qualifier: "",
        hall: "I -4 2",
        basisop_idx: 0,
        ops: OPS_397,
    },
    BioSpaceGroup {
        number: 122,
        ccp4: 122,
        hm: "I -4 2 d",
        ext: 0,
        qualifier: "",
        hall: "I -4 2bw",
        basisop_idx: 0,
        ops: OPS_398,
    },
    BioSpaceGroup {
        number: 123,
        ccp4: 123,
        hm: "P 4/m m m",
        ext: 0,
        qualifier: "",
        hall: "-P 4 2",
        basisop_idx: 0,
        ops: OPS_399,
    },
    BioSpaceGroup {
        number: 124,
        ccp4: 124,
        hm: "P 4/m c c",
        ext: 0,
        qualifier: "",
        hall: "-P 4 2c",
        basisop_idx: 0,
        ops: OPS_400,
    },
    BioSpaceGroup {
        number: 125,
        ccp4: 125,
        hm: "P 4/n b m",
        ext: 49,
        qualifier: "",
        hall: "P 4 2 -1ab",
        basisop_idx: 21,
        ops: OPS_401,
    },
    BioSpaceGroup {
        number: 125,
        ccp4: 0,
        hm: "P 4/n b m",
        ext: 50,
        qualifier: "",
        hall: "-P 4a 2b",
        basisop_idx: 0,
        ops: OPS_402,
    },
    BioSpaceGroup {
        number: 126,
        ccp4: 126,
        hm: "P 4/n n c",
        ext: 49,
        qualifier: "",
        hall: "P 4 2 -1n",
        basisop_idx: 20,
        ops: OPS_403,
    },
    BioSpaceGroup {
        number: 126,
        ccp4: 0,
        hm: "P 4/n n c",
        ext: 50,
        qualifier: "",
        hall: "-P 4a 2bc",
        basisop_idx: 0,
        ops: OPS_404,
    },
    BioSpaceGroup {
        number: 127,
        ccp4: 127,
        hm: "P 4/m b m",
        ext: 0,
        qualifier: "",
        hall: "-P 4 2ab",
        basisop_idx: 0,
        ops: OPS_405,
    },
    BioSpaceGroup {
        number: 128,
        ccp4: 128,
        hm: "P 4/m n c",
        ext: 0,
        qualifier: "",
        hall: "-P 4 2n",
        basisop_idx: 0,
        ops: OPS_406,
    },
    BioSpaceGroup {
        number: 129,
        ccp4: 129,
        hm: "P 4/n m m",
        ext: 49,
        qualifier: "",
        hall: "P 4ab 2ab -1ab",
        basisop_idx: 29,
        ops: OPS_407,
    },
    BioSpaceGroup {
        number: 129,
        ccp4: 0,
        hm: "P 4/n m m",
        ext: 50,
        qualifier: "",
        hall: "-P 4a 2a",
        basisop_idx: 0,
        ops: OPS_408,
    },
    BioSpaceGroup {
        number: 130,
        ccp4: 130,
        hm: "P 4/n c c",
        ext: 49,
        qualifier: "",
        hall: "P 4ab 2n -1ab",
        basisop_idx: 29,
        ops: OPS_409,
    },
    BioSpaceGroup {
        number: 130,
        ccp4: 0,
        hm: "P 4/n c c",
        ext: 50,
        qualifier: "",
        hall: "-P 4a 2ac",
        basisop_idx: 0,
        ops: OPS_410,
    },
    BioSpaceGroup {
        number: 131,
        ccp4: 131,
        hm: "P 42/m m c",
        ext: 0,
        qualifier: "",
        hall: "-P 4c 2",
        basisop_idx: 0,
        ops: OPS_411,
    },
    BioSpaceGroup {
        number: 132,
        ccp4: 132,
        hm: "P 42/m c m",
        ext: 0,
        qualifier: "",
        hall: "-P 4c 2c",
        basisop_idx: 0,
        ops: OPS_412,
    },
    BioSpaceGroup {
        number: 133,
        ccp4: 133,
        hm: "P 42/n b c",
        ext: 49,
        qualifier: "",
        hall: "P 4n 2c -1n",
        basisop_idx: 32,
        ops: OPS_413,
    },
    BioSpaceGroup {
        number: 133,
        ccp4: 0,
        hm: "P 42/n b c",
        ext: 50,
        qualifier: "",
        hall: "-P 4ac 2b",
        basisop_idx: 0,
        ops: OPS_414,
    },
    BioSpaceGroup {
        number: 134,
        ccp4: 134,
        hm: "P 42/n n m",
        ext: 49,
        qualifier: "",
        hall: "P 4n 2 -1n",
        basisop_idx: 33,
        ops: OPS_415,
    },
    BioSpaceGroup {
        number: 134,
        ccp4: 0,
        hm: "P 42/n n m",
        ext: 50,
        qualifier: "",
        hall: "-P 4ac 2bc",
        basisop_idx: 0,
        ops: OPS_416,
    },
    BioSpaceGroup {
        number: 135,
        ccp4: 135,
        hm: "P 42/m b c",
        ext: 0,
        qualifier: "",
        hall: "-P 4c 2ab",
        basisop_idx: 0,
        ops: OPS_417,
    },
    BioSpaceGroup {
        number: 136,
        ccp4: 136,
        hm: "P 42/m n m",
        ext: 0,
        qualifier: "",
        hall: "-P 4n 2n",
        basisop_idx: 0,
        ops: OPS_418,
    },
    BioSpaceGroup {
        number: 137,
        ccp4: 137,
        hm: "P 42/n m c",
        ext: 49,
        qualifier: "",
        hall: "P 4n 2n -1n",
        basisop_idx: 32,
        ops: OPS_419,
    },
    BioSpaceGroup {
        number: 137,
        ccp4: 0,
        hm: "P 42/n m c",
        ext: 50,
        qualifier: "",
        hall: "-P 4ac 2a",
        basisop_idx: 0,
        ops: OPS_420,
    },
    BioSpaceGroup {
        number: 138,
        ccp4: 138,
        hm: "P 42/n c m",
        ext: 49,
        qualifier: "",
        hall: "P 4n 2ab -1n",
        basisop_idx: 33,
        ops: OPS_421,
    },
    BioSpaceGroup {
        number: 138,
        ccp4: 0,
        hm: "P 42/n c m",
        ext: 50,
        qualifier: "",
        hall: "-P 4ac 2ac",
        basisop_idx: 0,
        ops: OPS_422,
    },
    BioSpaceGroup {
        number: 139,
        ccp4: 139,
        hm: "I 4/m m m",
        ext: 0,
        qualifier: "",
        hall: "-I 4 2",
        basisop_idx: 0,
        ops: OPS_423,
    },
    BioSpaceGroup {
        number: 140,
        ccp4: 140,
        hm: "I 4/m c m",
        ext: 0,
        qualifier: "",
        hall: "-I 4 2c",
        basisop_idx: 0,
        ops: OPS_424,
    },
    BioSpaceGroup {
        number: 141,
        ccp4: 141,
        hm: "I 41/a m d",
        ext: 49,
        qualifier: "",
        hall: "I 4bw 2bw -1bw",
        basisop_idx: 34,
        ops: OPS_425,
    },
    BioSpaceGroup {
        number: 141,
        ccp4: 0,
        hm: "I 41/a m d",
        ext: 50,
        qualifier: "",
        hall: "-I 4bd 2",
        basisop_idx: 0,
        ops: OPS_426,
    },
    BioSpaceGroup {
        number: 142,
        ccp4: 142,
        hm: "I 41/a c d",
        ext: 49,
        qualifier: "",
        hall: "I 4bw 2aw -1bw",
        basisop_idx: 35,
        ops: OPS_427,
    },
    BioSpaceGroup {
        number: 142,
        ccp4: 0,
        hm: "I 41/a c d",
        ext: 50,
        qualifier: "",
        hall: "-I 4bd 2c",
        basisop_idx: 0,
        ops: OPS_428,
    },
    BioSpaceGroup {
        number: 143,
        ccp4: 143,
        hm: "P 3",
        ext: 0,
        qualifier: "",
        hall: "P 3",
        basisop_idx: 0,
        ops: OPS_429,
    },
    BioSpaceGroup {
        number: 144,
        ccp4: 144,
        hm: "P 31",
        ext: 0,
        qualifier: "",
        hall: "P 31",
        basisop_idx: 0,
        ops: OPS_430,
    },
    BioSpaceGroup {
        number: 145,
        ccp4: 145,
        hm: "P 32",
        ext: 0,
        qualifier: "",
        hall: "P 32",
        basisop_idx: 0,
        ops: OPS_431,
    },
    BioSpaceGroup {
        number: 146,
        ccp4: 146,
        hm: "R 3",
        ext: 72,
        qualifier: "",
        hall: "R 3",
        basisop_idx: 0,
        ops: OPS_432,
    },
    BioSpaceGroup {
        number: 146,
        ccp4: 1146,
        hm: "R 3",
        ext: 82,
        qualifier: "",
        hall: "P 3*",
        basisop_idx: 36,
        ops: OPS_433,
    },
    BioSpaceGroup {
        number: 147,
        ccp4: 147,
        hm: "P -3",
        ext: 0,
        qualifier: "",
        hall: "-P 3",
        basisop_idx: 0,
        ops: OPS_434,
    },
    BioSpaceGroup {
        number: 148,
        ccp4: 148,
        hm: "R -3",
        ext: 72,
        qualifier: "",
        hall: "-R 3",
        basisop_idx: 0,
        ops: OPS_435,
    },
    BioSpaceGroup {
        number: 148,
        ccp4: 1148,
        hm: "R -3",
        ext: 82,
        qualifier: "",
        hall: "-P 3*",
        basisop_idx: 36,
        ops: OPS_436,
    },
    BioSpaceGroup {
        number: 149,
        ccp4: 149,
        hm: "P 3 1 2",
        ext: 0,
        qualifier: "",
        hall: "P 3 2",
        basisop_idx: 0,
        ops: OPS_437,
    },
    BioSpaceGroup {
        number: 150,
        ccp4: 150,
        hm: "P 3 2 1",
        ext: 0,
        qualifier: "",
        hall: "P 3 2\"",
        basisop_idx: 0,
        ops: OPS_438,
    },
    BioSpaceGroup {
        number: 151,
        ccp4: 151,
        hm: "P 31 1 2",
        ext: 0,
        qualifier: "",
        hall: "P 31 2 (0 0 4)",
        basisop_idx: 0,
        ops: OPS_439,
    },
    BioSpaceGroup {
        number: 152,
        ccp4: 152,
        hm: "P 31 2 1",
        ext: 0,
        qualifier: "",
        hall: "P 31 2\"",
        basisop_idx: 0,
        ops: OPS_440,
    },
    BioSpaceGroup {
        number: 153,
        ccp4: 153,
        hm: "P 32 1 2",
        ext: 0,
        qualifier: "",
        hall: "P 32 2 (0 0 2)",
        basisop_idx: 0,
        ops: OPS_441,
    },
    BioSpaceGroup {
        number: 154,
        ccp4: 154,
        hm: "P 32 2 1",
        ext: 0,
        qualifier: "",
        hall: "P 32 2\"",
        basisop_idx: 0,
        ops: OPS_442,
    },
    BioSpaceGroup {
        number: 155,
        ccp4: 155,
        hm: "R 3 2",
        ext: 72,
        qualifier: "",
        hall: "R 3 2\"",
        basisop_idx: 0,
        ops: OPS_443,
    },
    BioSpaceGroup {
        number: 155,
        ccp4: 1155,
        hm: "R 3 2",
        ext: 82,
        qualifier: "",
        hall: "P 3* 2",
        basisop_idx: 36,
        ops: OPS_444,
    },
    BioSpaceGroup {
        number: 156,
        ccp4: 156,
        hm: "P 3 m 1",
        ext: 0,
        qualifier: "",
        hall: "P 3 -2\"",
        basisop_idx: 0,
        ops: OPS_445,
    },
    BioSpaceGroup {
        number: 157,
        ccp4: 157,
        hm: "P 3 1 m",
        ext: 0,
        qualifier: "",
        hall: "P 3 -2",
        basisop_idx: 0,
        ops: OPS_446,
    },
    BioSpaceGroup {
        number: 158,
        ccp4: 158,
        hm: "P 3 c 1",
        ext: 0,
        qualifier: "",
        hall: "P 3 -2\"c",
        basisop_idx: 0,
        ops: OPS_447,
    },
    BioSpaceGroup {
        number: 159,
        ccp4: 159,
        hm: "P 3 1 c",
        ext: 0,
        qualifier: "",
        hall: "P 3 -2c",
        basisop_idx: 0,
        ops: OPS_448,
    },
    BioSpaceGroup {
        number: 160,
        ccp4: 160,
        hm: "R 3 m",
        ext: 72,
        qualifier: "",
        hall: "R 3 -2\"",
        basisop_idx: 0,
        ops: OPS_449,
    },
    BioSpaceGroup {
        number: 160,
        ccp4: 1160,
        hm: "R 3 m",
        ext: 82,
        qualifier: "",
        hall: "P 3* -2",
        basisop_idx: 36,
        ops: OPS_450,
    },
    BioSpaceGroup {
        number: 161,
        ccp4: 161,
        hm: "R 3 c",
        ext: 72,
        qualifier: "",
        hall: "R 3 -2\"c",
        basisop_idx: 0,
        ops: OPS_451,
    },
    BioSpaceGroup {
        number: 161,
        ccp4: 1161,
        hm: "R 3 c",
        ext: 82,
        qualifier: "",
        hall: "P 3* -2n",
        basisop_idx: 36,
        ops: OPS_452,
    },
    BioSpaceGroup {
        number: 162,
        ccp4: 162,
        hm: "P -3 1 m",
        ext: 0,
        qualifier: "",
        hall: "-P 3 2",
        basisop_idx: 0,
        ops: OPS_453,
    },
    BioSpaceGroup {
        number: 163,
        ccp4: 163,
        hm: "P -3 1 c",
        ext: 0,
        qualifier: "",
        hall: "-P 3 2c",
        basisop_idx: 0,
        ops: OPS_454,
    },
    BioSpaceGroup {
        number: 164,
        ccp4: 164,
        hm: "P -3 m 1",
        ext: 0,
        qualifier: "",
        hall: "-P 3 2\"",
        basisop_idx: 0,
        ops: OPS_455,
    },
    BioSpaceGroup {
        number: 165,
        ccp4: 165,
        hm: "P -3 c 1",
        ext: 0,
        qualifier: "",
        hall: "-P 3 2\"c",
        basisop_idx: 0,
        ops: OPS_456,
    },
    BioSpaceGroup {
        number: 166,
        ccp4: 166,
        hm: "R -3 m",
        ext: 72,
        qualifier: "",
        hall: "-R 3 2\"",
        basisop_idx: 0,
        ops: OPS_457,
    },
    BioSpaceGroup {
        number: 166,
        ccp4: 1166,
        hm: "R -3 m",
        ext: 82,
        qualifier: "",
        hall: "-P 3* 2",
        basisop_idx: 36,
        ops: OPS_458,
    },
    BioSpaceGroup {
        number: 167,
        ccp4: 167,
        hm: "R -3 c",
        ext: 72,
        qualifier: "",
        hall: "-R 3 2\"c",
        basisop_idx: 0,
        ops: OPS_459,
    },
    BioSpaceGroup {
        number: 167,
        ccp4: 1167,
        hm: "R -3 c",
        ext: 82,
        qualifier: "",
        hall: "-P 3* 2n",
        basisop_idx: 36,
        ops: OPS_460,
    },
    BioSpaceGroup {
        number: 168,
        ccp4: 168,
        hm: "P 6",
        ext: 0,
        qualifier: "",
        hall: "P 6",
        basisop_idx: 0,
        ops: OPS_461,
    },
    BioSpaceGroup {
        number: 169,
        ccp4: 169,
        hm: "P 61",
        ext: 0,
        qualifier: "",
        hall: "P 61",
        basisop_idx: 0,
        ops: OPS_462,
    },
    BioSpaceGroup {
        number: 170,
        ccp4: 170,
        hm: "P 65",
        ext: 0,
        qualifier: "",
        hall: "P 65",
        basisop_idx: 0,
        ops: OPS_463,
    },
    BioSpaceGroup {
        number: 171,
        ccp4: 171,
        hm: "P 62",
        ext: 0,
        qualifier: "",
        hall: "P 62",
        basisop_idx: 0,
        ops: OPS_464,
    },
    BioSpaceGroup {
        number: 172,
        ccp4: 172,
        hm: "P 64",
        ext: 0,
        qualifier: "",
        hall: "P 64",
        basisop_idx: 0,
        ops: OPS_465,
    },
    BioSpaceGroup {
        number: 173,
        ccp4: 173,
        hm: "P 63",
        ext: 0,
        qualifier: "",
        hall: "P 6c",
        basisop_idx: 0,
        ops: OPS_466,
    },
    BioSpaceGroup {
        number: 174,
        ccp4: 174,
        hm: "P -6",
        ext: 0,
        qualifier: "",
        hall: "P -6",
        basisop_idx: 0,
        ops: OPS_467,
    },
    BioSpaceGroup {
        number: 175,
        ccp4: 175,
        hm: "P 6/m",
        ext: 0,
        qualifier: "",
        hall: "-P 6",
        basisop_idx: 0,
        ops: OPS_468,
    },
    BioSpaceGroup {
        number: 176,
        ccp4: 176,
        hm: "P 63/m",
        ext: 0,
        qualifier: "",
        hall: "-P 6c",
        basisop_idx: 0,
        ops: OPS_469,
    },
    BioSpaceGroup {
        number: 177,
        ccp4: 177,
        hm: "P 6 2 2",
        ext: 0,
        qualifier: "",
        hall: "P 6 2",
        basisop_idx: 0,
        ops: OPS_470,
    },
    BioSpaceGroup {
        number: 178,
        ccp4: 178,
        hm: "P 61 2 2",
        ext: 0,
        qualifier: "",
        hall: "P 61 2 (0 0 5)",
        basisop_idx: 0,
        ops: OPS_471,
    },
    BioSpaceGroup {
        number: 179,
        ccp4: 179,
        hm: "P 65 2 2",
        ext: 0,
        qualifier: "",
        hall: "P 65 2 (0 0 1)",
        basisop_idx: 0,
        ops: OPS_472,
    },
    BioSpaceGroup {
        number: 180,
        ccp4: 180,
        hm: "P 62 2 2",
        ext: 0,
        qualifier: "",
        hall: "P 62 2 (0 0 4)",
        basisop_idx: 0,
        ops: OPS_473,
    },
    BioSpaceGroup {
        number: 181,
        ccp4: 181,
        hm: "P 64 2 2",
        ext: 0,
        qualifier: "",
        hall: "P 64 2 (0 0 2)",
        basisop_idx: 0,
        ops: OPS_474,
    },
    BioSpaceGroup {
        number: 182,
        ccp4: 182,
        hm: "P 63 2 2",
        ext: 0,
        qualifier: "",
        hall: "P 6c 2c",
        basisop_idx: 0,
        ops: OPS_475,
    },
    BioSpaceGroup {
        number: 183,
        ccp4: 183,
        hm: "P 6 m m",
        ext: 0,
        qualifier: "",
        hall: "P 6 -2",
        basisop_idx: 0,
        ops: OPS_476,
    },
    BioSpaceGroup {
        number: 184,
        ccp4: 184,
        hm: "P 6 c c",
        ext: 0,
        qualifier: "",
        hall: "P 6 -2c",
        basisop_idx: 0,
        ops: OPS_477,
    },
    BioSpaceGroup {
        number: 185,
        ccp4: 185,
        hm: "P 63 c m",
        ext: 0,
        qualifier: "",
        hall: "P 6c -2",
        basisop_idx: 0,
        ops: OPS_478,
    },
    BioSpaceGroup {
        number: 186,
        ccp4: 186,
        hm: "P 63 m c",
        ext: 0,
        qualifier: "",
        hall: "P 6c -2c",
        basisop_idx: 0,
        ops: OPS_479,
    },
    BioSpaceGroup {
        number: 187,
        ccp4: 187,
        hm: "P -6 m 2",
        ext: 0,
        qualifier: "",
        hall: "P -6 2",
        basisop_idx: 0,
        ops: OPS_480,
    },
    BioSpaceGroup {
        number: 188,
        ccp4: 188,
        hm: "P -6 c 2",
        ext: 0,
        qualifier: "",
        hall: "P -6c 2",
        basisop_idx: 0,
        ops: OPS_481,
    },
    BioSpaceGroup {
        number: 189,
        ccp4: 189,
        hm: "P -6 2 m",
        ext: 0,
        qualifier: "",
        hall: "P -6 -2",
        basisop_idx: 0,
        ops: OPS_482,
    },
    BioSpaceGroup {
        number: 190,
        ccp4: 190,
        hm: "P -6 2 c",
        ext: 0,
        qualifier: "",
        hall: "P -6c -2c",
        basisop_idx: 0,
        ops: OPS_483,
    },
    BioSpaceGroup {
        number: 191,
        ccp4: 191,
        hm: "P 6/m m m",
        ext: 0,
        qualifier: "",
        hall: "-P 6 2",
        basisop_idx: 0,
        ops: OPS_484,
    },
    BioSpaceGroup {
        number: 192,
        ccp4: 192,
        hm: "P 6/m c c",
        ext: 0,
        qualifier: "",
        hall: "-P 6 2c",
        basisop_idx: 0,
        ops: OPS_485,
    },
    BioSpaceGroup {
        number: 193,
        ccp4: 193,
        hm: "P 63/m c m",
        ext: 0,
        qualifier: "",
        hall: "-P 6c 2",
        basisop_idx: 0,
        ops: OPS_486,
    },
    BioSpaceGroup {
        number: 194,
        ccp4: 194,
        hm: "P 63/m m c",
        ext: 0,
        qualifier: "",
        hall: "-P 6c 2c",
        basisop_idx: 0,
        ops: OPS_487,
    },
    BioSpaceGroup {
        number: 195,
        ccp4: 195,
        hm: "P 2 3",
        ext: 0,
        qualifier: "",
        hall: "P 2 2 3",
        basisop_idx: 0,
        ops: OPS_488,
    },
    BioSpaceGroup {
        number: 196,
        ccp4: 196,
        hm: "F 2 3",
        ext: 0,
        qualifier: "",
        hall: "F 2 2 3",
        basisop_idx: 0,
        ops: OPS_489,
    },
    BioSpaceGroup {
        number: 197,
        ccp4: 197,
        hm: "I 2 3",
        ext: 0,
        qualifier: "",
        hall: "I 2 2 3",
        basisop_idx: 0,
        ops: OPS_490,
    },
    BioSpaceGroup {
        number: 198,
        ccp4: 198,
        hm: "P 21 3",
        ext: 0,
        qualifier: "",
        hall: "P 2ac 2ab 3",
        basisop_idx: 0,
        ops: OPS_491,
    },
    BioSpaceGroup {
        number: 199,
        ccp4: 199,
        hm: "I 21 3",
        ext: 0,
        qualifier: "",
        hall: "I 2b 2c 3",
        basisop_idx: 0,
        ops: OPS_492,
    },
    BioSpaceGroup {
        number: 200,
        ccp4: 200,
        hm: "P m -3",
        ext: 0,
        qualifier: "",
        hall: "-P 2 2 3",
        basisop_idx: 0,
        ops: OPS_493,
    },
    BioSpaceGroup {
        number: 201,
        ccp4: 201,
        hm: "P n -3",
        ext: 49,
        qualifier: "",
        hall: "P 2 2 3 -1n",
        basisop_idx: 20,
        ops: OPS_494,
    },
    BioSpaceGroup {
        number: 201,
        ccp4: 0,
        hm: "P n -3",
        ext: 50,
        qualifier: "",
        hall: "-P 2ab 2bc 3",
        basisop_idx: 0,
        ops: OPS_495,
    },
    BioSpaceGroup {
        number: 202,
        ccp4: 202,
        hm: "F m -3",
        ext: 0,
        qualifier: "",
        hall: "-F 2 2 3",
        basisop_idx: 0,
        ops: OPS_496,
    },
    BioSpaceGroup {
        number: 203,
        ccp4: 203,
        hm: "F d -3",
        ext: 49,
        qualifier: "",
        hall: "F 2 2 3 -1d",
        basisop_idx: 27,
        ops: OPS_497,
    },
    BioSpaceGroup {
        number: 203,
        ccp4: 0,
        hm: "F d -3",
        ext: 50,
        qualifier: "",
        hall: "-F 2uv 2vw 3",
        basisop_idx: 0,
        ops: OPS_498,
    },
    BioSpaceGroup {
        number: 204,
        ccp4: 204,
        hm: "I m -3",
        ext: 0,
        qualifier: "",
        hall: "-I 2 2 3",
        basisop_idx: 0,
        ops: OPS_499,
    },
    BioSpaceGroup {
        number: 205,
        ccp4: 205,
        hm: "P a -3",
        ext: 0,
        qualifier: "",
        hall: "-P 2ac 2ab 3",
        basisop_idx: 0,
        ops: OPS_500,
    },
    BioSpaceGroup {
        number: 206,
        ccp4: 206,
        hm: "I a -3",
        ext: 0,
        qualifier: "",
        hall: "-I 2b 2c 3",
        basisop_idx: 0,
        ops: OPS_501,
    },
    BioSpaceGroup {
        number: 207,
        ccp4: 207,
        hm: "P 4 3 2",
        ext: 0,
        qualifier: "",
        hall: "P 4 2 3",
        basisop_idx: 0,
        ops: OPS_502,
    },
    BioSpaceGroup {
        number: 208,
        ccp4: 208,
        hm: "P 42 3 2",
        ext: 0,
        qualifier: "",
        hall: "P 4n 2 3",
        basisop_idx: 0,
        ops: OPS_503,
    },
    BioSpaceGroup {
        number: 209,
        ccp4: 209,
        hm: "F 4 3 2",
        ext: 0,
        qualifier: "",
        hall: "F 4 2 3",
        basisop_idx: 0,
        ops: OPS_504,
    },
    BioSpaceGroup {
        number: 210,
        ccp4: 210,
        hm: "F 41 3 2",
        ext: 0,
        qualifier: "",
        hall: "F 4d 2 3",
        basisop_idx: 0,
        ops: OPS_505,
    },
    BioSpaceGroup {
        number: 211,
        ccp4: 211,
        hm: "I 4 3 2",
        ext: 0,
        qualifier: "",
        hall: "I 4 2 3",
        basisop_idx: 0,
        ops: OPS_506,
    },
    BioSpaceGroup {
        number: 212,
        ccp4: 212,
        hm: "P 43 3 2",
        ext: 0,
        qualifier: "",
        hall: "P 4acd 2ab 3",
        basisop_idx: 0,
        ops: OPS_507,
    },
    BioSpaceGroup {
        number: 213,
        ccp4: 213,
        hm: "P 41 3 2",
        ext: 0,
        qualifier: "",
        hall: "P 4bd 2ab 3",
        basisop_idx: 0,
        ops: OPS_508,
    },
    BioSpaceGroup {
        number: 214,
        ccp4: 214,
        hm: "I 41 3 2",
        ext: 0,
        qualifier: "",
        hall: "I 4bd 2c 3",
        basisop_idx: 0,
        ops: OPS_509,
    },
    BioSpaceGroup {
        number: 215,
        ccp4: 215,
        hm: "P -4 3 m",
        ext: 0,
        qualifier: "",
        hall: "P -4 2 3",
        basisop_idx: 0,
        ops: OPS_510,
    },
    BioSpaceGroup {
        number: 216,
        ccp4: 216,
        hm: "F -4 3 m",
        ext: 0,
        qualifier: "",
        hall: "F -4 2 3",
        basisop_idx: 0,
        ops: OPS_511,
    },
    BioSpaceGroup {
        number: 217,
        ccp4: 217,
        hm: "I -4 3 m",
        ext: 0,
        qualifier: "",
        hall: "I -4 2 3",
        basisop_idx: 0,
        ops: OPS_512,
    },
    BioSpaceGroup {
        number: 218,
        ccp4: 218,
        hm: "P -4 3 n",
        ext: 0,
        qualifier: "",
        hall: "P -4n 2 3",
        basisop_idx: 0,
        ops: OPS_513,
    },
    BioSpaceGroup {
        number: 219,
        ccp4: 219,
        hm: "F -4 3 c",
        ext: 0,
        qualifier: "",
        hall: "F -4a 2 3",
        basisop_idx: 0,
        ops: OPS_514,
    },
    BioSpaceGroup {
        number: 220,
        ccp4: 220,
        hm: "I -4 3 d",
        ext: 0,
        qualifier: "",
        hall: "I -4bd 2c 3",
        basisop_idx: 0,
        ops: OPS_515,
    },
    BioSpaceGroup {
        number: 221,
        ccp4: 221,
        hm: "P m -3 m",
        ext: 0,
        qualifier: "",
        hall: "-P 4 2 3",
        basisop_idx: 0,
        ops: OPS_516,
    },
    BioSpaceGroup {
        number: 222,
        ccp4: 222,
        hm: "P n -3 n",
        ext: 49,
        qualifier: "",
        hall: "P 4 2 3 -1n",
        basisop_idx: 20,
        ops: OPS_517,
    },
    BioSpaceGroup {
        number: 222,
        ccp4: 0,
        hm: "P n -3 n",
        ext: 50,
        qualifier: "",
        hall: "-P 4a 2bc 3",
        basisop_idx: 0,
        ops: OPS_518,
    },
    BioSpaceGroup {
        number: 223,
        ccp4: 223,
        hm: "P m -3 n",
        ext: 0,
        qualifier: "",
        hall: "-P 4n 2 3",
        basisop_idx: 0,
        ops: OPS_519,
    },
    BioSpaceGroup {
        number: 224,
        ccp4: 224,
        hm: "P n -3 m",
        ext: 49,
        qualifier: "",
        hall: "P 4n 2 3 -1n",
        basisop_idx: 30,
        ops: OPS_520,
    },
    BioSpaceGroup {
        number: 224,
        ccp4: 0,
        hm: "P n -3 m",
        ext: 50,
        qualifier: "",
        hall: "-P 4bc 2bc 3",
        basisop_idx: 0,
        ops: OPS_521,
    },
    BioSpaceGroup {
        number: 225,
        ccp4: 225,
        hm: "F m -3 m",
        ext: 0,
        qualifier: "",
        hall: "-F 4 2 3",
        basisop_idx: 0,
        ops: OPS_522,
    },
    BioSpaceGroup {
        number: 226,
        ccp4: 226,
        hm: "F m -3 c",
        ext: 0,
        qualifier: "",
        hall: "-F 4a 2 3",
        basisop_idx: 0,
        ops: OPS_523,
    },
    BioSpaceGroup {
        number: 227,
        ccp4: 227,
        hm: "F d -3 m",
        ext: 49,
        qualifier: "",
        hall: "F 4d 2 3 -1d",
        basisop_idx: 27,
        ops: OPS_524,
    },
    BioSpaceGroup {
        number: 227,
        ccp4: 0,
        hm: "F d -3 m",
        ext: 50,
        qualifier: "",
        hall: "-F 4vw 2vw 3",
        basisop_idx: 0,
        ops: OPS_525,
    },
    BioSpaceGroup {
        number: 228,
        ccp4: 228,
        hm: "F d -3 c",
        ext: 49,
        qualifier: "",
        hall: "F 4d 2 3 -1ad",
        basisop_idx: 37,
        ops: OPS_526,
    },
    BioSpaceGroup {
        number: 228,
        ccp4: 0,
        hm: "F d -3 c",
        ext: 50,
        qualifier: "",
        hall: "-F 4ud 2vw 3",
        basisop_idx: 0,
        ops: OPS_527,
    },
    BioSpaceGroup {
        number: 229,
        ccp4: 229,
        hm: "I m -3 m",
        ext: 0,
        qualifier: "",
        hall: "-I 4 2 3",
        basisop_idx: 0,
        ops: OPS_528,
    },
    BioSpaceGroup {
        number: 230,
        ccp4: 230,
        hm: "I a -3 d",
        ext: 0,
        qualifier: "",
        hall: "-I 4bd 2c 3",
        basisop_idx: 0,
        ops: OPS_529,
    },
    BioSpaceGroup {
        number: 5,
        ccp4: 5005,
        hm: "I 1 21 1",
        ext: 0,
        qualifier: "b4",
        hall: "I 2yb",
        basisop_idx: 38,
        ops: OPS_530,
    },
    BioSpaceGroup {
        number: 5,
        ccp4: 3005,
        hm: "C 1 21 1",
        ext: 0,
        qualifier: "b5",
        hall: "C 2yb",
        basisop_idx: 14,
        ops: OPS_531,
    },
    BioSpaceGroup {
        number: 18,
        ccp4: 1018,
        hm: "P 21212(a)",
        ext: 0,
        qualifier: "",
        hall: "P 2ab 2a",
        basisop_idx: 14,
        ops: OPS_532,
    },
    BioSpaceGroup {
        number: 20,
        ccp4: 1020,
        hm: "C 2 2 21a)",
        ext: 0,
        qualifier: "",
        hall: "C 2ac 2",
        basisop_idx: 39,
        ops: OPS_533,
    },
    BioSpaceGroup {
        number: 21,
        ccp4: 1021,
        hm: "C 2 2 2a",
        ext: 0,
        qualifier: "",
        hall: "C 2ab 2b",
        basisop_idx: 14,
        ops: OPS_534,
    },
    BioSpaceGroup {
        number: 22,
        ccp4: 1022,
        hm: "F 2 2 2a",
        ext: 0,
        qualifier: "",
        hall: "F 2 2c",
        basisop_idx: 40,
        ops: OPS_535,
    },
    BioSpaceGroup {
        number: 23,
        ccp4: 1023,
        hm: "I 2 2 2a",
        ext: 0,
        qualifier: "",
        hall: "I 2ab 2bc",
        basisop_idx: 33,
        ops: OPS_536,
    },
    BioSpaceGroup {
        number: 94,
        ccp4: 1094,
        hm: "P 42 21 2a",
        ext: 0,
        qualifier: "",
        hall: "P 4bc 2a",
        basisop_idx: 20,
        ops: OPS_537,
    },
    BioSpaceGroup {
        number: 197,
        ccp4: 1197,
        hm: "I 2 3a",
        ext: 0,
        qualifier: "",
        hall: "I 2ab 2bc 3",
        basisop_idx: 30,
        ops: OPS_538,
    },
    BioSpaceGroup {
        number: 1,
        ccp4: 0,
        hm: "A 1",
        ext: 0,
        qualifier: "",
        hall: "A 1",
        basisop_idx: 41,
        ops: OPS_539,
    },
    BioSpaceGroup {
        number: 1,
        ccp4: 0,
        hm: "B 1",
        ext: 0,
        qualifier: "",
        hall: "B 1",
        basisop_idx: 42,
        ops: OPS_540,
    },
    BioSpaceGroup {
        number: 1,
        ccp4: 0,
        hm: "C 1",
        ext: 0,
        qualifier: "",
        hall: "C 1",
        basisop_idx: 43,
        ops: OPS_541,
    },
    BioSpaceGroup {
        number: 1,
        ccp4: 0,
        hm: "F 1",
        ext: 0,
        qualifier: "",
        hall: "F 1",
        basisop_idx: 44,
        ops: OPS_542,
    },
    BioSpaceGroup {
        number: 1,
        ccp4: 0,
        hm: "I 1",
        ext: 0,
        qualifier: "",
        hall: "I 1",
        basisop_idx: 45,
        ops: OPS_543,
    },
    BioSpaceGroup {
        number: 2,
        ccp4: 0,
        hm: "A -1",
        ext: 0,
        qualifier: "",
        hall: "-A 1",
        basisop_idx: 41,
        ops: OPS_544,
    },
    BioSpaceGroup {
        number: 2,
        ccp4: 0,
        hm: "B -1",
        ext: 0,
        qualifier: "",
        hall: "-B 1",
        basisop_idx: 42,
        ops: OPS_545,
    },
    BioSpaceGroup {
        number: 2,
        ccp4: 0,
        hm: "C -1",
        ext: 0,
        qualifier: "",
        hall: "-C 1",
        basisop_idx: 43,
        ops: OPS_546,
    },
    BioSpaceGroup {
        number: 2,
        ccp4: 0,
        hm: "F -1",
        ext: 0,
        qualifier: "",
        hall: "-F 1",
        basisop_idx: 44,
        ops: OPS_547,
    },
    BioSpaceGroup {
        number: 2,
        ccp4: 0,
        hm: "I -1",
        ext: 0,
        qualifier: "",
        hall: "-I 1",
        basisop_idx: 45,
        ops: OPS_548,
    },
    BioSpaceGroup {
        number: 3,
        ccp4: 0,
        hm: "B 1 2 1",
        ext: 0,
        qualifier: "b1",
        hall: "B 2y",
        basisop_idx: 46,
        ops: OPS_549,
    },
    BioSpaceGroup {
        number: 3,
        ccp4: 0,
        hm: "C 1 1 2",
        ext: 0,
        qualifier: "c1",
        hall: "C 2",
        basisop_idx: 47,
        ops: OPS_550,
    },
    BioSpaceGroup {
        number: 4,
        ccp4: 0,
        hm: "B 1 21 1",
        ext: 0,
        qualifier: "b1",
        hall: "B 2yb",
        basisop_idx: 46,
        ops: OPS_551,
    },
    BioSpaceGroup {
        number: 4,
        ccp4: 0,
        hm: "C 1 1 21",
        ext: 0,
        qualifier: "c2",
        hall: "C 2c",
        basisop_idx: 47,
        ops: OPS_552,
    },
    BioSpaceGroup {
        number: 5,
        ccp4: 0,
        hm: "F 1 2 1",
        ext: 0,
        qualifier: "b6",
        hall: "F 2y",
        basisop_idx: 48,
        ops: OPS_553,
    },
    BioSpaceGroup {
        number: 8,
        ccp4: 0,
        hm: "F 1 m 1",
        ext: 0,
        qualifier: "b4",
        hall: "F -2y",
        basisop_idx: 48,
        ops: OPS_554,
    },
    BioSpaceGroup {
        number: 9,
        ccp4: 0,
        hm: "F 1 d 1",
        ext: 0,
        qualifier: "b4",
        hall: "F -2yuw",
        basisop_idx: 49,
        ops: OPS_555,
    },
    BioSpaceGroup {
        number: 12,
        ccp4: 0,
        hm: "F 1 2/m 1",
        ext: 0,
        qualifier: "b4",
        hall: "-F 2y",
        basisop_idx: 48,
        ops: OPS_556,
    },
    BioSpaceGroup {
        number: 64,
        ccp4: 0,
        hm: "A b a m",
        ext: 0,
        qualifier: "",
        hall: "-A 2 2ab",
        basisop_idx: 3,
        ops: OPS_557,
    },
    BioSpaceGroup {
        number: 89,
        ccp4: 0,
        hm: "C 4 2 2",
        ext: 0,
        qualifier: "",
        hall: "C 4 2",
        basisop_idx: 50,
        ops: OPS_558,
    },
    BioSpaceGroup {
        number: 90,
        ccp4: 0,
        hm: "C 4 2 21",
        ext: 0,
        qualifier: "",
        hall: "C 4a 2",
        basisop_idx: 50,
        ops: OPS_559,
    },
    BioSpaceGroup {
        number: 97,
        ccp4: 0,
        hm: "F 4 2 2",
        ext: 0,
        qualifier: "",
        hall: "F 4 2",
        basisop_idx: 50,
        ops: OPS_560,
    },
    BioSpaceGroup {
        number: 115,
        ccp4: 0,
        hm: "C -4 2 m",
        ext: 0,
        qualifier: "",
        hall: "C -4 2",
        basisop_idx: 50,
        ops: OPS_561,
    },
    BioSpaceGroup {
        number: 117,
        ccp4: 0,
        hm: "C -4 2 b",
        ext: 0,
        qualifier: "",
        hall: "C -4 2ya",
        basisop_idx: 50,
        ops: OPS_562,
    },
    BioSpaceGroup {
        number: 139,
        ccp4: 0,
        hm: "F 4/m m m",
        ext: 0,
        qualifier: "",
        hall: "-F 4 2",
        basisop_idx: 50,
        ops: OPS_563,
    },
];

// Gemmi✔️✔️: const SpaceGroupAltName spacegroup_tables::alt_names[28] = {
static GEMMI_ALT_NAMES: &[BioAltName] = &[
    BioAltName {
        hm: "A e m 2",
        ext: 0,
        pos: 190,
    },
    BioAltName {
        hm: "B m e 2",
        ext: 0,
        pos: 191,
    },
    BioAltName {
        hm: "B 2 e m",
        ext: 0,
        pos: 192,
    },
    BioAltName {
        hm: "C 2 m e",
        ext: 0,
        pos: 193,
    },
    BioAltName {
        hm: "C m 2 e",
        ext: 0,
        pos: 194,
    },
    BioAltName {
        hm: "A e 2 m",
        ext: 0,
        pos: 195,
    },
    BioAltName {
        hm: "A e a 2",
        ext: 0,
        pos: 202,
    },
    BioAltName {
        hm: "B b e 2",
        ext: 0,
        pos: 203,
    },
    BioAltName {
        hm: "B 2 e b",
        ext: 0,
        pos: 204,
    },
    BioAltName {
        hm: "C 2 c e",
        ext: 0,
        pos: 205,
    },
    BioAltName {
        hm: "C c 2 e",
        ext: 0,
        pos: 206,
    },
    BioAltName {
        hm: "A e 2 a",
        ext: 0,
        pos: 207,
    },
    BioAltName {
        hm: "C m c e",
        ext: 0,
        pos: 303,
    },
    BioAltName {
        hm: "C c m e",
        ext: 0,
        pos: 304,
    },
    BioAltName {
        hm: "A e m a",
        ext: 0,
        pos: 305,
    },
    BioAltName {
        hm: "A e a m",
        ext: 0,
        pos: 306,
    },
    BioAltName {
        hm: "B b e m",
        ext: 0,
        pos: 307,
    },
    BioAltName {
        hm: "B m e b",
        ext: 0,
        pos: 308,
    },
    BioAltName {
        hm: "C m m e",
        ext: 0,
        pos: 315,
    },
    BioAltName {
        hm: "A e m m",
        ext: 0,
        pos: 317,
    },
    BioAltName {
        hm: "B m e m",
        ext: 0,
        pos: 319,
    },
    BioAltName {
        hm: "C c c e",
        ext: 49,
        pos: 321,
    },
    BioAltName {
        hm: "C c c e",
        ext: 50,
        pos: 322,
    },
    BioAltName {
        hm: "A e a a",
        ext: 49,
        pos: 325,
    },
    BioAltName {
        hm: "A e a a",
        ext: 50,
        pos: 326,
    },
    BioAltName {
        hm: "B b e b",
        ext: 49,
        pos: 329,
    },
    BioAltName {
        hm: "B b e b",
        ext: 50,
        pos: 330,
    },
    BioAltName {
        hm: "P 21 21 2a",
        ext: 0,
        pos: 532,
    },
];

#[cfg(test)]
mod tests {
    use super::super::{
        BioCrystalCell, BioCrystalInfo, BioNcsOperator, BioTransform, find_nearest_image,
    };
    use super::{
        BioSpaceGroup, BioSymOp, GEMMI_ALT_NAMES, GEMMI_OP_DEN, GEMMI_SPACEGROUPS,
        NormalizedSpaceGroupName, ParsedSpaceGroupName, add_ncs_images_to_crystal,
        find_main_spacegroup_by_name, find_spacegroup_by_name, find_spacegroup_by_number,
        find_structure_spacegroup, parse_spacegroup_name, set_crystal_symmetry_images,
        setup_cell_images,
    };
    use crate::BioAsu;

    fn find_numeric_spacegroup_name(name: &str) -> Option<&'static BioSpaceGroup> {
        match parse_spacegroup_name(name)? {
            ParsedSpaceGroupName::Numeric { ccp4 } => find_spacegroup_by_number(ccp4),
            ParsedSpaceGroupName::Named(_) => None,
        }
    }

    fn find_main_named_spacegroup(
        name: &str,
        alpha: f64,
        gamma: f64,
    ) -> Option<&'static BioSpaceGroup> {
        match parse_spacegroup_name(name)? {
            ParsedSpaceGroupName::Named(normalized) => {
                find_main_spacegroup_by_name(&normalized, alpha, gamma, false, false)
            }
            ParsedSpaceGroupName::Numeric { .. } => None,
        }
    }

    fn crystal_with_group_and_images(
        group_name: &str,
        cs_count: i16,
        images: Vec<BioTransform>,
    ) -> BioCrystalInfo {
        BioCrystalInfo::new(
            BioCrystalCell {
                a: 10.0,
                b: 10.0,
                c: 10.0,
                ..BioCrystalCell::default()
            },
            Some(group_name.to_owned()),
            None,
            BioTransform::new(
                [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
                [0.0; 3],
            ),
            BioTransform::new(
                [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]],
                [0.0; 3],
            ),
            false,
            cs_count,
            images,
        )
    }

    #[test]
    fn lattice_spacegroup_numeric_names_follow_ccp4_and_strtol_boundaries() {
        let group_1003 = find_numeric_spacegroup_name("1003").unwrap();
        assert_eq!(
            (group_1003.ccp4, group_1003.number, group_1003.hm),
            (1003, 3, "P 1 1 2")
        );
        assert_eq!(
            GEMMI_SPACEGROUPS
                .iter()
                .position(|group| group.ccp4 == 1003),
            Some(3)
        );

        let group_1004 = find_numeric_spacegroup_name("1004").unwrap();
        assert_eq!(
            (group_1004.ccp4, group_1004.number, group_1004.hm),
            (1004, 4, "P 1 1 21")
        );
        let group_19 = find_numeric_spacegroup_name("19").unwrap();
        assert_eq!(
            (group_19.ccp4, group_19.number, group_19.hm),
            (19, 19, "P 21 21 21")
        );

        // Gemmi maps numeric zero to main[0], not to a table row whose ccp4
        // field happens to be zero.
        let group_zero = find_numeric_spacegroup_name("0").unwrap();
        assert_eq!(
            (group_zero.ccp4, group_zero.number, group_zero.hm),
            (1, 1, "P 1")
        );
        assert!(GEMMI_SPACEGROUPS.iter().any(|group| group.ccp4 == 0));

        assert_eq!(
            find_numeric_spacegroup_name(" \t_1003").unwrap().hm,
            "P 1 1 2"
        );
        assert!(find_numeric_spacegroup_name("999999").is_none());
        assert!(find_numeric_spacegroup_name("19x").is_none());
        assert!(find_numeric_spacegroup_name("19 ").is_none());
        assert!(find_numeric_spacegroup_name("+19").is_none());
        assert!(find_numeric_spacegroup_name("-19").is_none());
        assert!(find_numeric_spacegroup_name("2147483648").is_none());
        assert!(
            find_numeric_spacegroup_name("999999999999999999999999999999999999999999999999")
                .is_none()
        );
        assert_eq!(
            find_numeric_spacegroup_name("19\0ignored").map(|group| (group.ccp4, group.hm)),
            Some((19, "P 21 21 21"))
        );

        // The pinned LP64 Gemmi/GCC reference narrows long to int modulo 2^32;
        // ILP32 strtol saturates at LONG_MAX before the same int parameter.
        let wrapped_long = if std::mem::size_of::<std::ffi::c_long>() > std::mem::size_of::<i32>() {
            Some((1003, 3, "P 1 1 2"))
        } else {
            None
        };
        assert_eq!(
            find_numeric_spacegroup_name("4294968299").map(|group| (
                group.ccp4,
                group.number,
                group.hm
            )),
            wrapped_long
        );
    }

    #[test]
    fn lattice_spacegroup_name_normalization_keeps_source_byte_rules() {
        let parse_named = |name: &str| match parse_spacegroup_name(name) {
            Some(ParsedSpaceGroupName::Named(normalized)) => normalized,
            other => panic!("expected normalized name, got {other:?}"),
        };

        let lowercase_initial = parse_named("p 21 21 21");
        assert_eq!(lowercase_initial.first, b'P');
        assert_eq!(lowercase_initial.bytes, b"p 21 21 21");
        assert_eq!(lowercase_initial.start, 2);
        assert_eq!(lowercase_initial.c_string_end, b"p 21 21 21".len());

        let source_whitespace = "\t_ P\t21_21 21";
        let spaced = parse_named(source_whitespace);
        assert_eq!(spaced.first, b'P');
        assert_eq!(spaced.bytes, source_whitespace.as_bytes());
        assert_eq!(spaced.start, 5);
        assert_eq!(spaced.c_string_end, source_whitespace.len());

        assert_eq!(
            parse_named("r 3:h"),
            NormalizedSpaceGroupName {
                bytes: b"r 3:H".to_vec(),
                first: b'R',
                start: 2,
                c_string_end: 5,
            }
        );
        assert_eq!(
            parse_named("R3R"),
            NormalizedSpaceGroupName {
                bytes: b"R3:R".to_vec(),
                first: b'R',
                start: 1,
                c_string_end: 4,
            }
        );
        assert_eq!(
            parse_named("R3r ").bytes,
            b"R3r ",
            "the source suffix rule checks the actual final byte"
        );
    }

    #[test]
    fn lattice_spacegroup_main_name_lookup_matches_pinned_gemmi() {
        let summary = |name: &str, alpha: f64, gamma: f64| {
            find_main_named_spacegroup(name, alpha, gamma)
                .map(|group| (group.number, group.ccp4, group.ext, group.hm))
        };

        // Pinned Gemmi 5cc1c23: when alpha is exactly zero, the source's
        // ternary selects the null-preference flag rather than comparing
        // gamma against alpha, so the default H setting remains first.
        assert_eq!(summary("R 3", 0.0, 120.0), Some((146, 146, b'H', "R 3")));
        assert_eq!(
            summary("R 3", 90.0, 101.249),
            Some((146, 1146, b'R', "R 3"))
        );
        assert_eq!(
            summary("R 3", 90.0, 101.25),
            Some((146, 146, b'H', "R 3")),
            "the source comparison is strict at 1.125 * alpha"
        );
        assert_eq!(summary("R 3", 90.0, 101.251), Some((146, 146, b'H', "R 3")));

        assert_eq!(summary("R 3:H", 90.0, 100.0), Some((146, 146, b'H', "R 3")));
        assert_eq!(
            summary("R 3:R", 90.0, 120.0),
            Some((146, 1146, b'R', "R 3"))
        );
        assert_eq!(summary("R3R", 90.0, 120.0), Some((146, 1146, b'R', "R 3")));

        for cubic_name in ["P m -3 m", "P m 3 m"] {
            assert_eq!(
                summary(cubic_name, 90.0, 90.0),
                Some((221, 221, 0, "P m -3 m")),
                "both the current and legacy cubic spelling resolve to the first matching main row"
            );
        }
        for cubic_name in ["F d -3 m", "F d 3 m"] {
            assert_eq!(
                summary(cubic_name, 90.0, 90.0),
                Some((227, 227, b'1', "F d -3 m"))
            );
        }
        assert_eq!(
            summary("P 21 21 21", 90.0, 90.0),
            Some((19, 19, 0, "P 21 21 21"))
        );
        assert!(summary("Q 1", 90.0, 90.0).is_none());
        assert!(summary("R 3:Q", 90.0, 90.0).is_none());
    }

    #[test]
    fn lattice_spacegroup_monoclinic_short_names_preserve_source_order_and_b_exception() {
        let summary = |name: &str| {
            find_spacegroup_by_name(name, 90.0, 90.0)
                .map(|group| (group.number, group.ccp4, group.ext, group.hm))
        };

        // Pinned Gemmi's short-name loop returns the first main-table row
        // whose remaining symbol matches; the B exception mutates its source
        // cursor twice and makes "B 2" select B 1 1 2 rather than B 1 2 1.
        assert_eq!(summary("P2"), Some((3, 3, 0, "P 1 2 1")));
        assert_eq!(summary("P 2"), Some((3, 3, 0, "P 1 2 1")));
        assert_eq!(summary("P 2 1 1"), Some((3, 0, 0, "P 2 1 1")));
        assert_eq!(summary("P 1 2 1"), Some((3, 3, 0, "P 1 2 1")));
        assert_eq!(summary("B2"), Some((5, 1005, 0, "B 1 1 2")));
        assert_eq!(summary("B 1 1 2"), Some((5, 1005, 0, "B 1 1 2")));
        assert_eq!(summary("B 1 2 1"), Some((3, 0, 0, "B 1 2 1")));
    }

    #[test]
    fn lattice_spacegroup_alternate_names_preserve_source_extensions_and_tail_rules() {
        let summary = |name: &str| {
            find_spacegroup_by_name(name, 90.0, 90.0)
                .map(|group| (group.number, group.ccp4, group.ext, group.hm))
        };

        assert_eq!(summary("A e m 2"), Some((39, 39, 0, "A b m 2")));
        assert_eq!(summary("AEM2"), Some((39, 39, 0, "A b m 2")));
        assert_eq!(summary("A e m 2 "), Some((39, 39, 0, "A b m 2")));
        assert_eq!(summary("A e m 2x"), None);
        assert_eq!(summary("C c c e"), Some((68, 68, b'1', "C c c a")));
        assert_eq!(summary("C c c e:1"), Some((68, 68, b'1', "C c c a")));
        assert_eq!(summary("C c c e:2"), Some((68, 0, b'2', "C c c a")));
        assert_eq!(summary("C c c e:3"), None);
        // The pinned source checks the first extension byte but does not
        // require the input to end immediately after that byte.
        assert_eq!(summary("C c c e:1junk"), Some((68, 68, b'1', "C c c a")));
        assert_eq!(summary("C c c e:2junk"), Some((68, 0, b'2', "C c c a")));
        assert_eq!(summary("P 21 21 2a"), Some((18, 1018, 0, "P 21212(a)")));
        assert_eq!(summary("Q 2"), None);
    }

    #[test]
    fn lattice_spacegroup_structure_dispatch_suppresses_noncrystal_lookup() {
        let cell = BioCrystalCell::default();
        let noncrystal = BioCrystalInfo::new(
            cell,
            Some("P2".to_owned()),
            None,
            BioTransform::identity(),
            BioTransform::identity(),
            false,
            0,
            Vec::new(),
        );
        assert!(!noncrystal.is_crystal());
        assert!(find_spacegroup_by_name("P2", cell.alpha, cell.gamma).is_some());
        assert!(find_structure_spacegroup(&noncrystal).is_none());

        let crystal_cell = BioCrystalCell {
            a: 10.0,
            b: 10.0,
            c: 10.0,
            ..BioCrystalCell::default()
        };
        let crystal = BioCrystalInfo::new(
            crystal_cell,
            Some("P2".to_owned()),
            None,
            BioTransform::new(
                [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
                [0.0; 3],
            ),
            BioTransform::new(
                [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]],
                [0.0; 3],
            ),
            false,
            0,
            Vec::new(),
        );
        assert!(crystal.is_crystal());
        assert_eq!(
            find_structure_spacegroup(&crystal).map(|group| (group.ccp4, group.hm)),
            Some((3, "P 1 2 1"))
        );
    }

    #[test]
    fn lattice_crystal_images_identity_only_clears_previous_state() {
        let previous = BioTransform::new(
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [9.0, 8.0, 7.0],
        );
        let mut crystal = crystal_with_group_and_images("P 1", 17, vec![previous]);

        set_crystal_symmetry_images(&mut crystal);

        assert_eq!(crystal.cs_count(), 0);
        assert!(crystal.symmetry_images().is_empty());
    }

    #[test]
    fn lattice_crystal_images_centered_group_preserves_order_and_rebuilds() {
        let previous = BioTransform::new(
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [9.0, 8.0, 7.0],
        );
        let mut crystal = crystal_with_group_and_images("C 2 2 2", 17, vec![previous]);
        let expected = [
            BioTransform::new(
                [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]],
                [0.0; 3],
            ),
            BioTransform::new(
                [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                [0.0; 3],
            ),
            BioTransform::new(
                [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                [0.0; 3],
            ),
            BioTransform::new(
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [0.5, 0.5, 0.0],
            ),
            BioTransform::new(
                [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]],
                [0.5, 0.5, 0.0],
            ),
            BioTransform::new(
                [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                [0.5, 0.5, 0.0],
            ),
            BioTransform::new(
                [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                [0.5, 0.5, 0.0],
            ),
        ];

        assert_eq!(
            find_structure_spacegroup(&crystal).map(|group| (group.number, group.hm)),
            Some((21, "C 2 2 2"))
        );
        set_crystal_symmetry_images(&mut crystal);
        assert_eq!(crystal.cs_count(), 7);
        assert_eq!(crystal.symmetry_images(), expected);

        set_crystal_symmetry_images(&mut crystal);
        assert_eq!(crystal.cs_count(), 7);
        assert_eq!(crystal.symmetry_images(), expected);
    }

    #[test]
    fn lattice_crystal_images_unknown_group_clears_previous_state() {
        let previous = BioTransform::new(
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [9.0, 8.0, 7.0],
        );
        let mut crystal = crystal_with_group_and_images("Q 9", 17, vec![previous]);

        assert!(crystal.is_crystal());
        assert!(find_structure_spacegroup(&crystal).is_none());
        set_crystal_symmetry_images(&mut crystal);

        assert_eq!(crystal.cs_count(), 0);
        assert!(crystal.symmetry_images().is_empty());
    }

    #[test]
    fn lattice_ncs_images_expand_fixed_order_from_original_crystal_images() {
        let mut crystal = crystal_with_group_and_images("C 2 2 2", 0, Vec::new());
        set_crystal_symmetry_images(&mut crystal);
        assert_eq!(crystal.cs_count(), 7);

        crystal.orthogonal = BioTransform::new(
            [[2.0, 1.0, 0.0], [0.0, 3.0, 1.0], [0.0, 0.0, 4.0]],
            [5.0, 6.0, 7.0],
        );
        crystal.fractional = BioTransform::new(
            [[0.5, 0.0, 0.0], [0.0, 0.25, 0.0], [0.0, 0.0, 0.125]],
            [-1.0, -2.0, -3.0],
        );

        let ncs = [
            BioNcsOperator::new(
                "given-before".to_owned(),
                true,
                BioTransform::new(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                    [101.0, 102.0, 103.0],
                ),
            ),
            BioNcsOperator::new(
                "first-expanded".to_owned(),
                false,
                BioTransform::new(
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                    [1.0, 2.0, 3.0],
                ),
            ),
            BioNcsOperator::new(
                "given-between".to_owned(),
                true,
                BioTransform::new(
                    [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]],
                    [-101.0, -102.0, -103.0],
                ),
            ),
            BioNcsOperator::new(
                "second-expanded".to_owned(),
                false,
                BioTransform::new(
                    [[1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                    [-2.0, 1.0, 0.0],
                ),
            ),
        ];

        let first_ncs = BioTransform::new(
            [[1.0, 0.5, 0.0], [0.0, 0.75, 0.25], [0.0, 0.0, 0.5]],
            [2.0, 0.0, -1.75],
        );
        let first_combinations = [
            BioTransform::new(
                [[-1.0, -0.5, 0.0], [0.0, -0.75, -0.25], [0.0, 0.0, 0.5]],
                [-2.0, 0.0, -1.75],
            ),
            BioTransform::new(
                [[1.0, 0.5, 0.0], [0.0, -0.75, -0.25], [0.0, 0.0, -0.5]],
                [2.0, 0.0, 1.75],
            ),
            BioTransform::new(
                [[-1.0, -0.5, 0.0], [0.0, 0.75, 0.25], [0.0, 0.0, -0.5]],
                [-2.0, 0.0, 1.75],
            ),
            BioTransform::new(
                [[1.0, 0.5, 0.0], [0.0, 0.75, 0.25], [0.0, 0.0, 0.5]],
                [2.5, 0.5, -1.75],
            ),
            BioTransform::new(
                [[-1.0, -0.5, 0.0], [0.0, -0.75, -0.25], [0.0, 0.0, 0.5]],
                [-1.5, 0.5, -1.75],
            ),
            BioTransform::new(
                [[1.0, 0.5, 0.0], [0.0, -0.75, -0.25], [0.0, 0.0, -0.5]],
                [2.5, 0.5, 1.75],
            ),
            BioTransform::new(
                [[-1.0, -0.5, 0.0], [0.0, 0.75, 0.25], [0.0, 0.0, -0.5]],
                [-1.5, 0.5, 1.75],
            ),
        ];
        let second_ncs = BioTransform::new(
            [[1.0, 2.0, 0.5], [0.0, 0.75, 0.25], [0.0, 0.0, 0.5]],
            [3.5, -0.25, -2.125],
        );
        let second_combinations = [
            BioTransform::new(
                [[-1.0, -2.0, -0.5], [0.0, -0.75, -0.25], [0.0, 0.0, 0.5]],
                [-3.5, 0.25, -2.125],
            ),
            BioTransform::new(
                [[1.0, 2.0, 0.5], [0.0, -0.75, -0.25], [0.0, 0.0, -0.5]],
                [3.5, 0.25, 2.125],
            ),
            BioTransform::new(
                [[-1.0, -2.0, -0.5], [0.0, 0.75, 0.25], [0.0, 0.0, -0.5]],
                [-3.5, -0.25, 2.125],
            ),
            BioTransform::new(
                [[1.0, 2.0, 0.5], [0.0, 0.75, 0.25], [0.0, 0.0, 0.5]],
                [4.0, 0.25, -2.125],
            ),
            BioTransform::new(
                [[-1.0, -2.0, -0.5], [0.0, -0.75, -0.25], [0.0, 0.0, 0.5]],
                [-3.0, 0.75, -2.125],
            ),
            BioTransform::new(
                [[1.0, 2.0, 0.5], [0.0, -0.75, -0.25], [0.0, 0.0, -0.5]],
                [4.0, 0.75, 2.125],
            ),
            BioTransform::new(
                [[-1.0, -2.0, -0.5], [0.0, 0.75, 0.25], [0.0, 0.0, -0.5]],
                [-3.0, 0.25, 2.125],
            ),
        ];
        let expected = [
            first_ncs,
            first_combinations[0],
            first_combinations[1],
            first_combinations[2],
            first_combinations[3],
            first_combinations[4],
            first_combinations[5],
            first_combinations[6],
            second_ncs,
            second_combinations[0],
            second_combinations[1],
            second_combinations[2],
            second_combinations[3],
            second_combinations[4],
            second_combinations[5],
            second_combinations[6],
        ];
        let original_crystal_images = crystal.symmetry_images().to_vec();

        add_ncs_images_to_crystal(&mut crystal, &ncs);

        assert_eq!(crystal.cs_count(), 7);
        assert_eq!(crystal.symmetry_images().len(), 23);
        assert_eq!(
            &crystal.symmetry_images()[..7],
            original_crystal_images.as_slice()
        );
        assert_eq!(
            &crystal.symmetry_images()[7..],
            expected,
            "NCS rows stay in source order; each combines only with the original seven crystal images"
        );
    }

    #[test]
    fn lattice_setup_cell_images_rebuilds_group_then_appends_ncs_once() {
        let stale = BioTransform::new(
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [9.0, 8.0, 7.0],
        );
        let mut crystal = crystal_with_group_and_images("C 2 2 2", 19, vec![stale]);
        let given = BioNcsOperator::new(
            "given".to_owned(),
            true,
            BioTransform::new(
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [9.0, 9.0, 9.0],
            ),
        );
        let generated = BioNcsOperator::new(
            "generated".to_owned(),
            false,
            BioTransform::new(
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [2.5, 0.0, 0.0],
            ),
        );
        let ncs = [given, generated];

        setup_cell_images(&mut crystal, &ncs);

        let expected = [
            BioTransform::new(
                [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]],
                [0.0; 3],
            ),
            BioTransform::new(
                [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                [0.0; 3],
            ),
            BioTransform::new(
                [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                [0.0; 3],
            ),
            BioTransform::new(
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [0.5, 0.5, 0.0],
            ),
            BioTransform::new(
                [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]],
                [0.5, 0.5, 0.0],
            ),
            BioTransform::new(
                [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                [0.5, 0.5, 0.0],
            ),
            BioTransform::new(
                [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                [0.5, 0.5, 0.0],
            ),
            BioTransform::new(
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [0.25, 0.0, 0.0],
            ),
            BioTransform::new(
                [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]],
                [-0.25, 0.0, 0.0],
            ),
            BioTransform::new(
                [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                [0.25, 0.0, 0.0],
            ),
            BioTransform::new(
                [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                [-0.25, 0.0, 0.0],
            ),
            BioTransform::new(
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [0.75, 0.5, 0.0],
            ),
            BioTransform::new(
                [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]],
                [0.25, 0.5, 0.0],
            ),
            BioTransform::new(
                [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
                [0.75, 0.5, 0.0],
            ),
            BioTransform::new(
                [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
                [0.25, 0.5, 0.0],
            ),
        ];
        assert_eq!(crystal.cs_count(), 7);
        assert_eq!(crystal.symmetry_images(), expected);

        let once = crystal.clone();
        setup_cell_images(&mut crystal, &ncs);
        assert_eq!(
            crystal, once,
            "setup replaces rather than accumulates images"
        );
    }

    #[test]
    fn lattice_setup_cell_images_clears_unknown_and_noncrystal_state() {
        let stale = BioTransform::new(
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            [9.0, 8.0, 7.0],
        );
        let mut unknown = crystal_with_group_and_images("Q 9", 1, vec![stale]);
        setup_cell_images(&mut unknown, &[]);
        assert_eq!(unknown.cs_count(), 0);
        assert!(unknown.symmetry_images().is_empty());

        let mut noncrystal = BioCrystalInfo::new(
            BioCrystalCell::default(),
            Some("P 1".to_owned()),
            None,
            BioTransform::identity(),
            BioTransform::identity(),
            false,
            1,
            vec![stale],
        );
        assert!(!noncrystal.is_crystal());
        setup_cell_images(&mut noncrystal, &[]);
        assert_eq!(noncrystal.cs_count(), 0);
        assert!(noncrystal.symmetry_images().is_empty());
    }

    #[test]
    fn lattice_ncs_setup_changes_nearest_image_selection() {
        let mut crystal = BioCrystalInfo::new(
            BioCrystalCell {
                a: 8.0,
                b: 8.0,
                c: 8.0,
                ..BioCrystalCell::default()
            },
            Some("P 1".to_owned()),
            None,
            BioTransform::new(
                [[8.0, 0.0, 0.0], [0.0, 8.0, 0.0], [0.0, 0.0, 8.0]],
                [0.0; 3],
            ),
            BioTransform::new(
                [[0.125, 0.0, 0.0], [0.0, 0.125, 0.0], [0.0, 0.0, 0.125]],
                [0.0; 3],
            ),
            false,
            3,
            Vec::new(),
        );
        setup_cell_images(&mut crystal, &[]);
        let without_ncs =
            find_nearest_image(&crystal, [0.0, 0.0, 0.0], [5.0, 0.0, 0.0], BioAsu::Any);
        assert_eq!(without_ncs.dist_sq().to_bits(), 9.0_f64.to_bits());
        assert_eq!(without_ncs.pbc_shift(), &[-1, 0, 0]);
        assert_eq!(without_ncs.sym_idx(), 0);

        let ncs = [BioNcsOperator::new(
            "translated".to_owned(),
            false,
            BioTransform::new(
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [2.0, 0.0, 0.0],
            ),
        )];
        setup_cell_images(&mut crystal, &ncs);
        let with_ncs = find_nearest_image(&crystal, [0.0, 0.0, 0.0], [5.0, 0.0, 0.0], BioAsu::Any);
        assert_eq!(crystal.cs_count(), 0);
        assert_eq!(crystal.symmetry_images().len(), 1);
        assert_eq!(with_ncs.dist_sq().to_bits(), 1.0_f64.to_bits());
        assert_eq!(with_ncs.pbc_shift(), &[-1, 0, 0]);
        assert_eq!(with_ncs.sym_idx(), 1);
    }

    #[test]
    fn lattice_spacegroup_table_retains_all_rows_and_ccp4_mapping_data() {
        assert_eq!(GEMMI_SPACEGROUPS.len(), 564);
        assert_eq!(GEMMI_ALT_NAMES.len(), 28);
        assert_eq!(
            GEMMI_SPACEGROUPS
                .iter()
                .map(|group| group.ops.len())
                .sum::<usize>(),
            7_740
        );

        let ccp4_1003 = GEMMI_SPACEGROUPS
            .iter()
            .position(|group| group.ccp4 == 1003);
        assert_eq!(ccp4_1003, Some(3));
        let ccp4_1003_group = &GEMMI_SPACEGROUPS[3];
        assert_eq!(ccp4_1003_group.number, 3);
        assert_eq!(ccp4_1003_group.hm, "P 1 1 2");
        assert_eq!(ccp4_1003_group.qualifier, "c");
        assert_eq!(ccp4_1003_group.hall, "P 2");
        assert_eq!(ccp4_1003_group.basisop_idx, 1);

        // This row has ccp4 == 0, while Gemmi's numeric lookup treats input 0
        // specially as main[0]. The data must not be conflated with `number`.
        assert_eq!(GEMMI_SPACEGROUPS[0].ccp4, 1);
        assert_eq!(GEMMI_SPACEGROUPS[4].number, 3);
        assert_eq!(GEMMI_SPACEGROUPS[4].ccp4, 0);
        assert_eq!(GEMMI_SPACEGROUPS[4].hm, "P 2 1 1");
    }

    #[test]
    fn lattice_spacegroup_alternate_names_preserve_targets_and_extensions() {
        let actual = GEMMI_ALT_NAMES
            .iter()
            .map(|entry| (entry.hm, entry.ext, entry.pos))
            .collect::<Vec<_>>();
        let expected: [(&str, u8, usize); 28] = [
            ("A e m 2", 0, 190),
            ("B m e 2", 0, 191),
            ("B 2 e m", 0, 192),
            ("C 2 m e", 0, 193),
            ("C m 2 e", 0, 194),
            ("A e 2 m", 0, 195),
            ("A e a 2", 0, 202),
            ("B b e 2", 0, 203),
            ("B 2 e b", 0, 204),
            ("C 2 c e", 0, 205),
            ("C c 2 e", 0, 206),
            ("A e 2 a", 0, 207),
            ("C m c e", 0, 303),
            ("C c m e", 0, 304),
            ("A e m a", 0, 305),
            ("A e a m", 0, 306),
            ("B b e m", 0, 307),
            ("B m e b", 0, 308),
            ("C m m e", 0, 315),
            ("A e m m", 0, 317),
            ("B m e m", 0, 319),
            ("C c c e", b'1', 321),
            ("C c c e", b'2', 322),
            ("A e a a", b'1', 325),
            ("A e a a", b'2', 326),
            ("B b e b", b'1', 329),
            ("B b e b", b'2', 330),
            ("P 21 21 2a", 0, 532),
        ];
        assert_eq!(actual, expected);
        assert!(
            GEMMI_ALT_NAMES
                .iter()
                .all(|entry| entry.pos < GEMMI_SPACEGROUPS.len())
        );

        assert_eq!(GEMMI_SPACEGROUPS[190].hm, "A b m 2");
        assert_eq!(GEMMI_SPACEGROUPS[321].hm, "C c c a");
        assert_eq!(GEMMI_SPACEGROUPS[321].ext, b'1');
        assert_eq!(GEMMI_SPACEGROUPS[322].hm, "C c c a");
        assert_eq!(GEMMI_SPACEGROUPS[322].ext, b'2');
    }

    #[test]
    fn lattice_spacegroup_operations_retain_denominator_order_and_settings() {
        assert_eq!(GEMMI_OP_DEN, 24);

        let p212121: &BioSpaceGroup = &GEMMI_SPACEGROUPS[114];
        assert_eq!(p212121.number, 19);
        assert_eq!(p212121.ccp4, 19);
        assert_eq!(p212121.hm, "P 21 21 21");
        let expected_p212121 = [
            BioSymOp {
                rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
                tran: [0, 0, 0],
            },
            BioSymOp {
                rot: [[-24, 0, 0], [0, -24, 0], [0, 0, 24]],
                tran: [12, 0, 12],
            },
            BioSymOp {
                rot: [[24, 0, 0], [0, -24, 0], [0, 0, -24]],
                tran: [12, 12, 0],
            },
            BioSymOp {
                rot: [[-24, 0, 0], [0, 24, 0], [0, 0, -24]],
                tran: [0, 12, 12],
            },
        ];
        assert_eq!(p212121.ops, expected_p212121);

        let r3_hexagonal = &GEMMI_SPACEGROUPS[432];
        assert_eq!(r3_hexagonal.ccp4, 146);
        assert_eq!(r3_hexagonal.ext, b'H');
        assert_eq!(r3_hexagonal.hall, "R 3");
        assert_eq!(r3_hexagonal.ops.len(), 9);
        assert_eq!(r3_hexagonal.ops[3].tran, [16, 8, 8]);
        assert_eq!(r3_hexagonal.ops[6].tran, [8, 16, 16]);

        let r3_rhombohedral = &GEMMI_SPACEGROUPS[433];
        assert_eq!(r3_rhombohedral.ccp4, 1146);
        assert_eq!(r3_rhombohedral.ext, b'R');
        assert_eq!(r3_rhombohedral.hall, "P 3*");
        assert_eq!(r3_rhombohedral.basisop_idx, 36);
        assert_eq!(
            r3_rhombohedral.ops,
            [
                BioSymOp {
                    rot: [[24, 0, 0], [0, 24, 0], [0, 0, 24]],
                    tran: [0, 0, 0],
                },
                BioSymOp {
                    rot: [[0, 0, 24], [24, 0, 0], [0, 24, 0]],
                    tran: [0, 0, 0],
                },
                BioSymOp {
                    rot: [[0, 24, 0], [0, 0, 24], [24, 0, 0]],
                    tran: [0, 0, 0],
                },
            ]
        );
    }
}
