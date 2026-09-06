//! Source-shaped hash helpers shared by Avalon fingerprint feature families.

pub(super) fn next_hash(mut hash: u64, data: u64) -> u64 {
    // Avalon❗✔️: uint64_t next_hash(uint64_t hash, uint64_t data)
    // Avalon❗✔️: {
    // Avalon❗✔️:     hash += data;
    hash = hash.wrapping_add(data);
    // Avalon❗✔️:     hash += (uint64_t)(hash << (uint64_t)10);
    hash = hash.wrapping_add(hash.wrapping_shl(10));
    // Avalon❗✔️:     hash ^= (uint64_t)(hash >> (uint64_t)6);
    hash ^= hash.wrapping_shr(6);
    // Avalon❗✔️:     return hash;
    hash
    // Avalon❗✔️: }
}

pub(super) fn hash_string(value: &str) -> u64 {
    // Avalon❗✔️: uint64_t hash_string(char *str)
    // Avalon❗✔️: {
    // Avalon❗✔️:    uint64_t hash = 1001;
    let mut hash = 1001_u64;
    // Avalon❗✔️:    for (; (*str); str++)
    // Avalon❗✔️:       hash = next_hash(hash, *str);
    for byte in value.bytes() {
        // The pinned Linux C ABI uses signed `char`; conversion to uint64_t
        // therefore sign-extends bytes above 0x7f before wrapping.
        let data = u64::from_ne_bytes(i64::from(byte as i8).to_ne_bytes());
        hash = next_hash(hash, data);
    }
    // Avalon❗✔️:    return hash;
    hash
    // Avalon❗✔️: }
}

pub(super) fn hash_position(mut hash: u64, slot_count: usize) -> usize {
    // Avalon❗✔️: uint64_t hash_position(uint64_t hash, int nslots)
    // Avalon❗✔️: {
    // Avalon❗✔️:     hash += (uint64_t)(hash << (uint64_t)3);
    hash = hash.wrapping_add(hash.wrapping_shl(3));
    // Avalon❗✔️:     hash ^= (uint64_t)(hash >> (uint64_t)11);
    hash ^= hash.wrapping_shr(11);
    // Avalon❗✔️:     hash += (uint64_t)(hash << (uint64_t)15);
    hash = hash.wrapping_add(hash.wrapping_shl(15));
    // Avalon❗✔️:     return (hash%(uint64_t)nslots);
    (hash % slot_count as u64) as usize
    // Avalon❗✔️: }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shortcut_label_hash_matches_native_ala_color_source_value() {
        assert_eq!((hash_string("ALA") & 0x00ff_ff00) | 119, 14_602_103);
    }

    #[test]
    fn next_hash_wraps_like_source_uint64_arithmetic() {
        assert_eq!(next_hash(u64::MAX, u64::MAX), 0xfc00_0000_0000_0821);
    }

    #[test]
    fn hash_position_uses_source_finalization_and_modulo() {
        assert_eq!(next_hash(17, 16), 34_353);
        assert_eq!(hash_position(next_hash(17, 16), 64), 47);
        assert_eq!(hash_position(u64::MAX, 257), 238);
    }
}
