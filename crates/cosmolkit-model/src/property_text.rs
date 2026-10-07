//! Canonical byte text carried by source properties and chemistry notation.

use std::borrow::Borrow;

/// Owned source text. Every byte, including NUL and non-UTF-8 bytes, is data.
///
/// Chemistry reads and writes bytes. Language projections perform their own
/// explicit text decoding after chemical storage and notation output.
#[doc(hidden)]
#[derive(Debug, Clone, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PropertyText(Vec<u8>);

impl PropertyText {
    #[must_use]
    pub const fn new() -> Self {
        Self(Vec::new())
    }

    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self(Vec::with_capacity(capacity))
    }

    pub fn insert_byte(&mut self, index: usize, byte: u8) {
        self.0.insert(index, byte);
    }

    /// Copy a borrowed byte sequence into the canonical owning text value.
    #[must_use]
    pub fn from_bytes(bytes: &[u8]) -> Self {
        // RDKit✔️✔️: inline RDValue(const std::string &v)
        // RDKit✔️✔️:       : value(new std::string(v)), type(RDTypeTag::StringTag) {}
        // Behavior: the source String payload owns a byte-for-byte copy; no
        // byte is a Unicode scalar or a terminator inside this counted value.
        // PropertyValue owns the source tag, this type owns its text payload.
        // Complexity: one owned buffer copy, O(bytes), matching std::string's
        // counted copy. No validity scan, secondary representation or cache.
        Self(bytes.to_vec())
    }

    #[must_use]
    pub const fn as_bytes(&self) -> &[u8] {
        self.0.as_slice()
    }

    #[must_use]
    pub fn into_bytes(self) -> Vec<u8> {
        self.0
    }

    #[must_use]
    pub const fn len(&self) -> usize {
        self.0.len()
    }

    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn push_byte(&mut self, byte: u8) {
        self.0.push(byte);
    }

    pub fn extend_bytes(&mut self, bytes: &[u8]) {
        self.0.extend_from_slice(bytes);
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }
}

impl AsRef<[u8]> for PropertyText {
    fn as_ref(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl Borrow<[u8]> for PropertyText {
    fn borrow(&self) -> &[u8] {
        self.as_bytes()
    }
}

impl From<Vec<u8>> for PropertyText {
    fn from(bytes: Vec<u8>) -> Self {
        // Owned Rust inputs transfer the same counted buffer. This is a data
        // transport adapter, not a second string tag or a text conversion.
        Self(bytes)
    }
}

impl From<String> for PropertyText {
    fn from(text: String) -> Self {
        Self(text.into_bytes())
    }
}

impl From<&str> for PropertyText {
    fn from(text: &str) -> Self {
        Self::from_bytes(text.as_bytes())
    }
}

impl From<&String> for PropertyText {
    fn from(text: &String) -> Self {
        Self::from_bytes(text.as_bytes())
    }
}

impl From<&[u8]> for PropertyText {
    fn from(bytes: &[u8]) -> Self {
        Self::from_bytes(bytes)
    }
}

impl<const N: usize> From<&[u8; N]> for PropertyText {
    fn from(bytes: &[u8; N]) -> Self {
        Self::from_bytes(bytes)
    }
}

impl From<&PropertyText> for PropertyText {
    fn from(text: &PropertyText) -> Self {
        text.clone()
    }
}

// Rust formatting writes already-valid generated text into the one counted
// byte buffer. Raw property payloads use extend_bytes directly; this adapter
// neither decodes bytes nor defines a Display implementation for them.
impl std::fmt::Write for PropertyText {
    fn write_str(&mut self, text: &str) -> std::fmt::Result {
        self.extend_bytes(text.as_bytes());
        Ok(())
    }
}
