//! String scan cursor corresponding to musl's string pseudo-FILE.
//!
//! Source: `third_party/musl/src/internal/shgetc.h` and `shgetc.c` (1.2.5),
//! string-only path. `sh_fromstring` sets `buf = rpos = s`, `rend = (void*)-1`
//! and `shlim(f, 0)` makes `shend = rend`; `shgetc` then reads `*rpos++` and
//! `shunget` decrements `rpos` while `shlim >= 0`. `shcnt` reports
//! `shcnt + (rpos - buf)`; for the string path `shcnt` starts as `buf - rpos`
//! so the count is the number of bytes advanced, and `shlim(f, 0)` resets it.
//!
//! Only the string path is ported. FILE stream IO, `__uflow`, locking and the
//! limit machinery for nonzero limits are not reachable here.

/// A byte cursor over a NUL-terminated string view.
///
/// Reading at or past the end yields the virtual terminating NUL (`0`) and
/// reproduces the source's `shunget` undo behavior for that terminator. The
/// consumed count starts at zero and `shlim(f, 0)` is modeled by
/// [`ScanCursor::reset_count`].
pub(crate) struct ScanCursor<'a> {
    bytes: &'a [u8],
    pos: usize,
    start: usize,
}

impl<'a> ScanCursor<'a> {
    pub(crate) const fn new(bytes: &'a [u8]) -> Self {
        // musl✔️✔️: #define sh_fromstring(f, s) \
        // musl✔️✔️: 	((f)->buf = (f)->rpos = (void *)(s), (f)->rend = (void*)-1)
        // String-only initialization; no allocation or FILE stream machinery.
        Self {
            bytes,
            pos: 0,
            start: 0,
        }
    }

    /// `shgetc(f)`: return the current byte and advance, or the virtual NUL.
    pub(crate) fn getc(&mut self) -> i32 {
        // musl✔️✔️: #define shgetc(f) (((f)->rpos != (f)->shend) ? *(f)->rpos++ : __shgetc(f))
        // On the string path shend is (void*)-1, so the buffered branch is
        // always used. Bounds checks provide a virtual NUL instead of reading
        // outside the slice; no conversion path needs bytes beyond that NUL.
        let byte = self.bytes.get(self.pos).copied().unwrap_or(0);
        if self.pos <= self.bytes.len() {
            self.pos += 1;
        }
        i32::from(byte)
    }

    /// `shunget(f)`: move the read position back one byte.
    pub(crate) fn ungetc(&mut self) {
        // musl✔️✔️: #define shunget(f) ((f)->shlim>=0 ? (void)(f)->rpos-- : (void)0)
        // The string specialization always has shlim == 0.
        if self.pos > 0 {
            self.pos -= 1;
        }
    }

    /// `shlim(f, 0)`: reset the consumed-byte count without moving the cursor.
    pub(crate) fn reset_count(&mut self) {
        // musl✔️✔️: void __shlim(FILE *f, off_t lim)
        // musl✔️✔️: {
        // musl✔️✔️: 	f->shlim = lim;
        // musl✔️✔️: 	f->shcnt = f->buf - f->rpos;
        // musl✔️✔️: 	/* If lim is nonzero, rend must be a valid pointer. */
        // musl✔️✔️: 	if (lim && f->rend - f->rpos > lim)
        // musl✔️✔️: 		f->shend = f->rpos + lim;
        // musl✔️✔️: 	else
        // musl✔️✔️: 		f->shend = f->rend;
        // musl✔️✔️: }
        // Only lim == 0 is reachable; pos - start is the source offset sum.
        self.start = self.pos;
    }

    /// `shcnt(f)`: bytes consumed since construction or the last reset.
    pub(crate) fn consumed(&self) -> usize {
        // musl✔️✔️: #define shcnt(f) ((f)->shcnt + ((f)->rpos - (f)->buf))
        self.pos - self.start
    }
}

#[cfg(test)]
mod tests {
    use super::ScanCursor;

    #[test]
    fn atof_port_cursor_reads_bytes_then_virtual_nul() {
        let mut cursor = ScanCursor::new(b"1e");
        assert_eq!(cursor.getc(), i32::from(b'1'));
        assert_eq!(cursor.getc(), i32::from(b'e'));
        assert_eq!(cursor.getc(), 0);
        assert_eq!(cursor.getc(), 0);
        // The first virtual NUL advances to len+1; the clamp stops there.
        assert_eq!(cursor.consumed(), 3);
        cursor.ungetc();
        cursor.ungetc();
        assert_eq!(cursor.consumed(), 1);
    }

    #[test]
    fn atof_port_cursor_ungetc_after_virtual_nul_does_not_advance() {
        let mut cursor = ScanCursor::new(b"x");
        assert_eq!(cursor.getc(), i32::from(b'x'));
        assert_eq!(cursor.getc(), 0);
        cursor.ungetc();
        assert_eq!(cursor.getc(), 0);
        assert_eq!(cursor.consumed(), 2);
    }

    #[test]
    fn atof_port_cursor_embedded_nul_stops_at_zero_byte() {
        let mut cursor = ScanCursor::new(b"1\0002");
        assert_eq!(cursor.getc(), i32::from(b'1'));
        assert_eq!(cursor.getc(), 0);
        assert_eq!(cursor.consumed(), 2);
    }

    #[test]
    fn atof_port_cursor_reset_count_matches_shlim_zero() {
        let mut cursor = ScanCursor::new(b"abc");
        let _ = cursor.getc();
        let _ = cursor.getc();
        assert_eq!(cursor.consumed(), 2);
        cursor.reset_count();
        assert_eq!(cursor.consumed(), 0);
    }

    #[test]
    fn atof_port_cursor_empty_input_is_immediately_terminal() {
        let mut cursor = ScanCursor::new(b"");
        assert_eq!(cursor.getc(), 0);
        assert_eq!(cursor.consumed(), 1);
    }
}
