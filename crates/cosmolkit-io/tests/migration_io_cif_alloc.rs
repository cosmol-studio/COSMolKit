//! BIO-CID C20: single-final-allocation proof for the composed CIF number
//! formatter on the RUST path (exactly one heap allocation per call: the
//! returned String; no intermediate strings). This is NOT a cross-language
//! parity measurement: the source to_str returns std::string(buf, len),
//! whose allocation count is ABI/SSO-dependent (host libstdc++15 inline
//! capacity 15 — basic_string.h:218 / basic_string.tcc:233 — so 0
//! allocations for outputs up to 15 bytes, 1 beyond) and was not measured.
//! Count only the calling thread during each formatter call. Even one test
//! runs alongside libtest's main thread, whose bookkeeping may allocate.

use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;
use std::sync::atomic::{AtomicBool, Ordering};

use cosmolkit_io::cif::format_cif_f64;

thread_local! {
    // Const initialization needs no allocation or destructor. None disables
    // measurement; Some(count) records this thread's active call only.
    static ALLOCATION_COUNT: Cell<Option<usize>> = const { Cell::new(None) };
}

fn record_allocation() {
    // Allocator hooks must not allocate or unwind, including during TLS teardown.
    let _ = ALLOCATION_COUNT.try_with(|count| {
        if let Some(value) = count.get() {
            count.set(Some(value.saturating_add(1)));
        }
    });
}

struct CountingAllocator;

unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        record_allocation();
        unsafe { System.alloc(layout) }
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        record_allocation();
        unsafe { System.alloc_zeroed(layout) }
    }
    unsafe fn realloc(&self, pointer: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        record_allocation();
        unsafe { System.realloc(pointer, layout, new_size) }
    }
    unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
        unsafe { System.dealloc(pointer, layout) }
    }
}

#[global_allocator]
static GLOBAL: CountingAllocator = CountingAllocator;

#[test]
fn bio_cid_c20_allocation_count() {
    // One row per dispatch branch: zero, negative zero, integer fixed,
    // fractional fixed, fixed with zero-fill, exponent small/large,
    // subnormal, f64::MAX, and each special spelling. Each RUST call must
    // make EXACTLY one heap allocation (the returned String; Rust String
    // has no SSO): the significand carrier is stack state and both notation
    // sinks are bounded stack buffers (BIO-CID C16-C19). No claim is made
    // about the source std::string allocation count (ABI/SSO-dependent).
    let rows: [(f64, &str); 13] = [
        (0.0, "0"),
        (-0.0, "-0"),
        (1023.0, "1023"),
        (1234.5, "1234.5"),
        (100.0, "100"),
        (1.23456789, "1.23456789"),
        (1e-5, "1e-05"),
        (1.23456789e18, "1.23456789e+18"),
        (5e-324, "4.94065646e-324"),
        (f64::MAX, "1.79769313e+308"),
        (f64::INFINITY, "Inf"),
        (f64::NEG_INFINITY, "-Inf"),
        (f64::NAN, "NaN"),
    ];
    for (value, expected) in rows {
        ALLOCATION_COUNT.set(Some(0));
        let output = format_cif_f64(std::hint::black_box(value));
        let delta = ALLOCATION_COUNT.take().unwrap();
        assert_eq!(delta, 1, "allocation delta for {value:?} -> {output}");
        assert_eq!(output, expected);
    }
}

#[test]
fn allocation_counter_ignores_other_threads_and_inactive_scope() {
    std::hint::black_box(vec![0_u8; 64]);
    assert_eq!(ALLOCATION_COUNT.get(), None);

    let start = AtomicBool::new(false);
    let finished = AtomicBool::new(false);
    std::thread::scope(|scope| {
        let worker = scope.spawn(|| {
            while !start.load(Ordering::Acquire) {
                std::hint::spin_loop();
            }
            std::hint::black_box(vec![0_u8; 64]);
            finished.store(true, Ordering::Release);
        });

        ALLOCATION_COUNT.set(Some(0));
        start.store(true, Ordering::Release);
        while !finished.load(Ordering::Acquire) {
            std::hint::spin_loop();
        }
        // The background allocation happened inside this measurement window,
        // but only this formatter's returned String belongs in the count.
        let output = format_cif_f64(std::hint::black_box(0.0));
        let count = ALLOCATION_COUNT.take().unwrap();
        assert_eq!(count, 1);
        assert_eq!(output, "0");
        worker.join().unwrap();
    });
    assert_eq!(ALLOCATION_COUNT.get(), None);
}
