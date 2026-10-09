//! Sole reset-only RDKit ControlC state shared by native source algorithms.
//! No Molecule operation, runtime storage or constructor/destructor policy.
//! A second SIGINT after the first uses SIG_DFL, as in the source reset-only
//! path whose d_prev_handler is zero initialized. Nested reset is global.

static SOURCE_INTERRUPTED: std::sync::atomic::AtomicBool =
    std::sync::atomic::AtomicBool::new(false);
pub fn got_signal() -> bool {
    // BEGIN RDKIT CPP FUNCTION source_got_signal (RDGeneral/ControlCHandler.h)
    // RDKit❗✔️:   static bool getGotSignal() { return d_gotSignal; }
    // END RDKIT CPP FUNCTION source_got_signal

    SOURCE_INTERRUPTED.load(std::sync::atomic::Ordering::SeqCst)
}
#[cfg(not(target_arch = "wasm32"))]
extern "C" fn source_signal_handler(signal_number: libc::c_int) {
    // BEGIN RDKIT CPP FUNCTION source_signal_handler (RDGeneral/ControlCHandler.h)
    // RDKit❗✔️:   static void signalHandler(int signalNumber) {
    // RDKit❗✔️:     if (signalNumber == SIGINT) {
    // RDKit❗✔️:       d_gotSignal = true;
    // RDKit❗✔️:       std::signal(SIGINT, d_prev_handler);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION source_signal_handler

    if signal_number == libc::SIGINT {
        SOURCE_INTERRUPTED.store(true, std::sync::atomic::Ordering::SeqCst);
        // Safety: this is the C signal() entrypoint with SIGINT and SIG_DFL, the
        // zero-initialized prior handler in the source reset-only embedding path.
        unsafe {
            libc::signal(libc::SIGINT, libc::SIG_DFL);
        }
    }
}
#[cfg(not(target_arch = "wasm32"))]
pub fn reset() -> bool {
    // BEGIN RDKIT CPP FUNCTION RDKit::ControlCHandler::reset_native
    // RDKit❗✔️:   static void reset() {
    // RDKit❗✔️:     d_gotSignal = false;
    // RDKit❗✔️:     std::signal(SIGINT, signalHandler);
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION RDKit::ControlCHandler::reset_native
    SOURCE_INTERRUPTED.store(false, std::sync::atomic::Ordering::SeqCst);
    // Safety: static C ABI callback, the existing lock-free flag and signal reset.
    let previous = unsafe {
        libc::signal(
            libc::SIGINT,
            source_signal_handler as *const () as libc::sighandler_t,
        )
    };
    // Windows declares SIG_ERR as c_int, while signal() returns sighandler_t.
    // Convert the -1 sentinel to the return type without narrowing its bits.
    previous != libc::SIG_ERR as libc::sighandler_t
}
#[cfg(target_arch = "wasm32")]
pub fn reset() -> bool {
    // BEGIN RDKIT CPP FUNCTION RDKit::ControlCHandler::reset_wasm
    // RDKit❗✔️:   static void reset() {
    // RDKit❗✔️:     d_gotSignal = false;
    // RDKit❗✔️:     std::signal(SIGINT, signalHandler);
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION RDKit::ControlCHandler::reset_wasm
    // Approved Web adaptation: browser/Node WASM has no process SIGINT handler.
    // Reset the same interruption state without installing an OS handler;
    // computation and deadline checks remain enabled. Host cancellation (e.g.
    // terminating a Worker) is separate from synchronous chemistry execution.
    SOURCE_INTERRUPTED.store(false, std::sync::atomic::Ordering::SeqCst);
    true
}
