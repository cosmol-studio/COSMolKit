//! Sole RDKit ControlC state shared by native source algorithms.
//! No Molecule operation or runtime storage. Native scoped and reset-only
//! consumers share the same source static previous-handler semantics.

static SOURCE_INTERRUPTED: std::sync::atomic::AtomicBool =
    std::sync::atomic::AtomicBool::new(false);
#[cfg(not(target_arch = "wasm32"))]
static PREVIOUS_HANDLER: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(libc::SIG_DFL as usize);

/// Internal source-lifetime guard; Web targets retain the accepted no-OS-signal policy.
#[doc(hidden)]
pub struct ControlCHandler;

impl ControlCHandler {
    pub fn new() -> Option<Self> {
        // RDKit✔️✔️:   ControlCHandler() { d_prev_handler = std::signal(SIGINT, signalHandler); }
        // One process-global previous slot, exactly as in upstream (not a
        // per-instance stack). Atomic storage keeps Rust signal access sound.
        // No allocation or mutex in signal context. Constructor does not reset
        // an existing flag; the source only clears it in reset/destruction.
        #[cfg(not(target_arch = "wasm32"))]
        {
            // SAFETY: static C ABI callback; signal returns the opaque previous
            // handler token, which is stored/reinstalled without dereferencing.
            let previous = unsafe {
                libc::signal(
                    libc::SIGINT,
                    source_signal_handler as *const () as libc::sighandler_t,
                )
            };
            if previous == libc::SIG_ERR as libc::sighandler_t {
                return None;
            }
            PREVIOUS_HANDLER.store(previous as usize, std::sync::atomic::Ordering::SeqCst);
        }
        // WASM has no process SIGINT: only scoped flag lifetime is projected.
        Some(Self)
    }
}

impl Drop for ControlCHandler {
    fn drop(&mut self) {
        // RDKit✔️✔️:   ~ControlCHandler() {
        // RDKit✔️✔️:     std::signal(SIGINT, d_prev_handler);
        // RDKit✔️✔️:     d_gotSignal = false;
        // RDKit✔️✔️:   }
        #[cfg(not(target_arch = "wasm32"))]
        // SAFETY: this is the unchanged opaque token returned by signal(), not
        // a Rust data pointer or a dynamically allocated callback.
        unsafe {
            libc::signal(
                libc::SIGINT,
                PREVIOUS_HANDLER.load(std::sync::atomic::Ordering::SeqCst) as libc::sighandler_t,
            );
        }
        SOURCE_INTERRUPTED.store(false, std::sync::atomic::Ordering::SeqCst);
    }
}
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
        // Safety: re-install the source static previous handler token.
        unsafe {
            libc::signal(
                libc::SIGINT,
                PREVIOUS_HANDLER.load(std::sync::atomic::Ordering::SeqCst) as libc::sighandler_t,
            );
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

#[cfg(all(test, not(target_arch = "wasm32")))]
mod tests {
    #[test]
    fn scoped_handler_restores_native_handler_and_clears_signal() {
        const CHILD: &str = "CK_SCOPED_SIGNAL_REGRESSION_CHILD";
        if std::env::var_os(CHILD).is_none() {
            let status = std::process::Command::new(std::env::current_exe().unwrap())
                .args(["--exact", "source_control_c::tests::scoped_handler_restores_native_handler_and_clears_signal", "--nocapture"])
                .env(CHILD, "1").status().unwrap();
            assert!(status.success());
            return;
        }
        // Isolated child: never change the test runner's global SIGINT policy.
        unsafe {
            libc::signal(libc::SIGINT, libc::SIG_IGN);
        }
        {
            let _handler = super::ControlCHandler::new().unwrap();
            assert!(!super::got_signal());
            assert_eq!(unsafe { libc::raise(libc::SIGINT) }, 0);
            assert!(super::got_signal());
        }
        assert!(!super::got_signal());
        assert_eq!(
            unsafe { libc::signal(libc::SIGINT, libc::SIG_IGN) },
            libc::SIG_IGN
        );
    }
}
