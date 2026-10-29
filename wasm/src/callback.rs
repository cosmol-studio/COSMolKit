//! Synchronous native callback boundary. JS values never enter Send closures.
use std::{
    cell::RefCell,
    collections::HashMap,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
    thread::ThreadId,
};
use wasm_bindgen::{JsCast, JsValue};

struct Entry {
    function: js_sys::Function,
    native: Option<usize>,
}
#[derive(Default)]
struct Registry {
    next: u64,
    entries: HashMap<u64, Entry>,
    native: HashMap<usize, u64>,
}
struct Frame {
    error: Option<JsValue>,
    targets: Option<Vec<cosmolkit::Molecule>>,
}
thread_local! {
    static REGISTRY: RefCell<Registry> = RefCell::new(Registry::default());
    static FRAMES: RefCell<Vec<Frame>> = const { RefCell::new(Vec::new()) };
    static WRONG_THREAD: Arc<AtomicBool> = Arc::new(AtomicBool::new(false));
}
/// The callback's native owner holds this Send handle, never the JS function.
pub(crate) struct Handle {
    id: u64,
    owner: ThreadId,
    wrong_thread: Arc<AtomicBool>,
}
impl Drop for Handle {
    fn drop(&mut self) {
        if self.owner != std::thread::current().id() {
            return;
        }
        REGISTRY.with(|registry| {
            let mut registry = registry.borrow_mut();
            if let Some(entry) = registry.entries.remove(&self.id) {
                if let Some(native) = entry.native {
                    registry.native.remove(&native);
                }
            }
        });
    }
}
impl Handle {
    pub(crate) fn new(raw: &JsValue) -> Result<Arc<Self>, JsValue> {
        let function = raw
            .dyn_ref::<js_sys::Function>()
            .ok_or_else(|| js_sys::TypeError::new("expected a synchronous callback"))?
            .clone();
        let id = REGISTRY.with(|registry| {
            let mut registry = registry.borrow_mut();
            registry.next = registry
                .next
                .checked_add(1)
                .ok_or_else(|| js_sys::Error::new("callback identifiers exhausted"))?;
            let id = registry.next;
            registry.entries.insert(
                id,
                Entry {
                    function,
                    native: None,
                },
            );
            Ok::<_, JsValue>(id)
        })?;
        Ok(Arc::new(Self {
            id,
            owner: std::thread::current().id(),
            wrong_thread: WRONG_THREAD.with(Arc::clone),
        }))
    }
    pub(crate) fn bind<T: ?Sized>(&self, native: &Arc<T>) {
        let pointer = Arc::as_ptr(native).cast::<()>() as usize;
        REGISTRY.with(|registry| {
            let mut registry = registry.borrow_mut();
            registry.native.insert(pointer, self.id);
            registry
                .entries
                .get_mut(&self.id)
                .expect("live callback handle")
                .native = Some(pointer);
        });
    }
    pub(crate) fn ensure_thread(&self) -> Result<(), JsValue> {
        if self.owner != std::thread::current().id() {
            self.wrong_thread.store(true, Ordering::Relaxed);
            // No JS value may be accessed from a different native thread.
            return Err(JsValue::UNDEFINED);
        }
        Ok(())
    }
    pub(crate) fn call(
        &self,
        arguments: impl FnOnce() -> Result<js_sys::Array, JsValue>,
    ) -> Result<JsValue, JsValue> {
        self.ensure_thread()?;
        let failure =
            FRAMES.with(|frames| frames.borrow().iter().rev().find_map(|f| f.error.clone()));
        if let Some(error) = failure {
            return Err(error);
        }
        let function = REGISTRY
            .with(|registry| {
                registry
                    .borrow()
                    .entries
                    .get(&self.id)
                    .map(|e| e.function.clone())
            })
            .ok_or_else(|| js_sys::Error::new("callback owner is unavailable"))?;
        function.apply(&JsValue::UNDEFINED, &arguments()?)
    }
    pub(crate) fn boolean(
        &self,
        arguments: impl FnOnce() -> Result<js_sys::Array, JsValue>,
    ) -> bool {
        match self
            .call(arguments)
            .and_then(|v| crate::input::read::<bool>(tsify::Ts::new_unchecked(v)))
        {
            Ok(value) => value,
            Err(error) => {
                record(error);
                false
            }
        }
    }
}
pub(crate) fn original<T: ?Sized>(native: &Arc<T>) -> Result<JsValue, JsValue> {
    let pointer = Arc::as_ptr(native).cast::<()>() as usize;
    REGISTRY.with(|registry| {
        let registry = registry.borrow();
        let entry = registry
            .native
            .get(&pointer)
            .and_then(|id| registry.entries.get(id));
        entry
            .map(|entry| JsValue::from(entry.function.clone()))
            .ok_or_else(|| js_sys::Error::new("native callback has no JavaScript owner").into())
    })
}
pub(crate) fn record(error: JsValue) {
    FRAMES.with(|frames| {
        if let Some(frame) = frames.borrow_mut().last_mut() {
            if frame.error.is_none() {
                frame.error = Some(error);
            }
        }
    });
}
struct Scope;
impl Drop for Scope {
    fn drop(&mut self) {
        FRAMES.with(|frames| {
            frames.borrow_mut().pop();
        });
    }
}
pub(crate) fn run<R>(call: impl FnOnce() -> Result<R, JsValue>) -> Result<R, JsValue> {
    run_with_targets(None, call)
}
pub(crate) fn run_with_targets<R>(
    targets: Option<Vec<cosmolkit::Molecule>>,
    call: impl FnOnce() -> Result<R, JsValue>,
) -> Result<R, JsValue> {
    FRAMES.with(|frames| {
        frames.borrow_mut().push(Frame {
            error: None,
            targets,
        })
    });
    let scope = Scope;
    let result = call();
    let error = FRAMES.with(|frames| frames.borrow().last().and_then(|frame| frame.error.clone()));
    drop(scope);
    if WRONG_THREAD.with(|flag| flag.swap(false, Ordering::Relaxed)) {
        return Err(js_sys::Error::new(
            "JavaScript callbacks require the owning synchronous thread",
        )
        .into());
    }
    match error {
        Some(error) => Err(error),
        None => result,
    }
}
pub(crate) fn target(index: Option<usize>) -> Result<cosmolkit::Molecule, JsValue> {
    FRAMES.with(|frames| {
        let frames = frames.borrow();
        let targets = frames
            .iter()
            .rev()
            .find_map(|frame| frame.targets.as_ref())
            .ok_or_else(|| js_sys::Error::new("matcher callback has no active target scope"))?;
        let index = if targets.len() == 1 {
            0
        } else {
            index.ok_or_else(|| {
                js_sys::Error::new("native multi-target callback omitted its input index")
            })?
        };
        targets.get(index).cloned().ok_or_else(|| {
            js_sys::Error::new("native callback target index is out of range").into()
        })
    })
}
