//! Checked native field ownership. Views resolve a path on every access.
use std::{cell::RefCell, rc::Rc};
use wasm_bindgen::JsValue;

/// Shared immutable native value, optionally projected from an owning parent.
/// Retains non-Clone values without reconstructing or copying their contents.
/// Unlike NativeField, this container provides no mutable access.
pub(crate) struct ReadonlyValue<T: ?Sized>(Rc<dyn AsRef<T>>);
impl<T: ?Sized> Clone for ReadonlyValue<T> {
    fn clone(&self) -> Self {
        Self(self.0.clone())
    }
}
impl<T: ?Sized> AsRef<T> for ReadonlyValue<T> {
    fn as_ref(&self) -> &T {
        self.0.as_ref().as_ref()
    }
}
impl<T: ?Sized> std::ops::Deref for ReadonlyValue<T> {
    type Target = T;
    fn deref(&self) -> &T {
        self.as_ref()
    }
}
struct ReadonlyOwned<T>(T);
impl<T> AsRef<T> for ReadonlyOwned<T> {
    fn as_ref(&self) -> &T {
        &self.0
    }
}
struct ReadonlyChild<P: ?Sized, T: ?Sized> {
    parent: ReadonlyValue<P>,
    get: Box<dyn for<'a> Fn(&'a P) -> &'a T>,
}
impl<P: ?Sized, T: ?Sized> AsRef<T> for ReadonlyChild<P, T> {
    fn as_ref(&self) -> &T {
        (self.get)(self.parent.as_ref())
    }
}
impl<T: 'static> ReadonlyValue<T> {
    pub(crate) fn owned(value: T) -> Self {
        Self(Rc::new(ReadonlyOwned(value)))
    }
}
impl<T: ?Sized + 'static> ReadonlyValue<T> {
    pub(crate) fn child<U: ?Sized + 'static>(
        &self,
        get: impl for<'a> Fn(&'a T) -> &'a U + 'static,
    ) -> ReadonlyValue<U> {
        ReadonlyValue(Rc::new(ReadonlyChild {
            parent: self.clone(),
            get: Box::new(get),
        }))
    }
}

/// Borrow a mutable exported owner through its checked class ABI. This is not
/// a configuration conversion: unique native executor/output state stays put.
macro_rules! mutable_owner {
    ($wrapper:ident, $native:ty, $input:ident, $ts:literal) => {
        #[wasm_bindgen::prelude::wasm_bindgen]
        extern "C" {
            #[wasm_bindgen(typescript_type = $ts)]
            pub type $input;
        }
        const _: () = {
            #[wasm_bindgen::prelude::wasm_bindgen]
            extern "C" {
                #[wasm_bindgen(catch, js_namespace = Reflect, js_name = apply)]
                fn visit_mut(
                    visit: &mut dyn FnMut(&$wrapper),
                    receiver: &wasm_bindgen::JsValue,
                    args: &js_sys::Array,
                ) -> Result<wasm_bindgen::JsValue, wasm_bindgen::JsValue>;
            }
            impl $wrapper {
                pub(crate) fn with_mut<R>(
                    value: &$input,
                    apply: impl FnOnce(&mut $native) -> Result<R, wasm_bindgen::JsValue>,
                ) -> Result<R, wasm_bindgen::JsValue> {
                    crate::input::with_instance(
                        value.as_ref(),
                        $ts,
                        |visit, receiver, args| {
                            visit_mut(&mut |owner| visit(owner), receiver, args)
                        },
                        |owner| owner.0.write(apply)?,
                    )
                }
            }
        };
    };
}
pub(crate) use mutable_owner;

trait Storage<T> {
    fn read(&self, visit: &mut dyn FnMut(&T)) -> Result<(), JsValue>;
    fn write(&self, visit: &mut dyn FnMut(&mut T)) -> Result<(), JsValue>;
}

struct Owned<T>(RefCell<T>);
impl<T> Storage<T> for Owned<T> {
    fn read(&self, visit: &mut dyn FnMut(&T)) -> Result<(), JsValue> {
        let value = self
            .0
            .try_borrow()
            .map_err(|_| js_sys::TypeError::new("configuration is currently being modified"))?;
        visit(&value);
        Ok(())
    }
    fn write(&self, visit: &mut dyn FnMut(&mut T)) -> Result<(), JsValue> {
        let mut value = self
            .0
            .try_borrow_mut()
            .map_err(|_| js_sys::TypeError::new("configuration is currently borrowed"))?;
        visit(&mut value);
        Ok(())
    }
}

pub(crate) struct NativeField<T>(Rc<dyn Storage<T>>);
impl<T> Clone for NativeField<T> {
    fn clone(&self) -> Self {
        Self(self.0.clone())
    }
}

impl<T: 'static> NativeField<T> {
    pub fn owned(value: T) -> Self {
        Self(Rc::new(Owned(RefCell::new(value))))
    }

    pub fn read<R>(&self, visit: impl FnOnce(&T) -> R) -> Result<R, JsValue> {
        let mut visit = Some(visit);
        let mut result = None;
        self.0
            .read(&mut |value| result = Some(visit.take().unwrap()(value)))?;
        Ok(result.unwrap())
    }

    pub fn write<R>(&self, visit: impl FnOnce(&mut T) -> R) -> Result<R, JsValue> {
        let mut visit = Some(visit);
        let mut result = None;
        self.0
            .write(&mut |value| result = Some(visit.take().unwrap()(value)))?;
        Ok(result.unwrap())
    }

    pub fn child<U: 'static>(
        &self,
        get: fn(&T) -> Option<&U>,
        get_mut: fn(&mut T) -> Option<&mut U>,
    ) -> NativeField<U> {
        NativeField(Rc::new(Child {
            parent: self.clone(),
            get,
            get_mut: Some(get_mut),
        }))
    }

    /// Preserve a native shared borrow without granting a mutable child path.
    pub fn readonly_child<U: 'static>(&self, get: fn(&T) -> Option<&U>) -> NativeField<U> {
        NativeField(Rc::new(Child {
            parent: self.clone(),
            get,
            get_mut: None,
        }))
    }
}

struct Child<P, T> {
    parent: NativeField<P>,
    get: fn(&P) -> Option<&T>,
    get_mut: Option<fn(&mut P) -> Option<&mut T>>,
}
impl<P: 'static, T> Storage<T> for Child<P, T> {
    fn read(&self, visit: &mut dyn FnMut(&T)) -> Result<(), JsValue> {
        self.parent.read(|parent| {
            let child = (self.get)(parent).ok_or_else(|| {
                js_sys::ReferenceError::new("the optional configuration field is absent")
            })?;
            visit(child);
            Ok(())
        })?
    }
    fn write(&self, visit: &mut dyn FnMut(&mut T)) -> Result<(), JsValue> {
        let get_mut = self.get_mut.ok_or_else(|| {
            js_sys::TypeError::new("this native view is readonly; copy it before editing")
        })?;
        self.parent.write(|parent| {
            let child = get_mut(parent).ok_or_else(|| {
                js_sys::ReferenceError::new("the optional configuration field is absent")
            })?;
            visit(child);
            Ok(())
        })?
    }
}
