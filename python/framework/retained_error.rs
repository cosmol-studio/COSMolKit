//! A typed view of an actual retained native error, including non-Clone causes.
use crate::io_source_errors::ErrorRef;
use std::error::Error;
use std::marker::PhantomData;
use std::ops::Deref;
use std::sync::Arc;

pub(crate) struct RetainedError<T> {
    node: ErrorRef,
    marker: PhantomData<T>,
}
impl<T> Clone for RetainedError<T> {
    fn clone(&self) -> Self {
        Self {
            node: self.node.clone(),
            marker: PhantomData,
        }
    }
}
impl<T: Error + Send + Sync + 'static> RetainedError<T> {
    pub(crate) fn new(inner: T) -> Self {
        Self::from_ref(ErrorRef::new(Arc::new(inner)))
    }
    pub(crate) fn from_ref(node: ErrorRef) -> Self {
        assert!(node.error().is::<T>(), "concrete native error view type");
        Self {
            node,
            marker: PhantomData,
        }
    }
    pub(crate) fn node(&self) -> &ErrorRef {
        &self.node
    }
}
impl<T: Error + 'static> Deref for RetainedError<T> {
    type Target = T;
    fn deref(&self) -> &T {
        self.node
            .error()
            .downcast_ref()
            .expect("retained concrete native error")
    }
}
impl<T: Error + 'static> AsRef<T> for RetainedError<T> {
    fn as_ref(&self) -> &T {
        self
    }
}
