//! Thin projection of registered molecule component operations.
use crate::Molecule;
impl Molecule {
    pub fn fragments(&self) -> Result<Vec<Self>, cosmolkit::OperationError> {
        self.inner.borrow().fragments().map(|values| {
            values
                .into_iter()
                .map(|inner| Self {
                    inner: std::cell::RefCell::new(inner),
                })
                .collect()
        })
    }
    pub fn largest_fragment(&self) -> Result<Self, cosmolkit::OperationError> {
        self.inner.borrow().largest_fragment().map(|inner| Self {
            inner: std::cell::RefCell::new(inner),
        })
    }
}
