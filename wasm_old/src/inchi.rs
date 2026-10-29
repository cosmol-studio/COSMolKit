//! InChI transport delegates every conversion to the public facade.
impl crate::Molecule {
    pub fn to_inchi_result(&self) -> Result<cosmolkit::InchiResult, cosmolkit::InchiError> {
        self.inner.borrow().to_inchi_result()
    }
    pub fn to_inchi_result_with_params(
        &self,
        params: &cosmolkit::InchiWriteParams,
    ) -> Result<cosmolkit::InchiResult, cosmolkit::InchiError> {
        self.inner.borrow().to_inchi_result_with_params(params)
    }
    pub fn from_inchi(text: &str) -> Result<Self, cosmolkit::InchiError> {
        cosmolkit::Molecule::from_inchi(text).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn from_inchi_with_params(
        text: &str,
        params: &cosmolkit::InchiReadParams,
    ) -> Result<Self, cosmolkit::InchiError> {
        cosmolkit::Molecule::from_inchi_with_params(text, params).map(|inner| Self {
            inner: inner.into(),
        })
    }
    pub fn to_inchi(&self) -> Result<String, cosmolkit::InchiError> {
        self.inner.borrow().to_inchi()
    }
    pub fn to_inchi_with_params(
        &self,
        params: &cosmolkit::InchiWriteParams,
    ) -> Result<String, cosmolkit::InchiError> {
        self.inner.borrow().to_inchi_with_params(params)
    }
    pub fn to_inchi_key(&self) -> Result<String, cosmolkit::InchiError> {
        self.inner.borrow().to_inchi_key()
    }
    pub fn to_inchi_key_with_params(
        &self,
        params: &cosmolkit::InchiWriteParams,
    ) -> Result<String, cosmolkit::InchiError> {
        self.inner.borrow().to_inchi_key_with_params(params)
    }
}
