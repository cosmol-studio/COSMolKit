//! Canonical lazy stereoisomer transport, including decimal big-integer values.
use crate::Molecule;
use cosmolkit as ck;
impl Molecule {
    pub fn stereoisomer_count_with_options(
        &self,
        options: &ck::StereoisomerOptions,
    ) -> Result<String, ck::EnumerationError> {
        self.inner
            .borrow()
            .stereoisomer_count_with_options(options)
            .map(|v| v.to_string())
    }
    pub fn enumerate_stereoisomers_with_options(
        &self,
        options: &ck::StereoisomerOptions,
    ) -> Result<ck::StereoisomerIterator, ck::OperationError> {
        self.inner
            .borrow()
            .enumerate_stereoisomers_with_options(options)
    }
    pub fn enumerate_stereoisomers_with_random_bits(
        &self,
        options: &ck::StereoisomerOptions,
        mut callback: Box<dyn FnMut(usize) -> Result<String, String> + Send + Sync + 'static>,
    ) -> Result<ck::StereoisomerIterator, ck::OperationError> {
        self.inner
            .borrow()
            .enumerate_stereoisomers_with_random_bits(
                options,
                Box::new(move |width| {
                    callback(width)?
                        .parse()
                        .map_err(|e| format!("invalid unsigned random bits: {e}"))
                }),
            )
    }
}
