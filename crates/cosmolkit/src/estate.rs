//! Molecular electrotopological state and its surface-area projections.
use crate::{DescriptorReadError, Molecule};
use cosmolkit_descriptors::{EStateExtrema, EStateVsaFamily};

impl Molecule {
    /// Return electrotopological state indices in current graph atom order.
    pub fn estate_indices(&self) -> Result<Vec<f64>, DescriptorReadError> {
        cosmolkit_descriptors::estate_indices(
            self.topology(),
            self.derived_cache_runtime().valence_assignment(),
        )
        .map_err(|source| DescriptorReadError::Algorithm { source })
    }

    fn estate_extrema(&self) -> Result<EStateExtrema, DescriptorReadError> {
        EStateExtrema::from_indices(&self.estate_indices()?)
            .map_err(|source| DescriptorReadError::Algorithm { source })
    }

    /// Return the maximum EState index; fail if the graph is empty.
    pub fn max_estate_index(&self) -> Result<f64, DescriptorReadError> {
        Ok(self.estate_extrema()?.maximum)
    }

    /// Return the minimum EState index; fail if the graph is empty.
    pub fn min_estate_index(&self) -> Result<f64, DescriptorReadError> {
        Ok(self.estate_extrema()?.minimum)
    }

    /// Return the maximum absolute atom EState index; fail if the graph is empty.
    pub fn max_abs_estate_index(&self) -> Result<f64, DescriptorReadError> {
        Ok(self.estate_extrema()?.maximum_absolute)
    }

    /// Return the minimum absolute atom EState index; fail if the graph is empty.
    pub fn min_abs_estate_index(&self) -> Result<f64, DescriptorReadError> {
        Ok(self.estate_extrema()?.minimum_absolute)
    }
}

macro_rules! numbered_bins {
    ($family:ident; $($name:ident => $index:literal),+ $(,)?) => {$(
        #[doc = concat!("Return the selected bin of the ", stringify!($family), " projection.")]
        pub fn $name(&self) -> Result<f64, DescriptorReadError> {
            self.project_estate_vsa(EStateVsaFamily::$family, false, |bins| bins[$index])
        }
    )+};
}

impl Molecule {
    fn project_estate_vsa<T>(
        &self,
        family: EStateVsaFamily,
        refresh: bool,
        project: impl FnOnce(&[f64]) -> T,
    ) -> Result<T, DescriptorReadError> {
        let mut state = self.descriptor_queries_runtime()?;
        let bins = family
            .bins(
                self.topology(),
                self.derived_cache_runtime().valence_assignment(),
                &mut state,
                refresh,
            )
            .map_err(|source| DescriptorReadError::Algorithm { source })?;
        Ok(project(bins))
    }

    /// Return the hydrogen-inclusive Labute atom areas summed in eleven EState bins.
    pub fn estate_vsa(&self) -> Result<Vec<f64>, DescriptorReadError> {
        self.project_estate_vsa(EStateVsaFamily::AreaByEState, true, <[f64]>::to_vec)
    }

    /// Return EState indices summed in ten Labute atom-area bins.
    pub fn vsa_estate(&self) -> Result<Vec<f64>, DescriptorReadError> {
        self.project_estate_vsa(EStateVsaFamily::EStateByArea, true, <[f64]>::to_vec)
    }

    numbered_bins!(AreaByEState; estate_vsa_1 => 0, estate_vsa_2 => 1, estate_vsa_3 => 2, estate_vsa_4 => 3, estate_vsa_5 => 4, estate_vsa_6 => 5, estate_vsa_7 => 6, estate_vsa_8 => 7, estate_vsa_9 => 8, estate_vsa_10 => 9, estate_vsa_11 => 10);
    numbered_bins!(EStateByArea; vsa_estate_1 => 0, vsa_estate_2 => 1, vsa_estate_3 => 2, vsa_estate_4 => 3, vsa_estate_5 => 4, vsa_estate_6 => 5, vsa_estate_7 => 6, vsa_estate_8 => 7, vsa_estate_9 => 8, vsa_estate_10 => 9);
}
