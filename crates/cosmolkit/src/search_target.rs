use crate::Molecule;

impl Molecule {
    pub(crate) fn detached_search_target(&self) -> cosmolkit_search::SearchTarget<'_> {
        let cache = self.derived_cache_runtime();
        cosmolkit_search::SearchTarget::new(
            self.topology(),
            self.coordinate_block_runtime(),
            &self.topology().stereo_groups,
            cache.valid_ring_info(),
            cache.valence_assignment(),
        )
    }
}
