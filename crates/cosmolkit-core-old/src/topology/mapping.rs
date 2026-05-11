#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct IndexMapping {
    old_to_new: Vec<Option<usize>>,
    new_to_old: Vec<Option<usize>>,
}

impl IndexMapping {
    pub(crate) fn new(old_to_new: Vec<Option<usize>>, new_to_old: Vec<Option<usize>>) -> Self {
        Self {
            old_to_new,
            new_to_old,
        }
    }

    pub(crate) fn from_kept_indices(old_len: usize, kept_old_indices: &[usize]) -> Self {
        let mut old_to_new = vec![None; old_len];
        let mut new_to_old = Vec::with_capacity(kept_old_indices.len());

        for (new_idx, &old_idx) in kept_old_indices.iter().enumerate() {
            old_to_new[old_idx] = Some(new_idx);
            new_to_old.push(Some(old_idx));
        }

        Self {
            old_to_new,
            new_to_old,
        }
    }

    pub(crate) fn old_to_new(&self, old_idx: usize) -> Option<usize> {
        self.old_to_new.get(old_idx).copied().flatten()
    }

    pub(crate) fn new_to_old(&self, new_idx: usize) -> Option<usize> {
        self.new_to_old.get(new_idx).copied().flatten()
    }

    pub(crate) fn old_len(&self) -> usize {
        self.old_to_new.len()
    }

    pub(crate) fn new_len(&self) -> usize {
        self.new_to_old.len()
    }

    pub(crate) fn is_deleted_old(&self, old_idx: usize) -> bool {
        self.old_to_new(old_idx).is_none()
    }

    pub(crate) fn is_newly_created(&self, new_idx: usize) -> bool {
        self.new_to_old(new_idx).is_none()
    }
}

pub(crate) type AtomMapping = IndexMapping;
pub(crate) type BondMapping = IndexMapping;
