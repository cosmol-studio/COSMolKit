//! Shared RDKit conformer-ID selection.

use crate::Conformer3D;

pub(crate) fn resolve_3d_conformer_index(conformers: &[Conformer3D], id: i64) -> Option<usize> {
    // BEGIN RDKIT CPP FUNCTION RDKit::ROMol::getConformer (ROMol.cpp)
    // RDKit✔️✔️:   if (d_confs.size() == 0) {
    // RDKit✔️✔️:     throw ConformerException("No conformations available on the molecule");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (id < 0) {
    // RDKit✔️✔️:     return *(d_confs.front());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto cid = (unsigned int)id;
    // RDKit✔️✔️:   for (auto conf : d_confs) {
    // RDKit✔️✔️:     if (conf->getId() == cid) { return *conf; }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION RDKit::ROMol::getConformer
    if conformers.is_empty() {
        return None;
    }
    if id < 0 {
        return Some(0);
    }
    let requested = usize::try_from(id).ok()?;
    conformers.iter().position(|conformer| conformer.id() == requested)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn negative_ids_select_first_and_sparse_ids_use_stored_identity() {
        let conformers = [Conformer3D::new(7, vec![], true), Conformer3D::new(17, vec![], true)];
        assert_eq!(resolve_3d_conformer_index(&conformers, -9), Some(0));
        assert_eq!(resolve_3d_conformer_index(&conformers, 17), Some(1));
        assert_eq!(resolve_3d_conformer_index(&conformers, 1), None);
        assert_eq!(resolve_3d_conformer_index(&[], -1), None);
    }
}
