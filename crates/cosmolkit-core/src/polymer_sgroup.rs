//! Source-shaped topology-only polymer SGroup algorithms.

use cosmolkit_model::{AtomId, BondId, SGroupConnection, SubstanceGroup};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PolymerSGroupError {
    #[error("no atoms in polymer sgroup")]
    EmptyAtoms,
    #[error("polymer SGroup atom {atom} is outside the topology with {atom_count} atoms")]
    AtomOutOfRange { atom: AtomId, atom_count: usize },
    #[error("polymer SGroup bond {bond} is outside the topology with {bond_count} bonds")]
    BondOutOfRange { bond: BondId, bond_count: usize },
}

/// Append inferred head and tail crossing bonds in source adjacency order.
///
/// `neighbors_of` exposes only the ordered incident atom and bond IDs for one
/// atom. This keeps both concrete and query lowerers on the same topology-only
/// algorithm without cloning either graph representation.
pub fn setup_unmarked_polymer_sgroup<F, I>(
    atom_count: usize,
    group_atoms: &[AtomId],
    mut neighbors_of: F,
    head_crossings: &mut Vec<BondId>,
    tail_crossings: &mut Vec<BondId>,
) -> Result<(), PolymerSGroupError>
where
    F: FnMut(AtomId) -> I,
    I: IntoIterator<Item = (AtomId, BondId)>,
{
    // RDKit❗✔️: Source anchor copied verbatim from CXSmilesOps.cpp::setupUnmarkedPolymerSGroup.
    /*
    void setupUnmarkedPolymerSGroup(RWMol &mol, SubstanceGroup &sgroup,
                                    std::vector<unsigned int> &headCrossings,
                                    std::vector<unsigned int> &tailCrossings) {
      const auto &atoms = sgroup.getAtoms();
      if (atoms.empty()) {
        throw SmilesParseException("no atoms in polymer sgroup");
      }
      const auto firstAtom = mol.getAtomWithIdx(atoms.front());
      for (auto nbr : boost::make_iterator_range(mol.getAtomNeighbors(firstAtom))) {
        const auto nbrAtom = mol[nbr];
        if (std::find(atoms.begin(), atoms.end(), nbrAtom->getIdx()) ==
            atoms.end()) {
          // in most cases we just add this to the set of headCrossings.
          // The exception occurs when there's only one atom in the SGroup and
          //  we already have a headCrossing, in which case we may put this one
          //  as a tailCrossing
          if (atoms.size() > 1 || headCrossings.empty()) {
            headCrossings.push_back(
                mol.getBondBetweenAtoms(firstAtom->getIdx(), nbrAtom->getIdx())
                    ->getIdx());
          } else if (atoms.size() == 1) {
            if (tailCrossings.empty()) {
              tailCrossings.push_back(
                  mol.getBondBetweenAtoms(firstAtom->getIdx(), nbrAtom->getIdx())
                      ->getIdx());
            } else {
              BOOST_LOG(rdWarningLog)
                  << " single atom polymer Sgroup has more than two bonds to "
                     "external atoms. Ignoring all bonds after the first two."
                  << std::endl;
            }
          }
        }
      }
      if (atoms.size() > 1) {
        const auto lastAtom = mol.getAtomWithIdx(atoms.back());
        for (auto nbr :
             boost::make_iterator_range(mol.getAtomNeighbors(lastAtom))) {
          const auto nbrAtom = mol[nbr];
          if (std::find(atoms.begin(), atoms.end(), nbrAtom->getIdx()) ==
              atoms.end()) {
            tailCrossings.push_back(
                mol.getBondBetweenAtoms(lastAtom->getIdx(), nbrAtom->getIdx())
                    ->getIdx());
          }
        }
      }
    }
    */
    if group_atoms.is_empty() {
        return Err(PolymerSGroupError::EmptyAtoms);
    }
    let first_atom = group_atoms[0];
    let last_atom = *group_atoms.last().expect("nonempty polymer group");
    for atom in [first_atom, last_atom] {
        if atom.index() >= atom_count {
            return Err(PolymerSGroupError::AtomOutOfRange { atom, atom_count });
        }
    }

    for (neighbor_atom, bond) in neighbors_of(first_atom) {
        if group_atoms.contains(&neighbor_atom) {
            continue;
        }
        if group_atoms.len() > 1 || head_crossings.is_empty() {
            head_crossings.push(bond);
        } else if tail_crossings.is_empty() {
            tail_crossings.push(bond);
        }
    }

    if group_atoms.len() > 1 {
        for (neighbor_atom, bond) in neighbors_of(last_atom) {
            if !group_atoms.contains(&neighbor_atom) {
                tail_crossings.push(bond);
            }
        }
    }

    Ok(())
}

/// Normalize polymer CONNECT state and install ordered crossing references.
///
/// Explicit crossings are already translated to detached `BondId` values by
/// their notation owner. When both sides are empty, this helper infers them
/// through the same ordered adjacency view used by concrete and query graphs.
pub fn finalize_polymer_sgroup<F, I>(
    group: &mut SubstanceGroup,
    source_connect: Option<&str>,
    source_head_crossings: &[BondId],
    source_tail_crossings: &[BondId],
    atom_count: usize,
    bond_count: usize,
    neighbors_of: F,
) -> Result<(), PolymerSGroupError>
where
    F: FnMut(AtomId) -> I,
    I: IntoIterator<Item = (AtomId, BondId)>,
{
    // RDKit❗✔️: Source anchor copied verbatim from CXSmilesOps.cpp::finalizePolymerSGroup.
    /*
    // deal with setting up the crossing bonds, etc.
    void finalizePolymerSGroup(RWMol &mol, SubstanceGroup &sgroup) {
      bool isFlipped = false;
      std::string connect = "EU";
      if (sgroup.getPropIfPresent("CONNECT", connect)) {
        if (connect.find(",f") != std::string::npos) {
          isFlipped = true;
          boost::replace_all(connect, ",f", "");
        }
      }
      if (connect == "hh") {
        connect = "HH";
      } else if (connect == "ht") {
        connect = "HT";
      } else if (connect == "eu") {
        connect = "EU";
      } else {
        BOOST_LOG(rdWarningLog) << "unrecognized CXSMILES CONNECT value: '"
                                << connect << "'. Assuming 'eu'" << std::endl;
        connect = "EU";
      }
      sgroup.setProp("CONNECT", connect);

      std::vector<unsigned int> headCrossings;
      std::vector<unsigned int> tailCrossings;
      sgroup.getPropIfPresent(_headCrossings, headCrossings);
      sgroup.clearProp(_headCrossings);
      sgroup.getPropIfPresent(_tailCrossings, tailCrossings);
      sgroup.clearProp(_tailCrossings);
      if (headCrossings.empty() && tailCrossings.empty()) {
        setupUnmarkedPolymerSGroup(mol, sgroup, headCrossings, tailCrossings);
      }
      if (headCrossings.empty() && tailCrossings.empty()) {
        // we tried... nothing more we can do
        return;
      }

      for (auto &bondIdx : headCrossings) {
        sgroup.addBondWithIdx(bondIdx);
      }
      sgroup.setProp("XBHEAD", headCrossings);

      for (auto &bondIdx : tailCrossings) {
        sgroup.addBondWithIdx(bondIdx);
      }

      // now we can setup XBCORR
      std::vector<unsigned int> xbcorr;
      for (unsigned int i = 0;
           i < std::min(headCrossings.size(), tailCrossings.size()); ++i) {
        unsigned headIdx = headCrossings[i];
        unsigned tailIdx = tailCrossings[i];
        if (isFlipped) {
          tailIdx = tailCrossings[tailCrossings.size() - i - 1];
        }
        xbcorr.push_back(headIdx);
        xbcorr.push_back(tailIdx);
      }
      sgroup.setProp("XBCORR", xbcorr);
    }
    */
    let mut connect = source_connect.unwrap_or("EU").to_owned();
    let flipped = connect.contains(",f");
    if flipped {
        connect = connect.replace(",f", "");
    }
    let (connection, normalized_connect) = match connect.as_str() {
        "hh" => (SGroupConnection::HeadToHead, "HH"),
        "ht" => (SGroupConnection::HeadToTail, "HT"),
        "eu" => (SGroupConnection::Either, "EU"),
        _ => (SGroupConnection::Either, "EU"),
    };
    group.set_connection(connection);
    group.set_prop("CONNECT", normalized_connect);

    let mut head = source_head_crossings.to_vec();
    let mut tail = source_tail_crossings.to_vec();
    if head.is_empty() && tail.is_empty() {
        setup_unmarked_polymer_sgroup(
            atom_count,
            group.atoms(),
            neighbors_of,
            &mut head,
            &mut tail,
        )?;
    }
    if head.is_empty() && tail.is_empty() {
        return Ok(());
    }

    for &bond in head.iter().chain(&tail) {
        if bond.index() >= bond_count {
            return Err(PolymerSGroupError::BondOutOfRange { bond, bond_count });
        }
    }
    for bond in head.iter().chain(&tail) {
        group.push_bond(*bond);
    }
    for bond in &head {
        group.push_head_crossing_bond(*bond);
    }
    for index in 0..head.len().min(tail.len()) {
        group.push_crossing_bond_correspondence(head[index]);
        let tail_index = if flipped {
            tail.len() - index - 1
        } else {
            index
        };
        group.push_crossing_bond_correspondence(tail[tail_index]);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{PolymerSGroupError, finalize_polymer_sgroup, setup_unmarked_polymer_sgroup};
    use cosmolkit_model::{
        AtomId, BondId, SGroupConnection, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind,
    };

    fn polymer_group(atoms: &[AtomId]) -> SubstanceGroup {
        SubstanceGroup::new(
            SubstanceGroupId::new(0),
            SubstanceGroupKind::StructuralRepeatUnit,
        )
        .with_atoms(atoms.to_vec())
    }

    #[test]
    fn query_sgroups_unmarked_polymer_rejects_empty_group() {
        let mut head = Vec::new();
        let mut tail = Vec::new();
        let error =
            setup_unmarked_polymer_sgroup(0, &[], |_| std::iter::empty(), &mut head, &mut tail)
                .expect_err("source rejects an empty polymer group");

        assert_eq!(error, PolymerSGroupError::EmptyAtoms);
        assert!(head.is_empty());
        assert!(tail.is_empty());
    }

    #[test]
    fn query_sgroups_unmarked_polymer_one_atom_keeps_first_two_external_bonds() {
        let adjacency = [
            vec![
                (AtomId::new(1), BondId::new(4)),
                (AtomId::new(2), BondId::new(2)),
                (AtomId::new(3), BondId::new(7)),
            ],
            Vec::new(),
            Vec::new(),
            Vec::new(),
        ];
        let mut head = Vec::new();
        let mut tail = Vec::new();

        setup_unmarked_polymer_sgroup(
            adjacency.len(),
            &[AtomId::new(0)],
            |atom| adjacency[atom.index()].iter().copied(),
            &mut head,
            &mut tail,
        )
        .expect("one atom uses its first two external bonds");

        assert_eq!(head, [BondId::new(4)]);
        assert_eq!(tail, [BondId::new(2)]);
    }

    #[test]
    fn query_sgroups_unmarked_polymer_multi_atom_uses_first_and_last_members() {
        let adjacency = [
            Vec::new(),
            vec![
                (AtomId::new(0), BondId::new(0)),
                (AtomId::new(2), BondId::new(1)),
                (AtomId::new(4), BondId::new(2)),
            ],
            Vec::new(),
            vec![
                (AtomId::new(2), BondId::new(3)),
                (AtomId::new(5), BondId::new(4)),
                (AtomId::new(4), BondId::new(5)),
            ],
            Vec::new(),
            Vec::new(),
        ];
        let mut head = Vec::new();
        let mut tail = Vec::new();

        setup_unmarked_polymer_sgroup(
            adjacency.len(),
            &[AtomId::new(1), AtomId::new(2), AtomId::new(3)],
            |atom| adjacency[atom.index()].iter().copied(),
            &mut head,
            &mut tail,
        )
        .expect("multi-atom group uses ordered first and last members");

        assert_eq!(head, [BondId::new(0), BondId::new(2)]);
        assert_eq!(tail, [BondId::new(4), BondId::new(5)]);
    }

    #[test]
    fn query_sgroups_unmarked_polymer_preserves_duplicate_group_membership_effects() {
        let adjacency = [vec![(AtomId::new(1), BondId::new(6))], Vec::new()];
        let mut head = Vec::new();
        let mut tail = Vec::new();

        setup_unmarked_polymer_sgroup(
            adjacency.len(),
            &[AtomId::new(0), AtomId::new(0)],
            |atom| adjacency[atom.index()].iter().copied(),
            &mut head,
            &mut tail,
        )
        .expect("the source preserves repeated group member rows");

        assert_eq!(head, [BondId::new(6)]);
        assert_eq!(tail, [BondId::new(6)]);
    }

    #[test]
    fn query_sgroups_polymer_connect_accepts_only_source_lowercase_values() {
        for (source, expected_connection, expected_prop) in [
            ("hh", SGroupConnection::HeadToHead, "HH"),
            ("ht", SGroupConnection::HeadToTail, "HT"),
            ("eu", SGroupConnection::Either, "EU"),
            ("HH", SGroupConnection::Either, "EU"),
            ("unknown", SGroupConnection::Either, "EU"),
        ] {
            let mut group = polymer_group(&[AtomId::new(0)]);
            finalize_polymer_sgroup(
                &mut group,
                Some(source),
                &[BondId::new(0)],
                &[],
                1,
                1,
                |_| std::iter::empty::<(AtomId, BondId)>(),
            )
            .expect("source CONNECT fallback is non-failing");

            assert_eq!(group.connection(), Some(&expected_connection), "{source}");
            assert_eq!(
                group.props().get("CONNECT").map(String::as_str),
                Some(expected_prop)
            );
        }
    }

    #[test]
    fn query_sgroups_polymer_flip_removes_every_marker_and_reverses_tail_pairs() {
        let mut group = polymer_group(&[AtomId::new(0), AtomId::new(1)]);
        finalize_polymer_sgroup(
            &mut group,
            Some("hh,f,f"),
            &[BondId::new(1), BondId::new(1)],
            &[BondId::new(2), BondId::new(2), BondId::new(3)],
            2,
            4,
            |_| std::iter::empty::<(AtomId, BondId)>(),
        )
        .expect("repeated flip markers use source replace-all");

        assert_eq!(group.connection(), Some(&SGroupConnection::HeadToHead));
        assert_eq!(group.props().get("CONNECT").map(String::as_str), Some("HH"));
        assert_eq!(
            group.bonds(),
            &[
                BondId::new(1),
                BondId::new(1),
                BondId::new(2),
                BondId::new(2),
                BondId::new(3),
            ]
        );
        assert_eq!(
            group.head_crossing_bonds(),
            &[BondId::new(1), BondId::new(1)]
        );
        assert_eq!(
            group.crossing_bond_correspondence(),
            &[
                BondId::new(1),
                BondId::new(3),
                BondId::new(1),
                BondId::new(2),
            ]
        );
    }

    #[test]
    fn query_sgroups_polymer_correspondence_uses_minimum_unflipped_length() {
        let mut group = polymer_group(&[AtomId::new(0), AtomId::new(1)]);
        finalize_polymer_sgroup(
            &mut group,
            Some("ht"),
            &[BondId::new(0), BondId::new(1), BondId::new(2)],
            &[BondId::new(3)],
            2,
            4,
            |_| std::iter::empty::<(AtomId, BondId)>(),
        )
        .expect("source pairs only as many crossings as both sides provide");

        assert_eq!(
            group.crossing_bond_correspondence(),
            &[BondId::new(0), BondId::new(3)]
        );
    }

    #[test]
    fn query_sgroups_polymer_absent_crossings_keep_connect_without_inventing_bonds() {
        let mut group = polymer_group(&[AtomId::new(0)]);
        finalize_polymer_sgroup(&mut group, None, &[], &[], 1, 0, |_| {
            std::iter::empty::<(AtomId, BondId)>()
        })
        .expect("no inferred crossings is a source early return");

        assert_eq!(group.connection(), Some(&SGroupConnection::Either));
        assert_eq!(group.props().get("CONNECT").map(String::as_str), Some("EU"));
        assert!(group.bonds().is_empty());
        assert!(group.head_crossing_bonds().is_empty());
        assert!(group.crossing_bond_correspondence().is_empty());
    }

    #[test]
    fn query_sgroups_polymer_infers_only_when_both_crossing_arrays_are_empty() {
        let adjacency = [vec![(AtomId::new(1), BondId::new(0))], Vec::new()];
        let mut group = polymer_group(&[AtomId::new(0)]);
        finalize_polymer_sgroup(
            &mut group,
            Some("eu"),
            &[],
            &[BondId::new(1)],
            adjacency.len(),
            2,
            |atom| adjacency[atom.index()].iter().copied(),
        )
        .expect("one explicit side prevents inference of the other");

        assert!(group.head_crossing_bonds().is_empty());
        assert_eq!(group.bonds(), &[BondId::new(1)]);
        assert!(group.crossing_bond_correspondence().is_empty());
    }
}
