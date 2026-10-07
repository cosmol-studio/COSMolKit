//! Source-shaped topology-only polymer SGroup algorithms.

use cosmolkit_model::{AtomId, BondId, SGroupConnection, SubstanceGroup};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PolymerSGroupError {
    #[error("polymer SGroup property operation failed: {0}")]
    Property(#[from] cosmolkit_model::MoleculePropertyError),
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
    // BEGIN COMPLETE PINNED SF182
    // RDKit✔️✔️: void setupUnmarkedPolymerSGroup(RWMol &mol, SubstanceGroup &sgroup,
    // RDKit✔️✔️:                                 std::vector<unsigned int> &headCrossings,
    // RDKit✔️✔️:                                 std::vector<unsigned int> &tailCrossings) {
    // RDKit✔️✔️:   const auto &atoms = sgroup.getAtoms();
    // RDKit✔️✔️:   if (atoms.empty()) {
    // RDKit✔️✔️:     throw SmilesParseException("no atoms in polymer sgroup");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const auto firstAtom = mol.getAtomWithIdx(atoms.front());
    // RDKit✔️✔️:   for (auto nbr : boost::make_iterator_range(mol.getAtomNeighbors(firstAtom))) {
    // RDKit✔️✔️:     const auto nbrAtom = mol[nbr];
    // RDKit✔️✔️:     if (std::find(atoms.begin(), atoms.end(), nbrAtom->getIdx()) ==
    // RDKit✔️✔️:         atoms.end()) {
    // RDKit✔️✔️:       // in most cases we just add this to the set of headCrossings.
    // RDKit✔️✔️:       // The exception occurs when there's only one atom in the SGroup and
    // RDKit✔️✔️:       //  we already have a headCrossing, in which case we may put this one
    // RDKit✔️✔️:       //  as a tailCrossing
    // RDKit✔️✔️:       if (atoms.size() > 1 || headCrossings.empty()) {
    // RDKit✔️✔️:         headCrossings.push_back(
    // RDKit✔️✔️:             mol.getBondBetweenAtoms(firstAtom->getIdx(), nbrAtom->getIdx())
    // RDKit✔️✔️:                 ->getIdx());
    // RDKit✔️✔️:       } else if (atoms.size() == 1) {
    // RDKit✔️✔️:         if (tailCrossings.empty()) {
    // RDKit✔️✔️:           tailCrossings.push_back(
    // RDKit✔️✔️:               mol.getBondBetweenAtoms(firstAtom->getIdx(), nbrAtom->getIdx())
    // RDKit✔️✔️:                   ->getIdx());
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:               << " single atom polymer Sgroup has more than two bonds to "
    // RDKit✔️✔️:                  "external atoms. Ignoring all bonds after the first two."
    // RDKit✔️✔️:               << std::endl;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (atoms.size() > 1) {
    // RDKit✔️✔️:     const auto lastAtom = mol.getAtomWithIdx(atoms.back());
    // RDKit✔️✔️:     for (auto nbr :
    // RDKit✔️✔️:          boost::make_iterator_range(mol.getAtomNeighbors(lastAtom))) {
    // RDKit✔️✔️:       const auto nbrAtom = mol[nbr];
    // RDKit✔️✔️:       if (std::find(atoms.begin(), atoms.end(), nbrAtom->getIdx()) ==
    // RDKit✔️✔️:           atoms.end()) {
    // RDKit✔️✔️:         tailCrossings.push_back(
    // RDKit✔️✔️:             mol.getBondBetweenAtoms(lastAtom->getIdx(), nbrAtom->getIdx())
    // RDKit✔️✔️:                 ->getIdx());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END COMPLETE PINNED SF182
    // Complete reached ROMol::getAtomWithIdx(unsigned int)
    // RDKit✔️✔️: Atom *ROMol::getAtomWithIdx(unsigned int idx) {
    // RDKit✔️✔️:   URANGE_CHECK(idx, getNumAtoms());
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto vd = boost::vertex(idx, d_graph);
    // RDKit✔️✔️:   auto res = d_graph[vd];
    // RDKit✔️✔️:   POSTCONDITION(res, "");
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Complete reached ROMol::getAtomWithIdx(unsigned int) const
    // RDKit✔️✔️: const Atom *ROMol::getAtomWithIdx(unsigned int idx) const {
    // RDKit✔️✔️:   URANGE_CHECK(idx, getNumAtoms());
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto vd = boost::vertex(idx, d_graph);
    // RDKit✔️✔️:   const auto res = d_graph[vd];
    // RDKit✔️✔️:
    // RDKit✔️✔️:   POSTCONDITION(res, "");
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Complete reached ROMol::getBondBetweenAtoms mutable forwarding
    // RDKit✔️✔️: Bond *ROMol::getBondBetweenAtoms(unsigned int idx1, unsigned int idx2) {
    // RDKit✔️✔️:   return const_cast<Bond *>(
    // RDKit✔️✔️:       static_cast<const ROMol *>(this)->getBondBetweenAtoms(
    // RDKit✔️✔️:           idx1, idx2));  // avoid code duplication
    // RDKit✔️✔️: }
    // Complete reached ROMol::getBondBetweenAtoms const
    // RDKit✔️✔️: const Bond *ROMol::getBondBetweenAtoms(unsigned int idx1,
    // RDKit✔️✔️:                                        unsigned int idx2) const {
    // RDKit✔️✔️:   URANGE_CHECK(idx1, getNumAtoms());
    // RDKit✔️✔️:   URANGE_CHECK(idx2, getNumAtoms());
    // RDKit✔️✔️:   const Bond *res = nullptr;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto [edge, found] = boost::edge(boost::vertex(idx1, d_graph),
    // RDKit✔️✔️:                                    boost::vertex(idx2, d_graph), d_graph);
    // RDKit✔️✔️:   if (found) {
    // RDKit✔️✔️:     res = d_graph[edge];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Complete reached ROMol::getAtomNeighbors
    // RDKit✔️✔️: ROMol::ADJ_ITER_PAIR ROMol::getAtomNeighbors(Atom const *at) const {
    // RDKit✔️✔️:   PRECONDITION(at, "no atom");
    // RDKit✔️✔️:   PRECONDITION(&at->getOwningMol() == this,
    // RDKit✔️✔️:                "atom not associated with this molecule");
    // RDKit✔️✔️:   return boost::adjacent_vertices(at->getIdx(), d_graph);
    // RDKit✔️✔️: };
    // Complete reached SubstanceGroup.h
    // RDKit✔️✔️: const std::vector<unsigned int> &getAtoms() const { return d_atoms; }
    // Complete reached Atom.h
    // RDKit✔️✔️: unsigned int getIdx() const { return d_index; }
    // Complete reached Bond.h
    // RDKit✔️✔️: unsigned int getIdx() const { return d_index; }
    // Behavior: ordered borrowed members/adjacency preserve source row order,
    // duplicates, existing head/tail arrays and warning branch. Bounds checks
    // occur only when source getAtomWithIdx reaches the corresponding endpoint.
    // Complexity: two degree-bounded adjacency passes, linear member lookup,
    // amortized vector append and O(1) extra state, no group/graph clones.
    // The already validated adjacency tuple supplies the source bond ID without
    // an additional getBondBetweenAtoms lookup or a second topology index.
    if group_atoms.is_empty() {
        return Err(PolymerSGroupError::EmptyAtoms);
    }
    let first_atom = group_atoms[0];
    if first_atom.index() >= atom_count {
        return Err(PolymerSGroupError::AtomOutOfRange {
            atom: first_atom,
            atom_count,
        });
    }

    for (neighbor_atom, bond) in neighbors_of(first_atom) {
        if group_atoms.contains(&neighbor_atom) {
            continue;
        }
        if group_atoms.len() > 1 || head_crossings.is_empty() {
            head_crossings.push(bond);
        } else if tail_crossings.is_empty() {
            tail_crossings.push(bond);
        } else {
            eprintln!(
                " single atom polymer Sgroup has more than two bonds to external atoms. Ignoring all bonds after the first two."
            );
        }
    }

    if group_atoms.len() > 1 {
        let last_atom = *group_atoms.last().expect("nonempty polymer group");
        if last_atom.index() >= atom_count {
            return Err(PolymerSGroupError::AtomOutOfRange {
                atom: last_atom,
                atom_count,
            });
        }
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
    source_connect: Option<&[u8]>,
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
    // BEGIN COMPLETE PINNED SF183
    // RDKit✔️✔️: void finalizePolymerSGroup(RWMol &mol, SubstanceGroup &sgroup) {
    // RDKit✔️✔️:   bool isFlipped = false;
    // RDKit✔️✔️:   std::string connect = "EU";
    // RDKit✔️✔️:   if (sgroup.getPropIfPresent("CONNECT", connect)) {
    // RDKit✔️✔️:     if (connect.find(",f") != std::string::npos) {
    // RDKit✔️✔️:       isFlipped = true;
    // RDKit✔️✔️:       boost::replace_all(connect, ",f", "");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (connect == "hh") {
    // RDKit✔️✔️:     connect = "HH";
    // RDKit✔️✔️:   } else if (connect == "ht") {
    // RDKit✔️✔️:     connect = "HT";
    // RDKit✔️✔️:   } else if (connect == "eu") {
    // RDKit✔️✔️:     connect = "EU";
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog) << "unrecognized CXSMILES CONNECT value: '"
    // RDKit✔️✔️:                             << connect << "'. Assuming 'eu'" << std::endl;
    // RDKit✔️✔️:     connect = "EU";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   sgroup.setProp("CONNECT", connect);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<unsigned int> headCrossings;
    // RDKit✔️✔️:   std::vector<unsigned int> tailCrossings;
    // RDKit✔️✔️:   sgroup.getPropIfPresent(_headCrossings, headCrossings);
    // RDKit✔️✔️:   sgroup.clearProp(_headCrossings);
    // RDKit✔️✔️:   sgroup.getPropIfPresent(_tailCrossings, tailCrossings);
    // RDKit✔️✔️:   sgroup.clearProp(_tailCrossings);
    // RDKit✔️✔️:   if (headCrossings.empty() && tailCrossings.empty()) {
    // RDKit✔️✔️:     setupUnmarkedPolymerSGroup(mol, sgroup, headCrossings, tailCrossings);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (headCrossings.empty() && tailCrossings.empty()) {
    // RDKit✔️✔️:     // we tried... nothing more we can do
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto &bondIdx : headCrossings) {
    // RDKit✔️✔️:     sgroup.addBondWithIdx(bondIdx);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   sgroup.setProp("XBHEAD", headCrossings);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto &bondIdx : tailCrossings) {
    // RDKit✔️✔️:     sgroup.addBondWithIdx(bondIdx);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // now we can setup XBCORR
    // RDKit✔️✔️:   std::vector<unsigned int> xbcorr;
    // RDKit✔️✔️:   for (unsigned int i = 0;
    // RDKit✔️✔️:        i < std::min(headCrossings.size(), tailCrossings.size()); ++i) {
    // RDKit✔️✔️:     unsigned headIdx = headCrossings[i];
    // RDKit✔️✔️:     unsigned tailIdx = tailCrossings[i];
    // RDKit✔️✔️:     if (isFlipped) {
    // RDKit✔️✔️:       tailIdx = tailCrossings[tailCrossings.size() - i - 1];
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     xbcorr.push_back(headIdx);
    // RDKit✔️✔️:     xbcorr.push_back(tailIdx);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   sgroup.setProp("XBCORR", xbcorr);
    // RDKit✔️✔️: }
    // END COMPLETE PINNED SF183
    // Complete reached SubstanceGroup::addBondWithIdx
    // RDKit✔️✔️: void SubstanceGroup::addBondWithIdx(unsigned int idx) {
    // RDKit✔️✔️:   PRECONDITION(dp_mol, "bad mol");
    // RDKit✔️✔️:   if (idx >= dp_mol->getNumBonds()) {
    // RDKit✔️✔️:     throw ValueErrorException("Bond index out of range");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   d_bonds.push_back(idx);
    // RDKit✔️✔️: }
    // Complete reached RDProps helper
    // RDKit✔️✔️: void setProp(const std::string_view key, T val, bool computed = false) const {
    // RDKit✔️✔️:     if(key.empty()) {
    // RDKit✔️✔️:       throw ValueErrorException("Cannot set property with empty key");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (computed) {
    // RDKit✔️✔️:       STR_VECT compLst;
    // RDKit✔️✔️:       getPropIfPresent(RDKit::detail::computedPropName, compLst);
    // RDKit✔️✔️:       if (std::find(compLst.begin(), compLst.end(), key) == compLst.end()) {
    // RDKit✔️✔️:         compLst.emplace_back(key);
    // RDKit✔️✔️:         d_props.setVal(RDKit::detail::computedPropName, compLst);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     d_props.setVal(key, val);
    // RDKit✔️✔️:   }
    // Complete reached RDProps helper
    // RDKit✔️✔️: bool getPropIfPresent(const std::string_view key, T &res) const {
    // RDKit✔️✔️:     return d_props.getValIfPresent(key, res);
    // RDKit✔️✔️:   }
    // Complete reached RDProps helper
    // RDKit✔️✔️: void clearProp(const std::string_view key) const {
    // RDKit✔️✔️:     STR_VECT compLst;
    // RDKit✔️✔️:     if (getPropIfPresent(RDKit::detail::computedPropName, compLst)) {
    // RDKit✔️✔️:       auto svi = std::find(compLst.begin(), compLst.end(), key);
    // RDKit✔️✔️:       if (svi != compLst.end()) {
    // RDKit✔️✔️:         compLst.erase(svi);
    // RDKit✔️✔️:         d_props.setVal(RDKit::detail::computedPropName, compLst);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     d_props.clearVal(key);
    // RDKit✔️✔️:   }
    // Complete reached pinned Boost1.81 replace.hpp
    // Boost✔️✔️: inline void replace_all(
    // Boost✔️✔️:             SequenceT& Input,
    // Boost✔️✔️:             const Range1T& Search,
    // Boost✔️✔️:             const Range2T& Format )
    // Boost✔️✔️:         {
    // Boost✔️✔️:             ::boost::algorithm::find_format_all(
    // Boost✔️✔️:                 Input,
    // Boost✔️✔️:                 ::boost::algorithm::first_finder(Search),
    // Boost✔️✔️:                 ::boost::algorithm::const_formatter(Format) );
    // Boost✔️✔️:         }
    // Complete reached pinned Boost1.81 find_format.hpp
    // Boost✔️✔️: inline void find_format_all(
    // Boost✔️✔️:             SequenceT& Input,
    // Boost✔️✔️:             FinderT Finder,
    // Boost✔️✔️:             FormatterT Formatter )
    // Boost✔️✔️:         {
    // Boost✔️✔️:             // Concept check
    // Boost✔️✔️:             BOOST_CONCEPT_ASSERT((
    // Boost✔️✔️:                 FinderConcept<
    // Boost✔️✔️:                     FinderT,
    // Boost✔️✔️:                     BOOST_STRING_TYPENAME range_const_iterator<SequenceT>::type>
    // Boost✔️✔️:                 ));
    // Boost✔️✔️:             BOOST_CONCEPT_ASSERT((
    // Boost✔️✔️:                 FormatterConcept<
    // Boost✔️✔️:                     FormatterT,
    // Boost✔️✔️:                     FinderT,BOOST_STRING_TYPENAME range_const_iterator<SequenceT>::type>
    // Boost✔️✔️:                 ));
    // Boost✔️✔️:
    // Boost✔️✔️:             detail::find_format_all_impl(
    // Boost✔️✔️:                 Input,
    // Boost✔️✔️:                 Finder,
    // Boost✔️✔️:                 Formatter,
    // Boost✔️✔️:                 Finder(::boost::begin(Input), ::boost::end(Input)));
    // Boost✔️✔️:
    // Boost✔️✔️:         }
    // Complete reached pinned Boost1.81 detail/find_format_all.hpp
    // Boost✔️✔️: inline void find_format_all_impl2(
    // Boost✔️✔️:                 InputT& Input,
    // Boost✔️✔️:                 FinderT Finder,
    // Boost✔️✔️:                 FormatterT Formatter,
    // Boost✔️✔️:                 FindResultT FindResult,
    // Boost✔️✔️:                 FormatResultT FormatResult)
    // Boost✔️✔️:             {
    // Boost✔️✔️:                 typedef BOOST_STRING_TYPENAME
    // Boost✔️✔️:                     range_iterator<InputT>::type input_iterator_type;
    // Boost✔️✔️:                 typedef find_format_store<
    // Boost✔️✔️:                         input_iterator_type,
    // Boost✔️✔️:                         FormatterT,
    // Boost✔️✔️:                         FormatResultT > store_type;
    // Boost✔️✔️:
    // Boost✔️✔️:                 // Create store for the find result
    // Boost✔️✔️:                 store_type M( FindResult, FormatResult, Formatter );
    // Boost✔️✔️:
    // Boost✔️✔️:                 // Instantiate replacement storage
    // Boost✔️✔️:                 std::deque<
    // Boost✔️✔️:                     BOOST_STRING_TYPENAME range_value<InputT>::type> Storage;
    // Boost✔️✔️:
    // Boost✔️✔️:                 // Initialize replacement iterators
    // Boost✔️✔️:                 input_iterator_type InsertIt=::boost::begin(Input);
    // Boost✔️✔️:                 input_iterator_type SearchIt=::boost::begin(Input);
    // Boost✔️✔️:
    // Boost✔️✔️:                 while( M )
    // Boost✔️✔️:                 {
    // Boost✔️✔️:                     // process the segment
    // Boost✔️✔️:                     InsertIt=process_segment(
    // Boost✔️✔️:                         Storage,
    // Boost✔️✔️:                         Input,
    // Boost✔️✔️:                         InsertIt,
    // Boost✔️✔️:                         SearchIt,
    // Boost✔️✔️:                         M.begin() );
    // Boost✔️✔️:
    // Boost✔️✔️:                     // Adjust search iterator
    // Boost✔️✔️:                     SearchIt=M.end();
    // Boost✔️✔️:
    // Boost✔️✔️:                     // Copy formatted replace to the storage
    // Boost✔️✔️:                     ::boost::algorithm::detail::copy_to_storage( Storage, M.format_result() );
    // Boost✔️✔️:
    // Boost✔️✔️:                     // Find range for a next match
    // Boost✔️✔️:                     M=Finder( SearchIt, ::boost::end(Input) );
    // Boost✔️✔️:                 }
    // Boost✔️✔️:
    // Boost✔️✔️:                 // process the last segment
    // Boost✔️✔️:                 InsertIt=::boost::algorithm::detail::process_segment(
    // Boost✔️✔️:                     Storage,
    // Boost✔️✔️:                     Input,
    // Boost✔️✔️:                     InsertIt,
    // Boost✔️✔️:                     SearchIt,
    // Boost✔️✔️:                     ::boost::end(Input) );
    // Boost✔️✔️:
    // Boost✔️✔️:                 if ( Storage.empty() )
    // Boost✔️✔️:                 {
    // Boost✔️✔️:                     // Truncate input
    // Boost✔️✔️:                     ::boost::algorithm::detail::erase( Input, InsertIt, ::boost::end(Input) );
    // Boost✔️✔️:                 }
    // Boost✔️✔️:                 else
    // Boost✔️✔️:                 {
    // Boost✔️✔️:                     // Copy remaining data to the end of input
    // Boost✔️✔️:                     ::boost::algorithm::detail::insert( Input, ::boost::end(Input), Storage.begin(), Storage.end() );
    // Boost✔️✔️:                 }
    // Boost✔️✔️:             }
    // Complete reached pinned Boost1.81 detail/find_format_all.hpp
    // Boost✔️✔️: inline void find_format_all_impl(
    // Boost✔️✔️:                 InputT& Input,
    // Boost✔️✔️:                 FinderT Finder,
    // Boost✔️✔️:                 FormatterT Formatter,
    // Boost✔️✔️:                 FindResultT FindResult)
    // Boost✔️✔️:             {
    // Boost✔️✔️:                 if( ::boost::algorithm::detail::check_find_result(Input, FindResult) ) {
    // Boost✔️✔️:                     ::boost::algorithm::detail::find_format_all_impl2(
    // Boost✔️✔️:                         Input,
    // Boost✔️✔️:                         Finder,
    // Boost✔️✔️:                         Formatter,
    // Boost✔️✔️:                         FindResult,
    // Boost✔️✔️:                         Formatter(FindResult) );
    // Boost✔️✔️:                 }
    // Boost✔️✔️:             }
    // Complete reached pinned Boost1.81 first_finderF::operator()
    // Boost✔️✔️: operator()(
    // Boost✔️✔️:                     ForwardIteratorT Begin,
    // Boost✔️✔️:                     ForwardIteratorT End ) const
    // Boost✔️✔️:                 {
    // Boost✔️✔️:                     typedef iterator_range<ForwardIteratorT> result_type;
    // Boost✔️✔️:                     typedef ForwardIteratorT input_iterator_type;
    // Boost✔️✔️:
    // Boost✔️✔️:                     // Outer loop
    // Boost✔️✔️:                     for(input_iterator_type OuterIt=Begin;
    // Boost✔️✔️:                         OuterIt!=End;
    // Boost✔️✔️:                         ++OuterIt)
    // Boost✔️✔️:                     {
    // Boost✔️✔️:                         // Sanity check
    // Boost✔️✔️:                         if( boost::empty(m_Search) )
    // Boost✔️✔️:                             return result_type( End, End );
    // Boost✔️✔️:
    // Boost✔️✔️:                         input_iterator_type InnerIt=OuterIt;
    // Boost✔️✔️:                         search_iterator_type SubstrIt=m_Search.begin();
    // Boost✔️✔️:                         for(;
    // Boost✔️✔️:                             InnerIt!=End && SubstrIt!=m_Search.end();
    // Boost✔️✔️:                             ++InnerIt,++SubstrIt)
    // Boost✔️✔️:                         {
    // Boost✔️✔️:                             if( !( m_Comp(*InnerIt,*SubstrIt) ) )
    // Boost✔️✔️:                                 break;
    // Boost✔️✔️:                         }
    // Boost✔️✔️:
    // Boost✔️✔️:                         // Substring matching succeeded
    // Boost✔️✔️:                         if ( SubstrIt==m_Search.end() )
    // Boost✔️✔️:                             return result_type( OuterIt, InnerIt );
    // Boost✔️✔️:                     }
    // Boost✔️✔️:
    // Boost✔️✔️:                     return result_type( End, End );
    // Boost✔️✔️:                 }
    // Complete reached pinned Boost1.81 detail/formatter.hpp
    // Boost✔️✔️: const result_type& operator()(const Range2T&) const
    // Boost✔️✔️:                 {
    // Boost✔️✔️:                     return m_Format;
    // Boost✔️✔️:                 }
    // Complete reached pinned Boost1.81 detail/find_format_store.hpp
    // Boost✔️✔️: find_format_store& operator=( FindResultT FindResult )
    // Boost✔️✔️:                 {
    // Boost✔️✔️:                     iterator_range<ForwardIteratorT>::operator=(FindResult);
    // Boost✔️✔️:                     if( !this->empty() ) {
    // Boost✔️✔️:                         m_FormatResult=m_Formatter(FindResult);
    // Boost✔️✔️:                     }
    // Boost✔️✔️:
    // Boost✔️✔️:                     return *this;
    // Boost✔️✔️:                 }
    // Complete reached pinned Boost1.81 detail/find_format_store.hpp
    // Boost✔️✔️: const format_result_type& format_result()
    // Boost✔️✔️:                 {
    // Boost✔️✔️:                     return m_FormatResult;
    // Boost✔️✔️:                 }
    // Complete reached pinned Boost1.81 detail/find_format_store.hpp
    // Boost✔️✔️: bool check_find_result(InputT&, FindResultT& FindResult)
    // Boost✔️✔️:             {
    // Boost✔️✔️:                 typedef BOOST_STRING_TYPENAME
    // Boost✔️✔️:                     range_const_iterator<InputT>::type input_iterator_type;
    // Boost✔️✔️:                 iterator_range<input_iterator_type> ResultRange(FindResult);
    // Boost✔️✔️:                 return !ResultRange.empty();
    // Boost✔️✔️:             }
    // Complete reached pinned Boost1.81 detail/replace_storage.hpp
    // Boost✔️✔️: inline OutputIteratorT move_from_storage(
    // Boost✔️✔️:                 StorageT& Storage,
    // Boost✔️✔️:                 OutputIteratorT DestBegin,
    // Boost✔️✔️:                 OutputIteratorT DestEnd )
    // Boost✔️✔️:             {
    // Boost✔️✔️:                 OutputIteratorT OutputIt=DestBegin;
    // Boost✔️✔️:
    // Boost✔️✔️:                 while( !Storage.empty() && OutputIt!=DestEnd )
    // Boost✔️✔️:                 {
    // Boost✔️✔️:                     *OutputIt=Storage.front();
    // Boost✔️✔️:                     Storage.pop_front();
    // Boost✔️✔️:                     ++OutputIt;
    // Boost✔️✔️:                 }
    // Boost✔️✔️:
    // Boost✔️✔️:                 return OutputIt;
    // Boost✔️✔️:             }
    // Complete reached pinned Boost1.81 detail/replace_storage.hpp
    // Boost✔️✔️: inline void copy_to_storage(
    // Boost✔️✔️:                 StorageT& Storage,
    // Boost✔️✔️:                 const WhatT& What )
    // Boost✔️✔️:             {
    // Boost✔️✔️:                 Storage.insert( Storage.end(), ::boost::begin(What), ::boost::end(What) );
    // Boost✔️✔️:             }
    // Complete reached pinned Boost1.81 detail/replace_storage.hpp
    // Boost✔️✔️: ForwardIteratorT operator()(
    // Boost✔️✔️:                     StorageT& Storage,
    // Boost✔️✔️:                     InputT& /*Input*/,
    // Boost✔️✔️:                     ForwardIteratorT InsertIt,
    // Boost✔️✔️:                     ForwardIteratorT SegmentBegin,
    // Boost✔️✔️:                     ForwardIteratorT SegmentEnd )
    // Boost✔️✔️:                 {
    // Boost✔️✔️:                     // Copy data from the storage until the beginning of the segment
    // Boost✔️✔️:                     ForwardIteratorT It=::boost::algorithm::detail::move_from_storage( Storage, InsertIt, SegmentBegin );
    // Boost✔️✔️:
    // Boost✔️✔️:                     // 3 cases are possible :
    // Boost✔️✔️:                     //   a) Storage is empty, It==SegmentBegin
    // Boost✔️✔️:                     //   b) Storage is empty, It!=SegmentBegin
    // Boost✔️✔️:                     //   c) Storage is not empty
    // Boost✔️✔️:
    // Boost✔️✔️:                     if( Storage.empty() )
    // Boost✔️✔️:                     {
    // Boost✔️✔️:                         if( It==SegmentBegin )
    // Boost✔️✔️:                         {
    // Boost✔️✔️:                             // Case a) everything is grand, just return end of segment
    // Boost✔️✔️:                             return SegmentEnd;
    // Boost✔️✔️:                         }
    // Boost✔️✔️:                         else
    // Boost✔️✔️:                         {
    // Boost✔️✔️:                             // Case b) move the segment backwards
    // Boost✔️✔️:                             return std::copy( SegmentBegin, SegmentEnd, It );
    // Boost✔️✔️:                         }
    // Boost✔️✔️:                     }
    // Boost✔️✔️:                     else
    // Boost✔️✔️:                     {
    // Boost✔️✔️:                         // Case c) -> shift the segment to the left and keep the overlap in the storage
    // Boost✔️✔️:                         while( It!=SegmentEnd )
    // Boost✔️✔️:                         {
    // Boost✔️✔️:                             // Store value into storage
    // Boost✔️✔️:                             Storage.push_back( *It );
    // Boost✔️✔️:                             // Get the top from the storage and put it here
    // Boost✔️✔️:                             *It=Storage.front();
    // Boost✔️✔️:                             Storage.pop_front();
    // Boost✔️✔️:
    // Boost✔️✔️:                             // Advance
    // Boost✔️✔️:                             ++It;
    // Boost✔️✔️:                         }
    // Boost✔️✔️:
    // Boost✔️✔️:                         return It;
    // Boost✔️✔️:                     }
    // Boost✔️✔️:                 }
    // Complete reached pinned Boost1.81 detail/replace_storage.hpp
    // Boost✔️✔️: inline ForwardIteratorT process_segment(
    // Boost✔️✔️:                 StorageT& Storage,
    // Boost✔️✔️:                 InputT& Input,
    // Boost✔️✔️:                 ForwardIteratorT InsertIt,
    // Boost✔️✔️:                 ForwardIteratorT SegmentBegin,
    // Boost✔️✔️:                 ForwardIteratorT SegmentEnd )
    // Boost✔️✔️:             {
    // Boost✔️✔️:                 return
    // Boost✔️✔️:                     process_segment_helper<
    // Boost✔️✔️:                         has_stable_iterators<InputT>::value>()(
    // Boost✔️✔️:                                 Storage, Input, InsertIt, SegmentBegin, SegmentEnd );
    // Boost✔️✔️:             }
    // Complete reached pinned Boost1.81 detail/sequence.hpp
    // Boost✔️✔️: inline typename InputT::iterator erase(
    // Boost✔️✔️:                 InputT& Input,
    // Boost✔️✔️:                 BOOST_STRING_TYPENAME InputT::iterator From,
    // Boost✔️✔️:                 BOOST_STRING_TYPENAME InputT::iterator To )
    // Boost✔️✔️:             {
    // Boost✔️✔️:                 return Input.erase( From, To );
    // Boost✔️✔️:             }
    // Behavior: byte-exact lowercase CONNECT dispatch, every non-overlapping
    // flip marker removed, source warning/fallback and mutation order retained.
    // Detached head/tail arguments are the canonical typed transient vectors;
    // XBHEAD/XBCORR are represented by their existing typed MODEL fields.
    // Logging control/timestamps are independent unmodeled RDLog capabilities;
    // this owner preserves the source warning branch and counted payload bytes.
    // Complexity: one O(L) string copy and linear in-place byte compaction,
    // source-sized head/tail copies, O(H+T) appends and O(min(H,T)) pairing.
    // Existing value builders replace typed vectors through O(1) owned moves,
    // without cloning any graph/group state or introducing a second interface.
    use std::io::Write as _;

    let mut connect = source_connect.unwrap_or(b"EU").to_vec();
    let flipped = connect.windows(2).any(|pair| pair == b",f");
    if flipped {
        let mut read = 0;
        let mut write = 0;
        while read < connect.len() {
            if connect[read..].starts_with(b",f") {
                read += 2;
            } else {
                connect[write] = connect[read];
                read += 1;
                write += 1;
            }
        }
        connect.truncate(write);
    }
    let (connection, normalized_connect) = match connect.as_slice() {
        b"hh" => (SGroupConnection::HeadToHead, "HH"),
        b"ht" => (SGroupConnection::HeadToTail, "HT"),
        b"eu" => (SGroupConnection::Either, "EU"),
        _ => {
            // std::ostream's default exception mask does not turn warning sink
            // failures into chemistry failures. Preserve raw counted bytes;
            // the source message is not a UTF-8 or C-string conversion.
            let mut warning = std::io::stderr().lock();
            let _ = warning
                .write_all(b"unrecognized CXSMILES CONNECT value: '")
                .and_then(|()| warning.write_all(&connect))
                .and_then(|()| warning.write_all(b"'. Assuming 'eu'\n"))
                .and_then(|()| warning.flush());
            (SGroupConnection::Either, "EU")
        }
    };
    group.set_connection(connection);
    group.set_prop("CONNECT", normalized_connect)?;

    let mut head = source_head_crossings.to_vec();
    group.clear_prop("_headCrossings")?;
    let mut tail = source_tail_crossings.to_vec();
    group.clear_prop("_tailCrossings")?;
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

    for &bond in &head {
        if bond.index() >= bond_count {
            return Err(PolymerSGroupError::BondOutOfRange { bond, bond_count });
        }
        group.push_bond(bond);
    }
    // The temporary value is inaccessible under this exclusive borrow and
    // never interpreted or committed. Existing builders move all original
    // fields back and replace only the canonical XBHEAD vector.
    let scratch = SubstanceGroup::new(
        group.id(),
        cosmolkit_model::SubstanceGroupKind::StructuralRepeatUnit,
    );
    let owned_group = std::mem::replace(group, scratch);
    *group = owned_group.with_head_crossing_bonds(head.clone());

    for &bond in &tail {
        if bond.index() >= bond_count {
            return Err(PolymerSGroupError::BondOutOfRange { bond, bond_count });
        }
        group.push_bond(bond);
    }
    let mut correspondence = Vec::new();
    for index in 0..head.len().min(tail.len()) {
        correspondence.push(head[index]);
        let tail_index = if flipped {
            tail.len() - index - 1
        } else {
            index
        };
        correspondence.push(tail[tail_index]);
    }
    let scratch = SubstanceGroup::new(
        group.id(),
        cosmolkit_model::SubstanceGroupKind::StructuralRepeatUnit,
    );
    let owned_group = std::mem::replace(group, scratch);
    *group = owned_group.with_crossing_bond_correspondence(correspondence);
    Ok(())
}

#[cfg(test)]
fn fixed_property_text(value: &cosmolkit_model::PropertyText) -> &str {
    std::str::from_utf8(value.as_bytes()).expect("fixed fixture text is UTF8")
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
                Some(source.as_bytes()),
                &[BondId::new(0)],
                &[],
                1,
                1,
                |_| std::iter::empty::<(AtomId, BondId)>(),
            )
            .expect("source CONNECT fallback is non-failing");

            assert_eq!(group.connection(), Some(&expected_connection), "{source}");
            assert_eq!(
                group
                    .props()
                    .get("CONNECT".as_bytes())
                    .map(|v| v.as_string().expect("fixed StringTag"))
                    .map(super::fixed_property_text),
                Some(expected_prop)
            );
        }
    }

    #[test]
    fn query_sgroups_polymer_flip_removes_every_marker_and_reverses_tail_pairs() {
        let mut group = polymer_group(&[AtomId::new(0), AtomId::new(1)]);
        finalize_polymer_sgroup(
            &mut group,
            Some("hh,f,f".as_bytes()),
            &[BondId::new(1), BondId::new(1)],
            &[BondId::new(2), BondId::new(2), BondId::new(3)],
            2,
            4,
            |_| std::iter::empty::<(AtomId, BondId)>(),
        )
        .expect("repeated flip markers use source replace-all");

        assert_eq!(group.connection(), Some(&SGroupConnection::HeadToHead));
        assert_eq!(
            group
                .props()
                .get("CONNECT".as_bytes())
                .map(|v| v.as_string().expect("fixed StringTag"))
                .map(super::fixed_property_text),
            Some("HH")
        );
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
            Some("ht".as_bytes()),
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
        assert_eq!(
            group
                .props()
                .get("CONNECT".as_bytes())
                .map(|v| v.as_string().expect("fixed StringTag"))
                .map(super::fixed_property_text),
            Some("EU")
        );
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
            Some("eu".as_bytes()),
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
