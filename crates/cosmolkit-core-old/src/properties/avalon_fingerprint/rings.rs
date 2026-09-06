//! Avalon basis-ring implementation used by fingerprint preprocessing.

const BASE_BITS: usize = 16;
const UNLINKED: isize = -1;
const NO_COLOR: i32 = 0;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct BondSet {
    max_member: usize,
    bit_array: Vec<u32>,
}

impl BondSet {
    fn new(max_member: usize) -> Self {
        // Avalon❗✔️: bit_set_t *NewSet(unsigned int max_member)
        // Avalon❗✔️: /*
        // Avalon❗✔️:  * Allocate a new set. The set is initialized to cover values from
        // Avalon❗✔️:  * 0 to max_member.
        // Avalon❗✔️:  */
        // Avalon❗✔️: {
        // Avalon❗✔️:    bit_set_t *result;
        // Avalon❗✔️:
        // Avalon❗✔️:    result = TypeAlloc(1, bit_set_t);
        // Avalon❗✔️:    result->max_member = max_member;
        // Avalon❗✔️:    result->bit_array = TypeAlloc((max_member/BASE_BITS)+1, set_base_t);
        // Avalon❗✔️:    return(result);
        // Avalon❗✔️: }
        Self {
            max_member,
            bit_array: vec![0; max_member / BASE_BITS + 1],
        }
    }

    fn clear(&mut self) {
        // Avalon❗✔️: bit_set_t *ClearSet(bit_set_t *set)
        // Avalon❗✔️: /*
        // Avalon❗✔️:  * Clears all bits in *set.
        // Avalon❗✔️:  */
        // Avalon❗✔️: {
        // Avalon❗✔️:    register int i;
        // Avalon❗✔️:    register set_base_t *ip;
        // Avalon❗✔️:
        // Avalon❗✔️:    if (set == (bit_set_t *)NULL)
        // Avalon❗✔️:       ShowMessage("globbered set pointer","ClearSet");
        // Avalon❗✔️:    else
        // Avalon❗✔️:       for (i=set->max_member/BASE_BITS, ip=set->bit_array;
        // Avalon❗✔️:            i>=0;
        // Avalon❗✔️:            i--, ip++)
        // Avalon❗✔️:          *ip = 0;
        // Avalon❗✔️:    return(set);
        // Avalon❗✔️: }
        self.bit_array.fill(0);
    }

    fn put(&mut self, member: usize) {
        // Avalon❗✔️: bit_set_t *PutMember(bit_set_t *set, unsigned int member)
        // Avalon❗✔️: /*
        // Avalon❗✔️:  * Adds member to *set.
        // Avalon❗✔️:  */
        // Avalon❗✔️: {
        // Avalon❗✔️:    bit_set_t *sp;
        // Avalon❗✔️:
        // Avalon❗✔️:    if (set == (bit_set_t *)NULL || set->bit_array == (set_base_t *)NULL)
        // Avalon❗✔️:       ShowMessage("globbered set pointer","PutMember");
        // Avalon❗✔️:    else if (set->max_member < member)
        // Avalon❗✔️:    {
        // Avalon❗✔️:       sp = NewSet(member);
        // Avalon❗✔️:       sp = CopySet(sp, set);
        // Avalon❗✔️:       DisposeSet(set);
        // Avalon❗✔️:       set = sp;
        if member > self.max_member {
            self.bit_array.resize(member / BASE_BITS + 1, 0);
            self.max_member = member;
        }
        // Avalon❗✔️:       set->bit_array[member/BASE_BITS] |= ((set_base_t)1<<(member%BASE_BITS));
        // Avalon❗✔️:    }
        // Avalon❗✔️:    else
        // Avalon❗✔️:       set->bit_array[member/BASE_BITS] |= ((set_base_t)1<<(member%BASE_BITS));
        self.bit_array[member / BASE_BITS] |= 1_u32 << (member % BASE_BITS);
        // Avalon❗✔️:    return(set);
        // Avalon❗✔️: }
    }

    pub(super) fn contains(&self, member: usize) -> bool {
        // Avalon❗✔️: int IsMember(bit_set_t *set, unsigned member)
        // Avalon❗✔️: /*
        // Avalon❗✔️:  * TRUE if member is in *set
        // Avalon❗✔️:  */
        // Avalon❗✔️: {
        // Avalon❗✔️:    int result;
        // Avalon❗✔️:
        // Avalon❗✔️:    if (set == (bit_set_t *)NULL) {
        // Avalon❗✔️:       ShowMessage("globbered set pointer", "IsMember");
        // Avalon❗✔️:       return (FALSE);
        // Avalon❗✔️:    }
        // Avalon❗✔️:    else if (set->max_member < member  ||  set->bit_array[member/BASE_BITS] == 0)
        // Avalon❗✔️:       return (FALSE);
        // Avalon❗✔️:    else
        // Avalon❗✔️:       return (0 != (set->bit_array[member/BASE_BITS] & ((set_base_t)1<<(member%BASE_BITS))));
        // Avalon❗✔️: }
        member <= self.max_member
            && self.bit_array[member / BASE_BITS] & (1_u32 << (member % BASE_BITS)) != 0
    }

    pub(super) const fn max_member(&self) -> usize {
        self.max_member
    }

    fn copy_from(&mut self, source: &Self) {
        // Avalon❗✔️: bit_set_t *CopySet(bit_set_t *dest, bit_set_t *src)
        // Avalon❗✔️: /*
        // Avalon❗✔️:  * Copies *src onto *dest. However, it does _not_ allocate the
        // Avalon❗✔️:  * set, it only checks if destination is large enough.
        // Avalon❗✔️:  */
        // Avalon❗✔️: {
        // Avalon❗✔️:    register int i;
        // Avalon❗✔️:    register set_base_t *ip1, *ip2;
        // Avalon❗✔️:
        // Avalon❗✔️:    if (dest == (bit_set_t *)NULL || src == (bit_set_t *)NULL)
        // Avalon❗✔️:       ShowMessage("globbered set pointer","CopySet");
        // Avalon❗✔️:    else if (dest->max_member < src->max_member)
        // Avalon❗✔️:       ShowMessage("destination set size < source set size","CopySet");
        if self.max_member < source.max_member {
            return;
        }
        // Avalon❗✔️:    else
        // Avalon❗✔️:       for (i=src->max_member/BASE_BITS,
        // Avalon❗✔️:              ip1=dest->bit_array, ip2=src->bit_array;
        // Avalon❗✔️:            i>=0;
        // Avalon❗✔️:            i--, ip1++, ip2++)
        // Avalon❗✔️:          *ip1 = *ip2;
        let source_words = source.max_member / BASE_BITS + 1;
        self.bit_array[..source_words].copy_from_slice(&source.bit_array[..source_words]);
        // Avalon❗✔️:    return(dest);
        // Avalon❗✔️: }
    }

    fn exclusive_union_with(&mut self, other: &Self) {
        // Avalon❗✔️: bit_set_t *SetExclusiveUnion(bit_set_t *set1, bit_set_t *set2)
        // Avalon❗✔️: /*
        // Avalon❗✔️:  * Computes the exclusive set union of *set1 and *set2. The result is
        // Avalon❗✔️:  * a pointer to the changed *set1. The result is not allocated.
        // Avalon❗✔️:  */
        // Avalon❗✔️: {
        // Avalon❗✔️:    register int i;
        // Avalon❗✔️:
        // Avalon❗✔️:    if (set1 == (bit_set_t *)NULL || set2 == (bit_set_t *)NULL)
        // Avalon❗✔️:       ShowMessage("globbered set pointer","SetExclusiveUnion");
        // Avalon❗✔️:    else if (set1->max_member < set2->max_member)
        // Avalon❗✔️:       ShowMessage("destination set size < source set size",
        // Avalon❗✔️:                   "SetExclusiveUnion");
        if self.max_member < other.max_member {
            return;
        }
        // Avalon❗✔️:    else
        // Avalon❗✔️:       for (i=0; i<(set2->max_member/BASE_BITS)+1; i++)
        // Avalon❗✔️:          set1->bit_array[i] ^= set2->bit_array[i];
        let words = other.max_member / BASE_BITS + 1;
        for index in 0..words {
            self.bit_array[index] ^= other.bit_array[index];
        }
        // Avalon❗✔️:    return(set1);
        // Avalon❗✔️: }
    }

    fn intersection_with(&mut self, other: &Self) {
        // Avalon❗✔️: bit_set_t *SetIntersection(bit_set_t *set1, bit_set_t *set2)
        // Avalon❗✔️: /*
        // Avalon❗✔️:  * Computes the set intersection of *set1 and *set2. The result is a
        // Avalon❗✔️:  * pointer to the changed *set1. The result is not allocated.
        // Avalon❗✔️:  */
        // Avalon❗✔️: {
        // Avalon❗✔️:    register int i;
        // Avalon❗✔️:
        // Avalon❗✔️:    if (set1 == (bit_set_t *)NULL || set2 == (bit_set_t *)NULL)
        // Avalon❗✔️:       ShowMessage("globbered set pointer", "SetIntersection");
        // Avalon❗✔️:    else if (set1->max_member < set2->max_member)
        // Avalon❗✔️:       ShowMessage("destination set size < source set size","SetIntersection");
        if self.max_member < other.max_member {
            return;
        }
        // Avalon❗✔️:    else
        // Avalon❗✔️:    {
        // Avalon❗✔️:       for (i=0; i<(set2->max_member/BASE_BITS)+1; i++)
        // Avalon❗✔️:          set1->bit_array[i] &= set2->bit_array[i];
        let other_words = other.max_member / BASE_BITS + 1;
        for index in 0..other_words {
            self.bit_array[index] &= other.bit_array[index];
        }
        // Avalon❗✔️:       for (i=(set2->max_member/BASE_BITS)+1;
        // Avalon❗✔️:            i<(set1->max_member/BASE_BITS)+1; i++)
        // Avalon❗✔️:          set1->bit_array[i] = 0;
        self.bit_array[other_words..].fill(0);
        // Avalon❗✔️:    }
        // Avalon❗✔️:    return(set1);
        // Avalon❗✔️: }
    }

    fn intersection_is_empty(&self, other: &Self) -> bool {
        // Avalon❗✔️: int IntersectionIsEmpty(bit_set_t *set1, bit_set_t *set2)
        // Avalon❗✔️: /*
        // Avalon❗✔️:  * Tests if set1 and set2 don't share any bits.
        // Avalon❗✔️:  */
        // Avalon❗✔️: {
        // Avalon❗✔️:    register int i, n;
        // Avalon❗✔️:
        // Avalon❗✔️:    if (set1 == (bit_set_t *)NULL || set2 == (bit_set_t *)NULL)
        // Avalon❗✔️:       ShowMessage("globbered set pointer", "IntersectionIsEmpty");
        // Avalon❗✔️:    else
        // Avalon❗✔️:    {
        // Avalon❗✔️:       n = set1->max_member/BASE_BITS;
        // Avalon❗✔️:       if ((set2->max_member/BASE_BITS)>n) n = set2->max_member/BASE_BITS;
        // Avalon❗✔️:       n++;
        // Avalon❗✔️:
        // Avalon❗✔️:       for (i=0; i<n; i++)
        // Avalon❗✔️:           if (set1->bit_array[i] & set2->bit_array[i]) return (FALSE);
        let words = self.bit_array.len().max(other.bit_array.len());
        let result = (0..words).all(|index| {
            self.bit_array.get(index).copied().unwrap_or(0)
                & other.bit_array.get(index).copied().unwrap_or(0)
                == 0
        });
        // Avalon❗✔️:    }
        // Avalon❗✔️:    return(TRUE);
        // Avalon❗✔️: }
        result
    }

    fn cardinality(&self) -> usize {
        // Avalon❗✔️: int Cardinality(bit_set_t *set)
        // Avalon❗✔️: /*
        // Avalon❗✔️:  * Returns the number of elements in *set.
        // Avalon❗✔️:  */
        // Avalon❗✔️: {
        // Avalon❗✔️:    register unsigned int i, n;
        // Avalon❗✔️:
        // Avalon❗✔️:    for (i=0, n=0; i<=set->max_member; i++)
        // Avalon❗✔️:       if (ISMEMBER(set, i)) n++;
        // Avalon❗✔️:
        // Avalon❗✔️:    return(n);
        // Avalon❗✔️: }
        (0..=self.max_member)
            .filter(|&member| self.contains(member))
            .count()
    }

    #[cfg(test)]
    fn members(&self) -> Vec<usize> {
        (0..=self.max_member)
            .filter(|&member| self.contains(member))
            .collect()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct BondSetNode {
    pub(super) cardinality: usize,
    pub(super) bond_set: BondSet,
}

#[derive(Debug, Clone, Copy)]
struct TreeCell {
    color: i32,
    link: isize,
}

pub(super) fn ring_list(bonds: &[[usize; 2]]) -> Vec<BondSetNode> {
    // Avalon❗✔️: bond_set_node *RingList(unsigned bonds[][2], unsigned nbonds)
    // Avalon❗✔️: /*
    // Avalon❗✔️:  * Returns a list of basis rings of the graph defined by
    // Avalon❗✔️:  * bonds[0..nbonds-1][0..1]
    // Avalon❗✔️:  *
    // Avalon❗✔️:  * Generalized to ignore any bonds that contain nodes with number 0.
    // Avalon❗✔️:  * This can be used to perform ring analysis only on a selected subset
    // Avalon❗✔️:  * of the bonds of the source graph.
    // Avalon❗✔️:  */
    // Avalon❗✔️: {
    // Avalon❗✔️:    int i;
    // Avalon❗✔️:    unsigned int b;
    // Avalon❗✔️:    int at1, at2;
    // Avalon❗✔️:    int level1, level2;
    // Avalon❗✔️:    int trace, tmp;
    // Avalon❗✔️:    int color, old_color, new_color;
    // Avalon❗✔️:    int tmp_link, new_link;
    // Avalon❗✔️:    bond_set_node *result, *p;
    // Avalon❗✔️:    unsigned natoms;
    // Avalon❗✔️:    struct tree_cell *tree;
    // Avalon❗✔️:
    // Avalon❗✔️:    result = (bond_set_node *)NULL;
    let mut result = Vec::new();
    // Avalon❗✔️:
    // Avalon❗✔️:    for (i=0,natoms=0; i<nbonds; i++)         /* make space for spanning tree */
    // Avalon❗✔️:    {
    // Avalon❗✔️:       if (natoms < bonds[i][0]) natoms = bonds[i][0];
    // Avalon❗✔️:       if (natoms < bonds[i][1]) natoms = bonds[i][1];
    // Avalon❗✔️:    }
    let natoms = bonds
        .iter()
        .flat_map(|bond| bond.iter().copied())
        .max()
        .unwrap_or(0);
    // Avalon❗✔️:    tree = TypeAlloc(natoms+1,struct tree_cell);
    // Avalon❗✔️:
    // Avalon❗✔️:    for (i=0; i<=natoms; i++)
    // Avalon❗✔️:    {
    // Avalon❗✔️:       tree[i].color = NO_COLOR; tree[i].link = UNLINKED;
    // Avalon❗✔️:    }
    let mut tree = vec![
        TreeCell {
            color: NO_COLOR,
            link: UNLINKED
        };
        natoms + 1
    ];
    // Avalon❗✔️:
    // Avalon❗✔️:    /* Add bonds to spanning tree one at a time. If a bond doesn't link */
    // Avalon❗✔️:    /* to an old component, label the tree cells for both atoms with a  */
    // Avalon❗✔️:    /* new (unique) label. If one of the atoms corresponds to a cell    */
    // Avalon❗✔️:    /* that already has a label, label the cell corresponding to the    */
    // Avalon❗✔️:    /* other atom with the same label. If both atoms correspond to      */
    // Avalon❗✔️:    /* different labels, relabel the cells with the higher label to the */
    // Avalon❗✔️:    /* lower one to show, that they belong to the same component. If    */
    // Avalon❗✔️:    /* both cells are label the same, a ring is found. Then the links   */
    // Avalon❗✔️:    /* of the tree cells are followed to trace back to the common parent*/
    // Avalon❗✔️:    /* and the new ring, i.e. the set of bonds, is added to the result. */
    // Avalon❗✔️:
    // Avalon❗✔️:    for (b=0,color=NO_COLOR; b<nbonds; b++)
    // Avalon❗✔️:    {
    let mut color = NO_COLOR;
    for (bond_index, bond) in bonds.iter().enumerate() {
        // Avalon❗✔️:       at1 = bonds[b][0]; at2 = bonds[b][1];
        let (mut at1, mut at2) = (bond[0], bond[1]);
        // Avalon❗✔️:       if (at1 == 0  ||  at2 == 0) continue;  /* ignore bonds with atom 0 */
        if at1 == 0 || at2 == 0 {
            continue;
        }
        // Avalon❗✔️:       if (tree[at1].color == NO_COLOR &&     /* new component */
        // Avalon❗✔️:          tree[at2].color == NO_COLOR)
        // Avalon❗✔️:       {
        if tree[at1].color == NO_COLOR && tree[at2].color == NO_COLOR {
            // Avalon❗✔️:          color++; tree[at2].link = b;
            color += 1;
            tree[at2].link = bond_index as isize;
            // Avalon❗✔️:          tree[at1].color = tree[at2].color = color;
            tree[at1].color = color;
            tree[at2].color = color;
            // Avalon❗✔️:       }
            // Avalon❗✔️:       else if (tree[at1].color == NO_COLOR)           /* link first atom */
            // Avalon❗✔️:       {
        } else if tree[at1].color == NO_COLOR {
            // Avalon❗✔️:          tree[at1].color = tree[at2].color;
            tree[at1].color = tree[at2].color;
            // Avalon❗✔️:          tree[at1].link  = b;
            tree[at1].link = bond_index as isize;
            // Avalon❗✔️:       }
            // Avalon❗✔️:       else if (tree[at2].color == NO_COLOR)         /* link second atom */
            // Avalon❗✔️:       {
        } else if tree[at2].color == NO_COLOR {
            // Avalon❗✔️:          tree[at2].color = tree[at1].color;
            tree[at2].color = tree[at1].color;
            // Avalon❗✔️:          tree[at2].link  = b;
            tree[at2].link = bond_index as isize;
            // Avalon❗✔️:       }
            // Avalon❗✔️:       else if (tree[at1].color != tree[at2].color)  /* link two compnts. */
            // Avalon❗✔️:       {
        } else if tree[at1].color != tree[at2].color {
            // Avalon❗✔️:          new_color = tree[at1].color; old_color = tree[at2].color;
            let new_color = tree[at1].color;
            let old_color = tree[at2].color;
            // Avalon❗✔️:          for (i=0; i<=natoms; i++)
            // Avalon❗✔️:             if (tree[i].color == old_color) tree[i].color = new_color;
            for cell in &mut tree {
                if cell.color == old_color {
                    cell.color = new_color;
                }
            }
            // Avalon❗✔️:
            // Avalon❗✔️:          tmp_link = tree[at2].link;     /* trace the links of component 2 */
            let mut tmp_link = tree[at2].link;
            // Avalon❗✔️:          tree[at2].link = b;            /* and revert linkage             */
            tree[at2].link = bond_index as isize;
            // Avalon❗✔️:          while (tmp_link != UNLINKED)
            // Avalon❗✔️:          {
            while tmp_link != UNLINKED {
                // Avalon❗✔️:             at2 = TheOtherAtom(bonds[tmp_link],at2);
                at2 = other_atom(bonds[tmp_link as usize], at2);
                // Avalon❗✔️:             new_link = tmp_link;
                let new_link = tmp_link;
                // Avalon❗✔️:             tmp_link = tree[at2].link;
                tmp_link = tree[at2].link;
                // Avalon❗✔️:             tree[at2].link = new_link;
                tree[at2].link = new_link;
                // Avalon❗✔️:          }
            }
            // Avalon❗✔️:       }
            // Avalon❗✔️:       else                            /* ring found -> add it to list */
            // Avalon❗✔️:       {                               /* trace back to root from both atoms */
        } else {
            // Avalon❗✔️:          for (trace=at1,level1=0; tree[trace].link != UNLINKED; level1++)
            // Avalon❗✔️:             trace = TheOtherAtom(bonds[tree[trace].link],trace);
            let mut trace = at1;
            let mut level1 = 0;
            while tree[trace].link != UNLINKED {
                trace = other_atom(bonds[tree[trace].link as usize], trace);
                level1 += 1;
            }
            // Avalon❗✔️:          for (trace=at2,level2=0; tree[trace].link != UNLINKED; level2++)
            // Avalon❗✔️:             trace = TheOtherAtom(bonds[tree[trace].link],trace);
            trace = at2;
            let mut level2 = 0;
            while tree[trace].link != UNLINKED {
                trace = other_atom(bonds[tree[trace].link as usize], trace);
                level2 += 1;
            }
            // Avalon❗✔️:
            // Avalon❗✔️:          if (level1 > level2)   /* make path 1 the shorter one of the two */
            // Avalon❗✔️:          {
            if level1 > level2 {
                // Avalon❗✔️:             tmp = level1; level1 = level2; level2 = tmp;
                std::mem::swap(&mut level1, &mut level2);
                // Avalon❗✔️:             tmp = at1; at1 = at2; at2 = tmp;
                std::mem::swap(&mut at1, &mut at2);
                // Avalon❗✔️:          }
            }
            // Avalon❗✔️:
            // Avalon❗✔️:          p = NewBondSetNode(nbonds); p->next = result; result = p;
            // Avalon❗✔️: bond_set_node *NewBondSetNode(unsigned max_member)
            // Avalon❗✔️: /*
            // Avalon❗✔️:  * Allocates a new node.
            // Avalon❗✔️:  */
            // Avalon❗✔️: {
            // Avalon❗✔️:    bond_set_node *result;
            // Avalon❗✔️:
            // Avalon❗✔️:    result = TypeAlloc(1,bond_set_node);
            // Avalon❗✔️:    result->next = (bond_set_node *)NULL;
            // Avalon❗✔️:    result->cardinality = 0;
            // Avalon❗✔️:    result->bond_set = NewSet(max_member);
            // Avalon❗✔️:
            // Avalon❗✔️:    return(result);
            // Avalon❗✔️: }
            // Avalon❗✔️:          p->bond_set = ClearSet(p->bond_set); p->cardinality = 1;
            let mut bond_set = BondSet::new(bonds.len());
            bond_set.clear();
            // Avalon❗✔️:          p->bond_set = PutMember(p->bond_set,b);
            bond_set.put(bond_index);
            let mut cardinality = 1;
            // Avalon❗✔️:
            // Avalon❗✔️:          for (i=0; i<level2-level1; i++) /* trace back excess of long path */
            // Avalon❗✔️:          {
            for _ in 0..level2 - level1 {
                // Avalon❗✔️:             p->bond_set = PutMember(p->bond_set,(unsigned)tree[at2].link);
                bond_set.put(tree[at2].link as usize);
                // Avalon❗✔️:             p->cardinality++;
                cardinality += 1;
                // Avalon❗✔️:             at2 = TheOtherAtom(bonds[tree[at2].link],at2);
                at2 = other_atom(bonds[tree[at2].link as usize], at2);
                // Avalon❗✔️:          }
            }
            // Avalon❗✔️:
            // Avalon❗✔️:          while (at1 != at2)     /* simultaneously trace back both paths */
            // Avalon❗✔️:          {
            while at1 != at2 {
                // Avalon❗✔️:             p->bond_set = PutMember(p->bond_set,(unsigned)tree[at1].link);
                bond_set.put(tree[at1].link as usize);
                // Avalon❗✔️:             p->cardinality++;
                cardinality += 1;
                // Avalon❗✔️:             at1 = TheOtherAtom(bonds[tree[at1].link],at1);
                at1 = other_atom(bonds[tree[at1].link as usize], at1);
                // Avalon❗✔️:             p->bond_set = PutMember(p->bond_set,(unsigned)tree[at2].link);
                bond_set.put(tree[at2].link as usize);
                // Avalon❗✔️:             p->cardinality++;
                cardinality += 1;
                // Avalon❗✔️:             at2 = TheOtherAtom(bonds[tree[at2].link],at2);
                at2 = other_atom(bonds[tree[at2].link as usize], at2);
                // Avalon❗✔️:          }
            }
            // The C linked list prepends in O(1). Appending then reversing once
            // preserves its exact order while avoiding repeated Vec front insertion.
            result.push(BondSetNode {
                cardinality,
                bond_set,
            });
            // Avalon❗✔️:       }       /* else ring found */
        }
        // Avalon❗✔️:    }       /* for all bonds */
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    free((char *)tree);
    // Rust releases the tree vector on scope exit.
    // Avalon❗✔️:    return(result);
    // Avalon❗✔️: }
    result.reverse();
    result
}

fn other_atom(bond: [usize; 2], atom: usize) -> usize {
    bond[0] + bond[1] - atom
}

fn sort_rings(rings: &mut [BondSetNode]) {
    // Avalon❗✔️: bond_set_node *SortRings(bond_set_node *list)
    // Avalon❗✔️: /*
    // Avalon❗✔️:  * Sorts *list into descending order with respect to cardinality.
    // Avalon❗✔️:  */
    // Avalon❗✔️: {
    // Avalon❗✔️:    bond_set_node *p1, *p2;
    // Avalon❗✔️:    bit_set_t *set;
    // Avalon❗✔️:    int size;
    // Avalon❗✔️:
    // Avalon❗✔️:    for (p1 = list; p1 != (bond_set_node *)NULL; p1 = p1->next)
    // Avalon❗✔️:    {
    for first in 0..rings.len() {
        // Avalon❗✔️:       for (p2 = p1->next; p2 != (bond_set_node *)NULL; p2 = p2->next)
        // Avalon❗✔️:          if (p1->cardinality > p2->cardinality)  /* out of order -> exchange */
        // Avalon❗✔️:          {
        for second in first + 1..rings.len() {
            if rings[first].cardinality > rings[second].cardinality {
                // Avalon❗✔️:             size            = p1->cardinality;
                // Avalon❗✔️:             p1->cardinality = p2->cardinality;
                // Avalon❗✔️:             p2->cardinality = size;
                // Avalon❗✔️:             set          = p1->bond_set;
                // Avalon❗✔️:             p1->bond_set = p2->bond_set;
                // Avalon❗✔️:             p2->bond_set = set;
                rings.swap(first, second);
                // Avalon❗✔️:          }
            }
        }
        // Avalon❗✔️:    }
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    return(list);
    // Avalon❗✔️: }
}

pub(super) fn combine_rings(rings: &mut Vec<BondSetNode>) {
    // Avalon❗✔️: bond_set_node *CombineRings(bond_set_node *list)
    // Avalon❗✔️: /*
    // Avalon❗✔️:  * Combines pairs of rings until selfconsistency to get a list of
    // Avalon❗✔️:  * smaller basis rings.
    // Avalon❗✔️:  */
    // Avalon❗✔️: {
    // Avalon❗✔️:    bond_set_node *p1, *p2;
    // Avalon❗✔️:    int size;
    // Avalon❗✔️:    bit_set_t *set,*tmp;
    // Avalon❗✔️:    int changed;
    // Avalon❗✔️:    int ntoggle;
    // Avalon❗✔️:
    // Avalon❗✔️:    if (list == (bond_set_node *)NULL) return(list);
    if rings.is_empty() {
        return;
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    srand(1);         /* make sure that algorithm works reproducibly */
    let mut random = GlibcRandom::seeded(1);
    // Avalon❗✔️:    ntoggle = 0;      /* safeguard against infinite looping */
    let mut ntoggle = 0;
    // Avalon❗✔️:    set = NewSet(MaxMember(list->bond_set));  /* assumes allocated size is */
    let mut scratch = BondSet::new(rings[0].bond_set.max_member);
    // Avalon❗✔️:    do                                        /* the same for all rings!!! */
    // Avalon❗✔️:    {
    loop {
        // Avalon❗✔️:       changed = FALSE;
        let mut changed = false;
        // Avalon❗✔️:       list = SortRings(list);  /* loop over all pairs of different rings */
        sort_rings(rings);
        // Avalon❗✔️:       for (p1=list; p1!=(bond_set_node *)NULL; p1=p1->next)
        // Avalon❗✔️:          for (p2=p1->next; p2!=(bond_set_node *)NULL; p2=p2->next)
        // Avalon❗✔️:          {
        for first in 0..rings.len() {
            for second in first + 1..rings.len() {
                // Avalon❗✔️:             // check first if there is any overlap
                // Avalon❗✔️:             if (IntersectionIsEmpty(p1->bond_set,p2->bond_set)) continue;;
                if rings[first]
                    .bond_set
                    .intersection_is_empty(&rings[second].bond_set)
                {
                    continue;
                }
                // Avalon❗✔️:             set = SetExclusiveUnion(CopySet(set,p1->bond_set),p2->bond_set);
                scratch.copy_from(&rings[first].bond_set);
                scratch.exclusive_union_with(&rings[second].bond_set);
                // Avalon❗✔️:             size = Cardinality(set);
                let size = scratch.cardinality();
                // Avalon❗✔️:             if (size > 0 &&
                // Avalon❗✔️:                 (size <= p1->cardinality || size <= p2->cardinality))
                // Avalon❗✔️:             {
                if size > 0
                    && (size <= rings[first].cardinality || size <= rings[second].cardinality)
                {
                    // Avalon❗✔️:                if (p1->cardinality > p2->cardinality)
                    // Avalon❗✔️:                {
                    let target = if rings[first].cardinality > rings[second].cardinality {
                        // Avalon❗✔️:                   if (p1->cardinality > size  || (rand() / 10) % 2)
                        // Avalon❗✔️:                   {
                        if rings[first].cardinality > size || (random.next() / 10) % 2 != 0 {
                            Some(first)
                        } else {
                            None
                        }
                        // Avalon❗✔️:                }
                        // Avalon❗✔️:                else
                        // Avalon❗✔️:                {
                    } else {
                        // Avalon❗✔️:                   if (p2->cardinality > size  || (rand() / 10) % 2)
                        // Avalon❗✔️:                   {
                        if rings[second].cardinality > size || (random.next() / 10) % 2 != 0 {
                            Some(second)
                        } else {
                            None
                        }
                        // Avalon❗✔️:                   }
                        // Avalon❗✔️:                }
                    };
                    if let Some(target) = target {
                        // Avalon❗✔️:                      if (p1->cardinality == size)
                        // Avalon❗✔️:                         ntoggle++;
                        // Avalon❗✔️:                      else
                        // Avalon❗✔️:                         ntoggle = 0;
                        // Avalon❗✔️:                      changed = TRUE; p1->cardinality = size;
                        // Avalon❗✔️:                      tmp = p1->bond_set; p1->bond_set = set; set = tmp;
                        // Avalon❗✔️:                      if (p2->cardinality == size)
                        // Avalon❗✔️:                         ntoggle++;
                        // Avalon❗✔️:                      else
                        // Avalon❗✔️:                         ntoggle = 0;
                        // Avalon❗✔️:                      changed = TRUE; p2->cardinality = size;
                        // Avalon❗✔️:                      tmp = p2->bond_set; p2->bond_set = set; set = tmp;
                        if rings[target].cardinality == size {
                            ntoggle += 1;
                        } else {
                            ntoggle = 0;
                        }
                        changed = true;
                        rings[target].cardinality = size;
                        std::mem::swap(&mut rings[target].bond_set, &mut scratch);
                    }
                    // The source duplicates the preceding update block for p2;
                    // `target` selects the same branch without changing call order.
                    // Avalon❗✔️:             }
                }
                // Avalon❗✔️:          }
            }
        }
        // Avalon❗✔️:       if (ntoggle > 4) changed = FALSE;          /* limit to 4 toggles */
        if ntoggle > 4 {
            changed = false;
        }
        // Avalon❗✔️:    } while(changed);
        if !changed {
            break;
        }
    }
    // Avalon❗✔️:
    // Avalon❗✔️:    DisposeSet(set);
    // Rust releases the scratch set on scope exit.
    // Avalon❗✔️:    return(list);
    // Avalon❗✔️: }
}

pub(super) fn proper_ring_pairs(
    base_rings: &[BondSetNode],
    maxnode: usize,
    bonds: &[[usize; 2]],
) -> Vec<BondSetNode> {
    // Avalon❗✔️: bond_set_node *ProperRingPairs(bond_set_node *base_rings, int maxnode,
    // Avalon❗✔️:                                unsigned bonds[][2])
    // Avalon❗✔️: /*
    // Avalon❗✔️:  * Returns the list of all proper rings that can be generated by XORing
    // Avalon❗✔️:  * the base_rings. Note: It may happen that the XOR of two proper rings
    // Avalon❗✔️:  * yields a disconnected set of more than one ring!
    // Avalon❗✔️:  * maxnode is needed to speed up internal operations.
    // Avalon❗✔️:  */
    // Avalon❗✔️: {
    // Avalon❗✔️:    bond_set_node *p1, *p2, *result, *ph;
    // Avalon❗✔️:    int i, b;
    // Avalon❗✔️:    int nterminal;
    // Avalon❗✔️:    bit_set_t *set;
    // Avalon❗✔️:    int *atom_touched;
    // Avalon❗✔️:
    // Avalon❗✔️:    result = (bond_set_node *)NULL;
    let mut result = Vec::new();
    // Avalon❗✔️:    if (base_rings == (bond_set_node *)NULL) return result;
    let Some(first_ring) = base_rings.first() else {
        return result;
    };
    // Avalon❗✔️:    atom_touched = TypeAlloc(maxnode+1, int);
    let mut atom_touched = vec![0_i32; maxnode + 1];
    // Avalon❗✔️:
    // Avalon❗✔️:    /* assumes allocated size is */
    // Avalon❗✔️:    set = NewSet(MaxMember(base_rings->bond_set));
    let mut intersection = BondSet::new(first_ring.bond_set.max_member);
    // Avalon❗✔️:    for (p1=base_rings; p1 != (bond_set_node *)NULL; p1=p1->next)
    // Avalon❗✔️:       for (p2=p1->next; p2 != (bond_set_node *)NULL; p2=p2->next)
    // Avalon❗✔️:       {
    for first in 0..base_rings.len() {
        for second in first + 1..base_rings.len() {
            // Avalon❗✔️:          // check first if there is any overlap
            // Avalon❗✔️:          if (IntersectionIsEmpty(p1->bond_set,p2->bond_set)) continue;;
            if base_rings[first]
                .bond_set
                .intersection_is_empty(&base_rings[second].bond_set)
            {
                continue;
            }
            // Avalon❗✔️:          SetIntersection(CopySet(set,p1->bond_set),p2->bond_set);
            intersection.copy_from(&base_rings[first].bond_set);
            intersection.intersection_with(&base_rings[second].bond_set);
            // Avalon❗✔️:          /* Now, we test if there is only one connected path in common */
            // Avalon❗✔️:          /* between p1 and p2. This is the case if exactly 2 atoms referred */
            // Avalon❗✔️:          /* to by the AND of the bond sets bonds are listed once. */
            // Avalon❗✔️:          for (i=0; i<maxnode+1; i++) atom_touched[i] = 0;
            atom_touched.fill(0);
            // Avalon❗✔️:          for (b=0; b<MaxMember(set); b++)
            // Avalon deliberately excludes the allocated maximum member here.
            for bond_index in 0..intersection.max_member {
                // Avalon❗✔️:             if (IsMember(set,b))
                // Avalon❗✔️:             {
                if intersection.contains(bond_index) {
                    // Avalon❗✔️:                atom_touched[bonds[b][0]]++;
                    atom_touched[bonds[bond_index][0]] += 1;
                    // Avalon❗✔️:                atom_touched[bonds[b][1]]++;
                    atom_touched[bonds[bond_index][1]] += 1;
                    // Avalon❗✔️:             }
                }
            }
            // Avalon❗✔️:          nterminal = 0;
            // Avalon❗✔️:          for (i=0; i<maxnode+1; i++)
            // Avalon❗✔️:             if (atom_touched[i]==1) nterminal++;
            let nterminal = atom_touched.iter().filter(|&&count| count == 1).count();
            // Avalon❗✔️:          if (nterminal == 2)
            // Avalon❗✔️:          {
            if nterminal == 2 {
                // Avalon❗✔️:             ph = NewBondSetNode(MaxMember(p1->bond_set));
                let mut bond_set = BondSet::new(base_rings[first].bond_set.max_member);
                // Avalon❗✔️:             SetExclusiveUnion(CopySet(ph->bond_set,p1->bond_set),p2->bond_set);
                bond_set.copy_from(&base_rings[first].bond_set);
                bond_set.exclusive_union_with(&base_rings[second].bond_set);
                // Avalon❗✔️:             ph->cardinality = Cardinality(ph->bond_set);
                let cardinality = bond_set.cardinality();
                // Avalon❗✔️:             ph->next = result;
                // Avalon❗✔️:             result = ph;
                result.push(BondSetNode {
                    cardinality,
                    bond_set,
                });
                // Avalon❗✔️:          }
            }
            // Avalon❗✔️:       }
        }
    }
    // `push` records enumeration order; one linear reversal reproduces the
    // source's constant-time linked-list prepends without quadratic inserts.
    result.reverse();
    // Avalon❗✔️:
    // Avalon❗✔️:    DisposeSet(set);
    // Avalon❗✔️:    MyFree((char *)atom_touched);
    // Rust releases both scratch allocations on scope exit.
    // Avalon❗✔️:    return (result);
    // Avalon❗✔️: }
    result
}

pub(super) fn prepend_fused_ring_pairs(rings: &mut Vec<BondSetNode>) {
    if rings.is_empty() {
        return;
    }
    // Avalon❗✔️:    /* add fused pairs that could be also aromatic */
    // Avalon❗✔️:    set = NewSet(MaxMember(ring_list->bond_set));
    let mut scratch = BondSet::new(rings[0].bond_set.max_member);
    // Avalon❗✔️:    aromatic_candidates = ring_list;
    let original_len = rings.len();
    let mut fused = Vec::new();
    // Avalon❗✔️:    for (plist=ring_list; plist!=(bond_set_node *)NULL; plist=plist->next)
    // Avalon❗✔️:       for (plisth=plist->next;
    // Avalon❗✔️:            plisth!=(bond_set_node *)NULL;
    // Avalon❗✔️: 	   plisth=plisth->next)
    // Avalon❗✔️:       {
    for first in 0..original_len {
        for second in first + 1..original_len {
            // Avalon❗✔️:          // check first if there is any overlap
            // Avalon❗✔️:          if (IntersectionIsEmpty(plist->bond_set,plisth->bond_set)) continue;;
            if rings[first]
                .bond_set
                .intersection_is_empty(&rings[second].bond_set)
            {
                continue;
            }
            // Avalon❗✔️:          /* assumes allocated size is large enough */
            // Avalon❗✔️:          set = SetExclusiveUnion(CopySet(set,plist->bond_set),plisth->bond_set);
            scratch.copy_from(&rings[first].bond_set);
            scratch.exclusive_union_with(&rings[second].bond_set);
            // Avalon❗✔️: 	 /* only fused rings are considered */
            // Avalon❗✔️: 	 if (Cardinality(set) == plist->cardinality +
            // Avalon❗✔️: 	                         plisth->cardinality - 2)
            // Avalon❗✔️: 	 {
            let cardinality = scratch.cardinality();
            if cardinality
                == rings[first]
                    .cardinality
                    .wrapping_add(rings[second].cardinality)
                    .wrapping_sub(2)
            {
                // Avalon❗✔️: 	    p = NewBondSetNode(MaxMember(plist->bond_set));
                // Avalon❗✔️: 	    p->next = aromatic_candidates;
                // Avalon❗✔️: 	    aromatic_candidates = p;
                // Avalon❗✔️: 	    CopySet(p->bond_set, set);
                // Avalon❗✔️:             p->cardinality = Cardinality(set);
                fused.push(BondSetNode {
                    cardinality,
                    bond_set: scratch.clone(),
                });
                // Avalon❗✔️: 	 }
            }
            // Avalon❗✔️:       }
        }
    }
    // Source prepends each generated node. Reversing once preserves that
    // linked-list order without quadratic front insertion into the vector.
    fused.reverse();
    fused.append(rings);
    *rings = fused;
    // Avalon❗✔️:    DisposeSet(set);
    // Rust releases the reusable XOR scratch set on scope exit.
}

/// Local form of glibc's TYPE_3 additive generator used by `rand()`.
///
/// Avalon reseeds the process-global generator for every ring combination.
/// Keeping the identical state per invocation preserves that sequence while
/// making concurrent fingerprint calls independent.
struct GlibcRandom {
    state: [i32; 31],
    front: usize,
    rear: usize,
}

impl GlibcRandom {
    fn seeded(seed: u32) -> Self {
        // glibc❗✔️:   state = buf->state;
        // glibc❗✔️:   /* We must make sure the seed is not 0.  Take arbitrarily 1 in this case.  */
        // glibc❗✔️:   if (seed == 0)
        // glibc❗✔️:     seed = 1;
        let seed = if seed == 0 { 1 } else { seed };
        // glibc❗✔️:   write_state (state, 0, seed);
        let mut state = [0_i32; 31];
        state[0] = seed as i32;
        // glibc❗✔️:   if (type == TYPE_0)
        // glibc❗✔️:     goto done;
        // Avalon always reaches this helper through glibc's default TYPE_3 state.
        // glibc❗✔️:
        // glibc❗✔️:   dst = state;
        // glibc❗✔️:   word = seed;
        let mut word = i64::from(state[0]);
        // glibc❗✔️:   kc = buf->rand_deg;
        // glibc❗✔️:   for (i = 1; i < kc; ++i)
        // glibc❗✔️:     {
        // glibc❗✔️:       /* This does:
        // glibc❗✔️: 	   state[i] = (16807 * state[i - 1]) % 2147483647;
        // glibc❗✔️: 	 but avoids overflowing 31 bits.  */
        for slot in state.iter_mut().skip(1) {
            // glibc❗✔️:       long int hi = word / 127773;
            let high = word / 127_773;
            // glibc❗✔️:       long int lo = word % 127773;
            let low = word % 127_773;
            // glibc❗✔️:       word = 16807 * lo - 2836 * hi;
            word = 16_807 * low - 2_836 * high;
            // glibc❗✔️:       if (word < 0)
            // glibc❗✔️: 	word += 2147483647;
            if word < 0 {
                word += 2_147_483_647;
            }
            // glibc❗✔️:       write_state (++dst, 0, word);
            *slot = word as i32;
            // glibc❗✔️:     }
        }
        // glibc❗✔️:
        // glibc❗✔️:   buf->fptr = &state[buf->rand_sep];
        // glibc❗✔️:   buf->rptr = &state[0];
        let mut result = Self {
            state,
            front: 3,
            rear: 0,
        };
        // glibc❗✔️:   kc *= 10;
        // glibc❗✔️:   while (--kc >= 0)
        // glibc❗✔️:     {
        // glibc❗✔️:       int32_t discard;
        // glibc❗✔️:       (void) __random_r (buf, &discard);
        // glibc❗✔️:     }
        for _ in 0..310 {
            result.next();
        }
        result
    }

    fn next(&mut self) -> u32 {
        // glibc❗✔️:   else
        // glibc❗✔️:     {
        // glibc❗✔️:       int32_t *fptr = buf->fptr;
        // glibc❗✔️:       int32_t *rptr = buf->rptr;
        // glibc❗✔️:       int32_t *end_ptr = buf->end_ptr;
        // glibc❗✔️:       uint32_t val;
        // glibc❗✔️:
        // glibc❗✔️:       /* Avoid integer overflow with uint32_t arihmetic.  */
        // glibc❗✔️:       val = read_state (fptr, 0);
        // glibc❗✔️:       val += read_state (rptr, 0);
        let value = (self.state[self.front] as u32).wrapping_add(self.state[self.rear] as u32);
        // glibc❗✔️:       write_state (fptr, 0, val);
        self.state[self.front] = value as i32;
        // glibc❗✔️:       /* Chucking least random bit.  */
        // glibc❗✔️:       *result = val >> 1;
        let result = (value >> 1) & 0x7fff_ffff;
        // glibc❗✔️:       ++fptr;
        self.front += 1;
        // glibc❗✔️:       if (fptr >= end_ptr)
        // glibc❗✔️: 	{
        // glibc❗✔️: 	  fptr = state;
        // glibc❗✔️: 	  ++rptr;
        // glibc❗✔️: 	}
        if self.front == self.state.len() {
            self.front = 0;
            self.rear += 1;
            // glibc❗✔️:       else
            // glibc❗✔️: 	{
            // glibc❗✔️: 	  ++rptr;
            // glibc❗✔️: 	  if (rptr >= end_ptr)
            // glibc❗✔️: 	    rptr = state;
            // glibc❗✔️: 	}
        } else {
            self.rear += 1;
            if self.rear == self.state.len() {
                self.rear = 0;
            }
        }
        // glibc❗✔️:       buf->fptr = fptr;
        // glibc❗✔️:       buf->rptr = rptr;
        // glibc❗✔️:     }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn member_lists(rings: &[BondSetNode]) -> Vec<Vec<usize>> {
        rings.iter().map(|ring| ring.bond_set.members()).collect()
    }

    #[test]
    fn local_rng_matches_glibc_rand_after_srand_one() {
        let mut random = GlibcRandom::seeded(1);
        let values: Vec<u32> = (0..5).map(|_| random.next()).collect();
        assert_eq!(
            values,
            [
                1_804_289_383,
                846_930_886,
                1_681_692_777,
                1_714_636_915,
                1_957_747_793
            ]
        );
    }

    #[test]
    fn ring_list_preserves_source_prepend_order_for_disconnected_cycles() {
        let bonds = [[1, 2], [2, 3], [3, 1], [4, 5], [5, 6], [6, 4]];
        let rings = ring_list(&bonds);
        assert_eq!(member_lists(&rings), vec![vec![3, 4, 5], vec![0, 1, 2]]);
        assert_eq!(
            rings
                .iter()
                .map(|ring| ring.cardinality)
                .collect::<Vec<_>>(),
            vec![3, 3]
        );
    }

    #[test]
    fn ring_list_ignores_bonds_with_zero_endpoints() {
        let rings = ring_list(&[[1, 2], [2, 3], [3, 1], [0, 4], [4, 1]]);
        assert_eq!(member_lists(&rings), vec![vec![0, 1, 2]]);
    }

    #[test]
    fn combine_rings_reduces_overlapping_basis_cycles() {
        let bonds = [[1, 2], [2, 3], [3, 4], [4, 1], [1, 3]];
        let mut rings = ring_list(&bonds);
        combine_rings(&mut rings);

        assert_eq!(rings.len(), 2);
        assert_eq!(
            rings
                .iter()
                .map(|ring| ring.cardinality)
                .collect::<Vec<_>>(),
            vec![3, 3]
        );
        assert_eq!(member_lists(&rings), vec![vec![0, 1, 4], vec![2, 3, 4]]);
    }

    #[test]
    fn native_fused_and_bridged_basis_member_order_is_locked() {
        let fused = [[1, 2], [2, 3], [3, 4], [4, 1], [3, 5], [5, 6], [6, 4]];
        let mut rings = ring_list(&fused);
        assert_eq!(
            member_lists(&rings),
            vec![vec![2, 4, 5, 6], vec![0, 1, 2, 3]]
        );
        combine_rings(&mut rings);
        assert_eq!(
            member_lists(&rings),
            vec![vec![2, 4, 5, 6], vec![0, 1, 2, 3]]
        );

        let bridged = [[1, 2], [2, 3], [3, 1], [3, 4], [4, 5], [5, 6], [6, 4]];
        let mut rings = ring_list(&bridged);
        assert_eq!(member_lists(&rings), vec![vec![4, 5, 6], vec![0, 1, 2]]);
        combine_rings(&mut rings);
        assert_eq!(member_lists(&rings), vec![vec![4, 5, 6], vec![0, 1, 2]]);
    }

    #[test]
    fn proper_ring_pairs_preserves_source_prepend_order() {
        let bonds = [[1, 2], [2, 3], [3, 4], [4, 1], [3, 5], [5, 6], [6, 4]];
        let mut rings = ring_list(&bonds);
        combine_rings(&mut rings);

        let pairs = proper_ring_pairs(&rings, 6, &bonds);

        assert_eq!(member_lists(&pairs), vec![vec![0, 1, 3, 4, 5, 6]]);
        assert_eq!(pairs[0].cardinality, 6);
    }

    #[test]
    fn proper_ring_pairs_keeps_source_max_member_exclusion() {
        let bonds = [[1, 2], [2, 3], [3, 1]];
        let mut first = BondSet::new(2);
        first.put(0);
        first.put(2);
        let mut second = BondSet::new(2);
        second.put(1);
        second.put(2);
        let rings = vec![
            BondSetNode {
                cardinality: 2,
                bond_set: first,
            },
            BondSetNode {
                cardinality: 2,
                bond_set: second,
            },
        ];

        // The only shared bond is member 2, equal to MaxMember(set), so the
        // source endpoint-count loop does not see it and rejects the pair.
        assert!(proper_ring_pairs(&rings, 3, &bonds).is_empty());
    }

    #[test]
    fn multiple_proper_ring_pairs_match_native_prepend_order() {
        let bonds = [
            [1, 2],
            [2, 3],
            [3, 4],
            [4, 1],
            [3, 5],
            [5, 6],
            [6, 4],
            [5, 7],
            [7, 8],
            [8, 6],
        ];
        let mut rings = ring_list(&bonds);
        combine_rings(&mut rings);

        let pairs = proper_ring_pairs(&rings, 8, &bonds);

        assert_eq!(
            member_lists(&pairs),
            vec![vec![0, 1, 3, 4, 5, 6], vec![2, 4, 6, 7, 8, 9]]
        );
    }

    #[test]
    fn fused_ring_pair_duplicates_and_prepend_order_match_native() {
        let bonds = [
            [1, 2],
            [2, 3],
            [3, 4],
            [4, 1],
            [3, 5],
            [5, 6],
            [6, 4],
            [5, 7],
            [7, 8],
            [8, 6],
        ];
        let mut rings = ring_list(&bonds);
        combine_rings(&mut rings);
        let mut proper = proper_ring_pairs(&rings, 8, &bonds);
        proper.reverse();
        proper.append(&mut rings);
        rings = proper;

        prepend_fused_ring_pairs(&mut rings);

        assert_eq!(
            member_lists(&rings),
            vec![
                vec![0, 1, 3, 4, 5, 6],
                vec![2, 4, 6, 7, 8, 9],
                vec![0, 1, 3, 4, 6, 7, 8, 9],
                vec![0, 1, 3, 4, 6, 7, 8, 9],
                vec![2, 4, 6, 7, 8, 9],
                vec![0, 1, 3, 4, 5, 6],
                vec![5, 7, 8, 9],
                vec![2, 4, 5, 6],
                vec![0, 1, 2, 3],
            ]
        );
    }
}
