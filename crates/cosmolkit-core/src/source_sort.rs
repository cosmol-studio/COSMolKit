//! Private source-defined libstdc++ ordering shared by ring and wedge owners.
//! This generalizes the existing ring-size introsort; no chemistry is implemented here.

pub(crate) fn sort_by<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], mut less: F) {
    // libstdc++❗✔️: __sort(_RandomAccessIterator __first, _RandomAccessIterator __last,
    // libstdc++❗✔️: 	   _Compare __comp)
    // libstdc++❗✔️:     {
    // libstdc++❗✔️:       if (__first != __last)
    // libstdc++❗✔️: 	{
    // libstdc++❗✔️: 	  std::__introsort_loop(__first, __last,
    // libstdc++❗✔️: 				std::__lg(__last - __first) * 2,
    // libstdc++❗✔️: 				__comp);
    // libstdc++❗✔️: 	  std::__final_insertion_sort(__first, __last, __comp);
    // libstdc++❗✔️: 	}
    // libstdc++❗✔️:     }
    // Behavior: preserve the pinned source's equivalent-key permutations.
    // Complexity: O(n log n), logarithmic stack, no new storage or element cloning.
    if v.len() < 2 {
        return;
    }
    let depth = (usize::BITS - 1 - v.len().leading_zeros()) as usize * 2;
    introsort(v, 0, v.len(), depth, &mut less);
    final_insertion(v, 0, v.len(), &mut less);
}

fn introsort<T, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    mut last: usize,
    mut depth: usize,
    less: &mut F,
) {
    // libstdc++❗✔️: __introsort_loop(_RandomAccessIterator __first,
    // libstdc++❗✔️: 		     _RandomAccessIterator __last,
    // libstdc++❗✔️: 		     _Size __depth_limit, _Compare __comp)
    // libstdc++❗✔️:     {
    // libstdc++❗✔️:       while (__last - __first > int(_S_threshold))
    // libstdc++❗✔️: 	{
    // libstdc++❗✔️: 	  if (__depth_limit == 0)
    // libstdc++❗✔️: 	    {
    // libstdc++❗✔️: 	      std::__partial_sort(__first, __last, __last, __comp);
    // libstdc++❗✔️: 	      return;
    // libstdc++❗✔️: 	    }
    // libstdc++❗✔️: 	  --__depth_limit;
    // libstdc++❗✔️: 	  _RandomAccessIterator __cut =
    // libstdc++❗✔️: 	    std::__unguarded_partition_pivot(__first, __last, __comp);
    // libstdc++❗✔️: 	  std::__introsort_loop(__cut, __last, __depth_limit, __comp);
    // libstdc++❗✔️: 	  __last = __cut;
    // libstdc++❗✔️: 	}
    // libstdc++❗✔️:     }
    while last - first > 16 {
        if depth == 0 {
            heap_sort(&mut v[first..last], less);
            return;
        }
        depth -= 1;
        let cut = partition_pivot(v, first, last, less);
        introsort(v, cut, last, depth, less);
        last = cut;
    }
}

fn partition_pivot<T, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    less: &mut F,
) -> usize {
    // libstdc++❗✔️: __unguarded_partition_pivot(_RandomAccessIterator __first,
    // libstdc++❗✔️: 				_RandomAccessIterator __last, _Compare __comp)
    // libstdc++❗✔️:     {
    // libstdc++❗✔️:       _RandomAccessIterator __mid = __first + (__last - __first) / 2;
    // libstdc++❗✔️:       std::__move_median_to_first(__first, __first + 1, __mid, __last - 1,
    // libstdc++❗✔️: 				  __comp);
    // libstdc++❗✔️:       return std::__unguarded_partition(__first + 1, __last, __first, __comp);
    // libstdc++❗✔️:     }
    let mid = first + (last - first) / 2;
    median(v, first, first + 1, mid, last - 1, less);
    partition(v, first + 1, last, first, less)
}

fn median<T, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    result: usize,
    a: usize,
    b: usize,
    c: usize,
    less: &mut F,
) {
    // libstdc++❗✔️: __move_median_to_first(_Iterator __result,_Iterator __a, _Iterator __b,
    // libstdc++❗✔️: 			   _Iterator __c, _Compare __comp)
    // libstdc++❗✔️:     {
    // libstdc++❗✔️:       if (__comp(__a, __b))
    // libstdc++❗✔️: 	{
    // libstdc++❗✔️: 	  if (__comp(__b, __c))
    // libstdc++❗✔️: 	    std::iter_swap(__result, __b);
    // libstdc++❗✔️: 	  else if (__comp(__a, __c))
    // libstdc++❗✔️: 	    std::iter_swap(__result, __c);
    // libstdc++❗✔️: 	  else
    // libstdc++❗✔️: 	    std::iter_swap(__result, __a);
    // libstdc++❗✔️: 	}
    // libstdc++❗✔️:       else if (__comp(__a, __c))
    // libstdc++❗✔️: 	std::iter_swap(__result, __a);
    // libstdc++❗✔️:       else if (__comp(__b, __c))
    // libstdc++❗✔️: 	std::iter_swap(__result, __c);
    // libstdc++❗✔️:       else
    // libstdc++❗✔️: 	std::iter_swap(__result, __b);
    // libstdc++❗✔️:     }
    if less(&v[a], &v[b]) {
        if less(&v[b], &v[c]) {
            v.swap(result, b);
        } else if less(&v[a], &v[c]) {
            v.swap(result, c);
        } else {
            v.swap(result, a);
        }
    } else if less(&v[a], &v[c]) {
        v.swap(result, a);
    } else if less(&v[b], &v[c]) {
        v.swap(result, c);
    } else {
        v.swap(result, b);
    }
}

fn partition<T, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    mut first: usize,
    mut last: usize,
    pivot: usize,
    less: &mut F,
) -> usize {
    // libstdc++❗✔️: __unguarded_partition(_RandomAccessIterator __first,
    // libstdc++❗✔️: 			  _RandomAccessIterator __last,
    // libstdc++❗✔️: 			  _RandomAccessIterator __pivot, _Compare __comp)
    // libstdc++❗✔️:     {
    // libstdc++❗✔️:       while (true)
    // libstdc++❗✔️: 	{
    // libstdc++❗✔️: 	  while (__comp(__first, __pivot))
    // libstdc++❗✔️: 	    ++__first;
    // libstdc++❗✔️: 	  --__last;
    // libstdc++❗✔️: 	  while (__comp(__pivot, __last))
    // libstdc++❗✔️: 	    --__last;
    // libstdc++❗✔️: 	  if (!(__first < __last))
    // libstdc++❗✔️: 	    return __first;
    // libstdc++❗✔️: 	  std::iter_swap(__first, __last);
    // libstdc++❗✔️: 	  ++__first;
    // libstdc++❗✔️: 	}
    // libstdc++❗✔️:     }
    loop {
        while less(&v[first], &v[pivot]) {
            first += 1;
        }
        last -= 1;
        while less(&v[pivot], &v[last]) {
            last -= 1;
        }
        if first >= last {
            return first;
        }
        v.swap(first, last);
        first += 1;
    }
}

fn final_insertion<T, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    first: usize,
    last: usize,
    less: &mut F,
) {
    // libstdc++❗✔️: __final_insertion_sort(_RandomAccessIterator __first,
    // libstdc++❗✔️: 			   _RandomAccessIterator __last, _Compare __comp)
    // libstdc++❗✔️:     {
    // libstdc++❗✔️:       if (__last - __first > int(_S_threshold))
    // libstdc++❗✔️: 	{
    // libstdc++❗✔️: 	  std::__insertion_sort(__first, __first + int(_S_threshold), __comp);
    // libstdc++❗✔️: 	  std::__unguarded_insertion_sort(__first + int(_S_threshold), __last,
    // libstdc++❗✔️: 					  __comp);
    // libstdc++❗✔️: 	}
    // libstdc++❗✔️:       else
    // libstdc++❗✔️: 	std::__insertion_sort(__first, __last, __comp);
    // libstdc++❗✔️:     }
    if last - first > 16 {
        insertion(v, first, first + 16, less);
        for i in first + 16..last {
            linear_insert(v, i, less);
        }
    } else {
        insertion(v, first, last, less);
    }
}

fn insertion<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], first: usize, last: usize, less: &mut F) {
    // libstdc++❗✔️: __insertion_sort(_RandomAccessIterator __first,
    // libstdc++❗✔️: 		     _RandomAccessIterator __last, _Compare __comp)
    // libstdc++❗✔️:     {
    // libstdc++❗✔️:       if (__first == __last) return;
    // libstdc++❗✔️:
    // libstdc++❗✔️:       for (_RandomAccessIterator __i = __first + 1; __i != __last; ++__i)
    // libstdc++❗✔️: 	{
    // libstdc++❗✔️: 	  if (__comp(__i, __first))
    // libstdc++❗✔️: 	    {
    // libstdc++❗✔️: 	      typename iterator_traits<_RandomAccessIterator>::value_type
    // libstdc++❗✔️: 		__val = _GLIBCXX_MOVE(*__i);
    // libstdc++❗✔️: 	      _GLIBCXX_MOVE_BACKWARD3(__first, __i, __i + 1);
    // libstdc++❗✔️: 	      *__first = _GLIBCXX_MOVE(__val);
    // libstdc++❗✔️: 	    }
    // libstdc++❗✔️: 	  else
    // libstdc++❗✔️: 	    std::__unguarded_linear_insert(__i,
    // libstdc++❗✔️: 				__gnu_cxx::__ops::__val_comp_iter(__comp));
    // libstdc++❗✔️: 	}
    // libstdc++❗✔️:     }
    if first == last {
        return;
    }
    for i in first + 1..last {
        if less(&v[i], &v[first]) {
            v[first..=i].rotate_right(1);
        } else {
            linear_insert(v, i, less);
        }
    }
}

fn linear_insert<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], mut last: usize, less: &mut F) {
    // libstdc++❗✔️: __unguarded_linear_insert(_RandomAccessIterator __last,
    // libstdc++❗✔️: 			      _Compare __comp)
    // libstdc++❗✔️:     {
    // libstdc++❗✔️:       typename iterator_traits<_RandomAccessIterator>::value_type
    // libstdc++❗✔️: 	__val = _GLIBCXX_MOVE(*__last);
    // libstdc++❗✔️:       _RandomAccessIterator __next = __last;
    // libstdc++❗✔️:       --__next;
    // libstdc++❗✔️:       while (__comp(__val, __next))
    // libstdc++❗✔️: 	{
    // libstdc++❗✔️: 	  *__last = _GLIBCXX_MOVE(*__next);
    // libstdc++❗✔️: 	  __last = __next;
    // libstdc++❗✔️: 	  --__next;
    // libstdc++❗✔️: 	}
    // libstdc++❗✔️:       *__last = _GLIBCXX_MOVE(__val);
    // libstdc++❗✔️:     }
    while less(&v[last], &v[last - 1]) {
        v.swap(last, last - 1);
        last -= 1;
    }
}

fn heap_sort<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], less: &mut F) {
    // libstdc++❗✔️: __partial_sort(_RandomAccessIterator __first,
    // libstdc++❗✔️: 		   _RandomAccessIterator __middle,
    // libstdc++❗✔️: 		   _RandomAccessIterator __last,
    // libstdc++❗✔️: 		   _Compare __comp)
    // libstdc++❗✔️:     {
    // libstdc++❗✔️:       std::__heap_select(__first, __middle, __last, __comp);
    // libstdc++❗✔️:       std::__sort_heap(__first, __middle, __comp);
    // libstdc++❗✔️:     }
    // The source partial_sort(first,last,last) selects the entire range.
    make_heap(v, less);
    sort_heap(v, less);
}

fn make_heap<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], less: &mut F) {
    // libstdc++❗✔️: __make_heap(_RandomAccessIterator __first, _RandomAccessIterator __last,
    // libstdc++❗✔️: 		_Compare& __comp)
    // libstdc++❗✔️:     {
    // libstdc++❗✔️:       typedef typename iterator_traits<_RandomAccessIterator>::value_type
    // libstdc++❗✔️: 	  _ValueType;
    // libstdc++❗✔️:       typedef typename iterator_traits<_RandomAccessIterator>::difference_type
    // libstdc++❗✔️: 	  _DistanceType;
    // libstdc++❗✔️:
    // libstdc++❗✔️:       if (__last - __first < 2)
    // libstdc++❗✔️: 	return;
    // libstdc++❗✔️:
    // libstdc++❗✔️:       const _DistanceType __len = __last - __first;
    // libstdc++❗✔️:       _DistanceType __parent = (__len - 2) / 2;
    // libstdc++❗✔️:       while (true)
    // libstdc++❗✔️: 	{
    // libstdc++❗✔️: 	  _ValueType __value = _GLIBCXX_MOVE(*(__first + __parent));
    // libstdc++❗✔️: 	  std::__adjust_heap(__first, __parent, __len, _GLIBCXX_MOVE(__value),
    // libstdc++❗✔️: 			     __comp);
    // libstdc++❗✔️: 	  if (__parent == 0)
    // libstdc++❗✔️: 	    return;
    // libstdc++❗✔️: 	  __parent--;
    // libstdc++❗✔️: 	}
    // libstdc++❗✔️:     }
    if v.len() < 2 {
        return;
    }
    let mut parent = (v.len() - 2) / 2;
    loop {
        adjust_heap(v, parent, v.len(), less);
        if parent == 0 {
            return;
        }
        parent -= 1;
    }
}

fn sort_heap<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], less: &mut F) {
    // libstdc++❗✔️: __sort_heap(_RandomAccessIterator __first, _RandomAccessIterator __last,
    // libstdc++❗✔️: 		_Compare& __comp)
    // libstdc++❗✔️:     {
    // libstdc++❗✔️:       while (__last - __first > 1)
    // libstdc++❗✔️: 	{
    // libstdc++❗✔️: 	  --__last;
    // libstdc++❗✔️: 	  std::__pop_heap(__first, __last, __last, __comp);
    // libstdc++❗✔️: 	}
    // libstdc++❗✔️:     }
    let mut last = v.len();
    while last > 1 {
        last -= 1;
        v.swap(0, last);
        adjust_heap(v, 0, last, less);
    }
}

fn adjust_heap<T, F: FnMut(&T, &T) -> bool>(
    v: &mut [T],
    mut hole: usize,
    len: usize,
    less: &mut F,
) {
    // libstdc++❗✔️: __adjust_heap(_RandomAccessIterator __first, _Distance __holeIndex,
    // libstdc++❗✔️: 		  _Distance __len, _Tp __value, _Compare __comp)
    // libstdc++❗✔️:     {
    // libstdc++❗✔️:       const _Distance __topIndex = __holeIndex;
    // libstdc++❗✔️:       _Distance __secondChild = __holeIndex;
    // libstdc++❗✔️:       while (__secondChild < (__len - 1) / 2)
    // libstdc++❗✔️: 	{
    // libstdc++❗✔️: 	  __secondChild = 2 * (__secondChild + 1);
    // libstdc++❗✔️: 	  if (__comp(__first + __secondChild,
    // libstdc++❗✔️: 		     __first + (__secondChild - 1)))
    // libstdc++❗✔️: 	    __secondChild--;
    // libstdc++❗✔️: 	  *(__first + __holeIndex) = _GLIBCXX_MOVE(*(__first + __secondChild));
    // libstdc++❗✔️: 	  __holeIndex = __secondChild;
    // libstdc++❗✔️: 	}
    // libstdc++❗✔️:       if ((__len & 1) == 0 && __secondChild == (__len - 2) / 2)
    // libstdc++❗✔️: 	{
    // libstdc++❗✔️: 	  __secondChild = 2 * (__secondChild + 1);
    // libstdc++❗✔️: 	  *(__first + __holeIndex) = _GLIBCXX_MOVE(*(__first
    // libstdc++❗✔️: 						     + (__secondChild - 1)));
    // libstdc++❗✔️: 	  __holeIndex = __secondChild - 1;
    // libstdc++❗✔️: 	}
    // libstdc++❗✔️:       __decltype(__gnu_cxx::__ops::__iter_comp_val(_GLIBCXX_MOVE(__comp)))
    // libstdc++❗✔️: 	__cmp(_GLIBCXX_MOVE(__comp));
    // libstdc++❗✔️:       std::__push_heap(__first, __holeIndex, __topIndex,
    // libstdc++❗✔️: 		       _GLIBCXX_MOVE(__value), __cmp);
    // libstdc++❗✔️:     }
    // Swaps carry the saved source value down the hole, avoiding Clone/unsafe.
    let top = hole;
    let mut child = hole;
    while child < (len - 1) / 2 {
        child = 2 * (child + 1);
        if less(&v[child], &v[child - 1]) {
            child -= 1;
        }
        v.swap(hole, child);
        hole = child;
    }
    if len % 2 == 0 && child == (len - 2) / 2 {
        child = 2 * (child + 1);
        v.swap(hole, child - 1);
        hole = child - 1;
    }
    push_heap(v, hole, top, less);
}

fn push_heap<T, F: FnMut(&T, &T) -> bool>(v: &mut [T], mut hole: usize, top: usize, less: &mut F) {
    // libstdc++❗✔️: __push_heap(_RandomAccessIterator __first,
    // libstdc++❗✔️: 		_Distance __holeIndex, _Distance __topIndex, _Tp __value,
    // libstdc++❗✔️: 		_Compare& __comp)
    // libstdc++❗✔️:     {
    // libstdc++❗✔️:       _Distance __parent = (__holeIndex - 1) / 2;
    // libstdc++❗✔️:       while (__holeIndex > __topIndex && __comp(__first + __parent, __value))
    // libstdc++❗✔️: 	{
    // libstdc++❗✔️: 	  *(__first + __holeIndex) = _GLIBCXX_MOVE(*(__first + __parent));
    // libstdc++❗✔️: 	  __holeIndex = __parent;
    // libstdc++❗✔️: 	  __parent = (__holeIndex - 1) / 2;
    // libstdc++❗✔️: 	}
    // libstdc++❗✔️:       *(__first + __holeIndex) = _GLIBCXX_MOVE(__value);
    // libstdc++❗✔️:     }
    while hole > top {
        let parent = (hole - 1) / 2;
        if !less(&v[parent], &v[hole]) {
            break;
        }
        v.swap(hole, parent);
        hole = parent;
    }
}

#[cfg(test)]
mod proposals {
    use super::*;
    #[test]
    fn io44_equal_score_order_matches_cpp_row3500() {
        let scores: Vec<i32> = (0..33)
            .map(|i| if (24..32).contains(&i) { -3 } else { 100 })
            .collect();
        let mut indices: Vec<_> = (0..33).collect();
        sort_by(&mut indices, |a, b| scores[*a] < scores[*b]);
        assert_eq!(
            indices,
            vec![
                31, 30, 29, 28, 27, 26, 25, 24, 17, 16, 18, 19, 20, 21, 22, 23, 32, 0, 15, 14, 13,
                12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1
            ]
        );
    }
}
