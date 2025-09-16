//! Place a sequence of iterators into a *heap of iterators* (HIT); the
//! result is a new iterator that returns items in sorted order, *provided 
//! that the original iterators are sorted.*
//! 
//! This data structure has (we suspect) performance advantages over a more naive
//! implementation on a vector of iterators, that searches each entry in the vector
//! for the minimum value.  In particular, if the vector is very long one may need to
//! compare many elements.
//!
//!
//! # Three commands for 95% of uses
//!
//! Most users only need two functions to work with HIT's:
//! - [`hit_merge`] creates a HIT
//! - [`bulk_insert`](IteratorsMergedInSortedOrder::bulk_insert) modifies the HIT by inserting a new collection of iterators.
//! - [`hit_merge_by_predicate`] is both a function and (for convenience) a trait method, which allows
//! the user to specify their own order on entries.  This can be useful, for example,
//! if you wish to merge several iterators into a single iterator that returns values in *descending* order. 
//! 
//! The other items in the module are primarily just "internal machinery." 
//! 
//! # See also
//! 
//! OAT also offers specialized utilities for intersection of sets represented by sorted
//! iterators; see [`IntersectOrderedIterators`](crate::utilities::iterators::general).
//! 
//! # Examples
//! 
//! ```
//! // Import the definition of the `hit_merge` function; 
//! use oat_rust::utilities::iterators::merge::heap_of_iterators::hit_merge;
//! 
//! // The `bulk_insert` function is a method on the `IteratorsMergedInSortedOrder` struct; 
//! // so to import the function, we import the struct.
//! use oat_rust::utilities::iterators::merge::heap_of_iterators::IteratorsMergedInSortedOrder;
//! 
//! // Define two iterators
//! let iter_1 = vec![0, 2, ];
//! let iter_2 = vec![1, 3, ];
//! 
//! // Merge the iterators into a HIT
//! let mut hit = hit_merge( vec![ iter_1, iter_2 ] );
//! 
//! // The HIT returns entries from both iterators, in ascending order
//! assert_eq!( hit.next(),  Some( 0 ) );
//! assert_eq!( hit.next(),  Some( 1 ) );
//! assert_eq!( hit.next(),  Some( 2 ) );
//! assert_eq!( hit.next(),  Some( 3 ) );
//! 
//! // Now the HIT is empty, but we cah refill it with new iterators
//! let iter_3 = vec![0, 2, ];
//! let iter_4 = vec![1, 3, ];
//! hit.bulk_insert( vec![ iter_3, iter_4 ] );
//! 
//! assert_eq!( hit.next(),  Some( 0 ) );
//! assert_eq!( hit.next(),  Some( 1 ) );
//! assert_eq!( hit.next(),  Some( 2 ) );
//! assert_eq!( hit.next(),  Some( 3 ) );
//! ```
//! 
//! # Adapted from itertools
//! 
//! This module is simlar to (and adapted from) the `kmerge_by` module from 
//! itertools.  The key difference is that OAT allows you to *modify* a
//! HIT by adding new iterators, after it has been created.



use crate::utilities::heaps::heap::{ heapify, heapify_tail, sift_down };
use crate::utilities::order::{JudgePartialOrder, OrderOperatorAuto, OrderOperatorAutoReverse, OrderOperatorByLessThan};



// ----------------------------------------------------------------------------
// COPIED MACROS    
// ----------------------------------------------------------------------------


// Implementation's internal macros

macro_rules! debug_fmt_fields {
    ($tyname:ident, $($($field:ident).+),*) => {
        fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
            f.debug_struct(stringify!($tyname))
                $(
              .field(stringify!($($field).+), &self.$($field).+)
              )*
              .finish()
        }
    }
}

macro_rules! clone_fields {
    ($($field:ident),*) => {
        fn clone(&self) -> Self {
            Self {
                $($field: self.$field.clone(),)*
            }
        }
    }
}


// ----------------------------------------------------------------------------
// COPIED SIZE HINT FUNCTIONS (from other itertools file)
// ----------------------------------------------------------------------------


/// Add **x** correctly to a **SizeHint**.
#[inline]
fn hacked_size_hint_add_scalar(sh: (usize, Option<usize>), x: usize) -> (usize, Option<usize>)
{
    let (mut low, mut hi) = sh;
    low = low.saturating_add(x);
    hi = hi.and_then(|elt| elt.checked_add(x));
    (low, hi)
}

/// Add **SizeHint** correctly.
#[inline]
fn hacked_size_hint_add(a: (usize, Option<usize>), b: (usize, Option<usize>)) -> (usize, Option<usize>)
{
    let min = a.0.checked_add(b.0).unwrap_or(usize::MAX);
    let max = match (a.1, b.1) {
        (Some(x), Some(y)) => x.checked_add(y),
        _ => None,
    };

    (min, max)
}


// ----------------------------------------------------------------------------
// COPIED ORIGINAL FILE CONTENT
// ----------------------------------------------------------------------------


use std::mem::replace;
use std::fmt;

use super::super::general::PeekUnqualified;



//  HEAD/TAIL 
//  ---------------------------------------------------------------------------

/// Iterator wrapper having two fields: `head` (an item)  and `tail` (an iterator). 
///
/// If an iterator represents a sequence of elements, then `head` corresponds
/// to the first element, and `tail` corresponds to the sequence of all
/// remaining elements.
///
/// `PartialEq`, `Eq`, `PartialOrd` and `Ord` are implemented by comparing sequences based on
/// first items (which are guaranteed to exist).
///
#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct HeadTailHit<I>
    where I: Iterator
{
    pub head: I::Item,
    pub tail: I,
}

impl<I> HeadTailHit<I>
    where I: Iterator
{
    /// Constructs a `HeadTailHit` from an `Iterator`. Returns `None` if the `Iterator` is empty.
    pub fn new(mut it: I) -> Option<HeadTailHit<I>> {
        let head = it.next();
        head.map(|h| {
            HeadTailHit {
                head: h,
                tail: it,
            }
        })
    }

    /// Get the next element and update `head`, returning the old head in `Some`.
    ///
    /// Returns `None` when the tail is exhausted (only `head` then remains).
    fn next(&mut self) -> Option<I::Item> {
        if let Some(next) = self.tail.next() {
            Some(replace(&mut self.head, next))
        } else {
            None
        }
    }

    // ADDAPTED FROM:
    // /// Hints at the size of the sequence, same as the `Iterator` method.
    // fn size_hint(&self) -> (usize, Option<usize>) {
    //     size_hint::add_scalar(self.tail.size_hint(), 1)
    // }

    /// Hints at the size of the sequence, same as the `Iterator` method.
    fn size_hint(&self) -> (usize, Option<usize>) {
        hacked_size_hint_add_scalar(self.tail.size_hint(), 1)
    }

}


impl<I> Clone for HeadTailHit<I>
    where I: Iterator + Clone,
          I::Item: Clone
{
    clone_fields!(head, tail);
}


//  HEAP
//  ---------------------------------------------------------------------------

//  NB: this was original content of the file; we extracted and moved to a
//      separate file/folder dedicated to heaps.


// /// Make `data` a heap (min-heap w.r.t the sorting).
// fn heapify<T, S>(data: &mut [T], mut less_than: S)
//     where S: FnMut(&T, &T) -> bool
// {
//     for i in (0..data.len() / 2).rev() {
//         sift_down(data, i, &mut less_than);
//     }
// }
// 
// /// Sift down element at `index` (`heap` is a min-heap wrt the ordering)
// pub fn sift_down<T, S>(heap: &mut [T], index: usize, mut less_than: S)
//     where S: FnMut(&T, &T) -> bool
// {
//     debug_assert!(index <= heap.len());
//     let mut pos = index;
//     let mut child = 2 * pos + 1;
//     // the `pos` conditional is to avoid a bounds check
//     while pos < heap.len() && child < heap.len() {
//         let right = child + 1;
// 
//         // pick the smaller of the two children
//         if right < heap.len() && less_than(&heap[right], &heap[child]) {
//             child = right;
//         }
// 
//         // sift down is done if we are already in order
//         if !less_than(&heap[child], &heap[pos]) {
//             return;
//         }
//         heap.swap(pos, child);
//         pos = child;
//         child = 2 * pos + 1;
//     }
// }


//  IteratorsMergedInSortedOrder object
//  ---------------------------------------------------------------------------


/// An iterator adaptor that merges an abitrary number of base iterators
/// according to an ordering function.
///
/// Iterator element type is `I::Item`.
///
/// See [`hit_merge_by_fnmut`] for more information.
#[must_use = "iterator adaptors are lazy and do nothing unless consumed"]
pub struct IteratorsMergedInSortedOrder<I, F>
    where I: Iterator,
{
    heap: Vec<HeadTailHit<I>>,
    less_than: F,
}

impl <I, F>  IteratorsMergedInSortedOrder < I, F > 
    where 
        I: Iterator,
        F: JudgePartialOrder<  <I as IntoIterator>::Item >
{    
    /// Create an empty `IteratorsMergedInSortedOrder`.
    pub fn new( less_than: F ) -> IteratorsMergedInSortedOrder< I, F > { IteratorsMergedInSortedOrder{ heap: Vec::new(), less_than: less_than } }

    /// Returns `true` iff `self.next() = None`.
    /// 
    /// For details, see the implementation of the `Iterator` trait for [`IteratorsMergedInSortedOrder`].
    /// 
    /// ```
    /// use oat_rust::utilities::iterators::merge::heap_of_iterators::{hit_merge};
    /// 
    /// // Create a HIT iterator, and pop off some elements
    /// let mut hit = hit_merge( vec![ vec![1, 2], vec![0, 3] ] );
    /// assert!( ! hit.is_empty() );
    /// for step in 0..4 { let _ = hit.next(); }
    /// assert!( hit.is_empty() );
    /// 
    /// let mut hit = hit_merge( vec![ Vec::<usize>::new() ] );
    /// assert!( hit.is_empty() );
    /// ```    
    pub fn is_empty( &self ) -> bool { self.heap.is_empty() }

    /// Get the number of iterators in the heap.
    pub fn len( &self ) -> usize { self.heap.len() }

    /// Clears the heap, removing all iterators.
    /// 
    /// Concretely, this is achieved by calling `self.heap.clear()` on the vector `self.heap`.
    /// Thus it does not change the capacity of the vector that stores the heap.
    pub fn clear( &mut self )  { self.heap.clear() }

    /// Append a new seqeunce of iterators to the merge heap.
    ///
    /// ```
    /// use oat_rust::utilities::iterators::merge::heap_of_iterators::{hit_merge};
    /// 
    /// // Create a HIT iterator, and pop off some elements
    /// let ordered_sequences = vec![ vec![1, 2], vec![0, 3] ];
    /// let mut hit = hit_merge( ordered_sequences );
    /// assert_eq!( Some(0), hit.next() );
    /// assert_eq!( Some(1), hit.next() );
    /// 
    /// // Insert new iterators into the heap.
    /// hit.bulk_insert( vec![ vec![ 4, 5 ], vec![ 6 ] ] );
    /// let vec : Vec<usize> = hit.collect();
    /// assert_eq!( vec, vec![ 2, 3, 4, 5, 6 ] )
    /// ```    
    pub fn bulk_insert< J >( &mut self, iterables_to_insert: J )
        where 
            J:          IntoIterator,
            J::Item:    IntoIterator< IntoIter = I >,            
    {   
        // this is where we'll start the bulk heapify
        let tail_base = self.heap.len();     

        // push the new iterators onto the heap
        let iter = iterables_to_insert.into_iter();
        self.heap.extend(iter.filter_map(|it| HeadTailHit::new(it.into_iter())));
        
        // heapify
        let less_than = &mut self.less_than;
        heapify_tail(&mut self.heap, |a, b| less_than.judge_lt( &a.head, &b.head),
        & tail_base);        
    }  

    /// Append a new seqeunce of iterators to the merge heap.
    ///
    /// ```
    /// use oat_rust::utilities::iterators::merge::heap_of_iterators::{hit_merge};
    /// 
    /// // Create a HIT iterator, and pop off some elements
    /// let ordered_sequences = vec![ vec![1, 3], vec![0, 4] ];
    /// let mut hit = hit_merge( ordered_sequences );
    /// assert_eq!( Some(0), hit.next() );
    /// assert_eq!( Some(1), hit.next() );
    /// 
    /// // Insert new iterators into the heap.
    /// hit.insert_one_iter( vec![ 2, 5 ] );
    /// let vec : Vec<usize> = hit.collect();
    /// assert_eq!( vec, vec![ 2, 3, 4, 5 ] )
    /// ```    
    pub fn insert_one_iter< J >( &mut self, iterable: J )
        where 
            J:      IntoIterator< IntoIter = I >,            
    {   
        // wrap the item in a `Once` iterator
        let iterables_to_insert = std::iter::once( iterable );
        // call bulk_insert
        self.bulk_insert(iterables_to_insert);
    }      
 
}

impl<I, F> fmt::Debug for IteratorsMergedInSortedOrder<I, F>
    where I: Iterator + fmt::Debug,
          I::Item: fmt::Debug,
{
    debug_fmt_fields!(IteratorsMergedInSortedOrder, heap);
}


impl<I, F> Clone for IteratorsMergedInSortedOrder<I, F>
    where I: Iterator + Clone,
          I::Item: Clone,
          F: Clone,
{
    clone_fields!(heap, less_than);
}

impl<I, F> Iterator for IteratorsMergedInSortedOrder<I, F>
    where I: Iterator,
          F: JudgePartialOrder< I::Item>
{
    type Item = I::Item;

    fn next(&mut self) -> Option<Self::Item> {
        if self.heap.is_empty() {
            return None;
        }
        let result = if let Some(next) = self.heap[0].next() {
            next
        } else {
            self.heap.swap_remove(0).head
        };
        let less_than = &mut self.less_than;
        sift_down(&mut self.heap, 0, |a, b| less_than.judge_lt(&a.head, &b.head));
        Some(result)
    }

//    Original version, which uses private methods from itertools library.    
//    fn size_hint(&self) -> (usize, Option<usize>) {
//        self.heap.iter()
//                 .map(|i| i.size_hint())
//                 .fold1(size_hint::add)
//                 .unwrap_or((0, Some(0)))
//    }  
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.heap.iter()
                 .map(|i| i.size_hint())
                 .reduce(hacked_size_hint_add)
                 .unwrap_or((0, Some(0)))
    }  
}      

impl < I, F > PeekUnqualified for IteratorsMergedInSortedOrder< I, F > 
        where I: Iterator,
        F: JudgePartialOrder< I::Item>,
{
    fn peek_unqualified( &mut self ) -> Option< & <Self as Iterator>::Item > {
        match self.heap.is_empty() {
            true    =>  { None },
            false   =>  { Some( & self.heap[0].head ) }
        }        
    }
}





//  IteratorsMergedInSortedOrder makers
//  ---------------------------------------------------------------------------


/// Merge a sequence of iterators into a single iterator; result is sorted by
/// `less_than` if each iterator in the original sequences is sorted by 
/// `less_than`.
///
/// This is very similar to `hit_merge_by_fnmut` but the type constraint on 
/// the third parameter is different.  Moreover, `hit_merge_by_fnmut` is 
/// esseentially just a wrapper around this function.  We believe that
/// the original reason for splitting it into two had to do with 
/// engineering type constraints.
pub fn hit_merge_by_predicate<I, F>(iterable: I, less_than: F)
    -> IteratorsMergedInSortedOrder<<I::Item as IntoIterator>::IntoIter, F>
    where I: IntoIterator,
          I::Item: IntoIterator,
          F: JudgePartialOrder< <<I as IntoIterator>::Item as IntoIterator>::Item>,
{
    let iter = iterable.into_iter();
    let (lower, _) = iter.size_hint();
    let mut heap: Vec<_> = Vec::with_capacity(lower);
    heap.extend(iter.filter_map(|it| HeadTailHit::new(it.into_iter())));
    heapify(&mut heap, |a, b| less_than.judge_lt(&a.head, &b.head));
    IteratorsMergedInSortedOrder { heap, less_than }
}


/// Merge a sequence of iterators into a single iterator; result is sorted by
/// `less_than` if each iterator in the original sequence is sorted by 
/// `less_than`.

/// ```
/// use oat_rust::utilities::iterators::merge::heap_of_iterators::hit_merge_by_fnmut;
/// use num_traits::sign::Signed;
/// 
/// // Result may not respect the order function if the input sequences do not.
/// let unordered_sequences =  vec![ vec![ 2, 0, -4 ] ];
/// let x : Vec<_> = hit_merge_by_fnmut( unordered_sequences, |a, b| &a.abs() < &b.abs() ).collect();
/// assert_eq!( x , vec![ 2, 0, -4 ] );
///
/// // Result *will* respect the order function if the input sequences do.
/// let ordered_sequences = vec![ vec![1, -2], vec![0, -3] ];
/// let y : Vec<_> = hit_merge_by_fnmut( ordered_sequences, |a, b| &a.abs() < &b.abs() ).collect();
/// assert_eq!( y, vec![ 0, 1, -2, -3 ] )
/// ```
pub fn hit_merge_by_fnmut<I, F>(iter: I, less_than: F)
    -> IteratorsMergedInSortedOrder<<I::Item as IntoIterator>::IntoIter, OrderOperatorByLessThan< F, <I::Item as IntoIterator>::Item > >
    where I: Sized + IntoIterator,
          I::Item: IntoIterator,
          F: Fn(&<I::Item as IntoIterator>::Item,
                &<I::Item as IntoIterator>::Item) -> bool
{
    let order_comparator: OrderOperatorByLessThan<F, <I::Item as IntoIterator>::Item > = OrderOperatorByLessThan::new( less_than );
    hit_merge_by_predicate(iter, order_comparator )
}


/// Merge a sequence of iterators into a single iterator; result is sorted in
/// ascending order 
/// if each iterator in the original sequence is sorted in ascending order.
///
/// ```
/// use oat_rust::utilities::iterators::merge::heap_of_iterators::hit_merge;
/// 
/// // Result may not respect order if the input sequences do not.
/// let data_unordered = vec![ vec![2, 0, 4]];
/// let x : Vec<usize> = hit_merge( data_unordered ).collect();
/// assert_eq!( x, vec![ 2, 0, 4 ] );
///
/// let data_ordered = vec![ vec![1, 2], vec![0, 3] ];
/// let y : Vec<usize> = hit_merge( data_ordered ).collect();
/// assert_eq!( y, vec![ 0, 1, 2, 3 ] )
/// ```
pub fn hit_merge<I>(iterable: I) 
    -> IteratorsMergedInSortedOrder<<I::Item as IntoIterator>::IntoIter, OrderOperatorAuto>
    
    where I: IntoIterator,
          I::Item: IntoIterator,
          <<I as IntoIterator>::Item as IntoIterator>::Item: PartialOrd
{
    hit_merge_by_predicate(iterable, OrderOperatorAuto)
}

/// Merge a sequence of iterators into a single iterator; result is sorted in
/// descending order 
/// if each iterator in the original sequence is sorted in descending order.
///
/// ```
/// use oat_rust::utilities::iterators::merge::heap_of_iterators::hit_merge_descend;
/// 
/// // Result may not respect order if the input sequences do not.
/// let data_unordered = vec![ vec![2, 0, 4]];
/// let merged_unordered : Vec<usize> = hit_merge_descend( data_unordered ).collect();
/// assert_eq!( merged_unordered, vec![ 2, 0, 4 ] );
///
/// let data_ordered = vec![ vec![6, 4], vec![5, 3] ];
/// let merged_ordered : Vec<usize> = hit_merge_descend( data_ordered ).collect();
/// assert_eq!( merged_ordered, vec![ 6, 5, 4, 3 ] )
/// ```
pub fn hit_merge_descend<I>(iterable: I) 
    -> IteratorsMergedInSortedOrder<<I::Item as IntoIterator>::IntoIter, OrderOperatorAutoReverse>
    
    where I: IntoIterator,
          I::Item: IntoIterator,
          <<I as IntoIterator>::Item as IntoIterator>::Item: PartialOrd
{
    hit_merge_by_predicate(iterable, OrderOperatorAutoReverse)
}


//  ---------------------------------------------------------------------------
//  NEW CODE: MODIFY HEAP POST-HOC
//  ---------------------------------------------------------------------------


/// Append a new seqeunce of iterators to the merge heap.
///
/// ```
/// use oat_rust::utilities::iterators::merge::heap_of_iterators::{hit_merge, hit_bulk_insert};
/// 
/// // Create a HIT iterator, and pop off some elements
/// let ordered_sequences = vec![ vec![1, 2], vec![0, 3] ];
/// let mut hit = hit_merge( ordered_sequences );
/// assert_eq!( Some(0), hit.next() );
/// assert_eq!( Some(1), hit.next() );
/// 
/// // Insert new iterators into the heap.
/// hit_bulk_insert( &mut hit, vec![ vec![ 4, 5 ], vec![ 6 ] ] );
/// let vec : Vec<usize> = hit.collect();
/// assert_eq!( vec, vec![ 2, 3, 4, 5, 6 ] )
/// ```
pub fn hit_bulk_insert< I, F >( 
            merged : &mut IteratorsMergedInSortedOrder<<I::Item as IntoIterator>::IntoIter, F>, 
            iterable: I,
        )
    where I: IntoIterator,
          I::Item: IntoIterator,
          F: JudgePartialOrder< <<I as IntoIterator>::Item as IntoIterator>::Item>
{

    merged.bulk_insert( iterable )

}






pub trait HitMergeByPredicateTrait{
    /// Merge an iterable of iterators
    /// 
    /// Equivalent to `hit_merge_by_predicate( self, less_than )`
    fn hit_merge_by_predicate< F >( self, less_than: F)
    -> IteratorsMergedInSortedOrder<<Self::Item as IntoIterator>::IntoIter, F>
    where Self: IntoIterator + Sized,
          Self::Item: IntoIterator,
          F: JudgePartialOrder< <<Self as IntoIterator>::Item as IntoIterator>::Item>
    {
        hit_merge_by_predicate( self, less_than )
    }
}

impl < T >HitMergeByPredicateTrait for

    T

    where
          T: IntoIterator + Sized,
          T::Item: IntoIterator,
{}






    
#[cfg(test)]
mod tests {
    use ordered_float::OrderedFloat;

    use super::*;
    use crate::{utilities::{iterators::general::PeekUnqualified, }, rings::types::field_prime_order::{PrimeOrderField}, vectors::operations::{Transform}};    

    /// Verify that cloning works.
    #[test]
    fn test_heap_clone() {
        let ordered_sequences = vec![ vec![1, -2], vec![0, -3] ];
        let merged = hit_merge_by_fnmut( ordered_sequences, |a: &i64, b: &i64| &a.abs() < &b.abs() );
        let _ = merged.clone();
    }

    #[test]
    fn test_heap_functions() {

        // use oat_rust::utilities::iterators::merge::heap_of_iterators::hit_merge_by_fnmut;
        // use num_traits::sign::Signed;

        
        let ordered_sequences = vec![ vec![1, -2], vec![0, -3] ];
        let mut merged = hit_merge_by_fnmut( ordered_sequences, |a: &i64, b: &i64| &a.abs() < &b.abs() );
        while let Some( x ) = merged.peek_unqualified() {
            assert_eq!( Some(x.clone()), merged.next() )
        }

        let ordered_sequences: Vec<Vec<usize>> = vec![ vec![1, 2, 3], vec![] ];
        let mut merged = hit_merge_by_fnmut( ordered_sequences, |a: &usize, b: &usize| &a < &b );
        while let Some( x ) = merged.peek_unqualified() {
            assert_eq!( Some(x.clone()), merged.next() )
        }    
        
        let ordered_sequences: Vec<Vec<i64>> = vec![ vec![], vec![] ];
        let mut merged = hit_merge_by_fnmut( ordered_sequences, |a: &i64, b: &i64| &a.abs() < &b.abs() );
        while let Some( x ) = merged.peek_unqualified() {
            assert_eq!( Some(x.clone()), merged.next() )
        }            
    }

    #[test]    
    fn test_heap_allocations() {
        let nvertices = 800;
        // let ordered_sequence: Vec<_> = (0..nvertices).map(|x| (  Simplex{ filvalue: OrderedFloat( x as f64 ), vertices: vec![x,x+1,x+2] }  ) ).collect();

        let ring_operator   =   PrimeOrderField::new(3);        
        let ordered_sequences: Vec< Vec< _ > > = vec![ 
                (0..nvertices)
                    .map(
                        |x|  
                        ( 
                            ( OrderedFloat( x as f64 ), vec![x,x+1,x+2] ),
                            1usize
                        )
                    ).collect() 
            ; nvertices as usize ];
        let input_1 = ordered_sequences.iter().map(|x| x.iter().cloned().peekable());
        let _input_2 = ordered_sequences.iter().map(|x| x.iter().cloned().peekable());

        let mut h = hit_merge_by_predicate( input_1, OrderOperatorAuto  ).simplify( ring_operator.clone() );
        
        for vec_num in 0 .. nvertices {
            hit_bulk_insert( &mut h.unsimplified, vec![ ordered_sequences[ vec_num as usize ].iter().cloned().peekable() ] );
            let _ = h.next();
        }

        let input_2 = ordered_sequences.iter().map(|x| x.iter().cloned().peekable());        
        hit_bulk_insert( &mut h.unsimplified, input_2 );


        println!("henry");
    }

}