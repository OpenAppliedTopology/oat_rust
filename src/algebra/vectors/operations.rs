//! Operations on sparse vectors: [add](crate::algebra::vectors::operations::VectorOperations::add), [subtract](crate::algebra::vectors::operations::VectorOperations::subtract), [scale](crate::algebra::vectors::operations::VectorOperations::scale), [simplify](crate::algebra::vectors::operations::VectorOperations::add), [combine](crate::algebra::vectors::operations::Combine::combine), etc.
//!
//! 
//! # Convenient methods
//! 
//! The [VectorOperations] trait provides convenient methods for many of the most common vector operations, including
//! - [linear combination](crate::algebra::vectors::operations::VectorOperations::combine), including [sums](crate::algebra::vectors::operations::VectorOperations::sum) of any number of vectors
//! - [pairwise addition](crate::algebra::vectors::operations::VectorOperations::add) and [subtraction](crate::algebra::vectors::operations::VectorOperations::subtract)
//! - [scaling](crate::algebra::vectors::operations::VectorOperations::scale)
//! 
//! For a **full list of available methods**, see the lefthand menu bar for the [VectorOperations] page.
//!
//! # Additional functionality
//! 
//! There are also structs that perform operations not provided by the traits, see below for a complete list.
//! 
//! 
//! 
// //! # Examples
// //! 
// //! 
// //! Let `iter_a`, `iter_b`, and `iter_c` be sparse vectors, i.e. iterators that run over 
// //! sparse matrix entries.  For example, we could define `iter_a`, `iter_b`, `iter_c` as follows
// //! 
// //! ```
// //! // First define the entries
// //! // Note that `vec!` creates a standard rust vector, which is different 
// //! // from the sort of vector we care about)
// //! let entries_a   =   vec![ (1, 1.), (4, 4.) ];
// //! let entries_b   =   vec![ (2, 2.), (3, 3.) ];
// //! let entries_c   =   vec![ (1, 1.), (2, 2.), (3, 3.), (3, 3.), (4, 0.) ];
// //! 
// //! // Now define the sparse vector iterators
// //! // Note that `iter()` creates an interator, and `cloned()` reformats 
// //! // the entries of each iterator.
// //! let iter_a      =   entries_a.iter().cloned(); 
// //! let iter_b      =   entries_b.iter().cloned();
// //! let iter_c      =   entries_c.iter().cloned();
// //! ```
// //! 
// //! Let's also define the operator of the coefficient ring we want to work with.
// //!     
// //! ```
// //! // Load the module that allows us to define our coefficient ring.
// //! use oat_rust::algebra::rings::types::native::*;
// 
// //! // Define the operator of a coefficient ring (in this case, floating point real numbers)
// //! let ring_operator = RingOperatorForNativeRustNumberType::<f64>::new();   
// //! ```
// //!     
// //! * **Scale, drop zeros, gather, simplify**
// //! 
// //!     We can scale, drop zero entries, and gather terms as follows
// //!     
// //!     ```
// //!     use oat_rust::algebra::vectors::operations::*;
// //!     use oat_rust::algebra::rings::types::native::*;
// //! 
// //!     # // Define the vector
// //!     # let entries_a   =   vec![ (1, 1.), (4, 4.) ];
// //!     # let entries_c   =   vec![ (1, 1.), (2, 2.), (3, 3.), (3, 3.), (4, 0.) ];
// //!     # let iter_a      =   entries_a.iter().cloned();
// //!     # let iter_c      =   entries_c.iter().cloned();
// //!     #
// //!     # // Define the operator of the coefficient ring (in this case, floating point real numbers)
// //!     # let ring_operator = RingOperatorForNativeRustNumberType::<f64>::new();        
// //!       
// //!     // SCALE A VECTOR BY 2.
// //!     // Example: convert [ (1, 1.), (4, 4.) ] into [ (1, 2.), (4, 8.) ]
// //!     let scaled : Vec<_> = iter_a
// //!                             .clone() // this makes a copy of the iterator, so the original stays unchanged
// //!                             .scale( 2., ring_operator.clone() )
// //!                             .collect(); // this collects the entries of the iterator into a standard Rust vector
// //!     assert_eq!( scaled, vec![ (1, 2.), (4, 8.) ]);
// //!       
// //!     // DROP ZERO ENTRIES
// //!     // Example: convert [ (1, 1.), (2, 2.), (3, 3.), (3, 3.), (4, 0.) ] into [ (1, 1.), (2, 2.), (3, 3.), (3, 3.) ]
// //!     let dropped : Vec<_> = iter_c
// //!                             .clone() // this makes a copy of the iterator, so the original stays unchanged
// //!                             .drop_zeros( ring_operator.clone() )
// //!                             .collect(); // this collects the entries of the iterator into a standard Rust vector
// //!     
// //!     assert_eq!( dropped, vec![ (1, 1.), (2, 2.), (3, 3.), (3, 3.) ]);
// //!       
// //!     // GATHER CONSECUTIVE ENTRIES THAT SHARE THE SAME INDEX
// //!     // The resulting vector has no repeating consecutive indices; each index gets 
// //!     // the sum of the corresponding coefficients.
// //!     // Example: convert [(1,1.), (1,0.5), (2,0.), (1,0.)] into [(1,1.5), (2,0.), (1,0.)]
// //!     let gathered : Vec<_> = iter_c
// //!                             .clone() // this makes a copy of the iterator, so the original stays unchanged
// //!                             .peekable() // this puts the iterator in a slightly different form, which is compatible with gather
// //!                             .gather( ring_operator.clone() )
// //!                             .collect(); // this collects the entries of the iterator into a standard Rust vector
// //!     assert_eq!( gathered, vec![ (1, 1.), (2, 2.), (3, 6.), (4, 0.) ]);   
// //! 
// //!     // GATHER CONSECUTIVE ENTRIES THAT SHARE THE SAME INDEX, AND DROP RESULTING ZEROS
// //!     let simplified : Vec<_> = iter_c
// //!                                 .clone() // this makes a copy of the iterator, so the original stays unchanged
// //!                                 .peekable() // this puts the iterator in a slightly different form, which is compatible with gather
// //!                                 .simplify( ring_operator.clone() )
// //!                                 .collect(); // this collects the entries of the iterator into a standard Rust vector
// //!     assert_eq!( simplified, vec![ (1, 1.), (2, 2.), (3, 6.), ]);  
// //!     ```
// //! * **Combine iterators in sorted order** (basic)
// //! 
// //!   We can combine two iterators, `A` and `B`, into a single iterator `C` using the 
// //!     [merge](https://docs.rs/itertools/0.7.2/itertools/fn.merge.html) function from [itertools](https://docs.rs/itertools/latest/itertools/). 
// //!     The resulting iterator, `C`, is a 
// //!     [Merge struct](https://docs.rs/itertools/0.7.8/itertools/structs/struct.Merge.html).
// //!     Iterator `C` will iterate over all the entries in `A` and `B`.
// //!     If the items of `A` and `B` appear in sorted order, then the items of `C` will also 
// //!     appear in sorted order.
// //!     ```
// //!     # let entries_a   =   vec![ (1, 1.), (4, 4.) ];
// //!     # let entries_b   =   vec![ (2, 2.), (3, 3.) ];
// //!     # let iter_a      =   entries_a.iter().cloned(); 
// //!     # let iter_b      =   entries_b.iter().cloned();
// //! 
// //!     use itertools::merge;
// //!     use std::iter::FromIterator;
// //!     
// //!     // Merge [ (1, 1.), (4, 4.) ] and [ (2, 2.), (3, 3.) ].
// //!     // The entries in these vectors are in sorted order, so the resulting iterator will 
// //!     // also produce items in sorted order.
// //!     let iter_merged   =   merge(iter_a, iter_b);
// //!     let entries_mrgd  =   Vec::from_iter(iter_merged);
// //!     assert_eq!( entries_mrgd, vec![ (1, 1.), (2, 2.), (3, 3.), (4, 4.) ])
// //!     ```
// //! 
// //!     We can merge any `k` iterators of the same kind using the 
// //!     [merge](https://docs.rs/itertools/0.7.2/itertools/fn.merge.html) 
// //!     function from 
// //!     [itertools](https://docs.rs/itertools/latest/itertools/).  
// //! 
// //! * **Combine iterators in sorted order** (advanced)
// //!  
// //!     For advanced usage (eg matrix reduction), we also provide a
// //!     customized merge process in the [hit](crate::utilities::iterators::merge::hit) module.
// //! 
// //! * **Add**
// //! 
// //!     We can add the vectors represented by `iter_a` and `iter_b` by 
// //!     first combining (e.g., with the `merge` function discussed above), 
// //!     then applying the `gather` method.
// //! 
// //! * **Subtract**
// //! 
// //!     We can subtract  `iter_a` from `iter_b` by first scaling `iter_a` by `-1`, then adding.



use crate::algebra::matrices::query::{ MatrixAlgebra, MatrixOracle, };
use crate::algebra::matrices::operations::combine_rows_and_columns::{LinearCombinationOfColumns, LinearCombinationOfColumnsReverse, LinearCombinationOfRows, LinearCombinationOfRowsReverse};
use crate::utilities::functions::evaluate::EvaluateFunction;
use crate::utilities::iterators::general::{PeekUnqualified};
use crate::algebra::vectors::entries::{KeyValGet, KeyValSet, KeyValNew};
use crate::algebra::rings::traits::{SemiringOperations, RingOperations};
use crate::utilities::iterators::merge::two_type::MergeTwoIteratorsByOrderOperator;
use crate::utilities::iterators::merge::hit::{IteratorsMergedInSortedOrder, HitMergeByPredicateTrait};
use crate::utilities::order::{JudgeOrder, JudgePartialOrder, ReverseOrder};
use crate::utilities::sets::MapHasKeyOrSequenceHasElement;

use std::collections::HashMap;
use std::fmt::{Debug};
use std::cmp::Ordering;
use std::hash::Hash;
use std::iter::Peekable;
use std::marker::PhantomData;


use derive_getters::{Getters, Dissolve};
use derive_new::new;



//  ---------------------------------------------------------------------------
//  ITERATOR WRAPPER / TRANSFORMATIONS 
//  ---------------------------------------------------------------------------



//  ---------------------------------------------------------------------------
//  Linear combinations


/// A wrapper struct representing a (simplified) linear combination of vectors.  If each
/// vector in the combination returns elements in (strictly or nonstrictly) ascending 
/// order of index, then this struct returns entries in **non-strictly** ascending order
/// of index (repeat indices may occurr).
/// 
/// This struct is intended simply to simplify the expression of a more complicated type
pub type LinearCombinationUnsimplified
            < SparseVector, RingOperator, OrderOperator >  
        =   
            IteratorsMergedInSortedOrder<   
                    // type parameter #1 = a scalar multiple of a vector
                    Scale<
                        SparseVector,
                        RingOperator,
                    >,
                    // type parameter #2 = the object that determines order on entries
                    OrderOperator,
                >;         


/// Type alias for a (simplified) linear combination of vectors.  If each
/// vector in the combination returns elements in (strictly or nonstrictly) ascending 
/// order of index, then this struct returns entries in **strictly** ascending order
/// of index (no repeat indices).
/// 
/// This struct is intended simply to simplify the expression of a more complicated type
pub type LinearCombinationSimplified
            < SparseVector, RingOperator, OrderOperator >  
        =   
            Simplify<
                    IteratorsMergedInSortedOrder<   
                            // type parameter #1 = a scalar multiple of a vector
                            Scale<
                                SparseVector,
                                RingOperator,
                            >,
                            // type parameter #2 = the object that determines order on entries
                            OrderOperator,
                        >,
                    RingOperator,
                >;   




//  ---------------------------------------------------------------------------
//  DROP ZEROS


/// Iterates over the same items as `self.undropped`, skipping any items with coefficient 0.
/// 
/// Formally, we skip any item `x` such that `self.ring.is_0( x.val() )==true`.
/// 
#[derive(new,Clone,Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct DropZeros 
            < SparseVector, RingOperator >             
{
    undropped:              SparseVector,
    ring_operator:          RingOperator,
}

impl    < SparseVector, RingOperator > 
        
        Iterator for 
        
        DropZeros 
            < SparseVector, RingOperator > 

    where       RingOperator:           SemiringOperations,
                SparseVector:           Iterator,
                SparseVector::Item:     KeyValGet< Val = RingOperator::Element >,

{
    type Item = SparseVector::Item;

    fn next( &mut self) -> Option< Self::Item > 
    {
        let mut next = self.undropped.next();

        while let Some( ref x ) = next {
            if self.ring_operator.is_0( x.val() ) { next = self.undropped.next(); }
            else {break} 
        }
        next 
    }
}

//  ---------------------------------------------------------------------------
//  SCALE    


/// Iterates over the same items as `self.unscaled`, with all coefficients scaled by `self.scalar`.
#[derive(new,Clone,Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct Scale
    < SparseVector, RingOperator > 

    where   RingOperator:               SemiringOperations,
{
    unscaled:                   SparseVector,
    ring_operator:              RingOperator,    
    scaling_coefficient:        RingOperator::Element,    
}


impl    < SparseVector, RingOperator, > 
        
    Iterator for 
    
    Scale
        < SparseVector, RingOperator, > 
   
    where   SparseVector:               Iterator,
            SparseVector::Item:         KeyValSet< Val = RingOperator::Element >,
            RingOperator:               SemiringOperations,

{
    type Item = SparseVector::Item;

    fn next( &mut self) -> Option< Self::Item > 
        {
            if let Some( mut x ) = self.unscaled.next() { 
                x.set_val( 
                    self.ring_operator.multiply( 
                        x.val(), 
                        self.scaling_coefficient.clone(), 
                    )
                );
                Some(x)
            }
            else { None }
        }
}


//  ---------------------------------------------------------------------------
//  GATHER COEFFICIENTS 


/// Iterates over the same items as `self.ungathered`, except that 
/// consecutive entries with equal indices are merged into a single entry whose
/// coefficient is the sum of the coefficients.
#[derive(new,Clone,Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct Gather
            < SparseVector, RingOperator > 

{
    ungathered:         SparseVector,
    ring_operator:      RingOperator,
}



impl    < SparseVector, RingOperator > 

        Iterator for Gather
    
        < SparseVector, RingOperator > 
   
    where   SparseVector:               PeekUnqualified + 
                                        Iterator<
                                            Item:       KeyValSet< 
                                                            Key:        PartialEq,
                                                            Val     =   RingOperator::Element 
                                                        >,
                                        >,
            RingOperator:               SemiringOperations,
{
    type Item = SparseVector::Item;

    fn next( &mut self) -> Option< Self::Item > 
    {
        if let Some( mut x ) = self.ungathered.next() {
            while let Some( peek ) = self.ungathered.peek_unqualified() {
                if peek.key() == x.key() { 
                    x.set_val(
                        self.ring_operator.add( 
                            x.val(), 
                            peek.val() 
                        )
                    );
                    let _ = self.ungathered.next(); // we have already gotten what we need
                }
                else { break }
            }
            Some( x )
        }
        else 
        { None }
    }
}


//  ---------------------------------------------------------------------------
//  SIMPLIFY


/// Combines consecutive elements with equal index into a single entry; drops this entry, if it is zero.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::rings::types::native::RingOperatorForNativeRustNumberType;
/// use oat_rust::algebra::vectors::operations::Simplify;
/// use itertools::Itertools;
/// 
/// // Define the operator of the coefficient ring_operator.
/// let ring_operator = RingOperatorForNativeRustNumberType::<f64>::new();  
/// 
/// // Simplify some iterators:
/// let vec: Vec<(i32,f64)> = vec![ (0, 0.), (0, 1.), (1, 1.), (1, -1.), (2, 1.)];
/// let unsimplified = vec.iter().cloned().peekable();
/// let simplified = Simplify::new( unsimplified, ring_operator.clone() );       
/// assert_eq!( simplified.collect_vec(), vec![ (0,1.), (2, 1.) ] );
/// 
/// let vec: Vec<(usize, f64)> = vec![];
/// let unsimplified = vec.iter().cloned().peekable();            
/// let simplified = Simplify::new( unsimplified, ring_operator.clone() );       
/// assert_eq!( simplified.collect_vec(), vec![] );
/// 
/// let vec = vec![ (0, 0.), (0, 0.) ];
/// let unsimplified = vec.iter().cloned().peekable();            
/// let simplified = Simplify::new( unsimplified, ring_operator.clone() );              
/// assert_eq!( simplified.collect_vec(), vec![] );
/// ```
#[derive(new,Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd,Hash)]
pub struct Simplify
    
    < SparseVector, RingOperator > 
{
    pub unsimplified:           SparseVector,
    pub ring_operator:          RingOperator,
}


impl    < SparseVector, RingOperator > 

        Iterator for 
        
        Simplify
            < SparseVector, RingOperator > 
   
    where   SparseVector:                                   Iterator<
                                                                Item:   KeyValSet< 
                                                                            Key:        PartialEq,
                                                                            Val     =   RingOperator::Element,
                                                                        >
                                                            >,
            SparseVector:                                   PeekUnqualified,
            RingOperator:                                   SemiringOperations,
{
    type Item = SparseVector::Item;

    fn next( &mut self) -> Option< Self::Item > 
    {
        while let Some( mut x ) = self.unsimplified.next() {

            while let Some( peek ) = self.unsimplified.peek_unqualified() {

                if peek.key() == x.key() { 
                    x.set_val(
                        self.ring_operator.add( 
                            x.val(), 
                            peek.val() 
                        )
                    );
                    let _ = self.unsimplified.next(); // we have already gotten what we need
                }
                else { break }
            }
            match self.ring_operator.is_0( x.val() ){
                false   => { return Some( x ) },    //  return the gathered entry, if its coefficient is nonzero                
                true    => { continue },         //  otherwise skip to the next entry
            }
        }
        // else 
        // { println!("returning a simplified term"); None }        
        // println!("returning a simplified term"); 
        None
    }
}


    
impl < SparseVector, RingOperator > 

    Clone for

    Simplify
        < SparseVector, RingOperator > 

    where   SparseVector:               Clone,
            RingOperator:               Clone,

{
    fn clone(&self) -> Self {
        Simplify { unsimplified: self.unsimplified.clone(), ring_operator: self.ring_operator.clone() }
    }
}










//  ---------------------------------------------------------------------------
//  NEGATE


/// Negates the coefficient of each entry, i.e., swaps `( i, a )` with `(i, -a)`.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::vectors::operations::Negate;
/// use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;
/// 
/// let vec = vec![ (0, 0usize), (1, 1usize) ];
/// let ring_operator = PrimeOrderField::new(7);
/// let negated = Negate::new( vec.iter().cloned(), ring_operator );
/// itertools::assert_equal( negated, vec![ (0, 0usize), (1, 6usize) ] );
/// ```
#[derive(new,Clone,Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct Negate 
            < SparseVector, RingOperator > 
{
    unnegated:              SparseVector,
    ring_operator:          RingOperator,   
}
          

impl    < SparseVector, RingOperator > 
        
        Iterator for 
        
        Negate 
            < SparseVector, RingOperator > 

    where   SparseVector:               Iterator,
            SparseVector::Item:         PartialEq + KeyValSet< Val = RingOperator::Element >,
            RingOperator:               RingOperations,

{
    type Item = SparseVector::Item;

    fn next( &mut self) -> Option< Self::Item > 
    {
        match self.unnegated.next() {
            Some( mut next_item ) => {
                next_item.set_val( 
                        self.ring_operator.negate( next_item.val() )
                    );
                Some( next_item )
            },
            None => { None }
        }
    }
}




//  ---------------------------------------------------------------------------
//  MASK



/// A "mask" over the iterator; it returns an index-coefficient pair iff the index is a key in the
/// gien hashmap.
/// 
/// # Design notes
/// 
/// It would be possible to construct an object with the same behavior via Rust's [filter](https://doc.rust-lang.org/stable/std/iter/trait.Iterator.html#method.filter)
/// method.  However, to achieve the desired result one needs access to the hashmap; the most natural way
/// to provide such acceess to a [`Filter`](https://doc.rust-lang.org/stable/std/iter/struct.Filter.html) object is
/// to use a closure.  However, closures give rise to a host of troubles vis-a-vis type inferrence.  The object
/// that we construct here is intended to serve as a basic primitive; a building block that will be used 
/// throughout the library and combined in complicated ways with other structs.  Therefore we prefer to define
/// a new object with better properties vis-a-vis type inference.

#[derive(new,Clone,Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct OnlyIndicesInsideCollection
                < SparseVector, CollectionToInclude, >
{
    entry_iter:                 SparseVector,  
    collection_to_include:      CollectionToInclude,  
}

impl < SparseVector, CollectionToInclude, > 
    
    Iterator for

    OnlyIndicesInsideCollection
        < SparseVector, CollectionToInclude > 

    where   SparseVector:               Iterator,
            SparseVector::Item:         KeyValGet,
            CollectionToInclude:        MapHasKeyOrSequenceHasElement< <SparseVector::Item as KeyValGet>::Key >,
{
    type Item = SparseVector::Item;

    fn next( &mut self ) -> Option< Self::Item > {
        // println!("OnlyIndicesInsideCollection: next(): BEGINNING");
        for entry in self.entry_iter.by_ref() {
            // println!("OnlyIndicesInsideCollection: next(): BEGINNING WHILE LOOP");
            // println!("IT APPEARS THAT THE ERROR OCCURS WHEN WE TRY TO CALL self.collection_to_include.map_has_key_or_sequence_has_element ");
            match self.collection_to_include.map_has_key_or_sequence_has_element( & entry.key() ) {
                true    =>  { return Some( entry ) },
                false   =>  { continue },
            }
        }
        // println!("OnlyIndicesInsideCollection: next(): ENDING");        
        None
    }
}


//  ---------------------------------------------------------------------------
//  OTHER MASK



/// A "mask" over an iterator; it returns an index-coefficient pair iff the index is **not** a key in the
/// given hashmap or an element of the given vector.
/// 
/// # Design notes
/// 
/// It would be possible to construct an object with the same behavior via Rust's [filter](https://doc.rust-lang.org/stable/std/iter/trait.Iterator.html#method.filter)
/// method.  However, to achieve the desired result one needs access to the hashmap; the most natural way
/// to provide such acceess to a [`Filter`](https://doc.rust-lang.org/stable/std/iter/struct.Filter.html) object is
/// to use a closure.  However, closures give rise to a host of troubles vis-a-vis type inferrence.  The object
/// that we construct here is intended to serve as a basic primitive; a building block that will be used 
/// throughout the library and combined in complicated ways with other structs.  Therefore we prefer to define
/// a new object with better properties vis-a-vis type inference.
#[derive(new,Clone,Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct OnlyIndicesOutsideCollection
            < SparseVector, CollectionToExclude >
{
    entry_iter:                 SparseVector,
    collection_to_exclude:      CollectionToExclude, 
}


impl < SparseVector, CollectionToExclude > 
    
    Iterator for

    OnlyIndicesOutsideCollection
        < SparseVector, CollectionToExclude > 

    where   SparseVector:               Iterator,
            SparseVector::Item:         KeyValGet,
            CollectionToExclude:        MapHasKeyOrSequenceHasElement< <SparseVector::Item as KeyValGet>::Key >,
{
    type Item = SparseVector::Item;

    fn next( &mut self ) -> Option< Self::Item > {
        for entry in self.entry_iter.by_ref() {
            match self.collection_to_exclude.map_has_key_or_sequence_has_element( & entry.key() ) {
                false   =>  return Some( entry ),
                true    =>  continue,
            }
        }
        None
    }
}


/// Transforms a tuple entry by remapping its index and val independently, then forming a new entry type
/// 
/// Concretely, this struct implements an instance of `EvaluateFunction` that performs this transformation.
#[derive(Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct RemapEntryTuple< KeyOld, KeyNew, KeyMap, ValOld, ValNew, ValMap, EntryNew > {
    pub keymap:         KeyMap,
    pub valmap:         ValMap,
    phantom_keyold:     PhantomData<KeyOld>,
    phantom_keynew:     PhantomData<KeyNew>,
    phantom_valold:     PhantomData<ValOld>,
    phantom_valnew:     PhantomData<ValNew>,  
    phantom_entrynew:   PhantomData<EntryNew>          
}

impl < KeyOld, KeyNew, KeyMap, ValOld, ValNew, ValMap, EntryNew >

    RemapEntryTuple
        < KeyOld, KeyNew, KeyMap, ValOld, ValNew, ValMap, EntryNew > 
{
    /// Form a new instance from an objection`keymap` that remaps keys, and an object `valmap` that remaps values.
    pub fn new( keymap: KeyMap, valmap: ValMap ) -> Self {
        RemapEntryTuple{ keymap, valmap, phantom_entrynew: PhantomData, phantom_keynew: PhantomData, phantom_keyold: PhantomData, phantom_valnew: PhantomData, phantom_valold: PhantomData}
    }
}


impl < KeyOld, KeyNew, KeyMap, ValOld, ValNew, ValMap, EntryNew >
    
    EvaluateFunction
        < (KeyOld,ValOld), EntryNew > for 

    RemapEntryTuple
        < KeyOld, KeyNew, KeyMap, ValOld, ValNew, ValMap, EntryNew > where 

    KeyMap:     EvaluateFunction< KeyOld, KeyNew >,
    ValMap:     EvaluateFunction< ValOld, ValNew >,
    EntryNew:   KeyValNew < Key = KeyNew, Val = ValNew >
{
    /// Transforms a tuple entry by remapping its index and val independently, then forming a new entry type
    fn evaluate_function( &self, input: (KeyOld,ValOld) ) -> EntryNew {
        EntryNew::new( self.keymap.evaluate_function(input.0), self.valmap.evaluate_function(input.1) )
    }
}


impl < KeyOld, KeyNew, KeyMap, ValOld, ValNew, ValMap, EntryNew >
    
    Clone for 

    RemapEntryTuple
        < KeyOld, KeyNew, KeyMap, ValOld, ValNew, ValMap, EntryNew > where 

    KeyMap:     Clone,
    ValMap:     Clone,
{
    /// Transforms a tuple entry by remapping its index and val independently, then forming a new entry type
    fn clone( &self) -> Self {
        RemapEntryTuple::new( self.keymap.clone(), self.valmap.clone() )
    }
}




//  ---------------------------------------------------------------------------
//  MAP BY INDEX


/// Similar to map for iterators, however the only thing we map is the index.
/// 
/// This object always iteratos over tuples.  If you want to interate over some other kind of object,
/// it may be worth considering the [`MapByTransform`](oat_rust::utilities::iterators::general) struct.
///
/// # Examples
/// 
/// ```
/// use std::collections::HashMap;
/// use oat_rust::algebra::vectors::operations::ChangeIndexSimple;
/// 
/// // create iterator
/// let entry_iter_data: Vec< (usize, usize) > = vec![(1,1), (2,2)];
/// let entry_iter = entry_iter_data.iter();
/// 
/// // hashmaps implement `EvaluateFunction` automatically; see documentaion for `EvaluateFunction` for details
/// let mut hash = HashMap::new();
/// hash.insert( 1usize, -1i32 );
/// hash.insert( 2usize, -2i32 );
/// 
/// // create a filter/mapped iterator        
/// let reindexed: ChangeIndexSimple<_,_,i32> = ChangeIndexSimple::new( entry_iter, &hash );
/// 
/// // check that it is correct        
/// itertools::assert_equal( reindexed, vec![(-1,1usize), (-2,2usize)] );
/// ```
#[derive(new,Clone,Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct ChangeIndexSimple
                < SparseVector, IndexChanger, NewIndex, >
{
    entry_iter:                 SparseVector,    
    index_changer:              IndexChanger, 
    phantom_indexnew:           PhantomData< NewIndex >,    
}


impl < SparseVector, IndexChanger, NewIndex, Coefficient > 
    
    Iterator for

    ChangeIndexSimple
        < SparseVector, IndexChanger, NewIndex, > 

    where   
            SparseVector:               Iterator,
            SparseVector::Item:         KeyValGet< Val = Coefficient >,
            IndexChanger:               EvaluateFunction< <SparseVector::Item as KeyValGet>::Key, NewIndex >,
{
    type Item = ( NewIndex, Coefficient );

    fn next( &mut self ) -> Option< Self::Item > {
        match self.entry_iter.next() {
            Some( entry ) => {
                Some(  (  self.index_changer.evaluate_function( entry.key() ), entry.val()  )  )
            },
            None => { None }
        }
    }
}



//  ---------------------------------------------------------------------------
//  CHANGE ENTRY TYPE


/// Wrapper around an entry iterator; the wrapper iterates over the same sequence of index-coefficient pairs, but changes the type of each entry to `EntryNew`
#[derive(new,Clone,Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct ChangeEntryType
                < SparseVector, EntryNew >                  
{
    entry_iter:                 SparseVector,
    phantom_entrynew:           PhantomData< EntryNew >,
}


impl < SparseVector, EntryNew, Index, Coefficient >
    
    Iterator for

    ChangeEntryType
        < SparseVector, EntryNew >

    where   SparseVector:               Iterator,
            SparseVector::Item:         KeyValGet< Key = Index, Val = Coefficient >,
            EntryNew:                   KeyValNew < 
                                            Key = Index, 
                                            Val = Coefficient 
                                        >,            
{
    type Item = EntryNew;

    fn next( &mut self ) -> Option< Self::Item > {
        self.entry_iter.next().map(|entry| EntryNew::new(  entry.key(), entry.val()  ))
    }
}


//  ---------------------------------------------------------------------------
//  FILTER CHANGE INDEX


/// Similar to filter-map for iterators, however the only thing we filter and map is the index.
///
/// # Examples
/// 
/// ```
/// use std::collections::HashMap;
/// use oat_rust::algebra::vectors::operations::{FilterChangeIndex, };
/// 
/// // create iterator
/// let entry_iter_data: Vec< (usize, usize) > = vec![(1,1), (2,2), (3,3)];
/// let entry_iter = entry_iter_data.iter();
/// 
/// // hashmaps implement `EvaluateFunction` automatically; see documentaion for `EvaluateFunction` for details
/// let mut hash = HashMap::new();
/// hash.insert( 1usize, -1i32 );
/// hash.insert( 2usize, -2i32 );
/// 
/// // create a filter/mapped iterator        
/// let iter_filtermapped_by_index: FilterChangeIndex<
///                                     std::slice::Iter<'_, (usize, usize)>,
///                                     HashMap<usize,i32>,
///                                     i32
///                                 > 
///     = FilterChangeIndex::new( entry_iter, hash );
/// 
/// // check that it is correct        
/// itertools::assert_equal( iter_filtermapped_by_index, vec![(-1,1usize), (-2,2usize)] );
/// ```
#[derive(new,Clone,Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct FilterChangeIndex
                < SparseVector, IndexFilterChanger, IndexNew, >
{
    entry_iter:                 SparseVector,    
    index_filter_changer:       IndexFilterChanger, 
    phantom_indexnew:           PhantomData< IndexNew >,    
}


impl < SparseVector, IndexFilterChanger, IndexNew, Coefficient, > 
    
    Iterator for

    FilterChangeIndex
        < SparseVector, IndexFilterChanger, IndexNew, > 

    where   SparseVector:               Iterator,
            SparseVector::Item:         KeyValGet< Val = Coefficient >,
            IndexFilterChanger:         EvaluateFunction< <SparseVector::Item as KeyValGet>::Key, Option<IndexNew> >,
{
    type Item = ( IndexNew, Coefficient );

    fn next( &mut self ) -> Option< Self::Item > {
        for entry in self.entry_iter.by_ref() {
            match self.index_filter_changer.evaluate_function( entry.key() ) {
                None    =>  return None,
                Some( index_new)    => return Some( (index_new, entry.val() ) ),
            }
        }
        None
    }
}



//  ---------------------------------------------------------------------------
//  ENTRYWISE PRODUCT


/// Entrywise product of two vectors
/// 
/// Concretely, if `u = [u0 .. um]` and `v = [v0 .. vm ]` are vectors, then this returns
/// a vector of form `[u0 * v0 .. um * vm]`.
/// 
/// **Assumes that the entries of both `u` and `v` are sorted in strictly ascending order of index, with respect to the [order operator](crate::utilities::order).**
/// 
/// # See also
/// 
/// The [VectorOperations::multiply_entrywise] method.
/// 
/// # Examples
/// 
/// See the [VectorOperations::multiply_entrywise] method.
#[derive(new,Clone,Copy,Debug,Dissolve,Eq,PartialEq,Getters,Ord,PartialOrd)]
pub struct EntrywiseProduct      
    
    < VectorA, VectorB, RingOperator, OrderOperator > 
    
    where   
            VectorA:            Iterator,
            VectorB:            Iterator,            
{
    vec_a:                  VectorA,
    vec_b:                  VectorB,
    next_a:                 Option< VectorA::Item >,
    next_b:                 Option< VectorB::Item >,
    ring_operator:          RingOperator,
    order_operator:         OrderOperator,
}

impl    < VectorA, VectorB, RingOperator, OrderOperator, Index > 
        
        Iterator for 
            
        EntrywiseProduct
        
        < VectorA, VectorB, RingOperator, OrderOperator > 
   
    where   RingOperator:       SemiringOperations,
            OrderOperator:      JudgeOrder< Index >,
            VectorA:            Iterator,
            VectorB:            Iterator,            
            VectorA::Item:      KeyValSet<
                                    Key     =   Index,
                                    Val     =   RingOperator::Element,
                                >,
            VectorB::Item:      KeyValGet<
                                    Key     =   Index,
                                    Val     =   RingOperator::Element,
                                >,                                

{
    type Item = VectorA::Item;

    fn next( &mut self) -> Option< Self::Item > 
        {
            loop {
                if self.next_a.is_none() || self.next_b.is_none() { return None }

                let a_key = self.next_a.as_ref().unwrap().key();
                let b_key = self.next_b.as_ref().unwrap().key();

                match self
                        .order_operator
                        .judge_cmp( 
                            &a_key, 
                            &b_key 
                        )  {
                    // if next_a < next_b, then delete next_a
                    Ordering::Less => {
                        self.next_a = self.vec_a.next(); 
                    } 
                    // if next_b < next_a, then delete next_b
                    Ordering::Greater => {
                        self.next_b = self.vec_b.next();
                    }
                    // otherwise multiply the coefficients and return
                    Ordering::Equal => {
                        let next_a  =   self.vec_a.next();
                        let next_b = self.vec_b.next();
                        let next_a_old  =   std::mem::replace(&mut self.next_a, next_a).unwrap();
                        let next_b_old  =   std::mem::replace(&mut self.next_b, next_b).unwrap();

                        let mut return_value    =   next_a_old;
                        return_value.set_val(
                            self.ring_operator.multiply(
                                return_value.val(),
                                next_b_old.val(),
                            )
                        );
                        return Some( return_value )
                    }
                }
            }
        }
}

















//  ===========================================================================
//  SINGLE-VECTOR OPERATIONS
//  ===========================================================================








/// Sort the structural nonzero entries by index, using an order operator
/// 
/// This is equivalent to `vec.sort_by(|a,b| order_operator.judge_cmp(& a.key(),& b.key()) )`.
/// 
/// # Example
/// 
/// ```
/// use oat_rust::algebra::vectors::operations::sort_by_order_operator_on_indices;
/// use oat_rust::utilities::order::OrderOperatorAuto; // this implements the default order on Rust objects
/// 
/// let mut vec = vec![ (3, 1.0), (2, 1.0), (1, 1.0) ];
/// let order_operator = OrderOperatorAuto;
/// 
/// sort_by_order_operator_on_indices( &mut vec, order_operator );
/// assert_eq!( vec, vec![ (1, 1.0), (2, 1.0), (3, 1.0) ] );
/// ```
pub fn sort_by_order_operator_on_indices< T, OrderOperator >( vec: &mut Vec< T >, order_operator: OrderOperator ) 
    where
        T:                  KeyValGet,
        OrderOperator:      JudgeOrder< T::Key >,
{
    vec.sort_by(|a,b| order_operator.judge_cmp(& a.key(),& b.key()) )
}


/// Sort the structural nonzero entries by index, using an order operator
/// 
/// This is equivalent to `vec.sort_unstable_by(|a,b| order_operator.judge_cmp(& a.key(),& b.key()) )`.
/// 
/// # Example
/// 
/// ```
/// use oat_rust::algebra::vectors::operations::sort_unstable_by_order_operator_on_indices;
/// use oat_rust::utilities::order::OrderOperatorAuto; // this implements the default order on Rust objects
/// 
/// let mut vec = vec![ (3, 1.0), (2, 1.0), (1, 1.0) ];
/// let order_operator = OrderOperatorAuto;
/// 
/// sort_unstable_by_order_operator_on_indices( &mut vec, order_operator );
/// assert_eq!( vec, vec![ (1, 1.0), (2, 1.0), (3, 1.0) ] );
/// ```
pub fn sort_unstable_by_order_operator_on_indices< T, OrderOperator >( vec: &mut Vec< T >, order_operator: OrderOperator ) 
    where
        T:                  KeyValGet,
        OrderOperator:      JudgeOrder< T::Key >,
{
    vec.sort_unstable_by(|a,b| order_operator.judge_cmp(& a.key(),& b.key()) )
}











/// Convenient methods to transform a vector
///
/// See the lefthand menu for a list of available methods.
/// 
/// These methods are available for structs that
/// implement `Iterator< Item : KeyValGet >`.
/// 
/// Most methods require the entry iterator to be sized; *however* reuse that method
/// this trait is just a convenience.  You can usually perform the desired
/// transformation on you iterator directly, even if you cannot perform it
/// with this trait because the iterator violates the `Sized` criterion.
pub trait VectorOperations< Index, RingElement>

    where   Self:               IntoIterator,
            Self::Item:         KeyValGet< Key = Index, Val = RingElement >,

{

    /// Returns an interator that iterates over the same items as `self`, 
    /// skipping any items with coefficient 0.
    /// 
    /// # Example
    /// 
    /// To drop zero entries, we use a method from the [VectorOperations] trait.
    /// 
    /// ```
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// 
    /// // define inputs
    /// let ring_operator   =   FieldFloat64::new();
    /// let v               =   vec![ (0,  0.), (1,  1.), (2,  0.) ];
    /// let d               =   vec![ (1,  1.) ];
    /// 
    /// // negate
    /// let p               =   v.drop_zeros( ring_operator );
    ///
    /// // verify
    /// assert!( p.eq( d ) );
    /// ```
    fn drop_zeros< RingOperator >
        ( self, ring_operator: RingOperator ) 
        -> 
        DropZeros< Self::IntoIter, RingOperator > 
        
        where   Self:                           Sized,
                RingOperator:                   SemiringOperations< Element = RingElement >,
    {
        DropZeros{ undropped: self.into_iter(), ring_operator, }
    }

    /// Returns an interator that iterates over the same items as `self`, 
    /// with all coefficients scaled by `scalar`.
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// 
    /// // define inputs
    /// let ring_operator   =   FieldFloat64::new();
    /// 
    /// // negate
    /// let p               =   vec![ (0,  0.), (1,  1.) ].scale_by( 2.0, ring_operator );
    ///
    /// // verify
    /// assert!( p.eq( vec![ (0,  0.), (1,  2.) ] ) );
    /// ```
    fn scale_by 
        < RingOperator > 
        ( self, scalar: RingElement, ring_operator: RingOperator )
        -> 
        Scale < Self::IntoIter, RingOperator >
        
        where   Self:                           Sized,
                < Self as IntoIterator>::Item:  KeyValSet + PartialEq,
                RingOperator:                   SemiringOperations< Element = RingElement >,
    {
        Scale{ unscaled: self.into_iter(), scaling_coefficient: scalar, ring_operator, }
    }

    /// Returns an interator that iterates over the same items as `self`, except that 
    /// consecutive entries with equal indices are merged into a single entry whose
    /// coefficient is the sum of the coefficients.  
    /// 
    /// # Example
    /// 
    /// To combine connsecutive entries with equal indices by adding their coefficients, we use a method from the [VectorOperations] trait.
    /// 
    /// ```
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// 
    /// // define inputs
    /// let ring_operator   =   FieldFloat64::new();
    /// let v               =   vec![ (0,  0.), (1,  1.), (1, 1.), (2, 2.), (1, 1.) ].into_iter().peekable();
    /// let g               =   vec![ (0,  0.), (1, 2.), (2, 2.), (1, 1.) ];
    /// 
    /// // negate
    /// let h               =   v.gather( ring_operator );
    ///
    /// // verify
    /// assert!( h.eq( g ) );
    /// ```
    fn gather < RingOperator > ( self, ring_operator: RingOperator )
        -> Gather < Self::IntoIter, RingOperator >

        where   Self:                               Sized,
                Self::IntoIter:                     PeekUnqualified,
                <Self as IntoIterator>::Item:       KeyValGet< Val = RingElement > + KeyValSet < Val = RingOperator::Element > + PartialEq,
                RingOperator:                       SemiringOperations< Element = RingElement >,
                Index:                              PartialEq,              
    {
        Gather{ ungathered: self.into_iter(), ring_operator,  } 
    }

    /// Returns an interator that iterates over the same items as `self`, except that 
    /// (i) consecutive entries with equal indices are merged into a single entry whose
    /// coefficient is the sum of the coefficients, (ii) if this entry has coefficient 0,
    /// then it is omitted.  
    /// 
    /// # Example
    /// 
    /// The `.simplify` method from the [VectorOperations] trait first gathers consecutive entries, then drops zeros.  **Note** The output may have two or more consecutive entries with identical indices, as a result of dropping zeros.  
    /// 
    /// ```
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// 
    /// // define inputs
    /// let ring_operator   =   FieldFloat64::new();
    /// let v               =   vec![ (0,  1.), (1,  1.), (1, -1.), (0,  2.), (3, 3.) ].into_iter().peekable();
    /// let s               =   vec![ (0,  1.), (0,  2.), (3, 3.), ];
    /// 
    /// // negate
    /// let t               =   v.simplify( ring_operator );
    ///
    /// // verify
    /// assert!( t.eq( s ) );
    /// ```
    fn simplify < RingOperator > ( self, ring_operator: RingOperator )
        -> Simplify < Self::IntoIter, RingOperator >

        where   Self:                               Sized,
                Self::IntoIter:                     PeekUnqualified,
                <Self as IntoIterator>::Item:       KeyValSet < Val = RingElement >,
                RingOperator:                       SemiringOperations< Element = RingElement >,
                Index:                              PartialEq,             
    {
        Simplify{ unsimplified: self.into_iter(), ring_operator, } 
    }        

    /// Flip the sign on every entry.
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// 
    /// // define inputs
    /// let ring_operator   =   FieldFloat64::new();
    /// let v               =   vec![ (0,  0.), (1,  1.) ].into_iter().peekable();
    /// let m               =   vec![ (0, -0.), (1, -1.) ];
    /// 
    /// // negate
    /// let n               =   v.negate( ring_operator );
    ///
    /// // verify
    /// assert!( n.eq( m ) );
    /// ```
    fn negate < RingOperator > (self, ring_operator: RingOperator )
        -> Negate< Self::IntoIter, RingOperator >
        where 
            Self:           Sized,
            RingOperator:   RingOperations< Element = RingElement >,
            Self::Item:     KeyValSet + PartialEq,
    {
        Negate::new( self.into_iter(), ring_operator )
    }


    /// Multiply with a matrix, but do not simplify
    /// 
    /// # Arguments
    /// 
    /// - `matrix`: a function that consumes an index and returns a vector
    /// - `ring_operator`: determines the ring operations on the coefficients
    /// - `orger_operator`: specifies the order on indices; 
    /// 
    /// **Each vector returned by `matrix` must have entries sorted in ascending order, acording to `order_operator`**
    /// 
    /// # Example 
    /// 
    /// ```
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::utilities::order::OrderOperatorAuto;
    /// 
    /// // define inputs
    /// let v               =   vec![     (0, 1.), (1, 1.)    ];
    /// let data            =   vec![                         
    ///                             vec![ (0, 1.), (1, 1.)    ],
    ///                             vec![          (1, 1.)    ],
    ///                             vec![                     ],  
    ///                         ];
    /// let matrix          =   |i| { let row: &Vec<_> = &data[i]; row.clone() };
    /// let ring_operator   =   FieldFloat64::new();
    /// let order_operator  =   OrderOperatorAuto::new();
    /// 
    /// // multiply
    /// let u               =   v.multiply_with_matrix_fnmut_unsimplified( matrix, ring_operator, order_operator );
    ///
    /// // verify
    /// assert!(    u.eq( vec![ (0, 1.), (1, 1.), (1, 1.) ] )     );
    /// ```
    fn multiply_with_matrix_fnmut_unsimplified < RingOperator, OrderOperator, Matrix, MatrixRowOrColumn, MatrixIndex, > (
            self, 
            mut matrix: Matrix, 
            ring_operator: RingOperator,
            order_operator: OrderOperator,
        )
        ->  IteratorsMergedInSortedOrder<
                    Scale< MatrixRowOrColumn::IntoIter, RingOperator >,
                    OrderOperator,
                >
        where 
            Self:                                           Sized,    
            Self::Item:                                     KeyValGet < Val = RingElement >,                 
            Matrix:                                         FnMut( Index ) -> MatrixRowOrColumn,
            MatrixRowOrColumn:                              IntoIterator,
            MatrixRowOrColumn::Item:                        KeyValSet< Key = MatrixIndex, Val = RingElement >,
            OrderOperator:                                  JudgePartialOrder< < MatrixRowOrColumn as IntoIterator >::Item >,
            MatrixIndex:                                    PartialEq,
            RingOperator:                                   Clone + SemiringOperations< Element = RingElement >,
    {
        self.into_iter()
                    .map(   |x| 
                            ( x.val(), matrix(x.key()).into_iter() ) 
                        )
                    .linearly_combine_scalar_vector_pairs_without_symplifying( ring_operator, order_operator )
    }  



    /// Multiply with a matrix, but do not simplify;
    /// return `None` if one of the index-coefficient pairs in `self` contains an index which is invalid for the matrix
    /// 
    /// # Arguments
    /// 
    /// - `matrix`: a function that consumes an index and returns an `Option< V >` where `V` is any type that can be treated as a vector
    /// - `ring_operator`: determines the ring operations on the coefficients
    /// - `orger_operator`: specifies the order on matrix entries; 
    /// 
    /// **Each vector returned by `matrix` must have entries sorted in strictly ascending order, acording to index.**
    /// In effect, this means that there needs to be an order on the indices of the matrix, and `order_operator` must
    /// assert that `(i,a) <= (j,b)` iff `i <= j`.
    /// 
    /// # Example 
    /// 
    /// ```
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::utilities::order::OrderOperatorAuto;
    /// 
    /// // define inputs
    /// let v               =   vec![     (0, 1.), (1, 1.)    ];
    /// let data            =   vec![                         
    ///                             vec![ (0, 1.), (1, 1.)    ],
    ///                             vec![          (1, 1.)    ],
    ///                             vec![                     ],  
    ///                         ];
    /// let matrix          =   |i| { let row: &Vec<_> = &data[i]; row.clone() };
    /// let ring_operator   =   FieldFloat64::new();
    /// let order_operator  =   OrderOperatorAuto::new();
    /// 
    /// // multiply
    /// let u               =   v.multiply_with_matrix_fnmut_unsimplified( matrix, ring_operator, order_operator );
    ///
    /// // verify
    /// assert!(    u.eq( vec![ (0, 1.), (1, 1.), (1, 1.) ] )     );
    /// ```
    fn multiply_with_matrix_fnmut_unsimplified_result < RingOperator, OrderOperator, Matrix, MatrixView, MatrixIndex, > (
            self, 
            mut matrix: Matrix, 
            ring_operator: RingOperator,
            order_operator: OrderOperator,
        )
        ->  Result< 
                IteratorsMergedInSortedOrder<
                    Scale< MatrixView::IntoIter, RingOperator, >,
                    OrderOperator,
                >,
                Index, 
            >
        where 
            Self:                                   Sized,    
            Self::Item:                             KeyValGet< Val = RingElement >,                 
            Matrix:                                 FnMut( Index ) -> Result< MatrixView, Index >,            
            MatrixView:                             IntoIterator< 
                                                        Item:   KeyValSet< 
                                                                    Key:    PartialEq,
                                                                    Val=    RingElement 
                                                                >    
                                                    >,            
            OrderOperator:                          JudgePartialOrder< < MatrixView as IntoIterator>::Item >, // JudgePartialOrder< < MatrixView as SparseVector >::Index >,
            MatrixIndex:                            PartialEq,
            RingOperator:                           Clone + SemiringOperations< Element = RingElement >,
    {
        let mut summands = Vec::new();
        for x in self.into_iter() {
            let slice = matrix( x.key() )?;
            summands.push( (x.val(), slice ) );
        }
        Ok( summands
                .linearly_combine_scalar_vector_pairs_without_symplifying(
                    ring_operator, 
                    order_operator,
                ) 
        )
    }      


    /// Multiply with a matrix
    /// 
    /// # Arguments
    /// 
    /// - `matrix`: a function that consumes an index and returns a `V`, where `V` is any type that can be treated as a vector
    /// - `ring_operator`: determines the ring operations on the coefficients
    /// - `orger_operator`: specifies the order on indices; 
    /// 
    /// **Each vector returned by `matrix` must have entries sorted in ascending order, acording to `order_operator`**
    /// 
    /// # Example 
    /// 
    /// ```
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::utilities::order::OrderOperatorAuto;
    /// 
    /// // define inputs
    /// let v               =   vec![     (0, 1.), (1, 1.), (2, 1.)     ];
    /// let data            =   vec![   
    ///                             vec![ (0, 1.), (1, 1.), (2, 1.)     ],
    ///                             vec![          (1, 1.), (2, 1.)     ],
    ///                             vec![                   (2, 1.)     ],  
    ///                         ];
    /// let matrix          =   |i| { let row: &Vec<_> = &data[i]; row.clone() };
    /// let ring_operator   =   FieldFloat64::new();
    /// let order_operator  =   OrderOperatorAuto::new();
    /// 
    /// // multiply
    /// let u               =   v.multiply_matrix_fnmut( matrix, ring_operator, order_operator );
    ///
    /// // verify
    /// assert!( u.eq(  vec![ (0, 1.), (1, 2.), (2, 3.) ]  ) );
    /// ```
    fn multiply_matrix_fnmut < RingOperator, OrderOperator, Matrix, MatrixView, MatrixIndex, > (
            self, 
            matrix: Matrix, 
            ring_operator: RingOperator,
            order_operator: OrderOperator,
        )
        ->  Simplify<
                    IteratorsMergedInSortedOrder<
                            Scale< MatrixView::IntoIter, RingOperator, >,
                            OrderOperator,
                        >,
                    RingOperator,
                >
        where 
            Self:                                   Sized,                   
            Matrix:                                 FnMut( Index ) -> MatrixView,
            MatrixView:                             IntoIterator<
                                                        Item:   KeyValSet<
                                                                    Key = MatrixIndex, 
                                                                    Val = RingOperator::Element 
                                                                >,
                                                    >,        
            OrderOperator:                          JudgePartialOrder< < MatrixView as IntoIterator >::Item >,
            MatrixIndex:                            PartialEq,
            RingOperator:                           Clone + SemiringOperations< Element = RingElement >,
    {
        self.multiply_with_matrix_fnmut_unsimplified(matrix, ring_operator.clone(), order_operator)
            .simplify( ring_operator )
    }  


    /// Multiply with a matrix, returning `None` if an index is invalid
    /// 
    /// # Arguments
    /// 
    /// - `matrix`: a function that consumes an index and returns an `Option< V >` where `V` is any type that can be treated as a vector
    /// - `ring_operator`: determines the ring operations on the coefficients
    /// - `orger_operator`: specifies the order on matrix entries; 
    /// 
    /// **Each vector returned by `matrix` must have entries sorted in strictly ascending order, acording to index. This means
    /// that there must be a total order on indices, and `(i,a) <= (j,b)` if and only if `i <= j`.**
    /// 
    /// # Example 
    /// 
    /// ```
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::utilities::order::OrderOperatorByKey;
    /// 
    /// // define inputs
    /// let v               =   vec![     (0, 1.), (1, 1.), (2, 1.)     ];
    /// let data: Vec<Vec<(usize,f64)>>            =   vec![   
    ///                             vec![ (0, 1.), (1, 1.), (2, 1.)     ],
    ///                             vec![          (1, 1.), (2, 1.)     ],
    ///                             vec![                   (2, 1.)     ],  
    ///                         ];
    /// let matrix          =   |i: usize| { 
    ///                             match i < data.len() {
    ///                                 true    =>  {
    ///                                     let row = data[i].clone();
    ///                                     Ok(row)
    ///                                 },
    ///                                 false   =>  Err( i ),
    ///                             }
    ///                         };
    /// let ring_operator   =   FieldFloat64::new();
    /// let order_operator  =   OrderOperatorByKey::new();
    /// 
    /// // multiply
    /// let u               =   v.multiply_matrix_fnmut_result( matrix, ring_operator, order_operator ).ok().unwrap();
    /// // verify
    /// assert!( u.eq(  vec![ (0, 1.), (1, 2.), (2, 3.) ]  ) );  
    /// ```
    fn multiply_matrix_fnmut_result < RingOperator, OrderOperator, Matrix, MatrixView, MatrixIndex, > (
        self, 
        matrix: Matrix, 
        ring_operator: RingOperator,
        order_operator: OrderOperator,
    )
    ->  Result<
            Simplify<
                IteratorsMergedInSortedOrder<
                    Scale< MatrixView::IntoIter, RingOperator, >,
                    OrderOperator,
                >,
                RingOperator,
            >,
            Index,
        >
    where 
        Self:                                   Sized,                     
        Matrix:                                 FnMut( Index ) -> Result< MatrixView, Index >,
        MatrixView:                             IntoIterator<
                                                    Item:   KeyValSet< 
                                                                Key = MatrixIndex, 
                                                                Val = RingOperator::Element 
                                                            >,    
                                                >,
        OrderOperator:                          JudgePartialOrder< MatrixView::Item >,
        MatrixIndex:                            PartialEq,
        RingOperator:                           Clone + SemiringOperations< Element = RingElement >,
    {
        let unsimplified    =   self.multiply_with_matrix_fnmut_unsimplified_result::<
                                                                        RingOperator, 
                                                                        OrderOperator, 
                                                                        Matrix, 
                                                                        MatrixView,  
                                                                        MatrixIndex
                                                                    >
                                    (
                                        matrix, 
                                        ring_operator.clone(), 
                                        order_operator
                                    )?;
        return Ok( unsimplified.simplify( ring_operator ) )
    }      

    /// Multiply self as a row vector with a matrix
    /// 
    /// # Arguments
    /// 
    /// - `matrix`: a struct that implements [MatrixOracle]
    /// 
    /// # Returns
    /// 
    /// A [simplified](Simplify) vector representing a linear combination of ascending rows, where the combination is specified by the vector.
    /// 
    /// 
    /// # Example 
    /// 
    /// **Note** In particular, entries of the product appear in strictly ascending order.
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket};
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::utilities::order::OrderOperatorAuto;
    /// 
    /// // define inputs
    /// let v               =   vec![     (0, 1.), (1, 1.), (2, 1.)     ];
    /// let data            =   vec![   
    ///                             vec![ (0, 1.), (1, 1.), (2, 1.)     ],
    ///                             vec![          (1, 1.), (2, 1.)     ],
    ///                             vec![                   (2, 1.)     ],  
    ///                         ];
    /// let matrix          =   VecOfVec::new( data ).ok().unwrap();
    /// 
    /// // multiply
    /// let u               =   v.multiply_self_as_a_row_vector_with_matrix_custom( 
    ///                             & matrix,
    ///                             FieldFloat64::new(),
    ///                             OrderOperatorAuto,
    ///                         );
    ///
    /// // verify
    /// assert!( u.eq(  vec![ (0, 1.), (1, 2.), (2, 3.) ]  ) );
    /// ```
    fn multiply_self_as_a_row_vector_with_matrix_custom < Matrix, RingOperator, OrderOperatorForRowEntries, > (
        self, 
        matrix: Matrix,
        ring_operator:  RingOperator,
        order_operator_for_row_entries: OrderOperatorForRowEntries,
    )
    ->  Simplify<
                IteratorsMergedInSortedOrder<
                        Scale< Matrix::Row, RingOperator >,
                        OrderOperatorForRowEntries,
                    >,
                RingOperator,
            >
    where 
        Self:                                   Sized,                 
        Matrix:                                 MatrixOracle< RowIndex=Index >,
        Matrix::RowEntry:                       KeyValSet< Key = Matrix::ColumnIndex, Val = RingOperator::Element >,
        OrderOperatorForRowEntries:                 Clone + JudgePartialOrder< Matrix::RowEntry >,
        Matrix::ColumnIndex:                       PartialEq,
        RingElement:                            Clone,
        RingOperator:                           Clone + SemiringOperations< Element = RingElement >,        
    {
        let matrix = |i| { matrix.row(&i) };
        self.multiply_matrix_fnmut( matrix, ring_operator, order_operator_for_row_entries )
    }
    

    /// Multiply self as a column vector with a matrix, and return entries in reverse order
    /// 
    /// # Arguments
    /// 
    /// - `matrix`: a struct that implements [MatrixOracle]
    /// 
    /// # Returns
    /// 
    /// A [simplified](Simplify) vector representing a linear combination of descending columns, where the combination is specified by the vector.
    /// 
    /// 
    /// # Example 
    /// 
    /// **Note** In particular, entries of the product appear in strictly descending order.
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket};
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::utilities::order::OrderOperatorAuto;
    /// 
    /// // define inputs
    /// let v               =   vec![     (0, 1.), (1, 1.), (2, 1.)     ];
    /// let data            =   vec![   
    ///                             vec![ (0, 1.), (1, 1.), (2, 1.)     ],
    ///                             vec![          (1, 1.), (2, 1.)     ],
    ///                             vec![                   (2, 1.)     ],  
    ///                         ];
    /// let matrix          =   VecOfVec::new( data ).ok().unwrap();
    /// 
    /// // multiply
    /// let u               =   v.multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order_custom( 
    ///                             & matrix,
    ///                             FieldFloat64::new(),
    ///                             OrderOperatorAuto,
    ///                         );
    /// let u: Vec<_> = u.collect();
    /// println!("{:?}", &u );
    /// let u = u.into_iter();
    ///
    /// // verify
    /// assert!( u.eq(  vec![ (2, 1.), (1, 2.), (0, 3.) ]  ) );
    /// ```
    fn multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order_custom < Matrix, RingOperator, OrderOperatorForColumnEntries > (
        self,
        matrix: Matrix,
        ring_operator:  RingOperator,
        order_operator_for_column_entries: OrderOperatorForColumnEntries,
    )
    ->  Simplify<
                IteratorsMergedInSortedOrder<
                        Scale< Matrix::ColumnReverse, RingOperator >,
                        ReverseOrder< OrderOperatorForColumnEntries >,
                    >,
                RingOperator,
            >
    where 
        Self:                                   Sized,                 
        Matrix:                                 MatrixOracle< ColumnIndex = Index >,
        Matrix::ColumnEntry:                    PartialEq + KeyValSet < Key = Matrix::RowIndex, Val = RingOperator::Element >,
        OrderOperatorForColumnEntries:          Clone + JudgePartialOrder< Matrix::ColumnEntry >,
        Matrix::RowIndex:                       PartialEq,
        RingElement:                            Clone,
        RingOperator:                           Clone + SemiringOperations< Element = RingElement >,        
    {
        let matrix_function = |i| { matrix.column_reverse(&i) };
        self.multiply_matrix_fnmut( matrix_function, ring_operator, ReverseOrder::new( order_operator_for_column_entries ) )
    }

    /// Returns `self * M`, where `M` is a matrix, `*` is vector-matrix multiplication, and we regard `self` as a row vector
    /// 
    /// - The vector returned has entries sorted in strictly ascending order, according to the
    ///   order operator of `M`
    /// - The entries of `self` do not have to be sorted in any order
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket};
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// 
    /// // define inputs
    /// let v               =   vec![     (0, 1.), (1, 1.), (2, 1.)     ];
    /// let data            =   vec![   
    ///                             vec![ (0, 1.), (1, 1.), (2, 1.)     ],
    ///                             vec![          (1, 1.), (2, 1.)     ],
    ///                             vec![                   (2, 1.)     ],  
    ///                         ];
    /// let data            =   VecOfVec::new( data ).ok().unwrap();
    /// let matrix          =   MatrixAlgebraPacket::with_default_order_and_f64_coefficients( &data );
    /// 
    /// // multiply
    /// let u               =   v.multiply_self_as_a_row_vector_with_matrix( matrix );
    ///
    /// // verify
    /// assert!( u.eq(  vec![ (0, 1.), (1, 2.), (2, 3.) ]  ) );
    /// ```
    fn multiply_self_as_a_row_vector_with_matrix< Matrix >(
        self,
        matrix: Matrix,
    )
    ->  LinearCombinationOfRows< Matrix >
    where 
        Self:                                   Sized,     
        Matrix:                                 MatrixAlgebra< Coefficient = RingElement, RowIndex = Index >,
        Matrix::RowEntry:                       KeyValSet< Key = Matrix::ColumnIndex, Val = Matrix::Coefficient >,
        Matrix::Coefficient:                    Clone,
        Matrix::RingOperator:                   Clone,
        Matrix::ColumnIndex:                    PartialEq,
    {
        let matrix_function = |i| { matrix.row(&i) };
        let order_operator = matrix.order_operator_for_row_entries();
        self.multiply_matrix_fnmut( matrix_function, matrix.ring_operator(),  order_operator )
    }  

    /// Similar to [multiply_self_as_a_row_vector_with_matrix](VectorOperations::multiply_self_as_a_row_vector_with_matrix), but returns `None` if `self` contains an entry `(k,v)` for any `k` that is an invalid row index of `M`
    /// 
    /// - The vector returned has entries sorted in strictly ascending order, according to the
    ///   order operator of `M`
    /// - The entries of `self` do not have to be sorted in any order
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket};
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::utilities::order::{OrderOperatorAuto, OrderOperatorByKey};
    /// 
    /// // define inputs
    /// let v               =   vec![     (0, 1.), (1, 1.), (2, 1.)     ];
    /// let data            =   vec![   
    ///                             vec![ (0, 1.), (1, 1.), (2, 1.)     ],
    ///                             vec![          (1, 1.), (2, 1.)     ],
    ///                             vec![                   (2, 1.)     ],  
    ///                         ];
    /// let data            =   VecOfVec::new( data ).ok().unwrap();
    /// let matrix          =   MatrixAlgebraPacket{
    ///                             matrix:                             & data,
    ///                             ring_operator:                      FieldFloat64::new(),
    ///                             order_operator_for_row_entries:     OrderOperatorByKey::new(),
    ///                             order_operator_for_row_indices:     OrderOperatorAuto,    
    ///                             order_operator_for_column_entries:  OrderOperatorByKey::new(),
    ///                             order_operator_for_column_indices:  OrderOperatorAuto,
    ///                         };
    /// 
    /// // multiply
    /// let u               =   v.multiply_self_as_a_row_vector_with_matrix_result( matrix );
    ///
    /// // verify
    /// assert!( u.is_ok() );
    /// assert!( u.unwrap().eq( vec![ (0, 1.), (1, 2.), (2, 3.) ] ) );
    /// ```
    fn multiply_self_as_a_row_vector_with_matrix_result< Matrix >(
        self,
        matrix: Matrix,
    )
    ->  Result< 
            LinearCombinationOfRows< Matrix >,
            Index,
        >
    where 
        Self:                                   Sized,     
        Matrix:                                 MatrixAlgebra< Coefficient = RingElement, RowIndex = Index >,
        Matrix::RowEntry:                       KeyValSet< Key = Matrix::ColumnIndex, Val = Matrix::Coefficient >,
        Matrix::Coefficient:                    Clone,        
        Matrix::RingOperator:                   Clone,
        Matrix::ColumnIndex:                    PartialEq,        
    {
        let matrix_function = |i| { matrix.row_result(&i) };
        let order_operator = matrix.order_operator_for_row_entries();
        self.multiply_matrix_fnmut_result( matrix_function, matrix.ring_operator(),  order_operator )
    }      

    /// Returns `self * M`, where `M` is a matrix, `*` is vector-matrix multiplication, and we regard `self` as a row vector;
    /// entries appear in **reverse order**
    /// 
    /// - The vector returned has entries sorted in strictly **descending** order, according to the
    ///   order operator of `M`
    /// - The entries of `self` do not have to be sorted in any order
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket};
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::utilities::order::{OrderOperatorAuto, OrderOperatorByKey};
    /// 
    /// // define inputs
    /// let v               =   vec![     (0, 1.), (1, 1.), (2, 1.)     ];
    /// let data            =   vec![   
    ///                             vec![ (0, 1.), (1, 1.), (2, 1.)     ],
    ///                             vec![          (1, 1.), (2, 1.)     ],
    ///                             vec![                   (2, 1.)     ],  
    ///                         ];
    /// let data            =   VecOfVec::new( data ).ok().unwrap();
    /// let matrix          =   MatrixAlgebraPacket{
    ///                             matrix:                             & data,
    ///                             ring_operator:                      FieldFloat64::new(),
    ///                             order_operator_for_row_entries:     OrderOperatorByKey::new(),
    ///                             order_operator_for_row_indices:     OrderOperatorAuto,    
    ///                             order_operator_for_column_entries:  OrderOperatorByKey::new(),
    ///                             order_operator_for_column_indices:  OrderOperatorAuto,
    ///                         };
    /// 
    /// // multiply
    /// let u               =   v.multiply_self_as_a_row_vector_with_matrix_and_return_entries_in_reverse_order( matrix );
    ///
    /// // verify
    /// assert!( u.eq(  vec![ (2, 3.), (1, 2.), (0, 1.) ]  ) );
    /// ```
    fn multiply_self_as_a_row_vector_with_matrix_and_return_entries_in_reverse_order< Matrix >(
        self,
        matrix: Matrix,
    )
    ->  LinearCombinationOfRowsReverse< Matrix >
    where 
        Self:                                   Sized,     
        Matrix:                                 MatrixAlgebra< Coefficient = RingElement, RowIndex = Index >,
        Matrix::RowEntry:                       KeyValSet< Key = Matrix::ColumnIndex, Val = Matrix::Coefficient >,
        Matrix::Coefficient:                    Clone,       
        Matrix::RingOperator:                   Clone,
        Matrix::ColumnIndex:                    PartialEq,         
    {
        let matrix_function = |i| { matrix.row_reverse(&i) };
        let order_operator = ReverseOrder::new(matrix.order_operator_for_row_entries());
        self.multiply_matrix_fnmut( matrix_function, matrix.ring_operator(),  order_operator )
    }   

    /// Similar to [multiply_self_as_a_row_vector_with_matrix_and_return_entries_in_reverse_order](VectorOperations::multiply_self_as_a_row_vector_with_matrix_and_return_entries_in_reverse_order), but returns `Err(k)` if `self` contains an entry `(k,v)` for any `k` that is an invalid row index of `M`
    /// 
    /// - The vector returned has entries sorted in strictly **descending** order, according to the
    ///   order operator of `M`
    /// - The entries of `self` do not have to be sorted in any order
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket};
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::utilities::order::{OrderOperatorAuto, OrderOperatorByKey};
    /// 
    /// // define inputs
    /// let v               =   vec![     (0, 1.), (1, 1.), (2, 1.)     ];
    /// let data            =   vec![   
    ///                             vec![ (0, 1.), (1, 1.), (2, 1.)     ],
    ///                             vec![          (1, 1.), (2, 1.)     ],
    ///                             vec![                   (2, 1.)     ],  
    ///                         ];
    /// let data            =   VecOfVec::new( data ).ok().unwrap();
    /// let matrix          =   MatrixAlgebraPacket{
    ///                             matrix:                             & data,
    ///                             ring_operator:                      FieldFloat64::new(),
    ///                             order_operator_for_row_entries:     OrderOperatorByKey::new(),
    ///                             order_operator_for_row_indices:     OrderOperatorAuto,    
    ///                             order_operator_for_column_entries:  OrderOperatorByKey::new(),
    ///                             order_operator_for_column_indices:  OrderOperatorAuto,
    ///                         };
    /// 
    /// // multiply
    /// let u               =   v.multiply_self_as_a_row_vector_with_matrix_and_return_entries_in_reverse_order_result( matrix );
    ///
    /// // verify
    /// assert!( u.is_ok() );
    /// assert!( u.unwrap().eq( vec![ (2, 3.), (1, 2.), (0, 1.) ] ) );
    /// ```
    fn multiply_self_as_a_row_vector_with_matrix_and_return_entries_in_reverse_order_result< Matrix >(
        self,
        matrix: Matrix,
    )
    ->  Result< LinearCombinationOfRowsReverse< Matrix >, Index >
    where 
        Self:                                   Sized,     
        Matrix:                                 MatrixAlgebra< Coefficient = RingElement, RowIndex = Index >,
        Matrix::RowEntry:                       KeyValSet< Key = Matrix::ColumnIndex, Val = Matrix::Coefficient >,
        Matrix::Coefficient:                    Clone,       
        Matrix::RingOperator:                   Clone,
        Matrix::ColumnIndex:                    PartialEq,         
    {
        let matrix_function = |i| { matrix.row_reverse_result(&i) };
        let order_operator = ReverseOrder::new(matrix.order_operator_for_row_entries());
        self.multiply_matrix_fnmut_result( matrix_function, matrix.ring_operator(),  order_operator )
    }          


    


    /// Returns `M * self`, where `M` is a matrix, `*` is vector-matrix multiplication, and we regard `self` as a column vector
    /// 
    /// - The vector returned has entries sorted in strictly ascending order, according to the
    ///   order operator of `M`
    /// - The entries of `self` do not have to be sorted in any order
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket};
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// 
    /// // define inputs
    /// let v               =   vec![     (0, 1.), (1, 1.), (2, 1.)     ];
    /// let data            =   vec![   
    ///                             vec![ (0, 1.), (1, 1.), (2, 1.)     ],
    ///                             vec![          (1, 1.), (2, 1.)     ],
    ///                             vec![                   (2, 1.)     ],  
    ///                         ];
    /// let data            =   VecOfVec::new( data ).ok().unwrap();
    /// let matrix          =   MatrixAlgebraPacket::with_default_order_and_f64_coefficients(&data);
    ///        
    /// // multiply
    /// let u               =   v.multiply_self_as_a_column_vector_with_matrix( matrix );
    /// 
    /// // check answer
    /// assert!( u.eq( vec![ (0, 3.), (1, 2.), (2, 1.) ]  ) );
    /// ```
    fn multiply_self_as_a_column_vector_with_matrix< Matrix >(
        self,
        matrix: Matrix,
    )
    ->  LinearCombinationOfColumns< Matrix >
    where 
        Self:                                   Sized,     
        Matrix:                                 MatrixAlgebra< Coefficient=RingElement, ColumnIndex=Index>,
        Matrix::ColumnEntry:                    KeyValSet < Key=Matrix::RowIndex, Val=Matrix::Coefficient >,
        Matrix::Coefficient:                    Clone,        
        Matrix::RingOperator:                   Clone,
        Matrix::RowIndex:                       PartialEq,        
    {
        let matrix_function = |i| { matrix.column(&i) };
        let order_operator = matrix.order_operator_for_column_entries();
        self.multiply_matrix_fnmut( matrix_function, matrix.ring_operator(),  order_operator )
    }  


    /// Similar to [multiply_self_as_a_column_vector_with_matrix](VectorOperations::multiply_self_as_a_column_vector_with_matrix), but returns `Err(k)` if `self` contains an entry `(k,v)` for any `k` that is an invalid column index of `M`
    /// 
    /// - The vector returned has entries sorted in strictly ascending order, according to the
    ///   order operator of `M`
    /// - The entries of `self` do not have to be sorted in any order
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket};
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// 
    /// // define inputs
    /// let v               =   vec![     (0, 1.), (1, 1.), (2, 1.)     ];
    /// let data            =   vec![   
    ///                             vec![ (0, 1.), (1, 1.), (2, 1.)     ],
    ///                             vec![          (1, 1.), (2, 1.)     ],
    ///                             vec![                   (2, 1.)     ],  
    ///                         ];
    /// let data            =   VecOfVec::new( data ).ok().unwrap();
    /// let matrix          =   MatrixAlgebraPacket::with_default_order_and_f64_coefficients(&data);
    /// 
    /// // multiply
    /// let u               =   v.multiply_self_as_a_column_vector_with_matrix_result( matrix );
    ///
    /// // verify
    /// assert!( u.is_ok() );
    /// assert!( u.ok().unwrap().eq( vec![ (0, 3.), (1, 2.), (2, 1.) ] ) );
    /// ```
    fn multiply_self_as_a_column_vector_with_matrix_result< Matrix >(
        self,
        matrix: Matrix,
    )
    ->  Result< LinearCombinationOfColumns< Matrix >, Index >
    where 
        Self:                                   Sized,     
        Matrix:                                 MatrixAlgebra< Coefficient=RingElement, ColumnIndex=Index>,
        Matrix::ColumnEntry:                    KeyValSet < Key=Matrix::RowIndex, Val=Matrix::Coefficient >,
        Matrix::Coefficient:                    Clone,
        Matrix::RingOperator:                   Clone,   
        Matrix::RowIndex:                       PartialEq,                   
    {
        let matrix_function = |i| { matrix.column_result(&i) };
        let order_operator = matrix.order_operator_for_column_entries();
        self.multiply_matrix_fnmut_result( matrix_function, matrix.ring_operator(),  order_operator )
    }  





    /// Similar to [multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order](VectorOperations::multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order), but returns `Err(k)` if `self` contains an entry `(k,v)` for any `k` that is an invalid column index of `M`
    /// 
    /// - The vector returned has entries sorted in strictly **descending** order, according to the
    ///   order operator of `M`
    /// - The entries of `self` do not have to be sorted in any order
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket};
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::utilities::order::{OrderOperatorAuto, OrderOperatorByKey};
    /// 
    /// // define inputs
    /// let v               =   vec![     (0, 1.), (1, 1.), (2, 1.)     ];
    /// let data            =   vec![   
    ///                             vec![ (0, 1.), (1, 1.), (2, 1.)     ],
    ///                             vec![          (1, 1.), (2, 1.)     ],
    ///                             vec![                   (2, 1.)     ],  
    ///                         ];
    /// let data            =   VecOfVec::new( data ).ok().unwrap();
    /// let matrix          =   MatrixAlgebraPacket{
    ///                             matrix:                             & data,
    ///                             ring_operator:                      FieldFloat64::new(),
    ///                             order_operator_for_row_entries:     OrderOperatorByKey::new(),
    ///                             order_operator_for_row_indices:     OrderOperatorAuto,    
    ///                             order_operator_for_column_entries:  OrderOperatorByKey::new(),
    ///                             order_operator_for_column_indices:  OrderOperatorAuto,
    ///                         };
    /// 
    /// // multiply
    /// let u               =   v.multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order_result( matrix );
    ///
    /// // verify
    /// assert!( u.is_ok() );
    /// assert!( u.unwrap().eq( vec![ (2, 1.), (1, 2.), (0, 3.) ] ) );
    /// ```
    fn multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order_result< Matrix >(
        self,
        matrix: Matrix,
    )
    ->  Result< LinearCombinationOfColumnsReverse< Matrix >, Index >
    where 
        Self:                                   Sized,     
        Matrix:                                 MatrixAlgebra< Coefficient=RingElement, ColumnIndex=Index >, 
        Matrix::ColumnEntry:                    KeyValSet < Key=Matrix::RowIndex, Val=Matrix::Coefficient >,
        Matrix::Coefficient:                    Clone,   
        Matrix::RingOperator:                   Clone,
        Matrix::RowIndex:                       PartialEq,             
    {
        let matrix_function = |i| { matrix.column_reverse_result(&i) };
        let order_operator = ReverseOrder::new(matrix.order_operator_for_column_entries());
        self.multiply_matrix_fnmut_result( matrix_function, matrix.ring_operator(),  order_operator )
    }  


    /// Returns `M * self`, where `M` is a matrix, `*` is vector-matrix multiplication, and we regard `self` as a column vector;
    /// entries appear in **reverse order**
    /// 
    /// - The vector returned has entries sorted in strictly **descending** order, according to the
    ///   order operator of `M`
    /// - The entries of `self` do not have to be sorted in any order
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket};
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::utilities::order::{OrderOperatorAuto, OrderOperatorByKey};
    /// 
    /// // define inputs
    /// let v               =   vec![     (0, 1.), (1, 1.), (2, 1.)     ];
    /// let data            =   vec![   
    ///                             vec![ (0, 1.), (1, 1.), (2, 1.)     ],
    ///                             vec![          (1, 1.), (2, 1.)     ],
    ///                             vec![                   (2, 1.)     ],  
    ///                         ];
    /// let data            =   VecOfVec::new( data ).ok().unwrap();
    /// let matrix          =   MatrixAlgebraPacket{
    ///                             matrix:                             & data,
    ///                             ring_operator:                      FieldFloat64::new(),
    ///                             order_operator_for_row_entries:     OrderOperatorByKey::new(),
    ///                             order_operator_for_row_indices:     OrderOperatorAuto,    
    ///                             order_operator_for_column_entries:  OrderOperatorByKey::new(),
    ///                             order_operator_for_column_indices:  OrderOperatorAuto,
    ///                         };
    /// 
    /// // multiply
    /// let u               =   v.multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order( matrix );
    ///
    /// // verify
    /// assert!( u.eq( vec![ (2, 1.), (1, 2.), (0, 3.) ] ) );
    /// ```
    fn multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order< Matrix >(
        self,
        matrix: Matrix,
    )
    ->  LinearCombinationOfColumnsReverse< Matrix >
    where 
        Self:                                   Sized,     
        Matrix:                                 MatrixAlgebra< Coefficient=RingElement, ColumnIndex=Index >,     
        Matrix::ColumnEntry:                    KeyValSet < Key=Matrix::RowIndex, Val=Matrix::Coefficient >,
        Matrix::Coefficient:                    Clone,   
        Matrix::RingOperator:                   Clone,
        Matrix::RowIndex:                       PartialEq,   
    {
        let matrix_function = |i| { matrix.column_reverse(&i) };
        let order_operator = ReverseOrder::new(matrix.order_operator_for_column_entries());
        self.multiply_matrix_fnmut( matrix_function, matrix.ring_operator(),  order_operator )
    }          

    
             


             





    /// Multiply two vectors entry-wise, returning a new vector
    /// 
    /// Concretely, if `u = [u0 .. um]` and `v = [v0 .. vm ]` are vectors, then this returns
    /// a vector of form `[u0 * v0 .. um * vm]`.
    /// 
    /// # Caution
    /// 
    /// - Assumes that the entries of both `u` and `v` are sorted in strictly ascending order of index, with respect to the [order operator](crate::utilities::order).
    /// - The resulting vector will store an explicit zero entry at index `i`
    ///   if `u` and `v` both store explicit entries for index `i`, and at least one of those entries is zero.
    ///   **Therefore the output is only guaranteed to be simplified if the two inputs are simplified**.
    /// 
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::utilities::order::OrderOperatorAuto;
    /// 
    /// // define inputs
    /// let u                   =   vec![     (0, 0.),          (2, 2.)  ];
    /// let v                   =   vec![     (0, 0.), (1, 1.), (2, 2.)  ];
    /// 
    /// let product             =   u.multiply_entrywise(
    ///                                 v,
    ///                                 FieldFloat64::new(),
    ///                                 OrderOperatorAuto,
    ///                             );
    /// let product: Vec<_>     =   product.collect(); // collect the entries into a vector
    /// assert_eq!( product, vec![ (0,0.), (2,4.) ] )
    /// ```
    fn multiply_entrywise
        < RingOperator, OrderOperator, OtherIterable > 
        ( self, other_iterable: OtherIterable, ring_operator: RingOperator, order_operator: OrderOperator )
        -> 
        EntrywiseProduct< Self::IntoIter, OtherIterable::IntoIter, RingOperator, OrderOperator >
        
        where   Self:                           Sized,
                OtherIterable:                  IntoIterator,
                OtherIterable::Item:            KeyValGet < Key = Index, Val = RingElement >,
                RingOperator:                   SemiringOperations< Element = RingElement >,
                OrderOperator:                  JudgeOrder< Index >,
    {
        let mut vec_a   =   self.into_iter();
        let mut vec_b   =   other_iterable.into_iter();
        let next_a = vec_a.next();
        let next_b = vec_b.next();
        EntrywiseProduct{
            vec_a,
            vec_b,
            next_a,
            next_b,
            ring_operator,
            order_operator,
        }
    } 
  


    /// Computes the dot product with another vector
    /// 
    /// Concretely, if `u = [u0 .. um]` and `v = [v0 .. vm ]` are vectors, then this returns
    ///  `(u0 * v0)  + ..  + (um * vm)`.
    /// 
    /// Returns zero if the vector stores no explict entries.
    /// 
    /// # Caution 
    /// 
    /// **Assumes that the entries of both `u` and `v` are sorted in strictly ascending order of index, with respect to the [order operator](crate::utilities::order).**
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::utilities::order::OrderOperatorAuto;
    /// 
    /// // define inputs
    /// let v                   =   vec![     (0, 0.), (1, 1.), (2, 2.)  ];
    /// let u                   =   vec![     (0, 0.),          (2, 2.)  ];
    /// 
    /// let dot                 =   v.dot(
    ///                                 u,
    ///                                 FieldFloat64::new(),
    ///                                 OrderOperatorAuto,
    ///                             );
    /// assert_eq!( dot, 4.0);
    /// ```
    /// 
    /// # See also
    /// 
    /// The [VectorOperations::dot] method is relatively efficient, because it uses an order operator to reduce the time and
    /// memory used in the calculation. If your vectors are not sorted, or if you don't have an order operator, then you can
    /// use the [VectorOperations::dot_slow] method instead.
    fn dot
        < RingOperator, OrderOperator, OtherIterable > 
        ( self, other_iterable: OtherIterable, ring_operator: RingOperator, order_operator: OrderOperator )
        -> 
        RingElement
        
        where   Self:                           Sized,
                <Self as IntoIterator>::Item:   KeyValSet < Val = RingElement >,
                OtherIterable:                  IntoIterator,
                OtherIterable::Item:            KeyValGet< Key = Index, Val = RingElement >,
                RingOperator:                   Clone + SemiringOperations< Element = RingElement >,
                OrderOperator:                  JudgeOrder< Index >,
    {
        self.multiply_entrywise(
                other_iterable,
                ring_operator.clone(),
                order_operator,
            )
            .sum_coefficients(ring_operator)
        // let vec_a = self.into_iter();
        // let vec_b = other_iterator.into_iter();
        // let next_a = vec_a.next();
        // let next_b = vec_b.next();
        // EntrywiseProduct{
        //     vec_a,
        //     vec_b,
        //     next_a,
        //     next_b,
        //     ring_operator:  ring_operator.clone(),
        //     order_operator,
        //     phantom_coefficient: PhantomData,
        //     phantom_index: PhantomData,
        // }.sum_coefficients( ring_operator )
    }    


    /// Computes the dot product with another vector
    /// 
    /// Concretely, if `u = [u0 .. um]` and `v = [v0 .. vm ]` are vectors, then this returns
    ///  `(u0 * v0)  + ..  + (um * vm)`.
    /// 
    /// Returns zero if the vector `self` stores no explict entries.
    /// 
    /// This method works by emptying all the key-value pairs from `u` into a `HashMap< Key, Val >`; it then
    /// iterates over `v`, using the hashmap to look up corresponding entries of `u`.
    /// 
    /// # Caution 
    /// 
    /// **Assumes that the entries of both `u` and `v` are sorted in strictly ascending order of index, with respect to the [order operator](crate::utilities::order).**
    /// 
    /// # See also
    /// 
    /// The [VectorOperations::dot] method is at least as efficient as [VectorOperations::dot_slow], and sometimes much more
    /// efficient in terms of time and memory. This is because [VectorOperations::dot] never has to empty its entries into
    /// a hash map.
    /// So, if your vectors are sorted and you have an [order operator](crate::utilities::order), 
    /// the [VectorOperations::dot] method will serve you better.
    /// 
    /// # Examples
    /// 
    /// ```
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// 
    /// // define inputs
    /// let v                   =   vec![     (0, 0.), (1, 1.), (2, 2.)  ];
    /// let u                   =   vec![     (0, 0.),          (2, 2.)  ];
    /// 
    /// let dot                 =   v.dot_slow(
    ///                                 u,
    ///                                 FieldFloat64::new(),
    ///                             );
    /// assert_eq!( dot, 4.0);
    /// ```
    fn dot_slow
        < RingOperator, OtherIterable > 
        ( self, other_iterator: OtherIterable, ring_operator: RingOperator, )
        -> 
        RingElement
        
        where   Self:                           Sized,
                <Self as IntoIterator>::Item:   KeyValSet < Val = RingElement >,
                OtherIterable:                  IntoIterator,
                OtherIterable::Item:            KeyValGet< Key = Index, Val = RingElement >,
                RingOperator:                   Clone + SemiringOperations< Element= RingElement >,
                Index:                          Hash + std::cmp::Eq,
    {
        let mut u: HashMap< Index, RingElement > =  self.into_iter()
                                                    .map(|x| (x.key(), x.val()) )
                                                    .collect();
        let mut dot             =   RingOperator::zero();
        for entry_a in other_iterator.into_iter() {
            if let Some( entry_b_coefficient )           =   u.remove( &entry_a.key() ) {
                let product     =   ring_operator.multiply( entry_a.val(), entry_b_coefficient );
                dot                         = < RingOperator as SemiringOperations >::add( & ring_operator, dot, product );
            }
        }
        return dot
    }         


    /// Add a vector, without simplifying
    fn add_unsimplified < OrderOperator, Other > (self, other: Other, order_operator: OrderOperator )
        -> MergeTwoIteratorsByOrderOperator< Self::IntoIter, Other::IntoIter, OrderOperator,>
        where 
            Self:                               Sized,
            < Self as IntoIterator>::IntoIter:  PeekUnqualified,                    
            Other:                              IntoIterator< Item = < Self as IntoIterator >::Item >,
            Other::IntoIter:                    PeekUnqualified,
            OrderOperator:                      JudgePartialOrder< < Self as IntoIterator >::Item >,
            Index:                              PartialEq,
    {
        MergeTwoIteratorsByOrderOperator::new( self.into_iter(), other.into_iter(), order_operator )
    }   

    /// Add another vector, and simplify
    /// 
    /// See also the [VectorOperations::sum] method in the [VectorOperations] trait, for adding `k ≥ 3` vectors.
    /// 
    /// # Examples
    /// 
    /// They require the user to specify an order operator (which says which entries 
    /// are supposed to precede others) and a ring operator (which specifies addition, subtraction, etc. of coefficients).
    /// Note that they require iterators to implement the `PeekUnqualified` trait;
    /// if your iterator doesn't satisfy this condition, you can transform it into
    /// one that does by calling `iter.peekable()`.
    /// 
    /// ```
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::utilities::order::OrderOperatorAuto;
    /// 
    /// // define inputs
    /// let order_operator  =   OrderOperatorAuto::new();
    /// let ring_operator   =   FieldFloat64::new();
    /// let a               =   vec![ (1,1.), (2,2.) ].into_iter().peekable();
    /// let b               =   vec![ (2,2.), (3,3.) ].into_iter().peekable();
    /// 
    /// // add
    /// let sum             =   a.add( b, ring_operator, order_operator );
    /// 
    /// // verify
    /// assert!( sum.eq( vec![ (1,1.), (2,4.), (3,3.) ] ) ) 
    /// ```
    /// 
    /// 
    /// #### Alternate addition strategies
    /// 
    /// It's also possible to implement addition using lower level operations.  The standard procedure is currently as follows:
    /// 
    /// - First combine the entries of two or more sparse vector iterators into a single iterator, ensuring that entries with equal indices group together (for example, `[ (0,1), (1,1), (1,2) ]` and `[ (1,2), (1,1), (0,1) ]` are ok but `[ (1,1), (0,1), (1,2) ]` is not).  [There are many tools to make merging sorted iterators easier.](crate::utilities::iterators::merge).
    /// - Then simplify the merged iterator by adding terms
    /// with equal entries, and dropping zeros.
    /// 
    /// Add two sorted vectors using the merge method from the popular Itertools library:
    /// 
    /// ```
    /// //  import definitions
    /// use itertools::Itertools;
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// 
    /// // add
    /// let sum         =   vec![ (1,1.), (2,2.) ]
    ///                         .into_iter() // convert the vector to an iterator
    ///                         .merge( vec![ (2,2.), (3,3.) ] ) // merge with another vector
    ///                         .peekable() // make the first entry of the merged vector peekable
    ///                         .simplify( FieldFloat64::new() ); // combine entries with equal indices, and drop zeros
    /// 
    /// // verify
    /// assert!( sum.eq( vec![ (1,1.), (2,4.), (3,3.) ] ) )
    /// ```
    /// 
    /// 
    /// 
    fn add < RingOperator, OrderOperator, Other > (self, other: Other, ring_operator: RingOperator, order_operator: OrderOperator )
        ->  Simplify< 
                    MergeTwoIteratorsByOrderOperator< 
                            Self::IntoIter, 
                            Other::IntoIter,
                            OrderOperator,
                        >, 
                    RingOperator, 
                >
        where 
            Self:                               Sized,
            < Self as IntoIterator>::IntoIter:  PeekUnqualified,
            < Self as IntoIterator >::Item:     PartialEq + KeyValSet < Val = RingElement >,            
            Other:                              IntoIterator< Item = < Self as IntoIterator >::Item >,
            Other::IntoIter:                    PeekUnqualified,
            OrderOperator:                      JudgePartialOrder< < Self as IntoIterator >::Item >,
            RingOperator:                       SemiringOperations< Element = RingElement >,
            Index:                              PartialEq,
    {
        MergeTwoIteratorsByOrderOperator::new( self.into_iter(), other.into_iter(), order_operator ).simplify(ring_operator)
    }    

    /// Subtract a vector, without simplifying
    fn subtract_unsimplified < RingOperator, OrderOperator, Other > (self, other: Other, ring_operator: RingOperator, order_operator: OrderOperator )
        -> MergeTwoIteratorsByOrderOperator< Self::IntoIter, Peekable< Negate<Other::IntoIter,RingOperator> >,OrderOperator,>
        where 
            Self:                               Sized, // peeking required for MergeTwoIteratorsByOrderOperator
            < Self as IntoIterator>::IntoIter:  PeekUnqualified, // peeking required for MergeTwoIteratorsByOrderOperator
            < Self as IntoIterator >::Item:     PartialEq + KeyValSet < Val = RingElement >,              
            Other:                              IntoIterator< Item = < Self as IntoIterator >::Item >,
            Other::IntoIter:                    PeekUnqualified,            
            OrderOperator:                      JudgePartialOrder< < Self as IntoIterator >::Item >,
            RingOperator:                       Clone + RingOperations< Element = RingElement >,            
            Index:                              PartialEq,
    {
        let pos = self.into_iter();
        let neg = Negate::new( other.into_iter(), ring_operator ).peekable();
        
        MergeTwoIteratorsByOrderOperator::new( pos, neg, order_operator )
    }   

    /// Subtract a vector, and simplify
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::utilities::order::OrderOperatorAuto;
    /// 
    /// // define inputs
    /// let order_operator  =   OrderOperatorAuto::new();
    /// let ring_operator   =   FieldFloat64::new();
    /// let a               =   vec![ (1,1.), (2,2.) ].into_iter().peekable();
    /// let b               =   vec![ (2,2.), (3,3.) ].into_iter().peekable();
    /// 
    /// // add
    /// let sum             =   a.subtract( b, ring_operator, order_operator );
    /// 
    /// // verify
    /// assert!( sum.eq( vec![ (1,1.), (3,-3.) ] ) )  
    /// ```
    /// 
    /// #### Low-level implemntation
    /// 
    /// It's also possible to implement subtraction using lower-level operations:
    /// 
    /// ```
    /// //  import definitions
    /// use itertools::Itertools;
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// 
    /// //  define the terms
    /// let a   =   vec![ (1,1.), (2,2.) ];
    /// let b   =   vec![ (2,2.), (3,3.) ];
    /// 
    /// //  subtract b from a
    /// let sum         =   a
    ///                         .into_iter() // convert the vector to an iterator
    ///                         .merge( // merge with negative b
    ///                             b.into_iter() // convert b to an iterator
    ///                                 .negate( FieldFloat64::new() ) // convert b to negative b
    ///                         ) 
    ///                         .peekable() // make the first entry of the merged vector peekable (a technical requirement)
    ///                         .simplify( FieldFloat64::new() ); // combine entries with equal indices, and drop zeros
    /// 
    /// //  verify
    /// assert!( sum.eq( vec![ (1,1.), (3,-3.) ] ) )   
    /// ```
    fn subtract < RingOperator, OrderOperator, Other > (self, other: Other, ring_operator: RingOperator, order_operator: OrderOperator )
        ->  Simplify< 
                    MergeTwoIteratorsByOrderOperator< 
                            Self::IntoIter, 
                            Peekable< Negate< Other::IntoIter,RingOperator > >,
                            OrderOperator,
                        >, 
                    RingOperator, 
                >
        where 
            Self:                               Sized,
            < Self as IntoIterator>::IntoIter:  PeekUnqualified,
            < Self as IntoIterator >::Item:     PartialEq + KeyValSet < Val = RingElement >,              
            Other:                              IntoIterator< Item = < Self as IntoIterator >::Item >,
            Other::IntoIter:                    PeekUnqualified,
            OrderOperator:                      JudgePartialOrder< < Self as IntoIterator >::Item >,
            RingOperator:                       Clone + RingOperations< Element = RingElement >,            
            Index:                              PartialEq,
    {
        self.into_iter().subtract_unsimplified(other, ring_operator.clone(), order_operator).simplify(ring_operator)
    }    


    /// Remove entries whose indices fall **outside** of a given collection.
    ///
    ///  # Example
    /// 
    /// ```
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// 
    /// // define vector
    /// let v                       =   vec![ (1,1.), (3,3.), (5,5.) ];
    /// 
    /// // define indices to exclude
    /// let indices_to_include      =   vec![ 1, 2 ];
    /// 
    /// // exclude entries outside the collection
    /// let v               =   v.include_collection( indices_to_include );
    /// 
    /// // verify
    /// assert!( v.eq( vec![ (1,1.) ] ) )
    /// ```
    fn include_collection
        < CollectionToInclude: MapHasKeyOrSequenceHasElement<Index> >
        ( self, collection_to_include: CollectionToInclude ) 
        -> 
        OnlyIndicesInsideCollection< Self::IntoIter, CollectionToInclude, > 
            where
                CollectionToInclude:    MapHasKeyOrSequenceHasElement<Index>,
                Self:                   Sized,
    {
        OnlyIndicesInsideCollection::new( self.into_iter(), collection_to_include ) 
    }

    /// Remove entries whose indices fall **inside** of a given collection.
    ///
    ///  # Example
    /// 
    /// ```
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// 
    /// // define vector
    /// let v                       =   vec![ (1,1.), (3,3.), (5,5.) ];
    /// 
    /// // define indices to exclude
    /// let indices_to_include      =   vec![ 1, 2 ];
    /// 
    /// // exclude entries outside the collection
    /// let v               =   v.exclude_collection( indices_to_include );
    /// 
    /// // verify
    /// assert!( v.eq( vec![ (3,3.), (5,5.) ] ) )
    /// ```
    fn exclude_collection
        < CollectionToExclude: MapHasKeyOrSequenceHasElement<Index> >
        ( self, collection_to_exclude: CollectionToExclude ) 
        -> 
        OnlyIndicesOutsideCollection< Self::IntoIter, CollectionToExclude, > 
            where
                Self:           Sized,
    {
        OnlyIndicesOutsideCollection::new( self.into_iter(), collection_to_exclude ) 
    }





    /// Sum the coefficients in a vector
    /// 
    /// This is equivalent to `self.hit_merge_by_predicate( order_operator ).simplify( ring_operator )`, and you can use that method in order to avoid the `Sized` requirement.
    /// 
    ///  # Example
    /// 
    /// ```
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::utilities::order::OrderOperatorAuto;
    /// 
    /// // define inputs
    /// let ring_operator   =   FieldFloat64::new();
    /// let v               =   vec![ (1,1.), (2,2.) ].into_iter().peekable();
    /// 
    /// // add
    /// let sum             =   v.sum_coefficients( ring_operator );
    /// 
    /// // verify
    /// assert_eq!( sum, 3.0 )  
    /// ```
    fn sum_coefficients< RingOperator >( self, ring_operator: RingOperator ) -> RingElement

        where 
            Self:                               Sized,
            // Self::Item:                         IntoIterator,        
            < Self as IntoIterator >::Item:     KeyValGet< Val = RingElement >,
            RingOperator:                       SemiringOperations< Element = RingElement >,
    {
        let mut sum     =   RingOperator::zero();
        for entry in self { sum = ring_operator.add( sum, entry.val() ) }
        return sum
    }



        
}

// We implement this trait automatically on all iterators.
impl    < T, Index, RingElement >

        VectorOperations < Index, RingElement >
        
        for T

        where 
            T:          IntoIterator,
            T::Item:    KeyValGet< Key = Index, Val = RingElement >,
  
{} // everything implemented automatically




















//  ===========================================================================
//  MULTI-VECTOR OPERATIONS
//  ===========================================================================






/// Convenient methods to transform a vector
///
/// See the lefthand menu for a list of available methods.
/// 
/// **This trait is auto-implemented for every struct; however, the methods it provides are
/// only available for iterables that run over [sparse vectors](crate::algebra::vectors), ordered pairs of form (scalar, sparse_vector), or similar.**
/// 
/// Most methods require the entry iterator to be sized; *however* reuse that method
/// this trait is just a convenience.  You can usually perform the desired
/// transformation on you iterator directly, even if you cannot perform it
/// with this trait because the iterator violates the `Sized` criterion.
pub trait MultiVectorOperations
{

    /// Returns a linear combination, without simplifying
    /// 
    /// This method is called on an iterable that runs over `(scalar, vector)` pairs.
    /// It scales each vector with the corresponding scalar, then merges the scaled
    /// vectors into a single vector. Entries with the same index are *not* combined, so
    /// multiple entries may have the same index.
    /// 
    /// **Caution** In order to ensure that the resulting vector has entries in sorted order
    /// each of the input vectors must have entries in sorted order, with 
    /// respect to the [order operator](crate::utilities::order)
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::vectors::operations::MultiVectorOperations;
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::utilities::order::OrderOperatorAuto;;
    /// 
    /// // define inputs
    /// let ring_operator   =   FieldFloat64::new();
    /// let order_operator  =   OrderOperatorAuto::new();
    /// let u               =   vec![ (0usize,  1.),                ].into_iter().peekable();
    /// let v               =   vec![ (0usize,  1.), (1usize,  1.), ].into_iter().peekable();
    /// 
    /// // combine
    /// let combination     =   [ (0.0, u), (7.0, v) ]
    ///                             .linearly_combine_scalar_vector_pairs_without_symplifying( 
    ///                                 ring_operator, 
    ///                                 order_operator 
    ///                             );
    ///
    /// // verify
    /// assert!( combination.eq(  vec![ (0, 0.), (0, 7.), (1, 7.) ] ) );
    /// ```
    fn linearly_combine_scalar_vector_pairs_without_symplifying< RingOperator, OrderOperator, Vector, Index, RingElement >( 
            self, 
            ring_operator: RingOperator, 
            order_operator: OrderOperator 
        )
            -> 
            IteratorsMergedInSortedOrder<
                    Scale< Vector::IntoIter, RingOperator >,
                    OrderOperator,
                >
        where
            Self:                   Sized,        
            Self:                   IntoIterator< Item = (RingOperator::Element, Vector ) >,

            // ------------------------------------------------------
            Vector:                 IntoIterator<
                                        IntoIter:   Iterator<
                                                        Item:   KeyValSet< 
                                                                    Key = Index, 
                                                                    Val = RingElement 
                                                                >,
                                                    >
                                    >,
            Index:                  PartialEq,
            // ------------------------------------------------------            

            OrderOperator:          JudgePartialOrder< Vector::Item >,
            RingOperator:           Clone + SemiringOperations< Element = RingElement >,
    {
        self.into_iter().map(
                    | (scale, unscaled) |
                    Scale{ unscaled: unscaled.into_iter(), scaling_coefficient: scale, ring_operator: ring_operator.clone(), }
                )
            .hit_merge_by_predicate(order_operator)
    }

    /// Returns a simplified linear combination
    /// 
    /// This method is called on an iterable that runs over `(scalar, vector)` pairs.
    /// It scales each vector with the corresponding scalar, then merges the scaled
    /// vectors into a single vector. To wit, this method returns the linear combination
    /// of the given vectors with their respective scalars.
    /// 
    /// **Caution** Each of the input vectors must have entries in sorted order, with 
    /// respect to the [order operator](crate::utilities::order).
    /// 
    /// # Example
    /// 
    /// ```
    /// use oat_rust::algebra::vectors::operations::MultiVectorOperations;
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::utilities::order::OrderOperatorAuto;;
    /// 
    /// // define inputs
    /// let ring_operator   =   FieldFloat64::new();
    /// let order_operator  =   OrderOperatorAuto::new();
    /// let u               =   vec![ (0,  0.), (1,  1.), (2,  2.) ].into_iter().peekable();
    /// let v               =   vec![ (1,  1.), (2,  2.), (3,  3.) ].into_iter().peekable();
    /// let w               =   vec![ (1, -1.), (2, -2.), (3, -6.) ];
    /// 
    /// // combine
    /// let comb            =   [ (1.0, u), (-2.0, v) ].linearly_combine_scalar_vector_pairs( ring_operator, order_operator );
    ///
    /// // verify
    /// assert!( comb.eq( w ) );
    /// ```
    fn linearly_combine_scalar_vector_pairs< RingOperator, OrderOperator, Vector, RingElement >( 
        self, 
        ring_operator:      RingOperator, 
        order_operator:     OrderOperator 
    )
        -> 
        Simplify< 
            IteratorsMergedInSortedOrder<
                Scale< Vector::IntoIter, RingOperator, >,
                OrderOperator,
            >,
            RingOperator,        
        >

        where
            Self:                                   Sized,        
            Self:                                   IntoIterator< Item = (RingElement, Vector ) >,
            Vector:                                 IntoIterator,
            Vector::Item:                           PartialEq + KeyValSet < Val = RingElement >,
            < Vector::Item as KeyValGet >::Key:     PartialEq,
            OrderOperator:                          JudgePartialOrder< Vector::Item >,
            RingOperator:                           Clone + SemiringOperations< Element = RingElement >,
    {
        self.into_iter().map(
                    | (scale, unscaled) |
                    Scale{ unscaled: unscaled.into_iter(), scaling_coefficient: scale, ring_operator: ring_operator.clone(), }
                )
            .hit_merge_by_predicate(order_operator)
            .simplify( ring_operator )
    }    







    /// Sum two or more vectors without simplifying
    /// 
    /// The resulting iterator may have **multiple entries with the same index**.
    /// 
    /// This is equivalent to `self.hit_merge_by_predicate( order_operator )` (and you can use that method in order to avoid the `Sized` requirement)
    /// 
    /// # See also
    /// 
    /// The method [VectorOperations::sum_unsimplified] calls on a iterator of vectors and returns a vector. Contrast this with
    /// 
    /// - [VectorOperations::add_unsimplified], which is a method called on a vector, and returns a vector
    /// - [VectorOperations::sum_coefficients], which is a method called on a vector, and returns a scalar
    fn sum_vectors_unsimplified< OrderOperator >( self, order_operator: OrderOperator )
            -> IteratorsMergedInSortedOrder< <Self::Item as IntoIterator >::IntoIter, OrderOperator >

        where 
            Self:               Sized,
            Self:               IntoIterator,
            Self::Item:         IntoIterator,  
            OrderOperator:      JudgePartialOrder<
                                        < <Self as IntoIterator >::Item as IntoIterator >::Item 
                                    >,
    {
        self.hit_merge_by_predicate( order_operator )
    }    



    /// Sum a collection of vectors, and simplify
    /// 
    /// This is equivalent to `self.hit_merge_by_predicate( order_operator ).simplify( ring_operator )`; you can use that method instead if you need to avoid the `Sized` requirement.
    /// 
    ///  # Example
    /// 
    /// ```
    /// use oat_rust::algebra::vectors::operations::MultiVectorOperations;
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::utilities::order::OrderOperatorAuto;
    /// 
    /// // define inputs
    /// let order_operator  =   OrderOperatorAuto::new();
    /// let ring_operator   =   FieldFloat64::new();
    /// let a               =   vec![ (1,1.), (2,2.) ].into_iter().peekable();
    /// let b               =   vec![ (2,2.), (3,3.) ].into_iter().peekable();
    /// let c               =   vec![ (3,3.), (4,4.) ].into_iter().peekable();
    /// 
    /// // add
    /// let sum             =   [a,b,c].sum_vectors( ring_operator, order_operator );
    /// 
    /// // verify
    /// assert!( sum.eq( vec![ (1,1.), (2,4.), (3,6.), (4,4.) ] ) )  
    /// ```
    /// 
    /// # Alternate approaches
    /// 
    /// There are many alternate strategies to add vectors using low-level operations:
    /// 
    /// Add three sorted vectors using the kmerge method from the popular Itertools library. 
    /// 
    /// ```
    /// //  import definitions
    /// use itertools::Itertools;
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// 
    /// // add
    /// let vectors     =   vec![ vec![ (1,1.), (2,2.) ], vec![ (2,2.), (3,3.) ], vec![ (3,3.), (4,4.) ] ];
    /// let sum         =   vectors
    ///                         .into_iter() // convert the vector to an iterator
    ///                         .kmerge() // merge the iterators into a single sorted iterator
    ///                         .peekable() // make the first entry of the merged vector peekable
    ///                         .simplify( FieldFloat64::new() ); // combine entries with equal indices, and drop zeros
    /// 
    /// // verify
    /// assert!( sum.eq( vec![ (1,1.), (2,4.), (3,6.), (4,4.) ] ) )  
    /// ```
    /// 
    /// Add three vectors **sorted in descending order** using the the [hit_merge_by_predicate](crate::utilities::iterators::merge::hit::HitMergeByPredicateTrait::hit_merge_by_predicate) method of the [HitMergeByPredicateTrait](crate::utilities::iterators::merge::hit::HitMergeByPredicateTrait).
    /// 
    /// 
    /// ```
    /// //  import definitions
    /// use oat_rust::utilities::iterators::merge::hit::{HitMergeByPredicateTrait};
    /// use oat_rust::algebra::vectors::operations::VectorOperations;
    /// use oat_rust::algebra::rings::types::native::FieldFloat64;
    /// use oat_rust::utilities::order::{OrderOperatorAuto, ReverseOrder};
    /// 
    /// //  define a "predicate" object to judge the order of elements
    /// let order_operator = OrderOperatorAuto::new(); // implments the default order on Rust objects
    /// let order_operator = ReverseOrder::new( order_operator ); // reverses the order
    /// 
    /// // add
    /// let vectors     =   vec![ vec![ (2,2.), (1,1.) ], vec![ (3,3.), (2,2.) ], vec![ (4,4.), (3,3.) ] ];
    /// let sum         =   vectors
    ///                         .into_iter() // convert the vector to an iterator
    ///                         .hit_merge_by_predicate(order_operator) // merge the iterators into a single sorted iterator
    ///                         .peekable() // make the first entry of the merged vector peekable
    ///                         .simplify( FieldFloat64::new() ); // combine entries with equal indices, and drop zeros
    /// 
    /// // verify
    /// assert!( sum.eq( vec![ (4,4.), (3,6.), (2,4.), (1,1.) ] ) )  
    /// ```
    /// 
    /// For more examples on how to combine iterators, check out the [merge](crate::utilities::iterators::merge) module.
    /// 
    /// # See also
    /// 
    /// The method [VectorOperations::sum] calls on a iterator of vectors and returns a vector. Contrast this with
    /// 
    /// - [VectorOperations::add], which is a method called on a vector, and returns a vector
    /// - [VectorOperations::sum_coefficients], which is a method called on a vector, and returns a scalar
    fn sum_vectors< RingOperator, OrderOperator, RingElement, Index >( self, ring_operator: RingOperator, order_operator: OrderOperator )
            -> 
            Simplify<
                    IteratorsMergedInSortedOrder< <Self::Item as IntoIterator >::IntoIter, OrderOperator >,
                    RingOperator,
                >

        where 
            Self:                                   Sized,
            Self:                                   IntoIterator,
            Self::Item:                             IntoIterator,        
            < Self::Item as IntoIterator >::Item:   PartialEq + KeyValSet < Key = Index, Val = RingElement >,
            OrderOperator:      JudgePartialOrder<
                                        < <Self as IntoIterator >::Item as IntoIterator >::Item 
                                    >,
            RingOperator:       SemiringOperations< Element = RingElement >,
            Index:              PartialEq,
    {
        self.hit_merge_by_predicate( order_operator ).simplify( ring_operator )
    }    

}


// We implement this trait automatically on all iterators.
impl    < T >

    MultiVectorOperations for T  
{} // everything implemented automatically






































//  ---------------------------------------------------------------------------
//  TEST
//  ---------------------------------------------------------------------------


#[cfg(test)]
mod tests {
    use crate::algebra::vectors::operations::MultiVectorOperations;

    

    #[test]
    fn test_subtract_with_trait() {
        use crate::algebra::vectors::operations::VectorOperations;
        use crate::algebra::rings::types::native::FieldFloat64;
        use crate::utilities::order::OrderOperatorAuto;

        // define inputs
        let order_operator  =   OrderOperatorAuto::new();
        let ring_operator   =   FieldFloat64::new();
        let a               =   vec![ (1,1.), (2,2.) ].into_iter().peekable();
        let b               =   vec![ (2,2.), (3,3.) ].into_iter().peekable();

        // add
        let sum             =   a.subtract( b, ring_operator, order_operator );

        // verify
        assert!( sum.eq( vec![ (1,1.), (3,-3.) ] ) )        
    }    
    
    #[test]
    fn test_add_2_with_trait() {
        use crate::algebra::vectors::operations::VectorOperations;
        use crate::algebra::rings::types::native::FieldFloat64;
        use crate::utilities::order::OrderOperatorAuto;

        // define inputs
        let order_operator  =   OrderOperatorAuto::new();
        let ring_operator   =   FieldFloat64::new();
        let a               =   vec![ (1,1.), (2,2.) ].into_iter().peekable();
        let b               =   vec![ (2,2.), (3,3.) ].into_iter().peekable();

        // add
        let sum             =   a.add( b, ring_operator, order_operator );

        // verify
        assert!( sum.eq( vec![ (1,1.), (2,4.), (3,3.) ] ) )        
    }

    #[test]
    fn test_add_k_with_trait() {
        
        use crate::algebra::rings::types::native::FieldFloat64;
        use crate::utilities::order::OrderOperatorAuto;

        // define inputs
        let order_operator  =   OrderOperatorAuto::new();
        let ring_operator   =   FieldFloat64::new();
        let a               =   vec![ (1,1.), (2,2.) ].into_iter().peekable();
        let b               =   vec![ (2,2.), (3,3.) ].into_iter().peekable();
        let c               =   vec![ (3,3.), (4,4.) ].into_iter().peekable();

        // add
        let sum             =   [a,b,c].sum_vectors( ring_operator, order_operator );

        // verify
        assert!( sum.eq( vec![ (1,1.), (2,4.), (3,6.), (4,4.) ] ) )        
    }    

    #[test]
    fn test_add_2() {

        //  import definitions
        use itertools::Itertools;
        use crate::algebra::vectors::operations::VectorOperations;
        use crate::algebra::rings::types::native::FieldFloat64;

        // add
        let sum         =   vec![ (1,1.), (2,2.) ]
                                .into_iter() // convert the vector to an iterator
                                .merge( vec![ (2,2.), (3,3.) ] ) // merge with another vector
                                .peekable() // make the first entry of the merged vector peekable
                                .simplify( FieldFloat64::new() ); // combine entries with equal indices, and drop zeros
        
        // verify
        assert!( sum.eq( vec![ (1,1.), (2,4.), (3,3.) ] ) )
    }

    #[test]
    fn test_add_k() {
        //  import definitions
        use itertools::Itertools;
        use crate::algebra::vectors::operations::VectorOperations;
        use crate::algebra::rings::types::native::FieldFloat64;        

        // add
        let vectors     =   vec![ vec![ (1,1.), (2,2.) ], vec![ (2,2.), (3,3.) ], vec![ (3,3.), (4,4.) ] ];
        let sum         =   vectors
                                .into_iter() // convert the vector to an iterator
                                .kmerge() // merge the iterators into a single sorted iterator
                                .peekable() // make the first entry of the merged vector peekable
                                .simplify( FieldFloat64::new() ); // combine entries with equal indices, and drop zeros

        // verify
        assert!( sum.eq( vec![ (1,1.), (2,4.), (3,6.), (4,4.) ] ) )  
    }

    #[test]
    fn test_add_k_with_hit() {
        //  import definitions
        use crate::utilities::iterators::merge::hit::HitMergeByPredicateTrait;
        use crate::algebra::vectors::operations::VectorOperations;
        use crate::algebra::rings::types::native::FieldFloat64;
        use crate::utilities::order::{OrderOperatorAuto, ReverseOrder};

        //  define a "predicate" object to judge the order of elements
        let order_operator = OrderOperatorAuto::new(); // implments the default order on Rust objects
        let order_operator = ReverseOrder::new( order_operator ); // reverses the order

        // add
        let vectors     =   vec![ vec![ (2,2.), (1,1.) ], vec![ (3,3.), (2,2.) ], vec![ (4,4.), (3,3.) ] ];
        let sum         =   vectors
                                .into_iter() // convert the vector to an iterator
                                .hit_merge_by_predicate(order_operator) // merge the iterators into a single sorted iterator
                                .peekable() // make the first entry of the merged vector peekable
                                .simplify( FieldFloat64::new() ); // combine entries with equal indices, and drop zeros

        // verify
        assert!( sum.eq( vec![ (4,4.), (3,6.), (2,4.), (1,1.) ] ) )  
    }    
  
    #[test]
    fn test_subtract() {
        //  import definitions
        use itertools::Itertools;
        use crate::algebra::vectors::operations::VectorOperations;
        use crate::algebra::rings::types::native::FieldFloat64;

        //  define the terms
        let a   =   vec![ (1,1.), (2,2.) ];
        let b   =   vec![ (2,2.), (3,3.) ];

        //  subtract b from a
        let sum         =   a
                                .into_iter() // convert the vector to an iterator
                                .merge( // merge with negative b
                                    b.into_iter() // convert b to an iterator
                                        .negate( FieldFloat64::new() ) // convert b to negative b
                                ) 
                                .peekable() // make the first entry of the merged vector peekable (a technical requirement)
                                .simplify( FieldFloat64::new() ); // combine entries with equal indices, and drop zeros

        //  verify
        assert!( sum.eq( vec![ (1,1.), (3,-3.) ] ) )        
    }

    #[test]
    fn test_basic_transforms() {

        use crate::algebra::rings::types::native::FieldFloat64;
        use crate::algebra::vectors::operations::VectorOperations;

        // DEFINE A COEFFICIENT RING
        let ring_operator = FieldFloat64::new(); // the ring of real numbers, represented by floats

        // DEFINE A SEQUENCE OF VECTOR ENTRIES
        let entry_data = vec![ (1, 1.), (2, 2.), (3, 3.), (3, 3.), (4, 0.) ];
    
        // DEFINE A SPARSE VECTOR ITERATOR (i.e., an iterator that runs over entries)       
        let sparse_vec = entry_data.iter().cloned();
    
        // SCALE THE VECTOR BY 2.
        let scaled : Vec<_> = sparse_vec
                                .clone() // this makes a copy of the iterator, so the original stays unchanged
                                .scale_by( 2., ring_operator )
                                .collect(); // this collects the entries of the iterator into a standard Rust vector
        assert_eq!( scaled, vec![ (1, 2.), (2, 4.), (3, 6.), (3, 6.), (4, 0.) ]);
    
        // DROP ZERO ENTRIES
        let dropped : Vec<_> = sparse_vec
                                .clone() // this makes a copy of the iterator, so the original stays unchanged
                                .drop_zeros( ring_operator )
                                .collect(); // this collects the entries of the iterator into a standard Rust vector
        
        assert_eq!( dropped, vec![ (1, 1.), (2, 2.), (3, 3.), (3, 3.) ]);
    
        // MERGE CONSECUTIVE ENTRIES THAT SHARE THE SAME INDEX
        let gathered : Vec<_> = sparse_vec
                                .clone() // this makes a copy of the iterator, so the original stays unchanged
                                .peekable() // this puts the iterator in a slightly different form, which is compatible with gather
                                .gather( ring_operator )
                                .collect(); // this collects the entries of the iterator into a standard Rust vector
        assert_eq!( gathered, vec![ (1, 1.), (2, 2.), (3, 6.), (4, 0.) ]);  
        
        // MERGE CONSECUTIVE ENTRIES THAT SHARE THE SAME INDEX, AND DROP RESULTING ZEROS
        let gathered : Vec<_> = sparse_vec
                                .clone() // this makes a copy of the iterator, so the original stays unchanged
                                .peekable() // this puts the iterator in a slightly different form, which is compatible with gather
                                .simplify( ring_operator )
                                .collect(); // this collects the entries of the iterator into a standard Rust vector
        assert_eq!( gathered, vec![ (1, 1.), (2, 2.), (3, 6.), ]);          
        
    }     
    

    #[test]
    fn doctest_draft_test_filter_map() {    
        use std::collections::HashMap;
        use crate::algebra::vectors::operations::{FilterChangeIndex, };

        // create iterator
        let entry_iter_data: Vec< (usize, usize) > = vec![(1,1), (2,2), (3,3)];
        let entry_iter = entry_iter_data.iter();
        
        // hashmaps implement `EvaluateFunction` automatically; see documentaion for `EvaluateFunction` for details
        let mut hash = HashMap::new();
        hash.insert( 1usize, -1i32 );
        hash.insert( 2usize, -2i32 );

        // create a filter/mapped iterator        
        let iter_filtermapped_by_index: FilterChangeIndex<
                                            std::slice::Iter<'_, (usize, usize)>,
                                            HashMap<usize,i32>,
                                            i32
                                        > 
            = FilterChangeIndex::new( entry_iter, hash );

        // check that it is correct        
        itertools::assert_equal( iter_filtermapped_by_index, vec![(-1,1usize), (-2,2usize)] );
    }

    #[test]
    fn doctest_draft_test_negate() {
        use crate::algebra::vectors::operations::Negate;
        use crate::algebra::rings::types::field_prime_order::PrimeOrderField;

        let vec = vec![ (0, 0usize), (1, 1usize) ];
        let ring_operator = PrimeOrderField::new(7);
        let negated = Negate::new( vec.iter().cloned(), ring_operator );
        itertools::assert_equal( negated, vec![ (0, 0usize), (1, 6usize) ] );
    }

    #[test]
    fn test_multiply_matrix() {

        use crate::algebra::rings::types::native::FieldFloat64;
        use crate::algebra::vectors::operations::VectorOperations;
        use crate::utilities::order::OrderOperatorAuto;
        
        // define inputs
        let v               =   vec![     (0, 1.), (1, 1.)    ];
        let data =   vec![                         
                                    vec![ (0, 1.), (1, 1.)    ],
                                    vec![          (1, 1.)    ],
                                    vec![                     ],  
                                ];
        let matrix          =   |i| { let row: &Vec<_> = &data[i]; row.clone() };
        let ring_operator   =   FieldFloat64::new();
        let order_operator  =   OrderOperatorAuto::new();
        
        // multiply
        let u: Vec<_>               =   v.multiply_with_matrix_fnmut_unsimplified( matrix, ring_operator, order_operator ).collect();

        println!("{:?}", & u );
        assert!(    u.into_iter().eq( vec![ (0, 1.), (1, 1.), (1, 1.) ] )       );        
    }




    #[test]
    fn test_misc_in_operations() {   
    }



}


