//! Simplices represented by Rust vectors

use derive_new::new;


use crate::utilities::iterators::general::TwoTypeIterator;
use crate::utilities::iterators::merge::hit::{HitMergeByPredicateTrait, IteratorsMergedInSortedOrder};
use crate::utilities::order::{OrderOperatorAutoReverse, JudgeOrder, JudgePartialOrder};
use crate::utilities::sequences_and_ordinals::{CombinationsReverse, SortedVec};
use itertools::{Dedup, KMerge, Itertools, Combinations};
use std::cmp::Ordering;
use std::iter::{Flatten, Cloned};
use std::slice::Iter;

use std::fmt::Debug;



//  ================================================================================
//  TWIST ORDER OPERATOR - UNFILTERED SIMPLEX
//  ================================================================================

/// Orders simplices first by dimension (descending) then by lexicographic order (ascending)
/// 
/// This order is specially designed to facilitate the `twist` optimization in 
/// persistent homology calculations.
/// 
/// This object implements the `JudgeOrder` and `JudgePartialOrder` traits.
#[derive(Clone,Debug,new)]
pub struct OrderOperatorTwistSimplex;

impl < Vertex: Ord >

    JudgePartialOrder
        < Vec< Vertex > > for 

    OrderOperatorTwistSimplex
{
    fn judge_partial_cmp
        ( & self, lhs: & Vec< Vertex >, rhs: & Vec< Vertex > ) -> Option<Ordering>
    {   Some( self.judge_cmp(lhs,rhs) )      }
}


impl < Vertex: Ord >

    JudgeOrder
        < Vec< Vertex > > for 

    OrderOperatorTwistSimplex
{
    fn judge_cmp
        ( & self, lhs: & Vec< Vertex >, rhs: & Vec< Vertex > ) -> Ordering
    {
        // compare simplex dimensions
        let comp = rhs.len().cmp( & lhs.len() ); // note we put RHS on the LEFT here
        if comp != Ordering::Equal { return comp }
        // then compare simplices lexicographically
        lhs.cmp( rhs )        
    }
}


//  ================================================================================
//  Insertion
//  ================================================================================


/// Inserts a vertex into a simplex, and returns a new simplex and the index of the inserted vertex.
/// 
/// For example, in the simplex is `[0, 2, 3]` and the vertex is `1`, the result will be `([0, 1, 2, 3], 1)`.
/// 
/// **Note** This function does not check that the simplex is strictly sorted or that the vertex is not already in the simplex.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::topology::simplicial::simplices::vector::insert_vertex;
/// let simplex = vec![0, 2, 3];
/// let (new_simplex, index) = insert_vertex(1, &simplex);
/// assert_eq!(new_simplex, vec![0, 1, 2, 3]);
/// assert_eq!(index, 1);
/// ```
/// 
/// ```
/// use oat_rust::topology::simplicial::simplices::vector::insert_vertex;
/// let simplex = vec![];
/// let (new_simplex, index) = insert_vertex(1, &simplex);
/// assert_eq!(new_simplex, vec![1]);
/// assert_eq!(index, 0);
/// ```
pub fn insert_vertex< Vertex >( vertex: Vertex, simplex: & Vec< Vertex > )
       -> ( Vec< Vertex >, usize )
    where Vertex: Ord + Clone
{
    // if the simplex is empty, return a new simplex with the vertex
    if simplex.is_empty() { 
        return ( vec![ vertex ], 0 ) 
    }

    // otherwise, find the insertion point
    let mut i = 0;
    while i < simplex.len() && simplex[ i ] < vertex {
        i += 1;
    }
    
    // create a new simplex with the vertex inserted at position i
    let mut new_simplex = Vec::with_capacity( simplex.len() + 1 );
    new_simplex.extend_from_slice( &simplex[ 0 .. i ] ); // copy the first i elements
    new_simplex.push( vertex ); // insert the new vertex
    new_simplex.extend_from_slice( &simplex[ i .. ] ); // copy the rest of the elements
    
    ( new_simplex, i )
}



//  ================================================================================
//  FACETS OF SIMPLICES
//  ================================================================================



//  ---------------------------------------------------------------------------
//  FACES FROM FACETS-OF-THE-COMPLEX (NEW)
//  ---------------------------------------------------------------------------



/// Assuming the user-provided facets have vertices sorted in ascending order, returns
/// an iterator that runs over `dim`-dimensional subsimplices in **descending** lexicographic 
/// order.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::topology::simplicial::simplices::vector::dimension_d_simplices_in_reverse_lexicographic_order_iter;
/// use oat_rust::utilities::sequences_and_ordinals::SortedVec;
/// 
/// let dowker_simplices    =   vec![ SortedVec::new( vec![0,1,2] ).unwrap(), 
///                                   SortedVec::new( vec![0,3]   ).unwrap(),  ];
/// let one_simplices       =   dimension_d_simplices_in_reverse_lexicographic_order_iter( &dowker_simplices, 1 );
/// 
/// itertools::assert_equal( one_simplices, vec![ vec![1,2], vec![0,3], vec![0,2], vec![0,1] ] );
/// ```
pub fn dimension_d_simplices_in_reverse_lexicographic_order_iter< Vertex >(
    facets: &Vec< SortedVec< Vertex >>, 
    dim: isize 
) 
    -> 
    TwoTypeIterator<
        std::iter::Empty< Vec<Vertex> >,
        Dedup< 
            IteratorsMergedInSortedOrder< 
                CombinationsReverse< Vertex, &Vec< Vertex > >,
                OrderOperatorAutoReverse,
            >,         
        >,
    >
    
    where Vertex: Ord + Clone
{
    // if dimemsion is < -1 return an empty sequence
    if dim < -1 { 
        return  TwoTypeIterator::Version1(
                    std::iter::empty()
                )
    }

    // otherwise return a sequence of faces   
    let iter = 
            facets
                .iter()
                .map( |vertices| CombinationsReverse::from_vec( dim as usize +1, vertices.vec() ) )
                .hit_merge_by_predicate( OrderOperatorAutoReverse )
                .dedup();
    return TwoTypeIterator::Version2(iter)
}

/// Returns
/// an iterator that runs over `dim`-dimensional subsimplices in **ascending** lexicographic 
/// order.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::topology::simplicial::simplices::vector::dimension_d_simplices_in_lexicographic_order_iter;
/// use oat_rust::utilities::sequences_and_ordinals::SortedVec;
/// 
/// let dowker_simplices    =   vec![ SortedVec::new( vec![0,1,2] ).unwrap(), 
///                                   SortedVec::new( vec![0,3]   ).unwrap(),  ];
/// let one_simplices       =   dimension_d_simplices_in_lexicographic_order_iter( &dowker_simplices, 1 );
/// 
/// itertools::assert_equal( one_simplices, vec![ vec![0,1], vec![0,2], vec![0,3], vec![1,2], ] );
/// ```
pub fn dimension_d_simplices_in_lexicographic_order_iter< Vertex >(
    facets: & Vec< SortedVec< Vertex >>, 
    dim: isize 
) 
    -> 
    TwoTypeIterator<
        std::iter::Empty< Vec<Vertex> >,
        Dedup< KMerge<  Combinations<Cloned<Iter<Vertex>>> > >,
    >
    where Vertex: Ord + Clone
{
    // if dimension is less than -1, return an emtpy iterator
    if dim < -1 { 
        return  TwoTypeIterator::Version1(
                    std::iter::empty()
                ) 
    }    
    // otherwise return 
    let iter    =   
            facets
                .iter()
                .map( |x| x.vec().iter().cloned().combinations( (dim + 1isize) as usize )  )
                .kmerge()
                .dedup();
    return TwoTypeIterator::Version2(
        iter
    )
}



/// Returns
/// an iterator that runs over all dimension-k subsimplices, for 0 ≤ k ≤ d,
/// sorted first by dimension, then lexicographically, in ascending order.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::topology::simplicial::simplices::vector::dimension_0_through_d_simplices_in_dimensionwise_lexicographic_order_iter;
/// use oat_rust::utilities::sequences_and_ordinals::SortedVec;
/// 
/// let dowker_simplices    =   vec![ SortedVec::new( vec![0,1,2] ).unwrap(), 
///                                   SortedVec::new( vec![0,3]   ).unwrap(),  ];
/// let all_simplices       =   dimension_0_through_d_simplices_in_dimensionwise_lexicographic_order_iter( &dowker_simplices, 1 );
/// 
/// itertools::assert_equal( all_simplices, vec![ vec![0], vec![1], vec![2], vec![3], vec![0,1], vec![0,2], vec![0,3], vec![1,2], ] );
/// ```
pub fn dimension_0_through_d_simplices_in_dimensionwise_lexicographic_order_iter< Vertex >( 
    facets: & Vec< SortedVec< Vertex >>, 
    max_dim: isize 
) 
    -> 
    Flatten<
        std::vec::IntoIter<
            TwoTypeIterator<
                std::iter::Empty< Vec<Vertex> >,
                Dedup< KMerge<  Combinations<Cloned<Iter<Vertex>>> > >,
            >
        >
    >
    where
        Vertex:         Ord + Clone,
{
    (0 .. max_dim + 1)
        .map(   |dim|
                dimension_d_simplices_in_lexicographic_order_iter(facets, dim)
            )
        .collect_vec()
        .into_iter()
        .flatten()
}

/// Returns
/// an iterator that runs over all dimension-k subsimplices, for 0 ≤ k ≤ d,
/// sorted first by dimension (strictly descending order), then lexicographically (strictly descending order).
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::topology::simplicial::simplices::vector::dimension_0_through_d_simplices_in_reverse_dimensionwise_lexicographic_order_iter;
/// use oat_rust::utilities::sequences_and_ordinals::SortedVec;
/// 
/// let dowker_simplices    =   vec![ SortedVec::new( vec![0,1,2] ).unwrap(), 
///                                   SortedVec::new( vec![0,3]   ).unwrap(),  ];
/// let all_simplices       =   dimension_0_through_d_simplices_in_reverse_dimensionwise_lexicographic_order_iter( &dowker_simplices, 1 );
/// 
/// itertools::assert_equal( all_simplices, vec![ vec![1,2], vec![0,3], vec![0,2], vec![0,1], vec![3], vec![2], vec![1], vec![0], ] );
/// ```
pub fn dimension_0_through_d_simplices_in_reverse_dimensionwise_lexicographic_order_iter< Vertex >( 
    facets: &Vec< SortedVec< Vertex >>, 
    max_dim: isize 
) 
    -> 
    Flatten< std::vec::IntoIter<
        TwoTypeIterator<
            std::iter::Empty< Vec<Vertex> >,
            Dedup< 
                IteratorsMergedInSortedOrder< 
                    CombinationsReverse< Vertex, &Vec< Vertex > >,
                    OrderOperatorAutoReverse,
                >,  
            >        
        >,   
    >>
    where
        Vertex:         Ord + Clone,
{
    (0 .. max_dim + 1)
        .rev()
        .map(   |dim|
                dimension_d_simplices_in_reverse_lexicographic_order_iter( facets, dim )
            )
        .collect_vec()
        .into_iter()
        .flatten()
}


/// Returns
/// an iterator that runs over all dimension-k subsimplices, for 0 ≤ k ≤ d,
/// sorted first by dimension (strictly ascending order), then lexicographically (strictly descending order).
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::topology::simplicial::simplices::vector::dimension_0_through_d_simplices_in_ascending_dimension_descending_lexicographic_order_iter;
/// use oat_rust::utilities::sequences_and_ordinals::SortedVec;
/// 
/// let dowker_simplices    =   vec![ SortedVec::new( vec![0,1,2] ).unwrap(), 
///                                   SortedVec::new( vec![0,3]   ).unwrap(),  ];
/// let all_simplices       =   dimension_0_through_d_simplices_in_ascending_dimension_descending_lexicographic_order_iter( &dowker_simplices, 1 );
/// 
/// itertools::assert_equal( all_simplices, vec![ vec![3], vec![2], vec![1], vec![0], vec![1,2], vec![0,3], vec![0,2], vec![0,1], ] );
/// ```
pub fn dimension_0_through_d_simplices_in_ascending_dimension_descending_lexicographic_order_iter< Vertex >( 
        facets: &Vec< SortedVec< Vertex >>, 
        max_dim: isize 
    ) 
    -> 
    Flatten< std::vec::IntoIter<
        TwoTypeIterator<
            std::iter::Empty< Vec<Vertex> >,
            Dedup< 
                IteratorsMergedInSortedOrder< 
                    CombinationsReverse< Vertex, &Vec< Vertex > >,
                    OrderOperatorAutoReverse,
                >,  
            >        
        >,   
    >>
    where
        Vertex:         Ord + Clone,
{
    (0 .. max_dim + 1)
        .map(   |dim|
                dimension_d_simplices_in_reverse_lexicographic_order_iter( facets, dim )
            )
        .collect_vec()
        .into_iter()
        .flatten()
}



/// Returns
/// a vector of vectors whose kth entry contains all dimension-k subsimplices,
/// sorted in in lexicographic order.
pub fn dimension_0_through_d_simplices_in_lexicographic_order_binned_by_dimension< Vertex >( 
            facets: & Vec< SortedVec< Vertex >>, 
            max_dim: isize 
        ) 
        -> 
        Vec< Vec< Vec< Vertex >>> 
        where Vertex: Ord + Clone
{
    if max_dim < 0 { return Vec::new() } // return the empty set if max_dim < 0
    let mut seq             =   Vec::with_capacity( max_dim as usize );
    for dim in 0 .. max_dim + 1  {
        let vec: Vec<_>     =   dimension_d_simplices_in_lexicographic_order_iter(
                                    facets,
                                    dim
                                )
                                .collect();
        seq.push( vec );
    }
    seq
}

/// Similar to `dimension_0_through_d_simplices_in_lexicographic_order_binned_by_dimension`, but stores all simplices in a single vector, sorted first by dimension and second by lexicographic order.
pub fn dimension_0_through_d_simplices_in_dimensionwise_lexicographic_order_list< Vertex >( 
    facets: & Vec< SortedVec< Vertex >>, 
    max_dim: isize 
) 
-> 
Vec< Vec< Vertex >>
    where Vertex: Ord + Clone
{
    let mut a = dimension_0_through_d_simplices_in_lexicographic_order_binned_by_dimension( facets, max_dim );
    let mut b = Vec::new();
    for i in 0 .. a.len() {
        b.append( &mut a[ i ] );
    }
    b
}




//  ---------------------------------------------------------------------------
//  FACETS-OF-A-SIMPLEX
//  ---------------------------------------------------------------------------

/// Maintains an "inner state" that steps through the facets of a simplex in 
/// ascending lexicographic order; only returns `Some(())` or `None`.
/// 
/// This object can be useful in applications that may require access to a sequence of
/// simplices by reference, but don't need to own the referenced simplex.
/// 
/// # Examples
/// 
/// ```

/// use oat_rust::topology::simplicial::simplices::vector::FacetIteratorNoReturnAscendingLex;

/// // Create the iterator
/// let mut facet_iterator_noreturn     
///     =   FacetIteratorNoReturnAscendingLex::new( vec![0, 1, 2], None );
///
/// // Test it                                                
/// let mut answers = vec![
///     FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2] , facet: vec![0, 1], deleted_vertex_index: Some(2) },
///     FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2] , facet: vec![0, 2], deleted_vertex_index: Some(1) },
///     FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2] , facet: vec![1, 2], deleted_vertex_index: Some(0) },
///     FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2] , facet: vec![]    , deleted_vertex_index: None    },            
///     FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2] , facet: vec![0, 1], deleted_vertex_index: Some(2) },                                                                        
/// ];
///
/// for i in 0..5 {
///     facet_iterator_noreturn.next();
///     assert_eq!( &facet_iterator_noreturn, &answers[ i ] )    
/// }      
//
/// // Re-initialize with a new simplex
///
/// facet_iterator_noreturn.reinitialize_with_simplex( vec![0 ,3] );
///
/// // Test again        
///
/// answers = vec![
///     FacetIteratorNoReturnAscendingLex { simplex: vec![0, 3], facet: vec![0], deleted_vertex_index: Some(1) },
///     FacetIteratorNoReturnAscendingLex { simplex: vec![0, 3], facet: vec![3], deleted_vertex_index: Some(0) },
///     FacetIteratorNoReturnAscendingLex { simplex: vec![0, 3], facet: vec![] , deleted_vertex_index: None    },
///     FacetIteratorNoReturnAscendingLex { simplex: vec![0, 3], facet: vec![0], deleted_vertex_index: Some(1) },                                                                                                      
/// ];    
///
/// for i in 0..4 {
///     facet_iterator_noreturn.next();
///     assert_eq!( &facet_iterator_noreturn, &answers[ i ] )    
/// }   
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct FacetIteratorNoReturnAscendingLex< Vertex: Ord >
{
    pub simplex: Vec< Vertex > ,
    pub facet: Vec< Vertex >,
    pub deleted_vertex_index: Option<usize>
}

impl < Vertex: Ord > FacetIteratorNoReturnAscendingLex < Vertex > 
    where Vertex: Clone
{

    /// Initialize a no-return facet iterator.  
    /// 
    /// When it is first constructed, the internal state of the iterator does not represent a facet.
    /// The internal state will be updated to represent the first facet when `next()` is called for the first time.
    pub fn new( simplex: Vec< Vertex >, buffer: Option< Vec<Vertex> > ) -> Self {
        let buffer = buffer.unwrap_or( Vec::with_capacity( simplex.len()-1 ) );
        
        FacetIteratorNoReturnAscendingLex {
            simplex,
            facet: buffer,
            deleted_vertex_index: None
        }
    }

    /// Reinitialize a no-return facet iterator with a new simplex.
    /// 
    /// The result is essentially the same as replacing our iterator with `FacetIteratorNoReturnAscendingLex::new( simplex )`.
    pub fn reinitialize_with_simplex( &mut self, simplex: Vec< Vertex > ) {
        // if necessary, expand the capacity of the facet vector
        if simplex.len() > 1 + self.facet.capacity() { 
            self.facet.reserve_exact(
                simplex.len() - 1 - self.facet.capacity()
            ) 
        }
        // replace the old simplex with the new
        self.simplex    =   simplex;
        // update the state to indicate that it does not represent a facet
        self.facet.clear();
        self.deleted_vertex_index = None;
    }
}


impl < Vertex: Ord >
    Iterator for 
    FacetIteratorNoReturnAscendingLex < Vertex >     
    where Vertex : Clone
{
    type Item    =   ();

    fn next( &mut self ) -> Option<()> {

        if let Some( deleted_index ) = self.deleted_vertex_index {

            if deleted_index == 0 {
                // if we start from the facet obtained by deleting vertex 0, then the 
                // next state should **not** represent a facet
                self.deleted_vertex_index   =   None;
                self.facet.clear();
                None
                
            } else {
                // if we start from the facet obtained by deleting vertex k > 0, then 
                // the next state should represent the facet obtained by deleting vertex k-1
                let next_deleted_index  =   deleted_index - 1;
                self.facet[ next_deleted_index ] = self.simplex[ deleted_index ].clone(); // replace the deleted index and remove the next one
                self.deleted_vertex_index = Some( next_deleted_index );
                Some( () )
            }
        
        } else {

            self.facet.clear();
            for i in 0..self.simplex.len()-1 {   // dim = 1 less than num_vertices = INDEX OF LAST VERTEX IN SIMPLEX
                self.facet.push( self.simplex[ i ].clone() ) 
            }      
            // set deleted vertex equal to last
            self.deleted_vertex_index = Some( self.simplex.len()-1 );  // dim = 1 less than num_vertices = INDEX OF LAST VERTEX IN SIMPLEX
            Some(())            
        }
    }
}









//  ===========================================================================
//  ===========================================================================
//  TESTS
//  ===========================================================================
//  ===========================================================================





#[cfg(test)]
mod tests {
    use crate::utilities::order::{is_sorted_strictly, OrderOperatorAuto};

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;


    //  =======================================================================   
    //  SMALL EXAMPLES / EXACT RESULTS
    //  =======================================================================     


    //  ENUMERATION METHODS
    //  =======================================================================

    #[test]
    fn test_ordered_subsimplices_up_thru_dim() {

        let facets          =   vec![ 
                                                                SortedVec::new(vec![0, 1, 2]).unwrap()
                                                            ];

        assert_eq!(         dimension_0_through_d_simplices_in_lexicographic_order_binned_by_dimension( & facets, 2),
                            vec![
                                vec![   vec![0],     vec![1],    vec![2]         ],                                
                                vec![   vec![0,1],   vec![0,2],  vec![1,2]       ],
                                vec![   vec![0,1,2]                              ]
                            ]                            
        );

        assert_eq!(         dimension_0_through_d_simplices_in_dimensionwise_lexicographic_order_list( & facets, 2),
                            vec![
                                        vec![0],     vec![1],    vec![2],                                         
                                        vec![0,1],   vec![0,2],  vec![1,2],       
                                        vec![0,1,2]                              
                            ]
        ) ;       


    }

    //  SIMPLEX ITERATOR
    //  =======================================================================    

    #[test]
    fn test_ascending_facet_iterator_no_return()
    {

        // Create the iterator
        let mut facet_iterator_noreturn     =   FacetIteratorNoReturnAscendingLex::new(
                                                    vec![0, 1, 2],
                                                    None
                                                );

        // Test it                                                
        let mut answers = vec![
            FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2], facet: vec![0, 1], deleted_vertex_index: Some(2) },
            FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2], facet: vec![0, 2], deleted_vertex_index: Some(1) },
            FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2], facet: vec![1, 2], deleted_vertex_index: Some(0) },
            FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2], facet: vec![]    , deleted_vertex_index: None    },            
            FacetIteratorNoReturnAscendingLex { simplex: vec![0, 1, 2], facet: vec![0, 1], deleted_vertex_index: Some(2) },                                                                        
        ];

        for i in 0..5 {
            facet_iterator_noreturn.next();
            assert_eq!( &facet_iterator_noreturn, &answers[ i ] )    
        }      

        // Re-initialize with a new simplex

        facet_iterator_noreturn.reinitialize_with_simplex( vec![0 ,3] );

        // Test again        

        answers = vec![
            FacetIteratorNoReturnAscendingLex { simplex: vec![0, 3], facet: vec![0], deleted_vertex_index: Some(1) },
            FacetIteratorNoReturnAscendingLex { simplex: vec![0, 3], facet: vec![3], deleted_vertex_index: Some(0) },
            FacetIteratorNoReturnAscendingLex { simplex: vec![0, 3], facet: vec![] , deleted_vertex_index: None    },
            FacetIteratorNoReturnAscendingLex { simplex: vec![0, 3], facet: vec![0], deleted_vertex_index: Some(1) },                                                                                                      
        ];    

        for i in 0..4 {
            facet_iterator_noreturn.next();
            assert_eq!( &facet_iterator_noreturn, &answers[ i ] )    
        }     
                
    }    






    //  =======================================================================   
    //  AUTOMATED TEST OF ALL ENUMERATION TECHNIQUES (EXCLUDING ITERATORS)
    //  =======================================================================   


    //  EXECUTE THE TESTS
    //  =======================================================================   

    fn test_all_enumeration_techniques() {

        use crate::utilities::random::random_sequences;

        let size_of_ambient_set             =   12;
        let cardinalities                               =   0 .. 13;
        let max_dim                                     =   5;

        for probability in [0.0, 0.01, 0.05 ] {
            let facets                                  =   random_sequences( 20, cardinalities.clone(), probability);
            verify_simplex_enumeartion_methods_are_consistent( &facets, max_dim )
        }

    }    


    //  HELPER FUNCTIONS
    //  =======================================================================      


    //  VALIDATE RESULTS
    //  -----------------------------------------------------------------------    

    /// Given a collection of facets, verifies that all the functions in this module which produce
    /// lists of simplices return the same results.
    fn verify_simplex_enumeartion_methods_are_consistent( 
            facets: &Vec< SortedVec< usize >>, 
            max_dim: isize,            
        ) 
    {

        // Check mutual consistency
        // --------------------------------------------------------------------

        let simplices: Vec<_>                           =   dimension_0_through_d_simplices_in_dimensionwise_lexicographic_order_iter( facets, max_dim ).collect();        

        // test 1
        let test                                        =   dimension_0_through_d_simplices_in_dimensionwise_lexicographic_order_list(facets, max_dim);
        assert_eq!( & simplices, & test );

        // test 2
        let test                                        =   dimension_0_through_d_simplices_in_lexicographic_order_binned_by_dimension(facets, max_dim)
                                                                .into_iter()
                                                                .flatten()
                                                                .collect::<Vec<_>>();
        assert_eq!( & simplices, & test );

        // test 3
        let mut test                                    =   dimension_0_through_d_simplices_in_ascending_dimension_descending_lexicographic_order_iter(facets, max_dim)
                                                                .collect::<Vec<_>>();

        ( &mut test ).sort_by(| a,b|  b.len().cmp(  &a.len()  )  ); // a stable sorting, so that simplices go into reverse order by length but otherwise remain in the same order
        ( &mut test).reverse();
        assert_eq!( & simplices, & test );        

        // test 4
        let mut test                                    =   dimension_0_through_d_simplices_in_reverse_dimensionwise_lexicographic_order_iter(facets, max_dim).collect::<Vec<_>>();
        ( &mut test).reverse();
        assert_eq!( & simplices, & test );

        // test 5 and 6
        let simplices_binned_by_dimension               =   dimension_0_through_d_simplices_in_lexicographic_order_binned_by_dimension(facets, max_dim);
        for dim in 0 .. max_dim + 1 {
            let test                                    =   dimension_d_simplices_in_lexicographic_order_iter( facets, dim ).collect_vec();
            assert_eq!( & simplices_binned_by_dimension[ dim as usize ], & test  );
            let mut test                                =   dimension_d_simplices_in_reverse_lexicographic_order_iter( facets, dim ).collect_vec();            
            ( &mut test).reverse();
            assert_eq!( & simplices_binned_by_dimension[ dim as usize ], & test  );            
        }

        // Check ordering of simplices
        // --------------------------------------------------------------------        

        let binned_simplices                            =   dimension_0_through_d_simplices_in_lexicographic_order_binned_by_dimension(facets, max_dim);
        for ( dim, simplices ) in binned_simplices.into_iter().enumerate() {
            
            // test 7: within a dimension, simplices are sorted lexicographically
            assert!( simplices.is_sorted() );

            // test 8: each simplex in the bin for dimension dim does indeed have dimension dim
            for simplex in simplices.iter() { 
                assert!( simplex.len() == dim + 1 ) 
            }        
        }

        // Check the simplices themselves
        // --------------------------------------------------------------------

        // test 9: each simplex is sorted strictly
        for simplex in simplices.iter() { 
            assert!( is_sorted_strictly( & simplex, & OrderOperatorAuto::new() ) ) 
        }   

        // test 10: each simplex is a subset of one of the facets
        for simplex in simplices.iter() { 

            let order_certified_simplex         =   SortedVec::new( simplex.clone() ).unwrap();
            let is_ok                                   =   facets
                                                                .iter()
                                                                .any( 
                                                                    |facet| 
                                                                    facet.contains_subset( & order_certified_simplex ) 
                                                                );
            assert!( is_ok );
        }            

        // test 11: every length ≤ max_dim + 1 subset of the facets is a simplex

        let simplices_var                               =   subsequences_up_to_length_m_multi_source( facets, (max_dim + 1) as usize );
        assert_eq!( & simplices, & simplices_var );

    }



    //  ENUMERATE SUBSETS (ALTERNATE METHOD)
    //  -----------------------------------------------------------------------

    /// This is a helper function which generates a list of all length ≤ m subsequences
    /// of multiple sequences.
    /// It uses a different technique from our other enumeration methods, so we can use
    /// it to test that our methods are working.
    fn subsequences_up_to_length_m_multi_source<T: Clone + Ord >(sequences: &Vec< SortedVec<T> >, m: usize) -> Vec<Vec<T>> {


        let mut simplices                               =   Vec::new();
        for sequence in sequences {
            let subsequences    =   subsequences_up_to_length_m_single_source(sequence.vec(), m);
            simplices.extend( subsequences.into_iter() );
        }
        simplices.sort();
        simplices.dedup();

        simplices
    }



    /// This is a helper function which generates a list o fall length ≤ m subsequences.
    /// It uses a different techqique from our other enumeration methods, so we can use
    /// it to test that our methods are working.
    fn subsequences_up_to_length_m_single_source<T: Clone>(vec: &Vec<T>, m: usize) -> Vec<Vec<T>> {
        let mut result = Vec::new();
        let mut current = Vec::new();
    
        generate_subsequences(vec, m, 0, &mut current, &mut result);
    
        result
    }


    /// A subroutine for the alternate function that generates a list of all subsequences
    fn generate_subsequences<T: Clone>(
        vec: &Vec<T>, 
        max_len: usize, 
        start: usize, 
        current: &mut Vec<T>, 
        result: &mut Vec<Vec<T>>
    ) {
        // If the current subsequence is not empty, add it to the result
        if !current.is_empty() {
            result.push(current.clone());
        }
    
        // Explore further only if the length of the current subsequence is less than `max_len`
        if current.len() < max_len {
            for i in start..vec.len() {
                // Include the current element in the subsequence
                current.push(vec[i].clone());
    
                // Recurse to generate subsequences that include the current element
                generate_subsequences(vec, max_len, i + 1, current, result);
    
                // Backtrack (remove the last element) to explore other possibilities
                current.pop();
            }
        }
    }
    







}