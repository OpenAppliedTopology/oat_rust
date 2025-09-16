//! The Dowker complex of a set relation.
//! 
//! Equivalently, a simplicial complex `K` represented by a 
//! subset of simplices `Q` containing every maximal simplex of `K`.
//! 
//! 
//! 
//! # Example
//! 
//! Here we compute a basis for homology in dimension 1,
//! for the simplicial complex on vertex set `{0,1,2,3}`
//! whose maximal faces are `{0,1,2}, {0,3}, {1,3}, {2,3}`.
//! It "looks" like a tetrahedron where the 3-face and 
//! all but one of the 2-faces have been removed.
//! The coefficient field is the finite field of order 3.
//! 
//! 
//! ```
//! use oat_rust::algebra::matrices::operations::umatch::differential::DifferentialUmatch;   
//! use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;
//!             
//! use oat_rust::topology::simplicial::from::relation::DowkerComplex;
//! use oat_rust::topology::simplicial::simplices::vector::{dimension_d_simplices_in_lexicographic_order_iter, dimension_d_simplices_in_reverse_lexicographic_order_iter}; 
//! use oat_rust::topology::simplicial::simplices::vector::{dimension_0_through_d_simplices_in_ascending_dimension_descending_lexicographic_order_iter}; 
//!     
//! use oat_rust::utilities::order::OrderOperatorAuto;        
//! use oat_rust::utilities::sequences_and_ordinals::SortedVec;
//! 
//! use itertools::Itertools;
//!     
//! // Parameters
//! // ----------
//!     
//! // Define the maximum homology dimensino we want to compute.
//! let min_homology_dimension                  =   0;
//! let max_homology_dimension                  =   2;
//!     
//! // Define the ring operator for the finite field of order 3.
//! // You can use this object to perform arithmetic operations, e.g., 
//! // to add 1 and 1, you can run `ring_operator.add(1,1)`.
//! let ring_operator          =   PrimeOrderField::new(3);
//!     
//! // We will build a dowker complex.
//! // A dowker complex is defined by a vertex set V and a family S
//! // of subsets of V.  A subset of V forms a simplex iff it is 
//! // a subset of some element of S.  We refer to the elements 
//! // of S as "dowker simplices".
//!     
//! // Each dowker simplex is represented by a SortedVec of vertices.
//! // We store the list of all such simplices inside a larger vector.
//! let dowker_simplices: Vec<SortedVec<usize>>
//!     =   vec![    
//!                 vec![0,1,2], 
//!                 vec![0,3], 
//!                 vec![1,3], 
//!                 vec![2,3]  
//!             ]
//!             .into_iter()
//!             .map( |x| SortedVec::new(x).unwrap() )  // we unwrap because `new` can return an error
//!             .collect_vec();
//!     
//! //  Boundary matrix
//! //  ---------------
//!     
//! // This is a lazy object that generates rows/columns of the boundary matrix, on demand.
//! let boundary_matrix = DowkerComplex::new( dowker_simplices.clone(), ring_operator.clone() );
//!     
//! //  Simplex iterators
//! //  -----------------
//!     
//! // An iterator that runs over all triangles in the complex, in ascending 
//! // lexicographic order
//! let triangles_ascending_order = dimension_d_simplices_in_lexicographic_order_iter( &dowker_simplices, 2);
//!     
//! // An iterator that runs over all edges in the complex, in descending 
//! // lexicographic order
//! let triangles_descending_order = dimension_d_simplices_in_reverse_lexicographic_order_iter( &dowker_simplices, 2);   
//!     
//! // An iterator that runs over simplices of dimension 0 through max_homology_dimension,
//! // ordered first by dimension (ascending) and second by lexicographic order (descending)
//! let row_indices = boundary_matrix.simplices_in_row_reduction_order( max_homology_dimension );
//! 
//! // row_indices contains a reference to boundary_matrix, which will cause problems. Instead,
//! // we can construct an iterator that runs over the same sequence of simplices, like so:
//! let row_indices     =   dimension_0_through_d_simplices_in_ascending_dimension_descending_lexicographic_order_iter(&dowker_simplices, max_homology_dimension);
//!     
//! //  Matrix factorization (used to compute homology)
//! //  -----------------------------------------------
//!     
//! // Factor the boundary matrix
//! let factored    =   DifferentialUmatch::new( 
//!                         boundary_matrix, 
//!                         min_homology_dimension,
//!                         max_homology_dimension,
//!                     );
//!                   
//! // Betti numbers  
//! // -------------
//! // 
//! // To extract betti numbers, the user has to provide a function that assigns a
//! // dimension to each index (i.e. to each simplex). This is just the length of 
//! // the simplex minus 1.
//! // The resulting object, `betti_numbers`, is a hashmap that maps dimensions to betti numbers.
//! let betti_numbers   =   factored.betti_numbers(); 
//! 
//! // This loop prints the betti numbers in dimensions 0 and 1
//! for dim in 0 .. 2 {
//!     println!(
//!             // we'll insert two values into this string
//!             "The betti number in dimension {:?} is {:?}.",
//!             dim,                        
//!             betti_numbers.get( & dim )  
//!         );            
//! }
//! println!(""); // insert line break
//! 
//! // This should print the following:
//! // The betti number in dimension 0 is 1.
//! // The betti number in dimension 1 is 2.
//! 
//! // Cycle representatives for homology
//! // ----------------------------------
//! 
//! println!(
//!         "The following are basis vectors for homology in dimensions 0 through {:?}",
//!         max_homology_dimension,
//!     );
//! for (cycle_number, cycle) in factored.homology_basis().enumerate() {
//!     // `cycle` is an iterator.  For convenience, collect the elements of the
//!     // iterator into a Rust vector.
//!     let cycle: Vec<_> = cycle.collect();
//!     println!("Basis vector {:?} = {:?}", cycle_number, cycle);
//! }
//! 
//! // This should print the following:
//! //
//! // The following are basis vectors for homology in dimensions 0 through 2
//! // Basis vector 0 = [([0], 1)]
//! // Basis vector 1 = [([0, 2], 1), ([0, 3], 2), ([2, 3], 1)]
//! // Basis vector 2 = [([0, 1], 1), ([0, 3], 2), ([1, 3], 1)]
//! ```
//! 
//! ### Try changing the coefficient ring
//! 
//! OAT has a number of different [predefined coefficient rings](crate::algebra::rings::types), which you can substitute into
//! the example above in order to calculate homology with different coefficients.  Simply replace `PrimeOrderField::new(3)` 
//! in the following line
//! ```
//! # use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;       
//! let ring_operator   =   PrimeOrderField::new(3);
//! ``` 
//! with the [ring operator](crate::algebra::rings) of your choice.
//! 
//! ### Try changing the complex
//! 
//! These are the maximal faces in a (minimal) triangulation of the real projective plane:
//! 
//! ```
//! # use oat_rust::utilities::sequences_and_ordinals::SortedVec;
//! # use std::iter::FromIterator;
//! let dowker_simplices    =   Vec::from_iter(
//!     vec![    
//!             vec![1,2,6], vec![1,3,6], vec![1,4,5], 
//!             vec![1,2,5], vec![2,3,5], vec![2,4,6], 
//!             vec![1,3,4], vec![2,3,4], vec![3,5,6], 
//!             vec![4,5,6],  
//!         ]
//!         .into_iter()
//!         .map( |x| SortedVec::new(x) ) 
//! );
//! ```
//! 
//! # Implementation notes
//! 
//! This module uses the [SortedVec] data structure
//! for a variety of internal opertaions, e.g. binary search, set union, set difference, evaluation of
//! set containment, etc.

use crate::algebra::chain_complexes::ChainComplex;
use crate::algebra::matrices::query::{MatrixAlgebra, MatrixOracle, };
use crate::algebra::matrices::operations::MatrixOracleOperations;
use crate::algebra::rings::traits::RingOperations;
use crate::utilities::iterators::general::{symmetric_difference_of_ordered_iterators, IntersectOrderedIterators, TwoTypeIterator};
use crate::utilities::iterators::merge::hit::IteratorsMergedInSortedOrder;
use crate::utilities::sequences_and_ordinals::{SortedVec, CombinationsReverse};

use std::fmt::Debug;
use std::hash::Hash;
use std::iter::{ Cloned, Flatten};
use std::slice::Iter;

use itertools::{Combinations, Dedup, Itertools, KMerge};


use crate::topology::simplicial::boundary::{SimplexBoundaryAscend, SimplexBoundaryDescend};
use crate::topology::simplicial::simplices::vector::{dimension_0_through_d_simplices_in_ascending_dimension_descending_lexicographic_order_iter, dimension_0_through_d_simplices_in_dimensionwise_lexicographic_order_iter, dimension_d_simplices_in_lexicographic_order_iter, dimension_d_simplices_in_reverse_lexicographic_order_iter};
use crate::utilities::order::{ LexicographicOrderDominatedByReverselength, OrderOperatorAuto, OrderOperatorAutoReverse, OrderOperatorByKey, OrderOperatorByKeyCustom};        








//  ===================================================================================
//  DOWKER COMPLEX
//  ===================================================================================


/// Represents the Dowker complex of a binary relation
/// 
/// Here with think of a binary relation between sets `S` and `T` as a binary matrix
/// with rows labeled by the elements of `S` and columns labeled by the elements of `T`.
/// 
/// In order to assist with indexing, we constrain the elements of `T` as follows. 
/// Note that **we call the elements of `T` vertices**.
/// - They should implement [Ord].
/// - It should be possible to convert them to `usize`; this requirement is formalized
///   by requiring `usize` to implement `From<Vertex>`.
/// - The order on vertices should be the same as the order obtained by converting
///   vertices to usize.
/// 
/// We represent the relation as a  `Vec< SortedVec< usize > >`. If `r` has this type
/// then we regard `r` as the relation that contains `(s,t)` iff `r[s]` contains `t`.
#[derive(Clone)]
pub struct DowkerComplex
            < Vertex, RingOperator >
    where
        Vertex:             Clone + Ord + Hash,
        RingOperator:       Clone + RingOperations,
        usize:              From< Vertex >,                
{
    relation_rows:          Vec< SortedVec< Vertex > >,
    relation_columns:       Vec< SortedVec< usize > >,    
    ring_operator:          RingOperator,
}


impl < Vertex, RingOperator >
    
    DowkerComplex
        < Vertex, RingOperator >
    where
        Vertex:             Clone + Debug + Ord + Hash,
        RingOperator:       Clone + RingOperations,
        usize:              From< Vertex >,                
{
    /// Construct a new Dowker complex
    pub fn new( relation_rows: Vec< SortedVec< Vertex > >, ring_operator: RingOperator ) -> Self {

        let max_column_index: Option<usize>        =   relation_rows
                                            .iter()
                                            .map(|simplex| simplex.vec().iter().cloned().map(|x| x.into())  )
                                            .flatten()
                                            .max();

        // if no dowker simplex is nonempty then we don't need to do any more work
        if max_column_index.is_none() {
            return DowkerComplex{ relation_rows, relation_columns: Vec::with_capacity(0), ring_operator, }
        }

        // otherwise let's construct the dual dowker sets
        // first we will count how much space to allocate to each row
        let mut row_counts          =   vec![ 0usize; max_column_index.unwrap() + 1 ];
        for dowker_set in relation_rows.iter() {
            for vertex in dowker_set.vec() {
                let vertex_usize: usize     =   vertex.clone().into();
                row_counts[ vertex_usize ] += 1; // the `into` method convertex the vertex to usize
            }
        }
        // now allocate the dual
        let mut transpose: Vec<Vec<usize>>    =   Vec::with_capacity( max_column_index.unwrap() );
        for row_count in row_counts { transpose.push( Vec::with_capacity(row_count) ) }

        // now construct the dual
        for ( row_index, dowker_set ) in relation_rows.iter().enumerate() {
            for vertex in dowker_set.vec() {
                let vertex_usize: usize     =   vertex.clone().into();
                let row_of_transpose: &mut Vec<usize> = transpose.get_mut( vertex_usize ).unwrap();
                row_of_transpose.push( row_index ); // the `into` method convertex the vertex to usize
            }
        }      
        let mut transpose_safe =   Vec::with_capacity( transpose.len() );
        for sorted_set in transpose.into_iter() {
            transpose_safe.push( SortedVec::new(sorted_set).ok().unwrap() );
        }

        DowkerComplex{ relation_rows, relation_columns: transpose_safe, ring_operator, }
    }

    /// Construct a new Dowker complex from vectors; panics if one of the vectors is not sorted in strictly ascending order.
    pub fn from_vectors( dowker_simplices: Vec< Vec< Vertex > >, ring_operator: RingOperator ) -> Result< Self, Vec< Vertex > > {
        let mut relation_rows = Vec::with_capacity( dowker_simplices.len() );
        for ( counter, simplex ) in dowker_simplices.into_iter().enumerate() {
            match SortedVec::new(simplex) {
                Err( vec ) => { 
                    println!("Error: attempted to create a Dowker boundary matrix from a sequence of simplices, but simplex number {:?}, which is {:?}, is not sorted.", counter, &vec); 
                    return Err( vec );
                } Ok(set) => { 
                    relation_rows.push( set ); 
                }
            }            
        }
        // let relation_rows = dowker_simplices.into_iter().map(|sequence| SortedVec::new(sequence)? ).collect();
        Ok( DowkerComplex::new( relation_rows, ring_operator ) )
    }    
    
    /// The rows of the relation
    /// 
    /// Each row is formatted as a [SortedVec] of elements (or vertices).
    pub fn relation_rows( &self ) -> & Vec< SortedVec< Vertex > > { & self.relation_rows }  

    /// The columns of the relation
    /// 
    /// Each column is formatted as a [SortedVec] of `usize` integers
    pub fn relation_columns( &self ) -> & Vec< SortedVec< usize > > { & self.relation_columns }      


    /// Returns the maximum index of any vertex (convertex to `usize`)
    fn max_vertex( &self ) -> Option< Vertex > {
        self.relation_rows
            .iter()
            .map(|simplex| simplex.vec().iter().cloned()  )
            .flatten()
            .max()
    }    

    /// Returns the sequence of simplices
    /// contained in the Dowker complex that have dimension ≤ `max_simplex_dimension`, ordered (first) 
    /// in ascending order of dimension, and (second) in descending lexicographic
    /// order (excluding simplices of dimension `> max_simplex_dimension`).
    /// 
    /// This is the same order in which we visit rows of the bounday matrix during the cohomology reduction
    /// algorithm.
    /// 
    /// 
    /// **There are many other tools to enumerate the simplices in the Dowker complex, in addition.** See the
    /// methods in [crate::topology::simplicial::simplices::vector] for examples.
    pub fn simplices_in_row_reduction_order( 
                    &self, 
                    max_simplex_dimension: isize 
                ) -> 
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
    {
        return dimension_0_through_d_simplices_in_ascending_dimension_descending_lexicographic_order_iter( 
                    self.relation_rows(), 
                    max_simplex_dimension 
                )

    } 

    /// Returns the simplices of the Dowker comples (up to dimension `max_simplex_dimension`) in lexicographic order
    /// 
    /// **There are many other tools to enumerate the simplices in the Dowker complex, in addition.** See the
    /// methods in [crate::topology::simplicial::simplices::vector] for examples. Most of these methods take a
    /// `& Vec< SortedVec< Vertex > >` as input; you can obtain this from a [DowkerComplex] using the [DowkerComplex::relation_rows]
    /// method.
    pub fn simplices_in_lexicographic_order( &self, max_simplex_dimension: isize )
        ->
        Flatten<
            std::vec::IntoIter<
                TwoTypeIterator<
                    std::iter::Empty< Vec<Vertex> >,
                    Dedup< KMerge<  Combinations<Cloned<Iter<Vertex>>> > >,
                >
            >
        >        
    {
        dimension_0_through_d_simplices_in_dimensionwise_lexicographic_order_iter( & self.relation_rows, max_simplex_dimension )
    }


}






impl < Vertex, RingOperator >

    MatrixOracle for 

    DowkerComplex
        < Vertex, RingOperator >

    where
        Vertex:                     Clone + Debug + Ord + Hash,
        RingOperator:               Clone + RingOperations,
        usize:                      From< Vertex >,
{
    type Coefficient            =   RingOperator::Element;

    type RowIndex               =   Vec< Vertex >;

    type ColumnIndex            =   Vec< Vertex >;

    type RowEntry               =   ( Vec< Vertex >, Self::Coefficient );

    type ColumnEntry            =   ( Vec< Vertex >, Self::Coefficient );

    type Row                    =   DowkerBoundaryMatrixRow< Vertex, RingOperator >;

    type RowReverse             =   DowkerBoundaryMatrixRowReverse< Vertex, RingOperator >;

    type Column                 =   SimplexBoundaryAscend< Vertex, RingOperator >;

    type ColumnReverse          =   SimplexBoundaryDescend< Vertex, RingOperator >;

    fn row(                     &   self, index: & Self::RowIndex   )   -> Self::Row {
        DowkerBoundaryMatrixRow::from_vec_of_dowker_sets( index.clone(), & self.relation_rows, self.ring_operator.clone() ).unwrap()
    }

    fn row_reverse(             &   self, index: & Self::RowIndex   )   -> Self::RowReverse {
        DowkerBoundaryMatrixRowReverse::from_vec_of_dowker_sets( index.clone(), & self.relation_rows, self.ring_operator.clone() ).unwrap()
    }

    fn column(                  &   self, index: & Self::ColumnIndex)   -> Self::Column {
        SimplexBoundaryAscend::new( index.clone(), self.ring_operator.clone() )
    }

    fn column_reverse(          &   self, index: & Self::ColumnIndex)   -> Self::ColumnReverse {
        SimplexBoundaryDescend::new( index.clone(), self.ring_operator.clone() )
    }

    fn has_row_for_index(     &   self, index: & Self::RowIndex   )   -> bool {

        // the vector must be strictly sorted
        if ! index.iter().is_sorted_by( |a,b| a.lt(b) ) {
            return false
        }

        // we have an `index` which is a collection of integers indexing the columns of the binary relation
        // matrix. the question is whether there exists a row index where all of these rows are
        // nonzero. we can answer this question by intersecting the column indices of all the rows
        IntersectOrderedIterators::new(
            index.iter().map(|v| self.relation_columns[ usize::from(v.clone()) ].vec().iter().cloned() )               
        )
        .next()
        .is_some()
    }

    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool {
        // the test for containing a row is the same as the test for containing a column
        self.has_row_for_index(index)
    }

    fn structural_nonzero_entry(&   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient > {

        // both indices must be
        if !  self.has_row_for_index(row)  {
            panic!("Attempted to look up the entry in row {:?}, column {:?} of a Dowker boundary matrix, but {:?} is either an invalid simplex (not sorted or has a repeat entry) or it is not contained in the Dowker compelx", row, column, row );
        }
        if !  self.has_column_for_index(row)  {
            panic!("Attempted to look up the entry in row {:?}, column {:?} of a Dowker boundary matrix, but {:?} is either an invalid simplex (not sorted or has a repeat entry) or it is not contained in the Dowker compelx", row, column, column );
        }        

        // this criterion has to be satisfied to have a facet-cofacet pair
        if column.len() != row.len() + 1 { 
            return None 
        }

        // get the symmetric difference
        let mut symmetric_difference                        =   symmetric_difference_of_ordered_iterators( row.iter(), column.iter() );
        let first_different_vertex                      =   symmetric_difference.next();  
        
        if symmetric_difference.next().is_some() {
            // in this case the row and column indices differ by more than one vertex, so they cannot be a facet-cofacet pair
            return None
        } else {
            // in this case we do have a facet-cofacet pair, so we have to calculate the corresponding coefficient by finding 
            // which vertex in column is not contained in row, and then calculating the corresponding coefficient
            for (counter, vertex) in column.iter().enumerate() {
                if Some(vertex) == first_different_vertex {
                    let coefficient         =   self.ring_operator.minus_one_to_power( counter );
                    return Some( coefficient )
                }
            }
        }

        panic!("Error retreiving structural nonzero entry.");
    }
}









impl < Vertex, RingOperator >

    MatrixAlgebra for 

    DowkerComplex
        < Vertex, RingOperator >

    where
        Vertex:                     Clone + Debug + Ord + Hash,
        RingOperator:               Clone + RingOperations,
        usize:                      From< Vertex >,
{
    type RingOperator                                   =   RingOperator;

    type OrderOperatorForRowEntries                     =   OrderOperatorByKeyCustom<LexicographicOrderDominatedByReverselength>;

    type OrderOperatorForRowIndices                     =   LexicographicOrderDominatedByReverselength;

    type OrderOperatorForColumnEntries                  =   OrderOperatorByKeyCustom<LexicographicOrderDominatedByReverselength>;

    type OrderOperatorForColumnIndices                  =   LexicographicOrderDominatedByReverselength;

    fn ring_operator( &self ) -> Self::RingOperator {
        self.ring_operator.clone()
    }

    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries {
        OrderOperatorByKeyCustom::new(LexicographicOrderDominatedByReverselength::new())
    }

    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices {
        LexicographicOrderDominatedByReverselength::new()
    }

    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries {
        OrderOperatorByKeyCustom::new(LexicographicOrderDominatedByReverselength::new())
    }

    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices {
        LexicographicOrderDominatedByReverselength::new()
    }
}










impl < Vertex, RingOperator >

    MatrixOracleOperations for 

    DowkerComplex
        < Vertex, RingOperator >

    where // these are the requirements to implement the `MatrixOracle` trait
        Vertex:                     Clone + Debug + Ord + Hash,
        RingOperator:               Clone + RingOperations,
        usize:                      From< Vertex >,

{}





















impl < Vertex, RingOperator >

    ChainComplex for 

    DowkerComplex
        < Vertex, RingOperator >

    where
        Vertex:                     Clone + Debug + Ord + Hash,
        RingOperator:               Clone + RingOperations,
        usize:                      From< Vertex >,
{
    type BasisVectorIndicesIterable         =   Vec< Vec< Vertex > >;

    /// Returns simplices of dimension `dimension` in lexicographic order.
    fn basis_vector_indices_for_dimension( &self, dimension: isize ) -> Self::BasisVectorIndicesIterable {
        dimension_d_simplices_in_lexicographic_order_iter(
            & self.relation_rows, 
            dimension
        ).collect()
    }

    fn dimension_for_basis_vector_with_index( &self, index: & Self::RowIndex ) -> Result<isize, Self::RowIndex> {
        if ! self.has_row_for_index(index) {
            return Err( index.clone() )
        }
        Ok( (index.len() as isize) - 1 )
    }
}
















//  ===================================================================================
//  ROWS AND COLUMNS
//  ===================================================================================


//  COUBOUNDARY VECTOR DESCENDING
//  -----------------------------------------------------------------------------------


/// Iterates over the terms of the coboundary of a simplex in descending lexicographic order.
/// 
/// The coefficient is calculated as `(-1)^k`, where `xk` is the vertex which appears twice in the template.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::topology::simplicial::from::relation::DowkerBoundaryMatrixRowReverse;
/// use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;
/// use oat_rust::utilities::sequences_and_ordinals::SortedVec;
/// 
/// let ring_operator = PrimeOrderField::new(3);
/// let simplex = vec![1,3];
/// let relation_rows = vec![ SortedVec::from_iter( vec![0,1,2,3,4] ).unwrap() ];
/// // note we have to unwrap, because the constructor returns a Result
/// let coboundary = DowkerBoundaryMatrixRowReverse::from_vec_of_dowker_sets( simplex, &relation_rows, ring_operator ).unwrap();
/// 
/// itertools::assert_equal( coboundary, vec![ (vec![1,3,4], 1), (vec![1,2,3], 2), (vec![0,1,3], 1) ]);
/// ```
#[derive(Clone, Debug)]
pub struct DowkerBoundaryMatrixRowReverse< Vertex, RingOperator >
    where 
        RingOperator:       RingOperations,
        Vertex:             Ord + Clone,
{
    next_cofacet_opt:           Option< Vec< Vertex > >,
    next_coefficient:           RingOperator::Element,
    vertices_to_insert:         Vec< Vertex >,              // should be sorted in ascending order
    retrieval_locus:            usize,                      // the place where the vertex inserted into `next_cofacet` was drawn from 
    insertion_locus:            usize,                      // this points to the place in `next_cofacet` where we inserted a vertex    
    ring_operator:              RingOperator,
}   


impl < Vertex, RingOperator >

DowkerBoundaryMatrixRowReverse
        < Vertex, RingOperator >

    where 
        RingOperator:       RingOperations,
        Vertex:             Clone + Debug + Hash + Ord,
{
    /// Generates a [DowkerBoundaryMatrixRow] for a simplex, given a CSR representation of the binary dowker matrix.    
    pub fn from_vec_of_dowker_sets( facet: Vec< Vertex >, relation_rows: &Vec< SortedVec< Vertex > >, ring_operator: RingOperator ) -> Result< Self, Vec< Vertex > > {

        //  the coboundary of the empty simplex is zero;
        //  ---------------------------------------------
        if facet.is_empty() {
            // this iterator will be identified as empty, because `vertices_to_insert` is empty
            return Ok( DowkerBoundaryMatrixRowReverse{
                next_cofacet_opt:           None,
                next_coefficient:           RingOperator::one(),   // these are arbitrary values            
                vertices_to_insert:         vec![],                // these are arbitrary values  
                retrieval_locus:            0,                     // these are arbitrary values  
                insertion_locus:            0,                     // these are arbitrary values  
                ring_operator,         // these are arbitrary values                    
            } )                   
        }

        let facet = SortedVec::new( facet ); // place the facet in a safe wrapper; this will panic if the facet isn't strictly sorted
        if let Err(vec) = facet { 
            println!("Error: attempted to compute the coboundary of simplex {:?}, but the vertices of the simplex are not sorted in strictly ascending order", &vec );
            return Err( vec );
        }
        let facet = facet.unwrap(); // we have handled the error, so it's safe to unwrap
        let vertices_to_insert: Vec< Vertex >   =   relation_rows
                                                        .iter()
                                                        .filter( |x|  x.contains_subset( &facet ) )
                                                        .map(|x| x.vec().iter())
                                                        .kmerge()
                                                        .dedup()
                                                        .filter(|x| ! facet.contains(x) )
                                                        .cloned()
                                                        .collect();


        // //  obtain a list of vertices to insert
        // //  ---------------------------------------------        
        // let facet_vertex_set = SortedVec::from_iter( facet.iter().cloned() );
        // let mut vertices_to_insert = SortedVec::new();

        // // take the union of all dowker super-simplices
        // for dowker_set in relation_rows
        //                                 .iter() // iterate over dowker siplices
        //                                 .filter( |x| x.is_superset( &facet_vertex_set) ) // only keep those that contain the facet
        //     {
        //         vertices_to_insert.extend( dowker_set.iter().cloned() )
        //     }

        // // take a set difference: vertices_to_insert \ facet
        // for vertex in facet.iter() { vertices_to_insert.remove( vertex ); }

        // // convert to a vector
        // let mut vertices_to_insert: Vec< Vertex > = vertices_to_insert.into_iter().collect();

        // // sort the vector
        // vertices_to_insert.sort();
        

        //  if there are no vertices to insert, then the coboundary is zero
        //  ---------------------------------------------
        if vertices_to_insert.is_empty() {
            // this iterator is empty, because `next_cofacet_opt` is None
            return Ok( DowkerBoundaryMatrixRowReverse{
                next_cofacet_opt:           None,
                next_coefficient:           RingOperator::one(),     // these are abitrary values            
                vertices_to_insert:         vec![],                  // these are abitrary values  
                retrieval_locus:            0,                       // these are abitrary values  
                insertion_locus:            0,                       // these are abitrary values  
                ring_operator,           // these are abitrary values                    
            } )                  
        }


        //  generate the rest of the initial data
        //  ---------------------------------------------                

        let mut coefficient = ring_operator.minus_one_to_power( facet.len() );
        let mut next_cofacet = facet.into_vec();
        let mut insertion_locus = next_cofacet.len();
        let retrieval_locus = vertices_to_insert.len() - 1;        
        let inserted_vertex = vertices_to_insert[retrieval_locus].clone();
        while   ( insertion_locus > 0 ) 
                && 
                ( inserted_vertex < next_cofacet[ insertion_locus - 1 ] ) {
            insertion_locus -= 1;
            coefficient = ring_operator.negate( coefficient );
        }
        next_cofacet.insert( insertion_locus, inserted_vertex );

        Ok( DowkerBoundaryMatrixRowReverse{
            next_cofacet_opt:           Some( next_cofacet ),
            next_coefficient:           coefficient,
            vertices_to_insert,         // should be sorted in ascending order
            retrieval_locus,            // the place where the inserted vertex was drawn from 
            insertion_locus,            // this points to the place in `next_cofacet` where we inserted a vertex
            ring_operator,             
        } )

    }

    // /// Generates a [DowkerBoundaryMatrixRow] for a simplex, given CSR and CSC representations of the binary dowker matrix.
    // fn from_csr_and_csc_matrices( facet: Vec< Vertex >, dowker_matrix_csr: Vec< SortedVec< Vertex > >, dowker_matrix_csc: Vec< SortedVec< Vertex > > ) {        
    // }
}


impl < Vertex, RingOperator >

    Iterator for

    DowkerBoundaryMatrixRowReverse
        < Vertex, RingOperator >

    where 
        RingOperator:       RingOperations,
        Vertex:             Ord + Clone,   
{
    type Item = (Vec<Vertex>, RingOperator::Element);

    fn next( &mut self ) -> Option< Self::Item >{

        match self.next_cofacet_opt {
            None => { None }
            Some( ref mut next_cofacet ) => {

                // make a copy of the return value, before we change the internal state of the iterator
                let return_value    =   ( next_cofacet.clone(), self.next_coefficient.clone() );

                // if there are more vertices to insert, update the internal state with data for the next entry
                if self.retrieval_locus > 0 {
                    
                    // grab the next vertex
                    self.retrieval_locus -= 1; 
                    let inserted_vertex     =   self.vertices_to_insert[ self.retrieval_locus ].clone();
                    
                    // update pointers and coefficients
                    while   ( self.insertion_locus > 0 ) 
                            && 
                            ( inserted_vertex < next_cofacet[ self.insertion_locus - 1 ] ) {
                        next_cofacet[ self.insertion_locus ] = next_cofacet[ self.insertion_locus - 1 ].clone();
                        self.insertion_locus -= 1;
                        self.next_coefficient = self.ring_operator.negate( self.next_coefficient.clone() );
                    }
                    // update the cofacet
                    next_cofacet[ self.insertion_locus ] = inserted_vertex;
                } 
                // otherwise flag the fact that there will be no more entries, after we return the one that we have already extracted
                else {
                    self.next_cofacet_opt = None;
                }

                Some( return_value )
            }
        }

    }
}       



//  COUBOUNDARY VECTOR ASCENDING
//  -----------------------------------------------------------------------------------

/// Iterates over the terms of the coboundary of a simplex in ascending lexicographic order.
/// 
/// The coefficient is calculated as `(-1)^k`, where `xk` is the vertex which appears twice in the template.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::topology::simplicial::from::relation::DowkerBoundaryMatrixRow;
/// use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;
/// use oat_rust::utilities::sequences_and_ordinals::SortedVec;
/// 
/// let ring_operator = PrimeOrderField::new(3);
/// let simplex = vec![1,3];
/// let relation_rows = vec![ SortedVec::from_iter( vec![0,1,2,3,4] ).unwrap() ];
/// 
/// // note we have to unwrap, because the constructor returns a Result
/// let coboundary = DowkerBoundaryMatrixRow::from_vec_of_dowker_sets( simplex, &relation_rows, ring_operator ).unwrap();
/// 
/// itertools::assert_equal( coboundary, vec![ (vec![0,1,3], 1), (vec![1,2,3], 2), (vec![1,3,4], 1) ]);
/// ```
#[derive(Clone, Debug)]
pub struct DowkerBoundaryMatrixRow< Vertex, RingOperator >
    where 
        RingOperator:       RingOperations,
        Vertex:             Ord + Clone,
{
    next_cofacet_opt:           Option< Vec< Vertex > >,
    next_coefficient:           RingOperator::Element,
    vertices_to_insert:         Vec< Vertex >,              // should be sorted in ascending order
    retrieval_locus:            usize,                      // the place where the vertex inserted into `next_cofacet` was drawn from 
    insertion_locus:            usize,                      // this points to the place in `next_cofacet` where we inserted a vertex    
    ring_operator:              RingOperator,
}   


impl < Vertex, RingOperator >

    DowkerBoundaryMatrixRow
        < Vertex, RingOperator >

    where 
        RingOperator:       RingOperations,
        Vertex:             Clone + Debug + Hash + Ord,
{
    /// Generates a [DowkerBoundaryMatrixRow] for a simplex, given a CSR representation of the binary dowker matrix.    
    pub fn from_vec_of_dowker_sets( facet: Vec< Vertex >, relation_rows: &Vec< SortedVec< Vertex > >, ring_operator: RingOperator ) -> Result< Self, Vec< Vertex > > {

        //  the coboundary of the empty simplex is zero;
        //  ---------------------------------------------
        if facet.is_empty() {
            // this iterator will be identified as empty, because `vertices_to_insert` is empty
            return Ok( DowkerBoundaryMatrixRow{
                next_cofacet_opt:           None,
                next_coefficient:           RingOperator::one(),   // these are arbitrary values            
                vertices_to_insert:         vec![],                // these are arbitrary values  
                retrieval_locus:            0,                     // these are arbitrary values  
                insertion_locus:            0,                     // these are arbitrary values  
                ring_operator,         // these are arbitrary values                    
            } )                
        }


        let facet = SortedVec::new( facet ); // place the facet in a safe wrapper; this will panic if the facet isn't strictly sorted
        if let Err(vec) = facet { 
            println!("Error: attempted to compute the coboundary of simplex {:?}, but the vertices of the simplex are not sorted in strictly ascending order", &vec );
            return Err( vec );
        }
        let facet = facet.unwrap(); // we have handled the error, so it's safe to unwrap
        let vertices_to_insert: Vec< Vertex >   =   relation_rows
                                                        .iter()
                                                        .filter( |x|  x.contains_subset( &facet ) )
                                                        .map(|x| x.vec().iter())
                                                        .kmerge()
                                                        .dedup()
                                                        .filter(|x| ! facet.contains(x) )
                                                        .cloned()
                                                        .collect();
        // //  obtain a list of vertices to insert
        // //  ---------------------------------------------        
        // let facet_vertex_set = SortedVec::from_iter( facet.iter().cloned() );
        // let mut vertices_to_insert = SortedVec::new();

        // // take the union of all dowker super-simplices
        // for dowker_simplex in relation_rows
        //                                 .iter() // iterate over dowker siplices
        //                                 .filter( |x| x.is_superset( &facet_vertex_set) ) // only keep those that contain the facet
        //     {
        //         vertices_to_insert.extend( dowker_simplex.iter().cloned() )
        //     }

        // // take a set difference: vertices_to_insert \ facet
        // for vertex in facet.iter() { vertices_to_insert.remove( vertex ); }

        // // convert to a vector
        // let mut vertices_to_insert: Vec< Vertex > = vertices_to_insert.into_iter().collect();

        // // sort the vector
        // vertices_to_insert.sort();
        

        //  the coboundary of the empty simplex is zero;
        //  ---------------------------------------------
        if vertices_to_insert.is_empty() {
            // this iterator is empty, because `vertices_to_insert` is empty
            return Ok( DowkerBoundaryMatrixRow{
                next_cofacet_opt:           None,
                next_coefficient:           RingOperator::one(),
                vertices_to_insert:         vec![],     
                retrieval_locus:            0,          
                insertion_locus:            0,          
                ring_operator,                
            } )                   
        }


        //  generate the rest of the initial data
        //  ---------------------------------------------                

        let mut coefficient = RingOperator::one();
        let mut next_cofacet = facet.into_vec();
        let mut insertion_locus = 0;
        let retrieval_locus = 0;        
        let inserted_vertex = vertices_to_insert[retrieval_locus].clone();
        while   ( insertion_locus < next_cofacet.len() ) 
                && 
                ( inserted_vertex > next_cofacet[ insertion_locus ] ) {
            insertion_locus += 1;
            coefficient = ring_operator.negate( coefficient );
        }
        next_cofacet.insert( insertion_locus, inserted_vertex );

        Ok( DowkerBoundaryMatrixRow{
            next_cofacet_opt:           Some( next_cofacet ),
            next_coefficient:           coefficient,
            vertices_to_insert,         // should be sorted in ascending order
            retrieval_locus,            // the place where the inserted vertex was drawn from 
            insertion_locus,            // this points to the place in `next_cofacet` where we inserted a vertex
            ring_operator,             
        } )

    }

    // /// Generates a [DowkerBoundaryMatrixRow] for a simplex, given CSR and CSC representations of the binary dowker matrix.
    // fn from_csr_and_csc_matrices( facet: Vec< Vertex >, dowker_matrix_csr: Vec< SortedVec< Vertex > >, dowker_matrix_csc: Vec< SortedVec< Vertex > > ) {        
    // }
}


impl < Vertex, RingOperator >

    Iterator for

    DowkerBoundaryMatrixRow
        < Vertex, RingOperator >

    where 
        RingOperator:       RingOperations,
        Vertex:             Ord + Clone,   
{
    type Item = (Vec<Vertex>, RingOperator::Element);

    fn next( &mut self ) -> Option< Self::Item >{

        // println!("{:?} -- DELETE THIS AND ALL DEBUG REQUIREMENTS ON THIS IMPLEMENTATION OF ITER", &self);

        match self.next_cofacet_opt {
            None => { None }
            Some( ref mut next_cofacet ) => {

                // make a copy of the return value, before we change the internal state of the iterator
                let return_value    =   ( next_cofacet.clone(), self.next_coefficient.clone() );

                // if there are more vertices to insert, update the internal state with data for the next entry
                if self.retrieval_locus + 1 < self.vertices_to_insert.len() {
                    
                    // grab the next vertex
                    self.retrieval_locus += 1; 
                    let inserted_vertex     =   self.vertices_to_insert[ self.retrieval_locus ].clone();
                    
                    // update pointers and coefficients
                    while   ( self.insertion_locus + 1 < next_cofacet.len() ) 
                            && 
                            ( inserted_vertex > next_cofacet[ self.insertion_locus + 1 ] ) {
                        next_cofacet[ self.insertion_locus ] = next_cofacet[ self.insertion_locus +1 ].clone();
                        self.insertion_locus += 1;
                        self.next_coefficient = self.ring_operator.negate( self.next_coefficient.clone() );
                    }
                    // update the cofacet
                    next_cofacet[ self.insertion_locus ] = inserted_vertex;
                } 
                // otherwise flag the fact that there will be no more entries, after we return the one that we have already extracted
                else {
                    self.next_cofacet_opt = None;
                }

                Some( return_value )
            }
        }

    }
}    















//  ===================================================================================
//  CONSTRUCTORS
//  ===================================================================================



/// Construct a "sideways ladder".  
/// 
/// # Examples
/// 
/// The following is created by calling
/// `sideways_ladder_edges( 0, 2 )`
/// 
/// ```text
/// 1 ---- 3 ---- 5
/// |      |      |
/// |      |      |
/// 0 ---- 2 ---- 4
/// ```
/// 
/// The following is created by calling
/// `sideways_ladder( 1,3 )`
/// 
/// ```text
/// 3 ---- 5 ---- 7 ---- 9
/// |      |      |      |
/// |      |      |      |
/// 2 ---- 4 ---- 6 ---- 8
/// ```
/// # Code example
/// 
/// ```
/// use oat_rust::topology::simplicial::from::relation::sideways_ladder_edges;
/// 
/// // construct the edges
/// let number_of_holes         =   1;
/// let offset_from_left        =   1;
/// let mut edges           =   sideways_ladder_edges(number_of_holes, offset_from_left);
/// edges.sort();
/// 
/// // this is the ground truth
/// let ground_truth            =   vec![ 
///                                     vec![2,3], 
///                                     vec![2,4],
///                                     vec![3,5],
///                                     vec![4,5],                                            
///                                 ];
/// 
/// assert_eq!( edges, ground_truth );
/// ```
pub fn sideways_ladder_edges( offset_from_left: usize, number_of_holes: usize ) -> Vec< Vec< usize > > {
    let mut hyperedges                              =   Vec::with_capacity( 1 + 3 * number_of_holes );

    // vertical parts
    for p in offset_from_left .. offset_from_left + number_of_holes + 1 {
        hyperedges.push( vec![2*p, 2*p + 1] )
    }
        
    // horizontal parts
    for p in offset_from_left .. offset_from_left + number_of_holes {
        hyperedges.push( vec![ 2*p,     2*(p+1)     ] );
        hyperedges.push( vec![ 2*p + 1, 2*(p+1) + 1 ] );
    }

    return hyperedges
}






















//  ===================================================================================
//  TESTING
//  ===================================================================================




use crate::algebra::matrices::debug::{matrix_oracle_is_internally_consistent, matrix_order_operators_are_internally_consistent };
use crate::algebra::rings::types::field_prime_order::PrimeOrderField;


/// Builds a Dowker boundary matrix for the user-provided list of relation rows, and validates it with [matrix_oracle_is_internally_consistent] and [matrix_order_operators_are_internally_consistent]
/// 
/// The `relation_rows` argument is a list of sorted lists, where each sorted list records the column indices of the nonzero entries in a given row of the 
/// binary relation matrix.
pub fn validate_dowker_boundary_matrix<Vertex>( 
        relation_rows: Vec<SortedVec<Vertex>>, 
        max_dim: isize, 
    )  

    where 
        Vertex:     Clone + Hash + Debug + Ord,
        usize:      From< Vertex >,
{
    
    // define the boundary matrix
    // --------------------------

    let ring_operator                               =   PrimeOrderField::new(47);    
    let boundary_matrix                             =   DowkerComplex::new( relation_rows.clone(), ring_operator );        

    // verify that the matrix lookup operations are internally consistent
    // ------------------------------------------------------------------

    // get row indices in order
    let mut sorted_row_indices: Vec<_>               =   boundary_matrix.simplices_in_row_reduction_order( max_dim ).collect();
    ( &mut sorted_row_indices ).reverse(); // row reduction order is the REVERSE of the actual order on rows

    // get column indices in order    
    let mut sorted_column_indices: Vec<_>               =   boundary_matrix.simplices_in_row_reduction_order( max_dim + 1 ).collect();
    ( &mut sorted_column_indices ).reverse(); // row reduction order is the REVERSE of the actual order on rows    


    assert!(
        matrix_oracle_is_internally_consistent(
            & boundary_matrix,
            sorted_row_indices.clone(),
            sorted_column_indices.clone(),
        )
    );

    // verify that the matrix order operators are internally consistent
    // ----------------------------------------------------------------
    assert!(
        matrix_order_operators_are_internally_consistent(
            & boundary_matrix,
            sorted_row_indices.clone(),
            sorted_column_indices.clone(),
        )
        .is_ok()
    ); 
}




#[cfg(test)]
mod tests {
    

use crate::utilities::random::random_sequences;

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_dowker_boundary_small(){
        
        // used for enumerating simplices
        let relation_rows: Vec<Vec<usize>> =  
                vec![ 
                        vec![0,1,2],
                        vec![0,3],                      
                        vec![1,3], 
                        vec![2,3]                                           
                    ];     
        let relation_rows =     relation_rows.into_iter()
                                    .map(|v| SortedVec::new(v).ok().unwrap() )
                                    .collect::<Vec<_>>();

        validate_dowker_boundary_matrix( relation_rows, 2 );
    }



    #[test]
    fn test_dowker_boundary_big(){
        
        for _ in 0..10 {
            
            let relation_rows    =   random_sequences(7, (0..7).step_by(2), 0.1 );
            let _verbose = false;
            validate_dowker_boundary_matrix( relation_rows, 2 );
        }

    }    



}   



#[cfg(test)]
mod docstring_tests {
    use crate::topology::simplicial::simplices::vector::{dimension_d_simplices_in_lexicographic_order_iter, dimension_d_simplices_in_reverse_lexicographic_order_iter, dimension_0_through_d_simplices_in_ascending_dimension_descending_lexicographic_order_iter};

    



    #[test]
    fn docstring_test_dowker_homology() {
        use itertools::Itertools;
            
        use crate::topology::simplicial::from::relation::DowkerComplex;
            
        use crate::algebra::matrices::operations::umatch::differential::DifferentialUmatch;     
        use crate::algebra::rings::types::field_prime_order::PrimeOrderField;        
      
        use crate::utilities::sequences_and_ordinals::SortedVec;
            
        // Parameters
        // ----------
            
        // Define the maximum homology dimensino we want to compute.
        let min_homology_dimension                         =   0;
        let max_homology_dimension                  =   2;
            
        // Define the ring operator for the finite field of order 3.
        // You can use this object to perform arithmetic operations, e.g., 
        // to add 1 and 1, you can run `ring_operator.add(1,1)`.
        let ring_operator          =   PrimeOrderField::new(3);
            
        // We will build a dowker complex.
        // A dowker complex is defined by a vertex set V and a family S
        // of subsets of V.  A subset of V forms a simplex iff it is 
        // a subset of some element of S.  We refer to the elements 
        // of S as "dowker simplices".
            
        // Each dowker simplex is represented by a SortedVec of vertices.
        // We store the list of all such simplices inside a larger vector.
        let dowker_simplices: Vec< SortedVec< usize > > 
            =   vec![    
                        vec![0,1,2], 
                        vec![0,3], 
                        vec![1,3], 
                        vec![2,3]  
                    ]
                    .into_iter()
                    .map( |x| SortedVec::new(x).unwrap() ) 
                    .collect_vec();
            
        //  Boundary matrix
        //  ---------------
            
        // This is a lazy object that generates rows/columns of the boundary matrix, on demand.
        let boundary_matrix = DowkerComplex::new( dowker_simplices.clone(), ring_operator );
            
        //  Simplex iterators
        //  -----------------
            
        // An iterator that runs over all triangles in the complex, in ascending 
        // lexicographic order
        let _triangles_ascending_order = dimension_d_simplices_in_lexicographic_order_iter( &dowker_simplices, 2);
            
        // An iterator that runs over all edges in the complex, in descending 
        // lexicographic order
        let _triangles_descending_order = dimension_d_simplices_in_reverse_lexicographic_order_iter( &dowker_simplices, 2);   
            
        // An iterator that runs over simplices of dimension 0 through max_homology_dimension,
        // ordered first by dimension (ascending) and second by lexicographic order (descending)
        let _row_indices = boundary_matrix.simplices_in_row_reduction_order( max_homology_dimension );

        // row_indices contains a reference to boundary_matrix, which will cause problems. Instead,
        // we can construct an iterator that runs over the same sequence of simplices, like so:
        let row_indices     =   dimension_0_through_d_simplices_in_ascending_dimension_descending_lexicographic_order_iter(&dowker_simplices, max_homology_dimension);
            
        //  Homology computation (by matrix factorization)
        //  ----------------------------------------------
            
        // Factor the boundary matrix
        let factored    =   DifferentialUmatch::new( 
                                boundary_matrix, 
                                min_homology_dimension,
                                max_homology_dimension,
                            );
                        
        //  Printing results
        //  ----------------
                        
        // Betti numbers.  For this computation we have to provide a
        // function that assigns a dimension to each index (i.e. to each simplex)
        let homology_dimensions   =   factored.betti_numbers(); 
        for dim in 0 .. 2 {
            println!(
                    // we'll insert two values into this string
                    "The betti number in dimension {:?} is {:?}.",
                    dim,                  // the dimension
                    homology_dimensions         // and the betti number
                        .get( & dim )     // this looks up a value in the hashmap
                        .unwrap_or( & 0)  // if the hashmap doesn't have the value, then use 0
                );            
        }
        println!(); // insert line break
        
        // Cycle representatives for homology
        println!(
                "The following are basis vectors for homology in dimensions 0 through {:?}",
                max_homology_dimension,
            );
        for (cycle_number, cycle) in factored.homology_basis().enumerate() {
            // `cycle` is an iterator.  For convenience, collect the elements of the
            // iterator into a Rust vector.
            let cycle: Vec<_> = cycle.collect();
            println!("Basis vector {:?} = {:?}", cycle_number, cycle);
        }
    }  


    #[test]
    fn test() {
        use crate::topology::simplicial::from::relation::DowkerBoundaryMatrixRow;
        use crate::algebra::rings::types::field_prime_order::PrimeOrderField;
        use crate::utilities::sequences_and_ordinals::SortedVec;

        let ring_operator = PrimeOrderField::new(3);
        let simplex = vec![1,3];
        let relation_rows = vec![ SortedVec::from_iter( vec![0,1,2,3,4] ).unwrap() ];
        let coboundary = DowkerBoundaryMatrixRow::from_vec_of_dowker_sets( simplex, &relation_rows, ring_operator ).unwrap();

        itertools::assert_equal( coboundary, vec![ (vec![0,1,3], 1), (vec![1,2,3], 2), (vec![1,3,4], 1) ]);        
    }    


    #[test]
    fn doc_test_sideways_ladder() {
        use crate::topology::simplicial::from::relation::sideways_ladder_edges;

        // construct the edges
        let number_of_holes         =   1;
        let offset_from_left        =   1;
        let mut edges           =   sideways_ladder_edges(number_of_holes, offset_from_left);
        edges.sort();

        // this is the ground truth
        let ground_truth            =   vec![ 
                                            vec![2,3], 
                                            vec![2,4], 
                                            vec![3,5], 
                                            vec![4,5],
                                        ];

        assert_eq!( edges, ground_truth );
    }   




    #[test]
    fn doc_test_misc_REDUNDANT_OK_TO_DELETE() {
        use crate::algebra::matrices::operations::umatch::differential::DifferentialUmatch;   
        use crate::algebra::rings::types::field_prime_order::PrimeOrderField;
        use crate::topology::simplicial::from::relation::DowkerComplex;
        use crate::topology::simplicial::simplices::vector::{dimension_d_simplices_in_lexicographic_order_iter, dimension_d_simplices_in_reverse_lexicographic_order_iter}; 
        use crate::topology::simplicial::simplices::vector::dimension_0_through_d_simplices_in_ascending_dimension_descending_lexicographic_order_iter;        
        use crate::utilities::sequences_and_ordinals::SortedVec;

        use itertools::Itertools;

        // Parameters
        // ----------

        // Define the maximum homology dimensino we want to compute.
        let min_homology_dimension                  =   0;
        let max_homology_dimension                  =   2;        

        // Define the ring operator for the finite field of order 3.
        // You can use this object to perform arithmetic operations, e.g., 
        // to add 1 and 1, you can run `ring_operator.add(1,1)`.
        let ring_operator          =   PrimeOrderField::new(3);

        // We will build a dowker complex.
        // A dowker complex is defined by a vertex set V and a family S
        // of subsets of V.  A subset of V forms a simplex iff it is 
        // a subset of some element of S.  We refer to the elements 
        // of S as "dowker simplices".

        // Each dowker simplex is represented by a SortedVec of vertices.
        // We store the list of all such simplices inside a larger vector.
        let dowker_simplices: Vec<SortedVec<usize>>
            =   vec![    
                        vec![0,1,2], 
                        vec![0,3], 
                        vec![1,3], 
                        vec![2,3]  
                    ]
                    .into_iter()
                    .map( |x| SortedVec::new(x).unwrap() )  // we unwrap because `new` can return an error
                    .collect_vec();

        //  Boundary matrix
        //  ---------------

        // This is a lazy object that generates rows/columns of the boundary matrix, on demand.
        let boundary_matrix = DowkerComplex::new( dowker_simplices.clone(), ring_operator.clone() );

        //  Simplex iterators
        //  -----------------

        // An iterator that runs over all triangles in the complex, in ascending 
        // lexicographic order
        let triangles_ascending_order = dimension_d_simplices_in_lexicographic_order_iter( &dowker_simplices, 2);

        // An iterator that runs over all edges in the complex, in descending 
        // lexicographic order
        let triangles_descending_order = dimension_d_simplices_in_reverse_lexicographic_order_iter( &dowker_simplices, 2);   

        // An iterator that runs over simplices of dimension 0 through max_homology_dimension,
        // ordered first by dimension (ascending) and second by lexicographic order (descending)
        let row_indices = boundary_matrix.simplices_in_row_reduction_order( max_homology_dimension );

        // row_indices contains a reference to boundary_matrix, which will cause problems. Instead,
        // we can construct an iterator that runs over the same sequence of simplices, like so:
        let row_indices     =   dimension_0_through_d_simplices_in_ascending_dimension_descending_lexicographic_order_iter(&dowker_simplices, max_homology_dimension);

        //  Homology computation (by matrix factorization)
        //  ----------------------------------------------

        // Factor the boundary matrix
        let factored    =   DifferentialUmatch::new( 
                                boundary_matrix, 
                                min_homology_dimension,
                                max_homology_dimension,
                            );
                        
        //  Printing results
        //  ----------------
                        
        // Betti numbers.  For this computation we have to provide a
        // function that assigns a dimension to each index (i.e. to each simplex)
        let homology_dimensions   =   factored.betti_numbers(); 
        for dim in 0 .. 2 {
            println!(
                    // we'll insert two values into this string
                    "The betti number in dimension {:?} is {:?}.",
                    dim,                  // the dimension
                    homology_dimensions         // and the betti number
                        .get( & dim )     // this looks up a value in the hashmap
                        .unwrap_or( & 0)  // if the hashmap doesn't have the value, then use 0
                );            
        }
        println!(""); // insert line break

        // Cycle representatives for homology
        println!(
                "The following are basis vectors for homology in dimensions 0 through {:?}",
                max_homology_dimension,
            );
        for (cycle_number, cycle) in factored.homology_basis().enumerate() {
            // `cycle` is an iterator.  For convenience, collect the elements of the
            // iterator into a Rust vector.
            let cycle: Vec<_> = cycle.collect();
            println!("Basis vector {:?} = {:?}", cycle_number, cycle);
        }
    }     

}





