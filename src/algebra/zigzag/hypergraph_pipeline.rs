use std::collections::HashSet;
use std::time::Instant;


use itertools::Itertools;
use num::Integer;
use serde::{Deserialize, Serialize};


use crate::algebra::zigzag::{cospan_pipeline, span_pipeline};
use crate::algebra::zigzag::spans::induced_span;
use crate::algebra::rings::traits::DivisionRingOperations; 
use crate::algebra::rings::types::field_prime_order::BooleanField;
use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
use crate::utilities::order::OrderOperatorByKey;



use super::decompose::Diagonalization;
use super::{cospans::{factor_dowker_complex, induced_cospan}, decompose::{QuiverReprsentation, SingleBarBasisVectorIndexLedger}};


use derive_getters::{Getters, Dissolve};
use derive_new::new;

// THIS JUST MAKES THE CODE EASIER TO READ
type Simplex        =   Vec<usize>;
type Chain<RingElement>          =   Vec< (Simplex, RingElement) >;



/// The dimension of a linear combination of simplices
/// 
/// Assumes that the linear combination is homogeous with respect to simplex dimension
pub fn chain_dimension< RingElement >( chain: & Vec< ( Simplex, RingElement ) > ) -> usize {
    chain[0].0.len() - 1
}









#[derive(new,Getters,Clone,Debug,Dissolve,Eq,PartialEq,Ord,PartialOrd,Serialize,Deserialize)]
pub struct IntervalSubmoduleBasis< RingElement > {
    pub leftmost_vertex:                                    usize,
    pub basis_vectors:                                      Vec< Chain< RingElement > >,
}



impl < RingElement >

    IntervalSubmoduleBasis< RingElement >
{

    pub fn interval_endpoints( &self ) -> (usize, usize) {
        let left    =   self.leftmost_vertex.clone();
        let right   =   left + self.basis_vectors.len();
        return ( left, right )
    }
}    








/// Returns an interval decomposition of the zigzag module of a sequence of hypergraphs connected by their unions.
/// 
/// Concretely, we take a sequence of hypergraphs `A, B, C, ..` and decompose the zigzag module for 
/// `A --> A u B <-- B --> B u C <-- ..`.
/// 
/// The coefficient field is determined by the choice of `ring_operator`.
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::rings::types::field_prime_order::BooleanField;
/// use oat_rust::topology::simplicial::from::relation::sideways_ladder_edges;
/// use oat_rust::algebra::zigzag::hypergraph_pipeline::interval_decomposition_for_zigzag_of_hypgeraph_unions;
/// 
/// 
/// // Construct length-2 sequence of hypergraphs of the following form. (in this
/// // case there are actually just simple undirected graphs)
/// //
/// // hypergraph 0:
/// //
/// // 1 ---- 3 ---- 5
/// // |      |      |
/// // |      |      |
/// // 0 ---- 2 ---- 4
/// //
/// // hypergraph 1:
/// //
/// //        3 ---- 5 ---- 7
/// //        |      |      |
/// //        |      |      |
/// //        2 ---- 4 ---- 6
/// //
/// // --------------------------------------------------------------------
/// 
/// let number_of_ladders                           =   2;
/// let holes_per_ladder                            =   2;
/// 
/// let mut hypergraphs                             =   Vec::new();
/// for offset_from_left in 0 .. number_of_ladders {
///     hypergraphs.push( 
///         sideways_ladder_edges(offset_from_left, holes_per_ladder)  
///     );   
/// }
/// 
/// // Compute the barcode for the corresponding zigzag:
/// //
/// // hypergraph_0 ----> ( hypergraph_0  U  hypergraph_1 ) <---- hypergraph_1
/// //
/// // --------------------------------------------------------------------
/// 
/// let ring_operator                               =   BooleanField::new();
/// let max_homology_dimension                      =   2;
/// let interval_decomposition                      =   interval_decomposition_for_zigzag_of_hypgeraph_unions(hypergraphs, ring_operator, max_homology_dimension);
/// let barcode                                     =   interval_decomposition.iter()
///                                                         .map(| (dim, submodule)| ( dim.clone(), submodule.interval_endpoints() ) )
///                                                         .collect::<Vec<_>>();
/// 
/// // Verify that the barcode is correct
/// // Note that each pair (a,b) represents a half-open interval of form [a,b)
/// // The actual barcode is
/// //
/// //             
/// // dimension 0:
/// //       
/// //      0 ---- 1 ---- 2             
/// //
/// // dimension 1  
/// //      
/// //      0 ---- 1
/// //      0 ---- 1 ---- 2
/// //             1 ---- 2
/// //
/// // --------------------------------------------------------------------
///         
/// let barcode_ground_truth: Vec< (usize,(usize,usize)) >                        =   vec![
///                                                         ( 0, (0,3) ),
///                                                         ( 1, (0,2) ),
///                                                         ( 1, (0,3) ),
///                                                         ( 1, (1,3) ),                                                                                                                                                                                               
///                                                     ];
/// 
/// assert_eq!( barcode, barcode_ground_truth );
/// ```
pub fn interval_decomposition_for_zigzag_of_hypgeraph_unions< RingOperator >( 
        hypergraphs:                Vec< Vec< Vec< usize > > >, 
        ring_operator:              RingOperator,
        max_homology_dimension:     usize,
    )
    -> 
    Vec< ( usize, IntervalSubmoduleBasis< RingOperator::Element > ) >

    where
        RingOperator:               Clone +  DivisionRingOperations + std::fmt::Debug,
        RingOperator::Element:      Clone + std::fmt::Debug + Eq + Ord,
{
    // compute number of vertices and arrows
    let n_hypergraphs                               =   hypergraphs.len();
    

    // handle the edge case where the number of hypergraphs is 0 or 1
    // ========================================================================    
    if n_hypergraphs == 0 { // if there are no hypergraphs then short circuit
        return Vec::with_capacity(0);
    } else if n_hypergraphs == 1 {
        let factored                                    =   factor_dowker_complex(hypergraphs[0].clone(), ring_operator.clone(), max_homology_dimension );
        let mut submodules  =   Vec::with_capacity(factored.homology_indices().len());
        for basis_vector in factored.homology_basis() {
            let basis_vector: Vec<_> = basis_vector.collect();
            let dimension = chain_dimension(&basis_vector);
            let leftmost_vertex = 0;
            let submodule   =   IntervalSubmoduleBasis::new( leftmost_vertex , vec![basis_vector] );
            submodules.push( (dimension,submodule) );
        }
        return submodules 
    }

    // calculate the cospans
    // ========================================================================

    let start = Instant::now();
    
    let mut cospans                                     =   Vec::with_capacity( n_hypergraphs - 1 );
    for ( hypergraph_a, hypergraph_b ) in hypergraphs.into_iter().tuple_windows() {   
        cospans.push(
            induced_cospan( hypergraph_a, hypergraph_b, ring_operator.clone(), max_homology_dimension )
        );         
    }

    let duration = start.elapsed();                     // Stop the timer
    println!("Time to build cospans: {:?}", duration);     // Print the elapsed time    

    // return the decomposition
    // ========================================================================    

    cospan_pipeline::interval_decomposition_for_zigzag_of_cospans(cospans, ring_operator)
}







pub fn interval_decomposition_for_zigzag_of_hypgeraph_unions_WITH_SPANS< RingOperator >( 
        hypergraphs:                Vec< Vec< Vec< usize > > >, 
        ring_operator:              RingOperator,
        max_homology_dimension:     usize,
    )
    -> 
    Vec< ( usize, IntervalSubmoduleBasis< RingOperator::Element > ) >

    where
        RingOperator:               Clone +  DivisionRingOperations + std::fmt::Debug,
        RingOperator::Element:      Clone + std::fmt::Debug + Eq + Ord,
{
    // compute number of vertices and arrows
    let n_hypergraphs                               =   hypergraphs.len();


    // handle the edge case where the number of hypergraphs is 0 or 1
    // ========================================================================    
    if n_hypergraphs == 0 { // if there are no hypergraphs then short circuit
        return Vec::with_capacity(0);
    } else if n_hypergraphs == 1 {
        let factored                                    =   factor_dowker_complex(hypergraphs[0].clone(), ring_operator.clone(), max_homology_dimension );
        let mut submodules  =   Vec::with_capacity(factored.homology_indices().len());
        for basis_vector in factored.homology_basis() {
            let basis_vector: Vec<_> = basis_vector.collect();
            let dimension = chain_dimension(&basis_vector);
            let leftmost_vertex = 0;
            let submodule   =   IntervalSubmoduleBasis::new( leftmost_vertex , vec![basis_vector] );
            submodules.push( (dimension,submodule) );
        }
        return submodules 
    }

    // calculate the cospans
    // ========================================================================

    let start = Instant::now();

    let mut spans                                     =   Vec::with_capacity( n_hypergraphs - 1 );
    for ( counter, ( hypergraph_a, hypergraph_b ) ) in hypergraphs.into_iter().tuple_windows().enumerate() {   
        let return_left_basis       =   counter == 0;
        spans.push(
            induced_span( hypergraph_a, hypergraph_b, ring_operator.clone(), max_homology_dimension, return_left_basis )
        );         
    }

    let duration = start.elapsed();                     // Stop the timer
    println!("Time to build spans: {:?}", duration);     // Print the elapsed time    

    // return the decomposition
    // ========================================================================    

    span_pipeline::interval_decomposition_for_zigzag_of_spans(spans, ring_operator)
}
































#[cfg(test)]
mod tests {



    

    #[test]
    fn doc_test_sideways_ladder() {
        use crate::algebra::rings::types::field_prime_order::BooleanField;
        use crate::topology::simplicial::from::relation::sideways_ladder_edges;
        use crate::algebra::zigzag::hypergraph_pipeline::interval_decomposition_for_zigzag_of_hypgeraph_unions;


        // Construct length-2 sequence of hypergraphs of the following form. (in this
        // case there are actually just simple undirected graphs)
        //
        // hypergraph 0:
        //
        // 1 ---- 3 ---- 5
        // |      |      |
        // |      |      |
        // 0 ---- 2 ---- 4
        //
        // hypergraph 1:
        //
        //        3 ---- 5 ---- 7
        //        |      |      |
        //        |      |      |
        //        2 ---- 4 ---- 6
        //
        // --------------------------------------------------------------------

        let number_of_ladders                           =   2;
        let holes_per_ladder                            =   2;

        let mut hypergraphs                             =   Vec::new();
        for offset_from_left in 0 .. number_of_ladders {
            hypergraphs.push( 
                sideways_ladder_edges(offset_from_left, holes_per_ladder)  
            );   
        }

        // Compute the barcode for the corresponding zigzag:
        //
        // hypergraph_0 ----> ( hypergraph_0  U  hypergraph_1 ) <---- hypergraph_1
        //
        // --------------------------------------------------------------------

        let ring_operator                               =   BooleanField::new();
        let max_homology_dimension                      =   2;
        let interval_decomposition                      =   interval_decomposition_for_zigzag_of_hypgeraph_unions(hypergraphs, ring_operator, max_homology_dimension);
        let barcode                                     =   interval_decomposition.iter()
                                                                .map(| (dim, submodule)| ( dim.clone(), submodule.interval_endpoints() ) )
                                                                .collect::<Vec<_>>();
        
        // Verify that the barcode is correct
        // Note that each pair (a,b) represents a half-open interval of form [a,b)
        // The actual barcode is
        //
        //             
        // dimension 0:
        //       
        //      0 ---- 1 ---- 2             
        //
        // dimension 1  
        //      
        //      0 ---- 1
        //      0 ---- 1 ---- 2
        //             1 ---- 2
        //
        // --------------------------------------------------------------------
                
        let barcode_ground_truth: Vec< (usize,(usize,usize)) >                        =   vec![
                                                                ( 0, (0,3) ),
                                                                ( 1, (0,2) ),
                                                                ( 1, (0,3) ),
                                                                ( 1, (1,3) ),                                                                                                                                                                                               
                                                            ];

        assert_eq!( barcode, barcode_ground_truth );

    }





    #[test]
    fn test_sideways_ladder() {
        use crate::algebra::rings::types::field_prime_order::BooleanField;
        use crate::topology::simplicial::from::relation::sideways_ladder_edges;
        use crate::algebra::zigzag::hypergraph_pipeline::interval_decomposition_for_zigzag_of_hypgeraph_unions;


        let number_of_ladders                           =   3;
        let holes_per_ladder                            =   2;

        let mut hypergraphs                             =   Vec::new();
        for offset_from_left in 0 .. number_of_ladders {
            hypergraphs.push( 
                sideways_ladder_edges(offset_from_left, holes_per_ladder)  
            );   
        }

        // Compute the barcode for the corresponding zigzag:
        // hypergraph_0 ----> ( hypergraph_0  U  hypergraph_1 ) <---- hypergraph_1
        let ring_operator                               =   BooleanField::new();
        let max_homology_dimension                      =   2;        
        let interval_decomposition                      =   interval_decomposition_for_zigzag_of_hypgeraph_unions(hypergraphs, ring_operator, max_homology_dimension);
        let barcode                                     =   interval_decomposition.iter()
                                                                .map(| (dim, submodule)| ( dim.clone(), submodule.interval_endpoints() ) )
                                                                .collect::<Vec<_>>();

        // Verify that the barcode is correct
        // Note that each pair (a,b) represents a half-open interval of form [a,b)
        //
        // Here's a "diagram" of the ladders, where each * represents a hole
        //
        //  laddder 0        : * *
        //  laddder 1 (union): * * *
        //  laddder 2        :   * *
        //  laddder 3 (union):   * * *
        //  laddder 4        :     * *
        //
        // By in specting this we can see that the barcode should be
        //
        // Dimension 0:
        //      0 ---- 1 ---- 2 ---- 3 ---- 4
        //
        // Dimension 2:
        //      0 ---- 1 
        //      0 ---- 1 ---- 2 ---- 3
        //             1 ---- 2 ---- 3 ---- 4
        //                           3 ---- 4
                                     
        let barcode_ground_truth                        =   vec![
                                                                ( 0, (0,5) ),
                                                                ( 1, (0,2) ),
                                                                ( 1, (0,4) ),
                                                                ( 1, (1,5) ),
                                                                ( 1, (3,5) ),
                                                            ];

        assert_eq!( barcode, barcode_ground_truth );

    }    


    /// Similar to preceding test, but for every hypergraph we add a hyperdge that contains all the odd vertices.
    #[test]
    fn test_sideways_ladder_with_top_blob() {
        use num::Integer;
        use std::collections::HashSet;
        use crate::algebra::rings::types::field_prime_order::BooleanField;
        use crate::topology::simplicial::from::relation::sideways_ladder_edges;
        use crate::algebra::zigzag::hypergraph_pipeline::interval_decomposition_for_zigzag_of_hypgeraph_unions;


        let number_of_ladders                           =   3;
        let holes_per_ladder                            =   2;

        let mut hypergraphs                             =   Vec::new();


        // build the ladder hypergraphs
        for offset_from_left in 0 .. number_of_ladders {
            hypergraphs.push( 
                sideways_ladder_edges(offset_from_left, holes_per_ladder)  
            );   
        }

        // add the blob
        let mut blob                                    =   HashSet::new();
        blob.extend( 
            hypergraphs.iter()
                .map(|x| x.iter().map(|y| y.iter().cloned() ).flatten() )
                .flatten()
                .filter(|x| x.is_odd() )
            );
        let mut blob: Vec<_>                                    =   blob.into_iter().collect();
        blob.sort();

        for hypergraph in hypergraphs.iter_mut() {
            hypergraph.push( blob.clone() );
        }            

        // Compute the barcode for the corresponding zigzag:
        // hypergraph_0 ----> ( hypergraph_0  U  hypergraph_1 ) <---- hypergraph_1
        let ring_operator                               =   BooleanField::new();
        let max_homology_dimension                      =   2;        
        let interval_decomposition                      =   interval_decomposition_for_zigzag_of_hypgeraph_unions(hypergraphs, ring_operator, max_homology_dimension);
        let barcode                                     =   interval_decomposition.iter()
                                                                .map(| (dim, submodule)| ( dim.clone(), submodule.interval_endpoints() ) )
                                                                .collect::<Vec<_>>();

        // Verify that the barcode is correct
        // Note that each pair (a,b) represents a half-open interval of form [a,b)
        //
        // Here's a "diagram" of the ladders, where each * represents a hole
        //
        //  laddder 0        : * *
        //  laddder 1 (union): * * *
        //  laddder 2        :   * *
        //  laddder 3 (union):   * * *
        //  laddder 4        :     * *
        //
        // By in specting this we can see that the barcode should be
        //
        // Dimension 0:
        //      0 ---- 1 ---- 2 ---- 3 ---- 4
        //
        // Dimension 2:
        //      0 ---- 1 
        //      0 ---- 1 ---- 2 ---- 3
        //             1 ---- 2 ---- 3 ---- 4
        //                           3 ---- 4
                                     

        let barcode_ground_truth                        =   vec![
                                                                ( 0, (0,5) ),
                                                                ( 1, (0,2) ),
                                                                ( 1, (0,4) ),
                                                                ( 1, (1,5) ),
                                                                ( 1, (3,5) ),
                                                            ];                                     

        assert_eq!( barcode, barcode_ground_truth );

    }        

}
    




