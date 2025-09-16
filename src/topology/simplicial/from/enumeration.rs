//! Simplicial complex represented by the set of all its simplices
//! 
//! This module contains a function that accepts the sequence of all simplices as input, and returns an integer-indexed matrix, as output.

use itertools::Itertools;
use crate::algebra::rings::traits::RingOperations;
use crate::utilities::sequences_and_ordinals::{BijectiveSequence};
use std::fmt::Debug;
use std::hash::Hash;



//  ---------------------------------------------------------------------------
//  SIMPLEX BIMAP TO BOUNDARY - SIMPLICES AS VECTORS
//  ---------------------------------------------------------------------------

/// Returns a (integer indexed) boundary matrix of a simplicial complex, given (i) a bijection between the set of simplices `K` and the set of integers `{0, .., K-1}`, and (ii) a ring operator.
/// 
/// The `k`th row and column of the matrix correspond to the `k`th simplex, as enumerated by the bimap.
/// 
/// The constructed matrix has rows and columns indexed by `usize`, and coefficients of type `RingElt`.
pub fn boundary_matrix_from_simplex_bimap< Vertex, RingOperator >( 
            simplex_bimap:      & BijectiveSequence< Vec < Vertex > >,
            ring_operator:      RingOperator
        ) 
        ->
        Vec< Vec < (usize, RingOperator::Element) >>

        where   Vertex:             Ord + Hash + Clone + Debug,      
                RingOperator:       RingOperations,
{
    if simplex_bimap.is_empty() { return vec![] }

    let mut boundary            =   Vec::with_capacity( simplex_bimap.len() );  

    let mut simplex_dim;
    let mut simplex_num_verts;

    for simplex in simplex_bimap.vec_elements_in_order() {

        simplex_num_verts       =   simplex.len();
        simplex_dim             =   simplex_num_verts - 1;

        // no need to calculate boundaries of dim-0 cells
        if simplex_dim == 0 {
            boundary.push( Vec::with_capacity(0) );
            continue;
        }  

        let mut vec             =   Vec::with_capacity( simplex_num_verts );    // num_vertices = NUMBER OF FACETS

        for (facet_count, facet)  in simplex.iter().cloned().combinations( simplex_dim ).enumerate() {
            vec.push( 
                (
                    simplex_bimap.ordinal_for_element( &facet ).unwrap(),
                    ring_operator.minus_one_to_power( simplex_dim - facet_count )
                ) 
            )            
        }
        boundary.push( vec );
    }

    boundary

}




//  ---------------------------------------------------------------------------
//  SIMPLEX BIMAP TO BOUNDARY - SIMPLICES AS SIMPLEX STRUCTS
//  ---------------------------------------------------------------------------




// /// Similar to [boundary_matrix_from_simplex_bimap], but simplices are represented by objects of type `Simplex`, not `Vec`.
// pub fn boundary_matrix_from_complex_facets_simplexform< Vertex, RingOperator, RingElt >( 
//     simplex_bimap:      BijectiveSequence< Simplex< Vertex > >,
//     ring_operator:      RingOperator
// ) 
// ->
// Vec< Vec < (usize, RingElt) >>

// where   Vertex:    Ord + Hash + Clone + Debug,      
//         RingOperator:     Semiring< RingElt > + Ring< RingElt >,
// {
//     if simplex_bimap.is_empty() { return vec![] }

//     let mut boundary            =   Vec::with_capacity( simplex_bimap.len() );  

//     let mut state_iter          =   FacetIteratorNoReturnAscendingLex{
//                                     simplex: Simplex{ vertices: SortedVec::new(vec![]) },
//                                     facet: vec![],
//                                     deleted_vertex_index: None
//                                 };

//     let mut global_int_index;
//     let mut simplex_dim;
//     let mut simplex_num_verts;

//     for simplex in simplex_bimap.vec_elements_in_order().iter().cloned() {

//         simplex_dim             =   simplex.dimension();
//         simplex_num_verts       =   simplex.number_of_vertices();
//         state_iter.reinitialize_with_simplex( simplex );

//         let mut vec             =   Vec::with_capacity( simplex_num_verts );    // num_vertices = NUMBER OF FACETS

//         for i in 0 .. simplex_num_verts {
//             state_iter.next();
            
//             println!("{:?}", &state_iter);
//             println!("{:?}", &simplex_bimap);            

//             global_int_index    =   simplex_bimap.ordinal_for_element( & Simplex::from_vec( state_iter.facet.clone() ) ).unwrap();
//             vec.push( 
//                 (
//                     global_int_index.clone(),
//                     ring_operator.minus_one_to_power( simplex_dim - i )
//                 ) 
//             )
//         }
//         boundary.push( vec );
//     }

// boundary

// }



//  ===========================================================================
//  ===========================================================================
//  TESTS
//  ===========================================================================
//  ===========================================================================




#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use crate::{utilities::sequences_and_ordinals::SortedVec, topology::simplicial::simplices::vector::dimension_0_through_d_simplices_in_dimensionwise_lexicographic_order_list};

    #[test]
    fn test_bimap_to_boundary () {

        let ring                    =   crate::algebra::rings::types::native::RingOperatorForNativeRustNumberType::< f64 >::new();
        let complex_facets          =   vec![  SortedVec::new( vec![0,1,2] ).ok().unwrap() ];


        let bimap_sequential        =   BijectiveSequence::from_vec(
                                            dimension_0_through_d_simplices_in_dimensionwise_lexicographic_order_list( & complex_facets, 2 )
                                        ).unwrap();  

        let boundary                =   boundary_matrix_from_simplex_bimap( & bimap_sequential, ring );

        assert_eq!(     &   boundary,
                        &   vec![
                                    vec![],
                                    vec![],
                                    vec![],
                                    vec![(0, -1.0), (1, 1.0)],
                                    vec![(0, -1.0), (2, 1.0)],
                                    vec![(1, -1.0), (2, 1.0)],
                                    vec![(3, 1.0), (4, -1.0), (5, 1.0)]
                            ]
        )
    }    


}    