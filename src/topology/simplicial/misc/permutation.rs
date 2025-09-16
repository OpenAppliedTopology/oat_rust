//! Permutations on simplices induced by permutations on vertices

use crate::utilities::indexing_and_bijection::{compose_f_after_g, sort_perm_by_key};




//  ===========================================================================
//  ===========================================================================
//  SIMPLEX - AS - SORTED VECTOR
//  ===========================================================================
//  ===========================================================================

//  ---------------------------------------------------------------------------
//  PERMUTATION ON SIMPLICES INDUCED BY PERMUTATION ON VERTICES
//  ---------------------------------------------------------------------------


/// Obtain an argsort permutation on simplices from a permutation on vertices.
/// 
/// Given a vector f = [f0, .., fn] representing a function of form 
/// f: old_vertex_number -> new_vertex_number, obtain the vector 
/// g: old_simplex_number -> new_simplex_number; here "number" refers to 
/// lexicographic order. 
/// 
/// This function does not assume that the `simplex_sequence` has lexicogrphic 
/// order, but it **does** assume that the new simplex sequence has lexicographic
/// order.
pub fn simplex_perm_o2n_from_vertex_perm_o2n( 
    simplex_sequence:           &   Vec< Vec< usize > >,
    vertex_perm_old_to_new:     &   Vec< usize >
    ) 
    ->
    Vec< usize >
{
    // let mut get_key = 
    //     | x | 
    //     {
    //         let mut new = compose_f_after_g( &vertex_perm_old_to_new, x );
    //         new.sort();
    //         new
    //     };

    sort_perm_by_key( simplex_sequence, 
            | x | 
            {
                let mut new = compose_f_after_g( vertex_perm_old_to_new, x );
                new.sort();
                new
            } 
        )
    // // Create vector of new simplices
    // let mut new_simplex_sequence =  Vec::from_iter(
    //                                     simplex_sequence
    //                                         .iter()
    //                                         .cloned()
    //                                         .map(
    //                                             |x|
    //                                             Simplex{ 
    //                                                 vertices:  compose_f_after_g( &vertex_perm_old_to_new, &x )
    //                                             }
    //                                         )
    //                                 );

    // // We must remember to sort the new vertices                                    
    // for simplex in new_simplex_sequence.iter_mut() { simplex.vertices.sort()}

    // // Obtain the sort permutation
    // sort_perm( &new_simplex_sequence )
}







#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use crate::topology::simplicial::simplices::unweighted::Simplex;
    use crate::topology::simplicial::simplices::vector::dimension_0_through_d_simplices_in_dimensionwise_lexicographic_order_list;
    use crate::utilities::indexing_and_bijection::{inverse_perm};
    use crate::utilities::sequences_and_ordinals::SortedVec;
    use std::iter::FromIterator;

    #[test]
    fn test_simplex_perm_o2n_from_vertex_perm_o2n() {

        // sequence_old:          [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 2]]
        // sequence_old_permuted: [[0], [1], [3], [2], [0, 1], [0, 3], [0, 2], [1, 2]]
        // new_sequence:          [[0], [1], [2], [3], [0, 1], [0, 2], [0, 3], [1, 3]]
        // permutation: simplex old -> new 
        //                        [ 0,   1,   3,   2,   4,      6,      5,      7]


        let complex_facets          =   vec![  vec![0,1,2], vec![0, 3] ];
        let complex_facets  =   complex_facets.into_iter().map(|x| SortedVec::new(x).unwrap() ).collect();
        let simplex_sequence_old    =   dimension_0_through_d_simplices_in_dimensionwise_lexicographic_order_list( &complex_facets, 1);   
        let perm_v_o2n              =   vec![0, 1, 3, 2];

        let mut simplex_sequence_new =  Vec::from_iter(
            simplex_sequence_old
                .iter()
                .cloned()
                .map(
                    |x|
                    Simplex{ 
                        vertices:  compose_f_after_g( &perm_v_o2n, &x )
                    }
                )
        );        

        // We must remember to sort the new vertices                                    
        for simplex in simplex_sequence_new.iter_mut() { simplex.vertices.sort() }        

        // perm: simplex OLD -> NEW
        let perm_s_o2n              =   simplex_perm_o2n_from_vertex_perm_o2n( &simplex_sequence_old, &perm_v_o2n );
        // perm: simplex NEW -> OLD
        let perm_s_n2o               =   inverse_perm( &perm_s_o2n );

        let simplex_sequence_permuted    =   compose_f_after_g( &simplex_sequence_old, &perm_s_n2o );

        let mut simplex_sequence_permuted_vertex_translated     =   simplex_sequence_permuted.clone();
        for i in 0..simplex_sequence_permuted_vertex_translated.len() { simplex_sequence_permuted_vertex_translated[i] = compose_f_after_g( & perm_v_o2n, & simplex_sequence_permuted[i]) };
        for i in 0..simplex_sequence_permuted_vertex_translated.len() { simplex_sequence_permuted_vertex_translated[i].sort() };        
        

        println!("sequence_old:          {:?}",     & simplex_sequence_old );
        println!("sequence_old_permuted: {:?}",     & simplex_sequence_permuted );        
        println!("new_sequence:          {:?}",     & simplex_sequence_permuted_vertex_translated );     
        println!("permutation: simplex old -> new {:?}", & perm_s_o2n);           

    }


}    
