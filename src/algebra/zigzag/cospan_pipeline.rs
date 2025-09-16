use std::collections::HashSet;
use std::time::Instant;


use itertools::Itertools;
use num::Integer;


use crate::algebra::matrices::operations::MatrixOracleOperations;
use crate::algebra::matrices::types::packet::MatrixAlgebraPacket;
use crate::algebra::zigzag::hypergraph_pipeline::chain_dimension;
use crate::algebra::rings::traits::DivisionRingOperations; 
use crate::algebra::rings::types::field_prime_order::BooleanField;
use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
use crate::utilities::order::OrderOperatorByKey;



use super::cospans::Cospan;
use super::decompose::Diagonalization;
use super::hypergraph_pipeline::IntervalSubmoduleBasis;
use super::{cospans::{factor_dowker_complex, induced_cospan}, decompose::{QuiverReprsentation, SingleBarBasisVectorIndexLedger}};


use derive_getters::{Getters, Dissolve};
use derive_new::new;

// THIS JUST MAKES THE CODE EASIER TO READ
type Simplex        =   Vec<usize>;
type Chain<RingElement>          =   Vec< (Simplex, RingElement) >;








/// Returns an interval decomposition of the zigzag module of a sequence of cospans, where each cospan has an associated basis of simplicial chains
/// 
/// **This method consumes the cospans**, because we place them in a quiver representation.
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
pub fn interval_decomposition_for_zigzag_of_cospans< RingOperator >( 
        cospans:                    Vec< Cospan< RingOperator::Element > >, 
        ring_operator:              RingOperator,
    )
    -> 
    Vec< ( usize, IntervalSubmoduleBasis< RingOperator::Element > ) >

    where
        RingOperator:       Clone +  DivisionRingOperations + std::fmt::Debug,
        RingOperator::Element:        Clone + std::fmt::Debug + Eq + Ord,
{
    // compute number of vertices and arrows    
    if cospans.len() == 0 { // if there are no hypergraphs then short circuit
        return Vec::with_capacity(0);
    }    
    let n_cospans                                       =   cospans.len();
    let n_arrows                                        =   2 * n_cospans;
    let n_vertices                                      =   2 * n_cospans + 1.min( cospans.len() ); // the term 1.min( cospans.len() adds 1 if there is at least one cospan in the list; technically we don't need to worry about this because we can't get to this poitn in the code unless there is at least one cospan ... but it makes me feel better to have formulae that work everywhere

    // build quiver representation
    // ---------------------------

    let start = Instant::now();

    // define arrow directions
    let arrow_directions                                =   ( 0 .. n_arrows ).map(|i| i.is_even() ).collect::<Vec<_>>();    



    // initialize matrices and cycle bases
    let mut matrices                                    =   Vec::with_capacity(n_arrows);  
    let mut homology_cycle_bases                        =   Vec::with_capacity(n_vertices);

    // insert the basis for the leftmost vertex
    homology_cycle_bases.push( cospans[0].left_basis.clone() );
    
    
    // fill in matrices, vector space grades, and gradings
    for cospan in cospans {      
        matrices.push( cospan.left_morphism );
        matrices.push( cospan.right_morphism );
        homology_cycle_bases.push( cospan.center_basis );
        homology_cycle_bases.push( cospan.right_basis );            
    }

    let vector_space_dimensions                         =   homology_cycle_bases.iter().map(|basis| basis.number_of_rows() ).collect::<Vec<_>>();
    let quiver_representation                           =   QuiverReprsentation::new(
                                                                arrow_directions, 
                                                                matrices, 
                                                                vector_space_dimensions,
                                                                ring_operator.clone(),
                                                            );

                
    
    let duration = start.elapsed();                     // Stop the timer
    println!("Time to build quiver: {:?}", duration);     // Print the elapsed time

    // diagonalize the representation
    // ----------------------------------
    let diagonalization                                 =   quiver_representation.diagonalize().unwrap();



    // compute submodules

    let mut submodules                                  =   Vec::new();
    for single_bar_ledger in diagonalization.list_of_single_bar_basis_vector_index_ledgers() {

        let mut chains                                  =   Vec::with_capacity(single_bar_ledger.bar_length());
        for (vertex, basis_vector_index) in single_bar_ledger.iter() {
            let row       =   diagonalization.bases()[ vertex ].row_ref( basis_vector_index );

            let homology_cycle_basis_matrix   =     MatrixAlgebraPacket::with_default_order(
                                                                & homology_cycle_bases[ vertex ], 
                                                                ring_operator.clone()
                                                            );

            let chain   =   homology_cycle_basis_matrix.multiply_with_row_vector(row).collect();
            
            chains.push(chain);
        }

        let dimension                               =   chain_dimension( & chains[0] );        
        let submodule                               =   IntervalSubmoduleBasis{ 
                                                            leftmost_vertex:    single_bar_ledger.leftmost_vertex(),
                                                            basis_vectors:      chains,
                                                        };

        submodules.push( ( dimension, submodule ) );        
    }

    submodules.sort_by_key( |(dim,submodule)|  ( dim.clone(), submodule.interval_endpoints() )  );
    submodules
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

        // println!("max grade with bars: {:?}", diagonalization.max_nontrivial_grade() );
        // println!("diagonalization: {:#?}", & diagonalization );
        // println!("barcode: {:?}", barcode);

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

        // println!("max grade with bars: {:?}", diagonalization.max_nontrivial_grade() );
        // println!("diagonalization: {:#?}", & diagonalization );
        // println!("barcode: {:?}", barcode);

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
    




