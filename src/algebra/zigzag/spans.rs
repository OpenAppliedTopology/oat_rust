use std::cell::Cell;
use std::collections::HashMap;
use std::hash::Hash;
use std::convert::TryInto;


use crate::{algebra::{matrices::{operations::solve::triangle::{TriangularSolveForColumnVectorReverse, TriangularSolveForRowVector}, types::vec_of_vec::sorted::VecOfVec}, rings::traits::DivisionRingOperations}, utilities::sequences_and_ordinals::SortedVec};
use crate::algebra::rings::types::field_prime_order::BooleanField;


use itertools::{Diff, Itertools};
use serde::{Serialize, Deserialize};

use crate::{algebra::matrices::operations::umatch::differential::DifferentialUmatch, topology::simplicial::from::relation::DowkerComplex};













/// Takes a collection of hyperedges and a ring operator as input; returns a factored boundary matrix as output.
pub fn factor_dowker_complex< RingOperator >( 
        dowker_simplices:           Vec<Vec<usize>>,
        ring_operator:              RingOperator,
        max_homology_dimension:     usize,    
    )  
    -> 
    DifferentialUmatch<
        DowkerComplex<usize, RingOperator>,
    > 

    where
        RingOperator:       Clone + DivisionRingOperations,
{

    // Parameters
    // ----------
        
    // We will build a dowker complex.
    // A dowker complex is defined by a vertex set V and a family S
    // of subsets of V.  A subset of V forms a simplex iff it is 
    // a subset of some element of S.  We refer to the elements 
    // of S as "dowker simplices".
        
    // Each dowker simplex is represented by a SortedVec of vertices.
    // We store the list of all such simplices inside a larger vector.
    let dowker_simplices
        =   dowker_simplices
                .into_iter()
                .map( |x| SortedVec::new(x).unwrap() )  // we unwrap because `new` can return an error
                .collect_vec();
        
    //  Boundary matrix
    //  ---------------
        
    // This is a lazy object that generates rows/columns of the boundary matrix, on demand.
    let boundary_matrix = DowkerComplex::new( dowker_simplices.clone(), ring_operator.clone() );
        
        
    // An iterator that runs over simplices of dimension 0 through max_homology_dimension,
    // ordered first by dimension (ascending) and second by lexicographic order (descending)
    let row_indices = boundary_matrix.simplices_in_row_reduction_order( max_homology_dimension as isize ).collect::<Vec<_>>();
        
    //  Homology computation (by matrix factorization)
    //  ----------------------------------------------
        
    // Return the factored boundary matrix
    DifferentialUmatch::new( 
        boundary_matrix, 
        0, // min homology dimension
        max_homology_dimension as isize,
    )
    
}








#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Span< Coefficient > {
    pub left_morphism:      VecOfVec< usize, Coefficient >,
    pub right_morphism:     VecOfVec< usize, Coefficient >,
    pub left_basis_opt:     Option< VecOfVec< Vec<usize>, Coefficient > >,
    pub center_basis:       VecOfVec< Vec<usize>, Coefficient >,// Vec< ( Vec<usize>, Coefficient) >,
    pub right_basis:        VecOfVec< Vec<usize>, Coefficient >,// Vec< ( Vec<usize>, Coefficient) >,        
}








/// Returns matrices which represent linear operations on **ROWS**
pub fn induced_span< RingOperator >( 
        hypergraph_a:               Vec<Vec<usize>>, 
        hypergraph_b:               Vec<Vec<usize>>,
        ring_operator:              RingOperator,
        max_homology_dimension:     usize,
        return_left_basis:          bool,
    )
    ->
    Span< RingOperator::Element >
    where
        RingOperator:       Clone + DivisionRingOperations,


{

    let mut hypergraph_u                                =   hypergraph_a.clone();
    hypergraph_u.extend( hypergraph_b.iter().cloned() );


    // factor the boundary matrices
    let factored_boundary_matrix_a          =   factor_dowker_complex(hypergraph_a, ring_operator.clone(), max_homology_dimension );
    let factored_boundary_matrix_b          =   factor_dowker_complex(hypergraph_b, ring_operator.clone(), max_homology_dimension );
    let factored_boundary_matrix_u          =   factor_dowker_complex(hypergraph_u, ring_operator.clone(), max_homology_dimension );




    for (counter, decomposition) in vec![ & factored_boundary_matrix_a, & factored_boundary_matrix_b, & factored_boundary_matrix_u ].into_iter().enumerate() {
        println!("Homology basis for space number {}", counter);
        decomposition.print_cohomology_basis();
    }

    

    let homology_basis_in_u                             =   VecOfVec::from_iterable_of_iterables(factored_boundary_matrix_u.cohomology_basis()).ok().unwrap();



    let mut induced_maps        =   Vec::new();
    // let mut homology_cycle_bases                        =   Vec::new();

    // for each of the two bookend spaces (a and b), compute the induced map into u
    for bookend_boundary_matrix_decomposition in vec![ & factored_boundary_matrix_a, & factored_boundary_matrix_b ] {



        // place the homology basis vectors (or rather their indices) for space u into bijection with some integers
        //
        // NB: when OAT enumerates these indices it does an exhaustive search over the indices of the rows of the boundary matrix.
        //     this search visits simplices in the same order that we use when we factor the boundary matrix using the 
        //     cohomology algorithm. this means that we visit simplices in ascending order of dimension but *descending*
        //     lexicographic order within a dimension.
        //     therefore if simplices S and T have equal dimension and S ≤ T lexicographically, then  we will havee
        //     harmonic_basis_vector_index_in_bookend_to_ordinal( T ) ≤ harmonic_basis_vector_index_in_bookend_to_ordinal ( S )
        let mut harmonic_basis_vector_index_in_bookend_to_ordinal               =   HashMap::new();        
        
        for ( counter, index ) in bookend_boundary_matrix_decomposition.homology_indices().into_iter().enumerate() {
            harmonic_basis_vector_index_in_bookend_to_ordinal.insert( index, counter  );
        }




        let cohomology_basis_in_u                       =     VecOfVec::from_iterable_of_iterables(
            factored_boundary_matrix_u.cohomology_basis()
        ).ok().unwrap();


        // initialize the matrix that represents the induced map
        let mut induced_map                                 =   Vec::new();
    
        // iterate over basis vectors for the homology of the bookend space
        for cohomology_basis_vector_in_u        in    cohomology_basis_in_u.inner_vec_of_vec_ref() {

            let restriction_to_bookend                  =   cohomology_basis_vector_in_u
                                                                    .iter()
                                                                    .filter(
                                                                        |(simplex, _coefficient)|
                                                                        {
                                                                            // this is the simplex index of the entry
                                                                            let simplex     =   SortedVec::new(simplex.clone() ).unwrap();
                                                                            // check to see if any hyperedge in the bookend space contains this simplex
                                                                            let hypergraph  =   bookend_boundary_matrix_decomposition
                                                                                                    .boundary_matrix()
                                                                                                    .relation_rows();
                                                                            for hyperedge in hypergraph {
                                                                                if hyperedge.contains_subset( &simplex ) { return true }
                                                                            }
                                                                            // if it does not then we reject the entry
                                                                            return false
                                                                        }
                                                                    );

            // find the family of coefficients a_i such that harmonic_basis_vector_in_bookend = sum_i ( a_i  * the_ith_jordan_basis_vector_in_u )
            let corresponding_linear_combination_in_bookend       =     TriangularSolveForRowVector::solve(
                                                                            restriction_to_bookend.into_iter().cloned(), 
                                                                            bookend_boundary_matrix_decomposition.differential_comb_inverse()
                                                                        ).ok().unwrap();

            // remove the coefficients that correspond to non-harmonic basis vectors (i.e. basis vectors in the coboundary)
            let mut projection                              =   corresponding_linear_combination_in_bookend
                                                                    .filter_map( 
                                                                        | (simplex,coefficient) |
                                                                        {   
                                                                            harmonic_basis_vector_index_in_bookend_to_ordinal
                                                                                .get( & simplex )
                                                                                .cloned()
                                                                                .map( | ordinal | (ordinal, coefficient) )
                                                                        }
    
                                                                    )
                                                                    .collect::<Vec<_>>(); 

            projection.reverse(); // reversing entries places them in sorted order                                                          
            

            // push the combination as a row to the induced map
            // 
            // NB: if we had not re-indexed the projection with integers then we sould have had to reverse the vector, because TriangularsolveMinorDescend
            //     returns entries in descending order of index. however, as we noted in the comments above, the map harmonic_basis_vector_index_in_bookend_to_ordinal
            //     reverses the order of indices within a fixed dimension
            induced_map.push( projection );

        }

        // convert the induced map matrix into a VecOfVec, and add it to our list
        let induced_map                                 =   VecOfVec::new(induced_map).ok().unwrap();
        induced_maps.push( induced_map );
    }


    let right_morphism                                  =   induced_maps.remove(1);
    let left_morphism                                   =   induced_maps.remove(0);   
    // let right_basis                                     =   homology_cycle_bases.remove(1);
    let left_basis_opt                                  =   if return_left_basis {  
                                                                                                    let matrix = factored_boundary_matrix_a.cohomology_basis();
                                                                                                    let matrix = VecOfVec::from_iterable_of_iterables(matrix).ok().unwrap();
                                                                                                    Some( matrix )
                                                                                            } else {
                                                                                                None
                                                                                            }; 
    let right_basis                                     =    VecOfVec::from_iterable_of_iterables( factored_boundary_matrix_b.cohomology_basis() ).ok().unwrap();


    Span{
        left_morphism,
        right_morphism,
        left_basis_opt,
        center_basis:       homology_basis_in_u,
        right_basis,
    }

}






#[cfg(test)]
mod tests {    

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn doc_test_induced_cospan_two_cycles() {


        //
        // 3               5
        // 2       1       4
        //         0    
        //
        // the underlying simplicial complex we want to build looks like the result of
        // glueing two length-4 cycle graphs together by identifying a single edge
        // from each cycle

        // the "harmonic" basis of cycle representatives for the union (with vectors 
        // appearing in the same order in which they are returned by OAT iterators) is 
        // (omitting coefficients since we are working over the integers modulo 2):
        //
        // basis vector 0:  [0] + [1] + [2] + [3]    // a sum of four 0-simplices
        // basis vector 1:  [2,3]                    // this is a single 1-simplex
        // basis vector 2:  [4,5]                    // this is a single 1-simplex
        //
        // the "harmonic indices" associated with each vector is the last index of
        // any nonzero entry in the vector. so they are:
        //
        // basis vector 0:  [0] + [1] + [2] + [3] 
        // basis vector 1:  [2,3]
        // basis vector 2:  [4,5]
        //
        // the basis of cycle representatives for the other hypergraphs are
        //
        // hypergraph a:
        //   basis vector 0:    [0] + [1] + [2] + [3] 
        //   basis vector 1:    [2,3]
        // hypergraph b:
        //   basis vector 0:    [0] + [1] + [2] + [3] 
        //   basis vector 1:    [4,5]
        //
        // therefore the matrix representations for the maps induced by inclusion are
        //
        // hypergraph a:
        //   1 0
        //   0 1
        //   0 0 
        // hypergraph b:
        //   1 0 
        //   0 0
        //   0 1        

    


        // first half-filled cycle
        let hypergraph_a    =   vec![
                                    vec![ 0, 1 ],
                                    vec![ 0, 2 ],
                                    vec![ 1, 3 ],
                                    vec![ 2, 3 ],                                    
                                ];

        // second half-filled cycle                                
        let hypergraph_b    =   vec![
                                    vec![ 0, 1 ],
                                    vec![ 0, 4 ],
                                    vec![ 1, 5 ],
                                    vec![ 4, 5 ],
                                ];               

        // work over the two element field
        let ring_operator       =   BooleanField::new();

        // set max homology dimension
        let max_homology_dimension                      =   2;

        // compute the span
        let return_left_basis               =   false;
        let span              =   induced_span( hypergraph_a, hypergraph_b, ring_operator, max_homology_dimension, return_left_basis );

        // check the matrices
        let matrix_a            =   & span.left_morphism;
        let matrix_b            =   & span.right_morphism;        
        let matrix_a_ground_truth: VecOfVec< usize, bool >   =   VecOfVec::new( 
                                            vec![
                                                vec![ (0,true)             ], // this is for the 0-dimensional homology class
                                                vec![           (1,true )  ], // this is for the "lefthand" cycle  
                                                vec![                      ], // this is for the "righthand" cycle                                                                                                
                                            ]     
                                        ).ok().unwrap();
        let matrix_b_ground_truth: VecOfVec< usize, bool >   =   VecOfVec::new( 
                                            vec![
                                                vec![ (0,true)             ], // this is for the 0-dimensional homology class
                                                vec![                      ], // this is for the "lefthand" cycle  
                                                vec![           (1,true )  ], // this is for the "righthand" cycle                                                                                                
                                            ]     
                                        ).ok().unwrap();                                   
        
        assert_eq!(  matrix_a, & matrix_a_ground_truth );
        assert_eq!(  matrix_b, & matrix_b_ground_truth );       
    }

}