//! Example homology calculations.
//! 
//! This section is under construction.  See the [Dowker complex module](oat_rust::topology::simplicial::from::relation) for a worked example.


// //! An example homology calculation.
// //! 
// //! Here we compute a basis for homology in dimension 1,
// //! for the simplicial complex on vertex set `{0,1,2,3}`
// //! whose maximal faces are `{0,1,2}, {0,3}, {1,3}, {2,3}`.
// //! The coefficient field is the finite field of order 3.
// //! 
// //! **To run this example on your desktop computer** first checkout the
// //! [quick start tutorial in OAT]() for instructions
// //! on installing Rust and running a program.  As part of this process,
// //! you'll create a new folder that contains a file called `main.rs`.  Inside
// //! `main.rs` is some text that reads `fn main{ .. }`.  Delete everything
// //! between `{` and `}`, and paste in the following:
// //! 
// //!
// //! ```
// //! use std::collections::HashSet;
// //! use std::iter::FromIterator;
// //! use oat_rust::algebra::matrices::query::MatrixOracle;
// //! use oat_rust::algebra::matrices::operations::umatch::row_major::{Umatch::new_with_compression};
// //! use oat_rust::topology::simplicial::from::relation::DowkerComplex;
// //! use oat_rust::topology::simplicial::simplices::vector::{dimension_0_through_d_simplices_in_reverse_dimensionwise_lexicographic_order_iter, dimension_d_simplices_in_lexicographic_order_iter};
// //! use oat_rust::utilities::order::OrderOperatorAutoLt;        
// //! use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;
// //! 
// //! // Define the ring operator for the finite field of order 3.
// //! // You can use this object to perform arithmetic operations, e.g., 
// //! // to add 1 and 1, you can run `ring_operator.add(1,1)`.
// //! let ring_operator          =   PrimeOrderField::new(3);
// //! 
// //! // We will build a dowker complex.
// //! // A dowker complex is defined by a vertex set V and a family S
// //! // of subsets of V.  A subset of V forms a simplex iff it is 
// //! // a subset of some element of S.  We refer to the elements 
// //! // of S as "dowker simplices".
// //! 
// //! // Each dowker simplex is represented by a vector, and we store the
// //! // list of all such simplices inside a larger vector.
// //! let dowker_simplices_vec_format =   vec![    
// //!                                         vec![0,1,2], 
// //!                                         vec![0,3], 
// //!                                         vec![1,3], 
// //!                                         vec![2,3]  
// //!                                     ];                                
// //! 
// //! // For certain calculations it's better to store a simplex as a *set*,
// //! // rather than a vector.  This object stores the same collection of 
// //! // simplices, represented by sets.
// //! let dowker_simplices_set_format: Vec<_>  =   dowker_simplices_vec_format
// //!                                                         .iter()
// //!                                                         .cloned()
// //!                                                         .map( |x| HashSet::from_iter( x ) )
// //!                                                         .collect();
// //! 
// //! // Build the boundary matrix.
// //! // This is a lazy object that generates rows/columns of the boundary matrix, on demand.
// //! let boundary_matrix = DowkerComplex::new( dowker_simplices_set_format, ring_operator.clone() );
// //! 
// //! // This iterates over simplices in descending order of 
// //! // dimension (first) and descending lexicographic order (second).
// //! // When computing homology of dimension d, we only need to iterate
// //! // over simplices of dimension d and below.
// //! let iter_row_index = dimension_0_through_d_simplices_in_reverse_dimensionwise_lexicographic_order_iter( &dowker_simplices_vec_format, 1 );
// //! 
// //! // Compute a umatch factorization of the boundary matrix.
// //! // For details on what this factorization entails, see the paper 
// //! // "U-match factorization: sparse homological algebra, lazy cycle representatives, and dualities in persistent (co) homology"
// //! // by Hang, Giusti, Ziegelmeier, and Henselman-Petrusek.  You can also check out the
// //! // OAT documentation for `Umatch`.
// //! let umatch = Umatch::new_with_compression(
// //!         boundary_matrix, 
// //!         iter_row_index, 
// //!         ring_operator.clone(), 
// //!         OrderOperatorAutoLt::new(), 
// //!         OrderOperatorAutoLt::new(),
// //!     );
// //! 
// //! // Get the source COMB (cf the  paper on umatch factorization)
// //! let comb_source = umatch.comb_source();
// //! 
// //! // Get the matching array (cf the  paper on umatch factorization)
// //! let matching = umatch.generalized_matching_matrix_ref();
// //! 
// //! // The set {columns of the source comb that are not matched upward or downward} 
// //! // forms a basis for homology
// //! let dim = 1;
// //! let mut betti = 0;
// //! 
// //! // Print the basis vectors
// //! println!(""); // an empty line, for spacing        
// //! println!("Each of the following lines represents a basis vector for homology in dimension {:?}", dim);
// //! for simplex in dimension_d_simplices_in_lexicographic_order_iter( &dowker_simplices_vec_format, dim ).unwrap() {
// //!     if matching.has_a_match_for_row_index( &simplex ) { continue }
// //!     if matching.has_a_match_for_column_index( &simplex ) { continue }        
// //!     let basis_vec = comb_source.column_reverse( simplex );
// //!     let basis_vec: Vec<_> = basis_vec.collect();
// //!     println!("basis vector {:?}: {:?}", betti, basis_vec);
// //!     betti += 1;
// //! }
// //! 
// //! // Print the betti number
// //! println!(""); // an empty line, for spacing
// //! println!("The betti number in dimension {:?} is {:?}.", dim, betti);
// //! println!(""); // an empty line, for spacing    
// //! ```    
// //! 
// //! Make sure that all changes are saved to `main.rs`, then run the
// //! program as described in the [quick start tutorial in OAT]().
// //! This should print the following:
// //! 
// //! ```bash
// //! $ Each of the following lines represents a basis vector for homology in dimension 1
// //! $ basis vector 0: [([1, 3], 1), ([0, 3], 2), ([0, 1], 2)]
// //! $ basis vector 1: [([2, 3], 1), ([0, 3], 2), ([0, 2], 2)]
// //! $ 
// //! $ The betti number in dimension 1 is 2
// //! ```
// //! 
// //! 
// //! # Change the coefficient ring
// //! 
// //! OAT has a number of different [predefined coefficient rings](oat_rust::rings), which you can substitute into
// //! the example above in order to calculate homology with different coefficients.  Simply replace the
// //! line 
// //! ```
// //! use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;       
// //! let ring_operator   =   PrimeOrderField::new(3);
// //! ``` 
// //! with one of the `let ring_operator = ...` lines listed under *Predefined Rings*, [here](oat_rust::rings).
// //! 
// //! 
// //! 
// //! 
// // //! 
// // //! In the example above, we defined the ring operator in two steps: first, importing the definition of the operator
// // //! using a `use` statement, then building the operator itself using a `new` function:
// // //! 
// // //! ```
// // //! // import the definition
// // //! use oat_rust::rings::types::field_prime_order::PrimeOrderField; 
// // //! // build the operator
// // //! let ring_operator   =   PrimeOrderField::new(3); 
// // //! ```
// // //! 
// // //! If we preffered, we could have combined these steps:
// // //! 
// // //! ```
// // //! let ring_operator   =   oat_rust::rings::types::field_prime_order::PrimeOrderField::new(3);
// // //! ```
// // //! 
// // //! OAT has a number of different [predefined coefficient rings](oat_rust::rings) (see the section titled *Predefined rings*), which you can substitute into
// // //! the example above, to calculate homology with different coefficients.



// #[cfg(test)]
// mod doc_test_drafts {
//     use itertools::Itertools;
//     use num::rational::Ratio;
//     use oat_rust::utilities::sequences_and_ordinals::SortedVec;

//     use crate::topology::simplicial::simplices::vector::{dimension_0_through_d_simplices_in_reverse_dimensionwise_lexicographic_order_iter, dimension_d_simplices_in_lexicographic_order_iter, dimension_0_through_d_simplices_in_dimensionwise_lexicographic_order_iter, dimension_d_simplices_in_reverse_lexicographic_order_iter};        


//     #[test]
//     fn compute_homology_var() {

//         use itertools::Itertools;
//         use std::collections::HashSet;
//         use std::iter::FromIterator;
//         use crate::topology::simplicial::from::relation::DowkerComplex;
//         use crate::topology::simplicial::simplices::vector::{dimension_0_through_d_simplices_in_reverse_dimensionwise_lexicographic_order_iter, dimension_d_simplices_in_lexicographic_order_iter};  
//         use oat_rust::chains::factored::factor_boundary_matrix;     
//         use crate::algebra::matrices::query::MatrixOracle;
//         use crate::algebra::matrices::operations::umatch::row_major::{Umatch::new_with_compression};
//         use oat_rust::utilities::order::OrderOperatorAutoLt;        
//         use oat_rust::rings::types::field_prime_order::PrimeOrderField;

//         // Parameters
//         // ----------

//         // Define the maximum homology dimensino we want to compute.
//         let max_homology_dimension                  =   2;

//         // Define the ring operator for the finite field of order 3.
//         // You can use this object to perform arithmetic operations, e.g., 
//         // to add 1 and 1, you can run `ring_operator.add(1,1)`.
//         let ring_operator          =   PrimeOrderField::new(3);

//         // We will build a dowker complex.
//         // A dowker complex is defined by a vertex set V and a family S
//         // of subsets of V.  A subset of V forms a simplex iff it is 
//         // a subset of some element of S.  We refer to the elements 
//         // of S as "dowker simplices".

//         // Each dowker simplex is represented by a SortedVec of vertices.
//         // We store the list of all such simplices inside a larger vector.
//         let dowker_simplices 
//             =   vec![    
//                         vec![0,1,2], 
//                         vec![0,3], 
//                         vec![1,3], 
//                         vec![2,3]  
//                     ]
//                     .into_iter()
//                     .map( |x| SortedVec::new(x) ) 
//                     .collect_vec();

//         //  Boundary matrix
//         //  ---------------

//         // This is a lazy object that generates rows/columns of the boundary matrix, on demand.
//         let boundary_matrix = DowkerComplex::new( dowker_simplices, ring_operator.clone() );

//         //  Simplex iterators
//         //  -----------------

//         // An iterator that runs over all triangles in the complex, in ascending 
//         // lexicographic order
//         let triangles_ascending_order = dimension_d_simplices_in_lexicographic_order_iter( &dowker_simplices, 2);

//         // An iterator that runs over all edges in the complex, in descending 
//         // lexicographic order
//         let triangles_descending_order = dimension_d_simplices_in_reverse_lexicographic_order_iter( &dowker_simplices, 2);   

//         // An iterator that runs over simplices of dimension 0 through max_homology_dimension,
//         // ordered first by dimension (ascending) and second by lexicographic order (descending)
//         let row_indices = boundary_matrix.simplices_in_row_reduction_order( max_homology_dimension );

//         //  Homology computation (by matrix factorization)
//         //  ----------------------------------------------

//         // Factor the boundary matrix
//         let factored    =   DifferentialUmatch::new( 
//                                 boundary_matrix, 
//                                 ring_operator, 
//                                 OrderOperatorAutoLt::new(), 
//                                 row_indices,
//                             );

//         //  Printing results
//         //  ----------------

//         // Betti numbers.  For this computation we have to provide a
//         // function that assigns a dimension to each index (i.e. to each simplex)
//         let homology_dimensions   =   factored.homology_dimensions(|x| x.len() as isize -1 ); 
//         for dim in 0 .. 2 {
//             println!(
//                     // we'll insert two values into this string
//                     "The betti number in dimension {:?} is {:?}.",
//                     dim,                  // the dimension
//                     homology_dimensions         // and the betti number
//                         .get( & dim )     // this looks up a value in the hashmap
//                         .unwrap_or( & 0)  // if the hashmap doesn't have the value, then use 0
//                 )
//         }

//         // Cycle representatives for homology
//         println!(
//                 "The following are basis vectors for homology in dimensions {:?} through {:?}",
//                 0,
//                 max_homology_dimension,
//             );
//         for (cycle_number, cycle) in factored.homology_basis().enumerate() {
//             // `cycle` is an iterator.  For convenience, collect the elements of the
//             // iterator into a Rust vector.
//             let cycle: Vec<_> = cycle.collect();
//             println!("Cycle number {:?} is:", cycle_number);
//             println!("{:?}", cycle );          
//         }
//     }




//     #[test]
//     fn compute_homology_projective() {

//         use itertools::Itertools;
//         use std::collections::HashSet;
//         use std::iter::FromIterator;
//         use crate::topology::simplicial::from::relation::DowkerComplex;
//         use crate::topology::simplicial::simplices::vector::{dimension_0_through_d_simplices_in_reverse_dimensionwise_lexicographic_order_iter, dimension_d_simplices_in_lexicographic_order_iter};  
//         use oat_rust::chains::factored::factor_boundary_matrix;     
//         use crate::algebra::matrices::query::MatrixOracle;
//         use crate::algebra::matrices::operations::umatch::row_major::{Umatch::new_with_compression};
//         use oat_rust::utilities::order::OrderOperatorAutoLt;        
//         use oat_rust::rings::types::field_prime_order::PrimeOrderField;

//         // Parameters
//         // ----------

//         // Define the maximum homology dimensino we want to compute.
//         let max_homology_dimension                  =   2;

//         // Define the ring operator for the finite field of order 3.
//         // You can use this object to perform arithmetic operations, e.g., 
//         // to add 1 and 1, you can run `ring_operator.add(1,1)`.
//         let ring_operator          =   PrimeOrderField::new(3);

//         // We will build a dowker complex.
//         // A dowker complex is defined by a vertex set V and a family S
//         // of subsets of V.  A subset of V forms a simplex iff it is 
//         // a subset of some element of S.  We refer to the elements 
//         // of S as "dowker simplices".

//         // Each dowker simplex is represented by a SortedVec of vertices.
//         // We store the list of all such simplices inside a larger vector.
//         let dowker_simplices 
//             =   vec![ 
//                     vec![0, 1, 2], vec![0, 3, 4], vec![1, 3, 5], vec![2, 4, 5], vec![0, 2, 3], vec![2, 3, 5], vec![1, 2, 4], vec![0, 1, 5], vec![1, 3, 4], vec![0, 4, 5]                                       
//                 ]
//                 .into_iter()
//                 .map( |x| SortedVec::new(x) ) 
//                 .collect_vec();

//         //  Boundary matrix
//         //  ---------------

//         // This is a lazy object that generates rows/columns of the boundary matrix, on demand.
//         let boundary_matrix = DowkerComplex::new( dowker_simplices, ring_operator.clone() );

//         //  Simplex iterators
//         //  -----------------

//         // An iterator that runs over all triangles in the complex, in ascending 
//         // lexicographic order
//         let triangles_ascending_order = dimension_d_simplices_in_lexicographic_order_iter( &dowker_simplices, 2);

//         // An iterator that runs over all edges in the complex, in descending 
//         // lexicographic order
//         let triangles_descending_order = dimension_d_simplices_in_reverse_lexicographic_order_iter( &dowker_simplices, 2);   

//         // An iterator that runs over simplices of dimension 0 through max_homology_dimension,
//         // ordered first by dimension (ascending) and second by lexicographic order (descending)
//         let row_indices = boundary_matrix.simplices_in_row_reduction_order( max_homology_dimension );

//         //  Homology computation (by matrix factorization)
//         //  ----------------------------------------------

//         // Factor the boundary matrix
//         let factored    =   DifferentialUmatch::new( 
//                                 boundary_matrix, 
//                                 ring_operator, 
//                                 OrderOperatorAutoLt::new(), 
//                                 row_indices,
//                             );

//         //  Printing results
//         //  ----------------

//         // Betti numbers.  For this computation we have to provide a
//         // function that assigns a dimension to each index (i.e. to each simplex)
//         let homology_dimensions   =   factored.homology_dimensions(|x| x.len() as isize -1 ); 
//         for dim in 0 .. 2 {
//             println!(
//                     // we'll insert two values into this string
//                     "The betti number in dimension {:?} is {:?}.",
//                     dim,                  // the dimension
//                     homology_dimensions         // and the betti number
//                         .get( & dim )     // this looks up a value in the hashmap
//                         .unwrap_or( & 0)  // if the hashmap doesn't have the value, then use 0
//                 )
//         }

//         // Cycle representatives for homology
//         println!(
//                 "The following are basis vectors for homology in dimensions {:?} through {:?}",
//                 0,
//                 max_homology_dimension,
//             );
//         for (cycle_number, cycle) in factored.homology_basis().enumerate() {
//             // `cycle` is an iterator.  For convenience, collect the elements of the
//             // iterator into a Rust vector.
//             let cycle: Vec<_> = cycle.collect();
//             println!("Cycle number {:?} is:", cycle_number);
//             println!("{:?}", cycle );          
//         }
//     }
// }