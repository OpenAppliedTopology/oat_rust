//! Solve `Ax = b` for `x`, where `A` is a matrix partially in echelon form.
//! 
//! - [Background](#background)
//! - [Solve xA = b](#solvexab), with several [alternate formulations](#alternate-formulations)
//! - [Solve Ay = c](#solveayc)
//! - [Quotients, remainders, and solutions](#qrsolutions)
//! - [Panics](#panics) will (eventually) occur when `b` or `c` aren't properly sorted.
//! 
//! # Background
//! 
//! Suppose we have an invertible upper triangular matrix `A` and a vector `c`.
//! 
//! We can solve `xA = c` for `c` using [back substitution](https://algowiki-project.org/en/Backward_substitution).
//! This process adds linear multiples of the rows of `A` to `c` iteratively, to eliminate leading entries, until
//! the vector has been reduced to zero. 
//! 
//! If we permute the rows of `A` we can still solve `xA = c` via back-substitution.  To do this efficiently, we
//! just need an efficient way to look up the row with leading entry `i`, for each column index `i = 0, 1, 2, ..`.
//! We refer to the map sending column indices to the corresponding row indices the **leading entry bijection**.
//! 
//! If we delete some rows from the permuted matrix `A` then it may no longer be possible to solve `xA = c` for x.
//! However, we can still apply the back substitution algorithm using just the rows that haven't been deleted --
//! stopping when we encounter a leading entry that can't be eliminated with the leading entry of a row from `A`.
//! This algorithm will produce a vector `q`.  We call
//! 
//! - `q` the **quotient**, and
//! - `c - Aq` the **remainder**
//! 
//! Since we've deleted some rows, there's no longer a nice bijection from column indices to row indices.  But there's
//! still a bijection from a subset of column indices to the row indices that haven't been deleted.  We call this the
//! 
//! - **partially defined leading entry bijection**
//! 
//! 
// //! We start with a matrix `A` and a vector `b`, were `b` is an iterator that runs over the
// //! structural nonzero entries of a sparse vector.  
// //! 
// //! Here is a conceptual example of what is meant by "echelon solve".  Suppose we have
// //! a matrix `A` and vector `b` with the sparse structure show below (`*` indicates a nonzero entry).
// // //! It can be shown that `b` lies in the row space of `A` iff `b[1:3]` lies in the
// // //! row space of the submatrix `A[2:4, 1:3]`.  
// //! If there exists a 
// //! vector `x` such that `xA = b`, then we can find such an `x` via "back substituion": first we add a scalar
// //! multiple of `A[2,:]` to clear the first entry in `b`, resulting in a vector `c`.
// //! Then we add a scalar multiple of `A[3,:]` to clear the first nonzero entry in `c` (which occurs in column 2).
// //! This process continues until we have eliminated all nonzero entries.  The
// //! solution `x` can then be constructed from the collection of scalar multipliers
// //! used to perform the elimination.
// //! 
// //! 
// //! ```text
// //!              1  2  3  4  5  6  7
// //!           -------------------------
// //!         1 |                        |
// //!         2 |  *  *  *  *  *  *  *   |
// //!   A  =  3 |     *  *  *  *  *  *   |
// //!         4 |        *  *  *  *  *   |
// //!         5 |                        |
// //!           -------------------------
// //! 
// //!              1  2  3  4  5  6  7
// //!            -------------------------
// //!   b  =    |  *  *  *  *  *  *  *   |
// //!            -------------------------
// //! ```
// //! 
// //! We can use an analogous method whenever we have the following data
// //! 
// //! - a matrix `A`
// //! - a vector `b`
// //! - a set of row indices `I` and a set of column indices `J` such that `A[I,J]` is upper triangular and invertible.
// //! 
// //! we require, in addition, that 
// //! 
// //! - if `[i,j]` is a diagonal entry of `A[I,J]`, then all entries to the left of `[i,j]` in `A` vanish.
// //!   - this applies only in cases where we we add rows of `A` to `b`
// //! - if `[i,j]` is a diagonal entry of `A[I,J]`, then all entries below `[i,j]`  in `A` vanish.
// //!   - this applies only in cases where we we add columns of `A` to `b`
// //! 
// //! In practice, when we use OAT to solve `xA = b`, we do not provide `I` and `J` explicitly;
// //! rather we note that there are (mutually inverse) bijections `f: I -> J` and `g: J -> I` such that
// //! `[i, f(i)]` and `[g(j), j]` lie on the diagonal of `A[I,J]`.  We pass these bijections to the solver, instead.
//! 
//! 
//! # <a name="solvexab">Solve xA = b for y</a>
//! 
//! OAT's echelon solvers will compute a quotient and remainder, given the following inputs
//! 
//! - a matrix `A` where no two rows share a leading nonzero entry in the same column
//! - a partially defined leading entry bijection
//! - an order operator (used to specify the order of row entries)
//! 
//! 
//!  
//! ```
//! use oat_rust::algebra::matrices::types::{vec_of_vec::sorted::VecOfVec, packet::MatrixAlgebraPacket};
//! use oat_rust::algebra::matrices::operations::solve::echelon::{RowEchelonSolver, ColumnEchelonSolverReverse};
//! use oat_rust::algebra::matrices::operations::multiply::{multiply_row_vector_with_matrix};        
//! use oat_rust::algebra::rings::types::field_prime_order::BooleanField;
//! use oat_rust::algebra::vectors::operations::VectorOperations;
//! use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
//! use oat_rust::utilities::order::{OrderOperatorByKey, OrderOperatorByKeyReverse};
//! use itertools::Itertools;
//! use assert_panic::assert_panic;
//! 
//! // define the ring operator
//! let ring_operator = BooleanField::new();
//! 
//! // a matrix A with an invertible upper triangular submatrix
//! let matrix  =   VecOfVec::new(
//!                     vec![ 
//!                         vec![                                     ], 
//!                         vec![                                     ],                                     
//!                         vec![ (0,true),  (1,true), (2,true)       ], 
//!                         vec![            (1,true), (2,true)       ],
//!                         vec![                      (2,true)       ],                                       
//!                     ],
//!                 ).ok().unwrap();    
//! let matrix  =   MatrixAlgebraPacket::with_default_order_and_boolean_coefficients( &matrix );
//! 
//! // partial bijection from the column indices of the diagonal elements to their row indices
//! let partial_bijection = |x: usize| { if x < 3 { Some(x+2) } else { None } };
//! 
//! // the problem vector b
//! let b = vec![ (0,true), (1,true) ];
//! 
//! // solver; attempts to solve xA = b for x
//! let division =  RowEchelonSolver::solve(
//!                                                                 b.clone(), // entries must appear in strictly ASCENDING order
//!                                                                 & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
//!                                                                 EvaluateFunctionFnMutWrapper::new( partial_bijection ), // maps *matched* column indices to matched row indices
//!                                                                 BooleanField::new(), // defines the ring operations
//!                                                                 OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
//!                                                             );
//! 
//! // multiply the solution with A, and check that the product equals b                               
//! let product     =   division
//!                         .solution()
//!                         .unwrap()
//!                         .multiply_self_as_a_row_vector_with_matrix( &matrix );
//! assert_eq!( product.collect_vec(), b );   // check the solution, i.e. check that xA = b  
//! ```
//! 
//! # Alternate formulations
//! 
//! There are many ways to specify the partially defined leading entry bijection.
//! You can use any object that implements the right version of the [EvaluateFunction] trait.  For example,
//! you can replace `let partial_bijection = |x: usize| { if x < 3 { Some(x+2) } else { None } };` in the preceding
//! example with a `Vec`:
//! 
//! ```
//! # use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
//! # use oat_rust::algebra::matrices::types::packet::MatrixAlgebraPacket;
//! # use oat_rust::algebra::matrices::operations::solve::echelon::{RowEchelonSolver, ColumnEchelonSolverReverse};
//! # use oat_rust::algebra::matrices::operations::multiply::{multiply_row_vector_with_matrix};        
//! # use oat_rust::algebra::rings::types::field_prime_order::BooleanField;
//! # use oat_rust::algebra::vectors::operations::VectorOperations;
//! # use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
//! # use oat_rust::utilities::order::{OrderOperatorByKey, OrderOperatorByKeyReverse};
//! # use std::collections::HashMap;
//! # use itertools::Itertools;
//! # use assert_panic::assert_panic;
//! # 
//! # // a matrix A with an invertible upper triangular submatrix
//! # let matrix  =   VecOfVec::new(
//! #                         vec![ 
//! #                             vec![                                     ], 
//! #                             vec![                                     ],                                     
//! #                             vec![ (0,true),  (1,true), (2,true)       ], 
//! #                             vec![            (1,true), (2,true)       ],
//! #                             vec![                      (2,true)       ],                                       
//! #                         ],
//! #                     ).ok().unwrap(); 
//! # let matrix  =   MatrixAlgebraPacket::with_default_order_and_boolean_coefficients( &matrix ); 
//! # 
//! # // the problem vector b
//! # let b = vec![ (0,true), (1,true) ];
//! # 
//! # // the partial bijection
//! let partial_bijection = vec![2,3,4]; // or as a vector
//! # 
//! # // create a solver to solve xA = b for x
//! # let division =  RowEchelonSolver::solve(
//! #                                                                 b.clone(), // entries must appear in strictly ASCENDING order
//! #                                                                 & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
//! #                                                                 partial_bijection, // maps *matched* column indices to matched row indices
//! #                                                                 BooleanField::new(), // defines the ring operations
//! #                                                                 OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
//! #                                                             );
//! #                                                         
//! # // check the solution, i.e. check that xA = b
//! # let product = division
//!                     .solution()
//!                     .unwrap()
//!                     .multiply_self_as_a_row_vector_with_matrix( &matrix );
//! # assert_eq!( product.collect_vec(), b.clone() );  
//! ```
//! 
//! or a `HashMap`: 
//! 
//! ```
//! # use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
//! # use oat_rust::algebra::matrices::types::packet::MatrixAlgebraPacket;
//! # use oat_rust::algebra::matrices::operations::solve::echelon::{RowEchelonSolver, ColumnEchelonSolverReverse};
//! # use oat_rust::algebra::matrices::operations::multiply::{multiply_row_vector_with_matrix};        
//! # use oat_rust::algebra::rings::types::field_prime_order::BooleanField;
//! # use oat_rust::algebra::vectors::operations::VectorOperations;
//! # use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
//! # use oat_rust::utilities::order::{OrderOperatorByKey, OrderOperatorByKeyReverse};
//! # use std::collections::HashMap;
//! # use itertools::Itertools;
//! # use assert_panic::assert_panic;
//! # 
//! # // a matrix A with an invertible upper triangular submatrix
//! # let matrix  =   VecOfVec::new(
//! #                         vec![ 
//! #                             vec![                                     ], 
//! #                             vec![                                     ],                                     
//! #                             vec![ (0,true),  (1,true), (2,true)       ], 
//! #                             vec![            (1,true), (2,true)       ],
//! #                             vec![                      (2,true)       ],                                       
//! #                         ],
//! #                     ).ok().unwrap(); 
//! # let matrix  =   MatrixAlgebraPacket::with_default_order_and_boolean_coefficients( &matrix );
//! # 
//! # // the problem vector b
//! # let b = vec![ (0,true), (1,true) ];
//! # 
//! # // the partial bijection
//! let partial_bijection = HashMap::from( [(0,2),(1,3),(2,4)] ); // or as a hashmap
//! # 
//! # // create a solver to solve xA = b for x
//! # let division =  RowEchelonSolver::solve(
//! #                                                                 b.clone(), // entries must appear in strictly ASCENDING order
//! #                                                                 & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
//! #                                                                 partial_bijection, // maps *matched* column indices to matched row indices
//! #                                                                 BooleanField::new(), // defines the ring operations
//! #                                                                 OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
//! #                                                             );
//! #                                                         
//! # // check the solution, i.e. check that xA = b
//! # let product = division
//!                     .quotient()
//!                     .multiply_self_as_a_row_vector_with_matrix( & matrix );
//! # assert_eq!( product.collect_vec(), b.clone() );  
//! ```
//! 
//! 
//! 
//! # <a name="solveayc">Solve Ay = c for y</a>
//! 
//! When solving for columns we have to eliminate **trailing** entries, rather than leading entries, in order to take advantage
//! of the upper triangular structure.  Therefore the entries of `c` have to appear in **strictly descending order**.
//! 
//! ```
//! use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
//! use oat_rust::algebra::matrices::operations::solve::echelon::{RowEchelonSolver, ColumnEchelonSolverReverse};
//! use oat_rust::algebra::matrices::operations::multiply::{multiply_row_vector_with_matrix};        
//! use oat_rust::algebra::rings::types::field_prime_order::BooleanField;
//! use oat_rust::algebra::vectors::operations::VectorOperations;
//! use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
//! use oat_rust::utilities::order::{OrderOperatorByKey, OrderOperatorByKeyReverse};
//! use std::collections::HashMap;
//! use itertools::Itertools;
//! use assert_panic::assert_panic;
//! 
//! let matrix  =   VecOfVec::new(
//!                     vec![ 
//!                         vec![                                     ], 
//!                         vec![                                     ],                                     
//!                         vec![ (0,true),  (1,true), (2,true)       ], 
//!                         vec![            (1,true), (2,true)       ],
//!                         vec![                      (2,true)       ],                                       
//!                     ],
//!                 ).ok().unwrap();  
//! 
//! let c = vec![ (3,true), (2,true) ];   // define a sparse vector c                                                                                                                          
//! let division =  ColumnEchelonSolverReverse::solve(
//!                                                                 c.clone(), // entries must appear in strictly DESCENDING order
//!                                                                 & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
//!                                                                 HashMap::from( [(2,0),(3,1),(4,2)] ), // or as a hashmap, // maps *matched* row indices to matched column indices
//!                                                                 BooleanField::new(), // defines the ring operations
//!                                                                 OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever j ≤ i
//!                                                             );
//!                                                         
//! let product = division
//!                 .quotient()
//!                 .multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order_custom( 
//!                     &matrix, 
//!                     BooleanField::new(), 
//!                     OrderOperatorByKey::new() 
//!                 );
//! 
//! // solution satisfies Ay = c      
//! assert_eq!( 
//!     product.collect_vec(), 
//!     c.clone() 
//! );   
//! ```
//! 
//! # <a name="qrsolutions">Quotients, remainders, and solutions</a>
//! 
//! The solvers return a [QuotientRemainderSolver] object, that can be used to extract quotients, remainders, and solutions.
//! 
//! **Fact** a quotient is a solution iff the corresponding remainder is zero.
//! 
//! Calling `.solution()` on a [QuotientRemainderSolver] object will therefore return `Some(quotient)` if the remainer is zero, and `None`
//! otherwise.
//! 
//! 
//! #### Case 1: There exists a solution
//! 
//! ```
//! use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
//! use oat_rust::algebra::matrices::operations::solve::echelon::{RowEchelonSolver, ColumnEchelonSolverReverse};
//! use oat_rust::algebra::matrices::operations::multiply::{multiply_row_vector_with_matrix, };        
//! use oat_rust::algebra::rings::types::field_prime_order::BooleanField;
//! use oat_rust::algebra::vectors::operations::VectorOperations;
//! use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
//! use oat_rust::utilities::order::{OrderOperatorByKey, OrderOperatorByKeyReverse};
//! use std::collections::HashMap;
//! use itertools::Itertools;
//! use assert_panic::assert_panic;
//! 
//! // define the ring operator
//! let ring_operator = BooleanField::new();
//! 
//! // a matrix A with an invertible upper triangular submatrix
//! let matrix  =   VecOfVec::new(
//!                     vec![ 
//!                         vec![                                     ], 
//!                         vec![                                     ],                                     
//!                         vec![ (0,true),  (1,true), (2,true)       ], 
//!                         vec![            (1,true), (2,true)       ],
//!                         vec![                      (2,true)       ],                                       
//!                     ],
//!                 ).ok().unwrap();  
//! 
//! // create a solver to solve xA = b for x
//! let division =  RowEchelonSolver::solve(
//!                                                                 vec![ (0,true), (1,true) ], // entries must appear in strictly ASCENDING order
//!                                                                 & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
//!                                                                 vec![2,3,4], // maps *matched* column indices to matched row indices
//!                                                                 BooleanField::new(), // defines the ring operations
//!                                                                 OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
//!                                                             );
//! 
//! // the remainder is zero
//! assert_eq!( 
//!     division.clone().remainder().collect_vec(),
//!     vec![ ]
//! );
//!
//! // therefore the solution equals the quotient
//! assert_eq!( 
//!     division.clone().solution(), 
//!     Some( division.clone().quotient().collect_vec() ) 
//! ); 
//! 
//! // we can also split the division object into a separate quotient and remainder
//! let qr = division.clone().quotient_remainder();
//! assert_eq!( 
//!     qr.remainder.clone().collect_vec(), 
//!     vec![] 
//! );
//! assert_eq!( 
//!     division.clone().solution(), 
//!     Some( qr.quotient.clone() ),
//! );
//! ```
//! 
//! 
//! #### Case 2: There exists no solution
//! 
//! ```
//! use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
//! use oat_rust::algebra::matrices::operations::solve::echelon::{RowEchelonSolver, ColumnEchelonSolverReverse};
//! use oat_rust::algebra::matrices::operations::multiply::{multiply_row_vector_with_matrix};        
//! use oat_rust::algebra::rings::types::field_prime_order::BooleanField;
//! use oat_rust::algebra::vectors::operations::VectorOperations;
//! use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
//! use oat_rust::utilities::order::{OrderOperatorByKey, OrderOperatorByKeyReverse};
//! use std::collections::HashMap;
//! use itertools::Itertools;
//! use assert_panic::assert_panic;
//! 
//! // define the ring operator
//! let ring_operator = BooleanField::new();
//! 
//! // a matrix A with an invertible upper triangular submatrix
//! let matrix  =   VecOfVec::new(
//!                     vec![ 
//!                         vec![                                     ], 
//!                         vec![                                     ],                                     
//!                         vec![ (0,true),  (1,true), (2,true)       ], 
//!                         vec![            (1,true), (2,true)       ],
//!                         vec![                      (2,true)       ],                                       
//!                     ],
//!                 ).ok().unwrap();  
//! 
//! // create a solver to solve xA = b for x
//! let division =  RowEchelonSolver::solve(
//!                                                                 vec![ (0,true), (1,true), (2,true), (3,true) ], // entries must appear in strictly ASCENDING order
//!                                                                 & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
//!                                                                 vec![2,3,4], // maps *matched* column indices to matched row indices
//!                                                                 BooleanField::new(), // defines the ring operations
//!                                                                 OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
//!                                                             );
//! 
//! // the remainder is nonzero
//! assert_eq!( 
//!     division.clone().remainder().collect_vec(),
//!     vec![ (3,true) ]
//! );
//!
//! // therefore the solution equals None
//! assert_eq!( 
//!     division.clone().solution(), 
//!     None
//! );                                     
//! ```
//! 
//! 
//! # Panics
//! 
//! 
//! The solver will panic if we pass an unsorted input, and pull all entries from the quotient and remainder
//! 
//! ```
//! use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
//! use oat_rust::algebra::matrices::operations::solve::echelon::{RowEchelonSolver, ColumnEchelonSolverReverse};
//! use oat_rust::algebra::matrices::operations::multiply::{multiply_row_vector_with_matrix};        
//! use oat_rust::algebra::rings::types::field_prime_order::BooleanField;
//! use oat_rust::algebra::vectors::operations::VectorOperations;
//! use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
//! use oat_rust::utilities::order::{OrderOperatorByKey, OrderOperatorByKeyReverse};
//! use itertools::Itertools;
//! use assert_panic::assert_panic;
//! 
//! // define the ring operator
//! let ring_operator = BooleanField::new();
//! 
//! // a matrix A with an invertible upper triangular submatrix
//! let matrix  =   VecOfVec::new(
//!                     vec![ 
//!                         vec![                                     ], 
//!                         vec![                                     ],                                     
//!                         vec![ (0,true),  (1,true), (2,true)       ], 
//!                         vec![            (1,true), (2,true)       ],
//!                         vec![                      (2,true)       ],                                       
//!                     ],
//!                 ).ok().unwrap();  
//! 
//! let d = vec![ (1,true), (0,true), (3,true), (2,true), ];
//! let division =  RowEchelonSolver::solve(
//!                                                                 d.clone(), // entries must appear in strictly ASCENDING order
//!                                                                 & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
//!                                                                 vec![2,3,4], // maps *matched* column indices to matched row indices
//!                                                                 BooleanField::new(), // defines the ring operations
//!                                                                 OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
//!                                                             );                                                                
//! assert_panic!( for _ in division.solver() {} );   // the elimination procedure clears the first three columns, but cannot clearn the fourth
//! ```
//! 
//! The solver *may possibly* compute a quotient without throwing an error,
//! if it finds an entry it can't eliminate before hitting the first consecutive out-of-order pair
//! 
//! 
//! ```
//! use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
//! use oat_rust::algebra::matrices::operations::solve::echelon::{RowEchelonSolver, ColumnEchelonSolverReverse};
//! use oat_rust::algebra::matrices::operations::multiply::{multiply_row_vector_with_matrix};        
//! use oat_rust::algebra::rings::types::field_prime_order::BooleanField;
//! use oat_rust::algebra::vectors::operations::VectorOperations;
//! use oat_rust::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
//! use oat_rust::utilities::order::{OrderOperatorByKey, OrderOperatorByKeyReverse};
//! use itertools::Itertools;
//! use assert_panic::assert_panic;
//! 
//! // define the ring operator
//! let ring_operator = BooleanField::new();
//! 
//! // a matrix A with an invertible upper triangular submatrix
//! let matrix  =   VecOfVec::new(
//!                     vec![ 
//!                         vec![                                     ], 
//!                         vec![                                     ],                                     
//!                         vec![ (0,true),  (1,true), (2,true)       ], 
//!                         vec![            (1,true), (2,true)       ],
//!                         vec![                      (2,true)       ],                                       
//!                     ],
//!                 ).ok().unwrap();  
//! 
//! 
//! let d = vec![ (0,true), (1,true), (2,true), (3,true), (4,true), (2,true), ];
//! let division =  RowEchelonSolver::solve(
//!                                                                 d.clone(), // entries must appear in strictly ASCENDING order
//!                                                                 & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
//!                                                                 vec![2,3,4], // maps *matched* column indices to matched row indices
//!                                                                 BooleanField::new(), // defines the ring operations
//!                                                                 OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
//!                                                             );
//! 
//! let mut qr = division.quotient_remainder();
//! 
//! // no panic when computing the quotient
//! assert_eq!( 
//!     qr.quotient.clone(), 
//!     vec![ (2,true) ] 
//! ); 
//! 
//! // no panic when computing the first element of the remainder 
//! assert_eq!( 
//!     qr.remainder.next(), 
//!     Some( (3,true) ) 
//! ); 
//! 
//! // panics once we exhaust the problem vector
//! assert_panic!( 
//!     for _ in qr.remainder {} 
//! );   
//! ```



//  ===========================================================================================
//  ===========================================================================================
//  ===========================================================================================
//  ===========================================================================================



use derive_getters::Dissolve;
use derive_new::new;


use crate::algebra::matrices::types::transpose::OrderAntiTranspose;
use crate::algebra::matrices::query::{MatrixOracle};
use crate::algebra::rings::traits::DivisionRingOperations;
use crate::algebra::vectors::entries::{KeyValSet, KeyValGet};
use crate::algebra::vectors::operations::{VectorOperations, LinearCombinationSimplified, Scale, Simplify};

use crate::utilities::iterators::merge::hit::{hit_bulk_insert, hit_merge_by_predicate};
use crate::utilities::order::{JudgePartialOrder, ReverseOrder};
use crate::utilities::functions::evaluate::EvaluateFunction;
use crate::utilities::iterators::general::{HeadTail, TwoTypeIterator, RequireStrictAscentWithPanic, TransformIter};
use crate::utilities::iterators::merge::hit::IteratorsMergedInSortedOrder;
// use debugit::DebugIt as D;


/// Standardizes the method of extracting a remainder from a solver
/// 
/// See the documentaiton for [QuotientRemainderSolver], for an explanation.
pub trait IntoRemainder{
    type Remainder;
    fn into_remainder( self ) -> Self::Remainder;
}

/// Wrapper for the quotient and remainder of a division b/A
/// 
/// See the documentaiton for [QuotientRemainderSolver], for an explanation of terms.
#[derive(Clone,Copy,Debug)]
pub struct QuotientRemainderOutput< Q, R >{
    pub quotient: Q,
    pub remainder: R,
}

// /// Splits a solver into a quotient and a solution
// /// 
// /// See [QuotientRemainderOutput] for the full interpretation.
// pub trait IntoQuotientRemainderOutput< Remainder >
//     where 
//         Self:       Iterator + IntoRemainder,
// {
//     /// Splits the solution into a quotient and a remainder
//     /// 
//     /// Operates by draining the iterator 
//     /// 
//     /// See [QuotientRemainderOutput] for the full interpretation.
//     fn into_quotient_remainder( self ) -> QuotientRemainderOutput< Vec< Self::Item >, Remainder > {
//         // let mut q = Vec::new();
//         // while let Some(i) = self.next() { q.push(i) }
//         // q.shrink_to_fit();
//         // return (q, self)
//     }
// }


/// Safe wrapper for iterators that solve linear systems, with quotient and remainder parts
/// 
/// Echelon solvers use an incremental approach that involves iterative
/// elimination of leading entries of a vector `r` by adding scalar multiples of rows or columns from a matrix `A`.  The solver terminates when either
/// 
/// - no nonzero entries remain, OR
/// - the solver cannot eleiminate a leading entry
/// 
/// We call the linear coefficients used to combine rows/columns of `A` the **quotient**, and we call whatever remains after the elimination the **remainder**
/// 
/// A solution exists iff the remainder equals zero; in this case the quotient is a solution.
/// 
/// This struct provides methods for extracting quotients, remainders, and solutions from a solver.

#[derive(Copy,Debug,new)]
pub struct QuotientRemainderSolver< I >
    where I:    Iterator + IntoRemainder,
{
    solver:  I
}

impl < I > QuotientRemainderSolver< I >
    where 
        I:              Iterator + IntoRemainder,
        I::Remainder:   Iterator,

{
    /// Splits apart the quotient and remainder of an echelon solver
    /// 
    /// The quotient is a bona fide solution to `Ax = b` or `xA = b` iff the remainder is empty.
    /// 
    /// See [QuotientRemainderSolver] for details.
    pub fn quotient_remainder( mut self ) -> QuotientRemainderOutput< Vec<I::Item>, I::Remainder > {
        let mut quotient = Vec::new();
        for i in self.solver.by_ref() { quotient.push(i) }
        quotient.shrink_to_fit();
        let remainder = self.solver.into_remainder();
        QuotientRemainderOutput{ quotient, remainder }       
    }

    /// Returns the remainder
    /// 
    /// See [QuotientRemainderSolver] for details.
    pub fn remainder( mut self ) -> I::Remainder {
        for _ in self.solver.by_ref() {} // drain the quotient
        self.solver.into_remainder() // what remains is the remainder
    }    

    /// Returns the quotient, which is a solution iff the remainder is zero
    /// 
    /// See [QuotientRemainderSolver] for furhter details.
    /// 
    /// This is equivalent to calling `self.solver()`, because the solver iterates over elements of the quotient.
    pub fn quotient( self ) -> I { self.solver }

    /// Returns the solver, which iterates over elements of the quotient
    /// 
    /// See [QuotientRemainderSolver] for furhter details.
    pub fn solver( self ) -> I { self.solver }    

    /// Returns a solution formated as `Vec<T>`, if a solution exists
    /// 
    /// Concretely, returns the quotient of `b/A` if the remainder is zero, and returns `None` otherwise.
    /// 
    /// Alternate methods
    /// 
    /// - [QuotientRemainderSolver::quotient] provides an iterator for the quotient; this can be faster and more memory efficient, however **the quotient is only a bona fide solution when the remainder is zero**
    /// - [QuotientRemainderSolver::quotient_remainder], provides the same information, plus information about the remainder.  However, take care to remember that the quotient returned by that method is only a bona fide solution if the remainder is empty.
    /// 
    /// See [QuotientRemainderSolver] for furhter details.
    pub fn solution( self ) -> Option< Vec< I::Item > > {   
        let mut qr = self.quotient_remainder();
        match qr.remainder.next().is_none() {
            true    =>  Some( qr.quotient ),
            false   =>  None,
        }
    }    
}


impl < I >

    Clone for 
    
    QuotientRemainderSolver< I >

    where I:    Clone + Iterator + IntoRemainder,
{
    fn clone(&self) -> Self {
        QuotientRemainderSolver { solver: self.solver.clone() }
    }
}


//  ======================================================================================
//  ======================================================================================
//  ======================================================================================
//  ASCENDING MAJOR SOLVE
//  ======================================================================================
//  ======================================================================================
//  ======================================================================================



//  ======================================================================================
//  ECHELON SOLVE ASCEND WITH MINOR KEYS
//  ======================================================================================




// ECHELON MAJOR ASCEND SOLVE W/ MINOR KEYS (ITERATOR)
// ---------------------------------------------------------------------------

/// Iterates over the entries of the solution `x` to a matrix equation `xA = b`, where
/// (i) each entry of `x` representing a pair of form `(row_index, coefficient)` is replaced by an entry representing `(match(row_index), coefficient)`, where
/// `match(row_index)` is the matched column index, and
/// (ii) entries appear in ascending order, according to column index
/// 
/// # Assumptions
/// 
/// - `b` is a sparse vector iterator indexed by the column indices of `A`, and the entries of `b` appear in ascending order of index
/// - There is a partial bijection `match: column_index -> row_index` from the column indices of `A` to the row indices of `A`; concretely, this means a bijection from a subset of the column indices to a subset of the row indices
/// - Whenever `row_index = match( column_index )`, the leading entry of `A.row( &row_index )` has index `column_index`
/// - A solution to `xA = b` exists
/// 
/// # What happens when an assumption is violated
/// 
/// This iterator will generate entries in the manner described below (under heading "Implementation"), whether or not there exists a solution to `Ax = b`.
/// After the iterator is exhausted, the user can call `x.remainder()`; the resulting iterator, `y`, represents what is "left over" from the vector `b` after some of its entries are eliminated by adding
/// scalar multiples of rows of `A`.  If `y` is empty, then a solution exists; in particular, `x` is a feasible solution.  On the other hand, if `y` is
/// nonempty then no solution exists, and `x` can be regarded merely as a failed attempt.
/// 
/// # Implementation
/// 
/// The iterator proceeds by repeating the following steps: (i) if `b = 0` has no leading entry, return `None`; we have found a valid solution, `x`, (ii) otherwise `b` has a 
/// leading (nonzero) entry.  Let `column_index` be the index of this entry.  (iii) if `column_index` is unmatched, return `None`; there exists no solution.  (iv) otherwise, let `row_index` be the row index matched to `column_index`.  Add a scalar multiple of `A.row( &row_index )` to `b` to 
/// eliminate its leading entry.  Let `alpha` denote the scalar by which we multiply, and return `Some( (column_index, -alpha) )`.
#[derive(Dissolve)]
pub struct RowEchelonSolverReindexed< 
                    ProblemVector,
                    MatchingFromColumnIndicesToRowIndices,
                    Matrix,                   
                    RingOperator,
                    OrderOperator,
                > 
    where 
        ProblemVector:                              IntoIterator< Item = Matrix::RowEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromColumnIndicesToRowIndices:      EvaluateFunction< Matrix::ColumnIndex, Option< Matrix::RowIndex > >,
        Matrix:                                     MatrixOracle,
        Matrix::RowEntry:                           KeyValSet,      
        RingOperator:                               Clone + DivisionRingOperations< Element = Matrix::Coefficient >,
        OrderOperator:                              Clone + JudgePartialOrder <  Matrix::RowEntry >,                   

{   
    ring_operator:                      RingOperator, // the operator for the coefficient ring_operator    
    matrix:                             Matrix, // the `A` in `Ax = b`
    entries_to_eliminate:               HeadTail<
                                                LinearCombinationSimplified<   
                                                        TwoTypeIterator< // this enum allows us to treat two different iterator types as a single type
                                                                Matrix::Row,  // the iterators returned by the matrix                                                                
                                                                RequireStrictAscentWithPanic< ProblemVector::IntoIter, OrderOperator >,
                                                            >,
                                                        RingOperator,
                                                        OrderOperator,
                                                    >
                                            >,
    matching_from_column_index_to_row_index:     MatchingFromColumnIndicesToRowIndices,
    // phantom_lifetime:                   PhantomData< Matrix::RowIndex > 
} 




impl    < 
            ProblemVector,
            MatchingFromColumnIndicesToRowIndices,
            Matrix,                   
            RingOperator,
            OrderOperator,
        >  

        RowEchelonSolverReindexed<
                ProblemVector,
                MatchingFromColumnIndicesToRowIndices,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                  IntoIterator< Item = Matrix::RowEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromColumnIndicesToRowIndices:     EvaluateFunction< Matrix::ColumnIndex, Option< Matrix::RowIndex > >,
        Matrix:                         MatrixOracle,
        Matrix::RowEntry:               KeyValSet< Key=Matrix::ColumnIndex, Val=Matrix::Coefficient >,      
        RingOperator:                   Clone + DivisionRingOperations< Element = Matrix::Coefficient >,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::RowEntry >,  

{

    /// Attempts to solve `xA = b` for `x`; see [`EchelonSolverMajorAscendWithMinorKeys`].
    pub fn solve(
                b:                              ProblemVector,                
                a:                              Matrix,
                matching_from_column_index_to_row_index: MatchingFromColumnIndicesToRowIndices,
                ring_operator:                  RingOperator,
                order_operator:                 OrderOperator,
            )
            ->
            QuotientRemainderSolver<
                    RowEchelonSolverReindexed<
                            ProblemVector,
                            MatchingFromColumnIndicesToRowIndices,
                            Matrix,                    
                            RingOperator,
                            OrderOperator,
                        >
                >
    {
        // we will eliminate entries by adding other vectors; addition is performed by merging in new iterators;
        // because the vectors we merge in may have a different type than the problem vector b, we need to do a
        // "type union" by wrapping `b` inside an `TwoTypeIterator`
        let entries_to_eliminate
                =   hit_merge_by_predicate( 
                            std::iter::once(    
                                    TwoTypeIterator::Version2( b.into_iter().require_strict_ascent_with_panic( order_operator.clone() ) )   // wrap b in an TwoTypeIterator enum so that it can be easily combined with  
                                        .scale_by( ring_operator.minus_one(), ring_operator.clone() )
                                ),
                            order_operator,
                        )
                        .simplify( ring_operator.clone() );
        // load the entries in the tail of a HeadTail
        let entries_to_eliminate = HeadTail { head: None, tail: entries_to_eliminate };
        QuotientRemainderSolver::new(
            RowEchelonSolverReindexed{
                    ring_operator,
                    matching_from_column_index_to_row_index,
                    matrix:                             a,
                    entries_to_eliminate, 
                    // phantom_lifetime:                   PhantomData,
                }
        )
    }
  
}



impl    < 
            ProblemVector,
            MatchingFromColumnIndicesToRowIndices,
            Matrix,                   
            RingOperator,
            OrderOperator,
        >  

        IntoRemainder for

        RowEchelonSolverReindexed<
                ProblemVector,
                MatchingFromColumnIndicesToRowIndices,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                  IntoIterator< Item = Matrix::RowEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromColumnIndicesToRowIndices:     EvaluateFunction< Matrix::ColumnIndex, Option< Matrix::RowIndex > >,
        Matrix:                         MatrixOracle,
        Matrix::RowEntry:               KeyValSet< Key=Matrix::ColumnIndex, Val=Matrix::Coefficient >,      
        RingOperator:                   Clone + DivisionRingOperations< Element = Matrix::Coefficient >,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::RowEntry >,  


{
    type Remainder = 
            HeadTail<
                    LinearCombinationSimplified<   
                            TwoTypeIterator< // this enum allows us to treat two different iterator types as a single type
                                    < Matrix::Row as IntoIterator >::IntoIter,  // the iterators returned by the matrix                                                                
                                    RequireStrictAscentWithPanic< ProblemVector::IntoIter, OrderOperator >, // the vector we have tried to eliminate
                                >,
                            RingOperator, OrderOperator 
                        >        
                >;


        /// Returns an iterator that runs over all the entries that have not yet been eliminated, consuming the solver
        /// in the process.
        /// 
        /// Recall that the solver works by iteratively adding linear multiples rows or columns of matrix `A` to the
        /// problem vector `b`, in order to solve `xA=b`.  If `v_k` is the vector we add on iteration `k`, then the remainder after iteration `k` is
        /// `b - v_0 - v_1 ... - v_k`.  The leading entry of the remainder is the entry we must try to eliminate on 
        /// iteration `v_{k+1}`.  Note that the remainder can be nonzero even when `self.next() = None`.  This happens
        /// if and only if there exists no solution to `xA = b`.
    fn into_remainder( self ) -> Self::Remainder {
        self.entries_to_eliminate
    }

}




impl    <
            ProblemVector,
            MatchingFromColumnIndicesToRowIndices,
            Matrix,                   
            RingOperator,
            OrderOperator,
        > 

        Iterator for    

        RowEchelonSolverReindexed<
                ProblemVector,
                MatchingFromColumnIndicesToRowIndices,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                  IntoIterator< Item = Matrix::RowEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromColumnIndicesToRowIndices:     EvaluateFunction< Matrix::ColumnIndex, Option< Matrix::RowIndex > >,
        Matrix:                         MatrixOracle,
        Matrix::RowEntry:               KeyValSet< Key=Matrix::ColumnIndex, Val=Matrix::Coefficient >,      
        RingOperator:                   Clone + DivisionRingOperations< Element = Matrix::Coefficient >,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::RowEntry >, 

{
    type Item = Matrix::RowEntry;

    fn next( &mut self ) -> Option<Self::Item> { 

        // THIS ITERATOR PROCEEDS BY ELIMINATING ENTRIES IN `self.entries_to_eliminate`

        //  IF THIS ITERATOR IS NONEMPTY THEN WE WILL GENERATE THE NEXT ENTRY OF THE SOLUTION VECTOR; 
        //  note that `entry_to_eliminate` will always have a nonzero coefficient, because `self.entries_to_eliminate` is a struct of type `Simplify< .. >`; structs of this type only return nonzero entries.
        if let Some( mut entry_to_eliminate )       =   self.entries_to_eliminate.next() {

            // println!("EchelonSolverMajorAscendWithMinorKeys.next(): BRANCH TO ELIMINATE ENTRY: STARTING");

            // KEY TO LOOK UP THE VECTOR THAT WILL EXECUTE THE ELIMINATION
            // let row_index_of_eliminating_row = self.matching_from_column_index_to_row_index.evaluate_function( entry_to_eliminate.key() );

            match self.matching_from_column_index_to_row_index.evaluate_function( entry_to_eliminate.key() ) {
                None =>  {
                    // put back the entry
                    self.entries_to_eliminate.head = Some( entry_to_eliminate );
                    None
                } Some( row_index_of_eliminating_row ) => {
                    // println!("EchelonSolverMajorAscendWithMinorKeys.next(): BRANCH TO ELIMINATE ENTRY: WAYPOINT 1");            

                    // TO ELIMINATE `entry_to_eliminate`, WE ADD A SCALAR MULTIPLE OF A VECTOR TAKEN FROM THE MATRIX; DENOTE THE UNSCLAED VECTOR BY v
                    let mut seed_of_eliminating_iterator    =   self.matrix.row( &row_index_of_eliminating_row ).into_iter();
                    
                    // println!("EchelonSolverMajorAscendWithMinorKeys.next(): BRANCH TO ELIMINATE ENTRY: WAYPOINT 2");            

                    // POP THE LEADING ENTRY OUT OF THE UNSCALED VECTOR v; USE THIS ENTRY TO DETERMINE THE SCALE FACTOR BY WHICH WE'LL MULTIPLY THE REST OF v
                    let eliminating_entry           =   seed_of_eliminating_iterator.next().unwrap();

                    // THE SCALE FACTOR IS EQUAL TO -(entry_to_eliminate)/(eliminating_entry)
                    let scale_factor     =    self.ring_operator.negate(
                                                            self.ring_operator.divide( entry_to_eliminate.val(), eliminating_entry.val() )
                                                        );
                    
                    // SCALE v SO THAT ITS LEADING ENTRY WILL CANCEL `entry_to_eliminate`
                    // note that it won't *actually* cancel that entry in practice, because we already popped the leading entry off of v and off of entries_to_eliminate
                    let eliminating_iterator       
                            =   TwoTypeIterator::Version1( seed_of_eliminating_iterator ) 
                                    .scale_by( scale_factor.clone(), self.ring_operator.clone() );

                    // MERGE THE (BEHEADED) SCALAR MULTIPLE OF v INTO THE HEAP, NAMELY `entries_to_eliminate`
                    // println!("EchelonSolverMajorAscendWithMinorKeys.next(): hit_bulk_insert: STARTING");
                    hit_bulk_insert( &mut self.entries_to_eliminate.tail.unsimplified, vec![eliminating_iterator] ); 
                    // println!("EchelonSolverMajorAscendWithMinorKeys.next(): hit_bulk_insert: ENDING");

                    // NOW RE-USE `entry_to_eliminate` AS A CONTAINER FOR THE NEXT ENTRY IN THE INVERSE MATRIX
                    entry_to_eliminate.set_val( scale_factor );

                    Some( entry_to_eliminate ) // recall that we have re-purposed `entry_to_eliminate` to hold the output of this function
                }
            }

        } else {           
            None
        }
    }
}  



impl    < 
            ProblemVector,
            MatchingFromColumnIndicesToRowIndices,
            Matrix,                   
            RingOperator,
            OrderOperator,
        > 

        Clone for    

        RowEchelonSolverReindexed<
                ProblemVector,
                MatchingFromColumnIndicesToRowIndices,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                              IntoIterator< 
                                                        Item = Matrix::RowEntry, // the iterator runs over entries of the same type as the matrix
                                                        IntoIter: Clone 
                                                    >, 
        MatchingFromColumnIndicesToRowIndices:      Clone + EvaluateFunction< Matrix::ColumnIndex, Option< Matrix::RowIndex > >,
        Matrix:                                     Clone + MatrixOracle< Row: Clone >,
        Matrix::RowEntry:                           Clone + KeyValSet,      
        RingOperator:                               Clone + DivisionRingOperations< Element = Matrix::Coefficient >,
        OrderOperator:                              Clone + JudgePartialOrder <  Matrix::RowEntry >,             
{
    fn clone(&self) -> Self {
        RowEchelonSolverReindexed{ 
                ring_operator: self.ring_operator.clone(),
                matrix: self.matrix.clone(),
                entries_to_eliminate: self.entries_to_eliminate.clone(),
                matching_from_column_index_to_row_index: self.matching_from_column_index_to_row_index.clone(),
            }
    }
}



















//  ======================================================================================
//  ECHELON MAJOR ASCEND SOLVE ASCEND WITH ROW INDEXS
//  ======================================================================================




// ECHELON SOLVE W/ ROW INDEXS (ITERATOR)
// ---------------------------------------------------------------------------

/// Iterates over the entries of the solution `x` to a matrix equation `xA = b`, where entries appear in ascending order, according to row index
/// 
/// # Assumptions
/// 
/// - `b` is a sparse vector iterator indexed by the column indices of `A`, and the entries of `b` appear in ascending order of index
/// - There is a partial bijection `match: column_index -> row_index` from the column indices of `A` to the row indices of `A`; concretely, this means a bijection from a subset of the column indices to a subset of the row indices
/// - Whenever `row_index = match( column_index )`, the leading entry of `A.row( &row_index )` has index `column_index`
/// - A solution to `xA = b` exists
/// 
/// # What happens when an assumption is violated
/// 
/// This iterator will generate entries in the manner described below (under heading "Implementation"), whether or not there exists a solution to `Ax = b`.
/// After the iterator is exhausted, the user can call `x.remainder()`; the resulting iterator, `y`, represents what is "left over" from the vector `b` after some of its entries are eliminated by adding
/// scalar multiples of rows of `A`.  If `y` is empty, then a solution exists; in particular, `x` is a feasible solution.  On the other hand, if `y` is
/// nonempty then no solution exists, and `x` can be regarded merely as a failed attempt.
/// 
/// # Implementation
/// 
/// The iterator proceeds by repeating the following steps: (i) if `b = 0` has no leading entry, return `None`; we have found a valid solution, `x`, (ii) otherwise `b` has a 
/// leading (nonzero) entry.  Let `column_index` be the index of this entry.  (iii) if `column_index` is unmatched, return `None`; there exists no solution.  (iv) otherwise, let `row_index` be the row index matched to `column_index`.  Add a scalar multiple of `A.row( &row_index )` to `b` to 
/// eliminate its leading entry.  Let `alpha` denote the scalar by which we multiply, and return `Some( (row_index, -alpha) )`.
#[derive(Dissolve)]
pub struct RowEchelonSolver<
                    ProblemVector,
                    MatchingFromColumnIndicesToRowIndices,
                    Matrix,                    
                    RingOperator,
                    OrderOperator,
                > 
    where 
        ProblemVector:                              IntoIterator< Item = Matrix::RowEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromColumnIndicesToRowIndices:      EvaluateFunction< Matrix::ColumnIndex, Option< Matrix::RowIndex > >,
        Matrix:                                     MatrixOracle,
        Matrix::RowEntry:                           KeyValSet< Key=Matrix::ColumnIndex, Val=Matrix::Coefficient >,      
        RingOperator:                               DivisionRingOperations< Element = Matrix::Coefficient >,
        OrderOperator:                              JudgePartialOrder <  Matrix::RowEntry >,                   

{   
    ring_operator:                      RingOperator, // the operator for the coefficient ring_operator    
    matrix:                             Matrix, // the `A` in `Ax = b`
    entries_to_eliminate:               HeadTail<
                                                LinearCombinationSimplified<   
                                                        TwoTypeIterator< // this enum allows us to treat two different iterator types as a single type
                                                                < Matrix::Row as IntoIterator >::IntoIter,  // the iterators returned by the matrix                                                                
                                                                RequireStrictAscentWithPanic< ProblemVector::IntoIter, OrderOperator >,
                                                                // Matrix::RowEntry,
                                                            >,
                                                        RingOperator, 
                                                        OrderOperator,
                                                    >
                                            >,
    matching_from_column_index_to_row_index:     MatchingFromColumnIndicesToRowIndices,
} 


impl    < 
            ProblemVector,
            MatchingFromColumnIndicesToRowIndices,
            Matrix,                   
            RingOperator,
            OrderOperator,
        >  

        RowEchelonSolver<
                ProblemVector,
                MatchingFromColumnIndicesToRowIndices,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                              IntoIterator< Item = Matrix::RowEntry >, // the iterator runs over entries of the same type as the matrix        
        MatchingFromColumnIndicesToRowIndices:      EvaluateFunction< Matrix::ColumnIndex, Option< Matrix::RowIndex > >,        
        Matrix:                                     MatrixOracle,
        Matrix::ColumnIndex:                        Clone + PartialEq,
        Matrix::RowIndex:                           Clone,
        RingOperator:                               Clone + DivisionRingOperations< Element = Matrix::Coefficient >,
        OrderOperator:                              Clone + JudgePartialOrder <  Matrix::RowEntry >,               
        Matrix::RowEntry:                           KeyValSet,  

{

    /// Attempts to solve `xA = b` for `x`; see [`RowEchelonSolver`].
    pub fn solve(
                b:                                          ProblemVector,                
                a:                                          Matrix,
                matching_from_column_index_to_row_index:    MatchingFromColumnIndicesToRowIndices,
                ring_operator:                              RingOperator,
                order_operator:                             OrderOperator,
            )
            ->
            QuotientRemainderSolver<
                    RowEchelonSolver<
                            ProblemVector,
                            MatchingFromColumnIndicesToRowIndices,
                            Matrix,                    
                            RingOperator,
                            OrderOperator,
                        >  
                >
    {
        // we will eliminate entries by adding other vectors; addition is performed by merging in new iterators;
        // because the vectors we merge in may have a different type than the problem vector b, we need to do a
        // "type union" by wrapping `b` inside an `TwoTypeIterator`
        let entries_to_eliminate
                =   hit_merge_by_predicate( 
                            std::iter::once(    
                                    TwoTypeIterator::Version2( b.into_iter().require_strict_ascent_with_panic(order_operator.clone()) )   // wrap b in an TwoTypeIterator enum so that it can be easily combined with  
                                        .scale_by( ring_operator.minus_one(), ring_operator.clone() ) 
                                ),
                            order_operator,
                        )
                        .simplify( ring_operator.clone() );
        // load the entries in the tail of a HeadTail
        let entries_to_eliminate = HeadTail { head: None, tail: entries_to_eliminate }; 
        QuotientRemainderSolver::new(
            RowEchelonSolver{
                ring_operator,
                matching_from_column_index_to_row_index,
                matrix:                             a,
                entries_to_eliminate, 
            } 
        )
    }
}





impl    < 
            ProblemVector,
            MatchingFromColumnIndicesToRowIndices,
            Matrix,                   
            RingOperator,
            OrderOperator,
        >  

        IntoRemainder for

        RowEchelonSolver<
                ProblemVector,
                MatchingFromColumnIndicesToRowIndices,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                  IntoIterator< Item = Matrix::RowEntry >, // the iterator runs over entries of the same type as the matrix        
        MatchingFromColumnIndicesToRowIndices:     EvaluateFunction< Matrix::ColumnIndex, Option< Matrix::RowIndex > >,        
        Matrix:                         MatrixOracle,
        Matrix::ColumnIndex:                         Clone + PartialEq,
        RingOperator:                   Clone + DivisionRingOperations< Element = Matrix::Coefficient >,
        Matrix::Coefficient:                         Clone,
        OrderOperator:                  Clone + JudgePartialOrder <  Matrix::RowEntry >,               
        Matrix::RowEntry:               Clone + KeyValSet< Key=Matrix::ColumnIndex, Val=Matrix::Coefficient >,  

{

    type Remainder = 
        HeadTail<
            LinearCombinationSimplified<   
                    TwoTypeIterator< // this enum allows us to treat two different iterator types as a single type
                            Matrix::Row,  // the iterators returned by the matrix                                                                
                            RequireStrictAscentWithPanic< ProblemVector::IntoIter, OrderOperator >, // the vector we have tried to eliminate
                        >,
                    RingOperator, 
                    OrderOperator,
                >
            >;


    /// Returns an iterator that runs over all the entries that have not yet been eliminated, consuming the solver
    /// in the process.
    /// 
    /// Recall that the solver works by iteratively adding linear multiples rows or columns of matrix `A` to the
    /// problem vector `b`, in order to solve `xA=b`.  If `v_k` is the vector we add on iteration `k`, then the remainder after iteration `k` is
    /// `b - v_0 - v_1 ... - v_k`.  The leading entry of the remainder is the entry we must try to eliminate on 
    /// iteration `v_{k+1}`.  Note that the remainder can be nonzero even when `self.next() = None`.  This happens
    /// if and only if there exists no solution to `xA = b`.
    fn into_remainder( self ) -> Self::Remainder
    {
        self.entries_to_eliminate
    }   
}








impl    < 
            ProblemVector,
            MatchingFromColumnIndicesToRowIndices,
            Matrix,                   
            RingOperator,
            OrderOperator,
        > 

        Iterator for    

        RowEchelonSolver<
                ProblemVector,
                MatchingFromColumnIndicesToRowIndices,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                              IntoIterator< Item = Matrix::RowEntry >, // the iterator runs over entries of the same type as the matrix        
        MatchingFromColumnIndicesToRowIndices:      EvaluateFunction< Matrix::ColumnIndex, Option< Matrix::RowIndex > >,        
        Matrix:                                     MatrixOracle,
        RingOperator:                               Clone + DivisionRingOperations< Element = Matrix::Coefficient >, // division is required to perform one of the solve operations
        OrderOperator:                              Clone + JudgePartialOrder <  Matrix::RowEntry >,               
        Matrix::RowEntry:                           KeyValSet< Key=Matrix::ColumnIndex, Val=Matrix::Coefficient >,

{
    type Item = (Matrix::RowIndex, Matrix::Coefficient);

    fn next( &mut self ) -> Option<Self::Item> { 

    // unsafe {
    //     println!("about to invoke RowEchelonSolver.nex(); problem vector =  {:?}", D(&self.entries_to_eliminate) );
    //     std::thread::sleep(std::time::Duration::from_secs(1));

    // }
    

    // THIS ITERATOR PROCEED BY ELIMINATING ENTRIES IN `self.entries_to_eliminate`

        //  IF THIS ITERATOR IS NONEMPTY THEN WE WILL GENERATE THE NEXT ENTRY OF THE SOLUTION VECTOR; 
        //  note that `entry_to_eliminate` will always have a nonzero coefficient, because `self.entries_to_eliminate` is a struct of type `Simplify< .. >`; structs of this type only return nonzero entries.
        if let Some( entry_to_eliminate )       =   self.entries_to_eliminate.next() {

            // unsafe {
            //     println!("entry to eliminate =  {:?}", D(&entry_to_eliminate) );
            // }
                            

            // INDEX OF THE VECTOR THAT WILL EXECUTE THE ELIMINATION
            // let row_index_of_eliminating_row = self.matching_from_column_index_to_row_index.evaluate_function( entry_to_eliminate.key() );

            match self.matching_from_column_index_to_row_index.evaluate_function( entry_to_eliminate.key() ) {
                None => {
                    self.entries_to_eliminate.head = Some( entry_to_eliminate ); // put back the entry; we cannot clear it
                    None
                }, Some( row_index_of_eliminating_row ) => {
                    // TO ELIMINATE `entry_to_eliminate`, WE ADD A SCALAR MULTIPLE OF A VECTOR TAKEN FROM THE MATRIX; DENOTE THE UNSCLAED VECTOR BY v
                    let mut seed_of_eliminating_iterator    =   self.matrix.row( & row_index_of_eliminating_row ).into_iter();
                    
                    // POP THE LEADING ENTRY OUT OF THE UNSCALED VECTOR v; USE THIS ENTRY TO DETERMINE THE SCALE FACTOR BY WHICH WE'LL MULTIPLY THE REST OF v
                    let eliminating_entry           =   seed_of_eliminating_iterator.next().unwrap();

                    // unsafe {
                    //     println!("row_index_of_eliminating_row =  {:?}", D(&row_index_of_eliminating_row) );                
                    //     println!("eliminating entry =  {:?}", D(&eliminating_entry) );
                    //     let seed_of_eliminating_iterator_2    =   self.matrix.row( row_index_of_eliminating_row.clone() ).into_iter().collect_vec();
                    //     println!("eliminating iterator entries =  {:?}", D(&seed_of_eliminating_iterator_2) );
                    // }            

                    // THE SCALE FACTOR IS EQUAL TO -(entry_to_eliminate)/(eliminating_entry)
                    let scale_factor     =    self.ring_operator.negate(
                                                            self.ring_operator.divide( entry_to_eliminate.val(), eliminating_entry.val() )
                                                        );
                    
                    // SCALE v SO THAT ITS LEADING ENTRY WILL CANCEL `entry_to_eliminate`
                    // note that it won't *actually* cancel that entry in practice, because we already popped the leading entry off of v and off of entries_to_eliminate
                    let eliminating_iterator       
                            =   TwoTypeIterator::Version1( seed_of_eliminating_iterator ) 
                                    .scale_by( scale_factor.clone(), self.ring_operator.clone() );


                    // MERGE THE (BEHEADED) SCALAR MULTIPLE OF v INTO THE HEAP, NAMELY `entries_to_eliminate`
                    hit_bulk_insert( &mut self.entries_to_eliminate.tail.unsimplified, vec![eliminating_iterator] ); 

                    Some( (row_index_of_eliminating_row, scale_factor) ) // recall that we have re-purposed `entry_to_eliminate` to hold the output of this function
                }
            }

        } else {
            None
        }
    }
}  






impl    < 
            ProblemVector,
            MatchingFromColumnIndicesToRowIndices,
            Matrix,                   
            RingOperator,
            OrderOperator,
        > 

        Clone for    

        RowEchelonSolver<
                ProblemVector,
                MatchingFromColumnIndicesToRowIndices,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                              IntoIterator< 
                                                        Item        =   Matrix::RowEntry,    // the iterator runs over entries of the same type as the matrix        
                                                        IntoIter:       Clone,
                                                    >, 
        MatchingFromColumnIndicesToRowIndices:      Clone + EvaluateFunction< Matrix::ColumnIndex, Option< Matrix::RowIndex > >,        
        Matrix:                                     Clone + MatrixOracle,
        RingOperator:                               Clone + DivisionRingOperations< Element = Matrix::Coefficient >,
        OrderOperator:                              Clone + JudgePartialOrder <  Matrix::RowEntry >,               
        Matrix::RowEntry:                           KeyValSet,    
        Simplify<IteratorsMergedInSortedOrder<Scale<TwoTypeIterator<Matrix::Row, RequireStrictAscentWithPanic<ProblemVector::IntoIter, OrderOperator>>, RingOperator>, OrderOperator>, RingOperator>: Clone,     
{
    fn clone(&self) -> Self {
        RowEchelonSolver{ 
                ring_operator: self.ring_operator.clone(),
                matrix: self.matrix.clone(),
                entries_to_eliminate: self.entries_to_eliminate.clone(),
                matching_from_column_index_to_row_index: self.matching_from_column_index_to_row_index.clone(),
            }
    }
}










//  ======================================================================================
//  ======================================================================================
//  ======================================================================================
//  DESCENDING MINOR SOLVE
//  ======================================================================================
//  ======================================================================================
//  ======================================================================================







//  ======================================================================================
//  ECHELON MINOR DESCEND SOLVE WITH ROW INDEXS
//  ======================================================================================




// ECHELON SOLVE W/ ROW INDEXS (ITERATOR)
// ---------------------------------------------------------------------------


/// Iterates over the entries of the solution `x` to a matrix equation `Ax = b`, where
/// (i) each entry of `x`, which represents a pair of form `(column_index, coefficient)`, is replaced by an entry representing `(match(column_index), coefficient)`, where
/// `match(column_index)` is the matched row index, and
/// (ii) entries appear in descending order, according to row index
/// 
/// # Assumptions
/// 
/// - `b` is a sparse vector iterator indexed by the row indices of `A`, and the entries of `b` appear in descending order of index
/// - There is a partial bijection `match: row_index -> column_index` from the row indices of `A` to the column indices of `A`; concretely, this means a bijection from a subset of the row indices to a subset of the column indices
/// - Whenever `column_index = match( row_index )`, the leading entry of `A.column_reverse( column_index )` has index `row_index`
/// - A solution to `Ax = b` exists
/// 
/// # What happens when an assumption is violated
/// 
/// This iterator will generate entries in the manner described below (under heading "Implementation"), whether or not there exists a solution to `Ax = b`.
/// After the iterator is exhausted, the user can call `x.remainder()`; the resulting iterator, `y`, represents what is "left over" from the vector `b` after some of its entries are eliminated by adding
/// scalar multiples of rows of `A`.  If `y` is empty, then a solution exists; in particular, `x` is a feasible solution.  On the other hand, if `y` is
/// nonempty then no solution exists, and `x` can be regarded merely as a failed attempt.
/// 
/// # Implementation
/// 
/// The iterator proceeds by repeating the following steps: (i) if `b = 0` has no leading entry, return `None`; we have found a valid solution, `x`, (ii) otherwise `b` has a 
/// leading (nonzero) entry.  Let `column_index` be the index of this entry.  (iii) if `column_index` is unmatched, return `None`; there exists no solution.  (iv) otherwise, let `row_index` be the row index matched to `column_index`.  Add a scalar multiple of `A.row( &row_index )` to `b` to 
/// eliminate its leading entry.  Let `alpha` denote the scalar by which we multiply, and return `Some( (row_index, -alpha) )`.
/// 
/// Under the hood, `ColumnEchelonSolverReverseReindexed` is simply a wrapper for a struct `EchelonSolverMajorAscendWithMinorKeys` which
/// contains a reference to the antitranspose of the matrix; see source code for full details.
#[derive(Dissolve)]
pub struct ColumnEchelonSolverReverseReindexed<
                    ProblemVector,
                    MatchingFromRowIndicesToColumnIndices,
                    Matrix,                     
                    RingOperator,
                    OrderOperatorUnreversed,
                > 
    where 
        ProblemVector:                              IntoIterator< Item = Matrix::ColumnEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromRowIndicesToColumnIndices:      EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColumnIndex > >,
        Matrix:                                     MatrixOracle,
        Matrix::ColumnEntry:                        KeyValSet < Key=Matrix::RowIndex, Val=Matrix::Coefficient >,      
        RingOperator:                               Clone + DivisionRingOperations< Element = Matrix::Coefficient >,
        OrderOperatorUnreversed:                    Clone + JudgePartialOrder <  Matrix::ColumnEntry >,                   

{   
    antitransposed_solver:      RowEchelonSolverReindexed< // the column indices of the antitransposed matrix are the row indices of the un-antitransposed matrix
                                        ProblemVector,
                                        MatchingFromRowIndicesToColumnIndices,
                                        OrderAntiTranspose< Matrix >,
                                        RingOperator,
                                        ReverseOrder< OrderOperatorUnreversed >,
                                    >,
} 

impl    < 
            ProblemVector,
            MatchingFromRowIndicesToColumnIndices,
            Matrix,                   
            RingOperator,
            OrderOperatorUnreversed,
        >  

        ColumnEchelonSolverReverseReindexed<
                ProblemVector,
                MatchingFromRowIndicesToColumnIndices,
                Matrix,                   
                RingOperator,
                OrderOperatorUnreversed,
            > 

    where 
        ProblemVector:                              IntoIterator< Item = Matrix::ColumnEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromRowIndicesToColumnIndices:      EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColumnIndex > >,
        Matrix:                                     MatrixOracle,
        Matrix::ColumnEntry:                        KeyValSet < Key=Matrix::RowIndex, Val=Matrix::Coefficient >,      
        RingOperator:                               Clone + DivisionRingOperations< Element = Matrix::Coefficient >,
        OrderOperatorUnreversed:                    Clone + JudgePartialOrder <  Matrix::ColumnEntry >,     

{

    /// Attempts to solve `xA = b` for `x`; see [`EchelonSolverMajorAscendWithMinorKeys`].
    pub fn solve(
                b:                              ProblemVector,                
                a:                              Matrix,
                matching_from_row_index_to_column_index: MatchingFromRowIndicesToColumnIndices,
                ring_operator:                  RingOperator,
                order_operator_unreversed:    OrderOperatorUnreversed,
            )
            ->
            QuotientRemainderSolver<
                    ColumnEchelonSolverReverseReindexed< 
                            ProblemVector,
                            MatchingFromRowIndicesToColumnIndices,
                            Matrix,                     
                            RingOperator,
                            OrderOperatorUnreversed,
                        >             
                > 
    {
        QuotientRemainderSolver::new(
            ColumnEchelonSolverReverseReindexed{
                    // Note that we use EchelonSolverMajorAscendWithMinorKeys, since the column indices of the antitranspose are the row indices of the original matrix
                    antitransposed_solver:   RowEchelonSolverReindexed::solve(
                                                    b,
                                                    OrderAntiTranspose::new( a ),
                                                    matching_from_row_index_to_column_index,
                                                    ring_operator,
                                                    ReverseOrder::new( order_operator_unreversed ),
                                                ).quotient()
                }
        )
    }

}








impl    < 
            ProblemVector,
            MatchingFromRowIndicesToColumnIndices,
            Matrix,                   
            RingOperator,
            OrderOperatorUnreversed,
        >  

        IntoRemainder for

        ColumnEchelonSolverReverseReindexed<
                ProblemVector,
                MatchingFromRowIndicesToColumnIndices,
                Matrix,                   
                RingOperator,
                OrderOperatorUnreversed,
            > 

    where 
        ProblemVector:                              IntoIterator< Item = Matrix::ColumnEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromRowIndicesToColumnIndices:      EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColumnIndex > >,
        Matrix:                                     MatrixOracle,      
        Matrix::ColumnEntry:                        KeyValSet < Key=Matrix::RowIndex, Val=Matrix::Coefficient >,      
        RingOperator:                               Clone + DivisionRingOperations< Element = Matrix::Coefficient >,
        OrderOperatorUnreversed:                    Clone + JudgePartialOrder <  Matrix::ColumnEntry >,   

{

    type Remainder = 
        HeadTail<
                LinearCombinationSimplified<   
                        TwoTypeIterator< // this enum allows us to treat two different iterator types as a single type
                                < Matrix::ColumnReverse as IntoIterator >::IntoIter,  // the iterators returned by the matrix                                                                
                                RequireStrictAscentWithPanic< ProblemVector::IntoIter, ReverseOrder< OrderOperatorUnreversed > >, // the vector we have tried to eliminate
                            >,
                        RingOperator, ReverseOrder< OrderOperatorUnreversed >
                    >
            >;

    /// Returns an iterator that runs over all the entries that have not yet been eliminated, consuming the solver
    /// in the process.
    /// 
    /// Recall that the solver works by iteratively adding linear multiples rows or columns of matrix `A` to the
    /// problem vector `b`, in order to solve `xA=b`.  If `v_k` is the vector we add on iteration `k`, then the remainder after iteration `k` is
    /// `b - v_0 - v_1 ... - v_k`.  The leading entry of the remainder is the entry we must try to eliminate on 
    /// iteration `v_{k+1}`.  Note that the remainder can be nonzero even when `self.next() = None`.  This happens
    /// if and only if there exists no solution to `xA = b`.
    fn into_remainder( self ) -> Self::Remainder                 
    {
        self.antitransposed_solver.into_remainder()
    }    
}









impl    <
            ProblemVector,
            MatchingFromRowIndicesToColumnIndices,
            Matrix,                   
            RingOperator,
            OrderOperatorUnreversed,
        > 

        Iterator for    

        ColumnEchelonSolverReverseReindexed<
                ProblemVector,
                MatchingFromRowIndicesToColumnIndices,
                Matrix,                   
                RingOperator,
                OrderOperatorUnreversed,
            > 

    where 
        ProblemVector:                              IntoIterator< Item = Matrix::ColumnEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromRowIndicesToColumnIndices:      EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColumnIndex > >,
        Matrix:                                     MatrixOracle,      
        Matrix::ColumnEntry:                        KeyValSet < Key=Matrix::RowIndex, Val=Matrix::Coefficient >,      
        RingOperator:                               Clone + DivisionRingOperations< Element = Matrix::Coefficient >,
        OrderOperatorUnreversed:                              Clone + JudgePartialOrder <  Matrix::ColumnEntry >,    

{
    type Item = Matrix::ColumnEntry;

    fn next( &mut self ) -> Option<Self::Item> { 
        self.antitransposed_solver.next()
    }
}  






impl    < 
            ProblemVector,
            MatchingFromRowIndicesToColumnIndices,
            Matrix,                   
            RingOperator,
            OrderOperatorUnreversed,
        > 

        Clone for    

ColumnEchelonSolverReverseReindexed<
                    ProblemVector,
                    MatchingFromRowIndicesToColumnIndices,
                    Matrix,                     
                    RingOperator,
                    OrderOperatorUnreversed,
                > 
    where 
        ProblemVector:                              IntoIterator< Item = Matrix::ColumnEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromRowIndicesToColumnIndices:      EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColumnIndex > >,
        Matrix:                                     MatrixOracle,      
        Matrix::ColumnEntry:                        KeyValSet < Key=Matrix::RowIndex, Val=Matrix::Coefficient >,      
        RingOperator:                               Clone + DivisionRingOperations< Element = Matrix::Coefficient >,
        OrderOperatorUnreversed:                              Clone + JudgePartialOrder <  Matrix::ColumnEntry >,  
        
        // our struct contains a single field; the following is the type of that field; we simply require this type to implement clone!
        RowEchelonSolverReindexed< 
                                        ProblemVector,
                                        MatchingFromRowIndicesToColumnIndices,
                                        OrderAntiTranspose< Matrix >,
                                        RingOperator,
                                        ReverseOrder< OrderOperatorUnreversed >,
                                    >:              Clone,                     
{
    fn clone(&self) -> Self {
        ColumnEchelonSolverReverseReindexed{ 
                antitransposed_solver: self.antitransposed_solver.clone()
            }
    }
}
















//  ======================================================================================
//  ECHELON MINOR DESCEND SOLVE WITH MINOR KEYS
//  ======================================================================================




// ECHELON SOLVE W/ MINOR KEYS (ITERATOR)
// ---------------------------------------------------------------------------


/// Iterates over the entries of the solution `x` to a matrix equation `Ax = b`, where entries appear in descending order, according to  index
/// 
/// # Assumptions
/// 
/// - `b` is a sparse [vector](oat_rust::algebra::vectors) indexed by the row indices of `A`, and the entries of `b` appear in descending order of row index
/// - There is a partial bijection `match: row_index -> column_index` from the row indices of `A` to the column indices of `A`; concretely, this means a bijection from a subset of the row indices to a subset of the column indices
/// - Whenever `column_index = match( row_index )`, the leading entry of `A.column_reverse( column_index )` has index `row_index`
/// - A solution to `Ax = b` exists
/// 
/// # What happens when an assumption is violated
/// 
/// This iterator will generate entries in the manner described below (under heading "Implementation"), whether or not there exists a solution to `Ax = b`.
/// After the iterator is exhausted, the user can call `x.remainder()`; the resulting iterator, `y`, represents what is "left over" from the vector `b` after some of its entries are eliminated by adding
/// scalar multiples of rows of `A`.  If `y` is empty, then a solution exists; in particular, `x` is a feasible solution.  On the other hand, if `y` is
/// nonempty then no solution exists, and `x` can be regarded merely as a failed attempt.
/// 
/// # Implementation
/// 
/// The iterator proceeds by repeating the following steps: (i) If `b` has no leading entry then `b=0`. In this case return `None`, because we have found a valid solution, `x`. (ii) Otherwise `b` has a 
/// leading (nonzero) entry.  Let `column_index` be the index of this entry.  (iii) If `column_index` is unmatched, return `None`; there exists no solution.  (iv) Otherwise, let `row_index` be the row index matched to `column_index`.  Add a scalar multiple of `A.row( &row_index )` to `b` to 
/// eliminate its leading entry.  Let `alpha` denote the scalar by which we multiply, and return `Some( (row_index, -alpha) )`.
pub struct ColumnEchelonSolverReverse<
                    ProblemVector,
                    MatchingFromRowIndicesToColumnIndices,
                    Matrix,                     
                    RingOperator,
                    OrderOperator,
                > 
    where 
        ProblemVector:                              IntoIterator< Item = Matrix::ColumnEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromRowIndicesToColumnIndices:      EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColumnIndex > >,
        Matrix:                                     MatrixOracle,  
        Matrix::ColumnEntry:                        KeyValSet < Key=Matrix::RowIndex, Val=Matrix::Coefficient >,      
        RingOperator:                               Clone + DivisionRingOperations< Element = Matrix::Coefficient >,
        OrderOperator:                              Clone + JudgePartialOrder <  Matrix::ColumnEntry >,                   

{   
    antitransposed_solver:      RowEchelonSolver< // the row indices of the antitransposed matrix are the column indices of the un-antitransposed matrix
                                        ProblemVector,
                                        MatchingFromRowIndicesToColumnIndices,
                                        OrderAntiTranspose< Matrix >,                 
                                        RingOperator,
                                        ReverseOrder< OrderOperator >,
                                    >,
} 

impl    < 
            ProblemVector,
            MatchingFromRowIndicesToColumnIndices,
            Matrix,                   
            RingOperator,
            OrderOperator,
        >  

        ColumnEchelonSolverReverse<
                ProblemVector,
                MatchingFromRowIndicesToColumnIndices,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                              IntoIterator< Item = Matrix::ColumnEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromRowIndicesToColumnIndices:      EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColumnIndex > >,
        Matrix:                                     MatrixOracle,  
        Matrix::ColumnEntry:                        KeyValSet < Key=Matrix::RowIndex, Val=Matrix::Coefficient >,      
        RingOperator:                               Clone + DivisionRingOperations< Element = Matrix::Coefficient >,
        OrderOperator:                              Clone + JudgePartialOrder <  Matrix::ColumnEntry >,    

{

    /// Attempts to solve `xA = b` for `x`; see [`EchelonSolverMajorAscendWithMinorKeys`].
    pub fn solve(
                b:                              ProblemVector,                
                a:                              Matrix,
                matching_from_row_index_to_column_index: MatchingFromRowIndicesToColumnIndices,
                ring_operator:                  RingOperator,
                order_operator:               OrderOperator,
            )
            ->
            QuotientRemainderSolver<
                    ColumnEchelonSolverReverse< 
                            ProblemVector,
                            MatchingFromRowIndicesToColumnIndices,
                            Matrix,                     
                            RingOperator,
                            OrderOperator,
                        >             
                > 
    {
        QuotientRemainderSolver::new(
            ColumnEchelonSolverReverse{
                // Note that we use RowEchelonSolver, since the row indices of the antitranspose are the column indices of the original matrix
                antitransposed_solver:   RowEchelonSolver::solve(
                                                b,
                                                OrderAntiTranspose::new( a ),
                                                matching_from_row_index_to_column_index,
                                                ring_operator,
                                                ReverseOrder::new( order_operator ),
                                            ).quotient()
            }
        )
    }
}





impl    < 
            ProblemVector,
            MatchingFromRowIndicesToColumnIndices,
            Matrix,                   
            RingOperator,
            OrderOperator,
        >  

        IntoRemainder for

        ColumnEchelonSolverReverse<
                ProblemVector,
                MatchingFromRowIndicesToColumnIndices,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                              IntoIterator< Item = Matrix::ColumnEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromRowIndicesToColumnIndices:      EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColumnIndex > >,
        Matrix:                                     MatrixOracle,    
        Matrix::ColumnEntry:                        KeyValSet < Key=Matrix::RowIndex, Val=Matrix::Coefficient >,      
        RingOperator:                               Clone + DivisionRingOperations< Element = Matrix::Coefficient >,
        OrderOperator:                              Clone + JudgePartialOrder <  Matrix::ColumnEntry >,   

{

    type Remainder = 
            HeadTail<
                    LinearCombinationSimplified<   
                            TwoTypeIterator< // this enum allows us to treat two different iterator types as a single type
                                    < Matrix::ColumnReverse as IntoIterator >::IntoIter,  // the iterators returned by the matrix                                                                
                                    RequireStrictAscentWithPanic< ProblemVector::IntoIter, ReverseOrder< OrderOperator > >, // the vector we have tried to eliminate
                                >,
                            RingOperator, ReverseOrder< OrderOperator >
                        >
                >;

    /// Returns an iterator that runs over all the entries that have not yet been eliminated, consuming the solver
    /// in the process.
    /// 
    /// Recall that the solver works by iteratively adding linear multiples rows or columns of matrix `A` to the
    /// problem vector `b`, in order to solve `xA=b`.  If `v_k` is the vector we add on iteration `k`, then the remainder after iteration `k` is
    /// `b - v_0 - v_1 ... - v_k`.  The leading entry of the remainder is the entry we must try to eliminate on 
    /// iteration `v_{k+1}`.  Note that the remainder can be nonzero even when `self.next() = None`.  This happens
    /// if and only if there exists no solution to `xA = b`.
    fn into_remainder( self ) -> Self::Remainder     
    {
        self.antitransposed_solver.into_remainder()
    }    
}





impl    <
            ProblemVector,
            MatchingFromRowIndicesToColumnIndices,
            Matrix,                   
            RingOperator,
            OrderOperator,
        > 

        Iterator for    

        ColumnEchelonSolverReverse<
                ProblemVector,
                MatchingFromRowIndicesToColumnIndices,
                Matrix,                   
                RingOperator,
                OrderOperator,
            > 

    where 
        ProblemVector:                              IntoIterator< Item = Matrix::ColumnEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromRowIndicesToColumnIndices:      EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColumnIndex > >,
        Matrix:                                     MatrixOracle,    
        Matrix::ColumnEntry:                        KeyValSet < Key=Matrix::RowIndex, Val=Matrix::Coefficient >,      
        RingOperator:                               Clone + DivisionRingOperations< Element = Matrix::Coefficient >,
        OrderOperator:                              Clone + JudgePartialOrder <  Matrix::ColumnEntry >,   

{
    type Item = (Matrix::ColumnIndex, Matrix::Coefficient);

    fn next( &mut self ) -> Option<Self::Item> { 
        self.antitransposed_solver.next()
    }
}  





impl    < 
            ProblemVector,
            MatchingFromRowIndicesToColumnIndices,
            Matrix,                   
            RingOperator,
            OrderOperator,
        > 

        Clone for    

        ColumnEchelonSolverReverse<
                ProblemVector,
                MatchingFromRowIndicesToColumnIndices,
                Matrix,                     
                RingOperator,
                OrderOperator,
            > 
    where 
        ProblemVector:                              IntoIterator< Item = Matrix::ColumnEntry >, // the iterator runs over entries of the same type as the matrix
        MatchingFromRowIndicesToColumnIndices:      EvaluateFunction< Matrix::RowIndex, Option< Matrix::ColumnIndex > >,
        Matrix:                                     MatrixOracle,    
        Matrix::ColumnEntry:                        KeyValSet < Key=Matrix::RowIndex, Val=Matrix::Coefficient >,      
        RingOperator:                               Clone + DivisionRingOperations< Element = Matrix::Coefficient >,
        OrderOperator:                              Clone + JudgePartialOrder <  Matrix::ColumnEntry >,   
        
        // The struct ColumnEchelonSolverReverse contains only one field; below is the name of that field, and the requirement that it implement Clone
        RowEchelonSolver< // the row indices of the antitransposed matrix are the column indices of the un-antitransposed matrix
            ProblemVector,
            MatchingFromRowIndicesToColumnIndices,
            OrderAntiTranspose< Matrix >,                 
            RingOperator,
            ReverseOrder< OrderOperator >,
        >:                                          Clone,                       
{
    fn clone(&self) -> Self {
        ColumnEchelonSolverReverse{ 
                antitransposed_solver: self.antitransposed_solver.clone()
            }
    }
}












//  =========================================================================
//  =========================================================================
//  =========================================================================
//  IMPLEMENTATIONS OF CLONE
//  =========================================================================
//  =========================================================================
//  =========================================================================












































//  TESTS
//  ===========================================================================


#[cfg(test)]
mod doctring_tests {



    #[test]
    fn doc_test_module() {
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::operations::solve::echelon::{RowEchelonSolver, ColumnEchelonSolverReverse};
                
        use crate::algebra::rings::types::field_prime_order::BooleanField;
        use crate::algebra::vectors::operations::VectorOperations;
        use crate::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;
        use crate::utilities::order::{OrderOperatorByKey};
        use itertools::Itertools;
        use assert_panic::assert_panic;

        // define the ring operator
        let ring_operator = BooleanField::new();

        // build a matrix A with an invertible upper triangular submatrix
        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![                                     ], 
                                    vec![                                     ],                                     
                                    vec![ (0,true),  (1,true), (2,true)       ], 
                                    vec![            (1,true), (2,true)       ],
                                    vec![                      (2,true)       ],                                       
                                ],
                            ).ok().unwrap();    

        // ======================================================================================                        
        // SOLVE xA = b FOR x
        // ======================================================================================
                        
        let b = vec![ (0,true), (1,true) ];   // define a sparse vector b           

        // define the key-matrix_to_factor with a wrapper around a closure operator (|x| -> expression(x))      
        // --------------------------------------------------------------------------------------
                        
        let partial_bijection = |x: usize| { if x < 3 { Some(x+2) } else { None } };
        let division =  RowEchelonSolver::solve(
                                                                        b.clone(), // entries must appear in strictly ASCENDING order
                                                                        & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
                                                                        EvaluateFunctionFnMutWrapper::new( partial_bijection ), // maps *matched* column indices to matched row indices
                                                                        BooleanField::new(), // defines the ring operations
                                                                        OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
                                                                    );
                                                                
        let product = division
                            .quotient()
                            .multiply_self_as_a_row_vector_with_matrix_custom( 
                                &matrix, 
                                ring_operator, 
                                OrderOperatorByKey::new() 
                            );
        assert_eq!( product.collect_vec(), b );   // check the solution, i.e. check that xA = b  


        // define the key-matrix_to_factor with vector
        // --------------------------------------------------------------------------------------   
                        
        // create a solver to solve xA = b for x
        let division =  RowEchelonSolver::solve(
                                                                        b.clone(), // entries must appear in strictly ASCENDING order
                                                                        & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
                                                                        vec![2,3,4], // maps *matched* column indices to matched row indices
                                                                        BooleanField::new(), // defines the ring operations
                                                                        OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
                                                                    );
                                                                
        // check the solution, i.e. check that xA = b
        let product = division
                            .clone()
                            .quotient()
                            .multiply_self_as_a_row_vector_with_matrix_custom( 
                                &matrix, 
                                ring_operator, 
                                OrderOperatorByKey::new() 
                            );
        assert_eq!( product.collect_vec(), b );  

        // different look-up methods
        // --------------------------------------------------------------------------------------             

        assert_eq!( division.clone().solution(), Some( division.clone().quotient().collect_vec() ) );        

        // ======================================================================================
        // SOLVE Ay = c FOR y
        // ======================================================================================

        let c = vec![ (3,true), (2,true) ];   // define a sparse vector c  

        // define the key-matrix_to_factor with a wrapper around a closure operator (|x| -> expression(x))      
        // --------------------------------------------------------------------------------------        
                                                                                                                                
        let partial_bijection = |x: usize| { if (2..=4).contains(&x) { Some(x-2) } else { None } };
        let division =  ColumnEchelonSolverReverse::solve(
                                                                        c.clone(), // entries must appear in strictly DESCENDING order
                                                                        & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
                                                                        EvaluateFunctionFnMutWrapper::new( partial_bijection ), // maps *matched* row indices to matched column indices
                                                                        BooleanField::new(), // defines the ring operations
                                                                        OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever j ≤ i
                                                                    );
                                                                
        let product = division
                            .quotient()
                            .multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order_custom( 
                                &matrix, 
                                ring_operator, 
                                OrderOperatorByKey::new()
                            );
        assert_eq!( product.collect_vec(), c );   // check the solution, i.e. check that Ay = c      

        // ======================================================================================                        
        // GET A REMAINDER
        // ======================================================================================


        // quotienting out columns
        // --------------------------------------------------------------------------------------   

        let d = vec![ (3,true), (2,true), (1,true), (0,true)  ];

        let partial_bijection = |x: usize| { if (2..=4).contains(&x) { Some(x-2) } else { None } };
        let division =  ColumnEchelonSolverReverse::solve(
            d, // entries must appear in strictly DESCENDING order
            & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
            EvaluateFunctionFnMutWrapper::new( partial_bijection ), // maps *matched* row indices to matched column indices
            BooleanField::new(), // defines the ring operations
            OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever j ≤ i
        );    
        let remainder = division.remainder().collect_vec();                                                        
        assert!( remainder.eq( &vec![ (1,true), (0,true) ] ) );   // the elimination procedure clears all but the first two entries

        // quotienting out rows
        // --------------------------------------------------------------------------------------   

        let d = vec![ (0,true), (1,true), (2,true), (3,true), ];
        let division =  RowEchelonSolver::solve(
                                                                        d, // entries must appear in strictly ASCENDING order
                                                                        & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
                                                                        vec![2,3,4], // maps *matched* column indices to matched row indices
                                                                        BooleanField::new(), // defines the ring operations
                                                                        OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
                                                                    );                                                                
        assert!( division.clone().remainder().eq( vec![ (3,true) ] ) );   // the elimination procedure clears the first three columns, but cannot clearn the fourth

        // alternate look-up method
        // --------------------------------------------------------------------------------------   
                                                          
        assert!( division.clone().quotient_remainder().remainder.eq( vec![ (3,true) ] ) );   // the elimination procedure clears all but the first two entries        
         

        // ======================================================================================                        
        // PANIC FOR UNSORTED INPUTS
        // ======================================================================================

        // panics if we pass an unsorted input, and pull all entries from the quotient and remainder
        // --------------------------------------------------------------------------------------   

        let d = vec![ (1,true), (0,true), (3,true), (2,true), ];
        let division =  RowEchelonSolver::solve(
                                                                        d, // entries must appear in strictly ASCENDING order
                                                                        & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
                                                                        vec![2,3,4], // maps *matched* column indices to matched row indices
                                                                        BooleanField::new(), // defines the ring operations
                                                                        OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
                                                                    );                                                                
        assert_panic!( for _ in division.solver() {} );   // the elimination procedure clears the first three columns, but cannot clearn the fourth


        // the solver *may possibly* compute a quotient without throwing an error, if it finds an entry it can't eliminate before hitting the first consecutive out-of-order pair
        // --------------------------------------------------------------------------------------   

        let d = vec![ (0,true), (1,true), (2,true), (3,true), (4,true), (2,true), ];
        let division =  RowEchelonSolver::solve(
                                                                        d, // entries must appear in strictly ASCENDING order
                                                                        & matrix, // the matrix oracle traits are only implemented for *references* to VecOfVec, not on VecOfVec itself
                                                                        vec![2,3,4], // maps *matched* column indices to matched row indices
                                                                        BooleanField::new(), // defines the ring operations
                                                                        OrderOperatorByKey::new(), // an object that asserts (i, s) ≤ (j, t) whenever i ≤ j
                                                                    );
        
        let mut qr = division.quotient_remainder();
        assert!( qr.quotient.clone().eq( &vec![ (2,true) ] ) ); // no panic when computing the quotient
        assert_eq!( qr.remainder.next(), Some( (3,true) ) ); // no panic when computing the first element of the remainder 
        assert_panic!( for _ in qr.remainder {} );   // panics once we exhaust the problem vector

    }

}










//  TESTS
//  ===========================================================================


#[cfg(test)]
mod tests {

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    
    use itertools::Itertools;
    use crate::{algebra::matrices::types::vec_of_vec::sorted::VecOfVec, algebra::rings::types::field_prime_order::PrimeOrderField};


    // Invert upper triangular matrices
    // =====================================================================================================================

    /// [this is a subroutine / a helper function] 
    /// Verifies that the output of a solver actually solves the problem Ax = b
    /// **IT IS IMPORTANT THAT THE MATRIX BE IN PARTIAL ECHECLON FORM**
    #[cfg(test)] // because this is a helper function, we mark it with the macro shown left to prevent "warning: unused function" from showing up
    fn test_echelon_solve_on_invertible_uppertriangular_matrix< 'a, MatchingFromColumnIndicesToRowIndices, MatchingFromRowIndicesToColumnIndices, Coefficient, RingOperator >( 
                matrix:                             & VecOfVec< usize, Coefficient >, 
                matching_from_column_index_to_row_index:     MatchingFromColumnIndicesToRowIndices,
                matching_from_row_index_to_column_index:     MatchingFromRowIndicesToColumnIndices,
                ring_operator:                      RingOperator, 
                matrix_size:                        usize 
            ) 
        where   Coefficient:                                Clone + PartialEq + Ord + std::fmt::Debug,
                RingOperator:                               DivisionRingOperations< Element = Coefficient > + Clone,
                MatchingFromColumnIndicesToRowIndices:      Clone + EvaluateFunction< usize, Option< usize > >,
                MatchingFromRowIndicesToColumnIndices:      Clone + EvaluateFunction< usize, Option< usize > >,                
    {
        use crate::{algebra::matrices::operations::multiply::{multiply_row_vector_with_matrix, multiply_column_vector_with_matrix_and_return_reversed}, utilities::order::OrderOperatorAuto};

        


        // generate random vectors of 0's and 1's; try to get one vector each with 0, 1, 2, 3, and `matrix_size` nonzero entries in total
        for num_nonzeros in [0, 1, 2, 3, matrix_size].iter().map(|x| std::cmp::min(x, &matrix_size) ) {
            for vector_support in (0 .. matrix_size).combinations( *num_nonzeros ) {
                let mut vec = vector_support
                                .iter()
                                .cloned()
                                .map(|x| ( x, RingOperator::one() ))
                                .collect_vec();
                vec.sort();


                //  row vector solutions (corresponding to xA = b)
                //  --------------------------------------------------------------------                

                // compute a row vector solution with column indices
                let solution_major_with_minor_keys = RowEchelonSolverReindexed::solve(
                                        vec.iter().cloned(),                    
                                        matrix,
                                        matching_from_column_index_to_row_index.clone(),
                                        ring_operator.clone(),
                                        OrderOperatorAuto,                    
                                    )
                                    .quotient()
                                    .peekable()
                                    .simplify( ring_operator.clone() );
                // verify the solution with column indices
                itertools::assert_equal(
                    vec.iter().cloned().peekable().simplify( ring_operator.clone() ),
                    multiply_row_vector_with_matrix( 
                            solution_major_with_minor_keys, 
                            matrix, 
                            ring_operator.clone(), 
                            OrderOperatorAuto
                        )
                );

                // compute a row vector solution with row indices
                let solution_major_with_major_keys = RowEchelonSolver::solve(
                                        vec.iter().cloned(),                    
                                        matrix,
                                        matching_from_column_index_to_row_index.clone(),
                                        ring_operator.clone(),
                                        OrderOperatorAuto,                    
                                    )
                                    .quotient()
                                    .peekable()
                                    .simplify( ring_operator.clone() );
                // verify the solution with row indices
                itertools::assert_equal(
                    vec.iter().cloned().peekable().simplify( ring_operator.clone() ),
                    multiply_row_vector_with_matrix( 
                            solution_major_with_major_keys, 
                            matrix, 
                            ring_operator.clone(), 
                            OrderOperatorAuto
                        )
                );    


                //  column vector solutions (corresponding to Ax = b)
                //  --------------------------------------------------------------------

                // REVERSE THE ORDER OF ENTRIES, SINCE WE ARE NOW DEALING WITH DESCENDING VIEWS
                vec.reverse(); 

                // compute a column vector solution with column indices
                let solution_minor_with_minor_keys = ColumnEchelonSolverReverse::solve(
                                        vec.iter().cloned(),                    
                                        matrix,
                                        matching_from_row_index_to_column_index.clone(),
                                        ring_operator.clone(),
                                        OrderOperatorAuto,                    
                                    )
                                    .quotient()
                                    .peekable()
                                    .simplify( ring_operator.clone() );
                // verify the solution with column indices
                itertools::assert_equal(
                    vec.iter().cloned().peekable().simplify( ring_operator.clone() ),
                    multiply_column_vector_with_matrix_and_return_reversed( 
                            solution_minor_with_minor_keys, 
                            matrix, 
                            ring_operator.clone(), 
                            OrderOperatorAuto
                        )
                );

                // compute a column vector solution with row indices
                let solution_major_with_major_keys = ColumnEchelonSolverReverseReindexed::solve(
                                        vec.iter().cloned(),                    
                                        matrix,
                                        matching_from_row_index_to_column_index.clone(),
                                        ring_operator.clone(),
                                        OrderOperatorAuto,                    
                                    )
                                    .quotient()
                                    .peekable()
                                    .simplify( ring_operator.clone() );
                // verify the solution with row indices
                itertools::assert_equal(
                    vec.iter().cloned().peekable().simplify( ring_operator.clone() ),
                    multiply_column_vector_with_matrix_and_return_reversed( 
                            solution_major_with_major_keys, 
                            matrix, 
                            ring_operator.clone(), 
                            OrderOperatorAuto
                        )
                );                                
                

            }
        }                                                               
    }


    
    #[test]
    fn test_echelon_solve_on_specific_matrices() {
        use num::rational::Ratio;        
        // use crate::algebra::matrices::query::MajorDimension;
        use crate::algebra::rings::types::native::RingOperatorForNativeRustNumberType;
        use crate::utilities::functions::evaluate::EvaluateFunctionFnMutWrapper;

        // Define the integer type for the rational numbers we'll use
        type IntegerType = isize;
        let modulus = 1049;

        // Define shorthand functions for creating new rational numbers
        let r = Ratio::new;    
        let q = Ratio::from_integer;   

        // Define matching function from Matrix::ColumnIndex (=usize) to Matrix::RowIndex (=usize)
        let matching_column_index_to_row_index_closure = | x: usize | -> Option< usize > { Some( x ) };
        let matching_column_index_to_row_index = EvaluateFunctionFnMutWrapper::new( matching_column_index_to_row_index_closure );
        
        // Define the ring operators
        let ring_operator_q  =   RingOperatorForNativeRustNumberType::< Ratio<IntegerType> >::new();  
        let ring_operator_p  =   PrimeOrderField::new(modulus);                  

        // Test individual matrices        
        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![ (0,q(1)) ]
                                ],
                        ).ok().unwrap();
        test_echelon_solve_on_invertible_uppertriangular_matrix( &matrix, matching_column_index_to_row_index.clone(), matching_column_index_to_row_index.clone(), ring_operator_q, 1 );       

        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![ (0,q(1)),  (1,q(1)), ], 
                                    vec![            (1,q(1)), ],
                                ],
                        ).ok().unwrap();
        test_echelon_solve_on_invertible_uppertriangular_matrix( &matrix, matching_column_index_to_row_index.clone(), matching_column_index_to_row_index.clone(), ring_operator_q, 2 );       

        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![ (0,q(1)),  (1,q(-1)), (2,q(3))       ], 
                                    vec![            (1,q(-1)), (2,q(4))       ],
                                    vec![                       (2,q(5))       ],                                    
                                ],
                        ).ok().unwrap(); 
        test_echelon_solve_on_invertible_uppertriangular_matrix( &matrix, matching_column_index_to_row_index.clone(), matching_column_index_to_row_index.clone(), ring_operator_q, 3 );        
        
        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![ (0,r(4,1)),  (1,r(-6,1))                 ], 
                                    vec![              (1,r(4,1)),   (2,r(-2,1))   ],
                                    vec![                            (2,r(7,1))    ],                                    
                                ],
                        ).ok().unwrap();      
        test_echelon_solve_on_invertible_uppertriangular_matrix( &matrix, matching_column_index_to_row_index.clone(), matching_column_index_to_row_index.clone(), ring_operator_q, 3 );           

        let matrix  =   VecOfVec::new(
                                vec![ 
                                    vec![   (0,r(2,1)),     (1,r(6,1)),                     (3,r(8,1)),  ], 
                                    vec![                   (1,r(2,1)),     (2,r(9,1)),                  ],
                                    vec![                                   (2,r(2,1)),     (3,r(-4,1)), ],  
                                    vec![                                                   (3,r(4,1)),  ],                                  
                                ],
                        ).ok().unwrap();
        test_echelon_solve_on_invertible_uppertriangular_matrix( &matrix, matching_column_index_to_row_index.clone(), matching_column_index_to_row_index.clone(), ring_operator_q, 4 );         
        
        // MUCH LARGER randomized matrix to invert, with nonzero diagonal entries
        use rand::Rng;        // we use this module to generate random elements
        let matrix_size =   20;

        let mut rng = rand::thread_rng(); // this generates random integers
        let mut vec_of_vec = vec![];
        for row_index in 0 .. matrix_size {
            let coefficient_leading         =   rng.gen_range( 1 .. modulus );
            let mut new_vec     =   vec![ (row_index, coefficient_leading) ]; // start with a nonzero diagonal element            
            for q in row_index+1 .. matrix_size { // fill out the rest of this row of the matrix
                let coefficient   = rng.gen_range( 0 .. modulus );
                let flag = rng.gen_range(0usize .. 3usize);
                if      flag == 0 { new_vec.push( ( q, 0 )           ) }
                else if flag == 1 { new_vec.push( ( q, coefficient ) ) }
                else              { continue }
            }
            vec_of_vec.push( new_vec );  // the row is now filled out; push into the matrix
        }

        let matrix  =   VecOfVec::new(vec_of_vec).ok().unwrap(); // formally wrap the matrix in a VecOfVec struct
        test_echelon_solve_on_invertible_uppertriangular_matrix( &matrix, matching_column_index_to_row_index.clone(), matching_column_index_to_row_index, ring_operator_p, matrix_size );                 

    }

}
