//! Matrix-matrix and matrix-vector multiplication.
//! 
//! 
//! To multiply a matrix with a matrix, use one of the following
//! - [MatrixOracleOperations](crate::algebra::matrices::operations::MatrixOracleOperations), or
//! - [ProductMatrix::new]
//! 
//! To multiply a matrix with a vector, use one of the following
//! - [MatrixOracleOperations](crate::algebra::matrices::operations::MatrixOracleOperations)
//! - [VectorOperations](crate::algebra::vectors::operations::VectorOperations)
//! - the functions listed at the bottom of this page
//! 
//! **[MatrixOracleOperations](crate::algebra::matrices::operations::MatrixOracleOperations) and [VectorOperations] are only available
//! for matrix structs that implement `Sized`. 
//! If your matrices do not implement `Sized`, then you can use the functions listed at the bottom of this page
//! (for matrix-vector multiplication) or [ProductMatrix::new] (for matrix-matrix multiplication).
//! **
//! 
//! # See also
//! 
//! The [product matrix module](crate::algebra::matrices::types::Product), which defines the [ProductMatrix] type.


use crate::algebra::matrices::query::{ MatrixAlgebra, MatrixOracle };
use crate::algebra::matrices::operations::MatrixOracleOperations;
use crate::algebra::matrices::types::product::ProductMatrix;

use crate::algebra::rings::traits::SemiringOperations;

use crate::algebra::vectors::operations::{LinearCombinationSimplified, LinearCombinationUnsimplified, VectorOperations};
use crate::algebra::vectors::entries::{KeyValGet, KeyValSet};

use crate::utilities::iterators::merge::hit::hit_merge_by_predicate;
use crate::utilities::order::{JudgePartialOrder, ReverseOrder};



//  MATRIX - VECTOR MULTIPLICATION
//  ===========================================================================




//  ROW
//  ---------------------------------------------------------------------------

/// Calculates `v * M` (unsimplified)
/// 
/// Calculates `v * M`, where `v` is a sparse vector and `M` is a matrix.
/// 
/// The result is an iterator that runs over nonzero entries in ascending order of index.
/// 
/// **Note** the resulting iterator is not simplified -- it may return several entries with the *same index.*
/// 
/// # Examples
/// 
/// If `A` is a **row**-major representation of a matrix, and `v` is a sparse vector, then this function returns `vA`, a linear combination of **rows** of `A`.
/// If instead the representation is column-major, then this function returns `Av`.
/// 
/// ```
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;
/// use oat_rust::utilities::order::OrderOperatorByKey;
/// use oat_rust::algebra::matrices::operations::multiply::multiply_row_vector_with_matrix_unsimplified;
/// 
/// // a row-major matrix representation with usize row indices and isize column indices
/// let matrix      =   VecOfVec::new(   
///                         vec![   vec![ (0isize, 1), (1isize, 6) ], 
///                                 vec![ (0isize, 1), (1isize, 1) ],     ]     
///                     ).ok().unwrap();
/// // the vector [1, 1]
/// let vector      =   vec![ (0usize, 1), (1usize, 1) ];
/// 
/// // the product
/// let product     =   multiply_row_vector_with_matrix_unsimplified(
///                         vector,
///                         & matrix,
///                         PrimeOrderField::new(7), // the finite field of order 7
///                         OrderOperatorByKey::new(), // this compares tuples according to their lexicographic order
///                     );
/// // the product should return the following sequence of entries: (0,1), (0,1), (1,1), (1,6)   
/// itertools::assert_equal( product, vec![(0isize,1usize), (0,1), (1,1), (1,6)] )
/// ```
pub fn multiply_row_vector_with_matrix_unsimplified 
            < Matrix, RingOperator, SparseVecIter, OrderOperator> 
            ( 
                sparse_vec:         SparseVecIter,
                matrix:             Matrix, 
                ring_operator:      RingOperator,
                order_operator:     OrderOperator,
            ) 
            ->
            LinearCombinationUnsimplified
                < Matrix::Row, RingOperator, OrderOperator >            
    where 
        Matrix:                         MatrixOracle,
        Matrix::RowEntry:               KeyValSet,
        RingOperator:                   Clone + SemiringOperations< Element = Matrix::Coefficient >, // ring operators must typically implement Clone for operations such as simplification and dropping zeros
        OrderOperator:                  Clone + JudgePartialOrder<  Matrix::RowEntry >, // order comparators must often implement Clone when one wishes to construct new, ordered objects
        SparseVecIter:                  IntoIterator,
        SparseVecIter::Item:            KeyValGet < Key=Matrix::RowIndex, Val=Matrix::Coefficient >,
{
    // LinearCombinationUnsimplified{
    //     linear_combination_unsimplified:
            hit_merge_by_predicate(
                sparse_vec
                        .into_iter()
                        .map(   |x| 
                                matrix
                                    .row( & x.key() )
                                    .scale_by( x.val(), ring_operator.clone() )
                        ),
                    order_operator
            )        
    // }
}


/// Calculates `v * M` (unsimplified)
/// 
/// Calculates `v * M`, where `v` is a sparse vector and `M` is a matrix.
/// 
/// The result is an iterator of type [LinearCombinationSimplified], which represents a linear combination of the rows of `M`
/// - only entries with nonzero coefficients are returned
/// - entries are sorted in ascending order, according to index
/// - every entry has a distinct index (no repeat indices)
/// 
/// 
/// # Examples
/// 
/// 
/// ```
/// use oat_rust::algebra::matrices::operations::multiply::multiply_row_vector_with_matrix;
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;
/// use oat_rust::utilities::order::OrderOperatorByKey;
/// 
/// // a row-major matrix representation with usize row indices and isize column indices
/// let matrix      =   VecOfVec::new(   
///                         vec![   vec![ (0isize, 1), (1isize, 6) ], 
///                                 vec![ (0isize, 1), (1isize, 1) ],     ]     
///                     ).ok().unwrap();
/// // the vector [1, 1]
/// let vector      =   vec![ (0usize, 1), (1usize, 1) ];
/// 
/// // the product
/// let product     =   multiply_row_vector_with_matrix(
///                         vector,
///                         & matrix,
///                         PrimeOrderField::new(7), // the finite field of order 7
///                         OrderOperatorByKey::new(), // this compares tuples according to their lexicographic order
///                     );
/// // the product should equal vector [2 0]        
/// itertools::assert_equal( product, vec![(0isize, 2)] )
/// ```
pub fn multiply_row_vector_with_matrix 
            < Matrix, RingOperator, SparseVecIter, OrderOperator> 
            ( 
                sparse_vec:         SparseVecIter,
                matrix:             Matrix, 
                ring_operator:      RingOperator,
                order_operator:     OrderOperator,
            ) 
            ->
            LinearCombinationSimplified
                < Matrix::Row, RingOperator, OrderOperator >            
    where 
        Matrix:                         MatrixOracle,    
        Matrix::RowEntry:               KeyValSet,
        RingOperator:                   Clone + SemiringOperations< Element = Matrix::Coefficient >, // ring operators must typically implement Clone for operations such as simplification and dropping zeros        
        OrderOperator:                  Clone + JudgePartialOrder<  Matrix::RowEntry >, // order comparators must often implement Clone when one wishes to construct new, ordered objects
        SparseVecIter:                  IntoIterator,
        SparseVecIter::Item:            KeyValGet < Key=Matrix::RowIndex, Val=Matrix::Coefficient >,
{
    multiply_row_vector_with_matrix_unsimplified( 
            sparse_vec,
            matrix,
            ring_operator.clone(),
            order_operator,
        )
        // .linear_combination_unsimplified
        .simplify(ring_operator)
}        



//  COLUMN
//  ---------------------------------------------------------------------------


/// Returns the product of a matrix with a vector; the resulting iterator is an unsimplified linear combination of the matrix's columns.
/// 
/// The resulting iterator returns entries in descending (but **non-strictly** descending) order of index, provided that
/// [matrix.column_reverse(&column_index)](crate::algebra::matrices::query::MatrixOracle) returns entries in 
/// descending order according to index, consistent with the user-provided [order operator](crate::utilities::order).
/// 
/// **Note** the resulting iterator is not simplified -- it may return multiple entries with the same index.
/// 
/// ```
/// use oat_rust::algebra::matrices::operations::multiply::multiply_column_vector_with_matrix_and_return_unsimplified_and_reversed;
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;
/// use oat_rust::utilities::order::OrderOperatorAuto;
/// 
/// // a row-major matrix representation with usize row indices and isize column indices
/// let matrix      =   VecOfVec::new(   
///                         vec![   vec![ (0isize, 1), (1isize, 6) ], 
///                                 vec![ (0isize, 1), (1isize, 1) ],     ]     
///                     ).ok().unwrap();   
/// // the vector [1, 1]
/// let vector      =   vec![ (0isize, 1), (1isize, 1) ];
/// 
/// // the product
/// let mut product     =   multiply_column_vector_with_matrix_and_return_unsimplified_and_reversed(
///                         vector,
///                         & matrix,
///                         PrimeOrderField::new(7), // the finite field of order 7
///                         OrderOperatorAuto, // this compares tuples according to their lexicographic order
///                     );
/// // the product should return the following entries in sequence: (1,1), (1,1), (0,1), (0,6)       
/// itertools::assert_equal( product, vec![(1,1), (1,1), (0,6), (0,1)] )
/// ```
pub fn multiply_column_vector_with_matrix_and_return_unsimplified_and_reversed 
            < Matrix, RingOperator, SparseVecIter, OrderOperatorForRowIndices> 
            ( 
                sparse_vec:                         SparseVecIter,
                matrix:                             Matrix, 
                ring_operator:                      RingOperator,
                order_operator_for_row_indices:     OrderOperatorForRowIndices,
            ) 
            ->
            LinearCombinationUnsimplified
                < Matrix::ColumnReverse, RingOperator, ReverseOrder< OrderOperatorForRowIndices > >            
    where 
        Matrix:                             MatrixOracle,
        Matrix::ColumnEntry:                KeyValSet,
        RingOperator:                       Clone + SemiringOperations< Element = Matrix::Coefficient >, // ring operators must typically implement Clone for operations such as simplification and dropping zeros
        OrderOperatorForRowIndices:                      Clone + JudgePartialOrder<  Matrix::ColumnEntry >, // order comparators must often implement Clone when one wishes to construct new, ordered objects
        SparseVecIter:                      IntoIterator,
        SparseVecIter::Item:                KeyValGet < Key = Matrix::ColumnIndex, Val = Matrix::Coefficient >,
{
    // LinearCombinationUnsimplified{
    //     linear_combination_unsimplified:
            hit_merge_by_predicate(
                sparse_vec
                        .into_iter()
                        .map(   |x| 
                                matrix
                                    .column_reverse(& x.key() )
                                    .scale_by( x.val(), ring_operator.clone() )
                        ),
                    ReverseOrder::new( order_operator_for_row_indices )
            )        
    // }
}


/// Returns the product of a matrix with a vector; the result is a linear combination of the matrix's columns.
/// 
/// The resulting iterator is a simplified linear combination of columns.  It returns entries in **strictly** descending order of index, provided that
/// [matrix.column_reverse(&column_index)](crate::algebra::matrices::query::MatrixOracle) returns entries in 
/// ascending order according to index, consistent with the user-provided [order operator](crate::utilities::order).
/// 
/// In particular, no two consequtive indices have the same
/// index.
/// 
/// # Examples
/// 
/// If `A` is a **row**-major representation of a matrix, and `v` is a sparse vector, then this function returns `A v`, a linear combination of **columns** of `A`.
/// If instead the representation is column-major, then this function returns `vA`.
/// 
/// ```
/// use oat_rust::algebra::matrices::operations::multiply::multiply_column_vector_with_matrix_and_return_reversed;
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::rings::types::field_prime_order::PrimeOrderField;
/// use oat_rust::utilities::order::OrderOperatorByKey;
/// 
/// // a row-major matrix representation with usize row indices and isize column indices
/// let matrix      =   VecOfVec::new(   
///                         vec![   vec![ (0isize, 1), (1isize, 6) ], 
///                                 vec![ (0isize, 1), (1isize, 1) ],     ]     
///                     ).ok().unwrap();
/// // the vector [1, 1]
/// let vector      =   vec![ (0isize, 1), (1isize, 1) ];
/// 
/// // the product
/// let product     =   multiply_column_vector_with_matrix_and_return_reversed(
///                         vector,
///                         & matrix,
///                         PrimeOrderField::new(7), // the finite field of order 7
///                         OrderOperatorByKey::new(), // this compares order of entries according to their indices
///                     );
/// // the product should equal vector [0 1]        
/// itertools::assert_equal( product, vec![(1usize, 2)] )
/// ```
pub fn multiply_column_vector_with_matrix_and_return_reversed 
            < Matrix, RingOperator, SparseVecIter, OrderOperatorForRowIndices> 
            ( 
                sparse_vec:                         SparseVecIter,
                matrix:                             Matrix, 
                ring_operator:                      RingOperator,
                order_operator_row_row_indices:     OrderOperatorForRowIndices,
            ) 
            ->
            LinearCombinationSimplified
                < Matrix::ColumnReverse, RingOperator, ReverseOrder< OrderOperatorForRowIndices > >

    where 
        Matrix:                             MatrixOracle<
                                                ColumnEntry:                KeyValSet,
                                                // ColumnReverse:              PeekUnqualified,
                                            >,
        RingOperator:                       Clone + SemiringOperations< Element = Matrix::Coefficient >, // ring operators must typically implement Clone for operations such as simplification and dropping zeros
        OrderOperatorForRowIndices:                      Clone + JudgePartialOrder<  Matrix::ColumnEntry >, // order comparators must often implement Clone when one wishes to construct new, ordered objects
        SparseVecIter:                      IntoIterator,
        SparseVecIter::Item:                KeyValGet < Key = Matrix::ColumnIndex, Val = Matrix::Coefficient >,
{

    multiply_column_vector_with_matrix_and_return_unsimplified_and_reversed( 
            sparse_vec,
            matrix,
            ring_operator.clone(),
            order_operator_row_row_indices,
        )
        // .linear_combination_unsimplified
        .simplify(ring_operator)
}    







































//  =========================================================================================================
//  TESTS
//  =========================================================================================================

//  ---------------------------------------------------------------------------
//  DOC-TEST DRAFTS
//  ---------------------------------------------------------------------------

//  We use the following module to draft doc tests, which are easier to debug here than in doc strings.

#[cfg(test)]
mod doc_tests {
    

    use crate::algebra::matrices::operations::MatrixOracleOperations;


    #[test]
    fn doc_test_vector_matrix_product_major_ascend_unsimplified() {
        use crate::algebra::matrices::operations::multiply::multiply_row_vector_with_matrix_unsimplified;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::rings::types::field_prime_order::PrimeOrderField;
        use crate::utilities::order::OrderOperatorByKey;

        // a row-major matrix representation with usize row indices and isize column indices
        let matrix      =   VecOfVec::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
                                                           vec![ (0isize, 1), (1isize, 1) ],     ]     ).unwrap();
        // the vector [1, 1]
        let vector      =   vec![ (0usize, 1), (1usize, 1) ];
        
        // the product
        let product     =   multiply_row_vector_with_matrix_unsimplified(
                                vector,
                                (& matrix).into_peekable_matrix(),
                                PrimeOrderField::new(7), // the finite field of order 7
                                OrderOperatorByKey::new(), // this compares tuples according to their lexicographic order
                            );
        // the product should return the following sequence of entries: (0,1), (0,1), (1,1), (1,6)   
        itertools::assert_equal( product, vec![(0isize,1usize), (0,1), (1,1), (1,6)] )
    }     

    #[test]
    fn doc_test_vector_matrix_product_major_ascend_simplified() {
        use crate::algebra::matrices::operations::multiply::multiply_row_vector_with_matrix;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::rings::types::field_prime_order::PrimeOrderField;
        use crate::utilities::order::OrderOperatorByKey;

        // a row-major matrix representation with usize row indices and isize column indices
        let matrix      =   VecOfVec::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
                                                           vec![ (0isize, 1), (1isize, 1) ],     ]     ).unwrap();
        // the vector [1, 1]
        let vector      =   vec![ (0usize, 1), (1usize, 1) ];
        
        // the product
        let product     =   multiply_row_vector_with_matrix(
                                vector,
                                ( & matrix ).into_peekable_matrix(),
                                PrimeOrderField::new(7), // the finite field of order 7
                                OrderOperatorByKey::new(), // this compares tuples according to their lexicographic order
                            );
        // the product should equal vector [2 0]        
        itertools::assert_equal( product, vec![(0isize, 2)] )
    }    

    #[test]
    fn doc_test_vector_matrix_product_minor_descend_unsimplified() {
        use crate::algebra::matrices::operations::multiply::multiply_column_vector_with_matrix_and_return_unsimplified_and_reversed;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::rings::types::field_prime_order::PrimeOrderField;
        use crate::utilities::order::OrderOperatorAuto;

        // a row-major matrix representation with usize row indices and isize column indices
        let matrix      =   VecOfVec::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
                                                           vec![ (0isize, 1), (1isize, 1) ],     ]     ).unwrap();;
        // the vector [1, 1]
        let vector      =   vec![ (0isize, 1), (1isize, 1) ];
        
        // the product
        let product     =   multiply_column_vector_with_matrix_and_return_unsimplified_and_reversed(
                                vector,
                                ( & matrix ).into_peekable_matrix(),
                                PrimeOrderField::new(7), // the finite field of order 7
                                OrderOperatorAuto, // this compares tuples according to their lexicographic order
                            );
        // the product should return the following entries in sequence: (1,1), (1,1), (0,1), (0,6)       
        itertools::assert_equal( product, vec![(1,1), (1,1), (0,6), (0,1)] )
    }

    #[test]
    fn doc_test_vector_matrix_product_minor_descend_simplified() {
        use crate::algebra::matrices::operations::multiply::multiply_column_vector_with_matrix_and_return_reversed;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::rings::types::field_prime_order::PrimeOrderField;
        use crate::utilities::order::OrderOperatorByKey;

        // a row-major matrix representation with usize row indices and isize column indices
        let matrix      =   VecOfVec::new(   vec![   vec![ (0isize, 1), (1isize, 6) ], 
                                                           vec![ (0isize, 1), (1isize, 1) ],     ]     ).unwrap();;
        // the vector [1, 1]
        let vector      =   vec![ (0isize, 1), (1isize, 1) ];
        
        // the product
        let product     =   multiply_column_vector_with_matrix_and_return_reversed(
                                vector,
                                ( & matrix ).into_peekable_matrix(),
                                PrimeOrderField::new(7), // the finite field of order 7
                                OrderOperatorByKey::new(), // this compares tuples according to their lexicographic order
                            );
        // the product should equal vector [0 2]        
        itertools::assert_equal( product, vec![(1usize, 2)] )
    }
    
    







    #[test]
    fn doc_test_matrix_product_row_unsimplified() {

        // Import the necessary traits and structs
        use crate::algebra::rings::types::native::RingOperatorForNativeRustNumberType;
        use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
        use crate::algebra::matrices::query::MatrixOracle;
        use crate::algebra::matrices::types::product::ProductMatrix;
        use crate::algebra::matrices::types::packet::MatrixAlgebraPacket;
        use std::iter::FromIterator;

        // Define the operator for the coefficient ring
        let ring_operator = RingOperatorForNativeRustNumberType::<f64>::new();     

        // Define matrix A, a 2x2 matrix with 1's above the diagonal and 0's below
        let matrix_1   =   VecOfVec::new(    vec![ vec![ (0,1.), (1,1.) ], vec![ (1,1.) ] ]    ).unwrap();
        // Define matrix B, a 2x2 matrix with 1's above the diagonal and 0's below
        let matrix_2   =   VecOfVec::new(    vec![ vec![ (0,1.), (1,1.) ], vec![ (1,1.) ] ]    ).unwrap();
        // Define references to the two arrays 
        let matrix_1_packet    =   MatrixAlgebraPacket::with_default_order_and_f64_coefficients( &matrix_1 );
        let matrix_2_packet    =   MatrixAlgebraPacket::with_default_order_and_f64_coefficients( &matrix_2 );         
        // Define the lazy product of A and B
        let product     =   ProductMatrix::new( 
                                    matrix_1_packet,                 // matrix A
                                    matrix_2_packet,                 // matrix B
                                );
                            
        // First, check the output iterators themselves (these are wrapped in `Simplify` structs).
        // ---------------------------------------------------------------------------------------
                            
        // Use this object to create a vector-of-vectors that represents the product
        let output =   Vec::from_iter(
                                                (0..2)
                                                    .map(   |x| 
                                                            Vec::from_iter(
                                                                    product
                                                                    .row(&x)
                                                                )
                                                        )
                                            );
                                        
        // Check that the answer is correct                                    
        assert_eq!(     
                output,
                vec![ vec![(0,1.), (1,2.)], vec![(1,1.)]  ]
            ); 
        
        // Second, check the underlying, unsimplified iterators (these are are contained inside the `Simplify` structs).
        // -------------------------------------------------------------------------------------------------------------
        
        // Use this object to create a vector-of-vectors that represents the product
        let output =   Vec::from_iter(
                                                (0..2)
                                                    .map(   |x| 
                                                            Vec::from_iter(
                                                                    product
                                                                    .row(&x)
                                                                    .unsimplified // this unwraps the inner iterator
                                                                )
                                                        )
                                            );
                                        
        // Check that the answer is correct                                    
        assert_eq!(     
                output,
                vec![ vec![(0,1.), (1,1.), (1,1.)], vec![(1,1.)]  ]
            ); 
    }

}



//  ---------------------------------------------------------------------------
//  TEST
//  ---------------------------------------------------------------------------


#[cfg(test)]
mod tests {
    use std::iter::FromIterator;

    use itertools::Itertools;

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use crate::algebra::matrices::types::packet::MatrixAlgebraPacket;
    use crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
    use crate::algebra::rings::types::native::{RingOperatorForNativeRustNumberType};
    use crate::utilities::order::OrderOperatorByLessThan;


    // Test that lazy MATRIX-MATRIX multiplication works properly.
    #[test]
    pub fn matrix_by_matrix_multiply_major_ascend_test_1() {

        // Define the operator for the coefficient ring
        let ring_operator = RingOperatorForNativeRustNumberType::<f64>::new();     
        
        // Define the factor A, a 2x2 matrix with 1's above the diagonal and 0's below.
        let matrix_1   =   VecOfVec::new(
                                                    // MajorDimension::Row,
                                             vec![ vec![ (0,1.), (1,1.) ], vec![ (1,1.) ] ]
                                                ).ok().unwrap();

        // Define the factor B, a 2x2 matrix with 1's above the diagonal and 0's below.
        let matrix_2   =   VecOfVec::new(
                                                    // MajorDimension::Row,
                                             vec![ vec![ (0,1.), (1,1.) ], vec![ (1,1.) ] ]
                                                ).ok().unwrap();
        let matrix_1_packet    =   MatrixAlgebraPacket::with_default_order_and_f64_coefficients( &matrix_1 );
        let matrix_2_packet    =   MatrixAlgebraPacket::with_default_order_and_f64_coefficients( &matrix_2 );        

        // Define the lazy product of A and B
        let product     =   ProductMatrix::new( 
                                                                        matrix_1_packet, // first borrow is to create something that implements the oracle trait, second borrow is needed to meet requirements of the function syntax
                                                                        matrix_2_packet,
                                                                    );

        // CHECK THE MAJOR ITERATOR ITSELF

        // Use this object to create a vector-of-vectors that represents the product 
        let output =   Vec::from_iter(
                                                (0..2)
                                                    .map(   |x| 
                                                            Vec::from_iter(
                                                                    product
                                                                    .row(&x)
                                                                )
                                                        )
                                            );
        
        // Check that the answer is correct                                    
        assert_eq!(     
                output,
                vec![ vec![(0,1.), (1,2.)], vec![(1,1.)]  ]
            );     

        // CHECK THE UNDERLYING **UNSIMPLIFIED** ITERATOR   
        
        // Use this object to create a vector-of-vectors that represents the product 
        let output =   Vec::from_iter(
                                                (0..2)
                                                    .map(   |x| 
                                                            Vec::from_iter(
                                                                    product
                                                                    .row(&x)
                                                                    .unsimplified
                                                                )
                                                        )
                                            );
        
        // Check that the answer is correct                                    
        assert_eq!(     
                output,
                vec![ vec![(0,1.), (1,1.), (1,1.)], vec![(1,1.)]  ]
            );           
        
    }


    // Test that lazy MATRIX-MATRIX multiplication works properly - another way.
    #[test]
    pub fn matrix_by_matrix_multiply_major_ascend_test_2() { 
        
        // Define factor A 
        let matrix_1 =   
            VecOfVec::<usize,i64>::new(
                    vec![ 
                            vec![                  ],
                            vec![ (0,0)            ],                                                            
                            vec![ (0,4),   (1,-1)  ], 
                            vec![ (0,5),   (1,7)   ], 
                            vec![ (0,11),  (1,13)  ],                                                         
                    ]
                ).ok().unwrap();
        // Define factor B 
        let matrix_2 =    
            VecOfVec::<usize,i64>::new(
                    vec![                                                                                                                           
                            vec![ (0,1),   (1,2),   (2,3) ], 
                            vec![ (0,4),   (1,5),   (2,6) ],
                            vec![                         ],  
                            vec![ (0,0)                   ],                                                              
                        ]
                ).ok().unwrap();
        let matrix_1_packet    =   MatrixAlgebraPacket::with_default_order_and_i64_coefficients( &matrix_1 );
        let matrix_2_packet    =   MatrixAlgebraPacket::with_default_order_and_i64_coefficients( &matrix_2 );        

        // This is the *actual* product A*B, if A and B are row-major.
        let product_true: Vec<Vec<(usize,i64)>>   =   
                    vec![ 
                            vec![                                   ], 
                            vec![                                   ],                                                             
                            vec![          (1,3),     (2,6)         ], 
                            vec![ (0,33),  (1,45),    (2,57)        ],
                            vec![ (0,63),  (1,87),    (2,111)       ],                                                             
                        ];                                              

        // Define the lazy product of A and B
        let product     =   ProductMatrix::new( 
                                                                matrix_1_packet,  // first borrow is to create something that implements the oracle trait, second borrow is needed to meet requirements of the function syntax
                                                                matrix_2_packet, 
                                                        );


        // Check that each row of the lazy MATRIX product agrees with each row of the true product.
        for row_index in 0 .. 5 {
            let vec_lazy   =   product
                                                .row( &row_index )
                                                .collect_vec();
            let vec_true =  product_true[ row_index ].clone();
            assert_eq!( vec_lazy, vec_true )
        }   
        
    }  
    

   // Test that lazy MATRIX-VECTOR multiplication works properly.
   #[test]
   pub fn matrix_by_vector_multiply_major_ascend_test_2() {

       // Define the operator for the coefficient ring
       let ring_operator = RingOperatorForNativeRustNumberType::<i32>::new();     
       
       // Define factor A 
       let matrix_1 =   
           VecOfVec::new(
                   vec![ 
                           vec![                  ],
                           vec![ (0,0)            ],                                                            
                           vec![ (0,4),   (1,-1)  ], 
                           vec![ (0,5),   (1,7)   ], 
                           vec![ (0,11),  (1,13)  ],                                                         
                   ]
               ).ok().unwrap();
       // Define factor B 
       let matrix_2 =    
           VecOfVec::new(
                   vec![                                                                                                                           
                           vec![ (0,1),   (1,2),   (2,3) ], 
                           vec![ (0,4),   (1,5),   (2,6) ],
                           vec![                         ],  
                           vec![ (0,0)                   ],                                                              
                       ]
                ).ok().unwrap();
       // This is the *actual* product A*B, if A and B are row-major.
       let product_true   =   
            vec![ 
                    vec![                                   ], 
                    vec![                                   ],                                                             
                    vec![          (1,3),     (2,6)         ], 
                    vec![ (0,33),  (1,45),    (2,57)        ],
                    vec![ (0,63),  (1,87),    (2,111)       ],                                                             
                ];                                                                             


       // Check that each row of the lazy MATRIX product agrees with each row of the true product.
       for row_index in 0 .. 5 {
           let vec_lazy   =   multiply_row_vector_with_matrix(
                                                    (& matrix_1).row( &row_index ).clone(),               
                                                    & matrix_2,  // first borrow is to create something that implements the oracle trait, second borrow is needed to meet requirements of the function syntax
                                                    ring_operator,
                                                    OrderOperatorByLessThan::new( |x:&(i32, i32), y:&(i32, i32)| x.key() < y.key()  )
                                                )
                                               .collect_vec();
           let vec_true =  product_true[ row_index ].clone();
           assert_eq!( vec_lazy, vec_true )
       }   
       
   } 

}