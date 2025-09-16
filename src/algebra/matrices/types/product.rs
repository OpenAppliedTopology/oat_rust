//! The (lazy) product of two matrices
//! 
//! This module provides a [ProductMatrix] struct, which represents the product of two matrices.
//! This object is lazy, in the sense that it stores a copy of each matrix, and only calcuates
//! entries of the product when they are requested.
//! 
//! # See also
//! 
//! The [matrix multiplication module](crate::algebra::matrices::operations::multiply).



//  MATRIX - MATRIX MULTIPLICATION
//  ===========================================================================

use derive_new::new;

use crate::algebra::matrices::query::{MatrixAlgebra, MatrixOracle};
use crate::algebra::matrices::operations::MatrixOracleOperations;
use crate::algebra::rings::traits::SemiringOperations;
use crate::algebra::vectors::entries::{KeyValGet, KeyValSet};
use crate::algebra::vectors::operations::{sort_unstable_by_order_operator_on_indices, LinearCombinationSimplified, VectorOperations};

use crate::utilities::order::ReverseOrder;

use derive_getters::Dissolve;



/// Product of two matrices `A * B`, computed in lazy fashion
/// 
/// - `A` and `B` do not have to have the same type. They can be any types that implement the [MatrixAlgebra] trait.
/// - The column indices of `A` and the row indices of `B` must have the same type.
/// - Check the source code to see how rows, columns, and entries are computed.
///
/// # What data is stored
/// 
/// This struct contains only `A` and `B`, no extra information is stored.
///
/// # Structural nonzeros
/// 
/// If `P` is the product `A * B`, then `P[i,j]` is considered structurally nonzero if and only if `P[i,j]`
/// is actually nonzero. Concretely, we use the `is_zero` method from the ring operator associated with matrix `B`
/// to determine whether each entry is zero.
/// 
/// # Note on trait bounds
/// 
/// This trait requires the row entries of the lefthand matrix to implement [crate::algebra::vectors::entries::KeyValSet]. 
/// This requirement exists to facilitate computation of dot-products of rows of the lefthand matrix with columns of the righthand matrix.
/// It's likely that we could find a way to remove this requirement.
/// 
/// # Performance
/// 
/// We could achieve better performance for the method `self.structural_nonzero_entry(i,j)` if we required the user to
/// ensure that the order operator for the column indices of `A` implemented the same linear order on indices as the 
/// order operator for the row indices of `B`. If we did make that requirement, then we could compute `A.row(i).dot( B.column(j), .. )`
/// relatively efficiently, using the `dot` operation from the [crate::algebra::vectors::operations::VectorOperations] trait.
/// 
/// Instead, we first collect the entries of `B.column(j)` in a vector `v`, then sort `v` using the order operator on column indices
/// provided by `A`. This ensures that `v` and `A.row(i)` are sorted in a consistent order. Then we compute `A.row(i).dot( v, .. )`.
/// 
/// So, the real cost here comes in the extra time and space it takes to dump then entries of `B.col(j)` into a vector, and sort it.
/// 
/// The user can of course achieve better performance by calling `A.row(i).dot( B.column(j), .. )` directly
/// 
/// # Examples
/// 
/// ```
/// use oat_rust::algebra::rings::types::native::RingOperatorForNativeRustNumberType;
/// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;
/// use oat_rust::algebra::matrices::types::packet::MatrixAlgebraPacket;
/// use oat_rust::algebra::matrices::query::{MatrixOracle, MatrixAlgebra};
/// use oat_rust::algebra::matrices::types::product::ProductMatrix;
/// use oat_rust::algebra::vectors::entries::KeyValGet;
/// use oat_rust::utilities::order::OrderOperatorByLessThan;
/// use std::iter::FromIterator;
/// use itertools::Itertools;
/// 
/// // Initialize variables
/// // ---------------------------------------------------------------------------------------
/// 
/// // Define matrix A, a 2x2 matrix with 1's above the diagonal and 0's below
/// let matrix_1        =   VecOfVec::new(      vec![ 
///                                                     vec![   (0,1.),   (1,1.)    ], 
///                                                     vec![             (1,1.)    ],
///                                             ]    
///                         ).ok().unwrap();
/// // Define matrix B, a 2x2 matrix with 1's above the diagonal and 0's below
/// let matrix_2        =   VecOfVec::new(      vec![ 
///                                                     vec![   (0,1.),   (1,1.)    ], 
///                                                     vec![             (1,1.)    ],
///                                             ]       
///                         ).ok().unwrap();
/// // Place the matrices inside "packets" that contain extra info about how to add/multiply coefficients and order nonzero entries
/// let matrix_1_pac    =   MatrixAlgebraPacket::with_default_order_and_f64_coefficients( &matrix_1 );
/// let matrix_2_pac    =   MatrixAlgebraPacket::with_default_order_and_f64_coefficients( &matrix_2 );
///                     
/// // Define the lazy product of A and B
/// // ---------------------------------------------------------------------------------------
/// 
/// let product         =   ProductMatrix::new( 
///                             matrix_1_pac,                 // matrix A
///                             matrix_2_pac,                 // matrix B
///                         );
///                     
/// // Check that the rows of the product are correct.
/// // ---------------------------------------------------------------------------------------
///                     
/// // Use this object to create a vector-of-vectors that represents the product
/// let output          =   vec![
///                             product.row( & 0 ).collect_vec(),
///                             product.row( & 1 ).collect_vec(),
///                         ];
/// 
/// // Check that the answer is correct                                    
/// assert_eq!(     
///         output,
///         vec![ 
///                 vec![   (0,1.),   (1,2.)    ], 
///                 vec![             (1,1.)    ],
///         ]    
///     ); 
/// ```
#[derive(new,Clone,Copy,Debug,Dissolve,Eq,PartialEq,Ord,PartialOrd)]
pub struct ProductMatrix < MatrixLeft, MatrixRight > 
{
    matrix_left:    MatrixLeft,     
    matrix_right:   MatrixRight,        
}


impl < MatrixLeft, MatrixRight > 
    
    ProductMatrix
        < MatrixLeft, MatrixRight > 

    where
        MatrixLeft:     MatrixAlgebra<
                            RowEntry:                           KeyValSet, // this "assymetric" requirement exists to support a single dot product computation `vec_a.dot( .. )`, see below
                            ColumnEntry:                        KeyValSet,      
                        >,
        MatrixRight:    MatrixAlgebra< 
                            Coefficient                     =   MatrixLeft::Coefficient, // matrices must have same coefficients  
                            RowIndex                        =   MatrixLeft::ColumnIndex, //    
                            RowEntry:                           KeyValSet,                 
                        >,
{
    /// Computes the structural nonzero entry in `(row,column)`
    /// 
    /// This method is provided for users who seek a (slight) performance advantage over the the method `structural_nonzero_entry` implemented on [ProductMatrix] for the [MatrixOracle] trait.
    /// The user **must ensure that the order operators for the column indices of the left matrix and the row indices of the right matrix induces the same order**.
    /// If they do not, then the returned value may be incorrect. 
    /// 
    /// By contrast, the value returned by the `structural_nonzero_entry` is always correct,
    /// provided that there are no errors in the implementation of the [MatrixOracle] and [MatrixAlgebra] traits for both matrices.
    /// 
    /// See the documentation on [MatrixOracle::structural_nonzero_entry] for [ProductMatrix] for a discussion on performance considerations.
    pub fn structural_nonzero_entry_fast_unsafe(&   self, row:   & MatrixLeft::RowIndex, column: & MatrixRight::ColumnIndex ) ->  Option< MatrixLeft::Coefficient > {
        let ring_operator       =   self.matrix_left.ring_operator();
        let order_operator  =   self.matrix_left.order_operator_for_column_indices();        
        let vec_a   =   self.matrix_left.row( row );
        let vec_b   =   self.matrix_right.column( column );
        let dot                 =   vec_a.dot( 
                                        vec_b, 
                                        ring_operator.clone(), 
                                        order_operator 
                                    );
        if ring_operator.is_0( dot.clone() ) {
            return None
        } else {
            Some( dot )
        }
    }
}






impl < MatrixLeft, MatrixRight > 

    MatrixOracle for 
    
    ProductMatrix
        < MatrixLeft, MatrixRight > 

    where
        MatrixLeft:     MatrixAlgebra<
                            RowEntry:                           KeyValSet, // this "assymetric" requirement exists to support a single dot product computation `vec_a.dot( .. )`, see below
                            ColumnEntry:                        KeyValSet,      
                        >,
        MatrixRight:    MatrixAlgebra< 
                            Coefficient                     =   MatrixLeft::Coefficient, // matrices must have same coefficients  
                            RowIndex                        =   MatrixLeft::ColumnIndex, //    
                            RowEntry:                           KeyValSet,                 
                        >,
{
    type Coefficient            =   MatrixRight::Coefficient;    // The type of coefficient stored in each entry of the matrix    
    type RowIndex               =   MatrixLeft::RowIndex;    // The type key used to look up rows.  Matrices can have rows indexed by anything: integers, simplices, strings, etc.
    type ColumnIndex            =   MatrixRight::ColumnIndex;    // The type of column indices    
    type RowEntry               =   MatrixRight::RowEntry; // The type of entries in each row; these are essentially pairs of form `(column_index, coefficient)`
    type ColumnEntry            =   MatrixLeft::ColumnEntry; // The type of entries in each column; these are essentially pairs of form `(row_index, coefficient)`
    type Row                    =   LinearCombinationSimplified
                                        < MatrixRight::Row, MatrixRight::RingOperator, MatrixRight::OrderOperatorForRowEntries >;  // What you get when you ask for a row.
    type RowReverse             =   LinearCombinationSimplified
                                        < MatrixRight::RowReverse, MatrixRight::RingOperator, ReverseOrder< MatrixRight::OrderOperatorForRowEntries > >;  // What you get when you ask for a row with the order of entries reversed
    type Column                 =   LinearCombinationSimplified
                                        < MatrixLeft::Column, MatrixLeft::RingOperator, MatrixLeft::OrderOperatorForColumnEntries >;  // What you get when you ask for a column
    type ColumnReverse          =   LinearCombinationSimplified
                                        < MatrixLeft::ColumnReverse, MatrixRight::RingOperator, ReverseOrder< MatrixLeft::OrderOperatorForColumnEntries > >;  // What you get when you ask for a column with the order of entries reversed                             

    fn structural_nonzero_entry(&   self, row:   & Self::RowIndex, column: & Self::ColumnIndex ) ->  Option< Self::Coefficient > {
        let ring_operator       =   self.matrix_left.ring_operator();
        let order_operator  =   self.matrix_left.order_operator_for_column_indices();        
        let vec_a   =   self.matrix_left.row( row );
        let vec_b   =   self.matrix_right.column( column );
        let mut vec_b = vec_b.collect::<Vec<_>>();
        sort_unstable_by_order_operator_on_indices( &mut vec_b, order_operator.clone() ); // ensure that vec_b is sorted in an order consistent with a

        let dot                 =   vec_a.dot( 
                                        vec_b, 
                                        ring_operator.clone(), 
                                        order_operator 
                                    );
        if ring_operator.is_0( dot.clone() ) {
            return None
        } else {
            Some( dot )
        }
    }

    fn has_column_for_index(  &   self, index: & Self::ColumnIndex)   -> bool {
        self.matrix_right.has_column_for_index( index )
    }
    fn has_row_for_index(     &   self, index: & Self::RowIndex   )   -> bool {
        self.matrix_left.has_row_for_index( index )
    }

    fn row(                     &   self, index: &Self::RowIndex   )   -> Self::Row {
        let coefficients    =   self.matrix_left.row( index );
        coefficients.multiply_self_as_a_row_vector_with_matrix(
            &self.matrix_right, 
        )
    }
    fn row_reverse(             &   self, index: &Self::RowIndex   )   -> Self::RowReverse          { 
        let coefficients    =   self.matrix_left.row( index );
        coefficients.multiply_self_as_a_row_vector_with_matrix_and_return_entries_in_reverse_order(
            &self.matrix_right, 
        )    
    }
    fn column(                  &   self, index: &Self::ColumnIndex)   -> Self::Column              { 
        let coefficients    =   self.matrix_right.column( index );
        coefficients.multiply_self_as_a_column_vector_with_matrix(
            &self.matrix_left, 
        )
    }
    fn column_reverse(              &   self, index: &Self::ColumnIndex)   -> Self::ColumnReverse              { 
        let coefficients    =   self.matrix_right.column( index );
        coefficients.multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order_custom(
            &self.matrix_left, 
            self.matrix_right.ring_operator(),
            self.matrix_left.order_operator_for_column_entries(),
        )

        
    } 

}




impl < MatrixLeft, MatrixRight > 

    MatrixAlgebra for 
    
    ProductMatrix
        < MatrixLeft, MatrixRight > 

    where
        MatrixLeft:     MatrixAlgebra<
                            RowEntry:                           KeyValSet,
                            ColumnEntry:                        KeyValSet,      
                        >,
        MatrixRight:    MatrixAlgebra< 
                            Coefficient                     =   MatrixLeft::Coefficient, // matrices must have same coefficients  
                            RowIndex                        =   MatrixLeft::ColumnIndex, //    
                            RowEntry:                           KeyValSet,                 
                        >,
{
    type OrderOperatorForColumnEntries      =   MatrixLeft::OrderOperatorForColumnEntries;
    type OrderOperatorForColumnIndices      =   MatrixRight::OrderOperatorForColumnIndices;
    type OrderOperatorForRowEntries         =   MatrixRight::OrderOperatorForRowEntries;
    type OrderOperatorForRowIndices         =   MatrixLeft::OrderOperatorForRowIndices;
    type RingOperator                       =   MatrixRight::RingOperator;

    fn order_operator_for_column_entries( &self ) -> Self::OrderOperatorForColumnEntries { 
        self.matrix_left.order_operator_for_column_entries()
    }
    fn order_operator_for_column_indices( &self ) -> Self::OrderOperatorForColumnIndices { 
        self.matrix_right.order_operator_for_column_indices()
    }   
    fn order_operator_for_row_entries( &self ) -> Self::OrderOperatorForRowEntries { 
        self.matrix_right.order_operator_for_row_entries()
    }
    fn order_operator_for_row_indices( &self ) -> Self::OrderOperatorForRowIndices { 
        self.matrix_left.order_operator_for_row_indices()
    }        
    fn ring_operator( &self ) -> Self::RingOperator {
        self.matrix_right.ring_operator()
    }
}





impl < MatrixLeft, MatrixRight > 

    MatrixOracleOperations for 
    
    ProductMatrix
        < MatrixLeft, MatrixRight > 

{}        