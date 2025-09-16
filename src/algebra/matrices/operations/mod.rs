//! Matrix multiplication, inversion, factorization, etc.
//! 
//!
//! - [accessing entries in rows and columns](crate::algebra::matrices::query)
//! - [matrix and vector multiplication](crate::algebra::matrices::operations::multiply)
//! - [solving systems of equations](crate::algebra::matrices::operations::solve)
//! - [matrix factorization](crate::algebra::matrices::operations::umatch)
//! - [matrix inversion](crate::algebra::matrices::operations::invert)
//! - [kernels](crate::algebra::matrices::operations::umatch::row_major::Umatch::kernel)
//! - [images](crate::algebra::matrices::operations::umatch::row_major::Umatch::image)
//! - [transforming and combining vectors (scale, add, drop zeros, reindex, etc.)](crate::algebra::vectors::operations)

pub mod multiply; 
pub mod invert;
pub mod vec_of_vec_reduction;
pub mod solve;
pub mod transform_entry_wise;
pub mod transform_vector_wise;
pub mod umatch;
// pub mod display;
pub mod transpose;
pub mod combine_rows_and_columns;
pub mod iterate_rows_and_columns;


// =========================================================================

use iterate_rows_and_columns::{SequenceOfReverseRows, SequenceOfRows};

use crate::algebra::matrices::operations::transform_vector_wise::PutbackIteratorMatrix;
use crate::algebra::rings::traits::DivisionRingOperations;
use crate::algebra::vectors::entries::{KeyValGet, KeyValPair, KeyValSet};

use crate::algebra::matrices::types::peekable::PeekableMatrix;
use crate::algebra::vectors::operations::VectorOperations;
// use crate::utilities::iterators::general::PeekUnqualified;
use crate::utilities::order::{JudgeOrder, OrderOperatorByKeyCustom};

use self::umatch::row_major::Umatch;
use super::types::product::ProductMatrix;

use super::query::{ MatrixAlgebra, MatrixOracle, };
use super::operations::combine_rows_and_columns::{LinearCombinationOfColumns, LinearCombinationOfColumnsReverse, LinearCombinationOfRows, LinearCombinationOfRowsReverse};
use super::operations::iterate_rows_and_columns::{SequenceOfColumns, SequenceOfReverseColumns,};
use super::types::packet::MatrixAlgebraPacket;
use super::types::transpose::{Transpose, OrderAntiTranspose};
use super::types::reverse::ReverseMatrix;

use std::hash::Hash;





/// Convenient methods for matrix operations, including multiplication, decomposition, etc.
pub trait MatrixOracleOperations {




    /// Lefthand multiplication with another matrix
    /// 
    /// Returns `self * other_matrix`
    fn multiply_on_the_left_of< Othr >( 
                self, 
                other: Othr,
            )
        ->
        ProductMatrix< Self, Othr > 

        where 
            Self:               Sized + MatrixAlgebra,
            Othr:               MatrixAlgebra< Coefficient = Self::Coefficient, RowIndex = Self::ColumnIndex, OrderOperatorForRowIndices = Self::OrderOperatorForColumnIndices >,
    {
        ProductMatrix::new( self, other )
    }

    /// Righthand matrix multiplication with another matrix
    /// 
    /// Returns `other_matrix * self`
    fn multiply_on_the_right_of< Othr >( 
                self, 
                other: Othr,
            )
        ->
        ProductMatrix< Othr, Self > 

        where 
            Self:               Sized + MatrixAlgebra,
            Othr:               MatrixAlgebra< Coefficient = Self::Coefficient, ColumnIndex = Self::RowIndex, OrderOperatorForColumnIndices = Self::OrderOperatorForRowIndices >,
    {
        ProductMatrix::new( other, self )
    } 

    /// Returns a [Umatch] decomposition of `self`
    /// 
    /// The argument `iter_row_index` must run over all row indices in *strictly descending order*,
    /// as determined by `self.order_operator_for_row_indices()`.
    /// 
    /// # Example
    /// 
    /// In this example we compute a U-match decomposition of a matrix `D`, then use
    /// the decomposition to solve `Dx = b` for `x`.
    /// 
    /// ```
    /// use crate::oat_rust::algebra::matrices::operations::MatrixOracleOperations;
    /// use oat_rust::algebra::matrices::operations::umatch::row_major::Umatch;
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;   
    /// use oat_rust::algebra::matrices::types::packet::MatrixAlgebraPacket;     
    /// use itertools::Itertools;
    /// 
    /// // DEFINE THE MATRIX
    /// // ===============================
    /// let matrix          =   VecOfVec::new( 
    ///                                         vec![   
    ///                                                     vec![(0,true), (1,true), (2,true)],
    ///                                                     vec![                            ], 
    ///                                                     vec![                    (2,true)], 
    ///                                         ] 
    ///                                     ).ok().unwrap();
    /// let matrix          =   MatrixAlgebraPacket::with_default_order_and_boolean_coefficients(
    ///                             &matrix 
    ///                         );
    ///                                 
    /// // COMPUTE U-MATCH
    /// // ===============================
    ///                                 
    /// let umatch
    ///     =   matrix.into_umatch( 
    ///             (0..3).rev(), // an iterator that runs over all row indices, from bottom to top
    ///         );        
    /// 
    /// // SOLVE Dx = b FOR x
    /// // ===============================
    /// 
    /// let b   =   [ (0,true), (2,true) ];
    /// let x   =   umatch.solve_dx_equals_b( b.clone() ).unwrap();
    /// let dx  =   umatch.multiply_dx(x);
    /// assert!( dx.eq( b ) );   
    /// ```
    fn into_umatch< IterRowIndex >(
                self, 
                iter_row_index:                            IterRowIndex, 
            )
        ->
        Umatch< Self >  
            where   Self:                   Sized + 
                                            MatrixAlgebra<
                                                Row:            Clone,
                                                RowEntry:       KeyValPair,
                                                ColumnEntry:    KeyValPair,
                                                ColumnIndex:    Hash,
                                                RowIndex:       Hash,
                                                RingOperator:   DivisionRingOperations,
                                            >, 
                    IterRowIndex:           Iterator < Item = Self::RowIndex >,
    {
        Umatch::new( self, iter_row_index, )
    }

    /// Returns a [Umatch] decomposition of `self`
    /// 
    /// - The argument `iter_row_index` must run over all row indices in *strictly descending order*, as determined by `order_operator_for_row_indices`.
    /// - This function constructs a [MatrixAlgebraPacket] from `self` together with the ring operator and order operators for row and column indices.
    ///   Order operators for row and column entries
    ///   are obtained from the user-provided order operators for row and column indices by wrapping them in [OrderOperatorByKeyCustom].
    ///   The give order operators **must conform to the conditions specified in the documentation for [MatrixAlgebraPacket]**. Otherwise the U-match factorization may be
    ///   calculated incorrectly, **without any warnings given**.
    /// 
    /// This function is similar to [into_umatch](Self::into_umatch), 
    /// but allows the user to specify a custom ring operator and custom order operators for row
    /// and column indices. Unlike [into_umatch](Self::into_umatch), this function
    /// does not require `Self` to implement [MatrixAlgebra](crate::algebra::matrices::query::MatrixAlgebra),
    /// only [MatrixOracle](crate::algebra::matrices::query::MatrixOracle).
    /// 
    /// # Example
    /// 
    /// In this example we compute a U-match decomposition of a matrix `D`, then use
    /// the decomposition to solve `Dx = b` for `x`.
    /// 
    /// ```
    /// use crate::oat_rust::algebra::matrices::operations::MatrixOracleOperations;
    /// use oat_rust::algebra::matrices::operations::umatch::row_major::Umatch;
    /// use oat_rust::algebra::matrices::types::vec_of_vec::sorted::VecOfVec;   
    /// use oat_rust::algebra::matrices::types::packet::MatrixAlgebraPacket;  
    /// use oat_rust::algebra::rings::types::field_prime_order::BooleanField;
    /// use oat_rust::utilities::order::OrderOperatorAuto;
    /// use itertools::Itertools;
    /// 
    /// // DEFINE THE MATRIX
    /// // ===============================
    /// let matrix          =   & VecOfVec::new( 
    ///                                         vec![   
    ///                                                     vec![(0,true), (1,true), (2,true)],
    ///                                                     vec![                            ], 
    ///                                                     vec![                    (2,true)], 
    ///                                         ] 
    ///                                     ).ok().unwrap();
    ///                                 
    /// // COMPUTE U-MATCH
    /// // ===============================
    ///                                 
    /// let umatch
    ///     =   matrix.into_umatch_custom( 
    ///             (0..3).rev(),           // an iterator that runs over all row indices, from bottom to top
    ///             BooleanField::new(),    // ring operator
    ///             OrderOperatorAuto,      // order operator for row indices
    ///             OrderOperatorAuto,      // order operator for column indices
    ///         );        
    /// 
    /// // SOLVE Dx = b FOR x
    /// // ===============================
    /// 
    /// let b   =   [ (0,true), (2,true) ];
    /// let x   =   umatch.solve_dx_equals_b( b.clone() ).unwrap();
    /// let dx  =   umatch.multiply_dx(x);
    /// assert!( dx.eq( b ) );   
    /// ```
    fn into_umatch_custom< RingOperator, OrderOperatorForRowIndices, OrderOperatorForColumnIndices, IterRowIndex >(
                self, 
                iter_row_index:                            IterRowIndex, 
                ring_operator:                          RingOperator, 
                order_operator_for_row_indices:         OrderOperatorForRowIndices, 
                order_operator_for_column_indices:      OrderOperatorForColumnIndices,            
            )
        ->
        Umatch 
            <
                MatrixAlgebraPacket<
                    Self,
                    RingOperator,
                    OrderOperatorByKeyCustom < OrderOperatorForColumnIndices, >, // recall that row entries have column indices
                    OrderOperatorForRowIndices,
                    OrderOperatorByKeyCustom < OrderOperatorForRowIndices, >, // recall that column entries have row indices
                    OrderOperatorForColumnIndices,
                > 
            >  
        where   Self:                               Sized + 
                                                    MatrixOracle<
                                                        Row:            Clone,
                                                        RowEntry:       KeyValPair,
                                                        ColumnEntry:    KeyValPair,
                                                        ColumnIndex:    Hash,
                                                        RowIndex:       Hash,
                                                    >, 
                IterRowIndex:                       Iterator < Item = Self::RowIndex >,
                RingOperator:                       Clone + DivisionRingOperations< Element = Self::Coefficient >,
                OrderOperatorForColumnIndices:      Clone + JudgeOrder <  Self::ColumnIndex >,
                OrderOperatorForRowIndices:         Clone + JudgeOrder <  Self::RowIndex >,    

        {
            let packet  =   MatrixAlgebraPacket{
                    matrix:                             self,
                    ring_operator:                      ring_operator,
                    order_operator_for_row_indices:     order_operator_for_row_indices.clone(),
                    order_operator_for_column_indices:  order_operator_for_column_indices.clone(),
                    order_operator_for_row_entries:     OrderOperatorByKeyCustom::<    OrderOperatorForColumnIndices, >::new( order_operator_for_column_indices ),
                    order_operator_for_column_entries:  OrderOperatorByKeyCustom::< OrderOperatorForRowIndices,    >::new( order_operator_for_row_indices    ),
                };
            Umatch::new( packet, iter_row_index, )
        }   

    /// Wraps `self` in a [crate::algebra::matrices::types::transpose::Transpose] struct.
    fn transpose( self ) -> Transpose<Self> 
        where
            Self:   Sized,
    { Transpose::new( self )  } 

    /// Wraps `& self` in a [crate::algebra::matrices::types::transpose::Transpose] struct.
    fn transpose_of_ref( &self ) -> Transpose<&Self> 
        where
    { Transpose::new( & self )  }         

    /// Wraps `self` in an [OrderAntiTranspose] struct.
    /// 
    /// # Caution
    /// 
    /// There are three important differences between [OrderAntiTranspose] and the matrix returned by [antitranspose_deep](crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec::antitranspose_deep), which is available for [VecOfVec](crate::algebra::matrices::types::vec_of_vec::sorted::VecOfVec) matrices. See the [OrderAntiTranspose] documentation for details.
    fn order_antitranspose( self ) -> OrderAntiTranspose<Self> 
        where
            Self:   Sized,
    { OrderAntiTranspose::new( self )  } 

    /// Wraps `& self` in an [OrderAntiTranspose] struct.
    fn order_antitranspose_of_ref( &self ) -> OrderAntiTranspose<&Self> 
        where
            Self:   Sized,    
    { OrderAntiTranspose::new( &self ) }    

    /// Wraps `self` in a [Reverse] struct.
    fn reverse_order_of_rows_and_columns( self ) -> ReverseMatrix<Self> 
        where
            Self:   Sized,
    { ReverseMatrix::new( self )  } 

    /// Wraps `& self` in a [Reverse] struct.
    fn reverse_of_ref( &self ) -> ReverseMatrix<&Self> 
        where
    { ReverseMatrix::new( & self )  }


    /// Multiply with a row vector
    /// 
    /// Calling `self.multiply_with_row_vector( v )` is the same as calling `v.multiply_self_as_a_row_vector_with_matrix( &self )`. For details see
    /// the documentation for [multiply_self_as_a_row_vector_with_matrix](crate::algebra::vectors::operations::VectorOperations::multiply_self_as_a_row_vector_with_matrix).
    fn multiply_with_row_vector< V >( &self, vector: V ) -> LinearCombinationOfRows< Self >
        where
            Self:               Sized + MatrixAlgebra,
            V:                  IntoIterator,
            V::Item:            KeyValGet < Key = Self::RowIndex, Val = Self::Coefficient >,
            Self::RowEntry:     KeyValSet < Key = Self::ColumnIndex, Val = Self::Coefficient >,
            // Self::Row:          PeekUnqualified,
    {
        vector.multiply_self_as_a_row_vector_with_matrix( &self )
    }
    /// Returns `Err(k)` if the vector contains an entry `(k,v)`, where `k` is an invaled row index for `self`
    /// 
    /// Calling `self.multiply_with_row_vector_result( v )` is the same as calling `v.matrix_multiply_as_row( &self )`. For details see
    /// the documentation for [matrix_multiply_as_row_opt](crate::algebra::vectors::operations::VectorOperations::matrix_multiply_as_row_opt).
    fn multiply_with_row_vector_result< V >( &self, vector: V ) -> Result< LinearCombinationOfRows< Self >, Self::RowIndex >
        where
            Self:               Sized + MatrixAlgebra,
            V:                  IntoIterator,
            V::Item:            KeyValGet < Key = Self::RowIndex, Val = Self::Coefficient >,
            Self::RowEntry:     KeyValSet < Key = Self::ColumnIndex, Val = Self::Coefficient >,
            // Self::Row:          PeekUnqualified
    {
        vector.multiply_self_as_a_row_vector_with_matrix_result( &self )
    }   
    /// Returns entries in reverse order
    /// 
    /// Calling `self.multiply_vector_as_row_reverse( v )` is the same as calling `v.matrix_multiply_as_row_reverse( &self )`. For details see
    /// the documentation for [matrix_multiply_as_row_reverse](crate::algebra::vectors::operations::VectorOperations::matrix_multiply_as_row_reverse).
    fn multiply_with_row_vector_and_return_entries_in_reverse_order< V >( &self, vector: V ) -> LinearCombinationOfRowsReverse< Self >
        where
            Self:               Sized + MatrixAlgebra,
            V:                  IntoIterator,
            V::Item:            KeyValGet < Key = Self::RowIndex, Val = Self::Coefficient >,
            Self::RowEntry:     KeyValSet < Key = Self::ColumnIndex, Val = Self::Coefficient >,
            // Self::RowReverse:   PeekUnqualified,  
    {
        vector.multiply_self_as_a_row_vector_with_matrix_and_return_entries_in_reverse_order( &self )
    }
    /// Returns `None` if the vector contains an entry `(k,v)`, where `k` is an invaled row index for `self`
    /// 
    /// Calling `self.multiply_with_row_vector_and_return_entries_in_reverse_order_result( v )` is the same as calling `v.multiply_self_as_a_row_vector_with_matrix_and_return_entries_in_reverse_order_result( &self )`. For details see
    /// the documentation for [multiply_self_as_a_row_vector_with_matrix_and_return_entries_in_reverse_order_result](crate::algebra::vectors::operations::VectorOperations::multiply_self_as_a_row_vector_with_matrix_and_return_entries_in_reverse_order_result).
    fn multiply_with_row_vector_and_return_entries_in_reverse_order_result< V >( &self, vector: V ) -> Result< LinearCombinationOfRowsReverse< Self >, Self::RowIndex >
        where
            Self:               Sized + MatrixAlgebra,
            V:                  IntoIterator,
            V::Item:            KeyValGet < Key = Self::RowIndex, Val = Self::Coefficient >,
            Self::RowEntry:     KeyValSet < Key = Self::ColumnIndex, Val = Self::Coefficient >,
            // Self::RowReverse:   PeekUnqualified,
    {
        vector.multiply_self_as_a_row_vector_with_matrix_and_return_entries_in_reverse_order_result( &self )
    }       


    /// Multiply with a column vector
    /// 
    /// Calling `self.multiply_with_column_vector( v )` is the same as calling `v.multiply_self_as_a_column_vector_with_matrix( &self )`. For details see
    /// the documentation for [matrix_multiply_as_column](crate::algebra::vectors::operations::VectorOperations::matrix_multiply_as_column).
    fn multiply_with_column_vector< V >( &self, vector: V ) -> LinearCombinationOfColumns< Self >
        where
            Self:               Sized + MatrixAlgebra,
            V:                  IntoIterator,
            V::Item:            KeyValGet < Key = Self::ColumnIndex, Val = Self::Coefficient >,
            Self::ColumnEntry:  KeyValSet < Key = Self::RowIndex, Val = Self::Coefficient >,
            // Self::Column:       PeekUnqualified,             
    {
        vector.multiply_self_as_a_column_vector_with_matrix( &self )
    }
    /// Returns `Err(k)` if the vector contains an entry `(k,v)`, where `k` is an invaled column index for `self`
    /// 
    /// Calling `self.multiply_with_column_vector_result( v )` is the same as calling `v.multiply_self_as_a_column_vector_with_matrix_result( &self )`. For details see
    /// the documentation for [multiply_self_as_a_column_vector_with_matrix_result](crate::algebra::vectors::operations::VectorOperations::multiply_self_as_a_column_vector_with_matrix_result).
    fn multiply_with_column_vector_result< V >( &self, vector: V ) -> Result< LinearCombinationOfColumns< Self >, Self::ColumnIndex >
        where
            Self:               Sized + MatrixAlgebra,
            V:                  IntoIterator,
            V::Item:            KeyValGet < Key = Self::ColumnIndex, Val = Self::Coefficient >,
            Self::ColumnEntry:  KeyValSet < Key = Self::RowIndex, Val = Self::Coefficient >,
            // Self::Column:       PeekUnqualified,                
    {
        vector.multiply_self_as_a_column_vector_with_matrix_result( &self )
    }        
    /// Returns entries in reverse order
    /// 
    /// Calling `self.multiply_with_column_vector_reverse( v )` is the same as calling `v.multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order( &self )`. For details see
    /// the documentation for [multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order](crate::algebra::vectors::operations::VectorOperations::multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order).
    fn multiply_with_column_vector_reverse< V >( &self, vector: V ) -> LinearCombinationOfColumnsReverse< Self >
        where
            Self:               Sized + MatrixAlgebra,
            V:                  IntoIterator,
            V::Item:            KeyValGet < Key = Self::ColumnIndex, Val = Self::Coefficient >,
            Self::ColumnEntry:  KeyValSet < Key = Self::RowIndex, Val = Self::Coefficient >,
            // Self::ColumnReverse:PeekUnqualified,                
    {
        vector.multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order( &self )
    }
    /// Returns `None` if the vector contains an entry `(k,v)`, where `k` is an invaled column index for `self`
    /// 
    /// Calling `self.multiply_with_column_vector_reverse_result( v )` is the same as calling `v.multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order( &self )`. For details see
    /// the documentation for [multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order_result](crate::algebra::vectors::operations::VectorOperations::multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order_result).
    fn multiply_with_column_vector_reverse_result< V >( &self, vector: V ) -> Result< LinearCombinationOfColumnsReverse< Self >, Self::ColumnIndex >
        where
            Self:               Sized + MatrixAlgebra,
            V:                  IntoIterator,
            V::Item:            KeyValGet < Key = Self::ColumnIndex, Val = Self::Coefficient >,
            Self::ColumnEntry:  KeyValSet < Key = Self::RowIndex, Val = Self::Coefficient >,
            // Self::ColumnReverse:PeekUnqualified,                
    {
        vector.multiply_self_as_a_column_vector_with_matrix_and_return_entries_in_reverse_order_result( &self )
    }      



    /// Returns the structural nonzero entry with the largest row index in the specified column
    /// 
    /// Returns `None` if the column is empty.
    fn bottom_entry_for_column(&self, column_index: &Self::ColumnIndex) 
        -> 
    Option<Self::ColumnEntry> 
    
    where
        Self:               MatrixAlgebra
    {
        self.column_reverse( column_index ).next()
    }



    /// Wraps `self` in [PutbackIteratorMatrix], which makes rows and columns peekable
    /// 
    /// See the documentation for [PeekableMatrix] for details.
    fn into_put_backable_matrix( self ) -> PutbackIteratorMatrix< Self >
        where
            Self:               Sized
    {
        PutbackIteratorMatrix::new( self )
    }


    /// Wraps `self` in [PeekableMatrix], which makes rows and columns peekable
    /// 
    /// See the documentation for [PeekableMatrix] for details.
    fn into_peekable_matrix( self ) -> PeekableMatrix< Self >
        where
            Self:               Sized
    {
        PeekableMatrix::new( self )
    }

    /// Iterates over columns
    /// 
    /// This function returns an iterator which runs over columns of the matrix in the order specified by the user-provided iterator of column indices.
    fn sequence_of_columns< ColumnIndexIterator >( &self, indices: ColumnIndexIterator ) 
        -> SequenceOfColumns< & Self, ColumnIndexIterator > {
        SequenceOfColumns::new( self, indices )
    }
    
    /// Iterates over columns (order of entries in each column is reversed)
    /// 
    /// This function returns an iterator which runs over columns of the matrix in the order specified by the user-provided iterator of column indices.
    /// The order of entries in each column is reversed.
    fn sequence_of_reverse_columns< ColumnIndexIterator >( &self, indices: ColumnIndexIterator ) 
        -> SequenceOfReverseColumns< & Self, ColumnIndexIterator > {
        SequenceOfReverseColumns::new( self, indices )
    }

    /// Iterates over rows
    /// 
    /// This function returns an iterator which runs over rows of the matrix in the order specified by the user-provided iterator of row indices.    
    fn sequence_of_rows< RowIndexIterator >( &self, indices: RowIndexIterator ) 
        -> SequenceOfRows< & Self, RowIndexIterator > {
        SequenceOfRows::new( self, indices )
    }

    /// Iterates over rows (order of entries in each row is reversed)
    /// 
    /// This function returns an iterator which runs over rows of the matrix in the order specified by the user-provided iterator of row indices.
    /// The order of entries in each row is reversed.    
    fn sequence_of_reverse_rows< RowIndexIterator >( &self, indices: RowIndexIterator ) 
        -> SequenceOfReverseRows< & Self, RowIndexIterator > {
        SequenceOfReverseRows::new( self, indices )
    }




}



// //  Auto-implement this trait on all types
// //  --------------------------------------------------------------
// impl < Matrix > 

//     MatrixOracleOperations for 
    
//     Matrix

// {}    